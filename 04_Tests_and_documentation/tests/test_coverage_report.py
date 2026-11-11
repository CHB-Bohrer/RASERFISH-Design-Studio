"""Coverage-only tests: saved results and synthetic intervals, never new QC runs."""
import copy
import gzip
import hashlib
import json
import re
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import coverage_report as coverage
import studio

INSTALLED = Path(__file__).resolve().parents[1]
EXPERIMENT = '20261019_111448_six_10kb'
AUDIT = Path(__file__).resolve().parents[2] / 'outputs/final_six_target_review_20261021'


def fixture(intervals=((101, 136), (1101, 1136), (2500, 2535)), strand='+'):
    candidates = [dict(id='c' + str(i), name='Probe ' + str(i), bin_id='r_B001', reference='r',
        imager='DP1', status='selected', bound_strand=strand,
        sites=[dict(reference='r', start=a, end=b, strand=strand)]) for i, (a, b) in enumerate(intervals)]
    return dict(status='PROVISIONAL_DRAFT_CROSS_PRIMARY_DISABLED', title='Synthetic saved result',
        policy=dict(bound_strand=strand, screen_cross_primary=False),
        references=[dict(id='r', length=3000, mapping=dict(verified=True, chromosome='chr2',
            start=10001, end=13000, orientation='+', assembly='GRCh38'))],
        bins=[dict(id='r_B001', reference='r', imager='DP1', start=101, end=2600, length=2500)],
        candidates=candidates, selected_ids=[c['id'] for c in candidates],
        counts=dict(selected=len(candidates), exact_sites=len(candidates), candidates=len(candidates), regions=1))


def feature(kind, ident, start, end, strand='+', parent=''):
    return dict(type=kind, feature_id=ident, label=ident, start=start, end=end,
                strand=strand, parent=parent, source='Synthetic annotations')


def add_annotation(result, features):
    result['references'][0].update(features=features, gene_annotation=dict(status='FETCHED', provider='Synthetic'))
    return result


class GeometryTests(unittest.TestCase):
    def test_partial_sectors_gap_edges_and_offsets(self):
        r = fixture()
        before = copy.deepcopy(r)
        report = coverage.build_report(r)
        region = report['regions'][0]
        self.assertEqual(region['sectors'], [dict(start=101,end=1100,count=1),
            dict(start=1101,end=2100,count=1), dict(start=2101,end=2600,count=1)])
        self.assertEqual(region['metrics']['largest_gap_bp'], 1363)
        self.assertEqual(region['metrics']['left_edge_gap_bp'], 0)
        self.assertEqual(region['metrics']['right_edge_gap_bp'], 65)
        self.assertEqual(region['metrics']['left_half'], 2)
        self.assertEqual(region['metrics']['right_half'], 1)
        self.assertAlmostEqual(region['metrics']['mean_offset_bp'], (118.5+1118.5+2517.5)/3-1350.5)
        self.assertEqual((region['genomic_start'],region['genomic_end']), (10101,12600))
        self.assertEqual(region['sites'][0]['start'], 101)
        self.assertEqual(r, before)

    def test_missing_annotation_is_unknown_but_fetched_empty_is_zero(self):
        unknown = coverage.build_report(fixture())
        self.assertEqual(unknown['annotation']['status'], 'UNKNOWN')
        self.assertIsNone(unknown['regions'][0]['metrics']['complete_intron_sites'])
        empty = coverage.build_report(add_annotation(fixture(), []))
        self.assertEqual(empty['annotation']['status'], 'AVAILABLE')
        self.assertEqual(empty['regions'][0]['metrics']['complete_intron_sites'], 0)
        self.assertIn('does not establish', empty['annotation']['note'])

    def test_empty_selection_has_unknown_mean_and_no_internal_gap(self):
        r = coverage.build_report(add_annotation(fixture(()), []))['regions'][0]
        self.assertEqual(r['metrics']['sites'], 0)
        self.assertIsNone(r['metrics']['mean_offset_bp'])
        self.assertIsNone(r['metrics']['largest_gap_bp'])
        self.assertEqual(r['metrics']['left_edge_gap_bp'], 2500)
        self.assertEqual(r['metrics']['right_edge_gap_bp'], 2500)
        self.assertEqual(sum(s['count'] for s in r['sectors']), 0)

    def test_introns_are_per_isoform_not_gene_union(self):
        r = fixture(((200,235),(400,435),(800,835)))
        f = [feature('gene','g',1,900), feature('transcript','t1',1,900,parent='g'),
             feature('exon','e1',1,100,parent='t1'), feature('exon','e2',700,900,parent='t1'),
             feature('transcript','t2',1,900,parent='g'), feature('exon','e3',1,300,parent='t2'),
             feature('exon','e4',350,500,parent='t2'), feature('exon','e5',700,900,parent='t2')]
        row = coverage.build_report(add_annotation(r,f))['regions'][0]
        self.assertEqual(row['metrics']['complete_intron_sites'], 2)
        self.assertEqual(len(row['genes'][0]['transcripts']), 2)
        self.assertEqual(row['genes'][0]['transcripts'][0]['exons'], [[700,900]])
        self.assertEqual(row['genes'][0]['exons'], [[101,300],[350,500],[700,900]])

    def test_minus_mapping_preserves_local_coordinates_and_flips_genomic_fallback(self):
        r = fixture(((200,235),), strand='-')
        r['references'][0]['mapping']['orientation'] = '-'
        # These genomic plus-strand features become reference-local minus-strand.
        f = []
        for kind,ident,a,b,parent in [('gene','g',1,900,''),('transcript','t',1,900,'g'),
                                     ('exon','e1',1,100,'t'),('exon','e2',700,900,'t')]:
            f.append(dict(type=kind,feature_id=ident,parent=parent,strand='+',
                          genomic_start=13000-b+1,genomic_end=13000-a+1))
        row = coverage.build_report(add_annotation(r,f))['regions'][0]
        self.assertEqual((row['genomic_start'],row['genomic_end']), (10401,12900))
        self.assertEqual(row['mapping_orientation'],'-')
        self.assertEqual((row['sites'][0]['start'],row['sites'][0]['strand']), (200,'-'))
        self.assertEqual(row['genes'][0]['strand'], '-')
        self.assertEqual(row['metrics']['complete_intron_sites'], 1)
        # Opposite local transcript strand is not complementary to these arms.
        for x in f:x['strand']='-'
        self.assertEqual(coverage.build_report(r)['regions'][0]['metrics']['complete_intron_sites'],0)

    def test_local_minus_annotations_not_flipped_twice(self):
        r=fixture(((200,235),),strand='-');r['references'][0]['mapping']['orientation']='-'
        f=[feature('gene','g',1,900,'-'),feature('transcript','t',1,900,'-',parent='g'),
           feature('exon','e1',1,100,'-',parent='t'),feature('exon','e2',700,900,'-',parent='t')]
        row=coverage.build_report(add_annotation(r,f))['regions'][0]
        self.assertEqual(row['genes'][0]['strand'],'-')
        self.assertEqual(row['metrics']['complete_intron_sites'],1)

    def test_missing_isoform_exons_or_unknown_strand_are_unknown_not_zero(self):
        cases=[
            [feature('gene','g',1,900)],
            [feature('gene','g',1,900),feature('transcript','t',1,900,parent='g')],
            [feature('gene','g',1,900,'?'),feature('transcript','t',1,900,'?',parent='g'),
             feature('exon','e1',1,100,'?',parent='t'),feature('exon','e2',700,900,'?',parent='t')]]
        for features in cases:
            row=coverage.build_report(add_annotation(fixture(((200,235),)),features))['regions'][0]
            self.assertEqual(row['annotation']['status'],'PARTIAL')
            self.assertIsNone(row['metrics']['complete_intron_sites'])
            self.assertEqual(len(row['genes']),1)
        # A known gene wholly outside the target is not a missing local model.
        row=coverage.build_report(add_annotation(fixture(),[feature('gene','g',1,50)]))['regions'][0]
        self.assertEqual(row['annotation']['status'],'AVAILABLE')
        self.assertEqual(row['metrics']['complete_intron_sites'],0)

    def test_duplicate_missing_and_outside_sites_are_rejected(self):
        for mutate in [lambda r:r['selected_ids'].append('missing'),
                       lambda r:r['selected_ids'].append('c0'),
                       lambda r:r['candidates'][0]['sites'][0].update(start=0),
                       lambda r:r['candidates'][0]['sites'][0].update(strand='-')]:
            r=fixture();mutate(r)
            with self.assertRaises(ValueError):coverage.build_report(r)

    def test_source_label_counts_sites_and_not_primaries(self):
        records=[(dict(start=1,end=10000),Path('/unused'),dict(sites=3,selected=2)) for _ in range(6)]
        self.assertEqual(coverage._experiment_label({},records,'20261026_x'),'Six 10 kb targets · 18 sites')

    def test_truncated_transcript_edges_are_uncertain_per_site_and_strand(self):
        r=fixture(((200,235),(800,835),(420,435)))
        r['references'][0]['length']=1000;r['references'][0]['mapping']['end']=11000
        r['bins'][0].update(start=1,end=1000,length=1000)
        features=[feature('gene','g',-200,1500),feature('transcript','t',-200,1500,parent='g'),
                  feature('exon','e1',100,150,parent='t'),feature('exon','e2',400,450,parent='t')]
        row=coverage.build_report(add_annotation(r,features))['regions'][0]
        self.assertIsNone(row['metrics']['complete_intron_sites'])
        self.assertEqual(row['metrics']['known_complete_intron_sites'],1)
        self.assertEqual(row['metrics']['uncertain_intron_sites'],1)
        self.assertEqual([s['intron_status'] for s in row['sites']],
                         ['COMPLETE_INTRON_MATCH','NO_COMPLETE_INTRON_MATCH','UNKNOWN'])
        # A truncated reverse-strand model is not uncertainty for Top+ arms.
        for f in features:f['strand']='-'
        row=coverage.build_report(r)['regions'][0]
        self.assertEqual(row['metrics']['complete_intron_sites'],0)
        self.assertEqual(row['metrics']['uncertain_intron_sites'],0)

    def test_known_other_isoform_match_resolves_truncated_edge_uncertainty(self):
        r=fixture(((200,235),(800,835)))
        r['references'][0]['length']=1000;r['references'][0]['mapping']['end']=11000
        r['bins'][0].update(start=1,end=1000,length=1000)
        features=[feature('gene','g',-200,1500),feature('transcript','t1',-200,1500,parent='g'),
                  feature('exon','e1',100,150,parent='t1'),feature('exon','e2',400,450,parent='t1'),
                  feature('transcript','t2',700,900,parent='g'),feature('exon','e3',700,750,parent='t2'),
                  feature('exon','e4',850,900,parent='t2')]
        row=coverage.build_report(add_annotation(r,features))['regions'][0]
        self.assertEqual(row['metrics']['complete_intron_sites'],2)
        self.assertEqual(row['metrics']['known_complete_intron_sites'],2)
        self.assertEqual(row['metrics']['uncertain_intron_sites'],0)

    def test_no_large_evidence_is_traversed(self):
        class Explosive:
            def __iter__(self):raise AssertionError('Large scientific evidence traversed')
        r=fixture();r['candidates'][0]['offtarget_hits']=Explosive();r['blast_runs']=Explosive()
        self.assertEqual(coverage.build_report(r)['regions'][0]['metrics']['sites'],3)


class SavedSourceTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.base=Path(self.temp.name)
        self.jobs=self.base/'design_runs';self.jobs.mkdir()
        self.patch=patch.object(studio,'JOBS',self.jobs);self.patch.start()
        self.ident='20261026_test';self.root=self.write_run(self.ident)

    def tearDown(self):
        self.patch.stop();self.temp.cleanup()

    def write_run(self, ident, result=None):
        root=self.jobs/ident;(root/'result').mkdir(parents=True)
        (root/'result/design_result.json').write_text(json.dumps(result or fixture()))
        (root/'analysis_status.json').write_text(json.dumps(dict(status='complete',analysis_complete=True,job_id=ident)))
        (root/'job_status.json').write_text(json.dumps(dict(status='complete',id=ident)))
        (self.jobs/(ident+'.zip')).write_bytes(b'committed archive fixture')
        return root

    def experiment(self):
        root=self.base/'experiments'/'six_targets';root.mkdir(parents=True)
        targets=[]
        for i in range(6):
            ident='run_'+str(i);r=fixture();r['bins'][0]['id']='region_'+str(i)
            for c in r['candidates']:c['bin_id']=r['bins'][0]['id']
            source=self.write_run(ident,r)
            targets.append(dict(label='Target '+str(i),studio_run_id=ident,source_run_dir=str(source),
                status='COMPLETE',analysis_complete=True,imager='DP1',chromosome='chr2',start=10101,end=12600))
        manifest=dict(analysis_complete=True,queue_status='COMPLETE',targets=targets)
        (root/'experiment_manifest.json').write_text(json.dumps(manifest))
        return root,manifest

    def test_completed_run_and_experiment_order(self):
        root,manifest=self.experiment()
        report=coverage.load_source(self.base,self.jobs,studio.completed_run,'experiment','six_targets')
        self.assertEqual([r['label'] for r in report['regions']],[t['label'] for t in manifest['targets']])
        sources=coverage.list_sources(self.base,self.jobs,studio.completed_run)['sources']
        self.assertEqual(sources[0]['kind'],'experiment')
        self.assertIn(self.ident,[s['id'] for s in sources])
        self.assertEqual(set(sources[0]),{'id','label','kind'})

    def test_incomplete_and_uncommitted_sources_rejected_before_result_read(self):
        (self.root/'analysis_status.json').write_text(json.dumps(dict(status='incomplete',analysis_complete=False,job_id=self.ident)))
        (self.root/'result/design_result.json').write_text('unreadable result should never be read')
        with self.assertRaisesRegex(ValueError,'incomplete'):
            coverage.load_source(self.base,self.jobs,studio.completed_run,'run',self.ident)
        self.assertEqual(coverage.list_sources(self.base,self.jobs,studio.completed_run)['sources'],[])
        (self.root/'analysis_status.json').write_text(json.dumps(dict(status='complete',analysis_complete=True,job_id=self.ident)))
        (self.jobs/(self.ident+'.zip')).unlink()
        with self.assertRaisesRegex(ValueError,'incomplete'):
            coverage.load_source(self.base,self.jobs,studio.completed_run,'run',self.ident)

    def test_experiment_requires_every_source_and_exact_provenance(self):
        root,m=self.experiment();m['targets'][5]['source_run_dir']=str(self.root)
        (root/'experiment_manifest.json').write_text(json.dumps(m))
        with self.assertRaisesRegex(ValueError,'provenance'):
            coverage.load_source(self.base,self.jobs,studio.completed_run,'experiment',root.name)
        m['targets'][5]['source_run_dir']=str(self.jobs/'run_5')
        (root/'experiment_manifest.json').write_text(json.dumps(m));(self.jobs/'run_5.zip').unlink()
        with self.assertRaisesRegex(ValueError,'incomplete'):
            coverage.load_source(self.base,self.jobs,studio.completed_run,'experiment',root.name)

    def test_independent_run_passes_never_become_cross_run_pass(self):
        root,_=self.experiment()
        for source in ['run_'+str(i) for i in range(6)]:
            p=self.jobs/source/'result/design_result.json';r=json.loads(p.read_text())
            r['policy']['screen_cross_primary']=True;r['cross_primary_screening']={'status':'PASS','enabled':True}
            p.write_text(json.dumps(r))
        report=coverage.load_source(self.base,self.jobs,studio.completed_run,'experiment',root.name)
        self.assertEqual(report['cross_primary_status'],'BETWEEN SEPARATE RUNS NOT ASSESSED')
        self.assertTrue(all(r['cross_primary_status']=='PASS' for r in report['regions']))
        self.assertTrue(any('without a combined PASS' in w for w in report['warnings']))
        _,small=coverage._experiment(self.base,self.jobs,root.name,studio.completed_run)
        self.assertTrue(all('result' not in metadata for _,_,metadata in small))

    def test_paths_symlinked_results_and_annotation_directories_rejected(self):
        for ident in ['../escape','/tmp/outside','a/b','a\\b','',self.ident+'%2f..']:
            with self.assertRaises(ValueError):coverage.load_source(self.base,self.jobs,studio.completed_run,'run',ident)
        (self.jobs/'alias').symlink_to(self.root,target_is_directory=True)
        with self.assertRaises(ValueError):coverage.load_source(self.base,self.jobs,studio.completed_run,'run','alias')
        outside=self.base/'outside';outside.mkdir()
        (self.root/'result/inputs').symlink_to(outside,target_is_directory=True)
        (self.root/'result/run_config.json').write_text(json.dumps(dict(configuration=dict(references=[]))))
        with self.assertRaisesRegex(ValueError,'symlinks'):
            coverage.load_source(self.base,self.jobs,studio.completed_run,'run',self.ident)

    def annotation_assets(self, features):
        out=self.root/'result';(out/'inputs').mkdir(exist_ok=True)
        ref=copy.deepcopy(fixture()['references'][0]);ref.update(features={'$annotation':'/references/0/features'},gene_annotation={'status':'FETCHED'})
        path=out/'inputs/annotations.json.gz';path.write_bytes(gzip.compress(json.dumps({'/references/0/features':features}).encode()))
        config=dict(configuration=dict(references=[ref]),assets={'inputs/annotations.json.gz':{'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}})
        (out/'run_config.json').write_text(json.dumps(config));return path

    def test_compressed_annotation_integrity_and_bounded_expansion(self):
        path=self.annotation_assets([])
        report=coverage.load_source(self.base,self.jobs,studio.completed_run,'run',self.ident)
        self.assertEqual(report['annotation']['status'],'AVAILABLE')
        path.write_bytes(gzip.compress(b'{}'))
        with self.assertRaisesRegex(ValueError,'checksum'):
            coverage.load_source(self.base,self.jobs,studio.completed_run,'run',self.ident)
        path=self.annotation_assets([feature('gene','g',1,100)])
        with patch.object(coverage,'MAX_ANNOTATION_BYTES',10):
            with self.assertRaisesRegex(ValueError,'expanded read limit'):
                coverage.load_source(self.base,self.jobs,studio.completed_run,'run',self.ident)

    def test_directory_and_report_limits(self):
        with patch.object(coverage,'MAX_DIRECTORY_ENTRIES',0):
            with self.assertRaisesRegex(ValueError,'Too many'):
                coverage.list_sources(self.base,self.jobs,studio.completed_run)
        with patch.object(coverage,'MAX_REPORT_BYTES',10):
            with self.assertRaisesRegex(ValueError,'response limit'):coverage.build_report(fixture())

    def request(self,path):
        # Exercise the actual GET dispatch without opening a server or changing app state.
        handler=object.__new__(studio.Handler);handler.path=path;handler.safe_host=lambda:True
        responses=[];handler.send=lambda data,code=200,**kw:responses.append((code,data,kw))
        with patch.object(studio,'BASE',self.base),patch.object(studio,'STATE',{}):handler.do_GET()
        return responses[0]

    def test_api_result_attaches_small_report_without_rewriting(self):
        path=self.root/'result/design_result.json';before=path.read_bytes()
        code,body,_=self.request('/api/result?id='+self.ident)
        self.assertEqual(code,200);self.assertEqual(body['coverage_report']['schema'],coverage.SCHEMA)
        self.assertEqual(path.read_bytes(),before)
        code,body,_=self.request('/api/coverage?kind=run&id='+self.ident)
        self.assertEqual(code,200);self.assertEqual(body['regions'][0]['metrics']['sites'],3)
        self.assertEqual(self.request('/api/coverage?kind=run&id=../escape')[0],400)

    def test_dynamic_viewer_download_keeps_historical_viewer(self):
        path=self.root/'result/viewer.html';path.write_text('old historical viewer')
        with patch.object(studio,'standalone',side_effect=lambda result:'new viewer '+result['coverage_report']['schema']):
            code,body,meta=self.request('/download?id='+self.ident+'&file=viewer.html')
        self.assertEqual(code,200);self.assertIn(coverage.SCHEMA,body);self.assertIn('viewer.html',meta['filename'])
        self.assertEqual(path.read_text(),'old historical viewer')

    def test_visual_failure_does_not_replace_scientific_completion(self):
        r=fixture()
        with patch.object(coverage,'build_report',side_effect=ValueError('annotation unavailable')):
            enriched=studio.result_with_coverage(r,self.root/'result')
        self.assertEqual(enriched['status'],r['status']);self.assertNotIn('coverage_report',enriched)
        self.assertIn('annotation unavailable',enriched['coverage_report_error']);self.assertNotIn('coverage_report_error',r)

    def test_offline_viewer_embeds_coverage_and_both_modules(self):
        enriched=studio.result_with_coverage(fixture(),self.root/'result')
        html=studio.standalone(enriched)
        boot=json.loads(re.search(r'window.STUDIO_BOOT=(.*?);</script>',html,re.S).group(1))
        self.assertTrue(boot['offline'])
        self.assertEqual(boot['demo']['coverage_report'],enriched['coverage_report'])
        self.assertNotIn('src="/web/probe_coverage.js"',html)
        self.assertNotIn('src="/web/coverage_panel.js"',html)
        for name in ('probe_coverage.js','coverage_panel.js'):
            self.assertIn((studio.BASE/'web'/name).read_text(),html)


@unittest.skipUnless((INSTALLED/'experiments'/EXPERIMENT/'experiment_manifest.json').is_file() and
                     (AUDIT/'coverage_review_summary.json').is_file(), 'Saved six-target audit fixture is not installed')
class SavedSixTargetParity(unittest.TestCase):
    def test_six_saved_run_parity_and_read_only_experiment(self):
        manifest=json.loads((INSTALLED/'experiments'/EXPERIMENT/'experiment_manifest.json').read_text())
        baseline=json.loads((AUDIT/'coverage_review_summary.json').read_text())['targets']
        tracks=json.loads((AUDIT/'transcript_tracks.json').read_text())['targets']
        paths=[Path(t['source_run_dir'])/'result/design_result.json' for t in manifest['targets']]
        before={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
        with patch.object(studio,'JOBS',INSTALLED/'design_runs'):
            report=coverage.load_source(INSTALLED,INSTALLED/'design_runs',studio.completed_run,'experiment',EXPERIMENT)
        self.assertEqual(report['cross_primary_status'],coverage.DISABLED)
        self.assertEqual(report['order_status'],'DRAFT')
        self.assertEqual([r['metrics']['sites'] for r in report['regions']],[103,80,89,80,98,144])
        self.assertEqual([r['metrics']['complete_intron_sites'] for r in report['regions']],[0,80,89,13,0,16])
        for row,expected,track in zip(report['regions'],baseline,tracks):
            self.assertEqual([s['count'] for s in row['sectors']],expected['sector_counts_1kb'])
            self.assertEqual(row['metrics']['left_half'],expected['left_half'])
            self.assertEqual(row['metrics']['right_half'],expected['right_half'])
            self.assertAlmostEqual(row['metrics']['mean_offset_bp'],expected['equal_weight_mean_genomic_offset_bp'])
            self.assertEqual(row['metrics']['largest_gap_bp'],expected['largest_internal_arm_gap_bp'])
            self.assertEqual(sum(len(g['transcripts']) for g in row['genes']),sum(len(g['transcripts']) for g in track['genes']))
            self.assertEqual([[s['start']+7820000,s['end']+7820000] for s in row['sites']],sorted(expected['selected_site_intervals']))
            self.assertTrue(all(row['start']<=s['start']<=s['end']<=row['end'] for s in row['sites']))
        self.assertEqual(before,{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})


if __name__=='__main__':unittest.main()
