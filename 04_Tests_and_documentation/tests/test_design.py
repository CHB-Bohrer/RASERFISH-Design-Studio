import copy
import io
import json
import random
import shutil
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from sequence_tools import *
from design_engine import *
import annotation_sources as ann

BASE=Path(__file__).resolve().parents[1]
IMS=validate_imagers(json.loads((BASE/'imagers.json').read_text()))
REF={'id':'X','sequence':'A'*100,'length':100,'mapping':None,'features':[]}
BINS=[{'id':'X_B001','reference':'X','number':1,'start':1,'end':100,'length':100,'center':50.5,'imager':'DP7'}]

def cand(name='A',start=5,end=40,tm=47,bin_id='X_B001'):
    arm='ACGTTGCA'*4+'ACGT';full=IMS[0]['docking']+'TT'+arm+'TT'+IMS[0]['docking']
    return dict(id=name,bin_id=bin_id,reference='X',imager='DP7',arm=arm,full=full,prefix=IMS[0]['docking']+'TT',suffix='TT'+IMS[0]['docking'],
      start=start,end=end,sites=[dict(reference='X',start=start,end=end,extended_start=start,extended_end=end,strand='-')],
      tm=tm,gc=50,fold={'mean_arm_unpaired':.8},specificity={'risk_penalty':0},reject_reasons=[],status='candidate')

def hit(**kwargs):
    d=dict(reference='chr7',source='human',start=100,end=139,bound_strand='-',qseq='GCGCATGC'*5,sseq='GCGCATGC'*5,
      query_start=1,query_end=40,query_coverage=100.,identity=100.,alignment_length=40,mismatches=0,gap_opens=0,
      rna=False,allowed=False,full_primary=False)
    d.update(kwargs);return d

class SequenceTests(unittest.TestCase):
    def test_inventory_handles(self):
        by={i['id']:i for i in IMS};self.assertEqual(by['DP7']['docking'],'TATTCGTATCGG');self.assertEqual(by['DP3']['docking'],'CCGCGATATTAA')
    def test_inventory_variable_length(self):
        by={i['id']:i for i in IMS};self.assertEqual(len(by['DP2']['sequence']),13);self.assertEqual(len(by['DP5']['sequence']),13)
    def test_duplicate_imager(self):
        with self.assertRaises(ValueError):validate_imagers([IMS[0],dict(IMS[0],id='NEW')])
    def test_bad_sequence_not_cleaned_away(self):
        with self.assertRaises(ValueError):normalize('ACGT 123 ACTG')
    def test_iupac_coordinates_retained(self):
        self.assertEqual(parse_fasta('>X\nacgtnry')[0]['length'],7)
    def test_unique_ids(self):
        with self.assertRaises(ValueError):parse_fasta('>X\nACGT\n>X\nTTTT')
    def test_ref_aliases(self):
        self.assertEqual(canonical_ref('chr7'),canonical_ref('ref|NC_000007.14|'))
        self.assertNotEqual(canonical_ref('chr7'),canonical_ref('NC_000007.13'))
    def test_alt_not_collapsed(self):
        self.assertEqual(canonical_ref('ref|NT_123.1|'),'NT_123.1')
    def test_bins_equal_adjacent(self):
        r=dict(REF,length=307,sequence='A'*307);bs,ts=make_bins([r],100)
        self.assertEqual([b['length'] for b in bs],[100]*3);self.assertEqual([b['start'] for b in bs],[1,101,201]);self.assertEqual(ts[0]['length'],7)
    def test_no_short_final_bin(self):
        bs,ts=make_bins([dict(REF,length=199)],100);self.assertEqual(len(bs),1);self.assertEqual(ts[0]['length'],99)
    def test_exact_division_required(self):
        with self.assertRaises(ValueError):make_bins([dict(REF,length=307)],100,'strict')
    def test_bounds(self):
        bs,ts=make_bins([dict(REF,length=350)],100,'strict',bounds={'X':[51,350]});self.assertEqual(bs[0]['start'],51);self.assertEqual(bs[-1]['end'],350)
    def test_features_not_lost_at_boundary(self):
        f={'start':80,'end':120,'label':'gene'}
        a=clip_features([f],1,100)[0];b=clip_features([f],101,200)[0]
        self.assertEqual((a['view_start'],a['view_end']),(80,100));self.assertEqual((b['view_start'],b['view_end']),(101,120));self.assertEqual(f['end'],120)
    def test_reverse_mapping(self):
        r=dict(REF,mapping=dict(start=1001,end=1100,chromosome='chr7',orientation='-',assembly='GRCh38'))
        g=genomic_interval(r,10,20,'-');self.assertEqual((g['start'],g['end'],g['strand']),(1081,1091,'+'))
    def test_genbank_features(self):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature,FeatureLocation,CompoundLocation
        rec=SeqRecord(Seq('A'*100),id='test',annotations={'molecule_type':'DNA'})
        rec.features=[SeqFeature(CompoundLocation([FeatureLocation(3,8),FeatureLocation(40,60)]),type='CDS',qualifiers={'gene':['demo']})]
        stream=io.StringIO();SeqIO.write(rec,stream,'genbank');r=parse_annotated(stream.getvalue().encode(),'x.gb')[0]
        self.assertEqual([(f['start'],f['end']) for f in r['features']],[(4,8),(41,60)])
    def test_coverage_union(self):
        b={'reference':'X','start':1,'end':100,'length':100,'center':50.5};a=cand(start=5,end=40);d=cand('B',35,70)
        m=region_metrics(b,[a,d]);self.assertEqual(m['covered_bases'],66);self.assertEqual(m['probes'],2)

class ThermoSpecificityTests(unittest.TestCase):
    def setUp(self):self.p=Policy(mode='local_preview')
    def test_perfect_offtarget_Tm_equals_same_arm_model(self):
        s='GCGCATGC'*5;h=enrich_hit(hit(),cand(),self.p)
        self.assertAlmostEqual(h['offtarget_tm'],calculate_tm(s,PROFILES[self.p.profile])['adjusted_screen_C'],places=5)
    def test_mismatched_hit_not_on_target_Tm(self):
        q='GCGCATGC'*5;s=q[:15]+'A'+q[16:]
        a=duplex_estimate(q,q,'oligominer1_2nM')[0];b=duplex_estimate(q,s,'oligominer1_2nM')[0]
        if b is not None:self.assertNotEqual(a,b)
    def test_gapped_is_unknown(self):self.assertIsNone(duplex_estimate('AAC-GTT','AACC GTT'.replace(' ',''),'oligominer1_2nM')[0])
    def test_rna_not_DNA_tm(self):self.assertEqual(duplex_estimate('A'*40,'A'*40,'oligominer1_2nM',True),(None,'RNA_SEQUENCE_ONLY'))
    def test_coordinate_intended_exemption(self):
        r=dict(REF,mapping=dict(start=100,end=199,chromosome='chr7',orientation='+',assembly='GRCh38'))
        self.assertTrue(expected_hit(hit(),cand(),{'X_B001':BINS[0]},{'X':r},[]))
    def test_other_human_location_not_exempt(self):
        r=dict(REF,mapping=dict(start=100,end=199,chromosome='chr7',orientation='+',assembly='GRCh38'))
        self.assertFalse(expected_hit(hit(start=200,end=239),cand(),{'X_B001':BINS[0]},{'X':r},[]))
    def test_wrong_strand_not_exempt(self):
        r=dict(REF,mapping=dict(start=100,end=199,chromosome='chr7',orientation='+',assembly='GRCh38'))
        self.assertFalse(expected_hit(hit(bound_strand='+'),cand(),{'X_B001':BINS[0]},{'X':r},[]))
    def test_RNA_not_exempt_by_DNA_allow(self):
        self.assertFalse(expected_hit(hit(rna=True),cand(),{'X_B001':BINS[0]},{'X':REF},[dict(reference='chr7',start=1,end=500,strand='-')]))
    def test_duplicate_hsps_merge(self):
        self.assertEqual(len(cluster_hits([hit(),hit(start=101,end=140,full_primary=True)])),1)
    def test_disjoint_sites_separate(self):self.assertEqual(len(cluster_hits([hit(),hit(start=500,end=540)])),2)
    def test_distinct_chromosomes_separate(self):self.assertEqual(len(cluster_hits([hit(),hit(reference='chr8')])),2)
    def test_high_hit_rejects(self):
        c=cand();c['offtarget_hits']=[enrich_hit(hit(),c,self.p)];summarize_specificity(c,self.p);self.assertTrue(c['reject_reasons'])
    def test_allowed_high_hit_not_rejected(self):
        c=cand();c['offtarget_hits']=[enrich_hit(hit(allowed=True),c,self.p)];summarize_specificity(c,self.p);self.assertEqual(c['specificity']['high_tm_sites'],0)
    def test_report_only_override_retains_evidence(self):
        c=cand();c['offtarget_hits']=[enrich_hit(hit(),c,self.p)];p=Policy(mode='local_preview',human_policy='report_only');summarize_specificity(c,p)
        self.assertEqual(c['reject_reasons'],[]);self.assertEqual(c['specificity']['ignored_human_review_hsps'],1)
    def test_report_only_does_not_hide_local_hit(self):
        c=cand();c['offtarget_hits']=[enrich_hit(hit(source='targets'),c,self.p)];summarize_specificity(c,Policy(mode='local_preview',human_policy='report_only'));self.assertTrue(c['reject_reasons'])
    def test_many_moderate_sites_taken_seriously(self):
        p=Policy(mode='local_preview',max_moderate_offtarget_sites=2,max_high_offtarget_sites=10)
        c=cand();c['offtarget_hits']=[dict(hit(start=100*i,end=100*i+39),review=True,high_tm=False,moderate_tm=True,unknown_long=False,strong_rna=False) for i in (1,2,3)]
        summarize_specificity(c,p);self.assertTrue(c['reject_reasons'])
    def test_no_missing_Tm_as_zero(self):
        c=cand();h=enrich_hit(hit(qseq='GC-GCATGC'*5,sseq='GCTGCATGC'*5),c,self.p);self.assertIsNone(h['offtarget_tm']);self.assertTrue(h['unknown_long'])
    def test_allowlist_rejects_RNA_exemption(self):
        with self.assertRaises(ValueError):validate_allow([dict(reference='rna',start=1,end=100,strand='-',source='rna_human')],BINS)

class SpacingSelectionTests(unittest.TestCase):
    def test_five_bp_gap(self):
        p=Policy();a=cand(start=5,end=40);b=cand('B',46,81);self.assertFalse(footprints_conflict(a['sites'],b['sites'],p));b['sites'][0]['extended_start']=45;self.assertTrue(footprints_conflict(a['sites'],b['sites'],p))
    def test_flank_extension(self):
        seq='CCCCTT'+'GCGC'+'TTAAAA';a,b=footprint(seq,7,10,'AAATT','TTGG','-');self.assertEqual((a,b),(5,12))
    def test_reverse_flank_orientation(self):
        seq='CC'+'ACGT'+'GG';a,b=footprint(seq,3,6,'CC','GG','+');self.assertEqual((a,b),(1,8))
    def test_footprints_across_bins(self):
        a=cand(start=65,end=99);b=cand('B',103,140,bin_id='X_B002');self.assertTrue(footprints_conflict(a['sites'],b['sites'],Policy()))
    def test_local_seed_search_matches_bruteforce(self):
        rand=random.Random(5);seq=''.join(rand.choice('ACGT') for _ in range(150));q=seq[44:84];q=q[:12]+('C' if q[12]!='C' else 'A')+q[13:]
        r={'id':'X','sequence':seq};ix=LocalIndex([r]);found=ix.full_matches(q,r,3);brute=[(i,sum(a!=b for a,b in zip(q,seq[i:i+len(q)]))) for i in range(len(seq)-len(q)+1)];self.assertEqual(found,[h for h in brute if h[1]<=3])
    def test_repeat_cross_bin_rejected(self):
        self.assertFalse(same_bin_site({'reference':'X','start':105,'end':140,'strand':'-'},BINS[0],Policy()))
    def test_selection_excludes_QC_failures(self):
        class Fake:
            def pair(self,a,b):return {'blocked':False}
        a=cand('A');a['reject_reasons']=['off-target'];b=cand('B',51,86);out=select_candidates([a,b],copy.deepcopy(BINS),Policy(mode='local_preview',target_per_bin=2,minimum_per_bin=1),Fake(),lambda *a,**k:None)
        self.assertEqual([x['id'] for x in out],['B'])
    def test_pair_affects_selection_not_just_label(self):
        class Fake:
            def pair(self,a,b):return {'blocked':a['id']!=b['id']}
        a=cand('A');b=cand('B',51,86);out=select_candidates([a,b],copy.deepcopy(BINS),Policy(mode='local_preview',screen_cross_primary=True,target_per_bin=2,minimum_per_bin=1),Fake(),lambda *a,**k:None)
        self.assertEqual(len(out),1);self.assertTrue(any(c['status']=='not_selected_pair' for c in (a,b)))
    def test_no_silent_threshold_relaxation(self):
        class Fake:
            def pair(self,a,b):return {'blocked':False}
        bins=copy.deepcopy(BINS);a=cand('A');a['reject_reasons']=['folding'];out=select_candidates([a],bins,Policy(mode='local_preview'),Fake(),lambda *a,**k:None)
        self.assertEqual(out,[]);self.assertTrue(bins[0]['metrics']['below_minimum']);self.assertEqual(bins[0]['metrics']['shortfall_to_target'],Policy().target_per_bin)
    def test_full_mode_does_not_fallback(self):
        with patch.dict(sys.modules,{'RNA':None}):
            with self.assertRaises(RuntimeError):Folding(Policy(mode='full'))
    def test_invalid_thermo_bounds(self):
        with self.assertRaises(ValueError):Policy.from_dict({'tm_min':55,'tm_max':42})
    def test_input_sequence_preserved_in_demo(self):
        r=json.loads((BASE/'examples/local_design/design_result.json').read_text());self.assertEqual(r['reference_hashes_before'],r['reference_hashes_after']);self.assertEqual(r['status'],'LOCAL_ONLY_DESIGN_PREVIEW')
    def test_demo_has_optimized_new_arms(self):
        r=json.loads((BASE/'examples/local_design/design_result.json').read_text());by={c['id']:c for c in r['candidates']};self.assertTrue(all(by[i]['status']=='selected' for i in r['selected_ids']));self.assertGreater(r['counts']['QC_rejected'],0)
    def test_demo_final_global_spacing_all_copies(self):
        r=json.loads((BASE/'examples/local_design/design_result.json').read_text());p=Policy(**r['policy']);ids=set(r['selected_ids']);sel=[c for c in r['candidates'] if c['id'] in ids]
        for i,a in enumerate(sel):
            for b in sel[i+1:]:self.assertFalse(footprints_conflict(a['sites'],b['sites'],p))
    def test_demo_no_reused_imager(self):
        r=json.loads((BASE/'examples/local_design/design_result.json').read_text());self.assertEqual(len(r['bins']),len(set(b['imager'] for b in r['bins'])))
    def test_demo_Tm_gate(self):
        r=json.loads((BASE/'examples/local_design/design_result.json').read_text());sel=set(r['selected_ids']);self.assertTrue(all(r['policy']['tm_min']<c['tm']<r['policy']['tm_max'] for c in r['candidates'] if c['id'] in sel))

class AnnotationTests(unittest.TestCase):
    def exp(self,cell='HCT116',assembly='GRCh38',status='released',modified=False):
        return dict(status=status,accession='ENCSRTEST',biosample_ontology={'term_name':cell},target={'label':'H3K27ac'},
          assay_term_name='ChIP-seq',biosample_summary='Homo sapiens '+cell+(' genetically modified' if modified else ''),replicates=[],
          files=[dict(status='released',assembly=assembly,file_format='bed',file_format_type='narrowPeak',output_type='replicated peaks',
                    accession='ENCFFTEST',href='/files/ENCFFTEST/@@download/ENCFFTEST.bed.gz',biological_replicates=[1,2])])
    def test_cell_aliases(self):self.assertEqual(ann.cell_name('IMR90'),'IMR-90');self.assertEqual(ann.cell_name('hct-116'),'HCT116')
    def test_genomic_retrieval_count(self):
        ref=ann.get_genomic_reference('chr7:101-110',lambda url:{'dna':'ACGTACGTAC'});self.assertEqual(ref['mapping']['start'],101);self.assertEqual(ref['length'],10)
    def test_sequence_length_mismatch(self):
        with self.assertRaises(RuntimeError):ann.get_genomic_reference('chr7:101-110',lambda url:{'dna':'AAA'})
    def test_no_fake_gene_mapping(self):
        with self.assertRaises(ValueError):ann.fetch_genes(REF,lambda url:[])
    def test_genes_keep_true_span(self):
        ref=dict(REF,mapping=dict(verified=True,start=1001,end=1100,chromosome='chr7'))
        d=ann.fetch_genes(ref,lambda u:[dict(start=990,end=1050,strand=1,feature_type='gene',id='GENE',assembly_name='GRCh38')]);self.assertEqual((d['features'][0]['start'],d['features'][0]['end']),(-10,50))
    def test_genes_assembly_mismatch(self):
        ref=dict(REF,mapping=dict(verified=True,start=1001,end=1100,chromosome='chr7'))
        with self.assertRaises(RuntimeError):ann.fetch_genes(ref,lambda u:[dict(start=1001,end=1010,assembly_name='GRCh37')])
    def test_requested_cell_first(self):
        got=ann.discover_h3k27ac('HCT116',True,lambda url:{'@graph':[self.exp()]});self.assertFalse(got['is_proxy']);self.assertEqual(got['used_cell'],'HCT116')
    def test_fallback_to_HCT116(self):
        def get(url):return {'@graph':[self.exp()]} if 'HCT116' in url else {'@graph':[]}
        got=ann.discover_h3k27ac('MyCells',True,get);self.assertTrue(got['is_proxy']);self.assertEqual(got['used_cell'],'HCT116')
    def test_fallback_to_IMR90(self):
        def get(url):return {'@graph':[self.exp('IMR-90')]} if 'IMR-90' in url else {'@graph':[]}
        got=ann.discover_h3k27ac('MyCells',True,get);self.assertEqual(got['used_cell'],'IMR-90');self.assertEqual(len(got['attempts']),3)
    def test_network_failure_not_no_data(self):
        with self.assertRaises(RuntimeError):ann.discover_h3k27ac('MyCells',True,lambda url:(_ for _ in ()).throw(RuntimeError('network')))
    def test_hg19_not_overlayed(self):
        got=ann.discover_h3k27ac('HCT116',False,lambda u:{'@graph':[self.exp(assembly='hg19')]});self.assertEqual(got['choices'],[])
    def test_modified_or_treated_excluded(self):
        self.assertFalse(ann.eligible_experiment(self.exp(modified=True),'HCT116'));e=self.exp();e['treatments']=['drug'];self.assertFalse(ann.eligible_experiment(e,'HCT116'))
    def test_bed_coordinates(self):
        ref=dict(REF,mapping=dict(start=101,end=200,chromosome='chr7'))
        got=ann.parse_bed_track('chr7\t109\t120\tpeak\t0\t.\t12',ref,{'kind':'peaks','cell':'HCT116'});self.assertEqual((got[0]['start'],got[0]['end']),(10,20));self.assertEqual(got[0]['signal'],12.)
    def test_local_bed_reference_id(self):
        got=ann.parse_bed_track('X\t0\t5\tfeature\nY\t0\t50\tother',REF,{'kind':'annotation'},True);self.assertEqual(len(got),1);self.assertEqual((got[0]['start'],got[0]['end']),(1,5))
    def test_wrong_BED_chromosome(self):
        ref=dict(REF,mapping=dict(start=101,end=200,chromosome='chr7'));self.assertEqual(ann.parse_bed_track('chr8\t109\t120',ref,{}),[])

class RealEngineTests(unittest.TestCase):
    @unittest.skipUnless(__import__('importlib').util.find_spec('RNA'),'ViennaRNA not installed here; run this test in the existing local .qc_env')
    def test_real_full_length_DNA(self):
        f=Folding(Policy(mode='full'));a=cand();result=f.fold(a);self.assertEqual(len(result['base_unpaired']),len(a['full']));self.assertGreaterEqual(result['mean_arm_unpaired'],0);self.assertLessEqual(result['mean_arm_unpaired'],1);self.assertEqual(f.pair(a,a)['status'],'RUN_AQUEOUS_MODEL')

if __name__=='__main__':unittest.main()
