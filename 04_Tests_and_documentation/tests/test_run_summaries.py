import csv
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from run_summaries import human_readable_run_name, summarize_completed_result, write_run_summary


def fixture():
    return dict(title='ERRFI1 / DP1 — 20 kb', status='PROVISIONAL_DRAFT_CROSS_PRIMARY_DISABLED',
        policy={'screen_cross_primary': False}, generated_utc='2026-10-17T16:00:00+00:00',
        bins=[dict(id='r_B001', reference='r', start=366425, end=386424, length=20000,
            imager='DP3', chromosome='chr1', genomic_start=8186425, genomic_end=8206424,
            metrics={'probes': 1, 'exact_sites': 1, 'site_centers': [366450], 'largest_gap_bp': 1200})],
        assignments={'stale_B002': 'DP11'}, selected_ids=['C1'], counts={'selected': 1, 'candidates': 600, 'QC_rejected': 334},
        candidates=[dict(id='C1',name='DP3_001',bin_id='r_B001',imager='DP3',status='selected',
            tm=46.5,gc=50.0,arm_length=38,bound_strand='+',sites=[{'start':366425}],reject_reasons=[],
            offtarget_hits=Explosive(),specificity=dict(high_tm_sites=0,moderate_tm_sites=2,unknown_long_sites=0,strong_rna_sites=0,rna_risk_sites=0),
            fold=dict(status='RUN_AQUEOUS_MODEL',mean_arm_unpaired=.625,mfe_kcal_mol=-2.1,base_unpaired=Explosive()),
            self_pair=dict(status='RUN_AQUEOUS_MODEL',blocked=False,max_exact_stem=8,mfe_gain_vs_separate=-8.7,intermolecular_pairs=[[1,2],[3,4]]))],
        checks=[dict(name='Cross-primary interactions between regions',status='NOT RUN - USER DISABLED')])


class Explosive:
    def __iter__(self):
        raise AssertionError('Large evidence collection was traversed')


class RunSummaryTests(unittest.TestCase):
    def test_actual_bins_override_stale_title_and_assignments(self):
        self.assertEqual(human_readable_run_name(fixture()), 'DP3 · chr1:8186425-8206424')

    def test_configuration_names_only_current_regions(self):
        config=dict(title='ERRFI1 / DP1', bin_width=20000, references=[dict(id='r',length=400000,
            mapping=dict(chromosome='chr1',start=7820001,end=8220000))], bounds={'r':[366425,386424]},
            assignments={'r_B001':'DP3','r_B002':'DP11'})
        self.assertEqual(human_readable_run_name(config), 'DP3 · chr1:8186425-8206424')

    def test_reverse_mapped_config_uses_physical_coordinates(self):
        c=dict(bin_width=100,references=[dict(id='r',length=200,mapping={'chromosome':'chr5','start':1001,'end':1200,'orientation':'-'})],bounds={'r':[1,100]},assignments={'r_B001':'DP11'})
        self.assertEqual(human_readable_run_name(c),'DP11 · chr5:1101-1200')

    def test_user_name_supported_without_using_old_title(self):
        r=fixture();r['run_name']='  My   DP3\n design  '
        self.assertEqual(human_readable_run_name(r),'My DP3 design')

    def test_large_evidence_never_read_and_scores_exact(self):
        s=summarize_completed_result(fixture(),'run1'); row=s['selected_primary_scores'][0]
        self.assertEqual(row['moderate_tm_offtarget_sites'],2)
        self.assertEqual(row['arm_accessibility'],.625)
        self.assertEqual(row['self_dimer_intermolecular_pair_count'],2)
        self.assertEqual(row['self_dimer_gain_kcal_mol'],-8.7)
        self.assertNotIn('site_centers',s['regions'][0]['metrics'])
        self.assertNotIn('base_unpaired',json.dumps(s));self.assertNotIn('offtarget_hits',json.dumps(s))

    def test_disabled_within_and_between_never_pass(self):
        r=fixture();r['cross_primary_screening']={'status':'PASS','enabled':True,'evaluated':5}
        cross=summarize_completed_result(r,'r')['cross_primary_screening']
        self.assertEqual(cross['evaluated'],0)
        for key in ('status','within_region_status','between_region_status'):
            self.assertEqual(cross[key],'NOT RUN - USER DISABLED')

    def test_missing_measurement_not_zero_or_pass(self):
        r=fixture();r['candidates'][0]['fold']={};r['candidates'][0]['self_pair']={};r['candidates'][0]['specificity']={}
        row=summarize_completed_result(r,'r')['selected_primary_scores'][0]
        self.assertIsNone(row['high_tm_offtarget_sites']);self.assertIsNone(row['arm_accessibility'])
        self.assertEqual(row['individual_folding_status'],'NOT RECORDED');self.assertEqual(row['self_dimer_status'],'NOT RECORDED')

    def test_incomplete_results_refused(self):
        for status in ['STOPPED - INCOMPLETE','FAILED','CANCELLED']:
            r=fixture();r['status']=status
            with self.assertRaisesRegex(ValueError,'incomplete'):summarize_completed_result(r,'r')
        r=fixture();r['analysis_complete']=False
        with self.assertRaisesRegex(ValueError,'incomplete'):summarize_completed_result(r,'r')

    def test_legacy_status_not_promoted_to_all_checks_passed(self):
        s=summarize_completed_result(fixture(),'run1')
        self.assertEqual(s['status'],'complete');self.assertEqual(s['analysis_status'],'PROVISIONAL_DRAFT_CROSS_PRIMARY_DISABLED')
        self.assertEqual(s['order_status'],'DRAFT');self.assertNotIn('score',s)

    def test_writes_small_scores_and_name(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp);s=write_run_summary(fixture(),out,'r')
            self.assertLess((out/'run_summary.json').stat().st_size,10000)
            self.assertEqual(json.loads((out/'run_summary.json').read_text()),s)
            self.assertEqual((out/'RUN_NAME.txt').read_text(),'DP3 · chr1:8186425-8206424\n')
            rows=list(csv.DictReader(io.StringIO((out/'primary_scores.csv').read_text())))
            self.assertEqual(rows[0]['cross_primary_status'],'NOT RUN - USER DISABLED')
            self.assertEqual(rows[0]['self_dimer_intermolecular_pair_count'],'2')


if __name__=='__main__':unittest.main()
