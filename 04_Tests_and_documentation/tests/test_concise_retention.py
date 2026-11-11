import copy
import gzip
import json
import sys
import tempfile
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from compact_outputs import concise_saved_result, finalize_outputs, write_input_manifest
from test_compact_outputs import fixture_config,fixture_result,fixture_bins


class Unreadable:
    def __iter__(self):raise AssertionError('Omitted detailed data traversed')


class ConciseRetentionTests(unittest.TestCase):
    def fixture(self):
        config=fixture_config();r=fixture_result(config)
        r['candidates'][0]['fold'].update(base_unpaired=[.1,.9],structure='..',mfe_kcal_mol=-2.0)
        r['candidates'][0]['self_pair'].update(intermolecular_pairs=[[1,2],[3,4]],structure='(())')
        r['bins']=[dict(id='b',metrics={'probes':1,'site_centers':[10,20]})]
        r['candidate_summary']={'total':600,'statuses':{'selected':1,'rejected_by_QC':599}}
        r['counts']={'candidates':600,'selected':1,'QC_rejected':599}
        r['imagers']=[dict(id='DP3',sequence='ACGT',docking='ACGT')]
        r['blast_runs']=[dict(source='human',raw_alignments=20,retained_alignments=10,output_sha256='original',
            output_sha256_scope='ORDERED_BATCH_SHA256_AND_LENGTH_MANIFEST',command=['blastn','-word_size','7'],
            database_identity={'metadata_sha256':'dbhash'},search_manifest='blast_attempts/source/search_manifest.json',batches=[dict(batch=1,status='COMPLETE',output_sha256='abc',command=['blastn','1']),dict(batch=2,status='COMPLETE',output_sha256='def',command=['blastn','2'])])]
        return r

    def test_preserves_actual_counts_scores_sequences_and_sites(self):
        r=self.fixture();old=copy.deepcopy(r);got=concise_saved_result(r)
        for key in ('candidate_summary','counts','selected_ids','policy','imagers'):
            self.assertEqual(got[key],r[key])
        self.assertEqual(got['candidates'][0]['specificity'],r['candidates'][0]['specificity'])
        for key in ('arm','full','sites'):
            self.assertEqual(got['candidates'][0][key],r['candidates'][0][key])
        self.assertEqual(got['candidates'][0]['fold']['mfe_kcal_mol'],-2.0)
        self.assertEqual(got['candidates'][0]['self_pair']['intermolecular_pairs'],2)
        self.assertNotIn('structure',got['candidates'][0]['self_pair']);self.assertNotIn('base_unpaired',got['candidates'][0]['fold'])
        self.assertNotIn('site_centers',got['bins'][0]['metrics']);self.assertEqual(r,old)

    def test_omitted_arrays_and_hit_payloads_not_read(self):
        r=self.fixture();r['candidates'][0]['fold']['base_unpaired']=Unreadable();r['candidates'][0]['offtarget_hits']=Unreadable()
        got=concise_saved_result(r);self.assertNotIn('offtarget_hits',got['candidates'][0])

    def test_batch_summary_retains_all_receipt_digest_and_scientific_options(self):
        r=self.fixture();got=concise_saved_result(r)['blast_runs'][0]
        self.assertNotIn('batches',got);self.assertIsNone(got['search_manifest']);self.assertEqual(got['command'],r['blast_runs'][0]['command'])
        self.assertEqual(got['output_sha256'],'original');self.assertEqual(got['raw_alignments'],20)
        self.assertEqual(got['batch_receipt_summary']['count'],2)
        self.assertEqual(got['batch_receipt_summary']['statuses'],{'COMPLETE':2})
        r['blast_runs'][0]['batches'][1]['output_sha256']='changed'
        self.assertNotEqual(got['batch_receipt_summary']['ordered_receipts_sha256'],concise_saved_result(r)['blast_runs'][0]['batch_receipt_summary']['ordered_receipts_sha256'])

    def test_idempotent_conversion_keeps_original_accounting(self):
        first=concise_saved_result(self.fixture());self.assertEqual(concise_saved_result(first),first)

    def test_standard_export_omits_detail_full_audit_preserves_debug_detail(self):
        for mode in ('standard','full_audit'):
            with self.subTest(mode=mode),tempfile.TemporaryDirectory() as tmp:
                config=fixture_config(mode);r=fixture_result(config)
                r['candidates'][0]['fold'].update(base_unpaired=[.1,.9],structure='..')
                r['candidates'][0]['self_pair'].update(intermolecular_pairs=[[1,2],[3,4]],structure='(())')
                root=Path(tmp);manifest=write_input_manifest(config,fixture_bins(),root)
                got=finalize_outputs(r,manifest,root)
                if mode=='standard':
                    self.assertNotIn('base_unpaired',got['candidates'][0]['fold'])
                    self.assertNotIn('structure',got['candidates'][0]['self_pair'])
                    self.assertEqual(got['candidates'][0]['self_pair']['intermolecular_pairs'],2)
                else:
                    with gzip.open(root/'debug/all_candidates.jsonl.gz','rt') as f:detail=json.loads(next(f))
                    self.assertEqual(detail['fold']['base_unpaired'],[.1,.9]);self.assertEqual(detail['self_pair']['structure'],'(())')

if __name__=='__main__':unittest.main()
