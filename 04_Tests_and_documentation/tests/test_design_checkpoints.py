"""Interrupted work can be resumed only from authenticated committed batches."""
import copy
import json
import os
import random
import signal
import sqlite3
import subprocess
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
sys.path.insert(0,str(Path(__file__).resolve().parent))
import design_engine as engine
from specificity_store import HitStore
from design_checkpoints import read_checkpoint_summary
from test_blast_batches import fixture_candidates,fixture_output,expected_batch_environment
from test_design import IMS,hit,cand


def quiet(*args,**kwargs):pass


class TransactionalBatchTests(unittest.TestCase):
    def run_batches(self,folder,store_path,*,resume=False,fail_at=None,invalid_at=None,cleanup_failure=False):
        folder.mkdir()
        candidates,refs,bins=fixture_candidates(3)
        store=HitStore(store_path,resume=resume)
        for c in candidates:c['offtarget_hits']=store.hits(c['id'])
        calls=[]
        def execute(command,outlog,cancel=None,resource_guard=None,*,env_overrides=None):
            self.assertEqual(env_overrides,expected_batch_environment(command))
            calls.append(command)
            data=fixture_output(command)
            if len(calls)==invalid_at:data=data.splitlines(keepends=True)[0]+b'not a valid alignment\n'
            Path(command[command.index('-out')+1]).write_bytes(data)
            outlog.write_text('Synthetic checkpoint fixture\n')
            if len(calls)==fail_at:raise RuntimeError('Forced failed subprocess')
        checkpoint=types.SimpleNamespace(update_progress=quiet)
        try:
            with patch.object(engine,'BLAST_BATCH_CANDIDATES',1),patch.object(engine.shutil,'which',return_value='/fixture/blastn'),patch.object(engine.subprocess,'check_output',return_value='fixture'),patch.object(engine,'run_process',side_effect=execute):
                if cleanup_failure:
                    with patch.object(engine,'archive_transient',side_effect=RuntimeError('Forced cleanup interruption')):
                        result=engine.blast_candidates(candidates,refs,bins,[],{'human':'human_fixture','rna_human':'rna_fixture'},engine.Policy(),folder,quiet,checkpoint_session=checkpoint)
                else:result=engine.blast_candidates(candidates,refs,bins,[],{'human':'human_fixture','rna_human':'rna_fixture'},engine.Policy(),folder,quiet,checkpoint_session=checkpoint)
            snapshots=[]
            for c in candidates:
                engine.summarize_specificity(c,engine.Policy())
                snapshots.append(dict(hits=list(c['offtarget_hits']),count=len(c['offtarget_hits']),specificity=c['specificity'],reasons=c['reject_reasons']))
            return result,snapshots,calls
        finally:store.close()

    def test_resume_skips_successful_batch_preserves_every_hit_and_cap_subject(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);db=root/'saved.sqlite'
            with self.assertRaisesRegex(RuntimeError,'Forced failed'):
                self.run_batches(root/'stopped',db,fail_at=2)
            with sqlite3.connect(db) as con:
                self.assertEqual(con.execute('SELECT COUNT(*) FROM checkpoints').fetchone()[0],1)
            resumed,hits,calls=self.run_batches(root/'resumed',db,resume=True)
            fresh,fresh_hits,fresh_calls=self.run_batches(root/'fresh',root/'fresh.sqlite')
            self.assertEqual(len(calls),5);self.assertEqual(len(fresh_calls),6)
            self.assertEqual(hits,fresh_hits)
            for actual,expected in zip(resumed,fresh):
                for key in ('raw_alignments','retained_alignments','query_sha256','output_sha256','output_sha256_scope','capped_candidate_ids'):
                    self.assertEqual(actual[key],expected[key],key)
            self.assertEqual(resumed[0]['output_sha256_scope'],'ORDERED_BATCH_SHA256_AND_LENGTH_MANIFEST')
            self.assertTrue(resumed[0]['batches'][0]['reused_from_checkpoint'])
            final,final_hits,final_calls=self.run_batches(root/'already_saved',db,resume=True)
            self.assertEqual(final_calls,[]);self.assertEqual(final_hits,fresh_hits)

    def test_interrupted_parser_rolls_back_partial_hit_and_subject_insertions(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);db=root/'saved.sqlite'
            with self.assertRaises((TypeError,KeyError,ValueError)):
                self.run_batches(root/'stopped',db,invalid_at=2)
            with sqlite3.connect(db) as con:
                self.assertEqual(con.execute('SELECT DISTINCT candidate FROM hits').fetchall(),[('candidate_0000',)])
                self.assertEqual(con.execute('SELECT DISTINCT candidate FROM subjects').fetchall(),[('candidate_0000',)])
            _,actual,calls=self.run_batches(root/'resumed',db,resume=True)
            _,expected,_=self.run_batches(root/'fresh',root/'fresh.sqlite')
            self.assertEqual(actual,expected);self.assertEqual(len(calls),5)

    def test_commit_before_raw_cleanup_survives_interruption_without_researching(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);db=root/'saved.sqlite'
            with self.assertRaisesRegex(RuntimeError,'cleanup'):
                self.run_batches(root/'stopped',db,cleanup_failure=True)
            raw=next((root/'stopped').rglob('*.tsv'));digest=engine.file_sha(raw)
            _,_,calls=self.run_batches(root/'resumed',db,resume=True)
            self.assertEqual(len(calls),5);self.assertEqual(engine.file_sha(raw),digest)

    def test_same_row_count_index_flag_corruption_refuses_resume(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);db=root/'saved.sqlite'
            self.run_batches(root/'fresh',db)
            with sqlite3.connect(db) as con:con.execute('UPDATE hits SET high=1-high WHERE rowid=1')
            with self.assertRaisesRegex(ValueError,'checksum'):HitStore(db,resume=True)

    def test_subject_corruption_refuses_resume(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);db=root/'saved.sqlite'
            self.run_batches(root/'fresh',db)
            with sqlite3.connect(db) as con:con.execute("UPDATE subjects SET subject='altered' WHERE candidate='candidate_0000' AND source='human'")
            with self.assertRaisesRegex(ValueError,'subject evidence'):HitStore(db,resume=True)

    def test_missing_row_and_corrupt_metadata_refuse_resume(self):
        for mutation in ('DELETE FROM hits WHERE rowid=1',"UPDATE checkpoints SET payload='{}' WHERE key='blast:human:1'"):
            with self.subTest(mutation=mutation),tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp);db=root/'saved.sqlite';self.run_batches(root/'fresh',db)
                with sqlite3.connect(db) as con:con.execute(mutation)
                with self.assertRaises(ValueError):HitStore(db,resume=True)

    def test_deleted_batch_marker_cannot_leave_unauthenticated_evidence(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);db=root/'saved.sqlite';self.run_batches(root/'fresh',db)
            with sqlite3.connect(db) as con:con.execute("DELETE FROM checkpoints WHERE key='blast:human:1'")
            with self.assertRaisesRegex(ValueError,'subject evidence|gap|coverage'):HitStore(db,resume=True)

    def test_valid_database_from_different_identity_cannot_be_substituted(self):
        with tempfile.TemporaryDirectory() as tmp:
            db=Path(tmp)/'saved.sqlite';store=HitStore(db,identity_digest='original-identity')
            store.commit_checkpoint('local',dict(status='COMPLETE'));store.close()
            with self.assertRaisesRegex(ValueError,'different design inputs'):
                HitStore(db,resume=True,identity_digest='different-identity')
            verified=HitStore(db,resume=True,identity_digest='original-identity');verified.close()

    def test_sigkill_mid_transaction_keeps_only_previously_committed_batch(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);db=root/'saved.sqlite';store=HitStore(db)
            row=engine.enrich_hit(hit(),cand(),engine.Policy())
            store.hits('saved').append(row);store.record_subject('human','saved','chr7')
            store.commit_checkpoint('blast:human:1',dict(status='COMPLETE'));store.close()
            script="""import json,sys,os,signal
from pathlib import Path
from specificity_store import HitStore
s=HitStore(Path(sys.argv[1]),resume=True)
h=json.loads(sys.argv[2]);hits=s.hits('partial')
for i in range(12000):
 h['start']=i+200;h['end']=i+239;hits.append(h)
s.record_subject('human','partial','chr8')
os.kill(os.getpid(),signal.SIGKILL)
"""
            result=subprocess.run([sys.executable,'-c',script,str(db),json.dumps(row)],cwd=Path(engine.__file__).parent,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            self.assertEqual(result.returncode,-signal.SIGKILL,result.stderr.decode())
            reopened=HitStore(db,resume=True)
            try:
                self.assertEqual(len(reopened.hits('saved')),1);self.assertEqual(len(reopened.hits('partial')),0)
                self.assertEqual(reopened.connection.execute('SELECT COUNT(*) FROM subjects').fetchone()[0],1)
            finally:reopened.close()


class IndividualCheckpointTests(unittest.TestCase):
    def config(self):
        rng=random.Random(771)
        seq=''.join(rng.choice('ACGT') for _ in range(800))
        return dict(title='Checkpoint implementation fixture',references=[dict(id='TINY',sequence=seq,length=len(seq),features=[],mapping=None)],
            imagers=IMS,bin_width=400,tail_policy='strict',assignments={'TINY_B001':'DP1','TINY_B002':'DP2'},
            policy=dict(mode='local_preview',shortlist_per_bin=4,extension_shortlist_per_bin=2,screen_cross_primary=False),output_detail='standard')

    def test_folding_resume_exact_scientific_parity_no_repeated_local_or_saved_fold_calls(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);config=self.config();calls=[];original_fold=engine.Folding.fold
            def stop_after_two(fold,c):
                calls.append(c['id'])
                if len(calls)==3:raise RuntimeError('Forced folding interruption')
                return original_fold(fold,c)
            with patch.object(engine.Folding,'fold',stop_after_two):
                with self.assertRaisesRegex(RuntimeError,'folding interruption'):
                    engine.design(copy.deepcopy(config),root/'stopped',quiet,checkpoint_dir=root/'checkpoints')
            saved=set(calls[:2]);resumed_calls=[];original_pair=engine.Folding.pair
            def track_fold(fold,c):
                self.assertNotIn(c['id'],saved);resumed_calls.append(c['id']);return original_fold(fold,c)
            def track_pair(fold,a,b):
                self.assertEqual(a['id'],b['id'],'Disabled cross-primary check was called')
                return original_pair(fold,a,b)
            with patch.object(engine,'local_qc',side_effect=AssertionError('Saved local QC must not run again')),patch.object(engine.Folding,'fold',track_fold),patch.object(engine.Folding,'pair',track_pair):
                actual=engine.design(copy.deepcopy(config),root/'resumed',quiet,checkpoint_dir=root/'checkpoints')
            expected=engine.design(copy.deepcopy(config),root/'fresh',quiet,checkpoint_dir=root/'fresh_checkpoints')
            for key in ('selected_ids','counts','candidate_summary','cross_primary_screening','self_dimer_screening','candidates','selection_search'):
                self.assertEqual(actual.get(key),expected.get(key),key)
            self.assertGreater(len(resumed_calls),0)
            with patch.object(engine,'local_qc',side_effect=AssertionError('Local repeated')),patch.object(engine.Folding,'fold',side_effect=AssertionError('Fold repeated')),patch.object(engine,'exact_pair_segments',side_effect=AssertionError('Stem repeated')):
                complete=engine.design(copy.deepcopy(config),root/'all_saved',quiet,checkpoint_dir=root/'checkpoints')
            self.assertEqual(complete['selected_ids'],expected['selected_ids'])
            progress=read_checkpoint_summary(root/'checkpoints')
            self.assertEqual(progress['status'],'INCOMPLETE');self.assertEqual(progress['folded_candidates'],progress['total_candidates'])
            self.assertFalse((root/'stopped'/'design_result.json').exists())

    def test_local_preview_ignores_unused_empty_database_paths_and_resumes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);config=self.config();config['databases']={'human':'','rna_human':''}
            first=engine.design(copy.deepcopy(config),root/'fresh',quiet,checkpoint_dir=root/'checkpoints')
            with patch.object(engine,'blast_candidates',side_effect=AssertionError('Preview must not run BLAST')),patch.object(engine,'local_qc',side_effect=AssertionError('Saved local checks repeated')),patch.object(engine.Folding,'fold',side_effect=AssertionError('Saved folding checks repeated')):
                resumed=engine.design(copy.deepcopy(config),root/'resumed',quiet,checkpoint_dir=root/'checkpoints')
            self.assertEqual(first['selected_ids'],resumed['selected_ids'])
            identity=json.loads((root/'checkpoints'/'checkpoint_manifest.json').read_text())['identity']
            self.assertEqual(identity['databases'],{})
            self.assertEqual(read_checkpoint_summary(root/'checkpoints')['total_batches'],0)


if __name__=='__main__':unittest.main()
