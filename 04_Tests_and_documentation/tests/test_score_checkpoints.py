"""Concise Standard evidence must preserve the full scientific decisions."""
import copy
import hashlib
import json
import os
from pathlib import Path
import signal
import sqlite3
import subprocess
import sys
import tempfile
import threading
import types
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
sys.path.insert(0,str(Path(__file__).resolve().parent))
import design_engine as engine
from specificity_store import ScoreStore, HitStore, ScoreHits
from test_design import hit, cand
from test_blast_batches import fixture_candidates, fixture_output, expected_batch_environment
import test_design_checkpoints as checkpoint_fixtures


def quiet(*a,**k):pass

def enriched(**kwargs):
    value=engine.enrich_hit(hit(**kwargs),cand(),engine.Policy())
    return value


class ScoreParityTests(unittest.TestCase):
    def compare(self,rows,policy=None):
        policy=policy or engine.Policy()
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);score=ScoreStore(root/'scores.sqlite');raw=HitStore(root/'hits.sqlite')
            try:
                for row in rows:score.hits('c').append(row);raw.hits('c').append(row)
                score.stage_scores('local',['c'],policy);score.commit_checkpoint('local',{'status':'COMPLETE'})
                a=cand('c');b=cand('c');a['offtarget_hits']=score.hits('c');b['offtarget_hits']=raw.hits('c')
                engine.summarize_specificity(a,policy);engine.summarize_specificity(b,policy)
                self.assertEqual(a['specificity'],b['specificity'])
                self.assertEqual(a['reject_reasons'],b['reject_reasons'])
                self.assertEqual(a['offtarget_alignment_count'],b['offtarget_alignment_count'])
                return a['specificity']
            finally:score.close();raw.close()

    def test_unsorted_bridge_merges_both_queries_but_not_adjacent_or_different_sources(self):
        rows=[enriched(start=100,end=120),enriched(start=140,end=160,full_primary=True),
              enriched(start=115,end=145,full_primary=True),enriched(start=161,end=180),
              enriched(start=100,end=160,source='construct'),enriched(start=100,end=160,bound_strand='+'),
              enriched(start=100,end=160,reference='NT_ALT.1')]
        self.assertEqual(self.compare(rows)['review_site_clusters'],5)

    def test_aliases_merge_but_versioned_noncanonical_contigs_do_not(self):
        rows=[enriched(),enriched(reference='ref|NC_000007.14|'),enriched(reference='NC_000007.13')]
        self.assertEqual(self.compare(rows)['review_site_clusters'],2)

    def test_zero_risk_reviewed_clusters_and_allowed_duplicates_count_exactly(self):
        safe=enriched();safe.update(high_tm=False,moderate_tm=False,unknown_long=False,strong_rna=False,rna_exact_risk=False)
        result=self.compare([safe,enriched(allowed=True),enriched(allowed=True)])
        self.assertEqual(result['allowed_hsps'],2);self.assertEqual(result['high_tm_sites'],0)
        self.assertEqual(result['review_site_clusters'],1);self.assertGreater(result['risk_penalty'],0)

    def test_report_only_keeps_ignored_hsp_count_and_other_sources(self):
        result=self.compare([enriched(),enriched(full_primary=True),enriched(source='construct')],engine.Policy(human_policy='report_only'))
        self.assertEqual(result['ignored_human_review_hsps'],2);self.assertEqual(result['review_site_clusters'],1)

    def test_rna_exact_and_strong_flags_are_policy_conditional(self):
        exact=enriched(source='rna_human',rna=True,bound_strand='RNA');exact.update(strong_rna=False,rna_exact_risk=True)
        strong=dict(exact,start=500,end=539,strong_rna=True,rna_exact_risk=False)
        self.assertEqual(self.compare([exact,strong],engine.Policy(avoid_gene_rna=False))['rna_risk_sites'],1)
        self.assertEqual(self.compare([exact,strong],engine.Policy(avoid_gene_rna=True))['rna_risk_sites'],2)

    def test_expected_target_requires_perfect_allowed_human_arm(self):
        variants=[{},dict(full_primary=True),dict(identity=99.9),dict(query_coverage=99.9),dict(allowed=False),dict(source='construct')]
        for index,variant in enumerate(variants):
            with self.subTest(variant=variant),tempfile.TemporaryDirectory() as tmp:
                store=ScoreStore(Path(tmp)/'score.sqlite')
                try:
                    store.hits('c').append(enriched(**dict(dict(allowed=True),**variant)))
                    store.stage_scores('local',['c'],engine.Policy())
                    self.assertEqual(store.hits('c').has_expected_genomic_target(),index==0)
                finally:store.close()

    def test_local_database_source_collision_unions_intervals_without_double_count(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);store=ScoreStore(root/'score.sqlite',deferred_local_sources=['construct']);raw=HitStore(root/'raw.sqlite')
            try:
                local=[enriched(source='construct',start=100,end=120),enriched(source='construct',start=140,end=160)]
                blast=[enriched(source='construct',start=115,end=145,full_primary=True)]
                for h in local:store.hits('c').append(h);raw.hits('c').append(h)
                store.stage_scores('local',['c'],engine.Policy());store.commit_checkpoint('local',dict(status='COMPLETE'))
                store.restore_local_source('construct','c')
                for h in blast:store.hits('c').append(h);raw.hits('c').append(h)
                store.stage_scores('blast:construct:1',['c'],engine.Policy(),'construct');store.commit_checkpoint('blast:construct:1',dict(status='COMPLETE'))
                self.assertEqual(store.hits('c').specificity_counts(engine.Policy()),raw.hits('c').specificity_counts(engine.Policy()))
                self.assertEqual(len(store.hits('c')),len(raw.hits('c')))
                self.assertEqual(store.hits('c').specificity_counts(engine.Policy())['review_site_clusters'],1)
            finally:store.close();raw.close()

    def test_distinct_subject_count_is_union_across_both_queries_even_without_retained_hits(self):
        with tempfile.TemporaryDirectory() as tmp:
            store=ScoreStore(Path(tmp)/'score.sqlite')
            try:
                for subject in ['NM_1','NM_1','NM_2','NM_3']:store.record_subject('rna_human','c',subject)
                store.stage_scores('blast:rna_human:1',['c'],engine.Policy(),'rna_human')
                store.commit_checkpoint('blast:rna_human:1',dict(status='COMPLETE'))
                self.assertEqual(store.capped_candidates('rna_human',3),['c'])
                self.assertEqual(store.capped_candidates('rna_human',4),[])
                self.assertEqual(len(store.hits('c')),0)
            finally:store.close()


class ScoreDurabilityTests(unittest.TestCase):
    def prepared(self,path):
        store=ScoreStore(path,identity_digest='fixture')
        store.hits('c').append(enriched())
        store.stage_scores('blast:human:1',['c'],engine.Policy(),'human')
        store.commit_checkpoint('blast:human:1',dict(status='COMPLETE'))
        store.close()

    def test_scores_have_no_sequence_or_alignment_payload_and_work_file_reclaims_each_batch(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'score.sqlite';store=ScoreStore(path)
            try:
                work_sizes=[]
                for batch in range(20):
                    for n in range(500):store.hits(str(batch)).append(enriched(start=100+n,end=140+n))
                    self.assertGreater(store._work.execute('SELECT COUNT(*) FROM hits').fetchone()[0],0)
                    store.stage_scores(f'blast:human:{batch}',[str(batch)],engine.Policy(),'human')
                    store.commit_checkpoint(f'blast:human:{batch}',dict(status='COMPLETE'))
                    self.assertEqual(store._work.execute('SELECT COUNT(*) FROM hits').fetchone()[0],0)
                    work_sizes.append(store.work_path.stat().st_size)
                self.assertEqual(min(work_sizes),max(work_sizes))
                self.assertLess(path.stat().st_size,100000)
                blob=path.read_bytes()
                for forbidden in (b'qseq',b'sseq',b'GCGCATGCGCGCATGC',b'deferred_intervals'):self.assertNotIn(forbidden,blob)
                tables={r[0] for r in store.connection.execute("SELECT name FROM sqlite_master WHERE type='table'")}
                self.assertEqual(tables,{'scores','checkpoints','checkpoint_store','checkpoint_identity'})
                with self.assertRaisesRegex(RuntimeError,'Full audit'):list(store.hits('0'))
            finally:store.close()
            self.assertFalse(path.with_name('score.working.sqlite').exists())

    def test_corrupt_or_uncovered_scores_refuse_resume(self):
        for mutation in ("UPDATE scores SET payload='{}'",'DELETE FROM scores','DELETE FROM checkpoints'):
            with self.subTest(mutation=mutation),tempfile.TemporaryDirectory() as tmp:
                path=Path(tmp)/'score.sqlite';self.prepared(path)
                with sqlite3.connect(path) as db:db.execute(mutation)
                with self.assertRaises(ValueError):ScoreStore(path,resume=True,identity_digest='fixture')

    def test_wrong_identity_or_legacy_raw_store_refuse_resume(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'score.sqlite';self.prepared(path)
            with self.assertRaisesRegex(ValueError,'different design inputs'):ScoreStore(path,resume=True,identity_digest='changed')
            legacy=Path(tmp)/'legacy.sqlite';HitStore(legacy).close()
            with self.assertRaisesRegex(ValueError,'Unsupported score checkpoint schema'):ScoreStore(legacy,resume=True)

    def test_sigkill_mid_uncommitted_reduction_keeps_only_prior_scores(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'score.sqlite';self.prepared(path)
            program='''import os,signal,sys,json
from specificity_store import ScoreStore
from design_engine import Policy
s=ScoreStore(sys.argv[1],resume=True,identity_digest='fixture')
s.hits('partial').append(json.loads(sys.argv[2]))
s.stage_scores('blast:human:2',['partial'],Policy(),'human')
os.kill(os.getpid(),signal.SIGKILL)
'''
            proc=subprocess.run([sys.executable,'-c',program,str(path),json.dumps(enriched())],cwd=Path(engine.__file__).parent,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            self.assertEqual(proc.returncode,-signal.SIGKILL,proc.stderr.decode())
            store=ScoreStore(path,resume=True,identity_digest='fixture')
            try:self.assertEqual(len(store.hits('c')),1);self.assertEqual(len(store.hits('partial')),0);self.assertEqual(store.checkpoint_count('blast:'),1)
            finally:store.close()


class ScoreBatchTests(unittest.TestCase):
    def test_streamed_partial_process_failure_guard_and_cancel_rollback_then_recompute(self):
        from runtime_safety import MemoryLimitExceeded
        native_process=engine.run_process;native_append=ScoreHits.append;native_which=engine.shutil.which
        for failure in ('nonzero','guard','cancel'):
            with self.subTest(failure=failure),tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp);path=root/'scores.sqlite';partial_seen=threading.Event()
                calls=[]
                def attempt(folder,*,resume=False,fail=False):
                    folder.mkdir();candidates,refs,bins=fixture_candidates(1)
                    store=ScoreStore(path if folder.name!='fresh' else root/'fresh.sqlite',resume=resume,identity_digest='fixture')
                    cancel=threading.Event()
                    for candidate in candidates:candidate['offtarget_hits']=store.hits(candidate['id'])
                    class Guard:
                        def checkpoint(self,stage,force=False):
                            if fail and failure=='guard' and stage=='BLAST subprocess' and partial_seen.is_set():
                                raise MemoryLimitExceeded('Synthetic streaming resource stop',{'status':'STOP_TEST'})
                    def append(hits,hit):
                        native_append(hits,hit)
                        if fail:
                            partial_seen.set()
                            if failure=='cancel':cancel.set()
                    def execute(command,outlog,current_cancel=None,resource_guard=None,*,env_overrides=None):
                        self.assertEqual(env_overrides,expected_batch_environment(command))
                        calls.append(command)
                        data=fixture_output(command)
                        if fail:data=data.splitlines(keepends=True)[0]
                        delay=.1 if fail and failure=='nonzero' else 30 if fail else 0
                        code=7 if fail and failure=='nonzero' else 0
                        source='import sys,time; f=open(sys.argv[1],"wb",buffering=0); f.write(bytes.fromhex(sys.argv[2])); time.sleep(float(sys.argv[3])); f.close(); sys.exit(int(sys.argv[4]))'
                        producer=[sys.executable,'-c',source,command[command.index('-out')+1],data.hex(),str(delay),str(code)]
                        return native_process(producer,outlog,current_cancel,resource_guard,env_overrides=env_overrides)
                    try:
                        with patch.object(engine,'BLAST_BATCH_CANDIDATES',1),patch.object(engine.shutil,'which',side_effect=lambda name,*args,**kwargs:'/fixture/blastn' if name=='blastn' else native_which(name,*args,**kwargs)),patch.object(engine.subprocess,'check_output',return_value='fixture'),patch.object(engine,'run_process',side_effect=execute),patch.object(ScoreHits,'append',append):
                            result=engine.blast_candidates(candidates,refs,bins,[],{'human':'human_fixture'},engine.Policy(),folder,quiet,cancel=cancel,resource_guard=Guard(),checkpoint_session=types.SimpleNamespace(update_progress=quiet))
                        engine.summarize_specificity(candidates[0],engine.Policy())
                        return result,candidates[0]['specificity']
                    finally:store.close()
                expected_error=MemoryLimitExceeded if failure=='guard' else engine.Cancelled if failure=='cancel' else RuntimeError
                with self.assertRaises(expected_error):attempt(root/'failed',fail=True)
                self.assertTrue(partial_seen.is_set(),'Fixture failed before any alignment reached the score store')
                store=ScoreStore(path,resume=True,identity_digest='fixture')
                try:
                    self.assertEqual(store.checkpoint_count('blast:'),0)
                    self.assertEqual(store.connection.execute('SELECT COUNT(*) FROM scores').fetchone()[0],0)
                    self.assertEqual(len(store.hits('candidate_0000')),0)
                finally:store.close()
                manifest=json.loads(next((root/'failed').rglob('search_manifest.json')).read_text())
                self.assertEqual(manifest['status'],'INCOMPLETE');self.assertEqual(manifest['completed_batches'],0)
                self.assertFalse(list((root/'failed').rglob('*_complete.json')))
                self.assertFalse(list((root/'failed').rglob('*_candidates_blast.tsv')))
                actual,score=attempt(root/'resumed',resume=True)
                self.assertEqual(len(calls),2,'Interrupted batch was not recomputed')
                expected,baseline=attempt(root/'fresh')
                self.assertEqual(score,baseline)
                self.assertEqual(actual[0]['batches'][0]['output_sha256'],expected[0]['batches'][0]['output_sha256'])

    def run_fixture(self,root,store_path,*,resume=False,fail_at=None,invalid_at=None,kind=ScoreStore):
        root.mkdir();candidates,refs,bins=fixture_candidates(3)
        store=kind(store_path,resume=resume,identity_digest='fixture')
        calls=[]
        for c in candidates:c['offtarget_hits']=store.hits(c['id'])
        def execute(command,outlog,cancel=None,resource_guard=None,*,env_overrides=None):
            self.assertEqual(env_overrides,expected_batch_environment(command))
            calls.append(command);data=fixture_output(command)
            if len(calls)==invalid_at:data=data.splitlines(keepends=True)[0]+b'invalid\trow\n'
            Path(command[command.index('-out')+1]).write_bytes(data);outlog.write_text('fixture')
            if len(calls)==fail_at:raise RuntimeError('Deliberate failed BLAST')
            return {'status':'COMPLETE','returncode':0}
        try:
            with patch.object(engine,'BLAST_BATCH_CANDIDATES',1),patch.object(engine.shutil,'which',return_value='/fixture/blastn'),patch.object(engine.subprocess,'check_output',return_value='fixture'),patch.object(engine,'run_process',side_effect=execute):
                result=engine.blast_candidates(candidates,refs,bins,[],{'human':'human_fixture','rna_human':'rna_fixture'},engine.Policy(),root,quiet,checkpoint_session=types.SimpleNamespace(update_progress=quiet))
            snapshots=[]
            for c in candidates:
                engine.summarize_specificity(c,engine.Policy())
                expected=(c['offtarget_hits'].has_expected_genomic_target() if isinstance(c['offtarget_hits'],ScoreHits) else any(h['source']=='human' and h['allowed'] and not h['full_primary'] and h['identity']==100 and h['query_coverage']==100 for h in c['offtarget_hits']))
                snapshots.append(dict(specificity=c['specificity'],reasons=c['reject_reasons'],count=c['offtarget_alignment_count'],expected=expected))
            return result,snapshots,calls
        finally:store.close()

    def test_resumed_scores_exactly_match_full_retained_alignment_evidence(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);path=root/'score.sqlite'
            with self.assertRaisesRegex(RuntimeError,'failed BLAST'):self.run_fixture(root/'failed',path,fail_at=2)
            self.assertEqual(list((root/'failed').rglob('*.tsv')),[])
            resumed,scores,calls=self.run_fixture(root/'resumed',path,resume=True)
            raw,baseline,_=self.run_fixture(root/'raw',root/'raw.sqlite',kind=HitStore)
            self.assertEqual(scores,baseline);self.assertEqual(len(calls),5)
            for actual,expected in zip(resumed,raw):
                for key in ('raw_alignments','retained_alignments','output_sha256','capped_candidate_ids'):self.assertEqual(actual[key],expected[key])
            _,saved,calls=self.run_fixture(root/'all_saved',path,resume=True)
            self.assertEqual(calls,[]);self.assertEqual(saved,baseline)

    def test_partial_parser_does_not_commit_or_leave_raw_alignment_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);path=root/'score.sqlite'
            with self.assertRaises((ValueError,TypeError,KeyError)):self.run_fixture(root/'failed',path,invalid_at=2)
            self.assertEqual(list((root/'failed').rglob('*.tsv')),[])
            store=ScoreStore(path,resume=True,identity_digest='fixture')
            try:self.assertEqual(store.checkpoint_count('blast:'),1);self.assertEqual(len(store.hits('candidate_0001')),0)
            finally:store.close()
            _,actual,_=self.run_fixture(root/'resumed',path,resume=True)
            _,expected,_=self.run_fixture(root/'fresh',root/'fresh.sqlite')
            self.assertEqual(actual,expected)

    def test_post_commit_cleanup_failure_preserves_saved_score_and_removes_raw_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);path=root/'score.sqlite'
            with patch.object(engine,'archive_transient',side_effect=RuntimeError('Deliberate cleanup failure')):
                with self.assertRaisesRegex(RuntimeError,'cleanup failure'):
                    self.run_fixture(root/'failed',path)
            self.assertEqual(list((root/'failed').rglob('*.tsv')),[])
            _,actual,calls=self.run_fixture(root/'resumed',path,resume=True)
            _,expected,_=self.run_fixture(root/'fresh',root/'fresh.sqlite')
            self.assertEqual(len(calls),5);self.assertEqual(actual,expected)

    def test_secondary_rollback_failure_never_masks_primary_blast_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);path=root/'score.sqlite'
            with patch.object(ScoreStore,'rollback',side_effect=RuntimeError('Secondary rollback error')):
                with self.assertRaisesRegex(RuntimeError,'Deliberate failed BLAST'):
                    self.run_fixture(root/'failed',path,fail_at=2)
            self.assertEqual(list((root/'failed').rglob('*.tsv')),[])
            manifest=json.loads(next((root/'failed').rglob('search_manifest.json')).read_text())
            self.assertEqual(manifest['error'],'Deliberate failed BLAST')
            self.assertIn('Secondary rollback error',manifest['batches'][-1]['cleanup_errors'][0])
            _,actual,calls=self.run_fixture(root/'resumed',path,resume=True)
            _,expected,_=self.run_fixture(root/'fresh',root/'fresh.sqlite')
            self.assertEqual(len(calls),5);self.assertEqual(actual,expected)

    def test_integrated_standard_and_full_audit_have_equal_selected_scientific_results(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);config=checkpoint_fixtures.IndividualCheckpointTests().config()
            compact=engine.design(copy.deepcopy(config),root/'standard',quiet,checkpoint_dir=root/'checkpoints')
            config['output_detail']='full_audit'
            audit=engine.design(copy.deepcopy(config),root/'audit',quiet)
            # Full audit additionally exports structures/per-base arrays. Its
            # scientific data must match after applying only the documented
            # concise serialization; independent count/cluster parity tests
            # above compare ScoreStore directly with retained raw evidence.
            from compact_outputs import concise_saved_result
            audit_scores=concise_saved_result(audit)
            for key in ('selected_ids','counts','candidates','selection_search','cross_primary_screening','self_dimer_screening'):
                self.assertEqual(compact[key],audit_scores[key],key)
            self.assertTrue((root/'checkpoints'/'specificity_scores.sqlite').is_file())
            self.assertFalse((root/'checkpoints'/'retained_hits.sqlite').exists())
            self.assertEqual(compact['cross_primary_screening']['evaluated'],0)
            self.assertGreater(compact['counts']['monomer_evaluations'],0)
            self.assertGreater(compact['counts']['self_dimer_evaluations'],0)

if __name__=='__main__':unittest.main()
