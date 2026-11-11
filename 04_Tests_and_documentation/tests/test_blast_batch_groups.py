"""Small candidate groups preserve every query and transaction boundary."""
import hashlib
import json
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
sys.path.insert(0,str(Path(__file__).resolve().parent))
import design_engine as engine
from specificity_store import ScoreStore,HitStore,ScoreHits,DiskHits
from test_blast_batches import fixture_candidates,fixture_output,query_records,expected_batch_environment


def quiet(*args,**kwargs):pass


class GroupCheckpointTests(unittest.TestCase):
    def attempt(self,root,path,*,kind=ScoreStore,count=10,batch_size=4,resume=False,
                fail_at=None,invalid_at=None,order=None,checkpoint=True):
        root.mkdir()
        candidates,refs,bins=fixture_candidates(count)
        if order is not None:candidates=[candidates[i] for i in order]
        store=kind(path,resume=resume,identity_digest='group-fixture')
        for c in candidates:c['offtarget_hits']=store.hits(c['id'])
        calls=[];raw_queries={};seen=[]
        hits_class=ScoreHits if kind is ScoreStore else DiskHits
        native_append=hits_class.append
        def append(hits,hit):
            seen.append(hits.candidate)
            return native_append(hits,hit)
        def execute(command,outlog,cancel=None,resource_guard=None,*,env_overrides=None):
            self.assertEqual(env_overrides,expected_batch_environment(command))
            records=query_records(command)
            calls.append([ident[:-2] for ident,_ in records if ident.endswith('_A')])
            data=fixture_output(command)
            if len(calls)==fail_at or len(calls)==invalid_at:
                # Enough valid output from the first candidate to reach the
                # store before the malformed row or unsuccessful child exit.
                data=b''.join(data.splitlines(keepends=True)[:4])
                if len(calls)==invalid_at:data+=b'invalid\trow\n'
            else:
                source='rna_human' if 'rna' in Path(command[command.index('-db')+1]).name else 'human'
                for line in data.splitlines(keepends=True):
                    ident=line.split(b'\t',1)[0].decode()
                    key=(source,ident)
                    entry=raw_queries.setdefault(key,[hashlib.sha256(),0])
                    entry[0].update(line);entry[1]+=1
            Path(command[command.index('-out')+1]).write_bytes(data)
            outlog.write_text('Synthetic grouped-checkpoint fixture; no real genome search.\n')
            if len(calls)==fail_at:raise RuntimeError('Deliberate partial grouped process failure')
            return {'status':'COMPLETE','returncode':0}
        try:
            with patch.object(engine,'BLAST_BATCH_CANDIDATES',batch_size),patch.object(engine.shutil,'which',return_value='/fixture/blastn'),patch.object(engine.subprocess,'check_output',return_value='fixture'),patch.object(engine,'run_process',side_effect=execute),patch.object(hits_class,'append',append):
                result=engine.blast_candidates(candidates,refs,bins,[],{'human':'human_fixture','rna_human':'rna_fixture'},engine.Policy(),root,quiet,
                    detail='standard' if kind is ScoreStore else 'full_audit',
                    checkpoint_session=types.SimpleNamespace(update_progress=quiet) if checkpoint else None)
            scores=[]
            for c in candidates:
                engine.summarize_specificity(c,engine.Policy())
                scores.append(dict(id=c['id'],specificity=c['specificity'],reasons=c['reject_reasons'],alignments=c['offtarget_alignment_count']))
            return result,scores,calls,{k:(v[0].hexdigest(),v[1]) for k,v in raw_queries.items()}
        except BaseException as error:
            error.fixture_calls=calls
            error.fixture_seen=seen
            raise
        finally:store.close()

    def test_partial_second_group_rolls_back_every_candidate_and_retries_same_group(self):
        for kind,batch_size in ((ScoreStore,2),(HitStore,2),(ScoreStore,4),(HitStore,4)):
            for failure in ('process','parser'):
                with self.subTest(kind=kind.__name__,batch_size=batch_size,failure=failure),tempfile.TemporaryDirectory() as tmp:
                    root=Path(tmp);path=root/'saved.sqlite'
                    kwargs={'fail_at':2} if failure=='process' else {'invalid_at':2}
                    with self.assertRaises((RuntimeError,ValueError,TypeError,KeyError)) as stopped:
                        self.attempt(root/'failed',path,kind=kind,batch_size=batch_size,**kwargs)
                    error=stopped.exception
                    self.assertEqual(error.fixture_calls,[[f'candidate_{i:04d}' for i in range(start,start+batch_size)] for start in (0,batch_size)])
                    # Full audit deliberately verifies a subprocess before
                    # parsing its file, while Standard streams partial rows.
                    if kind is ScoreStore or failure=='parser':self.assertIn(f'candidate_{batch_size:04d}',error.fixture_seen)
                    store=kind(path,resume=True,identity_digest='group-fixture')
                    try:
                        self.assertEqual(store.checkpoint_count('blast:'),1)
                        saved=store.checkpoint('blast:human:1')
                        self.assertEqual(saved['candidate_ids'],error.fixture_calls[0])
                        self.assertIsNone(store.checkpoint('blast:human:2'))
                        self.assertGreater(len(store.hits('candidate_0000')),0)
                        for index in range(batch_size,10):self.assertEqual(len(store.hits(f'candidate_{index:04d}')),0)
                        if kind is ScoreStore:self.assertEqual(store.connection.execute('SELECT DISTINCT candidate FROM scores ORDER BY candidate').fetchall(),[(f'candidate_{i:04d}',) for i in range(batch_size)])
                    finally:store.close()
                    manifest=json.loads(next((root/'failed').rglob('search_manifest.json')).read_text())
                    self.assertEqual(manifest['status'],'INCOMPLETE')
                    self.assertEqual(manifest['completed_batches'],1)
                    self.assertEqual(manifest['batches'][1]['status'],'INCOMPLETE')
                    self.assertEqual(len(list((root/'failed').rglob('*_complete.json'))),1)
                    actual,scores,calls,_=self.attempt(root/'resumed',path,kind=kind,batch_size=batch_size,resume=True)
                    expected,baseline,fresh_calls,_=self.attempt(root/'fresh',root/'fresh.sqlite',kind=kind,batch_size=batch_size)
                    self.assertEqual(calls,fresh_calls[1:])
                    self.assertEqual(calls[0],error.fixture_calls[1])
                    self.assertEqual(scores,baseline)
                    for a,b in zip(actual,expected):
                        for key in ('raw_alignments','retained_alignments','output_sha256','query_sha256','capped_candidate_ids'):self.assertEqual(a[key],b[key],key)
                    self.assertTrue(actual[0]['batches'][0]['reused_from_checkpoint'])
                    final,final_scores,final_calls,_=self.attempt(root/'all_saved',path,kind=kind,batch_size=batch_size,resume=True)
                    self.assertEqual(final_calls,[]);self.assertEqual(final_scores,baseline)

    def test_last_partial_group_covers_both_queries_and_every_candidate_once(self):
        for kind,batch_size in ((ScoreStore,2),(HitStore,2),(ScoreStore,4),(HitStore,4)):
            with self.subTest(kind=kind.__name__,batch_size=batch_size),tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp);runs,_,calls,raw=self.attempt(root/'completed',root/'scores.sqlite',kind=kind,count=9,batch_size=batch_size)
                groups=[batch_size]*(9//batch_size)+[9%batch_size]
                self.assertEqual([len(group) for group in calls],groups*2)
                self.assertEqual(len(raw),36)
                for run in runs:
                    self.assertEqual(run['batch_count'],len(groups))
                    self.assertEqual([r['query_count'] for r in run['batches']],[2*n for n in groups])
                    self.assertEqual([i for r in run['batches'] for i in r['candidate_ids']],[f'candidate_{i:04d}' for i in range(9)])

    def test_reordered_candidates_cannot_reuse_group_checkpoint(self):
        for kind in (ScoreStore,HitStore):
            with self.subTest(kind=kind.__name__),tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp);path=root/'saved.sqlite'
                self.attempt(root/'saved',path,kind=kind,count=5)
                before=hashlib.sha256(path.read_bytes()).hexdigest()
                with self.assertRaisesRegex(ValueError,'does not match the current candidate batch') as error:
                    self.attempt(root/'reordered',path,kind=kind,count=5,resume=True,order=[1,0,2,3,4])
                self.assertEqual(error.exception.fixture_calls,[])
                self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(),before)
                self.assertFalse(list((root/'reordered').rglob('*_complete.json')))
                manifest=json.loads(next((root/'reordered').rglob('search_manifest.json')).read_text())
                self.assertEqual(manifest['status'],'INCOMPLETE')

    def test_one_two_and_four_candidate_batches_match_raw_queries_and_exact_scores_in_both_modes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);results=[]
            for kind in (ScoreStore,HitStore):
                for size in (1,2,4):
                    label=f'{kind.__name__}_{size}'
                    results.append(self.attempt(root/label,root/(label+'.sqlite'),kind=kind,count=10,batch_size=size,checkpoint=False))
            baseline=results[0]
            for result in results[1:]:
                self.assertEqual(result[1],baseline[1])
                self.assertEqual(result[3],baseline[3])
                for actual,expected in zip(result[0],baseline[0]):
                    for key in ('raw_alignments','retained_alignments','query_sha256','output_sha256','capped_candidate_ids'):self.assertEqual(actual[key],expected[key],key)
            # Scientific checks include allowed intended hits, several additional
            # loci on that chromosome, full-primary junctions and RNA orientation.
            scores=baseline[1][0]['specificity']
            self.assertGreater(scores['allowed_hsps'],0)
            self.assertGreater(scores['review_site_clusters'],1)
            self.assertGreater(scores['strong_rna_sites'],0)


if __name__=='__main__':unittest.main()
