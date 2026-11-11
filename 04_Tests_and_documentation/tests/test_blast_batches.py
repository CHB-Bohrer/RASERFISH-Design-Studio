"""Bounded BLAST dispatch must preserve every scientific search and hit rule."""
import copy
import gzip
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
import threading
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
sys.path.insert(0,str(Path(__file__).resolve().parent))
import design_engine as engine
from specificity_store import HitStore
from test_blast_adapter import inputs


def fixture_candidates(count):
    candidate,reference,bins=inputs()
    candidates=[]
    for index in range(count):
        item=copy.deepcopy(candidate);item['id']=f'candidate_{index:04d}'
        candidates.append(item)
    return candidates,[reference],bins


def query_records(command):
    lines=Path(command[command.index('-query')+1]).read_text().splitlines()
    return [(lines[i][1:],lines[i+1]) for i in range(0,len(lines),2)]


def expected_batch_environment(command):
    records=query_records(command)
    lengths={ident:len(sequence) for ident,sequence in records}
    arms=[ident for ident,_ in records if ident.endswith('_A')]
    limit=max(lengths[ident] for ident in arms)+1
    if limit>min(lengths[ident]+lengths[ident[:-2]+'_F'] for ident in arms):limit=1
    return {'BATCH_SIZE':str(limit)}


def fixture_output(command):
    """Intended locus, two extra sites on the SAME chromosome, full junctions,
    and both RNA orientations. Multiple sites must never be collapsed to one.
    """
    candidate,_,_=inputs();rows=[]
    rna='rna' in Path(command[command.index('-db')+1]).name
    for ident,sequence in query_records(command):
        full=ident.endswith('_F')
        qstart=len(candidate['prefix'])+1 if full else 1
        for start in (100,1000,2000):
            for reverse in ([False,True] if rna else [False]):
                end=start+len(candidate['arm'])-1
                subject_start,subject_end=(end,start) if reverse else (start,end)
                rows.append([ident,'NM_FIXTURE.1' if rna else 'NC_000007.14',100,len(candidate['arm']),0,0,
                    qstart,qstart+len(candidate['arm'])-1,subject_start,subject_end,1e-10,80,len(sequence),
                    candidate['arm'],candidate['arm']])
        if full:
            q=sequence[:24]
            rows.append([ident,'NM_FIXTURE.1' if rna else 'NC_000007.14',100,len(q),0,0,1,len(q),
                3023 if rna else 3000,3000 if rna else 3023,1e-6,40,len(sequence),q,q])
    return ''.join('\t'.join(map(str,row))+'\n' for row in rows).encode()


class BatchTests(unittest.TestCase):
    def run_fixture(self,root,count=21,batch_size=10,detail='standard',fail_at=None,databases=None):
        candidates,refs,bins=fixture_candidates(count)
        store=HitStore(root/'hits.sqlite')
        for candidate in candidates:candidate['offtarget_hits']=store.hits(candidate['id'])
        calls=[];outputs=[]
        def fake(command,outlog,cancel=None,resource_guard=None,*,env_overrides=None):
            self.assertEqual(env_overrides,expected_batch_environment(command))
            calls.append(command)
            if len(calls)>1 and detail=='standard':
                previous=Path(calls[-2][calls[-2].index('-out')+1])
                self.assertFalse(previous.exists(),'Previous successful batch must be processed and removed before launching another')
                self.assertGreater(store.connection.execute('SELECT COUNT(*) FROM hits').fetchone()[0],0)
            output=fixture_output(command);outputs.append(output)
            Path(command[command.index('-out')+1]).write_bytes(output)
            outlog.write_text('Explicit fixture, not a genomic result.\n')
            if len(calls)==fail_at:raise RuntimeError('Fixture BLAST terminated by signal 9 (SIGKILL)')
        try:
            with patch.object(engine,'BLAST_BATCH_CANDIDATES',batch_size),patch.object(engine.shutil,'which',return_value='/fixture/blastn'),patch.object(engine,'run_process',side_effect=fake),patch.object(engine.subprocess,'check_output',return_value='fixture BLAST version'):
                runs=engine.blast_candidates(candidates,refs,bins,[],databases or {'human':'human_fixture','rna_human':'rna_fixture'},engine.Policy(),root,lambda *a,**k:None,detail=detail)
            snapshots=[]
            for candidate in candidates:
                engine.summarize_specificity(candidate,engine.Policy())
                snapshots.append(dict(id=candidate['id'],hits=list(candidate['offtarget_hits']),specificity=candidate['specificity'],reasons=candidate['reject_reasons']))
            return runs,snapshots,calls,outputs
        finally:store.close()

    def test_batched_and_single_search_have_exact_hit_and_site_parity(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);(root/'many').mkdir();(root/'one').mkdir()
            batched,batched_hits,calls,_=self.run_fixture(root/'many')
            single,single_hits,_,_=self.run_fixture(root/'one',batch_size=1000)
            self.assertEqual(batched_hits,single_hits)
            self.assertEqual(len(calls),6)
            for result,baseline in zip(batched,single):
                for key in ('raw_alignments','retained_alignments','query_sha256','output_sha256','capped_candidate_ids'):
                    self.assertEqual(result[key],baseline[key],key)
            first=batched_hits[0]
            self.assertTrue(any(hit['allowed'] and hit['start']==100 for hit in first['hits']))
            for start in (1000,2000):
                self.assertTrue(any(hit['source']=='human' and not hit['allowed'] and hit['start']==start for hit in first['hits']))
            self.assertTrue(any(hit['full_primary'] and hit['query_start']==1 and hit['start']==3000 for hit in first['hits']))
            self.assertTrue(any(hit['source']=='rna_human' for hit in first['hits']))
            self.assertTrue(all(hit['bound_strand']=='RNA' for hit in first['hits'] if hit['source']=='rna_human'))

    def test_all_600_candidates_and_1200_queries_are_searched_with_unchanged_flags(self):
        with tempfile.TemporaryDirectory() as temp:
            self.assertEqual(engine.BLAST_BATCH_CANDIDATES,4)
            runs,_,calls,_=self.run_fixture(Path(temp),count=600,batch_size=engine.BLAST_BATCH_CANDIDATES,databases={'human':'human_fixture'})
            self.assertEqual(len(calls),150)
            self.assertEqual(sum(batch['query_count'] for batch in runs[0]['batches']),1200)
            ids=[ident for batch in runs[0]['batches'] for ident in batch['candidate_ids']]
            self.assertEqual(len(set(ids)),600)
            expected={'-task':'blastn-short','-strand':'both','-dust':'no','-soft_masking':'false','-word_size':'7',
                '-evalue':'1000','-reward':'1','-penalty':'-3','-gapopen':'5','-gapextend':'2','-max_target_seqs':'100000','-num_threads':'2',
                '-outfmt':'6 '+' '.join(engine.BLAST_FIELDS)}
            for command in calls:
                self.assertEqual(command[0],'/fixture/blastn')
                self.assertNotIn('-max_hsps',command)
                for flag,value in expected.items():self.assertEqual(command[command.index(flag)+1],value)

    def test_success_receipts_contain_hashes_counts_and_complete_status(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);runs,_,calls,outputs=self.run_fixture(root,databases={'human':'human_fixture'})
            manifest=json.loads((root/runs[0]['search_manifest']).read_text())
            self.assertEqual(manifest['status'],'COMPLETE')
            self.assertEqual(manifest['completed_batches'],3)
            self.assertEqual(manifest['query_count'],42)
            self.assertFalse(manifest['reuse_allowed'])
            self.assertEqual(manifest['output_sha256'],hashlib.sha256(b''.join(outputs)).hexdigest())
            for receipt,output,command in zip(manifest['batches'],outputs,calls):
                self.assertEqual(receipt['environment_overrides'],{'BATCH_SIZE':str(len(inputs()[0]['arm'])+1)})
                self.assertEqual(receipt['native_batch_strategy'],'ONE_CANDIDATE_PAIR_OR_SMALLER')
                self.assertEqual(receipt['status'],'COMPLETE')
                self.assertEqual(receipt['output_sha256'],hashlib.sha256(output).hexdigest())
                self.assertIsNone(receipt['raw_file']);self.assertIsNone(receipt['query_file'])
            self.assertEqual(len(list(root.rglob('*_complete.json'))),3)

    def test_full_audit_archives_are_unique_and_reproduce_each_batch(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);runs,_,_,outputs=self.run_fixture(root,detail='full_audit')
            receipts=[batch for run in runs for batch in run['batches']]
            self.assertEqual(len(set(batch['raw_file'] for batch in receipts)),len(receipts))
            for batch,expected in zip(receipts,outputs):
                with gzip.open(root/batch['raw_file'],'rb') as stream:self.assertEqual(stream.read(),expected)
                with gzip.open(root/batch['query_file'],'rb') as stream:self.assertEqual(hashlib.sha256(stream.read()).hexdigest(),batch['query_sha256'])
            self.assertFalse(list(root.rglob('*_candidates_blast.tsv')))

    def test_failure_keeps_partial_output_and_never_reuses_it_on_next_attempt(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            with self.assertRaisesRegex(RuntimeError,'SIGKILL'):self.run_fixture(root,fail_at=2,databases={'human':'human_fixture'})
            old_manifest_path=next(root.rglob('search_manifest.json'))
            old=json.loads(old_manifest_path.read_text())
            self.assertEqual(old['status'],'INCOMPLETE');self.assertEqual(old['completed_batches'],1)
            self.assertEqual(old['failed_batch'],2);self.assertEqual(len(list(root.rglob('*_complete.json'))),1)
            interrupted=root/old['batches'][1]['raw_file'];digest=engine.file_sha(interrupted)
            # The next call uses a new scratch store, as every real run does.
            (root/'hits.sqlite').rename(root/'old_hits.sqlite')
            runs,_,calls,_=self.run_fixture(root,databases={'human':'human_fixture'})
            self.assertEqual(len(calls),3)
            self.assertNotEqual(root/runs[0]['search_manifest'],old_manifest_path)
            self.assertEqual(engine.file_sha(interrupted),digest)
            self.assertEqual(json.loads(old_manifest_path.read_text())['status'],'INCOMPLETE')

    def test_full_audit_fresh_attempt_does_not_replace_prior_archives(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            first,_,_,_=self.run_fixture(root,count=1,detail='full_audit')
            original={batch['raw_file']:engine.file_sha(root/batch['raw_file']) for run in first for batch in run['batches']}
            (root/'hits.sqlite').rename(root/'old_hits.sqlite')
            second,_,_,_=self.run_fixture(root,count=1,detail='full_audit')
            later={batch['raw_file'] for run in second for batch in run['batches']}
            self.assertFalse(set(original)&later)
            self.assertEqual(original,{path:engine.file_sha(root/path) for path in original})

    def test_tsv_parsing_does_not_use_whole_file_reads(self):
        original_read_bytes=Path.read_bytes;original_read_text=Path.read_text
        def safe_bytes(path,*args,**kwargs):
            self.assertNotEqual(path.suffix,'.tsv','BLAST TSV must be streamed')
            return original_read_bytes(path,*args,**kwargs)
        def safe_text(path,*args,**kwargs):
            self.assertNotEqual(path.suffix,'.tsv','BLAST TSV must be streamed')
            return original_read_text(path,*args,**kwargs)
        with tempfile.TemporaryDirectory() as temp,patch.object(Path,'read_bytes',safe_bytes),patch.object(Path,'read_text',safe_text):
            self.run_fixture(Path(temp))

    def test_cancel_before_batch_does_not_launch_or_mark_complete(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);candidates,refs,bins=fixture_candidates(1);cancel=threading.Event();cancel.set()
            with patch.object(engine.shutil,'which',return_value='/fixture/blastn'),patch.object(engine.subprocess,'check_output',return_value='fixture'),patch.object(engine,'run_process') as launch:
                with self.assertRaises(engine.Cancelled):engine.blast_candidates(candidates,refs,bins,[],{'human':'fixture'},engine.Policy(),root,lambda *a,**k:None,cancel=cancel)
                launch.assert_not_called()
            self.assertEqual(json.loads(next(root.rglob('search_manifest.json')).read_text())['status'],'INCOMPLETE')
            self.assertFalse(list(root.rglob('*_complete.json')))

    def test_cancel_after_empty_subprocess_output_is_still_incomplete(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);candidates,refs,bins=fixture_candidates(1);cancel=threading.Event()
            def fake(command,outlog,unused_cancel=None,resource_guard=None,*,env_overrides=None):
                self.assertEqual(env_overrides,expected_batch_environment(command))
                Path(command[command.index('-out')+1]).write_text('');cancel.set()
            with patch.object(engine.shutil,'which',return_value='/fixture/blastn'),patch.object(engine.subprocess,'check_output',return_value='fixture'),patch.object(engine,'run_process',side_effect=fake):
                with self.assertRaises(engine.Cancelled):engine.blast_candidates(candidates,refs,bins,[],{'human':'fixture'},engine.Policy(),root,lambda *a,**k:None,cancel=cancel)
            self.assertEqual(json.loads(next(root.rglob('search_manifest.json')).read_text())['status'],'INCOMPLETE')
            self.assertFalse(list(root.rglob('*_complete.json')))


@unittest.skipUnless(shutil.which('blastn') and shutil.which('makeblastdb'),'Requires existing BLAST executables')
class RealBatchParityTests(unittest.TestCase):
    def test_tiny_real_batched_matches_single_for_genome_and_rna_with_multiple_same_chromosome_sites(self):
        candidate,ref,bins=inputs()
        genome=ref['sequence']+'N'*90+candidate['arm']+'N'*80+candidate['full']
        rna=engine.rc(candidate['full'])+'N'*80+engine.rc(candidate['arm'])
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);databases={}
            for source,seq,ident in [('human',genome,'NC_000007.14'),('rna_human',rna,'NM_FIXTURE.1')]:
                fasta=root/(source+'.fa');fasta.write_text('>'+ident+'\n'+seq+'\n')
                subprocess.run(['makeblastdb','-in',str(fasta),'-dbtype','nucl','-parse_seqids','-blastdb_version','4','-out',str(root/source)],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
                databases[source]=str(root/source)
            results=[]
            for label,batch_size in [('batched',engine.BLAST_BATCH_CANDIDATES),('single',1000)]:
                folder=root/label;folder.mkdir();candidates,refs,bins=fixture_candidates(11)
                with patch.object(engine,'BLAST_BATCH_CANDIDATES',batch_size):
                    runs=engine.blast_candidates(candidates,refs,bins,[],databases,engine.Policy(),folder,lambda *a,**k:None)
                hits=[]
                for c in candidates:
                    engine.summarize_specificity(c,engine.Policy())
                    hits.append(dict(hits=sorted(c['offtarget_hits'],key=lambda h:json.dumps(h,sort_keys=True)),specificity=c['specificity'],reasons=c['reject_reasons']))
                results.append((runs,hits))
            self.assertEqual(results[0][1],results[1][1])
            hits=results[0][1][0]['hits']
            self.assertTrue(any(h['source']=='human' and h['allowed'] for h in hits))
            self.assertGreaterEqual(len({h['start'] for h in hits if h['source']=='human' and not h['allowed'] and h['query_coverage']==100 and not h['full_primary']}),2)
            self.assertTrue(any(h['source']=='rna_human' and h['full_primary'] and h['query_coverage']==100 for h in hits))
            for run in results[0][0]:self.assertEqual(run['batch_count'],3)


if __name__=='__main__':unittest.main()
