"""v1.6.1 scope and bounded-memory regressions, plus a tiny real-engine run.

Synthetic references test implementation, not human-genome suitability.
"""
import copy
import csv
import gzip
import json
import os
import random
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import design_engine as engine
from design_engine import Folding, Policy, design, select_candidates, summarize_specificity
from specificity_store import HitStore
from test_design import cand, hit, IMS
from test_spatial_search import candidate, region


class CrossPrimaryScopeTests(unittest.TestCase):
    def test_default_off_and_saved_thread_request_is_capped(self):
        self.assertFalse(Policy().screen_cross_primary)
        self.assertEqual(Policy().threads,2)
        self.assertEqual(Policy.from_dict({'threads':4}).threads,2)
        with self.assertRaises(ValueError):Policy.from_dict({'screen_cross_primary':'false'})
        with self.assertRaises(ValueError):Policy.from_dict({'threads':True})
        with self.assertRaises(ValueError):Policy(threads=3).validate()

    def test_disabled_pair_function_refuses_before_stem_or_native_call(self):
        f=Folding(Policy(mode='local_preview'))
        with patch('design_engine.exact_pair_segments',side_effect=AssertionError('must not run')):
            with self.assertRaisesRegex(RuntimeError,'NOT RUN - USER DISABLED'):
                f.pair(candidate('a',1),candidate('b',100))
        self.assertEqual(dict(f.evaluated_counts),{})

    def test_off_skips_cross_calls_in_selection_repair_and_final_checks(self):
        calls=[]
        class SelfOnly:
            def pair(self,a,b):
                if a['id']!=b['id']:raise AssertionError('cross calculation was called while disabled')
                calls.append(a['id'])
                return {'blocked':a['id']=='self_bad'}
        bins=[region(500),region(500,start=501,number=2)]
        cs=[candidate('a',40),candidate('alternate',42),candidate('b',130),
            candidate('self_bad',240),candidate('boundary_left',465),
            candidate('boundary_right',501,'X_B002'),candidate('d',650,'X_B002')]
        cs.append(candidate('individual_bad',780,'X_B002'));cs[-1]['reject_reasons']=['Individual accessibility below threshold']
        reports=[]
        out=select_candidates(cs,bins,Policy(mode='local_preview'),SelfOnly(),lambda *a,**k:None,selection_report=reports)
        self.assertGreater(len(out),2)
        self.assertEqual({c['bin_id'] for c in out},{'X_B001','X_B002'})
        self.assertNotIn('self_bad',{c['id'] for c in out})
        self.assertNotIn('individual_bad',{c['id'] for c in out})
        self.assertFalse({'boundary_left','boundary_right'}<={c['id'] for c in out})
        self.assertGreater(next(r for r in reports if r['chosen'])['repair_seed_checks'],0)
        self.assertTrue(calls)

    def test_enabled_cross_checks_still_change_selection(self):
        calls=[]
        class BlockCross:
            def pair(self,a,b):
                calls.append((a['id'],b['id']))
                return {'blocked':a['id']!=b['id']}
        out=select_candidates([candidate('a',40),candidate('b',150)],[region(500)],
            Policy(mode='local_preview',screen_cross_primary=True),BlockCross(),lambda *a,**k:None)
        self.assertEqual(len(out),1)
        self.assertTrue(any(a!=b for a,b in calls))

    def test_self_exact_stem_threshold_is_retained(self):
        f=Folding(Policy(mode='local_preview'))
        a=candidate('self',10);a['full']='A'*12+'T'*12
        result=f.pair(a,a)
        self.assertTrue(result['blocked'])
        self.assertEqual(result['scope'],'SELF_DIMER')
        self.assertGreater(result['max_exact_stem'],11)


class DiskSpecificityTests(unittest.TestCase):
    def synthetic_hits(self):
        rng=random.Random(482)
        rows=[]
        for i in range(6000):
            source=rng.choice(['human','targets','rna_human','rna_local','construct'])
            start=rng.randrange(0,50)*100+rng.randrange(35)
            rows.append(dict(hit(source=source,reference=rng.choice(['chr7','ref|NC_000007.14|','chr8','NT_123.1']),
                start=start,end=start+rng.randrange(10,75),bound_strand=rng.choice(['+','-','RNA']),allowed=rng.random()<.1),
                review=rng.random()<.8,high_tm=rng.random()<.08,moderate_tm=rng.random()<.3,
                unknown_long=rng.random()<.04,strong_rna=source.startswith('rna') and rng.random()<.2,
                rna_exact_risk=source.startswith('rna') and rng.random()<.3))
        return rows

    def test_disk_clusters_and_every_gate_match_existing_list_algorithm(self):
        rows=self.synthetic_hits()
        with tempfile.TemporaryDirectory() as tmp:
            store=HitStore(Path(tmp)/'hits.sqlite')
            try:
                disk=store.hits('c')
                for h in rows:disk.append(h)
                store.flush()
                for human in ('intended_only','report_only'):
                    for rna in (False,True):
                        policy=Policy(human_policy=human,avoid_gene_rna=rna)
                        a=dict(offtarget_hits=rows,reject_reasons=[])
                        b=dict(offtarget_hits=disk,reject_reasons=[])
                        summarize_specificity(a,policy);summarize_specificity(b,policy)
                        self.assertEqual(a,b|{'offtarget_hits':rows})
                self.assertEqual(len(disk),len(rows))
                self.assertEqual(list(disk),rows)
                self.assertEqual(store.connection.execute('PRAGMA cache_size').fetchone()[0],-4096)
                self.assertEqual(store.connection.execute('PRAGMA temp_store').fetchone()[0],1)
            finally:store.close()

    def test_cap_detection_uses_distinct_subjects_and_does_not_truncate_hits(self):
        with tempfile.TemporaryDirectory() as tmp:
            store=HitStore(Path(tmp)/'hits.sqlite')
            try:
                for subject in ('chr1','chr1','chr2','chr3'):store.record_subject('human','c',subject)
                self.assertEqual(store.capped_candidates('human',3),['c'])
                self.assertEqual(store.capped_candidates('human',4),[])
                self.assertEqual(store.capped_candidates('rna_human',1),[])
            finally:store.close()

    def test_duplex_cache_is_bounded_without_changing_results(self):
        self.assertEqual(engine.duplex_estimate.cache_parameters()['maxsize'],8192)
        q='GCGCATGC'*5
        first=engine.duplex_estimate(q,q,Policy().profile)
        engine.duplex_estimate.cache_clear()
        self.assertEqual(first,engine.duplex_estimate(q,q,Policy().profile))

    def test_fold_result_caches_are_bounded_and_eviction_keeps_same_answer(self):
        fold=Folding(Policy(mode='local_preview'))
        first=candidate('c0',10);first['full']='ACGTACGT'
        expected=fold.pair(first,first)
        for n in range(1,4200):
            c=candidate('c'+str(n),10)
            c['full']=''.join('ACGT'[(n//(4**k))%4] for k in range(7))
            fold.fold(c);fold.pair(c,c)
        self.assertLessEqual(len(fold.monomers),1024)
        self.assertLessEqual(len(fold.pairs),4096)
        self.assertEqual(fold.pair(first,first),expected)


class ProcessCleanupTests(unittest.TestCase):
    def test_memory_guard_failure_terminates_and_reaps_native_child(self):
        created=[];original=subprocess.Popen
        class Guard:
            def checkpoint(self,*args,**kwargs):raise RuntimeError('Conservative memory stop test')
        def launch(*args,**kwargs):
            child=original(*args,**kwargs);created.append(child);return child
        with tempfile.TemporaryDirectory() as tmp,patch('design_engine.subprocess.Popen',launch),patch.dict(os.environ,{'RASERFISH_DESIGN_WORKER':'0'}):
            with self.assertRaisesRegex(RuntimeError,'memory stop'):
                engine.run_process([engine.sys.executable,'-c','import time; time.sleep(10)'],Path(tmp)/'native.log',resource_guard=Guard())
        self.assertEqual(len(created),1)
        self.assertIsNotNone(created[0].returncode)
        with self.assertRaises(ProcessLookupError):os.kill(created[0].pid,0)

    def test_worker_child_keeps_supervised_group_and_inherited_lock(self):
        with tempfile.TemporaryDirectory() as tmp,patch('design_engine.subprocess.Popen') as launch,patch.dict(os.environ,{'RASERFISH_DESIGN_WORKER':'1','RASERFISH_DESIGN_LOCK_FD':'123'}):
            launch.return_value.poll.return_value=0
            launch.return_value.returncode=0
            engine.run_process(['test-tool'],Path(tmp)/'native.log')
        self.assertFalse(launch.call_args.kwargs['start_new_session'])
        self.assertEqual(launch.call_args.kwargs['pass_fds'],(123,))


@unittest.skipUnless(shutil.which('blastn') and shutil.which('makeblastdb'),'Requires installed BLAST+ and ViennaRNA')
class RealEngineDisabledTests(unittest.TestCase):
    def test_standard_streaming_matches_full_audit_with_retained_native_qc(self):
        import RNA
        from compact_outputs import concise_saved_result
        rng=random.Random(771)
        sequence=''.join(rng.choice('ACGT') for _ in range(800))
        reference=dict(id='TINY',sequence=sequence,length=len(sequence),features=[],
            mapping=dict(verified=True,start=1,end=len(sequence),chromosome='chr7',orientation='+',assembly='GRCh38'))
        config=dict(title='Tiny Standard streaming parity fixture, not a genomic validation',references=[reference],imagers=IMS,
            bin_width=400,tail_policy='strict',assignments={'TINY_B001':'DP1','TINY_B002':'DP2'},
            policy=dict(shortlist_per_bin=4,extension_shortlist_per_bin=2,threads=2,screen_cross_primary=False))
        original_pair=Folding.pair;original_fc=RNA.fold_compound;original_stems=engine.exact_pair_segments
        original_search=engine.blast_candidates;search_evidence={}
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            for source,seq in [('human',sequence),('rna_human','A'*120)]:
                fasta=root/(source+'.fa');fasta.write_text(('>NC_000007.14\n' if source=='human' else '>NM_SYNTHETIC.1\n')+seq+'\n')
                subprocess.run(['makeblastdb','-in',str(fasta),'-dbtype','nucl','-parse_seqids','-blastdb_version','4','-out',str(root/source)],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
            config['databases']={source:str(root/source) for source in ('human','rna_human')}
            results={}
            for detail in ('full_audit','standard'):
                pair_calls=[];native_sequences=[];stem_calls=[];guard_stages=[]
                def pair_spy(fold,a,b):
                    self.assertEqual(a['id'],b['id'],'Disabled cross-primary calculation ran')
                    pair_calls.append((a['id'],b['id']));return original_pair(fold,a,b)
                def native_spy(seq,*args,**kwargs):
                    native_sequences.append(seq)
                    if '&' in seq:
                        a,b=seq.split('&');self.assertEqual(a,b,'Disabled cross-primary cofold ran')
                    return original_fc(seq,*args,**kwargs)
                def stem_spy(a,b,*args,**kwargs):
                    self.assertEqual(a,b,'Disabled cross-primary exact-stem calculation ran')
                    stem_calls.append((a,b));return original_stems(a,b,*args,**kwargs)
                class Guard:
                    def checkpoint(self,stage,force=False):guard_stages.append(stage)
                def search_spy(*args,**kwargs):
                    value=original_search(*args,**kwargs)
                    search_evidence[detail]=copy.deepcopy(value)
                    return value
                current=copy.deepcopy(config);current['output_detail']=detail
                kwargs={'checkpoint_dir':root/'standard_checkpoints'} if detail=='standard' else {}
                with patch.object(Folding,'pair',pair_spy),patch.object(RNA,'fold_compound',native_spy),patch.object(engine,'exact_pair_segments',stem_spy),patch.object(engine,'blast_candidates',search_spy):
                    result=design(current,root/detail,lambda *a,**k:None,resource_guard=Guard(),**kwargs)
                count=result['counts']['candidates']
                self.assertGreater(count,0);self.assertGreater(result['counts']['selected'],0)
                self.assertEqual(sum('&' not in seq for seq in native_sequences),count)
                self.assertEqual(sum('&' in seq for seq in native_sequences),count)
                self.assertTrue(pair_calls and stem_calls)
                self.assertEqual(result['counts']['cross_primary_evaluations'],0)
                self.assertEqual(result['cross_primary_screening']['within_region_status'],'NOT RUN - USER DISABLED')
                self.assertEqual(result['cross_primary_screening']['between_region_status'],'NOT RUN - USER DISABLED')
                self.assertEqual(result['self_dimer_screening']['status'],'RUN_AQUEOUS_MODEL')
                self.assertEqual(result['order_status'],'DRAFT')
                self.assertEqual({run['source'] for run in result['blast_runs']},{'human','rna_human'})
                self.assertIn('Individual-primary folding',guard_stages)
                self.assertIn('BLAST subprocess',guard_stages)
                self.assertIn('Final DRAFT export',guard_stages)
                self.assertFalse(list((root/detail).rglob('*_candidates_blast.tsv')))
                results[detail]=result if detail=='standard' else concise_saved_result(result)
            for key in ('selected_ids','counts','candidates','selection_search','cross_primary_screening','self_dimer_screening'):
                self.assertEqual(results['standard'][key],results['full_audit'][key],key)
            for actual,expected in zip(search_evidence['standard'],search_evidence['full_audit']):
                for key in ('source','raw_alignments','retained_alignments','capped_candidate_ids'):
                    self.assertEqual(actual[key],expected[key],key)
                self.assertEqual([(b['output_sha256'],b['output_bytes']) for b in actual['batches']],
                                 [(b['output_sha256'],b['output_bytes']) for b in expected['batches']])
            self.assertTrue((root/'standard_checkpoints'/'specificity_scores.sqlite').is_file())

    def test_tiny_full_engine_retains_native_monomer_self_and_both_reference_searches(self):
        import RNA
        rng=random.Random(771)
        sequence=''.join(rng.choice('ACGT') for _ in range(800))
        reference=dict(id='TINY',sequence=sequence,length=len(sequence),features=[],
            mapping=dict(verified=True,start=1,end=len(sequence),chromosome='chr7',orientation='+',assembly='GRCh38'))
        config=dict(title='Tiny real-engine implementation test, not a genome validation',references=[reference],imagers=IMS,
            bin_width=400,tail_policy='strict',assignments={'TINY_B001':'DP1','TINY_B002':'DP2'},
            policy=dict(shortlist_per_bin=4,extension_shortlist_per_bin=2,threads=2,screen_cross_primary=False),
            output_detail='full_audit')
        pair_calls=[];native_sequences=[];stem_calls=[];guard_stages=[]
        original_pair=Folding.pair;original_fc=RNA.fold_compound;original_stems=engine.exact_pair_segments
        def pair_spy(fold,a,b):
            pair_calls.append((a['id'],b['id']))
            if a['id']!=b['id']:raise AssertionError('Disabled cross-primary computation ran')
            return original_pair(fold,a,b)
        def native_spy(sequence,*args,**kwargs):
            native_sequences.append(sequence)
            if '&' in sequence:
                left,right=sequence.split('&');self.assertEqual(left,right)
            return original_fc(sequence,*args,**kwargs)
        def stem_spy(a,b,*args,**kwargs):
            stem_calls.append((a,b));self.assertEqual(a,b)
            return original_stems(a,b,*args,**kwargs)
        class Guard:
            def checkpoint(self,stage,force=False):guard_stages.append(stage)
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            for source,seq in [('human',sequence),('rna_human','A'*120)]:
                fasta=root/(source+'.fa');fasta.write_text(('>NC_000007.14\n' if source=='human' else '>NM_SYNTHETIC.1\n')+seq+'\n')
                subprocess.run(['makeblastdb','-in',str(fasta),'-dbtype','nucl','-parse_seqids','-blastdb_version','4','-out',str(root/source)],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
            config['databases']={k:str(root/k) for k in ('human','rna_human')}
            with patch.object(Folding,'pair',pair_spy),patch.object(RNA,'fold_compound',native_spy),patch('design_engine.exact_pair_segments',stem_spy):
                result=design(copy.deepcopy(config),root/'result',lambda *a,**k:None,resource_guard=Guard())
            self.assertGreater(result['counts']['candidates'],0)
            self.assertGreater(result['counts']['selected'],0)
            self.assertEqual(result['counts']['cross_primary_evaluations'],0)
            self.assertGreater(result['counts']['monomer_evaluations'],0)
            self.assertGreater(result['counts']['self_dimer_evaluations'],0)
            self.assertTrue(pair_calls and stem_calls)
            self.assertTrue(any('&' not in s for s in native_sequences))
            self.assertTrue(any('&' in s for s in native_sequences))
            self.assertEqual({r['source'] for r in result['blast_runs']},{'human','rna_human'})
            for run in result['blast_runs']:
                self.assertEqual(run['command'][run['command'].index('-num_threads')+1],'2')
            self.assertEqual(result['cross_primary_screening']['status'],'NOT RUN - USER DISABLED')
            self.assertEqual(result['cross_primary_screening']['within_region_status'],'NOT RUN - USER DISABLED')
            self.assertEqual(result['cross_primary_screening']['between_region_status'],'NOT RUN - USER DISABLED')
            self.assertEqual(result['self_dimer_screening']['status'],'RUN_AQUEOUS_MODEL')
            self.assertEqual(result['order_status'],'DRAFT')
            self.assertIn('PROVISIONAL_DRAFT',result['status'])
            self.assertFalse(result['all_scientific_checks_complete'])
            self.assertTrue((root/'result'/'DRAFT_IDT_oPools.csv').is_file())
            with gzip.open(root/'result'/'debug'/'all_selected_pair_interactions.csv.gz','rt') as stream:
                rows=list(csv.DictReader(stream));self.assertTrue(rows)
                self.assertTrue(all(row['a']==row['b'] and row['scope']=='SELF_DIMER' for row in rows))
            self.assertIn('Individual-primary folding',guard_stages)
            self.assertIn('BLAST subprocess',guard_stages)
            self.assertIn('Final DRAFT export',guard_stages)


if __name__=='__main__':unittest.main()
