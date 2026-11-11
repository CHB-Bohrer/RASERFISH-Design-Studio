"""Checkpoint discovery/resume HTTP and file boundaries; no scientific run starts."""
import copy
import hashlib
import io
import json
import os
from pathlib import Path
import select
import subprocess
import sys
import tempfile
import threading
import unittest
from unittest.mock import patch
import urllib.error
import urllib.request
import zipfile

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import design_checkpoints
import runtime_safety
import studio
from compact_outputs import write_input_manifest
from design_engine import Policy
from sequence_tools import make_bins,validate_imagers
BASE=Path(__file__).resolve().parents[1]

def small_config():
    refs=[dict(id='target',description='Small DP3 checkpoint fixture',sequence='ACGT'*100,length=400,
               mapping=dict(chromosome='chr1',start=8186425,end=8186824),features=[])]
    imagers=validate_imagers([im for im in json.loads((BASE/'imagers.json').read_text()) if im['id']=='DP3'])
    bins,_=make_bins(refs,400,'strict',None,{'target':[1,400]})
    return dict(title='Saved DP3 fixture',output_detail='standard',references=refs,imagers=imagers,
        policy=vars(Policy.from_dict({})),bin_width=400,tail_policy='strict',bounds={'target':[1,400]},
        assignments={bins[0]['id']:'DP3'},allowed_intervals=[],transcription={},annotation_provenance=[],
        databases={'human':'/existing/human','rna_human':'/existing/rna'},
        background=[dict(id='construct',description='original saved construct',sequence='GATTACA')],
        transcripts=[dict(id='rna',description='original saved RNA',sequence='ACGTTGCA')])

class NoRunThread:
    def __init__(self,**kwargs):self.kwargs=kwargs;self.started=False
    def start(self):self.started=True

class CheckpointServerTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.base=Path(self.temp.name).resolve();self.jobs=self.base/'runs';self.jobs.mkdir()
        self.lock_path=self.base/'private_design.lock';self.state={};self.summaries={};self.server=None;self.server_thread=None
        self.patches=[patch.object(studio,'JOBS',self.jobs),patch.object(studio,'STATE',self.state),
            patch.object(runtime_safety,'default_lock_path',return_value=self.lock_path),
            patch.object(design_checkpoints,'read_checkpoint_summary',side_effect=lambda directory:copy.deepcopy(self.summaries.get(str(directory))))]
        for p in self.patches:p.start()
        self.config=small_config()
    def tearDown(self):
        if self.server:
            self.server.shutdown();self.server.server_close();self.server_thread.join(timeout=3)
        for item in self.state.values():
            lease=item.get('design_lock')
            if lease:lease.release()
        for p in reversed(self.patches):p.stop()
        self.temp.cleanup()
    def saved(self,ident='20261016_100000_origin',status='incomplete',origin=None,with_checkpoint=True):
        root=self.jobs/ident;root.mkdir()
        (root/'analysis_status.json').write_text(json.dumps(dict(status=status,analysis_complete=status=='complete',job_id=ident)))
        if origin:
            (root/'checkpoint_link.json').write_text(json.dumps(dict(schema='raserfish.checkpoint-link.v1',origin_run_id=origin,attempt_run_id=ident)))
        else:
            result=root/'result';result.mkdir();bins,_=make_bins(self.config['references'],400,'strict',self.config['assignments'],self.config['bounds'])
            write_input_manifest(self.config,bins,result,application_version='1.6.1')
            if with_checkpoint:
                directory=root/'checkpoints';directory.mkdir()
                self.summaries[str(directory)]=dict(summary='Saved DP3 fixture',completed_batches=2,total_batches=8,
                    folded_candidates=1,total_candidates=4,created_utc='2026-10-16T12:00:00Z',status='INCOMPLETE',resumable=True)
        return root
    def snapshot(self,root):
        return {p.relative_to(root).as_posix():(p.stat().st_size,p.stat().st_mtime_ns,hashlib.sha256(p.read_bytes()).hexdigest())
                for p in root.rglob('*') if p.is_file() and not p.is_symlink()}
    def start_http(self):
        if self.server:return
        self.server=studio.ThreadingHTTPServer(('127.0.0.1',0),studio.Handler)
        self.server_thread=threading.Thread(target=self.server.serve_forever,daemon=True);self.server_thread.start()
    def request(self,path,body=None,headers=None):
        self.start_http();headers=dict(headers or {})
        payload=None if body is None else json.dumps(body).encode()
        request=urllib.request.Request(f'http://127.0.0.1:{self.server.server_port}'+path,data=payload,headers=headers)
        try:
            with urllib.request.urlopen(request,timeout=4) as response:return response.status,json.load(response)
        except urllib.error.HTTPError as response:return response.code,json.load(response)
    def test_discovery_deduplicates_newest_incomplete_attempt_and_never_changes_old_files(self):
        root=self.saved();new=self.saved('20261016_120000_new',origin=root.name);before=self.snapshot(self.jobs)
        reply=studio.saved_checkpoint_runs()
        self.assertEqual([r['id'] for r in reply['runs']],[new.name]);self.assertEqual(reply['runs'][0]['origin_run_id'],root.name)
        self.assertEqual(reply['runs'][0]['status'],'incomplete');self.assertFalse(reply['runs'][0]['analysis_complete'])
        self.assertFalse(reply['active_design']);self.assertEqual(self.snapshot(self.jobs),before)
    def test_completed_newest_attempt_hides_older_and_direct_old_id_cannot_resume(self):
        root=self.saved();new=self.saved('20261016_120000_complete',origin=root.name,status='complete')
        self.assertEqual(studio.saved_checkpoint_runs()['runs'],[])
        with self.assertRaisesRegex(ValueError,'already completed.*'+new.name):studio.prepare_saved_resume(root.name)
        with self.assertRaisesRegex(ValueError,'already completed'):studio.prepare_saved_resume(new.name)
    def test_old_run_without_checkpoint_cannot_be_resumed_or_listed(self):
        root=self.saved(with_checkpoint=False)
        (root/'result/human_candidates_blast.tsv').write_text('Interrupted raw search is not a checkpoint\n')
        self.assertEqual(studio.saved_checkpoint_runs()['runs'],[])
        with self.assertRaisesRegex(ValueError,'Older interrupted BLAST files cannot be resumed'):studio.prepare_saved_resume(root.name)
        self.assertEqual(self.state,{})
    def test_unsafe_summary_is_hidden_and_resume_preserves_reason(self):
        root=self.saved();self.summaries[str(root/'checkpoints')].update(resumable=False,resume_reason='Retained-hit store is missing or unsafe.')
        self.assertEqual(studio.saved_checkpoint_runs()['runs'],[])
        with self.assertRaisesRegex(ValueError,'store is missing or unsafe'):studio.prepare_saved_resume(root.name)
    def test_saved_configuration_roundtrips_original_inputs_without_mutation(self):
        root=self.saved();attempt=self.saved('20261016_120000_new',origin=root.name);before=self.snapshot(self.jobs)
        config,origin=studio.prepare_saved_resume(attempt.name)
        self.assertEqual(config,self.config);self.assertEqual(origin,root.name);self.assertEqual(self.snapshot(self.jobs),before)
        self.assertEqual(config['imagers'][0]['id'],'DP3');self.assertFalse(config['policy']['screen_cross_primary'])
    def test_changed_compressed_input_asset_is_rejected_before_any_attempt(self):
        root=self.saved();manifest=json.loads((root/'result/run_config.json').read_text())
        asset=root/'result'/next(iter(manifest['assets']));asset.write_bytes(asset.read_bytes()+b'changed')
        with self.assertRaisesRegex(ValueError,'SHA-256 mismatch'):studio.prepare_saved_resume(root.name)
        self.assertEqual(self.state,{})
    def test_changed_manifest_scientific_settings_fail_original_configuration_hash(self):
        root=self.saved();path=root/'result/run_config.json';manifest=json.loads(path.read_text());manifest['configuration']['policy']['tm_min']=44
        path.write_text(json.dumps(manifest))
        with self.assertRaisesRegex(ValueError,'configuration SHA-256 mismatch'):studio.prepare_saved_resume(root.name)
    def test_manifest_size_and_invalid_json_rejected_with_bounded_read(self):
        root=self.saved();path=root/'result/run_config.json';path.write_text('{broken')
        with self.assertRaises(ValueError):studio.prepare_saved_resume(root.name)
        with path.open('wb') as f:f.truncate(16*1024*1024+1)
        with self.assertRaisesRegex(ValueError,'safe read limit'):studio.prepare_saved_resume(root.name)
    def test_full_audit_resume_is_rejected_to_preserve_original_raw_evidence_scope(self):
        self.config['output_detail']='full_audit';root=self.saved()
        with self.assertRaisesRegex(ValueError,'Full audit requires a fresh run'):studio.prepare_saved_resume(root.name)
    def test_run_directory_and_checkpoint_symlink_escape_are_rejected(self):
        outside=self.base/'outside';outside.mkdir();(self.jobs/'escaped').symlink_to(outside,target_is_directory=True)
        for ident in ('../outside','/tmp','escaped',None,17):
            with self.subTest(ident=ident),self.assertRaises(ValueError):studio.saved_run_root(ident)
        root=self.saved();(root/'checkpoints').rmdir();(root/'checkpoints').symlink_to(outside,target_is_directory=True)
        with self.assertRaisesRegex(ValueError,'Invalid checkpoint directory'):studio.prepare_saved_resume(root.name)
    def test_manifest_and_asset_symlink_escape_are_rejected(self):
        root=self.saved();manifest=root/'result/run_config.json';saved=manifest.read_bytes();outside=self.base/'outside_manifest.json';outside.write_bytes(saved)
        manifest.unlink();manifest.symlink_to(outside)
        with self.assertRaisesRegex(ValueError,'Invalid saved input manifest'):studio.prepare_saved_resume(root.name)
        manifest.unlink();manifest.write_bytes(saved);record=json.loads(saved);asset=root/'result'/next(iter(record['assets']))
        outside_asset=self.base/'outside.gz';outside_asset.write_bytes(asset.read_bytes());asset.unlink();asset.symlink_to(outside_asset)
        with self.assertRaisesRegex(ValueError,'escapes the run directory'):studio.prepare_saved_resume(root.name)
    def test_link_id_cannot_escape_and_small_status_read_is_bounded(self):
        root=self.saved();link=root/'checkpoint_link.json';link.write_text(json.dumps(dict(schema='raserfish.checkpoint-link.v1',origin_run_id='../outside')))
        with self.assertRaises(ValueError):studio.prepare_saved_resume(root.name)
        with link.open('wb') as f:f.truncate(1024*1024+1)
        with self.assertRaisesRegex(ValueError,'safe read limit'):studio.checkpoint_origin(root)
    def test_lock_probe_is_readonly_for_active_and_stale_metadata(self):
        self.assertIsNone(studio.active_design_snapshot());self.assertFalse(self.lock_path.exists())
        lease=runtime_safety.acquire_design_lock(dict(job_id='original',app_pid=os.getpid()))
        try:
            before=self.lock_path.read_bytes();info=self.lock_path.stat()
            active=studio.active_design_snapshot();self.assertEqual(active['job_id'],'original')
            self.assertEqual(self.lock_path.read_bytes(),before);self.assertEqual(self.lock_path.stat().st_mtime_ns,info.st_mtime_ns)
        finally:lease.release()
        self.assertIsNone(studio.active_design_snapshot());self.assertEqual(self.lock_path.read_bytes(),before)
        self.assertEqual(self.lock_path.stat().st_mtime_ns,info.st_mtime_ns)
    def test_actual_separate_process_lease_blocks_fresh_and_resumed_attempts(self):
        root=self.saved();code='''import json,sys\nfrom runtime_safety import acquire_design_lock\nlock=acquire_design_lock({'job_id':'other_instance'},path=sys.argv[1])\nprint(json.dumps(lock.metadata),flush=True)\nsys.stdin.readline()\nlock.release()\n'''
        proc=subprocess.Popen([sys.executable,'-c',code,str(self.lock_path)],cwd=BASE,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
        try:
            ready,_,_=select.select([proc.stdout],[],[],8);self.assertTrue(ready,'Lock holder did not become ready')
            self.assertEqual(json.loads(proc.stdout.readline())['job_id'],'other_instance')
            self.assertEqual(studio.active_design_snapshot()['job_id'],'other_instance')
            before=self.lock_path.read_bytes();reply=studio.saved_checkpoint_runs();self.assertTrue(reply['active_design']);self.assertEqual(reply['active_run']['id'],'other_instance')
            for origin in (None,root.name):
                with self.subTest(origin=origin),self.assertRaises(runtime_safety.DesignAlreadyRunning):studio.start_design_job(self.config,resume_from=origin)
            self.assertEqual(self.state,{});self.assertEqual(self.lock_path.read_bytes(),before)
        finally:
            if proc.poll() is None:proc.communicate('\n',timeout=5)
            if proc.poll() is None:proc.kill();proc.wait()
        self.assertEqual(proc.returncode,0)
    def test_start_resumed_attempt_creates_new_id_and_does_not_write_old_run(self):
        root=self.saved();before=self.snapshot(self.jobs);config,origin=studio.prepare_saved_resume(root.name)
        with patch.object(studio.threading,'Thread',NoRunThread):ident,restored=studio.start_design_job(config,resume_from=origin)
        self.assertNotEqual(ident,root.name);self.assertEqual(restored,self.config);self.assertEqual(config,self.config)
        self.assertEqual(self.state[ident]['checkpoint_origin'],root.name);self.assertTrue(self.state[ident]['thread'].started)
        self.assertEqual(self.state[ident]['thread'].kwargs['args'],(ident,restored));self.assertEqual({k:v for k,v in self.snapshot(self.jobs).items() if '/automatic_recovery.' not in k and '/automatic_recovery_' not in k},before)
        self.assertFalse((self.jobs/ident).exists(),'The stub must not execute run_job or scientific work')
    def test_resume_http_auth_token_host_origin_and_exact_saved_config(self):
        root=self.saved();self.start_http();before=self.snapshot(self.jobs);real_thread=threading.Thread
        def thread_factory(*args,**kwargs):
            return NoRunThread(**kwargs) if kwargs.get('target') is studio.run_job else real_thread(*args,**kwargs)
        with patch.object(studio.threading,'Thread',side_effect=thread_factory):
            for headers in ({},{'X-Studio-Token':'wrong'},{'X-Studio-Token':studio.TOKEN,'Origin':'https://unrelated.invalid'},
                            {'X-Studio-Token':studio.TOKEN,'Host':'unrelated.invalid'}):
                with self.subTest(headers=headers):self.assertEqual(self.request('/api/resume',{'id':root.name},headers)[0],403)
            self.assertEqual(self.state,{})
            code,body=self.request('/api/resume',{'id':root.name,'config':{'title':'Untrusted current edits','policy':{'tm_min':0}}},{'X-Studio-Token':studio.TOKEN})
        self.assertEqual(code,200);self.assertNotEqual(body['id'],root.name);self.assertEqual(body['config'],self.config)
        self.assertEqual(self.state[body['id']]['checkpoint_origin'],root.name);self.assertEqual({k:v for k,v in self.snapshot(self.jobs).items() if '/automatic_recovery.' not in k and '/automatic_recovery_' not in k},before)
    def test_checkpoint_listing_http_converts_unsafe_lock_to_useful_json_error(self):
        self.lock_path.write_text('{}');self.lock_path.chmod(0o644)
        code,body=self.request('/api/checkpoints');self.assertEqual(code,400);self.assertIn('Unsafe design lock file',body['error'])
    def test_result_packaging_prunes_checkpoint_tree_without_opening_or_traversing_it(self):
        root=self.saved();checkpoint=root/'checkpoints';huge=checkpoint/'retained_hits.sqlite'
        with huge.open('wb') as f:f.truncate(5*1024**3)
        nested=root/'result/checkpoints';nested.mkdir();(nested/'receipt.json').write_text('{}')
        outside=self.base/'external';outside.mkdir();(outside/'secret').write_text('must not archive')
        (root/'outside_link').symlink_to(outside,target_is_directory=True);(root/'file_link').symlink_to(outside/'secret')
        original=os.scandir
        def bounded_scandir(path):
            self.assertNotIn('checkpoints',Path(path).parts,'Packaging traversed durable working evidence')
            return original(path)
        with patch.object(studio.os,'scandir',side_effect=bounded_scandir):
            files=list(studio.result_archive_files(root))
            data=io.BytesIO()
            with zipfile.ZipFile(data,'w') as archive:
                for file in files:archive.write(file,file.relative_to(root).as_posix())
        names=[p.relative_to(root).as_posix() for p in files]
        self.assertIn('result/run_config.json',names);self.assertFalse(any('checkpoints' in Path(name).parts for name in names))
        self.assertNotIn('file_link',names);self.assertFalse(any('outside_link' in name for name in names));self.assertLess(len(data.getvalue()),100000)
        self.assertEqual(huge.stat().st_size,5*1024**3)

if __name__=='__main__':unittest.main()
