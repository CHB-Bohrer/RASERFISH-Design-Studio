"""Recovery is status-only; tiny HTTP servers and owned subprocesses, no design."""
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import threading
import time
import unittest
from unittest.mock import patch
import urllib.error
import urllib.request

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import studio
BASE=Path(__file__).resolve().parents[1]

class PersistedStatusTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.jobs=Path(self.temp.name)/'runs';self.jobs.mkdir()
        self.ident='20261016_111134_cb49ce';self.root=self.jobs/self.ident;self.root.mkdir()
        self.state={};self.patches=[patch.object(studio,'JOBS',self.jobs),patch.object(studio,'STATE',self.state)]
        for p in self.patches:p.start()
        self.server=studio.ThreadingHTTPServer(('127.0.0.1',0),studio.Handler)
        self.thread=threading.Thread(target=self.server.serve_forever,daemon=True);self.thread.start()
    def tearDown(self):
        self.server.shutdown();self.server.server_close();self.thread.join()
        for p in reversed(self.patches):p.stop()
        self.temp.cleanup()
    def write(self,name,record):
        (self.root/name).write_text(json.dumps(record))
    def get(self,path):
        try:
            with urllib.request.urlopen('http://127.0.0.1:'+str(self.server.server_port)+path,timeout=3) as r:return r.status,json.load(r)
        except urllib.error.HTTPError as r:return r.code,json.load(r)
    def test_latest_orphan_fixture_recovers_actual_incomplete_without_state_or_download(self):
        error='Design supervisor exited before its worker; analysis is incomplete.'
        self.write('analysis_status.json',dict(status='incomplete',analysis_complete=False,error=error))
        self.write('worker_status.json',dict(status='incomplete',error=error))
        (self.root/'run.log').write_text('BLAST: human; batch 348/600 (1 candidates).\n')
        before={p.name:p.read_bytes() for p in self.root.iterdir()}
        status,body=self.get('/api/job?id='+self.ident)
        self.assertEqual(status,200);self.assertEqual(body['status'],'incomplete');self.assertEqual(body['error'],error)
        self.assertFalse(body['analysis_complete']);self.assertTrue(body['recovered_from_disk'])
        self.assertIn('348/600',body['messages'][0]);self.assertEqual(self.state,{})
        self.assertEqual(before,{p.name:p.read_bytes() for p in self.root.iterdir()})
        for endpoint in ('/download','/api/result'):
            status,body=self.get(endpoint+'?id='+self.ident);self.assertEqual(status,400)
    def test_completion_claim_without_live_export_transaction_is_not_recovered(self):
        for name in ('analysis_status.json','worker_status.json','job_status.json'):
            self.write(name,dict(status='complete',analysis_complete=True))
        (self.jobs/(self.ident+'.zip')).write_bytes(b'not a verified archive')
        code,body=self.get('/api/job?id='+self.ident)
        self.assertEqual(code,400);self.assertIn('saved results archive',body['error'])
        self.assertEqual(self.get('/download?id='+self.ident)[0],400)
    def test_incomplete_evidence_overrides_other_completion_claims(self):
        self.write('analysis_status.json',dict(status='incomplete',error='Export did not finish'))
        self.write('worker_status.json',dict(status='complete'))
        self.assertEqual(self.get('/api/job?id='+self.ident)[1]['status'],'incomplete')
    def test_stale_running_is_unverified_not_complete_or_resumed(self):
        self.write('analysis_status.json',dict(status='running',analysis_complete=False))
        code,body=self.get('/api/job?id='+self.ident)
        self.assertEqual(code,400);self.assertIn('unverified',body['error']);self.assertEqual(self.state,{})
    def test_oversized_status_is_rejected_and_log_tail_is_bounded(self):
        self.write('analysis_status.json',dict(status='incomplete',error='Stopped'))
        with (self.root/'run.log').open('wb') as f:
            for _ in range(20):f.write(b'A'*65536)
            f.write(b'\nFinal progress\n')
        body=self.get('/api/job?id='+self.ident)[1]
        self.assertLess(sum(len(x) for x in body['messages']),70000)
        (self.root/'analysis_status.json').write_bytes(b' '*1048577)
        self.assertEqual(self.get('/api/job?id='+self.ident)[0],400)
    def test_malformed_mismatched_and_symlinked_status_rejected(self):
        for content in ('[]','not json','{"status":"incomplete","job_id":"elsewhere"}'):
            (self.root/'analysis_status.json').write_text(content)
            self.assertEqual(self.get('/api/job?id='+self.ident)[0],400)
        (self.root/'analysis_status.json').unlink()
        outside=self.jobs/'outside.json';outside.write_text('{"status":"incomplete"}')
        (self.root/'analysis_status.json').symlink_to(outside)
        self.assertEqual(self.get('/api/job?id='+self.ident)[0],400)
        self.assertEqual(self.get('/api/job?id=../outside')[0],400)
    def test_server_identity_contains_no_token(self):
        code,body=self.get('/api/server')
        self.assertEqual(code,200);self.assertEqual(body['pid'],os.getpid());self.assertEqual(body['port'],self.server.server_port)
        self.assertEqual(set(body),{'version','pid','instance_id','app_path','port'})
    def test_live_state_takes_precedence(self):
        self.write('analysis_status.json',dict(status='incomplete',error='old fixture'))
        self.state[self.ident]=dict(id=self.ident,status='running',messages=['new live fixture'],cancel=threading.Event())
        code,body=self.get('/api/job?id='+self.ident)
        self.assertEqual(body['status'],'running');self.assertNotIn('cancel',body)

SERVER_FIXTURE=r'''
import json,os,subprocess,sys,threading,time
from pathlib import Path
import studio
root=Path(sys.argv[1]);studio.BASE=root;studio.JOBS=root/'runs';studio.JOBS.mkdir()
job=studio.JOBS/'owned_fixture';job.mkdir()
worker=subprocess.Popen([sys.executable,'-c','import time; time.sleep(30)'],start_new_session=True)
worker._raserfish_owned_pgid=worker.pid
item=dict(id='owned_fixture',status='running',cancel=threading.Event(),process=worker,messages=[])
studio.STATE['owned_fixture']=item
(root/'worker.pid').write_text(str(worker.pid))
def own_worker():
    while not item['cancel'].wait(.02):
        if worker.poll() is not None:return
    studio.stop_worker(worker,job,reason=item['cancel_reason'])
    item['status']='incomplete'
    studio.analysis_status(job,'incomplete',error=item['cancel_reason'])
thread=threading.Thread(target=own_worker);item['thread']=thread;thread.start()
raise SystemExit(studio.serve_studio(0,True))
'''

class ServerLifecycleTests(unittest.TestCase):
    def test_shutdown_record_write_failure_still_stops_only_owned_worker(self):
        injection="""
def fail_recovery_record(reason):
    raise OSError('Synthetic recovery record write failure')
studio.stop_owned_recovery_chains=fail_recovery_record
"""
        fixture=SERVER_FIXTURE.replace('raise SystemExit(studio.serve_studio(0,True))',
            injection+'\nraise SystemExit(studio.serve_studio(0,True))')
        unrelated=subprocess.Popen([sys.executable,'-c','import time; time.sleep(40)'])
        try:
            with tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp)
                with (root/'console.log').open('w') as output:
                    proc=subprocess.Popen([sys.executable,'-c',fixture,tmp],cwd=BASE,
                        stdout=output,stderr=output,start_new_session=True)
                    worker=None
                    try:
                        journal=root/'logs/server_lifecycle.jsonl';deadline=time.monotonic()+8
                        while time.monotonic()<deadline:
                            if journal.exists() and 'SERVER_STARTED' in journal.read_text():break
                            if proc.poll() is not None:self.fail((root/'console.log').read_text())
                            time.sleep(.02)
                        self.assertIn('SERVER_STARTED',journal.read_text())
                        worker=int((root/'worker.pid').read_text())
                        os.kill(proc.pid,signal.SIGTERM)
                        self.assertEqual(proc.wait(timeout=12),128+signal.SIGTERM)
                        events=[json.loads(line) for line in journal.read_text().splitlines()]
                        failure=next(event for event in events if event['event']=='RECOVERY_SHUTDOWN_RECORD_FAILED')
                        self.assertIn('Synthetic recovery record write failure',failure['error'])
                        self.assertEqual(events[-1]['event'],'SERVER_STOPPED')
                        self.assertEqual(events[-1]['unfinished_supervisor_threads'],[])
                        status=json.loads((root/'runs/owned_fixture/analysis_status.json').read_text())
                        self.assertEqual(status['status'],'incomplete')
                        self.assertIn('SIGTERM',status['error'])
                        with self.assertRaises(ProcessLookupError):os.kill(worker,0)
                        worker=None
                        self.assertIsNone(unrelated.poll())
                    finally:
                        if proc.poll() is None:proc.kill();proc.wait()
                        # If the regression recurs, clean up only this fixture's
                        # explicitly recorded child process group.
                        if worker is not None:
                            try:os.killpg(worker,signal.SIGKILL)
                            except ProcessLookupError:pass
        finally:unrelated.terminate();unrelated.wait()

    def test_catchable_signals_are_journalled_and_stop_only_owned_worker(self):
        unrelated=subprocess.Popen([sys.executable,'-c','import time; time.sleep(40)'])
        try:
            for sig in (signal.SIGTERM,signal.SIGHUP,signal.SIGINT):
                with self.subTest(signal=sig),tempfile.TemporaryDirectory() as tmp:
                    root=Path(tmp)
                    with (root/'console.log').open('w') as output:
                        proc=subprocess.Popen([sys.executable,'-c',SERVER_FIXTURE,tmp],cwd=BASE,stdout=output,stderr=output,start_new_session=True)
                        try:
                            journal=root/'logs/server_lifecycle.jsonl';deadline=time.monotonic()+8
                            while time.monotonic()<deadline:
                                if journal.exists() and 'SERVER_STARTED' in journal.read_text():break
                                if proc.poll() is not None:self.fail((root/'console.log').read_text())
                                time.sleep(.02)
                            self.assertIn('SERVER_STARTED',journal.read_text())
                            worker=int((root/'worker.pid').read_text())
                            os.kill(proc.pid,sig);self.assertEqual(proc.wait(timeout=12),128+sig)
                            events=[json.loads(x) for x in journal.read_text().splitlines()]
                            received=next(e for e in events if e['event']=='SIGNAL_RECEIVED')
                            self.assertEqual(received['signal'],signal.Signals(sig).name)
                            self.assertEqual(received['pid'],proc.pid);self.assertEqual(received['session_id'],proc.pid)
                            self.assertEqual(events[-1]['event'],'SERVER_STOPPED');self.assertEqual(events[-1]['unfinished_supervisor_threads'],[])
                            status=json.loads((root/'runs/owned_fixture/analysis_status.json').read_text())
                            self.assertEqual(status['status'],'incomplete');self.assertIn(signal.Signals(sig).name,status['error'])
                            with self.assertRaises(ProcessLookupError):os.kill(worker,0)
                            self.assertIsNone(unrelated.poll())
                        finally:
                            if proc.poll() is None:proc.kill();proc.wait()
        finally:unrelated.terminate();unrelated.wait()
    def test_normal_return_and_unhandled_exception_are_durable(self):
        for error in (None,RuntimeError('Synthetic server failure')):
            with self.subTest(error=error),tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp)
                class FixtureServer:
                    server_port=12345
                    def __init__(self,*a):pass
                    def serve_forever(self):
                        if error:raise error
                    def server_close(self):pass
                with patch.object(studio,'BASE',root),patch.object(studio,'STATE',{}),patch.object(studio,'ThreadingHTTPServer',FixtureServer):
                    if error:
                        with self.assertRaisesRegex(RuntimeError,'Synthetic'):studio.serve_studio(0,True)
                    else:self.assertEqual(studio.serve_studio(0,True),0)
                events=[json.loads(x) for x in (root/'logs/server_lifecycle.jsonl').read_text().splitlines()]
                self.assertEqual(events[-1]['event'],'SERVER_STOPPED')
                self.assertEqual(events[-1]['exit_code'],1 if error else 0)
                self.assertEqual(any(e['event']=='UNHANDLED_EXCEPTION' for e in events),bool(error))

if __name__=='__main__':unittest.main()
