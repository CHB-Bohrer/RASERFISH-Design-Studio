"""Tiny real subprocesses exercise supervision; no scientific design is started."""
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
import zipfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import runtime_safety
import studio

BASE = Path(__file__).resolve().parents[1]
REAL_POPEN = subprocess.Popen
REAL_ZIPFILE = zipfile.ZipFile
REAL_REPLACE = Path.replace

WORKER = r'''
import json, os, subprocess, sys, time
from pathlib import Path
from runtime_safety import inherited_design_lock
root, fd, lock_path, mode = Path(sys.argv[1]), int(sys.argv[2]), sys.argv[3], sys.argv[4]
lease = inherited_design_lock(fd, lock_path)
config = json.load(sys.stdin)
(root/'result').mkdir(exist_ok=True)
(root/'result'/'DRAFT_IDT_oPools.csv').write_text('DRAFT fixture,not an order\n')
if mode == 'waiting':
    code = "import os,signal,sys,time; from pathlib import Path; signal.signal(signal.SIGTERM,signal.SIG_IGN); p=Path(sys.argv[1]); p.write_text(str(os.getpid()));\nwhile True:\n p.with_suffix('.heartbeat').write_text(str(time.monotonic())); time.sleep(.05)"
    child = subprocess.Popen([sys.executable, '-c', code, str(root/'child.pid')])
    deadline=time.monotonic()+5
    while not (root/'child.pid').exists() and time.monotonic()<deadline: time.sleep(.01)
    (root/'ready.json').write_text(json.dumps(dict(worker_pid=os.getpid(),child_pid=child.pid,worker_group=os.getpgrp(),child_group=os.getpgid(child.pid))))
    print(json.dumps(dict(message='Tiny worker and BLAST-like child ready',stage='Fixture native work')),flush=True)
    time.sleep(20)
    raise SystemExit(99)
if mode == 'crash':
    (root/'worker_status.json').write_text(json.dumps(dict(status='incomplete',error='Tiny fixture worker failed.')))
    raise SystemExit(7)
(root/'result'/'design_result.json').write_text(json.dumps(dict(counts={'selected':1},warnings=[],status='PROVISIONAL_DRAFT_CROSS_PRIMARY_DISABLED',order_status='DRAFT')))
(root/'worker_status.json').write_text(json.dumps(dict(status='complete')))
lease.release()
'''


def running(pid):
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    # A reparented zombie is already stopped, even if a container's PID 1 has not
    # reaped it yet. macOS launchd normally reaps it immediately.
    status = Path(f'/proc/{pid}/stat')
    if status.exists():
        try:
            if status.read_text().split(')', 1)[1].strip().split()[0] == 'Z':
                return False
        except (OSError, IndexError):
            pass
    return True


class WorkerSupervisionTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.jobs = Path(self.tmp.name) / 'runs'
        self.jobs.mkdir()
        self.lock_path = Path(self.tmp.name) / 'design.lock'
        self.job_id = 'fixture_job'
        self.root = self.jobs / self.job_id
        self.state = {}
        self.children = []
        self.leases = []
        self.patchers = [patch.object(studio, 'JOBS', self.jobs), patch.object(studio, 'STATE', self.state)]
        for patcher in self.patchers:
            patcher.start()

    def tearDown(self):
        for child in self.children:
            studio.stop_worker(child)
            for stream in (child.stdin, child.stdout, child.stderr):
                if stream is not None and not stream.closed:
                    stream.close()
        for lease in self.leases:
            lease.release()
        for patcher in reversed(self.patchers):
            patcher.stop()
        self.tmp.cleanup()

    def prepare(self):
        lease = runtime_safety.acquire_design_lock(dict(job_id=self.job_id), self.lock_path)
        self.leases.append(lease)
        self.state[self.job_id] = dict(id=self.job_id, status='running', stage='Starting', messages=[],
                                       cancel=threading.Event(), design_lock=lease)
        return lease

    def factory(self, mode):
        def launch(command, **kwargs):
            self.assertTrue(kwargs.get('start_new_session'))
            self.assertEqual(kwargs.get('pass_fds'), (self.state[self.job_id]['design_lock'].fd,))
            root = command[command.index('--worker-output') + 1]
            fd = command[command.index('--lock-fd') + 1]
            worker = REAL_POPEN([sys.executable, '-c', WORKER, root, fd, str(self.lock_path), mode],
                                cwd=BASE, **kwargs)
            self.children.append(worker)
            return worker
        return launch

    def guard(self, action=None):
        outer = self
        class FixtureGuard:
            def __init__(self, log_path, root_pid=None, logger=None):
                self.root_pid = root_pid
                self.log_path = log_path
                self.done = False
                self.deadline = time.monotonic() + 8
            def check(self, stage='working', starting=False):
                self.log_path.parent.mkdir(parents=True, exist_ok=True)
                with self.log_path.open('a') as out:
                    out.write(json.dumps(dict(stage=stage,starting=starting,synthetic=True)) + '\n')
                return {'status':'OK'}
            def checkpoint(self, stage='working', force=False):
                if action and not self.done:
                    if (outer.root/'ready.json').exists():
                        self.done = True
                        with outer.assertRaises(runtime_safety.DesignAlreadyRunning):
                            runtime_safety.acquire_design_lock(path=outer.lock_path)
                        if action == 'pressure':
                            raise runtime_safety.MemoryLimitExceeded('Synthetic memory pressure: analysis incomplete.')
                        if action == 'disk':
                            partial=outer.root/'result/blast_attempts/human-fixture/0001_candidates_blast.tsv'
                            partial.parent.mkdir(parents=True,exist_ok=True);partial.write_text('uncommitted output')
                            raise runtime_safety.DiskSpaceLimitExceeded('Synthetic disk reserve stop: analysis incomplete.')
                        if action == 'cancel':
                            outer.state[outer.job_id]['cancel'].set()
                    elif time.monotonic() > self.deadline:
                        raise RuntimeError('Tiny child startup timed out.')
                return self.check(stage)
        return FixtureGuard

    def assert_incomplete(self):
        self.assertEqual(self.state[self.job_id]['status'], 'incomplete')
        self.assertEqual(json.loads((self.root/'job_status.json').read_text())['status'], 'incomplete')
        persisted = json.loads((self.root/'analysis_status.json').read_text())
        self.assertEqual(persisted['status'], 'incomplete')
        self.assertFalse(persisted['analysis_complete'])
        self.assertFalse((self.jobs/(self.job_id+'.zip')).exists())
        process_record=json.loads((self.root/'supervisor_process.json').read_text())
        self.assertEqual(process_record['status'],'INCOMPLETE')
        self.assertFalse(process_record['analysis_complete'])
        self.assertIsNotNone(process_record['returncode'])
        self.assertTrue((self.root/'result'/'DRAFT_IDT_oPools.csv').exists())
        with runtime_safety.acquire_design_lock(path=self.lock_path):
            pass
        for route in ('/download?id=fixture_job', '/api/result?id=fixture_job'):
            handler = studio.Handler.__new__(studio.Handler)
            handler.path, handler.headers = route, {'Host': '127.0.0.1'}
            calls = []
            handler.send = lambda data, status=200, **kwargs: calls.append((data, status))
            handler.do_GET()
            self.assertEqual(calls[-1][1], 400)
            self.assertIn('incomplete', calls[-1][0]['error'].lower())

    def assert_stopped_group(self):
        ready = json.loads((self.root/'ready.json').read_text())
        self.assertEqual(ready['worker_group'], ready['worker_pid'])
        self.assertEqual(ready['child_group'], ready['worker_group'])
        self.assertIsNotNone(self.children[0].poll())
        deadline = time.monotonic() + 3
        while running(ready['child_pid']) and time.monotonic() < deadline:
            time.sleep(.03)
        self.assertFalse(running(ready['child_pid']), 'Owned BLAST-like child survived stop.')

    def test_pressure_stops_worker_and_term_resistant_child(self):
        self.prepare()
        with patch('runtime_safety.MemoryGuard', self.guard('pressure')), patch('studio.subprocess.Popen', self.factory('waiting')):
            studio.run_job(self.job_id, {'fixture': True})
        self.assert_stopped_group()
        self.assert_incomplete()
        self.assertIn('Synthetic memory pressure', self.state[self.job_id]['error'])

    def test_disk_stop_terminates_owned_group_and_removes_only_uncommitted_payload(self):
        self.prepare()
        with patch('runtime_safety.MemoryGuard', self.guard('disk')), patch('studio.subprocess.Popen', self.factory('waiting')):
            studio.run_job(self.job_id, {'fixture': True})
        self.assert_stopped_group()
        self.assert_incomplete()
        self.assertIn('disk reserve', self.state[self.job_id]['error'])
        self.assertFalse(list((self.root/'result/blast_attempts').rglob('*_candidates_blast.tsv')))
        self.assertTrue((self.root/'storage_events.jsonl').exists())
        self.assertTrue((self.root/'supervisor_process.json').exists())

    def test_cancel_stops_worker_and_term_resistant_child(self):
        self.prepare()
        with patch('runtime_safety.MemoryGuard', self.guard('cancel')), patch('studio.subprocess.Popen', self.factory('waiting')):
            studio.run_job(self.job_id, {'fixture': True})
        self.assert_stopped_group()
        self.assert_incomplete()
        self.assertIn('Cancelled', self.state[self.job_id]['error'])

    def test_worker_failure_is_persistently_incomplete(self):
        self.prepare()
        with patch('runtime_safety.MemoryGuard', self.guard()), patch('studio.subprocess.Popen', self.factory('crash')):
            studio.run_job(self.job_id, {'fixture': True})
        self.assert_incomplete()
        self.assertIn('Tiny fixture worker failed', self.state[self.job_id]['error'])

    def test_zip_close_failure_never_publishes_complete_journal(self):
        self.prepare()
        observed = []
        outer = self
        class FailZipClose(REAL_ZIPFILE):
            failed = False
            def close(self):
                if self.fp is not None and not self.failed:
                    self.failed = True
                    observed.append(json.loads((outer.root/'analysis_status.json').read_text())['status'])
                    super().close()
                    raise OSError('Synthetic ZIP close failure')
                return super().close()
        with patch('runtime_safety.MemoryGuard', self.guard()), patch('studio.subprocess.Popen', self.factory('success')), patch('studio.standalone', return_value='<html>DRAFT fixture</html>'), patch('studio.zipfile.ZipFile', FailZipClose):
            studio.run_job(self.job_id, {'fixture': True})
        self.assert_incomplete()
        self.assertEqual(observed, ['running'], 'Local complete status was published before ZIP closed.')

    def test_zip_rename_failure_never_publishes_complete_journal(self):
        self.prepare()
        observed = []
        outer = self
        def fail_replace(path, target):
            if path.name.endswith('.incomplete.zip'):
                observed.append(json.loads((outer.root/'analysis_status.json').read_text())['status'])
                raise OSError('Synthetic ZIP rename failure')
            return REAL_REPLACE(path, target)
        with patch('runtime_safety.MemoryGuard', self.guard()), patch('studio.subprocess.Popen', self.factory('success')), patch('studio.standalone', return_value='<html>DRAFT fixture</html>'), patch.object(Path, 'replace', fail_replace):
            studio.run_job(self.job_id, {'fixture': True})
        self.assert_incomplete()
        self.assertTrue(observed)
        self.assertEqual(set(observed), {'running'}, 'Local complete status was published before ZIP renamed.')

    def test_guard_signal_journal_precedes_kill_and_unrelated_process_survives(self):
        self.prepare()
        other=REAL_POPEN([sys.executable,'-c','import time; time.sleep(20)'],start_new_session=True)
        real_killpg=os.killpg
        observed=[]
        def signal_owned(pgid,sig):
            if sig:
                events=[json.loads(line) for line in (self.root/'termination_events.jsonl').read_text().splitlines()]
                self.assertEqual(events[-1]['signal'],signal.Signals(sig).name)
                self.assertIn('Synthetic memory pressure',events[-1]['reason'])
                self.assertNotEqual(pgid,other.pid)
                observed.append(sig)
            return real_killpg(pgid,sig)
        try:
            with patch('runtime_safety.MemoryGuard',self.guard('pressure')), patch('studio.subprocess.Popen',self.factory('waiting')), patch('studio.os.killpg',signal_owned):
                studio.run_job(self.job_id,{'fixture':True})
            self.assertTrue(observed)
            self.assertIsNone(other.poll())
            self.assert_stopped_group()
            self.assert_incomplete()
        finally:
            other.terminate();other.wait(timeout=3)

    def test_completed_worker_group_is_retired_before_later_cleanup(self):
        self.prepare()
        with patch('runtime_safety.MemoryGuard',self.guard()), patch('studio.subprocess.Popen',self.factory('success')), patch('studio.standalone',return_value='<html>DRAFT</html>'):
            studio.run_job(self.job_id,{'fixture':True})
        child=self.children[0]
        self.assertTrue(child._raserfish_group_retired)
        with patch('studio.os.killpg') as signal_group:
            self.assertEqual(studio.stop_worker(child,self.root),[])
            signal_group.assert_not_called()

    def test_termination_journal_failure_still_stops_owned_group(self):
        self.prepare()
        other=REAL_POPEN([sys.executable,'-c','import time; time.sleep(20)'],start_new_session=True)
        try:
            with patch('runtime_safety.MemoryGuard',self.guard('pressure')), patch('studio.subprocess.Popen',self.factory('waiting')), patch('studio.append_event',side_effect=OSError('Fixture journal disk failure')), patch('studio.os.write',return_value=0):
                studio.run_job(self.job_id,{'fixture':True})
            self.assert_stopped_group()
            self.assert_incomplete()
            self.assertIsNone(other.poll())
        finally:
            other.terminate();other.wait(timeout=3)

    def test_failure_after_archive_rename_quarantines_unpublished_zip(self):
        self.prepare()
        original_atomic=studio.atomic_json
        failed=[False]
        def write(path,payload):
            if Path(path).name=='job_status.json' and payload.get('status')=='complete' and not failed[0]:
                failed[0]=True
                self.assertTrue((self.jobs/(self.job_id+'.zip')).exists())
                raise OSError('Fixture complete journal failure after archive rename')
            return original_atomic(path,payload)
        with patch('runtime_safety.MemoryGuard',self.guard()), patch('studio.subprocess.Popen',self.factory('success')), patch('studio.standalone',return_value='<html>DRAFT</html>'), patch('studio.atomic_json',write):
            studio.run_job(self.job_id,{'fixture':True})
        self.assertTrue(failed[0]);self.assert_incomplete()
        quarantined=self.jobs/(self.job_id+'.INCOMPLETE_EXPORT')
        self.assertTrue(quarantined.exists())
        marker=json.loads((self.jobs/(self.job_id+'.INCOMPLETE_EXPORT.json')).read_text())
        self.assertEqual(marker['status'],'INCOMPLETE')
        self.assertIn('invalid',marker['note'])
        self.assertFalse((self.jobs/(self.job_id+'.zip')).exists())

    def test_completed_archive_and_journals_agree(self):
        self.prepare()
        with patch('runtime_safety.MemoryGuard', self.guard()), patch('studio.subprocess.Popen', self.factory('success')), patch('studio.standalone', return_value='<html>DRAFT fixture</html>'):
            studio.run_job(self.job_id, {'fixture': True})
        self.assertEqual(self.state[self.job_id]['status'], 'complete')
        self.assertTrue(json.loads((self.root/'analysis_status.json').read_text())['analysis_complete'])
        with REAL_ZIPFILE(self.jobs/(self.job_id+'.zip')) as archive:
            self.assertEqual(json.loads(archive.read(self.job_id+'/job_status.json'))['status'], 'complete')
            self.assertTrue(json.loads(archive.read(self.job_id+'/analysis_status.json'))['analysis_complete'])
            self.assertIn('DRAFT', archive.read(self.job_id+'/result/DRAFT_IDT_oPools.csv').decode())
        with runtime_safety.acquire_design_lock(path=self.lock_path):
            pass


if __name__ == '__main__':
    unittest.main()
