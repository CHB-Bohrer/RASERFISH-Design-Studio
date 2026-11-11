"""Small real subprocesses; no reference databases or design searches are used."""
import io
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import threading
import unittest
from unittest.mock import patch

import process_diagnostics as diagnostics
import runtime_safety
import studio


class ProcessDiagnosticsTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.log = self.root / 'human_blast.log'
        self.environment = patch.dict(os.environ, {'RASERFISH_DESIGN_WORKER':'0'})
        self.environment.start()

    def tearDown(self):
        self.environment.stop()
        self.temp.cleanup()

    def record(self):
        return json.loads(Path(str(self.log)+'.process.json').read_text())

    def test_success_records_exact_command_and_distinct_disk_streams(self):
        command=[sys.executable,'-c','import sys; print("stdout only"); print("stderr only",file=sys.stderr)']
        diagnostics.run_process(command,self.log)
        record=self.record()
        self.assertEqual(record['argv'],command)
        self.assertEqual(record['cwd'],os.getcwd())
        self.assertEqual(record['status'],'COMPLETE')
        self.assertEqual(record['returncode'],0)
        self.assertTrue(record['pid'])
        self.assertTrue(record['pgid'])
        self.assertEqual(Path(record['stdout_path']).read_text(),'stdout only\n')
        self.assertEqual(Path(record['stderr_path']).read_text(),'stderr only\n')
        self.assertIn('return code 0',self.log.read_text())

    def test_explicit_environment_is_child_only_and_journal_excludes_inherited_values(self):
        command=[sys.executable,'-c','import json,os; print(json.dumps({k:os.environ[k] for k in ["BATCH_SIZE","FIXTURE_INHERITED"]}))']
        overrides={'BATCH_SIZE':'91'}
        with patch.dict(os.environ,{'BATCH_SIZE':'parent-value','FIXTURE_INHERITED':'inherited-value'}):
            before=dict(os.environ)
            diagnostics.run_process(command,self.log,env_overrides=overrides)
            self.assertEqual(dict(os.environ),before)
        self.assertEqual(overrides,{'BATCH_SIZE':'91'})
        record=self.record()
        self.assertEqual(record['argv'],command)
        self.assertEqual(record['executable'],sys.executable)
        self.assertEqual(record['environment_overrides'],overrides)
        self.assertNotIn('FIXTURE_INHERITED',record['environment_overrides'])
        self.assertEqual(json.loads(Path(record['stdout_path']).read_text()),
                         {'BATCH_SIZE':'91','FIXTURE_INHERITED':'inherited-value'})
        self.assertIn('Environment overrides: {"BATCH_SIZE": "91"}',self.log.read_text())

    def test_failed_child_retains_explicit_environment_and_original_exit(self):
        command=[sys.executable,'-c','import os,sys; print(os.environ["BATCH_SIZE"],file=sys.stderr); sys.exit(7)']
        with self.assertRaisesRegex(RuntimeError,'return code 7'):
            diagnostics.run_process(command,self.log,env_overrides={'BATCH_SIZE':'91'})
        record=self.record()
        self.assertEqual(record['environment_overrides'],{'BATCH_SIZE':'91'})
        self.assertEqual(record['argv'],command)
        self.assertEqual(record['status'],'INCOMPLETE')
        self.assertEqual(record['returncode'],7)
        self.assertEqual(record['stderr_tail'],'91\n')
        self.assertFalse(record['analysis_complete'])

    def test_guard_stop_retains_explicit_environment_before_owned_signal(self):
        class Guard:
            def checkpoint(self,*_):
                raise runtime_safety.MemoryLimitExceeded('Fixture override guard',{'fixture_limit':123})
        original_stop=diagnostics._stop_owned_child
        seen=[]
        def checked_stop(proc,own_group,record,diagnostic,reason):
            persisted=json.loads(Path(diagnostic).read_text())
            self.assertEqual(persisted['environment_overrides'],{'BATCH_SIZE':'91'})
            self.assertEqual(persisted['status'],'STOP_REQUESTED')
            seen.append(proc.pid)
            return original_stop(proc,own_group,record,diagnostic,reason)
        with patch.object(diagnostics,'_stop_owned_child',side_effect=checked_stop):
            with self.assertRaisesRegex(runtime_safety.MemoryLimitExceeded,'Fixture override guard'):
                diagnostics.run_process([sys.executable,'-c','import time; time.sleep(20)'],self.log,
                                        resource_guard=Guard(),env_overrides={'BATCH_SIZE':'91'})
        record=self.record()
        self.assertEqual(seen,[record['pid']])
        self.assertEqual(record['environment_overrides'],{'BATCH_SIZE':'91'})
        self.assertEqual(record['status'],'INCOMPLETE')
        self.assertFalse(record['analysis_complete'])
        self.assertTrue(record['termination_actions'])
        self.assertTrue(all(action['target_pid']==record['pid'] for action in record['termination_actions']))
        self.assertEqual(record['memory_sample'],{'fixture_limit':123})

    def test_exit7_with_empty_streams_is_useful_and_incomplete(self):
        with self.assertRaisesRegex(RuntimeError,'return code 7') as context:
            diagnostics.run_process([sys.executable,'-c','raise SystemExit(7)'],self.log)
        self.assertIn('No stdout/stderr was emitted',str(context.exception))
        self.assertEqual(self.record()['status'],'INCOMPLETE')
        self.assertFalse(self.record()['analysis_complete'])
        self.assertEqual(self.record()['returncode'],7)
        self.assertIn('return code 7',self.log.read_text())
        self.assertIn('.process.json',str(context.exception))

    def test_terminating_signal_is_observed_not_guessed(self):
        with self.assertRaisesRegex(RuntimeError,'SIGTERM'):
            diagnostics.run_process([sys.executable,'-c','import os,signal; os.kill(os.getpid(),signal.SIGTERM)'],self.log)
        self.assertEqual(self.record()['returncode'],-signal.SIGTERM)
        self.assertEqual(self.record()['terminating_signal'],'SIGTERM')
        self.assertEqual(self.record()['termination_actions'],[])

    def test_missing_executable_records_spawn_failure(self):
        missing=str(self.root/'absent_blastn')
        with self.assertRaisesRegex(RuntimeError,'could not start'):
            diagnostics.run_process([missing,'-query','fixture.fa'],self.log)
        self.assertIsNone(self.record()['returncode'])
        self.assertIn('FileNotFoundError',self.record()['spawn_error'])
        self.assertEqual(self.record()['executable'],missing)
        self.assertEqual(self.record()['argv'][1:],['-query','fixture.fa'])

    def test_large_stdout_stays_on_disk_with_bounded_tail(self):
        diagnostics.run_process([sys.executable,'-c','import sys; sys.stdout.write("A"*2000000)'],self.log)
        record=self.record()
        self.assertEqual(record['stdout_bytes'],2000000)
        self.assertEqual(len(record['stdout_tail']),diagnostics.TAIL_BYTES)
        self.assertLess(Path(str(self.log)+'.process.json').stat().st_size,12000)

    def test_memory_stop_is_recorded_before_signal_and_other_process_survives(self):
        original_killpg=os.killpg
        unrelated=subprocess.Popen([sys.executable,'-c','import time; time.sleep(20)'],start_new_session=True)
        seen=[]
        class Guard:
            def checkpoint(self,*_):
                raise runtime_safety.MemoryLimitExceeded('Fixture guard: limit exceeded',{'fixture_limit':123})
        def checked_signal(pid,sig):
            if sig:
                record=self.record()
                self.assertEqual(record['status'],'STOP_REQUESTED')
                self.assertEqual(record['memory_sample'],{'fixture_limit':123})
                self.assertEqual(record['termination_actions'][-1]['signal'],signal.Signals(sig).name)
                self.assertNotEqual(pid,unrelated.pid)
                seen.append(sig)
            return original_killpg(pid,sig)
        try:
            with patch('process_diagnostics.os.killpg',checked_signal):
                with self.assertRaises(runtime_safety.MemoryLimitExceeded):
                    diagnostics.run_process([sys.executable,'-c','import time; time.sleep(20)'],self.log,resource_guard=Guard())
            self.assertTrue(seen)
            self.assertEqual(self.record()['status'],'INCOMPLETE')
            self.assertIsNone(unrelated.poll())
            self.assertIn('Fixture guard',self.log.read_text())
        finally:
            unrelated.terminate();unrelated.wait(timeout=3)

    def test_cancellation_retains_distinct_exception_and_stop_reason(self):
        class Cancelled(RuntimeError):pass
        event=threading.Event();event.set()
        with self.assertRaises(Cancelled):
            diagnostics.run_process([sys.executable,'-c','import time; time.sleep(20)'],self.log,event,cancellation_class=Cancelled)
        self.assertEqual(self.record()['stop_exception'],'Cancelled')
        self.assertIn('Cancelled',self.record()['stop_reason'])
        self.assertEqual(self.record()['status'],'INCOMPLETE')

    def test_journal_failure_after_spawn_terminates_owned_child(self):
        original_popen=subprocess.Popen
        original_write=diagnostics.atomic_json
        children=[];failed=[False]
        def spawn(*args,**kwargs):
            child=original_popen(*args,**kwargs);children.append(child);return child
        def write(path,record):
            if record.get('status')=='RUNNING' and not failed[0]:
                failed[0]=True
                raise OSError('Fixture diagnostic disk write failed after spawn')
            return original_write(path,record)
        with patch('process_diagnostics.subprocess.Popen',spawn),patch('process_diagnostics.atomic_json',write):
            with self.assertRaisesRegex(OSError,'disk write failed'):
                diagnostics.run_process([sys.executable,'-c','import time; time.sleep(20)'],self.log)
        self.assertTrue(failed[0])
        self.assertEqual(len(children),1)
        self.assertIsNotNone(children[0].poll())
        self.assertEqual(self.record()['status'],'INCOMPLETE')

    def test_summary_failure_after_spawn_terminates_owned_child(self):
        original_popen=subprocess.Popen
        original_write=diagnostics.write_summary
        children=[];failed=[False]
        def spawn(*args,**kwargs):
            child=original_popen(*args,**kwargs);children.append(child);return child
        def write(path,record,diagnostic):
            if record.get('status')=='RUNNING' and not failed[0]:
                failed[0]=True
                raise OSError('Fixture summary write failed after spawn')
            return original_write(path,record,diagnostic)
        with patch('process_diagnostics.subprocess.Popen',spawn),patch('process_diagnostics.write_summary',write):
            with self.assertRaisesRegex(OSError,'summary write failed'):
                diagnostics.run_process([sys.executable,'-c','import time; time.sleep(20)'],self.log)
        self.assertEqual(len(children),1)
        self.assertIsNotNone(children[0].poll())
        self.assertEqual(self.record()['status'],'INCOMPLETE')

    def test_supervisor_finalizes_inflight_without_inventing_child_returncode(self):
        record=diagnostics.command_record(['blastn','-query','fixture.fa'])
        stdout=self.root/'stdout';stderr=self.root/'stderr'
        stdout.write_text('');stderr.write_text('')
        record.update(status='RUNNING',stdout_path=str(stdout),stderr_path=str(stderr),summary_path=str(self.log))
        diagnostic=Path(str(self.log)+'.process.json')
        diagnostics.atomic_json(diagnostic,record)
        diagnostics.finalize_interrupted_commands(self.root,'Parent memory guard stop',[{'signal':'SIGTERM'}])
        record=self.record()
        self.assertEqual(record['status'],'INCOMPLETE')
        self.assertIsNone(record['returncode'])
        self.assertIsNone(record['terminating_signal'])
        self.assertFalse(record['child_returncode_observed'])
        self.assertEqual(record['supervisor_actions'],[{'signal':'SIGTERM'}])


class WorkerPipeTests(unittest.TestCase):
    def test_broken_supervisor_pipe_does_not_replace_original_failure(self):
        class Lease:
            def release(self):pass
        class Guard:
            def __init__(self,*_,**__):pass
            def checkpoint(self,*_,**__):pass
        class Broken(io.StringIO):
            def write(self,*_):raise BrokenPipeError('fixture closed supervisor pipe')
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            with patch('runtime_safety.inherited_design_lock',return_value=Lease()), patch('runtime_safety.MemoryGuard',Guard), patch('studio.sys.stdin',io.StringIO('{}')), patch('studio.sys.stdout',Broken()), patch('studio.sys.stderr',io.StringIO()), patch('design_engine.design',side_effect=RuntimeError('Original BLAST return code 7')):
                status=studio.worker_main(root,1)
            self.assertEqual(status,1)
            persisted=json.loads((root/'worker_status.json').read_text())
            self.assertEqual(persisted['status'],'incomplete')
            self.assertEqual(persisted['error'],'Original BLAST return code 7')
            self.assertIn('SUPERVISOR_PIPE_LOST',(root/'worker_events.jsonl').read_text())

    def test_failed_progress_report_leaves_worker_incomplete(self):
        class Lease:
            def release(self):pass
        class Guard:
            def __init__(self,*_,**__):pass
            def checkpoint(self,*_,**__):pass
        class Broken(io.StringIO):
            def write(self,*_):raise BrokenPipeError('fixture closed supervisor pipe')
        def design(_config,_output,report,**kwargs):
            report('Working')
            return {'counts':{},'status':'would be complete'}
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            with patch('runtime_safety.inherited_design_lock',return_value=Lease()), patch('runtime_safety.MemoryGuard',Guard), patch('studio.sys.stdin',io.StringIO('{}')), patch('studio.sys.stdout',Broken()), patch('studio.sys.stderr',io.StringIO()), patch('design_engine.design',side_effect=design):
                self.assertEqual(studio.worker_main(root,1),1)
            self.assertEqual(json.loads((root/'worker_status.json').read_text())['status'],'incomplete')


class ParentLossTests(unittest.TestCase):
    def test_lost_supervisor_stops_only_owned_worker_family_and_marks_incomplete(self):
        self.check_parent_loss(False)

    def test_lost_supervisor_still_stops_family_when_journal_disk_fails(self):
        self.check_parent_loss(True)

    def check_parent_loss(self,fail_journal):
        import time
        worker_code=r'''
import io,json,os,subprocess,sys,time
from pathlib import Path
import runtime_safety,design_engine,studio
root=Path(sys.argv[1])
if sys.argv[2]=='1':
    def fail(*a,**kw):raise OSError('Fixture journal disk failure')
    studio.atomic_json=fail;studio.append_event=fail
class Lease:
    def release(self):pass
class Guard:
    def __init__(self,*a,**k):pass
    def checkpoint(self,*a,**k):pass
runtime_safety.inherited_design_lock=lambda fd:Lease()
runtime_safety.MemoryGuard=Guard
sys.stdin=io.StringIO('{}')
def design(config,out,report,**kw):
    code='import signal,time; signal.signal(signal.SIGTERM,signal.SIG_IGN); time.sleep(20)'
    child=subprocess.Popen([sys.executable,'-c',code])
    (root/'ready.json').write_text(json.dumps(dict(worker_pid=os.getpid(),child_pid=child.pid)))
    time.sleep(20)
    return dict(counts={},status='fixture')
design_engine.design=design
raise SystemExit(studio.worker_main(root,1))
'''
        parent_code=r'''
import os,subprocess,sys,time
from pathlib import Path
root=Path(sys.argv[1])
with (root/'fixture_stderr.log').open('w') as err:
    child=subprocess.Popen([sys.executable,'-c',sys.argv[2],str(root),sys.argv[3]],start_new_session=True,stdout=subprocess.DEVNULL,stderr=err)
    while True:time.sleep(1)
'''
        def alive(pid):
            try:os.kill(pid,0)
            except ProcessLookupError:return False
            stat=Path(f'/proc/{pid}/stat')
            if stat.exists():
                try:
                    if stat.read_text().split(')',1)[1].strip().split()[0]=='Z':return False
                except (OSError,IndexError):pass
            return True
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary)
            parent=subprocess.Popen([sys.executable,'-c',parent_code,str(root),worker_code,'1' if fail_journal else '0'],cwd=Path(__file__).resolve().parents[1],start_new_session=True)
            other=subprocess.Popen([sys.executable,'-c','import time; time.sleep(20)'],start_new_session=True)
            ready=None
            try:
                deadline=time.monotonic()+8
                while not (root/'ready.json').exists() and time.monotonic()<deadline:time.sleep(.03)
                self.assertTrue((root/'ready.json').exists())
                ready=json.loads((root/'ready.json').read_text())
                parent.terminate();parent.wait(timeout=3)
                while (alive(ready['worker_pid']) or alive(ready['child_pid'])) and time.monotonic()<deadline:time.sleep(.05)
                self.assertFalse(alive(ready['worker_pid']))
                self.assertFalse(alive(ready['child_pid']))
                self.assertIsNone(other.poll())
                if fail_journal:
                    self.assertIn('Fixture journal disk failure',(root/'fixture_stderr.log').read_text())
                    self.assertIn('SUPERVISOR_LOST_STOP',(root/'fixture_stderr.log').read_text())
                    self.assertFalse((root/'worker_status.json').exists())
                else:
                    self.assertEqual(json.loads((root/'worker_status.json').read_text())['status'],'incomplete')
                    self.assertFalse(json.loads((root/'analysis_status.json').read_text())['analysis_complete'])
                    event=json.loads((root/'termination_events.jsonl').read_text().splitlines()[-1])
                    self.assertEqual(event['action'],'SUPERVISOR_LOST_STOP')
                    self.assertEqual(event['signal'],'SIGKILL')
                    self.assertEqual(event['target_pgid'],ready['worker_pid'])
            finally:
                if parent.poll() is None:parent.terminate();parent.wait(timeout=3)
                if other.poll() is None:other.terminate();other.wait(timeout=3)
                if ready and alive(ready['worker_pid']):
                    os.killpg(ready['worker_pid'],signal.SIGKILL)

if __name__=='__main__':unittest.main()
