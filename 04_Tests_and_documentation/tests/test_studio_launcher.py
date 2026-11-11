"""Owned tiny HTTP workers exercise server recovery; no design/BLAST is run."""
import json
import os
from pathlib import Path
import signal
import socket
import subprocess
import sys
import tempfile
import time
import unittest
from unittest.mock import patch

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE))
import studio_launcher as launcher

FAKE_SERVER = r'''
import argparse,json,os,signal,sys,threading,time
from pathlib import Path
from http.server import BaseHTTPRequestHandler,ThreadingHTTPServer
root=Path(__file__).resolve().parent
parser=argparse.ArgumentParser();parser.add_argument('--qc-root');parser.add_argument('--port',type=int);parser.add_argument('--no-open',action='store_true');args=parser.parse_args()
settings=json.loads((root/'fixture.json').read_text())
if settings.get('startup_hang'):time.sleep(30)
with (root/'starts.jsonl').open('a') as handle:handle.write(json.dumps(dict(pid=os.getpid(),port=args.port,pgid=os.getpgrp(),session_id=os.getsid(0)))+'\n')
count=len((root/'starts.jsonl').read_text().splitlines())
if settings.get('output_bytes'):sys.stdout.write('x'*settings['output_bytes']);sys.stdout.flush()
class Handler(BaseHTTPRequestHandler):
 def do_GET(self):
  body=json.dumps(dict(version='1.6.1',pid=os.getpid(),instance_id='fixture',app_path=str(root),port=self.server.server_port)).encode()
  self.send_response(200);self.send_header('Content-Length',str(len(body)));self.end_headers();self.wfile.write(body)
 def log_message(self,*args):pass
server=ThreadingHTTPServer(('127.0.0.1',args.port),Handler)
exit_codes=settings.get('exit_codes',[])
if count<=len(exit_codes):
 code=exit_codes[count-1]
 def exit_later():
  time.sleep(settings.get('exit_delay',.20))
  if code<0:os.kill(os.getpid(),-code)
  else:os._exit(code)
 threading.Thread(target=exit_later,daemon=True).start()
server.serve_forever()
'''

LAUNCH_HELPER = '''
import sys
from pathlib import Path
sys.path.insert(0,sys.argv[1])
from studio_launcher import LauncherConfig,StudioLauncher
root=Path(sys.argv[2])
c=LauncherConfig(root,root,port=int(sys.argv[3]),no_open=True,backoff=(.05,.10),poll_seconds=.01,startup_timeout=3,stop_grace=.5)
settings=__import__('json').loads((root/'fixture.json').read_text())
if 'launcher_backoff' in settings:c.backoff=tuple(settings['launcher_backoff'])
if 'launcher_timeout' in settings:c.startup_timeout=settings['launcher_timeout']
raise SystemExit(StudioLauncher(c).run())
'''


class LauncherTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name).resolve()
        (self.root / 'studio.py').write_text(FAKE_SERVER)
        (self.root / 'fixture.json').write_text('{}')
        self.processes = []
        self.outputs = []

    def tearDown(self):
        for proc in reversed(self.processes):
            if proc.poll() is None:
                proc.terminate()
                try: proc.wait(timeout=3)
                except subprocess.TimeoutExpired:
                    proc.kill(); proc.wait(timeout=3)
        for output in self.outputs: output.close()
        self.tmp.cleanup()

    def settings(self, **kwargs):
        (self.root / 'fixture.json').write_text(json.dumps(kwargs))

    def start(self, port=0):
        output = (self.root / ('launcher-output-' + str(len(self.outputs)) + '.txt')).open('wb')
        self.outputs.append(output)
        proc = subprocess.Popen([sys.executable, '-c', LAUNCH_HELPER, str(BASE), str(self.root), str(port)],
                                stdout=output, stderr=output, start_new_session=True)
        self.processes.append(proc)
        return proc

    def starts(self):
        try: return [json.loads(line) for line in (self.root / 'starts.jsonl').read_text().splitlines()]
        except FileNotFoundError: return []

    def state(self):
        try: return json.loads((self.root / 'logs/studio_launcher.json').read_text())
        except (FileNotFoundError, ValueError): return {}

    def wait_for(self, predicate, timeout=8):
        deadline=time.monotonic()+timeout
        while time.monotonic()<deadline:
            value=predicate()
            if value:return value
            time.sleep(.02)
        self.fail('Timed out. State: '+str(self.state())+'; starts: '+str(self.starts()))

    def ready(self):
        state=self.state()
        return state if state.get('event')=='SERVER_READY' else None

    def events(self):
        return [json.loads(line) for line in (self.root/'logs/studio_launcher.jsonl').read_text().splitlines()]

    def test_crash_restarts_once_at_same_port_then_ready(self):
        self.settings(exit_codes=[7])
        proc=self.start();self.wait_for(lambda:len(self.starts())==2 and self.ready())
        starts=self.starts();self.assertEqual(len(set(row['port'] for row in starts)),1)
        self.assertEqual(len(set(row['pid'] for row in starts)),2)
        self.assertTrue(all(row['pid']==row['pgid']==row['session_id'] for row in starts))
        self.assertEqual(self.state()['restarts'],1)
        self.assertTrue(any(event.get('returncode')==7 for event in self.events()))
        proc.terminate();self.assertEqual(proc.wait(timeout=5),143)

    def test_repeated_crash_stops_after_two_restarts(self):
        self.settings(exit_codes=[9,9,9,9])
        proc=self.start();self.assertEqual(proc.wait(timeout=8),9)
        self.assertEqual(len(self.starts()),3)
        self.assertEqual(self.state()['event'],'RESTART_LIMIT_REACHED')
        self.assertEqual(self.state()['status'],'failed')
        self.assertIn('incomplete',self.state()['reason'])

    def test_unexpected_sigkill_restarts(self):
        self.settings(exit_codes=[-signal.SIGKILL])
        proc=self.start();self.wait_for(lambda:len(self.starts())==2 and self.ready())
        exited=[event for event in self.events() if event['event']=='SERVER_EXITED']
        self.assertEqual(exited[0]['terminating_signal'],'SIGKILL')
        proc.terminate();proc.wait(timeout=5)

    def test_child_explicit_signals_do_not_restart(self):
        for sig in (signal.SIGTERM,signal.SIGINT,signal.SIGHUP):
            with self.subTest(signal=sig):
                # Use a fresh app root to keep fixture exit sequence independent.
                self.settings(exit_codes=[-sig],exit_delay=.15)
                (self.root/'starts.jsonl').unlink(missing_ok=True)
                proc=self.start();self.assertEqual(proc.wait(timeout=5),128+sig)
                self.assertEqual(len(self.starts()),1)
                self.assertEqual(self.state()['event'],'LAUNCHER_STOPPED')

    def test_child_handled_signal_exit_codes_do_not_restart(self):
        for code in (0,129,130,143):
            with self.subTest(code=code):
                self.settings(exit_codes=[code],exit_delay=.15)
                (self.root/'starts.jsonl').unlink(missing_ok=True)
                proc=self.start();self.assertEqual(proc.wait(timeout=5),code)
                self.assertEqual(len(self.starts()),1)

    def test_launcher_ctrl_c_stops_only_its_child_and_no_respawn(self):
        unrelated=subprocess.Popen([sys.executable,'-c','import time;time.sleep(30)'])
        self.processes.append(unrelated)
        proc=self.start();state=self.wait_for(self.ready);child_pid=state['server_pid']
        proc.send_signal(signal.SIGINT);self.assertEqual(proc.wait(timeout=5),130)
        self.assertEqual(len(self.starts()),1);self.assertIsNone(unrelated.poll())
        with self.assertRaises(ProcessLookupError):os.kill(child_pid,0)
        stop=[event for event in self.events() if event['event']=='OWNED_SERVER_STOP_REQUESTED']
        self.assertEqual([event['server_pid'] for event in stop],[child_pid])
        self.assertEqual(self.state()['signal'],'SIGINT')

    def test_launcher_hangup_stops_without_restart(self):
        proc=self.start();self.wait_for(self.ready)
        proc.send_signal(signal.SIGHUP);self.assertEqual(proc.wait(timeout=5),129)
        self.assertEqual(len(self.starts()),1)

    def test_duplicate_launcher_reuses_server_without_spawning(self):
        proc=self.start();ready=self.wait_for(self.ready)
        duplicate=self.start();self.assertEqual(duplicate.wait(timeout=5),0)
        self.assertEqual(len(self.starts()),1);self.assertIsNone(proc.poll())
        self.assertEqual(self.state()['launcher_pid'],proc.pid)
        self.assertIsNotNone(launcher.server_identity(ready['port'],self.root))
        proc.terminate();proc.wait(timeout=5)

    def test_duplicate_during_restart_backoff_never_spawns(self):
        first=launcher.StudioLauncher(launcher.LauncherConfig(self.root,self.root,0,True))
        self.assertTrue(first.acquire())
        try:
            duplicate=self.start();self.assertEqual(duplicate.wait(timeout=5),0)
            self.assertEqual(self.starts(),[])
        finally:os.close(first.lock_fd)

    def test_existing_direct_server_reused_without_taking_ownership(self):
        port=launcher.choose_port(0)
        proc=subprocess.Popen([sys.executable,str(self.root/'studio.py'),'--port',str(port)],
                              stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,start_new_session=True)
        self.processes.append(proc)
        self.wait_for(lambda:launcher.server_identity(port,self.root))
        duplicate=self.start(port);self.assertEqual(duplicate.wait(timeout=5),0)
        self.assertEqual(len(self.starts()),1);self.assertIsNone(proc.poll())
        self.assertFalse((self.root/'logs/studio_launcher.json').exists())

    def test_occupied_foreign_port_uses_free_port_without_touching_owner(self):
        with socket.socket() as occupied:
            occupied.bind(('127.0.0.1',0));occupied.listen()
            port=occupied.getsockname()[1]
            proc=self.start(port);ready=self.wait_for(self.ready)
            self.assertNotEqual(port,ready['port'])
            self.assertEqual(occupied.getsockname()[1],port)
            proc.terminate();proc.wait(timeout=5)

    def test_large_console_output_goes_directly_to_file(self):
        self.settings(output_bytes=2*1024*1024)
        proc=self.start();state=self.wait_for(self.ready)
        self.assertGreaterEqual(Path(state['console_log']).stat().st_size,2*1024*1024)
        self.assertEqual(state['argv'][2],str(self.root/'studio.py'))
        self.assertEqual(state['argv'][-1],'--no-open')
        proc.terminate();proc.wait(timeout=5)

    def test_launcher_releases_lock_after_exit(self):
        self.settings(exit_codes=[0])
        proc=self.start();self.assertEqual(proc.wait(timeout=5),0)
        self.settings()
        again=self.start();self.wait_for(lambda:len(self.starts())==2 and self.ready())
        self.assertIsNone(again.poll());again.terminate();again.wait(timeout=5)

    def test_server_identity_rejects_wrong_app_and_oversized_record(self):
        proc=self.start();state=self.wait_for(self.ready)
        self.assertIsNone(launcher.server_identity(state['port'],self.root/'elsewhere'))
        self.assertIsNone(launcher.server_identity(0,self.root))
        proc.terminate();proc.wait(timeout=5)
        class Response:
            def __enter__(self):return self
            def __exit__(self,*args):pass
            def read(self,limit):self.limit=limit;return b'x'*limit
        response=Response()
        with patch.object(launcher.urllib.request,'build_opener') as build:
            build.return_value.open.return_value=response
            self.assertIsNone(launcher.server_identity(12345,self.root))
            self.assertEqual(response.limit,16385)

    def test_launcher_lock_symlink_rejected(self):
        (self.root/'logs').mkdir()
        victim=self.root/'untouched';victim.write_text('untouched')
        (self.root/'logs/studio_launcher.lock').symlink_to(victim)
        instance=launcher.StudioLauncher(launcher.LauncherConfig(self.root,self.root))
        with self.assertRaises(OSError):instance.acquire()
        self.assertEqual(victim.read_text(),'untouched')

    def test_startup_identity_timeout_stops_owned_child_without_retry(self):
        self.settings(startup_hang=True,launcher_timeout=.2)
        proc=self.start();self.assertEqual(proc.wait(timeout=5),1)
        self.assertEqual(self.state()['event'],'LAUNCHER_STOPPED')
        self.assertIn('no automatic retry',self.state()['reason'])
        spawned=[row for row in self.events() if row['event']=='SERVER_SPAWNED']
        self.assertEqual(len(spawned),1)
        with self.assertRaises(ProcessLookupError):os.kill(spawned[0]['server_pid'],0)

    def test_stop_during_backoff_cancels_pending_restart(self):
        self.settings(exit_codes=[7],launcher_backoff=[2,2])
        proc=self.start();self.wait_for(lambda:self.state().get('event')=='SERVER_RESTART_SCHEDULED')
        proc.terminate();self.assertEqual(proc.wait(timeout=5),143)
        self.assertEqual(len(self.starts()),1)
        self.assertEqual(self.state()['status'],'stopped')

    def test_shell_launch_without_no_open_works_on_macos_bash(self):
        # Bash 3.2 treats an empty array expansion as unbound with set -u.
        envbin=self.root/'.qc_env/bin';envbin.mkdir(parents=True)
        python=envbin/'python';python.write_text('#!/bin/bash\nprintf "%s\\n" "$@"\n')
        python.chmod(0o755)
        refs=self.root/'references/GCF_000001405.40';refs.mkdir(parents=True)
        (refs/'reference_manifest.json').write_text('{}')
        script=self.root/'start_studio.sh';script.write_text((BASE/'start_studio.sh').read_text())
        proc=subprocess.run(['bash',str(script)],capture_output=True,text=True)
        self.assertEqual(proc.returncode,0,proc.stderr)
        self.assertIn('studio_launcher.py',proc.stdout)
        self.assertNotIn('--no-open',proc.stdout)
        again=subprocess.run(['bash',str(script),'--no-open','--port','0'],capture_output=True,text=True)
        self.assertEqual(again.returncode,0,again.stderr)
        self.assertIn('--no-open',again.stdout)
        self.assertIn('--port\n0',again.stdout)

    def test_script_uses_existing_environment_and_supervisor(self):
        script=(BASE/'start_studio.sh').read_text()
        self.assertIn('studio_launcher.py',script)
        self.assertNotIn('RASERFISH_Design_Studio_v1_5',script)
        self.assertNotIn('curl ',script);self.assertNotIn('pip install',script)
        self.assertIn('set -- "$QC_ROOT/.qc_env/bin/python"',script)
        self.assertEqual(subprocess.run(['bash','-n',str(BASE/'start_studio.sh')]).returncode,0)


if __name__=='__main__':unittest.main()
