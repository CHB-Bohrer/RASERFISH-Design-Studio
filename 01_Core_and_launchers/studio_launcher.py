#!/usr/bin/env python3
"""Keep the local Studio server available after at most two unexpected exits.

This launcher never starts a design itself. The server alone decides whether an
explicitly enrolled checkpoint chain can resume. Ctrl-C and intentional signals
stop the one child owned by this launcher and never cause a restart.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import fcntl
import json
import os
from pathlib import Path
import shlex
import signal
import socket
import stat
import subprocess
import sys
import threading
import time
import urllib.error
import urllib.request
import webbrowser

from process_diagnostics import atomic_json, append_event, bounded_tail, returncode_fields, utc_now

BASE = Path(__file__).resolve().parent
STOP_SIGNALS = (signal.SIGINT, signal.SIGTERM, signal.SIGHUP)
INTENTIONAL_EXITS = {0, *(128 + s for s in STOP_SIGNALS), *(-s for s in STOP_SIGNALS)}


@dataclass
class LauncherConfig:
    app_root: Path
    qc_root: Path
    port: int = 8769
    no_open: bool = False
    max_restarts: int = 2
    backoff: tuple[float, ...] = (5.0, 15.0)
    poll_seconds: float = 0.2
    startup_timeout: float = 20.0
    stop_grace: float = 25.0


def server_identity(port, app_root, timeout=0.5):
    """Read only the bounded, token-free identity of this exact installed app."""
    if not isinstance(port, int) or not 0 < port < 65536:
        return None
    try:
        opener = urllib.request.build_opener(urllib.request.ProxyHandler({}))
        with opener.open(f'http://127.0.0.1:{port}/api/server', timeout=timeout) as response:
            raw = response.read(16385)
        if len(raw) > 16384:
            return None
        record = json.loads(raw)
        if (not isinstance(record, dict) or record.get('port') != port or
                not isinstance(record.get('pid'), int) or record['pid'] <= 0 or
                not isinstance(record.get('app_path'), str) or
                Path(record['app_path']).resolve() != Path(app_root).resolve()):
            return None
        return record
    except (OSError, ValueError, urllib.error.URLError):
        return None


def candidate_ports(app_root, requested):
    """Use this app's own recent records; never enumerate other applications."""
    ports = [requested, 8769]
    logs = Path(app_root) / 'logs'
    try:
        with (logs / 'studio_launcher.json').open('rb') as handle:
            raw = handle.read(16385)
        if len(raw) <= 16384:
            record = json.loads(raw)
            if record.get('app_path') == str(Path(app_root).resolve()):
                ports.append(record.get('port'))
    except (OSError, ValueError, AttributeError):
        pass
    for line in reversed(bounded_tail(logs / 'server_lifecycle.jsonl', 65536).splitlines()):
        try:
            record = json.loads(line)
            if (record.get('event') == 'SERVER_STARTED' and
                    record.get('app_path') == str(Path(app_root).resolve())):
                ports.append(record.get('port'))
                break
        except (ValueError, AttributeError):
            continue
    return list(dict.fromkeys(port for port in ports if isinstance(port, int) and 0 < port < 65536))


def choose_port(requested):
    if not 0 <= requested < 65536:
        raise ValueError('Local port must be between 0 and 65535.')
    with socket.socket() as handle:
        # HTTPServer also enables address reuse. TIME_WAIT is not a live server.
        handle.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        try:
            handle.bind(('127.0.0.1', requested))
        except OSError:
            if requested == 0:
                raise
            handle.bind(('127.0.0.1', 0))
        return handle.getsockname()[1]


class StudioLauncher:
    def __init__(self, config):
        self.config = config
        self.app = Path(config.app_root).resolve()
        self.qc_root = Path(config.qc_root).resolve()
        self.logs = self.app / 'logs'
        self.state_path = self.logs / 'studio_launcher.json'
        self.event_path = self.logs / 'studio_launcher.jsonl'
        self.stop_event = threading.Event()
        self.stop_signal = None
        self.child = None
        self.lock_fd = None
        self.port = None
        self.record = dict(schema='raserfish.studio-launcher.v1', app_path=str(self.app),
                           qc_root=str(self.qc_root), launcher_pid=os.getpid(),
                           started_utc=utc_now(), max_restarts=config.max_restarts,
                           restarts=0, server_pid=None)

    def journal(self, event, **fields):
        self.record.update(fields, event=event, updated_utc=utc_now())
        atomic_json(self.state_path, self.record)
        append_event(self.event_path, self.record)

    def acquire(self):
        self.logs.mkdir(parents=True, exist_ok=True)
        path = self.logs / 'studio_launcher.lock'
        fd = os.open(path, os.O_CREAT | os.O_RDWR | getattr(os, 'O_NOFOLLOW', 0), 0o600)
        info = os.fstat(fd)
        if not stat.S_ISREG(info.st_mode) or info.st_uid != os.getuid() or info.st_nlink != 1:
            os.close(fd)
            raise RuntimeError('Unsafe Studio launcher lock file.')
        try:
            fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            os.close(fd)
            return False
        self.lock_fd = fd
        os.ftruncate(fd, 0)
        os.write(fd, json.dumps(dict(pid=os.getpid(), app_path=str(self.app))).encode())
        os.fsync(fd)
        return True

    def existing(self):
        for port in candidate_ports(self.app, self.config.port):
            identity = server_identity(port, self.app)
            if identity:
                return identity
        return None

    def show_existing(self, identity):
        address = f"http://127.0.0.1:{identity['port']}/"
        print(f'Studio is already open: {address}\nUsing the existing server; no duplicate was started.', flush=True)
        if not self.config.no_open:
            webbrowser.open(address)

    def request_stop(self, signum, frame=None):
        # Signal handlers only set state. Durable logging and owned-child cleanup
        # happen in the main loop; a second signal cannot interrupt that cleanup.
        if self.stop_signal is None:
            self.stop_signal = signum
        self.stop_event.set()

    def stop_owned_child(self, reason):
        child = self.child
        if child is None or child.poll() is not None:
            return
        self.journal('OWNED_SERVER_STOP_REQUESTED', reason=reason, server_pid=child.pid,
                     signal='SIGTERM')
        try:
            child.send_signal(signal.SIGTERM)
        except ProcessLookupError:
            return
        try:
            child.wait(timeout=self.config.stop_grace)
        except subprocess.TimeoutExpired:
            self.journal('OWNED_SERVER_FORCE_STOP_REQUESTED', reason=reason,
                         server_pid=child.pid, signal='SIGKILL')
            try:
                child.kill()
            except ProcessLookupError:
                pass
            child.wait(timeout=5)

    def wait_ready(self):
        deadline = time.monotonic() + self.config.startup_timeout
        while not self.stop_event.is_set() and self.child.poll() is None:
            identity = server_identity(self.port, self.app)
            if identity and identity['pid'] == self.child.pid:
                return True
            if time.monotonic() >= deadline:
                return False
            self.stop_event.wait(self.config.poll_seconds)
        return False

    def supervise(self):
        self.port = choose_port(self.config.port)
        address = f'http://127.0.0.1:{self.port}/'
        command = [sys.executable, '-u', str(self.app / 'studio.py'), '--qc-root',
                   str(self.qc_root), '--port', str(self.port), '--no-open']
        log_path = self.logs / (time.strftime('studio_auto_%Y%m%d_%H%M%S_') + str(os.getpid()) + '.log')
        self.record.update(port=self.port, address=address, argv=command,
                           command_display=shlex.join(command), console_log=str(log_path))
        self.journal('LAUNCHER_STARTED', status='starting')
        print(f'Studio address: {address}\nConsole log: {log_path}\n'
              'Keep this Terminal open. Stop with Ctrl-C. Unexpected server exits retry at most twice.', flush=True)
        opened = False
        # Child output is written directly to disk. There is no unread pipe or
        # retained stream whose size can grow in launcher memory.
        with log_path.open('ab', buffering=0) as console:
            for attempt in range(self.config.max_restarts + 1):
                if self.stop_event.is_set():
                    break
                self.journal('SERVER_START_REQUESTED', status='starting', restarts=attempt)
                env = os.environ.copy()
                env['PYTHONNOUSERSITE'] = '1'
                env['PYTHONDONTWRITEBYTECODE'] = '1'
                env['RASERFISH_SUPERVISED_SERVER'] = '1'
                env['PATH'] = str(self.qc_root / '.qc_env/bin') + os.pathsep + env.get('PATH', '')
                env.pop('PYTHONPATH', None)
                env.pop('PYTHONHOME', None)
                self.child = subprocess.Popen(command, cwd=self.app, env=env,
                                              stdin=subprocess.DEVNULL, stdout=console, stderr=console,
                                              start_new_session=True, close_fds=True)
                self.journal('SERVER_SPAWNED', server_pid=self.child.pid)
                ready = self.wait_ready()
                if ready:
                    self.journal('SERVER_READY', status='running')
                    if not opened and not self.config.no_open:
                        webbrowser.open(address)
                        opened = True
                    while self.child.poll() is None and not self.stop_event.wait(self.config.poll_seconds):
                        pass
                elif self.child.poll() is None and not self.stop_event.is_set():
                    self.journal('SERVER_STARTUP_UNCONFIRMED', status='failed',
                                 reason='The owned server did not confirm its identity before the startup timeout.')
                    self.stop_owned_child('Server startup identity was not confirmed.')
                    # A nonresponsive startup is not a verified crash; do not keep
                    # killing and restarting it while it could be doing work.
                    self.journal('LAUNCHER_STOPPED', status='failed',
                                 reason='Server startup identity was not confirmed; no automatic retry.')
                    return 1
                if self.stop_event.is_set():
                    self.stop_owned_child('Launcher received ' + signal.Signals(self.stop_signal).name)
                    break
                code = self.child.wait()
                self.journal('SERVER_EXITED', **returncode_fields(code))
                if code in INTENTIONAL_EXITS:
                    self.journal('LAUNCHER_STOPPED', status='stopped',
                                 reason='Server exited normally or received an intentional stop signal.')
                    return code if code >= 0 else 128 - code
                if attempt >= self.config.max_restarts:
                    self.journal('RESTART_LIMIT_REACHED', status='failed',
                                 reason='Automatic server restart limit reached; saved analysis remains incomplete.')
                    print('Studio stopped after repeated failures. Saved analysis remains incomplete. '
                          f'Details: {log_path}', flush=True)
                    return code if code > 0 else 1
                delay = self.config.backoff[min(attempt, len(self.config.backoff) - 1)]
                self.journal('SERVER_RESTART_SCHEDULED', status='retry_pending', delay_seconds=delay,
                             reason='Unexpected server exit; the next server verifies saved progress before any recovery.')
                print(f'Studio exited unexpectedly; restarting in {delay:g} seconds at {address}', flush=True)
                if self.stop_event.wait(delay):
                    break
        self.journal('LAUNCHER_STOPPED', status='stopped',
                     reason='Launcher received ' + signal.Signals(self.stop_signal).name,
                     signal=signal.Signals(self.stop_signal).name)
        return 128 + self.stop_signal

    def run(self):
        old_handlers = {}
        if not self.acquire():
            identity = self.existing()
            if identity:
                self.show_existing(identity)
            else:
                print('The Studio launcher is already running and may be restarting its server. '
                      'No duplicate was started.', flush=True)
            return 0
        try:
            identity = self.existing()
            if identity:
                self.show_existing(identity)
                return 0
            for sig in STOP_SIGNALS:
                old_handlers[sig] = signal.signal(sig, self.request_stop)
            return self.supervise()
        except BaseException as exc:
            try:
                self.journal('LAUNCHER_ERROR', status='failed', error=str(exc),
                             exception_type=type(exc).__name__)
                self.stop_owned_child('Launcher error: ' + str(exc))
            except Exception as cleanup_error:
                print('Studio launcher cleanup error: ' + str(cleanup_error), file=sys.stderr, flush=True)
            raise
        finally:
            for sig, old in old_handlers.items():
                signal.signal(sig, old)
            if self.lock_fd is not None:
                os.close(self.lock_fd)
                self.lock_fd = None


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--qc-root', type=Path, default=BASE)
    parser.add_argument('--port', type=int, default=8769)
    parser.add_argument('--no-open', action='store_true')
    args = parser.parse_args(argv)
    if not 0 <= args.port < 65536:
        parser.error('--port must be between 0 and 65535')
    return StudioLauncher(LauncherConfig(BASE, args.qc_root, args.port, args.no_open)).run()


if __name__ == '__main__':
    raise SystemExit(main())
