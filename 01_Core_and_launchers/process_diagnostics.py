"""Bounded subprocess diagnostics. Command output goes to disk, never PIPE buffers.

A diagnostic with status other than COMPLETE is not evidence of a successful
search. Only the owning Popen child (or its newly created process group) is
signalled here. Supervisor interruption is recorded separately from an observed
child return code: a requested signal is never reported as a known exit cause.
"""
from __future__ import annotations
import json
import os
from pathlib import Path
import shlex
import shutil
import signal
import subprocess
import time
import uuid
from datetime import datetime, timezone

TAIL_BYTES = 4096


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def atomic_json(path, payload):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name('.' + path.name + '.' + uuid.uuid4().hex + '.tmp')
    try:
        with temporary.open('w', encoding='utf-8') as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write('\n')
            handle.flush()
            os.fsync(handle.fileno())
        temporary.replace(path)
        directory = os.open(str(path.parent), os.O_RDONLY)
        try: os.fsync(directory)
        finally: os.close(directory)
    finally:
        if temporary.exists(): temporary.unlink()


def append_event(path, event):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('a', encoding='utf-8') as handle:
        handle.write(json.dumps(dict(timestamp=utc_now(), **event), sort_keys=True) + '\n')
        handle.flush()
        os.fsync(handle.fileno())


def bounded_tail(path, maximum=TAIL_BYTES):
    try:
        with Path(path).open('rb') as handle:
            handle.seek(0, os.SEEK_END)
            handle.seek(max(0, handle.tell() - maximum))
            return handle.read(maximum).decode('utf-8', errors='replace')
    except (FileNotFoundError, IsADirectoryError):
        return ''


def command_record(command, cwd=None, *, environment=None):
    requested = [str(value) for value in command]
    if not requested: raise ValueError('Cannot execute an empty command.')
    working_directory = str(Path(cwd or os.getcwd()).resolve())
    executable = shutil.which(requested[0], path=(environment if environment is not None else os.environ).get('PATH'))
    if executable is None and '/' in requested[0]:
        candidate = Path(requested[0])
        executable = str((Path(working_directory) / candidate).absolute())
    argv = list(requested)
    if executable:
        executable = str(Path(executable).absolute())
        argv[0] = executable
    return dict(requested_argv=requested, argv=argv, executable=executable,
                command_display=shlex.join(argv), cwd=working_directory,
                started_utc=utc_now(), launcher_pid=os.getpid(), pid=None, pgid=None,
                returncode=None, terminating_signal=None, termination_actions=[],
                status='STARTING', analysis_complete=False)


def returncode_fields(code):
    if code is None: return dict(returncode=None, terminating_signal=None)
    name = None
    if code < 0:
        try: name = signal.Signals(-code).name
        except ValueError: name = 'SIGNAL_' + str(-code)
    return dict(returncode=code, terminating_signal=name)


def exit_description(record):
    if record.get('spawn_error'):
        return 'could not start: ' + record['spawn_error']
    code = record.get('returncode')
    if code is None:
        return 'did not have a verified successful exit (child return code unavailable)'
    if code < 0:
        return 'terminated by ' + str(record.get('terminating_signal')) + f' (return code {code})'
    return f'exited with return code {code}'


def failure_message(record, diagnostic_path):
    name = Path(record.get('executable') or record['argv'][0]).name
    message = name + ' ' + exit_description(record) + '.'
    if record.get('stop_reason'): message += ' Stop reason: ' + record['stop_reason']
    if record.get('stderr_tail'): message += '\nStderr: ' + record['stderr_tail'].strip()
    elif record.get('stdout_tail'): message += '\nStdout: ' + record['stdout_tail'].strip()
    else: message += ' No stdout/stderr was emitted.'
    return message + '\nDetails: ' + str(diagnostic_path)


def write_summary(path, record, diagnostic_path):
    text = (f"Status: {record['status']}\nCommand: {record['command_display']}\n"
            f"Working directory: {record['cwd']}\nPID: {record.get('pid')}  PGID: {record.get('pgid')}\n"
            f"Started: {record['started_utc']}\nEnded: {record.get('ended_utc', 'not yet finished')}\n"
            f"Outcome: {exit_description(record)}\n")
    if record.get('environment_overrides'): text += 'Environment overrides: ' + json.dumps(record['environment_overrides'], sort_keys=True) + '\n'
    if record.get('stop_reason'): text += 'Stop reason: ' + record['stop_reason'] + '\n'
    text += f"Stdout: {record.get('stdout_path')}\nStderr: {record.get('stderr_path')}\nDiagnostics: {diagnostic_path}\n"
    if record.get('stderr_tail'): text += '\nStderr tail:\n' + record['stderr_tail']
    if record.get('stdout_tail'): text += '\nStdout tail:\n' + record['stdout_tail']
    Path(path).write_text(text, encoding='utf-8')


def _stop_owned_child(proc, own_group, record, diagnostic, reason):
    for sig, timeout in ((signal.SIGTERM, 3), (signal.SIGKILL, 3)):
        if proc.poll() is not None: break
        action = dict(timestamp=utc_now(), actor_pid=os.getpid(), target_pid=proc.pid,
                      target_pgid=proc.pid if own_group else None,
                      signal=signal.Signals(sig).name, reason=reason, action='SIGNAL_REQUESTED')
        record['termination_actions'].append(action)
        try: atomic_json(diagnostic, record)  # durable BEFORE signalling
        except Exception as io_error:
            # If the journal itself is unwritable, stop instead of allowing an
            # unmonitored child to continue. The supervisor stderr is the fallback.
            record['diagnostic_write_error']=str(io_error)
            fallback=dict(action,diagnostic_write_error=str(io_error),diagnostic_path=str(diagnostic))
            try: os.write(2,(json.dumps(fallback)+'\n').encode('utf-8'))
            except OSError: pass
        try:
            if own_group: os.killpg(proc.pid, sig)
            else: proc.send_signal(sig)
        except ProcessLookupError: pass
        try: proc.wait(timeout=timeout)
        except subprocess.TimeoutExpired: continue
    if proc.poll() is None:
        raise RuntimeError('Owned subprocess did not stop after SIGKILL; analysis incomplete.')


def run_process(command, outlog, cancel=None, resource_guard=None, *, cancellation_class=RuntimeError, env_overrides=None):
    """Run one owned child, preserving precise failures even with empty stderr.

    Explicit environment overrides apply only to the child. Journals retain
    only those explicit values, never the rest of the inherited environment.
    """
    outlog = Path(outlog)
    outlog.parent.mkdir(parents=True, exist_ok=True)
    diagnostic = Path(str(outlog) + '.process.json')
    stdout_path, stderr_path = Path(str(outlog) + '.stdout.log'), Path(str(outlog) + '.stderr.log')
    explicit_environment = dict(env_overrides or {})
    child_environment = os.environ.copy()
    child_environment.update(explicit_environment)
    record = command_record(command, environment=child_environment)
    record['environment_overrides'] = explicit_environment
    own_group = os.environ.get('RASERFISH_DESIGN_WORKER') != '1'
    record.update(stdout_path=str(stdout_path), stderr_path=str(stderr_path),
                  summary_path=str(outlog), owns_process_group=own_group)
    atomic_json(diagnostic, record)
    write_summary(outlog, record, diagnostic)
    proc = None
    caught = None
    start = time.monotonic()
    try:
        lock_fd = os.environ.get('RASERFISH_DESIGN_LOCK_FD')
        inherited_fds = (int(lock_fd),) if lock_fd is not None else ()
        with stdout_path.open('wb') as stdout, stderr_path.open('wb') as stderr:
            proc = subprocess.Popen(record['argv'], cwd=record['cwd'], stdout=stdout, stderr=stderr,
                                    start_new_session=own_group, pass_fds=inherited_fds, env=child_environment)
            # Non-integer PIDs occur only in explicitly mocked subprocess tests.
            record['pid'] = proc.pid if isinstance(proc.pid, int) else None
            try: record['pgid'] = os.getpgid(proc.pid) if record['pid'] is not None else None
            except ProcessLookupError: record['pgid'] = proc.pid if own_group else os.getpgrp()
            record['status'] = 'RUNNING'
            atomic_json(diagnostic, record)
            write_summary(outlog, record, diagnostic)
            try:
                while proc.poll() is None:
                    if cancel and cancel.is_set():
                        raise cancellation_class('Cancelled by user. Analysis is incomplete; partial files retained.')
                    if resource_guard: resource_guard.checkpoint('BLAST subprocess')
                    time.sleep(.2)
            except BaseException as exc:
                caught = exc
                record.update(status='STOP_REQUESTED', stop_reason=str(exc),
                              stop_exception=type(exc).__name__, memory_sample=getattr(exc, 'sample', None))
                atomic_json(diagnostic, record)
                _stop_owned_child(proc, own_group, record, diagnostic, str(exc))
    except BaseException as exc:
        if caught is None: caught = exc
        if proc is None: record['spawn_error'] = f'{type(exc).__name__}: {exc}'
        elif proc.poll() is None:
            # This includes journal failures immediately after Popen, before the
            # normal wait loop has begun. Never leak a child on diagnostic I/O failure.
            record.update(status='STOP_REQUESTED')
            record.setdefault('stop_reason',str(exc));record.setdefault('stop_exception',type(exc).__name__)
            _stop_owned_child(proc, own_group, record, diagnostic, record['stop_reason'])
    finally:
        record.update(returncode_fields(proc.returncode if proc is not None else None))
        record.update(ended_utc=utc_now(), elapsed_seconds=time.monotonic()-start,
                      stdout_tail=bounded_tail(stdout_path), stderr_tail=bounded_tail(stderr_path))
        record['stdout_bytes'] = stdout_path.stat().st_size if stdout_path.exists() else 0
        record['stderr_bytes'] = stderr_path.stat().st_size if stderr_path.exists() else 0
        record['status'] = 'COMPLETE' if caught is None and record['returncode'] == 0 else 'INCOMPLETE'
        record['analysis_complete'] = record['status'] == 'COMPLETE'
        atomic_json(diagnostic, record)
        write_summary(outlog, record, diagnostic)
    if caught is not None:
        # Preserve cancellation and resource-guard exception identity for callers.
        if hasattr(caught, 'add_note'): caught.add_note(failure_message(record, diagnostic))
        if record.get('stop_reason'): raise caught
        raise RuntimeError(failure_message(record, diagnostic)) from caught
    if record['returncode'] != 0: raise RuntimeError(failure_message(record, diagnostic))
    return record


def finalize_interrupted_commands(root, reason, actions=None):
    """Finalize only small, still-running journals after the owner group stops.

    A supervisor cannot reap its grandchild, so the grandchild's return code
    remains unknown instead of attributing SIGKILL or a memory kill without proof.
    """
    finalized = []
    for path in Path(root).rglob('*.process.json'):
        if path.stat().st_size > 1024 * 1024: continue
        try: record = json.loads(path.read_text())
        except (OSError, ValueError): continue
        if record.get('status') not in ('STARTING', 'RUNNING', 'STOP_REQUESTED'): continue
        record.update(status='INCOMPLETE', analysis_complete=False, ended_utc=utc_now(),
                      stop_reason=reason, supervisor_interrupted=True,
                      child_returncode_observed=False, supervisor_actions=actions or [],
                      stdout_tail=bounded_tail(record.get('stdout_path', '')),
                      stderr_tail=bounded_tail(record.get('stderr_path', '')))
        atomic_json(path, record)
        if record.get('summary_path'): write_summary(record['summary_path'], record, path)
        finalized.append(str(path))
    return finalized
