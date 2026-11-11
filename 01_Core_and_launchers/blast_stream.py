"""Bounded Standard-mode BLAST output streaming with verified child completion.

Only a private FIFO transports alignment bytes. The caller still owns parsing,
hashing, exact scoring and the final transaction. Iterator exhaustion means the
entire output was consumed AND the supplied owned-process runner verified exit
zero. Any partial parse, subprocess failure or cancellation is incomplete.
"""
from __future__ import annotations

from contextlib import contextmanager
import os
from pathlib import Path
import select
import stat
import threading


READ_CHUNK_BYTES = 64 * 1024
MAX_BUFFER_BYTES = 1024 * 1024
POLL_SECONDS = 0.05
# process_diagnostics allows up to three seconds each for TERM and KILL.
JOIN_TIMEOUT_SECONDS = 8.0


class _CombinedCancellation:
    def __init__(self, external, internal):
        self.external = external
        self.internal = internal

    def is_set(self):
        return self.internal.is_set() or bool(self.external and self.external.is_set())


@contextmanager
def stream_command(command, outlog, output_path, cancel=None, resource_guard=None,
                   *, run_process, cancellation_class=RuntimeError):
    """Yield LF-terminated raw lines from one monitored owned subprocess.

    ``command`` must already name ``output_path`` as its output destination.
    ``run_process`` has the design-engine runner signature and must return its
    successful process diagnostic (status COMPLETE, returncode 0), or raise.
    Diagnostics and stdout/stderr remain handled by that existing runner.
    """
    output_path = Path(output_path)
    internal_cancel = threading.Event()
    done = threading.Event()
    combined = _CombinedCancellation(cancel, internal_cancel)
    outcome = {}
    state = {'complete': False}
    descriptor = None
    thread = None
    fifo_identity = None
    iterator = None
    primary_error = None

    def execute():
        try:
            outcome['record'] = run_process(command, outlog, combined, resource_guard)
        except BaseException as error:
            outcome['error'] = error
        finally:
            done.set()

    def check_completion():
        if 'error' in outcome:
            raise outcome['error']
        record = outcome.get('record')
        if not isinstance(record, dict) or record.get('status') != 'COMPLETE' or record.get('returncode') != 0:
            raise RuntimeError('Streamed BLAST did not have a verified successful exit; analysis incomplete.')

    def lines():
        buffer = bytearray()
        while True:
            if cancel and cancel.is_set():
                raise cancellation_class('Cancelled while reading BLAST output; analysis incomplete.')
            # Snapshot completion BEFORE reading. If the producer finishes
            # between an empty read and the next check, repeat the read; bytes
            # written in that interval must not be mistaken for empty output.
            completed_before_read = done.is_set()
            raw = None
            remaining = MAX_BUFFER_BYTES - len(buffer)
            if remaining <= 0:
                raise RuntimeError('BLAST output line exceeds the 1 MiB streaming buffer; analysis incomplete.')
            try:
                raw = os.read(descriptor, min(READ_CHUNK_BYTES, remaining))
            except BlockingIOError:
                pass
            if raw:
                buffer.extend(raw)
                while True:
                    end = buffer.find(b'\n')
                    if end < 0:
                        break
                    line = bytes(buffer[:end + 1])
                    del buffer[:end + 1]
                    yield line
                if len(buffer) >= MAX_BUFFER_BYTES:
                    raise RuntimeError('BLAST output line exceeds the 1 MiB streaming buffer; analysis incomplete.')
                continue
            if completed_before_read:
                # FIFO data is drained before inspecting the result. A child
                # that exits before opening the FIFO also reaches this branch.
                check_completion()
                if buffer:
                    raise RuntimeError('BLAST output ended with an unterminated line; analysis incomplete.')
                state['complete'] = True
                return
            if raw == b'':
                # No writer yet, or the writer closed before its diagnostics
                # finished. EOF alone is never evidence of successful exit.
                done.wait(POLL_SECONDS)
            else:
                select.select([descriptor], [], [], POLL_SECONDS)

    try:
        # Refuse to replace any existing path. Only this exact FIFO inode may
        # be removed on exit, even if another actor replaces its pathname.
        os.mkfifo(output_path, 0o600)
        info = output_path.lstat()
        fifo_identity = (info.st_dev, info.st_ino)
        descriptor = os.open(output_path, os.O_RDONLY | os.O_NONBLOCK)
        thread = threading.Thread(target=execute, name='raserfish-blast-stream', daemon=True)
        thread.start()
        iterator = lines()
        yield iterator
        if not state['complete']:
            raise RuntimeError('BLAST output was not fully consumed; analysis incomplete.')
    except BaseException as error:
        primary_error = error
        raise
    finally:
        internal_cancel.set()
        cleanup_errors = []
        if iterator is not None:
            try:
                iterator.close()
            except BaseException as error:
                cleanup_errors.append('Closing output iterator: ' + str(error))
        if descriptor is not None:
            try:
                os.close(descriptor)
            except OSError as error:
                cleanup_errors.append('Closing FIFO reader: ' + str(error))
        if thread is not None and thread.ident is not None:
            thread.join(JOIN_TIMEOUT_SECONDS)
            if thread.is_alive():
                cleanup_errors.append('Owned BLAST runner did not stop within the bounded cleanup interval; analysis incomplete.')
        if fifo_identity is not None:
            try:
                info = output_path.lstat()
                if (info.st_dev, info.st_ino) == fifo_identity and stat.S_ISFIFO(info.st_mode):
                    output_path.unlink()
                else:
                    cleanup_errors.append('Output pathname was replaced; the replacement was preserved.')
            except FileNotFoundError:
                pass
            except OSError as error:
                cleanup_errors.append('Removing FIFO: ' + str(error))
        if cleanup_errors:
            message = '\n'.join(cleanup_errors)
            if primary_error is not None:
                if hasattr(primary_error, 'add_note'):
                    primary_error.add_note(message)
            else:
                raise RuntimeError(message)
