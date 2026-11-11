"""Local single-design locking and conservative, observable process-memory limits.

Only standard-library facilities are used. The UI may start without these metrics;
calling ``MemoryGuard.check`` before starting a design fails closed if they cannot
be obtained. A supervisor must keep checking every ``POLL_INTERVAL_SECONDS`` and
stop only its own worker process group on failure, including during native calls.
"""
from __future__ import annotations

import ctypes
import errno
import fcntl
from functools import lru_cache
import json
import os
import platform
import re
import stat
import shutil
import tempfile
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path

GIB = 1024 ** 3
POLL_INTERVAL_SECONDS = 5.0
MINIMUM_FREE_DISK_BYTES = GIB


class DesignAlreadyRunning(RuntimeError):
    """Another process holds the shared local design lock."""

    def __init__(self, path, metadata=None):
        self.path = str(path)
        self.metadata = metadata or {}
        details = self.metadata.get('job_id') or self.metadata.get('pid')
        suffix = f' (job/process {details})' if details else ''
        super().__init__('Another design is already active on this Mac' + suffix +
                         '. Wait for it to finish or stop it in its own app before starting another design.')


class MemoryLimitExceeded(RuntimeError):
    """The design must stop incomplete; this never means QC passed."""

    def __init__(self, message, sample=None):
        self.sample = sample or {}
        super().__init__(message)


class DiskSpaceLimitExceeded(MemoryLimitExceeded):
    """Owned work must stop incomplete before filling the output/scratch disk."""


class MemoryMetricsUnavailable(MemoryLimitExceeded):
    """Required memory metrics could not be verified, so execution is blocked."""


def default_lock_path():
    # /tmp resolves to /private/tmp on macOS. Independent of app location or port.
    return Path('/tmp') / f'raserfish-design-{os.getuid()}.lock'


def _read_lock_metadata(fd):
    try:
        value = json.loads(os.pread(fd, 16384, 0).decode('utf-8'))
        return value if isinstance(value, dict) else {}
    except (OSError, ValueError, UnicodeError):
        return {}


def _validate_lock_fd(fd, path):
    info = os.fstat(fd)
    linked = os.stat(path, follow_symlinks=False)
    if (not stat.S_ISREG(info.st_mode) or not stat.S_ISREG(linked.st_mode)
            or info.st_uid != os.getuid() or linked.st_uid != os.getuid()
            or (info.st_dev, info.st_ino) != (linked.st_dev, linked.st_ino)
            or info.st_mode & 0o077):
        raise RuntimeError('Unsafe design lock file: expected a private regular file owned by this user.')


class DesignLock:
    """A lock descriptor transferable through subprocess ``pass_fds``.

    Release closes this process's descriptor; it deliberately does NOT call
    LOCK_UN. A worker inheriting the same open-file description must continue to
    hold the lock if its supervisor exits. Never unlink the shared lock file.
    """

    def __init__(self, fd, path):
        self.fd = fd
        self.path = Path(path)

    @property
    def metadata(self):
        return _read_lock_metadata(self.fd) if self.fd is not None else {}

    def update_metadata(self, metadata=None, **fields):
        if self.fd is None:
            raise RuntimeError('The design lock has already been released.')
        value = self.metadata
        value.update(metadata or {})
        value.update(fields)
        value.setdefault('pid', os.getpid())
        value.setdefault('started_at', datetime.now(timezone.utc).isoformat())
        data = json.dumps(value, sort_keys=True).encode('utf-8')
        if len(data) > 16384:
            raise ValueError('Design lock metadata is too large.')
        os.ftruncate(self.fd, 0)
        os.pwrite(self.fd, data, 0)
        os.fsync(self.fd)
        return value

    def release(self):
        if self.fd is not None:
            fd, self.fd = self.fd, None
            os.close(fd)

    close = release

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.release()


def acquire_design_lock(metadata=None, path=None):
    path = Path(path) if path is not None else default_lock_path()
    flags = os.O_CREAT | os.O_RDWR | getattr(os, 'O_CLOEXEC', 0) | getattr(os, 'O_NOFOLLOW', 0)
    fd = os.open(path, flags, 0o600)
    try:
        _validate_lock_fd(fd, path)
        try:
            fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise DesignAlreadyRunning(path, _read_lock_metadata(fd)) from None
        lock = DesignLock(fd, path)
        # Stale metadata is diagnostic only and must not survive a new acquisition.
        os.ftruncate(fd, 0)
        lock.update_metadata(metadata)
        return lock
    except BaseException:
        os.close(fd)
        raise


def inherited_design_lock(fd, path=None):
    """Adopt the descriptor passed by the supervisor without creating a new lock.

    Validation prevents accepting an unrelated descriptor or a replaced path.
    Calling flock on the inherited open-file description preserves its lock.
    """
    if not isinstance(fd, int) or isinstance(fd, bool) or fd < 0:
        raise ValueError('A valid inherited design lock descriptor is required.')
    path = Path(path) if path is not None else default_lock_path()
    _validate_lock_fd(fd, path)
    try:
        fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        raise DesignAlreadyRunning(path, _read_lock_metadata(fd)) from None
    os.set_inheritable(fd, False)  # BLAST does not need to inherit this descriptor.
    return DesignLock(fd, path)


def _command(args):
    done = subprocess.run(args, check=True, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE, text=True, timeout=3)
    return done.stdout


def parse_vm_stat(text):
    """Parse vm_stat without assuming Intel's page size on Apple Silicon."""
    page = re.search(r'page size of\s+(\d+)\s+bytes', text)
    if not page:
        raise ValueError('vm_stat did not report its page size.')
    page_size = int(page.group(1))
    if not 1024 <= page_size <= 1048576:
        raise ValueError('Invalid vm_stat page size.')
    fields = {}
    for line in text.splitlines():
        match = re.match(r'\s*([^:]+):\s*(\d+)\.?\s*$', line)
        if match:
            fields[match.group(1).strip()] = int(match.group(2))
    required = ('Pages free', 'Pages inactive', 'File-backed pages')
    if any(key not in fields for key in required):
        raise ValueError('vm_stat is missing free, inactive, or file-backed pages.')
    # This is an estimate, not a promise that all reclaimable memory is instantly
    # free. Cap inactive pages by file-backed pages; do not add speculative or
    # purgeable counters, which could double-count reclaimable pages.
    available = (fields['Pages free'] + min(fields['Pages inactive'], fields['File-backed pages'])) * page_size
    return dict(available_bytes=available,
                available_method='macOS conservative estimate: free + min(inactive, file-backed) pages',
                page_size_bytes=page_size,
                vm_pages={key: fields[key] for key in required},
                compressed_bytes=fields.get('Pages occupied by compressor', 0) * page_size)


@lru_cache(maxsize=1)
def _mac_physical_memory_bytes():
    """Read installed RAM once through macOS's public, read-only sysctl API.

    Installed RAM is static for this process's lifetime. Avoid launching a
    sysctl subprocess on every safety poll; cache only a successfully validated
    value. Available memory, process accounting, and disk space remain fresh.
    """
    try:
        function = ctypes.CDLL('/usr/lib/libSystem.B.dylib', use_errno=True).sysctlbyname
    except AttributeError as exc:
        raise RuntimeError('macOS sysctlbyname physical-memory accounting is unavailable.') from exc
    function.argtypes = (ctypes.c_char_p, ctypes.c_void_p,
                         ctypes.POINTER(ctypes.c_size_t), ctypes.c_void_p,
                         ctypes.c_size_t)
    function.restype = ctypes.c_int
    physical = ctypes.c_uint64()
    size = ctypes.c_size_t(ctypes.sizeof(physical))
    ctypes.set_errno(0)
    if function(b'hw.memsize', ctypes.byref(physical), ctypes.byref(size), None, 0) != 0:
        code = ctypes.get_errno() or errno.EIO
        raise OSError(code, f'sysctlbyname hw.memsize: {os.strerror(code)}')
    if size.value != ctypes.sizeof(physical) or physical.value <= 0:
        raise ValueError('Invalid macOS sysctlbyname hw.memsize response.')
    return int(physical.value)


def system_memory_snapshot():
    system = platform.system()
    if system == 'Darwin':
        physical = _mac_physical_memory_bytes()
        result = parse_vm_stat(_command(['/usr/bin/vm_stat']))
    elif system == 'Linux':
        fields = {}
        for line in Path('/proc/meminfo').read_text().splitlines():
            key, _, value = line.partition(':')
            fields[key] = int(value.strip().split()[0]) * 1024
        physical = fields['MemTotal']
        result = dict(available_bytes=fields['MemAvailable'],
                      available_method='Linux /proc/meminfo MemAvailable')
    else:
        raise RuntimeError(f'Memory monitoring is unavailable on {system}.')
    if physical <= 0 or not 0 <= result['available_bytes'] <= physical:
        raise ValueError('Invalid physical/available memory metrics.')
    result['physical_bytes'] = physical
    return result


def parse_process_table(text):
    rows = []
    for line in text.splitlines():
        parts = line.strip().split(None, 3)
        if len(parts) < 3:
            continue
        try:
            pid, ppid, rss_kib = map(int, parts[:3])
        except ValueError:
            continue
        if pid <= 0 or ppid < 0 or rss_kib < 0:
            continue
        rows.append(dict(pid=pid, ppid=ppid, rss_bytes=rss_kib * 1024,
                         name=parts[3] if len(parts) > 3 else ''))
    return rows


class _MacRusageInfoV4(ctypes.Structure):
    """Public macOS sys/resource.h rusage_info_v4 ABI (RUSAGE_INFO_V4=4)."""

    _fields_ = [('ri_uuid', ctypes.c_uint8 * 16)] + [
        (name, ctypes.c_uint64) for name in (
            'ri_user_time', 'ri_system_time', 'ri_pkg_idle_wkups',
            'ri_interrupt_wkups', 'ri_pageins', 'ri_wired_size',
            'ri_resident_size', 'ri_phys_footprint', 'ri_proc_start_abstime',
            'ri_proc_exit_abstime', 'ri_child_user_time', 'ri_child_system_time',
            'ri_child_pkg_idle_wkups', 'ri_child_interrupt_wkups',
            'ri_child_pageins', 'ri_child_elapsed_abstime', 'ri_diskio_bytesread',
            'ri_diskio_byteswritten', 'ri_cpu_time_qos_default',
            'ri_cpu_time_qos_maintenance', 'ri_cpu_time_qos_background',
            'ri_cpu_time_qos_utility', 'ri_cpu_time_qos_legacy',
            'ri_cpu_time_qos_user_initiated', 'ri_cpu_time_qos_user_interactive',
            'ri_billed_system_time', 'ri_serviced_system_time', 'ri_logical_writes',
            'ri_lifetime_max_phys_footprint', 'ri_instructions', 'ri_cycles',
            'ri_billed_energy', 'ri_serviced_energy',
            'ri_interval_max_phys_footprint', 'ri_runnable_time')]


@lru_cache(maxsize=1)
def _mac_proc_pid_rusage():
    # Read-only local API, with no external dependency or privileged attachment.
    try:
        function = ctypes.CDLL('/usr/lib/libproc.dylib', use_errno=True).proc_pid_rusage
    except AttributeError as exc:
        raise RuntimeError('macOS proc_pid_rusage memory accounting is unavailable.') from exc
    function.argtypes = (ctypes.c_int, ctypes.c_int, ctypes.c_void_p)
    function.restype = ctypes.c_int
    return function


def mac_process_memory(pid):
    """Account for charged compressed/swapped pages that ps RSS can omit.

    Public proc_pid_rusage's physical footprint is the kernel memory-accounting
    metric, not virtual address-space size. Lifetime peak is diagnostic only;
    peaks from different processes need not have occurred at the same time.
    """
    usage = _MacRusageInfoV4()
    ctypes.set_errno(0)
    if _mac_proc_pid_rusage()(int(pid), 4, ctypes.byref(usage)) != 0:
        code = ctypes.get_errno() or errno.EIO
        raise OSError(code, f'proc_pid_rusage PID {pid}: {os.strerror(code)}')
    return dict(physical_footprint_bytes=int(usage.ri_phys_footprint),
                lifetime_peak_physical_footprint_bytes=int(usage.ri_lifetime_max_phys_footprint),
                kernel_resident_bytes=int(usage.ri_resident_size),
                process_start_abstime=int(usage.ri_proc_start_abstime))


def _process_lifecycle_snapshot(pid):
    """Recheck one denied descendant; never infer exit from permission errors.

    macOS ps is setuid on some systems. A supervisor can see a worker's short
    lived ps helper while libproc denies its footprint. A fresh lifecycle query
    distinguishes that race from an unmeasurable, still-running design process.
    """
    done = subprocess.run(['/bin/ps', '-p', str(pid), '-o', 'pid=,ppid=,stat=,comm='],
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                          text=True, timeout=3)
    if done.returncode == 1 and not done.stdout.strip() and not done.stderr.strip():
        return dict(pid=pid, state='absent')
    if done.returncode != 0:
        raise RuntimeError(f'Cannot verify lifecycle of PID {pid}: ps returned '
                           f'{done.returncode}: {done.stderr.strip()[:1000]}')
    parts = done.stdout.strip().split(None, 3)
    if (len(parts) != 4 or not parts[0].isdigit() or int(parts[0]) != pid
            or not parts[1].isdigit() or not parts[2]):
        raise RuntimeError(f'Cannot verify lifecycle of PID {pid}: invalid ps response.')
    return dict(pid=pid, ppid=int(parts[1]), process_status=parts[2], name=parts[3],
                state='zombie' if parts[2].startswith('Z') else 'live')


def _recover_child_memory_sample(row, original_error):
    """Bounded retry; a live child with persistently denied metrics still stops.

    Neither names nor zero RSS grant an exemption. Only a fresh absent/zombie
    observation (or libproc ESRCH) permits excluding a descendant from totals.
    """
    pid = row['pid']
    event = dict(pid=pid, ppid=row['ppid'], name=row['name'],
                 initial_errno=original_error.errno, lifecycle_checks=[])
    error = original_error
    for attempt in range(3):
        try:
            state = _process_lifecycle_snapshot(pid)
        except (OSError, RuntimeError, subprocess.SubprocessError) as exc:
            event['outcome'] = 'lifecycle_verification_failed'
            exc.process_snapshot_error = event
            raise
        event['lifecycle_checks'].append(state)
        if state['state'] != 'absent' and state['ppid'] != row['ppid']:
            event['outcome'] = 'process_identity_changed'
            exc = RuntimeError(f'PID {pid} parent changed during memory sampling; '
                               'cannot verify the original descendant.')
            exc.process_snapshot_error = event
            raise exc
        if state['state'] in ('absent', 'zombie'):
            event['outcome'] = 'verified_exited_after_permission_error'
            return None, event
        try:
            memory = mac_process_memory(pid)
        except OSError as exc:
            if exc.errno == errno.ESRCH:
                event['outcome'] = 'verified_exited_during_retry'
                return None, event
            if exc.errno not in (errno.EPERM, errno.EACCES):
                event['outcome'] = 'memory_retry_failed'
                exc.process_snapshot_error = event
                raise
            error = exc
        else:
            event['outcome'] = 'memory_read_recovered'
            return memory, event
        if attempt < 2:
            time.sleep(.02)
    event['outcome'] = 'live_process_memory_denied'
    exc = OSError(error.errno,
                  f'proc_pid_rusage PID {pid} ({row["name"]}, parent {row["ppid"]}): '
                  f'{os.strerror(error.errno)}; descendant remains live after '
                  '3 lifecycle checks; memory safety remains enforced')
    exc.process_snapshot_error = event
    raise exc from original_error


def process_tree_snapshot(root_pid):
    rows = parse_process_table(_command(['/bin/ps', '-axo', 'pid=,ppid=,rss=,comm=']))
    by_pid = {row['pid']: row for row in rows}
    if root_pid not in by_pid:
        raise RuntimeError(f'Design worker PID {root_pid} is absent from the process table.')
    children = {}
    for row in rows:
        children.setdefault(row['ppid'], []).append(row['pid'])
    selected, stack = set(), [root_pid]
    while stack:
        pid = stack.pop()
        if pid in selected:
            continue
        selected.add(pid)
        stack.extend(children.get(pid, ()))
    processes = []
    exited = []
    events = []
    mac = platform.system() == 'Darwin'
    for pid in sorted(selected):
        row = by_pid[pid]
        if mac:
            try:
                row.update(mac_process_memory(pid))
            except OSError as exc:
                # ESRCH is an explicit kernel exit observation. A permission
                # failure alone is not; verify lifecycle and retry narrowly.
                if pid != root_pid and exc.errno == errno.ESRCH:
                    exited.append(pid)
                    continue
                if pid != root_pid and exc.errno in (errno.EPERM, errno.EACCES):
                    memory, event = _recover_child_memory_sample(row, exc)
                    events.append(event)
                    if memory is None:
                        exited.append(pid)
                        continue
                    row.update(memory)
                else:
                    raise
        processes.append(row)
    result = dict(root_pid=root_pid, process_count=len(processes),
                  process_family_rss_bytes=sum(row['rss_bytes'] for row in processes),
                  processes=processes)
    if mac:
        result.update(process_family_physical_footprint_bytes=sum(
            row['physical_footprint_bytes'] for row in processes),
            physical_footprint_method='macOS libproc proc_pid_rusage RUSAGE_INFO_V4; includes charged compressed/swapped memory',
            processes_exited_during_snapshot=exited,
            process_snapshot_events=events)
    return result


def memory_limits(physical_bytes):
    """Safety limits, independent of all scientific screening thresholds."""
    if not isinstance(physical_bytes, (int, float)) or physical_bytes <= 0:
        raise ValueError('Physical memory must be known before deriving safety limits.')
    return dict(available_memory_enforcement='disabled_by_user',
                maximum_process_family_rss_bytes=min(8 * GIB, int(physical_bytes * .50)),
                maximum_process_family_physical_footprint_bytes=min(8 * GIB, int(physical_bytes * .50)),
                poll_interval_seconds=POLL_INTERVAL_SECONDS)


def _gib(value):
    return f'{value / GIB:.2f} GiB'


class MemoryGuard:
    """Log one JSONL observation and enforce limits at every check.

    ``root_pid`` should be the dedicated design worker, so native BLAST children
    are included. Start with ``check('preflight', starting=True)`` before expensive
    work. The supervisor must then call ``check`` at most five seconds apart;
    cooperative engine checkpoints alone cannot interrupt a long native call.
    """

    def __init__(self, log_path, root_pid=None, logger=None, log_callback=None):
        self.log_path = Path(log_path)
        self.root_pid = os.getpid() if root_pid is None else int(root_pid)
        self.logger = logger if logger is not None else log_callback
        self.limits = None
        self.last_sample = None
        self._last_checked_at = None

    def sample(self):
        record = dict(timestamp=datetime.now(timezone.utc).isoformat(),
                      monotonic_seconds=time.monotonic(), root_pid=self.root_pid)
        errors = []
        try:
            record.update(system_memory_snapshot())
        except (OSError, ValueError, KeyError, RuntimeError, subprocess.SubprocessError) as exc:
            errors.append('system memory: ' + str(exc))
        try:
            record.update(process_tree_snapshot(self.root_pid))
            if platform.system() == 'Darwin' and 'process_family_physical_footprint_bytes' not in record:
                raise ValueError('macOS physical-footprint metric is missing.')
        except (OSError, ValueError, KeyError, RuntimeError, subprocess.SubprocessError) as exc:
            errors.append('process family memory: ' + str(exc))
            if hasattr(exc, 'process_snapshot_error'):
                record['process_snapshot_error'] = exc.process_snapshot_error
        record['metrics_available'] = not errors
        if errors:
            record['metric_errors'] = errors
        if 'physical_bytes' in record:
            self.limits = memory_limits(record['physical_bytes'])
            record['limits'] = dict(self.limits)
        return record

    def check(self, stage='working', starting=False):
        record = self.sample()
        record.update(stage=str(stage), starting=bool(starting), status='OK')
        error_type = MemoryLimitExceeded
        message = None
        if not record['metrics_available']:
            record['status'] = 'METRICS_UNAVAILABLE'
            error_type = MemoryMetricsUnavailable
            message = ('Design stopped incomplete: memory safety metrics are unavailable. ' +
                       '; '.join(record.get('metric_errors', [])))
        else:
            limits = record['limits']
            # Available memory is diagnostic only, at the user's request. A
            # separate running reserve would recreate the removed startup block.
            if record['process_family_rss_bytes'] > limits['maximum_process_family_rss_bytes']:
                record['status'] = 'STOP_PROCESS_MEMORY_LIMIT'
                message = ('Design stopped incomplete: worker and child process RSS is ' +
                           _gib(record['process_family_rss_bytes']) + ', above the ' +
                           _gib(limits['maximum_process_family_rss_bytes']) + ' memory limit.')
            elif record.get('process_family_physical_footprint_bytes', 0) > limits['maximum_process_family_physical_footprint_bytes']:
                record['status'] = 'STOP_PROCESS_PHYSICAL_FOOTPRINT_LIMIT'
                message = ('Design stopped incomplete: worker and child process physical memory footprint ' +
                           '(including charged compressed/swapped memory) is ' +
                           _gib(record['process_family_physical_footprint_bytes']) + ', above the ' +
                           _gib(limits['maximum_process_family_physical_footprint_bytes']) +
                           ' memory limit; RSS is ' + _gib(record['process_family_rss_bytes']) + '.')
        # Check both output and scratch volumes. This is disk space, not the
        # available-RAM startup reserve removed at the user's request.
        disk = []
        try:
            for path in (self.log_path.parent, Path(tempfile.gettempdir())):
                while not path.exists():path = path.parent
                free = shutil.disk_usage(path).free
                disk.append(dict(path=str(path), free_bytes=free))
            record['disk_space'] = disk
            record['minimum_free_disk_bytes'] = MINIMUM_FREE_DISK_BYTES
            if not message and min(x['free_bytes'] for x in disk) < MINIMUM_FREE_DISK_BYTES:
                record['status'] = 'STOP_FREE_DISK_LIMIT'
                error_type = DiskSpaceLimitExceeded
                message = ('Design stopped incomplete: available disk space is below the ' +
                           _gib(MINIMUM_FREE_DISK_BYTES) + ' disk reserve. Free disk space before resuming.')
        except OSError as exc:
            if not message:
                record['status'] = 'DISK_METRICS_UNAVAILABLE'
                error_type = DiskSpaceLimitExceeded
                message = 'Design stopped incomplete: available disk space could not be checked: ' + str(exc)
        if message:
            record['reason'] = message
            if self.logger is not None:
                try:
                    self.logger(message)
                except Exception as exc:
                    # A disconnected supervisor may close worker stdout. The
                    # original stop and its sample must remain authoritative.
                    record['logger_error'] = type(exc).__name__ + ': ' + str(exc)
        self.last_sample = record
        self._last_checked_at = time.monotonic()
        self.log_path.parent.mkdir(parents=True, exist_ok=True)
        with self.log_path.open('a', encoding='utf-8') as handle:
            handle.write(json.dumps(record, sort_keys=True) + '\n')
        if message:
            raise error_type(message, record)
        return record

    def preflight(self, stage='preflight'):
        return self.check(stage, starting=True)

    def poll(self, stage='working'):
        return self.check(stage)

    def checkpoint(self, stage='working', force=False):
        """Cheap cooperative checkpoint; the independent supervisor still polls."""
        now = time.monotonic()
        if (not force and self._last_checked_at is not None
                and now - self._last_checked_at < POLL_INTERVAL_SECONDS):
            return self.last_sample
        return self.check(stage)
