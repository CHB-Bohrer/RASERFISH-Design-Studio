"""Memory-pressure fixtures never allocate large memory or run a full design."""
import json
import errno
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import unittest
from unittest.mock import patch
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import runtime_safety as safety

BASE = Path(__file__).resolve().parents[1]
VM_STAT = '''Mach Virtual Memory Statistics: (page size of 16384 bytes)
Pages free: 100.
Pages inactive: 500.
Pages speculative: 200.
Pages purgeable: 75.
File-backed pages: 300.
Pages occupied by compressor: 250.
'''


class LockTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.path = Path(self.tmp.name) / 'test.lock'

    def tearDown(self):
        self.tmp.cleanup()

    def test_shared_path_independent_of_app_and_port(self):
        self.assertEqual(safety.default_lock_path(), Path('/tmp') / f'raserfish-design-{os.getuid()}.lock')

    def test_second_acquisition_blocked_with_metadata(self):
        with safety.acquire_design_lock({'job_id': 'first'}, self.path) as lock:
            self.assertEqual(lock.metadata['job_id'], 'first')
            with self.assertRaises(safety.DesignAlreadyRunning) as ctx:
                safety.acquire_design_lock(path=self.path)
            self.assertEqual(ctx.exception.metadata['job_id'], 'first')
            inode = self.path.stat().st_ino
        with safety.acquire_design_lock({'job_id': 'second'}, self.path) as lock:
            self.assertEqual(lock.metadata['job_id'], 'second')
            self.assertEqual(self.path.stat().st_ino, inode)
        self.assertTrue(self.path.exists(), 'The shared lock inode must never be unlinked.')

    def test_release_is_idempotent(self):
        lock = safety.acquire_design_lock(path=self.path)
        lock.release()
        lock.release()
        with self.assertRaises(RuntimeError):
            lock.update_metadata({'x': 1})

    def test_separate_process_cannot_acquire_held_lock(self):
        with safety.acquire_design_lock(path=self.path):
            script = ('import sys; from runtime_safety import *\n'
                      'try: acquire_design_lock(path=sys.argv[1])\n'
                      'except DesignAlreadyRunning: sys.exit(23)\n')
            result = subprocess.run([sys.executable, '-c', script, str(self.path)], cwd=BASE, timeout=10)
            self.assertEqual(result.returncode, 23)

    def test_inherited_worker_holds_lock_after_supervisor_exit(self):
        worker_script = (
            'import sys; from runtime_safety import *\n'
            'lock=inherited_design_lock(int(sys.argv[2]),sys.argv[1])\n'
            'print("ready",flush=True)\n'
            'sys.stdin.readline()\n'
            'lock.release()\n')
        lock = safety.acquire_design_lock(path=self.path)
        child = subprocess.Popen([sys.executable, '-c', worker_script, str(self.path), str(lock.fd)],
                                 cwd=BASE, pass_fds=(lock.fd,), stdin=subprocess.PIPE,
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        try:
            self.assertEqual(child.stdout.readline().strip(), 'ready')
            lock.release()  # Same descriptor closure as supervisor process termination.
            with self.assertRaises(safety.DesignAlreadyRunning):
                safety.acquire_design_lock(path=self.path)
            child.communicate('\n', timeout=10)
            self.assertEqual(child.returncode, 0)
            with safety.acquire_design_lock(path=self.path):
                pass
        finally:
            lock.release()
            if child.poll() is None:
                child.kill()
                child.wait(timeout=10)
            for handle in (child.stdin, child.stdout, child.stderr):
                if handle:
                    handle.close()

    def test_interrupted_process_releases_lock(self):
        script = ('import sys,time; from runtime_safety import *\n'
                  'lock=acquire_design_lock(path=sys.argv[1]); print("ready",flush=True); time.sleep(30)')
        child = subprocess.Popen([sys.executable, '-c', script, str(self.path)], cwd=BASE,
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        try:
            self.assertEqual(child.stdout.readline().strip(), 'ready')
            child.terminate()
            child.wait(timeout=10)
            with safety.acquire_design_lock(path=self.path):
                pass
        finally:
            if child.poll() is None:
                child.kill()
                child.wait(timeout=10)
            child.stdout.close()
            child.stderr.close()

    def test_unsafe_symlink_and_permissions_rejected(self):
        target = Path(self.tmp.name) / 'other'
        target.write_text('unchanged')
        self.path.symlink_to(target)
        with self.assertRaises(OSError):
            safety.acquire_design_lock(path=self.path)
        self.assertEqual(target.read_text(), 'unchanged')
        self.path.unlink()
        self.path.touch(mode=0o644)
        with self.assertRaises(RuntimeError):
            safety.acquire_design_lock(path=self.path)

    def test_inherited_wrong_file_rejected(self):
        with safety.acquire_design_lock(path=self.path):
            other = Path(self.tmp.name) / 'other'
            fd = os.open(other, os.O_CREAT | os.O_RDWR, 0o600)
            try:
                with self.assertRaises(RuntimeError):
                    safety.inherited_design_lock(fd, self.path)
            finally:
                os.close(fd)


class MemoryTests(unittest.TestCase):
    def test_mac_page_size_and_no_double_counting(self):
        result = safety.parse_vm_stat(VM_STAT)
        self.assertEqual(result['available_bytes'], (100 + 300) * 16384)
        self.assertEqual(result['compressed_bytes'], 250 * 16384)
        result = safety.parse_vm_stat(VM_STAT.replace('500.', '100.'))
        self.assertEqual(result['available_bytes'], (100 + 100) * 16384)

    def test_missing_mac_metrics_rejected(self):
        with self.assertRaises(ValueError):
            safety.parse_vm_stat('Pages free: 100.')
        with self.assertRaises(ValueError):
            safety.parse_vm_stat(VM_STAT.replace('File-backed pages: 300.', ''))

    def test_mac_snapshot_has_no_external_hw_memsize_command(self):
        calls = []
        def run(args):
            calls.append(args)
            return VM_STAT
        with patch('runtime_safety.platform.system', return_value='Darwin'), patch('runtime_safety._command', side_effect=run), patch('runtime_safety._mac_physical_memory_bytes', return_value=16 * safety.GIB):
            result = safety.system_memory_snapshot()
        self.assertEqual(result['physical_bytes'], 16 * safety.GIB)
        self.assertEqual(calls, [['/usr/bin/vm_stat']])

    def test_process_tree_includes_native_children_excludes_unrelated(self):
        table = '11 1 100 /python\n12 11 200 /blastn\n13 12 300 /child\n14 1 999 /other\n'
        with patch('runtime_safety._command', return_value=table), patch('runtime_safety.platform.system', return_value='Linux'):
            result = safety.process_tree_snapshot(11)
        self.assertEqual([r['pid'] for r in result['processes']], [11, 12, 13])
        self.assertEqual(result['process_family_rss_bytes'], 600 * 1024)
        with patch('runtime_safety._command', return_value=table):
            with self.assertRaises(RuntimeError):
                safety.process_tree_snapshot(404)

    def test_memory_limits_for_this_16gb_mac(self):
        limits = safety.memory_limits(16 * safety.GIB)
        self.assertEqual(limits['available_memory_enforcement'], 'disabled_by_user')
        self.assertNotIn('minimum_start_available_bytes', limits)
        self.assertNotIn('minimum_available_bytes', limits)
        self.assertEqual(limits['maximum_process_family_rss_bytes'], 8 * safety.GIB)
        self.assertEqual(limits['maximum_process_family_physical_footprint_bytes'], 8 * safety.GIB)
        self.assertEqual(limits['poll_interval_seconds'], 5)

    def fixture_guard(self, folder, available=6 * safety.GIB, rss=1 * safety.GIB, footprint=None):
        footprint = rss if footprint is None else footprint
        system = dict(physical_bytes=16 * safety.GIB, available_bytes=available, available_method='test fixture')
        process = dict(root_pid=11, process_count=2, process_family_rss_bytes=rss,
                       process_family_physical_footprint_bytes=footprint,
                       processes=[dict(pid=11, ppid=1, rss_bytes=rss // 2, name='test worker'),
                                  dict(pid=12, ppid=11, rss_bytes=rss // 2, name='test blast')])
        return safety.MemoryGuard(Path(folder) / 'memory.jsonl', root_pid=11), system, process

    def test_memory_log_contains_metrics_processes_limits_stage(self):
        with tempfile.TemporaryDirectory() as folder:
            guard, system, process = self.fixture_guard(folder)
            with patch('runtime_safety.system_memory_snapshot', return_value=system), patch('runtime_safety.process_tree_snapshot', return_value=process):
                guard.preflight()
                guard.poll('genomic BLAST')
            records = [json.loads(line) for line in guard.log_path.read_text().splitlines()]
        self.assertEqual(len(records), 2)
        self.assertEqual(records[0]['starting'], True)
        self.assertEqual(records[1]['stage'], 'genomic BLAST')
        self.assertEqual(records[1]['status'], 'OK')
        self.assertEqual(records[1]['process_count'], 2)
        self.assertIn('maximum_process_family_rss_bytes', records[1]['limits'])

    def test_reported_216_gib_allows_start_and_running_checks(self):
        with tempfile.TemporaryDirectory() as folder:
            guard, system, process = self.fixture_guard(folder, available=int(2.16 * safety.GIB))
            with patch('runtime_safety.system_memory_snapshot', return_value=system), patch('runtime_safety.process_tree_snapshot', return_value=process):
                self.assertEqual(guard.preflight()['status'], 'OK')
                self.assertEqual(guard.poll()['status'], 'OK')
            records = [json.loads(line) for line in guard.log_path.read_text().splitlines()]
        self.assertEqual(len(records), 2)
        for record in records:
            self.assertEqual(record['available_bytes'], int(2.16 * safety.GIB))
            self.assertEqual(record['limits']['available_memory_enforcement'], 'disabled_by_user')
            self.assertNotIn('reason', record)

    def test_zero_available_memory_is_logged_without_reserve_stop(self):
        with tempfile.TemporaryDirectory() as folder:
            guard, system, process = self.fixture_guard(folder, available=0)
            with patch('runtime_safety.system_memory_snapshot', return_value=system), patch('runtime_safety.process_tree_snapshot', return_value=process):
                guard.poll('individual folding')
            record = json.loads(guard.log_path.read_text())
        self.assertEqual(record['status'], 'OK')
        self.assertEqual(record['available_bytes'], 0)
        self.assertEqual(record['stage'], 'individual folding')

    def test_family_cap_stops_even_with_plenty_of_available_memory(self):
        with tempfile.TemporaryDirectory() as folder:
            guard, system, process = self.fixture_guard(folder, available=10 * safety.GIB, rss=9 * safety.GIB)
            with patch('runtime_safety.system_memory_snapshot', return_value=system), patch('runtime_safety.process_tree_snapshot', return_value=process):
                with self.assertRaises(safety.MemoryLimitExceeded) as ctx:
                    guard.poll()
        self.assertEqual(ctx.exception.sample['status'], 'STOP_PROCESS_MEMORY_LIMIT')

    def test_low_rss_large_compressed_footprint_stops_and_preserves_sample(self):
        with tempfile.TemporaryDirectory() as folder:
            guard, system, process = self.fixture_guard(folder, rss=6 * safety.GIB, footprint=41_010_000_000)
            with patch('runtime_safety.system_memory_snapshot', return_value=system), patch('runtime_safety.process_tree_snapshot', return_value=process):
                with self.assertRaises(safety.MemoryLimitExceeded) as ctx:
                    guard.poll('genomic BLAST')
            logged = json.loads(guard.log_path.read_text())
        self.assertEqual(ctx.exception.sample['status'], 'STOP_PROCESS_PHYSICAL_FOOTPRINT_LIMIT')
        self.assertEqual(logged, ctx.exception.sample)
        self.assertIn('compressed/swapped', str(ctx.exception))
        self.assertEqual(logged['process_family_rss_bytes'], 6 * safety.GIB)
        self.assertEqual(logged['process_family_physical_footprint_bytes'], 41_010_000_000)

    def test_rss_cap_independent_of_lower_footprint(self):
        with tempfile.TemporaryDirectory() as folder:
            guard, system, process = self.fixture_guard(folder, rss=9 * safety.GIB, footprint=1 * safety.GIB)
            with patch('runtime_safety.system_memory_snapshot', return_value=system), patch('runtime_safety.process_tree_snapshot', return_value=process):
                with self.assertRaises(safety.MemoryLimitExceeded) as ctx:
                    guard.poll()
        self.assertEqual(ctx.exception.sample['status'], 'STOP_PROCESS_MEMORY_LIMIT')

    def test_mac_missing_footprint_fails_closed(self):
        with tempfile.TemporaryDirectory() as folder:
            guard, system, process = self.fixture_guard(folder)
            process.pop('process_family_physical_footprint_bytes')
            with patch('runtime_safety.platform.system', return_value='Darwin'), patch('runtime_safety.system_memory_snapshot', return_value=system), patch('runtime_safety.process_tree_snapshot', return_value=process):
                with self.assertRaises(safety.MemoryMetricsUnavailable) as ctx:
                    guard.poll()
        self.assertIn('physical-footprint metric is missing', str(ctx.exception))

    def test_failed_logger_does_not_replace_memory_stop(self):
        with tempfile.TemporaryDirectory() as folder:
            guard, system, process = self.fixture_guard(folder, footprint=41 * safety.GIB)
            def broken_logger(message):
                raise BrokenPipeError('supervisor disconnected')
            guard.logger = broken_logger
            with patch('runtime_safety.system_memory_snapshot', return_value=system), patch('runtime_safety.process_tree_snapshot', return_value=process):
                with self.assertRaises(safety.MemoryLimitExceeded) as ctx:
                    guard.poll()
            record = json.loads(guard.log_path.read_text())
        self.assertEqual(record, ctx.exception.sample)
        self.assertEqual(record['logger_error'], 'BrokenPipeError: supervisor disconnected')
        self.assertEqual(record['status'], 'STOP_PROCESS_PHYSICAL_FOOTPRINT_LIMIT')

    def test_mac_tree_sums_only_owned_footprints_and_tolerates_exited_child(self):
        table = '11 1 100 /python\n12 11 200 /blastn\n13 11 300 /ps\n14 1 999 /other\n'
        def memory(pid):
            if pid == 13:
                raise OSError(errno.ESRCH, 'exited')
            return dict(physical_footprint_bytes=pid * 1000, lifetime_peak_physical_footprint_bytes=pid * 2000)
        with patch('runtime_safety._command', return_value=table), patch('runtime_safety.platform.system', return_value='Darwin'), patch('runtime_safety.mac_process_memory', side_effect=memory) as get_memory:
            result = safety.process_tree_snapshot(11)
        self.assertEqual([r['pid'] for r in result['processes']], [11, 12])
        self.assertEqual(result['process_family_physical_footprint_bytes'], 23000)
        self.assertEqual(result['process_family_rss_bytes'], 300 * 1024)
        self.assertEqual(result['processes_exited_during_snapshot'], [13])
        self.assertEqual([c.args[0] for c in get_memory.call_args_list], [11, 12, 13])

    def test_mac_root_exit_and_child_permission_failure_are_not_ignored(self):
        table = '11 1 100 /python\n12 11 200 /blastn\n'
        with patch('runtime_safety._command', return_value=table), patch('runtime_safety.platform.system', return_value='Darwin'):
            with patch('runtime_safety.mac_process_memory', side_effect=OSError(errno.ESRCH, 'root gone')):
                with self.assertRaises(OSError):
                    safety.process_tree_snapshot(11)
            with patch('runtime_safety.mac_process_memory', side_effect=[{'physical_footprint_bytes':100}] + [OSError(errno.EPERM, 'denied')] * 4), patch('runtime_safety._process_lifecycle_snapshot', return_value=dict(pid=12, ppid=11, name='/blastn', state='live', process_status='R')), patch('runtime_safety.time.sleep'):
                with self.assertRaises(OSError):
                    safety.process_tree_snapshot(11)

    def test_mac_libproc_success_and_failure_preserve_errno(self):
        def fake_call(pid, flavor, buffer):
            self.assertEqual((pid, flavor), (123, 4))
            buffer._obj.ri_phys_footprint = 987654321
            buffer._obj.ri_resident_size = 1234
            buffer._obj.ri_lifetime_max_phys_footprint = 999999999
            return 0
        with patch('runtime_safety._mac_proc_pid_rusage', return_value=fake_call):
            sample = safety.mac_process_memory(123)
        self.assertEqual(sample['physical_footprint_bytes'], 987654321)
        self.assertEqual(sample['kernel_resident_bytes'], 1234)
        self.assertEqual(sample['lifetime_peak_physical_footprint_bytes'], 999999999)
        def missing(*args):
            safety.ctypes.set_errno(errno.ESRCH)
            return -1
        with patch('runtime_safety._mac_proc_pid_rusage', return_value=missing):
            with self.assertRaises(OSError) as ctx:
                safety.mac_process_memory(123)
        self.assertEqual(ctx.exception.errno, errno.ESRCH)

    @unittest.skipUnless(sys.platform == 'darwin', 'macOS read-only libproc smoke check')
    def test_real_mac_current_process_footprint(self):
        sample = safety.mac_process_memory(os.getpid())
        self.assertGreater(sample['physical_footprint_bytes'], 0)
        self.assertGreater(sample['lifetime_peak_physical_footprint_bytes'], 0)
        self.assertGreater(sample['kernel_resident_bytes'], 0)

    def test_metrics_unavailable_fails_closed_and_logs_error(self):
        with tempfile.TemporaryDirectory() as folder:
            guard, system, process = self.fixture_guard(folder)
            with patch('runtime_safety.system_memory_snapshot', side_effect=OSError('denied')), patch('runtime_safety.process_tree_snapshot', return_value=process):
                with self.assertRaises(safety.MemoryMetricsUnavailable):
                    guard.preflight()
            record = json.loads(guard.log_path.read_text())
        self.assertEqual(record['status'], 'METRICS_UNAVAILABLE')
        self.assertFalse(record['metrics_available'])

    def test_missing_worker_rss_fails_closed(self):
        with tempfile.TemporaryDirectory() as folder:
            guard, system, process = self.fixture_guard(folder)
            with patch('runtime_safety.system_memory_snapshot', return_value=system), patch('runtime_safety.process_tree_snapshot', side_effect=RuntimeError('no worker')):
                with self.assertRaises(safety.MemoryMetricsUnavailable):
                    guard.poll()

    def test_cooperative_checkpoint_throttles_without_disabling_supervisor(self):
        with tempfile.TemporaryDirectory() as folder:
            guard, system, process = self.fixture_guard(folder)
            with patch('runtime_safety.system_memory_snapshot', return_value=system) as get_system, patch('runtime_safety.process_tree_snapshot', return_value=process):
                guard.checkpoint('before folding')
                guard.checkpoint('before next folding')
                self.assertEqual(get_system.call_count, 1)
                guard.checkpoint('new stage', force=True)
                self.assertEqual(get_system.call_count, 2)
                guard.poll('supervisor poll')
                self.assertEqual(get_system.call_count, 3)


class NativePhysicalMemoryTests(unittest.TestCase):
    """Static RAM may be cached; every changing safety metric must stay fresh."""

    def setUp(self):
        safety._mac_physical_memory_bytes.cache_clear()

    def tearDown(self):
        safety._mac_physical_memory_bytes.cache_clear()

    def fake_native(self, value=16 * safety.GIB, length=8, code=0):
        def read(name, buffer, size, new_value, new_size):
            self.assertEqual(name, b'hw.memsize')
            self.assertIsNone(new_value, 'The sysctl call must be read-only.')
            self.assertEqual(new_size, 0)
            self.assertEqual(size._obj.value, 8)
            buffer._obj.value = value
            size._obj.value = length
            safety.ctypes.set_errno(code)
            return -1 if code else 0
        return read

    def test_native_read_is_readonly_validated_and_cached_without_subprocess(self):
        native = self.fake_native()
        with patch('runtime_safety.ctypes.CDLL', return_value=SimpleNamespace(sysctlbyname=native)) as load, patch('runtime_safety.subprocess.run', side_effect=AssertionError('No external process permitted')):
            self.assertEqual(safety._mac_physical_memory_bytes(), 16 * safety.GIB)
            self.assertEqual(safety._mac_physical_memory_bytes(), 16 * safety.GIB)
        load.assert_called_once_with('/usr/lib/libSystem.B.dylib', use_errno=True)
        self.assertEqual(native.restype, safety.ctypes.c_int)
        self.assertEqual(native.argtypes, (safety.ctypes.c_char_p, safety.ctypes.c_void_p, safety.ctypes.POINTER(safety.ctypes.c_size_t), safety.ctypes.c_void_p, safety.ctypes.c_size_t))

    def test_native_error_preserves_errno_and_is_not_cached(self):
        with patch('runtime_safety.ctypes.CDLL', side_effect=[SimpleNamespace(sysctlbyname=self.fake_native(code=errno.EPERM)), SimpleNamespace(sysctlbyname=self.fake_native())]) as load:
            with self.assertRaises(OSError) as failure:
                safety._mac_physical_memory_bytes()
            self.assertEqual(failure.exception.errno, errno.EPERM)
            self.assertIn('hw.memsize', str(failure.exception))
            self.assertEqual(safety._mac_physical_memory_bytes(), 16 * safety.GIB)
            self.assertEqual(load.call_count, 2)

    def test_invalid_size_or_zero_ram_is_rejected_and_not_cached(self):
        for value, length in ((0, 8), (16 * safety.GIB, 4), (16 * safety.GIB, 0), (16 * safety.GIB, 16)):
            with self.subTest(value=value, length=length):
                safety._mac_physical_memory_bytes.cache_clear()
                with patch('runtime_safety.ctypes.CDLL', side_effect=[SimpleNamespace(sysctlbyname=self.fake_native(value, length)), SimpleNamespace(sysctlbyname=self.fake_native())]):
                    with self.assertRaises(ValueError):
                        safety._mac_physical_memory_bytes()
                    self.assertEqual(safety._mac_physical_memory_bytes(), 16 * safety.GIB)

    def test_missing_native_function_or_library_fails_closed(self):
        for result in (SimpleNamespace(), OSError(errno.ENOENT, 'library unavailable')):
            with self.subTest(result=result):
                mocked = {'side_effect': result} if isinstance(result, Exception) else {'return_value': result}
                with patch('runtime_safety.ctypes.CDLL', **mocked):
                    with self.assertRaises((OSError, RuntimeError)):
                        safety._mac_physical_memory_bytes()

    def test_system_snapshots_refresh_vm_stat_after_static_ram_cached(self):
        with patch('runtime_safety.ctypes.CDLL', return_value=SimpleNamespace(sysctlbyname=self.fake_native())) as load, patch('runtime_safety.platform.system', return_value='Darwin'), patch('runtime_safety._command', side_effect=[VM_STAT, VM_STAT.replace('Pages free: 100.', 'Pages free: 50.')]) as command:
            first = safety.system_memory_snapshot()
            second = safety.system_memory_snapshot()
        self.assertEqual(load.call_count, 1)
        self.assertEqual([c.args[0] for c in command.call_args_list], [['/usr/bin/vm_stat'], ['/usr/bin/vm_stat']])
        self.assertEqual(first['available_bytes'] - second['available_bytes'], 50 * 16384)

    def test_vm_stat_timeout_after_cache_stops_incomplete_and_records_error(self):
        process = dict(process_family_rss_bytes=safety.GIB, process_family_physical_footprint_bytes=safety.GIB)
        timeout = subprocess.TimeoutExpired(['/usr/bin/vm_stat'], 3)
        with tempfile.TemporaryDirectory() as folder, patch('runtime_safety.ctypes.CDLL', return_value=SimpleNamespace(sysctlbyname=self.fake_native())) as load, patch('runtime_safety.platform.system', return_value='Darwin'), patch('runtime_safety._command', side_effect=[VM_STAT, timeout]), patch('runtime_safety.process_tree_snapshot', return_value=process) as processes:
            guard = safety.MemoryGuard(Path(folder) / 'memory.jsonl')
            self.assertEqual(guard.check()['status'], 'OK')
            with self.assertRaises(safety.MemoryMetricsUnavailable) as failure:
                guard.check('RNA BLAST')
            records = [json.loads(line) for line in guard.log_path.read_text().splitlines()]
            self.assertEqual(records[-1], failure.exception.sample)
        self.assertEqual(load.call_count, 1)
        self.assertEqual(processes.call_count, 2)
        self.assertEqual(records[-1]['status'], 'METRICS_UNAVAILABLE')
        self.assertIn('vm_stat', records[-1]['reason'])
        self.assertIn('stopped incomplete', records[-1]['reason'])

    def test_cached_ram_does_not_hide_new_process_memory_or_disk_limit(self):
        process = dict(process_family_rss_bytes=safety.GIB, process_family_physical_footprint_bytes=safety.GIB)
        for failure_kind in ('footprint', 'process_metrics', 'disk'):
            with self.subTest(failure_kind=failure_kind):
                safety._mac_physical_memory_bytes.cache_clear()
                second_process = dict(process, process_family_physical_footprint_bytes=9 * safety.GIB) if failure_kind == 'footprint' else (RuntimeError('unavailable live worker') if failure_kind == 'process_metrics' else process)
                second_free = safety.GIB // 2 if failure_kind == 'disk' else 2 * safety.GIB
                disk = [SimpleNamespace(free=free) for free in (2 * safety.GIB, 2 * safety.GIB, second_free, second_free)]
                with tempfile.TemporaryDirectory() as folder, patch('runtime_safety.ctypes.CDLL', return_value=SimpleNamespace(sysctlbyname=self.fake_native())) as load, patch('runtime_safety.platform.system', return_value='Darwin'), patch('runtime_safety._command', return_value=VM_STAT) as command, patch('runtime_safety.process_tree_snapshot', side_effect=[process, second_process]) as processes, patch('runtime_safety.shutil.disk_usage', side_effect=disk) as disks:
                    guard = safety.MemoryGuard(Path(folder) / 'memory.jsonl')
                    self.assertEqual(guard.check()['status'], 'OK')
                    with self.assertRaises(safety.MemoryLimitExceeded) as failure:
                        guard.check()
                self.assertEqual(load.call_count, 1)
                self.assertEqual(command.call_count, 2)
                self.assertEqual(processes.call_count, 2)
                self.assertEqual(disks.call_count, 4)
                self.assertEqual(failure.exception.sample['status'], {'footprint':'STOP_PROCESS_PHYSICAL_FOOTPRINT_LIMIT', 'process_metrics':'METRICS_UNAVAILABLE', 'disk':'STOP_FREE_DISK_LIMIT'}[failure_kind])

    @unittest.skipUnless(sys.platform == 'darwin', 'macOS read-only sysctl smoke check')
    def test_real_mac_native_physical_memory_read_without_subprocess(self):
        with patch('runtime_safety.subprocess.run', side_effect=AssertionError('No external process permitted')):
            physical = safety._mac_physical_memory_bytes()
        self.assertIsInstance(physical, int)
        self.assertGreater(physical, 0)
        self.assertEqual(safety._mac_physical_memory_bytes(), physical)


class ChildMemoryRetryTests(unittest.TestCase):
    """A denied monitor helper can exit; a live unmeasurable worker cannot pass."""

    table = '11 1 100 /python\n12 11 200 /bin/ps\n14 1 999 /unrelated\n'
    root_memory = dict(physical_footprint_bytes=1000, lifetime_peak_physical_footprint_bytes=2000)
    live = dict(pid=12, ppid=11, name='/bin/ps', state='live', process_status='R')
    row = dict(pid=12, ppid=11, name='/bin/ps', rss_bytes=200 * 1024)

    def tree(self, memory, lifecycle):
        with patch('runtime_safety._command', return_value=self.table), patch('runtime_safety.platform.system', return_value='Darwin'), patch('runtime_safety.mac_process_memory', side_effect=memory) as metrics, patch('runtime_safety._process_lifecycle_snapshot', side_effect=lifecycle) as status, patch('runtime_safety.time.sleep') as pause:
            result = safety.process_tree_snapshot(11)
        return result, metrics, status, pause

    def test_denied_child_verified_absent_or_zombie_is_audited_and_excluded(self):
        for code in (errno.EPERM, errno.EACCES):
            for state in (dict(pid=12, state='absent'), dict(self.live, state='zombie', process_status='Z')):
                with self.subTest(code=code, state=state['state']):
                    result, metrics, status, pause = self.tree([self.root_memory, OSError(code, 'denied')], [state])
                    self.assertEqual(result['process_count'], 1)
                    self.assertEqual(result['process_family_physical_footprint_bytes'], 1000)
                    self.assertEqual(result['process_family_rss_bytes'], 100 * 1024)
                    self.assertEqual(result['processes_exited_during_snapshot'], [12])
                    event = result['process_snapshot_events'][0]
                    self.assertEqual(event['outcome'], 'verified_exited_after_permission_error')
                    self.assertEqual(event['initial_errno'], code)
                    self.assertEqual(event['lifecycle_checks'], [state])
                    self.assertEqual([call.args[0] for call in metrics.call_args_list], [11, 12])
                    status.assert_called_once_with(12)
                    pause.assert_not_called()

    def test_live_child_recovers_full_memory_accounting(self):
        recovered = dict(physical_footprint_bytes=123456, lifetime_peak_physical_footprint_bytes=200000, kernel_resident_bytes=9000, process_start_abstime=123)
        result, metrics, status, pause = self.tree([self.root_memory, OSError(errno.EPERM, 'denied'), recovered], [self.live])
        self.assertEqual(result['process_family_physical_footprint_bytes'], 124456)
        self.assertEqual(result['process_family_rss_bytes'], 300 * 1024)
        self.assertEqual(result['processes'][1]['kernel_resident_bytes'], 9000)
        self.assertEqual(result['process_snapshot_events'][0]['outcome'], 'memory_read_recovered')
        self.assertEqual(result['processes_exited_during_snapshot'], [])
        self.assertEqual([call.args[0] for call in metrics.call_args_list], [11, 12, 12])
        pause.assert_not_called()

    def test_live_then_gone_child_is_verified_before_exclusion(self):
        gone = dict(pid=12, state='absent')
        result, metrics, status, pause = self.tree([self.root_memory, OSError(errno.EPERM, 'denied'), OSError(errno.EPERM, 'still denied')], [self.live, gone])
        event = result['process_snapshot_events'][0]
        self.assertEqual(event['lifecycle_checks'], [self.live, gone])
        self.assertEqual(event['outcome'], 'verified_exited_after_permission_error')
        self.assertEqual(status.call_count, 2)
        self.assertEqual(metrics.call_count, 3)
        pause.assert_called_once_with(.02)

    def test_retry_esrch_is_explicit_exit_evidence(self):
        result, metrics, status, pause = self.tree([self.root_memory, OSError(errno.EPERM, 'denied'), OSError(errno.ESRCH, 'gone')], [self.live])
        self.assertEqual(result['processes_exited_during_snapshot'], [12])
        self.assertEqual(result['process_snapshot_events'][0]['outcome'], 'verified_exited_during_retry')
        pause.assert_not_called()

    def test_persistently_live_denied_child_stops_after_three_checks(self):
        for code in (errno.EPERM, errno.EACCES):
            with self.subTest(code=code), patch('runtime_safety._process_lifecycle_snapshot', return_value=self.live) as status, patch('runtime_safety.mac_process_memory', side_effect=OSError(code, 'denied')) as metrics, patch('runtime_safety.time.sleep') as pause:
                with self.assertRaises(OSError) as failure:
                    safety._recover_child_memory_sample(self.row, OSError(code, 'initial denial'))
                self.assertEqual(failure.exception.errno, code)
                self.assertIn('PID 12', str(failure.exception))
                self.assertIn('/bin/ps', str(failure.exception))
                self.assertIn('parent 11', str(failure.exception))
                self.assertEqual(status.call_count, 3)
                self.assertEqual(metrics.call_count, 3)
                self.assertEqual(pause.call_count, 2)
                event = failure.exception.process_snapshot_error
                self.assertEqual(event['outcome'], 'live_process_memory_denied')
                self.assertEqual(event['lifecycle_checks'], [self.live] * 3)

    def test_denied_root_never_uses_child_lifecycle_exception(self):
        for code in (errno.EPERM, errno.EACCES, errno.ESRCH):
            with self.subTest(code=code), patch('runtime_safety._command', return_value=self.table), patch('runtime_safety.platform.system', return_value='Darwin'), patch('runtime_safety.mac_process_memory', side_effect=OSError(code, 'root unavailable')), patch('runtime_safety._process_lifecycle_snapshot') as status:
                with self.assertRaises(OSError) as failure:
                    safety.process_tree_snapshot(11)
                self.assertEqual(failure.exception.errno, code)
                status.assert_not_called()

    def test_failed_lifecycle_lookup_does_not_authorize_exclusion(self):
        for error in (RuntimeError('invalid process table'), PermissionError(errno.EPERM, 'ps blocked'), subprocess.TimeoutExpired('/bin/ps', 3)):
            with self.subTest(error=type(error).__name__), patch('runtime_safety._process_lifecycle_snapshot', side_effect=error), patch('runtime_safety.mac_process_memory') as metrics:
                with self.assertRaises(type(error)) as failure:
                    safety._recover_child_memory_sample(self.row, OSError(errno.EPERM, 'denied'))
                metrics.assert_not_called()
                self.assertEqual(failure.exception.process_snapshot_error['outcome'], 'lifecycle_verification_failed')

    def test_nonpermission_retry_error_is_not_hidden(self):
        with patch('runtime_safety._process_lifecycle_snapshot', return_value=self.live), patch('runtime_safety.mac_process_memory', side_effect=OSError(errno.EIO, 'accounting failed')):
            with self.assertRaises(OSError) as failure:
                safety._recover_child_memory_sample(self.row, OSError(errno.EPERM, 'denied'))
        self.assertEqual(failure.exception.errno, errno.EIO)
        self.assertEqual(failure.exception.process_snapshot_error['outcome'], 'memory_retry_failed')

    def test_changed_parent_is_not_treated_as_original_descendant(self):
        for process_state in ('live', 'zombie'):
            changed = dict(self.live, ppid=999, state=process_state)
            with self.subTest(state=process_state), patch('runtime_safety._process_lifecycle_snapshot', return_value=changed), patch('runtime_safety.mac_process_memory') as metrics:
                with self.assertRaises(RuntimeError) as failure:
                    safety._recover_child_memory_sample(self.row, OSError(errno.EPERM, 'denied'))
                metrics.assert_not_called()
                self.assertEqual(failure.exception.process_snapshot_error['outcome'], 'process_identity_changed')
                self.assertEqual(failure.exception.process_snapshot_error['lifecycle_checks'], [changed])

    def test_persistent_live_denial_is_recorded_in_guard_stop_sample(self):
        system = dict(physical_bytes=16 * safety.GIB, available_bytes=6 * safety.GIB)
        memory = [self.root_memory] + [OSError(errno.EPERM, 'denied')] * 4
        with tempfile.TemporaryDirectory() as folder, patch('runtime_safety.system_memory_snapshot', return_value=system), patch('runtime_safety._command', return_value=self.table), patch('runtime_safety.platform.system', return_value='Darwin'), patch('runtime_safety.mac_process_memory', side_effect=memory), patch('runtime_safety._process_lifecycle_snapshot', return_value=self.live), patch('runtime_safety.time.sleep'):
            guard = safety.MemoryGuard(Path(folder) / 'memory.jsonl', root_pid=11)
            with self.assertRaises(safety.MemoryMetricsUnavailable) as failure:
                guard.check('live denied worker')
            sample = failure.exception.sample
            self.assertEqual(sample['status'], 'METRICS_UNAVAILABLE')
            self.assertEqual(sample['process_snapshot_error']['pid'], 12)
            self.assertEqual(sample['process_snapshot_error']['ppid'], 11)
            self.assertEqual(sample['process_snapshot_error']['name'], '/bin/ps')
            self.assertEqual(sample['process_snapshot_error']['outcome'], 'live_process_memory_denied')
            self.assertEqual(sample['process_snapshot_error']['lifecycle_checks'], [self.live] * 3)
            self.assertEqual(json.loads(guard.log_path.read_text()), sample)

    def test_lifecycle_query_accepts_only_clean_missing_pid_return(self):
        done = subprocess.CompletedProcess([], 1, stdout='', stderr='')
        with patch('runtime_safety.subprocess.run', return_value=done) as run:
            self.assertEqual(safety._process_lifecycle_snapshot(12), dict(pid=12, state='absent'))
        self.assertEqual(run.call_args.args[0], ['/bin/ps', '-p', '12', '-o', 'pid=,ppid=,stat=,comm='])
        for done in (subprocess.CompletedProcess([], 1, stdout='', stderr='permission denied'), subprocess.CompletedProcess([], 2, stdout='', stderr='usage error'), subprocess.CompletedProcess([], -15, stdout='', stderr=''), subprocess.CompletedProcess([], 1, stdout='12 11 R /bin/ps', stderr='')):
            with self.subTest(returncode=done.returncode, stderr=done.stderr), patch('runtime_safety.subprocess.run', return_value=done):
                with self.assertRaises(RuntimeError):
                    safety._process_lifecycle_snapshot(12)

    def test_lifecycle_query_parses_live_and_zombie_status(self):
        for process_status, state in (('R', 'live'), ('S+', 'live'), ('Z', 'zombie'), ('Z+', 'zombie')):
            done = subprocess.CompletedProcess([], 0, stdout=f'  12  11  {process_status} /path with spaces/ps\n', stderr='')
            with self.subTest(process_status=process_status), patch('runtime_safety.subprocess.run', return_value=done):
                self.assertEqual(safety._process_lifecycle_snapshot(12), dict(pid=12, ppid=11, process_status=process_status, name='/path with spaces/ps', state=state))

    def test_lifecycle_query_rejects_invalid_or_mismatched_response(self):
        for output in ('', '12 11', '13 11 R /bin/ps', '12 nope R /bin/ps', '12 11 R', 'garbage'):
            done = subprocess.CompletedProcess([], 0, stdout=output, stderr='')
            with self.subTest(output=output), patch('runtime_safety.subprocess.run', return_value=done):
                with self.assertRaises(RuntimeError):
                    safety._process_lifecycle_snapshot(12)

    def test_recovered_oversized_footprint_still_stops_guard(self):
        system = dict(physical_bytes=16 * safety.GIB, available_bytes=6 * safety.GIB)
        memory = [self.root_memory, OSError(errno.EPERM, 'transient'), dict(physical_footprint_bytes=9 * safety.GIB)]
        with tempfile.TemporaryDirectory() as folder, patch('runtime_safety.system_memory_snapshot', return_value=system), patch('runtime_safety._command', return_value=self.table), patch('runtime_safety.platform.system', return_value='Darwin'), patch('runtime_safety.mac_process_memory', side_effect=memory), patch('runtime_safety._process_lifecycle_snapshot', return_value=self.live):
            guard = safety.MemoryGuard(Path(folder) / 'memory.jsonl', root_pid=11)
            with self.assertRaises(safety.MemoryLimitExceeded) as failure:
                guard.check('recovered helper')
            self.assertEqual(failure.exception.sample['status'], 'STOP_PROCESS_PHYSICAL_FOOTPRINT_LIMIT')
            self.assertEqual(failure.exception.sample['process_snapshot_events'][0]['outcome'], 'memory_read_recovered')


if __name__ == '__main__':
    unittest.main()
