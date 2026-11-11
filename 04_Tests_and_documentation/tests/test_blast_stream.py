"""Real tiny producers exercise streamed output and owned-process cleanup."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import threading
import time
import unittest
from unittest.mock import patch

from blast_stream import stream_command, MAX_BUFFER_BYTES
from process_diagnostics import run_process


class Cancelled(RuntimeError):
    pass


class BlastStreamTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.output = self.root / 'alignments.fifo'
        self.log = self.root / 'blast.log'

    def tearDown(self):
        self.temporary.cleanup()

    @staticmethod
    def runner(command, outlog, cancel, guard):
        return run_process(command, outlog, cancel, guard, cancellation_class=Cancelled)

    def stream(self, source, cancel=None, guard=None, runner=None):
        command = [sys.executable, '-c', source, str(self.output)]
        return stream_command(command, self.log, self.output, cancel, guard,
                              run_process=runner or self.runner, cancellation_class=Cancelled)

    def diagnostic(self):
        return json.loads(Path(str(self.log) + '.process.json').read_text())

    def test_fragmented_success_preserves_every_raw_byte(self):
        source = "import sys,time; f=open(sys.argv[1],'wb',buffering=0); f.write(b'q1\\tA'); time.sleep(.02); f.write(b'CG\\nq2\\tT\\n'); f.close()"
        with self.stream(source) as lines:
            self.assertEqual(list(lines), [b'q1\tACG\n', b'q2\tT\n'])
        self.assertFalse(self.output.exists())
        self.assertEqual(self.diagnostic()['returncode'], 0)

    def test_immediate_producer_completion_does_not_drop_buffered_bytes(self):
        data = b'query\tACGT\n' * 10000
        def runner(*args):
            with self.output.open('wb') as stream:
                stream.write(data)
            return {'status': 'COMPLETE', 'returncode': 0}
        for _ in range(10):
            with self.stream('pass', runner=runner) as lines:
                self.assertEqual(b''.join(lines), data)

    def test_empty_success_with_and_without_opening_output(self):
        for source in ("import sys; open(sys.argv[1],'wb').close()", 'pass'):
            with self.subTest(source=source):
                with self.stream(source) as lines:
                    self.assertEqual(list(lines), [])
                self.assertFalse(self.output.exists())
                self.assertEqual(self.diagnostic()['status'], 'COMPLETE')

    def test_failure_before_open_does_not_wait_forever(self):
        before = time.monotonic()
        with self.assertRaisesRegex(RuntimeError, 'return code 17'):
            with self.stream("import sys; print('missing database',file=sys.stderr); sys.exit(17)") as lines:
                list(lines)
        self.assertLess(time.monotonic() - before, 5)
        self.assertFalse(self.output.exists())
        self.assertEqual(self.diagnostic()['returncode'], 17)

    def test_complete_lines_do_not_hide_nonzero_exit(self):
        with self.assertRaisesRegex(RuntimeError, 'return code 5'):
            with self.stream("import sys; f=open(sys.argv[1],'wb'); f.write(b'partial\\n'); f.close(); sys.exit(5)") as lines:
                list(lines)
        self.assertEqual(self.diagnostic()['status'], 'INCOMPLETE')
        self.assertFalse(self.output.exists())

    def test_successful_exit_with_unterminated_line_is_incomplete(self):
        with self.assertRaisesRegex(RuntimeError, 'unterminated line'):
            with self.stream("import sys; open(sys.argv[1],'wb').write(b'truncated')") as lines:
                list(lines)
        self.assertFalse(self.output.exists())

    def test_oversized_line_stops_producer_with_bounded_buffer(self):
        self.assertEqual(MAX_BUFFER_BYTES, 1024 * 1024)
        source = "import sys,time; f=open(sys.argv[1],'wb',buffering=0); f.write(b'x'*(1024*1024+1)); time.sleep(30)"
        with self.assertRaisesRegex(RuntimeError, '1 MiB streaming buffer'):
            with self.stream(source) as lines:
                list(lines)
        self.assertFalse(self.output.exists())
        self.assertIsNotNone(self.diagnostic()['returncode'])

    def test_parser_exception_is_preserved_and_only_owned_child_stops(self):
        unrelated = subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(30)'])
        sentinel = ValueError('deliberate parser failure')
        try:
            source = "import sys,time; f=open(sys.argv[1],'wb',buffering=0); f.write(b'row\\n'); time.sleep(30)"
            with self.assertRaises(ValueError) as caught:
                with self.stream(source) as lines:
                    self.assertEqual(next(lines), b'row\n')
                    raise sentinel
            self.assertIs(caught.exception, sentinel)
            record = self.diagnostic()
            self.assertIsNotNone(record['returncode'])
            self.assertTrue(all(action['target_pid'] == record['pid'] for action in record['termination_actions']))
            self.assertIsNone(unrelated.poll())
            self.assertFalse(self.output.exists())
        finally:
            unrelated.terminate()
            unrelated.wait(timeout=5)

    def test_external_cancellation_stops_silent_hanging_producer(self):
        cancel = threading.Event()
        timer = threading.Timer(.2, cancel.set)
        timer.start()
        try:
            source = "import sys,time; f=open(sys.argv[1],'wb'); time.sleep(30)"
            before = time.monotonic()
            with self.assertRaises(Cancelled):
                with self.stream(source, cancel) as lines:
                    list(lines)
            self.assertLess(time.monotonic() - before, 5)
            self.assertFalse(self.output.exists())
            self.assertIsNotNone(self.diagnostic()['returncode'])
        finally:
            timer.cancel()
            timer.join()

    def test_resource_stop_keeps_original_exception(self):
        sentinel = RuntimeError('synthetic resource stop')
        class Guard:
            def checkpoint(self, *args, **kwargs):
                raise sentinel
        with self.assertRaises(RuntimeError) as caught:
            with self.stream('import time; time.sleep(30)', guard=Guard()) as lines:
                list(lines)
        self.assertIs(caught.exception, sentinel)
        record = self.diagnostic()
        self.assertEqual(record['status'], 'INCOMPLETE')
        self.assertEqual(record['stop_reason'], str(sentinel))
        self.assertIsNotNone(record['returncode'])

    def test_early_consumer_exit_is_not_success(self):
        source = "import sys,time; f=open(sys.argv[1],'wb',buffering=0); f.write(b'row\\n'); time.sleep(30)"
        with self.assertRaisesRegex(RuntimeError, 'not fully consumed'):
            with self.stream(source) as lines:
                self.assertEqual(next(lines), b'row\n')
        self.assertIsNotNone(self.diagnostic()['returncode'])
        self.assertFalse(self.output.exists())

    def test_runner_must_supply_verified_zero_exit(self):
        for result in (None, {'status': 'INCOMPLETE', 'returncode': 0}, {'status': 'COMPLETE', 'returncode': 1}):
            with self.subTest(result=result):
                with self.assertRaisesRegex(RuntimeError, 'verified successful exit'):
                    with self.stream('pass', runner=lambda *a: result) as lines:
                        list(lines)
                self.assertFalse(self.output.exists())

    def test_existing_output_path_is_never_overwritten(self):
        self.output.write_text('preserve me')
        with self.assertRaises(FileExistsError):
            with self.stream('pass') as lines:
                list(lines)
        self.assertEqual(self.output.read_text(), 'preserve me')

    def test_replaced_output_path_is_preserved(self):
        def runner(*args):
            self.output.unlink()
            self.output.write_text('replacement')
            return {'status': 'COMPLETE', 'returncode': 0}
        with self.assertRaisesRegex(RuntimeError, 'replacement was preserved'):
            with self.stream('pass', runner=runner) as lines:
                list(lines)
        self.assertEqual(self.output.read_text(), 'replacement')

    def test_cleanup_wait_is_bounded_for_noncooperative_test_runner(self):
        release = threading.Event()
        def runner(*args):
            release.wait(5)
            return {'status': 'COMPLETE', 'returncode': 0}
        try:
            with patch('blast_stream.JOIN_TIMEOUT_SECONDS', .05):
                with self.assertRaisesRegex(RuntimeError, 'not fully consumed') as caught:
                    with self.stream('pass', runner=runner):
                        pass
            self.assertTrue(any('did not stop' in note for note in getattr(caught.exception, '__notes__', [])))
        finally:
            release.set()


if __name__ == '__main__':
    unittest.main()
