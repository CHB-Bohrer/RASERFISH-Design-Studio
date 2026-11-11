"""Saved checkpoints bind immutable inputs; incomplete evidence is never a pass."""
import copy
import json
import os
from pathlib import Path
import sqlite3
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import design_checkpoints as checkpoints


class CheckpointIdentityTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name).resolve()
        self.path = self.root / 'checkpoints'
        self.config = dict(output_detail='standard', title='DP3',
                           references=[dict(id='target', sequence='ACGT' * 20)],
                           policy=dict(mode='local_preview', screen_cross_primary=False))
        self.candidates = [dict(id='candidate_1', arm='ACGT', full='TTACGTTT',
                                reference='target', bin_id='region_1', prefix='TT', suffix='TT',
                                sites=[dict(reference='target', start=2, end=5, strand='+')]),
                           dict(id='candidate_2', arm='GCTA', full='TTGCTATT',
                                reference='target', bin_id='region_1', prefix='TT', suffix='TT',
                                sites=[dict(reference='target', start=8, end=11, strand='+')])]
        self.policy = dict(mode='local_preview', screen_cross_primary=False, max_high_offtarget_sites=0)
        self.engine = dict(sources={'design_engine.py': 'fixture1'},
                           tools={'blastn': {'sha256': 'fixture2'}},
                           dependencies={'RNA': {'version': 'fixture3'}})

    def session(self, **overrides):
        kwargs = dict(path=self.path, config=self.config, candidates=self.candidates,
                      databases={}, policy=self.policy, engine_identity=self.engine)
        kwargs.update(overrides)
        return checkpoints.CheckpointSession(**kwargs)

    def initialized(self, **overrides):
        session = self.session(**overrides)
        with sqlite3.connect(session.store_path) as connection:
            connection.execute('CREATE TABLE fixture (value TEXT)')
        return session

    def test_identical_identity_resumes_without_rewriting_manifest(self):
        first = self.initialized()
        first.update_progress(completed_batches=2, total_batches=4, folded_candidates=1,
                              summary='Two batches committed.')
        before = first.manifest_path.read_bytes(), first.manifest_path.stat().st_mtime_ns
        second = self.session(config=dict(reversed(list(self.config.items()))))
        self.assertTrue(second.is_resume)
        self.assertEqual(first.identity, second.identity)
        self.assertEqual(before, (second.manifest_path.read_bytes(), second.manifest_path.stat().st_mtime_ns))
        self.assertEqual(second.progress['completed_batches'], 2)

    def test_config_sequence_change_fails_closed(self):
        self.initialized()
        changed = copy.deepcopy(self.config)
        changed['references'][0]['sequence'] += 'A'
        with self.assertRaisesRegex(ValueError, 'identity mismatch.*config_sha256'):
            self.session(config=changed)

    def test_candidate_order_sequence_and_intended_sites_bind_identity(self):
        self.initialized()
        changed_sets = [list(reversed(self.candidates))]
        for field, value in [('arm', 'AAAA'), ('full', 'TTAAAATT'), ('prefix', 'AA'),
                             ('reference', 'other'), ('bin_id', 'region_2'),
                             ('sites', [dict(reference='target', start=3, end=6, strand='+')])]:
            changed = copy.deepcopy(self.candidates)
            changed[0][field] = value
            changed_sets.append(changed)
        for changed in changed_sets:
            with self.subTest(changed=changed), self.assertRaisesRegex(ValueError, 'identity mismatch.*candidates'):
                self.session(candidates=changed)

    def test_mutable_qc_results_do_not_change_input_identity(self):
        self.initialized()
        changed = copy.deepcopy(self.candidates)
        changed[0].update(offtarget_hits=object(), reject_reasons=['partial'], fold={'status': 'partial'})
        self.assertTrue(self.session(candidates=changed).is_resume)

    def test_policy_and_code_and_tool_changes_fail_closed(self):
        self.initialized()
        changed_policy = dict(self.policy, max_high_offtarget_sites=1)
        with self.assertRaisesRegex(ValueError, 'identity mismatch.*policy'):
            self.session(policy=changed_policy)
        for section in ('sources', 'tools', 'dependencies'):
            changed = copy.deepcopy(self.engine)
            changed[section]['new_version'] = 'different'
            with self.subTest(section=section), self.assertRaisesRegex(ValueError, 'identity mismatch.*engine'):
                self.session(engine_identity=changed)

    def test_database_metadata_change_fails_without_reading_genome_contents(self):
        prefix = self.root / 'human'
        index = self.root / 'human.nsq'
        index.write_bytes(b'fixture-index')
        self.initialized(databases={'human': str(prefix)})
        original = Path.open
        def guarded_open(path, *args, **kwargs):
            if path == index:
                raise AssertionError('Genome index contents must not be read for checkpoint identity')
            return original(path, *args, **kwargs)
        with patch.object(Path, 'open', guarded_open):
            self.assertTrue(self.session(databases={'human': str(prefix)}).is_resume)
        info = index.stat()
        os.utime(index, ns=(info.st_atime_ns, info.st_mtime_ns + 1000000))
        with self.assertRaisesRegex(ValueError, 'identity mismatch.*databases'):
            self.session(databases={'human': str(prefix)})

    def test_alias_database_identity_includes_external_volume(self):
        volume = self.root / 'volumes'
        volume.mkdir()
        (volume / 'piece.nsq').write_bytes(b'genome-volume')
        (self.root / 'human.nal').write_text('TITLE test\nDBLIST "volumes/piece"\n')
        identity = checkpoints.database_identity(self.root / 'human')
        paths = [item['path'] for item in identity['index_files']]
        self.assertEqual(paths, sorted([str(self.root / 'human.nal'), str(volume / 'piece.nsq')]))
        self.assertTrue(all({'bytes', 'mtime_ns', 'ctime_ns', 'inode', 'device'} <= set(item)
                            for item in identity['index_files']))

    def test_missing_database_index_is_not_a_valid_identity(self):
        with self.assertRaisesRegex(ValueError, 'No local BLAST index'):
            self.session(databases={'human': str(self.root / 'missing')})

    def test_missing_store_and_missing_manifest_refuse_resume(self):
        self.session()
        with self.assertRaisesRegex(ValueError, 'store is missing'):
            self.session()
        other = self.root / 'other'
        other.mkdir()
        (other / 'retained_hits.sqlite').write_bytes(b'old-unidentified-evidence')
        with self.assertRaisesRegex(ValueError, 'without its identity'):
            self.session(path=other)

    def test_tampered_manifest_or_progress_refuses_resume(self):
        first = self.initialized()
        manifest = json.loads(first.manifest_path.read_text())
        original = first.manifest_path.read_bytes()
        manifest['identity']['config_sha256'] = 'changed'
        first.manifest_path.write_text(json.dumps(manifest))
        with self.assertRaisesRegex(ValueError, 'identity metadata is damaged'):
            self.session()
        first.manifest_path.write_bytes(original)
        progress = json.loads(first.progress_path.read_text())
        progress['identity_sha256'] = 'different'
        first.progress_path.write_text(json.dumps(progress))
        with self.assertRaisesRegex(ValueError, 'progress does not match'):
            self.session()

    def test_progress_always_incomplete_even_if_caller_claims_complete(self):
        first = self.initialized()
        first.update_progress(status='COMPLETE', completed_batches=4, total_batches=4,
                              folded_candidates=2, total_candidates=2,
                              summary='Individual QC saved; final export still requires verification.')
        self.assertEqual(first.progress['status'], 'INCOMPLETE')
        summary = checkpoints.read_checkpoint_summary(self.path)
        self.assertEqual(summary['status'], 'INCOMPLETE')
        self.assertEqual(summary['completed_batches'], 4)
        self.assertTrue(summary['resumable'])

    def test_full_audit_preserves_artifacts_and_explicitly_refuses_resume(self):
        config = dict(self.config, output_detail='full_audit')
        first = self.initialized(config=config)
        self.assertFalse(first.progress['resumable'])
        before = first.manifest_path.read_bytes()
        with self.assertRaisesRegex(ValueError, 'Full-audit checkpoints cannot be resumed'):
            self.session(config=config)
        self.assertEqual(before, first.manifest_path.read_bytes())

    def test_summary_reads_only_bounded_json_and_never_store(self):
        first = self.initialized()
        original = Path.open
        def guarded_open(path, *args, **kwargs):
            if path == first.store_path:
                raise AssertionError('Listing checkpoints must not read the retained-hit store')
            return original(path, *args, **kwargs)
        with patch.object(Path, 'open', guarded_open):
            self.assertTrue(checkpoints.read_checkpoint_summary(self.path)['resumable'])
        first.progress_path.write_bytes(b' ' * (checkpoints.MAX_METADATA_BYTES + 1))
        summary = checkpoints.read_checkpoint_summary(self.path)
        self.assertEqual(summary['status'], 'INCOMPLETE')
        self.assertFalse(summary['resumable'])
        self.assertIn('bounded size', summary['resume_reason'])

    def test_summary_absent_unknown_or_missing_store_is_not_resumable(self):
        self.assertIsNone(checkpoints.read_checkpoint_summary(self.path))
        first = self.session()
        self.assertFalse(checkpoints.read_checkpoint_summary(self.path)['resumable'])
        first.manifest_path.write_text('{"schema":"old-unverified-format"}')
        self.assertIsNone(checkpoints.read_checkpoint_summary(self.path))

    def test_summary_rejects_corrupt_identity_and_false_completion_claim(self):
        first = self.initialized()
        progress = json.loads(first.progress_path.read_text())
        progress['status'] = 'COMPLETE'
        first.progress_path.write_text(json.dumps(progress))
        summary = checkpoints.read_checkpoint_summary(self.path)
        self.assertFalse(summary['resumable'])
        self.assertEqual(summary['status'], 'INCOMPLETE')
        first.update_progress()
        manifest = json.loads(first.manifest_path.read_text())
        manifest['identity_sha256'] = 'invalid'
        first.manifest_path.write_text(json.dumps(manifest))
        self.assertFalse(checkpoints.read_checkpoint_summary(self.path)['resumable'])

    def test_symlinked_metadata_or_store_is_refused(self):
        first = self.initialized()
        external = self.root / 'external.sqlite'
        first.store_path.rename(external)
        first.store_path.symlink_to(external)
        with self.assertRaisesRegex(ValueError, 'store is missing or unsafe'):
            self.session()
        self.assertFalse(checkpoints.read_checkpoint_summary(self.path)['resumable'])

    def test_new_checkpoint_refuses_a_dangling_store_symlink(self):
        self.path.mkdir()
        target = self.root / 'not-created.sqlite'
        (self.path / 'retained_hits.sqlite').symlink_to(target)
        with self.assertRaisesRegex(ValueError, 'without its identity'):
            self.session()
        self.assertFalse(target.exists())
        self.assertFalse((self.path / 'checkpoint_manifest.json').exists())

    def test_atomic_writes_sync_file_and_parent_and_leave_no_temporary_files(self):
        with patch.object(checkpoints.os, 'fsync', wraps=os.fsync) as fsync:
            first = self.initialized()
            first.update_progress(completed_batches=1)
        self.assertEqual(fsync.call_count, 6)
        self.assertEqual(list(self.path.glob('*.tmp')), [])

    def test_bad_counts_and_nonfinite_inputs_are_rejected(self):
        first = self.initialized()
        before = first.progress_path.read_bytes()
        for value in (-1, True, 1.5):
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, 'nonnegative integer'):
                first.update_progress(completed_batches=value)
        self.assertEqual(before, first.progress_path.read_bytes())
        with self.assertRaises(ValueError):
            self.session(config=dict(self.config, threshold=float('nan')))


if __name__ == '__main__':
    unittest.main()
