import io,json,sys,tempfile,unittest,zipfile
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import studio, runtime_safety as safety, automatic_recovery as recovery
class CompactBackendTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.base=Path(self.temp.name);self.jobs=self.base/'runs';self.jobs.mkdir()
        self.patch=patch.object(studio,'JOBS',self.jobs);self.patch.start()
    def tearDown(self):self.patch.stop();self.temp.cleanup()
    def completed(self):
        ident='20261018_000000_example';root=self.jobs/ident;(root/'result').mkdir(parents=True)
        (root/'analysis_status.json').write_text(json.dumps(dict(status='complete',analysis_complete=True,job_id=ident)))
        (root/'job_status.json').write_text(json.dumps(dict(status='complete',id=ident)))
        (root/'result/design_result.json').write_text(json.dumps(dict(status='DRAFT',counts=dict(selected=1),bins=[dict(id='x',reference='chr1',start=1,end=50,imager='DP3',chromosome='chr1',genomic_start=1,genomic_end=50)],policy=dict(screen_cross_primary=False))))
        with zipfile.ZipFile(self.jobs/(ident+'.zip'),'w') as z:z.writestr('score','DRAFT')
        return root
    def test_completed_runs_survive_empty_server_state_and_use_actual_imager(self):
        root=self.completed()
        with patch.object(studio,'STATE',{}):
            rows=studio.completed_runs()['runs'];self.assertEqual(len(rows),1);self.assertIn('DP3',rows[0]['name'])
            self.assertEqual(studio.completed_run(root.name,True)['result']['run_name'],rows[0]['name'])
    def test_incomplete_or_uncommitted_archives_never_list(self):
        root=self.completed();(self.jobs/(root.name+'.INCOMPLETE_EXPORT.json')).write_text('{}')
        self.assertEqual(studio.completed_runs()['runs'],[])
        with self.assertRaises(ValueError):studio.completed_run(root.name)
    def test_status_and_run_identity_must_match(self):
        root=self.completed();(root/'analysis_status.json').write_text(json.dumps(dict(status='complete',analysis_complete=True,job_id='different')))
        self.assertEqual(studio.completed_runs()['runs'],[])
    def test_saved_path_cannot_escape_runs(self):
        with self.assertRaises(ValueError):studio.completed_run('../outside')
    def test_partial_payload_cleanup_preserves_error_diagnostics_and_symlinks(self):
        root=self.completed();attempt=root/'result/blast_attempts/human-a';attempt.mkdir(parents=True)
        for name in ('0001_candidates_blast.tsv','0001_queries.fa','command.json','human_blast.log'):(attempt/name).write_text('fixture')
        outside=self.base/'outside';outside.write_text('keep');(attempt/'link_candidates_blast.tsv').symlink_to(outside)
        deleted=studio.prune_interrupted_standard_batches(root)
        self.assertEqual(len(deleted),2);self.assertTrue((attempt/'command.json').exists());self.assertTrue((attempt/'human_blast.log').exists());self.assertEqual(outside.read_text(),'keep')
    def test_memory_logs_reduce_to_peak_scores(self):
        root=self.completed();path=root/'memory_usage.jsonl'
        path.write_text('\n'.join(json.dumps(dict(status='OK',process_family_rss_bytes=x)) for x in (1,500,20)))
        summary=studio.compact_memory_logs(root)
        self.assertEqual(summary['memory_usage.jsonl']['peaks']['process_family_rss_bytes'],500)
        self.assertFalse(path.exists());self.assertTrue((root/'result/runtime_summary.json').exists())
    def test_standard_archive_omits_runtime_snapshots(self):
        root=self.completed();(root/'memory_usage.jsonl').write_text('many samples')
        self.assertNotIn(root/'memory_usage.jsonl',list(studio.result_archive_files(root)))
    def test_full_audit_archive_retains_opt_in_runtime_evidence(self):
        root=self.completed();(root/'result/run_config.json').write_text(json.dumps(dict(output_detail='full_audit')));path=root/'memory_usage.jsonl';path.write_text('many samples')
        self.assertIn(path,list(studio.result_archive_files(root)))
    def test_disk_guard_stops_and_records_without_starting_work(self):
        record=dict(metrics_available=True,limits=safety.memory_limits(32*safety.GIB),process_family_rss_bytes=1000)
        usage=type('Usage',(),{'free':500_000_000})()
        guard=safety.MemoryGuard(self.base/'memory.jsonl')
        with patch.object(guard,'sample',return_value=record),patch.object(safety.shutil,'disk_usage',return_value=usage):
            with self.assertRaises(safety.DiskSpaceLimitExceeded) as caught:guard.check('fixture')
        self.assertIn('disk reserve',str(caught.exception));self.assertEqual(json.loads((self.base/'memory.jsonl').read_text())['status'],'STOP_FREE_DISK_LIMIT')
        self.assertIsNotNone(recovery.retry_block_reason(str(caught.exception),type(caught.exception).__name__))
    def test_disk_guard_does_not_restore_removed_ram_reserve(self):
        record=dict(metrics_available=True,limits=safety.memory_limits(32*safety.GIB),process_family_rss_bytes=1000,available_estimated_bytes=1)
        usage=type('Usage',(),{'free':10*safety.GIB})()
        guard=safety.MemoryGuard(self.base/'memory.jsonl')
        with patch.object(guard,'sample',return_value=record),patch.object(safety.shutil,'disk_usage',return_value=usage):self.assertEqual(guard.check('fixture')['status'],'OK')
if __name__=='__main__':unittest.main()
