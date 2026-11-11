"""Automatic retry policy and real local HTTP/lease integration, with no full design."""
import copy
import json
import os
from pathlib import Path
import sys
import threading
import time
import unittest
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import automatic_recovery as recovery
import runtime_safety
import studio
import test_checkpoint_server as fixtures
NoRunThread=fixtures.NoRunThread

class AutomaticRecoveryTests(unittest.TestCase):
    setUp=fixtures.CheckpointServerTests.setUp
    saved=fixtures.CheckpointServerTests.saved
    snapshot=fixtures.CheckpointServerTests.snapshot
    start_http=fixtures.CheckpointServerTests.start_http
    request=fixtures.CheckpointServerTests.request
    def tearDown(self):
        studio.SERVER_STOPPING.clear();studio.RECOVERY_OWNED.clear()
        fixtures.CheckpointServerTests.tearDown(self)
    def enroll(self):
        root=self.saved();recovery.enroll(root,root.name);return root
    def pending(self):
        root=self.enroll();recovery.schedule(root,root.name,'Worker interrupted',now=0);return root
    def test_only_explicitly_enrolled_runs_recover(self):
        root=self.saved();before=self.snapshot(self.jobs)
        with patch.object(studio,'start_design_job') as start:studio.recovery_tick();start.assert_not_called()
        self.assertEqual(before,self.snapshot(self.jobs))
    def test_backoff_is_durable_and_total_retries_bounded_even_after_reloading(self):
        root=self.enroll();record=recovery.schedule(root,root.name,'exit 7',now=100)
        self.assertEqual((record['state'],record['next_retry_at']),( 'pending',130))
        with self.assertRaisesRegex(ValueError,'not yet due'):recovery.claim(root,'attempt1',now=129)
        self.assertEqual(recovery.claim(root,'attempt1',now=130)['attempts'],1)
        self.assertEqual(recovery.schedule(root,'attempt1','exit 7',now=200)['next_retry_at'],320)
        self.assertEqual(recovery.claim(root,'attempt2',now=320)['attempts'],2)
        record=recovery.schedule(root,'attempt2','exit 7',now=400)
        self.assertEqual(record['state'],'stopped');self.assertIn('limit reached',record['reason'])
        with self.assertRaises(ValueError):recovery.claim(root,'attempt3',now=1000)
        self.assertEqual(recovery.read_record(root)['attempts'],2)
    def test_duplicate_completion_cannot_schedule_or_advance_twice(self):
        root=self.enroll();recovery.schedule(root,root.name,'failed',now=100)
        record=recovery.schedule(root,root.name,'duplicate',now=200)
        self.assertEqual(record['next_retry_at'],130)
        recovery.claim(root,'retry',now=130)
        record=recovery.schedule(root,root.name,'stale callback',now=999)
        self.assertEqual(record['state'],'active');self.assertEqual(record['latest_run_id'],'retry')
    def test_cancellation_durable_and_pending_attempt_cannot_be_claimed(self):
        root=self.pending();studio.cancel_design_chain(root.name)
        self.assertEqual(recovery.read_record(root)['state'],'stopped')
        with self.assertRaises(ValueError):recovery.claim(root,'forbidden',now=999999)
        with patch.object(studio,'start_design_job') as start:studio.recovery_tick();start.assert_not_called()
    def test_cancel_old_id_targets_active_retry_and_preserves_failed_status(self):
        root=self.pending();retry=self.saved('retry',origin=root.name)
        recovery.claim(root,retry.name,now=999999)
        event=threading.Event();self.state[retry.name]=dict(id=retry.name,status='running',cancel=event,checkpoint_origin=root.name)
        studio.cancel_design_chain(root.name)
        self.assertTrue(event.is_set());self.assertEqual(recovery.read_record(root)['state'],'stopped')
        self.assertEqual(json.loads((root/'analysis_status.json').read_text())['status'],'incomplete')
    def test_resource_safety_validation_permission_and_storage_stops_are_not_retried(self):
        cases=[('limit','MemoryLimitExceeded',{}),('metrics','MemoryMetricsUnavailable',{}),
            ('bad checkpoint','ValueError',None),('Permission denied','RuntimeError',None),
            ('No space left on device','OSError',None),('cannot allocate memory','RuntimeError',None),
            ('worker cleanup failed','RuntimeError',None),('failed','MemoryError',None)]
        for error,kind,sample in cases:
            with self.subTest(error=error):self.assertIsNotNone(recovery.retry_block_reason(error,kind,sample))
        self.assertIsNone(recovery.retry_block_reason('Subprocess returned exit code 7','RuntimeError'))
    def test_restart_waits_for_owned_worker_lease_and_preserves_state(self):
        root=self.enroll();lease=runtime_safety.acquire_design_lock(dict(job_id=root.name))
        try:
            before=(root/'automatic_recovery.json').read_bytes()
            studio.recovery_tick();self.assertEqual((root/'automatic_recovery.json').read_bytes(),before)
        finally:lease.release()
        studio.recovery_tick();record=recovery.read_record(root)
        self.assertEqual(record['state'],'pending');self.assertEqual(record['attempts'],0)
    def test_missing_checkpoint_stops_incomplete_without_new_work(self):
        root=self.enroll();self.summaries.clear()
        studio.recovery_tick();record=recovery.read_record(root)
        self.assertEqual(record['state'],'stopped');self.assertIn('No verified-format',record['reason'])
        self.assertEqual(self.state,{})
    def test_automatic_claim_uses_original_configuration_and_same_cross_instance_lease(self):
        root=self.pending();original=self.snapshot(root/'result')
        with patch.object(studio.threading,'Thread',NoRunThread):studio.recovery_tick()
        record=recovery.read_record(root);self.assertEqual(record['attempts'],1)
        ident=record['latest_run_id'];self.assertNotEqual(ident,root.name)
        self.assertEqual(self.state[ident]['configuration'],self.config)
        self.assertTrue(self.state[ident]['automatic_attempt']);self.assertEqual(original,self.snapshot(root/'result'))
        with self.assertRaises(runtime_safety.DesignAlreadyRunning):runtime_safety.acquire_design_lock()
        reply=studio.job_with_recovery(root.name)
        self.assertEqual(reply['status'],'incomplete');self.assertEqual(reply['auto_recovery']['latest_run_id'],ident)
    def test_fresh_design_refused_while_original_awaits_retry(self):
        root=self.pending()
        with self.assertRaisesRegex(ValueError,'awaiting automatic recovery'):studio.start_design_job(self.config)
        self.assertEqual(self.state,{})
    def test_input_hash_change_stops_before_spawn(self):
        root=self.pending();path=root/'result/run_config.json';m=json.loads(path.read_text())
        m['configuration']['policy']['tm_min']=44;path.write_text(json.dumps(m))
        with patch.object(studio,'start_design_job') as start:studio.recovery_tick();start.assert_not_called()
        record=recovery.read_record(root);self.assertEqual(record['state'],'stopped');self.assertIn('SHA-256 mismatch',record['reason'])
    def test_memory_worker_record_is_honored_after_server_restart(self):
        root=self.enroll();(root/'worker_status.json').write_text(json.dumps(dict(status='incomplete',error='guard fired',exception_type='MemoryLimitExceeded')))
        studio.recovery_tick();record=recovery.read_record(root)
        self.assertEqual(record['state'],'stopped');self.assertIn('Resource safety',record['reason'])
    def test_http_discovery_does_not_launch_and_cancel_auth_is_required(self):
        root=self.pending();before=recovery.read_record(root)
        with patch.object(studio,'start_design_job') as start:
            code,reply=self.request('/api/active-design');self.assertEqual(code,200)
            self.assertEqual(reply['id'],root.name);self.assertEqual(reply['config'],self.config);self.assertEqual(reply['token'],studio.TOKEN)
            code,listing=self.request('/api/checkpoints');self.assertEqual(code,200)
            self.assertTrue(listing['active_design']);self.assertEqual(listing['active_run']['status'],'retry_wait')
            self.assertEqual(self.request('/api/cancel',{'id':root.name})[0],403)
            self.assertEqual(recovery.read_record(root),before);start.assert_not_called()
            self.assertEqual(self.request('/api/cancel',{'id':root.name},{'X-Studio-Token':studio.TOKEN})[0],200)
            self.assertEqual(recovery.read_record(root)['state'],'stopped')
    def test_shutdown_blocks_automatic_claim(self):
        self.pending();studio.SERVER_STOPPING.set()
        with patch.object(studio,'start_design_job') as start:studio.recovery_tick();start.assert_not_called()
        with self.assertRaisesRegex(ValueError,'shutting down'):studio.start_design_job(self.config)
    def test_completed_status_without_export_commit_never_becomes_complete(self):
        root=self.enroll();(root/'analysis_status.json').write_text(json.dumps(dict(status='complete',analysis_complete=True,job_id=root.name)))
        (root/'job_status.json').write_text(json.dumps(dict(status='complete',id=root.name)))
        self.assertFalse(studio.restore_completed_recovery(root,recovery.read_record(root)))
        self.assertEqual(recovery.read_record(root)['state'],'active')
    def test_committed_final_export_survives_server_loss_without_rerunning(self):
        root=self.enroll();(root/'analysis_status.json').write_text(json.dumps(dict(status='complete',analysis_complete=True,job_id=root.name)))
        (root/'job_status.json').write_text(json.dumps(dict(status='complete',id=root.name,messages=[])))
        (self.jobs/(root.name+'.zip')).write_bytes(b'fixture: committed archive')
        studio.recovery_tick()
        self.assertEqual(recovery.read_record(root)['state'],'complete');self.assertEqual(self.state[root.name]['status'],'complete')
    def test_cancel_immediately_after_manual_resume_targets_new_attempt(self):
        root=self.pending()
        with patch.object(studio.threading,'Thread',NoRunThread):ident,_=studio.start_design_job(self.config,resume_from=root.name)
        self.assertEqual(recovery.read_record(root)['latest_run_id'],ident)
        studio.cancel_design_chain(root.name)
        self.assertTrue(self.state[ident]['cancel'].is_set())
        recovery.enroll(root,ident)
        self.assertEqual(recovery.read_record(root)['state'],'stopped')
    def test_completed_recovery_is_readable_after_another_server_restart(self):
        root=self.enroll();(root/'analysis_status.json').write_text(json.dumps(dict(status='complete',analysis_complete=True,job_id=root.name)))
        (root/'job_status.json').write_text(json.dumps(dict(status='complete',id=root.name,messages=[])))
        (self.jobs/(root.name+'.zip')).write_bytes(b'committed archive')
        recovery.stop(root,'Complete',complete=True)
        reply=studio.job_with_recovery(root.name)
        self.assertEqual(reply['status'],'complete');self.assertEqual(reply['auto_recovery']['state'],'complete')
    def test_symlink_or_oversized_record_is_rejected(self):
        root=self.enroll();path=root/'automatic_recovery.json';outside=self.base/'outside';outside.write_bytes(path.read_bytes())
        path.unlink();path.symlink_to(outside)
        with self.assertRaisesRegex(ValueError,'Unsafe'):recovery.read_record(root)
        path.unlink();path.write_bytes(b' '*65537)
        with self.assertRaisesRegex(ValueError,'safe read limit'):recovery.read_record(root)

    def test_cancel_racing_claim_targets_new_attempt_returned_by_atomic_stop(self):
        root=self.pending();retry=self.saved('racing_retry',origin=root.name)
        event=threading.Event();self.state[retry.name]=dict(id=retry.name,status='running',cancel=event)
        original_stop=recovery.stop
        def race(*args,**kwargs):
            recovery.claim(root,retry.name,now=999999)
            return original_stop(*args,**kwargs)
        with patch.object(recovery,'stop',side_effect=race):studio.cancel_design_chain(root.name)
        self.assertTrue(event.is_set());self.assertEqual(recovery.read_record(root)['state'],'stopped')
    def test_cancellation_from_another_instance_is_observed_by_owned_worker(self):
        root=self.enroll();event=threading.Event();item=dict(id=root.name,cancel=event)
        recovery.stop(root,'Cancelled from another browser')
        with self.assertRaisesRegex(RuntimeError,'Cancelled from another browser'):
            studio.observe_recovery_cancellation(item,root)
        self.assertTrue(event.is_set());self.assertTrue(item['user_cancelled'])
    def test_obsolete_retry_failure_does_not_stop_newer_attempt(self):
        root=self.pending();before=recovery.read_record(root)
        recovery.claim(root,'newer_attempt',now=999999)
        recovery.stop(root,'obsolete exception',expected=before)
        record=recovery.read_record(root);self.assertEqual(record['state'],'active');self.assertEqual(record['latest_run_id'],'newer_attempt')
    def test_invalid_pending_timestamps_refused(self):
        root=self.pending();path=root/'automatic_recovery.json';original=json.loads(path.read_text())
        for when in (None,'later',float('nan'),float('inf'),-1,True):
            with self.subTest(when=when):
                path.write_text(json.dumps(dict(original,next_retry_at=when)))
                with self.assertRaisesRegex(ValueError,'retry time'):recovery.read_record(root)
    def test_supervisor_memory_stop_is_honored_after_server_loss(self):
        root=self.enroll();(root/'supervisor_process.json').write_text(json.dumps(dict(stop_reason='supervisor guard stop',memory_sample={'status':'STOP'})))
        studio.recovery_tick();self.assertEqual(recovery.read_record(root)['state'],'stopped')
    def test_completion_cannot_be_borrowed_from_an_unrelated_chain(self):
        root=self.enroll();other=self.saved('unrelated',status='complete');record=recovery.read_record(root)
        record['latest_run_id']=other.name
        (other/'analysis_status.json').write_text(json.dumps(dict(status='complete',analysis_complete=True,job_id=other.name)))
        (other/'job_status.json').write_text(json.dumps(dict(status='complete',id=other.name)))
        (self.jobs/(other.name+'.zip')).write_bytes(b'committed')
        self.assertFalse(studio.restore_completed_recovery(root,record));self.assertEqual(self.state,{})

    def test_idle_second_server_shutdown_does_not_cancel_other_attempt(self):
        root=self.enroll()
        studio.stop_owned_recovery_chains('Idle server intentionally stopped')
        self.assertEqual(recovery.read_record(root)['state'],'active')
    def test_old_server_shutdown_does_not_cancel_newly_claimed_attempt(self):
        root=self.pending();studio.RECOVERY_OWNED.add((root.name,root.name))
        recovery.claim(root,'new_owner_attempt',now=999999)
        studio.stop_owned_recovery_chains('Old server intentionally stopped')
        self.assertEqual(recovery.read_record(root)['state'],'active')
    def test_shutdown_stop_is_atomic_against_a_racing_new_claim(self):
        root=self.pending();studio.RECOVERY_OWNED.add((root.name,root.name))
        original=recovery.stop
        def race(*args,**kwargs):
            recovery.claim(root,'new_owner_attempt',now=999999)
            return original(*args,**kwargs)
        with patch.object(recovery,'stop',side_effect=race):studio.stop_owned_recovery_chains('Old server stopped')
        record=recovery.read_record(root);self.assertEqual(record['state'],'active');self.assertEqual(record['latest_run_id'],'new_owner_attempt')
    def test_shutdown_stops_its_own_pending_attempt(self):
        root=self.pending();studio.RECOVERY_OWNED.add((root.name,root.name))
        studio.stop_owned_recovery_chains('Designer intentionally closed')
        self.assertEqual(recovery.read_record(root)['state'],'stopped')

if __name__=='__main__':unittest.main()
