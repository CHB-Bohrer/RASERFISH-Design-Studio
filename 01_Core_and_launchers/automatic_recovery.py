"""Bounded, durable automatic resume policy. Scientific evidence stays in checkpoints."""
from contextlib import contextmanager
import errno
import fcntl
import json
import math
import os
from pathlib import Path
import re
import stat
import time
from process_diagnostics import atomic_json, append_event, utc_now

SCHEMA='raserfish.automatic-recovery.v1'
MAX_ATTEMPTS=2
RETRY_DELAYS=(30,120)
ACTIVE_STATES={'active','pending'}
ID_PATTERN=re.compile(r'[a-z0-9_-]{1,100}')


def read_record(origin):
    path=Path(origin)/'automatic_recovery.json'
    if not path.exists():return None
    if path.is_symlink():raise ValueError('Unsafe automatic recovery record')
    with path.open('rb') as stream:blob=stream.read(65537)
    if len(blob)>65536:raise ValueError('Automatic recovery record exceeds its safe read limit')
    value=json.loads(blob)
    if (not isinstance(value,dict) or value.get('schema')!=SCHEMA or value.get('origin_run_id')!=Path(origin).name
        or not ID_PATTERN.fullmatch(str(value.get('latest_run_id','')))
        or value.get('state') not in {'active','pending','stopped','complete'}
        or type(value.get('attempts')) is not int or not 0<=value['attempts']<=MAX_ATTEMPTS):
        raise ValueError('Invalid automatic recovery record')
    when=value.get('next_retry_at')
    if value['state']=='pending' and (type(when) not in (int,float) or not math.isfinite(when) or when<0):
        raise ValueError('Invalid automatic recovery retry time')
    return value


@contextmanager
def record_lock(origin):
    """Serialize scheduling, cancellation and claims across local Studio instances."""
    origin=Path(origin)
    if origin.is_symlink() or not origin.is_dir():raise ValueError('Unsafe recovery run directory')
    path=origin/'automatic_recovery.lock'
    fd=os.open(path,os.O_CREAT|os.O_RDWR|getattr(os,'O_NOFOLLOW',0),0o600)
    try:
        info=os.fstat(fd);linked=path.stat(follow_symlinks=False)
        if not stat.S_ISREG(info.st_mode) or info.st_uid!=os.getuid() or info.st_mode&0o077 or (info.st_ino,info.st_dev)!=(linked.st_ino,linked.st_dev):
            raise ValueError('Unsafe recovery lock')
        fcntl.flock(fd,fcntl.LOCK_EX)
        yield
    finally:os.close(fd)


def write_record(origin,record,event,**updates):
    record=dict(record,**updates,updated_utc=utc_now())
    atomic_json(Path(origin)/'automatic_recovery.json',record)
    append_event(Path(origin)/'automatic_recovery_events.jsonl',dict(event=event,**record))
    return record


def enroll(origin,job_id,reset=False):
    with record_lock(origin):
        old=read_record(origin)
        if old and not reset:return old
        return write_record(origin,dict(schema=SCHEMA,origin_run_id=Path(origin).name,latest_run_id=job_id,
            state='active',attempts=0,max_attempts=MAX_ATTEMPTS,next_retry_at=None,
            reason='Automatic checkpoint recovery enabled for this explicitly started design.'),'ENROLLED')


def stop(origin,reason,job_id=None,complete=False,expected=None):
    with record_lock(origin):
        record=read_record(origin)
        if not record or (job_id and record['latest_run_id']!=job_id):return record
        if expected and any(record.get(k)!=expected.get(k) for k in ('latest_run_id','attempts','state','next_retry_at')):return record
        if record['state']=='complete':return record
        return write_record(origin,record,'COMPLETED' if complete else 'STOPPED',
            state='complete' if complete else 'stopped',next_retry_at=None,reason=str(reason)[:4000])


def retry_block_reason(error,exception_type='',memory_sample=None,cancelled=False):
    if cancelled:return 'Cancelled by user; automatic recovery is stopped.'
    if memory_sample or exception_type in ('MemoryLimitExceeded','MemoryMetricsUnavailable','DiskSpaceLimitExceeded','MemoryError'):
        return 'Resource safety stop; automatic recovery is stopped. Review resource diagnostics before resuming.'
    if exception_type in ('ValueError','TypeError','AssertionError','PermissionError','FileNotFoundError','DatabaseError','IntegrityError','OperationalError'):
        return 'Inputs, evidence or local files require review before resuming: '+str(error)[:3000]
    message=str(error).lower()
    if any(word in message for word in ('disk reserve','disk space could not be checked','memory safety','memory limit','memory metrics','memory pressure','memoryguard','memory guard','out of memory','cannot allocate memory','bad_alloc')):
        return 'Resource safety stop; automatic recovery is stopped. '+str(error)[:3000]
    if any(word in message for word in ('checkpoint identity mismatch','checkpoint identity metadata is damaged','checksum','sha-256 mismatch','database disk image is malformed','database or disk is full','no space left','permission denied','operation not permitted','cleanup failed','did not exit after sigkill','cannot preserve worker output')):
        return 'Saved evidence or local resources require review before resuming: '+str(error)[:3000]
    return None


def schedule(origin,job_id,reason,*,blocked=None,now=None):
    now=time.time() if now is None else now
    with record_lock(origin):
        record=read_record(origin)
        if not record or record['latest_run_id']!=job_id or record['state']!='active':return record
        if blocked or record['attempts']>=MAX_ATTEMPTS:
            return write_record(origin,record,'STOPPED',state='stopped',next_retry_at=None,
                reason=blocked or 'Automatic retry limit reached (2 retries). Analysis remains incomplete. '+str(reason)[:3000])
        return write_record(origin,record,'RETRY_SCHEDULED',state='pending',
            next_retry_at=now+RETRY_DELAYS[record['attempts']],reason=str(reason)[:4000])


def claim(origin,new_id,*,now=None):
    """Caller must already own the global single-design lease."""
    now=time.time() if now is None else now
    with record_lock(origin):
        record=read_record(origin)
        if not record or record['state']!='pending' or record['next_retry_at']>now:
            raise ValueError('Automatic recovery is no longer pending or not yet due')
        if record['attempts']>=MAX_ATTEMPTS:raise ValueError('Automatic recovery retry limit reached')
        return write_record(origin,record,'RETRY_STARTED',state='active',latest_run_id=new_id,
            attempts=record['attempts']+1,next_retry_at=None,
            reason='Automatically resuming verified checkpoints; analysis remains incomplete.')
