#!/usr/bin/env python3
"""Local browser application for QC-guided RASER-FISH selection. Bind loopback only."""
from __future__ import annotations
import argparse
import base64
import copy
import gzip
import io
import importlib.util
import json
import mimetypes
import os
import re
import secrets
import signal
import shutil
import subprocess
import sys
import threading
import tempfile
import traceback
import time
import urllib.parse
import uuid
import webbrowser
import zipfile
from datetime import datetime,timezone
from http.server import BaseHTTPRequestHandler,ThreadingHTTPServer
from pathlib import Path
from sequence_tools import parse_fasta,parse_annotated,validate_imagers,make_bins
import annotation_sources as annot
import coverage_report
from design_presets import PRESETS, PARAMETER_NOTES, ERRFI1_DEFAULTS
from rna_avoidance import feature_intervals
from process_diagnostics import (atomic_json, append_event, bounded_tail, command_record,
    returncode_fields, failure_message, finalize_interrupted_commands, utc_now)

BASE=Path(__file__).resolve().parent
JOBS=BASE/'design_runs'
TOKEN=secrets.token_urlsafe(32)
SERVER_INSTANCE=uuid.uuid4().hex
STATE={};LOCK=threading.Lock();ROOT=BASE
SERVER_STOPPING=threading.Event()
RECOVERY_OWNED=set()
INTERNAL_STATE_FIELDS={'cancel','thread','design_lock','process','configuration'}

def public_job(item):
    return {key:value for key,value in item.items() if key not in INTERNAL_STATE_FIELDS}

def saved_run_root(ident):
    if not isinstance(ident,str) or not re.fullmatch(r'[a-z0-9_-]{1,100}',ident):
        raise ValueError('Unknown saved run')
    root=(JOBS/ident).resolve()
    if root.parent!=JOBS.resolve() or not root.is_dir():raise ValueError('Unknown saved run')
    return root

def small_saved_json(path,limit=1024*1024):
    path=Path(path)
    if path.is_symlink():raise ValueError('Saved progress must not use symlinks')
    with path.open('rb') as handle:data=handle.read(limit+1)
    if len(data)>limit:raise ValueError('Saved progress record exceeds the safe read limit')
    record=json.loads(data)
    if not isinstance(record,dict):raise ValueError('Invalid saved progress record')
    return record

def completed_run(ident, include_result=False):
    """Only committed on-disk analyses may be opened after a server restart."""
    root=saved_run_root(ident)
    if (JOBS/ident).is_symlink() or (root/'result').is_symlink():raise ValueError('Unsafe saved result')
    status=small_saved_json(root/'analysis_status.json');job=small_saved_json(root/'job_status.json')
    archive=JOBS/(ident+'.zip')
    if (status.get('status')!='complete' or status.get('analysis_complete') is not True
        or status.get('job_id')!=ident or job.get('status')!='complete' or job.get('id')!=ident
        or not archive.is_file() or archive.is_symlink() or (JOBS/(ident+'.INCOMPLETE_EXPORT.json')).exists()):
        raise ValueError('Analysis is incomplete; no completed result is available.')
    result=small_saved_json(root/'result/design_result.json',limit=32*1024*1024)
    from run_summaries import summarize_completed_result
    summary=summarize_completed_result(result,ident)
    if (root/'result/run_summary.json').exists():
        saved=small_saved_json(root/'result/run_summary.json',limit=4*1024*1024)
        if saved.get('run_id')==ident and saved.get('status')=='complete':summary=saved
    value=dict(id=ident,name=summary['name'],status='complete',summary=dict(counts=summary['counts']))
    if include_result:value['result']=dict(result,title=summary['name'],run_name=summary['name'])
    return value


def completed_runs():
    rows=[]
    if JOBS.exists():
        for root in sorted(JOBS.iterdir(),key=lambda p:p.name,reverse=True):
            if not root.is_dir() or root.is_symlink():continue
            try:rows.append(completed_run(root.name))
            except (OSError,ValueError,TypeError,KeyError):continue
    return dict(runs=rows)


def result_with_coverage(result, result_dir):
    """Enrich a response/new export without rewriting historical result files.

    Presentation failures must never turn a completed scientific analysis into a
    failure. The dedicated coverage endpoint reports such failures explicitly.
    """
    value=dict(result)
    try:
        value['coverage_report']=coverage_report.build_report(result,result_dir)
        value.pop('coverage_report_error',None)
    except (OSError,ValueError,TypeError,KeyError,OverflowError) as exc:
        value.pop('coverage_report',None)
        value['coverage_report_error']='Saved coverage could not be displayed: '+str(exc)
    return value


def prune_interrupted_standard_batches(root):
    """Called with the global lease after the owned worker has stopped.

    Failed TSV/query files are never checkpoint evidence. Preserve command,
    return-code, stderr and stop records; only remove disposable search payloads.
    """
    folder=Path(root)/'result/blast_attempts'
    removed=[]
    if not folder.is_dir() or folder.is_symlink():return removed
    for parent,dirs,names in os.walk(folder,followlinks=False):
        dirs[:]=[d for d in dirs if not (Path(parent)/d).is_symlink()]
        for name in names:
            path=Path(parent)/name
            if not path.is_symlink() and (name.endswith('_candidates_blast.tsv') or name.endswith('_queries.fa')):
                size=path.stat().st_size;path.unlink();removed.append(dict(path=str(path.relative_to(root)),bytes=size))
    if removed:record_event_safely(Path(root)/'storage_events.jsonl',dict(event='DISCARDED_UNCOMMITTED_SEARCH_PAYLOADS',files=removed))
    return removed


def compact_memory_logs(root):
    """Keep peak observations, not thousands of repeated process snapshots."""
    summary={}
    paths=[]
    for name in ('memory_usage.jsonl','worker_memory_usage.jsonl'):
        path=Path(root)/name
        if not path.is_file() or path.is_symlink():continue
        peaks={};samples=0;stops=[]
        with path.open() as handle:
            for line in handle:
                try:sample=json.loads(line)
                except ValueError:continue
                samples+=1
                for key in ('process_family_rss_bytes','process_family_physical_footprint_bytes'):
                    if isinstance(sample.get(key),(int,float)):peaks[key]=max(peaks.get(key,0),sample[key])
                if sample.get('status') not in (None,'OK'):stops.append({k:sample[k] for k in ('timestamp','status','reason') if k in sample})
        summary[name]=dict(samples=samples,peaks=peaks,stops=stops[-20:])
        paths.append(path)
    if summary:
        atomic_json(Path(root)/'result/runtime_summary.json',dict(process_memory=summary,note='Peak measurements only; detailed process snapshots require Full audit output.'))
        for path in paths:path.unlink()
    return summary

def checkpoint_origin(root):
    """Resolve only a direct local run ID, never an arbitrary caller path."""
    root=Path(root).resolve()
    link=root/'checkpoint_link.json'
    origin=root
    if link.exists():
        record=small_saved_json(link)
        if record.get('schema')!='raserfish.checkpoint-link.v1':raise ValueError('Unsupported checkpoint link')
        origin=saved_run_root(record.get('origin_run_id'))
    directory=origin/'checkpoints'
    if directory.is_symlink() or directory.resolve().parent!=origin:
        raise ValueError('Invalid checkpoint directory')
    return origin,directory

def recovery_records():
    """Only explicitly enrolled new/resumed runs can restart automatically."""
    import automatic_recovery as recovery
    if not JOBS.exists():return []
    records=[]
    for path in sorted(JOBS.iterdir(),key=lambda p:p.name,reverse=True):
        if not path.is_dir() or path.is_symlink() or not (path/'automatic_recovery.json').exists():continue
        try:
            origin=saved_run_root(path.name);record=recovery.read_record(origin)
            if record:records.append((origin,record))
        except (OSError,ValueError,TypeError) as exc:
            record_event_safely(BASE/'logs/automatic_recovery_errors.jsonl',dict(event='INVALID_RECOVERY_RECORD',run=path.name,error=str(exc)))
    return records


def recovery_for_job(ident):
    import automatic_recovery as recovery
    try:
        origin,_=checkpoint_origin(saved_run_root(ident))
        return recovery.read_record(origin)
    except (OSError,ValueError):
        return next((record for _,record in recovery_records() if record['latest_run_id']==ident),None)


def job_with_recovery(ident):
    with LOCK:item=STATE.get(ident)
    recovery=recovery_for_job(ident)
    if item is None and recovery and recovery['state']=='complete' and recovery['latest_run_id']==ident:
        try:
            if restore_completed_recovery(saved_run_root(recovery['origin_run_id']),recovery):item=STATE.get(ident)
        except (OSError,ValueError):pass
    if item is not None:result=public_job(item)
    else:
        try:result=persisted_incomplete_job(ident)
        except ValueError:
            if not recovery:raise
            result=dict(id=ident,status='incomplete',stage='Interrupted - incomplete',analysis_complete=False,
                messages=['The previous server stopped; saved progress is being checked.'])
    if recovery:result['auto_recovery']=recovery
    return result


def recovery_active_design():
    """Read-only browser attachment; never launches a run."""
    with LOCK:
        item=next((j for j in STATE.values() if j['status']=='running'),None)
        if item:
            return dict(id=item['id'],config=copy.deepcopy(item.get('configuration')),
                        auto_recovery=recovery_for_job(item['id']),token=TOKEN)
    for origin,record in recovery_records():
        if record['state'] not in ('active','pending'):continue
        try:
            config,_=prepare_saved_resume(origin.name)
            return dict(id=record['latest_run_id'],config=config,auto_recovery=record,token=TOKEN)
        except (OSError,ValueError):continue
    return dict(id=None,config=None,auto_recovery=None,token=TOKEN)


def cancel_design_chain(ident,reason='Cancelled by user; automatic recovery is stopped.'):
    import automatic_recovery as recovery
    record=recovery_for_job(ident)
    target=record['latest_run_id'] if record else ident
    if record:
        stopped=recovery.stop(saved_run_root(record['origin_run_id']),reason)
        target=stopped['latest_run_id'] if stopped else target
    with LOCK:
        item=STATE.get(target)
        if item:
            item['user_cancelled']=True;item['cancel_reason']=reason;item['cancel'].set()
    return target


def stop_owned_recovery_chains(reason):
    """Intentional shutdown cannot cancel a retry claimed by another server."""
    import automatic_recovery as recovery
    with LOCK:
        owned={(j.get('checkpoint_origin',j.get('id')),j.get('id')) for j in STATE.values()}
    owned.update(RECOVERY_OWNED)
    for origin,record in recovery_records():
        if record['state'] not in ('active','pending') or (origin.name,record['latest_run_id']) not in owned:continue
        # A claim changes latest_run_id. A stale shutdown must not follow it.
        recovery.stop(origin,reason,expected=record)


def observe_recovery_cancellation(item,origin):
    import automatic_recovery as recovery
    record=recovery.read_record(origin)
    if record and record['latest_run_id']==item['id'] and record['state']=='stopped':
        item['user_cancelled']=True;item['cancel_reason']=record['reason'];item['cancel'].set()
    if item['cancel'].is_set():raise RuntimeError(item.get('cancel_reason','Cancelled by user. Analysis is incomplete.'))


def finish_recovery_attempt(job_id,exc=None):
    """Called after owned-worker cleanup, before releasing the global design lease."""
    import automatic_recovery as recovery
    from design_checkpoints import read_checkpoint_summary
    item=STATE[job_id];origin,_=checkpoint_origin(JOBS/job_id)
    record=recovery.read_record(origin)
    if not record:return
    if item.get('status')=='complete':
        recovery.stop(origin,'All requested analyses and final export finished.',job_id=job_id,complete=True);return
    details={}
    try:details=small_saved_json(JOBS/job_id/'worker_status.json')
    except (OSError,ValueError):pass
    cancelled=item.get('user_cancelled') or (item['cancel'].is_set() and not item.get('shutdown_interruption'))
    blocked=recovery.retry_block_reason(str(exc or item.get('error','Interrupted attempt')),
        exception_type=details.get('exception_type',type(exc).__name__ if exc else ''),
        memory_sample=details.get('memory_sample') or getattr(exc,'sample',None),cancelled=cancelled)
    summary=read_checkpoint_summary(origin/'checkpoints')
    if not summary or not summary.get('resumable',True):
        blocked=blocked or (summary or {}).get('resume_reason') or 'No verified-format checkpoint is available; automatic recovery stopped.'
    record=recovery.schedule(origin,job_id,str(exc or item.get('error','Interrupted attempt')),blocked=blocked)
    if record:item['auto_recovery']=record


def restore_completed_recovery(origin,record):
    """Require all export commit records; an interrupted archive is never complete."""
    root=saved_run_root(record['latest_run_id']);ident=root.name
    linked,_=checkpoint_origin(root)
    if linked!=origin.resolve():return False
    if (root/'checkpoint_link.json').exists():
        link=small_saved_json(root/'checkpoint_link.json')
        if link.get('attempt_run_id')!=ident:return False
    status=small_saved_json(root/'analysis_status.json');job=small_saved_json(root/'job_status.json')
    archive=JOBS/(ident+'.zip')
    if (status.get('status')!='complete' or not status.get('analysis_complete') or status.get('job_id')!=ident
        or job.get('id')!=ident or job.get('status')!='complete' or not archive.is_file() or archive.is_symlink()
        or (JOBS/(ident+'.INCOMPLETE_EXPORT.json')).exists()):return False
    with LOCK:
        if ident not in STATE:STATE[ident]=dict(job)
    import automatic_recovery as recovery
    recovery.stop(origin,'All requested analyses and final export finished.',job_id=ident,complete=True)
    return True


def recovery_tick():
    """Reconcile durable attempts and start at most one due retry with the global lease."""
    import automatic_recovery as recovery
    from design_checkpoints import read_checkpoint_summary
    from runtime_safety import DesignAlreadyRunning,acquire_design_lock
    if SERVER_STOPPING.is_set():return
    # The worker retains this lock after server loss until its own cleanup ends.
    if active_design_snapshot() is not None:return
    with LOCK:
        if any(j.get('status')=='running' for j in STATE.values()):return
    for origin,record in recovery_records():
        if record['state'] in ('stopped','complete'):continue
        if record['state']=='active':
            # Serialize reconciliation with fresh starts and other Studio instances.
            try:lease=acquire_design_lock(dict(app_pid=os.getpid(),app=str(BASE),maintenance='Automatic recovery reconciliation'))
            except DesignAlreadyRunning:return
            try:
                record=recovery.read_record(origin)
                if record['state']!='active':continue
                try:
                    if restore_completed_recovery(origin,record):continue
                except (OSError,ValueError):pass
                ident=record['latest_run_id'];details={};status={}
                try:details=small_saved_json(JOBS/ident/'worker_status.json')
                except (OSError,ValueError):pass
                try:status=small_saved_json(JOBS/ident/'analysis_status.json')
                except (OSError,ValueError):pass
                supervisor={}
                try:supervisor=small_saved_json(JOBS/ident/'supervisor_process.json')
                except (OSError,ValueError):pass
                error=details.get('error') or supervisor.get('stop_reason') or status.get('error') or 'The designer or its worker was interrupted before completion.'
                blocked=recovery.retry_block_reason(error,details.get('exception_type',''),details.get('memory_sample') or supervisor.get('memory_sample'))
                summary=read_checkpoint_summary(origin/'checkpoints')
                if not summary or not summary.get('resumable',True):
                    blocked=blocked or (summary or {}).get('resume_reason') or 'No verified-format checkpoint is available; automatic recovery stopped.'
                record=recovery.schedule(origin,ident,error,blocked=blocked)
                RECOVERY_OWNED.add((origin.name,record['latest_run_id']))
            finally:lease.release()
        if record['state']!='pending' or record['next_retry_at']>time.time():continue
        try:
            config,_=prepare_saved_resume(origin.name)
            start_design_job(config,resume_from=origin.name,automatic=True)
            return
        except DesignAlreadyRunning:return
        except Exception as exc:
            recovery.stop(origin,'Automatic recovery stopped: '+str(exc),expected=record)
            record_event_safely(BASE/'logs/automatic_recovery_errors.jsonl',dict(event='RESUME_REFUSED',run=origin.name,error=str(exc)))


def automatic_recovery_loop(stop):
    while not stop.is_set():
        try:recovery_tick()
        except Exception as exc:
            record_event_safely(BASE/'logs/automatic_recovery_errors.jsonl',dict(event='RECOVERY_CHECK_FAILED',error=str(exc)))
        stop.wait(1)


def active_design_snapshot():
    """Observe the shared design lock without overwriting its metadata."""
    import fcntl
    from runtime_safety import default_lock_path,_validate_lock_fd,_read_lock_metadata
    path=default_lock_path()
    try:fd=os.open(path,os.O_RDONLY|getattr(os,'O_NOFOLLOW',0))
    except FileNotFoundError:return None
    try:
        _validate_lock_fd(fd,path)
        try:fcntl.flock(fd,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:return _read_lock_metadata(fd) or {'active':True}
        return None
    finally:os.close(fd)

def saved_checkpoint_runs():
    """List progress hints only. The worker verifies evidence before reusing it."""
    from design_checkpoints import read_checkpoint_summary
    active=active_design_snapshot()
    with LOCK:live=next((public_job(j) for j in STATE.values() if j['status']=='running'),None)
    seen=set();runs=[]
    if JOBS.exists():
        for path in sorted(JOBS.iterdir(),key=lambda p:p.name,reverse=True):
            if not path.is_dir() or path.is_symlink():continue
            try:
                root=saved_run_root(path.name);origin,directory=checkpoint_origin(root)
                if origin.name in seen:continue
                summary=read_checkpoint_summary(directory)
                if summary is None:continue
                # The newest attempt represents this checkpoint chain, including
                # a completed attempt which must hide earlier incomplete ones.
                seen.add(origin.name)
                status=small_saved_json(root/'analysis_status.json')
                if status.get('status')=='complete':continue
                if not summary.get('resumable',True):continue
                runs.append(dict(summary,id=root.name,origin_run_id=origin.name,status='incomplete',
                    analysis_complete=False,auto_recovery=recovery_for_job(root.name),note='Saved progress; integrity and matching inputs are verified on resume.'))
            except (OSError,ValueError,KeyError):continue
    pending=next((r for _,r in recovery_records() if r['state'] in ('active','pending')),None)
    return dict(runs=runs,active_design=active is not None or live is not None or pending is not None,
                active_run=live or ({'id':active.get('job_id'),'status':'running'} if active else
                    ({'id':pending['latest_run_id'],'status':'retry_wait','auto_recovery':pending} if pending else None)))

def completed_checkpoint_attempt(origin):
    """A completed attempt closes its checkpoint chain, including older IDs."""
    for path in sorted(JOBS.iterdir(),key=lambda p:p.name,reverse=True):
        if not path.is_dir() or path.is_symlink():continue
        try:
            root=saved_run_root(path.name);candidate,_=checkpoint_origin(root)
            if candidate!=origin:continue
            status=small_saved_json(root/'analysis_status.json')
            if str(status.get('status','')).lower() in ('complete','completed'):return root.name
        except (OSError,ValueError,KeyError):continue
    return None

def prepare_saved_resume(ident):
    from compact_outputs import load_input_manifest,output_detail
    from design_checkpoints import read_checkpoint_summary
    root=saved_run_root(ident);origin,directory=checkpoint_origin(root)
    summary=read_checkpoint_summary(directory)
    if summary is None:raise ValueError('This run has no verified-format checkpoints. Older interrupted BLAST files cannot be resumed.')
    if not summary.get('resumable',True):raise ValueError(summary.get('resume_reason','This checkpoint cannot be resumed.'))
    status=small_saved_json(root/'analysis_status.json')
    if str(status.get('status','')).lower() in ('complete','completed'):raise ValueError('This run already completed; it cannot be resumed.')
    completed=completed_checkpoint_attempt(origin)
    if completed:raise ValueError('This checkpoint chain already completed in run '+completed+'; older attempts cannot be resumed.')
    manifest=origin/'result/run_config.json'
    if manifest.resolve().parent!=origin/'result' or manifest.is_symlink():raise ValueError('Invalid saved input manifest')
    if manifest.stat().st_size>16*1024*1024:raise ValueError('Saved input manifest exceeds the safe read limit')
    config=load_input_manifest(manifest,aggregate_decompressed_limit=128*1024*1024)
    if output_detail(config)!='standard':raise ValueError('Resume currently supports Standard output. Full audit requires a fresh run to preserve its complete raw evidence.')
    return config,origin.name

def start_design_job(config,resume_from=None,automatic=False):
    """Acquire the same cross-instance lease for fresh and resumed attempts."""
    from design_engine import Policy
    from runtime_safety import acquire_design_lock
    import automatic_recovery as recovery
    if SERVER_STOPPING.is_set():raise ValueError('The designer is shutting down.')
    config=copy.deepcopy(config)
    config['policy']=vars(Policy.from_dict(config.get('policy',{})))
    config['imagers']=validate_imagers(config['imagers'])
    if config.get('background_text'):config['background']=parse_fasta(config.pop('background_text'))
    else:config.pop('background_text',None)
    if config.get('transcript_text'):
        config['transcripts']=parse_fasta('\n'.join(line if line.startswith('>') else line.upper().replace('U','T') for line in config.pop('transcript_text').splitlines()))
    else:config.pop('transcript_text',None)
    with LOCK:
        if any(j['status']=='running' for j in STATE.values()):raise ValueError('A design is already running. Wait or cancel it before starting another.')
        if any(r['state'] in ('active','pending') and r['origin_run_id']!=resume_from for _,r in recovery_records()):
            raise ValueError('A design is awaiting automatic recovery. Cancel it before starting another design.')
        ident=datetime.now().strftime('%Y%m%d_%H%M%S')+'_'+uuid.uuid4().hex[:6]
        lease=acquire_design_lock(dict(job_id=ident,app_pid=os.getpid(),app=str(BASE)))
        try:
            if automatic:recovery.claim(saved_run_root(resume_from),ident)
            elif resume_from:
                from compact_outputs import output_detail
                if output_detail(config)=='standard':recovery.enroll(saved_run_root(resume_from),ident,reset=True)
        except BaseException:
            lease.release();raise
        STATE[ident]=dict(id=ident,status='running',stage='Verifying saved progress' if resume_from else 'Starting',
            messages=[],cancel=threading.Event(),design_lock=lease,checkpoint_origin=resume_from or ident,
            configuration=config,automatic_attempt=automatic)
    try:
        thread=threading.Thread(target=run_job,args=(ident,config),daemon=False);STATE[ident]['thread']=thread;thread.start()
    except BaseException:
        lease.release();STATE.pop(ident,None);raise
    return ident,config

def persisted_incomplete_job(ident):
    """Read bounded status evidence after restart; never revive a job or export.

    Only an explicit unfinished status is recovered. Old completion claims do
    not confer download access: the original export transaction is not replayed.
    """
    if not re.fullmatch(r'[a-z0-9_-]{1,100}',ident):raise ValueError('Unknown run')
    root=(JOBS/ident).resolve()
    if root.parent!=JOBS.resolve() or not root.is_dir():raise ValueError('Unknown run')
    records=[]
    for name in ('analysis_status.json','worker_status.json','job_status.json'):
        path=root/name
        if not path.exists():continue
        if path.resolve().parent!=root:raise ValueError('Invalid saved run status path')
        with path.open('rb') as handle:blob=handle.read(1024*1024+1)
        if len(blob)>1024*1024:raise ValueError('Saved run status exceeds the safe read limit')
        try:record=json.loads(blob)
        except (ValueError,UnicodeError) as exc:raise ValueError('Saved run status is unreadable; completion is unverified') from exc
        if not isinstance(record,dict) or not isinstance(record.get('status'),str):
            raise ValueError('Invalid saved run status; completion is unverified')
        if record.get('job_id',record.get('id',ident))!=ident:raise ValueError('Saved run status ID does not match')
        records.append(record)
    if not records:raise ValueError('Unknown run; no saved status was found')
    unfinished=[r for r in records if r['status'].lower() in ('incomplete','cancelled','failed','error','stopped')]
    if not unfinished:
        if any(r['status'].lower() in ('running','starting') for r in records):
            raise ValueError('This run belongs to another or previous server. Its completion is unverified; no analysis was resumed.')
        raise ValueError('This saved run is not loaded by this server. Open its saved results archive to inspect it; no analysis was resumed.')
    message=next((str(r['error'])[:4000] for r in unfinished if r.get('error')),
                 'This analysis stopped or did not finish. No completed result is available.')
    log_path=root/'run.log'
    if log_path.exists() and log_path.resolve().parent!=root:raise ValueError('Invalid saved run log path')
    messages=bounded_tail(log_path,65536).splitlines()[-120:]
    if not messages or message not in messages[-1]:messages.append('INCOMPLETE: '+message)
    return dict(id=ident,status='incomplete',stage='Stopped - incomplete',analysis_complete=False,
        error=message,messages=messages,recovered_from_disk=True,
        note='Saved status only. No work was resumed and interrupted search output is not completed evidence.')

def analysis_status(root, status, **details):
    payload=dict(status=status, analysis_complete=status=='complete',
                 updated_utc=datetime.now(timezone.utc).isoformat(), **details)
    atomic_json(root/'analysis_status.json',payload)
    return payload

def record_event_safely(path,event):
    """Failure to write diagnostics must never leave an owned worker running."""
    try:append_event(path,event)
    except Exception as exc:
        event=dict(event,diagnostic_write_error=str(exc),diagnostic_path=str(path))
        try:os.write(2,(json.dumps(event)+'\n').encode('utf-8')[-8192:])
        except OSError:pass
    return event

def stop_worker(process, root=None, reason='Owned worker cleanup', sample=None):
    """Record the reason before signalling only this job's owned process group."""
    if process is None or getattr(process,'_raserfish_group_retired',False):return []
    pgid=getattr(process,'_raserfish_owned_pgid',None)
    if pgid is None:
        try:pgid=os.getpgid(process.pid)
        except ProcessLookupError:
            process._raserfish_group_retired=True
            return []
        if pgid!=process.pid or pgid==os.getpgrp():
            raise RuntimeError('Refusing to stop a process group not owned by this design supervisor.')
    if pgid!=process.pid or pgid==os.getpgrp():
        raise RuntimeError('Refusing unsafe worker process-group termination.')
    actions=[]
    def signal_owned(sig):
        try:os.killpg(pgid,0)
        except ProcessLookupError:return False
        event=dict(action='SIGNAL_REQUESTED',actor_pid=os.getpid(),worker_pid=process.pid,
                   target_pgid=pgid,signal=signal.Signals(sig).name,reason=str(reason),memory_sample=sample)
        if root is not None:event=record_event_safely(Path(root)/'termination_events.jsonl',event)
        actions.append(dict(timestamp=utc_now(),**event))
        try:os.killpg(pgid,sig)
        except ProcessLookupError:pass
        return True
    signal_owned(signal.SIGTERM)
    try:process.wait(timeout=3)
    except subprocess.TimeoutExpired:pass
    # A TERM-resistant child may outlive the worker. Its owned group must end.
    signal_owned(signal.SIGKILL)
    try:process.wait(timeout=3)
    except subprocess.TimeoutExpired:
        raise RuntimeError('Owned design worker did not exit after SIGKILL; analysis incomplete.')
    # Never revisit this numeric PGID after packaging or other long work. It was
    # either absent or its owned members have now been sent SIGKILL and reaped.
    process._raserfish_group_retired=True
    if root is not None:
        finalize_interrupted_commands(Path(root)/'result',str(reason),actions)
    return actions

class SupervisorUnavailable(RuntimeError):
    pass

def worker_main(root, lock_fd):
    """Persist failures even when the supervisor pipe disappears."""
    from runtime_safety import inherited_design_lock,MemoryGuard
    from design_engine import design
    lease=inherited_design_lock(lock_fd)
    parent_pid=os.getppid();watcher_done=threading.Event()
    def report(message,**kw):
        event=dict(message=str(message),stage=kw.get('stage','Working'))
        try:print(json.dumps(event),flush=True)
        except (BrokenPipeError,OSError) as exc:
            failure='Design supervisor connection was lost; analysis is incomplete.'
            record_event_safely(root/'worker_events.jsonl',dict(action='SUPERVISOR_PIPE_LOST',reason=failure,
                original_event=event,error=str(exc),supervisor_pid=parent_pid))
            raise SupervisorUnavailable(failure) from exc
    def watch_parent():
        while not watcher_done.wait(1):
            if os.getppid()==parent_pid and parent_pid>1:continue
            failure='Design supervisor exited before its worker; analysis is incomplete.'
            event=dict(action='SUPERVISOR_LOST_STOP',actor_pid=os.getpid(),target_pgid=os.getpgrp(),
                       signal='SIGKILL',reason=failure)
            try:
                atomic_json(root/'worker_status.json',dict(status='incomplete',error=failure,updated_utc=utc_now()))
                analysis_status(root,'incomplete',error=failure)
                record_event_safely(root/'termination_events.jsonl',event)
                finalize_interrupted_commands(root/'result',failure,[event])
            except Exception as exc:
                record_event_safely(root/'termination_events.jsonl',dict(event,diagnostic_error=str(exc)))
            finally:
                # Only the explicitly created worker session is eligible. Even
                # a full/unwritable disk must not disable orphan termination.
                if os.getpgrp()==os.getpid():os.killpg(os.getpgrp(),signal.SIGKILL)
                os._exit(1)
    watcher=threading.Thread(target=watch_parent,daemon=True);watcher.start()
    try:
        config=json.load(sys.stdin)
        guard=MemoryGuard(root/'worker_memory_usage.jsonl',root_pid=os.getpid(),logger=lambda m:report(m,stage='Memory guard'))
        guard.checkpoint('Worker startup',force=True)
        _,checkpoint_dir=checkpoint_origin(root)
        result=design(config,root/'result',report,resource_guard=guard,checkpoint_dir=checkpoint_dir)
        if os.getppid()!=parent_pid:raise SupervisorUnavailable('Design supervisor exited; analysis is incomplete.')
        atomic_json(root/'worker_status.json',dict(status='complete',updated_utc=utc_now(),
            counts=result['counts'],result_status=result['status']))
        return 0
    except BaseException as exc:
        atomic_json(root/'worker_status.json',dict(status='incomplete',error=str(exc),
            exception_type=type(exc).__name__,memory_sample=getattr(exc,'sample',None),
            exception_notes=getattr(exc,'__notes__',[]),updated_utc=utc_now()))
        traceback.print_exc()
        try:report('INCOMPLETE: '+str(exc),stage='Stopped - incomplete')
        except SupervisorUnavailable:pass  # Never replace the original failure.
        return 1
    finally:
        watcher_done.set();lease.release()

def environment(root:Path):
    manifest=root/'references'/'GCF_000001405.40'/'reference_manifest.json'
    databases={};reference_info=None
    if manifest.exists():
        reference_info=json.loads(manifest.read_text())
        for label,v in reference_info.get('databases',{}).items():
            key='rna_human' if label.startswith('rna') else 'human'
            databases[key]=str(manifest.parent/v['database_prefix_relative'])
    return dict(python=sys.executable,python_version=sys.version.split()[0],qc_root=str(root),
       biopython=bool(importlib.util.find_spec('Bio')),vienna=bool(importlib.util.find_spec('RNA')),
       blastn=shutil.which('blastn'),pybigwig=bool(importlib.util.find_spec('pyBigWig')),
       databases=databases,reference_manifest=str(manifest) if manifest.exists() else None,
       assembly=reference_info.get('assembly_accession') if reference_info else None,
       note='Existing environment and databases are reused; nothing is downloaded during startup.')

def load_demo():
    path=BASE/'examples'/'compact_demo'/'design_result.json'
    if not path.exists():path=BASE/'examples'/'local_design'/'design_result.json'
    return json.loads(path.read_text()) if path.exists() else None

def bootstrap():
    with gzip.open(BASE/'inputs/errfi1_reference.json.gz','rt',encoding='utf-8') as f:
        references=json.load(f)
    return dict(token=TOKEN,environment=environment(ROOT),imagers=validate_imagers(json.loads((BASE/'imagers.json').read_text())),
      glossary=json.loads((BASE/'glossary.json').read_text()),references=references,
      default_config=copy.deepcopy(ERRFI1_DEFAULTS),demo=load_demo(),version='1.6.1',presets=PRESETS,parameter_notes=PARAMETER_NOTES)

def standalone(result:dict,preview=False,initial_step=None):
    boot=bootstrap();boot['token']='';boot['environment']={'offline':True};boot['demo']=result
    if not preview:
        boot['references']=[]
        boot.pop('default_config',None)
    boot['offline']=True;boot['preview']=preview
    if initial_step is not None:
        if initial_step not in range(1,7):raise ValueError('Initial step must be between 1 and 6.')
        boot['initial_step']=initial_step
    if not preview and result.get('imagers'):
        boot['imagers']=validate_imagers(result['imagers'])
    data=json.dumps(boot,separators=(',',':')).replace('<','\\u003c').replace('>','\\u003e').replace('&','\\u0026')
    html=(BASE/'web'/'index.html').read_text()
    html=html.replace('<link rel="stylesheet" href="/web/styles.css">','<style>'+(BASE/'web'/'styles.css').read_text()+'</style>')
    html=html.replace('<script src="/bootstrap.js"></script>','<script>window.STUDIO_BOOT='+data+';</script>')
    html=html.replace('<script src="/web/imager_binding.js"></script>','<script>'+(BASE/'web'/'imager_binding.js').read_text()+'</script>')
    for module in ('probe_coverage.js','coverage_panel.js'):
        tag='<script src="/web/'+module+'"></script>'
        if tag in html:html=html.replace(tag,'<script>'+(BASE/'web'/module).read_text()+'</script>')
    html=html.replace('<script src="/web/app.js"></script>','<script>'+(BASE/'web'/'app.js').read_text()+'</script>')
    return html

def quarantine_unfinished_export(job_id,reason):
    """Preserve a new failed export without publishing its prepared journals.

    Rename is constant-memory even for a large full-audit archive. A clear marker
    invalidates any completion journals prepared before the transaction finished.
    This never operates on old runs: callers pass only their newly created job.
    """
    final=JOBS/(job_id+'.zip');partial=JOBS/(job_id+'.incomplete.zip')
    source=final if final.exists() else partial
    if not source.exists():return None
    target=JOBS/(job_id+'.INCOMPLETE_EXPORT')
    source.replace(target)
    marker=dict(status='INCOMPLETE',analysis_complete=False,archive=str(target),reason=str(reason),
        note='Export finalization failed. Any completion journals prepared inside this unpublished archive are invalid. Do not use as completed results.',
        updated_utc=utc_now())
    atomic_json(JOBS/(job_id+'.INCOMPLETE_EXPORT.json'),marker)
    return target

def run_job(job_id:str,config:dict):
    from runtime_safety import MemoryGuard,MemoryMetricsUnavailable
    item=STATE[job_id];root=JOBS/job_id
    process=None;reader=None;lease=item['design_lock'];supervisor_record=None;reader_errors=[];archive_started=False;failure=None
    def log(message,**kw):
        with LOCK:
            item['messages'].append(message);item['stage']=kw.get('stage','Working')
            item['messages']=item['messages'][-2000:]
        with (root/'run.log').open('a') as f:f.write(datetime.now(timezone.utc).isoformat()+' '+message+'\n')
    try:
        root.mkdir(parents=True,exist_ok=True)
        origin=item.get('checkpoint_origin',job_id)
        saved_run_root(origin)
        from compact_outputs import output_detail
        if output_detail(config)=='standard':
            import automatic_recovery as recovery
            recovery.enroll(JOBS/origin,job_id)
            # The lease proves no previous attempt is still executing. Dispose
            # of native output left by a hard interruption before retrying.
            for old in JOBS.iterdir():
                if not old.is_dir() or old.is_symlink():continue
                try:linked,_=checkpoint_origin(old)
                except (OSError,ValueError):continue
                if linked==saved_run_root(origin):prune_interrupted_standard_batches(old)
        atomic_json(root/'checkpoint_link.json',dict(schema='raserfish.checkpoint-link.v1',
            origin_run_id=origin,attempt_run_id=job_id,
            note='Local resumable work only; a checkpoint is not a completed analysis or an order.'))
        analysis_status(root,'running',job_id=job_id,note='Unfinished analysis. Draft files are not a completed result.')
        guard=MemoryGuard(root/'memory_usage.jsonl',root_pid=os.getpid(),logger=lambda m:log(m,stage='Memory guard'))
        observe_recovery_cancellation(item,JOBS/origin)
        guard.check('Before starting design',starting=True)
        worker_env=os.environ.copy();worker_env.update(RASERFISH_DESIGN_WORKER='1',
            RASERFISH_DESIGN_LOCK_FD=str(lease.fd),PYTHONDONTWRITEBYTECODE='1')
        worker_command=[sys.executable,str(BASE/'studio.py'),'--worker-output',str(root),'--lock-fd',str(lease.fd)]
        supervisor_record=command_record(worker_command)
        supervisor_record.update(stdout_path=str(root/'worker_stdout.log'),stderr_path=str(root/'worker_stderr.log'))
        atomic_json(root/'supervisor_process.json',supervisor_record)
        with (root/'worker_stderr.log').open('w') as errors:
            process=subprocess.Popen(worker_command,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=errors,
                text=True,env=worker_env,start_new_session=True,pass_fds=(lease.fd,))
            process._raserfish_owned_pgid=process.pid
            supervisor_record.update(status='RUNNING',pid=process.pid,pgid=process.pid)
            atomic_json(root/'supervisor_process.json',supervisor_record)
            item['process']=process
            lease.update_metadata(dict(job_id=job_id,worker_pid=process.pid,app_pid=os.getpid(),app=str(BASE)))
            def receive():
                try:
                    with (root/'worker_stdout.log').open('w') as output:
                        while True:
                            line=process.stdout.readline(65536)
                            if not line:break
                            output.write(line);output.flush()
                            try:
                                event=json.loads(line)
                                log(str(event['message']),stage=event.get('stage','Working'))
                            except (ValueError,KeyError,TypeError):log(line.rstrip()[:4000],stage='Worker output')
                except Exception as exc:
                    reader_errors.append(str(exc))
            reader=threading.Thread(target=receive,daemon=True);reader.start()
            process.stdin.write(json.dumps(config));process.stdin.close()
            guard.root_pid=process.pid
            while process.poll() is None:
                observe_recovery_cancellation(item,JOBS/origin)
                if reader_errors:raise RuntimeError('Cannot preserve worker output: '+reader_errors[0])
                if item['cancel'].is_set():raise RuntimeError(item.get('cancel_reason','Cancelled by user. Analysis is incomplete; partial files retained.'))
                try:guard.checkpoint(item.get('stage','Working'))
                except MemoryMetricsUnavailable:
                    # Normal worker exit can race the process-table snapshot.
                    # Final export is checked again against the living supervisor.
                    if process.poll() is not None:break
                    raise
                time.sleep(.25)
            reader.join(timeout=3)
            cleanup_actions=stop_worker(process,root,reason='Owned child cleanup immediately after worker exit')
            if cleanup_actions and process.returncode==0:
                raise RuntimeError('Design worker exited before its owned children finished; analysis is incomplete.')
            supervisor_record.update(returncode_fields(process.returncode))
            supervisor_record.update(ended_utc=utc_now(),status='COMPLETE' if process.returncode==0 else 'INCOMPLETE',
                analysis_complete=process.returncode==0,stdout_tail=bounded_tail(root/'worker_stdout.log'),
                stderr_tail=bounded_tail(root/'worker_stderr.log'))
            atomic_json(root/'supervisor_process.json',supervisor_record)
            if reader_errors:raise RuntimeError('Cannot preserve worker output: '+reader_errors[0])
            if process.returncode:
                record=root/'worker_status.json'
                detail=json.loads(record.read_text()).get('error') if record.exists() else None
                outcome=failure_message(supervisor_record,root/'supervisor_process.json')
                raise RuntimeError((detail+'\n' if detail else '')+outcome)
        record=json.loads((root/'worker_status.json').read_text())
        if record.get('status')!='complete':raise RuntimeError('Worker did not finish all requested analyses.')
        observe_recovery_cancellation(item,JOBS/origin)
        if item['cancel'].is_set():raise RuntimeError('Cancelled before export. Analysis is incomplete.')
        guard.root_pid=os.getpid();guard.checkpoint('Final export',force=True)
        result=json.loads((root/'result'/'design_result.json').read_text())
        from run_summaries import write_run_summary
        summary=write_run_summary(result,root/'result',job_id,config)
        result.update(title=summary['name'],run_name=summary['name'])
        result=result_with_coverage(result,root/'result')
        atomic_json(root/'result/design_result.json',result)
        if output_detail(config)=='standard':compact_memory_logs(root)
        (root/'result'/'viewer.html').write_text(standalone(result),encoding='utf-8')
        partial=JOBS/(job_id+'.incomplete.zip')
        archive_started=True
        with zipfile.ZipFile(partial,'w',zipfile.ZIP_DEFLATED) as archive:
            for file in result_archive_files(root):
                if file.is_file() and file.name not in ('job_status.json','analysis_status.json'):
                    observe_recovery_cancellation(item,JOBS/origin)
                    guard.checkpoint('Packaging results')
                    if item['cancel'].is_set():raise RuntimeError('Cancelled during export. Analysis is incomplete.')
                    archive.write(file,file.relative_to(JOBS).as_posix())
            guard.checkpoint('Export finished',force=True)
            status=public_job(item)
            status.update(status='complete',counts=result['counts'],warnings=result['warnings'],stage='Complete')
            completed=dict(status='complete',analysis_complete=True,updated_utc=datetime.now(timezone.utc).isoformat(),
                job_id=job_id,result_status=result['status'],cross_primary_screening=result.get('cross_primary_screening'),
                note='Requested analyses finished. Omitted analyses remain NOT RUN; all orders are DRAFT.')
            archive.writestr(job_id+'/job_status.json',json.dumps(status,indent=2)+'\n')
            archive.writestr(job_id+'/analysis_status.json',json.dumps(completed,indent=2)+'\n')
        partial.replace(JOBS/(job_id+'.zip'))
        atomic_json(root/'job_status.json',status)
        atomic_json(root/'analysis_status.json',completed)
        with LOCK:item.update(status='complete',counts=result['counts'],warnings=result['warnings'],stage='Complete')
    except Exception as e:
        failure=e
        with LOCK:item.update(status='incomplete',error=str(e),stage='Stopped - incomplete')
        if archive_started:
            try:quarantine_unfinished_export(job_id,str(e))
            except Exception as export_error:
                record_event_safely(root/'termination_events.jsonl',dict(action='EXPORT_FINALIZATION_FAILED',
                    reason=str(e),quarantine_error=str(export_error),analysis_complete=False))
        try:actions=stop_worker(process,root,reason=str(e),sample=getattr(e,'sample',None))
        except Exception as stop_error:
            actions=[]
            e=RuntimeError(str(e)+'; owned-worker cleanup failed: '+str(stop_error))
            failure=e
            record_event_safely(root/'termination_events.jsonl',dict(action='STOP_FAILED',reason=str(e)))
        if supervisor_record is not None:
            supervisor_record.update(status='INCOMPLETE',analysis_complete=False,ended_utc=utc_now(),
                stop_reason=str(e),memory_sample=getattr(e,'sample',None),termination_actions=actions,
                stdout_tail=bounded_tail(root/'worker_stdout.log'),stderr_tail=bounded_tail(root/'worker_stderr.log'))
            supervisor_record.update(returncode_fields(process.returncode if process is not None else None))
            if process is None:supervisor_record['spawn_error']=str(e)
            atomic_json(root/'supervisor_process.json',supervisor_record)
        with LOCK:item.update(status='incomplete',error=str(e),stage='Stopped - incomplete')
        analysis_status(root,'incomplete',job_id=job_id,error=str(e),note='Stopped or unfinished analyses are not PASS. Partial outputs remain DRAFT and incomplete.')
        atomic_json(root/'job_status.json',public_job(item))
        (root/'ERROR.txt').write_text(traceback.format_exc())
        log('STOPPED - INCOMPLETE: '+str(e),stage='Stopped - incomplete')
    finally:
        try:
            if process is not None:stop_worker(process,root,reason='Final owned worker cleanup')
            if reader is not None:reader.join(timeout=3)
            if config.get('output_detail','standard')=='standard':
                try:prune_interrupted_standard_batches(root)
                except OSError as exc:record_event_safely(root/'storage_events.jsonl',dict(event='PAYLOAD_CLEANUP_FAILED',error=str(exc)))
            if item.get('status')=='complete' and config.get('output_detail','standard')=='standard':
                attempts=root/'result/blast_attempts'
                if attempts.is_dir() and not attempts.is_symlink():
                    try:shutil.rmtree(attempts)
                    except OSError as exc:record_event_safely(root/'storage_events.jsonl',dict(event='COMPLETED_DIAGNOSTIC_CLEANUP_FAILED',error=str(exc)))
            try:finish_recovery_attempt(job_id,failure)
            except Exception as exc:
                record_event_safely(root/'automatic_recovery_errors.jsonl',dict(event='RECOVERY_FINALIZE_FAILED',error=str(exc)))
        finally:lease.release()

def result_archive_files(root):
    """Standard exports contain result scores/inputs and final status journals."""
    manifest=Path(root)/'result/run_config.json'
    detail=small_saved_json(manifest).get('output_detail','standard') if manifest.exists() else 'standard'
    for folder,dirs,files in os.walk(root,followlinks=False):
        dirs[:]=sorted(name for name in dirs if name!='checkpoints' and (detail!='standard' or name!='blast_attempts') and not (Path(folder)/name).is_symlink())
        for name in sorted(files):
            file=Path(folder)/name
            if detail=='standard' and Path(folder)==Path(root) and name not in ('job_status.json','analysis_status.json'):continue
            if not file.is_symlink():yield file

def normalize_config_inputs(config):
    config=copy.deepcopy(config)
    for text_key,records_key in (('background_text','background'),('transcript_text','transcripts')):
        text=config.pop(text_key,'')
        if text:
            if records_key=='transcripts':
                text='\n'.join(line if line.startswith('>') else line.upper().replace('U','T') for line in text.splitlines())
            config[records_key]=parse_fasta(text)
    return config

def export_settings(config):
    """Portable settings: one manifest and compressed input assets, no database copy."""
    from compact_outputs import write_input_manifest
    from design_engine import Policy
    config=normalize_config_inputs(config)
    config['policy']=vars(Policy.from_dict(config.get('policy',{})))
    config['imagers']=validate_imagers(config['imagers'])
    bins,_=make_bins(config['references'],int(config.get('bin_width',20000)),config.get('tail_policy','strict'),config.get('assignments'),config.get('bounds'))
    with tempfile.TemporaryDirectory(prefix='raserfish-settings-') as temp:
        root=Path(temp);write_input_manifest(config,bins,root,application_version='1.6.1')
        data=io.BytesIO()
        with zipfile.ZipFile(data,'w',zipfile.ZIP_DEFLATED) as archive:
            for file in sorted(root.rglob('*')):
                if file.is_file():archive.write(file,file.relative_to(root).as_posix())
        return data.getvalue()

def import_settings(data):
    from compact_outputs import load_input_manifest
    from design_engine import Policy
    if len(data)>25_000_000:raise ValueError('Settings archive exceeds 25 MB.')
    with tempfile.TemporaryDirectory(prefix='raserfish-settings-') as temp:
        root=Path(temp)
        with zipfile.ZipFile(io.BytesIO(data)) as archive:
            entries=archive.infolist()
            if len(entries)>100 or sum(x.file_size for x in entries)>50_000_000:
                raise ValueError('Settings archive has too many files or exceeds the expanded-size limit.')
            seen=set()
            for info in entries:
                path=Path(info.filename)
                if path.is_absolute() or '..' in path.parts or '\\' in info.filename or info.filename in seen:
                    raise ValueError('Unsafe or duplicate settings archive path.')
                seen.add(info.filename)
                if info.is_dir():continue
                dest=root/path
                dest.parent.mkdir(parents=True,exist_ok=True)
                dest.write_bytes(archive.read(info))
        config=load_input_manifest(root/'run_config.json',aggregate_decompressed_limit=50*1024*1024)
    Policy.from_dict(config.get('policy',{}))
    config['imagers']=validate_imagers(config['imagers'])
    return config

class Handler(BaseHTTPRequestHandler):
    server_version='RASERFISH-Design-Studio/1.6.1'
    def log_message(self,fmt,*args):
        if '/api/job?' not in str(args):super().log_message(fmt,*args)
    def safe_host(self):
        return self.headers.get('Host','').split(':')[0] in ('127.0.0.1','localhost')
    def send(self,data,status=200,kind='application/json',filename=None):
        blob=(json.dumps(data,allow_nan=False).encode() if kind=='application/json' else data.encode() if isinstance(data,str) else data)
        self.send_response(status);self.send_header('Content-Type',kind);self.send_header('Content-Length',str(len(blob)))
        self.send_header('X-Content-Type-Options','nosniff');self.send_header('Cache-Control','no-store')
        self.send_header('Content-Security-Policy',"default-src 'self'; script-src 'self' 'unsafe-inline'; style-src 'self' 'unsafe-inline'; img-src 'self' data:; connect-src 'self'; frame-ancestors 'none'")
        if filename:self.send_header('Content-Disposition',f'attachment; filename="{filename}"')
        self.end_headers();self.wfile.write(blob)
    def send_file(self,path,kind='application/octet-stream',filename=None):
        """Stream large exports using one 1 MiB buffer, including full-audit ZIPs."""
        with Path(path).open('rb') as handle:
            length=os.fstat(handle.fileno()).st_size
            self.send_response(200);self.send_header('Content-Type',kind)
            self.send_header('Content-Length',str(length))
            self.send_header('X-Content-Type-Options','nosniff');self.send_header('Cache-Control','no-store')
            self.send_header('Content-Security-Policy',"default-src 'self'; script-src 'self' 'unsafe-inline'; style-src 'self' 'unsafe-inline'; img-src 'self' data:; connect-src 'self'; frame-ancestors 'none'")
            if filename:self.send_header('Content-Disposition',f'attachment; filename="{filename}"')
            self.end_headers()
            try:
                while True:
                    chunk=handle.read(1024*1024)
                    if not chunk:break
                    self.wfile.write(chunk)
            except (BrokenPipeError,ConnectionResetError):return

    def do_GET(self):
        if not self.safe_host():return self.send({'error':'Invalid host'},403)
        u=urllib.parse.urlparse(self.path);q=urllib.parse.parse_qs(u.query)
        try:
            if u.path=='/api/server':
                return self.send(dict(version='1.6.1',pid=os.getpid(),instance_id=SERVER_INSTANCE,
                    app_path=str(BASE),port=self.server.server_port))
            if u.path=='/api/runs':return self.send(completed_runs())
            if u.path=='/api/coverage-sources':
                return self.send(coverage_report.list_sources(BASE,JOBS,completed_run))
            if u.path=='/api/coverage':
                return self.send(coverage_report.load_source(BASE,JOBS,completed_run,
                    q.get('kind',[''])[0],q.get('id',[''])[0]))
            if u.path=='/api/checkpoints':return self.send(saved_checkpoint_runs())
            if u.path=='/api/active-design':return self.send(recovery_active_design())
            if u.path=='/':return self.send((BASE/'web/index.html').read_text(),kind='text/html; charset=utf-8')
            if u.path=='/bootstrap.js':return self.send('window.STUDIO_BOOT='+json.dumps(bootstrap()).replace('<','\\u003c')+';',kind='application/javascript')
            if u.path.startswith('/web/'):
                path=(BASE/u.path.lstrip('/')).resolve()
                if BASE/'web' not in path.parents:raise ValueError('Invalid path')
                return self.send(path.read_bytes(),kind=mimetypes.guess_type(path)[0] or 'application/octet-stream')
            if u.path=='/api/job':
                ident=q.get('id',[''])[0]
                return self.send(job_with_recovery(ident))
            if u.path in ('/api/result','/download'):
                ident=q.get('id',[''])[0]
                if ident in STATE and STATE[ident].get('status')!='complete':raise ValueError('Analysis is incomplete; no completed result is available.')
                name=q.get('file',['results.zip'])[0]
                needs_result=u.path=='/api/result' or name=='viewer.html'
                saved=completed_run(ident,include_result=needs_result)
                root=JOBS/ident
                if needs_result:saved['result']=result_with_coverage(saved['result'],root/'result')
                if u.path=='/api/result':return self.send(saved['result'])
                if name=='viewer.html':
                    return self.send(standalone(saved['result']),kind='text/html; charset=utf-8',
                        filename=re.sub(r'[^A-Za-z0-9._-]+','_',saved['name'])[:150]+'_'+ident[:8]+'_viewer.html')
                allowed={'results.zip':JOBS/(ident+'.zip'),'viewer.html':root/'result/viewer.html',
                   'selected_primaries.csv':root/'result/selected_primaries.csv','region_summary.csv':root/'result/region_summary.csv',
                   'candidate_decisions.csv.gz':root/'result/debug/candidate_decisions.csv.gz','DRAFT_IDT_oPools.csv':root/'result/DRAFT_IDT_oPools.csv',
                   'primary_scores.csv':root/'result/primary_scores.csv','run_summary.json':root/'result/run_summary.json',
                   'imager_inventory.csv':root/'result/imager_inventory.csv','readout_assignments.csv':root/'result/readout_assignments.csv'}
                if name not in allowed:raise ValueError('Unknown export')
                return self.send_file(allowed[name],kind='application/octet-stream',filename=re.sub(r'[^A-Za-z0-9._-]+','_',saved['name'])[:150]+'_'+ident[:8]+'_'+name)
            return self.send({'error':'Not found'},404)
        except (ValueError,OSError,RuntimeError,TypeError,KeyError,OverflowError) as e:return self.send({'error':str(e)},400)
    def do_POST(self):
        if not self.safe_host() or self.headers.get('X-Studio-Token')!=TOKEN:return self.send({'error':'Invalid local session token; reload the page.'},403)
        origin=self.headers.get('Origin')
        if origin and origin!=f'http://{self.headers.get("Host")}':return self.send({'error':'Cross-origin request blocked'},403)
        try:
            n=int(self.headers.get('Content-Length','0'))
            if not 0<n<=100_000_000:raise ValueError('Empty or oversized request; do not upload whole genomes.')
            body=json.loads(self.rfile.read(n));path=urllib.parse.urlparse(self.path).path
            if path=='/api/settings/export':
                return self.send(export_settings(body['config']),kind='application/zip',filename='RASERFISH_design_settings.zip')
            if path=='/api/settings/import':
                return self.send({'config':import_settings(base64.b64decode(body['base64'],validate=True))})
            if path=='/api/references':
                mode=body['mode']
                if mode=='bundled':refs=parse_fasta((BASE/'inputs/target_sequences.fasta').read_text())
                elif mode=='fasta':refs=parse_fasta(body['text'])
                elif mode=='annotated':refs=parse_annotated(base64.b64decode(body['base64'],validate=True),body['filename'])
                elif mode=='genomic':refs=[annot.get_genomic_reference(body['locus'])]
                else:raise ValueError('Unknown input mode')
                return self.send({'references':refs})
            if path=='/api/imagers/resolve':
                return self.send({'imagers':validate_imagers(body['imagers'])})
            if path=='/api/bins':
                bins,tails=make_bins(body['references'],int(body['width']),body.get('tail_policy','strict'),body.get('assignments'),body.get('bounds'))
                return self.send({'bins':bins,'tails':tails})
            if path=='/api/genes':return self.send(annot.fetch_genes(body['reference']))
            if path=='/api/transcription/preview':
                from rna_avoidance import prepare_transcription
                _,rna,models=prepare_transcription(body['references'],body.get('transcription',{}),body.get('enabled',True))
                return self.send({'models':models,'derived_rna_count':len(rna)})
            if path=='/api/h3k/discover':return self.send(annot.discover_h3k27ac(body.get('cell',''),body.get('fallback',True)))
            if path=='/api/h3k/load':return self.send(annot.load_h3k27ac(body['reference'],body['choice'],body.get('requested_cell','')))
            if path=='/api/features/import':
                choice=dict(kind=body.get('kind','annotation'),label=body.get('label','User annotation'),accession=body.get('filename','user BED'),cell=body.get('cell',''))
                return self.send({'features':annot.parse_bed_track(body['text'],body['reference'],choice,local=body.get('local',True))})
            if path=='/api/design':
                ident,_=start_design_job(body['config'])
                return self.send({'id':ident})
            if path=='/api/resume':
                config,origin=prepare_saved_resume(body.get('id'))
                ident,config=start_design_job(config,resume_from=origin)
                return self.send({'id':ident,'config':config})
            if path=='/api/cancel':
                cancel_design_chain(body['id'])
                return self.send({'message':'Cancellation requested. The owned design process and its BLAST child will stop; analysis remains incomplete.'})
            return self.send({'error':'Not found'},404)
        except Exception as e:return self.send({'error':str(e)},400)

class ServerStop(BaseException):
    def __init__(self,signum):self.signum=signum

def serve_studio(port,no_open=False):
    """Journal catchable server exits and let each owning supervisor stop work.

    SIGKILL and OS disappearance cannot run Python cleanup; the worker's
    independent parent watcher still leaves those runs incomplete.
    """
    server=None;exit_code=0;reason='Server loop returned normally';old_handlers={}
    SERVER_STOPPING.clear();recovery_stop=threading.Event();recovery_thread=None;intentional_stop=False
    lifecycle=BASE/'logs'/'server_lifecycle.jsonl'
    context=dict(pid=os.getpid(),parent_pid=os.getppid(),pgid=os.getpgrp(),
        session_id=os.getsid(0),instance_id=SERVER_INSTANCE,app_path=str(BASE),
        port=port,argv=list(sys.argv))
    def journal(event,**details):
        record_event_safely(lifecycle,dict(context,event=event,**details))
    def caught_signal(signum,frame):
        # A second Ctrl-C/hangup must not interrupt owned-worker cleanup.
        for candidate in old_handlers:signal.signal(candidate,signal.SIG_IGN)
        journal('SIGNAL_RECEIVED',signal=signal.Signals(signum).name)
        raise ServerStop(signum)
    journal('SERVER_STARTING')
    try:
        for candidate in (signal.SIGTERM,signal.SIGHUP,signal.SIGINT):
            old_handlers[candidate]=signal.signal(candidate,caught_signal)
        server=ThreadingHTTPServer(('127.0.0.1',port),Handler)
        context['port']=server.server_port
        address=f'http://127.0.0.1:{server.server_port}'
        journal('SERVER_STARTED',address=address)
        print(f'RASER-FISH Design Studio 1.6.1: {address}\nKeep this Terminal open. Stop with Ctrl-C.',flush=True)
        recovery_thread=threading.Thread(target=automatic_recovery_loop,args=(recovery_stop,),daemon=True)
        recovery_thread.start()
        if not no_open:webbrowser.open(address)
        server.serve_forever()
    except ServerStop as exc:
        intentional_stop=True
        exit_code=128+exc.signum;reason='Server received '+signal.Signals(exc.signum).name
    except KeyboardInterrupt:
        intentional_stop=True
        exit_code=130;reason='Server interrupted by keyboard'
        journal('KEYBOARD_INTERRUPT')
    except BaseException as exc:
        exit_code=1;reason='Server exception: '+str(exc)
        journal('UNHANDLED_EXCEPTION',exception_type=type(exc).__name__,error=str(exc),
            traceback=traceback.format_exc()[-16384:])
        raise
    finally:
        SERVER_STOPPING.set();recovery_stop.set()
        if recovery_thread:recovery_thread.join(timeout=5)
        if intentional_stop:
            try:stop_owned_recovery_chains('Designer intentionally closed; automatic recovery is stopped.')
            except Exception as exc:
                journal('RECOVERY_SHUTDOWN_RECORD_FAILED',error=str(exc))
        journal('SERVER_STOPPING',reason=reason,exit_code=exit_code)
        with LOCK:items=list(STATE.values())
        for item in items:
            if item.get('status')=='running':
                item['shutdown_interruption']=not intentional_stop
                item['cancel_reason']=reason+'. Analysis is incomplete; partial files retained.'
                item['cancel'].set()
        if server is not None:server.server_close()
        for item in items:
            if item.get('thread'):item['thread'].join(timeout=10)
        # Normally run_job already reaped its child. Use only its explicit owned
        # Popen handle if an unresponsive supervisor exceeded the join window.
        for item in items:
            thread=item.get('thread')
            if thread and thread.is_alive() and item.get('process') is not None:
                try:stop_worker(item['process'],JOBS/item['id'],reason=reason)
                except Exception as exc:journal('WORKER_CLEANUP_ERROR',job_id=item.get('id'),error=str(exc))
                thread.join(timeout=3)
        pending=[item.get('id') for item in items if item.get('thread') and item['thread'].is_alive()]
        journal('SERVER_STOPPED',reason=reason,exit_code=exit_code,unfinished_supervisor_threads=pending)
        for candidate,old in old_handlers.items():signal.signal(candidate,old)
    return exit_code

def main():
    global ROOT
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--qc-root',type=Path,default=BASE)
    parser.add_argument('--port',type=int,default=0);parser.add_argument('--no-open',action='store_true')
    parser.add_argument('--worker-output',type=Path,help=argparse.SUPPRESS)
    parser.add_argument('--lock-fd',type=int,help=argparse.SUPPRESS)
    parser.add_argument('--build-preview',action='store_true');args=parser.parse_args();ROOT=args.qc_root.expanduser().resolve()
    if args.worker_output is not None:
        if args.lock_fd is None or os.environ.get('RASERFISH_DESIGN_WORKER')!='1':
            raise SystemExit('Design workers require the application single-design lock.')
        raise SystemExit(worker_main(args.worker_output.resolve(),args.lock_fd))
    if args.build_preview:
        result=load_demo()
        if not result:raise SystemExit('Generate examples/local_design first.')
        path=BASE/'OPEN_STUDIO_PREVIEW.html';path.write_text(standalone(result,True),encoding='utf-8');print(path);return
    JOBS.mkdir(exist_ok=True)
    return serve_studio(args.port,args.no_open)
if __name__=='__main__':raise SystemExit(main())
