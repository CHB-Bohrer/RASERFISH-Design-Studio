"""Durable checkpoint identity; saved partial work is never a completed design.

The score store validates its own transaction receipts. This module binds those
receipts to the exact inputs, candidate order, scientific implementation, tools,
and local database index metadata. Genome volumes are never read into memory or
hashed here. Metadata matching detects replacement, not arbitrary undetectable
in-place modification of database files; databases must remain immutable.
"""
from __future__ import annotations

import dataclasses
import datetime as dt
import hashlib
import importlib
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tempfile


SCHEMA = 'raserfish.design-checkpoints.v2'
BASE = Path(__file__).resolve().parent
MAX_METADATA_BYTES = 4 * 1024 * 1024
SOURCE_FILES = (
    'design_engine.py', 'specificity_store.py', 'design_checkpoints.py', 'blast_stream.py', 'process_diagnostics.py',
    'sequence_tools.py', 'rna_avoidance.py', 'strand_guard.py',
    'engine/qc_core.py', 'engine/oligominer_thermo_qc.py',
    'engine/probe_binding_viewer.py', 'design_presets.py', 'requirements.txt',
)
CANDIDATE_FIELDS = (
    'id', 'arm', 'full', 'reference', 'bin_id', 'prefix', 'suffix', 'sites',
    'start', 'end', 'strand', 'arm_length', 'tm', 'gc',
)


def _utc_now():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def _json_value(value):
    if dataclasses.is_dataclass(value):
        value = dataclasses.asdict(value)
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        if any(not isinstance(key, str) for key in value):
            raise ValueError('Checkpoint identity requires string dictionary keys.')
        return {key: _json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_value(item) for item in value]
    return value


def _digest(value):
    digest = hashlib.sha256()
    encoder = json.JSONEncoder(sort_keys=True, separators=(',', ':'),
                               ensure_ascii=False, allow_nan=False)
    for part in encoder.iterencode(_json_value(value)):
        digest.update(part.encode('utf-8'))
    return digest.hexdigest()


def _file_identity(path, *, hash_contents=False):
    path = Path(path).expanduser().resolve(strict=True)
    if not path.is_file():
        raise ValueError('Checkpoint identity expected a file: ' + str(path))
    before = path.stat()
    result = dict(path=str(path), bytes=before.st_size,
                  mtime_ns=before.st_mtime_ns, ctime_ns=before.st_ctime_ns,
                  inode=before.st_ino, device=before.st_dev)
    if hash_contents:
        digest = hashlib.sha256()
        with path.open('rb') as stream:
            for chunk in iter(lambda: stream.read(1024 * 1024), b''):
                digest.update(chunk)
        after = path.stat()
        if (before.st_size, before.st_mtime_ns, before.st_ctime_ns, before.st_ino,
            before.st_dev) != (after.st_size, after.st_mtime_ns, after.st_ctime_ns,
                              after.st_ino, after.st_dev):
            raise ValueError('A tool or source file changed while its checkpoint identity was read: ' + str(path))
        result['sha256'] = digest.hexdigest()
    return result


def database_identity(database):
    """Include local alias targets and volume metadata without reading genomes."""
    root = Path(database).expanduser().resolve()
    files, visited = {}, set()

    def visit(prefix):
        prefix = prefix.resolve()
        if str(prefix) in visited:
            return
        visited.add(str(prefix))
        matches = sorted(item for item in prefix.parent.glob(prefix.name + '.*') if item.is_file())
        if not matches:
            raise ValueError('No local BLAST index files found for checkpoint identity: ' + str(prefix))
        for item in matches:
            identity = _file_identity(item)
            files[identity['path']] = identity
            if item.suffix == '.nal':
                if item.stat().st_size > MAX_METADATA_BYTES:
                    raise ValueError('BLAST alias metadata is unexpectedly large: ' + str(item))
                text = item.read_text().replace('\\\n', ' ')
                for line in text.splitlines():
                    fields = shlex.split(line, comments=True)
                    if fields and fields[0].upper() == 'DBLIST':
                        for name in fields[1:]:
                            target = Path(name).expanduser()
                            visit(target if target.is_absolute() else item.parent / target)

    visit(root)
    indices = [files[key] for key in sorted(files)]
    return dict(path=str(root), index_files=indices, metadata_sha256=_digest(indices),
                validation='LOCAL_INDEX_METADATA; genome volumes are not content-hashed')


def _module_identity(name):
    module = importlib.import_module(name)
    path = getattr(module, '__file__', None)
    if not path:
        raise ValueError('Cannot identify QC dependency ' + name)
    return dict(version=str(getattr(module, '__version__', 'not exposed')),
                file=_file_identity(path, hash_contents=True))


def scientific_engine_identity(policy):
    policy = _json_value(policy)
    sources = {name: _file_identity(BASE / name, hash_contents=True)['sha256']
               for name in SOURCE_FILES}
    dependencies = {name: _module_identity(name) for name in (
        'Bio', 'Bio.Seq', 'Bio.SeqUtils', 'Bio.SeqUtils.MeltingTemp',
        'numpy', 'numpy._core._multiarray_umath',
    )}
    tools = {}
    if policy.get('mode', 'full') == 'full':
        executable = shutil.which('blastn')
        if not executable:
            raise ValueError('BLAST executable is unavailable; checkpoint identity cannot be verified.')
        binary = _file_identity(executable, hash_contents=True)
        version = subprocess.run([binary['path'], '-version'], check=True,
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                 text=True, timeout=10)
        tools['blastn'] = dict(binary=binary, version=version.stdout.strip(),
                              version_stderr=version.stderr.strip())
        dependencies['RNA'] = _module_identity('RNA')
        dependencies['RNA._RNA'] = _module_identity('RNA._RNA')
    return dict(sources=sources, dependencies=dependencies, tools=tools,
                python=dict(version=sys.version, implementation=sys.implementation.name,
                            executable=_file_identity(sys.executable, hash_contents=True)),
                folding_model=dict(parameters='DNA_Mathews2004', dangles=2, noGU=1,
                                   gquad=0, temperature=policy.get('folding_temperature'),
                                   sodium_M=policy.get('sodium_M')))


def _read_json(path):
    path = Path(path)
    if path.is_symlink():
        raise ValueError('Checkpoint metadata must not be a symbolic link: ' + str(path))
    with path.open('rb') as stream:
        raw = stream.read(MAX_METADATA_BYTES + 1)
    if len(raw) > MAX_METADATA_BYTES:
        raise ValueError('Checkpoint metadata exceeds its bounded size limit.')
    value = json.loads(raw)
    if not isinstance(value, dict):
        raise ValueError('Checkpoint metadata is not an object.')
    return value


def _atomic_json(path, value):
    path = Path(path)
    if path.is_symlink():
        raise ValueError('Checkpoint metadata must not be a symbolic link: ' + str(path))
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode='w', encoding='utf-8',
                                         prefix='.' + path.name + '.', suffix='.tmp',
                                         dir=path.parent, delete=False) as stream:
            temporary = Path(stream.name)
            json.dump(value, stream, sort_keys=True, separators=(',', ':'), allow_nan=False)
            stream.write('\n')
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
        descriptor = os.open(str(path.parent), os.O_RDONLY)
        try:
            os.fsync(descriptor)
        finally:
            os.close(descriptor)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


class CheckpointSession:
    """Validate an exact resume identity before the evidence store is opened."""

    def __init__(self, path, config, candidates, databases, policy, engine_identity=None):
        self.path = Path(path)
        if self.path.is_symlink():
            raise ValueError('Checkpoint directory must not be a symbolic link.')
        self.path.mkdir(parents=True, exist_ok=True)
        self.manifest_path = self.path / 'checkpoint_manifest.json'
        self.progress_path = self.path / 'checkpoint_progress.json'
        self.is_resume = self.manifest_path.exists()
        detail = config.get('output_detail', 'standard')
        if detail not in ('standard', 'full_audit'):
            raise ValueError('Unknown checkpoint output detail.')
        self.store_path = self.path / ('specificity_scores.sqlite' if detail=='standard' else 'retained_hits.sqlite')
        self.resumable = detail == 'standard'
        self.resume_reason = ('Saved checkpoints require matching inputs, tools and databases.'
                              if self.resumable else
                              'Full-audit checkpoints cannot be resumed; to preserve every raw audit artifact, start a new design.')
        if self.is_resume and not self.resumable:
            raise ValueError(self.resume_reason)
        self.identity = dict(
            config_sha256=_digest(config), policy=_json_value(policy),
            candidates=[dict(id=candidate['id'], sha256=_digest({key: candidate[key]
                         for key in CANDIDATE_FIELDS if key in candidate})) for candidate in candidates],
            databases={source: database_identity(database) for source, database in sorted(databases.items())},
            engine=_json_value(engine_identity) if engine_identity is not None else scientific_engine_identity(policy),
            output_detail=detail,
        )
        self.identity_digest = _digest(self.identity)
        if self.is_resume:
            saved = _read_json(self.manifest_path)
            if saved.get('schema') != SCHEMA:
                raise ValueError('Checkpoint format is not recognized; start a new design.')
            if saved.get('identity_sha256') != _digest(saved.get('identity')):
                raise ValueError('Checkpoint identity metadata is damaged; saved evidence was not reused.')
            if saved['identity_sha256'] != self.identity_digest:
                changed = [key for key in self.identity if saved['identity'].get(key) != self.identity[key]]
                raise ValueError('Checkpoint identity mismatch (' + ', '.join(changed) +
                                 '); inputs, candidate order, scientific code, tools and databases must match. Start a new design.')
            if self.store_path.is_symlink() or not self.store_path.is_file():
                raise ValueError('Checkpoint score store is missing or unsafe; interrupted evidence cannot be reused.')
            self.created_utc = saved['created_utc']
            self.progress = _read_json(self.progress_path)
            if self.progress.get('schema') != SCHEMA or self.progress.get('identity_sha256') != self.identity_digest:
                raise ValueError('Checkpoint progress does not match its immutable identity.')
            if self.progress.get('status') != 'INCOMPLETE':
                raise ValueError('Checkpoint progress must remain INCOMPLETE until final export verification.')
        else:
            known_stores=(self.path/'specificity_scores.sqlite',self.path/'retained_hits.sqlite')
            if any(item.exists() or item.is_symlink() for item in known_stores) or self.progress_path.exists():
                raise ValueError('Checkpoint evidence exists without its identity; start a new design.')
            self.created_utc = _utc_now()
            _atomic_json(self.manifest_path, dict(schema=SCHEMA, identity=self.identity,
                         identity_sha256=self.identity_digest, created_utc=self.created_utc))
            self.progress = {}
            self.update_progress(completed_batches=0, total_batches=0, folded_candidates=0,
                                 total_candidates=len(candidates), summary='Checkpoint initialized; analysis is incomplete.')

    def update_progress(self, **fields):
        """Publish bounded informational progress after successful DB commits."""
        progress = dict(self.progress, **_json_value(fields))
        progress.update(schema=SCHEMA, status='INCOMPLETE', identity_sha256=self.identity_digest,
                        created_utc=self.created_utc, updated_utc=_utc_now(),
                        resumable=self.resumable, resume_reason=self.resume_reason)
        for name in ('completed_batches', 'total_batches', 'folded_candidates', 'total_candidates'):
            value = progress.get(name, 0)
            if type(value) is not int or value < 0:
                raise ValueError('Checkpoint progress count must be a nonnegative integer: ' + name)
            progress[name] = value
        if len(json.dumps(progress, ensure_ascii=False).encode('utf-8')) > MAX_METADATA_BYTES:
            raise ValueError('Checkpoint progress exceeds its bounded size limit.')
        _atomic_json(self.progress_path, progress)
        self.progress = progress
        return dict(progress)


def read_checkpoint_summary(path):
    """Read only small progress metadata; never open or validate a large store."""
    path = Path(path)
    manifest_path = path / 'checkpoint_manifest.json'
    if not manifest_path.exists():
        return None
    try:
        if path.is_symlink():
            raise ValueError('Checkpoint directory is a symbolic link.')
        manifest = _read_json(manifest_path)
        if manifest.get('schema') != SCHEMA:
            return None
        if manifest.get('identity_sha256') != _digest(manifest.get('identity')):
            raise ValueError('Checkpoint identity metadata is damaged.')
        progress = _read_json(path / 'checkpoint_progress.json')
        if progress.get('schema') != SCHEMA or progress.get('identity_sha256') != manifest.get('identity_sha256'):
            raise ValueError('Checkpoint metadata does not match.')
        if progress.get('status') != 'INCOMPLETE':
            raise ValueError('Checkpoint progress has an invalid completion claim.')
        summary = {key: progress[key] for key in (
            'completed_batches', 'total_batches', 'folded_candidates', 'total_candidates',
            'created_utc', 'updated_utc', 'summary', 'resumable', 'resume_reason') if key in progress}
        for name in ('completed_batches', 'total_batches', 'folded_candidates', 'total_candidates'):
            if type(summary.get(name)) is not int or summary[name] < 0:
                raise ValueError('Invalid checkpoint progress count.')
        summary.update(schema=SCHEMA, status='INCOMPLETE')
        if manifest['identity'].get('output_detail') != 'standard':
            summary.update(resumable=False, resume_reason='Full-audit checkpoints cannot be resumed; to preserve every raw audit artifact, start a new design.')
        store = path / ('specificity_scores.sqlite' if manifest['identity'].get('output_detail')=='standard' else 'retained_hits.sqlite')
        if store.is_symlink() or not store.is_file():
            summary.update(resumable=False, resume_reason='Specificity score store is missing or unsafe.')
        return summary
    except (OSError, ValueError, TypeError, KeyError) as exc:
        return dict(schema=SCHEMA, status='INCOMPLETE', resumable=False,
                    completed_batches=0, total_batches=0, folded_candidates=0,
                    total_candidates=0, summary='Checkpoint metadata cannot be verified.',
                    resume_reason=str(exc))
