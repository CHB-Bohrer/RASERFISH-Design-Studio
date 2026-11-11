#!/usr/bin/env python3
"""Same-folder setup helper for the consolidated Design Studio 1.6.

Only software checks and PUBLIC reference preparation. No design is run, and
no existing probe, target, scientific policy, or historical order is edited.
Called by the bundled install.sh; no old folder or patch is required.
"""
from __future__ import annotations
import argparse
import datetime as dt
import hashlib
import json
import math
import os
from pathlib import Path
import shlex
import shutil
import stat
import subprocess
import sys
import tempfile
import uuid
import zipfile

ACCESSION = 'GCF_000001405.40'
LABELS = ('human_GRCh38', 'rna_human_GRCh38')


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as fh:
        for block in iter(lambda: fh.read(4 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def save_json(path: Path, obj: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.tmp')
    tmp.write_text(json.dumps(obj, indent=2) + '\n', encoding='utf-8')
    tmp.replace(path)


def within(base: Path, relative: str) -> Path:
    """Never follow a reference-manifest entry outside its reference directory."""
    rel = Path(relative)
    if rel.is_absolute() or '..' in rel.parts:
        raise ValueError('Expected a safe relative reference path: ' + relative)
    path = (base / rel).resolve()
    if path != base.resolve() and base.resolve() not in path.parents:
        raise ValueError('Reference path escapes this folder: ' + relative)
    return path


def command(cmd: list[str], root: Path, log: Path) -> None:
    """Stream progress without a shell; preserve command and failures in the log."""
    log.parent.mkdir(parents=True, exist_ok=True)
    line = '$ ' + shlex.join([str(x) for x in cmd])
    print('\n' + line, flush=True)
    with log.open('a', encoding='utf-8') as fh:
        fh.write('\n' + line + '\n')
        fh.flush()
        proc = subprocess.Popen(cmd, cwd=root, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, text=True,
                                encoding='utf-8', errors='replace', bufsize=1)
        try:
            assert proc.stdout is not None
            for text in proc.stdout:
                print(text, end='', flush=True)
                fh.write(text)
                fh.flush()
            code = proc.wait()
        except BaseException:
            proc.terminate()
            try:
                proc.wait(timeout=10)
            except subprocess.TimeoutExpired:
                proc.kill()
                proc.wait()
            raise
        finally:
            if proc.stdout is not None:
                proc.stdout.close()
    if code:
        raise RuntimeError(f'Command failed (exit {code}). See {log}')


def extract_zip_safely(path: Path, destination: Path) -> None:
    """Validate member destinations and CRC while extracting NCBI data."""
    destination.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(path) as z:
        for item in z.infolist():
            within(destination, item.filename)
            if stat.S_ISLNK(item.external_attr >> 16):
                raise ValueError('Unexpected symbolic link in reference ZIP: ' + item.filename)
        z.extractall(destination)


def db_info(prefix: Path) -> str:
    return subprocess.check_output(['blastdbcmd', '-db', str(prefix), '-info'],
                                   text=True, stderr=subprocess.STDOUT)


def protect_inputs(root: Path) -> dict:
    """Hash the application and scientific inputs, ignoring runtime/cache folders."""
    paths = list(root.glob('*.py')) + list(root.glob('*.json'))
    paths += [root / 'start_studio.sh']
    for name in ('inputs', 'existing_design', 'current_design', 'engine', 'web'):
        folder = root / name
        if folder.is_dir():
            paths.extend(p for p in folder.rglob('*') if p.is_file()
                         and '__pycache__' not in p.parts and p.suffix != '.pyc')
    return {str(p.relative_to(root)): sha256(p) for p in sorted(set(paths)) if p.is_file()}


def verify_inputs(root: Path, before: dict) -> None:
    changes = [name for name, value in before.items()
               if not (root / name).is_file() or sha256(root / name) != value]
    if changes:
        raise RuntimeError('Application/input preservation check failed: ' + ', '.join(changes))


def check_engines(root: Path) -> dict:
    """Exercise the existing designer's actual folding and BLAST adapters."""
    sys.path.insert(0, str(root))
    import Bio
    import RNA
    import numpy
    from design_engine import Policy, Folding, blast_candidates
    from sequence_tools import rc
    env = root / '.qc_env'
    if Path(sys.prefix).resolve() != env.resolve():
        raise RuntimeError('Not using the Python environment in THIS Studio folder.')
    for module in (Bio, RNA, numpy):
        location = Path(module.__file__).resolve()
        if env.resolve() not in location.parents:
            raise RuntimeError('A dependency loaded from outside the local .qc_env: ' + str(location))
    for tool in ('blastn', 'makeblastdb', 'blastdbcmd', 'datasets'):
        found = shutil.which(tool)
        if not found or env.resolve() not in Path(found).resolve().parents:
            raise RuntimeError(tool + ' must be installed inside this folder\'s .qc_env.')
    if Bio.__version__ != '1.86' or RNA.__version__ != '2.7.2':
        raise RuntimeError('Expected Biopython 1.86 and ViennaRNA 2.7.2. No scientific setting was relaxed.')
    # Synthetic positive controls, never an analysis of the user's probes.
    arm = 'AGCTTGACGTCAGCTAGGCTCATGACCTAGTCGGATCGTA'
    docking = 'TATTCGTATCGG'
    pre, suf = docking + 'TT', 'TT' + docking
    candidate = dict(id='setup_positive', arm=arm, full=pre + arm + suf,
                     prefix=pre, suffix=suf, bin_id='T_B001', sites=[dict(strand='-')],
                     tm=48.0, offtarget_hits=[], reject_reasons=[])
    other = dict(candidate, id='setup_partner', arm=rc(arm), full=pre + rc(arm) + suf)
    policy = Policy(mode='full', threads=1, bound_strand='-')
    folding = Folding(policy)
    monomer = folding.fold(candidate)
    if monomer['status'] != 'RUN_AQUEOUS_MODEL' or not math.isfinite(monomer['mfe_kcal_mol']):
        raise RuntimeError('Real DNA folding did not complete.')
    if len(monomer['base_unpaired']) != len(candidate['full']):
        raise RuntimeError('Real full-primary folding returned an incomplete sequence.')
    for pair in (folding.pair(candidate, candidate), folding.pair(candidate, other)):
        if pair['status'] != 'RUN_AQUEOUS_MODEL' or not math.isfinite(pair['mfe_gain_vs_separate']):
            raise RuntimeError('Real full-primary cofolding did not complete.')
    ref = dict(id='T', sequence='A' * 99 + arm + 'C' * 100, length=199 + len(arm),
               mapping=dict(verified=True, start=1, end=199 + len(arm), chromosome='chr7',
                            orientation='+', assembly='GRCh38'))
    bins = [dict(id='T_B001', reference='T', start=1, end=ref['length'])]
    parent = root / '.tools' / 'studio_setup' / 'temporary_tests'
    parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='engines_', dir=parent) as tmp:
        folder = Path(tmp)
        fa = folder / 'tiny.fa'
        fa.write_text('>NC_000007.14\n' + ref['sequence'] + '\n')
        subprocess.run(['makeblastdb', '-in', str(fa), '-dbtype', 'nucl', '-parse_seqids',
                        '-blastdb_version', '4', '-out', str(folder / 'db')],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        runs = blast_candidates([candidate], [ref], bins, [], {'human': str(folder / 'db')},
                                policy, folder, lambda *a, **k: None)
        hits = candidate['offtarget_hits']
        if not any(h['allowed'] and not h['full_primary'] and h['query_coverage'] == 100
                   and h['identity'] == 100 for h in hits):
            raise RuntimeError('Real BLAST/parser missed its known exact intended control.')
        if not any(h['full_primary'] for h in hits) or runs[0]['status'] != 'RUN_HEURISTIC_SEARCH':
            raise RuntimeError('Real complete-primary BLAST adapter test failed.')
    return {'status': 'PASS_WITH_REAL_ENGINES', 'Biopython': Bio.__version__,
            'ViennaRNA': RNA.__version__, 'NumPy': numpy.__version__, 'Python': sys.executable,
            'BLAST': subprocess.check_output(['blastn', '-version'], text=True).strip(),
            'datasets': subprocess.check_output(['datasets', '--version'], text=True).strip(),
            'scope': 'Synthetic software tests only; not a probe-QC run.'}


def verify_references(refs: Path) -> dict:
    manifest = json.loads((refs / 'reference_manifest.json').read_text())
    if manifest.get('assembly_accession') != ACCESSION:
        raise RuntimeError('Reference manifest is not for the expected assembly ' + ACCESSION)
    if set(manifest.get('databases', {})) != set(LABELS):
        raise RuntimeError('Both human DNA and annotated RNA databases are required.')
    for label in LABELS:
        database = manifest['databases'][label]
        fasta = within(refs, database['fasta_relative'])
        if not fasta.is_file() or sha256(fasta) != database['fasta_sha256']:
            raise RuntimeError('Missing or changed reference FASTA: ' + str(fasta))
        prefix = within(refs, database['database_prefix_relative'])
        print(db_info(prefix), flush=True)
        if not database.get('index_files'):
            raise RuntimeError('Missing reference index integrity records for ' + label)
        for recorded in database['index_files']:
            index = within(refs, recorded['path_relative'])
            if not index.is_file() or sha256(index) != recorded['sha256']:
                raise RuntimeError('Missing or changed BLAST index: ' + str(index))
    return manifest


def prepare_references(root: Path) -> dict:
    refs = root / 'references' / ACCESSION
    if refs.exists() and root not in refs.resolve().parents:
        raise RuntimeError('References are linked outside this Studio folder; expected real local files.')
    refs.mkdir(parents=True, exist_ok=True)
    ready = refs / 'reference_manifest.json'
    if ready.is_file():
        print('\nChecking the completed local references. No new download is needed.', flush=True)
        return verify_references(refs)
    free = shutil.disk_usage(root).free / 1024 ** 3
    if free < 15:
        raise RuntimeError(f'Only {free:.1f} GiB free. Free space before preparing references (15 GiB minimum check; plan about 30 GB).')
    archive = refs / 'GRCh38_genome_and_RNA.zip'
    complete = refs / 'download_complete.json'
    log = refs / 'PREPARE_LOG.txt'
    cmd = ['datasets', 'download', 'genome', 'accession', ACCESSION,
           '--include', 'genome,rna,seq-report', '--filename', str(archive)]
    if not complete.exists():
        if archive.exists():
            previous = archive.with_name('previous_incomplete_' + uuid.uuid4().hex[:8] + '.zip')
            archive.replace(previous)  # Preserve interrupted downloads; never delete user data.
            print('Retained an incomplete download as ' + str(previous), flush=True)
        command(cmd, root, log)
        save_json(complete, {'assembly': ACCESSION, 'command': cmd, 'archive_sha256': sha256(archive)})
    else:
        info = json.loads(complete.read_text())
        if info.get('assembly') != ACCESSION or not archive.is_file() or sha256(archive) != info['archive_sha256']:
            raise RuntimeError('Completed download record or archive changed. Do not reuse it silently.')
    unpacked = refs / 'downloaded'
    print('\nUnpacking the public genome and annotated RNA files.', flush=True)
    extract_zip_safely(archive, unpacked)
    assembly = unpacked / 'ncbi_dataset' / 'data' / ACCESSION
    genomes, rnas = list(assembly.glob('*_genomic.fna')), list(assembly.glob('rna.fna'))
    if len(genomes) != 1 or len(rnas) != 1:
        raise RuntimeError('Expected one genomic.fna and one rna.fna in the NCBI assembly package.')
    database_root = refs / 'blast'
    database_root.mkdir(exist_ok=True)
    databases = {}
    for label, fasta in zip(LABELS, (genomes[0], rnas[0])):
        prefix = database_root / label
        command(['makeblastdb', '-in', str(fasta), '-dbtype', 'nucl', '-parse_seqids',
                 '-blastdb_version', '4', '-out', str(prefix), '-title', ACCESSION + ' ' + label], root, log)
        info = db_info(prefix)
        print(info, flush=True)
        (refs / (label + '_database_info.txt')).write_text(info)
        indices = sorted(prefix.parent.glob(prefix.name + '.*'))
        if not indices:
            raise RuntimeError('No local index files were produced for ' + label)
        databases[label] = {'fasta_relative': str(fasta.relative_to(refs)), 'fasta_sha256': sha256(fasta),
                            'fasta_bytes': fasta.stat().st_size, 'database_prefix_relative': str(prefix.relative_to(refs)),
                            'database_info': info,
                            'index_files': [{'path_relative': str(p.relative_to(refs)), 'sha256': sha256(p),
                                             'bytes': p.stat().st_size} for p in indices if p.is_file()]}
    result = {'assembly_accession': ACCESSION, 'prepared_utc': dt.datetime.now(dt.timezone.utc).isoformat(),
              'datasets_version': subprocess.check_output(['datasets', '--version'], text=True).strip(),
              'makeblastdb_version': subprocess.check_output(['makeblastdb', '-version'], text=True).strip(),
              'download': json.loads(complete.read_text()), 'databases': databases,
              'scope': 'NCBI GRCh38.p14 genome and available annotated RNA for this download; records metadata and hashes. '
                       'Not a full custom construct or every possible nascent/novel transcript. '
                       'A fresh annotated RNA download can differ from an older download of the same assembly.'}
    save_json(ready, result)
    return result


def verify_studio_discovery(root: Path) -> dict:
    sys.path.insert(0, str(root))
    import studio
    status = studio.environment(root)
    if not status.get('vienna') or not status.get('blastn') or not status.get('biopython'):
        raise RuntimeError('Studio could not discover its locally installed engines.')
    if set(status.get('databases', {})) != {'human', 'rna_human'}:
        raise RuntimeError('Studio did not discover both local reference databases.')
    for prefix in status['databases'].values():
        if (root / 'references').resolve() not in Path(prefix).resolve().parents:
            raise RuntimeError('Studio discovered a reference outside this folder: ' + prefix)
    return status


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=('snapshot', 'engines', 'references', 'finish'))
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--run', type=Path, required=True)
    args = parser.parse_args()
    root, log = args.root.resolve(), args.run.resolve()
    log.mkdir(parents=True, exist_ok=True)
    try:
        if args.stage == 'snapshot':
            save_json(log / 'INPUTS_BEFORE.json', protect_inputs(root))
        elif args.stage == 'engines':
            save_json(log / 'ENGINE_CHECK.json', check_engines(root))
            print('\nReal DNA folding and BLAST adapter tests passed.', flush=True)
        elif args.stage == 'references':
            manifest = prepare_references(root)
            save_json(log / 'REFERENCE_MANIFEST.json', manifest)
        elif args.stage == 'finish':
            before = json.loads((log / 'INPUTS_BEFORE.json').read_text())
            verify_inputs(root, before)
            status = verify_studio_discovery(root)
            save_json(log / 'INSTALLATION_READY.json', {
                'status': 'LOCAL_INSTALLATION_READY_NOT_A_PROBE_QC_RESULT', 'studio_folder': str(root),
                'application_and_scientific_inputs_unchanged': True, 'checked_files': len(before),
                'studio_environment': status, 'engine_check': json.loads((log / 'ENGINE_CHECK.json').read_text()),
                'reference_manifest': str(root / 'references' / ACCESSION / 'reference_manifest.json')})
            print('\nStudio found BOTH local databases. Application and input hashes are unchanged.', flush=True)
        return 0
    except (Exception, KeyboardInterrupt) as exc:
        print('\nSETUP STOPPED: ' + str(exc), file=sys.stderr, flush=True)
        print('Keep the setup_logs directory. No missing analysis was labeled as passed.', file=sys.stderr)
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
