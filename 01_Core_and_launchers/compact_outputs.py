"""Small, reproducible run bundles; output choices never change QC decisions.

Standard bundles retain every actual input, selected-probe QC, rejection counts,
and pair warning evidence. Full audit additionally keeps compressed detailed
candidate, alignment, and search artifacts. No large table is assembled in RAM
merely to serialize it.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import shutil
from collections import Counter
from pathlib import Path, PurePosixPath

SCHEMA = 'raserfish.run-manifest.v1'
DETAIL_MODES = ('standard', 'full_audit')
MAX_ASSET_BYTES = 256 * 1024 * 1024


def digest(data):
    return hashlib.sha256(data.encode('utf-8') if isinstance(data, str) else data).hexdigest()


def file_digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def output_detail(config):
    detail = config.get('output_detail', 'standard')
    if detail not in DETAIL_MODES:
        raise ValueError('Output detail must be standard or full_audit.')
    return detail


def _write_json(path, value, compressed=False):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if compressed else open
    with opener(path, 'wt', encoding='utf-8') as stream:
        json.dump(value, stream, separators=(',', ':'), allow_nan=False)


def _write_jsonl(path, records):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, 'wt', encoding='utf-8') as stream:
        for record in records:
            json.dump(record, stream, separators=(',', ':'), allow_nan=False)
            stream.write('\n')


def _write_fasta(path, records):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, 'wt', encoding='ascii', newline='\n') as stream:
        for name, sequence in records:
            stream.write('>' + name + '\n')
            for i in range(0, len(sequence), 80):
                stream.write(sequence[i:i + 80] + '\n')


def write_input_manifest(config, bins, out, application_version=None):
    """Save complete inputs once, and exact designed slices with original bounds.

    Full context is retained when a target slice does not cover the loaded
    reference: exact-site/local-QC results depend on that surrounding context.
    The small manifest never embeds reference sequence or feature annotations.
    """
    out = Path(out)
    detail = output_detail(config)
    assets = {}
    annotations = {}
    context_sequences = []
    extra_sequences = []
    slice_sequences = []
    targets = []
    sources = []
    exact_cover = {}
    lookup = {r['id']: r for r in config['references']}
    for n, region in enumerate(bins):
        reference = lookup[region['reference']]
        sequence = reference['sequence'][region['start'] - 1:region['end']]
        record = 'region_' + str(n + 1)
        slice_sequences.append((record, sequence))
        targets.append(dict(id=region['id'], reference=reference['id'],
                            start=region['start'], end=region['end'],
                            length=len(sequence), sha256=digest(sequence),
                            asset='inputs/target_region.fa.gz', record=record))
        if region['start'] == 1 and region['end'] == len(reference['sequence']):
            exact_cover[reference['id']] = record

    def externalize(value, path, sequence_group=None):
        if isinstance(value, list):
            return [externalize(v, path + '/' + str(i), sequence_group) for i, v in enumerate(value)]
        if not isinstance(value, dict):
            return value
        result = {}
        for key, item in value.items():
            pointer = path + '/' + key
            if key == 'features' and isinstance(item, list):
                annotations[pointer] = item
                result[key] = {'$annotation': pointer}
            elif key == 'sequence' and sequence_group and isinstance(item, str):
                # Top-level target sequence can reuse an exact full-length slice.
                reference_id = value.get('id')
                is_reference = path.startswith('/references/') and path.count('/') == 2
                if is_reference and reference_id in exact_cover:
                    asset, record = 'inputs/target_region.fa.gz', exact_cover[reference_id]
                else:
                    records = context_sequences if is_reference else extra_sequences
                    asset = 'inputs/reference_context.fa.gz' if is_reference else 'inputs/additional_sequences.fa.gz'
                    record = ('context_' if is_reference else 'input_') + str(len(records) + 1)
                    records.append((record, item))
                result[key] = {'$sequence': record, 'asset': asset, 'length': len(item), 'sha256': digest(item)}
            else:
                result[key] = externalize(item, pointer, sequence_group)
        return result

    configuration = {}
    for key, value in config.items():
        configuration[key] = externalize(value, '/' + key, key in ('references', 'background', 'transcripts'))
    for reference in config['references']:
        sources.append(dict(id=reference['id'], description=reference.get('description', ''),
                            length=len(reference['sequence']), sha256=digest(reference['sequence']),
                            mapping=reference.get('mapping'), feature_count=len(reference.get('features', []))))
    for name, records in (('inputs/target_region.fa.gz', slice_sequences),
                          ('inputs/reference_context.fa.gz', context_sequences),
                          ('inputs/additional_sequences.fa.gz', extra_sequences)):
        if records:
            _write_fasta(out / name, records)
            assets[name] = dict(sha256=file_digest(out / name), format='fasta.gz', records=len(records))
    if annotations:
        name = 'inputs/annotations.json.gz'
        _write_json(out / name, annotations, compressed=True)
        assets[name] = dict(sha256=file_digest(out / name), format='json.gz')
    manifest = dict(schema=SCHEMA, application_version=application_version, output_detail=detail, configuration=configuration,
                    input_config_sha256=digest(json.dumps(config, sort_keys=True)),
                    sources=sources, target_regions=targets, assets=assets,
                    context_note='Original local coordinates and complete loaded context are retained for exact-site and local QC reproducibility. External database paths identify dependencies; databases are not bundled.')
    _write_json(out / 'run_config.json', manifest)
    return manifest


def _safe_asset(root, name):
    if not isinstance(name, str) or '\\' in name:
        raise ValueError('Manifest asset path must be a safe relative path.')
    path = PurePosixPath(name)
    if path.is_absolute() or '..' in path.parts or not path.parts or ':' in path.parts[0]:
        raise ValueError('Manifest asset path must be a safe relative path.')
    resolved = (root / Path(*path.parts)).resolve()
    if not resolved.is_relative_to(root):
        raise ValueError('Manifest asset escapes the run directory.')
    return resolved


def _read_limited_gzip(path, remaining_bytes=None):
    limit = MAX_ASSET_BYTES if remaining_bytes is None else min(MAX_ASSET_BYTES, remaining_bytes)
    with gzip.open(path, 'rb') as stream:
        data = stream.read(limit + 1)
    if len(data) > limit:
        raise ValueError('Manifest input asset exceeds the supported decompressed size.')
    return data.decode('utf-8'), len(data)


def load_input_manifest(path, aggregate_decompressed_limit=None):
    """Restore the original configuration only after validating all input hashes.

    ``path`` is the run_config.json file, with its inputs/ directory alongside.
    Plain legacy configurations are accepted unchanged for backward compatibility.
    No path in a manifest can leave that directory, including via symlinks.
    Web uploads can supply an aggregate decompression budget in addition to the
    per-asset limit, preventing many small nested gzip files from expanding huge.
    """
    path = Path(path).resolve()
    if aggregate_decompressed_limit is not None and (isinstance(aggregate_decompressed_limit, bool)
            or not isinstance(aggregate_decompressed_limit, int) or aggregate_decompressed_limit < 0):
        raise ValueError('Aggregate decompressed input limit must be a nonnegative integer.')
    manifest = json.loads(path.read_text())
    if manifest.get('schema') != SCHEMA:
        if 'schema' in manifest:
            raise ValueError('Unsupported run manifest schema.')
        output_detail(manifest)
        return manifest
    root = path.parent
    fasta = {}
    annotations = {}
    expanded_bytes = 0
    for name, metadata in manifest['assets'].items():
        asset = _safe_asset(root, name)
        if file_digest(asset) != metadata['sha256']:
            raise ValueError('Input asset SHA-256 mismatch: ' + name)
        remaining = None if aggregate_decompressed_limit is None else aggregate_decompressed_limit - expanded_bytes
        text, asset_bytes = _read_limited_gzip(asset, remaining)
        expanded_bytes += asset_bytes
        if metadata['format'] == 'fasta.gz':
            records = {}
            current = None
            for line in text.splitlines():
                if line.startswith('>'):
                    current = line[1:]
                    if current in records or not current:
                        raise ValueError('Duplicate or empty manifest FASTA record.')
                    records[current] = []
                elif current is None:
                    raise ValueError('Invalid manifest FASTA asset.')
                else:
                    records[current].append(line)
            fasta[name] = {key: ''.join(lines) for key, lines in records.items()}
        elif metadata['format'] == 'json.gz' and name == 'inputs/annotations.json.gz':
            annotations = json.loads(text)
        else:
            raise ValueError('Unsupported manifest input asset format.')

    def restore(value):
        if isinstance(value, list):
            return [restore(item) for item in value]
        if not isinstance(value, dict):
            return value
        if '$sequence' in value:
            try:
                sequence = fasta[value['asset']][value['$sequence']]
            except KeyError as exc:
                raise ValueError('Missing manifest sequence record.') from exc
            if len(sequence) != value['length'] or digest(sequence) != value['sha256']:
                raise ValueError('Manifest sequence SHA-256 or length mismatch.')
            return sequence
        if '$annotation' in value:
            if value['$annotation'] not in annotations:
                raise ValueError('Missing manifest annotation record.')
            return annotations[value['$annotation']]
        return {key: restore(item) for key, item in value.items()}

    configuration = restore(manifest['configuration'])
    if digest(json.dumps(configuration, sort_keys=True)) != manifest['input_config_sha256']:
        raise ValueError('Restored input configuration SHA-256 mismatch.')
    # Also verify the independently exported target slices against restored context.
    refs = {r['id']: r for r in configuration['references']}
    for region in manifest['target_regions']:
        sequence = refs[region['reference']]['sequence'][region['start'] - 1:region['end']]
        if (sequence != fasta[region['asset']][region['record']] or
                len(sequence) != region['length'] or digest(sequence) != region['sha256']):
            raise ValueError('Target slice differs from restored source context.')
    output_detail(configuration)
    return configuration


_LARGE_DATA_KEYS = frozenset(('sequence', 'features', 'annotated_rna_sequences', 'offtarget_hits', 'qseq', 'sseq'))
DISABLED_CROSS_PRIMARY = 'NOT RUN - USER DISABLED'


def _without_large_data(value):
    # Filter keys before descending: never traverse or deepcopy the excluded HSP
    # collections, which can be much larger than the remaining result.
    if isinstance(value, list):
        return [_without_large_data(item) for item in value]
    if not isinstance(value, dict):
        return value
    return {key: _without_large_data(item) for key, item in value.items()
            if key not in _LARGE_DATA_KEYS}


def pair_flags(pair, policy):
    """Reporting margins only. These flags never participate in acceptance."""
    flags = []
    if pair.get('blocked'):
        flags.append('FAIL_COMPATIBILITY')
    stem = pair.get('max_exact_stem')
    if stem is not None and stem >= policy['maximum_pair_stem'] - 1:
        flags.append('NEAR_EXACT_STEM_LIMIT')
    gain = pair.get('mfe_gain_vs_separate')
    if gain is not None and pair.get('intermolecular_pairs') and gain <= policy['minimum_cofold_gain'] + 1.0:
        flags.append('NEAR_COFOLD_LIMIT')
    if str(pair.get('status', '')).startswith(('ERROR', 'FAILED', 'UNRESOLVED')):
        flags.append('UNRESOLVED_PAIR_CHECK')
    return flags


def _pair_scope(pair, by_id):
    if pair.get('a') == pair.get('b'):
        return 'SELF_DIMER'
    if pair.get('scope') in ('WITHIN_REGION', 'BETWEEN_REGION'):
        return pair['scope']
    a, b = by_id.get(pair.get('a'), {}), by_id.get(pair.get('b'), {})
    if a.get('bin_id') is not None and b.get('bin_id') is not None:
        return 'WITHIN_REGION' if a['bin_id'] == b['bin_id'] else 'BETWEEN_REGION'
    return 'CROSS_PRIMARY'


def _pair_summary(pairs, flagged):
    gains = [p['mfe_gain_vs_separate'] for p in pairs if p.get('mfe_gain_vs_separate') is not None]
    return dict(evaluated=len(pairs), blocked=sum(bool(p.get('blocked')) for p in pairs),
                flagged=len(flagged), max_exact_stem=max((p.get('max_exact_stem', 0) for p in pairs), default=0),
                minimum_mfe_gain_vs_separate=min(gains, default=None))


_CONCISE_OMITTED_KEYS = (_LARGE_DATA_KEYS - {'sequence'}) | frozenset(('base_unpaired', 'structure', 'site_centers'))


def concise_saved_result(result):
    """Slim an already-finalized result without recomputing its QC accounting.

    This also handles existing compact-result files. Candidate totals, rejection
    counts and actual checks are copied unchanged; only detailed display payloads
    and repeated per-batch command receipts are reduced. No raw hit collection is
    traversed. Full audit calls should retain the original detailed result.
    """
    def clean(value):
        if isinstance(value, list):
            return [clean(item) for item in value]
        if not isinstance(value, dict):
            return value
        output = {}
        for key, item in value.items():
            if key in _CONCISE_OMITTED_KEYS:
                continue
            if key == 'intermolecular_pairs' and isinstance(item, (list, tuple)):
                output[key] = len(item)
            elif key == 'blast_runs' and isinstance(item, list):
                output[key] = [summarize_search(run) for run in item]
            else:
                output[key] = clean(item)
        return output

    def summarize_search(run):
        if not isinstance(run, dict):
            return clean(run)
        compact = {key: clean(value) for key, value in run.items() if key != 'batches'}
        if 'search_manifest' in compact:
            compact['search_manifest'] = None
            compact['batch_detail_retention'] = 'NOT RETAINED - CONCISE OUTPUT; select Full audit before a new run for individual batch evidence.'
        batches = run.get('batches')
        if isinstance(batches, list):
            digest = hashlib.sha256()
            encoder = json.JSONEncoder(sort_keys=True, separators=(',', ':'), allow_nan=False)
            statuses = Counter()
            for batch in batches:
                for chunk in encoder.iterencode(batch):
                    digest.update(chunk.encode('utf-8'))
                digest.update(b'\n')
                statuses[str(batch.get('status', 'NOT RECORDED'))] += 1
            compact['batch_receipt_summary'] = dict(count=len(batches), statuses=dict(statuses),
                ordered_receipts_sha256=digest.hexdigest(),
                digest_scope='CANONICAL_JSON_RECEIPTS_NEWLINE_DELIMITED',
                detail='Individual batch receipts omitted from concise result; all alignments were processed before scores were saved.')
            if compact.get('command'):
                compact['command_scope'] = 'FIRST_BATCH_EXAMPLE; common scientific search options retained; individual command receipts omitted from concise result'
        return compact

    compact = clean(result)
    compact['output_detail'] = 'standard'
    compact['retained_detail'] = 'CONCISE_SCORES'
    compact['detail_note'] = 'Scalar QC scores, selected sequences and binding locations retained. Per-base folding values, structures, pair-coordinate lists and repeated batch command receipts require Full audit selected before a run.'
    return compact


def compact_result(result, manifest):
    policy = result['policy']
    selected_ids = set(result['selected_ids'])
    candidates = result['candidates']
    selected = [dict(_without_large_data(c), offtarget_alignment_count=c['offtarget_alignment_count'] if 'offtarget_alignment_count' in c else len(c.get('offtarget_hits', [])))
                for c in candidates if c['id'] in selected_ids]
    statuses = Counter(c.get('status', 'UNKNOWN') for c in candidates)
    rejections = Counter(reason for c in candidates for reason in set(c.get('reject_reasons', [])))
    by_id = {c['id']: c for c in candidates}
    pairs = [dict(pair, scope=_pair_scope(pair, by_id)) for pair in result.get('selected_pair_interactions', [])]
    # Missing setting means a historical result with the original enabled checks.
    # This does not apply the new default retroactively to old saved evidence.
    cross_enabled = policy.get('screen_cross_primary', True)
    self_pairs = [p for p in pairs if p['scope'] == 'SELF_DIMER']
    cross_pairs = [p for p in pairs if p['scope'] != 'SELF_DIMER']
    if not cross_enabled and cross_pairs:
        raise ValueError('Disabled cross-primary screening cannot export evaluated cross-primary pairs.')
    flagged = [dict(pair, report_flags=flags) for pair in pairs if (flags := pair_flags(pair, policy))]
    self_summary = _pair_summary(self_pairs, [p for p in flagged if p['scope'] == 'SELF_DIMER'])
    self_summary.update(result.get('self_dimer_screening', {}))
    self_summary['selected_evaluated'] = len(self_pairs)
    self_summary.setdefault('status', 'RECORDED_IN_LEGACY_RESULT' if self_pairs else 'NOT RECORDED')
    cross_summary = _pair_summary(cross_pairs, [p for p in flagged if p['scope'] != 'SELF_DIMER'])
    cross_summary.update(result.get('cross_primary_screening', {}))
    cross_summary.update(enabled=cross_enabled, selected_evaluated=len(cross_pairs), scope=['within_region', 'between_region'])
    if not cross_enabled:
        cross_summary.update(status=DISABLED_CROSS_PRIMARY, within_region_status=DISABLED_CROSS_PRIMARY,
                             between_region_status=DISABLED_CROSS_PRIMARY, evaluated=0, selected_evaluated=0,
                             max_exact_stem=None, minimum_mfe_gain_vs_separate=None)
    else:
        cross_summary.setdefault('status', 'RECORDED_IN_LEGACY_RESULT')
        cross_summary.setdefault('within_region_status', cross_summary['status'])
        cross_summary.setdefault('between_region_status', cross_summary['status'])
    summary = _pair_summary(pairs, flagged)
    summary.update(status='SELF-DIMER CHECKS ONLY' if not cross_enabled else 'SELF-DIMER AND CROSS-PRIMARY CHECKS',
                   self_dimer=self_summary, cross_primary=cross_summary,
                   cross_primary_status=cross_summary['status'],
                   within_region_status=cross_summary['within_region_status'],
                   between_region_status=cross_summary['between_region_status'],
                   reporting_margins=dict(exact_stem_bases=1, cofold_kcal_per_mol=1.0),
                   note='Counts cover only analyses actually run. Cross-primary checks are separate from self-dimer checks. Flag margins are display-only.')
    compact = {key: _without_large_data(value) for key, value in result.items()
               if key not in ('candidates', 'selected_pair_interactions', 'references')}
    compact['references'] = [dict(source, sequence_asset='run_config.json') for source in manifest['sources']]
    compact['imagers'] = result.get('imagers', [])
    compact.update(candidates=selected, selected_pair_interactions=flagged,
                   selected_pair_summary=summary, cross_primary_screening=cross_summary,
                   self_dimer_screening=self_summary, order_status='DRAFT',
                   candidate_summary=dict(total=len(candidates), statuses=dict(statuses), rejection_reasons=dict(rejections)),
                   output_detail=manifest['output_detail'], run_manifest='run_config.json',
                   result_schema='raserfish.compact-result.v1')
    for candidate in selected:
        flags = []
        if candidate.get('reject_reasons'):
            flags.append('FAIL_INDIVIDUAL_QC')
        if candidate.get('specificity', {}).get('moderate_tm_sites', 0):
            flags.append('MODERATE_OFFTARGETS_PRESENT')
        mean = candidate.get('fold', {}).get('mean_arm_unpaired')
        if mean is not None and mean <= policy.get('minimum_arm_unpaired', 0.4) + 0.05:
            flags.append('NEAR_UNPAIRED_LIMIT')
        flags.extend(pair_flags(candidate.get('self_pair', {}), policy))
        candidate['report_flags'] = flags
    return concise_saved_result(compact) if manifest['output_detail'] == 'standard' else compact


def _csv(path, rows, fields, compressed=False):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if compressed else open
    with opener(path, 'wt', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction='ignore')
        writer.writeheader()
        for row in rows:
            cleaned = {}
            for key in fields:
                value = row.get(key, '')
                if isinstance(value, (dict, list, tuple)):
                    value = json.dumps(value, separators=(',', ':'))
                if isinstance(value, str) and value.startswith(('=', '+', '@')):
                    value = "'" + value
                cleaned[key] = value
            writer.writerow(cleaned)


def archive_transient(path, out, detail):
    """Compress or discard a finalized search artifact, using bounded memory."""
    path, out = Path(path), Path(out)
    if not path.exists():
        return None
    if detail == 'full_audit':
        target = out / 'debug' / (path.name + '.gz')
        target.parent.mkdir(parents=True, exist_ok=True)
        with path.open('rb') as source, gzip.open(target, 'wb') as dest:
            shutil.copyfileobj(source, dest, 1024 * 1024)
        path.unlink()
        return str(target.relative_to(out))
    path.unlink()
    return None


def finalize_outputs(result, manifest, out):
    out = Path(out)
    detail = manifest['output_detail']
    compact = compact_result(result, manifest)
    fields = ['scope', 'a', 'b', 'max_exact_stem', 'mfe_gain_vs_separate', 'structure',
              'intermolecular_pairs', 'blocked', 'status', 'report_flags']
    _csv(out / 'flagged_pair_interactions.csv', compact['selected_pair_interactions'], fields)
    _write_json(out / 'pair_interaction_summary.json', compact['selected_pair_summary'])
    _csv(out / 'pair_interaction_summary.csv', [compact['selected_pair_summary']],
         ['status', 'evaluated', 'blocked', 'flagged', 'max_exact_stem', 'minimum_mfe_gain_vs_separate', 'cross_primary_status', 'within_region_status', 'between_region_status', 'self_dimer', 'cross_primary', 'reporting_margins', 'note'])
    _write_json(out / 'cross_primary_screening.json', compact['cross_primary_screening'])
    _write_json(out / 'self_dimer_screening.json', compact['self_dimer_screening'])
    _write_json(out / 'candidate_summary.json', compact['candidate_summary'])
    _csv(out / 'selected_qc_flags.csv', compact['candidates'], ['id', 'name', 'bin_id', 'report_flags', 'specificity', 'fold', 'self_pair', 'reject_reasons'])
    if detail == 'full_audit':
        _write_jsonl(out / 'debug' / 'all_candidates.jsonl.gz',
                     ({key: value for key, value in candidate.items() if key != 'offtarget_hits'} for candidate in result['candidates']))
        _write_jsonl(out / 'debug' / 'candidate_review_alignments.jsonl.gz',
                     (dict(candidate_id=c['id'], bin_id=c['bin_id'], hit=h) for c in result['candidates'] for h in c.get('offtarget_hits', [])))
        by_id = {c['id']: c for c in result['candidates']}
        _csv(out / 'debug' / 'all_selected_pair_interactions.csv.gz', (dict(p, scope=_pair_scope(p, by_id)) for p in result.get('selected_pair_interactions', [])), fields, compressed=True)
    for run in compact.get('blast_runs', []):
        original = run.get('raw_file')
        if original:
            # A search implementation may already archive immediately after parsing.
            if original.startswith('debug/') and (out / original).is_file():
                archived = original
            else:
                archived = archive_transient(out / original, out, detail)
            run['raw_file'] = archived
        run['raw_retention'] = 'COMPRESSED_FULL_AUDIT' if detail == 'full_audit' else 'DISCARDED_STANDARD'
        run['note'] = 'Search counts, command and SHA-256 retained. Raw HSPs are compressed only in full audit mode; heuristic searches do not prove exhaustive absence of binding.'
    for path in list(out.glob('*_candidates_blast.tsv')) + list(out.glob('*_blast.log')) + list(out.glob('candidate_queries.fa')):
        archive_transient(path, out, detail)
    _write_json(out / 'design_result.json', compact)
    return compact
