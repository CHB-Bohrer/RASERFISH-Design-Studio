"""Read-only coverage summaries from committed results and bounded saved annotation.

No sequence assets, alignments, BLAST databases or checkpoint stores are opened.
Coordinates and transcript strands in this API are reference-local, 1-based inclusive.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path

SCHEMA = 'raserfish.coverage-report.v1'
DISABLED = 'NOT RUN - USER DISABLED'
MAX_SOURCES = 200
MAX_DIRECTORY_ENTRIES = 2000
MAX_REGIONS = 64
MAX_FEATURES = 50000
MAX_CANDIDATES = 50000
MAX_SITES = 100000
MAX_JSON_BYTES = 4 * 1024 * 1024
MAX_ANNOTATION_BYTES = 8 * 1024 * 1024
MAX_REPORT_BYTES = 4 * 1024 * 1024
MAX_RESULT_BYTES = 16 * 1024 * 1024
NOTE = ('Saved annotation describes reference transcript models, not expression or RNA binding. '
        'Introns are inferred within each isoform; collapsed gene tracks are not single transcripts. '
        'RNA with the same reference-local strand as the bound DNA is complementary to the arm. '
        'No annotated overlap does not establish absence of transcription.')


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _integer(value, label):
    _require(isinstance(value, int) and not isinstance(value, bool), 'Invalid ' + label)
    return value


def _identifier(value):
    _require(isinstance(value, str) and re.fullmatch(r'[A-Za-z0-9_-]{1,100}', value), 'Invalid coverage source ID')
    return value


def _root(parent, ident):
    parent = Path(parent)
    _identifier(ident)
    _require(not parent.is_symlink(), 'Coverage sources must not use symlinks')
    root = parent / ident
    _require(not root.is_symlink() and root.is_dir() and root.resolve().parent == parent.resolve(),
             'Unknown coverage source')
    return root


def _asset(root, relative):
    root = Path(root)
    _require(not root.is_symlink(), 'Saved coverage assets must not use symlinks')
    path = root
    for part in Path(relative).parts:
        _require(part not in ('..', '.', '') and not Path(part).is_absolute(), 'Unsafe coverage asset')
        path /= part
        _require(not path.is_symlink(), 'Saved coverage assets must not use symlinks')
    _require(path.resolve().is_relative_to(root.resolve()), 'Unsafe coverage asset')
    return path


def _bytes(path, limit=MAX_JSON_BYTES):
    _require(not path.is_symlink() and path.is_file(), 'Missing or unsafe saved coverage asset')
    _require(path.stat().st_size <= limit, 'Saved coverage asset exceeds the read limit')
    with path.open('rb') as stream:
        data = stream.read(limit + 1)
    _require(len(data) <= limit, 'Saved coverage asset exceeds the read limit')
    return data


def _json(path, limit=MAX_JSON_BYTES):
    value = json.loads(_bytes(path, limit))
    _require(isinstance(value, dict), 'Invalid saved coverage JSON')
    return value


def _merge(intervals):
    merged = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return merged


def _clip(start, end, lo, hi):
    start, end = max(start, lo), min(end, hi)
    return [start, end] if start <= end else None


def _annotation(result, result_dir):
    """Resolve only features, never externalized sequence or alignment assets."""
    refs = result.get('references', [])
    by_id = {r['id']: r for r in refs}
    config_refs = {}
    annotations = {}
    if result_dir is not None:
        config_path = _asset(result_dir, 'run_config.json')
        if config_path.is_file():
            saved = _json(config_path)
            configuration = saved.get('configuration', saved)
            _require(isinstance(configuration, dict), 'Invalid saved annotation configuration')
            config_refs = {r['id']: r for r in configuration.get('references', [])}
            annotation_path = _asset(result_dir, 'inputs/annotations.json.gz')
            if annotation_path.is_file():
                compressed = _bytes(annotation_path)
                expected = saved.get('assets', {}).get('inputs/annotations.json.gz', {}).get('sha256')
                if expected:
                    _require(hashlib.sha256(compressed).hexdigest() == expected, 'Saved annotation checksum mismatch')
                with gzip.open(annotation_path, 'rb') as stream:
                    expanded = stream.read(MAX_ANNOTATION_BYTES + 1)
                _require(len(expanded) <= MAX_ANNOTATION_BYTES, 'Saved annotation exceeds the expanded read limit')
                annotations = json.loads(expanded)
                _require(isinstance(annotations, dict), 'Invalid saved annotation')
    output = {}
    feature_count = 0
    for rid, ref in by_id.items():
        configured = config_refs.get(rid, {})
        raw = configured.get('features', ref.get('features'))
        if isinstance(raw, dict) and set(raw) == {'$annotation'}:
            raw = annotations.get(raw['$annotation'])
        metadata = configured.get('gene_annotation', ref.get('gene_annotation', {}))
        metadata = metadata if isinstance(metadata, dict) else {}
        available = isinstance(raw, list) and (bool(raw) or metadata.get('status') == 'FETCHED')
        features = raw if isinstance(raw, list) else []
        feature_count += len(features)
        _require(feature_count <= MAX_FEATURES, 'Saved annotation contains too many features')
        # Peak/signal features alone do not establish transcript annotation availability.
        if features and not any(f.get('type') in ('gene', 'transcript', 'exon', 'mRNA') for f in features):
            available = metadata.get('status') == 'FETCHED'
        source = str(metadata.get('provider') or next((f.get('source') for f in features if f.get('source')), 'Not recorded'))
        output[rid] = dict(features=features if available else [], status='AVAILABLE' if available else 'UNKNOWN',
                          source=source, note=NOTE if available else 'Saved transcript annotation is unavailable; overlap counts are unknown. ' + NOTE)
    return output


def _feature_interval(feature, reference):
    # Stored local coordinates/strands are already reversed by annotation_sources
    # for minus-oriented references. Genomic coordinates are only a fallback.
    start, end = feature.get('start'), feature.get('end')
    strand = feature.get('strand', '?')
    if start is None or end is None:
        mapping = reference.get('mapping') or {}
        gs, ge = feature.get('genomic_start'), feature.get('genomic_end')
        _require(mapping.get('verified') is True and gs is not None and ge is not None,
                 'Annotation feature has no local coordinates')
        if mapping.get('orientation', '+') == '-':
            start, end = mapping['end'] - ge + 1, mapping['end'] - gs + 1
            strand = {'+': '-', '-': '+'}.get(strand, '?')
        else:
            start, end = gs - mapping['start'] + 1, ge - mapping['start'] + 1
    start, end = _integer(start, 'feature start'), _integer(end, 'feature end')
    _require(start <= end, 'Invalid annotation interval')
    return start, end, strand if strand in ('+', '-') else '?'


def _tracks(features, reference, start, end):
    genes, txs, exons = {}, {}, defaultdict(list)
    for feature in features:
        kind = feature.get('type')
        if kind not in ('gene', 'transcript', 'mRNA', 'exon'):
            continue
        a, b, strand = _feature_interval(feature, reference)
        fid = str(feature.get('feature_id') or feature.get('id') or '')
        item = dict(id=fid, label=str(feature.get('label') or fid), strand=strand,
                    start=a, end=b, parent=str(feature.get('parent') or ''))
        if kind == 'gene' and fid:
            genes[fid] = item
        elif kind in ('transcript', 'mRNA') and fid:
            txs[fid] = item
        elif kind == 'exon' and item['parent']:
            exons[item['parent']].append((a, b, strand))
    grouped = defaultdict(list)
    full_introns = []
    uncertain_intervals = []
    for tid, tx in txs.items():
        span = _clip(tx['start'], tx['end'], start, end)
        if span is None:
            continue
        blocks = []
        for a, b, strand in exons[tid]:
            _require(tx['start'] <= a <= b <= tx['end'] and strand == tx['strand'],
                     'Exon disagrees with its transcript model')
            blocks.append([a, b])
        merged = _merge(blocks)
        if not merged:
            uncertain_intervals.append((tx['start'], tx['end'], tx['strand']))
        else:
            # Overlap/region snapshots can omit exons outside the reference.
            # A transcript's unbounded terminal stretches therefore cannot be
            # called exon-free introns, even when some internal exons were saved.
            if tx['start'] < merged[0][0]:
                uncertain_intervals.append((tx['start'], merged[0][0] - 1, tx['strand']))
            if merged[-1][1] < tx['end']:
                uncertain_intervals.append((merged[-1][1] + 1, tx['end'], tx['strand']))
        # Gaps between saved exons of the same isoform remain positive evidence,
        # including when its distant terminal exons were outside the snapshot.
        for left, right in zip(merged, merged[1:]):
            if left[1] + 1 <= right[0] - 1:
                full_introns.append((left[1] + 1, right[0] - 1, tx['strand']))
                if tx['strand'] not in ('+', '-'):
                    uncertain_intervals.append((left[1] + 1, right[0] - 1, '?'))
        visible = [p for a, b in merged if (p := _clip(a, b, start, end)) is not None]
        gid = tx['parent'] or 'transcript:' + tid
        grouped[gid].append(dict(id=tid, label=tx['label'], strand=tx['strand'], start=span[0], end=span[1], exons=visible))
    out = []
    for gid, transcripts in grouped.items():
        gene = genes.get(gid, {})
        strands = {t['strand'] for t in transcripts}
        strand = gene.get('strand', next(iter(strands)) if len(strands) == 1 else '?')
        _require(all(t['strand'] == strand for t in transcripts), 'Gene and transcript strands disagree')
        out.append(dict(id=gid, label=gene.get('label') or gid.removeprefix('transcript:'), strand=strand,
                        spans=_merge([[t['start'], t['end']] for t in transcripts]),
                        exons=_merge([e for t in transcripts for e in t['exons']]),
                        transcripts=sorted(transcripts, key=lambda t: (t['start'], t['end'], t['id']))))
    for gid, gene in genes.items():
        span = _clip(gene['start'], gene['end'], start, end)
        if span and gid not in grouped:
            uncertain_intervals.append((gene['start'], gene['end'], gene['strand']))
            out.append(dict(id=gid, label=gene['label'], strand=gene['strand'],
                            spans=[span], exons=[], transcripts=[]))
    return sorted(out, key=lambda g: (g['spans'][0][0], g['label'], g['id'])), full_introns, uncertain_intervals


def _mapping(reference, start, end):
    mapping = reference.get('mapping') or {}
    if mapping.get('verified') is not True or not mapping.get('chromosome'):
        return None, None, None, 'Reference-local coordinates (1-based inclusive)'
    _require(mapping.get('orientation', '+') in ('+', '-'), 'Invalid reference mapping orientation')
    if mapping.get('orientation', '+') == '-':
        gs, ge = mapping['end'] - end + 1, mapping['end'] - start + 1
    else:
        gs, ge = mapping['start'] + start - 1, mapping['start'] + end - 1
    chromosome = str(mapping['chromosome'])
    label = f"{mapping.get('assembly', 'Genome')} {chromosome}:{gs:,}-{ge:,} · reference {mapping.get('orientation', '+')}"
    return chromosome, gs, ge, label


def build_report(result, result_dir=None, title=None):
    """Derive presentation data without changing a result or recalculating QC."""
    _require(isinstance(result, dict), 'Invalid coverage result')
    _require(result.get('analysis_complete') is not False, 'Coverage requires completed results')
    _require(not any(s in str(result.get('status', '')).upper() for s in ('INCOMPLETE', 'FAILED', 'CANCELLED', 'STOPPED')),
             'Coverage requires completed results')
    refs = result.get('references', [])
    bins = result.get('bins', [])
    candidates = result.get('candidates', [])
    _require(isinstance(refs, list) and isinstance(bins, list) and isinstance(candidates, list), 'Invalid saved result collections')
    _require(0 < len(bins) <= MAX_REGIONS and len(refs) <= MAX_REGIONS and len(candidates) <= MAX_CANDIDATES,
             'Coverage result exceeds the collection limit')
    by_ref = {r['id']: r for r in refs}
    _require(len(by_ref) == len(refs), 'Duplicate reference IDs')
    _require(len({r['id'] for r in bins}) == len(bins), 'Duplicate region IDs')
    selected_ids = result.get('selected_ids', [])
    _require(isinstance(selected_ids, list) and len(selected_ids) <= MAX_CANDIDATES and len(set(selected_ids)) == len(selected_ids),
             'Invalid selected primary IDs')
    selected_set = set(selected_ids)
    selected = [c for c in candidates if c.get('id') in selected_set]
    _require(len(selected) == len(selected_ids) and len({c['id'] for c in selected}) == len(selected), 'Selected primary records are missing or duplicated')
    if 'selected' in result.get('counts', {}):
        _require(result['counts']['selected'] == len(selected), 'Saved selected count differs from primary records')
    annotations = _annotation(result, result_dir)
    by_bin = defaultdict(list)
    for candidate in selected:
        by_bin[candidate.get('bin_id')].append(candidate)
    regions = []
    site_count = 0
    warnings = ['Coverage is a read-only view of saved selections; no design or QC was rerun.',
                'Reference-local mean offsets describe sequence positions, not measured three-dimensional centroid errors.', NOTE]
    for region in bins:
        rid, refid = region['id'], region['reference']
        _require(refid in by_ref, 'Region reference is missing')
        ref = by_ref[refid]
        start, end = _integer(region['start'], 'region start'), _integer(region['end'], 'region end')
        _require(1 <= start <= end <= ref['length'], 'Invalid region bounds')
        sites = []
        bound = result.get('policy', {}).get('bound_strand')
        for c in by_bin[rid]:
            for index, site in enumerate(c.get('sites', [])):
                a, b = _integer(site['start'], 'site start'), _integer(site['end'], 'site end')
                strand = site.get('strand', c.get('bound_strand'))
                _require(site.get('reference', c.get('reference')) == refid and start <= a <= b <= end,
                         'Selected site is outside its assigned region')
                _require(strand in ('+', '-'), 'Selected site strand is unknown')
                if bound not in ('+', '-'):
                    bound = strand
                _require(strand == bound, 'Selected sites do not share the bound DNA strand')
                sites.append(dict(id=c['id'] + ':' + str(index + 1), name=str(c.get('name') or c['id']),
                                  start=a, end=b, strand=strand))
        site_count += len(sites)
        _require(site_count <= MAX_SITES, 'Coverage contains too many sites')
        sites.sort(key=lambda s: (s['start'], s['end'], s['id']))
        centers = [(s['start'] + s['end']) / 2 for s in sites]
        midpoint = (start + end) / 2
        _require((end - start) // 1000 + 1 <= 10000, 'Region exceeds the coverage sector limit')
        sectors = [dict(start=a, end=min(end, a + 999), count=0) for a in range(start, end + 1, 1000)]
        for center in centers:
            sectors[int((center - start) // 1000)]['count'] += 1
        covered = _merge([[s['start'], s['end']] for s in sites])
        ann = annotations[refid]
        genes, introns, uncertain_intervals = _tracks(ann['features'], ref, start, end) if ann['status'] == 'AVAILABLE' else ([], [], [])
        known_introns = 0
        uncertain_sites = 0
        for site in sites:
            known = any(a <= site['start'] <= site['end'] <= b and strand == site['strand']
                        for a, b, strand in introns)
            uncertain = (ann['status'] != 'AVAILABLE' or
                         any(a <= site['start'] <= site['end'] <= b and strand in (site['strand'], '?')
                             for a, b, strand in uncertain_intervals))
            # Any proven complementary isoform match settles this site's count;
            # another truncated isoform cannot erase that positive evidence.
            site['intron_status'] = 'COMPLETE_INTRON_MATCH' if known else 'UNKNOWN' if uncertain else 'NO_COMPLETE_INTRON_MATCH'
            known_introns += int(known)
            uncertain_sites += int(uncertain and not known)
        complete_introns = known_introns if ann['status'] == 'AVAILABLE' and not uncertain_sites else None
        region_annotation = {k: ann[k] for k in ('status', 'source', 'note')}
        if uncertain_sites and ann['status'] == 'AVAILABLE':
            region_annotation.update(status='PARTIAL', note='Some selected sites have uncertain complementary-strand intron classification because transcript, exon or strand data are missing, or the saved model ends before the transcript boundary. Confirmed interior intron matches are retained separately. ' + NOTE)
        chromosome, gs, ge, coordinate_label = _mapping(ref, start, end)
        regions.append(dict(id=str(rid), label=str(region.get('label') or rid), imager=region.get('imager'),
            reference=refid, start=start, end=end, length=end - start + 1,
            chromosome=chromosome, genomic_start=gs, genomic_end=ge, coordinate_label=coordinate_label,
            mapping_orientation=(ref.get('mapping') or {}).get('orientation', '+') if chromosome else None,
            bound_strand=bound if bound in ('+', '-') else None, sites=sites, sectors=sectors,
            metrics=dict(primaries=len(by_bin[rid]), sites=len(sites),
                left_half=sum(c < midpoint for c in centers), right_half=sum(c >= midpoint for c in centers),
                mean_offset_bp=sum(c - midpoint for c in centers) / len(centers) if centers else None,
                largest_gap_bp=max((b[0] - a[1] - 1 for a, b in zip(covered, covered[1:])), default=0) if covered else None,
                left_edge_gap_bp=covered[0][0] - start if covered else end - start + 1,
                right_edge_gap_bp=end - covered[-1][1] if covered else end - start + 1,
                complete_intron_sites=complete_introns,
                known_complete_intron_sites=known_introns if ann['status'] == 'AVAILABLE' else None,
                uncertain_intron_sites=uncertain_sites), genes=genes,
            annotation=region_annotation))
    _require(set(by_bin) <= {b['id'] for b in bins}, 'Selected primary has no assigned region')
    if 'exact_sites' in result.get('counts', {}):
        _require(result['counts']['exact_sites'] == site_count, 'Saved site count differs from binding-site records')
    statuses = {r['annotation']['status'] for r in regions}
    status = next(iter(statuses)) if len(statuses) == 1 else 'PARTIAL'
    if status != 'AVAILABLE':
        warnings.append('Some saved transcript annotation is unavailable; unknown overlap counts are shown explicitly.')
    cross = result.get('cross_primary_screening', {})
    cross_status = DISABLED if result.get('policy', {}).get('screen_cross_primary') is False else cross.get('status', 'NOT RECORDED')
    for region in regions:
        region['cross_primary_status'] = cross_status
    if cross_status == DISABLED:
        warnings.append('Cross-primary interactions within and between regions: ' + DISABLED + '.')
    report = dict(schema=SCHEMA, title=str(title or result.get('run_name') or result.get('title') or 'Saved probe coverage'),
        order_status='DRAFT', cross_primary_status=cross_status, regions=regions,
        annotation=dict(status=status, source='; '.join(sorted({r['annotation']['source'] for r in regions})), note=NOTE),
        warnings=warnings)
    _require(len(json.dumps(report, separators=(',', ':')).encode()) <= MAX_REPORT_BYTES, 'Coverage report exceeds the response limit')
    return report


def _completed(jobs, ident, completed_run):
    root = _root(jobs, ident)
    result_path = _asset(root, 'result/design_result.json')
    _require(result_path.is_file() and result_path.stat().st_size <= MAX_RESULT_BYTES,
             'Saved result exceeds the coverage read limit')
    saved = completed_run(ident, include_result=True)
    _require(saved.get('status') == 'complete' and isinstance(saved.get('result'), dict), 'Coverage source is incomplete')
    return root, saved


def _experiment(base, jobs, ident, completed_run):
    root = _root(Path(base) / 'experiments', ident)
    manifest = _json(_asset(root, 'experiment_manifest.json'))
    _require(manifest.get('analysis_complete') is True and manifest.get('queue_status') == 'COMPLETE',
             'Coverage experiment is incomplete')
    targets = manifest.get('targets', [])
    _require(isinstance(targets, list) and 0 < len(targets) <= MAX_REGIONS, 'Invalid experiment target list')
    completed = []
    seen = set()
    for target in targets:
        _require(target.get('status') == 'COMPLETE' and target.get('analysis_complete') is True, 'Experiment target is incomplete')
        run_id = _identifier(target.get('studio_run_id'))
        _require(run_id not in seen, 'Experiment repeats a saved run')
        seen.add(run_id)
        expected = Path(jobs) / run_id
        _require(target.get('source_run_dir') == str(expected), 'Experiment saved-run provenance does not match')
        run_root, saved = _completed(jobs, run_id, completed_run)
        result = saved['result']
        if target.get('summary', {}).get('result_sha256'):
            actual = hashlib.sha256(_bytes(_asset(run_root, 'result/design_result.json'), 32 * 1024 * 1024)).hexdigest()
            _require(actual == target['summary']['result_sha256'], 'Experiment result checksum mismatch')
        _require(len(result.get('bins', [])) == 1, 'Experiment target must identify one designed region')
        region = result['bins'][0]
        reference = next((r for r in result.get('references', []) if r.get('id') == region.get('reference')), {})
        chromosome, gs, ge, _ = _mapping(reference, region['start'], region['end'])
        _require((region.get('imager'), gs, ge) == (target.get('imager'), target.get('start'), target.get('end'))
                 and chromosome == target.get('chromosome', chromosome), 'Experiment target differs from its saved run')
        # Do not retain full scientific results across runs. Listing an experiment
        # and aggregating its compact plots must remain bounded to one run at a time.
        completed.append((target, run_root, dict(id=saved['id'], name=saved.get('name'),
            sites=result.get('counts', {}).get('exact_sites'))))
        del saved, result, reference, region
    return manifest, completed


def _experiment_label(manifest, completed, ident):
    if manifest.get('title'):
        return str(manifest['title'])
    if len(completed) == 6 and all(t.get('end', 0) - t.get('start', 0) + 1 == 10000 for t, _, _ in completed):
        counts = [meta['sites'] for _, _, meta in completed]
        if all(isinstance(n, int) and n >= 0 for n in counts):
            return 'Six 10 kb targets · ' + str(sum(counts)) + ' sites'
    date = ident[:8]
    date = date[:4] + '-' + date[4:6] + '-' + date[6:] if date.isdigit() and len(date) == 8 else ident
    return str(len(completed)) + ' target experiment · ' + date


def load_source(base, jobs, completed_run, kind, ident):
    _identifier(ident)
    if kind == 'run':
        root, saved = _completed(jobs, ident, completed_run)
        return build_report(saved['result'], root / 'result', saved.get('name'))
    _require(kind == 'experiment', 'Unknown coverage source kind')
    manifest, completed = _experiment(base, jobs, ident, completed_run)
    reports = []
    for target, root, metadata in completed:
        _, saved = _completed(jobs, metadata['id'], completed_run)
        reports.append(build_report(saved['result'], root / 'result', target.get('label')))
        del saved
    regions = []
    for (target, _, _), report in zip(completed, reports):
        for region in report['regions']:
            region['label'] = str(target.get('label') or region['label'])
            regions.append(region)
    _require(len({r['id'] for r in regions}) == len(regions), 'Experiment region IDs collide')
    statuses = {r['annotation']['status'] for r in reports}
    cross = {r['cross_primary_status'] for r in reports}
    warnings = list(dict.fromkeys(w for r in reports for w in r['warnings']))
    aggregate_cross = DISABLED if cross == {DISABLED} else 'BETWEEN SEPARATE RUNS NOT ASSESSED'
    if aggregate_cross != DISABLED:
        warnings.append('Recorded interaction checks apply to each saved run separately. Interactions between separate runs were not assessed; individual region statuses are shown without a combined PASS claim.')
    report = dict(schema=SCHEMA, title=_experiment_label(manifest, completed, ident),
        order_status='DRAFT', cross_primary_status=aggregate_cross,
        regions=regions, annotation=dict(status=next(iter(statuses)) if len(statuses) == 1 else 'PARTIAL',
            source='; '.join(sorted({r['annotation']['source'] for r in reports})), note=NOTE),
        warnings=warnings)
    _require(len(json.dumps(report, separators=(',', ':')).encode()) <= MAX_REPORT_BYTES, 'Coverage report exceeds the response limit')
    return report


def _children(parent):
    parent = Path(parent)
    if not parent.is_dir() or parent.is_symlink():
        return []
    entries = []
    # Stop before directory enumeration or reading result bodies becomes unbounded.
    for index, child in enumerate(parent.iterdir()):
        _require(index < MAX_DIRECTORY_ENTRIES, 'Too many saved coverage sources to list safely')
        if child.is_dir() and not child.is_symlink() and re.fullmatch(r'[A-Za-z0-9_-]{1,100}', child.name):
            entries.append(child)
    return sorted(entries, key=lambda p: p.name, reverse=True)[:MAX_SOURCES]


def list_sources(base, jobs, completed_run):
    sources = []
    for root in _children(Path(base) / 'experiments'):
        try:
            manifest, completed = _experiment(base, jobs, root.name, completed_run)
            sources.append(dict(id=root.name, label=_experiment_label(manifest, completed, root.name), kind='experiment'))
        except (OSError, ValueError, TypeError, KeyError):
            continue
    for root in _children(jobs):
        try:
            _, saved = _completed(jobs, root.name, completed_run)
            sources.append(dict(id=root.name, label=str(saved.get('name') or root.name), kind='run'))
            del saved
        except (OSError, ValueError, TypeError, KeyError):
            continue
    return dict(sources=sources[:MAX_SOURCES])
