"""Small, human-readable result scores; never recompute or relax scientific QC.

The caller verifies that the run finished. Saved values are copied into explicit
columns, not combined into an invented probability or quality score. Large
alignment, per-base, and pair-coordinate collections are never traversed here.
"""
from __future__ import annotations

import csv
import json
import math
import statistics
from pathlib import Path

SCHEMA = 'raserfish.concise-scores.v1'
DISABLED = 'NOT RUN - USER DISABLED'
PRIMARY_SCORE_FIELDS = (
    'primary_id', 'name', 'region', 'imager', 'decision', 'reject_reasons',
    'tm_C', 'gc_percent', 'arm_length_nt', 'bound_strand', 'exact_sites',
    'high_tm_offtarget_sites', 'moderate_tm_offtarget_sites',
    'unresolved_long_sites', 'strong_rna_risk_sites', 'rna_risk_sites',
    'individual_folding_status', 'arm_accessibility', 'fold_mfe_kcal_mol',
    'self_dimer_status', 'self_dimer_blocked', 'self_dimer_max_exact_stem_bp',
    'self_dimer_gain_kcal_mol', 'self_dimer_intermolecular_pair_count',
    'cross_primary_status', 'order_status',
)


def _clean_name(value):
    return ' '.join(str(value or '').split())[:180]


def _regions(value):
    """Use actual bins; never include stale assignments outside designed bounds."""
    if value.get('bins'):
        return value['bins']
    if value.get('configuration'):
        value = value['configuration']
    try:
        from sequence_tools import make_bins
        refs = [dict(r, length=r.get('length') or len(r.get('sequence', '')))
                for r in value.get('references', [])]
        return make_bins(refs, int(value.get('bin_width', 20000)),
                         value.get('tail_policy', 'strict'), value.get('assignments'),
                         value.get('bounds'))[0]
    except (TypeError, ValueError, KeyError):
        return []


def human_readable_run_name(result_or_config, run_id=None):
    """Name from designed readouts/coordinates, ignoring old preset titles."""
    value = result_or_config.get('configuration', result_or_config)
    custom = _clean_name(value.get('run_name'))
    if custom:
        return custom
    bins = _regions(result_or_config)
    imagers = list(dict.fromkeys(_clean_name(b.get('imager')) for b in bins if b.get('imager')))
    locations = []
    for region in bins:
        if region.get('chromosome') and region.get('genomic_start') is not None and region.get('genomic_end') is not None:
            location = f"{region['chromosome']}:{region['genomic_start']}-{region['genomic_end']}"
        elif region.get('reference') and region.get('start') is not None and region.get('end') is not None:
            location = f"{region['reference']}:{region['start']}-{region['end']} (local)"
        else:
            continue
        if location not in locations:
            locations.append(location)
    if not locations:
        return 'RASERFISH design' + (' · ' + str(run_id) if run_id else '')
    readout = ', '.join(imagers[:4]) + (f' +{len(imagers)-4}' if len(imagers) > 4 else '')
    where = '; '.join(locations[:2]) + (f'; +{len(locations)-2} regions' if len(locations) > 2 else '')
    return _clean_name((readout + ' · ' if readout else '') + where)


def _scalar_fields(value, fields):
    if not isinstance(value, dict):
        return {}
    return {key: value[key] for key in fields if key in value and
            isinstance(value[key], (str, int, float, bool, type(None)))}


def _range(values):
    numbers = [x for x in values if isinstance(x, (float, int)) and not isinstance(x, bool) and math.isfinite(x)]
    return dict(minimum=min(numbers), median=statistics.median(numbers), maximum=max(numbers)) if numbers else None


def summarize_completed_result(result, run_id, config=None):
    """Extract concise scores from a verified completed result (including legacy).

    Known incomplete results are refused. Legacy completed results can lack the
    modern analysis_complete field; completion must be established by the caller
    from the run's durable result/export evidence, not inferred from this summary.
    """
    status = str(result.get('status', ''))
    if result.get('analysis_complete') is False or any(word in status.upper() for word in ('INCOMPLETE', 'STOPPED', 'FAILED', 'CANCELLED')):
        raise ValueError('An incomplete analysis cannot be summarized as a completed run.')
    selected_ids = set(result.get('selected_ids', []))
    candidates = [c for c in result.get('candidates', []) if c.get('id') in selected_ids or (not selected_ids and c.get('status') == 'selected')]
    policy = result.get('policy', {})
    cross_recorded = result.get('cross_primary_screening') or result.get('selected_pair_summary', {}).get('cross_primary', {})
    disabled = policy.get('screen_cross_primary') is False or cross_recorded.get('enabled') is False
    cross = dict(enabled=False, status=DISABLED, within_region_status=DISABLED,
                 between_region_status=DISABLED, evaluated=0) if disabled else _scalar_fields(cross_recorded, (
                     'enabled', 'status', 'within_region_status', 'between_region_status', 'evaluated', 'selected_evaluated', 'blocked'))
    cross.setdefault('status', 'NOT RECORDED')
    rows = []
    for candidate in candidates:
        specificity = candidate.get('specificity') or {}
        fold = candidate.get('fold') or {}
        pair = candidate.get('self_pair') or {}
        inter = pair.get('intermolecular_pairs')
        count = len(inter) if isinstance(inter, (list, tuple)) else inter if isinstance(inter, int) else None
        rows.append(dict(zip(PRIMARY_SCORE_FIELDS, (
            candidate.get('id'), candidate.get('name') or candidate.get('id'), candidate.get('bin_id'),
            candidate.get('imager'), candidate.get('status', 'NOT RECORDED'), candidate.get('reject_reasons', []),
            candidate.get('tm'), candidate.get('gc'), candidate.get('arm_length'), candidate.get('bound_strand'),
            len(candidate.get('sites', [])), specificity.get('high_tm_sites'), specificity.get('moderate_tm_sites'),
            specificity.get('unknown_long_sites'), specificity.get('strong_rna_sites'), specificity.get('rna_risk_sites'),
            fold.get('status', 'NOT RECORDED'), fold.get('mean_arm_unpaired'), fold.get('mfe_kcal_mol'),
            pair.get('status', 'NOT RECORDED'), pair.get('blocked'), pair.get('max_exact_stem'),
            pair.get('mfe_gain_vs_separate'), count, cross['status'], 'DRAFT',
        ))))
    regions = []
    for region in result.get('bins', []):
        item = _scalar_fields(region, ('id', 'reference', 'imager', 'start', 'end', 'length', 'chromosome', 'genomic_start', 'genomic_end'))
        item['metrics'] = _scalar_fields(region.get('metrics', {}), (
            'probes', 'exact_sites', 'coverage_percent', 'largest_gap_bp', 'median_spacing_bp', 'median_tm',
            'shortfall_to_target', 'below_minimum', 'extended_primaries', 'min_arm_length', 'max_arm_length',
            'top_strand_sites', 'bottom_strand_sites'))
        regions.append(item)
    counts = _scalar_fields(result.get('counts', {}), ('candidates', 'selected', 'regions', 'exact_sites',
            'extended_primaries', 'QC_rejected', 'monomer_evaluations', 'self_dimer_evaluations', 'cross_primary_evaluations'))
    summary = dict(schema=SCHEMA, run_id=run_id,
        name=human_readable_run_name({**result, 'run_name': (config or {}).get('run_name') or result.get('run_name')}, run_id),
        generated_utc=result.get('generated_utc'), status='complete', analysis_status=status,
        order_status='DRAFT', counts=counts, regions=regions,
        selected_score_ranges={key: _range(row[key] for row in rows) for key in ('tm_C', 'gc_percent', 'arm_accessibility', 'self_dimer_gain_kcal_mol')},
        cross_primary_screening=cross,
        self_dimer_screening=_scalar_fields(result.get('self_dimer_screening') or result.get('selected_pair_summary', {}).get('self_dimer', {}),
            ('status', 'evaluated', 'selected_evaluated', 'blocked', 'max_exact_stem', 'minimum_mfe_gain_vs_separate')),
        checks=[_scalar_fields(check, ('name', 'status')) for check in result.get('checks', [])],
        selected_primary_scores=rows,
        note='Existing QC values are reported separately, without an invented composite score. Omitted analysis is NOT RUN. Order remains DRAFT.')
    return summary


def write_run_summary(result, run_dir, run_id=None, config=None):
    """Write score JSON/CSV and a readable label beside an existing result."""
    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    summary = summarize_completed_result(result, run_id or run_dir.parent.name, config)
    (run_dir / 'run_summary.json').write_text(json.dumps(summary, separators=(',', ':'), allow_nan=False) + '\n')
    with (run_dir / 'primary_scores.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=PRIMARY_SCORE_FIELDS)
        writer.writeheader()
        for row in summary['selected_primary_scores']:
            flat = {key: json.dumps(value, separators=(',', ':')) if isinstance(value, (list, dict)) else value for key, value in row.items()}
            for key, value in flat.items():
                if isinstance(value, str) and value.startswith(('=', '+', '@', '\t', '\r')):
                    flat[key] = "'" + value
            writer.writerow(flat)
    (run_dir / 'RUN_NAME.txt').write_text(summary['name'] + '\n')
    return summary
