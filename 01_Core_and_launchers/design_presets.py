"""Starting values with field-level provenance, not a claimed full paper pipeline."""
from __future__ import annotations

LEGACY_PRESETS = {
    'adaptive_45_55': {
        'label': 'Adaptive high-yield / 45-55 C (user settings)',
        'description': 'Prefer short arms in 36-42 nt, but evaluate longer alternatives up to the adjustable 90-nt ceiling. '
                       'Every retained alternative must meet the SAME RNA, specificity, structure and spacing criteria. '
                       'Select more compatible binding sites; the count target is a benchmark, not a stopping cap. '
                       'These are user-requested departures from the published 42-50 C / 36-42 nt window.',
        'values': dict(arm_min=36, preferred_arm_max=42, arm_max=90, length_mode='adaptive',
                       extension_shortlist_per_bin=80, selection_goal='maximize_sites',
                       gc_min=0.0, gc_max=76.0, tm_min=45.0, tm_max=55.0, tm_strict=True,
                       profile='oligominer1_2nM', homopolymer_max=5, minimum_gap=5, footprint='extended',
                       folding_temperature=42.0, sodium_M=0.390),
    },
    'ellenberg_informed': {
        'label': 'Ellenberg / LoopTrace-informed starting parameters',
        'description': 'Published primary-arm length and Tm window. Unreported details use clearly labeled project settings. '
                       'Your five-base footprint spacing and separately chosen RNA policy are retained; this is not an exact reproduction of the paper.',
        'values': dict(arm_min=36, arm_max=42, length_mode='fixed', preferred_arm_max=42, extension_shortlist_per_bin=80, selection_goal='target_count', gc_min=0.0, gc_max=76.0,
                       tm_min=42.0, tm_max=50.0, tm_strict=True, profile='oligominer1_2nM',
                       homopolymer_max=5, minimum_gap=5, footprint='extended',
                       folding_temperature=42.0, sodium_M=0.390),
    },
    'expanded_g5': {
        'label': 'Expanded G5-length alternative (36-90 nt; not the initial pilot)',
        'description': 'Optional broader-length project alternative. Uses the current nearest-neighbor model during selection, '
                       'not the original approximate GC/length Tm formula. No silent extension of the paper-length preset.',
        'values': dict(arm_min=36, arm_max=90, length_mode='fixed', preferred_arm_max=42, extension_shortlist_per_bin=80, selection_goal='target_count', gc_min=0.0, gc_max=76.0,
                       tm_min=42.0, tm_max=50.0, tm_strict=True, profile='oligominer1_2nM',
                       homopolymer_max=5, minimum_gap=5, footprint='extended',
                       folding_temperature=42.0, sodium_M=0.390),
    },
}

# All active controls share the verified ERRFI1 run; old IDs remain importable.
import json
from pathlib import Path
ERRFI1_DEFAULTS = json.loads((Path(__file__).resolve().parent/'inputs/errfi1_defaults.json').read_text())
PRESETS = {'errfi1': {
    'label': 'ERRFI1 / DP1 — verified run defaults',
    'description': 'Saved ERRFI1 scientific thresholds unchanged. Top (+), transcript screening ON, whole-gene exclusion OFF; preferred 36–42 nt, adaptive through 90 nt; strict 45–55 C. Search budget increased to 400 preferred + 200 extended candidates in 1 kb sectors.',
    'values': ERRFI1_DEFAULTS['policy'],
}}

SOURCES = {
    'errfi1_run': 'Verified saved ERRFI1 / DP1 run 20261015_100112_62f902; see verification/errfi1_baseline.json. v1.6 changes only search budgets/methods and output defaults.',
    'user_update': 'User request, 13-14 October 2026: fixed trace-wide strand, transcript-only RNA policy available/default in 1.5; permit longer arms when useful to meet remaining criteria and increase compatible binding-site count; change selection Tm window to 45-55 C.',
    'paper': 'Beckwith et al., Nature Communications 16:6673 (2025), Primary FISH probe design and Non-denaturing FISH; DOI 10.1038/s41467-025-61689-y.',
    'initial_record': 'RASERFISH_primary_probe_design_parameters_and_logic.pdf, sections 3-6: initial 68-oligo pilot used mostly 36-42 nt candidates; four long R2 additions were exceptions.',
    'project_audit': 'REVIEW_README.txt, G5 arm-selection audit: 36-90 nt, GC <=76%, no explicit lower GC cutoff, no homopolymer >=6; inherited expanded-project settings, not claimed original paper settings.',
}

PARAMETER_NOTES = [
    dict(setting='Preferred DNA-binding-arm length', value='36-42 nt, with adaptive longer alternatives', origin='PAPER BASE RANGE / USER ADAPTATION',
         explanation='The preferred short range is retained. New candidates may be extended and rechecked to recover useful sites; an existing ordered sequence is never silently lengthened.'),
    dict(setting='Maximum adaptive arm length', value='90 nt, adjustable', origin='PROJECT SEARCH CEILING, NOT PAPER LIMIT',
         explanation='This is the maximum genomic arm, not the full primary. The program checks lengths one base at a time and favors short feasible alternatives; it does not turn every probe into a 90-mer.'),
    dict(setting='Current on-target Tm selection window', value='45-55 C; strict endpoints', origin='NEW USER REQUIREMENT',
         explanation='This replaces 42-50 C for the active default. Salt, concentration, formamide correction, aqueous folding temperature and off-target thresholds are unchanged.'),
    dict(setting='Count-selection goal', value='More compatible exact binding sites', origin='NEW USER REQUIREMENT',
         explanation='Compare a short-arm baseline, longer infill and two joint greedy searches, followed by bounded spatial repair. Keep the highest-site-count feasible result. Count targets are benchmarks unless capped-count mode is selected; no global-optimum or brighter-signal guarantee.'),
    dict(setting='Published comparison window (not the current user default)', value='42-50 C; strict software comparisons', origin='PAPER WINDOW / SOFTWARE BOUNDS',
         explanation='The legacy Ellenberg-informed comparison reports 42-50 C accounting for 2xSSC and 50% formamide. Strict endpoints are inherited OligoMiner-style implementation details.'),
    dict(setting='Tm model and concentration', value='DNA_NN3; 390 mM sodium; 2 nM probe / trace target; 50% FA x 0.65 C/%', origin='EXPLICIT IMPLEMENTATION / PROJECT ASSUMPTIONS',
         explanation='The paper does not give its complete thermodynamic command line. 2 nM per primary is a reported experimental concentration; using it with trace target is our modeling choice, not proof of the paper\'s calculation. Default 25/25 nM software comparison remains selectable.'),
    dict(setting='Arm GC limits', value='No independent minimum (0%); maximum 76%', origin='INHERITED PROJECT FALLBACK',
         explanation='The paper does not report primary-arm GC cutoffs. The 76% cap is documented in our expanded G5 audit; a precise cap was not established for the earliest pilot. The paper\'s 40-50% GC rule applies to short imagers, NOT primary binding arms. Tm still constrains low-GC candidates.'),
    dict(setting='Primary footprint spacing', value='At least 5 bp; extended footprint', origin='USER REQUIREMENT, STRICTER THAN PAPER',
         explanation='The paper used at least 2 nt spacing. We retain your requested 5 bp and include contiguous spacer/handle matching. Tracing regions themselves remain adjacent.'),
    dict(setting='Identical-base / simple-repeat filters', value='At most 5 consecutive identical bases; exclude 8 dinucleotide units', origin='INHERITED PROJECT HEURISTICS',
         explanation='These are project sequence-quality filters. The paper\'s separate short-imager filters are not substituted for primary-arm rules.'),
    dict(setting='Region width and target oligo count', value='User-chosen separately in steps 2 and 4', origin='USER LAYOUT / PAPER EXAMPLE AVAILABLE',
         explanation='The paper reports up to 100 primaries in 10 kb per position. A separate button applies that planning example. It is not forced on short construct inputs or on a user-chosen tracing resolution.'),
    dict(setting='Trace-wide strand and RNA policy', value='One fixed DNA strand; transcript-sequence screening ON; whole-gene exclusion OFF by default', origin='USER REQUIREMENT / REVISED PROJECT POLICY',
         explanation='All tracing probes keep the same chosen DNA-bound strand. Gene spans do not automatically exclude candidates in the default mode. Long/high-identity complementary RNA matches still affect selection; isolated 18-nt hits are review evidence. Strict span exclusion remains optional. No RNA-free or allele-identity guarantee.'),
    dict(setting='Off-target / folding / pair thresholds', value='Displayed project QC heuristics retained', origin='PROJECT HEURISTICS, NOT PAPER REPRODUCTION',
         explanation='The paper reports a 42 C LDA model at 0.8 and 5 or 10 off-target 18-mers. This program does not implement that complete LDA/k-mer pipeline. Its BLAST, off-target Tm, aqueous folding and pair gates are labeled separately.'),
]


def preset_provenance(policy: dict) -> dict:
    ident=policy.get('preset_id','custom')
    preset=PRESETS.get(ident,LEGACY_PRESETS.get(ident,{}))
    changes={k:{'preset':v,'actual':policy.get(k)} for k,v in preset.get('values',{}).items() if policy.get(k)!=v}
    return dict(preset_id=ident,label=preset.get('label','Custom settings'),modified_fields=changes,
                notes=PARAMETER_NOTES,sources=SOURCES,
                explicit_nonimplementation=['Complete Ellenberg command line not recovered','OligoMiner LDA not implemented','Whole-genome k-mer pipeline not implemented'],
                statement='Reported paper values are distinguished from project/user choices. Selecting a preset is not experimental validation.')
