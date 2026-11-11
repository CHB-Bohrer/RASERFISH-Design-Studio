"""Sequence-orientation safeguard for one RASER-FISH tracing set.

A single explicit physical DNA-bound strand is required across the complete run.
RNA annotations reject candidates; they never change this strand. This is a
sequence-design invariant, not certification of chromatid or homolog identity.
"""
from collections import Counter, defaultdict
from sequence_tools import canonical_ref, genomic_interval

ERROR = ("RASER-FISH tracing requires ONE fixed DNA-bound strand for the entire "
         "tracing set. In Step 4 choose Top (+) or Bottom (-) for ALL regions. "
         "Automatic per-gene strand switching is not allowed. RNA conflicts "
         "exclude candidates; they do not switch the target strand.")


def require_fixed_strand(policy):
    strand = policy.get('bound_strand') if isinstance(policy, dict) else policy.bound_strand
    if strand not in ('+', '-'):
        raise ValueError(ERROR)
    return strand


def check_reference_polarity(refs, policy):
    """Reject differently oriented fragments of one mapped chromosome.

    Local references use their supplied orientation. Without a verified mapping
    the program cannot establish how separate construct fragments are connected.
    """
    strand = require_fixed_strand(policy)
    orientations = defaultdict(set)
    for ref in refs:
        mapping = ref.get('mapping')
        if not mapping:
            continue
        orientation = mapping.get('orientation', '+')
        if orientation not in ('+', '-'):
            raise ValueError('Unknown genomic reference orientation for ' + ref['id'])
        if not mapping.get('verified'):
            raise ValueError('Unverified genomic reference orientation for ' + ref['id'])
        chromosome = canonical_ref(mapping['chromosome'])
        orientations[chromosome].add(strand if orientation == '+' else ('-' if strand == '+' else '+'))
    for chromosome, actual in orientations.items():
        if len(actual) > 1:
            raise ValueError('References on ' + chromosome + ' have inconsistent sequence orientations. '
                             'Normalize the fragments to one genomic orientation before tracing; '
                             'the same local Top/Bottom label would otherwise bind opposite genomic strands.')
    return {chromosome: next(iter(actual)) for chromosome, actual in orientations.items()}


def audit_selected_strands(selected, refs, bins, policy):
    """Independently check every selected primary and EVERY exact binding copy."""
    locked = require_fixed_strand(policy)
    genomic = check_reference_polarity(refs, policy)
    ref_lookup = {ref['id']: ref for ref in refs}
    region_lookup = {region['id']: region for region in bins}
    by_region = {region['id']: Counter() for region in bins}
    site_counts = Counter()
    for candidate in selected:
        if candidate.get('bound_strand') != locked:
            raise AssertionError('Single-strand safeguard: wrong/missing primary strand for ' + candidate['id'])
        region_id = candidate['bin_id']
        if region_id not in region_lookup:
            raise AssertionError('Single-strand safeguard: candidate names an unknown region')
        if not candidate.get('sites'):
            raise AssertionError('Single-strand safeguard: selected primary has no mapped sites')
        for site in candidate['sites']:
            if site.get('strand') != locked:
                raise AssertionError('Single-strand safeguard: wrong/missing strand at an exact repeat copy of ' + candidate['id'])
            if site['reference'] not in ref_lookup:
                raise AssertionError('Single-strand safeguard: unmapped reference')
            ref = ref_lookup[site['reference']]
            mapping = genomic_interval(ref, site['start'], site['end'], site['strand'])
            if mapping and mapping['strand'] != genomic[canonical_ref(mapping['reference'])]:
                raise AssertionError('Single-strand safeguard: inconsistent genomic target strand')
            by_region[region_id][locked] += 1
            site_counts[locked] += 1
    rows = []
    for region in bins:
        counts = by_region[region['id']]
        rows.append(dict(Region=region['id'], Reference=region['reference'],
                         Required_bound_strand=locked, Top_sites=counts['+'], Bottom_sites=counts['-'],
                         Status='PASS_SEQUENCE_ORIENTATION' if sum(counts.values()) else 'NO_SELECTED_SITES'))
    return dict(version='1.0', status='PASS_SEQUENCE_ORIENTATION' if selected else 'NO_SELECTED_PROBES',
                required_bound_strand=locked, selected_primaries=len(selected),
                top_sites=site_counts['+'], bottom_sites=site_counts['-'],
                opposite_strand_sites=site_counts['-' if locked == '+' else '+'],
                genomic_strand_by_chromosome=genomic, per_region=rows,
                scope='One fixed reference-relative target strand across this run. '
                      'Mapped fragments on the same chromosome are checked for a consistent genomic strand. '
                      'Unmapped construct-fragment connectivity is not inferred. '
                      'Other independent runs must use the same global strand and be checked together. '
                      'This does not identify maternal/paternal homologs or guarantee one chromatid experimentally.')
