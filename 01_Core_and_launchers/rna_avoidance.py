"""Gene-direction-aware DNA probe selection.

Coordinates and strands here are relative to the supplied forward reference.
A plus-sense transcript requires a minus-bound DNA primary (sense arm).
A minus-sense transcript requires a plus-bound DNA primary (reverse-complement arm).
This is a conservative sequence-design guard, never a guarantee of no RNA binding.
"""
from __future__ import annotations
from collections import defaultdict
from sequence_tools import rc, sha, normalize

TRANSCRIBED_TYPES = {
    'gene', 'transcript', 'mrna', 'ncrna', 'lncrna', 'rrna', 'trna', 'misc_rna',
    'pre_rna', 'primary_transcript', 'precursor_rna', 'cds', 'exon', 'utr',
    '5\'utr', '3\'utr', '5utr', '3utr', 'five_prime_utr', 'three_prime_utr',
}
MODES = ('unknown', 'annotations', 'forward', 'reverse', 'both', 'manual', 'none')


def feature_intervals(ref: dict) -> list[dict]:
    """Gene/transcript spans include introns; histone peaks are never transcripts.

    Joined feature parts with a common identifier are spanned to avoid treating
    an intron as RNA-free. Explicit circular-origin features retain their parts.
    Unknown strand is preserved and blocks both polarities where overlapping.
    """
    groups = defaultdict(list)
    for i, f in enumerate(ref.get('features', [])):
        if str(f.get('type', '')).lower() not in TRANSCRIBED_TYPES:
            continue
        ident = str(f.get('feature_id') or f.get('label') or f'feature_{i}')
        st = f.get('strand') if f.get('strand') in ('+', '-') else '?'
        # Type is included so gene/CDS entries with different semantics do not
        # turn into one invented transcript.
        key = (ident, str(f.get('type', '')).lower(), st)
        groups[key].append(f)
    out = []
    for (ident, kind, st), features in groups.items():
        groups_to_span = [[f] for f in features] if any(f.get('wraps_origin') for f in features) else [features]
        for fs in groups_to_span:
            a = max(1, min(int(f['start']) for f in fs))
            b = min(len(ref['sequence']), max(int(f['end']) for f in fs))
            if b < a:
                continue
            out.append(dict(start=a, end=b, strand=st, name=str(fs[0].get('label') or ident),
                            feature_id=ident, source=str(fs[0].get('source', 'input annotation')),
                            kind=kind, scope='gene/transcript span; introns included where spanned'))
    # Preserve semantic records but eliminate identical spans/strand from redundant annotations.
    seen = set(); unique = []
    for item in out:
        key = (item['start'], item['end'], item['strand'])
        if key not in seen:
            unique.append(item); seen.add(key)
    return sorted(unique, key=lambda x: (x['start'], x['end'], x['strand']))


def _manual_intervals(ref: dict, rows: list[dict]) -> list[dict]:
    out = []
    for n, row in enumerate(rows, 1):
        a, b = row.get('start'), row.get('end')
        if isinstance(a, bool) or isinstance(b, bool) or not isinstance(a, int) or not isinstance(b, int):
            raise ValueError(f'{ref["id"]}: transcript coordinates must be whole numbers.')
        if not 1 <= a <= b <= len(ref['sequence']):
            raise ValueError(f'{ref["id"]}: transcript {n} extends outside the reference.')
        strand = row.get('strand')
        if strand not in ('+', '-', 'both'):
            raise ValueError(f'{ref["id"]}: choose a transcription direction for transcript {n}.')
        for st in ('+', '-') if strand == 'both' else (strand,):
            out.append(dict(start=a, end=b, strand=st, name=str(row.get('name') or f'Transcript {n}'),
                            feature_id=f'manual:{n}:{st}', source='user-confirmed transcription interval',
                            kind='transcript', scope='user-supplied transcript span; includes introns'))
    return out


def prepare_transcription(refs: list[dict], declarations: dict | None, enabled: bool = True):
    """Validate the declared/annotated RNA model and make local RNA-sense references.

    Raw/construct inputs require an affirmative direction review. A successfully
    fetched empty genomic gene annotation can be used, but is not expression data.
    No feature label, FASTA name or H3K27ac track is used to infer direction.
    """
    declarations = declarations or {}
    if not isinstance(declarations, dict):
        raise ValueError('transcription must be a per-reference object.')
    if not enabled:
        return {}, [], [dict(reference=r['id'], status='GENE_SPANS_NOT_USED_FOR_EXCLUSION', intervals=[],
                            loaded_feature_count=len(r.get('features',[])),
                            note='Annotations remain visible; no whole-span masking or inferred unspliced RNA was generated. Supplied RNA sequences and requested RNA database searches are independent.') for r in refs]
    model = {}; summaries = []; rna_refs = []; seen_rna = set()
    for ref in refs:
        rid = ref['id']; dec = declarations.get(rid, {})
        mapped = bool(ref.get('mapping', {}).get('verified')) if ref.get('mapping') else False
        mode = dec.get('mode', 'annotations' if mapped else 'unknown')
        if mode not in MODES:
            raise ValueError(f'{rid}: unknown transcription mode {mode!r}.')
        annotated = feature_intervals(ref)
        metadata = ref.get('gene_annotation') or {}
        if mode == 'unknown':
            raise ValueError(f'{rid}: gene direction has not been confirmed. In step 1 choose forward, reverse, '
                             'both, reviewed annotations, manual transcript intervals, or explicitly no known transcription. '
                             'DNA targeting is retained; this question is about avoiding its RNA.')
        if not mapped and dec.get('confirmed') is not True:
            raise ValueError(f'{rid}: confirm the construct transcription directions in step 1 before designing.')
        if mapped and mode != 'annotations' and dec.get('confirmed') is not True:
            raise ValueError(f'{rid}: a manual genomic transcription declaration requires confirmation.')
        if mode == 'annotations':
            if mapped and metadata.get('status') not in ('FETCHED', 'USER_VERIFIED_COMPLETE'):
                raise ValueError(f'{rid}: fetch gene/transcript annotations first (step 1). Missing annotations are not RNA-free DNA.')
            if not mapped and not annotated:
                raise ValueError(f'{rid}: this file has no usable gene/transcript annotations. Enter the gene direction or manual intervals.')
        intervals = list(annotated)
        if mode in ('forward', 'reverse', 'both'):
            # Conservatively cover the whole reference when only direction, not
            # transcript boundaries, is known. This can lower probe yield.
            for st in ('+', '-') if mode == 'both' else ('+',) if mode == 'forward' else ('-',):
                intervals.append(dict(start=1, end=len(ref['sequence']), strand=st, name='Declared transcription '+st,
                    feature_id=f'declared:{rid}:{st}', source='user-confirmed whole-reference direction', kind='transcript',
                    scope='whole supplied sequence is conservatively treated as transcribed'))
        if mode == 'manual':
            rows = dec.get('intervals', [])
            if not isinstance(rows, list) or not rows:
                raise ValueError(f'{rid}: manual mode needs at least one transcript interval.')
            intervals += _manual_intervals(ref, rows)
        if mode == 'none' and annotated:
            raise ValueError(f'{rid}: annotated gene/transcript spans are present; they cannot be erased by selecting no transcription. '
                             'Use reviewed annotations or correct the source file.')
        if mode == 'none' and dec.get('confirmed') is not True:
            raise ValueError(f'{rid}: confirm the explicit no-known-transcription declaration.')
        # Every record remains attributable. Identical model spans can be merged.
        unique = {}; 
        for it in intervals:
            unique.setdefault((it['start'], it['end'], it['strand']), it)
        intervals = sorted(unique.values(), key=lambda x:(x['start'],x['end'],x['strand']))
        model[rid] = intervals
        unknown = sum(it['strand'] == '?' for it in intervals)
        summaries.append(dict(reference=rid, mode=mode, status='DIRECTION_MODEL_APPLIED', confirmed=dec.get('confirmed', mapped),
            interval_count=len(intervals), unknown_strand_intervals=unknown, intervals=intervals,
            annotation=metadata, note='RNA-avoidance relative to supplied annotations/declarations, not an RNA-free guarantee. '
            'Opposite-strand overlaps are skipped; unspliced local gene spans are screened as RNA-sense sequences.'))
        for it in intervals:
            fragment = ref['sequence'][it['start']-1:it['end']]
            for st in ('+', '-') if it['strand'] == '?' else (it['strand'],):
                seq = fragment if st == '+' else rc(fragment)
                if seq in seen_rna: continue
                seen_rna.add(seq)
                rna_refs.append(dict(id='__rna_span_'+sha(rid+'|'+str(it['start'])+'|'+str(it['end'])+'|'+st)[:16],
                    sequence=seq, length=len(seq), features=[], mapping=None,
                    source_reference=rid, source_start=it['start'], source_end=it['end'], transcription_strand=st,
                    description=it['name']+' | RNA-sense, unspliced local span (not complete mature transcript)',
                    inferred_gene_span=True))
        for record in ref.get('annotated_rna_sequences', []):
            seq=normalize(str(record['sequence']).upper().replace('U','T'))
            if seq in seen_rna: continue
            seen_rna.add(seq)
            rna_refs.append(dict(id='__rna_feature_'+sha(rid+'|'+str(record['id']))[:16],sequence=seq,length=len(seq),
                features=[],mapping=None,source_reference=rid,description=str(record.get('name','annotated RNA')),
                source_scope=record.get('scope','annotated RNA-sense sequence'),inferred_gene_span=False))
    return model, rna_refs, summaries


def permitted_strands(intervals: list[dict], start: int, end: int, requested: str = 'auto') -> list[str]:
    """Forbid the DNA strand with the SAME sign as each overlapping RNA sense.

    This ensures the ordered arm is sense to that transcript, not its complement.
    Any overlap is treated conservatively; opposite-strand genes may exclude both.
    Outside the known intervals, auto chooses bottom only; this is a convention,
    not proof that intergenic sequence cannot be transcribed.
    """
    forbidden = set()
    for it in intervals:
        if it['start'] <= end and it['end'] >= start:
            forbidden.update(('+', '-') if it['strand'] == '?' else (it['strand'],))
    allowed = [s for s in ('-', '+') if s not in forbidden]
    if requested != 'auto':
        return [requested] if requested in allowed else []
    if not forbidden: return ['-']
    return allowed


def candidate_rna_safety(sites: list[dict], model: dict) -> dict:
    evidence = []; unsafe = []
    for site in sites:
        intervals = model.get(site['reference'], [])
        overlap = [it for it in intervals if it['start'] <= site['end'] and it['end'] >= site['start']]
        ok = bool(permitted_strands(intervals, site['start'], site['end'], site['strand']))
        item = dict(reference=site['reference'], start=site['start'], end=site['end'],
                    bound_strand=site['strand'], genes=[it['name'] for it in overlap],
                    transcription_strands=sorted(set(it['strand'] for it in overlap)), safe_direction=ok)
        evidence.append(item)
        if not ok: unsafe.append(item)
    return dict(status='DIRECTION_COMPATIBLE' if not unsafe else 'RNA_DIRECTION_CONFLICT',
                sites=evidence, conflicts=unsafe, note='Sense-arm design relative to declared/annotated RNA; independent sequence screening also required.')
