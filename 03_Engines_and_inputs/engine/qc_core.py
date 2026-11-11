"""Input-preserving local QC and adapters; see METHODS.md for scope/limitations.

Coordinates in JSON/CSV are 1-based inclusive. Input arms are synthesized 5'->3'.
Matching a forward reference means the physically bound DNA strand is MINUS.
A target allowlist is separate from the search references, including repeat sites.
"""
from __future__ import annotations
import csv
import hashlib
import json
import math
import re
import shutil
import subprocess
from collections import defaultdict
from dataclasses import asdict, dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np
from Bio.SeqUtils import MeltingTemp as mt
from oligominer_thermo_qc import PROFILES, calculate_tm
from probe_binding_viewer import reverse_complement as rc

VERSION = '2.0.1-setup'
DNA = set('ACGT')

@dataclass(frozen=True)
class Settings:
    minimum_gap: int = 5
    max_mismatches: int = 3
    exact_kmer_length: int = 18
    stem_min: int = 6
    stem_review_length: int = 12
    hairpin_review_length: int = 8
    hairpin_min_loop: int = 3
    tm_low: float = 42.0
    tm_high: float = 50.0
    folding_temperature: float = 42.0
    sodium_M: float = 0.390
    max_full_local_bases: int = 2_000_000


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_csv(path: Path, rows: list[dict], columns: list[str] | None = None) -> None:
    if columns is None:
        columns = list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w', encoding='utf-8', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=columns or ['status'], extrasaction='ignore')
        w.writeheader()
        for row in rows:
            clean = {}
            for k in columns or ['status']:
                v = row.get(k, '')
                if isinstance(v, (dict, list, tuple)):
                    v = json.dumps(v, separators=(',', ':'))
                if isinstance(v, str) and v.startswith(('=', '+', '@')):
                    v = "'" + v
                clean[k] = v
            w.writerow(clean)


def load_allowlist(path: Path | None, data: dict) -> list[dict]:
    """Default: each probe's assigned region and observed exact bound strand.

    Full same-strand assigned intervals are allowed, NOT all database matches.
    With no known exact site, default to bottom strand, flagged by mapper as absent.
    User allowlist CSV: Probe_Name or Pool or Readout,Reference,Start,End,Bound_strand.
    Wildcard '*' is supported only in selector columns, never reference/strand.
    Explicit rows extend (rather than replace) the assigned local-region allowance.
    """
    rules = []
    refs = {r['id']: r for r in data['references']}
    for p in data['probes']:
        rid = p['expected_region']
        if rid not in refs:
            raise ValueError(f"{p['name']}: Region must name one of the target FASTA records")
        strands = {data['hits'][hid]['bound_strand'] for hid in p['hits']
                   if data['hits'][hid]['reference_id'] == rid} or {'-'}
        for strand in sorted(strands):
            rules.append(dict(probe=p['name'], reference=rid, start=1,
                              end=refs[rid]['length'], strand=strand, origin='assigned_target_region'))
    if path:
        with path.open(encoding='utf-8-sig', newline='') as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                rid = row.get('Reference', '').strip()
                strand = row.get('Bound_strand', '').strip()
                start, end = int(row['Start']), int(row['End'])
                if not rid or strand not in ['+', '-'] or not 1 <= start <= end:
                    raise ValueError('Invalid allowlist interval/strand')
                matches = [p for p in data['probes'] if
                    (row.get('Probe_Name', '*') or '*') in ('*', p['name']) and
                    (row.get('Pool', '*') or '*') in ('*', p['pool']) and
                    (row.get('Readout', '*') or '*') in ('*', p['readout'])]
                if not matches:
                    raise ValueError(f'Allowlist row matches no probe: {row}')
                rules.extend(dict(probe=p['name'], reference=rid, start=start,
                                  end=end, strand=strand, origin='explicit_allowlist') for p in matches)
    return rules


def allowed(rules: list[dict], name: str, reference: str, start: int,
            end: int, bound_strand: str) -> bool:
    return any(r['probe'] == name and r['reference'] == reference and
               r['strand'] == bound_strand and r['start'] <= start <= end <= r['end'] for r in rules)


def initialize_qc(data: dict, settings: Settings, active: str) -> None:
    if active not in PROFILES:
        raise ValueError('Unknown temperature profile')
    for p in data['probes']:
        profiles = {}
        for name, prof in PROFILES.items():
            vals = calculate_tm(p['arm'], prof)
            profiles[name] = {**vals, 'in_window': settings.tm_low < vals['adjusted_screen_C'] < settings.tm_high}
        seq = p['arm']
        longest = max(len(m.group()) for m in re.finditer('A+|C+|G+|T+', seq))
        seq_checks = {'length_36_90': 36 <= len(seq) <= 90,
                      'GC_at_most_76': p['gc_percent'] <= 76,
                      'no_homopolymer_ge6': longest < 6,
                      'no_dinucleotide_8units': not bool(re.search(r'([ACGT]{2})\1{7}', seq)),
                      'longest_homopolymer': longest}
        p['qc'] = dict(thermo=profiles, composition=seq_checks, off_targets=[], local_expected=[],
                       hairpin_stems=[], partners=[], spacing_flags=[], full_fold=None,
                       full_fold_status='NOT_RUN')
    data['qc'] = {'version': VERSION, 'active_profile': active,
                  'profiles': {name: asdict(prof) for name, prof in PROFILES.items()},
                  'settings': asdict(settings), 'allowlist': [], 'checks': [], 'sources': [],
                  'no_sequences_changed': True, 'sequence_sha256': hashlib.sha256(
                      '\n'.join(p['name'] + '\t' + p['full'] for p in data['probes']).encode()).hexdigest(),
                  'full_folding_scope': 'NOT_RUN; exact-stem geometry only',
                  'no_experimental_validation': True}


def flank_spacing(data: dict, settings: Settings) -> None:
    refs = {r['id']: r['sequence'] for r in data['references']}
    for h in data['hits']:
        p = data['probes'][h['probe_index']]
        seq = refs[h['reference_id']]
        prefix, suffix = p['prefix'], p['suffix']
        if h['bound_strand'] == '+':
            prefix, suffix = rc(suffix), rc(prefix)
        # Identity to the top reference is equivalent to pairing with bottom DNA.
        left = right = 0
        start, end = h['start'] - 1, h['end']
        while left < len(prefix) and start - left - 1 >= 0 and prefix[-1-left] == seq[start-left-1]:
            left += 1
        while right < len(suffix) and end + right < len(seq) and suffix[right] == seq[end+right]:
            right += 1
        h.update(flank_left=left, flank_right=right,
                 extended_start=h['start'] - left, extended_end=h['end'] + right)
    for item in data['spacing_audit']:
        a, b = data['hits'][item['left_hit_id']], data['hits'][item['right_hit_id']]
        ext_gap = b['extended_start'] - a['extended_end'] - 1
        item['potential_extended_gap'] = ext_gap
        item['potential_flank_warning'] = ext_gap < settings.minimum_gap
        item['same_strand'] = a['bound_strand'] == b['bound_strand']
        if item['potential_flank_warning']:
            for h, other in [(a, b), (b, a)]:
                data['probes'][h['probe_index']]['qc']['spacing_flags'].append({
                    'site_id': h['id'], 'partner_index': other['probe_index'],
                    'partner_name': data['probes'][other['probe_index']]['name'],
                    'reference': h['reference_id'], 'start': h['start'], 'end': h['end'],
                    'arm_gap': item['length'], 'potential_gap': ext_gap,
                    'same_strand': item['same_strand']})


def full_hamming_matches(query: str, reference: str, max_mismatches: int) -> list[tuple[int, int]]:
    """Exhaustive full-length ungapped matches, including overlapping positions.
    Ambiguous reference bases never count as matches. No seed/aligner hit cap.
    """
    n = len(query)
    if len(reference) < n:
        return []
    q = np.frombuffer(query.encode('ascii'), dtype=np.uint8)
    ref = np.frombuffer(reference.encode('ascii'), dtype=np.uint8)
    result = []
    # Chunk to keep memory bounded, even for larger local construct references.
    chunk_size = 20000
    for base in range(0, len(reference)-n+1, chunk_size):
        raw = ref[base:min(len(ref), base+chunk_size+n-1)]
        windows = np.lib.stride_tricks.sliding_window_view(raw, n)
        nm = (windows != q).sum(axis=1)
        result.extend((base + int(i), int(nm[i])) for i in np.flatnonzero(nm <= max_mismatches))
    return result


def maximal_exact_segments(query: str, reference: str, k: int) -> list[tuple[int, int, int]]:
    """All maximal exact tracts >=k, returns query start, ref start, length (0-based)."""
    if min(len(query), len(reference)) < k:
        return []
    idx = defaultdict(list)
    for i in range(len(reference)-k+1):
        word = reference[i:i+k]
        if not set(word)-DNA:
            idx[word].append(i)
    out = set()
    for i in range(len(query)-k+1):
        for j in idx.get(query[i:i+k], []):
            # Emit only at maximal left boundary, then extend right.
            if i and j and query[i-1] == reference[j-1]:
                continue
            n = k
            while i+n < len(query) and j+n < len(reference) and query[i+n] == reference[j+n]:
                n += 1
            out.add((i, j, n))
    return sorted(out)


def local_screen(data: dict, search_refs: list[dict], rules: list[dict],
                 settings: Settings, rna_sense: bool = False) -> list[dict]:
    """Full-arm mismatch screen + exact >=18-bp partial screen, no gapped alignment.

    DNA: both reference orientations, complementarity-based strand assignment.
    RNA: only rc(arm) matches the 5'->3' RNA sequence; U normalized to T by caller.
    full-primary exact partial screen also includes spacer/handle junctions.
    """
    if sum(len(r['sequence']) for r in search_refs) > settings.max_full_local_bases:
        raise ValueError('Local references exceed configured size; use --blast-db for a genome/transcriptome index')
    all_rows = []
    for p in data['probes']:
        for r in search_refs:
            reference, rid = r['sequence'], r['id']
            full_found = []
            dirs = [('+', rc(p['arm']))] if rna_sense else [('-', p['arm']), ('+', rc(p['arm']))]
            for strand, query in dirs:
                for pos, mismatch in full_hamming_matches(query, reference, settings.max_mismatches):
                    ok = False if rna_sense else allowed(rules, p['name'], rid, pos+1, pos+len(query), strand)
                    row = dict(probe=p['name'], probe_index=p['index'], reference=rid,
                        start=pos+1, end=pos+len(query), bound_strand='RNA' if rna_sense else strand,
                        kind='full_arm_ungapped', mismatches=mismatch, alignment_length=len(query),
                        query_start=1, query_end=len(query), query_coverage=100.0,
                        identity=100.0*(len(query)-mismatch)/len(query), allowed=ok,
                        source=r.get('source_label', 'local'), full_primary=False,
                        qseq=p['arm'],sseq=reference[pos:pos+len(query)] if strand=='-' else rc(reference[pos:pos+len(query)]))
                    full_found.append(row)
                    all_rows.append(row)
            seq = p['full'] or p['arm']
            dirs = [('+', rc(seq))] if rna_sense else [('-', seq), ('+', rc(seq))]
            for strand, query in dirs:
                for qi, pos, n in maximal_exact_segments(query, reference, settings.exact_kmer_length):
                    # Do not duplicate a partial hit represented by a full-arm hit.
                    if any(h['bound_strand'] == ('RNA' if rna_sense else strand) and
                           h['start'] <= pos+1 <= pos+n <= h['end'] for h in full_found):
                        continue
                    ok = False if rna_sense else allowed(rules, p['name'], rid, pos+1, pos+n, strand)
                    start_q = qi+1 if strand == '-' else len(seq)-qi-n+1
                    row = dict(probe=p['name'], probe_index=p['index'], reference=rid,
                        start=pos+1, end=pos+n, bound_strand='RNA' if rna_sense else strand,
                        kind='exact_partial_full_primary', mismatches=0, alignment_length=n,
                        query_start=start_q, query_end=start_q+n-1,
                        query_coverage=round(100.0*n/len(seq), 3), identity=100.0, allowed=ok,
                        source=r.get('source_label', 'local'), full_primary=True,
                        qseq=seq[start_q-1:start_q+n-1],sseq=reference[pos:pos+n] if strand=='-' else rc(reference[pos:pos+n]))
                    all_rows.append(row)
    for row in all_rows:
        q = data['probes'][row['probe_index']]['qc']
        q['local_expected' if row['allowed'] else 'off_targets'].append(row)
    return all_rows


def exact_pair_segments(a: str, b: str, minimum: int) -> list[dict]:
    """Maximal Watson-Crick complementary tracts on full-length primary sequences."""
    out = []
    for i, j, n in maximal_exact_segments(a, rc(b), minimum):
        out.append({'a_start': i+1, 'a_end': i+n,
                    'b_start': len(b)-j-n+1, 'b_end': len(b)-j,
                    'length': n, 'sequence_a': a[i:i+n]})
    return out


@lru_cache(maxsize=10000)
def isolated_stem_tm(seq: str, sodium_mM: float = 390.0) -> float:
    # Both interacting primary species at 2 nM, not a 2/0 target calculation.
    # This is ONLY a short isolated duplex proxy, NOT a full dimer or hairpin Tm.
    temp = mt.Tm_NN(seq, nn_table=mt.DNA_NN3, Na=sodium_mM, Mg=0, K=0, Tris=0,
                   dnac1=2, dnac2=2, dNTPs=0, selfcomp=(seq == rc(seq)), saltcorr=5)
    return float(temp-0.65*50)


def exact_stem_screen(data: dict, settings: Settings) -> list[dict]:
    probes = data['probes']
    pairs = []
    for p in probes:
        seq = p['full'] or p['arm']
        hp = []
        seen = set()
        for hit in exact_pair_segments(seq, seq, settings.stem_min):
            # Truncate overlapping outer tract to respect minimum hairpin loop.
            left, right = hit['a_start'], hit['b_end']
            if left >= right:
                continue
            max_n = (right-left+1-settings.hairpin_min_loop)//2
            n = min(hit['length'], max_n)
            if n < settings.stem_min:
                continue
            arm_right = right-n+1
            loop = arm_right-(left+n)
            coords = (left,left+n-1,arm_right,right)
            if coords in seen:
                continue
            seen.add(coords)
            hp.append(dict(a_start=left,a_end=left+n-1,b_start=arm_right,b_end=right,
                           length=n,loop_bases=loop,sequence_a=seq[left-1:left-1+n]))
        hp.sort(key=lambda x: (-x['length'], x['loop_bases'], x['a_start']))
        p['qc']['hairpin_stems'] = hp
        p['qc']['max_hairpin_stem'] = max((h['length'] for h in hp), default=0)
    for i, a in enumerate(probes):
        seq_a = a['full'] or a['arm']
        for j in range(i, len(probes)):
            b = probes[j]
            segs = exact_pair_segments(seq_a, b['full'] or b['arm'], settings.stem_min)
            for s in segs:
                s['isolated_duplex_Tm_proxy_C'] = round(isolated_stem_tm(s['sequence_a']), 4)
            segs.sort(key=lambda h: (-h['length'], -h['isolated_duplex_Tm_proxy_C']))
            best = segs[0] if segs else None
            entry = dict(a_index=i,b_index=j,a_name=a['name'],b_name=b['name'],
                         is_homodimer=i==j, max_exact_stem=best['length'] if best else 0,
                         max_isolated_duplex_Tm_proxy_C=max((s['isolated_duplex_Tm_proxy_C'] for s in segs),default=None),
                         best_segment=best,segments=segs, exact_review=bool(best and best['length']>=settings.stem_review_length),
                         full_cofold=None)
            pairs.append(entry)
            if best:
                a['qc']['partners'].append(dict(partner_index=j,partner_name=b['name'],
                                                max_stem=best['length'],pair_index=len(pairs)-1))
                if j != i:
                    b['qc']['partners'].append(dict(partner_index=i,partner_name=a['name'],
                                                    max_stem=best['length'],pair_index=len(pairs)-1))
    for p in probes:
        p['qc']['partners'].sort(key=lambda q:(-q['max_stem'],q['partner_index']))
        p['qc']['max_dimer_stem'] = max((r['max_stem'] for r in p['qc']['partners']),default=0)
    data['qc']['pair_interactions'] = pairs
    return pairs


def _dot_pairs(structure: str) -> list[tuple[int,int]]:
    stack, pairs = [], []
    for i,c in enumerate(structure.replace('&','')):
        if c == '(':
            stack.append(i)
        elif c == ')':
            if not stack:
                raise ValueError('Unbalanced folding structure')
            pairs.append((stack.pop(), i))
    if stack:
        raise ValueError('Unbalanced folding structure')
    return sorted(pairs)


def run_vienna(data: dict, settings: Settings, mode: str = 'auto') -> str:
    """Optional full-length monomer MFE/PF and all-pair MFE cofold adapter.

    DNA Mathews 2004, aqueous Na=0.390 M, T=42 C. No formamide model,
    no pseudoknots/G-quadruplexes, no all-pool concentration equilibrium.
    Does NOT report a pass/fail probability or a concentration-corrected dimer Tm.
    """
    if mode == 'stems':
        return 'NOT_RUN: exact-stem mode requested'
    try:
        import RNA
    except ImportError:
        if mode == 'vienna':
            raise RuntimeError('Full folding requested but ViennaRNA is missing. Install ViennaRNA==2.7.2.')
        return 'NOT_RUN: ViennaRNA unavailable; exact-stem checks were run instead'
    if not hasattr(RNA, 'params_load_DNA_Mathews2004'):
        raise RuntimeError('ViennaRNA lacks DNA Mathews2004 parameters; refusing RNA defaults')
    RNA.params_load_DNA_Mathews2004()
    md = RNA.md()
    md.temperature = settings.folding_temperature
    md.salt = settings.sodium_M
    md.dangles = 2
    md.noGU = 1
    md.gquad = 0
    for fold_index, p in enumerate(data['probes'], 1):
        if fold_index == 1 or fold_index % 25 == 0 or fold_index == len(data['probes']):
            print(f'  DNA monomer folding {fold_index}/{len(data["probes"])}', flush=True)
        seq = p['full'] or p['arm']
        fc = RNA.fold_compound(seq, md)
        structure, energy = fc.mfe()
        fc.exp_params_rescale(energy)
        _, ensemble = fc.pf()
        matrix = fc.bpp()
        unpaired = []
        for i in range(1,len(seq)+1):
            unpaired.append(max(0.0, min(1.0, 1.0 - sum(matrix[min(i,j)][max(i,j)]
                                               for j in range(1,len(seq)+1) if i != j))))
        a0 = len(p['prefix'])
        a1 = a0 + len(p['arm'])
        p['qc']['full_fold'] = dict(structure=structure,mfe_kcal_mol=float(energy),
            ensemble_kcal_mol=float(ensemble),mean_arm_unpaired=sum(unpaired[a0:a1])/len(p['arm']),
            base_unpaired=unpaired,paired_positions=_dot_pairs(structure),engine='ViennaRNA '+RNA.__version__)
        p['qc']['full_fold_status'] = 'RUN_AQUEOUS_MODEL'
    for pair_index, pair in enumerate(data['qc']['pair_interactions'], 1):
        if pair_index == 1 or pair_index % 500 == 0 or pair_index == len(data['qc']['pair_interactions']):
            print(f'  DNA pair cofolding {pair_index}/{len(data["qc"]["pair_interactions"])}', flush=True)
        a,b = [data['probes'][pair[k]] for k in ['a_index','b_index']]
        sa,sb = a['full'] or a['arm'], b['full'] or b['arm']
        fc = RNA.fold_compound(sa+'&'+sb,md)
        structure, energy = fc.mfe()
        inter = [(i,j-len(sa)) for i,j in _dot_pairs(structure) if i<len(sa)<=j]
        delta = energy-a['qc']['full_fold']['mfe_kcal_mol']-b['qc']['full_fold']['mfe_kcal_mol']
        pair['full_cofold'] = dict(structure=structure,mfe_kcal_mol=float(energy),
            mfe_gain_vs_separate_kcal_mol=float(delta),intermolecular_pairs=inter)
    return 'RUN_AQUEOUS_MODEL: ViennaRNA '+RNA.__version__+'; no formamide correction or pool equilibrium'


BLAST_FIELDS = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend',
                'sstart','send','evalue','bitscore','qlen','qseq','sseq']


def parse_blast(path: Path, data: dict, rules: list[dict], source: str,
                rna: bool = False) -> list[dict]:
    """Parse tabular outfmt 6. Preserve weak HSPs too; flag, do not filter probes.
    DNA query = synthesized arm. Query plus orientation is forced by runner.
    RNA subject is transcript sense; only reverse-strand HSPs count as binding.
    """
    pmap = {p['name']:p for p in data['probes']}
    result = []
    with path.open(encoding='utf-8') as fh:
        for num,line in enumerate(fh,1):
            if not line.strip() or line.startswith('#'):
                continue
            vals = line.rstrip('\r\n').split('\t')
            if len(vals)!=len(BLAST_FIELDS):
                raise ValueError(f'BLAST row {num} must contain {len(BLAST_FIELDS)} fields; see HOW_TO_RUN.txt')
            raw = dict(zip(BLAST_FIELDS,vals))
            if raw['qseqid'] not in pmap:
                raise ValueError('BLAST includes unknown probe: '+raw['qseqid'])
            p=pmap[raw['qseqid']]
            q1,q2,s1,s2=[int(raw[k]) for k in ['qstart','qend','sstart','send']]
            if q1>q2:
                raise ValueError('BLAST import requires query-plus coordinates. Generate output with supplied runner')
            if int(raw['qlen'])!=len(p['arm']):
                raise ValueError('BLAST query length differs from current probe arm')
            if p['arm'][q1-1:q2] != raw['qseq'].replace('-','').upper():
                raise ValueError('BLAST query sequence differs from current probe input')
            refstrand='+' if s1<=s2 else '-'
            if rna and refstrand!='-':
                continue
            bound='-' if refstrand=='+' else '+'
            identical=[aa==bb and aa!='-' for aa,bb in zip(raw['qseq'].upper(),raw['sseq'].upper())]
            runs=''.join('1' if x else '0' for x in identical).split('0')
            exact=max(map(len,runs), default=0)
            n=int(raw['length']); ident=float(raw['pident'])
            qcov=100.0*(q2-q1+1)/len(p['arm'])
            review=exact>=18 or (ident>=80 and (qcov>=70 or n>=30))
            start,end=sorted((s1,s2))
            ok=False if rna else allowed(rules,p['name'],raw['sseqid'],start,end,bound)
            row=dict(probe=p['name'],probe_index=p['index'],reference=raw['sseqid'],start=start,end=end,
                bound_strand='RNA' if rna else bound,kind='blast_local_alignment',mismatches=int(raw['mismatch']),
                gap_opens=int(raw['gapopen']),alignment_length=n,query_start=q1,query_end=q2,
                query_coverage=qcov,identity=ident,allowed=ok,source=source,exact_tract=exact,
                review=review,evalue=float(raw['evalue']),bitscore=float(raw['bitscore']),
                qseq=raw['qseq'],sseq=raw['sseq'],full_primary=False)
            result.append(row)
            p['qc']['local_expected' if ok else 'off_targets'].append(row)
    return result


def run_blast(db: str, query: Path, out: Path, threads: int=4, executable: str='blastn') -> dict:
    binary=shutil.which(executable)
    if not binary:
        raise RuntimeError('blastn is unavailable. Install NCBI BLAST+ or supply --blast-tsv plus its run metadata')
    cmd=[binary,'-query',str(query),'-db',db,'-task','blastn-short','-strand','both',
         '-dust','no','-soft_masking','false','-word_size','7','-evalue','1000',
         '-reward','1','-penalty','-3','-gapopen','5','-gapextend','2',
         '-max_target_seqs','100000','-num_threads',str(threads),
         '-outfmt','6 '+' '.join(BLAST_FIELDS),'-out',str(out)]
    # BLASTN normalizes its query to plus and encodes minus orientation on subject.
    proc=subprocess.run(cmd,text=True,capture_output=True)
    if proc.returncode:
        raise RuntimeError('BLAST failed: '+proc.stderr[-5000:])
    version=subprocess.run([binary,'-version'],text=True,capture_output=True).stdout.strip()
    return dict(command=cmd,version=version,stderr=proc.stderr,database=db,
                output_sha256=sha256(out),heuristic_search=True,max_target_seqs=100000,
                cap_note='100000 subject-sequence limit; absence of a cap warning is not proof of exhaustive alignments')


def finish(data: dict, settings: Settings) -> None:
    for p in data['probes']:
        q=p['qc']
        flagged=[h for h in q['off_targets'] if h.get('review',True)]
        q['off_target_review_count']=len(flagged)
        q['local_status']='REVIEW' if flagged else 'NO_FLAGS_IN_SEARCHED_REFERENCES'
        q['exact_stem_review']=q['max_dimer_stem']>=settings.stem_review_length or q['max_hairpin_stem']>=settings.hairpin_review_length
        q['not_full_validation']=True
    data['qc']['counts']={
        'probes':len(data['probes']),'exact_sites':len(data['hits']),
        'tm_in_window':{name:sum(p['qc']['thermo'][name]['in_window'] for p in data['probes']) for name in PROFILES},
        'local_offtarget_flagged_probes':sum(p['qc']['off_target_review_count']>0 for p in data['probes']),
        'exact_pair_comparisons':len(data['qc']['pair_interactions']),
        'exact_pair_flags':sum(p['exact_review'] for p in data['qc']['pair_interactions']),
        'hairpin_geometry_flags':sum(p['qc']['max_hairpin_stem']>=settings.hairpin_review_length for p in data['probes']),
        'flank_spacing_pairs':sum(r['potential_flank_warning'] for r in data['spacing_audit']),
    }
