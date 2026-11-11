"""Reference handling. Public coordinates are always explicit; local display is 1-based inclusive."""
from __future__ import annotations
import base64
import hashlib
import io
import json
import math
import re
from functools import lru_cache
from pathlib import Path
from typing import Any
from Bio import SeqIO

ALPHABET = set('ACGTRYSWKMBDHVN')
COMPLEMENT = str.maketrans('ACGTRYSWKMBDHVN', 'TGCAYRSWMKVHDBN')
# GRCh38 primary chromosome RefSeq accessions; versions are intentional.
CHROMOSOMES = dict(zip(['chr'+str(i) for i in range(1,23)]+['chrX','chrY','chrM'],
 ['NC_000001.11','NC_000002.12','NC_000003.12','NC_000004.12','NC_000005.10','NC_000006.12',
  'NC_000007.14','NC_000008.11','NC_000009.12','NC_000010.11','NC_000011.10','NC_000012.12',
  'NC_000013.11','NC_000014.9','NC_000015.10','NC_000016.10','NC_000017.11','NC_000018.10',
  'NC_000019.10','NC_000020.11','NC_000021.9','NC_000022.11','NC_000023.11','NC_000024.10','NC_012920.1']))

def rc(s: str) -> str:
    return s.translate(COMPLEMENT)[::-1]

def sha(s: str | bytes) -> str:
    return hashlib.sha256(s.encode('utf-8') if isinstance(s,str) else s).hexdigest()

def canonical_ref(s: str) -> str:
    """Aliases for the same primary GRCh38 sequence only; do not merge alt/patch loci."""
    return _canonical_ref_cached(s.strip())

@lru_cache(maxsize=4096)
def _canonical_ref_cached(s: str) -> str:
    # CHROMOSOMES is the fixed GRCh38 accession table above. Do not cache any
    # coordinate, strand, intended-site allowance, or caller-owned mapping.
    if s.startswith('ref|') and s.endswith('|'): s=s[4:-1]
    if s in CHROMOSOMES: return CHROMOSOMES[s]
    if 'chr'+s in CHROMOSOMES: return CHROMOSOMES['chr'+s]
    return s

def normalize(s: str, strict=False) -> str:
    s=''.join(s.split()).upper()
    if not s or set(s)-(set('ACGT') if strict else ALPHABET):
        raise ValueError('Expected DNA sequence (A/C/G/T or IUPAC reference bases). No annotations or numbers in sequence text.')
    return s

def parse_fasta(text: str) -> list[dict]:
    text=text.lstrip('\ufeff').strip()
    if not text: raise ValueError('No DNA sequence provided.')
    if not text.startswith('>'): text='>sequence\n'+text
    records=[]
    for record in SeqIO.parse(io.StringIO(text),'fasta'):
        seq=normalize(str(record.seq))
        rid=record.id
        if not re.fullmatch(r'[A-Za-z0-9_.:-]+',rid): raise ValueError('Use simple FASTA IDs: letters, numbers, underscore, dot, colon, or dash.')
        records.append(dict(id=rid,description=record.description,sequence=seq,length=len(seq),sha256=sha(seq),features=[],mapping=None))
    if not records or len({r['id'] for r in records})!=len(records): raise ValueError('FASTA IDs must be nonempty and unique.')
    return records

def parse_annotated(blob: bytes, filename: str) -> list[dict]:
    fmt='snapgene' if filename.lower().endswith('.dna') else 'genbank'
    stream=io.BytesIO(blob) if fmt=='snapgene' else io.StringIO(blob.decode('utf-8-sig'))
    records=[]
    for i,rec in enumerate(SeqIO.parse(stream,fmt)):
        seq=normalize(str(rec.seq)); rid=re.sub(r'[^A-Za-z0-9_.:-]','_',rec.id) or f'seq{i+1}'
        feats=[]; annotated_rna=[]
        for j,ft in enumerate(rec.features):
            if ft.type=='source' or ft.location is None: continue
            quals=ft.qualifiers
            label=next((str(quals[k][0]) for k in ('label','gene','locus_tag','product','note') if quals.get(k)),ft.type)
            # Preserve annotated spliced mRNA/ncRNA/CDS sequences as additional
            # RNA-sense references. A CDS is labeled partial, not a full transcript.
            if ft.type.lower() in ('mrna','transcript','ncrna','rrna','trna','misc_rna','cds') and ft.location.strand in (-1,1):
                rna_seq=str(ft.extract(rec.seq)).upper()
                if rna_seq:
                    annotated_rna.append(dict(id=f'{rid}:{j}',name=label,sequence=rna_seq,
                        source=filename,type=ft.type,scope='annotated spliced feature; CDS-only records omit untranslated sequence'))
            # Store multipart joins separately, not as a misleading uninterrupted exon.
            for part in ft.location.parts:
                a,b=int(part.start)+1,int(part.end)
                if b>=a:
                    feats.append(dict(type=ft.type,label=label,start=a,end=b,strand='-' if part.strand==-1 else '+' if part.strand==1 else '?',
                        source=filename,feature_id=f'{rid}:{j}',original_start=a,original_end=b,
                        wraps_origin=bool(rec.annotations.get('topology')=='circular' and len(ft.location.parts)>1 and
                          ((ft.location.strand==1 and int(ft.location.parts[0].start)>int(ft.location.parts[-1].start)) or
                           (ft.location.strand==-1 and int(ft.location.parts[0].start)<int(ft.location.parts[-1].start))))))
        records.append(dict(id=rid,description=rec.description,sequence=seq,length=len(seq),sha256=sha(seq),
            features=feats,mapping=None,annotated_rna_sequences=annotated_rna,topology=rec.annotations.get('topology','unknown')))
    if not records: raise ValueError('No record found in annotated file.')
    if len({r['id'] for r in records})!=len(records): raise ValueError('Annotated records need unique IDs.')
    return records

def resolve_imager_entry(entry: dict) -> dict:
    """Resolve two explicit input conventions to one canonical imager/docking pair.

    All stored sequences are 5'-to-3'. `sequence` ALWAYS means the fluorescent
    imaging oligo, even for a record originally entered as a docking site.
    New input uses input_kind='imager'/'docking' + input_sequence_5to3. Existing
    v1.2 records with only `sequence` retain their original interpretation.
    Redundant sequence fields are checked rather than silently overwritten.
    """
    if not isinstance(entry, dict):
        raise ValueError('An imaging-oligo entry must be an object.')
    if entry.get('sequence_orientation', '5to3') != '5to3':
        raise ValueError('Enter and store each sequence 5-prime to 3-prime; 3-prime-to-5-prime display strings are not accepted.')
    ident = str(entry.get('id', '')).strip().upper()
    if not re.fullmatch('[A-Z0-9_-]+', ident):
        raise ValueError('Imager IDs must be simple unique names.')

    def read_seq(key: str):
        value = entry.get(key)
        if value is None or value == '':
            return None
        if not isinstance(value, str):
            raise ValueError(key + ' must be DNA text written 5-prime to 3-prime.')
        seq = normalize(value, True)
        if not 8 <= len(seq) <= 40:
            raise ValueError('Imager/docking input must be 8-40 DNA bases, 5-prime to 3-prime, without dye notation or TT spacers.')
        return seq

    def agree(first, second, label):
        if first is not None and second is not None and first != second:
            raise ValueError('Conflicting ' + label + ' fields. No sequence was changed; check the input type and orientation.')
        return first if first is not None else second

    imager = agree(read_seq('sequence'), read_seq('imager_sequence_5to3'), 'imaging-oligo')
    docking = agree(read_seq('docking'), read_seq('docking_sequence_5to3'), 'docking-site')
    entered = read_seq('input_sequence_5to3')
    kind = entry.get('input_kind')
    if kind is None:
        if entered is not None:
            raise ValueError('Choose input_kind=imager or docking for input_sequence_5to3. Sequence polarity cannot be guessed.')
        # Legacy `sequence` is a fluorescent imager. Explicit docking-only
        # records are also unambiguous; the generic sequence field never flips.
        kind = 'imager' if imager is not None else 'docking' if docking is not None else None
    if kind not in ('imager', 'docking'):
        raise ValueError('Choose whether the entered sequence is the imaging oligo or its docking/binding site on the primary.')
    if entered is None:
        entered = imager if kind == 'imager' else docking
    if entered is None:
        raise ValueError('Missing entered sequence for the selected input type. Use input_sequence_5to3, not an ambiguously labeled sequence.')
    resolved_imager = entered if kind == 'imager' else rc(entered)
    resolved_docking = rc(entered) if kind == 'imager' else entered
    agree(imager, resolved_imager, 'imaging-oligo')
    agree(docking, resolved_docking, 'docking-site / reverse-complement')
    return dict(id=ident, sequence=resolved_imager, docking=resolved_docking,
        imager_sequence_5to3=resolved_imager, docking_sequence_5to3=resolved_docking,
        input_kind=kind, input_sequence_5to3=entered, sequence_orientation='5to3',
        dye=str(entry.get('dye', 'not specified')),
        source=entry.get('source', 'user-defined'))


def validate_imagers(imagers: list[dict]) -> list[dict]:
    out = [resolve_imager_entry(entry) for entry in imagers]
    if len(set(x['id'] for x in out)) != len(out):
        raise ValueError('Duplicate imager ID.')
    if len(set(x['sequence'] for x in out)) != len(out):
        raise ValueError('Two entries resolve to the same imaging-oligo/docking pair and cannot distinguish bins, even when entered using different input types.')
    return out

def make_bins(refs: list[dict], width: int, tail_policy='trim', assignments=None, bounds=None) -> tuple[list[dict],list[dict]]:
    if not isinstance(width,int) or width<50 or width>1_000_000: raise ValueError('Region width must be an integer from 50 to 1,000,000 bp.')
    if tail_policy not in ('trim','strict'): raise ValueError('Unknown remainder policy.')
    bins=[]; tails=[]; assignments=assignments or {}; bounds=bounds or {}
    for r in refs:
        a,b=bounds.get(r['id'],[1,r['length']])
        a,b=int(a),int(b)
        if not 1<=a<=b<=r['length']: raise ValueError('Design bounds are outside '+r['id'])
        n,rem=divmod(b-a+1,width)
        if n==0: raise ValueError(f'{r["id"]}: interval is shorter than one {width}-bp region.')
        if rem and tail_policy=='strict': raise ValueError(f'{r["id"]}: {rem} bp remain. Choose a divisible interval or explicitly exclude the remainder.')
        for j in range(n):
            start=a+j*width; end=start+width-1; ident=f'{r["id"]}_B{j+1:03d}'
            mapping=r.get('mapping'); gs=ge=None
            if mapping:
                if mapping.get('orientation','+')=='+': gs=mapping['start']+start-1;ge=mapping['start']+end-1
                else: gs=mapping['end']-end+1;ge=mapping['end']-start+1
            bins.append(dict(id=ident,reference=r['id'],number=j+1,start=start,end=end,length=width,
                             center=(start+end)/2,imager=assignments.get(ident,''),genomic_start=gs,genomic_end=ge,
                             chromosome=mapping.get('chromosome') if mapping else None))
        if rem: tails.append(dict(reference=r['id'],start=b-rem+1,end=b,length=rem,reason='Excluded to retain equal region widths'))
    return bins,tails

def clip_features(features: list[dict], start:int, end:int) -> list[dict]:
    return [dict(f,view_start=max(start,f['start']),view_end=min(end,f['end']),
                 clipped=f['start']<start or f['end']>end) for f in features if f['start']<=end and f['end']>=start]

def genomic_interval(ref:dict, start:int,end:int,strand:str) -> dict | None:
    m=ref.get('mapping')
    if not m: return None
    if m.get('orientation','+')=='+': a=m['start']+start-1;b=m['start']+end-1;s=strand
    else: a=m['end']-end+1;b=m['end']-start+1;s='+' if strand=='-' else '-'
    return dict(reference=canonical_ref(m['chromosome']),start=a,end=b,strand=s,assembly=m.get('assembly','GRCh38'))

def interval_union(intervals):
    out=[]
    for a,b in sorted(intervals):
        if out and a<=out[-1][1]+1: out[-1][1]=max(b,out[-1][1])
        else: out.append([a,b])
    return out

def region_metrics(b:dict, candidates:list[dict]) -> dict:
    sites=[h for c in candidates for h in c['sites'] if h['reference']==b['reference'] and b['start']<=h['start']<=h['end']<=b['end']]
    intervals=interval_union([(h['start'],h['end']) for h in sites]);covered=sum(e-a+1 for a,e in intervals)
    gaps=[];last=b['start']-1
    for a,e in intervals:
        if a-last-1>0: gaps.append(a-last-1)
        last=e
    if b['end']-last>0:gaps.append(b['end']-last)
    centers=sorted((s['start']+s['end'])/2 for s in sites)
    hist=[0]*10
    for center in centers:hist[min(9,int((center-b['start'])*10/b['length']))]+=1
    sector_count=max(1,(b['length']+999)//1000);sectors=[0]*sector_count
    for center in centers:sectors[min(sector_count-1,int((center-b['start'])//1000))]+=1
    spacings=[right-left for left,right in zip(centers,centers[1:])]
    centroid=sum(centers)/len(centers) if centers else None
    return dict(probes=len(candidates),exact_sites=len(sites),covered_bases=covered,coverage_percent=100*covered/b['length'],
      longest_gap=max(gaps,default=0),centroid=centroid,centroid_offset=centroid-b['center'] if centroid is not None else None,
      largest_gap_bp=max(gaps,default=0),median_spacing_bp=__import__('statistics').median(spacings) if spacings else None,
      sector_width_bp=1000,sector_count=sector_count,sector_counts=sectors,
      distribution=hist,site_centers=centers,median_tm=__import__('statistics').median(c['tm'] for c in candidates) if candidates else None,
      probe_density_per_kb=len(candidates)*1000/b['length'])
