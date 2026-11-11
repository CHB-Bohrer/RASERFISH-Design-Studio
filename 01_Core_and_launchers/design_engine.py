"""QC-guided candidate selection, not post-hoc flagging.

A transparent, balanced greedy search. It is NOT the complete OligoMiner/LDA
pipeline, an exhaustive combinatorial optimizer, or an experimental certificate.
All selection thresholds are explicit, adjustable heuristics saved with the run.
"""
from __future__ import annotations
import csv
import datetime as dt
import io
import hashlib
import json
import math
import os
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
from collections import Counter,defaultdict,OrderedDict
from contextlib import contextmanager
from dataclasses import dataclass,asdict
from functools import lru_cache
from pathlib import Path
from typing import Callable

BASE=Path(__file__).resolve().parent
sys.path.insert(0,str(BASE/'engine'))
from oligominer_thermo_qc import PROFILES,calculate_tm
from Bio.SeqUtils import MeltingTemp as mt
from sequence_tools import rc,sha,normalize,canonical_ref,make_bins,genomic_interval,region_metrics,validate_imagers
from qc_core import exact_pair_segments,_dot_pairs
from rna_avoidance import prepare_transcription, permitted_strands, candidate_rna_safety
from design_presets import preset_provenance
from strand_guard import require_fixed_strand, check_reference_polarity, audit_selected_strands
from compact_outputs import write_input_manifest, finalize_outputs, archive_transient, output_detail
from specificity_store import HitStore, DiskHits, ScoreStore, ScoreHits

VERSION='1.6.1'
BLAST_FIELDS='qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qseq sseq'.split()
# Resource partition only: every shortlisted candidate still contributes both
# queries, and every HSP has the same retention and intended-site rules.
# Reuse each BLAST process for a small checkpoint group. Native input batches
# are bounded separately below so repeat-heavy candidates are not pooled.
BLAST_BATCH_CANDIDATES=4

@dataclass
class Policy:
    profile:str='oligominer1_2nM'
    mode:str='full'
    arm_min:int=36
    arm_max:int=90
    length_mode:str='adaptive'
    preferred_arm_max:int=42
    extension_shortlist_per_bin:int=200
    selection_goal:str='maximize_sites'
    gc_min:float=0.0
    gc_max:float=76.0
    tm_min:float=45.0
    tm_max:float=55.0
    tm_strict:bool=True
    candidate_step:int=1
    variants_per_start:int=2
    shortlist_per_bin:int=400
    max_total_candidates:int=5000
    target_per_bin:int=100
    minimum_per_bin:int=50
    minimum_gap:int=5
    footprint:str='extended'
    bound_strand:str='+'
    avoid_gene_rna:bool=False
    rna_exact_min:int=18
    preset_id:str='errfi1'
    homopolymer_max:int=5
    human_policy:str='intended_only'
    repeat_policy:str='within_bin'
    high_offtarget_tm:float=42.0
    moderate_offtarget_tm:float=35.0
    max_high_offtarget_sites:int=0
    max_moderate_offtarget_sites:int=10
    reject_unknown_long:bool=True
    screen_rna:bool=True
    screen_cross_primary:bool=False
    max_strong_rna_sites:int=0
    minimum_arm_unpaired:float=0.40
    minimum_cofold_gain:float=-12.0
    maximum_pair_stem:int=11
    folding_temperature:float=42.0
    sodium_M:float=0.390
    allow_unresolved_bins:bool=False
    threads:int=2

    @classmethod
    def from_dict(cls,d:dict):
        known=set(cls.__dataclass_fields__)
        extra=set(d)-known
        if extra:raise ValueError('Unknown policy fields: '+', '.join(sorted(extra)))
        d=dict(d)
        # Older saved configurations may request four threads. Enforce the
        # v1.6.1 resource limit without changing their scientific thresholds.
        if type(d.get('threads')) is int and d['threads']>2:d['threads']=2
        p=cls(**d);p.validate();return p

    def validate(self):
        for key in ('arm_min','arm_max','minimum_per_bin','minimum_gap','homopolymer_max','max_high_offtarget_sites',
                    'max_moderate_offtarget_sites','max_strong_rna_sites','maximum_pair_stem','rna_exact_min','preferred_arm_max','extension_shortlist_per_bin'):
            value=getattr(self,key)
            if isinstance(value,bool) or not isinstance(value,int) or value<0:
                raise ValueError(key+' must be a nonnegative integer.')
        for key in ('tm_strict','reject_unknown_long','screen_rna','screen_cross_primary','allow_unresolved_bins','avoid_gene_rna'):
            if not isinstance(getattr(self,key),bool):raise ValueError(key+' must be true or false.')
        if self.homopolymer_max<1 or self.maximum_pair_stem<1:raise ValueError('Repeat and stem lengths must be positive.')
        if self.profile not in PROFILES:raise ValueError('Unknown thermodynamic profile.')
        if self.mode not in ('full','local_preview'):raise ValueError('Choose full selection or local-only preview.')
        if not 18<=self.arm_min<=self.arm_max<=120:raise ValueError('Arm lengths must satisfy 18 <= minimum <= maximum <=120.')
        if self.length_mode not in ('adaptive','fixed'):raise ValueError('Choose adaptive or fixed length search.')
        if self.length_mode=='adaptive' and not self.arm_min<=self.preferred_arm_max<=self.arm_max:raise ValueError('Adaptive lengths require minimum <= preferred maximum <= extension ceiling.')
        if not 1<=self.extension_shortlist_per_bin<=2000:raise ValueError('Extension shortlist must be 1-2000 candidates per region.')
        if self.selection_goal not in ('maximize_sites','target_count'):raise ValueError('Choose maximum compatible sites or a target-count cap.')
        if not 0<=self.gc_min<=self.gc_max<=100:raise ValueError('Invalid GC window.')
        if not self.tm_min<self.tm_max:raise ValueError('Minimum Tm must be below maximum Tm.')
        require_fixed_strand(self)
        if not 12<=self.rna_exact_min<=60:raise ValueError('RNA exact-complement threshold must be 12-60 bases.')
        if self.avoid_gene_rna and not self.screen_rna:raise ValueError('Gene-RNA avoidance requires complementary RNA screening. Disable the avoidance mode explicitly to omit it.')
        if self.preset_id not in ('errfi1','adaptive_45_55','ellenberg_informed','expanded_g5','custom'):raise ValueError('Unknown parameter preset.')
        if self.footprint not in ('extended','arm'):raise ValueError('Unknown footprint model.')
        if self.human_policy not in ('intended_only','report_only'):raise ValueError('Unknown human-hit policy.')
        if self.repeat_policy not in ('within_bin','unique'):raise ValueError('Unknown repeat policy.')
        for k in ('candidate_step','variants_per_start','shortlist_per_bin','max_total_candidates','target_per_bin','threads'):
            v=getattr(self,k)
            if not isinstance(v,int) or not 1<=v<=100000:raise ValueError(k+' must be a positive integer.')
        if self.candidate_step>50 or self.variants_per_start>10:raise ValueError('Excessive candidate stride/length variants.')
        if self.threads>2 or isinstance(self.threads,bool):raise ValueError('BLAST threads must be 1 or 2; v1.6.1 limits concurrent memory use.')
        if self.shortlist_per_bin>2000 or self.max_total_candidates>20000:raise ValueError('This version supports <=2000 candidates/bin and <=20000/run.')
        if not 0<=self.minimum_per_bin<=self.target_per_bin:raise ValueError('Minimum count must not exceed target count.')
        if not 0<=self.minimum_arm_unpaired<=1:raise ValueError('Unpaired probability must be 0-1.')
        if self.minimum_cofold_gain>0:raise ValueError('Cofold gain threshold must be <=0.')
        if not 0<=self.minimum_gap<=100:raise ValueError('Minimum oligo gap must be 0-100 bases.')
        for k in ('tm_min','tm_max','gc_min','gc_max','high_offtarget_tm','moderate_offtarget_tm','folding_temperature','sodium_M','minimum_cofold_gain'):
            if not math.isfinite(getattr(self,k)):raise ValueError(k+' must be finite.')
        if self.moderate_offtarget_tm>self.high_offtarget_tm:raise ValueError('Moderate off-target Tm must be <= high off-target Tm.')
        if not 0.01<=self.sodium_M<=1.5 or not 0<=self.folding_temperature<=95:raise ValueError('Invalid aqueous folding conditions.')

class Cancelled(RuntimeError):pass

def log_default(message:str,**kw):print(message,flush=True)

def file_sha(path:Path):
    import hashlib
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()

def occurrences(pattern:str,sequence:str):
    pos=sequence.find(pattern)
    while pos>=0:
        yield pos
        pos=sequence.find(pattern,pos+1)

def footprint(refseq:str,arm_start:int,arm_end:int,prefix:str,suffix:str,strand:str):
    if strand=='+':prefix,suffix=rc(suffix),rc(prefix)
    left=right=0
    while left<len(prefix) and arm_start-left-2>=0 and prefix[-1-left]==refseq[arm_start-left-2]:left+=1
    while right<len(suffix) and arm_end+right<len(refseq) and suffix[right]==refseq[arm_end+right]:right+=1
    return arm_start-left,arm_end+right

def exact_sites(arm:str,full:str,prefix:str,suffix:str,refs:list[dict]):
    sites=[]
    for r in refs:
        for strand,q in (('-',arm),('+',rc(arm))):
            for pos in occurrences(q,r['sequence']):
                a,b=pos+1,pos+len(arm);ea,eb=footprint(r['sequence'],a,b,prefix,suffix,strand)
                sites.append(dict(reference=r['id'],start=a,end=b,strand=strand,extended_start=ea,extended_end=eb))
    return sites

def same_bin_site(h:dict,b:dict,p:Policy,expected_strand=None):
    return h['reference']==b['reference'] and h['strand']==(expected_strand or p.bound_strand) and b['start']<=h['start']<=h['end']<=b['end']

def footprints_conflict(sites1:list[dict],sites2:list[dict],p:Policy):
    ak,bk=('extended_start','extended_end') if p.footprint=='extended' else ('start','end')
    for a in sites1:
        for b in sites2:
            if a['reference']==b['reference'] and not(a[bk]+p.minimum_gap<b[ak] or b[bk]+p.minimum_gap<a[ak]):return True
    return False

def _spatial_shortlist(pool:list[dict], b:dict, limit:int, shortest_first:bool, center:float):
    """Balance fixed 1-kb sectors, then spread starts within each sector.

    Midpoint/farthest-point sampling prevents adjacent, highly ranked length
    alternatives from spending a sector's entire expensive-QC budget. Sequence
    quality breaks spatial ties. This changes search order, never a QC gate.
    """
    nb=max(1,math.ceil(b['length']/1000));buckets=defaultdict(list)
    for c in pool:
        idx=max(0,min(nb-1,(c['start']-b['start'])//1000))
        buckets[idx].append(c)
    retained=[];anchors={}
    for k in buckets:
        # Virtual boundary anchors place the first choice near the center,
        # followed by choices in the largest remaining coordinate gaps.
        anchors[k]=[b['start']+1000*k-1,min(b['end'],b['start']+1000*(k+1)-1)+1]
    while len(retained)<limit:
        active=sorted(k for k in buckets if buckets[k])
        if not active:break
        take=min(limit-len(retained),len(active))
        # A partial round covers the full interval instead of favoring its left end.
        sectors=[active[min(len(active)-1,int((j+.5)*len(active)/take))] for j in range(take)]
        for k in sectors:
            def rank(c):
                distance=min(abs(c['start']-a) for a in anchors[k])
                quality=(c['arm_length'],abs(c['tm']-center)) if shortest_first else (abs(c['tm']-center),c['arm_length'])
                return (-distance,*quality,c['start'],c['id'])
            chosen=min(buckets[k],key=rank)
            retained.append(chosen);buckets[k].remove(chosen);anchors[k].append(chosen['start'])
    return retained


def generate_candidates(refs:list[dict],bins:list[dict],imagers:list[dict],policy:Policy,log:Callable=log_default,rna_model:dict|None=None):
    require_fixed_strand(policy)
    check_reference_polarity(refs,policy)
    if policy.avoid_gene_rna and rna_model is None:
        raise ValueError('Prepare the transcription model before RNA-aware candidate generation; missing directions are not a pass.')
    rna_model=rna_model or {}
    lookup={r['id']:r for r in refs};ilook={i['id']:i for i in imagers};all_candidates=[];counts=[]
    profile=PROFILES[policy.profile];center=(policy.tm_min+policy.tm_max)/2
    adaptive=policy.length_mode=='adaptive'
    for b in bins:
        sequence=lookup[b['reference']]['sequence'];imager=ilook[b['imager']];dock=imager['docking'];pre=dock+'TT';suf='TT'+dock
        stages=Counter();pool=[];seen=set()
        for pos in range(b['start']-1,b['end']-policy.arm_min+1,policy.candidate_step):
            variants=[]
            for size in range(policy.arm_min,min(policy.arm_max,b['end']-pos)+1):
                window=sequence[pos:pos+size];stages['windows_enumerated']+=1
                extended=adaptive and size>policy.preferred_arm_max
                stages['extension_windows_enumerated' if extended else 'preferred_windows_enumerated']+=1
                if set(window)-set('ACGT'):stages['ambiguous_base']+=1;continue
                gc=100*(window.count('G')+window.count('C'))/size
                if not policy.gc_min<=gc<=policy.gc_max:stages['GC_outside']+=1;continue
                if re.search(r'(A{%d,}|C{%d,}|G{%d,}|T{%d,})'%((policy.homopolymer_max+1,)*4),window) or re.search(r'([ACGT]{2})\1{7}',window):
                    stages['simple_repeat']+=1;continue
                strands=permitted_strands(rna_model.get(b['reference'],[]),pos+1,pos+size,policy.bound_strand) if policy.avoid_gene_rna else [policy.bound_strand]
                if not strands:stages['RNA_direction_conflict']+=1;continue
                for strand in strands:
                    arm=window if strand=='-' else rc(window)
                    value=calculate_tm(arm,profile)['adjusted_screen_C']
                    inside=policy.tm_min<value<policy.tm_max
                    if policy.tm_strict and not inside:stages['Tm_outside']+=1;continue
                    stages['composition_Tm_pass']+=1
                    variants.append((abs(value-center),size,arm,value,gc,pos+1,pos+size,strand))
            # Try the shortest passing arms in each band first. Longer alternatives
            # are independently screened; extension is not a waiver of any QC rule.
            bands=([v for v in variants if v[1]<=policy.preferred_arm_max],
                   [v for v in variants if v[1]>policy.preferred_arm_max]) if adaptive else (variants,)
            preferred_cheap_pass=bool(bands[0])
            for band_no,band in enumerate(bands):
                accepted=0
                ordered=sorted(band,key=lambda v:(v[1],v[0],v[2])) if adaptive else sorted(band)
                for _,size,arm,value,gc,a,e,strand in ordered:
                    if accepted>=policy.variants_per_start:break
                    if arm in seen:continue
                    seen.add(arm);full=pre+arm+suf
                    sites=exact_sites(arm,full,pre,suf,refs)
                    if not sites:raise AssertionError('Generated candidate has no intended exact site')
                    outside=[h for h in sites if not same_bin_site(h,b,policy,strand)]
                    if outside:stages['cross_bin_or_other_reference_exact']+=1;continue
                    rna_evidence=candidate_rna_safety(sites,rna_model) if policy.avoid_gene_rna else {'status':'GENE_SPANS_NOT_USED_FOR_EXCLUSION','sites':[], 'note':'Transcript-sequence screening is separate. Gene annotations were not used to reject this candidate by span.'}
                    if rna_evidence.get('conflicts'):stages['RNA_direction_conflict_at_repeat_copy']+=1;continue
                    if policy.repeat_policy=='unique' and len(sites)!=1:stages['repeat_excluded']+=1;continue
                    if any(footprints_conflict([x],[y],policy) for n,x in enumerate(sites) for y in sites[n+1:]):
                        stages['self_repeat_footprint_conflict']+=1;continue
                    if policy.footprint=='extended' and any(h['extended_start']<b['start'] or h['extended_end']>b['end'] for h in sites):
                        stages['flank_crosses_bin_boundary']+=1;continue
                    is_extended=adaptive and size>policy.preferred_arm_max
                    note=('No preferred-length window at this start passed the composition/Tm/direction screen.' if not preferred_cheap_pass else
                          'Longer alternative for specificity, folding or packing; not necessarily required by Tm.') if is_extended else 'Preferred-length candidate.'
                    pool.append(dict(id='C_'+sha(b['id']+'|'+imager['id']+'|'+arm)[:16],bin_id=b['id'],reference=b['reference'],
                      imager=imager['id'],imager_sequence=imager['sequence'],docking_sequence=dock,
                      imager_input_kind=imager.get('input_kind','imager'),imager_entered_sequence_5to3=imager.get('input_sequence_5to3',imager['sequence']),
                      dye=imager['dye'],arm=arm,full=full,prefix=pre,suffix=suf,
                      bound_strand=strand,rna_avoidance=rna_evidence,start=min(h['start'] for h in sites),end=max(h['end'] for h in sites),sites=sites,gc=gc,tm=value,tm_in_window=policy.tm_min<value<policy.tm_max,
                      arm_length=size,full_length=len(full),length_band='extended' if is_extended else ('preferred' if adaptive else 'fixed_range'),
                      extension_bases=max(0,size-policy.preferred_arm_max) if adaptive else 0,extension_note=note,
                      off_target_recheck_required=is_extended,offtarget_hits=[],fold=None,reject_reasons=[],status='candidate',
                      generation_start=a,generation_end=e))
                    accepted+=1
                stages['length_alternatives_not_retained']+=max(0,len(band)-accepted)
        stages['unique_locally_compatible_candidates']=len(pool)
        if adaptive:
            short=[c for c in pool if c['length_band']=='preferred'];long=[c for c in pool if c['length_band']=='extended']
            retained_short=_spatial_shortlist(short,b,policy.shortlist_per_bin,True,center)
            retained_long=_spatial_shortlist(long,b,policy.extension_shortlist_per_bin,True,center)
            retained=retained_short+retained_long
            stages['preferred_shortlisted_for_QC']=len(retained_short);stages['extension_shortlisted_for_QC']=len(retained_long)
            stages['preferred_locally_compatible']=len(short);stages['extension_locally_compatible']=len(long)
        else:
            retained=_spatial_shortlist(pool,b,policy.shortlist_per_bin,False,center)
            stages['preferred_shortlisted_for_QC']=len(retained);stages['extension_shortlisted_for_QC']=0
        stages['shortlisted_for_QC']=len(retained);stages['not_shortlisted']=len(pool)-len(retained)
        counts.append(dict(bin_id=b['id'],**dict(stages)));all_candidates+=retained
        log(f'{b["id"]}: {len(retained)} detailed-QC candidates ({stages["extension_shortlisted_for_QC"]} extended); all lengths use the same QC rules.',stage='Candidates')
    if len(all_candidates)>policy.max_total_candidates:
        raise ValueError(f'{len(all_candidates)} shortlisted candidates exceed the configured {policy.max_total_candidates} limit. Reduce normal/extension shortlist or region count, or explicitly raise the limit.')
    return all_candidates,counts

class LocalIndex:
    def __init__(self,refs:list[dict]):
        self.refs=refs;self.indices={}
    def index(self,r:dict,k:int):
        key=(r['id'],k)
        if key not in self.indices:
            idx=defaultdict(list);s=r['sequence']
            for i in range(len(s)-k+1):
                word=s[i:i+k]
                if set(word)<=set('ACGT'):idx[word].append(i)
            self.indices[key]=idx
        return self.indices[key]
    def full_matches(self,query:str,r:dict,mismatches=3):
        n=len(query);s=r['sequence'];blocks=mismatches+1;seeds=[]
        for k in range(blocks):
            a=k*n//blocks;b=(k+1)*n//blocks;seeds.append((a,query[a:b]))
        starts=set()
        for offset,seed in seeds:
            for j in self.index(r,len(seed)).get(seed,[]):
                st=j-offset
                if 0<=st<=len(s)-n:starts.add(st)
        hits=[]
        for st in sorted(starts):
            t=s[st:st+n];mm=sum(a!=b for a,b in zip(query,t))
            if mm<=mismatches:hits.append((st,mm))
        return hits
    def partials(self,query:str,r:dict,k=18):
        s=r['sequence'];idx=self.index(r,k);found=set()
        for i in range(len(query)-k+1):
            for j in idx.get(query[i:i+k],[]):
                if i and j and query[i-1]==s[j-1]:continue
                n=k
                while i+n<len(query) and j+n<len(s) and query[i+n]==s[j+n]:n+=1
                found.add((i,j,n))
        return sorted(found)

def expected_hit(h:dict,c:dict,binlookup:dict,reflookup:dict,extra_allow:list[dict]):
    if h.get('rna'):return False
    b=binlookup[c['bin_id']]
    if h['source']=='targets':
        return h['reference']==b['reference'] and h['bound_strand']==c['sites'][0]['strand'] and b['start']<=h['start']<=h['end']<=b['end']
    g=genomic_interval(reflookup[b['reference']],b['start'],b['end'],c['sites'][0]['strand'])
    if h['source']=='human' and g and g['reference']==canonical_ref(h['reference']) and h['bound_strand']==g['strand'] and g['start']<=h['start']<=h['end']<=g['end']:
        return True
    for rule in extra_allow:
        if rule.get('bin_id','*') not in ('*',c['bin_id']):continue
        if rule.get('source','human')!=h['source']:continue
        if canonical_ref(rule['reference'])==canonical_ref(h['reference']) and rule['strand']==h['bound_strand'] and rule['start']<=h['start']<=h['end']<=rule['end']:
            return True
    return False

@lru_cache(maxsize=8192)
def _perfect_tract_metrics(q:str,s:str):
    """Immutable sequence-only evidence; never cache a mutable hit/candidate.

    The original tract algorithm (including its handling of unequal lengths)
    is retained. Repeated genomic loci often share the same aligned strings.
    """
    runs=[];start=None
    for i,(a,b) in enumerate(zip(q,s)):
        if a==b and a in 'ACGT':
            if start is None:start=i
        elif start is not None:runs.append(q[start:i]);start=None
    if start is not None:runs.append(q[start:])
    return max(map(len,runs),default=0),tuple(runs)

def perfect_tracts(q:str,s:str):
    # Public callers still receive their own list, as before.
    return list(_perfect_tract_metrics(q,s)[1])

@lru_cache(maxsize=8192)
def duplex_estimate(q:str,s:str,profile:str,rna:bool=False) -> tuple[float|None,str]:
    """BLAST qseq/sseq have the same alignment orientation; physical target is complement(sseq).
    Unsupported mismatches/gaps yield UNKNOWN, never a made-up low Tm.
    RNA alignments remain sequence-only: no transplant of DNA parameters.
    """
    if rna:return None,'RNA_SEQUENCE_ONLY'
    if '-' in q+s or len(q)!=len(s):return None,'GAPPED_DUPLEX_NOT_MODELED'
    if len(q)<8 or set(q+s)-set('ACGT'):return None,'TOO_SHORT_OR_AMBIGUOUS'
    p=PROFILES[profile]
    try:
        val=mt.Tm_NN(q,c_seq=s.translate(str.maketrans('ACGT','TGCA')),nn_table=getattr(mt,p.nn_table),
             dnac1=p.higher_strand_nM,dnac2=p.lower_strand_nM,Na=p.sodium_mM,K=0,Tris=0,Mg=0,dNTPs=0,saltcorr=p.saltcorr,strict=True)
        if p.round_before_formamide:val=round(val,2)
        val-=p.formamide_percent*p.formamide_factor
        if not math.isfinite(val):raise ValueError('Nonfinite Tm')
        return float(val),'EXACT_SEGMENT_DNA_NN' if q==s else 'SUPPORTED_MISMATCH_DNA_NN'
    except (ValueError,KeyError,IndexError):return None,'UNSUPPORTED_MISMATCH_CONTEXT'

@lru_cache(maxsize=8192)
def _exact_tract_tm(q:str,s:str,profile:str):
    """Only sequence/profile-dependent DNA evidence; limits stay per hit."""
    runs=_perfect_tract_metrics(q,s)[1]
    exact_tms=[duplex_estimate(run,run,profile,False)[0] for run in runs if len(run)>=8]
    return max((v for v in exact_tms if v is not None),default=None)

def clear_alignment_caches():
    """Release bounded sequence evidence when a design finishes or stops."""
    _perfect_tract_metrics.cache_clear()
    _exact_tract_tm.cache_clear()
    duplex_estimate.cache_clear()

def enrich_hit(h:dict,c:dict,p:Policy):
    q=h['qseq'].upper();s=h['sseq'].upper()
    h['exact_tract']=_perfect_tract_metrics(q,s)[0]
    h['offtarget_tm'],h['tm_model']=duplex_estimate(q,s,p.profile,h.get('rna',False))
    # Complementary exact stems can still flag risk when the whole alignment is unsupported.
    # Preserve assignment/error order and keep all candidate-specific decisions
    # below outside these bounded caches.
    h['exact_tract_tm']=_exact_tract_tm(q,s,p.profile) if not h.get('rna') else None
    values=[h[k] for k in ('offtarget_tm','exact_tract_tm') if h[k] is not None]
    h['risk_tm']=max(values,default=None)
    h['on_minus_off_tm']=c['tm']-h['offtarget_tm'] if h['offtarget_tm'] is not None else None
    h['review']=h['exact_tract']>=18 or (h['identity']>=80 and (h['query_coverage']>=70 or h['alignment_length']>=30))
    h['high_tm']=h['risk_tm'] is not None and h['risk_tm']>=p.high_offtarget_tm
    h['moderate_tm']=h['risk_tm'] is not None and h['risk_tm']>=p.moderate_offtarget_tm
    h['unknown_long']=h['offtarget_tm'] is None and not h.get('rna') and h['identity']>=90 and (h['alignment_length']>=30 or h['query_coverage']>=90)
    h['strong_rna']=bool(h.get('rna') and h['identity']>=90 and (h['alignment_length']>=30 or h['query_coverage']>=90))
    h['rna_exact_risk']=bool(h.get('rna') and h['exact_tract']>=p.rna_exact_min)
    if h['rna_exact_risk']:h['review']=True
    return h

def local_qc(candidates,refs,bins,background,transcripts,allow,p,log,resource_guard=None):
    lookup={b['id']:b for b in bins};rl={r['id']:r for r in refs}
    searches=[(r,'targets',False) for r in refs]+[(r,'construct',False) for r in background]+[(r,'rna_local',True) for r in transcripts]
    idx=LocalIndex([r for r,_,_ in searches])
    for n,c in enumerate(candidates,1):
        if resource_guard:resource_guard.checkpoint('Local specificity')
        for r,source,rna in searches:
            full_hits=[]
            for strand,q in ([('+',rc(c['arm']))] if rna else [('-',c['arm']),('+',rc(c['arm']))]):
                for st,mm in idx.full_matches(q,r):
                    ss=r['sequence'][st:st+len(q)];ss=ss if strand=='-' else rc(ss)
                    h=dict(reference=r['id'],source=source,start=st+1,end=st+len(q),bound_strand='RNA' if rna else strand,
                        qseq=c['arm'],sseq=ss,alignment_length=len(q),query_start=1,query_end=len(q),query_coverage=100.,
                        identity=100*(len(q)-mm)/len(q),mismatches=mm,gap_opens=0,full_primary=False,rna=rna)
                    h['allowed']=expected_hit(h,c,lookup,rl,allow);full_hits.append(h)
                    if not h['allowed']:c['offtarget_hits'].append(enrich_hit(h,c,p))
            for strand,q in ([('+',rc(c['full']))] if rna else [('-',c['full']),('+',rc(c['full']))]):
                for qi,st,size in idx.partials(q,r,k=min(18,p.rna_exact_min) if rna else 18):
                    if any(h['bound_strand']==('RNA' if rna else strand) and h['start']<=st+1<=st+size<=h['end'] for h in full_hits):continue
                    qa=qi+1 if strand=='-' else len(q)-qi-size+1
                    seq=c['full'][qa-1:qa-1+size]
                    h=dict(reference=r['id'],source=source,start=st+1,end=st+size,bound_strand='RNA' if rna else strand,
                        qseq=seq,sseq=seq,alignment_length=size,query_start=qa,query_end=qa+size-1,query_coverage=100*size/len(q),
                        identity=100.,mismatches=0,gap_opens=0,full_primary=True,rna=rna)
                    h['allowed']=expected_hit(h,c,lookup,rl,allow)
                    if not h['allowed']:c['offtarget_hits'].append(enrich_hit(h,c,p))
        if isinstance(c['offtarget_hits'],ScoreHits):
            c['offtarget_hits'].store.stage_scores('local',[c['id']],p)
        if n==1 or n%50==0 or n==len(candidates):log(f'Local full-arm and junction checks: {n}/{len(candidates)}',stage='Local specificity')


def run_process(command:list[str],outlog:Path,cancel=None,resource_guard=None,*,env_overrides=None):
    from process_diagnostics import run_process as execute
    return execute(command,outlog,cancel,resource_guard,cancellation_class=Cancelled,env_overrides=env_overrides)


def _blast_internal_batch_letters(candidates):
    """Bound BLAST's native query batch to one candidate pair or smaller.

    BLAST reads whole queries until the base count reaches BATCH_SIZE. Queries
    are written arm then full primary. A threshold greater than every arm and
    no greater than any pair holds exactly one pair per native batch. Unusual
    length combinations fall back to one query per native batch. This affects
    resource partitioning only; no search flag or alignment limit is changed.
    """
    if not candidates:return 1
    letters=max(len(c['arm']) for c in candidates)+1
    return letters if letters<=min(len(c['arm'])+len(c['full']) for c in candidates) else 1


@contextmanager
def _blast_alignment_stream(command,outlog,path,cancel,resource_guard,stream_output,*,internal_batch_letters=None):
    """Standard consumes bounded pipe output; audit retains its original file.

    Successful iterator exhaustion also requires successful child completion.
    The caller cannot commit partial pipe output as completed evidence.
    """
    def execute(command,outlog,cancel=None,resource_guard=None):
        if internal_batch_letters is None:
            return run_process(command,outlog,cancel,resource_guard)
        return run_process(command,outlog,cancel,resource_guard,
                           env_overrides={'BATCH_SIZE':str(internal_batch_letters)})
    if stream_output:
        from blast_stream import stream_command
        with stream_command(command,outlog,path,cancel,resource_guard,
                            run_process=execute,cancellation_class=Cancelled) as stream:
            yield stream
    else:
        execute(command,outlog,cancel,resource_guard)
        with path.open('rb') as stream:yield stream


def _blast_database_identity(database):
    """Record local index identity without loading or hashing genome volumes.

    Outputs are never reused, even when this identity matches a previous run.
    The file metadata is provenance, not a claim of full database validation.
    """
    path=Path(database).expanduser().resolve()
    files=[]
    for item in sorted(path.parent.glob(path.name+'.*')):
        if item.is_file():
            info=item.stat()
            files.append(dict(path=str(item),bytes=info.st_size,mtime_ns=info.st_mtime_ns))
    return dict(path=str(path),index_files=files,
                metadata_sha256=sha(json.dumps(files,sort_keys=True)),reuse_allowed=False)


def _blast_write_receipt(path,record):
    temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(record,indent=2,sort_keys=True)+'\n')
    temporary.replace(path)

def blast_candidates(candidates,refs,bins,allow,databases,p,out,log,cancel=None,detail='standard',resource_guard=None,checkpoint_session=None):
    executable=shutil.which('blastn')
    if not executable:raise RuntimeError('blastn not found; launch with your existing .qc_env.')
    executable=str(Path(executable).resolve())
    qpath=out/'candidate_queries.fa';by={c['id']:c for c in candidates}
    with qpath.open('w') as f:
        for c in candidates:
            for tag,seq in (('A',c['arm']),('F',c['full'])):
                qid=c['id']+'_'+tag;f.write('>'+qid+'\n'+seq+'\n')
    runs=[];binlookup={b['id']:b for b in bins};reflookup={r['id']:r for r in refs}
    version=subprocess.check_output([executable,'-version'],text=True).strip()
    query_digest=file_sha(qpath)
    disk_store=next((c['offtarget_hits'].store for c in candidates if isinstance(c['offtarget_hits'],DiskHits)),None)
    attempts=out/'blast_attempts';attempts.mkdir(exist_ok=True)
    for source,db in databases.items():
        # A fresh directory is mandatory: an old interrupted TSV is never
        # evidence, even if its inputs appear identical or its size is nonzero.
        safe_source=re.sub(r'[^A-Za-z0-9_.-]','_',source)
        attempt=Path(tempfile.mkdtemp(prefix=safe_source+'-',dir=attempts))
        manifest_path=attempt/'search_manifest.json';rna=source.startswith('rna')
        database_identity=_blast_database_identity(db)
        batch_count=math.ceil(len(candidates)/BLAST_BATCH_CANDIDATES)
        manifest=dict(schema='raserfish.blast-batches.v1',status='INCOMPLETE',source=source,
            database=str(db),database_identity=database_identity,blast_version=version,
            candidate_count=len(candidates),query_count=2*len(candidates),query_sha256=query_digest,
            batch_candidate_limit=BLAST_BATCH_CANDIDATES,batch_count=batch_count,completed_batches=0,
            raw_alignments=0,retained_alignments=0,reuse_allowed=bool(checkpoint_session),batches=[])
        _blast_write_receipt(manifest_path,manifest)
        retained=total=0;subjects=defaultdict(set);all_output_digest=hashlib.sha256()
        log(f'BLAST: {source}; {len(candidates)} candidates, {2*len(candidates)} arm/full-primary queries in {batch_count} batches. All search settings and binding sites are evaluated; Standard saves exact scores, Full audit retains alignments.',stage='Reference searches')
        for batch_index,start in enumerate(range(0,len(candidates),BLAST_BATCH_CANDIDATES),1):
            batch_candidates=candidates[start:start+BLAST_BATCH_CANDIDATES]
            internal_batch_letters=_blast_internal_batch_letters(batch_candidates)
            checkpoint_key=f'blast:{source}:{batch_index}'
            saved=disk_store.checkpoint(checkpoint_key) if checkpoint_session else None
            if saved is not None:
                if saved.get('status')!='COMPLETE' or saved.get('candidate_ids')!=[c['id'] for c in batch_candidates]:
                    raise ValueError('Saved BLAST checkpoint does not match the current candidate batch.')
                receipt=dict(saved,reused_from_checkpoint=True)
                manifest['batches'].append(receipt)
                total+=receipt['raw_alignments'];retained+=receipt['retained_alignments']
                manifest.update(completed_batches=batch_index,raw_alignments=total,retained_alignments=retained)
                _blast_write_receipt(manifest_path,manifest)
                log(f'BLAST: {source}; restored verified batch {batch_index}/{batch_count} from saved progress.',stage='Reference searches')
                if cancel and cancel.is_set():raise Cancelled('Cancelled while restoring saved BLAST progress; analysis incomplete.')
                if resource_guard:resource_guard.checkpoint('Restoring verified BLAST checkpoints')
                continue
            stem=attempt.name+f'_{batch_index:04d}'
            batch_query=attempt/(stem+'_queries.fa');path=attempt/(stem+'_candidates_blast.tsv')
            outlog=attempt/(stem+'_blast.log')
            batch_lookup={}
            with batch_query.open('w') as stream:
                for c in batch_candidates:
                    for tag,seq in (('A',c['arm']),('F',c['full'])):
                        qid=c['id']+'_'+tag;batch_lookup[qid]=(c,tag,seq)
                        stream.write('>'+qid+'\n'+seq+'\n')
            cmd=[executable,'-query',str(batch_query),'-db',database_identity['path'],'-task','blastn-short','-strand','both',
                 '-dust','no','-soft_masking','false','-word_size','7','-evalue','1000','-reward','1','-penalty','-3',
                 '-gapopen','5','-gapextend','2','-max_target_seqs','100000','-num_threads',str(p.threads),
                 '-outfmt','6 '+' '.join(BLAST_FIELDS),'-out',str(path)]
            receipt=dict(batch=batch_index,status='INCOMPLETE',candidate_ids=[c['id'] for c in batch_candidates],
                query_count=len(batch_lookup),query_sha256=file_sha(batch_query),command=cmd,
                query_file=str(batch_query.relative_to(out)),raw_file=str(path.relative_to(out)),
                log_file=str(outlog.relative_to(out)),database_metadata_sha256=database_identity['metadata_sha256'])
            receipt.update(environment_overrides={'BATCH_SIZE':str(internal_batch_letters)},
                           native_batch_strategy='ONE_CANDIDATE_PAIR_OR_SMALLER')
            manifest['batches'].append(receipt);_blast_write_receipt(manifest_path,manifest)
            try:
                if cancel and cancel.is_set():raise Cancelled('Cancelled before BLAST batch; analysis incomplete.')
                if resource_guard:resource_guard.checkpoint(f'BLAST {source} batch {batch_index}/{batch_count}',force=True)
                log(f'BLAST: {source}; batch {batch_index}/{batch_count} ({len(batch_candidates)} candidates).',stage='Reference searches')
                if resource_guard:resource_guard.checkpoint('Parsing '+source+' BLAST alignments',force=True)
                if isinstance(disk_store,ScoreStore):
                    for candidate in batch_candidates:disk_store.restore_local_source(source,candidate['id'])
                batch_total=batch_retained=batch_output_bytes=0;batch_digest=hashlib.sha256()
                # Hash and parse one line at a time; do not read a complete
                # BLAST output, or retain raw HSP rows, in Python memory.
                stream_output=detail=='standard' and isinstance(disk_store,ScoreStore)
                with _blast_alignment_stream(cmd,outlog,path,cancel,resource_guard,stream_output,
                                             internal_batch_letters=internal_batch_letters) as stream:
                    def lines():
                        nonlocal batch_output_bytes
                        for raw in stream:
                            batch_output_bytes+=len(raw)
                            batch_digest.update(raw);all_output_digest.update(raw)
                            yield raw.decode('utf-8')
                    for row in csv.DictReader(lines(),delimiter='\t',fieldnames=BLAST_FIELDS):
                        batch_total+=1;c,tag,seq=batch_lookup[row['qseqid']]
                        if batch_total==1 or batch_total%10000==0:
                            if cancel and cancel.is_set():raise Cancelled('Cancelled during BLAST parsing; analysis incomplete.')
                            # During streaming the owned subprocess runner samples
                            # this guard; avoid racing its non-thread-safe state.
                            # Forced checks surround the stream, and the separate
                            # worker supervisor remains active throughout.
                            if resource_guard and not stream_output:resource_guard.checkpoint('Parsing '+source+' BLAST alignments')
                        q1,q2=int(row['qstart']),int(row['qend']);s1,s2=int(row['sstart']),int(row['send'])
                        if q1>q2 or int(row['qlen'])!=len(seq) or row['qseq'].replace('-','').upper()!=seq[q1-1:q2]:
                            raise ValueError('BLAST query differs from the current candidate. Aborting.')
                        if disk_store:disk_store.record_subject(source,c['id'],row['sseqid'])
                        else:subjects[c['id']].add(row['sseqid'])
                        if rna and s1<=s2:continue
                        h=dict(reference=row['sseqid'],source=source,start=min(s1,s2),end=max(s1,s2),bound_strand='RNA' if rna else '-' if s1<=s2 else '+',
                              qseq=row['qseq'].upper(),sseq=row['sseq'].upper(),alignment_length=int(row['length']),query_start=q1,query_end=q2,
                              query_coverage=100*(q2-q1+1)/len(seq),identity=float(row['pident']),mismatches=int(row['mismatch']),
                              gap_opens=int(row['gapopen']),full_primary=tag=='F',rna=rna,evalue=float(row['evalue']),bitscore=float(row['bitscore']))
                        h['allowed']=expected_hit(h,c,binlookup,reflookup,allow)
                        # Exactly the previous scientific retention rules.
                        ex=_perfect_tract_metrics(h['qseq'],h['sseq'])[0]
                        potential=ex>=(min(18,p.rna_exact_min) if rna else 18) or (h['identity']>=80 and (h['query_coverage']>=70 or h['alignment_length']>=30))
                        if potential or h['allowed']:
                            c['offtarget_hits'].append(enrich_hit(h,c,p));batch_retained+=1
                if cancel and cancel.is_set():raise Cancelled('Cancelled during BLAST parsing; analysis incomplete.')
                if resource_guard:resource_guard.checkpoint('Finalizing '+source+' BLAST batch',force=True)
                total+=batch_total;retained+=batch_retained
                receipt.update(raw_alignments=batch_total,retained_alignments=batch_retained,
                    output_sha256=batch_digest.hexdigest(),output_bytes=batch_output_bytes,
                    output_transport='BOUNDED_PIPE' if stream_output else 'FILE')
                if isinstance(disk_store,ScoreStore):
                    disk_store.stage_scores(checkpoint_key,[c['id'] for c in batch_candidates],p,source)
                if checkpoint_session:
                    # The database is the authority. Its completion marker and
                    # All new scores (Standard) or raw evidence (Full audit)
                    # commit together before raw cleanup.
                    committed=dict(receipt,status='COMPLETE',raw_file=None,query_file=None,
                        raw_retention='COMPRESSED_FULL_AUDIT' if detail=='full_audit' else 'DISCARDED_STANDARD',
                        evidence_origin=str(out.resolve()),reused_from_checkpoint=False)
                    disk_store.commit_checkpoint(checkpoint_key,committed)
                    checkpoint_session.update_progress(completed_batches=disk_store.checkpoint_count('blast:'),
                        summary=f'Saved BLAST {source} batch {batch_index}/{batch_count}; analysis remains incomplete.')
                elif disk_store:disk_store.flush()
                # Standard discards raw alignments after their exact scores
                # commit. Full audit explicitly retains detailed evidence.
                receipt['raw_file']=archive_transient(path,out,detail)
                receipt['query_file']=archive_transient(batch_query,out,detail)
                receipt['raw_retention']='COMPRESSED_FULL_AUDIT' if detail=='full_audit' else 'DISCARDED_STANDARD'
                receipt['status']='COMPLETE'
                _blast_write_receipt(attempt/(stem+'_complete.json'),receipt)
                manifest.update(completed_batches=batch_index,raw_alignments=total,retained_alignments=retained)
                _blast_write_receipt(manifest_path,manifest)
            except BaseException as error:
                cleanup_errors=[]
                if checkpoint_session or isinstance(disk_store,ScoreStore):
                    try:disk_store.rollback()
                    except BaseException as cleanup:cleanup_errors.append('Rollback: '+str(cleanup))
                if detail=='standard' and isinstance(disk_store,ScoreStore):
                    # Partial raw output is never reusable checkpoint evidence.
                    # Preserve diagnostics, but do not strand large TSVs after
                    # failure. Cleanup errors must never mask its actual cause.
                    try:receipt['discarded_partial_bytes']=path.stat().st_size if path.exists() else 0
                    except OSError as cleanup:cleanup_errors.append('Partial output size: '+str(cleanup))
                    for transient in (path,batch_query):
                        try:
                            if transient.is_file() and not transient.is_symlink():transient.unlink()
                        except OSError as cleanup:cleanup_errors.append('Partial output cleanup: '+str(cleanup))
                    receipt.update(raw_file=None,query_file=None,raw_retention='DISCARDED_INCOMPLETE_STANDARD')
                receipt.update(status='INCOMPLETE',error=str(error),error_type=type(error).__name__)
                if cleanup_errors:receipt['cleanup_errors']=cleanup_errors
                manifest.update(status='INCOMPLETE',failed_batch=batch_index,error=str(error))
                try:_blast_write_receipt(manifest_path,manifest)
                except OSError:
                    # The worker's outer error propagation retains the primary
                    # exception even if a full disk prevents this extra receipt.
                    pass
                raise
        cap=disk_store.capped_candidates(source,100000) if disk_store else [c['id'] for c in candidates if len(subjects[c['id']])>=100000]
        for ident in cap:by[ident]['reject_reasons'].append('BLAST output may be capped; specificity unresolved')
        output_digest=all_output_digest.hexdigest()
        output_scope='CONCATENATED_BATCH_OUTPUT_BYTES'
        if checkpoint_session:
            # Raw files are already discarded in Standard mode. This ordered
            # digest is reproducible on resume and explicitly has a new scope.
            output_digest=sha(json.dumps([dict(batch=b['batch'],sha256=b['output_sha256'],bytes=b['output_bytes'])
                for b in manifest['batches']],sort_keys=True,separators=(',',':')))
            output_scope='ORDERED_BATCH_SHA256_AND_LENGTH_MANIFEST'
        manifest.update(status='COMPLETE',output_sha256=output_digest,
                        output_sha256_scope=output_scope,capped_candidate_ids=cap)
        _blast_write_receipt(manifest_path,manifest)
        runs.append(dict(source=source,database=str(db),database_identity=database_identity,
            command=manifest['batches'][0]['command'] if manifest['batches'] else [],command_scope='FIRST_BATCH; all exact commands are in batches',
            blast_version=version,raw_alignments=total,retained_alignments=retained,status='RUN_HEURISTIC_SEARCH',
            query_sha256=query_digest,output_sha256=output_digest,output_sha256_scope=output_scope,
            capped_candidate_ids=cap,raw_file=None,batches=manifest['batches'],batch_count=batch_count,
            search_manifest=str(manifest_path.relative_to(out)),reuse_allowed=bool(checkpoint_session),
            raw_retention='COMPRESSED_FULL_AUDIT' if detail=='full_audit' else 'DISCARDED_STANDARD',
            note='A heuristic search cannot prove exhaustive absence of binding. Alt/patch records are not automatically independent loci. All batches must complete; only identity-validated transactional checkpoints may be reused. Interrupted outputs are never evidence.'))
    return runs


def cluster_hits(hits:list[dict]):
    """Merge overlapping HSPs at the same reference/strand, including arm/full duplicates.
    Alternative contigs are NOT merged with a main chromosome without a mapping.
    Counts are site clusters in this reference representation, not physical copies.
    """
    groups=defaultdict(list)
    for h in hits:groups[(h['source'],canonical_ref(h['reference']),h['bound_strand'])].append(h)
    clusters=[]
    for (source,ref,strand),ar in groups.items():
        for h in sorted(ar,key=lambda x:(x['start'],x['end'])):
            if clusters and (clusters[-1]['source'],clusters[-1]['reference'],clusters[-1]['strand'])==(source,ref,strand) and h['start']<=clusters[-1]['end']:
                cl=clusters[-1];cl['end']=max(cl['end'],h['end']);cl['alignments'].append(h)
            else:
                clusters.append(dict(source=source,reference=ref,strand=strand,start=h['start'],end=h['end'],alignments=[h]))
    return clusters

def summarize_specificity(c:dict,p:Policy):
    if isinstance(c['offtarget_hits'],DiskHits):
        counts=c['offtarget_hits'].specificity_counts(p)
        high,moderate,unknown,rna,rna_risk=(counts[k] for k in ('high_tm_sites','moderate_tm_sites','unknown_long_sites','strong_rna_sites','rna_risk_sites'))
        cluster_count=counts['review_site_clusters'];allowed_count=counts['allowed_hsps'];ignored_count=counts['ignored_human_review_hsps']
    else:
        unexpected=[h for h in c['offtarget_hits'] if not h['allowed'] and h['review']]
        ignored=[h for h in unexpected if h['source']=='human' and p.human_policy=='report_only']
        active=[h for h in unexpected if not(h['source']=='human' and p.human_policy=='report_only')]
        clusters=cluster_hits(active)
        high=sum(any(h['high_tm'] for h in cl['alignments']) for cl in clusters)
        moderate=sum(any(h['moderate_tm'] for h in cl['alignments']) for cl in clusters)
        unknown=sum(any(h['unknown_long'] for h in cl['alignments']) for cl in clusters)
        rna=sum(any(h['strong_rna'] for h in cl['alignments']) for cl in clusters)
        rna_risk=sum(any(h['strong_rna'] or (p.avoid_gene_rna and h.get('rna_exact_risk',False)) for h in cl['alignments']) for cl in clusters)
        cluster_count=len(clusters);allowed_count=sum(h['allowed'] for h in c['offtarget_hits']);ignored_count=len(ignored)
    if high>p.max_high_offtarget_sites:c['reject_reasons'].append(f'{high} unintended high-Tm site cluster(s)')
    if moderate>p.max_moderate_offtarget_sites:c['reject_reasons'].append(f'{moderate} unintended moderate/high-Tm site clusters')
    if unknown and p.reject_unknown_long:c['reject_reasons'].append(f'{unknown} long strong alignment(s) with unresolved duplex Tm')
    if rna_risk>p.max_strong_rna_sites:c['reject_reasons'].append(f'{rna_risk} complementary RNA risk cluster(s): long/near-full match or exact tract meeting the RNA threshold')
    # Explicit heuristic weight; cannot be interpreted as expected background or probability.
    penalty=4*high+1.5*moderate+3*unknown+3*rna_risk+0.12*math.log1p(cluster_count)
    c['offtarget_alignment_count']=len(c['offtarget_hits'])
    c['specificity']=dict(review_site_clusters=cluster_count,high_tm_sites=high,moderate_tm_sites=moderate,unknown_long_sites=unknown,
       strong_rna_sites=rna,rna_risk_sites=rna_risk,allowed_hsps=allowed_count,ignored_human_review_hsps=ignored_count,
       risk_penalty=penalty,genome_scope='RUN' if p.mode=='full' else 'NOT_RUN',
       site_count_note='Overlapping HSP clusters in searched references, not proven independent physical copies.')

class Folding:
    def __init__(self,p:Policy,resource_guard=None):
        self.p=p;self.monomers=OrderedDict();self.pairs=OrderedDict();self.RNA=None;self.md=None
        self.resource_guard=resource_guard;self.evaluated_counts=Counter();self.monomer_evaluations=0
        if p.mode=='full':
            try:import RNA
            except ImportError:raise RuntimeError('ViennaRNA is required for full QC-guided selection. Use your installed .qc_env; no fallback was applied.')
            if not hasattr(RNA,'params_load_DNA_Mathews2004'):raise RuntimeError('DNA parameters unavailable; RNA defaults will not be used.')
            self.RNA=RNA;RNA.params_load_DNA_Mathews2004();md=RNA.md();md.temperature=p.folding_temperature
            md.salt=p.sodium_M;md.dangles=2;md.noGU=1;md.gquad=0;self.md=md
    def fold(self,c:dict):
        seq=c['full']
        if seq in self.monomers:return self.monomers[seq]
        if self.resource_guard:self.resource_guard.checkpoint('Individual-primary folding',force=True)
        if self.RNA is None:
            result=dict(status='NOT_RUN',engine='exact-stem geometry only; local preview',mean_arm_unpaired=None,mfe_kcal_mol=None)
        else:
            fc=self.RNA.fold_compound(seq,self.md);st,en=fc.mfe();fc.exp_params_rescale(en);_,pf=fc.pf();matrix=fc.bpp()
            unpaired=[max(0.,min(1.,1-sum(matrix[min(i,j)][max(i,j)] for j in range(1,len(seq)+1) if i!=j))) for i in range(1,len(seq)+1)]
            a=len(c['prefix']);b=a+len(c['arm'])
            result=dict(status='RUN_AQUEOUS_MODEL',engine='ViennaRNA '+self.RNA.__version__,structure=st,mfe_kcal_mol=float(en),ensemble_kcal_mol=float(pf),
               mean_arm_unpaired=sum(unpaired[a:b])/(b-a),base_unpaired=unpaired,
               mean_docking_unpaired=(sum(unpaired[:a-2])+sum(unpaired[b+2:]))/(len(seq)-len(c['arm'])-4))
            # Do not retain native fold compounds or the quadratic BPP matrix.
            del matrix,fc
        self.monomer_evaluations+=1
        self.monomers[seq]=result
        if len(self.monomers)>1024:self.monomers.popitem(last=False)
        return result
    def pair(self,a:dict,b:dict):
        if a['id']!=b['id'] and not self.p.screen_cross_primary:
            raise RuntimeError('Cross-primary interaction analysis is NOT RUN - USER DISABLED; a disabled calculation was requested.')
        key=tuple(sorted((a['id'],b['id'])))
        if key in self.pairs:return self.pairs[key]
        scope='SELF_DIMER' if a['id']==b['id'] else 'WITHIN_REGION' if a['bin_id']==b['bin_id'] else 'BETWEEN_REGION'
        if self.resource_guard:self.resource_guard.checkpoint('Self-dimer folding' if scope=='SELF_DIMER' else 'Cross-primary folding',force=True)
        segs=exact_pair_segments(a['full'],b['full'],6)
        stem=max((s['length'] for s in segs),default=0);gain=None;st=None;inter=[]
        if self.RNA is not None:
            fa=self.fold(a);fb=self.fold(b);fc=self.RNA.fold_compound(a['full']+'&'+b['full'],self.md);st,en=fc.mfe()
            gain=float(en-fa['mfe_kcal_mol']-fb['mfe_kcal_mol'])
            inter=[(i,j-len(a['full'])) for i,j in _dot_pairs(st) if i<len(a['full'])<=j]
            del fc
        blocked=stem>self.p.maximum_pair_stem or (gain is not None and inter and gain<self.p.minimum_cofold_gain)
        result=dict(a=a['id'],b=b['id'],a_name=a.get('name',a['id']),b_name=b.get('name',b['id']),
            max_exact_stem=stem,mfe_gain_vs_separate=gain,structure=st,intermolecular_pairs=inter,
            blocked=bool(blocked),scope=scope,status='RUN_AQUEOUS_MODEL' if gain is not None else 'GEOMETRY_ONLY',
            reason='pair exceeds the configured heuristic compatibility limit' if blocked else '')
        self.evaluated_counts[scope]+=1
        self.pairs[key]=result
        if len(self.pairs)>4096:self.pairs.popitem(last=False)
        return result


def _selection_pass(candidates,bins,p,fold,style,initial=None,cancel=None):
    """One feasible greedy solution; no gate or quota is relaxed in this helper."""
    if cancel and cancel.is_set():raise Cancelled('Cancelled by user; original inputs are unchanged.')
    selected=list(initial or []);selected_ids={c['id'] for c in selected}
    selected_by={b['id']:[c for c in selected if c['bin_id']==b['id']] for b in bins}
    available={b['id']:[c for c in candidates if c['bin_id']==b['id'] and not c['reject_reasons'] and c['id'] not in selected_ids] for b in bins}
    decisions={c['id']:('selected','') for c in selected};finished=set()
    while len(finished)<len(bins):
        eligible=[b for b in bins if b['id'] not in finished]
        b=min(eligible,key=lambda x:(len(selected_by[x['id']])/p.target_per_bin,x['reference'],x['number']))
        ident=b['id'];current=selected_by[ident]
        if (p.selection_goal=='target_count' and len(current)>=p.target_per_bin) or not available[ident]:
            finished.add(ident);continue
        if cancel and cancel.is_set():raise Cancelled('Cancelled by user; original inputs are unchanged.')
        metrics=region_metrics(b,current);hist=metrics['sector_counts']
        centers=[b['start']-1,*metrics['site_centers'],b['end']+1]
        def rank(c):
            size=c.get('arm_length',len(c.get('arm','')))
            if style=='packing':
                # Earlier finishing single sites pack efficiently; repeated intended
                # sites carry their actual nonoverlapping site count, not duplicate rows.
                return (-len(c['sites']),c['end'],size,abs(c['tm']-(p.tm_min+p.tm_max)/2),c['id'])
            sub={min(len(hist)-1,max(0,int(((h['start']+h['end'])/2-b['start'])//1000))) for h in c['sites']}
            diversity=sum(1/(1+hist[i]) for i in sub)
            gap_fill=sum(min(abs((h['start']+h['end'])/2-x) for x in centers) for h in c['sites'])
            unpaired=c['fold'].get('mean_arm_unpaired')
            quality=(unpaired if unpaired is not None else 0.5)-0.15*abs(c['tm']-(p.tm_min+p.tm_max)/2)
            return (-len(c['sites']),-diversity,-gap_fill,-quality+c['specificity']['risk_penalty'],size,c['start'],c['id'])
        chosen=None
        for c in sorted(available[ident],key=rank):
            if cancel and cancel.is_set():raise Cancelled('Cancelled by user; original inputs are unchanged.')
            conflict=next((q for q in selected if footprints_conflict(c['sites'],q['sites'],p)),None)
            if conflict:
                decisions[c['id']]=('not_selected_spacing','Conflicts with '+conflict['id']);available[ident].remove(c);continue
            if fold.pair(c,c)['blocked']:
                decisions[c['id']]=('not_selected_pair','Self-interaction exceeds the configured compatibility limit');available[ident].remove(c);continue
            pairblock=None
            for q in selected if p.screen_cross_primary else ():
                if cancel and cancel.is_set():raise Cancelled('Cancelled by user; original inputs are unchanged.')
                if fold.pair(c,q)['blocked']:pairblock=q;break
            if pairblock:
                decisions[c['id']]=('not_selected_pair','Incompatible with '+pairblock['id']);available[ident].remove(c);continue
            chosen=c;break
        if chosen:
            decisions[chosen['id']]=('selected','');available[ident].remove(chosen);selected.append(chosen);selected_by[ident].append(chosen)
        else:finished.add(ident)
    return selected,decisions


def select_candidates(candidates:list[dict],bins:list[dict],p:Policy,fold:Folding,log:Callable=log_default,cancel=None,selection_report:list|None=None):
    """Compare feasible heuristics, then attempt bounded swaps and infill.

    Exact-site count is the primary objective. At equal count, minimize the
    largest uncovered gap (including region edges), then regional gaps and
    1-kb-sector unevenness. This is not a certified global optimum.
    """
    if cancel and cancel.is_set():raise Cancelled('Cancelled by user; original inputs are unchanged.')
    bins_by={b['id']:b for b in bins};reports=[];solutions=[]
    def score(items):
        groups={b['id']:[c for c in items if c['bin_id']==b['id']] for b in bins}
        metrics=[region_metrics(b,groups[b['id']]) for b in bins]
        site_count=sum(len(c['sites']) for c in items)
        filled=sum(len(group)>=p.minimum_per_bin for group in groups.values())
        floor=min((len(group) for group in groups.values()),default=0)
        gaps=[m['largest_gap_bp'] for m in metrics]
        # Scale the expected occupancy by sector length, including short final sectors.
        uneven=sum(sum((n-m['exact_sites']*min(1000,b['length']-1000*i)/b['length'])**2
                       for i,n in enumerate(m['sector_counts'])) for b,m in zip(bins,metrics))
        bases=sum(c.get('arm_length',len(c.get('arm',''))) *len(c['sites']) for c in items)
        return site_count,-max(gaps,default=0),-sum(gaps),-uneven,filled,floor,len(items),-bases
    def run(label,items,style,initial=None):
        log('Trying '+label+'; retained individual QC and spacing rules remain active. Cross-primary screening: '+('ON' if p.screen_cross_primary else 'OFF'),stage='Joint selection')
        chosen,decisions=_selection_pass(items,bins,p,fold,style,initial,cancel)
        report=dict(strategy=label,selected_primaries=len(chosen),exact_sites=sum(len(c['sites']) for c in chosen),
                    extended_primaries=sum(c.get('length_band')=='extended' for c in chosen),
                    per_region={b['id']:sum(c['bin_id']==b['id'] for c in chosen) for b in bins},chosen=False,
                    largest_gap_bp=-score(chosen)[1],search_method='deterministic feasible heuristic; no global-optimum certificate',
                    repair_attempts=0,repair_accepted=0)
        reports.append(report);solutions.append((chosen,decisions,score(chosen)))
        return chosen
    if p.selection_goal=='maximize_sites':
        if p.length_mode=='adaptive':
            base=run('preferred_length_baseline',[c for c in candidates if c.get('length_band')!='extended'],'packing')
            run('preferred_baseline_plus_longer_infill',candidates,'packing',base)
        run('all_lengths_site_packing',candidates,'packing')
        run('all_lengths_spatial_balance',candidates,'distribution')
    else:
        run('balanced_target_count',candidates,'distribution')
    best=max(range(len(solutions)),key=lambda i:solutions[i][2]);selected,decisions,best_score=solutions[best];reports[best]['chosen']=True
    report=reports[best]
    report['pre_repair_exact_sites']=best_score[0];report['pre_repair_largest_gap_bp']=-best_score[1]
    # Bound both seed scans and full greedy refills. Each trial removes at most
    # two blockers, inserts one QC-passing alternate, and attempts full infill.
    # Only strict lexicographic improvements replace the best feasible solution.
    # The original baselines remain available; scientific gates are never waived.
    limits=dict(max_rounds=2,max_seed_checks=128,max_refills=32,max_removals=2)
    report['repair_limits']=limits;seed_checks=attempts=accepted=0
    if p.selection_goal=='maximize_sites':
        log('Repairing sparse sectors with bounded swaps and infill; exact-site count cannot decrease.',stage='Joint selection')
        for _ in range(limits['max_rounds']):
            selected_ids={c['id'] for c in selected}
            metrics={b['id']:region_metrics(b,[c for c in selected if c['bin_id']==b['id']]) for b in bins}
            def seed_rank(c):
                b=bins_by[c['bin_id']];m=metrics[b['id']]
                anchors=[b['start']-1,*m['site_centers'],b['end']+1]
                need=sum(min(abs((h['start']+h['end'])/2-a) for a in anchors) for h in c['sites'])
                return (-len(c['sites']),-need,c.get('arm_length',len(c.get('arm',''))),c['id'])
            alternatives=sorted((c for c in candidates if not c['reject_reasons'] and c['id'] not in selected_ids),key=seed_rank)
            improved=False
            for seed in alternatives:
                if cancel and cancel.is_set():raise Cancelled('Cancelled by user; original inputs are unchanged.')
                if attempts>=limits['max_refills'] or seed_checks>=limits['max_seed_checks']:break
                seed_checks+=1
                if seed['id'] in {c['id'] for c in selected} or fold.pair(seed,seed)['blocked']:continue
                blockers=[]
                for q in selected:
                    if cancel and cancel.is_set():raise Cancelled('Cancelled by user; original inputs are unchanged.')
                    if footprints_conflict(seed['sites'],q['sites'],p) or (p.screen_cross_primary and fold.pair(seed,q)['blocked']):
                        blockers.append(q['id'])
                        if len(blockers)>limits['max_removals']:break
                if len(blockers)>limits['max_removals']:continue
                initial=[c for c in selected if c['id'] not in blockers]+[seed]
                attempts+=1
                trial,trial_decisions=_selection_pass(candidates,bins,p,fold,'packing',initial,cancel)
                trial_score=score(trial)
                if trial_score>best_score:
                    selected,decisions,best_score=trial,trial_decisions,trial_score
                    accepted+=1;improved=True
            if not improved or attempts>=limits['max_refills'] or seed_checks>=limits['max_seed_checks']:break
    report.update(selected_primaries=len(selected),exact_sites=best_score[0],largest_gap_bp=-best_score[1],
                  extended_primaries=sum(c.get('length_band')=='extended' for c in selected),
                  per_region={b['id']:sum(c['bin_id']==b['id'] for c in selected) for b in bins},
                  repair_attempts=attempts,repair_accepted=accepted,repair_seed_checks=seed_checks)
    if selection_report is not None:selection_report.extend(reports)
    selected_ids={c['id'] for c in selected}
    for c in candidates:
        if c['reject_reasons']:c['status']='rejected_by_QC';c.pop('selection_note',None)
        elif c['id'] in selected_ids:c['status']='selected';c.pop('selection_note',None)
        else:
            status,note=decisions.get(c['id'],('not_selected_quota_or_rank','Not used in the winning feasible selection; no automatic QC failure.'))
            c['status']=status;c['selection_note']=note
    for i,c in enumerate(selected):
        if cancel and cancel.is_set():raise Cancelled('Cancelled by user; original inputs are unchanged.')
        if c['reject_reasons']:raise AssertionError('QC-rejected candidate selected')
        if not p.arm_min<=c.get('arm_length',len(c.get('arm','')))<=p.arm_max:raise AssertionError('Selected arm exceeds configured length limits')
        if p.tm_strict and not p.tm_min<c['tm']<p.tm_max:raise AssertionError('Selected candidate violates Tm window')
        if fold.pair(c,c)['blocked']:raise AssertionError('Self-interaction violation')
        for other in selected[i+1:]:
            if cancel and cancel.is_set():raise Cancelled('Cancelled by user; original inputs are unchanged.')
            if footprints_conflict(c['sites'],other['sites'],p):raise AssertionError('Final footprint-spacing violation')
            if p.screen_cross_primary and fold.pair(c,other)['blocked']:raise AssertionError('Final pair-compatibility violation')
    selected.sort(key=lambda c:(bins_by[c['bin_id']]['reference'],bins_by[c['bin_id']]['number'],c['start']))
    numbers=Counter()
    for c in selected:
        numbers[c['bin_id']]+=1;c['name']=f'{c["bin_id"]}_{c["imager"]}_{numbers[c["bin_id"]]:03d}'
    for b in bins:
        group=[c for c in selected if c['bin_id']==b['id']]
        b['metrics']=region_metrics(b,group);b['metrics']['shortfall_to_target']=max(0,p.target_per_bin-b['metrics']['probes'])
        b['metrics']['bottom_strand_sites']=sum(h['strand']=='-' for c in group for h in c['sites'])
        b['metrics']['top_strand_sites']=sum(h['strand']=='+' for c in group for h in c['sites'])
        b['metrics']['below_minimum']=b['metrics']['probes']<p.minimum_per_bin
        b['metrics']['extended_primaries']=sum(c.get('length_band')=='extended' for c in group)
        sizes=[c.get('arm_length',len(c.get('arm',''))) for c in group]
        b['metrics']['min_arm_length']=min(sizes,default=None);b['metrics']['max_arm_length']=max(sizes,default=None)
    log(f'Kept {len(selected)} primaries / {sum(len(c["sites"]) for c in selected)} compatible exact sites using {reports[best]["strategy"]}.',stage='Joint selection')
    return selected


def validate_allow(rows:list[dict],bins:list[dict]):
    names={b['id'] for b in bins};out=[]
    for r in rows:
        a=int(r['start']);b=int(r['end']);ident=r.get('bin_id','*')
        if ident not in names and ident!='*':raise ValueError('Allowed interval names an unknown region: '+ident)
        if a<1 or b<a or r.get('strand') not in ('+','-'):raise ValueError('Allowed intervals need positive 1-based bounds and a physical target strand.')
        if r.get('source','human') not in ('human','construct'):raise ValueError('Allowed sites may not silently exempt RNA targets.')
        if not r.get('reference'):raise ValueError('Allowed interval needs an explicit chromosome/reference.')
        out.append(dict(r,start=a,end=b,bin_id=ident))
    return out

def write_csv(path:Path,rows:list[dict],fields:list[str]):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');w.writeheader()
        for row in rows:
            clean={}
            for k in fields:
                v=row.get(k,'')
                if isinstance(v,(dict,list,tuple)):v=json.dumps(v,separators=(',',':'))
                if isinstance(v,str) and v.startswith(('=','+','@')):v="'"+v
                clean[k]=v
            w.writerow(clean)


def design(config:dict,out:Path,log:Callable=log_default,cancel=None,resource_guard=None,*,checkpoint_dir=None) -> dict:
    """Use bounded disk-backed evidence, optionally saving durable progress.

    A resumed design gets a fresh output directory and opens its original
    checkpoint store in place. Old outputs and incomplete attempts are kept.
    """
    if checkpoint_dir is not None:
        context=dict(path=Path(checkpoint_dir),store=None,session=None)
        try:return _design(config,out,log,cancel,resource_guard,None,context)
        finally:
            if context['store'] is not None:context['store'].close()
            clear_alignment_caches()
    with tempfile.TemporaryDirectory(prefix='raserfish_hsps_') as scratch:
        store=(ScoreStore(Path(scratch)/'specificity_scores.sqlite',deferred_local_sources=config.get('databases',{}))
               if output_detail(config)=='standard' else HitStore(Path(scratch)/'retained_hits.sqlite'))
        try:return _design(config,out,log,cancel,resource_guard,store)
        finally:
            store.close()
            clear_alignment_caches()


def _design(config:dict,out:Path,log:Callable,cancel,resource_guard,hit_store,checkpoint_context=None) -> dict:
    detail=output_detail(config)
    p=Policy.from_dict(config.get('policy',{}));refs=config['references'];imagers=validate_imagers(config['imagers'])
    if resource_guard:resource_guard.checkpoint('Validating design inputs',force=True)
    if config.get('policy',{}).get('threads',2)!=p.threads:
        log('BLAST thread request reduced to 2 by the v1.6.1 resource limit.',stage='Resource policy')
    if not refs:raise ValueError('No target references.')
    for r in refs:
        r['sequence']=normalize(r['sequence']);r['length']=len(r['sequence']);r['sha256']=sha(r['sequence'])
        if r.get('mapping') and not r['mapping'].get('verified'):raise ValueError('Unverified genomic mapping must be removed or verified before selection.')
    check_reference_polarity(refs,p)
    if len(set(r['id'] for r in refs))!=len(refs):raise ValueError('Duplicate reference IDs.')
    if sum(r['length'] for r in refs)>1_000_000:raise ValueError('This first version limits local target sequence to 1 Mb per run.')
    bins,tails=make_bins(refs,int(config.get('bin_width',20000)),config.get('tail_policy','strict'),config.get('assignments'),config.get('bounds'))
    ilook={i['id']:i for i in imagers}
    if any(b['imager'] not in ilook for b in bins):raise ValueError('Choose a valid imager for EVERY region before starting.')
    repeated=[i for i,n in Counter(b['imager'] for b in bins).items() if n>1]
    if repeated and not p.allow_unresolved_bins:
        raise ValueError('Imager reused across regions: '+', '.join(repeated)+'. The same docking barcode cannot resolve those regions, even in different rounds. Add imagers or explicitly permit merged-label groups.')
    allow=validate_allow(config.get('allowed_intervals',[]),bins)
    background=config.get('background',[]);transcripts=list(config.get('transcripts',[]))
    for rr in background+transcripts:
        rr['sequence']=normalize(rr['sequence'].upper().replace('U','T'));rr['length']=len(rr['sequence'])
    rna_model,derived_rna,transcription_summary=prepare_transcription(refs,config.get('transcription',{}),p.avoid_gene_rna)
    transcripts+=derived_rna
    all_ids=[r['id'] for r in refs+background+transcripts]
    if len(set(all_ids))!=len(all_ids):raise ValueError('Target, construct and transcript reference IDs must be unique.')
    databases=config.get('databases',{})
    if p.mode=='full':
        if not databases.get('human'):raise RuntimeError('Full mode requires a locally installed human DNA BLAST database, even when human hits are report-only.')
        if p.screen_rna and not databases.get('rna_human'):raise RuntimeError('RNA screening was requested but its database is missing.')
        if not shutil.which('blastn'):raise RuntimeError('BLAST executable missing. Use your existing .qc_env.')
    if not p.screen_rna:databases={k:v for k,v in databases.items() if not k.startswith('rna')}
    if any(not str(r.get('id','')).strip() for r in refs+background+transcripts):raise ValueError('All reference IDs must be nonempty.')
    for r in refs+background+transcripts:
        if r.get('mapping') and (r['mapping'].get('assembly')!='GRCh38' or r['mapping']['end']-r['mapping']['start']+1!=len(r['sequence'])):
            raise ValueError('Genomic coordinate mapping must match the sequence length and assembly.')
    if out.exists() and any(out.iterdir()):raise ValueError('Refusing to overwrite a previous run.')
    out.mkdir(parents=True,exist_ok=True)
    effective_config=dict(config,policy=asdict(p),output_detail=detail,
                          bin_width=int(config.get('bin_width',20000)),tail_policy=config.get('tail_policy','strict'))
    input_hash=sha(json.dumps(effective_config,sort_keys=True))
    manifest=write_input_manifest(effective_config,bins,out,application_version=VERSION)
    before={r['id']:r['sha256'] for r in refs}
    log('Generating candidate probes before QC.',stage='Candidates')
    log('Probe-probe interaction screening: '+('ON' if p.screen_cross_primary else 'OFF — NOT RUN - USER DISABLED')+'. Individual-primary and self-dimer checks are retained.',stage='QC scope')
    fold=Folding(p,resource_guard)
    candidates,generation=generate_candidates(refs,bins,imagers,p,log,rna_model)
    if not candidates:raise ValueError('No candidates passed the configured sequence/region gates on the locked strand '+p.bound_strand+'. RNA conflicts were excluded, not rescued by switching strands. No constraints were relaxed.')
    checkpoint=None
    if checkpoint_context is not None:
        from design_checkpoints import CheckpointSession
        checkpoint=CheckpointSession(checkpoint_context['path'],effective_config,candidates,databases if p.mode=='full' else {},p)
        checkpoint_context['session']=checkpoint
        log('Validating saved progress and specificity scores.' if checkpoint.is_resume else 'Saving progress after each completed search batch and individual QC check.',stage='Saved progress')
        if cancel and cancel.is_set():raise Cancelled('Cancelled before checkpoint validation; analysis incomplete.')
        if resource_guard:resource_guard.checkpoint('Validating saved progress',force=True)
        hit_store=(ScoreStore(checkpoint.store_path,resume=checkpoint.is_resume,identity_digest=checkpoint.identity_digest,
                             deferred_local_sources=databases if p.mode=='full' else {})
                   if detail=='standard' else HitStore(checkpoint.store_path,resume=checkpoint.is_resume,identity_digest=checkpoint.identity_digest))
        checkpoint_context['store']=hit_store
        if checkpoint.is_resume and not hit_store.checkpoint('local') and (
                hit_store.checkpoint_count('blast:') or hit_store.checkpoint_count('fold:')):
            raise ValueError('Saved progress is missing its completed local-QC checkpoint; evidence is incomplete.')
        if cancel and cancel.is_set():raise Cancelled('Cancelled after checkpoint validation; analysis incomplete.')
        checkpoint.update_progress(total_candidates=len(candidates),
            total_batches=math.ceil(len(candidates)/BLAST_BATCH_CANDIDATES)*len(databases) if p.mode=='full' else 0,
            completed_batches=hit_store.checkpoint_count('blast:'),folded_candidates=hit_store.checkpoint_count('fold:'),
            summary='Verified saved progress; remaining analysis has not completed.' if checkpoint.is_resume else 'No completed batches yet; analysis remains incomplete.')
    for c in candidates:
        prior=c['offtarget_hits'];c['offtarget_hits']=hit_store.hits(c['id'])
        if not (checkpoint and hit_store.checkpoint('local')):
            for hit in prior:c['offtarget_hits'].append(hit)
    if checkpoint and hit_store.checkpoint('local'):
        log('Restored completed local full-arm and junction checks from verified saved progress.',stage='Local specificity')
    else:
        local_qc(candidates,refs,bins,background,transcripts,allow,p,log,resource_guard)
        if cancel and cancel.is_set():raise Cancelled('Cancelled during local specificity; analysis incomplete.')
        if checkpoint:
            hit_store.commit_checkpoint('local',dict(status='COMPLETE',candidate_ids=[c['id'] for c in candidates]))
            checkpoint.update_progress(summary='Local full-arm and junction checks saved; analysis remains incomplete.')
        else:hit_store.flush()
    blast_runs=[]
    if p.mode=='full':blast_runs=blast_candidates(candidates,refs,bins,allow,databases,p,out,log,cancel,detail,resource_guard,checkpoint)
    for n,c in enumerate(candidates,1):
        if resource_guard:resource_guard.checkpoint('Individual QC '+str(n)+'/'+str(len(candidates)),force=True)
        if cancel and cancel.is_set():raise Cancelled('Cancelled by user.')
        saved=hit_store.checkpoint('fold:'+c['id']) if checkpoint else None
        if saved is not None:
            if saved.get('status')!='COMPLETE' or saved.get('candidate_id')!=c['id']:
                raise ValueError('Saved individual-QC checkpoint does not match its candidate.')
            c.update(saved['candidate_qc'])
            fold.monomers[c['full']]=c['fold']
            fold.pairs[(c['id'],c['id'])]=c['self_pair']
            if len(fold.monomers)>1024:fold.monomers.popitem(last=False)
            if len(fold.pairs)>4096:fold.pairs.popitem(last=False)
            # Counts describe completed scientific checks across attempts, not
            # the number of engine invocations in this resumed worker alone.
            fold.monomer_evaluations=saved['monomer_evaluations']
            fold.evaluated_counts=Counter(saved['evaluated_counts'])
            if n==1 or n%50==0 or n==len(candidates):log(f'Restored monomer/self-pair checks: {n}/{len(candidates)}',stage='Folding')
            continue
        mapped=next(r for r in refs if r['id']==c['reference']).get('mapping')
        if p.mode=='full' and mapped:
            recovered_target=(c['offtarget_hits'].has_expected_genomic_target() if isinstance(c['offtarget_hits'],ScoreHits) else
                              any(h['source']=='human' and h['allowed'] and not h['full_primary'] and h['identity']==100 and h['query_coverage']==100 for h in c['offtarget_hits']))
            if not recovered_target:
                c['reject_reasons'].append('Expected perfect genomic target not recovered from the supplied database; reference identity unresolved')
        if cancel and cancel.is_set():raise Cancelled('Cancelled by user.')
        summarize_specificity(c,p)
        # Analyze every shortlisted candidate even if another test rejects it, for interpretable evidence.
        c['fold']=fold.fold(c)
        if c['fold'].get('mean_arm_unpaired') is not None and c['fold']['mean_arm_unpaired']<p.minimum_arm_unpaired:
            c['reject_reasons'].append('Mean arm unpaired probability below aqueous-model threshold')
        self_pair=fold.pair(c,c);c['self_pair']=self_pair
        if self_pair['blocked']:c['reject_reasons'].append('Self-pair exceeds compatibility threshold')
        if checkpoint:
            if cancel and cancel.is_set():raise Cancelled('Cancelled before saving individual QC; analysis incomplete.')
            hit_store.commit_checkpoint('fold:'+c['id'],dict(status='COMPLETE',candidate_id=c['id'],
                candidate_qc={key:c[key] for key in ('specificity','reject_reasons','fold','self_pair','offtarget_alignment_count')},
                monomer_evaluations=fold.monomer_evaluations,evaluated_counts=dict(fold.evaluated_counts)))
            checkpoint.update_progress(folded_candidates=n,summary=f'Saved individual folding and self-dimer checks: {n}/{len(candidates)}. Analysis remains incomplete.')
        if n==1 or n%50==0 or n==len(candidates):log(f'Monomer/self-pair checks: {n}/{len(candidates)}',stage='Folding')
    selection_search=[]
    selected=select_candidates(candidates,bins,p,fold,log,cancel,selection_search)
    if checkpoint:checkpoint.update_progress(summary='All requested search and individual-QC checkpoints saved. Selection/export must finish before this analysis is complete.')
    strand_audit=audit_selected_strands(selected,refs,bins,p)
    if p.avoid_gene_rna:
        for c in selected:
            if candidate_rna_safety(c['sites'],rna_model)['conflicts']:
                raise AssertionError('Selected primary conflicts with an annotated/declared transcription direction')
    for c in selected:
        if c['specificity']['rna_risk_sites']>p.max_strong_rna_sites:
            raise AssertionError('Selected primary exceeds the active RNA-complementarity limit')
    status='LOCAL_ONLY_DESIGN_PREVIEW' if p.mode=='local_preview' else 'QC_GUIDED_CANDIDATE_DESIGN' if p.screen_cross_primary else 'PROVISIONAL_DRAFT_CROSS_PRIMARY_DISABLED'
    warnings=[]
    if not p.screen_cross_primary:warnings.append('Probe-probe interaction screening: OFF. Cross-primary exact-complementarity/stem searches and thermodynamic cofolding within and between regions are NOT RUN - USER DISABLED. Individual-primary accessibility and self-dimer checks remain separate. Order files are DRAFT.')
    if not p.avoid_gene_rna:
        warnings.append('Whole-gene RNA-span exclusion is OFF. ' + ('Annotated-transcript sequence screening remains requested; short partial RNA hits are review evidence, not automatic rejections. ' if p.screen_rna else 'Public annotated-RNA database screening is also OFF. ') + 'This is not an RNA-free guarantee and does not assume RNase treatment was performed.')
    else:warnings.append('RNA avoidance applies to provided/annotated directions and searched RNA sequences, not to every possible unannotated or antisense RNA.')
    if any(x.get('unknown_strand_intervals',0) for x in transcription_summary):warnings.append('Some annotated transcript strands were unknown; overlapping candidates on both strands were excluded.')
    if p.length_mode=='adaptive':warnings.append(f'Adaptive lengths: prefer {p.arm_min}-{p.preferred_arm_max} nt; consider extensions through {p.arm_max} nt under unchanged QC. Finite shortlists and greedy searches do not guarantee a global maximum.')
    if p.selection_goal=='maximize_sites':warnings.append('Maximum compatible-site mode: the per-region target count is a comparison benchmark, not a stopping cap. Region widths remain fixed.')
    if tails:warnings.append('Trailing sequence was excluded to keep every tracing region the same length. See excluded tails.')
    warnings.append('All selected primaries must bind the locked '+p.bound_strand+' reference-relative DNA strand across ALL tracing regions. RNA criteria operate within the locked strand; there is never strand switching. Other runs must use the same trace-wide strand. This is not a guarantee of chromatid/homolog identity.')
    if repeated:warnings.append('Some regions share the same imager and are merged readout groups, not independently resolved tracing positions.')
    if p.human_policy=='report_only':warnings.append('Human genomic off-targets were reported but did NOT affect selection. This is a user override, not a specificity pass.')
    if not background:warnings.append('Complete construct/vector NOT PROVIDED; public human references do not cover a custom integration.')
    if p.mode=='local_preview':warnings.append('Genome/RNA BLAST and full folding NOT RUN. This preview must not be treated as a genome-screened order.')
    if any(b['metrics']['below_minimum'] for b in bins):warnings.append('One or more regions have fewer primaries than the requested minimum. No QC thresholds were silently relaxed.')
    if p.footprint=='arm':warnings.append('Spacing checks only designated arms; TT spacers/handles can extend a duplex.')
    # Explicitly collect retained self checks. A disabled cross-primary branch
    # never calls pair(), including during final export.
    pair_results=[fold.pair(c,c) for c in selected]
    if p.screen_cross_primary:
        pair_results.extend(fold.pair(c,other) for i,c in enumerate(selected) for other in selected[i+1:])
    selected_cross=sum(x['a']!=x['b'] for x in pair_results)
    cross_status=('RUN_AQUEOUS_MODEL' if p.mode=='full' else 'EXACT_STEMS_ONLY') if p.screen_cross_primary else 'NOT RUN - USER DISABLED'
    cross_summary=dict(enabled=p.screen_cross_primary,status=cross_status,within_region_status=cross_status,
        between_region_status=cross_status,scope=['within_region','between_region'],
        evaluated=fold.evaluated_counts['WITHIN_REGION']+fold.evaluated_counts['BETWEEN_REGION'],selected_evaluated=selected_cross)
    self_summary=dict(status='RUN_AQUEOUS_MODEL' if p.mode=='full' else 'GEOMETRY_ONLY',
        evaluated=fold.evaluated_counts['SELF_DIMER'],selected_evaluated=len(selected))
    checks=[dict(name='One fixed DNA-bound strand across ALL tracing regions',status=strand_audit['status']),
        dict(name='Whole-gene RNA-span exclusion',status='RUN_ANNOTATION_SCOPED' if p.avoid_gene_rna else 'NOT_USED_TRANSCRIPT_POLICY' if p.screen_rna else 'DISABLED_BY_USER'),
        dict(name='Derived unspliced local RNA-span search',status='RUN' if derived_rna else 'NO_DERIVED_SEQUENCES'),
        dict(name='Thermodynamics and composition used in selection',status='RUN'),
        dict(name='Local target/construct/RNA sequence specificity',status='RUN'),
        dict(name='Human DNA search',status='RUN_HEURISTIC' if p.mode=='full' else 'NOT_RUN'),
        dict(name='Complementary annotated RNA search',status='RUN_HEURISTIC' if p.mode=='full' and p.screen_rna else 'NOT_RUN'),
        dict(name='Individual-primary folding and arm accessibility',status='RUN_AQUEOUS_MODEL' if p.mode=='full' else 'NOT_RUN'),
        dict(name='Self-dimer compatibility (each primary with itself)',status=self_summary['status']),
        dict(name='Cross-primary interactions within regions',status=cross_status),
        dict(name='Cross-primary interactions between regions',status=cross_status),
        dict(name='Complete construct',status='SUPPLIED_LOCAL_SEARCH' if background else 'NOT_PROVIDED'),
        dict(name='OligoMiner LDA and whole-genome k-mer pipeline',status='NOT_IMPLEMENTED')]
    result=dict(version=VERSION,status=status,order_status='DRAFT',all_scientific_checks_complete=p.mode=='full' and p.screen_cross_primary,
       cross_primary_screening=cross_summary,self_dimer_screening=self_summary,
       title=config.get('title','RASER-FISH QC-guided design'),generated_utc=dt.datetime.now(dt.timezone.utc).isoformat(),
       policy=asdict(p),parameter_provenance=preset_provenance(asdict(p)),transcription_summary=transcription_summary,derived_rna_reference_count=len(derived_rna),references=refs,bins=bins,excluded_tails=tails,imagers=imagers,allowed_intervals=allow,
       rna_policy_summary=dict(mode='strict_gene_spans_and_sequences' if p.avoid_gene_rna else 'transcript_sequences_only' if p.screen_rna else 'public_rna_search_disabled',
           whole_gene_exclusion=p.avoid_gene_rna, annotated_rna_search_requested=p.screen_rna,
           exact_tract_is_hard_gate=p.avoid_gene_rna, strong_rule='identity >=90% AND (aligned length >=30 nt OR query coverage >=90%)',
           rnase_treatment='NOT_ASSUMED_OR_VERIFIED', note='Sequence-heuristic policy, not RNA-DNA Tm or evidence of expression. See execution checks for what actually ran.'),
       strand_audit=strand_audit,selection_search=selection_search,candidates=candidates,selected_ids=list(c['id'] for c in selected),generation=generation,selected_pair_interactions=pair_results,
       warnings=warnings,checks=checks,blast_runs=blast_runs,annotation_provenance=config.get('annotation_provenance',[]),
       counts=dict(candidates=len(candidates),selected=len(selected),regions=len(bins),exact_sites=sum(len(c['sites']) for c in selected),
         extended_primaries=sum(c.get('length_band')=='extended' for c in selected),
         QC_rejected=sum(bool(c['reject_reasons']) for c in candidates),pairs_evaluated_during_selection=sum(fold.evaluated_counts.values()),
         cross_primary_evaluations=cross_summary['evaluated'],self_dimer_evaluations=self_summary['evaluated'],monomer_evaluations=fold.monomer_evaluations,
         selected_pairs_rechecked=len(pair_results)),reference_hashes_before=before,
       reference_hashes_after={r['id']:sha(r['sequence']) for r in refs},input_config_sha256=input_hash,
       algorithm='QC-gated adaptive candidate bands; compare preferred-only baseline, longer infill, site-packing and spatial-balance greedy solutions. Maximize compatible exact sites among tried solutions, not a certified global optimum.',
       scope='DNA duplex Tm is model-based. Folding is aqueous without formamide/pool equilibrium. Counts are not experimental brightness. '
             'Full mode searches arms and complete primaries against the supplied genome/RNA databases; local preview does not. Gene-RNA avoidance is annotation- and declaration-scoped.  Only explicit intended intervals are exempt. No experimental validation.')
    if result['reference_hashes_after']!=before:raise AssertionError('Reference changed during design')
    write_csv(out/'selected_primaries.csv',[dict(Probe_Name=c['name'],Region=c['bin_id'],Reference=c['reference'],Readout=c['imager'],
       Pool=c['bin_id']+'_'+c['imager'],Imager_Sequence_5to3=c['imager_sequence'],Docking_Sequence_in_Primary_5to3=c['docking_sequence'],
       Imager_Input_Type=c['imager_input_kind'],Imager_Entered_Sequence_5to3=c['imager_entered_sequence_5to3'],DNA_Binding_Arm_5to3=c['arm'],Full_Primary_5to3=c['full'],Tm_C=c['tm'],GC_percent=c['gc'],
       Exact_sites=len(c['sites']),Arm_length_nt=c['arm_length'],Length_band=c['length_band'],Extension_bases=c['extension_bases'],Bound_DNA_strand=c['bound_strand'],RNA_direction_status=c['rna_avoidance']['status'],QC_scope=status) for c in selected],
       ['Probe_Name','Region','Reference','Readout','Pool','Imager_Sequence_5to3','Docking_Sequence_in_Primary_5to3','Imager_Input_Type','Imager_Entered_Sequence_5to3','DNA_Binding_Arm_5to3','Full_Primary_5to3','Tm_C','GC_percent','Exact_sites','Arm_length_nt','Length_band','Extension_bases','Bound_DNA_strand','RNA_direction_status','QC_scope'])
    # Explicit convention/provenance exports; neither file is an oligo order.
    inventory_rows=[dict(Imager_ID=i['id'],Imager_Sequence_5to3=i['sequence'],
        Docking_Sequence_in_Primary_5to3=i['docking'],Entered_Type=i['input_kind'],
        Entered_Sequence_5to3=i['input_sequence_5to3'],Dye_or_note=i['dye'],Source=i['source']) for i in imagers]
    fields=['Imager_ID','Imager_Sequence_5to3','Docking_Sequence_in_Primary_5to3',
        'Entered_Type','Entered_Sequence_5to3','Dye_or_note','Source']
    write_csv(out/'imager_inventory.csv', inventory_rows, fields)
    inventory_by_id={i['Imager_ID']:i for i in inventory_rows}
    write_csv(out/'readout_assignments.csv', [dict(Region=b['id'],Reference=b['reference'],
        **inventory_by_id[b['imager']]) for b in bins], ['Region','Reference']+fields)
    write_csv(out/'DRAFT_IDT_oPools.csv',[{'Pool name':c['bin_id']+'_'+c['imager'],'Sequence':c['full']} for c in selected],['Pool name','Sequence'])
    write_csv(out/'region_summary.csv',[dict(Region=b['id'],Reference=b['reference'],Start=b['start'],End=b['end'],Width_bp=b['length'],Imager=b['imager'],**b['metrics']) for b in bins],
       ['Region','Reference','Start','End','Width_bp','Imager','probes','exact_sites','coverage_percent','longest_gap','centroid','centroid_offset','distribution','shortfall_to_target','below_minimum','median_tm','bottom_strand_sites','top_strand_sites','extended_primaries','min_arm_length','max_arm_length','sector_width_bp','sector_counts','largest_gap_bp','median_spacing_bp','sector_count'])
    from compact_outputs import _csv as write_compact_csv
    write_compact_csv(out/'debug'/'candidate_decisions.csv.gz',(dict(Candidate=c['id'],Region=c['bin_id'],Imager=c['imager'],Imager_Sequence_5to3=c['imager_sequence'],Docking_Sequence_in_Primary_5to3=c['docking_sequence'],Imager_Input_Type=c['imager_input_kind'],Status=c['status'],Tm_C=c['tm'],GC_percent=c['gc'],
       Arm_length_nt=c['arm_length'],Length_band=c['length_band'],Extension_bases=c['extension_bases'],Extension_note=c['extension_note'],Reject_reasons=c['reject_reasons'],Selection_note=c.get('selection_note',''),Arm=c['arm'],Full_primary=c['full'],**c['specificity'],
       mean_arm_unpaired=c['fold'].get('mean_arm_unpaired')) for c in candidates),
       ['Candidate','Region','Imager','Imager_Sequence_5to3','Docking_Sequence_in_Primary_5to3','Imager_Input_Type','Status','Tm_C','GC_percent','Arm_length_nt','Length_band','Extension_bases','Extension_note','Reject_reasons','Selection_note','Arm','Full_primary',
        'review_site_clusters','high_tm_sites','moderate_tm_sites','unknown_long_sites','strong_rna_sites','rna_risk_sites','ignored_human_review_hsps','mean_arm_unpaired'],compressed=True)
    write_csv(out/'selected_binding_sites.csv',[dict(Probe=c['name'],Region=c['bin_id'],**h) for c in selected for h in c['sites']],
       ['Probe','Region','reference','start','end','strand','extended_start','extended_end'])
    write_csv(out/'selection_search.csv',selection_search,['strategy','selected_primaries','exact_sites','extended_primaries','per_region','chosen','repair_attempts','repair_accepted','pre_repair_exact_sites','pre_repair_largest_gap_bp','largest_gap_bp','search_method'])
    write_csv(out/'generation_counts.csv',generation,sorted(set(k for d in generation for k in d)))
    write_csv(out/'transcription_model.csv',[dict(Reference=x['reference'],Mode=x.get('mode','disabled'),Status=x['status'],**it) for x in transcription_summary for it in x.get('intervals',[])],
       ['Reference','Mode','Status','name','start','end','strand','source','scope'])
    (out/'strand_consistency.json').write_text(json.dumps(strand_audit,indent=2))
    write_csv(out/'strand_consistency.csv',strand_audit['per_region'],['Region','Reference','Required_bound_strand','Top_sites','Bottom_sites','Status'])
    (out/'parameter_provenance.json').write_text(json.dumps(result['parameter_provenance'],indent=2))
    if resource_guard:resource_guard.checkpoint('Final DRAFT export',force=True)
    result=finalize_outputs(result,manifest,out)
    (out/'READ_ME_RESULT.txt').write_text(status+'\n'+result['scope']+'\n\n'+'\n'.join(warnings)+'\n\nOutput detail: '+detail+'\n'
       'run_config.json is a reproducible input manifest. Keep its inputs/ directory alongside it.\n'
       'debug/candidate_decisions.csv.gz records every shortlisted candidate. design_result.json and the viewer contain selected primaries only.\n'
       'Self-dimer status: '+self_summary['status']+'. Cross-primary analysis: '+cross_status+'.\n'
       'Self-dimers and cross-primary interactions are reported separately. Omitted cross-primary analysis is never a PASS.\n'
       'Standard output saves exact specificity counts and scores only; per-batch intervals are discarded after summarizing. No per-alignment sequence payload is kept. Full audit adds compressed detailed candidates, alignments, selected pairs, BLAST queries and search files under debug/.\n'
       'External genome/RNA databases are not bundled. Reproduction also requires matching installed tools and databases.\n'
       '\nDRAFT_IDT_oPools.csv is a candidate order, not a validation certificate.\n')
    log(f'Done: {len(selected)} selected primaries across {len(bins)} equal regions. Historical orders unchanged.',stage='Complete')
    return result
