"""Public annotation adapters. Network errors are not evidence of absent biology.

Coordinates: UCSC/BED are 0-based half-open at input; UI is 1-based inclusive.
Only public coordinates/cell names are sent remotely, never private probe DNA.
"""
from __future__ import annotations
import gzip
import hashlib
import io
import json
import math
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable
from sequence_tools import normalize, sha, rc

BASE=Path(__file__).resolve().parent
CACHE=BASE/'cache'/'annotations'
_LAST=0.0

def fetch_json(url: str, refresh=False) -> dict | list:
    global _LAST
    CACHE.mkdir(parents=True,exist_ok=True)
    key=sha(url); dest=CACHE/(key+'.json')
    if dest.exists() and not refresh:
        return json.loads(dest.read_text())['response']
    if time.monotonic()-_LAST<1.05: time.sleep(1.05-(time.monotonic()-_LAST))
    req=urllib.request.Request(url,headers={'Accept':'application/json','User-Agent':'RASERFISH-Design-Studio/1.0 public-annotation-client'})
    try:
        try:
            response=urllib.request.urlopen(req,timeout=45)
        except urllib.error.HTTPError as exc:
            if exc.code==404 and '/search/' in url:
                body=json.loads(exc.read())
                if body.get('total')==0 or (body.get('@graph')==[] and 'no results' in body.get('notification','').lower()):
                    return body
            raise
        with response as fh:
            data=fh.read(40_000_001)
            if len(data)>40_000_000: raise RuntimeError('Public metadata response exceeds 40 MB.')
        result=json.loads(data)
    except Exception as e:
        raise RuntimeError(f'Annotation request failed ({type(e).__name__}): {url}. This is NOT a finding that no data exist.') from e
    finally:_LAST=time.monotonic()
    dest.write_text(json.dumps(dict(url=url,retrieved_utc=datetime.now(timezone.utc).isoformat(),response=result)))
    return result

def locus(text: str) -> tuple[str,int,int]:
    match=re.fullmatch(r'\s*(?:chr)?([0-9]+|X|Y|M):([\d,]+)-([\d,]+)\s*',text,re.I)
    if not match: raise ValueError('Use chr7:117475001-117485000 (1-based inclusive).')
    chrom='chr'+match[1].upper();a=int(match[2].replace(',',''));b=int(match[3].replace(',',''))
    if not 1<=a<=b or b-a+1>1_000_000: raise ValueError('Select a valid interval up to 1 Mb in this first version.')
    return chrom,a,b

def get_genomic_reference(text:str, getter:Callable=fetch_json) -> dict:
    chrom,a,b=locus(text)
    url='https://api.genome.ucsc.edu/getData/sequence?'+urllib.parse.urlencode(dict(genome='hg38',chrom=chrom,start=a-1,end=b)).replace('&',';')
    data=getter(url);seq=normalize(data.get('dna',''))
    if len(seq)!=b-a+1:raise RuntimeError('Sequence length does not match the requested interval.')
    return dict(id=f'{chrom}_{a}_{b}',description=f'GRCh38 {chrom}:{a}-{b}',sequence=seq,length=len(seq),sha256=sha(seq),
      features=[],mapping=dict(assembly='GRCh38',chromosome=chrom,start=a,end=b,orientation='+',verified=True),
      source=dict(provider='UCSC',url=url,loaded_utc=datetime.now(timezone.utc).isoformat()))

def fetch_genes(ref:dict, getter:Callable=fetch_json) -> dict:
    m=ref.get('mapping')
    if not m or not m.get('verified'):raise ValueError('Public gene tracks need a verified genomic coordinate mapping. A reporter filename is not a mapping.')
    chrom=m['chromosome'].removeprefix('chr');a,b=m['start'],m['end']
    url=f'https://rest.ensembl.org/overlap/region/human/{chrom}:{a}-{b}?feature=gene;feature=transcript;feature=exon;content-type=application/json'
    raw=getter(url)
    if not isinstance(raw,list):raise RuntimeError('Unexpected Ensembl feature response.')
    features=[]
    for f in raw:
        if f.get('assembly_name') not in (None,'GRCh38'):raise RuntimeError('Ensembl assembly does not match GRCh38; refusing to mix coordinates.')
        ga,gb=int(f['start']),int(f['end'])
        if m.get('orientation','+')=='+':la,lb=ga-a+1,gb-a+1;st='+' if f.get('strand')==1 else '-' if f.get('strand')==-1 else '?'
        else:la,lb=b-gb+1,b-ga+1;st='-' if f.get('strand')==1 else '+' if f.get('strand')==-1 else '?'
        features.append(dict(type=f.get('feature_type','gene'),label=f.get('external_name',f.get('id','feature')),
          start=la,end=lb,strand=st,source='Ensembl GRCh38',feature_id=f.get('id',''),
          parent=f.get('Parent',''),genomic_start=ga,genomic_end=gb,biotype=f.get('biotype','')))
    return dict(features=features,status='FETCHED',provider='Ensembl',url=url,assembly='GRCh38',
                loaded_utc=datetime.now(timezone.utc).isoformat(),cell_specific=False)

def cell_name(s:str) -> str:
    val=s.strip(); key=re.sub(r'[^A-Z0-9]','',val.upper())
    return {'IMR90':'IMR-90','HCT116':'HCT116'}.get(key,val)

def eligible_experiment(e:dict,cell:str) -> bool:
    if e.get('status')!='released':return False
    if e.get('biosample_ontology',{}).get('term_name')!=cell:return False
    if e.get('target',{}).get('label')!='H3K27ac':return False
    if e.get('assay_term_name')!='ChIP-seq':return False
    if e.get('treatments') or e.get('genetic_modifications'):return False
    if 'genetically modified' in e.get('biosample_summary','').lower() or ' treated with ' in e.get('biosample_summary','').lower():return False
    for r in e.get('replicates',[]):
        bios=r.get('library',{}).get('biosample',{}) if isinstance(r,dict) else {}
        if isinstance(bios,dict) and (bios.get('treatments') or bios.get('genetic_modifications')):return False
    return True

def discover_h3k27ac(requested:str, fallback=True,getter:Callable=fetch_json) -> dict:
    """Return source choices. A successful empty search allows fallback; an error does not."""
    req=cell_name(requested)
    cells=[req] if req else []
    if fallback:
        for c in ('HCT116','IMR-90'):
            if c not in cells:cells.append(c)
    if not cells:raise ValueError('Enter a cell line or allow HCT116/IMR-90 fallback.')
    attempts=[]
    for cell in cells:
        query=dict(type='Experiment',status='released',assay_term_name='ChIP-seq',**{'target.label':'H3K27ac',
          'biosample_ontology.term_name':cell},format='json',frame='embedded',limit='all')
        url='https://www.encodeproject.org/search/?'+urllib.parse.urlencode(query)
        raw=getter(url)
        exps=raw.get('@graph',[])
        choices=[]
        for e in exps:
            if not eligible_experiment(e,cell):continue
            files=e.get('files',[])
            for f in files:
                if isinstance(f,str):f=getter('https://www.encodeproject.org'+f+'?format=json')
                if f.get('status')!='released' or f.get('assembly') not in ('GRCh38','GRCh38.p14','hg38'):continue
                fmt=f.get('file_format','');sub=f.get('file_format_type','');ot=f.get('output_type','')
                peaks=fmt in ('bed','bigBed') and ('peak' in ot.lower() or 'Peak' in sub)
                signal=fmt=='bigWig' and ('signal' in ot.lower())
                if not (peaks or signal):continue
                audit={k:len(v) for k,v in (e.get('audit') or {}).items() if isinstance(v,list)}
                href=f.get('href') or ''
                if not href:continue
                choices.append(dict(accession=f['accession'],experiment=e['accession'],cell=cell,assembly=f['assembly'],
                   format=fmt,subtype=sub,output_type=ot,kind='peaks' if peaks else 'signal',
                   url=urllib.parse.urljoin('https://www.encodeproject.org',href),
                   md5=f.get('md5sum'),replicates=f.get('biological_replicates',[]),audit=audit,
                   experiment_url='https://www.encodeproject.org/experiments/'+e['accession']+'/',
                   file_url='https://www.encodeproject.org/files/'+f['accession']+'/',
                   laboratory=e.get('lab',{}).get('title',''),biosample_summary=e.get('biosample_summary','')))
        choices.sort(key=lambda c:(c['kind']!='peaks',c['audit'].get('ERROR',0),
                      c['audit'].get('NOT_COMPLIANT',0),not any(s in c['output_type'].lower() for s in ('optimal','replicated')),
                      c['format']!='bed',-len(c['replicates']),c['accession']))
        attempts.append(dict(cell=cell,experiments_returned=len(exps),compatible_files=len(choices),query=url))
        if choices:return dict(requested_cell=req or '(not specified)',used_cell=cell,
           is_proxy=not req or req!=cell,choices=choices,attempts=attempts,status='FOUND',
           note='Cell-line proxy; not a measurement in your cells.' if not req or req!=cell else 'Exact cell-name match; inspect experiment and biosample metadata.')
    return dict(requested_cell=req,used_cell=None,is_proxy=False,choices=[],attempts=attempts,status='NO_COMPATIBLE_RELEASED_FILE',
                note='No eligible untreated/unmodified GRCh38 H3K27ac file found. This does not mean H3K27ac is biologically absent.')

def get_track_bytes(choice:dict) -> tuple[bytes,str]:
    url=choice['url'];parsed=urllib.parse.urlparse(url)
    if parsed.scheme!='https' or parsed.hostname!='www.encodeproject.org':raise ValueError('Only ENCODE-discovered HTTPS download URLs are accepted.')
    CACHE.mkdir(parents=True,exist_ok=True)
    path=CACHE/(choice['accession']+'.bin')
    if not path.exists():
        req=urllib.request.Request(url,headers={'User-Agent':'RASERFISH-Design-Studio/1.0'})
        with urllib.request.urlopen(req,timeout=90) as f:
            content=f.read(120_000_001)
        if len(content)>120_000_000:raise RuntimeError('Track exceeds 120 MB download limit; choose another file or a bigWig track.')
        if choice.get('md5') and hashlib.md5(content).hexdigest()!=choice['md5']:raise RuntimeError('ENCODE file MD5 mismatch; download not accepted.')
        path.write_bytes(content)
    data=path.read_bytes()
    if choice.get('md5') and hashlib.md5(data).hexdigest()!=choice['md5']:raise RuntimeError('Cached track MD5 mismatch.')
    return data,sha(data)

def parse_bed_track(text:str,ref:dict,choice:dict,local=False) -> list[dict]:
    features=[];mapping=ref.get('mapping')
    for line in text.splitlines():
        if not line.strip() or line.startswith(('#','track','browser')):continue
        f=line.split('\t')
        if len(f)<3:f=line.split()
        if len(f)<3:raise ValueError('BED track requires chromosome/reference, start, and end.')
        if local:
            if f[0]!=ref['id']:continue
            a,b=int(f[1])+1,int(f[2])
        else:
            if not mapping:raise ValueError('No genomic mapping for the public track.')
            if f[0].removeprefix('chr')!=mapping['chromosome'].removeprefix('chr'):continue
            ga,gb=int(f[1])+1,int(f[2])
            if ga>mapping['end'] or gb<mapping['start']:continue
            if mapping.get('orientation','+')=='+':a,b=ga-mapping['start']+1,gb-mapping['start']+1
            else:a,b=mapping['end']-gb+1,mapping['end']-ga+1
        if a>ref['length'] or b<1:continue
        value=None
        if len(f)>6:
            try:value=float(f[6])
            except ValueError:pass
        if choice.get('kind')=='bedgraph' and len(f)>3:value=float(f[3])
        features.append(dict(type='H3K27ac' if choice.get('kind')!='annotation' else 'feature',
          label=f[3] if len(f)>3 and choice.get('kind')!='bedgraph' else choice.get('label','H3K27ac'),
          start=a,end=b,signal=value,source=choice.get('accession','user BED'),cell=choice.get('cell'),
          is_proxy=choice.get('is_proxy',False),score_type='peak signalValue' if choice.get('kind')!='bedgraph' else 'user bedGraph value'))
    return features

def load_h3k27ac(ref:dict,choice:dict,requested_cell:str) -> dict:
    if not ref.get('mapping',{}).get('verified'):raise ValueError('Public H3K27ac needs a verified genomic mapping.')
    if choice.get('assembly') not in ('GRCh38','GRCh38.p14','hg38'):raise ValueError('Track assembly differs; no automatic liftOver.')
    choice=dict(choice,is_proxy=cell_name(requested_cell)!=choice['cell'])
    if choice['format']=='bed':
        blob,checksum=get_track_bytes(choice)
        raw=gzip.decompress(blob) if blob[:2]==b'\x1f\x8b' else blob
        features=parse_bed_track(raw.decode('utf-8'),ref,choice)
        return dict(features=features,signal=[],status='FETCHED_PEAKS',metadata=dict(choice,sha256=checksum,
              loaded_utc=datetime.now(timezone.utc).isoformat()))
    try:import pyBigWig
    except ImportError:raise RuntimeError('This file requires optional pyBigWig. Choose a BED peak file, or install pyBigWig in the existing .qc_env.')
    if not pyBigWig.remote:raise RuntimeError('This pyBigWig build has no remote support. Choose a BED peak file.')
    m=ref['mapping'];bw=pyBigWig.open(choice['url'])
    if bw is None:raise RuntimeError('Unable to open ENCODE bigWig/bigBed.')
    try:
        chrom=m['chromosome'];a=m['start']-1;b=m['end']
        if chrom not in bw.chroms():raise ValueError('Track chromosome naming does not match the request.')
        if choice['format']=='bigBed':
            entries=bw.entries(chrom,a,b) or []
            text='\n'.join(f'{chrom}\t{x}\t{y}\t{rest}' for x,y,rest in entries)
            features=parse_bed_track(text,ref,choice)
            return dict(features=features,signal=[],status='FETCHED_PEAKS',metadata=choice)
        n=min(800,b-a);values=bw.stats(chrom,a,b,nBins=n,exact=True)
        signal=[]
        for i,v in enumerate(values):
            s=1+(i*(b-a))//n;e=((i+1)*(b-a))//n
            if m.get('orientation','+')=='-':s,e=ref['length']-e+1,ref['length']-s+1
            signal.append(dict(start=s,end=e,value=float(v) if v is not None and math.isfinite(v) else None))
        return dict(features=[],signal=signal,status='FETCHED_SIGNAL',metadata=dict(choice,
          note='Exact mean signal per display interval; absent values remain missing, not zero.',
          loaded_utc=datetime.now(timezone.utc).isoformat()))
    finally:bw.close()
