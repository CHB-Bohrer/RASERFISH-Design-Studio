"""Exact parity and isolation for bounded alignment-computation caches."""
import copy
from dataclasses import replace
from pathlib import Path
import random
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import design_engine as engine
import sequence_tools


def original_tracts(q,s):
    runs=[];start=None
    for i,(a,b) in enumerate(zip(q,s)):
        if a==b and a in 'ACGT':
            if start is None:start=i
        elif start is not None:runs.append(q[start:i]);start=None
    if start is not None:runs.append(q[start:])
    return runs


def original_enrich(h,c,p):
    runs=original_tracts(h['qseq'].upper(),h['sseq'].upper())
    h['exact_tract']=max(map(len,runs),default=0)
    h['offtarget_tm'],h['tm_model']=engine.duplex_estimate(h['qseq'].upper(),h['sseq'].upper(),p.profile,h.get('rna',False))
    exact_tms=[]
    if not h.get('rna'):
        exact_tms=[engine.duplex_estimate(s,s,p.profile,False)[0] for s in runs if len(s)>=8]
    h['exact_tract_tm']=max((v for v in exact_tms if v is not None),default=None)
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


def hit(q,s,**kw):
    return dict(qseq=q,sseq=s,rna=False,identity=100.,query_coverage=100.,alignment_length=len(q),**kw)


PAIRS=[
    ('GTAATCTCAGCACTTTGGGAGGCCAAGGCGGGCAGA',)*2,
    ('GCGCATGC'*5,'GCGCATGC'*5),
    ('GCGCATGC'*5,'GCGCATGT'+'GCGCATGC'*4),
    ('GCGCATGC'*5,'GCGCATGG'+'GCGCATGC'*4),
    ('GCGCATGC'*5,'GCGCATG-'+'GCGCATGC'*4),
    ('GCGCATGC'*5,'GCGCATGN'+'GCGCATGC'*4),
    ('NACGTACGTACGTNACGTACGT','NACGTACGTACGTNACGTACGT'),
    ('acgtacgtacgt','ACGTACGTACGT'),
    ('ACGTACGTACGT','ACGTACGT'),
    ('ACGTACGT','ACGTACGTACGT'),
    ('',''),('A','A'),('ACG-T','ACGGT'),('NNNNNNNN','NNNNNNNN'),
]


class AlignmentPerformanceTests(unittest.TestCase):
    def setUp(self):
        engine.clear_alignment_caches()
        sequence_tools._canonical_ref_cached.cache_clear()

    def test_tract_exact_parity_and_public_list_isolation(self):
        rng=random.Random(7)
        pairs=PAIRS+[( ''.join(rng.choices('ACGTN-',k=rng.randrange(0,121))),
                       ''.join(rng.choices('ACGTN-',k=rng.randrange(0,121)))) for _ in range(500)]
        for q,s in pairs:
            expected=original_tracts(q,s)
            self.assertEqual(engine.perfect_tracts(q,s),expected)
            got=engine.perfect_tracts(q,s);got.append('MUTATED')
            self.assertEqual(engine.perfect_tracts(q,s),expected)
            metrics=engine._perfect_tract_metrics(q,s)
            self.assertIsInstance(metrics,tuple);self.assertIsInstance(metrics[1],tuple)
            self.assertEqual(metrics[0],max(map(len,expected),default=0))

    def test_enriched_dictionaries_match_exactly_for_all_profiles_and_rna(self):
        for profile in engine.PROFILES:
            for rna in (False,True):
                for q,s in PAIRS:
                    for identity,coverage,length in [(79.,69.,29),(80.,70.,30),(89.9,89.9,29),(90.,90.,30),(100.,100.,len(q))]:
                        h=hit(q,s);h.update(rna=rna,identity=identity,query_coverage=coverage,alignment_length=length)
                        p=engine.Policy(profile=profile)
                        self.assertEqual(engine.enrich_hit(copy.deepcopy(h),{'tm':49.25},p),
                                         original_enrich(copy.deepcopy(h),{'tm':49.25},p))

    def test_candidate_and_policy_values_are_not_cached(self):
        q='GCGCATGC'*5;s=q;h=hit(q,s)
        base=engine.Policy()
        for p in [base,replace(base,high_offtarget_tm=-100,moderate_offtarget_tm=-99),
                  replace(base,high_offtarget_tm=1000,moderate_offtarget_tm=999),
                  replace(base,rna_exact_min=6),replace(base,rna_exact_min=200)]:
            for rna in (False,True):
                for tm in (10,45,55,100):
                    for identity,coverage,length in [(0,0,0),(89,89,29),(90,90,30)]:
                        h.update(rna=rna,identity=identity,query_coverage=coverage,alignment_length=length)
                        self.assertEqual(engine.enrich_hit(dict(h),{'tm':tm},p),original_enrich(dict(h),{'tm':tm},p))
        a=engine.enrich_hit(dict(h),{'tm':50},base);a['risk_tm']='MUTATED'
        self.assertEqual(engine.enrich_hit(dict(h),{'tm':50},base),original_enrich(dict(h),{'tm':50},base))

    def test_cache_eviction_is_bounded_and_same_answer_returns(self):
        q,s=PAIRS[0];expected=engine.enrich_hit(hit(q,s),{'tm':49},engine.Policy())
        for i in range(8300):
            # Ambiguous separators keep tract thermodynamics cheap while creating
            # distinct valid cache keys. This exercises eviction, not NN runtime.
            seq='N'.join('ACGT'[(i//(4**k))%4] for k in range(7))
            engine.enrich_hit(hit(seq,seq),{'tm':49},engine.Policy())
        self.assertEqual(engine._perfect_tract_metrics.cache_info().currsize,8192)
        self.assertEqual(engine._exact_tract_tm.cache_info().currsize,8192)
        self.assertLessEqual(engine.duplex_estimate.cache_info().currsize,8192)
        self.assertEqual(engine.enrich_hit(hit(q,s),{'tm':49},engine.Policy()),expected)

    def test_cache_clear_releases_all_sequence_evidence(self):
        q,s=PAIRS[0];engine.enrich_hit(hit(q,s),{'tm':49},engine.Policy())
        engine.clear_alignment_caches()
        for func in [engine._perfect_tract_metrics,engine._exact_tract_tm,engine.duplex_estimate]:
            self.assertEqual(func.cache_info().currsize,0)

    def test_design_exit_releases_sequence_caches_after_success_and_failure(self):
        for checkpointed in (False,True):
            for fails in (False,True):
                with tempfile.TemporaryDirectory() as directory:
                    root=Path(directory)
                    def stub(*args):
                        q,s=PAIRS[0]
                        engine.enrich_hit(hit(q,s),{'tm':49},engine.Policy())
                        self.assertGreater(engine._perfect_tract_metrics.cache_info().currsize,0)
                        if fails:raise ValueError('synthetic design failure')
                        return {'status':'stub completed'}
                    with patch.object(engine,'_design',side_effect=stub):
                        opts={'checkpoint_dir':root/'checkpoint'} if checkpointed else {}
                        if fails:
                            with self.assertRaisesRegex(ValueError,'synthetic design failure'):
                                engine.design({},root/'result',**opts)
                        else:self.assertEqual(engine.design({},root/'result',**opts),{'status':'stub completed'})
                    for func in [engine._perfect_tract_metrics,engine._exact_tract_tm,engine.duplex_estimate]:
                        self.assertEqual(func.cache_info().currsize,0)

    def test_error_type_message_and_partial_assignments_are_unchanged(self):
        for q,s in [PAIRS[0],('ACGTACGT-NACGTACGT','ACGTACGT-NACGTACGT'),('A','A')]:
            old=hit(q,s);new=hit(q,s);p=engine.Policy(profile='unknown-profile')
            errors=[]
            for fn,h in [(original_enrich,old),(engine.enrich_hit,new)]:
                try:fn(h,{'tm':49},p)
                except Exception as exc:errors.append((type(exc),str(exc)))
            self.assertEqual(len(errors)%2,0)
            if errors:self.assertEqual(errors[0],errors[1])
            self.assertEqual(old,new)

    def test_intended_location_classification_remains_per_hit(self):
        ref={'id':'r','sequence':'A'*100,'mapping':{'chromosome':'chr1','start':1001,'end':1100,'strand':'+','verified':True}}
        bins={'b':{'reference':'r','start':1,'end':100}}
        c={'bin_id':'b','sites':[{'strand':'+'}]}
        for reference in ['chr1','1','ref|NC_000001.11|','NC_000001.11','NT_187693.1']:
            for strand in ['+','-']:
                for a,b in [(1001,1040),(1061,1100),(1000,1039),(1062,1101)]:
                    h={'rna':False,'source':'human','reference':reference,'bound_strand':strand,'start':a,'end':b}
                    expected=reference!='NT_187693.1' and strand=='+' and a>=1001 and b<=1100
                    self.assertEqual(engine.expected_hit(h,c,bins,{'r':ref},[]),expected)
                    h['rna']=True;self.assertFalse(engine.expected_hit(h,c,bins,{'r':ref},[]))

    def test_canonical_cache_aliases_alt_loci_and_eviction(self):
        def original(s):
            s=s.strip()
            if s.startswith('ref|') and s.endswith('|'):s=s[4:-1]
            if s in sequence_tools.CHROMOSOMES:return sequence_tools.CHROMOSOMES[s]
            if 'chr'+s in sequence_tools.CHROMOSOMES:return sequence_tools.CHROMOSOMES['chr'+s]
            return s
        names=['chr1','1','ref|NC_000001.11|','NC_000001.10','NT_187693.1','chrM','M',' MT ','ref|NT_187693.1|','ref|x','']
        for s in names:self.assertEqual(sequence_tools.canonical_ref(s),original(s))
        for i in range(4200):self.assertEqual(sequence_tools.canonical_ref('alt_'+str(i)),'alt_'+str(i))
        self.assertEqual(sequence_tools._canonical_ref_cached.cache_info().currsize,4096)
        for s in names:self.assertEqual(sequence_tools.canonical_ref(s),original(s))
        with self.assertRaises(AttributeError):sequence_tools.canonical_ref(None)


if __name__=='__main__':unittest.main()
