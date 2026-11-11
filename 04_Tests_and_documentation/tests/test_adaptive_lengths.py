"""Regression tests for the user's adaptive-length / 45-55 C request.
Synthetic fixtures only. Engine-independent tests never claim full genome QC.
"""
import copy
import csv
import json
import random
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from design_engine import (Policy, generate_candidates, design, select_candidates,
                           summarize_specificity, enrich_hit, footprints_conflict)
from design_presets import PRESETS, LEGACY_PRESETS
from rna_avoidance import prepare_transcription, candidate_rna_safety
from sequence_tools import make_bins, validate_imagers
from test_design import cand, hit
BASE=Path(__file__).resolve().parents[1]
IMS=validate_imagers(json.loads((BASE/'imagers.json').read_text()))

def fixture(gc=0.32,n=650,seed=127):
    rng=random.Random(seed)
    seq=''.join(rng.choices('ACGT',weights=[(1-gc)/2,gc/2,gc/2,(1-gc)/2],k=n))
    r=dict(id='ADAPT',sequence=seq,length=n,mapping=None,features=[])
    cfg=dict(title='Synthetic adaptive test',references=[r],imagers=IMS,bin_width=n,
             assignments={'ADAPT_B001':'DP7'},transcription={'ADAPT':dict(mode='forward',confirmed=True)},
             policy={'bound_strand':'-','mode':'local_preview','shortlist_per_bin':40,'extension_shortlist_per_bin':35,
                     'target_per_bin':3,'minimum_per_bin':1})
    return cfg

def mined(cfg,updates=None):
    p=Policy.from_dict(dict(cfg['policy'],**(updates or {})))
    model,_,_=prepare_transcription(cfg['references'],cfg['transcription'])
    bins,_=make_bins(cfg['references'],cfg['bin_width'],assignments=cfg['assignments'])
    cs,counts=generate_candidates(cfg['references'],bins,IMS,p,lambda *a,**k:None,model)
    return cs,counts,p,model,bins

class AdaptivePolicyTests(unittest.TestCase):
    def test_user_requested_window(self):self.assertEqual((Policy().tm_min,Policy().tm_max),(45,55))
    def test_requested_defaults(self):
        p=Policy();self.assertEqual((p.length_mode,p.preferred_arm_max,p.arm_max),('adaptive',42,90));self.assertEqual(p.selection_goal,'maximize_sites')
    def test_other_scientific_gates_unchanged(self):
        p=Policy();self.assertEqual((p.high_offtarget_tm,p.moderate_offtarget_tm,p.folding_temperature),(42,35,42));self.assertEqual(p.minimum_gap,5);self.assertFalse(p.avoid_gene_rna);self.assertEqual(p.max_strong_rna_sites,0)
    def test_paper_preset_is_not_relabelled(self):
        p=LEGACY_PRESETS['ellenberg_informed']['values'];self.assertEqual((p['arm_min'],p['arm_max'],p['tm_min'],p['tm_max'],p['length_mode']),(36,42,42,50,'fixed'))
    def test_preferred_limit_must_fit_ceiling(self):
        with self.assertRaises(ValueError):Policy.from_dict({'arm_max':41})
    def test_bad_length_mode(self):
        with self.assertRaises(ValueError):Policy.from_dict({'length_mode':'unlimited'})
    def test_bad_selection_goal(self):
        with self.assertRaises(ValueError):Policy.from_dict({'selection_goal':'ignore_qc'})
    def test_shortlist_limit_counts_extensions(self):
        cfg=fixture()
        with self.assertRaisesRegex(ValueError,'exceed'):mined(cfg,{'max_total_candidates':1})
    def test_ceiling_is_enforced(self):
        cs,_,p,_,_=mined(fixture(),{'arm_max':60});self.assertTrue(cs);self.assertTrue(all(36<=c['arm_length']<=60 for c in cs))
    def test_extension_recovers_low_gc_candidates(self):
        cfg=fixture();cs,_,_,_,_=mined(cfg);short,_,_,_,_=mined(cfg,{'length_mode':'fixed','arm_max':42})
        self.assertTrue(any(c['arm_length']>42 for c in cs));self.assertGreater(len(cs),len(short))
    def test_shortest_rescue_retained(self):
        cfg=fixture(n=90,seed=128)
        cs,_,p,_,_=mined(cfg,{'shortlist_per_bin':160,'extension_shortlist_per_bin':160,'variants_per_start':1})
        from engine.oligominer_thermo_qc import calculate_tm,PROFILES
        from design_engine import exact_sites
        for c in cs:
            if c['generation_start']==1 and c['length_band']=='extended':
                # Every tested shorter extended window is out of the requested
                # temperature window or fails another explicit composition rule.
                for n in range(43,c['arm_length']):
                    w=cfg['references'][0]['sequence'][:n]
                    val=calculate_tm(w,PROFILES[p.profile])['adjusted_screen_C']
                    self.assertLessEqual(val,p.tm_min)
                return
        self.fail('Expected synthetic start-1 length rescue missing')
    def test_extension_obeys_new_tm_window(self):
        cs,_,p,_,_=mined(fixture());self.assertTrue(all(p.tm_min<c['tm']<p.tm_max for c in cs))
    def test_fixed_mode_never_extends(self):
        cs,counts,_,_,_=mined(fixture(gc=.5),{'length_mode':'fixed','arm_max':42});self.assertTrue(cs);self.assertFalse(any(c['length_band']=='extended' for c in cs));self.assertEqual(counts[0]['extension_shortlisted_for_QC'],0)
    def test_short_candidates_not_displaced_by_extended_band(self):
        cfg=fixture(gc=.5)
        cs,_,_,_,_=mined(cfg,{'arm_max':90});base,_,_,_,_=mined(cfg,{'arm_max':42})
        self.assertEqual({c['id'] for c in cs if c['length_band']=='preferred'},{c['id'] for c in base})
    def test_extended_candidates_keep_correct_RNA_strand(self):
        cs,_,_,model,_=mined(fixture());self.assertTrue(cs)
        self.assertTrue(all(c['bound_strand']=='-' and not candidate_rna_safety(c['sites'],model)['conflicts'] for c in cs))
    def test_ceiling_does_not_enlarge_region(self):
        cs,_,_,_,bins=mined(fixture(n=200));self.assertTrue(all(h['extended_start']>=1 and h['extended_end']<=200 for c in cs for h in c['sites']));self.assertEqual(bins[0]['length'],200)
    def test_extension_metadata_identifies_bases(self):
        cs,_,p,_,_=mined(fixture());long=[c for c in cs if c['length_band']=='extended'];self.assertTrue(long);self.assertTrue(all(c['extension_bases']==c['arm_length']-p.preferred_arm_max and c['extension_note'] for c in long))

class CountSelectionTests(unittest.TestCase):
    class Fold:
        def pair(self,a,b):return {'blocked':False}
    def setUp(self):
        self.bins=[dict(id='X_B001',reference='X',number=1,start=1,end=300,length=300,center=150.5,imager='DP7')]
    def test_maximize_mode_not_capped_by_benchmark(self):
        cs=[cand(str(i),5+i*50,40+i*50) for i in range(5)];p=Policy(mode='local_preview',target_per_bin=1,minimum_per_bin=1)
        sel=select_candidates(cs,self.bins,p,self.Fold(),lambda *a,**k:None);self.assertEqual(len(sel),5)
    def test_target_count_mode_remains_capped(self):
        cs=[cand(str(i),5+i*50,40+i*50) for i in range(5)];p=Policy(mode='local_preview',target_per_bin=1,minimum_per_bin=1,selection_goal='target_count')
        sel=select_candidates(cs,self.bins,p,self.Fold(),lambda *a,**k:None);self.assertEqual(len(sel),1)
    def test_QC_rejected_long_arm_cannot_rescue(self):
        a=cand('base');b=cand('extension',101,150);b.update(length_band='extended',reject_reasons=['RNA risk'])
        sel=select_candidates([a,b],self.bins,Policy(mode='local_preview'),self.Fold(),lambda *a,**k:None)
        self.assertEqual([c['id'] for c in sel],['base']);self.assertEqual(b['status'],'rejected_by_QC')
    def test_multiple_searches_preserve_baseline_site_count(self):
        a=cand('base');b=cand('long',5,90);b['length_band']='extended';c=cand('later',51,86)
        reports=[];sel=select_candidates([a,b,c],self.bins,Policy(mode='local_preview'),self.Fold(),lambda *a,**k:None,selection_report=reports)
        self.assertEqual(len(sel),2);self.assertEqual(len(reports),4);self.assertEqual(sum(x['chosen'] for x in reports),1)
        self.assertGreaterEqual(next(x['exact_sites'] for x in reports if x['chosen']),reports[0]['exact_sites'])
    def test_long_arm_still_fails_high_offtarget(self):
        c=cand('long');c.update(length_band='extended',arm_length=60);p=Policy(mode='local_preview')
        c['offtarget_hits']=[enrich_hit(hit(),c,p)];summarize_specificity(c,p);self.assertTrue(c['reject_reasons'])
        sel=select_candidates([c],self.bins,p,self.Fold(),lambda *a,**k:None);self.assertFalse(sel)
    def test_pair_block_still_enforced(self):
        class Block:
            def pair(self,a,b):return {'blocked':a['id']!=b['id']}
        cs=[cand('a'),dict(cand('b',100,135),length_band='extended')]
        sel=select_candidates(cs,self.bins,Policy(mode='local_preview',screen_cross_primary=True),Block(),lambda *a,**k:None);self.assertEqual(len(sel),1)
    def test_real_local_design_has_longer_selected_and_saved_evidence(self):
        cfg=fixture()
        with tempfile.TemporaryDirectory() as t:
            out=Path(t)/'result';r=design(copy.deepcopy(cfg),out,lambda *a,**k:None)
            self.assertGreater(r['counts']['selected'],0);self.assertGreater(r['counts']['extended_primaries'],0)
            self.assertEqual(r['reference_hashes_before'],r['reference_hashes_after'])
            self.assertEqual(r['status'],'LOCAL_ONLY_DESIGN_PREVIEW');self.assertEqual(r['policy']['tm_min'],45)
            self.assertTrue((out/'selection_search.csv').exists());self.assertIn('Arm_length_nt',(out/'selected_primaries.csv').read_text().splitlines()[0])
            selected=[c for c in r['candidates'] if c['status']=='selected']
            for i,c in enumerate(selected):
                self.assertFalse(c['reject_reasons']);self.assertEqual(c['specificity']['rna_risk_sites'],0)
                for other in selected[i+1:]:self.assertFalse(footprints_conflict(c['sites'],other['sites'],Policy.from_dict(r['policy'])))

if __name__=='__main__':unittest.main()
