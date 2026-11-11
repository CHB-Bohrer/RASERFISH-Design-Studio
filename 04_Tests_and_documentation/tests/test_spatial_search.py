"""Spatial-search regressions; synthetic fixtures do not certify genome QC."""
import copy
import random
import threading
import unittest
from unittest.mock import patch

from design_engine import (Cancelled, Folding, Policy, _selection_pass,
                           _spatial_shortlist, footprints_conflict, select_candidates)
from sequence_tools import region_metrics
from test_design import cand


def region(length=20000,start=1,number=1):
    return dict(id=f'X_B{number:03d}',reference='X',number=number,start=start,
                end=start+length-1,length=length,center=start+(length-1)/2,imager='DP1')


def candidate(ident,start,bin_id='X_B001'):
    c=cand(ident,start,start+35,tm=49,bin_id=bin_id)
    c.update(arm_length=36,length_band='preferred')
    return c


class CompatibleFold:
    def pair(self,a,b):
        return {'blocked':False}


class SpatialShortlistTests(unittest.TestCase):
    def test_twenty_fixed_one_kb_sectors_receive_equal_budgets(self):
        b=region();pool=[candidate(f'{sector}-{j}',1+1000*sector+20*j)
                         for sector in range(20) for j in range(48)]
        out=_spatial_shortlist(pool,b,400,True,50)
        counts=[sum((c['start']-1)//1000==i for c in out) for i in range(20)]
        self.assertEqual(counts,[20]*20)

    def test_short_final_sector_is_real_not_rescaled(self):
        b=region(2500,start=189019)
        pool=[candidate(str(j),b['start']+j) for j in range(0,2460,20)]
        out=_spatial_shortlist(pool,b,6,True,50)
        counts=[sum((c['start']-b['start'])//1000==i for c in out) for i in range(3)]
        self.assertEqual(counts,[2,2,2])

    def test_near_identical_ranked_alternates_do_not_crowd_out_spatial_options(self):
        b=region(1000)
        pool=[candidate(f'cluster-{j:03d}',490+j%10) for j in range(150)]
        pool += [candidate(str(pos),pos) for pos in (50,250,750,950)]
        # The dispersed options have worse Tm distance and a longer arm. They
        # still get search coverage; they must subsequently pass unchanged QC.
        for c in pool[-4:]:c.update(tm=45.1,arm_length=42)
        out=_spatial_shortlist(pool,b,5,True,50)
        self.assertEqual({c['start'] for c in out if not c['id'].startswith('cluster')},{50,250,750,950})

    def test_budget_smaller_than_sector_count_covers_entire_interval(self):
        b=region();pool=[candidate(str(i),1000*i+500) for i in range(20)]
        out=_spatial_shortlist(pool,b,4,True,50)
        self.assertEqual([(c['start']-1)//1000 for c in out],[2,7,12,17])

    def test_shortlist_is_deterministic_under_input_reordering(self):
        b=region(3000);pool=[candidate(str(j),j) for j in range(1,2960,11)]
        a=[c['id'] for c in _spatial_shortlist(pool,b,31,True,50)]
        random.Random(17).shuffle(pool)
        self.assertEqual(a,[c['id'] for c in _spatial_shortlist(pool,b,31,True,50)])


class DistributionMetricsTests(unittest.TestCase):
    def test_gap_includes_edges_and_spacing_uses_site_centers(self):
        metrics=region_metrics(region(2000),[candidate('a',101),candidate('b',1101)])
        self.assertEqual(metrics['sector_counts'],[1,1])
        self.assertEqual(metrics['sector_width_bp'],1000)
        self.assertEqual(metrics['sector_count'],2)
        self.assertEqual(metrics['largest_gap_bp'],964)
        self.assertEqual(metrics['largest_gap_bp'],metrics['longest_gap'])
        self.assertEqual(metrics['median_spacing_bp'],1000)
        self.assertEqual(len(metrics['distribution']),10)

    def test_empty_region_has_full_gap_and_no_invented_spacing(self):
        metrics=region_metrics(region(2500),[])
        self.assertEqual(metrics['largest_gap_bp'],2500)
        self.assertIsNone(metrics['median_spacing_bp'])
        self.assertEqual(metrics['sector_counts'],[0,0,0])

    def test_repeat_sites_count_independently(self):
        a=candidate('repeat',101);a['sites'].append(dict(a['sites'][0],start=1101,end=1136,extended_start=1101,extended_end=1136))
        metrics=region_metrics(region(2000),[a])
        self.assertEqual((metrics['probes'],metrics['exact_sites']),(1,2))
        self.assertEqual(metrics['sector_counts'],[1,1])


class RepairTests(unittest.TestCase):
    class HubFold:
        def pair(self,a,b):
            return {'blocked':a['id']!=b['id'] and 'hub' in (a['id'],b['id'])}

    def fixture(self):
        hub=candidate('hub',1)
        hub['sites'].append(dict(hub['sites'][0],start=400,end=435,extended_start=400,extended_end=435))
        hub['end']=435
        return [hub]+[candidate(name,pos) for name,pos in zip(('b','c','d'),(60,160,260))]

    def select(self,cs,bins=None,fold=None,**policy):
        reports=[]
        selected=select_candidates(cs,bins or [region(500)],Policy(mode='local_preview',screen_cross_primary=True,**policy),
                                   fold or self.HubFold(),lambda *a,**k:None,selection_report=reports)
        return selected,reports

    def test_repair_replaces_two_site_hub_with_three_compatible_sites(self):
        selected,reports=self.select(self.fixture())
        chosen=next(r for r in reports if r['chosen'])
        self.assertEqual({c['id'] for c in selected},{'b','c','d'})
        self.assertEqual(chosen['pre_repair_exact_sites'],2)
        self.assertEqual(chosen['exact_sites'],3)
        self.assertGreater(chosen['repair_accepted'],0)
        self.assertLessEqual(chosen['repair_attempts'],32)
        self.assertLessEqual(chosen['repair_seed_checks'],128)
        self.assertIn('no global-optimum',chosen['search_method'])

    def test_gap_improvement_cannot_cost_exact_sites(self):
        selected,reports=self.select(self.fixture()[:2])
        self.assertEqual([c['id'] for c in selected],['hub'])
        self.assertEqual(next(r for r in reports if r['chosen'])['exact_sites'],2)

    def test_QC_rejected_alternate_cannot_be_used_for_repair(self):
        cs=self.fixture();cs[-1]['reject_reasons']=['High-Tm off-target']
        selected,_=self.select(cs)
        self.assertNotIn('d',{c['id'] for c in selected})
        self.assertEqual(cs[-1]['status'],'rejected_by_QC')
        self.assertTrue(all(not c['reject_reasons'] for c in selected))

    def test_result_and_repairs_are_deterministic(self):
        cs=self.fixture();a,ra=self.select(copy.deepcopy(cs))
        random.Random(98).shuffle(cs);b,rb=self.select(copy.deepcopy(cs))
        self.assertEqual([c['id'] for c in a],[c['id'] for c in b])
        self.assertEqual(ra,rb)

    def test_target_count_cap_is_preserved(self):
        selected,_=self.select([candidate(str(i),40+60*i) for i in range(6)],fold=CompatibleFold(),
                               selection_goal='target_count',target_per_bin=2,minimum_per_bin=1)
        self.assertEqual(len(selected),2)

    def test_repair_keeps_global_cross_region_pair_and_spacing_constraints(self):
        bins=[region(500),region(500,start=501,number=2)]
        cs=self.fixture()
        cs += [candidate('adjacent_left',465),candidate('adjacent_right',501,'X_B002'),candidate('remote',700,'X_B002')]
        calls=[]
        class Fold(self.HubFold):
            def pair(self,a,b):
                calls.append((a['bin_id'],b['bin_id']))
                blocked=super().pair(a,b)['blocked'] or {a['id'],b['id']}=={'b','remote'}
                return {'blocked':blocked}
        fold=Fold();selected,_=self.select(cs,bins,fold)
        for i,c in enumerate(selected):
            for other in selected[i+1:]:
                self.assertFalse(footprints_conflict(c['sites'],other['sites'],Policy()))
                self.assertFalse(fold.pair(c,other)['blocked'])
        self.assertTrue(any(a!=b for a,b in calls))


class GateAndCancellationTests(unittest.TestCase):
    def test_actual_cofold_gate_keeps_equality_and_blocks_exceedance(self):
        fold=Folding(Policy(mode='local_preview',screen_cross_primary=True,minimum_cofold_gain=-12))
        a=candidate('a',1);b=candidate('b',100)
        a['full']='AAAA';b['full']='TTTT'
        class Compound:
            energy=-12.0
            def mfe(self):return '((....))',self.energy
        class Engine:
            def fold_compound(self,*args):return Compound()
        fold.RNA=Engine();fold.fold=lambda c:{'mfe_kcal_mol':0.0}
        with patch('design_engine.exact_pair_segments',return_value=[]):
            self.assertFalse(fold.pair(a,b)['blocked'])
            fold.pairs.clear();Compound.energy=-12.0001
            self.assertTrue(fold.pair(a,b)['blocked'])

    def test_actual_exact_stem_gate_keeps_equality_and_blocks_exceedance(self):
        fold=Folding(Policy(mode='local_preview',screen_cross_primary=True,maximum_pair_stem=11))
        a=candidate('a',1);b=candidate('b',100)
        a['full']='A'*11;b['full']='T'*11
        self.assertFalse(fold.pair(a,b)['blocked'])
        a=dict(a,id='a12',full='A'*12);b=dict(b,id='b12',full='T'*12)
        self.assertTrue(fold.pair(a,b)['blocked'])

    def test_extended_footprint_five_base_gap_is_unchanged(self):
        a=candidate('a',5);b=candidate('b',46)
        self.assertFalse(footprints_conflict(a['sites'],b['sites'],Policy()))
        b['sites'][0]['extended_start']=45
        self.assertTrue(footprints_conflict(a['sites'],b['sites'],Policy()))

    def test_cancelled_empty_search_does_not_report_success(self):
        event=threading.Event();event.set()
        with self.assertRaises(Cancelled):
            select_candidates([],[region()],Policy(mode='local_preview'),CompatibleFold(),cancel=event)

    def test_cancellation_between_pair_checks_is_honored(self):
        event=threading.Event();calls=[]
        class Fold:
            def pair(self,a,b):
                calls.append((a['id'],b['id']))
                if len(calls)==3:event.set()
                return {'blocked':False}
        with self.assertRaises(Cancelled):
            _selection_pass([candidate(str(i),40+70*i) for i in range(6)],[region(500)],
                            Policy(mode='local_preview'),Fold(),'packing',cancel=event)
        self.assertLessEqual(len(calls),3)

    def test_cancellation_entering_repair_stops_before_any_swap(self):
        event=threading.Event();fixture=RepairTests().fixture()
        def log(message,**kwargs):
            if message.startswith('Repairing'):event.set()
        with self.assertRaises(Cancelled):
            select_candidates(fixture,[region(500)],Policy(mode='local_preview'),
                              RepairTests.HubFold(),log,cancel=event)


if __name__=='__main__':unittest.main()
