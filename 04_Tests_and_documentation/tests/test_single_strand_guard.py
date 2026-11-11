"""Fixed-strand safety checks. Synthetic fixtures, not experimental QC."""
import copy
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from design_engine import Policy, generate_candidates, design
from strand_guard import require_fixed_strand, check_reference_polarity, audit_selected_strands
from design_presets import PRESETS
from sequence_tools import make_bins
from test_rna_direction import ref, genes, declaration, IMAGERS
from rna_avoidance import prepare_transcription


def primary(ident='P', strand='+', rid='Test', bid='Test_B001', start=20, end=55):
    return dict(id=ident, bound_strand=strand, bin_id=bid,
                sites=[dict(reference=rid, start=start, end=end, strand=strand)])


def mapped(rid, chromosome='chr1', orientation='+', start=1001):
    r=ref(n=100);r['id']=rid
    r['mapping']=dict(verified=True,assembly='GRCh38',chromosome=chromosome,
                      start=start,end=start+99,orientation=orientation)
    return r


class FixedPolicyTests(unittest.TestCase):
    def test_errfi1_top_is_default(self):
        self.assertEqual(Policy().bound_strand,'+')
    def test_missing_choice_uses_verified_default(self):
        self.assertEqual(Policy.from_dict({}).bound_strand,'+')
    def test_auto_is_rejected(self):
        with self.assertRaisesRegex(ValueError,'ONE fixed'):Policy.from_dict({'bound_strand':'auto'})
    def test_unknown_both_and_blank_rejected(self):
        for st in (None,'both','','top','unselected'):
            with self.subTest(strand=st),self.assertRaises(ValueError):Policy.from_dict({'bound_strand':st})
    def test_both_explicit_choices_allowed(self):
        for st in ('+','-'):self.assertEqual(Policy.from_dict({'bound_strand':st}).bound_strand,st)
    def test_disabling_rna_never_enables_auto(self):
        with self.assertRaises(ValueError):Policy.from_dict({'bound_strand':'auto','avoid_gene_rna':False,'screen_rna':False})
    def test_preset_restores_verified_strand(self):
        for preset in PRESETS.values():self.assertEqual(preset['values']['bound_strand'],'+')
    def test_direct_mining_cannot_bypass_guard(self):
        r=ref();bins,_=make_bins([r],600,assignments={'Test_B001':'DP7'})
        for st in ('auto','unselected'):
            with self.subTest(strand=st),self.assertRaises(ValueError):
                generate_candidates([r],bins,IMAGERS,Policy(mode='local_preview',bound_strand=st),lambda *a,**k:None,{})
    def test_api_style_request_stops_before_creating_output(self):
        for st in ('auto','unselected'):
            with tempfile.TemporaryDirectory() as t:
                out=Path(t)/'result'
                with self.assertRaises(ValueError):design({'policy':{'bound_strand':st}},out)
                self.assertFalse(out.exists())


class FixedMiningTests(unittest.TestCase):
    def mine(self,st):
        r=genes(ref(n=500),[('plus',1,190,'+'),('minus',251,500,'-')])
        model,_,_=prepare_transcription([r],declaration('annotations'))
        bins,_=make_bins([r],500,assignments={'Test_B001':'DP7'})
        p=Policy(mode='local_preview',bound_strand=st,shortlist_per_bin=20,extension_shortlist_per_bin=10,avoid_gene_rna=True)
        cs,counts=generate_candidates([r],bins,IMAGERS,p,lambda *a,**k:None,model)
        return cs,counts
    def test_top_lock_rejects_plus_genes_without_switching(self):
        cs,counts=self.mine('+');self.assertTrue(cs)
        self.assertEqual({c['bound_strand'] for c in cs},{'+'})
        self.assertTrue(all(c['generation_start']>190 for c in cs))
        self.assertGreater(counts[0]['RNA_direction_conflict'],0)
    def test_bottom_lock_rejects_minus_genes_without_switching(self):
        cs,counts=self.mine('-');self.assertTrue(cs)
        self.assertEqual({c['bound_strand'] for c in cs},{'-'})
        self.assertTrue(all(c['generation_end']<251 for c in cs))
        self.assertGreater(counts[0]['RNA_direction_conflict'],0)


class StrandAuditTests(unittest.TestCase):
    def setUp(self):
        self.refs=[ref()];self.bins=[dict(id='Test_B001',reference='Test'),dict(id='Test_B002',reference='Test')]
    def test_uniform_top_passes(self):
        a=audit_selected_strands([primary()],self.refs,self.bins,Policy(bound_strand='+'))
        self.assertEqual((a['status'],a['top_sites'],a['opposite_strand_sites']),('PASS_SEQUENCE_ORIENTATION',1,0))
    def test_uniform_bottom_passes(self):
        a=audit_selected_strands([primary(strand='-')],self.refs,self.bins,Policy(bound_strand='-'))
        self.assertEqual(a['bottom_sites'],1)
    def test_mixed_regions_rejected_even_if_each_region_is_uniform(self):
        cs=[primary(),primary('Q','-',bid='Test_B002',start=250,end=285)]
        with self.assertRaises(AssertionError):audit_selected_strands(cs,self.refs,self.bins,Policy(bound_strand='+'))
    def test_wrong_repeat_copy_rejected(self):
        c=primary();c['sites'].append(dict(reference='Test',start=250,end=285,strand='-'))
        with self.assertRaisesRegex(AssertionError,'repeat copy'):audit_selected_strands([c],self.refs,self.bins,Policy(bound_strand='+'))
    def test_empty_not_called_pass(self):
        a=audit_selected_strands([],self.refs,self.bins,Policy(bound_strand='+'))
        self.assertEqual(a['status'],'NO_SELECTED_PROBES')
        self.assertTrue(all(r['Status']=='NO_SELECTED_SITES' for r in a['per_region']))
    def test_missing_sites_not_called_pass(self):
        c=primary();c['sites']=[]
        with self.assertRaises(AssertionError):audit_selected_strands([c],self.refs,self.bins,Policy(bound_strand='+'))
    def test_reference_orientations_consistent(self):
        self.assertEqual(check_reference_polarity([mapped('A'),mapped('B',start=2001)],Policy(bound_strand='+')),{'NC_000001.11':'+'})
    def test_reversed_reference_changes_physical_strand_label(self):
        self.assertEqual(check_reference_polarity([mapped('A',orientation='-')],Policy(bound_strand='+')),{'NC_000001.11':'-'})
    def test_mixed_orientations_same_chromosome_rejected(self):
        with self.assertRaisesRegex(ValueError,'inconsistent sequence orientations'):
            check_reference_polarity([mapped('A'),mapped('B',orientation='-')],Policy(bound_strand='+'))
    def test_unknown_orientation_not_guessed(self):
        with self.assertRaises(ValueError):check_reference_polarity([mapped('A',orientation='?')],Policy(bound_strand='+'))
    def test_independent_chromosomes_do_not_claim_same_molecule(self):
        a=check_reference_polarity([mapped('A'),mapped('B',chromosome='chr2',orientation='-')],Policy(bound_strand='+'))
        self.assertEqual(a,{'NC_000001.11':'+','NC_000002.12':'-'})
    def test_correctly_mapped_sites_pass(self):
        r=mapped('Test',orientation='-');c=primary()
        a=audit_selected_strands([c],[r],self.bins,Policy(bound_strand='+'))
        self.assertEqual(a['genomic_strand_by_chromosome'],{'NC_000001.11':'-'})
    def test_final_audit_runs_before_order_export(self):
        cfg=self.config()
        def corrupt(*args,**kwargs):return [primary(strand='-')]
        with tempfile.TemporaryDirectory() as t,patch('design_engine.select_candidates',side_effect=corrupt):
            with self.assertRaisesRegex(AssertionError,'wrong/missing primary strand'):
                design(cfg,Path(t)/'result',lambda *a,**k:None)
            self.assertFalse((Path(t)/'result/DRAFT_IDT_oPools.csv').exists())
    def config(self):
        r=ref(n=260)
        return dict(references=[r],imagers=IMAGERS,assignments={'Test_B001':'DP7'},bin_width=260,
                    transcription=declaration('reverse'),policy=dict(mode='local_preview',bound_strand='+',
                    shortlist_per_bin=8,extension_shortlist_per_bin=6,variants_per_start=1))
    def test_real_local_run_writes_matching_json_and_csv(self):
        cfg=self.config()
        with tempfile.TemporaryDirectory() as t:
            r=design(cfg,Path(t)/'result',lambda *a,**k:None)
            a=json.loads((Path(t)/'result/strand_consistency.json').read_text())
            self.assertEqual(a,r['strand_audit']);self.assertTrue(a['selected_primaries'])
            self.assertEqual(a['bottom_sites'],0)
            self.assertTrue((Path(t)/'result/strand_consistency.csv').exists())
            self.assertTrue(any(x['name']=='One fixed DNA-bound strand across ALL tracing regions' for x in r['checks']))


if __name__=='__main__':unittest.main()
