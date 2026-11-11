"""Sequence convention, round-trip, assembly and export regression tests.
All new named sequences in these tests are fixtures, not recommended imagers.
"""
import copy
import csv
import json
import random
import sys
import tempfile
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from sequence_tools import resolve_imager_entry, validate_imagers, rc
from design_engine import design, Policy
from test_adaptive_lengths import fixture

BASE=Path(__file__).resolve().parents[1]
I='CCGATACGAATA'
D='TATTCGTATCGG'

def entry(kind='imager', value=None, **extra):
    return dict(id='TEST_DP7', input_kind=kind,
                input_sequence_5to3=value or (I if kind=='imager' else D), **extra)

class ImagerConventionTests(unittest.TestCase):
    def test_imager_input(self):
        r=resolve_imager_entry(entry());self.assertEqual(r['sequence'],I);self.assertEqual(r['docking'],D)
    def test_docking_input(self):
        r=resolve_imager_entry(entry('docking'));self.assertEqual(r['sequence'],I);self.assertEqual(r['docking'],D)
    def test_explicit_sequence_names(self):
        r=resolve_imager_entry(entry('docking'));self.assertEqual(r['imager_sequence_5to3'],I);self.assertEqual(r['docking_sequence_5to3'],D)
    def test_provenance_docking(self):
        r=resolve_imager_entry(entry('docking'));self.assertEqual((r['input_kind'],r['input_sequence_5to3']),('docking',D))
    def test_provenance_imager(self):
        r=resolve_imager_entry(entry());self.assertEqual((r['input_kind'],r['input_sequence_5to3']),('imager',I))
    def test_legacy_sequence(self):
        r=resolve_imager_entry({'id':'old','sequence':I});self.assertEqual(r['docking'],D);self.assertEqual(r['input_kind'],'imager')
    def test_docking_only_explicit_key(self):
        r=resolve_imager_entry({'id':'dock','docking_sequence_5to3':D});self.assertEqual(r['sequence'],I);self.assertEqual(r['input_kind'],'docking')
    def test_imager_only_explicit_key(self):
        r=resolve_imager_entry({'id':'im','imager_sequence_5to3':I});self.assertEqual(r['sequence'],I)
    def test_missing_type_on_unlabeled_entered_sequence_rejected(self):
        with self.assertRaises(ValueError):resolve_imager_entry({'id':'X','input_sequence_5to3':I})
    def test_unknown_type_rejected(self):
        with self.assertRaises(ValueError):resolve_imager_entry(entry('target_dna'))
    def test_blank_type_not_inferred(self):
        with self.assertRaises(ValueError):resolve_imager_entry(entry(''))
    def test_no_sequence(self):
        with self.assertRaises(ValueError):resolve_imager_entry({'id':'X','input_kind':'imager'})
    def test_reversed_direction_metadata_rejected(self):
        with self.assertRaises(ValueError):resolve_imager_entry(entry(sequence_orientation='3to5'))
    def test_bad_name(self):
        bad=entry();bad['id']='not a valid name'
        with self.assertRaises(ValueError):resolve_imager_entry(bad)
    def test_normalization(self):
        r=resolve_imager_entry(entry('docking',' tat tcg\ntatcgg '));self.assertEqual(r['docking'],D)
    def test_reject_invalid_characters(self):
        for value in ('ACGTNNNNACGT','ACGTUACGTACG',"5'-"+I+"-3'",I+'123',I+'/ATTO647/',I+'-',I+'\u200b'):
            with self.subTest(value=value),self.assertRaises(ValueError):resolve_imager_entry(entry('imager',value))
    def test_length_bounds(self):
        for n in (7,41):
            with self.subTest(n=n),self.assertRaises(ValueError):resolve_imager_entry(entry('imager','A'*n))
        for n in (8,40):self.assertEqual(len(resolve_imager_entry(entry('imager','A'*n))['docking']),n)
    def test_wrong_cached_docking_rejected(self):
        with self.assertRaises(ValueError):resolve_imager_entry(dict(id='X',sequence=I,docking='AAAAAAAAAAAA'))
    def test_wrong_cached_imager_rejected(self):
        with self.assertRaises(ValueError):resolve_imager_entry(entry('docking',sequence='AAAAAAAAAAAA'))
    def test_duplicate_pair_across_conventions(self):
        a=entry();b=entry('docking');b['id']='OTHER'
        with self.assertRaises(ValueError):validate_imagers([a,b])
    def test_duplicate_id(self):
        with self.assertRaises(ValueError):validate_imagers([entry(),entry('imager','ACGTACGTACGT')])
    def test_not_same_as_reverse_imager(self):
        a=entry();b=entry('imager',D);b['id']='OPPOSITE'
        out=validate_imagers([a,b]);self.assertEqual(len(out),2);self.assertNotEqual(out[0]['sequence'],out[1]['sequence'])
    def test_no_in_place_input_mutation(self):
        original=entry('docking');snap=copy.deepcopy(original);resolve_imager_entry(original);self.assertEqual(original,snap)
    def test_idempotent_serialization(self):
        for kind in ('imager','docking'):
            once=resolve_imager_entry(entry(kind));self.assertEqual(resolve_imager_entry(json.loads(json.dumps(once))),once)
    def test_many_pairs_same_canonical_result(self):
        rng=random.Random(1213)
        for _ in range(200):
            seq=''.join(rng.choice('ACGT') for _ in range(rng.randint(8,40)))
            a=resolve_imager_entry(entry('imager',seq));b=resolve_imager_entry(entry('docking',rc(seq)))
            self.assertEqual((a['sequence'],a['docking']),(b['sequence'],b['docking']))
    def test_all_inventory_pairs_both_modes(self):
        for original in json.loads((BASE/'imagers.json').read_text()):
            a=resolve_imager_entry(original);b=resolve_imager_entry(dict(id=original['id'],input_kind='docking',input_sequence_5to3=rc(original['sequence'])))
            self.assertEqual((a['sequence'],a['docking']),(b['sequence'],b['docking']))
    def test_defaults_not_changed(self):
        p=Policy();self.assertEqual((p.tm_min,p.tm_max,p.preferred_arm_max,p.arm_max,p.minimum_gap,p.length_mode),(45,55,42,90,5,'adaptive'))
        self.assertFalse(p.avoid_gene_rna)

class DesignPolarityIntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp=tempfile.TemporaryDirectory();cls.root=Path(cls.tmp.name);cls.runs=[]
        cfg=fixture(n=260);cfg['policy'].update(shortlist_per_bin=8,extension_shortlist_per_bin=8,variants_per_start=1)
        # Hold all scientific settings and the physical imager constant.
        for kind in ('imager','docking'):
            c=copy.deepcopy(cfg)
            c['imagers']=[dict(id=i['id'],input_kind=kind,input_sequence_5to3=i['sequence'] if kind=='imager' else rc(i['sequence']),dye=i.get('dye','not specified')) for i in cfg['imagers']]
            r=design(c,cls.root/kind,lambda *a,**k:None);cls.runs.append(r)
    @classmethod
    def tearDownClass(cls):cls.tmp.cleanup()
    def test_same_generated_physical_primaries(self):
        self.assertEqual([c['full'] for c in self.runs[0]['candidates']],[c['full'] for c in self.runs[1]['candidates']])
    def test_same_selected_pool(self):
        self.assertTrue(self.runs[0]['selected_ids']);self.assertEqual(self.runs[0]['selected_ids'],self.runs[1]['selected_ids'])
    def test_same_draft_order_bytes(self):
        self.assertEqual((self.root/'imager/DRAFT_IDT_oPools.csv').read_bytes(),(self.root/'docking/DRAFT_IDT_oPools.csv').read_bytes())
    def test_canonical_used_for_real_assembly(self):
        for r in self.runs:
            for c in r['candidates']:
                d=rc(c['imager_sequence']);self.assertEqual(c['full'],d+'TT'+c['arm']+'TT'+d);self.assertEqual(c['docking_sequence'],d)
    def test_input_provenance_recorded_in_candidates(self):
        for r,kind in zip(self.runs,('imager','docking')):
            for c in r['candidates']:
                self.assertEqual(c['imager_input_kind'],kind)
                self.assertEqual(c['imager_entered_sequence_5to3'],c['imager_sequence'] if kind=='imager' else c['docking_sequence'])
    def test_unambiguous_selected_export(self):
        with (self.root/'docking/selected_primaries.csv').open() as f:rows=list(csv.DictReader(f))
        self.assertTrue(rows)
        for row in rows:
            self.assertEqual(row['Imager_Input_Type'],'docking')
            self.assertEqual(row['Docking_Sequence_in_Primary_5to3'],rc(row['Imager_Sequence_5to3']))
    def test_inventory_export_has_both_sequences(self):
        with (self.root/'docking/imager_inventory.csv').open() as f:rows=list(csv.DictReader(f))
        self.assertEqual(rows[0]['Entered_Type'],'docking');self.assertEqual(rows[0]['Entered_Sequence_5to3'],rows[0]['Docking_Sequence_in_Primary_5to3'])
    def test_readout_assignment_export(self):
        with (self.root/'docking/readout_assignments.csv').open() as f:rows=list(csv.DictReader(f))
        self.assertEqual(len(rows),self.runs[1]['counts']['regions']);self.assertIn('Imager_Sequence_5to3',rows[0])
    def test_legacy_inputs_supported(self):
        for r in self.runs:self.assertEqual(r['reference_hashes_before'],r['reference_hashes_after'])

if __name__=='__main__':unittest.main()
