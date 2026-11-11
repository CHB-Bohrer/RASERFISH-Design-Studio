"""Portable saved-run and scientific-code regression checks; no network or full QC."""
import ast
import copy
import gzip
import hashlib
import io
import json
import sys
import unittest
import zipfile
from dataclasses import asdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from design_engine import Policy
from design_presets import PRESETS, preset_provenance
from sequence_tools import make_bins, rc, validate_imagers


def read_json(relative):
    return json.loads((ROOT / relative).read_text())


class ERRFI1DefaultsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.baseline = read_json('verification/errfi1_baseline.json')
        cls.saved = read_json('inputs/errfi1_defaults.json')
        with gzip.open(ROOT / 'inputs/errfi1_reference.json.gz', 'rt') as handle:
            cls.references = json.load(handle)
        cls.expected_policy = dict(cls.baseline['policy_exact'])
        cls.expected_policy.update(shortlist_per_bin=400,
                                   extension_shortlist_per_bin=200,
                                   preset_id='errfi1', threads=2, screen_cross_primary=False)

    def test_backend_defaults_match_verified_run_and_only_authorized_changes(self):
        self.assertEqual(asdict(Policy()), self.expected_policy)
        self.assertEqual(asdict(Policy.from_dict({})), self.expected_policy)

    def test_saved_defaults_and_current_preset_are_complete_and_identical(self):
        self.assertEqual(self.saved['policy'], self.expected_policy)
        self.assertIn('errfi1', PRESETS)
        for name, preset in PRESETS.items():
            with self.subTest(preset=name):
                self.assertEqual(preset['values'], self.expected_policy)
        self.assertEqual(preset_provenance(self.expected_policy)['modified_fields'], {})

    def test_every_scientific_threshold_and_mode_matches_run(self):
        for field, value in self.baseline['scientific_policy_exact'].items():
            with self.subTest(field=field):
                self.assertEqual(getattr(Policy(), field), value)
                self.assertEqual(self.saved['policy'][field], value)

    def test_bundled_reference_sequence_and_mapping_match_saved_run(self):
        metadata = {r['id']: r for r in self.baseline['reference_metadata']}
        self.assertEqual({r['id'] for r in self.references}, set(metadata))
        for reference in self.references:
            expected = metadata[reference['id']]
            with self.subTest(reference=reference['id']):
                self.assertEqual(len(reference['sequence']), expected['length'])
                digest = hashlib.sha256(reference['sequence'].encode()).hexdigest()
                self.assertEqual(digest, expected['sha256'])
                self.assertEqual(reference['sha256'], digest)
                self.assertEqual(reference['mapping'], expected['mapping'])

    def test_layout_reproduces_exact_errfi1_window_and_dp1_assignment(self):
        for field in ('bin_width', 'tail_policy', 'bounds', 'allowed_intervals', 'transcription'):
            with self.subTest(field=field):
                self.assertEqual(self.saved[field], self.baseline['layout_exact'][field])
        bins, tails = make_bins(self.references, self.saved['bin_width'],
                               self.saved['tail_policy'], self.saved['assignments'],
                               self.saved['bounds'])
        self.assertEqual(len(bins), 1)
        self.assertEqual(tails, [])
        region = bins[0]
        self.assertEqual((region['start'], region['end'], region['length']),
                         (189019, 209018, 20000))
        self.assertEqual((region['chromosome'], region['genomic_start'], region['genomic_end']),
                         ('chr1', 8009019, 8029018))
        self.assertEqual(region['imager'], 'DP1')

    def test_dp1_imager_and_docking_polarity_are_verified(self):
        dp1 = next(i for i in validate_imagers(self.saved['imagers']) if i['id'] == 'DP1')
        expected = self.baseline['dp1_inventory_exact']
        self.assertEqual(dp1['imager_sequence_5to3'], expected['imager_sequence_5to3'])
        self.assertEqual(dp1['docking_sequence_5to3'], expected['docking_sequence_5to3'])
        self.assertEqual(dp1['imager_sequence_5to3'], 'CTAACGAATCGG')
        self.assertEqual(dp1['docking_sequence_5to3'], 'CCGATTCGTTAG')
        self.assertEqual(rc(dp1['imager_sequence_5to3']), dp1['docking_sequence_5to3'])
        self.assertEqual(dp1['input_kind'], 'imager')

    def test_saved_baseline_distinguishes_selected_tables_from_export_completion(self):
        observations = self.baseline['prior_result_observations']
        self.assertEqual(observations['selected_primaries_count'], 41)
        self.assertEqual(observations['selected_readout_counts'], {'DP1': 41})
        self.assertEqual(observations['strand_consistency']['top_sites'], 41)
        self.assertEqual(observations['strand_consistency']['bottom_sites'], 0)
        self.assertEqual(observations['strand_consistency']['opposite_strand_sites'], 0)
        self.assertFalse(observations['final_export_files_presence']['design_result.json'])

    def test_ui_bootstrap_uses_verified_reference_and_same_defaults(self):
        import studio
        boot = studio.bootstrap()
        self.assertEqual(boot['default_config'], self.saved)
        self.assertEqual(boot['default_config']['policy'], self.expected_policy)
        self.assertEqual(boot['references'], self.references)
        self.assertEqual(boot['version'], '1.6.1')

    def test_portable_settings_export_import_restores_actual_errfi1_inputs(self):
        import studio
        config = copy.deepcopy(self.saved)
        config['references'] = copy.deepcopy(self.references)
        config['background_text'] = '>construct_background\nACGTCGATCGATGCATGCAT\n'
        config['transcript_text'] = '>test_transcript\nACGUACGUACGUACGUACGU\n'
        expected = studio.normalize_config_inputs(config)
        expected['imagers'] = validate_imagers(expected['imagers'])
        original = copy.deepcopy(config)
        archive = studio.export_settings(config)
        restored = studio.import_settings(archive)
        self.assertEqual(restored, expected)
        self.assertEqual(config, original, 'Export must not mutate active settings.')
        self.assertEqual(restored['transcripts'][0]['sequence'], 'ACGTACGTACGTACGTACGT')
        self.assertEqual(restored['references'][0]['sequence'], self.references[0]['sequence'])
        with zipfile.ZipFile(io.BytesIO(archive)) as bundle:
            manifest = json.loads(bundle.read('run_config.json'))
            sequence_pointer = manifest['configuration']['references'][0]['sequence']
            self.assertIsInstance(sequence_pointer, dict)
            self.assertEqual(sequence_pointer['sha256'], self.references[0]['sha256'])
            self.assertIn(sequence_pointer['asset'], bundle.namelist())
            self.assertIn('inputs/target_region.fa.gz', bundle.namelist())
            self.assertLess(len(bundle.read('run_config.json')), 100000)


class ScientificImplementationPreservationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.manifest = read_json('verification/source_manifest.json')

    def test_specificity_folding_exact_matching_and_spacing_definitions_unchanged(self):
        source = (ROOT / 'design_engine.py').read_text()
        definitions = {n.name: n for n in ast.parse(source).body
                       if isinstance(n, (ast.FunctionDef, ast.ClassDef))}
        # v1.6.1 intentionally changes folding dispatch/cache management and
        # adds disk-backed specificity aggregation; dedicated equivalence tests
        # cover these. Retain the original hashes for all other scientific code.
        changed_for_v161 = {'Folding', 'summarize_specificity'}
        for name, expected in self.manifest['design_engine_definition_sha256'].items():
            if name in changed_for_v161:continue
            with self.subTest(definition=name):
                self.assertIn(name, definitions)
                text = ast.get_source_segment(source, definitions[name])
                # November 3 bounded memoization moved sequence-only work into
                # helpers. Verify their exact identities before normalizing the
                # extraction, then retain the original scientific baseline hash.
                # test_alignment_performance independently checks output parity,
                # thresholds, errors, mutation isolation and cache bounds.
                if name == 'perfect_tracts':
                    self.assertEqual(hashlib.sha256(text.encode()).hexdigest(),
                                     '8a19ace5e1a17711818f49a6a656126fa6bf77604d97374b791a9defa69aad10')
                    helper = ast.get_source_segment(source, definitions['_perfect_tract_metrics'])
                    self.assertEqual(hashlib.sha256(helper.encode()).hexdigest(),
                                     '5dc09625f45bf0230840e9118d4ab0aa084f07357b08fe14ae6db1a0e29fab7f')
                    start = '    runs=[];start=None'
                    returned = '    return max(map(len,runs),default=0),tuple(runs)'
                    self.assertEqual(helper.count(start), 1)
                    self.assertEqual(helper.count(returned), 1)
                    text = 'def perfect_tracts(q:str,s:str):\n' + helper[helper.index(start):].replace(returned, '    return runs')
                if name == 'enrich_hit':
                    helper = ast.get_source_segment(source, definitions['_exact_tract_tm'])
                    self.assertEqual(hashlib.sha256(helper.encode()).hexdigest(),
                                     'c6f2f93085cf11878646bb11aa734988c9104babcf48476122bcec8a1ca574c0')
                    cached = ("    q=h['qseq'].upper();s=h['sseq'].upper()\n"
                              "    h['exact_tract']=_perfect_tract_metrics(q,s)[0]\n"
                              "    h['offtarget_tm'],h['tm_model']=duplex_estimate(q,s,p.profile,h.get('rna',False))\n")
                    original = ("    runs=perfect_tracts(h['qseq'].upper(),h['sseq'].upper())\n"
                                "    h['exact_tract']=max(map(len,runs),default=0)\n"
                                "    h['offtarget_tm'],h['tm_model']=duplex_estimate(h['qseq'].upper(),h['sseq'].upper(),p.profile,h.get('rna',False))\n")
                    self.assertEqual(text.count(cached), 1)
                    text = text.replace(cached, original)
                    cached = ("    # Preserve assignment/error order and keep all candidate-specific decisions\n"
                              "    # below outside these bounded caches.\n"
                              "    h['exact_tract_tm']=_exact_tract_tm(q,s,p.profile) if not h.get('rna') else None\n")
                    original = ("    exact_tms=[]\n"
                                "    if not h.get('rna'):\n"
                                "        exact_tms=[duplex_estimate(s,s,p.profile,False)[0] for s in runs if len(s)>=8]\n"
                                "    h['exact_tract_tm']=max((v for v in exact_tms if v is not None),default=None)\n")
                    self.assertEqual(text.count(cached), 1)
                    text = text.replace(cached, original)
                if name == 'local_qc':
                    # Resource checks and the post-candidate score-storage
                    # hook are operational additions. Removing exactly those
                    # hooks must leave the original local specificity algorithm
                    # byte-for-byte unchanged. Score/raw parity, source unions,
                    # resume and disk bounds are tested in test_score_checkpoints.
                    text = text.replace('allow,p,log,resource_guard=None):', 'allow,p,log):')
                    text = text.replace("        if resource_guard:resource_guard.checkpoint('Local specificity')\n", '')
                    score_hook=("        if isinstance(c['offtarget_hits'],ScoreHits):\n"
                                "            c['offtarget_hits'].store.stage_scores('local',[c['id']],p)\n")
                    self.assertEqual(text.count(score_hook),1)
                    text = text.replace(score_hook,'')
                self.assertEqual(hashlib.sha256(text.encode()).hexdigest(), expected,
                                 name + ' differs from the verified v1.5 scientific implementation')

    def test_thermodynamics_rna_and_strand_modules_unchanged(self):
        for name, expected in self.manifest['unchanged_file_sha256'].items():
            with self.subTest(module=name):
                self.assertEqual(hashlib.sha256((ROOT / name).read_bytes()).hexdigest(), expected)


if __name__ == '__main__':
    unittest.main()
