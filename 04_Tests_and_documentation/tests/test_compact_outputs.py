"""Output regression fixtures are not a new genome-screened ERRFI1 design."""
import copy
import gzip
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from compact_outputs import (write_input_manifest, load_input_manifest, compact_result,
                             finalize_outputs, pair_flags, digest, file_digest)


def fixture_config(detail='standard'):
    return dict(title='compact fixture', output_detail=detail, bin_width=200,
                policy=dict(maximum_pair_stem=11, minimum_cofold_gain=-12.0,
                            minimum_arm_unpaired=.4, bound_strand='+'),
                references=[dict(id='X', sequence='ACGT' * 100, length=400, mapping=None,
                                 description='context', features=[dict(label='unique-feature-label', start=1, end=400, strand='+')],
                                 annotated_rna_sequences=[dict(id='spliced', sequence='GATTACA' * 6)])],
                background=[dict(id='vector', sequence='ACCTGG' * 4, features=[])],
                transcripts=[dict(id='RNA1', sequence='GGAGAA' * 5)],
                imagers=[dict(id='DP1', sequence='CTAGATGTAT')],
                bounds={'X': [101, 300]}, assignments={'X_B001': 'DP1'},
                databases={'human': '/external/database/human'})


def fixture_bins():
    return [dict(id='X_B001', reference='X', start=101, end=300, length=200)]


def fixture_result(config):
    pair = dict(a='C1', b='C1', max_exact_stem=10, mfe_gain_vs_separate=-11.75,
                intermolecular_pairs=4, blocked=False, status='RUN_AQUEOUS_MODEL')
    good = dict(id='C1', name='P1', bin_id='X_B001', status='selected', arm='GATTACA' * 6,
                full='TT' + 'GATTACA' * 6 + 'TT', reject_reasons=[], sites=[],
                specificity=dict(high_tm_sites=0, moderate_tm_sites=1),
                fold=dict(mean_arm_unpaired=.44), self_pair=pair,
                offtarget_hits=[dict(qseq='GATTACA', sseq='GATTACA', allowed=False)])
    bad = dict(id='C2', bin_id='X_B001', status='QC_rejected', reject_reasons=['high off-target'],
               offtarget_hits=[dict(qseq='GGG', sseq='GGG')])
    other = dict(a='C1', b='C1', max_exact_stem=2, mfe_gain_vs_separate=-2,
                 intermolecular_pairs=1, blocked=False, status='RUN_AQUEOUS_MODEL')
    return dict(version='1.6.0', policy=config['policy'], references=config['references'],
                imagers=config['imagers'], bins=fixture_bins(), candidates=[good, bad],
                selected_ids=['C1'], counts=dict(candidates=2, selected=1),
                selected_pair_interactions=[pair, other], blast_runs=[], warnings=['fixture only'])


class ManifestTests(unittest.TestCase):
    def test_all_inputs_restored_including_outside_context_and_spliced_rna(self):
        config = fixture_config()
        before = copy.deepcopy(config)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = write_input_manifest(config, fixture_bins(), root, '1.6.0')
            self.assertEqual(load_input_manifest(root / 'run_config.json'), config)
            self.assertEqual(config, before)
            self.assertTrue((root / 'inputs/reference_context.fa.gz').exists())
            self.assertEqual(manifest['target_regions'][0]['start'], 101)
            manifest_text = (root / 'run_config.json').read_text()
            self.assertNotIn(config['references'][0]['sequence'], manifest_text)
            self.assertNotIn('unique-feature-label', manifest_text)
            self.assertEqual(manifest['application_version'], '1.6.0')
            self.assertEqual(len(list(root.rglob('*annotations*'))), 1)

    def test_whole_context_target_reuses_one_sequence_asset(self):
        config = fixture_config()
        bins = [dict(id='X_B001', reference='X', start=1, end=400)]
        with tempfile.TemporaryDirectory() as tmp:
            manifest = write_input_manifest(config, bins, Path(tmp))
            self.assertNotIn('inputs/reference_context.fa.gz', manifest['assets'])
            self.assertEqual(load_input_manifest(Path(tmp) / 'run_config.json'), config)

    def test_damaged_asset_rejected_before_import(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            write_input_manifest(fixture_config(), fixture_bins(), root)
            with (root / 'inputs/reference_context.fa.gz').open('ab') as stream:
                stream.write(b'changed')
            with self.assertRaisesRegex(ValueError, 'SHA-256'):
                load_input_manifest(root / 'run_config.json')

    def test_changed_configuration_digest_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = write_input_manifest(fixture_config(), fixture_bins(), root)
            manifest['configuration']['policy']['maximum_pair_stem'] = 12
            (root / 'run_config.json').write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError, 'SHA-256'):
                load_input_manifest(root / 'run_config.json')

    def test_malicious_paths_rejected(self):
        for bad in ('../outside.gz', '/tmp/outside.gz', 'inputs/../../outside.gz', 'C:\\outside.gz'):
            with self.subTest(path=bad), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                manifest = write_input_manifest(fixture_config(), fixture_bins(), root)
                manifest['assets'][bad] = dict(sha256='0' * 64, format='fasta.gz')
                (root / 'run_config.json').write_text(json.dumps(manifest))
                with self.assertRaisesRegex(ValueError, 'safe relative path'):
                    load_input_manifest(root / 'run_config.json')

    def test_symlink_escape_rejected(self):
        with tempfile.TemporaryDirectory() as tmp, tempfile.TemporaryDirectory() as other:
            root = Path(tmp)
            manifest = write_input_manifest(fixture_config(), fixture_bins(), root)
            external = Path(other) / 'out.gz'
            external.write_bytes(b'outside')
            (root / 'escape.gz').symlink_to(external)
            manifest['assets']['escape.gz'] = dict(sha256=file_digest(external), format='fasta.gz')
            (root / 'run_config.json').write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError, 'escapes'):
                load_input_manifest(root / 'run_config.json')

    def test_target_slice_content_verified_against_source(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = write_input_manifest(fixture_config(), fixture_bins(), root)
            name = 'inputs/target_region.fa.gz'
            with gzip.open(root / name, 'wt') as stream:
                stream.write('>region_1\n' + 'A' * 200 + '\n')
            manifest['assets'][name]['sha256'] = file_digest(root / name)
            (root / 'run_config.json').write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError, 'Target slice'):
                load_input_manifest(root / 'run_config.json')

    def test_aggregate_nested_gzip_budget_rejects_combined_expansion(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = write_input_manifest(fixture_config(), fixture_bins(), root)
            sizes = [len(gzip.decompress((root / name).read_bytes())) for name in manifest['assets']]
            # Each individual asset fits; their combined decoded contents do not.
            budget = max(sizes)
            self.assertGreater(sum(sizes), budget)
            with self.assertRaisesRegex(ValueError, 'decompressed'):
                load_input_manifest(root / 'run_config.json', aggregate_decompressed_limit=budget)
            self.assertEqual(load_input_manifest(root / 'run_config.json', aggregate_decompressed_limit=sum(sizes)), fixture_config())


class CompactResultTests(unittest.TestCase):
    def test_selected_qc_and_original_counts_retained_without_raw_hits(self):
        config = fixture_config()
        result = fixture_result(config)
        before = copy.deepcopy(result)
        with tempfile.TemporaryDirectory() as tmp:
            manifest = write_input_manifest(config, fixture_bins(), Path(tmp))
            compact = compact_result(result, manifest)
        self.assertEqual(result, before)
        self.assertEqual([c['id'] for c in compact['candidates']], ['C1'])
        self.assertEqual(compact['counts'], result['counts'])
        self.assertEqual(compact['candidates'][0]['specificity'], result['candidates'][0]['specificity'])
        self.assertEqual(compact['candidates'][0]['offtarget_alignment_count'], 1)
        self.assertEqual(compact['imagers'], result['imagers'])
        text = json.dumps(compact)
        for forbidden in ('offtarget_hits', 'qseq', 'sseq', 'unique-feature-label'):
            self.assertNotIn(forbidden, text)
        self.assertNotIn(config['references'][0]['sequence'], text)
        self.assertEqual(compact['candidate_summary']['rejection_reasons'], {'high off-target': 1})
        self.assertEqual(compact['selected_pair_summary']['evaluated'], 2)
        self.assertEqual(len(compact['selected_pair_interactions']), 1)

    def test_failures_and_near_threshold_flags_are_only_reporting(self):
        policy = fixture_config()['policy']
        bad = dict(blocked=True, max_exact_stem=12, mfe_gain_vs_separate=-12.5, intermolecular_pairs=1)
        self.assertIn('FAIL_COMPATIBILITY', pair_flags(bad, policy))
        self.assertIn('NEAR_COFOLD_LIMIT', pair_flags(bad, policy))
        self.assertEqual(pair_flags(dict(blocked=False, max_exact_stem=9, mfe_gain_vs_separate=-10.9, intermolecular_pairs=1), policy), [])
        self.assertEqual(policy['minimum_cofold_gain'], -12.0)

    def test_standard_discards_raw_but_full_audit_preserves_compressed_once(self):
        for detail in ('standard', 'full_audit'):
            with self.subTest(detail=detail), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                config = fixture_config(detail)
                manifest = write_input_manifest(config, fixture_bins(), root)
                result = fixture_result(config)
                raw = b'fake raw alignment\n' * 100
                (root / 'human_candidates_blast.tsv').write_bytes(raw)
                (root / 'candidate_queries.fa').write_text('>query\nACGT\n')
                result['blast_runs'] = [dict(raw_file='human_candidates_blast.tsv', output_sha256=digest(raw))]
                compact = finalize_outputs(result, manifest, root)
                self.assertFalse((root / 'human_candidates_blast.tsv').exists())
                self.assertFalse((root / 'candidate_queries.fa').exists())
                self.assertEqual(json.loads((root / 'design_result.json').read_text()), compact)
                if detail == 'full_audit':
                    with gzip.open(root / 'debug/human_candidates_blast.tsv.gz', 'rb') as stream:
                        self.assertEqual(stream.read(), raw)
                    with gzip.open(root / 'debug/all_candidates.jsonl.gz', 'rt') as stream:
                        rows = [json.loads(line) for line in stream]
                    self.assertEqual([row['id'] for row in rows], ['C1', 'C2'])
                    self.assertNotIn('offtarget_hits', rows[0])
                    with gzip.open(root / 'debug/candidate_review_alignments.jsonl.gz', 'rt') as stream:
                        hits = [json.loads(line) for line in stream]
                    self.assertEqual(len(hits), 2)
                    self.assertEqual(hits[0]['hit'], result['candidates'][0]['offtarget_hits'][0])
                else:
                    self.assertFalse((root / 'debug').exists())
                    self.assertIsNone(compact['blast_runs'][0]['raw_file'])


if __name__ == '__main__':
    unittest.main()
