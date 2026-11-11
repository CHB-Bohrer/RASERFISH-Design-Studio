"""Consolidated edition regression checks. Synthetic data only, not FISH QC."""
import copy
import importlib.util
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch
from design_engine import Policy, enrich_hit, summarize_specificity, generate_candidates, design
from design_presets import PRESETS
from rna_avoidance import prepare_transcription
from sequence_tools import make_bins, rc
from test_rna_direction import ref, genes, declaration, IMAGERS

BASE=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('setup_helper', BASE/'setup/studio_install_helper.py')
helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper)


def scored_rna(length=18,coverage=45,identity=100,strict=False,source='rna_human'):
    seq=('ACGTTAGCCTAGGTCGACTGCAATGCTCGATGCATCGTACGATGC')[:length]
    p=Policy(bound_strand='+',avoid_gene_rna=strict,mode='local_preview')
    c={'tm':50,'offtarget_hits':[],'reject_reasons':[]}
    h=dict(source=source,reference='tx',start=1,end=length,bound_strand='RNA',qseq=seq,sseq=seq,
           identity=identity,query_coverage=coverage,alignment_length=length,allowed=False,rna=True)
    c['offtarget_hits']=[enrich_hit(h,c,p)];summarize_specificity(c,p)
    return c

class RNAPolicyTests(unittest.TestCase):
    def test_default_is_transcript_search_without_span_mask(self):
        p=Policy.from_dict({'bound_strand':'+'})
        self.assertFalse(p.avoid_gene_rna);self.assertTrue(p.screen_rna);self.assertEqual(p.max_strong_rna_sites,0)
    def test_no_automatic_strand_at_new_default(self):
        self.assertEqual(Policy.from_dict({}).bound_strand,'+')
    def test_presets_restore_errfi1_rna_policy(self):
        for p in PRESETS.values():
            for key,value in dict(avoid_gene_rna=False,screen_rna=True,rna_exact_min=18,max_strong_rna_sites=0,bound_strand='+').items():
                self.assertEqual(p['values'][key],value)
    def test_18nt_default_is_review_not_rejection(self):
        c=scored_rna();self.assertFalse(c['reject_reasons']);self.assertTrue(c['offtarget_hits'][0]['review'])
        self.assertTrue(c['offtarget_hits'][0]['rna_exact_risk']);self.assertEqual(c['specificity']['rna_risk_sites'],0)
    def test_18nt_strict_is_rejection(self):self.assertTrue(scored_rna(strict=True)['reject_reasons'])
    def test_30nt_default_still_rejected(self):self.assertTrue(scored_rna(length=30,coverage=60)['reject_reasons'])
    def test_near_full_default_still_rejected(self):self.assertTrue(scored_rna(length=29,coverage=92)['reject_reasons'])
    def test_supplied_rna_not_exempted(self):self.assertTrue(scored_rna(length=35,coverage=95,source='rna_local')['reject_reasons'])
    def test_strict_requires_public_rna_search(self):
        with self.assertRaises(ValueError):Policy.from_dict({'bound_strand':'+','avoid_gene_rna':True,'screen_rna':False})
    def test_relaxed_does_not_generate_unspliced_references(self):
        r=genes(ref(n=240),[('PLUS',1,240,'+')]);before=copy.deepcopy(r)
        model,tx,summary=prepare_transcription([r],{},False)
        self.assertFalse(model);self.assertFalse(tx);self.assertEqual(r,before)
        self.assertEqual(summary[0]['status'],'GENE_SPANS_NOT_USED_FOR_EXCLUSION')
    def test_opposing_annotation_does_not_flip_or_mask_relaxed_candidates(self):
        r=genes(ref(n=240),[('PLUS',1,240,'+')]);bins,_=make_bins([r],240,assignments={'Test_B001':'DP7'})
        p=Policy(bound_strand='+',mode='local_preview',shortlist_per_bin=10,extension_shortlist_per_bin=5)
        cs,_=generate_candidates([r],bins,IMAGERS,p,lambda *a,**k:None,{})
        self.assertTrue(cs);self.assertEqual({c['bound_strand'] for c in cs},{'+'})
        for c in cs:
            self.assertEqual(c['arm'],rc(r['sequence'][c['generation_start']-1:c['generation_end']]))
    def test_strict_still_excludes_conflicting_whole_span(self):
        r=genes(ref(n=240),[('PLUS',1,240,'+')]);bins,_=make_bins([r],240,assignments={'Test_B001':'DP7'})
        p=Policy(bound_strand='+',mode='local_preview',avoid_gene_rna=True,shortlist_per_bin=10,extension_shortlist_per_bin=5)
        model,_,_=prepare_transcription([r],declaration('annotations'),True)
        cs,stats=generate_candidates([r],bins,IMAGERS,p,lambda *a,**k:None,model)
        self.assertFalse(cs);self.assertGreater(stats[0]['RNA_direction_conflict'],0)
    def test_saved_tutorial_policy_and_one_strand(self):
        d=json.loads((BASE/'examples/local_design/design_result.json').read_text())
        self.assertEqual(d['version'],'1.5.0');self.assertEqual(d['rna_policy_summary']['mode'],'transcript_sequences_only')
        self.assertEqual(d['rna_policy_summary']['rnase_treatment'],'NOT_ASSUMED_OR_VERIFIED')
        self.assertEqual(d['strand_audit']['status'],'PASS_SEQUENCE_ORIENTATION')
    def test_full_mode_requires_rna_database_in_relaxed_mode(self):
        r=ref(n=100);cfg=dict(references=[r],imagers=IMAGERS,assignments={'Test_B001':'DP7'},bin_width=100,
            policy=dict(bound_strand='+',mode='full'),databases={'human':'nonexistent'})
        with tempfile.TemporaryDirectory() as t:
            with self.assertRaisesRegex(RuntimeError,'RNA screening'):design(cfg,Path(t)/'result',lambda *a,**k:None)

class SameFolderInstallTests(unittest.TestCase):
    def skeleton(self,root):
        for name in ('studio.py','design_engine.py','sequence_tools.py','engine/qc_core.py','start_studio.sh','setup/studio_install_helper.py'):
            f=root/name;f.parent.mkdir(parents=True,exist_ok=True);f.write_text('# fixture only\n')
        (root/'install.sh').write_text((BASE/'install.sh').read_text())
    def test_dry_run_uses_script_folder_not_current_directory(self):
        with tempfile.TemporaryDirectory() as t:
            base=Path(t);root=base/'fresh';root.mkdir();self.skeleton(root)
            p=subprocess.run(['bash',str(root/'install.sh'),'--dry-run'],cwd=base,text=True,capture_output=True)
            self.assertEqual(p.returncode,0,p.stderr);self.assertIn(str(root/'.qc_env'),p.stdout)
            self.assertFalse((root/'.tools').exists());self.assertFalse((base/'.qc_env').exists())
    def test_symbolic_reference_directory_rejected(self):
        with tempfile.TemporaryDirectory() as t:
            base=Path(t);root=base/'fresh';root.mkdir();self.skeleton(root);(root/'references').symlink_to(base)
            p=subprocess.run(['bash',str(root/'install.sh'),'--dry-run'],text=True,capture_output=True)
            self.assertNotEqual(p.returncode,0);self.assertIn('symbolic link',p.stderr)
    def test_unknown_installer_flag_fails(self):
        p=subprocess.run(['bash',str(BASE/'install.sh'),'--erase-old'],text=True,capture_output=True)
        self.assertNotEqual(p.returncode,0);self.assertIn('Unknown option',p.stderr)
    def test_launcher_does_not_require_old_folder(self):
        text=(BASE/'start_studio.sh').read_text();self.assertNotIn('Documents/Codex',text)
        self.assertNotIn('RASERFISH_QC_ROOT',text);self.assertIn('First-time setup: bash install.sh',text)
    def test_helper_safe_relative_path(self):
        with tempfile.TemporaryDirectory() as t:
            base=Path(t)
            for name in ('../outside','/absolute'):
                with self.assertRaises(ValueError):helper.within(base,name)
    def test_archive_traversal_rejected(self):
        import zipfile
        with tempfile.TemporaryDirectory() as t:
            p=Path(t);z=p/'bad.zip'
            with zipfile.ZipFile(z,'w') as f:f.writestr('../not_here','x')
            with self.assertRaises(ValueError):helper.extract_zip_safely(z,p/'out')
            self.assertFalse((p/'not_here').exists())
    def test_reference_missing_index_integrity_rejected(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t);f=p/'seq.fa';f.write_text('>fixture\nACGT\n')
            d=dict(fasta_relative='seq.fa',fasta_sha256=helper.sha256(f),database_prefix_relative='db')
            (p/'reference_manifest.json').write_text(json.dumps(dict(assembly_accession=helper.ACCESSION,databases={x:d for x in helper.LABELS})))
            with patch.object(helper,'db_info',return_value='fixture: not a real BLAST database'):
                with self.assertRaisesRegex(RuntimeError,'integrity records'):helper.verify_references(p)
    def test_no_probe_selection_in_installer(self):
        text=(BASE/'install.sh').read_text();self.assertNotIn('apply_single_strand_guard',text)
        self.assertNotIn('old/.qc_env',text);self.assertIn('REGRESSION_TESTS.txt',text)
    def test_bash_syntax(self):
        for f in ('install.sh','start_studio.sh','Install.command','Open_Studio.command'):
            p=subprocess.run(['bash','-n',str(BASE/f)],capture_output=True,text=True)
            self.assertEqual(p.returncode,0,p.stderr)

if __name__=='__main__':unittest.main()
