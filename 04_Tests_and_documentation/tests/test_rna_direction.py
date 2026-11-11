"""Tests of source-labeled presets, RNA direction, and selection (synthetic data only)."""
import copy
import io
import json
import random
import tempfile
import unittest
from pathlib import Path
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio import SeqIO
from design_engine import Policy, generate_candidates, design, local_qc, summarize_specificity, enrich_hit
from design_presets import PRESETS, LEGACY_PRESETS, preset_provenance
from sequence_tools import parse_fasta, parse_annotated, make_bins, validate_imagers, rc
from rna_avoidance import prepare_transcription, permitted_strands, feature_intervals, candidate_rna_safety
import annotation_sources as ann

BASE=Path(__file__).resolve().parents[1]
IMAGERS=validate_imagers(json.loads((BASE/'imagers.json').read_text()))

def ref(seed=113,n=600):
    rand=random.Random(seed)
    seq=''.join(rand.choices('ACGT',weights=[2,3,3,2],k=n))
    return dict(id='Test',sequence=seq,length=n,mapping=None,features=[])

def declaration(mode,**kw):return {'Test':dict(mode=mode,confirmed=True,**kw)}
def genes(r,features):
    r=copy.deepcopy(r);r['features']=[dict(type='gene',label=name,start=a,end=b,strand=s,feature_id=name,source='synthetic fixture') for name,a,b,s in features];return r

def candidates(r,mode='forward',forced=None,**kw):
    if forced is None:forced='+' if mode=='reverse' else '-'
    p=Policy(mode='local_preview',bound_strand=forced,shortlist_per_bin=70,avoid_gene_rna=True,**kw)
    model,tx,summary=prepare_transcription([r],declaration(mode))
    bins,_=make_bins([r],r['length'],assignments={'Test_B001':'DP7'})
    cs,count=generate_candidates([r],bins,IMAGERS,p,lambda *a,**k:None,model)
    return cs,count,model,tx,summary,bins,p

class PresetTests(unittest.TestCase):
    def test_default_preferred_and_extended_lengths(self):self.assertEqual((Policy().arm_min,Policy().preferred_arm_max,Policy().arm_max),(36,42,90))
    def test_default_tm_window(self):self.assertEqual((Policy().tm_min,Policy().tm_max),(45,55))
    def test_extended_length_is_explicit_alternative(self):self.assertEqual(LEGACY_PRESETS['expanded_g5']['values']['arm_max'],90)
    def test_no_false_imager_GC_rule(self):self.assertEqual((Policy().gc_min,Policy().gc_max),(0,76))
    def test_user_spacing_retained(self):self.assertEqual((Policy().minimum_gap,Policy().footprint),(5,'extended'))
    def test_default_requires_explicit_global_direction(self):self.assertEqual(Policy().bound_strand,'+');self.assertFalse(Policy().avoid_gene_rna);self.assertTrue(Policy().screen_rna)
    def test_preset_changes_recorded(self):self.assertIn('arm_max',preset_provenance({'preset_id':'ellenberg_informed','arm_max':90})['modified_fields'])
    def test_reject_nonboolean_avoidance(self):
        with self.assertRaises(ValueError):Policy.from_dict({'avoid_gene_rna':'yes'})
    def test_cannot_omit_RNA_search_under_avoidance(self):
        with self.assertRaises(ValueError):Policy.from_dict({'bound_strand':'+','avoid_gene_rna':True,'screen_rna':False})
    def test_explicit_opt_out_is_possible_and_labeled(self):
        p=Policy.from_dict({'bound_strand':'-','avoid_gene_rna':False,'screen_rna':False});self.assertFalse(p.avoid_gene_rna)
    def test_default_preset_matches_backend(self):
        p=Policy()
        for k,v in PRESETS['errfi1']['values'].items():self.assertEqual(getattr(p,k),v,k)

class DirectionTests(unittest.TestCase):
    def test_plus_gene_requires_bottom_DNA(self):self.assertEqual(permitted_strands([dict(start=1,end=100,strand='+')],20,60),['-'])
    def test_minus_gene_requires_top_DNA(self):self.assertEqual(permitted_strands([dict(start=1,end=100,strand='-')],20,60),['+'])
    def test_overlap_opposite_genes_excludes_both(self):self.assertEqual(permitted_strands([dict(start=1,end=80,strand='+'),dict(start=50,end=120,strand='-')],60,90),[])
    def test_unknown_gene_strand_excludes_both(self):self.assertEqual(permitted_strands([dict(start=1,end=100,strand='?')],20,60),[])
    def test_any_partial_gene_overlap_is_guarded(self):self.assertEqual(permitted_strands([dict(start=60,end=100,strand='+')],20,60,'+'),[])
    def test_intergenic_convention_is_bottom(self):self.assertEqual(permitted_strands([dict(start=1,end=100,strand='-')],101,160),['-'])
    def test_forced_wrong_strand_is_not_silently_flipped(self):self.assertEqual(permitted_strands([dict(start=1,end=100,strand='+')],20,60,'+'),[])
    def test_raw_sequence_requires_direction(self):
        with self.assertRaisesRegex(ValueError,'not been confirmed'):prepare_transcription([ref()],{})
    def test_construct_requires_confirmation(self):
        with self.assertRaisesRegex(ValueError,'confirm'):prepare_transcription([ref()],{'Test':{'mode':'forward'}})
    def test_fetched_empty_gene_list_not_missing(self):
        r=ref();r['mapping']={'verified':True};r['gene_annotation']={'status':'FETCHED'}
        model,tx,s=prepare_transcription([r],{});self.assertEqual(model['Test'],[]);self.assertEqual(tx,[])
    def test_missing_genomic_annotations_block(self):
        r=ref();r['mapping']={'verified':True}
        with self.assertRaisesRegex(ValueError,'fetch gene'):prepare_transcription([r],{})
    def test_H3K27ac_does_not_infer_transcription(self):
        r=ref();r['features']=[dict(type='H3K27ac',start=1,end=600,strand='+')]
        self.assertEqual(feature_intervals(r),[])
        with self.assertRaises(ValueError):prepare_transcription([r],{})
    def test_no_known_transcription_does_not_erase_gene(self):
        with self.assertRaisesRegex(ValueError,'cannot be erased'):prepare_transcription([genes(ref(),[('X',1,600,'+')])],declaration('none'))
    def test_manual_coordinates_are_checked(self):
        with self.assertRaises(ValueError):prepare_transcription([ref()],declaration('manual',intervals=[dict(start=0,end=500,strand='+')]))
    def test_manual_float_rejected(self):
        with self.assertRaises(ValueError):prepare_transcription([ref()],declaration('manual',intervals=[dict(start=1.1,end=500,strand='+')]))
    def test_no_manual_rows_rejected(self):
        with self.assertRaises(ValueError):prepare_transcription([ref()],declaration('manual',intervals=[]))
    def test_both_declared_conservatively_excludes_all(self):
        r=ref();model,_,_=prepare_transcription([r],declaration('both'));self.assertEqual(permitted_strands(model['Test'],10,70),[])
    def test_plus_RNA_sequence_is_reference_sense(self):
        r=ref();_,tx,_=prepare_transcription([r],declaration('forward'));self.assertEqual(tx[0]['sequence'],r['sequence'])
    def test_minus_RNA_sequence_is_reverse_complement(self):
        r=ref();_,tx,_=prepare_transcription([r],declaration('reverse'));self.assertEqual(tx[0]['sequence'],rc(r['sequence']))
    def test_repeat_copy_wrong_direction_rejected(self):
        sites=[dict(reference='Test',start=10,end=49,strand='-'),dict(reference='Test',start=200,end=239,strand='-')]
        result=candidate_rna_safety(sites,{'Test':[dict(start=1,end=80,strand='+',name='plus'),dict(start=190,end=300,strand='-',name='minus')]})
        self.assertEqual(result['status'],'RNA_DIRECTION_CONFLICT');self.assertEqual(len(result['conflicts']),1)
    def test_unknown_annotation_strand_preserved(self):
        rec=SeqRecord(Seq('A'*100),id='Test',annotations={'molecule_type':'DNA'});rec.features=[SeqFeature(FeatureLocation(10,50,strand=None),type='gene')]
        buf=io.StringIO();SeqIO.write(rec,buf,'genbank');r=parse_annotated(buf.getvalue().encode(),'test.gb')[0]
        # GenBank serializes an unqualified location as forward; parser fixture
        # below explicitly tests backend handling of a missing direction.
        f=dict(start=10,end=50,type='gene',label='unknown',strand='?');self.assertEqual(feature_intervals(dict(ref(),features=[f]))[0]['strand'],'?')
    def test_compound_mRNA_protects_introns_and_splice(self):
        r=ref(n=100);rec=SeqRecord(Seq(r['sequence']),id='Test',annotations={'molecule_type':'DNA'})
        rec.features=[SeqFeature(CompoundLocation([FeatureLocation(10,30,strand=-1),FeatureLocation(50,75,strand=-1)]),type='mRNA',qualifiers={'gene':['TestGene']})]
        buf=io.StringIO();SeqIO.write(rec,buf,'genbank');rr=parse_annotated(buf.getvalue().encode(),'test.gb')[0]
        spans=feature_intervals(rr);self.assertEqual((spans[0]['start'],spans[0]['end']),(11,75));self.assertEqual(len(rr['annotated_rna_sequences'][0]['sequence']),45)
        _,tx,_=prepare_transcription([rr],declaration('annotations'));self.assertEqual(len(tx),2)
    def test_genomic_reverse_view_flips_gene_direction(self):
        r=ref(n=100);r['mapping']=dict(verified=True,start=1001,end=1100,chromosome='chr7',orientation='-')
        d=ann.fetch_genes(r,lambda url:[dict(start=1010,end=1080,strand=1,feature_type='gene',id='X',assembly_name='GRCh38')]);self.assertEqual(d['features'][0]['strand'],'-')
    def test_gene_annotation_request_includes_transcripts(self):
        r=ref(n=100);r['mapping']=dict(verified=True,start=1001,end=1100,chromosome='chr7')
        seen=[];ann.fetch_genes(r,lambda url:(seen.append(url) or []));self.assertIn('feature=transcript',seen[0])

class SelectionTests(unittest.TestCase):
    def test_forward_candidates_are_bottom_and_sense(self):
        r=ref();cs,*_=candidates(r);self.assertTrue(cs)
        for c in cs:
            self.assertEqual(c['bound_strand'],'-');a=c['generation_start'];e=c['generation_end'];self.assertEqual(c['arm'],r['sequence'][a-1:e]);self.assertLessEqual(len(c['arm']),Policy().arm_max)
    def test_reverse_candidates_are_top_and_RNA_sense(self):
        r=ref();cs,*_=candidates(r,'reverse');self.assertTrue(cs)
        for c in cs:
            self.assertEqual(c['bound_strand'],'+');a=c['generation_start'];e=c['generation_end'];self.assertEqual(c['arm'],rc(r['sequence'][a-1:e]))
    def test_wrong_forced_strand_yields_no_candidates(self):
        cs,count,*_=candidates(ref(),'forward','+');self.assertEqual(cs,[]);self.assertGreater(count[0]['RNA_direction_conflict'],0)
    def test_unknown_annotation_blocks_only_overlap(self):
        r=genes(ref(),[('unknown',200,400,'?')]);cs,count,*_=candidates(r,'annotations');self.assertTrue(cs);self.assertTrue(all(c['generation_end']<200 or c['generation_start']>400 for c in cs))
    def test_opposing_gene_excludes_candidates_without_switching_strand(self):
        r=genes(ref(),[('p',1,240,'+'),('m',301,600,'-')]);cs,*_=candidates(r,'annotations');self.assertEqual({c['bound_strand'] for c in cs},{'-'});self.assertTrue(all(c['generation_end']<301 for c in cs))
    def test_RNA18_rejected_even_when_human_DNA_is_ignored(self):
        p=Policy(mode='local_preview',human_policy='report_only',avoid_gene_rna=True);c={'tm':47,'offtarget_hits':[],'reject_reasons':[]}
        seq='ACGTTAGCCTAGGTCGAC';h=dict(source='rna_human',reference='tx',start=1,end=18,bound_strand='RNA',qseq=seq,sseq=seq,identity=100,query_coverage=45,alignment_length=18,allowed=False,rna=True)
        c['offtarget_hits']=[enrich_hit(h,c,p)];summarize_specificity(c,p);self.assertEqual(c['specificity']['rna_risk_sites'],1);self.assertTrue(c['reject_reasons'])
    def test_same_sense_not_counted_as_RNA_complement(self):
        r=ref();cs,_,model,tx,_,bins,p=candidates(r);local_qc(cs[:4],[r],bins,[],tx,[],p,lambda *a,**k:None)
        self.assertFalse(any(h.get('rna') and h['query_coverage']==100 and h['identity']==100 for c in cs[:4] for h in c['offtarget_hits']))
    def test_wrong_RNA_source_exact_arm_is_rejected(self):
        r=ref();cs,_,model,tx,_,bins,p=candidates(r);c=cs[0];tx=[dict(id='bad_rna',sequence=rc(c['arm']))]
        local_qc([c],[r],bins,[],tx,[],p,lambda *a,**k:None);summarize_specificity(c,p);self.assertTrue(c['reject_reasons'])
    def test_missing_transcription_stops_before_output(self):
        r=ref();cfg=dict(references=[r],imagers=IMAGERS,assignments={'Test_B001':'DP7'},bin_width=600,policy={'bound_strand':'-','mode':'local_preview','avoid_gene_rna':True})
        with tempfile.TemporaryDirectory() as t:
            with self.assertRaisesRegex(ValueError,'gene direction'):design(cfg,Path(t)/'out',lambda *a,**k:None)
            self.assertFalse((Path(t)/'out').exists())
    def test_real_de_novo_selection_records_new_rules(self):
        r=genes(ref(),[('p',1,240,'+'),('m',301,600,'-')]);cfg=dict(references=[r],imagers=IMAGERS,assignments={'Test_B001':'DP7'},bin_width=600,transcription=declaration('annotations'),policy={'bound_strand':'-','mode':'local_preview','avoid_gene_rna':True,'shortlist_per_bin':30,'target_per_bin':6,'minimum_per_bin':2})
        with tempfile.TemporaryDirectory() as t:
            result=design(cfg,Path(t)/'out',lambda *a,**k:None)
            self.assertTrue(result['selected_ids']);self.assertEqual(result['version'],'1.6.1');self.assertTrue((Path(t)/'out'/'transcription_model.csv').exists())
            self.assertTrue(result['transcription_summary']);self.assertTrue(result['parameter_provenance'])
            for c in result['candidates']:
                if c['id'] in result['selected_ids']:self.assertEqual(c['rna_avoidance']['status'],'DIRECTION_COMPATIBLE');self.assertEqual(c['specificity']['rna_risk_sites'],0)

if __name__=='__main__':unittest.main()
