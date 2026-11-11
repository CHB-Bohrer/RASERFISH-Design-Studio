"""Fixture and optional real-engine checks; fixtures are not genomic results."""
import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from design_engine import *
from sequence_tools import validate_imagers
BASE=Path(__file__).resolve().parents[1]
IMS=validate_imagers(json.loads((BASE/'imagers.json').read_text()))

def inputs():
    arm='AGCTTGACGTCAGCTAGGCTCATGACCTAGTCGGATCGTA'
    pre=IMS[0]['docking']+'TT';suf='TT'+IMS[0]['docking']
    c=dict(id='fixture',arm=arm,full=pre+arm+suf,prefix=pre,suffix=suf,bin_id='T_B001',sites=[dict(strand='-')],
      tm=48.0,offtarget_hits=[],reject_reasons=[])
    ref=dict(id='T',sequence='A'*99+arm+'C'*100,length=199+len(arm),mapping=dict(verified=True,start=1,end=199+len(arm),chromosome='chr7',orientation='+',assembly='GRCh38'))
    bins=[dict(id='T_B001',reference='T',start=1,end=ref['length'])]
    return c,ref,bins

class BlastFixtureTests(unittest.TestCase):
    def execute(self,rna=False,wrong=False):
        c,ref,bins=inputs()
        def fake(command,log,cancel=None,resource_guard=None,*,env_overrides=None):
            self.assertEqual(env_overrides,{'BATCH_SIZE':str(len(c['arm'])+1)})
            qpath=Path(command[command.index('-query')+1]);records=qpath.read_text().splitlines()
            self.assertEqual(len(records),4)
            rows=[]
            for index in (0,2):
                ident=records[index][1:];seq=records[index+1]
                qstart=1 if index==0 else len(c['prefix'])+1
                q=c['arm'] if not wrong else 'C'*len(c['arm'])
                for rev in ([False,True] if rna else [False]):
                    a,e=100,99+len(c['arm']);a,e=(e,a) if rev else (a,e)
                    fields=[ident,'ref|NM_TEST.1|' if rna else 'ref|NC_000007.14|',100,len(q),0,0,qstart,qstart+len(q)-1,a,e,1e-10,80,len(seq),q,c['arm']]
                    rows.append('\t'.join(map(str,fields)))
            Path(command[command.index('-out')+1]).write_text('\n'.join(rows)+'\n')
            log.write_text('Explicit mock BLAST fixture\n')
        with tempfile.TemporaryDirectory() as temp,patch('design_engine.shutil.which',return_value='/fixture/blastn'),patch('design_engine.run_process',side_effect=fake),patch('design_engine.subprocess.check_output',return_value='fixture version'):
            runs=blast_candidates([c],[ref],bins,[],{'rna_human' if rna else 'human':'fixture'},Policy(mode='local_preview'),Path(temp),lambda *a,**k:None)
        return c,runs
    def test_searches_both_full_and_arm(self):
        c,r=self.execute();self.assertEqual({h['full_primary'] for h in c['offtarget_hits']},{True,False});self.assertEqual(len(c['offtarget_hits']),2)
    def test_expected_native_target_retained_as_allowed(self):
        c,r=self.execute();self.assertTrue(all(h['allowed'] for h in c['offtarget_hits']));summarize_specificity(c,Policy());self.assertEqual(c['specificity']['high_tm_sites'],0)
    def test_rna_only_complementary_orientation_kept(self):
        c,r=self.execute(True);self.assertEqual(len(c['offtarget_hits']),2);self.assertTrue(all(h['bound_strand']=='RNA' and not h['allowed'] for h in c['offtarget_hits']))
    def test_wrong_query_sequence_aborts(self):
        with self.assertRaises(ValueError):self.execute(wrong=True)
    def test_empty_search_does_not_claim_exhaustiveness(self):
        c,r=self.execute();self.assertIn('heuristic',r[0]['note'])
    def test_policy_rejects_fractional_integer_fields(self):
        with self.assertRaises(ValueError):Policy.from_dict({'minimum_gap':3.5})
    def test_policy_rejects_negative_site_limit(self):
        with self.assertRaises(ValueError):Policy.from_dict({'max_high_offtarget_sites':-1})
    def test_policy_rejects_false_string_boolean(self):
        with self.assertRaises(ValueError):Policy.from_dict({'tm_strict':'false'})

@unittest.skipUnless(shutil.which('blastn') and shutil.which('makeblastdb'),'Real BLAST+ not installed; explicit fixtures above are not engine validation')
class RealBlastTests(unittest.TestCase):
    def test_tiny_real_native_reference(self):
        c,ref,bins=inputs()
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);fa=p/'tiny.fa';fa.write_text('>NC_000007.14\n'+ref['sequence']+'\n')
            subprocess.run(['makeblastdb','-in',str(fa),'-dbtype','nucl','-parse_seqids','-blastdb_version','4','-out',str(p/'db')],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
            runs=blast_candidates([c],[ref],bins,[],{'human':str(p/'db')},Policy(mode='full',threads=1),p,lambda *a,**k:None)
            self.assertTrue(any(h['allowed'] and not h['full_primary'] and h['query_coverage']==100 and h['identity']==100 for h in c['offtarget_hits']))
            self.assertTrue(any(h['full_primary'] for h in c['offtarget_hits']))
            self.assertEqual(runs[0]['status'],'RUN_HEURISTIC_SEARCH')
