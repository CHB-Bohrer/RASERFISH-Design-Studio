"""Reporting verifies omitted cross-primary checks never become a PASS."""
import copy
import csv
import gzip
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'tests'))
from compact_outputs import compact_result, finalize_outputs, write_input_manifest
from test_compact_outputs import fixture_config, fixture_bins, fixture_result

DISABLED = 'NOT RUN - USER DISABLED'


class UnreadableLargePayload(list):
    def __iter__(self):
        raise AssertionError('Excluded heavy evidence must not be traversed')

    def __len__(self):
        raise AssertionError('Recorded alignment count must avoid touching heavy evidence')

    def __deepcopy__(self, memo):
        raise AssertionError('Excluded heavy evidence must not be copied')


class ReportingTests(unittest.TestCase):
    def make_result(self, enabled=False, detail='standard'):
        config = fixture_config(detail)
        config['policy']['screen_cross_primary'] = enabled
        result = fixture_result(config)
        result['version'] = '1.6.1'
        result['cross_primary_screening'] = dict(enabled=enabled, status='RUN_AQUEOUS_MODEL' if enabled else DISABLED,
                                                 within_region_status='RUN_AQUEOUS_MODEL' if enabled else DISABLED,
                                                 between_region_status='RUN_AQUEOUS_MODEL' if enabled else DISABLED,
                                                 evaluated=0)
        result['self_dimer_screening'] = dict(status='RUN_AQUEOUS_MODEL', evaluated=2)
        return config, result

    def compact(self, config, result):
        with tempfile.TemporaryDirectory() as tmp:
            manifest = write_input_manifest(config, fixture_bins(), tmp)
            return compact_result(result, manifest)

    def test_disabled_cross_primary_scope_explicit_with_retained_self_evidence(self):
        config, result = self.make_result()
        compact = self.compact(config, result)
        cross = compact['cross_primary_screening']
        self.assertEqual(cross['status'], DISABLED)
        self.assertEqual(cross['within_region_status'], DISABLED)
        self.assertEqual(cross['between_region_status'], DISABLED)
        self.assertEqual(cross['evaluated'], 0)
        self.assertIsNone(cross['max_exact_stem'])
        self.assertEqual(compact['order_status'], 'DRAFT')
        self.assertEqual(compact['self_dimer_screening']['status'], 'RUN_AQUEOUS_MODEL')
        self.assertEqual(compact['candidates'][0]['fold'], result['candidates'][0]['fold'])
        self.assertEqual(compact['candidates'][0]['self_pair'], result['candidates'][0]['self_pair'])
        self.assertEqual({p['scope'] for p in compact['selected_pair_interactions']}, {'SELF_DIMER'})
        self.assertNotIn('PASS', json.dumps(cross))

    def test_enabled_restored_and_regions_are_separately_identified(self):
        config, result = self.make_result(True)
        for cid, region in [('C3', 'X_B001'), ('C4', 'X_B002')]:
            c = copy.deepcopy(result['candidates'][0]); c.update(id=cid, bin_id=region)
            result['candidates'].append(c); result['selected_ids'].append(cid)
            pair = dict(result['selected_pair_interactions'][0], a='C1', b=cid)
            result['selected_pair_interactions'].append(pair)
        compact = self.compact(config, result)
        self.assertTrue(compact['cross_primary_screening']['enabled'])
        self.assertEqual(compact['cross_primary_screening']['selected_evaluated'], 2)
        self.assertEqual({p['scope'] for p in compact['selected_pair_interactions']}, {'SELF_DIMER', 'WITHIN_REGION', 'BETWEEN_REGION'})
        self.assertNotIn(DISABLED, json.dumps(compact['cross_primary_screening']))

    def test_historical_missing_setting_is_not_relabelled_disabled(self):
        config, result = self.make_result()
        del config['policy']['screen_cross_primary']
        del result['cross_primary_screening']
        compact = self.compact(config, result)
        self.assertTrue(compact['cross_primary_screening']['enabled'])
        self.assertEqual(compact['cross_primary_screening']['status'], 'RECORDED_IN_LEGACY_RESULT')
        self.assertNotIn(DISABLED, json.dumps(compact['cross_primary_screening']))

    def test_inconsistent_disabled_result_with_cross_rows_cannot_export(self):
        config, result = self.make_result()
        result['selected_pair_interactions'].append(dict(result['selected_pair_interactions'][0], a='C1', b='C2'))
        with self.assertRaisesRegex(ValueError, 'Disabled cross-primary'):
            self.compact(config, result)

    def test_heavy_alignment_lists_not_traversed_or_copied(self):
        config, result = self.make_result()
        result['candidates'][0]['offtarget_hits'] = UnreadableLargePayload()
        result['candidates'][0]['offtarget_alignment_count'] = 10000000
        compact = self.compact(config, result)
        self.assertEqual(compact['candidates'][0]['offtarget_alignment_count'], 10000000)
        self.assertNotIn('offtarget_hits', compact['candidates'][0])

    def test_standard_and_full_audit_export_disabled_scope_and_self_tag(self):
        for detail in ('standard', 'full_audit'):
            with self.subTest(detail=detail), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                config, result = self.make_result(detail=detail)
                manifest = write_input_manifest(config, fixture_bins(), root)
                finalize_outputs(result, manifest, root)
                self.assertEqual(json.loads((root/'cross_primary_screening.json').read_text())['status'], DISABLED)
                self.assertEqual(json.loads((root/'self_dimer_screening.json').read_text())['status'], 'RUN_AQUEOUS_MODEL')
                with (root/'pair_interaction_summary.csv').open() as stream:
                    rows = list(csv.DictReader(stream))
                self.assertEqual(rows[0]['cross_primary_status'], DISABLED)
                self.assertEqual(rows[0]['within_region_status'], DISABLED)
                self.assertEqual(rows[0]['between_region_status'], DISABLED)
                if detail == 'full_audit':
                    with gzip.open(root/'debug/all_selected_pair_interactions.csv.gz', 'rt') as stream:
                        self.assertEqual({r['scope'] for r in csv.DictReader(stream)}, {'SELF_DIMER'})

    def test_ui_uses_saved_scope_and_preserves_explicit_enabled_import(self):
        node = shutil.which('node')
        if not node or not Path(node).exists():
            self.skipTest('Node runtime unavailable')
        script = r"""
const fs=require('fs'), vm=require('vm'), assert=require('assert');
const src=fs.readFileSync('web/app.js','utf8');
const begin=src.indexOf('function crossPrimaryResultStatus('), end=src.indexOf('function renderPairs(',begin);
const context={};vm.createContext(context);vm.runInContext(src.slice(begin,end),context);
const no={policy:{screen_cross_primary:false},cross_primary_screening:{enabled:false,status:'PASS'}};
assert.equal(context.crossPrimaryResultStatus(no).status,'NOT RUN - USER DISABLED');
const text=context.crossPrimaryResultText(no);assert(text.includes('Within-region: NOT RUN - USER DISABLED'));assert(text.includes('Between-region: NOT RUN - USER DISABLED'));assert(text.includes('Order status: DRAFT'));assert(!text.includes('PASS'));
const old={policy:{mode:'full'}};assert.equal(context.crossPrimaryResultStatus(old).enabled,true);assert(!context.crossPrimaryResultText(old).includes('USER DISABLED'));
const yes={policy:{screen_cross_primary:true},cross_primary_screening:{enabled:true,status:'RUN_AQUEOUS_MODEL',within_region_status:'RUN_AQUEOUS_MODEL',between_region_status:'RUN_AQUEOUS_MODEL'}};assert(context.crossPrimaryResultText(yes).includes('screening: ON'));assert(!context.crossPrimaryResultText(yes).includes('USER DISABLED'));
assert.equal(context.resultPairScope({a:'A',b:'A'},{}),'SELF_DIMER');
const candidates=[{id:'A',bin_id:'R1'},{id:'B',bin_id:'R1'},{id:'C',bin_id:'R2'}];assert.equal(context.resultPairScope({a:'A',b:'B'},{candidates}),'WITHIN_REGION');assert.equal(context.resultPairScope({a:'A',b:'C'},{candidates}),'BETWEEN_REGION');
// Execute the actual import function with only its UI dependencies stubbed.
const elements=new Map();const el=id=>{if(!elements.has(id))elements.set(id,{value:'',checked:false,type:id==='p_screen_cross_primary'?'checkbox':'text',classList:{toggle(){}}});return elements.get(id);};
Object.assign(context,{$:el,S:{},B:{presets:{errfi1:{description:'default'}}},copy:x=>JSON.parse(JSON.stringify(x)),resolveImagerEntry:x=>x,loadedRefs(){},makeBins(){},renderAllowed(){},updateAdaptiveControls(){},markStale(){},toast(){}});
const rb=src.indexOf('function restoreConfiguration('),re=src.indexOf("$('saveConfig').onclick",rb);vm.runInContext(src.slice(rb,re),context);
const config={references:[{id:'X'}],imagers:[{id:'DP1'}],policy:{screen_cross_primary:true,threads:4}};context.restoreConfiguration(config,false);assert.equal(el('p_screen_cross_primary').checked,true);assert.equal(Number(el('p_threads').value),2);
config.policy={};context.restoreConfiguration(config,false);assert.equal(el('p_screen_cross_primary').checked,false);assert.equal(Number(el('p_threads').value),2);
console.log('UI saved-result and import behavior checks passed');
"""
        run = subprocess.run([node, '-e', script], cwd=ROOT, text=True, capture_output=True)
        self.assertEqual(run.returncode, 0, run.stdout + run.stderr)

    def test_ui_stopped_analysis_does_not_request_result_or_order(self):
        node = shutil.which('node')
        if not node or not Path(node).exists():
            self.skipTest('Node runtime unavailable')
        script = r"""
const fs=require('fs'), vm=require('vm'), assert=require('assert');
const src=fs.readFileSync('web/app.js','utf8');
const begin=src.indexOf('let pollFailures=0;'),end=src.indexOf("$('cancelRun').onclick",begin);
(async()=>{for(const status of ['failed','incomplete','stopped','cancelled']){
 const elements=new Map(),requests=[];const el=id=>{if(!elements.has(id))elements.set(id,{textContent:'',classList:{add(){}}});return elements.get(id);};
 const context={$:el,S:{job:'fixture'},fetch:async path=>{requests.push(path);return {ok:true,json:async()=>({status,stage:'Stopped',messages:['Memory guard stopped run'],error:'Available memory below reserve'})};},toast(){},refreshCheckpoints(){},setTimeout(){throw Error('Incomplete job must not retry result retrieval');},encodeURIComponent};
 vm.createContext(context);vm.runInContext(src.slice(begin,end),context);await context.pollJob();
 assert.deepEqual(requests,['/api/job?id=fixture']);assert(el('runState').textContent.startsWith('INCOMPLETE'));assert.equal(context.S.result,undefined);assert.equal(el('runDesign').disabled,false);
}})().catch(e=>{console.error(e);process.exit(1);});
"""
        run = subprocess.run([node, '-e', script], cwd=ROOT, text=True, capture_output=True)
        self.assertEqual(run.returncode, 0, run.stdout + run.stderr)


if __name__ == '__main__':
    unittest.main()
