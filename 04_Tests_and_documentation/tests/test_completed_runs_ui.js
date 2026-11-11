'use strict';
const {test}=require('node:test');
const assert=require('node:assert/strict');
const fs=require('node:fs');
const path=require('node:path');
const vm=require('node:vm');
const source=fs.readFileSync(path.join(__dirname,'../web/app.js'),'utf8');
const html=fs.readFileSync(path.join(__dirname,'../web/index.html'),'utf8');
const block=source.slice(source.indexOf('// Completed-run browsing'),source.indexOf('// End completed-run browsing.'));
const names=source.slice(source.indexOf('// Readable run names'),source.indexOf('// End readable run names.'));
const response=data=>({ok:true,json:async()=>data});
function setup(opts={}){
 const nodes={};const node=id=>nodes[id]||(nodes[id]={value:'',textContent:'',innerHTML:'',disabled:false,classList:{toggle(){}}});
 const requests=[],initialized=[],rendered=[],toasts=[];
 const run={id:'run_dp3',name:'DP3 · chr1:8186425-8206424',status:'complete',summary:{counts:{selected:129}}};
 const result={title:run.name,status:'PROVISIONAL_DRAFT_CROSS_PRIMARY_DISABLED',candidates:[],bins:[]};
 const context={B:{offline:!!opts.offline},S:{designActive:!!opts.active},$:node,fetch:async(url,...args)=>{requests.push(url);return opts.fetch?opts.fetch(url,...args):response(url==='/api/runs'?{runs:[run]}:result);},
  initResult:()=>initialized.push(1),renderResults:()=>rendered.push(1),toast:x=>toasts.push(x),esc:x=>String(x??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]))};
 vm.createContext(context);vm.runInContext(block,context);
 return {context,nodes,requests,initialized,rendered,toasts,refresh:()=>vm.runInContext('refreshCompletedRuns()',context),open:()=>vm.runInContext('openCompletedRun()',context),initialize:()=>vm.runInContext('initializeCompletedRuns()',context)};
}
test('completed run names load without starting designs',async()=>{const t=setup();t.initialize();assert.equal(t.requests.length,0);await t.refresh();assert.deepEqual(t.requests,['/api/runs']);assert.match(t.nodes.completedRunChoice.innerHTML,/DP3/);assert.match(t.nodes.completedRunChoice.innerHTML,/129 primaries/);assert.equal(t.nodes.openCompletedRun.disabled,false);});
test('only explicit open fetches and displays selected saved result',async()=>{const t=setup();await t.refresh();await t.open();assert.deepEqual(t.requests,['/api/runs','/api/result?id=run_dp3']);assert.equal(t.context.S.job,'run_dp3');assert.equal(t.initialized.length,1);assert.equal(t.rendered.length,1);assert.match(t.nodes.completedRunsStatus.textContent,/DRAFT/);});
test('active run blocks browsing result that would replace polling job',async()=>{const t=setup({active:true});await t.refresh();await t.open();assert.equal(t.requests.length,1);assert.equal(t.nodes.openCompletedRun.disabled,true);assert.equal(t.initialized.length,0);});
test('offline browsing has no network access',async()=>{const t=setup({offline:true});t.initialize();await t.refresh();await t.open();assert.equal(t.requests.length,0);assert.equal(t.nodes.openCompletedRun.disabled,true);});
test('failed and incomplete entries excluded and labels escaped',async()=>{const t=setup({fetch:async()=>response({runs:[{id:'bad',status:'incomplete'},{id:'done',name:'<script>name</script>',status:'complete'}]})});await t.refresh();assert.ok(!t.nodes.completedRunChoice.innerHTML.includes('value="bad"'));assert.match(t.nodes.completedRunChoice.innerHTML,/&lt;script&gt;/);});
test('failed lookup clears old selectable results',async()=>{let ok=true;const t=setup({fetch:async()=>{if(ok)return response({runs:[{id:'a',status:'complete'}]});throw Error('Disconnected');}});await t.refresh();ok=false;await t.refresh();await t.open();assert.equal(t.nodes.completedRunChoice.innerHTML,'');assert.equal(t.nodes.openCompletedRun.disabled,true);assert.match(t.nodes.completedRunsStatus.textContent,/Disconnected/);});
test('incomplete result never replaces current saved result',async()=>{const t=setup({fetch:async url=>response(url==='/api/runs'?{runs:[{id:'a',status:'complete'}]}:{status:'STOPPED - INCOMPLETE',candidates:[],bins:[]})});await t.refresh();await t.open();assert.equal(t.initialized.length,0);assert.equal(t.context.S.result,undefined);assert.match(t.nodes.completedRunsStatus.textContent,/incomplete/);});
test('malformed listing fails without guessing completion',async()=>{const t=setup({fetch:async()=>response({other:[]})});await t.refresh();assert.match(t.nodes.completedRunsStatus.textContent,/unreadable/);assert.equal(t.nodes.openCompletedRun.disabled,true);});
test('run name derives actual bins rather than stale preset or extra assignment',()=>{const context={S:{title:'ERRFI1 DP1',bins:[{id:'b1',imager:'DP3',chromosome:'chr1',genomic_start:8186425,genomic_end:8206424}],assignments:{b1:'DP3',b2:'DP11'}},$:()=>({value:''})};vm.createContext(context);vm.runInContext(names,context);assert.equal(vm.runInContext('proposedRunName()',context),'DP3 · chr1:8186425-8206424');});
test('concise default retains stable config key and opt-in detail',()=>{assert.match(html,/<option value="standard">Concise scores \(default\)<\/option>/);assert.match(html,/<option value="full_audit">Full audit — retain detailed evidence \(large\)<\/option>/);assert.match(html,/Both modes apply the same scientific checks/);assert.match(html,/id="downloadScores"/);});
