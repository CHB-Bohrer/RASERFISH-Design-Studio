'use strict';
const {test}=require('node:test');
const assert=require('node:assert/strict');
const fs=require('node:fs');
const path=require('node:path');
const vm=require('node:vm');
const source=fs.readFileSync(path.join(__dirname,'../web/app.js'),'utf8');
const checkpointSource=source.slice(source.indexOf('// Saved checkpoints are listed'),source.indexOf('// End saved-checkpoint controls.'));
const html=fs.readFileSync(path.join(__dirname,'../web/index.html'),'utf8');
const run={id:'20261017_saved',summary:'DP3 — chr1:8186425–8206424',status:'incomplete',created_utc:'2026-10-17T12:00:00Z',completed_batches:348,total_batches:1200,folded_candidates:4,total_candidates:600};
const response=data=>({ok:true,json:async()=>data});
function setup(options={}){
 const nodes={},requests=[],resumed=[],restored=[],polls=[],toasts=[];
 const context={B:{offline:!!options.offline},S:{job:null,designActive:false},pollFailures:0,
  $:id=>nodes[id]||(nodes[id]={textContent:'',innerHTML:'',value:'',disabled:false,classList:{hidden:false,toggle(name,on){this[name]=on;},add(name){this[name]=true;},remove(name){this[name]=false;}}}),
  fetch:async(...args)=>{requests.push(args);return options.fetch?options.fetch(...args):response({runs:[run],active_run:null,active_design:false});},
  api:async(...args)=>{resumed.push(args);return options.api?options.api(...args):{id:'new_attempt',config:{references:[{id:'saved_reference'}],policy:{screen_cross_primary:false}}};},
  restoreConfiguration:c=>{restored.push(c);if(options.restoreError)throw Error('Display error');},renderRunSummary(){},
  pollJob(){polls.push(context.S.job);},toast:s=>toasts.push(s),
  fmt:String,esc:x=>String(x??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]))};
 vm.createContext(context);vm.runInContext(checkpointSource,context);
 const call=name=>vm.runInContext(name+'()',context);
 return {context,nodes,requests,resumed,restored,polls,toasts,refresh:()=>call('refreshCheckpoints'),resume:()=>call('resumeCheckpoint'),initialize:()=>call('initializeCheckpoints')};
}
test('page reload lists stopped checkpoints without launching an automatic retry',async()=>{
 const t=setup();t.initialize();await new Promise(resolve=>setImmediate(resolve));
 assert.equal(t.requests.length,1);assert.equal(t.requests[0][0],'/api/checkpoints');assert.equal(t.resumed.length,0);assert.equal(t.polls.length,0);
 assert.equal(t.nodes.checkpointPanel.classList.hidden,false);assert.equal(t.nodes.resumeCheckpoint.disabled,false);
 assert.match(t.nodes.checkpointStatus.textContent,/Automatic recovery is managed by the local program/);
 assert.match(t.nodes.checkpointDetail.textContent,/INCOMPLETE/);assert.match(t.nodes.checkpointDetail.textContent,/348 \/ 1200/);assert.match(t.nodes.checkpointDetail.textContent,/4 \/ 600/);
 assert.equal(typeof t.nodes.resumeCheckpoint.onclick,'function');
});
test('offline results hide checkpoint controls and perform no network requests',async()=>{
 const t=setup({offline:true});t.initialize();await t.refresh();await t.resume();
 assert.equal(t.nodes.checkpointPanel.classList.hidden,true);assert.equal(t.nodes.resumeCheckpoint.disabled,true);assert.equal(t.requests.length,0);assert.equal(t.resumed.length,0);
});
test('empty checkpoint list is visible and cannot be resumed',async()=>{
 const t=setup({fetch:async()=>response({runs:[]})});await t.refresh();await t.resume();
 assert.match(t.nodes.checkpointStatus.textContent,/No resumable saved checkpoints/);assert.equal(t.nodes.resumeCheckpoint.disabled,true);assert.equal(t.resumed.length,0);
});
test('missing endpoint produces a useful restart message without a resume request',async()=>{
 const t=setup({fetch:async()=>({ok:false,status:404,json:async()=>({error:'Not found'})})});await t.refresh();await t.resume();
 assert.match(t.nodes.checkpointStatus.textContent,/Restart the updated Studio/);assert.equal(t.nodes.resumeCheckpoint.disabled,true);assert.equal(t.resumed.length,0);
});
test('server listing failure is shown without treating checkpoints as complete',async()=>{
 const t=setup({fetch:async()=>({ok:false,status:500,json:async()=>({error:'Checkpoint index unavailable'})})});await t.refresh();
 assert.equal(t.nodes.checkpointStatus.textContent,'Checkpoint index unavailable');assert.equal(t.nodes.resumeCheckpoint.disabled,true);assert.equal(t.nodes.checkpointDetail.textContent,'');
});
test('connection failure clears a stale selectable checkpoint',async()=>{
 let connected=true;const t=setup({fetch:async()=>{if(connected)return response({runs:[run]});throw Error('Connection unavailable');}});
 await t.refresh();connected=false;await t.refresh();await t.resume();
 assert.equal(t.nodes.checkpointChoice.innerHTML,'');assert.equal(t.nodes.resumeCheckpoint.disabled,true);assert.equal(t.resumed.length,0);assert.equal(t.nodes.checkpointStatus.textContent,'Connection unavailable');
});
test('malformed listing fails visibly rather than guessing checkpoint eligibility',async()=>{
 const t=setup({fetch:async()=>response({not_runs:[run]})});await t.refresh();
 assert.match(t.nodes.checkpointStatus.textContent,/unreadable/);assert.equal(t.nodes.resumeCheckpoint.disabled,true);
});
test('only incomplete runs are offered and names are escaped',async()=>{
 const t=setup({fetch:async()=>response({runs:[null,{id:'complete',status:'complete'},{id:'running',status:'running'},{...run,summary:'<script>malicious</script>'}]})});await t.refresh();
 assert.ok(!t.nodes.checkpointChoice.innerHTML.includes('<script>'));assert.match(t.nodes.checkpointChoice.innerHTML,/&lt;script&gt;/);assert.ok(!t.nodes.checkpointChoice.innerHTML.includes('value="complete"'));assert.ok(!t.nodes.checkpointChoice.innerHTML.includes('value="running"'));
});
test('explicit resume sends only saved id, restores returned settings, and polls new attempt',async()=>{
 const t=setup();await t.refresh();await t.resume();
 assert.equal(t.resumed.length,1);assert.equal(t.resumed[0][0],'/api/resume');assert.equal(JSON.stringify(t.resumed[0][1]),JSON.stringify({id:run.id}));
 assert.equal(t.restored.length,1);assert.equal(t.restored[0].references[0].id,'saved_reference');assert.equal(t.restored[0].policy.screen_cross_primary,false);
 assert.deepEqual(t.polls,['new_attempt']);assert.equal(t.context.S.job,'new_attempt');assert.equal(t.context.S.designActive,true);assert.equal(t.nodes.runDesign.disabled,true);
 assert.match(t.nodes.runLog.textContent,/original inputs and scientific settings/);assert.match(t.nodes.checkpointStatus.textContent,/remains incomplete/);assert.equal(t.nodes.resumeCheckpoint.disabled,true);
});
test('active local or cross-instance design blocks resume',async()=>{
 for(const active of [{active_run:{id:'running'}},{active_design:true}]){
  const t=setup({fetch:async()=>response({runs:[run],...active})});await t.refresh();await t.resume();
  assert.equal(t.nodes.resumeCheckpoint.disabled,true);assert.equal(t.resumed.length,0);assert.match(t.nodes.checkpointStatus.textContent,/already active/);
 }
 const t=setup();t.context.S.designActive=true;await t.refresh();await t.resume();assert.equal(t.resumed.length,0);assert.equal(t.nodes.resumeCheckpoint.disabled,true);
});
test('pending resume disables duplicate requests',async()=>{
 let resolve;const t=setup({api:()=>new Promise(r=>resolve=r)});await t.refresh();const first=t.resume();await t.resume();
 assert.equal(t.resumed.length,1);assert.equal(t.nodes.resumeCheckpoint.disabled,true);assert.equal(t.nodes.refreshCheckpoints.disabled,true);
 resolve({id:'new_attempt',config:{references:[{id:'r'}]}});await first;assert.equal(t.polls.length,1);
});
test('failed or unconfirmed resume shows error and requires refresh before retry',async()=>{
 const t=setup({api:async()=>{throw Error('Saved input fingerprint mismatch');}});await t.refresh();await t.resume();await t.resume();
 assert.equal(t.resumed.length,1);assert.equal(t.polls.length,0);assert.equal(t.restored.length,0);assert.equal(t.nodes.runDesign.disabled,false);assert.equal(t.nodes.resumeCheckpoint.disabled,true);
 assert.match(t.nodes.checkpointStatus.textContent,/Resume request was not confirmed/);assert.match(t.nodes.checkpointStatus.textContent,/fingerprint mismatch/);
 await t.refresh();assert.equal(t.nodes.resumeCheckpoint.disabled,false);
});
test('missing resume identifier never polls or calls restoreConfiguration',async()=>{
 const t=setup({api:async()=>({config:{references:[]}})});await t.refresh();await t.resume();
 assert.equal(t.polls.length,0);assert.equal(t.restored.length,0);assert.match(t.nodes.checkpointStatus.textContent,/did not return a run identifier/);
});
test('display-only restore failure still monitors the already accepted saved run',async()=>{
 const t=setup({restoreError:true});await t.refresh();await t.resume();
 assert.deepEqual(t.polls,['new_attempt']);assert.equal(t.context.S.designActive,true);assert.match(t.toasts[0],/saved run has started/);
});
test('markup explains saved settings replace current edits and does not imply a completed pass',()=>{
 assert.match(html,/id="checkpointPanel" class="notice hidden"/);assert.match(html,/Resume saved progress/);
 assert.match(html,/Current edits in the controls are replaced/);assert.match(html,/not a PASS or a completed order/);
});
