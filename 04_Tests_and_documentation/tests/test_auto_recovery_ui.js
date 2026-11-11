'use strict';
const {test}=require('node:test');
const assert=require('node:assert/strict');
const fs=require('node:fs');
const path=require('node:path');
const vm=require('node:vm');
const source=fs.readFileSync(path.join(__dirname,'../web/app.js'),'utf8');
const checkpointSource=source.slice(source.indexOf('// Saved checkpoints are listed'),source.indexOf('// End saved-checkpoint controls.'));
const pollSource=source.slice(source.indexOf('let pollFailures=0;'),source.indexOf('// SVG drawing:'));
const html=fs.readFileSync(path.join(__dirname,'../web/index.html'),'utf8');
const answer=data=>({ok:true,json:async()=>data});
const recovery=(state,extra={})=>({state,attempts:0,max_attempts:2,next_retry_at:Date.now()/1000+30,reason:'Worker interrupted.',latest_run_id:'first',origin_run_id:'first',...extra});
const status=(extra={})=>({status:'incomplete',stage:'Stopped - incomplete',messages:['BLAST batch 348 saved.'],...extra});
function setup(options={}){
 const nodes={},requests=[],posts=[],timers=[],toasts=[],restored=[],views=[],results=[];
 const node=id=>nodes[id]||(nodes[id]={textContent:'',innerHTML:'',value:'',disabled:false,scrollTop:0,scrollHeight:1,classList:{hidden:false,toggle(name,on){this[name]=on;},add(name){this[name]=true;},remove(name){this[name]=false;}}});
 const context={B:{offline:!!options.offline,token:'old_token'},S:{job:options.initialJob===undefined?'first':options.initialJob,designActive:false},$ :node,
  fetch:async(url,opts)=>{requests.push({url,opts});return options.fetch?options.fetch(url,opts):answer(status());},
  api:async(url,body)=>{posts.push({url,body,token:context.B.token});return options.api?options.api(url,body):{};},
  setTimeout:(fn,delay)=>{timers.push({fn,delay});return timers.length;},
  restoreConfiguration:c=>restored.push(c),renderRunSummary(){},go:n=>views.push(n),initResult:()=>results.push(context.S.result),toast:x=>toasts.push(x),fmt:String,esc:String,
  Date,Math,encodeURIComponent};
 vm.createContext(context);vm.runInContext(checkpointSource+'\n'+pollSource,context);
 const call=code=>vm.runInContext(code,context);
 return {nodes,context,requests,posts,timers,toasts,restored,views,results,
  poll:()=>call('pollJob()'),refresh:()=>call('refreshCheckpoints()'),cancel:()=>node('cancelRun').onclick(),initialize:()=>call('initializeCheckpoints()')};
}
test('pending retry stays incomplete, preserves cancel, and never submits a browser retry',async()=>{
 const t=setup({fetch:async()=>answer(status({auto_recovery:recovery('pending')}))});await t.poll();
 assert.equal(t.context.S.designActive,true);assert.equal(t.nodes.runDesign.disabled,true);assert.equal(t.nodes.cancelRun.classList.hidden,false);
 assert.match(t.nodes.runState.textContent,/INCOMPLETE.*automatic retry/);assert.match(t.nodes.autoRecoveryStatus.textContent,/retry 1 of 2/);assert.match(t.nodes.autoRecoveryStatus.textContent,/Worker interrupted/);
 assert.equal(t.timers[0].delay,1500);assert.equal(t.posts.length,0);assert.equal(t.results.length,0);assert.equal(t.toasts.length,0);
});
test('pending recovery takes precedence over persisted running status after a server restart',async()=>{
 const t=setup({fetch:async()=>answer(status({status:'running',auto_recovery:recovery('pending',{attempts:1})}))});await t.poll();
 assert.match(t.nodes.runState.textContent,/INCOMPLETE.*automatic retry/);assert.match(t.nodes.autoRecoveryStatus.textContent,/retry 2 of 2/);assert.equal(t.results.length,0);
});
test('failed worker remains watched during cleanup before its retry decision is persisted',async()=>{
 let polls=0;
 const t=setup({fetch:async()=>answer(status({auto_recovery:recovery(polls++?'pending':'active')}))});
 await t.poll();assert.match(t.nodes.runState.textContent,/INCOMPLETE.*checking automatic recovery/);assert.equal(t.context.S.designActive,true);
 assert.equal(t.nodes.cancelRun.classList.hidden,false);assert.equal(t.nodes.runDesign.disabled,true);assert.equal(t.timers[0].delay,1500);assert.equal(t.toasts.length,0);
 await t.timers[0].fn();assert.match(t.nodes.runState.textContent,/waiting for automatic retry/);assert.equal(t.timers.length,2);assert.equal(t.posts.length,0);assert.equal(t.results.length,0);
});
test('browser follows latest active attempt without requesting a new design',async()=>{
 const t=setup({fetch:async url=>answer(status(url.includes('id=first')?{auto_recovery:recovery('active',{attempts:1,latest_run_id:'retry_1'})}:{status:'running',messages:['Resumed at batch 349'],auto_recovery:recovery('active',{attempts:1,latest_run_id:'retry_1'})}))});
 await t.poll();assert.equal(t.context.S.job,'retry_1');assert.equal(t.timers[0].delay,0);await t.timers[0].fn();
 assert.equal(t.requests[1].url,'/api/job?id=retry_1');assert.equal(t.nodes.runLog.textContent,'Resumed at batch 349');assert.equal(t.posts.length,0);assert.equal(t.results.length,0);
});
test('completed recovery loads only the verified latest result',async()=>{
 const t=setup({fetch:async url=>{
  if(url==='/api/job?id=first')return answer(status({auto_recovery:recovery('complete',{attempts:1,latest_run_id:'retry_1'})}));
  if(url==='/api/job?id=retry_1')return answer(status({status:'complete',auto_recovery:recovery('complete',{attempts:1,latest_run_id:'retry_1'})}));
  if(url==='/api/result?id=retry_1')return answer({result_id:'retry_1',order_status:'DRAFT'});
  if(url==='/api/checkpoints')return answer({runs:[],active_design:false});
  throw Error('Unexpected URL '+url);
 }});await t.poll();assert.equal(t.results.length,0);await t.timers[0].fn();
 assert.equal(t.context.S.result.result_id,'retry_1');assert.equal(t.results.length,1);assert.deepEqual(t.views,[6]);assert.equal(t.posts.length,0);
 assert.ok(!t.requests.some(x=>x.url==='/api/result?id=first'));
});
test('recovery complete flag alone cannot turn an incomplete analysis into a pass',async()=>{
 const t=setup({fetch:async url=>answer(url==='/api/checkpoints'?{runs:[]}:status({auto_recovery:recovery('complete')}))});await t.poll();
 assert.match(t.nodes.runState.textContent,/INCOMPLETE/);assert.equal(t.results.length,0);assert.ok(!t.requests.some(x=>x.url.startsWith('/api/result')));
});
test('resource stop terminates UI polling and displays the stop reason without retry',async()=>{
 const t=setup({fetch:async url=>answer(url==='/api/checkpoints'?{runs:[]}:status({auto_recovery:recovery('stopped',{reason:'Memory safety guard stopped the worker.'})}))});await t.poll();
 assert.match(t.nodes.autoRecoveryStatus.textContent,/Memory safety guard/);assert.match(t.nodes.runState.textContent,/INCOMPLETE/);assert.equal(t.timers.length,0);assert.equal(t.posts.length,0);assert.equal(t.nodes.cancelRun.classList.hidden,true);
});
test('stopped recovery cannot display persisted pre-crash running state as active',async()=>{
 const t=setup({fetch:async url=>answer(url==='/api/checkpoints'?{runs:[]}:status({status:'running',auto_recovery:recovery('stopped',{reason:'Automatic retry limit reached.'})}))});await t.poll();
 assert.match(t.nodes.runState.textContent,/INCOMPLETE.*stopped/);assert.equal(t.context.S.designActive,false);assert.equal(t.timers.length,0);assert.equal(t.results.length,0);
});
test('exhausted recovery follows the last failed attempt so its actual error is shown',async()=>{
 const t=setup({fetch:async url=>answer(url==='/api/checkpoints'?{runs:[]}:status({error:url.includes('retry_2')?'BLAST exited 2.':undefined,auto_recovery:recovery('stopped',{attempts:2,latest_run_id:'retry_2',reason:'Automatic retry limit reached.'})}))});
 await t.poll();await t.timers[0].fn();assert.equal(t.context.S.job,'retry_2');assert.equal(t.toasts[0],'BLAST exited 2.');assert.equal(t.timers.length,1);assert.equal(t.posts.length,0);
});
test('reload reattaches to retry-wait, restores original DP3 config, and refreshes token read-only',async()=>{
 const config={assignments:{DP3_region:'DP3'},policy:{screen_cross_primary:false,shortlist_per_bin:400,extension_shortlist_per_bin:200}};
 const t=setup({initialJob:null,fetch:async url=>{
  if(url==='/api/checkpoints')return answer({runs:[],active_design:true,active_run:{id:'first',status:'retry_wait'}});
  if(url==='/api/active-design')return answer({id:'first',config,token:'restarted_token'});
  return answer(status({auto_recovery:recovery('pending')}));
 }});await t.refresh();await new Promise(resolve=>setImmediate(resolve));
 assert.equal(t.context.S.job,'first');assert.equal(t.context.S.designActive,true);assert.equal(t.context.B.token,'restarted_token');assert.deepEqual(t.restored,[config]);assert.deepEqual(t.views,[5]);assert.equal(t.posts.length,0);assert.match(t.nodes.runState.textContent,/automatic retry/);
});
test('refresh with an existing watched job cannot restore different settings or duplicate its poll',async()=>{
 const t=setup({fetch:async()=>answer({runs:[],active_design:true,active_run:{id:'retry_1'}})});await t.refresh();
 assert.equal(t.context.S.job,'first');assert.equal(t.restored.length,0);assert.equal(t.requests.length,1);assert.equal(t.posts.length,0);assert.equal(t.timers.length,0);
});
test('connection recovery refreshes server token without restarting a worker',async()=>{
 let interrupted=true;
 const t=setup({fetch:async url=>{
  if(interrupted){interrupted=false;throw Error('server restarting');}
  if(url==='/api/active-design')return answer({id:'first',token:'new_token'});
  return answer(status({auto_recovery:recovery('pending')}));
 }});await t.poll();await t.timers[0].fn();
 assert.equal(t.context.B.token,'new_token');assert.ok(t.requests.some(x=>x.url==='/api/active-design'));assert.equal(t.posts.length,0);assert.match(t.nodes.runState.textContent,/automatic retry/);
});
test('cancel during retry-wait refreshes token and cancels the whole existing chain',async()=>{
 const t=setup({fetch:async url=>answer(url==='/api/active-design'?{id:'retry_1',token:'new_token'}:status({auto_recovery:recovery('pending')}))});await t.poll();await t.cancel();
 assert.equal(t.posts.length,1);assert.equal(t.posts[0].url,'/api/cancel');assert.equal(t.posts[0].body.id,'first');assert.equal(t.posts[0].token,'new_token');assert.equal(t.nodes.cancelRun.disabled,false);
 assert.match(t.nodes.runState.textContent,/automatic recovery will stop/);assert.match(t.toasts[0],/pending automatic retries/);
});
test('unconfirmed cancellation remains visible and never claims the worker stopped',async()=>{
 const t=setup({fetch:async()=>{throw Error('server offline');}});await t.cancel();
 assert.equal(t.posts.length,0);assert.equal(t.nodes.cancelRun.disabled,false);assert.match(t.toasts[0],/Cancellation was not confirmed.*server offline/);
});
test('duplicate cancellation clicks cannot submit concurrent requests',async()=>{
 let done;const t=setup({fetch:()=>new Promise(resolve=>done=resolve)});const first=t.cancel();await t.cancel();assert.equal(t.requests.length,1);
 done(answer({id:'first',token:'current'}));await first;assert.equal(t.posts.length,1);
});
test('stale scheduled polls cannot reopen or overwrite a newer design',async()=>{
 const t=setup({fetch:async()=>answer(status({auto_recovery:recovery('pending')}))});await t.poll();t.context.S.job='new_design';await t.timers[0].fn();
 assert.equal(t.requests.length,1);assert.equal(t.context.S.job,'new_design');
});
test('stale result response cannot replace the newer design result',async()=>{
 let done;const t=setup({fetch:async url=>url.startsWith('/api/result')?new Promise(resolve=>done=resolve):answer(status({status:'complete'}))});
 const waiting=t.poll();await new Promise(resolve=>setImmediate(resolve));t.context.S.job='new_design';done(answer({result_id:'first'}));await waiting;
 assert.equal(t.context.S.result,undefined);assert.equal(t.results.length,0);
});
test('offline saved result makes no recovery calls',async()=>{
 const t=setup({offline:true,initialJob:null});t.initialize();await t.refresh();assert.equal(t.requests.length,0);assert.equal(t.posts.length,0);
});
test('page explains two retries and the conditions that stop automatic recovery',()=>{
 assert.match(html,/Automatic restart: ON for new Standard runs/);assert.match(html,/at most twice: after 30 seconds, then after 2 minutes/);
 assert.match(html,/Cancellation, memory\/resource stops, and damaged or incompatible checkpoints stop automatic recovery/);
 assert.match(html,/analysis remains INCOMPLETE/);assert.match(html,/manual resume option/);
});
