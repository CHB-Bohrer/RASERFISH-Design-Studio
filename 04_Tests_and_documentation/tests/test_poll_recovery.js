'use strict';
const {test}=require('node:test');
const assert=require('node:assert/strict');
const fs=require('node:fs');
const vm=require('node:vm');
const source=fs.readFileSync(require('node:path').join(__dirname,'../web/app.js'),'utf8');
const pollSource=source.slice(source.indexOf('let pollFailures=0;'),source.indexOf("$('cancelRun').onclick="));
function setup(fetch){
 const nodes={};const timers=[];const toasts=[];
 const context={S:{job:'saved_run'},fetch,encodeURIComponent,Error,Math,
  $:id=>nodes[id]||(nodes[id]={textContent:'',disabled:false,classList:{add(){},remove(){}},scrollTop:0,scrollHeight:1}),
  setTimeout:(fn,delay)=>timers.push({fn,delay}),toast:s=>toasts.push(s),refreshCheckpoints(){},initResult(){throw Error('No completion expected');},go(){throw Error('No completion expected');}};
 vm.createContext(context);vm.runInContext(pollSource,context);
 context.$('runLog').textContent='BLAST: human; batch 348/600 (1 candidates).';context.$('progressStage').textContent='Reference searches';
 return {context,nodes,timers,toasts,poll:()=>vm.runInContext('pollJob()',context)};
}
test('connection failures preserve last progress without repeated log lines and back off',async()=>{
 const t=setup(async()=>{throw Error('offline');});const before=t.nodes.runLog.textContent;
 for(let i=0;i<15;i++)await t.poll();
 assert.equal(t.nodes.runLog.textContent,before);assert.equal(t.nodes.progressStage.textContent,'Reference searches');
 assert.match(t.nodes.runState.textContent,/completion is unverified/);assert.equal(t.nodes.runDesign.disabled,true);
 assert.equal(t.timers[0].delay,1500);assert.equal(t.timers.at(-1).delay,15000);
});
test('retry recovers saved incomplete status without loading completed results',async()=>{
 let calls=0;const t=setup(async()=>{if(calls++===0)throw Error('offline');return {ok:true,json:async()=>({status:'incomplete',stage:'Stopped - incomplete',messages:['Saved batch 348'],error:'Design supervisor exited.'})};});
 await t.poll();await t.timers[0].fn();
 assert.equal(calls,2);assert.match(t.nodes.runState.textContent,/INCOMPLETE/);assert.equal(t.nodes.runLog.textContent,'Saved batch 348');
 assert.equal(t.nodes.runDesign.disabled,false);assert.equal(t.toasts[0],'Design supervisor exited.');assert.equal(t.timers.length,1);
});
test('definitive status error is visible without indefinite retry or changing progress',async()=>{
 const t=setup(async()=>({ok:false,status:400,json:async()=>({error:'This run is unverified.'})}));await t.poll();
 assert.equal(t.timers.length,0);assert.match(t.nodes.runState.textContent,/This run is unverified/);assert.match(t.nodes.runLog.textContent,/348\/600/);
});
test('stale response cannot overwrite a newer run',async()=>{
 let resolve;const t=setup(()=>new Promise(r=>resolve=r));const waiting=t.poll();t.context.S.job='new_run';
 resolve({ok:true,json:async()=>({status:'complete'})});await waiting;
 assert.match(t.nodes.runLog.textContent,/348\/600/);assert.equal(t.timers.length,0);
});
