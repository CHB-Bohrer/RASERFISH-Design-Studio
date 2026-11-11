'use strict';
const {test} = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const binding = require('../web/imager_binding.js');
const rc = seq => [...seq].reverse().map(c=>({A:'T',T:'A',G:'C',C:'G'}[c])).join('');
const I='CCGATACGAATA', D='TATTCGTATCGG';
function entry(seq=I){return {id:'DP7',sequence:seq,docking:rc(seq),dye:'Atto565N'};}
test('DP7 actual canonical sequences and reversed display',()=>{
 const p=binding.pairData(entry());
 assert.equal(p.sequence,I); assert.equal(p.docking,D);
 assert.equal(p.imagerDisplay,'ATAAGCATAGCC');
 assert.equal(p.columns[0].dockingBase,'T'); assert.equal(p.columns[0].imagerBase,'A');
 assert.equal(p.columns[0].imagerPosition,12); assert.equal(p.columns.at(-1).imagerPosition,1);
});
test('all six inventory sequences get one correct column per base',()=>{
 const inv=JSON.parse(fs.readFileSync(path.join(__dirname,'..','imagers.json')));
 assert.equal(inv.length,6);
 for(const im of inv){const p=binding.pairData({...im,docking:rc(im.sequence)});
  assert.equal(p.columns.length,im.sequence.length);
  for(const c of p.columns){assert.equal(c.complementary,true);assert.equal(im.sequence[c.imagerPosition-1],c.imagerBase);}
 }
});
test('every supported length 8 through 40 with mixed nonpalindromic input',()=>{
 for(let n=8;n<=40;n++){const seq='ACCTATGGTACAGTCCGATAGGCTCAGTCAACGGTTGCACT'.slice(0,n);const p=binding.pairData(entry(seq));
  assert.equal(p.columns.length,n); assert.equal(p.imagerDisplay,[...seq].reverse().join(''));
  assert.ok(p.columns.every(x=>x.complementary));
  const h=binding.render(entry(seq));assert.equal((h.match(/data-base-index=/g)||[]).length,n);
 }
});
test('canonical input object is not changed',()=>{
 const e=entry();const before=JSON.stringify(e);binding.render(e,{showOverview:true});assert.equal(JSON.stringify(e),before);
});
test('inconsistent pairs show a clear error and no duplex',()=>{
 const im={...entry(),docking:'A'.repeat(12)};assert.throws(()=>binding.pairData(im));
 const h=binding.render(im);assert.ok(h.includes('data-binding-error'));assert.ok(!h.includes('data-binding-svg'));
});
test('conflicting cached aliases are rejected',()=>{
 assert.throws(()=>binding.pairData({...entry(),imager_sequence_5to3:'T'.repeat(12)}));
});
test('missing sequence is not invented from name',()=>{
 assert.throws(()=>binding.pairData({id:'DP7'}));
});
test('ambiguous, modified and invalid bases fail visibly',()=>{
 for(const sequence of ['ACGTNNNNACGT',I+'/dye/',"5'-"+I+"-3'",'A'.repeat(7),'A'.repeat(41)]){
  assert.ok(binding.render({...entry(),sequence}).includes('data-binding-error'));
 }
});
test('8 and 40 nt edges are included',()=>{for(const n of [8,40])assert.equal(binding.pairData(entry('A'.repeat(n))).columns.length,n);});
test('lower-case and whitespace normalize without mutating source',()=>{
 const im={sequence:'ccgat acgaata',docking:'tattcgtatcgg'};const p=binding.pairData(im);assert.equal(p.sequence,I);assert.equal(im.sequence,'ccgat acgaata');
});
test('actual primary assembly is checked before displaying architecture',()=>{
 const c={prefix:D+'TT',suffix:'TT'+D,arm:'ACGTACCTGAGC'};c.full=c.prefix+c.arm+c.suffix;
 const h=binding.render(entry(),{primary:c,showOverview:true});assert.ok(!h.includes('data-binding-error'));assert.ok(h.includes(c.arm.length+' bases, from this candidate'));
 assert.ok(binding.render(entry(),{primary:{...c,full:c.full+'T'}}).includes('data-binding-error'));
});
test('opposite DNA arm orientations do not reverse docking/imager pairing',()=>{
 for(const arm of ['ACCTGTAGGCAT',rc('ACCTGTAGGCAT')]){
  const c={prefix:D+'TT',suffix:'TT'+D,arm};c.full=c.prefix+arm+c.suffix;
  assert.ok(binding.render(entry(),{primary:c}).includes('data-docking-5to3="'+D+'"'));
 }
});
test('user names and dye notes are escaped',()=>{
 const h=binding.render({...entry(),id:'<img src=x onerror=alert(1)>',dye:'<script>bad()</script>'});
 assert.ok(!h.includes('<img'));assert.ok(!h.includes('<script>'));assert.ok(h.includes('&lt;script&gt;'));
});
test('exported pairing has self-contained colors, labels and identity',()=>{
 const h=binding.render(entry());assert.ok(h.includes('xmlns="http://www.w3.org/2000/svg"'));
 assert.ok(h.includes('data-imager-5to3="'+I+'"'));assert.ok(h.includes('data-docking-5to3="'+D+'"'));
 assert.ok(h.includes('data-binding-export="duplex"'));assert.ok(h.includes('not a folding'));
 assert.ok(h.includes('3\u2032 to 5\u2032'));assert.ok(h.includes('5\u2032 to 3\u2032'));
});
test('hover/focus groups have correct numbering and accessible descriptions',()=>{
 const h=binding.render(entry());assert.equal((h.match(/tabindex="0"/g)||[]).length,12);
 assert.ok(h.includes('data-docking-position="1" data-imager-position="12"'));
 assert.ok(h.includes('data-docking-position="12" data-imager-position="1"'));
 assert.ok(h.includes('aria-label='));
});
test('full-primary cartoon depicts two separate imagers and distinct central arm',()=>{
 const h=binding.render(entry(),{showOverview:true});assert.ok(h.includes('Imager copy 1'));assert.ok(h.includes('Imager copy 2'));
 assert.ok(h.includes('DNA-binding arm'));assert.ok(h.includes('Dye attachment end is not inferred'));
});
