'use strict';
const {test}=require('node:test');
const assert=require('node:assert/strict');
const fs=require('node:fs');
const path=require('node:path');
const coverage=require('../web/probe_coverage.js');
function fixture() {
  return {title:'Saved probe coverage',order_status:'DRAFT',cross_primary_status:'NOT RUN - USER DISABLED',annotation:{status:'FETCHED'},regions:[{
    id:'r1',label:'Region one',imager:'DP1',reference:'context',start:201,end:2200,length:2000,chromosome:'chr1',genomic_start:8000201,genomic_end:8002200,mapping_orientation:'+',bound_strand:'+',
    sites:[{id:'c1',name:'Probe one',start:201,end:240,strand:'+'},{id:'c2',name:'Probe two',start:2161,end:2200,strand:'+'}],
    sectors:[{start:201,end:1200,count:1},{start:1201,end:2200,count:1}],metrics:{primaries:2,sites:2,left_half:1,right_half:1,mean_offset_bp:0,largest_gap_bp:1920,complete_intron_sites:1},
    genes:[{id:'g1',label:'Gene plus',strand:'+',spans:[[101,2300]],exons:[[101,260],[2100,2300]],transcripts:[{id:'tx1',label:'Isoform',strand:'+',start:101,end:2300,exons:[[101,260],[2100,2300]]}]},
      {id:'g2',label:'Gene minus',strand:'-',spans:[[1000,1600]],exons:[[1000,1100]],transcripts:[{id:'tx2',label:'Other isoform',strand:'-',start:1000,end:1600,exons:[[1000,1100]]}]}]
  }]};
}
test('inclusive intervals and base centres respect a nonzero local origin',()=>{
  const r={start:201,end:2200};
  assert.deepEqual(coverage.intervalToPixels([201,2200],r,100,1000),[100,1100]);
  assert.deepEqual(coverage.intervalToPixels([1,2500],r,100,1000),[100,1100]);
  assert.equal(coverage.intervalToPixels([1,200],r),null);
  assert.equal(coverage.pointToPixel(1200.5,r,100,1000),600);
  assert.equal(coverage.pointToPixel(201,r,100,1000),100.25);
});
test('rendered clipping keeps exon and span edges inside the shared region axis',()=>{
  const h=coverage.render(fixture());
  assert.ok(h.includes('x1="328"'));assert.ok(h.includes('x2="1288"'));
  assert.ok(h.includes('x="328"')&&h.includes('width="28.8"'));
  assert.ok(!h.includes('NaN'));assert.ok(h.includes('local 201–240; chr1:8,000,201–8,000,240'));
  assert.ok(h.includes('data-coverage-mean="r1" d="M808,303'));
});
test('gene arrows point in the correct direction without reversing coordinates',()=>{
  const h=coverage.render(fixture());
  assert.match(h, /data-transcript-strand="\+"[\s\S]*?d="M340,329L346,332L340,335"/);
  assert.match(h, /data-transcript-strand="-"[\s\S]*?d="M729\.52,354L723\.52,357L729\.52,360"/);
});
test('genomic tooltips follow explicit reverse mapping and omit unknown orientation',()=>{
  const f=fixture();f.regions[0].mapping_orientation='-';
  assert.ok(coverage.render(f).includes('local 201–240; chr1:8,002,161–8,002,200'));
  delete f.regions[0].mapping_orientation;
  assert.ok(!coverage.render(f).includes('local 201–240; chr1:'));
});
test('mean and RNA orientation labels remain reference-local for reverse references',()=>{
  const f=fixture();f.regions[0].mapping_orientation='-';f.regions[0].bound_strand='-';f.regions[0].metrics.mean_offset_bp=-40;
  const h=coverage.render(f);
  assert.ok(h.includes('Mean local offset  -40 bp'));assert.ok(h.includes('reference-local mean offset: -40 bp'));
  assert.ok(h.includes('DNA target - · complementary RNA -'));assert.ok(!h.includes('mean genomic offset'));
});
test('count scales stay common across regions, including a filtered detail view',()=>{
  const f=fixture(), second=structuredClone(f.regions[0]);second.id='r2';second.sectors[0].count=19;f.regions.push(second);
  assert.equal(coverage.prepare(f).countMaximum,20);
  assert.equal(coverage.prepare(f,{regionId:'r1',mode:'isoforms'}).countMaximum,20);
  const h=coverage.render(f,{regionId:'r1',mode:'isoforms'});
  assert.ok(h.includes('data-count-maximum="20"'));assert.ok(!h.includes('data-region-id="r2"'));
});
test('all dynamic labels and identifiers are XML escaped',()=>{
  const f=fixture();f.title='<script>bad()</script>';f.regions[0].id='" onclick="bad()';
  f.regions[0].sites[0].name='<img src=x onerror=bad()>';f.regions[0].genes[0].label='A&B <gene>';f.warnings=['<script>alert(1)</script>'];
  const h=coverage.render(f);
  assert.ok(!h.includes('<script>'));assert.ok(!h.includes('<img'));assert.ok(!h.includes(' onclick="'));
  assert.ok(h.includes('&lt;script&gt;'));assert.ok(h.includes('A&amp;B &lt;gene&gt;'));assert.ok(h.includes('&quot; onclick=&quot;'));
});
test('unavailable annotation is unknown rather than zero transcript overlap',()=>{
  const f=fixture();f.annotation={status:'NOT_LOADED'};f.regions[0].genes=[];f.regions[0].metrics.complete_intron_sites=0;
  const h=coverage.render(f);
  assert.ok(h.includes('Annotation unavailable — transcript overlap is unknown.'));
  assert.ok(h.includes('Complete-arm intron overlaps: unknown'));assert.ok(!h.includes('No overlapping transcripts in the saved annotation.'));
});
test('available empty annotation does not claim absence of transcription',()=>{
  const f=fixture();f.regions[0].genes=[];
  const h=coverage.render(f);assert.ok(h.includes('No overlapping transcripts in the saved annotation.'));
  assert.ok(h.includes('No annotated overlap is not an RNA-free guarantee.'));
});
test('partial exon models retain tracks without presenting line segments as certain introns',()=>{
  const f=fixture();f.regions[0].annotation={status:'PARTIAL'};f.regions[0].metrics.complete_intron_sites=null;
  const h=coverage.render(f,{mode:'isoforms'});
  assert.ok(h.includes('Annotation incomplete — some transcript or exon models are unavailable.'));
  assert.ok(h.includes('exon/intron identity may be incomplete'));assert.ok(h.includes('Complete-arm intron overlaps: unknown'));
  assert.ok(!h.includes('Connecting lines: introns.'));
  assert.equal((h.match(/data-transcript-id=/g)||[]).length,2);
});
test('small empty regions render a finite nonzero count scale without a mean marker',()=>{
  const f=fixture();Object.assign(f.regions[0],{start:1,end:1,sites:[],sectors:[],genes:[],metrics:{sites:0,primaries:0}});
  const h=coverage.render(f);assert.ok(h.includes('No selected binding sites'));assert.ok(h.includes('data-count-maximum="1"'));
  assert.ok(h.includes('Position within region (bp)'));assert.ok(!h.includes('NaN'));assert.ok(!h.includes('data-coverage-mean='));
});
test('every individual isoform is rendered, with optional transcript hiding',()=>{
  const f=fixture();const h=coverage.render(f,{mode:'isoforms'});
  assert.equal((h.match(/data-transcript-id=/g)||[]).length,2);
  assert.ok(h.includes('data-transcript-id="tx1"'));assert.ok(h.includes('data-transcript-id="tx2"'));
  assert.equal((coverage.render(f,{showTranscripts:false}).match(/data-transcript-strand=/g)||[]).length,0);
});
test('standalone output carries accessibility, white background, draft and scientific limits',()=>{
  const h=coverage.render(fixture());assert.ok(h.startsWith('<svg xmlns="http://www.w3.org/2000/svg"'));
  for(const token of ['role="img"','<title>','<desc>','fill="#ffffff"','font-family="Arial,Helvetica,sans-serif"','data-site-id="c1"','data-region-id="r1"','NOT RUN - USER DISABLED','not measured 3D centroid errors','annotation does not establish expression'])assert.ok(h.includes(token),token);
  assert.ok(!h.includes('foreignObject'));assert.ok(!h.includes('<image'));assert.ok(!h.includes('http://',h.indexOf('xmlns=')+50));
});
test('rendering does not mutate saved report, and malformed input fails visibly',()=>{
  const f=fixture(),before=JSON.stringify(f);coverage.render(f);assert.equal(JSON.stringify(f),before);
  for(const invalid of [null,{}, {regions:[{start:0,end:0}]}])assert.ok(coverage.render(invalid).includes('data-coverage-error="true"'));
});
test('known general warnings are consolidated while specific warnings survive',()=>{
  const f=fixture();f.warnings=['Coverage is a read-only view of saved selections; no design or QC was rerun.',
    'Cross-primary interactions within and between regions: NOT RUN - USER DISABLED.',
    'Custom saved evidence limitation.','Custom saved evidence limitation.'];
  const h=coverage.render(f);
  assert.ok(!h.includes('Coverage is a read-only view of saved selections;'));
  // One visible warning, and one copy in the accessible SVG description.
  assert.equal((h.match(/Custom saved evidence limitation\./g)||[]).length,2);
});
test('saved six-target annotation fixture retains all 64 target–isoform overlaps',t=>{
  const folder=path.resolve(__dirname,'../../outputs/final_six_target_review_20261021');
  const tf=path.join(folder,'transcript_tracks.json'),sf=path.join(folder,'coverage_review_summary.json');
  if(!fs.existsSync(tf)||!fs.existsSync(sf)){t.skip('Optional saved experiment fixture is not available in this checkout.');return;}
  const tracks=JSON.parse(fs.readFileSync(tf)),summary=JSON.parse(fs.readFileSync(sf));
  const regions=tracks.targets.map((r,i)=>{const s=summary.targets[i];return {id:r.target,label:r.target,start:r.start,end:r.end,imager:s.imager,
    sites:s.selected_site_intervals.map((p,j)=>({id:`${r.target}_${j}`,name:`Site ${j+1}`,start:p[0],end:p[1]})),
    sectors:s.sector_counts_1kb.map((count,j)=>({start:r.start+j*1000,end:r.start+(j+1)*1000-1,count})),metrics:{sites:s.sites,primaries:s.primaries,left_half:s.left_half,right_half:s.right_half,mean_offset_bp:s.equal_weight_mean_genomic_offset_bp,largest_gap_bp:s.largest_internal_arm_gap_bp,complete_intron_sites:s.potential_complete_arm_unspliced_rna_sites},
    genes:r.genes.map(g=>({...g,id:g.gene_id,spans:g.transcript_spans,exons:g.exon_union,transcripts:g.transcripts.map(x=>({...x,start:x.span[0],end:x.span[1]}))}))};});
  const report={title:'Six-target coverage',regions,annotation:{status:'FETCHED'}};
  const h=coverage.render(report,{mode:'isoforms'}),overview=coverage.render(report);
  assert.equal((h.match(/data-transcript-id=/g)||[]).length,64);
  assert.equal((h.match(/data-site-id=/g)||[]).length,594);
  assert.deepEqual(regions.map(r=>r.genes.reduce((sum,g)=>sum+g.transcripts.length,0)),[22,18,8,1,0,15]);
  assert.ok(Buffer.byteLength(overview)<200000,`Overview is ${Buffer.byteLength(overview)} bytes`);
});
