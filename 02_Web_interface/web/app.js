'use strict';
const B=window.STUDIO_BOOT||{};
const $=id=>document.getElementById(id);
const copy=x=>JSON.parse(JSON.stringify(x));
const esc=x=>String(x??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const fmt=x=>x==null?'not available':Number(x).toLocaleString(undefined,{maximumFractionDigits:1});
const rc=s=>s.split('').reverse().map(c=>({A:'T',T:'A',C:'G',G:'C',N:'N',R:'Y',Y:'R',S:'S',W:'W',K:'M',M:'K',B:'V',V:'B',D:'H',H:'D'}[c]||c)).join('');
const complement=s=>rc(s).split('').reverse().join('');
const COLORS=['#65d6c4','#bba0f0','#e7b16d','#85b9f4','#a3d38f','#f29aba','#c5b495','#90c9df'];
const color=id=>COLORS[Math.max(0,S.imagers.findIndex(x=>x.id===id))%COLORS.length];
const S={refs:copy(B.references||[]),imagers:copy(B.imagers||[]),bins:[],tails:[],assignments:{},allowed:[],transcription:{},preset:'errfi1',
 tracks:null,annotationProvenance:[],result:B.demo?copy(B.demo):null,job:null,step:1,region:null,viewRef:null,viewStart:1,viewEnd:1,selected:null,stale:false};
S.imagers=S.imagers.map(resolveImagerEntry);
const STEPS=['Targets & features','Tracing regions','Imaging probes','QC & selection','Run design','Compare results'];
$('steps').innerHTML=STEPS.map((s,i)=>`<button data-step-nav="${i+1}" ${i===0?'class="active"':''}>${i+1}. ${s}</button>`).join('');
function toast(m){$('toast').textContent=m;$('toast').classList.remove('hidden');setTimeout(()=>$('toast').classList.add('hidden'),7000);}
function markStale(){S.stale=true;if(S.result)$('staleResult').classList.remove('hidden');}
function help(key){const g=B.glossary||{};const items=key&&g[key]?[[key,g[key]]]:Object.entries(g);$('helpTitle').textContent=key&&g[key]?g[key].label:'Terms explained';$('helpText').innerHTML=items.map(([k,v])=>`<section class="term-block"><h3>${esc(v.label)}</h3><p>${esc(v.text)}</p></section>`).join('');$('helpDialog').showModal();}
document.addEventListener('click',e=>{const h=e.target.closest('[data-help]');if(h){e.preventDefault();help(h.dataset.help);}const n=e.target.closest('[data-step-nav],[data-next]');if(n){e.preventDefault();go(Number(n.dataset.stepNav||n.dataset.next));}});
$('helpAll').onclick=()=>help();$('closeHelp').onclick=()=>$('helpDialog').close();
function go(n){S.step=n;document.querySelectorAll('.step').forEach(e=>e.classList.toggle('hidden',Number(e.dataset.step)!==n));document.querySelectorAll('[data-step-nav]').forEach(e=>e.classList.toggle('active',Number(e.dataset.stepNav)===n));if(n===2){renderBins();renderBinMap();}if(n===3)renderImagers();if(n===5){refreshCheckpoints();renderRunSummary();const issues=transcriptionIssues();$('runSummary').innerHTML+='<div class="notice '+(issues.length?'danger':'')+'"><b>RNA policy</b><br>'+esc(issues.length?issues.join(' '):rnaPolicyText())+'</div>';}if(n===6){refreshCompletedRuns();renderResults();coveragePanel?.refreshSources();}window.scrollTo({top:0,behavior:'smooth'});}
async function api(path,body){if(B.offline)throw Error('This is an offline preview/results file. Start the local Design Studio program to fetch annotations or run QC-guided selection.');const r=await fetch(path,{method:'POST',headers:{'Content-Type':'application/json','X-Studio-Token':B.token},body:JSON.stringify(body)});let d;try{d=await r.json();}catch{throw Error('The local program returned an unreadable response. Keep its Terminal open.');}if(!r.ok)throw Error(d.error||'Request failed');return d;}
async function task(button,status,fn){const btn=$(button),el=$(status);btn.disabled=true;el.classList.remove('error');el.textContent='Working...';try{await fn();}catch(e){el.classList.add('error');el.textContent=e.message;toast(e.message);}finally{btn.disabled=false;}}
function parseFasta(text){text=text.trim();if(!text.startsWith('>'))text='>sequence\n'+text;let refs=[],r=null;for(const line of text.split(/\r?\n/)){if(!line.trim())continue;if(line.startsWith('>')){if(r)refs.push(r);let h=line.slice(1).trim();r={id:h.split(/\s/)[0],description:h,sequence:'',mapping:null,features:[]};}else{if(!r)throw Error('Missing FASTA header');r.sequence+=line.replace(/\s/g,'').toUpperCase();}}if(r)refs.push(r);if(!refs.length||new Set(refs.map(r=>r.id)).size!==refs.length)throw Error('FASTA needs unique record IDs.');for(const r of refs){if(!r.sequence||/[^ACGTRYSWKMBDHVN]/.test(r.sequence))throw Error('Invalid reference letters. Do not put notes in DNA text.');r.length=r.sequence.length;}return refs;}
function refOptions(select,refs=S.refs){const old=$(select).value;$(select).innerHTML=refs.map(r=>`<option value="${esc(r.id)}">${esc(r.id)} (${fmt(r.length)} bp)</option>`).join('');if(refs.some(r=>r.id===old))$(select).value=old;}
function loadedRefs(){refOptions('annotationRef');refOptions('contextRef');$('referenceCards').innerHTML=S.refs.map(r=>`<div class="refcard"><b>${esc(r.id)}</b>${fmt(r.length)} bp<br>${r.mapping?'GRCh38 '+esc(r.mapping.chromosome)+':'+fmt(r.mapping.start)+'-'+fmt(r.mapping.end):'Local sequence; no genomic coordinates inferred'}</div>`).join('');$('boundsEditor').innerHTML=S.refs.map(r=>`<div class="grid3"><div>${esc(r.id)}</div><label>Start<input type="number" id="boundA_${esc(r.id)}" value="1" min="1" max="${r.length}"></label><label>End<input type="number" id="boundB_${esc(r.id)}" value="${r.length}" min="1" max="${r.length}"></label></div>`).join('');updateAnnotationAvailability();renderTranscription();contextMap();makeBins();}
$('inputMode').onchange=()=>{for(const [id,mode] of [['locusBox','genomic'],['fastaBox','fasta'],['annotatedBox','annotated']])$(id).classList.toggle('hidden',$('inputMode').value!==mode);};
$('fastaFile').onchange=async()=>{const f=$('fastaFile').files[0];if(f)$('fastaText').value=await f.text();};
$('loadSequence').onclick=()=>task('loadSequence','loadStatus',async()=>{const mode=$('inputMode').value;let d;
 if(mode==='bundled')d={references:copy(B.references)};
 else if(mode==='legacy')d=await api('/api/references',{mode:'bundled'});
 else if(mode==='fasta')d=B.offline?{references:parseFasta($('fastaText').value)}:await api('/api/references',{mode,text:$('fastaText').value});
 else if(mode==='genomic')d=await api('/api/references',{mode,locus:$('locus').value});
 else {const f=$('annotatedFile').files[0];if(!f)throw Error('Choose an annotated file first.');const a=new Uint8Array(await f.arrayBuffer());let s='';for(let i=0;i<a.length;i+=8192)s+=String.fromCharCode(...a.slice(i,i+8192));d=await api('/api/references',{mode,filename:f.name,base64:btoa(s)});}
 S.refs=d.references;S.assignments={};S.annotationProvenance=[];S.allowed=[];S.transcription={};
 if(mode==='genomic'&&!B.offline){for(const r of S.refs){try{const g=await api('/api/genes',{reference:r});r.features=(r.features||[]).concat(g.features);r.gene_annotation={...g,features:undefined};S.annotationProvenance.push({...g,features:undefined});}catch(e){r.gene_annotation={status:'FAILED',error:e.message};toast('DNA loaded, but gene directions are not available yet: '+e.message);}}}
 loadedRefs();if(mode==='bundled')restoreConfiguration({...copy(B.default_config),references:copy(B.references)},false);renderAllowed();markStale();$('loadStatus').textContent=`Loaded ${S.refs.length} reference(s). Historical order files were not changed.`;
});
function updateAnnotationAvailability(){const r=S.refs.find(r=>r.id===$('annotationRef').value);const mapped=!!r?.mapping?.verified;$('annotationAvailability').textContent=mapped?'Public tracks can be loaded on this verified GRCh38 interval.':'This is a local/construct sequence. Import its actual annotations, or use genomic-coordinate input. Public tracks will not be pasted onto an unverified coordinate map.';$('fetchGenes').disabled=!mapped;$('discoverTracks').disabled=!mapped;$('loadTrack').disabled=!mapped||!S.tracks?.choices?.length;}
$('annotationRef').onchange=updateAnnotationAvailability;$('contextRef').onchange=contextMap;
$('fetchGenes').onclick=()=>task('fetchGenes','trackStatus',async()=>{let r=S.refs.find(r=>r.id===$('annotationRef').value);const d=await api('/api/genes',{reference:r});r.features=(r.features||[]).filter(f=>f.source!=='Ensembl GRCh38').concat(d.features);r.gene_annotation={...d,features:undefined};S.annotationProvenance.push({...d,features:undefined});$('trackStatus').textContent=`Loaded ${d.features.length} gene/transcript/exon features from Ensembl. Only strict RNA mode uses whole spans to reject candidates.`;renderTranscription();contextMap();markStale();});
$('discoverTracks').onclick=()=>task('discoverTracks','trackStatus',async()=>{S.tracks=await api('/api/h3k/discover',{cell:$('cellLine').value,fallback:$('allowFallback').checked});const t=S.tracks;$('trackChoice').innerHTML=t.choices.map((x,i)=>`<option value="${i}">${esc(x.cell)} | ${esc(x.accession)} | ${esc(x.output_type)} | ${esc(x.format)} | ${esc(x.experiment)}</option>`).join('')||'<option>No compatible files found</option>';$('trackStatus').textContent=t.status+'\n'+t.note+(t.is_proxy?'\nPROXY: using '+t.used_cell+' rather than the requested cell.':'');S.annotationProvenance.push(t);updateAnnotationAvailability();});
$('loadTrack').onclick=()=>task('loadTrack','trackStatus',async()=>{if(!S.tracks?.choices?.length)throw Error('Find a compatible track first.');const choice=S.tracks.choices[Number($('trackChoice').value)];let r=S.refs.find(r=>r.id===$('annotationRef').value);const d=await api('/api/h3k/load',{reference:r,choice,requested_cell:$('cellLine').value});r.features=(r.features||[]).filter(f=>f.type!=='H3K27ac').concat(d.features);r.h3k_signal=d.signal;r.h3k_metadata=d.metadata;S.annotationProvenance.push(d);$('trackStatus').textContent=d.status+': '+choice.accession+' / '+choice.cell+(d.metadata.is_proxy?' (PROXY)':'');contextMap();markStale();});
$('importBed').onclick=()=>task('importBed','trackStatus',async()=>{let r=S.refs.find(r=>r.id===$('annotationRef').value);const f=$('bedFile').files[0];if(!f)throw Error('Select a BED/bedGraph file.');const d=await api('/api/features/import',{reference:r,text:await f.text(),filename:f.name,kind:$('bedKind').value,local:$('bedLocal').checked,cell:$('cellLine').value});r.features=(r.features||[]).concat(d.features);S.annotationProvenance.push({source:f.name,kind:$('bedKind').value,cell:$('cellLine').value});$('trackStatus').textContent='Imported '+d.features.length+' interval(s).';contextMap();markStale();});
const RNA_TYPES=new Set(['gene','transcript','mrna','ncrna','lncrna','rrna','trna','misc_rna','pre_rna','primary_transcript','precursor_rna','cds','exon','utr',"5'utr","3'utr",'5utr','3utr','five_prime_utr','three_prime_utr']);
function geneFeatures(r){return (r.features||[]).filter(f=>RNA_TYPES.has(String(f.type||'').toLowerCase()));}
function declaration(r){if(!S.transcription[r.id])S.transcription[r.id]={mode:r.mapping?.verified?'annotations':'unknown',confirmed:false,intervals:[]};return S.transcription[r.id];}
function transcriptionIssues(){
 if(!$('rnaAvoidanceEnabled').checked)return [];
 const issues=[];
 if($('p_screen_rna')&&!$('p_screen_rna').checked)issues.push('RNA avoidance requires the complementary RNA search to remain enabled.');
 for(const r of S.refs){const d=declaration(r),mapped=!!r.mapping?.verified,fs=geneFeatures(r);
  if(d.mode==='unknown'){issues.push(r.id+': choose and confirm the gene/transcription direction in step 1.');continue;}
  if(!mapped&&!d.confirmed)issues.push(r.id+': confirm the construct transcription directions in step 1.');
  if(mapped&&d.mode!=='annotations'&&!d.confirmed)issues.push(r.id+': confirm the manual transcription declaration.');
  if(d.mode==='annotations'&&mapped&&!['FETCHED','USER_VERIFIED_COMPLETE'].includes(r.gene_annotation?.status))issues.push(r.id+': fetch gene/transcript annotations before design.');
  if(d.mode==='annotations'&&!mapped&&!fs.length)issues.push(r.id+': no gene/transcript features are loaded; specify directions instead.');
  if(d.mode==='none'&&fs.length)issues.push(r.id+': known gene annotations cannot be discarded by choosing no transcription.');
  if(d.mode==='none'&&!d.confirmed)issues.push(r.id+': confirm the no-known-transcription declaration.');
  if(d.mode==='manual'){
   if(!d.intervals.length)issues.push(r.id+': add at least one transcript interval.');
   for(const it of d.intervals)if(!Number.isInteger(it.start)||!Number.isInteger(it.end)||it.start<1||it.end>r.length||it.end<it.start||!['+','-','both'].includes(it.strand))issues.push(r.id+': check transcript bounds and direction.');
  }
 }
 return issues;
}
function rnaPolicyText(){const strict=$('rnaAvoidanceEnabled').checked,screen=$('p_screen_rna')?.checked!==false;return strict?'STRICT: whole-gene span exclusion and transcript-sequence screening. RNA conflicts reject candidates; the DNA-bound strand never changes.':screen?'TRANSCRIPT-SEQUENCE SCREENING: whole-gene exclusion OFF; public annotated-RNA search ON in full mode. Supplied RNA sequences are also checked. Short partial RNA hits are review evidence; long/high-identity hits still affect selection. This does not establish absence of RNA binding.':'PUBLIC RNA DATABASE SEARCH OFF: whole-gene exclusion is also off. Supplied RNA references, if present, are still checked locally. This is an explicit reduced-scope setting, not an RNA pass.';}
function showRNAReadiness(){const x=transcriptionIssues(),strict=$('rnaAvoidanceEnabled').checked,screen=$('p_screen_rna')?.checked!==false;$('rnaPolicyMode').value=screen?(strict?'strict':'transcripts'):'custom';$('rnaReadiness').classList.toggle('error',!!x.length);$('rnaReadiness').textContent=x.length?x.join('\n'):rnaPolicyText();}
$('rnaPolicyMode').onchange=()=>{const value=$('rnaPolicyMode').value;if(value==='custom')return;$('rnaAvoidanceEnabled').checked=value==='strict';if($('p_screen_rna'))$('p_screen_rna').checked=true;renderTranscription();$('strictDirectionDetails').open=value==='strict';markStale();};
function renderTranscription(){
 const modes=[['unknown','Choose the transcription direction...'],['annotations','Use annotated gene/transcript directions (review below)'],['forward','Forward (+): RNA follows the supplied sequence ->'],['reverse','Reverse (-): RNA follows its reverse complement <-'],['both','Both directions across the whole reference (usually excludes all)'],['manual','Multiple genes / known transcript boundaries: enter intervals'],['none','Explicitly no known transcription in this sequence']];
 $('transcriptionEditors').innerHTML=S.refs.map(r=>{const d=declaration(r),fs=geneFeatures(r),mapped=!!r.mapping?.verified;
  return `<section class="transcription-editor"><div class="row between"><h4>${esc(r.id)} <span class="hint">${fmt(r.length)} bp / ${mapped?'genomic interval':'local or synthetic DNA'}</span></h4><span class="pill ${d.mode==='unknown'?'warning':''}">${!$('rnaAvoidanceEnabled').checked?'NOT USED FOR SPAN EXCLUSION':mapped?esc(r.gene_annotation?.status||'GENES NOT LOADED'):d.confirmed?'DIRECTIONS CONFIRMED':'REVIEW REQUIRED'}</span></div>
  <label>Which direction does the gene / RNA run?<select data-transcription-ref="${esc(r.id)}" data-transcription-key="mode">${modes.map(([v,l])=>`<option value="${v}" ${v===d.mode?'selected':''}>${esc(l)}</option>`).join('')}</select></label>
  ${d.mode==='both'?'<div class="notice danger">Whole-sequence transcription on both strands leaves no RNA-avoiding DNA strand. For genes in separate intervals, use manual intervals instead.</div>':''}
  ${['forward','reverse'].includes(d.mode)?`<p class="hint">Declared ${d.mode} transcription: the <b>whole supplied sequence</b> is conservatively protected. Strict RNA compatibility requires ${d.mode==='forward'?'BOTTOM (-)':'TOP (+)'} DNA; if the trace-wide strand differs, strict mode rejects candidates. Use manual intervals to limit the modeled transcript.</p>`:''}
  ${d.mode==='manual'?`<div class="table-wrap"><table><thead><tr><th>Gene / transcript</th><th>Start (1-based)</th><th>End (inclusive)</th><th>RNA direction</th><th></th></tr></thead><tbody>${d.intervals.map((it,i)=>`<tr><td><input data-rna-ref="${esc(r.id)}" data-rna-index="${i}" data-rna-key="name" value="${esc(it.name)}"></td><td><input type="number" data-rna-ref="${esc(r.id)}" data-rna-index="${i}" data-rna-key="start" min="1" max="${r.length}" value="${it.start}"></td><td><input type="number" data-rna-ref="${esc(r.id)}" data-rna-index="${i}" data-rna-key="end" min="1" max="${r.length}" value="${it.end}"></td><td><select data-rna-ref="${esc(r.id)}" data-rna-index="${i}" data-rna-key="strand">${[['+','Forward (+) ->'],['-','<- Reverse (-)'],['both','Both directions']].map(([v,l])=>`<option value="${v}" ${it.strand===v?'selected':''}>${esc(l)}</option>`).join('')}</select></td><td><button data-rna-remove-ref="${esc(r.id)}" data-rna-remove-index="${i}">Remove</button></td></tr>`).join('')}</tbody></table></div><button data-add-transcript="${esc(r.id)}">Add gene/transcript interval</button><p class="hint">Start/end are local sequence coordinates. Known imported genes are also retained. Enter the full transcribed span including UTRs and introns. Unknown opposite-strand overlap is not treated as safe.</p>`:''}
  ${(!mapped||d.mode!=='annotations')?`<label class="check"><input type="checkbox" data-transcription-ref="${esc(r.id)}" data-transcription-key="confirmed" ${d.confirmed?'checked':''}> I reviewed this sequence's transcription directions and the stated interval assumptions.</label>`:'<p class="hint">Direction is read from the fetched annotation, not the H3K27ac cell line. Opposite-strand gene overlaps are exclusions only in strict mode.</p>'}
  <details><summary>Review loaded gene/transcript directions (${fs.length} feature parts)</summary>${fs.length?`<div class="table-wrap"><table><thead><tr><th>Name / type</th><th>Span in reference</th><th>Transcription</th><th>RNA-compatible strand (not an instruction to switch)</th></tr></thead><tbody>${fs.slice(0,150).map(f=>`<tr><td>${esc(f.label)} / ${esc(f.type)}</td><td>${fmt(f.start)}-${fmt(f.end)}</td><td>${f.strand==='+'?'Forward (+) ->':f.strand==='-'?'<- Reverse (-)':'UNKNOWN'}</td><td>${f.strand==='+'?'Bottom (-)':f.strand==='-'?'Top (+)':'Exclude both within overlap'}</td></tr>`).join('')}</tbody></table></div>`:'<p class="hint">No gene/transcript features loaded. This alone does not establish absence of transcription.</p>'}${fs.length>150?'<p class="hint">First 150 parts displayed; all remain in the saved annotations. Strict mode uses all for whole-span exclusions.</p>':''}</details></section>`;
 }).join('');showRNAReadiness();
}
$('rnaAvoidanceEnabled').onchange=()=>{if($('rnaAvoidanceEnabled').checked&&$('p_screen_rna'))$('p_screen_rna').checked=true;showRNAReadiness();markStale();};
document.addEventListener('change',e=>{
 const el=e.target;if(el.dataset.transcriptionRef){const d=S.transcription[el.dataset.transcriptionRef],k=el.dataset.transcriptionKey;d[k]=k==='confirmed'?el.checked:el.value;if(k==='mode')d.confirmed=false;renderTranscription();markStale();}
 if(el.dataset.rnaRef){const d=S.transcription[el.dataset.rnaRef],it=d.intervals[Number(el.dataset.rnaIndex)],k=el.dataset.rnaKey;it[k]=['start','end'].includes(k)?Number(el.value):el.value;d.confirmed=false;renderTranscription();markStale();}
 if(el.id==='p_screen_rna'){showRNAReadiness();markStale();}
});
document.addEventListener('click',e=>{const add=e.target.closest('[data-add-transcript]'),remove=e.target.closest('[data-rna-remove-ref]');if(add){const rid=add.dataset.addTranscript,r=S.refs.find(r=>r.id===rid);S.transcription[rid].intervals.push({name:'Transcript '+(S.transcription[rid].intervals.length+1),start:1,end:r.length,strand:'+'});S.transcription[rid].confirmed=false;renderTranscription();markStale();}if(remove){const d=S.transcription[remove.dataset.rnaRemoveRef];d.intervals.splice(Number(remove.dataset.rnaRemoveIndex),1);d.confirmed=false;renderTranscription();markStale();}});
function initializePresets(){
 $('parameterPreset').innerHTML=Object.entries(B.presets||{}).map(([id,p])=>`<option value="${esc(id)}">${esc(p.label)}</option>`).join('');
 $('parameterProvenance').innerHTML='<table><thead><tr><th>Setting</th><th>Starting value</th><th>Where it comes from</th><th>What to know</th></tr></thead><tbody>'+(B.parameter_notes||[]).map(x=>`<tr><td>${esc(x.setting)}</td><td>${esc(x.value)}</td><td>${esc(x.origin)}</td><td>${esc(x.explanation)}</td></tr>`).join('')+'</tbody></table>';
 $('applyPreset').onclick=()=>applyPreset($('parameterPreset').value,true);
 $('parameterPreset').onchange=()=>{const p=B.presets?.[$('parameterPreset').value];$('presetSummary').textContent=(p?.description||'')+' Click Apply to change the live controls.';};
 $('restoreErrfi1').onclick=()=>{restoreConfiguration({...copy(B.default_config),references:copy(B.references)},false);toast('Restored ERRFI1 / DP1 defaults.');};
 $('paperLayout').onclick=()=>{if(!confirm('Set all tracing-region widths to 10,000 bp and the primary-count target to 100, as in the paper example? Shorter sequences will not form a full region. Your existing result/order will not change.'))return;$('binWidth').value=10000;$('p_target_per_bin').value=100;$('p_selection_goal').value='target_count';updateAdaptiveControls();makeBins();markStale();toast('Applied the 10 kb / 100-primary planning example. Review regions in step 2.');};
 applyPreset('errfi1',false);
}
function lengthSummary(p){return p.length_mode==='adaptive'?`${p.arm_min}-${p.preferred_arm_max} nt preferred; extend through ${p.arm_max} nt when useful`:`${p.arm_min}-${p.arm_max} nt fixed search range`;}
function applyPreset(id,notify=true){const pr=B.presets?.[id];if(!pr)return;S.preset=id;$('parameterPreset').value=id;for(const [k,v] of Object.entries(pr.values)){const el=$(k==='avoid_gene_rna'?'rnaAvoidanceEnabled':k==='human_policy'?'humanPolicy':k==='allow_unresolved_bins'?'allowGrouped':'p_'+k);if(!el)continue;if(el.type==='checkbox')el.checked=v;else el.value=v;}$('presetSummary').innerHTML='<b>'+esc(pr.label)+'</b><br>'+esc(pr.description)+'<br><strong>'+esc(lengthSummary(pr.values))+'; '+pr.values.tm_min+'-'+pr.values.tm_max+' C; '+pr.values.minimum_gap+'-bp footprint gap.</strong>';updateAdaptiveControls();if(S.refs.length)renderTranscription();if(notify){markStale();toast('ERRFI1 QC defaults applied, including Top (+) and transcript-only RNA screening. Review the current regions and assignments.');}}
function updateAdaptiveControls(){const adaptive=$('p_length_mode').value==='adaptive';$('p_preferred_arm_max').disabled=!adaptive;$('p_extension_shortlist_per_bin').disabled=!adaptive;const p=getPolicy();$('adaptiveSummary').innerHTML='<b>'+esc(lengthSummary(p))+'.</b> On-target Tm: '+p.tm_min+'-'+p.tm_max+' C. '+(adaptive?'Short alternatives and longer rescue alternatives have separate QC budgets. Longer arms receive the same genome/RNA and individual-primary checks; retained QC thresholds are unchanged.':'Only the configured fixed range is searched.')+'<br><b>'+ (p.selection_goal==='maximize_sites'?'Keep adding compatible sites; the per-region target is a benchmark, not a cap.':'Stop each region at the selected primary-count target.')+'</b> Tracing regions remain the same size. Folding temperature and off-target thresholds do not follow the on-target Tm window.';const interaction=$('interactionPolicyStatus');if(interaction){interaction.textContent=interactionPolicyText(p);interaction.classList.toggle('danger',!p.screen_cross_primary);}}


function getBounds(){const o={};for(const r of S.refs)o[r.id]=[Number($('boundA_'+r.id)?.value||1),Number($('boundB_'+r.id)?.value||r.length)];return o;}
function makeBins(){try{const width=Number($('binWidth').value),bounds=getBounds();if(!Number.isInteger(width)||width<50||width>1000000)throw Error('Region width must be 50-1,000,000 bp.');let bins=[],tails=[];
 for(const r of S.refs){let [a,b]=bounds[r.id];if(!Number.isInteger(a)||!Number.isInteger(b)||a<1||b>r.length||b<a)throw Error('Design bounds are outside '+r.id);let count=Math.floor((b-a+1)/width),rem=(b-a+1)%width;if(!count)throw Error(r.id+' is shorter than one region.');if(rem&&$('tailPolicy').value==='strict')throw Error(r.id+': '+rem+' bp remain. Choose an exactly divisible interval or exclude the remainder.');
 for(let j=0;j<count;j++){let start=a+j*width,end=start+width-1,id=r.id+'_B'+String(j+1).padStart(3,'0'),m=r.mapping;let genomic_start=null,genomic_end=null;
 if(m){genomic_start=m.orientation==='-'?m.end-end+1:m.start+start-1;genomic_end=m.orientation==='-'?m.end-start+1:m.start+end-1;}
 bins.push({id,reference:r.id,number:j+1,start,end,length:width,center:(start+end)/2,imager:S.assignments[id]||'',genomic_start,genomic_end,chromosome:m?.chromosome||null});}
 if(rem)tails.push({reference:r.id,start:b-rem+1,end:b,length:rem});}
 for(const b of bins){if(!b.imager){const used=new Set(bins.map(x=>x.imager));let im=S.imagers.find(x=>!used.has(x.id));if(im)b.imager=im.id;}S.assignments[b.id]=b.imager;}
 S.bins=bins;S.tails=tails;$('binStatus').classList.remove('error');$('binStatus').textContent=bins.length+' adjacent '+fmt(width)+'-bp tracing regions. Individual oligo gaps are configured in step 4.';renderBins();renderImagers();renderBinMap();renderAllowed();
 }catch(e){$('binStatus').textContent=e.message;$('binStatus').classList.add('error');S.bins=[];renderBins();toast(e.message);}}
$('splitRegions').onclick=()=>{makeBins();markStale();};
function renderBins(){$('binRows').innerHTML=S.bins.map(b=>`<tr class="clickable" data-bin-preview="${esc(b.id)}"><td><b>${esc(b.id)}</b></td><td>${esc(b.reference)}</td><td>${fmt(b.start)}</td><td>${fmt(b.end)}</td><td>${fmt(b.length)} bp</td><td>${b.chromosome?esc(b.chromosome)+':'+fmt(b.genomic_start)+'-'+fmt(b.genomic_end):'Local sequence coordinates'}</td></tr>`).join('');$('tailWarning').classList.toggle('hidden',!S.tails.length);$('tailWarning').innerHTML='<b>Excluded trailing sequence (not shorter final regions):</b> '+S.tails.map(t=>esc(t.reference)+': '+fmt(t.start)+'-'+fmt(t.end)+' ('+fmt(t.length)+' bp)').join('; ');}
function renderImagers(){$('imagerRows').innerHTML=S.bins.map(b=>{let im=S.imagers.find(x=>x.id===S.assignments[b.id]);return `<tr><td><b>${esc(b.id)}</b></td><td><select data-imager-bin="${esc(b.id)}"><option value="">Choose an imager</option>${S.imagers.map(i=>`<option value="${esc(i.id)}" ${im?.id===i.id?'selected':''}>${esc(i.id)}</option>`).join('')}</select></td><td class="mono">${esc(im?.sequence||'not assigned')}</td><td class="mono">${esc(im?.docking||'not assigned')}</td><td>${esc(im?.dye||'')}</td><td>${im?`<button type="button" class="binding-preview-button" data-preview-imager="${esc(im.id)}" data-preview-region="${esc(b.id)}">Show binding</button>`:''}</td></tr>`;}).join('');renderImagerInventory();imagerWarning();renderBindingPreview();}
function imagerWarning(){const ids=S.bins.map(b=>S.assignments[b.id]),dups=[...new Set(ids.filter((v,i)=>v&&ids.indexOf(v)!==i))],missing=ids.filter(x=>!x).length;let s=[];if(missing)s.push(missing+' regions still need an imager; add new verified imager sequences or reduce the number of regions.');if(dups.length)s.push('Repeated imager(s): '+dups.join(', ')+'. Those regions cannot be independently decoded.');$('imagerWarning').classList.toggle('hidden',!s.length);$('imagerWarning').textContent=s.join(' ');}
document.addEventListener('change',e=>{if(e.target.dataset.imagerBin){const id=e.target.dataset.imagerBin;S.assignments[id]=e.target.value;S.bins.find(b=>b.id===id).imager=e.target.value;S.bindingImager=e.target.value;S.bindingRegion=id;renderImagers();markStale();}});
// An explicit convention is required for newly entered DNA. The legacy
// `sequence` field always means the imaging oligo; it never silently flips.
function resolveImagerEntry(entry){
 if(entry.sequence_orientation!=null&&entry.sequence_orientation!=='5to3')throw Error('Both sequences must be entered 5-prime to 3-prime.');
 const id=String(entry.id||'').trim().toUpperCase();
 if(!/^[A-Z0-9_-]+$/.test(id))throw Error('Use a name containing letters, numbers, underscores or dashes.');
 function read(key){const value=entry[key];if(value==null||value==='')return null;if(typeof value!=='string')throw Error(key+' must be DNA text.');const seq=value.replace(/\s/g,'').toUpperCase();if(!/^[ACGT]{8,40}$/.test(seq))throw Error('Enter 8-40 A/C/G/T bases written 5′ to 3′. Remove dye notation, direction labels and TT spacers; put dye details in Dye / note.');return seq;}
 function agree(a,b,label){if(a!=null&&b!=null&&a!==b)throw Error('Conflicting '+label+' fields. Check the sequence type and orientation.');return a??b;}
 const im=agree(read('sequence'),read('imager_sequence_5to3'),'imaging-oligo'),dock=agree(read('docking'),read('docking_sequence_5to3'),'docking-site');
 let entered=read('input_sequence_5to3'),kind=entry.input_kind;
 if(kind==null){if(entered!=null)throw Error('Choose whether these bases are the imaging oligo or its docking site.');kind=im!=null?'imager':dock!=null?'docking':null;}
 if(!['imager','docking'].includes(kind))throw Error('Choose which sequence you are entering: imaging oligo or docking site on the primary.');
 entered=entered??(kind==='imager'?im:dock);if(entered==null)throw Error('Enter the sequence for the selected input type.');
 const sequence=kind==='imager'?entered:rc(entered),docking=kind==='docking'?entered:rc(entered);
 agree(im,sequence,'imaging-oligo');agree(dock,docking,'docking-site / reverse-complement');
 return {id,sequence,docking,imager_sequence_5to3:sequence,docking_sequence_5to3:docking,input_kind:kind,input_sequence_5to3:entered,sequence_orientation:'5to3',dye:entry.dye||'not specified',source:entry.source||'user-defined'};
}
const inputKindLabel=kind=>kind==='docking'?'Binding-site / docking sequence':kind==='imager'?'Imaging oligo itself':'Not recorded in older result';
function imagerPairMarkup(im,showAlignment=true){return `<div class="sequence-pair"><div class="sequence-pair-row"><span class="sequence-kind">Imaging oligo itself <button class="help" data-help="imager">?</button></span><code data-sequence-role="imager">5′-${esc(im.sequence)}-3′</code><span class="hint">Separate fluorescent oligo; NOT appended to the primary.</span></div><div class="sequence-pair-row docking-row"><span class="sequence-kind">Binding-site / docking sequence <button class="help" data-help="docking">?</button></span><code data-sequence-role="docking">5′-${esc(im.docking)}-3′</code><span class="hint">Placed on BOTH ends of the primary, outside its TT spacers.</span></div></div>${showAlignment?ImagerBinding.render(im):''}`;}
function renderBindingPreview(id, region){
 if(id!==undefined){S.bindingImager=id;S.bindingRegion=region||null;}
 let im=S.imagers.find(i=>i.id===S.bindingImager)||S.imagers[0];
 if(!im){$('assignedBindingPreview').innerHTML='<p class="hint">Add an imager or docking site first.</p>';return;}
 S.bindingImager=im.id;
 $('bindingImagerChoice').innerHTML=S.imagers.map(i=>`<option value="${esc(i.id)}" ${i.id===im.id?'selected':''}>${esc(i.id)} - ${esc(i.dye)}</option>`).join('');
 const assigned=S.bins.filter(b=>S.assignments[b.id]===im.id).map(b=>b.id);
 $('bindingPreviewContext').textContent=(S.bindingRegion?'Viewing '+S.bindingRegion+'. ':assigned.length?'Assigned to '+assigned.join(', ')+'. ':'Not assigned to a region. ')+'Preview only; use the region dropdown above to change an assignment.';
 $('assignedBindingPreview').innerHTML=ImagerBinding.render(im,{showOverview:true})+`<details class="binding-canonical"><summary>Sequences to copy / order: both written 5′ to 3′</summary>${imagerPairMarkup(im,false)}</details>`;
}
$('bindingImagerChoice').onchange=()=>renderBindingPreview($('bindingImagerChoice').value);
document.addEventListener('click',e=>{
 const view=e.target.closest('[data-preview-imager]');
 if(view){renderBindingPreview(view.dataset.previewImager,view.dataset.previewRegion);$('assignedBindingPreview').scrollIntoView({behavior:'smooth',block:'center'});}
 const exp=e.target.closest('[data-binding-export]');
 if(exp){const figure=exp.closest('[data-binding-figure]'),kind=exp.dataset.bindingExport,svg=figure.querySelector('svg[data-binding-svg="'+kind+'"]');
   if(svg)saveBlob((figure.dataset.bindingImager||'imager').replace(/[^A-Za-z0-9_-]/g,'_')+'_'+kind+'.svg',new XMLSerializer().serializeToString(svg),'image/svg+xml');}
});
function describeBasePair(e){const col=e.target.closest?.('[data-base-index]');if(!col)return;const note=col.closest('[data-binding-figure]')?.querySelector('.binding-hover-help');if(note)note.textContent=col.dataset.pairDescription;}
document.addEventListener('pointerover',describeBasePair);
document.addEventListener('focusin',describeBasePair);
let imagerDraftVersion=0,imagerAdding=false;
function draftImager(){return resolveImagerEntry({id:$('newImagerName').value.trim()||'PREVIEW',input_kind:$('newImagerType').value,input_sequence_5to3:$('newImagerSequence').value,dye:$('newImagerDye').value.trim()||'not specified',source:'user-defined'});}
function updateImagerPreview(reset=true){
 if(reset){imagerDraftVersion++;$('newImagerConfirmed').checked=false;}
 const kind=$('newImagerType').value;
 $('newImagerSequenceLabel').textContent=kind==='imager'?"Imaging oligo itself, 5′ to 3′":kind==='docking'?"Binding-site / docking sequence ON THE PRIMARY, 5′ to 3′":"Entered DNA sequence, 5′ to 3′ (choose its type above)";
 $('newImagerTypeHelp').textContent=kind==='docking'?'These entered bases will be used directly as the primary docking site. The imaging oligo is calculated as their reverse complement. This is NOT the genomic target DNA.':kind==='imager'?'These entered bases are the fluorescent imaging oligo. Its reverse complement will be used as the docking site on each primary.':'This choice changes how the entered bases are interpreted. A DNA sequence alone cannot reveal which molecule it represents.';
 $('newImagerError').textContent='';$('newImagerConfirmed').disabled=true;$('addImager').disabled=true;
 if(!kind||!$('newImagerSequence').value.trim()){$('newImagerPreview').innerHTML='<p class="hint">Choose a sequence type and enter its bases to preview both molecules.</p>';return;}
 let validPair=false;try{const entry=draftImager();validPair=true;$('newImagerPreview').innerHTML=imagerPairMarkup(entry)+`<p class="hint"><b>Entered as:</b> ${esc(inputKindLabel(entry.input_kind))}. Both sequences above are written 5′ to 3′.</p>`;
   const conflict=S.imagers.find(i=>i.id===entry.id||i.sequence===entry.sequence);
   if(conflict)throw Error('This name or resolved sequence pair already exists as '+conflict.id+'. Select that existing imager instead of adding a second copy.');
   if(!$('newImagerName').value.trim()){ $('newImagerError').textContent='Add a unique name before saving this pair.';return; }
   $('newImagerConfirmed').disabled=false;$('addImager').disabled=imagerAdding||!$('newImagerConfirmed').checked;
 }catch(e){if(!validPair)$('newImagerPreview').innerHTML='<p class="hint">No valid sequence pair to display. Correct the entry above.</p>';$('newImagerError').textContent=e.message;}
}
for(const id of ['newImagerName','newImagerSequence','newImagerDye'])$(id).addEventListener('input',()=>updateImagerPreview());
$('newImagerType').addEventListener('change',()=>updateImagerPreview());
$('newImagerConfirmed').addEventListener('change',()=>updateImagerPreview(false));
function renderImagerInventory(){$('imagerInventoryRows').innerHTML=S.imagers.map(i=>`<tr><td><b>${esc(i.id)}</b></td><td class="mono">${esc(i.sequence)}</td><td class="mono">${esc(i.docking)}</td><td>${esc(inputKindLabel(i.input_kind))}</td><td>${esc(i.dye)}</td><td><button type="button" class="binding-preview-button" data-preview-imager="${esc(i.id)}">Show binding</button></td></tr>`).join('');}
$('addImager').onclick=async()=>{const version=imagerDraftVersion;try{
 if(!$('newImagerConfirmed').checked||$('newImagerConfirmed').disabled)throw Error('Check the two labeled sequence previews before adding this pair.');
 const entry=draftImager();if(!$('newImagerName').value.trim())throw Error('Enter a unique imager name.');
 if(S.imagers.some(i=>i.id===entry.id||i.sequence===entry.sequence))throw Error('That name or resolved imaging/docking pair already exists.');
 imagerAdding=true;$('addImager').disabled=true;
 const candidates=[...copy(S.imagers),entry];
 const resolved=B.offline?candidates:(await api('/api/imagers/resolve',{imagers:candidates})).imagers;
 if(version!==imagerDraftVersion)throw Error('The entry changed during validation. Check the new preview and add again.');
 S.imagers=resolved;renderImagers();markStale();$('newImagerName').value='';$('newImagerSequence').value='';updateImagerPreview();
 toast('Added '+entry.id+' from '+inputKindLabel(entry.input_kind).toLowerCase()+'. Assign it to a region. The primary uses '+entry.docking+'.');
 }catch(e){$('newImagerError').textContent=e.message;toast(e.message);}finally{imagerAdding=false;$('addImager').disabled=!$('newImagerConfirmed').checked||$('newImagerConfirmed').disabled;}
};
$('exportImagerInventory').onclick=()=>saveBlob('Imager_and_docking_inventory.csv',csv([['Imager_ID','Imager_Sequence_5to3','Docking_Sequence_in_Primary_5to3','Entered_Type','Entered_Sequence_5to3','Dye_or_note'],...S.imagers.map(i=>[i.id,i.sequence,i.docking,i.input_kind,i.input_sequence_5to3,i.dye])]),'text/csv');

const POLICY_GROUPS=[
 ['Analysis mode and target counts','Full mode uses installed BLAST and full DNA folding before selecting. A local preview is clearly labeled and cannot establish human-genome specificity.',[
  ['mode','Mode','select','full',[['full','Full QC-guided selection'],['local_preview','LOCAL-ONLY preview; missing genome and full folding']],'unknown'],
  ['selection_goal','How many primaries should be selected?','select','maximize_sites',[['maximize_sites','Maximize compatible sites (no count-target cap)'],['target_count','Stop at target primary count per region']],'selection_goal'],
  ['target_per_bin','Primary-count benchmark / cap in target mode','number',100,null,'shortfall'],['minimum_per_bin','Minimum count for a region to be considered filled','number',50,null,'shortfall'],
  ['threads','BLAST threads (maximum 2)','number',2,null,'offtarget']]],
 ['Binding-arm thermodynamics and sequence','Uses the same temperature profile to screen candidates and report the selected set. Thresholds apply to the arm, not the complete primary.',[
  ['profile','Temperature profile','select','oligominer1_2nM',[['oligominer1_2nM','OligoMiner-style: 2 nM / trace target'],['oligominer1_defaults','OligoMiner source defaults: 25 / 25 nM'],['previous_NN_audit','Earlier NN4 audit: 2 nM / trace target']],'profile'],
  ['tm_min','Minimum arm Tm (C)','number',45,null,'tm'],['tm_max','Maximum arm Tm (C)','number',55,null,'tm'],['tm_strict','Enforce the Tm window (strict bounds)','check',true,null,'tm'],
  ['gc_min','Minimum arm GC (%) (0 = no separate lower cutoff)','number',0,null,'gc'],['gc_max','Maximum arm GC (%)','number',76,null,'gc'],
  ['length_mode','Binding-arm length strategy','select','adaptive',[['adaptive','Adaptive: short preferred, longer alternatives allowed'],['fixed','Fixed allowed range only']],'adaptive_length'],
  ['arm_min','Minimum arm length (nt)','number',36,null,'adaptive_length'],
  ['preferred_arm_max','Preferred short-arm maximum (nt)','number',42,null,'adaptive_length'],
  ['arm_max','Maximum allowed arm length (extension ceiling, nt)','number',90,null,'adaptive_length'],
  ['bound_strand','ONE DNA-bound strand for ALL tracing regions','select','+',[['unselected','Choose a fixed strand before running'],['+','Top (+) for ALL regions: arm is the reverse complement'],['-','Bottom (-) for ALL regions: arm equals the forward sequence']],'single_strand'],
  ['homopolymer_max','Longest permitted identical-base run','number',5,null,'gc']]],
 ['Candidate search and spatial layout','Sampling settings control workload. Counts of generated alternatives and those receiving detailed QC are reported. All exact repeat copies are considered.',[
  ['candidate_step','Candidate start stride (bp)','number',1,null,'shortlist'],['variants_per_start','Alternatives per start, per length band','number',2,null,'shortlist'],
  ['shortlist_per_bin','Preferred / fixed-range QC shortlist per region','number',400,null,'shortlist'],['extension_shortlist_per_bin','Additional longer-arm QC candidates per region','number',200,null,'adaptive_length'],['max_total_candidates','Maximum detailed-QC candidates per run','number',5000,null,'shortlist'],
  ['minimum_gap','Minimum primary footprint gap (bp)','number',5,null,'footprint'],['footprint','Spacing footprint','select','extended',[['extended','Arm plus contiguous spacer/handle matching'],['arm','Arm only (less strict; extensions possible)']],'footprint'],
  ['repeat_policy','Repeat policy','select','within_bin',[['within_bin','Allow compatible repeats inside the same region'],['unique','Require exactly one local full-arm site']],'repeats']]],
 ['Off-target evidence used during selection','High predicted stability at multiple unintended sites is penalized. Unsupported mismatch estimates remain unknown; they are not silently treated as weak binding.',[
  ['high_offtarget_tm','High off-target Tm threshold (C)','number',42,null,'high_hit'],['max_high_offtarget_sites','Maximum high-Tm unintended site clusters','number',0,null,'hit_count'],
  ['moderate_offtarget_tm','Moderate off-target Tm threshold (C)','number',35,null,'offtarget_tm'],['max_moderate_offtarget_sites','Maximum moderate/high-Tm site clusters','number',10,null,'hit_count'],
  ['reject_unknown_long','Reject long, strong DNA alignments with unknown Tm','check',true,null,'unknown'],
  ['screen_rna','Search annotated human RNA (keep ON in full mode)','check',true,null,'rna_screen'],['rna_exact_min','Exact RNA tract threshold in STRICT mode only (nt)','number',18,null,'rna_screen'],['max_strong_rna_sites','Maximum strong RNA-risk clusters (0 = reject any)','number',0,null,'rna_screen']]],
 ['Individual-primary folding and interaction screening','Individual-primary folding/accessibility and self-dimer checks remain active in full mode. Cross-primary screening is a separate explicit setting; OFF skips its exact-stem searches and pairwise cofolding. Aqueous cofold energy is not an experimental dimer fraction.',[
  ['screen_cross_primary','Probe-probe interaction screening','check',false,null,'cofold'],
  ['minimum_arm_unpaired','Minimum mean arm unpaired probability (0-1)','number',0.4,null,'unpaired'],
  ['minimum_cofold_gain','Self-dimer / enabled cross-primary gain limit (kcal/mol)','number',-12,null,'cofold'],
  ['maximum_pair_stem','Self-dimer / enabled cross-primary stem limit (bp)','number',11,null,'exact_stem'],
  ['folding_temperature','Aqueous folding temperature (C)','number',42,null,'folding'],['sodium_M','Aqueous folding sodium (M)','number',0.390,null,'folding']]]
];
for(const [,,fields] of POLICY_GROUPS)for(const field of fields){if(Object.hasOwn(B.presets?.errfi1?.values||{},field[0]))field[3]=B.presets.errfi1.values[field[0]];}
$('policyFields').innerHTML=POLICY_GROUPS.map(([title,description,fields])=>`<article class="card"><h3>${title}</h3><p class="hint">${description}</p><div class="grid4">${fields.map(([key,label,type,value,options,term])=>`<div><label ${type==='check'?'class="check"':''}>${type==='check'?`<input id="p_${key}" type="checkbox" ${value?'checked':''}>`:''}${label} <button class="help" data-help="${term}">?</button>${type==='select'?`<select id="p_${key}">${options.map(([v,l])=>`<option value="${v}" ${v===value?'selected':''}>${l}</option>`).join('')}</select>`:type==='number'?`<input id="p_${key}" type="number" value="${value}" step="any">`:''}</label><div class="setting-description">${esc(B.glossary?.[term]?.label||'')}</div></div>`).join('')}</div></article>`).join('');
document.addEventListener('input',e=>{if(e.target.id?.startsWith('p_')){updateAdaptiveControls();markStale();}});
document.addEventListener('change',e=>{if(e.target.id?.startsWith('p_')){updateAdaptiveControls();markStale();}});
function interactionPolicyText(p){return 'Probe-probe interaction screening: '+(p.screen_cross_primary?'ON — within-region and between-region cross-primary checks requested.':'OFF — within-region and between-region: NOT RUN - USER DISABLED. No cross-primary stem searches or thermodynamic pairwise cofolds are requested.')+' Individual-primary folding/accessibility and self-dimer checks remain requested in full mode. Self-dimer thresholds still apply. Order files remain DRAFT.';}
function getPolicy(){const p={};for(const [, ,fields] of POLICY_GROUPS)for(const [key,,type] of fields){p[key]=type==='check'?$('p_'+key).checked:type==='number'?Number($('p_'+key).value):$('p_'+key).value;}p.threads=Math.max(1,Math.min(2,Math.floor(p.threads)||2));p.avoid_gene_rna=$('rnaAvoidanceEnabled').checked;p.preset_id=S.preset;p.human_policy=$('humanPolicy').value;p.allow_unresolved_bins=$('allowGrouped').checked;return p;}
$('humanPolicy').onchange=()=>{$('humanOverrideWarning').classList.toggle('hidden',$('humanPolicy').value!=='report_only');markStale();};
function renderAllowed(){$('allowRows').innerHTML=S.allowed.map((r,i)=>`<tr><td><select data-allow="${i}" data-key="bin_id"><option value="*">All regions (use carefully)</option>${S.bins.map(b=>`<option ${r.bin_id===b.id?'selected':''}>${esc(b.id)}</option>`).join('')}</select></td><td><select data-allow="${i}" data-key="source"><option value="human" ${r.source==='human'?'selected':''}>Human DNA</option><option value="construct" ${r.source==='construct'?'selected':''}>Construct</option></select></td><td><input data-allow="${i}" data-key="reference" value="${esc(r.reference)}" placeholder="chr7"></td><td><input data-allow="${i}" data-key="start" value="${r.start}" type="number"></td><td><input data-allow="${i}" data-key="end" value="${r.end}" type="number"></td><td><select data-allow="${i}" data-key="strand"><option value="-" ${r.strand==='-'?'selected':''}>Bottom (-)</option><option value="+" ${r.strand==='+'?'selected':''}>Top (+)</option></select></td><td><button data-remove-allow="${i}">Remove</button></td></tr>`).join('');}
$('addAllowed').onclick=()=>{S.allowed.push({bin_id:S.bins[0]?.id||'*',source:'human',reference:'',start:1,end:1,strand:'-'});renderAllowed();markStale();};
document.addEventListener('change',e=>{if(e.target.dataset.allow!=null){const r=S.allowed[Number(e.target.dataset.allow)],k=e.target.dataset.key;r[k]=['start','end'].includes(k)?Number(e.target.value):e.target.value;markStale();}});
document.addEventListener('click',e=>{const b=e.target.closest('[data-remove-allow]');if(b){S.allowed.splice(Number(b.dataset.removeAllow),1);renderAllowed();markStale();}});
function showEnvironment(){let e=B.environment||{};$('environment').innerHTML=[['Python',e.python_version||'offline'],['BLAST',e.blastn?'found':'NOT AVAILABLE'],['DNA folding',e.vienna?'found':'NOT AVAILABLE'],['Public databases',Object.keys(e.databases||{}).length?'found':'NOT LOADED'],['bigWig reader (optional)',e.pybigwig?'found':'BED peaks still supported']].map(([a,b])=>`<span class="pill ${b.includes('NOT')?'warning':''}">${esc(a)}: ${esc(b)}</span>`).join('');$('humanDB').value=e.databases?.human||'';$('rnaDB').value=e.databases?.rna_human||'';}
// Readable run names use only the regions actually being designed.
function proposedRunName(){
 const custom=$('runName').value.trim().replace(/\s+/g,' ').slice(0,180);if(custom)return custom;
 const imagers=[...new Set(S.bins.map(b=>S.assignments[b.id]||b.imager).filter(Boolean))];
 const locations=[...new Set(S.bins.map(b=>b.chromosome&&b.genomic_start!=null&&b.genomic_end!=null?`${b.chromosome}:${b.genomic_start}-${b.genomic_end}`:`${b.reference}:${b.start}-${b.end} (local)`))];
 if(!locations.length)return 'RASERFISH design';
 return (imagers.slice(0,4).join(', ')+(imagers.length>4?' +'+(imagers.length-4):'')+(imagers.length?' · ':'')+locations.slice(0,2).join('; ')+(locations.length>2?'; +'+(locations.length-2)+' regions':'')).slice(0,180);
}
// End readable run names.
function configuration(){const c={title:proposedRunName(),run_name:$('runName').value.trim(),output_detail:$('outputDetail').value,references:copy(S.refs),imagers:copy(S.imagers),assignments:copy(S.assignments),allowed_intervals:copy(S.allowed),bin_width:Number($('binWidth').value),tail_policy:$('tailPolicy').value,bounds:getBounds(),policy:getPolicy(),transcription:copy(S.transcription),databases:{human:$('humanDB').value.trim(),rna_human:$('rnaDB').value.trim()},background_text:$('constructFasta').value.trim(),transcript_text:$('transcriptFasta').value.trim(),annotation_provenance:copy(S.annotationProvenance)};for(const [key,textKey] of [['background','background_text'],['transcripts','transcript_text']]){if(S.importedInputs?.[key]&&c[textKey]===S.importedInputs[key].text){c[key]=copy(S.importedInputs[key].records);delete c[textKey];}}return c;}
function renderRunSummary(){$('runName').placeholder=proposedRunName();const p=getPolicy();$('runSummary').innerHTML=`<div class="notice ${['+','-'].includes(p.bound_strand)?'':'danger'}"><b>RASER-FISH strand lock: ${p.bound_strand==='+'?'TOP (+) for ALL regions':p.bound_strand==='-'?'BOTTOM (-) for ALL regions':'NOT SELECTED - choose in Step 4'}</b><br>The chosen RNA policy can reject candidates, but never flips a primary to the opposite strand. Use this same strand in every separate run belonging to the trace.</div><div class="kpis"><div class="kpi"><strong>${S.refs.length}</strong><span>reference sequences</span></div><div class="kpi"><strong>${S.bins.length}</strong><span>equal regions</span></div><div class="kpi"><strong>${fmt($('binWidth').value)}</strong><span>bp per region</span></div><div class="kpi"><strong>${p.target_per_bin}</strong><span>${p.selection_goal==='maximize_sites'?'benchmark (not a cap)':'target primaries / region'}</span></div><div class="kpi"><strong>${p.shortlist_per_bin+(p.length_mode==='adaptive'?p.extension_shortlist_per_bin:0)}</strong><span>maximum QC candidates / region</span></div></div><p class="hint">${esc(lengthSummary(p))}; ${p.tm_min}-${p.tm_max} C. ${p.selection_goal==='maximize_sites'?'Compatible-site count search':'Target-count capped search'}; retained individual-primary and specificity QC gates stay active.</p><div class="notice ${p.screen_cross_primary?'':'danger'}">${esc(interactionPolicyText(p))}</div><div class="notice ${p.mode==='local_preview'?'danger':''}"><b>${p.mode==='full'?'FULL QC-GUIDED SELECTION':'LOCAL-ONLY PREVIEW'}</b><br>${p.mode==='full'?'Requires real DNA folding and the selected local BLAST databases. No software fallback is allowed.':'Genome/RNA searches and full thermodynamic folding will NOT run. A preview is not a genome-screened probe order.'}</div>`;}
function saveBlob(name,text,type='text/plain'){const a=document.createElement('a');a.href=URL.createObjectURL(new Blob([text],{type}));a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(a.href),3000);}
function restoreConfiguration(c,keepDatabases=true){
 if(!c?.references?.length)throw Error('Settings need their sequence inputs. Import the complete settings ZIP.');
 S.refs=copy(c.references);S.imagers=copy(c.imagers).map(resolveImagerEntry);S.assignments=copy(c.assignments||{});S.allowed=copy(c.allowed_intervals||[]);S.transcription=copy(c.transcription||{});S.annotationProvenance=copy(c.annotation_provenance||[]);S.title=c.title;$('runName').value=c.run_name||'';S.bindingImager=Object.values(S.assignments)[0]||S.imagers[0]?.id;S.bindingRegion=Object.keys(S.assignments)[0]||null;
 $('binWidth').value=c.bin_width??20000;$('tailPolicy').value=c.tail_policy||'strict';$('outputDetail').value=c.output_detail||'standard';
 S.preset=c.policy?.preset_id||'custom';
 $('p_screen_cross_primary').checked=c.policy?.screen_cross_primary===true;
 for(const [k,v] of Object.entries(c.policy||{})){const el=$(k==='avoid_gene_rna'?'rnaAvoidanceEnabled':k==='human_policy'?'humanPolicy':k==='allow_unresolved_bins'?'allowGrouped':'p_'+k);if(el){if(el.type==='checkbox')el.checked=v;else el.value=v;}}
 $('p_threads').value=Math.max(1,Math.min(2,Math.floor(Number(c.policy?.threads))||2));
 if(Number(c.policy?.threads)>2)toast('Imported BLAST threads reduced to 2 for the v1.6.1 memory limit.');
 $('parameterPreset').value=B.presets?.[S.preset]?S.preset:'errfi1';
 $('presetSummary').textContent=S.preset==='errfi1'?B.presets.errfi1.description:'Imported settings preserved; applying ERRFI1 will reset the QC controls.';
 loadedRefs();for(const r of S.refs){const bounds=c.bounds?.[r.id];if(bounds){$('boundA_'+r.id).value=bounds[0];$('boundB_'+r.id).value=bounds[1];}}
 S.assignments=copy(c.assignments||{});makeBins();renderAllowed();updateAdaptiveControls();
 const fasta=rows=>(rows||[]).map(r=>'>'+(r.description||r.id)+'\n'+r.sequence).join('\n');
 $('constructFasta').value=c.background_text||fasta(c.background);$('transcriptFasta').value=c.transcript_text||fasta(c.transcripts);S.importedInputs={background:Array.isArray(c.background)&&!c.background_text?{records:copy(c.background),text:$('constructFasta').value.trim()}:null,transcripts:Array.isArray(c.transcripts)&&!c.transcript_text?{records:copy(c.transcripts),text:$('transcriptFasta').value.trim()}:null};
 if(keepDatabases){$('humanDB').value=c.databases?.human||B.environment?.databases?.human||'';$('rnaDB').value=c.databases?.rna_human||B.environment?.databases?.rna_human||'';}
 $('humanOverrideWarning').classList.toggle('hidden',$('humanPolicy').value!=='report_only');markStale();
}
$('saveConfig').onclick=async()=>{try{
 if(B.offline)throw Error('Open the local Studio to save a portable settings ZIP.');
 const response=await fetch('/api/settings/export',{method:'POST',headers:{'Content-Type':'application/json','X-Studio-Token':B.token},body:JSON.stringify({config:configuration()})});
 if(!response.ok)throw Error((await response.json()).error||'Settings export failed.');
 saveBlob('RASERFISH_design_settings.zip',await response.blob(),'application/zip');
}catch(e){toast(e.message);}};
$('loadConfig').onchange=async()=>{const f=$('loadConfig').files[0];if(!f)return;try{let c;if(f.name.toLowerCase().endsWith('.json')){c=JSON.parse(await f.text());if(!c.references?.[0]?.sequence)throw Error('Import the settings ZIP with its inputs, not the manifest alone.');}else{const bytes=new Uint8Array(await f.arrayBuffer());let raw='';for(let i=0;i<bytes.length;i+=8192)raw+=String.fromCharCode(...bytes.slice(i,i+8192));c=(await api('/api/settings/import',{base64:btoa(raw)})).config;}restoreConfiguration(c);toast('Imported settings and inputs. Review database paths before running.');}catch(e){toast(e.message);}finally{$('loadConfig').value='';}};
// Completed-run browsing is read-only and never starts or resumes a design.
let completedRuns=[],completedRunsLoading=false,completedRunOpening=false;
function selectedCompletedRun(){return completedRuns.find(r=>r.id===$('completedRunChoice').value);}
function updateCompletedRunControls(){
 const busy=B.offline||completedRunsLoading||completedRunOpening;
 $('refreshCompletedRuns').disabled=!!busy;
 $('completedRunChoice').disabled=!!(busy||!completedRuns.length);
 $('openCompletedRun').disabled=!!(busy||S.designActive||!selectedCompletedRun());
}
async function refreshCompletedRuns(){
 if(B.offline||completedRunsLoading||completedRunOpening)return;
 completedRunsLoading=true;$('completedRunsStatus').textContent='Reading completed run names…';updateCompletedRunControls();
 try{
  const response=await fetch('/api/runs',{cache:'no-store'});const data=await response.json();
  if(!response.ok)throw Error(data.error||(response.status===404?'Restart the updated Studio to browse completed runs.':'Could not read completed runs.'));
  if(!Array.isArray(data.runs))throw Error('The completed-run list is unreadable.');
  const previous=$('completedRunChoice').value;
  completedRuns=data.runs.filter(r=>r&&typeof r.id==='string'&&r.status==='complete');
  $('completedRunChoice').innerHTML=completedRuns.map(r=>`<option value="${esc(r.id)}">${esc(r.name||r.id)}${Number.isSafeInteger(r.summary?.counts?.selected)?' — '+r.summary.counts.selected+' primaries':''} · ${esc(r.id)}</option>`).join('');
  $('completedRunChoice').value=completedRuns.some(r=>r.id===previous)?previous:(completedRuns[0]?.id||'');
  $('completedRunsStatus').textContent=S.designActive?'A design is active. Saved results can be opened after it stops.':completedRuns.length?completedRuns.length+' completed runs available.':'No completed runs are available.';
 }catch(e){completedRuns=[];$('completedRunChoice').innerHTML='';$('completedRunsStatus').textContent=e.message||'Could not connect to the local program.';}
 finally{completedRunsLoading=false;updateCompletedRunControls();}
}
async function openCompletedRun(){
 const run=selectedCompletedRun();if(!run||B.offline||completedRunsLoading||completedRunOpening||S.designActive)return;
 completedRunOpening=true;updateCompletedRunControls();$('completedRunsStatus').textContent='Opening saved result…';
 try{
  const response=await fetch('/api/result?id='+encodeURIComponent(run.id),{cache:'no-store'});const result=await response.json();
  if(!response.ok)throw Error(result.error||'Saved result is not available.');
  if(result.analysis_complete===false||/INCOMPLETE|STOPPED|FAILED|CANCELLED/.test(result.status||''))throw Error('This analysis is incomplete; no completed result was opened.');
  if(!Array.isArray(result.candidates)||!Array.isArray(result.bins))throw Error('Saved result data are incomplete.');
  S.result=result;S.job=run.id;S.stale=false;S.viewRef=null;initResult();renderResults();
  $('completedRunsStatus').textContent='Opened '+(run.name||run.id)+'. Order status: DRAFT.';
 }catch(e){$('completedRunsStatus').textContent=e.message;toast(e.message);}
 finally{completedRunOpening=false;updateCompletedRunControls();}
}
function initializeCompletedRuns(){
 $('completedRunsPanel').classList.toggle('hidden',!!B.offline);
 $('refreshCompletedRuns').onclick=refreshCompletedRuns;$('openCompletedRun').onclick=openCompletedRun;
 $('completedRunChoice').onchange=updateCompletedRunControls;updateCompletedRunControls();
}
// End completed-run browsing.
// Saved checkpoints are listed on reload; automatic retries are controlled by the server.
let checkpointRuns=[],checkpointLoading=false,checkpointBusy=false,checkpointActiveElsewhere=false,checkpointResumeUncertain=false,activeDesignLookup=null;
function selectedCheckpoint(){return checkpointRuns.find(r=>r.id===$('checkpointChoice').value);}
function updateCheckpointControls(){
 const blocked=B.offline||checkpointLoading||checkpointBusy||S.designActive||checkpointActiveElsewhere||checkpointResumeUncertain;
 $('refreshCheckpoints').disabled=!!(B.offline||checkpointLoading||checkpointBusy);
 $('checkpointChoice').disabled=!!(blocked||!checkpointRuns.length);
 $('resumeCheckpoint').disabled=!!(blocked||!selectedCheckpoint());
 if(S.designActive||checkpointActiveElsewhere)$('runDesign').disabled=true;
}
function renderCheckpointDetail(){
 const r=selectedCheckpoint();
 const count=n=>Number.isSafeInteger(n)&&n>=0?fmt(n):'not recorded';
 $('checkpointDetail').textContent=r?`${r.summary||r.id}\nSaved run: ${r.id}${r.created_utc?' · '+r.created_utc:''}\nReference search batches saved: ${count(r.completed_batches)} / ${count(r.total_batches)}. Individual-primary checks saved: ${count(r.folded_candidates)} / ${count(r.total_candidates)}.\nINCOMPLETE — saved evidence will be verified before reuse. This is not a completed QC result.${r.auto_recovery?'\nAutomatic recovery: '+r.auto_recovery.state+(r.auto_recovery.reason?' — '+r.auto_recovery.reason:''):''}`:'';
 updateCheckpointControls();
}
async function refreshCheckpoints(){
 if(B.offline||checkpointLoading||checkpointBusy)return;
 checkpointLoading=true;$('checkpointStatus').textContent='Checking for saved progress…';updateCheckpointControls();
 try{
  const response=await fetch('/api/checkpoints',{cache:'no-store'});
  if(!response.ok){let message='Cannot list saved progress.';try{message=(await response.json()).error||message;}catch{}if(response.status===404)message='Saved-progress controls are unavailable in this running program. Restart the updated Studio to use them.';throw Error(message);}
  const data=await response.json();if(!Array.isArray(data.runs))throw Error('The saved-progress list is unreadable. Refresh to try again.');
  const previous=$('checkpointChoice').value;
  checkpointRuns=data.runs.filter(r=>r&&typeof r.id==='string'&&r.status==='incomplete');
  checkpointActiveElsewhere=!!data.active_design||!!data.active_run;checkpointResumeUncertain=false;
  $('checkpointChoice').innerHTML=checkpointRuns.map(r=>`<option value="${esc(r.id)}">${esc(typeof r.summary==='string'?r.summary:r.id)} — ${esc(r.id)}</option>`).join('');
  $('checkpointChoice').value=checkpointRuns.some(r=>r.id===previous)?previous:(checkpointRuns[0]?.id||'');
  $('checkpointStatus').textContent=checkpointActiveElsewhere?'A design is already active or waiting for an automatic retry. Its saved progress is managed by the local program.':checkpointRuns.length?'These runs are incomplete. Automatic recovery is managed by the local program; Resume saved progress remains available when recovery has stopped.':'No resumable saved checkpoints are available. Older runs may not contain the checkpoints needed for a safe resume.';
  renderCheckpointDetail();
  if(checkpointActiveElsewhere&&!S.job)await refreshActiveDesign(true);
 }catch(e){checkpointRuns=[];$('checkpointChoice').innerHTML='';$('checkpointDetail').textContent='';$('checkpointStatus').textContent=e.message||'Cannot connect to the local program. Saved progress has not been checked.';}
 finally{checkpointLoading=false;updateCheckpointControls();}
}
// Read-only reconnect: the browser never launches an automatic retry.
async function refreshActiveDesign(attach=false){
 if(B.offline)return null;
 if(!activeDesignLookup)activeDesignLookup=(async()=>{
  const response=await fetch('/api/active-design',{cache:'no-store'});
  if(!response.ok){if(response.status===404)return null;let message='Cannot reconnect to the active design.';try{message=(await response.json()).error||message;}catch{}throw Error(message);}
  const data=await response.json();if(typeof data.token==='string'&&data.token)B.token=data.token;
  return data;
 })();
 let data;try{data=await activeDesignLookup;}finally{activeDesignLookup=null;}
 if(attach&&data&&typeof data.id==='string'&&data.id&&!S.job){
  S.job=data.id;S.designActive=true;pollFailures=0;$('runDesign').disabled=true;$('cancelRun').classList.remove('hidden');
  $('runState').textContent='Reconnected to saved design';
  try{restoreConfiguration(data.config);renderRunSummary();}catch(e){toast('The active design is being monitored, but its settings could not be displayed: '+e.message);}
  updateCheckpointControls();go(5);pollJob();
 }
 return data;
}
async function resumeCheckpoint(){
 const run=selectedCheckpoint();if(!run||B.offline||checkpointLoading||checkpointBusy||S.designActive||checkpointActiveElsewhere||checkpointResumeUncertain)return;
 checkpointBusy=true;const wasDisabled=$('runDesign').disabled;$('runDesign').disabled=true;updateCheckpointControls();
 $('checkpointStatus').textContent='Verifying saved inputs and checkpoints…';
 try{
  const d=await api('/api/resume',{id:run.id});
  if(!d||typeof d.id!=='string'||!d.id)throw Error('The local program did not return a run identifier. Refresh before trying again.');
  S.job=d.id;S.designActive=true;pollFailures=0;
  $('cancelRun').classList.remove('hidden');$('runState').textContent='Resuming saved progress';
  $('progressStage').textContent='Continuing verified checkpoints';
  $('runLog').textContent='Resuming saved run '+run.id+' with its original inputs and scientific settings.';
  $('checkpointStatus').textContent='Resuming '+run.id+' in a new attempt. The saved inputs and settings are used; this analysis remains incomplete until every required check finishes.';
  try{restoreConfiguration(d.config);renderRunSummary();}catch(e){toast('The saved run has started, but its settings could not be displayed: '+e.message);}
  pollJob();
 }catch(e){checkpointResumeUncertain=true;$('runDesign').disabled=wasDisabled;$('checkpointStatus').textContent='Resume request was not confirmed: '+e.message+' Refresh saved progress before trying again.';toast(e.message);}
 finally{checkpointBusy=false;updateCheckpointControls();}
}
function initializeCheckpoints(){
 $('checkpointPanel').classList.toggle('hidden',!!B.offline);
 $('refreshCheckpoints').onclick=refreshCheckpoints;$('resumeCheckpoint').onclick=resumeCheckpoint;$('checkpointChoice').onchange=renderCheckpointDetail;
 updateCheckpointControls();if(!B.offline)refreshCheckpoints();
}
// End saved-checkpoint controls.
$('runDesign').onclick=async()=>{try{if(!S.bins.length)throw Error('Create tracing regions first.');if(S.bins.some(b=>!S.assignments[b.id]))throw Error('Every region needs an imager.');const p=getPolicy();if(!['+','-'].includes(p.bound_strand))throw Error('Choose ONE fixed DNA-bound strand in Step 4 for ALL tracing regions. Per-gene strand switching is not permitted for RASER-FISH tracing.');const missing=transcriptionIssues();if(missing.length)throw Error(missing.join(' '));if(!p.screen_rna&&!confirm('The public annotated-RNA database search is OFF. Continue with reduced RNA QC scope?'))return;if(p.human_policy==='report_only'&&!confirm('Human genomic matches will NOT influence selection, even if there are many high-Tm hits. Continue with this explicit override?'))return;const d=await api('/api/design',{config:configuration()});S.job=d.id;S.designActive=true;pollFailures=0;updateCheckpointControls();$('runDesign').disabled=true;$('cancelRun').classList.remove('hidden');$('runState').textContent='Running';pollJob();}catch(e){$('runLog').textContent=e.message;toast(e.message);}};
let pollFailures=0;
function scheduleJobPoll(job,delay){setTimeout(()=>{if(S.job===job)return pollJob();},delay);}
function showAutoRecovery(recovery){
 const el=$('autoRecoveryStatus');if(!el)return;
 if(!recovery){el.textContent='Automatic restart: ON for new Standard runs. Verified checkpoints are reused after an interruption.';return;}
 const count=Number.isSafeInteger(recovery.attempts)?recovery.attempts:0,limit=Number.isSafeInteger(recovery.max_attempts)?recovery.max_attempts:2;
 const reason=typeof recovery.reason==='string'&&recovery.reason?' '+recovery.reason:'';
 if(recovery.state==='pending'){
  const delay=Math.max(0,Math.ceil((Number(recovery.next_retry_at)||0)-Date.now()/1000));
  el.textContent=`INCOMPLETE — automatic retry ${Math.min(count+1,limit)} of ${limit}${delay?' in '+delay+' seconds':' is waiting'}. Verified progress will be reused.`+reason;
 }else if(recovery.state==='active')el.textContent=count?`Automatic recovery attempt ${count} of ${limit} is running. The analysis remains INCOMPLETE until every required check finishes.`:'Automatic restart is enabled. Verified checkpoints will be reused if this run is interrupted.';
 else if(recovery.state==='stopped')el.textContent='Automatic recovery stopped.'+reason+' This analysis remains INCOMPLETE.';
 else if(recovery.state==='complete')el.textContent='Automatic recovery finished. The completed result is verified separately below.';
 else el.textContent='Automatic restart is enabled. Verified checkpoints will be reused if this run is interrupted.';
}
async function pollJob(){const job=S.job;try{
 const response=await fetch('/api/job?id='+encodeURIComponent(job));if(job!==S.job)return;
 if(!response.ok){if(response.status>=400&&response.status<500){const failure=await response.json();$('runState').textContent='Status unavailable — '+(failure.error||'This run cannot be verified.');$('runDesign').disabled=true;$('cancelRun').classList.add('hidden');return;}throw Error('Cannot read run status.');}
 const d=await response.json();if(job!==S.job)return;
 if(pollFailures){try{await refreshActiveDesign(false);}catch{}if(job!==S.job)return;}pollFailures=0;
 $('progressStage').textContent=d.stage;$('runLog').textContent=(d.messages||[]).slice(-120).join('\n');$('runLog').scrollTop=$('runLog').scrollHeight;$('runState').textContent=d.status;
 const recovery=d.auto_recovery;showAutoRecovery(recovery);
 if(['active','complete','stopped'].includes(recovery?.state)&&typeof recovery.latest_run_id==='string'&&recovery.latest_run_id&&recovery.latest_run_id!==job){
  S.job=recovery.latest_run_id;S.designActive=true;$('runDesign').disabled=true;$('cancelRun').classList.remove('hidden');$('runState').textContent='Following automatic recovery';scheduleJobPoll(S.job,0);return;
 }
 if(recovery?.state==='pending'){
  S.designActive=true;$('runDesign').disabled=true;$('cancelRun').classList.remove('hidden');$('runState').textContent='INCOMPLETE — waiting for automatic retry';scheduleJobPoll(job,1500);return;
 }
 // Worker cleanup precedes the durable retry decision; keep observing that transition.
 if(recovery?.state==='active'&&!['running','complete','completed'].includes(d.status)){
  S.designActive=true;$('runDesign').disabled=true;$('cancelRun').classList.remove('hidden');
  $('runState').textContent='INCOMPLETE — checking automatic recovery';$('autoRecoveryStatus').textContent='The worker stopped. Waiting for cleanup and the automatic recovery decision; saved progress remains incomplete.';
  scheduleJobPoll(job,1500);return;
 }
 if(d.status==='running'&&recovery?.state!=='stopped'){S.designActive=true;$('runDesign').disabled=true;$('cancelRun').classList.remove('hidden');scheduleJobPoll(job,1500);return;}
 if(!['complete','completed'].includes(d.status)){S.designActive=false;refreshCheckpoints();$('runDesign').disabled=false;$('cancelRun').classList.add('hidden');$('runState').textContent='INCOMPLETE — '+(d.status==='running'?'automatic recovery stopped':d.status);toast(d.error||recovery?.reason||'This analysis is incomplete. No completed result or order was produced.');return;}
 const responseResult=await fetch('/api/result?id='+encodeURIComponent(job));if(job!==S.job)return;if(!responseResult.ok)throw Error('Result is not available yet.');const result=await responseResult.json();if(job!==S.job)return;S.result=result;
 S.designActive=false;refreshCheckpoints();$('runDesign').disabled=false;$('cancelRun').classList.add('hidden');S.stale=false;initResult();go(6);toast('New selection completed.');
 }catch(e){if(job!==S.job)return;pollFailures++;$('runState').textContent='Connection unavailable — retrying. Last progress is shown below; completion is unverified.';$('runDesign').disabled=true;scheduleJobPoll(job,Math.min(15000,1500*pollFailures));}}
$('cancelRun').onclick=async()=>{
 if(!S.job||$('cancelRun').disabled)return;const job=S.job;$('cancelRun').disabled=true;
 try{await refreshActiveDesign(false);await api('/api/cancel',{id:job});$('runState').textContent='Cancellation requested — automatic recovery will stop.';toast('Cancellation requested, including pending automatic retries.');}
 catch(e){toast('Cancellation was not confirmed: '+e.message);}
 finally{$('cancelRun').disabled=false;}
};
// SVG drawing: coordinates are always calculated from data, never screenshots.
function axis(ref,a,b,w=1100){const x=p=>90+(p-a)/(b-a+1)*(w-115);return {x,w,plot:w-115};}
function renderMap(ref,bins=[],candidates=[],a=1,b=ref.length,selected=null){const {x,w,plot}=axis(ref,a,b);const visible=(f)=>f.start<=b&&f.end>=a;let svg=[],y=28;const txt=(px,py,s,extra='')=>`<text x="${px}" y="${py}" font-size="11" fill="var(--muted)" ${extra}>${esc(s)}</text>`;
 svg.push(txt(5,y,ref.id),`<line x1="90" x2="${w-25}" y1="${y+5}" y2="${y+5}" stroke="#aebfc7"/>`);for(let i=0;i<=5;i++){let pos=a+i*(b-a)/5;svg.push(txt(x(pos),y-2,fmt(Math.round(pos)),'text-anchor="middle"'));}
 y+=35;const features=(ref.features||[]).filter(visible),genes=features.filter(f=>f.type!=='H3K27ac');const laneEnds=[-Infinity,-Infinity,-Infinity,-Infinity];
 svg.push(txt(4,y+15,'Features'));
 for(const [i,f] of genes.slice(0,500).entries()){let la=0;while(la<3&&laneEnds[la]>=Math.max(a,f.start))la++;laneEnds[la]=Math.min(b,f.end)+(b-a)/50;const fy=y+la*23,fx=x(Math.max(a,f.start)),width=Math.max(2,x(Math.min(b,f.end)+1)-fx),h=f.type==='exon'?9:5;svg.push(`<g data-feature="${esc(i)}"><title>${esc(f.label+' | '+f.type+' | '+f.start+'-'+f.end+' | '+(f.source||''))}</title><rect x="${fx}" y="${fy}" width="${width}" height="${h}" fill="${f.type==='exon'?'#88bad8':'#78a4b9'}"/>${f.type!=='exon'?txt(fx,fy-3,f.label+(f.strand==='-'?' <':f.strand==='+'?' >':' ?')):''}</g>`);}
 if(!genes.length)svg.push(txt(100,y+13,'No gene/features loaded for this coordinate system.'));
 y+=105;const peaks=features.filter(f=>f.type==='H3K27ac'),sig=ref.h3k_signal||[],meta=ref.h3k_metadata;svg.push(txt(4,y+17,'H3K27ac'));
 const values=peaks.map(f=>f.signal||1).concat(sig.map(s=>s.value||0)),max=Math.max(1,...values);let dataPresent=false;
 for(const f of peaks){const fx=x(Math.max(a,f.start)),width=Math.max(2,x(Math.min(b,f.end)+1)-fx),h=10+30*Math.max(0,f.signal||1)/max;svg.push(`<rect x="${fx}" y="${y+43-h}" width="${width}" height="${h}" fill="#bd8240" opacity=".65"><title>${esc(f.label+' | '+(f.cell||'user track')+' | '+(f.source||''))}</title></rect>`);dataPresent=true;}
 for(const f of sig.filter(visible)){if(f.value==null)continue;const fx=x(Math.max(a,f.start)),width=Math.max(1,x(Math.min(b,f.end)+1)-fx),h=40*Math.max(0,f.value)/max;svg.push(`<rect x="${fx}" y="${y+43-h}" width="${width}" height="${h}" fill="#b47537"/>`);dataPresent=true;}
 if(!dataPresent)svg.push(txt(100,y+23,meta?'Track loaded; no peak/signal interval in this view.':'Not loaded. No cell-line data are inferred.'));
 if(meta)svg.push(txt(95,y+58,(meta.is_proxy?'PROXY ':'')+(meta.cell||'')+' | '+(meta.accession||'')+' | '+(meta.output_type||'')));
 y+=83;svg.push(txt(4,y+16,'Regions'));
 for(const bn of bins.filter(bn=>bn.reference===ref.id&&visible(bn))){let fx=x(Math.max(a,bn.start)),width=Math.max(2,x(Math.min(b,bn.end)+1)-fx);svg.push(`<g data-bin="${esc(bn.id)}"><rect x="${fx}" y="${y}" width="${width}" height="28" fill="${color(bn.imager)}" opacity="${selected===bn.id?'.65':'.25'}" stroke="var(--paper)"/><title>${esc(bn.id+' | '+bn.length+' bp | '+bn.imager)}</title>${width>75?txt(fx+5,y+19,bn.id):''}</g>`);}
 y+=50;svg.push(txt(4,y+10,'Primaries'));
 for(const c of candidates){for(const h of c.sites||[]){if(h.reference!==ref.id||!visible(h))continue;let fx=x(Math.max(a,h.start)),width=Math.max(2,x(Math.min(b,h.end)+1)-fx);let fy=y+(c.id===S.selected?15:0);svg.push(`<g data-candidate="${esc(c.id)}"><rect x="${fx}" y="${fy}" width="${width}" height="12" fill="${color(c.imager)}" ${c.id===S.selected?'stroke="#bf8725" stroke-width="3"':''}/><title>${esc((c.name||c.id)+' | '+h.start+'-'+h.end+' | binds '+h.strand)}</title></g>`);}}
 if(!candidates.length)svg.push(txt(100,y+12,'No selected probe sites in this panel yet.'));
 y+=45;
 if(b-a<=95&&ref.sequence){const seq=ref.sequence.slice(a-1,b);svg.push(txt(3,y,'TOP 5\u2032'),txt(3,y+24,'BOTTOM 3\u2032'));for(let i=0;i<seq.length;i++){let xx=x(a+i+.5);svg.push(txt(xx,y,seq[i],'text-anchor="middle"'),txt(xx,y+24,complement(seq[i]),'text-anchor="middle"'));}y+=49;}
 svg.push(txt(92,y,'Positions are 1-based inclusive. Regions are adjacent; individual primary footprints retain their configured minimum gap.'));return `<svg viewBox="0 0 ${w} ${y+15}" role="img" aria-label="Sequence, features, H3K27ac, regions and probe binding sites">${svg.join('')}</svg>`;}
function contextMap(){const r=S.refs.find(r=>r.id===$('contextRef').value)||S.refs[0];if(!r)return;$('contextMap').innerHTML=renderMap(r,[],[]);$('contextSources').textContent=(r.source?'Sequence source: '+r.source.provider+'; ':'')+(r.h3k_metadata?JSON.stringify(r.h3k_metadata):'Public cell-line tracks have not been assigned to this sequence.');}
function renderBinMap(rid){const r=S.refs.find(r=>r.id===(rid||S.refs[0]?.id));if(r)$('binMap').innerHTML=renderMap(r,S.bins,[]);}
document.addEventListener('click',e=>{const b=e.target.closest('[data-bin-preview]');if(b){let bin=S.bins.find(x=>x.id===b.dataset.binPreview);renderBinMap(bin.reference);}});
function initResult(){const r=S.result;if(!r)return;S.viewRef=r.references[0].id;S.viewStart=1;S.viewEnd=r.references[0].length;S.region=null;S.selected=r.selected_ids[0]||null;refOptions('resultReference',r.references);$('resultReference').value=S.viewRef;}
function resultCandidates(){return S.result?.candidates||[];}
function selectedCandidates(){const ids=new Set(S.result?.selected_ids||[]);return resultCandidates().filter(c=>ids.has(c.id));}
function miniHist(hist){let max=Math.max(1,...hist);return `<div class="mini-distribution" title="${hist.join(', ')}">${hist.map(n=>`<span style="height:${Math.max(2,23*n/max)}px"></span>`).join('')}</div>`;}
function resultStrandReview(r){
 const cs=(r?.candidates||[]).filter(c=>c.status==='selected'),local=new Set(),byChr={};let valid=true;
 for(const c of cs){const st=c.bound_strand||c.sites?.[0]?.strand;if(!['+','-'].includes(st)||!c.sites?.length)valid=false;local.add(st);
  for(const h of c.sites||[]){if(h.strand!==st)valid=false;local.add(h.strand);const ref=(r.references||[]).find(x=>x.id===h.reference),m=ref?.mapping;
   if(m){const ori=m.orientation||'+',gs=ori==='+'?h.strand:(h.strand==='+'?'-':'+');if(!['+','-'].includes(ori)||!m.verified)valid=false;const k=m.chromosome;(byChr[k]??=new Set()).add(gs);}}}
 const mixed=local.size>1||Object.values(byChr).some(x=>x.size>1),ok=cs.length>0&&valid&&!mixed;
 const text=mixed?'MIXED STRANDS: this saved design is not suitable as one continuous RASER-FISH trace. Do not combine it without fixed-strand reselection.':!cs.length?'No selected sites to check.':!valid?'Incomplete strand evidence: review before treating this as one trace.':`All saved selected sites bind ${[...local][0]} DNA. ${r.strand_audit?'Fixed-strand safeguard recorded.':'Older result: same sequence orientation only; no new safeguarded run was performed.'}`;
 return {ok,mixed,text};
}
function renderResults(){let r=S.result;coveragePanel?.setCurrentReport(r?.coverage_report,r?.coverage_report_error);if(r){$('selectionSearchReport').textContent=JSON.stringify(r.selection_search||'Not recorded in this older result',null,2);const ss=r.candidates.filter(c=>c.status==='selected'),ll=ss.map(c=>c.arm_length);$('lengthResults').textContent=ss.length?`${ss.filter(c=>c.length_band==='extended').length} extended primaries; selected arm lengths ${Math.min(...ll)}-${Math.max(...ll)} nt. ${r.policy.selection_goal==='maximize_sites'?'Count benchmark was NOT used as a stopping cap.':'Target-count / older selection mode.'}`:'No primaries selected.';}$('resultEmpty').classList.toggle('hidden',!!r);$('resultBody').classList.toggle('hidden',!r);if(!r)return;if(!S.viewRef)initResult();$('resultTitle').textContent=r.run_name||r.title;$('resultScope').classList.toggle('danger',r.status.includes('PREVIEW')||r.policy?.screen_cross_primary===false||r.analysis_complete===false);$('resultScope').innerHTML='<b>'+esc(r.status)+'</b><br>'+esc((r.warnings||[]).join(' '));$('staleResult').classList.toggle('hidden',!S.stale);const strandReview=resultStrandReview(r);$('resultScope').innerHTML+='<div class="notice '+(strandReview.ok?'':'danger')+'"><b>Trace-wide strand check</b><br>'+esc(strandReview.text)+'</div>';$('downloadOrder').disabled=!strandReview.ok||r.analysis_complete===false||/INCOMPLETE|STOPPED|FAILED|CANCELLED/.test(r.status);
 $('resultKpis').innerHTML=[[r.counts.regions,'equal regions'],[r.counts.candidates,'candidates given detailed QC'],[r.counts.selected,'selected primary oligos'],[r.counts.exact_sites,'exact target sites'],[r.counts.QC_rejected,'candidates rejected by QC']].map(([a,b])=>`<div class="kpi"><strong>${fmt(a)}</strong><span>${b}</span></div>`).join('');
 $('resultChecks').innerHTML='<div class="checks-grid">'+(r.checks||[]).map(c=>`<div>${esc(c.name)}<br><b>${esc(c.status)}</b></div>`).join('')+'</div>';$('resultSettings').textContent=JSON.stringify(r.policy,null,2);$('resultTranscription').textContent=JSON.stringify(r.transcription_summary||'NOT ASSESSED IN THIS OLDER RESULT',null,2);$('resultParameterSources').textContent=JSON.stringify(r.parameter_provenance||'No preset provenance in this older result',null,2);$('resultSources').textContent=r.scope+' '+(r.annotation_provenance?.length?'Annotation metadata are saved in the result JSON.':'No public annotation track was used in this run.');
 $('resultRegions').innerHTML=r.bins.map(b=>{const m=b.metrics;return `<tr class="clickable ${S.region===b.id?'selected':''}" data-result-bin="${esc(b.id)}"><td><b>${esc(b.id)}</b></td><td>${esc(b.imager)}</td><td>${fmt(b.length)}</td><td>${m.probes}</td><td>${m.min_arm_length==null?"not recorded":m.min_arm_length+"-"+m.max_arm_length+" nt"}</td><td>${m.extended_primaries??"not recorded"}</td><td>${m.exact_sites}</td><td>${fmt(m.coverage_percent)}%</td><td>${fmt(m.longest_gap)} bp</td><td>${fmt(m.centroid_offset)} bp</td><td>${miniHist(m.distribution)}</td><td><span class="tag ${m.below_minimum?'reject':''}">${m.shortfall_to_target} below target${m.below_minimum?'; BELOW MINIMUM':''}</span></td></tr>`;}).join('');
 const max=Math.max(1,...r.bins.map(b=>b.metrics.exact_sites),r.policy.target_per_bin),height=65+r.bins.length*36;let s=[];
 s.push('<text x="230" y="20" font-size="12" fill="var(--muted)">Solid = distinct primaries; outline = exact sites including repeats; dashed line = count benchmark / cap in target mode</text>');
 r.bins.forEach((b,i)=>{let y=37+i*36;s.push(`<g data-result-bin="${esc(b.id)}" style="cursor:pointer"><text x="8" y="${y+17}" font-size="12">${esc(b.id)} / ${esc(b.imager)}</text><rect x="230" y="${y}" width="${b.metrics.exact_sites/max*690}" height="22" fill="none" stroke="${color(b.imager)}"/><rect x="230" y="${y+3}" width="${b.metrics.probes/max*690}" height="16" fill="${color(b.imager)}"/><text x="${240+b.metrics.exact_sites/max*690}" y="${y+16}" font-size="12">${b.metrics.probes} primaries / ${b.metrics.exact_sites} sites</text></g>`);});s.push(`<line x1="${230+r.policy.target_per_bin/max*690}" x2="${230+r.policy.target_per_bin/max*690}" y1="30" y2="${height-10}" stroke="var(--muted)" stroke-dasharray="5,4"/>`);$('countChart').innerHTML=`<svg viewBox="0 0 1100 ${height}" role="img" aria-label="Per-region primary and exact site counts">${s.join('')}</svg>`;
 $('generationReport').textContent=r.algorithm+'\n\n'+JSON.stringify({generation:r.generation,candidate_summary:r.candidate_summary,rejection_summary:r.rejection_summary},null,2);$('spatialSummary').textContent=r.bins.map(b=>b.id+' — 1 kb sector counts: '+(b.metrics.sector_counts||[]).join(' ')+'; largest uncovered gap: '+fmt(b.metrics.largest_gap_bp??b.metrics.longest_gap)+' bp; median center spacing: '+fmt(b.metrics.median_spacing_bp)+' bp').join('\n');$('pairSummary').textContent=JSON.stringify({self_dimer:r.self_dimer_screening||r.selected_pair_summary?.self_dimer||'See saved self-pair evidence',cross_primary:crossPrimaryResultStatus(r),selected_interactions:r.selected_pair_summary||'Not summarized in this older result'},null,2);$('crossPrimaryResultStatus').textContent=crossPrimaryResultText(r);$('interactionResultNote').textContent='Compact view lists only failed or near-threshold interactions among the analyses actually run. Self-dimer checks are identified separately. Near-threshold flags do not alter scientific acceptance. This is not an experimental dimer fraction.';$('candidateFilter').disabled=r.result_schema==='raserfish.compact-result.v1';if($('candidateFilter').disabled)$('candidateFilter').value='selected';$('downloadDecisions').disabled=B.offline||!S.job;$('downloadResults').disabled=!S.job||B.offline;$('downloadScores').disabled=!S.job||B.offline;renderResultMap();renderProbeTable();renderInspector();renderPairs();}
function renderResultMap(){if(!S.result)return;const ref=S.result.references.find(r=>r.id===S.viewRef)||S.result.references[0];S.viewStart=Math.max(1,Math.min(S.viewStart,ref.length));S.viewEnd=Math.max(S.viewStart,Math.min(S.viewEnd,ref.length));$('viewStart').value=S.viewStart;$('viewEnd').value=S.viewEnd;$('resultMap').innerHTML=renderMap(ref,S.result.bins,selectedCandidates(),S.viewStart,S.viewEnd,S.region);$('selectedRegionNote').textContent=S.region?'Inspecting '+S.region+'. Select Whole sequence to compare its neighbors.':'Whole-reference view. Click a region to inspect its binding-site distribution.';}
function chooseBin(id){const b=S.result?.bins.find(b=>b.id===id);if(!b)return;S.region=id;S.viewRef=b.reference;S.viewStart=b.start;S.viewEnd=b.end;$('resultReference').value=b.reference;renderResults();}
document.addEventListener('click',e=>{let b=e.target.closest('[data-result-bin]');if(b){chooseBin(b.dataset.resultBin);return;}b=e.target.closest('[data-bin]');if(b&&S.step===6){chooseBin(b.dataset.bin);return;}b=e.target.closest('[data-candidate]');if(b){S.selected=b.dataset.candidate;renderProbeTable();renderInspector();renderResultMap();}});
$('resultReference').onchange=()=>{S.viewRef=$('resultReference').value;let r=S.result.references.find(r=>r.id===S.viewRef);S.viewStart=1;S.viewEnd=r.length;S.region=null;renderResults();};
$('applyView').onclick=()=>{S.viewStart=Math.floor(Number($('viewStart').value));S.viewEnd=Math.floor(Number($('viewEnd').value));if(!Number.isFinite(S.viewStart)||!Number.isFinite(S.viewEnd)||S.viewEnd<S.viewStart){toast('Invalid map interval.');return;}renderResultMap();renderProbeTable();};
$('wholeView').onclick=()=>{const r=S.result.references.find(r=>r.id===S.viewRef);S.viewStart=1;S.viewEnd=r.length;S.region=null;renderResults();};
function zoom(f){const ref=S.result?.references.find(r=>r.id===S.viewRef);if(!ref)return;const mid=(S.viewStart+S.viewEnd)/2,w=Math.min(ref.length,Math.max(30,(S.viewEnd-S.viewStart+1)*f));S.viewStart=Math.max(1,Math.floor(mid-w/2));S.viewEnd=Math.min(ref.length,S.viewStart+Math.ceil(w)-1);renderResultMap();renderProbeTable();}
$('zoomIn').onclick=()=>zoom(.5);$('zoomOut').onclick=()=>zoom(2);
function shownCandidates(){const search=$('probeSearch').value.toUpperCase(),filter=$('candidateFilter').value;return resultCandidates().filter(c=>c.reference===S.viewRef&&(!S.region||c.bin_id===S.region)&&(!search||((c.name||c.id)+' '+c.arm).toUpperCase().includes(search))&&(filter==='all'||filter==='selected'&&c.status==='selected'||filter==='rejected'&&c.status==='rejected_by_QC'||filter==='notselected'&&c.status.startsWith('not_selected'))).sort((a,b)=>a.start-b.start||a.id.localeCompare(b.id));}
function renderProbeTable(){const cs=shownCandidates();$('probeCount').textContent=`${cs.length} candidates shown. Display filters do not rerun selection or modify the saved order.`;$('probeRows').innerHTML=cs.slice(0,1500).map(c=>`<tr class="clickable ${S.selected===c.id?'selected':''}" data-candidate="${esc(c.id)}"><td>${esc(c.name||c.id)}</td><td>${esc(c.bin_id)}</td><td>${c.arm_length} nt${c.length_band==='extended'?' <span class="tag">extended</span>':''}</td><td>${fmt(c.tm)}</td><td>${fmt(c.gc)}%</td><td>${c.specificity.genome_scope==='NOT_RUN'?'local only: ':''}${c.specificity.high_tm_sites}</td><td>${c.fold?.mean_arm_unpaired==null?'NOT RUN':fmt(100*c.fold.mean_arm_unpaired)+'%'}</td><td><span class="tag ${c.status==='selected'?'selected':c.status==='rejected_by_QC'?'reject':''}">${esc(c.status.replaceAll('_',' '))}</span></td></tr>`).join('');}
$('probeSearch').oninput=renderProbeTable;$('candidateFilter').onchange=renderProbeTable;
function renderInspector(){const c=resultCandidates().find(c=>c.id===S.selected);if(!c){$('probeInspector').textContent='Select a candidate or primary.';return;}const q=c.specificity,ff=c.fold,other=c.offtarget_hits||[];let html=`<h3 style="overflow-wrap:anywhere">${esc(c.name||c.id)}</h3><span class="tag ${c.status==='selected'?'selected':'reject'}">${esc(c.status.replaceAll('_',' '))}</span><p class="hint">${esc(c.bin_id)} / ${esc(c.imager)} / binds ${(c.bound_strand||c.sites?.[0]?.strand)==='-'?'BOTTOM (-)':'TOP (+)'} strand</p>`;
 html+=`<div class="facts"><div>Arm Tm <button class="help" data-help="tm">?</button><strong>${fmt(c.tm)} C</strong></div><div>GC <button class="help" data-help="gc">?</button><strong>${fmt(c.gc)}%</strong></div><div>Arm / complete primary<strong>${c.arm_length} / ${c.full_length} nt</strong></div><div>Exact sites in its region<strong>${c.sites.length}</strong></div><div>Unintended high-Tm site clusters <button class="help" data-help="high_hit">?</button><strong>${q.high_tm_sites}</strong></div><div>Mean arm unpaired <button class="help" data-help="unpaired">?</button><strong>${ff?.mean_arm_unpaired==null?'NOT RUN':fmt(100*ff.mean_arm_unpaired)+'%'}</strong></div></div>`;
 html+=`<p class="hint"><b>Length selection:</b> ${esc(c.extension_note||'Fixed-range / older result.')} ${c.extension_bases?`+${c.extension_bases} bases beyond the preferred maximum.`:''} <button class="help" data-help="adaptive_length">?</button></p>`;
 if(c.reject_reasons.length)html+='<div class="notice danger"><b>Not selected because:</b><br>'+c.reject_reasons.map(esc).join('<br>')+'</div>';if(c.selection_note)html+='<p class="hint">'+esc(c.selection_note)+'</p>';
 html+=`<details open><summary>Gene/RNA avoidance <button class="help" data-help="rna_direction">?</button></summary><p class="hint">${esc(c.rna_avoidance?.status||'NOT ASSESSED IN THIS OLDER RESULT')}.</p>${(c.rna_avoidance?.sites||[]).map(h=>`<p class="hint">${esc(h.reference)}:${fmt(h.start)}-${fmt(h.end)} binds ${esc(h.bound_strand)} DNA; overlapping transcript directions: ${esc(h.transcription_strands.join(', ')||'none in supplied model')}. ${esc(h.genes.join('; '))}</p>`).join('')}<p class="hint">Complementary RNA-risk clusters: ${q.rna_risk_sites??q.strong_rna_sites??'not calculated'}. Same-sense RNA is not a hybridization match; annotation coverage limits this assessment.</p></details>`;
 const savedImager=(S.result.imagers||[]).find(i=>i.id===c.imager),savedDock=c.docking_sequence||c.prefix.slice(0,-2),savedSequence=c.imager_sequence||savedImager?.sequence||rc(savedDock);
 html+=`<details open class="saved-readout-pair"><summary>Imaging oligo versus its docking site <button class="help" data-help="imager_input_type">?</button></summary>${imagerPairMarkup({sequence:savedSequence,docking:savedDock},false)}${ImagerBinding.render({id:c.imager,sequence:savedSequence,docking:savedDock,dye:savedImager?.dye},{primary:c,showOverview:true})}<p class="hint">From this SAVED candidate, not any unsaved dropdown change. Input convention: ${esc(inputKindLabel(c.imager_input_kind||savedImager?.input_kind))}.</p></details>`;
 html+=`<h4>Complete primary, 5' to 3'</h4><div class="seq"><span class="handle">${esc(c.prefix.slice(0,-2))}</span>TT<b>${esc(c.arm)}</b>TT<span class="handle">${esc(c.suffix.slice(2))}</span></div><p class="hint">Purple: imager-binding docking sites, NOT the fluorescent imager itself. Bold: DNA-binding arm that recognizes the genomic/construct target. TT: spacers. These complete primaries are used for folding.</p><details open><summary>Exact binding sites and extended footprints <button class="help" data-help="footprint">?</button></summary>${c.sites.map(h=>`<p class="hint">${esc(h.reference)}:${fmt(h.start)}-${fmt(h.end)}; footprint ${fmt(h.extended_start)}-${fmt(h.extended_end)}</p>`).join('')}</details>`;
 html+=`<details><summary>Specificity summary (${c.offtarget_alignment_count??other.length} evaluated alignments) <button class="help" data-help="offtarget">?</button></summary><p class="hint">Allowed sites are shown, not counted as unwanted. High/unknown site counts are applied according to the saved human-hit policy. The compact view shows QC summaries; raw BLAST HSPs are saved only in Full audit/debug mode.</p><pre>${esc(JSON.stringify(q,null,2))}</pre>`;
 const ordered=[...other].sort((a,b)=>Number(a.allowed)-Number(b.allowed)||(b.risk_tm??-999)-(a.risk_tm??-999)||b.alignment_length-a.alignment_length);
 for(const h of ordered.slice(0,100))html+=`<details><summary>${h.allowed?'ALLOWED':'REVIEW'} ${esc(h.source)} ${esc(h.reference)}:${fmt(h.start)}-${fmt(h.end)}</summary><p class="hint">${h.alignment_length} aligned bases; ${fmt(h.identity)}% identity; ${fmt(h.query_coverage)}% of ${h.full_primary?'complete primary':'arm'}. Physical strand: ${esc(h.bound_strand)}.</p><p class="hint">Off-target Tm: <b>${h.offtarget_tm==null?'UNKNOWN / sequence-only':fmt(h.offtarget_tm)+' C'}</b> (${esc(h.tm_model)}) <button class="help" data-help="offtarget_tm">?</button><br>Exact-subsegment proxy: ${h.exact_tract_tm==null?'not available':fmt(h.exact_tract_tm)+' C'}; on-minus-off Tm: ${fmt(h.on_minus_off_tm)}</p><pre>${esc(h.qseq)}\n${esc(h.sseq)}</pre></details>`;
 if(ordered.length>100)html+='<p class="hint">First 100 by risk shown here; full-audit alignment evidence is under debug/.</p>';html+='</details>';
 html+=`<details open><summary>Complete-primary folding <button class="help" data-help="folding">?</button></summary><p class="hint">${esc(ff?.status||'NOT_RUN')}. Aqueous model; no formamide or whole-pool equilibrium.</p>${ff?.structure?`<pre>${esc(c.full)}\n${esc(ff.structure)}</pre>`:ff?.status?.startsWith('RUN')?'<p class="hint">Detailed structure is not retained in concise output.</p>':''}${ff?.mfe_kcal_mol!=null?`<p class="hint">MFE: ${fmt(ff.mfe_kcal_mol)} kcal/mol; ensemble free energy: ${fmt(ff.ensemble_kcal_mol)} kcal/mol; mean docking unpaired: ${ff.mean_docking_unpaired==null?'not recorded':fmt(100*ff.mean_docking_unpaired)+'%'}.</p>`:''}</details>`;
 html+=`<details><summary>Self-dimer check (same primary sequence) <button class="help" data-help="cofold">?</button></summary><pre>${esc(JSON.stringify(c.self_pair,null,2))}</pre></details>`;$('probeInspector').innerHTML=html;}
function crossPrimaryResultStatus(r){
 const recorded=r.cross_primary_screening||r.selected_pair_summary?.cross_primary;
 if(r.policy?.screen_cross_primary===false||recorded?.enabled===false)return {...recorded,enabled:false,status:'NOT RUN - USER DISABLED',within_region_status:'NOT RUN - USER DISABLED',between_region_status:'NOT RUN - USER DISABLED',evaluated:0,selected_evaluated:0};
 return recorded||{enabled:true,status:'RECORDED IN LEGACY RESULT — see saved checks',within_region_status:'See saved checks',between_region_status:'See saved checks'};
}
function crossPrimaryResultText(r){const x=crossPrimaryResultStatus(r);return 'Probe-probe interaction screening: '+(x.enabled?'ON':'OFF')+' — '+x.status+'. Within-region: '+x.within_region_status+'. Between-region: '+x.between_region_status+'. Self-dimer evidence is separate. Order status: DRAFT.';}
function resultPairScope(pair,r){if(pair.a===pair.b)return 'SELF_DIMER';if(pair.scope)return pair.scope;const by=new Map((r.candidates||[]).map(c=>[c.id,c]));const a=by.get(pair.a),b=by.get(pair.b);return a&&b?(a.bin_id===b.bin_id?'WITHIN_REGION':'BETWEEN_REGION'):'CROSS_PRIMARY';}
function renderPairs(){const r=S.result;if(!r)return;const cl=new Map(r.candidates.map(c=>[c.id,c])),pairs=[...(r.selected_pair_interactions||[])].sort((a,b)=>(a.mfe_gain_vs_separate??999)-(b.mfe_gain_vs_separate??999)||b.max_exact_stem-a.max_exact_stem);S.pairView=pairs;$('pairRows').innerHTML=pairs.slice(0,100).map((p,i)=>`<tr class="clickable" data-pair="${i}"><td>${esc(resultPairScope(p,r).replaceAll('_',' '))}</td><td>${esc(cl.get(p.a)?.name||p.a)}</td><td>${esc(cl.get(p.b)?.name||p.b)}</td><td>${p.mfe_gain_vs_separate==null?'NOT RUN':fmt(p.mfe_gain_vs_separate)}</td><td>${p.max_exact_stem==null?'NOT RUN':p.max_exact_stem+' bp'}</td><td>${esc(p.status)}</td></tr>`).join('');if(!pairs.length)$('pairRows').innerHTML='<tr><td colspan="6">No flagged interactions among the checks actually run. This is not evidence for an omitted analysis.</td></tr>';}
document.addEventListener('click',e=>{const el=e.target.closest('[data-pair]');if(!el)return;const p=S.pairView[Number(el.dataset.pair)],by=new Map(S.result.candidates.map(c=>[c.id,c])),a=by.get(p.a),b=by.get(p.b);$('pairDetail').innerHTML=`<div class="pair-detail"><b>${esc(a?.name||p.a)} + ${esc(b?.name||p.b)} (${esc(resultPairScope(p,S.result).replaceAll('_',' '))})</b><pre>${esc(a?.full||'')} &amp; ${esc(b?.full||'')}\n${esc(p.structure||(p.status?.startsWith('RUN')?'Detailed structure is not retained in concise output.':'Full cofold not run'))}</pre><p class="hint">Gain relative to separately folded primaries: ${fmt(p.mfe_gain_vs_separate)} kcal/mol. ${Array.isArray(p.intermolecular_pairs)?p.intermolecular_pairs.length:(p.intermolecular_pairs??0)} modeled intermolecular pairs. Not a dimer fraction at the experimental concentration.</p></div>`;});
function csv(rows){return rows.map(row=>row.map(v=>'"'+String(v??'').replace(/"/g,'""')+'"').join(',')).join('\r\n')+'\r\n';}
function dlFile(name){if(S.job&&!B.offline){window.location.href='/download?id='+encodeURIComponent(S.job)+'&file='+encodeURIComponent(name);return true;}return false;}
$('downloadResults').onclick=()=>dlFile('results.zip');
$('downloadScores').onclick=()=>dlFile('primary_scores.csv');
$('downloadOrder').onclick=()=>{if(!S.result)return;if(S.result.analysis_complete===false||/INCOMPLETE|STOPPED|FAILED|CANCELLED/.test(S.result.status)){toast('This analysis is incomplete; no draft order is available.');return;}if(!resultStrandReview(S.result).ok){toast('A draft order for one RASER trace requires complete, consistent single-strand evidence. Reselect on one fixed strand.');return;}if(!confirm('Export the SAVED selected candidate sequences? This is a DRAFT order, not experimental validation. '+crossPrimaryResultText(S.result)+' Historical orders are unchanged.'))return;if(!dlFile('DRAFT_IDT_oPools.csv'))saveBlob('DRAFT_'+S.result.status+'_oPools.csv',csv([['Pool name','Sequence'],...selectedCandidates().map(c=>[c.bin_id+'_'+c.imager,c.full])]),'text/csv');};
$('downloadRegion').onclick=()=>{if(!dlFile('region_summary.csv'))saveBlob('region_summary.csv',csv([['Region','Reference','Start','End','Width','Imager','Primaries','Exact sites','Coverage percent','Longest gap','Centroid offset','Count shortfall'],...S.result.bins.map(b=>[b.id,b.reference,b.start,b.end,b.length,b.imager,b.metrics.probes,b.metrics.exact_sites,b.metrics.coverage_percent,b.metrics.longest_gap,b.metrics.centroid_offset,b.metrics.shortfall_to_target])]),'text/csv');};
$('downloadDecisions').onclick=()=>{if(!dlFile('candidate_decisions.csv.gz'))saveBlob('candidate_decisions.csv',csv([['Candidate','Selected name','Region','Imager','Imager_Sequence_5to3','Docking_Sequence_in_Primary_5to3','Imager_Input_Type','Arm','Full primary','Tm','Decision','Reasons','Selection note'],...resultCandidates().map(c=>[c.id,c.name,c.bin_id,c.imager,c.imager_sequence,c.docking_sequence||c.prefix.slice(0,-2),c.imager_input_kind||'not recorded',c.arm,c.full,c.tm,c.status,c.reject_reasons.join('; '),c.selection_note])]),'text/csv');};
$('exportResultHtml').onclick=async()=>{try{if(S.job&&!B.offline){dlFile('viewer.html');return;}const clone=document.documentElement.cloneNode(true),boot=copy(B);boot.demo=S.result;boot.imagers=copy(S.result?.imagers||S.imagers);boot.offline=true;boot.preview=false;boot.initial_step=6;boot.token='';boot.environment={offline:true};boot.references=[];delete boot.default_config;let found=false;for(const s of clone.querySelectorAll('script')){if(s.src?.includes('bootstrap.js')||s.textContent.startsWith('window.STUDIO_BOOT=')){s.removeAttribute('src');s.textContent='window.STUDIO_BOOT='+JSON.stringify(boot).replace(/</g,'\\u003c')+';';found=true;}else if(['/web/app.js','/web/imager_binding.js','/web/probe_coverage.js','/web/coverage_panel.js'].includes(s.getAttribute('src'))){const url=s.getAttribute('src');s.removeAttribute('src');s.textContent=await (await fetch(url)).text();}}if(!found)throw Error('Embedded bootstrap not found.');const l=clone.querySelector('link[rel="stylesheet"]');if(l){const style=document.createElement('style');style.textContent=await (await fetch('/web/styles.css')).text();l.replaceWith(style);}saveBlob('RASERFISH_Design_Results.html','<!doctype html>\n'+clone.outerHTML,'text/html');}catch(e){toast(e.message);}};
$('showDemo').onclick=()=>{if(!B.demo){toast('No example is bundled.');return;}S.result=copy(B.demo);S.job=null;S.stale=false;initResult();go(6);};
document.querySelectorAll('.step[data-step="4"] input,.step[data-step="4"] select,.step[data-step="3"] input').forEach(el=>{el.addEventListener('change',markStale);el.addEventListener('input',markStale);});
$('connection').textContent=B.offline?'Offline preview / saved result':'Local computer';$('connection').classList.toggle('warning',!!B.offline);
if(B.offline){$('globalNotice').classList.remove('hidden');$('globalNotice').innerHTML='<b>Offline '+(B.preview?'interface preview':'saved result')+'.</b> Region planning, dropdowns and result inspection work here. Live annotation retrieval and QC-guided computation require the local Studio program. Check the saved execution checks for the scope of this result. Raw alignments are intentionally excluded from compact viewers.';}
const coveragePanel=window.CoveragePanel?.create({offline:!!B.offline});
initializeCompletedRuns();initializeCheckpoints();initializePresets();showEnvironment();loadedRefs();if(B.default_config&&(!B.offline||B.preview))restoreConfiguration({...copy(B.default_config),references:copy(B.references)},false);renderImagerInventory();updateImagerPreview();if(S.result)initResult();go(B.offline&&Number.isInteger(B.initial_step)&&B.initial_step>=1&&B.initial_step<=6?B.initial_step:B.offline&&!B.preview?6:1);
window.designStudio={state:S,restoreConfiguration,go,makeBins,configuration,renderResults,help,renderTranscription,transcriptionIssues,applyPreset,resolveImagerEntry,updateImagerPreview,renderImagers,renderBindingPreview};
