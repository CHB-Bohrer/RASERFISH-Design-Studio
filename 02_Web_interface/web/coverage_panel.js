/* Saved-data visualization only: never submits a design or changes its policy. */
(function (root, factory) {
  const api = factory();
  if (typeof module === 'object' && module.exports) module.exports = api;
  if (root) root.CoveragePanel = api;
})(typeof window !== 'undefined' ? window : null, function () {
  'use strict';
  const escape = value => String(value ?? '').replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  function filename(title) { return String(title || 'probe_coverage').replace(/[^a-zA-Z0-9._-]+/g, '_').slice(0, 100) || 'probe_coverage'; }
  function create(options = {}) {
    const doc = options.document || document, get = id => doc.getElementById(id);
    const renderer = options.renderer || window.ProbeCoverage;
    const request = options.fetch || window.fetch.bind(window);
    const offline = !!options.offline;
    let sources = [], current = null, report = null, requestId = 0, loaded = false, pending = false, svg = '';
    const source = get('coverageSource'), region = get('coverageRegion'), mode = get('coverageMode'), status = get('coverageStatus');
    if (!source || !renderer) return {refreshSources(){}, setCurrentReport(){}};
    function controls() {
      source.disabled = pending || (!sources.length && !current);
      get('refreshCoverage').disabled = pending || offline;
      for (const id of ['coverageRegion','coverageMode','coverageTranscripts']) get(id).disabled = pending || !report;
      for (const id of ['coverageSvg','coveragePng','coverageHtml']) get(id).disabled = pending || !svg;
      get('coverageRegionLabel').classList.toggle('hidden', !report);
    }
    function sourceOptions(preferred) {
      source.innerHTML = (current ? '<option value="current">Opened result</option>' : '') +
        sources.map((s,i) => `<option value="saved:${i}">${escape(s.label)}</option>`).join('');
      const valid = preferred === 'current' ? !!current : /^saved:\d+$/.test(preferred || '') && !!sources[Number(preferred.split(':')[1])];
      source.value = valid ? preferred : sources.length ? 'saved:0' : current ? 'current' : '';
      controls();
    }
    function draw() {
      if (!report) {svg = ''; get('coverageFigure').innerHTML = ''; controls(); return;}
      try {
        svg = renderer.render(report, {mode:mode.value,regionId:region.value || undefined,showTranscripts:get('coverageTranscripts').checked});
        if(svg.includes('data-coverage-error="true"'))throw Error('Saved coordinates could not be drawn. Review the saved coverage data.');
        get('coverageFigure').innerHTML = svg;
        get('coverageProvenance').textContent = [report.annotation?.source, report.annotation?.note, ...(report.warnings || [])].filter(Boolean).join('\n');
        get('coverageScope').textContent = 'DRAFT · Cross-primary interactions: ' + (report.cross_primary_status || 'NOT RECORDED');
        get('coverageTitle').textContent = report.title || 'Saved probe coverage';
        const shown = report.regions.filter(r => !region.value || r.id === region.value);
        const sites = shown.reduce((n,r) => n+(Number(r.metrics?.sites)||0),0);
        status.textContent = `${shown.length}${shown.length < report.regions.length ? ' of '+report.regions.length : ''} region${shown.length === 1 && shown.length === report.regions.length ? '' : 's'} · ${sites} selected binding sites · drawn from saved results`;
      } catch (e) {
        svg = ''; get('coverageFigure').innerHTML = ''; status.textContent = 'Could not draw this saved view: '+e.message;
      }
      controls();
      if (!svg) for (const id of ['coverageSvg','coveragePng','coverageHtml']) get(id).disabled = true;
    }
    function display(value) {
      if (!value || !Array.isArray(value.regions)) throw Error('The saved coverage data are incomplete.');
      report = value;
      const previous = region.value;
      region.innerHTML = '<option value="">All regions</option>'+report.regions.map(r => `<option value="${escape(r.id)}">${escape(r.label || r.id)} / ${escape(r.imager)}</option>`).join('');
      region.value = report.regions.some(r => r.id === previous) ? previous : '';
      draw();
    }
    async function chooseSource() {
      const generation = ++requestId, selected = source.value;
      report = null; svg = ''; get('coverageFigure').innerHTML = ''; get('coverageProvenance').textContent = ''; get('coverageScope').textContent = ''; get('coverageTitle').textContent = 'Saved probe coverage';
      if (selected === 'current') {pending=false; if(current)display(current); else controls(); return;}
      const item = sources[Number(selected.split(':')[1])];
      if (offline || !item) {pending=false; status.textContent='Open a completed result to see its coverage.'; controls(); return;}
      pending=true; controls(); status.textContent='Loading saved coverage and transcript annotations…';
      try {
        const response=await request('/api/coverage?kind='+encodeURIComponent(item.kind)+'&id='+encodeURIComponent(item.id),{cache:'no-store'});
        const data=await response.json(); if(generation!==requestId)return;
        if(!response.ok)throw Error(data.error || 'Could not read saved coverage.');
        display(data);
      } catch(e) {if(generation===requestId)status.textContent=e.message;}
      finally {if(generation===requestId){pending=false;controls();}}
    }
    async function refreshSources(force=false) {
      if(offline || pending || (loaded && !force))return;
      const prior=sources[Number(source.value.split(':')[1])], keepCurrent=loaded && source.value==='current';
      pending=true; controls();
      try {
        const response=await request('/api/coverage-sources',{cache:'no-store'}), data=await response.json();
        if(!response.ok)throw Error(data.error || 'Restart the updated designer to load coverage views.');
        if(!Array.isArray(data.sources))throw Error('The saved coverage list is unreadable.');
        sources=data.sources.filter(s=>s && ['run','experiment'].includes(s.kind) && typeof s.id==='string');
        loaded=true;
        const index=prior ? sources.findIndex(s=>s.kind===prior.kind && s.id===prior.id) : -1;
        sourceOptions(keepCurrent?'current':index>=0?'saved:'+index:undefined);
        pending=false; await chooseSource();
      }catch(e){status.textContent=e.message;}
      finally{pending=false;controls();}
    }
    function setCurrentReport(value, error) {
      if(!value && error && (offline || source.value==='current' || !sources.length)){
        ++requestId;current=null;report=null;svg='';pending=false;get('coverageFigure').innerHTML='';
        status.textContent=error;controls();return;
      }
      if(value===current)return;
      const follow=source.value==='current' || (!loaded && !sources.length);
      current=value || null; sourceOptions(follow&&current?'current':source.value);
      if(follow)chooseSource();
    }
    function download(blob,name) {
      const url=URL.createObjectURL(blob), a=doc.createElement('a');
      a.href=url; a.download=name; doc.body.appendChild(a); a.click(); a.remove();
      setTimeout(()=>URL.revokeObjectURL(url),1000);
    }
    const exportName=()=>filename(report?.title)+(region.value?'_'+filename(region.value):'')+'_'+mode.value;
    source.onchange=chooseSource;
    get('refreshCoverage').onclick=()=>refreshSources(true);
    for(const el of [region,mode,get('coverageTranscripts')])el.onchange=draw;
    get('coverageSvg').onclick=()=>{if(svg)download(new Blob([svg],{type:'image/svg+xml;charset=utf-8'}),exportName()+'.svg');};
    get('coverageHtml').onclick=()=>{if(svg)download(new Blob(['<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width"><title>'+escape(report.title)+'</title><style>body{margin:24px;font:14px system-ui;color:#182637;background:#fff}svg{display:block;width:100%;height:auto}pre{white-space:pre-wrap}</style>'+svg+'<pre>'+escape(get('coverageProvenance').textContent)+'</pre></html>'],{type:'text/html;charset=utf-8'}),exportName()+'.html');};
    get('coveragePng').onclick=async()=>{
      if(!svg)return;
      const content=svg,name=exportName(), button=get('coveragePng');button.disabled=true;
      try {
        const image=new Image();
        await new Promise((resolve,reject)=>{image.onload=resolve;image.onerror=()=>reject(Error('Could not prepare PNG. The SVG download remains available.'));image.src='data:image/svg+xml;charset=utf-8,'+encodeURIComponent(content);});
        const box=get('coverageFigure').querySelector('svg')?.viewBox.baseVal;
        const w=box?.width || image.naturalWidth, h=box?.height || image.naturalHeight;
        if(!w||!h)throw Error('This view has no drawable dimensions.');
        const scale=Math.min(2,Math.sqrt(24000000/(w*h)),16000/w,16000/h);
        const canvas=doc.createElement('canvas');canvas.width=Math.round(w*scale);canvas.height=Math.round(h*scale);
        const context=canvas.getContext('2d');context.fillStyle='#ffffff';context.fillRect(0,0,canvas.width,canvas.height);context.drawImage(image,0,0,canvas.width,canvas.height);
        const blob=await new Promise(resolve=>canvas.toBlob(resolve,'image/png'));
        if(!blob)throw Error('PNG export failed. Use SVG for this view.');download(blob,name+'.png');
      }catch(e){status.textContent=e.message;}finally{controls();}
    };
    get('coverageFigure').onclick=e=>{
      const hit=e.target.closest('[data-site-id]');
      if(hit)status.textContent=(hit.getAttribute('data-probe-name') || hit.getAttribute('data-site-id'))+' · '+(hit.querySelector('title')?.textContent || 'Selected binding site');
    };
    get('refreshCoverage').classList.toggle('hidden',offline);
    controls();
    return {refreshSources,setCurrentReport};
  }
  return {create,escape,filename};
});
