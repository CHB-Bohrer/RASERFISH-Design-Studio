/* Saved-result coverage plots. Rendering only: no network, design or QC work.
 * Input intervals are local, 1-based and inclusive. SVG edges use [start,end+1).
 * The exported SVG includes its own styles and a white background.
 */
(function (root, factory) {
  'use strict';
  const api = factory();
  if (typeof module === 'object' && module.exports) module.exports = api;
  else root.ProbeCoverage = api;
})(typeof globalThis !== 'undefined' ? globalThis : this, function () {
  'use strict';
  const C = Object.freeze({ink:'#182637', muted:'#647383', blue:'#2b678c',
    plus:'#b23d60', minus:'#7952a3', grid:'#dce3e9', pale:'#f4f7fa', white:'#ffffff'});
  const W = 1320, X = 328, PW = 960, BAR_H = 74;
  const esc = value => String(value ?? '').replace(/[\x00-\x08\x0b\x0c\x0e-\x1f]/g, '').replace(/[&<>"']/g,
    c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  const number = value => value !== null && value !== '' && Number.isFinite(Number(value)) ? Number(value) : null;
  const fmt = value => number(value) === null ? 'unknown' : Number(value).toLocaleString('en-US', {maximumFractionDigits:1});
  const short = (value, max=40) => { const s=String(value ?? ''); return s.length > max ? s.slice(0,max-1)+'…' : s; };
  const n = value => Math.round(value*100)/100;
  const text = (x,y,label,size=12,color=C.ink,extra='') => `<text x="${n(x)}" y="${n(y)}" font-family="Arial,Helvetica,sans-serif" font-size="${size}" fill="${color}" ${extra}>${esc(label)}</text>`;
  const line = (x1,y1,x2,y2,color=C.grid,extra='') => `<line x1="${n(x1)}" y1="${n(y1)}" x2="${n(x2)}" y2="${n(y2)}" stroke="${color}" ${extra}/>`;
  const rect = (x,y,w,h,color,extra='') => `<rect x="${n(x)}" y="${n(y)}" width="${n(w)}" height="${n(h)}" fill="${color}" ${extra}/>`;
  const title = value => `<title>${esc(value)}</title>`;
  const strandColor = strand => strand === '+' ? C.plus : strand === '-' ? C.minus : C.muted;
  function wrap(value, width=150) {
    const result=[]; let row='';
    for (const word of String(value ?? '').split(/\s+/)) {
      if (row && row.length+word.length+1>width) { result.push(row); row=''; }
      row+=(row?' ':'')+word;
    }
    if (row) result.push(row);
    return result;
  }
  function clipped(interval, region) {
    if (!Array.isArray(interval) || interval.length < 2) return null;
    const a=number(interval[0]), b=number(interval[1]);
    if (a===null || b===null || a>b || b<region.start || a>region.end) return null;
    return [Math.max(a,region.start),Math.min(b,region.end)];
  }
  function intervals(items, region) { return (Array.isArray(items)?items:[]).map(x=>clipped(x,region)).filter(Boolean); }
  function intervalToPixels(interval, region, x=0, width=1) {
    const c=clipped(interval,region); if (!c) return null;
    const length=region.end-region.start+1;
    return [x+(c[0]-region.start)/length*width,x+(c[1]+1-region.start)/length*width];
  }
  // A base coordinate denotes its centre; adding half a base aligns it to SVG edges.
  function pointToPixel(coordinate, region, x=0, width=1) {
    return x+(coordinate-region.start+0.5)/(region.end-region.start+1)*width;
  }
  function available(annotation) {
    const status=String(annotation?.status || '').toUpperCase();
    return ['FETCHED','LOADED','AVAILABLE','COMPLETE','CACHED','PRESENT','OK'].includes(status);
  }
  function niceMaximum(value) {
    if (value<=5) return Math.max(1,Math.ceil(value));
    const step=Math.pow(10,Math.floor(Math.log10(value))-1);
    return Math.ceil(value/(step*5))*step*5;
  }
  function normalize(region, annotation) {
    const start=number(region.start), end=number(region.end);
    if (start===null || end===null || start<1 || start>end || !Number.isInteger(start) || !Number.isInteger(end)) {
      throw Error('A saved region has invalid local coordinates.');
    }
    const r={...region,start,end,length:end-start+1,id:String(region.id ?? region.label ?? region.reference ?? ''),
      label:String(region.label ?? region.id ?? 'Region'),annotation:region.annotation || annotation || {},metrics:{...(region.metrics || {})}};
    r.sites=(Array.isArray(region.sites)?region.sites:[]).filter(s=>clipped([s.start,s.end],r)).map((s,i)=>({...s,
      id:String(s.id ?? s.name ?? i),name:String(s.name ?? s.id ?? `Site ${i+1}`),start:Number(s.start),end:Number(s.end)}));
    r.sectors=(Array.isArray(region.sectors)?region.sectors:[]).map(s=>{
      const interval=clipped([s.start,s.end],r), count=number(s.count);
      return interval && count!==null && count>=0 ? {start:interval[0],end:interval[1],count} : null;
    }).filter(Boolean);
    if (!r.sectors.length) {
      const step=Math.max(1000,Math.ceil(r.length/100/1000)*1000);
      for (let a=start;a<=end;a+=step) {
        const b=Math.min(a+step-1,end);
        r.sectors.push({start:a,end:b,count:r.sites.filter(s=>{const m=(s.start+s.end)/2;return m>=a&&m<b+1;}).length});
      }
    }
    r.genes=(Array.isArray(region.genes)?region.genes:[]).map(g=>({...g,
      id:String(g.id ?? g.label ?? ''),label:String(g.label ?? g.id ?? 'Gene'),
      spans:intervals(g.spans,r),exons:intervals(g.exons,r),
      transcripts:(Array.isArray(g.transcripts)?g.transcripts:[]).filter(t=>clipped([t.start,t.end],r)).map(t=>({...t,
        id:String(t.id ?? t.label ?? ''),label:String(t.label ?? t.id ?? 'Transcript'),
        start:Number(t.start),end:Number(t.end),exons:intervals(t.exons,r)}))
    }));
    return r;
  }
  function prepare(report, options={}) {
    if (!report || !Array.isArray(report.regions)) throw Error('Saved coverage data is unavailable.');
    const all=report.regions.map(r=>normalize(r,report.annotation));
    const regions=options.regionId && options.regionId!=='all' ? all.filter(r=>r.id===String(options.regionId)) : all;
    const maximum=all.reduce((max,r)=>r.sectors.reduce((v,s)=>Math.max(v,s.count),max),0);
    return {regions,allRegions:all,countMaximum:niceMaximum(maximum),mode:options.mode==='isoforms'?'isoforms':'overview',showTranscripts:options.showTranscripts!==false};
  }
  function track(spans,exons,strand,r,y,description) {
    const color=strandColor(strand), out=[`<g data-transcript-strand="${esc(strand || '?')}" stroke="${color}">`,title(description)];
    spans.forEach(span=>{
      const p=intervalToPixels(span,r,X,PW); if (!p) return;
      out.push(line(p[0],y,p[1],y,color,'stroke-width="1.3"'));
      if (strand==='+' || strand==='-') {
        for(let a=p[0]+15;a<p[1]-6;a+=52) {
          const tip=a+(strand==='+'?3:-3), tail=a+(strand==='+'?-3:3);
          out.push(`<path d="M${n(tail)},${n(y-3)}L${n(tip)},${n(y)}L${n(tail)},${n(y+3)}" fill="none" stroke-width="1.2"/>`);
        }
      }
    });
    exons.forEach(exon=>{const p=intervalToPixels(exon,r,X,PW);if(p)out.push(rect(p[0],y-4.5,p[1]-p[0],9,color));});
    out.push('</g>'); return out.join('');
  }
  function annotationRows(r,model) {
    if (!model.showTranscripts) return [];
    const rows=[];
    if (!available(r.annotation)) rows.push({type:'message',label:r.annotation?.status==='PARTIAL'?
      'Annotation incomplete — some transcript or exon models are unavailable.':'Annotation unavailable — transcript overlap is unknown.'});
    if (!r.genes.length) {
      if (available(r.annotation)) rows.push({type:'message',label:'No overlapping transcripts in the saved annotation.'});
      return rows;
    }
    r.genes.forEach(g=>{
      if (model.mode==='overview') rows.push({type:'gene',gene:g});
      else {
        rows.push({type:'gene-heading',gene:g});
        if (!g.transcripts.length) rows.push({type:'message',label:'Individual isoform models are unavailable for this gene.'});
        g.transcripts.forEach(t=>rows.push({type:'transcript',gene:g,transcript:t}));
      }
    });
    return rows;
  }
  function coordinateLabel(r) {
    if (r.coordinate_label) return String(r.coordinate_label);
    if (number(r.genomic_start)!==null && number(r.genomic_end)!==null) return `${r.chromosome || r.reference || 'Genomic'}:${fmt(r.genomic_start)}–${fmt(r.genomic_end)}`;
    return `${r.reference || 'Reference'}:${fmt(r.start)}–${fmt(r.end)} (local)`;
  }
  function sectorLabel(r) {
    const first=r.sectors[0], width=first.end-first.start+1;
    const standard=r.sectors.every((s,i)=>s.end-s.start+1===width || i===r.sectors.length-1);
    if (!standard) return 'Sites per sector';
    const unit=width>=1000 && width%1000===0 ? `${fmt(width/1000)} kb` : `${fmt(width)} bp`;
    return `Sites per ${unit}${r.sectors.at(-1).end-r.sectors.at(-1).start+1<width?' · last sector shorter':''}`;
  }
  function drawRegion(r,model,top) {
    const rows=annotationRows(r,model), tracksY=top+170, axisY=tracksY+Math.max(0,rows.length-1)*25+25,
      bottom=axisY+51, baseline=top+99, rugY=top+133, out=[];
    out.push(`<g data-region-id="${esc(r.id)}" data-coverage-region="${esc(r.id)}">`,title(`${r.label}; ${coordinateLabel(r)}`));
    out.push(text(24,top+19,short(`${r.label}${r.imager?' / '+r.imager:''}`,32),17,C.ink,'font-weight="700"'));
    const coordinateRows=wrap(coordinateLabel(r),45).slice(0,2);
    out.push(`<g>${title(coordinateLabel(r))}${coordinateRows.map((label,i)=>text(24,top+40+i*14,short(label,48),11,C.muted)).join('')}</g>`);
    out.push(text(24,top+76,`${fmt(r.metrics.sites ?? r.sites.length)} sites · ${fmt(r.metrics.primaries ?? r.sites.length)} primaries`,12));
    out.push(text(24,top+95,`Left : right  ${fmt(r.metrics.left_half)} : ${fmt(r.metrics.right_half)}`,12));
    const mean=number(r.metrics.mean_offset_bp);
    out.push(text(24,top+114,`Mean local offset  ${mean===null?'unknown':(mean>0?'+':'')+fmt(mean)+' bp'}`,12));
    out.push(`<g>${title(`Unprobed arm intervals: largest internal gap ${fmt(r.metrics.largest_gap_bp)} bp; left edge ${fmt(r.metrics.left_edge_gap_bp)} bp; right edge ${fmt(r.metrics.right_edge_gap_bp)} bp.`)}`,
      text(24,top+133,`Largest gap  ${fmt(r.metrics.largest_gap_bp)} bp`,12),'</g>');
    const bound=r.bound_strand;
    out.push(text(24,top+152,bound==='+'||bound==='-'?`DNA target ${bound} · complementary RNA ${bound}`:'DNA/RNA orientation unknown',10,C.muted));
    r.sectors.forEach((s,i)=>{
      const p=intervalToPixels([s.start,s.end],r,X,PW);
      if(i%2===0) out.push(rect(p[0],top+22,p[1]-p[0],axisY-top-22,C.pale));
    });
    [0,model.countMaximum/2,model.countMaximum].forEach(v=>{
      const y=baseline-v/model.countMaximum*BAR_H;
      out.push(line(X,y,X+PW,y,C.grid,'stroke-width="0.7"'));
      out.push(text(X-10,y+3,fmt(v),10,C.muted,'text-anchor="end"'));
    });
    out.push(text(X,top+11,sectorLabel(r),11,C.muted));
    r.sectors.forEach(s=>{
      const p=intervalToPixels([s.start,s.end],r,X,PW), full=p[1]-p[0], gap=Math.min(8,full*0.12), height=s.count/model.countMaximum*BAR_H;
      out.push(`<g class="pc-sector">${title(`${s.count} exact sites; local ${fmt(s.start)}–${fmt(s.end)}; ${fmt(s.end-s.start+1)} bp`)}`,
        rect(p[0]+gap/2,baseline-height,Math.max(0.2,full-gap),height,C.blue));
      if (full>=23) out.push(text((p[0]+p[1])/2,baseline-height-5,fmt(s.count),11,C.blue,'text-anchor="middle"'));
      out.push('</g>');
    });
    out.push(line(X,rugY,X+PW,rugY,C.grid),text(X,rugY+21,'Exact-site midpoints',10,C.muted));
    out.push(`<g stroke="${C.blue}" stroke-width="1" fill="none">`);
    r.sites.forEach(s=>{
      const x=pointToPixel((s.start+s.end)/2,r,X,PW);
      if (x<X || x>X+PW) return;
      let genomic='';
      if (r.mapping_orientation==='+' && number(r.genomic_start)!==null) {
        genomic=`; ${r.chromosome || 'genomic'}:${fmt(r.genomic_start+s.start-r.start)}–${fmt(r.genomic_start+s.end-r.start)}`;
      } else if (r.mapping_orientation==='-' && number(r.genomic_end)!==null) {
        genomic=`; ${r.chromosome || 'genomic'}:${fmt(r.genomic_end-(s.end-r.start))}–${fmt(r.genomic_end-(s.start-r.start))}`;
      }
      out.push(`<path data-site-id="${esc(s.id)}" data-probe-name="${esc(s.name)}" tabindex="0" d="M${n(x)},${rugY-5}v10">${title(`${s.name}; local ${fmt(s.start)}–${fmt(s.end)}${genomic}`)}</path>`);
    });
    out.push('</g>');
    out.push(line(X+PW/2,top+19,X+PW/2,axisY,C.muted,'stroke-width="0.8" stroke-dasharray="3 4"'));
    if (mean!==null && r.sites.length) {
      const x=pointToPixel((r.start+r.end)/2+mean,r,X,PW);
      if(x>=X && x<=X+PW)out.push(`<path data-coverage-mean="${esc(r.id)}" d="M${n(x)},${rugY+8}l-4,4l4,4l4,-4Z" fill="${C.ink}">${title(`Equal-weight reference-local mean offset: ${fmt(mean)} bp from region midpoint. Not a measured 3D error.`)}</path>`);
    }
    if (!r.sites.length) out.push(text(X+PW/2,top+67,'No selected binding sites',13,C.muted,'text-anchor="middle"'));
    rows.forEach((row,i)=>{
      const y=tracksY+i*25;
      if(row.type==='message') {out.push(text(X,y+4,row.label,12,C.muted));return;}
      const g=row.gene, color=strandColor(g.strand);
      if(row.type==='gene-heading') {
        out.push(text(24,y+4,short(`${g.label} (${g.strand || '?'}; ${g.transcripts.length} isoforms)`,43),12,color,'font-weight="700"'));
      } else if(row.type==='gene') {
        const label=`${g.label} (${g.strand || '?'}; ${g.transcripts.length} isoform${g.transcripts.length===1?'':'s'})`;
        out.push(`<g>${title(label)}${text(24,y+4,short(label,43),12,color,'font-weight="700"')}</g>`);
        out.push(track(g.spans,g.exons,g.strand,r,y,`${label}. Collapsed union of transcript spans and exons; not a single isoform.`));
      } else {
        const t=row.transcript, label=t.label===t.id?t.id:`${t.label} · ${t.id}`;
        out.push(`<g data-transcript-id="${esc(t.id)}">${title(label)}`,text(34,y+4,short(label,44),11,C.ink),
          track([[t.start,t.end]],t.exons,t.strand,r,y,`${label}; strand ${t.strand || '?'}; local ${fmt(t.start)}–${fmt(t.end)}. ${available(r.annotation)?'Exon blocks and connecting introns of this isoform.':'Saved exon blocks and transcript span; the exon model may be incomplete.'}`),'</g>');
      }
    });
    out.push(line(X,axisY,X+PW,axisY,C.muted,'stroke-width="0.8"'));
    const kb=r.length>=1000, divisions=Math.min(10,Math.max(1,Math.floor(r.length/100)));
    for(let i=0;i<=divisions;i++) {
      const x=X+PW*i/divisions, offset=r.length*i/divisions;
      out.push(line(x,axisY,x,axisY+4,C.muted),text(x,axisY+17,fmt(offset/(kb?1000:1)),10,C.muted,`text-anchor="${i===0?'start':i===divisions?'end':'middle'}"`));
    }
    out.push(text(X,axisY+35,`Position within region (${kb?'kb':'bp'}) · shared axis for sites and transcripts`,11,C.muted));
    const introns=number(r.metrics.complete_intron_sites), annotationKnown=available(r.annotation);
    out.push(text(24,axisY+16,`Complete-arm intron overlaps: ${annotationKnown&&introns!==null?fmt(introns):'unknown'}`,11,C.muted));
    out.push(line(24,bottom,W-24,bottom,C.grid),'</g>');
    return {svg:out.join(''),bottom};
  }
  function shell(body,height,label,description,extra='') {
    return `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${W} ${height}" width="${W}" height="${height}" role="img" aria-label="${esc(label)}" data-probe-coverage-svg="true" ${extra}>`+
      title(label)+`<desc>${esc(description)}</desc><style>[data-probe-coverage-svg] [data-site-id]{cursor:pointer;stroke-linecap:round}[data-probe-coverage-svg] [data-site-id]:hover,[data-probe-coverage-svg] [data-site-id]:focus{stroke:#182637;stroke-width:3;outline:none}.pc-sector:hover{opacity:.8}</style>`+
      rect(0,0,W,height,C.white)+body+'</svg>';
  }
  function render(report,options={}) {
    try {
      const model=prepare(report,options), mode=model.mode, label=`${report.title || 'Probe coverage'}${mode==='isoforms'?' · individual transcript isoforms':''}`;
      const incomplete=model.regions.some(r=>!available(r.annotation));
      const out=[text(24,38,short(label,83),25,C.ink,'font-weight="700"'),
        text(24,63,`${model.regions.length} region${model.regions.length===1?'':'s'} · saved selected binding sites · ${String(report.order_status || 'DRAFT')}${report.annotation?.source?' · Annotation: '+short(report.annotation.source,85):''}`,12,C.muted),
        text(24,88,'Bars: site counts. Ticks: exact-site midpoints. Diamond: mean sequence position. Dashed line: region midpoint.',12),
        text(24,110,model.showTranscripts?(mode==='isoforms'?(incomplete?'Boxes: saved exons. Lines: transcript spans; exon/intron identity may be incomplete. Arrows: transcript direction.':'Boxes: exons of each individual isoform. Connecting lines: introns. Arrows: transcript direction.'):'Boxes: exon union. Arrows: transcript direction. Each gene combines its overlapping isoforms.'):'Transcript tracks are hidden. Coverage counts describe saved binding sites.',12),
        text(24,135,'Forward (+) →',12,C.plus,'font-weight="700"'),text(165,135,'← Reverse (−)',12,C.minus,'font-weight="700"'),
        text(310,135,'Reference-local transcript direction; colors do not measure RNA risk or expression.',12,C.muted)];
      let y=162;
      if (!model.regions.length) {out.push(text(24,y+30,'No saved regions match this view.',15,C.muted));y+=70;}
      model.regions.forEach(r=>{const block=drawRegion(r,model,y);out.push(block.svg);y=block.bottom+26;});
      const cross=report.cross_primary_status || 'UNKNOWN — no saved cross-primary status';
      const suppliedWarnings=(Array.isArray(report.warnings)?report.warnings:[]).filter(w=>
        !/^(Coverage is a read-only view of saved selections|Genomic offsets describe sequence positions|Reference-local mean offsets describe sequence positions|Saved annotation describes reference transcript models|Cross-primary interactions within and between regions:)/.test(String(w)));
      const notes=[`Count axes share a common maximum of ${fmt(model.countMaximum)} sites across all saved regions. Region lengths and sector widths may differ.`,
        'Transcript tracks are clipped to each region. Collapsed exon unions are not single isoform models.',
        'RNA matching the bound DNA strand is complementary to the arm. Intron overlaps indicate potential unspliced RNA.',
        'Saved annotation does not establish expression, accessibility or binding. This view does not rerun design or QC.',
        'No annotated overlap is not an RNA-free guarantee. Reference-local mean offsets are not measured 3D centroid errors.',
        `${String(report.order_status || 'DRAFT')} · Cross-primary interactions: ${String(cross)}`,
        ...new Set(suppliedWarnings.map(String))];
      notes.forEach(note=>wrap(note,168).forEach(row=>{out.push(text(24,y+14,row,11,C.muted));y+=18;}));
      const description='Scientific coverage view using saved selected sites and local 1-based inclusive coordinates. '+notes.join(' ');
      return shell(out.join(''),Math.ceil(y+18),label,description,`data-coverage-mode="${mode}" data-count-maximum="${model.countMaximum}"`);
    } catch(error) {
      const message=error instanceof Error?error.message:'Saved coverage data could not be displayed.';
      return shell(text(24,40,'Coverage unavailable',22,C.ink,'font-weight="700"')+text(24,70,message,13,C.muted),112,'Coverage unavailable',message,'data-coverage-error="true"');
    }
  }
  return Object.freeze({render,prepare,intervalToPixels,pointToPixel});
});
