/* Sequence-derived imaging-oligo / primary-docking visualization.
 * Display layer only. Does not mutate input sequences or run thermodynamic QC.
 * Both stored sequences are 5'->3'; only the imager's display row is reversed.
 * No external packages, fonts, images, or network requests are required.
 */
(function (root, factory) {
  'use strict';
  const api = factory();
  if (typeof module === 'object' && module.exports) module.exports = api;
  else root.ImagerBinding = api;
})(typeof globalThis !== 'undefined' ? globalThis : this, function () {
  'use strict';
  const C = Object.freeze({A:'T', T:'A', C:'G', G:'C'});
  const escape = x => String(x ?? '').replace(/[&<>"']/g, c =>
    ({'&':'&amp;', '<':'&lt;', '>':'&gt;', '"':'&quot;', "'":'&#39;'}[c]));
  const reverse = s => [...s].reverse().join('');
  const rc = s => reverse(s).split('').map(c => C[c]).join('');
  function dna(value, name) {
    if (typeof value !== 'string') throw Error(name + ' sequence is missing.');
    const seq = value.replace(/\s/g, '').toUpperCase();
    if (!/^[ACGT]{8,40}$/.test(seq)) throw Error(name + ' must contain 8-40 A/C/G/T bases.');
    return seq;
  }
  function readPair(im) {
    if (!im || typeof im !== 'object') throw Error('Choose an imaging oligo first.');
    const sequence = dna(im.sequence ?? im.imager_sequence_5to3, 'Imager');
    const docking = dna(im.docking ?? im.docking_sequence_5to3, 'Docking');
    if (rc(sequence) !== docking) throw Error('Imager and docking sequences are not reverse complements. No paired diagram is drawn.');
    for (const [k, val] of [['imager_sequence_5to3',sequence],['docking_sequence_5to3',docking]]) {
      if (im[k] != null && dna(im[k], k) !== val) throw Error('Conflicting saved sequence fields.');
    }
    const display = reverse(sequence);
    const columns = [...docking].map((base, i) => ({
      index:i, dockingPosition:i+1, dockingBase:base,
      imagerPosition:sequence.length-i, imagerBase:display[i],
      complementary:C[base] === display[i]
    }));
    return {id:im.id || 'Imaging oligo', dye:im.dye || 'not specified',
      sequence, docking, imagerDisplay:display, columns};
  }
  function checkPrimary(primary, pair) {
    if (!primary) return null;
    const {arm, prefix, suffix, full} = primary;
    if (typeof arm !== 'string' || !/^[ACGT]+$/.test(arm) ||
        prefix !== pair.docking+'TT' || suffix !== 'TT'+pair.docking ||
        full !== prefix+arm+suffix) {
      throw Error('Saved primary does not match docking-TT-arm-TT-docking. Check the saved sequence; the diagram will not invent an assembly.');
    }
    return primary;
  }
  const text = (x,y,label,fill='#dce7f5',size=13,extra='') =>
    `<text x="${x}" y="${y}" fill="${fill}" font-family="system-ui,-apple-system,sans-serif" font-size="${size}" ${extra}>${escape(label)}</text>`;
  const rect = (x,y,w,h,fill,stroke,extra='') =>
    `<rect x="${x}" y="${y}" width="${w}" height="${h}" rx="5" fill="${fill}" stroke="${stroke}" ${extra}/>`;
  function duplex(pair) {
    const step=24, n=pair.docking.length, w=Math.max(590, n*step+112), h=238;
    const left=(w-n*step)/2, y1=68, y2=157;
    const out=[`<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${w} ${h}" width="${w}" height="${h}" role="img" aria-label="${escape(pair.id)} imaging oligo base-paired antiparallel to its docking site. Imager drawn 3-prime to 5-prime; primary docking drawn 5-prime to 3-prime." style="min-width:${Math.round(w*.8)}px;max-width:${Math.round(w*1.38)}px" data-binding-svg="duplex" data-imager-id="${escape(pair.id)}" data-imager-5to3="${pair.sequence}" data-docking-5to3="${pair.docking}">`,
      '<title>Imaging oligo and primary docking site: actual sequence pairing</title>',
      `<desc>Imaging oligo for ordering, 5-prime ${pair.sequence} 3-prime. Primary docking sequence, 5-prime ${pair.docking} 3-prime. The top row is reversed for antiparallel display only. This is not a folding or occupancy prediction.</desc>`,
      '<style>.binding-column:focus{outline:none}.binding-column .column-highlight{opacity:0}.binding-column:hover .column-highlight,.binding-column:focus .column-highlight{opacity:1}.binding-column:hover .base-cell,.binding-column:focus .base-cell{stroke:#ffffff;stroke-width:2}</style>',
      rect(1,1,w-2,h-2,'#101925','#324356'),
      text(22,28,'FLUORESCENT IMAGING OLIGO  |  3\u2032 to 5\u2032 in this view','#ffd18c',12,'font-weight="700"'),
      text(left-26,y1+5,'3\u2032','#ffd18c',16),text(left+n*step+10,y1+5,'5\u2032','#ffd18c',16),
      text(left-26,y2+5,'5\u2032','#cdb6ff',16),text(left+n*step+10,y2+5,'3\u2032','#cdb6ff',16)];
    pair.columns.forEach(col => {
      const x=left+col.index*step+step/2;
      const label=`Docking base ${col.dockingPosition}: ${col.dockingBase} pairs with imager base ${col.imagerPosition}: ${col.imagerBase}. Imager positions use its stored 5-prime-to-3-prime sequence.`;
      out.push(`<g class="binding-column" data-base-index="${col.index}" data-docking-base="${col.dockingBase}" data-imager-base="${col.imagerBase}" data-docking-position="${col.dockingPosition}" data-imager-position="${col.imagerPosition}" data-pair-description="${escape(label)}" tabindex="0" role="group" aria-label="${escape(label)}"><title>${escape(label)}</title>`,
        rect(x-12,y1-25,24,125,'#23344b','#819db3','class="column-highlight"'),
        `<line x1="${x}" y1="${y1+20}" x2="${x}" y2="${y2-20}" stroke="#7599a5" stroke-width="1.5"/>`,
        rect(x-10,y1-17,20,29,'#392f24','#8e744d','class="base-cell"'),
        rect(x-10,y2-17,20,29,'#2c2442','#84709e','class="base-cell"'),
        text(x,y1+4,col.imagerBase,'#ffe0a3',17,'text-anchor="middle" style="font-family:ui-monospace,SFMono-Regular,Consolas,monospace;font-weight:700" data-letter="imager"'),
        text(x,y2+4,col.dockingBase,'#e2d3ff',17,'text-anchor="middle" style="font-family:ui-monospace,SFMono-Regular,Consolas,monospace;font-weight:700" data-letter="docking"'),
        text(x,y2+31,col.dockingPosition,'#96aabd',10,'text-anchor="middle"'),'</g>');
    });
    out.push(text(22,214,'DOCKING SITE ON THE PRIMARY  |  5\u2032 to 3\u2032','#cdb6ff',12,'font-weight="700"'),'</svg>');
    return out.join('');
  }
  function overview(pair, primary) {
    const w=730, h=203, top=57, bottom=129;
    const out=[`<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${w} ${h}" width="${w}" height="${h}" role="img" aria-label="Two separate copies of ${escape(pair.id)} on the two identical docking sites of one primary oligo. The middle DNA-binding arm is separate." data-binding-svg="overview">`,
      '<title>Two imaging molecules on one primary probe, schematic not to scale</title>',
      '<desc>The continuous primary runs from 5-prime on the left to 3-prime on the right: docking site, TT spacer, DNA-binding arm, TT spacer, docking site. Each imaging molecule is antiparallel to its docking site. Two possible occupied sites are illustrated, not measured occupancy. Dye attachment end is unspecified.</desc>',
      rect(1,1,w-2,h-2,'#101925','#324356'),
      text(20,25,'WHERE THE IMAGER BINDS  |  one primary shown 5\u2032 to 3\u2032','#dce7f5',12,'font-weight="700"'),
      `<line x1="42" y1="${bottom}" x2="688" y2="${bottom}" stroke="#6a7e96" stroke-width="2"/>`,
      text(18,bottom+5,'5\u2032','#dce7f5',15),text(697,bottom+5,'3\u2032','#dce7f5',15)];
    [[60,'1'],[540,'2']].forEach(([x,number])=>{
      out.push(rect(x,top-14,130,30,'#392f24','#ae8850'),text(x+65,top+5,'Imager copy '+number,'#ffe0a3',12,'text-anchor="middle"'),
        text(x-23,top+5,'3\u2032','#ffd18c',12),text(x+137,top+5,'5\u2032','#ffd18c',12),
        `<line x1="${x+117}" y1="${top-24}" x2="${x+13}" y2="${top-24}" stroke="#ffd18c"/><path d="M ${x+13} ${top-24} l 6 -3 v 6 z" fill="#ffd18c"/>`,
        rect(x,bottom-16,130,32,'#2c2442','#9780b6'),text(x+65,bottom+5,'Docking site '+number,'#e2d3ff',12,'text-anchor="middle"'));
      for(let j=0;j<7;j++)out.push(`<line x1="${x+17+j*16}" y1="${top+17}" x2="${x+17+j*16}" y2="${bottom-18}" stroke="#7599a5" stroke-dasharray="3 3"/>`);
    });
    out.push(rect(209,bottom-16,35,32,'#26303b','#8190a0'),text(226,bottom+5,'TT','#dce7f5',12,'text-anchor="middle"'),
      rect(262,bottom-20,207,40,'#183b3e','#65d6c4'),text(365,bottom+1,'DNA-binding arm','#a9f1e7',14,'text-anchor="middle" font-weight="700"'),
      text(365,bottom+15,primary?primary.arm.length+' bases, from this candidate':'Designed against your target DNA','#a9f1e7',10,'text-anchor="middle"'),
      rect(487,bottom-16,35,32,'#26303b','#8190a0'),text(504,bottom+5,'TT','#dce7f5',12,'text-anchor="middle"'),
      text(20,181,'Purple = docking handles; teal = target-binding arm; TT = spacers. Not to scale.','#a6bacd',11),'</svg>');
    return out.join('');
  }
  function render(im, options={}) {
    try {
      const pair=readPair(im), primary=checkPrimary(options.primary, pair);
      const open=options.showOverview ? ' open' : '';
      return `<figure class="binding-figure" data-binding-figure data-binding-imager="${escape(pair.id)}">`+
        `<figcaption class="binding-heading"><strong>${escape(pair.id)}: imager bound to its docking site</strong><span class="binding-badge">${pair.docking.length} paired bases</span></figcaption>`+
        `<div class="binding-scroll">${duplex(pair)}</div>`+
        `<p class="binding-hover-help" role="status">Point at or tab to a base pair to highlight it and see both positions.</p>`+
        `<p class="binding-note"><strong>Display orientation only:</strong> the gold imager row runs 3\u2032 to 5\u2032 so paired bases line up. Copy/order from the labeled 5\u2032-to-3\u2032 fields, not the reversed gold row.</p>`+
        `<details class="binding-overview"${open}><summary>See both docking sites on one complete primary</summary><div class="binding-scroll">${overview(pair,primary)}</div><p class="binding-note">Both docking sites carry the sequence shown above. Two separate imager molecules are illustrated. Dye / note: <strong>${escape(pair.dye)}</strong>. Dye attachment end is not inferred. This is a pairing guide, not a folding, binding-strength, or occupancy prediction.</p></details>`+
        '<div class="binding-actions"><button type="button" data-binding-export="duplex">Save pairing SVG</button><button type="button" data-binding-export="overview">Save primary schematic SVG</button></div></figure>';
    } catch (err) {
      return `<div class="notice danger" role="alert" data-binding-error>Binding diagram unavailable: ${escape(err.message)}</div>`;
    }
  }
  return Object.freeze({pairData:readPair, render});
});
