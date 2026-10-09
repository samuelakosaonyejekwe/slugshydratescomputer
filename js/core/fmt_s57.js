// IHO S-57 electronic navigational chart reader (base cells, ".000"): an ISO/IEC 8211 file of feature and vector records.
//
// Read:   the ISO 8211 container (leaders, directories, the field descriptions of the data descriptive record with
//         their subfield labels and format controls: binary b1x / b2x integers, A and A(n) text, I and R numbers, B(n) bit
//         strings, repeating groups); the S-57 edition 3 records: DSID / DSPM (coordinate and sounding multiplication
//         factors, compilation scale, datums), vector records (isolated and connected nodes, edges with their SG2D
//         points and VRPT end nodes, SG3D soundings) and feature records (object class, primitive, ATTF attributes, FSPT
//         pointers with orientation). Line features come as polylines of their edges chained end to end, area features
//         as their closed boundary rings, soundings as points (longitude, latitude, elevation = -depth in metres),
//         other point features as markers; the layer of each is the object-class acronym (DEPCNT, PIPSOL, CBLSUB …).
// Not read: update cells (.001 …: they are not applied to their base cell), the Unicode national attributes (NATF), feature
//         relations and collections, masked-edge flags, files in the ASCII ("level 1") implementation as far as they use
//         formats other than those named, and S-101 data sets.
// File content is untrusted: record lengths, directory entries and counts are bounds-checked and capped.

function fail(msg) { const e = new Error(msg); e.user = true; throw e; }
// object classes of the S-57 catalogue, in code order
const OBJ = ('ADMARE AIRARE ACHBRT ACHARE BCNCAR BCNISD BCNLAT BCNSAW BCNSPP BERTHS BRIDGE BUISGL BUAARE BOYCAR BOYINB BOYISD BOYLAT BOYSAW BOYSPP CBLARE CBLOHD CBLSUB CANALS CANBNK CTSARE CAUSWY CTNARE CHKPNT CGUSTA COALNE CONZNE COSARE CTRPNT CONVYR CRANES CURENT CUSZNE DAMCON DAYMAR DWRTCL DWRTPT DEPARE DEPCNT DISMAR DOCARE DRGARE DRYDOC DMPGRD DYKCON EXEZNE FAIRWY FNCLNE FERYRT FSHZNE FSHFAC FSHGRD FLODOC FOGSIG FORSTC FRPARE GATCON GRIDRN HRBARE HRBFAC HULKES ICEARE ICNARE ISTZNE LAKARE LAKSHR LNDARE LNDELV LNDRGN LNDMRK LIGHTS LITFLT LITVES LOCMAG LOKBSN LOGPON MAGVAR MARCUL MIPARE MORFAC NAVLNE OBSTRN OFSPLF OSPARE OILBAR PILPNT PILBOP PIPARE PIPOHD PIPSOL PONTON PRCARE PRDARE PYLONS RADLNE RADRNG RADRFL RADSTA RTPBCN RDOCAL RDOSTA RAILWY RAPIDS RCRTCL RECTRC RCTLPT RSCSTA RESARE RETRFL RIVERS RIVBNK ROADWY RUNWAY SNDWAV SEAARE SPLARE SBDARE SLCONS SISTAT SISTAW SILTNK SLOTOP SLOGRD SMCFAC SOUNDG SPRING SQUARE STSLNE SUBTLN SWPARE TESARE TS_PRH TS_PNH TS_PAD TS_TIS T_HMON T_NHMN T_TIMS TIDEWY TOPMAR TSELNE TSSBND TSSCRS TSSLPT TSSRON TSEZNE TUNNEL TWRTPT UWTROC UNSARE VEGATN WATTUR WATFAL WEDKLP WRECKS').split(' ');
const META = { 300: 'M_ACCY', 301: 'M_CSCL', 302: 'M_COVR', 303: 'M_HDAT', 304: 'M_HOPA', 305: 'M_NPUB', 306: 'M_NSYS', 307: 'M_PROD', 308: 'M_QUAL', 309: 'M_SDAT', 310: 'M_SREL', 311: 'M_UNIT', 312: 'M_VDAT', 400: 'C_AGGR', 401: 'C_ASSO', 402: 'C_STAC' };
const ATT = { 3: 'BURDEP', 87: 'DRVAL1', 88: 'DRVAL2', 116: 'OBJNAM', 174: 'VALDCO', 179: 'VALSOU' };
const className = (c) => OBJ[c - 1] || META[c] || `class ${c}`;

/** True for an ISO 8211 file whose first record is a data descriptive record. */
export const isS57 = (u8) => u8.length > 24 && /^\d{5}[ 123]L/.test(String.fromCharCode(...u8.subarray(0, 7)));

/** Format controls "(b11,2b24,A(2),*…)" -> flat list of [kind, width]: b unsigned / s signed binary, A text, B bits. */
function formats(ctl) {
  const out = [];
  const walk = (s, depth) => {
    let i = 0;
    while (i < s.length && out.length < 4096) {
      if (s[i] === ',') { i++; continue; }
      let rep = 0;
      while (s[i] >= '0' && s[i] <= '9') rep = rep * 10 + +s[i++];
      rep = Math.min(rep || 1, 1024);
      if (s[i] === '(') { let d = 1, j = i + 1; while (j < s.length && d) { if (s[j] === '(') d++; else if (s[j] === ')') d--; j++; } if (depth < 8) for (let k = 0; k < rep; k++) walk(s.slice(i + 1, j - 1), depth + 1); i = j; continue; }
      const kind = s[i++];
      let item;
      if (kind === 'b') { const sign = s[i++], w = +s[i++]; item = [sign === '2' ? 's' : 'b', w]; }
      else { let w = 0; if (s[i] === '(') { const j = s.indexOf(')', i); w = +s.slice(i + 1, j < 0 ? s.length : j); i = j < 0 ? s.length : j + 1; } item = [kind === 'B' ? 'B' : 'A', kind === 'B' ? w >> 3 : w]; }
      for (let k = 0; k < rep; k++) out.push(item);
    }
  };
  walk(ctl.replace(/^\(|\)$/g, ''), 0);
  return out;
}

/**
 * Parse an S-57 base cell. opts: { maxRecords = 2e6 }. Returns
 * { polylines: [{ x (longitude), y (latitude), closed, layer (class acronym), type: 'line' | 'area', attrs }],
 *   points: flat lon, lat, -depth of the soundings, markers: [{ x, y, layer }] (point features),
 *   classes: { acronym: feature count }, scale, comf, somf, hdat, vdat, counts, warnings }
 */
export function parseS57(u8, opts = {}) {
  if (!isS57(u8)) fail('Not an S-57 / ISO 8211 file (the record leader is missing).');
  const N = u8.length, dv = new DataView(u8.buffer, u8.byteOffset, N), maxRecords = opts.maxRecords ?? 2e6, warnings = [];
  const txt = (p, n) => { let s = ''; for (let k = 0; k < n; k++) s += String.fromCharCode(u8[p + k]); return s; }, bad = (w) => fail(`The S-57 file is truncated or corrupt (${w}).`);
  const defs = new Map();                                   // tag -> { labels, fmt, rep: index of the first repeating subfield or -1 }
  let comf = 1e7, somf = 10, scale = 0, hdat = 0, vdat = 0, nRec = 0, aall = 0;
  const nodes = new Map(), edges = new Map(), feats = [], soundings = new Map();
  /** Rows of a field: one object per repetition, by subfield label. */
  const rows = (tag, p, end) => {
    const d = defs.get(tag), out = [];
    if (!d || !d.fmt.length) return out;
    const fixed = d.rep < 0 ? [] : d.fmt.slice(0, d.rep), fl = d.rep < 0 ? [] : d.labels.slice(0, d.rep), rf = d.rep < 0 ? d.fmt : d.fmt.slice(d.rep), rl = d.rep < 0 ? d.labels : d.labels.slice(d.rep);
    const take = (fm, lb, row) => {
      for (let k = 0; k < lb.length && k < fm.length; k++) {
        const [kind, w] = fm[k];
        let v;
        if (kind === 'A') { if (w) { if (p + w > end) return false; v = txt(p, w); p += w; } else { let q = p; while (q < end && u8[q] !== 0x1f && u8[q] !== 0x1e) q++; v = txt(p, q - p); p = q < end ? q + 1 : q; } }
        else { if (!(w > 0) || p + w > end) return false; if (kind === 'B') v = u8.subarray(p, p + w); else if (w === 1) v = kind === 's' ? dv.getInt8(p) : u8[p]; else if (w === 2) v = kind === 's' ? dv.getInt16(p, true) : dv.getUint16(p, true); else if (w === 4) v = kind === 's' ? dv.getInt32(p, true) : dv.getUint32(p, true); else v = 0; p += w; }
        row[lb[k]] = v;
      }
      return true;
    };
    const head = {};
    if (fixed.length && !take(fixed, fl, head)) return out;
    if (!rl.length) { out.push(head); return out; }
    // rows of fixed binary width are counted off (their bytes may equal the field terminator); text rows end at the terminator
    const size = rf.slice(0, rl.length).every((q) => q[1] > 0) ? rf.slice(0, rl.length).reduce((t, q) => t + q[1], 0) : 0;
    while ((size > 1 ? p + size <= end : p < end && u8[p] !== 0x1e) && out.length < 4e6) { const row = { ...head }; if (!take(rf, rl, row)) break; out.push(row); }
    if (!out.length && fixed.length) out.push(head);
    return out;
  };
  const name = (b) => (b && b.length >= 5 ? b[0] * 4294967296 + ((b[1] | (b[2] << 8) | (b[3] << 16) | (b[4] << 24)) >>> 0) : -1);
  let lenSize = 0, posSize = 0, tagSize = 0;
  for (let p = 0; p + 24 <= N; ) {
    let rl = +txt(p, 5);
    const lid = txt(p + 6, 1), base = +txt(p + 12, 5);
    if (lid === 'L' || !lenSize || txt(p + 20, 1).trim()) { lenSize = +txt(p + 20, 1) || lenSize; posSize = +txt(p + 21, 1) || posSize; tagSize = +txt(p + 23, 1) || tagSize; }
    if (!(base >= 24) || !(lenSize > 0 && posSize > 0 && tagSize > 0)) bad('record leader');
    if (++nRec > maxRecords) fail('The S-57 file holds too many records.');
    // directory
    const dir = [];
    let q = p + 24, last = 0;
    for (; q + tagSize + lenSize + posSize <= p + base && u8[q] !== 0x1e && dir.length < 1e4; q += tagSize + lenSize + posSize) { const len = +txt(q + tagSize, lenSize), pos = +txt(q + tagSize + lenSize, posSize); if (!(len >= 0 && pos >= 0)) bad('directory'); dir.push([txt(q, tagSize), p + base + pos, len]); last = Math.max(last, base + pos + len); }
    if (!rl) rl = last;                                     // a zero length stands for "as long as the directory says"
    if (!(rl >= 24) || p + rl > N) { warnings.push('The file ends inside a record; the records read so far are kept.'); break; }
    if (lid === 'L') {
      const fcl = +txt(p + 10, 2) || 9;
      for (const [tag, at, len] of dir) {
        const parts = txt(at + fcl, Math.max(0, len - fcl)).replace(/\x1e$/, '').split('\x1f'), labels = (parts[1] || '').split('!'), rep = labels.findIndex((l) => l[0] === '*');
        defs.set(tag, { labels: labels.map((l) => l.replace(/^\*/, '')).filter(Boolean), fmt: formats(parts[2] || ''), rep: rep < 0 && (parts[1] || '')[0] === '*' ? 0 : rep });
      }
    } else {
      let cur = null;
      for (const [tag, at, len] of dir) {
        const end = Math.min(at + len, p + rl);
        if (at > end) continue;
        if (tag === 'DSPM') { const r = rows(tag, at, end)[0] || {}; if (r.COMF > 0) comf = r.COMF; if (r.SOMF > 0) somf = r.SOMF; scale = r.CSCL || 0; hdat = r.HDAT || 0; vdat = r.VDAT || 0; }
        else if (tag === 'DSSI') { const r = rows(tag, at, end)[0] || {}; aall = r.AALL || 0; }
        else if (tag === 'VRID') { const r = rows(tag, at, end)[0] || {}; cur = { kind: 'v', rcnm: r.RCNM, id: r.RCNM * 4294967296 + (r.RCID >>> 0), pts: [], ends: [] }; if (r.RCNM === 130) edges.set(cur.id, cur); else nodes.set(cur.id, cur); }
        else if (tag === 'FRID') { const r = rows(tag, at, end)[0] || {}; cur = { kind: 'f', prim: r.PRIM, objl: r.OBJL, attrs: {}, refs: [] }; feats.push(cur); }
        else if (!cur) continue;
        else if (tag === 'SG2D' && cur.kind === 'v') for (const r of rows(tag, at, end)) cur.pts.push([r.XCOO / comf, r.YCOO / comf]);
        else if (tag === 'SG3D' && cur.kind === 'v') { const l = []; for (const r of rows(tag, at, end)) l.push([r.XCOO / comf, r.YCOO / comf, r.VE3D / somf]); soundings.set(cur.id, l); }
        else if (tag === 'VRPT' && cur.kind === 'v') for (const r of rows(tag, at, end)) cur.ends.push([name(r.NAME), r.TOPI]);
        else if (tag === 'ATTF' && cur.kind === 'f') for (const r of rows(tag, at, end)) { const k = ATT[r.ATTL]; if (k) { const v = String(r.ATVL ?? '').trim(); cur.attrs[k] = v !== '' && Number.isFinite(+v) && k !== 'OBJNAM' ? +v : v; } }
        else if (tag === 'FSPT' && cur.kind === 'f') for (const r of rows(tag, at, end)) cur.refs.push([name(r.NAME), r.ORNT, r.USAG]);
      }
    }
    p += rl;
  }
  if (!defs.size) bad('no data descriptive record');
  if (!feats.length && !edges.size) fail('The S-57 file holds no feature or vector records (an update cell or a catalogue file, not a base cell).');
  // edge -> points from its start node through its own points to its end node
  const nodePt = (id) => { const n = nodes.get(id); return n && n.pts.length ? n.pts[0] : null; };
  const edgePts = (id) => { const e = edges.get(id); if (!e) return null; if (!e.line) { const a = e.ends.find((q) => q[1] === 1), b = e.ends.find((q) => q[1] === 2), A = a && nodePt(a[0]), B = b && nodePt(b[0]); e.line = [...(A ? [A] : []), ...e.pts, ...(B ? [B] : [])]; } return e.line; };
  const polylines = [], points = [], markers = [], classes = {}, same = (a, b) => a[0] === b[0] && a[1] === b[1];
  let nVert = 0;
  for (const f of feats) {
    const layer = className(f.objl);
    classes[layer] = (classes[layer] || 0) + 1;
    if (f.prim === 1) {
      for (const [id] of f.refs) { const s = soundings.get(id), q = nodePt(id); if (s) for (const v of s) points.push(v[0], v[1], -v[2]); else if (q && markers.length < 2e5) markers.push({ x: q[0], y: q[1], layer }); }
      continue;
    }
    if (f.prim !== 2 && f.prim !== 3) continue;
    // chain the edges in the order given; a new piece starts where the next edge does not continue the last point
    let cur = null;
    const flush = () => { if (cur && cur.length > 1) { const closed = cur.length > 3 && same(cur[0], cur[cur.length - 1]); if (closed) cur.pop(); nVert += cur.length; polylines.push({ x: cur.map((q) => q[0]), y: cur.map((q) => q[1]), closed, layer, type: f.prim === 3 ? 'area' : 'line', ...(Object.keys(f.attrs).length ? { attrs: f.attrs } : {}) }); } cur = null; };
    for (const [id, ornt] of f.refs) {
      const l = edgePts(id);
      if (!l || !l.length) continue;
      const seq = ornt === 2 ? l.slice().reverse() : l;
      if (cur && same(cur[cur.length - 1], seq[0]) && !(cur.length > 3 && same(cur[0], cur[cur.length - 1]))) for (let k = 1; k < seq.length; k++) cur.push(seq[k]);
      else { flush(); cur = seq.slice(); }
      if (nVert + cur.length > 2e7) fail('The S-57 cell holds too many vertices.');
    }
    flush();
  }
  if (aall > 1) warnings.push('National-language attributes are not read.');
  return { polylines, points, markers, classes, scale, comf, somf, hdat, vdat, counts: { records: nRec, features: feats.length, edges: edges.size, nodes: nodes.size, soundings: points.length / 3 }, warnings };
}
