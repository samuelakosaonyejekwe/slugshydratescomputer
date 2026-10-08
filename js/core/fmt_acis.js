// ACIS save-file reader (text .sat) producing the neutral B-rep model of fmt_brep.js.
//
// Read:   the header (save version, product, date, millimetres per unit, tolerances) of releases 1.x to 3x (and the
//         ShapeManager "ASM" flavour written by Autodesk products), the entity records with implicit or explicit ("-n")
//         indices, the topology body → lump → shell → face → loop → coedge → edge → vertex → point, the body transform,
//         the curves straight-curve and ellipse-curve (circles and ellipses, sampled between the edge vertices), and the
//         surfaces plane-surface, cone-surface (cylinders, cones, elliptical ones flagged by their ratio), sphere-surface
//         and torus-surface with the face and surface senses.
// Not read: spline-surface faces (counted in `skipped`; their boundary is still delivered), intcurve-curve and other
//         procedural curves (edges on them become the chord between their vertices), wires, attributes, history, and
//         binary .sab files, whose tagged encoding is recognised and rejected with an explanation.
// File content is untrusted: record and token counts are capped and every pointer is range-checked.
import { sub, cross, vlen, unit, perp, ellipseArc } from './fmt_brep.js';

function fail(msg) { const e = new Error(msg); e.user = true; throw e; }

/** True for a binary ACIS / ShapeManager file ("ACIS BinaryFile", "ASM BinaryFile4/8"). */
export const isSAB = (u8) => { const h = String.fromCharCode(...u8.subarray(0, 16)); return /^(ACIS|ASM) BinaryFile/.test(h); };
export function parseSAB() { return fail('This is a binary ACIS file (.sab): its tagged binary records are not read. Save the model as text ACIS (.sat) or as STEP AP214 / AP242 from the CAD package.'); }

/** Parse a text SAT file into the neutral B-rep model (see fmt_brep.js). */
export function parseSAT(text) {
  const n = text.length;
  let p = 0;
  const ws = () => { while (p < n && text.charCodeAt(p) <= 32) p++; };
  const tok = () => { ws(); const s = p; while (p < n && text.charCodeAt(p) > 32) p++; return text.slice(s, p); };
  const version = +tok();
  tok(); tok(); tok();
  if (!(version >= 100 && version < 1e6) || n < 20) fail('Not an ACIS SAT text file (the version header is missing).');
  let unitMm = NaN, product = '';
  if (version >= 200) {
    for (let k = 0; k < 3; k++) { const t = tok(), len = +(t[0] === '@' ? t.slice(1) : t); if (!(len >= 0 && len < 4096)) fail('The ACIS header is malformed.'); p++; if (k === 0) product = text.substr(p, len); p += len; }
    unitMm = +tok(); tok(); tok();
  }
  // records: whitespace-separated tokens up to '#'; "@n" introduces a counted string
  const recs = [], warnings = [];
  let cur = [], seq = 0;
  while (p < n) {
    const c = text.charCodeAt(p);
    if (c <= 32) { p++; continue; }
    if (c === 35) {
      p++;
      if (cur.length) { let idx = seq; if (/^-\d+$/.test(cur[0])) { idx = +cur[0].slice(1); cur.shift(); } if (idx > 4e6 || recs.length > 2e6) fail('The ACIS file holds too many entities.'); recs[idx] = cur; seq = idx + 1; cur = []; }
      continue;
    }
    if (c === 64) { let q = p + 1; while (q < n && text.charCodeAt(q) > 32) q++; const len = +text.slice(p + 1, q); if (len >= 0 && q + 1 + len <= n) { cur.push('@'); p = q + 1 + len; continue; } }
    const s = p;
    while (p < n && text.charCodeAt(p) > 32 && text.charCodeAt(p) !== 35) p++;
    const t = text.slice(s, p);
    if (/^End-of-(ACIS|ASM)-data/.test(t)) break;
    if (cur.length < 1e5) cur.push(t);
  }
  if (!recs.length) fail('The ACIS file holds no entity records.');
  const B = version >= 700 ? 4 : 2, rec = (i) => (i >= 0 && i < recs.length ? recs[i] : undefined) || null, ptr = (t) => (t && t[0] === '$' ? +t.slice(1) : -1);
  const kindOf = (r) => (r ? r[0] : ''), num3 = (r, o) => [+r[o], +r[o + 1], +r[o + 2]], rev = (t) => t === 'reversed' || t === '1' || t === 'reversed_v';
  const counts = { body: 0, lump: 0, shell: 0, face: 0, loop: 0, coedge: 0, edge: 0, vertex: 0, surfaces: {}, curves: {} }, skipped = {};
  for (const r of recs) if (r && counts[r[0]] !== undefined && typeof counts[r[0]] === 'number') counts[r[0]]++;
  // body transform: p' = scale · (p · M) + t  (row-vector convention of the save file)
  let M = null, T = [0, 0, 0], S = 1;
  for (const r of recs) if (r && r[0] === 'body') { const tr = rec(ptr(r[B + 2])); if (tr && tr[0] === 'transform' && !M) { const o = tr.length >= B + 13 ? B : 2, v = tr.slice(o, o + 13).map(Number); if (v.every(Number.isFinite)) { M = v.slice(0, 9); T = v.slice(9, 12); S = v[12] || 1; } } }
  if (counts.body > 1 && M) warnings.push('The file holds several bodies; the transform of the first one is applied to all.');
  const ident = !M || (M.join() === '1,0,0,0,1,0,0,0,1' && S === 1 && !T.some((v) => v));
  const txv = (v) => (M ? [v[0] * M[0] + v[1] * M[3] + v[2] * M[6], v[0] * M[1] + v[1] * M[4] + v[2] * M[7], v[0] * M[2] + v[1] * M[5] + v[2] * M[8]] : v);
  const tx = (q) => { if (ident) return q; const v = txv(q); return [v[0] * S + T[0], v[1] * S + T[1], v[2] * S + T[2]]; };
  const pointOf = (vi) => { const v = rec(vi), pt = v && v[0] === 'vertex' ? rec(ptr(v[B + 1])) : null; if (!pt || pt[0] !== 'point') return null; const q = num3(pt, B); return q.every(Number.isFinite) ? tx(q) : null; };
  // edges, sampled once in their own direction
  const edgeCache = new Map(), edges = [];
  const edgePts = (ei) => {
    let e = edgeCache.get(ei);
    if (e) return e;
    const r = rec(ei);
    e = { pts: [], kind: 'other', closed: false };
    edgeCache.set(ei, e);
    if (!r || r[0] !== 'edge') return e;
    // start vertex [start parameter] end vertex [end parameter] coedge curve sense
    let o = B;
    const va = ptr(r[o++]);
    if (r[o] && r[o][0] !== '$') o++;
    const vb = ptr(r[o++]);
    if (r[o] && r[o][0] !== '$') o++;
    const cv = rec(ptr(r[o + 1])), fwd = !rev(r[o + 2]), A = pointOf(va), Bp = pointOf(vb), ck = kindOf(cv);
    if (!A || !Bp) return e;
    const closed = va === vb || vlen(sub(A, Bp)) === 0;
    e.kind = ck === 'straight-curve' ? 'line' : ck === 'ellipse-curve' ? 'ellipse' : ck ? 'spline' : 'other';
    counts.curves[e.kind] = (counts.curves[e.kind] || 0) + 1;
    if (ck === 'ellipse-curve') {
      const c = tx(num3(cv, B)), nrm = unit(txv(num3(cv, B + 3))), Mv = txv(num3(cv, B + 6)).map((q) => q * S), ratio = +cv[B + 9];
      if (c.every(Number.isFinite) && vlen(Mv) > 0 && ratio > 0) { if (Math.abs(ratio - 1) < 1e-12) e.kind = 'circle'; e.pts = ellipseArc(c, Mv, cross(nrm, Mv).map((q) => q * ratio), A, closed ? null : Bp, fwd); e.closed = closed; if (closed) e.pts.pop(); return e; }
    }
    if (ck && ck !== 'straight-curve' && ck !== 'ellipse-curve') skipped[ck] = (skipped[ck] || 0) + 1;
    e.pts = closed ? [A] : [A, Bp]; e.closed = false; e.null = !ck;
    return e;
  };
  const surfOf = (si) => {
    const r = rec(si), k = kindOf(r);
    if (!r) return null;
    if (k === 'plane-surface') { const z = unit(txv(num3(r, B + 3))); return { type: 'plane', o: tx(num3(r, B)), z, x: perp(z, txv(num3(r, B + 6))), hollow: false }; }
    if (k === 'cone-surface') {
      const z = unit(txv(num3(r, B + 3))), Mv = txv(num3(r, B + 6)).map((q) => q * S), ratio = +r[B + 9];
      let o = B + 10;
      while (o < r.length && (r[o] === 'I' || r[o] === 'F' || r[o] === 'T')) o += r[o] === 'I' ? 1 : 2;      // range of the base ellipse
      const sin = +r[o], cos = +r[o + 1], R = vlen(Mv);
      if (!(R > 0) || !Number.isFinite(sin) || !cos) return { type: 'other' };
      // the cross-section grows along the axis when sine and cosine share their sign; a negative cosine marks a hollow cone
      const s = { type: Math.abs(sin) < 1e-12 ? 'cylinder' : 'cone', o: tx(num3(r, B)), z, x: unit(Mv), r: R, tanA: Math.abs(sin) < 1e-12 ? 0 : sin / cos, hollow: cos < 0 };
      if (Math.abs(ratio - 1) > 1e-9) s.ratio = ratio;
      return s;
    }
    if (k === 'sphere-surface') { const rad = +r[B + 3] * S, x = unit(txv(num3(r, B + 4))), z = unit(txv(num3(r, B + 7))); return { type: 'sphere', o: tx(num3(r, B)), z: vlen(z) ? z : [0, 0, 1], x: vlen(x) ? x : [1, 0, 0], r: Math.abs(rad), hollow: rad < 0 }; }
    if (k === 'torus-surface') { const z = unit(txv(num3(r, B + 3))), R = +r[B + 6] * S, rr = +r[B + 7] * S; return R > 0 && Math.abs(rr) < R ? { type: 'torus', o: tx(num3(r, B)), z, x: perp(z, txv(num3(r, B + 8))), R, r: Math.abs(rr), hollow: rr < 0 } : { type: 'other' }; }
    return { type: k === 'spline-surface' ? 'spline' : 'other' };
  };
  const faces = [];
  let guardTotal = 0;
  recs.forEach((r, fi) => {
    if (!r || r[0] !== 'face') return;
    const surf = surfOf(ptr(r[B + 4])) || { type: 'other' }, loops = [];
    counts.surfaces[surf.type] = (counts.surfaces[surf.type] || 0) + 1;
    if (surf.type === 'spline' || surf.type === 'other') { const nm = kindOf(rec(ptr(r[B + 4]))) || 'no surface'; skipped[nm] = (skipped[nm] || 0) + 1; }
    for (let li = ptr(r[B + 1]), nl = 0; li >= 0 && nl < 1e5; nl++) {
      const lp = rec(li);
      if (!lp || lp[0] !== 'loop') break;
      const first = ptr(lp[B + 1]), pts = [];
      let single = false;
      for (let ci = first, nc = 0; ci >= 0 && nc < 1e5; nc++) {
        const ce = rec(ci);
        if (!ce || ce[0] !== 'coedge' || ++guardTotal > 2e7) break;
        const e = edgePts(ptr(ce[B + 3])), seq = rev(ce[B + 4]) && !e.closed ? e.pts.slice().reverse() : rev(ce[B + 4]) ? [e.pts[0], ...e.pts.slice(1).reverse()] : e.pts;
        if (e.null && seq.length === 1) single = true;
        for (let k = 0; k < seq.length - (e.closed ? 0 : 1) || (k === 0 && seq.length === 1); k++) pts.push(seq[k]);
        ci = ptr(ce[B]);
        if (ci === first) break;
      }
      if (pts.length) loops.push({ pts, single: single && pts.length === 1 });
      li = ptr(lp[B]);
    }
    const { hollow, ...s } = surf;
    faces.push({ surf: s, sense: !rev(r[B + 5]) !== !!hollow, loops, id: fi });
  });
  recs.forEach((r, i) => { if (r && r[0] === 'edge') { const e = edgePts(i); if (e.pts.length > 1) edges.push({ pts: e.pts, kind: e.kind, closed: e.closed }); } });
  if (!faces.length && !edges.length) fail('The ACIS file holds no faces or edges (it may contain attributes or history only).');
  if (!ident && M && M.join() !== '1,0,0,0,1,0,0,0,1') warnings.push('The body carries a rotation; it was applied to every point and direction.');
  return { format: 'ACIS SAT', version: version >= 20000 ? `ASM ${Math.floor(version / 100)}` : `${Math.floor(version / 100)}.${version % 100}`, saveVersion: version, product: product.trim(), unitScale: unitMm > 0 ? unitMm / 1000 : null, faces, edges, vertices: counts.vertex, counts, skipped, warnings };
}
