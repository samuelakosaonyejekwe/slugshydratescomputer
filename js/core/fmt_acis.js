// ACIS save-file reader (text .sat and binary .sab) producing the neutral B-rep model of fmt_brep.js.
//
// Read:   the header (save version, product, date, millimetres per unit, tolerances) of releases 1.x to 3x and the
//         ShapeManager "ASM" flavour written by Autodesk products; text records with implicit or explicit ("-n") indices
//         and the tagged binary records of "ACIS BinaryFile" and "ASM BinaryFile4" streams (checked on ASM 208, 212, 218
//         and 223 as embedded in DWG 2004 to 2018 and on ezdxf-written files), both feeding one model builder; the
//         topology body → lump → shell → face → loop → coedge → edge → vertex → point, the body transform, the curves
//         straight-curve and ellipse-curve (circles and ellipses, sampled between the edge vertices) and intcurve-curve
//         through the B-spline saved with it (exact for "exactcur", the fitted approximation for intersection and
//         other procedural curves), "{ ref n }" subtype references, and the surfaces plane-surface, cone-surface
//         (cylinders, cones, elliptical ones flagged by their ratio), sphere-surface and torus-surface with the face and
//         surface senses.
// Not read: spline-surface faces other than exact B-splines (counted in `skipped`; their boundary is still delivered);
//         the exact form ("exactsur") is decoded like the curves, but no file holding one was at hand, so that path is
//         unverified; procedural curves saved without their B-spline (chord between the vertices), wires, attributes,
//         history, and "ASM BinaryFile8" streams with 64-bit integers (decoded by the same rules, unverified).
// File content is untrusted: record, token and string sizes are capped and every pointer is range-checked.
import { nurbsCurve, nurbsSurfValid, sub, cross, vlen, unit, perp, ellipseArc } from './fmt_brep.js';

function fail(msg) { const e = new Error(msg); e.user = true; throw e; }

/** True for a binary ACIS / ShapeManager file ("ACIS BinaryFile", "ASM BinaryFile4/8"). */
export const isSAB = (u8) => { const h = String.fromCharCode(...u8.subarray(0, 16)); return /^(ACIS|ASM) BinaryFile/.test(h); };

/**
 * Parse a binary SAB file into the neutral B-rep model (see fmt_brep.js). The tagged values are turned into the token
 * records of the text form: logicals become "#t" / "#f", enumeration values "e<n>", pointers "$n", strings "@…".
 */
export function parseSAB(u8) {
  if (!isSAB(u8)) fail('Not a binary ACIS file (the "ACIS BinaryFile" / "ASM BinaryFile" header is missing).');
  const N = u8.length, dv = new DataView(u8.buffer, u8.byteOffset, N), W = u8[14] === 0x38 ? 8 : 4, bad = () => fail('The binary ACIS file is truncated or corrupt.');
  let p = 15;
  const need = (n) => { if (p + n > N) bad(); }, int = () => { need(W); const v = W === 8 ? Number(dv.getBigInt64(p, true)) : dv.getInt32(p, true); p += W; return v; };
  const f64 = () => { need(8); const v = dv.getFloat64(p, true); p += 8; return v; }, str = (n) => { need(n); let t = ''; for (let k = 0; k < n; k++) t += String.fromCharCode(u8[p + k]); p += n; return t; };
  const lenOf = (t) => { let n; if (t === 7 || t === 13 || t === 14) { need(1); n = u8[p++]; } else if (t === 8) { need(2); n = dv.getUint16(p, true); p += 2; } else { need(4); n = dv.getUint32(p, true); p += 4; } if (n > 1e6) bad(); return n; };
  const version = int();
  int(); int(); int();
  if (!(version >= 100 && version < 1e6)) fail('Not a binary ACIS file (the save version is out of range).');
  let unitMm = NaN, product = '';
  if (version >= 200) {
    for (let k = 0; k < 3; k++) { need(1); const t = u8[p++]; if (t < 7 || t > 9) bad(); const v = str(lenOf(t)); if (k === 0) product = v; }
    for (let k = 0; k < 3; k++) { need(1); if (u8[p++] !== 6) bad(); const v = f64(); if (k === 0) unitMm = v; }
  }
  const recs = [];
  let cur = [], name = [], nTok = 0;
  while (p < N) {
    const t = u8[p++];
    if (t === 17) { if (/^End-of-(ACIS|ASM)-data/.test(cur[0])) break; if (cur.length) recs.push(cur); if (recs.length > 2e6) fail('The ACIS file holds too many entities.'); cur = []; name = []; continue; }
    if (++nTok > 5e7) fail('The ACIS file holds too many values.');
    if (t === 2) { need(1); cur.push(dv.getInt8(p++)); }
    else if (t === 3) { need(2); cur.push(dv.getInt16(p, true)); p += 2; }
    else if (t === 4) cur.push(int());
    else if (t === 5) { need(4); cur.push(dv.getFloat32(p, true)); p += 4; }
    else if (t === 6) cur.push(f64());
    else if (t >= 7 && t <= 9) { const n = lenOf(t); need(n); p += n; cur.push('@'); }
    else if (t === 18) { for (const w of str(lenOf(t)).split(/\s+/)) if (w) cur.push(w); }   // literal text, e.g. the numbers of a transform
    else if (t === 10 || t === 11) cur.push(t === 10 ? '#t' : '#f');
    else if (t === 12) cur.push('$' + int());
    else if (t === 13 || t === 14) { const v = str(lenOf(t)); if (cur.length) cur.push(v); else if (t === 14) name.push(v); else cur.push([...name, v].join('-')); }
    else if (t === 15 || t === 16) cur.push(t === 15 ? '{' : '}');
    else if (t === 19 || t === 20) cur.push(f64(), f64(), f64());
    else if (t === 21) cur.push('e' + int());
    else if (t === 22) cur.push(f64(), f64());
    else fail(`The binary ACIS file uses an unknown value tag (${t}); save the model as text ACIS (.sat) or STEP.`);
  }
  if (cur.length && !/^End-of-/.test(cur[0])) recs.push(cur);
  return build(version, product, unitMm, recs, 'ACIS SAB');
}

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
  const recs = [];
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
  return build(version, product, unitMm, recs, 'ACIS SAT');
}

/** Token records (entity name first) -> neutral model. */
function build(version, product, unitMm, recs, format) {
  const warnings = [];
  if (!recs.length) fail('The ACIS file holds no entity records.');
  const B = version >= 700 ? 4 : 2, rec = (i) => (i >= 0 && i < recs.length ? recs[i] : undefined) || null, ptr = (t) => (t && t[0] === '$' ? +t.slice(1) : -1);
  const kindOf = (r) => (r ? r[0] : ''), num3 = (r, o) => [+r[o], +r[o + 1], +r[o + 2]], rev = (t) => t === 'reversed' || t === '1' || t === 'reversed_v' || t === '#t' || t === 'e1';
  const counts = { body: 0, lump: 0, shell: 0, face: 0, loop: 0, coedge: 0, edge: 0, vertex: 0, surfaces: {}, curves: {} }, skipped = {};
  for (const r of recs) if (r && counts[r[0]] !== undefined && typeof counts[r[0]] === 'number') counts[r[0]]++;
  // body transform: p' = scale · (p · M) + t  (row-vector convention of the save file)
  let M = null, T = [0, 0, 0], S = 1;
  for (const r of recs) if (r && r[0] === 'body') { const tr = rec(ptr(r[B + 2])); if (tr && tr[0] === 'transform' && !M) { const o = tr.length < B + 13 ? 2 : String(tr[B - 1])[0] === '$' ? B : B - 1, v = tr.slice(o, o + 13).map(Number); if (v.every(Number.isFinite)) { M = v.slice(0, 9); T = v.slice(9, 12); S = v[12] || 1; } } }
  if (counts.body > 1 && M) warnings.push('The file holds several bodies; the transform of the first one is applied to all.');
  const ident = !M || (M.join() === '1,0,0,0,1,0,0,0,1' && S === 1 && !T.some((v) => v));
  const txv = (v) => (M ? [v[0] * M[0] + v[1] * M[3] + v[2] * M[6], v[0] * M[1] + v[1] * M[4] + v[2] * M[7], v[0] * M[2] + v[1] * M[5] + v[2] * M[8]] : v);
  const tx = (q) => { if (ident) return q; const v = txv(q); return [v[0] * S + T[0], v[1] * S + T[1], v[2] * S + T[2]]; };
  // a vertex names its point last; ASM 212 and later put a count before it
  const pointOf = (vi) => { const v = rec(vi), pt = v && v[0] === 'vertex' ? rec(ptr(String(v[B + 1])[0] === '$' ? v[B + 1] : v[B + 2])) : null; if (!pt || pt[0] !== 'point') return null; const q = num3(pt, B); return q.every(Number.isFinite) ? tx(q) : null; };
  // subtype objects "{ name … }" in file order, for "{ ref n }"
  const subtypes = [];
  recs.forEach((r) => { if (r) for (let i = 0; i + 1 < r.length; i++) if (r[i] === '{' && r[i + 1] !== 'ref' && subtypes.length < 1e6) subtypes.push([r, i]); });
  /**
   * B-spline saved at the start of the first subtype of record r: { name [full] nubs | nurbs, degree(s), closure(s), … }.
   * dim 1: curve -> { name, deg, knots, cps, w }; dim 2: surface -> { name, nurbs } (see nurbsSurfPoint). Knots are saved
   * with their multiplicities and without the first and last repetition.
   */
  const bspline = (r, dim) => {
    let i = r.indexOf('{'), src = r;
    if (i < 0) return null;
    if (src[i + 1] === 'ref') { const t = subtypes[+src[i + 2]]; if (!t) return null; [src, i] = t; }
    const name = src[++i];
    i++;
    if (/^(full|e0)$/.test(src[i])) i++; else if (/^(none|summary|e\d+)$/.test(src[i])) return null;
    const kind = src[i++];
    if (kind !== 'nubs' && kind !== 'nurbs') return null;
    const deg = [], nk = [], knots = [];
    for (let d = 0; d < dim; d++) deg.push(+src[i++]);
    i += dim;                                              // closure in each direction
    if (dim === 2) i += 2;                                 // singularities
    for (let d = 0; d < dim; d++) nk.push(+src[i++]);
    for (let d = 0; d < dim; d++) {
      const k = [];
      if (!(deg[d] >= 1 && deg[d] <= 25 && nk[d] >= 2 && nk[d] <= 1e5)) return null;
      for (let j = 0; j < nk[d]; j++) { const v = +src[i++], m = +src[i++]; if (!(m >= 1 && m <= deg[d] + 1) || k.length > 2e5) return null; for (let q = 0; q < m; q++) k.push(v); }
      knots.push([k[0], ...k, k[k.length - 1]]);
    }
    const n = knots.map((k, d) => k.length - deg[d] - 1), total = n.reduce((a, b) => a * b, 1), w = kind === 'nurbs' ? [] : null, cps = [];
    if (n.some((v) => v <= 0) || total > 4e5 || i + total * (w ? 4 : 3) > src.length) return null;
    for (let j = 0; j < total; j++) { cps.push([+src[i], +src[i + 1], +src[i + 2]]); i += 3; if (w) w.push(+src[i++]); }
    if (dim === 1) return { name, deg: deg[0], knots: knots[0], cps, w };
    // saved with v running fastest: reorder so that u does
    const g = { degU: deg[0], degV: deg[1], nU: n[0], nV: n[1], knotsU: knots[0], knotsV: knots[1], cps: [], w: w ? [] : null };
    for (let b = 0; b < n[1]; b++) for (let a = 0; a < n[0]; a++) { g.cps.push(cps[a * n[1] + b]); if (w) g.w.push(w[a * n[1] + b]); }
    return nurbsSurfValid(g) ? { name, nurbs: g } : null;
  };
  let approx = 0;
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
    if (r[o] !== undefined && String(r[o])[0] !== '$') o++;
    const vb = ptr(r[o++]);
    if (r[o] !== undefined && String(r[o])[0] !== '$') o++;
    const cv = rec(ptr(r[o + 1])), fwd = !rev(r[o + 2]), A = pointOf(va), Bp = pointOf(vb), ck = kindOf(cv);
    if (!A || !Bp) return e;
    const closed = va === vb || vlen(sub(A, Bp)) === 0;
    e.kind = ck === 'straight-curve' ? 'line' : ck === 'ellipse-curve' ? 'ellipse' : ck ? 'spline' : 'other';
    counts.curves[e.kind] = (counts.curves[e.kind] || 0) + 1;
    if (ck === 'ellipse-curve') {
      const c = tx(num3(cv, B)), nrm = unit(txv(num3(cv, B + 3))), Mv = txv(num3(cv, B + 6)).map((q) => q * S), ratio = +cv[B + 9];
      if (c.every(Number.isFinite) && vlen(Mv) > 0 && ratio > 0) { if (Math.abs(ratio - 1) < 1e-12) e.kind = 'circle'; e.pts = ellipseArc(c, Mv, cross(nrm, Mv).map((q) => q * ratio), A, closed ? null : Bp, fwd); e.closed = closed; if (closed) e.pts.pop(); return e; }
    }
    if (ck === 'intcurve-curve') {
      // the B-spline saved with the curve: its exact form (exactcur) or the fitted approximation of a procedural curve
      const bs = bspline(cv, 1);
      let smp = bs ? nurbsCurve(bs.deg, bs.knots, bs.cps, bs.w, NaN, NaN, Math.min(1000, 16 * bs.cps.length)) : null;
      if (smp) {
        smp = smp.map(tx);
        const near = (P) => { let bi = 0, bd = Infinity; smp.forEach((q, i) => { const d = vlen(sub(q, P)); if (d < bd) { bd = d; bi = i; } }); return bi; };
        let ia = near(A), ib = closed ? ia : near(Bp);
        if (bs.name !== 'exactcur') approx++;
        if (closed) { const ring = vlen(sub(smp[0], smp[smp.length - 1])) <= 1e-6 * (vlen(sub(smp[0], smp[smp.length >> 1])) || 1); if (ring) { e.pts = [A, ...smp.slice(ia + 1, -1), ...smp.slice(0, ia)]; e.closed = true; return e; } }
        else if (ia !== ib) { const cut = ia < ib ? smp.slice(ia + 1, ib) : smp.slice(ib + 1, ia).reverse(); e.pts = [A, ...cut, Bp]; return e; }
      }
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
      while (o < r.length && typeof r[o] === 'string' && /^(I|F|T|#t|#f)$/.test(r[o])) o += r[o] === 'I' || r[o] === '#f' ? 1 : 2;      // range of the base ellipse
      const sin = +r[o], cos = +r[o + 1], R = vlen(Mv);
      if (!(R > 0) || !Number.isFinite(sin) || !cos) return { type: 'other' };
      // the cross-section grows along the axis when sine and cosine share their sign; a negative cosine marks a hollow cone
      const s = { type: Math.abs(sin) < 1e-12 ? 'cylinder' : 'cone', o: tx(num3(r, B)), z, x: unit(Mv), r: R, tanA: Math.abs(sin) < 1e-12 ? 0 : sin / cos, hollow: cos < 0 };
      if (Math.abs(ratio - 1) > 1e-9) s.ratio = ratio;
      return s;
    }
    if (k === 'sphere-surface') { const rad = +r[B + 3] * S, x = unit(txv(num3(r, B + 4))), z = unit(txv(num3(r, B + 7))); return { type: 'sphere', o: tx(num3(r, B)), z: vlen(z) ? z : [0, 0, 1], x: vlen(x) ? x : [1, 0, 0], r: Math.abs(rad), hollow: rad < 0 }; }
    if (k === 'torus-surface') { const z = unit(txv(num3(r, B + 3))), R = +r[B + 6] * S, rr = +r[B + 7] * S; return R > 0 && Math.abs(rr) < R ? { type: 'torus', o: tx(num3(r, B)), z, x: perp(z, txv(num3(r, B + 8))), R, r: Math.abs(rr), hollow: rr < 0 } : { type: 'other' }; }
    if (k === 'spline-surface') {
      // only an exact B-spline surface is delivered; the fitted approximations of procedural surfaces are not
      const bs = bspline(r, 2);
      if (bs && bs.name === 'exactsur') { const g = bs.nurbs; g.cps = g.cps.map(tx); return { type: 'spline', nurbs: g }; }
      return { type: 'spline' };
    }
    return { type: 'other' };
  };
  const faces = [];
  let guardTotal = 0;
  recs.forEach((r, fi) => {
    if (!r || r[0] !== 'face') return;
    const surf = surfOf(ptr(r[B + 4])) || { type: 'other' }, loops = [];
    counts.surfaces[surf.type] = (counts.surfaces[surf.type] || 0) + 1;
    if ((surf.type === 'spline' && !surf.nurbs) || surf.type === 'other') { const nm = kindOf(rec(ptr(r[B + 4]))) || 'no surface'; skipped[nm] = (skipped[nm] || 0) + 1; }
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
  if (approx) warnings.push(`${approx} edges lie on procedural curves (intersections, offsets); they follow the B-spline approximation saved with them.`);
  if (!ident && M && M.join() !== '1,0,0,0,1,0,0,0,1') warnings.push('The body carries a rotation; it was applied to every point and direction.');
  return { format, version: version >= 20000 ? `ASM ${Math.floor(version / 100)}` : `${Math.floor(version / 100)}.${version % 100}`, saveVersion: version, product: product.trim(), unitScale: unitMm > 0 ? unitMm / 1000 : null, faces, edges, vertices: counts.vertex, counts, skipped, warnings };
}
