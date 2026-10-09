// MicroStation / Intergraph design-file reader: DGN V7 (ISFF, Intergraph Standard File Formats).
//
// Read:   2-D and 3-D design files: the TCB (type 9: sub-units per master unit, units of resolution per sub-unit, unit
//         names, global origin), element headers (level, complex and deleted bits, words-to-follow), lines (3), line
//         strings (4), shapes (6), curves (11, evaluated with the ISFF slope-weighted cubic), ellipses (15) and arcs (16)
//         with their 2-D rotation or 3-D quaternion orientation, complex chains (12) and complex shapes (14) joined into
//         one polyline, B-spline poles (21) as their control polygon, point strings (22) and zero-length lines as points.
//         Components of cells (2) are read in place (they carry absolute coordinates). Coordinates are returned in master
//         units with the global origin applied.
// Not read: text (7, 17), tags, dimensions, raster and application elements (counted in `skipped`); the placement of
//         shared-cell instances (35) is not decoded: when a file uses shared cells, each definition (34) is delivered
//         once in its own coordinates; B-spline knots and weights are ignored. Cell libraries are
//         recognised and rejected with an explanation. DGN V8 files (compound files) are read by parseDGN8 below for their
//         basic elements; see its comment for what that covers and how far it was verified.
// File content is untrusted: element lengths and vertex counts are bounds-checked and element / vertex totals are capped.

import { openCFB } from './fmt_cfb.js';
import { inflateSync } from './hdf5.js';

function fail(msg) { const e = new Error(msg); e.user = true; throw e; }
const TAU = 2 * Math.PI, ARC_N = 48;
const TYPE_NAME = { 1: 'cell library header', 2: 'cell header', 5: 'group data', 7: 'text node', 8: 'digitizer setup', 9: 'design file header', 10: 'level symbology', 17: 'text', 18: 'surface', 19: 'solid', 23: 'cone', 24: 'B-spline surface', 25: 'B-spline surface boundary', 26: 'B-spline knots', 27: 'B-spline curve header', 28: 'B-spline weights', 33: 'dimension', 36: 'multiline', 37: 'tag value', 66: 'application data', 87: 'raster header', 88: 'raster data' };

/** 'v7' | 'v8' | 'cell' (V7 cell library) | null. */
export function dgnKind(u8) {
  if (u8.length >= 8 && u8[0] === 0xd0 && u8[1] === 0xcf && u8[2] === 0x11 && u8[3] === 0xe0) return 'v8';
  if (u8.length < 4) return null;
  if (u8[0] === 0x08 && u8[1] === 0x05 && u8[2] === 0x17 && u8[3] === 0x00) return 'cell';
  return (u8[0] === 0x08 || u8[0] === 0xc8) && u8[1] === 0x09 && u8[2] === 0xfe && u8[3] === 0x02 ? 'v7' : null;
}

/**
 * Parse a V7 design file. opts: { maxElements = 2e6, maxVertices = 5e6 }. Returns
 * { version: 'V7', is3d, units: { master, sub, subPerMaster, uorPerSub }, origin: [x, y, z] (master units),
 *   polylines: [{ x, y, z? (3-D files), closed, level, type }], points: flat x, y, z, counts: { type: n },
 *   levels: number[], skipped: { name: n }, warnings: string[] }.
 */
export function parseDGN(u8, opts = {}) {
  const kind = dgnKind(u8);
  if (kind === 'v8') return parseDGN8(u8, opts);
  if (kind === 'cell') fail('This is a MicroStation V7 cell library, not a design file: it holds cell definitions without a model. Place the cells in a design file and save that, or export DXF.');
  if (kind !== 'v7') fail('Not a MicroStation V7 design file (the design-file header element is missing).');
  const maxElements = opts.maxElements ?? 2e6, maxVertices = opts.maxVertices ?? 5e6, N = u8.length;
  const i32 = (p) => (u8[p + 2] | (u8[p + 3] << 8) | (u8[p] << 16) | (u8[p + 1] << 24)), u16 = (p) => u8[p] | (u8[p + 1] << 8);
  /** VAX D-float (word-swapped) → number. */
  const dbl = (p) => {
    const hi = (u8[p] | (u8[p + 1] << 8)) >>> 0, lo = ((u8[p + 2] | (u8[p + 3] << 8)) * 65536 + (u8[p + 4] | (u8[p + 5] << 8))) * 65536 + (u8[p + 6] | (u8[p + 7] << 8)), e = (hi >> 7) & 255;
    if (!e) return 0;
    const v = (1 + ((hi & 127) * 2 ** 48 + lo) / 2 ** 55) * 2 ** (e - 129);
    return hi & 0x8000 ? -v : v;
  };
  const is3d = u8[0] === 0xc8, warnings = [], main = { polylines: [], points: [] }, defs = { polylines: [], points: [] }, counts = {}, skipped = {}, levels = new Set();
  let sink = main;
  let scale = 1, origin = [0, 0, 0], units = { master: '', sub: '', subPerMaster: 0, uorPerSub: 0 }, gotTcb = false, nVert = 0, nElem = 0, cut = false;
  let chain = null, inSharedDef = false;
  const P = (p) => (is3d ? [i32(p) * scale - origin[0], i32(p + 4) * scale - origin[1], i32(p + 8) * scale - origin[2]] : [i32(p) * scale - origin[0], i32(p + 4) * scale - origin[1], 0]);
  const emit = (pts, closed, level, type) => {
    if (pts.length === 2 && pts[0][0] === pts[1][0] && pts[0][1] === pts[1][1] && pts[0][2] === pts[1][2]) pts = [pts[0]];
    if (pts.length === 1) { sink.points.push(pts[0][0], pts[0][1], pts[0][2]); return; }
    if (pts.length < 2) return;
    if (closed && pts.length > 2) { const a = pts[0], b = pts[pts.length - 1]; if (a[0] === b[0] && a[1] === b[1] && a[2] === b[2]) pts.pop(); }
    nVert += pts.length;
    sink.polylines.push({ x: pts.map((q) => q[0]), y: pts.map((q) => q[1]), ...(is3d ? { z: pts.map((q) => q[2]) } : {}), closed, level, type });
  };
  const d2 = (a, b) => (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2;
  /** Append a component to the open complex chain / shape, turning it round when its far end is the nearer one. */
  const link = (pts) => {
    const c = chain.pts;
    if (!c.length) { chain.first = pts.length; c.push(...pts); return; }
    const end = c[c.length - 1];
    if (chain.parts === 1 && chain.first > 1) { const s = c[0]; if (Math.min(d2(s, pts[0]), d2(s, pts[pts.length - 1])) < Math.min(d2(end, pts[0]), d2(end, pts[pts.length - 1]))) c.reverse(); }
    const e2 = c[c.length - 1], seq = d2(e2, pts[pts.length - 1]) < d2(e2, pts[0]) ? pts.slice().reverse() : pts;
    for (let k = d2(e2, seq[0]) <= 1e-18 * (1 + d2(e2, [0, 0, 0])) ? 1 : 0; k < seq.length; k++) c.push(seq[k]);
  };
  const flushChain = () => { if (chain && chain.pts.length) emit(chain.pts, chain.closed, chain.level, chain.type); chain = null; };
  /** ISFF curve: the first and last two vertices only define the end slopes; each span is a cubic in chord length. */
  const curve = (v) => {
    const n = v.length;
    if (n < 6) return v;
    const D = [], M = [], T = [], out = [];
    for (let k = 0; k < n - 1; k++) { const d = Math.hypot(v[k + 1][0] - v[k][0], v[k + 1][1] - v[k][1], v[k + 1][2] - v[k][2]); D.push(d || 1e-4); M.push(d ? [0, 1, 2].map((c) => (v[k + 1][c] - v[k][c]) / d) : [0, 0, 0]); }
    for (let k = 2; k < n - 2; k++) T[k] = [0, 1, 2].map((c) => { const a = Math.abs(M[k + 1][c] - M[k][c]), b = Math.abs(M[k - 1][c] - M[k - 2][c]); return a + b === 0 ? (M[k][c] + M[k - 1][c]) / 2 : (M[k - 1][c] * a + M[k][c] * b) / (a + b); });
    for (let k = 2; k < n - 3; k++) {
      const d = D[k];
      for (let s = 0; s < 8; s++) { const t = (d * s) / 8; out.push([0, 1, 2].map((c) => { const dv = (v[k + 1][c] - v[k][c]) / d, B = (3 * dv - 2 * T[k][c] - T[k + 1][c]) / d, A = (T[k][c] + T[k + 1][c] - 2 * dv) / (d * d); return ((A * t + B) * t + T[k][c]) * t + v[k][c]; })); }
    }
    out.push(v[n - 3]);
    return out;
  };
  /** Ellipse / arc sampled in its own plane, then turned by the 2-D rotation or the 3-D quaternion. */
  const arc = (p, isArc) => {
    let q = p + 36, start = 0, sweep = 360;
    if (isArc) {
      start = i32(q) / 360000;
      const neg = u8[q + 5] & 0x80, raw = ((u8[q + 6] | (u8[q + 7] << 8)) + (u8[q + 4] << 16) + ((u8[q + 5] & 0x7f) * 16777216)) / 360000;
      sweep = (neg ? -raw : raw) || 360;
      q += 8;
    }
    const a = dbl(q) * scale, b = dbl(q + 8) * scale;
    let R, o;
    if (is3d) {
      const w = i32(q + 16) / 2147483648, x = i32(q + 20) / 2147483648, y = i32(q + 24) / 2147483648, z = i32(q + 28) / 2147483648;
      // rows = design-space directions of the element's own x and y axes (checked on wireframes whose arcs must meet their lines)
      R = [[w * w + x * x - y * y - z * z, 2 * (x * y - w * z), 2 * (x * z + w * y)], [2 * (x * y + w * z), w * w - x * x + y * y - z * z, 2 * (y * z - w * x)], [2 * (x * z - w * y), 2 * (y * z + w * x), w * w - x * x - y * y + z * z]];
      o = [dbl(q + 32) * scale - origin[0], dbl(q + 40) * scale - origin[1], dbl(q + 48) * scale - origin[2]];
    } else {
      const r = (i32(q + 16) / 360000) * (Math.PI / 180), c = Math.cos(r), s = Math.sin(r);
      R = [[c, s, 0], [-s, c, 0], [0, 0, 1]];
      o = [dbl(q + 20) * scale - origin[0], dbl(q + 28) * scale - origin[1], 0];
    }
    if (!(a > 0 && b > 0) || ![...o, start, sweep].every(Number.isFinite)) return null;
    const n = Math.max(2, Math.ceil((Math.abs(sweep) / 360) * ARC_N - 1e-9)), pts = [];
    for (let i = 0; i <= n; i++) { const t = ((start + (sweep * i) / n) * Math.PI) / 180, u = a * Math.cos(t), v = b * Math.sin(t); pts.push([o[0] + u * R[0][0] + v * R[1][0], o[1] + u * R[0][1] + v * R[1][1], o[2] + u * R[0][2] + v * R[1][2]]); }
    if (Math.abs(sweep) >= 360 - 1e-9) pts.pop();
    return pts;
  };
  let p = 0;
  for (; p + 4 <= N; nElem++) {
    if (u8[p] === 0xff && u8[p + 1] === 0xff) break;
    const size = 4 + 2 * u16(p + 2), type = u8[p + 1] & 0x7f, level = u8[p] & 0x3f, complex = !!(u8[p] & 0x80), deleted = !!(u8[p + 1] & 0x80);
    if (p + size > N) { warnings.push('The design file is truncated; the elements before the break were read.'); break; }
    if (nElem >= maxElements || nVert > maxVertices) { cut = true; break; }
    const e = p;
    p += size;
    if (chain && (!complex || chain.left <= 0)) flushChain();
    if (!complex) inSharedDef = false;
    if (type === 9 && !gotTcb && size >= 1264) {
      const spm = i32(e + 1112), ups = i32(e + 1116), chars = (o) => String.fromCharCode(...[u8[e + o], u8[e + o + 1]].filter((c) => c >= 32 && c < 127)).trim();
      units = { master: chars(1120), sub: chars(1122), subPerMaster: spm, uorPerSub: ups };
      if (spm > 0 && ups > 0) { scale = 1 / (spm * ups); origin = [dbl(e + 1240) * scale, dbl(e + 1248) * scale, dbl(e + 1256) * scale]; if (!origin.every(Number.isFinite)) origin = [0, 0, 0]; }
      else warnings.push('The design-file header gives no working units; coordinates are returned in units of resolution.');
      gotTcb = true;
      continue;
    }
    if (deleted) { skipped.deleted = (skipped.deleted || 0) + 1; if (chain) chain.left--; continue; }
    if (type === 34) { inSharedDef = true; counts[34] = (counts[34] || 0) + 1; continue; }
    if (type === 35) { counts[35] = (counts[35] || 0) + 1; continue; }
    sink = inSharedDef && complex ? defs : main;             // components of a shared-cell definition are kept apart
    let pts = null, closed = false;
    if (type === 3 && size >= (is3d ? 60 : 52)) pts = [P(e + 36), P(e + (is3d ? 48 : 44))];
    else if ((type === 4 || type === 6 || type === 11 || type === 21 || type === 22) && size >= 38) {
      const n = u16(e + 36), w = is3d ? 12 : 8;
      if (38 + n * w <= size) { pts = []; for (let k = 0; k < n; k++) pts.push(P(e + 38 + k * w)); if (type === 11) pts = curve(pts); closed = type === 6; }
    } else if ((type === 15 && size >= (is3d ? 92 : 72)) || (type === 16 && size >= (is3d ? 100 : 80))) { pts = arc(e, type === 16); closed = type === 15; }
    else if ((type === 12 || type === 14) && size >= 40) {
      flushChain();
      chain = { pts: [], left: u16(e + 38), closed: type === 14, level, type, parts: 0, first: 0 };
      counts[type] = (counts[type] || 0) + 1; levels.add(level);
      continue;
    } else { const nm = TYPE_NAME[type] || `type ${type}`; if (type !== 9 && type !== 10 && type !== 8 && type !== 5 && type !== 2) skipped[nm] = (skipped[nm] || 0) + 1; else if (type === 2) counts[2] = (counts[2] || 0) + 1; if (chain && complex) chain.left--; continue; }
    if (!pts || !pts.length || pts.some((q) => !Number.isFinite(q[0] + q[1] + q[2]))) { skipped.malformed = (skipped.malformed || 0) + 1; if (chain && complex) chain.left--; continue; }
    counts[type] = (counts[type] || 0) + 1; levels.add(level);
    if (chain && complex) { link(pts); chain.parts++; chain.left--; }
    else if (type === 22) for (const q of pts) sink.points.push(q[0], q[1], q[2]);
    else emit(pts, closed, level, type);
  }
  flushChain();
  if (cut) warnings.push('The design file holds more elements than can be shown; the remaining ones were not read.');
  if (!gotTcb) warnings.push('No design-file header (TCB) was found; coordinates are returned in units of resolution.');
  if (counts[21]) warnings.push(`${counts[21]} B-spline curves are shown as their control polygons.`);
  // shared cells: the placement of an instance (type 35) is not decoded, so a definition is shown once, as defined
  const useDefs = defs.polylines.length + defs.points.length > 0 && (counts[35] > 0 || !(main.polylines.length + main.points.length));
  if (useDefs) warnings.push(`${counts[34] || 0} shared-cell definitions are shown once each in their own coordinates; the placements of their ${counts[35] || 0} instances are not applied.`);
  else if (counts[35]) warnings.push(`${counts[35]} shared-cell instances were not expanded.`);
  const polylines = useDefs ? main.polylines.concat(defs.polylines) : main.polylines, points = useDefs ? main.points.concat(defs.points) : main.points;
  return { version: 'V7', is3d, units, origin, polylines, points, counts, levels: [...levels].sort((a, b) => a - b), skipped, warnings };
}

/**
 * MicroStation V8 design file: a compound file whose model storages ("Dgn-Md/#000000" …) hold the graphic elements as one
 * zlib-packed stream ("Dgn^G/$1": 16-byte block header, then the packed element records). No specification is published;
 * the record layout below was established on a file written through the Open Design Alliance library and checked element
 * by element against GDAL's reading of it: a 32-byte header (type, size in 16-bit words behind the first four bytes, offset
 * of the linkages, level, element id, time stamp), 24 bytes of symbology whose property word marks 3-D elements, a 48-byte
 * integer range, then from byte 104 the geometry as IEEE doubles in units of resolution: lines (3), line strings, shapes,
 * curves and point strings (4, 6, 11, 22: a count, then the points), ellipses (15: axes, rotation or quaternion, centre) and
 * arcs (16: start, sweep, axes, rotation or quaternion, centre); complex chains (12), complex shapes (14) and cells (2)
 * are followed by their components. Text, B-splines, shared cells, references, further models and element blocks beyond
 * the first one are not read. Same result structure as for V7, with version 'V8'.
 */
function parseDGN8(u8, opts = {}) {
  const cfb = openCFB(u8), maxElements = opts.maxElements ?? 2e6, maxVertices = opts.maxVertices ?? 5e6, warnings = [], counts = {}, skipped = {}, levels = new Set(), polylines = [], points = [];
  const models = cfb.find(/^Dgn-Md\/#[0-9A-Fa-f]+$/).filter((e) => e.type === 'storage').map((e) => e.path).sort(), how = 'In MicroStation use File > Save As and choose "MicroStation V7 DGN", or File > Export > DXF.';
  if (!models.length) fail(`This compound file holds no MicroStation V8 model (streams: ${cfb.entries.filter((e) => e.type === 'stream').slice(0, 6).map((e) => e.path).join(', ') || 'none'}).`);
  const stream = cfb.find((e) => e.type === 'stream' && e.path.startsWith(models[0] + '/Dgn^G/$'))[0];
  if (!stream) fail(`The default model of this MicroStation V8 design file holds no graphic-element stream. ${how}`);
  const unpack = (b, skip) => { let k = skip; while (k + 2 < b.length && !(b[k] === 0x78 && (b[k + 1] === 0x9c || b[k + 1] === 0x01 || b[k + 1] === 0xda || b[k + 1] === 0x5e))) k++; if (k + 2 >= b.length) fail(`The element stream of this MicroStation V8 design file is not zlib-packed as expected, so it is not read. ${how}`); return inflateSync(b.subarray(k), 256e6); };
  const d = unpack(cfb.read(stream), 16), N = d.length, dv = new DataView(d.buffer, d.byteOffset, N), f64 = (p) => dv.getFloat64(p, true);
  // model header: units of resolution per master unit and the unit labels (UTF-16 linkages)
  let uor = 0, master = '', sub = '', name = '';
  const mh = cfb.find((e) => e.type === 'stream' && e.path === models[0] + '/Dgn~Mh')[0];
  if (mh) {
    try {
      const h = unpack(cfb.read(mh), 0), hv = new DataView(h.buffer, h.byteOffset, h.length), at = 4100;
      if (h.length >= at + 248) { const v = hv.getFloat64(at + 224, true); if (v > 0 && Number.isFinite(v)) uor = v; }
      for (let k = at + 240; k + 16 <= h.length; k += 4) {
        if (hv.getUint16(k + 2, true) !== 0x56d2) continue;
        const id = hv.getUint16(k, true), len = hv.getUint32(k + 8, true);
        if (h[k + 12] !== 0xff || h[k + 13] !== 0xfd || len > 512 || k + 12 + len > h.length) continue;
        let t = '';
        for (let q = k + 14; q + 1 < k + 12 + len; q += 2) { const c = hv.getUint16(q, true); if (!c) break; t += String.fromCharCode(c); }
        if (id === 0x100f) name = t; else if (id === 0x1007) master = t; else if (id === 0x100b) sub = t;
      }
    } catch (e) { if (!(e && e.user) && !(e instanceof RangeError)) throw e; }
  }
  if (!(uor > 0)) { uor = 1; warnings.push('The model header gives no resolution; coordinates are returned in units of resolution.'); }
  let nVert = 0, nEl = 0, is3d = false, group = null, p = 4;
  const count = (t) => { counts[t] = (counts[t] || 0) + 1; };
  const emit = (pts, closed, level, type) => {
    if (group) { for (const q of pts) group.pts.push(q); if (--group.left <= 0) { const g = group; group = null; emit(g.pts.filter((q, i) => !i || q.some((v, k) => v !== g.pts[i - 1][k])), g.type === 14, g.level, g.type); } return; }
    if (pts.length === 1 || (pts.length === 2 && pts[0].every((v, k) => v === pts[1][k]))) { points.push(pts[0][0], pts[0][1], pts[0][2] || 0); return; }
    if (pts.length < 2) return;
    if (closed && pts.length > 2 && pts[0].every((v, k) => v === pts[pts.length - 1][k])) pts = pts.slice(0, -1);
    nVert += pts.length;
    polylines.push({ x: pts.map((q) => q[0]), y: pts.map((q) => q[1]), ...(is3d ? { z: pts.map((q) => q[2] || 0) } : {}), closed, level, type });
  };
  for (; p + 104 <= N; ) {
    const type = d[p], words = dv.getUint32(p + 4, true), attr = dv.getUint32(p + 8, true), level = dv.getUint32(p + 12, true), size = 4 + 2 * words, end = p + Math.min(size, 4 + 2 * Math.max(attr, 52));
    if (!(words >= 50) || p + size > N) { if (p + 4 < N && words) warnings.push('The element stream ends inside an element; the elements read so far are kept.'); break; }
    if (++nEl > maxElements || nVert > maxVertices) { warnings.push('The design file holds more elements than can be shown; the remaining ones were not read.'); break; }
    const three = !!(dv.getUint32(p + 40, true) & 0x800), w = three ? 24 : 16, g = p + 104, pt = (q) => (three ? [f64(q) / uor, f64(q + 8) / uor, f64(q + 16) / uor] : [f64(q) / uor, f64(q + 8) / uor, 0]);
    if (three) is3d = true;
    if (type === 3 && g + 2 * w <= end) { count(3); levels.add(level); emit([pt(g), pt(g + w)], false, level, 3); }
    else if ((type === 4 || type === 6 || type === 11 || type === 22) && g + 8 <= end) {
      const n = dv.getUint32(g, true);
      if (n > 0 && g + 8 + n * w <= end) { count(type); levels.add(level); const pts = []; for (let k = 0; k < n; k++) pts.push(pt(g + 8 + k * w)); if (type === 22 && !group) for (const q of pts) points.push(q[0], q[1], q[2]); else emit(pts, type === 6, level, type); }
      else skipped['damaged element'] = (skipped['damaged element'] || 0) + 1;
    } else if ((type === 15 || type === 16) && g + (type === 16 ? 16 : 0) + (three ? 72 : 40) <= end) {
      const a = type === 16 ? g + 16 : g, start = type === 16 ? f64(g) : 0, sweep = type === 16 ? f64(g + 8) : TAU, r1 = f64(a) / uor, r2 = f64(a + 8) / uor, c = three ? pt(a + 48) : pt(a + 24), n = Math.max(2, Math.ceil((Math.abs(sweep) / TAU) * ARC_N - 1e-9)), pts = [];
      // orientation: a rotation about z (radians) in 2-D, a quaternion (w, x, y, z) in 3-D
      let X = [1, 0, 0], Y = [0, 1, 0];
      if (three) { const qw = f64(a + 16), qx = f64(a + 24), qy = f64(a + 32), qz = f64(a + 40); X = [1 - 2 * (qy * qy + qz * qz), 2 * (qx * qy - qw * qz), 2 * (qx * qz + qw * qy)]; Y = [2 * (qx * qy + qw * qz), 1 - 2 * (qx * qx + qz * qz), 2 * (qy * qz - qw * qx)]; }
      else { const r = f64(a + 16); X = [Math.cos(r), Math.sin(r), 0]; Y = [-Math.sin(r), Math.cos(r), 0]; }
      if ([start, sweep, r1, r2, ...c, ...X, ...Y].every(Number.isFinite)) {
        count(type); levels.add(level);
        const full = type === 15 || Math.abs(Math.abs(sweep) - TAU) < 1e-9;
        for (let k = 0; k <= n - (full && !group ? 1 : 0); k++) { const t = start + (sweep * k) / n, u = r1 * Math.cos(t), v = r2 * Math.sin(t); pts.push([c[0] + u * X[0] + v * Y[0], c[1] + u * X[1] + v * Y[1], c[2] + u * X[2] + v * Y[2]]); }
        emit(pts, full, level, type);
      } else skipped['damaged element'] = (skipped['damaged element'] || 0) + 1;
    } else if ((type === 12 || type === 14) && g + 4 <= end && !group) { const n = dv.getUint32(g, true); count(type); levels.add(level); if (n > 0 && n < 1e6) group = { left: n, pts: [], type, level }; }
    else if (type === 2) count(2);
    else { if (group) { group.left--; if (group.left <= 0) { const q = group; group = null; emit(q.pts, q.type === 14, q.level, q.type); } } const nm = TYPE_NAME[type] || `type ${type}`; skipped[nm] = (skipped[nm] || 0) + 1; }
    p += size;
  }
  if (group && group.pts.length > 1) { const q = group; group = null; emit(q.pts, q.type === 14, q.level, q.type); }
  if (models.length > 1) warnings.push(`The design file holds ${models.length} models; only the first one is read.`);
  warnings.push('MicroStation V8 files have no published specification: the element layout was established on files written through the Open Design Alliance library; check the result against the drawing.');
  return { version: 'V8', is3d, units: { master, sub, subPerMaster: 0, uorPerSub: 0, uorPerMaster: uor }, origin: [0, 0, 0], model: name, polylines, points, counts, levels: [...levels].sort((a, b) => a - b), skipped, warnings };
}
