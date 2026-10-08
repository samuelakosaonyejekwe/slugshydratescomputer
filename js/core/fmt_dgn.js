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
//         once in its own coordinates; B-spline knots and weights are ignored. DGN V8 files (OLE2 compound documents with an unpublished element stream) and
//         cell libraries are recognised and rejected with an explanation.
// File content is untrusted: element lengths and vertex counts are bounds-checked and element / vertex totals are capped.

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
  if (kind === 'v8') fail('This is a MicroStation V8 design file: an OLE2 compound document whose element stream is not publicly documented, so it is not read (V7 design files are). In MicroStation use File > Save As and choose "MicroStation V7 DGN", or File > Export > DXF.');
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
