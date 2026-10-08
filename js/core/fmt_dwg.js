// AutoCAD DWG reader, written against the Open Design Alliance "Open Design Specification for .dwg files".
//
// Containers: R13 / R14 / 2000 (AC1012, AC1014, AC1015: section locator, plain sections) and 2004 / 2010 / 2013 / 2018
//         (AC1018, AC1024, AC1027, AC1032: XOR-masked file header, section page map, section map, LZ-compressed data
//         pages). The 2007 container (AC1021: Reed-Solomon coded pages, its own compression) and releases before R13
//         are recognised and rejected with the way out.
// Read:   the object map (handle → offset) and the bit-coded objects: LINE, POINT, CIRCLE, ARC, ELLIPSE, LWPOLYLINE
//         (bulges as arcs, elevation, closed flag), POLYLINE_2D / _3D with their VERTEX chains, SPLINE (returned as its
//         definition for the caller to evaluate), 3DFACE, INSERT and MINSERT of block definitions (scale, rotation,
//         extrusion, rows and columns; nested, depth-capped), the LAYER table (entity → layer name), the object
//         coordinate system of planar entities and, for 2000 drawings, the drawing unit $INSUNITS from the header
//         variables. Only model space is delivered; invisible entities are left out.
// Not read: text, dimensions, hatches, leaders, solids, regions and ACIS bodies, images, attributes, polyface and
//         polygon meshes (counted in `skipped` where their type is known), line types, colours, paper space, the other
//         header variables (so $INSUNITS of 2004 and later drawings is not known); checksums are not verified.
// File content is untrusted: every offset is bounds-checked, decompressed sizes, object, vertex and nesting counts are
// capped, and chains of linked entities are walked with a visit guard.

function fail(msg) { const e = new Error(msg); e.user = true; throw e; }
const TAU = 2 * Math.PI, ARC_N = 48;
const RELEASE = { 'MC0.0': 'R1.1', 'AC1.2': 'R1.2', 'AC1.4': 'R1.4', 'AC1.50': 'R2.0', 'AC2.10': 'R2.10', AC1002: 'R2.5', AC1003: 'R2.6', AC1004: 'R9', AC1006: 'R10', AC1009: 'R11 / R12', AC1012: 'R13', AC1014: 'R14', AC1015: '2000', AC1018: '2004', AC1021: '2007', AC1024: '2010', AC1027: '2013', AC1032: '2018' };
const TYPE = { 1: 'TEXT', 2: 'ATTRIB', 3: 'ATTDEF', 4: 'BLOCK', 5: 'ENDBLK', 6: 'SEQEND', 7: 'INSERT', 8: 'MINSERT', 10: 'VERTEX_2D', 11: 'VERTEX_3D', 12: 'VERTEX_MESH', 13: 'VERTEX_PFACE', 14: 'VERTEX_PFACE_FACE', 15: 'POLYLINE_2D', 16: 'POLYLINE_3D', 17: 'ARC', 18: 'CIRCLE', 19: 'LINE', 20: 'DIMENSION', 21: 'DIMENSION', 22: 'DIMENSION', 23: 'DIMENSION', 24: 'DIMENSION', 25: 'DIMENSION', 26: 'DIMENSION', 27: 'POINT', 28: '3DFACE', 29: 'POLYLINE_PFACE', 30: 'POLYLINE_MESH', 31: 'SOLID', 32: 'TRACE', 33: 'SHAPE', 34: 'VIEWPORT', 35: 'ELLIPSE', 36: 'SPLINE', 37: 'REGION', 38: '3DSOLID', 39: 'BODY', 40: 'RAY', 41: 'XLINE', 43: 'OLEFRAME', 44: 'MTEXT', 45: 'LEADER', 46: 'TOLERANCE', 47: 'MLINE', 49: 'BLOCK_HEADER', 51: 'LAYER', 74: 'OLE2FRAME', 77: 'LWPOLYLINE', 78: 'HATCH' };
const isEntityType = (t) => (t >= 1 && t <= 41 && t !== 9) || (t >= 43 && t <= 47) || t === 74 || t === 77 || t === 78;

/** { code: 'AC1015', release: '2000', readable } for a DWG file, or null. */
export function dwgVersion(u8) {
  const code = String.fromCharCode(...u8.subarray(0, 6));
  if (!/^(AC\d\.\d\d?|AC\d{4}|MC0\.0)/.test(code)) return null;
  const c = /^(AC\d{4}|MC0\.0)/.test(code) ? code.slice(0, 6).replace(/\0.*$/, '') : code.replace(/[^\x20-\x7e].*$/, '');
  return { code: c, release: RELEASE[c] || 'unknown', readable: ['AC1012', 'AC1014', 'AC1015', 'AC1018', 'AC1024', 'AC1027', 'AC1032'].includes(c) };
}

/** Bit-stream reader over the file (most significant bit first; multi-byte raw values little-endian). */
class Bits {
  constructor(u8, bit, end) { this.u8 = u8; this.p = bit; this.end = end; this.f8 = new DataView(new ArrayBuffer(8)); }
  b() { if (this.p >= this.end) throw new RangeError('eof'); const v = (this.u8[this.p >> 3] >> (7 - (this.p & 7))) & 1; this.p++; return v; }
  bb() { return (this.b() << 1) | this.b(); }
  rc() { if (this.p + 8 > this.end) throw new RangeError('eof'); const i = this.p >> 3, s = this.p & 7; this.p += 8; return s ? ((this.u8[i] << s) | (this.u8[i + 1] >> (8 - s))) & 255 : this.u8[i]; }
  rs() { const a = this.rc(); return a | (this.rc() << 8); }
  rl() { const a = this.rs(); return (a | (this.rs() << 16)) >>> 0; }
  rd() { for (let k = 0; k < 8; k++) this.f8.setUint8(k, this.rc()); return this.f8.getFloat64(0, true); }
  bs() { const c = this.bb(); return c === 0 ? this.rs() : c === 1 ? this.rc() : c === 2 ? 0 : 256; }
  bl() { const c = this.bb(); return c === 0 ? this.rl() : c === 1 ? this.rc() : 0; }
  bd() { const c = this.bb(); return c === 0 ? this.rd() : c === 1 ? 1 : 0; }
  /** double with a default: 0, 4, 6 or 8 bytes patched into the default value. */
  dd(def) {
    const c = this.bb();
    if (c === 0) return def;
    if (c === 3) return this.rd();
    this.f8.setFloat64(0, def, true);
    if (c === 2) { this.f8.setUint8(4, this.rc()); this.f8.setUint8(5, this.rc()); }
    for (let k = 0; k < 4; k++) this.f8.setUint8(k, this.rc());
    return this.f8.getFloat64(0, true);
  }
  h() { const code = (this.bb() << 2) | this.bb(), n = (this.bb() << 2) | this.bb(); let v = 0; for (let k = 0; k < n; k++) v = v * 256 + this.rc(); return { code, v }; }
  t() { const n = this.bs(); if (n > 4096) throw new RangeError('text'); let s = ''; for (let k = 0; k < n; k++) { const c = this.rc(); if (c) s += String.fromCharCode(c); } return s; }
  p3() { return [this.bd(), this.bd(), this.bd()]; }
  /** UTF-16 text of the 2007+ string stream. */
  tu() { const n = this.bs(); if (n > 4096) throw new RangeError('text'); let s = ''; for (let k = 0; k < n; k++) { const c = this.rs(); if (c) s += String.fromCharCode(c); } return s; }
  mc() { let v = 0, sh = 0; for (let k = 0; k < 5; k++) { const c = this.rc(); v += (c & 127) * 2 ** sh; sh += 7; if (!(c & 128)) break; } return v; }
}

/** DWG 2004 section decompression (an LZ77 variant with literal runs and back references). */
function unLZ(src, size) {
  const out = new Uint8Array(size), n = src.length;
  let p = 0, o = 0, op = 0;
  const rd = () => { if (p >= n) throw new RangeError('lz'); return src[p++]; };
  const lit = () => { const b = rd(); if (b & 0xf0) { op = b; return 0; } if (b) return b + 3; let t = 0x0f, c = rd(); while (!c) { t += 0xff; c = rd(); } return t + c + 3; };
  const long = () => { let t = 0, c = rd(); while (!c) { t += 0xff; c = rd(); } return t + c; };
  const copyLit = (k) => { if (p + k > n || o + k > size) throw new RangeError('lz'); out.set(src.subarray(p, p + k), o); p += k; o += k; };
  copyLit(lit());
  for (let guard = 0; guard < 4e8; guard++) {
    if (!op) op = rd();
    const c = op;
    op = 0;
    if (c === 0x11) break;
    let len, off, nl;
    if (c >= 0x40) { len = (c >> 4) - 1; off = ((c >> 2) & 3) | (rd() << 2); nl = c & 3; }
    else {
      if (c < 0x10) throw new RangeError('lz');
      len = c === 0x10 ? long() + 9 : c < 0x20 ? (c & 15) + 2 : c === 0x20 ? long() + 0x21 : c - 0x1e;
      const b0 = rd();
      off = (b0 >> 2) | (rd() << 6); nl = b0 & 3;
      if (c < 0x20) off += 0x3fff;
    }
    const from = o - off - 1;
    if (from < 0 || o + len > size) throw new RangeError('lz');
    for (let k = 0; k < len; k++) out[o + k] = out[from + k];
    o += len;
    copyLit(nl || lit());
  }
  return out;
}
/** Logical sections of a 2004-layout file: name → bytes, for the names asked for. */
function sections2004(u8, want) {
  const N = u8.length, dv = new DataView(u8.buffer, u8.byteOffset, N), bad = (w) => fail(`The DWG file is truncated or corrupt (${w}).`);
  if (N < 0x100) bad('file header');
  const h = u8.slice(0x80, 0x80 + 0x6c), hd = new DataView(h.buffer);
  for (let k = 0, seed = 1; k < 0x6c; k++) { seed = (Math.imul(seed, 0x343fd) + 0x269ec3) | 0; h[k] ^= (seed >>> 16) & 255; }
  if (String.fromCharCode(...h.subarray(0, 11)) !== 'AcFssFcAJMB') bad('encrypted file header');
  const system = (addr, type) => { if (!(addr >= 0x100 && addr + 20 <= N) || dv.getUint32(addr, true) !== type) bad('system section'); const ds = dv.getUint32(addr + 4, true), cs = dv.getUint32(addr + 8, true), ct = dv.getUint32(addr + 12, true); if (ds > 64e6 || addr + 20 + cs > N) bad('system section'); return ct === 2 ? unLZ(u8.subarray(addr + 20, addr + 20 + cs), ds) : u8.subarray(addr + 20, addr + 20 + ds); };
  const pm = system(hd.getUint32(0x54, true) + hd.getUint32(0x58, true) * 4294967296 + 0x100, 0x41630e3b), pmv = new DataView(pm.buffer, pm.byteOffset, pm.byteLength), pages = new Map();
  for (let p = 0, addr = 0x100; p + 8 <= pm.length; ) { const id = pmv.getInt32(p, true), size = pmv.getUint32(p + 4, true); p += 8; if (id >= 0) pages.set(id, addr); else p += 16; addr += size; }
  const sm = system(pages.get(hd.getUint32(0x5c, true)) ?? -1, 0x4163003b), smv = new DataView(sm.buffer, sm.byteOffset, sm.byteLength), nd = smv.getUint32(0, true), out = {};
  for (let p = 20, k = 0; k < nd && k < 1000 && p + 96 <= sm.length; k++) {
    const size = smv.getUint32(p, true) + smv.getUint32(p + 4, true) * 4294967296, np = smv.getUint32(p + 8, true), maxPage = smv.getUint32(p + 12, true), comp = smv.getUint32(p + 20, true), enc = smv.getUint32(p + 28, true), name = String.fromCharCode(...sm.subarray(p + 32, p + 96)).replace(/\0[\s\S]*$/, '');
    p += 96;
    if (np > 1e6 || p + 16 * np > sm.length) bad('section map');
    if (want.includes(name)) {
      if (size > 512e6 || enc === 1) fail('A section of the DWG file is too large or encrypted.');
      const buf = new Uint8Array(size);
      for (let q = 0; q < np; q++) {
        const e = p + 16 * q, pa = pages.get(smv.getUint32(e, true)), cs = smv.getUint32(e + 4, true), start = smv.getUint32(e + 8, true) + smv.getUint32(e + 12, true) * 4294967296;
        if (pa === undefined || pa + 32 + cs > N || start > size) bad('data page');
        const ds = Math.min(maxPage, size - start), raw = u8.subarray(pa + 32, pa + 32 + cs);   // a page always unpacks to the full page size
        if (maxPage > 16e6) bad('data page');
        buf.set((comp === 2 ? unLZ(raw, maxPage) : raw).subarray(0, ds), start);
      }
      out[name] = buf;
    }
    p += 16 * np;
  }
  return out;
}
const ref = (h, base) => (h.code === 6 ? base + 1 : h.code === 8 ? base - 1 : h.code === 10 ? base + h.v : h.code === 12 ? base - h.v : h.v);
/** Axes of the object coordinate system of extrusion direction n (the "arbitrary axis algorithm"). */
function ocs(n) {
  const l = Math.hypot(n[0], n[1], n[2]) || 1, z = [n[0] / l, n[1] / l, n[2] / l];
  if (Math.abs(z[0]) < 1e-12 && Math.abs(z[1]) < 1e-12 && z[2] > 0) return null;
  const a = Math.abs(z[0]) < 1 / 64 && Math.abs(z[1]) < 1 / 64 ? [z[2], 0, -z[0]] : [-z[1], z[0], 0], al = Math.hypot(a[0], a[1], a[2]), x = [a[0] / al, a[1] / al, a[2] / al];
  return { x, y: [z[1] * x[2] - z[2] * x[1], z[2] * x[0] - z[0] * x[2], z[0] * x[1] - z[1] * x[0]], z };
}
const toWcs = (o, p) => (o ? [p[0] * o.x[0] + p[1] * o.y[0] + p[2] * o.z[0], p[0] * o.x[1] + p[1] * o.y[1] + p[2] * o.z[1], p[0] * o.x[2] + p[1] * o.y[2] + p[2] * o.z[2]] : p);
/** Arc of a polyline bulge between two points (same z): intermediate points only. */
function bulgePts(a, b, bulge) {
  const ang = 4 * Math.atan(bulge), n = Math.min(4 * ARC_N, Math.ceil((Math.abs(ang) / TAU) * ARC_N - 1e-9)), out = [];
  if (!(n > 1) || !Number.isFinite(ang)) return out;
  const dx = b[0] - a[0], dy = b[1] - a[1], k = (1 / bulge - bulge) / 2, cx = (a[0] + b[0]) / 2 - (k * dy) / 2, cy = (a[1] + b[1]) / 2 + (k * dx) / 2, r = Math.hypot(a[0] - cx, a[1] - cy), t0 = Math.atan2(a[1] - cy, a[0] - cx);
  for (let i = 1; i < n; i++) { const t = t0 + (ang * i) / n; out.push([cx + r * Math.cos(t), cy + r * Math.sin(t), a[2]]); }
  return out;
}

/**
 * Parse a DWG file. opts: { maxEntities = 500000, maxVertices = 5e6, maxInsertDepth = 8 }. Returns
 * { version: 'AC1015', release: '2000', insunits: $INSUNITS code (0 none, 1 in, 2 ft, 4 mm, 5 cm, 6 m …; 2000 drawings only) or null,
 *   polylines: [{ x, y, z? (only when some z ≠ 0), closed, layer, type }]   world coordinates, curves flattened (48 segments
 *     per full turn); a SPLINE comes as { spline: { degree, knots, points, weights, fit, closed }, layer, type: 'SPLINE' }
 *     for the caller to evaluate,
 *   points: flat x, y, z, faces: flat triangles (3DFACE), counts: { TYPE: n } (entities emitted, block contents included),
 *   layers: string[], skipped: { TYPE: n }, warnings: string[] }
 */
export function parseDWG(u8, opts = {}) {
  const ver = dwgVersion(u8);
  if (!ver) fail('Not a DWG file (the "AC10xx" version signature is missing).');
  const dxf = 'Save the drawing as ASCII DXF, which is read: in AutoCAD use SAVEAS and choose DXF, or convert it with the free ODA File Converter (output "ASCII DXF").';
  if (!ver.readable) fail(ver.code === 'AC1021' ? `This is an AutoCAD 2007 drawing (AC1021): its pages are Reed-Solomon coded and compressed in a layout of their own, which is not read (R13, R14, 2000, 2004, 2010, 2013 and 2018 drawings are). ${dxf}` : `This is an AutoCAD ${ver.release} drawing (${ver.code}), a layout older than R13 that is not read. ${dxf}`);
  const V = +ver.code.slice(4), r2000 = V >= 15, r2004 = V >= 18, r2007 = V >= 21, r2010 = V >= 24, r2013 = V >= 27;
  let file = u8, mapBytes = null;
  if (r2004) {
    let secs;
    try { secs = sections2004(u8, ['AcDb:Handles', 'AcDb:AcDbObjects']); } catch (e) { if (e instanceof RangeError) fail('The DWG file is truncated or corrupt (compressed section).'); throw e; }
    file = secs['AcDb:AcDbObjects']; mapBytes = secs['AcDb:Handles'];
    if (!file || !mapBytes) fail('The DWG file is truncated or corrupt (object sections are missing).');
    u8 = file;
  }
  const N = u8.length, dv = new DataView(u8.buffer, u8.byteOffset, N), maxEntities = opts.maxEntities ?? 5e5, maxVertices = opts.maxVertices ?? 5e6, maxDepth = opts.maxInsertDepth ?? 8;
  const bad = (what) => fail(`The DWG file is truncated or corrupt (${what}).`);
  const sec = {};
  if (!r2004) {
    if (N < 0x40) bad('header');
    const nsec = dv.getUint32(0x15, true);
    if (nsec > 16 || 0x19 + 9 * nsec > N) bad('section locator');
    for (let k = 0; k < nsec; k++) { const o = 0x19 + 9 * k; sec[u8[o]] = [dv.getUint32(o + 1, true), dv.getUint32(o + 5, true)]; }
    if (!sec[2] || sec[2][0] + sec[2][1] > N) bad('object map');
    mapBytes = u8.subarray(sec[2][0], sec[2][0] + sec[2][1]);
  }
  // classes: number → DXF name / entity flag (variable object types start at 500)
  const classes = [];
  if (!r2004 && sec[1] && sec[1][0] + 20 < N && sec[1][1] > 20) {
    try {
      const size = dv.getUint32(sec[1][0] + 16, true), b = new Bits(u8, (sec[1][0] + 20) * 8, Math.min(N, sec[1][0] + 20 + size) * 8);
      while (b.p + 40 < b.end && classes.length < 4096) { const num = b.bs(); b.bs(); b.t(); b.t(); const name = b.t(); b.b(); const id = b.bs(); classes[num - 500] = { name, entity: id === 0x1f2 }; }
    } catch (e) { if (!(e instanceof RangeError)) throw e; }
  }
  // R2000 header variables up to $INSUNITS: d BD, s BS, l BL, b B, t text, h handle, 3 3BD, 2 2RD (the list of the ODA specification)
  let insunits = null;
  if (V === 15 && sec[0] && sec[0][0] + 20 < N) {
    try {
      const hb = new Bits(u8, (sec[0][0] + 20) * 8, Math.min(N, sec[0][0] + 20 + dv.getUint32(sec[0][0] + 16, true)) * 8), first = hb.bd();
      const seq = 'dddttttllhbbbbbbbbbbbbbbbbbbbbsssssssssssssssssssssssssssdddddddddddddddddddddtllllllllshhhhhhd33322d333hhsh33333333322d333hhsh333333ttdddddddddbbbbbbsssdddddddddbsbbbbssssssssssssssbbssssbshhhhhsshhhhhhhhhhhhhsstthhhl';
      for (const c of seq) { if (c === 'd') hb.bd(); else if (c === 's') hb.bs(); else if (c === 'l') hb.bl(); else if (c === 'b') hb.b(); else if (c === 't') hb.t(); else if (c === 'h') hb.h(); else if (c === '3') hb.p3(); else { hb.rd(); hb.rd(); } }
      const v = hb.bs();
      if (first === 412148564080 && v >= 0 && v <= 24) insunits = v;
    } catch (e) { if (!(e instanceof RangeError)) throw e; }
  }
  // object map: sections of (handle delta, offset delta) pairs
  const map = new Map();
  for (let p = 0, end = mapBytes.length, guard = 0; p + 2 <= end && guard < 1e5; guard++) {
    const size = (mapBytes[p] << 8) | mapBytes[p + 1], se = p + size;
    if (size <= 2 || se > end) break;
    let q = p + 2, handle = 0, loc = 0;
    while (q < se) {
      let v = 0, sh = 0, c;
      do { c = mapBytes[q++]; v += (c & 127) * 2 ** sh; sh += 7; } while (c & 128 && q < se);
      handle += v;
      let o = 0, neg = false;
      sh = 0;
      for (;;) { c = mapBytes[q++]; if (c & 128 && q < se) { o += (c & 127) * 2 ** sh; sh += 7; } else { o += (c & 63) * 2 ** sh; neg = !!(c & 64); break; } }
      loc += neg ? -o : o;
      if (loc > 0 && loc < N) map.set(handle, loc);
      if (map.size > 4e6) fail('The DWG file holds too many objects.');
    }
    p = se + 2;
  }
  if (!map.size) bad('object map is empty');
  const warnings = [], skipped = {}, counts = {}, layerName = new Map(), blocks = new Map(), ents = new Map(), byOwner = new Map(), model = [];
  let nDecoded = 0, nBad = 0;
  /** Decode the object at a file offset: common data, the supported entity bodies, then the handle stream. */
  const decode = (off) => {
    let size = 0, sh = 0, q = off;
    for (let k = 0; k < 2; k++) { if (q + 2 > N) return null; const w = dv.getUint16(q, true); q += 2; size |= (w & 0x7fff) << sh; sh += 15; if (!(w & 0x8000)) break; }
    if (!size || q + size > N) return null;
    let start = q * 8, bitsize = 0;
    const b = new Bits(u8, start, (q + size) * 8);
    let type;
    if (r2010) { const hbits = b.mc(); start = b.p; b.end = start + size * 8; if (b.end > N * 8) return null; bitsize = size * 8 - hbits; const c = b.bb(); type = c === 0 ? b.rc() : c === 1 ? b.rc() + 0x1f0 : b.rs(); }   // the size counts from behind the handle-stream size
    else { type = b.bs(); if (r2000) bitsize = b.rl(); }
    const handle = b.h().v, o = { type, handle };
    // 2007+: texts live in a string stream at the end of the object data, located through a flag and a size word
    let str = null;
    if (r2007 && bitsize > 17 && bitsize <= size * 8) {
      let sp = start + bitsize - 1;
      if ((u8[sp >> 3] >> (7 - (sp & 7))) & 1) { const w = new Bits(u8, sp - 16, sp); let ss = w.rs(); sp -= 16; if (ss & 0x8000) { const w2 = new Bits(u8, sp - 16, sp); ss = (ss & 0x7fff) | (w2.rs() << 15); sp -= 16; } if (ss > 0 && sp - ss >= start) str = new Bits(u8, sp - ss, sp); }
    }
    const T = () => (r2007 ? (str ? str.tu() : '') : b.t());
    for (let es = b.bs(), g = 0; es && g < 1000; g++, es = b.bs()) { b.h(); b.p += 8 * es; }
    const cls = type >= 500 ? classes[type - 500] : null, name = cls ? cls.name : TYPE[type], entity = cls ? cls.entity : isEntityType(type);
    o.name = name || `TYPE_${type}`;
    let mode = 0, nreact = 0, nolinks = true, ltByLayer = true, ltf = 0, psf = 0, matf = 0, vsf = 0, xdicMissing = false, colorRef = false, owned = 0;
    if (entity) {
      if (b.b()) { const gs = b.rl(); if (gs > size) throw new RangeError('graphics'); b.p += 8 * gs; }
      if (!r2000) bitsize = b.rl();
      mode = b.bb(); nreact = b.bl();
      if (!r2000) ltByLayer = !!b.b();
      if (r2004) { xdicMissing = !!b.b(); if (r2013) b.b(); const cf = b.bs(); if (cf & 0x2000) b.bl(); if (cf & 0x8000) b.bl(); colorRef = !!(cf & 0x4000); }
      else { nolinks = !!b.b(); b.bs(); }
      b.bd();
      if (r2000) { ltf = b.bb(); psf = b.bb(); }
      if (r2007) { matf = b.bb(); b.rc(); }
      if (r2010) { vsf = b.b() + b.b() + b.b(); }
      o.invisible = b.bs() & 1;
      if (r2000) b.rc();
      o.entity = true; o.mode = mode;
    } else { if (!r2000) bitsize = b.rl(); nreact = b.bl(); if (r2004) { xdicMissing = !!b.b(); if (r2013) b.b(); } }
    if (nreact > 1e5) throw new RangeError('reactors');
    const p2 = () => [b.rd(), b.rd()], be = () => (r2000 && b.b() ? [0, 0, 1] : b.p3()), bt = () => (r2000 && b.b() ? 0 : b.bd());
    const nm = o.name;
    if (nm === 'LINE') {
      if (r2000) { const zf = b.b(), x1 = b.rd(), x2 = b.dd(x1), y1 = b.rd(), y2 = b.dd(y1), z1 = zf ? 0 : b.rd(), z2 = zf ? 0 : b.dd(z1); o.a = [x1, y1, z1]; o.b = [x2, y2, z2]; }
      else { o.a = b.p3(); o.b = b.p3(); }
    } else if (nm === 'POINT') o.pt = b.p3();
    else if (nm === 'CIRCLE' || nm === 'ARC') { o.c = b.p3(); o.r = b.bd(); bt(); o.n = be(); if (nm === 'ARC') { o.a0 = b.bd(); o.a1 = b.bd(); } }
    else if (nm === 'ELLIPSE') { o.c = b.p3(); o.major = b.p3(); o.n = b.p3(); o.ratio = b.bd(); o.a0 = b.bd(); o.a1 = b.bd(); }
    else if (nm === 'LWPOLYLINE') {
      const fl = b.bs();
      if (fl & 4) b.bd();
      o.elev = fl & 8 ? b.bd() : 0;
      if (fl & 2) b.bd();
      o.n = fl & 1 ? b.p3() : [0, 0, 1];
      const np = b.bl(), nb = fl & 16 ? b.bl() : 0, nid = r2010 && fl & 1024 ? b.bl() : 0, nw = fl & 32 ? b.bl() : 0;
      if (np > maxVertices || nb > np || nid > np) throw new RangeError('lwpolyline');
      o.closed = !!(fl & 512); o.pts = []; o.bulges = [];
      for (let k = 0; k < np; k++) { const pr = o.pts[k - 1]; o.pts.push(r2000 && k ? [b.dd(pr[0]), b.dd(pr[1])] : p2()); }
      for (let k = 0; k < nb; k++) o.bulges.push(b.bd());
      for (let k = 0; k < nid; k++) b.bl();
      for (let k = 0; k < nw; k++) { b.bd(); b.bd(); }
    } else if (nm === 'POLYLINE_2D') { o.closed = !!(b.bs() & 1); b.bs(); b.bd(); b.bd(); bt(); o.elev = b.bd(); o.n = be(); o.poly = 2; if (r2004) owned = b.bl(); }
    else if (nm === 'POLYLINE_3D') { b.rc(); o.closed = !!(b.rc() & 1); o.poly = 3; if (r2004) owned = b.bl(); }
    else if (nm === 'VERTEX_2D') { o.vflags = b.rc(); o.pt = b.p3(); const sw = b.bd(); if (sw >= 0) b.bd(); o.bulge = b.bd(); }
    else if (nm === 'VERTEX_3D') { o.vflags = b.rc(); o.pt = b.p3(); }
    else if (nm === '3DFACE') {
      if (r2000) { const noFlags = b.b(), z0 = b.b(), c1 = [b.rd(), b.rd(), z0 ? 0 : b.rd()], c2 = [b.dd(c1[0]), b.dd(c1[1]), b.dd(c1[2])], c3 = [b.dd(c2[0]), b.dd(c2[1]), b.dd(c2[2])], c4 = [b.dd(c3[0]), b.dd(c3[1]), b.dd(c3[2])]; o.corners = [c1, c2, c3, c4]; if (!noFlags) b.bs(); }
      else o.corners = [b.p3(), b.p3(), b.p3(), b.p3()];
    } else if (nm === 'SPLINE') {
      let scen = b.bl();
      if (r2013) { const f1 = b.bl(), kp = b.bl(); if (f1 & 1) scen = 2; if (kp === 15) scen = 1; }
      const s = { degree: b.bl(), knots: [], points: [], weights: null, fit: [], closed: false };
      if (scen === 2) { b.bd(); b.p3(); b.p3(); const nf = b.bl(); if (nf > 1e5) throw new RangeError('spline'); for (let k = 0; k < nf; k++) s.fit.push(b.p3()); }
      else if (scen === 1) {
        b.b(); s.closed = !!b.b(); b.b(); b.bd(); b.bd();
        const nk = b.bl(), nc = b.bl(), wt = b.b();
        if (nk > 1e5 || nc > 1e5) throw new RangeError('spline');
        for (let k = 0; k < nk; k++) s.knots.push(b.bd());
        if (wt) s.weights = [];
        for (let k = 0; k < nc; k++) { s.points.push(b.p3()); if (wt) s.weights.push(b.bd()); }
      }
      o.spline = s;
    } else if (nm === 'INSERT' || nm === 'MINSERT') {
      o.ins = b.p3();
      if (r2000) { const f = b.bb(); if (f === 3) o.scale = [1, 1, 1]; else if (f === 1) { const y = b.dd(1); o.scale = [1, y, b.dd(1)]; } else if (f === 2) { const x = b.rd(); o.scale = [x, x, x]; } else { const x = b.rd(), y = b.dd(x); o.scale = [x, y, b.dd(x)]; } }
      else o.scale = b.p3();
      o.rot = b.bd(); o.n = b.p3(); o.attribs = !!b.b();
      if (r2004 && o.attribs) owned = b.bl();
      if (nm === 'MINSERT') { o.cols = b.bs(); o.rows = b.bs(); o.dx = b.bd(); o.dy = b.bd(); }
    } else if (nm === 'BLOCK_HEADER') {
      o.bname = T();
      if (r2007) b.bs(); else { b.b(); b.bs(); b.b(); }
      b.b(); b.b(); o.xref = !!b.b(); o.overlay = !!b.b();
      if (r2000) b.b();
      if (r2004) owned = b.bl();
      o.base = b.p3(); T();
      if (r2000) { o.ninsert = 0; for (let g = 0; g < 1e6 && b.rc(); g++) o.ninsert++; T(); const ps = b.bl(); if (ps > size) throw new RangeError('preview'); }
    } else if (nm === 'LAYER') o.lname = T();
    if (owned > maxVertices) throw new RangeError('owned objects');
    // the handle stream follows the data, at the bit position given by the object size
    if (!(bitsize > 0 && bitsize <= size * 8)) return o;
    const hs = new Bits(u8, start + bitsize, b.end), H = () => ref(hs.h(), handle);
    try {
      if (entity) {
        if (mode === 0) o.owner = H();
        for (let k = 0; k < nreact; k++) hs.h();
        if (!xdicMissing) hs.h();
        if (!r2000) { o.layer = H(); if (!ltByLayer) hs.h(); }
        if (r2004) { if (colorRef) hs.h(); } else if (!nolinks) { o.prev = H(); o.next = H(); } else o.next = handle + 1;
        if (r2000) { o.layer = H(); if (ltf === 3) hs.h(); if (matf === 3) hs.h(); if (psf === 3) hs.h(); for (let k = 0; k < vsf; k++) hs.h(); }
        if (o.poly && r2004) { o.owned = []; for (let k = 0; k < owned; k++) o.owned.push(H()); }
        else if (o.poly) { o.first = H(); o.last = H(); }
        else if (nm === 'INSERT' || nm === 'MINSERT') o.block = H();
      } else if (nm === 'BLOCK_HEADER') {
        hs.h();
        for (let k = 0; k < nreact; k++) hs.h();
        hs.h(); hs.h(); o.blockEnt = H();
        if (!o.xref && !o.overlay) { o.first = H(); o.last = H(); }
      }
    } catch (e) { if (!(e instanceof RangeError)) throw e; }
    return o;
  };
  for (const [handle, off] of map) {
    let o = null;
    try { o = decode(off); } catch (e) { if (!(e instanceof RangeError)) throw e; nBad++; }
    if (!o) continue;
    if (++nDecoded > 4e6) break;
    if (o.name === 'LAYER') layerName.set(o.handle || handle, o.lname);
    else if (o.name === 'BLOCK_HEADER') blocks.set(o.handle || handle, o);
    else if (o.entity) {
      ents.set(o.handle, o);
      if (o.mode === 2) model.push(o);
      else if (o.mode === 0 && o.owner !== undefined) { const l = byOwner.get(o.owner); if (l) l.push(o); else byOwner.set(o.owner, [o]); }
    }
  }
  // geometry
  const polylines = [], points = [], faces = [];
  let nVert = 0, nEnt = 0, cut = false;
  const ident = (p) => p;
  const emit = (pts, closed, o, T, type) => {
    if (pts.length < 2) return;
    const w = T === ident ? pts : pts.map(T);
    if (closed && w.length > 2) { const a = w[0], z = w[w.length - 1]; if (a[0] === z[0] && a[1] === z[1] && a[2] === z[2]) w.pop(); }
    nVert += w.length;
    const has3 = w.some((q) => q[2] !== 0);
    polylines.push({ x: w.map((q) => q[0]), y: w.map((q) => q[1]), ...(has3 ? { z: w.map((q) => q[2]) } : {}), closed, layer: layerName.get(o.layer) || '', type });
  };
  const arc = (c, r, n, a0, sweep) => { const ax = ocs(n), k = Math.max(2, Math.ceil((Math.abs(sweep) / TAU) * ARC_N - 1e-9)), out = []; for (let i = 0; i <= k; i++) { const t = a0 + (sweep * i) / k; out.push(toWcs(ax, [c[0] + r * Math.cos(t), c[1] + r * Math.sin(t), c[2]])); } return out; };
  /** Vertices of an old-style polyline: the linked chain first → … → last. */
  const vertices = (pl) => {
    const out = [];
    const ok = (v) => v && v.pt && !(v.vflags & 16) && (v.name === 'VERTEX_2D' || v.name === 'VERTEX_3D');
    if (pl.owned) { for (const h of pl.owned) { const v = ents.get(h); if (ok(v)) out.push(v); } return out; }
    for (let h = pl.first, g = 0; h !== undefined && g < maxVertices; g++) { const v = ents.get(h); if (!v) break; if (ok(v)) out.push(v); if (h === pl.last) break; h = v.next; }
    return out;
  };
  const count = (nm) => { counts[nm] = (counts[nm] || 0) + 1; };
  const draw = (o, T, depth, seen) => {
    if (cut) return;
    if (++nEnt > maxEntities || nVert > maxVertices) { cut = true; return; }
    const nm = o.name;
    if (nm === 'LINE') { count(nm); if (o.a[0] === o.b[0] && o.a[1] === o.b[1] && o.a[2] === o.b[2]) { const q = T(o.a); points.push(q[0], q[1], q[2]); } else emit([o.a, o.b], false, o, T, nm); }
    else if (nm === 'POINT') { count(nm); const q = T(o.pt); points.push(q[0], q[1], q[2]); }
    else if (nm === 'CIRCLE') { count(nm); if (o.r > 0) { const pts = arc(o.c, o.r, o.n, 0, TAU); pts.pop(); emit(pts, true, o, T, nm); } }
    else if (nm === 'ARC') { count(nm); let sw = (o.a1 - o.a0) % TAU; if (sw <= 0) sw += TAU; if (o.r > 0) emit(arc(o.c, o.r, o.n, o.a0, sw), false, o, T, nm); }
    else if (nm === 'ELLIPSE') {
      count(nm);
      const M = o.major, m = [(o.n[1] * M[2] - o.n[2] * M[1]) * o.ratio, (o.n[2] * M[0] - o.n[0] * M[2]) * o.ratio, (o.n[0] * M[1] - o.n[1] * M[0]) * o.ratio];
      let sw = o.a1 - o.a0;
      if (sw <= 0 || sw > TAU) sw = ((sw % TAU) + TAU) % TAU || TAU;
      const k = Math.max(2, Math.ceil((sw / TAU) * ARC_N - 1e-9)), pts = [], full = sw >= TAU - 1e-9;
      for (let i = 0; i <= k; i++) { const t = o.a0 + (sw * i) / k, c = Math.cos(t), s = Math.sin(t); pts.push([o.c[0] + c * M[0] + s * m[0], o.c[1] + c * M[1] + s * m[1], o.c[2] + c * M[2] + s * m[2]]); }
      if (full) pts.pop();
      emit(pts, full, o, T, nm);
    } else if (nm === 'LWPOLYLINE' || o.poly === 2) {
      count(nm === 'LWPOLYLINE' ? nm : 'POLYLINE');
      const ax = ocs(o.n), src = nm === 'LWPOLYLINE' ? o.pts.map((q, k) => ({ p: [q[0], q[1], o.elev], bulge: o.bulges[k] || 0 })) : vertices(o).map((v) => ({ p: [v.pt[0], v.pt[1], o.elev], bulge: v.bulge || 0 })), pts = [];
      src.forEach((v, k) => { pts.push(v.p); const nx = src[k + 1] || (o.closed ? src[0] : null); if (nx && v.bulge) pts.push(...bulgePts(v.p, nx.p, v.bulge)); });
      emit(ax ? pts.map((q) => toWcs(ax, q)) : pts, o.closed, o, T, nm === 'LWPOLYLINE' ? nm : 'POLYLINE');
    } else if (o.poly === 3) { count('POLYLINE'); emit(vertices(o).map((v) => v.pt), o.closed, o, T, 'POLYLINE'); }
    else if (nm === '3DFACE') { count(nm); const c = o.corners.map(T); faces.push(...c[0], ...c[1], ...c[2]); if (c[3].some((v, k) => v !== c[2][k])) faces.push(...c[0], ...c[2], ...c[3]); }
    else if (nm === 'SPLINE') { count(nm); const s = o.spline; polylines.push({ spline: { ...s, points: s.points.map(T), fit: s.fit.map(T) }, layer: layerName.get(o.layer) || '', type: nm }); nVert += s.points.length + s.fit.length; }
    else if (nm === 'INSERT' || nm === 'MINSERT') {
      const blk = blocks.get(o.block);
      if (!blk || depth >= maxDepth || seen.has(o.block)) { skipped[blk ? 'INSERT (nested too deeply)' : 'INSERT (block not found)'] = (skipped[blk ? 'INSERT (nested too deeply)' : 'INSERT (block not found)'] || 0) + 1; return; }
      count(nm);
      const list = byOwner.get(o.block) || [], ax = ocs(o.n), c = Math.cos(o.rot), s = Math.sin(o.rot), base = blk.base || [0, 0, 0], inner = new Set(seen).add(o.block);
      for (let r = 0; r < Math.min(o.rows || 1, 1000); r++) for (let k = 0; k < Math.min(o.cols || 1, 1000); k++) {
        const ox = k * (o.dx || 0), oy = r * (o.dy || 0);
        const T2 = (q) => { const u = (q[0] - base[0]) * o.scale[0], v = (q[1] - base[1]) * o.scale[1], w = (q[2] - base[2]) * o.scale[2]; return T(toWcs(ax, [o.ins[0] + (u + ox) * c - (v + oy) * s, o.ins[1] + (u + ox) * s + (v + oy) * c, o.ins[2] + w])); };
        for (const e of list) if (!e.invisible) draw(e, T2, depth + 1, inner);
      }
    } else if (!['VERTEX_2D', 'VERTEX_3D', 'VERTEX_MESH', 'VERTEX_PFACE', 'VERTEX_PFACE_FACE', 'SEQEND', 'BLOCK', 'ENDBLK', 'ATTRIB', 'VIEWPORT'].includes(nm)) skipped[nm] = (skipped[nm] || 0) + 1;
  };
  for (const o of model) if (!o.invisible) draw(o, ident, 0, new Set());
  if (cut) warnings.push('The drawing holds more entities than can be shown; the remaining ones were not read.');
  if (nBad) warnings.push(`${nBad} objects could not be decoded and were skipped.`);
  if (!model.length) warnings.push('No model-space entities were found in the drawing.');
  return { version: ver.code, release: ver.release, insunits, polylines, points, faces, counts, layers: [...new Set(layerName.values())].filter(Boolean), skipped, warnings };
}
