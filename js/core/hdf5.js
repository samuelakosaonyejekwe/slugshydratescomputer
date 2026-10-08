// Read-only HDF5 container reader in plain JavaScript (browser and Node), written against the HDF5 File Format
// Specification version 3.0. It serves the geometry importer (CGNS, MED, Exodus II, NetCDF-4, MATLAB v7.3, plain .h5).
//
// Read:   superblock versions 0-3 (with a user block), version 1 and 2 object headers with continuation blocks,
//         old-style groups (symbol table: v1 B-tree + local heap), new-style groups (link messages; dense link storage
//         through the fractal heap and its v2 B-tree name index), hard and soft links, committed datatypes,
//         datatypes fixed-point, floating-point (16/32/64 bit), fixed and variable-length strings, enum, array,
//         object reference, opaque / bitfield (raw) and compounds of those, simple / scalar / null dataspaces,
//         compact, contiguous and chunked layouts (v1 B-tree; version 4 layouts with single-chunk, implicit, fixed-array,
//         extensible-array and v2 B-tree chunk indexes), filters deflate, shuffle and fletcher32 (checksum stripped, not
//         verified), fill values, attributes in the object header and dense attribute storage, hyperslab selection.
// Not read: external links and external raw-data files, virtual datasets, szip / n-bit / scale-offset / third-party
//         filters, filtered fractal heaps, "huge" fractal-heap objects, region references, time datatypes, shared-message
//         heaps (SOHM), free-space and driver blocks; metadata checksums are not verified; nothing is ever written.
// File content is untrusted: every address and length is bounds-checked, object, node and recursion counts are capped.

function fail(msg) { const e = new Error(msg); e.user = true; throw e; }
const latin1 = new TextDecoder('latin1'), utf8 = new TextDecoder();
const HOST_LE = new Uint8Array(new Uint16Array([1]).buffer)[0] === 1;
const FILTER_NAME = { 1: 'deflate', 2: 'shuffle', 3: 'fletcher32', 4: 'szip', 5: 'n-bit', 6: 'scale-offset', 307: 'bzip2', 32000: 'LZF', 32001: 'Blosc', 32004: 'LZ4', 32008: 'Bitshuffle', 32013: 'zfp', 32015: 'Zstandard' };

/** Offset of the HDF5 signature (0, or 512 · 2^k behind a user block), or -1. */
export function hdf5Offset(u8) {
  for (let o = 0; o + 8 <= u8.length && o <= (1 << 26); o = o ? o * 2 : 512) if (u8[o] === 0x89 && u8[o + 1] === 0x48 && u8[o + 2] === 0x44 && u8[o + 3] === 0x46 && u8[o + 4] === 0x0d && u8[o + 5] === 0x0a && u8[o + 6] === 0x1a && u8[o + 7] === 0x0a) return o;
  return -1;
}

/**
 * Synchronous inflate of a zlib (RFC 1950) or raw deflate (RFC 1951) stream into at most `max` bytes. Small buffers (HDF5
 * chunks, record blobs) are decoded far faster this way than through one DecompressionStream each.
 */
export function inflateSync(src, max = 512e6, raw = false) {
  const n = src.length, bad = () => fail('Compressed data in the file is corrupt or truncated.');
  let p = raw ? 0 : 2, bb = 0, bn = 0, out = new Uint8Array(Math.min(max, Math.max(1024, n * 4))), o = 0;
  if (!raw && (n < 2 || (src[0] & 15) !== 8 || ((src[0] << 8) | src[1]) % 31 || src[1] & 32)) bad();
  const bits = (k) => { while (bn < k) { if (p >= n) bad(); bb |= src[p++] << bn; bn += 8; } const v = bb & ((1 << k) - 1); bb >>>= k; bn -= k; return v; };
  const grow = (k) => { if (o + k > max) fail('Decompressed data is too large.'); if (o + k > out.length) { const t = new Uint8Array(Math.min(max, Math.max(out.length * 2, o + k))); t.set(out.subarray(0, o)); out = t; } };
  const build = (len, off, cnt) => {                       // canonical Huffman table: [count per length, symbols sorted by code]
    const c = new Uint16Array(16), s = new Uint16Array(cnt), offs = new Uint16Array(16);
    for (let i = 0; i < cnt; i++) c[len[off + i]]++;
    c[0] = 0;
    for (let i = 1; i < 16; i++) offs[i] = offs[i - 1] + c[i - 1];
    for (let i = 0; i < cnt; i++) if (len[off + i]) s[offs[len[off + i]]++] = i;
    return [c, s];
  };
  const sym = (t) => { let code = 0, first = 0, idx = 0; for (let l = 1; l < 16; l++) { code |= bits(1); const c = t[0][l]; if (code - c < first) return t[1][idx + code - first]; idx += c; first = (first + c) << 1; code <<= 1; } return bad(); };
  const LB = [3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 23, 27, 31, 35, 43, 51, 59, 67, 83, 99, 115, 131, 163, 195, 227, 258], LE = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0];
  const DB = [1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 257, 385, 513, 769, 1025, 1537, 2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577], DE = [0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13];
  const ORD = [16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15];
  let fixed = null;
  for (let last = 0; !last;) {
    last = bits(1);
    const type = bits(2);
    if (type === 0) {
      bb = 0; bn = 0;
      if (p + 4 > n) bad();
      const len = src[p] | (src[p + 1] << 8);
      p += 4;
      if (p + len > n) bad();
      grow(len); out.set(src.subarray(p, p + len), o); o += len; p += len;
      continue;
    }
    if (type === 3) bad();
    let lt, dt;
    if (type === 1) {
      if (!fixed) { const l = new Uint8Array(320); l.fill(8, 0, 144); l.fill(9, 144, 256); l.fill(7, 256, 280); l.fill(8, 280, 288); l.fill(5, 288, 320); fixed = [build(l, 0, 288), build(l, 288, 30)]; }
      [lt, dt] = fixed;
    } else {
      const nl = bits(5) + 257, nd = bits(5) + 1, nc = bits(4) + 4, cl = new Uint8Array(19), l = new Uint8Array(nl + nd);
      if (nl > 286 || nd > 30) bad();
      for (let i = 0; i < nc; i++) cl[ORD[i]] = bits(3);
      const ct = build(cl, 0, 19);
      for (let i = 0; i < nl + nd;) {
        const s = sym(ct);
        if (s < 16) { l[i++] = s; continue; }
        const prev = s === 16 ? (i ? l[i - 1] : bad()) : 0, rep = s === 16 ? 3 + bits(2) : s === 17 ? 3 + bits(3) : 11 + bits(7);
        if (i + rep > nl + nd) bad();
        l.fill(prev, i, i + rep); i += rep;
      }
      lt = build(l, 0, nl); dt = build(l, nl, nd);
    }
    for (;;) {
      const s = sym(lt);
      if (s < 256) { if (o >= out.length) grow(1); out[o++] = s; continue; }
      if (s === 256) break;
      if (s > 285) bad();
      const len = LB[s - 257] + bits(LE[s - 257]), ds = sym(dt);
      if (ds > 29) bad();
      const dist = DB[ds] + bits(DE[ds]);
      if (dist > o) bad();
      grow(len);
      if (dist >= len) out.copyWithin(o, o - dist, o - dist + len); else for (let i = 0; i < len; i++) out[o + i] = out[o - dist + i];
      o += len;
    }
  }
  return out.subarray(0, o);
}

async function inflateZlib(bytes, max) {
  if (bytes.length < 1 << 17) return inflateSync(bytes, max);
  if (typeof DecompressionStream !== 'function') fail('This browser cannot decompress data (DecompressionStream is missing).');
  const chunks = [];
  let n = 0;
  try {
    const reader = new Blob([bytes]).stream().pipeThrough(new DecompressionStream('deflate')).getReader();
    for (;;) { const { done, value } = await reader.read(); if (done) break; n += value.length; if (n > max) { await reader.cancel().catch(() => {}); fail('Decompressed data is too large.'); } chunks.push(value); }
  } catch (e) { if (e && e.user) throw e; fail('Compressed data in the HDF5 file is corrupt or truncated.'); }
  const out = new Uint8Array(n);
  let p = 0;
  for (const c of chunks) { out.set(c, p); p += c.length; }
  return out;
}

/** Readable name of a datatype descriptor ('float64', 'int32', 'string', 'compound{…}' …). */
export function typeName(t) {
  if (!t) return 'unknown';
  switch (t.cls) {
    case 0: return (t.signed ? 'int' : 'uint') + t.size * 8;
    case 1: return 'float' + t.size * 8;
    case 3: return 'string';
    case 4: return 'bitfield' + t.size * 8;
    case 5: return 'opaque';
    case 6: return `compound{${t.members.map((m) => m.name).slice(0, 6).join(',')}}`;
    case 7: return 'reference';
    case 8: return 'enum';
    case 9: return t.vlen === 1 ? 'string' : `vlen<${typeName(t.base)}>`;
    case 10: return `${typeName(t.base)}[${t.dims.join('×')}]`;
    default: return 'class' + t.cls;
  }
}
/** True for datatypes that decode to plain numbers. */
export const isNumeric = (t) => !!t && (t.cls === 0 || t.cls === 1 || t.cls === 8 || (t.cls === 10 && isNumeric(t.base)));

/**
 * Open an HDF5 file held in memory. opts: { inflate(bytes, 'deflate', max) → Promise<Uint8Array>, maxBytes (512e6: cap on
 * one dataset read), maxObjects (2e5) }. Returns:
 *   root                      address of the root group
 *   kids(addr)                Map name → { addr } | { soft: path } | { external: true }
 *   info(addr)                { kind: 'group' | 'dataset' | 'datatype', shape, maxShape, type, layout, chunk, filters, bytes }
 *   attrs(addr)               { name: number | string | typed array | string[] | … }
 *   resolve(path, from)       address of the object at a '/'-separated path (soft links followed) or -1
 *   read(addr, { start, count, maxBytes })   Promise<{ data, shape, type }>; data is a typed array (64-bit integers as
 *                             doubles), string[], an object of member columns for compounds, or Uint8Array for raw classes
 *   walk(fn(path, addr, info), { maxDepth, maxNodes })   depth-first over the link graph, every object once
 */
export function openHDF5(u8, opts = {}) {
  const maxBytes = opts.maxBytes ?? 512e6, maxObjects = opts.maxObjects ?? 2e5, inflate = opts.inflate || ((b, _f, m) => inflateZlib(b, m));
  const N = u8.length, dv = new DataView(u8.buffer, u8.byteOffset, u8.byteLength), sb = hdf5Offset(u8);
  if (sb < 0) fail('Not an HDF5 file (signature missing).');
  const bad = (what) => fail(`The HDF5 file is truncated or corrupt (${what}).`);
  const need = (p, n) => { if (!(p >= 0) || !(n >= 0) || p + n > N) bad('a structure points beyond the end of the file'); };
  const u1 = (p) => { need(p, 1); return u8[p]; }, u2 = (p) => { need(p, 2); return dv.getUint16(p, true); }, u4 = (p) => { need(p, 4); return dv.getUint32(p, true); };
  /** n-byte little-endian unsigned integer; all bits set (undefined address / unlimited) → -1. */
  const uN = (p, n) => { need(p, n); let v = 0, m = 1, ff = n >= 2; for (let k = 0; k < n; k++) { const b = u8[p + k]; if (b !== 255) ff = false; v += b * m; m *= 256; } return ff ? -1 : v; };
  const sig = (p, s) => p >= 0 && p + 4 <= N && u8[p] === s.charCodeAt(0) && u8[p + 1] === s.charCodeAt(1) && u8[p + 2] === s.charCodeAt(2) && u8[p + 3] === s.charCodeAt(3);
  const cstr = (p, max = 65536) => { need(p, 0); let e = p; while (e < N && e - p < max && u8[e]) e++; return utf8.decode(u8.subarray(p, e)); };
  const log2 = (v) => Math.floor(Math.log2(v)), encSize = (v) => (v > 0 ? Math.floor(log2(v) / 8) + 1 : 1);

  const ver = u1(sb + 8), base = sb;
  let O, Lz, root;
  if (ver <= 1) { O = u1(sb + 13); Lz = u1(sb + 14); }
  else if (ver <= 3) { O = u1(sb + 9); Lz = u1(sb + 10); }
  else fail(`HDF5 superblock version ${ver} is not supported.`);
  if (![2, 4, 8].includes(O) || ![2, 4, 8].includes(Lz)) bad('superblock');
  const A = (p) => { const a = uN(p, O); return a < 0 ? -1 : a + base; }, LN = (p) => uN(p, Lz);
  root = ver <= 1 ? A(sb + (ver === 0 ? 24 : 28) + 5 * O) : A(sb + 12 + 3 * O);
  if (root < 0) bad('root group');
  let budget = 4e6;                                           // metadata nodes visited over the lifetime of this handle
  const tick = (n = 1) => { if ((budget -= n) < 0) fail('The HDF5 file has too much metadata to index.'); };

  // ---- object headers -------------------------------------------------------------------------------------------
  const headers = new Map();
  function header(addr) {
    let msgs = headers.get(addr);
    if (msgs) return msgs;
    if (headers.size >= maxObjects) fail('The HDF5 file holds too many objects.');
    msgs = [];
    const blocks = [];
    need(addr, 16);
    const v2 = sig(addr, 'OHDR');
    let hflags = 0, left = Infinity;
    if (v2) {
      if (u8[addr + 4] !== 2) bad('object header version');
      hflags = u8[addr + 5];
      let p = addr + 6 + (hflags & 0x20 ? 16 : 0) + (hflags & 0x10 ? 4 : 0);
      const szb = 1 << (hflags & 3), size0 = uN(p, szb);
      blocks.push([p + szb, p + szb + size0]);
    } else {
      if (u8[addr] !== 1) bad('object header');
      left = u2(addr + 2);
      blocks.push([addr + 16, addr + 16 + u4(addr + 8)]);
    }
    for (let b = 0; b < blocks.length && b < 4096; b++) {
      let [s, e] = blocks[b];
      if (e > N) e = N;
      const hs = v2 ? 4 + (hflags & 4 ? 2 : 0) : 8;
      while (s + hs <= e && left > 0 && msgs.length < 1e5) {
        const type = v2 ? u8[s] : dv.getUint16(s, true), size = dv.getUint16(s + (v2 ? 1 : 2), true), flags = u8[s + (v2 ? 3 : 4)];
        s += hs; left--;
        if (s + size > e) break;
        if (type === 0x10) {
          const ca = A(s), cl = LN(s + O);
          if (ca >= 0 && cl > 0 && ca + cl <= N) { if (v2) { if (sig(ca, 'OCHK')) blocks.push([ca + 4, ca + cl - 4]); } else blocks.push([ca, ca + cl]); }
        } else if (type !== 0) msgs.push({ type, flags, p: s, n: size });
        s += size;
        tick();
      }
    }
    headers.set(addr, msgs);
    return msgs;
  }
  const msg = (addr, type) => header(addr).find((m) => m.type === type);

  // ---- message decoders -------------------------------------------------------------------------------------------
  function dataspace(p) {
    const v = u1(p), rank = u1(p + 1), fl = u1(p + 2);
    let q, type = rank ? 1 : 0;
    if (v === 1) q = p + 8; else if (v === 2) { type = u1(p + 3); q = p + 4; } else bad('dataspace version');
    if (rank > 32) bad('dataspace rank');
    const dims = [], max = [];
    for (let k = 0; k < rank; k++) dims.push(LN(q + k * Lz));
    if (fl & 1) for (let k = 0; k < rank; k++) max.push(LN(q + (rank + k) * Lz));
    if (dims.some((d) => d < 0)) bad('dataspace extent');
    return { rank, dims, max: max.length ? max : dims.slice(), type };
  }
  function datatype(p, depth = 0) {
    if (depth > 8) fail('An HDF5 datatype is nested too deeply.');
    const cv = u1(p), cls = cv & 15, tv = cv >> 4, b0 = u1(p + 1), b1 = u1(p + 2), size = u4(p + 4), t = { cls, size };
    let q = p + 8;
    const padName = () => { let e = q; while (u1(e)) e++; const s = utf8.decode(u8.subarray(q, e)); q = tv < 3 ? q + (((e - q + 8) >> 3) << 3) : e + 1; return s; };
    if (cls === 0 || cls === 4) { t.le = !(b0 & 1); t.signed = cls === 0 && !!(b0 & 8); q += 4; }
    else if (cls === 1) { t.le = !(b0 & 1); q += 12; }
    else if (cls === 2) q += 2;
    else if (cls === 3) { t.pad = b0 & 15; t.utf8 = (b0 >> 4) === 1; }
    else if (cls === 5) q += b0;
    else if (cls === 6) {
      const n = b0 | (b1 << 8), ob = size < 256 ? 1 : size < 65536 ? 2 : size < 16777216 ? 3 : 4;
      t.members = [];
      for (let k = 0; k < n; k++) {
        const name = padName();
        let off;
        if (tv < 3) { off = u4(q); q += tv === 1 ? 32 : 4; } else { off = uN(q, ob); q += ob; }
        const mt = datatype(q, depth + 1);
        q += mt.len;
        if (off < 0 || off + mt.size > size) bad('compound member');
        t.members.push({ name, off, type: mt });
      }
    } else if (cls === 7) t.ref = b0 & 15;
    else if (cls === 8) {
      t.base = datatype(q, depth + 1); q += t.base.len;
      const n = b0 | (b1 << 8);
      t.names = [];
      for (let k = 0; k < n; k++) t.names.push(padName());
      q += n * t.base.size;
      t.le = t.base.le; t.signed = t.base.signed;
    } else if (cls === 9) { t.vlen = b0 & 15; t.utf8 = (b1 & 15) === 1; t.base = datatype(q, depth + 1); q += t.base.len; }
    else if (cls === 10) {
      const nd = u1(q);
      q += tv < 3 ? 4 : 1;
      if (nd > 32) bad('array datatype');
      t.dims = [];
      for (let k = 0; k < nd; k++) { t.dims.push(u4(q)); q += 4; }
      if (tv < 3) q += 4 * nd;
      t.base = datatype(q, depth + 1); q += t.base.len;
    } else bad('datatype class');
    t.len = q - p;
    return t;
  }
  /** Datatype of a message that may be a pointer to a committed (named) datatype. */
  function typeOf(m, depth = 0) {
    if (!(m.flags & 2)) return datatype(m.p);
    const v = u1(m.p), st = u1(m.p + 1), a = v === 1 ? A(m.p + 8) : v === 2 || (v === 3 && st !== 0) ? A(m.p + 2) : -1;
    const tm = a >= 0 && depth < 4 ? msg(a, 3) : null;
    if (!tm) fail('The HDF5 file uses a shared-message heap for its datatypes, which is not read.');
    return typeOf(tm, depth + 1);
  }
  function layout(p) {
    const v = u1(p), out = { addr: -1, size: -1 };
    if (v < 3) {
      const nd = u1(p + 1), cls = u1(p + 2);
      let q = p + 8;
      if (nd > 33) bad('layout');
      if (cls !== 0) { out.addr = A(q); q += O; }
      const d = [];
      for (let k = 0; k < nd; k++) { d.push(u4(q)); q += 4; }
      if (cls === 0) { out.cls = 'compact'; out.size = u4(q); out.addr = q + 4; need(out.addr, out.size); }
      else if (cls === 1) out.cls = 'contiguous';
      else if (cls === 2) { out.cls = 'chunked'; out.idx = 'btree1'; out.chunk = d.slice(0, -1); out.esize = d[d.length - 1]; }
      else bad('layout class');
      return out;
    }
    const cls = u1(p + 1);
    if (cls === 0) { out.cls = 'compact'; out.size = u2(p + 2); out.addr = p + 4; need(out.addr, out.size); }
    else if (cls === 1) { out.cls = 'contiguous'; out.addr = A(p + 2); out.size = LN(p + 2 + O); }
    else if (cls === 2 && v === 3) {
      const nd = u1(p + 2);
      if (nd > 33 || nd < 1) bad('layout');
      out.cls = 'chunked'; out.idx = 'btree1'; out.addr = A(p + 3);
      const d = [];
      for (let k = 0; k < nd; k++) d.push(u4(p + 3 + O + 4 * k));
      out.chunk = d.slice(0, -1); out.esize = d[nd - 1];
    } else if (cls === 2 && (v === 4 || v === 5)) {            // version 5 (HDF5 2.x) only widens the stored chunk sizes, which the index entries describe themselves
      const fl = u1(p + 2), nd = u1(p + 3), eb = u1(p + 4);
      if (nd > 33 || nd < 1 || eb < 1 || eb > 8) bad('layout');
      let q = p + 5;
      const d = [];
      for (let k = 0; k < nd; k++) { d.push(uN(q, eb)); q += eb; }
      out.cls = 'chunked'; out.chunk = d.slice(0, -1); out.esize = d[nd - 1];
      const it = u1(q++);
      out.idx = { 1: 'single', 2: 'implicit', 3: 'farray', 4: 'earray', 5: 'btree2' }[it] || bad('chunk index type');
      if (it === 1) { if (fl & 2) { out.fsize = LN(q); out.fmask = u4(q + Lz); q += Lz + 4; } }
      else if (it === 3) q += 1; else if (it === 4) q += 5; else if (it === 5) q += 6;
      out.addr = A(q);
    } else if (cls === 3) out.cls = 'virtual';
    else bad('layout class');
    return out;
  }
  function filters(p) {
    const v = u1(p), n = u1(p + 1), out = [];
    let q = p + (v === 1 ? 8 : 2);
    if (v !== 1 && v !== 2) bad('filter pipeline');
    for (let k = 0; k < n && k < 32; k++) {
      const id = u2(q);
      let nl = 0;
      q += 2;
      if (v === 1 || id >= 256) { nl = u2(q); q += 2; }
      const ncd = u2(q + 2);
      q += 4;
      q += v === 1 ? (nl + 7) & ~7 : nl;
      const cd = [];
      for (let c = 0; c < ncd && c < 64; c++) cd.push(u4(q + 4 * c));
      q += 4 * ncd + (v === 1 && ncd % 2 ? 4 : 0);
      out.push({ id, cd });
    }
    return out;
  }
  function fillValue(addr) {
    const m5 = msg(addr, 5), m4 = msg(addr, 4);
    if (m5) {
      const v = u1(m5.p);
      if (v === 3) { const fl = u1(m5.p + 1); if (fl & 0x20) { const n = u4(m5.p + 2); need(m5.p + 6, n); return u8.subarray(m5.p + 6, m5.p + 6 + n); } return null; }
      if (v === 1 || (v === 2 && u1(m5.p + 3))) { const n = u4(m5.p + 4); if (n > 0 && n <= 65536) { need(m5.p + 8, n); return u8.subarray(m5.p + 8, m5.p + 8 + n); } }
      return null;
    }
    if (m4) { const n = u4(m4.p); if (n > 0 && n <= 65536 && n + 4 <= m4.n) return u8.subarray(m4.p + 4, m4.p + 4 + n); }
    return null;
  }

  // ---- heaps and B-trees ---------------------------------------------------------------------------------------------
  /** v2 B-tree: visit(recordOffset) for every record; returns the record size. */
  function btree2(addr, visit) {
    if (!sig(addr, 'BTHD')) bad('v2 B-tree header');
    const nodeSize = u4(addr + 6), recSize = u2(addr + 10), depth = u2(addr + 12), rootA = A(addr + 16), nroot = u2(addr + 16 + O);
    if (rootA < 0 || !nroot) return recSize;
    if (!(recSize > 0) || depth > 32 || nodeSize < 16) bad('v2 B-tree header');
    const maxLeaf = Math.floor((nodeSize - 10) / recSize), nrecSize = encSize(maxLeaf), cum = [maxLeaf], cumSize = [0];
    for (let u = 1; u <= depth; u++) { const ptr = O + nrecSize + cumSize[u - 1], mx = Math.floor((nodeSize - (10 + ptr)) / (recSize + ptr)); cum[u] = (mx + 1) * cum[u - 1] + mx; cumSize[u] = encSize(cum[u]); }
    const node = (a, n, d) => {
      tick(n + 1);
      if (!sig(a, d ? 'BTIN' : 'BTLF')) bad('v2 B-tree node');
      need(a + 6, n * recSize);
      for (let k = 0; k < n; k++) visit(a + 6 + k * recSize);
      if (!d) return;
      let p = a + 6 + n * recSize;
      for (let k = 0; k <= n; k++) { const ca = A(p), cn = uN(p + O, nrecSize); p += O + nrecSize + (d > 1 ? cumSize[d - 1] : 0); if (ca >= 0 && cn > 0) node(ca, cn, d - 1); }
    };
    node(rootA, nroot, depth);
    return recSize;
  }
  /** Fractal heap accessor: get(idOffset) → bytes of a managed or tiny object, or null (huge / filtered). */
  function fheap(addr) {
    if (!sig(addr, 'FRHP')) bad('fractal heap');
    const idLen = u2(addr + 5), filtLen = u2(addr + 7), maxMan = u4(addr + 10);
    let p = addr + 14 + 10 * Lz + 2 * O;
    const width = u2(p), start = LN(p + 2), maxDirect = LN(p + 2 + Lz), bits = u2(p + 2 + 2 * Lz), rootA = A(p + 6 + 2 * Lz), curRows = u2(p + 6 + 2 * Lz + O);
    if (!(width > 0 && start > 0 && maxDirect >= start) || bits > 64) bad('fractal heap');
    const offSize = (bits + 7) >> 3, lenSize = Math.min((log2(maxDirect) + 7) >> 3, encSize(maxMan)), maxDRows = log2(maxDirect) - log2(start) + 2;
    const locate = (ia, nrows, boff, target, depth) => {
      if (depth > 16 || !sig(ia, 'FHIB')) bad('fractal heap block');
      tick();
      const q = ia + 5 + O + offSize, nd = Math.min(nrows, maxDRows) * width;
      let off = boff;
      for (let r = 0; r < nrows && r < 64; r++) {
        const size = r < 2 ? start : start * 2 ** (r - 1);
        if (target >= off + size * width) { off += size * width; continue; }
        const c = Math.floor((target - off) / size);
        off += c * size;
        if (r < maxDRows) return [A(q + (r * width + c) * O), off, size];
        return locate(A(q + (nd + (r - maxDRows) * width + c) * O), log2(size) - log2(start * width) + 1, off, target, depth + 1);
      }
      return [-1, 0, 0];
    };
    return {
      idLen,
      get(id) {
        const b = u1(id), kind = (b >> 4) & 3;
        if (kind === 2) { const n = (b & 15) + 1; need(id + 1, n); return u8.subarray(id + 1, id + 1 + n); }
        if (kind !== 0 || filtLen || rootA < 0) return null;
        const off = uN(id + 1, offSize), len = uN(id + 1 + offSize, lenSize), [ba, boff, bsize] = curRows ? locate(rootA, curRows, 0, off, 0) : [rootA, 0, start];
        if (ba < 0 || off < boff || off - boff + len > bsize) return null;
        need(ba + off - boff, len);
        return u8.subarray(ba + off - boff, ba + off - boff + len);
      },
    };
  }
  const gcols = new Map();
  /** Global heap object (collection address, index) → bytes or null. */
  function gheap(addr, idx) {
    let col = gcols.get(addr);
    if (!col) {
      col = new Map();
      if (sig(addr, 'GCOL')) {
        const end = Math.min(N, addr + LN(addr + 8));
        for (let p = addr + 16, k = 0; p + 8 + Lz <= end && k < 1e6; k++) {
          const id = u2(p), n = LN(p + 8);
          if (!id || n < 0 || p + 8 + Lz + n > end) break;
          col.set(id, [p + 8 + Lz, n]);
          p += 8 + Lz + ((n + 7) & ~7);
        }
      }
      if (gcols.size < 1e5) gcols.set(addr, col);
    }
    const e = col.get(idx);
    return e ? u8.subarray(e[0], e[0] + e[1]) : null;
  }

  // ---- groups ------------------------------------------------------------------------------------------------------
  function linkMsg(b, o, out) {
    const d = new DataView(b.buffer, b.byteOffset, b.byteLength), n = b.length;
    if (o + 3 > n || b[o] !== 1) return;
    const fl = b[o + 1], ls = 1 << (fl & 3);
    let q = o + 2, lt = 0;
    if (fl & 8) lt = b[q++];
    if (fl & 4) q += 8;
    if (fl & 16) q++;
    if (q + ls > n) return;
    let nl = 0;
    for (let k = ls - 1; k >= 0; k--) nl = nl * 256 + b[q + k];
    q += ls;
    if (q + nl > n) return;
    const name = utf8.decode(b.subarray(q, q + nl));
    q += nl;
    if (lt === 0 && q + O <= n) { let a = 0, ff = true; for (let k = O - 1; k >= 0; k--) { a = a * 256 + b[q + k]; if (b[q + k] !== 255) ff = false; } if (!ff) out.set(name, { addr: a + base }); }
    else if (lt === 1 && q + 2 <= n) { const sl = d.getUint16(q, true); if (q + 2 + sl <= n) out.set(name, { soft: utf8.decode(b.subarray(q + 2, q + 2 + sl)) }); }
    else out.set(name, { external: true });
  }
  const kidCache = new Map();
  function kids(addr) {
    let out = kidCache.get(addr);
    if (out) return out;
    out = new Map();
    for (const m of header(addr)) {
      if (out.size > 1e6) fail('An HDF5 group holds too many links.');
      if (m.type === 6) linkMsg(u8, m.p, out);
      else if (m.type === 0x11) {
        const bt = A(m.p), hp = A(m.p + O);
        if (bt < 0 || !sig(hp, 'HEAP')) continue;
        const hd = A(hp + 8 + 2 * Lz);
        const walk = (b, depth) => {
          if (depth > 32 || !sig(b, 'TREE') || u1(b + 4) !== 0) bad('group B-tree');
          tick();
          const level = u1(b + 5), n = u2(b + 6), p = b + 8 + 2 * O;
          for (let k = 0; k < n; k++) {
            const c = A(p + Lz + k * (Lz + O));
            if (c < 0) continue;
            if (level > 0) { walk(c, depth + 1); continue; }
            if (!sig(c, 'SNOD')) bad('symbol table node');
            const ns = u2(c + 6);
            tick(ns);
            for (let s = 0, e = c + 8; s < ns; s++, e += 2 * O + 24) {
              const name = cstr(hd + uN(e, O)), oa = A(e + O), ct = u4(e + 2 * O);
              if (ct === 2) out.set(name, { soft: cstr(hd + u4(e + 2 * O + 8)) }); else if (oa >= 0) out.set(name, { addr: oa });
            }
          }
        };
        walk(bt, 0);
      } else if (m.type === 2) {
        const fl = u1(m.p + 1), q = m.p + 2 + (fl & 1 ? 8 : 0), ha = A(q), ba = A(q + O);
        if (ha < 0 || ba < 0) continue;
        const heap = fheap(ha);
        btree2(ba, (r) => { const o = heap.get(r + 4); if (o) linkMsg(o, 0, out); });
      }
    }
    kidCache.set(addr, out);
    return out;
  }
  function resolve(path, from = root, depth = 0) {
    let cur = String(path).startsWith('/') ? root : from;
    for (const part of String(path).split('/')) {
      if (!part || part === '.') continue;
      const l = kids(cur).get(part);
      if (!l) return -1;
      if (l.addr !== undefined) cur = l.addr;
      else if (l.soft !== undefined && depth < 16) { cur = resolve(l.soft, cur, depth + 1); if (cur < 0) return -1; }
      else return -1;
    }
    return cur;
  }

  // ---- value decoding ---------------------------------------------------------------------------------------------
  const swap = (b, sz) => { const c = b.slice(); for (let i = 0; i + sz <= c.length; i += sz) for (let a = 0, z = sz - 1; a < z; a++, z--) { const t = c[i + a]; c[i + a] = c[i + z]; c[i + z] = t; } return c; };
  const aligned = (b, sz, le) => (le !== HOST_LE && sz > 1 ? swap(b, sz) : b.byteOffset % sz || !(b.buffer instanceof ArrayBuffer) ? b.slice() : b);
  /** raw bytes of n elements of type t → values. */
  function decode(raw, t, n, depth = 0) {
    if (depth > 8) fail('An HDF5 datatype is nested too deeply.');
    const sz = t.size;
    if (raw.length < n * sz) bad('dataset shorter than its extent');
    if (t.cls === 0 || t.cls === 4 || t.cls === 8) {
      const s = t.cls === 8 ? t.base.signed : t.signed, C = { 1: s ? Int8Array : Uint8Array, 2: s ? Int16Array : Uint16Array, 4: s ? Int32Array : Uint32Array, 8: s ? BigInt64Array : BigUint64Array }[sz];
      if (!C) { const out = new Float64Array(n); for (let i = 0; i < n; i++) { let v = 0; for (let k = 0; k < sz; k++) v += raw[i * sz + (t.le ? k : sz - 1 - k)] * 256 ** k; out[i] = v; } return out; }
      const b = aligned(raw.subarray(0, n * sz), sz, t.le), a = new C(b.buffer, b.byteOffset, n);
      if (sz < 8) return a;
      const f = new Float64Array(n);
      for (let i = 0; i < n; i++) f[i] = Number(a[i]);
      return f;
    }
    if (t.cls === 1) {
      if (sz === 4 || sz === 8) { const b = aligned(raw.subarray(0, n * sz), sz, t.le); return new (sz === 4 ? Float32Array : Float64Array)(b.buffer, b.byteOffset, n); }
      if (sz !== 2) fail(`HDF5 floating-point numbers of ${sz * 8} bits are not supported.`);
      const out = new Float32Array(n);
      for (let i = 0; i < n; i++) { const h = t.le ? raw[2 * i] | (raw[2 * i + 1] << 8) : raw[2 * i + 1] | (raw[2 * i] << 8), e = (h >> 10) & 31, f = h & 1023, v = e === 0 ? f * 2 ** -24 : e === 31 ? (f ? NaN : Infinity) : (1 + f / 1024) * 2 ** (e - 15); out[i] = h & 0x8000 ? -v : v; }
      return out;
    }
    if (t.cls === 3) {
      const out = new Array(n), dec = t.utf8 ? utf8 : latin1;
      for (let i = 0; i < n; i++) { let e = i * sz + sz; const s0 = i * sz; if (t.pad === 2) { while (e > s0 && (raw[e - 1] === 32 || raw[e - 1] === 0)) e--; } else { const z = raw.indexOf(0, s0); if (z >= 0 && z < e) e = z; } out[i] = dec.decode(raw.subarray(s0, e)); }
      return out;
    }
    if (t.cls === 9) {
      const out = new Array(n), d = new DataView(raw.buffer, raw.byteOffset, raw.byteLength);
      if (sz < 8 + O) bad('variable-length datatype');
      for (let i = 0; i < n; i++) {
        const len = d.getUint32(i * sz, true);
        let a = 0, ff = true;
        for (let k = O - 1; k >= 0; k--) { const b = raw[i * sz + 4 + k]; a = a * 256 + b; if (b !== 255) ff = false; }
        const obj = ff || !len ? null : gheap(a + base, d.getUint32(i * sz + 4 + O, true));
        tick();
        if (t.vlen === 1) out[i] = obj ? (t.utf8 ? utf8 : latin1).decode(obj.subarray(0, Math.min(len, obj.length))).replace(/\0[\s\S]*$/, '') : '';
        else out[i] = obj && obj.length >= len * t.base.size ? decode(obj, t.base, len, depth + 1) : [];
      }
      return out;
    }
    if (t.cls === 7) {
      const out = new Float64Array(n);
      for (let i = 0; i < n; i++) { let a = 0, ff = true; for (let k = Math.min(O, sz) - 1; k >= 0; k--) { const b = raw[i * sz + k]; a = a * 256 + b; if (b !== 255) ff = false; } out[i] = ff || t.ref !== 0 ? -1 : a + base; }
      return out;
    }
    if (t.cls === 10) return decode(raw, t.base, n * t.dims.reduce((a, b) => a * b, 1), depth + 1);
    if (t.cls === 6) {
      const out = {};
      for (const m of t.members.slice(0, 256)) {
        const ms = m.type.size, col = new Uint8Array(n * ms);
        for (let i = 0; i < n; i++) col.set(raw.subarray(i * sz + m.off, i * sz + m.off + ms), i * ms);
        try { out[m.name] = decode(col, m.type, n, depth + 1); } catch (e) { if (!e || !e.user) throw e; }
      }
      return out;
    }
    return raw.slice(0, n * sz);
  }
  const scalarise = (v, shape) => (shape.length === 0 && v && typeof v === 'object' && v.length === 1 ? v[0] : v);

  // ---- attributes -------------------------------------------------------------------------------------------------
  function attrMsg(b, o, out) {
    const d = new DataView(b.buffer, b.byteOffset, b.byteLength), n = b.length;
    if (o + 8 > n) return;
    const v = b[o], fl = v > 1 ? b[o + 1] : 0, ns = d.getUint16(o + 2, true), ts = d.getUint16(o + 4, true), ss = d.getUint16(o + 6, true), pad = v === 1 ? (x) => (x + 7) & ~7 : (x) => x;
    if (v < 1 || v > 3 || fl & 3) return;                    // shared datatype / dataspace of an attribute: skipped
    let q = o + 8 + (v === 3 ? 1 : 0);
    if (q + pad(ns) + pad(ts) + pad(ss) > n) return;
    const name = utf8.decode(b.subarray(q, q + ns)).replace(/\0[\s\S]*$/, '');
    q += pad(ns);
    // the embedded messages are decoded through the file-level readers, so they need absolute offsets
    const abs = b === u8 ? 0 : b.byteOffset - u8.byteOffset;
    let t, s;
    try { t = datatype(abs + q); s = dataspace(abs + q + pad(ts)); } catch (e) { if (e && e.user) return; throw e; }
    q += pad(ts) + pad(ss);
    const cnt = s.type === 2 ? 0 : s.dims.reduce((a, c) => a * c, 1);
    if (!(cnt * t.size <= n - q) || cnt > 4e6) return;
    try { out[name] = cnt ? scalarise(decode(b.subarray(q, q + cnt * t.size), t, cnt), s.dims) : null; } catch (e) { if (!e || !e.user) throw e; }
  }
  const attrCache = new Map();
  function attrs(addr) {
    let out = attrCache.get(addr);
    if (out) return out;
    out = Object.create(null);
    for (const m of header(addr)) {
      if (m.type === 0x0c) attrMsg(u8, m.p, out);
      else if (m.type === 0x15) {
        const fl = u1(m.p + 1), q = m.p + 2 + (fl & 1 ? 2 : 0), ha = A(q), ba = A(q + O);
        if (ha < 0 || ba < 0) continue;
        const heap = fheap(ha);
        btree2(ba, (r) => { const o = heap.get(r); if (o) attrMsg(o, 0, out); });
      }
    }
    attrCache.set(addr, out);
    return out;
  }

  // ---- objects and data -----------------------------------------------------------------------------------------------
  const infoCache = new Map();
  function info(addr) {
    let o = infoCache.get(addr);
    if (o) return o;
    const ms = msg(addr, 1), mt = msg(addr, 3), ml = msg(addr, 8), mf = msg(addr, 0x0b);
    if (ms && mt && ml) {
      const s = dataspace(ms.p), lay = layout(ml.p);
      let type = null, note;
      try { type = typeOf(mt); } catch (e) { if (!e || !e.user) throw e; note = e.message; }
      o = { kind: 'dataset', shape: s.type === 2 ? [0] : s.dims, maxShape: s.max, space: s, type, layout: lay.cls, chunk: lay.chunk || null, lay, filters: mf ? filters(mf.p) : [], external: !!msg(addr, 7), note };
      o.count = o.shape.reduce((a, b) => a * b, 1);
      o.bytes = type ? o.count * type.size : 0;
    } else if (mt && !ms) o = { kind: 'datatype', type: typeOf(mt) };
    else o = { kind: 'group' };
    infoCache.set(addr, o);
    return o;
  }
  /** Chunk records { addr, size, mask, off[] (element offsets) } of a chunked dataset. */
  function chunkList(inf) {
    const { lay, shape, maxShape } = inf, rank = shape.length, cd = lay.chunk, out = [], cbytes = cd.reduce((a, b) => a * b, 1) * lay.esize;
    if (cd.length !== rank || cd.some((c) => !(c > 0))) bad('chunk dimensions');
    const push = (addr, size, mask, off) => { if (out.length >= 4e6) fail('The HDF5 dataset has too many chunks.'); if (addr >= 0) out.push({ addr, size, mask, off }); };
    if (lay.addr < 0) return out;
    // linear chunk index → element offsets, by the "maximal extent" chunk counts (an unlimited dimension runs slowest)
    const nmax = maxShape.map((m, d) => (m < 0 ? -1 : Math.ceil(m / cd[d]))), unl = nmax.indexOf(-1), order = unl > 0 ? [unl, ...nmax.map((_, d) => d).filter((d) => d !== unl)] : nmax.map((_, d) => d);
    const fromIndex = (i) => { const off = new Array(rank); for (let k = rank - 1; k >= 0; k--) { const d = order[k], m = nmax[d]; if (m < 0) { off[d] = i * cd[d]; i = 0; } else { off[d] = (i % m) * cd[d]; i = Math.floor(i / m); } } return off; };
    const elem = (p, k, filtered, es) => (filtered ? push(A(p), uN(p + O, es - O - 4), u4(p + es - 4), fromIndex(k)) : push(A(p), cbytes, 0, fromIndex(k)));
    if (lay.idx === 'btree1') {
      const ks = 8 + 8 * (rank + 1);
      const walk = (b, depth) => {
        if (depth > 32 || !sig(b, 'TREE') || u1(b + 4) !== 1) bad('chunk B-tree');
        const level = u1(b + 5), n = u2(b + 6), p = b + 8 + 2 * O;
        tick(n + 1);
        for (let k = 0; k < n; k++) {
          const e = p + k * (ks + O), c = A(e + ks);
          if (c < 0) continue;
          if (level > 0) { walk(c, depth + 1); continue; }
          const off = [];
          for (let d = 0; d < rank; d++) off.push(uN(e + 8 + 8 * d, 8));
          push(c, u4(e), u4(e + 4), off);
        }
      };
      walk(lay.addr, 0);
    } else if (lay.idx === 'single') push(lay.addr, lay.fsize ?? cbytes, lay.fmask || 0, new Array(rank).fill(0));
    else if (lay.idx === 'implicit') {
      const n = shape.reduce((a, s, d) => a * Math.ceil(s / cd[d]), 1);
      for (let i = 0; i < n; i++) push(lay.addr + i * cbytes, cbytes, 0, fromIndex(i));
    } else if (lay.idx === 'farray') {
      const h = lay.addr;
      if (!sig(h, 'FAHD')) bad('fixed-array chunk index');
      const filtered = u1(h + 5) === 1, es = u1(h + 6), pb = u1(h + 7), n = LN(h + 8), db = A(h + 8 + Lz);
      if (db < 0) return out;
      if (!sig(db, 'FADB') || es < O || pb > 30 || !(n >= 0 && n <= 4e6)) bad('fixed-array chunk index');
      const pn = 2 ** pb;
      let p = db + 6 + O;
      tick(n);
      if (n <= pn) for (let k = 0; k < n; k++) elem(p + k * es, k, filtered, es);
      else {
        const np = Math.ceil(n / pn), bm = p;
        p += ((np + 7) >> 3) + 4;
        for (let g = 0; g < np; g++) if (u1(bm + (g >> 3)) & (0x80 >> (g & 7))) for (let k = g * pn, e = Math.min(n, k + pn), q = p + g * (pn * es + 4); k < e; k++, q += es) elem(q, k, filtered, es);
      }
    } else if (lay.idx === 'earray') {
      const h = lay.addr;
      if (!sig(h, 'EAHD')) bad('extensible-array chunk index');
      const filtered = u1(h + 5) === 1, es = u1(h + 6), maxBits = u1(h + 7), nIdx = u1(h + 8), minElm = u1(h + 9), minPtr = u1(h + 10), pageBits = u1(h + 11), ib = A(h + 12 + 6 * Lz);
      if (ib < 0) return out;
      if (!sig(ib, 'EAIB') || es < O || !(minElm > 0 && minPtr > 1) || maxBits > 64 || pageBits > 30) bad('extensible-array chunk index');
      const offSize = (maxBits + 7) >> 3, nsb = 1 + maxBits - log2(minElm), ibSb = 2 * log2(minPtr), nDb = 2 * (minPtr - 1), pageN = 2 ** pageBits;
      let p = ib + 6 + O, k = 0;
      for (let i = 0; i < nIdx; i++, k++) elem(p + i * es, k, filtered, es);
      p += nIdx * es;
      const dblock = (a, ne, k0, bitmap) => {
        if (a < 0) return;
        if (!sig(a, 'EADB')) bad('extensible-array data block');
        tick(ne);
        const q = a + 6 + O + offSize;
        if (ne <= pageN) { for (let i = 0; i < ne; i++) elem(q + i * es, k0 + i, filtered, es); return; }
        for (let g = 0, np = Math.ceil(ne / pageN); g < np; g++) if (bitmap < 0 || u1(bitmap[0] + ((bitmap[1] + g) >> 3)) & (0x80 >> ((bitmap[1] + g) & 7))) for (let i = g * pageN, e = Math.min(ne, i + pageN), r = q + 4 + g * (pageN * es + 4); i < e; i++, r += es) elem(r, k0 + i, filtered, es);
      };
      const info2 = (u) => [2 ** (u >> 1), 2 ** ((u + 1) >> 1) * minElm];
      let di = 0;
      for (let u = 0; u < Math.min(ibSb, nsb); u++) { const [nd, ne] = info2(u); for (let j = 0; j < nd && di < nDb; j++, di++, k += ne) dblock(A(p + di * O), ne, k, -1); }
      p += nDb * O;
      for (let u = ibSb; u < nsb && u < 128; u++) {
        const [nd, ne] = info2(u), sa = A(p + (u - ibSb) * O);
        if (k > 4e6 * 64) break;
        if (sa >= 0) {
          if (!sig(sa, 'EASB')) bad('extensible-array super block');
          const paged = ne > pageN, bmSize = paged ? (Math.ceil(ne / pageN) + 7) >> 3 : 0, q = sa + 6 + O + offSize;
          tick(nd);
          for (let j = 0; j < nd; j++) dblock(A(q + nd * bmSize + j * O), ne, k + j * ne, paged ? [q, j * Math.ceil(ne / pageN)] : -1);   // page-initialised bits run on from one data block to the next
        }
        k += nd * ne;
      }
    } else if (lay.idx === 'btree2') {
      let rs = 0;
      const recs = [];
      rs = btree2(lay.addr, (r) => recs.push(r));
      const filtered = u1(lay.addr + 5) === 11, cs = rs - O - 4 - 8 * rank;
      if (filtered ? cs < 1 || cs > 8 : rs < O + 8 * rank) bad('v2 B-tree chunk index');
      for (const r of recs) { const q = r + O + (filtered ? cs + 4 : 0), off = []; for (let d = 0; d < rank; d++) off.push(uN(q + 8 * d, 8) * cd[d]); push(A(r), filtered ? uN(r + O, cs) : cbytes, filtered ? u4(r + O + cs) : 0, off); }
    }
    return out;
  }
  /** Copy the part of a block (origin bo, extent bd, elements of es bytes at src[sp…]) that falls in the selection. */
  function blit(src, sp, bo, bd, out, start, count, es) {
    const r = count.length, lo = new Array(r), hi = new Array(r);
    for (let d = 0; d < r; d++) { lo[d] = Math.max(start[d], bo[d]); hi[d] = Math.min(start[d] + count[d], bo[d] + bd[d]); if (lo[d] >= hi[d]) return; }
    const ss = new Array(r), os = new Array(r);
    ss[r - 1] = os[r - 1] = es;
    for (let d = r - 2; d >= 0; d--) { ss[d] = ss[d + 1] * bd[d + 1]; os[d] = os[d + 1] * count[d + 1]; }
    const i = lo.slice(), row = (hi[r - 1] - lo[r - 1]) * es;
    for (;;) {
      let a = sp, b = 0;
      for (let d = 0; d < r; d++) { a += (i[d] - bo[d]) * ss[d]; b += (i[d] - start[d]) * os[d]; }
      if (a + row > src.length) bad('dataset shorter than its extent');
      out.set(src.subarray(a, a + row), b);
      let d = r - 2;
      for (; d >= 0; d--) { if (++i[d] < hi[d]) break; i[d] = lo[d]; }
      if (d < 0) break;
    }
  }
  async function read(addr, sel = {}) {
    const inf = info(addr);
    if (inf.kind !== 'dataset') fail('The HDF5 object is not a dataset.');
    if (!inf.type) fail(inf.note || 'The HDF5 dataset has an unsupported datatype.');
    if (inf.external || inf.layout === 'virtual') fail('The HDF5 dataset keeps its values in other files (external or virtual storage), which is not read.');
    const t = inf.type, es = t.size, rank0 = inf.shape.length, shape = rank0 ? inf.shape : [1];
    const start = shape.map((_, d) => Math.max(0, Math.min(shape[d], (sel.start && sel.start[d]) | 0))), count = shape.map((s, d) => Math.max(0, Math.min(s - start[d], sel.count && sel.count[d] !== undefined ? sel.count[d] | 0 : s)));
    const n = inf.space.type === 2 ? 0 : count.reduce((a, b) => a * b, 1), total = n * es, cap = Math.min(maxBytes, sel.maxBytes ?? Infinity);
    if (!(total <= cap)) fail(`An HDF5 dataset of ${Math.round(total / 1e6)} MB is too large to read.`);
    const out = new Uint8Array(total), fv = fillValue(addr), lay = inf.lay;
    if (fv && fv.length === es && fv.some((b) => b)) for (let i = 0; i < total; i += es) out.set(fv, i);
    if (n && (lay.cls === 'compact' || lay.cls === 'contiguous')) {
      if (lay.addr >= 0) { const avail = lay.cls === 'compact' ? lay.size : N - lay.addr; if (inf.bytes > avail) bad('dataset shorter than its extent'); blit(u8, lay.addr, shape.map(() => 0), shape, out, start, count, es); }
    } else if (n) {
      const cd = rank0 ? lay.chunk : [1], cbytes = cd.reduce((a, b) => a * b, 1) * lay.esize;
      if (lay.esize !== es || !(cbytes > 0 && cbytes <= maxBytes)) bad('chunk size');
      for (const f of inf.filters) if (![1, 2, 3].includes(f.id)) fail(`The HDF5 dataset is compressed with the ${FILTER_NAME[f.id] || 'filter ' + f.id} filter, which is not read (deflate and shuffle are). Rewrite it with gzip compression: h5repack -f GZIP=4 in.h5 out.h5.`);
      for (const c of rank0 ? chunkList(inf) : [{ addr: lay.addr, size: cbytes, mask: 0, off: [0] }]) {
        if (c.off.some((o, d) => !(o >= 0) || o >= start[d] + count[d] || o + cd[d] <= start[d])) continue;
        need(c.addr, c.size);
        let b = u8.subarray(c.addr, c.addr + c.size);
        for (let k = inf.filters.length - 1; k >= 0; k--) {
          const f = inf.filters[k];
          if (c.mask & (1 << k)) continue;
          if (f.id === 1) b = await inflate(b, 'deflate', cbytes + 64);
          else if (f.id === 3) b = b.subarray(0, Math.max(0, b.length - 4));
          else { const sz = f.cd[0] || es, m = Math.floor(b.length / sz), u = new Uint8Array(b.length); if (sz > 1) { for (let j = 0; j < sz; j++) for (let i = 0, o = j * m; i < m; i++) u[i * sz + j] = b[o + i]; u.set(b.subarray(m * sz), m * sz); b = u; } }
        }
        if (b.length < cbytes) bad('chunk shorter than its size');
        blit(b, 0, c.off, cd, out, start, count, es);
      }
    }
    return { data: decode(out, t, n), shape: rank0 ? count.concat(t.cls === 10 ? t.dims : []) : t.cls === 10 ? t.dims.slice() : [], type: t };
  }
  function walk(fn, { maxDepth = 24, maxNodes = 1e5 } = {}) {
    const seen = new Set([root]), stack = [['', root, 0]];
    let n = 0;
    while (stack.length) {
      const [path, addr, depth] = stack.pop();
      if (++n > maxNodes) fail('The HDF5 file holds too many objects.');
      let inf;
      try { inf = info(addr); } catch (e) { if (!e || !e.user) throw e; continue; }
      if (fn(path || '/', addr, inf) === false || inf.kind !== 'group' || depth >= maxDepth) continue;
      const list = [...kids(addr)].filter(([, l]) => l.addr !== undefined && !seen.has(l.addr));
      for (let k = list.length - 1; k >= 0; k--) { seen.add(list[k][1].addr); stack.push([`${path}/${list[k][0]}`, list[k][1].addr, depth + 1]); }
    }
  }
  return { version: ver, root, userBlock: sb, kids, info, attrs, resolve, read, walk };
}
