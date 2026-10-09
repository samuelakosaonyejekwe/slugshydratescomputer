// LAZ (LASzip-compressed LAS) reader: point coordinates.
//
// Read:   LAS 1.0-1.4 headers, the LASzip VLR (record 22204), the "pointwise" and "pointwise chunked" compressors with
//         the arithmetic coder and the version-2 item codecs POINT10, GPSTIME11, RGB12 and BYTE (point formats 0-3), and
//         the "layered chunked" compressor of LAS 1.4 (point formats 6-10, item version 3): per chunk the raw first
//         point, the point count and layer sizes, then the "channel, returns and XY" layer and the Z layer with one
//         context per scanner channel; the chunk table with fixed or variable chunk sizes (COPC files read as LAZ).
//         Only x, y, z are kept (scaled to real coordinates); the other fields are decoded as far as needed to stay in
//         step with the stream. Clouds above the cap are evenly sub-sampled.
// Not read: wave packets of formats 4, 5 (rejected), version-1 item codecs of LASzip releases before 2.0 (rejected), the
//         LASzip "compatibility mode" re-mapping, and every attribute other than the coordinates. Formats 9 and 10 share
//         the coordinate layers of 6-8 but were not checked against a file.
// File content is untrusted: counts and offsets are bounds-checked and reading past the end of the data stops the chunk.

function fail(msg) { const e = new Error(msg); e.user = true; throw e; }
const HOW = 'Decompress it first: laszip -i in.laz -o out.las, or pdal translate in.laz out.las (the LAS file is read directly).';

// ---- adaptive arithmetic decoder (the coder of LASzip) -------------------------------------------------------------
class SymbolModel {
  constructor(n) {
    this.n = n; this.dist = new Uint32Array(n); this.count = new Uint32Array(n); this.table = null; this.tshift = 0; this.tsize = 0;
    if (n > 16) { let tb = 3; while (n > 1 << (tb + 2)) tb++; this.tsize = 1 << tb; this.tshift = 15 - tb; this.table = new Uint32Array(this.tsize + 2); }
    this.init();
  }
  init() { this.total = 0; this.cycle = this.n; this.count.fill(1); this.update(); this.left = this.cycle = (this.n + 6) >> 1; }
  update() {
    const { n, count, dist, table } = this;
    if ((this.total += this.cycle) > 32768) { this.total = 0; for (let k = 0; k < n; k++) this.total += count[k] = (count[k] + 1) >> 1; }
    const scale = Math.floor(0x80000000 / this.total);
    let sum = 0, s = 0;
    if (!table) for (let k = 0; k < n; k++) { dist[k] = (scale * sum) >>> 16; sum += count[k]; }
    else {
      for (let k = 0; k < n; k++) { dist[k] = (scale * sum) >>> 16; sum += count[k]; const w = dist[k] >>> this.tshift; while (s < w) table[++s] = k - 1; }
      table[0] = 0;
      while (s <= this.tsize) table[++s] = n - 1;
    }
    this.cycle = (5 * this.cycle) >> 2;
    const max = (n + 6) << 3;
    if (this.cycle > max) this.cycle = max;
    this.left = this.cycle;
  }
}
class BitModel {
  constructor() { this.c0 = 1; this.c = 2; this.p0 = 4096; this.cycle = this.left = 4; }
  update() {
    if ((this.c += this.cycle) > 8192) { this.c = (this.c + 1) >> 1; this.c0 = (this.c0 + 1) >> 1; if (this.c0 === this.c) this.c++; }
    this.p0 = (this.c0 * Math.floor(0x80000000 / this.c)) >>> 18;
    this.cycle = (5 * this.cycle) >> 2;
    if (this.cycle > 64) this.cycle = 64;
    this.left = this.cycle;
  }
}
class Decoder {
  constructor(u8, p, end) { this.u8 = u8; this.p = p + 4; this.end = end; this.over = 0; this.value = ((u8[p] << 24) | (u8[p + 1] << 16) | (u8[p + 2] << 8) | u8[p + 3]) >>> 0; this.length = 0xffffffff; }
  renorm() { do { this.value = ((this.value << 8) | (this.p < this.end ? this.u8[this.p++] : (this.over++, 0))) >>> 0; this.length = (this.length << 8) >>> 0; } while (this.length < 0x1000000); }
  bit(m) {
    const x = m.p0 * (this.length >>> 13), sym = this.value >= x ? 1 : 0;
    if (sym) { this.value -= x; this.length -= x; } else { this.length = x; m.c0++; }
    if (this.length < 0x1000000) this.renorm();
    if (--m.left === 0) m.update();
    return sym;
  }
  symbol(m) {
    const dist = m.dist;
    let sym, x, y = this.length, n;
    if (m.table) {
      const len = (this.length >>>= 15), dv = Math.floor(this.value / len), t = dv >>> m.tshift;
      sym = m.table[t]; n = m.table[t + 1] + 1;
      while (n > sym + 1) { const k = (sym + n) >>> 1; if (dist[k] > dv) n = k; else sym = k; }
      x = dist[sym] * len;
      if (sym !== m.n - 1) y = dist[sym + 1] * len;
    } else {
      const len = (this.length >>>= 15);
      x = sym = 0; n = m.n;
      let k = n >>> 1;
      do { const z = len * dist[k]; if (z > this.value) { n = k; y = z; } else { sym = k; x = z; } } while ((k = (sym + n) >>> 1) !== sym);
    }
    this.value -= x; this.length = y - x;
    if (this.length < 0x1000000) this.renorm();
    m.count[sym]++;
    if (--m.left === 0) m.update();
    return sym;
  }
  bits(b) {
    if (b > 19) { const lo = this.bits(16); return ((this.bits(b - 16) << 16) | lo) >>> 0; }
    const len = (this.length >>>= b), sym = Math.floor(this.value / len);
    this.value -= len * sym;
    if (this.length < 0x1000000) this.renorm();
    return sym;
  }
}
/** LASzip integer (de)compressor: corrector bit count from a per-context model, then the corrector itself. */
class IntCoder {
  constructor(dec, bits, contexts = 1) {
    this.dec = dec; this.bits = bits; this.k = 0; this.high = 8;
    this.range = bits < 32 ? 2 ** bits : 0; this.min = bits < 32 ? -(this.range / 2) : -2147483648;
    this.mBits = Array.from({ length: contexts }, () => new SymbolModel(bits + 1));
    this.mCorr = [new BitModel()];
    for (let i = 1; i <= bits; i++) this.mCorr.push(new SymbolModel(1 << Math.min(i, this.high)));
  }
  read(pred, ctx = 0) {
    const d = this.dec, k = (this.k = d.symbol(this.mBits[ctx]));
    let c;
    if (!k) c = d.bit(this.mCorr[0]);
    else if (k < 32) {
      if (k <= this.high) c = d.symbol(this.mCorr[k]);
      else { const k1 = k - this.high; c = d.symbol(this.mCorr[k]); c = c * 2 ** k1 + d.bits(k1); }
      if (c >= 2 ** (k - 1)) c += 1; else c -= 2 ** k - 1;
    } else c = this.min;
    let real = pred + c;
    if (!this.range) return real | 0;
    if (real < 0) real += this.range; else if (real >= this.range) real -= this.range;
    return real;
  }
}
class Median5 {
  constructor() { this.v = [0, 0, 0, 0, 0]; this.high = true; }
  add(x) {
    const v = this.v;
    if (this.high) {
      if (x < v[2]) { v[4] = v[3]; v[3] = v[2]; if (x < v[0]) { v[2] = v[1]; v[1] = v[0]; v[0] = x; } else if (x < v[1]) { v[2] = v[1]; v[1] = x; } else v[2] = x; }
      else { if (x < v[3]) { v[4] = v[3]; v[3] = x; } else v[4] = x; this.high = false; }
    } else if (v[2] < x) { v[0] = v[1]; v[1] = v[2]; if (v[4] < x) { v[2] = v[3]; v[3] = v[4]; v[4] = x; } else if (v[3] < x) { v[2] = v[3]; v[3] = x; } else v[2] = x; }
    else { if (v[1] < x) { v[0] = v[1]; v[1] = x; } else v[0] = x; this.high = true; }
  }
}
const RET_MAP = [[15, 14, 13, 12, 11, 10, 9, 8], [14, 0, 1, 3, 6, 10, 10, 9], [13, 1, 2, 4, 7, 11, 11, 10], [12, 3, 4, 5, 8, 12, 12, 11], [11, 6, 7, 8, 9, 13, 13, 12], [10, 10, 11, 12, 13, 14, 14, 13], [9, 10, 11, 12, 13, 14, 15, 14], [8, 9, 10, 11, 12, 13, 14, 15]];

// ---- item decoders (version 2); each keeps only the state that later symbols depend on ---------------------------------
function point10(dec, u8, p) {
  const dv = new DataView(u8.buffer, u8.byteOffset + p, 20), mChanged = new SymbolModel(64), icInt = new IntCoder(dec, 16, 4), mAngle = [new SymbolModel(256), new SymbolModel(256)], icSrc = new IntCoder(dec, 16);
  const mBit = [], mClass = [], mUser = [], icX = new IntCoder(dec, 32, 2), icY = new IntCoder(dec, 32, 22), icZ = new IntCoder(dec, 32, 20), mx = [], my = [], lastInt = new Uint16Array(16), lastH = new Int32Array(8);
  for (let i = 0; i < 16; i++) { mx.push(new Median5()); my.push(new Median5()); }
  const st = { x: dv.getInt32(0, true), y: dv.getInt32(4, true), z: dv.getInt32(8, true) };
  let bitByte = u8[p + 14], cls = u8[p + 15], user = u8[p + 17], src = dv.getUint16(18, true);
  st.read = () => {
    const ch = dec.symbol(mChanged);
    if (ch & 32) bitByte = dec.symbol(mBit[bitByte] || (mBit[bitByte] = new SymbolModel(256)));
    const r = bitByte & 7, n = (bitByte >> 3) & 7, m = RET_MAP[n][r], l = Math.abs(n - r);
    if (ch & 16) lastInt[m] = icInt.read(lastInt[m], m < 3 ? m : 3);
    if (ch & 8) cls = dec.symbol(mClass[cls] || (mClass[cls] = new SymbolModel(256)));
    if (ch & 4) dec.symbol(mAngle[(bitByte >> 6) & 1]);
    if (ch & 2) user = dec.symbol(mUser[user] || (mUser[user] = new SymbolModel(256)));
    if (ch & 1) src = icSrc.read(src);
    const n1 = n === 1 ? 1 : 0, dx = icX.read(mx[m].v[2], n1);
    st.x = (st.x + dx) | 0; mx[m].add(dx);
    const kx = icX.k, dy = icY.read(my[m].v[2], n1 + (kx < 20 ? kx & ~1 : 20));
    st.y = (st.y + dy) | 0; my[m].add(dy);
    const kb = (kx + icY.k) >> 1;
    st.z = lastH[l] = icZ.read(lastH[l], n1 + (kb < 18 ? kb & ~1 : 18));
  };
  return st;
}
function gpstime11(dec) {
  const mMulti = new SymbolModel(516), m0 = new SymbolModel(6), ic = new IntCoder(dec, 32, 9), diff = new Int32Array(4), extreme = new Int32Array(4);
  let last = 0, next = 0;
  const full = () => { next = (next + 1) & 3; ic.read(0, 8); dec.bits(32); last = next; diff[last] = 0; extreme[last] = 0; };
  const bump = (d) => { if (++extreme[last] > 3) { diff[last] = d; extreme[last] = 0; } };
  const read = (depth) => {
    if (diff[last] === 0) {
      const multi = dec.symbol(m0);
      if (multi === 1) { diff[last] = ic.read(0, 0); extreme[last] = 0; }
      else if (multi === 2) full();
      else if (multi > 2) { last = (last + multi - 2) & 3; if (depth < 8) read(depth + 1); }
    } else {
      let multi = dec.symbol(mMulti);
      if (multi === 1) { ic.read(diff[last], 1); extreme[last] = 0; }
      else if (multi < 511) {
        if (multi === 0) bump(ic.read(0, 7));
        else if (multi < 500) ic.read(Math.imul(multi, diff[last]), multi < 10 ? 2 : 3);
        else if (multi === 500) bump(ic.read(Math.imul(500, diff[last]), 4));
        else { multi = 500 - multi; if (multi > -10) ic.read(Math.imul(multi, diff[last]), 5); else bump(ic.read(Math.imul(-10, diff[last]), 6)); }
      } else if (multi === 512) full();
      else if (multi > 512) { last = (last + multi - 512) & 3; if (depth < 8) read(depth + 1); }
    }
  };
  return { read: () => read(0) };
}
function rgb12(dec) {
  const used = new SymbolModel(128), m = [0, 1, 2, 3, 4, 5].map(() => new SymbolModel(256));
  return { read() { const s = dec.symbol(used); if (s & 1) dec.symbol(m[0]); if (s & 2) dec.symbol(m[1]); if (s & 64) for (let k = 2; k < 6; k++) if (s & (1 << k)) dec.symbol(m[k]); } };
}
function bytes(dec, n) {
  const m = Array.from({ length: n }, () => new SymbolModel(256));
  return { read() { for (let k = 0; k < n; k++) dec.symbol(m[k]); } };
}

// ---- LAS 1.4 points (POINT14, item version 3): one arithmetic-coded layer per attribute ----------------------------------
// return-map context of (number of returns, return number); the return level is their difference capped at 7
const MAP6 = ['0123453445555555', '1013333333333333', '2124444444433333', '3345444444444444', '4344544444444444', '5344454444444444', '3344445444444444', '4344444544444444', '4344444454444444', '5344444445444444', '5344444444544444', '5334444444455444', '5334444444455544', '5334444444445554', '5334444444444555', '5334444444444455'];
/**
 * Decoder of the x, y, z of a layered chunk: the "channel, returns and XY" layer and the Z layer. The other layers
 * (classification, flags, intensity, scan angle, user data, point source, GPS time, colour, extra bytes) are skipped by
 * their sizes; of their content only the change flags carried by the first layer are needed. Each scanner channel keeps
 * a context of its own, started from the point at which the channel is first met.
 */
function point14(u8, p, decXY, decZ) {
  const dv = new DataView(u8.buffer, u8.byteOffset + p, 30), ctx = [null, null, null, null];
  const make = (X, Y, Z, n, r) => {
    const c = { X, Y, Z, n, r, gps: 0, mChanged: [], mChannel: new SymbolModel(3), mN: [], mR: [], mSame: new SymbolModel(13), icX: new IntCoder(decXY, 32, 2), icY: new IntCoder(decXY, 32, 22), icZ: decZ ? new IntCoder(decZ, 32, 20) : null, mx: [], my: [], lastZ: new Int32Array(8).fill(Z) };
    for (let i = 0; i < 8; i++) c.mChanged.push(new SymbolModel(128));
    for (let i = 0; i < 12; i++) { c.mx.push(new Median5()); c.my.push(new Median5()); }
    return c;
  };
  let cur = (u8[p + 15] >> 4) & 3, c = (ctx[cur] = make(dv.getInt32(0, true), dv.getInt32(4, true), dv.getInt32(8, true), u8[p + 14] >> 4, u8[p + 14] & 15));
  const st = { x: c.X, y: c.Y, z: c.Z };
  st.read = () => {
    const ch = decXY.symbol(c.mChanged[(c.r === 1 ? 1 : 0) + (c.r >= c.n ? 2 : 0) + (c.gps ? 4 : 0)]);
    if (ch & 64) { const to = (cur + decXY.symbol(c.mChannel) + 1) & 3; if (!ctx[to]) ctx[to] = make(c.X, c.Y, c.Z, c.n, c.r); cur = to; c = ctx[to]; }
    const gps = (ch >> 4) & 1;
    if (ch & 4) c.n = decXY.symbol(c.mN[c.n] || (c.mN[c.n] = new SymbolModel(16)));
    const how = ch & 3;
    if (how === 1) c.r = (c.r + 1) & 15; else if (how === 2) c.r = (c.r + 15) & 15;
    else if (how === 3) c.r = gps ? decXY.symbol(c.mR[c.r] || (c.mR[c.r] = new SymbolModel(16))) : (c.r + decXY.symbol(c.mSame) + 2) & 15;
    const n = c.n, r = c.r, mi = ((MAP6[n].charCodeAt(r) - 48) << 1) | gps, l = Math.min(7, Math.abs(n - r)), n1 = n === 1 ? 1 : 0;
    const dx = c.icX.read(c.mx[mi].v[2], n1);
    c.X = (c.X + dx) | 0; c.mx[mi].add(dx);
    const kx = c.icX.k, dy = c.icY.read(c.my[mi].v[2], n1 + (kx < 20 ? kx & ~1 : 20));
    c.Y = (c.Y + dy) | 0; c.my[mi].add(dy);
    if (c.icZ) { const kb = (kx + c.icY.k) >> 1; c.Z = c.lastZ[l] = c.icZ.read(c.lastZ[l], n1 + (kb < 18 ? kb & ~1 : 18)); }
    c.gps = gps; st.x = c.X; st.y = c.Y; st.z = c.Z;
  };
  return st;
}

/** True for a LAS file whose points are LASzip-compressed. */
export const isLAZ = (u8) => u8.length > 110 && u8[0] === 0x4c && u8[1] === 0x41 && u8[2] === 0x53 && u8[3] === 0x46 && (u8[104] & 0xc0) !== 0;

/**
 * Decode a LAZ file. opts: { maxPoints = 2e6, maxDecode = 40e6 }. Returns { xyz: Float64Array, total, kept, version,
 * pointFormat, compressor, lazVersion, chunks, warnings }.
 */
export function parseLAZ(u8, opts = {}) {
  const N = u8.length, dv = new DataView(u8.buffer, u8.byteOffset, N), maxPoints = opts.maxPoints ?? 2e6, maxDecode = opts.maxDecode ?? 40e6, warnings = [];
  if (N < 227 || !(u8[0] === 0x4c && u8[1] === 0x41 && u8[2] === 0x53 && u8[3] === 0x46)) fail('Not a LAS / LAZ file (LASF signature missing).');
  const major = u8[24], minor = u8[25], hsize = dv.getUint16(94, true), offPts = dv.getUint32(96, true), nVlr = dv.getUint32(100, true), fmt = u8[104] & 0x3f, recLen = dv.getUint16(105, true);
  let total = dv.getUint32(107, true);
  if (minor >= 4 && hsize >= 375 && N >= 255) { const n64 = Number(dv.getBigUint64(247, true)); if (n64 > 0) total = n64; }
  if (major !== 1 || minor > 4 || offPts < hsize || offPts + 8 > N) fail(`Unsupported or corrupt LAZ file (version ${major}.${minor}).`);
  const f = (o) => dv.getFloat64(o, true), sc = [f(131), f(139), f(147)], of = [f(155), f(163), f(171)];
  let vlr = -1;
  for (let p = hsize, k = 0; k < nVlr && k < 4096 && p + 54 <= offPts; k++) { const len = dv.getUint16(p + 20, true); if (dv.getUint16(p + 18, true) === 22204 && len >= 34) { vlr = p + 54; break; } p += 54 + len; }
  if (vlr < 0) fail('The LAZ file has no LASzip record describing its compression. ' + HOW);
  const compressor = dv.getUint16(vlr, true), coder = dv.getUint16(vlr + 2, true), lazVersion = `${u8[vlr + 4]}.${u8[vlr + 5]}.${dv.getUint16(vlr + 6, true)}`, chunkSize = dv.getUint32(vlr + 12, true), nItems = dv.getUint16(vlr + 32, true), items = [];
  if (vlr + 34 + 6 * nItems > offPts || nItems > 32) fail('The LASzip record is corrupt.');
  for (let k = 0; k < nItems; k++) items.push({ type: dv.getUint16(vlr + 34 + 6 * k, true), size: dv.getUint16(vlr + 36 + 6 * k, true), version: dv.getUint16(vlr + 38 + 6 * k, true) });
  const layered = compressor === 3;
  if (!layered && items.some((i) => i.type === 9)) fail(`This LAZ file holds wave-packet points (format ${fmt}), which are not read. ${HOW}`);
  if (layered ? coder !== 0 || !items.length || items[0].type !== 10 || items.some((i) => i.type < 10 || i.type > 14) : (compressor !== 1 && compressor !== 2) || coder !== 0 || !items.length || items[0].type !== 6 || items.some((i) => ![0, 6, 7, 8].includes(i.type))) fail(`This LAZ file uses an unsupported LASzip layout (compressor ${compressor}). ${HOW}`);
  if (items.some((i) => i.version !== (layered ? 3 : 2))) fail(`This LAZ file uses item codec version ${items.find((i) => i.version !== (layered ? 3 : 2)).version}${layered ? '' : ' of a very old LASzip release'}, which is not read. ${HOW}`);
  const rawLen = items.reduce((s, i) => s + i.size, 0);
  if (rawLen !== recLen || items[0].size !== (layered ? 30 : 20)) fail('The LASzip record does not match the point record length.');
  // layered chunks: point count, then one size word per layer (9 of the point item, 1 colour, 2 colour + NIR, 1 wave packet, 1 per extra byte)
  const nLayers = items.reduce((s, i) => s + (i.type === 10 ? 9 : i.type === 12 ? 2 : i.type === 14 ? i.size : 1), 0);
  if (!(total > 0)) fail('The LAZ file holds no points.');
  // chunk byte ranges
  const chunks = [];
  if (compressor === 1) chunks.push([offPts, N, total]);
  else {
    const variable = chunkSize === 0xffffffff;
    if (!(chunkSize > 0)) fail('The LASzip record is corrupt (chunk size 0).');
    let tp = Number(dv.getBigInt64(offPts, true));
    if (tp === -1 && N >= offPts + 16) tp = Number(dv.getBigInt64(N - 8, true));
    const nch = variable ? (tp + 8 <= N && tp >= 0 ? dv.getUint32(tp + 4, true) : 0) : Math.ceil(total / chunkSize);
    if (!(tp >= offPts + 8 && tp + 8 <= N) || dv.getUint32(tp, true) !== 0) fail('The LAZ chunk table is missing or corrupt; the file appears truncated.');
    const stored = dv.getUint32(tp + 4, true);
    if (stored > 1e7 || stored < nch - 1) fail('The LAZ chunk table is corrupt.');
    let start = offPts + 8;
    if (stored) {
      const d = new Decoder(u8, tp + 8, N), ic = new IntCoder(d, 32, 2);
      for (let k = 0, prev = 0, cnt = 0, done = 0; k < stored && k < nch; k++) {
        if (variable) cnt = ic.read(cnt, 0) >>> 0;        // variable chunks (as written for COPC) carry their point counts
        prev = ic.read(prev, 1) >>> 0;
        if (start + prev > tp) fail('The LAZ chunk table is corrupt.');
        const n = variable ? cnt : Math.min(chunkSize, total - k * chunkSize);
        chunks.push([start, start + prev, Math.min(n, total - done)]); start += prev; done += n;
      }
    }
    if (!variable && chunks.length < nch) chunks.push([start, tp, total - chunks.length * chunkSize]);
  }
  // even sub-sampling: every `step`-th point of every `cstep`-th chunk
  const cstep = Math.max(1, Math.ceil(total / maxDecode)), decoded = cstep === 1 ? total : chunks.reduce((s, c, k) => s + (k % cstep ? 0 : c[2]), 0), step = Math.max(1, Math.ceil(decoded / maxPoints));
  const xyz = new Float64Array(3 * (Math.ceil(decoded / step) + chunks.length));
  let kept = 0, seen = 0, truncated = false;
  const keep = (st) => { if (seen++ % step === 0) { xyz[3 * kept] = st.x * sc[0] + of[0]; xyz[3 * kept + 1] = st.y * sc[1] + of[1]; xyz[3 * kept + 2] = st.z * sc[2] + of[2]; kept++; } };
  for (let c = 0; c < chunks.length; c += cstep) {
    const [a, b, n] = chunks[c];
    if (!(n > 0) || a + rawLen + 4 > b || b > N) { truncated = true; break; }
    let dec, pt;
    const rest = [];
    if (layered) {
      const sizes = a + rawLen + 4, data = sizes + 4 * nLayers;
      if (data > b) { truncated = true; break; }
      const nxy = dv.getUint32(sizes, true), nz = dv.getUint32(sizes + 4, true);
      if (data + nxy + nz > b || (n > 1 && nxy < 4)) { truncated = true; break; }
      dec = new Decoder(u8, data, data + nxy);
      pt = point14(u8, a, dec, nz ? new Decoder(u8, data + nxy, data + nxy + nz) : null);
    } else {
      dec = new Decoder(u8, a + rawLen, b); pt = point10(dec, u8, a);
      for (let k = 1; k < items.length; k++) rest.push(items[k].type === 7 ? gpstime11(dec) : items[k].type === 8 ? rgb12(dec) : bytes(dec, items[k].size));
    }
    keep(pt);
    for (let i = 1; i < n; i++) {
      pt.read();
      for (let k = 0; k < rest.length; k++) rest[k].read();
      if (dec.over > 8) { truncated = true; break; }
      keep(pt);
    }
    if (truncated) break;
  }
  if (truncated) warnings.push(`The LAZ data ends early; ${seen} of ${decoded} points were decoded.`);
  if (cstep > 1) warnings.push(`Very large cloud: only every ${cstep}th chunk of ${chunkSize} points was decompressed.`);
  if (step > 1 || cstep > 1) warnings.push(`Point cloud of ${total} points sub-sampled to ${kept}.`);
  return { xyz: xyz.subarray(0, 3 * kept), total, kept, version: `${major}.${minor}`, pointFormat: fmt, compressor: compressor === 1 ? 'pointwise' : layered ? 'layered chunked' : 'pointwise chunked', lazVersion, chunks: chunks.length, warnings };
}
