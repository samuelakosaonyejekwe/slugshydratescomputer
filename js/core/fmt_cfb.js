// Compound File Binary (OLE2 structured storage, [MS-CFB]) reader: the container of DGN V8 and of many native CAD files.
//
// Read:   the header of version 3 (512-byte sectors) and version 4 (4096-byte sectors) files, the DIFAT chain, the FAT, the
//         mini-FAT with its mini stream, and the directory (storages and streams with their full paths).
// Not read: property-set contents, transactions; the red-black ordering of the directory is not checked.
// File content is untrusted: sector chains are bounds-checked and guarded against loops, entry counts and sizes are capped.

function fail(msg) { const e = new Error(msg); e.user = true; throw e; }

/** True for a compound file (signature D0 CF 11 E0 A1 B1 1A E1). */
export const isCFB = (u8) => u8.length >= 512 && u8[0] === 0xd0 && u8[1] === 0xcf && u8[2] === 0x11 && u8[3] === 0xe0 && u8[4] === 0xa1 && u8[5] === 0xb1 && u8[6] === 0x1a && u8[7] === 0xe1;

/**
 * Open a compound file. Returns { version, sectorSize, entries: [{ path, name, type: 'storage' | 'stream' | 'root', size, clsid }],
 * read(pathOrEntry, maxBytes = 256e6) → Uint8Array, find(test) → entries whose path matches a RegExp or predicate }.
 */
export function openCFB(u8) {
  if (!isCFB(u8)) fail('Not a compound file (the OLE2 signature is missing).');
  const N = u8.length, dv = new DataView(u8.buffer, u8.byteOffset, N), bad = (w) => fail(`The compound file is truncated or corrupt (${w}).`);
  const major = dv.getUint16(26, true), ss = 1 << dv.getUint16(30, true), mss = 1 << dv.getUint16(32, true);
  if ((major !== 3 && major !== 4) || (ss !== 512 && ss !== 4096) || mss !== 64) bad('header');
  const nFat = dv.getUint32(44, true), dirStart = dv.getUint32(48, true), cutoff = dv.getUint32(56, true), miniStart = dv.getUint32(60, true), difStart = dv.getUint32(68, true), nSec = Math.floor((N - ss) / ss) + 1;
  const off = (s) => { const o = (s + 1) * ss; if (!(s >= 0 && o + ss <= N + ss - 1 && o < N)) bad('sector'); return o; };
  // DIFAT: 109 entries in the header, then chained sectors
  const fatSecs = [];
  for (let k = 0; k < 109 && fatSecs.length < nFat; k++) { const s = dv.getUint32(76 + 4 * k, true); if (s < 0xfffffffa) fatSecs.push(s); }
  for (let s = difStart, g = 0; s < 0xfffffffa && fatSecs.length < nFat && g < 1e5; g++) { const o = off(s); for (let k = 0; k < ss / 4 - 1 && fatSecs.length < nFat; k++) { const v = dv.getUint32(o + 4 * k, true); if (v < 0xfffffffa) fatSecs.push(v); } s = o + ss <= N ? dv.getUint32(o + ss - 4, true) : 0xfffffffe; }
  if (fatSecs.length > nSec) bad('FAT');
  const fat = new Uint32Array(fatSecs.length * (ss / 4)).fill(0xffffffff);
  fatSecs.forEach((s, i) => { const o = off(s); for (let k = 0; k < ss / 4 && o + 4 * k + 4 <= N; k++) fat[i * (ss / 4) + k] = dv.getUint32(o + 4 * k, true); });
  /** Bytes of a sector chain (size < 0: to the end of the chain). */
  const chain = (start, size, max) => {
    const out = [], seen = new Set();
    let total = 0;
    for (let s = start; s < 0xfffffffa && (size < 0 || total < size); s = fat[s]) {
      if (s >= fat.length || seen.has(s) || total > max) bad('sector chain');
      seen.add(s);
      const o = off(s), n = Math.min(ss, N - o);
      out.push(u8.subarray(o, o + n)); total += n;
    }
    const buf = new Uint8Array(size < 0 ? total : size);
    let p = 0;
    for (const b of out) { const n = Math.min(b.length, buf.length - p); if (n <= 0) break; buf.set(b.subarray(0, n), p); p += n; }
    if (size >= 0 && p < size) bad('stream data');
    return buf;
  };
  const dir = chain(dirStart, -1, 64e6), dd = new DataView(dir.buffer), n = Math.min(Math.floor(dir.length / 128), 2e5), raw = [];
  for (let i = 0; i < n; i++) {
    const p = i * 128, len = Math.min(dd.getUint16(p + 64, true), 64), t = dir[p + 66];
    let name = '';
    for (let k = 0; k + 2 < len; k += 2) name += String.fromCharCode(dd.getUint16(p + k, true));
    raw.push({ name, t, left: dd.getUint32(p + 68, true), right: dd.getUint32(p + 72, true), child: dd.getUint32(p + 76, true), start: dd.getUint32(p + 116, true), size: dd.getUint32(p + 120, true) + (major === 4 ? dd.getUint32(p + 124, true) * 4294967296 : 0), clsid: Array.from(dir.subarray(p + 80, p + 96), (b) => b.toString(16).padStart(2, '0')).join('') });
  }
  if (!raw.length || raw[0].t !== 5) bad('root entry');
  const entries = [], seen = new Set();
  const walk = (i, parent, depth) => {
    for (const stack = [i]; stack.length; ) {
      const k = stack.pop();
      if (k >= raw.length || seen.has(k) || depth > 64) continue;
      seen.add(k);
      const r = raw[k], path = parent + r.name;
      if (r.t === 1 || r.t === 2) entries.push({ path, name: r.name, type: r.t === 1 ? 'storage' : 'stream', size: r.t === 2 ? r.size : 0, clsid: r.clsid, start: r.start });
      stack.push(r.left, r.right);
      if (r.t === 1 && r.child < raw.length) walk(r.child, path + '/', depth + 1);
    }
  };
  entries.push({ path: '', name: raw[0].name, type: 'root', size: 0, clsid: raw[0].clsid, start: raw[0].start });
  walk(raw[0].child, '', 0);
  let mini = null, miniFat = null;
  const read = (what, max = 256e6) => {
    const e = typeof what === 'string' ? entries.find((x) => x.path === what && x.type === 'stream') : what;
    if (!e || e.type !== 'stream') fail('The compound file holds no such stream.');
    if (e.size > max) fail('A stream of the compound file is too large.');
    if (e.size >= cutoff) return chain(e.start, e.size, max);
    if (!mini) { mini = chain(raw[0].start, raw[0].size, 256e6); const mf = chain(miniStart, -1, 64e6); miniFat = new Uint32Array(mf.buffer, 0, mf.length >> 2); if (new Uint8Array(new Uint16Array([1]).buffer)[0] !== 1) for (let k = 0; k < miniFat.length; k++) miniFat[k] = new DataView(mf.buffer).getUint32(4 * k, true); }
    const out = new Uint8Array(e.size), seenM = new Set();
    let p = 0;
    for (let s = e.start; p < e.size; s = miniFat[s]) { if (!(s < miniFat.length) || seenM.has(s) || (s + 1) * 64 > mini.length + 63) bad('mini stream'); seenM.add(s); const nb = Math.min(64, e.size - p); out.set(mini.subarray(s * 64, s * 64 + nb), p); p += nb; }
    return out;
  };
  return { version: major, sectorSize: ss, entries, read, find: (test) => entries.filter((e) => (test instanceof RegExp ? test.test(e.path) : test(e))) };
}
