// JT (ISO 14306) container reader: file header, table of contents, segment headers and the Parasolid B-rep segments.
//
// Read:   the 80-byte version string, the byte order, the TOC (segment id, offset, length, type) of JT 8, 9 and 10 files and
//         the header of every segment (its id, type and length are checked against the TOC); jtSummary() describes what
//         a file holds. XT B-rep segments (type 17) of JT 8 and 9 are unpacked (zlib) and their Parasolid neutral binary
//         transmit streams handed out by jtBrepStreams() for the Parasolid reader.
// Not read: the tessellation. Shape LOD segments keep their triangle strips in the JT codecs (bit-length, Huffman and
//         arithmetic coded index streams with predictors, quantised or topologically compressed vertex data, Deering
//         normals); JT B-rep segments, the logical scene graph (part names, placements), PMI and the LZMA-packed segments
//         of JT 10 are not decoded either. parseJT() ends in an error that says what the file contains and how to
//         export it.
// File content is untrusted: offsets and counts are bounds-checked and the TOC size is capped.

import { inflateSync } from './hdf5.js';

function fail(msg) { const e = new Error(msg); e.user = true; throw e; }
const SEG = { 1: 'logical scene graph', 2: 'JT B-rep', 3: 'PMI data', 4: 'meta data', 6: 'shape', 7: 'shape LOD 0', 8: 'shape LOD 1', 9: 'shape LOD 2', 10: 'shape LOD 3', 11: 'shape LOD 4', 12: 'shape LOD 5', 13: 'shape LOD 6', 14: 'shape LOD 7', 15: 'shape LOD 8', 16: 'shape LOD 9', 17: 'XT B-rep', 18: 'wireframe', 20: 'ULP', 24: 'LWPA' };
const hex = (u8, p) => Array.from(u8.subarray(p, p + 16), (b) => b.toString(16).padStart(2, '0')).join('');

export const isJT = (u8) => u8.length > 100 && /^Version \d+\.\d+/.test(String.fromCharCode(...u8.subarray(0, 16)));

/**
 * Describe a JT file: { version, major, minor, littleEndian, segments: [{ id, type, typeName, offset, length, compressed? }],
 * counts: { typeName: n }, shapes (number of shape LOD-0 segments), lods, brep: 'JT B-rep' | 'XT B-rep' | null }.
 */
export function jtSummary(u8) {
  if (!isJT(u8)) fail('Not a JT file (the "Version n.n JT" header is missing).');
  const N = u8.length, dv = new DataView(u8.buffer, u8.byteOffset, N), m = /^Version (\d+)\.(\d+)/.exec(String.fromCharCode(...u8.subarray(0, 80))), major = +m[1], minor = +m[2], le = u8[80] === 0, wide = major >= 10;
  if (major < 7 || major > 12) fail(`JT version ${major}.${minor} is not known to this reader.`);
  const u32 = (p) => dv.getUint32(p, le), toc = wide ? u32(le ? 85 : 89) + u32(le ? 89 : 85) * 4294967296 : u32(85), esz = wide ? 32 : 28;
  if (!(toc >= 85 && toc + 4 <= N)) fail('The JT file is truncated or corrupt (table of contents).');
  const n = dv.getInt32(toc, le), segments = [], counts = {};
  if (!(n >= 0 && n <= 2e6) || toc + 4 + n * esz > N) fail('The JT file is truncated or corrupt (table of contents).');
  for (let k = 0; k < n; k++) {
    const e = toc + 4 + k * esz, offset = wide ? u32(e + (le ? 16 : 20)) + u32(e + (le ? 20 : 16)) * 4294967296 : u32(e + 16), length = u32(e + esz - 8), type = u32(e + esz - 4) >>> 24, typeName = SEG[type] || `type ${type}`, s = { id: hex(u8, e), type, typeName, offset, length };
    // segment header: id, type, length, then for the compressible types a zlib flag
    s.valid = offset + 24 <= N && hex(u8, offset) === s.id && u32(offset + 16) === type && u32(offset + 20) === length;
    if (s.valid && offset + 33 <= N && [1, 2, 3, 4, 17, 18, 20, 24].includes(type)) s.compressed = u32(offset + 24) === 2 && major < 10 ? true : major >= 10 ? u32(offset + 24) >= 2 : false;
    counts[typeName] = (counts[typeName] || 0) + 1;
    segments.push(s);
  }
  return { version: `${major}.${minor}`, major, minor, littleEndian: le, segments, counts, shapes: counts['shape LOD 0'] || counts.shape || 0, lods: Object.keys(counts).filter((k) => /^shape LOD/.test(k)).length, brep: counts['XT B-rep'] ? 'XT B-rep' : counts['JT B-rep'] ? 'JT B-rep' : null };
}

/**
 * Parasolid streams of the XT B-rep segments of a JT file: { summary, streams: [Uint8Array …] (neutral binary transmit data,
 * starting at its "PS" flag), lzma: number of segments left out because JT 10 packs them with LZMA }.
 * A segment holds, behind its 24-byte header, a compression flag (2: zlib, 3: LZMA), the packed length, the algorithm byte
 * and the data; unpacked, the B-rep element carries the Parasolid stream after a short element header.
 */
export function jtBrepStreams(u8, maxBytes = 256e6) {
  const s = jtSummary(u8), dv = new DataView(u8.buffer, u8.byteOffset, u8.length), streams = [];
  let lzma = 0, total = 0;
  for (const g of s.segments) {
    if (g.type !== 17 || !g.valid || g.length < 40 || g.offset + g.length > u8.length) continue;
    const flag = dv.getUint32(g.offset + 24, s.littleEndian), alg = u8[g.offset + 32];
    let data = null;
    if (flag === 2 && alg === 2) { try { data = inflateSync(u8.subarray(g.offset + 33, g.offset + g.length), maxBytes); } catch (e) { if (e && e.user) continue; throw e; } }
    else if (flag === 3 || alg === 3) { lzma++; continue; }
    else data = u8.subarray(g.offset + 24, g.offset + g.length);
    // the transmit data begin at the "PS\0\0" flag followed by the length of the "… TRANSMIT FILE" line
    let at = -1;
    for (let k = 0; k + 12 < data.length && k < 4096; k++) if (data[k] === 0x50 && data[k + 1] === 0x53 && data[k + 2] === 0 && data[k + 3] === 0 && data[k + 6] === 0x3a && data[k + 7] === 0x20) { at = k; break; }
    if (at < 0 || (total += data.length - at) > maxBytes) continue;
    streams.push(data.subarray(at));
  }
  return { summary: s, streams, lzma };
}

/** Always fails: says what the JT file holds and why its geometry is not read. */
export function parseJT(u8, summary) {
  const s = summary || jtSummary(u8), parts = Object.entries(s.counts).map(([k, v]) => `${v} × ${k}`).join(', ');
  return fail(`This JT ${s.version} file holds ${s.segments.length} segments (${parts || 'none'}). Its tessellation is stored in the JT ${s.major >= 9 ? 'topological mesh' : 'vertex-shape'} codecs (bit-length, Huffman and arithmetic coded streams, quantised vertices)${s.brep ? ` and its exact shape as ${s.brep}` : ''}, which are not decoded. Export the model as STEP AP242 or STL from NX, Teamcenter Visualization or the authoring CAD package.`);
}
