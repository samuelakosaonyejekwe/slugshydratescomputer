// ADF (Advanced Data Format) container reader: the tree-of-nodes file that CGNS used before HDF5 and still writes on request.
//
// Read:   the file header (database version "A…" with ASCII-hex disk pointers and dimensions, version "B…" with binary
//         ones; the number format IEEE big- or little-endian), node headers (name, label, sub-node count and table, data
//         type, dimensions, data chunks), sub-node tables, data chunks and data-chunk tables, the data types I4, I8, U4,
//         U8, R4, R8, C1 and B1, and links (type LK) to a node of the same file. openADF() presents the tree through the
//         calls the CGNS reader makes on an HDF5 file (root, kids, attrs, resolve, info, read), so one CGNS reader serves
//         both containers; arrays come in file order, which is the order of the HDF5 flavour.
// Not read: links into other files (counted by the CGNS reader), the Cray number formats, complex and compound data
//         types (such nodes read as empty), the free-chunk tables.
// File content is untrusted: every pointer, tag and count is bounds-checked, node and byte totals are capped and the
// tree walk is guarded against cycles.

function fail(msg) { const e = new Error(msg); e.user = true; throw e; }
const BLOCK = 4096, DATA = 1e9;                             // disk block size; address offset of a node's data
const TYPES = { I4: [4, Int32Array, 'getInt32'], U4: [4, Uint32Array, 'getUint32'], R4: [4, Float32Array, 'getFloat32'], R8: [8, Float64Array, 'getFloat64'], I8: [8, Float64Array, 'getBigInt64'], U8: [8, Float64Array, 'getBigUint64'], C1: [1, Uint8Array], B1: [1, Uint8Array] };

/** True for an ADF file ("@(#)ADF Database Version …"). */
export const isADF = (u8) => u8.length > 266 && /^ADF Database Version [A-Z]\d/.test(String.fromCharCode(...u8.subarray(4, 32)));

/**
 * Open an ADF file held in memory. opts: { maxBytes (512e6: cap on one array), maxNodes (2e5) }. Returns
 *   version, format ('IEEE_BIG_32' …), root (0), nodes (count read so far)
 *   kids(a)        Map name → { addr }; a node with data also lists ' data' → { addr: a + 1e9 }
 *   attrs(a)       { name, label, type } (type: the ADF data type, 'LK' for a link)
 *   resolve(p, a)  address of the child path p below a (' data' and ' link' as in CGNS / HDF5) or -1
 *   info(a)        { kind: 'group' | 'dataset', shape, type }
 *   read(a)        Promise<{ data, shape, type }> for a data address: typed array (64-bit integers as doubles), C1 as bytes
 */
export function openADF(u8, opts = {}) {
  if (!isADF(u8)) fail('Not an ADF file (the "ADF Database Version" header is missing).');
  const N = u8.length, dv = new DataView(u8.buffer, u8.byteOffset, N), maxBytes = opts.maxBytes ?? 512e6, maxNodes = opts.maxNodes ?? 2e5;
  const txt = (p, n) => { let s = ''; for (let k = 0; k < n && p + k < N; k++) s += String.fromCharCode(u8[p + k]); return s; }, bad = (what) => fail(`The ADF file is truncated or corrupt (${what}).`);
  const version = txt(25, 6), old = u8[25] !== 0x42, fm = txt(100, 2), le = fm[0] === 'L';
  if (!/^[BL][LB]$/.test(fm)) fail(`This ADF file stores its numbers in the format "${fm}" (a Cray or machine-native layout), which is not read; IEEE big- and little-endian files are. Convert it with the CGNS tools: cgnsconvert -h in.cgns out.cgns.`);
  const hex = (p, n) => { const v = parseInt(txt(p, n), 16); if (!(v >= 0)) bad('number'); return v; };
  /** Disk pointer (block, offset) at p -> byte address, or -1 for the null pointer. */
  const ptr = (p) => { if (p + 12 > N) bad('pointer'); const b = old ? hex(p, 8) : Number(dv.getBigUint64(p, le)), o = old ? hex(p + 8, 4) : dv.getUint32(p + 8, le); if (o > BLOCK || !(b >= 0)) bad('pointer'); return o === BLOCK && b === 0 ? -1 : b * BLOCK + o; };
  const tag = (p, t) => p >= 0 && p + 4 <= N && txt(p, 4) === t;
  const nodes = [], byPos = new Map();
  /** Node header at byte address p -> node index. */
  const node = (p) => {
    let i = byPos.get(p);
    if (i !== undefined) return i;
    if (!tag(p, 'NoDe') || p + 246 > N || !tag(p + 242, 'TaiL')) bad('node header');
    if (nodes.length >= maxNodes) fail('The ADF file holds too many nodes.');
    const nd = { name: txt(p + 4, 32).replace(/[\s\0]+$/, ''), label: txt(p + 36, 32).replace(/[\s\0]+$/, ''), nsub: hex(p + 68, 8), table: ptr(p + 84), type: txt(p + 96, 32).replace(/[\s\0]+$/, ''), ndim: hex(p + 128, 2), dims: [], nchunk: hex(p + 226, 4), chunk: ptr(p + 230), kids: null };
    if (nd.ndim > 12) bad('dimensions');
    for (let k = 0; k < nd.ndim; k++) nd.dims.push(old ? hex(p + 130 + 8 * k, 8) : Number(dv.getBigUint64(p + 130 + 8 * k, le)));
    i = nodes.length; nodes.push(nd); byPos.set(p, i);
    return i;
  };
  const kids = (a) => {
    const nd = nodes[a];
    if (!nd) return new Map();
    if (!nd.kids) {
      nd.kids = new Map();
      if (nd.nsub > 0 && nd.table >= 0) {
        const t = nd.table;
        if (!tag(t, 'SNTb') || nd.nsub > 1e6 || t + 16 + 44 * nd.nsub > N) bad('sub-node table');
        for (let k = 0; k < nd.nsub; k++) { const e = t + 16 + 44 * k, name = txt(e, 32).replace(/[\s\0]+$/, ''), c = ptr(e + 32); if (c >= 0 && name) nd.kids.set(name, { addr: node(c) }); }
      }
      if (TYPES[nd.type] && nd.ndim > 0) nd.kids.set(' data', { addr: a + DATA });
    }
    return nd.kids;
  };
  /** Byte ranges of the data of a node, in order. */
  const ranges = (nd) => {
    const out = [], one = (p, end) => { if (!tag(p, 'DaTa') || !(end >= p + 16 && end + 4 <= N)) bad('data chunk'); out.push([p + 16, end]); };
    if (nd.nchunk === 1 && nd.chunk >= 0) one(nd.chunk, ptr(nd.chunk + 4));
    else if (nd.nchunk > 1 && nd.chunk >= 0) {
      const t = nd.chunk;
      if (!tag(t, 'DCtb') || t + 16 + 24 * nd.nchunk > N) bad('data-chunk table');
      for (let k = 0; k < nd.nchunk; k++) { const s = ptr(t + 16 + 24 * k), e = ptr(t + 28 + 24 * k); if (s < 0) bad('data chunk'); one(s, Math.min(e, ptr(s + 4))); }
    }
    return out;
  };
  const bytesOf = (nd, want) => {
    if (want > maxBytes) fail('An array of the ADF file is too large.');
    const out = new Uint8Array(want);
    let o = 0;
    for (const [s, e] of ranges(nd)) { const n = Math.min(e - s, want - o); if (n > 0) { out.set(u8.subarray(s, s + n), o); o += n; } }
    if (o < want) bad('array data');
    return out;
  };
  const link = (a, depth = 0) => {
    // link data: "file\0path" as C1; an empty file name points into this file
    const nd = nodes[a];
    if (!nd || nd.type !== 'LK' || depth > 20) return -1;
    const s = txt(0, 0) + String.fromCharCode(...bytesOf(nd, Math.min(4096, nd.dims[0] || 0))), cut = s.indexOf('\0'), file = cut >= 0 ? s.slice(0, cut) : '', path = (cut >= 0 ? s.slice(cut + 1) : s).replace(/\0[\s\S]*$/, '');
    if (file.trim()) return -1;
    let at = 0;
    for (const part of path.split('/').filter(Boolean)) { const k = kids(at).get(part); if (!k) return -1; at = nodes[k.addr].type === 'LK' ? link(k.addr, depth + 1) : k.addr; if (at < 0) return -1; }
    return at;
  };
  const root = node(ptr(134));
  return {
    version, format: `IEEE_${le ? 'LITTLE' : 'BIG'}_${fm[1] === 'B' ? 64 : 32}`, root, get nodes() { return nodes.length; }, kids,
    attrs: (a) => { const nd = nodes[a]; return nd ? { name: nd.name, label: nd.label, type: nd.type } : {}; },
    resolve: (path, from = root) => { if (path === ' link') return link(from); if (path === ' data') return nodes[from] && TYPES[nodes[from].type] && nodes[from].ndim > 0 ? from + DATA : -1; let at = from; for (const part of String(path).split('/').filter(Boolean)) { const k = kids(at).get(part); if (!k) return -1; at = k.addr; } return at; },
    info: (a) => { const nd = nodes[a >= DATA ? a - DATA : a]; return nd ? { kind: a >= DATA ? 'dataset' : 'group', shape: nd.dims.slice().reverse(), type: nd.type } : { kind: 'none' }; },
    read: async (a) => {
      const nd = nodes[a - DATA], ty = nd && TYPES[nd.type];
      if (!ty) fail('The ADF node holds no readable array.');
      const count = nd.dims.reduce((x, y) => x * y, 1);
      if (!(count >= 0) || count * ty[0] > maxBytes) fail('An array of the ADF file is too large.');
      const raw = bytesOf(nd, count * ty[0]);
      if (ty[0] === 1) return { data: raw, shape: nd.dims.slice().reverse(), type: nd.type };
      const out = new ty[1](count), rv = new DataView(raw.buffer), big = ty[2].includes('Big');
      for (let k = 0; k < count; k++) out[k] = big ? Number(rv[ty[2]](k * ty[0], le)) : rv[ty[2]](k * ty[0], le);
      return { data: out, shape: nd.dims.slice().reverse(), type: nd.type };
    },
  };
}
