// ASTM E57 (E2807) 3-D imaging file reader: point coordinates of every scan.
//
// Read:   the 48-byte file header, the paged layout (payload of each 1024-byte page without its trailing CRC-32C, which
//         is not verified), the XML section (data3D scans with name, pose and the points CompressedVector prototype),
//         CompressedVector binary sections (data packets; index and empty packets skipped) with the bit-pack codec for
//         Float (single / double), Integer and ScaledInteger fields, bit streams running on across packets;
//         cartesianX / Y / Z, or sphericalRange / Azimuth / Elevation converted to Cartesian; points flagged by
//         cartesianInvalidState / sphericalInvalidState are dropped; each scan's pose (rotation quaternion + translation)
//         is applied, so all scans share the file's global frame; clouds above the cap are evenly sub-sampled.
// Not read: images2D (counted), intensity / colour / row-column and time fields (their byte streams are skipped), index
//         packets, scans without Cartesian or spherical coordinates, integer fields wider than 53 bits.
// File content is untrusted: offsets, packet lengths and record counts are bounds-checked; XML size and nesting are capped.

function fail(msg) { const e = new Error(msg); e.user = true; throw e; }
const utf8 = new TextDecoder();

export const isE57 = (u8) => u8.length >= 48 && utf8.decode(u8.subarray(0, 8)) === 'ASTM-E57';

/** Tiny XML scanner → { name, attrs, kids, text }; no entities beyond the five predefined, no DTD. */
function xmlTree(s) {
  const root = { name: '#', attrs: {}, kids: [], text: '' }, stack = [root], re = /<!--[\s\S]*?-->|<\?[\s\S]*?\?>|<!\[CDATA\[([\s\S]*?)\]\]>|<(\/?)([^\s>\/]+)([^>]*?)(\/?)>|([^<]+)/g, ent = (t) => t.replace(/&(lt|gt|amp|quot|apos);/g, (_, c) => ({ lt: '<', gt: '>', amp: '&', quot: '"', apos: "'" })[c]);
  let m, count = 0;
  while ((m = re.exec(s))) {
    const top = stack[stack.length - 1];
    if (m[1] !== undefined) top.text += m[1];
    else if (m[6] !== undefined) { if (stack.length > 1) top.text += ent(m[6]); }
    else if (m[3]) {
      if (m[2]) { if (stack.length > 1) stack.pop(); continue; }
      if (++count > 2e6 || stack.length > 64) fail('The E57 XML section is too large or nested too deeply.');
      const node = { name: m[3].replace(/^.*:/, ''), attrs: {}, kids: [], text: '' }, ar = /([^\s=]+)\s*=\s*(?:"([^"]*)"|'([^']*)')/g;
      let a;
      while ((a = ar.exec(m[4]))) node.attrs[a[1]] = ent(a[2] ?? a[3]);
      top.kids.push(node);
      if (!m[5]) stack.push(node);
    }
  }
  return root;
}
const kid = (n, name) => (n ? n.kids.find((k) => k.name === name) : undefined);
const num = (n, name, dflt) => { const k = kid(n, name), v = k ? parseFloat(k.text) : NaN; return Number.isFinite(v) ? v : dflt; };

/**
 * Parse an E57 file. opts: { maxPoints = 2e6 }. Returns { xyz: Float64Array (x, y, z of the kept points), total (records in
 * the file), kept, invalid (records dropped by their invalid-state flag), scans: [{ name, points, fields, posed }],
 * images, version, warnings }.
 */
export function parseE57(u8, opts = {}) {
  if (!isE57(u8)) fail('Not an E57 file (signature "ASTM-E57" missing).');
  const dv = new DataView(u8.buffer, u8.byteOffset, u8.byteLength), N = u8.length, maxPoints = opts.maxPoints ?? 2e6, warnings = [];
  const u64 = (p) => dv.getUint32(p, true) + dv.getUint32(p + 4, true) * 4294967296;
  const major = dv.getUint32(8, true), minor = dv.getUint32(12, true), xmlPhys = u64(24), xmlLen = u64(32), ps = u64(40), body = ps - 4;
  if (major !== 1 || !(ps >= 64 && ps <= 1 << 20) || xmlPhys >= N || !(xmlLen > 0 && xmlLen <= 64e6)) fail('The E57 header is corrupt or of an unsupported version.');
  const toLog = (p) => p - 4 * Math.floor(p / ps), toPhys = (l) => l + 4 * Math.floor(l / body);
  /** n logical bytes from logical offset l (page checksums left out). */
  const read = (l, n) => {
    const out = new Uint8Array(n);
    for (let o = 0; o < n;) { const p = toPhys(l + o), room = body - (p % ps), take = Math.min(room, n - o); if (p + take > N) fail('The E57 file is truncated.'); out.set(u8.subarray(p, p + take), o); o += take; }
    return out;
  };
  const xml = xmlTree(utf8.decode(read(toLog(xmlPhys), xmlLen))), rootEl = kid(xml, 'e57Root');
  if (!rootEl) fail('The E57 XML section holds no e57Root element.');
  const scansEl = (kid(rootEl, 'data3D') || { kids: [] }).kids.slice(0, 10000), images = (kid(rootEl, 'images2D') || { kids: [] }).kids.length;
  if (!scansEl.length) fail(`The E57 file holds no 3-D scans${images ? ` (only ${images} images)` : ''}.`);
  const total = scansEl.reduce((s, sc) => s + (+((kid(sc, 'points') || { attrs: {} }).attrs.recordCount) || 0), 0), step = Math.max(1, Math.ceil(total / maxPoints));
  if (!(total > 0)) fail('The E57 scans hold no points.');
  if (total > 4e9) fail('The E57 file declares an impossible number of points.');
  const xyz = new Float64Array(3 * Math.ceil(total / step + scansEl.length)), scans = [];
  let kept = 0, invalid = 0, gidx = 0;
  for (const sc of scansEl) {
    const pts = kid(sc, 'points'), proto = kid(pts, 'prototype'), n = pts ? +pts.attrs.recordCount || 0 : 0, name = (kid(sc, 'name') || { text: '' }).text.trim();
    if (!pts || !proto || !n) { scans.push({ name, points: 0, fields: [], posed: false }); continue; }
    if (proto.kids.length > 256) fail('An E57 scan declares too many fields.');
    const fields = proto.kids.map((f) => {
      const t = f.attrs.type, a = f.attrs;
      if (t === 'Float') return { name: f.name, bits: a.precision === 'single' ? 32 : 64, float: true };
      const min = a.minimum !== undefined ? +a.minimum : -9223372036854775808, max = a.maximum !== undefined ? +a.maximum : 9223372036854775807, range = max - min;
      if (t !== 'Integer' && t !== 'ScaledInteger') fail(`E57 field "${f.name}" has the unsupported type "${String(t).slice(0, 30)}".`);
      const bits = range <= 0 ? 0 : Math.ceil(Math.log2(range + 1));
      return { name: f.name, bits, min, scale: t === 'ScaledInteger' && a.scale !== undefined ? +a.scale : 1, offset: t === 'ScaledInteger' && a.offset !== undefined ? +a.offset : 0 };
    });
    const idx = (nm) => fields.findIndex((f) => f.name === nm), cart = ['cartesianX', 'cartesianY', 'cartesianZ'].map(idx), sph = ['sphericalRange', 'sphericalAzimuth', 'sphericalElevation'].map(idx);
    const useCart = cart[0] >= 0 && cart[1] >= 0, useSph = !useCart && sph[0] >= 0 && sph[1] >= 0, want = useCart ? cart : sph, inv = idx(useCart ? 'cartesianInvalidState' : 'sphericalInvalidState');
    if (!useCart && !useSph) { warnings.push(`Scan "${name || scans.length + 1}" has neither Cartesian nor spherical coordinates and was skipped.`); scans.push({ name, points: 0, fields: fields.map((f) => f.name), posed: false }); gidx += n; continue; }
    for (const k of [...want, inv]) if (k >= 0 && !fields[k].float && fields[k].bits > 53) fail('E57 integer fields wider than 53 bits are not read.');
    const m = Math.ceil(n / step) + 1, cols = [0, 1, 2, 3].map((c) => ((c < 3 ? want[c] : inv) >= 0 ? new Float64Array(m) : null)), need = new Map();
    [...want, inv].forEach((k, c) => { if (k >= 0) need.set(k, { f: fields[k], col: cols[c], rec: 0, out: 0, buf: new Uint8Array(0), bit: 0 }); });
    // records kept: those whose global index is a multiple of `step`
    const first = (step - (gidx % step)) % step;
    const feed = (st, bytes) => {
      let b = bytes;
      if (st.buf.length) { b = new Uint8Array(st.buf.length + bytes.length); b.set(st.buf); b.set(bytes, st.buf.length); }
      const { f, col } = st, bd = new DataView(b.buffer, b.byteOffset, b.byteLength), nbits = b.length * 8;
      let bp = st.bit;
      if (!f.bits) { st.buf = new Uint8Array(0); return; }
      while (st.rec < n && bp + f.bits <= nbits) {
        if ((st.rec - first) % step === 0 && st.rec >= first) {
          let v;
          if (f.float) v = f.bits === 32 ? bd.getFloat32(bp >> 3, true) : bd.getFloat64(bp >> 3, true);
          else { v = 0; for (let got = 0, q = bp; got < f.bits;) { const off = q & 7, take = Math.min(8 - off, f.bits - got); v += ((b[q >> 3] >> off) & ((1 << take) - 1)) * 2 ** got; got += take; q += take; } v = (v + f.min) * f.scale + f.offset; }
          col[st.out++] = v;
        }
        st.rec++; bp += f.bits;
      }
      st.buf = b.subarray(bp >> 3); st.bit = bp & 7;
    };
    const sec = +pts.attrs.fileOffset;
    if (!(sec >= 48 && sec + 32 <= N) || u8[sec] !== 1) fail('An E57 point section is missing or corrupt.');
    const secLog = toLog(sec), secEnd = secLog + u64(sec + 8);
    let l = toLog(u64(sec + 16)), guard = 0;
    const states = [...need.values()], done = () => { for (const st of states) if (st.rec < n && st.f.bits) return false; return true; };
    while (l + 4 <= secEnd && !done()) {
      if (++guard > 5e7) fail('The E57 point section is corrupt.');
      const hp = toPhys(l), h = hp % ps <= body - 4 && hp + 4 <= N ? u8.subarray(hp, hp + 4) : read(l, 4), len = (h[2] | (h[3] << 8)) + 1;
      if (h[0] === 1) {
        const pk = read(l, len), nb = pk[4] | (pk[5] << 8);
        if (6 + 2 * nb > len) fail('An E57 data packet is corrupt.');
        let o = 6 + 2 * nb;
        for (let k = 0; k < nb; k++) { const bl = pk[6 + 2 * k] | (pk[7 + 2 * k] << 8); if (o + bl > len) fail('An E57 data packet is corrupt.'); const st = need.get(k); if (st) feed(st, pk.subarray(o, o + bl)); o += bl; }
      } else if (h[0] !== 0 && h[0] !== 2) fail('The E57 point section holds an unknown packet type.');
      l += len;
    }
    const got = Math.min(...states.filter((st) => st.f.bits).map((st) => st.out), m);
    if (states.some((st) => st.f.bits && st.rec < n)) warnings.push(`Scan "${name || scans.length + 1}" ends before its ${n} declared points; the points found were kept.`);
    // pose: unit quaternion (w, x, y, z) and translation
    const pose = kid(sc, 'pose'), rot = kid(pose, 'rotation'), tr = kid(pose, 'translation');
    let R = null;
    const T = [num(tr, 'x', 0), num(tr, 'y', 0), num(tr, 'z', 0)];
    if (rot) { const w = num(rot, 'w', 1), x = num(rot, 'x', 0), y = num(rot, 'y', 0), z = num(rot, 'z', 0), s = w * w + x * x + y * y + z * z || 1; R = [1 - (2 * (y * y + z * z)) / s, (2 * (x * y - w * z)) / s, (2 * (x * z + w * y)) / s, (2 * (x * y + w * z)) / s, 1 - (2 * (x * x + z * z)) / s, (2 * (y * z - w * x)) / s, (2 * (x * z - w * y)) / s, (2 * (y * z + w * x)) / s, 1 - (2 * (x * x + y * y)) / s]; }
    const before = kept;
    for (let i = 0; i < got; i++) {
      if (cols[3] && cols[3][i] !== 0) { invalid++; continue; }
      let x = cols[0][i], y = cols[1][i], z = cols[2] ? cols[2][i] : 0;
      if (useSph) { const r = x, az = y, el = z; x = r * Math.cos(el) * Math.cos(az); y = r * Math.cos(el) * Math.sin(az); z = r * Math.sin(el); }
      if (R) { const a = R[0] * x + R[1] * y + R[2] * z, b = R[3] * x + R[4] * y + R[5] * z, c = R[6] * x + R[7] * y + R[8] * z; x = a; y = b; z = c; }
      xyz[3 * kept] = x + T[0]; xyz[3 * kept + 1] = y + T[1]; xyz[3 * kept + 2] = z + T[2]; kept++;
    }
    scans.push({ name, points: n, kept: kept - before, fields: fields.map((f) => f.name), posed: !!pose });
    gidx += n;
  }
  if (step > 1) warnings.push(`Point cloud of ${total} points evenly sub-sampled to ${kept} (every ${step}th point).`);
  if (invalid) warnings.push(`${invalid}${step > 1 ? ' sampled' : ''} points flagged invalid by the scanner were dropped.`);
  if (images) warnings.push(`${images} images stored in the file were ignored.`);
  return { xyz: xyz.subarray(0, 3 * kept), total, kept, invalid, scans, images, version: `${major}.${minor}`, warnings };
}
