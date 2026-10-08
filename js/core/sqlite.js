// Minimal read-only SQLite 3 database-file reader (file format as documented at sqlite.org/fileformat2.html), written for
// GeoPackage import. No SQL engine: tables are scanned in rowid order and filtered by the caller.
//
// Read:   the 100-byte header (page size, reserved bytes, text encoding), table b-trees (interior and leaf pages) and
//         index b-trees (which also store WITHOUT ROWID tables), overflow page chains, the record format (all serial
//         types), sqlite_master / sqlite_schema, column names taken from the CREATE TABLE text (INTEGER PRIMARY KEY columns
//         aliasing the rowid), UTF-8, UTF-16LE and UTF-16BE text.
// Not read: a write-ahead log beside the file is NOT replayed (changes not yet checkpointed into the database file are
//         invisible; `wal` reports that the header announces WAL mode), rollback journals, encrypted databases, virtual
//         tables (their shadow tables are ordinary tables), views and triggers (listed, not evaluated), freelist pages.
// File content is untrusted: page numbers, cell offsets and payload sizes are bounds-checked; page visits are capped and
// each page is entered once per scan, so a cyclic or damaged tree ends in a clean error.

function fail(msg) { const e = new Error(msg); e.user = true; throw e; }
const utf8 = new TextDecoder('utf-8'), utf16le = new TextDecoder('utf-16le'), utf16be = new TextDecoder('utf-16be');

export const isSQLite = (u8) => u8.length >= 100 && utf8.decode(u8.subarray(0, 16)) === 'SQLite format 3\0';

/** Column names of a CREATE TABLE statement, the index of an INTEGER PRIMARY KEY (rowid alias) and the PRIMARY KEY columns. */
export function tableColumns(sql) {
  const s = String(sql || ''), open = s.indexOf('('), cols = [], pk = [];
  let rowidAlias = -1;
  if (open < 0) return { cols, rowidAlias, pk, withoutRowid: false };
  let depth = 0, q = '', start = open + 1, end = s.length;
  const parts = [];
  for (let i = open; i < s.length; i++) {
    const c = s[i];
    if (q) { if (c === q) q = ''; continue; }
    if (c === '"' || c === "'" || c === '`') q = c;
    else if (c === '[') q = ']';
    else if (c === '(') depth++;
    else if (c === ')') { if (--depth === 0) { end = i; break; } }
    else if (c === ',' && depth === 1) { parts.push(s.slice(start, i)); start = i + 1; }
  }
  parts.push(s.slice(start, end));
  const ident = (t) => { const m = /^\s*(?:"((?:[^"]|"")*)"|`([^`]*)`|\[([^\]]*)\]|'([^']*)'|([^\s(,]+))/.exec(t); return m ? (m[1] !== undefined ? m[1].replace(/""/g, '"') : m[2] ?? m[3] ?? m[4] ?? m[5]) : null; };
  for (const part of parts.slice(0, 4000)) {
    const head = /^\s*(CONSTRAINT|PRIMARY\s+KEY|UNIQUE|CHECK|FOREIGN\s+KEY)\b/i.exec(part);
    if (head) { const m = /PRIMARY\s+KEY\s*\(([^)]*)\)/i.exec(part); if (m) for (const c of m[1].split(',')) { const n = ident(c); if (n) pk.push(n); } continue; }
    const name = ident(part);
    if (name === null) continue;
    if (/\bPRIMARY\s+KEY\b/i.test(part)) { pk.push(name); if (/^\s*(?:"(?:[^"]|"")*"|`[^`]*`|\[[^\]]*\]|[^\s(,]+)\s+INTEGER\b(?!\s*\()/i.test(part) && !/\bPRIMARY\s+KEY\s+DESC\b/i.test(part)) rowidAlias = cols.length; }
    cols.push(name);
  }
  if (rowidAlias < 0 && pk.length === 1) { const k = cols.indexOf(pk[0]), def = parts.find((p) => ident(p) === pk[0]) || ''; if (k >= 0 && /^\s*\S+\s+INTEGER\s*(,|$|\b(?!\())/i.test(def + ',')) rowidAlias = k; }
  return { cols, rowidAlias, pk, withoutRowid: /\)\s*(?:STRICT\s*,\s*)?WITHOUT\s+ROWID/i.test(s.slice(end)) };
}

/**
 * Open a database held in memory. Returns:
 *   pageSize, pages, encoding ('utf-8' | 'utf-16le' | 'utf-16be'), wal (header announces write-ahead logging)
 *   schema                       [{ type, name, tbl_name, rootpage, sql }] (sqlite_master)
 *   tables                       Map name → { root, cols, rowidAlias, withoutRowid }
 *   scan(root, fn(values, rowid) → false to stop, { index })   visits every record of a b-tree in key order
 *   rows(name, { max, columns }) → [{ column: value }]  (INTEGER PRIMARY KEY filled from the rowid; BLOBs are Uint8Array)
 */
export function openSQLite(u8, opts = {}) {
  if (!isSQLite(u8)) fail('Not an SQLite database (header "SQLite format 3" missing).');
  const dv = new DataView(u8.buffer, u8.byteOffset, u8.byteLength), ps0 = dv.getUint16(16), pageSize = ps0 === 1 ? 65536 : ps0, U = pageSize - u8[20], enc = dv.getUint32(56) || 1;
  if (pageSize < 512 || pageSize & (pageSize - 1) || U < 480) fail('The SQLite header is corrupt (page size).');
  if (enc < 1 || enc > 3) fail('The SQLite header is corrupt (text encoding).');
  const pages = Math.floor(u8.length / pageSize), text = enc === 1 ? utf8 : enc === 2 ? utf16le : utf16be, maxVisits = opts.maxPages ?? 4e6;
  const bad = (what) => fail(`The SQLite database is truncated or corrupt (${what}).`);
  const pageOff = (n) => { if (!(n >= 1 && n <= pages)) bad(`page ${n} of ${pages}`); return (n - 1) * pageSize; };
  /** varint at p → [value, next offset] (9 bytes at most; values beyond 2^53 lose precision). */
  const varint = (p, end) => { let v = 0; for (let k = 0; k < 9; k++) { if (p >= end) bad('varint'); const b = u8[p++]; if (k === 8) return [v * 256 + b, p]; v = v * 128 + (b & 127); if (!(b & 128)) return [v, p]; } return [v, p]; };
  /** Payload of P bytes whose local part starts at p inside a page ending at pe; follows the overflow chain. */
  const payload = (p, pe, P, index) => {
    const X = index ? Math.floor(((U - 12) * 64) / 255) - 23 : U - 35, M = Math.floor(((U - 12) * 32) / 255) - 23;
    if (P <= X) { if (p + P > pe) bad('cell'); return u8.subarray(p, p + P); }
    if (P > 512e6) fail('A record in the SQLite database is too large.');
    const K = M + ((P - M) % (U - 4)), local = K <= X ? K : M;
    if (p + local + 4 > pe) bad('cell');
    const out = new Uint8Array(P);
    out.set(u8.subarray(p, p + local));
    let o = local, next = dv.getUint32(p + local);
    for (let guard = 0; o < P; guard++) {
      if (!next || guard > pages) bad('overflow chain');
      const q = pageOff(next), n = Math.min(P - o, U - 4);
      out.set(u8.subarray(q + 4, q + 4 + n), o);
      o += n; next = dv.getUint32(q);
    }
    return out;
  };
  /** Record bytes → array of values (null, number, string, Uint8Array). */
  const record = (b) => {
    const d = new DataView(b.buffer, b.byteOffset, b.byteLength), n = b.length, out = [];
    let h = 0, hs = 0, k = 0;
    for (; k < 9 && h < n; k++) { const c = b[h++]; if (k === 8) { hs = hs * 256 + c; break; } hs = hs * 128 + (c & 127); if (!(c & 128)) break; }
    if (hs > n || hs < h) bad('record header');
    let p = hs;
    while (h < hs && out.length < 32768) {
      let t = 0;
      for (k = 0; k < 9 && h < hs; k++) { const c = b[h++]; if (k === 8) { t = t * 256 + c; break; } t = t * 128 + (c & 127); if (!(c & 128)) break; }
      const len = t >= 12 ? Math.floor((t - 12) / 2) : [0, 1, 2, 3, 4, 6, 8, 8, 0, 0, 0, 0][t];
      if (p + len > n) bad('record body');
      if (t === 0 || t === 10 || t === 11) out.push(null);
      else if (t === 1) out.push(d.getInt8(p));
      else if (t === 2) out.push(d.getInt16(p));
      else if (t === 3) out.push((d.getInt8(p) << 16) | d.getUint16(p + 1));
      else if (t === 4) out.push(d.getInt32(p));
      else if (t === 5) out.push(d.getInt16(p) * 4294967296 + d.getUint32(p + 2));
      else if (t === 6) out.push(Number(d.getBigInt64(p)));
      else if (t === 7) out.push(d.getFloat64(p));
      else if (t === 8 || t === 9) out.push(t - 8);
      else if (t & 1) out.push(text.decode(b.subarray(p, p + len)));
      else out.push(b.subarray(p, p + len));
      p += len;
    }
    return out;
  };
  let visits = 0;
  function scan(root, fn) {
    const seen = new Set();
    const walk = (pg, depth) => {
      if (depth > 64 || seen.has(pg)) bad('b-tree loop');
      if (++visits > maxVisits) fail('The SQLite database has too many pages to scan.');
      seen.add(pg);
      const base = pageOff(pg), h = base + (pg === 1 ? 100 : 0), type = u8[h], pe = base + U, interior = type === 2 || type === 5, index = type === 2 || type === 10;
      if (![2, 5, 10, 13].includes(type)) bad(`page ${pg} is not a b-tree page`);
      const nc = dv.getUint16(h + 3), cp = h + (interior ? 12 : 8);
      if (cp + 2 * nc > pe) bad('cell pointers');
      for (let i = 0; i < nc; i++) {
        let c = base + dv.getUint16(cp + 2 * i);
        if (c < cp + 2 * nc || c + 4 > pe) bad('cell offset');
        if (interior) { if (walk(dv.getUint32(c), depth + 1) === false) return false; c += 4; if (!index) continue; }
        let P, rowid = null;
        [P, c] = varint(c, pe);
        if (!index) [rowid, c] = varint(c, pe);
        if (fn(record(payload(c, pe, P, index)), rowid) === false) return false;
      }
      return interior ? walk(dv.getUint32(h + 8), depth + 1) : true;
    };
    return walk(root, 0);
  }
  const schema = [], tables = new Map();
  scan(1, (v) => { if (schema.length >= 1e5) fail('The SQLite schema holds too many objects.'); schema.push({ type: String(v[0]), name: String(v[1]), tbl_name: String(v[2]), rootpage: +v[3] || 0, sql: v[4] === null || v[4] === undefined ? '' : String(v[4]) }); });
  for (const o of schema) if (o.type === 'table' && o.rootpage > 0) tables.set(o.name, { root: o.rootpage, ...tableColumns(o.sql) });
  function rows(name, { max = 1e6, where } = {}) {
    const t = tables.get(name) || [...tables].filter(([n]) => n.toLowerCase() === String(name).toLowerCase()).map(([, v]) => v)[0];
    if (!t) return [];
    const out = [], order = t.withoutRowid ? [...t.pk.map((c) => t.cols.indexOf(c)).filter((k) => k >= 0), ...t.cols.map((_, k) => k).filter((k) => !t.pk.includes(t.cols[k]))] : null;
    scan(t.root, (v, rowid) => {
      const r = {};
      if (order) order.forEach((ci, k) => { r[t.cols[ci]] = v[k] === undefined ? null : v[k]; });
      else for (let k = 0; k < t.cols.length; k++) r[t.cols[k]] = k === t.rowidAlias && (v[k] === null || v[k] === undefined) ? rowid : v[k] === undefined ? null : v[k];
      if (!order && t.rowidAlias < 0) r.rowid = rowid;
      if (!where || where(r)) out.push(r);
      return out.length < max;
    });
    return out;
  }
  return { pageSize, pages, encoding: text.encoding, wal: u8[18] === 2 || u8[19] === 2, schema, tables, scan, rows };
}
