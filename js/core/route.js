// Pipeline, riser, well and network geometry derived from imported Geometry objects (see geom.js): elevation profiles for
// 1-D multiphase models, draping of a route over bathymetry, minimum-curvature well paths, node-edge networks,
// wall-thickness / deposit maps and terrain grids. Pure functions: no DOM, no I/O; SI units (m, degrees) throughout.
import { gridOf, tableColumns, minimumCurvature as surveyPositions } from './geom.js';
import { clamp } from './num.js';

const D2R = Math.PI / 180, R_EARTH = 6371008.8, fin = Number.isFinite;
const UNIT = { m: 1, km: 1000, mm: 1e-3, cm: 0.01, ft: 0.3048, feet: 0.3048, usft: 1200 / 3937, in: 0.0254, inch: 0.0254 };
const isGeo = (g) => !!(g && (g.geographic || (g.grid && g.grid.geographic) || (g.stats && g.stats.geographic)));
const numOf = (v) => (typeof v === 'number' ? v : typeof v === 'string' && v.trim() !== '' ? Number(v) : NaN);
/** Great-circle distance (m) between two lon/lat points in degrees. */
export function haversine(lon1, lat1, lon2, lat2) {
  const a = Math.sin(((lat2 - lat1) * D2R) / 2) ** 2 + Math.cos(lat1 * D2R) * Math.cos(lat2 * D2R) * Math.sin(((lon2 - lon1) * D2R) / 2) ** 2;
  return 2 * R_EARTH * Math.asin(Math.min(1, Math.sqrt(a)));
}

/**
 * Minimum-curvature positions and dog-leg severity of a directional survey.
 * stations: [{ md, inc, azi }] | [[md, inc, azi]] | { md[], inc[], azi[] } (m, deg); start: { tvd, north, east } of station 0.
 * Returns { md, inc, azi, tvd, north, east, dogleg (deg per interval), dls (deg / 30 m), x (= east), y (= north), z (= -tvd) }.
 */
export function minimumCurvature(stations, start) {
  const s = surveyPositions(stations, start || {});
  return { ...s, x: s.east.slice(), y: s.north.slice(), z: s.tvd.map((v) => 0 - v) };
}

// ---- Profile object -------------------------------------------------------------------------------------------
/** x (cumulative horizontal distance) and z (elevation) -> profile; consecutive duplicate nodes are dropped. */
function build(x, z, source, warnings = []) {
  const X = [], Z = [];
  let dup = 0, bad = 0;
  for (let i = 0; i < x.length; i++) {
    if (!fin(x[i]) || !fin(z[i])) { bad++; continue; }
    if (X.length && x[i] === X[X.length - 1] && z[i] === Z[Z.length - 1]) { dup++; continue; }
    X.push(x[i]); Z.push(z[i]);
  }
  if (X.length < 2) throw new Error('A profile needs at least two distinct points.');
  if (bad) warnings.push(`${bad} points with non-numeric values were dropped.`);
  if (dup) warnings.push(`${dup} repeated points (zero-length elements) were dropped.`);
  const x0 = X[0], s = [0], incl = [];
  let zMin = Z[0], zMax = Z[0];
  for (let i = 0; i < X.length; i++) X[i] -= x0;
  for (let i = 1; i < X.length; i++) {
    const dx = X[i] - X[i - 1], dz = Z[i] - Z[i - 1];
    s.push(s[i - 1] + Math.hypot(dx, dz)); incl.push(Math.atan2(dz, dx) / D2R);
    if (Z[i] < zMin) zMin = Z[i]; if (Z[i] > zMax) zMax = Z[i];
  }
  return { x: X, z: Z, s, incl, length: s[s.length - 1], horizontalLength: X[X.length - 1], zMin, zMax, source, warnings };
}
/** 3-D path -> unrolled (x, z): x accumulates the horizontal distance between consecutive points. */
function unroll(X, Y, Z, geo) {
  const x = [0];
  for (let i = 1; i < X.length; i++) x.push(x[i - 1] + (geo ? haversine(X[i - 1], Y[i - 1], X[i], Y[i]) : Math.hypot(X[i] - X[i - 1], Y[i] - Y[i - 1])));
  return { x, z: Z };
}
/** Along-pipe lengths and elevations -> unrolled (x, z) that keeps every element length. */
function fromLengths(s, z, warnings, what = 'intervals') {
  const x = [0];
  let steep = 0;
  for (let i = 1; i < s.length; i++) { const ds = s[i] - s[i - 1], dz = z[i] - z[i - 1]; if (Math.abs(dz) > ds * (1 + 1e-9) + 1e-12) steep++; x.push(x[i - 1] + Math.sqrt(Math.max(0, ds * ds - dz * dz))); }
  if (steep) warnings.push(`${steep} ${what} change elevation by more than their length; they are drawn vertical and the profile is longer than the stated length there.`);
  return { x, z };
}

// ---- Douglas-Peucker in the x-z plane --------------------------------------------------------------------------
/** Weight of every vertex = the tolerance below which Douglas-Peucker keeps it (ends, highest and lowest point: Infinity). */
function dpWeights(x, z) {
  const n = x.length, w = new Float64Array(n);
  let iLo = 0, iHi = 0;
  for (let i = 1; i < n; i++) { if (z[i] < z[iLo]) iLo = i; if (z[i] > z[iHi]) iHi = i; }
  const forced = [...new Set([0, iLo, iHi, n - 1])].sort((a, b) => a - b), stack = [];
  for (const i of forced) w[i] = Infinity;
  for (let k = 0; k + 1 < forced.length; k++) stack.push([forced[k], forced[k + 1], Infinity]);
  while (stack.length) {
    const [a, b, cap] = stack.pop();
    if (b - a < 2) continue;
    const dx = x[b] - x[a], dz = z[b] - z[a], len = Math.hypot(dx, dz);
    let best = -1, bi = -1;
    for (let i = a + 1; i < b; i++) { const d = len > 0 ? Math.abs(dz * (x[i] - x[a]) - dx * (z[i] - z[a])) / len : Math.hypot(x[i] - x[a], z[i] - z[a]); if (d > best) { best = d; bi = i; } }
    const wt = Math.min(best, cap);
    w[bi] = wt;
    stack.push([a, bi, wt], [bi, b, wt]);
  }
  return w;
}
const pickNodes = (p, keep, note) => { const x = [], z = []; for (let i = 0; i < p.x.length; i++) if (keep(i)) { x.push(p.x[i]); z.push(p.z[i]); } return build(x, z, p.source, [...(p.warnings || []), ...(note ? [note(x.length)] : [])]); };
const asProfile = (p) => { if (!p || !Array.isArray(p.x) || !Array.isArray(p.z) || p.x.length !== p.z.length || p.x.length < 2) throw new Error('A profile { x[], z[] } with at least two points is required.'); return p; };

/** Douglas-Peucker simplification in x-z with tolerance (m); the end points and the highest and lowest point always stay. */
export function simplify(profile, tolerance = 0) {
  const p = asProfile(profile), tol = Math.max(0, +tolerance || 0), w = dpWeights(p.x, p.z);
  return pickNodes(p, (i) => w[i] > tol);
}
/** The n most significant points of a profile (Douglas-Peucker order); used to cap profile size. */
function simplifyTo(p, n, note) {
  if (p.x.length <= n) return p;
  const w = dpWeights(p.x, p.z), sorted = Array.from(w).sort((a, b) => b - a), thr = sorted[Math.max(1, n) - 1];
  let left = n - sorted.filter((v) => v > thr).length;
  return pickNodes(p, (i) => w[i] > thr || (w[i] === thr && left-- > 0), note);
}
/** n points equally spaced along the pipe (linear interpolation of x and z in arc length). */
export function resample(profile, n) {
  const p = asProfile(profile), src = p.s && p.s.length === p.x.length ? p : build(p.x, p.z, p.source, []), m = Math.max(2, Math.round(+n) || 2), x = [], z = [];
  for (let k = 0, i = 1; k < m; k++) {
    const sk = (src.length * k) / (m - 1);
    while (i < src.s.length - 1 && src.s[i] < sk) i++;
    const t = clamp((sk - src.s[i - 1]) / (src.s[i] - src.s[i - 1] || 1), 0, 1);
    x.push(src.x[i - 1] + t * (src.x[i] - src.x[i - 1])); z.push(src.z[i - 1] + t * (src.z[i] - src.z[i - 1]));
  }
  return build(x, z, src.source, [...(src.warnings || [])]);
}
/** Pipe elements of a profile: [{ x0, x1, z0, z1, length, incl }] with incl in degrees, positive uphill. */
export function segments(profile) {
  const p = asProfile(profile), out = [];
  for (let i = 1; i < p.x.length; i++) { const dx = p.x[i] - p.x[i - 1], dz = p.z[i] - p.z[i - 1]; out.push({ x0: p.x[i - 1], x1: p.x[i], z0: p.z[i - 1], z1: p.z[i], length: Math.hypot(dx, dz), incl: Math.atan2(dz, dx) / D2R }); }
  return out;
}

// ---- Centreline extraction ---------------------------------------------------------------------------------------
const lineLength = (p, geo) => { let l = 0; for (let i = 1; i < p.x.length; i++) l += geo ? haversine(p.x[i - 1], p.y[i - 1], p.x[i], p.y[i]) : Math.hypot(p.x[i] - p.x[i - 1], p.y[i] - p.y[i - 1], p.z ? p.z[i] - p.z[i - 1] : 0); return l; };
/** Longest polyline with every polyline that touches its ends chained on; returns { x, y, z | null, joined, left }. */
function mainLine(g) {
  const geo = isGeo(g), pls = (g.polylines || []).filter((p) => p && p.x && p.x.length >= 2).slice(0, 5000).map((p) => ({ x: p.x, y: p.y, z: p.z || null, closed: !!p.closed, len: lineLength(p, geo) })).sort((a, b) => b.len - a.len);
  if (!pls.length) throw new Error('The geometry holds no line with two or more points.');
  const first = pls[0], hasZ = pls.every((p) => p.z), cur = { x: first.x.slice(), y: first.y.slice(), z: hasZ || first.z ? first.z.slice() : null }, used = new Set([0]);
  let ext = 0;
  for (const p of pls) for (const a of [p.x, p.y]) for (const v of a) ext = Math.max(ext, Math.abs(v));
  const tol = 1e-7 * (ext || 1), near = (ax, ay, bx, by) => Math.abs(ax - bx) <= tol && Math.abs(ay - by) <= tol;
  for (let grown = true; grown && used.size < pls.length && !first.closed;) {
    grown = false;
    for (let k = 1; k < pls.length; k++) {
      if (used.has(k) || pls[k].closed) continue;
      const p = pls[k], n = p.x.length, m = cur.x.length, zOf = (rev) => (cur.z ? (p.z ? (rev ? p.z.slice().reverse() : p.z) : null) : null);
      let add = null, front = false, rev = false;
      if (near(cur.x[m - 1], cur.y[m - 1], p.x[0], p.y[0])) add = true;
      else if (near(cur.x[m - 1], cur.y[m - 1], p.x[n - 1], p.y[n - 1])) { add = true; rev = true; }
      else if (near(cur.x[0], cur.y[0], p.x[n - 1], p.y[n - 1])) { add = true; front = true; }
      else if (near(cur.x[0], cur.y[0], p.x[0], p.y[0])) { add = true; front = true; rev = true; }
      if (!add) continue;
      const px = rev ? p.x.slice().reverse() : p.x, py = rev ? p.y.slice().reverse() : p.y, pz = zOf(rev);
      if (cur.z && !pz) cur.z = null;
      if (front) { cur.x = px.slice(0, -1).concat(cur.x); cur.y = py.slice(0, -1).concat(cur.y); if (cur.z) cur.z = pz.slice(0, -1).concat(cur.z); }
      else { cur.x = cur.x.concat(px.slice(1)); cur.y = cur.y.concat(py.slice(1)); if (cur.z) cur.z = cur.z.concat(pz.slice(1)); }
      used.add(k); grown = true;
    }
  }
  return { ...cur, closed: first.closed, joined: used.size, left: pls.length - used.size, geo };
}
const monotone = (a) => { let up = true, dn = true; for (let i = 1; i < a.length; i++) { if (a[i] < a[i - 1]) up = false; if (a[i] > a[i - 1]) dn = false; } return up || dn; };
/** Order of scattered centreline points along their principal horizontal axis. */
function principalOrder(X, Y) {
  const n = X.length;
  let mx = 0, my = 0, sxx = 0, sxy = 0, syy = 0;
  for (let i = 0; i < n; i++) { mx += X[i] / n; my += Y[i] / n; }
  for (let i = 0; i < n; i++) { const a = X[i] - mx, b = Y[i] - my; sxx += a * a; sxy += a * b; syy += b * b; }
  const th = 0.5 * Math.atan2(2 * sxy, sxx - syy), c = Math.cos(th), s = Math.sin(th);
  return Array.from({ length: n }, (_, i) => i).sort((p, q) => (X[p] - X[q]) * c + (Y[p] - Y[q]) * s);
}
function pathLike(X, Y) {
  const n = X.length;
  let lo = [Infinity, Infinity], hi = [-Infinity, -Infinity], sum = 0;
  for (let i = 0; i < n; i++) { lo = [Math.min(lo[0], X[i]), Math.min(lo[1], Y[i])]; hi = [Math.max(hi[0], X[i]), Math.max(hi[1], Y[i])]; if (i) sum += Math.hypot(X[i] - X[i - 1], Y[i] - Y[i - 1]); }
  const diag = Math.hypot(hi[0] - lo[0], hi[1] - lo[1]);
  return n >= 2 && sum <= 3 * diag;
}

/** Longest route through a network: source-to-sink on an acyclic graph, otherwise the longest shortest path. */
function networkPath(g, W) {
  const { nodes, edges } = g.network, at = new Map(nodes.map((n) => [n.id, n])), E = [];
  let noLen = 0;
  for (const e of edges) {
    const a = at.get(e.from), b = at.get(e.to);
    if (!a || !b || e.from === e.to) continue;
    let len = fin(e.length) && e.length > 0 ? e.length : NaN;
    if (!fin(len) && [a.x, a.y, b.x, b.y].every(fin)) len = Math.hypot(b.x - a.x, b.y - a.y, (fin(b.z) ? b.z : 0) - (fin(a.z) ? a.z : 0));
    if (!fin(len) || len <= 0) { noLen++; continue; }
    E.push({ from: e.from, to: e.to, len });
  }
  if (!E.length) throw new Error('The network has no connection with a length or with node coordinates to measure one.');
  if (noLen) W.push(`${noLen} connections without a length or node coordinates were left out of the route.`);
  const out = new Map(), indeg = new Map(), und = new Map(), put = (m, k, v) => { const l = m.get(k); if (l) l.push(v); else m.set(k, [v]); };
  for (const e of E) { put(out, e.from, e); indeg.set(e.to, (indeg.get(e.to) || 0) + 1); if (!indeg.has(e.from)) indeg.set(e.from, indeg.get(e.from) || 0); put(und, e.from, [e.to, e.len]); put(und, e.to, [e.from, e.len]); }
  const ids = [...indeg.keys()], deg = new Map(indeg), ready = ids.filter((id) => !deg.get(id)), order = [];
  for (let q = 0; q < ready.length; q++) { const id = ready[q]; order.push(id); for (const e of out.get(id) || []) { deg.set(e.to, deg.get(e.to) - 1); if (!deg.get(e.to)) ready.push(e.to); } }
  const walk = (end, prev) => { const pid = [end], lens = []; for (let c = end; prev.has(c);) { const [q, l] = prev.get(c); lens.push(l); pid.push(q); c = q; } return { pid: pid.reverse(), lens: lens.reverse() }; };
  let path;                                         // { pid: node ids along the route, lens: length of each connection }
  if (order.length === ids.length) {
    const dist = new Map(ids.map((id) => [id, 0])), prev = new Map();
    for (const id of order) for (const e of out.get(id) || []) if (dist.get(id) + e.len > dist.get(e.to)) { dist.set(e.to, dist.get(id) + e.len); prev.set(e.to, [id, e.len]); }
    let end = ids[0];
    for (const id of ids) if (dist.get(id) > dist.get(end)) end = id;
    path = walk(end, prev);
  } else {
    W.push('The network has loops or no clear flow direction: the route is the longest of the shortest paths between two nodes, ignoring direction.');
    if (ids.length > 4000) throw new Error('The cyclic network is too large to search for a route.');
    const far = (src) => {
      const dist = new Map([[src, 0]]), prev = new Map(), done = new Set();
      for (;;) {
        let u = null, du = Infinity;
        for (const [k, d] of dist) if (!done.has(k) && d < du) { u = k; du = d; }
        if (u === null) break;
        done.add(u);
        for (const [v, l] of und.get(u) || []) if (!done.has(v) && du + l < (dist.get(v) ?? Infinity)) { dist.set(v, du + l); prev.set(v, [u, l]); }
      }
      let end = src;
      for (const [k, d] of dist) if (d > dist.get(end)) end = k;
      return { end, prev };
    };
    const a = far(ids[0]).end, b = far(a);
    path = walk(b.end, b.prev);
  }
  const s = [0], z = [];
  let noZ = 0;
  path.pid.forEach((id, i) => { const n = at.get(id); if (i) s.push(s[i - 1] + path.lens[i - 1]); if (fin(n.z)) z.push(n.z); else { z.push(0); noZ++; } });
  if (noZ) W.push(`${noZ} of ${path.pid.length} nodes on the route have no elevation; it is taken as 0 there.`);
  return { ...fromLengths(s, z, W, 'connections'), source: `network: longest route ${path.pid[0]} → ${path.pid[path.pid.length - 1]} over ${path.lens.length} connections (declared lengths, node elevations)`, nodes: path.pid };
}

/** Table columns -> unrolled path. */
function tablePath(g, o, W) {
  const H = g.headers || [], recs = g.records || [];
  if (!H.length || recs.length < 2) throw new Error('The table needs a header row and at least two rows.');
  const tc = tableColumns(H), c = tc.c, has = (k) => c[k] !== undefined, raw = (i) => recs.map((r) => numOf(r[H[i]])), col = (k) => raw(c[k]), U = (k, d) => UNIT[tc.unit[k]] ?? d, dpd = o.depthPositiveDown;
  const zCol = () => {
    if (has('z') && !(has('depth') && c.depth === c.z)) return col('z').map((v) => v * U('z', o.scale) * (dpd === true ? -1 : 1));
    if (has('depth')) { if (dpd !== false) W.push(`Column "${H[c.depth]}" is depth (positive down) and was converted to negative elevation.`); return col('depth').map((v) => v * U('depth', o.scale) * (dpd === false ? 1 : -1)); }
    return null;
  };
  if (tc.role === 'survey') {
    const um = U('md', o.scale), md = col('md').map((v) => v * um), rows = md.map((_, i) => i).filter((i) => fin(md[i]) && (has('inc') ? fin(numOf(recs[i][H[c.inc]])) : fin(numOf(recs[i][H[c.tvd]]))));
    let tvd;
    if (has('inc')) { const inc = col('inc'), azi = has('azi') ? col('azi') : null; tvd = surveyPositions(rows.map((i) => [md[i], inc[i], azi && fin(azi[i]) ? azi[i] : 0])).tvd; }
    else { const t = col('tvd'), ut = U('tvd', um); tvd = rows.map((i) => t[i] * ut); }
    return { ...fromLengths(rows.map((i) => md[i]), tvd.map((v) => (dpd === false ? v : -v)), W, 'survey intervals'), source: `table: well survey (${has('inc') ? 'MD + inclination' + (has('azi') ? ' + azimuth' : '') + ', minimum curvature' : 'MD + TVD'})`, well: true };
  }
  const z = zCol();
  if (has('chainage') && z && !(has('lon') && has('lat'))) {
    const kp = /^(kp|kilomet)/.test(String(H[c.chainage]).toLowerCase().replace(/[^a-z]/g, '')), uc = U('chainage', kp ? 1000 : o.scale), ch = col('chainage').map((v) => v * uc);
    if (kp && !tc.unit.chainage) W.push(`Column "${H[c.chainage]}" is a kilometre point and was converted to metres.`);
    let idx = ch.map((_, i) => i).filter((i) => fin(ch[i]) && fin(z[i]));
    if (!monotone(idx.map((i) => ch[i]))) { idx = idx.sort((p, q) => ch[p] - ch[q]); W.push('Rows were sorted by chainage.'); }
    return { x: idx.map((i) => Math.abs(ch[i] - ch[idx[0]])), z: idx.map((i) => z[i]), source: `table: ${H[c.chainage]} + ${H[has('z') ? c.z : c.depth]} (chainage taken as horizontal distance)` };
  }
  if (has('x') && has('y')) {
    const geo = has('lon') && has('lat'), X = col('x'), Y = col('y'), ux = geo ? 1 : U('x', o.scale), uy = geo ? 1 : U('y', o.scale), idx = X.map((_, i) => i).filter((i) => fin(X[i]) && fin(Y[i]) && (!z || fin(z[i])));
    if (!z) W.push('The table has no elevation or depth column; the route is flat (z = 0). Drape it over a bathymetry grid to get elevations.');
    return { ...unroll(idx.map((i) => X[i] * ux), idx.map((i) => Y[i] * uy), idx.map((i) => (z ? z[i] : 0)), geo), source: `table: ${H[c.x]}, ${H[c.y]}${z ? ', ' + H[has('z') ? c.z : c.depth] : ''}${geo ? ' (geographic, great-circle distances)' : ''}` };
  }
  const numeric = H.map((_, i) => i).filter((i) => { const v = raw(i); return v.filter(fin).length >= 0.9 * v.length; });
  if (numeric.length === 2) { const a = raw(numeric[0]), b = raw(numeric[1]); W.push(`Column names were not recognised: "${H[numeric[0]]}" is taken as distance and "${H[numeric[1]]}" as elevation.`); return { ...unroll(a.map((v) => v * o.scale), a.map(() => 0), b.map((v) => v * o.scale * (dpd === true ? -1 : 1)), false), source: 'table: two numeric columns (distance, elevation)' }; }
  if (numeric.length >= 3) { const [a, b, d] = numeric.slice(0, 3).map(raw); W.push(`Column names were not recognised: the first three numeric columns (${numeric.slice(0, 3).map((i) => H[i]).join(', ')}) are taken as x, y, z.`); return { ...unroll(a.map((v) => v * o.scale), b.map((v) => v * o.scale), d.map((v) => v * o.scale * (dpd === true ? -1 : 1)), false), source: 'table: first three numeric columns (x, y, z)' }; }
  throw new Error('The table has no columns that describe a route (chainage + elevation, x / y / z, lat / lon / depth, MD / TVD or MD / inclination / azimuth).');
}

/** Pipe-like surface mesh -> centreline: vertices sliced along the principal axis, one centroid per slice. */
function meshPath(g, o, W) {
  const t = g.triangles, nv = Math.floor(t.length / 3), m = [0, 0, 0], C = [[0, 0, 0], [0, 0, 0], [0, 0, 0]], step = Math.max(1, Math.floor(nv / 2e5));
  let cnt = 0;
  for (let i = 0; i < nv; i += step) { for (let k = 0; k < 3; k++) m[k] += t[3 * i + k]; cnt++; }
  for (let k = 0; k < 3; k++) m[k] /= cnt;
  for (let i = 0; i < nv; i += step) { const d = [t[3 * i] - m[0], t[3 * i + 1] - m[1], t[3 * i + 2] - m[2]]; for (let a = 0; a < 3; a++) for (let b = 0; b < 3; b++) C[a][b] += d[a] * d[b]; }
  let ax = [1, 0.37, 0.11];
  for (let it = 0; it < 80; it++) { const v = [0, 1, 2].map((a) => C[a][0] * ax[0] + C[a][1] * ax[1] + C[a][2] * ax[2]), l = Math.hypot(...v); if (!(l > 0)) break; ax = v.map((q) => q / l); }
  let lo = Infinity, hi = -Infinity;
  const proj = (i) => (t[3 * i] - m[0]) * ax[0] + (t[3 * i + 1] - m[1]) * ax[1] + (t[3 * i + 2] - m[2]) * ax[2];
  for (let i = 0; i < nv; i += step) { const p = proj(i); if (p < lo) lo = p; if (p > hi) hi = p; }
  if (!(hi > lo)) throw new Error('The mesh has no extent to take a centreline from.');
  const nb = clamp(Math.round(Math.sqrt(cnt) / 2), 4, Math.min(200, o.maxPoints)), S = Array.from({ length: nb }, () => [0, 0, 0, 0]);
  for (let i = 0; i < nv; i += step) { const b = Math.min(nb - 1, Math.floor(((proj(i) - lo) / (hi - lo)) * nb)), q = S[b]; q[0] += t[3 * i]; q[1] += t[3 * i + 1]; q[2] += t[3 * i + 2]; q[3]++; }
  const P = S.filter((q) => q[3]).map((q) => [(q[0] / q[3]) * o.scale, (q[1] / q[3]) * o.scale, (q[2] / q[3]) * o.scale * (o.depthPositiveDown === true ? -1 : 1)]);
  W.push(`Centreline ESTIMATED from the surface mesh: vertices were sliced into ${nb} stations along the principal axis and averaged. Bends are smoothed and the ends are shortened by up to half a slice; use the design centreline where one exists.`);
  return { ...unroll(P.map((q) => q[0]), P.map((q) => q[1]), P.map((q) => q[2]), false), source: `mesh: centreline estimated from ${nb} slices along the principal axis` };
}

function pathOf(g, o, W) {
  const sc = o.scale, zs = o.depthPositiveDown === true ? -1 : 1;
  if (g.kind === 'polylines') {
    if (g.survey && Array.isArray(g.survey.md) && g.survey.md.length > 1) return { ...fromLengths(g.survey.md, g.survey.tvd.map((v) => 0 - v), W, 'survey intervals'), source: 'well survey: measured depth and true vertical depth', well: true };
    const ml = mainLine(g), geo = ml.geo, sxy = geo ? 1 : sc;
    if (ml.closed) W.push('The line is a closed outline; it is followed once from its first point without the closing segment.');
    if (ml.left) W.push(`${ml.left} of ${ml.left + ml.joined} polylines do not connect to the main line and were ignored.`);
    if (ml.joined > 1) W.push(`${ml.joined} connected polylines were joined end to end.`);
    if (ml.z) return { ...unroll(ml.x.map((v) => v * sxy), ml.y.map((v) => v * sxy), ml.z.map((v) => v * sc * zs), geo), source: `3-D polyline${geo ? ' (geographic, great-circle distances)' : ''}` };
    const plane = o.plane === 'xz' || o.plane === 'xy' ? o.plane : geo || g.pathway === 'gis' || !monotone(ml.x) ? 'xy' : 'xz';
    if (plane === 'xz') {
      if (o.plane !== 'xz') W.push('The 2-D polyline is read as an elevation profile (x = horizontal distance, y = elevation) because x never turns back; pass plane: "xy" if it is a plan-view route.');
      return { ...unroll(ml.x.map((v) => v * sc), ml.x.map(() => 0), ml.y.map((v) => v * sc * zs), false), source: '2-D polyline read as an x–z elevation profile' };
    }
    W.push('The route has no elevations and is flat (z = 0); drape it over a bathymetry grid (drape(route, grid) or opts.grid) to get its profile.');
    return { ...unroll(ml.x.map((v) => v * sxy), ml.y.map((v) => v * sxy), ml.x.map(() => 0), geo), source: `plan-view polyline without elevations${geo ? ' (geographic)' : ''}` };
  }
  if (g.kind === 'points' || (g.kind === 'grid' && g.points && g.stats && g.stats.role === 'route3d')) {
    const n = Math.floor(g.points.length / 3), step = Math.max(1, Math.ceil(n / 50000)), geo = isGeo(g), sxy = geo ? 1 : sc;
    let X = [], Y = [], Z = [];
    for (let i = 0; i < n; i += step) { X.push(g.points[3 * i] * sxy); Y.push(g.points[3 * i + 1] * sxy); Z.push(g.points[3 * i + 2] * sc * zs); }
    if (step > 1) W.push(`Every ${step}th point of ${n} was used.`);
    let how = 'in file order';
    if (!((g.stats && g.stats.role === 'route3d') || pathLike(X, Y))) { const ord = principalOrder(X, Y); X = ord.map((i) => X[i]); Y = ord.map((i) => Y[i]); Z = ord.map((i) => Z[i]); how = 'sorted along their principal horizontal axis'; W.push('The points do not trace a line in file order; they were sorted along their principal horizontal axis. A scattered cloud is not a centreline.'); }
    return { ...unroll(X, Y, Z, geo), source: `point list ${how}${geo ? ' (geographic, great-circle distances)' : ''}` };
  }
  if (g.kind === 'network') return networkPath(g, W);
  if (g.kind === 'table') return tablePath(g, o, W);
  if (g.kind === 'mesh') return meshPath(g, o, W);
  if (g.kind === 'grid') throw new Error('An elevation grid is terrain, not a route. Drape a route over it with drape(route, grid).');
  throw new Error(`A ${g.kind || 'unknown'} geometry cannot give a pipeline profile.`);
}

/**
 * Elevation profile of a pipeline, riser, flowline or well for 1-D flow models.
 * Returns { x[] cumulative horizontal distance (m), z[] elevation (m, negative below sea level), s[] cumulative length along
 * the pipe (m), incl[] inclination of each element (deg, + uphill in the direction of the profile; one fewer than the nodes),
 * length, horizontalLength, zMin, zMax, source, warnings[] }.
 * Accepts 3-D polylines, 2-D x–z profiles, geographic polylines, ordered point lists, networks (longest route), tables
 * (chainage + elevation, x / y / z, lat / lon / depth, MD / TVD, MD / inc / azi), pipe-like meshes (estimate) and profiles.
 * opts: unit ('m' | 'km' | 'ft' | 'mm': unit of the coordinates when the file does not say), depthPositiveDown (true | false |
 * 'auto': only depth / TVD columns are flipped), reverse (default false; well surveys default to true so that the profile runs
 * from bottom-hole to wellhead), maxPoints (400), plane ('xz' | 'xy' for 2-D polylines), grid (bathymetry to drape a flat route on).
 */
export function profileFrom(g, opts = {}) {
  if (!g || typeof g !== 'object') throw new Error('profileFrom needs a geometry.');
  opts = opts && typeof opts === 'object' ? opts : {};
  const W = [], fileUnit = g.stats && typeof g.stats.units === 'string' ? g.stats.units.toLowerCase() : null, unit = opts.unit !== undefined ? String(opts.unit).toLowerCase() : fileUnit && UNIT[fileUnit] && !(g.stats && g.stats.wellSurvey) ? fileUnit : 'm';
  if (!UNIT[unit]) throw new Error(`Unknown length unit "${opts.unit}" (use m, km, ft or mm).`);
  const o = { scale: UNIT[unit], depthPositiveDown: opts.depthPositiveDown === true || opts.depthPositiveDown === false ? opts.depthPositiveDown : 'auto', maxPoints: Math.max(2, Math.round(+opts.maxPoints) || 400), plane: opts.plane };
  if (unit !== 'm' && opts.unit === undefined) W.push(`Coordinates converted from ${unit} to metres (unit given by the file).`);
  let p;
  if (!g.kind && Array.isArray(g.x) && Array.isArray(g.z)) p = { x: g.x.slice(), z: g.z.slice(), source: g.source || 'profile' };
  else {
    let src = g;
    if (opts.grid && (g.kind === 'polylines' || g.kind === 'points') && !(g.polylines || []).some((q) => q.z) && !(g.kind === 'points' && g.bbox && g.bbox.max[2] !== g.bbox.min[2])) { src = drape(g, opts.grid); W.push(...src.warnings); }
    p = pathOf(src, o, W);
    if (src !== g) p.source = `route draped over the elevation grid (${src.polylines[0].x.length} samples)`;
  }
  let { x, z } = p;
  const reverse = opts.reverse !== undefined ? !!opts.reverse : !!p.well;
  if (reverse) { const xe = x[x.length - 1]; x = x.map((v) => xe - v).reverse(); z = z.slice().reverse(); if (p.well && opts.reverse === undefined) W.push('Well survey reversed so that the profile runs from bottom-hole to wellhead (production direction); pass reverse: false to keep survey order.'); }
  const full = build(x, z, p.source, W);
  return simplifyTo(full, o.maxPoints, (n) => `Profile simplified from ${full.x.length} to ${n} points (length ${full.length.toPrecision(7)} m before); the end, highest and lowest points are kept.`);
}

// ---- Draping ---------------------------------------------------------------------------------------------------
const gridXYZ = (grid) => {
  const q = grid && grid.kind === 'grid' ? grid.grid : grid;
  const x = q && (q.x || q.lon), y = q && (q.y || q.lat), z = q && (q.z || q.elev);
  if (!x || !y || !z || x.length < 2 || y.length < 2 || z.length !== y.length) throw new Error('drape needs an elevation grid { x[], y[], z[][] } with at least 2 × 2 nodes.');
  return { x, y, z, geo: grid.kind === 'grid' ? !!grid.grid.geographic : !!(grid.geographic || grid.lon) };
};
/** Bilinear sample of a grid with ascending x and y; NaN cells are left out of the weighting; outside points are clamped. */
function sampleGrid(G, px, py) {
  const loc = (a, q) => { let lo = 0, hi = a.length - 1; while (hi - lo > 1) { const m = (lo + hi) >> 1; if (a[m] <= q) lo = m; else hi = m; } return [lo, hi, clamp((q - a[lo]) / (a[hi] - a[lo] || 1), 0, 1)]; };
  const [i0, i1, ti] = loc(G.x, px), [j0, j1, tj] = loc(G.y, py);
  let s = 0, w = 0;
  for (const [j, i, wt] of [[j0, i0, (1 - tj) * (1 - ti)], [j0, i1, (1 - tj) * ti], [j1, i0, tj * (1 - ti)], [j1, i1, tj * ti]]) { const v = G.z[j][i]; if (v === v && wt > 0) { s += wt * v; w += wt; } }
  return w > 0 ? s / w : NaN;
}
/**
 * Lay a plan-view route on an elevation grid. route: polylines or points Geometry, or { x[], y[] } (same coordinate system
 * as the grid: both geographic or both projected); gridGeometry: a 'grid' Geometry or { x[], y[], elev[][] | z[][], geographic }.
 * The route is densified to the grid spacing and each point gets the bilinearly interpolated elevation.
 * Returns a 'polylines' Geometry with one 3-D polyline, ready for profileFrom.
 */
export function drape(route, gridGeometry) {
  const G = gridXYZ(gridGeometry), warnings = [];
  let rx, ry, geo;
  if (route && route.kind === 'polylines') { const ml = mainLine(route); rx = ml.x; ry = ml.y; geo = ml.geo; if (ml.left) warnings.push(`${ml.left} polylines that do not connect to the main line were ignored.`); }
  else if (route && route.kind === 'points') { rx = []; ry = []; for (let i = 0; i + 2 < route.points.length; i += 3) { rx.push(route.points[i]); ry.push(route.points[i + 1]); } geo = isGeo(route); }
  else if (route && Array.isArray(route.x) && Array.isArray(route.y) && route.x.length === route.y.length) { rx = route.x; ry = route.y; geo = route.geographic !== undefined ? !!route.geographic : G.geo; }
  else throw new Error('drape needs a route: a polylines or points geometry, or { x[], y[] }.');
  if (rx.length < 2) throw new Error('The route needs at least two points.');
  if (geo !== G.geo) throw new Error(`The route is ${geo ? 'geographic (degrees)' : 'projected (metres)'} but the grid is ${G.geo ? 'geographic (degrees)' : 'projected (metres)'}; bring both into one coordinate system first.`);
  const dxg = Math.abs(G.x[G.x.length - 1] - G.x[0]) / (G.x.length - 1), dyg = Math.abs(G.y[G.y.length - 1] - G.y[0]) / (G.y.length - 1);
  let total = 0;
  for (let i = 1; i < rx.length; i++) total += Math.hypot(rx[i] - rx[i - 1], ry[i] - ry[i - 1]);
  const step = Math.max(Math.min(dxg, dyg) || total, total / 20000), x = [], y = [], z = [];
  let outside = 0, holes = 0;
  const put = (px, py) => {
    if (px < G.x[0] || px > G.x[G.x.length - 1] || py < G.y[0] || py > G.y[G.y.length - 1]) outside++;
    const v = sampleGrid(G, clamp(px, G.x[0], G.x[G.x.length - 1]), clamp(py, G.y[0], G.y[G.y.length - 1]));
    if (v !== v) { holes++; return; }
    x.push(px); y.push(py); z.push(v);
  };
  put(rx[0], ry[0]);
  for (let i = 1; i < rx.length; i++) { const d = Math.hypot(rx[i] - rx[i - 1], ry[i] - ry[i - 1]), m = Math.max(1, Math.ceil(d / step - 1e-9)); for (let k = 1; k <= m; k++) put(rx[i - 1] + ((rx[i] - rx[i - 1]) * k) / m, ry[i - 1] + ((ry[i] - ry[i - 1]) * k) / m); }
  if (x.length < 2) throw new Error('The route does not cross any valid cell of the elevation grid.');
  if (outside) warnings.push(`${outside} of ${x.length + holes} route samples lie outside the grid; they take the elevation of the nearest grid edge.`);
  if (holes) warnings.push(`${holes} route samples fall on cells without data and were left out.`);
  let lo = Infinity, hi = -Infinity;
  for (const v of z) { if (v < lo) lo = v; if (v > hi) hi = v; }
  return { kind: 'polylines', name: `${route.name || 'route'} on ${gridGeometry.name || 'grid'}`, format: 'Draped route', pathway: 'gis', polylines: [{ x, y, z, closed: false }], geographic: geo, bbox: { min: [Math.min(...x), Math.min(...y), lo], max: [Math.max(...x), Math.max(...y), hi] }, warnings, stats: { draped: true, polylines: 1, vertices: x.length, outside, routePoints: rx.length } };
}

// ---- Verification ------------------------------------------------------------------------------------------------
/**
 * Verification of a profile before it is used by a solver. opts.id = pipe inner diameter in metres (enables the volume check).
 * Returns [{ name, pass, got, expected, note }].
 */
export function checkProfile(profile, { id } = {}) {
  const p = asProfile(profile), seg = segments(p), n = p.x.length, L = p.s ? p.s[p.s.length - 1] : NaN, sum = seg.reduce((q, e) => q + e.length, 0), tolL = 1e-9 * Math.max(1, Math.abs(sum)), out = [];
  out.push({ name: 'Length conservation', pass: fin(L) && Math.abs(sum - L) <= tolL, got: sum, expected: L, note: 'Sum of element lengths equals the arc length at the last node.' });
  const zero = seg.filter((e) => !(e.length > 1e-12 * Math.max(1, sum))).length;
  out.push({ name: 'No zero-length elements', pass: zero === 0, got: zero, expected: 0, note: 'Every element has a positive length.' });
  const seen = new Set();
  let dup = 0;
  for (let i = 0; i < n; i++) { const k = p.x[i] + ',' + p.z[i]; if (seen.has(k)) dup++; else seen.add(k); }
  out.push({ name: 'No duplicate nodes', pass: dup === 0, got: dup, expected: 0, note: 'No two nodes share the same distance and elevation.' });
  let back = 0;
  for (let i = 1; i < n; i++) if (p.x[i] < p.x[i - 1] || (p.s && !(p.s[i] > p.s[i - 1]))) back++;
  out.push({ name: 'Monotonic chainage', pass: back === 0, got: back, expected: 0, note: 'Horizontal distance never decreases and arc length strictly increases.' });
  const bad = [...p.x, ...p.z, ...(p.s || []), ...(p.incl || [])].filter((v) => !fin(v)).length;
  out.push({ name: 'Finite values', pass: bad === 0 && (!p.s || p.s.length === n) && (!p.incl || p.incl.length === n - 1), got: bad, expected: 0, note: 'x, z, s and inclination are finite and of consistent length.' });
  if (fin(id) && id > 0) { const a = (Math.PI / 4) * id * id, vol = seg.reduce((q, e) => q + a * e.length, 0), exp = a * L; out.push({ name: 'Internal volume', pass: Math.abs(vol - exp) <= 1e-9 * Math.max(1, exp), got: vol, expected: exp, note: 'Sum of element volumes equals π/4 · id² · length (m³).' }); }
  const imax = seg.reduce((q, e) => Math.max(q, Math.abs(e.incl)), 0);
  out.push({ name: 'Inclination bounds', pass: imax <= 90 + 1e-9 && (!p.incl || p.incl.every((v, i) => Math.abs(v - seg[i].incl) <= 1e-6)), got: imax, expected: 90, note: 'Every element lies between −90° and +90° and matches the stored inclination.' });
  let rt = NaN;
  try { rt = profileFrom({ kind: 'table', headers: ['chainage', 'elevation'], records: p.x.map((x, i) => ({ chainage: x, elevation: p.z[i] })) }, { maxPoints: Math.max(n, 2) }).length; } catch { rt = NaN; }
  out.push({ name: 'Round trip', pass: fin(rt) && Math.abs(rt - sum) <= 1e-9 * Math.max(1, sum), got: rt, expected: sum, note: 'Profile → chainage / elevation table → profile reproduces the length.' });
  return out;
}

// ---- Networks ------------------------------------------------------------------------------------------------------
/**
 * Node-edge network of a geometry: a network as it is, polylines with every polyline an edge and coincident end points
 * welded into nodes (crossings between end points are not split), or a table with from / to columns.
 * Returns { nodes: [{ id, type, x, y, z }], edges: [{ from, to, type, length, diameter, name }],
 * issues: { disconnected: [node ids outside the largest connected part], duplicates: [edge names], zeroLength: [edge names], loops } }.
 */
export function networkFrom(g) {
  if (!g || typeof g !== 'object') throw new Error('networkFrom needs a geometry.');
  let nodes = [], edges = [];
  const num = (v) => (fin(v) ? v : null);
  if (g.kind === 'network') {
    nodes = g.network.nodes.map((n) => ({ id: String(n.id), type: n.type || 'node', x: num(n.x), y: num(n.y), z: num(n.z) }));
    edges = g.network.edges.map((e) => ({ from: String(e.from), to: String(e.to), type: e.type || 'pipe', length: num(e.length), diameter: num(e.diameter), name: e.name || `${e.from}-${e.to}` }));
  } else if (g.kind === 'polylines') {
    const geo = isGeo(g), u = geo ? 1 : UNIT[String((g.stats && g.stats.units) || 'm').toLowerCase()] ?? 1, pls = g.polylines.filter((p) => p.x.length >= 2).slice(0, 20000);
    let ext = 0;
    for (const p of pls) for (const a of [p.x, p.y]) for (const v of a) ext = Math.max(ext, Math.abs(v));
    const tol = 1e-7 * (ext || 1), cell = new Map();
    const nodeAt = (x, y, z) => {
      const kx = Math.round(x / tol / 4), ky = Math.round(y / tol / 4);
      for (let a = -1; a <= 1; a++) for (let b = -1; b <= 1; b++) { const q = cell.get(kx + a + ',' + (ky + b)); if (q && Math.abs(q.x0 - x) <= tol && Math.abs(q.y0 - y) <= tol) return q.id; }
      const n = { id: 'N' + (nodes.length + 1), type: 'node', x: geo ? x : x * u, y: geo ? y : y * u, z: z === null ? null : z * u };
      nodes.push(n); cell.set(kx + ',' + ky, { id: n.id, x0: x, y0: y });
      return n.id;
    };
    pls.forEach((p, k) => { const e = p.x.length - 1, a = nodeAt(p.x[0], p.y[0], p.z ? p.z[0] : null), b = p.closed ? a : nodeAt(p.x[e], p.y[e], p.z ? p.z[e] : null); edges.push({ from: a, to: b, type: 'pipe', length: (lineLength(p, geo) + (p.closed ? (geo ? haversine(p.x[e], p.y[e], p.x[0], p.y[0]) : Math.hypot(p.x[0] - p.x[e], p.y[0] - p.y[e], p.z ? p.z[0] - p.z[e] : 0)) : 0)) * (geo ? 1 : u), diameter: null, name: 'L' + (k + 1) }); });
  } else if (g.kind === 'table') {
    const H = g.headers || [], nrm = (h) => String(h).toLowerCase().replace(/[^a-z0-9]/g, ''), find = (...names) => H.find((h) => names.includes(nrm(h)));
    const hf = find('from', 'source', 'upstream', 'fromnode', 'start', 'inlet', 'node1'), ht = find('to', 'target', 'downstream', 'tonode', 'end', 'outlet', 'node2');
    if (!hf || !ht) throw new Error('The table needs from and to columns to describe a network.');
    const hl = find('length', 'l', 'len'), hd = find('diameter', 'd', 'dia', 'id', 'bore'), hty = find('type', 'kind', 'class'), hn = find('name', 'tag', 'label'), ids = new Map();
    const node = (id) => { id = String(id); if (!ids.has(id)) { ids.set(id, true); nodes.push({ id, type: 'node', x: null, y: null, z: null }); } return id; };
    for (const r of g.records || []) { if (r[hf] === '' || r[hf] === undefined || r[ht] === '' || r[ht] === undefined) continue; const a = node(r[hf]), b = node(r[ht]); edges.push({ from: a, to: b, type: hty ? String(r[hty] || 'pipe').toLowerCase() : 'pipe', length: hl ? num(numOf(r[hl])) : null, diameter: hd ? num(numOf(r[hd])) : null, name: hn && r[hn] !== '' ? String(r[hn]) : `${a}-${b}` }); }
  } else throw new Error(`A ${g.kind || 'unknown'} geometry does not describe a network.`);
  if (!nodes.length) throw new Error('The network holds no nodes.');
  // connectivity (union-find), duplicates, zero-length elements, independent loops
  const idx = new Map(nodes.map((n, i) => [n.id, i])), par = nodes.map((_, i) => i), root = (i) => { while (par[i] !== i) { par[i] = par[par[i]]; i = par[i]; } return i; };
  const pairs = new Map(), duplicates = [], zeroLength = [];
  let linked = 0;
  for (const e of edges) {
    const a = idx.get(e.from), b = idx.get(e.to);
    if (a === undefined || b === undefined) continue;
    linked++;
    if (a === b || e.length === 0) zeroLength.push(e.name);
    const k = a < b ? a + ',' + b : b + ',' + a;
    if (pairs.has(k)) duplicates.push(e.name); else pairs.set(k, 1);
    par[root(a)] = root(b);
  }
  const size = new Map();
  nodes.forEach((_, i) => size.set(root(i), (size.get(root(i)) || 0) + 1));
  let main = -1, best = 0;
  for (const [r, n] of size) if (n > best) { best = n; main = r; }
  return { nodes, edges, issues: { disconnected: nodes.filter((_, i) => root(i) !== main).map((n) => n.id), duplicates, zeroLength, loops: linked - nodes.length + size.size } };
}

// ---- Wall-thickness, corrosion and deposit maps ----------------------------------------------------------------------
/**
 * Map of a wall quantity over the pipe surface from a table (x | chainage, θ | clock position | y, value [, time]), a grid
 * (x = axial, y = circumferential) or points (x, θ, value).
 * Returns { x[] axial position (m when the unit is known), theta[] (deg; clock positions converted at 30° per hour),
 * t[theta][x] (NaN where there is no reading), kind: 'thickness' | 'deposit' | 'corrosion', min, max, mean, unit, missing, time, times }.
 * opts: kind (override), time (which time step of a time-resolved table; default the latest).
 */
export function mapFrom(g, opts = {}) {
  if (!g || typeof g !== 'object') throw new Error('mapFrom needs a geometry.');
  let x, theta, t, kind = opts.kind, unit = null, thetaUnit = 'deg', time = null, times = [];
  const lattice = (ax, an, val) => {
    const ux = [...new Set(ax)].sort((a, b) => a - b), ut = [...new Set(an)].sort((a, b) => a - b);
    if (ux.length * ut.length > 4e6) throw new Error('The map has too many distinct positions.');
    const xi = new Map(ux.map((v, k) => [v, k])), ti = new Map(ut.map((v, k) => [v, k])), m = ut.map(() => new Array(ux.length).fill(NaN));
    for (let i = 0; i < val.length; i++) m[ti.get(an[i])][xi.get(ax[i])] = val[i];
    x = ux; theta = ut; t = m;
  };
  if (g.kind === 'table') {
    const H = g.headers || [], tc = tableColumns(H), c = tc.c, ka = c.chainage ?? c.x, kt = c.theta ?? c.y, kv = c.deposit ?? c.thickness ?? c.loss;
    if (ka === undefined || kt === undefined || kv === undefined) throw new Error('The table needs an axial column (x or chainage), a circumferential column (theta, clock position or y) and a value column (thickness, corrosion depth or deposit).');
    kind = kind || (c.deposit !== undefined ? 'deposit' : c.thickness !== undefined ? 'thickness' : 'corrosion');
    let recs = g.records || [];
    if (c.time !== undefined) {
      const tv = recs.map((r) => (typeof r[H[c.time]] === 'number' ? r[H[c.time]] : Date.parse(r[H[c.time]])));
      times = [...new Set(tv.filter(fin))].sort((a, b) => a - b);
      if (times.length) { time = opts.time !== undefined && times.includes(+opts.time) ? +opts.time : times[times.length - 1]; recs = recs.filter((_, i) => tv[i] === time); }
    }
    const hn = String(H[kt]).toLowerCase(), clock = /clock/.test(hn), rad = (kt === c.theta ? tc.unit.theta : null) === 'rad', axU = ka === c.chainage ? tc.unit.chainage : tc.unit.x, ua = UNIT[axU] ?? (/^kp/.test(String(H[ka]).toLowerCase()) ? 1000 : 1);
    const ax = [], an = [], val = [];
    for (const r of recs) { const a = numOf(r[H[ka]]), b = numOf(r[H[kt]]), v = numOf(r[H[kv]]); if (fin(a) && fin(b) && fin(v)) { ax.push(a * ua); an.push(clock ? (b % 12) * 30 : rad ? b / D2R : b); val.push(v); } }
    if (!val.length) throw new Error('The map table holds no numeric rows.');
    lattice(ax, an, val);
    unit = (kv === c.deposit ? tc.unit.deposit : kv === c.thickness ? tc.unit.thickness : tc.unit.loss) || null;
    if (kt !== c.theta) thetaUnit = tc.unit.y || 'length';
  } else if (g.kind === 'grid') { x = Array.from(g.grid.x); theta = Array.from(g.grid.y); t = g.grid.z.map((r) => Array.from(r)); }
  else if (g.kind === 'points') {
    const n = Math.floor(g.points.length / 3), ax = [], an = [], val = [], sx = new Set(), sy = new Set();
    for (let i = 0; i < n; i++) { ax.push(g.points[3 * i]); an.push(g.points[3 * i + 1]); val.push(g.points[3 * i + 2]); if (sx.size <= 4096) sx.add(ax[i]); if (sy.size <= 4096) sy.add(an[i]); }
    if (sx.size * sy.size <= Math.max(4 * n, 64) && sx.size <= 4096 && sy.size <= 4096) lattice(ax, an, val);
    else { const gr = gridOf(g, 60, 36); x = Array.from(gr.x); theta = Array.from(gr.y); t = gr.z.map((r) => Array.from(r)); }
  } else throw new Error(`A ${g.kind || 'unknown'} geometry does not hold a wall map (use a table, grid or points with x, θ and a value).`);
  kind = ['thickness', 'deposit', 'corrosion'].includes(kind) ? kind : g.stats && g.stats.role === 'depositMap' ? 'deposit' : 'thickness';
  let min = Infinity, max = -Infinity, sum = 0, cnt = 0, missing = 0;
  for (const row of t) for (const v of row) { if (fin(v)) { if (v < min) min = v; if (v > max) max = v; sum += v; cnt++; } else missing++; }
  if (!cnt) throw new Error('The map holds no numeric values.');
  return { x, theta, t, kind, min, max, mean: sum / cnt, unit, thetaUnit, missing, time, times };
}

// ---- Terrain -----------------------------------------------------------------------------------------------------------
/**
 * Seabed / terrain grid for route planning: { x[], y[], elev[y][x], geographic, filled } from an elevation grid, soundings /
 * point cloud or a terrain surface mesh, at most n nodes per side; cells without data are filled from their neighbours.
 */
export function terrainFrom(g, n = 60) {
  if (!g || !['grid', 'points', 'mesh'].includes(g.kind)) throw new Error('terrainFrom needs an elevation grid, a point cloud or a surface mesh.');
  n = clamp(Math.round(+n) || 60, 2, 400);
  const small = g.kind === 'grid' && g.grid.x.length <= n && g.grid.y.length <= n, gr = small ? g.grid : gridOf(g, g.kind === 'grid' ? Math.min(n, g.grid.x.length) : n, g.kind === 'grid' ? Math.min(n, g.grid.y.length) : n);
  const x = Array.from(gr.x), y = Array.from(gr.y), elev = gr.z.map((r) => Array.from(r)), ny = y.length, nx = x.length;
  let holes = 0, filled = 0;
  for (const r of elev) for (const v of r) if (!fin(v)) holes++;
  if (holes === nx * ny) throw new Error('The terrain holds no valid elevations.');
  for (let pass = 0; holes > filled && pass < nx + ny; pass++) {
    const fill = [];
    for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {
      if (fin(elev[j][i])) continue;
      let s = 0, c = 0;
      for (const [a, b] of [[j - 1, i], [j + 1, i], [j, i - 1], [j, i + 1]]) if (a >= 0 && a < ny && b >= 0 && b < nx && fin(elev[a][b])) { s += elev[a][b]; c++; }
      if (c) fill.push([j, i, s / c]);
    }
    for (const [j, i, v] of fill) elev[j][i] = v;
    filled += fill.length;
    if (!fill.length) break;
  }
  return { x, y, elev, geographic: isGeo(g), filled, name: g.name };
}
