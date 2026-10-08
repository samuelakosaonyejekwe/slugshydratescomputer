// Routing of imported geometry onto suite inputs. Each rule looks at what was imported (pipeline route or
// profile, well survey, network, bathymetry, wall-thickness or deposit map, structural mesh, parameter
// table) and offers values only for inputs that the suite really has. Every offer names where the number came from.
import { sectionOf, gridOf } from './geom.js';
import { profileFrom, simplify, networkFrom, mapFrom, terrainFrom } from './route.js';

const UNIT = { mm: 1e-3, millimetre: 1e-3, cm: 1e-2, m: 1, metre: 1, km: 1000, in: 0.0254, inch: 0.0254, ft: 0.3048, um: 1e-6, micron: 1e-6 };
const toM = (g) => UNIT[String(g.stats?.units || 'm').toLowerCase()] ?? 1;
const fin = (x) => typeof x === 'number' && Number.isFinite(x);
const between = (x, a, b) => fin(x) && x >= a && x <= b;

/** Outlines of any geometry as an object the 2-D solvers accept. */
export function asOutline(g) {
  if (g.kind === 'mesh' || g.kind === 'polylines') return g;
  const segs = sectionOf(g) || [];
  if (!segs.length) return null;
  const lo = [Infinity, Infinity], hi = [-Infinity, -Infinity];
  for (const s of segs) { lo[0] = Math.min(lo[0], s[0], s[2]); hi[0] = Math.max(hi[0], s[0], s[2]); lo[1] = Math.min(lo[1], s[1], s[3]); hi[1] = Math.max(hi[1], s[1], s[3]); }
  return { kind: 'polylines', name: g.name, polylines: segs.map((s) => ({ x: [s[0], s[2]], y: [s[1], s[3]], closed: false })), bbox: { min: lo, max: hi } };
}

/** Elevation grid of a DEM, sounding cloud or terrain surface as { x | lon, y | lat, elev } (used by file inputs that take bathymetry). */
export function asBathy(g) {
  if (!['grid', 'points', 'mesh'].includes(g.kind)) return null;
  const gr = g.kind === 'grid' && g.grid.x.length <= 120 && g.grid.y.length <= 120 ? { x: g.grid.x, y: g.grid.y, z: g.grid.z } : gridOf(g, 60, 60);
  if (!gr) return null;
  const x = Array.from(gr.x), y = Array.from(gr.y), elev = gr.z.map((r) => Array.from(r, (v) => (Number.isFinite(v) ? v : 0)));
  const geo = g.kind === 'grid' && g.grid.geographic !== undefined ? g.grid.geographic : Math.abs(x[0]) <= 180 && Math.abs(y[0]) <= 90 && Math.abs(x[x.length - 1] - x[0]) < 5;
  return geo ? { lat: y, lon: x, elev, name: g.name } : { x, y, elev, name: g.name };
}

/** Name/value parameter tables and single-row tables whose columns match input keys. */
function fromTable(g, keys, fields) {
  if (g.kind !== 'table' || !g.records?.length) return [];
  const out = [], norm = (s) => String(s).toLowerCase().replace(/[^a-z0-9]/g, ''), byNorm = new Map();
  for (const f of fields) { if (f.type && f.type !== 'number') continue; byNorm.set(norm(f.key), f.key); byNorm.set(norm(f.label), f.key); }
  const hd = g.headers;
  if (hd.length >= 2 && g.records.every((r) => typeof r[hd[0]] === 'string')) { for (const r of g.records) { const k = byNorm.get(norm(r[hd[0]])); if (k && fin(r[hd[1]])) out.push({ key: k, value: r[hd[1]], from: `Parameter “${r[hd[0]]}” in ${g.name}` }); } }
  else for (const col of hd) { const k = byNorm.get(norm(col)); if (k && fin(g.records[0][col])) out.push({ key: k, value: g.records[0][col], from: `Column “${col}” in ${g.name}` }); }
  return out;
}

const tryIt = (fn) => { try { return fn(); } catch { return null; } };
const lazy = (d, key, fn) => { if (d && d[key] !== undefined) return d[key]; const v = tryIt(fn); if (d && typeof d === 'object') d[key] = v; return v; };
const isWell = (g) => !!(g.stats && g.stats.wellSurvey) || (g.kind === 'table' && g.stats && g.stats.role === 'survey');
/** Does this geometry describe one pipe route (as opposed to terrain, a wall map, a cloud or a solid)? */
function routeLike(g) {
  if (isWell(g)) return false;
  if (g.kind === 'polylines') return g.polylines.length > 0 && !g.polylines.every((p) => p.closed);
  if (g.kind === 'network') return true;
  if (g.kind === 'points') return !!(g.stats && g.stats.role === 'route3d');
  if (g.kind === 'table') return !!(g.stats && (g.stats.role === 'profile' || g.stats.role === 'route3d'));
  return false;
}
/** Pipeline profile of a route-like geometry (cached on the derived object). */
export const routeProfile = (g, d = {}) => (routeLike(g) ? lazy(d, 'profile', () => profileFrom(g)) : null);
/** Profile as the table value of a `profile` input: at most 60 points [{ x, z }] in metres. */
export function profileTable(p, max = 60) {
  let q = p;
  if (q.x.length > max) { let lo = 0, hi = Math.max(q.zMax - q.zMin, q.horizontalLength, 1); for (let it = 0; it < 40; it++) { const mid = 0.5 * (lo + hi); if (simplify(p, mid).x.length > max) lo = mid; else hi = mid; } q = simplify(p, hi); }
  return q.x.map((x, i) => ({ x: +x.toFixed(3), z: +q.z[i].toFixed(3) }));
}
/** Height of the final climb of a profile: the run of elements steeper than 20° that ends at the last node. */
function riserHeightOf(p) {
  let i = p.incl.length;
  while (i > 0 && p.incl[i - 1] >= 20) i--;
  return p.z[p.z.length - 1] - p.z[i];
}
const median = (a) => { const s = a.slice().sort((x, y) => x - y); return s.length ? s[Math.floor(s.length / 2)] : NaN; };
const netOf = (g, d) => (g.kind === 'network' || (g.kind === 'polylines' && g.polylines.length > 1 && !isWell(g)) ? lazy(d, 'network', () => networkFrom(g)) : null);
const mapOf = (g, d) => (g.kind === 'table' && g.stats && (g.stats.role === 'thicknessMap' || g.stats.role === 'depositMap') ? lazy(d, 'map', () => mapFrom(g)) : null);
const toMM = (m) => ({ m: 1000, cm: 10, in: 25.4, inch: 25.4, um: 1e-3, 'µm': 1e-3 })[m.unit] ?? 1;

/** Route length, depth, riser height and bore: offered wherever a suite declares those inputs. */
function pipeOffers(g, d, has, withProfile) {
  const out = [], p = routeProfile(g, d);
  if (p) {
    const what = `${g.format || g.kind}: ${p.source}`;
    if (withProfile && has('profile')) out.push({ key: 'profile', value: profileTable(p), from: `${what}, ${p.x.length} points${p.x.length > 60 ? ' simplified to at most 60' : ''}` });
    if (has('length') && between(p.length, 0.5, 2e6)) out.push({ key: 'length', value: p.length, from: `Length along the pipe of the imported route (${what})` });
    if (has('waterDepth') && between(-p.zMin, 1, 12000)) out.push({ key: 'waterDepth', value: -p.zMin, from: 'Deepest point of the imported route below datum' });
    const rh = riserHeightOf(p);
    if (has('riserHeight') && between(rh, 5, 12000)) out.push({ key: 'riserHeight', value: rh, from: 'Height of the final climb of the route (elements steeper than 20°)' });
  }
  if (g.kind === 'network' && has('id')) { const dia = g.network.edges.map((e) => e.diameter).filter((x) => between(x, 0.01, 5)); if (dia.length) out.push({ key: 'id', value: median(dia) * 1000, from: 'Median pipe diameter in the network (taken as metres, converted to mm)' }); }
  return out;
}
function networkOffer(g, d, has) {
  const n = has('network') ? netOf(g, d) : null;
  if (!n || !n.edges.length) return [];
  const bad = n.issues.disconnected.length + n.issues.zeroLength.length;
  return [{ key: 'network', value: n.edges.slice(0, 2000).map((e) => ({ from: e.from, to: e.to, type: e.type, length: e.length, diameter: e.diameter })), from: `${g.format || g.kind}: ${n.edges.length} connections between ${n.nodes.length} nodes${bad ? ` (${n.issues.disconnected.length} disconnected nodes, ${n.issues.zeroLength.length} zero-length connections: check before running)` : ''}` }];
}

const RULES = {
  net(g, d, has) {
    const out = pipeOffers(g, d, has, true);
    if (g.survey && has('survey')) { const s = g.survey, step = Math.max(1, Math.ceil(s.md.length / 400)), rows = []; for (let i = 0; i < s.md.length; i += step) rows.push({ md: s.md[i], inc: s.inc[i], azi: s.azi[i] }); if ((s.md.length - 1) % step) { const e = s.md.length - 1; rows.push({ md: s.md[e], inc: s.inc[e], azi: s.azi[e] }); } out.push({ key: 'survey', value: rows, from: `${g.format}: ${s.md.length} survey stations${step > 1 ? `, every ${step}th kept` : ''}` }); }
    if (isWell(g) && g.stats) {
      if (has('wellDepth') && between(g.stats.tvd, 1, 15000)) out.push({ key: 'wellDepth', value: g.stats.tvd, from: 'Deepest true vertical depth of the survey' });
      if (has('wellMD') && between(g.stats.md, 1, 20000)) out.push({ key: 'wellMD', value: g.stats.md, from: 'Measured depth of the last survey station' });
    }
    out.push(...networkOffer(g, d, has));
    if (has('terrain') && (g.kind === 'grid' || (g.kind === 'points' && !routeLike(g) && g.count >= 16 && !(g.stats && g.stats.role === 'particles')) || (g.kind === 'mesh' && (g.pathway === 'gis' || (g.stats && g.stats.surfaceOnly && g.bbox.max[2] - g.bbox.min[2] < 0.5 * Math.min(g.bbox.max[0] - g.bbox.min[0], g.bbox.max[1] - g.bbox.min[1])))))) {
      const t = lazy(d, 'terrain', () => terrainFrom(g, 60));
      if (t) { out.push({ key: 'terrain', value: t, from: `${g.format}: seabed / terrain grid ${t.x.length} × ${t.y.length}${t.geographic ? ' (geographic)' : ''}${t.filled ? `, ${t.filled} empty cells filled from neighbours` : ''}` }); const lo = Math.min(...t.elev.map((r) => Math.min(...r))); if (has('waterDepth') && between(-lo, 1, 12000) && !out.some((it) => it.key === 'waterDepth')) out.push({ key: 'waterDepth', value: -lo, from: 'Deepest point of the imported bathymetry' }); }
    }
    return out;
  },
  flow(g, d, has) { return pipeOffers(g, d, has, true); },
  solids(g, d, has) {
    const out = pipeOffers(g, d, has, false), m = mapOf(g, d);
    if (m && m.kind === 'deposit') {
      if (has('depositMap')) out.push({ key: 'depositMap', value: { ...m, name: g.name }, from: `${g.format}: deposit thickness over ${m.x.length} axial × ${m.theta.length} circumferential positions${m.times.length > 1 ? `, latest of ${m.times.length} time steps` : ''}` });
      if (has('depositMax') && between(m.max * toMM(m), 0, 5000)) out.push({ key: 'depositMax', value: m.max * toMM(m), from: `Largest deposit thickness in the imported map${m.unit ? '' : ' (values taken as mm)'}` });
    }
    return out;
  },
  ops(g, d, has) { return [...networkOffer(g, d, has), ...pipeOffers(g, d, has, false)]; },
  integ(g, d, has) {
    const out = pipeOffers(g, d, has, false), m = mapOf(g, d);
    if (m && m.kind !== 'deposit') {
      if (has('wtMap')) out.push({ key: 'wtMap', value: { ...m, name: g.name }, from: `${g.format}: ${m.kind === 'corrosion' ? 'metal-loss depth' : 'wall thickness'} over ${m.x.length} axial × ${m.theta.length} circumferential positions${m.missing ? `, ${m.missing} positions without a reading` : ''}` });
      if (m.kind === 'thickness' && has('minWt') && between(m.min * toMM(m), 0.1, 300)) out.push({ key: 'minWt', value: m.min * toMM(m), from: `Smallest wall thickness in the imported map${m.unit ? '' : ' (values taken as mm)'}` });
    }
    if (g.kind === 'mesh' && has('structMesh')) out.push({ key: 'structMesh', value: g, from: `${g.format}: ${g.count} surface triangles${g.stats && g.stats.cells ? ` from ${g.stats.cells} cells` : ''}` });
    if (has('spanLength') && g.kind === 'polylines' && g.polylines.length === 1 && !g.polylines[0].closed && !isWell(g)) { const p = routeProfile(g, d); if (p && between(p.length, 2, 300) && Math.abs(p.z[p.z.length - 1] - p.z[0]) <= 0.1 * p.length) out.push({ key: 'spanLength', value: p.horizontalLength, from: 'Horizontal length of the imported single span (ends at nearly the same level)' }); }
    return out;
  },
  econ(g, d, has) {
    const out = pipeOffers(g, d, has, false);
    if (isWell(g) && g.stats && has('wellDepth') && between(g.stats.tvd, 1, 15000)) out.push({ key: 'wellDepth', value: g.stats.tvd, from: 'Deepest true vertical depth of the survey' });
    return out;
  },
};

/** Offers for one suite. Returns [{ key, value, from }] restricted to inputs the suite declares. */
export function geometryLinks(suite, g, derived) {
  const fields = suite.inputs.flatMap((grp) => grp.fields), keys = new Set(fields.map((f) => f.key)), has = (k) => keys.has(k);
  const out = [...(RULES[suite.id]?.(g, derived || {}, has) || []), ...fromTable(g, keys, fields)];
  const custom = suite.geometry?.(g, derived) || [];
  const seen = new Set();
  return [...custom, ...out].filter((it) => it && keys.has(it.key) && it.value !== null && it.value !== undefined && !seen.has(it.key) && seen.add(it.key));
}
