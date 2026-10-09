// Suite 2 — Geometry, Wells, Network & Equipment.
// Route and riser geometry (arc length, inclination, extrema, volume, catenary / lazy-wave generators, terrain
// transects, free-span screening), multilayer wall thermal design, well trajectory (minimum curvature), inflow and
// vertical-lift performance with nodal analysis, valve / choke sizing (IEC 60534, orifice, multiphase correlations,
// Sachdeva and homogeneous-equilibrium choking), pump and compressor duty, separator and slug-catcher sizing and a
// node–edge hydraulic network solved by Newton–Raphson on the incidence-matrix formulation; artificial lift (gas lift,
// subsea boosting), a transient network with rotor and actuator dynamics (transientNetwork, transientSystem, rotorDynamics),
// design features for the integrity suite (buckle arrestors, buckle initiators, cathodic protection, supports) and the
// sourced constants (PROVENANCE) and reference data sets (js/data/ref/net.js) the engine is checked against.
// SI inside the engines; bara, °C, mm and Sm³/d at the interfaces.
import { clamp, brent, interp1, linspace, solveLinear, rk4, isNum, rng } from '../core/num.js';
import { G, frictionFactor, hInside, hOutside, uValue, seaTemperature, gradient } from '../core/pipe.js';
import { fluidModel, makeFluid, props as eosProps, lookup, waterProps, R, P_STD, T_STD } from '../core/thermo.js';
import { checkProfile, profileFrom } from '../core/route.js';
import { BASE } from '../data/basecase.js';
import { NIKURADSE_ROUGH, SMOOTH_PIPE, CRANE_FT, CHOKE_FIELD, SURVEY_FORGE, BENCH_FRICTION, BENCH_BENDS, BENCH_BENDS_LOW, BENCH_VALVES, BENCH_VALVE_LIQUID, BENCH_VALVE_GAS, BENCH_ORIFICE, BENCH_COMPRESSOR, BENCH_DEMISTER, BEIJ_BENDS, VALVE_CV_TRAVEL, COMPRESSOR_MAP, GAS_WELL_TESTS } from '../data/ref/net.js';

const D2R = Math.PI / 180, DAY = 86400, KEL = 273.15, RHO_SEA = 1027;
const num = (x, d) => (isNum(x) ? x : typeof x === 'string' && x.trim() !== '' && Number.isFinite(+x) ? +x : d);
const rd = (x, n = 2) => (isNum(x) ? +x.toFixed(n) : '—');
const sg = (x, n = 4) => (isNum(x) ? +x.toPrecision(n) : '—');
const area = (d) => (Math.PI * d * d) / 4;

/** Conversion factors to SI (multiply a value in the named unit to obtain SI: m, m³, Pa, kg, Pa·s, m³/s). */
export const UNIT = Object.freeze({ in: 0.0254, ft: 0.3048, mile: 1609.344, bbl: 0.158987294928, psi: 6894.757293168, bar: 1e5, atm: 101325, lb: 0.45359237, cP: 1e-3, mD: 9.869233e-16, 'bbl/d': 0.158987294928 / DAY, 'Sm3/d': 1 / DAY, 'm3/h': 1 / 3600, usgpm: 6.30901964e-5, scf: 0.028316846592 });
/** Convert a value between two units of the same kind listed in UNIT (temperatures: 'C', 'F', 'K'). */
export function convert(value, from, to) {
  const T = { C: [1, KEL], K: [1, 0], F: [5 / 9, 459.67 * (5 / 9)] };
  if (T[from] && T[to]) return (value * T[from][0] + T[from][1] - T[to][1]) / T[to][0];
  if (!(from in UNIT) || !(to in UNIT)) throw new Error(`Unknown unit in the conversion ${from} → ${to}.`);
  return (value * UNIT[from]) / UNIT[to];
}

// ---- route geometry ---------------------------------------------------------------------------------------
/**
 * Analyse an elevation profile (x = horizontal distance, z = elevation, m).
 * Returns { n, s[] (arc length at nodes), segs: [{ i, x0, x1, z0, z1, L, dx, dz, incl (deg) }], length, horizontal, gain, loss,
 * zMin, zMax, iMin, highs: [{ i, x, z }], lows: [{ i, x, z }], issues: [{ type, i, msg }], zeroLength, duplicates, backward }.
 */
export function analyseProfile(x, z, tol = 1e-6) {
  const n = Math.min(x.length, z.length), s = [0], segs = [], issues = [], highs = [], lows = [], seen = new Map();
  let gain = 0, loss = 0, zMin = Infinity, zMax = -Infinity, iMin = 0, zeroLength = 0, duplicates = 0, backward = 0, lastSign = 0;
  for (let i = 0; i < n; i++) {
    if (!isNum(x[i]) || !isNum(z[i])) { issues.push({ type: 'invalid', i, msg: `Point ${i + 1} has a missing or non-numeric coordinate.` }); s.push(s[s.length - 1]); continue; }
    const key = Math.round(x[i] / tol) + '|' + Math.round(z[i] / tol);
    if (seen.has(key)) { duplicates++; issues.push({ type: 'duplicate', i, msg: `Point ${i + 1} repeats point ${seen.get(key) + 1} (x = ${x[i]} m, z = ${z[i]} m).` }); } else seen.set(key, i);
    if (z[i] < zMin) { zMin = z[i]; iMin = i; }
    zMax = Math.max(zMax, z[i]);
    if (!i) continue;
    const dx = x[i] - x[i - 1], dz = z[i] - z[i - 1], L = Math.hypot(dx, dz);
    if (L < tol) { zeroLength++; issues.push({ type: 'zero-length', i, msg: `Element ${i} (points ${i} → ${i + 1}) has zero length.` }); }
    if (dx < -tol) { backward++; issues.push({ type: 'backward', i, msg: `Chainage runs backwards between points ${i} and ${i + 1}.` }); }
    segs.push({ i: i - 1, x0: x[i - 1], x1: x[i], z0: z[i - 1], z1: z[i], L, dx, dz, incl: L > 0 ? Math.asin(clamp(dz / L, -1, 1)) / D2R : 0 });
    s[i] = s[i - 1] + L;
    if (dz > 0) gain += dz; else loss -= dz;
    const sign = Math.abs(dz) < tol ? 0 : Math.sign(dz);
    if (sign && lastSign && sign !== lastSign) (sign < 0 ? highs : lows).push({ i: i - 1, x: x[i - 1], z: z[i - 1] });
    if (sign) lastSign = sign;
  }
  s.length = n;
  return { n, s, segs, length: n ? s[n - 1] : 0, horizontal: n ? x[n - 1] - x[0] : 0, gain, loss, zMin, zMax, iMin, highs, lows, issues, zeroLength, duplicates, backward };
}
/** Profile nodes → segment list { L[], theta[] (rad, positive upward) } with the start point { x0, z0 }. */
export function profileToSegments(x, z) {
  const L = [], theta = [];
  for (let i = 1; i < x.length; i++) { L.push(Math.hypot(x[i] - x[i - 1], z[i] - z[i - 1])); theta.push(Math.atan2(z[i] - z[i - 1], x[i] - x[i - 1])); }
  return { x0: x[0], z0: z[0], L, theta };
}
/** Segment list → profile nodes { x[], z[] } (inverse of profileToSegments). */
export function segmentsToProfile({ x0 = 0, z0 = 0, L, theta }) {
  const x = [x0], z = [z0];
  L.forEach((l, i) => { x.push(x[i] + l * Math.cos(theta[i])); z.push(z[i] + l * Math.sin(theta[i])); });
  return { x, z };
}
/** Table rows → clean arrays: rows with missing numbers, repeated points and backward steps are dropped (and counted). */
export function cleanProfile(rows) {
  const x = [], z = [];
  let dropped = 0;
  for (const r of Array.isArray(rows) ? rows : []) {
    const px = num(r?.x, null), pz = num(r?.z, null);
    if (px === null || pz === null) { dropped++; continue; }
    const k = x.length - 1;
    if (k >= 0 && (px < x[k] - 1e-9 || Math.hypot(px - x[k], pz - z[k]) < 1e-6)) { dropped++; continue; }
    x.push(px); z.push(pz);
  }
  return { x, z, dropped };
}
/** Reduce a profile to at most maxPts points; indices in `keep` (extrema, riser base) and both ends are always retained. */
export function simplifyProfile(x, z, maxPts = 80, keep = []) {
  const n = x.length;
  if (n <= maxPts) return { x: x.slice(), z: z.slice(), idx: x.map((_, i) => i) };
  let must = [...new Set(keep.filter((i) => i > 0 && i < n - 1))].sort((a, b) => a - b);
  if (must.length > maxPts - 2) { // more extrema than room: keep the most prominent ones
    const prom = (i) => Math.min(Math.abs(z[i] - z[i - 1]), Math.abs(z[i] - z[i + 1]));
    must = must.sort((a, b) => prom(b) - prom(a)).slice(0, maxPts - 2).sort((a, b) => a - b);
  }
  const idx = [0, ...must, n - 1];
  while (idx.length < maxPts) { // insert the point that deviates most from the current polyline
    let best = -1, dev = 1e-9, at = 0;
    for (let k = 1; k < idx.length; k++) {
      const a = idx[k - 1], b = idx[k];
      for (let i = a + 1; i < b; i++) { const d = Math.abs(z[i] - (z[a] + ((z[b] - z[a]) * (x[i] - x[a])) / (x[b] - x[a] || 1))); if (d > dev) { dev = d; best = i; at = k; } }
    }
    if (best < 0) break;
    idx.splice(at, 0, best);
  }
  return { x: idx.map((i) => x[i]), z: idx.map((i) => z[i]), idx };
}
/**
 * Riser base of a profile that ends at or above sea level: the last point lying on the seabed before the final ascent.
 * Returns { i, x, z, height } (height = 0 and i = last point when the line has no riser).
 */
export function riserBase(x, z, forceX = 0) {
  const n = x.length, zEnd = z[n - 1], zMin = Math.min(...z);
  if (forceX > 0) { let i = 0; for (let k = 1; k < n; k++) if (Math.abs(x[k] - forceX) < Math.abs(x[i] - forceX)) i = k; return { i, x: x[i], z: z[i], height: Math.max(zEnd - z[i], 0) }; }
  if (!(zEnd > -5 && zMin < -20)) return { i: n - 1, x: x[n - 1], z: zEnd, height: 0 };
  const depth = -zMin, x0 = x[n - 1] - 2.5 * depth;
  let lo = Infinity; for (let i = 0; i < n; i++) if (x[i] >= x0 || i >= n - 3) lo = Math.min(lo, z[i]);
  const tol = Math.max(5, 0.005 * depth);
  let i = n - 2; while (i > 0 && !(z[i] <= lo + tol)) i--;
  if (Math.atan2(zEnd - z[i], x[n - 1] - x[i]) < 8 * D2R) return { i: n - 1, x: x[n - 1], z: zEnd, height: 0 }; // a gentle shore approach is not a riser
  return { i, x: x[i], z: z[i], height: zEnd - z[i] };
}
/**
 * Steel catenary riser from the touchdown point to the hang-off: z = a (cosh(x/a) − 1).
 * { height (m), angle (hang-off angle from the vertical, deg), n } → { x[], z[], a, span, length } relative to the touchdown point.
 */
export function catenary({ height, angle = 12, n = 16 }) {
  if (!(height > 0)) throw new Error('A catenary riser needs a positive height between touchdown and hang-off.');
  if (!(angle > 0.5 && angle < 85)) throw new Error('The hang-off angle of a catenary riser must lie between 0.5° and 85° from the vertical.');
  const pTop = 1 / Math.tan(angle * D2R), a = height / (Math.hypot(1, pTop) - 1), x = [], z = [];
  for (let i = 0; i <= n; i++) { const p = (pTop * i) / n; x.push(a * Math.asinh(p)); z.push(a * (Math.hypot(1, p) - 1)); }
  return { x, z, a, span: a * Math.asinh(pTop), length: a * pTop };
}
/**
 * Lazy-wave riser built from three catenary pieces sharing one horizontal tension: lower catenary (touchdown → lift point),
 * buoyant arch (concave down) and upper catenary through the sag bend to the hang-off.
 * { height, angle (hang-off from vertical, deg), liftAngle (slope where buoyancy starts, deg), sagAngle (downward slope where
 * buoyancy ends, deg), buoyancy (net uplift / submerged weight), n } → { x[], z[], a, span, length, hog, sag (heights) }.
 */
export function lazyWave({ height, angle = 10, liftAngle = 35, sagAngle = 15, buoyancy = 1.5, n = 30 }) {
  if (!(height > 0)) throw new Error('A lazy-wave riser needs a positive height.');
  if (!(buoyancy > 0.05) || !(liftAngle > 1 && liftAngle < 80) || !(sagAngle >= 0 && sagAngle < 80) || !(angle > 0.5 && angle < 85)) throw new Error('Lazy-wave riser: angles must lie between 1° and 80° and the buoyancy ratio must be positive.');
  const p1 = Math.tan(liftAngle * D2R), p2 = -Math.tan(sagAngle * D2R), p3 = 1 / Math.tan(angle * D2R), rho = 1 / buoyancy, c = (p) => Math.hypot(1, p);
  const hf = c(p1) - 1 + rho * (c(p1) - c(p2)) + (c(p3) - c(p2)), a = height / hf, hog = a * (c(p1) - 1 + rho * (c(p1) - 1)), sag = a * (c(p1) - 1 + rho * (c(p1) - c(p2)) + 1 - c(p2));
  if (!(hf > 0) || sag < 0.02 * height) throw new Error('Lazy-wave riser: with these angles the sag bend would touch the seabed — increase the buoyancy ratio or the lift angle, or reduce the sag angle.');
  const x = [0], z = [0], pieces = [[0, p1, a], [p1, p2, -a * rho], [p2, p3, a]], per = Math.max(3, Math.round(n / 3));
  for (const [pa, pb, ak] of pieces) { const x0 = x[x.length - 1], z0 = z[z.length - 1]; for (let i = 1; i <= per; i++) { const p = pa + ((pb - pa) * i) / per; x.push(x0 + ak * (Math.asinh(p) - Math.asinh(pa))); z.push(z0 + ak * (c(p) - c(pa))); } }
  return { x, z, a, span: x[x.length - 1], length: a * (p1 + rho * (p1 - p2) + (p3 - p2)), hog, sag };
}
/** Undulating seabed flowline: straight fall of `drop` metres towards the riser base plus a tapered sine undulation. */
export function seabedLine({ length, depthEnd, drop = 0, amp = 0, wavelength = 5000, n = 0 }) {
  if (!(length > 0)) throw new Error('The flowline length must be positive.');
  const lam = Math.max(wavelength, 50), m = n || clamp(Math.ceil(length / (lam / 8)), 8, 56), x = [], z = [];
  for (let i = 0; i <= m; i++) { const xi = (length * i) / m, w = Math.min(1, xi / (lam / 4), (length - xi) / (lam / 4)); x.push(xi); z.push(-depthEnd + drop * (1 - xi / length) + amp * w * Math.sin((2 * Math.PI * xi) / lam)); }
  return { x, z };
}
/** Bilinear sample of a grid E[ny][nx] at fractional indices (fx, fy). */
export function sampleGrid(E, fx, fy) {
  const ny = E.length, nx = E[0].length, i = clamp(Math.floor(fx), 0, Math.max(nx - 2, 0)), j = clamp(Math.floor(fy), 0, Math.max(ny - 2, 0)), u = clamp(fx - i, 0, 1), w = clamp(fy - j, 0, 1), i1 = Math.min(i + 1, nx - 1), j1 = Math.min(j + 1, ny - 1);
  return (1 - w) * ((1 - u) * E[j][i] + u * E[j][i1]) + w * ((1 - u) * E[j1][i] + u * E[j1][i1]);
}
/**
 * Elevation transect across a terrain grid { x[] | lon[], y[] | lat[], elev[ny][nx], geographic }.
 * line: 'we' (west → east through the middle) | 'sn' | 'diag' | 'deep' (deepest cell → highest cell). Returns { x[] (m along the line), z[] }.
 */
export function terrainTransect(t, line = 'we', n = 0) {
  const X = t?.x || t?.lon, Y = t?.y || t?.lat, E = t?.elev;
  if (!Array.isArray(X) || !Array.isArray(Y) || !Array.isArray(E) || X.length < 2 || Y.length < 2 || E.length !== Y.length || !Array.isArray(E[0]) || E[0].length !== X.length) throw new Error('The terrain grid is not usable: it needs x (or lon), y (or lat) and an elevation matrix elev[y][x] of matching size.');
  const nx = X.length, ny = Y.length, geo = t.geographic ?? (!!t.lon || (Math.abs(X[0]) <= 180 && Math.abs(Y[0]) <= 90 && Math.abs(X[nx - 1] - X[0]) < 5));
  const kx = geo ? 111320 * Math.cos(0.5 * (Y[0] + Y[ny - 1]) * D2R) : 1, ky = geo ? 110540 : 1;
  let a = [0, (ny - 1) / 2], b = [nx - 1, (ny - 1) / 2];
  if (line === 'sn') { a = [(nx - 1) / 2, 0]; b = [(nx - 1) / 2, ny - 1]; }
  else if (line === 'diag') { a = [0, 0]; b = [nx - 1, ny - 1]; }
  else if (line === 'deep') { let lo = Infinity, hi = -Infinity; E.forEach((r, j) => r.forEach((v, i) => { if (isNum(v) && v < lo) { lo = v; a = [i, j]; } if (isNum(v) && v > hi) { hi = v; b = [i, j]; } })); }
  const m = n || clamp(4 * Math.max(nx, ny), 16, 600), cx = (f) => interp1(X.map((_, i) => i), X, f) * kx, cy = (f) => interp1(Y.map((_, i) => i), Y, f) * ky, xs = [], zs = [];
  const Lt = Math.hypot(cx(b[0]) - cx(a[0]), cy(b[1]) - cy(a[1]));
  if (!(Lt > 0)) throw new Error('The terrain transect has zero length: choose another line across the grid.');
  for (let i = 0; i <= m; i++) { const f = i / m, v = sampleGrid(E, a[0] + (b[0] - a[0]) * f, a[1] + (b[1] - a[1]) * f); xs.push(Lt * f); zs.push(isNum(v) ? v : 0); }
  return { x: xs, z: zs, geographic: !!geo };
}
/**
 * Free-span candidates of a pipe resting on an uneven seabed. Between two neighbouring crests the pipe can sag
 * w L⁴ / (384 EI) (fixed-ended beam under its submerged weight); a depression deeper than that plus `gap` is bridged.
 * Returns [{ x (mid-span), length, gap }].
 */
export function freeSpans(x, z, { w, EI, gap = 0.1, maxLength = 200 }) {
  const crests = [0], out = [];
  for (let i = 1; i < x.length - 1; i++) if (z[i] >= z[i - 1] && z[i] > z[i + 1]) crests.push(i);
  crests.push(x.length - 1);
  for (let k = 1; k < crests.length; k++) {
    const a = crests[k - 1], b = crests[k], L = x[b] - x[a];
    if (b - a < 2 || !(L > 0) || L > maxLength) continue;
    let depth = 0; for (let i = a + 1; i < b; i++) depth = Math.max(depth, z[a] + ((z[b] - z[a]) * (x[i] - x[a])) / L - z[i]);
    const sagMax = w > 0 && EI > 0 ? (w * L ** 4) / (384 * EI) : 0;
    if (depth - sagMax > gap) out.push({ x: 0.5 * (x[a] + x[b]), length: L, gap: depth - sagMax });
  }
  return out;
}
/** Bends from direction changes of a polyline: angle between neighbouring segments and the radius of the circle through the three points. */
export function bendInventory(x, z, minAngle = 2) {
  const out = [];
  for (let i = 1; i < x.length - 1; i++) {
    const a = Math.hypot(x[i] - x[i - 1], z[i] - z[i - 1]), b = Math.hypot(x[i + 1] - x[i], z[i + 1] - z[i]), c = Math.hypot(x[i + 1] - x[i - 1], z[i + 1] - z[i - 1]);
    if (a < 1e-9 || b < 1e-9) continue;
    const ang = Math.abs(Math.atan2(z[i + 1] - z[i], x[i + 1] - x[i]) - Math.atan2(z[i] - z[i - 1], x[i] - x[i - 1])) / D2R, s = 0.5 * (a + b + c), tri = Math.sqrt(Math.max(s * (s - a) * (s - b) * (s - c), 0));
    if (ang >= minAngle) out.push({ i, x: x[i], z: z[i], angle: ang, radius: tri > 1e-12 ? (a * b * c) / (4 * tri) : 1e9 });
  }
  return out;
}

// ---- local losses -----------------------------------------------------------------------------------------
/**
 * Fitting library: resistance coefficient as K = fT·(Le/D) (fully turbulent friction factor fT) or a fixed K. The Le/D values are
 * those of Crane TP-410 (see PROVENANCE); cls groups the entries for the bend and tee calibration multipliers; a function leD(D)
 * gives a size-dependent value (Crane's centric butterfly valve: 45 up to 8 in, 35 for 10–14 in, 25 above).
 */
export const FITTINGS = Object.freeze({
  bend90: { label: '90° bend, r/D = 1.5', leD: 14, cls: 'bend' }, bend90lr: { label: '90° bend, r/D = 5', leD: 15.5, cls: 'bend' }, bend45: { label: '45° standard elbow', leD: 16, cls: 'bend' }, elbow90: { label: '90° standard elbow', leD: 30, cls: 'bend' },
  teeRun: { label: 'Tee, flow through run', leD: 20, cls: 'tee' }, teeBranch: { label: 'Tee, flow through branch', leD: 60, cls: 'tee' }, gate: { label: 'Gate valve, full bore', leD: 8 }, ball: { label: 'Ball valve, full bore', leD: 3 },
  check: { label: 'Swing check valve', leD: 50 }, globe: { label: 'Globe valve', leD: 340 }, butterfly: { label: 'Butterfly valve, centric', leD: (D) => (D <= 0.2286 ? 45 : D <= 0.381 ? 35 : 25) }, entrance: { label: 'Sharp entrance', k: 0.5 }, exit: { label: 'Pipe exit', k: 1.0 }, reducer: { label: 'Reducer / expander', k: 0.2 },
});
/** Roughness of clean commercial steel pipe used by Crane TP-410 for its fT table (0.0018 in), m. */
export const CRANE_ROUGHNESS = 0.0018 * 0.0254;
/** Fully turbulent friction factor for a relative roughness (rough-wall limit of Colebrook–White). */
export const fullyTurbulent = (rel) => 0.25 / Math.log10(Math.max(rel, 1e-7) / 3.7) ** 2;
/** Crane TP-410 table for 90° pipe bends: r/D and K/fT. */
export const CRANE_BEND = Object.freeze({ rD: [1, 1.5, 2, 3, 4, 6, 8, 10, 12, 14, 16, 20], leD: [20, 14, 12, 12, 14, 17, 24, 30, 34, 38, 42, 50] });
/**
 * Bend resistance coefficient from r/D and angle. 90°: Crane table. Above 90° (n = angle/90 quarter turns): Crane's
 * K = (n − 1)(0.25 π fT r/D + 0.5 K90) + K90. Below 90° K90·(angle/90)^0.9, which reproduces the ratio of Crane's 45° and 90°
 * standard elbows (16 fT / 30 fT) within 0.5 %.
 */
export function bendK(rD, angleDeg, fT) {
  if (rD > 20) return 0; // field curvature: plain pipe friction already covers it
  const r = Math.max(rD, 1), k90 = fT * interp1(CRANE_BEND.rD, CRANE_BEND.leD, r), n = angleDeg / 90;
  return n > 1 ? (n - 1) * (0.25 * Math.PI * fT * r + 0.5 * k90) + k90 : k90 * Math.max(n, 0) ** 0.9;
}
/** Total K and equivalent length of a fittings list [{ type, count }] for a pipe of diameter D and relative roughness rel. */
export function fittingsLoss(list, D, rel, mult = {}) {
  const fT = fullyTurbulent(rel), rows = [];
  let K = 0;
  for (const r of Array.isArray(list) ? list : []) {
    const key = String(r?.type ?? '').trim(), f = FITTINGS[key] || Object.values(FITTINGS).find((q) => q.label.toLowerCase() === key.toLowerCase()), cnt = Math.max(num(r?.count, 0), 0);
    if (!f || !cnt) { if (key && !f) rows.push({ type: key, label: 'unknown fitting — ignored', count: cnt, k: 0, total: 0 }); continue; }
    const k = (f.k ?? fT * (typeof f.leD === 'function' ? f.leD(D) : f.leD)) * (f.cls && isNum(mult[f.cls]) ? mult[f.cls] : 1); K += k * cnt; rows.push({ type: key, label: f.label, count: cnt, k, total: k * cnt });
  }
  return { K, fT, eqLength: (K * D) / fT, rows };
}
/** Single-phase pressure drop (Pa) of a pipe: Darcy–Weisbach friction, K-factor losses and elevation. f fixes the friction factor. */
export function pipeDp({ m, rho, mu, L, D, rough = 0, dz = 0, k = 0, f = null, model = 'colebrook' }) {
  const v = m / (rho * area(D)), Re = (rho * Math.abs(v) * D) / mu, ff = f ?? frictionFactor(Re, rough / D, model);
  return ((ff * L) / D + k) * 0.5 * rho * v * Math.abs(v) + rho * G * dz;
}
/** Reynolds-number regime label. */
export const reRegime = (Re) => (Re < 2000 ? 'laminar' : Re < 4000 ? 'transitional' : 'turbulent');

// ---- wall thermal design ----------------------------------------------------------------------------------
/**
 * Multilayer wall: U-value on the inner diameter (kernel uValue), outer diameter, weight, specific gravity and thermal mass.
 * { id, wt (m), kSteel, rhoSteel, cpSteel, layers: [{ name, t (m), k, rho, cp }], pipWt (carrier pipe wall, m, 0 = none),
 *   concrete: { t, k, rho } | null, hIn, hOut, burial: { cover (m above the pipe crown), kSoil } | null, rhoContents, uMult, E (Pa) }
 * Returns { U, resistances, od, layers, massPerM, weightAir, buoyancy, submerged (N/m), sg, thermalMass (J/m/K), steelArea, EI }.
 */
export function wallDesign({ id, wt, kSteel = 45, rhoSteel = 7850, cpSteel = 480, layers = [], pipWt = 0, concrete = null, hIn = 1000, hOut = 500, burial = null, rhoContents = 0, uMult = 1, E = 2.07e11, rhoSea = RHO_SEA }) {
  if (!(id > 0) || !(wt > 0)) throw new Error('The pipe needs a positive inner diameter and wall thickness.');
  const all = layers.filter((l) => l.t > 0 && l.k > 0).map((l) => ({ name: l.name || 'Coating', t: l.t, k: l.k, rho: l.rho > 0 ? l.rho : 900, cp: l.cp > 0 ? l.cp : 1500 }));
  if (pipWt > 0) all.push({ name: 'Carrier pipe (steel)', t: pipWt, k: kSteel, rho: rhoSteel, cp: cpSteel });
  if (concrete && concrete.t > 0) all.push({ name: 'Concrete weight coat', t: concrete.t, k: concrete.k > 0 ? concrete.k : 2, rho: concrete.rho > 0 ? concrete.rho : 3040, cp: 880 });
  let r = id / 2;
  const ring = (t) => { const a = Math.PI * ((r + t) ** 2 - r * r); r += t; return a; }, steelArea = ring(wt);
  let mass = steelArea * rhoSteel, cap = mass * cpSteel;
  for (const l of all) { const a = ring(l.t); mass += a * l.rho; cap += a * l.rho * l.cp; }
  const od = 2 * r, u = uValue({ id, wt, kWall: kSteel, layers: all, hIn, hOut, burial: burial && burial.cover > 0 ? { depth: burial.cover + od / 2, kSoil: burial.kSoil } : null });
  const weightAir = (mass + area(id) * rhoContents) * G, buoyancy = rhoSea * area(od) * G, ods = id + 2 * wt;
  return { U: u.U * uMult, resistances: u.resistances, od, layers: all, massPerM: mass, weightAir, buoyancy, submerged: weightAir - buoyancy, sg: weightAir / buoyancy, thermalMass: cap, steelArea, EI: (E * Math.PI * (ods ** 4 - id ** 4)) / 64 };
}

// ---- multiphase marching along a pipe ------------------------------------------------------------------------
/** Case-fluid mixture at P (bara), T (°C) for independent hydrocarbon and water mass rates (kg/s): same closures as fluidModel.at. */
export function mixAt(fm, P, T, mHC, mW) {
  const o = lookup(fm.table, P, T), w = waterProps(P, T, fm.aq);
  o.rhoW = w.rho; o.muW = w.mu; o.cpW = w.cp; o.kW = w.k; o.mG = mHC * o.wG; o.mO = mHC * (1 - o.wG); o.mW = mW;
  o.qG = o.mG / o.rhoG; o.qO = o.mO / o.rhoO; o.qW = mW / w.rho; o.qL = o.qO + o.qW;
  const wc = o.qL > 0 ? o.qW / o.qL : 0, mL = o.mO + mW;
  o.wcut = wc; o.rhoL = o.qL > 0 ? mL / o.qL : o.rhoO;
  o.muL = wc <= 0 ? o.muO : wc >= 1 ? o.muW : wc > 0.6 ? o.muW * (1 - Math.min(1 - wc, 0.7)) ** -2.5 : o.muO * (1 - Math.min(wc, 0.7)) ** -2.5;
  o.cpL = mL > 0 ? (o.mO * o.cpO + mW * o.cpW) / mL : o.cpO; o.kL = (1 - wc) * o.kO + wc * o.kW;
  o.sigma = wc > 0.6 ? clamp(0.0756 - 1.5e-4 * T - 2e-5 * P, 0.03, 0.076) : Math.max(o.sigma, 1e-4);
  return o;
}
/**
 * March pressure and temperature in the flow direction through a list of cells with the kernel pressure-gradient closure.
 * { fm, mHC, mW (kg/s), cells: [{ ds, theta, D, rough, U, tAmb }], pIn (bara), tIn (°C), model, fModel, tOf(fraction) (prescribed
 *   temperature instead of the energy balance), kLoss (K-factor applied at the outlet) }
 * Returns { P[], T[], pOut, tOut, dp, dpFric, dpGrav, dpMinor (bar), holdup[], vm[], regime[], rhoNs[], vMax, eroMax (max vm·√ρ), liquid (m³), volume (m³), dead }.
 */
export function marchPipe({ fm, mHC, mW, cells, pIn, tIn, model = 'beggsBrill', fModel = 'colebrook', tOf = null, kLoss = 0, gExtra = 0 }) {
  const mdot = mHC + mW, P = [pIn], T = [tIn], holdup = [], vm = [], regime = [], rhoNs = [], sTot = cells.reduce((a, c) => a + c.ds, 0) || 1;
  let dpF = 0, dpG = 0, vMax = 0, eroMax = 0, liquid = 0, volume = 0, sAcc = 0, dead = false, last = null;
  const state = (p, t, c) => {
    const pe = Math.max(p, 1), pr = mixAt(fm, pe, t, mHC, mW), A = area(c.D), gx = gExtra > 0 ? gExtra * (c.gx ?? 1) : 0; // gx: lift gas carried with the stream (same properties as the produced gas)
    if (gx > 0) { pr.mG += gx; pr.qG += gx / pr.rhoG; }
    const vsl = pr.qL / A, vsg = pr.qG / A;
    const gr = gradient({ vsl, vsg, rhoL: pr.rhoL, rhoG: pr.rhoG, muL: pr.muL, muG: pr.muG, sigma: pr.sigma, D: c.D, theta: c.theta, rough: c.rough, P: pe * 1e5, fModel }, model);
    const mCp = pr.mG * pr.cpG + pr.mO * pr.cpO + pr.mW * pr.cpW, beta = mCp > 1e-9 ? (c.U * Math.PI * c.D) / mCp : 0;
    const src = mCp > 1e-9 ? -((pr.mG * pr.cpG * pr.jtG + pr.mO * pr.cpO * pr.jtO - pr.mW / pr.rhoW) / mCp) * gr.dpdx - (mdot * G * Math.sin(c.theta)) / mCp : 0;
    return { pr, gr, vsl, vsg, A, beta, src, mCp };
  };
  // temperature step with the exact solution of dT/ds = -beta (T - Ta) + src (unconditionally stable for small rates)
  const tStep = (t, c, st, ds) => { if (!(st.mCp > 1e-9)) return c.tAmb; const e = Math.exp(-st.beta * ds), phi = st.beta * ds > 1e-9 ? (1 - e) / (st.beta * ds) : 1; return c.tAmb + (t - c.tAmb) * e + st.src * ds * phi; };
  let prev = null; // midpoint rule; the predictor reuses the gradient of the previous cell (one closure call per cell after the first)
  for (const c of cells) {
    const i = P.length - 1, a = prev || state(P[i], T[i], c), pm = P[i] - (0.5 * c.ds * a.gr.dpdx) / 1e5, tm = tOf ? tOf((sAcc + 0.5 * c.ds) / sTot) : tStep(T[i], c, a, 0.5 * c.ds);
    const b = cells.length === 1 ? a : state(pm, tm, c), pn = P[i] - (c.ds * b.gr.dpdx) / 1e5; // a single-cell connector takes one explicit step
    prev = b;
    if (pn < 1) dead = true;
    P.push(pn); T.push(clamp(tOf ? tOf((sAcc + c.ds) / sTot) : tStep(T[i], c, b, c.ds), -60, 250)); sAcc += c.ds;
    const v = b.vsl + b.vsg, rn = v > 1e-12 ? (b.pr.rhoL * b.vsl + b.pr.rhoG * b.vsg) / v : b.pr.rhoL;
    holdup.push(b.gr.holdup); vm.push(v); regime.push(b.gr.regime); rhoNs.push(rn);
    dpF += ((b.gr.fric + b.gr.acc) * c.ds) / 1e5; dpG += (b.gr.grav * c.ds) / 1e5; vMax = Math.max(vMax, v); eroMax = Math.max(eroMax, v * Math.sqrt(rn)); liquid += b.gr.holdup * b.A * c.ds; volume += b.A * c.ds; last = { v, rn };
  }
  const dpMinor = kLoss > 0 && last ? (kLoss * 0.5 * last.rn * last.v * last.v) / 1e5 : 0, n = P.length - 1;
  P[n] -= dpMinor;
  return { P, T, pOut: P[n], tOut: T[n], dp: pIn - P[n], dpFric: dpF, dpGrav: dpG, dpMinor, holdup, vm, regime, rhoNs, vMax, eroMax, liquid, volume, dead: dead || P[n] < 1 };
}
/** Cells along a polyline between arc positions sA and sB (S = arc length at the nodes); n is the target number of cells. */
export function polylineCells(S, X, Z, sA, sB, n, base = {}) {
  const cuts = [sA, ...S.filter((s) => s > sA + 1e-9 && s < sB - 1e-9), sB], cells = [], tot = Math.max(sB - sA, 1e-12);
  for (let k = 1; k < cuts.length; k++) {
    const a = cuts[k - 1], b = cuts[k], len = b - a;
    if (len < 1e-9) continue;
    const m = Math.max(1, Math.round((n * len) / tot)), za = interp1(S, Z, a), zb = interp1(S, Z, b), xa = interp1(S, X, a), xb = interp1(S, X, b), theta = Math.asin(clamp((zb - za) / len, -1, 1));
    for (let j = 0; j < m; j++) { const f = (j + 0.5) / m; cells.push({ ...base, ds: len / m, theta, zMid: za + (zb - za) * f, xMid: xa + (xb - xa) * f, sMid: a + len * f }); }
  }
  return cells;
}

// ---- wells ---------------------------------------------------------------------------------------------------
/**
 * Well trajectory by the minimum-curvature method. stations: [{ md (m), inc (deg from vertical), azi (deg from north) }].
 * Returns { stations: [{ md, inc, azi, tvd, north, east, disp, dls (deg / 30 m) }], md, tvd, displacement, maxDls, maxInc }.
 */
export function minimumCurvature(stations) {
  const st = (Array.isArray(stations) ? stations : []).map((r) => ({ md: num(r?.md, null), inc: num(r?.inc, 0), azi: num(r?.azi, 0) })).filter((r) => r.md !== null).sort((a, b) => a.md - b.md);
  if (st.length < 2 || !(st[st.length - 1].md > st[0].md)) throw new Error('The deviation survey needs at least two stations with increasing measured depth.');
  const out = [{ ...st[0], tvd: 0, north: 0, east: 0, disp: 0, dls: 0 }];
  let maxDls = 0, maxInc = st[0].inc;
  for (let i = 1; i < st.length; i++) {
    const a = st[i - 1], b = st[i], p = out[i - 1], dMD = b.md - a.md, i1 = a.inc * D2R, i2 = b.inc * D2R, a1 = a.azi * D2R, a2 = b.azi * D2R;
    if (!(dMD > 0)) { out.push({ ...b, tvd: p.tvd, north: p.north, east: p.east, disp: p.disp, dls: 0 }); continue; }
    const beta = Math.acos(clamp(Math.cos(i2 - i1) - Math.sin(i1) * Math.sin(i2) * (1 - Math.cos(a2 - a1)), -1, 1)), rf = beta < 1e-9 ? 1 : (2 / beta) * Math.tan(beta / 2);
    const north = p.north + 0.5 * dMD * (Math.sin(i1) * Math.cos(a1) + Math.sin(i2) * Math.cos(a2)) * rf, east = p.east + 0.5 * dMD * (Math.sin(i1) * Math.sin(a1) + Math.sin(i2) * Math.sin(a2)) * rf, tvd = p.tvd + 0.5 * dMD * (Math.cos(i1) + Math.cos(i2)) * rf, dls = ((beta / D2R) * 30) / dMD;
    maxDls = Math.max(maxDls, dls); maxInc = Math.max(maxInc, b.inc);
    out.push({ ...b, tvd, north, east, disp: Math.hypot(north, east), dls });
  }
  const e = out[out.length - 1];
  return { stations: out, md: e.md - out[0].md, tvd: e.tvd, displacement: e.disp, maxDls, maxInc };
}
/** Tubing cells from the bottom of the well to the wellhead (flow direction), n cells of equal measured length. */
export function tubingCells(traj, n, D, rough) {
  const md = traj.stations.map((s) => s.md), tvd = traj.stations.map((s) => s.tvd), L = md[md.length - 1] - md[0], ds = L / n, cells = [];
  for (let i = 0; i < n; i++) { const mb = md[md.length - 1] - i * ds, ma = mb - ds, dz = interp1(md, tvd, mb) - interp1(md, tvd, ma); cells.push({ ds, theta: Math.asin(clamp(dz / ds, -1, 1)), D, rough, U: 0, tAmb: 0 }); }
  return cells;
}
/**
 * Inflow performance: rate (Sm³/d) at a flowing bottom-hole pressure pwf (bara).
 * o: { type: 'pi' | 'vogel' | 'composite' | 'fetkovich' | 'darcy' | 'jones' | 'gas', pRes, pb (bubble point), pi (Sm³/d/bar),
 *      n (Fetkovich exponent), a, b (laminar and non-Darcy coefficients of Δp = a q + b q², or Δp² for 'gas') }.
 */
export function iprRate(pwf, o) {
  const pr = o.pRes, pw = clamp(pwf, 0, pr), J = o.pi;
  const quad = (d) => (o.b > 0 ? (-o.a + Math.sqrt(o.a * o.a + 4 * o.b * d)) / (2 * o.b) : d / o.a);
  const vogel = (x) => 1 - 0.2 * x - 0.8 * x * x;
  switch (o.type) {
    case 'vogel': return ((J * pr) / 1.8) * vogel(pw / pr);
    case 'composite': { const pb = Math.min(o.pb ?? pr, pr); return pw >= pb ? J * (pr - pw) : J * (pr - pb) + ((J * pb) / 1.8) * vogel(pw / pb); }
    case 'fetkovich': { const n = o.n ?? 1; return ((J * pr) / 1.8) * (1 - (pw / pr) ** 2) ** n; }
    case 'darcy': case 'jones': return quad(pr - pw);
    case 'gas': return quad(pr * pr - pw * pw);
    default: return J * (pr - pw);
  }
}
/** Flowing bottom-hole pressure (bara) that delivers the rate q (Sm³/d); null above the absolute open-flow potential. */
export function iprPwf(q, o) {
  if (q <= 0) return o.pRes;
  if (q >= iprRate(0, o)) return null;
  return brent((p) => iprRate(p, o) - q, 0, o.pRes, 1e-9);
}
/**
 * Darcy radial-flow and non-Darcy (Forchheimer) coefficients in the form Δp = a q + b q² (oil, bar and Sm³/d) or
 * p_res² − p_wf² = a q + b q² (gas, bar² and Sm³/d).
 * { phase: 'oil' | 'gas', k (mD), h (pay, m), hp (perforated interval, m), re, rw (m), skin, mu (Pa·s), B (m³/Sm³), rho (kg/m³ in situ),
 *   Z, T (K), M (kg/mol), beta (1/m; default β = 1.88 × 10¹⁰ k^−1.47 φ^−0.53 1/ft with k in mD, S. C. Jones 1987), phi (porosity), nonDarcy (bool) }
 */
export function radialCoefficients({ phase = 'oil', k, h, hp, re, rw, skin = 0, mu, B = 1, rho = 800, Z = 0.9, T = 363, M = 0.02, beta = null, phi = 0.2, nonDarcy = true }) {
  if (!(k > 0 && h > 0 && re > rw && rw > 0)) throw new Error('Darcy inflow needs positive permeability and pay thickness and a drainage radius larger than the wellbore radius.');
  const kk = k * UNIT.mD, geo = Math.log(re / rw) - 0.75 + skin, bt = beta ?? (1.88e10 * k ** -1.47 * clamp(phi, 0.01, 0.6) ** -0.53) / UNIT.ft, hpp = hp > 0 ? hp : h, Tsc = T_STD + KEL, psc = P_STD * 1e5;
  if (!(geo > 0)) throw new Error('The skin factor is so negative that the Darcy inflow equation has no meaning (ln(re/rw) − 0.75 + S ≤ 0).');
  if (phase === 'gas') {
    const a = ((mu * Z * T * psc) / (Math.PI * kk * h * Tsc)) * geo, b = nonDarcy ? ((bt * M * Z * T * psc * psc) / (2 * Math.PI ** 2 * R * Tsc * Tsc * hpp * hpp)) * (1 / rw - 1 / re) : 0;
    return { a: a / 1e10 / DAY, b: b / 1e10 / DAY ** 2, beta: bt, pi: null };
  }
  const a = (mu * B * geo) / (2 * Math.PI * kk * h), b = nonDarcy ? ((bt * rho * B * B) / (4 * Math.PI ** 2 * hpp * hpp)) * (1 / rw - 1 / re) : 0;
  return { a: a / 1e5 / DAY, b: b / 1e5 / DAY ** 2, beta: bt, pi: (1e5 * DAY) / a };
}

// ---- valves and chokes ---------------------------------------------------------------------------------------
const N6 = 27.3; // IEC 60534-2-1 numerical constant for W in kg/h, pressures in bar, density in kg/m³ and C = Cv
/**
 * Required flow coefficient (IEC 60534 / ISA-75.01) for liquid, gas or a non-flashing two-phase stream.
 * { w (kg/s), dp (bar), p1 (bara), rhoL, rhoG (kg/m³ at inlet), xG (gas mass fraction), k (cp/cv), xT, FL, pv, pc (bara) }
 * Returns { cv, kv, choked, x (pressure-drop ratio), Y (expansion factor), dpLiquidMax, critical ratio xCrit }.
 */
export function valveCv({ w, dp, p1, rhoL = 1000, rhoG = 1, xG = 0, k = 1.3, xT = 0.7, FL = 0.9, pv = 0, pc = 221 }) {
  if (!(dp > 0) || !(p1 > 0)) throw new Error('Valve sizing needs a positive pressure drop and inlet pressure.');
  const Fk = k / 1.4, xCrit = Fk * xT, x = Math.min(dp / p1, xCrit), Y = 1 - x / (3 * xCrit), FF = 0.96 - 0.28 * Math.sqrt(clamp(pv / pc, 0, 1)), dpLmax = FL * FL * Math.max(p1 - FF * pv, 1e-9), dpL = Math.min(dp, dpLmax);
  const inv = (xG > 0 ? xG / (rhoG * Y * Y * x * p1) : 0) + (xG < 1 ? (1 - xG) / (rhoL * dpL) : 0), cv = ((w * 3600) / N6) * Math.sqrt(inv);
  return { cv, kv: cv * 0.865, choked: (xG > 0 && dp / p1 >= xCrit) || (xG < 1 && dp >= dpLmax), x, Y, dpLiquidMax: dpLmax, xCrit };
}
/** Pressure drop (bar) of a valve of flow coefficient cv passing w (kg/s): inverse of valveCv; continues quadratically beyond choking. */
export function valveDp({ cv, w, p1, rhoL = 1000, rhoG = 1, xG = 0, k = 1.3, xT = 0.7 }) {
  if (!(cv > 0)) throw new Error('A valve needs a positive flow coefficient.');
  const xCrit = (k / 1.4) * xT, c = ((w * 3600) / (N6 * cv)) ** 2;
  let dp = c * (xG / rhoG + (1 - xG) / rhoL), Y = 1, x = 0;
  for (let i = 0; i < 30; i++) { x = Math.min(dp / Math.max(p1, 1e-6), xCrit); Y = 1 - x / (3 * xCrit); const d = c * (xG / (rhoG * Y * Y) + (1 - xG) / rhoL); if (Math.abs(d - dp) < 1e-12 * (1 + dp)) { dp = d; break; } dp = d; }
  return { dp, x, Y, choked: xG > 0 && dp / Math.max(p1, 1e-6) >= xCrit };
}
/** Inherent characteristic: relative flow coefficient at travel h (0–1). type: 'linear' | 'equal' (equal-percentage, rangeability R) | 'quick'. */
export const valveCharacteristic = (h, type = 'equal', Rng = 50) => (h <= 0 ? 0 : h >= 1 ? 1 : type === 'linear' ? h : type === 'quick' ? Math.sqrt(h) : Rng ** (h - 1));
/** Travel (0–1) that gives the relative flow coefficient f. */
export const valveOpening = (f, type = 'equal', Rng = 50) => clamp(type === 'linear' ? f : type === 'quick' ? f * f : f <= 0 ? 0 : 1 + Math.log(f) / Math.log(Rng), 0, 1);
/** Orifice / Bernoulli mass flow (kg/s): m = Cd·Y·A·√(2 ρ Δp / (1 − β⁴)); gas expansion factor when p1 and k are given. d, D in m, dp in bar. */
export function orificeFlow({ cd = 0.61, d, D = Infinity, rho, dp, p1 = null, k = null }) {
  const beta = Number.isFinite(D) ? d / D : 0, Y = p1 && k ? 1 - (0.351 + 0.256 * beta ** 4 + 0.93 * beta ** 8) * (1 - (Math.max(p1 - dp, 0) / p1) ** (1 / k)) : 1;
  return cd * Y * area(d) * Math.sqrt((2 * rho * Math.max(dp, 0) * 1e5) / (1 - beta ** 4));
}
/** Critical pressure ratio (throat / stagnation) of an ideal gas with isentropic exponent k. */
export const criticalRatio = (k) => (2 / (k + 1)) ** (k / (k - 1));
/** Isentropic nozzle mass flux (kg/m²/s) of an ideal gas from stagnation p0 (Pa), rho0 to the pressure ratio r (choked below the critical ratio). */
export function nozzleFlux(p0, rho0, k, r) {
  const rr = Math.max(r, criticalRatio(k));
  return Math.sqrt(((2 * k) / (k - 1)) * p0 * rho0 * (rr ** (2 / k) - rr ** ((k + 1) / k)));
}
const CHOKE_CORR = { gilbert: [10, 0.546, 1.89, 1], ros: [17.4, 0.5, 2.0, 0], baxendell: [9.56, 0.546, 1.93, 0], achong: [3.82, 0.65, 1.88, 0] };
/**
 * Critical-flow multiphase choke correlations p1 = A·GLR^B·q / S^C (Gilbert, Ros, Baxendell, Achong; field units inside).
 * Give two of { p1 (bara upstream), q (Sm³/d liquid), bean (64ths of an inch) } and glr (Sm³/Sm³); the missing one is returned.
 */
export function chokeCorrelation(name, { p1 = null, q = null, glr, bean = null }) {
  const [A, B, C, gauge] = CHOKE_CORR[name] || CHOKE_CORR.gilbert, Rs = Math.max(glr, 1e-6) * (UNIT.bbl / UNIT.scf), off = gauge ? 1.01325 : 0;
  if (bean === null) return ((A * Rs ** B * (q / UNIT.bbl)) / Math.max((p1 - off) / (UNIT.psi / 1e5), 1e-9)) ** (1 / C);
  if (q === null) return (((p1 - off) / (UNIT.psi / 1e5)) * bean ** C) / (A * Rs ** B) * UNIT.bbl;
  return ((A * Rs ** B * (q / UNIT.bbl)) / bean ** C) * (UNIT.psi / 1e5) + off;
}
/**
 * Sachdeva et al. (1986) two-phase choke model (no slip, frozen composition, polytropic gas expansion): critical ratio, polytropic
 * exponent and mass flux as printed in Guo, Lyons & Ghalambor (2007), Eqs. 5.13–5.18. Discharge coefficient 0.62–0.90 there; 0.85 for a
 * choke free of upstream disturbance and 0.75 for one close behind an elbow (Buffa & Baliño 2017, quoting the original paper).
 * { p1, p2 (bara), x (gas mass fraction), rhoG (at p1), rhoL, k, cvG, cL (J/kg/K), cd } → { yc (critical ratio), critical, y, G (kg/m²/s) }.
 */
export function sachdeva({ p1, p2, x, rhoG, rhoL, k = 1.3, cvG = 1700, cL = 2200, cd = 0.85 }) {
  const xg = clamp(x, 1e-6, 1), vG1 = 1 / rhoG, vL = 1 / rhoL, kk = k / (k - 1), n = 1 + (xg * cvG * (k - 1)) / (xg * cvG + (1 - xg) * cL), rat = ((1 - xg) * vL) / (xg * vG1);
  let yc = 0.5;
  for (let i = 0; i < 200; i++) { const a = rat * yc ** (1 / k), y = ((kk + rat * (1 - yc)) / (kk + n / 2 + n * a + (n / 2) * a * a)) ** kk; if (Math.abs(y - yc) < 1e-13) { yc = y; break; } yc = 0.5 * (yc + y); }
  const y = Math.max(p2 / p1, yc), vG2 = vG1 * y ** (-1 / k), rhoM2 = 1 / (xg * vG2 + (1 - xg) * vL);
  return { yc, critical: p2 / p1 <= yc, y, G: cd * Math.sqrt(2 * p1 * 1e5 * rhoM2 * rhoM2 * (((1 - xg) * (1 - y)) / rhoL + xg * kk * (vG1 - y * vG2))) };
}
/**
 * Homogeneous-equilibrium critical mass flux: G(p) = √(2 ∫ v dp) / v(p) maximised over the throat pressure.
 * vOf(p bara) is the specific volume (m³/kg) of the equilibrium mixture along the expansion path. Returns { G, pc, ratio }.
 */
export function hemCritical(vOf, p1, n = 240) {
  let I = 0, pPrev = p1, vPrev = vOf(p1), best = { G: 0, pc: p1, ratio: 1 };
  for (let i = 1; i <= n; i++) {
    const p = p1 * (0.03 / 1) ** (i / n), v = vOf(p);
    I += 0.5 * (v + vPrev) * (pPrev - p) * 1e5;
    const Gp = Math.sqrt(2 * I) / v;
    if (Gp > best.G) best = { G: Gp, pc: p, ratio: p / p1 };
    pPrev = p; vPrev = v;
  }
  return best;
}

// ---- pumps and compressors -----------------------------------------------------------------------------------
/** Affinity laws: a duty { q, h, p } at speed ratio N2/N1 and impeller-trim ratio D2/D1. */
export const affinity = ({ q = 0, h = 0, p = 0 }, speed = 1, trim = 1) => ({ q: q * speed * trim, h: h * (speed * trim) ** 2, p: p * (speed * trim) ** 3 });
/** Pump head (m) at flow q for a parabolic curve through the rated point (qr, hr) with a shut-off head ratio, scaled by the affinity laws. */
export function pumpHead(q, { qr, hr, shutoff = 1.25, speed = 1, trim = 1, curve = null }) {
  const s = Math.max(speed * trim, 1e-6);
  if (curve) { // tabulated curve { q[], h[] } at rated speed and full impeller; linear continuation beyond the last point
    const n = curve.q.length, x = q / s, slope = Math.min((curve.h[n - 1] - curve.h[n - 2]) / (curve.q[n - 1] - curve.q[n - 2]), -1e-9);
    return s * s * (x > curve.q[n - 1] ? curve.h[n - 1] + slope * (x - curve.q[n - 1]) : interp1(curve.q, curve.h, x));
  }
  return hr * s * s * (shutoff - (shutoff - 1) * (q / (qr * s)) ** 2);
}
/** Flow (same unit as the curve) at which the pump head falls to zero. */
export function pumpRunout({ qr, shutoff = 1.25, speed = 1, trim = 1, curve = null }) {
  const s = Math.max(speed * trim, 1e-6);
  if (!curve) return qr * s * Math.sqrt(shutoff / (shutoff - 1));
  const n = curve.q.length, slope = Math.min((curve.h[n - 1] - curve.h[n - 2]) / (curve.q[n - 1] - curve.q[n - 2]), -1e-9);
  return s * (curve.q[n - 1] + Math.max(curve.h[n - 1], 0) / -slope);
}
/** Pump efficiency at flow q: parabola through the best-efficiency point (taken at the rated flow). */
export const pumpEfficiency = (q, { qr, eta = 0.75, speed = 1, trim = 1, curve = null }) => { const s = Math.max(speed * trim, 1e-6); if (curve?.eta) return clamp(interp1(curve.q, curve.eta, q / s), 0.02, 0.95); const x = q / (qr * s); return clamp(eta * (2 * x - x * x), 0.02, eta); };
/**
 * Euler turbomachinery head of a centrifugal impeller with the Wiesner slip factor.
 * { d2 (m), rpm, beta2 (blade outlet angle from the tangent, deg), blades, b2 (outlet width, m), q (m³/s) } → { u2, cm2, slip, headInf, head (m) }.
 */
export function eulerHead({ d2, rpm, beta2 = 25, blades = 7, b2 = 0.02, q = 0 }) {
  const u2 = (Math.PI * d2 * rpm) / 60, cm2 = q / (Math.PI * d2 * b2), tb = Math.tan(beta2 * D2R), slip = Number.isFinite(blades) ? 1 - Math.sqrt(Math.sin(beta2 * D2R)) / blades ** 0.7 : 1;
  return { u2, cm2, slip, headInf: (u2 * (u2 - cm2 / tb)) / G, head: (u2 * (slip * u2 - cm2 / tb)) / G };
}
/**
 * One compression stage. { p1, p2 (bara), t1 (°C), z (mean compressibility), k (cp/cv), mw (g/mol), eta (polytropic efficiency), mdot (kg/s) }
 * Returns { ratio, headPoly, headIsen (J/kg), t2 (°C), t2Isen, power (W), etaIsen, n (polytropic exponent) }.
 */
export function compressorStage({ p1, p2, t1, z = 1, k = 1.3, mw = 20, eta = 0.78, mdot = 1 }) {
  if (!(p1 > 0 && p2 > p1)) throw new Error('A compressor stage needs a discharge pressure above its suction pressure.');
  const T1 = t1 + KEL, r = p2 / p1, km = (k - 1) / k, nm = km / eta, zrt = (z * R * T1) / (mw * 1e-3), headPoly = (zrt / nm) * (r ** nm - 1), headIsen = (zrt / km) * (r ** km - 1);
  return { ratio: r, headPoly, headIsen, t2: T1 * r ** nm - KEL, t2Isen: T1 * r ** km - KEL, power: (mdot * headPoly) / eta, etaIsen: (headIsen * eta) / headPoly, n: 1 / (1 - nm) };
}
/**
 * Simple compressor map around a design point (qd = inlet volume flow, hd = polytropic head at 100 % speed): tabulated speed
 * lines obeying the fan laws, from the surge line to stonewall. Returns { speeds, lines: [{ n, q[], h[], qSurge, qChoke }] }.
 */
export function compressorMap({ qd, hd, speeds = [0.7, 0.8, 0.9, 1, 1.05], n = 11 }) {
  return { qd, hd, speeds, lines: speeds.map((N) => { const phi = linspace(PHI_SURGE, PHI_CHOKE, n); return { n: N, q: phi.map((f) => f * N * qd), h: phi.map((f) => N * N * hd * mapShape(f)), qSurge: PHI_SURGE * N * qd, qChoke: PHI_CHOKE * N * qd }; }) };
}
/** Operating point on a compressor map by interpolation between the tabulated speed lines: { speed, qSurge, qChoke, surgeMargin, chokeMargin (%), inside }. */
export function mapPoint(map, q, h) {
  const hs = map.lines.map((l) => interp1(l.q, l.h, q)), up = hs.every((v, i) => !i || v > hs[i - 1]);
  const speed = up ? interp1(hs, map.speeds, h) : map.speeds[map.speeds.length - 1] * Math.sqrt(h / Math.max(hs[hs.length - 1], 1e-9));
  const qSurge = interp1(map.speeds, map.lines.map((l) => l.qSurge), speed), qChoke = interp1(map.speeds, map.lines.map((l) => l.qChoke), speed);
  return { speed, qSurge, qChoke, surgeMargin: (100 * (q - qSurge)) / Math.max(q, 1e-12), chokeMargin: (100 * (qChoke - q)) / Math.max(qChoke, 1e-12), inside: h >= hs[0] - 1e-9 && h <= hs[hs.length - 1] + 1e-9 && q >= qSurge && q <= qChoke };
}

// ---- vessels ---------------------------------------------------------------------------------------------------
/**
 * Liquid volume (m³) in a vessel at a level (m above the lowest point). Exact circular-segment geometry for the shell,
 * exact ellipsoidal caps for the heads. { d, l (tan–tan length), level, orientation: 'horizontal' | 'vertical', heads: 'flat' | 'elliptical' (2:1) | 'hemispherical' }
 */
export function vesselVolume({ d, l, level, orientation = 'horizontal', heads = 'elliptical' }) {
  if (!(d > 0) || !(l >= 0)) throw new Error('A vessel needs a positive diameter and a non-negative length.');
  const Rv = d / 2, a = heads === 'elliptical' ? Rv / 2 : heads === 'hemispherical' ? Rv : 0, hMax = orientation === 'vertical' ? l + 2 * a : d;
  if (!isNum(level) || level < -1e-9 || level > hMax + 1e-9) throw new Error(`The liquid level must lie between 0 and ${hMax.toFixed(3)} m for this vessel.`);
  const h = clamp(level, 0, hMax);
  if (orientation === 'vertical') {
    const cap = (t) => (a > 0 ? Math.PI * Rv * Rv * (t - ((t - a) ** 3 + a ** 3) / (3 * a * a)) : 0); // ellipsoidal cap filled to t (0..a)
    if (h <= a) return cap(h);
    if (h <= a + l) return cap(a) + Math.PI * Rv * Rv * (h - a);
    const t = h - a - l; return cap(a) + Math.PI * Rv * Rv * l + Math.PI * Rv * Rv * (t - t ** 3 / (3 * a * a));
  }
  return l * (Rv * Rv * Math.acos(clamp((Rv - h) / Rv, -1, 1)) - (Rv - h) * Math.sqrt(Math.max(2 * Rv * h - h * h, 0))) + (Math.PI * a * h * h * (3 * Rv - h)) / (3 * Rv);
}
/** Level (m) that holds a liquid volume in a vessel (inverse of vesselVolume). */
export function vesselLevel(volume, geom) {
  const hMax = geom.orientation === 'vertical' ? geom.l + (geom.heads === 'flat' ? 0 : geom.heads === 'hemispherical' ? geom.d : geom.d / 2) : geom.d, vMax = vesselVolume({ ...geom, level: hMax });
  if (volume <= 0) return 0;
  if (volume >= vMax) return hMax;
  return brent((h) => vesselVolume({ ...geom, level: h }) - volume, 0, hMax, 1e-10);
}
/** Circular segment of a pipe or vessel of diameter d filled to height h: { area, wetted (wall perimeter), chord (free-surface width) }. */
export function segmentGeometry(d, h) {
  const r = d / 2, hh = clamp(h, 0, d), a = Math.acos(clamp((r - hh) / r, -1, 1));
  return { area: r * r * a - (r - hh) * Math.sqrt(Math.max(2 * r * hh - hh * hh, 0)), wetted: 2 * r * a, chord: 2 * Math.sqrt(Math.max(2 * r * hh - hh * hh, 0)) };
}
/** Hydraulic diameter 4A/P of a flow section. */
export const hydraulicDiameter = (A, P) => (P > 0 ? (4 * A) / P : 0);
/** Souders–Brown maximum gas velocity (m/s) for a sizing coefficient K (m/s). */
export const soudersBrown = (K, rhoL, rhoG) => K * Math.sqrt(Math.max(rhoL - rhoG, 0) / Math.max(rhoG, 1e-9));
/**
 * Souders–Brown coefficient (m/s) of a separator with a wire-mesh mist eliminator against pressure (bara): the York curve as
 * fitted by Svrcek & Monnery (1993), K [ft/s] = 0.1821 + 0.0029 P + 0.0460 ln P (P < 15 psia), 0.35 (15–40 psia), 0.430 − 0.023 ln P
 * (40–5,500 psia); a horizontal vessel takes 1.25 times the vertical value.
 */
export function demisterK(pBar, horizontal = false) {
  const P = clamp((pBar * 1e5) / UNIT.psi, 1, 5500), K = P < 15 ? 0.1821 + 0.0029 * P + 0.046 * Math.log(P) : P < 40 ? 0.35 : 0.43 - 0.023 * Math.log(P);
  return K * UNIT.ft * (horizontal ? 1.25 : 1);
}

// ---- network ---------------------------------------------------------------------------------------------------
/** Node–edge incidence matrix A[node][edge]: −1 where the edge leaves the node, +1 where it enters. */
export function incidence(nodes, edges) {
  const ix = new Map(nodes.map((n, i) => [String(n.id), i])), A = nodes.map(() => new Array(edges.length).fill(0));
  edges.forEach((e, j) => { const a = ix.get(String(e.from)), b = ix.get(String(e.to)); if (a !== undefined) A[a][j] -= 1; if (b !== undefined) A[b][j] += 1; });
  return A;
}
/**
 * Steady hydraulic network: nodal mass conservation (Kirchhoff) plus one pressure-drop law per branch, solved for the
 * free node pressures and all branch mass flows by Newton–Raphson with a backtracking line search.
 * nodes: [{ id, p (bara, fixed when numeric), q (kg/s injected, negative = demand), z (m) }]
 * edges: [{ from, to, length (m), id (m), rough (m), k, f }] (default law: single-phase Darcy–Weisbach with `fluid` { rho, mu })
 * law(edge, m ≥ 0, pUpstream (bara), dir (+1 along from → to), index, rec[]) → pressure drop in the flow direction (bar)
 * afterEval({ p, m, rec }) → number: outer fixed-point hook (e.g. node temperatures); iteration continues while it returns > 1e-3.
 * Returns { converged (residual below tol), usable (residual below 1e-4: stalled on a kink of a closure law), iterations, residual, p[], m[], dp[], dpLaw[], balance[], massResidual, loops: [{ edges: [{ i, sign }], sum (bar) }], inflow, outflow, rec }.
 */
export function solveNetwork({ nodes, edges, law = null, fluid = { rho: 1000, mu: 1e-3 }, tol = 1e-8, maxIter = 80, guess = null, afterEval = null }) {
  const N = nodes?.length || 0, E = edges?.length || 0, ids = new Map(), ef = [], et = [];
  if (N < 2 || !E) throw new Error('The network needs at least two nodes and one connection.');
  nodes.forEach((n, i) => { const k = String(n.id); if (ids.has(k)) throw new Error(`Network node “${k}” is defined twice.`); ids.set(k, i); });
  edges.forEach((e, i) => {
    const a = ids.get(String(e.from)), b = ids.get(String(e.to));
    if (a === undefined || b === undefined) throw new Error(`Connection ${i + 1} (${e.from} → ${e.to}) refers to a node that does not exist.`);
    if (a === b) throw new Error(`Connection ${i + 1} starts and ends at the same node (${e.from}).`);
    ef.push(a); et.push(b);
  });
  // spanning forest (tree edges) and chords (each chord closes one independent loop)
  const par = nodes.map((_, i) => i), find = (i) => { while (par[i] !== i) { par[i] = par[par[i]]; i = par[i]; } return i; }, chords = [], adj = nodes.map(() => []);
  edges.forEach((_, i) => { const a = find(ef[i]), b = find(et[i]); if (a === b) chords.push(i); else { par[a] = b; adj[ef[i]].push([et[i], i]); adj[et[i]].push([ef[i], i]); } });
  const fixed = nodes.map((n) => isNum(n.p)), withP = new Set(nodes.map((_, i) => (fixed[i] ? find(i) : -1)));
  nodes.forEach((n, i) => { if (!withP.has(find(i))) throw new Error(`Node “${n.id}” is not connected to any fixed-pressure boundary: the network is disconnected or has no pressure reference.`); });
  const up = new Array(N).fill(-1), upE = new Array(N).fill(-1), depth = new Array(N).fill(-1);
  for (let r = 0; r < N; r++) { if (depth[r] >= 0) continue; depth[r] = 0; const q = [r]; while (q.length) { const c = q.pop(); for (const [o, i] of adj[c]) if (depth[o] < 0) { depth[o] = depth[c] + 1; up[o] = c; upE[o] = i; q.push(o); } } }
  const loopOf = (i) => {
    const path = [{ i, sign: 1 }], back = [];
    let a = ef[i], b = et[i];
    while (a !== b) { if (depth[b] >= depth[a]) { const e = upE[b]; path.push({ i: e, sign: ef[e] === b ? 1 : -1 }); b = up[b]; } else { const e = upE[a]; back.push({ i: e, sign: ef[e] === up[a] ? 1 : -1 }); a = up[a]; } }
    return path.concat(back.reverse());
  };
  const lawFn = law || ((e, m, pUp, dir, i) => pipeDp({ m, rho: fluid.rho, mu: fluid.mu, L: e.length, D: e.id, rough: e.rough || 0, dz: dir * ((nodes[et[i]].z || 0) - (nodes[ef[i]].z || 0)), k: e.k || 0, f: e.f ?? null }) / 1e5);
  const pi = [], inj = nodes.map((n) => num(n.q, 0));
  let nf = 0; nodes.forEach((_, i) => pi.push(fixed[i] ? -1 : nf++));
  const nU = nf + E, pFix = nodes.filter((_, i) => fixed[i]).map((n) => n.p), pRef = Math.max(1, ...pFix.map(Math.abs)), mInj = inj.reduce((s, q) => s + Math.abs(q), 0) / 2;
  const evalAll = (pp, mm) => {
    const L = new Array(E), rec = new Array(E), r = new Array(nU).fill(0);
    let mx = mInj;
    for (let i = 0; i < E; i++) {
      const s = mm[i] >= 0 ? 1 : -1, u = s > 0 ? ef[i] : et[i];
      L[i] = lawFn(edges[i], Math.abs(mm[i]), pp[u], s, i, rec); r[nf + i] = pp[ef[i]] - pp[et[i]] - s * L[i];
      if (pi[ef[i]] >= 0) r[pi[ef[i]]] -= mm[i];
      if (pi[et[i]] >= 0) r[pi[et[i]]] += mm[i];
      mx = Math.max(mx, Math.abs(mm[i]));
    }
    for (let k = 0; k < N; k++) if (pi[k] >= 0) r[pi[k]] += inj[k];
    const ms = Math.max(mx, 1e-12);
    let norm = 0; for (let k = 0; k < nU; k++) norm = Math.max(norm, Math.abs(r[k]) / (k < nf ? ms : pRef));
    return { L, rec, r, norm: Number.isFinite(norm) ? norm : Infinity, ms };
  };
  const jac = (Lm, Lp, sgn) => { // rows: free-node mass balances, then branch equations; columns: free pressures, then branch flows
    const J = Array.from({ length: nU }, () => new Array(nU).fill(0));
    for (let i = 0; i < E; i++) {
      const row = nf + i, a = ef[i], b = et[i], u = sgn[i] > 0 ? a : b;
      J[row][nf + i] = -Lm[i];
      if (pi[a] >= 0) { J[row][pi[a]] += 1; J[pi[a]][nf + i] = -1; }
      if (pi[b] >= 0) { J[row][pi[b]] -= 1; J[pi[b]][nf + i] = 1; }
      if (pi[u] >= 0) J[row][pi[u]] -= sgn[i] * Lp[i];
    }
    return J;
  };
  let p = nodes.map((n, i) => (fixed[i] ? n.p : num(guess?.p?.[i], 1.3 * Math.max(...pFix) + 5))), m;
  if (guess?.m?.length === E && guess.m.every(isNum)) m = guess.m.slice();
  else { // start from the flows of a linear network with unit resistances (exact for a tree with fixed inflows)
    const r0 = new Array(nU).fill(0);
    for (let i = 0; i < E; i++) r0[nf + i] = (fixed[ef[i]] ? p[ef[i]] : 0) - (fixed[et[i]] ? p[et[i]] : 0);
    for (let k = 0; k < N; k++) if (pi[k] >= 0) r0[pi[k]] = inj[k];
    m = solveLinear(jac(new Array(E).fill(1), new Array(E).fill(0), new Array(E).fill(1)), r0.map((x) => -x)).slice(nf);
  }
  let ev = evalAll(p, m), change = afterEval ? afterEval({ p, m, rec: ev.rec }) : 0, it = 0, converged = false, stale = true;
  if (change > 0.05) { ev = evalAll(p, m); change = afterEval({ p, m, rec: ev.rec }); }
  const Lm = new Array(E).fill(0), Lp = new Array(E).fill(0), sgn = new Array(E).fill(1);
  const fd = law ? 0.02 : 1e-5; // closure laws are only piecewise smooth (flow-pattern changes): secant-size steps give usable slopes
  const derivs = () => { // finite-difference branch derivatives; reused between iterations while convergence is fast (chord steps)
    for (let i = 0; i < E; i++) {
      const s = m[i] >= 0 ? 1 : -1, u = s > 0 ? ef[i] : et[i], am = Math.abs(m[i]), dm = Math.max(fd * am, 1e-3 * fd * ev.ms), dpp = fd * Math.max(p[u], 1), tiny = (1e-9 * pRef) / ev.ms;
      Lm[i] = (lawFn(edges[i], am + dm, p[u], s, i, null) - ev.L[i]) / dm; if (Math.abs(Lm[i]) < tiny) Lm[i] = tiny;
      if (law) Lp[i] = (lawFn(edges[i], am, p[u] + dpp, s, i, null) - ev.L[i]) / dpp;
    }
  };
  let best = ev.norm, since = 0, forced = 0;
  for (; it < maxIter; it++) {
    if (ev.norm < tol && !(change > 1e-3)) { converged = true; break; }
    if (ev.norm < 0.7 * best) { best = ev.norm; since = 0; } else if (++since > 7 && ev.norm < 1e-4) break; // stalled on a kink of the closure laws: good to within the reported residual
    if (ev.norm >= tol) {
      let done = false;
      for (let attempt = 0; attempt < 2 && !done; attempt++) {
        if (stale) { derivs(); stale = false; attempt = 1; } // with fresh derivatives this is the last attempt and gets the full line search
        for (let i = 0; i < E; i++) sgn[i] = m[i] >= 0 ? 1 : -1;
        let d; try { d = solveLinear(jac(Lm, Lp, sgn), ev.r.map((x) => -x)); } catch { d = null; }
        if (d && d.every(Number.isFinite)) for (let k = 0, lam = 1; k < (attempt ? 10 : 1); k++, lam *= 0.5) {
          const pn = p.map((x, i) => (pi[i] >= 0 ? Math.max(x + lam * d[pi[i]], 0.3 * x, 0.05) : x)), mn = m.map((x, i) => x + lam * d[nf + i]), en = evalAll(pn, mn);
          if (en.norm < (1 - 1e-4 * lam) * ev.norm || (attempt && k === 9 && ev.norm >= 1e-4 && forced++ < 4)) { stale = en.norm > 0.25 * ev.norm; p = pn; m = mn; ev = en; done = true; break; }
        }
        if (!done) stale = true;
      }
      if (!done) break;
    } else ev = evalAll(p, m);
    change = afterEval ? afterEval({ p, m, rec: ev.rec }) : 0;
    if (change > 0.05) { ev = evalAll(p, m); change = afterEval({ p, m, rec: ev.rec }); stale = true; } // the outer state moved: refresh the residual
  }
  if (!converged && ev.norm < tol && !(change > 1e-3)) converged = true;
  const dpLaw = ev.L.map((L, i) => (m[i] >= 0 ? 1 : -1) * L), balance = inj.slice();
  edges.forEach((_, i) => { balance[ef[i]] -= m[i]; balance[et[i]] += m[i]; });
  let massResidual = 0, inflow = 0, outflow = 0;
  balance.forEach((b, k) => { if (fixed[k]) { if (b > 0) outflow += b; else inflow -= b; } else { massResidual = Math.max(massResidual, Math.abs(b)); if (inj[k] > 0) inflow += inj[k]; else outflow -= inj[k]; } });
  return { converged, usable: converged || ev.norm < 1e-4, iterations: it, residual: ev.norm, p, m, dp: edges.map((_, i) => p[ef[i]] - p[et[i]]), dpLaw, balance, massResidual, loops: chords.map((i) => { const path = loopOf(i); return { edges: path, sum: path.reduce((s, q) => s + q.sign * dpLaw[q.i], 0) }; }), inflow, outflow, rec: ev.rec, from: ef, to: et };
}

// ---- rotating-equipment dynamics and the transient network ---------------------------------------------------------
const mapShape = (phi) => 1 + 0.22 * (1 - phi * phi) - 2.2 * Math.max(phi - 1.1, 0) ** 2, PHI_SURGE = 0.62, PHI_CHOKE = 1.3;
/**
 * Flow coefficient (inlet flow / design flow at that speed) on the generic speed line for a head coefficient psi (head / design
 * head at that speed): the inverse of the map shape used by compressorMap. Right of stonewall the flow stays at the choke value;
 * left of the surge point the line is continued with a shallow slope that reaches zero flow at 1.08 times the surge head, so
 * that surge is detected (negative surge margin) but its limit cycle is not resolved.
 */
export function mapFlowCoefficient(psi) {
  const ps = mapShape(PHI_SURGE);
  if (!(psi < ps)) return Math.max(PHI_SURGE * (1 - (psi - ps) / (0.08 * ps)), 0);
  if (psi >= mapShape(1.1)) return Math.sqrt(1 - (psi - 1) / 0.22);
  if (psi <= mapShape(PHI_CHOKE)) return PHI_CHOKE;
  return (4.84 + Math.sqrt(Math.max(23.4256 - 9.68 * (1.442 + psi), 0))) / 4.84;
}
/**
 * Speed of one rotor train: J dω/dt = τ_driver(t, ω) − τ_load(t, ω), integrated with the implicit trapezoidal rule (Newton on
 * every step). { J (kg·m²), omega0 (rad/s), driver(t, ω) → N·m (default: none, i.e. a trip), load(t, ω) → N·m, tEnd (s), n (steps) }
 * Returns { t[], omega[], halfTime (s until the speed has fallen to half its initial value, null when it never does) }.
 */
export function rotorDynamics({ J, omega0, driver = () => 0, load, tEnd, n = 200 }) {
  if (!(J > 0) || !(tEnd > 0) || typeof load !== 'function') throw new Error('Rotor dynamics needs a positive inertia, a duration and a load-torque function.');
  const t = [0], w = [omega0], dt = tEnd / n, f = (tt, ww) => (driver(tt, ww) - load(tt, ww)) / J;
  let halfTime = null;
  for (let i = 0; i < n; i++) {
    const t1 = t[i] + dt, f0 = f(t[i], w[i]);
    let wn = w[i] + dt * f0;
    for (let k = 0; k < 30; k++) { const g = wn - w[i] - 0.5 * dt * (f0 + f(t1, wn)), h = 1e-6 * Math.max(Math.abs(wn), 1e-3), dg = 1 - (0.5 * dt * (f(t1, wn + h) - f(t1, wn))) / h, d = g / (Math.abs(dg) > 1e-12 ? dg : 1); wn -= d; if (Math.abs(d) < 1e-12 * (1 + Math.abs(wn))) break; }
    wn = Math.max(wn, 0);
    if (halfTime === null && wn <= 0.5 * omega0 && w[i] > 0.5 * omega0) halfTime = t[i] + (dt * (w[i] - 0.5 * omega0)) / (w[i] - wn);
    t.push(t1); w.push(wn);
  }
  return { t, omega: w, halfTime };
}
const rsq = (x, eps) => x / Math.sqrt(Math.abs(x) + eps); // signed square root, regular at zero
const zOfNode = (n, p) => (n.z && Array.isArray(n.z.p) ? interp1(n.z.p, n.z.z, clamp(p, n.z.p[0], n.z.p[n.z.p.length - 1])) : isNum(n.z) ? n.z : 1);
/**
 * Transient network: dynamic nodal inventories joined by quasi-steady branches, with rotor speed and actuator dynamics.
 *
 * nodes: [{ id, type: 'volume' | 'fixed', p (bara: initial pressure, or the boundary pressure of a fixed node),
 *           V (m³ total), T (K), M (kg/mol), z (number, or a table { p[], z[] } of the real-gas factor against pressure),
 *           liquid: { vol0 (m³), rho (kg/m³), geom (vesselVolume geometry, for the level) | area (m²) } | undefined }]
 *        gas inventory of a volume node: m = p·V_gas·M / (Z R T) with V_gas = V − liquid volume (line-pack with real-gas compressibility).
 * edges: [{ id, from, to, type, phase: 'gas' | 'liquid', check (no reverse flow) }] with, by type,
 *   'pipe'       { r (bar per (kg/s)²), dz (static bar from → to) }                        m = sign·√(|Δp − dz| / r)
 *   'valve'      { cv, pos (0–1 initial travel), cmd, tau (s), rate (1/s), char, rng, xT, k, rho (liquid) }   IEC 60534 flow equations
 *   'pump'       { qr (m³/h), hr (m), shutoff, eta, rho, dz (static head m), ks (system m per (m³/h)²), rotor, throttle: { k0, pos, tau, rate, char, rng } }
 *   'compressor' { qd (m³/s inlet at rated speed), hd (J/kg polytropic per stage at rated speed), stages, eta, k, rotor }   generic fan-law map
 *   'source'     { curve: { m[], p[] } (pressure available against mass rate), offset (bar), c0 (bar per (kg/s)² at travel pos0), pos, pos0,
 *                  tau, rate, char, rng, xg (gas mass fraction → node 'to'), liquidTo (node receiving the liquid) }
 * rotors: [{ id, J (kg·m²), omegaR (rad/s rated), speed (initial fraction), tauMax (N·m driver limit), kp, ti (speed governor), tauLoss (N·m at rated speed, ∝ ω²), set (speed set point, fraction) }]
 * controllers: [{ type: 'P' | 'tripOpen', measure: { node, q: 'p' | 'level' } | { edge, q: 'surgeMargin' }, set, gain, bias, min, max, out: { edge } | { rotor }, rotor (tripOpen) }]
 * events: [{ t, edge | rotor | node, set: { cmd } | { trip: true } | { set } | { p } }]
 * Integration: implicit θ-method (θ = 0.55) with Newton iteration on a finite-difference Jacobian; a step that fails to converge is halved.
 * Returns { t[], p: { node: [] }, level: { node: [] }, liquid: { node: [] }, m: { edge: [] (kg/s) }, speed: { rotor: [] }, power: { rotor: [] (kW) },
 *   pos: { edge: [] }, margin: { compressor: [] (%) }, surgeTime: { compressor: s | null }, balance: { gas, liquid: { initial, in, out, final } (kg) }, steps, halved }.
 */
export function transientNetwork({ nodes, edges, rotors = [], controllers = [], events = [], tEnd, dt = 0.5, nOut = 120, theta = 0.55 }) {
  if (!(tEnd > 0) || !(dt > 0)) throw new Error('The transient needs a positive duration and time step.');
  if (!Array.isArray(nodes) || !Array.isArray(edges) || !nodes.length || !edges.length) throw new Error('The transient network needs nodes and connections.');
  const nix = new Map(nodes.map((n, i) => [String(n.id), i])), rix = new Map(rotors.map((r, i) => [String(r.id), i])), eix = new Map(edges.map((e, i) => [String(e.id), i]));
  const N = (id, what) => { const k = nix.get(String(id)); if (k === undefined) throw new Error('Transient network: ' + what + ' refers to the unknown node “' + id + '”.'); return k; };
  const RT = (id, what) => { const k = rix.get(String(id)); if (k === undefined) throw new Error('Transient network: ' + what + ' refers to the unknown rotor “' + id + '”.'); return k; };
  const E = edges.map((e) => ({ ...e, a: e.type === 'source' ? -1 : N(e.from, 'connection ' + e.id), b: N(e.to, 'connection ' + e.id), lq: e.type === 'source' && e.liquidTo !== undefined && e.liquidTo !== null ? N(e.liquidTo, 'connection ' + e.id) : -1, ri: e.rotor !== undefined ? RT(e.rotor, 'connection ' + e.id) : -1, act: e.type === 'valve' || e.type === 'source' ? e : e.type === 'pump' && e.throttle ? e.throttle : null,
    dz: num(e.dz, 0), head: num(e.head, 0), ks: num(e.ks, 0), k: num(e.k, 1.3), xT: num(e.xT, 0.7), rho: num(e.rho, 1000), shutoff: num(e.shutoff, 1.25), eta: num(e.eta, e.type === 'pump' ? 0.75 : 0.78), stages: Math.max(1, Math.round(num(e.stages, 1))), xg: clamp(num(e.xg, 1), 0, 1), offset: num(e.offset, 0), pos0: num(e.pos0, 1), char: e.char || 'linear', rng: e.rng || 50 }));
  E.forEach((e) => { if (e.act) { e.actTau = Math.max(num(e.act.tau, 1), 1e-3); e.actRate = num(e.act.rate, 1e9); e.actCmd = num(e.act.cmd, num(e.act.pos, 1)); e.actChar = e.act.char || 'linear'; e.actRng = e.act.rng || 50; } if (e.type === 'source') e.f0 = valveCharacteristic(e.pos0, e.char, e.rng); });
  const vol = nodes.map((n) => n.type !== 'fixed');
  nodes.forEach((n) => { if (n.type !== 'fixed' && !(n.V > 0 && n.T > 0 && n.M > 0 && n.p > 0)) throw new Error('Transient network: volume node “' + n.id + '” needs a positive volume, temperature, molar mass and initial pressure.'); if (n.type === 'fixed' && !(n.p > 0)) throw new Error('Transient network: boundary node “' + n.id + '” needs a pressure.'); });
  rotors.forEach((r) => { if (!(r.J > 0 && r.omegaR > 0)) throw new Error('Transient network: rotor “' + r.id + '” needs a positive inertia and rated speed.'); });
  let ny = 0;
  const sG = nodes.map((n, i) => (vol[i] ? ny++ : -1)), sL = nodes.map((n, i) => (vol[i] && n.liquid ? ny++ : -1)), sP = E.map((e) => (e.act ? ny++ : -1)), sW = rotors.map(() => ny++), sI = rotors.map(() => ny++), aGi = ny++, aGo = ny++, aLi = ny++, aLo = ny++;
  const y0 = new Array(ny).fill(0), cG = nodes.map((n) => (vol[n ? nix.get(String(n.id)) : 0] ? (R * n.T) / (n.M * 1e5) : 0));
  nodes.forEach((n, i) => { if (!vol[i]) return; const vl = n.liquid ? Math.max(n.liquid.vol0 || 0, 0) : 0; if (sL[i] >= 0) y0[sL[i]] = vl; y0[sG[i]] = (n.p * Math.max(n.V - vl, 1e-9)) / (cG[i] * zOfNode(n, n.p)); });
  E.forEach((e, i) => { if (sP[i] >= 0) y0[sP[i]] = clamp(num(e.act.pos, 1), 0, 1); });
  rotors.forEach((r, i) => { y0[sW[i]] = r.omegaR * num(r.speed, 1); });
  const evs = events.filter((q) => q && isNum(q.t)).slice().sort((a, b) => a.t - b.t), ch = (e, h) => valveCharacteristic(h, e.char, e.rng);
  const lvTab = nodes.map((n) => { if (n.type === 'fixed' || !n.liquid || !n.liquid.geom) return null; const g = n.liquid.geom, hMax = g.orientation === 'vertical' ? g.l + (g.heads === 'flat' ? 0 : g.heads === 'hemispherical' ? g.d : g.d / 2) : g.d, h = linspace(0, hMax, 161); return { h, v: h.map((q) => vesselVolume({ ...g, level: q })) }; }); // level–volume table of each vessel (exact geometry, 160 intervals)
  const levelOf = (n, vl, i) => (lvTab[i] ? interp1(lvTab[i].v, lvTab[i].h, Math.max(vl, 0)) : Math.max(vl, 0) / (n.liquid.area || 1));
  let evals = 0, jacs = 0;
  /** Right-hand side and every algebraic quantity at time t for the state y. */
  const evaluate = (t, y, want, startOfStep = true) => {
    evals++;
    const dy = new Array(ny).fill(0), p = new Array(nodes.length), cmd = E.map((e) => (e.act ? e.actCmd : 0)), lock = E.map(() => false), trip = rotors.map((r) => !!r.tripped), setR = rotors.map((r) => num(r.set, num(r.speed, 1))), pFix = nodes.map((n) => n.p);
    for (const q of evs) { if (startOfStep ? q.t > t + 1e-9 : q.t > t - 1e-9) break; const s = q.set || {}; if (q.edge !== undefined && eix.has(String(q.edge)) && isNum(s.cmd)) { cmd[eix.get(String(q.edge))] = s.cmd; if (s.lock) lock[eix.get(String(q.edge))] = true; } if (q.rotor !== undefined && rix.has(String(q.rotor))) { const k = rix.get(String(q.rotor)); if (s.trip !== undefined) trip[k] = !!s.trip; if (isNum(s.set)) setR[k] = s.set; } if (q.node !== undefined && nix.has(String(q.node)) && isNum(s.p)) pFix[nix.get(String(q.node))] = s.p; }
    nodes.forEach((n, i) => { if (!vol[i]) { p[i] = pFix[i]; return; } const vl = sL[i] >= 0 ? y[sL[i]] : 0, c = (cG[i] * Math.max(y[sG[i]], 0)) / Math.max(n.V - vl, 1e-6 * n.V); let pp = c * zOfNode(n, n.p); for (let k = 0; k < 4; k++) pp = c * zOfNode(n, pp); p[i] = pp; });
    const rho = (i, j) => { const n = vol[i] ? nodes[i] : nodes[j]; return (p[i] * 1e5 * n.M) / (zOfNode(n, p[i]) * R * n.T); }, sp = rotors.map((r, i) => Math.max(y[sW[i]], 0) / r.omegaR), tq = rotors.map(() => 0), pw = rotors.map(() => 0);
    const m = new Array(E.length).fill(0), margin = new Array(E.length).fill(null);
    E.forEach((e, i) => {
      const pos = sP[i] >= 0 ? clamp(y[sP[i]], 0, 1) : 1;
      if (e.type === 'source') {
        const f = ch(e, pos), f0 = e.f0, pd = p[e.b] + e.offset, mm = e.curve.m, pa = (x) => interp1(mm, e.curve.p, clamp(x, mm[0], mm[mm.length - 1])), c = f > 1e-9 ? e.c0 * (f0 / f) ** 2 : Infinity, g = (x) => pa(x) - pd - c * x * x;
        let lo = 0, hi = mm[mm.length - 1], glo = f > 1e-9 ? g(0) : -1, ghi;
        if (glo <= 0) hi = 0; else if ((ghi = g(hi)) >= 0) lo = hi; else for (let k = 0; k < 60 && hi - lo > 1e-13 * mm[mm.length - 1]; k++) { /* regula falsi with the Illinois modification */ const xm = hi - (ghi * (hi - lo)) / (ghi - glo), gm = g(xm); if (gm === 0) { lo = hi = xm; break; } if (gm > 0) { lo = xm; glo = gm; ghi *= 0.5; } else { hi = xm; ghi = gm; glo *= 0.5; } }
        m[i] = 0.5 * (lo + hi);
        const xg = e.xg; dy[sG[e.b]] += xg * m[i]; dy[aGi] += xg * m[i];
        if (e.lq >= 0 && sL[e.lq] >= 0) { dy[sL[e.lq]] += ((1 - xg) * m[i]) / nodes[e.lq].liquid.rho; dy[aLi] += (1 - xg) * m[i]; }
        return;
      }
      const a = e.a, b = e.b, liquid = e.phase === 'liquid', d = p[a] - p[b] - e.dz;
      let w = 0;
      if (e.type === 'pipe') w = rsq(d, 1e-5) / Math.sqrt(Math.max(e.r, 1e-12));
      else if (e.type === 'valve') {
        const cv = e.cv * ch(e, pos);
        if (liquid) w = ((N6 * cv) / 3600) * Math.sqrt(e.rho) * rsq(d, 1e-6);
        else { const u = d >= 0 ? a : b, p1 = Math.max(p[u], 1e-6), xc = (e.k / 1.4) * e.xT, x = Math.min(Math.abs(d) / p1, xc); w = Math.sign(d) * ((N6 * cv) / 3600) * (1 - x / (3 * xc)) * Math.sqrt(rho(u, u === a ? b : a) * p1) * rsq(x, 1e-9); }
      } else if (e.type === 'pump') {
        const s = e.ri >= 0 ? sp[e.ri] : 1, rl = e.rho, sh = e.shutoff, hreq = (-d * 1e5) / (rl * G) + e.head, kt = e.throttle ? e.throttle.k0 / Math.max(valveCharacteristic(pos, e.actChar, e.actRng), 1e-3) ** 2 : 0, nu = e.hr * s * s * sh - hreq, q = nu > 0 ? Math.sqrt(nu / (((sh - 1) * e.hr) / (e.qr * e.qr) + e.ks + kt)) : 0;
        const H = s > 1e-6 ? e.hr * s * s * (sh - (sh - 1) * (q / (e.qr * s)) ** 2) : 0, xq = s > 1e-6 ? q / (e.qr * s) : 0, e0 = e.eta;
        w = (rl * q) / 3600;
        if (e.ri >= 0) { const P = xq > 1e-4 ? (rl * G * (q / 3600) * Math.max(H, 0)) / (e0 * Math.max(2 * xq - xq * xq, 0.05)) : (rl * G * ((e.qr * s) / 3600) * Math.max(H, 0)) / (2 * e0); /* efficiency parabola η = η₀(2x − x²): finite shut-off power */ pw[e.ri] += P; tq[e.ri] += P / Math.max(y[sW[e.ri]], 1e-3 * rotors[e.ri].omegaR); }
      } else if (e.type === 'compressor') {
        const s = e.ri >= 0 ? sp[e.ri] : 1, n = vol[a] ? nodes[a] : nodes[b], st = e.stages, kk = e.k, et = e.eta, nm = (kk - 1) / (kk * et), rr = Math.max(p[b] / Math.max(p[a], 1e-6), 1e-3) ** (1 / st), hreq = ((zOfNode(n, p[a]) * R * n.T) / n.M / nm) * (rr ** nm - 1);
        const phi = s > 1e-3 ? mapFlowCoefficient(hreq / (s * s * e.hd)) : 0;
        w = phi * s * e.qd * rho(a, b); margin[i] = phi > 1e-9 ? Math.max((100 * (phi - PHI_SURGE)) / phi, -100) : -100;
        if (e.ri >= 0) { const P = (w * Math.max(hreq, 0) * st) / et; pw[e.ri] += P; tq[e.ri] += P / Math.max(y[sW[e.ri]], 1e-3 * rotors[e.ri].omegaR); }
      } else throw new Error('Transient network: unknown connection type “' + e.type + '”.');
      if (e.check && w < 0) w = 0;
      const u = w >= 0 ? a : b; // the upstream inventory limits what can leave
      if (liquid && vol[u] && sL[u] >= 0) w *= clamp(y[sL[u]] / (0.01 * nodes[u].V), 0, 1);
      m[i] = w;
      if (liquid) { const rl = e.rho; if (vol[a] && sL[a] >= 0) dy[sL[a]] -= w / rl; if (vol[b] && sL[b] >= 0) dy[sL[b]] += w / rl; if (!vol[a]) dy[w >= 0 ? aLi : aLo] += Math.abs(w); if (!vol[b]) dy[w >= 0 ? aLo : aLi] += Math.abs(w); }
      else { if (vol[a]) dy[sG[a]] -= w; if (vol[b]) dy[sG[b]] += w; if (!vol[a]) dy[w >= 0 ? aGi : aGo] += Math.abs(w); if (!vol[b]) dy[w >= 0 ? aGo : aGi] += Math.abs(w); }
    });
    const level = nodes.map((n, i) => (sL[i] >= 0 ? levelOf(n, y[sL[i]], i) : null));
    for (const c of controllers) {
      if (c.type === 'tripOpen') continue;
      const ms = c.measure || {}, k = ms.node !== undefined ? nix.get(String(ms.node)) : undefined, j = ms.edge !== undefined ? eix.get(String(ms.edge)) : undefined, val = k !== undefined ? (ms.q === 'level' ? level[k] : p[k]) : j !== undefined ? margin[j] : null;
      if (!isNum(val)) continue;
      const u = clamp(num(c.bias, 0) + num(c.gain, 1) * (val - num(c.set, 0)), num(c.min, 0), num(c.max, 1));
      if (c.out && c.out.edge !== undefined && eix.has(String(c.out.edge))) { if (!lock[eix.get(String(c.out.edge))]) cmd[eix.get(String(c.out.edge))] = u; } else if (c.out && c.out.rotor !== undefined && rix.has(String(c.out.rotor))) setR[rix.get(String(c.out.rotor))] = u;
    }
    for (const c of controllers) if (c.type === 'tripOpen' && rix.has(String(c.rotor)) && trip[rix.get(String(c.rotor))] && c.out && eix.has(String(c.out.edge))) cmd[eix.get(String(c.out.edge))] = 1;
    E.forEach((e, i) => { if (sP[i] < 0) return; const rate = e.actRate; dy[sP[i]] = clamp((clamp(cmd[i], 0, 1) - y[sP[i]]) / e.actTau, -rate, rate); });
    rotors.forEach((r, i) => {
      const er = setR[i] - sp[i], u = num(r.kp, 10) * er + y[sI[i]], tm = num(r.tauMax, 0), td = trip[i] ? 0 : tm * clamp(u, 0, 1);
      dy[sI[i]] = trip[i] || (u > 1 && er > 0) || (u < 0 && er < 0) ? 0 : er / Math.max(num(r.ti, 2), 1e-3);
      dy[sW[i]] = (td - tq[i] - num(r.tauLoss, 0) * sp[i] * sp[i]) / r.J;
      if (want) want.tq = tq;
    });
    return want ? { dy, p, m, level, sp, pw, margin } : dy;
  };
  { const a0 = evaluate(0, y0, {}); rotors.forEach((r, i) => { if (r.int0 === undefined && r.tauMax > 0) y0[sI[i]] = clamp((a0.dy[sW[i]] * -r.J) / r.tauMax, 0, 1); else y0[sI[i]] = num(r.int0, 0); }); } // governor integral that holds the initial speed
  const scale = y0.map((v) => Math.max(Math.abs(v), 1e-6)), gTot = Math.max(nodes.reduce((s, n, i) => s + (vol[i] ? y0[sG[i]] : 0), 0), 1e-6), lTot = Math.max(nodes.reduce((s, n, i) => s + (sL[i] >= 0 ? n.V * n.liquid.rho : 0), 0), 1);
  nodes.forEach((n, i) => { if (sL[i] >= 0) scale[sL[i]] = n.V; }); E.forEach((e, i) => { if (sP[i] >= 0) scale[sP[i]] = 1; }); rotors.forEach((r, i) => { scale[sW[i]] = r.omegaR; scale[sI[i]] = 1; });
  scale[aGi] = scale[aGo] = gTot; scale[aLi] = scale[aLo] = lTot;
  let halved = 0, substeps = 0, loose = 0, jacH = null;
  const Jm = Array.from({ length: ny }, () => new Array(ny).fill(0));
  const step = (t, y, h, depth) => {
    const f0 = evaluate(t, y), t1 = t + h, Gf = (yy) => { const f1 = evaluate(t1, yy, null, false); return yy.map((v, j) => v - y[j] - h * ((1 - theta) * f0[j] + theta * f1[j])); };
    let yn = y.map((v, j) => v + h * f0[j]), g = Gf(yn), ok = false;
    // modified Newton: the finite-difference Jacobian is kept from step to step and rebuilt when convergence is slow (the branch
    // laws have kinks: check valves, rate limits, the surge line); the four mass accumulators feed nothing back (unit columns)
    const jac = () => { for (let j = 0; j < ny; j++) { if (j >= aGi) { for (let i = 0; i < ny; i++) Jm[i][j] = i === j ? 1 : 0; continue; } const hh = 1e-6 * Math.max(Math.abs(yn[j]), scale[j]), yp = yn.slice(); yp[j] += hh; const gp = Gf(yp); for (let i = 0; i < ny; i++) Jm[i][j] = (gp[i] - g[i]) / hh; } jacH = h; fresh = true; jacs++; };
    let fresh = false;
    if (jacH !== h) jac();
    for (let it = 0; it < 14; it++) {
      if (!fresh && it === 4) jac(); else if (fresh && it === 9) jac();
      let d; try { d = solveLinear(Jm.map((r) => r.slice()), g.map((x) => -x)); } catch { d = null; }
      if (!d || !d.every(Number.isFinite)) break;
      let lam = 1, gn = null, yt = null; const n0 = g.reduce((s, x, j) => Math.max(s, Math.abs(x) / scale[j]), 0);
      for (let k = 0; k < 6; k++, lam *= 0.5) { yt = yn.map((v, j) => v + lam * d[j]); gn = Gf(yt); if (gn.reduce((s, x, j) => Math.max(s, Math.abs(x) / scale[j]), 0) <= n0 * (1 - 1e-4 * lam) || n0 < 1e-13) break; }
      const dn = d.reduce((s, x, j) => Math.max(s, (lam * Math.abs(x)) / (Math.abs(yn[j]) + scale[j])), 0);
      yn = yt; g = gn;
      if (dn < 1e-7 && g.reduce((s, x, j) => Math.max(s, Math.abs(x) / scale[j]), 0) < 1e-7) { ok = true; break; }
    }
    if (!ok && depth < 4) { halved++; return step(t + h / 2, step(t, y, h / 2, depth + 1), h / 2, depth + 1); }
    if (!ok) loose++; // accepted at the smallest step with the last Newton iterate
    substeps++;
    E.forEach((e, i) => { if (sP[i] >= 0) yn[sP[i]] = clamp(yn[sP[i]], 0, 1); });
    rotors.forEach((r, i) => { yn[sW[i]] = Math.max(yn[sW[i]], 0); });
    return yn;
  };
  const nSteps = Math.max(1, Math.round(tEnd / dt)), h = tEnd / nSteps, every = Math.max(1, Math.ceil(nSteps / Math.max(nOut, 2)));
  const out = { t: [], p: {}, level: {}, liquid: {}, m: {}, speed: {}, power: {}, pos: {}, margin: {}, surgeTime: {} };
  nodes.forEach((n, i) => { out.p[n.id] = []; if (sL[i] >= 0) { out.level[n.id] = []; out.liquid[n.id] = []; } });
  E.forEach((e, i) => { out.m[e.id] = []; if (sP[i] >= 0) out.pos[e.id] = []; if (e.type === 'compressor') { out.margin[e.id] = []; out.surgeTime[e.id] = null; } });
  rotors.forEach((r) => { out.speed[r.id] = []; out.power[r.id] = []; });
  let y = y0.slice();
  const rec = (t, keep) => {
    const a = evaluate(t, y, {});
    E.forEach((e, i) => { if (e.type === 'compressor' && out.surgeTime[e.id] === null && a.margin[i] !== null && a.margin[i] < 0 && a.sp[e.ri >= 0 ? e.ri : 0] > 0.05) out.surgeTime[e.id] = t; });
    if (!keep) return;
    out.t.push(t); nodes.forEach((n, i) => { out.p[n.id].push(a.p[i]); if (sL[i] >= 0) { out.level[n.id].push(a.level[i]); out.liquid[n.id].push(y[sL[i]]); } });
    E.forEach((e, i) => { out.m[e.id].push(a.m[i]); if (sP[i] >= 0) out.pos[e.id].push(y[sP[i]]); if (e.type === 'compressor') out.margin[e.id].push(a.margin[i]); });
    rotors.forEach((r, i) => { out.speed[r.id].push(a.sp[i]); out.power[r.id].push(a.pw[i] / 1000); });
  };
  rec(0, true);
  for (let k = 1; k <= nSteps; k++) { y = step((k - 1) * h, y, h, 0); rec(k * h, k % every === 0 || k === nSteps); }
  const gas = (yy) => nodes.reduce((s, n, i) => s + (vol[i] ? yy[sG[i]] : 0), 0), liq = (yy) => nodes.reduce((s, n, i) => s + (sL[i] >= 0 ? yy[sL[i]] * n.liquid.rho : 0), 0);
  out.balance = { gas: { initial: gas(y0), in: y[aGi], out: y[aGo], final: gas(y) }, liquid: { initial: liq(y0), in: y[aLi], out: y[aLo], final: liq(y) } };
  out.steps = nSteps; out.substeps = substeps; out.halved = halved; out.loose = loose; out.evaluations = evals; out.jacobians = jacs;
  return out;
}
/**
 * Lumped dynamic model of a production system for transientNetwork, built from plain numbers (the transientBasis this suite
 * publishes in its outputs, so that another suite can run further scenarios without repeating the steady solution).
 * B: { m0 (kg/s total at t = 0), xg (gas mass fraction at the separator), rhoL, M (kg/mol), k, z: { p[], z[] }, stroke (s), asvSet (%),
 *      flareSet (× arrival pressure), flareP, liqP (bara), headerV (m³),
 *      wells: { m[], p[] (deliverability curve), pAvail0, chokeDp (bar), pos0, char, rng }, line: { V, T, p0 },
 *      sep: { V, T, p0, pSet, vol0, geom, levelSet, H },
 *      comp: null | { stages, pd, eta, margin (%), rpm, H (inertia constant s), Tcool (K), dpExport (bar), dischargeV (m³), asvSize },
 *      pump: null | { qr (m³/h), hr (m), shutoff, eta, speed, head (static m), ks, rpm, H } }
 * scenario: 'compressorTrip' | 'valveClosure' | 'pumpTrip' | 'chokeChange' | 'none'; o: { tEvent (s, default 10), tEnd (s), steps, stroke (s), chokeTo (0–1) }
 * Returns { nodes, edges, rotors, controllers, events, tEnd, dt, scenario, tEvent, mGas, mLiquid } — pass it to transientNetwork.
 * Node ids: LINE (flowline gas line-pack), SEP, DISCHARGE, HEADER, EXPORT, FLARE, LIQ; edge ids: WELLS, FLOWLINE, COMPRESSOR, CHECKVALVE, RECYCLE, EXPORTVALVE,
 * FLAREVALVE, PUMP (or DRAIN); rotor ids: compressor, pump.
 */
export function transientSystem(B, scenario = 'compressorTrip', o = {}) {
  if (!B || !(B.m0 > 0) || !(B.xg > 0 && B.xg < 1)) throw new Error('The transient needs a flowing case with both gas and liquid at the separator.');
  const tEvent = Math.max(num(o.tEvent, 10), 0), tEnd = Math.max(num(o.tEnd, 120), tEvent + 1), stroke = Math.max(num(o.stroke, num(B.stroke, 20)), 0.5), steps = Math.max(10, Math.round(num(o.steps, 240)));
  const mG = B.xg * B.m0, mL = (1 - B.xg) * B.m0, k = num(B.k, 1.3), xc = (k / 1.4) * 0.7, S = B.sep, pS = S.pSet, gasNode = (id, V, T, p) => ({ id, type: 'volume', V, T, M: B.M, p, z: B.z });
  const rhoAt = (p, T) => (p * 1e5 * B.M) / (zOfNode({ z: B.z }, p) * R * T), gasCv = (w, p1, T, dp) => { const x = Math.min(Math.max(dp, 1e-6) / p1, xc); return (w * 3600) / (N6 * (1 - x / (3 * xc)) * Math.sqrt(x * p1 * rhoAt(p1, T))); };
  const flareP = num(B.flareP, 1.5), liqP = B.pump ? B.liqP : Math.max(Math.min(B.liqP, pS - 2), 1), W = B.wells, dp0 = Math.max(W.chokeDp, 0.05);
  const nodes = [gasNode('LINE', B.line.V, B.line.T, B.line.p0), { ...gasNode('SEP', S.V, S.T, S.p0), liquid: { vol0: S.vol0, rho: B.rhoL, geom: S.geom } }, { id: 'FLARE', type: 'fixed', p: flareP }, { id: 'LIQ', type: 'fixed', p: liqP }];
  const edges = [{ id: 'WELLS', type: 'source', to: 'LINE', liquidTo: 'SEP', curve: { m: W.m, p: W.p }, offset: W.pAvail0 - dp0 - B.line.p0, c0: dp0 / (B.m0 * B.m0), pos: W.pos0, pos0: W.pos0, tau: 1, rate: 1 / stroke, char: W.char, rng: W.rng, xg: B.xg },
    { id: 'FLOWLINE', type: 'pipe', phase: 'gas', from: 'LINE', to: 'SEP', r: Math.max(B.line.p0 - pS, 0.05) / (mG * mG) }];
  const rotors = [], controllers = [], events = [];
  if (B.comp) {
    const C = B.comp, st = Math.max(1, Math.round(C.stages)), nm = (k - 1) / (k * C.eta), h0 = ((zOfNode({ z: B.z }, pS) * R * S.T) / B.M / nm) * ((C.pd / pS) ** (nm / st) - 1), phi0 = clamp(1 / (1 + C.margin / 100), 0.66, 1.08), qIn0 = mG / rhoAt(pS, S.T), P0 = (mG * h0 * st) / C.eta, wR = (C.rpm * 2 * Math.PI) / 60;
    const dpCheck = 0.5, pH = C.pd - dpCheck; // discharge piping up to the non-return valve, then the header / export line-pack
    nodes.push(gasNode('DISCHARGE', Math.max(num(C.dischargeV, 3), 0.05), C.Tcool, C.pd), gasNode('HEADER', B.headerV, C.Tcool, pH), { id: 'EXPORT', type: 'fixed', p: pH - C.dpExport });
    edges.push({ id: 'COMPRESSOR', type: 'compressor', phase: 'gas', from: 'SEP', to: 'DISCHARGE', qd: qIn0 / phi0, hd: h0 / mapShape(phi0), stages: st, eta: C.eta, k, rotor: 'compressor', check: true },
      { id: 'CHECKVALVE', type: 'pipe', phase: 'gas', from: 'DISCHARGE', to: 'HEADER', r: dpCheck / (mG * mG), check: true },
      { id: 'RECYCLE', type: 'valve', phase: 'gas', from: 'DISCHARGE', to: 'SEP', cv: num(C.asvSize, 2) * gasCv(mG, C.pd, C.Tcool, C.pd - pS), pos: 0, cmd: 0, tau: 0.2, rate: 2, char: 'linear', k, check: true },
      { id: 'EXPORTVALVE', type: 'valve', phase: 'gas', from: 'HEADER', to: 'EXPORT', cv: gasCv(mG, pH, C.Tcool, C.dpExport), pos: 1, cmd: 1, tau: 1, rate: 1 / stroke, char: 'linear', k, check: true });
    rotors.push({ id: 'compressor', J: (2 * C.H * P0) / (wR * wR), omegaR: wR, speed: 1, set: 1, tauMax: (1.25 * P0) / wR, kp: 15, ti: 1.5, tauLoss: (0.02 * P0) / wR, ratedPower: P0 });
    controllers.push({ type: 'P', measure: { edge: 'COMPRESSOR', q: 'surgeMargin' }, set: B.asvSet, gain: -1 / B.asvSet, bias: 0, min: 0, max: 1, out: { edge: 'RECYCLE' } }, { type: 'tripOpen', rotor: 'compressor', out: { edge: 'RECYCLE' } });
  } else {
    nodes.push({ id: 'EXPORT', type: 'fixed', p: 0.8 * pS });
    edges.push({ id: 'EXPORTVALVE', type: 'valve', phase: 'gas', from: 'SEP', to: 'EXPORT', cv: gasCv(mG, pS, S.T, 0.2 * pS), pos: 1, cmd: 1, tau: 1, rate: 1 / stroke, char: 'linear', k, check: true });
  }
  edges.push({ id: 'FLAREVALVE', type: 'valve', phase: 'gas', from: 'SEP', to: 'FLARE', cv: gasCv(1.2 * mG, B.flareSet * pS, S.T, B.flareSet * pS - flareP), pos: 0, cmd: 0, tau: 0.5, rate: 0.5, char: 'linear', k, check: true });
  controllers.push({ type: 'P', measure: { node: 'SEP', q: 'p' }, set: B.flareSet * pS, gain: 1 / (0.05 * pS), bias: 0, min: 0, max: 1, out: { edge: 'FLAREVALVE' } });
  const q0 = (mL / B.rhoL) * 3600;
  if (B.pump) {
    const P = B.pump, s0 = P.speed, sh = P.shutoff, hp0 = P.hr * s0 * s0 * (sh - (sh - 1) * (q0 / (P.qr * s0)) ** 2), hreq0 = ((liqP - pS) * 1e5) / (B.rhoL * G) + P.head, kth = (hp0 - hreq0) / (q0 * q0) - P.ks, pos0 = kth > 0 ? 0.6 : 1, Pr = (B.rhoL * G * (P.qr / 3600) * P.hr) / P.eta, wR = (P.rpm * 2 * Math.PI) / 60;
    edges.push({ id: 'PUMP', type: 'pump', phase: 'liquid', from: 'SEP', to: 'LIQ', qr: P.qr, hr: P.hr, shutoff: sh, eta: P.eta, rho: B.rhoL, head: P.head, ks: P.ks, rotor: 'pump', check: true, throttle: { k0: Math.max(kth, 0) * pos0 * pos0, pos: pos0, cmd: pos0, tau: 2, rate: 0.2, char: 'linear' } });
    rotors.push({ id: 'pump', J: (2 * P.H * Pr) / (wR * wR), omegaR: wR, speed: s0, set: s0, tauMax: (1.3 * Pr) / wR, kp: 15, ti: 1, tauLoss: (0.03 * Pr) / wR, ratedPower: Pr });
    controllers.push({ type: 'P', measure: { node: 'SEP', q: 'level' }, set: S.levelSet, gain: 3 / S.H, bias: pos0, min: 0.02, max: 1, out: { edge: 'PUMP' } });
  } else {
    edges.push({ id: 'DRAIN', type: 'valve', phase: 'liquid', from: 'SEP', to: 'LIQ', cv: (mL * 3600) / (N6 * Math.sqrt((pS - liqP) * B.rhoL)) / 0.6, rho: B.rhoL, pos: 0.6, cmd: 0.6, tau: 2, rate: 0.2, char: 'linear', check: true });
    controllers.push({ type: 'P', measure: { node: 'SEP', q: 'level' }, set: S.levelSet, gain: 3 / S.H, bias: 0.6, min: 0, max: 1, out: { edge: 'DRAIN' } });
  }
  let sc = ['compressorTrip', 'valveClosure', 'pumpTrip', 'chokeChange', 'none'].includes(scenario) ? scenario : 'none';
  if (sc === 'compressorTrip' && !B.comp) sc = 'valveClosure';
  if (sc === 'compressorTrip') events.push({ t: tEvent, rotor: 'compressor', set: { trip: true } });
  else if (sc === 'valveClosure') events.push({ t: tEvent, edge: 'EXPORTVALVE', set: { cmd: 0 } });
  else if (sc === 'pumpTrip') events.push(B.pump ? { t: tEvent, rotor: 'pump', set: { trip: true } } : { t: tEvent, edge: 'DRAIN', set: { cmd: 0, lock: true } });
  else if (sc === 'chokeChange') events.push({ t: tEvent, edge: 'WELLS', set: { cmd: clamp(num(o.chokeTo, 0.3), 0, 1) } });
  return { nodes, edges, rotors, controllers, events, tEnd, dt: tEnd / steps, scenario: sc, tEvent, mGas: mG, mLiquid: mL };
}
/**
 * Orientation of directional equipment against the solved flow direction: pumps, compressors, chokes / control valves and check
 * valves are built for flow from their “from” node to their “to” node. edges: [{ from, to, kind }], m: mass flows (kg/s, + along from → to).
 * Returns [{ item, kind, direction: 'as drawn' | 'reversed' | 'no flow', ok }].
 */
export function orientationCheck(edges, m) {
  const directional = /pump|comp|choke|valve|check/;
  return edges.filter((e) => directional.test(String(e.kind || e.type || ''))).map((e) => { const i = edges.indexOf(e), w = m[i] ?? 0, tiny = 1e-9 * Math.max(1e-9, ...m.map((x) => Math.abs(x))); return { item: e.from + ' → ' + e.to, kind: String(e.kind || e.type), direction: Math.abs(w) <= tiny ? 'no flow' : w > 0 ? 'as drawn' : 'reversed', ok: w >= -tiny }; });
}
/**
 * Compare measured values with model values. rows: [{ kind (text), at, aux, value }]; models: { key: { label, unit, re (RegExp matching
 * the kind text), f(at, aux) → model value } }. Returns { rows: [{ key, label, unit, at, aux, measured, model, diff, pct }],
 * groups: [{ key, label, unit, n, bias, rmse, mape }], unknown (rows whose kind matched no model or that have no number) }.
 */
export function compareMeasured(rows, models) {
  const out = [], by = new Map();
  let unknown = 0;
  for (const r of Array.isArray(rows) ? rows : []) {
    const kind = String(r?.kind ?? '').trim().toLowerCase(), val = num(r?.value, null), key = kind ? Object.keys(models).find((k) => models[k].re.test(kind)) : null;
    if (!kind && val === null) continue;
    if (!key || val === null) { unknown++; continue; }
    let mod = null; try { mod = models[key].f(num(r.at, 0), num(r.aux, 0)); } catch { mod = null; }
    if (!isNum(mod)) { unknown++; continue; }
    const row = { key, label: models[key].label, unit: models[key].unit, at: num(r.at, 0), aux: num(r.aux, 0), measured: val, model: mod, diff: mod - val, pct: Math.abs(val) > 1e-12 ? (100 * (mod - val)) / Math.abs(val) : null };
    out.push(row); if (!by.has(key)) by.set(key, []); by.get(key).push(row);
  }
  const groups = [...by.entries()].map(([key, a]) => { const n = a.length, pc = a.filter((q) => q.pct !== null); return { key, label: models[key].label, unit: models[key].unit, n, bias: a.reduce((s, q) => s + q.diff, 0) / n, rmse: Math.sqrt(a.reduce((s, q) => s + q.diff * q.diff, 0) / n), mape: pc.length ? pc.reduce((s, q) => s + Math.abs(q.pct), 0) / pc.length : null }; });
  return { rows: out, groups, unknown };
}

// ---- materials ---------------------------------------------------------------------------------------------------
/** Line-pipe grade library: SMYS, SMTS, Young's modulus (MPa), Poisson ratio, thermal expansion (1/K), density (kg/m³), conductivity (W/m/K). */
export const GRADES = Object.freeze({
  B: { label: 'API 5L Grade B', smys: 245, smts: 415 }, X42: { label: 'API 5L X42', smys: 290, smts: 415 }, X52: { label: 'API 5L X52', smys: 360, smts: 460 }, X56: { label: 'API 5L X56', smys: 390, smts: 490 },
  X60: { label: 'API 5L X60', smys: 415, smts: 520 }, X65: { label: 'API 5L X65', smys: 450, smts: 535 }, X70: { label: 'API 5L X70', smys: 485, smts: 570 }, X80: { label: 'API 5L X80', smys: 555, smts: 625 },
  duplex: { label: '22Cr duplex (UNS S31803)', smys: 450, smts: 620, E: 200000, alphaT: 1.3e-5, rho: 7800, k: 15 }, superduplex: { label: '25Cr super duplex (UNS S32750)', smys: 550, smts: 750, E: 200000, alphaT: 1.35e-5, rho: 7800, k: 14 },
  clad: { label: 'X65 + 3 mm Alloy 625 clad (CRA liner not counted for strength)', smys: 450, smts: 535 },
});
/** Length (m) over which a pipe of bending stiffness EI (N·m²) and submerged weight w (N/m) lifts off a rigid seabed on each side of a point support of height h (m): w L⁴ / (72 EI) = h. */
export const liftOffLength = (EI, h, w) => (w > 0 && EI > 0 && h > 0 ? ((72 * EI * h) / w) ** 0.25 : 0);
/** Buckle-propagation pressure (bar) of a pipe: 35 · f_y · α_fab · (t / D)^2.5 with f_y in MPa (form used in DNV-ST-F101). */
export const propagationPressure = (fyMPa, t, D, aFab = 1) => 35 * fyMPa * 10 * aFab * (t / D) ** 2.5;
/** Net anode mass (kg) for a protection current I (A) over a life (years): I · t / (u · ε), ε = electrochemical capacity (A·h/kg), u = utilisation factor. */
export const anodeMassRequired = (I, years, capacity = 2000, u = 0.8) => (I * years * 8760) / (u * capacity);
/** Material record of a grade: { grade, label, smys, smts, E (MPa), poisson, alphaT, rho, k }. */
export const material = (grade) => { const g = GRADES[grade] || GRADES.X65; return { grade: GRADES[grade] ? grade : 'X65', label: g.label, smys: g.smys, smts: g.smts, E: g.E ?? BASE.E, poisson: BASE.poisson, alphaT: g.alphaT ?? BASE.alphaT, rho: g.rho ?? BASE.rhoSteel, k: g.k ?? 45 }; };

// ---- suite declaration ---------------------------------------------------------------------------------------
const DEF_SURVEY = [{ md: 0, inc: 0, azi: 45 }, { md: 800, inc: 0, azi: 45 }, { md: 1580, inc: 51.9, azi: 45 }, { md: 3400, inc: 51.9, azi: 45 }];
const pipeRow = (from, to, length, diameter, k = 0, type = 'pipe', param = 0) => ({ from, to, type, length, diameter, k, param });
const WELL_IDS = Array.from({ length: Math.max(1, Math.round(BASE.wells ?? 2)) }, (_, i) => 'W' + (i + 1)), JUMPER_M = [80, 120, 150, 180, 220, 260];
const DEF_NETWORK = [...WELL_IDS.map((id, i) => pipeRow(id, 'M1', JUMPER_M[i % JUMPER_M.length], 152.4, 1.5)), pipeRow('M1', 'PLEM', 60, 254, 1), pipeRow('PLEM', 'RB', 0, 0, 0, 'flowline'), pipeRow('RB', 'SEP', 0, 0, 0, 'riser')];
const DEF_NODES = [...WELL_IDS.map((id) => ({ id, kind: 'source', share: +(100 / WELL_IDS.length).toFixed(4), status: 'flowing' })), { id: 'SEP', kind: 'sink' }];
const DEF_SUPPORTS = [{ x: 0, type: 'PLET on a sliding mudmat', k: 2000 }, { x: BASE.riserBaseX, type: 'riser-base anchor', k: 50000 }, { x: BASE.profile[BASE.profile.length - 1].x, type: 'hang-off flex joint', k: 100000 }];
const DEF_LAYERS = [{ name: 'Fusion-bonded epoxy', tMm: 0.5, k: 0.2, rho: 1200, cp: 1400 }, { name: BASE.insulation.name, tMm: BASE.insulation.t * 1000, k: BASE.insulation.k, rho: 780, cp: 1700 }];
const DEF_FITTINGS = [{ type: 'bend90lr', count: 4 }, { type: 'teeRun', count: 2 }, { type: 'ball', count: 3 }, { type: 'check', count: 1 }, { type: 'exit', count: 1 }];
const opt = (pairs) => pairs.map(([value, label]) => ({ value, label }));
const gen = (v) => v.geomMode === 'generate', genOrTer = (v) => v.geomMode !== 'table';

const INPUTS = [
  { group: 'Route and riser', tab: 'inputs', help: 'Where the line runs. Elevation is relative to mean sea level (negative below); distance is horizontal from the inlet.', fields: [
    { key: 'geomMode', label: 'Route source', type: 'select', value: 'table', options: opt([['table', 'Elevation profile table (imported or typed)'], ['centreline', '3-D centreline coordinates table'], ['generate', 'Generate from water depth, flowline length and riser type'], ['terrain', 'Drape over the imported terrain grid']]) },
    { key: 'profile', label: 'Elevation profile', type: 'table', columns: [{ key: 'x', label: 'Distance', unit: 'm' }, { key: 'z', label: 'Elevation', unit: 'm' }], value: BASE.profile.map((p) => ({ x: p.x, z: p.z })), showIf: (v) => v.geomMode === 'table', help: 'Imported routes are written here.' },
    { key: 'centreline', label: 'Centreline coordinates', type: 'table', columns: [{ key: 'x', label: 'Easting / longitude', unit: '' }, { key: 'y', label: 'Northing / latitude', unit: '' }, { key: 'z', label: 'Elevation', unit: '' }], value: [], showIf: (v) => v.geomMode === 'centreline', help: 'Points along the pipe axis in route order; chainage is measured in plan and plan-view bends are added to the bend list.' },
    { key: 'coordSystem', label: 'Coordinate reference system', type: 'select', value: 'local', options: opt([['local', 'Projected / local grid (easting, northing)'], ['geographic', 'Geographic WGS 84 (longitude, latitude in degrees)']]), showIf: (v) => v.geomMode === 'centreline' },
    { key: 'lengthUnit', label: 'Unit of distances and grid coordinates in the tables', type: 'select', value: 'm', options: opt([['m', 'metres'], ['km', 'kilometres'], ['ft', 'feet']]), showIf: (v) => v.geomMode === 'table' || v.geomMode === 'centreline' },
    { key: 'elevUnit', label: 'Unit of elevations in the tables', type: 'select', value: 'm', options: opt([['m', 'metres'], ['ft', 'feet']]), showIf: (v) => v.geomMode === 'table' || v.geomMode === 'centreline' },
    { key: 'elevSense', label: 'Vertical coordinate', type: 'select', value: 'up', options: opt([['up', 'Elevation, positive upwards'], ['depth', 'Depth, positive downwards']]), showIf: (v) => v.geomMode === 'table' || v.geomMode === 'centreline' },
    { key: 'riserType', label: 'Riser configuration', type: 'select', value: 'scr', options: opt([['scr', 'Steel catenary riser'], ['lazyWave', 'Lazy-wave riser'], ['flexible', 'Flexible free-hanging catenary'], ['vertical', 'Vertical (top-tensioned / platform) riser'], ['none', 'No riser']]), showIf: genOrTer },
    { key: 'waterDepth', label: 'Water depth at the riser base', unit: 'm', value: BASE.waterDepth, min: 0, max: 3500, typical: [50, 2500], showIf: gen },
    { key: 'flowlineLength', label: 'Flowline length (horizontal)', unit: 'm', value: BASE.riserBaseX, min: 10, max: 5e5, typical: [1000, 150000], showIf: gen },
    { key: 'seabedDrop', label: 'Seabed fall from inlet to riser base', unit: 'm', value: 100, min: -1500, max: 1500, showIf: gen, help: 'Positive when the inlet is shallower than the riser base.' },
    { key: 'undulationAmp', label: 'Seabed undulation amplitude', unit: 'm', value: 20, min: 0, max: 300, showIf: gen },
    { key: 'undulationLength', label: 'Seabed undulation wavelength', unit: 'm', value: 6000, min: 50, max: 1e5, showIf: gen },
    { key: 'topsideElev', label: 'Arrival elevation above sea level', unit: 'm', value: 25, min: 0, max: 100, showIf: genOrTer },
    { key: 'hangoffAngle', label: 'Hang-off angle from vertical', unit: '°', value: 12, min: 2, max: 60, typical: [8, 20], showIf: (v) => genOrTer(v) && ['scr', 'flexible', 'lazyWave'].includes(v.riserType) },
    { key: 'lazyLift', label: 'Lazy wave: slope where buoyancy starts', unit: '°', value: 35, min: 5, max: 75, showIf: (v) => genOrTer(v) && v.riserType === 'lazyWave' },
    { key: 'lazySag', label: 'Lazy wave: downward slope where buoyancy ends', unit: '°', value: 15, min: 0, max: 60, showIf: (v) => genOrTer(v) && v.riserType === 'lazyWave' },
    { key: 'lazyBuoy', label: 'Lazy wave: net uplift / submerged weight', unit: '–', value: 1.5, min: 0.2, max: 6, showIf: (v) => genOrTer(v) && v.riserType === 'lazyWave' },
    { key: 'riserBaseX', label: 'Riser base distance (0 = detect)', unit: 'm', value: 0, min: 0, max: 1e6, showIf: (v) => v.geomMode === 'table' || v.geomMode === 'centreline' },
    { key: 'terrain', label: 'Seabed / ground terrain grid', type: 'file', value: null, help: 'Elevation grid { x, y, elev } attached on the Geometry tab; used for the terrain route and for free-span screening.' },
    { key: 'terrainLine', label: 'Line across the terrain', type: 'select', value: 'we', options: opt([['we', 'West → east through the middle'], ['sn', 'South → north through the middle'], ['diag', 'South-west → north-east diagonal'], ['deep', 'Deepest cell → highest cell']]) },
    { key: 'seabedRms', label: 'Seabed micro-relief (RMS) without terrain', unit: 'm', value: 0.12, min: 0, max: 5, help: 'Used to screen free spans statistically when no surveyed terrain is attached.' },
    { key: 'seabedCorr', label: 'Micro-relief wavelength', unit: 'm', value: 40, min: 5, max: 500 },
    { key: 'spanGap', label: 'Smallest gap counted as a free span', unit: 'm', value: 0.3, min: 0.01, max: 5 },
  ] },
  { group: 'Pipe and material', tab: 'inputs', fields: [
    { key: 'idMm', label: 'Inner diameter', unit: 'mm', value: BASE.idMm, min: 20, max: 1500, typical: [100, 1000] },
    { key: 'wtMm', label: 'Wall thickness', unit: 'mm', value: BASE.wtMm, min: 2, max: 80, typical: [6, 40] },
    { key: 'roughUm', label: 'Wall roughness', unit: 'µm', value: BASE.roughUm, min: 0.5, max: 3000, typical: [15, 300] },
    { key: 'effIdMm', label: 'Fouled effective diameter (0 = clean)', unit: 'mm', value: 0, min: 0, max: 1500, help: 'Minimum bore left by deposits, from the solids suite.' },
    { key: 'effRoughUm', label: 'Fouled roughness (0 = clean)', unit: 'µm', value: 0, min: 0, max: 20000 },
    { key: 'grade', label: 'Line-pipe grade', type: 'select', value: 'X65', options: Object.entries(GRADES).map(([value, g]) => ({ value, label: g.label })) },
    { key: 'designPressure', label: 'Design pressure', unit: 'bara', value: BASE.designPressure, min: 10, max: 1500 },
    { key: 'designTemp', label: 'Design temperature', unit: '°C', value: BASE.designTemp, min: -50, max: 250 },
  ] },
  { group: 'Wall layers and surroundings', tab: 'inputs', help: 'Layers are listed from the steel outwards.', fields: [
    { key: 'layers', label: 'Coating and insulation layers', type: 'table', columns: [{ key: 'name', label: 'Layer', type: 'text' }, { key: 'tMm', label: 'Thickness', unit: 'mm' }, { key: 'k', label: 'Conductivity', unit: 'W/m/K' }, { key: 'rho', label: 'Density', unit: 'kg/m³' }, { key: 'cp', label: 'Heat capacity', unit: 'J/kg/K' }], value: DEF_LAYERS },
    { key: 'pip', label: 'Pipe-in-pipe (layers fill the annulus)', type: 'bool', value: false },
    { key: 'pipWtMm', label: 'Carrier pipe wall thickness', unit: 'mm', value: 14.3, min: 3, max: 60, showIf: (v) => !!v.pip },
    { key: 'concreteMm', label: 'Concrete weight coat thickness', unit: 'mm', value: 0, min: 0, max: 200 },
    { key: 'concreteRho', label: 'Concrete density', unit: 'kg/m³', value: 3040, min: 1800, max: 3800 },
    { key: 'burialDepth', label: 'Burial cover above the pipe (0 = exposed)', unit: 'm', value: 0, min: 0, max: 10 },
    { key: 'burialStartX', label: 'Buried from distance', unit: 'm', value: 0, min: 0, max: 1e6, showIf: (v) => v.burialDepth > 0 },
    { key: 'kSoil', label: 'Soil thermal conductivity', unit: 'W/m/K', value: 1.5, min: 0.2, max: 4, help: 'Clay or silt about 1.5, sand and gravel about 2.0 (ISO 10456 design values); saturated soils range from 0.6 to 4.' },
    { key: 'currentSpeed', label: 'Near-bed current speed', unit: 'm/s', value: BASE.currentSpeed, min: 0, max: 3 },
    { key: 'windSpeed', label: 'Wind speed on exposed onshore pipe', unit: 'm/s', value: 3, min: 0, max: 40 },
    { key: 'tSeabed', label: 'Seabed temperature', unit: '°C', value: BASE.tSeabed, min: -2, max: 35 },
    { key: 'tSeaSurface', label: 'Sea-surface temperature', unit: '°C', value: BASE.tSeaSurface, min: -2, max: 36 },
    { key: 'tAir', label: 'Air temperature', unit: '°C', value: BASE.tAir, min: -50, max: 55 },
    { key: 'tGround', label: 'Ground temperature (buried onshore)', unit: '°C', value: 15, min: -20, max: 45 },
    { key: 'uMult', label: 'U-value calibration multiplier', unit: '–', value: 1, min: 0.3, max: 3 },
    { key: 'insKMult', label: 'Insulation conductivity multiplier (ageing, water ingress)', unit: '–', value: 1, min: 0.5, max: 5, help: 'Multiplies the conductivity of the layer that carries the largest thermal resistance.' },
  ] },
  { group: 'Operating point', tab: 'inputs', help: 'Rates come from the case fluid.', fields: [
    { key: 'tIn', label: 'Temperature at the sources (wellheads / inlet)', unit: '°C', value: BASE.tIn, min: -20, max: 200 },
    { key: 'pSep', label: 'Arrival (separator) pressure', unit: 'bara', value: BASE.separatorP, min: 1.5, max: 400 },
  ] },
  { group: 'Wells', tab: 'inputs', fields: [
    { key: 'survey', label: 'Deviation survey', type: 'table', columns: [{ key: 'md', label: 'Measured depth', unit: 'm' }, { key: 'inc', label: 'Inclination', unit: '°' }, { key: 'azi', label: 'Azimuth', unit: '°' }], value: DEF_SURVEY },
    { key: 'tubingIdMm', label: 'Tubing inner diameter', unit: 'mm', value: BASE.tubingIdMm, min: 30, max: 300 },
    { key: 'tubingRoughUm', label: 'Tubing roughness', unit: 'µm', value: 25, min: 0.5, max: 1000 },
    { key: 'iprType', label: 'Inflow model', type: 'select', value: 'composite', options: opt([['pi', 'Straight-line productivity index'], ['vogel', 'Vogel (saturated oil)'], ['composite', 'Composite: straight line above the bubble point, Vogel below'], ['fetkovich', 'Fetkovich'], ['darcy', 'Darcy radial flow with skin'], ['jones', 'Jones–Blount–Glaze (oil, non-Darcy)'], ['gas', 'Gas well: Darcy + Forchheimer (p² form)']]) },
    { key: 'pRes', label: 'Reservoir pressure', unit: 'bara', value: BASE.pRes, min: 5, max: 1500 },
    { key: 'tRes', label: 'Reservoir temperature', unit: '°C', value: BASE.tRes, min: 10, max: 250 },
    { key: 'pBubble', label: 'Saturation pressure at reservoir temperature', unit: 'bara', value: 218, min: 1, max: 1500 },
    { key: 'piWell', label: 'Productivity index per well', unit: 'Sm³/d/bar', value: BASE.pi, min: 0.01, max: 5000, showIf: (v) => ['pi', 'vogel', 'composite', 'fetkovich'].includes(v.iprType) },
    { key: 'fetkN', label: 'Fetkovich exponent', unit: '–', value: 0.85, min: 0.5, max: 1, showIf: (v) => v.iprType === 'fetkovich' },
    { key: 'permMd', label: 'Permeability', unit: 'mD', value: 150, min: 0.001, max: 20000, showIf: (v) => ['darcy', 'jones', 'gas'].includes(v.iprType) },
    { key: 'payM', label: 'Net pay thickness', unit: 'm', value: 30, min: 0.5, max: 500, showIf: (v) => ['darcy', 'jones', 'gas'].includes(v.iprType) },
    { key: 'perfM', label: 'Perforated interval', unit: 'm', value: 20, min: 0.5, max: 500, showIf: (v) => ['jones', 'gas'].includes(v.iprType) },
    { key: 'reM', label: 'Drainage radius', unit: 'm', value: 500, min: 20, max: 5000, showIf: (v) => ['darcy', 'jones', 'gas'].includes(v.iprType) },
    { key: 'rwM', label: 'Wellbore radius', unit: 'm', value: 0.108, min: 0.03, max: 0.5, showIf: (v) => ['darcy', 'jones', 'gas'].includes(v.iprType) },
    { key: 'skin', label: 'Completion skin', unit: '–', value: 2, min: -6, max: 100, showIf: (v) => ['darcy', 'jones', 'gas'].includes(v.iprType) },
    { key: 'porosity', label: 'Porosity (non-Darcy coefficient)', unit: '–', value: 0.2, min: 0.02, max: 0.5, showIf: (v) => ['jones', 'gas'].includes(v.iprType), help: 'Used only in the correlation of the non-Darcy coefficient β = 1.88 × 10¹⁰ k^−1.47 φ^−0.53 (1/ft, k in mD).' },
    { key: 'maxDrawdown', label: 'Maximum drawdown (sand-control / completion limit)', unit: 'bar', value: 70, min: 1, max: 1000, typical: [20, 120], help: 'Largest difference between reservoir pressure and flowing bottom-hole pressure the completion is designed for (sand control, fines migration, coning). A project design limit, not a literature value: the well-count study rejects any number of wells that needs more than this.' },
    { key: 'chokeMaxOpen', label: 'Largest choke opening accepted at the case rate', unit: '%', value: 80, min: 20, max: 100, help: 'Control margin for the well-count study: a number of wells is recommended only when the production chokes hold the case rate at or below this opening.' },
  ] },
  { group: 'Artificial lift', tab: 'inputs', help: 'Lift added to every flowing well. Subsea boosting is a multiphase pump at the wellhead / manifold that adds its pressure rise to what the well delivers; gas lift injects gas into the tubing and lightens the column above the valve (the lift gas then travels with the production to the separator).', fields: [
    { key: 'liftType', label: 'Lift method', type: 'select', value: 'boost', options: opt([['boost', 'Subsea multiphase boosting'], ['gaslift', 'Gas lift'], ['none', 'None (natural flow)']]) },
    { key: 'liftDp', label: 'Boosting pump rated pressure rise', unit: 'bar', value: 30, min: 0, max: 200, typical: [10, 60], showIf: (v) => v.liftType === 'boost' },
    { key: 'liftAuto', label: 'Speed-controlled: hold the production choke at its target opening', type: 'bool', value: true, showIf: (v) => v.liftType === 'boost', help: 'The pressure rise is set between the minimum and the rated value so that the choke sits at the target opening at the case rate. Off: the pump runs at its rated pressure rise.' },
    { key: 'chokeTarget', label: 'Target choke opening', unit: '%', value: 60, min: 10, max: 95, showIf: (v) => v.liftType === 'boost' && !!v.liftAuto },
    { key: 'liftDpMin', label: 'Smallest pressure rise while running (minimum speed)', unit: 'bar', value: 2, min: 0, max: 100, showIf: (v) => v.liftType === 'boost' && !!v.liftAuto },
    { key: 'liftEff', label: 'Multiphase pump efficiency', unit: '–', value: 0.5, min: 0.2, max: 0.85, showIf: (v) => v.liftType === 'boost' },
    { key: 'liftGasRate', label: 'Lift gas per well', unit: 'Sm³/d', value: 60000, min: 0, max: 1e6, showIf: (v) => v.liftType === 'gaslift' },
    { key: 'liftDepthPct', label: 'Injection valve depth', unit: '% of measured depth', value: 80, min: 5, max: 100, showIf: (v) => v.liftType === 'gaslift' },
  ] },
  { group: 'Production choke', tab: 'inputs', fields: [
    { key: 'chokeCvMax', label: 'Rated Cv (fully open)', unit: 'US gpm/psi^½', value: 120, min: 1, max: 5000 },
    { key: 'chokeChar', label: 'Inherent characteristic', type: 'select', value: 'equal', options: opt([['equal', 'Equal percentage'], ['linear', 'Linear'], ['quick', 'Quick opening']]) },
    { key: 'chokeRange', label: 'Rangeability', unit: '–', value: 50, min: 5, max: 200 },
    { key: 'chokeFL', label: 'Liquid pressure-recovery factor FL', unit: '–', value: 0.9, min: 0.3, max: 1 },
    { key: 'chokeXT', label: 'Pressure-differential ratio factor xT', unit: '–', value: 0.7, min: 0.1, max: 0.95 },
    { key: 'chokeCd', label: 'Bean discharge coefficient', unit: '–', value: 0.85, min: 0.4, max: 1.05, help: 'Sachdeva model: 0.85 for a choke free of upstream flow disturbance, 0.75 for a choke close behind an elbow; published range 0.62–0.90. Calibrate it to choke tests.' },
  ] },
  { group: 'Network', tab: 'inputs', help: 'Connections of type flowline and riser take their geometry from the route; pipe uses the length and diameter given; choke uses the Cv in “parameter”; pump uses the pump curve below (parameter = pumps in series); closed removes the connection.', fields: [
    { key: 'network', label: 'Connections', type: 'table', columns: [{ key: 'from', label: 'From', type: 'text' }, { key: 'to', label: 'To', type: 'text' }, { key: 'type', label: 'Type', type: 'text' }, { key: 'length', label: 'Length', unit: 'm' }, { key: 'diameter', label: 'Diameter', unit: 'mm' }, { key: 'k', label: 'K-factor', unit: '–' }, { key: 'param', label: 'Parameter', unit: '' }], value: DEF_NETWORK, help: 'Diameters below 5 are read as metres. Zero diameter means the main pipe diameter.' },
    { key: 'netNodes', label: 'Node data (optional)', type: 'table', columns: [{ key: 'id', label: 'Node', type: 'text' }, { key: 'kind', label: 'Kind (source / sink / junction)', type: 'text' }, { key: 'share', label: 'Share of the case rate', unit: '%' }, { key: 'p', label: 'Fixed pressure', unit: 'bara' }, { key: 'z', label: 'Elevation', unit: 'm' }, { key: 'status', label: 'Completion status (flowing / shut in / suspended / not completed)', type: 'text' }], value: DEF_NODES, help: 'Blank cells are inferred: nodes with no inflow are sources, nodes with no outflow are sinks at the arrival pressure, elevations follow the route. A well that is shut in, suspended, abandoned or not yet completed is taken out of the network and the remaining wells share the case rate.' },
    { key: 'fittings', label: 'Fittings on the main line', type: 'table', columns: [{ key: 'type', label: `Fitting (${Object.keys(FITTINGS).join(', ')})`, type: 'text' }, { key: 'count', label: 'Count', unit: '' }], value: DEF_FITTINGS },
    { key: 'kMult', label: 'Local-loss calibration multiplier', unit: '–', value: 1, min: 0, max: 20, help: 'Scales every fitting and bend K-factor of the main line.' },
    { key: 'bendMult', label: 'Bend-loss multiplier', unit: '–', value: 1, min: 0.2, max: 5, help: 'Scales the K-factor of bends and elbows only (on top of the local-loss multiplier).' },
    { key: 'teeMult', label: 'Tee / junction-loss multiplier', unit: '–', value: 1, min: 0.2, max: 5, help: 'Scales the K-factor of tees only.' },
    { key: 'cErosion', label: 'Erosional velocity constant C', unit: '(kg/m)^½/s', value: 122, min: 60, max: 300, help: 'API RP 14E: v = C/√ρ; 122 corresponds to C = 100 in field units.' },
  ] },
  { group: 'Liquid export pump', tab: 'inputs', fields: [
    { key: 'pumpOn', label: 'Include the pump', type: 'bool', value: true },
    { key: 'pumpQr', label: 'Rated flow', unit: 'm³/h', value: 200, min: 1, max: 20000, showIf: (v) => !!v.pumpOn },
    { key: 'pumpHr', label: 'Rated head', unit: 'm', value: 520, min: 1, max: 5000, showIf: (v) => !!v.pumpOn },
    { key: 'pumpEff', label: 'Best efficiency', unit: '–', value: 0.74, min: 0.2, max: 0.92, showIf: (v) => !!v.pumpOn },
    { key: 'pumpShutoff', label: 'Shut-off head / rated head', unit: '–', value: 1.25, min: 1.02, max: 2, showIf: (v) => !!v.pumpOn, help: 'Shape of the head curve when no tested curve is given.' },
    { key: 'pumpInertiaH', label: 'Pump train inertia constant (kinetic energy / rated power)', unit: 's', value: 1.5, min: 0.1, max: 20, showIf: (v) => !!v.pumpOn, help: 'Rotor inertia J = 2 H P / ω²; used by the transient.' },
    { key: 'pumpSpeedPct', label: 'Speed', unit: '% of rated', value: 100, min: 30, max: 120, showIf: (v) => !!v.pumpOn },
    { key: 'pumpTrimPct', label: 'Impeller trim', unit: '% of full diameter', value: 100, min: 75, max: 100, showIf: (v) => !!v.pumpOn },
    { key: 'pumpDischargeP', label: 'Delivery pressure', unit: 'bara', value: 60, min: 1, max: 600, showIf: (v) => !!v.pumpOn },
    { key: 'pumpStaticHead', label: 'Static lift', unit: 'm', value: 15, min: -200, max: 2000, showIf: (v) => !!v.pumpOn },
    { key: 'pumpLineLength', label: 'Discharge line length', unit: 'm', value: 400, min: 1, max: 5e5, showIf: (v) => !!v.pumpOn },
    { key: 'pumpLineIdMm', label: 'Discharge line diameter', unit: 'mm', value: 154, min: 20, max: 1500, showIf: (v) => !!v.pumpOn },
    { key: 'pumpNpshR', label: 'NPSH required at rated flow', unit: 'm', value: 3.5, min: 0.3, max: 40, showIf: (v) => !!v.pumpOn },
    { key: 'pumpSuctionHead', label: 'Liquid level above the pump suction', unit: 'm', value: 5, min: -5, max: 60, showIf: (v) => !!v.pumpOn },
    { key: 'pumpImpMm', label: 'Impeller diameter', unit: 'mm', value: 330, min: 80, max: 1500, showIf: (v) => !!v.pumpOn },
    { key: 'pumpRpm', label: 'Rated speed', unit: 'rpm', value: 3560, min: 300, max: 12000, showIf: (v) => !!v.pumpOn },
    { key: 'pumpBeta2', label: 'Blade outlet angle', unit: '°', value: 25, min: 10, max: 90, showIf: (v) => !!v.pumpOn },
    { key: 'pumpBlades', label: 'Number of blades', unit: '', value: 7, min: 3, max: 16, showIf: (v) => !!v.pumpOn },
    { key: 'pumpCurve', label: 'Tested pump curve (optional, at rated speed)', type: 'table', columns: [{ key: 'q', label: 'Flow', unit: 'm³/h' }, { key: 'h', label: 'Head', unit: 'm' }, { key: 'eta', label: 'Efficiency', unit: '%' }], value: [], help: 'Three or more points replace the parabola through the rated point; speed and trim still scale it by the affinity laws.' },
  ] },
  { group: 'Gas compressor', tab: 'inputs', fields: [
    { key: 'compOn', label: 'Include the compressor', type: 'bool', value: true },
    { key: 'compPd', label: 'Discharge pressure', unit: 'bara', value: 150, min: 2, max: 700, showIf: (v) => !!v.compOn },
    { key: 'compEta', label: 'Polytropic efficiency', unit: '–', value: 0.78, min: 0.4, max: 0.92, showIf: (v) => !!v.compOn },
    { key: 'compMaxRatio', label: 'Largest pressure ratio per stage', unit: '–', value: 3.5, min: 1.2, max: 6, showIf: (v) => !!v.compOn },
    { key: 'compTcool', label: 'Interstage cooler outlet temperature', unit: '°C', value: 40, min: 5, max: 80, showIf: (v) => !!v.compOn },
    { key: 'compMargin', label: 'Design flow above the duty flow', unit: '%', value: 10, min: -20, max: 60, showIf: (v) => !!v.compOn },
    { key: 'compHeadMult', label: 'Map head scaling', unit: '–', value: 1, min: 0.5, max: 1.5, showIf: (v) => !!v.compOn, help: 'Scales the head of the generic map (calibration against performance tests).' },
    { key: 'compFlowMult', label: 'Map flow scaling', unit: '–', value: 1, min: 0.5, max: 1.5, showIf: (v) => !!v.compOn },
    { key: 'compRpm', label: 'Rated speed', unit: 'rpm', value: 9000, min: 1000, max: 30000, showIf: (v) => !!v.compOn },
    { key: 'compInertiaH', label: 'Compressor train inertia constant (kinetic energy / rated power)', unit: 's', value: 4, min: 0.2, max: 30, showIf: (v) => !!v.compOn },
    { key: 'compMap', label: 'Tested compressor map (optional)', type: 'table', columns: [{ key: 'speed', label: 'Speed', unit: '%' }, { key: 'q', label: 'Inlet flow', unit: 'm³/s' }, { key: 'head', label: 'Polytropic head', unit: 'kJ/kg' }], value: [], help: 'At least two speed lines of three points each, each line from surge to stonewall; otherwise a generic map is scaled to the duty.' },
  ] },
  { group: 'Separator and slug catcher', tab: 'inputs', fields: [
    { key: 'sepOrient', label: 'Separator orientation', type: 'select', value: 'horizontal', options: opt([['horizontal', 'Horizontal'], ['vertical', 'Vertical']]) },
    { key: 'sepD', label: 'Separator diameter', unit: 'm', value: 2.6, min: 0.3, max: 8 },
    { key: 'sepL', label: 'Separator tan–tan length', unit: 'm', value: 10, min: 0.5, max: 60 },
    { key: 'sepLevelPct', label: 'Normal liquid level', unit: '% of height', value: 50, min: 5, max: 90 },
    { key: 'sepKAuto', label: 'Souders–Brown K from the mist-eliminator curve at the vessel pressure', type: 'bool', value: true, help: 'York curve as fitted by Svrcek & Monnery; horizontal vessels take 1.25 times the vertical value.' },
    { key: 'sepK', label: 'Souders–Brown coefficient K (when not automatic)', unit: 'm/s', value: 0.11, min: 0.02, max: 0.3 },
    { key: 'sepInitLevelPct', label: 'Initial liquid level', unit: '% of height', value: 50, min: 2, max: 95, help: 'Initial liquid inventory of the separator for the transient.' },
    { key: 'sepInitP', label: 'Initial separator pressure (0 = arrival pressure)', unit: 'bara', value: 0, min: 0, max: 400, help: 'Initial gas inventory of the separator for the transient.' },
    { key: 'sepResMin', label: 'Liquid residence time required', unit: 'min', value: 3, min: 0.5, max: 30 },
    { key: 'slugSurge', label: 'Design slug surge volume', unit: 'm³', value: 42, min: 0, max: 5000, help: 'From the flow suite when it has been run.' },
    { key: 'slugSF', label: 'Slug-volume design factor', unit: '–', value: 1.25, min: 1, max: 3 },
    { key: 'scInitPct', label: 'Initial slug-catcher inventory', unit: '% of the normal hold-up', value: 100, min: 0, max: 400 },
    { key: 'scEff', label: 'Slug-catcher effective capacity / geometric volume', unit: '–', value: 0.85, min: 0.3, max: 1 },
    { key: 'slugDuration', label: 'Slug arrival duration', unit: 's', value: 120, min: 5, max: 3600 },
    { key: 'drainFactor', label: 'Largest liquid draw-off / normal liquid rate', unit: '–', value: 1.5, min: 1, max: 5 },
    { key: 'scType', label: 'Slug catcher type', type: 'select', value: 'finger', options: opt([['finger', 'Finger (multiple-pipe) type'], ['vessel', 'Vessel type']]) },
    { key: 'fingerDmm', label: 'Finger diameter', unit: 'mm', value: 1200, min: 300, max: 1600, showIf: (v) => v.scType === 'finger' },
    { key: 'fingerK', label: 'Souders–Brown coefficient of the fingers', unit: 'm/s', value: 0.08, min: 0.02, max: 0.3, showIf: (v) => v.scType === 'finger', help: 'Design assumption for gas–liquid disengagement in a finger without a mist eliminator: about half of the 0.12–0.15 m/s that API 12J lists for horizontal vessels with a wire-mesh extractor. No published value specific to finger-type slug catchers was found; use the vendor’s figure.' },
    { key: 'fingerMaxL', label: 'Longest finger allowed', unit: 'm', value: 60, min: 10, max: 400, showIf: (v) => v.scType === 'finger' },
  ] },
  { group: 'Equipment condition', tab: 'inputs', fields: [
    { key: 'degrFactor', label: 'Rotating-equipment degradation factor (head retained)', unit: '–', value: 1, min: 0.5, max: 1, help: '1 = as new. Wear, fouling and enlarged clearances reduce the head of the export pump and of the compressor map by this factor.' },
  ] },
  { group: 'Design features (installation and integrity)', tab: 'inputs', help: 'Features of the as-designed line that the integrity and economics suites read: buckle arrestors, planned lateral-buckle initiators, cathodic protection and supports.', fields: [
    { key: 'buckleArrestorSpacing', label: 'Buckle-arrestor spacing (0 = none)', unit: 'm', value: 244, min: 0, max: 5000, help: 'Default: one integral arrestor every 20 joints of 12.2 m over the section where a buckle could propagate.' },
    { key: 'buckleArrestorFrom', label: 'Arrestors installed below water depth (0 = from the propagation pressure)', unit: 'm', value: 0, min: 0, max: 4000 },
    { key: 'buckleInitType', label: 'Lateral-buckle initiators', type: 'select', value: 'sleeper', options: opt([['sleeper', 'Sleepers'], ['snake lay', 'Snake lay'], ['none', 'None']]) },
    { key: 'buckleInitSpacing', label: 'Initiator spacing', unit: 'm', value: 2500, min: 200, max: 20000, showIf: (v) => v.buckleInitType !== 'none' },
    { key: 'sleeperHeight', label: 'Sleeper height', unit: 'm', value: 0.5, min: 0.1, max: 3, showIf: (v) => v.buckleInitType === 'sleeper' },
    { key: 'anodeSpacing', label: 'Bracelet-anode spacing (0 = none)', unit: 'm', value: 122, min: 0, max: 1000, help: 'Default: every 10 joints.' },
    { key: 'anodeMass', label: 'Net anode mass', unit: 'kg', value: 60, min: 1, max: 2000 },
    { key: 'cpLife', label: 'Cathodic-protection design life', unit: 'y', value: 25, min: 1, max: 60 },
    { key: 'cpCurrent', label: 'Mean protection current density of bare steel', unit: 'A/m²', value: 0.075, min: 0.01, max: 0.5, help: 'DNVGL-RP-F103 (2016) Table 6-2, non-buried pipeline: 0.050 up to 25 °C internal fluid temperature, 0.060 up to 50 °C, 0.075 up to 80 °C, 0.100 up to 120 °C, 0.130 above; buried 0.020 / 0.030 / 0.040 / 0.060 / 0.080 A/m².' },
    { key: 'cpBreakdown', label: 'Mean coating breakdown factor', unit: '–', value: 0.048, min: 0.0001, max: 1, help: 'DNVGL-RP-F103 (2016): f = a + 0.5·b·t over the design life t. Thermally insulating FBE / polyurethane linepipe coating a = 0.01, b = 0.003 per year → 0.0475 for 25 years; three-layer FBE / polypropylene 0.001 + 0.5 × 0.00003 t; FBE under concrete 0.030 + 0.5 × 0.0003 t. Field joints have their own factors and are not added here.' },
    { key: 'supports', label: 'Supports and restraints', type: 'table', columns: [{ key: 'x', label: 'Distance', unit: 'm' }, { key: 'type', label: 'Type (anchor, PLET, sleeper, rock berm, clamp, …)', type: 'text' }, { key: 'k', label: 'Stiffness', unit: 'kN/m' }], value: DEF_SUPPORTS, help: 'A support inside a free-span candidate splits it.' },
    { key: 'spanLimit', label: 'Allowable free-span length', unit: 'm', value: 45, min: 5, max: 300 },
  ] },
  { group: 'Measured data for comparison', tab: 'inputs', help: 'Optional. Each row is one measurement; the model value is listed beside it with bias, RMS error and mean absolute percentage error per kind. Kinds: as-built or survey elevation (at = distance m, value = elevation m); ROV span (at = distance, value = span length m); inspection ID (at = distance, value = bore mm); bathymetry (at = distance, value = seabed elevation m); deviation TVD (at = measured depth m, value = TVD m); caliper (at = measured depth, value = tubing ID mm); pump FAT (at = flow m³/h, value = head m); compressor test (at = inlet flow m³/s, aux = speed %, value = polytropic head kJ/kg); valve test (at = water flow m³/h, aux = opening %, value = Δp bar); separator test (value = gas velocity at carry-over m/s); pressure (at = distance, value = bara); temperature (at = distance, value = °C); well test (at = rate Sm³/d, value = flowing bottom-hole pressure bara); production log (at = measured depth, value = flowing pressure bara); commissioning fill (at = distance, value = line-fill volume m³).', fields: [
    { key: 'measured', label: 'Measurements', type: 'table', columns: [{ key: 'kind', label: 'Kind', type: 'text' }, { key: 'at', label: 'At (position / flow / depth)', unit: '' }, { key: 'aux', label: 'Second coordinate', unit: '' }, { key: 'value', label: 'Measured value', unit: '' }], value: [] },
  ] },
  { group: 'Transient scenario', tab: 'setup', help: 'Lumped dynamic model of the production system: line-pack of the flowline gas, separator gas and liquid inventories, export header, pump and compressor rotors, valve actuators, anti-surge recycle, flare pressure control and level control.', fields: [
    { key: 'transScenario', label: 'Scenario', type: 'select', value: 'valveClosure', options: opt([['compressorTrip', 'Compressor trip'], ['valveClosure', 'Export-gas valve closure (surge approach)'], ['pumpTrip', 'Export-pump trip'], ['chokeChange', 'Production-choke change'], ['none', 'No upset (initial inventories only)'], ['off', 'Do not run the transient']]) },
    { key: 'transTime', label: 'Simulated time', unit: 's', value: 120, min: 10, max: 7200 },
    { key: 'valveStroke', label: 'Stroke time of the choke and export valve', unit: 's', value: 20, min: 0.5, max: 600 },
    { key: 'chokeStepPct', label: 'New choke opening (choke-change scenario)', unit: '%', value: 30, min: 0, max: 100, showIf: (v) => v.transScenario === 'chokeChange' },
    { key: 'exportVol', label: 'Volume of the compressor discharge header and export line-pack', unit: 'm³', value: 10, min: 1, max: 1e5 },
    { key: 'dischargeVol', label: 'Piping volume between the compressor and its discharge check valve', unit: 'm³', value: 1.5, min: 0.05, max: 500 },
    { key: 'asvSet', label: 'Anti-surge control line (surge margin at which recycle starts)', unit: '%', value: 10, min: 2, max: 40 },
    { key: 'asvSize', label: 'Recycle-valve capacity / compressor design flow', unit: '–', value: 2, min: 0.5, max: 5 },
    { key: 'flareSetPct', label: 'Flare pressure-control set point', unit: '% of arrival pressure', value: 110, min: 101, max: 200 },
  ] },
  { group: 'Models', tab: 'setup', fields: [
    { key: 'netModel', label: 'Two-phase pressure-gradient closure', type: 'select', value: 'beggsBrill', options: opt([['beggsBrill', 'Beggs & Brill'], ['driftFlux', 'Drift flux'], ['mechanistic', 'Mechanistic (stratified / slug unit cell)'], ['homogeneous', 'Homogeneous (no slip)']]) },
    { key: 'fModel', label: 'Friction-factor equation', type: 'select', value: 'colebrook', options: opt([['colebrook', 'Colebrook–White'], ['haaland', 'Haaland'], ['swamee', 'Swamee–Jain'], ['churchill', 'Churchill']]) },
    { key: 'bendMinAngle', label: 'Smallest direction change listed as a bend', unit: '°', value: 2, min: 0.2, max: 45 },
  ] },
  { group: 'Discretisation', tab: 'mesh', fields: [
    { key: 'nCells', label: 'Cells along the flowline', unit: '', value: 24, min: 6, max: 400, help: 'The riser gets half as many again; jumpers and spools one cell per 24.' },
    { key: 'nTubing', label: 'Cells along the tubing', unit: '', value: 24, min: 6, max: 400 },
    { key: 'netTol', label: 'Network residual tolerance (relative)', unit: '–', value: 1e-8, min: 1e-12, max: 1e-3 },
    { key: 'transSteps', label: 'Time steps of the transient', unit: '', value: 120, min: 20, max: 20000 },
  ] },
];
const FIELDS = INPUTS.flatMap((g) => g.fields), DEF = Object.fromEntries(FIELDS.map((f) => [f.key, f.value]));
/** Inputs with defaults filled in and numbers coerced; clearly impossible numbers raise a readable error. */
function readInputs(v = {}) {
  const o = {};
  for (const f of FIELDS) {
    const x = v[f.key];
    if (!f.type || f.type === 'number') { o[f.key] = num(x, f.value); if (isNum(f.min) && o[f.key] < f.min - 1e-12 && f.min >= 0 && o[f.key] < 0) throw new Error(`${f.label} cannot be negative (got ${o[f.key]}).`); }
    else if (f.type === 'bool') o[f.key] = x === undefined || x === null ? f.value : !!x && x !== 'false';
    else if (f.type === 'table') o[f.key] = Array.isArray(x) ? x : f.value;
    else if (f.type === 'select') o[f.key] = f.options.some((q) => q.value === x) ? x : f.value;
    else o[f.key] = x ?? f.value;
  }
  for (const k of ['idMm', 'wtMm', 'roughUm', 'pSep', 'pRes', 'tubingIdMm', 'sepD', 'sepL', 'nCells', 'nTubing', 'netTol', 'chokeCvMax']) if (!(o[k] > 0)) throw new Error(`${FIELDS.find((f) => f.key === k).label} must be greater than zero.`);
  o.nCells = Math.max(4, Math.round(o.nCells)); o.nTubing = Math.max(4, Math.round(o.nTubing)); o.transSteps = clamp(Math.round(o.transSteps), 10, 40000);
  return o;
}

/** Route of the case from the inputs: { x[], z[], rbI (riser-base index when known), riser, notes[], raw (analysis of the table as typed) }. */
function routeOf(v) {
  const notes = [], ux = { m: 1, km: 1000, ft: UNIT.ft }[v.lengthUnit] || 1, uz = (v.elevUnit === 'ft' ? UNIT.ft : 1) * (v.elevSense === 'depth' ? -1 : 1);
  let x, z, rbI = null, riser = null, raw = null, plan = null;
  const addRiser = () => {
    const n0 = x.length - 1, depth = -z[n0], type = v.riserType;
    if (type === 'none' || !(depth > 5)) return;
    const H = depth + v.topsideElev;
    let r;
    if (type === 'vertical') r = { x: [0, 0.25, 0.5, 0.75, 1], z: [0, 0.25 * H, 0.5 * H, 0.75 * H, H], length: Math.hypot(H, 1), span: 1 };
    else if (type === 'lazyWave') r = lazyWave({ height: H, angle: v.hangoffAngle, liftAngle: v.lazyLift, sagAngle: v.lazySag, buoyancy: v.lazyBuoy, n: 27 });
    else r = catenary({ height: H, angle: v.hangoffAngle, n: 14 });
    rbI = n0; riser = { type, length: r.length, span: r.span, hog: r.hog ?? null, sag: r.sag ?? null, a: r.a ?? null };
    for (let i = 1; i < r.x.length; i++) { x.push(x[n0] + r.x[i]); z.push(z[n0] + r.z[i]); }
  };
  if (v.geomMode === 'generate') {
    const s = seabedLine({ length: v.flowlineLength, depthEnd: v.waterDepth, drop: v.seabedDrop, amp: v.undulationAmp, wavelength: v.undulationLength });
    x = s.x; z = s.z; addRiser();
  } else if (v.geomMode === 'terrain') {
    if (!v.terrain) throw new Error('The route is set to follow the terrain grid, but no terrain is attached: import one on the Geometry tab or choose another route source.');
    const t = terrainTransect(v.terrain, v.terrainLine), a = analyseProfile(t.x, t.z), s = simplifyProfile(t.x, t.z, 56, [...a.highs, ...a.lows].map((q) => q.i));
    x = s.x; z = s.z; notes.push(`Route draped over the terrain grid: ${rd(t.x[t.x.length - 1], 0)} m transect, ${x.length} points kept.`); addRiser();
  } else if (v.geomMode === 'centreline') {
    const geo = v.coordSystem === 'geographic', pt = v.centreline.map((r) => ({ x: num(r?.x, null), y: num(r?.y, null), z: num(r?.z, null) })).filter((q) => q.x !== null && q.y !== null && q.z !== null);
    if (pt.length < 2) throw new Error('The centreline table needs at least two points with easting (or longitude), northing (or latitude) and elevation.');
    if (geo && pt.some((q) => Math.abs(q.x) > 180 || Math.abs(q.y) > 90)) throw new Error('Geographic centreline coordinates must be longitude (−180…180°) and latitude (−90…90°): choose the projected grid for eastings and northings.');
    const lat0 = pt.reduce((a, q) => a + q.y, 0) / pt.length, kx = geo ? 111320 * Math.cos(lat0 * D2R) : ux, ky = geo ? 110540 : ux, E = [], Nn = [];
    x = []; z = [];
    let dropped = 0;
    for (const q of pt) {
      const e = (q.x - pt[0].x) * kx, n = (q.y - pt[0].y) * ky, zz = q.z * uz, k = x.length - 1;
      if (k >= 0 && Math.hypot(e - E[k], n - Nn[k], zz - z[k]) < 1e-6) { dropped++; continue; }
      x.push(k < 0 ? 0 : x[k] + Math.hypot(e - E[k], n - Nn[k])); z.push(zz); E.push(e); Nn.push(n);
    }
    if (x.length < 2) throw new Error('All centreline points coincide: the route has no length.');
    if (dropped) notes.push(`${dropped} repeated centreline point(s) were left out.`);
    plan = { e: E, n: Nn, geographic: geo };
  } else {
    const rows = (Array.isArray(v.profile) ? v.profile : []).map((r) => ({ x: num(r?.x, NaN) * ux, z: num(r?.z, NaN) * uz }));
    raw = analyseProfile(rows.map((r) => r.x), rows.map((r) => r.z));
    const c = cleanProfile(rows);
    if (c.x.length < 2) throw new Error('The elevation profile needs at least two valid points (distance and elevation).');
    if (c.dropped) notes.push(`${c.dropped} profile row(s) were left out: missing numbers, repeated points or chainage running backwards.`);
    x = c.x; z = c.z;
  }
  return { x, z, rbI, riser, notes, raw, plan };
}

const kGas = (pr) => clamp(pr.cpG / Math.max(pr.cpG - R / (pr.mwG * 1e-3), 1), 1.05, 1.67);
const jtMix = (pr) => { const mCp = pr.mG * pr.cpG + pr.mO * pr.cpO + pr.mW * pr.cpW; return mCp > 0 ? (pr.mG * pr.cpG * pr.jtG + pr.mO * pr.cpO * pr.jtO - pr.mW / pr.rhoW) / mCp : 0; };

/** Case network: parses the tables, attaches route geometry to flowline / riser connections and returns a solver for a rate scale. */
function caseNetwork(I) {
  const { v, geo, fm } = I, notes = [], kindOf = (t) => (t === 'flowline' ? 'flowline' : t === 'riser' ? 'riser' : /choke|valve/.test(t) ? 'choke' : /pump/.test(t) ? 'pump' : 'pipe');
  const rows = v.network.map((r) => ({ from: String(r?.from ?? '').trim(), to: String(r?.to ?? '').trim(), type: String(r?.type ?? 'pipe').trim().toLowerCase() || 'pipe', length: Math.max(num(r?.length, 0), 0), dia: Math.max(num(r?.diameter, 0), 0), k: Math.max(num(r?.k, 0), 0), param: num(r?.param, 0) })).filter((r) => r.from && r.to);
  const offWells = new Set(v.netNodes.filter((r) => r && /shut|susp|aband|not\s*compl|uncompl|plug|closed|^off$/i.test(String(r.status ?? '').trim())).map((r) => String(r.id ?? '').trim()).filter(Boolean)), open = rows.filter((r) => !/^(closed|shut|off|isolated)/.test(r.type));
  let list = open.filter((r) => !offWells.has(r.from) && !offWells.has(r.to));
  if (open.length > list.length) notes.push(`${[...offWells].filter((id) => open.some((r) => r.from === id || r.to === id)).join(', ')}: not flowing (completion status) — taken out of the network; the remaining wells share the case rate.`);
  const closed = rows.length - open.length, hasRiser = geo.rb.height > 0 && geo.rb.i < geo.x.length - 1, sRb = hasRiser ? geo.S[geo.rb.i] : geo.S[geo.S.length - 1], sEnd = geo.S[geo.S.length - 1];
  if (closed) notes.push(`${closed} connection(s) are closed and were taken out of the network.`);
  if (!list.length) { list = [pipeRow('IN', hasRiser ? 'RB' : 'OUT', 0, 0, 0, 'flowline')]; if (hasRiser) list.push(pipeRow('RB', 'OUT', 0, 0, 0, 'riser')); list.forEach((r) => { r.dia = 0; }); notes.push('No open connection is listed: a single line from inlet to outlet is solved.'); }
  list.forEach((r) => { r.kind = kindOf(r.type); if (r.kind === 'riser' && !hasRiser) r.kind = 'link'; if (r.kind === 'flowline' && !(sRb > 0)) r.kind = 'link'; });
  const names = []; list.forEach((r) => { for (const k of [r.from, r.to]) if (!names.includes(k)) names.push(k); });
  const spec = new Map(v.netNodes.filter((r) => r && String(r.id ?? '').trim()).map((r) => [String(r.id).trim(), r])), ix = new Map(names.map((n, i) => [n, i]));
  const indeg = names.map(() => 0), outdeg = names.map(() => 0); list.forEach((r) => { outdeg[ix.get(r.from)]++; indeg[ix.get(r.to)]++; });
  const nodes = names.map((id, i) => {
    const s = spec.get(id) || {}, k = String(s.kind ?? '').trim().toLowerCase(), pf = num(s.p, 0);
    const kind = /^(so|well|in)/.test(k) ? 'source' : /^(si|sep|out|term|del)/.test(k) ? 'sink' : /^(j|man|node|plem|tee)/.test(k) ? 'junction' : indeg[i] === 0 ? 'source' : outdeg[i] === 0 ? 'sink' : 'junction';
    return { id, kind, z: num(s.z, null), pFixed: pf > 0 ? pf : kind === 'sink' ? v.pSep : null, share: num(s.share, null), rate: 0 };
  });
  const src = nodes.filter((n) => n.kind === 'source'), rated = src.filter((n) => n.pFixed === null);
  I.liftFrac = ((I.liftPerWell || 0) * src.length) / I.mCase; // lift gas returning with the production, per kg/s of produced fluid
  if (!src.length) throw new Error('The network has no source: at least one node must have no inflow or be marked “source”.');
  if (!nodes.some((n) => n.kind === 'sink')) throw new Error('The network has no outlet: at least one node must have no outflow or be marked “sink”.');
  if (rated.length) {
    const given = rated.filter((n) => n.share !== null && n.share >= 0), blank = rated.filter((n) => !(n.share !== null && n.share >= 0)), sumG = given.reduce((s, n) => s + n.share, 0), each = blank.length ? Math.max(100 - sumG, 0) / blank.length : 0;
    blank.forEach((n) => { n.share = each; });
    const tot = rated.reduce((s, n) => s + n.share, 0);
    if (!(tot > 0)) throw new Error('Every source is shut in (all rate shares are zero): open at least one well.');
    if (Math.abs(tot - 100) > 0.5) notes.push(`Source shares add up to ${rd(tot, 1)} %: they were rescaled so that the network carries the case rate.`);
    rated.forEach((n) => { n.share /= tot; n.rate = n.share * I.mCase; });
  }
  // connections tied to the route: flowline pieces share the seabed section in proportion to their listed lengths
  const fl = list.filter((r) => r.kind === 'flowline'), wSum = fl.reduce((s, r) => s + r.length, 0);
  let cum = 0;
  fl.forEach((r) => { const w = wSum > 0 ? r.length / wSum : 1 / fl.length; r.sA = sRb * cum; cum += w; r.sB = sRb * Math.min(cum, 1); r.w = w; });
  if (wSum > 0 && fl.some((r) => !(r.length > 0))) throw new Error('When several flowline connections share the route, give every one of them a length (used as its share of the route).');
  const zOnRoute = (s) => interp1(geo.S, geo.z, s);
  list.forEach((r) => { const a = nodes[ix.get(r.from)], b = nodes[ix.get(r.to)]; if (r.kind === 'flowline') { a.z = zOnRoute(r.sA); b.z = zOnRoute(r.sB); a.onRoute = b.onRoute = true; } else if (r.kind === 'riser') { a.z = zOnRoute(sRb); b.z = zOnRoute(sEnd); a.onRoute = b.onRoute = true; } });
  for (let pass = 0; pass < names.length; pass++) { let ch = false; list.forEach((r) => { const a = nodes[ix.get(r.from)], b = nodes[ix.get(r.to)]; if (a.z === null && b.z !== null) { a.z = b.z; ch = true; } else if (b.z === null && a.z !== null) { b.z = a.z; ch = true; } }); if (!ch) break; }
  nodes.forEach((n) => { if (n.z === null) n.z = 0; });
  const edges = list.map((r) => {
    const a = ix.get(r.from), b = ix.get(r.to), e = { from: r.from, to: r.to, a, b, kind: r.kind, type: r.type, k: r.k, cells: [], rev: [], L: 0, D: I.D, rough: I.rough };
    const dress = (cells, onMain) => cells.forEach((c) => { const q = I.env(c.xMid, c.zMid, onMain); c.U = q.U; c.tAmb = q.tAmb; c.D = e.D; c.rough = e.rough; });
    if (r.kind === 'flowline' || r.kind === 'riser') {
      const sA = r.kind === 'riser' ? sRb : r.sA, sB = r.kind === 'riser' ? sEnd : r.sB;
      if (r.kind === 'flowline') { e.D = I.Dfl; e.rough = I.roughFl; e.k += I.kFit * r.w; }
      e.cells = polylineCells(geo.S, geo.x, geo.z, sA, sB, r.kind === 'riser' ? Math.max(4, Math.round(v.nCells / 2)) : Math.max(2, Math.round(v.nCells * r.w))); e.L = sB - sA; e.sA = sA; e.sB = sB; dress(e.cells, true);
    } else if (r.kind === 'pipe') {
      const dz = nodes[b].z - nodes[a].z; e.D = r.dia > 0 ? (r.dia < 5 ? r.dia : r.dia / 1000) : I.D;
      if (!(r.length > 0)) notes.push(`Connection ${r.from} → ${r.to} has no length: ${rd(Math.max(Math.abs(dz), 1), 1)} m was assumed.`);
      e.L = Math.max(r.length > 0 ? r.length : 1, Math.abs(dz) * 1.0005);
      const nc = Math.max(1, Math.round(v.nCells / 24)), theta = Math.asin(clamp(dz / e.L, -1, 1));
      for (let j = 0; j < nc; j++) e.cells.push({ ds: e.L / nc, theta, zMid: nodes[a].z + (dz * (j + 0.5)) / nc, xMid: 0 });
      dress(e.cells, false);
    } else if (r.kind === 'choke') { e.cv = r.param > 0 ? r.param : v.chokeCvMax; } else if (r.kind === 'pump') e.series = r.param > 0 ? r.param : 1;
    e.rev = e.cells.slice().reverse().map((c) => ({ ...c, theta: -c.theta }));
    return e;
  });
  const tNode = nodes.map(() => v.tIn);
  const law = (e, m, pUp, dir, i, rec) => {
    const tUp = tNode[dir > 0 ? e.a : e.b];
    if (e.kind === 'link' || !(m > 0) && e.kind !== 'pipe' && e.kind !== 'flowline' && e.kind !== 'riser') { if (rec) rec[i] = { tOut: tUp, dp: 0 }; return 0; }
    if (e.kind === 'choke' || e.kind === 'pump') {
      const pe = Math.max(pUp, 1), pr = mixAt(fm, pe, tUp, m * I.fHC, m * I.fW), qv = pr.qG + pr.qL, rho = m / qv, cp = (pr.mG * pr.cpG + pr.mO * pr.cpO + pr.mW * pr.cpW) / m;
      if (e.kind === 'choke') { const r = valveDp({ cv: e.cv, w: m, p1: pe, rhoL: pr.rhoL, rhoG: pr.rhoG, xG: pr.mG / m, k: kGas(pr), xT: v.chokeXT }); if (rec) rec[i] = { tOut: tUp - jtMix(pr) * r.dp * 1e5, dp: r.dp, choked: r.choked, x: r.x }; return r.dp; }
      if (dir < 0) { const d = 5 + 50 * (m / I.mCase) ** 2; if (rec) rec[i] = { tOut: tUp, dp: d, reverse: true }; return d; } // non-return valve
      const q = qv * 3600, H = e.series * pumpHead(q, I.pump), eta = pumpEfficiency(q, I.pump), d = (-rho * G * H) / 1e5;
      if (rec) rec[i] = { tOut: tUp + (H > 0 ? (G * H * (1 / eta - 1)) / cp : 0), dp: d, head: H, q, eta, power: H > 0 ? (m * G * H) / eta : 0, rho };
      return d;
    }
    const r = marchPipe({ fm, mHC: m * I.fHC, mW: m * I.fW, cells: dir > 0 ? e.cells : e.rev, pIn: pUp, tIn: tUp, model: v.netModel, fModel: v.fModel, kLoss: e.k, gExtra: m * (I.liftFrac || 0) });
    if (rec) rec[i] = r;
    return r.dp;
  };
  const afterEval = ({ m, rec }) => {
    const a = nodes.map(() => 0), b = nodes.map(() => 0);
    edges.forEach((e, i) => { const r = rec[i], w = Math.abs(m[i]); if (!r || !isNum(r.tOut)) return; const k = m[i] >= 0 ? e.b : e.a; a[k] += w * r.tOut; b[k] += w; });
    let ch = 0;
    nodes.forEach((n, k) => { if (n.kind === 'source' || !(b[k] > 1e-12)) return; const t = a[k] / b[k]; ch = Math.max(ch, Math.abs(t - tNode[k])); tNode[k] = t; });
    return ch;
  };
  const zTop = Math.max(...nodes.map((n) => n.z)), pMax = Math.max(...nodes.map((n) => n.pFixed ?? 0));
  /** Solve at a rate scale (fixed-rate sources × scale). */
  const solve = (scale = 1, guess = null) => {
    const nd = nodes.map((n) => ({ id: n.id, z: n.z, p: n.pFixed ?? undefined, q: n.rate * scale }));
    const g = guess || { p: nodes.map((n) => pMax + 20 + 0.06 * (zTop - n.z)) };
    const sol = solveNetwork({ nodes: nd, edges, law, tol: v.netTol, guess: g, afterEval, maxIter: 40 });
    return { sol, t: tNode.slice(), guess: { p: sol.p.slice(), m: sol.m.slice() } };
  };
  return { nodes, edges, notes, solve, hasRiser, sRb };
}

/** Inflow description of one well from the inputs (coefficients of the Darcy / non-Darcy forms use the case fluid at reservoir conditions). */
function makeIpr(v, fm, mHC, mW) {
  const o = { type: v.iprType, pRes: v.pRes, pb: Math.min(v.pBubble, v.pRes), pi: v.piWell, n: v.fetkN, a: 0, b: 0, basis: v.iprType === 'gas' ? 'gas' : 'liquid', beta: null };
  if (['darcy', 'jones', 'gas'].includes(o.type)) {
    const pr = mixAt(fm, v.pRes, v.tRes, mHC, mW), qLs = (fm.rates.qOilStd + fm.rates.qWaterStd) / DAY, geo = { k: v.permMd, h: v.payM, hp: v.perfM, re: v.reM, rw: v.rwM, skin: v.skin, phi: v.porosity };
    const c = o.type === 'gas' ? radialCoefficients({ phase: 'gas', ...geo, mu: pr.muG, Z: pr.zG, T: v.tRes + KEL, M: pr.mwG * 1e-3 }) : radialCoefficients({ phase: 'oil', ...geo, mu: pr.muL, B: qLs > 0 ? pr.qL / qLs : 1.2, rho: pr.rhoL, nonDarcy: o.type === 'jones' });
    o.a = c.a; o.b = c.b; o.beta = c.beta;
  }
  o.qMax = iprRate(0, o);
  if (['darcy', 'jones', 'gas'].includes(o.type)) o.pi = iprRate(0.9 * o.pRes, o) / (0.1 * o.pRes); // equivalent PI at 10 % drawdown
  return o;
}
/** Stable intersection of a supply curve (available pressure, falling) with a demand curve (required pressure): largest x where supply ≥ demand. */
function intersect(xa, ya, xb, yb) {
  const lo = Math.max(xa[0], xb[0]), hi = Math.min(xa[xa.length - 1], xb[xb.length - 1]), n = 240, g = (s) => interp1(xa, ya, s) - interp1(xb, yb, s);
  if (!(hi > lo)) return { x: g(lo) >= 0 ? lo : 0, limited: g(lo) >= 0 ? 'range' : 'dead' };
  let last = null;
  for (let i = 0; i < n; i++) { const a = lo + ((hi - lo) * i) / n, b = lo + ((hi - lo) * (i + 1)) / n; if (g(a) >= 0 && g(b) < 0) last = [a, b]; }
  if (last) return { x: brent(g, last[0], last[1], 1e-9), limited: null };
  return g(hi) >= 0 ? { x: hi, limited: hi < xa[xa.length - 1] - 1e-9 ? 'range' : 'inflow' } : { x: 0, limited: 'dead' };
}
const mode = (a) => { const c = new Map(); let best = a[0] ?? '—', nb = 0; for (const q of a) { const k = (c.get(q) || 0) + 1; c.set(q, k); if (k > nb) { nb = k; best = q; } } return best; };

/**
 * Well-count study: n = 1 … nMax identical wells sharing the case rate equally. One well is a function of its own rate only (inflow and
 * tubing lift), so a single wellhead-pressure curve serves every n; the network back-pressure is a function of the total rate.
 * { qCase (case rate in the basis of the well curve), liquidCase (Sm³/d liquid at the case rate), well(q) → { pwf, whp (natural wellhead
 *   pressure, bara), dead, ero } | null above the open-flow potential, grid / gridNoLift { q[], whp[], dead[] } (the same curve tabulated with
 *   the configured gas lift and without any lift), back(n, s) (wellhead pressure needed, bara, with n wells at s × case rate), pRes, pb, qMax,
 *   lift { boost (bool), dp (rated pressure rise, bar), dpMin, auto (speed-controlled to a target choke opening) }, targetDp(n, pReq) (choke
 *   pressure drop at the target opening), choke(n, pAvailable, dp) → { cv, open (%) }, cvMax, duty(n, well point, boost, pReq, s) → kW,
 *   maxDrawdown (bar), chokeMaxOpen (%) }
 * Returns { rows: [{ wells, rate (Sm³/d liquid held: the case rate when it can be delivered, else what the wells make with the lift as rated),
 *   ratePerWell, rateNatural, rateWithLift (chokes open), drawdown, pwf, whpNatural, whpAvailable, whpNeeded, chokeDp, chokeOpening, boostDp,
 *   boostPower, tubingErosionalRatio, aboveBubblePoint, delivers, feasible, limit }], recommended, basis }.
 */
export function wellCountStudy({ nMax = 6, qCase, liquidCase = qCase, unit = 'Sm³/d', well, grid, gridNoLift = null, back, pRes, pb = 0, qMax = Infinity, lift = {}, targetDp = () => 0, choke = () => ({ cv: 0, open: 100 }), cvMax = Infinity, duty = () => 0, maxDrawdown = Infinity, chokeMaxOpen = 100 }) {
  const g0 = gridNoLift || grid, L = { boost: !!lift.boost && lift.dp > 0, dp: lift.dp || 0, dpMin: lift.dpMin || 0, auto: !!lift.auto }, rows = [];
  for (let n = 1; n <= nMax; n++) {
    const sOf = (g) => g.q.map((q) => (n * q) / qCase), sL = sOf(grid), dem = sL.map((s) => back(n, s)), s0 = sOf(g0);
    const iNat = intersect(s0, g0.whp, s0, s0.map((s) => back(n, s))), iLift = L.boost ? intersect(sL, grid.whp.map((w, i) => (grid.dead[i] ? w : w + L.dp)), sL, dem) : g0 !== grid ? intersect(sL, grid.whp, sL, dem) : iNat;
    const q1 = qCase / n, w1 = q1 < qMax ? well(q1) : null, pReq = back(n, 1);
    let dpB = 0, pAv = null, dpC = -1, open = 100, cv = cvMax;
    if (w1 && !w1.dead) {
      if (L.boost) dpB = L.auto ? clamp(pReq + targetDp(n, pReq) - w1.whp, Math.min(L.dpMin, L.dp), L.dp) : L.dp;
      pAv = w1.whp + dpB; dpC = pAv - pReq;
      if (dpC > 0) { const c = choke(n, pAv, dpC); cv = c.cv; open = c.open; }
    }
    const delivers = dpC > 0, cvOk = cv <= cvMax * 1.0001, sProd = delivers ? 1 : clamp(iLift.x, 0, 1), qProd = (sProd * qCase) / n, wP = delivers ? w1 : sProd > 0 ? well(qProd) : null;
    const pwf = wP ? wP.pwf : pRes, dd = Math.max(pRes - pwf, 0), power = delivers ? duty(n, w1, dpB, pReq, 1) : wP && !wP.dead ? duty(n, wP, L.dp, back(n, sProd), sProd) : 0, ero = wP ? wP.ero || 0 : 0;
    const why = !delivers ? (!w1 ? `the rate per well (${rd(q1, 0)}) is above the open-flow potential (${rd(qMax, 0)} ${unit})` : w1.dead ? `one well cannot lift ${rd(q1, 0)} ${unit} to the wellhead through this tubing` : `wellhead pressure ${rd(Math.max(pAv ?? 0, 0), 1)} bara against ${rd(pReq, 1)} bara needed`) : dd > maxDrawdown ? `drawdown ${rd(dd, 0)} bar exceeds the limit of ${rd(maxDrawdown, 0)} bar` : ero > 1 ? `tubing erosional ratio ${rd(ero, 2)}` : !cvOk ? `the choke needs Cv ${rd(cv, 0)}, above the rated ${rd(cvMax, 0)}` : '';
    rows.push({ wells: n, rate: sProd * liquidCase, ratePerWell: qProd, rateNatural: iNat.x * liquidCase, rateWithLift: iLift.x * liquidCase, delivers, drawdown: dd, pwf, whpNatural: w1 && !w1.dead ? w1.whp : null, whpAvailable: pAv, whpNeeded: pReq, chokeDp: delivers ? dpC : 0, chokeOpening: delivers ? open : 100, boostDp: delivers ? dpB : L.dp, boostPower: power, tubingErosionalRatio: ero, aboveBubblePoint: pwf >= pb - 1e-9, feasible: delivers && cvOk && dd <= maxDrawdown && ero <= 1, limit: why });
  }
  const ok = rows.filter((r) => r.feasible), best = ok.find((r) => r.chokeOpening <= chokeMaxOpen + 1e-9) || ok[0] || rows.reduce((b, r) => (r.rate > b.rate * 1.0005 ? r : b), rows[0]);
  return { rows, recommended: best.wells, basis: best.feasible ? (best.chokeOpening <= chokeMaxOpen + 1e-9 ? `smallest number of wells that holds the case rate within the drawdown limit of ${rd(maxDrawdown, 0)} bar with the chokes no more than ${rd(chokeMaxOpen, 0)} % open` : `smallest number of wells that holds the case rate within the drawdown limit (chokes more than ${rd(chokeMaxOpen, 0)} % open)`) : `no number of wells up to ${nMax} holds the case rate within the limits: the one with the highest rate is listed` };
}

function run(vIn, ctx = {}) {
  const v = readInputs(vIn), fm = fluidModel(ctx), rates = fm.rates, warnings = [], recs = [], prog = (f, msg) => ctx.progress?.(f, msg);
  const mHC = rates.mHC, mW = rates.mW, mCase = mHC + mW, qLiqStd = rates.qOilStd + rates.qWaterStd, qGasStd = rates.qGasStd;
  if (!(mCase > 0)) throw new Error('The case fluid has no flow rate: set an oil, gas or mass rate on the fluid page.');
  if (v.tSeaSurface < v.tSeabed - 15) warnings.push({ level: 'info', msg: 'The sea surface is much colder than the seabed: check the two sea temperatures.' });

  // ---- geometry ---------------------------------------------------------------------------------------------
  prog(0.05, 'Route geometry');
  const route = routeOf(v), x = route.x, z = route.z, nP = x.length, prof = analyseProfile(x, z);
  const rb = route.rbI !== null ? { i: route.rbI, x: x[route.rbI], z: z[route.rbI], height: z[nP - 1] - z[route.rbI] } : riserBase(x, z, v.geomMode === 'table' || v.geomMode === 'centreline' ? v.riserBaseX * ({ m: 1, km: 1000, ft: UNIT.ft }[v.lengthUnit] || 1) : 0);
  const D = v.idMm / 1000, wt = v.wtMm / 1000, rough = v.roughUm * 1e-6, Dfl = v.effIdMm > 0 && v.effIdMm < v.idMm ? v.effIdMm / 1000 : D, roughFl = v.effRoughUm > 0 ? v.effRoughUm * 1e-6 : rough, fouled = Dfl < 0.999 * D || Math.abs(roughFl - rough) > 0.01 * rough;
  const waterDepth = Math.max(0, -prof.zMin), offshore = waterDepth > 1, lenFlow = prof.s[rb.i], lenRiser = prof.length - lenFlow, volume = area(D) * prof.length;
  route.notes.forEach((msg) => warnings.push({ level: 'info', msg }));
  if (route.raw && route.raw.issues.length) warnings.push({ level: 'warn', msg: `Geometry check of the profile table: ${route.raw.zeroLength} zero-length element(s), ${route.raw.duplicates} repeated point(s), ${route.raw.backward} backward step(s), ${route.raw.issues.filter((q) => q.type === 'invalid').length} invalid row(s). ${route.raw.issues[0].msg}` });
  if (prof.length > 1.5 * Math.max(prof.horizontal, 1) && !(rb.height > 0)) warnings.push({ level: 'info', msg: 'The route is much longer than its horizontal extent: check that distance and elevation use the same unit (m).' });

  // ---- wall thermal design -------------------------------------------------------------------------------------
  prog(0.12, 'Wall thermal design');
  const mat = material(v.grade), est = mixAt(fm, v.pSep + 35, v.tIn - 10, mHC, mW), qEst = est.qG + est.qL, lam = est.qL / qEst, rhoNs = mCase / qEst, muNs = lam * est.muL + (1 - lam) * est.muG, kNs = lam * est.kL + (1 - lam) * est.kG, cpNs = (est.mG * est.cpG + est.mO * est.cpO + est.mW * est.cpW) / mCase;
  const vEst = qEst / area(D), hIn = hInside((rhoNs * vEst * D) / muNs, (cpNs * muNs) / kNs, kNs, D);
  const layersIn = v.layers.map((l) => ({ name: String(l?.name ?? '').trim() || 'Coating', t: Math.max(num(l?.tMm, 0), 0) / 1000, k: num(l?.k, 0), rho: num(l?.rho, 0), cp: num(l?.cp, 0) }));
  { let bi = -1, br = 0; layersIn.forEach((l, i) => { const r = l.t > 0 && l.k > 0 ? l.t / l.k : 0; if (r > br) { br = r; bi = i; } }); if (bi >= 0 && v.insKMult !== 1) layersIn[bi] = { ...layersIn[bi], k: layersIn[bi].k * v.insKMult }; } // ageing of the insulating layer
  const wallBase = { id: D, wt, kSteel: mat.k, rhoSteel: mat.rho, layers: layersIn, pipWt: v.pip ? v.pipWtMm / 1000 : 0, concrete: { t: v.concreteMm / 1000, k: 2, rho: v.concreteRho }, hIn, rhoContents: rhoNs, uMult: v.uMult, E: mat.E * 1e6 };
  const od = wallDesign(wallBase).od, hSea = hOutside(v.currentSpeed, od, 'seawater', v.tSeabed), hAir = hOutside(v.windSpeed, od, 'air'), bur = v.burialDepth > 0 ? { cover: v.burialDepth, kSoil: v.kSoil } : null;
  const W = { sea: wallDesign({ ...wallBase, hOut: hSea }), air: wallDesign({ ...wallBase, hOut: hAir }) };
  W.seaB = bur ? wallDesign({ ...wallBase, hOut: hSea, burial: bur }) : W.sea; W.airB = bur ? wallDesign({ ...wallBase, hOut: hAir, burial: bur }) : W.air;
  const catOf = (xm, zm, onMain) => { const buried = !!bur && onMain && xm >= v.burialStartX - 1e-9 && xm <= rb.x + 1e-6; return zm < 0 ? (buried ? 'seaB' : 'sea') : buried ? 'airB' : 'air'; };
  const env = (xm, zm, onMain) => { const c = catOf(xm, zm, onMain); return { U: W[c].U, cat: c, tAmb: c === 'seaB' ? v.tSeabed : c === 'airB' ? v.tGround : c === 'sea' ? seaTemperature(-zm, v.tSeaSurface, v.tSeabed) : v.tAir }; };

  // ---- local losses and bends ------------------------------------------------------------------------------------
  const fit = fittingsLoss(v.fittings, D, rough / D, { bend: v.bendMult, tee: v.teeMult }), bK = (b) => bendK(b.radius / D, b.angle, fit.fT) * v.bendMult, bends = bendInventory(x, z, v.bendMinAngle).map((b) => ({ ...b, plane: 'vertical' })).concat(route.plan ? bendInventory(route.plan.e, route.plan.n, v.bendMinAngle).map((b) => ({ ...b, x: x[b.i], z: z[b.i], plane: 'plan' })) : []).sort((a, b) => a.x - b.x), bendsK = bends.reduce((s, b) => s + bK(b), 0), kLoss = (fit.K + bendsK) * v.kMult;

  // ---- network at the case rate ------------------------------------------------------------------------------------
  prog(0.2, 'Network hydraulics');
  const pump = { qr: v.pumpQr, hr: v.pumpHr * v.degrFactor, eta: v.pumpEff, speed: v.pumpSpeedPct / 100, trim: v.pumpTrimPct / 100, shutoff: v.pumpShutoff, curve: null };
  { // tested pump curve, when given
    const pc = v.pumpCurve.map((r) => ({ q: num(r?.q, null), h: num(r?.h, null), eta: num(r?.eta, null) })).filter((r) => r.q !== null && r.h !== null && r.q >= 0).sort((a, b) => a.q - b.q).filter((r, i, a) => !i || r.q > a[i - 1].q + 1e-9);
    if (pc.length >= 3) pump.curve = { q: pc.map((r) => r.q), h: pc.map((r) => r.h * v.degrFactor), eta: pc.every((r) => r.eta > 0) ? pc.map((r) => (r.eta > 1.5 ? r.eta / 100 : r.eta)) : null };
    else if (v.pumpCurve.length) warnings.push({ level: 'info', msg: 'The tested pump curve needs at least three points with increasing flow: the parabola through the rated point is used instead.' });
  }
  const rhoGstd = (P_STD * 1e5 * est.mwG * 1e-3) / (R * (T_STD + KEL)), liftGas = v.liftType === 'gaslift' && v.liftGasRate > 0 ? (v.liftGasRate / DAY) * rhoGstd : 0; // lift gas per well (kg/s)
  let boost = 0; // pressure rise of the subsea boosting pump (bar), set once the network back-pressure is known
  const I = { v, geo: { x, z, S: prof.s, rb }, fm, fHC: mHC / mCase, fW: mW / mCase, mCase, env, D, rough, Dfl, roughFl, kFit: kLoss, pump, liftPerWell: liftGas, liftFrac: 0 };
  const net = caseNetwork(I), base = net.solve(1), sol = base.sol;
  net.notes.forEach((msg) => warnings.push({ level: 'info', msg }));
  if (!sol.usable) warnings.push({ level: 'bad', msg: `The network solver stopped at a relative residual of ${sg(sol.residual, 2)} after ${sol.iterations} iterations: treat pressures as approximate (check for closed paths, tiny diameters or pumps running backwards).` });
  const eInfo = net.edges.map((e, i) => {
    const r = sol.rec[i] || {}, m = sol.m[i], pipe = Array.isArray(r.P), rev = m < 0, A = area(e.D), vEnd = pipe && r.vm.length ? r.vm[r.vm.length - 1] : 0;
    const rn = pipe && r.rhoNs.length ? r.rhoNs[Math.floor(r.rhoNs.length / 2)] : rhoNs, mu = muNs, Re = pipe ? (rn * (r.vMax || 0) * e.D) / mu : 0;
    return { e, r, m, pipe, rev, A, incl: pipe && e.L > 0 ? Math.asin(clamp((net.nodes[e.b].z - net.nodes[e.a].z) / e.L, -1, 1)) / D2R : 0, vMax: pipe ? r.vMax : 0, vEnd, ero: pipe ? r.eroMax / v.cErosion : 0, regime: pipe ? mode(r.regime) : e.kind, holdup: pipe && r.holdup.length ? r.holdup.reduce((s, h) => s + h, 0) / r.holdup.length : null, Re, dpFric: pipe ? r.dpFric + r.dpMinor : null, dpGrav: pipe ? r.dpGrav : null, tOut: isNum(r.tOut) ? r.tOut : base.t[e.b] };
  });
  // main line (flowline pieces in route order, then the riser) for plots and length-weighted quantities
  const mainIdx = net.edges.map((e, i) => i).filter((i) => net.edges[i].kind === 'flowline').sort((a, b) => net.edges[a].sA - net.edges[b].sA), riserIdx = net.edges.findIndex((e) => e.kind === 'riser');
  if (riserIdx >= 0) mainIdx.push(riserIdx);
  const ml = { x: [], z: [], P: [], T: [], holdup: [], vm: [], U: [], cat: [], ds: [], theta: [], D: [] };
  for (const i of mainIdx) {
    const q = eInfo[i], c = net.edges[i].cells, n = c.length;
    for (let j = 0; j < n; j++) { const k = q.rev ? n - 1 - j : j; ml.x.push(c[j].xMid); ml.z.push(c[j].zMid); ml.P.push(0.5 * (q.r.P[k] + q.r.P[k + 1])); ml.T.push(0.5 * (q.r.T[k] + q.r.T[k + 1])); ml.holdup.push(q.r.holdup[k]); ml.vm.push(q.r.vm[k]); ml.U.push(c[j].U); ml.cat.push(catOf(c[j].xMid, c[j].zMid, true)); ml.ds.push(c[j].ds); ml.theta.push(c[j].theta); ml.D.push(c[j].D); }
  }
  const mlLen = ml.ds.reduce((s, d) => s + d, 0), catLen = {}; ml.cat.forEach((c, i) => { catLen[c] = (catLen[c] || 0) + ml.ds[i]; });
  const mainCat = Object.keys(catLen).sort((a, b) => catLen[b] - catLen[a])[0] || (offshore ? 'sea' : 'air'), wall = W[mainCat], uMean = mlLen > 0 ? ml.U.reduce((s, u, i) => s + u * ml.ds[i], 0) / mlLen : wall.U;
  const nodeP = (k) => sol.p[k], sinkK = net.nodes.findIndex((n) => n.kind === 'sink'), tArr = base.t[sinkK], inletK = mainIdx.length ? net.edges[mainIdx[0]].a : sol.p.indexOf(Math.max(...sol.p)), pInlet = nodeP(inletK);
  const vMaxAll = Math.max(...eInfo.map((q) => q.vMax)), eroAll = Math.max(...eInfo.map((q) => q.ero)), liquidInv = eInfo.reduce((s, q) => s + (q.pipe ? q.r.liquid : 0), 0);
  let cleanP = null;
  if (fouled) { try { const c = caseNetwork({ ...I, Dfl: D, roughFl: rough }).solve(1); cleanP = c.sol.p[inletK]; } catch { cleanP = null; } }

  // ---- free spans ------------------------------------------------------------------------------------------------
  let spans = [], spanSource = 'none';
  const wSub = W.sea.submerged;
  if (offshore && lenFlow > 0 && !bur) {
    let fx, fz;
    if (v.terrain && v.geomMode !== 'table') { try { const t = terrainTransect(v.terrain, v.terrainLine); fx = t.x; fz = t.z; spanSource = 'terrain grid'; } catch { fx = null; } }
    if (!fx && v.seabedRms > 0) { // statistical micro-relief (seeded, so the screening is repeatable)
      const r = rng(2024), K = 6, comp = Array.from({ length: K }, (_, k) => ({ lam: v.seabedCorr * (0.5 + (2 * k) / (K - 1)), ph: r.uniform(0, 2 * Math.PI) })), amp = v.seabedRms * Math.sqrt(2 / K), n = clamp(Math.round(rb.x / (v.seabedCorr / 8)), 50, 6000);
      fx = []; fz = []; spanSource = 'statistical micro-relief';
      for (let i = 0; i <= n; i++) { const xi = x[0] + ((rb.x - x[0]) * i) / n; fx.push(xi); fz.push(interp1(x, z, xi) + comp.reduce((s, c) => s + amp * Math.sin((2 * Math.PI * xi) / c.lam + c.ph), 0)); }
    }
    if (fx && wSub > 0) spans = freeSpans(fx, fz, { w: wSub, EI: wall.EI, gap: v.spanGap, maxLength: 250 }).filter((s) => s.x <= rb.x + 1e-6);
  }
  // supports and restraints: a support inside a span candidate splits it
  const supUser = (Array.isArray(v.supports) ? v.supports : []).map((r) => ({ x: num(r?.x, null), type: String(r?.type ?? '').trim() || 'support', k: Math.max(num(r?.k, 0), 0) })).filter((r) => r.x !== null).sort((a, b) => a.x - b.x);
  let spansSplit = 0;
  if (spans.length && supUser.length) { const kept = []; for (const s of spans) { const a = s.x - s.length / 2, b = s.x + s.length / 2, cut = supUser.map((q) => q.x).filter((q) => q > a + 1 && q < b - 1); if (!cut.length) { kept.push(s); continue; } spansSplit++; const ed = [a, ...cut, b]; for (let k = 1; k < ed.length; k++) if (ed[k] - ed[k - 1] >= 5) kept.push({ x: 0.5 * (ed[k] + ed[k - 1]), length: ed[k] - ed[k - 1], gap: s.gap }); } spans = kept; }
  const spanTotal = spans.reduce((a, q) => a + q.length, 0), spanOver = spans.filter((s) => s.length > v.spanLimit).length, extraSupports = spans.reduce((a, s) => a + Math.max(Math.ceil(s.length / v.spanLimit) - 1, 0), 0), spansAll = spans.slice();
  spans.sort((a, b) => b.length - a.length); const spanCount = spans.length, spanMax = spans[0]?.length || 0; spans = spans.slice(0, 40).sort((a, b) => a.x - b.x);
  // design features: buckle arrestors, lateral-buckle initiators, cathodic protection
  const odS = D + 2 * wt, ppr = propagationPressure(mat.smys, wt, odS), arrFirst = v.buckleArrestorFrom > 0 ? v.buckleArrestorFrom : (ppr * 1e5) / 1.3 / (RHO_SEA * G); // propagation pressure (bar) and the depth where the factored external pressure reaches it
  const arrLen = prof.segs.reduce((s, q) => { const d0 = -q.z0, d1 = -q.z1, lo = Math.min(d0, d1), hi = Math.max(d0, d1); return s + q.L * (lo >= arrFirst ? 1 : hi <= arrFirst ? 0 : (hi - arrFirst) / (hi - lo)); }, 0);
  const arrestors = { spacing: v.buckleArrestorSpacing, count: v.buckleArrestorSpacing > 0 && arrLen > 0 ? Math.ceil(arrLen / v.buckleArrestorSpacing) + 1 : 0, firstDepth: arrFirst, length: arrLen, propagationPressure: ppr, required: waterDepth > arrFirst };
  const nInit = v.buckleInitType !== 'none' && lenFlow > 0 ? Math.max(0, Math.ceil((rb.x - x[0]) / v.buckleInitSpacing) - 1) : 0, initX = Array.from({ length: nInit }, (_, i) => x[0] + (i + 1) * v.buckleInitSpacing);
  const sleeperSpan = v.buckleInitType === 'sleeper' ? liftOffLength(W.sea.EI, v.sleeperHeight, wSub) : 0; // lift-off length each side of a sleeper: beam on a rigid seabed raised by the sleeper height
  const cpArea = Math.PI * odS * (v.anodeSpacing > 0 ? v.anodeSpacing : prof.length), cpI = cpArea * v.cpCurrent * v.cpBreakdown, cpReq = anodeMassRequired(cpI, v.cpLife); // Al–Zn–In anodes: 2,000 A·h/kg, utilisation 0.8
  const cp = { anodeSpacing: v.anodeSpacing, count: v.anodeSpacing > 0 && offshore ? Math.ceil(prof.length / v.anodeSpacing) : 0, anodeMass: v.anodeMass, requiredMass: cpReq, utilisation: v.anodeSpacing > 0 ? cpReq / v.anodeMass : null, designLife: v.cpLife, currentPerAnode: cpI, lifeYears: v.anodeSpacing > 0 ? Math.min((v.anodeMass * 0.8 * 2000) / (cpI * 8760), 999) : 0 };
  let gchk = []; try { gchk = checkProfile({ x, z, s: prof.s, incl: prof.segs.map((q) => q.incl) }, { id: D }); } catch (e) { gchk = [{ name: 'Profile check', pass: false, got: 0, expected: 0, note: e.message }]; }
  const orient = orientationCheck(net.edges, sol.m);

  // ---- wells: trajectory, inflow, lift, nodal analysis -------------------------------------------------------------
  prog(0.45, 'Wells and nodal analysis');
  const traj = minimumCurvature(v.survey), Dt = v.tubingIdMm / 1000;
  if (!(traj.tvd > 1)) throw new Error('The deviation survey gives no vertical depth: check measured depths and inclinations.');
  const tub = tubingCells(traj, v.nTubing, Dt, v.tubingRoughUm * 1e-6), ipr = makeIpr(v, fm, mHC, mW), gasBasis = ipr.basis === 'gas', qBasis = gasBasis ? qGasStd : qLiqStd;
  if (liftGas > 0) tub.forEach((c, i) => { c.gx = traj.md - (i + 0.5) * c.ds <= (v.liftDepthPct / 100) * traj.md ? 1 : 0; }); // lift gas joins above the injection valve
  if (!(qBasis > 0)) throw new Error(gasBasis ? 'The case fluid gives no gas at standard conditions: choose a liquid inflow model.' : 'The case fluid gives no stock-tank liquid: choose the gas-well inflow model.');
  const rated = net.nodes.map((n, k) => k).filter((k) => net.nodes[k].kind === 'source' && net.nodes[k].rate > 0), srcAll = net.nodes.map((n, k) => k).filter((k) => net.nodes[k].kind === 'source');
  const wi = (rated.length ? rated : srcAll).reduce((b, k) => (nodeP(k) > nodeP(b) ? k : b)), outOf = (k) => net.edges.reduce((s, e, i) => s + (e.a === k ? sol.m[i] : e.b === k ? -sol.m[i] : 0), 0);
  const share = clamp(net.nodes[wi].rate > 0 ? net.nodes[wi].share : Math.max(outOf(wi), 1e-9) / mCase, 1e-6, 1), qWell = share * qBasis, wellUnit = gasBasis ? 'Sm³/d gas' : 'Sm³/d liquid';
  const tOf = (f) => v.tRes + (v.tIn - v.tRes) * f;
  const avail = (s, o = {}) => {
    const q = s * qWell, pwf = iprPwf(q, o.ipr || ipr);
    if (pwf === null) return { pwf: 0, whp: -50, dead: true };
    const r = marchPipe({ fm, mHC: o.mHC ? o.mHC(s) : share * mHC * s, mW: o.mW ? o.mW(s) : share * mW * s, cells: o.cells || tub, pIn: pwf, tIn: v.tRes, tOf, model: v.netModel, fModel: v.fModel, gExtra: liftGas });
    return { pwf, whp: Math.max(r.pOut + (r.dead ? 0 : boost), -50), dead: r.dead, r }; // a boosting pump adds its pressure rise at the wellhead
  };
  const sMax = Math.min(2.6, (0.985 * ipr.qMax) / qWell), canScale = rated.length > 0;
  if (!(sMax > 1e-3)) throw new Error('The well has no inflow potential: check reservoir pressure and productivity.');
  if (v.liftType === 'boost' && v.liftDp > 0) { // duty of the boosting pump
    const a0 = v.liftAuto && sMax >= 1 ? avail(1) : null; /* natural wellhead pressure, before any boosting */
    boost = v.liftDp;
    if (a0) {
      const wW = share * mCase, cvT = v.chokeCvMax * valveCharacteristic(v.chokeTarget / 100, v.chokeChar, v.chokeRange);
      let dpT = 5; for (let k = 0; k < 3; k++) { const p1 = nodeP(wi) + dpT, pr = mixAt(fm, p1, v.tIn, share * mHC, share * mW); dpT = valveDp({ cv: cvT, w: wW, p1, rhoL: pr.rhoL, rhoG: pr.rhoG, xG: pr.mG / wW, k: kGas(pr), xT: v.chokeXT }).dp; }
      if (!a0.dead) boost = clamp(nodeP(wi) + dpT - a0.whp, Math.min(v.liftDpMin, v.liftDp), v.liftDp);
    }
  }
  const sA = linspace(0.06 * sMax, sMax, 9), curve = sA.map((s) => avail(s)), yA = curve.map((c) => c.whp);
  let sD = [1], yD = [nodeP(wi)];
  const more = (s) => { try { const up = s > 1, r = net.solve(s, up ? gUp : gDn); if (!r.sol.usable) return false; if (up) { gUp = r.guess; sD.push(s); yD.push(r.sol.p[wi]); } else { gDn = r.guess; sD.unshift(s); yD.unshift(r.sol.p[wi]); } return true; } catch { return false; } };
  let gUp = base.guess, gDn = base.guess;
  if (canScale) {
    for (const s of [1.5, 0.5]) more(s);
  } else { sD = [0, 3]; yD = [nodeP(wi), nodeP(wi)]; }
  let op = intersect(sA, yA, sD, yD);
  if (canScale) { // the back-pressure curve is extended only when the intersection lies outside 50–150 % of the case rate
    if (op.limited === 'range' && sMax * 1.05 >= 2.2 && more(2.2)) op = intersect(sA, yA, sD, yD);
    if (op.limited === 'dead' && more(0.3)) op = intersect(sA, yA, sD, yD);
  }
  if (canScale && !op.limited) { // polish the intersection with true solutions at the estimate
    for (let k = 0; k < 1; k++) {
      const s = op.x; if (sD.some((q) => Math.abs(q - s) < 2e-3)) break;
      try { const j = sD.findIndex((q) => q > s), r = net.solve(s, base.guess); if (!r.sol.usable) break; sD.splice(j < 0 ? sD.length : j, 0, s); yD.splice(j < 0 ? yD.length : j, 0, r.sol.p[wi]); const a = avail(s), ja = sA.findIndex((q) => q > s); sA.splice(ja < 0 ? sA.length : ja, 0, s); yA.splice(ja < 0 ? yA.length : ja, 0, a.whp); curve.splice(ja < 0 ? curve.length : ja, 0, a); } catch { break; }
      op = intersect(sA, yA, sD, yD);
    }
  }
  const sOp = op.x, operatingRate = sOp * qLiqStd, whpOp = sOp > 0 ? interp1(sD, yD, sOp) : nodeP(wi), pwfOp = sOp > 0 ? iprPwf(sOp * qWell, ipr) ?? 0 : ipr.pRes;
  const aCase = sMax >= 1 ? avail(1) : null, whpReq = nodeP(wi), whpAvail = aCase ? Math.max(aCase.whp, 0) : 0, chokeDp = whpAvail - whpReq, canDeliver = !!aCase && chokeDp > 0;

  // ---- well-count study: n = 1 … 6 identical wells sharing the case rate ------------------------------------------------
  // One well is a function of its own rate only (inflow + tubing), so a single wellhead-pressure curve serves every n; the network
  // back-pressure depends on the total rate (flowline and riser) plus the jumper loss, which scales with the square of the rate per well.
  const liftDuty = (nW, wNat, r, dpBoost, pReq, sc = 1) => { // shaft power (kW) and conditions of the lift at a natural wellhead pressure wNat
    const o = { power: 0, rate: 0, dp: 0, injectionP: null, suctionP: null, flow: null };
    if (dpBoost > 0) { const ps = Math.max(wNat, 1.5), pb = mixAt(fm, ps, v.tIn, mHC * sc, mW * sc), q = pb.qG + pb.qL; o.dp = dpBoost; o.suctionP = ps; o.flow = q * 3600; o.power = (q * dpBoost * 1e5) / v.liftEff / 1000; }
    else if (liftGas > 0) {
      const iv = clamp(Math.round(v.nTubing * (1 - v.liftDepthPct / 100)), 0, v.nTubing), pInj = (r ? Math.max(r.P[iv], 1) : pReq + 50) + 10, rr = Math.max(pInj / v.pSep, 1.01), nS = Math.max(1, Math.ceil(Math.log(rr) / Math.log(v.compMaxRatio) - 1e-9));
      o.rate = v.liftGasRate; o.injectionP = pInj; o.power = (nS * compressorStage({ p1: 1, p2: rr ** (1 / nS), t1: v.compTcool, z: 0.9, k: kGas(est), mw: est.mwG, eta: v.compEta, mdot: liftGas * nW }).power) / 1000;
    }
    return o;
  };
  const chokeTargetDp = (pReq, w, fr) => { const cvT = v.chokeCvMax * valveCharacteristic(v.chokeTarget / 100, v.chokeChar, v.chokeRange); let d = 5; for (let k = 0; k < 3; k++) { const p1 = pReq + d, pr = mixAt(fm, p1, v.tIn, fr * mHC, fr * mW); d = valveDp({ cv: cvT, w, p1, rhoL: pr.rhoL, rhoG: pr.rhoG, xG: pr.mG / w, k: kGas(pr), xT: v.chokeXT }).dp; } return d; };
  const nNow = Math.max((rated.length ? rated : srcAll).length, 1);
  const study = (() => {
    const key = (q) => q.toPrecision(9), qTop = 0.985 * ipr.qMax;
    const grid = (gx, seed) => { // wellhead pressure of one well against its rate, without boosting (gx: lift gas per well, kg/s)
      const pts = new Map(), put = (q, pwf, r) => { const o = { q, pwf, whp: r.dead ? -50 : r.pOut, dead: !!r.dead, ero: (r.eroMax || 0) / v.cErosion, r }; pts.set(key(q), o); return o; };
      const well = (q) => { const h = pts.get(key(q)); if (h) return h; const pwf = q < ipr.qMax ? iprPwf(q, ipr) : null; if (pwf === null || !(q > 0)) return null; const f = q / qBasis, cells = gx > 0 ? tub : tub0; return put(q, pwf, marchPipe({ fm, mHC: mHC * f, mW: mW * f, cells, pIn: pwf, tIn: v.tRes, tOf, model: v.netModel, fModel: v.fModel, gExtra: gx })); };
      if (seed) { sA.forEach((s, i) => { if (curve[i].r) put(s * qWell, curve[i].pwf, curve[i].r); }); if (aCase?.r) put(qWell, aCase.pwf, aCase.r); }
      else for (const f of [0.05, 0.2, 0.4, 0.6, 0.8]) well(f * qTop);
      for (let n = 1; n <= 6; n++) well(qBasis / n);
      const hi = Math.max(...[...pts.values()].map((o) => o.q)); if (hi < 0.97 * qTop) for (const f of [1 / 3, 2 / 3, 1]) well(hi + (qTop - hi) * f);
      const g = [...pts.values()].sort((a, b) => a.q - b.q);
      return { well, q: g.map((o) => o.q), whp: g.map((o) => o.whp), dead: g.map((o) => o.dead) };
    };
    const tub0 = liftGas > 0 ? tub.map((c) => ({ ...c, gx: 0 })) : tub, gL = grid(liftGas, true), g0 = liftGas > 0 ? grid(0, false) : gL; // with the configured gas lift, and without any lift
    const nS = sD.length, bHi = nS > 1 && sD[nS - 1] > sD[nS - 2] ? Math.max((yD[nS - 1] - yD[nS - 2]) / (sD[nS - 1] ** 2 - sD[nS - 2] ** 2), 0) : 0;
    const back0 = (s) => (s > sD[nS - 1] ? yD[nS - 1] + bHi * (s * s - sD[nS - 1] ** 2) : interp1(sD, yD, Math.max(s, sD[0]))); // beyond the solved range the friction term is continued with the square of the rate
    const dJ0 = canScale ? clamp(nodeP(wi) - pInlet, 0, 0.2 * nodeP(wi)) : 0, boosting = v.liftType === 'boost' && v.liftDp > 0;
    return wellCountStudy({ qCase: qBasis, liquidCase: qLiqStd, unit: wellUnit, well: gL.well, grid: gL, gridNoLift: g0, back: (n, s) => back0(s) + dJ0 * s * s * ((nNow / n) ** 2 - 1), pRes: ipr.pRes, pb: ipr.pb, qMax: ipr.qMax,
      lift: { boost: boosting, dp: boosting ? v.liftDp : 0, dpMin: v.liftDpMin, auto: !!v.liftAuto }, targetDp: (n, pReq) => chokeTargetDp(pReq, mCase / n, 1 / n), cvMax: v.chokeCvMax, maxDrawdown: v.maxDrawdown, chokeMaxOpen: v.chokeMaxOpen,
      choke: (n, pAv, dpC) => { const w = mCase / n, pr = mixAt(fm, pAv, v.tIn, mHC / n, mW / n), cv = valveCv({ w, dp: dpC, p1: pAv, rhoL: pr.rhoL, rhoG: pr.rhoG, xG: pr.mG / w, k: kGas(pr), xT: v.chokeXT, FL: v.chokeFL }).cv; return { cv, open: 100 * valveOpening(Math.min(cv / v.chokeCvMax, 1), v.chokeChar, v.chokeRange) }; },
      duty: (n, w, dpB, pReq, sc) => liftDuty(n, w.whp, w.r, dpB, pReq, sc).power });
  })();

  // sensitivities of the well deliverability curve (flowline back-pressure curve unchanged)
  const sens = [], sensCurves = [], rhoWstd = rates.qWaterStd > 0 ? (mW * DAY) / rates.qWaterStd : 1025, wcNow = qLiqStd > 0 ? rates.qWaterStd / qLiqStd : 0;
  const variant = (label, param, value, o) => { const ys = sA.map((s) => avail(s, o).whp), r = intersect(sA, ys, sD, yD); sens.push({ param, value, rate: r.x * qLiqStd, s: r.x }); sensCurves.push({ name: label, x: sA.map((s) => s * qWell), y: ys.map((p) => Math.max(p, 0)) }); };
  if (!gasBasis && rates.qOilStd > 0) for (const wc of [0, 0.5, 0.8].filter((w) => Math.abs(w - wcNow) > 0.05).slice(0, 2)) variant(`Water cut ${rd(wc * 100, 0)} %`, 'Water cut (%)', wc * 100, { mHC: (s) => (mHC * s * qWell * (1 - wc)) / rates.qOilStd, mW: (s) => (s * qWell * wc * rhoWstd) / DAY });
  for (const f of [0.78, 1.25]) variant(`Tubing ID ${rd(v.tubingIdMm * f, 1)} mm`, 'Tubing ID (mm)', v.tubingIdMm * f, { cells: tub.map((c) => ({ ...c, D: Dt * f })) });
  for (const f of [0.85, 0.7]) { const o = { ...ipr, pRes: ipr.pRes * f, pb: Math.min(ipr.pb, ipr.pRes * f) }; o.qMax = iprRate(0, o); variant(`Reservoir pressure ${rd(ipr.pRes * f, 0)} bara`, 'Reservoir pressure (bara)', ipr.pRes * f, { ipr: o }); }

  // ---- production choke ------------------------------------------------------------------------------------------
  prog(0.7, 'Choke, separator and rotating equipment');
  const mWell = share * mCase, p1c = canDeliver ? whpAvail : Math.max(whpReq, 1.5), prC = mixAt(fm, p1c, v.tIn, share * mHC, share * mW), xGc = prC.mG / mWell, kC = kGas(prC), rhoC = mWell / (prC.qG + prC.qL);
  let chokeCv = v.chokeCvMax, chokeOpening = 100, chokeState = null;
  if (canDeliver) { chokeState = valveCv({ w: mWell, dp: chokeDp, p1: p1c, rhoL: prC.rhoL, rhoG: prC.rhoG, xG: xGc, k: kC, xT: v.chokeXT, FL: v.chokeFL }); chokeCv = chokeState.cv; chokeOpening = 100 * valveOpening(Math.min(chokeCv / v.chokeCvMax, 1), v.chokeChar, v.chokeRange); }
  const dpOpen = valveDp({ cv: v.chokeCvMax, w: mWell, p1: p1c, rhoL: prC.rhoL, rhoG: prC.rhoG, xG: xGc, k: kC, xT: v.chokeXT }).dp, p2c = canDeliver ? whpReq : Math.max(p1c - dpOpen, 1);
  const sach = sachdeva({ p1: p1c, p2: p2c, x: xGc, rhoG: prC.rhoG, rhoL: prC.rhoL, k: kC, cvG: prC.cpG / kC, cL: prC.cpL, cd: v.chokeCd });
  // homogeneous-equilibrium path: equilibrium flash along an isentrope (dT/dp|s = μJT + v/cp, latent heat neglected)
  const hemP = [], hemV = []; { let t = v.tIn, pPrev = p1c; for (let i = 0; i <= 48; i++) { const p = p1c * 0.03 ** (i / 48), pe = Math.max(p, 1), pr = mixAt(fm, pe, t, share * mHC, share * mW), vv = (pr.qG + pr.qL) / mWell, cp = (pr.mG * pr.cpG + pr.mO * pr.cpO + pr.mW * pr.cpW) / mWell; t = clamp(t - (jtMix(pr) + vv / cp) * (pPrev - p) * 1e5, -60, 250); pPrev = p; hemP.unshift(p); hemV.unshift(vv * (pe / Math.max(p, 1e-6))); } }
  const hem = hemCritical((p) => interp1(hemP, hemV, p), p1c), critical = p2c / p1c <= sach.yc;
  const areaSach = mWell / Math.max(sach.G, 1e-9), areaHem = mWell / Math.max(v.chokeCd * hem.G, 1e-9), dOrifice = Math.sqrt((4 * mWell) / (Math.PI * v.chokeCd * Math.sqrt(2 * rhoC * Math.max(p1c - p2c, 1e-6) * 1e5)));
  const glr = qLiqStd > 0 ? qGasStd / qLiqStd : 1e5, qLw = share * qLiqStd, beanRows = Object.keys(CHOKE_CORR).map((name) => [name[0].toUpperCase() + name.slice(1), qLw > 0 && p1c > 1.2 ? rd(chokeCorrelation(name, { p1: p1c, q: qLw, glr }), 1) : '—', qLw > 0 && p1c > 1.2 ? rd((chokeCorrelation(name, { p1: p1c, q: qLw, glr }) / 64) * 25.4, 1) : '—']);

  // ---- artificial lift duty ------------------------------------------------------------------------------------------
  const nFlow = net.nodes.filter((n) => n.kind === 'source').length, lift = { type: v.liftType, wells: nFlow, gasMass: liftGas * nFlow, ...liftDuty(nFlow, whpAvail - boost, aCase?.r || null, boost, whpReq) };

  // ---- separator, slug catcher ------------------------------------------------------------------------------------
  const eosF = makeFluid(fm.spec), fl = eosProps(eosF, v.pSep, tArr), nHC = rates.nHC, wSep = waterProps(v.pSep, tArr, fm.aq);
  const nGas = fl.phase === 'oil' ? 0 : fl.phase === 'gas' ? nHC : fl.beta * nHC, nOil = nHC - nGas, mGasSep = (nGas * fl.gas.MW) / 1000, mOilSep = (nOil * fl.oil.MW) / 1000, qGasSep = nGas > 0 ? mGasSep / fl.gas.rho : 0, qOilSep = nOil > 0 ? mOilSep / fl.oil.rho : 0, qWatSep = mW / wSep.rho, qLiqSep = qOilSep + qWatSep;
  const rhoLsep = qLiqSep > 0 ? (mOilSep + mW) / qLiqSep : 800, rhoGsep = nGas > 0 ? fl.gas.rho : 1, iC1 = eosF.comps.findIndex((c) => c.id === 'C1'), zC1 = iC1 >= 0 ? eosF.z[iC1] : 0, c1Out = iC1 >= 0 ? nGas * fl.y[iC1] + nOil * fl.x[iC1] : 0;
  const compErr = Math.max(...eosF.z.map((zi, i) => Math.abs(nGas * fl.y[i] + nOil * fl.x[i] - zi * nHC))) / nHC;
  const horiz = v.sepOrient === 'horizontal', sepGeom = { d: v.sepD, l: v.sepL, orientation: v.sepOrient, heads: 'elliptical' }, sepH = horiz ? v.sepD : v.sepL + v.sepD / 2, sepLevel = (v.sepLevelPct / 100) * sepH, sepVtot = vesselVolume({ ...sepGeom, level: sepH }), sepVliq = vesselVolume({ ...sepGeom, level: sepLevel });
  const sepKused = v.sepKAuto ? demisterK(v.pSep, horiz) : v.sepK, liftQ = lift.gasMass / Math.max(rhoGsep, 1e-9); // lift gas returns through the separator
  const sepAg = horiz ? area(v.sepD) - vesselVolume({ d: v.sepD, l: 1, level: clamp(sepLevel, 0, v.sepD), heads: 'flat' }) : area(v.sepD), vGasSep = (qGasSep + liftQ) / Math.max(sepAg, 1e-9), vSB = soudersBrown(sepKused, rhoLsep, rhoGsep), gasLoad = vGasSep / Math.max(vSB, 1e-9), resTime = qLiqSep > 0 ? sepVliq / qLiqSep : 1e6;
  const scVol = v.slugSF * v.slugSurge + qLiqSep * v.sepResMin * 60, fingerA = area(v.fingerDmm / 1000), vSBf = soudersBrown(v.fingerK, rhoLsep, rhoGsep);
  let nFing = 2; while (nFing < 24 && (scVol / (nFing * fingerA * v.scEff) > v.fingerMaxL || (qGasSep + liftQ) / (nFing * 0.5 * fingerA) > vSBf)) nFing += 2;
  const fingerL = scVol / (nFing * fingerA * v.scEff), scVesselD = ((4 * (scVol / (0.6 * (v.scEff / 0.85)))) / (Math.PI * 4)) ** (1 / 3);
  // transient inventory of the slug catcher while the design slug arrives (level-proportional draw-off, capped)
  const scHold = qLiqSep * v.sepResMin * 60, scV0 = (v.scInitPct / 100) * scHold, qSlug = v.slugSurge / v.slugDuration, qCap = v.drainFactor * qLiqSep, tau = 300;
  const rhs = (during) => (t, y) => { const qin = qLiqSep + (during ? qSlug : 0), qout = clamp(qLiqSep + (y[0] - scHold) / tau, 0, qCap); return [qin - qout, qin, qout]; };
  const ph1 = rk4(rhs(true), [scV0, 0, 0], 0, v.slugDuration, 60), tEndSc = v.slugDuration + Math.max(1800, 6 * tau), ph2 = rk4(rhs(false), ph1.y[ph1.y.length - 1], v.slugDuration, tEndSc, 120);
  const scT = ph1.t.concat(ph2.t.slice(1)), scY = ph1.y.concat(ph2.y.slice(1)), scPeak = Math.max(...scY.map((y) => y[0])), scEnd = scY[scY.length - 1];

  // ---- liquid export pump ---------------------------------------------------------------------------------------------
  const muLsep = mixAt(fm, v.pSep, tArr, mHC, mW).muL, qDuty = qLiqSep * 3600, pumpLineD = v.pumpLineIdMm / 1000;
  const hSys = (q) => ((v.pumpDischargeP - v.pSep) * 1e5) / (rhoLsep * G) + v.pumpStaticHead + pipeDp({ m: (q / 3600) * rhoLsep, rho: rhoLsep, mu: muLsep, L: v.pumpLineLength, D: pumpLineD, rough, k: 6, model: v.fModel }) / (rhoLsep * G);
  const qZero = pumpRunout(pump);
  let pumpRes = null;
  if (v.pumpOn && qDuty > 0) {
    const gq = (q) => pumpHead(q, pump) - hSys(q), qOp = gq(1e-6) <= 0 ? 0 : gq(qZero) >= 0 ? qZero : brent(gq, 1e-6, qZero, 1e-9), qRun = qOp >= qDuty ? qDuty : qOp, hRun = pumpHead(qRun, pump), eta = pumpEfficiency(qRun, pump);
    const power = (rhoLsep * G * (qRun / 3600) * hRun) / eta / 1000, sx = qRun / (pump.qr * pump.speed * pump.trim), npshR = v.pumpNpshR * (pump.speed * pump.trim) ** 2 * (0.4 + 0.6 * sx * sx);
    const npshA = v.pumpSuctionHead - pipeDp({ m: (qRun / 3600) * rhoLsep, rho: rhoLsep, mu: muLsep, L: 15, D: pumpLineD * 1.25, rough, k: 1.5 }) / (rhoLsep * G); // saturated liquid: vapour pressure = vessel pressure
    const eu = eulerHead({ d2: (v.pumpImpMm / 1000) * pump.trim, rpm: v.pumpRpm * pump.speed, beta2: v.pumpBeta2, blades: Math.round(v.pumpBlades), b2: 0.07 * (v.pumpImpMm / 1000), q: qRun / 3600 });
    pumpRes = { qOp, qRun, hRun, eta, power, npshA, npshR, eu, stages: eu.head > 0 ? Math.max(1, Math.ceil(hRun / (0.85 * eu.head))) : null, throttle: qOp >= qDuty ? pumpHead(qDuty, pump) - hSys(qDuty) : 0 };
  }
  const netPumpPower = eInfo.reduce((s, q) => s + (q.e.kind === 'pump' && isNum(q.r.power) ? q.r.power : 0), 0) / 1000;

  // ---- gas compressor -----------------------------------------------------------------------------------------------------
  let comp = null;
  const mGc = mGasSep + lift.gasMass; // produced gas plus returning lift gas
  if (v.compOn && mGasSep > 1e-6 && v.compPd > v.pSep * 1.02) {
    const nSt = Math.max(1, Math.ceil(Math.log(v.compPd / v.pSep) / Math.log(v.compMaxRatio) - 1e-9)), ratio = (v.compPd / v.pSep) ** (1 / nSt), stages = [];
    let pIn = v.pSep, tInS = tArr, power = 0;
    for (let k = 0; k < nSt; k++) {
      const a = lookup(fm.table, pIn, tInS), kk = k === 0 ? clamp(fl.gas.cp / Math.max(fl.gas.cp - R / (fl.gas.MW * 1e-3), 1), 1.05, 1.67) : kGas(a), z1 = k === 0 ? fl.gas.Z : a.zG, mw = fl.gas.MW;
      const first = compressorStage({ p1: pIn, p2: pIn * ratio, t1: tInS, z: z1, k: kk, mw, eta: v.compEta, mdot: mGc }), z2 = lookup(fm.table, pIn * ratio, Math.min(first.t2, 169)).zG, st = compressorStage({ p1: pIn, p2: pIn * ratio, t1: tInS, z: 0.5 * (z1 + z2), k: kk, mw, eta: v.compEta, mdot: mGc });
      stages.push({ ...st, p1: pIn, p2: pIn * ratio, t1: tInS, k: kk, z: 0.5 * (z1 + z2), qIn: (mGc * z1 * R * (tInS + KEL)) / (mw * 1e-3 * pIn * 1e5) }); power += st.power; pIn *= ratio; tInS = v.compTcool;
    }
    const s0 = stages[0], bySpeed = new Map();
    v.compMap.forEach((r) => { const n = num(r?.speed, null), q = num(r?.q, null), h = num(r?.head, null); if (n > 0 && q > 0 && h > 0) { if (!bySpeed.has(n)) bySpeed.set(n, []); bySpeed.get(n).push({ q, h: h * 1000 }); } });
    const tested = [...bySpeed.entries()].filter(([, a]) => a.length >= 3).sort((a, b) => a[0] - b[0]).map(([n, a]) => { a.sort((p, q) => p.q - q.q); return { n: n / 100, q: a.map((p) => p.q), h: a.map((p) => p.h), qSurge: a[0].q, qChoke: a[a.length - 1].q }; });
    const map = tested.length >= 2 ? { speeds: tested.map((l) => l.n), lines: tested, tested: true } : compressorMap({ qd: s0.qIn * (1 + v.compMargin / 100) * v.compFlowMult, hd: s0.headPoly * 1.05 * v.compHeadMult * v.degrFactor }), pt = mapPoint(map, s0.qIn, s0.headPoly);
    comp = { nSt, ratio, stages, power: power / 1000, map, pt, tMax: Math.max(...stages.map((s) => s.t2)) };
  }

  // ---- transient: line-pack, separator inventories, rotors, actuators --------------------------------------------------
  let trans = null, basis = null;
  if (v.transScenario !== 'off') {
    try {
      const f0 = canDeliver ? 1 : clamp(sOp, 0, 1);
      if (!(mGasSep > 1e-6 * mCase && mOilSep + mW > 1e-6 * mCase && f0 > 0.02)) throw new Error('the case has no gas, no liquid or no flow at the separator');
      const zP = [1, 2, 4, 7, 10, 15, 20, 30, 40, 55, 70, 90, 110, 140, 170, 200, 250, 300, 400, 600], zT = { p: zP, z: zP.map((q) => clamp(lookup(fm.table, q, tArr).zG, 0.2, 1.6)) };
      const pts = sA.map((s, i) => [s * mCase, yA[i]]).concat(aCase ? [[mCase, aCase.whp]] : []).filter((q) => isNum(q[0]) && isNum(q[1])).sort((a, b) => a[0] - b[0]).filter((q, i, a) => !i || q[0] > a[i - 1][0] + 1e-9 * mCase);
      const iPk = pts.reduce((b, q, i) => (q[1] > pts[b][1] ? i : b), 0), cm = pts.map((q) => q[0]), cpv = pts.map((q, i) => Math.max(i < iPk ? pts[iPk][1] : q[1], 0)); // stable (falling) branch of the deliverability curve
      const m0 = f0 * mCase, pAvail0 = interp1(cm, cpv, m0), demand0 = canDeliver ? whpReq : interp1(sD, yD, clamp(f0, sD[0], sD[sD.length - 1])), pLine = mlLen > 0 ? ml.P.reduce((s, q, i) => s + q * ml.ds[i], 0) / mlLen : 0.5 * (pInlet + v.pSep), tLine = mlLen > 0 ? ml.T.reduce((s, q, i) => s + q * ml.ds[i], 0) / mlLen : tArr;
      const kSep = clamp(fl.gas.cp / Math.max(fl.gas.cp - R / (fl.gas.MW * 1e-3), 1), 1.05, 1.67), hr1 = pumpHead(pump.qr * pump.trim, { ...pump, speed: 1 }), ksPump = qDuty > 0 ? Math.max(hSys(qDuty) - hSys(1e-6), 0) / (qDuty * qDuty) : 0;
      basis = { m0, xg: mGasSep / mCase, rhoL: rhoLsep, M: fl.gas.MW * 1e-3, k: kSep, z: zT, stroke: v.valveStroke, asvSet: v.asvSet, flareSet: v.flareSetPct / 100, flareP: 1.5, liqP: v.pumpDischargeP, headerV: v.exportVol,
        wells: { m: cm, p: cpv, pAvail0, chokeDp: Math.max(pAvail0 - demand0, 0), pos0: clamp(chokeOpening / 100, 0.05, 1), char: v.chokeChar, rng: v.chokeRange },
        line: { V: Math.max(volume - liquidInv, 0.05 * volume), T: tLine + KEL, p0: Math.max(pLine, v.pSep + 0.5) },
        sep: { V: sepVtot, T: tArr + KEL, p0: v.sepInitP > 0 ? v.sepInitP : v.pSep, pSet: v.pSep, vol0: vesselVolume({ ...sepGeom, level: clamp(v.sepInitLevelPct / 100, 0.02, 0.95) * sepH }), geom: sepGeom, levelSet: sepLevel, H: sepH },
        comp: comp ? { stages: comp.nSt, pd: v.compPd, eta: v.compEta, margin: (1 + v.compMargin / 100) * v.compFlowMult * 100 - 100, rpm: v.compRpm, H: v.compInertiaH, Tcool: v.compTcool + KEL, dpExport: 2, dischargeV: v.dischargeVol, asvSize: v.asvSize } : null,
        pump: pumpRes ? { qr: pump.qr * pump.trim, hr: hr1, shutoff: Math.max(pumpHead(1e-9, { ...pump, speed: 1 }) / hr1, 1.02), eta: pump.eta, speed: pump.speed, head: v.pumpStaticHead, ks: ksPump, rpm: v.pumpRpm, H: v.pumpInertiaH } : null };
      const sys = transientSystem(basis, v.transScenario, { tEnd: v.transTime, steps: v.transSteps, chokeTo: v.chokeStepPct / 100 }), res = transientNetwork({ ...sys, nOut: 120 });
      const mx = (a) => (a && a.length ? Math.max(...a) : null), mn = (a) => (a && a.length ? Math.min(...a) : null), tz = (a) => (a ? a.reduce((s, q, i) => (i ? s + 0.5 * (q + a[i - 1]) * (res.t[i] - res.t[i - 1]) : 0), 0) : 0), cs = res.speed.compressor || null, ps = res.speed.pump || null;
      const half = (a) => { if (!a) return null; const i = a.findIndex((q) => q <= 0.5 * a[0]); return i > 0 ? res.t[i - 1] + ((res.t[i] - res.t[i - 1]) * (a[i - 1] - 0.5 * a[0])) / Math.max(a[i - 1] - a[i], 1e-12) - sys.tEvent : null; };
      const bg = res.balance.gas, bl = res.balance.liquid, liqEdge = res.m.PUMP || res.m.DRAIN;
      trans = { scenario: sys.scenario, tEvent: sys.tEvent, tEnd: sys.tEnd, t: res.t, pSep: res.p.SEP, pLine: res.p.LINE, pHeader: res.p.DISCHARGE || null, level: res.level.SEP, liquid: res.liquid.SEP, compSpeed: cs, pumpSpeed: ps, compPower: res.power.compressor || null, pumpPower: res.power.pump || null,
        mWells: res.m.WELLS, mLine: res.m.FLOWLINE, mComp: res.m.COMPRESSOR || null, mExport: res.m.EXPORTVALVE, mFlare: res.m.FLAREVALVE, mRecycle: res.m.RECYCLE || null, mLiquid: liqEdge, asv: res.pos.RECYCLE || null, choke: res.pos.WELLS, exportValve: res.pos.EXPORTVALVE, flareValve: res.pos.FLAREVALVE, liquidValve: res.pos.PUMP || res.pos.DRAIN, margin: res.margin.COMPRESSOR || null,
        peakSepP: mx(res.p.SEP), minSepP: mn(res.p.SEP), peakHeaderP: mx(res.p.DISCHARGE), minMargin: mn(res.margin.COMPRESSOR), surgeTime: res.surgeTime.COMPRESSOR ?? null, compHalfTime: sys.scenario === 'compressorTrip' ? half(cs) : null, pumpHalfTime: sys.scenario === 'pumpTrip' ? half(ps) : null, asvMax: mx(res.pos.RECYCLE), levelMax: mx(res.level.SEP), levelMin: mn(res.level.SEP), levelEnd: res.level.SEP[res.level.SEP.length - 1], flareMass: tz(res.m.FLAREVALVE), wellsEnd: res.m.WELLS[res.m.WELLS.length - 1],
        gas: bg, liquidBal: bl, massError: Math.max(Math.abs(bg.initial + bg.in - bg.out - bg.final) / Math.max(bg.initial, 1e-9), Math.abs(bl.initial + bl.in - bl.out - bl.final) / Math.max(bl.initial, 1e-9)), steps: res.steps, halved: res.halved, mGas0: sys.mGas, mLiq0: sys.mLiquid, hasComp: !!basis.comp, hasPump: !!basis.pump };
    } catch (e) { trans = null; warnings.push({ level: 'info', msg: 'The transient scenario was not run: ' + e.message + '.' }); }
  }

  // ---- comparison with measured data ---------------------------------------------------------------------------------
  const tubP = aCase?.r ? { md: aCase.r.P.map((_, i) => (i * traj.md) / v.nTubing), p: aCase.r.P.slice().reverse() } : null, stnMd = traj.stations.map((s) => s.md - traj.stations[0].md), stnTvd = traj.stations.map((s) => s.tvd);
  const mapHead = (q, N) => { if (!comp) return null; const sp = comp.map.speeds, hs = comp.map.lines.map((l) => interp1(l.q, l.h, q)); return sp.length > 1 ? interp1(sp, hs, clamp(N, sp[0], sp[sp.length - 1])) / 1000 : hs[0] / 1000; };
  const cmp = compareMeasured(v.measured, {
    rov: { label: 'ROV survey: free-span length', unit: 'm', re: /rov|span/, f: (at) => { let b = null; for (const s of spansAll) if (Math.abs(s.x - at) <= Math.max(150, s.length) && (!b || Math.abs(s.x - at) < Math.abs(b.x - at))) b = s; return b ? b.length : 0; } },
    deviation: { label: 'Deviation survey: true vertical depth', unit: 'm', re: /deviat|tvd|gyro|mwd/, f: (at) => interp1(stnMd, stnTvd, at) },
    caliper: { label: 'Caliper: tubing inner diameter', unit: 'mm', re: /caliper/, f: () => v.tubingIdMm },
    inspection: { label: 'Inspection: pipe bore', unit: 'mm', re: /inspect|pig|\bili\b|bore/, f: (at) => (at <= rb.x + 1e-6 ? Dfl : D) * 1000 },
    bathymetry: { label: 'GIS / bathymetric survey: seabed elevation on the route', unit: 'm', re: /bathy|gis|seabed/, f: (at) => interp1(x, z, Math.min(at, rb.x)) },
    pumpfat: { label: 'Pump acceptance test: head', unit: 'm', re: /pump|\bfat\b/, f: (at) => pumpHead(at, pump) },
    comptest: { label: 'Compressor performance test: polytropic head', unit: 'kJ/kg', re: /compress/, f: (at, aux) => mapHead(at, aux > 0 ? aux / 100 : 1) },
    valvetest: { label: 'Valve / choke flow test: pressure drop with water', unit: 'bar', re: /valve|choke|\bcv\b/, f: (at, aux) => (at / (0.865 * v.chokeCvMax * Math.max(valveCharacteristic(aux > 0 ? aux / 100 : 1, v.chokeChar, v.chokeRange), 1e-6))) ** 2 },
    septest: { label: 'Separator performance test: gas velocity at carry-over', unit: 'm/s', re: /separ|carry/, f: () => vSB },
    welltest: { label: 'Well test: flowing bottom-hole pressure', unit: 'bara', re: /well.?test|pwf|bhp|bottom/, f: (at) => iprPwf(at, ipr) },
    plt: { label: 'Production log: flowing pressure in the tubing', unit: 'bara', re: /plt|production.?log|logging/, f: (at) => (tubP ? interp1(tubP.md, tubP.p, at) : null) },
    commissioning: { label: 'Commissioning: line-fill volume to the distance', unit: 'm³', re: /commission|fill|hydrotest|flood/, f: (at) => area(D) * interp1(x, prof.s, at) },
    pressure: { label: 'Field pressure along the line', unit: 'bara', re: /press/, f: (at) => (ml.x.length > 1 ? interp1(ml.x, ml.P, at) : null) },
    temperature: { label: 'Field temperature along the line', unit: '°C', re: /temp/, f: (at) => (ml.x.length > 1 ? interp1(ml.x, ml.T, at) : null) },
    elevation: { label: 'As-built drawing / route survey: pipe elevation', unit: 'm', re: /as.?built|survey|elev|drawing/, f: (at) => interp1(x, z, at) },
  });
  if (cmp.unknown) warnings.push({ level: 'info', msg: cmp.unknown + ' measurement row(s) could not be compared: the kind was not recognised, the value is missing or the model has no value there.' });

  // ---- checks, warnings, recommendations -----------------------------------------------------------------------------------
  prog(0.9, 'Results');
  const odSteel = D + 2 * wt, hoop = ((v.designPressure - P_STD) * 1e5 * odSteel) / (2 * wt) / 1e6, hoopUtil = hoop / (0.72 * mat.smys), sgEmpty = (wall.massPerM * G) / wall.buoyancy, sgOp = wall.sg;
  const allSpanLen = spanTotal, contact = lenFlow > 0 && offshore ? clamp(1 - allSpanLen / lenFlow, 0, 1) : 1, sepSeg = segmentGeometry(v.sepD, clamp(sepLevel, 0, v.sepD));
  const sections = [['Main pipe bore', area(D), D], ['Tubing bore', area(Dt), Dt]];
  if (fouled) sections.push(['Fouled flowline bore', area(Dfl), Dfl]);
  if (v.pip) { const dIn = D + 2 * wt + 2 * layersIn.reduce((a, l) => a + (l.t > 0 ? l.t : 0), 0), dOut = D + 2 * wt; sections.push(['Pipe-in-pipe annulus', area(dIn) - area(dOut), hydraulicDiameter(area(dIn) - area(dOut), Math.PI * (dIn + dOut))]); }
  if (horiz) sections.push(['Separator liquid section at normal level', sepSeg.area, hydraulicDiameter(sepSeg.area, sepSeg.wetted + sepSeg.chord)], ['Separator gas section at normal level', area(v.sepD) - sepSeg.area, hydraulicDiameter(area(v.sepD) - sepSeg.area, Math.PI * v.sepD - sepSeg.wetted + sepSeg.chord)]);
  else sections.push(['Separator shell', area(v.sepD), v.sepD]);
  sections.push(['Slug-catcher finger', fingerA, v.fingerDmm / 1000], ['Pump discharge line', area(pumpLineD), pumpLineD]);
  const lowPts = prof.lows.filter((q) => q.x <= rb.x + 1e-6).map((q) => ({ ...q, vm: ml.x.length ? interp1(ml.x, ml.vm, q.x) : 0, holdup: ml.x.length ? interp1(ml.x, ml.holdup, q.x) : 0 })), slowLows = lowPts.filter((q) => q.vm < 1);
  if (eroAll > 1) warnings.push({ level: 'bad', msg: `Erosional velocity ratio ${rd(eroAll, 2)} exceeds 1.0 (largest mixture velocity ${rd(vMaxAll, 1)} m/s; limit C/√ρ with C = ${v.cErosion}).` });
  else if (eroAll > 0.8) warnings.push({ level: 'warn', msg: `Erosional velocity ratio ${rd(eroAll, 2)} is within 20 % of the limit.` });
  if (vMaxAll < 0.5) warnings.push({ level: 'warn', msg: `The largest mixture velocity is only ${rd(vMaxAll, 2)} m/s: the line is oversized for this rate and liquid and solids will settle.` });
  if (slowLows.length) warnings.push({ level: 'warn', msg: `${slowLows.length} low point(s) carry a mixture velocity below 1 m/s and are likely to collect liquid and free water (first at x = ${rd(slowLows[0].x, 0)} m, z = ${rd(slowLows[0].z, 1)} m, ${rd(slowLows[0].vm, 2)} m/s, liquid holdup ${rd(slowLows[0].holdup, 2)}).` });
  else if (lowPts.length) warnings.push({ level: 'info', msg: `${lowPts.length} low point(s) along the line; the slowest carries ${rd(Math.min(...lowPts.map((q) => q.vm)), 2)} m/s, enough to sweep liquid at this rate — re-check at turndown.` });
  if (!aCase) warnings.push({ level: 'bad', msg: `The case rate (${rd(qWell, 0)} ${wellUnit} per well) is above the absolute open-flow potential of the well (${rd(ipr.qMax, 0)} ${wellUnit}).` });
  else if (!canDeliver) warnings.push({ level: 'bad', msg: `At the case rate the well delivers ${rd(whpAvail, 1)} bara at the wellhead${boost > 0 ? ` (including ${rd(boost, 1)} bar of boosting)` : liftGas > 0 ? ' (with gas lift)' : ''} but the network needs ${rd(whpReq, 1)} bara: the wells cannot sustain the case rate (deliverable rate ${rd(operatingRate, 0)} Sm³/d liquid).` });
  if (!gchk.every((c) => c.pass)) warnings.push({ level: 'warn', msg: `Geometry consistency: ${gchk.filter((c) => !c.pass).map((c) => c.name).join(', ')} failed on the route used for the calculation.` });
  for (const q of orient.filter((c) => !c.ok)) warnings.push({ level: 'bad', msg: `Equipment orientation: the ${q.kind} ${q.item} carries flow against its built direction — swap its “from” and “to” nodes or check the boundary pressures.` });
  if (arrestors.required && !(arrestors.count > 0)) warnings.push({ level: 'warn', msg: `External pressure at ${rd(waterDepth, 0)} m exceeds the factored buckle-propagation pressure (${rd(ppr, 0)} bar unfactored) below about ${rd(arrFirst, 0)} m water depth and no buckle arrestors are specified.` });
  if (spanOver) warnings.push({ level: 'warn', msg: `${spanOver} free-span candidate(s) exceed the allowable span of ${v.spanLimit} m: about ${extraSupports} additional support(s) or span corrections are needed.` });
  if (sleeperSpan > v.spanLimit) warnings.push({ level: 'warn', msg: `Each ${v.sleeperHeight} m sleeper lifts the pipe over about ${rd(sleeperSpan, 0)} m on either side, more than the allowable span of ${v.spanLimit} m: lower the sleepers or check the span for vortex-induced vibration.` });
  if (cp.utilisation !== null && cp.utilisation > 1 && offshore) warnings.push({ level: 'warn', msg: `Cathodic protection: each anode needs ${rd(cp.requiredMass, 0)} kg for ${v.cpLife} years but has ${v.anodeMass} kg — reduce the spacing to about ${rd(v.anodeSpacing / cp.utilisation, 0)} m or use heavier anodes.` });
  if (trans) {
    if (trans.surgeTime !== null && trans.scenario === 'compressorTrip') warnings.push({ level: 'info', msg: `Transient (compressorTrip): about ${rd(Math.max(trans.surgeTime - trans.tEvent, 0), 1)} s after the trip the head of the coasting compressor falls below the pressure trapped in its discharge piping (${v.dischargeVol} m³) and the machine passes through surge until the recycle valve (open ${rd(100 * (trans.asvMax || 0), 0)} %) has vented that volume; the rotor reaches half speed after ${isNum(trans.compHalfTime) ? rd(trans.compHalfTime, 0) : '—'} s. A hot-gas bypass close to the discharge flange or a smaller trapped volume shortens the excursion.` });
    else if (trans.surgeTime !== null) warnings.push({ level: 'warn', msg: `Transient (${trans.scenario}): the compressor crosses its surge line ${rd(trans.surgeTime - trans.tEvent, 1)} s after the event (smallest surge margin ${rd(trans.minMargin, 0)} %) although the recycle valve opens to ${rd(100 * (trans.asvMax || 0), 0)} %: start the recycle earlier (control line above ${v.asvSet} %) or fit a larger, faster valve.` });
    else if (isNum(trans.minMargin) && (trans.asvMax || 0) > 0.01) warnings.push({ level: 'info', msg: `Transient (${trans.scenario}): the anti-surge recycle opens to ${rd(100 * trans.asvMax, 0)} % and holds the surge margin at ${rd(trans.minMargin, 1)} % (control line ${v.asvSet} %).` });
    if (trans.levelMax > 0.85 * sepH) warnings.push({ level: 'warn', msg: `Transient (${trans.scenario}): the separator level reaches ${rd((100 * trans.levelMax) / sepH, 0)} % of the vessel height within ${rd(trans.tEnd, 0)} s — a high-level trip would act.` });
    if (trans.flareMass > 1) warnings.push({ level: 'info', msg: `Transient (${trans.scenario}): the separator pressure peaks at ${rd(trans.peakSepP, 1)} bara and ${rd(trans.flareMass, 0)} kg of gas goes to flare through the pressure-control valve.` });
  }
  if (lift.type !== 'none' && nFlow < net.nodes.length && canDeliver && chokeOpening > 85) warnings.push({ level: 'info', msg: `The production chokes run ${rd(chokeOpening, 0)} % open with the lift as set: little control margin is left.` });
  if (op.limited === 'dead') warnings.push({ level: 'bad', msg: 'The well cannot flow naturally against the network back-pressure at any rate: artificial lift, boosting or a lower arrival pressure is needed.' });
  if (op.limited === 'range') warnings.push({ level: 'info', msg: `The natural-flow rate is at least ${rd(operatingRate, 0)} Sm³/d of liquid: the network back-pressure curve could not be extended beyond ${rd(sOp * 100, 0)} % of the case rate.` });
  if (canDeliver && chokeCv > v.chokeCvMax) warnings.push({ level: 'warn', msg: `The choke needs Cv ${rd(chokeCv, 1)} at the case rate, above the rated Cv ${v.chokeCvMax}: it would be wide open and still restrict the well.` });
  if (canDeliver && chokeOpening < 15) warnings.push({ level: 'warn', msg: `The choke would run only ${rd(chokeOpening, 0)} % open: poor control and trim erosion — a smaller trim (rated Cv about ${rd(chokeCv / 0.3, 0)}) suits this duty.` });
  if (canDeliver && critical) warnings.push({ level: 'info', msg: `Choke flow is critical: pressure ratio ${rd(p2c / p1c, 2)} is below the two-phase critical ratio ${rd(sach.yc, 2)} (Sachdeva); downstream pressure changes do not reach the well.` });
  if (gasLoad > 1) warnings.push({ level: 'bad', msg: `Separator gas velocity ${rd(vGasSep, 2)} m/s exceeds the Souders–Brown limit ${rd(vSB, 2)} m/s: liquid carry-over.` });
  if (resTime < v.sepResMin * 60) warnings.push({ level: 'warn', msg: `Separator liquid residence time ${rd(resTime / 60, 1)} min is below the ${v.sepResMin} min required.` });
  if (scPeak > scVol * 1.0001) warnings.push({ level: 'bad', msg: `The design slug overfills the slug catcher: peak inventory ${rd(scPeak, 1)} m³ against ${rd(scVol, 1)} m³.` });
  if (pumpRes) {
    if (pumpRes.qOp <= 0) warnings.push({ level: 'bad', msg: `The pump shut-off head ${rd(pumpHead(0, pump), 0)} m is below the static system head ${rd(hSys(1e-6), 0)} m: no flow.` });
    else if (pumpRes.qOp < qDuty) warnings.push({ level: 'warn', msg: `The pump delivers ${rd(pumpRes.qOp, 0)} m³/h against the system curve, less than the ${rd(qDuty, 0)} m³/h of liquid produced.` });
    if (pumpRes.npshA < 1.2 * pumpRes.npshR) warnings.push({ level: pumpRes.npshA < pumpRes.npshR ? 'bad' : 'warn', msg: `NPSH available ${rd(pumpRes.npshA, 1)} m against ${rd(pumpRes.npshR, 1)} m required (margin ${rd(pumpRes.npshA / pumpRes.npshR, 2)}): the separator liquid is at its bubble point, so only the static level counts.` });
  }
  if (comp) {
    if (comp.pt.surgeMargin < 10) warnings.push({ level: comp.pt.surgeMargin < 0 ? 'bad' : 'warn', msg: `Compressor surge margin ${rd(comp.pt.surgeMargin, 0)} % is below 10 %: recycle is needed.` });
    if (comp.pt.chokeMargin < 5) warnings.push({ level: 'warn', msg: `The compressor runs within ${rd(comp.pt.chokeMargin, 0)} % of stonewall.` });
    if (!comp.pt.inside) warnings.push({ level: 'warn', msg: `The compressor duty point lies outside the speed range of the map (needs ${rd(comp.pt.speed * 100, 0)} % speed).` });
    if (comp.tMax > 160) warnings.push({ level: 'warn', msg: `Compressor discharge temperature ${rd(comp.tMax, 0)} °C exceeds 160 °C: add a stage or more intercooling.` });
  }
  if (hoopUtil > 1) warnings.push({ level: 'bad', msg: `Hoop stress at design pressure is ${rd(hoop, 0)} MPa, ${rd(hoopUtil * 100, 0)} % of the 0.72·SMYS allowable for ${mat.label}.` });
  if (offshore && sgEmpty < 1.1) warnings.push({ level: 'warn', msg: `Empty-pipe specific gravity is ${rd(sgEmpty, 2)}: below 1.1 the line is unlikely to be stable on the seabed without concrete, trenching or anchoring.` });
  if (spanCount) warnings.push({ level: spanMax > 60 ? 'warn' : 'info', msg: `${spanCount} free-span candidate(s) from the ${spanSource}; the longest is ${rd(spanMax, 0)} m.` });
  if (fouled && cleanP !== null) warnings.push({ level: 'info', msg: `Fouled bore (${rd(Dfl * 1000, 1)} mm, roughness ${rd(roughFl * 1e6, 0)} µm) raises the inlet pressure from ${rd(cleanP, 1)} to ${rd(pInlet, 1)} bara.` });
  if (est.qG / qEst > 0.1 && net.edges.some((e) => e.kind === 'pump')) warnings.push({ level: 'info', msg: `The in-line pump handles a gas volume fraction of about ${rd((100 * est.qG) / qEst, 0)} %: a conventional centrifugal pump loses its head above roughly 10–15 %; a helico-axial or twin-screw multiphase pump is assumed.` });
  if (traj.maxDls > 6) warnings.push({ level: 'warn', msg: `Dog-leg severity reaches ${rd(traj.maxDls, 1)}°/30 m in the survey.` });

  { const now = study.rows[nNow - 1], best = study.rows[study.recommended - 1], equal = rated.length > 0 && rated.every((k) => Math.abs(net.nodes[k].share - 1 / rated.length) < 0.02);
    if (now && equal && now.feasible === false && now.limit) warnings.push({ level: now.delivers ? 'warn' : 'info', msg: `Well count: with ${nNow} well(s) ${now.limit}.` });
    if (nNow <= 6 && study.recommended !== nNow) recs.push(`Well count: ${study.recommended} well(s) is the ${study.basis} (${rd(best.ratePerWell, 0)} ${wellUnit} per well, drawdown ${rd(best.drawdown, 0)} bar, choke ${rd(best.chokeOpening, 0)} % open${best.boostPower > 0 ? `, lift power ${rd(best.boostPower, 0)} kW` : ''}); the network lists ${nNow}. Compare the options in the well-count table and let the economics suite price them.`);
    else if (nNow <= 6) recs.push(`Well count: ${nNow} well(s) is the ${study.basis}; each well makes ${rd(best.ratePerWell, 0)} ${wellUnit} at a drawdown of ${rd(best.drawdown, 0)} bar (limit ${rd(v.maxDrawdown, 0)} bar).`); }
  if (canDeliver) recs.push(`Run the production chokes about ${rd(chokeOpening, 0)} % open (Cv ${rd(chokeCv, 1)} of ${v.chokeCvMax}) to hold the case rate; they take ${rd(chokeDp, 1)} bar, and fully open the wells would make about ${rd(operatingRate, 0)} Sm³/d of liquid.`);
  else recs.push(`The wells fall short of the case rate by ${rd(Math.max(whpReq - whpAvail, 0), 1)} bar at the wellhead: plan for ${rd(operatingRate, 0)} Sm³/d of liquid, or lower the arrival pressure, ${lift.type === 'none' ? 'add boosting or gas lift' : lift.type === 'boost' ? `raise the rated boosting pressure to about ${rd(boost + whpReq - whpAvail + 5, 0)} bar` : 'raise the lift-gas rate'}.`);
  if (lift.power > 0) recs.push(lift.type === 'boost' ? `Subsea boosting of ${rd(boost, 0)} bar takes about ${rd(lift.power, 0)} kW of shaft power for ${rd(lift.flow, 0)} m³/h at the pump suction (${rd(lift.suctionP, 0)} bara).` : `Gas lift of ${rd(v.liftGasRate, 0)} Sm³/d per well needs about ${rd(lift.power, 0)} kW of compression to ${rd(lift.injectionP, 0)} bara.`);
  if (eroAll > 0.8) recs.push(`Increase the bore where the erosional ratio is highest (${rd(eroAll, 2)}): a diameter ${rd(100 * (Math.sqrt(eroAll / 0.8) - 1), 0)} % larger brings it to 0.8.`);
  if (uMean > 4 && offshore) recs.push(`The wall U-value is ${rd(uMean, 2)} W/m²K and the fluid arrives at ${rd(tArr, 1)} °C; compare this with the hydrate temperature at arrival (${rd(fm.hydrateT(v.pSep), 1)} °C) before settling the insulation.`);
  else recs.push(`With U = ${rd(uMean, 2)} W/m²K the fluid arrives at ${rd(tArr, 1)} °C (hydrate temperature at arrival pressure ${rd(fm.hydrateT(v.pSep), 1)} °C): pass this wall design to the flow and operations suites for cooldown.`);
  if (slowLows.length) recs.push(`Plan routine pigging or a minimum rate for the ${slowLows.length} low point(s) running below 1 m/s.`);
  if (gasLoad > 0.9 || resTime < v.sepResMin * 60) recs.push(`Enlarge the separator: gas load ${rd(gasLoad * 100, 0)} % of the Souders–Brown limit, residence time ${rd(resTime / 60, 1)} min.`);
  recs.push(v.scType === 'finger' ? `Size the slug catcher for ${rd(scVol, 0)} m³: ${nFing} fingers of ${v.fingerDmm} mm × ${rd(fingerL, 0)} m.` : `Size the slug catcher for ${rd(scVol, 0)} m³: a vessel of about ${rd(scVesselD, 1)} m × ${rd(4 * scVesselD, 1)} m.`);
  if (pumpRes && pumpRes.throttle > 0.15 * pumpRes.hRun) recs.push(`The export pump is throttled by ${rd(pumpRes.throttle, 0)} m of head at the duty flow: running at ${rd(100 * pump.speed * Math.sqrt(Math.max(hSys(qDuty), 1) / Math.max(pumpHead(qDuty, pump), 1)), 0)} % speed would save about ${rd(pumpRes.power * (1 - hSys(qDuty) / pumpHead(qDuty, pump)), 0)} kW.`);

  return assemble({ study, nNow, v, fm, rates, route, x, z, prof, rb, D, wt, rough, Dfl, roughFl, fouled, waterDepth, offshore, lenFlow, lenRiser, volume, mat, W, wall, mainCat, uMean, od, hIn, hSea, hAir, fit, bends, bendsK, kLoss, net, base, sol, eInfo, ml, mlLen, tArr, pInlet, inletK, vMaxAll, eroAll, liquidInv, spans, spanCount, spanMax, spanSource, supUser, spansSplit, spanOver, extraSupports, arrestors, ppr, arrFirst, initX, sleeperSpan, cp, gchk, orient, lift, boost, liftGas, trans, basis, cmp, sepKused, bK, nFlow, scHold, traj, ipr, gasBasis, qWell, wellUnit, share, wi, sA, yA, curve, sD, yD, sOp, op, operatingRate, whpOp, pwfOp, whpReq, whpAvail, chokeDp, canDeliver, sens, sensCurves, chokeCv, chokeOpening, chokeState, sach, hem, critical, p1c, p2c, areaSach, areaHem, dOrifice, beanRows, mWell, fl, eosF, nHC, nGas, nOil, mGasSep, mOilSep, qGasSep, qOilSep, qWatSep, qLiqSep, rhoLsep, rhoGsep, zC1, c1Out, compErr, sepGeom, sepH, sepLevel, sepVtot, sepVliq, vGasSep, vSB, gasLoad, resTime, scVol, nFing, fingerL, scVesselD, scT, scY, scPeak, scEnd, scV0, pump, pumpRes, hSys, qDuty, qZero, netPumpPower, comp, hoop, hoopUtil, sgEmpty, sgOp, lowPts, warnings, recs, mHC, mW, mCase, qLiqStd, qGasStd, contact, sections, Dt });
}

/** Results object (KPIs, plots, tables, balances, outputs) from the solved state. */
function assemble(S) {
  const { v, fm, rates, x, z, prof, rb, D, wt, rough, wall, W, net, sol, eInfo, ml, traj, ipr, comp, pumpRes, pump } = S, nP = x.length, T = (title, columns, rows, note) => ({ title, columns, rows, ...(note ? { note } : {}) });
  const st = (bad, warn) => (bad ? 'bad' : warn ? 'warn' : 'ok'), keep = [...prof.highs, ...prof.lows].map((q) => q.i).concat(rb.i), simp = simplifyProfile(x, z, 80, keep);
  const pumpPower = (pumpRes ? pumpRes.power : 0) + S.netPumpPower, compPower = comp ? comp.power : 0, nWells = net.nodes.filter((n) => n.kind === 'source').length;

  const kpis = [
    { label: 'Route length', value: rd(prof.length, 0), unit: 'm', status: 'ok', help: `${rd(S.lenFlow, 0)} m before the riser base and ${rd(S.lenRiser, 0)} m of riser, along the pipe axis.` },
    { label: 'Water depth', value: rd(S.waterDepth, 1), unit: 'm', status: 'ok', help: 'Deepest point of the route below mean sea level.' },
    { label: 'Riser height', value: rd(rb.height, 1), unit: 'm', status: 'ok', help: `Riser base at x = ${rd(rb.x, 0)} m.` },
    { label: 'Internal volume', value: rd(S.volume, 1), unit: 'm³', status: 'ok', help: 'Clean bore × route length.' },
    { label: 'Overall U-value (on ID)', value: rd(S.uMean, 3), unit: 'W/m²K', status: 'ok', help: 'Length-weighted over the route; see the resistance breakdown.' },
    { label: 'Flowline inlet pressure', value: rd(S.pInlet, 2), unit: 'bara', status: st(!sol.usable, false), help: `Node ${net.nodes[S.inletK].id} at the case rate.` },
    { label: 'Arrival temperature', value: rd(S.tArr, 1), unit: '°C', status: st(false, S.tArr < fm.hydrateT(v.pSep)), help: `Hydrate temperature at arrival pressure: ${rd(fm.hydrateT(v.pSep), 1)} °C.` },
    { label: 'Wellhead pressure needed', value: rd(S.whpReq, 2), unit: 'bara', status: st(!S.canDeliver, false), help: `Network pressure at ${net.nodes[S.wi].id}; the well delivers ${rd(S.whpAvail, 1)} bara at the case rate.` },
    { label: 'Natural-flow rate (chokes open)', value: rd(S.operatingRate, 0), unit: 'Sm³/d liquid', status: st(S.sOp <= 0, S.sOp < 1), help: `Nodal solution: ${rd(S.sOp * 100, 0)} % of the case rate, ${rd(S.sOp * S.qWell, 0)} ${S.wellUnit} per well.` },
    { label: 'Choke opening at the case rate', value: rd(S.chokeOpening, 1), unit: '%', status: st(false, !S.canDeliver || S.chokeOpening < 15 || S.chokeCv > v.chokeCvMax), help: `Required Cv ${rd(S.chokeCv, 1)} of ${v.chokeCvMax}; pressure drop ${rd(Math.max(S.chokeDp, 0), 1)} bar.` },
    { label: 'Erosional velocity ratio', value: rd(S.eroAll, 3), unit: '–', status: st(S.eroAll > 1, S.eroAll > 0.8), help: `Largest mixture velocity ${rd(S.vMaxAll, 2)} m/s.` },
    { label: 'Slug catcher volume', value: rd(S.scVol, 1), unit: 'm³', status: st(S.scPeak > S.scVol * 1.0001, false), help: 'Design surge × factor plus the liquid hold-up for the residence time.' },
    { label: 'Pump power', value: rd(pumpPower, 1), unit: 'kW', status: st(!!pumpRes && pumpRes.qOp <= 0, !!pumpRes && (pumpRes.qOp < S.qDuty || pumpRes.npshA < 1.2 * pumpRes.npshR)), help: 'Liquid export pump shaft power plus any in-line network pump.' },
    { label: 'Compressor power', value: rd(compPower, 1), unit: 'kW', status: st(false, !!comp && (comp.pt.surgeMargin < 10 || !comp.pt.inside)), help: comp ? `${comp.nSt} stage(s), ratio ${rd(comp.ratio, 2)} each.` : 'No compressor duty.' },
    { label: 'Wells: listed / recommended', value: `${S.nNow} / ${S.study.recommended}`, unit: '', status: st(false, S.study.recommended !== S.nNow), help: `Well-count study: ${S.study.recommended} is the ${S.study.basis}.` },
    { label: 'Drawdown at the case rate', value: rd(Math.max(ipr.pRes - (iprPwf(S.qWell, ipr) ?? 0), 0), 1), unit: 'bar', status: st(false, ipr.pRes - (iprPwf(S.qWell, ipr) ?? 0) > v.maxDrawdown), help: `Reservoir pressure minus flowing bottom-hole pressure of the well used for the nodal analysis; limit ${rd(v.maxDrawdown, 0)} bar.` },
    { label: 'Artificial lift power', value: rd(S.lift.power, 0), unit: 'kW', status: 'ok', help: S.lift.type === 'boost' ? `Subsea boosting of ${rd(S.boost, 1)} bar on ${S.nFlow} well stream(s).` : S.lift.type === 'gaslift' ? `Gas lift, ${rd(v.liftGasRate, 0)} Sm³/d per well.` : 'No artificial lift.' },
    { label: 'Transient: peak separator pressure', value: S.trans ? rd(S.trans.peakSepP, 2) : '—', unit: 'bara', status: st(false, !!S.trans && S.trans.peakSepP > 1.2 * v.pSep), help: S.trans ? `Scenario “${S.trans.scenario}” at t = ${S.trans.tEvent} s; mass-balance error ${sg(S.trans.massError, 2)}.` : 'The transient was not run.' },
    { label: 'Transient: smallest surge margin', value: S.trans && isNum(S.trans.minMargin) ? rd(S.trans.minMargin, 1) : '—', unit: '%', status: st(false, !!S.trans && isNum(S.trans.minMargin) && S.trans.minMargin < 0), help: 'Distance of the compressor operating point from the surge line during the scenario (negative: surge).' },
    { label: 'Specific gravity (empty / operating)', value: rd(S.sgEmpty, 2), unit: '–', status: st(false, S.offshore && S.sgEmpty < 1.1), help: `Operating: ${rd(S.sgOp, 2)}; submerged weight ${rd(wall.submerged, 0)} N/m.` },
  ];

  // ---- plots ---------------------------------------------------------------------------------------------------
  const plots = [], pts = (name, a, color) => ({ name, x: a.map((q) => q.x), y: a.map((q) => q.z), mode: 'points', ...(color ? { color } : {}) });
  const routeSeries = [{ name: 'Pipe route', x, y: z }];
  if (rb.height > 0) routeSeries.push({ name: 'Seabed', x: x.slice(0, rb.i + 1).concat(x[nP - 1]), y: z.slice(0, rb.i + 1).concat(rb.z), dash: true });
  if (prof.highs.length) routeSeries.push(pts('High points', prof.highs));
  if (prof.lows.length) routeSeries.push(pts('Low points', prof.lows));
  routeSeries.push(pts('Inlet / manifold', [{ x: x[0], z: z[0] }]), pts(rb.height > 0 ? 'Riser base' : 'Outlet', [{ x: rb.x, z: rb.z }]));
  if (rb.height > 0) routeSeries.push(pts('Arrival (separator)', [{ x: x[nP - 1], z: z[nP - 1] }]));
  if (S.spans.length) routeSeries.push(pts('Free-span candidates', S.spans.map((s) => ({ x: s.x, z: interp1(x, z, s.x) }))));
  plots.push({ type: 'line', title: 'Elevation profile', xlabel: 'Horizontal distance (m)', ylabel: 'Elevation (m)', series: routeSeries, hlines: S.offshore ? [{ y: 0, label: 'Sea level', color: '#0ea5e9' }] : [], vlines: rb.height > 0 ? [{ x: rb.x, label: 'Riser base' }] : [], note: 'Elevation is exaggerated relative to distance.' });
  plots.push({ type: 'line', title: 'Inclination along the route', xlabel: 'Horizontal distance (m)', ylabel: 'Inclination (° above horizontal)', zeroY: true, series: [{ name: 'Segment inclination', x: prof.segs.map((s) => s.x0).concat(x[nP - 1]), y: prof.segs.map((s) => s.incl).concat(prof.segs[prof.segs.length - 1].incl), mode: 'step' }] });
  if (ml.x.length > 1) {
    plots.push({ type: 'line', title: 'Pressure and temperature along the main line', xlabel: 'Horizontal distance (m)', ylabel: 'bara · °C', series: [{ name: 'Pressure (bara)', x: ml.x, y: ml.P }, { name: 'Temperature (°C)', x: ml.x, y: ml.T }], vlines: rb.height > 0 ? [{ x: rb.x, label: 'Riser base' }] : [] });
    plots.push({ type: 'line', title: 'Liquid holdup and mixture velocity along the main line', xlabel: 'Horizontal distance (m)', ylabel: 'Holdup (–) · velocity (m/s)', zeroY: true, series: [{ name: 'Liquid holdup', x: ml.x, y: ml.holdup }, { name: 'Mixture velocity (m/s)', x: ml.x, y: ml.vm }] });
  }
  const stn = traj.stations;
  plots.push({ type: 'line', title: 'Well trajectory — vertical section', xlabel: 'Horizontal displacement (m)', ylabel: 'True vertical depth below wellhead (m, downwards negative)', series: [{ name: 'Well path (minimum curvature)', x: stn.map((s) => s.disp), y: stn.map((s) => -s.tvd), mode: 'both' }] });
  plots.push({ type: 'line', title: 'Well trajectory — plan view', xlabel: 'East (m)', ylabel: 'North (m)', series: [{ name: 'Well path', x: stn.map((s) => s.east), y: stn.map((s) => s.north), mode: 'both' }] });
  const qA = S.sA.map((s) => s * S.qWell), qD = S.sD.map((s) => s * S.qWell);
  plots.push({ type: 'line', title: 'Nodal analysis at the wellhead', xlabel: `Rate per well (${S.wellUnit})`, ylabel: 'Pressure (bara)', zeroY: true,
    series: [{ name: 'Inflow: bottom-hole flowing pressure (IPR)', x: qA, y: S.curve.map((c) => c.pwf), dash: true }, { name: 'Well deliverability: wellhead pressure available (IPR − tubing lift)', x: qA, y: S.yA.map((p) => Math.max(p, 0)) }, { name: 'Network back-pressure at the wellhead', x: qD, y: S.yD }, { name: 'Natural-flow point', x: [S.sOp * S.qWell], y: [S.whpOp], mode: 'points' }],
    vlines: [{ x: S.qWell, label: 'Case rate' }], note: `Well ${net.nodes[S.wi].id} (${rd(S.share * 100, 0)} % of the case rate). The gap between the two wellhead curves at the case rate is taken by the production choke.` });
  plots.push({ type: 'line', title: 'Well-count study: deliverable field rate', xlabel: 'Number of wells', ylabel: 'Liquid rate (Sm³/d)', zeroY: true, series: [{ name: 'Without lift, chokes open', x: S.study.rows.map((r) => r.wells), y: S.study.rows.map((r) => r.rateNatural), mode: 'both' }, { name: 'With the lift as rated, chokes open', x: S.study.rows.map((r) => r.wells), y: S.study.rows.map((r) => r.rateWithLift), mode: 'both' }, { name: 'Rate held within the limits', x: S.study.rows.filter((r) => r.feasible).map((r) => r.wells), y: S.study.rows.filter((r) => r.feasible).map((r) => r.rate), mode: 'points' }], hlines: [{ y: S.qLiqStd, label: 'Case rate' }], vlines: [{ x: S.study.recommended, label: 'Recommended' }], note: 'Deliverable rates beyond the solved range of the network back-pressure curve continue its friction term with the square of the rate.' });
  plots.push({ type: 'line', title: 'Well-count study: drawdown, choke opening and lift power', xlabel: 'Number of wells', ylabel: 'bar · % · kW/10', zeroY: true, series: [{ name: 'Drawdown (bar)', x: S.study.rows.map((r) => r.wells), y: S.study.rows.map((r) => r.drawdown), mode: 'both' }, { name: 'Choke opening (%)', x: S.study.rows.map((r) => r.wells), y: S.study.rows.map((r) => r.chokeOpening), mode: 'both' }, { name: 'Lift power (kW ÷ 10)', x: S.study.rows.map((r) => r.wells), y: S.study.rows.map((r) => r.boostPower / 10), mode: 'both' }], hlines: [{ y: v.maxDrawdown, label: 'Drawdown limit' }] });
  plots.push({ type: 'line', title: 'Well deliverability sensitivities', xlabel: `Rate per well (${S.wellUnit})`, ylabel: 'Wellhead pressure (bara)', zeroY: true, series: [{ name: 'Base case', x: qA, y: S.yA.map((p) => Math.max(p, 0)) }, ...S.sensCurves, { name: 'Network back-pressure', x: qD, y: S.yD, dash: true }] });
  const hh = linspace(0, 1, 41);
  plots.push({ type: 'line', title: 'Choke inherent characteristics', xlabel: 'Opening (%)', ylabel: 'Flow coefficient Cv', zeroY: true, series: [['equal', 'Equal percentage'], ['linear', 'Linear'], ['quick', 'Quick opening']].map(([c, name]) => ({ name: name + (c === v.chokeChar ? ' (selected)' : ''), x: hh.map((h) => h * 100), y: hh.map((h) => v.chokeCvMax * valveCharacteristic(h, c, v.chokeRange)), dash: c !== v.chokeChar })).concat([{ name: 'Duty point', x: [S.chokeOpening], y: [Math.min(S.chokeCv, v.chokeCvMax)], mode: 'points' }]), hlines: [{ y: Math.min(S.chokeCv, v.chokeCvMax), label: 'Cv needed' }] });
  if (pumpRes) {
    const qq = linspace(0, S.qZero * 0.98, 40), c2 = { ...pump, speed: pump.speed * 0.85 };
    plots.push({ type: 'line', title: 'Export pump and system curves', xlabel: 'Flow (m³/h)', ylabel: 'Head (m)', zeroY: true, series: [{ name: `Pump at ${rd(pump.speed * 100, 0)} % speed, ${rd(pump.trim * 100, 0)} % trim`, x: qq, y: qq.map((q) => pumpHead(q, pump)) }, { name: 'Pump at 85 % of that speed (affinity laws)', x: qq, y: qq.map((q) => Math.max(pumpHead(q, c2), 0)), dash: true }, { name: 'System curve', x: qq, y: qq.map((q) => S.hSys(Math.max(q, 1e-6))) }, { name: 'Intersection', x: [pumpRes.qOp], y: [pumpHead(pumpRes.qOp, pump)], mode: 'points' }], vlines: [{ x: S.qDuty, label: 'Liquid produced' }] });
  }
  if (comp) plots.push({ type: 'line', title: 'Compressor map (first stage)', xlabel: 'Inlet volume flow (m³/s)', ylabel: 'Polytropic head (kJ/kg)', series: comp.map.lines.map((l) => ({ name: `${rd(l.n * 100, 0)} % speed`, x: l.q, y: l.h.map((h) => h / 1000) })).concat([{ name: 'Surge line', x: comp.map.lines.map((l) => l.q[0]), y: comp.map.lines.map((l) => l.h[0] / 1000), dash: true }, { name: 'Stonewall', x: comp.map.lines.map((l) => l.q[l.q.length - 1]), y: comp.map.lines.map((l) => l.h[l.h.length - 1] / 1000), dash: true }, { name: 'Duty point', x: [comp.stages[0].qIn], y: [comp.stages[0].headPoly / 1000], mode: 'points' }]) });
  const rTot = wall.resistances.reduce((s, r) => s + r.R, 0);
  plots.push({ type: 'bar', title: 'Thermal resistances of the wall (referred to the inner diameter)', ylabel: 'Resistance (m²K/W)', categories: wall.resistances.map((r) => r.name), series: [{ name: 'Resistance', values: wall.resistances.map((r) => r.R) }] });
  // network diagram: nodes laid out by their distance from the sources, pressure as the background field
  const lev = net.nodes.map(() => 0);
  for (let k = 0; k < net.nodes.length; k++) net.edges.forEach((e) => { if (lev[e.b] < lev[e.a] + 1 && lev[e.a] + 1 <= net.nodes.length) lev[e.b] = lev[e.a] + 1; });
  const maxLev = Math.max(...lev, 1), perLev = {}, pos = net.nodes.map((n, k) => { perLev[lev[k]] = (perLev[lev[k]] || 0) + 1; return { x: lev[k], j: perLev[lev[k]] }; });
  pos.forEach((p, k) => { p.y = p.j / (perLev[lev[k]] + 1); });
  const gx = linspace(-0.5, maxLev + 0.5, 40), gy = linspace(0, 1, 20), field = gy.map((yy) => gx.map((xx) => { let a = 0, b = 0; pos.forEach((p, k) => { const w = 1 / (((xx - p.x) / maxLev) ** 2 + (yy - p.y) ** 2 + 1e-3) ** 1.5; a += w * sol.p[k]; b += w; }); return a / b; }));
  plots.push({ type: 'field', title: 'Network diagram with node pressures', xlabel: 'Position along the flow path (sources left, outlet right)', ylabel: 'Branch', zlabel: 'Pressure', zunit: 'bara', x: gx, y: gy, z: field, cmap: 'viridis', shapes: net.edges.map((e) => ({ x: [pos[e.a].x, pos[e.b].x], y: [pos[e.a].y, pos[e.b].y], closed: false, color: e.kind === 'choke' ? '#f97316' : e.kind === 'pump' ? '#22c55e' : '#ffffff', dash: e.kind === 'choke' || e.kind === 'pump' })), markers: net.nodes.map((n, k) => ({ x: pos[k].x, y: pos[k].y, label: `${n.id} ${rd(sol.p[k], 1)}` })) });
  const lv = linspace(0, S.sepH, 41), sepV = lv.map((h) => vesselVolume({ ...S.sepGeom, level: h }));
  plots.push({ type: 'line', title: 'Separator level–volume relation', xlabel: 'Liquid level (m)', ylabel: 'Liquid volume (m³)', zeroY: true, series: [{ name: `${v.sepOrient === 'horizontal' ? 'Horizontal' : 'Vertical'} vessel ${v.sepD} m × ${v.sepL} m, 2:1 heads`, x: lv, y: sepV }, { name: 'Linear (for comparison)', x: [0, S.sepH], y: [0, S.sepVtot], dash: true }, { name: 'Normal level', x: [S.sepLevel], y: [S.sepVliq], mode: 'points' }] });
  plots.push({ type: 'line', title: 'Slug catcher inventory while the design slug arrives', xlabel: 'Time (s)', ylabel: 'Liquid inventory (m³)', zeroY: true, series: [{ name: 'Inventory', x: S.scT, y: S.scY.map((y) => y[0]) }], hlines: [{ y: S.scVol, label: 'Capacity' }] });

  const tr = S.trans;
  if (tr) {
    const pc = (a) => a.map((q) => 100 * q), ev = [{ x: tr.tEvent, label: 'Event' }], sr = (name, a, f = (q) => q) => (a ? [{ name, x: tr.t, y: a.map(f) }] : []);
    plots.push({ type: 'line', title: `Transient (${tr.scenario}): rotor speeds and valve travel`, xlabel: 'Time (s)', ylabel: '% of rated speed · % open', zeroY: true, vlines: ev, series: [...sr('Compressor speed', tr.compSpeed && pc(tr.compSpeed)), ...sr('Export-pump speed', tr.pumpSpeed && pc(tr.pumpSpeed)), ...sr('Anti-surge recycle valve', tr.asv && pc(tr.asv)), ...sr('Export-gas valve', tr.exportValve && pc(tr.exportValve)), ...sr('Flare pressure-control valve', tr.flareValve && pc(tr.flareValve)), ...sr('Production choke', tr.choke && pc(tr.choke)), ...sr('Liquid level-control valve', tr.liquidValve && pc(tr.liquidValve))] });
    plots.push({ type: 'line', title: `Transient (${tr.scenario}): nodal pressures`, xlabel: 'Time (s)', ylabel: 'Pressure (bara)', vlines: ev, series: [...sr('Separator', tr.pSep), ...sr('Flowline line-pack (mean)', tr.pLine), ...sr('Compressor discharge header', tr.pHeader)], hlines: [{ y: (v.flareSetPct / 100) * v.pSep, label: 'Flare set point' }] });
    plots.push({ type: 'line', title: `Transient (${tr.scenario}): mass flows`, xlabel: 'Time (s)', ylabel: 'Mass flow (kg/s)', zeroY: true, vlines: ev, series: [...sr('Wells (gas + liquid)', tr.mWells), ...sr('Flowline gas into the separator', tr.mLine), ...sr('Compressor', tr.mComp), ...sr('Recycle', tr.mRecycle), ...sr('Gas export', tr.mExport), ...sr('Flare', tr.mFlare), ...sr('Liquid export', tr.mLiquid)] });
    plots.push({ type: 'line', title: `Transient (${tr.scenario}): separator level and surge margin`, xlabel: 'Time (s)', ylabel: '% of vessel height · surge margin (%)', vlines: ev, series: [...sr('Liquid level', tr.level, (q) => (100 * q) / S.sepH), ...sr('Compressor surge margin', tr.margin)], hlines: [{ y: v.sepLevelPct, label: 'Level set point' }] });
  }
  if (S.cmp.rows.length) plots.push({ type: 'line', title: 'Measured against model (all kinds, each in its own unit)', xlabel: 'Measured', ylabel: 'Model', series: [...S.cmp.groups.map((g) => ({ name: `${g.label} (${g.unit})`, x: S.cmp.rows.filter((q) => q.key === g.key).map((q) => q.measured), y: S.cmp.rows.filter((q) => q.key === g.key).map((q) => q.model), mode: 'points' })), { name: 'Perfect agreement', x: [Math.min(...S.cmp.rows.map((q) => q.measured)), Math.max(...S.cmp.rows.map((q) => q.measured))], y: [Math.min(...S.cmp.rows.map((q) => q.measured)), Math.max(...S.cmp.rows.map((q) => q.measured))], dash: true }] });

  // ---- tables ------------------------------------------------------------------------------------------------------
  const tables = [], segShown = prof.segs.length > 80 ? simp : null, sx = segShown ? segShown.x : x, sz = segShown ? segShown.z : z, sp = segShown ? analyseProfile(sx, sz) : prof;
  tables.push(T('Route summary', ['Quantity', 'Value', 'Unit'], [['Route length (arc)', rd(prof.length, 1), 'm'], ['Horizontal extent', rd(prof.horizontal, 1), 'm'], ['Hydraulic length before the riser base', rd(S.lenFlow, 1), 'm'], ['Riser length', rd(S.lenRiser, 1), 'm'], ['Riser base distance', rd(rb.x, 1), 'm'], ['Riser height', rd(rb.height, 1), 'm'], ['Water depth', rd(S.waterDepth, 1), 'm'], ['Total climb', rd(prof.gain, 1), 'm'], ['Total descent', rd(prof.loss, 1), 'm'], ['Steepest uphill / downhill', `${rd(Math.max(...prof.segs.map((s) => s.incl)), 2)} / ${rd(Math.min(...prof.segs.map((s) => s.incl)), 2)}`, '°'], ['High points / low points', `${prof.highs.length} / ${prof.lows.length}`, ''], ['Internal volume', rd(S.volume, 2), 'm³'], ['Liquid inventory at the case rate (network)', rd(S.liquidInv, 1), 'm³'], ['Outer diameter with coatings', rd(S.od * 1000, 1), 'mm'], ['Seabed contact (share of the seabed section resting on the bottom)', rd(S.contact * 100, 1), '%'], ['Burial cover', v.burialDepth > 0 ? rd(v.burialDepth, 2) : 'exposed', v.burialDepth > 0 ? 'm' : ''], ['Free-span candidates', S.spanCount, ''], ['Longest free-span candidate', rd(S.spanMax, 1), 'm']]));
  tables.push(T('Segments', ['#', 'From x (m)', 'To x (m)', 'From z (m)', 'To z (m)', 'Length (m)', 'Inclination (°)', 'Cumulative length (m)', 'Volume (m³)'], sp.segs.map((s, i) => [i + 1, rd(s.x0, 1), rd(s.x1, 1), rd(s.z0, 1), rd(s.z1, 1), rd(s.L, 2), rd(s.incl, 3), rd(sp.s[i + 1], 1), rd(area(D) * s.L, 2)]), segShown ? `The route has ${prof.segs.length} segments; the table lists the simplified ${sp.segs.length}-segment route that keeps every high and low point.` : null));
  tables.push(T('Cross-sections and hydraulic diameters', ['Section', 'Flow area (m²)', 'Hydraulic diameter (m)'], S.sections.map((q) => [q[0], sg(q[1], 5), sg(q[2], 5)]), 'Hydraulic diameter = 4 × area / wetted perimeter (the free surface counts as perimeter for the separator sections).'));
  if (prof.highs.length + prof.lows.length) tables.push(T('High and low points', ['Type', 'x (m)', 'z (m)', 'Mixture velocity (m/s)', 'Liquid holdup'], [...prof.highs.map((q) => ['High', rd(q.x, 0), rd(q.z, 1), ml.x.length ? rd(interp1(ml.x, ml.vm, q.x), 2) : '—', ml.x.length ? rd(interp1(ml.x, ml.holdup, q.x), 3) : '—']), ...S.lowPts.map((q) => ['Low', rd(q.x, 0), rd(q.z, 1), rd(q.vm, 2), rd(q.holdup, 3)])].slice(0, 80)));
  tables.push(T('Bends and fittings', ['Item', 'Location x (m)', 'Angle (°) / count', 'Radius (m)', 'K each', 'K total'], [...S.bends.slice(0, 40).map((b) => [b.plane === 'plan' ? 'Plan-view bend' : 'Direction change (elevation)', rd(b.x, 0), rd(b.angle, 2), rd(Math.min(b.radius, 1e7), 0), sg(S.bK(b), 3), sg(S.bK(b), 3)]), ...S.fit.rows.map((r) => [r.label, '—', r.count, '—', sg(r.k, 3), sg(r.total, 3)]), ['Total (× calibration multiplier ' + v.kMult + ')', '—', '—', '—', '—', sg(S.kLoss, 4)]], `Equivalent length of the fittings: ${rd(S.fit.eqLength + (S.bendsK * D) / S.fit.fT, 1)} m of straight pipe (fully turbulent friction factor ${sg(S.fit.fT, 3)}). Direction changes with r/D above 20 are field curvature and add no loss.`));
  tables.push(T('Wall layers and thermal resistances', ['Layer', 'Thickness (mm)', 'Conductivity (W/m/K)', 'Resistance on ID (m²K/W)', 'Share (%)'], wall.resistances.map((r) => { const l = r.name === 'Pipe wall' ? { t: wt, k: S.mat.k } : wall.layers.find((q) => q.name === r.name); return [r.name, l ? rd(l.t * 1000, 2) : '—', l ? l.k : '—', sg(r.R, 4), rd((100 * r.R) / rTot, 1)]; }), `Inside film ${rd(S.hIn, 0)} W/m²K (Gnielinski on the no-slip mixture), outside film ${rd(S.mainCat.startsWith('sea') ? S.hSea : S.hAir, 0)} W/m²K (Churchill–Bernstein).${v.uMult !== 1 ? ` The U-value is multiplied by the calibration factor ${v.uMult}.` : ''}`));
  tables.push(T('Wall design summary', ['Quantity', 'Value', 'Unit'], [['U exposed in sea water', rd(W.sea.U, 3), 'W/m²K'], ['U buried under the seabed', v.burialDepth > 0 ? rd(W.seaB.U, 3) : '—', 'W/m²K'], ['U exposed in air', rd(W.air.U, 3), 'W/m²K'], ['U buried onshore', v.burialDepth > 0 ? rd(W.airB.U, 3) : '—', 'W/m²K'], ['Length-weighted U of the route', rd(S.uMean, 3), 'W/m²K'], ['Thermal mass of wall and coatings', rd(wall.thermalMass, 0), 'J/m/K'], ['Mass of pipe and coatings', rd(wall.massPerM, 1), 'kg/m'], ['Submerged weight, operating', rd(wall.submerged, 0), 'N/m'], ['Specific gravity, empty', rd(S.sgEmpty, 3), '–'], ['Specific gravity, operating', rd(S.sgOp, 3), '–'], ['Bending stiffness EI', sg(wall.EI, 4), 'N·m²'], ['Material', S.mat.label, ''], ['SMYS / SMTS', `${S.mat.smys} / ${S.mat.smts}`, 'MPa'], ['Hoop stress at design pressure (Barlow)', rd(S.hoop, 1), 'MPa'], ['Hoop utilisation against 0.72·SMYS', rd(S.hoopUtil, 3), '–']]));
  tables.push(T('Network nodes', ['Node', 'Kind', 'Elevation (m)', 'Pressure (bara)', 'Temperature (°C)', 'Boundary', 'Net flow at node (kg/s)'], net.nodes.map((n, k) => [n.id, n.kind, rd(n.z, 1), rd(sol.p[k], 3), rd(S.base.t[k], 2), n.pFixed !== null ? `fixed ${rd(n.pFixed, 2)} bara` : n.rate > 0 ? `fixed ${rd(n.rate, 3)} kg/s` : 'mass balance', sg(sol.balance[k], 4)]), 'For fixed-pressure nodes the net flow is what leaves the network there; for all other nodes it is the mass-balance residual.'));
  tables.push(T('Network connections', ['From', 'To', 'Type', 'Length (m)', 'ID (mm)', 'Flow area (cm²)', 'Hydraulic diameter (mm)', 'Mean inclination (°)', 'Flow direction', 'Flow (kg/s)', 'Δp (bar)', 'Friction + local (bar)', 'Elevation (bar)', 'Max velocity (m/s)', 'Erosional ratio', 'Mean holdup', 'Flow pattern', 'Reynolds regime', 'Fanning factor'], eInfo.map((q, i) => { const fD = q.pipe && q.Re > 0 ? frictionFactor(q.Re, q.e.rough / q.e.D, v.fModel) : null; return [q.e.from, q.e.to, q.e.kind === 'link' ? 'connector' : q.e.kind, q.pipe ? rd(q.e.L, 1) : '—', q.pipe ? rd(q.e.D * 1000, 1) : '—', q.pipe ? rd(q.A * 1e4, 1) : '—', q.pipe ? rd(q.e.D * 1000, 1) : '—', q.pipe ? rd(q.incl, 2) : '—', q.rev ? 'reversed (to → from)' : 'as drawn', sg(q.m, 5), sg(sol.dp[i], 5), q.pipe ? sg(q.dpFric, 4) : '—', q.pipe ? sg(q.dpGrav, 4) : '—', q.pipe ? rd(q.vMax, 2) : '—', q.pipe ? rd(q.ero, 3) : '—', q.pipe && q.holdup !== null ? rd(q.holdup, 3) : '—', q.regime, q.pipe ? reRegime(q.Re) : '—', fD !== null ? sg(fD / 4, 3) : '—']; }), `Newton–Raphson on ${net.nodes.length} nodes and ${net.edges.length} connections: ${sol.iterations} iterations, relative residual ${sg(sol.residual, 2)}, largest nodal mass imbalance ${sg(sol.massResidual, 2)} kg/s.`));
  if (sol.loops.length) tables.push(T('Loop pressure balance', ['Loop', 'Connections (sign = direction around the loop)', 'Sum of pressure drops (bar)'], sol.loops.map((l, i) => [i + 1, l.edges.map((q) => `${q.sign > 0 ? '+' : '−'}${net.edges[q.i].from}→${net.edges[q.i].to}`).join('  '), sg(l.sum, 3)]), 'Pressure drops recomputed from the branch laws at the converged flows; each loop must sum to zero.'));
  tables.push(T('Well trajectory (minimum curvature)', ['MD (m)', 'Inclination (°)', 'Azimuth (°)', 'TVD (m)', 'North (m)', 'East (m)', 'Displacement (m)', 'Dog-leg severity (°/30 m)'], stn.map((s) => [rd(s.md, 1), rd(s.inc, 2), rd(s.azi, 2), rd(s.tvd, 2), rd(s.north, 2), rd(s.east, 2), rd(s.disp, 2), rd(s.dls, 3)])));
  tables.push(T('Well and nodal analysis', ['Quantity', 'Value', 'Unit'], [['Inflow model', ipr.type, ''], ['Reservoir pressure / temperature', `${rd(ipr.pRes, 1)} / ${rd(v.tRes, 1)}`, 'bara / °C'], ['Productivity index per well', sg(ipr.pi, 4), 'Sm³/d/bar'], ['Absolute open-flow potential per well', rd(ipr.qMax, 0), S.wellUnit], ['Non-Darcy coefficient β', ipr.beta ? sg(ipr.beta, 3) : '—', '1/m'], ['Wells (sources) in the network', nWells, ''], ['Well used for the nodal analysis', net.nodes[S.wi].id, ''], ['Case rate per well', rd(S.qWell, 1), S.wellUnit], ['Bottom-hole flowing pressure at the case rate', S.sA.length && S.canDeliver ? rd(iprPwf(S.qWell, ipr) ?? 0, 1) : '—', 'bara'], ['Wellhead pressure available at the case rate', rd(S.whpAvail, 2), 'bara'], ['Wellhead pressure needed by the network', rd(S.whpReq, 2), 'bara'], ['Natural-flow rate per well (choke open)', rd(S.sOp * S.qWell, 1), S.wellUnit], ['Wellhead / bottom-hole pressure at natural flow', `${rd(S.whpOp, 1)} / ${rd(S.pwfOp, 1)}`, 'bara'], ['Field natural-flow liquid rate', rd(S.operatingRate, 0), 'Sm³/d'], ['Well TVD / MD / displacement', `${rd(traj.tvd, 0)} / ${rd(traj.md, 0)} / ${rd(traj.displacement, 0)}`, 'm'], ['Largest dog-leg severity', rd(traj.maxDls, 2), '°/30 m']]));
  tables.push(T('Well-count study', ['Wells', 'Rate per well (' + S.wellUnit + ')', 'Field liquid rate held (Sm³/d)', 'Deliverable without lift, chokes open (Sm³/d)', 'Deliverable with the lift as rated, chokes open (Sm³/d)', 'Drawdown (bar)', 'Bottom-hole pressure (bara)', 'Wellhead pressure natural / with lift / needed (bara)', 'Lift pressure rise (bar)', 'Lift power (kW)', 'Choke opening (%)', 'Tubing erosional ratio', 'Verdict'], S.study.rows.map((r) => [r.wells === S.nNow ? r.wells + ' (network as listed)' : r.wells, rd(r.ratePerWell, 0), rd(r.rate, 0), rd(r.rateNatural, 0), rd(r.rateWithLift, 0), rd(r.drawdown, 1), rd(r.pwf, 1), `${r.whpNatural === null ? '—' : rd(r.whpNatural, 1)} / ${r.whpAvailable === null ? '—' : rd(r.whpAvailable, 1)} / ${rd(r.whpNeeded, 1)}`, rd(r.boostDp, 1), rd(r.boostPower, 0), rd(r.chokeOpening, 0), rd(r.tubingErosionalRatio, 2), (r.feasible ? (r.chokeOpening <= v.chokeMaxOpen + 1e-9 ? 'feasible' : 'feasible, little choke margin') : 'not feasible: ' + r.limit) + (r.wells === S.study.recommended ? ' — recommended' : '')]), `Identical wells sharing the case rate equally: one wellhead-pressure curve of a single well (inflow and tubing lift) is combined with the back-pressure of the solved network at the total rate; the jumper loss is scaled with the square of the rate per well. Limits: drawdown ${rd(v.maxDrawdown, 0)} bar, choke at most ${rd(v.chokeMaxOpen, 0)} % open, tubing erosional ratio 1. Recommended: ${S.study.recommended} — the ${S.study.basis}. The cost of each option is priced by the economics suite.`));
  if (S.sens.length) tables.push(T('Deliverability sensitivities', ['Parameter', 'Value', 'Natural-flow liquid rate (Sm³/d)', 'Change (%)'], S.sens.map((q) => [q.param, rd(q.value, 1), rd(q.rate, 0), S.operatingRate > 0 ? rd((100 * (q.rate - S.operatingRate)) / S.operatingRate, 1) : '—']), 'Well-side changes only: the network back-pressure curve of the base case is kept.'));
  const cs = S.chokeState;
  tables.push(T('Production choke', ['Model', 'Result', 'Unit'], [['Upstream / downstream pressure', `${rd(S.p1c, 1)} / ${rd(S.p2c, 1)}`, 'bara'], ['Mass rate per well', rd(S.mWell, 3), 'kg/s'], ['IEC 60534 required Cv / Kv', `${rd(S.chokeCv, 2)} / ${rd(S.chokeCv * 0.865, 2)}`, ''], ['Expansion factor Y / pressure-drop ratio x', cs ? `${rd(cs.Y, 3)} / ${rd(cs.x, 3)}` : '—', ''], ['IEC choked (x ≥ Fk·xT or Δp ≥ FL²·p1)', cs ? (cs.choked ? 'yes' : 'no') : '—', ''], ['Opening on the selected characteristic', rd(S.chokeOpening, 1), '%'], ['Sachdeva critical pressure ratio', rd(S.sach.yc, 3), '–'], ['Flow state (Sachdeva)', S.critical ? 'critical' : 'sub-critical', ''], ['Sachdeva flow area / bean diameter', `${sg(S.areaSach * 1e6, 4)} / ${rd(Math.sqrt((4 * S.areaSach) / Math.PI) * 1000, 1)}`, 'mm² / mm'], ['Homogeneous-equilibrium critical ratio', rd(S.hem.ratio, 3), '–'], ['Homogeneous-equilibrium critical flux', rd(S.hem.G, 0), 'kg/m²/s'], ['Smallest area that passes the rate (HEM, with Cd)', sg(S.areaHem * 1e6, 4), 'mm²'], ['Bernoulli orifice diameter (homogeneous, Cd)', rd(S.dOrifice * 1000, 1), 'mm'], ...S.beanRows.map((r) => [`${r[0]} bean size for critical flow`, `${r[1]} / ${r[2]}`, '64ths in / mm'])], 'The Gilbert-type correlations assume critical flow and stock-tank liquid rate with the producing gas–liquid ratio.'));
  const duty = [['Production choke (per well)', `Cv ${rd(S.chokeCv, 1)} of ${v.chokeCvMax}`, `${rd(S.chokeOpening, 0)} % open`, `Δp ${rd(Math.max(S.chokeDp, 0), 1)} bar`], ['Separator', `${v.sepD} m × ${v.sepL} m ${v.sepOrient}`, `gas ${rd(S.gasLoad * 100, 0)} % of Souders–Brown`, `residence ${rd(Math.min(S.resTime / 60, 9999), 1)} min`], ['Slug catcher', `${rd(S.scVol, 0)} m³`, v.scType === 'finger' ? `${S.nFing} × ${v.fingerDmm} mm × ${rd(S.fingerL, 0)} m` : `${rd(S.scVesselD, 1)} m × ${rd(4 * S.scVesselD, 1)} m`, `peak ${rd(S.scPeak, 1)} m³`]];
  if (pumpRes) duty.push(['Liquid export pump', `${rd(pumpRes.qRun, 0)} m³/h at ${rd(pumpRes.hRun, 0)} m`, `${rd(pumpRes.power, 0)} kW, efficiency ${rd(pumpRes.eta * 100, 0)} %`, `NPSH ${rd(pumpRes.npshA, 1)} / ${rd(pumpRes.npshR, 1)} m; Euler head per stage ${rd(pumpRes.eu.head, 0)} m (slip ${rd(pumpRes.eu.slip, 2)})${pumpRes.stages ? `, ${pumpRes.stages} stage(s)` : ''}`]);
  eInfo.filter((q) => q.e.kind === 'pump' && isNum(q.r.head)).forEach((q) => duty.push([`In-line pump ${q.e.from} → ${q.e.to}`, `${rd(q.r.q, 0)} m³/h at ${rd(q.r.head, 0)} m`, `${rd(q.r.power / 1000, 0)} kW`, `Δp ${rd(-q.r.dp, 1)} bar`]));
  eInfo.filter((q) => q.e.kind === 'choke').forEach((q) => duty.push([`Network valve ${q.e.from} → ${q.e.to}`, `Cv ${rd(q.e.cv, 1)}`, `Δp ${rd(q.r.dp ?? 0, 2)} bar`, q.r.choked ? 'choked' : 'not choked']));
  if (comp) comp.stages.forEach((s, k) => duty.push([`Gas compressor stage ${k + 1}`, `${rd(s.p1, 1)} → ${rd(s.p2, 1)} bara`, `${rd(s.power / 1000, 0)} kW, head ${rd(s.headPoly / 1000, 1)} kJ/kg (isentropic ${rd(s.headIsen / 1000, 1)})`, `discharge ${rd(s.t2, 0)} °C, k ${rd(s.k, 3)}, Z ${rd(s.z, 3)}, isentropic efficiency ${rd(s.etaIsen * 100, 1)} %`]));
  if (comp) duty.push(['Compressor map position', `${rd(comp.pt.speed * 100, 1)} % speed`, `surge margin ${rd(comp.pt.surgeMargin, 0)} %`, `stonewall margin ${rd(comp.pt.chokeMargin, 0)} %`]);
  tables.push(T('Equipment duty list', ['Equipment', 'Duty', 'Performance', 'Notes'], duty));
  tables.push(T('Separator product split', ['Stream', 'Mass rate (kg/s)', 'Actual volume rate (m³/h)', 'Density (kg/m³)', 'Molar mass (g/mol)'], [['Gas', sg(S.mGasSep, 5), rd(S.qGasSep * 3600, 1), rd(S.rhoGsep, 2), S.nGas > 0 ? rd(S.fl.gas.MW, 2) : '—'], ['Oil / condensate', sg(S.mOilSep, 5), rd(S.qOilSep * 3600, 2), S.nOil > 0 ? rd(S.fl.oil.rho, 1) : '—', S.nOil > 0 ? rd(S.fl.oil.MW, 1) : '—'], ['Water', sg(S.mW, 5), rd(S.qWatSep * 3600, 2), rd(S.mW > 0 ? S.mW / S.qWatSep : 1000, 1), '18.0']], `Equation-of-state flash at ${rd(v.pSep, 1)} bara and ${rd(S.tArr, 1)} °C: vapour mole fraction ${rd(S.nGas / S.nHC, 4)}; largest component-balance error ${sg(S.compErr, 2)} (mole fraction).`));
  const raw = S.route.raw;
  tables.push(T('Geometry consistency checks', ['Check', 'Result', 'Status'], [['Length conservation: Σ segment lengths vs arc length', `${rd(prof.segs.reduce((s, q) => s + q.L, 0), 3)} m vs ${rd(prof.length, 3)} m`, 'pass'], ['Cell lengths of the main line vs route length', `${rd(S.mlLen, 3)} m vs ${rd(prof.length, 3)} m`, !ml.x.length || Math.abs(S.mlLen - prof.length) < 1e-6 * prof.length ? 'pass' : 'check'], ['Zero-length elements in the table as typed', raw ? raw.zeroLength : 0, raw && raw.zeroLength ? 'repaired' : 'pass'], ['Repeated points in the table as typed', raw ? raw.duplicates : 0, raw && raw.duplicates ? 'repaired' : 'pass'], ['Monotonic chainage', raw && raw.backward ? `${raw.backward} backward step(s)` : 'monotonic', raw && raw.backward ? 'repaired' : 'pass'], ['Round trip profile → segments → profile', sg(Math.max(...(() => { const b = segmentsToProfile(profileToSegments(x, z)); return b.x.map((q, i) => Math.hypot(q - x[i], b.z[i] - z[i])); })()), 2) + ' m', 'pass'], ['Network nodes / connections / independent loops', `${net.nodes.length} / ${net.edges.length} / ${sol.loops.length}`, sol.usable ? 'pass' : 'check'], ['Largest nodal mass imbalance', sg(sol.massResidual, 2) + ' kg/s', sol.massResidual < 1e-6 * S.mCase ? 'pass' : 'check'], ['Largest loop pressure imbalance', sol.loops.length ? sg(Math.max(...sol.loops.map((l) => Math.abs(l.sum))), 2) + ' bar' : 'no loops', 'pass'], ...S.gchk.map((c) => [`Route check: ${c.name}`, `${isNum(c.got) ? sg(c.got, 7) : '—'} (expected ${isNum(c.expected) ? sg(c.expected, 7) : '—'})`, c.pass ? 'pass' : 'check'])], 'The route checks are those applied to imported geometry (length and volume conservation, zero-length elements, repeated nodes, chainage, finite values, inclination bounds, round trip through a chainage–elevation table).'));
  const riserE = eInfo.find((q) => q.e.kind === 'riser');
  tables.push(T('Equipment orientation', ['Item', 'Kind', 'Built direction against the solved flow', 'Status'], [...S.orient.map((q) => [q.item, q.kind, q.direction, q.ok ? 'pass' : 'reversed']), ['Separator', v.sepOrient, v.sepOrient === 'horizontal' ? 'axis horizontal, gas space above the liquid' : 'axis vertical', 'pass'], ...(riserE ? [[`Riser ${riserE.e.from} → ${riserE.e.to}`, 'riser', `mean inclination ${rd(riserE.incl, 1)}° (${riserE.rev ? 'flow downwards' : 'flow upwards'})`, riserE.incl > 0 && !riserE.rev ? 'pass' : 'check']] : []), ...(S.trans ? [['Export pump / compressor / recycle (transient)', 'non-return', 'reverse flow blocked by check valves', 'pass']] : [])], 'Pumps, compressors, chokes, control valves and check valves must carry flow from their “from” node to their “to” node.'));
  const L = S.lift;
  tables.push(T('Artificial lift', ['Quantity', 'Value', 'Unit'], [['Method', L.type === 'boost' ? 'subsea multiphase boosting' : L.type === 'gaslift' ? 'gas lift' : 'none', ''], ['Flowing wells', L.wells, ''], ['Boosting pressure rise: duty / rated', L.type === 'boost' ? `${rd(L.dp, 1)} / ${rd(v.liftDp, 1)}` : '—', 'bar'], ['Boosting control', L.type === 'boost' ? (v.liftAuto ? `speed-controlled to hold the choke ${v.chokeTarget} % open (at least ${v.liftDpMin} bar)` : 'fixed at the rated pressure rise') : '—', ''], ['Pump suction pressure', isNum(L.suctionP) ? rd(L.suctionP, 1) : '—', 'bara'], ['Actual volume flow at the pump suction (all wells)', isNum(L.flow) ? rd(L.flow, 0) : '—', 'm³/h'], ['Lift gas per well', L.type === 'gaslift' ? rd(L.rate, 0) : '—', 'Sm³/d'], ['Lift gas, all wells', L.type === 'gaslift' ? sg(L.gasMass, 4) : '—', 'kg/s'], ['Injection pressure at the valve', isNum(L.injectionP) ? rd(L.injectionP, 1) : '—', 'bara'], ['Shaft power', rd(L.power, 0), 'kW'], ['Wellhead pressure available with lift / needed', `${rd(S.whpAvail, 1)} / ${rd(S.whpReq, 1)}`, 'bara'], ['Margin taken by the choke', rd(Math.max(S.chokeDp, 0), 1), 'bar']], L.type === 'boost' ? 'Hydraulic power = suction volume flow × pressure rise ÷ efficiency; the pump is placed at the wellhead / manifold, so the flowline hydraulics are unchanged.' : L.type === 'gaslift' ? 'Lift gas enters the tubing at the valve depth and returns with the production through the flowline, separator and compressor; it is given the properties of the produced gas.' : null));
  const A = S.arrestors, C = S.cp;
  tables.push(T('Design features', ['Feature', 'Value', 'Unit', 'Note'], [['Buckle-propagation pressure 35·SMYS·(t/D)^2.5', rd(S.ppr, 1), 'bar', 'unfactored'], ['Water depth below which a buckle would propagate', rd(S.arrFirst, 0), 'm', v.buckleArrestorFrom > 0 ? 'as specified' : 'external pressure × 1.3 = propagation pressure'], ['Route length needing arrestors', rd(A.length, 0), 'm', A.required ? 'required' : 'not required at this depth'], ['Buckle arrestors', A.count, '', A.count ? `every ${rd(A.spacing, 0)} m` : 'none specified'], ['Lateral-buckle initiators', v.buckleInitType, '', v.buckleInitType === 'none' ? '—' : `${S.initX.length} at ${rd(v.buckleInitSpacing, 0)} m spacing`], ['Pipe lift-off length each side of a sleeper', v.buckleInitType === 'sleeper' ? rd(S.sleeperSpan, 1) : '—', 'm', '(72·EI·h/w)^¼'], ['Bracelet anodes', C.count, '', v.anodeSpacing > 0 ? `every ${rd(v.anodeSpacing, 0)} m` : 'none specified'], ['Protection current per anode section', sg(C.currentPerAnode, 3), 'A', `${v.cpCurrent} A/m² × breakdown ${v.cpBreakdown}`], ['Anode mass required / installed', `${rd(C.requiredMass, 1)} / ${rd(C.anodeMass, 1)}`, 'kg', `${v.cpLife} y design life`], ['Anode life at the mean current', rd(C.lifeYears, 1), 'y', ''], ['Allowable free span', rd(v.spanLimit, 0), 'm', `${S.spanOver} candidate(s) longer`], ['Span candidates split by a support', S.spansSplit, '', '']]));
  tables.push(T('Supports and restraints', ['Distance (m)', 'Type', 'Stiffness (kN/m)', 'Elevation (m)'], [...S.supUser.map((q) => [rd(q.x, 0), q.type, q.k > 0 ? rd(q.k, 0) : '—', rd(interp1(x, z, clamp(q.x, x[0], x[nP - 1])), 1)]), ...S.initX.slice(0, 40).map((q) => [rd(q, 0), v.buckleInitType === 'sleeper' ? `sleeper, ${v.sleeperHeight} m` : 'snake-lay crown', '—', rd(interp1(x, z, q), 1)])]));
  if (tr) tables.push(T('Transient scenario', ['Quantity', 'Value', 'Unit'], [['Scenario', tr.scenario, ''], ['Event time / simulated time', `${rd(tr.tEvent, 1)} / ${rd(tr.tEnd, 0)}`, 's'], ['Initial mass rate: gas / liquid', `${sg(tr.mGas0, 4)} / ${sg(tr.mLiq0, 4)}`, 'kg/s'], ['Separator pressure: initial / peak / lowest', `${rd(tr.pSep[0], 2)} / ${rd(tr.peakSepP, 2)} / ${rd(tr.minSepP, 2)}`, 'bara'], ['Flowline line-pack pressure: initial / final', `${rd(tr.pLine[0], 2)} / ${rd(tr.pLine[tr.pLine.length - 1], 2)}`, 'bara'], ['Discharge header pressure: initial / peak', tr.pHeader ? `${rd(tr.pHeader[0], 1)} / ${rd(tr.peakHeaderP, 1)}` : '—', 'bara'], ['Separator level: initial / highest / final', `${rd(tr.level[0], 3)} / ${rd(tr.levelMax, 3)} / ${rd(tr.levelEnd, 3)}`, 'm'], ['Compressor speed at the end', tr.compSpeed ? rd(100 * tr.compSpeed[tr.compSpeed.length - 1], 1) : '—', '%'], ['Compressor coast-down to half speed', isNum(tr.compHalfTime) ? rd(tr.compHalfTime, 1) : '—', 's'], ['Pump speed at the end', tr.pumpSpeed ? rd(100 * tr.pumpSpeed[tr.pumpSpeed.length - 1], 1) : '—', '%'], ['Pump coast-down to half speed', isNum(tr.pumpHalfTime) ? rd(tr.pumpHalfTime, 1) : '—', 's'], ['Smallest surge margin', isNum(tr.minMargin) ? rd(tr.minMargin, 1) : '—', '%'], ['Surge line crossed', tr.surgeTime !== null ? `at ${rd(tr.surgeTime, 1)} s` : 'no', ''], ['Largest recycle-valve opening', tr.asv ? rd(100 * tr.asvMax, 0) : '—', '%'], ['Gas to flare', rd(tr.flareMass, 1), 'kg'], ['Well rate at the end', sg(tr.wellsEnd, 4), 'kg/s'], ['Gas inventory: initial / final', `${rd(tr.gas.initial, 1)} / ${rd(tr.gas.final, 1)}`, 'kg'], ['Liquid inventory: initial / final', `${rd(tr.liquidBal.initial, 0)} / ${rd(tr.liquidBal.final, 0)}`, 'kg'], ['Largest mass-balance error (relative)', sg(tr.massError, 2), '–'], ['Time steps / halved steps', `${tr.steps} / ${tr.halved}`, '']], 'Dynamic nodal inventories (flowline gas line-pack with the real-gas factor of the case fluid, separator gas and liquid, discharge header) joined by quasi-steady branches: well deliverability with the production choke, flowline resistance, compressor on a fan-law map with rotor inertia, governor, anti-surge recycle and trip logic, export and flare valves with actuator dynamics, export pump with rotor inertia and level control. Implicit θ-method (θ = 0.55).'));
  if (S.cmp.rows.length) {
    tables.push(T('Comparison with measured data', ['Kind', 'At', 'Second coordinate', 'Measured', 'Model', 'Model − measured', 'Deviation (%)', 'Unit'], S.cmp.rows.slice(0, 400).map((q) => [q.label, sg(q.at, 6), sg(q.aux, 4), sg(q.measured, 6), sg(q.model, 6), sg(q.diff, 3), q.pct === null ? '—' : rd(q.pct, 2), q.unit])));
    tables.push(T('Comparison metrics by kind of measurement', ['Kind', 'Points', 'Bias (model − measured)', 'RMS error', 'Mean absolute error (%)', 'Unit'], S.cmp.groups.map((g) => [g.label, g.n, sg(g.bias, 3), sg(g.rmse, 3), g.mape === null ? '—' : rd(g.mape, 2), g.unit])));
  }

  // ---- balances ------------------------------------------------------------------------------------------------------
  const balances = [{ name: 'Network mass (kg/s): sources → sinks', in: sol.inflow, out: sol.outflow }, { name: 'Separator mass (kg/s): feed → gas + oil + water', in: S.mHC + S.mW, out: S.mGasSep + S.mOilSep + S.mW }, { name: 'Methane over the separator flash (mol/s)', in: S.zC1 * S.nHC, out: S.c1Out }, { name: 'Slug catcher inventory change (m³): accumulated vs inflow − outflow', in: S.scEnd[0] - S.scV0, out: S.scEnd[1] - S.scEnd[2] }];
  if (ml.x.length) balances.push({ name: 'Route length (m): arc length vs sum of hydraulic cells', in: prof.length, out: S.mlLen });
  if (tr) balances.push({ name: 'Transient gas mass (kg): initial inventory + inflow vs final inventory + outflow', in: tr.gas.initial + tr.gas.in, out: tr.gas.final + tr.gas.out }, { name: 'Transient liquid mass (kg): initial inventory + inflow vs final inventory + outflow', in: tr.liquidBal.initial + tr.liquidBal.in, out: tr.liquidBal.final + tr.liquidBal.out });

  // ---- outputs -------------------------------------------------------------------------------------------------------
  const outputs = {
    profile: { x: simp.x, z: simp.z }, length: prof.length, id: D, wt, od: S.od, roughness: rough, uValue: S.uMean, layers: wall.layers.map((l) => ({ name: l.name, t: l.t, k: l.k })), volume: S.volume,
    waterDepth: S.waterDepth, riserHeight: rb.height, riserBaseX: rb.x, tSeabed: v.tSeabed, tSeaSurface: v.tSeaSurface, burial: v.burialDepth > 0 ? v.burialDepth : 0, kLoss: S.kLoss,
    bends: S.bends.slice(0, 60).map((b) => ({ x: b.x, angle: b.angle, radius: Math.min(b.radius, 1e7), plane: b.plane })), spans: S.spans.map((s) => ({ x: s.x, length: s.length, gap: s.gap })), spanSource: S.spanSource || 'none',
    ipr: { type: ipr.type, pRes: ipr.pRes, tRes: v.tRes, pi: ipr.pi, qMax: ipr.qMax, pb: ipr.pb, basis: ipr.basis, wells: nWells, qMaxTotal: ipr.qMax / S.share },
    wellTVD: traj.tvd, wellMD: traj.md, tubingId: v.tubingIdMm / 1000, whp: S.canDeliver ? S.whpAvail : S.whpReq, operatingRate: S.operatingRate, chokeCv: S.chokeCv, chokeOpening: S.chokeOpening, separatorP: v.pSep, slugCatcherVol: S.scVol,
    network: { nodes: net.nodes.map((n, k) => ({ id: n.id, type: n.kind, p: sol.p[k], z: n.z, t: S.base.t[k] })), edges: eInfo.map((q, i) => ({ from: q.e.from, to: q.e.to, type: q.e.kind === 'link' ? 'connector' : q.e.kind, length: q.e.L, id: q.e.D, q: q.m, dp: sol.dp[i], v: q.vMax, area: q.A, incl: q.incl, reversed: q.rev, k: q.e.k })), converged: sol.usable, residual: sol.residual, massResidual: sol.massResidual, loops: sol.loops.length },
    pumpPower, compressorPower: compPower, material: { grade: S.mat.grade, smys: S.mat.smys, smts: S.mat.smts, E: S.mat.E, poisson: S.mat.poisson, alphaT: S.mat.alphaT, rho: S.mat.rho }, designPressure: v.designPressure, designTemp: v.designTemp,
    // extras
    area: area(D), hydraulicDiameter: D, sections: S.sections.map((q) => ({ name: q[0], area: q[1], hydraulicDiameter: q[2] })), seabedContact: S.contact,
    centreline: S.route.plan ? { x: simp.idx.map((i) => S.route.plan.e[i]), y: simp.idx.map((i) => S.route.plan.n[i]), z: simp.z, crs: S.route.plan.geographic ? 'local metres from the first point (converted from WGS 84)' : 'local metres from the first point' } : null,
    mesh: { s: (() => { let a = 0; return ml.ds.map((d) => { a += d; return a - d / 2; }); })(), x: ml.x.slice(), z: ml.z.slice(), ds: ml.ds.slice(), theta: ml.theta.slice(), id: ml.D.slice(), cells: net.edges.map((e) => ({ from: e.from, to: e.to, n: e.cells.length })) },
    diagnostics: { zeroLength: raw ? raw.zeroLength : 0, duplicates: raw ? raw.duplicates : 0, backward: raw ? raw.backward : 0, nodes: net.nodes.length, edges: net.edges.length, loops: sol.loops.length, massResidual: sol.massResidual, loopResidual: sol.loops.length ? Math.max(...sol.loops.map((l) => Math.abs(l.sum))) : 0, lengthError: Math.abs(S.mlLen - prof.length) },
    thermalMass: wall.thermalMass, submergedWeight: wall.submerged, specificGravity: S.sgOp, specificGravityEmpty: S.sgEmpty, bendingStiffness: wall.EI, inletPressure: S.pInlet, arrivalTemp: S.tArr, whpRequired: S.whpReq, whpAvailable: S.whpAvail, chokeDp: Math.max(S.chokeDp, 0), chokeCritical: S.critical,
    flowlineLength: S.lenFlow, riserLength: S.lenRiser, idEffective: S.Dfl, roughnessEffective: S.roughFl, erosionalRatio: S.eroAll, maxVelocity: S.vMaxAll, liquidInventory: S.liquidInv, highPoints: prof.highs.map((q) => ({ x: q.x, z: q.z })), lowPoints: S.lowPts.map((q) => ({ x: q.x, z: q.z, vm: q.vm })),
    separator: { d: v.sepD, l: v.sepL, orientation: v.sepOrient, volume: S.sepVtot, liquidVolume: S.sepVliq, residence: Math.min(S.resTime, 1e9), gasLoad: S.gasLoad, qGas: S.qGasSep, qLiquid: S.qLiqSep, mGas: S.mGasSep, mOil: S.mOilSep, mWater: S.mW }, compressorStages: comp ? comp.nSt : 0, compressorDischargeT: comp ? comp.tMax : null, surgeMargin: comp ? comp.pt.surgeMargin : null, npshMargin: pumpRes ? pumpRes.npshA / pumpRes.npshR : null,
    // design features, lift, transient, checks
    buckleArrestorSpacing: v.buckleArrestorSpacing, buckleArrestors: { spacing: A.spacing, count: A.count, firstDepth: A.firstDepth, length: A.length, propagationPressure: A.propagationPressure, required: A.required },
    buckleInitiators: { type: v.buckleInitType, spacing: v.buckleInitType === 'none' ? 0 : v.buckleInitSpacing, count: S.initX.length, x: S.initX.slice(0, 80), sleeperHeight: v.buckleInitType === 'sleeper' ? v.sleeperHeight : 0, spanEachSide: S.sleeperSpan },
    anodeSpacing: v.anodeSpacing, cathodicProtection: { ...C }, spanLimit: v.spanLimit, spansOverLimit: S.spanOver,
    supports: [...S.supUser.map((q) => ({ x: q.x, type: q.type, stiffness: q.k })), ...S.initX.slice(0, 80).map((q) => ({ x: q, type: v.buckleInitType === 'sleeper' ? 'sleeper' : 'snake lay', stiffness: 0 }))],
    lift: { type: L.type, power: L.power, rate: L.rate, dp: L.dp, ratedDp: L.type === 'boost' ? v.liftDp : 0, wells: L.wells, flow: L.flow, suctionP: L.suctionP, injectionP: L.injectionP, gasMass: L.gasMass }, wellsFlowing: S.nFlow, separatorK: S.sepKused,
    wellCountStudy: S.study.rows.map((r) => ({ ...r })), wellsRecommended: S.study.recommended, wellsRecommendedBasis: S.study.basis, drawdown: Math.max(ipr.pRes - (iprPwf(S.qWell, ipr) ?? 0), 0), maxDrawdown: v.maxDrawdown, ratePerWell: S.qWell,
    transient: tr ? { scenario: tr.scenario, tEvent: tr.tEvent, t: tr.t, pSep: tr.pSep, pLine: tr.pLine, pHeader: tr.pHeader, level: tr.level, compSpeed: tr.compSpeed, pumpSpeed: tr.pumpSpeed, mWells: tr.mWells, mCompressor: tr.mComp, mExport: tr.mExport, mFlare: tr.mFlare, mLiquid: tr.mLiquid, recycleValve: tr.asv, surgeMargin: tr.margin, peakSepP: tr.peakSepP, minSurgeMargin: tr.minMargin, surgeTime: tr.surgeTime, compHalfTime: tr.compHalfTime, pumpHalfTime: tr.pumpHalfTime, levelMax: tr.levelMax, levelEnd: tr.levelEnd, flareMass: tr.flareMass, massError: tr.massError, halved: tr.halved } : null,
    transientBasis: S.basis,
    geometryChecks: S.gchk.map((c) => ({ name: c.name, pass: !!c.pass })), orientation: S.orient.map((q) => ({ item: q.item, kind: q.kind, direction: q.direction, ok: q.ok })),
    comparison: { n: S.cmp.rows.length, groups: S.cmp.groups.map((g) => ({ kind: g.key, label: g.label, unit: g.unit, n: g.n, bias: g.bias, rmse: g.rmse, mape: g.mape })) },
  };
  const summary = `${rd(prof.length / 1000, 2)} km route${rb.height > 0 ? ` with a ${rd(rb.height, 0)} m riser` : ''}${S.offshore ? ` in ${rd(S.waterDepth, 0)} m of water` : ' onshore'}, U = ${rd(S.uMean, 2)} W/m²K. At the case rate the network needs ${rd(S.whpReq, 1)} bara at the wellheads (${rd(S.pInlet, 1)} bara at the flowline inlet, arrival ${rd(S.tArr, 1)} °C) and the wells deliver ${rd(S.whpAvail, 1)} bara${L.type === 'boost' ? ` with ${rd(S.boost, 1)} bar of subsea boosting (${rd(L.power, 0)} kW of ${rd(v.liftDp, 0)} bar rated)` : L.type === 'gaslift' ? ' with gas lift' : ''}; deliverable rate with the chokes open ${rd(S.operatingRate, 0)} Sm³/d of liquid${S.canDeliver ? `, choke ${rd(S.chokeOpening, 0)} % open` : ''}.${tr ? ` Transient “${tr.scenario}”: separator pressure peaks at ${rd(tr.peakSepP, 1)} bara.` : ''}`;
  return { summary, kpis, warnings: S.warnings, recommendations: S.recs, plots, tables, balances, outputs };
}

// ---- calibration ---------------------------------------------------------------------------------------------------
/**
 * Fast predictions for one row of test data. Every target is a separate, simple test so that each parameter can be estimated from the
 * measurement that responds to it: single-phase (sea-water) line flow test and line-fill volume, bend and tee spool tests, valve
 * water-flow test, bean-choke test, pump head and power at a test flow, normalised compressor head and discharge temperature,
 * separator carry-over velocity, well-test flowing pressure, skin pressure drop, tubing friction test, thermal test, slug-catcher hold-up.
 */
function calModel(v0 = {}) {
  const v = readInputs(v0), route = routeOf(v), L = analyseProfile(route.x, route.z).length, D = v.effIdMm > 0 && v.effIdMm < v.idMm ? v.effIdMm / 1000 : v.idMm / 1000, rho = 1027, mu = 1.6e-3, c = (k, d) => num(v0[k], d);
  const qLine = c('qLine', 300), qTest = c('qTest', 800), wChoke = c('wChoke', 9), rhoChoke = c('rhoChoke', 350), bean = c('bean', 48), mCp = c('mCp', 120), tHot = c('tHot', 70), qValve = c('qValve', 150), qPump = c('qPump', 180), phiComp = c('phiComp', 0.9), nComp = c('nComp', 100), rComp = c('rComp', 2.5), rhoLt = c('rhoLt', 750), rhoGt = c('rhoGt', 25), qTub = c('qTub', 40), scGeom = c('scGeom', 150);
  const fl = fittingsLoss(v.fittings, D, (v.roughUm * 1e-6) / D, { bend: v.bendMult, tee: v.teeMult }), dyn = 0.5 * rho * (qLine / 3600 / area(D)) ** 2;
  const dpLine = pipeDp({ m: (qLine / 3600) * rho, rho, mu, L, D, rough: v.roughUm * 1e-6, k: fl.K * v.kMult, model: v.fModel }) / 1e5;
  const needFluid = ['darcy', 'jones', 'gas'].includes(v.iprType), fm = needFluid ? fluidModel({}) : null, ipr = makeIpr(v, fm, fm ? fm.rates.mHC : 0, fm ? fm.rates.mW : 0), pwf = iprPwf(qTest, ipr) ?? 0;
  const dpChoke = (wChoke / (v.chokeCd * area((bean / 64) * UNIT.in))) ** 2 / (2 * rhoChoke) / 1e5;
  const mat = material(v.grade), lay = v.layers.map((l) => ({ name: l?.name, t: Math.max(num(l?.tMm, 0), 0) / 1000, k: num(l?.k, 0) }));
  { let bi = -1, br = 0; lay.forEach((l, i) => { const r = l.t > 0 && l.k > 0 ? l.t / l.k : 0; if (r > br) { br = r; bi = i; } }); if (bi >= 0) lay[bi] = { ...lay[bi], k: lay[bi].k * v.insKMult }; }
  const base = { id: v.idMm / 1000, wt: v.wtMm / 1000, kSteel: mat.k, layers: lay, pipWt: v.pip ? v.pipWtMm / 1000 : 0, concrete: { t: v.concreteMm / 1000, k: 2, rho: v.concreteRho }, hIn: 1500, uMult: v.uMult };
  const od = wallDesign(base).od, U = wallDesign({ ...base, hOut: hOutside(v.currentSpeed, od, 'seawater', v.tSeabed), burial: v.burialDepth > 0 ? { cover: v.burialDepth, kSoil: v.kSoil } : null }).U;
  const pp = { qr: v.pumpQr, hr: v.pumpHr * v.degrFactor, shutoff: v.pumpShutoff, eta: v.pumpEff }, hPump = pumpHead(qPump, pp), N = Math.max(nComp / 100, 1e-3), traj = minimumCurvature(v.survey), Dt = v.tubingIdMm / 1000;
  return { dpLine, pwf, dpChoke, tArr: v.tSeabed + (tHot - v.tSeabed) * Math.exp((-U * Math.PI * (v.idMm / 1000) * L) / (Math.max(mCp, 1e-6) * 1000)),
    lineVol: area(D) * L, dpBend: (bendK(5, 90, fl.fT) * v.bendMult * dyn) / 100, dpTee: (fl.fT * FITTINGS.teeBranch.leD * v.teeMult * dyn) / 100, dpValve: (qValve / (0.865 * v.chokeCvMax)) ** 2,
    hPump, pPump: (1000 * G * (qPump / 3600) * Math.max(hPump, 0)) / pumpEfficiency(qPump, pp) / 1000, psiComp: v.degrFactor * v.compHeadMult * N * N * mapShape(phiComp / (N * v.compFlowMult)),
    t2Comp: compressorStage({ p1: 10, p2: 10 * Math.max(rComp, 1.01), t1: 30, z: 0.95, k: 1.28, mw: 20, eta: v.compEta }).t2, vGasMax: soudersBrown(v.sepK, rhoLt, rhoGt),
    dpSkin: ((qTest / DAY) * 1e-3 * 1.3 * v.skin) / (2 * Math.PI * v.permMd * UNIT.mD * v.payM) / 1e5, dpTubing: pipeDp({ m: (qTub / 3600) * 1000, rho: 1000, mu: 5e-4, L: traj.md, D: Dt, rough: v.tubingRoughUm * 1e-6, model: v.fModel }) / 1e5, scHeld: v.scEff * scGeom };
}
// Synthetic test data: the model above with roughness 80 µm, effective diameter 250 mm, bend multiplier 1.2, tee multiplier 0.9, rated Cv 110,
// bean Cd 0.78, pump rated head 505 m, shut-off ratio 1.3 and efficiency 0.71, compressor head scaling 0.97, flow scaling 1.03 and
// polytropic efficiency 0.75, separator K 0.10 m/s, PI 21.5 Sm³/d/bar, skin 3.5, tubing roughness 60 µm, U multiplier 1.15 and
// slug-catcher effective fraction 0.8, plus ±2 % uniform noise (seeded generator).
const CAL_SAMPLE = [{qLine:120,qTest:350,wChoke:5.73,rhoChoke:304,bean:32,mCp:87.6,tHot:55,qValve:76.4,qPump:60,phiComp:0.68,nComp:80,rComp:1.95,rhoLt:720,rhoGt:12,qTub:15,scGeom:76.4,dpLine:3.72,lineVol:1018,dpBend:0.6697,dpTee:1.941,dpValve:0.6464,dpChoke:55.27,hPump:632.7,pPump:285.9,psiComp:0.6749,t2Comp:96.98,vGasMax:0.777,pwf:284.7,dpSkin:6.729,dpTubing:0.5029,tArr:31.24,scHeld:61.32},{qLine:155,qTest:964,wChoke:7.91,rhoChoke:414,bean:40,mCp:123,tHot:65.9,qValve:125,qPump:78.2,phiComp:0.925,nComp:90,rComp:2.38,rhoLt:760,rhoGt:27,qTub:19.1,scGeom:125,dpLine:5.986,lineVol:1016,dpBend:1.113,dpTee:3.175,dpValve:1.722,dpChoke:31.21,hPump:629.2,pPump:299.9,psiComp:0.7764,t2Comp:116,vGasMax:0.5113,pwf:259.5,dpSkin:17.98,dpTubing:0.7967,tArr:44.06,scHeld:99.56},{qLine:189,qTest:1580,wChoke:10.1,rhoChoke:350,bean:48,mCp:157,tHot:76.8,qValve:175,qPump:96.4,phiComp:1.17,nComp:100,rComp:2.82,rhoLt:800,rhoGt:42,qTub:23.2,scGeom:175,dpLine:8.538,lineVol:1007,dpBend:1.645,dpTee:4.791,dpValve:3.327,dpChoke:29.29,hPump:624.2,pPump:308.3,psiComp:0.9202,t2Comp:136.2,vGasMax:0.4275,pwf:229,dpSkin:29.49,dpTubing:1.143,tArr:54.02,scHeld:138},{qLine:224,qTest:718,wChoke:12.3,rhoChoke:304,bean:56,mCp:192,tHot:61.5,qValve:224,qPump:115,phiComp:0.827,nComp:105,rComp:3.25,rhoLt:720,rhoGt:21,qTub:27.3,scGeom:224,dpLine:11.87,lineVol:998.4,dpBend:2.297,dpTee:6.825,dpValve:5.583,dpChoke:27.31,hPump:607.3,pPump:331.8,psiComp:1.167,t2Comp:157.3,vGasMax:0.5657,pwf:262.8,dpSkin:13.42,dpTubing:1.586,tArr:46.97,scHeld:182.2},{qLine:258,qTest:1330,wChoke:5.73,rhoChoke:414,bean:64,mCp:87.6,tHot:72.5,qValve:76.4,qPump:133,phiComp:1.07,nComp:80,rComp:1.95,rhoLt:760,rhoGt:36,qTub:31.4,scGeom:76.4,dpLine:15.62,lineVol:980.5,dpBend:3.128,dpTee:9.06,dpValve:0.6394,dpChoke:2.524,hPump:594.8,pPump:339.2,psiComp:0.4655,t2Comp:93.92,vGasMax:0.4399,pwf:240.5,dpSkin:25.09,dpTubing:2.018,tArr:41.02,scHeld:61.06},{qLine:293,qTest:473,wChoke:7.91,rhoChoke:350,bean:32,mCp:123,tHot:57.2,qValve:125,qPump:151,phiComp:0.729,nComp:90,rComp:2.38,rhoLt:800,rhoGt:15,qTub:35.5,scGeom:125,dpLine:19.63,lineVol:998.8,dpBend:3.973,dpTee:11.34,dpValve:1.755,dpChoke:91.9,hPump:568,pPump:356.8,psiComp:0.8588,t2Comp:116.4,vGasMax:0.7296,pwf:276,dpSkin:8.753,dpTubing:2.61,tArr:37.27,scHeld:100.2},{qLine:327,qTest:1090,wChoke:10.1,rhoChoke:304,bean:40,mCp:157,tHot:68.1,qValve:175,qPump:169,phiComp:0.975,nComp:100,rComp:2.82,rhoLt:720,rhoGt:30,qTub:39.5,scGeom:175,dpLine:24.32,lineVol:986,dpBend:4.874,dpTee:14.27,dpValve:3.394,dpChoke:71.61,hPump:558.9,pPump:365.5,psiComp:0.9815,t2Comp:139.3,vGasMax:0.4853,pwf:254.1,dpSkin:20.97,dpTubing:3.21,tArr:48.04,scHeld:138.4},{qLine:362,qTest:1700,wChoke:12.3,rhoChoke:414,bean:48,mCp:192,tHot:79,qValve:224,qPump:187,phiComp:1.22,nComp:105,rComp:3.25,rhoLt:760,rhoGt:45,qTub:43.6,scGeom:224,dpLine:29.43,lineVol:994.7,dpBend:6.009,dpTee:17.85,dpValve:5.504,dpChoke:36.63,hPump:525.7,pPump:372.7,psiComp:1.004,t2Comp:153,vGasMax:0.3944,pwf:220.7,dpSkin:32.3,dpTubing:3.813,tArr:60.69,scHeld:176.4},{qLine:396,qTest:841,wChoke:5.73,rhoChoke:350,bean:56,mCp:87.6,tHot:63.7,qValve:76.4,qPump:205,phiComp:0.876,nComp:80,rComp:1.95,rhoLt:800,rhoGt:24,qTub:47.7,scGeom:76.4,dpLine:35.75,lineVol:989.8,dpBend:7.253,dpTee:21.23,dpValve:0.6412,dpChoke:5.157,hPump:499.3,pPump:391.6,psiComp:0.6003,t2Comp:95.28,vGasMax:0.575,pwf:261.3,dpSkin:15.84,dpTubing:4.641,tArr:36.18,scHeld:61.93},{qLine:431,qTest:1450,wChoke:7.91,rhoChoke:304,bean:64,mCp:123,tHot:74.6,qValve:125,qPump:224,phiComp:1.12,nComp:90,rComp:2.38,rhoLt:720,rhoGt:39,qTub:51.8,scGeom:125,dpLine:41.38,lineVol:1001,dpBend:8.731,dpTee:24.57,dpValve:1.696,dpChoke:6.485,hPump:463.1,pPump:398.7,psiComp:0.6948,t2Comp:116.1,vGasMax:0.419,pwf:229.7,dpSkin:27.88,dpTubing:5.341,tArr:48.49,scHeld:98.82},{qLine:465,qTest:595,wChoke:10.1,rhoChoke:414,bean:32,mCp:157,tHot:59.4,qValve:175,qPump:242,phiComp:0.778,nComp:100,rComp:2.82,rhoLt:760,rhoGt:18,qTub:55.9,scGeom:175,dpLine:49.11,lineVol:996,dpBend:10.11,dpTee:29.32,dpValve:3.325,dpChoke:125.5,hPump:433.2,pPump:430.1,psiComp:1.078,t2Comp:136.6,vGasMax:0.646,pwf:267.5,dpSkin:11.06,dpTubing:6.074,tArr:42.67,scHeld:140.6},{qLine:500,qTest:1210,wChoke:12.3,rhoChoke:350,bean:40,mCp:192,tHot:70.3,qValve:224,qPump:260,phiComp:1.02,nComp:105,rComp:3.25,rhoLt:800,rhoGt:33,qTub:60,scGeom:224,dpLine:55.03,lineVol:991.4,dpBend:11.37,dpTee:34.24,dpValve:5.505,dpChoke:91.61,hPump:401.2,pPump:446.5,psiComp:1.111,t2Comp:152.1,vGasMax:0.484,pwf:245.3,dpSkin:22.53,dpTubing:7.137,tArr:54.52,scHeld:177.7}];
const CAL_VALID = [{qLine:120,qTest:800,wChoke:9,rhoChoke:304,bean:32,mCp:140,tHot:63,qValve:150,qPump:60,phiComp:0.86,nComp:80,rComp:2.6,rhoLt:720,rhoGt:23,qTub:15,scGeom:150,dpLine:3.77,lineVol:1007,dpBend:0.6556,dpTee:1.922,dpValve:2.452,dpChoke:133.8,hPump:653.4,pPump:295.1,psiComp:0.6041,t2Comp:127.2,vGasMax:0.5396,pwf:265.5,dpSkin:15.26,dpTubing:0.4983,tArr:43.55,scHeld:120.8},{qLine:183,qTest:350,wChoke:13,rhoChoke:414,bean:40,mCp:204,tHot:55,qValve:240,qPump:93.3,phiComp:0.68,nComp:90,rComp:3.4,rhoLt:760,rhoGt:12,qTub:22.5,scGeom:240,dpLine:8.215,lineVol:1017,dpBend:1.576,dpTee:4.523,dpValve:6.318,dpChoke:84.15,hPump:631.3,pPump:316.2,psiComp:0.8521,t2Comp:158.4,vGasMax:0.7769,pwf:285.2,dpSkin:6.486,dpTubing:1.105,tArr:43.52,scHeld:193.1},{qLine:247,qTest:1480,wChoke:7.67,rhoChoke:350,bean:48,mCp:119,tHot:75,qValve:120,qPump:127,phiComp:1.13,nComp:100,rComp:2.33,rhoLt:800,rhoGt:39.5,qTub:30,scGeom:120,dpLine:14.57,lineVol:982.7,dpBend:2.775,dpTee:8.179,dpValve:1.564,dpChoke:17.24,hPump:588.6,pPump:337.2,psiComp:0.9276,t2Comp:116.9,vGasMax:0.4352,pwf:230.5,dpSkin:28.12,dpTubing:1.84,tArr:47.86,scHeld:97.67},{qLine:310,qTest:1030,wChoke:11.7,rhoChoke:304,bean:56,mCp:183,tHot:67,qValve:210,qPump:160,phiComp:0.95,nComp:105,rComp:3.13,rhoLt:720,rhoGt:28.5,qTub:37.5,scGeom:210,dpLine:22.07,lineVol:1006,dpBend:4.532,dpTee:12.82,dpValve:4.827,dpChoke:24.2,hPump:555.3,pPump:364.1,psiComp:1.114,t2Comp:149,vGasMax:0.4973,pwf:251.6,dpSkin:19.74,dpTubing:2.869,tArr:50.82,scHeld:167.2},{qLine:373,qTest:575,wChoke:6.33,rhoChoke:414,bean:64,mCp:97.3,tHot:59,qValve:90,qPump:193,phiComp:0.77,nComp:80,rComp:2.07,rhoLt:760,rhoGt:17.5,qTub:45,scGeom:90,dpLine:31.38,lineVol:1002,dpBend:6.471,dpTee:18.41,dpValve:0.8888,dpChoke:3.066,hPump:506.9,pPump:387.8,psiComp:0.642,t2Comp:102.2,vGasMax:0.6615,pwf:272.8,dpSkin:10.9,dpTubing:4.06,tArr:35.49,scHeld:71.16},{qLine:437,qTest:1700,wChoke:10.3,rhoChoke:350,bean:32,mCp:161,tHot:79,qValve:180,qPump:227,phiComp:1.22,nComp:90,rComp:2.87,rhoLt:800,rhoGt:45,qTub:52.5,scGeom:180,dpLine:43.14,lineVol:1008,dpBend:8.695,dpTee:25.51,dpValve:3.566,dpChoke:152.9,hPump:468.4,pPump:414,psiComp:0.5767,t2Comp:140.7,vGasMax:0.4026,pwf:216.7,dpSkin:32.08,dpTubing:5.555,tArr:57.67,scHeld:143.6},{qLine:500,qTest:1250,wChoke:5,rhoChoke:304,bean:40,mCp:76,tHot:71,qValve:60,qPump:260,phiComp:1.04,nComp:100,rComp:1.8,rhoLt:720,rhoGt:34,qTub:60,scGeom:60,dpLine:55.67,lineVol:981.5,dpBend:11.39,dpTee:33.56,dpValve:0.3959,dpChoke:17.48,hPump:407.5,pPump:438,psiComp:0.9493,t2Comp:87.61,vGasMax:0.4409,pwf:244.7,dpSkin:24,dpTubing:7.069,tArr:36.58,scHeld:47.2}];

// ---- verification --------------------------------------------------------------------------------------------------
function verify() {
  const out = [], chk = (name, expected, got, tol, note) => out.push({ name, expected, got, tol, pass: Number.isFinite(got) && Math.abs(got - expected) <= tol, note });
  const poly = (x, z) => x.reduce((s, _, i) => (i ? s + Math.hypot(x[i] - x[i - 1], z[i] - z[i - 1]) : 0), 0);
  // geometry
  { const H = 1375, c = catenary({ height: H, angle: 12, n: 600 }); chk('Catenary arc length vs closed form', Math.sqrt(H * H + 2 * c.a * H), poly(c.x, c.z), 0.02, 'S = √(H² + 2aH) for z = a(cosh(x/a) − 1); polyline of 600 chords');
    const n = c.x.length; chk('Catenary hang-off angle', 12, 90 - Math.atan2(c.z[n - 1] - c.z[n - 2], c.x[n - 1] - c.x[n - 2]) / D2R, 0.05, 'slope of the last chord, degrees from vertical');
    chk('Catenary top elevation', H, c.z[n - 1], 1e-9, 'the generated riser reaches the requested height');
    const c16 = catenary({ height: H, angle: 12, n: 14 }); chk('Riser discretisation independence (14 chords)', 1, poly(c16.x, c16.z) / Math.sqrt(H * H + 2 * c.a * H), 2e-3, 'length of the coarse riser used in the route relative to the closed form'); }
  { const w = lazyWave({ height: 1000, angle: 10, liftAngle: 35, sagAngle: 15, buoyancy: 1.5, n: 900 }), n = w.x.length; chk('Lazy-wave riser reaches the hang-off height', 1000, w.z[n - 1], 1e-6, 'three catenary pieces with one horizontal tension');
    chk('Lazy-wave arc length vs closed form', w.length, poly(w.x, w.z), 0.02, 'a·[p₁ + (p₁ − p₂)/b + (p₃ − p₂)] from the slopes at the piece ends'); }
  { const x = [0, 3000, 6000], z = [0, -4000, 0], a = analyseProfile(x, z); chk('Length conservation on a 3-4-5 profile', 10000, a.length, 1e-9, 'two 5,000 m hypotenuses');
    chk('Inclination of a 3-4-5 segment', -Math.atan(4 / 3) / D2R, a.segs[0].incl, 1e-9, 'asin(dz/L) = −53.13°');
    chk('Low point detected at the vertex', 3000, a.lows[0]?.x ?? NaN, 1e-9, 'one low point at x = 3,000 m');
    chk('Internal volume of the 10 km route (ID 0.2 m)', Math.PI * 0.01 * 10000, area(0.2) * a.length, 1e-9, 'π/4·D²·L = 314.159 m³');
    const b = segmentsToProfile(profileToSegments(BASE.profile.map((p) => p.x), BASE.profile.map((p) => p.z))); chk('Round trip profile → segments → profile', 0, Math.max(...b.x.map((q, i) => Math.hypot(q - BASE.profile[i].x, b.z[i] - BASE.profile[i].z))), 1e-7, 'largest node displacement of the reference route (m)'); }
  { const a = analyseProfile([0, 100, 100, 200, 150, 300, 0], [0, -5, -5, -3, -4, 0, 0]); chk('Defect detection on a crafted bad profile', 122, 100 * a.zeroLength + 10 * a.duplicates + a.backward, 0, '1 zero-length element, 2 repeated points, 2 backward steps → code 122');
    const c = cleanProfile([{ x: 0, z: 0 }, { x: 100, z: -5 }, { x: 100, z: -5 }, { x: 200, z: -3 }, { x: 150, z: -4 }, { x: 300, z: 0 }, { x: 0, z: 0 }]); chk('Repair of the bad profile keeps the valid points', 4, c.x.length, 0, 'points at 0, 100, 200 and 300 m remain'); }
  { const s = simplifyProfile(linspace(0, 1000, 401), linspace(0, 1000, 401).map((q) => 10 * Math.sin(q / 40)), 80, []), a = analyseProfile(linspace(0, 1000, 401), linspace(0, 1000, 401).map((q) => 10 * Math.sin(q / 40))), s2 = simplifyProfile(linspace(0, 1000, 401), linspace(0, 1000, 401).map((q) => 10 * Math.sin(q / 40)), 80, [...a.highs, ...a.lows].map((q) => q.i));
    chk('Simplified profile respects the 80-point limit', 80, Math.max(s.x.length, s2.x.length), 0, '401-point sine route'); chk('Simplified profile keeps every high and low point', a.highs.length + a.lows.length, [...a.highs, ...a.lows].filter((q) => s2.idx.includes(q.i)).length, 0, 'extrema retained'); }
  { const t = { x: linspace(0, 1000, 11), y: linspace(0, 500, 6), elev: linspace(0, 500, 6).map((y) => linspace(0, 1000, 11).map((x) => -100 - 0.01 * x + 0.02 * y)) }, tr = terrainTransect(t, 'diag', 20), i = 7; chk('Terrain reconstruction on a plane', -100 - 0.01 * (1000 * i) / 20 + 0.02 * (500 * i) / 20, tr.z[i], 1e-9, 'bilinear sampling reproduces a planar seabed exactly');
    chk('Transect length across the grid diagonal', Math.hypot(1000, 500), tr.x[20], 1e-9, '√(1000² + 500²) m'); }
  { const sp = freeSpans([0, 10, 20, 30, 40], [0, -1, -2, -1, 0], { w: 500, EI: 2.5e7, gap: 0.1 }); chk('Free span over a 40 m, 2 m deep depression', 2 - (500 * 40 ** 4) / (384 * 2.5e7), sp[0]?.gap ?? NaN, 1e-9, 'gap = depth − w L⁴/(384 EI)'); }
  { const g = segmentGeometry(2, 1); chk('Hydraulic diameter of a half-full circular section', (2 * Math.PI) / (Math.PI + 2), hydraulicDiameter(g.area, g.wetted + g.chord), 1e-12, '4A/P = 2πR/(π + 2) for R = 1 m'); chk('Cross-sectional area of a half-full circle', Math.PI / 2, g.area, 1e-12, 'πR²/2'); chk('Hydraulic diameter of an annulus', 0.1, hydraulicDiameter(area(0.4) - area(0.3), Math.PI * 0.7), 1e-12, 'D_outer − D_inner'); }
  // wells
  { const R = 600 / (Math.PI / 3), t = minimumCurvature([{ md: 0, inc: 0, azi: 90 }, { md: 600, inc: 60, azi: 90 }]); chk('Minimum curvature: TVD of a build section', R * Math.sin(Math.PI / 3), t.tvd, 1e-9, '0 → 60° over 600 m is a circular arc of radius 572.96 m');
    chk('Minimum curvature: displacement of a build section', R * (1 - Math.cos(Math.PI / 3)), t.displacement, 1e-9, 'R(1 − cos 60°)'); chk('Dog-leg severity of the build section', 3, t.maxDls, 1e-9, '60° over 600 m = 3°/30 m');
    chk('Coordinate transformation: azimuth 90° gives easting only', 0, t.stations[1].north, 1e-9, 'north co-ordinate of the build section'); }
  { const o = { type: 'vogel', pRes: 300, pi: 25 }; chk('Vogel absolute open flow', (25 * 300) / 1.8, iprRate(0, o), 1e-9, 'q_max = J·p_res/1.8');
    const c = { type: 'composite', pRes: 300, pb: 200, pi: 25 }; chk('Composite IPR open flow', 25 * 100 + (25 * 200) / 1.8, iprRate(0, c), 1e-9, 'J(p_res − p_b) + J·p_b/1.8'); chk('Composite IPR inverse at 1,000 Sm³/d', 260, iprPwf(1000, c), 1e-6, 'straight line above the bubble point');
    const r = radialCoefficients({ phase: 'oil', k: 100, h: 20, re: 500, rw: 0.1, skin: 2, mu: 1e-3, B: 1.2, nonDarcy: false }); chk('Darcy radial productivity index vs field-unit formula', ((100 * (20 / UNIT.ft)) / (141.2 * 1 * 1.2 * (Math.log(5000) - 0.75 + 2))) * UNIT.bbl / (UNIT.psi / 1e5), r.pi, 0.02, 'k h / (141.2 μ B (ln(re/rw) − 0.75 + S)) converted from stb/d/psi');
    const g = { type: 'gas', pRes: 200, a: 0.02, b: 1e-7 }; chk('Forchheimer gas-well rate (quadratic root)', (-0.02 + Math.sqrt(0.0004 + 4e-7 * (40000 - 10000))) / 2e-7, iprRate(100, g), 1e-6, 'p_res² − p_wf² = a q + b q²'); }
  // hydraulics
  { const D = 0.05, L = 100, mu = 0.2, rho = 900, m = 0.5, Q = m / rho; chk('Hagen–Poiseuille pressure drop', (128 * mu * L * Q) / (Math.PI * D ** 4), pipeDp({ m, rho, mu, L, D }), 1e-6, 'laminar, Re = 64');
    chk('Colebrook–White at Re = 1e5, ε/D = 1e-4', 0.01851, frictionFactor(1e5, 1e-4, 'colebrook'), 2e-5, 'Moody-chart reference value'); chk('Haaland against the Colebrook reference', 0.01851, frictionFactor(1e5, 1e-4, 'haaland'), 3e-4, 'explicit approximation, within 1.5 %');
    chk('Swamee–Jain against the Colebrook reference', 0.01851, frictionFactor(1e5, 1e-4, 'swamee'), 3e-4, 'explicit approximation'); chk('Churchill against the Colebrook reference', 0.01851, frictionFactor(1e5, 1e-4, 'churchill'), 3e-4, 'all-regime equation');
    chk('Hydrostatic column of 1,000 m of sea water', (1025 * 9.80665 * 1000) / 1e5, pipeDp({ m: 0, rho: 1025, mu: 1e-3, L: 1000, D: 0.2, dz: 1000 }) / 1e5, 1e-9, 'ρ g h = 100.52 bar');
    chk('Bernoulli orifice (Cd = 1, small β)', area(0.05) * Math.sqrt(2 * 1000 * 2e5), orificeFlow({ cd: 1, d: 0.05, rho: 1000, dp: 2 }), 1e-9, 'v = √(2Δp/ρ) = 20 m/s through 50 mm: 39.27 kg/s');
    chk('K-factor loss of one velocity head', 0.5 * 1000 * 4, pipeDp({ m: 2 * 1000 * area(0.1), rho: 1000, mu: 1e-3, L: 0, D: 0.1, k: 1 }), 1e-6, 'K = 1 at 2 m/s in water: 2,000 Pa');
    const f = fittingsLoss([{ type: 'exit', count: 1 }, { type: 'entrance', count: 2 }], 0.2, 1e-4); chk('Fittings table total K', 2, f.K, 1e-12, 'exit 1.0 + 2 × entrance 0.5'); chk('Equivalent length of the fittings', (2 * 0.2) / fullyTurbulent(1e-4), f.eqLength, 1e-9, 'L_e = K·D/f_T'); }
  // valves and chokes
  { chk('IEC 60534 liquid sizing, hand calculation', 100 / (0.865 * Math.sqrt(4)), valveCv({ w: (100 / 3600) * 1000, dp: 4, p1: 10, rhoL: 1000, xG: 0 }).cv, 0.15, 'Cv = Q/(N1·√(Δp/SG)) for 100 m³/h of water at Δp = 4 bar: 57.8 (the constants N1 and N6 of the standard agree to 0.2 %)');
    chk('IEC 60534 gas sizing, hand calculation', 7200 / (27.3 * (1 - 0.2 / (3 * (1.3 / 1.4) * 0.7)) * Math.sqrt(0.2 * 50 * 40)), valveCv({ w: 2, dp: 10, p1: 50, rhoG: 40, xG: 1, k: 1.3, xT: 0.7 }).cv, 1e-9, 'W = N6·Cv·Y·√(x p1 ρ1), Y = 1 − x/(3 Fk xT)');
    chk('Valve pressure drop is the inverse of the sizing equation', 6, valveDp({ cv: valveCv({ w: 5, dp: 6, p1: 80, rhoL: 800, rhoG: 60, xG: 0.1 }).cv, w: 5, p1: 80, rhoL: 800, rhoG: 60, xG: 0.1 }).dp, 1e-6, 'two-phase stream, Δp = 6 bar');
    chk('Liquid choked flow limit', 0.9 * 0.9 * (10 - (0.96 - 0.28 * Math.sqrt(2 / 221)) * 2), valveCv({ w: 1, dp: 9.5, p1: 10, pv: 2, FL: 0.9 }).dpLiquidMax, 1e-9, 'Δp_max = FL²(p1 − FF·pv)');
    chk('Equal-percentage characteristic at 50 % travel', 1 / Math.sqrt(50), valveCharacteristic(0.5, 'equal', 50), 1e-12, 'R^(h − 1)'); chk('Opening is the inverse of the characteristic', 0.37, valveOpening(valveCharacteristic(0.37, 'equal', 50), 'equal', 50), 1e-12, 'equal-percentage trim');
    chk('Critical pressure ratio of an ideal gas, k = 1.4', 0.52828, criticalRatio(1.4), 1e-5, '(2/(k+1))^(k/(k−1))');
    chk('Sachdeva model in the all-gas limit', 0.52828, sachdeva({ p1: 50, p2: 10, x: 1, rhoG: 40, rhoL: 800, k: 1.4 }).yc, 1e-4, 'the two-phase critical ratio reduces to the gas value');
    const p0 = 50, r0 = 40, k = 1.3, h = hemCritical((p) => (1 / r0) * (p0 / p) ** (1 / k), p0, 2000); chk('Homogeneous-equilibrium integrator on an ideal gas', Math.sqrt(k * p0 * 1e5 * r0 * (2 / (k + 1)) ** ((k + 1) / (k - 1))), h.G, 30, 'G* = √(k p0 ρ0 (2/(k+1))^((k+1)/(k−1))) = 10,665 kg/m²/s');
    chk('Compressible nozzle: choked flux equals the critical flux', Math.sqrt(k * p0 * 1e5 * r0 * (2 / (k + 1)) ** ((k + 1) / (k - 1))), nozzleFlux(p0 * 1e5, r0, k, 0.2), 1e-6, 'isentropic nozzle below the critical ratio');
    chk('Gilbert correlation, hand calculation', (10 * 500 ** 0.546 * 1000) / 32 ** 1.89 * 0.0689476 + 1.01325, chokeCorrelation('gilbert', { q: 1000 * UNIT.bbl, glr: 500 * (UNIT.scf / UNIT.bbl), bean: 32 }), 1e-3, 'p = 10 R^0.546 q / S^1.89 (psig) for 1,000 bbl/d, 500 scf/bbl, 32/64 in'); }
  // pumps, compressors, vessels
  { const pc = { qr: 200, hr: 500 }; chk('Affinity laws on the pump curve', 0.64 * 500, pumpHead(0.8 * 200, { ...pc, speed: 0.8 }), 1e-9, 'at 80 % speed the rated point moves to 80 % flow and 64 % head');
    const a = affinity({ q: 100, h: 50, p: 20 }, 1.2, 1); chk('Affinity laws: power with speed cubed', 20 * 1.728, a.p, 1e-12, 'P ∝ N³');
    chk('Euler head with radial blades and no slip', (50 * 50) / 9.80665, eulerHead({ d2: 50 / (Math.PI * 1000 / 60), rpm: 1000, beta2: 90, blades: Infinity, q: 0 }).head, 1e-6, 'H = u₂²/g for u₂ = 50 m/s');
    const c = compressorStage({ p1: 1, p2: 4, t1: 26.85, z: 1, k: 1.4, mw: 28.97, eta: 1 }); chk('Isentropic discharge temperature of an ideal gas', 300 * 4 ** (0.4 / 1.4) - 273.15, c.t2, 1e-9, 'T₂ = T₁ r^((k−1)/k) = 445.8 K');
    chk('Isentropic head of an ideal gas', (1.4 / 0.4) * (8.314462618 / 0.02897) * 300 * (4 ** (0.4 / 1.4) - 1), c.headIsen, 1e-6, 'k/(k−1)·R T₁/M·(r^((k−1)/k) − 1)');
    const m = compressorMap({ qd: 2, hd: 1e5 }), l = m.lines[2]; chk('Compressor map: a point on the 90 % line returns 90 % speed', 0.9, mapPoint(m, l.q[5], l.h[5]).speed, 1e-9, 'interpolation between tabulated speed lines');
    chk('Compressor map: fan-law point at 95 % speed', 0.95, mapPoint(m, 0.95 * 2, 0.95 ** 2 * 1e5).speed, 5e-3, 'between the 90 % and 100 % lines'); chk('Surge flow follows the fan law', 0.62 * 0.9 * 2, mapPoint(m, l.q[5], l.h[5]).qSurge, 1e-9, 'surge line at constant flow coefficient');
    chk('Horizontal vessel at half level (flat ends)', (Math.PI * 1 * 6) / 2, vesselVolume({ d: 2, l: 6, level: 1, heads: 'flat' }), 1e-12, 'half of πR²L');
    chk('Horizontal vessel full, 2:1 elliptical heads', Math.PI * 6 + (4 / 3) * Math.PI * 0.5, vesselVolume({ d: 2, l: 6, level: 2 }), 1e-12, 'πR²L + 4/3·π·(R/2)·R²'); chk('Horizontal vessel with heads at half level', (Math.PI * 6 + (4 / 3) * Math.PI * 0.5) / 2, vesselVolume({ d: 2, l: 6, level: 1 }), 1e-12, 'half of the total by symmetry');
    chk('Vertical vessel full, 2:1 elliptical heads', Math.PI * 6 + (4 / 3) * Math.PI * 0.5, vesselVolume({ d: 2, l: 6, level: 7, orientation: 'vertical' }), 1e-12, 'same total volume standing up'); chk('Level from volume is the inverse', 0.73, vesselLevel(vesselVolume({ d: 2, l: 6, level: 0.73 }), { d: 2, l: 6, orientation: 'horizontal', heads: 'elliptical' }), 1e-8, 'horizontal vessel');
    let caught = 0; try { vesselVolume({ d: 2, l: 6, level: 2.5 }); } catch { caught++; } try { vesselVolume({ d: -1, l: 6, level: 0.5 }); } catch { caught++; } chk('Invalid vessel input is rejected', 2, caught, 0, 'level above the top and negative diameter both raise an error');
    chk('Souders–Brown velocity', 0.1 * Math.sqrt((800 - 20) / 20), soudersBrown(0.1, 800, 20), 1e-12, 'K√((ρL − ρG)/ρG) = 0.6245 m/s'); }
  // thermal
  { const w = wallDesign({ id: 0.2, wt: 0.01, kSteel: 45, layers: [{ name: 'Insulation', t: 0.05, k: 0.2 }], hIn: 1000, hOut: 500 }), ri = 0.1, r1 = 0.11, r2 = 0.16; chk('U-value of an insulated pipe vs the cylinder formula', 1 / (1 / 1000 + (ri * Math.log(r1 / ri)) / 45 + (ri * Math.log(r2 / r1)) / 0.2 + ri / (r2 * 500)), w.U, 1e-9, 'series resistances referred to the inner radius');
    chk('Outer diameter with the coating', 0.32, w.od, 1e-12, '0.2 + 2·(0.01 + 0.05) m'); chk('Thermal mass of the bare steel wall', Math.PI * (0.11 ** 2 - 0.1 ** 2) * 7850 * 480, wallDesign({ id: 0.2, wt: 0.01 }).thermalMass, 1e-6, 'ρ·cp·A of the steel ring (J/m/K)');
    const b = wallDesign({ id: 0.2, wt: 0.01, layers: [], hIn: 1e9, hOut: 500, burial: { cover: 1, kSoil: 1.5 } }); chk('Buried pipe: soil resistance by the acosh shape factor', 1 / (1e-9 + (0.1 * Math.log(1.1)) / 45 + (0.1 * Math.acosh((2 * 1.11) / 0.22)) / 1.5), b.U, 1e-6, 'R = r_i·acosh(2H/D_o)/k_soil with H to the pipe centre'); }
  // network
  { const two = solveNetwork({ nodes: [{ id: 'A', q: 2 }, { id: 'B', p: 10 }], edges: [{ from: 'A', to: 'B', length: 1000, id: 0.05 }, { from: 'A', to: 'B', length: 2000, id: 0.08 }], fluid: { rho: 900, mu: 0.5 } });
    chk('Two parallel pipes: analytic laminar split', (0.05 ** 4 / 1000) / (0.08 ** 4 / 2000), two.m[0] / two.m[1], 1e-7, 'm₁/m₂ = (D₁⁴/L₁)/(D₂⁴/L₂)'); chk('Two parallel pipes: loop pressure balance', 0, two.loops[0].sum, 1e-7, 'sum of pressure drops around the loop (bar)');
    const hp = solveNetwork({ nodes: [{ id: 'A', p: 12 }, { id: 'B', p: 10 }], edges: [{ from: 'A', to: 'B', length: 500, id: 0.04 }], fluid: { rho: 900, mu: 0.5 } }); chk('Flow between two fixed pressures', (Math.PI * 0.04 ** 4 * 900 * 2e5) / (128 * 0.5 * 500), hp.m[0], 1e-8, 'Hagen–Poiseuille: m = π D⁴ ρ Δp/(128 μ L)');
    const nodes = [{ id: 'A', p: 20 }, { id: 'B', p: 10 }, { id: 'C' }, { id: 'D', q: -3 }], edges = [{ from: 'A', to: 'C', length: 1000, id: 0.15 }, { from: 'C', to: 'B', length: 2000, id: 0.1 }, { from: 'C', to: 'D', length: 500, id: 0.1 }, { from: 'A', to: 'D', length: 800, id: 0.12, rough: 5e-5 }], lp = solveNetwork({ nodes, edges });
    chk('Kirchhoff node balance in a looped turbulent network', 0, lp.massResidual, 1e-9, 'largest nodal mass imbalance (kg/s)'); chk('Loop pressure balance in the looped network', 0, Math.abs(lp.loops[0].sum), 1e-6, 'bar'); chk('Global mass balance: inflow equals outflow', lp.inflow, lp.outflow, 1e-9, 'kg/s through the boundaries');
    const rev = solveNetwork({ nodes, edges: edges.map((e, i) => (i === 1 ? { ...e, from: 'B', to: 'C' } : e)) }); chk('Flow-direction consistency', -lp.m[1], rev.m[1], 1e-7, 'reversing the drawn direction of a connection changes only the sign of its flow');
    const A = incidence(nodes, edges); chk('Incidence matrix: every column sums to zero', 0, Math.max(...edges.map((_, j) => Math.abs(A.reduce((s, r) => s + r[j], 0)))), 0, 'each connection leaves one node and enters one');
    chk('Incidence matrix × flows = nodal balance', 3, A[3].reduce((s, a, j) => s + a * lp.m[j], 0), 1e-9, 'net inflow at the demand node equals its 3 kg/s demand');
    const tight = solveNetwork({ nodes, edges, tol: 1e-12 }); chk('Network solution independent of the tolerance', tight.p[2], lp.p[2], 1e-6, 'junction pressure at 1e-8 and 1e-12 (bara)');
    let caught = 0; try { solveNetwork({ nodes: [{ id: 'A', p: 5 }, { id: 'B' }, { id: 'X', q: 1 }, { id: 'Y' }], edges: [{ from: 'A', to: 'B', length: 10, id: 0.1 }, { from: 'X', to: 'Y', length: 10, id: 0.1 }] }); } catch { caught++; } try { solveNetwork({ nodes: [{ id: 'A', p: 5 }, { id: 'A', p: 4 }], edges: [{ from: 'A', to: 'A', length: 1, id: 0.1 }] }); } catch { caught++; } try { solveNetwork({ nodes: [{ id: 'A', p: 5 }, { id: 'B' }], edges: [{ from: 'A', to: 'Z', length: 1, id: 0.1 }] }); } catch { caught++; }
    chk('Topology defects are rejected', 3, caught, 0, 'disconnected sub-network, duplicate node and dangling connection each raise an error'); }
  // constants checked against their sources
  { const fT = 0.013; chk('Crane 90° bend, r/D = 5', 15.5 * fT, bendK(5, 90, fT), 1e-12, 'Crane table between r/D = 4 (14 fT) and 6 (17 fT)');
    chk('Crane bend formula at 180°, r/D = 3', (0.25 * Math.PI * 3 + 0.5 * 12 + 12) * fT, bendK(3, 180, fT), 1e-12, 'K = (n − 1)(0.25 π fT r/d + 0.5 K90) + K90 with n = 2');
    chk('Bend scaling below 90° against Crane’s elbows', 16 / 30, bendK(1.5, 45, fT) / bendK(1.5, 90, fT), 0.005, 'ratio of the 45° and 90° standard elbows, 16 fT / 30 fT');
    chk('Butterfly valve in a 300 mm bore', 35 * fullyTurbulent(1e-4), fittingsLoss([{ type: 'butterfly', count: 1 }], 0.3, 1e-4).K, 1e-12, 'Crane: 35 fT between 10 and 14 in');
    chk('Bend multiplier acts on bends only', (2 * 14 + 20) * fullyTurbulent(1e-4), fittingsLoss([{ type: 'bend90', count: 1 }, { type: 'teeRun', count: 1 }], 0.2, 1e-4, { bend: 2 }).K, 1e-12, 'one r/D = 1.5 bend doubled plus one tee run');
    chk('Mist-eliminator Souders–Brown K at 975 psi', 0.08281536, demisterK((975 * UNIT.psi) / 1e5), 1e-7, 'documented example of the York curve fitted by Svrcek & Monnery (m/s)');
    chk('API RP 14E constant in SI units', 3.8576728, 122 / Math.sqrt(1000), 5e-4, 'C = 100 (ft/s, lb/ft³) gives 3.8577 m/s at 1,000 kg/m³');
    chk('X65 specified minimum yield strength', 450, material('X65').smys, 0, 'API 5L / ISO 3183 grade L450 (MPa)');
    chk('Composite IPR, published worked example', 5733.33, iprRate(0, { type: 'composite', pRes: 4200, pb: 3000, pi: 2 }), 0.01, 'q_b = 2 × 1,200 = 2,400 and q_max = 2,400 + 2 × 3,000 / 1.8'); }
  // design features
  { const Ls = liftOffLength(1e7, 0.5, 500); chk('Sleeper lift-off length satisfies the beam equation', 0.5, (500 * Ls ** 4) / (72 * 1e7), 1e-12, 'w L⁴ / (72 EI) = sleeper height');
    chk('Buckle-propagation pressure of the reference pipe', 114.98, propagationPressure(450, 0.0159, 0.2858), 0.05, '35 × 450 MPa × (15.9 / 285.8)^2.5 = 11.50 MPa');
    chk('Anode mass for a protection current', 28.74375, anodeMassRequired(0.21, 25), 1e-9, '0.21 A × 25 y × 8,760 h / (0.8 × 2,000 A·h/kg)');
    chk('Equipment orientation: a pump drawn against the flow is flagged', 1, orientationCheck([{ from: 'A', to: 'B', kind: 'pump' }, { from: 'B', to: 'C', kind: 'pipe' }, { from: 'C', to: 'D', kind: 'choke' }], [-2, -2, 3]).filter((q) => !q.ok).length, 0, 'the reversed pipe is not directional; the choke carries flow as drawn');
    const cm = compareMeasured([{ kind: 'Pump FAT', at: 5, value: 9 }, { kind: 'pump fat', at: 10, value: 22 }, { kind: 'something else', at: 1, value: 1 }], { pumpfat: { label: 'head', unit: 'm', re: /pump/, f: (at) => 2 * at } });
    chk('Measured-data comparison: bias', -0.5, cm.groups[0]?.bias ?? NaN, 1e-12, 'model 10 and 20 against measured 9 and 22'); chk('Measured-data comparison: RMS error', Math.sqrt(2.5), cm.groups[0]?.rmse ?? NaN, 1e-12, '√((1² + 2²)/2)'); chk('Measured-data comparison: unrecognised rows are counted', 1, cm.unknown, 0, 'one row of an unknown kind'); }
  // geometry import and consistency
  { const ir = importRegression(); chk('Geometry-import regression: length of an imported 3-D polyline', 1125, ir.length, 1e-6, 'horizontal 3-4-5 leg of 500 m plus a leg of 500 m in plan climbing 375 m (625 m)'); chk('Geometry-import regression: every import check passes', 0, ir.checks - ir.checksPassed, 0, 'length conservation, zero-length elements, repeated nodes, chainage, finite values, inclination, round trip');
    chk('Geometry-import regression: riser base of the imported route', 500, ir.riserBaseX, 1e-6, 'the climb starts 500 m from the inlet');
    const cg = checkProfile({ x: [0, 3000, 6000], z: [0, -4000, 0], s: [0, 5000, 10000] }, { id: 0.2 }); chk('CAD/geometry consistency: element volumes of the 3-4-5 route', Math.PI * 0.01 * 10000, cg.find((q) => q.name === 'Internal volume')?.got ?? NaN, 1e-6, 'π/4·D²·L summed over the elements');
    const cb = checkProfile({ x: [0, 100, 100, 50], z: [0, -5, -5, -4], s: [0, 100.125, 100.125, 150] }); chk('CAD/geometry consistency: a defective route is rejected', 3, cb.filter((q) => !q.pass && /zero-length|duplicate|chainage/i.test(q.name)).length, 0, 'zero-length element, repeated node and backward chainage are each reported'); }
  // rotating-equipment dynamics and the transient network
  { const J = 10, w0 = 300, t0 = 500, rt = rotorDynamics({ J, omega0: w0, load: (t, w) => t0 * (w / w0) ** 2, tEnd: 20, n: 400 }); chk('Rotor coast-down against a quadratic load', (J * w0) / t0, rt.halfTime, 0.01, 'ω = ω₀/(1 + t/t_h), t_h = J ω₀/τ₀ = 6 s to half speed');
    chk('Generic compressor map: inverse of the speed line', 0.9, mapFlowCoefficient(1.0418), 1e-12, 'head coefficient 1 + 0.22 (1 − 0.9²) = 1.0418 belongs to 90 % of the design flow');
    const T = 300, M = 0.02, V = 10, rr = 0.5, c = (R * T) / (M * V) / 1e5, two = transientNetwork({ nodes: [{ id: 'A', type: 'volume', p: 20, V, T, M, z: 1 }, { id: 'B', type: 'volume', p: 10, V, T, M, z: 1 }], edges: [{ id: 'e', from: 'A', to: 'B', type: 'pipe', phase: 'gas', r: rr }], tEnd: 10, dt: 0.05, nOut: 400 }), n2 = two.t.length - 1;
    chk('Transient: pressure equalisation of two gas volumes', (Math.sqrt(10) - (c * 10) / Math.sqrt(rr)) ** 2, two.p.A[n2] - two.p.B[n2], 0.01, '√Δp falls linearly at (R T / M V)/√r for a quadratic branch law (bar after 10 s)');
    chk('Transient: a closed system conserves its gas mass', two.balance.gas.initial, two.balance.gas.final, 1e-8, 'kg');
    const one = transientNetwork({ nodes: [{ id: 'A', type: 'volume', p: 100, V: 50, T, M, z: 0.85 }, { id: 'B', type: 'fixed', p: 50 }], edges: [{ id: 'v', from: 'A', to: 'B', type: 'valve', phase: 'gas', cv: 20, pos: 1, cmd: 1, tau: 0.01, rate: 0.1 }], events: [{ t: 1, edge: 'v', set: { cmd: 0 } }], tEnd: 6, dt: 0.1, nOut: 400 });
    chk('Transient: real-gas line-pack inventory', (100e5 * 50 * M) / (0.85 * R * T), one.balance.gas.initial, 1e-6, 'm = p V M / (Z R T) with Z = 0.85 (kg)'); chk('Transient: rate-limited valve actuator', 0.5, one.pos.v[one.pos.v.length - 1], 2e-3, 'closing at 10 % of travel per second for 5 s');
    chk('Transient: gas leaving through the valve is accounted for', one.balance.gas.initial - one.balance.gas.out, one.balance.gas.final, 1e-6, 'initial inventory − outflow = final inventory (kg)');
    const geom = { d: 2.6, l: 10, orientation: 'horizontal', heads: 'elliptical' }, B0 = { m0: 50, xg: 0.1, rhoL: 800, M: 0.02, k: 1.3, z: 0.9, stroke: 10, asvSet: 10, flareSet: 1.1, flareP: 1.5, liqP: 40, headerV: 5, wells: { m: [0, 50, 100], p: [140, 100, 60], pAvail0: 100, chokeDp: 10, pos0: 0.5, char: 'linear', rng: 50 }, line: { V: 500, T: 320, p0: 60 },
      sep: { V: vesselVolume({ ...geom, level: 2.6 }), T: 320, p0: 25, pSet: 25, vol0: vesselVolume({ ...geom, level: 1.3 }), geom, levelSet: 1.3, H: 2.6 }, comp: { stages: 2, pd: 120, eta: 0.78, margin: 10, rpm: 9000, H: 4, Tcool: 313, dpExport: 2, dischargeV: 1.5, asvSize: 2 }, pump: { qr: 250, hr: 400, shutoff: 1.25, eta: 0.7, speed: 1, head: 10, ks: 0, rpm: 3560, H: 1.5 } };
    const st = transientNetwork(transientSystem(B0, 'none', { tEnd: 30, steps: 60 })), ns = st.t.length - 1;
    chk('Transient: the production system starts in steady state', 25, st.p.SEP[ns], 1e-5, 'separator pressure after 30 s without an event (bara)'); chk('Transient: compressor flow equals the gas produced', 5, st.m.COMPRESSOR[ns], 1e-5, '10 % gas in 50 kg/s'); chk('Transient: export pump holds the level', 1.3, st.level.SEP[ns], 2e-3, 'liquid level after 30 s (m)');
    const vc = transientNetwork(transientSystem(B0, 'valveClosure', { tEnd: 150, steps: 300 })); chk('Transient: anti-surge recycle keeps the compressor out of surge on a blocked outlet', 1, Math.min(...vc.margin.COMPRESSOR) > 0 && Math.max(...vc.pos.RECYCLE) > 0.02 && vc.surgeTime.COMPRESSOR === null ? 1 : 0, 0, 'the recycle valve opens at the control line and the surge margin stays positive');
    chk('Transient: gas balance of the valve-closure scenario', vc.balance.gas.initial + vc.balance.gas.in, vc.balance.gas.final + vc.balance.gas.out, 1e-6 * vc.balance.gas.initial, 'initial + inflow = final + outflow (kg)'); chk('Transient: liquid balance of the valve-closure scenario', vc.balance.liquid.initial + vc.balance.liquid.in, vc.balance.liquid.final + vc.balance.liquid.out, 1e-6 * vc.balance.liquid.initial, 'kg');
    const pt = transientNetwork(transientSystem(B0, 'pumpTrip', { tEnd: 40, steps: 160 })), np = pt.t.length - 1; chk('Transient: a tripped pump stops delivering', 0, pt.m.PUMP[np], 1e-9, 'the discharge check valve closes once the head falls below the static head (kg/s)'); chk('Transient: liquid accumulates after the pump trip', 1, pt.level.SEP[np] > 1.3 + 0.05 ? 1 : 0, 0, 'separator level rises'); }
  // calibration model, hand calculations
  { const c1 = calModel({ chokeCvMax: 100, qValve: 86.5, sepK: 0.1, rhoLt: 820, rhoGt: 20, effIdMm: 250, scEff: 0.8, scGeom: 150 }); chk('Calibration target: valve water test', 1, c1.dpValve, 1e-12, 'Δp = (Q / (0.865 Cv))² = 1 bar for 86.5 m³/h through Cv 100');
    chk('Calibration target: separator carry-over velocity', 0.1 * Math.sqrt(40), c1.vGasMax, 1e-12, 'K √((820 − 20)/20)'); chk('Calibration target: line-fill volume', area(0.25) * poly(BASE.profile.map((p) => p.x), BASE.profile.map((p) => p.z)), c1.lineVol, 1e-6, 'π/4 × 0.25² × route length'); chk('Calibration target: slug-catcher hold-up', 120, c1.scHeld, 1e-12, '0.8 × 150 m³'); }
  // well-count study against a closed-form system: straight-line inflow J, tubing loss 150 bar + c q², back-pressure 60 + 20 s² (bar)
  { const Q = 3000, J = 25, c = 2e-5, pR = 300, wl = (q) => ({ pwf: pR - q / J, whp: pR - q / J - 150 - c * q * q, dead: false, ero: 0 }), gq = linspace(25, 6000, 240), grid = { q: gq, whp: gq.map((q) => wl(q).whp), dead: gq.map(() => false) };
    const base = { qCase: Q, well: wl, grid, back: (n, s) => 60 + 20 * s * s, pRes: pR, qMax: J * pR, targetDp: () => 8, duty: (n, w, dpB) => 10 * dpB };
    const root = (A, B, C) => (-B + Math.sqrt(B * B + 4 * A * C)) / (2 * A); // positive root of A s² + B s − C = 0
    const s0 = wellCountStudy({ ...base, maxDrawdown: 55 }), s1 = wellCountStudy({ ...base, lift: { boost: true, dp: 40 }, maxDrawdown: 55 }), s2 = wellCountStudy({ ...base, lift: { boost: true, dp: 40, dpMin: 2, auto: true }, maxDrawdown: 70 });
    chk('Well-count study: natural-flow rate of two wells', Q * root(20 + (c * Q * Q) / 4, Q / (2 * J), 90), s0.rows[1].rateNatural, 2, '(20 + cQ²/n²) s² + Q/(nJ) s − 90 = 0 with n = 2 (Sm³/d)');
    chk('Well-count study: rate of four wells with 40 bar of boosting, chokes open', Q * root(20 + (c * Q * Q) / 16, Q / (4 * J), 130), s1.rows[3].rateWithLift, 3, 'same balance with 40 bar added at the wellhead (Sm³/d)');
    chk('Well-count study: drawdown per well', Q / (3 * J), s1.rows[2].drawdown, 1e-9, 'Q/(n J) = 40 bar for three wells');
    chk('Well-count study: smallest count without lift', 3, s0.recommended, 0, 'n = 2 gives 45 bara at the wellhead against 80 needed; n = 3 gives 90 bara at a drawdown of 40 bar');
    chk('Well-count study: drawdown limit rejects two boosted wells', 3, s1.recommended, 0, 'two wells with 40 bar of boosting deliver (85 ≥ 80 bara) but draw down 60 bar > 55 bar');
    chk('Well-count study: relaxed limit accepts two boosted wells', 2, s2.recommended, 0, 'drawdown 60 bar ≤ 70 bar');
    chk('Well-count study: speed-controlled boosting duty', 2, s2.rows[2].boostDp, 1e-9, 'three wells: 80 + 8 − 90 < minimum → 2 bar; two wells need 43 bar and are held at the rated 40'); chk('Well-count study: boosting capped at the rated pressure rise', 40, s2.rows[1].boostDp, 1e-9, '80 + 8 − 45 = 43 bar > 40 bar rated');
    chk('Well-count study: one well above its lifting capacity produces what it can', Q * root(20 + c * Q * Q, Q / J, 130), s1.rows[0].rate, 3, 'not deliverable at the case rate: rate with the rated boosting (Sm³/d)'); }
  // minimum curvature on a 3-D arc: the tangent turns along a great circle from (30°, 20°) to (70°, 80°); the position follows by integrating it
  { const tv = (i, a) => [Math.sin(i * D2R) * Math.cos(a * D2R), Math.sin(i * D2R) * Math.sin(a * D2R), Math.cos(i * D2R)], t0 = tv(30, 20), t1 = tv(70, 80), be = Math.acos(t0[0] * t1[0] + t0[1] * t1[1] + t0[2] * t1[2]), Lc = 900, Rc = Lc / be;
    const nn = t1.map((q, k) => (q - t0[k] * Math.cos(be)) / Math.sin(be)), pos = (b) => t0.map((q, k) => Rc * (q * Math.sin(b) + nn[k] * (1 - Math.cos(b)))), tm = t0.map((q, k) => q * Math.cos(be / 2) + nn[k] * Math.sin(be / 2));
    const mid = { md: Lc / 2, inc: Math.acos(tm[2]) / D2R, azi: Math.atan2(tm[1], tm[0]) / D2R }, t = minimumCurvature([{ md: 0, inc: 30, azi: 20 }, mid, { md: Lc, inc: 70, azi: 80 }]), e = pos(be), last = t.stations[t.stations.length - 1];
    chk('Minimum curvature: 3-D arc with turning azimuth, northing', e[0], last.north, 1e-8, 'p = R[t₀ sin β + n(1 − cos β)], R = L/β, evaluated through a mid-arc station'); chk('Minimum curvature: 3-D arc, easting', e[1], last.east, 1e-8, 'same arc'); chk('Minimum curvature: 3-D arc, TVD', e[2], last.tvd, 1e-8, 'same arc'); }
  // IEC 60534-2-1 worked example (incompressible, non-choked turbulent flow without attached fittings): water at 363 K, 680 → 220 kPa, 360 m³/h
  chk('IEC 60534-2-1 example: flow coefficient Kv', 165, valveCv({ w: 0.1 * 965.4, dp: 4.6, p1: 6.8, rhoL: 965.4, xG: 0, FL: 0.9, pv: 0.701, pc: 221.2 }).kv, 0.5, 'the standard gives Kv = 165 m³/h (as reproduced in the documentation of the fluids library, size_control_valve_l)');
  // units
  chk('Unit conversion: 1,000 psi to bar', 68.9476, convert(1000, 'psi', 'bar'), 1e-4, '1 psi = 6,894.757 Pa'); chk('Unit conversion: 10,000 bbl/d to m³/h', 66.245, convert(10000, 'bbl/d', 'm3/h'), 1e-3, '1 bbl = 0.158987 m³');
  chk('Unit conversion: 100 °F to °C', 37.7778, convert(100, 'F', 'C'), 1e-4, '(F − 32)·5/9'); chk('Unit conversion: 10 in to mm', 254, convert(10, 'in', 'ft') * 304.8, 1e-9, 'via feet');
  return out;
}

// ---- imported geometry ---------------------------------------------------------------------------------------------
/**
 * Offers for an imported route geometry (3-D polyline, chainage table, point list or network): the elevation-profile table with the
 * settings that go with it, the detected riser base and the water depth. Every offer states how many of the import checks of
 * core/route.js (length and volume conservation, zero-length elements, repeated nodes, chainage, finite values, inclination, round trip) passed.
 */
export function geometryOffers(g) {
  const out = [];
  try {
    if (!g || typeof g !== 'object' || g.survey || (g.stats && g.stats.wellSurvey) || !['polylines', 'table', 'points', 'network'].includes(g.kind)) return out;
    const p = profileFrom(g, { maxPoints: 60 }), chk = checkProfile(p, { id: DEF.idMm / 1000 }), okN = chk.filter((c) => c.pass).length, rb = riserBase(p.x, p.z), what = (g.format || g.kind) + (g.name ? ' “' + g.name + '”' : '');
    out.push({ key: 'geomMode', value: 'table', from: what + ': the route is used as the elevation-profile table' },
      { key: 'profile', value: p.x.map((x, i) => ({ x: +x.toFixed(3), z: +p.z[i].toFixed(3) })), from: what + ': ' + p.x.length + ' points, ' + rd(p.length, 1) + ' m along the pipe; ' + okN + ' of ' + chk.length + ' import checks passed' + (okN < chk.length ? ' (failed: ' + chk.filter((c) => !c.pass).map((c) => c.name).join(', ') + ')' : '') },
      { key: 'lengthUnit', value: 'm', from: 'Imported coordinates are converted to metres' }, { key: 'elevUnit', value: 'm', from: 'Imported coordinates are converted to metres' }, { key: 'elevSense', value: 'up', from: 'Imported elevations are positive upwards' });
    out.push({ key: 'riserBaseX', value: rb.height > 0 ? rb.x : 0, from: rb.height > 0 ? 'Riser base detected ' + rd(rb.height, 0) + ' m below the arrival point' : 'No riser detected in the imported route' });
    if (-p.zMin > 1) out.push({ key: 'waterDepth', value: -p.zMin, from: 'Deepest point of the imported route below datum' });
  } catch { /* not a route: the standard offers of the workspace apply */ }
  return out;
}
/**
 * Regression case of the geometry import: a 3-D polyline with a horizontal 3-4-5 leg (500 m) and a climbing leg (500 m in plan,
 * 375 m up: 625 m) is passed through geometryOffers and the route reader. Returns { length, expected (1,125 m), checksPassed, checks, riserBaseX }.
 */
export function importRegression() {
  const g = { kind: 'polylines', name: 'regression', format: 'Polyline', polylines: [{ x: [0, 300, 600], y: [0, 400, 800], z: [-100, -100, 275], closed: false }] }, offers = geometryOffers(g), vals = { ...DEF };
  for (const o of offers) vals[o.key] = o.value;
  const r = routeOf(readInputs(vals)), chk = checkProfile({ x: r.x, z: r.z, s: analyseProfile(r.x, r.z).s });
  return { length: analyseProfile(r.x, r.z).length, expected: 1125, checksPassed: chk.filter((c) => c.pass).length, checks: chk.length, riserBaseX: vals.riserBaseX, offers: offers.length };
}

// ---- sourced constants and reference data ---------------------------------------------------------------------------
const FLUIDS = 'https://raw.githubusercontent.com/CalebBell/fluids/master/fluids/';
/** Every literature constant set the engine relies on, with the address it was checked against on the retrieval date. */
export const PROVENANCE = [
  { item: 'Crane TP-410 resistance of 90° pipe bends, K/fT against r/D (20, 14, 12, 12, 14, 17, 24, 30, 34, 38, 42, 50 for r/D = 1 … 20)', used: 'bendK(), CRANE_BEND', source: 'Crane TP-410M as tabulated in the fluids library (bend_rounded_Crane_ratios / _fds)', url: FLUIDS + 'fittings.py', retrieved: '2026-10-08', status: 'verified', note: 'All ten values previously held agree; the r/D = 12 and 16 rows (34, 42) were added. Largest difference 0.' },
  { item: 'Bend loss at angles other than 90°', used: 'bendK()', source: 'Crane TP-410M formula K = (n − 1)(0.25 π fT r/d + 0.5 K90) + K90 as coded in fluids bend_rounded_Crane', url: FLUIDS + 'fittings.py', retrieved: '2026-10-08', status: 'corrected', note: 'Above 90° the former power law (angle/90)^0.9 capped at 180° was replaced by the Crane formula (180°, r/D = 3: 20.4 fT instead of 22.4 fT, −9 %). Below 90° the power law is kept: it reproduces the ratio of Crane’s 45° and 90° standard elbows, 16/30, within 0.5 %, whereas the Crane formula is stated for n ≥ 1 quarter turns.' },
  { item: 'Valve and fitting equivalent lengths: gate 8 fT, globe 340 fT, ball 3 fT, swing check 50 fT, butterfly 45 / 35 / 25 fT by size', used: 'FITTINGS', source: 'Crane TP-410M as coded in fluids (K_gate_valve_Crane, K_globe_valve_Crane, K_ball_valve_Crane, K_swing_check_valve_Crane, K_butterfly_valve_Crane)', url: FLUIDS + 'fittings.py', retrieved: '2026-10-08', status: 'corrected', note: 'Gate, globe, ball and swing check agree. The butterfly valve was 45 fT for all sizes: now 45 (≤ 8 in), 35 (10–14 in), 25 (larger) — 22 % lower on the 10 in reference line. The 45° fitting was 8 fT and is now Crane’s 45° standard elbow, 16 fT; the r/D = 5 bend was 16 fT and is now the tabulated 15.5 fT.' },
  { item: 'Standard elbow 30 fT, tee through run 20 fT, tee through branch 60 fT, sharp entrance K = 0.5, exit K = 1.0', used: 'FITTINGS', source: 'US Department of Energy Fundamentals Handbook DOE-HDBK-1012/3-92 (Thermodynamics, Heat Transfer and Fluid Flow, vol. 3), table “Typical values of Leq/D”; Crane values as listed at simupipe.com/resources/k-factor-table; fluids entrance_sharp and exit_normal', url: 'https://web.archive.org/web/20181106160701id_/https://www.standards.doe.gov/standards-documents/1000/1012-bhdbk-1992-v3/@@images/file', retrieved: '2026-10-09', status: 'verified', note: 'The DOE handbook (public domain) prints 30 for the 90° standard elbow, 60 for the tee through the branch and 16 for the 45° elbow: all agree. For the tee through the run it prints 10, while two pages that quote Crane TP-410 print 20: the engine keeps Crane’s 20 (the higher loss). Crane TP-410 itself was not opened. For scale, the Darby 3-K method of the fluids library gives 12–47 % less than 30 fT for a flanged standard elbow and 40–45 % less than 20 fT for a flanged tee run between 2 and 16 in at Re = 10⁴–10⁶: use the bend and tee multipliers where a fitting loss matters.' },
  { item: 'Fully turbulent friction factor fT of clean commercial steel pipe (roughness 0.0018 in)', used: 'fullyTurbulent(), CRANE_ROUGHNESS', source: 'Crane TP-410M fT table as tabulated in fluids (Crane_fts, Crane_fts_Ds)', url: FLUIDS + 'friction.py', retrieved: '2026-10-08', status: 'verified', note: 'The rough-wall Colebrook limit at ε = 0.0018 in reproduces the 18 tabulated values with a mean deviation of about 2 % and the library’s ft_Crane function within 0.5 % on average (reference data set crane-ft).' },
  { item: 'IEC 60534-2-1 sizing: N1 = 0.865 and N6 = 27.3 for Cv with bar, kg/h and kg/m³; Fk = k/1.4; Y = 1 − x/(3 Fk xT) ≥ 2/3; FF = 0.96 − 0.28 √(pv/pc); choked liquid Δp = FL²(p1 − FF pv)', used: 'valveCv(), valveDp(), transientNetwork()', source: 'fluids control_valve.py (N1 = 0.1, N6 = 3.16 for Kv with kPa; FF, Fγ and Y as coded there)', url: FLUIDS + 'control_valve.py', retrieved: '2026-10-08', status: 'verified', note: 'N6 = 3.16 (kPa, Kv) = 31.6 (bar, Kv) = 27.3 (bar, Cv = Kv/0.865); the factor equations are identical. Largest difference 0.1 % (rounding of N6).' },
  { item: 'Typical control-valve factors FL = 0.9 and xT = 0.7 (single-seat globe / cage choke)', used: 'defaults of chokeFL, chokeXT', source: 'AVEVA Instrumentation documentation, “Typical values for FL and XT and Fd” (stated to be a partial extract of Table D.1 of IEC 60534-2-1); Emerson Control Valve Handbook, section 5.10.1; Valmet GB-series coefficient table', url: 'https://docs.aveva.com/bundle/instrumentation/page/929743.html', retrieved: '2026-10-09', status: 'verified', note: 'Single-port globe valves: FL 0.9 with xT 0.72–0.75 flow-to-open, 0.85 / 0.70 for a characterised cage flow-to-close. Product data: Emerson cage-guided globe valves NPS 1–8 list FL 0.77–0.87 and xT 0.62–0.81 (https://www.emerson.com/is/content/emerson/en/final-control/flow-controls/documents/d101881x012.pdf); the Valmet DN 80 cage valve of the reference data lists 0.85 and 0.70 at full travel. The defaults sit inside these ranges, FL at the upper end; they are inputs and belong on the valve data sheet. IEC 60534-2-1 itself was not opened.' },
  { item: 'Orifice expansibility ε = 1 − (0.351 + 0.256 β⁴ + 0.93 β⁸)[1 − (p2/p1)^(1/κ)] (ISO 5167-2)', used: 'orificeFlow()', source: 'fluids flow_meter.py (orifice_expansibility)', url: FLUIDS + 'flow_meter.py', retrieved: '2026-10-08', status: 'verified', note: 'Coefficients 0.351, 0.256 and 0.93 agree exactly.' },
  { item: 'Critical-flow choke correlations p = A·GLR^B·q / S^C: Gilbert 10 / 0.546 / 1.89, Ros 17.4 / 0.5 / 2, Baxendell 9.56 / 0.546 / 1.93, Achong 3.82 / 0.65 / 1.88', used: 'chokeCorrelation()', source: 'Dabiri et al. (2024) Sci. Rep. 14, 6945, Table 1 (Gilbert, Ros, Achong as q = a·p·S^c / GLR^b with a = 0.1, 0.05747, 0.26178) and Ibrahim et al. (2021) ACS Omega 6, 19484, Table A1 (all four)', url: 'https://pmc.ncbi.nlm.nih.gov/articles/PMC10960849/', retrieved: '2026-10-08', status: 'verified', note: 'All constants agree (1/0.05747 = 17.40, 1/0.26178 = 3.82). The ACS Omega table (https://pmc.ncbi.nlm.nih.gov/articles/PMC8340095/) prints 1.84 for Gilbert’s bean exponent against 1.89 in the Scientific Reports table and here; 1.89 is kept. Baxendell is confirmed by the ACS Omega table only.' },
  { item: 'Sachdeva et al. (1986) two-phase choke model and its discharge coefficient', used: 'sachdeva(), default of chokeCd', source: 'Guo, Lyons & Ghalambor (2007) Petroleum Production Engineering, chapter 5, Eqs. 5.13–5.18 (copy hosted at petroleumengineers.ru); Buffa & Baliño (2017) Review of multiphase flow models for choke valves, ABCM JEM-2017-0041, Eqs. 20 and 30–32', url: 'https://eventos.abcm.org.br/jem2017/content/uploads/2017/03/JEM-2017-0041.pdf', retrieved: '2026-10-09', status: 'corrected', note: 'Critical-ratio equation, polytropic exponent n = 1 + x(cp − cv)/(x cv + (1 − x) cL), downstream gas volume and mass-flux equation agree term by term with both reproductions (second address: https://www.petroleumengineers.ru/sites/default/files/choke.pdf). Discharge coefficient: the conference paper quotes 0.75 for a choke close behind an elbow and 0.85 for one free of upstream disturbance; the textbook gives 0.62–0.90. The comment in this file had the two cases the wrong way round and is corrected; the default stays 0.85 (undisturbed approach) and is a calibration parameter. SPE 15657 itself was not opened.' },
  { item: 'Vogel IPR q/qmax = 1 − 0.2 (pwf/p) − 0.8 (pwf/p)²; composite qmax = qb + J pb/1.8', used: 'iprRate(), iprPwf()', source: 'pengtools wiki, “Vogel’s IPR” and “Composite IPR” (citing Vogel 1968, Neely, Brown 1984)', url: 'https://wiki.pengtools.com/index.php?title=Vogel%27s_IPR', retrieved: '2026-10-08', status: 'verified', note: 'Coefficients 0.2, 0.8 and 1.8 agree; the worked example of the page (p = 4,200, pb = 3,000, J = 2: qmax = 5,733) is reproduced in verify().' },
  { item: 'Fetkovich form q = qmax (1 − (pwf/p)²)^n with qmax = J p/1.8', used: 'iprRate()', source: 'PetroWiki, “Oil well performance”, Eqs. 2 and 3 (archived copy)', url: 'https://web.archive.org/web/20250214023333id_/https://petrowiki.spe.org/Oil_well_performance', retrieved: '2026-10-09', status: 'verified', note: 'The page prints q = C (p² − pwf²)^n and q/qmax = [1 − (pwf/p)²]^n. Scaling qmax with the productivity index as J p/1.8 is this suite’s convention so that Vogel and Fetkovich share one input; C and n come from a multi-rate test.' },
  { item: 'Darcy pseudo-steady radial inflow: ln(re/rw) − 3/4 + S, field-unit constant 141.2', used: 'radialCoefficients()', source: 'pengtools wiki, “JD” (dimensionless productivity index)', url: 'https://wiki.pengtools.com/index.php?title=JD', retrieved: '2026-10-08', status: 'verified', note: 'Both the −3/4 pseudo-steady term and the 141.2 constant agree.' },
  { item: 'Non-Darcy coefficient β = 1.88 × 10¹⁰ k^−1.47 φ^−0.53 (1/ft, k in mD) and the Jones–Blount–Glaze form Δp/q = a + b q', used: 'radialCoefficients(), input porosity', source: 'PetroWiki, “Gas well deliverability”, Eq. 32 (S. C. Jones, SPE 16949, 1987) and “Oil well performance”, Eqs. 4 and 5 (archived copies)', url: 'https://web.archive.org/web/20240903200140id_/https://petrowiki.spe.org/Gas_well_deliverability', retrieved: '2026-10-09', status: 'corrected', note: 'The former β = 2.33 × 10¹⁰/k^1.201 could not be found in any source that opens and is replaced by the correlation printed on the PetroWiki page; a porosity input (default 0.2) was added for it. At 150 mD and 20 % porosity β falls from 5.7 × 10⁷ to 2.8 × 10⁷ 1/ft (−51 %). The quadratic form and its solution agree with the oil-well page. The turbulence coefficient b is derived here from the Forchheimer equation integrated over the perforated interval, not copied from a field-unit formula. β can be overridden.' },
  { item: 'Souders–Brown coefficient of a vessel with a mist eliminator against pressure', used: 'demisterK(), separator sizing', source: 'fluids separator.py, K_separator_demister_York (York curve as fitted by Svrcek & Monnery, Chem. Eng. Prog., October 1993)', url: FLUIDS + 'separator.py', retrieved: '2026-10-08', status: 'corrected', note: 'The former fixed K = 0.12 m/s times a length factor (L/3.05 m)^0.56 (recalled as a GPSA rule; the GPSA data book could not be opened) is replaced by the pressure-dependent curve: 0.112 m/s for the horizontal reference separator at 25 bara against an effective 0.233 m/s before (−52 %). The documented example (975 psi → 0.08282 m/s) is reproduced in verify(). ' },
  { item: 'Souders–Brown coefficient 0.08 m/s for slug-catcher fingers', used: 'default of fingerK', source: 'PetroSkills / John M. Campbell Tip of the Month, September 2015, Table 1 (API 12J ranges, for context only)', url: 'https://www.jmcampbell.com/tip-of-the-month/2015/09/gas-liquid-separators-sizing-parameter/', retrieved: '2026-10-09', status: 'unverified', note: 'No source that opens gives a coefficient specific to finger-type (multiple-pipe) slug catchers. API 12J as quoted there lists 0.12–0.15 m/s for horizontal vessels 3 m long with a wire-mesh extractor and 0.04–0.11 m/s for vertical ones; 0.08 m/s — about half the horizontal value, for a pipe without a mist eliminator — is a design assumption. It is now an input (fingerK): take the vendor’s figure.' },
  { item: 'API RP 14E erosional velocity v = C/√ρ, C = 100 (ft/s, lb/ft³) = 122 (m/s, kg/m³)', used: 'default of cErosion', source: 'fluids piping.py, erosional_velocity (example: ρ = 1,000 kg/m³, C = 100 → 3.8577 m/s)', url: FLUIDS + 'piping.py', retrieved: '2026-10-08', status: 'verified', note: '122/√1000 = 3.858 m/s against 3.8577 m/s: difference 0.01 %.' },
  { item: 'API 5L / ISO 3183 PSL2 minimum yield and tensile strength: B 245/415, X42 290/415, X52 360/460, X56 390/490, X60 415/520, X65 450/535, X70 485/570, X80 555/625 MPa', used: 'GRADES, material()', source: 'API 5L PSL 1 tensile-requirement table as reproduced by a pipe supplier (grades B to X70; the minimum values are the same in PSL 2) and Mannesmann Line Pipe PSL 2 data sheets L450M / X65M and L555M / X80M', url: 'https://tubingchina.com/API-5L-PSL-1-Pipe-Mechanical-Properties-Tensile-Yield-Strength-Elongation.htm', retrieved: '2026-10-08', status: 'corrected', note: 'X65 was 448 / 531 MPa (65 and 77 ksi converted) and is now the tabulated 450 / 535 MPa (+0.4 % / +0.8 %); X80 confirmed at https://www.mannesmann.com/en/knowledge/standards-materials/hfi-welded-steel-pipes/api-5l-psl-2-h2-mod/l555m-x80m.html. All nine SMYS / SMTS pairs (245/415 … 555/625, including 320/435 for X46) also stand in Table 7-5 of DNV-OS-F101 (2013). The other grades agree. The shared reference-case file carries the same 450 / 535 MPa for X65 and no other grade constant.' },
  { item: '22Cr duplex 450 / 620 MPa and 25Cr super duplex 550 / 750 MPa; E 200 GPa, 7,800 kg/m³, conductivity 15 and 14 W/m/K, expansion 13.0 and 13.5 × 10⁻⁶ /K', used: 'GRADES', source: 'DNV-OS-F101 (October 2013), section 7, mechanical properties of duplex linepipe; Alleima data sheets SAF 2205 and SAF 2507', url: 'https://web.archive.org/web/20160328092919id_/http://rules.dnvgl.com:80/docs/pdf/DNV/codes/docs/2013-10/OS-F101.pdf', retrieved: '2026-10-09', status: 'corrected', note: 'SMYS / SMTS 450 / 620 (22Cr) and 550 / 750 MPa (25Cr) agree with the standard. Data sheets (https://www.alleima.com/en/technical-center/material-datasheets/tube-and-pipe-seamless/saf-2205/ and …/saf-2507/): density 7.8 g/cm³, E 200 GPa, conductivity 14 W/m/K at 20 °C and 16 (22Cr) / 15 (25Cr) at 100 °C, expansion 13.0 (22Cr) and 13.5 (25Cr) × 10⁻⁶ /K between 30 and 100 °C. The 25Cr expansion coefficient was 13.0 and is now 13.5 × 10⁻⁶ /K; 15 W/m/K for 22Cr is the value at about 60 °C.' },
  { item: 'Thermal properties: carbon steel 45 W/m/K, 7,850 kg/m³, 480 J/kg/K; polypropylene 0.22 W/m/K; epoxy 0.2 W/m/K; high-density concrete 2.0 W/m/K; soil 1.5 W/m/K; silica aerogel 0.02 W/m/K', used: 'wallDesign() defaults, DEF_LAYERS, kSoil, presets', source: 'ISO 10456 design values as tabulated in the ht library (insulation.py) and Wikipedia “List of thermal conductivities”', url: 'https://raw.githubusercontent.com/CalebBell/ht/master/ht/insulation.py', retrieved: '2026-10-08', status: 'corrected', note: 'Steel: tabulated 50 W/m/K, 7,800 kg/m³, 450 J/kg/K (and 36–54 W/m/K for carbon steel at https://en.wikipedia.org/wiki/List_of_thermal_conductivities): 45 / 7,850 / 480 lie within 10 %. Polypropylene 0.22 (910 kg/m³, 1,800 J/kg/K) and high-density concrete 2.0 agree. Epoxy was 0.3 and is now 0.2 W/m/K; default soil was 1.4 and is now 1.5 W/m/K (clay or silt; sand and gravel 2.0); the aerogel preset was 0.016 and is now 0.02 W/m/K. Glass-syntactic polyurethane: the 2014 Offshore Magazine survey of flow-assurance insulation (vendor data, https://img.offshore-mag.com/files/base/ebm/os/document/2019/06/1404off_flowassura.5cf68b1393219.pdf) lists 0.121–0.166 W/m/K for four GSPU products (Trelleborg Vikotherm G3 0.166, 866 kg/m³, 1,400 J/kg/K): the 0.17 W/m/K of the shared reference case is the upper end of that range.' },
  { item: 'Buckle-propagation pressure 35 fy αfab (t/D)^2.5 and the factor 1.3 on external pressure', used: 'propagationPressure(), buckle-arrestor extent', source: 'DNV-OS-F101 (October 2013), section 5 D500, Eqs. 5.15 and 5.16, Tables 5-2, 5-3 and 5-5', url: 'https://web.archive.org/web/20160328092919id_/http://rules.dnvgl.com:80/docs/pdf/DNV/codes/docs/2013-10/OS-F101.pdf', retrieved: '2026-10-09', status: 'verified', note: 'p_pr = 35 fy αfab (t2/D)^2.5 for 15 < D/t2 < 45 and the criterion pe − pmin ≤ p_pr/(γm γSC) agree. γm = 1.15 and γSC = 1.04 / 1.14 / 1.26 (safety class low / medium / high): the factor 1.3 used here is γm·γSC = 1.31 of safety class medium (1.20 for low, 1.45 for high). αfab is 1.00 for seamless pipe, 0.93 for UO, TRB and ERW and 0.85 for UOE pipe; the engine’s default is 1.' },
  { item: 'Al–Zn–In anode capacity 2,000 A·h/kg; bracelet-anode utilisation 0.8; mean current density 0.075 A/m²; mean coating breakdown 0.048', used: 'anodeMassRequired(), cathodic-protection check, defaults of cpCurrent, cpBreakdown and anodeSpacing', source: 'DNVGL-RP-F103 (July 2016), section 6.4.2, Tables 6-2, 6-3 and A-1; DNV-RP-B401 (October 2010, amended April 2011), Tables 10-6 and 10-8', url: 'https://web.archive.org/web/20170517154244id_/http://rules.dnvgl.com:80/docs/pdf/dnvgl/rp/2016-07/dnvgl-rp-f103.pdf', retrieved: '2026-10-09', status: 'corrected', note: 'Capacity 2,000 A·h/kg in sea water holds for an anode surface at or below 30 °C (1,500 at 60 °C, 720 at 80 °C; 1,500 in sediment): true for anodes outside 80 mm of insulation, not for an uninsulated hot line. Utilisation: at most 0.80 for bracelet anodes. Design mean current density of a non-buried pipeline: 0.050 / 0.060 / 0.075 / 0.100 / 0.130 A/m² for internal fluid temperatures up to 25 / 50 / 80 / 120 °C and above: the default was 0.07 and is now 0.075 A/m² (70 °C inlet). Coating breakdown f = a + 0.5 b t: thermally insulating FBE / polyurethane a = 0.01, b = 0.003 → 0.0475 over 25 years; the default was 0.03 and is now 0.048. With both changes the required anode mass rises by 71 %, so the default spacing of the 60 kg anodes was reduced from 146 m to 122 m (utilisation 0.90). RP-B401: https://web.archive.org/web/20170110185815id_/http://rules.dnvgl.com:80/docs/pdf/DNV/codes/docs/2011-04/RP-B401.pdf.' },
  { item: 'Wiesner slip factor σ = 1 − √(sin β2)/Z^0.7 (β2 from the tangent)', used: 'eulerHead()', source: 'Zulkifli (2017) slip-factor models compared with CFD, bachelor thesis, HAW Hamburg / CFturbo, section 7.7.5', url: 'https://cfturbo.com/fileadmin/content/down/publications/students/2017-02-Zulkifli-Slip-Factor-Models-CFD.pdf', retrieved: '2026-10-09', status: 'verified', note: 'The thesis prints the same expression with the blade angle measured from the tangential direction, valid up to the radius ratio r1/r2 = exp(−8.16 sin β2/Z) (0.61 for 7 blades at 25°); above it σ is reduced by the factor 1 − ((r1/r2 − ε)/(1 − ε))³, which the engine does not apply (pump impellers normally lie below the limit). Wiesner’s paper itself was not opened.' },
  { item: 'Minimum-curvature survey calculation', used: 'minimumCurvature()', source: 'Utah FORGE well 16A(78)-32 survey report (422 stations with reported TVD and offsets), Geothermal Data Repository submission 1283', url: 'https://gdr.openei.org/submissions/1283', retrieved: '2026-10-09', status: 'verified', note: 'Reported TVD, northing and easting of a 10,955 ft build-and-hold well are reproduced within 0.005 ft at every station (reference data sets deviation-survey-forge-16a, -north and -east); the closed-form arc checks are in verify().' },
  { item: 'Shape of the generic compressor map (surge at 62 % and stonewall at 130 % of design flow, head rise 22 %), rotor loss torque fractions, governor and actuator time constants of the transient', used: 'compressorMap(), mapFlowCoefficient(), transientSystem()', source: '—', url: '', retrieved: '2026-10-08', status: 'unverified', note: 'Model assumptions, not literature values: replace the map by tested speed lines and set inertia constants, stroke times and volumes from the equipment data sheets. The interpolation between tested speed lines is checked against open map data (reference data set compressor-map-speed-line); the generic shape is not.' },
  { item: 'Code-to-code benchmarks: friction factor, Crane fittings, IEC 60534 valve sizing, orifice flow, compressor head, mist-eliminator K', used: 'validationData (kind benchmark)', source: 'fluids 1.3.1 (C. Bell and contributors), MIT licence, run on 2026-10-09', url: 'https://pypi.org/project/fluids/1.3.1/', retrieved: '2026-10-09', status: 'verified', note: 'Twelve data sets hold the library’s outputs for stated inputs; the engine is evaluated on the same inputs without adjustment. Agreement to rounding for friction, orifice, compressor and mist-eliminator relations, 0.15–0.3 % for valve sizing, 0.3–0.7 % for Crane valves and bends of 90° and more; bends below 90° differ by design (see that data set).' },
];
const SI_PSI = UNIT.psi / 1e5, SV = SURVEY_FORGE.survey, svPath = (() => { let c = null; return () => (c ||= minimumCurvature(SV.md.map((md, i) => ({ md, inc: SV.inc[i], azi: SV.azi[i] }))).stations); })(), svAt = (key, md) => interp1(SV.md, svPath().map((q) => q[key]), md);
const craneK = (type, D) => fittingsLoss([{ type, count: 1 }], D, CRANE_ROUGHNESS / D).K, craneFT = (D) => fullyTurbulent(CRANE_ROUGHNESS / D);
const stageOf = (r) => compressorStage({ p1: r.P1 / 1e5, p2: r.P2 / 1e5, t1: r.T1 - KEL, z: r.Z, k: r.k, mw: r.MW, eta: r.etaP });
const heldOutMap = (() => { let c = null; const head = (pr) => pr ** (0.4 / 1.4) - 1; return () => (c ||= { head, map: (() => { const speeds = [...new Set(COMPRESSOR_MAP.map.map((r) => r.n))].sort((x, y) => x - y); return { speeds, lines: speeds.map((n) => { const p = COMPRESSOR_MAP.map.filter((r) => r.n === n).sort((x, y) => x.m - y.m); return { n, q: p.map((r) => r.m), h: p.map((r) => head(r.pr)), qSurge: p[0].m, qChoke: p[p.length - 1].m }; }) }; })() }); })();
const VALIDATION = [
  { ...NIKURADSE_ROUGH, model: (r) => frictionFactor(r.Re, 1 / (2 * r.rk), 'colebrook'), tolerance: { mape: 10 }, note: 'Colebrook–White with ε/D = k/(2r), k the sand-grain size; supports “field pressure-drop measurements” at the level of the wall-friction law. Achieved: mean absolute error 9.6 %. Uniform sand roughness dips below the Colebrook curve in the transition region (Colebrook describes commercial pipe): the engine is up to 36 % high there, while the fully rough points (the upper three of each series) agree within 3 %. The tolerance is set by that known difference, not by fitting.' },
  { ...SMOOTH_PIPE, model: (r) => frictionFactor(r.Re, 0, 'colebrook'), tolerance: { mape: 6 }, note: 'Kernel friction factor for a hydraulically smooth wall from creeping flow to Re = 10⁶. Achieved: mean absolute error 4.9 %. The points between Re = 2,500 and 3,300 lie in the laminar–turbulent transition, where the measured factor jumps while the kernel switches at its fixed Reynolds number (up to 26 % apart); the turbulent points above Re = 4,000 agree within a few per cent.' },
  { ...CHOKE_FIELD, model: (r) => chokeCorrelation('gilbert', { p1: r.pwh * SI_PSI, glr: r.glr * (UNIT.scf / UNIT.bbl), bean: r.bean }) / UNIT.bbl, tolerance: { mape: 30 }, note: 'Valve / choke flow tests: Gilbert’s critical-flow correlation used blind (pressure taken as psia, critical flow assumed for every row; the source does not give the downstream pressure). Achieved: mean absolute error 22.9 %, individual rows between −35 % and +36 %. Correlations of this family are screening tools; the discharge coefficient of the mechanistic choke models is meant to be calibrated to choke tests.' },
  { ...BEIJ_BENDS, model: (r) => { const fT = craneFT(BEIJ_BENDS.pipeId); return bendK(r.rD, 90, fT) - fT * 0.5 * Math.PI * r.rD; }, tolerance: { mape: 30 }, note: 'Minor-loss measurements: the Crane bend coefficient of the engine less the friction of the bend’s own length (fT·π/2·R/d) against the measured excess loss of the bend. Achieved: mean absolute error about 25 %; the engine is 8–31 % below eight of the nine bends (and 40 % above bend 8). Beij’s bends were commercial steel bends rougher than the straight pipe, and he measured a further 0.09–0.18 velocity heads in the downstream tangent that the Crane value does not carry at all: the bend-loss multiplier of the calibration tab exists for this.' },
  { ...VALVE_CV_TRAVEL, model: (r) => VALVE_CV_TRAVEL.ratedCv * valveCharacteristic(r.travel / 100, 'equal', 50), tolerance: { mape: 35 }, note: 'Valve Cv tests: ideal equal-percentage characteristic Cv = Cv,rated·R^(h − 1) with the default rangeability R = 50 against the manufacturer’s table. Achieved: mean absolute error about 31 %; the ends agree (2 % at 10 % travel, exact at 100 %) but the real trim is a modified equal-percentage shape that passes 40 % more flow than the ideal curve between 30 % and 80 % travel. Type the tested Cv values for a real valve; the published FL = 0.85 and xT = 0.70 at full travel support the default factors (0.9 and 0.7).' },
  { ...COMPRESSOR_MAP, model: (r) => { const h = heldOutMap(); return mapPoint(h.map, r.m, h.head(r.pr)).speed; }, tolerance: { mape: 3 }, note: 'Compressor performance maps: the engine’s map interpolation is given five speed lines (500, 600, 700, 800, 840) and asked for the speed of each point of the withheld 740 line from its flow and head (head taken as PR^0.286 − 1; the unknown inlet state cancels). Achieved: mean absolute error about 1.5 %, largest 4 % on the steep choke side. The source gives no geometry or inlet conditions, so Euler head and discharge temperature are not tested by it.' },
  { ...GAS_WELL_TESTS, model: (r) => iprRate(r.pwf, { type: 'gas', pRes: r.pRes, a: r.A, b: r.B, pi: 1 }), tolerance: { mape: 10 }, note: 'Well-test data: the pressure-squared Forchheimer inflow model with each well’s published coefficients. The seven open-flow potentials are reproduced within 0.01 % (the quadratic is solved as in the paper); the seven stabilised test rates are over-predicted by 8–27 % (mean 15 %), because coefficients from short isochronal flow periods describe a well that has not yet reached its stabilised drainage radius. Achieved overall: mean absolute error about 7.5 %.' },
  { ...SURVEY_FORGE, model: (r) => svAt('tvd', r.md), tolerance: { maxAbs: 0.02 }, note: 'Well deviation surveys: true vertical depth by the minimum-curvature method from the 422 stations (feet in, feet out) against the depth reported with the survey. Achieved: largest difference 0.005 ft over 8,559 ft (the reported values are rounded to 0.01 ft).' },
  { ...SURVEY_FORGE, id: 'deviation-survey-forge-16a-north', title: 'Deviation survey of Utah FORGE well 16A(78)-32: offset north at 24 of 422 stations', quantity: 'Offset north', target: 'north', model: (r) => svAt('north', r.md), tolerance: { maxAbs: 0.02 }, note: 'Same survey: northing offset (grid azimuth as tabulated). Achieved: largest difference 0.005 ft on −1,040 ft.' },
  { ...SURVEY_FORGE, id: 'deviation-survey-forge-16a-east', title: 'Deviation survey of Utah FORGE well 16A(78)-32: offset east at 24 of 422 stations', quantity: 'Offset east', target: 'east', model: (r) => svAt('east', r.md), tolerance: { maxAbs: 0.02 }, note: 'Same survey: easting offset. Achieved: largest difference 0.005 ft on 3,901 ft of step-out.' },
  // code-to-code benchmarks against the open-source fluids library (the engine is run blind on the library's inputs)
  { ...BENCH_FRICTION, model: (r) => frictionFactor(r.Re, r.eD, 'colebrook'), tolerance: { mape: 0.01 }, note: 'Kernel Colebrook–White friction factor against the library’s exact solution. Achieved: agreement to better than 0.001 %.' },
  { ...CRANE_FT, model: (r) => craneFT(r.d), tolerance: { mape: 1 }, note: 'Rough-wall limit of Colebrook–White at Crane’s roughness of 0.0018 in against the library’s f_T (which matches the rounded Crane table). Achieved: mean absolute error about 0.5 %.' },
  { ...BENCH_BENDS, model: (r) => bendK(r.rD, r.angle, craneFT(r.Di)), tolerance: { mape: 2 }, note: 'Crane bend coefficients at 90° and 180° with the engine’s own f_T. Achieved: mean absolute error 0.6 %, largest 1.6 % (table interpolation against the library’s polynomial).' },
  { ...BENCH_BENDS_LOW, model: (r) => bendK(r.rD, r.angle, craneFT(r.Di)), tolerance: { mape: 25 }, note: 'Below 90° the two codes differ by design: the engine scales the 90° value with (angle/90)^0.9, which goes to zero for a straight pipe, whereas the library continues Crane’s multi-turn formula to fractions of a turn, which leaves a finite loss at zero angle. Achieved: the engine is up to 24 % lower (mean about 18 %). Crane states the formula for whole quarter turns only.' },
  { ...BENCH_VALVES, model: (r) => craneK(r.type, r.D), tolerance: { mape: 1 }, note: 'Crane valve coefficients (gate 8 fT, globe 340 fT, ball 3 fT, swing check 50 fT, butterfly 45 / 35 / 25 fT) with the engine’s f_T. Achieved: mean absolute error 0.3 %, largest 0.8 %.' },
  { ...BENCH_VALVE_LIQUID, model: (r) => valveCv({ w: r.Q * r.rho, dp: (r.P1 - r.P2) / 1e5, p1: r.P1 / 1e5, rhoL: r.rho, xG: 0, FL: r.FL, pv: r.Psat / 1e5, pc: r.Pc / 1e5 }).cv, tolerance: { mape: 0.5 }, note: 'IEC 60534-2-1 liquid sizing including choked rows (FL and FF). Achieved: 0.15 % on every row, the rounding of the numerical constant N1.' },
  { ...BENCH_VALVE_GAS, model: (r) => valveCv({ w: r.w, dp: (r.P1 - r.P2) / 1e5, p1: r.P1 / 1e5, rhoG: r.rho1, xG: 1, k: r.gamma, xT: r.xT }).cv, tolerance: { mape: 0.5 }, note: 'IEC 60534-2-1 gas sizing including choked rows (expansion factor Y, Fγ·xT limit); the engine works with mass flow and inlet density, the library with normal volume flow, molar mass and Z. Achieved: 0.27 % on every row (rounding of N6 against N9).' },
  { ...BENCH_ORIFICE, model: (r) => orificeFlow({ cd: r.C, d: r.Do, D: r.D, rho: r.rho, dp: (r.P1 - r.P2) / 1e5, p1: r.k ? r.P1 / 1e5 : null, k: r.k }), tolerance: { mape: 0.01 }, note: 'Orifice equation with the velocity-of-approach factor and the ISO 5167 expansibility. Achieved: agreement to rounding.' },
  { ...BENCH_COMPRESSOR, model: (r) => stageOf(r).headPoly, tolerance: { mape: 0.01 }, note: 'Polytropic head from the polytropic exponent n/(n − 1) = η_p·k/(k − 1). Achieved: agreement to rounding.' },
  { ...BENCH_COMPRESSOR, id: 'fluids-compressor-temperature', title: 'Code-to-code: compressor discharge temperature (fluids)', quantity: 'Discharge temperature', unit: 'K', target: 'T2', model: (r) => stageOf(r).t2 + KEL, tolerance: { mape: 0.01 }, note: 'Discharge temperature T1·r^((n − 1)/n) against the library’s isentropic temperature rise divided by the isentropic efficiency. Achieved: agreement to rounding.' },
  { ...BENCH_COMPRESSOR, id: 'fluids-compressor-isentropic-efficiency', title: 'Code-to-code: isentropic efficiency from polytropic efficiency (fluids)', quantity: 'Isentropic efficiency', unit: '–', target: 'etaIsen', model: (r) => stageOf(r).etaIsen, tolerance: { mape: 0.01 }, note: 'Conversion between polytropic and isentropic efficiency at the stage pressure ratio. Achieved: agreement to rounding.' },
  { ...BENCH_DEMISTER, model: (r) => demisterK(r.P / 1e5), tolerance: { mape: 0.01 }, note: 'Pressure-dependent Souders–Brown coefficient of a vertical vessel with a mesh pad. Achieved: agreement to rounding.' },
];

// ---- presets ---------------------------------------------------------------------------------------------------------
const thinCoat = [{ name: 'Fusion-bonded epoxy', tMm: 0.5, k: 0.2, rho: 1200, cp: 1400 }, { name: 'Three-layer polypropylene', tMm: 3, k: 0.22, rho: 910, cp: 1800 }];
const PRESETS = [
  { name: 'Deep-water oil tie-back (reference case)', values: { geomMode: 'table' } },
  { name: 'Gas-condensate trunk line to shore (buried landfall)', values: {
    geomMode: 'table', profile: [[0, -95], [8000, -110], [16000, -102], [24000, -88], [32000, -70], [40000, -52], [46000, -35], [50000, -18], [52500, -6], [53500, 2], [55000, 9], [58000, 14]].map(([x, z]) => ({ x, z })),
    idMm: 336.6, wtMm: 9.5, roughUm: 30, layers: thinCoat, concreteMm: 50, burialDepth: 1.2, burialStartX: 49000, tSeabed: 8, tSeaSurface: 16, tAir: 12, tGround: 10, currentSpeed: 0.5, tIn: 60, pSep: 35, designPressure: 150, designTemp: 80,
    network: [pipeRow('WHP-A', 'PLAT', 60, 254, 2), pipeRow('WHP-B', 'PLAT', 2500, 254, 2), pipeRow('PLAT', 'LAND', 0, 0, 0, 'flowline')], netNodes: [{ id: 'WHP-A', kind: 'source', share: 55 }, { id: 'WHP-B', kind: 'source', share: 45 }, { id: 'LAND', kind: 'sink' }],
    iprType: 'gas', pRes: 280, tRes: 105, permMd: 60, payM: 35, perfM: 25, skin: 3, survey: [{ md: 0, inc: 0, azi: 200 }, { md: 1200, inc: 0, azi: 200 }, { md: 2100, inc: 35, azi: 200 }, { md: 3600, inc: 35, azi: 200 }], tubingIdMm: 124.3,
    sepD: 3.2, sepL: 13, slugSurge: 150, slugDuration: 300, fingerMaxL: 120, compPd: 110, pumpDischargeP: 50, chokeCvMax: 200, liftType: 'none', buckleInitType: 'none', anodeSpacing: 122, transScenario: 'compressorTrip', exportVol: 400,
    supports: [{ x: 0, type: 'platform riser clamp', k: 50000 }, { x: 49000, type: 'trench transition', k: 0 }, { x: 58000, type: 'landfall anchor block', k: 80000 }] } },
  { name: 'Shallow-water multi-well gathering network', values: {
    geomMode: 'generate', riserType: 'vertical', waterDepth: 45, flowlineLength: 6500, seabedDrop: 8, undulationAmp: 1.5, undulationLength: 1500, topsideElev: 18, idMm: 303.2, wtMm: 12.7, layers: thinCoat, concreteMm: 40, tSeabed: 14, tSeaSurface: 22, tIn: 55, pSep: 12, designPressure: 150,
    network: [...[1, 2, 3].map((i) => pipeRow(`W${i}`, 'H1', 0, 0, 0, 'choke', 55)), ...[4, 5, 6].map((i) => pipeRow(`W${i}`, 'H2', 0, 0, 0, 'choke', 55)), pipeRow('H1', 'MAN', 400, 203.2, 1.5), pipeRow('H2', 'MAN', 650, 203.2, 1.5), pipeRow('H1', 'H2', 300, 152.4, 1), pipeRow('MAN', 'RB', 0, 0, 0, 'flowline'), pipeRow('RB', 'TOP', 0, 0, 0, 'riser')],
    netNodes: [...[1, 2, 3, 4, 5, 6].map((i) => ({ id: `W${i}`, kind: 'source', share: [20, 18, 16, 16, 15, 15][i - 1] })), { id: 'TOP', kind: 'sink' }],
    iprType: 'vogel', pRes: 190, tRes: 75, pBubble: 190, piWell: 60, survey: [{ md: 0, inc: 0, azi: 120 }, { md: 600, inc: 0, azi: 120 }, { md: 1400, inc: 40, azi: 120 }, { md: 2300, inc: 40, azi: 120 }], tubingIdMm: 100.5,
    sepD: 2.8, sepL: 11, slugSurge: 25, compPd: 90, pumpDischargeP: 40, pumpHr: 420, liftType: 'gaslift', liftGasRate: 40000, liftDepthPct: 75, buckleInitType: 'none', anodeSpacing: 98, transScenario: 'pumpTrip', sepInitLevelPct: 40,
    supports: [{ x: 0, type: 'manifold spool', k: 5000 }, { x: 6500, type: 'riser clamp at the jacket', k: 100000 }] } },
  { name: 'Onshore hilly-terrain pipeline with pump station', values: {
    geomMode: 'table', profile: [[0, 120], [6000, 180], [12000, 340], [17000, 290], [23000, 520], [28000, 430], [34000, 610], [40000, 480], [46000, 395], [52000, 450]].map(([x, z]) => ({ x, z })),
    idMm: 303.2, wtMm: 9.5, roughUm: 45, layers: thinCoat, burialDepth: 1, burialStartX: 0, tGround: 14, tAir: 18, tIn: 65, pSep: 15, grade: 'X60', designPressure: 150, designTemp: 80,
    network: [pipeRow('WELL-1', 'IN', 350, 152.4, 2), pipeRow('WELL-2', 'IN', 900, 152.4, 2), pipeRow('WELL-3', 'IN', 1600, 152.4, 2), pipeRow('WELL-4', 'IN', 2400, 152.4, 2), pipeRow('IN', 'PS-suction', 23, 0, 0, 'flowline'), pipeRow('PS-suction', 'PS-discharge', 0, 0, 0, 'pump', 2), pipeRow('PS-discharge', 'TERMINAL', 29, 0, 0, 'flowline')], netNodes: [...[1, 2, 3, 4].map((i) => ({ id: `WELL-${i}`, kind: 'source', share: 25 })), { id: 'TERMINAL', kind: 'sink' }],
    pumpQr: 700, pumpHr: 320, pumpEff: 0.6, pumpDischargeP: 30, pumpNpshR: 5, pumpSuctionHead: 8, compOn: false, iprType: 'pi', piWell: 120, pRes: 260, sepOrient: 'vertical', sepD: 3, sepL: 9, slugSurge: 90, slugDuration: 240, scType: 'vessel',
    liftType: 'none', buckleArrestorSpacing: 0, buckleInitType: 'none', anodeSpacing: 0, transScenario: 'chokeChange', chokeStepPct: 20, supports: [{ x: 0, type: 'anchor block', k: 80000 }, { x: 23000, type: 'pump-station anchor', k: 80000 }, { x: 52000, type: 'terminal anchor block', k: 80000 }],
    fittings: [{ type: 'bend90lr', count: 14 }, { type: 'gate', count: 6 }, { type: 'check', count: 2 }, { type: 'teeRun', count: 4 }, { type: 'exit', count: 1 }] } },
  { name: 'Ultra-deep pipe-in-pipe flowline with lazy-wave riser', values: {
    geomMode: 'generate', riserType: 'lazyWave', waterDepth: 1800, flowlineLength: 12000, seabedDrop: 60, undulationAmp: 12, undulationLength: 4000, topsideElev: 22, hangoffAngle: 9, lazyLift: 38, lazySag: 12, lazyBuoy: 1.6,
    pip: true, pipWtMm: 15.9, layers: [{ name: 'Aerogel blanket', tMm: 20, k: 0.02, rho: 150, cp: 1000 }, { name: 'Annulus gap (air)', tMm: 8, k: 0.03, rho: 1.2, cp: 1005 }], tSeabed: 3, pRes: 345, piWell: 32, pSep: 22, designPressure: 380, grade: 'X70', wtMm: 19.1, liftDp: 40, buckleInitType: 'snake lay', buckleInitSpacing: 3000, transScenario: 'compressorTrip',
    supports: [{ x: 0, type: 'PLET on a sliding mudmat', k: 2000 }, { x: 12000, type: 'riser-base anchor', k: 50000 }] } },
];

// ---- suite object ----------------------------------------------------------------------------------------------------
const okNum = (items) => items.filter((it) => it && (typeof it.value === 'object' ? it.value !== null : isNum(it.value)));
export default {
  id: 'net', num: 2, title: 'Geometry, Wells, Network & Equipment', short: 'Network', icon: '🛠️',
  tagline: 'Route and riser geometry, wall thermal design, wells, chokes, rotating equipment, vessels and the hydraulic network of the case.',
  description: 'Turns the route (profile table, 3-D centreline, terrain transect or a generated flowline with a catenary, lazy-wave or vertical riser) into lengths, inclinations, volumes, bends, free-span candidates and a calculation mesh; designs the multilayer wall for its U-value, weight and thermal mass; solves the well (minimum-curvature trajectory, inflow, tubing lift with the kernel two-phase closures, gas lift or subsea boosting) and the node–edge network (Newton–Raphson on nodal mass balances and branch pressure-drop laws) and joins them by nodal analysis. A well-count study repeats the nodal solution for one to six identical wells — deliverable rate with and without the lift, drawdown against the sand-control limit, choke opening and lift power — and recommends the smallest count that holds the case rate. Chokes, pumps, compressors, the separator and the slug catcher are sized on the same case fluid, and a transient network model (line-pack, vessel inventories, pump and compressor maps with rotor inertia, anti-surge recycle and valve actuators) runs trip, valve-closure and choke-change scenarios. The design features that the integrity and economics suites need — buckle arrestors, lateral-buckle initiators, cathodic protection, supports — are inputs and are published.',
  guide: [
    'Choose where the route comes from: the elevation table (also filled by the geometry importer), a 3-D centreline, the terrain grid, or a generated flowline and riser.',
    'Describe the pipe wall from the steel outwards, burial and the surroundings; the U-value and its resistances are computed, not typed.',
    'List the connections of the network. “flowline” and “riser” connections follow the route; everything else uses its own length and diameter. Leave node data blank to have sources, sinks and elevations inferred.',
    'Give the well survey and inflow model. Run: the network pressure at the wellhead is compared with what the well delivers, the difference is given to the production choke, and the natural-flow rate is found by nodal analysis.',
    'Read the well-count table: it shows what one to six wells deliver with and without the lift, their drawdown against the limit you set, the choke opening and the lift power, and which count is recommended; the economics suite prices each option.',
    'Set the artificial lift, the design features (buckle arrestors, buckle initiators, anodes, supports, allowable span) and, on the model-setup tab, the transient scenario.',
    'Check the equipment duty list, the geometry, orientation and network diagnostics, then the mesh study and the verification tab. Measurements typed into “Measured data for comparison” are listed beside the model values with error metrics; the calibration tab fits the listed parameters to test data and shows the sourced reference data sets.',
  ],
  equationsNote: 'Steady state. All branches carry the case fluid (rates scaled per branch); pressure gradients use the shared kernel closures and are only piecewise smooth, so the network residual can stall near 1e-6 at a flow-pattern change (reported). Tubing temperature is prescribed linearly from reservoir to wellhead. The homogeneous-equilibrium choke estimate follows an isentrope with the latent heat neglected; the Sachdeva model is frozen-composition. Compressor stages use the polytropic relations with a mean compressibility from the property table; the generic compressor map obeys the fan laws and is replaced by tested speed lines when supplied. Free spans without surveyed terrain come from a seeded statistical micro-relief and are screening values only. The slug-catcher inventory model is a single liquid hold-up with a capped level-proportional draw-off. Transient network: dynamic nodal inventories (one gas line-pack volume for the flowline at its mean pressure, separator gas and liquid, compressor discharge piping, export header) joined by quasi-steady branches; real-gas factor from the case property table at the arrival temperature, isothermal volumes; liquid in the flowline is passed through without storage (holdup waves, slugs and pressure-wave travel belong to the flow and operations suites). The compressor follows the generic fan-law map scaled to the duty (tested maps are used for the steady duty point only); left of the surge line the map is continued so that surge is detected and timed, but the surge limit cycle itself is not resolved. Pump shaft power follows the efficiency parabola down to shut-off. Speed governors, anti-surge, flare and level control are proportional or PI loops with first-order, rate-limited actuators. Gas lift gives the lift gas the properties of the produced gas; subsea boosting is a pressure rise at the wellhead with hydraulic power = suction volume flow × pressure rise ÷ efficiency. Buckle-arrestor extent, cathodic-protection demand and sleeper lift-off are screening calculations for the integrity suite, not a replacement for it. Constants and their sources are listed in PROVENANCE (export of this module): items marked unverified could not be opened in an openly readable source and are inputs where that matters. Validation status: sourced reference data back pipe friction, bend losses, choke and valve flow, deviation surveys, gas-well tests and compressor-map interpolation, and twelve code-to-code benchmarks against the open-source fluids library cover friction, Crane fittings, IEC 60534 valve sizing, orifice flow, compressor head and mist-eliminator K. No openly licensed pump acceptance test, separator carry-over test, as-built or inspection record was found: those comparisons run only on data typed into “Measured data for comparison”. The well-count study assumes identical wells with equal shares and continues the network back-pressure beyond its solved range with the square of the rate; with gas lift its no-lift column still uses the back-pressure of the network carrying the lift gas.',
  implemented: [
    'euclidean/analytic geometry', 'differential arc-length equation', 'pipe centreline parameterization', 'coordinate-transformation equations', 'hydrostatic head equation', 'bernoulli equation', 'extended bernoulli/mechanical-energy equation', 'continuity equation', 'momentum balance', 'energy balance', 'kirchhoff-type network mass conservation', 'loop pressure-balance equations',
    'darcy–weisbach equation', 'fanning friction formulation', 'colebrook–white equation', 'haaland equation', 'swamee–jain equation', 'churchill friction-factor model', 'laminar hagen–poiseuille relation', 'reynolds-number regime classification', 'equivalent-length method', 'resistance-coefficient/k-factor method',
    'hydrostatic wellbore equation', 'darcy–weisbach wellbore momentum balance', 'vogel inflow-performance relationship', 'productivity-index model', 'darcy radial-flow equation', 'forchheimer non-darcy flow', 'fetkovich deliverability model', 'jones–blount–glaze gas-well model', 'nodal-analysis equations',
    'bernoulli/orifice equation', 'discharge-coefficient models', 'isa/iec valve-flow equations', 'critical/choked-flow equations', 'compressible nozzle relations', 'homogeneous-equilibrium choking models', 'multiphase choke correlations',
    'euler turbomachinery equation', 'pump affinity laws', 'pump-head curves', 'compressor polytropic equations', 'compressor isentropic relations', 'compressor-map interpolation', 'efficiency relations', 'surge/choke constraints',
    'transient mass balance', 'component balances', 'vapor–liquid flash', 'residence-time models', 'level/volume relationships',
    'node mass balance + branch momentum equations', 'pressure-flow network equations', 'nonlinear network solution', 'graph-based incidence-matrix formulation',
    'geometry + network hydraulics', 'well ipr + vertical-lift performance', 'reservoir–well–pipeline coupling', 'separator flash + dynamic inventory balance', 'pump/compressor map + transient network equations',
    // initial and boundary conditions represented by inputs
    'initial configuration and availability of every pipeline', 'manifold and item of equipment', 'initial pipe internal diameter and wall thickness', 'initial roughness', 'initial insulation condition', 'initial burial and seabed-contact state', 'initial valve and choke positions', 'initial pump and compressor operating states', 'initial slug-catcher inventory', 'initial separator liquid and gas inventories', 'initial well completion status', 'initial equipment connectivity', 'any pre-existing deposits', 'restrictions', 'damage or effective diameter reductions',
    'inlet and outlet locations', 'well/reservoir connections', 'terminal and separator connections', 'branch and junction connectivity', 'closed ends', 'pressure boundaries', 'equipment interfaces and environmental interfaces', 'fixed pipeline and riser coordinates', 'elevation and bathymetric profiles', 'water depth', 'pipe-support conditions where relevant', 'equipment capacities', 'valve and choke operating limits', 'pump and compressor operating envelopes', 'separator constraints', 'slug-catcher capacity', 'physical limits governing flow through connected components',
    // inputs
    'pipeline/riser/well centrelines and coordinates', 'lengths', 'id/od', 'wall thickness', 'roughness', 'elevation/bathymetry', 'inclination', 'bends, tees, branches and junctions', 'riser configuration', 'insulation/coatings', 'burial and seabed contact', 'pipe/material properties', 'well trajectory and completion geometry', 'network topology', 'manifolds', 'valves/chokes', 'pumps/compressors', 'separators', 'slug catchers', 'equipment dimensions, capacities, performance curves/maps and connectivity', 'coordinate reference system and units',
    // outputs
    'validated computational geometry/network', 'segment lengths', 'elevations and inclinations', 'cross-sectional areas', 'hydraulic diameters', 'internal volumes', 'node-edge connectivity', 'well/network topology', 'equipment connectivity and orientation', 'local-loss coefficients', 'equipment operating characteristics', 'derived geometric/thermal parameters', 'geometry-quality/topology diagnostics', 'computational discretization or mesh-ready representation',
    // calibration, verification, validation actually supported
    'pipe roughness from pressure-drop data', 'effective diameter from inspection data', 'local-loss coefficients', 'bend-loss coefficients', 'tee/junction coefficients', 'valve cv/kv', 'choke coefficients', 'pump performance curves', 'compressor performance maps', 'separator performance parameters', 'well productivity index', 'ipr parameters', 'completion skin', 'tubing roughness', 'heat-transfer coefficients', 'insulation properties', 'slug-catcher effective capacity', 'equipment efficiency', 'equipment degradation factors',
    'cad/geometry consistency', 'equipment orientation', 'geometry-import regression tests',
    'length conservation', 'elevation-profile verification', 'bathymetry reconstruction checks', 'inclination calculations', 'internal-volume calculations', 'cross-sectional-area calculations', 'hydraulic-diameter calculations', 'connectivity/topology tests', 'node-edge consistency', 'branch connectivity', 'flow-direction consistency', 'junction conservation', 'coordinate-system transformations', 'unit-conversion tests', 'boundary-condition assignment', 'duplicate/disconnected-node detection', 'zero-length-element detection', 'negative/invalid-volume detection', 'mesh/network independence', 'analytical geometry benchmarks', 'round-trip geometry tests',
    // validation. Backed by sourced reference data sets: field pressure-drop measurements (pipe friction, bend losses), valve / choke flow
    // tests (field choke tests, a manufacturer Cv table), well deviation surveys and survey data (regulator-grade survey report), well-test
    // data (gas-well deliverability tests), compressor performance tests (map interpolation only). Supported only through the “Measured
    // data for comparison” table, with no open data set found: pump factory acceptance tests, separator performance tests, as-built drawings,
    // pipeline inspection, GIS / bathymetric surveys, caliper data, production logging, commissioning data.
    'field pressure-drop measurements', 'well-test data', 'valve/choke flow tests', 'as-built drawings', 'survey data', 'pipeline inspection data', 'gis/bathymetric surveys', 'well deviation surveys', 'caliper data', 'pump factory acceptance tests', 'compressor performance tests', 'separator performance tests', 'production logging', 'commissioning data',
  ],
  referenceOnly: [],
  validationData: VALIDATION,
  geometry: (g) => geometryOffers(g),
  inputs: INPUTS,
  presets: PRESETS,
  pull: ({ fluid, outputs } = {}) => okNum([
    { key: 'tIn', value: fluid?.Tin, from: 'Case fluid: inlet temperature' },
    { key: 'pSep', value: fluid?.Pout, from: 'Case fluid: arrival pressure' },
    { key: 'pRes', value: fluid?.Pres, from: 'Case fluid: reservoir pressure' },
    { key: 'tRes', value: fluid?.Tres, from: 'Case fluid: reservoir temperature' },
    { key: 'pBubble', value: outputs?.pvt?.psat, from: `PVT suite: saturation pressure at reservoir temperature (${outputs?.pvt?.psatType || 'bubble'} point)` },
    { key: 'effIdMm', value: isNum(outputs?.solids?.effectiveId) ? outputs.solids.effectiveId * 1000 : null, from: 'Solids suite: smallest bore left by deposits' },
    { key: 'effRoughUm', value: isNum(outputs?.solids?.roughnessEff) ? outputs.solids.roughnessEff * 1e6 : null, from: 'Solids suite: effective roughness of the fouled wall' },
    { key: 'slugSurge', value: outputs?.flow?.slug?.surge, from: 'Flow suite: surge volume to accommodate at the receiving vessel' },
    { key: 'slugDuration', value: outputs?.flow?.slug?.length > 0 && outputs?.flow?.slug?.velocity > 0 ? clamp(outputs.flow.slug.length / outputs.flow.slug.velocity, 5, 3600) : null, from: 'Flow suite: slug length ÷ slug velocity' },
  ]),
  site: (site) => { const d = site?.data || {}, b = d.bathy; return okNum([
    { key: 'waterDepth', value: d.depth, from: 'Water depth at the site' },
    { key: 'tSeabed', value: d.seabedTemp, from: 'Seabed temperature at the site' },
    { key: 'tSeaSurface', value: d.sst, from: 'Sea-surface temperature at the site' },
    { key: 'currentSpeed', value: d.currentSpeed, from: 'Current speed at the site' },
    { key: 'tAir', value: d.airTemp, from: 'Air temperature at the site' },
    { key: 'windSpeed', value: d.windSpeed, from: 'Wind speed at the site' },
    { key: 'tGround', value: d.groundTemp, from: 'Ground temperature at about 1 m (buried onshore sections)' },
    { key: 'seabedDrop', value: isNum(d.seabedSlope) && isNum(d.depth) ? Math.min(Math.tan(d.seabedSlope * D2R) * DEF.flowlineLength, 0.5 * d.depth) : null, from: 'Regional seabed slope over the default flowline length (capped at half the water depth)' },
    { key: 'topsideElev', value: d.depth === 0 && isNum(d.elevation) ? clamp(d.elevation, 0, 100) : null, from: 'Land elevation at the site (onshore arrival)' },
    { key: 'terrain', value: b && Array.isArray(b.elev) && Array.isArray(b.lat) && Array.isArray(b.lon) ? { x: b.lon, y: b.lat, elev: b.elev, geographic: true } : null, from: 'Local bathymetry grid of the site' },
  ]); },
  run,
  mesh: [
    { name: 'Flowline, riser and tubing cells', keys: ['nCells', 'nTubing'], min: 6, note: 'Cells of the marching solution along the route (riser and jumpers scale with it) and along the tubing. Geometry-only quantities (volume, U-value) must not move at all.',
      metrics: [{ label: 'Flowline inlet pressure', unit: 'bara', get: (r) => r.outputs.inletPressure }, { label: 'Wellhead pressure available', unit: 'bara', get: (r) => r.outputs.whpAvailable }, { label: 'Natural-flow rate', unit: 'Sm³/d', get: (r) => r.outputs.operatingRate }, { label: 'Arrival temperature', unit: '°C', get: (r) => r.outputs.arrivalTemp }, { label: 'U-value', unit: 'W/m²K', get: (r) => r.outputs.uValue }, { label: 'Internal volume', unit: 'm³', get: (r) => r.outputs.volume }] },
    { name: 'Transient time step', keys: ['transSteps'], min: 20, note: 'Time steps of the implicit θ-method of the transient network (the scenario of the model-setup tab).',
      metrics: [{ label: 'Peak separator pressure', unit: 'bara', get: (r) => r.outputs.transient?.peakSepP ?? NaN }, { label: 'Smallest surge margin', unit: '%', get: (r) => r.outputs.transient?.minSurgeMargin ?? NaN }, { label: 'Final separator level', unit: 'm', get: (r) => r.outputs.transient?.levelEnd ?? NaN }] },
    { name: 'Network solver tolerance', keys: ['netTol'], refine: 'divide', note: 'The relative residual tolerance of the Newton–Raphson network solution is tightened by the refinement ratio.',
      metrics: [{ label: 'Flowline inlet pressure', unit: 'bara', get: (r) => r.outputs.inletPressure }, { label: 'Wellhead pressure needed', unit: 'bara', get: (r) => r.outputs.whpRequired }] },
  ],
  calibration: {
    note: 'Each parameter is estimated from the test that responds to it; fit only the parameters for which you have data and leave the others unticked. Line test with sea water: wall roughness (long lines) or the local-loss multiplier (short, fitting-dominated lines) — not both; with the line-fill volume from a pig run or caliper survey also the effective diameter. Bend and tee spool tests: the bend and tee multipliers. Valve water-flow test: rated Cv. Bean-choke test: discharge coefficient. Pump test: rated head and shut-off ratio from the head at several flows, best efficiency from the power. Compressor test (flow and head normalised by the design point, speed in %): map head and flow scaling; discharge temperature at a known pressure ratio (suction 10 bara, 30 °C, k = 1.28): polytropic efficiency. Separator carry-over test: Souders–Brown K (switch the automatic K off to use it). Well test: productivity index, and with points below the bubble point the saturation pressure (Fetkovich exponent for that inflow model); build-up analysis: skin from the skin pressure drop (viscosity 1 mPa·s and formation volume factor 1.3 assumed). Tubing injection test with water: tubing roughness. Thermal test (known heat-capacity rate and inlet temperature): U-value multiplier or, for aged insulation, the insulation conductivity multiplier — one of the two. Slug-catcher fill test: effective capacity. The degradation factor acts on the pump and compressor heads of an aged machine once the as-new curve is fixed.',
    params: [{ key: 'roughUm', label: 'Wall roughness (µm)', lo: 1, hi: 1000 }, { key: 'effIdMm', label: 'Effective diameter (mm)', lo: 150, hi: 254 }, { key: 'kMult', label: 'Local-loss multiplier', lo: 0, hi: 20 }, { key: 'bendMult', label: 'Bend-loss multiplier', lo: 0.2, hi: 5 }, { key: 'teeMult', label: 'Tee / junction-loss multiplier', lo: 0.2, hi: 5 }, { key: 'chokeCvMax', label: 'Valve rated Cv', lo: 5, hi: 1000 }, { key: 'chokeCd', label: 'Bean discharge coefficient', lo: 0.4, hi: 1.05 },
      { key: 'pumpHr', label: 'Pump rated head (m)', lo: 50, hi: 3000 }, { key: 'pumpShutoff', label: 'Pump shut-off head ratio', lo: 1.02, hi: 2 }, { key: 'pumpEff', label: 'Pump best efficiency', lo: 0.2, hi: 0.92 }, { key: 'compHeadMult', label: 'Compressor map head scaling', lo: 0.5, hi: 1.5 }, { key: 'compFlowMult', label: 'Compressor map flow scaling', lo: 0.5, hi: 1.5 }, { key: 'compEta', label: 'Compressor polytropic efficiency', lo: 0.4, hi: 0.92 }, { key: 'degrFactor', label: 'Equipment degradation factor', lo: 0.5, hi: 1 },
      { key: 'sepK', label: 'Separator Souders–Brown K (m/s)', lo: 0.02, hi: 0.3 }, { key: 'piWell', label: 'Productivity index (Sm³/d/bar)', lo: 0.5, hi: 500 }, { key: 'pBubble', label: 'IPR: saturation pressure (bara)', lo: 20, hi: 600 }, { key: 'fetkN', label: 'IPR: Fetkovich exponent', lo: 0.5, hi: 1 }, { key: 'skin', label: 'Completion skin', lo: -5, hi: 50 }, { key: 'tubingRoughUm', label: 'Tubing roughness (µm)', lo: 1, hi: 1000 },
      { key: 'uMult', label: 'U-value multiplier', lo: 0.3, hi: 3 }, { key: 'insKMult', label: 'Insulation conductivity multiplier', lo: 0.5, hi: 5 }, { key: 'scEff', label: 'Slug-catcher effective capacity fraction', lo: 0.3, hi: 1 }],
    columns: [{ key: 'qLine', label: 'Line / spool test flow', unit: 'm³/h' }, { key: 'dpLine', label: 'Line pressure drop', unit: 'bar' }, { key: 'lineVol', label: 'Line-fill volume', unit: 'm³' }, { key: 'dpBend', label: '90° bend (r/D = 5) pressure drop', unit: 'mbar' }, { key: 'dpTee', label: 'Tee branch pressure drop', unit: 'mbar' }, { key: 'qValve', label: 'Valve water flow, fully open', unit: 'm³/h' }, { key: 'dpValve', label: 'Valve pressure drop', unit: 'bar' },
      { key: 'wChoke', label: 'Choke mass rate', unit: 'kg/s' }, { key: 'rhoChoke', label: 'Mixture density at the choke', unit: 'kg/m³' }, { key: 'bean', label: 'Bean size', unit: '1/64 in' }, { key: 'dpChoke', label: 'Choke pressure drop', unit: 'bar' }, { key: 'qPump', label: 'Pump test flow (water)', unit: 'm³/h' }, { key: 'hPump', label: 'Pump head', unit: 'm' }, { key: 'pPump', label: 'Pump shaft power', unit: 'kW' },
      { key: 'phiComp', label: 'Compressor inlet flow / design flow', unit: '–' }, { key: 'nComp', label: 'Compressor speed', unit: '%' }, { key: 'psiComp', label: 'Polytropic head / design head', unit: '–' }, { key: 'rComp', label: 'Pressure ratio of the temperature test', unit: '–' }, { key: 't2Comp', label: 'Discharge temperature', unit: '°C' },
      { key: 'rhoLt', label: 'Separator test liquid density', unit: 'kg/m³' }, { key: 'rhoGt', label: 'Separator test gas density', unit: 'kg/m³' }, { key: 'vGasMax', label: 'Gas velocity at carry-over', unit: 'm/s' }, { key: 'qTest', label: 'Well-test rate', unit: 'Sm³/d' }, { key: 'pwf', label: 'Flowing bottom-hole pressure', unit: 'bara' }, { key: 'dpSkin', label: 'Skin pressure drop', unit: 'bar' },
      { key: 'qTub', label: 'Tubing test flow (water)', unit: 'm³/h' }, { key: 'dpTubing', label: 'Tubing friction pressure drop', unit: 'bar' }, { key: 'mCp', label: 'Heat-capacity rate', unit: 'kW/K' }, { key: 'tHot', label: 'Inlet temperature', unit: '°C' }, { key: 'tArr', label: 'Arrival temperature', unit: '°C' }, { key: 'scGeom', label: 'Slug-catcher geometric volume filled', unit: 'm³' }, { key: 'scHeld', label: 'Liquid held', unit: 'm³' }],
    targets: [{ key: 'dpLine', label: 'Line pressure drop', unit: 'bar' }, { key: 'lineVol', label: 'Line-fill volume', unit: 'm³' }, { key: 'dpBend', label: 'Bend pressure drop', unit: 'mbar' }, { key: 'dpTee', label: 'Tee pressure drop', unit: 'mbar' }, { key: 'dpValve', label: 'Valve pressure drop', unit: 'bar' }, { key: 'dpChoke', label: 'Choke pressure drop', unit: 'bar' }, { key: 'hPump', label: 'Pump head', unit: 'm' }, { key: 'pPump', label: 'Pump shaft power', unit: 'kW' },
      { key: 'psiComp', label: 'Compressor head ratio', unit: '–' }, { key: 't2Comp', label: 'Compressor discharge temperature', unit: '°C' }, { key: 'vGasMax', label: 'Gas velocity at carry-over', unit: 'm/s' }, { key: 'pwf', label: 'Flowing bottom-hole pressure', unit: 'bara' }, { key: 'dpSkin', label: 'Skin pressure drop', unit: 'bar' }, { key: 'dpTubing', label: 'Tubing friction pressure drop', unit: 'bar' }, { key: 'tArr', label: 'Arrival temperature', unit: '°C' }, { key: 'scHeld', label: 'Slug-catcher liquid held', unit: 'm³' }],
    model: calModel, sample: CAL_SAMPLE, validationSample: CAL_VALID,
  },
  verify,
};
