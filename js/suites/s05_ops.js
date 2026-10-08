// Suite 5 — Operations, Control & Flow-Assurance Management.
// Shutdown/cooldown (radial finite-volume conduction + lumped check), restart and ramp-up, blowdown, pigging,
// chemical injection, heating, riser-slugging control (PID / MPC / Kalman), operating logic (scheduler, alarms,
// state machine), operating envelope and optimisation, surrogate / residual correction and historical replay.
// SI inside the engines; bara, °C, h and mm at the interfaces.
import { clamp, linspace, interp1, brent, tridiag, solveLinear, rk45, nelderMead, lstsq, rng, metrics, mean, sum, isNum, fmt } from '../core/num.js';
import { fluidModel, inhibitorFor, hydrateDepression, INHIBITORS, R, VM_STD, makeFluid, phaseProps } from '../core/thermo.js';
import { G, uValue, hInside, hOutside, frictionFactor, marchSteady } from '../core/pipe.js';
import { flowPicture, caseLine, ambientAt } from '../core/caseflow.js';
import { BASE } from '../data/basecase.js';
import * as NET from './s02_net.js';
import { REF } from '../data/ref/ops.js';

const KEL = 273.15, HOUR = 3600, DAY = 86400;
const num = (x, d) => (isNum(+x) && x !== null && x !== '' ? +x : d);
const fin = (x, d = null) => (isNum(x) ? x : d);
const rd = (x, n = 2) => (isNum(x) ? +x.toFixed(n) : null);
const txt = (x, n = 2) => (isNum(x) ? +x.toFixed(n) : '—');

// ---- special functions and small dense linear algebra ----------------------------------------------------
/** Complementary error function (Chebyshev fit, relative error < 1.2e-7), valid for any real x. */
export function erfc(x) {
  const z = Math.abs(x), t = 1 / (1 + 0.5 * z);
  const r = t * Math.exp(-z * z - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 + t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))));
  return x >= 0 ? r : 2 - r;
}
const normCdf = (x) => 0.5 * erfc(-x / Math.SQRT2), normPdf = (x) => Math.exp(-0.5 * x * x) / Math.sqrt(2 * Math.PI);
const zeros = (n, m) => Array.from({ length: n }, () => new Array(m).fill(0));
const eye = (n) => zeros(n, n).map((r, i) => { r[i] = 1; return r; });
const mm = (X, Y) => X.map((r) => Y[0].map((_, j) => { let s = 0; for (let k = 0; k < r.length; k++) s += r[k] * Y[k][j]; return s; }));
const mv = (X, v) => X.map((r) => { let s = 0; for (let k = 0; k < r.length; k++) s += r[k] * v[k]; return s; });
const tr = (X) => X[0].map((_, j) => X.map((r) => r[j]));
const madd = (X, Y, f = 1) => X.map((r, i) => r.map((v, j) => v + f * Y[i][j]));
const dot = (a, b) => { let s = 0; for (let i = 0; i < a.length; i++) s += a[i] * b[i]; return s; };
const inv = (X) => tr(X.map((_, j) => solveLinear(X, X.map((__, i) => (i === j ? 1 : 0)))));
/** Matrix exponential by scaling and squaring of a Taylor series. */
export function expm(A) {
  const n = A.length, nrm = Math.max(...A.map((r) => sum(r.map(Math.abs)))), s = Math.max(0, Math.ceil(Math.log2(Math.max(nrm, 1e-12))) + 3), B = A.map((r) => r.map((v) => v / 2 ** s));
  let E = eye(n), T = eye(n);
  for (let k = 1; k <= 14; k++) { T = mm(T, B).map((r) => r.map((v) => v / k)); E = madd(E, T); }
  for (let k = 0; k < s; k++) E = mm(E, E);
  return E;
}
/** Zero-order-hold discretisation of dx/dt = A x + B u: returns { Ad, Bd } for a sample time Ts. */
export function c2d(A, B, Ts) {
  const n = A.length, m = B[0].length, M = zeros(n + m, n + m);
  for (let i = 0; i < n; i++) { for (let j = 0; j < n; j++) M[i][j] = A[i][j] * Ts; for (let j = 0; j < m; j++) M[i][n + j] = B[i][j] * Ts; }
  const E = expm(M);
  return { Ad: E.slice(0, n).map((r) => r.slice(0, n)), Bd: E.slice(0, n).map((r) => r.slice(n)) };
}
/** Roots of a monic polynomial x^n + c[1] x^(n-1) + … + c[n] (Durand–Kerner). Returns [[re, im], …]. */
export function polyRoots(c) {
  const n = c.length - 1, rad = 1 + Math.max(...c.slice(1).map((v, k) => Math.abs(v) ** (1 / (k + 1))));
  const z = Array.from({ length: n }, (_, k) => [rad * Math.cos(0.4 + (2 * Math.PI * k) / n), rad * Math.sin(0.4 + (2 * Math.PI * k) / n)]);
  const mul = (a, b) => [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]];
  for (let it = 0; it < 600; it++) {
    let ch = 0;
    for (let i = 0; i < n; i++) {
      let p = [1, 0], q = [1, 0];
      for (let k = 1; k <= n; k++) { p = mul(p, z[i]); p[0] += c[k]; }
      for (let j = 0; j < n; j++) if (j !== i) q = mul(q, [z[i][0] - z[j][0], z[i][1] - z[j][1]]);
      const d = q[0] * q[0] + q[1] * q[1] || 1e-300, dx = [(p[0] * q[0] + p[1] * q[1]) / d, (p[1] * q[0] - p[0] * q[1]) / d];
      z[i] = [z[i][0] - dx[0], z[i][1] - dx[1]]; ch = Math.max(ch, Math.hypot(dx[0], dx[1]) / (1 + Math.hypot(z[i][0], z[i][1])));
    }
    if (ch < 1e-14) break;
  }
  return z;
}
/** Characteristic polynomial (monic, descending powers) and adjugate series of A by Faddeev–LeVerrier. */
function charPoly(A) {
  const n = A.length, c = [1], N = [eye(n)];
  let M = A.map((r) => r.slice());
  for (let k = 1; k <= n; k++) { let t = 0; for (let i = 0; i < n; i++) t += M[i][i]; const ck = -t / k; c.push(ck); const Nk = M.map((r, i) => r.map((v, j) => v + (i === j ? ck : 0))); N.push(Nk); M = mm(A, Nk); }
  return { c, N }; // (sI − A)^-1 = Σ N[k] s^(n-1-k) / p(s), k = 0..n-1
}
/** Eigenvalues of a small real matrix as [[re, im], …] (balanced by a diagonal scaling first). */
export function eig(A) {
  const n = A.length, d = A.map((r, i) => clamp(Math.sqrt(Math.max(sum(r.map((v, j) => (j === i ? 0 : Math.abs(v)))), 1e-300) / Math.max(sum(A.map((q, j) => (j === i ? 0 : Math.abs(q[i])))), 1e-300)), 1e-8, 1e8));
  const B = A.map((r, i) => r.map((v, j) => (v * (Number.isFinite(d[j] / d[i]) ? d[j] / d[i] : 1)))), sc = Math.max(...B.map((r) => sum(r.map(Math.abs)))) || 1;
  return polyRoots(charPoly(B.map((r) => r.map((v) => v / sc))).c).map((z) => [z[0] * sc, z[1] * sc]);
}
/** Transfer function C (sI − A)^-1 B of a single-input single-output state-space model: { num[], den[] } in descending powers of s. */
export function ss2tf(A, B, C) {
  const { c, N } = charPoly(A), b = B.map((r) => r[0]);
  return { den: c, num: N.slice(0, A.length).map((Nk) => dot(C, mv(Nk, b))) };
}

// ---- shutdown: radial transient conduction and cooldown --------------------------------------------------------
/**
 * Radial finite-volume grid through the pipe wall and its coatings. layers: [{ t (m), k (W/m/K), rho (kg/m³), cp (J/kg/K) }] from
 * the bore outwards; nr cells are shared between the layers in proportion to thickness (at least one each).
 * Returns { rc[], C[] (J/m/K per cell), G[] (W/m/K: fluid→cell 1, cell i→i+1, …, cell N→ambient), ro, Rtot (m·K/W) }.
 */
export function radialGrid({ ri, layers, nr = 12, hIn = 200, hOut = 500 }) {
  const L = layers.filter((l) => l.t > 0 && l.k > 0), tot = sum(L.map((l) => l.t)), cells = [];
  if (!L.length) throw new Error('The wall needs at least one layer with a positive thickness and conductivity.');
  let r = ri;
  for (const l of L) { const n = Math.max(1, Math.round((nr * l.t) / tot)), dr = l.t / n; for (let i = 0; i < n; i++) { cells.push({ rw: r, re: r + dr, rc: Math.sqrt(0.5 * (r * r + (r + dr) ** 2)), k: l.k, C: (l.rho || 1000) * (l.cp || 1000) * Math.PI * ((r + dr) ** 2 - r * r) }); r += dr; } }
  const N = cells.length, Gs = [1 / (1 / (hIn * 2 * Math.PI * ri) + Math.log(cells[0].rc / ri) / (2 * Math.PI * cells[0].k))];
  for (let i = 0; i < N - 1; i++) Gs.push((2 * Math.PI) / (Math.log(cells[i].re / cells[i].rc) / cells[i].k + Math.log(cells[i + 1].rc / cells[i].re) / cells[i + 1].k));
  Gs.push(1 / (Math.log(r / cells[N - 1].rc) / (2 * Math.PI * cells[N - 1].k) + 1 / (hOut * 2 * Math.PI * r)));
  return { rc: cells.map((c) => c.rc), C: cells.map((c) => c.C), G: Gs, ro: r, Rtot: sum(Gs.map((g) => 1 / g)) };
}
/**
 * Cooldown of one cross-section after shut-in: lumped fluid node + radial finite-volume conduction (implicit, tridiagonal).
 * o: { T0 (°C fluid at shut-in), tAmb, ri, layers, hIn, hOut, cFluid (J/m/K), nr, dt (s), nSteps, theta (1 = backward Euler, 0.5 = Crank–Nicolson),
 *      qHeat (W/m into the fluid), init: 'steady' | 'uniform' | number[] (cell temperatures) }
 * Returns { t[] (s), Tf[], Twall[] (innermost cell), Tsurf[], U (W/m²K on the bore), tau (s, lumped time constant with the effective
 *           thermal mass), lumped(t), energy: { stored0, stored, lost, heat } (J/m), cells: final temperatures }.
 */
export function cooldown(o) {
  const q = { T0: 60, tAmb: 4, hIn: 200, hOut: 500, cFluid: 5e4, nr: 12, dt: 600, nSteps: 144, theta: 1, qHeat: 0, init: 'steady', ...o }, g = q.grid || radialGrid(q), N = g.C.length, n = N + 1, th = q.theta;
  const C = [q.cFluid, ...g.C], Gs = g.G, dT = q.T0 - q.tAmb;
  let T;
  if (Array.isArray(q.init)) T = [q.T0, ...q.init];
  else if (q.init === 'uniform') T = new Array(n).fill(q.T0);
  else { T = [q.T0]; let Rc = 0; for (let i = 0; i < N; i++) { Rc += 1 / Gs[i]; T.push(q.T0 - (dT * Rc) / g.Rtot); } }
  const w = T.map((v) => (dT ? (v - q.tAmb) / dT : 1)), cEff = sum(C.map((c, i) => c * w[i])), tau = cEff * g.Rtot, energyOf = (x) => sum(C.map((c, i) => c * (x[i] - q.tAmb)));
  const a = new Array(n), b = new Array(n), c = new Array(n), d = new Array(n), t = [0], Tf = [T[0]], Tw = [T[1]], Ts = [T[N]], E0 = energyOf(T);
  let lost = 0, heat = 0;
  for (let i = 0; i < n; i++) { const gl = i > 0 ? Gs[i - 1] : 0, gr = Gs[i]; a[i] = -th * gl; c[i] = i < n - 1 ? -th * gr : 0; b[i] = C[i] / q.dt + th * (gl + gr) - (i === 0 ? th * 0 : 0); }
  b[0] = C[0] / q.dt + th * Gs[0];
  for (let s = 1; s <= q.nSteps; s++) {
    for (let i = 0; i < n; i++) {
      const gl = i > 0 ? Gs[i - 1] : 0, gr = i === 0 ? Gs[0] : Gs[i], ex = (i > 0 ? gl * (T[i - 1] - T[i]) : 0) + (i < n - 1 ? gr * (T[i + 1] - T[i]) : Gs[N] * (q.tAmb - T[i]));
      d[i] = (C[i] / q.dt) * T[i] + (1 - th) * ex + (i === 0 ? q.qHeat : 0) + (i === n - 1 ? th * Gs[N] * q.tAmb : 0);
    }
    const Tn = tridiag(a, b, c, d);
    lost += q.dt * Gs[N] * (th * (Tn[N] - q.tAmb) + (1 - th) * (T[N] - q.tAmb)); heat += q.dt * q.qHeat; T = Tn;
    t.push(s * q.dt); Tf.push(T[0]); Tw.push(T[1]); Ts.push(T[N]);
  }
  return { t, Tf, Twall: Tw, Tsurf: Ts, U: 1 / (g.Rtot * 2 * Math.PI * q.ri), tau, cEff, lumped: (x) => q.tAmb + q.qHeat * g.Rtot + (dT - q.qHeat * g.Rtot) * Math.exp(-x / tau), energy: { stored0: E0, stored: energyOf(T), lost, heat }, cells: T.slice(1), grid: g };
}
/** First time (same unit as t) at which the series y falls to the series or constant lim; null when it never does. */
export function timeBelow(t, y, lim) {
  const L = (i) => (Array.isArray(lim) ? lim[i] : lim);
  if (y[0] <= L(0)) return 0;
  for (let i = 1; i < t.length; i++) { const d0 = y[i - 1] - L(i - 1), d1 = y[i] - L(i); if (d1 <= 0) return t[i - 1] + ((t[i] - t[i - 1]) * d0) / (d0 - d1 || 1e-300); }
  return null;
}
/**
 * Settle-out after shut-in: the liquid drains into the low points between successive crests, the gas equalises at one
 * pressure with its mass conserved. st: stations { s[], z[], ds, P[] (bara), T[] (°C), holdup[], A }; fm: fluid model.
 * Returns { pSettle (bara), holdup[] (settled), gasMass (kg), liquidVol (m³), gasVol (m³), headUphill (bar: liquid legs the restart must lift), levelRiser (m of liquid column at the line end),
 *   seals: { up[], down[] } (bar per basin: liquid head a gas pocket must lift to escape towards the outlet / towards the inlet) }.
 */
export function settleOut(st, fm) {
  const n = st.s.length, A = st.A, vol = st.ds * A, mG = sum(st.P.map((p, i) => fm.at(p, st.T[i]).rhoG * (1 - st.holdup[i]) * vol)), gasVol = sum(st.holdup.map((h) => (1 - h) * vol));
  const g = (P) => sum(st.T.map((T, i) => fm.at(P, T).rhoG * (1 - st.holdup[i]) * vol)) - mG, lo = Math.min(...st.P), hi = Math.max(...st.P), pSettle = hi - lo < 1e-6 ? lo : brent(g, lo, hi, 1e-6);
  // basins: split at local crests of the elevation profile
  const H = new Array(n).fill(0), cuts = [0];
  for (let i = 1; i < n - 1; i++) if (st.z[i] > st.z[i - 1] + 1e-6 && st.z[i] >= st.z[i + 1]) cuts.push(i);
  cuts.push(n);
  const half = Math.sqrt(A / Math.PI); // a cell spans its own elevation range ± the pipe radius
  for (let b = 0; b < cuts.length - 1; b++) {
    const i0 = cuts[b], i1 = cuts[b + 1], V = sum(st.holdup.slice(i0, i1)) * vol, zl = (i) => st.z[i] - Math.max(half, 0.5 * Math.abs((st.dz || [])[i] || 0)), zh = (i) => st.z[i] + Math.max(half, 0.5 * Math.abs((st.dz || [])[i] || 0));
    const fill = (lev) => { let v = 0; for (let i = i0; i < i1; i++) v += clamp((lev - zl(i)) / (zh(i) - zl(i)), 0, 1) * vol; return v; };
    let zmin = Infinity, zmax = -Infinity; for (let i = i0; i < i1; i++) { zmin = Math.min(zmin, zl(i)); zmax = Math.max(zmax, zh(i)); }
    const lev = V >= fill(zmax) ? zmax : brent((x) => fill(x) - V, zmin, zmax, 1e-9);
    for (let i = i0; i < i1; i++) H[i] = clamp((lev - zl(i)) / (zh(i) - zl(i)), 0, 1);
  }
  const rhoL = fm.at(pSettle, mean(st.T)).rhoL; let head = 0, lev = 0;
  for (let i = 0; i < n; i++) { const dz = (st.dz || [])[i] || 0; if (dz > 0) head += rhoL * G * dz * H[i]; }
  const up = [], down = [];
  for (let b = 0; b < cuts.length - 1; b++) { let u = 0, d = 0; for (let i = cuts[b]; i < cuts[b + 1]; i++) { const dz = (st.dz || [])[i] || 0; if (dz > 0) u += rhoL * G * dz * H[i]; else d -= rhoL * G * dz * H[i]; } up.push(u / 1e5); down.push(d / 1e5); }
  for (let i = n - 1; i >= 0 && st.z[i] > st.z[Math.max(i - 1, 0)] - 1e-9; i--) lev += H[i] * Math.max((st.dz || [])[i] || 0, 0);
  return { pSettle, holdup: H, gasMass: mG, gasVol, liquidVol: sum(st.holdup) * vol, headUphill: head / 1e5, rhoL, levelRiser: lev, seals: { up, down } };
}
/**
 * Warm-up after restart: plug-flow fluid energy equation with wall storage, implicit upwind marching (unconditionally stable).
 * o: { s[] (cell centres), ds, tAmb[], T0[] (fluid = wall at restart), mdot(t) (kg/s), cp (J/kg/K), cFluid[] , cWall[] (J/m/K), gIn, gOut (W/m/K), tIn(t) | number, dt, nSteps, ref[] }
 * Returns { t[], Tout[], Tmin[] (smallest T − ref along the line; ref defaults to 0), T (final fluid profile), field: [[…]] (every `every` steps), energy: { in, out, lost, stored } (J) }.
 */
export function warmUp(o) {
  const n = o.s.length, T = o.T0.slice(), Tw = (o.Tw0 || o.T0).slice(), ref = o.ref || new Array(n).fill(0), minOver = () => { let m = Infinity; for (let j = 0; j < n; j++) m = Math.min(m, T[j] - ref[j]); return m; }, t = [0], Tout = [T[n - 1]], Tmin = [minOver()], field = [T.slice()], tf = [0], every = o.every || Math.max(1, Math.round(o.nSteps / 40));
  const eOf = () => { let e = 0; for (let j = 0; j < n; j++) e += (o.cFluid[j] * T[j] + o.cWall[j] * Tw[j]) * o.ds; return e; }, E0 = eOf();
  let eIn = 0, eOut = 0, lost = 0;
  for (let k = 1; k <= o.nSteps; k++) {
    const tm = k * o.dt, m = Math.max(o.mdot(tm), 0), tin = typeof o.tIn === 'function' ? o.tIn(tm) : o.tIn, adv = (m * o.cp) / o.ds;
    let up = tin;
    for (let j = 0; j < n; j++) {
      // [a11 a12; a21 a22] [T; Tw] = [b1; b2]
      const a11 = o.cFluid[j] / o.dt + adv + o.gIn, a12 = -o.gIn, a21 = -o.gIn, a22 = o.cWall[j] / o.dt + o.gIn + o.gOut, b1 = (o.cFluid[j] / o.dt) * T[j] + adv * up, b2 = (o.cWall[j] / o.dt) * Tw[j] + o.gOut * o.tAmb[j], det = a11 * a22 - a12 * a21;
      T[j] = (b1 * a22 - a12 * b2) / det; Tw[j] = (a11 * b2 - a21 * b1) / det; up = T[j]; lost += o.gOut * (Tw[j] - o.tAmb[j]) * o.ds * o.dt;
    }
    eIn += m * o.cp * tin * o.dt; eOut += m * o.cp * T[n - 1] * o.dt;
    t.push(tm); Tout.push(T[n - 1]); Tmin.push(minOver()); if (k % every === 0) { field.push(T.slice()); tf.push(tm); }
  }
  return { t, Tout, Tmin, T, Tw, field, tField: tf, energy: { in: eIn, out: eOut, lost, stored: eOf() - E0 } };
}

// ---- depressurisation / blowdown ----------------------------------------------------------------------------------
/** Critical pressure ratio of an ideal gas with isentropic exponent k. */
export const criticalRatio = (k) => (2 / (k + 1)) ** (k / (k - 1));
/** Mass flux (kg/s/m²) of a gas through an orifice from P, T (Pa, K) to pBack: choked or sub-critical isentropic nozzle flow with compressibility Z. */
export function orificeFlux(P, T, pBack, { k = 1.3, mw = 0.02, Z = 1 } = {}) {
  if (!(P > pBack)) return 0;
  const r = Math.max(pBack / P, 0), rho = (P * mw) / (Z * R * T);
  if (r <= criticalRatio(k)) return Math.sqrt(k * P * rho) * (2 / (k + 1)) ** ((k + 1) / (2 * (k - 1)));
  return Math.sqrt(2 * P * rho * (k / (k - 1)) * (r ** (2 / k) - r ** ((k + 1) / k)));
}
/**
 * Homogeneous-equilibrium two-phase mass flux by the ω-method (Leung): P0 (Pa), rho0 (kg/m³ mixture), omega (compressibility
 * parameter, α0/k for a non-flashing gas–liquid mixture), pBack. Returns { G (kg/s/m²), etaC (critical pressure ratio), choked }.
 */
export function omegaFlux(P0, rho0, omega, pBack) {
  // the critical-ratio equation always has one root in (0, 1): f → −∞ as η → 0 and f(1) = 1
  const w = Math.max(omega, 1e-6), fc = (e) => e * e + (w * w - 2 * w) * (1 - e) ** 2 + 2 * w * w * Math.log(e) + 2 * w * w * (1 - e), etaC = brent(fc, 1e-12, 1, 1e-13), eta = pBack / P0;
  if (eta <= etaC) return { G: etaC * Math.sqrt((P0 * rho0) / w), etaC, choked: true };
  if (eta >= 1) return { G: 0, etaC, choked: false };
  return { G: (Math.sqrt(P0 * rho0) * Math.sqrt(-2 * (w * Math.log(eta) + (w - 1) * (1 - eta)))) / (w * (1 / eta - 1) + 1), etaC, choked: false };
}
/**
 * Blowdown of a gas inventory through a restriction: transient mass and energy balance of the gas with wall heat pick-up.
 * o: { V (m³ gas volume), P0 (Pa), T0 (K), pBack (Pa), area (m²), cd, k, mw (kg/mol), Z (number | (P, T) => Z), mode: 'isothermal' | 'adiabatic' | 'wall',
 *      wallC (J/K), wallUA (W/K gas ↔ wall), extUA (W/K wall ↔ ambient), tAmb (K), Tw0 (K), pEnd (Pa), pMark (Pa, extra pressure whose time is reported),
 *      n (time steps to pEnd, the step is fixed from the initial rate), jt (K/Pa across the valve),
 *      flash: (P, T) => kg of gas liberated per Pa of pressure drop (solution gas; cools the wall node by `latent` J/kg),
 *      hem: { mLiq (kg), rhoL, frac (liquid volume fraction at the valve inlet) } for two-phase discharge,
 *      relax: relaxation time of the gas liberation in s (homogeneous relaxation closure dm/dt = (m_eq − m)/Θ; 0 or absent = equilibrium),
 *             or 'dz' for the Downar-Zapolski correlation evaluated with the gas volume fraction of the inventory and the superheat pressure,
 *      liqVol (m³ liquid in the inventory, for the void fraction of the correlation) }
 * Returns { t[], P[], T[], Tw[], mdot[] (total), m[], tEnd (s to pEnd | null), tMark, minT, minTw, minTdown (K, downstream of the valve), peak (kg/s),
 *           discharged, flashed, mass: { initial, final }, liquidOut (kg) }.
 */
export function blowdown(o) {
  const q = { cd: 0.85, k: 1.3, mw: 0.02, Z: 1, mode: 'wall', wallC: 0, wallUA: 0, extUA: 0, tAmb: 277, n: 400, jt: 0, latent: 3e5, maxFactor: 8, ...o }, Zf = typeof q.Z === 'function' ? q.Z : () => q.Z;
  if (!(q.V > 0 && q.P0 > 0 && q.T0 > 0 && q.area > 0)) throw new Error('Blowdown needs a positive gas volume, pressure, temperature and orifice area.');
  const pEnd = Math.max(q.pEnd ?? q.pBack * 1.05, q.pBack * 1.0005), Tw0 = q.Tw0 ?? q.T0, rhoOf = (P, T) => (P * q.mw) / (Zf(P, T) * R * T);
  let mL = q.hem ? q.hem.mLiq : 0, Vg = q.V, m = rhoOf(q.P0, q.T0) * Vg, T = q.T0, Tw = Tw0, P = q.P0;
  const flow = (P, T, mLiq) => { // total mass rate and its gas fraction
    if (q.hem && mLiq > 0 && q.hem.frac > 0) { const a = 1 - q.hem.frac, rg = rhoOf(P, T), r0 = a * rg + (1 - a) * q.hem.rhoL, f = omegaFlux(P, r0, a / q.k, q.pBack); return [q.cd * q.area * f.G, (a * rg) / r0]; }
    return [q.cd * q.area * orificeFlux(P, T, q.pBack, { k: q.k, mw: q.mw, Z: Zf(P, T) }), 1];
  };
  const [w0] = flow(P, T, mL), tau0 = m / Math.max(w0, 1e-12), dt = (q.dt ?? (1.7 * tau0 * Math.log(q.P0 / pEnd) + 0.2 * tau0)) / q.n * (q.dt ? q.n : 1);
  const relaxOn = q.flash && (q.relax === 'dz' || q.relax > 0);
  const der = (s, flFix = null) => { // s = [m, T, Tw, mL]
    const Pn = (s[0] * Zf(P, s[1]) * R * s[1]) / (q.mw * Vg), [w, xg] = flow(Pn, s[1], s[3]), cv = (Zf(Pn, s[1]) * R) / (q.mw * (q.k - 1)), wg = w * xg, wl = w - wg;
    const FP = q.flash ? Math.max(0, q.flash(Pn, s[1])) * Pn : 0, fl = flFix ?? (wg * FP) / (Math.max(s[0], 1e-9) + FP); // liberation follows the actual pressure fall: dP/dt ≈ −P (wg − fl)/m
    const dT = q.mode === 'isothermal' ? 0 : (-(wg * Zf(Pn, s[1]) * R * s[1]) / q.mw + fl * cv * q.k * (s[2] - s[1])) / (Math.max(s[0], 1e-9) * cv);
    const dTw = q.mode === 'wall' && q.wallC > 0 ? (q.extUA * (q.tAmb - s[2]) - fl * q.latent) / q.wallC : 0;
    return { d: [-wg + fl, dT, dTw, -wl], w, wg, fl, P: Pn };
  };
  const out = { t: [0], P: [P], T: [T], Tw: [Tw], mdot: [w0], m: [m] }, m0 = m, mL0 = mL;
  let disc = 0, flashed = 0, tEnd = null, tMark = null, minT = T, minTw = Tw, minTd = T - q.jt * (P - q.pBack), peak = w0, time = 0, held = 0, thetaLast = 0;
  for (let k = 0; k < q.n * q.maxFactor && tEnd === null; k++) {
    const s0 = [m, T, Tw, mL]; let a = der(s0), flFix = null;
    if (relaxOn) { // gas the equilibrium would release in this step is held back and released with the relaxation time
      const coef = Math.max(q.flash(P, T), 1e-300), th = q.relax === 'dz' ? hrmRelaxationTime(Vg / (Vg + (q.liqVol ?? (q.hem ? mL / q.hem.rhoL : 0))), P + held / coef, P) : q.relax, rs = relaxStep(held, a.fl, th, dt);
      held = rs.E; flFix = rs.released / dt; thetaLast = th; a = der(s0, flFix);
    }
    const s1 = s0.map((v, i) => v + dt * a.d[i]); s1[0] = Math.max(s1[0], 1e-9); s1[1] = Math.max(s1[1], 20); s1[3] = Math.max(s1[3], 0);
    const b = der(s1, flFix), Pold = P;
    m = Math.max(m + 0.5 * dt * (a.d[0] + b.d[0]), 1e-9); T = Math.max(T + 0.5 * dt * (a.d[1] + b.d[1]), 20); Tw += 0.5 * dt * (a.d[2] + b.d[2]); mL = Math.max(mL + 0.5 * dt * (a.d[3] + b.d[3]), 0);
    disc += 0.5 * dt * (a.w + b.w); flashed += 0.5 * dt * (a.fl + b.fl); time += dt;
    if (q.mode === 'wall' && q.wallUA > 0) { // gas ↔ wall exchange integrated exactly over the step (stiff when little gas is left)
      const Cg = m * ((Zf(P, T) * R) / (q.mw * (q.k - 1))), Cw = q.wallC > 0 ? q.wallC : Infinity, Teq = Number.isFinite(Cw) ? (Cg * T + Cw * Tw) / (Cg + Cw) : Tw, f = Math.exp(-dt * q.wallUA * (1 / Cg + (Number.isFinite(Cw) ? 1 / Cw : 0)));
      T = Teq + (T - Teq) * f; if (Number.isFinite(Cw)) Tw = Teq + (Tw - Teq) * f;
    }
    if (q.hem) Vg = q.V + (mL0 - mL) / q.hem.rhoL;
    P = (m * Zf(P, T) * R * T) / (q.mw * Vg);
    const w = flow(P, T, mL)[0];
    out.t.push(time); out.P.push(P); out.T.push(T); out.Tw.push(Tw); out.mdot.push(w); out.m.push(m);
    minT = Math.min(minT, T); minTw = Math.min(minTw, Tw); minTd = Math.min(minTd, T - q.jt * Math.max(P - q.pBack, 0)); peak = Math.max(peak, w);
    if (q.pMark && tMark === null && P <= q.pMark) tMark = time - (dt * (q.pMark - P)) / (Pold - P || 1e-300);
    if (P <= pEnd) tEnd = time - (dt * (pEnd - P)) / (Pold - P || 1e-300);
  }
  return { ...out, dt, tEnd, tMark, minT, minTw, minTdown: minTd, peak, discharged: disc, flashed, unreleased: held, theta: thetaLast, liquidOut: mL0 - mL, mass: { initial: m0 + mL0, final: m + mL }, pFinal: P };
}

// ---- pigging --------------------------------------------------------------------------------------------------------
/**
 * Quasi-steady pig run along a line. o: { s[] (m, arc length at the nodes), z[], vm[] (mixture velocity behind the pig, m/s), holdup[], vsl[], rhoM[] (also the density of the fluid passing the bypass unless rhoBypass[] is given),
 *   rhoL, D, fric (Pa to keep the pig moving), mass (kg), bypass (bypass area / pipe area), cdBypass, leak (fraction of the swept liquid passing back through the pig),
 *   holdSlug (liquid fraction of the slug), qDrain (m³/s the receiving facility can process), fSlug (Darcy friction factor of the liquid slug),
 *   s0 (m: arc length at which a pig already in the line sits at the start; default the launcher), slug0 (m³ of liquid already ahead of it) }
 * The bypass discharge coefficient is Cd = 1/√K with K the pressure-loss coefficient on the bypass velocity (K = 1–1.5 is usual for plain bypass pigs).
 * Force balance on the pig: Δp A = friction + weight component; bypass leakage of the local mixture through the pig as an orifice under Δp; the liquid overtaken
 * collects as a slug ahead. Returns { t[], x[] (s of the pig), v[], slug[] (m³), transit (s), tFront (s), vMean, swept, leaked, received (m³), duration (s of liquid arrival),
 *   surge (m³ above the drain capacity), dpPig (Pa, mean), dpExtra (Pa, largest added line pressure drop), stalled }.
 */
export function pigRun(o) {
  const q = { fric: 1e5, mass: 80, bypass: 0.02, cdBypass: 0.9, leak: 0, holdSlug: 0.95, qDrain: 0.05, fSlug: 0.02, slug0: 0, ...o }, n = q.s.length, A = (Math.PI * q.D * q.D) / 4, L = q.s[n - 1], sStart = clamp(q.s0 ?? q.s[0], q.s[0], L - 1e-6);
  const out = { t: [0], x: [sStart], v: [], slug: [q.slug0] };
  let t = 0, Vs = q.slug0, swept = q.slug0, leaked = 0, received = 0, tFront = null, dpSum = 0, dpExtra = 0, stalled = false;
  for (let i = 0; i < n - 1; i++) {
    if (q.s[i + 1] <= sStart) continue;
    const ds = q.s[i + 1] - Math.max(q.s[i], sStart), sinT = clamp((q.z[i + 1] - q.z[i]) / (q.s[i + 1] - q.s[i]), -1, 1), dp = Math.max(q.fric + (q.mass * G * sinT) / A, 0.05 * q.fric), vLeak = q.cdBypass * q.bypass * Math.sqrt((2 * dp) / Math.max((q.rhoBypass || q.rhoM)[i], 0.5));
    const v = q.vm[i] - vLeak;
    if (!(v > 0.02)) { stalled = true; out.v.push(0); break; }
    const H = q.holdup[i], vL = H > 1e-6 ? Math.min(q.vsl[i] / H, v) : 0, dV = A * ds * H * (1 - vL / v), back = q.leak * dV;
    t += ds / v; swept += dV; leaked += back; Vs += dV - back; dpSum += dp * ds;
    const sPig = q.s[i + 1], Ls = Vs / (A * q.holdSlug);
    if (sPig + Ls >= L) { if (tFront === null) tFront = t - Math.min(ds, sPig + Ls - L) / v; const keep = A * q.holdSlug * (L - sPig); received += Vs - keep; Vs = keep; }
    // added pressure drop while the slug is in the line: liquid-full head and friction instead of the two-phase values
    let ex = 0; for (let j = i + 1; j < n - 1 && q.s[j] < sPig + Ls; j++) { const dz = q.z[j + 1] - q.z[j], dsj = q.s[j + 1] - q.s[j]; ex += (q.rhoL * q.holdSlug - q.rhoM[j]) * G * dz + ((q.fSlug * q.rhoL * q.vm[j] * q.vm[j]) / (2 * q.D)) * dsj * 0.5; }
    dpExtra = Math.max(dpExtra, ex + dp);
    out.t.push(t); out.x.push(sPig); out.v.push(v); out.slug.push(Vs + received);
  }
  out.v.push(out.v[out.v.length - 1] ?? 0);
  const duration = tFront === null ? 0 : Math.max(t - tFront, 1e-9);
  return { ...out, transit: stalled ? null : t, tFront, vMean: stalled ? 0 : (L - sStart) / t, swept, leaked, received, duration, surge: Math.max(0, received - q.qDrain * duration), dpPig: dpSum / Math.max(L - sStart, 1e-9), dpExtra, stalled, inPipe: Vs, start: sStart };
}

// ---- chemical injection ---------------------------------------------------------------------------------------------
/** Ogata–Banks solution of 1-D advection–dispersion for a step at the inlet: C/C0 at distance x, time t, velocity U, dispersion D. */
export function frontAnalytic(x, t, U, D) {
  if (t <= 0) return 0;
  const s = 2 * Math.sqrt(D * t), a = 0.5 * erfc((x - U * t) / s), pe = (U * x) / D;
  return a + (pe < 600 ? 0.5 * Math.exp(pe) * erfc((x + U * t) / s) : 0);
}
/**
 * Inhibitor front along the line: finite-volume solution on cells of equal transit time (advection is then an exact shift)
 * with implicit dispersion. o: { x[] (m, nodes), u[] (m/s liquid velocity at the nodes) | U (constant), L, D (m²/s), n (cells), tEnd (s), c0, cInit }
 * Returns { xc[] (cell centres), t[], c: [[…]] (rows = times), outlet[] (C at the last cell), tProtect (s until the whole line holds ≥ 95 % of c0 | null), transit (s) }.
 */
export function inhibitorFront(o) {
  const q = { D: 5, n: 120, c0: 1, cInit: 0, frac: 0.95, rows: 40, ...o }, L = q.L ?? q.x[q.x.length - 1], xs = q.x || [0, L], us = q.u || [q.U, q.U];
  const tt = [0]; for (let i = 1; i < xs.length; i++) tt.push(tt[i - 1] + (2 * (xs[i] - xs[i - 1])) / Math.max(us[i] + us[i - 1], 1e-6)); // transit time to each node
  const transit = tt[tt.length - 1], dtau = transit / q.n, edges = Array.from({ length: q.n + 1 }, (_, j) => interp1(tt, xs, j * dtau)), dx = edges.slice(1).map((e, j) => e - edges[j]), xc = dx.map((d, j) => edges[j] + 0.5 * d);
  const a = new Array(q.n), b = new Array(q.n), cc = new Array(q.n), steps = Math.max(1, Math.round((q.tEnd ?? 2 * transit) / dtau)), every = Math.max(1, Math.round(steps / q.rows));
  for (let j = 0; j < q.n; j++) { const gl = j > 0 ? (q.D * dtau) / (0.5 * (dx[j] + dx[j - 1])) : 0, gr = j < q.n - 1 ? (q.D * dtau) / (0.5 * (dx[j] + dx[j + 1])) : 0; a[j] = -gl; cc[j] = -gr; b[j] = dx[j] + gl + gr; }
  let c = new Array(q.n).fill(q.cInit), tProtect = null; const out = { xc, t: [0], c: [c.slice()], outlet: [c[q.n - 1]], tAll: [0] };
  for (let k = 1; k <= steps; k++) {
    for (let j = q.n - 1; j > 0; j--) c[j] = c[j - 1]; c[0] = q.c0;
    c = tridiag(a, b, cc, c.map((v, j) => v * dx[j]));
    const t = k * dtau, cmin = Math.min(...c); out.outlet.push(c[q.n - 1]); out.tAll.push(t);
    if (tProtect === null && cmin >= q.frac * q.c0) tProtect = t;
    if (k % every === 0 || k === steps) { out.t.push(t); out.c.push(c.slice()); }
  }
  return { ...out, tProtect, transit, dtau };
}
const PSAT_ETOH = (T) => 10 ** (8.20417 - 1642.89 / (T + 230.3)) * 133.322; // Antoine, mmHg and °C (not checked against a source: see PROVENANCE)
/**
 * Thermodynamic-inhibitor requirement with phase-partitioning losses. o: { dT (°C depression needed), inh ('MeOH' | 'MEG' | …), S (g/kg salinity),
 *   mWater (kg/s free water), lean (wt % purity of the injected chemical), P (bara), T (°C where the phases separate), qGasStd (Sm³/d), mOil (kg/s), rhoOil (kg/m³),
 *   eff (fraction of the injected chemical that reaches the aqueous phase before partitioning; injection efficiency) }
 * Partitioning: methanol to the gas by the K-value correlation methanolK (within its fitted range), 0.4 kg per m³ of liquid hydrocarbon;
 * glycols: gas loss neglected, 3.5 L per million Sm³ of gas to the hydrocarbon liquid; ethanol: modified Raoult's law (activity coefficient 1.6).
 * Returns { wt (wt % in the aqueous phase), mAq (kg/s inhibitor in the water), lossGas, lossOil (kg/s), mTotal (kg/s pure), qInject (m³/d of the lean chemical), rich (kg/s aqueous phase returned) }.
 */
export function inhibitorDose(o) {
  const q = { dT: 0, inh: 'MeOH', S: 0, mWater: 1, lean: 100, P: 25, T: 20, qGasStd: 0, mOil: 0, rhoOil: 800, eff: 1, ...o }, inh = INHIBITORS[q.inh] || INHIBITORS.MeOH, wt = q.dT > 0 ? inhibitorFor(q.dT, q.inh, q.S) : 0, w = wt / 100, lean = clamp(q.lean, 30, 100) / 100;
  if (w <= 0 || q.inh === 'none') return { wt: 0, mAq: 0, lossGas: 0, lossOil: 0, mTotal: 0, qInject: 0, rich: q.mWater, attainable: true };
  // water brought in by a lean (regenerated) chemical dilutes it: m_inh = w (mW + m_lean (1 − lean) + m_inh)
  const mAq = (w * q.mWater) / Math.max(1 - w / lean, 0.02), xAq = mAq / inh.MW / (mAq / inh.MW + (q.mWater + (mAq * (1 - lean)) / lean) / 18.015);
  const alcohol = q.inh === 'MeOH' || q.inh === 'EtOH', gasMol = q.qGasStd / DAY / VM_STD;
  const y = q.inh === 'MeOH' ? Math.min(methanolK(clamp(q.P, 6.9, 345), clamp(q.T, -23, 38)) * xAq, 0.2) : q.inh === 'EtOH' ? Math.min((1.6 * xAq * PSAT_ETOH(q.T)) / (q.P * 1e5), 0.2) : 0, lossGas = (y * gasMol * inh.MW) / 1000;
  const lossOil = alcohol ? (0.4 * q.mOil) / Math.max(q.rhoOil, 1) : 3.5e-9 * (q.qGasStd / DAY) * inh.rho, eff = clamp(q.eff, 0.05, 1), mTotal = (mAq + lossGas + lossOil) / eff;
  return { wt, mAq, lossGas, lossOil, mTotal, y, qInject: ((mTotal / lean) / (lean * inh.rho + (1 - lean) * 1000)) * DAY, rich: q.mWater + mAq / lean, attainable: wt < 93.9 && w / lean < 0.98 };
}

// ---- riser slugging: low-order four-state model ------------------------------------------------------------
/**
 * Four-state pipeline–riser model (gas and liquid mass in the feed pipeline and in the riser) with a low-point
 * orifice pair and a topside choke.
 * p: { D, Lp (feed length), Vp (feed volume), Lr (riser height), Vr (riser + topside volume), theta (rad, feed inclination at the low point),
 *      rhoL, mwG (kg/mol), Z, Tp, Tr (K), muL, wG, wL (kg/s inflow), Ps (Pa separator), Cv (choke), kH, kL, kG (optional), aLp (feed liquid fraction), rough,
 *      Dr (riser diameter when it differs from D), Lh (horizontal length between riser top and choke, added to the riser friction length), Kpc (valve constant in m², overrides Cv),
 *      chokeExp (valve characteristic f(z) = z^chokeExp; 1 = linear), fric: 'haaland' | 'dkm' (0.0056 + 0.5 Re^-0.32, as in the published model) }
 * Structure after Jahanshahi & Skogestad (2011): level at the low point h = kH hc ᾱL + (mL − ρL V ᾱL) sinθ / (A (1 − ᾱL) ρL), hc = D/cosθ, gas area A (1 − h/hc)²,
 * orifice equations for gas and liquid at the low point, top liquid fraction 2 ᾱLr − AL/A, valve w = Kpc f(z) √(ρt (Prt − Ps)).
 * Returns { p, alg(y, z, wG, wL, Ps), f(z)(t, y), steady(z, wG, wL), linearise(z), poles(z), critical() }.
 * y = [mGp, mLp, mGr, mLr] (kg), z = choke opening 0–1.
 */
export function slugModel(p) {
  const o = { D: 0.254, Lp: 5000, Lr: 300, theta: 0.02, rhoL: 800, mwG: 0.02, Z: 0.9, Tp: 320, Tr: 310, muL: 2e-3, wG: 1, wL: 20, Ps: 25e5, Cv: 400, kH: 0.7, kL: 0.3, aLp: 0.4, rough: 4.5e-5, chokeExp: 1, Lh: 0, fric: 'haaland', ...p };
  const r = o.D / 2, A = Math.PI * r * r, Dr = o.Dr || o.D, Ar = (Math.PI * Dr * Dr) / 4, th = Math.max(Math.abs(o.theta), 1e-3), hc = (2 * r) / Math.cos(th), sinT = Math.sin(th);
  o.Vp = o.Vp || A * o.Lp; o.Vr = o.Vr || Ar * (o.Lh > 0 ? o.Lr + o.Lh : o.Lr * 1.15);
  const RTp = (o.Z * R * o.Tp) / o.mwG, RTr = (o.Z * R * o.Tr) / o.mwG, aL = clamp(o.aLp, 0.05, 0.9), hbar = o.kH * hc * aL, Kc = o.Kpc > 0 ? o.Kpc : 2.403e-5 * o.Cv, fz = o.chokeExp === 1 ? (z) => z : (z) => Math.max(z, 0) ** o.chokeExp;
  const dmdh = (A * (1 - aL) * o.rhoL) / sinT, lamOf = (u, rho, D) => { const Re = Math.max((rho * u * D) / o.muL, 100); return o.fric === 'dkm' ? 0.0056 + 0.5 * Re ** -0.32 : frictionFactor(Re, o.rough / D, 'haaland'); };
  const areas = (h) => { const hh = clamp(h, 0, hc), AG = A * ((hc - hh) / hc) ** 2; return [AG, A - AG]; };
  // friction factors are frozen at the nominal superficial velocities (weak function of the state, large saving in the Jacobians)
  const rhoGref = (1.6 * o.Ps) / RTr, u0 = o.wL / (A * o.rhoL), um0 = o.wL / (Ar * o.rhoL) + o.wG / (rhoGref * Ar), lamP = lamOf(u0, o.rhoL, o.D), lamR = lamOf(um0, 0.5 * o.rhoL, Dr);
  const fricP = (wL) => { const u = wL / (A * o.rhoL); return (aL * lamP * o.rhoL * u * u * o.Lp) / (2 * o.D); };
  const fricR = (wG, wL, rhoM, aLr) => { const um = wL / (Ar * o.rhoL) + wG / (rhoGref * Ar); return (aLr * lamR * rhoM * um * um * (o.Lr + o.Lh)) / (2 * Dr); };
  if (!(o.kG > 0)) { // low-point gas coefficient from the nominal state: gas and liquid pass side by side at the mean level
    const [AG, AL] = areas(hbar), dPL = (o.wL / (o.kL * AL)) ** 2 / o.rhoL, dPG = Math.max(dPL - o.rhoL * G * hbar, 0.25 * dPL), rhoG = o.rhoGnom || 60;
    o.kG = o.wG / (AG * Math.sqrt(rhoG * dPG));
  }
  const alg = (y, z, wGin = o.wG, wLin = o.wL, Ps = o.Ps) => {
    const mGp = Math.max(y[0], 1e-6), mLp = y[1], mGr = Math.max(y[2], 1e-6), mLr = Math.max(y[3], 0);
    const rhoGp = mGp / Math.max(o.Vp - mLp / o.rhoL, 0.02 * o.Vp), Pp = rhoGp * RTp, h = hbar + (mLp - o.rhoL * o.Vp * aL) / dmdh;
    const VGr = Math.max(o.Vr - mLr / o.rhoL, 0.01 * o.Vr), rhoGr = mGr / VGr, aLr = clamp(mLr / (o.Vr * o.rhoL), 0, 1), rhoM = (mGr + mLr) / o.Vr;
    const Prt = rhoGr * RTr + 2e8 * Math.max(mLr / (o.Vr * o.rhoL) - 0.985, 0), Prb = Prt + rhoM * G * o.Lr + fricR(wGin, wLin, rhoM, aLr), [AG, AL] = areas(h);
    const dPG = Pp - fricP(wLin) - Prb, wGlp = dPG > 0 ? o.kG * AG * Math.sqrt(rhoGp * dPG) : 0;
    const dPL = dPG + o.rhoL * G * clamp(h, 0, 4 * hc), wLlp = dPL > 0 ? o.kL * AL * Math.sqrt(o.rhoL * dPL) : 0;
    const aLt = clamp(2 * aLr - AL / A, 0, 1), rhoT = aLt * o.rhoL + (1 - aLt) * rhoGr, xL = (aLt * o.rhoL) / Math.max(rhoT, 1e-9), dPc = Prt - Ps, w = dPc > 0 ? Kc * fz(z) * Math.sqrt(rhoT * dPc) : 0;
    return { Pp, Prt, Prb, h: h / hc, aLr, aLt, rhoT, w, wGout: (1 - xL) * w, wLout: xL * w, wGlp, wLlp, d: [wGin - wGlp, wLin - wLlp, wGlp - (1 - xL) * w, wLlp - xL * w] };
  };
  /** Equilibrium for a choke opening (exists for every opening: the level and the top pressure are found by 1-D root finding). */
  const steady = (z, wG = o.wG, wL = o.wL, Ps = o.Ps) => {
    const w = wG + wL, xL = wL / w, top = (Prt) => { const rg = Prt / RTr, aLt = (xL * rg) / (o.rhoL - xL * (o.rhoL - rg)); return { rg, aLt, rhoT: aLt * o.rhoL + (1 - aLt) * rg }; };
    const gT = (Prt) => Kc * fz(z) * Math.sqrt(top(Prt).rhoT * (Prt - Ps)) - w;
    let hi = Ps * 1.0001 + 10; while (gT(hi) < 0 && hi < 1e13) hi = Ps + (hi - Ps) * 2;
    const Prt = brent(gT, Ps, hi, 1e-10), t = top(Prt);
    const at = (h) => {
      const [AG, AL] = areas(h), aLr = clamp(0.5 * (t.aLt + AL / A), 0, 0.99), rhoM = aLr * o.rhoL + (1 - aLr) * t.rg, Prb = Prt + rhoM * G * o.Lr + fricR(wG, wL, rhoM, aLr);
      const dPG = (wL / (o.kL * AL)) ** 2 / o.rhoL - o.rhoL * G * h, Pp = Prb + fricP(wL) + dPG;
      return { AG, aLr, Prb, dPG, Pp, res: dPG > 0 ? o.kG * AG * Math.sqrt((Pp / RTp) * dPG) - wG : -wG };
    };
    const h = brent((x) => at(x).res, 1e-4 * hc, (1 - 1e-7) * hc, 1e-13), e = at(h), mLp = o.rhoL * o.Vp * aL + (h - hbar) * dmdh, mLr = e.aLr * o.Vr * o.rhoL;
    return { y: [(e.Pp / RTp) * (o.Vp - mLp / o.rhoL), mLp, t.rg * (o.Vr - mLr / o.rhoL), mLr], Pp: e.Pp, Prt, Prb: e.Prb, h: h / hc, aLr: e.aLr };
  };
  /** Jacobians at the equilibrium of opening z: dx/dt = A x + B u (u = opening), outputs [Pp, Prt] in Pa. */
  const linearise = (z) => {
    const s = steady(z), n = 4, Aj = zeros(n, n), C = zeros(2, n), out = (y, u) => { const a = alg(y, u); return [a.d, [a.Pp, a.Prt]]; };
    for (let j = 0; j < n; j++) { const e = Math.abs(s.y[j]) * 1e-6 + 1e-7, yp = s.y.slice(), ym = s.y.slice(); yp[j] += e; ym[j] -= e; const [dp, op] = out(yp, z), [dm, om] = out(ym, z); for (let i = 0; i < n; i++) Aj[i][j] = (dp[i] - dm[i]) / (2 * e); for (let i = 0; i < 2; i++) C[i][j] = (op[i] - om[i]) / (2 * e); }
    const e = 1e-5, [dp, op] = out(s.y, z + e), [dm, om] = out(s.y, z - e);
    return { A: Aj, B: dp.map((v, i) => [(v - dm[i]) / (2 * e)]), C, D: op.map((v, i) => (v - om[i]) / (2 * e)), ys: s.y, steady: s, z };
  };
  const poles = (z) => eig(linearise(z).A), growth = (z) => Math.max(...poles(z).map((e) => e[0]));
  /** Smallest opening at which the equilibrium loses stability (Hopf point), or null when stable up to fully open. */
  const critical = (zLo = 0.02, zHi = 1, n = 13) => {
    const zs = Array.from({ length: n }, (_, i) => zLo * (zHi / zLo) ** (i / (n - 1))), g = zs.map(growth), k = g.findIndex((v) => v > 0);
    if (k < 0) return null;
    if (k === 0) return zLo;
    return brent(growth, zs[k - 1], zs[k], 1e-4);
  };
  return { p: o, alg, f: (z) => (t, y) => alg(y, typeof z === 'function' ? z(t) : z).d, steady, linearise, poles, growth, critical, hc, A, dmdh, Kc, fz, zFloor: clamp(((o.wG + o.wL) / (Kc * Math.sqrt(o.rhoL * 300e5))) ** (1 / o.chokeExp), 0.02, 0.9) }; // zFloor: opening below which the choke alone would take more than about 300 bar
}

// ---- PID, FOPDT identification and tuning rules -----------------------------------------------------------------
/**
 * Discrete PID with derivative filter (on the measurement), back-calculation anti-windup, output limits, rate limit and
 * actuator dead time. c: { kc, ti (s, 0 = no integral), td (s), N (filter), dt, uMin, uMax, rate (per s), bias, dead (s), antiWindup }.
 * Returns { step(sp, pv, ff, track) -> u applied, state }.
 */
export function pidController(c) {
  const o = { kc: 1, ti: 0, td: 0, N: 10, dt: 1, uMin: 0, uMax: 1, rate: Infinity, bias: 0, dead: 0, antiWindup: true, ...c };
  const st = { I: 0, D: 0, pv: null, u: o.bias, sat: false, buf: new Array(Math.max(0, Math.round(o.dead / o.dt))).fill(o.bias) };
  const step = (sp, pv, ff = 0, track = null) => {
    const e = sp - pv, a = o.td > 0 ? o.td / (o.td + o.N * o.dt) : 0;
    if (st.pv !== null && o.td > 0) st.D = a * st.D - ((o.kc * o.td * o.N) / (o.td + o.N * o.dt)) * (pv - st.pv);
    st.pv = pv;
    const raw = o.bias + o.kc * e + st.I + st.D + ff;
    let u = clamp(raw, o.uMin, o.uMax);
    if (Number.isFinite(o.rate)) u = clamp(u, st.u - o.rate * o.dt, st.u + o.rate * o.dt);
    if (track !== null) u = Math.min(u, track); // override / selector output takes over
    st.sat = u !== raw;
    if (o.ti > 0) { st.I += (o.kc * e * o.dt) / o.ti; if (o.antiWindup && st.sat) st.I += ((u - raw) * o.dt) / Math.max(o.ti * 0.5, o.dt); }
    st.u = u;
    if (!st.buf.length) return u;
    st.buf.push(u); return st.buf.shift();
  };
  return { step, state: st, cfg: o };
}
/**
 * Closed-loop response of a first-order-plus-dead-time process K e^(-θs)/(τs+1) under PID (exact discretisation of the process).
 * o: { K, tau, theta, kc, ti, td, dt, tEnd, sp (step size), dist(t), uMin, uMax, rate, antiWindup, N }.
 * Returns { t, y, u, iae, overshoot (fraction), maxI }.
 */
export function pidLoop(o) {
  const q = { K: 1, tau: 10, theta: 0, kc: 1, ti: 0, td: 0, dt: 0.05, tEnd: 100, sp: 1, uMin: -1e9, uMax: 1e9, ...o }, n = Math.round(q.tEnd / q.dt), a = Math.exp(-q.dt / q.tau), nd = Math.round(q.theta / q.dt);
  const pid = pidController({ ...q, bias: 0, dead: 0 }), buf = new Array(nd).fill(0), t = [0], y = [0], u = [];
  let x = 0, iae = 0, peak = 0, maxI = 0;
  for (let k = 0; k < n; k++) {
    const sp = typeof q.sp === 'function' ? q.sp(k * q.dt) : q.sp, uk = pid.step(sp, x); u.push(uk); buf.push(uk);
    const ud = buf.shift() + (q.dist ? q.dist(k * q.dt) : 0);
    x = a * x + q.K * (1 - a) * ud; t.push((k + 1) * q.dt); y.push(x); iae += Math.abs(sp - x) * q.dt; peak = Math.max(peak, x); maxI = Math.max(maxI, Math.abs(pid.state.I));
  }
  u.push(u[u.length - 1]);
  const spEnd = typeof q.sp === 'function' ? q.sp(q.tEnd) : q.sp;
  return { t, y, u, iae, overshoot: spEnd ? Math.max(0, peak / spEnd - 1) : 0, maxI };
}
/** First-order-plus-dead-time fit of a step response (two-point start, least-squares refinement). du = input step. Returns { K, tau, theta, rmse }. */
export function identifyFOPDT(t, y, du = 1) {
  const y0 = y[0], yE = mean(y.slice(-Math.max(3, Math.floor(y.length / 20)))), K = (yE - y0) / du, frac = y.map((v) => (v - y0) / (yE - y0 || 1e-300));
  const cross = (f) => { for (let i = 1; i < t.length; i++) if (frac[i] >= f && frac[i - 1] < f) return t[i - 1] + ((f - frac[i - 1]) * (t[i] - t[i - 1])) / (frac[i] - frac[i - 1]); return t[t.length - 1]; };
  const t28 = cross(0.283) - t[0], t63 = cross(0.632) - t[0], tau0 = Math.max(1.5 * (t63 - t28), 1e-6), th0 = Math.max(t63 - tau0, 0);
  const model = (tau, th) => t.map((ti) => (ti - t[0] <= th ? 0 : 1 - Math.exp(-(ti - t[0] - th) / tau))), sse = (q) => { const m = model(Math.max(q[0], 1e-9), Math.max(q[1], 0)); let s = 0; for (let i = 0; i < m.length; i++) s += (m[i] - frac[i]) ** 2; return s; };
  const fit = nelderMead(sse, [tau0, th0], { lo: [tau0 * 0.05, 0], hi: [tau0 * 20 + 1e-9, t63 + 1e-9], maxIter: 150, tol: 1e-10 });
  return { K, tau: fit.x[0], theta: fit.x[1], rmse: Math.sqrt(fit.f / t.length) * Math.abs(yE - y0) };
}
/**
 * Controller settings from a first-order-plus-dead-time model: SIMC (Skogestad 2003: Kc = τ/(K (τc + θ)), τI = min(τ, 4 (τc + θ)); τc = closed-loop
 * time constant, default θ), the SIMC PID extension for a first-order process (Grimholt & Skogestad 2013: series form with τD = θ/3 and τc = θ/2 by
 * default, converted here to the ideal form used by pidController) and Ziegler–Nichols (ideal form) from the ultimate gain and period of the same model.
 * Returns rows { rule, mode, kc, ti, td } and { ku, pu }.
 */
export function tuningRules({ K, tau, theta }, tauC) {
  const th = Math.max(theta, 1e-9), tc = tauC ?? th, rows = [];
  rows.push({ rule: 'SIMC', mode: 'P', kc: tau / (K * (tc + th)), ti: 0, td: 0 });
  rows.push({ rule: 'SIMC', mode: 'PI', kc: tau / (K * (tc + th)), ti: Math.min(tau, 4 * (tc + th)), td: 0 });
  { const tcD = tauC ?? th / 2, kcS = tau / (K * (tcD + th)), tiS = Math.min(tau, 4 * (tcD + th)), tdS = th / 3, f = 1 + tdS / tiS; // series → ideal: Kc f, τI f, τD / f
    rows.push({ rule: 'SIMC', mode: 'PID', kc: kcS * f, ti: tiS * f, td: tdS / f, series: { kc: kcS, ti: tiS, td: tdS } }); }
  // ultimate point: −ωθ − atan(ωτ) = −π
  const wu = brent((w) => w * th + Math.atan(w * tau) - Math.PI, 1e-9 / th, Math.PI / th, 1e-12), ku = Math.sqrt(1 + (wu * tau) ** 2) / Math.abs(K), pu = (2 * Math.PI) / wu, sg = Math.sign(K) || 1;
  rows.push({ rule: 'Ziegler–Nichols', mode: 'P', kc: 0.5 * ku * sg, ti: 0, td: 0 });
  rows.push({ rule: 'Ziegler–Nichols', mode: 'PI', kc: 0.45 * ku * sg, ti: pu / 1.2, td: 0 });
  rows.push({ rule: 'Ziegler–Nichols', mode: 'PID', kc: 0.6 * ku * sg, ti: pu / 2, td: pu / 8 });
  return { rows, ku, pu };
}
/** Gain and phase margins of the loop PID × FOPDT from a frequency scan. Returns { gm (factor), pm (deg), wc, w180 }. */
export function loopMargins({ K, tau, theta }, { kc, ti, td = 0, N = 10 }) {
  const L = (w) => {
    // C(jw) = kc (1 + 1/(jw ti) + jw td / (1 + jw td/N)), G(jw) = K e^{-jwθ} / (1 + jwτ)
    const dre = td > 0 ? ((w * td) * (w * td / N)) / (1 + (w * td / N) ** 2) : 0, dim = td > 0 ? (w * td) / (1 + (w * td / N) ** 2) : 0, cre = kc * (1 + dre), cim = kc * ((ti > 0 ? -1 / (w * ti) : 0) + dim);
    const mag = (Math.hypot(cre, cim) * Math.abs(K)) / Math.sqrt(1 + (w * tau) ** 2), ph = Math.atan2(cim * Math.sign(K * kc || 1), cre * Math.sign(K * kc || 1)) - Math.atan(w * tau) - w * theta;
    return { mag, ph };
  };
  const ws = Array.from({ length: 900 }, (_, i) => 10 ** (-5 + (10 * i) / 899) / Math.max(tau, 1e-9));
  let gm = null, pm = null, wc = null, w180 = null;
  for (let i = 1; i < ws.length; i++) {
    const a = L(ws[i - 1]), b = L(ws[i]);
    if (pm === null && a.mag >= 1 && b.mag < 1) { const f = (a.mag - 1) / (a.mag - b.mag); wc = ws[i - 1] + f * (ws[i] - ws[i - 1]); pm = 180 + ((a.ph + f * (b.ph - a.ph)) * 180) / Math.PI; }
    if (gm === null && a.ph > -Math.PI && b.ph <= -Math.PI) { const f = (a.ph + Math.PI) / (a.ph - b.ph); w180 = ws[i - 1] + f * (ws[i] - ws[i - 1]); gm = 1 / (a.mag + f * (b.mag - a.mag)); }
  }
  return { gm, pm, wc, w180 };
}

// ---- linear MPC, LQ and Kalman filtering ------------------------------------------------------------------------
/** Box-constrained convex QP min ½ uᵀH u + fᵀu, lo ≤ u ≤ hi, by accelerated projected gradient. Returns { u, iterations }. */
export function boxQP(H, f, lo, hi, u0, { maxIter = 400, tol = 1e-9 } = {}) {
  const n = f.length; let v = new Array(n).fill(1), L = 1;
  for (let k = 0; k < 30; k++) { const w = mv(H, v), nw = Math.sqrt(dot(w, w)) || 1; L = nw; v = w.map((x) => x / nw); } // power iteration for the Lipschitz constant
  const proj = (u) => u.map((x, i) => clamp(x, lo[i], hi[i]));
  let u = proj(u0 || new Array(n).fill(0)), yk = u.slice(), tk = 1, it = 0;
  for (; it < maxIter; it++) {
    const g = mv(H, yk).map((x, i) => x + f[i]), un = proj(yk.map((x, i) => x - g[i] / (1.02 * L))), tn = 0.5 * (1 + Math.sqrt(1 + 4 * tk * tk)), ch = Math.max(...un.map((x, i) => Math.abs(x - u[i])));
    yk = un.map((x, i) => x + ((tk - 1) / tn) * (x - u[i])); u = un; tk = tn;
    if (ch < tol) break;
  }
  return { u, iterations: it };
}
/**
 * Linear MPC move for x⁺ = A x + B u, y = C x (single input, single output): minimises Σ q (y − r)² + rDu Δu² over the
 * prediction horizon np with nc free moves and bounds uMin ≤ u ≤ uMax and |Δu| ≤ duMax on every move (QP by projected gradient).
 * Returns { u (first move), seq, unconstrained (sequence without bounds), active (bool), Phi, H }.
 */
export function mpc({ A, B, C, x, uPrev = 0, r = 0, np = 20, nc = 5, q = 1, rDu = 1, uMin = -Infinity, uMax = Infinity, duMax = Infinity, cache }) {
  let m = cache && cache.H ? cache : null;
  if (!m) {
    const n = A.length, b = B.map((v) => (Array.isArray(v) ? v[0] : v)), F = [], Phi = zeros(np, nc);
    let Ak = eye(n); const imp = []; // impulse responses C A^k B
    for (let k = 0; k < np; k++) { imp.push(dot(C, mv(Ak, b))); Ak = mm(Ak, A); F.push(mv(tr(Ak), C)); } // F[k] = C A^(k+1)
    // u is parameterised by its moves Δu_0..Δu_{nc-1}; y_{k+1} = C A^{k+1} x + Σ_j step(k − j) Δu_j + step(k) uPrev
    const stp = []; let s = 0; for (let k = 0; k < np; k++) { s += imp[k]; stp.push(s); }
    for (let k = 0; k < np; k++) for (let j = 0; j < nc; j++) Phi[k][j] = k >= j ? stp[k - j] : 0;
    const H = zeros(nc, nc); for (let i = 0; i < nc; i++) for (let j = 0; j < nc; j++) { let v = 0; for (let k = 0; k < np; k++) v += Phi[k][i] * Phi[k][j]; H[i][j] = 2 * q * v + (i === j ? 2 * rDu : 0); }
    m = { F, Phi, H, stp }; if (cache) Object.assign(cache, m);
  }
  const free = m.F.map((Fk, k) => dot(Fk, x) + m.stp[k] * uPrev - (Array.isArray(r) ? r[Math.min(k, r.length - 1)] : r)), f = new Array(nc).fill(0);
  for (let j = 0; j < nc; j++) { let v = 0; for (let k = 0; k < np; k++) v += m.Phi[k][j] * free[k]; f[j] = 2 * q * v; }
  let unc; try { unc = solveLinear(m.H, f.map((v) => -v)); } catch { unc = f.map((v, j) => -v / (m.H[j][j] || 1)); }
  if (!unc.every(Number.isFinite)) unc = new Array(nc).fill(0);
  // bounds on the moves: rate limit directly, absolute limits through the running sum (tightened so every partial sum stays inside)
  const within = (d) => { let u = uPrev; for (const x of d) { if (Math.abs(x) > duMax + 1e-12) return false; u += x; if (u < uMin - 1e-12 || u > uMax + 1e-12) return false; } return true; };
  let seq = unc, active = false;
  if (!within(unc)) {
    active = true;
    // alternate projected-gradient solves with the box updated from the running input (exact for a single free move, tight in practice)
    let lo = new Array(nc).fill(-duMax), hi = new Array(nc).fill(duMax), d = unc.map((x) => clamp(x, -duMax, duMax));
    for (let pass = 0; pass < 6; pass++) {
      let u = uPrev;
      for (let j = 0; j < nc; j++) { lo[j] = Math.max(-duMax, uMin - u); hi[j] = Math.min(duMax, uMax - u); if (lo[j] > hi[j]) lo[j] = hi[j] = clamp(0, uMin - u, uMax - u); u += clamp(d[j], lo[j], hi[j]); }
      const nd = boxQP(m.H, f, lo, hi, d, { maxIter: 300 }).u, ch = Math.max(...nd.map((x, i) => Math.abs(x - d[i]))); d = nd;
      if (ch < 1e-10) break;
    }
    seq = d;
  }
  return { u: uPrev + seq[0], seq, unconstrained: unc, active };
}
/** Finite-horizon discrete LQ regulator by the Riccati recursion: cost Σ_{k=1..N} xᵀQx + Σ_{k=0..N-1} uᵀRu. Returns { K0 (first gain, u = −K0 x), P }. */
export function lqFinite(A, B, Q, Rm, N) {
  let P = Q.map((r) => r.slice()), K = null; const At = tr(A), Bt = tr(B);
  for (let k = N - 1; k >= 0; k--) {
    const S = madd(Rm, mm(mm(Bt, P), B)); K = mm(inv(S), mm(mm(Bt, P), A));
    if (k > 0) P = madd(Q, madd(mm(mm(At, P), A), mm(mm(mm(At, P), B), K), -1));
  }
  return { K0: K, P };
}
/**
 * Linear Kalman filter over a record. m: { A, B, C, Q, R, x0, P0 } (matrices as arrays of rows), u: inputs per step (arrays), y: measurements per step (arrays).
 * Returns { x: [[…]], K (last gain), P (last a-posteriori covariance), Pm (last a-priori covariance) }.
 */
export function kalman(m, u, y) {
  const n = m.A.length, At = tr(m.A), Ct = tr(m.C); let x = m.x0.slice(), P = m.P0.map((r) => r.slice()), K = null, Pm = P; const xs = [];
  for (let k = 0; k < y.length; k++) {
    const xm = mv(m.A, x).map((v, i) => v + (m.B && u ? dot(m.B[i], u[k]) : 0)); Pm = madd(mm(mm(m.A, P), At), m.Q);
    const S = madd(mm(mm(m.C, Pm), Ct), m.R); K = mm(mm(Pm, Ct), inv(S));
    const innov = y[k].map((v, i) => v - dot(m.C[i], xm)); x = xm.map((v, i) => v + dot(K[i], innov));
    const IKC = madd(eye(n), mm(K, m.C), -1); P = madd(mm(mm(IKC, Pm), tr(IKC)), mm(mm(K, m.R), tr(K))); // Joseph form
    xs.push(x.slice());
  }
  return { x: xs, K, P, Pm };
}
/**
 * Extended Kalman filter for x⁺ = F(x, u), y = h(x, u) with numerical Jacobians. o: { F, h, x0, P0, Q, R, u[], y[][] , scale[] (state scales for the finite differences) }.
 * Returns { x: [[…]], P }.
 */
export function ekf({ F, h, x0, P0, Q, R: Rn, u, y, scale }) {
  const n = x0.length; let x = x0.slice(), P = P0.map((r) => r.slice()); const xs = [];
  const jac = (fn, x) => { const f0 = fn(x), J = zeros(f0.length, n); for (let j = 0; j < n; j++) { const e = 1e-5 * (scale ? scale[j] : Math.abs(x[j]) + 1e-6), xp = x.slice(); xp[j] += e; const fp = fn(xp); for (let i = 0; i < f0.length; i++) J[i][j] = (fp[i] - f0[i]) / e; } return { f0, J }; };
  for (let k = 0; k < y.length; k++) {
    const fx = F(x, u ? u[k] : 0), pr = Array.isArray(fx) ? jac((s) => F(s, u ? u[k] : 0), x) : { f0: fx.x, J: fx.J }, Pm = madd(mm(mm(pr.J, P), tr(pr.J)), Q), me = jac((s) => h(s, u ? u[k] : 0), pr.f0), Ht = tr(me.J), S = madd(mm(mm(me.J, Pm), Ht), Rn), K = mm(mm(Pm, Ht), inv(S));
    const innov = y[k].map((v, i) => v - me.f0[i]); x = pr.f0.map((v, i) => v + dot(K[i], innov));
    const IKH = madd(eye(n), mm(K, me.J), -1); P = madd(mm(mm(IKH, Pm), tr(IKH)), mm(mm(K, Rn), tr(K))); P = P.map((r, i) => r.map((v, j) => 0.5 * (v + P[j][i])));
    xs.push(x.slice());
  }
  return { x: xs, P };
}

// ---- stiff integration and the closed-loop slugging / separator simulation ----------------------------------
/** One step of the L-stable two-stage Rosenbrock scheme ROS2 with a finite-difference Jacobian. Returns { y, err }. */
export function ros2Step(f, t, y, h) {
  const n = y.length, g = 1 + Math.SQRT1_2, f0 = f(t, y), M = new Array(n), A = new Array(n), yp = y.slice(), piv = [];
  for (let i = 0; i < n; i++) { M[i] = new Array(n); A[i] = new Array(n); }
  for (let j = 0; j < n; j++) { const e = 1e-7 * (Math.abs(y[j]) + 1e-3); yp[j] = y[j] + e; const fp = f(t, yp); yp[j] = y[j]; for (let i = 0; i < n; i++) A[i][j] = M[i][j] = (i === j ? 1 : 0) - (g * h * (fp[i] - f0[i])) / e; }
  // one LU factorisation (partial pivoting) serves both stages
  for (let k = 0; k < n; k++) { let p = k; for (let i = k + 1; i < n; i++) if (Math.abs(A[i][k]) > Math.abs(A[p][k])) p = i; piv.push(p); if (p !== k) { const tmp = A[k]; A[k] = A[p]; A[p] = tmp; } const d = A[k][k] || 1e-300; for (let i = k + 1; i < n; i++) { const l = (A[i][k] /= d); if (l !== 0) for (let j = k + 1; j < n; j++) A[i][j] -= l * A[k][j]; } }
  const lu = (b) => { const x = b.slice(); for (let k = 0; k < n; k++) { const p = piv[k]; if (p !== k) { const tmp = x[k]; x[k] = x[p]; x[p] = tmp; } } for (let k = 0; k < n; k++) for (let i = k + 1; i < n; i++) x[i] -= A[i][k] * x[k]; for (let i = n - 1; i >= 0; i--) { let s = x[i]; for (let j = i + 1; j < n; j++) s -= A[i][j] * x[j]; x[i] = s / (A[i][i] || 1e-300); } return x; };
  const k1 = lu(f0), f1 = f(t + h, y.map((v, i) => v + h * k1[i])), k2 = lu(f1.map((v, i) => v - 2 * k1[i]));
  return { y: y.map((v, i) => v + 1.5 * h * k1[i] + 0.5 * h * k2[i]), err: k1.map((v, i) => 0.5 * h * (v + k2[i])), M };
}
/** Adaptive ROS2 integration for stiff systems. Returns { t: [], y: [[]] }. opt: { rtol, atol, hInit, hMax, maxSteps }. */
export function integrateStiff(f, y0, t0, t1, { rtol = 1e-4, atol = 1e-6, hInit, hMax = Infinity, maxSteps = 20000 } = {}) {
  let t = t0, y = y0.slice(), h = Math.min(hInit ?? (t1 - t0) / 200, hMax); const ts = [t], ys = [y.slice()];
  for (let s = 0; s < maxSteps && t < t1 - 1e-12 * Math.abs(t1); s++) {
    if (t + h > t1) h = t1 - t;
    const st = ros2Step(f, t, y, h); let e = 0;
    for (let i = 0; i < y.length; i++) e = Math.max(e, Math.abs(st.err[i]) / (atol + rtol * Math.max(Math.abs(y[i]), Math.abs(st.y[i]))));
    if ((e <= 1 && st.y.every(Number.isFinite)) || h < 1e-9) { t += h; y = st.y; ts.push(t); ys.push(y.slice()); }
    h = Math.min(hMax, h * clamp(0.9 / Math.sqrt(Math.max(e, 1e-10)), 0.2, 4));
  }
  return { t: ts, y: ys };
}
/**
 * Closed-loop simulation of the slugging model with its topside choke and a receiving separator (liquid level and
 * pressure loops with first-order valve actuators). The run starts at the equilibrium of zTarget (slightly perturbed) under
 * the selected controller — 'pid' (inlet pressure → choke, optional cascade through a flow loop, high-level override),
 * 'mpc' (linear MPC on the model linearised at zTarget with a Kalman filter on the noisy inlet and topside pressures) or
 * 'open' — and the controller is switched to manual at tOff so that the open-loop behaviour of the same operating point shows.
 * c: { tEnd, dt, zTarget, tOff, tStep, dSp (bar), kc (opening per bar), ti, td, rate (1/s), dead (s), zMin, zMax, noise (bar), cascade, override, feedForward,
 *      sepV (m³), qDrain (m³/s), levelSp, levelHi, tauValve (s), np, nc, tsMpc, qY, rDu, y0, wOf(t) -> [wG, wL],
 *      tauSensor (s, first-order lag of the inlet-pressure transmitter), schedule(z) -> factor on kc (gain scheduling on the choke opening) }
 * Returns the time series and { ampClosed, ampOpen (bar peak-to-peak), suppressed, meanP, maxLevel, surge (m³), carryOver (m³), iae }.
 */
export function slugControl(sm, c) {
  const o = { tEnd: 12 * HOUR, dt: 30, zTarget: 0.3, tOff: 6 * HOUR, tStep: 2 * HOUR, dSp: 0, mode: 'pid', kc: -0.1, ti: 1800, td: 0, N: 8, rate: 0.01, dead: 0, zMin: 0.02, zMax: 1, noise: 0.05, seed: 11, cascade: false, override: true, feedForward: false, antiWindup: true,
    sepV: 60, qDrain: 0.08, levelSp: 0.5, levelHi: 0.8, tauValve: 8, np: 20, nc: 4, tsMpc: 120, qY: 1, rDu: 400, perturb: 1.0005, ...c };
  const p = sm.p, eq = sm.steady(o.zTarget), rand = rng(o.seed), wL0 = p.wL, wG0 = p.wG, RTs = (p.Z * R * p.Tr) / p.mwG;
  const xl0 = clamp(wL0 / p.rhoL / o.qDrain, 0.02, 0.95), wGmax = wG0 / 0.5, vg0 = o.sepV * (1 - o.levelSp);
  // separator: y[4] liquid volume, y[5] gas mass, y[6] liquid-valve position, y[7] gas-valve position
  const psep = (y) => (Math.max(y[5], 1e-6) * RTs) / Math.max(o.sepV - y[4], 0.03 * o.sepV), over = (y) => Math.max(0, y[4] - 0.95 * o.sepV) * 0.2;
  const cmd = { z: o.zTarget, xl: xl0, xg: 0.5 }, wOf = o.wOf || (() => [wG0, wL0]);
  const f = (t, y) => { const [wg, wl] = wOf(t), Ps = psep(y), a = sm.alg(y, cmd.z, wg, wl, Ps); return [a.d[0], a.d[1], a.d[2], a.d[3], a.wLout / p.rhoL - y[6] * o.qDrain * Math.min(1, y[4] / (0.02 * o.sepV)) - over(y), a.wGout - (y[7] * wGmax * Ps) / p.Ps - (20 * wGmax * Math.max(0, Ps - 1.25 * p.Ps)) / p.Ps, (cmd.xl - y[6]) / o.tauValve, (cmd.xg - y[7]) / o.tauValve]; };
  let y = o.y0 ? o.y0.slice() : [eq.y[0], eq.y[1] * o.perturb, eq.y[2], eq.y[3], o.sepV * o.levelSp, (p.Ps / RTs) * vg0, xl0, 0.5];
  const pid = pidController({ kc: o.kc, ti: o.ti, td: o.td, N: o.N, dt: o.dt, uMin: o.zMin, uMax: o.zMax, rate: o.rate, dead: o.dead, bias: o.zTarget, antiWindup: o.antiWindup }), w0 = wG0 + wL0;
  const flowPid = pidController({ kc: (0.6 * o.zTarget) / w0, ti: 3 * o.dt, dt: o.dt, uMin: o.zMin, uMax: o.zMax, rate: o.rate, bias: o.zTarget }), master = pidController({ kc: (o.kc * w0) / Math.max(o.zTarget, 0.02), ti: o.ti, dt: o.dt, uMin: 0, uMax: 4 * w0, bias: w0, antiWindup: o.antiWindup });
  const lvl = pidController({ kc: -(o.kcLevel ?? 2), ti: o.tiLevel ?? 900, dt: o.dt, uMin: 0, uMax: 1, bias: xl0 }), prs = pidController({ kc: -(o.kcPress ?? 4), ti: 120, dt: o.dt, uMin: 0, uMax: 1, bias: 0.5 });
  let mp = null; // MPC model: scaled deviation states, outputs in bar, input = opening deviation
  if (o.mode === 'mpc') {
    const lin = sm.linearise(o.zTarget), sc = lin.ys.map((v) => Math.abs(v) || 1), As = lin.A.map((r, i) => r.map((v, j) => (v * sc[j]) / sc[i])), Bs = lin.B.map((r, i) => [r[0] / sc[i]]), Cs = lin.C.map((r) => r.map((v, j) => (v * sc[j]) / 1e5)), d = c2d(As, Bs, o.tsMpc);
    const gr = Math.max(...eig(lin.A).map((e) => e[0])); if (gr > 0) o.np = Math.max(3, Math.min(o.np, Math.ceil(6 / (gr * o.tsMpc)))); o.nc = Math.min(o.nc, o.np); // keep the prediction of an unstable model well conditioned
    mp = { lin, sc, Ad: d.Ad, Bd: d.Bd, C: Cs, Dy: lin.D.map((v) => v / 1e5), x: y.slice(0, 4).map((v, i) => (v - lin.ys[i]) / sc[i]), P: eye(4).map((r) => r.map((v) => v * 1e-4)), cache: {}, every: Math.max(1, Math.round(o.tsMpc / o.dt)), u: 0, yacc: [0, 0], nacc: 0, active: 0, moves: 0 };
  }
  const n = Math.round(o.tEnd / o.dt), out = { t: [], pIn: [], pBase: [], pTop: [], z: [], wL: [], wG: [], level: [], pSep: [], mLr: [], mLrEst: [], qLout: [], sp: [] };
  let ovrCount = 0, maxLevel = 0, minLevel = 1, carry = 0, iae = 0, started = false, pf = null, maxLevCl = 0, maxLevOp = 0, carryCl = 0;
  for (let k = 0; k <= n; k++) {
    const t = k * o.dt, Ps = psep(y), a = sm.alg(y, cmd.z, ...wOf(t), Ps), level = y[4] / o.sepV, pRaw = a.Pp / 1e5 + (o.bias || 0) + o.noise * rand.normal(), ptm = a.Prt / 1e5 + o.noise * rand.normal();
    pf = pf === null || !(o.tauSensor > 0) ? pRaw : pf + (o.dt / (o.tauSensor + o.dt)) * (pRaw - pf); const pm = pf;
    if (o.schedule) { const g = clamp(o.schedule(cmd.z), 0.2, 5); pid.cfg.kc = o.kc * g; master.cfg.kc = ((o.kc * w0) / Math.max(o.zTarget, 0.02)) * g; }
    const sp = (o.sp ?? eq.Pp) / 1e5 + (t >= o.tStep ? o.dSp : 0), auto = t < o.tOff && o.mode !== 'open';
    out.t.push(t); out.pIn.push(a.Pp / 1e5); out.pTop.push(a.Prt / 1e5); out.pBase.push(a.Prb / 1e5); out.z.push(cmd.z); out.wL.push(a.wLout); out.wG.push(a.wGout); out.level.push(level); out.pSep.push(Ps / 1e5); out.mLr.push(y[3]); out.qLout.push(y[6] * o.qDrain); out.sp.push(sp);
    maxLevel = Math.max(maxLevel, level); minLevel = Math.min(minLevel, level); carry += over(y) * o.dt; if (auto) iae += Math.abs(sp - a.Pp / 1e5) * o.dt;
    if (auto) { maxLevCl = Math.max(maxLevCl, level); carryCl += over(y) * o.dt; } else maxLevOp = Math.max(maxLevOp, level);
    // separator loops (always in automatic; direct acting, so the gains are negative)
    cmd.xl = lvl.step(o.levelSp, level, o.feedForward ? a.wLout / p.rhoL / o.qDrain - xl0 : 0); cmd.xg = prs.step(1, Ps / p.Ps);
    if (auto) {
      const zOvr = o.override ? clamp(o.zMax - (o.overrideGain ?? 6) * (level - o.levelHi), o.zMin, o.zMax) : null; // high-level override on the choke (low select)
      if (o.override && level > o.levelHi) ovrCount++;
      if (o.mode === 'pid') {
        if (!started) { master.cfg.bias = a.w; started = true; }
        cmd.z = o.cascade ? flowPid.step(master.step(sp, pm), a.w, 0, zOvr) : pid.step(sp, pm, 0, zOvr);
      } else {
        mp.yacc[0] += pm; mp.yacc[1] += ptm; mp.nacc++;
        if (k % mp.every === 0) {
          const ym = [mp.yacc[0] / mp.nacc - eq.Pp / 1e5 - mp.Dy[0] * mp.u, mp.yacc[1] / mp.nacc - eq.Prt / 1e5 - mp.Dy[1] * mp.u]; mp.yacc = [0, 0]; mp.nacc = 0;
          const kf = kalman({ A: mp.Ad, B: mp.Bd, C: mp.C, Q: eye(4).map((r) => r.map((v) => v * 1e-8)), R: [[o.noise ** 2 + 1e-4, 0], [0, o.noise ** 2 + 1e-4]], x0: mp.xPrior || mp.x, P0: mp.P }, null, [ym]);
          mp.x = kf.x[0]; mp.P = kf.P;
          const mv1 = mpc({ A: mp.Ad, B: mp.Bd, C: mp.C[0], x: mp.x, uPrev: mp.u, r: sp - eq.Pp / 1e5, np: o.np, nc: o.nc, q: o.qY, rDu: o.rDu, uMin: o.zMin - o.zTarget, uMax: (zOvr ?? o.zMax) - o.zTarget, duMax: o.rate * o.tsMpc, cache: mp.cache });
          mp.u = mv1.u; mp.moves++; if (mv1.active) mp.active++; cmd.z = clamp(o.zTarget + mp.u, o.zMin, o.zMax);
          mp.xPrior = mv(mp.Ad, mp.x).map((v, i) => v + mp.Bd[i][0] * mp.u); mp.P = madd(mm(mm(mp.Ad, mp.P), tr(mp.Ad)), eye(4).map((r) => r.map((v) => v * 1e-8)));
        }
      }
    }
    out.mLrEst.push(mp ? mp.lin.ys[3] + mp.x[3] * mp.sc[3] : y[3]);
    if (k < n) { let st = ros2Step(f, t, y, o.dt); if (!st.y.every(Number.isFinite)) { const h = o.dt / 8; st = { y }; for (let j = 0; j < 8; j++) st = ros2Step(f, t + j * h, st.y, h); } y = st.y.map((v, i) => (i < 6 ? Math.max(v, 1e-9) : clamp(v, 0, 1))); }
  }
  const iOff = Math.min(n, Math.round(o.tOff / o.dt)), amp = (a, s, e) => { let lo = Infinity, hi = -Infinity; for (let i = s; i < e; i++) { lo = Math.min(lo, a[i]); hi = Math.max(hi, a[i]); } return e > s ? hi - lo : 0; };
  const cl0 = Math.floor(0.6 * iOff), ampCl = o.mode === 'open' ? amp(out.pIn, Math.floor(0.5 * n), n + 1) : amp(out.pIn, cl0, iOff), ampOl = iOff < n ? amp(out.pIn, iOff + Math.floor(0.3 * (n - iOff)), n + 1) : ampCl;
  return { ...out, yEnd: y, spBar: (o.sp ?? eq.Pp) / 1e5, ampOpen: ampOl, ampClosed: ampCl, meanP: mean(out.pIn.slice(cl0, Math.max(iOff, cl0 + 1))), meanZ: mean(out.z.slice(cl0, Math.max(iOff, cl0 + 1))), suppressed: Number.isFinite(ampCl) && ampCl < Math.max(1, 6 * o.noise + 0.02 * (o.sp ?? eq.Pp) / 1e5), maxLevel, minLevel, maxLevelClosed: maxLevCl, maxLevelOpen: maxLevOp, overrideSteps: ovrCount, surge: Math.max(0, maxLevel - o.levelSp) * o.sepV, surgeClosed: Math.max(0, maxLevCl - o.levelSp) * o.sepV, carryOver: carry, carryOverOpen: carry - carryCl, iae: iae / HOUR, mpcActive: mp ? mp.active / Math.max(mp.moves, 1) : 0 };
}
/** Growth rate (largest real part of the eigenvalues, 1/s) of a small matrix from repeated squaring of its exponential — robust for widely spread poles. */
export function growthRate(M, T0 = 300, m = 14) {
  let E = expm(M.map((r) => r.map((v) => v * T0))), ls = 0;
  for (let k = 0; k < m; k++) { E = mm(E, E); const nrm = Math.max(...E.map((r) => Math.max(...r.map(Math.abs)))) || 1e-300; E = E.map((r) => r.map((v) => v / nrm)); ls = 2 * ls + Math.log(nrm); }
  return ls / (T0 * 2 ** m);
}
/**
 * P / PI settings for an inlet-pressure loop on the linearised slugging model by a search over the closed-loop growth rate, with the
 * dead time (plus one sample) represented by a first-order Padé lag: the smallest gain (opening per bar) that reaches 60 % of the best
 * attainable decay rate, and its best integral time. Returns { kc, ti, decay (1/s, positive = stable), stable(kc, ti), poles2(kc, ti) }.
 */
export function tuneByPoles(lin, { kcLo = 1e-3, kcHi = 3, n = 15, tis = [0, 600, 1800, 5400, 14400], theta = 30 } = {}) {
  const sc = lin.ys.map((v) => Math.abs(v) || 1), A = lin.A.map((r, i) => r.map((v, j) => (v * sc[j]) / sc[i])), b = lin.B.map((r, i) => r[0] / sc[i]), cB = lin.C[0].map((v, j) => (v * sc[j]) / 1e5), sg = Math.sign(dot(cB, solveLinear(A, b.map((v) => -v)))) || -1, th = Math.max(theta, 1);
  const cl = (kc, ti) => { // states [x (4), w (Padé lag), q = ∫(−y)]: u = −kc y + (kc/ti) q, delayed input = 2w − u
    const N = 4, n2 = ti > 0 ? N + 2 : N + 1, M = zeros(n2, n2);
    for (let i = 0; i < N; i++) { for (let j = 0; j < N; j++) M[i][j] = A[i][j] + b[i] * kc * cB[j]; M[i][N] = 2 * b[i]; if (ti > 0) M[i][N + 1] = (-b[i] * kc) / ti; }
    for (let j = 0; j < N; j++) M[N][j] = (-2 * kc * cB[j]) / th; M[N][N] = -2 / th;
    if (ti > 0) { M[N][N + 1] = (2 * kc) / (ti * th); for (let j = 0; j < N; j++) M[N + 1][j] = -cB[j]; }
    return M;
  };
  const cand = [];
  for (let i = 0; i < n; i++) { const kc = sg * kcLo * (kcHi / kcLo) ** (i / (n - 1)); for (const ti of tis) cand.push({ kc, ti, decay: -growthRate(cl(kc, ti)) }); }
  const dBest = Math.max(...cand.map((c) => c.decay)), ok = cand.filter((c) => c.decay >= 0.6 * dBest && c.decay > 0), kMin = ok.length ? Math.min(...ok.map((c) => Math.abs(c.kc))) : null;
  const best = kMin === null ? cand.reduce((a, c) => (c.decay > a.decay ? c : a)) : ok.filter((c) => Math.abs(c.kc) <= kMin * 1.0001).reduce((a, c) => (c.decay > a.decay ? c : a));
  return { ...best, stable: (kc, ti) => growthRate(cl(kc, ti)) < 0, poles2: (kc, ti) => eig(cl(kc, ti)) };
}

/**
 * Nonlinear model-predictive control of the riser model by single shooting: at every sample the next nc choke openings (the last one held)
 * minimise Σ (P_in − set-point)² [bar²] + r Σ Δz² over np samples predicted with the nonlinear model (one L-stable Rosenbrock step per sample),
 * subject to zMin ≤ z ≤ zMax and a rate limit on the first move; the programme is solved by sqp(), warm-started from the shifted previous
 * solution. Full state feedback (no estimator). c: { z0, y0, sp (Pa), dSp (bar), tStep (s), np, nc, ts (s), tEnd, r, zMin, zMax, duMax, maxIter, sub }.
 * Returns { t, pIn (bar), z, sp, iae (bar·h), amp (bar peak-to-peak over the last third), evals, iterations (mean SQP iterations per move) }.
 */
export function nmpc(sm, c) {
  const o = { np: 6, nc: 2, ts: 120, tEnd: 2 * HOUR, r: 400, zMin: 0.02, zMax: 1, duMax: 0.2, maxIter: 3, dSp: 0, tStep: 0, sub: 2, ...c }, eq = sm.steady(o.z0), sp0 = (o.sp ?? eq.Pp) / 1e5;
  const stepF = (y, z, h) => { let st = ros2Step((t, s) => sm.alg(s, z).d, 0, y, h); if (!st.y.every(Number.isFinite)) { st = { y }; for (let j = 0; j < 8; j++) st = ros2Step((t, s) => sm.alg(s, z).d, 0, st.y, h / 8); } return st.y.map((x) => Math.max(x, 1e-9)); };
  let y = o.y0 ? o.y0.slice() : eq.y.map((x, i) => x * (i === 1 ? 1.0005 : 1)), zPrev = o.z0, seq = new Array(o.nc).fill(o.z0), evals = 0, its = 0, iae = 0;
  const n = Math.round(o.tEnd / o.ts), out = { t: [], pIn: [], z: [], sp: [] };
  for (let k = 0; k <= n; k++) {
    const t = k * o.ts, sp = sp0 + (t >= o.tStep ? o.dSp : 0), y0 = y, zp = zPrev;
    const cost = (u) => { let yy = y0, J = 0; evals++; for (let i = 0; i < o.np; i++) { const z = u[Math.min(i, o.nc - 1)]; yy = stepF(yy, z, o.ts); const e = sm.alg(yy, z).Pp / 1e5 - sp; J += e * e; } for (let j = 0; j < o.nc; j++) J += o.r * (u[j] - (j ? u[j - 1] : zp)) ** 2; return Number.isFinite(J) ? J : 1e12; };
    const lo = seq.map((_, j) => (j ? o.zMin : Math.max(o.zMin, zp - o.duMax))), hi = seq.map((_, j) => (j ? o.zMax : Math.min(o.zMax, zp + o.duMax))), res = sqp(cost, seq.map((x, j) => clamp(x, lo[j], hi[j])), { lo, hi, maxIter: o.maxIter, tol: 1e-6, h: 1e-4, forward: true });
    its += res.iterations; const z = clamp(res.x[0], lo[0], hi[0]); seq = [...res.x.slice(1), res.x[o.nc - 1]];
    const p = sm.alg(y, z).Pp / 1e5; out.t.push(t); out.pIn.push(p); out.z.push(z); out.sp.push(sp); iae += Math.abs(sp - p) * o.ts;
    if (k < n) for (let j = 0; j < o.sub; j++) y = stepF(y, z, o.ts / o.sub);
    zPrev = z;
  }
  const i0 = Math.floor((2 * out.pIn.length) / 3), tail = out.pIn.slice(i0);
  return { ...out, iae: iae / HOUR, amp: Math.max(...tail) - Math.min(...tail), evals, iterations: its / (n + 1), yEnd: y };
}

// ---- optimisation ---------------------------------------------------------------------------------------------------
/**
 * Linear programme max cᵀx subject to A x ≤ b, x ≥ 0 by the two-phase tableau simplex method with Bland's rule.
 * Returns { status: 'optimal' | 'infeasible' | 'unbounded', x, obj, iterations }.
 */
export function simplex(c, A, b) {
  const m = A.length, n = c.length, rows = [], basis = [], art = [];
  // columns: x (n) | slack/surplus (m) | artificial (one per negative-rhs row) | rhs
  const neg = b.map((v) => v < 0), nArt = neg.filter(Boolean).length, N = n + m + nArt;
  let ka = 0;
  for (let i = 0; i < m; i++) {
    const r = new Array(N + 1).fill(0), sg = neg[i] ? -1 : 1;
    for (let j = 0; j < n; j++) r[j] = sg * A[i][j];
    r[n + i] = sg; r[N] = sg * b[i];
    if (neg[i]) { r[n + m + ka] = 1; basis.push(n + m + ka); art.push(n + m + ka); ka++; } else basis.push(n + i);
    rows.push(r);
  }
  let it = 0;
  const pivot = (pr, pc) => { const pv = rows[pr][pc]; for (let j = 0; j <= N; j++) rows[pr][j] /= pv; for (let i = 0; i < m; i++) if (i !== pr) { const f = rows[i][pc]; if (f !== 0) for (let j = 0; j <= N; j++) rows[i][j] -= f * rows[pr][j]; } basis[pr] = pc; };
  const solve = (cost, allowed) => { // maximise costᵀ(all variables)
    for (; it < 2000; it++) {
      const red = (j) => { let z = 0; for (let i = 0; i < m; i++) z += cost[basis[i]] * rows[i][j]; return cost[j] - z; };
      let pc = -1; for (let j = 0; j < N; j++) if (allowed(j) && !basis.includes(j) && red(j) > 1e-9) { pc = j; break; }
      if (pc < 0) return 'optimal';
      let pr = -1, best = Infinity;
      for (let i = 0; i < m; i++) if (rows[i][pc] > 1e-10) { const ratio = rows[i][N] / rows[i][pc]; if (ratio < best - 1e-12 || (Math.abs(ratio - best) <= 1e-12 && basis[i] < basis[pr])) { best = ratio; pr = i; } }
      if (pr < 0) return 'unbounded';
      pivot(pr, pc);
    }
    return 'optimal';
  };
  if (nArt) {
    const c1 = new Array(N).fill(0); for (const j of art) c1[j] = -1;
    solve(c1, () => true);
    let inf = 0; for (let i = 0; i < m; i++) if (art.includes(basis[i])) inf += rows[i][N];
    if (inf > 1e-7) return { status: 'infeasible', x: new Array(n).fill(0), obj: null, iterations: it };
    for (let i = 0; i < m; i++) if (art.includes(basis[i])) { const pc = rows[i].slice(0, n + m).findIndex((v) => Math.abs(v) > 1e-9); if (pc >= 0) pivot(i, pc); } // drive degenerate artificials out
  }
  const c2 = new Array(N).fill(0); for (let j = 0; j < n; j++) c2[j] = c[j];
  const status = solve(c2, (j) => j < n + m), x = new Array(n).fill(0);
  for (let i = 0; i < m; i++) if (basis[i] < n) x[basis[i]] = rows[i][N];
  return { status, x, obj: status === 'optimal' ? dot(c, x) : null, iterations: it };
}
/** Mixed-integer LP (max cᵀx, A x ≤ b, x ≥ 0, x[i] integer for i in ints) by depth-first branch and bound on the simplex relaxation. Returns { status, x, obj, nodes }. */
export function branchBound(c, A, b, ints, { maxNodes = 400 } = {}) {
  let best = null, nodes = 0; const n = c.length, stack = [{ A, b }];
  while (stack.length && nodes < maxNodes) {
    const nd = stack.pop(), r = simplex(c, nd.A, nd.b); nodes++;
    if (r.status !== 'optimal' || (best && r.obj <= best.obj + 1e-9)) continue;
    const k = ints.find((i) => Math.abs(r.x[i] - Math.round(r.x[i])) > 1e-6);
    if (k === undefined) { best = { x: r.x.map((v, i) => (ints.includes(i) ? Math.round(v) : v)), obj: r.obj }; continue; }
    const e = new Array(n).fill(0); e[k] = 1;
    stack.push({ A: [...nd.A, e.map((v) => -v)], b: [...nd.b, -Math.ceil(r.x[k])] }, { A: [...nd.A, e], b: [...nd.b, Math.floor(r.x[k])] });
  }
  return best ? { status: 'optimal', ...best, nodes } : { status: 'infeasible', x: new Array(n).fill(0), obj: null, nodes };
}
/**
 * Primal log-barrier interior-point method for min f(x) subject to g_i(x) ≥ 0 (box bounds are added as constraints): damped Newton
 * steps on f − μ Σ ln g_i with finite-difference derivatives and a feasibility-preserving line search; μ is reduced geometrically.
 * x0 must be strictly feasible. Returns { x, f, iterations, history (objective per outer iteration), mu }.
 */
export function interiorPoint(f, cons, x0, { lo, hi, mu = 1, shrink = 0.2, outer = 9, inner = 25, h = 1e-4 } = {}) {
  const n = x0.length, gs = [...cons]; if (lo) lo.forEach((l, i) => { gs.push((x) => x[i] - l); gs.push((x) => hi[i] - x[i]); });
  const scale = lo ? lo.map((l, i) => hi[i] - l) : x0.map((v) => Math.abs(v) + 1), phi = (x, m) => { let s = f(x); for (const g of gs) { const v = g(x); if (!(v > 0)) return Infinity; s -= m * Math.log(v); } return s; };
  let x = x0.slice(), it = 0; const hist = [];
  if (!Number.isFinite(phi(x, mu))) throw new Error('The interior-point method needs a strictly feasible starting point.');
  for (let o = 0; o < outer; o++) {
    for (let k = 0; k < inner; k++, it++) {
      const p0 = phi(x, mu), gr = new Array(n), Hm = zeros(n, n), e = scale.map((s) => s * h), at = (d) => phi(x.map((v, i) => v + d[i]), mu), z = new Array(n).fill(0);
      const fp = [], fm = [];
      for (let i = 0; i < n; i++) { const d = z.slice(); d[i] = e[i]; fp.push(at(d)); d[i] = -e[i]; fm.push(at(d)); gr[i] = (fp[i] - fm[i]) / (2 * e[i]); Hm[i][i] = (fp[i] - 2 * p0 + fm[i]) / (e[i] * e[i]); }
      for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) { const d = z.slice(); d[i] = e[i]; d[j] = e[j]; const pp = at(d); d[j] = -e[j]; const pm = at(d); d[i] = -e[i]; const mmv = at(d); d[j] = e[j]; const mp = at(d); Hm[i][j] = Hm[j][i] = (pp - pm - mp + mmv) / (4 * e[i] * e[j]); }
      if (![...gr, ...Hm.flat()].every(Number.isFinite)) break;
      let step, reg = 0;
      for (let tries = 0; tries < 8; tries++) { try { step = solveLinear(Hm.map((r, i) => r.map((v, j) => v + (i === j ? reg : 0))), gr.map((v) => -v)); if (dot(step, gr) < 0) break; } catch { step = null; } reg = reg ? reg * 10 : 1e-6 * Math.max(...Hm.map((r, i) => Math.abs(r[i])), 1e-12); step = null; }
      if (!step) step = gr.map((v, i) => -v * scale[i] * scale[i]);
      let a = 1, moved = false;
      for (let ls = 0; ls < 30; ls++, a *= 0.5) { const xn = x.map((v, i) => v + a * step[i]), pn = phi(xn, mu); if (pn < p0 - 1e-4 * a * Math.abs(dot(step, gr))) { x = xn; moved = true; break; } }
      if (!moved || Math.max(...step.map((v, i) => Math.abs(a * v) / scale[i])) < 1e-7) break;
    }
    hist.push(f(x)); mu *= shrink;
  }
  return { x, f: f(x), iterations: it, history: hist, mu };
}
/** Real-coded genetic algorithm (tournament selection, blend crossover, Gaussian mutation, elitism); ints lists integer variables. Returns { x, f, history, evals }. */
export function geneticAlgorithm(f, lo, hi, { pop = 24, gens = 30, seed = 3, pm = 0.2, ints = [] } = {}) {
  const r = rng(seed), n = lo.length, fix = (x) => x.map((v, i) => { const c = clamp(v, lo[i], hi[i]); return ints.includes(i) ? Math.round(c) : c; });
  let P = Array.from({ length: pop }, () => fix(lo.map((l, i) => r.uniform(l, hi[i])))), F = P.map(f), evals = pop; const hist = [];
  const pick = () => { const a = r.int(pop), b = r.int(pop); return F[a] < F[b] ? P[a] : P[b]; };
  for (let g = 0; g < gens; g++) {
    const bi = F.indexOf(Math.min(...F)), Q = [P[bi].slice()];
    while (Q.length < pop) { const a = pick(), b = pick(); Q.push(fix(a.map((v, i) => { const u = r.uniform(-0.25, 1.25); let c = v + u * (b[i] - v); if (r.uniform() < pm) c += r.normal(0, 0.1 * (hi[i] - lo[i]) * (1 - g / gens) + 1e-12); return c; }))); }
    const FQ = Q.map((x, i) => (i === 0 ? F[bi] : f(x))); evals += pop - 1; P = Q; F = FQ; hist.push(Math.min(...F));
  }
  const bi = F.indexOf(Math.min(...F));
  return { x: P[bi], f: F[bi], history: hist, evals };
}
/** Particle-swarm optimisation with inertia damping and reflecting bounds. Returns { x, f, history, evals }. */
export function particleSwarm(f, lo, hi, { n = 20, iters = 40, seed = 5, w = 0.72, c1 = 1.5, c2 = 1.5 } = {}) {
  const r = rng(seed), d = lo.length, X = Array.from({ length: n }, () => lo.map((l, i) => r.uniform(l, hi[i]))), V = X.map((x) => x.map((_, i) => r.uniform(-0.1, 0.1) * (hi[i] - lo[i]))), pB = X.map((x) => x.slice()), pF = X.map(f);
  let g = pF.indexOf(Math.min(...pF)), gB = pB[g].slice(), gF = pF[g]; const hist = [];
  for (let it = 0; it < iters; it++) {
    for (let k = 0; k < n; k++) {
      for (let i = 0; i < d; i++) { V[k][i] = w * V[k][i] + c1 * r.uniform() * (pB[k][i] - X[k][i]) + c2 * r.uniform() * (gB[i] - X[k][i]); X[k][i] += V[k][i]; if (X[k][i] < lo[i]) { X[k][i] = lo[i]; V[k][i] *= -0.5; } if (X[k][i] > hi[i]) { X[k][i] = hi[i]; V[k][i] *= -0.5; } }
      const v = f(X[k]); if (v < pF[k]) { pF[k] = v; pB[k] = X[k].slice(); if (v < gF) { gF = v; gB = X[k].slice(); } }
    }
    hist.push(gF);
  }
  return { x: gB, f: gF, history: hist, evals: n * (iters + 1) };
}
/** Gaussian-process regression (squared-exponential kernel on inputs scaled to the unit box). Returns predict(x) -> { mean, sd }. */
export function gaussianProcess(X, y, { len = 0.25, noise = 1e-6 } = {}) {
  const n = X.length, my = mean(y), sy = Math.sqrt(mean(y.map((v) => (v - my) ** 2))) || 1, yn = y.map((v) => (v - my) / sy), kf = (a, b) => { let s = 0; for (let i = 0; i < a.length; i++) s += (a[i] - b[i]) ** 2; return Math.exp(-s / (2 * len * len)); };
  const K = X.map((a, i) => X.map((b, j) => kf(a, b) + (i === j ? noise : 0))), Lc = zeros(n, n);
  for (let i = 0; i < n; i++) for (let j = 0; j <= i; j++) { let s = K[i][j]; for (let k = 0; k < j; k++) s -= Lc[i][k] * Lc[j][k]; Lc[i][j] = i === j ? Math.sqrt(Math.max(s, 1e-12)) : s / Lc[j][j]; }
  const fwd = (bv) => { const z = new Array(n); for (let i = 0; i < n; i++) { let s = bv[i]; for (let k = 0; k < i; k++) s -= Lc[i][k] * z[k]; z[i] = s / Lc[i][i]; } return z; }, bwd = (z) => { const x = new Array(n); for (let i = n - 1; i >= 0; i--) { let s = z[i]; for (let k = i + 1; k < n; k++) s -= Lc[k][i] * x[k]; x[i] = s / Lc[i][i]; } return x; };
  const alpha = bwd(fwd(yn));
  return (x) => { const ks = X.map((a) => kf(a, x)), v = fwd(ks); return { mean: my + sy * dot(ks, alpha), sd: sy * Math.sqrt(Math.max(1 + noise - dot(v, v), 1e-12)) }; };
}
/** Bayesian optimisation (minimisation): Gaussian-process surrogate + expected improvement maximised over random candidates. Returns { x, f, history, evals }. */
export function bayesOpt(f, lo, hi, { n0 = 6, iters = 14, seed = 9, cand = 300, len = 0.25 } = {}) {
  const r = rng(seed), d = lo.length, toX = (u) => u.map((v, i) => lo[i] + v * (hi[i] - lo[i])), U = Array.from({ length: n0 }, (_, k) => lo.map(() => (k + r.uniform()) / n0)).map((u, k, all) => u.map((v, i) => all[(k * (i + 1)) % n0][i])), Y = U.map((u) => f(toX(u))), hist = [];
  for (let it = 0; it < iters; it++) {
    const gp = gaussianProcess(U, Y, { len }), best = Math.min(...Y); let bu = null, bEI = -1;
    for (let k = 0; k < cand; k++) { const u = lo.map(() => r.uniform()), p = gp(u), z = (best - p.mean) / p.sd, ei = (best - p.mean) * normCdf(z) + p.sd * normPdf(z); if (ei > bEI) { bEI = ei; bu = u; } }
    U.push(bu); Y.push(f(toX(bu))); hist.push(Math.min(...Y));
  }
  const bi = Y.indexOf(Math.min(...Y));
  return { x: toX(U[bi]), f: Y[bi], history: hist, evals: Y.length };
}
/** Least-squares polynomial response surface y(q) of the given degree. Returns { coef, predict(q), rmse, maxErr, loo (leave-one-out RMSE) }. */
export function responseSurface(q, y, degree = 3) {
  const deg = Math.min(degree, q.length - 1), row = (x) => Array.from({ length: deg + 1 }, (_, k) => x ** k), fit = (qq, yy) => lstsq(qq.map(row), yy), coef = fit(q, y), predict = (x) => dot(coef, row(x));
  const res = q.map((x, i) => predict(x) - y[i]); let loo = 0;
  if (q.length > deg + 2) for (let i = 0; i < q.length; i++) { const c = fit(q.filter((_, j) => j !== i), y.filter((_, j) => j !== i)); loo += (dot(c, row(q[i])) - y[i]) ** 2; }
  return { coef, predict, rmse: Math.sqrt(mean(res.map((v) => v * v))), maxErr: Math.max(...res.map(Math.abs)), loo: Math.sqrt(loo / q.length) };
}
/**
 * Operating window in rate from a scan. scan: { q[] (rate fractions, ascending) and one array per constraint }, cons: [{ key, name, type: 'min' | 'max', limit, unit }]
 * (type 'min' means the scanned value must stay ≥ limit). The window is the feasible interval around qRef (or the widest one).
 * Returns { qMin, qMax, feasible, limits: [{ name, bound: 'low' | 'high' | 'none', q, value at qRef, margin }], text[] }.
 */
export function operatingEnvelope(scan, cons, qRef = 1) {
  const n = 241, qs = linspace(scan.q[0], scan.q[scan.q.length - 1], n), ok = cons.map((c) => qs.map((x) => { const v = interp1(scan.q, scan[c.key], x); return c.type === 'min' ? v - c.limit : c.limit - v; }));
  const feas = qs.map((_, i) => ok.every((g) => g[i] >= 0));
  let i0 = -1, j0 = -1, bestLen = 0; // feasible run containing qRef, else the longest
  for (let i = 0; i < n;) { if (!feas[i]) { i++; continue; } let j = i; while (j + 1 < n && feas[j + 1]) j++; if (qs[i] <= qRef && qs[j] >= qRef) { i0 = i; bestLen = Infinity; j0 = j; } else if (j - i + 1 > bestLen) { bestLen = j - i + 1; i0 = i; j0 = j; } i = j + 1; }
  const cross = (g, i) => (i < 0 || i + 1 >= n ? null : qs[i] + ((qs[i + 1] - qs[i]) * g[i]) / (g[i] - g[i + 1] || 1e-300));
  const limits = cons.map((c, k) => { const g = ok[k], vRef = interp1(scan.q, scan[c.key], qRef); let bound = 'none', q = null; if (i0 >= 0) { if (i0 > 0 && g[i0 - 1] < 0) { bound = 'low'; q = cross(g, i0 - 1); } else if (j0 < n - 1 && g[j0 + 1] < 0) { bound = 'high'; q = cross(g, j0); } } else if (g.every((v) => v < 0)) bound = 'violated'; return { name: c.name, key: c.key, type: c.type, limit: c.limit, unit: c.unit || '', bound, q, value: vRef, margin: c.type === 'min' ? vRef - c.limit : c.limit - vRef }; });
  if (i0 < 0) { const never = limits.filter((l) => l.bound === 'violated'); return { qMin: null, qMax: null, feasible: false, limits, text: never.length ? never.map((l) => `${l.name} cannot be met at any scanned rate`) : ['the individual limits leave no common rate range'] }; }
  const lows = limits.filter((l) => l.bound === 'low'), highs = limits.filter((l) => l.bound === 'high'), qMin = lows.length ? Math.max(...lows.map((l) => l.q)) : qs[i0], qMax = highs.length ? Math.min(...highs.map((l) => l.q)) : qs[j0];
  const text = [...lows.map((l) => `minimum rate ${(100 * l.q).toFixed(0)} % set by ${l.name}`), ...highs.map((l) => `maximum rate ${(100 * l.q).toFixed(0)} % set by ${l.name}`)];
  if (!lows.length) text.push(`no lower limit down to ${(100 * qs[i0]).toFixed(0)} % (scan floor)`); if (!highs.length) text.push(`no upper limit up to ${(100 * qs[j0]).toFixed(0)} % (scan ceiling)`);
  return { qMin, qMax, feasible: true, limits, text };
}

// ---- operating logic: event scheduler, alarms and interlocks, state machine ------------------------------------------
/**
 * Discrete-event scheduler over a time march: the step is cut at every event so that events fire at their exact times.
 * events: [{ t (s), tag, action, set: { name: value } }]; opt: { tEnd, dt, state (initial variables), onStep(t0, t1, state) }.
 * Returns { fired: [{ t, tag, action }], steps, state }.
 */
export function eventScheduler(events, { tEnd, dt, state = {}, onStep } = {}) {
  const q = events.map((e, i) => ({ ...e, i })).filter((e) => e.t >= 0 && e.t <= tEnd).sort((a, b) => a.t - b.t || a.i - b.i), fired = [], st = { ...state };
  let t = 0, k = 0, steps = 0;
  const fire = () => { while (k < q.length && q[k].t <= t) { Object.assign(st, q[k].set || {}); fired.push({ t: q[k].t, tag: q[k].tag, action: q[k].action || '' }); k++; } };
  fire();
  while (t < tEnd) { const t1 = Math.min(tEnd, t + dt, k < q.length ? q[k].t : Infinity); if (onStep) onStep(t, t1, st); t = t1; steps++; fire(); }
  return { fired, steps, state: st };
}
/** k-out-of-n voting on boolean trip signals. */
export const vote = (signals, k) => signals.filter(Boolean).length >= k;
/**
 * Alarm and interlock evaluation. values: { key: number }, rules: [{ tag, key, type: 'low' | 'high', limit, level: 'alarm' | 'trip', msg, action }].
 * Returns the raised entries [{ tag, level, msg, value, limit, action }]; values that are missing or not finite raise nothing.
 */
export function evaluateAlarms(values, rules) {
  const out = [];
  for (const r of rules) { const v = values[r.key]; if (!isNum(v) || !isNum(+r.limit)) continue; if (r.type === 'low' ? v < +r.limit : v > +r.limit) out.push({ tag: String(r.tag), level: r.level === 'trip' ? 'trip' : 'alarm', msg: `${r.msg || r.key} ${r.type === 'low' ? 'below' : 'above'} ${+r.limit} (now ${+v.toFixed(2)})`, value: v, limit: +r.limit, action: r.action || '' }); }
  return out;
}
/** Operating states and the permitted transitions of the shutdown / preservation / restart sequence. FAILSAFE holds the fail-safe valve positions. */
export const OPS_STATES = Object.freeze({
  PRODUCING: { shutdown: 'SHUT_IN', trip: 'FAILSAFE', pig: 'PIGGING' },
  PIGGING: { pigReceived: 'PRODUCING', trip: 'FAILSAFE', shutdown: 'SHUT_IN' },
  SHUT_IN: { inhibit: 'PRESERVED', blowdown: 'DEPRESSURISED', heat: 'PRESERVED', restart: 'RESTARTING', trip: 'FAILSAFE' },
  PRESERVED: { blowdown: 'DEPRESSURISED', restart: 'RESTARTING', trip: 'FAILSAFE' },
  DEPRESSURISED: { repressurise: 'RESTARTING', restart: 'RESTARTING', trip: 'FAILSAFE' },
  RESTARTING: { rampDone: 'PRODUCING', shutdown: 'SHUT_IN', trip: 'FAILSAFE' },
  FAILSAFE: { reset: 'SHUT_IN' },
});
export const FAILSAFE_POSITIONS = Object.freeze({ 'Production choke': 'closed', 'Wing / master valves': 'closed', 'Riser ESD valve': 'closed', 'Blowdown valve': 'open on demand (fail-open)', 'Methanol injection valve': 'open (fail-open)', 'Heating': 'off', 'Pig launcher': 'isolated' });
/** Run a list of events through the state machine. Returns { state, log: [{ event, from, to, accepted }], rejected }. */
export function stateMachine(events, start = 'PRODUCING') {
  let s = start; const log = [];
  for (const e of events) { const to = OPS_STATES[s]?.[e]; log.push({ event: e, from: s, to: to || s, accepted: !!to }); if (to) s = to; }
  return { state: s, log, rejected: log.filter((l) => !l.accepted).length };
}
/**
 * Liquid surge at the receiving facility while the rate is ramped: the line inventory relaxes towards the steady inventory of the
 * current rate over the liquid residence time and the difference leaves (or is held back) at the outlet.
 * o: { nodes: [[t (s), q (fraction)], …] ramp schedule, inv(q) (m³), qLiq(q) (m³/s steady liquid outflow), qDrain (m³/s), inv0 (m³), dt, tEnd }
 * Returns { t[], q[], qOut[] (m³/s), V[] (m³ accumulated above the drain capacity), vMax, deferred (s of full-rate production lost), inventory[] }.
 */
export function rampSurge(o) {
  const tn = o.nodes.map((p) => p[0]), qn = o.nodes.map((p) => p[1]), q1 = qn[qn.length - 1], dt = o.dt || 60, tEnd = o.tEnd ?? tn[tn.length - 1] * 1.5 + 4 * HOUR, n = Math.ceil(tEnd / dt);
  let inv = o.inv0 ?? o.inv(qn[0]), V = 0, vMax = 0, def = 0; const out = { t: [], q: [], qOut: [], V: [], inventory: [] };
  for (let k = 0; k <= n; k++) {
    const t = k * dt, q = interp1(tn, qn, t), ql = o.qLiq(q), tau = clamp(o.inv(q) / Math.max(ql, 1e-6), 600, 12 * HOUR), dInv = (o.inv(q) - inv) / tau, qOut = Math.max(0, ql - dInv);
    out.t.push(t); out.q.push(q); out.qOut.push(qOut); out.V.push(V); out.inventory.push(inv);
    inv += (ql - qOut) * dt; V = Math.max(0, V + (qOut - o.qDrain) * dt); vMax = Math.max(vMax, V); def += Math.max(0, 1 - q / q1) * dt;
  }
  return { ...out, vMax, deferred: def };
}

// ---- frequency response, sensitivity peaks and robust PI tuning -----------------------------------------------
const cMul = (a, b) => [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]];
const cDiv = (a, b) => { const d = b[0] * b[0] + b[1] * b[1] || 1e-300; return [(a[0] * b[0] + a[1] * b[1]) / d, (a[1] * b[0] - a[0] * b[1]) / d]; };
const cPoly = (c, w) => { let p = [0, 0]; for (const a of c) p = [-p[1] * w + a, p[0] * w]; return p; }; // Horner at s = jw, descending powers
/**
 * Frequency response G(jw) of a plant { k, lags[] (τ of 1/(τs+1)), zeros[] (T of (Ts+1); negative = right-half-plane zero), ints (number of
 * integrators), delay (s), num[], den[] (optional polynomials in s, descending powers) }. Returns [re, im].
 */
export function plantResponse(p, w) {
  let g = [p.k ?? 1, 0];
  if (p.num && p.den) g = cMul(g, cDiv(cPoly(p.num, w), cPoly(p.den, w)));
  for (const T of p.zeros || []) g = cMul(g, [1, w * T]);
  for (const T of p.lags || []) g = cDiv(g, [1, w * T]);
  for (let i = 0; i < (p.ints || 0); i++) g = cDiv(g, [0, w]);
  if (p.delay) g = cMul(g, [Math.cos(w * p.delay), -Math.sin(w * p.delay)]);
  return g;
}
/** Frequency response of a controller { kc, ti, td, N, form: 'ideal' | 'series' } or a pure integrator { ki }. */
export function controllerResponse(c, w) {
  if (c.ki !== undefined && !(c.kc)) return [0, -c.ki / w];
  const I = c.ti > 0 ? [1, -1 / (w * c.ti)] : [1, 0], td = c.td || 0;
  if (c.form === 'series') return cMul([c.kc, 0], cMul(I, [1, w * td]));
  const N = c.N ?? 10, D = td > 0 ? cDiv([0, w * td], [1, (w * td) / N]) : [0, 0];
  return [c.kc * (I[0] + D[0]), c.kc * (I[1] + D[1])];
}
/**
 * Loop-shaping measures of L = C G from a frequency scan (logarithmic grid, peaks refined by golden section):
 * sensitivity peak Ms = max |1/(1+L)|, complementary sensitivity peak Mt = max |L/(1+L)|, gain margin, phase margin (deg) and
 * the crossover frequencies. The numbers describe the closed loop only when it is stable (check with the poles or a simulation).
 */
export function loopAnalysis(plant, ctrl, { wLo, wHi, n = 1200, fine = true } = {}) {
  const tRef = Math.max(plant.delay || 0, ...(plant.lags || []).map(Math.abs), ...(plant.zeros || []).map(Math.abs), plant.tRef || 0) || 1;
  const a = Math.log(wLo ?? 1e-4 / tRef), b = Math.log(wHi ?? 2e3 / tRef), L = (w) => cMul(controllerResponse(ctrl, w), plantResponse(plant, w));
  const sOf = (w) => { const l = L(w); return 1 / Math.hypot(1 + l[0], l[1]); }, tOf = (w) => { const l = L(w); return Math.hypot(l[0], l[1]) / Math.hypot(1 + l[0], l[1]); };
  const ws = [], S = [], T = [], mag = [], ph = [];
  let prev = null, acc = 0;
  for (let i = 0; i < n; i++) { const w = Math.exp(a + ((b - a) * i) / (n - 1)), l = L(w), m = Math.hypot(l[0], l[1]), d = Math.hypot(1 + l[0], l[1]), p = Math.atan2(l[1], l[0]); if (prev !== null) { let dp = p - prev; while (dp > Math.PI) dp -= 2 * Math.PI; while (dp < -Math.PI) dp += 2 * Math.PI; acc += dp; } else acc = p > 0.5 ? p - 2 * Math.PI : p; prev = p; ws.push(w); S.push(1 / d); T.push(m / d); mag.push(m); ph.push(acc); }
  const refine = (f, arr) => { let k = 0; for (let i = 1; i < n; i++) if (arr[i] > arr[k]) k = i; if (!fine) return [arr[k], ws[k]]; let lo = Math.log(ws[Math.max(k - 1, 0)]), hi = Math.log(ws[Math.min(k + 1, n - 1)]); const g = 0.381966; for (let it = 0; it < 40; it++) { const x1 = lo + g * (hi - lo), x2 = hi - g * (hi - lo); if (f(Math.exp(x1)) > f(Math.exp(x2))) hi = x2; else lo = x1; } const w = Math.exp(0.5 * (lo + hi)); return [Math.max(f(w), arr[k]), w]; };
  const [ms, wMs] = refine(sOf, S), [mt] = refine(tOf, T);
  let gm = null, pm = null, wc = null, w180 = null;
  for (let i = 1; i < n; i++) {
    if (pm === null && mag[i - 1] >= 1 && mag[i] < 1) { const f = (mag[i - 1] - 1) / (mag[i - 1] - mag[i]); wc = ws[i - 1] * (ws[i] / ws[i - 1]) ** f; pm = 180 + ((ph[i - 1] + f * (ph[i] - ph[i - 1])) * 180) / Math.PI; }
    if (gm === null && ph[i - 1] > -Math.PI && ph[i] <= -Math.PI) { const f = (ph[i - 1] + Math.PI) / (ph[i - 1] - ph[i]); w180 = ws[i - 1] * (ws[i] / ws[i - 1]) ** f; gm = 1 / (mag[i - 1] * (mag[i] / mag[i - 1]) ** f); }
  }
  return { ms, mt, wMs, gm, pm, wc, w180 };
}
/**
 * Robust PI tuning over a set of plant models (multi-model, no μ-analysis): the PI setting with the largest integral gain Kc/Ti whose
 * sensitivity peak stays at or below msMax for EVERY model of the set (the Åström–Schei criterion applied to the worst case).
 * models: plants for plantResponse; the first is the nominal one. Returns { kc, ti, ki, ms[] (per model), mt[] , worst }.
 */
export function robustPI(models, { msMax = 1.6, tis, n = 110 } = {}) {
  const m0 = models[0], sg = Math.sign(m0.k ?? 1) || 1, th = Math.max(m0.delay || 0, 1e-9), tau = Math.max(...(m0.lags || []), th), tiList = tis || Array.from({ length: 8 }, (_, i) => Math.min(tau, 2 * th) * (Math.max(tau, 16 * th) / Math.min(tau, 2 * th)) ** (i / 7));
  const wLo = 0.02 / Math.max(tau, th), wHi = 12 / th, worst = (kc, ti) => { let w = 0; for (const m of models) { w = Math.max(w, loopAnalysis(m, { kc, ti }, { n, fine: false, wLo, wHi }).ms); if (w > msMax) break; } return w; };
  let best = null; const k0 = (0.25 * tau) / (Math.abs(m0.k ?? 1) * th); // a quarter of the SIMC gain as the first trial
  for (const ti of tiList) {
    let lo = 0, hi = k0; for (let k = 0; k < 12 && worst(sg * hi, ti) < msMax; k++) { lo = hi; hi *= 2; }
    for (let k = 0; k < 9; k++) { const mid = 0.5 * (lo + hi); if (worst(sg * mid, ti) <= msMax) lo = mid; else hi = mid; }
    if (lo > 0 && (!best || lo / ti > best.ki)) best = { kc: sg * lo, ti, ki: lo / ti };
  }
  if (!best) return null;
  const an = models.map((m) => loopAnalysis(m, { kc: best.kc, ti: best.ti }, { n: 600 }));
  return { ...best, ms: an.map((x) => x.ms), mt: an.map((x) => x.mt), worst: Math.max(...an.map((x) => x.ms)), gm: an[0].gm, pm: an[0].pm };
}

// ---- sequential quadratic programming --------------------------------------------------------------------------
/** Exact active-set solution of the bound-constrained convex QP min ½ xᵀH x + fᵀx, lo ≤ x ≤ hi (H positive definite; ±Infinity bounds allowed). */
export function boxQPExact(H, f, lo, hi, x0) {
  const n = f.length, x = (x0 || new Array(n).fill(0)).map((v, i) => clamp(v, lo[i], hi[i])), act = x.map((v, i) => (v <= lo[i] ? -1 : v >= hi[i] ? 1 : 0));
  const ridge = 1e-12 * Math.max(...H.map((r, i) => Math.abs(r[i])), 1e-300);
  for (let it = 0; it < 12 * n + 40; it++) {
    const fr = []; for (let i = 0; i < n; i++) if (!act[i]) fr.push(i);
    let hit = -1;
    if (fr.length) {
      const A = fr.map((i) => fr.map((j) => H[i][j] + (i === j ? ridge : 0))), b = fr.map((i) => { let s = -f[i]; for (let j = 0; j < n; j++) if (act[j]) s -= H[i][j] * x[j]; return s; });
      let xt; try { xt = solveLinear(A, b); } catch { xt = fr.map((i) => x[i]); }
      let al = 1;
      fr.forEach((i, k) => { const d = xt[k] - x[i]; if (d < 0 && xt[k] < lo[i] - 1e-13) { const a2 = (lo[i] - x[i]) / d; if (a2 < al) { al = a2; hit = i; } } else if (d > 0 && xt[k] > hi[i] + 1e-13) { const a2 = (hi[i] - x[i]) / d; if (a2 < al) { al = a2; hit = i; } } });
      fr.forEach((i, k) => { x[i] += al * (xt[k] - x[i]); });
      if (hit >= 0) { act[hit] = x[hit] - lo[hit] < hi[hit] - x[hit] ? -1 : 1; x[hit] = act[hit] < 0 ? lo[hit] : hi[hit]; continue; }
    }
    let worst = 0, rel = -1; // multipliers of the active bounds: release the most wrong-signed one
    for (let i = 0; i < n; i++) if (act[i]) { let g = f[i]; for (let j = 0; j < n; j++) g += H[i][j] * x[j]; const v = act[i] < 0 ? -g : g; if (v > worst + 1e-11) { worst = v; rel = i; } }
    if (rel < 0) break;
    act[rel] = 0;
  }
  return x;
}
/**
 * Sequential quadratic programming for min f(x) subject to eq_i(x) = 0, ineq_j(x) ≥ 0 and lo ≤ x ≤ hi: damped BFGS approximation of
 * the Hessian of the Lagrangian (Powell), QP subproblem solved through its dual by an exact active-set method, ℓ1 merit function with a
 * backtracking line search, central-difference derivatives. Returns { x, f, iterations, violation, lambda, converged, history }.
 */
export function sqp(f, x0, { eq = [], ineq = [], lo, hi, maxIter = 80, tol = 1e-8, h = 1e-6, forward = false } = {}) {
  const n = x0.length, me = eq.length, mi = ineq.length, m = me + mi, L = lo || new Array(n).fill(-Infinity), U = hi || new Array(n).fill(Infinity);
  const cons = (x) => [...eq.map((c) => c(x)), ...ineq.map((c) => c(x))];
  const grad = forward ? (fn, x) => { const f0 = fn(x); return x.map((v, i) => { const e = h * (1 + Math.abs(v)), xp = x.slice(), s = Number.isFinite(U[i]) && v + e > U[i] ? -1 : 1; xp[i] = v + s * e; return (fn(xp) - f0) / (s * e); }); }
    : (fn, x) => x.map((v, i) => { const e = h * (1 + Math.abs(v)), xp = x.slice(), xm = x.slice(); xp[i] = v + e; xm[i] = v - e; return (fn(xp) - fn(xm)) / (2 * e); });
  const jac = (x) => [...eq, ...ineq].map((c) => grad(c, x)), viol = (c) => { let s = 0; for (let i = 0; i < m; i++) s += i < me ? Math.abs(c[i]) : Math.max(0, -c[i]); return s; };
  let x = x0.map((v, i) => clamp(v, L[i], U[i])), B = eye(n), fx = f(x), g = grad(f, x), c = cons(x), A = jac(x), mu = 1, lam = new Array(m).fill(0), it = 0, converged = false, stall = 0; const hist = [fx];
  for (; it < maxIter; it++) {
    // rows: general constraints, then the finite bounds as d_i ≥ lo_i − x_i and −d_i ≥ x_i − hi_i
    const rows = A.map((r) => r.slice()), rhs = c.slice(), free = []; for (let i = 0; i < me; i++) free.push(i);
    for (let i = 0; i < n; i++) { if (Number.isFinite(L[i])) { const r = new Array(n).fill(0); r[i] = 1; rows.push(r); rhs.push(x[i] - L[i]); } if (Number.isFinite(U[i])) { const r = new Array(n).fill(0); r[i] = -1; rows.push(r); rhs.push(U[i] - x[i]); } }
    const nr = rows.length; let d, lamQ = new Array(nr).fill(0);
    try {
      const Bi = inv(B), BiAt = rows.map((r) => mv(Bi, r)), Big = mv(Bi, g);
      if (nr) {
        const M = rows.map((r) => BiAt.map((q) => dot(r, q))), q = rows.map((r, i) => rhs[i] - dot(r, Big)), sc = Math.max(...M.map((r, i) => r[i]), 1e-300);
        for (let i = 0; i < nr; i++) M[i][i] += 1e-10 * sc;
        lamQ = boxQPExact(M, q, rows.map((_, i) => (i < me ? -1e9 : 0)), new Array(nr).fill(1e9), lamQ);
      }
      d = Big.map((v, i) => { let s = -v; for (let k = 0; k < nr; k++) s += lamQ[k] * BiAt[k][i]; return s; });
    } catch { d = null; }
    if (!d || !d.every(Number.isFinite)) { d = g.map((v) => -v); lamQ = new Array(nr).fill(0); }
    const dn = Math.max(...d.map(Math.abs)), xn = Math.max(...x.map(Math.abs), 1); if (dn > 1e3 * xn) d = d.map((v) => (v * 1e3 * xn) / dn);
    lam = lamQ.slice(0, m);
    const v0 = viol(c);
    if (Math.max(...d.map(Math.abs)) <= tol * xn && v0 <= 1e2 * tol) { converged = true; break; }
    mu = Math.max(mu, 1.5 * Math.max(0, ...lam.map(Math.abs)));
    const merit = (xx) => f(xx) + mu * viol(cons(xx)), p0 = fx + mu * v0, D = dot(g, d) - mu * v0;
    let al = 1, xNew = null;
    for (let ls = 0; ls < 30; ls++, al *= 0.5) { const xt = x.map((v, i) => clamp(v + al * d[i], L[i], U[i])), pt = merit(xt); if (Number.isFinite(pt) && pt <= p0 + 1e-4 * al * Math.min(D, 0)) { xNew = xt; break; } }
    if (!xNew) { converged = dn <= 1e3 * tol * xn && v0 <= 1e2 * tol; break; }
    const gOld = g, AOld = A, s = xNew.map((v, i) => v - x[i]), fOld = fx; x = xNew; fx = f(x); g = grad(f, x); c = cons(x); A = jac(x); hist.push(fx);
    stall = Math.abs(fx - fOld) <= 1e-13 * (1 + Math.abs(fx)) && viol(c) <= 1e2 * tol ? stall + 1 : 0; if (stall >= 3) { converged = true; break; }
    const gl = (gg, AA) => gg.map((v, i) => { let t = v; for (let k = 0; k < m; k++) t -= lam[k] * AA[k][i]; return t; }), y = gl(g, A).map((v, i) => v - gl(gOld, AOld)[i]), Bs = mv(B, s), sBs = dot(s, Bs), sy = dot(s, y);
    if (sBs > 1e-300) { const th = sy >= 0.2 * sBs ? 1 : (0.8 * sBs) / (sBs - sy), r = y.map((v, i) => th * v + (1 - th) * Bs[i]), sr = dot(s, r); if (sr > 1e-300) B = B.map((row, i) => row.map((v, j) => v - (Bs[i] * Bs[j]) / sBs + (r[i] * r[j]) / sr)); }
  }
  return { x, f: fx, iterations: it, violation: viol(c), lambda: lam, converged, history: hist };
}

// ---- unscented Kalman filter and recursive least squares -------------------------------------------------------
const cholesky = (P) => { const n = P.length, Lc = zeros(n, n); for (let i = 0; i < n; i++) for (let j = 0; j <= i; j++) { let s = P[i][j]; for (let k = 0; k < j; k++) s -= Lc[i][k] * Lc[j][k]; Lc[i][j] = i === j ? Math.sqrt(Math.max(s, 1e-300)) : s / Lc[j][j]; } return Lc; };
/**
 * Unscented Kalman filter for x⁺ = F(x, u), y = h(x, u): 2L + 1 sigma points x̂ ± α√κ · columns of chol(P) with the weights
 * W0ᵃ = (α²κ − L)/(α²κ), W0ᶜ = W0ᵃ + 1 − α² + β, Wj = 1/(2α²κ); the sigma points are redrawn from the predicted covariance for the update.
 * o: { F, h, x0, P0, Q, R, u[], y[][], alpha (1), beta (2), kappa (1.5 L) }. Returns { x: [[…]], P }.
 */
export function ukf({ F, h, x0, P0, Q, R: Rn, u, y, alpha = 1, beta = 2, kappa }) {
  const Ln = x0.length, kap = kappa ?? 1.5 * Ln, a2k = alpha * alpha * kap, w0a = (a2k - Ln) / a2k, w0c = w0a + 1 - alpha * alpha + beta, wj = 1 / (2 * a2k), sp = Math.sqrt(a2k);
  const sigma = (x, P) => { const C = cholesky(P), pts = [x.slice()]; for (let j = 0; j < Ln; j++) { pts.push(x.map((v, i) => v + sp * C[i][j])); pts.push(x.map((v, i) => v - sp * C[i][j])); } return pts; };
  const wm = (pts) => pts[0].map((_, i) => { let s = w0a * pts[0][i]; for (let j = 1; j < pts.length; j++) s += wj * pts[j][i]; return s; });
  const cov = (X, mx, Y, my) => { const C = zeros(mx.length, my.length); for (let j = 0; j < X.length; j++) { const w = j ? wj : w0c; for (let a = 0; a < mx.length; a++) { const da = X[j][a] - mx[a]; for (let b = 0; b < my.length; b++) C[a][b] += w * da * (Y[j][b] - my[b]); } } return C; };
  let x = x0.slice(), P = P0.map((r) => r.slice()); const xs = [];
  for (let k = 0; k < y.length; k++) {
    const uk = u ? u[k] : 0, Xp = sigma(x, P).map((s) => F(s, uk)), xm = wm(Xp), Pm = madd(cov(Xp, xm, Xp, xm), Q), Xs = sigma(xm, Pm.map((r, i) => r.map((v, j) => 0.5 * (v + Pm[j][i])))), Z = Xs.map((s) => h(s, uk)), zm = wm(Z);
    const S = madd(cov(Z, zm, Z, zm), Rn), Pxz = cov(Xs, xm, Z, zm), K = mm(Pxz, inv(S)), innov = y[k].map((v, i) => v - zm[i]);
    x = xm.map((v, i) => v + dot(K[i], innov)); P = madd(Pm, mm(mm(K, S), tr(K)), -1); P = P.map((r, i) => r.map((v, j) => 0.5 * (v + P[j][i])));
    xs.push(x.slice());
  }
  return { x: xs, P };
}
/** Recursive least squares with exponential forgetting: y_k = φ_kᵀθ. Returns { theta, P, history: [[θ…]] }. */
export function rls(phi, y, { lam = 1, p0 = 1e6, theta0 } = {}) {
  const n = phi[0].length; let th = theta0 ? theta0.slice() : new Array(n).fill(0), P = eye(n).map((r) => r.map((v) => v * p0)); const hist = [];
  for (let k = 0; k < y.length; k++) { const f = phi[k], Pf = mv(P, f), den = lam + dot(f, Pf), e = y[k] - dot(f, th); th = th.map((v, i) => v + (Pf[i] * e) / den); P = P.map((r, i) => r.map((v, j) => (v - (Pf[i] * Pf[j]) / den) / lam)); hist.push(th.slice()); }
  return { theta: th, P, history: hist };
}
/** First-order model K/(τs+1) from an input/output record sampled every dt by recursive least squares on y_k = a y_(k-1) + b u_(k-1-nd). */
export function rlsFirstOrder(u, y, dt, { nd = 0, lam = 1 } = {}) {
  const phi = [], yy = []; for (let k = nd + 1; k < y.length; k++) { phi.push([y[k - 1], u[k - 1 - nd]]); yy.push(y[k]); }
  const r = rls(phi, yy, { lam }), a = clamp(r.theta[0], 1e-9, 1 - 1e-12), b = r.theta[1];
  return { a, b, K: b / (1 - a), tau: -dt / Math.log(a) };
}
/**
 * Self-tuning PI on a first-order-plus-dead-time process whose gain changes with time: recursive least squares (forgetting factor)
 * identifies y_k = a y_(k-1) + b u_(k-1-nd) on line and the SIMC PI rule (τc = max(θ, dt)) is recomputed every sample from the estimate.
 * o: { K(t) | number, tau, theta, dt, tEnd, sp(t), adapt (bool), lam, K0 (design gain of the fixed controller) }. Returns { t, y, u, kHat[], iae, kcEnd }.
 */
export function selfTuningLoop(o) {
  const q = { tau: 20, theta: 2, dt: 0.5, tEnd: 600, adapt: true, lam: 0.98, ...o }, n = Math.round(q.tEnd / q.dt), nd = Math.round(q.theta / q.dt), a = Math.exp(-q.dt / q.tau), kOf = typeof q.K === 'function' ? q.K : () => q.K, K0 = q.K0 ?? kOf(0), tc = Math.max(q.theta, q.dt);
  const simc = (K, tau) => ({ kc: tau / (K * (tc + q.theta)), ti: Math.min(tau, 4 * (tc + q.theta)) });
  let set = simc(K0, q.tau), th = [a, K0 * (1 - a)], P = eye(2).map((r) => r.map((v) => v * 100)), x = 0, I = 0, iae = 0; const us = new Array(nd + 2).fill(0), ys = [0], out = { t: [0], y: [0], u: [], kHat: [K0] };
  for (let k = 0; k < n; k++) {
    const t = k * q.dt, sp = typeof q.sp === 'function' ? q.sp(t) : q.sp ?? 1, e = sp - x; I += (set.kc * e * q.dt) / set.ti; const u = set.kc * e + I; us.push(u);
    x = a * x + kOf(t) * (1 - a) * us[us.length - 1 - nd]; ys.push(x); iae += Math.abs(sp - x) * q.dt;
    if (q.adapt) { // regressor [y_(k), u_(k-nd)] predicts y_(k+1)
      const f = [ys[ys.length - 2], us[us.length - 1 - nd]], Pf = mv(P, f), den = q.lam + dot(f, Pf), err = x - dot(f, th);
      if (Math.abs(f[0]) + Math.abs(f[1]) > 1e-9) { th = th.map((v, i) => v + (Pf[i] * err) / den); P = P.map((r, i) => r.map((v, j) => clamp((v - (Pf[i] * Pf[j]) / den) / q.lam, -1e8, 1e8))); }
      const ah = clamp(th[0], 0.2, 0.99999), Kh = th[1] / (1 - ah); if (Math.abs(Kh) > 1e-3 * Math.abs(K0) && Math.sign(Kh) === Math.sign(K0)) set = simc(Kh, -q.dt / Math.log(ah));
      out.kHat.push(th[1] / (1 - clamp(th[0], 0.2, 0.99999)));
    } else out.kHat.push(K0);
    out.t.push(t + q.dt); out.y.push(x); out.u.push(u);
  }
  return { ...out, iae, kcEnd: set.kc, tiEnd: set.ti };
}
/**
 * Ratio control of an injected stream to a wild stream: set-point = ratio × measured wild flow (feed-forward) trimmed by a PI flow loop
 * on a pump with a first-order response and a capacity limit. o: { wild(t), ratio, tau (s), kc, ti (s), uMax, dt, tEnd, q0 }.
 * Returns { t, wild, q (injected), sp, ratioMin, ratioEnd, injected (∫q dt), required (∫ratio·wild dt), shortfall (∫max(0, sp − q) dt), saturated (s) }.
 */
export function ratioControl(o) {
  const q = { ratio: 1, tau: 20, kc: 0.5, ti: 30, uMax: Infinity, dt: 1, tEnd: 600, ...o }, n = Math.round(q.tEnd / q.dt), a = Math.exp(-q.dt / q.tau), out = { t: [], wild: [], q: [], sp: [] };
  let x = q.q0 ?? q.ratio * q.wild(0), I = 0, inj = 0, req = 0, short = 0, sat = 0, rMin = Infinity;
  for (let k = 0; k <= n; k++) {
    const t = k * q.dt, w = q.wild(t), sp = q.ratio * w, e = sp - x, raw = sp + q.kc * e + I, u = clamp(raw, 0, q.uMax);
    if (u === raw) I += (q.kc * e * q.dt) / q.ti; else sat += q.dt;
    out.t.push(t); out.wild.push(w); out.q.push(x); out.sp.push(sp); if (w > 1e-12) rMin = Math.min(rMin, x / w);
    if (k < n) { inj += x * q.dt; req += sp * q.dt; short += Math.max(0, sp - x) * q.dt; x = a * x + (1 - a) * u; }
  }
  return { ...out, ratioMin: rMin, ratioEnd: out.q[n] / Math.max(out.wild[n], 1e-12), injected: inj, required: req, shortfall: short, saturated: sat };
}

// ---- rotating equipment: variable-speed pump and compressor with anti-surge control ----------------------------
const ownPumpHead = (q, { qr, hr, shutoff = 1.25, speed = 1 }) => { const s = Math.max(speed, 1e-6); return hr * s * s * (shutoff - (shutoff - 1) * (q / (qr * s)) ** 2); };
const pumpHeadFn = typeof NET.pumpHead === 'function' ? NET.pumpHead : ownPumpHead;
/** Own fan-law compressor map (used when the network suite does not provide one): speed lines from surge to stonewall. */
const ownCompressorMap = ({ qd, hd, speeds = [0.7, 0.8, 0.9, 1, 1.05], n = 11 }) => { const shape = (f) => 1 + 0.22 * (1 - f * f) - 2.2 * Math.max(f - 1.1, 0) ** 2; return { qd, hd, speeds, lines: speeds.map((N) => { const phi = linspace(0.62, 1.3, n); return { n: N, q: phi.map((f) => f * N * qd), h: phi.map((f) => N * N * hd * shape(f)), qSurge: 0.62 * N * qd, qChoke: 1.3 * N * qd }; }) }; };
/**
 * Compressor performance from a map: the 100 % speed line { q[] (inlet volume flow, ascending from the surge point), h[] (polytropic head) }
 * scaled by the fan laws (Q ∝ N, H ∝ N²). table: optional rows [{ q, h }] at 100 % speed; otherwise the design point (qd, hd) and the
 * network suite's map generator (or the built-in one). Returns { q[], h[], qSurge, qChoke, head(Q, N), flow(H, N) -> { q, surge, choke }, map }.
 */
export function compressorCurve({ qd, hd, table }) {
  let q, hh, map = null;
  const rows = (Array.isArray(table) ? table : []).map((r) => ({ q: +r.q, h: +r.h })).filter((r) => r.q > 0 && r.h > 0).sort((a, b) => a.q - b.q);
  if (rows.length >= 3) { q = rows.map((r) => r.q); hh = rows.map((r) => r.h); }
  else { map = (typeof NET.compressorMap === 'function' ? NET.compressorMap : ownCompressorMap)({ qd, hd }); const l = map.lines.reduce((b, x) => (Math.abs(x.n - 1) < Math.abs(b.n - 1) ? x : b)); q = l.q.map((v) => v / l.n); hh = l.h.map((v) => v / (l.n * l.n)); }
  let kPk = 0; for (let i = 1; i < q.length; i++) if (hh[i] > hh[kPk]) kPk = i; // the stable branch starts at the head maximum
  const qs = q.slice(kPk), hs = hh.slice(kPk), nn = qs.length, qSurge = qs[0], qChoke = qs[nn - 1], hr = hs.slice().reverse(), qr = qs.slice().reverse();
  const head = (Q, N) => N * N * interp1(qs, hs, clamp(Q / Math.max(N, 1e-9), qSurge, qChoke));
  const flow = (H, N) => { const x = H / Math.max(N * N, 1e-12); if (x > hs[0]) return { q: 0, surge: true, choke: false }; if (x <= hs[nn - 1]) return { q: qChoke * N, surge: false, choke: true }; return { q: N * interp1(hr, qr, x), surge: false, choke: false }; };
  return { q: qs, h: hs, qSurge, qChoke, head, flow, map };
}
/** Surge margin (%) at a flow Q and speed N: 100 (Q − Q_surge)/Q on the same speed line. */
export const surgeMargin = (Q, N, qSurge) => (Q > 1e-12 ? (100 * (Q - qSurge * N)) / Q : -100);
/**
 * Variable-speed pump on a system curve H = hStatic + kSys Q²: rotor inertia with a torque-limited driver, flow PI acting on the driver
 * torque, quasi-steady hydraulics from the pump curve (affinity laws), shaft power P = P_r N³ (pShut + (1 − pShut) Q/(N Q_r)), trip and coast-down.
 * o: { qr (m³/s), hr (m), shutoff, eta, rho, J (kg·m²), rpm, hStatic (m), kSys (s²/m⁵; default puts the rated point on the system curve), curve,
 *      pShut (shut-off power fraction), tqMax (driver torque limit / rated torque), qSp(t) (m³/s), tTrip (s), kc, ti, nMin, nMax, dt, tEnd, n0 }
 * Returns { t, n[], q[], head[], power[] (W), torque: { rated }, pRated, tCoast50 (s from trip to half speed | null), qEnd, nEnd, iae, minFlowOk, analyticCoast(t) }.
 */
export function pumpSim(o) {
  const q = { shutoff: 1.25, eta: 0.75, rho: 800, J: 5, rpm: 3000, hStatic: 0, pShut: 0.4, tqMax: 1.5, kc: 2, ti: 4, nMin: 0, nMax: 1.05, dt: 0.02, tEnd: 60, tTrip: Infinity, n0: 1, ...o };
  const wr = (2 * Math.PI * q.rpm) / 60, pRated = (q.rho * G * q.qr * q.hr) / q.eta, tqR = pRated / wr, kSys = q.kSys ?? (q.hr - q.hStatic) / (q.qr * q.qr), hOf = (Q, N) => pumpHeadFn(Q, { qr: q.qr, hr: q.hr, shutoff: q.shutoff, speed: N, curve: q.curve || null });
  const flowAt = (N) => {
    if (N < 1e-6 || hOf(0, N) <= q.hStatic) return 0;
    if (!q.curve) return Math.sqrt(Math.max(q.hr * N * N * q.shutoff - q.hStatic, 0) / (kSys + (q.hr * (q.shutoff - 1)) / (q.qr * q.qr))); // parabolic curve: closed form
    const g = (Q) => hOf(Q, N) - q.hStatic - kSys * Q * Q; let hi = q.qr * Math.max(N, 0.05); for (let k = 0; k < 40 && g(hi) > 0; k++) hi *= 1.5; return brent(g, 0, hi, 1e-12);
  };
  const pwr = (Q, N) => pRated * N ** 3 * (q.pShut + (1 - q.pShut) * (Q / Math.max(N * q.qr, 1e-12)));
  const n = Math.round(q.tEnd / q.dt), out = { t: [], n: [], q: [], head: [], power: [] };
  let N = q.n0, I = null, iae = 0, tHalf = null, nTrip = null;
  for (let k = 0; k <= n; k++) {
    const t = k * q.dt, Q = flowAt(N), P = pwr(Q, N), tqL = N > 1e-6 ? P / (wr * N) : 0, tripped = t >= q.tTrip, sp = q.qSp ? q.qSp(t) : q.qr;
    if (I === null) I = clamp(tqL / (q.tqMax * tqR), 0, 1); // bumpless start at the present load
    let u = 0;
    if (!tripped) { const e = (sp - Q) / q.qr, raw = q.kc * e + I; u = clamp(raw, 0, 1); if (u === raw) I += (q.kc * e * q.dt) / q.ti; if (N >= q.nMax && u * q.tqMax * tqR > tqL) u = tqL / (q.tqMax * tqR); iae += Math.abs(e) * q.dt; }
    else { if (nTrip === null) nTrip = N; if (tHalf === null && N <= 0.5 * nTrip) tHalf = t - q.tTrip; }
    if (k % Math.max(1, Math.round(n / 600)) === 0 || k === n) { out.t.push(t); out.n.push(N); out.q.push(Q); out.head.push(hOf(Q, N)); out.power.push(P); }
    if (k < n) N = clamp(N + (q.dt * (u * q.tqMax * tqR - tqL)) / (q.J * wr), tripped ? 0 : q.nMin, 1.2);
  }
  const x0 = 1; // on a pure friction system the flow follows the speed (Q/(N Q_r) constant) and the torque falls with N²
  return { ...out, torque: { rated: tqR }, pRated, kSys, tCoast50: tHalf, qEnd: out.q[out.q.length - 1], nEnd: out.n[out.n.length - 1], iae, flowAt, power: out.power, powerAt: pwr, coastTime: (q.J * wr) / (tqR * (q.pShut + (1 - q.pShut) * x0)) };
}
/**
 * Centrifugal compressor between a suction drum and a discharge volume with speed and anti-surge control.
 * States: speed N (rotor inertia, torque- and power-limited driver), suction pressure, discharge pressure, recycle-valve position.
 * Flow through the machine follows the map at the head the pressure ratio demands (polytropic head with a constant efficiency); when the demanded
 * head exceeds the head at the surge point the flow collapses (a surge event is counted). The performance controller (PI) holds the suction
 * pressure with the speed; the anti-surge controller (PI on the surge margin with a fast-opening recycle valve) opens the cooled recycle when the
 * margin falls below the control line. Trip: driver torque to zero and recycle valve fully open.
 * o: { curve (compressorCurve), qd, hd, eta, ps (Pa), pExp (Pa export header), T1 (K), Z, k, mw (kg/mol), J, rpm, Vs, Vd (m³), feed(t) (kg/s), smCtl (%), kcAs, tiAs,
 *      tauOpen, tauClose (s), krFactor (recycle capacity / design flow), kcN, tiN, tqMax, pMax (W), nMin, nMax, tTrip, dt, tEnd, antiSurge (bool) }
 * Returns { t, n[], q[] (m³/s inlet), ps[], pd[] (bara), recycle[] (0–1), sm[] (%), power[] (W), design: { md, pd, power, torque, sm }, minSm (before the trip),
 *           surgeEvents, recycleMax, powerEnd, nEnd, tCoast50, energy (J) }.
 */
export function compressorSim(o) {
  const q = { eta: 0.78, Z: 0.9, k: 1.28, mw: 0.02, J: 60, rpm: 9000, Vs: 20, Vd: 10, smCtl: 10, kcAs: 0.08, tiAs: 6, tauOpen: 0.8, tauClose: 8, krFactor: 1.1, kcN: 4, tiN: 15, tqMax: 1.4, nMin: 0.7, nMax: 1.05, tTrip: Infinity, dt: 0.02, tEnd: 300, antiSurge: true, ...o };
  const cv = q.curve || compressorCurve({ qd: q.qd, hd: q.hd }), zrt = (q.Z * R * q.T1) / q.mw, nm = (q.k - 1) / (q.k * q.eta), rhoS0 = q.ps / zrt, md = rhoS0 * q.qd, hD = cv.head(q.qd, 1), ratio = (1 + (hD * nm) / zrt) ** (1 / nm), pd0 = q.ps * ratio;
  const wr = (2 * Math.PI * q.rpm) / 60, pD = (md * hD) / q.eta, tqR = pD / wr, pExp = q.pExp ?? 0.95 * pd0, Td = q.T1 * ratio ** nm, zrtD = (q.Z * R * Td) / q.mw, kd = md / Math.sqrt((pd0 / zrtD) * Math.max(pd0 - pExp, 1)), kr = (q.krFactor * md) / Math.sqrt((pd0 / zrtD) * (pd0 - q.ps));
  const headOf = (ps, pd) => (zrt / nm) * (Math.max(pd / ps, 1) ** nm - 1), pMax = q.pMax ?? 1.25 * pD, n = Math.round(q.tEnd / q.dt), every = Math.max(1, Math.round(n / 600));
  const out = { t: [], n: [], q: [], ps: [], pd: [], recycle: [], sm: [], power: [] };
  let N = clamp(q.n0 ?? 1, 0.3, 1.2), ps = q.ps, pd = pd0, xr = 0, In = tqR / (q.tqMax * tqR), Ias = 0, minSm = Infinity, surges = 0, wasSurge = false, rMax = 0, energy = 0, tHalf = null, nTrip = null;
  for (let k = 0; k <= n; k++) {
    const t = k * q.dt, tripped = t >= q.tTrip, H = headOf(ps, pd), fl = cv.flow(H, Math.max(N, 1e-6)), Q = fl.q, rhoS = ps / zrt, mc = rhoS * Q, sm = fl.surge ? -100 : surgeMargin(Q, N, cv.qSurge), P = (mc * H) / q.eta + 0.02 * pD * N ** 3;
    if (fl.surge && !wasSurge && !tripped) surges++; wasSurge = fl.surge;
    const rhoD = pd / zrtD, mr = xr * kr * Math.sqrt(rhoD * Math.max(pd - ps, 0)), mo = pd > pExp ? kd * Math.sqrt(rhoD * (pd - pExp)) : 0, feed = q.feed ? q.feed(t) : md;
    // performance controller: suction pressure above its set-point speeds the machine up
    let u = 0;
    if (!tripped) { const e = (ps - q.ps) / q.ps, raw = q.kcN * e + In; u = clamp(raw, 0, 1); if (u === raw || (raw > 1 && e < 0) || (raw < 0 && e > 0)) In += (q.kcN * e * q.dt) / q.tiN; minSm = Math.min(minSm, sm); energy += P * q.dt; }
    else { if (nTrip === null) nTrip = N; if (tHalf === null && N <= 0.5 * nTrip) tHalf = t - q.tTrip; }
    let tq = u * q.tqMax * tqR; if (tq * wr * N > pMax) tq = pMax / (wr * Math.max(N, 1e-6)); if (N >= q.nMax && tq > P / (wr * N)) tq = P / (wr * N);
    // anti-surge controller
    let xc = 0;
    if (tripped) xc = 1;
    else if (q.antiSurge) { const e = (q.smCtl - sm) / 100, raw = q.kcAs * e * 100 + Ias; xc = clamp(raw, 0, 1); if (xc === raw || (raw > 1 && e < 0) || (raw < 0 && e > 0)) Ias = Math.max(Ias + (q.kcAs * 100 * e * q.dt) / q.tiAs, -0.2); }
    if (k % every === 0 || k === n) { out.t.push(t); out.n.push(N); out.q.push(Q); out.ps.push(ps / 1e5); out.pd.push(pd / 1e5); out.recycle.push(xr); out.sm.push(sm); out.power.push(P); }
    if (!tripped) rMax = Math.max(rMax, xr);
    if (k < n) {
      const tqL = N > 1e-6 ? P / (wr * N) : 0;
      N = clamp(N + (q.dt * (tq - tqL)) / (q.J * wr), tripped ? 0 : q.nMin, 1.2);
      ps = Math.max(ps + ((q.dt * zrt) / q.Vs) * (feed + mr - mc), 0.05 * q.ps); pd = Math.max(pd + ((q.dt * zrtD) / q.Vd) * (mc - mr - mo), ps);
      xr += (q.dt * (xc - xr)) / (xc > xr ? q.tauOpen : q.tauClose); xr = clamp(xr, 0, 1);
    }
  }
  const last = out.t.length - 1;
  return { ...out, design: { md, pd: pd0 / 1e5, ratio, power: pD, torque: tqR, sm: surgeMargin(q.qd, 1, cv.qSurge), Td: Td - KEL, head: hD }, minSm, surgeEvents: surges, recycleMax: rMax, powerEnd: out.power[last], nEnd: out.n[last], tCoast50: tHalf, energy, curve: cv };
}

// ---- relaxation of flashing, chemical inventory, methanol partitioning -----------------------------------------
/**
 * Relaxation time (s) of the homogeneous relaxation model, dx/dt = (x_eq − x)/Θ (Downar-Zapolski et al. 1996, flashing water):
 * Θ = 6.51e-4 α^-0.257 ψ^-2.24 with ψ = (p_sat − p)/p_sat up to 10 bar, Θ = 3.84e-7 α^-0.54 φ^-1.76 with φ = (p_sat − p)/(p_crit − p_sat) above.
 */
export function hrmRelaxationTime(alpha, pSat, p, pCrit = 220.64e5) {
  if (!(pSat > p)) return Infinity;
  const a = clamp(alpha, 1e-6, 1);
  return pSat <= 10e5 ? 6.51e-4 * a ** -0.257 * ((pSat - p) / pSat) ** -2.24 : 3.84e-7 * a ** -0.54 * ((pSat - p) / Math.max(pCrit - pSat, 1)) ** -1.76;
}
/** Exact update of the relaxation closure over a step: unreleased mass E with equilibrium demand rate d (kg/s) and relaxation time Θ. Returns { E, released }. */
export function relaxStep(E, d, theta, dt) {
  if (!(theta > 0)) return { E: 0, released: E + d * dt };
  if (!Number.isFinite(theta)) return { E: E + d * dt, released: 0 };
  const f = Math.exp(-dt / theta), En = E * f + d * theta * (1 - f);
  return { E: En, released: E + d * dt - En };
}
/**
 * Chemical storage inventory over time: continuous use, batch withdrawals and deliveries triggered at the re-order level.
 * o: { V (tank m³), level0 (m³), use (m³/d), batches: [{ t (d), v (m³) }], reorder (m³; delivery ordered when the level falls to it), lead (d), delivery (m³), tEnd (d), dt (d) }
 * Returns { t[], level[], autonomy (d at the continuous use), tReorder (d | null), deliveries: [{ t, v }], minLevel, runOut (d | null), used, delivered, short (m³ that could not be supplied) }.
 */
export function chemicalInventory(o) {
  const q = { level0: o.V, use: 0, batches: [], reorder: 0, lead: 7, delivery: o.V, tEnd: 60, dt: 0.25, ...o }, n = Math.round(q.tEnd / q.dt), out = { t: [0], level: [q.level0] }, deliveries = [], bt = q.batches.map((b) => ({ ...b, done: false }));
  let lev = q.level0, due = null, tRe = null, minL = lev, runOut = null, used = 0, delivered = 0, short = 0;
  for (let k = 1; k <= n; k++) {
    const t = k * q.dt; let want = q.use * q.dt;
    for (const b of bt) if (!b.done && b.t <= t) { want += b.v; b.done = true; }
    const take = Math.min(want, lev); lev -= take; used += take; short += want - take; if (want > take + 1e-12 && runOut === null) runOut = t;
    if (due !== null && t >= due) { const v = Math.min(q.delivery, q.V - lev); lev += v; delivered += v; deliveries.push({ t, v }); due = null; }
    if (due === null && lev <= q.reorder && q.delivery > 0) { due = t + q.lead; if (tRe === null) tRe = t; }
    minL = Math.min(minL, lev); out.t.push(t); out.level.push(lev);
  }
  return { ...out, autonomy: q.use > 0 ? q.level0 / q.use : Infinity, tReorder: tRe, deliveries, minLevel: minL, runOut, used, delivered, short };
}
/**
 * Methanol vapour–aqueous K-value y/x (Moshfeghian's Wilson-type correlation fitted to data for −23…38 °C and 7…345 bar):
 * K = exp[5.37 (1 + ω*) (1 − 1/T*)] / P*, P* = P[psia]/35, T* = T[°R]/615, ω* = 2.95 − 0.02607 P* + 8.92828e-5 P*² − 0.851257/T*.
 */
export function methanolK(Pbar, Tc) {
  const Ps = (Pbar * 14.5038) / 35, Ts = ((Tc + KEL) * 1.8) / 615, w = 2.95 - 0.02607 * Ps + 8.92828e-5 * Ps * Ps - 0.851257 / Ts;
  return Math.exp(5.37 * (1 + w) * (1 - 1 / Ts)) / Ps;
}

// ---- suite-level assembly -------------------------------------------------------------------------------------------
const STEEL = { k: 45, rho: 7850, cp: 470 };
const INH_OPTS = ['MeOH', 'MEG', 'EtOH', 'DEG', 'TEG'];
/** Resample a flow picture to n cells of equal arc length. */
function buildStations(pic, n, A) {
  const S = [0]; for (let i = 1; i < pic.x.length; i++) S.push(S[i - 1] + Math.hypot(pic.x[i] - pic.x[i - 1], pic.z[i] - pic.z[i - 1]));
  const L = S[S.length - 1], ds = L / n, st = { s: [], x: [], z: [], dz: [], ds, A, L, S, n };
  const keys = ['P', 'T', 'holdup', 'tAmb', 'vm', 'vsl', 'vsg', 'rhoM', 'tHyd'];
  for (const k of keys) st[k] = [];
  for (let i = 0; i < n; i++) { const sc = (i + 0.5) * ds; st.s.push(sc); st.x.push(interp1(S, pic.x, sc)); st.z.push(interp1(S, pic.z, sc)); st.dz.push(interp1(S, pic.z, (i + 1) * ds) - interp1(S, pic.z, i * ds)); for (const k of keys) st[k].push(interp1(S, pic[k], sc)); }
  return st;
}
/** Wall and coating layers (bore outwards, SI) from the inputs, with the coating conductivity scaled to the stated U-value. */
function wallLayers(v, id, wt, hFlow, hOut, warnings) {
  const user = (Array.isArray(v.layers) ? v.layers : []).map((l) => ({ name: String(l.name || 'Layer'), t: num(l.t, 0) / 1000, k: num(l.k, 0), rho: num(l.rho, 900) || 900, cp: num(l.cp, 1700) || 1700 })).filter((l) => l.t > 0 && l.k > 0);
  const uOf = (f) => uValue({ id, wt, kWall: STEEL.k, layers: user.map((l) => ({ t: l.t, k: l.k * f })), hIn: hFlow, hOut }).U;
  let f = 1;
  if (v.uValue > 0 && user.length) { const g = (x) => uOf(x) - v.uValue; if (g(0.02) < 0 && g(50) > 0) f = brent(g, 0.02, 50, 1e-9); else warnings.push({ level: 'warn', msg: `The stated U-value of ${v.uValue} W/m²K cannot be reached with these layers; the layers are used as entered (U = ${uOf(1).toFixed(2)} W/m²K).` }); }
  f *= v.uMult;
  const layers = [{ name: 'Steel wall', t: wt, ...STEEL }, ...user.map((l) => ({ ...l, k: l.k * f }))];
  let r = id / 2; const capOf = (l) => { const c = l.rho * l.cp * Math.PI * ((r + l.t) ** 2 - r * r); r += l.t; return c; }, caps = layers.map(capOf), od = 2 * r;
  if (v.thermalMass > 0 && caps.length > 1) { const target = v.thermalMass * 1000, coat = sum(caps.slice(1)); if (target > caps[0] && coat > 0) layers.slice(1).forEach((l) => { l.rho *= (target - caps[0]) / coat; }); else warnings.push({ level: 'info', msg: 'The stated thermal mass is below that of the bare steel; the layer properties are used instead.' }); }
  layers.forEach((l) => { l.rho *= v.cMult; });
  if (v.burial > 0) { const depth = Math.max(v.burial, 0.55 * od); layers.push({ name: 'Soil (equivalent annulus)', t: (od / 2) * (Math.exp(Math.acosh((2 * depth) / od)) - 1), k: Math.max(v.kSoil, 0.1), rho: 1900, cp: 1300 }); }
  return { layers, factor: f, od, uFlow: v.burial > 0 ? uValue({ id, wt, kWall: STEEL.k, layers: layers.slice(1), hIn: hFlow, hOut: 1e4 }).U : uOf(f) };
}
/** Parameters of the slugging model from the steady picture at one rate. */
function slugParams(st, fm, line, v, rate, base) {
  const sRb = clamp(interp1(st.xNodes, st.S, line.riserBaseX), 0.3 * st.L, 0.98 * st.L), iF = st.s.map((s, i) => i).filter((i) => st.s[i] < sRb), iR = st.s.map((s, i) => i).filter((i) => st.s[i] >= sRb);
  const avg = (key, idx) => mean(idx.map((i) => st[key][i])), Pp = avg('P', iF), Tp = avg('T', iF), Tr = iR.length ? avg('T', iR) : Tp, k0 = iF[iF.length - 1], pr = fm.at(st.P[k0], st.T[k0], rate), pf = fm.at(Pp, Tp, rate);
  // feed inclination: mean downward slope over the last fifth of the flowline (never flatter than 0.1°)
  const j0 = iF[Math.floor(0.8 * iF.length)], slope = Math.abs(st.z[k0] - st.z[j0]) / Math.max(st.s[k0] - st.s[j0], 1), zTop = st.z[st.n - 1] + 0.5 * st.dz[st.n - 1], zRb = interp1(st.S, st.zNodes, sRb);
  return { D: Math.sqrt((4 * st.A) / Math.PI), Lp: sRb, Vp: st.A * sRb, Lr: Math.max(zTop - zRb, 5), Vr: st.A * (st.L - sRb + Math.max(v.topsideLen, 5)), theta: Math.max(slope, 1.75e-3), rhoL: pr.rhoL, mwG: pf.mwG / 1000, Z: pf.zG, Tp: Tp + KEL, Tr: Tr + KEL, muL: pr.muL,
    wG: Math.max(pr.mG, 1e-3), wL: Math.max(pr.mO + pr.mW, 1e-3), Ps: v.sepP * 1e5, Cv: v.chokeCv, chokeExp: clamp(v.chokeExp, 0.3, 4), kH: v.kH, kL: v.kL, aLp: clamp(avg('holdup', iF), 0.08, 0.85), rhoGnom: pf.rhoG, rough: line.roughness, ...(base ? { kG: base.kG } : {}) };
}
/** Replay of an operating log through the quasi-steady surrogate with a first-order thermal lag; optional ridge residual correction. */
function replayLog(rows, m) {
  const t = rows.map((r) => r.t * HOUR), n = rows.length, pIn = [], tArr = []; let T = null;
  for (let i = 0; i < n; i++) {
    const q = clamp(rows[i].rate / 100, 0, m.qHi), flowing = q > 0.02, z = clamp((rows[i].choke || 100) / 100, 0.02, 1), dt = i ? Math.max(t[i] - t[i - 1], 0) : 0;
    pIn.push(flowing ? m.pIn(Math.max(q, m.qLo)) + m.dpChoke(Math.max(q, m.qLo), z) : m.pSettle);
    const target = flowing ? m.tArr(Math.max(q, m.qLo)) : m.tAmbOut, tau = flowing ? m.tauWarm / Math.max(q, 0.1) : m.tauCool;
    T = T === null ? target : target + (T - target) * Math.exp(-dt / Math.max(tau, 60)); tArr.push(T);
  }
  return { t: rows.map((r) => r.t), pIn, tArr };
}
function ridgeResidual(rows, pred, key, lam = 1e-3) {
  const X = rows.map((r) => [1, r.rate / 100, (r.choke || 100) / 100]), y = rows.map((r, i) => r[key] - pred[i]), tr = X.map((_, i) => i).filter((i) => i % 2 === 0), te = X.map((_, i) => i).filter((i) => i % 2 === 1);
  if (tr.length < 4 || te.length < 2) return null;
  const Xt = tr.map((i) => X[i]), A = zeros(3, 3), g = [0, 0, 0];
  for (let i = 0; i < Xt.length; i++) for (let a = 0; a < 3; a++) { g[a] += Xt[i][a] * y[tr[i]]; for (let b = 0; b < 3; b++) A[a][b] += Xt[i][a] * Xt[i][b]; }
  for (let a = 0; a < 3; a++) A[a][a] += lam * Xt.length;
  let coef; try { coef = solveLinear(A, g); } catch { return null; }
  const corr = X.map((x) => dot(coef, x)), rm = (idx, c) => Math.sqrt(mean(idx.map((i) => (y[i] - (c ? corr[i] : 0)) ** 2)));
  return { coef, corr, before: rm(te, false), after: rm(te, true), trainBefore: rm(tr, false), trainAfter: rm(tr, true) };
}

async function run(v0, ctx = {}) {
  const v = { ...DEFAULTS, ...v0 }, warnings = [], rec = [], progress = (f, m) => ctx.progress?.(f, m), tick = async () => { if (ctx.tick) await ctx.tick(); };
  for (const k of NUMERIC) { v[k] = num(v[k], DEFAULTS[k]); }
  const need = (cond, msg) => { if (!cond) throw new Error(msg); };
  need(v.idMm >= 25 && v.idMm <= 1500, 'The inner diameter must lie between 25 and 1500 mm.');
  need(v.wtMm >= 1 && v.wtMm <= 100, 'The wall thickness must lie between 1 and 100 mm.');
  need(v.rateFrac >= 10 && v.rateFrac <= 200, 'The operating rate must lie between 10 and 200 % of the case rate.');
  need(v.tHorizon > 0.5 && v.tShut >= 0, 'The cooldown horizon must exceed half an hour and the shutdown duration cannot be negative.');
  need(v.orificeMm > 0.5 && v.cdBlow > 0.05 && v.cdBlow <= 1, 'The blowdown orifice needs a positive diameter and a discharge coefficient between 0.05 and 1.');
  need(v.sepP >= 1.2 && v.pBack >= 1, 'The separator pressure must be at least 1.2 bara and the flare back-pressure at least 1 bara.');
  need(v.slugCatcherVol > 0.5 && v.qDrainM3h > 0, 'The slug-catcher volume and the liquid handling capacity must be positive.');
  need(v.uMult > 0.05 && v.cMult > 0.05, 'The U-value and thermal-mass multipliers must be positive.');
  need(v.pigRatePct >= 10 && v.pigRatePct <= 150, 'The pigging rate must lie between 10 and 150 % of the case rate.');
  const nSt = Math.round(clamp(v.nStations, 8, 120)), nr = Math.round(clamp(v.nr, 3, 80)), ntCool = Math.round(clamp(v.ntCool, 12, 2000)), ntBlow = Math.round(clamp(v.ntBlow, 30, 6000)), nEnv = Math.round(clamp(v.nEnv, 4, 16)), dtCtl = clamp(v.dtCtl, 2, 120);
  const rate = v.rateFrac / 100, inhId = INH_OPTS.includes(v.inhibitor) ? v.inhibitor : 'MeOH', inh = INHIBITORS[inhId];

  // ---------- steady starting point ----------
  progress(0.02, 'Steady starting point');
  const fm = fluidModel(ctx), id = v.idMm / 1000, wt = v.wtMm / 1000, A = (Math.PI * id * id) / 4, over = { id, wt, ...(v.uValue > 0 ? { uValue: v.uValue * v.uMult } : {}), tSeabed: v.tSeabed, tSeaSurface: v.tSurface, pOut: v.sepP };
  const line = caseLine(ctx, over);
  // steady solution at another rate: kernel marching with a bracketed regula falsi on the inlet pressure (outlet met within 0.02 bar)
  let pGuess = null;
  const solve = (mScale, n = 36) => {
    const march = (pIn) => marchSteady({ fm, profile: line.profile, id: line.id, rough: line.roughness, U: line.uValue, tAmbOf: (s, z) => ambientAt(z, line), tIn: line.tIn, pIn, mScale, n }), res = (r) => (r.ok ? r.pOut - line.pOut : -1e3);
    let a = Math.max(pGuess ?? line.pOut + 60, line.pOut + 1), ra = march(a), fa = res(ra), b, rb, fb;
    if (Math.abs(fa) < 0.02) { pGuess = a; ra.line = line; ra.fm = fm; return ra; }
    for (let k = 0, step = 0.12 * a; ; k++, step *= 1.8) { b = fa < 0 ? a + step : Math.max(a - step, line.pOut + 0.2); rb = march(b); fb = res(rb); if (fa * fb < 0) break; if (k > 16 || b > 1400) throw new Error('no inlet pressure delivers this rate to the outlet'); a = b; ra = rb; fa = fb; }
    for (let k = 0; k < 40; k++) {
      const useSecant = fa > -900 && fb > -900, c = useSecant ? clamp(a - (fa * (b - a)) / (fb - fa), Math.min(a, b) + 0.02 * Math.abs(b - a), Math.max(a, b) - 0.02 * Math.abs(b - a)) : 0.5 * (a + b), rc = march(c), fc = res(rc);
      if (Math.abs(fc) < 0.02) { pGuess = c; rc.line = line; rc.fm = fm; return rc; }
      if (fa * fc < 0) { b = c; rb = rc; fb = fc; fa = useSecant ? fa * 0.5 : fa; } else { a = c; ra = rc; fa = fc; fb = useSecant ? fb * 0.5 : fb; }
    }
    throw new Error('the inlet-pressure iteration did not converge');
  };
  let pic;
  try { pic = Math.abs(rate - 1) < 1e-9 ? flowPicture(ctx, over) : { ...solve(rate, 100), source: 'kernel estimate' }; } catch (e) { throw new Error(`No steady flow solution at ${v.rateFrac} % of the case rate: ${e.message}`); }
  const st = buildStations(pic, nSt, A); st.xNodes = pic.x; st.zNodes = pic.z;
  const L = st.L, mdot = (fm.rates.mHC + fm.rates.mW) * rate, volScale = v.lineVolume > 0 ? clamp(v.lineVolume / (A * L), 0.3, 3) : 1, lineVol = A * L * volScale;
  const iMid = Math.floor(nSt / 2), pm = fm.at(st.P[iMid], st.T[iMid], rate), cpMix = (pm.mG * pm.cpG + pm.mO * pm.cpO + pm.mW * pm.cpW) / Math.max(pm.mG + pm.mO + pm.mW, 1e-9);
  const hFlow = hInside((Math.max(st.rhoM[iMid], 1) * Math.max(st.vm[iMid], 0.01) * id) / pm.muL, (pm.cpL * pm.muL) / pm.kL, pm.kL, id);
  const odGuess = id + 2 * wt + 2 * sum((v.layers || []).map((l) => num(l.t, 0) / 1000)), hSea = hOutside(v.currentSpeed, odGuess, 'seawater', v.tSeabed), hAir = hOutside(5, odGuess, 'air');
  const wall = wallLayers(v, id, wt, hFlow, hSea, warnings), buried = v.burial > 0;
  const pIn0 = pic.P[0], tArr0 = pic.T[pic.T.length - 1], liqInvSteady = sum(st.holdup) * st.ds * A * volScale;
  const tAmbOf = (i) => (st.z[i] >= 0 ? v.tAir : st.tAmb[i]);
  await tick();

  // ---------- A. shutdown and cooldown ----------
  progress(0.1, 'Settle-out and cooldown');
  const so = settleOut(st, fm), dtCool = (v.tHorizon * HOUR) / ntCool, grids = new Map();
  const gridOf = (h) => { const key = Math.round(h); if (!grids.has(key)) grids.set(key, radialGrid({ ri: id / 2, layers: wall.layers, nr, hIn: v.hInShut, hOut: h })); return grids.get(key); };
  const cool = [], eBal = { drop: 0, lost: 0 };
  for (let i = 0; i < nSt; i++) {
    const p = fm.at(so.pSettle, st.T[i]), H = so.holdup[i], cFluid = A * (H * p.rhoL * p.cpL + (1 - H) * p.rhoG * p.cpG);
    const c = cooldown({ grid: gridOf(buried ? 1e4 : st.z[i] >= 0 ? hAir : hSea), ri: id / 2, T0: st.T[i], tAmb: tAmbOf(i), cFluid, dt: dtCool, nSteps: ntCool });
    cool.push(c); eBal.drop += (c.energy.stored0 - c.energy.stored) * st.ds; eBal.lost += c.energy.lost * st.ds;
  }
  const tH = cool[0].t.map((t) => t / HOUR), vg = so.holdup.map((h) => (1 - h) * st.ds * A), vg0 = sum(vg.map((x, i) => x / (st.T[i] + KEL)));
  const pT = tH.map((_, k) => (vg0 > 0 ? (so.pSettle * vg0) / sum(vg.map((x, i) => x / (cool[i].Tf[k] + KEL))) : so.pSettle)), thT = pT.map((p) => fm.hydrateT(p)), limT = thT.map((x) => x + v.hydMargin);
  const tCool = cool.map((c) => timeBelow(tH, c.Tf, limT)), tWat = cool.map((c) => timeBelow(tH, c.Tf, v.wat)), tPour = cool.map((c) => timeBelow(tH, c.Tf, v.pourPoint));
  // special components: lumped capacitance with their own U-value and thermal mass
  const comps = (Array.isArray(v.components) ? v.components : []).filter((c) => num(c.u, 0) > 0).map((c) => {
    const i = clamp(Math.round((interp1(pic.x, st.S, clamp(num(c.x, 0), pic.x[0], pic.x[pic.x.length - 1])) / L) * nSt - 0.5), 0, nSt - 1), ta = tAmbOf(i), T0 = ta + clamp(num(c.t0f, 1), 0.05, 1) * (st.T[i] - ta), tau = (Math.max(num(c.mass, 1), 0.05) * cool[i].cEff) / (num(c.u, 1) * Math.PI * id);
    const series = tH.map((t) => ta + (T0 - ta) * Math.exp((-t * HOUR) / tau));
    return { name: String(c.name || 'Component'), x: st.x[i], u: num(c.u, 1), tau: tau / HOUR, T0, t: timeBelow(tH, series, limT), series };
  });
  const finite = (a) => a.filter((x) => x !== null), lineMin = finite(tCool).length ? Math.min(...finite(tCool)) : null, compMin = finite(comps.map((c) => c.t)).length ? Math.min(...finite(comps.map((c) => c.t))) : null;
  const iCold = lineMin === null ? tCool.reduce((b, _, i) => (cool[i].Tf[ntCool] - limT[ntCool] < cool[b].Tf[ntCool] - limT[ntCool] ? i : b), 0) : tCool.indexOf(lineMin);
  const coldComp = compMin !== null && (lineMin === null || compMin < lineMin) ? comps.find((c) => c.t === compMin) : null;
  const neverCools = lineMin === null && compMin === null, cooldownTime = neverCools ? v.tHorizon : Math.min(lineMin ?? Infinity, compMin ?? Infinity), coldSpotX = coldComp ? coldComp.x : st.x[iCold];
  const atTime = (c, h) => interp1(tH, c.Tf, clamp(h, 0, v.tHorizon)), tAtShut = cool.map((c) => atTime(c, v.tShut)), pAtShut = interp1(tH, pT, Math.min(v.tShut, v.tHorizon));
  const lumpedCold = cool[iCold].lumped(12 * HOUR), cold12 = atTime(cool[iCold], Math.min(12, v.tHorizon));
  if (v.tShut > v.tHorizon) warnings.push({ level: 'info', msg: `The shutdown (${v.tShut} h) is longer than the simulated cooldown (${v.tHorizon} h); the state at ${v.tHorizon} h is used for the restart.` });
  await tick();

  // ---------- envelope scan (steady solutions over rate) — also feeds the ramp-up, pigging and optimisation ----------
  progress(0.22, 'Rate scan');
  const scan = { q: [], pIn: [], tArr: [], margin: [], eros: [], qLiq: [], inv: [], vMean: [], zCrit: [], dpChoke: [], pReq: [], res: [] };
  const sm0 = slugModel(slugParams(st, fm, line, v, rate)), Kc = 2.403e-5 * v.chokeCv, qEnv = linspace(v.qLoPct / 100, v.qHiPct / 100, nEnv), qPig = v.pigRatePct / 100;
  { const k = qEnv.reduce((b, q, i) => (Math.abs(q - qPig) < Math.abs(qEnv[b] - qPig) ? i : b), 0); if (k > 0 && k < nEnv - 1 && qPig > qEnv[k - 1] && qPig < qEnv[k + 1]) qEnv[k] = qPig; } // the pigging rate becomes one of the scan points
  for (let k = 0; k < nEnv; k++) {
    let r; try { r = solve(qEnv[k]); } catch { continue; }
    const n = r.P.length, s2 = buildStations(r, Math.min(nSt, 30), A); s2.xNodes = r.x; s2.zNodes = r.z;
    let zc = null, dpc = 0;
    try { const smk = slugModel(slugParams(s2, fm, line, v, qEnv[k], sm0.p)); zc = smk.critical(smk.zFloor, 1, 9); const zUse = zc === null ? 1 : Math.min(1, v.slugControl ? 2 * zc : 0.9 * zc), eq = smk.steady(zUse); dpc = (eq.Prt - smk.p.Ps) / 1e5; } catch { zc = null; }
    scan.q.push(qEnv[k]); scan.pIn.push(r.pIn); scan.tArr.push(r.tOut); scan.margin.push(-Math.max(...r.subcooling)); scan.eros.push(Math.max(...r.vm.map((x, i) => x / (122 / Math.sqrt(Math.max(r.rhoM[i], 1))))));
    scan.qLiq.push(r.qL[n - 1] * 3600); scan.inv.push(r.liquidInventory * volScale); scan.vMean.push(r.length / r.residence); scan.zCrit.push(zc === null ? 1 : zc); scan.dpChoke.push(dpc); scan.pReq.push(r.pIn + dpc); scan.res.push(r);
    if (k % 3 === 2) await tick();
  }
  need(scan.q.length >= 3, 'Fewer than three rates of the scan have a steady solution; widen the scan range or check the case data.');
  const qLo = scan.q[0], qHi = scan.q[scan.q.length - 1], qClamp = (q) => clamp(q, qLo, qHi), at = (key, q) => interp1(scan.q, scan[key], qClamp(q));

  // ---------- solids state (solids suite) → action triggers ----------
  const depWax = Math.max(v.depWaxMm, 0), depHyd = Math.max(v.depHydMm, 0), depScale = Math.max(v.depScaleMm, 0), triggers = [], pigCycle = v.depRateMmD > 0 ? v.maxDepositMm / v.depRateMmD : null, daysToPig = v.depRateMmD > 0 ? Math.max(0, (v.maxDepositMm - depWax) / v.depRateMmD) : null;
  if (depWax >= v.maxDepositMm) triggers.push({ id: 'pig', tag: 'XA-150', msg: `Wax deposit ${depWax.toFixed(1)} mm has reached the ${v.maxDepositMm} mm limit`, action: 'Launch a pig as soon as the line is steady' });
  if (depScale >= v.maxDepositMm) triggers.push({ id: 'scale', tag: 'XA-151', msg: `Scale deposit ${depScale.toFixed(1)} mm has reached the ${v.maxDepositMm} mm limit`, action: 'Start scale-inhibitor injection and plan a squeeze or mechanical cleaning' });
  if (v.hydFrac >= v.maxHydFrac) triggers.push({ id: 'hydrate', tag: 'XA-152', msg: `Hydrate accumulation (${(100 * v.hydFrac).toFixed(1)} vol % in the liquid, ${depHyd.toFixed(1)} mm on the wall) has reached its limit of ${(100 * v.maxHydFrac).toFixed(0)} vol %`, action: 'Inject hydrate inhibitor; no shutdown or restart without preservation' });
  if (v.blockagePct >= v.maxBlockagePct) triggers.push({ id: 'blockage', tag: 'XA-153', msg: `Bore restriction ${v.blockagePct.toFixed(0)} % has reached the ${v.maxBlockagePct} % limit`, action: 'Reduce the rate and remediate by heating or by depressurising from both ends' });
  const trig = (id) => triggers.some((t) => t.id === id), pigIv = pigCycle !== null ? +(v.pigInterval > 0 ? Math.min(v.pigInterval, pigCycle) : pigCycle).toFixed(2) : v.pigInterval;

  // ---------- design surges and surge capacity (needed by the separator model of the control study) ----------
  const q0 = clamp(v.qStartPct / 100, 0.05, 1), rampS = Math.max(v.rampHours, 0.05) * HOUR, invOf = (q) => at('inv', q), qLiqOf = (q) => (at('qLiq', q) / 3600) * (q / qClamp(q)), qDrain = v.qDrainM3h / 3600;
  const surgeFor = (hours) => rampSurge({ nodes: [[0, q0 * rate], [hours * HOUR, rate]], inv: invOf, qLiq: qLiqOf, qDrain, inv0: liqInvSteady, dt: 120 });
  const ramp = surgeFor(v.rampHours);
  const pigDpFric = Math.max((4 * v.pigFric * v.pigContact * 1e5 * v.pigSealLen) / id, 100); // seal contact force: μ × contact pressure × seal area
  let rPig = scan.res[scan.q.indexOf(qPig)]; try { if (!rPig) rPig = solve(qPig); } catch (e) { throw new Error(`No steady flow solution at the pigging rate (${v.pigRatePct} %): ${e.message}`); }
  const pigS0 = v.pigX0 > 0 ? clamp(interp1(rPig.x, rPig.s, clamp(v.pigX0, rPig.x[0], rPig.x[rPig.x.length - 1])), 0, 0.999 * rPig.s[rPig.s.length - 1]) : 0, pigInLine = pigS0 > 0;
  const pigOf = (r, fromStart = false) => pigRun({ s: r.s, z: r.z, vm: r.vm, holdup: r.holdup, vsl: r.vsl, rhoG: r.rhoG, rhoM: r.rhoM, rhoL: mean(r.rhoL), D: id, fric: pigDpFric, mass: v.pigMass, bypass: v.pigBypass / 100, cdBypass: v.pigCd, leak: v.pigLeak / 100, qDrain, ...(pigInLine && !fromStart ? { s0: pigS0, slug0: Math.max(v.pigSlug0, 0) / volScale } : {}) });
  const pig = pigOf(rPig);
  const usable = clamp(v.catcherUsable, 5, 100) / 100, designSurge = Math.max(ramp.vMax, pig.surge * volScale, v.slugSurge), catReq = (1.1 * designSurge) / usable, catVol = v.catcherAuto ? Math.max(v.slugCatcherVol, catReq) : v.slugCatcherVol, surgeAllow = catVol * usable;
  if (v.catcherAuto && catReq > v.slugCatcherVol) warnings.push({ level: 'info', msg: `Surge capacity sized by the design: the largest design surge (${designSurge.toFixed(1)} m³) needs ${catReq.toFixed(0)} m³ with ${(100 * usable).toFixed(0)} % of the volume usable between the normal and the high level; the entered ${v.slugCatcherVol.toFixed(1)} m³ is replaced by that volume.` });

  // ---------- G. slugging and control ----------
  progress(0.4, 'Slugging model and control');
  const zCrit = sm0.critical(sm0.zFloor, 1), zAuto = zCrit === null ? 1 : Math.min(1, 2 * zCrit), sensorFailed = v.sensorState === 'failed', manual = sensorFailed || !v.ctlAuto, zCmd = clamp(v.chokePct > 0 ? v.chokePct / 100 : zAuto, 0.03, 1);
  // with the loop in manual (operator choice or a failed inlet-pressure transmitter) the choke is held at an opening that is stable without feedback
  const zTarget = clamp(manual && zCrit !== null ? Math.min(zCmd, 0.9 * zCrit) : zCmd, Math.max(v.zMinPct / 100, 0.03, sm0.zFloor), Math.max(v.zMaxPct / 100, 0.05, sm0.zFloor)), unstable = zCrit !== null && zTarget > zCrit;
  const lin = sm0.linearise(zTarget), eqT = lin.steady, polesOL = sm0.poles(zTarget), kStat = (sm0.steady(Math.min(1, zTarget * 1.02)).Pp - sm0.steady(zTarget * 0.98).Pp) / (Math.min(1, zTarget * 1.02) - zTarget * 0.98) / 1e5;
  // step test at a stable opening → first-order-plus-dead-time model → tuning rules
  const zId = Math.max(zCrit === null ? 0.5 * zTarget : 0.6 * zCrit, sm0.zFloor), eqId = sm0.steady(zId), decay = Math.max(-sm0.growth(zId), 1e-5), tId = clamp(6 / decay, 0.5 * HOUR, 16 * HOUR), dz = 0.1 * zId;
  const stepRun = integrateStiff(sm0.f(zId + dz), eqId.y, 0, tId, { rtol: 1e-4, atol: 1e-3, hInit: 5, hMax: tId / 150, maxSteps: 4000 });
  const loopDead = v.deadTime + v.actDead, stepY = stepRun.y.map((y) => sm0.alg(y, zId + dz).Pp / 1e5), fo = identifyFOPDT(stepRun.t, stepY, dz), foUse = { K: fo.K, tau: Math.max(fo.tau, 1), theta: Math.max(fo.theta, loopDead, dtCtl) };
  const rules = tuningRules(foUse, v.tauCFactor === 1 ? undefined : v.tauCFactor * foUse.theta), gainRatio = kStat !== 0 ? fo.K / kStat : 1, pole = tuneByPoles(lin, { tis: v.ctlMode === 'P' ? [0] : [1800, 5400, 14400], theta: loopDead + dtCtl });
  // robust PI: largest integral gain with the sensitivity peak below the limit for every model of the uncertainty set (gain and dead time)
  const gU = clamp(v.robGainPct, 0, 90) / 100, dU = clamp(v.robDelayPct, 0, 300) / 100, fopdt = (K, th) => ({ k: K, lags: [foUse.tau], delay: th }), robSet = [fopdt(foUse.K, foUse.theta), fopdt(foUse.K * (1 + gU), foUse.theta * (1 + dU)), fopdt(foUse.K * (1 - gU), foUse.theta * (1 + dU)), fopdt(foUse.K * (1 + gU), foUse.theta)];
  const rob = robustPI(robSet, { msMax: clamp(v.msMax, 1.1, 4) }); if (rob) rules.rows.push({ rule: `Robust multi-model (Ms ≤ ${clamp(v.msMax, 1.1, 4)})`, mode: 'PI', kc: rob.kc, ti: rob.ti, td: 0 });
  // recursive least squares on the step test → first-order model → self-tuned SIMC PI
  const nU = 160, tU = linspace(0, tId, nU), yU = tU.map((t) => interp1(stepRun.t, stepY, t) - stepY[0]), rl = rlsFirstOrder(new Array(nU).fill(dz), yU, tId / (nU - 1), { nd: Math.min(Math.round(Math.max(fo.theta, 0) / (tId / (nU - 1))), nU - 20) }), rlsRule = tuningRules({ K: rl.K, tau: Math.max(rl.tau, 1), theta: foUse.theta }).rows[1];
  if (isNum(rlsRule.kc) && isNum(rlsRule.ti)) rules.rows.push({ ...rlsRule, rule: 'Self-tuning (RLS model + SIMC)' });
  const tuneRows = rules.rows.map((r) => ({ ...r, kcT: r.kc * gainRatio, stable: pole.stable(r.kc * gainRatio, r.ti) }));
  let sel;
  if (v.tuning === 'manual') sel = { rule: 'Manual', kc: v.kcMan / 100, ti: v.ctlMode === 'P' ? 0 : v.tiMan, td: v.ctlMode === 'PID' ? v.tdMan : 0 };
  else if (v.tuning === 'robust' || v.tuning === 'rls') { const r = tuneRows.find((x) => x.rule.startsWith(v.tuning === 'robust' ? 'Robust' : 'Self-tuning')) || tuneRows[1]; sel = { rule: r.rule + ' (gain-scheduled)', kc: r.kcT, ti: v.ctlMode === 'P' ? 0 : r.ti, td: 0 }; }
  else if (v.tuning === 'simc' || v.tuning === 'zn') { const r = tuneRows.find((x) => x.mode === v.ctlMode && x.rule.startsWith(v.tuning === 'simc' ? 'SIMC' : 'Ziegler')); sel = { rule: r.rule + ' (gain-scheduled)', kc: r.kcT, ti: r.ti, td: r.td }; }
  else sel = { rule: 'Closed-loop pole search', kc: pole.kc, ti: v.ctlMode === 'P' ? 0 : pole.ti, td: v.ctlMode === 'PID' ? Math.max(v.deadTime, dtCtl) / 3 : 0 };
  let linStable = pole.stable(sel.kc, sel.ti);
  if (!linStable && v.tuning !== 'manual') { warnings.push({ level: 'warn', msg: `${sel.rule} settings do not stabilise the linearised model at ${(100 * zTarget).toFixed(0)} % opening; the pole-search settings are applied instead.` }); sel = { rule: 'Closed-loop pole search (fallback)', kc: pole.kc, ti: v.ctlMode === 'P' ? 0 : pole.ti, td: 0 }; linStable = pole.stable(sel.kc, sel.ti); }
  // gain scheduling on the choke opening: the controller gain follows the inverse of the local static gain of the riser (kept between 0.2 and 5 times the design value)
  const zSch = linspace(Math.max(sm0.zFloor, 0.03), 1, 8), kSch = zSch.map((z) => { try { return (sm0.steady(Math.min(1, z * 1.02)).Pp - sm0.steady(z * 0.98).Pp) / (Math.min(1, z * 1.02) - z * 0.98) / 1e5; } catch { return kStat; } }), schedule = v.adaptive === 'schedule' && kStat !== 0 ? (z) => { const k = interp1(zSch, kSch, clamp(z, zSch[0], 1)); return k * kStat > 0 ? kStat / k : 1; } : null;
  const tCtl = clamp(v.tCtl, 1, 72) * HOUR, ctlCfg = { dt: dtCtl, zTarget, tEnd: tCtl, tOff: 0.55 * tCtl, tStep: 0.2 * tCtl, dSp: v.spStep, kc: sel.kc, ti: sel.ti, td: sel.td, rate: v.rateLimit / 100, dead: loopDead, tauValve: Math.max(v.tauValve, 0.5), tauSensor: v.tauSensor, schedule, noise: v.noiseBar, cascade: !!v.cascade, override: !!v.override, feedForward: !!v.feedForward, antiWindup: true, sepV: catVol, levelSp: clamp(v.levelSp, 10, 90) / 100, levelHi: clamp(v.levelHi, 20, 99) / 100, bias: v.sensorState === 'bias' ? v.sensorBias : 0, ...(v.pSet > 0 ? { sp: v.pSet * 1e5 } : {}), qDrain: v.qDrainM3h / 3600, zMin: v.zMinPct / 100, zMax: v.zMaxPct / 100, np: Math.round(clamp(v.mpcNp, 3, 60)), nc: Math.round(clamp(v.mpcNc, 1, 12)), rDu: v.mpcR, tsMpc: Math.max(dtCtl, 4 * dtCtl) };
  const pidRun = slugControl(sm0, { ...ctlCfg, mode: manual ? 'open' : 'pid' });
  await tick();
  const mpcRun = slugControl(sm0, { ...ctlCfg, mode: 'mpc', tEnd: ctlCfg.tOff, tOff: ctlCfg.tOff + 1 });
  // extended Kalman filter on the nonlinear model: riser liquid mass from the two noisy pressures
  // the filters follow the closed-loop part of the record and the first oscillations after the switch to manual (at most 1000 samples)
  const stride = Math.max(1, Math.ceil(pidRun.t.length / 1500)), hK = dtCtl * stride, rnd = rng(23), idx = pidRun.t.map((_, i) => i).filter((i) => i % stride === 0 && i > 0).slice(0, 1000), g4 = 1 + Math.SQRT1_2;
  const sc4 = eqT.y.map((x) => Math.abs(x) || 1), sig4 = [0.01 * sc4[0], 0.05 * sm0.hc * sm0.dmdh, 0.2 * sc4[2], 0.25 * sc4[3]], // prior uncertainties sized to the sensitivity of each state
    F4 = (x, u) => { const s = ros2Step((t, y) => sm0.alg(y, u[0], undefined, undefined, u[1]).d, 0, x, hK); const J = s.M.map((r, i) => r.map((m, j) => ((i === j ? 1 : 0) - m) / (g4 * hK))), Phi = inv(madd(eye(4), J.map((r) => r.map((x2) => x2 * hK)), -1)); return { x: s.y.map((y, i) => clamp(y, 1e-6, 50 * sc4[i])), J: Phi }; };
  const ekfU = idx.map((i) => [pidRun.z[i], pidRun.pSep[i] * 1e5]), ekfY = idx.map((i) => [pidRun.pIn[i] + v.noiseBar * rnd.normal(), pidRun.pTop[i] + v.noiseBar * rnd.normal(), pidRun.pBase[i] + v.noiseBar * rnd.normal()]);
  const ekfRes = ekf({ F: F4, h: (x, u) => { const a = sm0.alg(x, u[0], undefined, undefined, u[1]); return [a.Pp / 1e5, a.Prt / 1e5, a.Prb / 1e5]; }, x0: eqT.y.map((x, i) => x * (i === 3 ? 1.25 : 1)), P0: eye(4).map((r, i) => r.map((x) => x * sig4[i] ** 2)), Q: eye(4).map((r, i) => r.map((x) => x * (0.03 * sig4[i]) ** 2)), R: [[v.noiseBar ** 2 + 1e-3, 0, 0], [0, v.noiseBar ** 2 + 1e-3, 0], [0, 0, v.noiseBar ** 2 + 1e-3]],
    u: ekfU, y: ekfY });
  const mTrue = idx.map((i) => pidRun.mLr[i]), mEst = ekfRes.x.map((x) => x[3]), half = Math.floor(idx.length / 4), ekfRmse = Math.sqrt(mean(mEst.slice(half).map((x, i) => (x - mTrue[half + i]) ** 2))), ekfRel = ekfRmse / Math.max(mean(mTrue), 1e-9);
  const tf = ss2tf(lin.A, lin.B, lin.C[0].map((x) => x / 1e5)), marg = loopAnalysis({ k: foUse.K, lags: [foUse.tau], delay: foUse.theta }, { kc: rules.rows[1].kc, ti: rules.rows[1].ti }, { n: 900 }), clPoles = pole.poles2(sel.kc, sel.ti);
  const scS = lin.ys.map((x) => Math.abs(x) || 1), dS = c2d(lin.A.map((r, i) => r.map((x, j) => (x * scS[j]) / scS[i])), lin.B.map((r, i) => [r[0] / scS[i]]), ctlCfg.tsMpc), cS = lin.C[0].map((x, j) => (x * scS[j]) / 1e5);
  let lqGain = null; try { lqGain = lqFinite(dS.Ad, dS.Bd, cS.map((a) => cS.map((b) => a * b)), [[v.mpcR]], 40).K0[0]; } catch { lqGain = null; }
  // unscented Kalman filter on the same model and measurements (first part of the record), beside the extended filter
  const nUk = Math.min(idx.length, 160), ukfRes = ukf({ F: (x, u) => ros2Step((t, y) => sm0.alg(y, u[0], undefined, undefined, u[1]).d, 0, x, hK).y.map((y, i) => clamp(y, 1e-6, 50 * sc4[i])), h: (x, u) => { const a = sm0.alg(x, u[0], undefined, undefined, u[1]); return [a.Pp / 1e5, a.Prt / 1e5, a.Prb / 1e5]; }, x0: eqT.y.map((x, i) => x * (i === 3 ? 1.25 : 1)), P0: eye(4).map((r, i) => r.map((x) => x * sig4[i] ** 2)), Q: eye(4).map((r, i) => r.map((x) => x * (0.03 * sig4[i]) ** 2)), R: [[v.noiseBar ** 2 + 1e-3, 0, 0], [0, v.noiseBar ** 2 + 1e-3, 0], [0, 0, v.noiseBar ** 2 + 1e-3]], u: ekfU.slice(0, nUk), y: ekfY.slice(0, nUk) });
  const h4 = Math.floor(nUk / 4), relOf = (est) => { const e = est.slice(h4, nUk).map((x, i) => (x[3] - mTrue[h4 + i]) ** 2); return e.length && e.every(Number.isFinite) ? Math.sqrt(mean(e)) / Math.max(mean(mTrue.slice(h4, nUk)), 1e-9) : null; }, ukfRel = relOf(ukfRes.x), ekfRelWin = relOf(ekfRes.x);
  // sensitivity peaks of the applied loop on the linearised riser model with the loop dead time (meaningful while the closed loop is stable)
  const loopLin = loopAnalysis({ k: 1, num: tf.num, den: tf.den, delay: loopDead + 0.5 * dtCtl, tRef: foUse.tau }, { kc: sel.kc, ti: sel.ti, td: sel.td }, { n: 900 });
  // nonlinear model-predictive control on the nonlinear riser model (single shooting, SQP), full state feedback
  let nm = null;
  if (v.nmpc && !manual) { try { nm = nmpc(sm0, { z0: zTarget, tEnd: HOUR, np: 5, maxIter: 2, ts: Math.max(ctlCfg.tsMpc, 120), dSp: v.spStep, tStep: 0.3 * HOUR, r: v.mpcR, zMin: Math.max(v.zMinPct / 100, sm0.zFloor), zMax: v.zMaxPct / 100, duMax: (v.rateLimit / 100) * ctlCfg.tsMpc }); } catch { nm = null; } }
  const slugSuppressed = !unstable || (pidRun.suppressed && !manual), chokeOut = slugSuppressed ? 100 * zTarget : 100 * (zCrit === null ? zTarget : 0.9 * zCrit), slugAmpOpen = unstable ? pidRun.ampOpen : 0;
  await tick();

  // ---------- B. restart and ramp-up ----------
  progress(0.55, 'Restart and ramp-up');
  const gelLen = sum(tAtShut.map((T) => (T < v.pourPoint ? st.ds : 0))), dpGel = (4 * v.yieldStress * gelLen) / id / 1e5, restartPressure = (v.sepP + so.headUphill + dpGel) * 1.05;
  const capSteel = STEEL.rho * STEEL.cp * Math.PI * ((id / 2 + wt) ** 2 - (id / 2) ** 2) * v.cMult, capCoat = sum(gridOf(hSea).C) - capSteel, cWall = Math.max(capSteel + 0.5 * Math.max(capCoat, 0), 1e3), gIn = hFlow * Math.PI * id, rTot = 1 / (wall.uFlow * Math.PI * id), gOut = 1 / Math.max(rTot - 1 / gIn, 1e-6);
  const cFl = st.holdup.map((H, i) => { const p = fm.at(st.P[i], st.T[i]); return A * (H * p.rhoL * p.cpL + (1 - H) * p.rhoG * p.cpG); });
  const tWarmEnd = Math.max(3 * rampS, 30 * HOUR), nWarm = 600, refHyd = st.P.map((p) => fm.hydrateT(p) + v.hydMargin);
  const warm = warmUp({ s: st.s, ds: st.ds, tAmb: st.s.map((_, i) => tAmbOf(i)), T0: tAtShut, mdot: (t) => mdot * (q0 + (1 - q0) * Math.min(t / rampS, 1)), cp: cpMix, cFluid: cFl, cWall: st.s.map(() => cWall), gIn, gOut, tIn: line.tIn, dt: tWarmEnd / nWarm, nSteps: nWarm, ref: refHyd });
  const tW = warm.t.map((t) => t / HOUR), tOutEnd = warm.Tout[nWarm], iSteady = warm.Tout.findIndex((T, k) => warm.Tout.slice(k).every((x) => Math.abs(x - tOutEnd) <= 1)), restartTime = tW[Math.max(iSteady, 0)];
  const iSafe = warm.Tmin.findIndex((_, k) => warm.Tmin.slice(k).every((x) => x >= 0)), tSafe = iSafe >= 0 ? tW[iSafe] : null;
  let rampReq = v.rampHours;
  if (ramp.vMax > surgeAllow) { let lo = v.rampHours, hi = v.rampHours * 2 + 1; for (let k = 0; k < 8 && surgeFor(hi).vMax > surgeAllow; k++) { lo = hi; hi *= 2; } if (surgeFor(hi).vMax <= surgeAllow) { for (let k = 0; k < 14; k++) { const m = 0.5 * (lo + hi); if (surgeFor(m).vMax > surgeAllow) lo = m; else hi = m; } rampReq = hi; } else rampReq = null; }
  else { let lo = 0.05, hi = v.rampHours; if (surgeFor(lo).vMax <= surgeAllow) rampReq = lo; else { for (let k = 0; k < 14; k++) { const m = 0.5 * (lo + hi); if (surgeFor(m).vMax > surgeAllow) lo = m; else hi = m; } rampReq = hi; } }
  // dynamic optimisation of a three-segment ramp: least deferred production with the surge inside the allowance
  const rampObj = (x) => { const T = Math.max(x[0], 0.05) * HOUR, f1 = clamp(x[1], 0, 1), f2 = clamp(Math.max(x[2], f1), 0, 1), r = rampSurge({ nodes: [[0, q0 * rate], [T / 3, rate * (q0 + (1 - q0) * f1)], [(2 * T) / 3, rate * (q0 + (1 - q0) * f2)], [T, rate]], inv: invOf, qLiq: qLiqOf, qDrain, inv0: liqInvSteady, dt: 240, tEnd: 60 * HOUR }); return r.deferred / HOUR + 50 * Math.max(0, r.vMax / Math.max(surgeAllow, 1e-6) - 1); };
  const rampT0 = rampReq ?? v.rampHours * 4, rampOpt = nelderMead(rampObj, [rampT0, 1 / 3, 2 / 3], { lo: [0.05, 0, 0], hi: [Math.max(4 * rampT0, 2), 1, 1], maxIter: 60, tol: 1e-5 }), rampLin = rampObj([rampT0, 1 / 3, 2 / 3]);
  await tick();

  // ---------- C. depressurisation ----------
  progress(0.63, 'Blowdown');
  const tB = clamp(v.tBlowStart, 0, v.tHorizon), TgB = sum(cool.map((c, i) => atTime(c, tB) * vg[i])) / Math.max(sum(vg), 1e-9), pB = interp1(tH, pT, tB), gB = fm.at(pB, TgB), tAmbMin = Math.min(...st.s.map((_, i) => tAmbOf(i)));
  // gas pockets behind liquid seals: a pocket at a crest keeps the vent pressure plus the liquid legs between it and the nearest open end
  const trapTop = (c) => sum(so.seals.up.slice(c)), trapIn = (c) => sum(so.seals.down.slice(0, c)), crests = so.seals.up.map((_, c) => c).concat(so.seals.up.length), bothEnds = v.bdRoute === 'both';
  const headTop = Math.max(...crests.map(trapTop)), headBoth = Math.max(...crests.map((c) => Math.min(trapTop(c), trapIn(c)))), headBd = bothEnds ? headBoth : headTop;
  const pSafeRaw = fm.hydrateP(tAmbMin - v.hydMargin), pSafe = pSafeRaw === null ? 700 : pSafeRaw, pEndAuto = Math.max(v.pBack * 1.1, pSafe - headBd), pEndB = v.pBlowEnd > 0 ? Math.max(v.pBlowEnd, v.pBack * 1.02) : pEndAuto;
  const zGrid = linspace(Math.log(1), Math.log(Math.max(pB, 2) * 1.05), 14), zVals = zGrid.map((lp) => fm.at(Math.exp(lp), TgB).zG), Zf = (P) => interp1(zGrid, zVals, Math.log(Math.max(P / 1e5, 1)));
  const mLiqLine = so.liquidVol * volScale * so.rhoL, steelMass = STEEL.rho * Math.PI * ((id / 2 + wt) ** 2 - (id / 2) ** 2) * L, hBar = clamp(so.liquidVol / (A * L), 0, 0.95), mHcLine = mLiqLine * (1 - gB.wcut) + so.gasMass * volScale;
  const flashRate = (P, T) => { const pb = P / 1e5; if (pb < 1.2) return 0; return (mHcLine * Math.max(0, fm.at(pb * 0.97, T - KEL).wG - fm.at(pb, T - KEL).wG)) / (0.03 * P); };
  const twoPhase = v.bdMode === 'hem' || v.bdMode === 'hrm', relaxing = v.bdMode === 'hrm';
  const bdCfg = { V: Math.max(so.gasVol * volScale, 1e-3), P0: pB * 1e5, T0: TgB + KEL, pBack: v.pBack * 1e5, area: (Math.PI * ((v.orificeMm / 1000) ** 2 + (bothEnds ? (v.orifice2Mm / 1000) ** 2 : 0))) / 4, cd: v.cdBlow, relax: relaxing ? (v.bdRelax > 0 ? v.bdRelax : 'dz') : 0, liqVol: so.liquidVol * volScale, k: v.kGas, mw: gB.mwG / 1000, Z: Zf, mode: 'wall', wallC: steelMass * STEEL.cp + mLiqLine * gB.cpL, wallUA: v.hGasWall * Math.PI * id * L * (1 - hBar), extUA: wall.uFlow * Math.PI * id * L, tAmb: mean(st.s.map((_, i) => tAmbOf(i))) + KEL, Tw0: TgB + KEL, pEnd: pEndB * 1e5, pMark: Math.max(pSafe - headBd, v.pBack * 1.02) * 1e5, n: ntBlow, jt: Math.max(gB.jtG, 0), flash: v.bdFlash || relaxing ? flashRate : null,
    hem: twoPhase ? { mLiq: mLiqLine, rhoL: so.rhoL, frac: clamp(v.bdLiquidFrac / 100, 0, 0.9) } : null };
  need(pB > v.pBack * 1.02, `The line pressure when the blowdown valve opens (${pB.toFixed(1)} bara) is not above the flare back-pressure (${v.pBack} bara).`);
  const bd = blowdown(bdCfg), bdReached = bd.tEnd !== null, blowdownTime = (bd.tEnd ?? bd.t[bd.t.length - 1]) / HOUR, bdMinT = bd.minT - KEL, bdMinTw = bd.minTw - KEL, bdMinTd = bd.minTdown - KEL, bdEndP = bd.pFinal / 1e5, seabedPAfter = bdEndP + headBd;
  const safeBd = pSafe - headBd > v.pBack * 1.02 && (!bdReached || seabedPAfter <= pSafe * 1.001 || v.pBlowEnd <= 0), safeTopOnly = pSafe - headTop > v.pBack * 1.02, peakStd = (bd.peak / (gB.mwG / 1000)) * VM_STD * DAY / 1e6; // million Sm³/d
  await tick();

  // ---------- D. pigging ----------
  progress(0.7, 'Pigging');
  const pigVsRate = scan.res.map((r, k) => { const p = pigOf(r); return { q: scan.q[k], v: p.vMean, transit: p.transit === null ? null : p.transit / HOUR, surge: p.surge }; }).filter((p) => p.transit !== null);
  const pigTransit = pig.transit === null ? null : pig.transit / HOUR, pigSurge = pig.surge * volScale, pigRuns = pigIv > 0 ? 365 / pigIv : 0;
  const waxLen = sum(rPig.T.slice(1).map((T, i) => (T - (wall.uFlow * (T - rPig.tAmb[i])) / hFlow < v.wat ? rPig.ds : 0))), waxVol = Math.PI * id * (v.waxThk / 1000) * waxLen;
  const pInPigMax = rPig.pIn + pig.dpExtra / 1e5, waxSeries = { t: [0], v: [0] };
  if (waxVol > 0 && pigIv > 0) { let w = 0; for (let k = 1; k <= 8; k++) { w += waxVol; waxSeries.t.push(k * pigIv); waxSeries.v.push(w); w *= 1 - clamp(v.pigEff, 0, 100) / 100; waxSeries.t.push(k * pigIv); waxSeries.v.push(w); } }
  const waxMax = Math.max(...waxSeries.v);

  // ---------- E. chemical injection ----------
  progress(0.75, 'Chemical injection');
  const S = fm.aq.S, tMinShut = Math.min(...tAtShut, ...comps.map((c) => interp1(tH, c.series, Math.min(v.tShut, v.tHorizon)))), dTshut = fm.hydrateT0(pAtShut) - tMinShut + v.inhMargin, dTsteady = Math.max(...st.P.map((p, i) => fm.hydrateT0(p) - st.T[i])) + v.inhMargin;
  const dTgov = v.dosingBasis === 'steady' ? dTsteady : v.dosingBasis === 'max' ? Math.max(dTshut, dTsteady) : dTshut, mWater = fm.rates.mW * rate, pOut = fm.at(v.sepP, tArr0, rate);
  const doseArgs = { inh: inhId, S, mWater, lean: v.leanWt, P: v.sepP, T: tArr0, qGasStd: fm.rates.qGasStd * rate, mOil: pOut.mO, rhoOil: isNum(pOut.rhoO) ? pOut.rhoO : pOut.rhoL, eff: clamp(v.injEff, 5, 100) / 100 };
  const doseGov = inhibitorDose({ ...doseArgs, dT: dTgov }), doseSteady = inhibitorDose({ ...doseArgs, dT: dTsteady }), doseWt = Math.max(doseGov.wt, v.inhRequired), steadyNeeded = dTsteady - v.inhMargin > hydrateDepression({ S, inhWt: 0, inh }) && (v.dosingBasis !== 'shutdown');
  const contRate = steadyNeeded ? (v.inhRequired > doseSteady.wt ? inhibitorDose({ ...doseArgs, dT: hydrateDepression({ S, inhWt: v.inhRequired, inh }) }).qInject : doseSteady.qInject) : 0;
  // batch treatment at shutdown: bring the settled water to the dose, plus the special components displaced 1.5 times
  const wDose = doseWt / 100, lean = clamp(v.leanWt, 30, 100) / 100, waterLine = sum(so.holdup.map((H, i) => H * fm.at(so.pSettle, st.T[i], rate).wcut)) * st.ds * A * volScale * 1020, batchVol = wDose > 0 ? (wDose * waterLine) / Math.max(1 - wDose / lean, 0.02) / lean / (lean * inh.rho + (1 - lean) * 1000) + v.bullheadVol : 0;
  const tBullhead = v.pumpMax > 0 ? batchVol / v.pumpMax : 0, restartInj = doseSteady.wt > 0 || dTshut > 0 ? inhibitorDose({ ...doseArgs, dT: dTshut, mWater: mWater * q0 }).qInject / 24 * (tSafe ?? restartTime) : 0;
  // inhibitor front along the flowing line (liquid velocity), analytic check with the mean velocity
  const uLiq = pic.x.map((_, i) => clamp(pic.vsl[i] / Math.max(pic.holdup[i], 0.02), 0.02, 50)), uMean = L / sum(st.s.map((_, i) => st.ds / clamp(st.vsl[i] / Math.max(st.holdup[i], 0.02), 0.02, 50)));
  const fDarcy = frictionFactor((pm.rhoL * uMean * id) / pm.muL, line.roughness / id), disp = Math.max(10.1 * (id / 2) * uMean * Math.sqrt(fDarcy / 8) * v.dispMult, 1e-4);
  const sInj = clamp(interp1(pic.x, st.S, clamp(v.injX, pic.x[0], pic.x[pic.x.length - 1])), 0, 0.95 * L), xF = [sInj, ...st.S.filter((s) => s > sInj + 1)], uF = xF.map((s) => interp1(st.S, uLiq, s));
  const front = inhibitorFront({ x: xF, u: uF, D: disp, n: Math.round(clamp(v.nxInh, 20, 600)), tEnd: 2.2 * ((L - sInj) / uMean) }), tProtect = (front.tProtect ?? front.transit * 2.2) / HOUR;
  // injection setting against requirement and pump capacity
  const contUsed = v.injRate > 0 ? v.injRate : contRate, mInj = (contUsed / DAY) * (lean * inh.rho + (1 - lean) * 1000) * lean * doseArgs.eff, doseAchieved = contUsed > 0 ? (100 * mInj) / (mWater + mInj / lean) : 0, doseShort = steadyNeeded ? Math.max(0, Math.max(doseSteady.wt, v.inhRequired) - doseAchieved) : 0, pumpUtil = contUsed / 24 / Math.max(v.pumpMax, 1e-9);
  const megLoop = inh.rho > 1000 ? { rich: doseSteady.rich, inventory: ((doseSteady.rich / 1050) * (liqInvSteady / Math.max(at('qLiq', rate) / 3600, 1e-6)) + contUsed * v.megStorageDays), duty: (mWater * (2.257e6 + 4186 * 90) + doseSteady.mAq * 2400 * 90) / 1000 } : null;
  const subGov = dTgov - v.inhMargin, wcStd = fm.rates.wc, ldhi = subGov <= 0 ? 'No hydrate driving force at the governing condition; no low-dosage inhibitor needed.' : `${subGov <= 10 ? 'A kinetic inhibitor (KHI) is a candidate: sub-cooling ' + subGov.toFixed(1) + ' °C ≤ 10 °C, provided the hold time stays below the tested induction time.' : 'Kinetic inhibitors are not suitable: sub-cooling ' + subGov.toFixed(1) + ' °C exceeds about 10 °C.'} ${wcStd <= 0.5 ? 'An anti-agglomerant (AA) is a candidate: water cut ' + (100 * wcStd).toFixed(0) + ' % ≤ 50 %, needs a liquid hydrocarbon phase and restart testing.' : 'Anti-agglomerants are doubtful: water cut ' + (100 * wcStd).toFixed(0) + ' % exceeds 50 %.'}`;

  // ---------- F. heating ----------
  const tHold = fm.hydrateT(so.pSettle) + v.hydMargin, heatLen = clamp(v.heatLengthPct, 0, 100) / 100, heatW = sum(st.s.map((_, i) => (st.s[i] <= heatLen * L ? Math.max(0, tHold - tAmbOf(i)) / gridOf(buried ? 1e4 : st.z[i] >= 0 ? hAir : hSea).Rtot * st.ds : 0))), heatingPower = heatW / 1000 / Math.max(v.heatEff / 100, 0.05);
  const hot = warmUp({ s: st.s, ds: st.ds, tAmb: st.s.map((_, i) => tAmbOf(i)), T0: st.s.map((_, i) => tAmbOf(i)), mdot: () => v.hotOilRate, cp: 2000, cFluid: st.s.map(() => A * 850 * 2000), cWall: st.s.map(() => cWall), gIn: 150 * Math.PI * id, gOut, tIn: v.hotOilT, dt: 240, nSteps: 450 }), iHot = hot.Tout.findIndex((T) => T >= tHold), hotOilTime = iHot >= 0 ? hot.t[iHot] / HOUR : null;
  // hot-oil circulation started while the line is still warm: steady return temperature and heater duty of the circulation
  const uL = (1 / (1 / gIn + 1 / gOut)) * L, tAmbMean = mean(st.s.map((_, i) => tAmbOf(i))), hotReturn = tAmbMean + (v.hotOilT - tAmbMean) * Math.exp(-uL / (v.hotOilRate * 2000)), hotOilDuty = (v.hotOilRate * 2000 * (v.hotOilT - hotReturn)) / 1000, hotOilHolds = hotReturn >= tHold;
  // cold start-up of the whole line from ambient (commissioning / first start) with the same ramp
  const cold = warmUp({ s: st.s, ds: st.ds, tAmb: st.s.map((_, i) => tAmbOf(i)), T0: st.s.map((_, i) => tAmbOf(i)), mdot: (t) => mdot * (q0 + (1 - q0) * Math.min(t / rampS, 1)), cp: cpMix, cFluid: cFl, cWall: st.s.map(() => cWall), gIn, gOut, tIn: line.tIn, dt: tWarmEnd / 300, nSteps: 300 });
  await tick();

  // ---------- I. envelope ----------
  progress(0.8, 'Operating envelope and optimisation');
  const inhibitedSteady = v.dosingBasis !== 'shutdown' && doseWt > 0, waxByPig = !!v.waxByPigging && pigIv > 0, cons = [
    { key: 'margin', name: 'hydrate margin', type: 'min', limit: inhibitedSteady ? -1e3 : v.hydMargin, unit: '°C' }, { key: 'tArr', name: 'arrival temperature above WAT', type: 'min', limit: waxByPig ? -1e3 : v.wat + v.watMargin, unit: '°C' },
    { key: 'eros', name: 'erosional velocity', type: 'max', limit: 1, unit: '–' }, { key: 'pReq', name: 'inlet pressure (incl. slug-stabilising choke Δp)', type: 'max', limit: v.pAvail, unit: 'bara' }, { key: 'qLiq', name: 'separator liquid capacity', type: 'max', limit: v.qDrainM3h, unit: 'm³/h' }];
  const turndown = num(v.turndown, 0); if (turndown > 0) { scan.qSelf = scan.q.slice(); cons.push({ key: 'qSelf', name: 'minimum stable rate (flow suite)', type: 'min', limit: turndown, unit: '–' }); }
  const env = operatingEnvelope(scan, cons, clamp(rate, qLo, qHi));
  if (waxByPig) env.text.push(`wax managed by pigging every ${pigIv} d (arrival may fall below the WAT)`);
  if (inhibitedSteady) env.text.push('hydrate margin provided by continuous inhibition');
  if (v.severeSlugging && !v.slugControl) env.text.push('severe slugging reported by the flow suite and no active slug control selected');

  // ---------- J. surrogate ----------
  const sur = {}; for (const k of ['pIn', 'tArr', 'margin', 'eros', 'qLiq', 'pReq', 'inv']) sur[k] = responseSurface(scan.q, scan[k], 3);

  // ---------- I. optimisation ----------
  const bblDay = fm.rates.qOilStd * 6.2898, mW1 = fm.rates.mW, depOf = (w) => hydrateDepression({ S, inhWt: clamp(w, 0, 90), inh }) - hydrateDepression({ S, inhWt: 0, inh });
  const chemCost = (q, w) => { const ww = clamp(w, 0, 90) / 100, vol = ((ww * mW1 * q) / Math.max(1 - ww / lean, 0.02) / lean / (lean * inh.rho + (1 - lean) * 1000)) * DAY; return megLoop ? vol * v.inhPrice * (v.megLossPct / 100) + (ww > 0 ? ((mW1 * q * 2.6e6) / 1000) * 24 * v.elecPrice : 0) : vol * v.inhPrice; }; // $/d; a glycol loop pays make-up and regeneration energy
  const heatCost = heatingPower * 24 * v.elecPrice, pigCostDay = (pigRuns * v.pigCost) / 365, needPig = (q) => sur.tArr.predict(q) < v.wat + v.watMargin;
  const gOf = (q, w, heat, pigs) => [sur.margin.predict(q) + depOf(w) + (heat ? 100 : 0) - v.hydMargin, heat || pigs || waxByPig ? 1 : sur.tArr.predict(q) - v.wat - v.watMargin, 1 - sur.eros.predict(q), (v.pAvail - sur.pReq.predict(q)) / 10, (v.qDrainM3h - sur.qLiq.predict(q)) / Math.max(v.qDrainM3h, 1)];
  const profit = (q, w, heat, pigs) => q * bblDay * v.oilPrice - chemCost(q, w) - (heat ? heatCost : 0) - (pigs ? pigCostDay : 0), scaleP = Math.max(bblDay * v.oilPrice, 1);
  const pen = (heat, pigs) => (x) => { const g = gOf(x[0], x[1], heat, pigs); return -profit(x[0], x[1], heat, pigs) / scaleP + 20 * sum(g.map((c) => Math.max(0, -c) ** 2)) + 2 * sum(g.map((c) => Math.max(0, -c))); };
  const lo2 = [qLo, 0], hi2 = [qHi, 60], combos = [[0, 0], [1, 0], [0, 1], [1, 1]].map(([h, p]) => { const r = nelderMead(pen(h, p), [clamp(rate, qLo, qHi), 10], { lo: lo2, hi: hi2, maxIter: 160, tol: 1e-9 }), g = gOf(r.x[0], r.x[1], h, p); return { heat: h, pigs: p, x: r.x, f: r.f, feasible: g.every((c) => c > -1e-3), profit: profit(r.x[0], r.x[1], h, p), evals: r.evals }; });
  const bestC = combos.filter((c) => c.feasible).sort((a, b) => b.profit - a.profit)[0] || combos.slice().sort((a, b) => a.f - b.f)[0], obj = pen(bestC.heat, bestC.pigs), methods = [{ name: 'Penalty + Nelder–Mead', x: bestC.x, f: obj(bestC.x), evals: bestC.evals, history: null }];
  const ga = geneticAlgorithm(obj, lo2, hi2, { pop: 20, gens: 25, seed: 3 }), pso = particleSwarm(obj, lo2, hi2, { n: 16, iters: 30, seed: 5 }), bo = bayesOpt(obj, lo2, hi2, { n0: 6, iters: 12, seed: 9, cand: 150 });
  methods.push({ name: 'Genetic algorithm', ...ga }, { name: 'Particle swarm', ...pso }, { name: 'Bayesian optimisation (GP + EI)', ...bo });
  try { // interior point from a strictly feasible start (inside the window, generous dose)
    const qs0 = env.feasible ? 0.5 * (env.qMin + env.qMax) : clamp(rate, qLo, qHi), gI = [0, 1, 2, 3, 4].map((k) => (x) => gOf(x[0], x[1], bestC.heat, bestC.pigs)[k] + 1e-9), start = [qs0, 45].map((x, i) => clamp(x, lo2[i] + 1e-3, hi2[i] - 1e-3));
    const ip = interiorPoint((x) => -profit(x[0], x[1], bestC.heat, bestC.pigs) / scaleP, gI, start, { lo: lo2, hi: hi2, mu: 0.1 }); methods.push({ name: 'Interior point (log barrier)', x: ip.x, f: obj(ip.x), evals: ip.iterations, history: ip.history.map((f) => f) });
  } catch { methods.push({ name: 'Interior point (log barrier)', x: null, f: null, evals: 0, history: null, note: 'no strictly feasible start' }); }
  try { // sequential quadratic programming on the same constrained problem (constraints handled exactly, no penalty)
    const sq = sqp((x) => -profit(x[0], x[1], bestC.heat, bestC.pigs) / scaleP, [clamp(rate, qLo, qHi), 10], { ineq: [0, 1, 2, 3, 4].map((k) => (x) => gOf(x[0], x[1], bestC.heat, bestC.pigs)[k]), lo: lo2, hi: hi2, maxIter: 40, tol: 1e-7, h: 1e-4 });
    if (sq.x.every(Number.isFinite)) methods.push({ name: 'Sequential quadratic programming (BFGS, active-set QP)', x: sq.x, f: obj(sq.x), evals: sq.iterations, history: sq.history.length > 1 ? sq.history : null });
  } catch { /* the other methods remain */ }
  // linear programme on the constraints linearised at the nonlinear optimum (variables: rate above the scan floor, dose) and its mixed-integer extension
  const xb = bestC.x, e1 = 1e-3, dgq = gOf(xb[0] + e1, xb[1], 0, 0).map((g, k) => (g - gOf(xb[0] - e1, xb[1], 0, 0)[k]) / (2 * e1)), dgw = gOf(xb[0], xb[1] + 0.05, 0, 0).map((g, k) => (g - gOf(xb[0], xb[1] - 0.05, 0, 0)[k]) / 0.1), g0 = gOf(xb[0], xb[1], 0, 0);
  const cq = bblDay * v.oilPrice - (chemCost(xb[0] + e1, xb[1]) - chemCost(xb[0] - e1, xb[1])) / (2 * e1), cw = (chemCost(xb[0], xb[1] + 0.05) - chemCost(xb[0], Math.max(xb[1] - 0.05, 0))) / (xb[1] >= 0.05 ? 0.1 : 0.05), bigM = 1e3;
  const rowsLP = [], rhsLP = []; for (let k = 0; k < 5; k++) { rowsLP.push([-dgq[k], -dgw[k], k === 0 || k === 1 ? -bigM : 0, k === 1 ? -bigM : 0]); rhsLP.push(g0[k] - dgq[k] * (xb[0] - qLo) - dgw[k] * xb[1]); }
  rowsLP.push([1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]); rhsLP.push(qHi - qLo, 60, 1, 1);
  const cLP = [cq / scaleP, -cw / scaleP, -heatCost / scaleP, -pigCostDay / scaleP], lp = simplex(cLP.slice(0, 2), rowsLP.filter((_, k) => k < 7).map((r) => r.slice(0, 2)), rhsLP.slice(0, 7)), milp = branchBound(cLP, rowsLP, rhsLP, [2, 3]);
  if (lp.status === 'optimal') methods.push({ name: 'Linear programme (simplex, linearised)', x: [qLo + lp.x[0], lp.x[1]], f: obj([qLo + lp.x[0], lp.x[1]]), evals: lp.iterations, history: null });
  const optBest = methods.filter((m) => m.x && isNum(m.f)).sort((a, b) => a.f - b.f)[0], qOpt = optBest.x[0], wOpt = optBest.x[1], profitOpt = profit(qOpt, wOpt, bestC.heat, bestC.pigs);
  // shutdown strategy: cheapest feasible preservation for the planned duration
  const gasValue = 0.2, strategies = [
    { id: 'none', phrase: 'restarting without preservation', name: 'No action (restart inside the cooldown time)', feasible: v.tShut <= cooldownTime || neverCools, cost: 0, lead: 0 },
    { id: 'inhibit', phrase: `${inhId} bullheading`, name: `${inhId} bullheading`, feasible: doseGov.attainable && tBullhead + v.tDecision <= cooldownTime, cost: (batchVol + restartInj) * v.inhPrice, lead: tBullhead },
    { id: 'blowdown', phrase: 'depressurising', name: bothEnds ? 'Depressurise from both ends (topside valve and inlet service line)' : 'Depressurise through the topside blowdown valve', feasible: safeBd && bdReached && blowdownTime + v.tDecision <= cooldownTime && Math.min(bdMinTw, bdMinTd) >= v.tMinDesign, cost: bd.discharged * gasValue + (blowdownTime / 24) * bblDay * rate * v.oilPrice * 0.25, lead: blowdownTime },
    { id: 'heat', phrase: 'switching on the heating', name: 'Electrical heating' + (v.heatInstalled ? '' : ' (not installed)'), feasible: !!v.heatInstalled && heatingPower > 0 && heatingPower <= v.heatMaxKw, cost: heatingPower * v.tShut * v.elecPrice, lead: 0.5 },
    { id: 'hotoil', phrase: 'circulating hot oil', name: 'Hot-oil circulation' + (v.hotOilInstalled ? '' : ' (no circulation loop)'), feasible: !!v.hotOilInstalled && hotOilHolds, cost: hotOilDuty * v.tShut * v.elecPrice, lead: 1 }];
  if (neverCools) strategies[0].feasible = true;
  if (trig('hydrate') && !neverCools && v.tShut > 0) { strategies[0].feasible = false; strategies[0].name += ' — not permitted: hydrate accumulation trigger'; }
  // lowest cost among the feasible strategies that leave the required no-touch time; if none does, lowest cost among the feasible ones
  const meets = (s) => s.feasible && (neverCools || s.id === 'none' || cooldownTime - v.tDecision - s.lead >= v.noTouchMin), byCost = (a, b) => a.cost - b.cost;
  const stratBest = strategies.filter(meets).sort(byCost)[0] || strategies.filter((s) => s.feasible).sort(byCost)[0] || null, stratSel = v.preserve === 'auto' ? stratBest || strategies[1] : strategies.find((s) => s.id === v.preserve) || strategies[0];
  await tick();

  // ---------- no-touch time, maximum shutdown ----------
  const noTouch = Math.max(0, cooldownTime - v.tDecision - stratSel.lead), gelLimit = (() => { // time at which the gelled-line restart pressure exceeds what is available
    if (!(v.yieldStress > 0)) return null;
    for (let k = 0; k < tH.length; k++) { const len = sum(cool.map((c) => (c.Tf[k] < v.pourPoint ? st.ds : 0))); if ((v.sepP + so.headUphill + (4 * v.yieldStress * len) / id / 1e5) * 1.05 > v.pAvail) return tH[k]; }
    return null; })();
  const preserved = stratSel.id !== 'none' && stratSel.feasible, maxShutdown = preserved ? (gelLimit ?? v.tHorizon) : Math.min(cooldownTime, gelLimit ?? Infinity);

  // ---------- availability ----------
  const rampLoss = ramp.deferred / HOUR, perPlanned = v.plannedHours + restartTime * 0 + rampLoss + (stratSel.id === 'blowdown' ? blowdownTime : 0), remediate = v.unplannedHours > cooldownTime && !neverCools ? Math.max(tBullhead, stratSel.id === 'blowdown' ? blowdownTime : 0) + v.remediationHours : 0;
  const perUnplanned = v.unplannedHours + rampLoss + remediate, pigLossH = pigTransit === null ? 0 : pigRuns * pigTransit * Math.max(0, 1 - qPig / Math.max(rate, 1e-9)), downH = v.plannedPerYear * perPlanned + v.unplannedPerYear * perUnplanned + pigLossH;
  const uptime = clamp(1 - downH / 8760, 0, 1), deferredVolume = (fm.rates.qOilStd * rate * downH) / 24, shutdownsYr = v.plannedPerYear + v.unplannedPerYear, blowdownsYr = stratSel.id === 'blowdown' ? shutdownsYr : v.unplannedPerYear * (remediate > 0 ? 1 : 0);
  const inhibitorRate = contUsed + (((stratSel.id === 'inhibit' ? batchVol : 0) + restartInj) * shutdownsYr) / 365, inhCostDay = megLoop ? inhibitorRate * v.inhPrice * (v.megLossPct / 100) + megLoop.duty * 24 * v.elecPrice * (contUsed > 0 ? 1 : 0) : inhibitorRate * v.inhPrice;

  // ---------- chemical inventory and ratio control of the injection ----------
  progress(0.86, 'Chemical inventory and rotating equipment');
  const rhoLean = lean * inh.rho + (1 - lean) * 1000, batchEvent = stratSel.id === 'inhibit' ? batchVol + restartInj : restartInj, tank0 = (v.tankVol * clamp(v.tankLevelPct, 0, 100)) / 100, invDays = 90;
  const nEv = Math.max(0, Math.round((shutdownsYr * invDays) / 365)), batches = Array.from({ length: nEv }, (_, k) => ({ t: ((k + 0.5) * invDays) / Math.max(nEv, 1), v: batchEvent })), reorder = Math.min(v.tankVol, v.resupplyLead * contUsed + (v.reorderBatches > 0 ? v.reorderBatches * batchEvent : 0));
  const stock = chemicalInventory({ V: v.tankVol, level0: tank0, use: contUsed, batches, reorder, lead: v.resupplyLead, delivery: v.resupplyVol, tEnd: invDays, dt: 0.25 }), autonomy = contUsed > 1e-9 ? tank0 / contUsed : null, treatments = batchEvent > 1e-9 ? tank0 / batchEvent : null;
  // ratio control during the restart: inhibitor in ratio to the water rate against a fixed injection sized for the full rate
  const wRest = dTshut > 0 ? inhibitorDose({ ...doseArgs, dT: dTshut }) : doseSteady, ratioKg = wRest.mTotal / Math.max(mWater, 1e-9), tRatio = Math.max((tSafe ?? restartTime) * HOUR, rampS, 600), rampUse = (rampReq ?? v.rampHours) * HOUR;
  const ratioRun = ratioControl({ wild: (t) => mWater * (q0 + (1 - q0) * Math.min(t / Math.max(rampUse, 1), 1)), ratio: ratioKg, tau: Math.max(v.tauPump, 1), kc: 0.5, ti: 3 * Math.max(v.tauPump, 1), uMax: (v.pumpMax / 3600) * rhoLean * lean, dt: Math.max(tRatio / 1500, 2), tEnd: tRatio });
  const ratioSaved = Math.max(0, ratioKg * mWater * tRatio - ratioRun.injected) / (rhoLean * lean), ratioShort = ratioRun.required > 0 ? ratioRun.shortfall / ratioRun.required : 0;

  // ---------- rotating equipment: liquid export pump and gas compressor ----------
  const qLiqArr = Math.max(at('qLiq', rate) / 3600, 1e-5), pumpQr = 1.15 * qLiqArr, pumpHr = (v.pumpDp * 1e5) / (pOut.rhoL * G), pumpTbl = (Array.isArray(v.pumpCurve) ? v.pumpCurve : []).map((r) => ({ q: num(r.q, 0) / 3600, h: num(r.h, 0) })).filter((r) => r.q >= 0 && r.h > 0).sort((a, b) => a.q - b.q);
  const pumpCv = pumpTbl.length >= 3 ? { q: pumpTbl.map((r) => r.q), h: pumpTbl.map((r) => r.h) } : null, pumpTd = clamp(v.pumpTurndownPct, 5, 100) / 100, pumpTrip = 80, pumpEnd = 110;
  const pump = pumpSim({ qr: pumpCv ? interp1(pumpCv.h.slice().reverse(), pumpCv.q.slice().reverse(), pumpHr) || pumpQr : pumpQr, hr: pumpHr, shutoff: v.pumpShutoff, eta: v.pumpEta / 100, rho: pOut.rhoL, J: v.pumpJ, rpm: v.pumpRpm, hStatic: (clamp(v.pumpStaticPct, 0, 95) / 100) * pumpHr, curve: pumpCv, qSp: (t) => (t < 30 ? qLiqArr : qLiqArr * pumpTd), tTrip: pumpTrip, tEnd: pumpEnd, dt: 0.02, n0: clamp(v.rotN0Pct, 20, 110) / 100, nMin: clamp(v.rotNminPct, 0, 100) / 100, nMax: clamp(v.rotNmaxPct, 100, 120) / 100 });
  const iP1 = pump.t.findIndex((t) => t >= 29), iP2 = pump.t.findIndex((t) => t >= pumpTrip - 1), pumpDuty = { n: pump.n[iP1], q: pump.q[iP1], p: pump.power[iP1] }, pumpLow = { n: pump.n[iP2], q: pump.q[iP2], p: pump.power[iP2] }, pumpMinFlow = (clamp(v.pumpMinFlowPct, 0, 100) / 100) * pump.q[0], pumpFlowMargin = pumpMinFlow > 0 ? 100 * (pumpLow.q / pumpMinFlow - 1) : 100;
  const cooled = v.coolerOutT > 0 && v.coolerOutT < tArr0, tSuc = cooled ? v.coolerOutT : tArr0, gS = fm.at(v.sepP, tSuc, rate), mGas = Math.max(pOut.mG, 1e-3), kPoly = (k1, eta) => (k1 - 1) / (k1 * eta);
  const compHeadOf = (T) => { const zrt = (gS.zG * R * (T + KEL)) / (gS.mwG / 1000), nmx = kPoly(v.kGas, v.compEta / 100); return (zrt / nmx) * (Math.max(v.compPd / v.sepP, 1.05) ** nmx - 1); };
  const compMd = 1.1 * mGas, rhoSuc = (v.sepP * 1e5 * (gS.mwG / 1000)) / (gS.zG * R * (tSuc + KEL)), compQd = compMd / rhoSuc, compHd = compHeadOf(tSuc), compTbl = (Array.isArray(v.compMap) ? v.compMap : []).map((r) => ({ q: num(r.q, 0), h: num(r.h, 0) * 1000 }));
  const compCurve = compressorCurve({ qd: compQd, hd: compHd, table: compTbl }), compTd = clamp(v.compTurndownPct, 5, 100) / 100, compTrip = 240, compFeed = (t) => mGas * (t < 60 ? 1 : t < 70 ? 1 - (1 - compTd) * ((t - 60) / 10) : compTd);
  const comp = compressorSim({ curve: compCurve, qd: compQd, hd: compHd, eta: v.compEta / 100, ps: v.sepP * 1e5, T1: tSuc + KEL, Z: gS.zG, k: v.kGas, mw: gS.mwG / 1000, J: v.compJ, rpm: v.compRpm, smCtl: v.compSmCtl, antiSurge: !!v.antiSurge, feed: compFeed, tTrip: compTrip, tEnd: 320, dt: 0.1, n0: clamp(v.rotN0Pct, 70, 110) / 100, nMin: clamp(v.rotNminPct, 0, 100) / 100, nMax: clamp(v.rotNmaxPct, 100, 120) / 100, pMax: (clamp(v.compPmaxPct, 50, 300) / 100) * ((compMd * compHd) / (v.compEta / 100)) });
  const iC1 = comp.t.findIndex((t) => t >= 58), iC2 = comp.t.findIndex((t) => t >= compTrip - 2), compDuty = { n: comp.n[iC1], sm: comp.sm[iC1], p: comp.power[iC1], r: comp.recycle[iC1] }, compLow = { n: comp.n[iC2], sm: comp.sm[iC2], p: comp.power[iC2], r: comp.recycle[iC2] };
  const compPowerNoCool = (mGas * compHeadOf(tArr0)) / (v.compEta / 100), compPowerDuty = (mGas * compHeadOf(tSuc)) / (v.compEta / 100), coolerDuty = cooled ? (mGas * gS.cpG * (tArr0 - tSuc)) / 1000 : 0, coolerHydT = fm.hydrateT(v.sepP), mapPt = typeof NET.mapPoint === 'function' && compCurve.map ? NET.mapPoint(compCurve.map, (mGas / rhoSuc), compHd * (compDuty.n ** 2 || 1)) : null;

  // ---------- H. alarms, interlocks, event sequence ----------
  progress(0.9, 'Operating logic');
  const pOperate = pIn0 + (eqT.Prt - sm0.p.Ps) / 1e5, steadyMargin = Math.min(...st.T.map((T, i) => T - fm.hydrateT(st.P[i]))), worstSurge = Math.max(ramp.vMax, pigSurge, v.slugSurge, slugSuppressed ? pidRun.surgeClosed : pidRun.surge), catcherLevel = clamp(clamp(v.levelSp, 10, 90) + (100 * worstSurge) / catVol, 0, 400);
  const almVals = { tArr: tArr0, pIn: pOperate, hydMargin: steadyMargin, catcherLevel, tBlowMin: Math.min(bdMinTw, bdMinTd), noTouch, erosion: at('eros', rate), restartP: restartPressure, pPig: pInPigMax, slugAmp: slugSuppressed ? pidRun.ampClosed * (unstable ? 1 : 0) : slugAmpOpen, cooldown: cooldownTime, doseShort, pumpUtil: 100 * pumpUtil, watMargin: tArr0 - v.wat, ...(v.plugTime > 0 ? { plugTime: v.plugTime } : {}),
    deposit: Math.max(depWax, depScale), hydFrac: 100 * v.hydFrac, blockage: v.blockagePct, surgeMargin: Math.max(comp.minSm, -100), pumpFlowMargin, compPower: (100 * compDuty.p) / Math.max(comp.design.power, 1e-9), ...(autonomy !== null ? { autonomy } : {}), ...(treatments !== null ? { treatments } : {}) };
  const rulesTbl = (Array.isArray(v.alarms) ? v.alarms : []).filter((r) => r && r.key && r.tag), raised = evaluateAlarms(almVals, rulesTbl), alarms = raised.map((a) => ({ tag: a.tag, level: a.level, msg: a.msg }));
  for (const l of env.limits) if (l.margin < 0) alarms.push({ tag: 'ENV-' + l.key.toUpperCase(), level: 'alarm', msg: `Operating point violates the ${l.name} limit (${txt(l.value)} against ${txt(l.limit)} ${l.unit}).` });
  for (const t of triggers) alarms.push({ tag: t.tag, level: 'alarm', msg: `${t.msg} → ${t.action}.` });
  if (comp.surgeEvents > 0) alarms.push({ tag: 'UA-160', level: 'alarm', msg: `Compressor surge in the turndown test (${comp.surgeEvents} event(s)): anti-surge control does not hold the margin.` });
  if (sensorFailed) alarms.push({ tag: 'PT-100', level: 'alarm', msg: `Inlet-pressure transmitter failed: slug controller forced to manual at ${(100 * zTarget).toFixed(0)} % opening.` });
  if (restartPressure > v.pAvail) alarms.push({ tag: 'RST-P', level: 'trip', msg: `Restart needs ${restartPressure.toFixed(0)} bara but only ${v.pAvail} bara is available.` });
  const trips = raised.filter((a) => a.level === 'trip'), safeState = trips.length ? stateMachine(['trip']).state : 'PRODUCING';
  // planned shutdown → preservation → restart sequence
  const tAct = Math.max(0, cooldownTime - stratSel.lead - v.tDecision), tRestart = Math.max(v.tShut, 0.1), seq = [{ t: 0, tag: 'XV-001', action: 'Close production choke and wing valves (planned shutdown)', ev: 'shutdown', set: { rate: 0 } }];
  if (inhibitedSteady === false && stratSel.id === 'inhibit' && tProtect > 0) seq.unshift({ t: 0, tag: 'P-201', action: `Inhibitor front already through the line: injection started ${tProtect.toFixed(1)} h before shut-in`, ev: null, set: {} });
  if (stratSel.id === 'inhibit' && tAct < tRestart) seq.push({ t: tAct, tag: 'P-201', action: `Start ${inh.name} bullheading at ${v.pumpMax} m³/h (${batchVol.toFixed(1)} m³)`, ev: 'inhibit', set: { inj: v.pumpMax } }, { t: Math.min(tAct + tBullhead, tRestart), tag: 'P-201', action: 'Stop bullheading: line inhibited', ev: null, set: { inj: 0 } });
  if (stratSel.id === 'blowdown' && tAct < tRestart) seq.push({ t: Math.min(tAct, tRestart), tag: 'BDV-301', action: `Open blowdown valve (${v.orificeMm} mm orifice) to ${pEndB.toFixed(1)} bara`, ev: 'blowdown', set: { bdv: 1 } }, { t: Math.min(tAct + blowdownTime, tRestart), tag: 'BDV-301', action: 'Close blowdown valve: line depressurised', ev: null, set: { bdv: 0 } });
  if (stratSel.id === 'heat') seq.push({ t: Math.min(tAct, tRestart), tag: 'DEH-401', action: `Switch on heating at ${heatingPower.toFixed(0)} kW`, ev: 'heat', set: { heat: 1 } });
  if (stratSel.id === 'hotoil') seq.push({ t: Math.min(tAct, tRestart), tag: 'P-601', action: `Start hot-oil circulation at ${v.hotOilRate} kg/s, ${v.hotOilT} °C (heater duty ${hotOilDuty.toFixed(0)} kW)`, ev: 'heat', set: { heat: 1 } });
  if (trig('hydrate')) seq.unshift({ t: 0, tag: 'XA-152', action: 'Hydrate accumulation trigger: inhibitor injection confirmed before shut-in', ev: null, set: {} });
  if (trig('blockage')) seq.unshift({ t: 0, tag: 'XA-153', action: 'Bore-restriction trigger: rate reduced, remediation planned', ev: null, set: {} });
  if (stock.tReorder !== null && stock.tReorder * 24 <= tRestart + 72) seq.push({ t: stock.tReorder * 24, tag: 'LAL-220', action: `Chemical storage at the re-order level: order ${v.resupplyVol} m³ (${v.resupplyLead} d lead time)`, ev: null, set: {} });
  seq.push({ t: tRestart, tag: 'XV-001', action: `Open wing valves; restart at ${v.qStartPct} % with ${inh.name} injection (needs ${restartPressure.toFixed(0)} bara)`, ev: 'restart', set: { rate: q0, heat: 0, inj: 0 } });
  if (tSafe !== null && tSafe > 0) seq.push({ t: tRestart + tSafe, tag: 'P-201', action: 'Stop restart inhibitor: whole line outside the hydrate region', ev: null, set: {} });
  seq.push({ t: tRestart + (rampReq ?? v.rampHours), tag: 'FIC-101', action: 'Ramp complete: full rate', ev: 'rampDone', set: { rate: 1 } }, { t: tRestart + Math.max(restartTime, rampReq ?? v.rampHours), tag: 'TI-102', action: 'Arrival temperature at steady state', ev: null, set: {} });
  const pigLaunchH = trig('pig') || pigInLine ? 0 : v.pigLaunch;
  if (pigLaunchH >= 0 && pigTransit !== null) { const tl = tRestart + Math.max(restartTime, rampReq ?? v.rampHours) + pigLaunchH; seq.push({ t: tl, tag: 'PL-501', action: pigInLine ? `Pig already in the line at ${(v.pigX0 / 1000).toFixed(1)} km: resume the run at ${v.pigRatePct} % rate` : `Launch pig at ${v.pigRatePct} % rate${trig('pig') ? ' (deposit trigger)' : ''}`, ev: 'pig', set: { rate: Math.min(1, qPig / Math.max(rate, 1e-9)) } }, { t: tl + pigTransit, tag: 'PR-502', action: `Receive pig; ${(pig.received * volScale).toFixed(0)} m³ of liquid ahead of it`, ev: 'pigReceived', set: { rate: 1 } }); }
  const seqS = seq.map((e) => ({ ...e, t: e.t * HOUR })), acc = { inj: 0, lost: 0 }, tSeqEnd = Math.max(...seqS.map((e) => e.t)) + HOUR;
  const sched = eventScheduler(seqS, { tEnd: tSeqEnd, dt: 600, state: { rate: 1, inj: 0, bdv: 0, heat: 0 }, onStep: (a, b, s) => { acc.inj += (s.inj * (b - a)) / HOUR; acc.lost += (1 - s.rate) * (b - a); } });
  const order = seq.map((e, i) => ({ e, i })).sort((a, b) => a.e.t - b.e.t || a.i - b.i).map((o) => o.e), fsm = stateMachine(order.filter((e) => e.ev).map((e) => e.ev)), stateAt = []; { let s = 'PRODUCING'; for (const e of order) { if (e.ev && OPS_STATES[s]?.[e.ev]) s = OPS_STATES[s][e.ev]; stateAt.push(s); } }

  // ---------- J. historical replay ----------
  const logRows = (Array.isArray(v.log) ? v.log : []).map((r) => ({ t: num(r.t, null), rate: num(r.rate, null), choke: num(r.choke, 100), pIn: num(r.pIn, null), tArr: num(r.tArr, null) })).filter((r) => r.t !== null && r.rate !== null).sort((a, b) => a.t - b.t);
  const iOut = nSt - 1, dpChokeOf = (q, z) => { const w = (sm0.p.wG + sm0.p.wL) * (q / Math.max(rate, 1e-9)); return (w / (Kc * Math.max(z, 0.02))) ** 2 / Math.max(eqT.y[3] / sm0.p.Vr * 0.6 + 30, 30) / 1e5; };
  const replayModel = { pIn: (q) => sur.pIn.predict(q), tArr: (q) => sur.tArr.predict(q), dpChoke: dpChokeOf, qLo, qHi, pSettle: so.pSettle, tAmbOut: tAmbOf(iOut), tauWarm: Math.max(restartTime, 0.5) * HOUR / 3, tauCool: cool[iOut].tau };
  let replay = null;
  if (logRows.length >= 3) {
    const rp = replayLog(logRows, replayModel), wp = logRows.map((r, i) => i).filter((i) => logRows[i].pIn !== null), wt2 = logRows.map((r, i) => i).filter((i) => logRows[i].tArr !== null);
    const mP = wp.length >= 3 ? metrics(wp.map((i) => logRows[i].pIn), wp.map((i) => rp.pIn[i])) : null, mT = wt2.length >= 3 ? metrics(wt2.map((i) => logRows[i].tArr), wt2.map((i) => rp.tArr[i])) : null;
    const full = wp.length === logRows.length && wt2.length === logRows.length && v.useResidual, rP = full ? ridgeResidual(logRows, rp.pIn, 'pIn') : null, rT = full ? ridgeResidual(logRows, rp.tArr, 'tArr') : null;
    replay = { rp, mP, mT, rP, rT };
  }

  // ---------- comparison with operating records (commissioning, start-up, shutdown, restart, ESD, blowdown, cooldown, pigging, tracer, …) ----------
  const recRows = (Array.isArray(v.records) ? v.records : []).map((r) => ({ kind: String(r.kind || '').trim().toLowerCase().replace(/[^a-z]/g, ''), t: num(r.t, 0), value: num(r.value, null), tag: String(r.tag || '').trim() })).filter((r) => r.kind && r.value !== null);
  const ser = (xs, ys) => (r) => interp1(xs, ys, clamp(r.t, xs[0], xs[xs.length - 1])), afterTrip = (m, tTrip) => { const i0 = m.t.findIndex((t) => t >= tTrip), n0 = m.n[Math.max(i0, 0)] || 1; return ser(m.t.slice(i0).map((t) => t - tTrip), m.n.slice(i0).map((x) => (100 * x) / n0)); };
  const tauEq = clamp((so.gasMass * volScale * Math.max(pIn0 - line.pOut, 0)) / (4 * Math.max(so.pSettle, 1) * Math.max(pm.mG, 1e-3)), 5, 6 * HOUR), pigRateT = pigVsRate.length > 1 ? (r) => interp1(pigVsRate.map((p) => 100 * p.q), pigVsRate.map((p) => p.transit), clamp(r.t, 100 * pigVsRate[0].q, 100 * pigVsRate[pigVsRate.length - 1].q)) : () => pigTransit ?? 0;
  const REC = {
    commissioning: { label: 'Commissioning: inlet pressure against rate', x: 'rate (% of case)', unit: 'bara', f: (r) => sur.pIn.predict(qClamp(r.t / 100)) },
    startup: { label: 'Start-up from cold: arrival temperature', x: 'h since start', unit: '°C', f: ser(cold.t.map((t) => t / HOUR), cold.Tout) },
    shutdown: { label: 'Shutdown: line pressure', x: 'h since shut-in', unit: 'bara', f: ser(tH, pT) },
    restart: { label: 'Restart: arrival temperature', x: 'h since restart', unit: '°C', f: ser(tW, warm.Tout) },
    esd: { label: 'Emergency shutdown: inlet pressure', x: 's since the trip', unit: 'bara', f: (r) => so.pSettle + (pIn0 - so.pSettle) * Math.exp(-Math.max(r.t, 0) / tauEq) },
    blowdown: { label: 'Blowdown test: line pressure', x: 'h since the valve opened', unit: 'bara', f: ser(bd.t.map((t) => t / HOUR), bd.P.map((p) => p / 1e5)) },
    cooldown: { label: 'Cooldown: cold-spot fluid temperature', x: 'h since shut-in', unit: '°C', f: ser(tH, cool[iCold].Tf) },
    pigging: { label: 'Pigging record: pig position', x: 'h since launch', unit: 'km', f: ser(pig.t.map((t) => t / HOUR), pig.x.map((x) => x / 1000)) },
    pigarrival: { label: 'Pig arrival time against rate', x: 'rate (% of case)', unit: 'h', f: pigRateT },
    surge: { label: 'Liquid surging: arrival rate during ramp-up', x: 'h since restart', unit: 'm³/h', f: ser(ramp.t.map((t) => t / HOUR), ramp.qOut.map((q) => q * 3600)) },
    tracer: { label: 'Chemical tracer: outlet concentration', x: 'h since injection', unit: 'C/C0', f: ser(front.tAll.map((t) => t / HOUR), front.outlet) },
    concentration: { label: 'MEG / methanol concentration at arrival', x: 'injection rate (m³/d)', unit: 'wt %', f: (r) => { const mI = (Math.max(r.t, 0) / DAY) * rhoLean * lean * doseArgs.eff; return (100 * mI) / Math.max(mWater + mI / lean, 1e-9); } },
    valve: { label: 'Valve response to a step command', x: 's since the command', unit: 'fraction of travel', f: (r) => (r.t <= v.actDead ? 0 : 1 - Math.exp(-(r.t - v.actDead) / Math.max(v.tauValve, 0.5))) },
    compressor: { label: 'Compressor transient: speed after a trip', x: 's since the trip', unit: '% of speed at trip', f: afterTrip(comp, compTrip) },
    pump: { label: 'Pump transient: speed after a trip', x: 's since the trip', unit: '% of speed at trip', f: afterTrip(pump, pumpTrip) },
    level: { label: 'Separator level history', x: 'h', unit: '%', f: ser(pidRun.t.map((t) => t / HOUR), pidRun.level.map((l) => 100 * l)) },
    alarm: { label: 'Field alarm / event history (1 = raised)', x: 'tag', unit: '–', f: (r) => (alarms.some((a) => a.tag === r.tag) ? 1 : 0) },
  };
  const ALIAS = { commission: 'commissioning', start: 'startup', esdtrip: 'esd', emergencyshutdown: 'esd', trip: 'esd', depressurisation: 'blowdown', depressurization: 'blowdown', pig: 'pigging', pigtime: 'pigarrival', surging: 'surge', meg: 'concentration', methanol: 'concentration', megconcentration: 'concentration', valveresponse: 'valve', separatorlevel: 'level', event: 'alarm', alarms: 'alarm' };
  const recCmp = recRows.map((r) => { const k = REC[r.kind] ? r.kind : ALIAS[r.kind]; if (!k) return { ...r, known: false, model: null }; let m = null; try { m = +REC[k].f(r); } catch { m = null; } return { ...r, kind: k, known: true, model: isNum(m) ? m : null }; });
  const recStats = Object.keys(REC).map((k) => { const rr = recCmp.filter((r) => r.kind === k && r.model !== null); if (!rr.length) return null; const mt = metrics(rr.map((r) => r.value), rr.map((r) => r.model)); return { kind: k, n: rr.length, bias: mt.bias, rmse: mt.rmse, mape: mt.mape, agree: k === 'alarm' ? mean(rr.map((r) => (Math.round(r.value) === r.model ? 1 : 0))) : null }; }).filter(Boolean);
  const recUnknown = [...new Set(recCmp.filter((r) => !r.known).map((r) => r.kind))];

  // ---------- results ----------
  progress(0.95, 'Assembling results');
  const km = (x) => x / 1000, thin = (a, n = 200) => { const s = Math.max(1, Math.ceil(a.length / n)); return a.filter((_, i) => i % s === 0 || i === a.length - 1); }, iRb = st.s.findIndex((s) => s >= sm0.p.Lp), iRbU = iRb < 0 ? nSt - 1 : Math.max(iRb - 1, 0);
  const plots = [], tables = [], kpis = [];
  { // cooldown curves
    const pickSt = [...new Set([0, iCold, iRbU, nSt - 1])], series = pickSt.map((i) => ({ name: `${i === iCold ? 'Cold spot, ' : ''}${km(st.x[i]).toFixed(1)} km`, x: thin(tH), y: thin(cool[i].Tf) }));
    series.push({ name: 'Hydrate temperature + margin', x: thin(tH), y: thin(limT), dash: true }, { name: 'Lumped-capacitance check (cold spot)', x: thin(tH), y: thin(tH.map((t) => cool[iCold].lumped(t * HOUR))), dash: true });
    if (coldComp) series.push({ name: coldComp.name, x: thin(tH), y: thin(coldComp.series) });
    plots.push({ type: 'line', title: 'Cooldown after shut-in', xlabel: 'Time since shut-in (h)', ylabel: 'Fluid temperature (°C)', series, hlines: [{ y: v.wat, label: 'WAT' }], vlines: neverCools ? [] : [{ x: cooldownTime, label: 'cooldown time' }, { x: noTouch, label: 'no-touch' }] });
    plots.push({ type: 'line', title: 'Cooldown time along the line (cold-spot map)', xlabel: 'Distance (km)', ylabel: 'Time to hydrate temperature + margin (h)', zeroY: true, series: [{ name: 'Line pipe', x: st.x.map(km), y: tCool.map((t) => t ?? v.tHorizon) }, ...(comps.length ? [{ name: 'Special components', x: comps.map((c) => km(c.x)), y: comps.map((c) => c.t ?? v.tHorizon), mode: 'points' }] : [])], hlines: [{ y: v.tHorizon, label: 'simulated horizon (never inside)' }], note: 'Points at the horizon never reach the hydrate temperature within the simulated time.' });
    const rowsT = Array.from({ length: Math.min(49, ntCool + 1) }, (_, r) => Math.round((r * ntCool) / (Math.min(49, ntCool + 1) - 1)));
    plots.push({ type: 'field', title: 'Temperature during cooldown (distance–time)', xlabel: 'Distance (km)', ylabel: 'Time since shut-in (h)', zlabel: 'Fluid temperature', zunit: '°C', x: st.x.map(km), y: rowsT.map((k) => tH[k]), z: rowsT.map((k) => cool.map((c) => c.Tf[k])), cmap: 'thermal', contours: 8, markers: neverCools ? [] : [{ x: km(coldSpotX), y: cooldownTime, label: 'first hydrate risk' }] });
  }
  plots.push({ type: 'line', title: 'Warm-up after restart', xlabel: 'Time since restart (h)', ylabel: 'Temperature (°C)', series: [{ name: 'Arrival temperature', x: thin(tW), y: thin(warm.Tout) }, { name: 'Smallest margin to hydrate temperature + margin', x: thin(tW), y: thin(warm.Tmin), dash: true }], hlines: [{ y: 0, label: 'hydrate-safe' }], vlines: [{ x: restartTime, label: 'steady' }, ...(tSafe !== null ? [{ x: tSafe, label: 'hydrate-safe' }] : [])] });
  plots.push({ type: 'line', title: 'Liquid surge during ramp-up', xlabel: 'Time since restart (h)', ylabel: 'm³/h · m³', series: [{ name: 'Liquid arriving (m³/h)', x: thin(ramp.t.map((t) => t / HOUR)), y: thin(ramp.qOut.map((q) => q * 3600)) }, { name: 'Volume above drain capacity (m³)', x: thin(ramp.t.map((t) => t / HOUR)), y: thin(ramp.V) }, { name: 'Rate (% of case)', x: thin(ramp.t.map((t) => t / HOUR)), y: thin(ramp.q.map((q) => 100 * q)), dash: true }, { name: 'Liquid inventory in the line (m³ / 10)', x: thin(ramp.t.map((t) => t / HOUR)), y: thin(ramp.inventory.map((x) => x / 10)), dash: true }], hlines: [{ y: v.qDrainM3h, label: 'drain capacity' }, { y: surgeAllow, label: 'surge allowance' }] });
  plots.push({ type: 'line', title: 'Blowdown: pressure and flare rate', xlabel: 'Time since the valve opened (h)', ylabel: 'bara · kg/s', series: [{ name: 'Line pressure (bara)', x: thin(bd.t.map((t) => t / HOUR)), y: thin(bd.P.map((p) => p / 1e5)) }, { name: 'Flare rate (kg/s)', x: thin(bd.t.map((t) => t / HOUR)), y: thin(bd.mdot) }], hlines: [{ y: Math.max(pSafe - headBd, 0), label: 'hydrate-safe at the top' }] });
  plots.push({ type: 'line', title: 'Blowdown: temperatures', xlabel: 'Time since the valve opened (h)', ylabel: 'Temperature (°C)', series: [{ name: 'Gas in the line', x: thin(bd.t.map((t) => t / HOUR)), y: thin(bd.T.map((T) => T - KEL)) }, { name: 'Wall and liquid', x: thin(bd.t.map((t) => t / HOUR)), y: thin(bd.Tw.map((T) => T - KEL)) }, { name: 'Downstream of the valve (Joule–Thomson)', x: thin(bd.t.map((t) => t / HOUR)), y: thin(bd.T.map((T, i) => T - KEL - bdCfg.jt * Math.max(bd.P[i] - bdCfg.pBack, 0))), dash: true }], hlines: [{ y: v.tMinDesign, label: 'minimum design temperature' }] });
  plots.push({ type: 'line', title: 'Pig position and liquid pushed ahead', xlabel: 'Time since launch (h)', ylabel: 'km · m³', series: [{ name: 'Pig position (km)', x: thin(pig.t.map((t) => t / HOUR)), y: thin(pig.x.map(km)) }, { name: 'Liquid collected ahead (m³)', x: thin(pig.t.map((t) => t / HOUR)), y: thin(pig.slug) }], hlines: [{ y: catVol, label: 'slug-catcher volume' }], vlines: pig.tFront !== null ? [{ x: pig.tFront / HOUR, label: 'slug front arrives' }] : [] });
  if (waxMax > 0) plots.push({ type: 'line', title: 'Wax inventory in the line between pig runs', xlabel: 'Time (d)', ylabel: 'Deposit volume (m³)', zeroY: true, series: [{ name: 'Wax in the line', x: waxSeries.t, y: waxSeries.v }], note: `Deposit grows to ${v.waxThk} mm over the pigging interval on the length below the WAT; each run removes ${v.pigEff} %.` });
  if (pigVsRate.length > 1) plots.push({ type: 'line', title: 'Pig velocity and receiver surge against rate', xlabel: 'Rate (% of case)', ylabel: 'm/s · m³', series: [{ name: 'Mean pig velocity (m/s)', x: pigVsRate.map((p) => 100 * p.q), y: pigVsRate.map((p) => p.v), mode: 'both' }, { name: 'Surge above drain capacity (m³ / 10)', x: pigVsRate.map((p) => 100 * p.q), y: pigVsRate.map((p) => (p.surge * volScale) / 10), mode: 'both' }], hlines: [{ y: 5, label: 'upper pig velocity guide (5 m/s)' }] });
  { const kk = [0.25, 0.5, 0.75, 1].map((f) => Math.min(front.t.length - 1, Math.max(1, Math.round(f * (front.t.length - 1) * 0.5)))), ser = []; for (const k of kk) { ser.push({ name: `t = ${(front.t[k] / HOUR).toFixed(2)} h`, x: front.xc.map(km), y: front.c[k].slice() }); } const kA = kk[1]; ser.push({ name: `Analytic (mean velocity), t = ${(front.t[kA] / HOUR).toFixed(2)} h`, x: front.xc.map(km), y: front.xc.map((x) => frontAnalytic(x, front.t[kA], uMean, disp)), dash: true });
    plots.push({ type: 'line', title: 'Inhibitor front along the line', xlabel: 'Distance (km)', ylabel: 'Concentration / injected concentration', ymin: 0, ymax: 1.05, series: ser, hlines: [{ y: 0.95, label: 'protected' }] }); }
  { const zs = linspace(Math.max(0.03, sm0.zFloor), 1, 24), gr = zs.map((z) => sm0.growth(z) * HOUR); plots.push({ type: 'line', title: 'Open-loop stability of the riser against choke opening', xlabel: 'Choke opening (%)', ylabel: 'Growth rate of the slowest mode (1/h)', series: [{ name: 'Largest real part of the poles', x: zs.map((z) => 100 * z), y: gr }], hlines: [{ y: 0, label: 'stability limit' }], vlines: [...(zCrit !== null ? [{ x: 100 * zCrit, label: 'limit cycle starts' }] : []), { x: 100 * zTarget, label: 'operating opening' }] }); }
  { const th = pidRun.t.map((t) => t / HOUR); plots.push({ type: 'line', title: 'Inlet pressure: closed loop, then controller in manual', xlabel: 'Time (h)', ylabel: 'Inlet pressure (bara)', series: [{ name: `${v.ctlMode} control, then open loop`, x: thin(th, 400), y: thin(pidRun.pIn, 400) }, { name: 'Set-point', x: thin(th, 400), y: thin(pidRun.sp, 400), dash: true }], vlines: [{ x: ctlCfg.tOff / HOUR, label: 'controller to manual' }] });
    plots.push({ type: 'line', title: 'Choke opening and separator level', xlabel: 'Time (h)', ylabel: '%', series: [{ name: 'Choke opening (%)', x: thin(th, 400), y: thin(pidRun.z.map((z) => 100 * z), 400) }, { name: 'Separator level (%)', x: thin(th, 400), y: thin(pidRun.level.map((l) => 100 * l), 400) }], hlines: [{ y: 80, label: 'high level' }], vlines: [{ x: ctlCfg.tOff / HOUR, label: 'controller to manual' }] });
    const tm = mpcRun.t.map((t) => t / HOUR), nC = Math.min(tm.length, Math.round(ctlCfg.tOff / dtCtl) + 1); plots.push({ type: 'line', title: 'PID against model-predictive control (set-point step)', xlabel: 'Time (h)', ylabel: 'Inlet pressure (bara)', series: [{ name: v.ctlMode, x: thin(th.slice(0, nC), 300), y: thin(pidRun.pIn.slice(0, nC), 300) }, { name: 'Linear MPC + Kalman filter', x: thin(tm, 300), y: thin(mpcRun.pIn, 300) }, { name: 'Set-point', x: thin(tm, 300), y: thin(mpcRun.sp, 300), dash: true }] });
    plots.push({ type: 'line', title: 'Extended Kalman filter: riser liquid mass from the inlet, riser-base and topside pressures', xlabel: 'Time (h)', ylabel: 'Liquid mass in the riser (t)', series: [{ name: 'Model (truth)', x: idx.map((i) => th[i]), y: mTrue.map((m) => m / 1000) }, { name: 'EKF estimate', x: idx.map((i) => th[i]), y: mEst.map((m) => m / 1000), dash: true }] }); }
  { const q100 = scan.q.map((q) => 100 * q), vl = env.feasible ? [{ x: 100 * env.qMin, label: 'min' }, { x: 100 * env.qMax, label: 'max' }] : [];
    plots.push({ type: 'line', title: 'Operating envelope: pressure and temperature limits', xlabel: 'Rate (% of case)', ylabel: 'bara · °C', series: [{ name: 'Inlet pressure (bara)', x: q100, y: scan.pIn, mode: 'both' }, { name: 'Inlet pressure incl. stabilising choke (bara)', x: q100, y: scan.pReq, mode: 'both' }, { name: 'Arrival temperature (°C)', x: q100, y: scan.tArr, mode: 'both' }, { name: 'Smallest hydrate margin (°C)', x: q100, y: scan.margin, mode: 'both' }], hlines: [{ y: v.pAvail, label: 'available pressure' }, { y: v.wat + v.watMargin, label: 'WAT + margin' }, { y: v.hydMargin, label: 'hydrate margin' }], vlines: vl });
    plots.push({ type: 'line', title: 'Operating envelope: utilisation of the capacity limits', xlabel: 'Rate (% of case)', ylabel: 'Fraction of limit', zeroY: true, series: [{ name: 'Erosional velocity ratio', x: q100, y: scan.eros, mode: 'both' }, { name: 'Liquid rate / separator capacity', x: q100, y: scan.qLiq.map((q) => q / v.qDrainM3h), mode: 'both' }, { name: 'Critical choke opening (fraction)', x: q100, y: scan.zCrit, mode: 'both', dash: true }], hlines: [{ y: 1, label: 'limit' }], vlines: vl }); }
  { const ser = methods.filter((m) => Array.isArray(m.history) && m.history.length > 1).map((m) => ({ name: m.name, x: m.history.map((_, i) => i + 1), y: m.history.map((f) => -f * scaleP / 1e6) })); if (ser.length) plots.push({ type: 'line', title: 'Optimisation convergence', xlabel: 'Iteration / generation', ylabel: 'Penalised daily margin (M$/d)', series: ser, hlines: [{ y: profitOpt / 1e6, label: 'best' }] }); }
  if (replay) { const tt = replay.rp.t, ser = [{ name: 'Model inlet pressure (bara)', x: tt, y: replay.rp.pIn }, { name: 'Model arrival temperature (°C)', x: tt, y: replay.rp.tArr }]; const wp = logRows.filter((r) => r.pIn !== null), wtm = logRows.filter((r) => r.tArr !== null); if (wp.length) ser.push({ name: 'Logged inlet pressure', x: wp.map((r) => r.t), y: wp.map((r) => r.pIn), mode: 'points' }); if (wtm.length) ser.push({ name: 'Logged arrival temperature', x: wtm.map((r) => r.t), y: wtm.map((r) => r.tArr), mode: 'points' });
    if (replay.rP) ser.push({ name: 'Inlet pressure with residual correction', x: tt, y: replay.rp.pIn.map((p, i) => p + replay.rP.corr[i]), dash: true }); plots.push({ type: 'line', title: 'Historical replay of the operating log', xlabel: 'Log time (h)', ylabel: 'bara · °C', series: ser }); }

  { // rotating equipment
    const sp = compCurve.map ? compCurve.map.lines : [0.7, 0.8, 0.9, 1, 1.05].map((N) => ({ n: N, q: compCurve.q.map((x) => x * N), h: compCurve.h.map((x) => x * N * N) }));
    const iT = comp.t.findIndex((t) => t >= compTrip), pathQ = comp.q.slice(0, iT), pathH = comp.q.slice(0, iT).map((Q, i) => compCurve.head(Q, comp.n[i]) / 1000);
    plots.push({ type: 'line', title: 'Compressor map with the turndown test', xlabel: 'Inlet volume flow (m³/s)', ylabel: 'Polytropic head (kJ/kg)', series: [...sp.map((l) => ({ name: `${(100 * l.n).toFixed(0)} % speed`, x: l.q, y: l.h.map((h) => h / 1000) })), { name: 'Surge line', x: sp.map((l) => l.q[0]), y: sp.map((l) => l.h[0] / 1000), dash: true }, { name: `Surge control line (${v.compSmCtl} % margin)`, x: sp.map((l) => l.q[0] / (1 - v.compSmCtl / 100)), y: sp.map((l) => compCurve.head(l.q[0] / (1 - v.compSmCtl / 100), l.n) / 1000), dash: true }, { name: 'Operating point during the test', x: thin(pathQ, 150), y: thin(pathH, 150), mode: 'both' }], note: 'Speed lines follow the fan laws; the operating point moves from the duty towards the surge control line when the gas rate falls and is held there by the recycle.' });
    plots.push({ type: 'line', title: 'Compressor: speed, surge margin and recycle (turndown, then trip)', xlabel: 'Time (s)', ylabel: '%', series: [{ name: 'Speed (%)', x: thin(comp.t, 300), y: thin(comp.n.map((x) => 100 * x), 300) }, { name: 'Surge margin (%)', x: thin(comp.t, 300), y: thin(comp.sm.map((x) => Math.max(x, -20)), 300) }, { name: 'Recycle valve (%)', x: thin(comp.t, 300), y: thin(comp.recycle.map((x) => 100 * x), 300) }, { name: 'Gas feed (% of duty)', x: thin(comp.t, 300), y: thin(comp.t.map((t) => (100 * compFeed(t)) / mGas), 300), dash: true }], hlines: [{ y: v.compSmCtl, label: 'surge control line' }, { y: 0, label: 'surge' }], vlines: [{ x: 60, label: 'turndown' }, { x: compTrip, label: 'trip' }] });
    plots.push({ type: 'line', title: 'Export pump: flow control and coast-down', xlabel: 'Time (s)', ylabel: '% · m³/h', series: [{ name: 'Speed (%)', x: thin(pump.t, 300), y: thin(pump.n.map((x) => 100 * x), 300) }, { name: 'Flow (m³/h)', x: thin(pump.t, 300), y: thin(pump.q.map((x) => 3600 * x), 300) }, { name: 'Shaft power (kW)', x: thin(pump.t, 300), y: thin(pump.power.map((x) => x / 1000), 300), dash: true }], hlines: [{ y: pumpMinFlow * 3600, label: 'minimum continuous flow' }], vlines: [{ x: 30, label: 'flow set-point step' }, { x: pumpTrip, label: 'trip' }] });
  }
  plots.push({ type: 'line', title: 'Chemical storage level', xlabel: 'Time (d)', ylabel: 'Stock (m³)', zeroY: true, series: [{ name: 'Stock', x: thin(stock.t, 300), y: thin(stock.level, 300) }], hlines: [{ y: reorder, label: 're-order level' }, { y: v.tankVol, label: 'tank volume' }], note: `Continuous use ${contUsed.toFixed(1)} m³/d and ${batchEvent.toFixed(0)} m³ per shutdown (${nEv} in ${invDays} d); deliveries of ${v.resupplyVol} m³ arrive ${v.resupplyLead} d after the re-order level is reached.` });
  plots.push({ type: 'line', title: 'Ratio control of the inhibitor during the restart', xlabel: 'Time since restart (h)', ylabel: 'kg/s', series: [{ name: 'Inhibitor injected (ratio control)', x: thin(ratioRun.t.map((t) => t / HOUR), 300), y: thin(ratioRun.q, 300) }, { name: 'Requirement (ratio × water rate)', x: thin(ratioRun.t.map((t) => t / HOUR), 300), y: thin(ratioRun.sp, 300), dash: true }, { name: 'Water rate / 10', x: thin(ratioRun.t.map((t) => t / HOUR), 300), y: thin(ratioRun.wild.map((w) => w / 10), 300), dash: true }], hlines: [{ y: ratioKg * mWater, label: 'fixed injection sized for full rate' }] });
  if (nm) plots.push({ type: 'line', title: 'Nonlinear model-predictive control on the nonlinear riser model', xlabel: 'Time (h)', ylabel: 'bara · %', series: [{ name: 'Inlet pressure (bara)', x: nm.t.map((t) => t / HOUR), y: nm.pIn }, { name: 'Set-point (bara)', x: nm.t.map((t) => t / HOUR), y: nm.sp, dash: true }, { name: 'Choke opening (%)', x: nm.t.map((t) => t / HOUR), y: nm.z.map((z) => 100 * z), mode: 'step' }], note: `Single shooting over 5 samples of ${Math.max(ctlCfg.tsMpc, 120)} s, two free moves, solved by sequential quadratic programming (${nm.iterations.toFixed(1)} iterations and ${(nm.evals / nm.t.length).toFixed(0)} model predictions per move on average).` });
  if (recStats.length) { const kinds = recStats.filter((s) => s.kind !== 'alarm').map((s) => s.kind); if (kinds.length) plots.push({ type: 'line', title: 'Operating records against the model (each scaled by its largest measured value)', xlabel: 'Measured / largest measured', ylabel: 'Model / largest measured', series: [...kinds.map((k) => { const rr = recCmp.filter((r) => r.kind === k && r.model !== null), sc = Math.max(...rr.map((r) => Math.abs(r.value)), 1e-12); return { name: REC[k].label, x: rr.map((r) => r.value / sc), y: rr.map((r) => r.model / sc), mode: 'points' }; }), { name: 'Perfect agreement', x: [0, 1], y: [0, 1], dash: true }] }); }

  // tables
  tables.push({ title: 'Shutdown – preservation – restart sequence', columns: ['Time (h)', 'Tag', 'Action', 'State after'], rows: order.map((e, i) => [rd(e.t, 2), e.tag, e.action, stateAt[i]]), note: `Scheduler fired ${sched.fired.length} events in ${sched.steps} steps; inhibitor injected in the sequence ${acc.inj.toFixed(1)} m³; production lost ${(acc.lost / HOUR).toFixed(1)} full-rate hours. State machine ended in ${fsm.state}${fsm.rejected ? ` (${fsm.rejected} event(s) not permitted in the current state)` : ''}.` });
  tables.push({ title: 'Alarms and interlocks', columns: ['Tag', 'Variable', 'Type', 'Limit', 'Value', 'Level', 'Status', 'Action'], rows: rulesTbl.map((r) => { const val = almVals[r.key], hit = raised.find((a) => a.tag === String(r.tag)); return [String(r.tag), String(r.key), String(r.type), txt(+r.limit), txt(val), String(r.level || 'alarm'), isNum(val) ? (hit ? 'RAISED' : 'normal') : 'not evaluated', String(r.action || '—')]; }), note: `Fail-safe state on a trip: ${safeState === 'FAILSAFE' ? 'FAILSAFE is active — ' : ''}${Object.entries(FAILSAFE_POSITIONS).map(([k, x]) => `${k}: ${x}`).join('; ')}.` });
  tables.push({ title: 'Controller tunings (inlet pressure → choke)', columns: ['Rule', 'Mode', 'Kc at test point (%/bar)', 'Kc at target (%/bar)', 'Ti (s)', 'Td (s)', 'Stabilises target'], rows: [...tuneRows.map((r) => [r.rule, r.mode, rd(100 * r.kc, 3), rd(100 * r.kcT, 3), rd(r.ti, 0), rd(r.td, 1), r.stable ? 'yes' : 'no']), ['Closed-loop pole search', pole.ti > 0 ? 'PI' : 'P', '—', rd(100 * pole.kc, 3), rd(pole.ti, 0), 0, pole.decay > 0 ? 'yes' : 'no'], ['APPLIED: ' + sel.rule, v.ctlMode, '—', rd(100 * sel.kc, 3), rd(sel.ti, 0), rd(sel.td, 1), linStable ? 'yes' : 'no']],
    note: `Step test at ${(100 * zId).toFixed(1)} % opening: K = ${fo.K.toFixed(1)} bar per unit opening, τ = ${(fo.tau / 60).toFixed(1)} min, θ = ${(foUse.theta).toFixed(0)} s (fit error ${fo.rmse.toFixed(3)} bar). Ultimate gain ${(100 * rules.ku).toFixed(2)} %/bar, period ${(rules.pu / 60).toFixed(1)} min. Gain scheduling by the static-gain ratio ${gainRatio.toFixed(2)} to the target opening. Loop margins of the SIMC PI settings on the identified model: gain margin ${marg.gm === null ? '∞' : marg.gm.toFixed(2)}, phase margin ${marg.pm === null ? '—' : marg.pm.toFixed(0) + '°'}, sensitivity peak Ms ${marg.ms.toFixed(2)}. The robust row is the PI setting with the largest integral gain whose Ms stays below the limit for the nominal model and for the models with the gain and dead time of the uncertainty set; the self-tuning row applies the SIMC rule to a first-order model identified by recursive least squares from the same step test (K = ${rl.K.toFixed(0)} bar per unit opening, τ = ${(rl.tau / 60).toFixed(1)} min).` });
  const fmtC = (e) => `${(e[0] * HOUR).toFixed(2)}${Math.abs(e[1]) > 1e-9 ? ' ± ' + Math.abs(e[1] * HOUR).toFixed(2) + 'j' : ''}`, uniq = (pl) => pl.filter((e) => e[1] >= -1e-12);
  tables.push({ title: 'Linear model at the operating opening', columns: ['Item', 'Value'], rows: [['Open-loop poles (1/h)', uniq(polesOL).map(fmtC).join(', ')], ['Closed-loop poles with the applied controller (1/h)', uniq(clPoles).map(fmtC).join(', ')], ['Transfer function numerator (bar per unit opening, descending powers of s)', tf.num.map((x) => x.toExponential(3)).join(', ')], ['Transfer function denominator', tf.den.map((x) => x.toExponential(3)).join(', ')], ['Static gain (bar per unit opening)', rd(kStat, 2)], ['Finite-horizon LQ gain on the scaled states', lqGain ? lqGain.map((x) => x.toExponential(2)).join(', ') : '—'], ['MPC moves with active constraints (%)', rd(100 * mpcRun.mpcActive, 1)], ['EKF error on riser liquid mass (% of mean)', rd(100 * ekfRel, 2)], [`EKF / UKF error over the first ${nUk} samples (% of mean)`, `${ekfRelWin === null ? '—' : (100 * ekfRelWin).toFixed(2)} / ${ukfRel === null ? '—' : (100 * ukfRel).toFixed(2)}`], ['Sensitivity peak Ms / complementary peak Mt of the applied loop', `${loopLin.ms.toFixed(2)} / ${loopLin.mt.toFixed(2)}`], ['Recursive-least-squares model of the step test (K bar per unit opening, τ min)', `${rl.K.toFixed(1)}, ${(rl.tau / 60).toFixed(1)}`], ['Nonlinear MPC: integral absolute error (bar·h) / remaining swing (bar)', nm ? `${nm.iae.toFixed(3)} / ${nm.amp.toFixed(2)}` : '—']], note: 'State vector: gas and liquid mass in the feed pipeline, gas and liquid mass in the riser; input: choke opening; outputs: inlet and topside pressure.' });
  tables.push({ title: 'Operating-envelope limits', columns: ['Constraint', 'Type', 'Limit', 'Value at the operating rate', 'Margin', 'Bounds the window', 'At rate (%)'], rows: env.limits.map((l) => [l.name, l.type === 'min' ? '≥' : '≤', txt(l.limit), txt(l.value), txt(l.margin), l.bound, l.q === null ? '—' : rd(100 * l.q, 0)]), note: env.text.join('; ') + '.' });
  tables.push({ title: 'Chemical injection summary', columns: ['Item', 'Value', 'Unit'], rows: [['Inhibitor', inh.name, ''], ['Governing condition', v.dosingBasis === 'steady' ? 'steady flow' : v.dosingBasis === 'max' ? 'worse of steady flow and shutdown' : 'shutdown cold spot', ''], ['Required depression (incl. margin)', rd(dTgov, 1), '°C'], ['Dose in the aqueous phase', rd(doseWt, 1), 'wt %'], ['Dose attainable with this chemical', doseGov.attainable ? 'yes' : 'no', ''], ['Continuous injection required at the operating rate', rd(contRate, 1), 'm³/d'], ['Injection setting used', rd(contUsed, 1), 'm³/d'], ['Concentration reached with that setting', rd(doseAchieved, 1), 'wt %'], ['Injection point', rd(sInj / 1000, 2), 'km from the inlet'], ['Pump utilisation', rd(100 * pumpUtil, 0), '%'], ['Loss to the gas phase', rd(doseSteady.lossGas * DAY / 1000, 3), 't/d'], ['Loss to the hydrocarbon liquid', rd(doseSteady.lossOil * DAY / 1000, 3), 't/d'], ['Batch for a shutdown (bullheading)', rd(batchVol, 1), 'm³'], ['Bullheading time at pump capacity', rd(tBullhead, 2), 'h'], ['Injection during restart until hydrate-safe', rd(restartInj, 1), 'm³'], ['Time for the front to protect the whole line', rd(tProtect, 2), 'h'], ['Axial dispersion coefficient', rd(disp, 2), 'm²/s'], ['Glycol loop: rich stream', megLoop ? rd(megLoop.rich, 2) : '—', 'kg/s'], ['Glycol loop: inventory', megLoop ? rd(megLoop.inventory, 0) : '—', 'm³'], ['Glycol loop: regeneration duty', megLoop ? rd(megLoop.duty, 0) : '—', 'kW']], note: ldhi });
  tables.push({ title: 'Cooldown of the line and special components', columns: ['Location', 'Distance (km)', 'Start temperature (°C)', 'Time to hydrate limit (h)', 'Time to WAT (h)', 'Temperature at restart (°C)'], rows: [...[...new Set([0, Math.floor(nSt / 4), iMid, iCold, iRbU, nSt - 1])].sort((a, b) => a - b).map((i) => [i === iCold ? 'Line (cold spot)' : 'Line', rd(km(st.x[i]), 2), rd(st.T[i], 1), tCool[i] === null ? `> ${v.tHorizon}` : rd(tCool[i], 1), tWat[i] === null ? `> ${v.tHorizon}` : rd(tWat[i], 1), rd(tAtShut[i], 1)]), ...comps.map((c) => [c.name, rd(km(c.x), 2), rd(c.T0, 1), c.t === null ? `> ${v.tHorizon}` : rd(c.t, 1), (() => { const t = timeBelow(tH, c.series, v.wat); return t === null ? `> ${v.tHorizon}` : rd(t, 1); })(), rd(interp1(tH, c.series, Math.min(v.tShut, v.tHorizon)), 1)])], note: `Settle-out pressure ${so.pSettle.toFixed(1)} bara; U-value ${cool[iCold].U.toFixed(2)} W/m²K at shut-in conditions; effective thermal mass at the cold spot ${(cool[iCold].cEff / 1000).toFixed(0)} kJ/m/K; lumped time constant ${(cool[iCold].tau / HOUR).toFixed(1)} h.` });
  tables.push({ title: 'Preservation strategies for this shutdown', columns: ['Strategy', 'Feasible', 'Lead time (h)', 'Cost per shutdown ($)', 'Selected', 'Lowest cost'], rows: strategies.map((s) => [s.name, s.feasible ? (meets(s) ? 'yes' : `yes, but leaves less than ${v.noTouchMin} h no-touch time`) : 'no', rd(s.lead, 2), rd(s.cost, 0), s.id === stratSel.id ? '●' : '', stratBest && s.id === stratBest.id ? '●' : '']) });
  tables.push({ title: 'Optimisation of the operating point', columns: ['Method', 'Rate (% of case)', 'Dose (wt %)', 'Daily margin (M$/d)', 'Feasible', 'Evaluations / iterations'], rows: [...methods.map((m) => (m.x ? [m.name, rd(100 * m.x[0], 1), rd(m.x[1], 1), rd(profit(m.x[0], m.x[1], bestC.heat, bestC.pigs) / 1e6, 4), gOf(m.x[0], m.x[1], bestC.heat, bestC.pigs).every((g) => g > -2e-2) ? 'yes' : 'no', m.evals] : [m.name, '—', '—', '—', m.note || '—', 0])), ...combos.map((c) => [`Enumeration: heating ${c.heat ? 'on' : 'off'}, pigging programme ${c.pigs ? 'on' : 'off'}`, rd(100 * c.x[0], 1), rd(c.x[1], 1), rd(c.profit / 1e6, 4), c.feasible ? 'yes' : 'no', c.evals]), ['Mixed-integer LP (branch and bound, linearised)', milp.status === 'optimal' ? rd(100 * (qLo + milp.x[0]), 1) : '—', milp.status === 'optimal' ? rd(milp.x[1], 1) : '—', milp.status === 'optimal' ? rd(profit(qLo + milp.x[0], milp.x[1], milp.x[2], milp.x[3]) / 1e6, 4) : '—', milp.status === 'optimal' ? `heating ${milp.x[2] ? 'on' : 'off'}, pigging ${milp.x[3] ? 'on' : 'off'}` : milp.status, milp.nodes]],
    note: `Objective: oil revenue minus inhibitor and heating cost per day, subject to hydrate margin, wax, erosion, inlet-pressure and separator limits evaluated on the response surfaces. Ramp optimisation (three segments): ${(rampOpt.x[0]).toFixed(1)} h ramp with ${(100 * rampOpt.x[1]).toFixed(0)} % and ${(100 * Math.max(rampOpt.x[2], rampOpt.x[1])).toFixed(0)} % of the rise after one and two thirds, objective ${rampOpt.f.toFixed(2)} against ${rampLin.toFixed(2)} for the linear ramp.` });
  tables.push({ title: 'Response-surface fit of the mechanistic model', columns: ['Quantity', 'RMS error', 'Largest error', 'Leave-one-out RMS error', 'Unit'], rows: [['Inlet pressure', rd(sur.pIn.rmse, 3), rd(sur.pIn.maxErr, 3), rd(sur.pIn.loo, 3), 'bar'], ['Arrival temperature', rd(sur.tArr.rmse, 3), rd(sur.tArr.maxErr, 3), rd(sur.tArr.loo, 3), '°C'], ['Hydrate margin', rd(sur.margin.rmse, 3), rd(sur.margin.maxErr, 3), rd(sur.margin.loo, 3), '°C'], ['Liquid inventory', rd(sur.inv.rmse, 2), rd(sur.inv.maxErr, 2), rd(sur.inv.loo, 2), 'm³'], ['Erosional ratio', rd(sur.eros.rmse, 4), rd(sur.eros.maxErr, 4), rd(sur.eros.loo, 4), '–']], note: `Cubic polynomials in rate fitted to ${scan.q.length} kernel solutions between ${(100 * qLo).toFixed(0)} and ${(100 * qHi).toFixed(0)} % of the case rate.` });
  if (replay) tables.push({ title: 'Historical replay: agreement with the operating log', columns: ['Quantity', 'Points', 'Bias', 'RMS error', 'R²', 'Hold-out RMS before correction', 'Hold-out RMS after correction'], rows: [['Inlet pressure (bar)', replay.mP ? replay.mP.n : 0, replay.mP ? txt(replay.mP.bias, 2) : '—', replay.mP ? txt(replay.mP.rmse, 2) : '—', replay.mP ? txt(replay.mP.r2, 3) : '—', replay.rP ? txt(replay.rP.before, 2) : '—', replay.rP ? txt(replay.rP.after, 2) : '—'], ['Arrival temperature (°C)', replay.mT ? replay.mT.n : 0, replay.mT ? txt(replay.mT.bias, 2) : '—', replay.mT ? txt(replay.mT.rmse, 2) : '—', replay.mT ? txt(replay.mT.r2, 3) : '—', replay.rT ? txt(replay.rT.before, 2) : '—', replay.rT ? txt(replay.rT.after, 2) : '—']], note: 'The model is driven by the logged rate and choke opening only. The residual correction is a ridge regression on rate and choke opening fitted on every second row and tested on the others.' });
  tables.push({ title: 'Blowdown and pigging summary', columns: ['Item', 'Value', 'Unit'], rows: [['Gas volume blown down', rd(bdCfg.V, 0), 'm³'], ['Start pressure / temperature', `${pB.toFixed(1)} / ${TgB.toFixed(1)}`, 'bara / °C'], ['Time to the target pressure', bdReached ? rd(blowdownTime, 2) : `> ${blowdownTime.toFixed(1)}`, 'h'], ['Time to the hydrate-safe pressure at the top', bd.tMark === null ? '—' : rd(bd.tMark / HOUR, 2), 'h'], ['Peak flare rate', `${bd.peak.toFixed(1)} kg/s (${peakStd.toFixed(2)} MSm³/d)`, ''], ['Gas discharged / liberated from the oil', `${(bd.discharged / 1000).toFixed(1)} / ${(bd.flashed / 1000).toFixed(1)}`, 't'], ['Liquid left in the line', rd((mLiqLine - bd.liquidOut) / so.rhoL, 0), 'm³'], ['Route', bothEnds ? 'both ends (topside valve + inlet service line)' : 'topside valve only', ''], ['Pressure kept by gas trapped behind liquid seals (topside only / both ends)', `${headTop.toFixed(1)} / ${headBoth.toFixed(1)}`, 'bar'], ['Discharge model', v.bdMode === 'hrm' ? `two-phase with relaxation of the gas liberation (Θ = ${isNum(bd.theta) && bd.theta < 1e6 ? bd.theta.toPrecision(3) + ' s at the end' : 'not reached'}; held back ${bd.unreleased.toFixed(0)} kg)` : v.bdMode === 'hem' ? 'two-phase, homogeneous equilibrium' : 'gas only', ''], ['Seabed pressure after blowdown', rd(seabedPAfter, 1), 'bara'], ['Hydrate-safe pressure at ambient', pSafeRaw === null ? 'no hydrate at ambient' : rd(pSafe, 1), 'bara'], ['Equivalent valve Cv', rd((v.cdBlow * bdCfg.area) / 1.7e-5, 1), 'US gpm/psi^0.5'], ['Pig differential pressure (seal friction ' + (pigDpFric / 1e5).toFixed(2) + ' bar)', rd(pig.dpPig / 1e5, 2), 'bar'], ['Pig mean velocity', rd(pig.vMean, 2), 'm/s'], ['Liquid swept / leaked past the pig', `${(pig.swept * volScale).toFixed(0)} / ${(pig.leaked * volScale).toFixed(1)}`, 'm³'], ['Liquid arrival duration', rd(pig.duration / 60, 1), 'min'], ['Peak inlet pressure while pigging', rd(pInPigMax, 1), 'bara'], ['Wax removed per run', rd(waxVol * clamp(v.pigEff, 0, 100) / 100, 2), 'm³'], ['Largest wax inventory in the line', rd(waxMax, 2), 'm³'], ['Hot-oil circulation time to the hold temperature', hotOilTime === null ? '> 30' : rd(hotOilTime, 1), 'h']] });

  tables.push({ title: 'Slug control: automatic against manual', columns: ['Case', 'Choke opening (%)', 'Inlet pressure (bara)', 'Pressure swing (bar)', 'Peak separator level (%)', 'Carry-over (m³)'], rows: [
    ...(manual ? [] : [[`Automatic (${v.ctlMode}, ${sel.rule})`, rd(100 * zTarget, 1), rd(pidRun.meanP, 1), rd(pidRun.ampClosed, 2), rd(100 * pidRun.maxLevelClosed, 0), rd(pidRun.carryOver - pidRun.carryOverOpen, 1)]]),
    [`Manual at the same opening${manual ? ' (selected)' : ' (comparison)'}`, rd(100 * zTarget, 1), rd(eqT.Pp / 1e5, 1), rd(manual ? pidRun.ampClosed : pidRun.ampOpen, 2), rd(100 * (manual ? pidRun.maxLevel : pidRun.maxLevelOpen), 0), rd(manual ? pidRun.carryOver : pidRun.carryOverOpen, 1)],
    ...(zCrit !== null ? [['Manual at 90 % of the critical opening (stable)', rd(90 * zCrit, 1), rd(sm0.steady(Math.max(0.9 * zCrit, sm0.zFloor)).Pp / 1e5, 1), 0, rd(clamp(v.levelSp, 10, 90), 0), 0]] : []),
    ...(nm ? [['Nonlinear MPC (state feedback, riser model only)', rd(100 * mean(nm.z.slice(-10)), 1), rd(mean(nm.pIn.slice(-10)), 1), rd(nm.amp, 2), '—', '—']] : [])],
    note: `Surge capacity ${catVol.toFixed(0)} m³ (${v.catcherAuto && catReq > v.slugCatcherVol ? 'sized by the design from ' + v.slugCatcherVol.toFixed(1) + ' m³ entered' : 'as entered'}; requirement ${catReq.toFixed(0)} m³ for a design surge of ${designSurge.toFixed(1)} m³ with ${(100 * usable).toFixed(0)} % usable). Loop-shaping check of the applied controller on the linearised riser with ${(loopDead + 0.5 * dtCtl).toFixed(0)} s dead time: Ms = ${loopLin.ms.toFixed(2)}, Mt = ${loopLin.mt.toFixed(2)}${linStable ? '' : ' (closed loop not stable: the peaks are not meaningful)'}; an open-loop unstable riser cannot reach the Ms of 1.2–2 usual for stable processes.${rob ? ` Robust PI on the step-test model: Kc ${(100 * rob.kc).toFixed(2)} %/bar, Ti ${rob.ti.toFixed(0)} s, worst-case Ms ${rob.worst.toFixed(2)} over gain ±${v.robGainPct} % and dead time +${v.robDelayPct} %.` : ''}${schedule ? ' Gain scheduling on the choke opening is active.' : ''}` });
  tables.push({ title: 'Rotating equipment: operating points, limits and control tests', columns: ['Item', 'Value', 'Limit', 'Unit', 'Status'], rows: [
    ['Compressor design flow / head', `${compQd.toFixed(3)} m³/s / ${(compHd / 1000).toFixed(1)} kJ/kg`, '—', '', compTbl.filter((r) => r.q > 0 && r.h > 0).length >= 3 ? 'map as entered' : compCurve.map && typeof NET.compressorMap === 'function' ? 'map from the network suite generator' : 'built-in fan-law map'],
    ['Compressor suction / discharge', `${v.sepP.toFixed(1)} bara, ${tSuc.toFixed(0)} °C → ${comp.design.pd.toFixed(1)} bara, ${comp.design.Td.toFixed(0)} °C`, '—', '', cooled ? `suction cooler ${coolerDuty.toFixed(0)} kW` : 'no suction cooler'],
    ['Compressor shaft power at duty', rd(compDuty.p / 1000, 0), rd((clamp(v.compPmaxPct, 50, 300) / 100) * comp.design.power / 1000, 0), 'kW', compDuty.p <= (clamp(v.compPmaxPct, 50, 300) / 100) * comp.design.power ? 'ok' : 'driver limit'],
    ['Compressor speed at duty / after turndown', `${(100 * compDuty.n).toFixed(1)} / ${(100 * compLow.n).toFixed(1)}`, `${v.rotNminPct}–${v.rotNmaxPct}`, '%', compLow.n <= v.rotNminPct / 100 + 1e-6 ? 'at minimum speed' : 'ok'],
    ['Surge margin at duty / minimum in the test', `${compDuty.sm.toFixed(1)} / ${Math.max(comp.minSm, -100).toFixed(1)}`, `≥ ${v.compSmCtl} (control line)`, '%', comp.surgeEvents ? 'SURGE' : comp.minSm >= 0.5 * v.compSmCtl ? 'ok' : 'low'],
    ['Recycle valve after turndown', rd(100 * compLow.r, 1), '100', '%', compLow.r < 0.95 ? 'ok' : 'fully open'],
    ['Surge events in the test', comp.surgeEvents, '0', '', comp.surgeEvents ? 'fail' : 'ok'],
    ['Compressor coast-down to half speed after trip', comp.tCoast50 === null ? `> ${320 - compTrip}` : rd(comp.tCoast50, 1), '—', 's', 'recycle fully open on trip'],
    ['Map point at duty (network suite)', mapPt ? `speed ${(100 * mapPt.speed).toFixed(0)} %, surge margin ${mapPt.surgeMargin.toFixed(0)} %` : '—', '—', '', mapPt ? (mapPt.inside ? 'inside the map' : 'outside the map') : 'not available'],
    ['Pump rated flow / head', `${(pump.q[0] * 3600).toFixed(0)} m³/h / ${pumpHr.toFixed(0)} m`, '—', '', pumpCv ? 'curve as entered' : typeof NET.pumpHead === 'function' ? 'curve from the network suite' : 'built-in parabolic curve'],
    ['Pump speed at duty / after turndown', `${(100 * pumpDuty.n).toFixed(1)} / ${(100 * pumpLow.n).toFixed(1)}`, `${v.rotNminPct}–${v.rotNmaxPct}`, '%', 'ok'],
    ['Pump shaft power at duty', rd(pumpDuty.p / 1000, 0), rd(pump.pRated / 1000, 0), 'kW', pumpDuty.p <= 1.05 * pump.pRated ? 'ok' : 'overload'],
    ['Pump flow after turndown', rd(pumpLow.q * 3600, 1), `≥ ${(pumpMinFlow * 3600).toFixed(0)}`, 'm³/h', pumpFlowMargin >= 0 ? 'ok' : 'below minimum flow'],
    ['Pump flow-control error (IAE)', rd(pump.iae, 3), '—', 's (fraction of rated flow)', 'ok'],
    ['Pump coast-down to half speed after trip', pump.tCoast50 === null ? `> ${pumpEnd - pumpTrip}` : rd(pump.tCoast50, 2), '—', 's', 'ok']],
    note: 'Tests: gas feed reduced to the turndown rate over 10 s at 60 s (anti-surge control), trip at 240 s (recycle opens, rotor coasts down); pump flow set-point stepped to the turndown rate at 30 s and tripped at 80 s. Speeds are relative to the rated speed.' });
  tables.push({ title: 'Chemical inventory, ratio control and solids triggers', columns: ['Item', 'Value', 'Unit'], rows: [['Stock at the start', rd(tank0, 1), 'm³'], ['Continuous use', rd(contUsed, 2), 'm³/d'], ['Use per shutdown event', rd(batchEvent, 1), 'm³'], ['Days of autonomy (continuous injection)', autonomy === null ? '—' : rd(autonomy, 1), 'd'], ['Shutdown treatments in stock', treatments === null ? '—' : rd(treatments, 2), ''], ['Re-order level', rd(reorder, 1), 'm³'], ['Re-order level reached after', stock.tReorder === null ? `> ${invDays}` : rd(stock.tReorder, 1), 'd'], ['Deliveries in the period', stock.deliveries.length, ''], ['Lowest stock', rd(stock.minLevel, 1), 'm³'], ['Inhibitor-to-water ratio at restart', rd(ratioKg, 4), 'kg/kg'], ['Lowest ratio achieved / required', rd(ratioRun.ratioMin / Math.max(ratioKg, 1e-12), 3), ''], ['Saved by ratio control per restart', rd(ratioSaved, 2), 'm³'],
    ['Wax deposit (now / limit)', `${depWax.toFixed(2)} / ${v.maxDepositMm}`, 'mm'], ['Scale deposit (now / limit)', `${depScale.toFixed(2)} / ${v.maxDepositMm}`, 'mm'], ['Hydrate in the liquid (now / limit)', `${(100 * v.hydFrac).toFixed(1)} / ${(100 * v.maxHydFrac).toFixed(0)}`, 'vol %'], ['Bore restriction (now / limit)', `${v.blockagePct.toFixed(1)} / ${v.maxBlockagePct}`, '%'], ['Days until the wax limit', daysToPig === null ? '—' : rd(daysToPig, 1), 'd'], ['Pigging interval used', pigIv > 0 ? rd(pigIv, 1) : '—', 'd'], ['Triggers active', triggers.length ? triggers.map((t) => t.tag).join(', ') : 'none', '']] });
  if (recRows.length) {
    tables.push({ title: 'Operating records against the model: metrics by record type', columns: ['Record type', 'Abscissa', 'Unit', 'Points', 'Bias', 'RMS error', 'MAPE (%)', 'Agreement (%)'], rows: recStats.map((s) => [REC[s.kind].label, REC[s.kind].x, REC[s.kind].unit, s.n, txt(s.bias, 3), txt(s.rmse, 3), txt(s.mape, 1), s.agree === null ? '—' : rd(100 * s.agree, 0)]), note: 'The model values are the predictions of this run (driven by the inputs only); nothing is fitted to the records here. Use the calibration tab to estimate parameters from such records.' });
    tables.push({ title: 'Operating records against the model: rows', columns: ['Record type', 'Abscissa / tag', 'Recorded', 'Model', 'Model − recorded'], rows: recCmp.slice(0, 200).map((r) => [r.known ? REC[r.kind].label : r.kind + ' (unknown type)', r.kind === 'alarm' ? r.tag : rd(r.t, 3), rd(r.value, 4), r.model === null ? '—' : rd(r.model, 4), r.model === null ? '—' : rd(r.model - r.value, 4)]) });
  }

  // KPIs
  const stt = (ok, warn) => (ok ? 'ok' : warn ? 'warn' : 'bad');
  kpis.push({ label: 'Cooldown time', value: rd(cooldownTime, 1), unit: neverCools ? 'h (not reached)' : 'h', status: stt(cooldownTime >= 12 || neverCools, cooldownTime >= 6), help: 'Time until the first point of the system reaches the hydrate temperature plus margin at the falling settle-out pressure.' });
  kpis.push({ label: 'No-touch time', value: rd(noTouch, 1), unit: 'h', status: stt(noTouch >= 4, noTouch >= 1.5), help: 'Cooldown time minus the decision allowance and the time the selected preservation needs.' });
  kpis.push({ label: 'Maximum shutdown', value: rd(maxShutdown, 1), unit: 'h', status: stt(maxShutdown >= v.tShut, maxShutdown >= 0.7 * v.tShut), help: preserved ? 'With the selected preservation in place, limited by the gelled-line restart pressure or the simulated horizon.' : 'Without preservation: the cooldown time.' });
  kpis.push({ label: 'Settle-out pressure', value: rd(so.pSettle, 1), unit: 'bara', status: 'ok' });
  kpis.push({ label: 'Restart pressure', value: rd(restartPressure, 1), unit: 'bara', status: stt(restartPressure <= 0.85 * v.pAvail, restartPressure <= v.pAvail), help: 'Separator pressure + settled liquid legs + gelled-crude yield term, with 5 % allowance.' });
  kpis.push({ label: 'Warm-up time', value: rd(restartTime, 1), unit: 'h', status: 'ok', help: 'Until the arrival temperature stays within 1 °C of its final value.' });
  kpis.push({ label: 'Ramp-up surge', value: rd(ramp.vMax, 1), unit: 'm³', status: stt(ramp.vMax <= surgeAllow, ramp.vMax <= catVol), help: `Allowance ${surgeAllow.toFixed(0)} m³.` });
  kpis.push({ label: 'Blowdown time', value: rd(blowdownTime, 2), unit: bdReached ? 'h' : 'h (target not reached)', status: stt(bdReached && blowdownTime <= Math.max(cooldownTime - v.tDecision, 0.1), bdReached) });
  kpis.push({ label: 'Blowdown minimum temperature', value: rd(Math.min(bdMinT, bdMinTd), 1), unit: '°C', status: stt(Math.min(bdMinTw, bdMinTd) >= v.tMinDesign + 10, Math.min(bdMinTw, bdMinTd) >= v.tMinDesign), help: `Gas ${bdMinT.toFixed(1)} °C, wall ${bdMinTw.toFixed(1)} °C, downstream of the valve ${bdMinTd.toFixed(1)} °C; design minimum ${v.tMinDesign} °C.` });
  kpis.push({ label: 'Pig transit', value: pigTransit === null ? 'stalled' : rd(pigTransit, 2), unit: 'h', status: stt(pigTransit !== null && pig.vMean >= 0.5 && pig.vMean <= 5, pigTransit !== null) });
  kpis.push({ label: 'Pig liquid surge', value: rd(pigSurge, 1), unit: 'm³', status: stt(pigSurge <= surgeAllow, pigSurge <= catVol), help: `Liquid arriving faster than the drain capacity; slug catcher ${catVol.toFixed(0)} m³.` });
  kpis.push({ label: `${inh.name} dose`, value: rd(doseWt, 1), unit: 'wt %', status: stt(doseGov.attainable && doseWt < 50, doseGov.attainable) });
  kpis.push({ label: 'Inhibitor use', value: rd(inhibitorRate, 2), unit: 'm³/d', status: 'ok', help: 'Continuous injection plus shutdown batches averaged over the year.' });
  kpis.push({ label: 'Heating power', value: rd(heatingPower, 0), unit: 'kW', status: stt(heatingPower <= v.heatMaxKw, heatingPower <= 1.2 * v.heatMaxKw), help: 'To hold the heated length at the hydrate temperature plus margin during shutdown.' });
  kpis.push({ label: 'Critical choke opening', value: zCrit === null ? 'stable' : rd(100 * zCrit, 1), unit: '%', status: stt(!unstable, slugSuppressed), help: 'Opening above which the riser limit cycle (severe slugging) starts in open loop.' });
  kpis.push({ label: 'Slugging amplitude under control', value: rd(pidRun.ampClosed, 2), unit: 'bar', status: stt(slugSuppressed, false), help: `Open loop at the same opening: ${pidRun.ampOpen.toFixed(1)} bar peak to peak.` });
  kpis.push({ label: 'Operating window', value: env.feasible ? `${(100 * env.qMin).toFixed(0)}–${(100 * env.qMax).toFixed(0)}` : 'none', unit: '% of case rate', status: stt(env.feasible && rate >= env.qMin && rate <= env.qMax, env.feasible) });
  kpis.push({ label: 'Uptime', value: rd(100 * uptime, 2), unit: '%', status: stt(uptime >= 0.95, uptime >= 0.9) });
  kpis.push({ label: 'Deferred production', value: rd(deferredVolume, 0), unit: 'Sm³/y', status: 'ok' });

  kpis.push({ label: 'Surge capacity', value: rd(catVol, 0), unit: 'm³', status: stt(catVol >= catReq * 0.999, catVol >= designSurge), help: `Required ${catReq.toFixed(0)} m³ for a design surge of ${designSurge.toFixed(1)} m³; peak level ${catcherLevel.toFixed(0)} %.` });
  kpis.push({ label: 'Compressor surge margin (minimum)', value: rd(Math.max(comp.minSm, -100), 1), unit: '%', status: stt(comp.surgeEvents === 0 && comp.minSm >= 0.5 * v.compSmCtl, comp.surgeEvents === 0), help: `Turndown to ${v.compTurndownPct} % of the gas rate; control line ${v.compSmCtl} %; ${comp.surgeEvents} surge event(s).` });
  kpis.push({ label: 'Compressor / pump power', value: `${(compDuty.p / 1000).toFixed(0)} / ${(pumpDuty.p / 1000).toFixed(0)}`, unit: 'kW', status: 'ok', help: 'Shaft power at the operating rate.' });
  kpis.push({ label: 'Chemical stock cover', value: treatments === null ? (autonomy === null ? 'not needed' : rd(autonomy, 0)) : rd(treatments, 1), unit: treatments === null ? (autonomy === null ? '' : 'd') : 'shutdowns', status: stt(stock.runOut === null && (treatments === null || treatments >= 2), stock.runOut === null), help: 'Shutdown treatments (or days of continuous injection) covered by the stock.' });
  kpis.push({ label: 'Blowdown leaves the line hydrate-safe', value: safeBd ? 'yes' : 'no', unit: bothEnds ? 'both ends' : 'topside only', status: stt(safeBd, stratSel.id !== 'blowdown'), help: `Trapped gas keeps ${headBd.toFixed(1)} bar above the vent pressure; hydrate-safe pressure ${pSafe.toFixed(1)} bara.` });

  // warnings and recommendations
  if (pic.source === 'kernel estimate') warnings.push({ level: 'info', msg: 'The flow suite has not published a profile for this case; the kernel steady solution is used as the starting point.' });
  if (!neverCools && v.tShut > cooldownTime && stratSel.id === 'none') warnings.push({ level: 'bad', msg: `The planned shutdown (${v.tShut} h) is longer than the cooldown time (${cooldownTime.toFixed(1)} h) and no preservation is selected: the line enters the hydrate region ${km(coldSpotX).toFixed(1)} km from the inlet.` });
  if (stratSel.id !== 'none' && !stratSel.feasible) warnings.push({ level: 'bad', msg: `The selected preservation (${stratSel.name}) is not feasible for this shutdown${stratBest ? `; ${stratBest.name} is` : ''}.` });
  if (!safeBd && pSafeRaw !== null) warnings.push({ level: stratSel.id === 'blowdown' ? 'warn' : 'info', msg: `${bothEnds ? 'Depressurising from both ends' : 'Topside blowdown alone'} cannot make the seabed hydrate-safe: gas pockets behind the liquid seals keep ${headBd.toFixed(1)} bar above the vent pressure against a hydrate-safe pressure of ${pSafe.toFixed(1)} bara at ${tAmbMin.toFixed(0)} °C.${!bothEnds && pSafe - headBoth > v.pBack * 1.02 ? ` Depressurising from both ends leaves ${headBoth.toFixed(1)} bar and is hydrate-safe.` : ''}` });
  if (cooled && tSuc < coolerHydT + v.hydMargin) warnings.push({ level: 'warn', msg: `The compressor suction cooler outlet (${tSuc.toFixed(0)} °C) is inside the hydrate margin at ${v.sepP} bara (hydrate temperature ${coolerHydT.toFixed(1)} °C): inject inhibitor upstream of the cooler or raise its outlet temperature.` });
  if (comp.surgeEvents > 0) warnings.push({ level: 'bad', msg: `The compressor surges ${comp.surgeEvents} time(s) when the gas rate falls to ${v.compTurndownPct} %: ${v.antiSurge ? 'the recycle valve is too slow or too small' : 'anti-surge control is switched off'}.` });
  else if (comp.minSm < 0.5 * v.compSmCtl && v.antiSurge) warnings.push({ level: 'warn', msg: `Surge margin falls to ${comp.minSm.toFixed(1)} % during the turndown (control line ${v.compSmCtl} %): open the recycle faster or raise the control margin.` });
  if (pumpFlowMargin < 0) warnings.push({ level: 'warn', msg: `At ${v.pumpTurndownPct} % rate the export pump runs at ${(pumpLow.q * 3600).toFixed(0)} m³/h, below its minimum continuous flow of ${(pumpMinFlow * 3600).toFixed(0)} m³/h: open the minimum-flow recycle.` });
  if (stock.runOut !== null) warnings.push({ level: 'bad', msg: `Chemical storage runs empty after ${stock.runOut.toFixed(0)} d (${v.tankVol} m³ tank, ${contUsed.toFixed(1)} m³/d and ${batchEvent.toFixed(0)} m³ per shutdown, ${v.resupplyLead} d resupply lead time).` });
  if (ratioShort > 0.02) warnings.push({ level: 'warn', msg: `During the restart the injection pump cannot hold the inhibitor-to-water ratio: ${(100 * ratioShort).toFixed(0)} % short of the requirement.` });
  for (const t of triggers) warnings.push({ level: 'warn', msg: `${t.tag}: ${t.msg} → ${t.action}.` });
  if (recUnknown.length) warnings.push({ level: 'info', msg: `Operating-record rows of unknown kind were ignored: ${recUnknown.join(', ')}.` });
  if (Math.min(bdMinTw, bdMinTd) < v.tMinDesign) warnings.push({ level: 'bad', msg: `Blowdown temperature ${Math.min(bdMinTw, bdMinTd).toFixed(0)} °C is below the minimum design metal temperature of ${v.tMinDesign} °C.` });
  if (restartPressure > v.pAvail) warnings.push({ level: 'bad', msg: `Restart needs ${restartPressure.toFixed(0)} bara (${so.headUphill.toFixed(0)} bar liquid legs, ${dpGel.toFixed(0)} bar gel) but ${v.pAvail} bara is available.` });
  if (rampReq === null) warnings.push({ level: 'warn', msg: 'No ramp duration up to 16 times the planned one keeps the surge inside the allowance; raise the drain capacity or the surge volume.' });
  if (pig.stalled) warnings.push({ level: 'bad', msg: 'The pig stalls: bypass leakage exceeds the gas velocity at the pigging rate.' });
  else if (pig.vMean > 5 || pig.vMean < 0.5) warnings.push({ level: 'warn', msg: `Mean pig velocity ${pig.vMean.toFixed(1)} m/s is outside the usual 0.5–5 m/s range.` });
  if (pInPigMax > v.pAvail) warnings.push({ level: 'warn', msg: `Inlet pressure while the pig-driven slug climbs the riser reaches ${pInPigMax.toFixed(0)} bara, above the ${v.pAvail} bara available.` });
  if (doseShort > 0.5) warnings.push({ level: 'bad', msg: `The injection setting of ${contUsed.toFixed(1)} m³/d gives ${doseAchieved.toFixed(0)} wt %, ${doseShort.toFixed(0)} wt % short of the requirement.` });
  if (pumpUtil > 1) warnings.push({ level: 'warn', msg: `Continuous injection of ${contUsed.toFixed(1)} m³/d exceeds the pump capacity of ${v.pumpMax} m³/h.` });
  if (sInj > 1 && (steadyNeeded || stratSel.id === 'inhibit')) warnings.push({ level: 'warn', msg: `The first ${(sInj / 1000).toFixed(1)} km upstream of the injection point receive no inhibitor.` });
  if (sensorFailed) warnings.push({ level: 'warn', msg: `Inlet-pressure transmitter failed: the slug loop is in manual at ${(100 * zTarget).toFixed(0)} % opening${zCrit !== null ? ', below the critical opening' : ''}.` });
  if (v.sensorState === 'bias' && Math.abs(v.sensorBias) > 0) warnings.push({ level: 'info', msg: `The inlet-pressure reading carries a bias of ${v.sensorBias} bar: the loop holds the true pressure ${(-v.sensorBias).toFixed(1)} bar away from the set-point.` });
  if (!doseGov.attainable) warnings.push({ level: 'bad', msg: `${inh.name} cannot provide ${dTgov.toFixed(0)} °C of depression; consider depressurisation or heating.` });
  if (unstable && !slugSuppressed) warnings.push({ level: 'bad', msg: `The ${v.ctlMode} loop does not hold ${(100 * zTarget).toFixed(0)} % opening (remaining oscillation ${pidRun.ampClosed.toFixed(1)} bar); operate at ≤ ${(90 * zCrit).toFixed(0)} % or retune.` });
  if (unstable && !manual && pidRun.carryOver - pidRun.carryOverOpen > 0.5) warnings.push({ level: 'warn', msg: `Even under automatic control the separator overfills (${(pidRun.carryOver - pidRun.carryOverOpen).toFixed(0)} m³ carried over): enlarge the surge volume or lower the opening.` });
  if (unstable && pidRun.carryOverOpen > 0.5) warnings.push({ level: manual ? 'bad' : 'info', msg: `${manual ? 'Controller in manual' : 'Comparison case, controller in manual'} at ${(100 * zTarget).toFixed(0)} % opening: ${pidRun.ampOpen.toFixed(0)} bar pressure swings and the slugs overfill the separator (${pidRun.carryOverOpen.toFixed(0)} m³ carried over in the simulated period)${manual ? '' : `; in automatic the level peaks at ${(100 * pidRun.maxLevelClosed).toFixed(0)} %`}.` });
  if (!env.feasible) warnings.push({ level: 'bad', msg: 'No rate in the scanned range satisfies all operating limits: ' + env.text.join('; ') + '.' });
  else if (rate < env.qMin || rate > env.qMax) warnings.push({ level: 'warn', msg: `The operating rate (${v.rateFrac} %) lies outside the window ${(100 * env.qMin).toFixed(0)}–${(100 * env.qMax).toFixed(0)} %.` });
  if (v.severeSlugging) warnings.push({ level: 'info', msg: 'The flow suite reports severe slugging for this case; the control section shows what the choke loop can do about it.' });
  for (const a of raised) warnings.push({ level: a.level === 'trip' ? 'bad' : 'warn', msg: `${a.tag}: ${a.msg}${a.action ? ' → ' + a.action : ''}` });

  if (neverCools) rec.push(`No point reaches the hydrate temperature within ${v.tHorizon} h of shut-in: no preservation is needed for shutdowns up to that duration.`);
  else rec.push(`Act within ${noTouch.toFixed(1)} h of shut-in: the cold spot (${coldComp ? coldComp.name + ', ' : ''}${km(coldSpotX).toFixed(1)} km) reaches the hydrate temperature + ${v.hydMargin} °C after ${cooldownTime.toFixed(1)} h, and ${stratSel.phrase} takes ${stratSel.lead.toFixed(1)} h plus ${v.tDecision} h to decide.`);
  if (stratBest && stratBest.id !== stratSel.id) rec.push(`For a ${v.tShut} h shutdown the lowest-cost feasible preservation is ${stratBest.phrase} (about $${fmt(stratBest.cost, 3)} per event against $${fmt(stratSel.cost, 3)} for the selected one${stratSel.feasible ? '' : ', which is not feasible'}).`);
  if (doseWt > 0) rec.push(`Dose ${inhId} to ${doseWt.toFixed(0)} wt % of the water phase (${dTgov.toFixed(0)} °C depression incl. ${v.inhMargin} °C margin): ${batchVol.toFixed(1)} m³ per shutdown${contRate > 0 ? ` and ${contRate.toFixed(1)} m³/d continuously` : ''}; start injection ${tProtect.toFixed(1)} h before a planned shut-in so the front covers the whole line.`);
  rec.push(rampReq === null ? `Restart at ${v.qStartPct} % with ${restartPressure.toFixed(0)} bara available at the inlet; the liquid surge cannot be kept inside ${surgeAllow.toFixed(0)} m³ by ramping alone.` : `Restart at ${v.qStartPct} % (needs ${restartPressure.toFixed(0)} bara) and ${rampReq > 0.99 * v.rampHours || rampReq > 0.5 ? `ramp at ≤ ${((100 - v.qStartPct) / Math.max(rampReq, 0.05)).toFixed(0)} %/h (${rampReq.toFixed(1)} h to full rate) to keep the liquid surge below ${surgeAllow.toFixed(0)} m³` : `ramp as planned over ${v.rampHours} h: the liquid surge (${ramp.vMax.toFixed(0)} m³) stays below ${surgeAllow.toFixed(0)} m³ even for a fast ramp`}${(tSafe ?? restartTime) > 0.05 ? `; keep ${inhId} on for ${(tSafe ?? restartTime).toFixed(1)} h until the whole line is outside the hydrate region` : ''}.`);
  rec.push(safeBd ? `Depressurise ${bothEnds ? `from both ends (${v.orificeMm} mm topside and ${v.orifice2Mm} mm at the inlet service line)` : `through the ${v.orificeMm} mm topside orifice`} to ${pEndB.toFixed(1)} bara: ${blowdownTime.toFixed(1)} h, peak flare ${bd.peak.toFixed(1)} kg/s, coldest metal ${Math.min(bdMinTw, bdMinTd).toFixed(0)} °C.` : `Do not rely on ${bothEnds ? 'depressurisation' : 'topside blowdown'} for hydrate protection: the liquid seals leave ${seabedPAfter.toFixed(0)} bara in the trapped gas pockets (hydrate-safe below ${pSafe.toFixed(0)} bara); ${!bothEnds && pSafe - headBoth > v.pBack * 1.02 ? 'depressurise from both ends' : 'plan inhibitor displacement or heating'} instead.`);
  if (pigTransit !== null) rec.push(`Pig at ${v.pigRatePct} % rate: ${pigTransit.toFixed(1)} h transit at ${pig.vMean.toFixed(1)} m/s; expect ${(pig.received * volScale).toFixed(0)} m³ of liquid over ${(pig.duration / 60).toFixed(0)} min${pigSurge > surgeAllow ? ` — ${pigSurge.toFixed(0)} m³ more than the drain can take, so lower the pigging rate or pre-drain the slug catcher` : ', inside the slug-catcher allowance'}${pigIv > 0 ? `; every ${pigIv} d (${pigRuns.toFixed(0)} runs a year)` : ''}.`);
  if (sensorFailed) rec.push(`Repair the inlet-pressure transmitter: until then keep the choke at ${(100 * zTarget).toFixed(0)} % in manual (costs ${Math.max(0, eqT.Pp / 1e5 - sm0.steady(Math.min(1, zAuto)).Pp / 1e5).toFixed(0)} bar of back-pressure against controlled operation at ${(100 * zAuto).toFixed(0)} %).`);
  else if (!v.ctlAuto && zCrit !== null) rec.push(`The slug controller is in manual: the choke is held at ${(100 * zTarget).toFixed(0)} % (stable without feedback). In automatic it could run at ${(100 * zAuto).toFixed(0)} % and the inlet pressure would be ${Math.max(0, eqT.Pp / 1e5 - sm0.steady(Math.min(1, zAuto)).Pp / 1e5).toFixed(0)} bar lower.`);
  else if (zCrit !== null) rec.push(slugSuppressed ? `Run the choke at ${(100 * zTarget).toFixed(0)} % under inlet-pressure control (Kc ${(100 * sel.kc).toFixed(1)} %/bar${sel.ti > 0 ? `, Ti ${(sel.ti / 60).toFixed(0)} min` : ''}): the open-loop limit cycle starting at ${(100 * zCrit).toFixed(0)} % (${pidRun.ampOpen.toFixed(0)} bar swings) is suppressed to ${pidRun.ampClosed.toFixed(2)} bar and the inlet pressure is ${(sm0.steady(zCrit).Pp / 1e5 - eqT.Pp / 1e5).toFixed(0)} bar lower than at the largest stable manual opening.` : `Keep the choke at or below ${(90 * zCrit).toFixed(0)} % in manual; the tested controller does not stabilise ${(100 * zTarget).toFixed(0)} %.`);
  else rec.push('The riser is stable at every choke opening for this rate: no slug control is needed.');
  rec.push(env.feasible ? `Keep the rate between ${(100 * env.qMin).toFixed(0)} and ${(100 * env.qMax).toFixed(0)} % of the case rate (${env.text.slice(0, 2).join('; ')}); the best daily margin is at ${(100 * qOpt).toFixed(0)} % with ${wOpt.toFixed(0)} wt % inhibitor, heating ${bestC.heat ? 'on' : 'off'}.` : `No feasible rate window: ${env.text.join('; ')}.`);
  if (heatingPower > 0) rec.push(`Electrical heating would need ${heatingPower.toFixed(0)} kW (${(heatingPower * 24 * v.elecPrice).toFixed(0)} $/d) to hold ${tHold.toFixed(0)} °C during a shutdown${hotOilTime !== null ? `; hot-oil circulation at ${v.hotOilRate} kg/s warms a cold line in ${hotOilTime.toFixed(1)} h` : ''}.`);

  rec.push(`Gas compressor: design ${(comp.design.power / 1000).toFixed(0)} kW at ${(comp.design.md).toFixed(1)} kg/s and a surge margin of ${comp.design.sm.toFixed(0)} %; when the gas rate falls to ${v.compTurndownPct} % the anti-surge loop ${comp.surgeEvents ? 'does NOT prevent surge' : `holds ${compLow.sm.toFixed(0)} % margin with the recycle ${(100 * compLow.r).toFixed(0)} % open at ${(100 * compLow.n).toFixed(0)} % speed`}${cooled ? `; the suction cooler (${tArr0.toFixed(0)} → ${tSuc.toFixed(0)} °C, ${coolerDuty.toFixed(0)} kW) saves ${((compPowerNoCool - compPowerDuty) / 1000).toFixed(0)} kW of shaft power` : ''}. Export pump: ${(pumpDuty.p / 1000).toFixed(0)} kW at ${(100 * pumpDuty.n).toFixed(0)} % speed for ${(pumpDuty.q * 3600).toFixed(0)} m³/h; it coasts to half speed in ${pump.tCoast50 === null ? 'more than ' + (pumpEnd - pumpTrip) : pump.tCoast50.toFixed(1)} s after a trip.`);
  rec.push(`Chemical storage: ${tank0.toFixed(0)} m³ in stock covers ${treatments === null ? 'no shutdown batches (none needed)' : treatments.toFixed(1) + ' shutdown treatment(s)'}${autonomy === null ? '' : ` and ${autonomy.toFixed(0)} d of continuous injection`}; re-order at ${reorder.toFixed(0)} m³${stock.tReorder === null ? ` (not reached in ${invDays} d)` : ` (reached after ${stock.tReorder.toFixed(0)} d)`}. Ratio control of the inhibitor to the water rate during the restart saves ${ratioSaved.toFixed(1)} m³ against a fixed full-rate injection.`);
  if (triggers.length) rec.push(`Solids triggers active: ${triggers.map((t) => t.action.toLowerCase()).join('; ')}.`); else if (daysToPig !== null) rec.push(`At ${v.depRateMmD} mm/d the wax deposit reaches the ${v.maxDepositMm} mm limit in ${daysToPig.toFixed(0)} d: pig at least every ${pigIv} d.`);
  const balances = [
    { name: 'Cooldown energy (stored heat released = heat lost to ambient, J)', in: eBal.drop, out: eBal.lost },
    { name: 'Blowdown mass (initial + liberated = remaining + discharged, kg)', in: bd.mass.initial + bd.flashed, out: bd.mass.final + bd.discharged },
    { name: 'Pig liquid (swept = received + still ahead + leaked, m³)', in: pig.swept, out: pig.received + pig.inPipe + pig.leaked },
    { name: 'Restart energy (inflow = outflow + losses + storage, J)', in: warm.energy.in, out: warm.energy.out + warm.energy.lost + warm.energy.stored },
    { name: 'Chemical inventory (start + delivered = end + used, m³)', in: tank0 + stock.delivered, out: stock.level[stock.level.length - 1] + stock.used },
    { name: 'Settle-out gas mass (before = after, kg)', in: so.gasMass, out: sum(st.T.map((T, i) => fm.at(so.pSettle, T).rhoG * (1 - st.holdup[i]) * st.ds * A)) },
  ];
  const outputs = {
    cooldownTime: rd(cooldownTime, 3), noTouchTime: rd(noTouch, 3), maxShutdown: rd(maxShutdown, 3), coldSpotX: rd(coldSpotX, 1), restartPressure: rd(restartPressure, 2), restartTime: rd(restartTime, 3),
    blowdownTime: rd(blowdownTime, 4), blowdownMinT: rd(Math.min(bdMinT, bdMinTd), 2), blowdownEndP: rd(bdEndP, 3), pigTransit: pigTransit === null ? null : rd(pigTransit, 3), pigSurge: rd(pigSurge, 2), pigDp: rd(pig.dpPig / 1e5, 3),
    inhibitorDose: rd(doseWt, 2), inhibitorRate: rd(inhibitorRate, 3), inhibitorCostPerDay: rd(inhCostDay, 0), heatingPower: rd(heatingPower, 1), uptime: rd(uptime, 5), deferredVolume: rd(deferredVolume, 0),
    envelope: { qMin: env.feasible ? rd(env.qMin, 3) : null, qMax: env.feasible ? rd(env.qMax, 3) : null, limits: env.text.slice() }, controller: { kc: rd(100 * sel.kc, 4), ti: rd(sel.ti, 1), td: rd(sel.td, 2), mode: v.ctlMode, rule: sel.rule, unit: '% opening per bar, s, s' },
    chokeOpening: rd(chokeOut, 1), slugSuppressed: !!slugSuppressed, alarms, eventsPerYear: { shutdowns: rd(shutdownsYr, 2), pigRuns: rd(pigRuns, 1), blowdowns: rd(blowdownsYr, 2) },
    // additional values
    cooldownReached: !neverCools, coldSpotT12: rd(cold12, 3), lumpedColdT12: rd(lumpedCold, 3), settleOutPressure: rd(so.pSettle, 2), watTime: tWat[iCold] === null ? null : rd(tWat[iCold], 2), hydrateSafeRestart: tSafe === null ? null : rd(tSafe, 2), rampTime: rampReq === null ? null : rd(rampReq, 2), rampSurge: rd(ramp.vMax, 2),
    blowdownMinMetalT: rd(Math.min(bdMinTw, bdMinTd), 2), blowdownPeakRate: rd(bd.peak, 3), blowdownSeabedP: rd(seabedPAfter, 2), blowdownReached: bdReached, pigVelocity: rd(pig.vMean, 3), pigPeakInletP: rd(pInPigMax, 2), waxRemoved: rd(waxVol * clamp(v.pigEff, 0, 100) / 100, 3), waxInventoryMax: rd(waxMax, 3), liquidInventory: rd(liqInvSteady, 1), liquidInventoryPeak: rd(Math.max(...ramp.inventory), 1), watMargin: rd(tArr0 - v.wat, 2), hydrateMargin: rd(steadyMargin, 2), doseAchieved: rd(doseAchieved, 2), warmUpTime: rd(restartTime, 3), restartRate: rd(v.qStartPct, 1), rampRate: rampReq === null ? null : rd((100 - v.qStartPct) / Math.max(rampReq, 0.05), 2),
    inhibitor: inhId, inhibitorBatch: rd(batchVol, 2), protectTime: rd(tProtect, 3), criticalChoke: zCrit === null ? null : rd(100 * zCrit, 2), slugAmplitudeOpen: rd(pidRun.ampOpen, 2), slugAmplitudeClosed: rd(pidRun.ampClosed, 3), ctlIae: rd(pidRun.iae, 4), ctlMeanP: rd(pidRun.meanP, 3), mpcIae: rd(mpcRun.iae, 4), ekfError: rd(ekfRel, 5),
    surgeCapacity: rd(catVol, 1), surgeCapacityRequired: rd(catReq, 1), designSurge: rd(designSurge, 2), controllerMode: manual ? 'manual' : 'automatic', sensitivityPeak: rd(loopLin.ms, 3), complementaryPeak: rd(loopLin.mt, 3), ukfError: ukfRel === null ? null : rd(ukfRel, 5), nmpcIae: nm ? rd(nm.iae, 4) : null, nmpcAmplitude: nm ? rd(nm.amp, 3) : null,
    blowdownRoute: bothEnds ? 'both ends' : 'topside', blowdownTrappedHead: rd(headBd, 2), blowdownHydrateSafe: !!safeBd, blowdownRelaxationTime: relaxing && isNum(bd.theta) ? rd(bd.theta, 4) : null,
    pump: { power: rd(pumpDuty.p / 1000, 1), speed: rd(100 * pumpDuty.n, 1), flow: rd(pumpDuty.q * 3600, 1), head: rd(pumpHr, 1), coast50: pump.tCoast50 === null ? null : rd(pump.tCoast50, 2), minFlowMargin: rd(pumpFlowMargin, 1) },
    compressor: { power: rd(compDuty.p / 1000, 1), designPower: rd(comp.design.power / 1000, 1), speed: rd(100 * compDuty.n, 1), surgeMargin: rd(compDuty.sm, 1), surgeMarginMin: rd(Math.max(comp.minSm, -100), 1), surgeEvents: comp.surgeEvents, recycle: rd(100 * compLow.r, 1), dischargeP: rd(comp.design.pd, 1), dischargeT: rd(comp.design.Td, 1), coast50: comp.tCoast50 === null ? null : rd(comp.tCoast50, 2), suctionT: rd(tSuc, 1), coolerDuty: rd(coolerDuty, 1) },
    pumpPower: rd(pumpDuty.p / 1000, 1), compressorPower: rd(compDuty.p / 1000, 1),
    chemical: { stock: rd(tank0, 1), autonomyDays: autonomy === null ? null : rd(autonomy, 1), treatments: treatments === null ? null : rd(treatments, 2), reorderLevel: rd(reorder, 1), reorderDay: stock.tReorder === null ? null : rd(stock.tReorder, 1), runOutDay: stock.runOut === null ? null : rd(stock.runOut, 1), ratioSaved: rd(ratioSaved, 2) },
    solidsTriggers: triggers.map((t) => ({ tag: t.tag, action: t.action })), pigInterval: rd(pigIv, 2), daysToPig: daysToPig === null ? null : rd(daysToPig, 1), pigInLine: pigInLine,
    records: recStats.map((s) => ({ kind: s.kind, n: s.n, bias: rd(s.bias, 4), rmse: rd(s.rmse, 4), mape: isNum(s.mape) ? rd(s.mape, 2) : null })),
    optimum: { rate: rd(qOpt, 3), dose: rd(wOpt, 2), heating: !!bestC.heat, pigging: !!bestC.pigs, marginPerDay: rd(profitOpt, 0) }, strategy: stratSel.id, bestStrategy: stratBest ? stratBest.id : null, heatingCostPerDay: rd(heatingPower * 24 * v.elecPrice, 0), state: safeState,
  };
  const summary = `${neverCools ? `No hydrate risk within ${v.tHorizon} h of shut-in` : `Cooldown time ${cooldownTime.toFixed(1)} h (no-touch ${noTouch.toFixed(1)} h) with the cold spot at ${km(coldSpotX).toFixed(1)} km`}; restart needs ${restartPressure.toFixed(0)} bara and ${restartTime.toFixed(1)} h to warm up; blowdown takes ${blowdownTime.toFixed(1)} h; ${zCrit === null ? 'the riser is stable at any choke opening' : `slugging starts above ${(100 * zCrit).toFixed(0)} % choke opening and is ${slugSuppressed ? 'suppressed' : 'not suppressed'} by the ${v.ctlMode} loop at ${(100 * zTarget).toFixed(0)} %`}; operating window ${env.feasible ? `${(100 * env.qMin).toFixed(0)}–${(100 * env.qMax).toFixed(0)} %` : 'not found'} of the case rate; uptime ${(100 * uptime).toFixed(1)} %.`;
  progress(1, 'Done');
  return { summary, kpis, warnings, recommendations: rec, plots, tables, balances, outputs };
}

// ---- declarations --------------------------------------------------------------------------------------------------
const sel = (key, label, value, options, help) => ({ key, label, type: 'select', value, options: options.map((o) => (Array.isArray(o) ? { value: o[0], label: o[1] } : { value: o, label: o })), help });
// synthetic 48 h historian extract: steady flow, turndown, a 4 h trip, restart and ramp, then 110 % (model + bias + noise)
const LOG_SAMPLE = [{ t: 0, rate: 100, choke: 25, pIn: 107.2, tArr: 41.8 }, { t: 2, rate: 100, choke: 25, pIn: 106.9, tArr: 41.5 }, { t: 4, rate: 100, choke: 25, pIn: 107.2, tArr: 42 }, { t: 6, rate: 100, choke: 25, pIn: 108, tArr: 41.9 }, { t: 8, rate: 100, choke: 25, pIn: 106.8, tArr: 41.9 }, { t: 10, rate: 100, choke: 25, pIn: 107.5, tArr: 42.1 }, { t: 11, rate: 80, choke: 22, pIn: 96.9, tArr: 41.4 }, { t: 12, rate: 80, choke: 22, pIn: 96.9, tArr: 40.2 }, { t: 14, rate: 80, choke: 22, pIn: 97.1, tArr: 39.4 }, { t: 16, rate: 0, choke: 2, pIn: 81.1, tArr: 36.3 }, { t: 17, rate: 0, choke: 2, pIn: 81.1, tArr: 34.5 }, { t: 18, rate: 0, choke: 2, pIn: 80.1, tArr: 32.7 }, { t: 19, rate: 0, choke: 2, pIn: 81.6, tArr: 30.7 }, { t: 20, rate: 30, choke: 12, pIn: 112.9, tArr: 29.2 }, { t: 21, rate: 45, choke: 15, pIn: 98.7, tArr: 29.9 }, { t: 22, rate: 60, choke: 18, pIn: 93.8, tArr: 29 }, { t: 23, rate: 80, choke: 22, pIn: 96.7, tArr: 30.2 }, { t: 24, rate: 100, choke: 25, pIn: 107.5, tArr: 33 }, { t: 26, rate: 100, choke: 25, pIn: 106.7, tArr: 36.3 }, { t: 28, rate: 100, choke: 25, pIn: 106.1, tArr: 37.9 }, { t: 30, rate: 100, choke: 25, pIn: 107.7, tArr: 39.2 }, { t: 33, rate: 100, choke: 25, pIn: 107.1, tArr: 40.9 }, { t: 36, rate: 110, choke: 28, pIn: 111.7, tArr: 41.9 }, { t: 39, rate: 110, choke: 28, pIn: 112.3, tArr: 43.2 }, { t: 42, rate: 110, choke: 28, pIn: 111.7, tArr: 43.8 }, { t: 45, rate: 110, choke: 28, pIn: 112.8, tArr: 43.9 }, { t: 48, rate: 110, choke: 28, pIn: 111.4, tArr: 43.6 }];
// placeholder rows that show the format of the operating-record table (synthetic: generated from this model for the reference case and rounded; not measurements)
const RECORD_SAMPLE = [{ kind: 'commissioning', t: 60, value: 84.05, tag: '' }, { kind: 'commissioning', t: 80, value: 81.33, tag: '' }, { kind: 'commissioning', t: 100, value: 90.82, tag: '' }, { kind: 'commissioning', t: 120, value: 95.82, tag: '' }, { kind: 'startup', t: 2, value: 4.702, tag: '' }, { kind: 'startup', t: 6, value: 4.388, tag: '' }, { kind: 'startup', t: 12, value: 38.56, tag: '' }, { kind: 'startup', t: 24, value: 48.35, tag: '' }, { kind: 'shutdown', t: 1, value: 84.54, tag: '' }, { kind: 'shutdown', t: 6, value: 81.03, tag: '' }, { kind: 'shutdown', t: 12, value: 76.18, tag: '' }, { kind: 'shutdown', t: 24, value: 69.12, tag: '' }, { kind: 'restart', t: 1, value: 13.53, tag: '' }, { kind: 'restart', t: 3, value: 16.44, tag: '' }, { kind: 'restart', t: 6, value: 13.33, tag: '' }, { kind: 'restart', t: 12, value: 42.33, tag: '' }, { kind: 'esd', t: 30, value: 96.03, tag: '' }, { kind: 'esd', t: 120, value: 95.26, tag: '' }, { kind: 'esd', t: 600, value: 94.1, tag: '' }, { kind: 'esd', t: 1800, value: 92.47, tag: '' }, { kind: 'blowdown', t: 0.25, value: 56.08, tag: '' }, { kind: 'blowdown', t: 0.5, value: 46.24, tag: '' }, { kind: 'blowdown', t: 1, value: 26.71, tag: '' }, { kind: 'blowdown', t: 2, value: 12.59, tag: '' }, { kind: 'cooldown', t: 2, value: 39.81, tag: '' }, { kind: 'cooldown', t: 6, value: 29.64, tag: '' }, { kind: 'cooldown', t: 12, value: 19.14, tag: '' }, { kind: 'cooldown', t: 24, value: 9.975, tag: '' }, { kind: 'pigging', t: 1, value: 3.907, tag: '' }, { kind: 'pigging', t: 2, value: 7.886, tag: '' }, { kind: 'pigging', t: 3, value: 11.85, tag: '' }, { kind: 'pigging', t: 4, value: 15.43, tag: '' }, { kind: 'pigarrival', t: 60, value: 8.448, tag: '' }, { kind: 'pigarrival', t: 80, value: 5.041, tag: '' }, { kind: 'pigarrival', t: 100, value: 4.009, tag: '' }, { kind: 'pigarrival', t: 120, value: 3.175, tag: '' }, { kind: 'surge', t: 0.5, value: 48.48, tag: '' }, { kind: 'surge', t: 2, value: 94.07, tag: '' }, { kind: 'surge', t: 4, value: 132.7, tag: '' }, { kind: 'surge', t: 8, value: 175.9, tag: '' }, { kind: 'tracer', t: 3, value: 1.002, tag: '' }, { kind: 'tracer', t: 4, value: 0.9548, tag: '' }, { kind: 'tracer', t: 5, value: 0.9864, tag: '' }, { kind: 'tracer', t: 6, value: 1.028, tag: '' }, { kind: 'concentration', t: 10, value: 0.9929, tag: '' }, { kind: 'concentration', t: 20, value: 2.092, tag: '' }, { kind: 'concentration', t: 40, value: 3.909, tag: '' }, { kind: 'concentration', t: 60, value: 5.947, tag: '' }, { kind: 'valve', t: 2, value: 0.2263, tag: '' }, { kind: 'valve', t: 5, value: 0.4654, tag: '' }, { kind: 'valve', t: 10, value: 0.7122, tag: '' }, { kind: 'valve', t: 20, value: 0.935, tag: '' }, { kind: 'compressor', t: 5, value: 96.58, tag: '' }, { kind: 'compressor', t: 10, value: 88.96, tag: '' }, { kind: 'compressor', t: 20, value: 77.77, tag: '' }, { kind: 'compressor', t: 40, value: 62.09, tag: '' }, { kind: 'pump', t: 1, value: 74.32, tag: '' }, { kind: 'pump', t: 2, value: 65.51, tag: '' }, { kind: 'pump', t: 4, value: 44.85, tag: '' }, { kind: 'pump', t: 8, value: 31.09, tag: '' }, { kind: 'level', t: 1, value: 48.46, tag: '' }, { kind: 'level', t: 3, value: 53.87, tag: '' }, { kind: 'level', t: 5, value: 49.35, tag: '' }, { kind: 'level', t: 7, value: 51.9, tag: '' }, { kind: 'alarm', t: 0, value: 0, tag: 'PAHH-100' }, { kind: 'alarm', t: 0, value: 0, tag: 'LAHH-200' }, { kind: 'alarm', t: 0, value: 0, tag: 'KAL-120' }, { kind: 'alarm', t: 0, value: 0, tag: 'XA-140' }];
const INPUTS = [
  { group: 'Line and facilities', tab: 'inputs', help: 'Geometry and fluid come from the case (network and flow suites when they have been run, otherwise the reference tie-back). These values can be linked from the other suites.', fields: [
    { key: 'rateFrac', label: 'Operating rate', unit: '% of case rate', value: 100, min: 10, max: 200, typical: [40, 110], help: 'Rate at which the line runs before the shutdown and for the control study.' },
    { key: 'idMm', label: 'Inner diameter', unit: 'mm', value: BASE.idMm, min: 25, max: 1500 },
    { key: 'wtMm', label: 'Wall thickness', unit: 'mm', value: BASE.wtMm, min: 1, max: 100 },
    { key: 'uValue', label: 'Overall U-value (0 = use the layers as entered)', unit: 'W/m²K on ID', value: BASE.U, min: 0, max: 200, typical: [1, 10], help: 'The coating conductivity is scaled so that the flowing U-value matches this number.' },
    { key: 'layers', label: 'Coating and insulation layers (bore outwards)', type: 'table', columns: [{ key: 'name', label: 'Layer', type: 'text' }, { key: 't', label: 'Thickness', unit: 'mm' }, { key: 'k', label: 'Conductivity', unit: 'W/m/K' }, { key: 'rho', label: 'Density', unit: 'kg/m³' }, { key: 'cp', label: 'Heat capacity', unit: 'J/kg/K' }], value: [{ name: BASE.insulation.name, t: BASE.insulation.t * 1000, k: BASE.insulation.k, rho: 900, cp: 1700 }] },
    { key: 'thermalMass', label: 'Wall + coating thermal mass (0 = from the layers)', unit: 'kJ/m/K', value: 0, min: 0, max: 5000 },
    { key: 'lineVolume', label: 'Line volume (0 = from geometry)', unit: 'm³', value: 0, min: 0, max: 1e6 },
    { key: 'burial', label: 'Burial depth to pipe centre (0 = exposed)', unit: 'm', value: 0, min: 0, max: 10 },
    { key: 'kSoil', label: 'Soil conductivity', unit: 'W/m/K', value: 1.5, min: 0.2, max: 4, showIf: (v) => v.burial > 0 },
    { key: 'tSeabed', label: 'Seabed temperature', unit: '°C', value: BASE.tSeabed, min: -2, max: 35 },
    { key: 'tSurface', label: 'Sea-surface temperature', unit: '°C', value: BASE.tSeaSurface, min: -2, max: 35 },
    { key: 'tAir', label: 'Air temperature', unit: '°C', value: BASE.tAir, min: -50, max: 55 },
    { key: 'currentSpeed', label: 'Seabed current', unit: 'm/s', value: BASE.currentSpeed, min: 0, max: 3 },
    { key: 'sepP', label: 'Separator pressure', unit: 'bara', value: BASE.separatorP, min: 1.2, max: 200 },
    { key: 'slugCatcherVol', label: 'Slug-catcher / separator volume', unit: 'm³', value: BASE.slugCatcherVol, min: 1, max: 5000 },
    { key: 'catcherUsable', label: 'Usable surge fraction of that volume', unit: '%', value: 30, min: 5, max: 100, help: 'Share of the volume between the normal level and the high-level alarm (50 % to 80 % by default).' },
    { key: 'catcherAuto', label: 'Size the surge capacity to the design surge when the entered volume is too small', type: 'bool', value: true, help: 'The design uses the larger of the entered volume and 1.1 × the largest design surge (ramp-up, pig, hydrodynamic slug) / usable fraction.' },
    { key: 'qDrainM3h', label: 'Liquid handling (drain) capacity', unit: 'm³/h', value: 230, min: 1, max: 20000, help: 'Largest liquid rate the facility can process continuously.' },
    { key: 'pAvail', label: 'Largest inlet pressure available / allowed', unit: 'bara', value: 180, min: 5, max: 1000, help: 'Shut-in pressure the wells can deliver at the flowline inlet, or the design limit if lower.' },
    { key: 'wat', label: 'Wax appearance temperature', unit: '°C', value: 30, min: -20, max: 90 },
    { key: 'watMargin', label: 'Margin above WAT at arrival', unit: '°C', value: 2, min: 0, max: 20 },
    { key: 'pourPoint', label: 'Pour point', unit: '°C', value: 9, min: -40, max: 60 },
    { key: 'yieldStress', label: 'Gel yield stress below the pour point', unit: 'Pa', value: 10, min: 0, max: 2000 },
    { key: 'slugSurge', label: 'Hydrodynamic slug surge (flow suite)', unit: 'm³', value: 0, min: 0, max: 5000 },
    { key: 'severeSlugging', label: 'Severe slugging reported by the flow suite', type: 'bool', value: false },
    { key: 'turndown', label: 'Minimum stable rate from the flow suite (0 = none)', unit: 'fraction', value: 0, min: 0, max: 1 },
    { key: 'plugTime', label: 'Hydrate plugging time from the solids suite (0 = none)', unit: 'h', value: 0, min: 0, max: 1e5 },
  ] },
  { group: 'Solids state and action triggers', tab: 'inputs', help: 'Deposit state from the solids suite (or from inspection) and the accumulation limits that trigger pigging, chemical injection, heating or a rate reduction in the operating logic.', fields: [
    { key: 'depWaxMm', label: 'Wax deposit thickness now', unit: 'mm', value: 0, min: 0, max: 100 },
    { key: 'depScaleMm', label: 'Scale deposit thickness now', unit: 'mm', value: 0, min: 0, max: 100 },
    { key: 'depHydMm', label: 'Hydrate deposit on the wall now', unit: 'mm', value: 0, min: 0, max: 200 },
    { key: 'hydFrac', label: 'Hydrate volume fraction in the liquid', unit: 'fraction', value: 0, min: 0, max: 1 },
    { key: 'blockagePct', label: 'Bore area lost to deposits', unit: '%', value: 0, min: 0, max: 100 },
    { key: 'depRateMmD', label: 'Wax deposition rate', unit: 'mm/d', value: 0, min: 0, max: 50 },
    { key: 'maxDepositMm', label: 'Largest wax or scale deposit allowed (pigging / chemical trigger)', unit: 'mm', value: 4, min: 0.1, max: 100 },
    { key: 'maxHydFrac', label: 'Largest hydrate fraction allowed (inhibitor trigger)', unit: 'fraction', value: 0.1, min: 0.001, max: 1 },
    { key: 'maxBlockagePct', label: 'Largest bore restriction allowed (remediation trigger)', unit: '%', value: 20, min: 1, max: 100 },
  ] },
  { group: 'Shutdown and cooldown', tab: 'inputs', fields: [
    { key: 'tShut', label: 'Planned shutdown duration', unit: 'h', value: 24, min: 0, max: 2000 },
    { key: 'tHorizon', label: 'Cooldown simulated', unit: 'h', value: 48, min: 1, max: 720 },
    { key: 'hydMargin', label: 'Hydrate safety margin', unit: '°C', value: 3, min: 0, max: 15 },
    { key: 'tDecision', label: 'Decision / mobilisation allowance', unit: 'h', value: 2, min: 0, max: 48 },
    { key: 'noTouchMin', label: 'No-touch time required', unit: 'h', value: 4, min: 0, max: 200, help: 'The automatic selection prefers strategies that leave at least this time before anything has to be done.' },
    sel('preserve', 'Preservation selected', 'auto', [['auto', 'Automatic: lowest-cost feasible strategy'], ['none', 'None (restart within the cooldown time)'], ['inhibit', 'Inhibitor bullheading'], ['blowdown', 'Depressurise'], ['heat', 'Electrical heating'], ['hotoil', 'Hot-oil circulation']]),
    { key: 'hInShut', label: 'Inside film coefficient after shut-in', unit: 'W/m²K', value: 120, min: 5, max: 2000, help: 'Natural convection of the settled fluids.' },
    { key: 'uMult', label: 'U-value multiplier (calibration)', unit: '–', value: 1, min: 0.2, max: 5 },
    { key: 'cMult', label: 'Thermal-mass multiplier (calibration)', unit: '–', value: 1, min: 0.2, max: 5 },
    { key: 'components', label: 'Special components (lower insulation, dead legs)', type: 'table', columns: [{ key: 'name', label: 'Component', type: 'text' }, { key: 'x', label: 'Distance', unit: 'm' }, { key: 'u', label: 'U-value', unit: 'W/m²K' }, { key: 'mass', label: 'Thermal mass / line pipe', unit: '–' }, { key: 't0f', label: 'Initial temperature fraction', unit: '–' }],
      value: [{ name: 'Well jumper', x: 0, u: 5, mass: 0.9, t0f: 1 }, { name: 'Manifold valve and connector', x: 50, u: 7, mass: 1.6, t0f: 1 }, { name: 'Chemical-injection dead leg', x: 9000, u: 5, mass: 0.7, t0f: 0.8 }, { name: 'Riser-base spool', x: BASE.riserBaseX, u: 4.5, mass: 1, t0f: 1 }], help: 'Initial temperature fraction: share of the local (fluid − ambient) difference present in the component at shut-in (below 1 for dead legs).' },
  ] },
  { group: 'Restart and ramp-up', tab: 'inputs', fields: [
    { key: 'qStartPct', label: 'Initial restart rate', unit: '% of operating rate', value: 30, min: 5, max: 100 },
    { key: 'rampHours', label: 'Planned ramp duration', unit: 'h', value: 6, min: 0.1, max: 200 },
  ] },
  { group: 'Depressurisation', tab: 'inputs', fields: [
    { key: 'orificeMm', label: 'Blowdown orifice diameter', unit: 'mm', value: 30, min: 1, max: 500 },
    { key: 'cdBlow', label: 'Discharge coefficient', unit: '–', value: 0.85, min: 0.1, max: 1 },
    { key: 'pBack', label: 'Flare back-pressure', unit: 'bara', value: 2, min: 1, max: 50 },
    { key: 'pBlowEnd', label: 'Target pressure (0 = hydrate-safe at ambient)', unit: 'bara', value: 0, min: 0, max: 500 },
    { key: 'tBlowStart', label: 'Valve opens after shut-in', unit: 'h', value: 4, min: 0, max: 500 },
    sel('bdRoute', 'Depressurisation route', 'both', [['both', 'Both ends: topside valve and inlet service line'], ['top', 'Topside valve only']], 'Gas trapped behind liquid seals in the low points keeps the pressure of the liquid legs between it and the nearest open end; venting from both ends halves that.'),
    { key: 'orifice2Mm', label: 'Orifice of the inlet service-line vent', unit: 'mm', value: 30, min: 1, max: 500, showIf: (v) => v.bdRoute === 'both' },
    sel('bdMode', 'Discharge model', 'gas', [['gas', 'Gas only (liquid stays behind)'], ['hem', 'Two-phase, homogeneous equilibrium (ω-method)'], ['hrm', 'Two-phase with relaxation of the gas liberation (homogeneous relaxation)']]),
    { key: 'bdLiquidFrac', label: 'Liquid volume fraction at the valve inlet', unit: '%', value: 5, min: 0, max: 90, showIf: (v) => v.bdMode === 'hem' || v.bdMode === 'hrm' },
    { key: 'bdRelax', label: 'Relaxation time of the gas liberation (0 = Downar-Zapolski correlation)', unit: 's', value: 0, min: 0, max: 1e5, showIf: (v) => v.bdMode === 'hrm' },
    { key: 'bdFlash', label: 'Include gas liberated from the oil', type: 'bool', value: true },
    { key: 'kGas', label: 'Isentropic exponent of the gas', unit: '–', value: 1.28, min: 1.05, max: 1.67 },
    { key: 'hGasWall', label: 'Gas-to-wall film coefficient', unit: 'W/m²K', value: 25, min: 0, max: 1000 },
    { key: 'tMinDesign', label: 'Minimum design metal temperature', unit: '°C', value: -29, min: -196, max: 20 },
  ] },
  { group: 'Pigging', tab: 'inputs', fields: [
    { key: 'pigRatePct', label: 'Rate while pigging', unit: '% of case rate', value: 80, min: 10, max: 150 },
    { key: 'pigFric', label: 'Seal friction coefficient', unit: '–', value: 0.3, min: 0.02, max: 2, help: 'Running differential pressure = 4 μ × contact pressure × seal length / bore.' },
    { key: 'pigSealLen', label: 'Total seal contact length', unit: 'm', value: 0.15, min: 0.01, max: 3 },
    { key: 'pigContact', label: 'Seal contact pressure (interference fit)', unit: 'bar', value: 1.4, min: 0.05, max: 30 },
    { key: 'pigEff', label: 'Wax removed per run', unit: '%', value: 90, min: 0, max: 100 },
    { key: 'pigLaunch', label: 'Launch after reaching steady state (negative = no pig in the sequence)', unit: 'h', value: 2, min: -1, max: 2000 },
    { key: 'pigMass', label: 'Pig mass', unit: 'kg', value: 80, min: 1, max: 5000 },
    { key: 'pigBypass', label: 'Bypass area', unit: '% of bore', value: 2, min: 0, max: 30 },
    { key: 'pigCd', label: 'Bypass discharge coefficient (1/√K)', unit: '–', value: 0.9, min: 0.3, max: 1, help: 'K = 1–1.5 on the bypass velocity is usual for plain bypass pigs (Cd 0.82–1); a deflector disk raises K to about 4 (Cd 0.5).' },
    { key: 'pigX0', label: 'Pig already in the line at (0 = at the launcher)', unit: 'm from the inlet', value: 0, min: 0, max: 1e6 },
    { key: 'pigSlug0', label: 'Liquid already ahead of that pig', unit: 'm³', value: 0, min: 0, max: 1e5, showIf: (v) => v.pigX0 > 0 },
    { key: 'pigLeak', label: 'Liquid leaking back past the pig', unit: '% of swept', value: 3, min: 0, max: 60 },
    { key: 'pigInterval', label: 'Pigging interval (0 = no routine pigging)', unit: 'd', value: 14, min: 0, max: 3650 },
    { key: 'waxByPigging', label: 'Wax managed by pigging (arrival may be below the WAT)', type: 'bool', value: false },
    { key: 'waxThk', label: 'Wax thickness at pigging', unit: 'mm', value: 2, min: 0, max: 50 },
    { key: 'pigCost', label: 'Cost per pig run', unit: '$', value: 15000, min: 0, max: 1e7 },
  ] },
  { group: 'Chemical injection and heating', tab: 'inputs', fields: [
    sel('inhibitor', 'Thermodynamic inhibitor', 'MeOH', INH_OPTS.map((k) => [k, INHIBITORS[k].name])),
    { key: 'leanWt', label: 'Purity of the injected chemical', unit: 'wt %', value: 100, min: 30, max: 100, help: 'About 90 wt % for regenerated (lean) MEG.' },
    sel('dosingBasis', 'Dose governed by', 'shutdown', [['shutdown', 'Shutdown cold spot'], ['steady', 'Steady flow'], ['max', 'Worse of the two (continuous injection)']]),
    { key: 'inhMargin', label: 'Dosing margin', unit: '°C', value: 3, min: 0, max: 15 },
    { key: 'inhRequired', label: 'Minimum dose from the solids suite (0 = none)', unit: 'wt %', value: 0, min: 0, max: 90 },
    { key: 'inhPrice', label: 'Chemical price', unit: '$/m³', value: 1150, min: 0, max: 1e5, help: 'Default: posted North-American methanol reference price of October 2026 (1,450 $/t × 0.792 t/m³).' },
    { key: 'injEff', label: 'Injection efficiency (chemical reaching the water phase)', unit: '%', value: 100, min: 5, max: 100 },
    { key: 'tankVol', label: 'Chemical storage volume', unit: 'm³', value: 250, min: 1, max: 1e5 },
    { key: 'tankLevelPct', label: 'Stock at the start', unit: '% of the tank', value: 80, min: 0, max: 100 },
    { key: 'resupplyLead', label: 'Resupply lead time', unit: 'd', value: 10, min: 0.25, max: 180 },
    { key: 'resupplyVol', label: 'Delivery size', unit: 'm³', value: 150, min: 0, max: 1e5 },
    { key: 'reorderBatches', label: 'Shutdown treatments kept in reserve at the re-order level', unit: '–', value: 1, min: 0, max: 20 },
    { key: 'tauPump', label: 'Injection pump response time', unit: 's', value: 20, min: 1, max: 3600 },
    { key: 'pumpMax', label: 'Injection pump capacity', unit: 'm³/h', value: 20, min: 0.01, max: 5000 },
    { key: 'injRate', label: 'Continuous injection setting (0 = as required)', unit: 'm³/d', value: 0, min: 0, max: 1e5 },
    { key: 'injX', label: 'Injection point', unit: 'm from the inlet', value: 0, min: 0, max: 1e6 },
    { key: 'bullheadVol', label: 'Extra volume for jumpers, trees and manifold', unit: 'm³', value: 8, min: 0, max: 5000 },
    { key: 'dispMult', label: 'Axial dispersion multiplier on the Taylor value', unit: '–', value: 150, min: 1, max: 1e5, help: 'Slug mixing spreads the front far more than single-phase turbulent dispersion.' },
    { key: 'megStorageDays', label: 'Glycol storage', unit: 'd', value: 3, min: 0, max: 60 },
    { key: 'megLossPct', label: 'Glycol make-up (losses)', unit: '% of circulation', value: 1, min: 0, max: 100 },
    { key: 'heatEff', label: 'Heating system efficiency', unit: '%', value: 70, min: 5, max: 100 },
    { key: 'heatLengthPct', label: 'Heated length', unit: '% of line', value: 100, min: 0, max: 100 },
    { key: 'heatInstalled', label: 'Electrical heating installed on the line', type: 'bool', value: false },
    { key: 'heatMaxKw', label: 'Heating power available', unit: 'kW', value: 3000, min: 0, max: 1e6 },
    { key: 'elecPrice', label: 'Electricity price', unit: '$/kWh', value: 0.092, min: 0, max: 5, help: 'Default: average US industrial price of June 2026 (9.17 ¢/kWh); offshore self-generation usually costs more.' },
    { key: 'hotOilInstalled', label: 'Hot-oil circulation loop available', type: 'bool', value: false },
    { key: 'hotOilRate', label: 'Hot-oil circulation rate', unit: 'kg/s', value: 30, min: 0.5, max: 2000 },
    { key: 'hotOilT', label: 'Hot-oil supply temperature', unit: '°C', value: 80, min: 20, max: 200 },
  ] },
  { group: 'Slugging and control', tab: 'setup', help: 'Low-order riser model (gas and liquid mass in the feed pipeline and in the riser) with a topside choke; controllers act on the choke from the inlet pressure.', fields: [
    { key: 'chokeCv', label: 'Topside choke Cv (fully open)', unit: 'US gpm/psi^0.5', value: 400, min: 5, max: 20000 },
    { key: 'chokePct', label: 'Operating choke opening (0 = twice the critical opening)', unit: '%', value: 0, min: 0, max: 100 },
    { key: 'slugControl', label: 'Active slug control available', type: 'bool', value: true, help: 'Used in the envelope: with control the choke may run at twice the critical opening, without it at 90 % of it.' },
    sel('ctlMode', 'Controller', 'PI', ['P', 'PI', 'PID']),
    { key: 'ctlAuto', label: 'Slug controller in automatic', type: 'bool', value: true, help: 'Off: the choke is held in manual at 90 % of the critical opening (stable without feedback).' },
    sel('tuning', 'Tuning', 'auto', [['auto', 'Closed-loop pole search on the linearised model'], ['simc', 'SIMC from the step test (gain-scheduled)'], ['zn', 'Ziegler–Nichols from the step test (gain-scheduled)'], ['robust', 'Robust multi-model PI: sensitivity peak below the limit for the whole uncertainty set'], ['rls', 'Self-tuning: recursive least squares model + SIMC'], ['manual', 'Manual']]),
    sel('adaptive', 'Adaptive element', 'schedule', [['schedule', 'Gain scheduling on the choke opening'], ['none', 'Fixed gain']]),
    { key: 'robGainPct', label: 'Uncertainty set: process gain', unit: '± %', value: 30, min: 0, max: 90 },
    { key: 'robDelayPct', label: 'Uncertainty set: dead time', unit: '+ %', value: 50, min: 0, max: 300 },
    { key: 'msMax', label: 'Largest sensitivity peak Ms allowed for the robust PI', unit: '–', value: 1.6, min: 1.1, max: 4 },
    { key: 'nmpc', label: 'Run the nonlinear MPC comparison', type: 'bool', value: true },
    { key: 'tauValve', label: 'Valve actuator time constant', unit: 's', value: 8, min: 0.5, max: 600 },
    { key: 'actDead', label: 'Actuator dead time', unit: 's', value: 0, min: 0, max: 600 },
    { key: 'tauSensor', label: 'Pressure transmitter time constant', unit: 's', value: 0, min: 0, max: 600 },
    { key: 'chokeExp', label: 'Choke characteristic exponent (flow ∝ opening^n; 1 = linear)', unit: '–', value: 1, min: 0.3, max: 4 },
    { key: 'procK', label: 'Process gain for loop calibration', unit: 'bar per % opening', value: -0.6, min: -100, max: 100, help: 'Used only by the calibration model of recorded step and closed-loop responses.' },
    { key: 'procTau', label: 'Process time constant for loop calibration', unit: 's', value: 600, min: 1, max: 1e6 },
    { key: 'kcMan', label: 'Manual gain', unit: '% opening per bar', value: -15, min: -1000, max: 1000, showIf: (v) => v.tuning === 'manual' },
    { key: 'tiMan', label: 'Manual integral time', unit: 's', value: 1800, min: 1, max: 1e6, showIf: (v) => v.tuning === 'manual' },
    { key: 'tdMan', label: 'Manual derivative time', unit: 's', value: 0, min: 0, max: 1e5, showIf: (v) => v.tuning === 'manual' },
    { key: 'tauCFactor', label: 'SIMC closed-loop time constant / dead time', unit: '–', value: 1, min: 0.2, max: 20 },
    { key: 'deadTime', label: 'Process and measurement dead time', unit: 's', value: 30, min: 0, max: 1800 },
    { key: 'rateLimit', label: 'Choke rate limit', unit: '%/s', value: 0.5, min: 0.001, max: 100 },
    { key: 'zMinPct', label: 'Smallest choke opening', unit: '%', value: 2, min: 0.5, max: 50 },
    { key: 'zMaxPct', label: 'Largest choke opening', unit: '%', value: 100, min: 5, max: 100 },
    { key: 'noiseBar', label: 'Pressure measurement noise (1σ)', unit: 'bar', value: 0.05, min: 0, max: 5 },
    { key: 'pSet', label: 'Inlet-pressure set-point (0 = equilibrium of the operating opening)', unit: 'bara', value: 0, min: 0, max: 1000 },
    { key: 'levelSp', label: 'Separator level set-point', unit: '%', value: 50, min: 10, max: 90 },
    { key: 'levelHi', label: 'High level that overrides the choke', unit: '%', value: 80, min: 20, max: 99 },
    sel('sensorState', 'Inlet-pressure transmitter', 'ok', [['ok', 'Healthy'], ['bias', 'Biased'], ['failed', 'Failed (loop to manual)']]),
    { key: 'sensorBias', label: 'Transmitter bias', unit: 'bar', value: 1, min: -20, max: 20, showIf: (v) => v.sensorState === 'bias' },
    { key: 'spStep', label: 'Set-point step in the test', unit: 'bar', value: -1, min: -20, max: 20 },
    { key: 'cascade', label: 'Cascade: pressure → flow → choke', type: 'bool', value: false },
    { key: 'override', label: 'High separator level overrides the choke', type: 'bool', value: true },
    { key: 'feedForward', label: 'Feed-forward of the liquid inflow to the level valve', type: 'bool', value: false },
    { key: 'tCtl', label: 'Control simulation length', unit: 'h', value: 12, min: 1, max: 72 },
    { key: 'kL', label: 'Low-point liquid orifice coefficient', unit: '–', value: 0.3, min: 0.01, max: 2 },
    { key: 'kH', label: 'Level correction factor', unit: '–', value: 0.7, min: 0.1, max: 1.5 },
    { key: 'topsideLen', label: 'Topside piping upstream of the choke', unit: 'm', value: 150, min: 5, max: 2000 },
    { key: 'mpcNp', label: 'MPC prediction horizon', unit: 'steps', value: 20, min: 3, max: 60 },
    { key: 'mpcNc', label: 'MPC control horizon', unit: 'moves', value: 4, min: 1, max: 12 },
    { key: 'mpcR', label: 'MPC move weight', unit: '–', value: 400, min: 0.01, max: 1e6 },
  ] },
  { group: 'Pump and compressor (maps, limits and control tests)', tab: 'setup', help: 'Liquid export pump on the separator liquid and gas compressor on the separator gas, both variable-speed. Leave the tables empty to use the curves of the network suite (or the built-in fan-law curves).', fields: [
    { key: 'pumpDp', label: 'Export pump differential pressure', unit: 'bar', value: 40, min: 1, max: 600 },
    { key: 'pumpShutoff', label: 'Pump shut-off head / rated head', unit: '–', value: 1.25, min: 1.02, max: 2 },
    { key: 'pumpEta', label: 'Pump efficiency at the rated point', unit: '%', value: 75, min: 20, max: 92 },
    { key: 'pumpStaticPct', label: 'Static part of the system head', unit: '% of rated head', value: 50, min: 0, max: 95 },
    { key: 'pumpJ', label: 'Pump and motor rotor inertia', unit: 'kg·m²', value: 3, min: 0.05, max: 5000 },
    { key: 'pumpRpm', label: 'Pump rated speed', unit: 'rpm', value: 3000, min: 300, max: 20000 },
    { key: 'pumpMinFlowPct', label: 'Minimum continuous flow', unit: '% of rated flow', value: 30, min: 0, max: 90 },
    { key: 'pumpTurndownPct', label: 'Pump control test: flow set-point step to', unit: '% of duty', value: 60, min: 5, max: 100 },
    { key: 'pumpCurve', label: 'Pump curve at rated speed (optional)', type: 'table', columns: [{ key: 'q', label: 'Flow', unit: 'm³/h' }, { key: 'h', label: 'Head', unit: 'm' }], value: [] },
    { key: 'compPd', label: 'Compressor discharge pressure', unit: 'bara', value: 60, min: 2, max: 700 },
    { key: 'compEta', label: 'Polytropic efficiency', unit: '%', value: 78, min: 40, max: 92 },
    { key: 'compJ', label: 'Compressor train rotor inertia', unit: 'kg·m²', value: 40, min: 0.5, max: 50000 },
    { key: 'compRpm', label: 'Compressor rated speed', unit: 'rpm', value: 10000, min: 1000, max: 40000 },
    { key: 'compSmCtl', label: 'Surge control line (margin to surge)', unit: '%', value: 10, min: 2, max: 40 },
    { key: 'antiSurge', label: 'Anti-surge control in service', type: 'bool', value: true },
    { key: 'compTurndownPct', label: 'Compressor control test: gas rate falls to', unit: '% of duty', value: 50, min: 5, max: 100 },
    { key: 'compPmaxPct', label: 'Driver power limit', unit: '% of design power', value: 125, min: 50, max: 300 },
    { key: 'rotN0Pct', label: 'Pump and compressor speed at the start of the tests (initial state)', unit: '% of rated', value: 100, min: 20, max: 110 },
    { key: 'rotNminPct', label: 'Lowest operating speed', unit: '% of rated', value: 70, min: 0, max: 100 },
    { key: 'rotNmaxPct', label: 'Highest operating speed', unit: '% of rated', value: 105, min: 100, max: 120 },
    { key: 'coolerOutT', label: 'Suction cooler outlet temperature (0 = no cooler)', unit: '°C', value: 0, min: 0, max: 150, help: 'Active cooling upstream of the compressor: lowers the suction volume flow, the head and the power; checked against the hydrate temperature.' },
    { key: 'compMap', label: 'Compressor map at 100 % speed, from the surge point (optional)', type: 'table', columns: [{ key: 'q', label: 'Inlet flow', unit: 'm³/s' }, { key: 'h', label: 'Polytropic head', unit: 'kJ/kg' }], value: [] },
  ] },
  { group: 'Alarms, interlocks and availability', tab: 'setup', fields: [
    { key: 'alarms', label: 'Alarm and trip thresholds', type: 'table', columns: [{ key: 'tag', label: 'Tag', type: 'text' }, { key: 'key', label: 'Variable', type: 'text' }, { key: 'type', label: 'low / high', type: 'text' }, { key: 'limit', label: 'Limit' }, { key: 'level', label: 'alarm / trip', type: 'text' }, { key: 'action', label: 'Action', type: 'text' }],
      help: 'Variables: tArr (arrival °C), pIn (inlet bara incl. choke), hydMargin (°C, steady), catcherLevel (%), tBlowMin (°C), noTouch (h), cooldown (h), erosion (ratio), restartP (bara), pPig (bara), slugAmp (bar), doseShort (wt %), pumpUtil (%), watMargin (°C), plugTime (h), deposit (mm wax or scale), hydFrac (vol %), blockage (%), surgeMargin (% minimum in the compressor test), pumpFlowMargin (% above minimum flow), compPower (% of design), autonomy (d), treatments (shutdown treatments in stock).',
      value: [{ tag: 'TAL-102', key: 'tArr', type: 'low', limit: 33, level: 'alarm', action: 'Raise rate or start wax inhibitor' }, { tag: 'PAH-100', key: 'pIn', type: 'high', limit: 150, level: 'alarm', action: 'Open choke / check for restriction' }, { tag: 'PAHH-100', key: 'pIn', type: 'high', limit: 185, level: 'trip', action: 'Shut in wells (ESD level 2)' },
        { tag: 'TDAL-110', key: 'hydMargin', type: 'low', limit: 3, level: 'alarm', action: 'Start hydrate inhibitor' }, { tag: 'LAH-200', key: 'catcherLevel', type: 'high', limit: 80, level: 'alarm', action: 'Slow the ramp / reduce pig speed' }, { tag: 'LAHH-200', key: 'catcherLevel', type: 'high', limit: 95, level: 'trip', action: 'Close inlet ESD valve' },
        { tag: 'TALL-301', key: 'tBlowMin', type: 'low', limit: -29, level: 'alarm', action: 'Throttle the blowdown' }, { tag: 'KAL-120', key: 'noTouch', type: 'low', limit: 4, level: 'alarm', action: 'Preserve immediately on shutdown' }, { tag: 'PAH-130', key: 'restartP', type: 'high', limit: 180, level: 'alarm', action: 'Displace / heat before restart' }, { tag: 'XA-140', key: 'slugAmp', type: 'high', limit: 5, level: 'alarm', action: 'Close choke to the stable opening' }, { tag: 'AAL-210', key: 'doseShort', type: 'high', limit: 1, level: 'alarm', action: 'Raise the inhibitor injection rate' },
        { tag: 'UAL-160', key: 'surgeMargin', type: 'low', limit: 3, level: 'alarm', action: 'Check the anti-surge valve and its tuning' }, { tag: 'FAL-170', key: 'pumpFlowMargin', type: 'low', limit: 0, level: 'alarm', action: 'Open the pump minimum-flow recycle' }, { tag: 'LAL-220', key: 'treatments', type: 'low', limit: 1, level: 'alarm', action: 'Order chemical: less than one shutdown treatment in stock' }, { tag: 'LAL-221', key: 'autonomy', type: 'low', limit: 7, level: 'alarm', action: 'Order chemical: less than a week of injection in stock' }] },
    { key: 'plannedPerYear', label: 'Planned shutdowns', unit: '1/y', value: 2, min: 0, max: 100 },
    { key: 'plannedHours', label: 'Planned shutdown duration', unit: 'h', value: 48, min: 0, max: 2000 },
    { key: 'unplannedPerYear', label: 'Unplanned shutdowns', unit: '1/y', value: 6, min: 0, max: 365 },
    { key: 'unplannedHours', label: 'Mean unplanned shutdown duration', unit: 'h', value: 10, min: 0, max: 2000 },
    { key: 'remediationHours', label: 'Extra time when a trip outlasts the cooldown time', unit: 'h', value: 12, min: 0, max: 2000 },
    { key: 'oilPrice', label: 'Oil price', unit: '$/bbl', value: BASE.oilPrice, min: 0, max: 500 },
    { key: 'qLoPct', label: 'Envelope scan: lowest rate', unit: '% of case rate', value: 25, min: 5, max: 95 },
    { key: 'qHiPct', label: 'Envelope scan: highest rate', unit: '% of case rate', value: 150, min: 100, max: 400 },
  ] },
  { group: 'Operating log (historical replay)', tab: 'inputs', help: 'Rows from a SCADA / historian export. The model is driven by the logged rate and choke opening and compared with the logged pressures and temperatures.', fields: [
    { key: 'log', label: 'Operating log', type: 'table', columns: [{ key: 't', label: 'Time', unit: 'h' }, { key: 'rate', label: 'Rate', unit: '% of case rate' }, { key: 'choke', label: 'Choke opening', unit: '%' }, { key: 'pIn', label: 'Inlet pressure', unit: 'bara' }, { key: 'tArr', label: 'Arrival temperature', unit: '°C' }], value: LOG_SAMPLE },
    { key: 'useResidual', label: 'Fit a residual correction on the log', type: 'bool', value: true },
  ] },
  { group: 'Operating records (commissioning, start-up, shutdown, ESD, blowdown, pigging, tracer, transients, alarms)', tab: 'inputs', help: 'One row per recorded value. Kind: commissioning (t = rate %, value = inlet bara), startup / restart (t = h, arrival °C), shutdown (h since shut-in, line bara), esd (s since the trip, inlet bara), blowdown (h, bara), cooldown (h, cold-spot °C), pigging (h since launch, pig position km), pigarrival (rate %, transit h), surge (h since restart, liquid arriving m³/h), tracer (h since injection, C/C0), concentration (injection m³/d, wt % at arrival), valve (s since a step command, fraction of travel), compressor / pump (s since a trip, % of the speed at trip), level (h, separator %), alarm (tag, 1 = raised in the field). The model prediction of this run is listed beside every row with metrics by kind.', fields: [
    { key: 'records', label: 'Operating records', type: 'table', columns: [{ key: 'kind', label: 'Kind', type: 'text' }, { key: 't', label: 'Abscissa' }, { key: 'value', label: 'Recorded value' }, { key: 'tag', label: 'Tag (alarms)', type: 'text' }], value: RECORD_SAMPLE },
  ] },
  { group: 'Numerical resolution', tab: 'mesh', fields: [
    { key: 'nStations', label: 'Axial stations', value: 30, min: 8, max: 120 },
    { key: 'nr', label: 'Radial cells through wall and coatings', value: 12, min: 3, max: 80 },
    { key: 'ntCool', label: 'Cooldown time steps', value: 144, min: 12, max: 2000 },
    { key: 'ntBlow', label: 'Blowdown time steps to the target', value: 300, min: 30, max: 6000 },
    { key: 'dtCtl', label: 'Control sample / integration step', unit: 's', value: 30, min: 2, max: 120 },
    { key: 'nEnv', label: 'Rates in the envelope scan', value: 7, min: 4, max: 16 },
    { key: 'nxInh', label: 'Cells for the inhibitor front', value: 120, min: 20, max: 600 },
  ] },
];
const FIELDS = INPUTS.flatMap((g) => g.fields), DEFAULTS = Object.fromEntries(FIELDS.map((f) => [f.key, f.value])), NUMERIC = FIELDS.filter((f) => !f.type || f.type === 'number').map((f) => f.key);

const CAL_MEMO = new Map();
/** Small memo for the parts of the calibration model: a part is recomputed only when one of the values it depends on changes. */
function memo(tag, key, fn) { const k = tag + '|' + key.join(','); if (CAL_MEMO.has(k)) return CAL_MEMO.get(k); if (CAL_MEMO.size > 4000) CAL_MEMO.clear(); const r = fn(); CAL_MEMO.set(k, r); return r; }
/**
 * Fast calibration model (reference-case fluid constants): cold-spot temperature after a shutdown, line pressure during a blowdown, pig arrival time and
 * seal differential pressure, valve travel after a step command, open-loop and closed-loop pressure responses, transmitter reading, choke pressure drop,
 * pump head and compressor head at a flow and speed, tracer concentration at the outlet, inhibitor concentration at arrival and restart pressure.
 */
function calModel(v0) {
  const v = { ...DEFAULTS, ...v0 }; for (const k of NUMERIC) v[k] = num(v[k], DEFAULTS[k]);
  const id = v.idMm / 1000, wt = v.wtMm / 1000, A = (Math.PI * id * id) / 4, line = memo('line', [], () => caseLine({})), L = line.length, tS = Math.max(num(v0.tShutH, 8), 0.01), t0 = num(v0.t0, 50), tBm = Math.max(num(v0.tBlowMin, 30), 0.01), p0 = Math.max(num(v0.p0, 80), v.pBack + 1), vGas = num(v0.vGas, 2.5);
  const c = memo('cool', [id, wt, v.uValue, v.uMult, v.cMult, v.thermalMass, v.burial, v.kSoil, v.hInShut, v.tSeabed, tS, t0, JSON.stringify(v.layers)], () => { const wall = wallLayers(v, id, wt, 800, 800, []); return cooldown({ ri: id / 2, layers: wall.layers, nr: 8, hIn: v.hInShut, hOut: 800, T0: t0, tAmb: v.tSeabed, cFluid: A * (0.45 * 800 * 2300 + 0.55 * 60 * 2600), dt: (tS * HOUR) / 48, nSteps: 48 }); });
  const b = memo('blow', [id, L, v.lineVolume, p0, v.pBack, v.orificeMm, v.cdBlow, v.kGas, tBm], () => blowdown({ V: (v.lineVolume > 0 ? v.lineVolume : A * L) * 0.55, P0: p0 * 1e5, T0: 290, pBack: v.pBack * 1e5, area: (Math.PI * (v.orificeMm / 1000) ** 2) / 4, cd: v.cdBlow, k: v.kGas, mw: 0.02, Z: 0.85, mode: 'isothermal', dt: (tBm * 60) / 60, n: 60, maxFactor: 1, pEnd: v.pBack * 1.001e5 }));
  const dpSeal = Math.max((4 * v.pigFric * v.pigContact * 1e5 * v.pigSealLen) / id, 100), leak = v.pigCd * (v.pigBypass / 100) * Math.sqrt((2 * dpSeal) / 450);
  // valve and loop records: actuator = dead time + first-order lag; process = gain, lag and dead time; transmitter = first-order lag and bias
  const tV = Math.max(num(v0.tValveS, 10), 0), tauV = Math.max(v.tauValve, 0.05), valvePos = tV <= v.actDead ? 0 : 1 - Math.exp(-(tV - v.actDead) / tauV);
  const tSt = Math.max(num(v0.tStepS, 600), 0), lags = [Math.max(v.procTau, 0.05), tauV, Math.max(v.tauSensor, 1e-3)], tEff = tSt - v.deadTime - v.actDead;
  let stepPv = 0; // unit step through three first-order lags in series (explicit sub-steps of the fastest lag)
  if (tEff > 0) { const n = 400, h = tEff / n, x = [0, 0, 0]; for (let k = 0; k < n; k++) { let u = 1; for (let j = 0; j < 3; j++) { x[j] = (x[j] + (h / lags[j]) * u) / (1 + h / lags[j]); u = x[j]; } } stepPv = v.procK * x[2]; }
  const tCl = Math.max(num(v0.tClS, 600), 1), clPv = memo('cl', [v.procK, v.procTau, v.deadTime, v.actDead, v.kcMan, v.tiMan, tCl], () => { const cl = pidLoop({ K: v.procK, tau: Math.max(v.procTau, 0.05), theta: v.deadTime + v.actDead, kc: v.kcMan, ti: v.tiMan, td: 0, dt: tCl / 300, tEnd: tCl, sp: 1 }); return cl.y[cl.y.length - 1]; });
  const zCh = clamp(num(v0.zChokePct, 30), 1, 100) / 100, wCh = Math.max(num(v0.wChoke, 40), 0), dpChoke = (wCh / (2.403e-5 * v.chokeCv * zCh ** v.chokeExp)) ** 2 / num(v0.rhoChoke, 300) / 1e5;
  const nP = clamp(num(v0.nPumpPct, 100), 5, 130) / 100, pumpHeadM = ownPumpHead(num(v0.qPumpM3h, 150) / 3600, { qr: num(v0.qPumpRated, 200) / 3600, hr: (v.pumpDp * 1e5) / (800 * G), shutoff: v.pumpShutoff, speed: nP });
  const nC = clamp(num(v0.nCompPct, 100), 30, 130) / 100, qdC = num(v0.qCompDesign, 0.25), hdC = (0.9 * R * 313) / 0.021 / ((v.kGas - 1) / (v.kGas * (v.compEta / 100))) * ((v.compPd / v.sepP) ** ((v.kGas - 1) / (v.kGas * (v.compEta / 100))) - 1), phi = clamp(num(v0.qCompM3s, 0.25) / (nC * qdC), 0.62, 1.3), compHeadKj = (nC * nC * hdC * (1 + 0.22 * (1 - phi * phi) - 2.2 * Math.max(phi - 1.1, 0) ** 2)) / 1000;
  const uL = Math.max(num(v0.uLiq, 1.2), 0.01), tracerC = frontAnalytic(L, Math.max(num(v0.tTracerH, 4), 1e-6) * HOUR, uL, Math.max(10.1 * (id / 2) * uL * Math.sqrt(0.02 / 8) * v.dispMult, 1e-4));
  const leanC = clamp(v.leanWt, 30, 100) / 100, rhoLeanC = leanC * 792 + (1 - leanC) * 1000, mI = (Math.max(num(v0.qInjM3d, 30), 0) / DAY) * rhoLeanC * leanC * (clamp(v.injEff, 5, 100) / 100), doseMeas = (100 * mI) / Math.max(num(v0.mWaterKgs, 2), 1e-6) / (1 + mI / leanC / Math.max(num(v0.mWaterKgs, 2), 1e-6));
  return { tCold: c.Tf[48], pBlow: b.pFinal / 1e5, pigArrival: L / Math.max(vGas - leak, 0.05) / HOUR, pigDp: dpSeal / 1e5, valvePos, stepPv, clPv, pvRead: num(v0.pTrue, 100) + v.sensorBias, dpChoke, pumpHeadM, compHeadKj, tracerC, doseMeas, restartP: (v.sepP + num(v0.headBar, 7) + (4 * v.yieldStress * Math.max(num(v0.gelLenM, 0), 0)) / id / 1e5) * 1.05 };
}

function verify() {
  const out = [], add = (name, expected, got, tol, note) => out.push({ name, expected, got, tol, pass: Number.isFinite(got) && Math.abs(got - expected) <= tol, note });
  const steel = { t: 0.004, k: 45, rho: 7850, cp: 470 }, ri = 0.1;
  { // 1 lumped limit: thin conductive wall, weak outside film
    const hOut = 12, cF = 4e4, c = cooldown({ ri, layers: [steel], nr: 4, hIn: 5000, hOut, T0: 60, tAmb: 5, cFluid: cF, dt: 60, nSteps: 600, theta: 0.5 }), Cw = steel.rho * steel.cp * Math.PI * ((ri + steel.t) ** 2 - ri * ri), tau = (cF + Cw) / (hOut * 2 * Math.PI * (ri + steel.t));
    add('Lumped cooldown against the analytic exponential', 5 + 55 * Math.exp(-36000 / tau), c.Tf[600], 0.05, 'Thin steel wall, Biot number ≪ 1: T = Ta + (T0 − Ta) exp(−t/τ), τ = ΣC / (h 2π ro)');
    add('Energy conservation during cooldown (Crank–Nicolson)', c.energy.stored0 - c.energy.stored, c.energy.lost, 1e-6 * c.energy.stored0, 'Stored heat released equals heat passed to the ambient, J/m');
  }
  const lay = [{ t: 0.0159, k: 45, rho: 7850, cp: 470 }, { t: 0.08, k: 0.17, rho: 900, cp: 1700 }], r0 = 0.127, r1 = r0 + 0.0159, r2 = r1 + 0.08;
  { // 2 steady multilayer cylinder with internal heating
    const q = 200, Rt = 1 / (300 * 2 * Math.PI * r0) + Math.log(r1 / r0) / (2 * Math.PI * 45) + Math.log(r2 / r1) / (2 * Math.PI * 0.17) + 1 / (600 * 2 * Math.PI * r2), c = cooldown({ ri: r0, layers: lay, nr: 14, hIn: 300, hOut: 600, T0: 4, tAmb: 4, cFluid: 5e4, dt: 7200, nSteps: 400, qHeat: q, init: 'uniform' });
    add('Radial conduction against the steady multilayer cylinder', 4 + q * Rt, c.Tf[400], 1e-3, 'Fluid temperature under a constant heat input: ΔT = q Σ ln(r₂/r₁)/(2πk) + films');
    add('U-value of the radial grid against the kernel resistance formula', uValue({ id: 2 * r0, wt: 0.0159, kWall: 45, layers: [{ t: 0.08, k: 0.17 }], hIn: 300, hOut: 600 }).U, c.U, 1e-9, 'W/m²K referred to the bore');
    const cs = cooldown({ ri: r0, layers: lay, nr: 12, hIn: 300, hOut: 600, T0: 50, tAmb: 4, cFluid: 5e4, dt: 3600, nSteps: 50, qHeat: 46 / Rt });
    add('Steady-to-transient initialisation holds the steady state', 50, cs.Tf[50], 1e-8, 'Starting from the steady radial profile with the steady heat input, nothing moves');
    const a = cooldown({ ri: r0, layers: lay, nr: 12, T0: 50, tAmb: 4, cFluid: 5e4, dt: 600, nSteps: 72 }), b = cooldown({ ri: r0, layers: lay, nr: 12, T0: 50, tAmb: 4, cFluid: 5e4, dt: 150, nSteps: 288 }), b2 = cooldown({ ri: r0, layers: lay, nr: 12, T0: 50, tAmb: 4, cFluid: 5e4, dt: 37.5, nSteps: 1152 });
    add('Time-step independence of the cooldown solver', b2.Tf[1152], a.Tf[72], 0.1, 'Temperature after 12 h with 10 min steps against 37.5 s steps');
    add('First-order convergence in time (error ratio of the two coarser steps)', (600 - 37.5) / (150 - 37.5), (a.Tf[72] - b2.Tf[1152]) / (b.Tf[288] - b2.Tf[1152]), 0.3, 'Backward Euler: the error falls in proportion to the step (600, 150 and 37.5 s)');
    add('Energy conservation during cooldown (backward Euler)', a.energy.stored0 - a.energy.stored, a.energy.lost, 1e-6 * a.energy.stored0, 'J/m');
  }
  const k = 1.3, mw = 0.02, T = 300, V = 500, area = 7e-4, cd = 0.85, cStar = cd * area * Math.sqrt((k * R * T) / mw) * (2 / (k + 1)) ** ((k + 1) / (2 * (k - 1)));
  { // blowdown
    const b = blowdown({ V, P0: 70e5, T0: T, pBack: 1e5, area, cd, k, mw, mode: 'isothermal', pEnd: 10e5, n: 600 });
    add('Isothermal choked blowdown against the exponential decay', (V / cStar) * Math.log(7), b.tEnd, 1, 'P = P0 exp(−t/τ), τ = V / (Cd A √(kRT/M) (2/(k+1))^((k+1)/(2(k−1)))), s to fall from 70 to 10 bara');
    const a = blowdown({ V, P0: 70e5, T0: T, pBack: 1e5, area, cd, k, mw, mode: 'adiabatic', pEnd: 10e5, n: 600 });
    add('Isentropic temperature of an ideal gas during adiabatic blowdown', T * (a.pFinal / 70e5) ** ((k - 1) / k), a.T[a.T.length - 1], 0.02, 'T/T0 = (P/P0)^((k−1)/k), K');
    add('Adiabatic choked blowdown time against the closed-form solution', ((7 ** ((k - 1) / (2 * k)) - 1) * 2 * V) / ((k - 1) * cStar), a.tEnd, 1, 'P = P0 [1 + (k−1) t / (2τ)]^(−2k/(k−1)), s');
    const h = blowdown({ V, P0: 70e5, T0: T, pBack: 1e5, area, cd, k, mw, mode: 'wall', pEnd: 3e5, n: 300, wallC: 5e8, wallUA: 5e4, extUA: 2e4, tAmb: 277, flash: () => 2e-4, hem: { mLiq: 2e5, rhoL: 800, frac: 0.08 } });
    add('Depressurisation mass conservation (two-phase discharge with gas liberation)', h.mass.initial + h.flashed, h.mass.final + h.discharged, 1e-6 * h.mass.initial, 'Initial + liberated = remaining + discharged, kg');
    add('Critical pressure ratio of the ω-method for ω = 1', 0.6065, omegaFlux(50e5, 50, 1, 1e5).etaC, 2e-3, 'Leung: η_c = 0.6065 for an isothermal ideal gas');
    add('Sub-critical orifice flux tends to the choked flux at the critical ratio', orificeFlux(50e5, 300, 0, { k, mw }), orificeFlux(50e5, 300, 50e5 * criticalRatio(k) * 1.000001, { k, mw }), 1, 'kg/s/m²');
  }
  { // advection–dispersion
    const U = 1.5, D = 4, f = inhibitorFront({ U, L: 10000, D, n: 500, tEnd: 4000, rows: 8 }), i = f.t.length - 1, j = f.xc.findIndex((x) => x >= U * f.t[i] - 100);
    add('Advection–dispersion against the Ogata–Banks (erfc) solution', frontAnalytic(f.xc[j], f.t[i], U, D), f.c[i][j], 5e-3, `C/C0 at x = ${f.xc[j].toFixed(0)} m, t = ${f.t[i].toFixed(0)} s`);
    add('Inhibitor mass balance before breakthrough', U * f.t[i], sum(f.c[i]) * (10000 / 500), 1e-6 * U * f.t[i], '∫C dx = U C0 t (per unit area)');
    const d = inhibitorDose({ dT: 18, inh: 'MEG', S: 20, mWater: 5, lean: 90 });
    add('Mixing with a lean glycol reaches the required aqueous concentration', d.wt, (100 * d.mAq) / (5 + d.mAq / 0.9), 1e-9, 'wt % in the rich phase after dilution by the water in the lean stream');
    add('Dose gives the required hydrate depression (Nielsen–Bucklin)', 18, hydrateDepression({ S: 20, inhWt: d.wt, inh: INHIBITORS.MEG }), 1e-5, '°C');
  }
  { // pigging
    const n = 101, s = Array.from({ length: n }, (_, i) => i * 100), flat = (x) => s.map(() => x), p = pigRun({ s, z: flat(0), vm: flat(2), holdup: flat(0.3), vsl: flat(0), rhoG: flat(60), rhoM: flat(280), rhoL: 800, D: 0.254, bypass: 0, fric: 1e5, qDrain: 0.01 });
    add('Pig transit time at constant velocity', 10000 / 2, p.transit, 1e-6, 'L / v, s');
    add('Liquid displaced equals the swept inventory', 0.3 * ((Math.PI * 0.254 ** 2) / 4) * 10000, p.received + p.inPipe, 1e-6, 'm³, no leakage, stationary liquid');
    add('Pig force balance on an incline', 1e5 + (80 * G * 0.1) / ((Math.PI * 0.254 ** 2) / 4), pigRun({ s, z: s.map((x) => 0.1 * x / Math.sqrt(1 + 0) * 1), vm: flat(2), holdup: flat(0.1), vsl: flat(0), rhoG: flat(60), rhoM: flat(100), rhoL: 800, D: 0.254, bypass: 0, fric: 1e5 }).dpPig, 5, 'Δp = friction + m g sinθ / A, Pa (sinθ ≈ 0.1)');
    const pb = pigRun({ s, z: flat(0), vm: flat(2), holdup: flat(0.3), vsl: flat(0), rhoG: flat(60), rhoM: flat(280), rhoL: 800, D: 0.254, bypass: 0.02, cdBypass: 0.7, fric: 1e5 });
    add('Pig slip from bypass leakage (orifice equation)', 2 - 0.7 * 0.02 * Math.sqrt(2e5 / 280), pb.vMean, 1e-9, 'v = v_gas − Cd (A_b/A) √(2Δp/ρ), m/s');
  }
  { // control
    const pl = pidLoop({ K: 2, tau: 10, theta: 0, kc: 1.5, ti: 0, dt: 0.002, tEnd: 10 });
    add('P-only step response against the first-order closed loop', ((2 * 1.5) / 4) * (1 - Math.exp(-4)), pl.y[pl.y.length - 1], 2e-4, 'y = KKc/(1+KKc) (1 − exp(−(1+KKc) t/τ)) at t = τ');
    const pi = pidLoop({ K: 2, tau: 10, theta: 0, kc: 2.5, ti: 10, dt: 0.002, tEnd: 2 });
    add('PI with Ti = τ gives a first-order closed loop', 1 - Math.exp(-1), pi.y[pi.y.length - 1], 5e-4, 'τc = τ/(K Kc) = 2 s: y(τc) = 1 − e⁻¹');
    const r = tuningRules({ K: 2, tau: 10, theta: 1 });
    add('SIMC PI gain (hand calculation)', 2.5, r.rows[1].kc, 1e-12, 'Kc = τ / (K (τc + θ)) with τc = θ: 10 / (2 × 2)');
    add('SIMC PI integral time (hand calculation)', 8, r.rows[1].ti, 1e-12, 'Ti = min(τ, 4 (τc + θ)) = min(10, 8)');
    const zi = tuningRules({ K: 1, tau: 1e6, theta: 2 });
    add('Ultimate period of an integrating process with dead time', 8, zi.pu, 1e-3, 'Pu = 4θ when τ ≫ θ');
    const aw = pidLoop({ K: 0.5, tau: 5, theta: 0, kc: 2, ti: 4, dt: 0.01, tEnd: 400, sp: 1, uMin: 0, uMax: 1 }), noAw = pidLoop({ K: 0.5, tau: 5, theta: 0, kc: 2, ti: 4, dt: 0.01, tEnd: 400, sp: 1, uMin: 0, uMax: 1, antiWindup: false });
    add('Anti-windup holds the integrator at its analytic bound', 1 - (2 * 0.5) / 2, aw.maxI, 2e-3, 'Back-calculation with Tt = Ti/2 at saturation: I = u_max − Kc e / 2, e = 0.5');
    add('Without anti-windup the integrator winds up (e Kc t / Ti)', (2 * 0.5 * 400) / 4, noAw.maxI, 2, 'Reference case showing what the anti-windup removes');
    const tt = Array.from({ length: 240 }, (_, i) => i * 0.5), fo = identifyFOPDT(tt, tt.map((t) => (t < 3 ? 1 : 1 + 4 * (1 - Math.exp(-(t - 3) / 12)))), 2);
    add('Step-test identification recovers the time constant', 12, fo.tau, 0.1, 'Synthetic first-order-plus-dead-time response, τ = 12 s');
    add('Step-test identification recovers the dead time', 3, fo.theta, 0.1, 'θ = 3 s');
    const m = loopMargins({ K: 1, tau: 1e4, theta: 1 }, { kc: 1e4 * 0.5, ti: 0 });
    add('Phase margin of an integrator with dead time', 90 - (0.5 * 180) / Math.PI, m.pm, 0.3, 'L = k e^(−θs)/s with kθ = 0.5: PM = 90° − kθ (rad)');
    const y = Array.from({ length: 300 }, (_, i) => [Math.sin(i)]), kf = kalman({ A: [[1]], C: [[1]], Q: [[0.01]], R: [[1]], x0: [0], P0: [[1]] }, null, y), Pm = (0.01 + Math.sqrt(1e-4 + 0.04)) / 2;
    add('Kalman gain on a scalar random walk against the algebraic Riccati solution', Pm / (Pm + 1), kf.K[0][0], 1e-9, 'P⁻ = (q + √(q² + 4qr))/2, K = P⁻/(P⁻ + r)');
    const Am = [[1, 0.1], [0, 0.9]], Bm = [[0.005], [0.1]], N = 12, mp = mpc({ A: Am, B: Bm, C: [1, 0], x: [1, 0.5], uPrev: 0.2, np: N, nc: N, q: 1, rDu: 0.1 }), lq = lqFinite([[1, 0.1, 0.005], [0, 0.9, 0.1], [0, 0, 1]], [[0.005], [0.1], [1]], [[1, 0, 0], [0, 0, 0], [0, 0, 0]], [[0.1]], N);
    add('MPC equals the finite-horizon LQ solution when no constraint is active', -(lq.K0[0][0] * 1 + lq.K0[0][1] * 0.5 + lq.K0[0][2] * 0.2), mp.seq[0], 1e-8, 'First move from the QP against the Riccati recursion on the augmented state');
    const mc = mpc({ A: Am, B: Bm, C: [1, 0], x: [1, 0.5], uPrev: 0.2, np: N, nc: N, q: 1, rDu: 0.1, duMax: 0.5 });
    add('MPC constraint handling: the first move sits on the rate limit', -0.5, mc.seq[0], 1e-9, 'Unconstrained move is −2.87; |Δu| ≤ 0.5');
    const ek = ekf({ F: (x) => [x[0]], h: (x) => [x[0]], x0: [0], P0: [[1]], Q: [[0.01]], R: [[1]], y }), kl = kalman({ A: [[1]], C: [[1]], Q: [[0.01]], R: [[1]], x0: [0], P0: [[1]] }, null, y);
    add('Extended Kalman filter reduces to the linear filter for a linear model', kl.x[299][0], ek.x[299][0], 1e-6, 'Same estimate after 300 updates');
    const e2 = eig([[0, 1], [-2, -3]]).map((e) => e[0]).sort((a, b) => a - b);
    add('Eigenvalues of a companion matrix', -2, e2[0], 1e-9, 'Poles of 1/(s² + 3s + 2): −1 and −2');
    const tf = ss2tf([[0, 1], [-2, -3]], [[0], [1]], [1, 0]);
    add('State space to transfer function', 2, tf.den[2] / tf.num[1], 1e-12, 'G(s) = 1/(s² + 3s + 2): denominator constant / numerator constant');
    add('Zero-order-hold discretisation against the scalar solution', Math.exp(-0.5), c2d([[-1]], [[1]], 0.5).Ad[0][0], 1e-12, 'Ad = exp(a Ts)');
    const sm = slugModel({ D: 0.2, Lp: 4000, Lr: 250, theta: 0.02, wG: 0.6, wL: 12, Ps: 20e5, Cv: 200, rhoGnom: 25 }), eq = sm.steady(0.2), res = sm.alg(eq.y, 0.2).d;
    add('Slugging model: mass balances vanish at the computed equilibrium', 0, Math.max(...res.map(Math.abs)), 1e-7, 'Largest residual of the four mass balances, kg/s');
    const a = integrateStiff(sm.f(0.2), eq.y.map((v2, i) => v2 * (i === 1 ? 1.0005 : 1)), 0, 600, { rtol: 1e-7, atol: 1e-6, hInit: 0.5 }), b = rk45(sm.f(0.2), eq.y.map((v2, i) => v2 * (i === 1 ? 1.0005 : 1)), 0, 600, { rtol: 1e-9, atol: 1e-9, maxSteps: 200000 });
    add('Stiff Rosenbrock integrator against adaptive Runge–Kutta 4(5)', sm.alg(b.y[b.y.length - 1], 0.2).Pp / 1e5, sm.alg(a.y[a.y.length - 1], 0.2).Pp / 1e5, 2e-3, 'Inlet pressure of the slugging model 10 min after a perturbation, bara');
    const zc = sm.critical();
    if (zc !== null && zc < 0.98) { const g1 = sm.growth(zc * 0.9), g2 = sm.growth(Math.min(1, zc * 1.1)); add('Hopf point separates decaying and growing oscillations', 1, g1 < 0 && g2 > 0 ? 1 : 0, 0, `Critical opening ${(100 * zc).toFixed(1)} %: real part of the slowest pole changes sign`); }
  }
  { // optimisation
    add('Simplex on a textbook LP', 36, simplex([3, 5], [[1, 0], [0, 2], [3, 2]], [4, 12, 18]).obj, 1e-9, 'max 3x + 5y, x ≤ 4, 2y ≤ 12, 3x + 2y ≤ 18 → (2, 6)');
    add('Two-phase simplex with a ≥ constraint', -9, simplex([-2, -3], [[-1, -1], [1, 0], [0, 1]], [-4, 3, 5]).obj, 1e-9, 'min 2x + 3y, x + y ≥ 4, x ≤ 3 → (3, 1)');
    add('Branch and bound on a small integer programme', 20, branchBound([5, 4], [[6, 4], [1, 2]], [24, 6], [0, 1]).obj, 1e-9, 'max 5x + 4y, 6x + 4y ≤ 24, x + 2y ≤ 6, integers → (4, 0); the relaxation gives 21');
    add('Interior-point method on a constrained quadratic', 2, interiorPoint((x) => (x[0] - 2) ** 2 + (x[1] - 2) ** 2, [(x) => 2 - x[0] - x[1]], [0.5, 0.5], { lo: [0, 0], hi: [3, 3] }).f, 1e-3, 'min (x−2)² + (y−2)², x + y ≤ 2 → (1, 1), f = 2');
    const sph = (x) => (x[0] - 1) ** 2 + (x[1] + 0.5) ** 2;
    add('Genetic algorithm finds a known minimum (seeded)', 1, geneticAlgorithm(sph, [-3, -3], [3, 3]).x[0], 0.01, 'Minimum of (x−1)² + (y+0.5)² at x = 1');
    add('Particle swarm finds a known minimum (seeded)', -0.5, particleSwarm(sph, [-3, -3], [3, 3]).x[1], 0.01, 'y = −0.5');
    add('Bayesian optimisation approaches a known minimum (seeded)', 0, bayesOpt(sph, [-3, -3], [3, 3], { iters: 25 }).f, 0.05, 'Gaussian-process surrogate with expected improvement, 31 evaluations');
    const gp = gaussianProcess([[0], [0.5], [1]], [1, 3, 2]);
    add('Gaussian process interpolates its training data', 3, gp([0.5]).mean, 1e-4, 'Posterior mean at a training point');
    const rs = responseSurface([0, 1, 2, 3, 4], [1, 3, 11, 31, 69], 3);
    add('Response surface reproduces a cubic exactly', 1 + 2.5 + 2.5 ** 3, rs.predict(2.5), 1e-6, 'y = 1 + x + x³ at x = 2.5');
    const env = operatingEnvelope({ q: [0.2, 0.6, 1, 1.4], a: [10, 30, 50, 70], b: [0.2, 0.6, 1.0, 1.4] }, [{ key: 'a', name: 'A', type: 'min', limit: 25 }, { key: 'b', name: 'B', type: 'max', limit: 1.2 }]);
    add('Operating-envelope lower boundary', 0.5, env.qMin, 1e-9, 'a = 50 q ≥ 25');
    add('Operating-envelope upper boundary', 1.2, env.qMax, 1e-9, 'b = q ≤ 1.2');
    const rp = rampSurge({ nodes: [[0, 1], [10, 1]], inv: () => 100, qLiq: () => 0.05, qDrain: 0.03, dt: 10, tEnd: 1000 });
    add('Surge accumulation at a constant excess rate', (0.05 - 0.03) * 1000, rp.V[rp.V.length - 1], 1e-9, 'Volume above the drain capacity = (q − q_drain) t, m³');
  }
  { // operating logic
    let vol = 0; const ev = eventScheduler([{ t: 250.5, tag: 'b', set: { q: 0 } }, { t: 100, tag: 'a', set: { q: 2 } }], { tEnd: 400, dt: 30, state: { q: 1 }, onStep: (a, b, s) => { vol += s.q * (b - a); } });
    add('Event scheduler fires at the exact event time', 250.5, ev.fired[1].t, 0, 'Events at 100 s and 250.5 s with a 30 s march');
    add('Mass balance across switching events', 100 + 2 * 150.5, vol, 1e-9, 'Tank filled at 1, then 2, then 0 m³/s: the steps are cut at the events');
    const rules = [{ tag: 'A', key: 'x', type: 'high', limit: 10, level: 'alarm' }, { tag: 'B', key: 'x', type: 'high', limit: 20, level: 'trip' }, { tag: 'C', key: 'y', type: 'low', limit: 5, level: 'alarm' }], count = (x, y) => evaluateAlarms({ x, y }, rules).length;
    add('Alarm logic truth table', 1213, count(5, 9) * 10000 + count(15, 9) * 1000 + count(25, 9) * 100 + count(5, 1) * 10 + count(25, 1), 0, 'Raised entries for (x, y) = (5, 9), (15, 9), (25, 9), (5, 1), (25, 1): 0, 1, 2, 1, 3 written as digits');
    add('Interlock voting two out of three', 1, (vote([true, false, true], 2) ? 1 : 0) * (vote([true, false, false], 2) ? 0 : 1), 0, 'Trips on two healthy signals, not on one');
    add('Missing measurement raises no alarm (fault handling)', 0, evaluateAlarms({ x: NaN }, rules).length, 0, 'A non-finite value is reported as not evaluated instead of tripping');
    const fsm = stateMachine(['shutdown', 'inhibit', 'restart', 'rampDone', 'restart']);
    add('State machine: planned sequence returns to PRODUCING with one rejected event', 1, fsm.state === 'PRODUCING' && fsm.rejected === 1 ? 1 : 0, 0, 'shutdown → inhibit → restart → rampDone; a second restart is not permitted');
    add('Fail-safe: a trip leads to FAILSAFE from every operating state', Object.keys(OPS_STATES).length - 1, Object.keys(OPS_STATES).filter((s) => s !== 'FAILSAFE' && stateMachine(['trip'], s).state === 'FAILSAFE').length, 0, 'And only a reset leaves it');
    add('Valve logic: FAILSAFE accepts only a reset', 1, stateMachine(['restart', 'blowdown', 'reset'], 'FAILSAFE').state === 'SHUT_IN' && stateMachine(['restart'], 'FAILSAFE').rejected === 1 ? 1 : 0, 0, 'Restart is rejected until the trip has been reset');
  }
  { // restart
    const n = 40, w = warmUp({ s: Array.from({ length: n }, (_, i) => i * 100 + 50), ds: 100, tAmb: new Array(n).fill(4), T0: new Array(n).fill(4), mdot: () => 30, cp: 2500, cFluid: new Array(n).fill(6e4), cWall: new Array(n).fill(1.2e5), gIn: 900, gOut: 2.4, tIn: 60, dt: 60, nSteps: 2400 });
    add('Restart energy conservation', w.energy.in, w.energy.out + w.energy.lost + w.energy.stored, 1e-9 * w.energy.in, 'Enthalpy in = enthalpy out + losses + storage in fluid and wall, J');
    add('Warm-up tends to the steady exponential profile', 4 + 56 * Math.exp((-n * 100) / ((30 * 2500) * (1 / 900 + 1 / 2.4))), w.Tout[w.Tout.length - 1], 0.3, 'T = Ta + (Tin − Ta) exp(−U′L / (ṁ cp)) after 40 h (first-order upwind in x)');
    const fmv = { at: (P, T2) => ({ rhoG: (P * 1e5 * 0.02) / (R * (T2 + KEL)), rhoL: 800 }) }, sv = { s: [50, 150, 250, 350], z: [0, -5, -5, 0], dz: [-5, 0, 0, 5], ds: 100, A: 0.05, P: [60, 60, 40, 40], T: [20, 20, 20, 20], holdup: [0.5, 0.5, 0.5, 0.5] };
    add('Settle-out pressure of an isothermal ideal gas', 50, settleOut(sv, fmv).pSettle, 1e-4, 'Equal gas volumes at 60 and 40 bara equalise at 50 bara');
  }
  { // tuning rules, loop shaping and robust tuning (published values)
    const r = tuningRules({ K: 2, tau: 10, theta: 1.5 }), f = 1 + 0.5 / 9;
    add('SIMC PID for a first-order process: series settings converted to the ideal form', (10 / (2 * 2.25)) * f, r.rows[2].kc, 1e-12, 'Series Kc = τ/(K (τc + θ)) with τc = θ/2, τI = min(τ, 4 (τc + θ)) = 9, τD = θ/3 = 0.5; ideal Kc = Kc (1 + τD/τI)');
    add('Ziegler–Nichols PID ratios', 0.6 * 8 + 1 / 2 + 1 / 8, r.rows[5].kc / r.ku * 8 + r.rows[5].ti / r.pu + r.rows[5].td / r.pu, 1e-9, 'Kc = 0.6 Ku, Ti = Pu/2, Td = Pu/8 (ideal form) written as 8 Kc/Ku + Ti/Pu + Td/Pu');
    const a = loopAnalysis({ k: 1, lags: [4], delay: 1 }, { kc: 2, ti: 4 }), b = loopAnalysis({ k: 1, ints: 1, delay: 1 }, { kc: 0.5, ti: 8 });
    add('Sensitivity peak of the SIMC PI loop on a first-order process with delay', REF.simcMargins.fopdt.ms, a.ms, 0.01, 'Published: Ms = 1.59 (τc = θ)');
    add('Gain margin of the same loop', Math.PI, a.gm, 0.01, 'L = 0.5 e^(−θs)/(θs): GM = π/2 / 0.5 = 3.14 (published 3.14)');
    add('Phase margin of the same loop', REF.simcMargins.fopdt.pm, a.pm, 0.1, 'Published 61.4°');
    add('Complementary sensitivity peak of the same loop', REF.simcMargins.fopdt.mt, a.mt, 0.01, 'Published 1.00');
    add('Sensitivity peak of the SIMC PI loop on an integrating process with delay', REF.simcMargins.integrating.ms, b.ms, 0.01, 'Published 1.70');
    add('Complementary sensitivity peak, integrating process', REF.simcMargins.integrating.mt, b.mt, 0.01, 'Published 1.30');
    add('Gain and phase margin, integrating process', REF.simcMargins.integrating.gm + REF.simcMargins.integrating.pm, b.gm + b.pm, 0.1, 'Published 2.96 and 46.9°');
    const g = plantResponse({ k: 1, lags: [1], delay: 1 }, 1);
    add('Frequency response of e^(−s)/(s + 1) at ω = 1', Math.SQRT1_2, Math.hypot(g[0], g[1]), 1e-12, '|G| = 1/√2');
    const set = [{ k: 1, lags: [4], delay: 1 }, { k: 1.3, lags: [4], delay: 1.5 }, { k: 0.7, lags: [4], delay: 1.5 }], rb = robustPI(set, { msMax: 1.6 }), nomOnly = robustPI([set[0]], { msMax: 1.6 });
    add('Robust PI: the worst-case sensitivity peak sits on the limit', 1.6, rb.worst, 0.03, 'Largest integral gain with Ms ≤ 1.6 for every model of the set (gain ±30 %, dead time +50 %)');
    add('Robust PI is detuned against the nominal-only design', 1, rb.ki < nomOnly.ki ? 1 : 0, 0, 'The uncertainty set costs integral gain');
  }
  { // optimisation: exact box QP and sequential quadratic programming
    const x = boxQPExact([[1, 0], [0, 1]], [-1, -2], [0, 0], [1.5, 1.5]);
    add('Active-set box QP against the analytic solution', 1 + 1.5 * 10, x[0] + 10 * x[1], 1e-10, 'min ½(x² + y²) − x − 2y, 0 ≤ x, y ≤ 1.5 → (1, 1.5)');
    const h71 = HS.HS71, s71 = sqp(h71.f, h71.x0, { eq: h71.eq, ineq: h71.ineq, lo: h71.lo, hi: h71.hi });
    add('SQP on Hock–Schittkowski problem 71', 17.0140173, s71.f, 1e-4, 'Published best value 17.0140173 (equality, inequality and bounds)');
    const h43 = HS.HS43, s43 = sqp(h43.f, h43.x0, { ineq: h43.ineq });
    add('SQP on the Rosen–Suzuki problem (Hock–Schittkowski 43)', -44, s43.f, 1e-5, 'Published best value −44 at (0, 1, 2, −1)');
    add('SQP solution of problem 43 against the published minimiser', 2, s43.x[2], 1e-3, 'x₃ = 2');
    const ros = sqp((p) => 100 * (p[1] - p[0] ** 2) ** 2 + (1 - p[0]) ** 2, [-1.2, 1], { maxIter: 200 });
    add('SQP (BFGS) on the unconstrained Rosenbrock function', 1, ros.x[0], 1e-4, 'Minimum at (1, 1)');
  }
  { // estimation and adaptation
    const y = Array.from({ length: 200 }, (_, i) => [Math.sin(0.3 * i)]), kl = kalman({ A: [[1]], C: [[1]], Q: [[0.01]], R: [[1]], x0: [0], P0: [[1]] }, null, y), uk = ukf({ F: (x) => [x[0]], h: (x) => [x[0]], x0: [0], P0: [[1]], Q: [[0.01]], R: [[1]], y });
    add('Unscented Kalman filter reduces to the linear filter for a linear model', kl.x[199][0], uk.x[199][0], 1e-9, 'Same estimate after 200 updates');
    const A2 = [[1, 0.1], [0, 0.95]], y2 = Array.from({ length: 80 }, (_, i) => [Math.cos(0.2 * i), 0.5 * Math.sin(0.1 * i)]), k2 = kalman({ A: A2, C: [[1, 0], [0, 1]], Q: [[0.01, 0], [0, 0.02]], R: [[0.5, 0], [0, 0.3]], x0: [0, 0], P0: [[1, 0], [0, 1]] }, null, y2), u2 = ukf({ F: (x) => mv(A2, x), h: (x) => x.slice(), x0: [0, 0], P0: [[1, 0], [0, 1]], Q: [[0.01, 0], [0, 0.02]], R: [[0.5, 0], [0, 0.3]], y: y2 });
    add('Unscented filter covariance equals the Riccati covariance (two states)', k2.P[0][1], u2.P[0][1], 1e-9, 'Off-diagonal element after 80 updates');
    const u = Array.from({ length: 60 }, (_, i) => (i % 14 < 7 ? 1 : -0.5)), ys = [0]; for (let k = 1; k < 60; k++) ys.push(0.9 * ys[k - 1] + 0.2 * u[k - 1]);
    const rl = rlsFirstOrder(u, ys, 2);
    add('Recursive least squares recovers the gain of a first-order process', 2, rl.K, 1e-6, 'y_k = 0.9 y_(k−1) + 0.2 u_(k−1): K = b/(1 − a)');
    add('Recursive least squares recovers the time constant', -2 / Math.log(0.9), rl.tau, 1e-4, 'τ = −Δt / ln a');
    const sq = (t) => (Math.floor(t / 60) % 2 ? 0 : 1), Kt = (t) => (t < 300 ? 1 : 4), ad = selfTuningLoop({ K: Kt, tau: 20, theta: 2, dt: 0.5, tEnd: 900, sp: sq, adapt: true }), fx = selfTuningLoop({ K: Kt, tau: 20, theta: 2, dt: 0.5, tEnd: 900, sp: sq, adapt: false });
    add('Self-tuning controller identifies the changed process gain', 4, ad.kHat[ad.kHat.length - 1], 0.2, 'The gain steps from 1 to 4 at t = 300 s');
    add('Self-tuning controller retunes: gain reduced by the gain ratio', 0.25, ad.kcEnd / (20 / (1 * 4)), 0.02, 'SIMC Kc = τ/(K (τc + θ)) follows 1/K');
    add('Adaptation lowers the control error after the gain change', 1, ad.iae < fx.iae ? 1 : 0, 0, `Integral absolute error ${ad.iae.toFixed(1)} against ${fx.iae.toFixed(1)} for the fixed controller`);
    const rc = ratioControl({ wild: (t) => (t < 100 ? 2 : 3), ratio: 0.4, tau: 10, kc: 0.5, ti: 30, dt: 0.5, tEnd: 600 });
    add('Ratio control holds the ratio after a change of the wild flow', 0.4, rc.ratioEnd, 1e-4, 'Injected / wild = 0.4 at steady state');
    const rs = ratioControl({ wild: () => 3, ratio: 0.4, tau: 10, uMax: 1, dt: 0.5, tEnd: 600, q0: 1 });
    add('Ratio control at the pump limit: shortfall equals the missing flow', (1.2 - 1) * 600, rs.shortfall, 1e-6, '(ratio × wild − capacity) × time');
  }
  { // rotating equipment: pump and compressor control tests
    const pp = { qr: 0.06, hr: 300, rho: 800, J: 4, rpm: 3000 }, wr = (2 * Math.PI * 3000) / 60, tqR = (800 * G * 0.06 * 300) / 0.75 / wr, p = pumpSim({ ...pp, tTrip: 0, tEnd: 12, dt: 0.002 });
    add('Pump coast-down against the analytic solution', (4 * wr) / tqR, p.tCoast50, 0.01 * ((4 * wr) / tqR), 'Torque ∝ N² on a friction system: N = N0/(1 + t/tc), tc = J ω / τ_rated; half speed at t = tc');
    add('Affinity laws: flow on a pure friction system follows the speed', 0.8 * 0.06, p.flowAt(0.8), 1e-12, 'Q ∝ N when the system curve is k Q²');
    const p2 = pumpSim({ ...pp, hStatic: 150, qSp: (t) => (t < 10 ? 0.06 : 0.04), tEnd: 60, dt: 0.01 });
    add('Pump flow control reaches the new set-point', 0.04, p2.qEnd, 1e-5, 'Variable-speed PI on the driver torque, static head 50 % of the rated head');
    add('Pump speed at the reduced flow from the curve and the system', Math.sqrt((150 + ((300 - 150) / 0.06 ** 2) * 0.04 ** 2 + ((300 * 0.25) / 0.06 ** 2) * 0.04 ** 2) / (300 * 1.25)), p2.nEnd, 1e-4, 'h_r N² s − h_r (s − 1) (Q/Q_r)² = h_static + k Q²');
    const cb = { qd: 0.5, hd: 110e3, ps: 25e5, T1: 313, Z: 0.9, k: 1.28, mw: 0.021, J: 40, rpm: 10000 }, own = { curve: compressorCurve({ qd: 0.5, hd: 110e3, table: [0.62, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3].map((f) => ({ q: 0.5 * f, h: 110e3 * (1 + 0.22 * (1 - f * f) - 2.2 * Math.max(f - 1.1, 0) ** 2) })) }) };
    add('Surge margin definition', 38, surgeMargin(0.5, 1, own.curve.qSurge), 1e-9, '100 (Q − Q_surge)/Q with the surge point at 62 % of the design flow');
    add('Fan laws on the compressor map', 0.64 * own.curve.head(0.5, 1), own.curve.head(0.4, 0.8), 1e-9 * 110e3, 'H(0.8 Q, 0.8 N) = 0.64 H(Q, N)');
    const c0 = compressorSim({ ...cb, ...own, tEnd: 60, dt: 0.05 }), rhoS = (25e5 * 0.021) / (0.9 * R * 313);
    add('Compressor steady shaft power at the design point', ((rhoS * 0.5 * 110e3) / 0.78) * 1.02, c0.powerEnd, 2e-3 * c0.powerEnd, 'ṁ H_p / η plus 2 % windage; the simulation stays at the design point');
    const md = c0.design.md, feed = (t) => md * (t < 30 ? 1 : 0.5), cA = compressorSim({ ...cb, ...own, tEnd: 260, dt: 0.05, feed }), cN = compressorSim({ ...cb, ...own, tEnd: 260, dt: 0.05, feed, antiSurge: false });
    add('Anti-surge control: no surge when the feed halves', 0, cA.surgeEvents, 0, 'Recycle opens on the surge control line');
    add('Anti-surge control holds the margin on the control line', 10, cA.sm[cA.sm.length - 1], 0.3, 'Integral action: steady margin = control line (10 %)');
    add('Without anti-surge control the same test surges', 1, cN.surgeEvents > 0 ? 1 : 0, 0, `${cN.surgeEvents} surge events`);
    add('Recycle makes up the flow to the control line', 1, Math.abs(cA.q[cA.q.length - 1] / (own.curve.qSurge * cA.nEnd) - 1 / 0.9) < 0.01 ? 1 : 0, 0, 'Q = Q_surge(N) / (1 − 0.10)');
  }
  { // relaxation, inventory, pig in the line, partitioning, seals
    let E = 0; for (let k = 0; k < 50; k++) E = relaxStep(E, 2, 5, 0.2).E;
    add('Relaxation closure against the analytic first-order lag', 2 * 5 * (1 - Math.exp(-10 / 5)), E, 1e-9, 'Held-back mass E = d Θ (1 − e^(−t/Θ)) for a constant demand d');
    add('Relaxation time of the low-pressure correlation (hand calculation)', 6.51e-4 * 10 ** 0.257 * 10 ** 2.24, hrmRelaxationTime(0.1, 5e5, 4.5e5), 1e-9, 'Θ = 6.51e-4 α^-0.257 ψ^-2.24 with α = 0.1, ψ = 0.1');
    const bb = { V: 500, P0: 70e5, T0: 300, pBack: 1e5, area: 7e-4, cd: 0.85, k: 1.3, mw: 0.02, mode: 'wall', pEnd: 5e5, n: 300, wallC: 5e8, wallUA: 5e4, extUA: 2e4, tAmb: 277, flash: () => 2e-4 }, eqb = blowdown(bb), fast = blowdown({ ...bb, relax: 1e-6 }), slow = blowdown({ ...bb, relax: 1e12 }), mid = blowdown({ ...bb, relax: 600 });
    add('Relaxation model tends to the equilibrium model for a short relaxation time', eqb.flashed, fast.flashed, 0.02 * eqb.flashed, 'Gas liberated during the blowdown, kg (the relaxation path takes the demand of each step explicitly: first-order difference to the Heun equilibrium path)');
    add('Relaxation model tends to frozen composition for a long relaxation time', 0, slow.flashed, 0.01, 'No gas is liberated, kg');
    add('Mass conservation with relaxation (liberated + held back = equilibrium demand)', mid.mass.initial + mid.flashed, mid.mass.final + mid.discharged, 1e-6 * mid.mass.initial, 'kg');
    const iv = chemicalInventory({ V: 250, level0: 200, use: 8, batches: [{ t: 5, v: 30 }], reorder: 60, lead: 4, delivery: 100, tEnd: 30, dt: 0.25 });
    add('Chemical inventory: days of autonomy', 25, iv.autonomy, 1e-12, 'Stock / continuous use = 200 / 8');
    add('Chemical inventory: re-order level reached at the expected time', (200 - 30 - 60) / 8, iv.tReorder, 0.25, '(stock − batch − re-order level) / use, within one time step');
    add('Chemical inventory balance', 200 + iv.delivered, iv.level[iv.level.length - 1] + iv.used, 1e-9, 'Start + delivered = end + used, m³');
    const n = 101, s = Array.from({ length: n }, (_, i) => i * 100), flat = (x) => s.map(() => x), pg = pigRun({ s, z: flat(0), vm: flat(2), holdup: flat(0.3), vsl: flat(0), rhoG: flat(60), rhoM: flat(280), rhoL: 800, D: 0.254, bypass: 0, fric: 1e5, s0: 4050, slug0: 5 });
    add('Pig already in the line: remaining transit time', (10000 - 4050) / 2, pg.transit, 1e-6, '(L − x0) / v, s');
    add('Pig already in the line: liquid ahead includes the initial slug', 5 + 0.3 * ((Math.PI * 0.254 ** 2) / 4) * (10000 - 4050), pg.received + pg.inPipe, 1e-6, 'm³');
    const ex = REF.methanolExample, Kx = methanolK(ex.psia / 14.5038, (ex.degF - 32) / 1.8);
    add('Methanol vapour–liquid K-value against the published worked example', ex.K, Kx, 4e-5, '1000 psia, 10 °F, 25 wt % methanol: published K = 0.00093 (worked with rounded intermediates; the unrounded correlation gives 0.00096)');
    add('Methanol in the gas for that example', 199, (1000 * Kx * ex.x * 32.042) / VM_STD, 8, 'Published 199 kg per million Sm³ (12.4 lbm/MMscf)');
    const fmv = { at: () => ({ rhoG: 50, rhoL: 800 }) }, sv = { s: [50, 150, 250, 350], z: [0, -5, -5, 0], dz: [-5, 0, 0, 5], ds: 100, A: 0.05, P: [50, 50, 50, 50], T: [20, 20, 20, 20], holdup: [0.75, 0.75, 0.75, 0.75] }, so = settleOut(sv, fmv);
    add('Liquid seal of a U-shaped dip: head a trapped gas pocket must lift', (800 * G * 5 * 0.5) / 1e5, so.seals.up[0], 1e-9, 'Dip full, both legs half full: ρ g Δz × 0.5, bar');
    add('Valve constant from the flow coefficient Cv', 6.30902e-5 * Math.sqrt(1000 / 6894.757), slugModel({ Cv: 1 }).Kc, 5e-9, 'w = Cv × 6.309e-5 m³/s √(Δp/psi × 1000/ρ) ρ = 2.4027e-5 Cv √(ρ Δp)');
  }
  { // nonlinear MPC and gain scheduling on the riser model
    const sm = slugModel({ D: 0.2, Lp: 4000, Lr: 250, theta: 0.02, wG: 0.6, wL: 12, Ps: 20e5, Cv: 200, rhoGnom: 25 }), zc = sm.critical();
    if (zc !== null && zc < 0.7) {
      const z0 = Math.min(1, 1.3 * zc), nm = nmpc(sm, { z0, tEnd: 5400, ts: 120, zMax: 1, duMax: 0.2 }), ol = integrateStiff(sm.f(z0), nm.yEnd.map((x, i) => x * (i === 1 ? 1.002 : 1)), 0, 1800, { rtol: 1e-5, atol: 1e-4, hInit: 1 });
      add('Nonlinear MPC holds the riser at an open-loop unstable opening', 0, nm.amp, 0.2, `Pressure swing over the last third of 1.5 h at ${(100 * z0).toFixed(0)} % opening (critical ${(100 * zc).toFixed(0)} %), bar`);
      add('The same operating point diverges without control', 1, sm.growth(z0) > 0 && ol.t.length > 2 ? 1 : 0, 0, 'Positive growth rate of the open-loop linearisation');
      add('Nonlinear MPC keeps the inlet pressure at its set-point', sm.steady(z0).Pp / 1e5, mean(nm.pIn.slice(-8)), 0.2, 'bara');
    }
  }
  { // critical-flow relations against the values printed in the fluids library (API 520 gas sizing)
    const rho = (50e5 * 0.02) / (R * 300);
    add('Choked nozzle flux against the API 520 coefficient C', 0.02669419967057233 / 0.03948, orificeFlux(50e5, 300, 0, { k: 1.35, mw: 0.02 }) / Math.sqrt(50e5 * rho), 1e-9, 'C(k = 1.35) = 0.0266942 = 0.03948 √(k (2/(k+1))^((k+1)/(k−1)))');
    add('Sub-critical nozzle flux against the API 520 coefficient F2', 0.8600724121105563, orificeFlux(50e5, 300, 35e5, { k: 1.8, mw: 0.02 }) / Math.sqrt(2 * 50e5 * rho * 0.3), 1e-9, 'F2(k = 1.8, r = 0.7) = 0.86007');
  }
  return out;
}

// ---- test problems, reference data sets and provenance ---------------------------------------------------------------
const S3 = Math.sqrt(3);
/** Hock–Schittkowski test problems (objective, constraints as ≥ 0 or = 0, bounds, standard starting point) as coded in the CUTEst SIF files. */
export const HS = {
  HS12: { f: (x) => 0.5 * x[0] ** 2 + x[1] ** 2 - x[0] * x[1] - 7 * x[0] - 7 * x[1], ineq: [(x) => 25 - 4 * x[0] ** 2 - x[1] ** 2], x0: [0, 0] },
  HS21: { f: (x) => 0.01 * x[0] ** 2 + x[1] ** 2 - 100, ineq: [(x) => 10 * x[0] - x[1] - 10], lo: [2, -50], hi: [50, 50], x0: [-1, -1] },
  HS22: { f: (x) => (x[0] - 2) ** 2 + (x[1] - 1) ** 2, ineq: [(x) => 2 - x[0] - x[1], (x) => x[1] - x[0] ** 2], x0: [2, 2] },
  HS24: { f: (x) => (((x[0] - 3) ** 2 - 9) * x[1] ** 3) / (27 * S3), ineq: [(x) => x[0] / S3 - x[1], (x) => x[0] + S3 * x[1], (x) => 6 - x[0] - S3 * x[1]], lo: [0, 0], hi: [Infinity, Infinity], x0: [1, 0.5] },
  HS29: { f: (x) => -x[0] * x[1] * x[2], ineq: [(x) => 48 - x[0] ** 2 - 2 * x[1] ** 2 - 4 * x[2] ** 2], x0: [1, 1, 1] },
  HS30: { f: (x) => x[0] ** 2 + x[1] ** 2 + x[2] ** 2, ineq: [(x) => x[0] ** 2 + x[1] ** 2 - 1], lo: [1, -10, -10], hi: [10, 10, 10], x0: [1, 1, 1] },
  HS34: { f: (x) => -x[0], ineq: [(x) => x[1] - Math.exp(x[0]), (x) => x[2] - Math.exp(x[1])], lo: [0, 0, 0], hi: [100, 100, 10], x0: [0, 1.05, 2.9] },
  HS35: { f: (x) => 9 - 8 * x[0] - 6 * x[1] - 4 * x[2] + 2 * x[0] ** 2 + 2 * x[1] ** 2 + x[2] ** 2 + 2 * x[0] * x[1] + 2 * x[0] * x[2], ineq: [(x) => 3 - x[0] - x[1] - 2 * x[2]], lo: [0, 0, 0], hi: [Infinity, Infinity, Infinity], x0: [0.5, 0.5, 0.5] },
  HS39: { f: (x) => -x[0], eq: [(x) => x[1] - x[0] ** 3 - x[2] ** 2, (x) => x[0] ** 2 - x[1] - x[3] ** 2], x0: [2, 2, 2, 2] },
  HS43: { f: (x) => x[0] ** 2 + x[1] ** 2 + 2 * x[2] ** 2 + x[3] ** 2 - 5 * x[0] - 5 * x[1] - 21 * x[2] + 7 * x[3], ineq: [(x) => 8 - x[0] ** 2 - x[1] ** 2 - x[2] ** 2 - x[3] ** 2 - x[0] + x[1] - x[2] + x[3], (x) => 10 - x[0] ** 2 - 2 * x[1] ** 2 - x[2] ** 2 - 2 * x[3] ** 2 + x[0] + x[3], (x) => 5 - 2 * x[0] ** 2 - x[1] ** 2 - x[2] ** 2 - 2 * x[0] + x[1] + x[3]], x0: [0, 0, 0, 0] },
  HS65: { f: (x) => (x[0] - x[1]) ** 2 + (x[0] + x[1] - 10) ** 2 / 9 + (x[2] - 5) ** 2, ineq: [(x) => 48 - x[0] ** 2 - x[1] ** 2 - x[2] ** 2], lo: [-4.5, -4.5, -5], hi: [4.5, 4.5, 5], x0: [-5, 5, 0] },
  HS71: { f: (x) => x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2], ineq: [(x) => x[0] * x[1] * x[2] * x[3] - 25], eq: [(x) => x[0] ** 2 + x[1] ** 2 + x[2] ** 2 + x[3] ** 2 - 40], lo: [1, 1, 1, 1], hi: [5, 5, 5, 5], x0: [1, 5, 5, 1] },
  HS100: { f: (x) => (x[0] - 10) ** 2 + 5 * (x[1] - 12) ** 2 + x[2] ** 4 + 3 * (x[3] - 11) ** 2 + 10 * x[4] ** 6 + 7 * x[5] ** 2 + x[6] ** 4 - 4 * x[5] * x[6] - 10 * x[5] - 8 * x[6], ineq: [(x) => 127 - 2 * x[0] ** 2 - 3 * x[1] ** 4 - x[2] - 4 * x[3] ** 2 - 5 * x[4], (x) => 282 - 7 * x[0] - 3 * x[1] - 10 * x[2] ** 2 - x[3] + x[4], (x) => 196 - 23 * x[0] - x[1] ** 2 - 6 * x[5] ** 2 + 8 * x[6], (x) => -4 * x[0] ** 2 - x[1] ** 2 + 3 * x[0] * x[1] - 2 * x[2] ** 2 - 5 * x[5] + 11 * x[6]], x0: [1, 2, 0, 4, 0, 1, 1] },
};
const once = (fn) => { let r, done = false; return () => { if (!done) { r = fn(); done = true; } return r; }; };
/** Blind prediction of the nitrogen vessel blowdown with the suite's blowdown model and its default gas-to-wall film coefficient (25 W/m²K). */
const n2Run = once(() => {
  const g = REF.n2.vessel, fl = makeFluid({ comp: { N2: 100 } }), Zc = new Map(), Zf = (P, T) => { const k = Math.round(P / 2e4) + ':' + Math.round(T); if (!Zc.has(k)) Zc.set(k, phaseProps(fl, [1], Math.max(P, 1e4) / 1e5, T - KEL, 'vapour').Z); return Zc.get(k); };
  const V = (Math.PI / 4) * g.diameter ** 2 * g.length, Ain = Math.PI * g.diameter * g.length + (Math.PI / 2) * g.diameter ** 2, mSteel = g.rho * ((Math.PI / 4) * ((g.diameter + 2 * g.thickness) ** 2 - g.diameter ** 2) * g.length + (Math.PI / 2) * g.diameter ** 2 * g.thickness);
  return blowdown({ V, P0: g.P0, T0: g.T0, pBack: g.pBack, area: (Math.PI / 4) * g.orifice ** 2, cd: g.cd, k: 1.4, mw: 0.028013, Z: Zf, mode: 'wall', wallC: mSteel * g.cp, wallUA: DEFAULTS.hGasWall * Ain, extUA: 5 * Ain, tAmb: g.tAmb, dt: 0.1, n: 1100, maxFactor: 1, pEnd: g.pBack * 1.0001 });
});
const atTime = (b, arr, t) => interp1(b.t, arr, clamp(t, 0, b.t[b.t.length - 1]));
/** Riser model with the published parameter set of the pipeline–riser test case: critical opening, period, steady state and slug-cycle extremes. */
const riserRun = once(() => {
  const o = REF.riser.olga, pNom = REF.riser.rows.find((r) => r.q === 'ssPin').value * 1e5, rg = (pNom * o.mwG) / (R * o.Tp), aLp = (rg * o.wL) / (rg * o.wL + o.rhoL * o.wG), sm = slugModel({ ...o, aLp }), zc = sm.critical(0.02, 1, 13);
  const pl = sm.poles(0.05).filter((e) => Math.abs(e[1]) > 1e-9).sort((a, b) => b[0] - a[0])[0], e1 = sm.steady(1), f = sm.f(1), dt = 0.25, n = Math.round((1.5 * HOUR) / dt), mn = { Pp: 1e99, Prb: 1e99, Prt: 1e99, w: 1e99 }, mx = { Pp: 0, Prb: 0, Prt: 0, w: 0 };
  let y = e1.y.map((x, i) => x * (i === 1 ? 1.05 : 1));
  for (let k = 0; k < n; k++) { y = ros2Step(f, 0, y, dt).y.map((x) => Math.max(x, 1e-9)); if (k * dt > 0.75 * HOUR) { const a = sm.alg(y, 1); for (const q of ['Pp', 'Prb', 'Prt', 'w']) { mn[q] = Math.min(mn[q], a[q]); mx[q] = Math.max(mx[q], a[q]); } } }
  // laboratory rig: flows from the stated volumetric rates (air taken at atmospheric conditions); the feed liquid fraction is made consistent with the model's own inlet pressure
  const g = REF.riser.rig, wL = (g.qWaterLmin / 60000) * g.rhoL, wG = (g.qAirLmin / 60000) * ((g.PsAtm * 0.029) / (R * g.Tp));
  let zr = null, pIn = 1.3 * g.PsAtm;
  for (let it = 0; it < 4; it++) { const rgr = (pIn * g.mwG) / (R * g.Tp), sr = slugModel({ ...g, wL, wG, Ps: g.PsAtm, aLp: (rgr * wL) / (rgr * wL + g.rhoL * wG) }); zr = sr.critical(0.03, 1, 13); try { pIn = sr.steady(zr ?? 0.2).Pp; } catch { break; } }
  return { zCrit: zc === null ? 100 : 100 * zc, period: pl ? (2 * Math.PI) / Math.abs(pl[1]) / 60 : 0, ssPin: e1.Pp / 1e5, ssPrb: e1.Prb / 1e5, ssPrt: e1.Prt / 1e5, minPin: mn.Pp / 1e5, minPrb: mn.Prb / 1e5, minPrt: mn.Prt / 1e5, minW: mn.w, maxPin: mx.Pp / 1e5, maxPrb: mx.Prb / 1e5, maxPrt: mx.Prt / 1e5, maxW: mx.w, zCritRig: zr === null ? 100 : 100 * zr };
});
/** Travelling front of the advection–dispersion solver far from the inlet: C/C0 against ξ = (x − U t)/(2 √(D t)). */
const frontRun = once(() => { const U = 1, D = 5, f = inhibitorFront({ U, L: 12000, D, n: 600, tEnd: 6000, rows: 2 }), i = f.t.length - 1; return { U, D, t: f.t[i], xc: f.xc, c: f.c[i] }; });
const hsRun = new Map();
const VALIDATION = [
  { id: REF.simc.id, title: REF.simc.title, quantity: 'Ms', unit: '–', kind: 'benchmark', source: REF.simc.source,
    columns: [{ key: 'id', label: 'Case' }, { key: 'process', label: 'Process' }, { key: 'kc', label: 'Kc' }, { key: 'ti', label: 'τI' }, { key: 'td', label: 'τD (series form)' }, { key: 'Ms', label: 'Published Ms' }],
    rows: REF.simc.rows, target: 'Ms',
    model: (r) => loopAnalysis({ k: r.k, lags: r.lags, zeros: r.zeros, ints: r.ints, delay: r.delay }, r.ki ? { ki: r.ki } : { kc: r.kc, ti: r.ti, td: r.td, form: 'series' }, { n: 2400 }).ms,
    tolerance: { mape: 1.5, maxAbs: 0.06 },
    note: 'Sensitivity peaks recomputed by the loop-analysis routine (frequency scan of the published process and controller) against the published values, which are printed to two decimals. Process E5 is entered as printed (smallest time constant 0.0008); the half-rule approximation printed in the same row corresponds to 0.008, with which the PID case gives Ms 1.83. Achieved: mean error 0.5 %, largest 0.04.' },
  { id: REF.n2.idP, title: 'Nitrogen vessel blowdown from 150 bara: pressure against time (experiment I1)', quantity: 'Pressure', unit: 'bara', kind: 'experiment', source: REF.n2.source,
    columns: [{ key: 't', label: 'Time', unit: 's' }, { key: 'P', label: 'Measured pressure', unit: 'bara' }], rows: REF.n2.pressure, target: 'P',
    model: (r) => atTime(n2Run(), n2Run().P, r.t) / 1e5, tolerance: { mape: 25, bias: 4 },
    note: 'Blind prediction with the lumped blowdown model (real-gas Z from the kernel equation of state, k = 1.4, discharge coefficient 0.8 from the cited input file, wall heat sink with the suite default film coefficient of 25 W/m²K). The model empties the vessel too fast: within 3 % for the first 10 s, then 10–28 % low from 15 s onwards (mean error 18 %), probably because the fixed film coefficient under-states the natural convection in dense cold nitrogen and the measured gas stays warmer. The tolerance states that miss; it is not a design-accuracy claim.' },
  { id: REF.n2.idT, title: 'Nitrogen vessel blowdown from 150 bara: gas temperature (upper and lower thermocouples)', quantity: 'Gas temperature', unit: 'K', kind: 'experiment', source: REF.n2.source,
    columns: [{ key: 't', label: 'Time', unit: 's' }, { key: 'sensor', label: 'Thermocouple' }, { key: 'T', label: 'Measured temperature', unit: 'K' }], rows: REF.n2.temperature, target: 'T',
    model: (r) => atTime(n2Run(), n2Run().T, r.t), tolerance: { mape: 12, bias: 30 },
    note: 'The model has one bulk gas temperature; the experiment shows 20–30 K of stratification between the two thermocouples. The predicted minimum is 15–35 K colder than the measured range (conservative for minimum-metal-temperature screening), for the same reason as the pressure miss.' },
  { id: REF.riser.id, title: 'Pipeline–riser slugging test case: stability limit, period and slug-cycle extremes (OLGA reference values; rig experiment)', quantity: 'Mixed quantities (see rows)', unit: 'as listed', kind: 'benchmark', source: REF.riser.source,
    columns: [{ key: 'quantity', label: 'Quantity' }, { key: 'unit', label: 'Unit' }, { key: 'value', label: 'Reference value' }], rows: REF.riser.rows, target: 'value',
    model: (r) => riserRun()[r.q], tolerance: { mape: 40 },
    note: 'The four-state riser model run with the published parameter set (orifice coefficients, level factor and valve constant of the cited authors; friction and feed liquid fraction by their relations). The first thirteen rows are OLGA simulation values, the last is the measured critical opening of the laboratory rig. For the OLGA case the critical opening (4.6 against 5 %), the period, the steady state and the pressure extremes agree within 12 %. Three rows are missed and dominate the mean error: the minimum outlet rate (0 against 0.79 kg/s), the maximum outlet rate (about 70 % too high; the published model itself reports 68 % and 32 % on these two) and the rig, where the model with the published rig parameters gives 35 % against the measured 15 % (the basis of the stated air rate and the rig valve characteristic are not given in the source). The model is therefore suited to locating the stability limit of a tuned case, not to predicting slug flow peaks.' },
  { id: REF.hs.id, title: 'Hock–Schittkowski constrained test problems: optimal objective values', quantity: 'f*', unit: '–', kind: 'benchmark', source: REF.hs.source,
    columns: [{ key: 'problem', label: 'Problem' }, { key: 'n', label: 'Variables' }, { key: 'm', label: 'Constraints' }, { key: 'fStar', label: 'Best known objective' }], rows: REF.hs.rows, target: 'fStar',
    model: (r) => { if (!hsRun.has(r.problem)) { const p = HS[r.problem]; hsRun.set(r.problem, sqp(p.f, p.x0, { eq: p.eq || [], ineq: p.ineq || [], lo: p.lo, hi: p.hi }).f); } return hsRun.get(r.problem); },
    tolerance: { mape: 0.01, maxAbs: 1e-3 },
    note: 'The SQP solver started from the standard starting points with finite-difference derivatives. All thirteen optima are reproduced to five significant figures or better.' },
  { id: REF.erfc.id, title: 'Advection–dispersion front against the tabulated complementary error function', quantity: 'erfc ξ = 2 C/C0', unit: '–', kind: 'benchmark', source: REF.erfc.source,
    columns: [{ key: 'xi', label: 'ξ = (x − U t) / (2 √(D t))' }, { key: 'erfc', label: 'Tabulated erfc ξ' }], rows: REF.erfc.rows, target: 'erfc',
    model: (r) => { const f = frontRun(); return 2 * interp1(f.xc, f.c, f.U * f.t + 2 * r.xi * Math.sqrt(f.D * f.t)); }, tolerance: { mape: 3, maxAbs: 0.01 },
    note: 'Concentration profile of the finite-volume inhibitor-front solver (600 cells, Péclet number 2400 at the outlet, so the second Ogata–Banks term is negligible) sampled at the tabulated ξ ahead of the front centre.' },
];
/** Where the literature constants and rules of this suite were checked. */
export const PROVENANCE = [
  { item: 'SIMC PI rule: Kc = τ/(K (τc + θ)), τI = min(τ, 4 (τc + θ)), τc = θ', used: 'tuningRules()', source: 'S. Skogestad (2003), Simple analytic rules for model reduction and PID controller tuning, J. Process Control 13, 291–309, eqs 23–25 and 28', url: 'https://folk.ntnu.no/skoge/publications/2003/tuningPID/finalpaper.pdf', retrieved: '2026-10-08', status: 'verified', note: 'Formulas identical. The published margins of the rule (GM 3.14, PM 61.4°, Ms 1.59, Mt 1.00; integrating process 2.96, 46.9°, 1.70, 1.30) are reproduced by loopAnalysis() within 0.01.' },
  { item: 'SIMC PID settings for a first-order process with delay', used: 'tuningRules()', source: 'C. Grimholt, S. Skogestad (2013), Optimal PID-control on first order plus time delay systems & verification of the SIMC rules, DYCOPS 2013; S. Skogestad, C. Grimholt (2012), The SIMC method for smooth PID controller tuning (improved PI rule)', url: 'https://folk.ntnu.no/skoge/publications/2013/grimholt-dycops/0122.pdf', retrieved: '2026-10-08', status: 'corrected', note: 'Old: Kc = (τ + θ/3)/(K (τc + θ)), τI = min(τ + θ/3, 4 (τc + θ)), τD = θ/3 applied directly in the ideal PID form with τc = θ — a mixture of the improved PI rule and the PID extension that appears in neither source. New: series-form Kc = τ/(K (τc + θ)), τI = min(τ, 4 (τc + θ)), τD = θ/3 with τc = θ/2, converted to the ideal form (Kc f, τI f, τD/f, f = 1 + τD/τI). For τ = 10, θ = 1 the ideal gain changes from 5.17/K to 7.04/K.' },
  { item: 'Ziegler–Nichols closed-loop settings (P 0.5 Ku; PI 0.45 Ku, Pu/1.2; PID 0.6 Ku, Pu/2, Pu/8, ideal form)', used: 'tuningRules()', source: 'Skogestad (2003) section 5.3, and Wikipedia, Ziegler–Nichols method', url: 'https://en.wikipedia.org/wiki/Ziegler%E2%80%93Nichols_method', retrieved: '2026-10-08', status: 'verified', note: 'All coefficients identical (0.8333 Tu = Pu/1.2). No difference.' },
  { item: 'Low-order riser slugging model: structure, level relation, orifice equations, top liquid fraction, valve equation; Kh = 0.7', used: 'slugModel()', source: 'E. Jahanshahi, S. Skogestad (2011), Simplified dynamical models for control of severe slugging in multiphase risers, 18th IFAC World Congress, eqs 1–43; E. Jahanshahi (2013), PhD thesis, NTNU, Table 2.1', url: 'https://folk.ntnu.no/skoge/publications/2011/jahanshahi_ifac-2011/0981.pdf', retrieved: '2026-10-08', status: 'verified', note: 'Equations 10–11, 31–37 and 41–43 match the code. Deliberate differences, now options: wall friction by Haaland with roughness instead of 0.0056 + 0.5 Re^-0.32 (fric: "dkm" reproduces the paper), friction factors frozen at the nominal velocities, riser friction length without the topside length unless Lh is given. Default low-point liquid coefficient 0.3 against the published fitted 0.281 for the OLGA case and 0.157 for the rig: a tuning parameter, kept as an input. With the published parameter set the model gives a critical opening of 4.6 % (published 5 %) and a period of 16.3 min (15.6).' },
  { item: 'Valve constant from Cv: w = 2.403e-5 Cv f(z) √(ρ Δp)', used: 'slugModel(), calModel()', source: 'Unit conversion of the Cv definition (US gal/min of water at 1 psi): 6.30902e-5 m³/s × √(1000 kg/m³ / 6894.757 Pa)', url: 'https://en.wikipedia.org/wiki/Flow_coefficient', retrieved: '2026-10-08', status: 'verified', note: 'Derived value 2.4027e-5; the code uses 2.403e-5 (0.01 % difference).' },
  { item: 'ω-method: mass flux of a non-flashing two-phase mixture and critical pressure ratio η² + (ω² − 2ω)(1 − η)² + 2ω² ln η + 2ω² (1 − η) = 0, ω = α/k', used: 'omegaFlux()', source: 'C. Hős, G. Burhani (2022), On the effect of mass fraction of frozen mixture flow on the dynamic performance of a direct spring operated safety valve, CMFF\'22 paper 114, eqs 13–15 (after J.C. Leung)', url: 'https://www.cmff.hu/papers/CMFF22_Final_Paper_PDF_114.pdf', retrieved: '2026-10-08', status: 'corrected', note: 'Flux and critical-ratio equations identical. Removed: an explicit fit η_c = 0.55 + 0.217 ln ω − 0.046 (ln ω)² used as a fallback, which could not be traced to an open source; the implicit equation always has a root in (0, 1) and is now solved in every case. The API 520 functions of the fluids library (safety_valve.py) were read: that library contains gas, steam and liquid sizing only and no ω-method, so it could not serve as the source for this item.' },
  { item: 'Critical pressure ratio and choked mass flux of an ideal gas', used: 'criticalRatio(), orificeFlux()', source: 'C. Bell, fluids (MIT licence), fluids/safety_valve.py, API520_C and API520_F2 docstrings (API 520 Part I gas sizing relations)', url: 'https://raw.githubusercontent.com/CalebBell/fluids/master/fluids/safety_valve.py', retrieved: '2026-10-08', status: 'verified', note: 'C ∝ √(k (2/(k + 1))^((k + 1)/(k − 1))) and the sub-critical factor F2 are the same isentropic nozzle relations as in the code; verify() checks the choked limit against the closed-form blowdown solutions.' },
  { item: 'Relaxation time of the homogeneous relaxation model', used: 'hrmRelaxationTime(), blowdown({ relax })', source: 'D. Schmidt, R. Maulik, K. Lyras (2021), Machine-learning accelerated turbulence modelling of transient flashing jets, arXiv:2109.15203, eqs 4–7 (correlation of Downar-Zapolski et al. 1996)', url: 'https://arxiv.org/pdf/2109.15203', retrieved: '2026-10-08', status: 'verified', note: 'Low-pressure branch (≤ 10 bar): Θ0 = 6.51e-4 s, exponents −0.257 and −2.24, identical. High-pressure branch: exponents −0.54 and −1.76 identical; the source prints Θ0 = 3.84·10^7 s, which is not physical, and 3.84e-7 s is used: that constant is NOT confirmed. The correlation was fitted to flashing water; its use for gas coming out of solution from oil is an assumption, so a user-set relaxation time is offered.' },
  { item: 'Methanol loss to the gas phase (vapour–aqueous K-value)', used: 'methanolK(), inhibitorDose()', source: 'M. Moshfeghian (2011), A simple model for estimation of methanol loss to vapor phase, Campbell Tip of the Month, August 2011, eqs 4–6 and Appendix A', url: 'https://jmcampbell.com/tip-of-the-month/2011/08/a-simple-model-for-estimation-of-methanol-loss-to-vapor-phase', retrieved: '2026-10-08', status: 'corrected', note: 'Old: modified Raoult\'s law with an activity coefficient of 1.6 recalled from memory; it gave y = 6.6e-5 for the published example (1000 psia, 10 °F, 25 wt %) against 1.47e-4, a factor 2.2 low. New: the published correlation; equation 3 is an image on the page and was reconstructed as the Wilson form K = exp[5.37 (1 + ω*)(1 − 1/T*)]/P*. With the rounded intermediates printed in the worked example it returns the published K = 0.00093; unrounded it gives 0.00096 (3 % higher, 206 against 199 kg per million Sm³). Applied inside −23…38 °C and 7…345 bar.' },
  { item: 'Methanol and glycol losses to the hydrocarbon liquid; glycol loss to the gas', used: 'inhibitorDose()', source: 'M. Moshfeghian (2010), Determination of traces of methanol in the TEG dehydrated gas, Campbell Tip of the Month, October 2010', url: 'https://www.jmcampbell.com/tip-of-the-month/2010/10/determination-of-traces-of-methanol-in-the-teg-dehydrated-gas/', retrieved: '2026-10-08', status: 'corrected', note: 'Old: distribution coefficients 6e-3 (alcohols) and 2e-4 (glycols) with an exponential temperature factor, recalled. New: methanol 0.4 kg per m³ of hydrocarbon liquid (planning value of the source); MEG 3.5 L per million Sm³ of gas to the hydrocarbon liquid and no loss to the gas ("negligible"). The same glycol values are applied to DEG and TEG and the methanol value to ethanol without a source of their own.' },
  { item: 'Antoine vapour pressure of ethanol (A 8.20417, B 1642.89, C 230.3; mmHg, °C)', used: 'inhibitorDose() for ethanol only', source: 'Dortmund Data Bank online Antoine calculation (methanol page read: 8.08097, 1582.27, 239.7 for 15–100 °C)', url: 'http://ddbonline.ddbst.com/AntoineCalculation/AntoineCalculationCGI.exe?component=Methanol', retrieved: '2026-10-08', status: 'unverified', note: 'The methanol constants formerly in the code agreed with the page read, but methanol now uses the K-value correlation. The ethanol page could not be retrieved, so the ethanol constants and its activity coefficient of 1.6 remain unchecked. The glycol Antoine sets were removed.' },
  { item: 'Bypass-pig discharge coefficient (pressure-loss coefficient of the bypass)', used: 'pigRun(), input pigCd', source: 'J.E. Azpiroz, M.H.W. Hendrix, W.P. Breugem, R.A.W.M. Henkes (2015), CFD modelling of bypass pigs with a deflector disk, 17th Int. Conf. Multiphase Technology, abstract', url: 'https://research.tudelft.nl/en/publications/cfd-modelling-of-bypass-pigs-with-a-deflector-disk/', retrieved: '2026-10-08', status: 'corrected', note: 'Old default Cd = 0.7 (K = 1/Cd² = 2.04). The source gives K = 1–1.5 as the range commonly used for plain bypass pigs and up to 4 with a deflector disk; new default Cd = 0.9 (K = 1.23). K is taken on the bypass velocity; the abstract does not state its reference velocity. The force balance (seal friction + weight component) and the orifice form are standard and were not compared with a specific source.' },
  { item: 'Laminar Taylor–Aris dispersion coefficient D (1 + Pe²/48)', used: 'theory reference for the dispersion input', source: 'Wikipedia, Taylor dispersion', url: 'https://en.wikipedia.org/wiki/Taylor_dispersion', retrieved: '2026-10-08', status: 'verified', note: 'Factor 48 confirmed; the laminar form is not used in the engine (pipeline flow is turbulent).' },
  { item: 'Turbulent axial dispersion coefficient 10.1 a u* (Taylor 1954)', used: 'run(): dispersion of the inhibitor front, before the user multiplier', source: 'G.I. Taylor (1954), Proc. R. Soc. A 223, 446–468 (not openly readable); course outline of MIT 2.27 lists the paper', url: 'https://ocw.mit.edu/courses/2-27-turbulent-flow-and-transport-spring-2002/5091c145ca7a3fe1d6d759cae9b803e2_9_Taylor_dispersion.pdf', retrieved: '2026-10-08', status: 'unverified', note: 'No openly readable page stating the coefficient 10.1 was found; the pages opened only cite the paper. The value is multiplied by a user factor (default 150) for slug mixing, which dominates the result and is itself an engineering assumption to be calibrated with tracer data.' },
  { item: 'Unscented transform: sigma points and weights', used: 'ukf()', source: 'Wikipedia, Kalman filter, section "Unscented Kalman filter" (after Wan & van der Merwe 2000)', url: 'https://en.wikipedia.org/wiki/Kalman_filter', retrieved: '2026-10-08', status: 'verified', note: 'W0ᵃ = (α²κ − L)/(α²κ), W0ᶜ = W0ᵃ + 1 − α² + β, Wj = 1/(2α²κ), points x̂ ± α√κ A_j: identical. Defaults α = 1, β = 2 (Gaussian), κ = 1.5 L (the page marks its κ recommendation as needing a citation).' },
  { item: 'Affinity (fan) laws Q ∝ N, H ∝ N², P ∝ N³ and the surge-margin definition 100 (ṁ − ṁ_surge)/ṁ', used: 'pumpSim(), compressorCurve(), surgeMargin()', source: 'Wikipedia, Affinity laws; Wikipedia, Compressor map', url: 'https://en.wikipedia.org/wiki/Compressor_map', retrieved: '2026-10-08', status: 'verified', note: 'Identical. The default surge control line of 10 % and the surge point at 62 % of the design flow on the generated map are design assumptions (inputs), not literature constants.' },
  { item: 'Wet-insulation conductivity (glass syntactic polyurethane, 0.17 W/m/K in the reference case)', used: 'default coating layer (from the base case)', source: '3M Glass Bubbles for insulation and buoyancy, sales card 2020, calculated conductivity of GSPU at 40 vol % loading', url: 'https://multimedia.3m.com/mws/media/1965520O/14473-glass-bubbles-insulation-and-buoyancy-sales-card-2020.pdf', retrieved: '2026-10-08', status: 'verified', note: 'The card lists 0.130–0.173 W/m/K depending on the bubble grade (0.19 for unfilled polyurethane); 0.17 lies inside that range, at its upper end.' },
  { item: 'Soil thermal conductivity default 1.5 W/m/K', used: 'input kSoil (buried lines)', source: 'World Oil, August 2000, Deepwater soil thermally insulates buried flowlines (article behind a subscription wall)', url: 'https://www.worldoil.com/magazine/2000/august-2000/special-report/deepwater-soil-thermally-insulates-buried-flowlines', retrieved: '2026-10-08', status: 'unverified', note: 'The page opened shows only the title; the numbers are not readable without a subscription. Search summaries quote about 1.0 W/m/K for deep-water Gulf of Mexico clay and 0.7–2.5 for other soils, which would bracket the default, but they were not read at the source. Not used in the reference case (line not buried).' },
  { item: 'Minimum design metal temperature default −29 °C (carbon steel)', used: 'input tMinDesign', source: 'ASME B31.3 Table A-1 / Fig. 323.2.2A (standard not openly readable)', url: 'https://amarineblog.com/2021/03/10/api-570-quiz-asme-31-3-mdmt-and-impact-test/', retrieved: '2026-10-08', status: 'unverified', note: 'The open page read discusses the impact-test curves but does not state the −29 °C (−20 °F) limit. The value is a user input and must be taken from the project material specification.' },
  { item: 'Methanol price default', used: 'input inhPrice', source: 'Methanex, posted regional contract prices: North America non-discounted reference price valid 1–31 October 2026', url: 'https://www.methanex.com/about-methanol/pricing/', retrieved: '2026-10-08', status: 'corrected', note: 'Old default 450 $/m³ (recalled). Posted price USD 1,450 per tonne = 4.36 $/gal → 1,148 $/m³ at 792 kg/m³; new default 1,150 $/m³. The MEG price used in one preset (1,100 $/m³) was not checked.' },
  { item: 'Electricity price default', used: 'input elecPrice', source: 'U.S. Energy Information Administration, Monthly Energy Review September 2026, Figure 9.2 (average price of electricity to ultimate customers, industrial sector, June 2026)', url: 'https://www.eia.gov/totalenergy/data/monthly/pdf/sec9_10.pdf', retrieved: '2026-10-08', status: 'corrected', note: 'Old default 0.12 $/kWh (recalled). The figure labels 9.17 ¢/kWh for the industrial sector in June 2026 (read from the chart labels, the lowest of the four sector values); new default 0.092 $/kWh. Offshore self-generated power usually costs more: set the site value.' },
  { item: 'Steel and coating thermal properties (steel 45 W/m/K, 7850 kg/m³, 470 J/kg/K; coating 900 kg/m³, 1700 J/kg/K; soil 1900 kg/m³, 1300 J/kg/K)', used: 'wallLayers(), cooldown thermal mass', source: 'none opened: generic handbook values', url: 'https://multimedia.3m.com/mws/media/1965520O/14473-glass-bubbles-insulation-and-buoyancy-sales-card-2020.pdf', retrieved: '2026-10-08', status: 'unverified', note: 'Only the coating conductivity was checked (see the wet-insulation entry; the card read gives no heat capacity). The thermal mass can be overridden by the network suite value or calibrated with the thermal-mass multiplier against cooldown records.' },
  { item: 'Erosional velocity 122/√ρ m/s (API RP 14E with C = 100 in US units) and the 0.5–5 m/s pig velocity guide', used: 'operating envelope; pigging warnings', source: 'API RP 14E (not openly readable); the factor 122 is the unit conversion of C = 100 (ft/s)(lb/ft³)^0.5', url: 'https://en.wikipedia.org/wiki/Flow_coefficient', retrieved: '2026-10-08', status: 'unverified', note: 'Neither the recommended practice nor a pigging guideline was opened; the address given is only the unit-definition page used for the conversions in this suite. Both limits are screening conventions and should be replaced by project criteria.' },
  { item: 'Selection hints for low-dosage hydrate inhibitors (kinetic inhibitors up to about 10 °C sub-cooling, anti-agglomerants up to about 50 % water cut)', used: 'note under the chemical-injection table', source: 'none opened', url: 'https://www.jmcampbell.com/tip-of-the-month/2010/10/determination-of-traces-of-methanol-in-the-teg-dehydrated-gas/', retrieved: '2026-10-08', status: 'unverified', note: 'Rules of thumb recalled, not found on the pages opened for the inhibitor-loss entries (the address given is that page); they affect only an advisory sentence, no computed number.' },
  { item: 'Hock–Schittkowski problem data and optimal values', used: 'HS, verify(), reference data set', source: 'CUTEst SIF files HS12 … HS100', url: 'https://bitbucket.org/optrove/sif/raw/HEAD/HS71.SIF', retrieved: '2026-10-08', status: 'verified', note: 'Objective, constraints, bounds and starting points of thirteen problems transcribed from the SIF files; optimal values from their LO SOLTN lines.' },
];
const CAL_SAMPLE = [{ tShutH: 2, t0: 52, tCold: 42.85, tBlowMin: 10, p0: 82, pBlow: 71.93, vGas: 1.6, pigArrival: 5.622, pigDp: 1.29, tValveS: 2, valvePos: -0.003671, tStepS: 30, stepPv: 0, tClS: 60, clPv: 0.2991, pTrue: 60, pvRead: 60.47, zChokePct: 15, wChoke: 30, dpChoke: 85.29, qPumpM3h: 40, nPumpPct: 100, pumpHeadM: 732, qCompM3s: 0.17, nCompPct: 100, compHeadKj: 123.9, tTracerH: 3, uLiq: 1.2, tracerC: 0.003308, qInjM3d: 10, mWaterKgs: 1.5, doseMeas: 4.987, gelLenM: 0, restartP: 33.55 }, { tShutH: 4, t0: 52, tCold: 37.49, tBlowMin: 20, p0: 82, pBlow: 65.5, vGas: 1.9, pigArrival: 4.328, pigDp: 1.332, tValveS: 4, valvePos: 0.07938, tStepS: 60, stepPv: -0.008529, tClS: 120, clPv: 1.127, pTrue: 70, pvRead: 71.22, zChokePct: 20, wChoke: 34, dpChoke: 49.63, qPumpM3h: 70, nPumpPct: 100, pumpHeadM: 732.4, qCompM3s: 0.19, nCompPct: 100, compHeadKj: 122.8, tTracerH: 3.4, uLiq: 1.2, tracerC: 0.006228, qInjM3d: 15, mWaterKgs: 1.5, doseMeas: 7.175, gelLenM: 500, restartP: 36.46 }, { tShutH: 6, t0: 52, tCold: 33.62, tBlowMin: 30, p0: 82, pBlow: 57.72, vGas: 2.2, pigArrival: 3.545, pigDp: 1.322, tValveS: 6, valvePos: 0.2179, tStepS: 120, stepPv: -0.09771, tClS: 200, clPv: 0.9797, pTrue: 80, pvRead: 79.88, zChokePct: 25, wChoke: 38, dpChoke: 33.15, qPumpM3h: 100, nPumpPct: 100, pumpHeadM: 692.6, qCompM3s: 0.21, nCompPct: 100, compHeadKj: 121.5, tTracerH: 3.8, uLiq: 1.2, tracerC: 0.004451, qInjM3d: 20, mWaterKgs: 2, doseMeas: 7.078, gelLenM: 1000, restartP: 37.42 }, { tShutH: 8, t0: 52, tCold: 28.72, tBlowMin: 45, p0: 82, pBlow: 49.31, vGas: 2.5, pigArrival: 3.021, pigDp: 1.336, tValveS: 8, valvePos: 0.332, tStepS: 200, stepPv: -0.2132, tClS: 300, clPv: 0.9069, pTrue: 90, pvRead: 91.43, zChokePct: 30, wChoke: 42, dpChoke: 24.56, qPumpM3h: 130, nPumpPct: 100, pumpHeadM: 669.5, qCompM3s: 0.23, nCompPct: 100, compHeadKj: 116, tTracerH: 4.1, uLiq: 1.2, tracerC: -0.005381, qInjM3d: 25, mWaterKgs: 2, doseMeas: 8.923, gelLenM: 2000, restartP: 41.73 }, { tShutH: 10, t0: 48, tCold: 23.21, tBlowMin: 60, p0: 78, pBlow: 38.71, vGas: 2.8, pigArrival: 2.623, pigDp: 1.341, tValveS: 10, valvePos: 0.4497, tStepS: 300, stepPv: -0.3361, tClS: 400, clPv: 0.9492, pTrue: 100, pvRead: 100.5, zChokePct: 35, wChoke: 46, dpChoke: 19.34, qPumpM3h: 160, nPumpPct: 100, pumpHeadM: 634.2, qCompM3s: 0.25, nCompPct: 100, compHeadKj: 111.7, tTracerH: 4.4, uLiq: 1.2, tracerC: 0.0315, qInjM3d: 30, mWaterKgs: 2.5, doseMeas: 8.624, gelLenM: 3000, restartP: 46.02 }, { tShutH: 12, t0: 48, tCold: 20.51, tBlowMin: 90, p0: 78, pBlow: 27.74, vGas: 3.1, pigArrival: 2.275, pigDp: 1.328, tValveS: 14, valvePos: 0.6048, tStepS: 450, stepPv: -0.4632, tClS: 500, clPv: 0.9471, pTrue: 110, pvRead: 110.1, zChokePct: 40, wChoke: 50, dpChoke: 15.28, qPumpM3h: 190, nPumpPct: 100, pumpHeadM: 587.4, qCompM3s: 0.27, nCompPct: 100, compHeadKj: 107.8, tTracerH: 4.7, uLiq: 1.2, tracerC: 0.4833, qInjM3d: 35, mWaterKgs: 2.5, doseMeas: 9.751, gelLenM: 4000, restartP: 49.34 }, { tShutH: 16, t0: 48, tCold: 16.32, tBlowMin: 120, p0: 78, pBlow: 19.56, vGas: 3.4, pigArrival: 2.01, pigDp: 1.296, tValveS: 18, valvePos: 0.6943, tStepS: 600, stepPv: -0.5611, tClS: 650, clPv: 0.9456, pTrue: 120, pvRead: 120.1, zChokePct: 50, wChoke: 54, dpChoke: 9.732, qPumpM3h: 210, nPumpPct: 100, pumpHeadM: 556.5, qCompM3s: 0.29, nCompPct: 100, compHeadKj: 102, tTracerH: 5, uLiq: 1.2, tracerC: 0.9263, qInjM3d: 40, mWaterKgs: 3, doseMeas: 9.282, gelLenM: 6000, restartP: 58.65 }, { tShutH: 20, t0: 45, tCold: 12.59, tBlowMin: 150, p0: 75, pBlow: 13.31, vGas: 2, pigArrival: 4.092, pigDp: 1.318, tValveS: 24, valvePos: 0.8302, tStepS: 800, stepPv: -0.6515, tClS: 800, clPv: 0.9369, pTrue: 130, pvRead: 127.7, zChokePct: 60, wChoke: 58, dpChoke: 6.646, qPumpM3h: 90, nPumpPct: 85, pumpHeadM: 496.5, qCompM3s: 0.31, nCompPct: 100, compHeadKj: 94.83, tTracerH: 5.4, uLiq: 1.2, tracerC: 1.011, qInjM3d: 50, mWaterKgs: 3, doseMeas: 11.53, gelLenM: 8000, restartP: 66.45 }, { tShutH: 24, t0: 45, tCold: 10.18, tBlowMin: 180, p0: 75, pBlow: 9.515, vGas: 2.6, pigArrival: 2.835, pigDp: 1.352, tValveS: 30, valvePos: 0.8839, tStepS: 1000, stepPv: -0.6939, tClS: 1000, clPv: 0.9528, pTrue: 95, pvRead: 94.6, zChokePct: 70, wChoke: 62, dpChoke: 4.882, qPumpM3h: 120, nPumpPct: 85, pumpHeadM: 481.4, qCompM3s: 0.18, nCompPct: 85, compHeadKj: 85.78, tTracerH: 5.8, uLiq: 1.2, tracerC: 1.001, qInjM3d: 60, mWaterKgs: 3.5, doseMeas: 11.71, gelLenM: 10000, restartP: 75.72 }, { tShutH: 30, t0: 45, tCold: 7.973, tBlowMin: 240, p0: 75, pBlow: 4.815, vGas: 3, pigArrival: 2.389, pigDp: 1.324, tValveS: 40, valvePos: 0.9738, tStepS: 1400, stepPv: -0.7452, tClS: 1300, clPv: 0.9507, pTrue: 85, pvRead: 86.49, zChokePct: 85, wChoke: 66, dpChoke: 3.215, qPumpM3h: 150, nPumpPct: 85, pumpHeadM: 440.6, qCompM3s: 0.2, nCompPct: 85, compHeadKj: 82.55, tTracerH: 6.3, uLiq: 1.2, tracerC: 0.983, qInjM3d: 70, mWaterKgs: 3.5, doseMeas: 13.51, gelLenM: 12000, restartP: 82.14 }, { tShutH: 14, t0: 50, tCold: 18.76, tBlowMin: 75, p0: 80, pBlow: 34.1, vGas: 2.3, pigArrival: 3.337, pigDp: 1.312, tValveS: 50, valvePos: 0.9902, tStepS: 1800, stepPv: -0.803, tClS: 1700, clPv: 0.9776, pTrue: 75, pvRead: 75.24, zChokePct: 100, wChoke: 70, dpChoke: 2.294, qPumpM3h: 180, nPumpPct: 90, pumpHeadM: 469, qCompM3s: 0.22, nCompPct: 90, compHeadKj: 90.5, tTracerH: 7, uLiq: 1.2, tracerC: 1.009, qInjM3d: 80, mWaterKgs: 4, doseMeas: 13.52, gelLenM: 15000, restartP: 95.95 }, { tShutH: 5, t0: 50, tCold: 34.64, tBlowMin: 25, p0: 80, pBlow: 59.09, vGas: 1.7, pigArrival: 5.325, pigDp: 1.328, tValveS: 12, valvePos: 0.524, tStepS: 250, stepPv: -0.2771, tClS: 2200, clPv: 0.9815, pTrue: 105, pvRead: 105.6, zChokePct: 45, wChoke: 48, dpChoke: 10.02, qPumpM3h: 60, nPumpPct: 80, pumpHeadM: 460.3, qCompM3s: 0.24, nCompPct: 95, compHeadKj: 98.3, tTracerH: 4.5, uLiq: 1.2, tracerC: 0.1175, qInjM3d: 45, mWaterKgs: 2.8, doseMeas: 11.02, gelLenM: 5000, restartP: 55.15 }];
const CAL_VALID = [{ tShutH: 3, t0: 50, tCold: 38.71, tBlowMin: 15, p0: 80, pBlow: 65.21, vGas: 1.8, pigArrival: 4.874, pigDp: 1.33, tValveS: 5, valvePos: 0.1588, tStepS: 90, stepPv: -0.05242, tClS: 150, clPv: 1.178, pTrue: 65, pvRead: 66.87, zChokePct: 18, wChoke: 32, dpChoke: 58.09, qPumpM3h: 50, nPumpPct: 100, pumpHeadM: 728.9, qCompM3s: 0.2, nCompPct: 100, compHeadKj: 122.1, tTracerH: 3.6, uLiq: 1.2, tracerC: 0.002727, qInjM3d: 12, mWaterKgs: 1.8, doseMeas: 4.899, gelLenM: 800, restartP: 36.56 }, { tShutH: 7, t0: 50, tCold: 30.17, tBlowMin: 40, p0: 80, pBlow: 48.93, vGas: 2.4, pigArrival: 3.195, pigDp: 1.296, tValveS: 9, valvePos: 0.3806, tStepS: 160, stepPv: -0.16, tClS: 350, clPv: 0.9565, pTrue: 88, pvRead: 89.3, zChokePct: 28, wChoke: 40, dpChoke: 26.48, qPumpM3h: 110, nPumpPct: 100, pumpHeadM: 683.5, qCompM3s: 0.24, nCompPct: 100, compHeadKj: 111, tTracerH: 4.2, uLiq: 1.2, tracerC: 0.002613, qInjM3d: 22, mWaterKgs: 2.2, doseMeas: 7.265, gelLenM: 2500, restartP: 44.15 }, { tShutH: 11, t0: 50, tCold: 22.64, tBlowMin: 75, p0: 80, pBlow: 34.02, vGas: 2.9, pigArrival: 2.455, pigDp: 1.314, tValveS: 16, valvePos: 0.6757, tStepS: 380, stepPv: -0.4093, tClS: 550, clPv: 0.9519, pTrue: 99, pvRead: 99.36, zChokePct: 38, wChoke: 44, dpChoke: 13.68, qPumpM3h: 140, nPumpPct: 90, pumpHeadM: 526.4, qCompM3s: 0.28, nCompPct: 100, compHeadKj: 104.8, tTracerH: 4.6, uLiq: 1.2, tracerC: 0.2693, qInjM3d: 33, mWaterKgs: 2.6, doseMeas: 8.988, gelLenM: 4500, restartP: 52.39 }, { tShutH: 15, t0: 46, tCold: 16.53, tBlowMin: 110, p0: 76, pBlow: 21.56, vGas: 3.2, pigArrival: 2.183, pigDp: 1.314, tValveS: 22, valvePos: 0.8008, tStepS: 700, stepPv: -0.61, tClS: 900, clPv: 0.94, pTrue: 115, pvRead: 115.1, zChokePct: 55, wChoke: 56, dpChoke: 8.135, qPumpM3h: 170, nPumpPct: 90, pumpHeadM: 482.6, qCompM3s: 0.19, nCompPct: 90, compHeadKj: 96.49, tTracerH: 5.2, uLiq: 1.2, tracerC: 0.9922, qInjM3d: 44, mWaterKgs: 3.2, doseMeas: 9.644, gelLenM: 7000, restartP: 62.34 }, { tShutH: 22, t0: 46, tCold: 11.42, tBlowMin: 160, p0: 76, pBlow: 12.14, vGas: 2.3, pigArrival: 3.281, pigDp: 1.326, tValveS: 35, valvePos: 0.9321, tStepS: 1200, stepPv: -0.7285, tClS: 1500, clPv: 0.957, pTrue: 125, pvRead: 129.8, zChokePct: 75, wChoke: 64, dpChoke: 4.169, qPumpM3h: 200, nPumpPct: 95, pumpHeadM: 512.6, qCompM3s: 0.21, nCompPct: 90, compHeadKj: 91.5, tTracerH: 6, uLiq: 1.2, tracerC: 0.9932, qInjM3d: 55, mWaterKgs: 3.6, doseMeas: 10.86, gelLenM: 9000, restartP: 69.81 }, { tShutH: 28, t0: 46, tCold: 8.853, tBlowMin: 210, p0: 76, pBlow: 6.897, vGas: 2.7, pigArrival: 2.705, pigDp: 1.302, tValveS: 45, valvePos: 0.9591, tStepS: 1600, stepPv: -0.7771, tClS: 2000, clPv: 0.9688, pTrue: 72, pvRead: 72.85, zChokePct: 90, wChoke: 68, dpChoke: 2.934, qPumpM3h: 80, nPumpPct: 85, pumpHeadM: 512, qCompM3s: 0.26, nCompPct: 95, compHeadKj: 95.29, tTracerH: 6.6, uLiq: 1.2, tracerC: 1.014, qInjM3d: 66, mWaterKgs: 3.9, doseMeas: 11.6, gelLenM: 13000, restartP: 85.91 }];

const suite = {
  id: 'ops', num: 5, title: 'Operations, Control & Flow-Assurance Management', short: 'Operations', icon: '🎛️',
  tagline: 'Shutdown, cooldown, restart, blowdown, pigging, chemicals, slug control, pumps and compressors, operating logic and the operating envelope in one study.',
  description: 'Starting from the steady flow picture of the case, the suite solves the cooldown of the line by radial finite-volume conduction, the settle-out, restart pressure, warm-up and ramp-up surge, the blowdown of the gas inventory, a pig run with the liquid it pushes, inhibitor dosing with the advection–dispersion of the front, and heating. A four-state riser model gives the choke opening where slugging starts; it is stabilised by PID (with gain scheduling, robust and self-tuned settings), by linear and nonlinear model-predictive control, with extended and unscented Kalman filters estimating the riser inventory. A variable-speed export pump and a gas compressor with anti-surge control are tested for turndown and trip. The rate window, an optimised operating point, a timed shutdown–restart sequence with alarms, the chemical stock, solids triggers and a comparison with operating records complete the study. Literature constants are listed with their sources, and published benchmarks and experiments are compared on the calibration tab.',
  guide: ['Run the fluid, network and flow suites first (or use the reference case) and link their values on the Inputs tab.', 'Set the shutdown duration, the preservation strategy and the restart ramp; check the cooldown and no-touch times.', 'Size the blowdown orifice and check the minimum temperatures and the seabed pressure left by the liquid head.', 'Review pigging surge against the slug-catcher and the inhibitor dose, volumes and front travel time.', 'On Model setup choose the controller and tuning; compare PID with linear and nonlinear MPC and read the critical choke opening; set the pump and compressor data (or maps) and read the surge margin in the turndown test.', 'Enter the chemical stock, the deposit state from the solids suite and its limits; paste commissioning, shutdown, blowdown, pigging or transient records into the operating-record table to see the model beside them.', 'Read the operating window, the optimised operating point and the event sequence; attach a historian export as the operating log to replay it.'],
  implemented: ['transient mass balance', 'momentum balance', 'energy balance', 'component balances', 'equipment inventory', 'tank/separator level', 'valve actuator dynamics', 'rotating-equipment dynamics', 'pid control', 'pi control', 'p control', 'feedforward control', 'feedback control', 'cascade control', 'ratio control', 'override/selective control', 'anti-windup', 'model predictive control', 'nonlinear mpc', 'adaptive control', 'robust control', 'optimal control', 'state-space model', 'transfer-function model', 'state observer', 'kalman filter', 'extended kalman filter', 'unscented kalman filter',
    'transient heat equation', 'fourier conduction', 'newton cooling', 'lumped-capacitance', 'multilayer cylindrical conduction', 'transient fluid energy equation', 'transient mass/energy balances', 'compressible-flow', 'critical-flow', 'homogeneous-equilibrium', 'homogeneous-relaxation model', 'joule–thomson cooling', 'pig force balance', 'differential-pressure equation', 'friction/contact-force', 'pig velocity equation', 'liquid inventory/displacement balance', 'bypass/leakage',
    'advection–diffusion equation', 'species conservation', 'mixing equations', 'partitioning models', 'inhibitor mass balance', 'linear programming', 'nonlinear programming', 'mixed-integer linear programming', 'mixed-integer nonlinear programming', 'dynamic optimization', 'sequential quadratic programming', 'interior-point', 'genetic algorithm', 'particle-swarm', 'bayesian optimization',
    'transient multiphase solver + mpc', 'hydrate-risk model + inhibitor optimizer', 'digital twin + state estimation', 'physics model + data-driven surrogate', 'mechanistic model + machine-learning residual correction', 'shutdown model + optimization', 'production optimization + flow-assurance constraints',
    'production rates', 'pressure and temperature', 'fluid inventories', 'current valve and choke openings', 'pump and compressor states and speeds', 'separator conditions', 'chemical inventories', 'inhibitor injection rates', 'heating status', 'pig location if a pig is already in the system', 'current alarm states', 'controller set points and operating mode', 'pressure', 'temperature', 'phase distribution', 'liquid accumulation',
    'production targets and permitted production ranges', 'startup and ramp-up profiles', 'shutdown sequences', 'restart schedules', 'minimum and maximum valve openings', 'choke limits', 'pump and compressor limits', 'separator pressure and level constraints', 'depressurization and blowdown schedules', 'inhibitor injection limits', 'heating limits', 'pig-launch and pig-receive conditions', 'emergency shutdown actions', 'set points', 'allowable deviations', 'actuator limits', 'operating constraints', 'alarm thresholds and safety limits', 'minimum acceptable hydrate safety margin', 'maximum shutdown duration', 'acceptable liquid surge', 'maximum solids accumulation and any conditions that trigger chemical injection', 'heating', 'depressurization or production-rate adjustment',
    'current state from modules 1-4', 'production targets and schedules', 'startup/shutdown/restart/ramp profiles', 'valve/choke commands and limits', 'depressurization/blowdown schedules', 'inhibitor/chemical injection rates, locations and capacity limits', 'pig geometry/friction/bypass and launch/receive schedule', 'heating/cooling strategy', 'pump/compressor controls and operating maps', 'separator/control settings', 'controller set points/gains', 'sensor states', 'alarm/interlock thresholds', 'operational and safety constraints',
    'startup/shutdown/restart trajectories', 'cooldown/warm-up time', 'restart pressure/rate requirements', 'depressurization/blowdown time and minimum temperature', 'transient liquid/solid inventory', 'pig position/velocity and generated liquid surge', 'inhibitor concentration/distribution and required dosage', 'hydrate/solids safety margin', 'equipment/control response', 'feasible operating envelope', 'maximum allowable shutdown duration', 'safe restart window', 'constraint violations, alarms/interlocks and recommended/optimized operating actions',
    // calibration: each quantity is a fitted parameter of the calibration model with a measured target that responds to it
    'valve actuator response', 'valve cv', 'choke characteristics', 'pump curves', 'compressor maps', 'controller gains', 'pid parameters', 'sensor dynamics', 'sensor bias', 'process dead time', 'actuator dead time', 'thermal time constants', 'inhibitor mixing/dispersion parameters', 'chemical-injection efficiency', 'pig friction', 'pig bypass', 'pig velocity parameters', 'cooldown parameters', 'restart friction/pressure parameters',
    'steady-to-transient initialization', 'event scheduler verification', 'valve opening/closing logic', 'controller logic', 'interlock logic', 'alarm logic', 'constraint handling', 'pid benchmark tests', 'pump/compressor control tests', 'mass balance during switching events', 'energy balance during shutdown', 'depressurization conservation', 'restart conservation', 'pig-tracking conservation', 'event-time accuracy', 'time-step independence', 'state-machine tests', 'fault-handling tests', 'fail-safe tests', 'optimization convergence', 'operating-envelope boundary verification',
    // validation: a tick means that this comparison is supported in the app — a sourced reference data set (blowdown) or a record type of the
    // operating-record table / operating log with the model prediction and error metrics beside it; it does not claim that field data are shipped
    'commissioning data', 'startup records', 'shutdown records', 'restart records', 'emergency shutdowns', 'blowdown/depressurization tests', 'cooldown measurements', 'pigging records', 'pig arrival times', 'liquid-surging measurements', 'chemical tracer measurements', 'meg/methanol concentration measurements', 'valve-response data', 'compressor transient data', 'pump transient data', 'separator-level histories', 'scada/historian data', 'field alarm/event histories'],
  referenceOnly: [],
  equationsNote: 'Screening-level operations models. Cooldown: radial conduction with a lumped fluid node at every axial station (no axial conduction, no natural-circulation redistribution of heat; liquid settles into the low points between crests). Settle-out and cooling pressure assume a fixed gas mass with no inter-phase mass transfer. Blowdown treats the line gas as one lumped volume with a lumped wall-and-liquid heat sink (no axial pressure gradient, so long lines blow down somewhat slower than predicted); against a published nitrogen vessel experiment the model empties the vessel too fast (see the reference data sets). Two-phase discharge is the ω-method homogeneous-equilibrium model with a fixed inlet liquid fraction; the homogeneous-relaxation option delays the liberation of dissolved gas with a relaxation time (dm/dt = (m_eq − m)/Θ, Downar-Zapolski correlation or a set value) — it is a lumped relaxation closure, not a one-dimensional nozzle solution. After blowdown, gas trapped behind liquid seals is assumed to keep the head of the liquid legs towards the nearest open end (upper bound). Pigging is quasi-steady on the steady profile. The slugging model is a four-state riser model tuned to the steady solution: it reproduces the onset and period of riser-induced slugging and its response to the choke, not hydrodynamic slug statistics. Robust control here means a multi-model PI design with a sensitivity-peak limit and a loop-shaping check (Ms, Mt); no H∞ synthesis or μ-analysis is performed. Adaptive control is gain scheduling on the choke opening and a recursive-least-squares self-tuner. Nonlinear MPC is single shooting on the riser model with full state feedback. Pump and compressor are quasi-steady machines on fan-law maps with rotor inertia, a torque-limited driver and lumped suction and discharge volumes; surge is counted when the demanded head exceeds the map, the surge cycle itself is not resolved. On the Equations tab a tick under calibration means the quantity can be estimated on the calibration tab; a tick under validation means that the comparison is supported in the app (sourced reference data for blowdown, control, optimisation and dispersion benchmarks; for field records the operating-record table and the operating log put the model prediction and error metrics beside your data) — no field data are shipped with the app.',
  inputs: INPUTS,
  presets: [
    { name: 'Planned shutdown and restart (reference tie-back)', values: { tShut: 24, preserve: 'inhibit', rampHours: 6, qStartPct: 30 } },
    { name: 'Emergency shutdown with blowdown', values: { tShut: 96, tHorizon: 96, preserve: 'blowdown', tBlowStart: 2, tDecision: 1, orificeMm: 45, pBack: 1.5, unplannedPerYear: 10, unplannedHours: 36, ntCool: 192 } },
    { name: 'Pigging campaign for wax management', values: { wat: 46, waxByPigging: true, pigRatePct: 90, pigInterval: 7, waxThk: 4, pigBypass: 4, pigFric: 0.45, pigLeak: 6, pigEff: 85, preserve: 'none', tShut: 6 } },
    { name: 'Severe-slugging control at low rate', values: { rateFrac: 40, chokePct: 0, ctlMode: 'PI', tuning: 'auto', cascade: true, tCtl: 24, dtCtl: 60, severeSlugging: true, qLoPct: 20, qHiPct: 120 } },
    { name: 'Continuous MEG injection with regeneration loop', values: { inhibitor: 'MEG', leanWt: 90, dosingBasis: 'max', inhPrice: 1100, megStorageDays: 5, megLossPct: 0.5, pumpMax: 40, preserve: 'inhibit' } },
    { name: 'Gas compression with a suction cooler and deep turndown', values: { coolerOutT: 25, compTurndownPct: 35, compSmCtl: 12, compPd: 75, rotNminPct: 75, pumpTurndownPct: 40, tuning: 'robust', msMax: 1.8 } },
    { name: 'Restart with a pig left in the line and a wax trigger', values: { pigX0: 9000, pigSlug0: 25, depWaxMm: 4.5, depRateMmD: 0.25, wat: 44, waxByPigging: true, preserve: 'inhibit', bdRoute: 'top', bdMode: 'hrm', tuning: 'rls', ctlAuto: false } },
    { name: 'Late-life turndown envelope', values: { rateFrac: 50, qLoPct: 15, qHiPct: 110, nEnv: 10, wat: 34, slugControl: false, pAvail: 140, plannedPerYear: 3, unplannedPerYear: 12, preserve: 'heat' } },
  ],
  pull: ({ fluid, outputs } = {}) => {
    const n = outputs?.net || {}, f = outputs?.flow || {}, s = outputs?.solids || {}, p = outputs?.pvt || {}, items = [], offer = (key, value, from, ok = isNum(value)) => { if (ok) items.push({ key, value, from }); };
    offer('sepP', fluid?.Pout, 'Case fluid: arrival pressure', isNum(fluid?.Pout) && fluid.Pout >= 1.2);
    offer('inhibitor', fluid?.inhibitor, 'Case fluid: inhibitor', INH_OPTS.includes(fluid?.inhibitor));
    offer('wat', p.wat, 'Fluid suite: wax appearance temperature');
    offer('idMm', n.id * 1000, 'Network suite: inner diameter', isNum(n.id) && n.id > 0.025);
    offer('wtMm', n.wt * 1000, 'Network suite: wall thickness', isNum(n.wt) && n.wt > 0.001);
    offer('uValue', n.uValue, 'Network suite: U-value', isNum(n.uValue) && n.uValue > 0);
    if (Array.isArray(n.layers) && n.layers.some((l) => l && l.t > 0 && l.k > 0)) items.push({ key: 'layers', value: n.layers.filter((l) => l && l.t > 0 && l.k > 0).map((l) => ({ name: String(l.name || 'Layer'), t: +(l.t * 1000).toFixed(2), k: l.k, rho: isNum(l.rho) ? l.rho : l.k > 5 ? 7850 : l.k > 0.6 ? 2400 : 900, cp: isNum(l.cp) ? l.cp : l.k > 5 ? 470 : l.k > 0.6 ? 900 : 1700 })), from: 'Network suite: coating layers' });
    offer('thermalMass', n.thermalMass / 1000, 'Network suite: thermal mass (J/m/K)', isNum(n.thermalMass) && n.thermalMass > 1e4 && n.thermalMass < 5e6);
    offer('lineVolume', n.volume, 'Network suite: line volume', isNum(n.volume) && n.volume > 0);
    offer('slugCatcherVol', n.slugCatcherVol, 'Network suite: slug-catcher volume', isNum(n.slugCatcherVol) && n.slugCatcherVol > 1);
    offer('sepP', n.separatorP, 'Network suite: separator pressure', isNum(n.separatorP) && n.separatorP >= 1.2);
    // the network suite's chokeCv is the Cv the WELLHEAD choke needs at the case rate, not the rated Cv of the topside choke used here: it is not linked
    { const dp = s.depositProfile, mx = (a) => (Array.isArray(a) && a.length ? Math.max(...a.filter(isNum), 0) * 1000 : null);
      if (dp) { offer('depWaxMm', mx(dp.wax), 'Solids suite: largest wax deposit'); offer('depScaleMm', mx(dp.scale), 'Solids suite: largest scale deposit'); offer('depHydMm', mx(dp.hydrate), 'Solids suite: largest hydrate deposit'); } }
    offer('hydFrac', s.hydrateFraction, 'Solids suite: hydrate fraction in the liquid', isNum(s.hydrateFraction) && s.hydrateFraction >= 0 && s.hydrateFraction <= 1);
    offer('blockagePct', 100 * s.blockage, 'Solids suite: bore area lost', isNum(s.blockage) && s.blockage >= 0 && s.blockage <= 1);
    offer('depRateMmD', s.waxRate, 'Solids suite: wax deposition rate', isNum(s.waxRate) && s.waxRate >= 0 && s.waxRate <= 50);
    offer('burial', n.burial, 'Network suite: burial depth', isNum(n.burial) && n.burial >= 0);
    offer('tSeabed', n.tSeabed, 'Network suite: seabed temperature'); offer('tSurface', n.tSeaSurface, 'Network suite: sea-surface temperature');
    offer('slugSurge', f.slug?.surge, 'Flow suite: slug surge volume', isNum(f.slug?.surge) && f.slug.surge >= 0);
    if (typeof f.severeSlugging === 'boolean') items.push({ key: 'severeSlugging', value: f.severeSlugging, from: 'Flow suite: severe slugging' });
    offer('turndown', f.turndownRate, 'Flow suite: minimum stable rate', isNum(f.turndownRate) && f.turndownRate >= 0 && f.turndownRate <= 1);
    offer('inhRequired', s.inhibitorRequired, 'Solids suite: inhibitor required', isNum(s.inhibitorRequired) && s.inhibitorRequired >= 0 && s.inhibitorRequired <= 90);
    offer('pigInterval', s.piggingInterval, 'Solids suite: pigging interval', isNum(s.piggingInterval) && s.piggingInterval >= 0);
    offer('plugTime', s.plugTime, 'Solids suite: time to plug', isNum(s.plugTime) && s.plugTime >= 0);
    offer('wat', s.wat, 'Solids suite: wax appearance temperature', isNum(s.wat) && !isNum(p.wat));
    return items;
  },
  site: (site) => { const d = site?.data || {}; return [['tSeabed', d.seabedTemp, 'Seabed temperature at site'], ['tSurface', d.sst, 'Sea-surface temperature at site'], ['tAir', d.airTemp, 'Air temperature at site'], ['currentSpeed', d.currentSpeed, 'Current speed at site'], ['elecPrice', d.electricityPrice, 'Electricity price at site']].filter((x) => isNum(x[1])).map(([key, value, from]) => ({ key, value, from })); },
  run,
  mesh: [
    { name: 'Cooldown solver: radial cells and time steps', keys: ['nr', 'ntCool'], min: 4, note: 'Radial finite-volume cells through wall and coatings and implicit time steps over the simulated cooldown.', metrics: [{ label: 'Cooldown time', unit: 'h', get: (r) => r.outputs.cooldownTime }, { label: 'Cold-spot temperature after 12 h', unit: '°C', get: (r) => r.outputs.coldSpotT12 }] },
    { name: 'Blowdown: time steps', keys: ['ntBlow'], min: 40, note: 'Second-order (Heun) steps from the opening of the valve to the target pressure.', metrics: [{ label: 'Blowdown time', unit: 'h', get: (r) => r.outputs.blowdownTime }, { label: 'Minimum temperature', unit: '°C', get: (r) => r.outputs.blowdownMinT }] },
    { name: 'Control loop: sample and integration step', keys: ['dtCtl'], refine: 'divide', min: 2, note: 'The controller sample time is also the step of the L-stable Rosenbrock integration of the riser and separator model.', metrics: [{ label: 'Mean inlet pressure under control', unit: 'bara', get: (r) => r.outputs.ctlMeanP }, { label: 'Integral absolute error', unit: 'bar·h', get: (r) => r.outputs.ctlIae }] },
  ],
  calibration: {
    note: 'Estimates the uncertain parameters of the operations models from records: U-value and thermal-mass multipliers (cooldown), blowdown discharge coefficient (equivalent valve Cv), pig seal friction, bypass area and bypass coefficient, valve actuator time constant and dead time, process gain, time constant and dead time, transmitter time constant and bias, controller gain and integral time (from a recorded closed-loop response), choke Cv and characteristic exponent, pump shut-off ratio and rated differential pressure, compressor efficiency (map head), dispersion multiplier, injection efficiency and gel yield stress. Each row may hold any subset of the measurements; untick the parameters your data cannot identify (bypass area and coefficient act as a product; U-value and thermal mass need records at several times).',
    params: [{ key: 'uMult', label: 'U-value multiplier (cooldown)', lo: 0.5, hi: 2 }, { key: 'cMult', label: 'Thermal-mass multiplier (thermal time constant)', lo: 0.5, hi: 2 }, { key: 'cdBlow', label: 'Blowdown discharge coefficient (valve Cv)', lo: 0.3, hi: 1 }, { key: 'pigFric', label: 'Pig seal friction coefficient', lo: 0.05, hi: 1.5 }, { key: 'pigBypass', label: 'Pig bypass area (% of bore)', lo: 0.1, hi: 15 }, { key: 'pigCd', label: 'Pig bypass discharge coefficient (velocity parameter)', lo: 0.4, hi: 1 },
      { key: 'tauValve', label: 'Valve actuator time constant (s)', lo: 0.5, hi: 120 }, { key: 'actDead', label: 'Actuator dead time (s)', lo: 0, hi: 60 }, { key: 'procK', label: 'Process gain (bar per % opening)', lo: -5, hi: -0.02 }, { key: 'procTau', label: 'Process time constant (s)', lo: 30, hi: 5000 }, { key: 'deadTime', label: 'Process dead time (s)', lo: 0, hi: 300 }, { key: 'tauSensor', label: 'Transmitter time constant (sensor dynamics, s)', lo: 0, hi: 120 }, { key: 'sensorBias', label: 'Transmitter bias (bar)', lo: -10, hi: 10 },
      { key: 'kcMan', label: 'Controller gain (% opening per bar)', lo: -100, hi: -0.5 }, { key: 'tiMan', label: 'Controller integral time (s)', lo: 30, hi: 20000 }, { key: 'chokeCv', label: 'Choke Cv (fully open)', lo: 50, hi: 3000 }, { key: 'chokeExp', label: 'Choke characteristic exponent', lo: 0.4, hi: 3.5 }, { key: 'pumpShutoff', label: 'Pump curve: shut-off head ratio', lo: 1.03, hi: 1.9 }, { key: 'pumpDp', label: 'Pump curve: rated differential pressure (bar)', lo: 5, hi: 200 },
      { key: 'compEta', label: 'Compressor map: polytropic efficiency (%)', lo: 45, hi: 92 }, { key: 'dispMult', label: 'Dispersion multiplier (inhibitor mixing)', lo: 1, hi: 5000 }, { key: 'injEff', label: 'Chemical-injection efficiency (%)', lo: 20, hi: 100 }, { key: 'yieldStress', label: 'Gel yield stress (restart pressure, Pa)', lo: 0, hi: 500 }],
    columns: [{ key: 'tShutH', label: 'Time since shut-in', unit: 'h' }, { key: 't0', label: 'Temperature at shut-in', unit: '°C' }, { key: 'tCold', label: 'Fluid temperature', unit: '°C' }, { key: 'tBlowMin', label: 'Time since blowdown valve opened', unit: 'min' }, { key: 'p0', label: 'Pressure at opening', unit: 'bara' }, { key: 'pBlow', label: 'Line pressure', unit: 'bara' }, { key: 'vGas', label: 'Mixture velocity behind the pig', unit: 'm/s' }, { key: 'pigArrival', label: 'Pig arrival time', unit: 'h' }, { key: 'pigDp', label: 'Pig differential pressure', unit: 'bar' },
      { key: 'tValveS', label: 'Time since valve step command', unit: 's' }, { key: 'valvePos', label: 'Valve travel', unit: 'fraction' }, { key: 'tStepS', label: 'Time since open-loop step of 1 % opening', unit: 's' }, { key: 'stepPv', label: 'Pressure reading change', unit: 'bar' }, { key: 'tClS', label: 'Time since 1 bar set-point step (closed loop)', unit: 's' }, { key: 'clPv', label: 'Pressure change in closed loop', unit: 'bar' }, { key: 'pTrue', label: 'Reference pressure (test gauge)', unit: 'bara' }, { key: 'pvRead', label: 'Transmitter reading', unit: 'bara' },
      { key: 'zChokePct', label: 'Choke opening', unit: '%' }, { key: 'wChoke', label: 'Mass rate through the choke', unit: 'kg/s' }, { key: 'dpChoke', label: 'Choke pressure drop', unit: 'bar' }, { key: 'qPumpM3h', label: 'Pump flow', unit: 'm³/h' }, { key: 'nPumpPct', label: 'Pump speed', unit: '%' }, { key: 'pumpHeadM', label: 'Pump head', unit: 'm' }, { key: 'qCompM3s', label: 'Compressor inlet flow', unit: 'm³/s' }, { key: 'nCompPct', label: 'Compressor speed', unit: '%' }, { key: 'compHeadKj', label: 'Compressor polytropic head', unit: 'kJ/kg' },
      { key: 'tTracerH', label: 'Time since tracer injection', unit: 'h' }, { key: 'uLiq', label: 'Liquid velocity', unit: 'm/s' }, { key: 'tracerC', label: 'Tracer at the outlet', unit: 'C/C0' }, { key: 'qInjM3d', label: 'Inhibitor injection rate', unit: 'm³/d' }, { key: 'mWaterKgs', label: 'Water rate', unit: 'kg/s' }, { key: 'doseMeas', label: 'Inhibitor in the water at arrival', unit: 'wt %' }, { key: 'gelLenM', label: 'Gelled length', unit: 'm' }, { key: 'restartP', label: 'Restart pressure', unit: 'bara' }],
    targets: [{ key: 'tCold', label: 'Fluid temperature', unit: '°C' }, { key: 'pBlow', label: 'Line pressure', unit: 'bara' }, { key: 'pigArrival', label: 'Pig arrival time', unit: 'h' }, { key: 'pigDp', label: 'Pig differential pressure', unit: 'bar' }, { key: 'valvePos', label: 'Valve travel', unit: 'fraction' }, { key: 'stepPv', label: 'Open-loop pressure change', unit: 'bar' }, { key: 'clPv', label: 'Closed-loop pressure change', unit: 'bar' }, { key: 'pvRead', label: 'Transmitter reading', unit: 'bara' }, { key: 'dpChoke', label: 'Choke pressure drop', unit: 'bar' },
      { key: 'pumpHeadM', label: 'Pump head', unit: 'm' }, { key: 'compHeadKj', label: 'Compressor polytropic head', unit: 'kJ/kg' }, { key: 'tracerC', label: 'Tracer at the outlet', unit: 'C/C0' }, { key: 'doseMeas', label: 'Inhibitor in the water', unit: 'wt %' }, { key: 'restartP', label: 'Restart pressure', unit: 'bara' }],
    model: calModel,
    sample: CAL_SAMPLE,
    validationSample: CAL_VALID,
  },
  verify,
  validationData: VALIDATION,
  live: { key: 'log', label: 'Operating log', help: 'Follow a SCADA / historian export (CSV with time in hours, rate in % of the case rate, choke opening in %, inlet pressure in bara and arrival temperature in °C); the replay and the residual correction update as rows arrive.' },
};
export default suite;
