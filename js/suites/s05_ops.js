// Suite 5 — Operations, Control & Flow-Assurance Management.
// Shutdown/cooldown (radial finite-volume conduction + lumped check), restart and ramp-up, blowdown, pigging,
// chemical injection, heating, riser-slugging control (PID / MPC / Kalman), operating logic (scheduler, alarms,
// state machine), operating envelope and optimisation, surrogate / residual correction and historical replay.
// SI inside the engines; bara, °C, h and mm at the interfaces.
import { clamp, linspace, interp1, brent, tridiag, solveLinear, rk45, nelderMead, lstsq, rng, metrics, mean, sum, isNum, fmt } from '../core/num.js';
import { fluidModel, inhibitorFor, hydrateDepression, INHIBITORS, R, VM_STD } from '../core/thermo.js';
import { G, uValue, hInside, hOutside, frictionFactor, marchSteady } from '../core/pipe.js';
import { flowPicture, caseLine, ambientAt } from '../core/caseflow.js';
import { BASE } from '../data/basecase.js';

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
 * Returns { pSettle (bara), holdup[] (settled), gasMass (kg), liquidVol (m³), gasVol (m³), headUphill (bar: liquid legs the restart must lift), levelRiser (m of liquid column at the line end) }.
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
  for (let i = n - 1; i >= 0 && st.z[i] > st.z[Math.max(i - 1, 0)] - 1e-9; i--) lev += H[i] * Math.max((st.dz || [])[i] || 0, 0);
  return { pSettle, holdup: H, gasMass: mG, gasVol, liquidVol: sum(st.holdup) * vol, headUphill: head / 1e5, rhoL, levelRiser: lev };
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
  const w = Math.max(omega, 1e-6), fc = (e) => e * e + (w * w - 2 * w) * (1 - e) ** 2 + 2 * w * w * Math.log(e) + 2 * w * w * (1 - e), etaC = fc(1e-6) * fc(0.999999) < 0 ? brent(fc, 1e-6, 0.999999, 1e-12) : 0.55 + 0.217 * Math.log(w) - 0.046 * Math.log(w) ** 2, eta = pBack / P0;
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
 *      hem: { mLiq (kg), rhoL, frac (liquid volume fraction at the valve inlet) } for two-phase discharge }
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
  const der = (s) => { // s = [m, T, Tw, mL]
    const Pn = (s[0] * Zf(P, s[1]) * R * s[1]) / (q.mw * Vg), [w, xg] = flow(Pn, s[1], s[3]), cv = (Zf(Pn, s[1]) * R) / (q.mw * (q.k - 1)), wg = w * xg, wl = w - wg;
    const FP = q.flash ? Math.max(0, q.flash(Pn, s[1])) * Pn : 0, fl = (wg * FP) / (Math.max(s[0], 1e-9) + FP); // liberation follows the actual pressure fall: dP/dt ≈ −P (wg − fl)/m
    const dT = q.mode === 'isothermal' ? 0 : (-(wg * Zf(Pn, s[1]) * R * s[1]) / q.mw + fl * cv * q.k * (s[2] - s[1])) / (Math.max(s[0], 1e-9) * cv);
    const dTw = q.mode === 'wall' && q.wallC > 0 ? (q.extUA * (q.tAmb - s[2]) - fl * q.latent) / q.wallC : 0;
    return { d: [-wg + fl, dT, dTw, -wl], w, wg, fl, P: Pn };
  };
  const out = { t: [0], P: [P], T: [T], Tw: [Tw], mdot: [w0], m: [m] }, m0 = m, mL0 = mL;
  let disc = 0, flashed = 0, tEnd = null, tMark = null, minT = T, minTw = Tw, minTd = T - q.jt * (P - q.pBack), peak = w0, time = 0;
  for (let k = 0; k < q.n * q.maxFactor && tEnd === null; k++) {
    const s0 = [m, T, Tw, mL], a = der(s0), s1 = s0.map((v, i) => v + dt * a.d[i]); s1[0] = Math.max(s1[0], 1e-9); s1[1] = Math.max(s1[1], 20); s1[3] = Math.max(s1[3], 0);
    const b = der(s1), Pold = P;
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
  return { ...out, dt, tEnd, tMark, minT, minTw, minTdown: minTd, peak, discharged: disc, flashed, liquidOut: mL0 - mL, mass: { initial: m0 + mL0, final: m + mL }, pFinal: P };
}

// ---- pigging --------------------------------------------------------------------------------------------------------
/**
 * Quasi-steady pig run along a line. o: { s[] (m, arc length at the nodes), z[], vm[] (mixture velocity behind the pig, m/s), holdup[], vsl[], rhoM[] (also the density of the fluid passing the bypass unless rhoBypass[] is given),
 *   rhoL, D, fric (Pa to keep the pig moving), mass (kg), bypass (bypass area / pipe area), cdBypass, leak (fraction of the swept liquid passing back through the pig),
 *   holdSlug (liquid fraction of the slug), qDrain (m³/s the receiving facility can process), fSlug (Darcy friction factor of the liquid slug) }
 * Force balance on the pig: Δp A = friction + weight component; bypass leakage of the local mixture through the pig as an orifice under Δp; the liquid overtaken
 * collects as a slug ahead. Returns { t[], x[] (s of the pig), v[], slug[] (m³), transit (s), tFront (s), vMean, swept, leaked, received (m³), duration (s of liquid arrival),
 *   surge (m³ above the drain capacity), dpPig (Pa, mean), dpExtra (Pa, largest added line pressure drop), stalled }.
 */
export function pigRun(o) {
  const q = { fric: 1e5, mass: 80, bypass: 0.02, cdBypass: 0.7, leak: 0, holdSlug: 0.95, qDrain: 0.05, fSlug: 0.02, ...o }, n = q.s.length, A = (Math.PI * q.D * q.D) / 4, L = q.s[n - 1];
  const out = { t: [0], x: [q.s[0]], v: [], slug: [0] };
  let t = 0, Vs = 0, swept = 0, leaked = 0, received = 0, tFront = null, dpSum = 0, dpExtra = 0, stalled = false;
  for (let i = 0; i < n - 1; i++) {
    const ds = q.s[i + 1] - q.s[i], sinT = clamp((q.z[i + 1] - q.z[i]) / ds, -1, 1), dp = Math.max(q.fric + (q.mass * G * sinT) / A, 0.05 * q.fric), vLeak = q.cdBypass * q.bypass * Math.sqrt((2 * dp) / Math.max((q.rhoBypass || q.rhoM)[i], 0.5));
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
  return { ...out, transit: stalled ? null : t, tFront, vMean: stalled ? 0 : (L - q.s[0]) / t, swept, leaked, received, duration, surge: Math.max(0, received - q.qDrain * duration), dpPig: dpSum / Math.max(L - q.s[0], 1e-9), dpExtra, stalled, inPipe: Vs };
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
const PSAT = { MeOH: (T) => 10 ** (8.08097 - 1582.271 / (T + 239.726)) * 133.322, EtOH: (T) => 10 ** (8.20417 - 1642.89 / (T + 230.3)) * 133.322, MEG: (T) => 10 ** (8.09083 - 2088.936 / (T + 203.454)) * 133.322, DEG: (T) => 10 ** (7.63666 - 1939.359 / (T + 162.714)) * 133.322, TEG: (T) => 10 ** (7.6302 - 2156.46 / (T + 165.92)) * 133.322 };
/**
 * Thermodynamic-inhibitor requirement with phase-partitioning losses. o: { dT (°C depression needed), inh ('MeOH' | 'MEG' | …), S (g/kg salinity),
 *   mWater (kg/s free water), lean (wt % purity of the injected chemical), P (bara), T (°C where the phases separate), qGasStd (Sm³/d), mOil (kg/s), mwOil (g/mol) }
 * Returns { wt (wt % in the aqueous phase), mAq (kg/s inhibitor in the water), lossGas, lossOil (kg/s), mTotal (kg/s pure), qInject (m³/d of the lean chemical), rich (kg/s aqueous phase returned) }.
 */
export function inhibitorDose(o) {
  const q = { dT: 0, inh: 'MeOH', S: 0, mWater: 1, lean: 100, P: 25, T: 20, qGasStd: 0, mOil: 0, mwOil: 150, ...o }, inh = INHIBITORS[q.inh] || INHIBITORS.MeOH, wt = q.dT > 0 ? inhibitorFor(q.dT, q.inh, q.S) : 0, w = wt / 100, lean = clamp(q.lean, 30, 100) / 100;
  if (w <= 0 || q.inh === 'none') return { wt: 0, mAq: 0, lossGas: 0, lossOil: 0, mTotal: 0, qInject: 0, rich: q.mWater, attainable: true };
  // water brought in by a lean (regenerated) chemical dilutes it: m_inh = w (mW + m_lean (1 − lean) + m_inh)
  const mAq = (w * q.mWater) / Math.max(1 - w / lean, 0.02), xAq = mAq / inh.MW / (mAq / inh.MW + (q.mWater + (mAq * (1 - lean)) / lean) / 18.015);
  const volatile = inh.rho < 1000, gamma = volatile ? 1.6 : 1, y = Math.min((gamma * xAq * (PSAT[q.inh] || PSAT.MEG)(q.T)) / (q.P * 1e5), 0.2), lossGas = (y * (q.qGasStd / DAY / VM_STD) * inh.MW) / 1000;
  const kHC = (volatile ? 6e-3 : 2e-4) * Math.exp(0.025 * (q.T - 20)), lossOil = (kHC * xAq * (q.mOil / (q.mwOil * 1e-3)) * inh.MW) / 1000, mTotal = mAq + lossGas + lossOil;
  return { wt, mAq, lossGas, lossOil, mTotal, qInject: ((mTotal / lean) / (lean * inh.rho + (1 - lean) * 1000)) * DAY, rich: q.mWater + mAq / lean, attainable: wt < 93.9 && w / lean < 0.98 };
}

// ---- riser slugging: low-order four-state model ------------------------------------------------------------
/**
 * Four-state pipeline–riser model (gas and liquid mass in the feed pipeline and in the riser) with a low-point
 * orifice pair and a topside choke.
 * p: { D, Lp (feed length), Vp (feed volume), Lr (riser height), Vr (riser + topside volume), theta (rad, feed inclination at the low point),
 *      rhoL, mwG (kg/mol), Z, Tp, Tr (K), muL, wG, wL (kg/s inflow), Ps (Pa separator), Cv (choke), kH, kL, kG (optional), aLp (feed liquid fraction), rough }
 * Returns { p, alg(y, z, wG, wL, Ps), f(z)(t, y), steady(z, wG, wL), linearise(z), poles(z), critical() }.
 * y = [mGp, mLp, mGr, mLr] (kg), z = choke opening 0–1.
 */
export function slugModel(p) {
  const o = { D: 0.254, Lp: 5000, Lr: 300, theta: 0.02, rhoL: 800, mwG: 0.02, Z: 0.9, Tp: 320, Tr: 310, muL: 2e-3, wG: 1, wL: 20, Ps: 25e5, Cv: 400, kH: 0.7, kL: 0.3, aLp: 0.4, rough: 4.5e-5, ...p };
  const r = o.D / 2, A = Math.PI * r * r, th = Math.max(Math.abs(o.theta), 1e-3), hc = (2 * r) / Math.cos(th), sinT = Math.sin(th);
  o.Vp = o.Vp || A * o.Lp; o.Vr = o.Vr || A * o.Lr * 1.15;
  const RTp = (o.Z * R * o.Tp) / o.mwG, RTr = (o.Z * R * o.Tr) / o.mwG, aL = clamp(o.aLp, 0.05, 0.9), hbar = o.kH * hc * aL, Kc = 2.403e-5 * o.Cv;
  const dmdh = (A * (1 - aL) * o.rhoL) / sinT, lamOf = (u, rho) => frictionFactor(Math.max((rho * u * o.D) / o.muL, 100), o.rough / o.D, 'haaland');
  const areas = (h) => { const hh = clamp(h, 0, hc), AG = A * ((hc - hh) / hc) ** 2; return [AG, A - AG]; };
  // friction factors are frozen at the nominal superficial velocities (weak function of the state, large saving in the Jacobians)
  const rhoGref = (1.6 * o.Ps) / RTr, u0 = o.wL / (A * o.rhoL), um0 = u0 + o.wG / (rhoGref * A), lamP = lamOf(u0, o.rhoL), lamR = lamOf(um0, 0.5 * o.rhoL);
  const fricP = (wL) => { const u = wL / (A * o.rhoL); return (aL * lamP * o.rhoL * u * u * o.Lp) / (2 * o.D); };
  const fricR = (wG, wL, rhoM, aLr) => { const um = wL / (A * o.rhoL) + wG / (rhoGref * A); return (aLr * lamR * rhoM * um * um * o.Lr) / (2 * o.D); };
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
    const aLt = clamp(2 * aLr - AL / A, 0, 1), rhoT = aLt * o.rhoL + (1 - aLt) * rhoGr, xL = (aLt * o.rhoL) / Math.max(rhoT, 1e-9), dPc = Prt - Ps, w = dPc > 0 ? Kc * z * Math.sqrt(rhoT * dPc) : 0;
    return { Pp, Prt, Prb, h: h / hc, aLr, aLt, rhoT, w, wGout: (1 - xL) * w, wLout: xL * w, wGlp, wLlp, d: [wGin - wGlp, wLin - wLlp, wGlp - (1 - xL) * w, wLlp - xL * w] };
  };
  /** Equilibrium for a choke opening (exists for every opening: the level and the top pressure are found by 1-D root finding). */
  const steady = (z, wG = o.wG, wL = o.wL, Ps = o.Ps) => {
    const w = wG + wL, xL = wL / w, top = (Prt) => { const rg = Prt / RTr, aLt = (xL * rg) / (o.rhoL - xL * (o.rhoL - rg)); return { rg, aLt, rhoT: aLt * o.rhoL + (1 - aLt) * rg }; };
    const gT = (Prt) => Kc * z * Math.sqrt(top(Prt).rhoT * (Prt - Ps)) - w;
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
  return { p: o, alg, f: (z) => (t, y) => alg(y, typeof z === 'function' ? z(t) : z).d, steady, linearise, poles, growth, critical, hc, A, dmdh, zFloor: clamp((o.wG + o.wL) / (Kc * Math.sqrt(o.rhoL * 300e5)), 0.02, 0.9) }; // zFloor: opening below which the choke alone would take more than about 300 bar
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
 * Controller settings from a first-order-plus-dead-time model: SIMC (Skogestad; τc = closed-loop time constant, default θ)
 * and Ziegler–Nichols from the ultimate gain and period of the same model. Returns rows { rule, mode, kc, ti, td } and { ku, pu }.
 */
export function tuningRules({ K, tau, theta }, tauC) {
  const th = Math.max(theta, 1e-9), tc = tauC ?? th, rows = [];
  rows.push({ rule: 'SIMC', mode: 'P', kc: tau / (K * (tc + th)), ti: 0, td: 0 });
  rows.push({ rule: 'SIMC', mode: 'PI', kc: tau / (K * (tc + th)), ti: Math.min(tau, 4 * (tc + th)), td: 0 });
  rows.push({ rule: 'SIMC', mode: 'PID', kc: (tau + th / 3) / (K * (tc + th)), ti: Math.min(tau + th / 3, 4 * (tc + th)), td: th / 3 });
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
  const n = y.length, g = 1 + Math.SQRT1_2, f0 = f(t, y), M = zeros(n, n);
  for (let j = 0; j < n; j++) { const e = 1e-7 * (Math.abs(y[j]) + 1e-3), yp = y.slice(); yp[j] += e; const fp = f(t, yp); for (let i = 0; i < n; i++) M[i][j] = (i === j ? 1 : 0) - (g * h * (fp[i] - f0[i])) / e; }
  // one LU factorisation (partial pivoting) serves both stages
  const A = M.map((r) => r.slice()), piv = [];
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
 *      sepV (m³), qDrain (m³/s), levelSp, levelHi, tauValve (s), np, nc, tsMpc, qY, rDu, y0, wOf(t) -> [wG, wL] }
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
  let ovrCount = 0, maxLevel = 0, minLevel = 1, carry = 0, iae = 0, started = false;
  for (let k = 0; k <= n; k++) {
    const t = k * o.dt, Ps = psep(y), a = sm.alg(y, cmd.z, ...wOf(t), Ps), level = y[4] / o.sepV, pm = a.Pp / 1e5 + (o.bias || 0) + o.noise * rand.normal(), ptm = a.Prt / 1e5 + o.noise * rand.normal();
    const sp = (o.sp ?? eq.Pp) / 1e5 + (t >= o.tStep ? o.dSp : 0), auto = t < o.tOff && o.mode !== 'open';
    out.t.push(t); out.pIn.push(a.Pp / 1e5); out.pTop.push(a.Prt / 1e5); out.pBase.push(a.Prb / 1e5); out.z.push(cmd.z); out.wL.push(a.wLout); out.wG.push(a.wGout); out.level.push(level); out.pSep.push(Ps / 1e5); out.mLr.push(y[3]); out.qLout.push(y[6] * o.qDrain); out.sp.push(sp);
    maxLevel = Math.max(maxLevel, level); minLevel = Math.min(minLevel, level); carry += over(y) * o.dt; if (auto) iae += Math.abs(sp - a.Pp / 1e5) * o.dt;
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
  return { ...out, yEnd: y, spBar: (o.sp ?? eq.Pp) / 1e5, ampOpen: ampOl, ampClosed: ampCl, meanP: mean(out.pIn.slice(cl0, Math.max(iOff, cl0 + 1))), meanZ: mean(out.z.slice(cl0, Math.max(iOff, cl0 + 1))), suppressed: Number.isFinite(ampCl) && ampCl < Math.max(1, 6 * o.noise + 0.02 * (o.sp ?? eq.Pp) / 1e5), maxLevel, minLevel, overrideSteps: ovrCount, surge: Math.max(0, maxLevel - o.levelSp) * o.sepV, carryOver: carry, iae: iae / HOUR, mpcActive: mp ? mp.active / Math.max(mp.moves, 1) : 0 };
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
    wG: Math.max(pr.mG, 1e-3), wL: Math.max(pr.mO + pr.mW, 1e-3), Ps: v.sepP * 1e5, Cv: v.chokeCv, kH: v.kH, kL: v.kL, aLp: clamp(avg('holdup', iF), 0.08, 0.85), rhoGnom: pf.rhoG, rough: line.roughness, ...(base ? { kG: base.kG } : {}) };
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

  // ---------- G. slugging and control ----------
  progress(0.4, 'Slugging model and control');
  const zCrit = sm0.critical(sm0.zFloor, 1), zAuto = zCrit === null ? 1 : Math.min(1, 2 * zCrit), sensorFailed = v.sensorState === 'failed', zCmd = clamp(v.chokePct > 0 ? v.chokePct / 100 : zAuto, 0.03, 1);
  // a failed inlet-pressure transmitter forces the loop to manual at an opening that is stable without feedback
  const zTarget = clamp(sensorFailed && zCrit !== null ? Math.min(zCmd, 0.9 * zCrit) : zCmd, Math.max(v.zMinPct / 100, 0.03, sm0.zFloor), Math.max(v.zMaxPct / 100, 0.05, sm0.zFloor)), unstable = zCrit !== null && zTarget > zCrit;
  const lin = sm0.linearise(zTarget), eqT = lin.steady, polesOL = sm0.poles(zTarget), kStat = (sm0.steady(Math.min(1, zTarget * 1.02)).Pp - sm0.steady(zTarget * 0.98).Pp) / (Math.min(1, zTarget * 1.02) - zTarget * 0.98) / 1e5;
  // step test at a stable opening → first-order-plus-dead-time model → tuning rules
  const zId = Math.max(zCrit === null ? 0.5 * zTarget : 0.6 * zCrit, sm0.zFloor), eqId = sm0.steady(zId), decay = Math.max(-sm0.growth(zId), 1e-5), tId = clamp(6 / decay, 0.5 * HOUR, 16 * HOUR), dz = 0.1 * zId;
  const stepRun = integrateStiff(sm0.f(zId + dz), eqId.y, 0, tId, { rtol: 1e-4, atol: 1e-3, hInit: 5, hMax: tId / 150, maxSteps: 4000 });
  const stepY = stepRun.y.map((y) => sm0.alg(y, zId + dz).Pp / 1e5), fo = identifyFOPDT(stepRun.t, stepY, dz), foUse = { K: fo.K, tau: Math.max(fo.tau, 1), theta: Math.max(fo.theta, v.deadTime, dtCtl) };
  const rules = tuningRules(foUse, Math.max(foUse.theta, v.tauCFactor * foUse.theta)), gainRatio = kStat !== 0 ? fo.K / kStat : 1, pole = tuneByPoles(lin, { tis: v.ctlMode === 'P' ? [0] : [1800, 5400, 14400], theta: v.deadTime + dtCtl });
  const tuneRows = rules.rows.map((r) => ({ ...r, kcT: r.kc * gainRatio, stable: pole.stable(r.kc * gainRatio, r.ti) }));
  let sel;
  if (v.tuning === 'manual') sel = { rule: 'Manual', kc: v.kcMan / 100, ti: v.ctlMode === 'P' ? 0 : v.tiMan, td: v.ctlMode === 'PID' ? v.tdMan : 0 };
  else if (v.tuning === 'simc' || v.tuning === 'zn') { const r = tuneRows.find((x) => x.mode === v.ctlMode && x.rule.startsWith(v.tuning === 'simc' ? 'SIMC' : 'Ziegler')); sel = { rule: r.rule + ' (gain-scheduled)', kc: r.kcT, ti: r.ti, td: r.td }; }
  else sel = { rule: 'Closed-loop pole search', kc: pole.kc, ti: v.ctlMode === 'P' ? 0 : pole.ti, td: v.ctlMode === 'PID' ? Math.max(v.deadTime, dtCtl) / 3 : 0 };
  let linStable = pole.stable(sel.kc, sel.ti);
  if (!linStable && v.tuning !== 'manual') { warnings.push({ level: 'warn', msg: `${sel.rule} settings do not stabilise the linearised model at ${(100 * zTarget).toFixed(0)} % opening; the pole-search settings are applied instead.` }); sel = { rule: 'Closed-loop pole search (fallback)', kc: pole.kc, ti: v.ctlMode === 'P' ? 0 : pole.ti, td: 0 }; linStable = pole.stable(sel.kc, sel.ti); }
  const tCtl = clamp(v.tCtl, 1, 72) * HOUR, ctlCfg = { dt: dtCtl, zTarget, tEnd: tCtl, tOff: 0.55 * tCtl, tStep: 0.2 * tCtl, dSp: v.spStep, kc: sel.kc, ti: sel.ti, td: sel.td, rate: v.rateLimit / 100, dead: v.deadTime, noise: v.noiseBar, cascade: !!v.cascade, override: !!v.override, feedForward: !!v.feedForward, antiWindup: true, sepV: v.slugCatcherVol, levelSp: clamp(v.levelSp, 10, 90) / 100, levelHi: clamp(v.levelHi, 20, 99) / 100, bias: v.sensorState === 'bias' ? v.sensorBias : 0, ...(v.pSet > 0 ? { sp: v.pSet * 1e5 } : {}), qDrain: v.qDrainM3h / 3600, zMin: v.zMinPct / 100, zMax: v.zMaxPct / 100, np: Math.round(clamp(v.mpcNp, 3, 60)), nc: Math.round(clamp(v.mpcNc, 1, 12)), rDu: v.mpcR, tsMpc: Math.max(dtCtl, 4 * dtCtl) };
  const pidRun = slugControl(sm0, { ...ctlCfg, mode: sensorFailed ? 'open' : 'pid' });
  await tick();
  const mpcRun = slugControl(sm0, { ...ctlCfg, mode: 'mpc', tEnd: ctlCfg.tOff, tOff: ctlCfg.tOff + 1 });
  // extended Kalman filter on the nonlinear model: riser liquid mass from the two noisy pressures
  const stride = Math.max(1, Math.ceil(pidRun.t.length / 1500)), hK = dtCtl * stride, rnd = rng(23), idx = pidRun.t.map((_, i) => i).filter((i) => i % stride === 0 && i > 0), g4 = 1 + Math.SQRT1_2;
  const sc4 = eqT.y.map((x) => Math.abs(x) || 1), sig4 = [0.01 * sc4[0], 0.05 * sm0.hc * sm0.dmdh, 0.2 * sc4[2], 0.25 * sc4[3]], // prior uncertainties sized to the sensitivity of each state
    F4 = (x, u) => { const s = ros2Step((t, y) => sm0.alg(y, u[0], undefined, undefined, u[1]).d, 0, x, hK); const J = s.M.map((r, i) => r.map((m, j) => ((i === j ? 1 : 0) - m) / (g4 * hK))), Phi = inv(madd(eye(4), J.map((r) => r.map((x2) => x2 * hK)), -1)); return { x: s.y.map((y, i) => clamp(y, 1e-6, 50 * sc4[i])), J: Phi }; };
  const ekfRes = ekf({ F: F4, h: (x, u) => { const a = sm0.alg(x, u[0], undefined, undefined, u[1]); return [a.Pp / 1e5, a.Prt / 1e5, a.Prb / 1e5]; }, x0: eqT.y.map((x, i) => x * (i === 3 ? 1.25 : 1)), P0: eye(4).map((r, i) => r.map((x) => x * sig4[i] ** 2)), Q: eye(4).map((r, i) => r.map((x) => x * (0.03 * sig4[i]) ** 2)), R: [[v.noiseBar ** 2 + 1e-3, 0, 0], [0, v.noiseBar ** 2 + 1e-3, 0], [0, 0, v.noiseBar ** 2 + 1e-3]],
    u: idx.map((i) => [pidRun.z[i], pidRun.pSep[i] * 1e5]), y: idx.map((i) => [pidRun.pIn[i] + v.noiseBar * rnd.normal(), pidRun.pTop[i] + v.noiseBar * rnd.normal(), pidRun.pBase[i] + v.noiseBar * rnd.normal()]) });
  const mTrue = idx.map((i) => pidRun.mLr[i]), mEst = ekfRes.x.map((x) => x[3]), half = Math.floor(idx.length / 4), ekfRmse = Math.sqrt(mean(mEst.slice(half).map((x, i) => (x - mTrue[half + i]) ** 2))), ekfRel = ekfRmse / Math.max(mean(mTrue), 1e-9);
  const tf = ss2tf(lin.A, lin.B, lin.C[0].map((x) => x / 1e5)), marg = loopMargins(foUse, rules.rows[1]), clPoles = pole.poles2(sel.kc, sel.ti);
  const scS = lin.ys.map((x) => Math.abs(x) || 1), dS = c2d(lin.A.map((r, i) => r.map((x, j) => (x * scS[j]) / scS[i])), lin.B.map((r, i) => [r[0] / scS[i]]), ctlCfg.tsMpc), cS = lin.C[0].map((x, j) => (x * scS[j]) / 1e5);
  let lqGain = null; try { lqGain = lqFinite(dS.Ad, dS.Bd, cS.map((a) => cS.map((b) => a * b)), [[v.mpcR]], 40).K0[0]; } catch { lqGain = null; }
  const slugSuppressed = !unstable || (pidRun.suppressed && !sensorFailed), chokeOut = slugSuppressed ? 100 * zTarget : 100 * (zCrit === null ? zTarget : 0.9 * zCrit), slugAmpOpen = unstable ? pidRun.ampOpen : 0;
  await tick();

  // ---------- B. restart and ramp-up ----------
  progress(0.55, 'Restart and ramp-up');
  const gelLen = sum(tAtShut.map((T) => (T < v.pourPoint ? st.ds : 0))), dpGel = (4 * v.yieldStress * gelLen) / id / 1e5, restartPressure = (v.sepP + so.headUphill + dpGel) * 1.05;
  const capSteel = STEEL.rho * STEEL.cp * Math.PI * ((id / 2 + wt) ** 2 - (id / 2) ** 2) * v.cMult, capCoat = sum(gridOf(hSea).C) - capSteel, cWall = Math.max(capSteel + 0.5 * Math.max(capCoat, 0), 1e3), gIn = hFlow * Math.PI * id, rTot = 1 / (wall.uFlow * Math.PI * id), gOut = 1 / Math.max(rTot - 1 / gIn, 1e-6);
  const cFl = st.holdup.map((H, i) => { const p = fm.at(st.P[i], st.T[i]); return A * (H * p.rhoL * p.cpL + (1 - H) * p.rhoG * p.cpG); }), q0 = clamp(v.qStartPct / 100, 0.05, 1), rampS = Math.max(v.rampHours, 0.05) * HOUR;
  const tWarmEnd = Math.max(3 * rampS, 30 * HOUR), nWarm = 600, refHyd = st.P.map((p) => fm.hydrateT(p) + v.hydMargin);
  const warm = warmUp({ s: st.s, ds: st.ds, tAmb: st.s.map((_, i) => tAmbOf(i)), T0: tAtShut, mdot: (t) => mdot * (q0 + (1 - q0) * Math.min(t / rampS, 1)), cp: cpMix, cFluid: cFl, cWall: st.s.map(() => cWall), gIn, gOut, tIn: line.tIn, dt: tWarmEnd / nWarm, nSteps: nWarm, ref: refHyd });
  const tW = warm.t.map((t) => t / HOUR), tOutEnd = warm.Tout[nWarm], iSteady = warm.Tout.findIndex((T, k) => warm.Tout.slice(k).every((x) => Math.abs(x - tOutEnd) <= 1)), restartTime = tW[Math.max(iSteady, 0)];
  const iSafe = warm.Tmin.findIndex((_, k) => warm.Tmin.slice(k).every((x) => x >= 0)), tSafe = iSafe >= 0 ? tW[iSafe] : null;
  const invOf = (q) => at('inv', q), qLiqOf = (q) => (at('qLiq', q) / 3600) * (q / qClamp(q)), qDrain = v.qDrainM3h / 3600, surgeAllow = (v.slugCatcherVol * v.catcherUsable) / 100;
  const surgeFor = (hours) => rampSurge({ nodes: [[0, q0 * rate], [hours * HOUR, rate]], inv: invOf, qLiq: qLiqOf, qDrain, inv0: liqInvSteady, dt: 120 });
  const ramp = surgeFor(v.rampHours);
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
  const pSafeRaw = fm.hydrateP(tAmbMin - v.hydMargin), pSafe = pSafeRaw === null ? 700 : pSafeRaw, headRiser = (so.rhoL * G * so.levelRiser) / 1e5, pEndAuto = Math.max(v.pBack * 1.1, pSafe - headRiser), pEndB = v.pBlowEnd > 0 ? Math.max(v.pBlowEnd, v.pBack * 1.02) : pEndAuto;
  const zGrid = linspace(Math.log(1), Math.log(Math.max(pB, 2) * 1.05), 14), zVals = zGrid.map((lp) => fm.at(Math.exp(lp), TgB).zG), Zf = (P) => interp1(zGrid, zVals, Math.log(Math.max(P / 1e5, 1)));
  const mLiqLine = so.liquidVol * volScale * so.rhoL, steelMass = STEEL.rho * Math.PI * ((id / 2 + wt) ** 2 - (id / 2) ** 2) * L, hBar = clamp(so.liquidVol / (A * L), 0, 0.95), mHcLine = mLiqLine * (1 - gB.wcut) + so.gasMass * volScale;
  const flashRate = (P, T) => { const pb = P / 1e5; if (pb < 1.2) return 0; return (mHcLine * Math.max(0, fm.at(pb * 0.97, T - KEL).wG - fm.at(pb, T - KEL).wG)) / (0.03 * P); };
  const bdCfg = { V: Math.max(so.gasVol * volScale, 1e-3), P0: pB * 1e5, T0: TgB + KEL, pBack: v.pBack * 1e5, area: (Math.PI * (v.orificeMm / 1000) ** 2) / 4, cd: v.cdBlow, k: v.kGas, mw: gB.mwG / 1000, Z: Zf, mode: 'wall', wallC: steelMass * STEEL.cp + mLiqLine * gB.cpL, wallUA: v.hGasWall * Math.PI * id * L * (1 - hBar), extUA: wall.uFlow * Math.PI * id * L, tAmb: mean(st.s.map((_, i) => tAmbOf(i))) + KEL, Tw0: TgB + KEL, pEnd: pEndB * 1e5, pMark: Math.max(pSafe - headRiser, v.pBack * 1.02) * 1e5, n: ntBlow, jt: Math.max(gB.jtG, 0), flash: v.bdFlash ? flashRate : null,
    hem: v.bdMode === 'hem' ? { mLiq: mLiqLine, rhoL: so.rhoL, frac: clamp(v.bdLiquidFrac / 100, 0, 0.9) } : null };
  need(pB > v.pBack * 1.02, `The line pressure when the blowdown valve opens (${pB.toFixed(1)} bara) is not above the flare back-pressure (${v.pBack} bara).`);
  const bd = blowdown(bdCfg), bdReached = bd.tEnd !== null, blowdownTime = (bd.tEnd ?? bd.t[bd.t.length - 1]) / HOUR, bdMinT = bd.minT - KEL, bdMinTw = bd.minTw - KEL, bdMinTd = bd.minTdown - KEL, bdEndP = bd.pFinal / 1e5, seabedPAfter = bdEndP + headRiser;
  const safeByTop = pSafe - headRiser > v.pBack * 1.02, peakStd = (bd.peak / (gB.mwG / 1000)) * VM_STD * DAY / 1e6; // million Sm³/d
  await tick();

  // ---------- D. pigging ----------
  progress(0.7, 'Pigging');
  const pigDpFric = Math.max((4 * v.pigFric * v.pigContact * 1e5 * v.pigSealLen) / id, 100); // seal contact force: μ × contact pressure × seal area
  let rPig = scan.res[scan.q.indexOf(qPig)]; try { if (!rPig) rPig = solve(qPig); } catch (e) { throw new Error(`No steady flow solution at the pigging rate (${v.pigRatePct} %): ${e.message}`); }
  const pigOf = (r) => pigRun({ s: r.s, z: r.z, vm: r.vm, holdup: r.holdup, vsl: r.vsl, rhoG: r.rhoG, rhoM: r.rhoM, rhoL: mean(r.rhoL), D: id, fric: pigDpFric, mass: v.pigMass, bypass: v.pigBypass / 100, leak: v.pigLeak / 100, qDrain });
  const pig = pigOf(rPig), pigVsRate = scan.res.map((r, k) => { const p = pigOf(r); return { q: scan.q[k], v: p.vMean, transit: p.transit === null ? null : p.transit / HOUR, surge: p.surge }; }).filter((p) => p.transit !== null);
  const pigTransit = pig.transit === null ? null : pig.transit / HOUR, pigSurge = pig.surge * volScale, pigRuns = v.pigInterval > 0 ? 365 / v.pigInterval : 0;
  const waxLen = sum(rPig.T.slice(1).map((T, i) => (T - (wall.uFlow * (T - rPig.tAmb[i])) / hFlow < v.wat ? rPig.ds : 0))), waxVol = Math.PI * id * (v.waxThk / 1000) * waxLen;
  const pInPigMax = rPig.pIn + pig.dpExtra / 1e5, waxSeries = { t: [0], v: [0] };
  if (waxVol > 0 && v.pigInterval > 0) { let w = 0; for (let k = 1; k <= 8; k++) { w += waxVol; waxSeries.t.push(k * v.pigInterval); waxSeries.v.push(w); w *= 1 - clamp(v.pigEff, 0, 100) / 100; waxSeries.t.push(k * v.pigInterval); waxSeries.v.push(w); } }
  const waxMax = Math.max(...waxSeries.v);

  // ---------- E. chemical injection ----------
  progress(0.75, 'Chemical injection');
  const S = fm.aq.S, tMinShut = Math.min(...tAtShut, ...comps.map((c) => interp1(tH, c.series, Math.min(v.tShut, v.tHorizon)))), dTshut = fm.hydrateT0(pAtShut) - tMinShut + v.inhMargin, dTsteady = Math.max(...st.P.map((p, i) => fm.hydrateT0(p) - st.T[i])) + v.inhMargin;
  const dTgov = v.dosingBasis === 'steady' ? dTsteady : v.dosingBasis === 'max' ? Math.max(dTshut, dTsteady) : dTshut, mWater = fm.rates.mW * rate, pOut = fm.at(v.sepP, tArr0, rate);
  const doseArgs = { inh: inhId, S, mWater, lean: v.leanWt, P: v.sepP, T: tArr0, qGasStd: fm.rates.qGasStd * rate, mOil: pOut.mO, mwOil: pOut.mwO };
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
  const contUsed = v.injRate > 0 ? v.injRate : contRate, mInj = (contUsed / DAY) * (lean * inh.rho + (1 - lean) * 1000) * lean, doseAchieved = contUsed > 0 ? (100 * mInj) / (mWater + mInj / lean) : 0, doseShort = steadyNeeded ? Math.max(0, Math.max(doseSteady.wt, v.inhRequired) - doseAchieved) : 0, pumpUtil = contUsed / 24 / Math.max(v.pumpMax, 1e-9);
  const megLoop = inh.rho > 1000 ? { rich: doseSteady.rich, inventory: ((doseSteady.rich / 1050) * (liqInvSteady / Math.max(at('qLiq', rate) / 3600, 1e-6)) + contUsed * v.megStorageDays), duty: (mWater * (2.257e6 + 4186 * 90) + doseSteady.mAq * 2400 * 90) / 1000 } : null;
  const subGov = dTgov - v.inhMargin, wcStd = fm.rates.wc, ldhi = subGov <= 0 ? 'No hydrate driving force at the governing condition; no low-dosage inhibitor needed.' : `${subGov <= 10 ? 'A kinetic inhibitor (KHI) is a candidate: sub-cooling ' + subGov.toFixed(1) + ' °C ≤ 10 °C, provided the hold time stays below the tested induction time.' : 'Kinetic inhibitors are not suitable: sub-cooling ' + subGov.toFixed(1) + ' °C exceeds about 10 °C.'} ${wcStd <= 0.5 ? 'An anti-agglomerant (AA) is a candidate: water cut ' + (100 * wcStd).toFixed(0) + ' % ≤ 50 %, needs a liquid hydrocarbon phase and restart testing.' : 'Anti-agglomerants are doubtful: water cut ' + (100 * wcStd).toFixed(0) + ' % exceeds 50 %.'}`;

  // ---------- F. heating ----------
  const tHold = fm.hydrateT(so.pSettle) + v.hydMargin, heatLen = clamp(v.heatLengthPct, 0, 100) / 100, heatW = sum(st.s.map((_, i) => (st.s[i] <= heatLen * L ? Math.max(0, tHold - tAmbOf(i)) / gridOf(buried ? 1e4 : st.z[i] >= 0 ? hAir : hSea).Rtot * st.ds : 0))), heatingPower = heatW / 1000 / Math.max(v.heatEff / 100, 0.05);
  const hot = warmUp({ s: st.s, ds: st.ds, tAmb: st.s.map((_, i) => tAmbOf(i)), T0: st.s.map((_, i) => tAmbOf(i)), mdot: () => v.hotOilRate, cp: 2000, cFluid: st.s.map(() => A * 850 * 2000), cWall: st.s.map(() => cWall), gIn: 150 * Math.PI * id, gOut, tIn: v.hotOilT, dt: 240, nSteps: 450 }), iHot = hot.Tout.findIndex((T) => T >= tHold), hotOilTime = iHot >= 0 ? hot.t[iHot] / HOUR : null;
  await tick();

  // ---------- I. envelope ----------
  progress(0.8, 'Operating envelope and optimisation');
  const inhibitedSteady = v.dosingBasis !== 'shutdown' && doseWt > 0, waxByPig = !!v.waxByPigging && v.pigInterval > 0, cons = [
    { key: 'margin', name: 'hydrate margin', type: 'min', limit: inhibitedSteady ? -1e3 : v.hydMargin, unit: '°C' }, { key: 'tArr', name: 'arrival temperature above WAT', type: 'min', limit: waxByPig ? -1e3 : v.wat + v.watMargin, unit: '°C' },
    { key: 'eros', name: 'erosional velocity', type: 'max', limit: 1, unit: '–' }, { key: 'pReq', name: 'inlet pressure (incl. slug-stabilising choke Δp)', type: 'max', limit: v.pAvail, unit: 'bara' }, { key: 'qLiq', name: 'separator liquid capacity', type: 'max', limit: v.qDrainM3h, unit: 'm³/h' }];
  const turndown = num(v.turndown, 0); if (turndown > 0) { scan.qSelf = scan.q.slice(); cons.push({ key: 'qSelf', name: 'minimum stable rate (flow suite)', type: 'min', limit: turndown, unit: '–' }); }
  const env = operatingEnvelope(scan, cons, clamp(rate, qLo, qHi));
  if (waxByPig) env.text.push(`wax managed by pigging every ${v.pigInterval} d (arrival may fall below the WAT)`);
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
    { id: 'blowdown', phrase: 'depressurising', name: 'Depressurise through the blowdown valve', feasible: safeByTop && bdReached && blowdownTime + v.tDecision <= cooldownTime && Math.min(bdMinTw, bdMinTd) >= v.tMinDesign, cost: bd.discharged * gasValue + (blowdownTime / 24) * bblDay * rate * v.oilPrice * 0.25, lead: blowdownTime },
    { id: 'heat', phrase: 'switching on the heating', name: 'Electrical heating', feasible: heatingPower > 0 && heatingPower <= v.heatMaxKw, cost: heatingPower * v.tShut * v.elecPrice, lead: 0.5 }];
  if (neverCools) strategies[0].feasible = true;
  const stratBest = strategies.filter((s) => s.feasible).sort((a, b) => a.cost - b.cost)[0] || null, stratSel = strategies.find((s) => s.id === v.preserve) || strategies[0];
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
  const inhibitorRate = contUsed + (stratSel.id === 'inhibit' ? ((batchVol + restartInj) * shutdownsYr) / 365 : 0), inhCostDay = megLoop ? inhibitorRate * v.inhPrice * (v.megLossPct / 100) + megLoop.duty * 24 * v.elecPrice * (contUsed > 0 ? 1 : 0) : inhibitorRate * v.inhPrice;

  // ---------- H. alarms, interlocks, event sequence ----------
  progress(0.9, 'Operating logic');
  const pOperate = pIn0 + (eqT.Prt - sm0.p.Ps) / 1e5, steadyMargin = Math.min(...st.T.map((T, i) => T - fm.hydrateT(st.P[i]))), worstSurge = Math.max(ramp.vMax, pigSurge, v.slugSurge, slugSuppressed ? 0 : pidRun.surge), catcherLevel = clamp(50 + (100 * worstSurge) / v.slugCatcherVol, 0, 400);
  const almVals = { tArr: tArr0, pIn: pOperate, hydMargin: steadyMargin, catcherLevel, tBlowMin: Math.min(bdMinTw, bdMinTd), noTouch, erosion: at('eros', rate), restartP: restartPressure, pPig: pInPigMax, slugAmp: slugSuppressed ? pidRun.ampClosed * (unstable ? 1 : 0) : slugAmpOpen, cooldown: cooldownTime, doseShort, pumpUtil: 100 * pumpUtil, watMargin: tArr0 - v.wat, ...(v.plugTime > 0 ? { plugTime: v.plugTime } : {}) };
  const rulesTbl = (Array.isArray(v.alarms) ? v.alarms : []).filter((r) => r && r.key && r.tag), raised = evaluateAlarms(almVals, rulesTbl), alarms = raised.map((a) => ({ tag: a.tag, level: a.level, msg: a.msg }));
  for (const l of env.limits) if (l.margin < 0) alarms.push({ tag: 'ENV-' + l.key.toUpperCase(), level: 'alarm', msg: `Operating point violates the ${l.name} limit (${txt(l.value)} against ${txt(l.limit)} ${l.unit}).` });
  if (sensorFailed) alarms.push({ tag: 'PT-100', level: 'alarm', msg: `Inlet-pressure transmitter failed: slug controller forced to manual at ${(100 * zTarget).toFixed(0)} % opening.` });
  if (restartPressure > v.pAvail) alarms.push({ tag: 'RST-P', level: 'trip', msg: `Restart needs ${restartPressure.toFixed(0)} bara but only ${v.pAvail} bara is available.` });
  const trips = raised.filter((a) => a.level === 'trip'), safeState = trips.length ? stateMachine(['trip']).state : 'PRODUCING';
  // planned shutdown → preservation → restart sequence
  const tAct = Math.max(0, cooldownTime - stratSel.lead - v.tDecision), tRestart = Math.max(v.tShut, 0.1), seq = [{ t: 0, tag: 'XV-001', action: 'Close production choke and wing valves (planned shutdown)', ev: 'shutdown', set: { rate: 0 } }];
  if (inhibitedSteady === false && stratSel.id === 'inhibit' && tProtect > 0) seq.unshift({ t: 0, tag: 'P-201', action: `Inhibitor front already through the line: injection started ${tProtect.toFixed(1)} h before shut-in`, ev: null, set: {} });
  if (stratSel.id === 'inhibit' && tAct < tRestart) seq.push({ t: tAct, tag: 'P-201', action: `Start ${inh.name} bullheading at ${v.pumpMax} m³/h (${batchVol.toFixed(1)} m³)`, ev: 'inhibit', set: { inj: v.pumpMax } }, { t: Math.min(tAct + tBullhead, tRestart), tag: 'P-201', action: 'Stop bullheading: line inhibited', ev: null, set: { inj: 0 } });
  if (stratSel.id === 'blowdown' && tAct < tRestart) seq.push({ t: Math.min(tAct, tRestart), tag: 'BDV-301', action: `Open blowdown valve (${v.orificeMm} mm orifice) to ${pEndB.toFixed(1)} bara`, ev: 'blowdown', set: { bdv: 1 } }, { t: Math.min(tAct + blowdownTime, tRestart), tag: 'BDV-301', action: 'Close blowdown valve: line depressurised', ev: null, set: { bdv: 0 } });
  if (stratSel.id === 'heat') seq.push({ t: Math.min(tAct, tRestart), tag: 'DEH-401', action: `Switch on heating at ${heatingPower.toFixed(0)} kW`, ev: 'heat', set: { heat: 1 } });
  seq.push({ t: tRestart, tag: 'XV-001', action: `Open wing valves; restart at ${v.qStartPct} % with ${inh.name} injection (needs ${restartPressure.toFixed(0)} bara)`, ev: 'restart', set: { rate: q0, heat: 0, inj: 0 } });
  if (tSafe !== null && tSafe > 0) seq.push({ t: tRestart + tSafe, tag: 'P-201', action: 'Stop restart inhibitor: whole line outside the hydrate region', ev: null, set: {} });
  seq.push({ t: tRestart + (rampReq ?? v.rampHours), tag: 'FIC-101', action: 'Ramp complete: full rate', ev: 'rampDone', set: { rate: 1 } }, { t: tRestart + Math.max(restartTime, rampReq ?? v.rampHours), tag: 'TI-102', action: 'Arrival temperature at steady state', ev: null, set: {} });
  if (v.pigLaunch >= 0 && pigTransit !== null) { const tl = tRestart + Math.max(restartTime, rampReq ?? v.rampHours) + v.pigLaunch; seq.push({ t: tl, tag: 'PL-501', action: `Launch pig at ${v.pigRatePct} % rate`, ev: 'pig', set: { rate: Math.min(1, qPig / Math.max(rate, 1e-9)) } }, { t: tl + pigTransit, tag: 'PR-502', action: `Receive pig; ${(pig.received * volScale).toFixed(0)} m³ of liquid ahead of it`, ev: 'pigReceived', set: { rate: 1 } }); }
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
  plots.push({ type: 'line', title: 'Blowdown: pressure and flare rate', xlabel: 'Time since the valve opened (h)', ylabel: 'bara · kg/s', series: [{ name: 'Line pressure (bara)', x: thin(bd.t.map((t) => t / HOUR)), y: thin(bd.P.map((p) => p / 1e5)) }, { name: 'Flare rate (kg/s)', x: thin(bd.t.map((t) => t / HOUR)), y: thin(bd.mdot) }], hlines: [{ y: Math.max(pSafe - headRiser, 0), label: 'hydrate-safe at the top' }] });
  plots.push({ type: 'line', title: 'Blowdown: temperatures', xlabel: 'Time since the valve opened (h)', ylabel: 'Temperature (°C)', series: [{ name: 'Gas in the line', x: thin(bd.t.map((t) => t / HOUR)), y: thin(bd.T.map((T) => T - KEL)) }, { name: 'Wall and liquid', x: thin(bd.t.map((t) => t / HOUR)), y: thin(bd.Tw.map((T) => T - KEL)) }, { name: 'Downstream of the valve (Joule–Thomson)', x: thin(bd.t.map((t) => t / HOUR)), y: thin(bd.T.map((T, i) => T - KEL - bdCfg.jt * Math.max(bd.P[i] - bdCfg.pBack, 0))), dash: true }], hlines: [{ y: v.tMinDesign, label: 'minimum design temperature' }] });
  plots.push({ type: 'line', title: 'Pig position and liquid pushed ahead', xlabel: 'Time since launch (h)', ylabel: 'km · m³', series: [{ name: 'Pig position (km)', x: thin(pig.t.map((t) => t / HOUR)), y: thin(pig.x.map(km)) }, { name: 'Liquid collected ahead (m³)', x: thin(pig.t.map((t) => t / HOUR)), y: thin(pig.slug) }], hlines: [{ y: v.slugCatcherVol, label: 'slug-catcher volume' }], vlines: pig.tFront !== null ? [{ x: pig.tFront / HOUR, label: 'slug front arrives' }] : [] });
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

  // tables
  tables.push({ title: 'Shutdown – preservation – restart sequence', columns: ['Time (h)', 'Tag', 'Action', 'State after'], rows: order.map((e, i) => [rd(e.t, 2), e.tag, e.action, stateAt[i]]), note: `Scheduler fired ${sched.fired.length} events in ${sched.steps} steps; inhibitor injected in the sequence ${acc.inj.toFixed(1)} m³; production lost ${(acc.lost / HOUR).toFixed(1)} full-rate hours. State machine ended in ${fsm.state}${fsm.rejected ? ` (${fsm.rejected} event(s) not permitted in the current state)` : ''}.` });
  tables.push({ title: 'Alarms and interlocks', columns: ['Tag', 'Variable', 'Type', 'Limit', 'Value', 'Level', 'Status', 'Action'], rows: rulesTbl.map((r) => { const val = almVals[r.key], hit = raised.find((a) => a.tag === String(r.tag)); return [String(r.tag), String(r.key), String(r.type), txt(+r.limit), txt(val), String(r.level || 'alarm'), isNum(val) ? (hit ? 'RAISED' : 'normal') : 'not evaluated', String(r.action || '—')]; }), note: `Fail-safe state on a trip: ${safeState === 'FAILSAFE' ? 'FAILSAFE is active — ' : ''}${Object.entries(FAILSAFE_POSITIONS).map(([k, x]) => `${k}: ${x}`).join('; ')}.` });
  tables.push({ title: 'Controller tunings (inlet pressure → choke)', columns: ['Rule', 'Mode', 'Kc at test point (%/bar)', 'Kc at target (%/bar)', 'Ti (s)', 'Td (s)', 'Stabilises target'], rows: [...tuneRows.map((r) => [r.rule, r.mode, rd(100 * r.kc, 3), rd(100 * r.kcT, 3), rd(r.ti, 0), rd(r.td, 1), r.stable ? 'yes' : 'no']), ['Closed-loop pole search', pole.ti > 0 ? 'PI' : 'P', '—', rd(100 * pole.kc, 3), rd(pole.ti, 0), 0, pole.decay > 0 ? 'yes' : 'no'], ['APPLIED: ' + sel.rule, v.ctlMode, '—', rd(100 * sel.kc, 3), rd(sel.ti, 0), rd(sel.td, 1), linStable ? 'yes' : 'no']],
    note: `Step test at ${(100 * zId).toFixed(1)} % opening: K = ${fo.K.toFixed(1)} bar per unit opening, τ = ${(fo.tau / 60).toFixed(1)} min, θ = ${(foUse.theta).toFixed(0)} s (fit error ${fo.rmse.toFixed(3)} bar). Ultimate gain ${(100 * rules.ku).toFixed(2)} %/bar, period ${(rules.pu / 60).toFixed(1)} min. Gain scheduling by the static-gain ratio ${gainRatio.toFixed(2)} to the target opening. Loop margins of the SIMC PI settings on the identified model: gain margin ${marg.gm === null ? '∞' : marg.gm.toFixed(2)}, phase margin ${marg.pm === null ? '—' : marg.pm.toFixed(0) + '°'}.` });
  const fmtC = (e) => `${(e[0] * HOUR).toFixed(2)}${Math.abs(e[1]) > 1e-9 ? ' ± ' + Math.abs(e[1] * HOUR).toFixed(2) + 'j' : ''}`, uniq = (pl) => pl.filter((e) => e[1] >= -1e-12);
  tables.push({ title: 'Linear model at the operating opening', columns: ['Item', 'Value'], rows: [['Open-loop poles (1/h)', uniq(polesOL).map(fmtC).join(', ')], ['Closed-loop poles with the applied controller (1/h)', uniq(clPoles).map(fmtC).join(', ')], ['Transfer function numerator (bar per unit opening, descending powers of s)', tf.num.map((x) => x.toExponential(3)).join(', ')], ['Transfer function denominator', tf.den.map((x) => x.toExponential(3)).join(', ')], ['Static gain (bar per unit opening)', rd(kStat, 2)], ['Finite-horizon LQ gain on the scaled states', lqGain ? lqGain.map((x) => x.toExponential(2)).join(', ') : '—'], ['MPC moves with active constraints (%)', rd(100 * mpcRun.mpcActive, 1)], ['EKF error on riser liquid mass (% of mean)', rd(100 * ekfRel, 2)]], note: 'State vector: gas and liquid mass in the feed pipeline, gas and liquid mass in the riser; input: choke opening; outputs: inlet and topside pressure.' });
  tables.push({ title: 'Operating-envelope limits', columns: ['Constraint', 'Type', 'Limit', 'Value at the operating rate', 'Margin', 'Bounds the window', 'At rate (%)'], rows: env.limits.map((l) => [l.name, l.type === 'min' ? '≥' : '≤', txt(l.limit), txt(l.value), txt(l.margin), l.bound, l.q === null ? '—' : rd(100 * l.q, 0)]), note: env.text.join('; ') + '.' });
  tables.push({ title: 'Chemical injection summary', columns: ['Item', 'Value', 'Unit'], rows: [['Inhibitor', inh.name, ''], ['Governing condition', v.dosingBasis === 'steady' ? 'steady flow' : v.dosingBasis === 'max' ? 'worse of steady flow and shutdown' : 'shutdown cold spot', ''], ['Required depression (incl. margin)', rd(dTgov, 1), '°C'], ['Dose in the aqueous phase', rd(doseWt, 1), 'wt %'], ['Dose attainable with this chemical', doseGov.attainable ? 'yes' : 'no', ''], ['Continuous injection required at the operating rate', rd(contRate, 1), 'm³/d'], ['Injection setting used', rd(contUsed, 1), 'm³/d'], ['Concentration reached with that setting', rd(doseAchieved, 1), 'wt %'], ['Injection point', rd(sInj / 1000, 2), 'km from the inlet'], ['Pump utilisation', rd(100 * pumpUtil, 0), '%'], ['Loss to the gas phase', rd(doseSteady.lossGas * DAY / 1000, 3), 't/d'], ['Loss to the hydrocarbon liquid', rd(doseSteady.lossOil * DAY / 1000, 3), 't/d'], ['Batch for a shutdown (bullheading)', rd(batchVol, 1), 'm³'], ['Bullheading time at pump capacity', rd(tBullhead, 2), 'h'], ['Injection during restart until hydrate-safe', rd(restartInj, 1), 'm³'], ['Time for the front to protect the whole line', rd(tProtect, 2), 'h'], ['Axial dispersion coefficient', rd(disp, 2), 'm²/s'], ['Glycol loop: rich stream', megLoop ? rd(megLoop.rich, 2) : '—', 'kg/s'], ['Glycol loop: inventory', megLoop ? rd(megLoop.inventory, 0) : '—', 'm³'], ['Glycol loop: regeneration duty', megLoop ? rd(megLoop.duty, 0) : '—', 'kW']], note: ldhi });
  tables.push({ title: 'Cooldown of the line and special components', columns: ['Location', 'Distance (km)', 'Start temperature (°C)', 'Time to hydrate limit (h)', 'Time to WAT (h)', 'Temperature at restart (°C)'], rows: [...[...new Set([0, Math.floor(nSt / 4), iMid, iCold, iRbU, nSt - 1])].sort((a, b) => a - b).map((i) => [i === iCold ? 'Line (cold spot)' : 'Line', rd(km(st.x[i]), 2), rd(st.T[i], 1), tCool[i] === null ? `> ${v.tHorizon}` : rd(tCool[i], 1), tWat[i] === null ? `> ${v.tHorizon}` : rd(tWat[i], 1), rd(tAtShut[i], 1)]), ...comps.map((c) => [c.name, rd(km(c.x), 2), rd(c.T0, 1), c.t === null ? `> ${v.tHorizon}` : rd(c.t, 1), (() => { const t = timeBelow(tH, c.series, v.wat); return t === null ? `> ${v.tHorizon}` : rd(t, 1); })(), rd(interp1(tH, c.series, Math.min(v.tShut, v.tHorizon)), 1)])], note: `Settle-out pressure ${so.pSettle.toFixed(1)} bara; U-value ${cool[iCold].U.toFixed(2)} W/m²K at shut-in conditions; effective thermal mass at the cold spot ${(cool[iCold].cEff / 1000).toFixed(0)} kJ/m/K; lumped time constant ${(cool[iCold].tau / HOUR).toFixed(1)} h.` });
  tables.push({ title: 'Preservation strategies for this shutdown', columns: ['Strategy', 'Feasible', 'Lead time (h)', 'Cost per shutdown ($)', 'Selected', 'Lowest cost'], rows: strategies.map((s) => [s.name, s.feasible ? 'yes' : 'no', rd(s.lead, 2), rd(s.cost, 0), s.id === stratSel.id ? '●' : '', stratBest && s.id === stratBest.id ? '●' : '']) });
  tables.push({ title: 'Optimisation of the operating point', columns: ['Method', 'Rate (% of case)', 'Dose (wt %)', 'Daily margin (M$/d)', 'Feasible', 'Evaluations / iterations'], rows: [...methods.map((m) => (m.x ? [m.name, rd(100 * m.x[0], 1), rd(m.x[1], 1), rd(profit(m.x[0], m.x[1], bestC.heat, bestC.pigs) / 1e6, 4), gOf(m.x[0], m.x[1], bestC.heat, bestC.pigs).every((g) => g > -2e-2) ? 'yes' : 'no', m.evals] : [m.name, '—', '—', '—', m.note || '—', 0])), ...combos.map((c) => [`Enumeration: heating ${c.heat ? 'on' : 'off'}, pigging programme ${c.pigs ? 'on' : 'off'}`, rd(100 * c.x[0], 1), rd(c.x[1], 1), rd(c.profit / 1e6, 4), c.feasible ? 'yes' : 'no', c.evals]), ['Mixed-integer LP (branch and bound, linearised)', milp.status === 'optimal' ? rd(100 * (qLo + milp.x[0]), 1) : '—', milp.status === 'optimal' ? rd(milp.x[1], 1) : '—', milp.status === 'optimal' ? rd(profit(qLo + milp.x[0], milp.x[1], milp.x[2], milp.x[3]) / 1e6, 4) : '—', milp.status === 'optimal' ? `heating ${milp.x[2] ? 'on' : 'off'}, pigging ${milp.x[3] ? 'on' : 'off'}` : milp.status, milp.nodes]],
    note: `Objective: oil revenue minus inhibitor and heating cost per day, subject to hydrate margin, wax, erosion, inlet-pressure and separator limits evaluated on the response surfaces. Ramp optimisation (three segments): ${(rampOpt.x[0]).toFixed(1)} h ramp with ${(100 * rampOpt.x[1]).toFixed(0)} % and ${(100 * Math.max(rampOpt.x[2], rampOpt.x[1])).toFixed(0)} % of the rise after one and two thirds, objective ${rampOpt.f.toFixed(2)} against ${rampLin.toFixed(2)} for the linear ramp.` });
  tables.push({ title: 'Response-surface fit of the mechanistic model', columns: ['Quantity', 'RMS error', 'Largest error', 'Leave-one-out RMS error', 'Unit'], rows: [['Inlet pressure', rd(sur.pIn.rmse, 3), rd(sur.pIn.maxErr, 3), rd(sur.pIn.loo, 3), 'bar'], ['Arrival temperature', rd(sur.tArr.rmse, 3), rd(sur.tArr.maxErr, 3), rd(sur.tArr.loo, 3), '°C'], ['Hydrate margin', rd(sur.margin.rmse, 3), rd(sur.margin.maxErr, 3), rd(sur.margin.loo, 3), '°C'], ['Liquid inventory', rd(sur.inv.rmse, 2), rd(sur.inv.maxErr, 2), rd(sur.inv.loo, 2), 'm³'], ['Erosional ratio', rd(sur.eros.rmse, 4), rd(sur.eros.maxErr, 4), rd(sur.eros.loo, 4), '–']], note: `Cubic polynomials in rate fitted to ${scan.q.length} kernel solutions between ${(100 * qLo).toFixed(0)} and ${(100 * qHi).toFixed(0)} % of the case rate.` });
  if (replay) tables.push({ title: 'Historical replay: agreement with the operating log', columns: ['Quantity', 'Points', 'Bias', 'RMS error', 'R²', 'Hold-out RMS before correction', 'Hold-out RMS after correction'], rows: [['Inlet pressure (bar)', replay.mP ? replay.mP.n : 0, replay.mP ? txt(replay.mP.bias, 2) : '—', replay.mP ? txt(replay.mP.rmse, 2) : '—', replay.mP ? txt(replay.mP.r2, 3) : '—', replay.rP ? txt(replay.rP.before, 2) : '—', replay.rP ? txt(replay.rP.after, 2) : '—'], ['Arrival temperature (°C)', replay.mT ? replay.mT.n : 0, replay.mT ? txt(replay.mT.bias, 2) : '—', replay.mT ? txt(replay.mT.rmse, 2) : '—', replay.mT ? txt(replay.mT.r2, 3) : '—', replay.rT ? txt(replay.rT.before, 2) : '—', replay.rT ? txt(replay.rT.after, 2) : '—']], note: 'The model is driven by the logged rate and choke opening only. The residual correction is a ridge regression on rate and choke opening fitted on every second row and tested on the others.' });
  tables.push({ title: 'Blowdown and pigging summary', columns: ['Item', 'Value', 'Unit'], rows: [['Gas volume blown down', rd(bdCfg.V, 0), 'm³'], ['Start pressure / temperature', `${pB.toFixed(1)} / ${TgB.toFixed(1)}`, 'bara / °C'], ['Time to the target pressure', bdReached ? rd(blowdownTime, 2) : `> ${blowdownTime.toFixed(1)}`, 'h'], ['Time to the hydrate-safe pressure at the top', bd.tMark === null ? '—' : rd(bd.tMark / HOUR, 2), 'h'], ['Peak flare rate', `${bd.peak.toFixed(1)} kg/s (${peakStd.toFixed(2)} MSm³/d)`, ''], ['Gas discharged / liberated from the oil', `${(bd.discharged / 1000).toFixed(1)} / ${(bd.flashed / 1000).toFixed(1)}`, 't'], ['Liquid left in the line', rd((mLiqLine - bd.liquidOut) / so.rhoL, 0), 'm³'], ['Residual liquid head on the low points', rd(headRiser, 1), 'bar'], ['Seabed pressure after blowdown', rd(seabedPAfter, 1), 'bara'], ['Hydrate-safe pressure at ambient', pSafeRaw === null ? 'no hydrate at ambient' : rd(pSafe, 1), 'bara'], ['Equivalent valve Cv', rd((v.cdBlow * bdCfg.area) / 1.7e-5, 1), 'US gpm/psi^0.5'], ['Pig differential pressure (seal friction ' + (pigDpFric / 1e5).toFixed(2) + ' bar)', rd(pig.dpPig / 1e5, 2), 'bar'], ['Pig mean velocity', rd(pig.vMean, 2), 'm/s'], ['Liquid swept / leaked past the pig', `${(pig.swept * volScale).toFixed(0)} / ${(pig.leaked * volScale).toFixed(1)}`, 'm³'], ['Liquid arrival duration', rd(pig.duration / 60, 1), 'min'], ['Peak inlet pressure while pigging', rd(pInPigMax, 1), 'bara'], ['Wax removed per run', rd(waxVol * clamp(v.pigEff, 0, 100) / 100, 2), 'm³'], ['Largest wax inventory in the line', rd(waxMax, 2), 'm³'], ['Hot-oil circulation time to the hold temperature', hotOilTime === null ? '> 30' : rd(hotOilTime, 1), 'h']] });

  // KPIs
  const stt = (ok, warn) => (ok ? 'ok' : warn ? 'warn' : 'bad');
  kpis.push({ label: 'Cooldown time', value: rd(cooldownTime, 1), unit: neverCools ? 'h (not reached)' : 'h', status: stt(cooldownTime >= 12 || neverCools, cooldownTime >= 6), help: 'Time until the first point of the system reaches the hydrate temperature plus margin at the falling settle-out pressure.' });
  kpis.push({ label: 'No-touch time', value: rd(noTouch, 1), unit: 'h', status: stt(noTouch >= 4, noTouch >= 1.5), help: 'Cooldown time minus the decision allowance and the time the selected preservation needs.' });
  kpis.push({ label: 'Maximum shutdown', value: rd(maxShutdown, 1), unit: 'h', status: stt(maxShutdown >= v.tShut, maxShutdown >= 0.7 * v.tShut), help: preserved ? 'With the selected preservation in place, limited by the gelled-line restart pressure or the simulated horizon.' : 'Without preservation: the cooldown time.' });
  kpis.push({ label: 'Settle-out pressure', value: rd(so.pSettle, 1), unit: 'bara', status: 'ok' });
  kpis.push({ label: 'Restart pressure', value: rd(restartPressure, 1), unit: 'bara', status: stt(restartPressure <= 0.85 * v.pAvail, restartPressure <= v.pAvail), help: 'Separator pressure + settled liquid legs + gelled-crude yield term, with 5 % allowance.' });
  kpis.push({ label: 'Warm-up time', value: rd(restartTime, 1), unit: 'h', status: 'ok', help: 'Until the arrival temperature stays within 1 °C of its final value.' });
  kpis.push({ label: 'Ramp-up surge', value: rd(ramp.vMax, 1), unit: 'm³', status: stt(ramp.vMax <= surgeAllow, ramp.vMax <= v.slugCatcherVol), help: `Allowance ${surgeAllow.toFixed(0)} m³.` });
  kpis.push({ label: 'Blowdown time', value: rd(blowdownTime, 2), unit: bdReached ? 'h' : 'h (target not reached)', status: stt(bdReached && blowdownTime <= Math.max(cooldownTime - v.tDecision, 0.1), bdReached) });
  kpis.push({ label: 'Blowdown minimum temperature', value: rd(Math.min(bdMinT, bdMinTd), 1), unit: '°C', status: stt(Math.min(bdMinTw, bdMinTd) >= v.tMinDesign + 10, Math.min(bdMinTw, bdMinTd) >= v.tMinDesign), help: `Gas ${bdMinT.toFixed(1)} °C, wall ${bdMinTw.toFixed(1)} °C, downstream of the valve ${bdMinTd.toFixed(1)} °C; design minimum ${v.tMinDesign} °C.` });
  kpis.push({ label: 'Pig transit', value: pigTransit === null ? 'stalled' : rd(pigTransit, 2), unit: 'h', status: stt(pigTransit !== null && pig.vMean >= 0.5 && pig.vMean <= 5, pigTransit !== null) });
  kpis.push({ label: 'Pig liquid surge', value: rd(pigSurge, 1), unit: 'm³', status: stt(pigSurge <= surgeAllow, pigSurge <= v.slugCatcherVol), help: `Liquid arriving faster than the drain capacity; slug catcher ${v.slugCatcherVol} m³.` });
  kpis.push({ label: `${inh.name} dose`, value: rd(doseWt, 1), unit: 'wt %', status: stt(doseGov.attainable && doseWt < 50, doseGov.attainable) });
  kpis.push({ label: 'Inhibitor use', value: rd(inhibitorRate, 2), unit: 'm³/d', status: 'ok', help: 'Continuous injection plus shutdown batches averaged over the year.' });
  kpis.push({ label: 'Heating power', value: rd(heatingPower, 0), unit: 'kW', status: stt(heatingPower <= v.heatMaxKw, heatingPower <= 1.2 * v.heatMaxKw), help: 'To hold the heated length at the hydrate temperature plus margin during shutdown.' });
  kpis.push({ label: 'Critical choke opening', value: zCrit === null ? 'stable' : rd(100 * zCrit, 1), unit: '%', status: stt(!unstable, slugSuppressed), help: 'Opening above which the riser limit cycle (severe slugging) starts in open loop.' });
  kpis.push({ label: 'Slugging amplitude under control', value: rd(pidRun.ampClosed, 2), unit: 'bar', status: stt(slugSuppressed, false), help: `Open loop at the same opening: ${pidRun.ampOpen.toFixed(1)} bar peak to peak.` });
  kpis.push({ label: 'Operating window', value: env.feasible ? `${(100 * env.qMin).toFixed(0)}–${(100 * env.qMax).toFixed(0)}` : 'none', unit: '% of case rate', status: stt(env.feasible && rate >= env.qMin && rate <= env.qMax, env.feasible) });
  kpis.push({ label: 'Uptime', value: rd(100 * uptime, 2), unit: '%', status: stt(uptime >= 0.95, uptime >= 0.9) });
  kpis.push({ label: 'Deferred production', value: rd(deferredVolume, 0), unit: 'Sm³/y', status: 'ok' });

  // warnings and recommendations
  if (pic.source === 'kernel estimate') warnings.push({ level: 'info', msg: 'The flow suite has not published a profile for this case; the kernel steady solution is used as the starting point.' });
  if (!neverCools && v.tShut > cooldownTime && stratSel.id === 'none') warnings.push({ level: 'bad', msg: `The planned shutdown (${v.tShut} h) is longer than the cooldown time (${cooldownTime.toFixed(1)} h) and no preservation is selected: the line enters the hydrate region ${km(coldSpotX).toFixed(1)} km from the inlet.` });
  if (stratSel.id !== 'none' && !stratSel.feasible) warnings.push({ level: 'bad', msg: `The selected preservation (${stratSel.name}) is not feasible for this shutdown${stratBest ? `; ${stratBest.name} is` : ''}.` });
  if (!safeByTop && pSafeRaw !== null) warnings.push({ level: 'warn', msg: `Topside blowdown alone cannot make the seabed hydrate-safe: ${headRiser.toFixed(0)} bar of settled liquid head remains against a hydrate-safe pressure of ${pSafe.toFixed(1)} bara at ${tAmbMin.toFixed(0)} °C.` });
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
  if (pidRun.carryOver > 0.5 && unstable) warnings.push({ level: 'warn', msg: `With the controller in manual the slugs overfill the separator (${pidRun.carryOver.toFixed(0)} m³ carried over in the simulated period).` });
  if (!env.feasible) warnings.push({ level: 'bad', msg: 'No rate in the scanned range satisfies all operating limits: ' + env.text.join('; ') + '.' });
  else if (rate < env.qMin || rate > env.qMax) warnings.push({ level: 'warn', msg: `The operating rate (${v.rateFrac} %) lies outside the window ${(100 * env.qMin).toFixed(0)}–${(100 * env.qMax).toFixed(0)} %.` });
  if (v.severeSlugging) warnings.push({ level: 'info', msg: 'The flow suite reports severe slugging for this case; the control section shows what the choke loop can do about it.' });
  for (const a of raised) warnings.push({ level: a.level === 'trip' ? 'bad' : 'warn', msg: `${a.tag}: ${a.msg}${a.action ? ' → ' + a.action : ''}` });

  if (neverCools) rec.push(`No point reaches the hydrate temperature within ${v.tHorizon} h of shut-in: no preservation is needed for shutdowns up to that duration.`);
  else rec.push(`Act within ${noTouch.toFixed(1)} h of shut-in: the cold spot (${coldComp ? coldComp.name + ', ' : ''}${km(coldSpotX).toFixed(1)} km) reaches the hydrate temperature + ${v.hydMargin} °C after ${cooldownTime.toFixed(1)} h, and ${stratSel.phrase} takes ${stratSel.lead.toFixed(1)} h plus ${v.tDecision} h to decide.`);
  if (stratBest && stratBest.id !== stratSel.id) rec.push(`For a ${v.tShut} h shutdown the lowest-cost feasible preservation is ${stratBest.phrase} (about $${fmt(stratBest.cost, 3)} per event against $${fmt(stratSel.cost, 3)} for the selected one${stratSel.feasible ? '' : ', which is not feasible'}).`);
  if (doseWt > 0) rec.push(`Dose ${inhId} to ${doseWt.toFixed(0)} wt % of the water phase (${dTgov.toFixed(0)} °C depression incl. ${v.inhMargin} °C margin): ${batchVol.toFixed(1)} m³ per shutdown${contRate > 0 ? ` and ${contRate.toFixed(1)} m³/d continuously` : ''}; start injection ${tProtect.toFixed(1)} h before a planned shut-in so the front covers the whole line.`);
  rec.push(rampReq === null ? `Restart at ${v.qStartPct} % with ${restartPressure.toFixed(0)} bara available at the inlet; the liquid surge cannot be kept inside ${surgeAllow.toFixed(0)} m³ by ramping alone.` : `Restart at ${v.qStartPct} % (needs ${restartPressure.toFixed(0)} bara) and ${rampReq > 0.99 * v.rampHours || rampReq > 0.5 ? `ramp at ≤ ${((100 - v.qStartPct) / Math.max(rampReq, 0.05)).toFixed(0)} %/h (${rampReq.toFixed(1)} h to full rate) to keep the liquid surge below ${surgeAllow.toFixed(0)} m³` : `ramp as planned over ${v.rampHours} h: the liquid surge (${ramp.vMax.toFixed(0)} m³) stays below ${surgeAllow.toFixed(0)} m³ even for a fast ramp`}${(tSafe ?? restartTime) > 0.05 ? `; keep ${inhId} on for ${(tSafe ?? restartTime).toFixed(1)} h until the whole line is outside the hydrate region` : ''}.`);
  rec.push(safeByTop ? `Depressurise through the ${v.orificeMm} mm orifice to ${pEndB.toFixed(1)} bara: ${blowdownTime.toFixed(1)} h, peak flare ${bd.peak.toFixed(1)} kg/s, coldest metal ${Math.min(bdMinTw, bdMinTd).toFixed(0)} °C.` : `Do not rely on topside blowdown for hydrate protection: the settled liquid leaves ${seabedPAfter.toFixed(0)} bara at the seabed (hydrate-safe below ${pSafe.toFixed(0)} bara); plan inhibitor displacement or heating instead.`);
  if (pigTransit !== null) rec.push(`Pig at ${v.pigRatePct} % rate: ${pigTransit.toFixed(1)} h transit at ${pig.vMean.toFixed(1)} m/s; expect ${(pig.received * volScale).toFixed(0)} m³ of liquid over ${(pig.duration / 60).toFixed(0)} min${pigSurge > surgeAllow ? ` — ${pigSurge.toFixed(0)} m³ more than the drain can take, so lower the pigging rate or pre-drain the slug catcher` : ', inside the slug-catcher allowance'}${v.pigInterval > 0 ? `; every ${v.pigInterval} d (${pigRuns.toFixed(0)} runs a year)` : ''}.`);
  if (sensorFailed) rec.push(`Repair the inlet-pressure transmitter: until then keep the choke at ${(100 * zTarget).toFixed(0)} % in manual (costs ${Math.max(0, eqT.Pp / 1e5 - sm0.steady(Math.min(1, zAuto)).Pp / 1e5).toFixed(0)} bar of back-pressure against controlled operation at ${(100 * zAuto).toFixed(0)} %).`);
  else if (zCrit !== null) rec.push(slugSuppressed ? `Run the choke at ${(100 * zTarget).toFixed(0)} % under inlet-pressure control (Kc ${(100 * sel.kc).toFixed(1)} %/bar${sel.ti > 0 ? `, Ti ${(sel.ti / 60).toFixed(0)} min` : ''}): the open-loop limit cycle starting at ${(100 * zCrit).toFixed(0)} % (${pidRun.ampOpen.toFixed(0)} bar swings) is suppressed to ${pidRun.ampClosed.toFixed(2)} bar and the inlet pressure is ${(sm0.steady(zCrit).Pp / 1e5 - eqT.Pp / 1e5).toFixed(0)} bar lower than at the largest stable manual opening.` : `Keep the choke at or below ${(90 * zCrit).toFixed(0)} % in manual; the tested controller does not stabilise ${(100 * zTarget).toFixed(0)} %.`);
  else rec.push('The riser is stable at every choke opening for this rate: no slug control is needed.');
  rec.push(env.feasible ? `Keep the rate between ${(100 * env.qMin).toFixed(0)} and ${(100 * env.qMax).toFixed(0)} % of the case rate (${env.text.slice(0, 2).join('; ')}); the best daily margin is at ${(100 * qOpt).toFixed(0)} % with ${wOpt.toFixed(0)} wt % inhibitor, heating ${bestC.heat ? 'on' : 'off'}.` : `No feasible rate window: ${env.text.join('; ')}.`);
  if (heatingPower > 0) rec.push(`Electrical heating would need ${heatingPower.toFixed(0)} kW (${(heatingPower * 24 * v.elecPrice).toFixed(0)} $/d) to hold ${tHold.toFixed(0)} °C during a shutdown${hotOilTime !== null ? `; hot-oil circulation at ${v.hotOilRate} kg/s warms a cold line in ${hotOilTime.toFixed(1)} h` : ''}.`);

  const balances = [
    { name: 'Cooldown energy (stored heat released = heat lost to ambient, J)', in: eBal.drop, out: eBal.lost },
    { name: 'Blowdown mass (initial + liberated = remaining + discharged, kg)', in: bd.mass.initial + bd.flashed, out: bd.mass.final + bd.discharged },
    { name: 'Pig liquid (swept = received + still ahead + leaked, m³)', in: pig.swept, out: pig.received + pig.inPipe + pig.leaked },
    { name: 'Restart energy (inflow = outflow + losses + storage, J)', in: warm.energy.in, out: warm.energy.out + warm.energy.lost + warm.energy.stored },
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
    { key: 'catcherUsable', label: 'Usable surge fraction of that volume', unit: '%', value: 60, min: 5, max: 100 },
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
  { group: 'Shutdown and cooldown', tab: 'inputs', fields: [
    { key: 'tShut', label: 'Planned shutdown duration', unit: 'h', value: 24, min: 0, max: 2000 },
    { key: 'tHorizon', label: 'Cooldown simulated', unit: 'h', value: 48, min: 1, max: 720 },
    { key: 'hydMargin', label: 'Hydrate safety margin', unit: '°C', value: 3, min: 0, max: 15 },
    { key: 'tDecision', label: 'Decision / mobilisation allowance', unit: 'h', value: 2, min: 0, max: 48 },
    sel('preserve', 'Preservation selected', 'inhibit', [['none', 'None (restart within the cooldown time)'], ['inhibit', 'Inhibitor bullheading'], ['blowdown', 'Depressurise'], ['heat', 'Electrical heating']]),
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
    sel('bdMode', 'Discharge model', 'gas', [['gas', 'Gas only (liquid stays behind)'], ['hem', 'Two-phase, homogeneous equilibrium (ω-method)']]),
    { key: 'bdLiquidFrac', label: 'Liquid volume fraction at the valve inlet', unit: '%', value: 5, min: 0, max: 90, showIf: (v) => v.bdMode === 'hem' },
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
    { key: 'inhPrice', label: 'Chemical price', unit: '$/m³', value: 450, min: 0, max: 1e5 },
    { key: 'pumpMax', label: 'Injection pump capacity', unit: 'm³/h', value: 20, min: 0.01, max: 5000 },
    { key: 'injRate', label: 'Continuous injection setting (0 = as required)', unit: 'm³/d', value: 0, min: 0, max: 1e5 },
    { key: 'injX', label: 'Injection point', unit: 'm from the inlet', value: 0, min: 0, max: 1e6 },
    { key: 'bullheadVol', label: 'Extra volume for jumpers, trees and manifold', unit: 'm³', value: 8, min: 0, max: 5000 },
    { key: 'dispMult', label: 'Axial dispersion multiplier on the Taylor value', unit: '–', value: 150, min: 1, max: 1e5, help: 'Slug mixing spreads the front far more than single-phase turbulent dispersion.' },
    { key: 'megStorageDays', label: 'Glycol storage', unit: 'd', value: 3, min: 0, max: 60 },
    { key: 'megLossPct', label: 'Glycol make-up (losses)', unit: '% of circulation', value: 1, min: 0, max: 100 },
    { key: 'heatEff', label: 'Heating system efficiency', unit: '%', value: 70, min: 5, max: 100 },
    { key: 'heatLengthPct', label: 'Heated length', unit: '% of line', value: 100, min: 0, max: 100 },
    { key: 'heatMaxKw', label: 'Heating power available', unit: 'kW', value: 3000, min: 0, max: 1e6 },
    { key: 'elecPrice', label: 'Electricity price', unit: '$/kWh', value: 0.12, min: 0, max: 5 },
    { key: 'hotOilRate', label: 'Hot-oil circulation rate', unit: 'kg/s', value: 30, min: 0.5, max: 2000 },
    { key: 'hotOilT', label: 'Hot-oil supply temperature', unit: '°C', value: 80, min: 20, max: 200 },
  ] },
  { group: 'Slugging and control', tab: 'setup', help: 'Low-order riser model (gas and liquid mass in the feed pipeline and in the riser) with a topside choke; controllers act on the choke from the inlet pressure.', fields: [
    { key: 'chokeCv', label: 'Topside choke Cv (fully open)', unit: 'US gpm/psi^0.5', value: 400, min: 5, max: 20000 },
    { key: 'chokePct', label: 'Operating choke opening (0 = twice the critical opening)', unit: '%', value: 0, min: 0, max: 100 },
    { key: 'slugControl', label: 'Active slug control available', type: 'bool', value: true, help: 'Used in the envelope: with control the choke may run at twice the critical opening, without it at 90 % of it.' },
    sel('ctlMode', 'Controller', 'PI', ['P', 'PI', 'PID']),
    sel('tuning', 'Tuning', 'auto', [['auto', 'Closed-loop pole search on the linearised model'], ['simc', 'SIMC from the step test (gain-scheduled)'], ['zn', 'Ziegler–Nichols from the step test (gain-scheduled)'], ['manual', 'Manual']]),
    { key: 'kcMan', label: 'Manual gain', unit: '% opening per bar', value: -15, min: -1000, max: 1000, showIf: (v) => v.tuning === 'manual' },
    { key: 'tiMan', label: 'Manual integral time', unit: 's', value: 1800, min: 1, max: 1e6, showIf: (v) => v.tuning === 'manual' },
    { key: 'tdMan', label: 'Manual derivative time', unit: 's', value: 0, min: 0, max: 1e5, showIf: (v) => v.tuning === 'manual' },
    { key: 'tauCFactor', label: 'SIMC closed-loop time constant / dead time', unit: '–', value: 1, min: 0.2, max: 20 },
    { key: 'deadTime', label: 'Actuator and measurement dead time', unit: 's', value: 30, min: 0, max: 1800 },
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
  { group: 'Alarms, interlocks and availability', tab: 'setup', fields: [
    { key: 'alarms', label: 'Alarm and trip thresholds', type: 'table', columns: [{ key: 'tag', label: 'Tag', type: 'text' }, { key: 'key', label: 'Variable', type: 'text' }, { key: 'type', label: 'low / high', type: 'text' }, { key: 'limit', label: 'Limit' }, { key: 'level', label: 'alarm / trip', type: 'text' }, { key: 'action', label: 'Action', type: 'text' }],
      help: 'Variables: tArr (arrival °C), pIn (inlet bara incl. choke), hydMargin (°C, steady), catcherLevel (%), tBlowMin (°C), noTouch (h), cooldown (h), erosion (ratio), restartP (bara), pPig (bara), slugAmp (bar), doseShort (wt %), pumpUtil (%), watMargin (°C), plugTime (h).',
      value: [{ tag: 'TAL-102', key: 'tArr', type: 'low', limit: 33, level: 'alarm', action: 'Raise rate or start wax inhibitor' }, { tag: 'PAH-100', key: 'pIn', type: 'high', limit: 150, level: 'alarm', action: 'Open choke / check for restriction' }, { tag: 'PAHH-100', key: 'pIn', type: 'high', limit: 185, level: 'trip', action: 'Shut in wells (ESD level 2)' },
        { tag: 'TDAL-110', key: 'hydMargin', type: 'low', limit: 3, level: 'alarm', action: 'Start hydrate inhibitor' }, { tag: 'LAH-200', key: 'catcherLevel', type: 'high', limit: 80, level: 'alarm', action: 'Slow the ramp / reduce pig speed' }, { tag: 'LAHH-200', key: 'catcherLevel', type: 'high', limit: 95, level: 'trip', action: 'Close inlet ESD valve' },
        { tag: 'TALL-301', key: 'tBlowMin', type: 'low', limit: -29, level: 'alarm', action: 'Throttle the blowdown' }, { tag: 'KAL-120', key: 'noTouch', type: 'low', limit: 4, level: 'alarm', action: 'Preserve immediately on shutdown' }, { tag: 'PAH-130', key: 'restartP', type: 'high', limit: 180, level: 'alarm', action: 'Displace / heat before restart' }, { tag: 'XA-140', key: 'slugAmp', type: 'high', limit: 5, level: 'alarm', action: 'Close choke to the stable opening' }, { tag: 'AAL-210', key: 'doseShort', type: 'high', limit: 1, level: 'alarm', action: 'Raise the inhibitor injection rate' }] },
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

/** Fast calibration model: cold-spot temperature after a shutdown, line pressure during a blowdown, pig arrival time (reference-case fluid constants). */
function calModel(v0) {
  const v = { ...DEFAULTS, ...v0 }; for (const k of NUMERIC) v[k] = num(v[k], DEFAULTS[k]);
  const id = v.idMm / 1000, wt = v.wtMm / 1000, A = (Math.PI * id * id) / 4, line = caseLine({}), L = line.length, tS = Math.max(num(v0.tShutH, 8), 0.01), t0 = num(v0.t0, 50), tBm = Math.max(num(v0.tBlowMin, 30), 0.01), p0 = Math.max(num(v0.p0, 80), v.pBack + 1), vGas = num(v0.vGas, 2.5);
  const wall = wallLayers(v, id, wt, 800, 800, []), c = cooldown({ ri: id / 2, layers: wall.layers, nr: 8, hIn: v.hInShut, hOut: 800, T0: t0, tAmb: v.tSeabed, cFluid: A * (0.45 * 800 * 2300 + 0.55 * 60 * 2600), dt: (tS * HOUR) / 48, nSteps: 48 });
  const b = blowdown({ V: (v.lineVolume > 0 ? v.lineVolume : A * L) * 0.55, P0: p0 * 1e5, T0: 290, pBack: v.pBack * 1e5, area: (Math.PI * (v.orificeMm / 1000) ** 2) / 4, cd: v.cdBlow, k: v.kGas, mw: 0.02, Z: 0.85, mode: 'isothermal', dt: (tBm * 60) / 60, n: 60, maxFactor: 1, pEnd: v.pBack * 1.001e5 });
  const leak = 0.7 * (v.pigBypass / 100) * Math.sqrt((2 * Math.max((4 * v.pigFric * v.pigContact * 1e5 * v.pigSealLen) / id, 100)) / 450);
  return { tCold: c.Tf[48], pBlow: b.pFinal / 1e5, pigArrival: L / Math.max(vGas - leak, 0.05) / HOUR };
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
  return out;
}

const suite = {
  id: 'ops', num: 5, title: 'Operations, Control & Flow-Assurance Management', short: 'Operations', icon: '🎛️',
  tagline: 'Shutdown, cooldown, restart, blowdown, pigging, chemicals, slug control, operating logic and the operating envelope in one study.',
  description: 'Starting from the steady flow picture of the case, the suite solves the cooldown of the line by radial finite-volume conduction, the settle-out, restart pressure, warm-up and ramp-up surge, the blowdown of the gas inventory, a pig run with the liquid it pushes, inhibitor dosing with the advection–dispersion of the front, and heating. A four-state riser model gives the choke opening where slugging starts and is stabilised by PID and by linear model-predictive control with Kalman filtering. The rate window, an optimised operating point, a timed shutdown–restart sequence with alarms, and a replay of the operating log complete the study.',
  guide: ['Run the fluid, network and flow suites first (or use the reference case) and link their values on the Inputs tab.', 'Set the shutdown duration, the preservation strategy and the restart ramp; check the cooldown and no-touch times.', 'Size the blowdown orifice and check the minimum temperatures and the seabed pressure left by the liquid head.', 'Review pigging surge against the slug-catcher and the inhibitor dose, volumes and front travel time.', 'On Model setup choose the controller and tuning; compare PID with MPC and read the critical choke opening.', 'Read the operating window, the optimised operating point and the event sequence; attach a historian export as the operating log to replay it.'],
  implemented: ['transient mass balance', 'momentum balance', 'energy balance', 'component balances', 'equipment inventory', 'tank/separator level', 'valve actuator dynamics', 'pid control', 'pi control', 'p control', 'feedforward control', 'feedback control', 'cascade control', 'override/selective control', 'anti-windup', 'model predictive control', 'optimal control', 'state-space model', 'transfer-function model', 'state observer', 'kalman filter', 'extended kalman filter',
    'transient heat equation', 'fourier conduction', 'newton cooling', 'lumped-capacitance', 'multilayer cylindrical conduction', 'transient fluid energy equation', 'transient mass/energy balances', 'compressible-flow', 'critical-flow', 'homogeneous-equilibrium', 'joule–thomson cooling', 'pig force balance', 'differential-pressure equation', 'friction/contact-force', 'pig velocity equation', 'liquid inventory/displacement balance', 'bypass/leakage',
    'advection–diffusion equation', 'species conservation', 'mixing equations', 'partitioning models', 'inhibitor mass balance', 'linear programming', 'nonlinear programming', 'mixed-integer linear programming', 'mixed-integer nonlinear programming', 'dynamic optimization', 'interior-point', 'genetic algorithm', 'particle-swarm', 'bayesian optimization',
    'transient multiphase solver + mpc', 'hydrate-risk model + inhibitor optimizer', 'digital twin + state estimation', 'physics model + data-driven surrogate', 'mechanistic model + machine-learning residual correction', 'shutdown model + optimization', 'production optimization + flow-assurance constraints',
    'production rates', 'pressure and temperature', 'fluid inventories', 'current valve and choke openings', 'separator conditions', 'inhibitor injection rates', 'heating status', 'current alarm states', 'controller set points and operating mode', 'pressure', 'temperature', 'phase distribution', 'liquid accumulation',
    'production targets and permitted production ranges', 'startup and ramp-up profiles', 'shutdown sequences', 'restart schedules', 'minimum and maximum valve openings', 'choke limits', 'separator pressure and level constraints', 'depressurization and blowdown schedules', 'inhibitor injection limits', 'heating limits', 'pig-launch and pig-receive conditions', 'emergency shutdown actions', 'set points', 'allowable deviations', 'actuator limits', 'operating constraints', 'alarm thresholds and safety limits', 'minimum acceptable hydrate safety margin', 'maximum shutdown duration', 'acceptable liquid surge', 'heating', 'depressurization or production-rate adjustment',
    'current state from modules 1-4', 'production targets and schedules', 'startup/shutdown/restart/ramp profiles', 'valve/choke commands and limits', 'depressurization/blowdown schedules', 'inhibitor/chemical injection rates, locations and capacity limits', 'pig geometry/friction/bypass and launch/receive schedule', 'separator/control settings', 'controller set points/gains', 'sensor states', 'alarm/interlock thresholds', 'operational and safety constraints',
    'startup/shutdown/restart trajectories', 'cooldown/warm-up time', 'restart pressure/rate requirements', 'depressurization/blowdown time and minimum temperature', 'transient liquid/solid inventory', 'pig position/velocity and generated liquid surge', 'inhibitor concentration/distribution and required dosage', 'hydrate/solids safety margin', 'equipment/control response', 'feasible operating envelope', 'maximum allowable shutdown duration', 'safe restart window', 'constraint violations, alarms/interlocks and recommended/optimized operating actions',
    'valve cv', 'pig friction', 'cooldown parameters', 'thermal time constants',
    'steady-to-transient initialization', 'event scheduler verification', 'valve opening/closing logic', 'controller logic', 'interlock logic', 'alarm logic', 'constraint handling', 'pid benchmark tests', 'mass balance during switching events', 'energy balance during shutdown', 'depressurization conservation', 'restart conservation', 'pig-tracking conservation', 'event-time accuracy', 'time-step independence', 'state-machine tests', 'fault-handling tests', 'fail-safe tests', 'optimization convergence', 'operating-envelope boundary verification',
    'shutdown records', 'blowdown/depressurization tests', 'cooldown measurements', 'pig arrival times', 'scada/historian data', 'emergency shutdowns'],
  referenceOnly: ['rotating-equipment', 'ratio control', 'heating/cooling strategy', 'nonlinear mpc', 'adaptive control', 'robust control', 'unscented kalman', 'homogeneous-relaxation', 'sequential quadratic programming', 'pump and compressor', 'pump/compressor', 'chemical inventories', 'pig location if a pig is already', 'maximum solids accumulation', 'valve actuator response', 'choke characteristics', 'pump curves', 'compressor maps', 'controller gains', 'pid parameters', 'sensor dynamics', 'sensor bias', 'process dead time', 'actuator dead time', 'inhibitor mixing/dispersion parameters', 'chemical-injection efficiency', 'pig bypass', 'pig velocity parameters', 'restart friction/pressure parameters',
    'commissioning data', 'startup records', 'restart records', 'pigging records', 'liquid-surging measurements', 'chemical tracer measurements', 'meg/methanol concentration measurements', 'valve-response data', 'compressor transient data', 'pump transient data', 'separator-level histories', 'field alarm/event histories'],
  equationsNote: 'Screening-level operations models. Cooldown: radial conduction with a lumped fluid node at every axial station (no axial conduction, no natural-circulation redistribution of heat; liquid settles into the low points between crests). Settle-out and cooling pressure assume a fixed gas mass with no inter-phase mass transfer. Blowdown treats the line gas as one lumped volume with a lumped wall-and-liquid heat sink (no axial pressure gradient, so long lines blow down somewhat slower than predicted); the two-phase option is the ω-method homogeneous-equilibrium model with a fixed inlet liquid fraction. Pigging is quasi-steady on the steady profile. The slugging model is a four-state riser model tuned to the steady solution: it reproduces the onset and period of riser-induced slugging and its response to the choke, not hydrodynamic slug statistics; its verdict should be confirmed with the transient flow suite. Pumps, compressors and their maps, ratio control, active cooling, nonlinear/adaptive/robust MPC, the unscented filter, relaxation (non-equilibrium) discharge and SQP are not solved.',
  inputs: INPUTS,
  presets: [
    { name: 'Planned shutdown and restart (reference tie-back)', values: { tShut: 24, preserve: 'inhibit', rampHours: 6, qStartPct: 30 } },
    { name: 'Emergency shutdown with blowdown', values: { tShut: 96, tHorizon: 96, preserve: 'blowdown', tBlowStart: 2, tDecision: 1, orificeMm: 45, pBack: 1.5, unplannedPerYear: 10, unplannedHours: 36, ntCool: 192 } },
    { name: 'Pigging campaign for wax management', values: { wat: 46, waxByPigging: true, pigRatePct: 90, pigInterval: 7, waxThk: 4, pigBypass: 4, pigFric: 0.45, pigLeak: 6, pigEff: 85, preserve: 'none', tShut: 6 } },
    { name: 'Severe-slugging control at low rate', values: { rateFrac: 40, chokePct: 0, ctlMode: 'PI', tuning: 'auto', cascade: true, tCtl: 24, dtCtl: 60, severeSlugging: true, qLoPct: 20, qHiPct: 120 } },
    { name: 'Continuous MEG injection with regeneration loop', values: { inhibitor: 'MEG', leanWt: 90, dosingBasis: 'max', inhPrice: 1100, megStorageDays: 5, megLossPct: 0.5, pumpMax: 40, preserve: 'inhibit' } },
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
    offer('chokeCv', n.chokeCv, 'Network suite: choke Cv', isNum(n.chokeCv) && n.chokeCv >= 5);
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
    note: 'Fits the U-value multiplier to cold-spot temperatures recorded during shutdowns, the blowdown discharge coefficient (equivalent valve Cv) to the pressure decay of a blowdown test and the pig friction to pig arrival times. The model uses the line data on the Inputs tab with reference-case fluid constants. Each row may hold any of the three measurements.',
    params: [{ key: 'uMult', label: 'U-value multiplier', lo: 0.5, hi: 2 }, { key: 'cdBlow', label: 'Blowdown discharge coefficient', lo: 0.3, hi: 1 }, { key: 'pigFric', label: 'Seal friction coefficient', lo: 0.05, hi: 1.5 }],
    columns: [{ key: 'tShutH', label: 'Time since shut-in', unit: 'h' }, { key: 't0', label: 'Temperature at shut-in', unit: '°C' }, { key: 'tCold', label: 'Fluid temperature', unit: '°C' }, { key: 'tBlowMin', label: 'Time since blowdown valve opened', unit: 'min' }, { key: 'p0', label: 'Pressure at opening', unit: 'bara' }, { key: 'pBlow', label: 'Line pressure', unit: 'bara' }, { key: 'vGas', label: 'Mixture velocity behind the pig', unit: 'm/s' }, { key: 'pigArrival', label: 'Pig arrival time', unit: 'h' }],
    targets: [{ key: 'tCold', label: 'Fluid temperature', unit: '°C' }, { key: 'pBlow', label: 'Line pressure', unit: 'bara' }, { key: 'pigArrival', label: 'Pig arrival time', unit: 'h' }],
    model: calModel,
    sample: [{ tShutH: 2, t0: 52, tCold: 44.6, tBlowMin: 10, p0: 82, pBlow: 71.6, vGas: 1.6, pigArrival: 4.53 }, { tShutH: 4, t0: 52, tCold: 38.8, tBlowMin: 20, p0: 82, pBlow: 64.9, vGas: 1.9, pigArrival: 3.54 }, { tShutH: 6, t0: 52, tCold: 33.8, tBlowMin: 30, p0: 82, pBlow: 58.3, vGas: 2.2, pigArrival: 3.07 }, { tShutH: 8, t0: 52, tCold: 29.5, tBlowMin: 45, p0: 82, pBlow: 49.8, vGas: 2.5, pigArrival: 2.63 }, { tShutH: 10, t0: 48, tCold: 24.6, tBlowMin: 60, p0: 78, pBlow: 38.8, vGas: 2.8, pigArrival: 2.24 }, { tShutH: 12, t0: 48, tCold: 21.5, tBlowMin: 90, p0: 78, pBlow: 27.1, vGas: 3.1, pigArrival: 2.05 }, { tShutH: 16, t0: 48, tCold: 16.8, tBlowMin: 120, p0: 78, pBlow: 18.8, vGas: 3.4, pigArrival: 1.87 }, { tShutH: 20, t0: 45, tCold: 13.1, tBlowMin: 150, p0: 75, pBlow: 12.9, vGas: 2, pigArrival: 3.44 }, { tShutH: 24, t0: 45, tCold: 11.2, tBlowMin: 180, p0: 75, pBlow: 9, vGas: 2.6, pigArrival: 2.49 }, { tShutH: 30, t0: 45, tCold: 8.1, tBlowMin: 240, p0: 75, pBlow: 4.4, vGas: 3, pigArrival: 2.13 }],
    validationSample: [{ tShutH: 3, t0: 50, tCold: 39.5, tBlowMin: 15, p0: 80, pBlow: 66.3, vGas: 1.8, pigArrival: 3.81 }, { tShutH: 7, t0: 50, tCold: 30, tBlowMin: 40, p0: 80, pBlow: 49.3, vGas: 2.4, pigArrival: 2.72 }, { tShutH: 11, t0: 50, tCold: 23.6, tBlowMin: 75, p0: 80, pBlow: 32.8, vGas: 2.9, pigArrival: 2.22 }, { tShutH: 15, t0: 46, tCold: 17.6, tBlowMin: 110, p0: 76, pBlow: 20.7, vGas: 3.2, pigArrival: 1.97 }, { tShutH: 22, t0: 46, tCold: 12.8, tBlowMin: 160, p0: 76, pBlow: 11.6, vGas: 2.3, pigArrival: 2.83 }, { tShutH: 28, t0: 46, tCold: 9.3, tBlowMin: 210, p0: 76, pBlow: 6.3, vGas: 2.7, pigArrival: 2.37 }],
  },
  verify,
  live: { key: 'log', label: 'Operating log', help: 'Follow a SCADA / historian export (CSV with time in hours, rate in % of the case rate, choke opening in %, inlet pressure in bara and arrival temperature in °C); the replay and the residual correction update as rows arrive.' },
};
export default suite;
