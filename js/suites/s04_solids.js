// Suite 4 — Hydrate & multiphase solids flow assurance.
// Hydrate thermodynamic driving force, classical nucleation, intrinsic / transfer-limited growth, a sectional
// population balance (with a quadrature-method-of-moments cross-check), slurry rheology, wall deposition and
// plugging marched in time on the case line; wax deposition, mineral-scale saturation, asphaltene screening and
// sand transport; one combined deposit profile for the backward coupling to the network and flow suites.
// SI units inside; bara and °C at the interfaces.
import { brent, clamp, interp1, linspace, logspace, rng, histogram, mean, quantile, isNum, rk45 } from '../core/num.js';
import { R, INHIBITORS, VM_STD, makeFluid, eosPhase, fluidModel, waterContent, hydrateDepression, inhibitorFor, hydrateT0 } from '../core/thermo.js';
import { G, gradient, frictionFactor, hInside } from '../core/pipe.js';
import { caseLine, steadyCase, flowPicture } from '../core/caseflow.js';
import { BASE } from '../data/basecase.js';

const KEL = 273.15, KB = 1.380649e-23, MW_W = 0.018015, PI = Math.PI;
/** Hydrate solid properties used throughout (structure-II natural-gas hydrate, rounded). */
export const HYDRATE = Object.freeze({ rho: 920, latent: 4.4e5, k: 0.6, cp: 2100 }); // kg/m³, J/kg, W/m/K, J/kg/K
const need = (cond, msg) => { if (!cond) throw new Error(msg); };
const num = (x, d) => (isNum(+x) ? +x : d);
const fin = (x, d = 0) => (Number.isFinite(x) ? x : d);
const nn = (x) => (Number.isFinite(x) ? x : null);

// =====================================================================================================
// 1. Hydrate thermodynamics: van der Waals–Platteeuw statistical model for methane structure I
// =====================================================================================================
let c1Fluid = null;
/** Fugacity (bar) of pure methane from the kernel Peng–Robinson model. */
export function methaneFugacity(Pbar, TK) { c1Fluid ||= makeFluid({ comp: { C1: 100 } }); return Pbar * Math.exp(eosPhase(c1Fluid, c1Fluid.z, Pbar, TK, 'vapour').lnphi[0]); }
/**
 * Langmuir cage occupancies of methane in structure I (Parrish–Prausnitz constants, C = A/T·exp(B/T) in 1/atm).
 * Returns { Cs, Cl (1/atm), thetaS, thetaL, hydrationNumber = 46 / (2θs + 6θl) }.
 */
export function langmuirOccupancy(TK, fBar) {
  const f = fBar / 1.01325, Cs = (3.7237e-3 / TK) * Math.exp(2708.8 / TK), Cl = (1.8372e-2 / TK) * Math.exp(2737.9 / TK);
  const thetaS = (Cs * f) / (1 + Cs * f), thetaL = (Cl * f) / (1 + Cl * f);
  return { Cs, Cl, thetaS, thetaL, hydrationNumber: 46 / Math.max(2 * thetaS + 6 * thetaL, 1e-9) };
}
/**
 * Three-phase (liquid water – hydrate – vapour) equilibrium pressure of methane hydrate from the equality of the
 * chemical potential of water in the hydrate lattice and in liquid water (van der Waals–Platteeuw, sI, T ≥ 0 °C).
 * Returns { P (bara) | null, thetaS, thetaL, hydrationNumber }.
 */
export function vdwpMethane(Tc) {
  const T = Tc + KEL, T0 = KEL, dmu0 = 1264, dh0 = -4858, dcp = -38.12, b = 0.141, dv = 4.6e-6; // J/mol, J/mol, J/mol/K, J/mol/K², m³/mol (empty lattice − liquid water)
  let I = 0; const n = 60; // ∫ Δh/(R T²) dT from T0 to T
  for (let k = 0; k < n; k++) { const t = T0 + ((T - T0) * (k + 0.5)) / n, dh = dh0 + dcp * (t - T0) + 0.5 * b * (t - T0) ** 2; I += ((dh / (R * t * t)) * (T - T0)) / n; }
  const g = (P) => {
    const f = methaneFugacity(P, T), o = langmuirOccupancy(T, f), xg = (f * 1e5) / (4.0e9 * Math.exp(-1700 * (1 / T - 1 / 298.15))); // dissolved methane lowers the water activity slightly
    const dmuH = -R * T * ((2 / 46) * Math.log(1 - o.thetaS) + (6 / 46) * Math.log(1 - o.thetaL));
    return dmuH - R * T * (dmu0 / (R * T0) - I + (dv * P * 1e5) / (R * T) - Math.log(1 - xg));
  };
  if (!(T >= 272.5) || g(1) > 0 || g(900) < 0) return { P: null, thetaS: null, thetaL: null, hydrationNumber: null };
  const P = brent(g, 1, 900, 1e-8), o = langmuirOccupancy(T, methaneFugacity(P, T));
  return { P, thetaS: o.thetaS, thetaL: o.thetaL, hydrationNumber: o.hydrationNumber };
}

// =====================================================================================================
// 2. Nucleation (classical nucleation theory) and induction time
// =====================================================================================================
/**
 * Classical nucleation rate. o: { TK, dT (K subcooling), TeqK, sigma (J/m², hydrate–water), theta (deg contact angle),
 * A (pre-exponential, 1/m³/s), het (true = heterogeneous) }. The volumetric driving force is ρ·L·ΔT/Teq.
 * Returns { J (1/m³ water/s), dG (J, barrier), rc (m, critical radius), f (contact-angle factor), exponent }.
 */
export function nucleationRate({ TK, dT, TeqK, sigma = 0.02, theta = 40, A = 3e7, het = true }) {
  const c = Math.cos((clamp(theta, 0, 180) * PI) / 180), f = het ? ((2 + c) * (1 - c) ** 2) / 4 : 1;
  if (!(dT > 0)) return { J: 0, dG: Infinity, rc: Infinity, f, exponent: Infinity };
  const dgv = (HYDRATE.rho * HYDRATE.latent * dT) / TeqK, dG = ((16 * PI * sigma ** 3) / (3 * dgv * dgv)) * f, ex = dG / (KB * TK);
  return { J: A * Math.exp(-Math.min(ex, 700)), dG, rc: (2 * sigma) / dgv, f, exponent: ex };
}
/** Mean induction time (s) of a water sample of volume V (m³) at a subcooling dT: 1 / (J·V), the mean of the exponential (Poisson) waiting time. */
export function inductionTime(dT, { TK = 277.15, TeqK, sigma, theta, A, het, V = 1e-3 } = {}) {
  const J = nucleationRate({ TK, dT, TeqK: TeqK ?? TK + dT, sigma, theta, A, het }).J;
  return J * V > 1e-30 ? 1 / (J * V) : 1e30;
}

// =====================================================================================================
// 3. Growth and dissociation kinetics
// =====================================================================================================
const fug = (Pbar, z) => Pbar * 1e5 * Math.exp(clamp(z - 1, -1.2, 0.3)); // Pa, first-order virial fugacity coefficient ln φ = Z − 1
const fugAt = (Pq, P, z) => Pq * 1e5 * Math.exp(clamp(((z - 1) * Pq) / Math.max(P, 1e-6), -1.2, 0.3));
/**
 * Hydrate growth flux on a particle surface with resistances in series (Kim–Bishnoi / Englezos intrinsic step on the
 * fugacity difference, liquid-film mass transfer, diffusion through a hydrate shell, heat removal from the particle).
 * o: { TK, P, Peq (bara), zG, kRef (mol/m²/Pa/s at 277.15 K), EaR (K), H (Pa·m³/mol, Henry constant of the gas in the
 *      continuous liquid), kFilm (m/s), kShell (m/s, Infinity = no shell), hPart (W/m²/K), dT (K), dHmol (J/mol gas) }
 * Returns { j (mol gas/m²/s, ≥ 0), jKin, jFilm, jShell, jHeat, kp (m/s, kinetic+film+shell conductance), dc (mol/m³), df (Pa), limiting }.
 */
export function hydrateGrowthRate({ TK, P, Peq, zG = 0.85, kRef = 1e-10, EaR = 13600, H = 2500, kFilm = Infinity, kShell = Infinity, hPart = Infinity, dT = 0, dHmol = 6e4 }) {
  const df = fug(P, zG) - fugAt(Peq, P, zG), z = { j: 0, jKin: 0, jFilm: 0, jShell: 0, jHeat: 0, kp: 0, dc: 0, df, limiting: 'none' };
  if (!(df > 0) || !(kRef > 0)) return z;
  const kStar = kRef * Math.exp(-EaR * (1 / TK - 1 / 277.15)), dc = df / H, kKin = kStar * H;
  const jKin = kKin * dc, jFilm = kFilm * dc, jShell = kShell * dc, jHeat = Number.isFinite(hPart) ? (hPart * Math.max(dT, 0)) / dHmol : Infinity;
  const kp = 1 / (1 / kKin + 1 / kFilm + 1 / kShell), j = 1 / (1 / (kp * dc) + 1 / jHeat);
  const m = Math.min(jKin, jFilm, jShell, jHeat);
  return { j, jKin, jFilm, jShell, jHeat, kp, dc, df, limiting: m === jKin ? 'intrinsic kinetics' : m === jFilm ? 'mass transfer' : m === jShell ? 'shell diffusion' : 'heat transfer' };
}
/**
 * Kim–Bishnoi dissociation flux (mol gas/m²/s) with an Arrhenius constant: K0·exp(−E/RT)·(f_eq − f).
 * Defaults are the methane values of Clarke & Bishnoi (K0 = 3.6e4 mol/m²/Pa/s, E = 81 kJ/mol).
 */
export function hydrateDissociationRate({ TK, P, Peq, zG = 0.85, K0 = 3.6e4, E = 81e3 }) {
  return K0 * Math.exp(-E / (R * TK)) * Math.max(fugAt(Peq, P, zG) - fug(P, zG), 0);
}
/** Shrinking-core shell conductance (m/s on the outer surface) for a converted fraction X of a droplet of radius Rd. */
export const shellConductance = (X, Rd, Dshell) => { const rc = Rd * Math.cbrt(clamp(1 - X, 0, 1)); return rc <= 0 ? 0 : Rd - rc < 1e-12 * Rd ? Infinity : (Dshell * rc) / (Rd * (Rd - rc)); };

// =====================================================================================================
// 4. Population balance: sectional (fixed pivot) and quadrature method of moments
// =====================================================================================================
/**
 * Geometric size grid with the fixed-pivot allocation tables for binary aggregation and binary breakage.
 * dCol: collision diameter floor (m) — particles smaller than this collide as if they had this size.
 * Returns { n, L[], v[], pk, pa (pair → lower pivot and number fraction assigned to it), bk, ba (breakage), gS, gD, gB (kernel geometry) }.
 */
export function pbeGrid(n = 14, Lmin = 2e-6, Lmax = 2e-2, dCol = 0) {
  n = Math.max(3, Math.round(n));
  const L = logspace(Lmin, Lmax, n), v = L.map((x) => (PI / 6) * x ** 3), find = (vs) => { let k = 0; while (k < n - 2 && v[k + 1] <= vs) k++; return k; };
  const pk = new Int32Array(n * n), pa = new Float64Array(n * n), pb = new Float64Array(n * n), gS = new Float64Array(n * n), gD = new Float64Array(n * n), gB = new Float64Array(n * n), bk = new Int32Array(n), ba = new Float64Array(n), bb = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      const vs = v[i] + v[j], q = i * n + j, a = Math.max(L[i], dCol), b = Math.max(L[j], dCol);
      if (vs >= v[n - 1]) { pk[q] = n - 2; pa[q] = 0; pb[q] = vs / v[n - 1]; } // overflow: volume kept in the top class
      else { const k = find(vs); pk[q] = k; pa[q] = (v[k + 1] - vs) / (v[k + 1] - v[k]); pb[q] = 1 - pa[q]; }
      gS[q] = (a + b) ** 3; gD[q] = (a + b) ** 2 * Math.abs(a * a - b * b); gB[q] = (a + b) ** 2 / (a * b);
    }
    const vd = v[i] / 2; // two equal daughters shared between the neighbouring pivots: conserves number (2) and volume
    if (vd < v[0]) { bk[i] = -1; } else { const k = find(vd); bk[i] = k; ba[i] = (2 * (v[k + 1] - vd)) / (v[k + 1] - v[k]); bb[i] = 2 - ba[i]; }
  }
  return { n, L, v, pk, pa, pb, bk, ba, bb, gS, gD, gB };
}
/**
 * Collision-frequency kernel (m³/s) between spheres of diameter Li and Lj and its parts.
 * e: { shear (1/s, laminar velocity gradient), eps (W/kg, turbulent dissipation), nu (m²/s), mu (Pa·s), TK, dRho (kg/m³), alpha (collision efficiency) }
 * Returns { total, shear (Smoluchowski orthokinetic), turbulent (Saffman–Turner), settling (differential Stokes), brownian }.
 */
export function aggregationKernel(Li, Lj, { shear = 0, eps = 0, nu = 1e-6, mu = 1e-3, TK = 277, dRho = 100, alpha = 1 } = {}) {
  const s = Li + Lj, sh = (shear / 6) * s ** 3, tu = 0.1618 * Math.sqrt(eps / nu) * s ** 3, se = (PI / 4) * s * s * ((Math.abs(dRho) * G) / (18 * mu)) * Math.abs(Li * Li - Lj * Lj), br = ((2 * KB * TK) / (3 * mu)) * (s * s) / (Li * Lj);
  return { total: alpha * (sh + tu + se + br), shear: sh, turbulent: tu, settling: se, brownian: br };
}
/**
 * Sectional population balance over a time t (s) in one well-mixed volume: aggregation (fixed pivot, explicit sub-steps),
 * binary breakage, growth/shrinkage in volume space (pivot shift with the exact volume rate and, for growth, exact number) and a source in one class.
 * N: number per class (any consistent basis). o: { beta: Float64Array n×n (m³/s on the same basis) | number (constant kernel),
 *   gv: [dv/dt per class, m³/s], S: [breakage frequency 1/s], src: { k, rate (1/s) }, frac (max fractional loss per sub-step), maxSub }.
 * Returns { N, sub, limited } — `limited` is true when the interval was too stiff for maxSub explicit sub-steps and
 * linearly implicit (modified Patankar–Euler) steps on the class volumes were used instead: unconditionally positive,
 * volume-conserving and convergent to the aggregation–breakage equilibrium.
 */
export function solvePBE(grid, N0, t, { beta = null, gv = null, S = null, src = null, frac = 0.25, maxSub = 400 } = {}) {
  const { n, v, pk, pa, pb, bk, ba, bb } = grid, ws = (grid.ws ||= { dN: new Float64Array(n), tmp: new Float64Array(n), dth: new Float64Array(n), lim: new Float64Array(n), V: new Float64Array(n), M: new Float64Array(n * n) });
  const N = new Float64Array(n), { dN, tmp, dth, lim, V, M } = ws; for (let i = 0; i < n; i++) N[i] = N0[i];
  const _c = 0, cst = typeof beta === 'number';
  const B = (i, j) => (cst ? beta : beta[i * n + j]);
  // death frequency of every class and the largest one among the classes that carry a noticeable share of the volume
  const rates = () => { let r = 0, vt = 0; for (let i = 0; i < n; i++) vt += N[i] * v[i]; for (let i = 0; i < n; i++) { let d = 0; if (beta !== null && N[i] > 0) for (let j = 0; j < n; j++) d += B(i, j) * N[j]; dth[i] = d; const tot = d + (S ? S[i] : 0); if (N[i] * v[i] > 1e-5 * vt && tot > r) r = tot; } return r; };
  let rate = rates(), sub = Math.max(1, Math.ceil((t * rate) / frac)), limited = false, dt = t / sub;
  if (sub > maxSub) { // too stiff for explicit sub-steps: linearly implicit, positivity-preserving and volume-conserving steps
    limited = true; sub = 0;
    const nImp = 3, h = t / nImp;
    const flow = (src, dst, fl) => { if (src === dst || !(fl > 0)) return; const c = fl / V[src]; M[dst * n + src] += c; M[src * n + src] -= c; };
    for (let st = 0; st < nImp; st++) {
      M.fill(0); for (let i = 0; i < n; i++) V[i] = N[i] * v[i];
      if (beta !== null) for (let i = 0; i < n; i++) {
        if (!(N[i] > 0)) continue;
        for (let j = i; j < n; j++) {
          if (!(N[j] > 0)) continue;
          const q = i * n + j, r = (i === j ? 0.5 : 1) * B(i, j) * N[i] * N[j], k = pk[q], vs = v[i] + v[j], wk = (pa[q] * v[k]) / vs, wk1 = (pb[q] * v[k + 1]) / vs;
          flow(i, k, r * v[i] * wk); flow(i, k + 1, r * v[i] * wk1); flow(j, k, r * v[j] * wk); flow(j, k + 1, r * v[j] * wk1);
        }
      }
      if (S) for (let i = 0; i < n; i++) { const k = bk[i]; if (k < 0 || !(S[i] > 0) || !(N[i] > 0)) continue; const r = S[i] * N[i]; flow(i, k, r * ba[i] * v[k]); flow(i, k + 1, r * bb[i] * v[k + 1]); }
      // (I − h·M) V' = V by elimination (column-diagonally-dominant M-matrix: no pivoting needed)
      for (let i = 0; i < n * n; i++) M[i] *= -h; for (let i = 0; i < n; i++) M[i * n + i] += 1;
      for (let c = 0; c < n; c++) { const pv = M[c * n + c]; for (let r2 = c + 1; r2 < n; r2++) { const m = M[r2 * n + c] / pv; if (m === 0) continue; for (let cc = c; cc < n; cc++) M[r2 * n + cc] -= m * M[c * n + cc]; V[r2] -= m * V[c]; } }
      for (let r2 = n - 1; r2 >= 0; r2--) { let x = V[r2]; for (let cc = r2 + 1; cc < n; cc++) x -= M[r2 * n + cc] * V[cc]; V[r2] = x / M[r2 * n + r2]; }
      for (let i = 0; i < n; i++) N[i] = Math.max(V[i], 0) / v[i];
    }
  }
  for (let s = 0; s < sub; s++) {
    if (s) rate = rates();
    if (rate > 0 && dt > 0) {
      dN.fill(0);
      // sparse classes with a very high death frequency are depleted exponentially, never below zero (pair-consistent, so volume is conserved)
      for (let i = 0; i < n; i++) { const x = dth[i] * dt; lim[i] = x > 1e-6 ? (1 - Math.exp(-x)) / x : 1; }
      if (beta !== null) for (let i = 0; i < n; i++) {
        if (!(N[i] > 0)) continue;
        for (let j = i; j < n; j++) {
          if (!(N[j] > 0)) continue;
          const q = i * n + j, r = (i === j ? 0.5 : 1) * B(i, j) * N[i] * N[j] * lim[i] * lim[j], k = pk[q];
          dN[i] -= r; dN[j] -= r; dN[k] += r * pa[q]; dN[k + 1] += r * pb[q];
        }
      }
      for (let i = 0; i < n; i++) N[i] = Math.max(N[i] + dt * dN[i], 0);
      if (S) { dN.fill(0); for (let i = 0; i < n; i++) { const k = bk[i]; if (k < 0 || !(S[i] > 0) || !(N[i] > 0)) continue; const r = N[i] * (1 - Math.exp(-S[i] * dt)); dN[i] -= r; dN[k] += r * ba[i]; dN[k + 1] += r * bb[i]; } for (let i = 0; i < n; i++) N[i] += dN[i]; }
    }
    if (src && src.rate > 0) N[src.k] += src.rate * dt;
  }
  if (gv) { // growth / shrinkage over the whole interval: every class is moved by its volume increment and shared between the neighbouring pivots
    tmp.fill(0);
    for (let i = 0; i < n; i++) {
      if (!(N[i] > 0)) continue;
      const dv = gv[i] * t, vn = v[i] + dv;
      if (dv === 0) { tmp[i] += N[i]; continue; }
      if (vn <= 0) continue; // dissolved completely
      if (vn <= v[0]) { tmp[0] += (N[i] * vn) / v[0]; continue; }
      if (vn >= v[n - 1]) { tmp[n - 1] += (N[i] * vn) / v[n - 1]; continue; }
      let k = dv > 0 ? i : 0; while (k < n - 2 && v[k + 1] <= vn) k++;
      const a = (v[k + 1] - vn) / (v[k + 1] - v[k]); tmp[k] += N[i] * a; tmp[k + 1] += N[i] * (1 - a);
    }
    N.set(tmp);
  }
  return { N, sub, limited };
}
/** Moments of a sectional distribution: { m0, m1, m2, m3, vol (Σ N v), d10, d32, d43 }. */
export function pbeMoments(grid, N) {
  let m0 = 0, m1 = 0, m2 = 0, m3 = 0, m4 = 0, vol = 0;
  for (let i = 0; i < grid.n; i++) { const L = grid.L[i], w = N[i]; m0 += w; m1 += w * L; m2 += w * L * L; m3 += w * L ** 3; m4 += w * L ** 4; vol += w * grid.v[i]; }
  return { m0, m1, m2, m3, vol, d10: m0 > 0 ? m1 / m0 : 0, d32: m2 > 0 ? m3 / m2 : 0, d43: m3 > 0 ? m4 / m3 : 0 };
}
function symEig(A) { // cyclic Jacobi for a small symmetric matrix: { val[], vec[][] (columns) }
  const n = A.length, a = A.map((r) => r.slice()), V = a.map((_, i) => a.map((__, j) => +(i === j)));
  for (let sweep = 0; sweep < 60; sweep++) {
    let off = 0; for (let p = 0; p < n; p++) for (let q = p + 1; q < n; q++) off += a[p][q] ** 2;
    if (off < 1e-26) break;
    for (let p = 0; p < n; p++) for (let q = p + 1; q < n; q++) {
      if (Math.abs(a[p][q]) < 1e-300) continue;
      const th = (a[q][q] - a[p][p]) / (2 * a[p][q]), t = Math.sign(th || 1) / (Math.abs(th) + Math.sqrt(th * th + 1)), c = 1 / Math.sqrt(t * t + 1), s = t * c;
      for (let k = 0; k < n; k++) { const x = a[k][p], y = a[k][q]; a[k][p] = c * x - s * y; a[k][q] = s * x + c * y; }
      for (let k = 0; k < n; k++) { const x = a[p][k], y = a[q][k]; a[p][k] = c * x - s * y; a[q][k] = s * x + c * y; }
      for (let k = 0; k < n; k++) { const x = V[k][p], y = V[k][q]; V[k][p] = c * x - s * y; V[k][q] = s * x + c * y; }
    }
  }
  return { val: a.map((r, i) => r[i]), vec: V };
}
/**
 * Quadrature nodes and weights from the first 2N moments of a size distribution (Wheeler algorithm).
 * Returns { L[], w[], ok } with Σ w L^k = m_k for k = 0 … 2N−1.
 */
export function qmomNodes(m) {
  const N = m.length >> 1, sc = m[1] / m[0], mm = m.map((x, k) => x / (m[0] * sc ** k)), a = new Array(N).fill(0), b = new Array(N).fill(0);
  let prev = new Array(2 * N).fill(0), cur = mm.slice();
  a[0] = mm[1];
  for (let k = 1; k < N; k++) {
    const nx = new Array(2 * N).fill(0);
    for (let l = k; l < 2 * N - k; l++) nx[l] = cur[l + 1] - a[k - 1] * cur[l] - b[k - 1] * prev[l];
    a[k] = nx[k + 1] / nx[k] - cur[k] / cur[k - 1]; b[k] = nx[k] / cur[k - 1];
    prev = cur; cur = nx;
  }
  if (b.slice(1).some((x) => !(x > 0)) || !a.every(Number.isFinite)) return { L: [m[1] / m[0]], w: [m[0]], ok: false };
  const J = a.map((_, i) => a.map((__, j) => (i === j ? a[i] : Math.abs(i - j) === 1 ? Math.sqrt(b[Math.max(i, j)]) : 0))), e = symEig(J);
  const L = e.val.map((x) => x * sc), w = e.val.map((_, j) => m[0] * e.vec[0][j] ** 2), ok = L.every((x) => x > 0 && Number.isFinite(x)) && w.every((x) => x >= 0);
  return { L, w, ok };
}
/**
 * Quadrature method of moments for the same processes as the sectional solver (length-based moments m0 … m(2N−1)).
 * o: { G(L) → dL/dt, beta(Li, Lj) → m³/s, S(L) → 1/s (binary equal-volume breakage), J (1/s source), L0 (m), steps }.
 * Returns { m[], nodes: { L[], w[] }, ok }.
 */
export function solveQMOM(m0, t, { G: Gf = null, beta = null, S = null, J = 0, L0 = 1e-6, steps = 200 } = {}) {
  const K = m0.length;
  let nodes = qmomNodes(m0), ok = nodes.ok;
  const rhs = (m) => {
    const q = qmomNodes(m); if (q.ok) nodes = q; else ok = false;
    const { L, w } = nodes, d = new Array(K).fill(0);
    for (let k = 0; k < K; k++) {
      let s = J * L0 ** k;
      for (let i = 0; i < L.length; i++) {
        if (Gf && k > 0) s += k * w[i] * Gf(L[i]) * L[i] ** (k - 1);
        if (S) s += w[i] * S(L[i]) * L[i] ** k * (2 ** (1 - k / 3) - 1);
        if (beta) for (let j = 0; j < L.length; j++) s += 0.5 * w[i] * w[j] * beta(L[i], L[j]) * ((L[i] ** 3 + L[j] ** 3) ** (k / 3) - L[i] ** k - L[j] ** k);
      }
      d[k] = s;
    }
    return d;
  };
  let m = m0.slice(); const dt = t / steps;
  for (let s = 0; s < steps; s++) { const k1 = rhs(m), mp = m.map((x, i) => x + dt * k1[i]), k2 = rhs(mp); m = m.map((x, i) => x + 0.5 * dt * (k1[i] + k2[i])); }
  const q = qmomNodes(m);
  return { m, nodes: q.ok ? { L: q.L, w: q.w } : { L: nodes.L, w: nodes.w }, ok: ok && q.ok };
}

// =====================================================================================================
// 5. Cohesion, slurry rheology, settling, particle momentum, wall deposition, porous plug
// =====================================================================================================
/**
 * Largest stable agglomerate from the Camargo–Palermo balance between the cohesive force and the shear stress.
 * { dp (m, primary particle), Fa (N, cohesive force between two primaries), mu0 (Pa·s), shear (1/s), phi (hydrate volume fraction), phiMax, fr (fractal dimension) }
 * Returns { ratio (dA/dp ≥ 1), dA (m), phiEff (effective volume fraction of the porous agglomerates) }.
 */
export function maxAgglomerateSize({ dp, Fa, mu0, shear, phi, phiMax = 4 / 7, fr = 2.5 }) {
  const p = clamp(phi, 0, 0.99 * phiMax), e = 3 - fr, cap = p > 1e-9 ? (phiMax / p) ** (1 / e) : 1e6;
  const g = (x) => x ** (4 - fr) - (Fa * (1 - (p / phiMax) * x ** e) ** 2) / (dp * dp * mu0 * Math.max(shear, 1e-9) * (1 - p * x ** e));
  let ratio = 1;
  if (g(1) < 0) { const hi = Math.min(cap * (1 - 1e-9), 1e6); ratio = g(hi) <= 0 ? hi : brent(g, 1, hi, 1e-9, 80); }
  return { ratio, dA: ratio * dp, phiEff: Math.min(p * ratio ** e, phiMax) };
}
/**
 * Relative viscosity of a suspension. model: 'mills' | 'krieger' (Krieger–Dougherty) | 'thomas' | 'einstein'.
 * phi is the effective volume fraction (agglomerates count with their trapped liquid). Capped at 1e4.
 */
export function slurryViscosity(phi, model = 'mills', { phiMax = 4 / 7, intrinsic = 2.5 } = {}) {
  const p = Math.max(phi, 0);
  if (model === 'einstein') return 1 + intrinsic * p;
  if (model === 'thomas') return Math.min(1 + 2.5 * p + 10.05 * p * p + 0.00273 * Math.exp(16.6 * p), 1e4);
  if (p >= phiMax * 0.9999) return 1e4;
  return Math.min(model === 'krieger' ? (1 - p / phiMax) ** (-intrinsic * phiMax) : (1 - p) / (1 - p / phiMax) ** 2, 1e4);
}
/**
 * Terminal settling velocity of a sphere. model: 'stokes' | 'schiller' (Schiller–Naumann drag, iterated).
 * Returns { v (m/s, positive = sinks, negative = rises), Re, Cd, n (Richardson–Zaki exponent), vHindered (at volume fraction phi) }.
 */
export function settlingVelocity(d, rhoP, rhoF, mu, { model = 'schiller', phi = 0 } = {}) {
  const dr = rhoP - rhoF, sgn = Math.sign(dr), vSt = (Math.abs(dr) * G * d * d) / (18 * mu);
  let v = vSt, Re = (rhoF * v * d) / mu, Cd = Re > 0 ? 24 / Re : Infinity;
  if (model !== 'stokes' && vSt > 0) for (let it = 0; it < 80; it++) {
    Re = (rhoF * v * d) / mu; Cd = Re < 1000 ? (24 / Re) * (1 + 0.15 * Re ** 0.687) : 0.44;
    const vn = Math.sqrt((4 * Math.abs(dr) * G * d) / (3 * Cd * rhoF));
    if (Math.abs(vn - v) < 1e-12 * v) { v = vn; break; } v = 0.5 * (v + vn);
  }
  Re = (rhoF * v * d) / mu;
  const n = Re < 0.2 ? 4.65 : Re < 1 ? 4.4 * Re ** -0.03 : Re < 500 ? 4.4 * Re ** -0.1 : 2.4;
  return { v: sgn * v, Re, Cd, n, vHindered: sgn * v * (1 - clamp(phi, 0, 0.99)) ** n };
}
/**
 * Particle momentum equation in a quiescent fluid: (ρp + ½ρf) dv/dt = (ρp − ρf) g − ¾ Cd ρf |v| v / d (drag by Schiller–Naumann,
 * buoyancy and added mass; history force neglected). Returns { t[], v[], vTerminal, tau (s, 63 % response time) }.
 */
export function particleRelaxation({ d, rhoP, rhoF, mu, tEnd = null }) {
  const vt = settlingVelocity(d, rhoP, rhoF, mu).v, tauS = ((rhoP + 0.5 * rhoF) * d * d) / (18 * mu), T = tEnd ?? 8 * tauS;
  const f = (t, y) => { const v = y[0], Re = Math.max((rhoF * Math.abs(v) * d) / mu, 1e-12), Cd = Re < 1000 ? (24 / Re) * (1 + 0.15 * Re ** 0.687) : 0.44; return [((rhoP - rhoF) * G - (0.75 * Cd * rhoF * Math.abs(v) * v) / d) / (rhoP + 0.5 * rhoF)]; };
  const r = rk45(f, [0], 0, T, { rtol: 1e-8, atol: 1e-14 }), v = r.y.map((y) => y[0]);
  let tau = T; for (let i = 1; i < v.length; i++) if (Math.abs(v[i]) >= 0.6321 * Math.abs(vt)) { const a = Math.abs(v[i - 1]), b = Math.abs(v[i]); tau = r.t[i - 1] + ((0.6321 * Math.abs(vt) - a) / (b - a || 1)) * (r.t[i] - r.t[i - 1]); break; }
  return { t: r.t, v, vTerminal: vt, tau };
}
/**
 * Turbulent deposition velocity of a particle onto a pipe wall (m/s): diffusion + eddy-impaction regimes capped at the
 * inertia-moderated plateau. V⁺ = min(0.14, 0.057 Sc^−2/3 + 4.5e-4 τ⁺²).
 */
export function depositionVelocity(d, rhoP, uStar, nu, rhoF, TK = 277) {
  if (!(uStar > 0)) return 0;
  const mu = nu * rhoF, tauP = (rhoP * d * d * uStar * uStar) / (18 * mu * nu), Sc = nu / (KB * TK / (3 * PI * mu * d));
  return uStar * Math.min(0.14, 0.057 * Sc ** (-2 / 3) + 4.5e-4 * tauP * tauP);
}
/** Kozeny–Carman permeability (m²) of a packed bed of porosity eps and grain size dp. */
export const kozenyCarman = (eps, dp) => (eps ** 3 * dp * dp) / (180 * (1 - eps) ** 2);
/**
 * Pressure gradient (Pa/m) for a superficial velocity v through a porous plug: Darcy term μv/k plus the Forchheimer
 * (Ergun inertial) term β ρ v². Returns { k (m²), beta (1/m), darcy, forchheimer, total }.
 */
export function porousGradient(v, mu, rho, eps, dp) {
  const k = kozenyCarman(eps, dp), beta = (1.75 * (1 - eps)) / (eps ** 3 * dp), darcy = (mu * v) / k, forchheimer = beta * rho * v * v;
  return { k, beta, darcy, forchheimer, total: darcy + forchheimer };
}
/** Superficial velocity (m/s) through a porous plug under a pressure gradient (Pa/m): root of the Darcy–Forchheimer quadratic. */
export function porousVelocity(dpdx, mu, rho, eps, dp) { const k = kozenyCarman(eps, dp), b = (1.75 * (1 - eps)) / (eps ** 3 * dp) * rho, a = mu / k; return b > 0 ? (-a + Math.sqrt(a * a + 4 * b * dpdx)) / (2 * b) : dpdx / a; }

// =====================================================================================================
// 6. Dissociation of a plug: Stefan moving boundary
// =====================================================================================================
/** Similarity constant λ of the one-phase Stefan problem: λ·exp(λ²)·erf(λ) = Ste/√π. */
export function stefanLambda(Ste) {
  const erf = (x) => { const t = 1 / (1 + 0.3275911 * x), y = 1 - ((((1.061405429 * t - 1.453152027) * t + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * Math.exp(-x * x); return y; };
  return brent((l) => l * Math.exp(l * l) * erf(l) - Ste / Math.sqrt(PI), 1e-9, 5, 1e-12);
}
/**
 * Radial melting of a hydrate plug that fills the bore after depressurisation (heat-transfer-controlled, quasi-steady
 * Stefan problem): heat flows from the surroundings through the overall coefficient U (on the bore radius) and the
 * annulus of released water to the dissociation front at Td.
 * { R (m, bore radius), U (W/m²/K), kW (W/m/K of the melt annulus), Tamb, Td (°C), eps (plug porosity), latent, rho, steps }
 * Returns { tAnalytic, tNumeric (s; Infinity when Tamb ≤ Td), r: [], t: [] (front radius against time) }.
 */
export function plugMeltTime({ R: Rb, U, kW = 0.58, Tamb, Td, eps = 0.4, latent = HYDRATE.latent, rho = HYDRATE.rho, steps = 400 }) {
  const dT = Tamb - Td, q = rho * (1 - eps) * latent;
  if (!(dT > 0)) return { tAnalytic: Infinity, tNumeric: Infinity, r: [Rb], t: [0] };
  const tAnalytic = (q * (Rb / (2 * U) + (Rb * Rb) / (4 * kW))) / dT, r = [Rb], t = [0];
  let tt = 0;
  for (let i = 0; i < steps; i++) { // front tracking: time for the front to move one radial increment at the local heat flux
    const s1 = Rb * (1 - i / steps), s2 = Rb * (1 - (i + 1) / steps), s = 0.5 * (s1 + s2), flux = dT / (1 / (U * Rb) + Math.log(Rb / s) / kW); // W/m per 2π
    tt += (q * s * (s1 - s2)) / flux; r.push(s2); t.push(tt);
  }
  return { tAnalytic, tNumeric: tt, r, t };
}

// =====================================================================================================
// 7. Wax
// =====================================================================================================
/**
 * Diffusivity of wax molecules in oil (m²/s). model: 'haydukMinhas' | 'wilkeChang'. mu in Pa·s, VA molar volume of the
 * wax (cm³/mol), MB solvent molar mass (g/mol).
 */
export function waxDiffusivity(TK, mu, { model = 'haydukMinhas', VA = 430, MB = 200, assoc = 1 } = {}) {
  const cP = mu * 1000;
  return model === 'wilkeChang' ? (7.4e-12 * Math.sqrt(assoc * MB) * TK) / (cP * VA ** 0.6) : (13.3e-12 * TK ** 1.47 * cP ** (10.2 / VA - 0.791)) / VA ** 0.71;
}
/**
 * Wax solubility curve: dissolved mass fraction falls exponentially below the wax appearance temperature.
 * Returns { dissolved, solid (mass fractions of the oil), dCdT (1/K) }.
 */
export function waxSolubility(T, wat, wTot, slope = 0.04) {
  if (T >= wat) return { dissolved: wTot, solid: 0, dCdT: 0 };
  const c = wTot * Math.exp(-slope * (wat - T));
  return { dissolved: c, solid: wTot - c, dCdT: slope * c };
}
/**
 * Local wax deposition rates on a cold wall.
 * o: { Tb, Tamb (°C), U (W/m²/K, clean overall coefficient on the bore), hIn (W/m²/K), kOil, rhoOil, muOil, wat, wTot (mass fraction),
 *      slope (1/K), delta (m, present thickness), Fw (wax fraction of the deposit), kDep (W/m/K), D (m), vL (m/s), gammaW (1/s wall shear rate),
 *      rhoMix, regime, mult (deposition multiplier), diffModel, wetFrac }
 * Returns { Ti (deposit surface °C), q (W/m²), dTdr (K/m), Dwo, jMol, jShear, jBrown (kg wax/m²/s), strip (shear-stripping factor 0–1),
 *           dDelta (m/s), dFw (1/s), Ueff }.
 */
export function waxDeposition({ Tb, Tamb, U, hIn, kOil = 0.13, rhoOil = 800, muOil = 3e-3, wat, wTot, slope = 0.04, delta = 0, Fw = 0.2, kDep = 0.25, D = 0.25, vL = 1, gammaW = 100, rhoMix = null, regime = '', mult = 1, diffModel = 'haydukMinhas', wetFrac = 1, rhoWax = 900 }) {
  const Ueff = 1 / (1 / U + delta / kDep), q = Ueff * (Tb - Tamb), Ti = Tb - q / Math.max(hIn, 1e-6), dTdr = q / kOil;
  const z = { Ti, q, dTdr, Dwo: 0, jMol: 0, jShear: 0, jBrown: 0, strip: 1, dDelta: 0, dFw: 0, Ueff };
  if (!(Ti < wat) || !(q > 0) || !(wTot > 0)) return z;
  const TK = Ti + KEL, Dwo = waxDiffusivity(TK, muOil, { model: diffModel }), sol = waxSolubility(Ti, wat, wTot, slope), bulk = waxSolubility(Tb, wat, wTot, slope);
  const jMol = rhoOil * Dwo * sol.dCdT * dTdr; // Fick's law with the solubility slope and the radial temperature gradient
  const dCr = 10e-6, sub = (5 * muOil) / (rhoOil * Math.max(Math.sqrt((gammaW * muOil) / rhoOil), 1e-6)); // crystal size, viscous sub-layer thickness
  const jShear = bulk.solid > 0 ? (rhoOil * 0.1 * (dCr / 2) ** 2 * gammaW * bulk.solid * bulk.solid) / sub : 0; // shear dispersion of precipitated crystals
  const jBrown = bulk.solid > 0 ? (rhoOil * ((KB * TK) / (3 * PI * muOil * dCr)) * bulk.solid) / sub : 0; // Brownian diffusion of crystals to the wall
  const nsr = ((/slug|bubble|churn/.test(regime) ? rhoMix ?? rhoOil : /annular/.test(regime) ? Math.sqrt((rhoMix ?? rhoOil) * rhoOil) : rhoOil) * vL * Math.max(delta, 1e-5)) / muOil;
  const strip = 1 / (1 + 0.055 * (nsr / 1000) ** 1.4); // shear stripping (Matzain form on a film Reynolds number in thousands)
  const j = mult * wetFrac * (jMol + jShear + jBrown), F = clamp(Fw, 0.02, 0.98), De = 1 / (1 + (64 * F * F) / (1 - F)), psi = De * (1 - F); // ageing: Cussler hindered diffusion into the gel
  const dDelta = ((j * (1 - psi)) / (rhoWax * F)) * strip, dFw = (j * psi) / (rhoWax * Math.max(delta, 2e-5));
  return { Ti, q, dTdr, Dwo, jMol, jShear, jBrown, strip, dDelta, dFw, Ueff };
}

// =====================================================================================================
// 8. Mineral scale
// =====================================================================================================
const IONS = { Na: { z: 1, M: 22.99, a: 4.0, b: 0.075 }, K: { z: 1, M: 39.098, a: 3.5, b: 0.015 }, Ca: { z: 2, M: 40.078, a: 5.0, b: 0.165 }, Mg: { z: 2, M: 24.305, a: 5.5, b: 0.2 }, Ba: { z: 2, M: 137.327, a: 5.0, b: 0 }, Sr: { z: 2, M: 87.62, a: 5.26, b: 0.121 }, Fe: { z: 2, M: 55.845, a: 6.0, b: 0 }, Cl: { z: -1, M: 35.453, a: 3.5, b: 0.015 }, SO4: { z: -2, M: 96.06, a: 5.0, b: -0.04 }, HCO3: { z: -1, M: 61.017, a: 5.4, b: 0 } };
export const ION_IDS = Object.freeze(Object.keys(IONS));
/** Standard seawater (mg/L) used for the injection-water mixing curve. */
export const SEAWATER = Object.freeze({ Na: 10781, K: 399, Ca: 412, Mg: 1284, Ba: 0.02, Sr: 7.9, Fe: 0.003, Cl: 19353, SO4: 2712, HCO3: 142 });
// log10 Ksp at 25 °C, ΔH (J/mol), ΔV of dissolution (cm³/mol), molar mass (g/mol), density (kg/m³), cation, anion
const MINERALS = [
  { id: 'calcite', name: 'Calcite (CaCO₃)', logK: -8.48, dH: -9610, dV: -58.4, M: 100.09, rho: 2710, cat: 'Ca', an: 'CO3' },
  { id: 'barite', name: 'Barite (BaSO₄)', logK: -9.97, dH: 26570, dV: -50.6, M: 233.39, rho: 4480, cat: 'Ba', an: 'SO4' },
  { id: 'celestite', name: 'Celestite (SrSO₄)', logK: -6.63, dH: -4340, dV: -49.7, M: 183.68, rho: 3960, cat: 'Sr', an: 'SO4' },
  { id: 'gypsum', name: 'Gypsum (CaSO₄·2H₂O)', logK: -4.58, dH: -456, dV: -42.4, M: 172.17, rho: 2320, cat: 'Ca', an: 'SO4' },
  { id: 'anhydrite', name: 'Anhydrite (CaSO₄)', logK: -4.36, dH: -7150, dV: -49.8, M: 136.14, rho: 2960, cat: 'Ca', an: 'SO4' },
  { id: 'siderite', name: 'Siderite (FeCO₃)', logK: -10.89, dH: -10380, dV: -55.3, M: 115.85, rho: 3870, cat: 'Fe', an: 'CO3' },
];
/** Ionic strength (mol/L) and molar concentrations of a water analysis given in mg/L: { I, m: { ion: mol/L }, tds (mg/L), balance (charge-balance error, fraction) }. */
export function ionicStrength(water) {
  const m = {}; let I = 0, tds = 0, cat = 0, an = 0;
  for (const id of ION_IDS) { const c = Math.max(num(water?.[id], 0), 0), ion = IONS[id]; m[id] = c / (ion.M * 1000); I += 0.5 * m[id] * ion.z * ion.z; tds += c; if (ion.z > 0) cat += m[id] * ion.z; else an -= m[id] * ion.z; }
  return { I, m, tds, balance: cat + an > 0 ? (cat - an) / (cat + an) : 0 };
}
/**
 * Single-ion activity coefficient. model: 'davies' | 'truesdellJones' (extended Debye–Hückel with ion-size and salting-out
 * parameters, usable to about 2 mol/L in chloride waters). ion: optional { a (Å), b }.
 */
export function activityCoefficient(z, I, Tc = 25, model = 'davies', ion = null) {
  const T = Tc + KEL, eps = 87.74 - 0.4008 * Tc + 9.398e-4 * Tc * Tc - 1.41e-6 * Tc ** 3, rho = 1 - 4.5e-6 * (Tc - 4) ** 2 * (Tc < 100 ? 1 : 1), A = (1.82483e6 * Math.sqrt(Math.max(rho, 0.9))) / (eps * T) ** 1.5, B = (50.2916 * Math.sqrt(Math.max(rho, 0.9))) / Math.sqrt(eps * T), s = Math.sqrt(I);
  if (model === 'truesdellJones' && ion) return 10 ** ((-A * z * z * s) / (1 + B * ion.a * s) + ion.b * I);
  return 10 ** (-A * z * z * (s / (1 + s) - 0.3 * I));
}
/** Mix two water analyses (mg/L) with a volume fraction f of the second one. */
export const mixWaters = (a, b, f) => Object.fromEntries(ION_IDS.map((id) => [id, (1 - f) * num(a?.[id], 0) + f * num(b?.[id], 0)]));
/**
 * Saturation indices of the common oilfield scales at T (°C), P (bara).
 * opt: { yCO2 (mole fraction of CO2 in the gas), model ('truesdellJones' | 'davies') }.
 * Carbonate system: pH from the CO2 fugacity and the bicarbonate alkalinity with temperature-dependent Henry and dissociation
 * constants (Plummer–Busenberg); solubility products by van 't Hoff with a volume-of-reaction pressure term.
 * Returns { I, tds, pH, fCO2 (bar), minerals: [{ id, name, SI, ptb (mg/L that can precipitate), logK }], max: { id, name, SI }, oddoTomson (calcite Is, Oddo–Tomson) }.
 */
export function scaleIndices(water, Tc, Pbar, { yCO2 = 0.03, model = 'truesdellJones' } = {}) {
  const { I, m, tds } = ionicStrength(water), T = Tc + KEL, lg = Math.log10, g = (id) => activityCoefficient(IONS[id].z, I, Tc, model, IONS[id]);
  const TF = Tc * 1.8 + 32, Ppsia = Pbar * 14.5038, phi = Math.exp(Ppsia * (2.84e-4 - 0.255 / (TF + 460))), fCO2 = Math.max(yCO2 * Pbar * phi, 1e-9);
  const lKH = 108.3865 + 0.01985076 * T - 6919.53 / T - 40.45154 * lg(T) + 669365 / T ** 2, lK1 = -356.3094 - 0.06091964 * T + 21834.37 / T + 126.8339 * lg(T) - 1684915 / T ** 2, lK2 = -107.8871 - 0.03252849 * T + 5151.79 / T + 38.92561 * lg(T) - 563713.9 / T ** 2;
  const gam = Object.fromEntries(ION_IDS.map((id) => [id, g(id)])), g2 = activityCoefficient(2, I, Tc, model, { a: 5.4, b: 0 });
  const aHCO3 = Math.max(gam.HCO3 * m.HCO3, 1e-12), aH = (10 ** (lK1 + lKH) * (fCO2 / 1.01325)) / aHCO3, pH = -lg(aH), aCO3 = (10 ** lK2 * aHCO3) / aH;
  const logKof = (mn) => (mn.id === 'calcite' ? -171.9065 - 0.077993 * T + 2839.319 / T + 71.595 * lg(T) : mn.logK - (mn.dH / (2.302585 * R)) * (1 / T - 1 / 298.15)) - (mn.dV * 1e-6 * (Pbar - 1) * 1e5) / (2.302585 * R * T);
  const minerals = MINERALS.map((mn) => {
    const logK = logKof(mn), K = 10 ** logK, carb = mn.an === 'CO3', aCat = gam[mn.cat] * m[mn.cat], SI = aCat > 0 && (carb ? aCO3 : m.SO4) > 0 ? lg((aCat * (carb ? aCO3 : gam.SO4 * m.SO4)) / K) : -99;
    let x = 0;
    if (SI > 0) {
      if (carb) { // M²⁺ + 2 HCO3⁻ → MCO3 + CO2 + H2O at constant CO2 fugacity
        const Kc = (K * 10 ** (lK1 + lKH) * (fCO2 / 1.01325)) / 10 ** lK2, fn = (y) => (m[mn.cat] - y) * (m.HCO3 - 2 * y) ** 2 * gam[mn.cat] * gam.HCO3 ** 2 - Kc, hi = Math.min(m[mn.cat], m.HCO3 / 2);
        x = fn(hi) >= 0 ? hi : brent(fn, 0, hi, 1e-14);
      } else { const a = m[mn.cat], b = m.SO4; x = 0.5 * (a + b - Math.sqrt((a - b) ** 2 + (4 * K) / (gam[mn.cat] * gam.SO4))); }
    }
    return { id: mn.id, name: mn.name, SI, ptb: Math.max(x, 0) * mn.M * 1000, logK, M: mn.M, rho: mn.rho };
  });
  const max = minerals.reduce((a, b) => (b.SI > a.SI ? b : a));
  // Oddo–Tomson calcite index (gas phase present), conditional constants in °F, psia and mol/L
  const oddoTomson = m.Ca > 0 && m.HCO3 > 0 ? lg((m.Ca * m.HCO3 ** 2) / ((fCO2 / phi) * 14.5038 * phi)) + 5.85 + 15.19e-3 * TF - 1.64e-6 * TF * TF - 5.27e-5 * Ppsia - 3.334 * Math.sqrt(I) + 1.431 * I : -99;
  return { I, tds, pH, fCO2, gamma: gam, gammaCO3: g2, minerals, max: { id: max.id, name: max.name, SI: max.SI }, oddoTomson };
}

// =====================================================================================================
// 9. Asphaltene
// =====================================================================================================
/**
 * de Boer screening: undersaturation (reservoir pressure − saturation pressure, bar) against the in-situ oil density (kg/m³).
 * The two boundaries are straight-line fits to the published plot and are for screening only.
 * Returns { cls: 'no problems' | 'slight problems' | 'severe problems', lower, upper (bar at this density) }.
 */
export function deBoer(rho, dP) {
  const lower = Math.max(0, 70 + 1.3 * (rho - 600)), upper = lower + 130;
  return { cls: dP <= lower ? 'no problems' : dP <= upper ? 'slight problems' : 'severe problems', lower, upper };
}
/** Colloidal instability index from a SARA analysis (wt %): (saturates + asphaltenes) / (aromatics + resins). */
export function colloidalInstability({ sat, aro, res, asp }) {
  const cii = (sat + asp) / Math.max(aro + res, 1e-9);
  return { cii, cls: cii < 0.7 ? 'stable' : cii <= 0.9 ? 'uncertain' : 'unstable' };
}
/**
 * Flory–Huggins (Hirschberg) maximum volume fraction of asphaltene soluble in a live oil.
 * { rhoL (kg/m³), mwL (g/mol), TK, deltaA (MPa^0.5 at 25 °C), vA (m³/kmol) }. The oil solubility parameter follows
 * δ = 17.347 ρ(g/cm³) + 2.904 MPa^0.5. Returns { phiMax, deltaL, deltaA }.
 */
export function asphalteneSolubility({ rhoL, mwL, TK, deltaA = 20.0, vA = 2.0 }) {
  const deltaL = 17.347 * (rhoL / 1000) + 2.904, dA = deltaA * (1 - 1.07e-3 * (TK - 298.15)), vL = mwL / rhoL; // m³/kmol
  const ex = (vA / vL) * (1 - vL / vA) - (vA * 1e-3 * ((dA - deltaL) * 1e3) ** 2) / (R * TK); // δ in Pa^0.5, v in m³/mol
  return { phiMax: Math.min(Math.exp(Math.min(ex, 0)), 1), deltaL, deltaA: dA };
}

// =====================================================================================================
// 10. Sand
// =====================================================================================================
/**
 * Minimum transport (critical) velocity of sand in a mostly horizontal pipe by three correlations.
 * { d (m), D (m), rhoP, rhoF (carrier liquid), mu (Pa·s), C (sand volume fraction), vsl, vm (m/s, for Salama's liquid fraction) }
 * Returns { oroskarTurian, salama, danielson (m/s; Salama and Danielson are mixture velocities), governing (the largest), settling }.
 */
export function sandCriticalVelocity({ d, D, rhoP = 2650, rhoF = 800, mu = 2e-3, C = 1e-4, vsl = 1, vm = 2 }) {
  const s = rhoP / rhoF, nu = mu / rhoF, c = clamp(C, 1e-7, 0.5), root = Math.sqrt(G * d * Math.max(s - 1, 1e-6));
  const oroskarTurian = root * 1.85 * c ** 0.1536 * (1 - c) ** 0.3564 * (d / D) ** -0.378 * ((D * rhoF * root) / mu) ** 0.09 * 0.96 ** 0.3;
  const salama = clamp(vsl / Math.max(vm, 1e-9), 0.01, 1) ** 0.53 * d ** 0.17 * nu ** -0.09 * Math.max(s - 1, 1e-6) ** 0.55 * D ** 0.47;
  const danielson = 0.23 * nu ** (-1 / 9) * d ** (1 / 9) * (G * D * Math.max(s - 1, 1e-6)) ** (5 / 9);
  return { oroskarTurian, salama, danielson, governing: Math.max(oroskarTurian, salama, danielson), settling: settlingVelocity(d, rhoP, rhoF, mu).v };
}
/** Screening erosion rate (mm/y) in a bend from sand: Salama (2000), E = W·V²·d / (Sm·D²·ρm) with W kg/d, d µm, D mm, Sm = 5.5. */
export const sandErosionScreen = (Wkgd, vm, dUm, Dmm, rhoM) => (Wkgd * vm * vm * dUm) / (5.5 * Dmm * Dmm * Math.max(rhoM, 1));

// =====================================================================================================
// 11. Case set-up on the solids grid
// =====================================================================================================
const WETGAS = { name: 'Lean wet gas', comp: { N2: 1, CO2: 2, H2S: 0, C1: 86, C2: 6, C3: 3, iC4: 0.5, nC4: 0.8, iC5: 0.2, nC5: 0.2, C6: 0.2, C7p: 0.1 }, c7MW: 120, c7SG: 0.76, rateBasis: 'gas', wc: 0 };
const C_STEEL = 470; // J/kg/K

/**
 * Build the line, fluid and reference flow pictures on n equal cells. Returns the set-up object used by the marching
 * solvers: geometry arrays, props(P, T) (memoised kernel properties at the case rate), prof(frac) (steady picture at a
 * rate fraction), grad(...) (kernel pressure gradient), the hydrate curve and its inverse, thermal mass per metre.
 */
export function buildSetup(v, ctx, n) {
  const c = { fluid: ctx?.fluid, site: ctx?.site, outputs: ctx?.outputs };
  const fo = v.fluidSystem === 'wetgas' ? { ...WETGAS, qGas: v.gasRate, qWater: v.gasWater } : v.fluidSystem === 'highwc' ? { rateBasis: 'oil', wc: v.highWc } : null;
  const line0 = caseLine(c), over = { id: v.idMm / 1000, roughness: v.roughUm * 1e-6, uValue: v.uValue, tSeabed: v.tSeabed };
  const differs = (a, b) => Math.abs(a - b) > 1e-6 * Math.max(Math.abs(b), 1e-9);
  const custom = !!fo || differs(over.id, line0.id) || differs(over.roughness, line0.roughness) || differs(over.uValue, line0.uValue) || differs(over.tSeabed, line0.tSeabed) || differs(v.lengthScale, 1);
  if (differs(v.lengthScale, 1)) over.profile = { x: line0.profile.x.map((x) => x * v.lengthScale), z: line0.profile.z.slice() };
  if (fo) over.override = fo;
  const fm = fluidModel(c, fo || {}), line = caseLine(c, over), D0 = line.id;
  const toCells = (pic) => {
    const sp = [0]; for (let i = 1; i < pic.x.length; i++) sp.push(sp[i - 1] + Math.max(Math.hypot(pic.x[i] - pic.x[i - 1], pic.z[i] - pic.z[i - 1]), 1e-9));
    const L = sp[sp.length - 1], ds = L / n, sc = Array.from({ length: n }, (_, i) => (i + 0.5) * ds), f = (a) => sc.map((s) => interp1(sp, a, s)), ze = Array.from({ length: n + 1 }, (_, i) => interp1(sp, pic.z, i * ds));
    return { L, ds, s: sc, x: f(pic.x), z: f(pic.z), theta: ze.slice(1).map((z, i) => Math.asin(clamp((z - ze[i]) / ds, -1, 1))), P: f(pic.P), T: f(pic.T), holdup: f(pic.holdup).map((h) => clamp(h, 0.01, 1)), dpdx: f(pic.dpdx), tauW: f(pic.tauW).map((t) => Math.max(t, 0)), rhoM: f(pic.rhoM), tAmb: f(pic.tAmb),
      regime: sc.map((s) => { let k = 0; while (k < sp.length - 1 && sp[k + 1] < s) k++; return String(pic.regime[Math.min(k, pic.regime.length - 1)] ?? ''); }), pOut: pic.P[pic.P.length - 1], pIn: pic.P[0], src: pic.source || 'kernel estimate' };
  };
  const memo = new Map();
  const prof = (frac) => {
    const key = Math.round(frac * 1e4);
    if (!memo.has(key)) {
      let pic;
      try { pic = key === 10000 ? flowPicture(c, custom ? { ...over, force: true } : {}) : steadyCase(c, { ...over, mScale: frac }); }
      catch (e) { throw new Error(`No steady flow solution at ${(frac * 100).toFixed(0)} % of the case rate: ${e.message}`); }
      memo.set(key, toCells(pic));
    }
    return memo.get(key);
  };
  const base = prof(1), pc = new Map();
  const props = (P, T) => {
    const kp = Math.round(48 * Math.log(clamp(P, 1, 900))), kt = Math.round(clamp(T, -30, 170)), key = kt * 4096 + kp;
    let o = pc.get(key); if (!o) { o = fm.at(Math.exp(kp / 48), kt); pc.set(key, o); } return o;
  };
  const aqS = fm.aq.S, inhId = v.inhibitor === 'case' ? fm.aq.inhId : v.inhibitor, inh = INHIBITORS[inhId] || INHIBITORS.none, inhWt = inhId === 'none' ? 0 : v.inhibitor === 'case' ? fm.aq.inhWt : v.inhWt;
  const Pg = logspace(1, 700, 90), Tg = Pg.map((p) => fm.hydrateT0(p)); for (let i = 1; i < Tg.length; i++) if (Tg[i] <= Tg[i - 1]) Tg[i] = Tg[i - 1] + 1e-9;
  const lnPg = Pg.map(Math.log), wt = line.wt;
  const S = { n, ds: base.ds, L: base.L, s: base.s, x: base.x, z: base.z, theta: base.theta, tAmb: base.tAmb, D0, rough0: line.roughness, U: line.uValue, wt, fm, line, props, prof, custom, src: base.src,
    hT0: (P) => fm.hydrateT0(P), peq: (Tfresh) => Math.exp(interp1(Tg, lnPg, Tfresh)), aq: { S: aqS, inh, inhId, inhWt }, pOut: base.pOut };
  S.grad = (i, P, T, D, rough, muFac, frac) => { const pr = props(P, T), A = (PI * D * D) / 4; return gradient({ vsl: (pr.qL * frac) / A, vsg: (pr.qG * frac) / A, rhoL: pr.rhoL, rhoG: pr.rhoG, muL: pr.muL * muFac, muG: pr.muG, sigma: pr.sigma, D, theta: S.theta[i], rough, P: P * 1e5 }); };
  const gm = new Map();
  S.gref = (frac) => { const key = Math.round(frac * 1e4); if (!gm.has(key)) { const r = prof(frac); gm.set(key, r.P.map((p, i) => S.grad(i, p, r.T[i], D0, S.rough0, 1, frac))); } return gm.get(key); };
  S.C = base.P.map((p, i) => { if (v.thermalMass > 0) return v.thermalMass * 1000; const pr = props(p, base.T[i]), A = (PI * D0 * D0) / 4, H = base.holdup[i]; return A * (H * pr.rhoL * pr.cpL + (1 - H) * pr.rhoG * pr.cpG) + BASE.rhoSteel * C_STEEL * PI * (D0 + wt) * wt; });
  return S;
}
/** Uniform laboratory-style set-up with constant fluid properties (used by the verification cases and the calibration model). */
export function labSetup({ n = 20, L = 2000, D = 0.1, T = 4, P = 80, U = 0, tAmb = 4, vsl = 1, vsg = 1, wcut = 0.2, sg = 0.7, C = 6e4 } = {}) {
  const A = (PI * D * D) / 4, qL = vsl * A, qG = vsg * A, rhoO = 800, rhoW = 1000, rhoG = 70, qW = qL * wcut, qO = qL - qW, rhoL = (qO * rhoO + qW * rhoW) / qL;
  const pr = { rhoG, rhoO, rhoW, rhoL, muG: 1.3e-5, muO: 3e-3, muW: 1.5e-3, muL: 3e-3 * (1 - Math.min(wcut, 0.7)) ** -2.5, cpG: 2500, cpO: 2100, cpW: 4100, cpL: (qO * rhoO * 2100 + qW * rhoW * 4100) / (qL * rhoL), kG: 0.04, kO: 0.14, kW: 0.58, kL: 0.2, sigma: 0.02, sigmaOW: 0.03, zG: 0.85, mwG: 18, mwO: 150, mG: qG * rhoG, mO: qO * rhoO, mW: qW * rhoW, qG, qO, qW, qL, wcut, phaseInv: wcut > 0.6 };
  const ds = L / n, arr = (x) => new Array(n).fill(x), Pg = logspace(1, 700, 90), Tg = Pg.map((p) => hydrateT0(p, sg)), lnPg = Pg.map(Math.log);
  const S = { n, ds, L, s: arr(0).map((_, i) => (i + 0.5) * ds), x: arr(0).map((_, i) => (i + 0.5) * ds), z: arr(0), theta: arr(0), tAmb: arr(tAmb), D0: D, rough0: 4.5e-5, U, wt: 0.008, props: () => pr, custom: true, src: 'laboratory set-up',
    hT0: (p) => hydrateT0(p, sg), peq: (Tf) => Math.exp(interp1(Tg, lnPg, Tf)), aq: { S: 0, inh: INHIBITORS.none, inhId: 'none', inhWt: 0 }, pOut: P, C: arr(C) };
  S.grad = (i, Pp, Tt, Dd, rough, muFac, frac) => { const a = (PI * Dd * Dd) / 4; return gradient({ vsl: (qL * frac) / a, vsg: (qG * frac) / a, rhoL, rhoG, muL: pr.muL * muFac, muG: pr.muG, sigma: pr.sigma, D: Dd, theta: 0, rough, P: Pp * 1e5 }); };
  const g1 = S.grad(0, P, T, D, 4.5e-5, 1, 1), pm = new Map();
  S.prof = (frac) => { const k = Math.round(frac * 1e4); if (!pm.has(k)) { const g = S.grad(0, P, T, D, 4.5e-5, 1, frac); pm.set(k, { L, ds, P: arr(P), T: arr(T), holdup: arr(g.holdup), dpdx: arr(0), tauW: arr(g.tauW), rhoM: arr(rhoL * g.holdup + rhoG * (1 - g.holdup)), tAmb: arr(tAmb), regime: arr(g.regime || 'slug'), pOut: P, pIn: P }); } return pm.get(k); };
  S.gref = (frac) => arr(S.grad(0, P, T, D, 4.5e-5, 1, frac)); S.g1 = g1;
  return S;
}

// =====================================================================================================
// 12. Transient hydrate march: transport, kinetics, population balance, deposition, bore feedback
// =====================================================================================================
/** Kinetic and deposition parameters of the hydrate march from the input values (SI). */
export function hydrateParams(v, S, over = {}) {
  const b = S.prof(1), mid = Math.floor(S.n / 2), pr = S.props(b.P[mid], b.T[mid]), lamL = pr.qL / Math.max(pr.qL + pr.qG, 1e-12);
  const mode = v.regime && v.regime !== 'auto' ? v.regime : lamL < 0.1 ? 'gas' : pr.wcut > 0.6 ? 'water' : 'oil';
  // primary particle: Boxall inertial droplet correlation d/D = 0.063 We^-3/5 on the mixture velocity
  const vm = (pr.qL + pr.qG) / ((PI * S.D0 * S.D0) / 4), We = ((mode === 'water' ? pr.rhoW : pr.rhoO) * vm * vm * S.D0) / Math.max(mode === 'oil' ? pr.sigmaOW : pr.sigma, 1e-4);
  const dAuto = clamp(0.063 * We ** -0.6 * S.D0, 10e-6, 400e-6), dPrim = v.primaryUm > 0 ? v.primaryUm * 1e-6 : dAuto;
  let hydN = v.hydNumber;
  if (!(hydN > 0)) { const Tm = mean(S.tAmb) + KEL, Pm = mean(b.P); hydN = clamp(langmuirOccupancy(Math.max(Tm, 273.2), methaneFugacity(Pm, Math.max(Tm, 273.2))).hydrationNumber, 5.75, 7.5); }
  const mwG = clamp(pr.mwG || 18, 16, 30) * 1e-3, mHyd = mwG + hydN * MW_W;
  const rhoH = v.rhoHyd > 0 ? v.rhoHyd : HYDRATE.rho;
  return { mode, dPrim, dAuto, hydN, mwG, rhoH, wfH: (hydN * MW_W) / mHyd, vmh: mHyd / rhoH, dHmol: HYDRATE.latent * mHyd,
    kinK: v.kinK * 1e-10, EaR: v.kinEa, shellD: v.shellD * 1e-13, H: mode === 'oil' ? 2500 : 7e4, mtMult: v.mtMult,
    nucA: 10 ** v.nucA, theta: v.contactAngle, sigma: v.sigmaHW * 1e-3, nucV: v.nucVolume * 1e-3, het: v.nucleation !== 'homogeneous', lamStar: 1,
    cohesion: v.cohesion * 1e-3, fr: v.fractal, phiMax: v.phiMax, viscModel: v.viscModel, aggEff: v.aggEff, kBreak: v.kBreak,
    inhEff: v.inhEff / 100, htMult: v.htMult, inletPhi: v.inletHydPct / 100,
    adhesion: v.adhesion, adhForce: v.adhForce * 1e-3, tauCrit: v.tauCrit, kRemove: v.kRemove / 3600, por0: v.porosity0, porInf: Math.min(v.porosityInf, v.porosity0), tAge: v.ageHours * 3600, filmMult: v.filmMult,
    plugBlock: v.plugBlockPct / 100, pInMax: v.pInMax, pShut: v.pShut, maxSub: 4, ...over };
}
/**
 * Transient hydrate solver on the line (a stepper: call step() until done, then result()).
 * Fluid-borne quantities (particle numbers per size class, primary-particle count, water removed to deposits, the
 * nucleation hazard integral) are carried as cell inventories with implicit upwind transport; in each cell the kinetics
 * act on the mixed state: nucleation hazard → onset → growth (resistances in series, limited by water, gas and the heat
 * that can be removed) → aggregation and breakage (sectional population balance) → wall capture, vapour-film growth,
 * shear removal and ageing of the deposit. The bore, the slurry viscosity and the wall shear feed back on the hydraulics.
 * o: { phases: [{ dur (s), dt (s), frac }], grid (pbeGrid), Dbase[] (bore before hydrate, m), roughBase[], dep0 (m), fracInit, cheap, rows }
 */
export function hydrateMarch(S, p, o) {
  const n = S.n, ds = S.ds, g = o.grid, nC = g.n, D0 = S.D0, RHO = p.rhoH, LAT = HYDRATE.latent, vp = (PI / 6) * p.dPrim ** 3, eexp = (S.props(S.pOut, 4).rhoW / RHO) / p.wfH; // hydrate volume per unit water volume
  const NP = Array.from({ length: n }, () => new Float64Array(nC)), PS = new Float64Array(n), WD = new Float64Array(n), LV = new Float64Array(n), HV = new Float64Array(n), film = new Float64Array(n), mDep = new Float64Array(n), por = new Float64Array(n).fill(p.por0);
  const Dbase = o.Dbase || new Array(n).fill(D0), roughBase = o.roughBase || new Array(n).fill(S.rough0), init = S.prof(o.fracInit ?? 1);
  for (let i = 0; i < n; i++) if (o.dep0 > 0) mDep[i] = RHO * (1 - p.por0) * (PI / 4) * (Dbase[i] ** 2 - Math.max(Dbase[i] - 2 * o.dep0, 0.05 * D0) ** 2) * ds;
  const mDep0 = mDep.reduce((a, b) => a + b, 0);
  const T = init.T.slice(), P = init.P.slice(), hold = init.holdup.slice(), dpdx = init.dpdx.slice(), tauW = init.tauW.slice(), muFac = new Array(n).fill(1), Dn = Dbase.slice();
  const cD = new Array(n).fill(-1), cMu = new Array(n).fill(1), cT = new Array(n).fill(0), cP = new Array(n).fill(0), cF = new Array(n).fill(0), cG = new Array(n).fill(null);
  const rec = { sub: new Array(n).fill(0), teq: new Array(n).fill(0), phi: new Array(n).fill(0), phiE: new Array(n).fill(0), d43: new Array(n).fill(0), dA: new Array(n).fill(0), J: new Array(n).fill(0), lam: new Array(n).fill(0), rate: new Array(n).fill(0), vL: new Array(n).fill(0), X: new Array(n).fill(0), tw: new Array(n).fill(0), lim: new Array(n).fill(''), depRate: new Array(n).fill(0), freeW: new Array(n).fill(1), subMax: new Array(n).fill(-99), phiMaxT: new Array(n).fill(0), Jmax: 0, gdot: new Array(n).fill(0), muC: new Array(n).fill(1e-3), rhoC: new Array(n).fill(800), exposure: new Array(n).fill(0), removal: new Array(n).fill(0), capture: new Array(n).fill(0) };
  const ser = { t: [], pIn: [], dp: [], susp: [], dep: [], blk: [], phi: [], rate: [], sub: [], frac: [], visc: [] }, fld = { t: [], phi: [], dep: [], sub: [] };
  const cpC = new Array(n).fill(null), L3 = g.L.map((L) => L ** 3), L4 = g.L.map((L) => L ** 4);
  const beta = new Float64Array(nC * nC), gv = new Float64Array(nC), Sb = new Float64Array(nC), FinN = new Float64Array(nC);
  const peak = { phi: 0, blk: 0, tPhi: 0, tBlk: 0, iPhi: 0, iBlk: 0, N: null, phiX: null, dHyd: null };
  const led = { inflow: 0, formed: 0, dissociated: 0, exported: 0, sloughed: 0, captured: 0, filmWall: 0, heat: 0 };
  const total = o.phases.reduce((a, b) => a + b.dur, 0), nSteps = o.phases.reduce((a, b) => a + Math.ceil(b.dur / b.dt - 1e-9), 0), stride = Math.max(1, Math.ceil(nSteps / (o.rows || 48)));
  let ph = 0, tPh = 0, t = 0, k = 0, done = nSteps === 0, ref = null, gRef = null, fe = 0, pst = null, tm0 = 1, onset = null, plug = null, limited = false, pInNow = P[0], released = true, newPhase = true;
  const dep0Of = (cf) => { const Sx = Math.min(S.aq.S * cf, 260), salt = hydrateDepression({ S: Sx, inhWt: 0, inh: S.aq.inh }); if (!(S.aq.inhWt > 0)) return salt; const w = S.aq.inhWt / 100; return salt + (p.inhEff ?? 1) * (hydrateDepression({ S: Sx, inhWt: (100 * w) / (w + (1 - w) / cf), inh: S.aq.inh }) - salt); };
  const depBase = dep0Of(1), kStarAt = (TK) => p.kinK * Math.exp(-p.EaR * (1 / TK - 1 / 277.15));
  const kP = (() => { let b = 0; for (let j = 1; j < nC; j++) if (Math.abs(Math.log(g.L[j] / p.dPrim)) < Math.abs(Math.log(g.L[b] / p.dPrim))) b = j; return b; })();
  const kSl = (L) => { let b = 0; for (let j = 1; j < nC; j++) if (Math.abs(Math.log(g.L[j] / L)) < Math.abs(Math.log(g.L[b] / L))) b = j; return b; };

  function step() {
    if (done) return;
    const phz = o.phases[ph], dt = Math.min(phz.dt, phz.dur - tPh), flowing = phz.frac > 0;
    if (newPhase) {
      newPhase = false;
      if (flowing) { fe = phz.frac; ref = S.prof(phz.frac); gRef = S.gref(phz.frac); released = false; cD.fill(-1); for (let i = 0; i < n; i++) if (o.cheap) { hold[i] = ref.holdup[i]; } }
      else { // static pressure: outlet pressure plus the head of the settled column, later scaled with the absolute gas temperature
        const b = S.prof(1), rm = mean(b.rhoM); pst = new Array(n); pst[n - 1] = p.pShut > 0 ? p.pShut : P[n - 1];
        for (let i = n - 2; i >= 0; i--) pst[i] = pst[i + 1] + (rm * G * (S.z[i + 1] - S.z[i])) / 1e5;
        tm0 = mean(T) + KEL;
      }
    }
    // ---- pressure along the line
    if (flowing) {
      P[n - 1] = S.pOut + (dpdx[n - 1] * ds) / 2e5; for (let i = n - 2; i >= 0; i--) P[i] = P[i + 1] + ((dpdx[i] + dpdx[i + 1]) * ds) / 2e5;
      pInNow = P[0] + (dpdx[0] * ds) / 2e5;
      if (pInNow > p.pInMax) { // the source cannot push harder: the rate falls and the profile is capped at the available pressure
        const r = (p.pInMax - S.pOut) / (pInNow - S.pOut); fe *= clamp(Math.sqrt(Math.max(r, 0)), 0.3, 1);
        for (let i = 0; i < n; i++) P[i] = S.pOut + (P[i] - S.pOut) * r; pInNow = p.pInMax;
        if (fe < 0.05 * phz.frac && !plug) { let im = 0; for (let i = 1; i < n; i++) if (dpdx[i] > dpdx[im]) im = i; plug = { t: t + dt, x: S.x[im], i: im, mech: 'flow stalled at the available inlet pressure' }; }
      }
      for (let i = 0; i < n; i++) P[i] = clamp(P[i], 1.05, 1500);
    } else { const f = (mean(T) + KEL) / tm0; for (let i = 0; i < n; i++) P[i] = Math.max(pst[i] * f, 1.05); pInNow = P[0]; }
    // ---- march downstream through the cells
    let thUp = 0, FinPS = 0, FinWD = 0, FinLV = 0, FinHV = 0, hydIn = 0, sumRate = 0, sumSusp = 0, maxBlk = 0, iBlk = 0, maxPhi = 0, iPhi = 0, maxVisc = 1, iVisc = 0, maxSub = -99;
    FinN.fill(0);
    if (flowing && p.inletPhi > 0) { const qin = S.props(P[0], T[0]).qL * fe * p.inletPhi; FinN[kP] = qin / g.v[kP]; FinHV = qin; FinPS = qin / (vp * eexp); led.inflow += qin * RHO * dt; }
    for (let i = 0; i < n; i++) {
      const pr0 = S.props(P[i], T[i]), mcp = flowing ? fe * (pr0.mG * pr0.cpG + pr0.mO * pr0.cpO + pr0.mW * pr0.cpW) : 0, Cds = (S.C[i] * ds) / dt, UA = S.U * PI * D0 * ds, den = Cds + mcp + UA;
      let Ti = flowing ? ref.T[i] + ((T[i] - ref.T[i]) * Cds + mcp * thUp) / den : (T[i] * Cds + UA * S.tAmb[i]) / den;
      const pr = S.props(P[i], Ti), Adep = mDep[i] / (RHO * (1 - por[i]) * ds), D = Math.sqrt(Math.max(Dbase[i] ** 2 - (4 * Adep) / PI, (0.03 * D0) ** 2)), A = (PI * D * D) / 4, dHyd = (Dbase[i] - D) / 2;
      Dn[i] = D;
      if (flowing) {
        const mu = muFac[i] * pr.muL;
        if (o.cheap) { const sc = (D0 / D) ** 4.8 * (mu / S.props(ref.P[i], ref.T[i]).muL) ** 0.2 * (fe / phz.frac) ** 1.8, fr0 = Math.max(gRef[i].fric, 1e-6); dpdx[i] = ref.dpdx[i] + fr0 * (sc - 1); tauW[i] = Math.max(ref.tauW[i], gRef[i].tauW) * sc * (D / D0); }
        else {
          // kernel gradient at anchor states; between anchors the friction follows the Blasius scaling in bore, viscosity and rate
          if (cD[i] < 0 || Math.abs(D / cD[i] - 1) > 0.06 || mu / cMu[i] > 2.5 || mu / cMu[i] < 0.4 || Math.abs(Ti - cT[i]) > 8 || Math.abs(P[i] / cP[i] - 1) > 0.15 || Math.abs(fe / cF[i] - 1) > 0.15) {
            const gk = S.grad(i, P[i], Ti, D, Math.min(roughBase[i] + 0.1 * dHyd, 0.05 * D), muFac[i], fe), g0 = gRef[i];
            cD[i] = D; cMu[i] = mu; cT[i] = Ti; cP[i] = P[i]; cF[i] = fe; cG[i] = { dp: ref.dpdx[i] + gk.dpdx - g0.dpdx, fric: Math.max(gk.fric, 0), tau: g0.tauW > 1e-9 ? (Math.max(ref.tauW[i], 1e-9) * gk.tauW) / g0.tauW : gk.tauW };
            hold[i] = clamp(ref.holdup[i] + gk.holdup - g0.holdup, 0.02, 1);
          }
          const sc = (cD[i] / D) ** 4.8 * (mu / cMu[i]) ** 0.2 * (fe / cF[i]) ** 1.8, an = cG[i];
          dpdx[i] = an.dp + an.fric * (sc - 1); tauW[i] = an.tau * sc * (D / cD[i]);
        }
      } else tauW[i] = 0;
      const HL = hold[i], VL = Math.max(A * HL * ds, 1e-12), qL = pr.qL * fe, tau = flowing ? VL / Math.max(qL, 1e-15) : Infinity, a = flowing ? dt / tau : 0, inv = 1 / (1 + a), h = dt * inv, q = 1 + a;
      const vsl = flowing ? qL / A : 0, vsg = flowing ? (pr.qG * fe) / A : 0, vm = vsl + vsg, vL = vsl / HL, wcut = pr.wcut, mW0 = pr.rhoW * wcut * VL;
      const Nc = NP[i], oil = p.mode === 'oil';
      let pVol = 0; // particle volume carried in the cell
      for (let c = 0; c < nC; c++) { Nc[c] = (Nc[c] + dt * FinN[c]) * inv; pVol += Nc[c] * g.v[c]; }
      // hydrate volume: a transported scalar for shelled droplets (oil-continuous), the particle volume itself for crystals
      let hydVol = oil ? (HV[i] + dt * FinHV) * inv : pVol, ps = (PS[i] + dt * FinPS) * inv, wd = (WD[i] + dt * FinWD) * inv, lv = (LV[i] + dt * FinLV) * inv;
      if (flowing && !released && film[i] > 0) { const fv = film[i] / q, c = kSl(Math.max(film[i] / (D * ds), p.dPrim)); Nc[c] += fv / g.v[c]; pVol += fv; ps += fv / (vp * eexp); hydVol += fv; film[i] = 0; } // the shut-in interface film breaks up into the stream
      const volBefore = hydVol, mWfree = Math.max(mW0 - (hydVol + film[i]) * RHO * p.wfH - wd, 0), cf = mW0 > 0 ? mW0 / Math.max(mWfree, 0.03 * mW0) : 1;
      const dep = cf > 1.0005 ? dep0Of(cf) : depBase, Teq = S.hT0(P[i]) - dep, dTs = Teq - Ti, TK = Ti + KEL, zG = pr.zG || 0.85;
      const Jn = dTs > 0 ? nucleationRate({ TK, dT: dTs, TeqK: Teq + KEL, sigma: p.sigma, theta: p.theta, A: p.nucA, het: p.het }).J : 0;
      if (mWfree > 0 && dTs > 0) lv += h * Jn * p.nucV * VL;
      const lam = lv / VL, rhoC = oil ? pr.rhoO : p.mode === 'water' ? pr.rhoW : pr.rhoL, muC = oil ? pr.muO : p.mode === 'water' ? pr.muW : pr.muL, nuC = muC / rhoC, kC = oil ? pr.kO : p.mode === 'water' ? pr.kW : pr.kL;
      const rhoM = pr.rhoL * HL + pr.rhoG * (1 - HL), epsT = (4 * tauW[i] * vm) / (rhoM * D), gdot = Math.max(Math.sqrt(epsT / nuC), (8 * vL) / D, 0.05);
      const Peq = S.peq(Ti + dep), Dg = (7.4e-12 * Math.sqrt(oil ? pr.mwO || 150 : 46.8) * TK) / (muC * 1000 * 37.7 ** 0.6); // methane in the continuous liquid (Wilke–Chang)
      let cpR = 1e9, dVol = 0, limTxt = dTs > 0 ? (lam >= p.lamStar ? 'no free water' : 'induction (no nuclei yet)') : 'outside the hydrate region', depVol = 0, mFilmWall = 0;
      if (!flowing) {
        // shut-in: phases are segregated, hydrate grows as a film at the water interface, fed by diffusion
        const Ai = D * Math.sqrt(Math.max(1 - (2 * HL - 1) ** 2, 0.05)) * ds;
        if (lam >= p.lamStar && dTs > 0 && mWfree > 0) {
          if (!onset) onset = { t: t + dt, x: S.x[i], i };
          const dc = (fug(P[i], zG) - fugAt(Peq, P[i], zG)) / p.H, j = dc > 0 ? dc / (1 / (kStarAt(TK) * p.H) + (0.25 * D * HL + 1e-3) / Dg + film[i] / Ai / p.shellD) : 0;
          dVol = Math.min(j * Ai * p.vmh * dt, (0.98 * mWfree) / (RHO * p.wfH), (0.9 * dTs * den * dt) / (RHO * LAT)); limTxt = 'diffusion through the interface film (shut-in)';
        } else if (dTs < 0 && film[i] > 0) { dVol = -Math.min(film[i], hydrateDissociationRate({ TK, P: P[i], Peq, zG }) * Ai * p.vmh * dt, (0.9 * -dTs * den * dt) / (RHO * LAT)); limTxt = 'dissociating'; }
        film[i] += dVol;
      } else {
        if (lam >= p.lamStar && dTs > 0 && mWfree > 0.02 * mW0) { // onset: the water phase is seeded as primary particles
          const target = (mWfree / pr.rhoW + hydVol / eexp) / vp;
          if (ps < target * 0.999) {
            const add = target - ps;
            if (oil) { Nc[kP] += (add * vp) / g.v[kP]; pVol += add * vp; hydVol += add * vp * eexp * 1e-3; } else { Nc[0] += add; hydVol += add * g.v[0]; pVol = hydVol; }
            ps = target; if (!onset) onset = { t: t + dt, x: S.x[i], i };
          }
        }
        let any = false, dV = 0;
        if (hydVol > 0 && dTs > 0 && mWfree > 0) {
          const dc = (fug(P[i], zG) - fugAt(Peq, P[i], zG)) / p.H, kKin = kStarAt(TK) * p.H;
          let rPart = 0, lim = 'intrinsic kinetics';
          if (dc > 0) {
            if (oil) {
              const X = clamp(hydVol / Math.max(ps * vp * eexp, 1e-300), 0, 1), Sh = 2 + 0.6 * Math.sqrt((gdot * p.dPrim * p.dPrim) / nuC) * Math.cbrt(nuC / Dg);
              const gr = hydrateGrowthRate({ TK, P: P[i], Peq, zG, kRef: p.kinK, EaR: p.EaR, H: p.H, kFilm: (Sh * Dg) / p.dPrim, kShell: shellConductance(X, p.dPrim / 2, p.shellD), hPart: (2 * kC) / p.dPrim, dT: dTs, dHmol: p.dHmol });
              rPart = X < 1 ? gr.j * ps * PI * p.dPrim * p.dPrim : 0; lim = gr.limiting;
            } else {
              for (let c = 0; c < nC; c++) { const L = g.L[c], j = 1 / (1 / (dc / (1 / kKin + L / (2 * Dg))) + p.dHmol / (((2 * kC) / L) * dTs)); gv[c] = PI * L * L * j; rPart += Nc[c] * gv[c]; }
              lim = 'surface growth (intrinsic + film)';
            }
            let s = 1;
            if (p.mode !== 'gas' && rPart > 0) { // gas must first dissolve in the liquid: absorption in series (Skovborg–Rasmussen)
              const kL = 0.4 * (Math.max(epsT, 1e-9) * nuC) ** 0.25 * Math.sqrt(Dg / nuC), aGL = ((D * Math.sqrt(Math.max(1 - (2 * HL - 1) ** 2, 0.05))) / (A * HL)) * (1 + (vm * vm) / (G * D)), rSup = kL * aGL * VL * dc * p.mtMult;
              s = 1 / (1 + rPart / Math.max(rSup, 1e-300)); if (s < 0.5) lim = 'gas absorption (mass transfer)';
            }
            const req = rPart * s * p.vmh * h, capW = (0.98 * mWfree) / (RHO * p.wfH), capG = Math.max(((0.9 * fe * pr.mG - hydIn * (1 - p.wfH)) * dt) / (q * RHO * (1 - p.wfH)), 0), capH = ((p.htMult ?? 1) * 0.9 * dTs * den * dt) / (RHO * LAT * q);
            const cap = Math.min(capW, capG, capH), sc = req > cap ? cap / req : 1;
            if (sc < 1) lim = cap === capH ? 'heat removal (fluid held at the hydrate temperature)' : cap === capW ? 'water availability' : 'gas availability';
            if (oil) dV = req * sc; else { const f = s * sc * p.vmh; for (let c = 0; c < nC; c++) gv[c] *= f; any = req > 0; }
            limTxt = lim;
          }
        } else if (hydVol > 0 && dTs < 0) {
          const jd = hydrateDissociationRate({ TK, P: P[i], Peq, zG }) * p.vmh, capH = (0.9 * -dTs * den * dt) / (RHO * LAT * q);
          if (oil) dV = -Math.min(hydVol, jd * ps * PI * p.dPrim * p.dPrim * h, capH);
          else { let tot = 0; for (let c = 0; c < nC; c++) { gv[c] = -PI * g.L[c] ** 2 * jd; tot -= Nc[c] * gv[c]; } const sc = tot * h > capH ? capH / (tot * h) : 1; for (let c = 0; c < nC; c++) gv[c] *= sc; any = tot > 0; }
          limTxt = 'dissociating';
        }
        if (pVol > 0) {
          const phiNow = hydVol / VL, cc = cpC[i]; let cp = cc && Math.abs(phiNow / cc.phi - 1) < 0.03 && Math.abs(gdot / cc.gdot - 1) < 0.03 && Math.abs(muC / cc.mu - 1) < 0.03 ? cc.cp : null;
          if (!cp) { cp = maxAgglomerateSize({ dp: p.dPrim, Fa: p.cohesion * p.dPrim, mu0: muC, shear: gdot, phi: phiNow, phiMax: p.phiMax, fr: p.fr }); cpC[i] = { phi: phiNow || 1e-300, gdot, mu: muC, cp }; }
          const agg = p.aggEff > 0 && dTs > 0;
          rec.dA[i] = cp.dA; cpR = cp.ratio;
          if (agg) {
            const cS = (p.aggEff * ((8 * vL) / D / 6 + 0.1618 * Math.sqrt(epsT / nuC))) / VL, cDf = (p.aggEff * (PI / 4) * Math.abs(RHO - rhoC) * G) / (18 * muC) / VL, cB = (p.aggEff * 2 * KB * TK) / (3 * muC) / VL;
            for (let c = 0; c < nC * nC; c++) beta[c] = cS * g.gS[c] + cDf * g.gD[c] + cB * g.gB[c];
          }
          { const kb = (p.kBreak * gdot) / cp.dA ** 3; for (let c = 0; c < nC; c++) Sb[c] = g.L[c] > 1.01 * p.dPrim ? kb * L3[c] : 0; } // agglomerates above the cohesive limit break even when nothing sticks any more
          const r = solvePBE(g, Nc, h, { beta: agg ? beta : null, gv: any ? gv : null, S: Sb, maxSub: p.maxSub, frac: 0.3 });
          if (r.limited) limited = true;
          pVol = 0; for (let c = 0; c < nC; c++) { Nc[c] = r.N[c]; pVol += Nc[c] * g.v[c]; }
        }
        hydVol = oil ? hydVol + dV : pVol;
        dVol = hydVol - volBefore;
        if (!oil && dTs < 0 && volBefore > 0) ps *= clamp(hydVol / volBefore, 0, 1);
        if (hydVol < 1e-30 || (dTs < 0 && hydVol < 1e-9 * pVol)) { Nc.fill(0); ps = 0; hydVol = 0; pVol = 0; dVol = -volBefore; }
        // ---- wall: particle capture on a sub-cooled wall, hydrate film from water vapour on the gas-wetted wall
        const hIn = hInside((pr.rhoL * vm * D) / pr.muL, (pr.cpL * pr.muL) / pr.kL, pr.kL, D), Tw = Ti - (S.U * (Ti - S.tAmb[i])) / Math.max(hIn, 1), cold = Tw < Teq, wet = /strat/.test(ref.regime[i]) ? Math.acos(clamp(1 - 2 * HL, -1, 1)) / PI : 1;
        rec.tw[i] = Tw;
        if (cold && pVol > 0 && p.adhesion > 0) {
          const uS = Math.sqrt(tauW[i] / pr.rhoL), geo = (PI * D * wet) / (A * HL); let lost = 0;
          for (let c = 0; c < nC; c++) {
            if (!(Nc[c] > 0)) continue;
            const L = g.L[c], lamK = depositionVelocity(L, RHO, uS, nuC, rhoC, TK) * p.adhesion * Math.min(1, p.adhForce / (8 * Math.max(tauW[i], 1e-6) * L)) * geo, loss = Nc[c] * (1 - Math.exp(-lamK * h));
            Nc[c] -= loss; lost += loss * g.v[c];
          }
          if (lost > 0) { const fP = lost / pVol; depVol = oil ? hydVol * fP : lost; ps *= 1 - fP; wd += depVol * RHO * p.wfH; hydVol -= depVol; pVol -= lost; }
        }
        if (cold && wet < 1 && vsg > 0 && mWfree > 0 && p.filmMult > 0) {
          const conv = pr.rhoG / ((pr.mwG * 1e-3) / VM_STD), dC = (waterContent(P[i], Ti) - waterContent(P[i], Math.max(Tw, -20))) * conv, Dwg = 2.2e-5 * (1.013 / P[i]) * (TK / 273.15) ** 1.75, Sh = 0.023 * ((pr.rhoG * vsg * D) / pr.muG) ** 0.83 * (pr.muG / (pr.rhoG * Dwg)) ** 0.33;
          mFilmWall = Math.min((p.filmMult * Math.max(dC, 0) * ((Sh * Dwg) / D) * PI * D * (1 - wet) * ds * dt) / p.wfH, (0.5 * mWfree * q) / p.wfH); wd += (mFilmWall * p.wfH) / q;
        }
      }
      // ---- ledger, heat of formation, deposit
      const mNew = dVol * RHO * (flowing ? q : 1);
      if (mNew >= 0) led.formed += mNew; else led.dissociated -= mNew;
      led.formed += mFilmWall; led.filmWall += mFilmWall; led.heat += mNew * LAT;
      Ti += (mNew * LAT) / dt / den;
      let m = mDep[i]; rec.removal[i] = 0;
      const tauC = p.tauCrit * ((1 - por[i]) / (1 - p.por0)) ** 2;
      if (m > 0 && flowing && tauW[i] > tauC && p.kRemove > 0) { // shear removal / sloughing back into the stream
        const dm = m * (1 - Math.exp(-p.kRemove * (tauW[i] / tauC - 1) * dt)), c = kSl(clamp(dHyd, p.dPrim, g.L[nC - 1])), fv = dm / (q * RHO);
        rec.removal[i] = dm / dt; m -= dm; Nc[c] += fv / g.v[c]; pVol += fv; ps += fv / (vp * eexp); wd -= (dm * p.wfH) / q; hydVol += fv; led.sloughed += dm;
      }
      if (m > 0 && dTs < 0) { // a deposit outside the hydrate region melts at the Kim–Bishnoi rate on its exposed surface, as fast as heat arrives
        const dm = Math.min(m, hydrateDissociationRate({ TK, P: P[i], Peq, zG }) * p.vmh * RHO * PI * D * ds * dt, (0.9 * Math.max(Teq - Ti, 0) * -1 + 0.9 * Math.max(Ti - Teq, 0)) * den * dt / LAT);
        if (dm > 0) { m -= dm; led.dissociated += dm; led.heat -= dm * LAT; Ti -= (dm * LAT) / dt / den; if (flowing) wd -= (dm * p.wfH) / q; }
      }
      const gain = depVol * RHO * q + mFilmWall;
      if (gain > 0) { const Vo = m / (RHO * (1 - por[i])), Vn = gain / (RHO * (1 - p.por0)); por[i] = (Vo * por[i] + Vn * p.por0) / (Vo + Vn); m += gain; led.captured += depVol * RHO * q; }
      por[i] = p.porInf + (por[i] - p.porInf) * Math.exp(-dt / p.tAge);
      rec.depRate[i] = (m - mDep[i]) / dt; mDep[i] = m;
      // ---- slurry state for the hydraulics of the next step
      const phiH = hydVol / VL; let m3 = 0, m4 = 0; for (let c = 0; c < nC; c++) { m3 += Nc[c] * L3[c]; m4 += Nc[c] * L4[c]; }
      const d43 = m3 > 0 ? m4 / m3 : 0, ratio = Math.max(1, Math.min(d43 / p.dPrim, 1.5 * cpR)), // the cohesive limit caps the size that enters the rheology
        phiE = Math.min(phiH * ratio ** (3 - p.fr), 1);
      muFac[i] = slurryViscosity(phiE, p.viscModel, { phiMax: p.phiMax });
      HV[i] = hydVol; PS[i] = ps; WD[i] = wd; LV[i] = lv; T[i] = Ti;
      if (flowing) { for (let c = 0; c < nC; c++) FinN[c] = Nc[c] / tau; FinPS = ps / tau; FinWD = wd / tau; FinLV = lv / tau; FinHV = hydVol / tau; hydIn = (hydVol * RHO) / tau; thUp = Ti - ref.T[i]; }
      const blk = 1 - (D * D) / (D0 * D0), rate = mNew / dt + mFilmWall / dt;
      rec.gdot[i] = gdot; rec.muC[i] = muC; rec.rhoC[i] = rhoC; if (dTs > 0) rec.exposure[i] += dt / 3600; rec.capture[i] = (depVol * RHO * q + mFilmWall) / dt;
      rec.sub[i] = dTs; rec.teq[i] = Teq; rec.phi[i] = phiH + film[i] / VL; rec.phiE[i] = phiE; rec.d43[i] = d43; rec.J[i] = Jn; if (Jn > rec.Jmax) rec.Jmax = Jn; rec.lam[i] = lam; rec.rate[i] = rate; rec.vL[i] = vL; rec.lim[i] = limTxt; rec.freeW[i] = mW0 > 0 ? mWfree / mW0 : 0;
      rec.X[i] = mW0 > 0 ? clamp(((hydVol + film[i]) * RHO * p.wfH) / mW0, 0, 1) : 0;
      if (dTs > rec.subMax[i]) rec.subMax[i] = dTs; if (rec.phi[i] > rec.phiMaxT[i]) rec.phiMaxT[i] = rec.phi[i];
      sumRate += rate; sumSusp += (hydVol + film[i]) * RHO;
      if (blk > maxBlk) { maxBlk = blk; iBlk = i; } if (rec.phi[i] > maxPhi) { maxPhi = rec.phi[i]; iPhi = i; } if (muFac[i] > maxVisc) { maxVisc = muFac[i]; iVisc = i; } if (dTs > maxSub) maxSub = dTs;
    }
    if (flowing) { led.exported += dt * hydIn; released = true; }
    t += dt; tPh += dt; k++;
    if (!plug && maxBlk >= p.plugBlock) plug = { t, x: S.x[iBlk], i: iBlk, mech: 'wall deposit closes the bore' };
    if (!plug && flowing && maxVisc >= 1000) plug = { t, x: S.x[iVisc], i: iVisc, mech: 'slurry jams (effective solids fraction at the packing limit)' };
    let depTot = 0; for (let i = 0; i < n; i++) depTot += mDep[i];
    if (maxPhi > peak.phi) Object.assign(peak, { phi: maxPhi, tPhi: t, iPhi, N: Array.from(NP[iPhi]), phiX: rec.phi.slice(), phiE: rec.phiE.slice(), d43: rec.d43.slice(), dA: rec.dA.slice(), gdot: rec.gdot[iPhi], muC: rec.muC[iPhi], rhoC: rec.rhoC[iPhi], T: T[iPhi], visc: muFac.slice() });
    if (maxBlk > peak.blk) Object.assign(peak, { blk: maxBlk, tBlk: t, iBlk, dHyd: Dn.map((d, i) => Math.max((Dbase[i] - d) / 2, 0)), por: Array.from(por) });
    ser.t.push(t / 3600); ser.pIn.push(pInNow); ser.dp.push(pInNow - P[n - 1]); ser.susp.push(sumSusp); ser.dep.push(depTot); ser.blk.push(maxBlk); ser.phi.push(maxPhi); ser.rate.push(sumRate); ser.sub.push(maxSub); ser.frac.push(flowing ? fe : 0); ser.visc.push(maxVisc);
    if (k % stride === 0 || plug || tPh >= phz.dur - 1e-9) { fld.t.push(t / 3600); fld.phi.push(rec.phi.slice()); fld.dep.push(Dn.map((d, i) => (Dbase[i] - d) / 2)); fld.sub.push(rec.sub.slice()); }
    if (tPh >= phz.dur - 1e-9) { ph++; tPh = 0; newPhase = true; }
    if (plug || ph >= o.phases.length) done = true;
  }
  return {
    step, get done() { return done; }, get progress() { return total > 0 ? t / total : 1; },
    result() {
      let susp = 0, depTot = 0; for (let i = 0; i < n; i++) { let hv = film[i] + HV[i]; susp += hv * RHO; depTot += mDep[i]; }
      const dHyd = Dn.map((d, i) => Math.max((Dbase[i] - d) / 2, 0));
      return { t: t / 3600, T: T.slice(), P: P.slice(), D: Dn.slice(), dHyd, por: Array.from(por), hold: hold.slice(), tauW: tauW.slice(), dpdx: dpdx.slice(), muFac: muFac.slice(), N: NP.map((a) => Array.from(a)), rec, ser, fld, peak, onset, plug, limited, pIn: pInNow,
        ledger: { ...led, suspended: susp, deposited: depTot, deposit0: mDep0, in: led.formed + mDep0 + led.inflow, out: susp + depTot + led.exported + led.dissociated } };
    },
  };
}
/** Run a hydrate march to completion (synchronous helper for the Monte Carlo loop, the verification cases and scripts). */
export function runHydrateMarch(S, p, o) { const m = hydrateMarch(S, p, o); while (!m.done) m.step(); return m.result(); }
/**
 * Implicit upwind transport of cell inventories E through cells of residence time tau (s) over a step dt with an inlet
 * flux fin (per second): E_i ← (E_i + dt·F_in)/(1 + dt/τ_i), F_out = E_i/τ_i. Returns the outlet flux. Conservative and
 * unconditionally stable; this is the transport operator of hydrateMarch.
 */
export function advectImplicit(E, tau, dt, fin = 0) { let F = fin; for (let i = 0; i < E.length; i++) { E[i] = (E[i] + dt * F) / (1 + dt / tau[i]); F = E[i] / tau[i]; } return F; }

// =====================================================================================================
// 13. Slow deposits in steady production: wax and mineral scale along the line
// =====================================================================================================
const WATER0 = Object.freeze({ Na: 12500, K: 250, Ca: 900, Mg: 150, Ba: 40, Sr: 90, Fe: 5, Cl: 21500, SO4: 15, HCO3: 650 });
const waterOf = (v) => { const row = Array.isArray(v.water) && v.water[0] ? v.water[0] : WATER0; return Object.fromEntries(ION_IDS.map((id) => [id, Math.max(num(row[id], 0), 0)])); };
/**
 * Wax and scale layers after `days` of steady production at a rate fraction.
 * Wax: Fick diffusion on the solubility slope and the wall heat flux, shear dispersion, Brownian flux, shear stripping,
 * ageing and the insulating feedback of the layer on the heat flux and on the bulk temperature downstream.
 * Scale: saturation indices at the local P, T, surface-reaction kinetics in series with ion mass transfer to the wall.
 * Returns { dWax[], Fw[], dScale[], si[] (max SI), mineral[], scaleRate[] (mm/y), ser: { t (d), dMax (mm), mass (kg) }, waxMass, waxRate0 (mm/d), piggingInterval (d) | null, ... }.
 */
export function slowDeposits(S, v, frac = 1, nt = 40) {
  const b = S.prof(frac), n = S.n, D0 = S.D0, ds = S.ds, days = v.depositDays, dt = (days * 86400) / nt, wTot = v.waxContent / 100, F0 = clamp(1 - v.waxOil / 100, 0.03, 0.95);
  const dWax = new Array(n).fill(0), Fw = new Array(n).fill(F0), water = mixWaters(waterOf(v), SEAWATER, clamp(v.swFrac / 100, 0, 1));
  const cells = b.P.map((P, i) => { const pr = S.props(P, b.T[i]), A = (PI * D0 * D0) / 4, vsl = (pr.qL * frac) / A, vsg = (pr.qG * frac) / A, vm = vsl + vsg, H = b.holdup[i];
    return { pr, vL: vsl / H, vm, hIn: hInside((pr.rhoL * vm * D0) / pr.muL, (pr.cpL * pr.muL) / pr.kL, pr.kL, D0), mcp: frac * (pr.mG * pr.cpG + pr.mO * pr.cpO + pr.mW * pr.cpW), gammaW: b.tauW[i] / pr.muL, oilWet: pr.phaseInv ? 0.15 : 1 - 0.5 * pr.wcut, waterWet: pr.phaseInv ? 1 : clamp(pr.wcut, 0, 1) }; });
  const ser = { t: [0], dMax: [0], mass: [0] }, last = new Array(n).fill(null), Tb = b.T.slice();
  let waxRate0 = 0;
  for (let k = 0; k < nt; k++) {
    let cum = 0, dMax = 0, mass = 0;
    for (let i = 0; i < n; i++) {
      const c = cells[i]; Tb[i] = S.tAmb[i] + (b.T[i] - S.tAmb[i]) * Math.exp(Math.min(cum, 3));
      const r = waxDeposition({ Tb: Tb[i], Tamb: S.tAmb[i], U: S.U, hIn: c.hIn, kOil: c.pr.kO, rhoOil: c.pr.rhoO, muOil: c.pr.muO, wat: v.wat, wTot, slope: v.waxSlope, delta: dWax[i], Fw: Fw[i], kDep: v.waxK, D: D0 - 2 * dWax[i], vL: c.vL, gammaW: c.gammaW, rhoMix: b.rhoM[i], regime: b.regime[i], mult: v.waxMult, diffModel: v.waxDiff, wetFrac: c.oilWet });
      if (k === 0 && r.dDelta * 86400e3 > waxRate0) waxRate0 = r.dDelta * 86400e3;
      cum += ((S.U - r.Ueff) * PI * D0 * ds) / Math.max(c.mcp, 1e-9); last[i] = r;
      dWax[i] = Math.min(dWax[i] + r.dDelta * dt, 0.4 * D0); Fw[i] = Math.min(Fw[i] + r.dFw * dt, 0.95);
      if (dWax[i] > dMax) dMax = dWax[i]; mass += 900 * Fw[i] * (PI / 4) * (D0 * D0 - (D0 - 2 * dWax[i]) ** 2) * ds;
    }
    ser.t.push(((k + 1) * dt) / 86400); ser.dMax.push(dMax * 1000); ser.mass.push(mass);
  }
  const lim = v.waxLimitMm, end = ser.dMax[nt]; let pig = null;
  if (end >= lim) { for (let k = 1; k <= nt; k++) if (ser.dMax[k] >= lim) { pig = ser.t[k - 1] + ((lim - ser.dMax[k - 1]) / (ser.dMax[k] - ser.dMax[k - 1] || 1)) * (ser.t[k] - ser.t[k - 1]); break; } }
  else if (end > 1e-6) { const rate = Math.max((end - ser.dMax[Math.floor(nt / 2)]) / (ser.t[nt] - ser.t[Math.floor(nt / 2)]), 1e-9); pig = Math.min(days + (lim - end) / rate, 3650); }
  // gel strength after a cold shutdown: yield stress from the wax precipitated at ambient temperature
  let restartDp = 0, gelLen = 0, tauY = 0;
  for (let i = 0; i < n; i++) { const sol = waxSolubility(S.tAmb[i], v.wat, wTot, v.waxSlope).solid * 100, ty = sol > 0 && S.tAmb[i] < v.wat - v.pourOffset ? v.gelCoef * sol * sol : 0; if (ty > 0) { restartDp += (4 * ty * ds) / D0; gelLen += ds; tauY = Math.max(tauY, ty); } }
  // scale
  const si = [], mineral = [], scaleRate = [], dScale = [], sis = [], opt = { yCO2: v.co2Pct / 100, model: v.actModel };
  for (let i = 0; i < n; i++) {
    const r = scaleIndices(water, b.T[i], b.P[i], opt), c = cells[i], T = b.T[i] + KEL, ReW = (c.pr.rhoL * c.vm * D0) / c.pr.muL, Sc = c.pr.muW / (c.pr.rhoW * 1e-9), km = (0.023 * ReW ** 0.83 * Sc ** 0.33 * 1e-9) / D0;
    let rate = 0;
    for (const m of r.minerals) { if (!(m.SI > 0) || !(c.pr.mW > 0)) continue; const kin = v.scaleK * 1e-8 * Math.exp((-45000 / R) * (1 / T - 1 / 298.15)) * (10 ** (m.SI / 2) - 1) ** 2, mt = km * (m.ptb / m.M), fl = 1 / (1 / Math.max(kin, 1e-300) + 1 / Math.max(mt, 1e-300)); rate += (c.waterWet * fl * m.M * 1e-3) / (m.rho * 0.8); } // m/s of a 20 % porous layer
    si.push(r.max.SI); mineral.push(r.max.SI > 0 ? r.max.name : 'none'); scaleRate.push(rate * 3.156e10); dScale.push(Math.min(rate * days * 86400, 0.2 * D0)); sis.push(r);
  }
  return { dWax, Fw, dScale, si, mineral, scaleRate, sis, ser, waxMass: ser.mass[nt], waxRate0, piggingInterval: pig, last, Tb, restartDp: restartDp / 1e5, gelLen, tauY, water, prof: b, cells };
}

// =====================================================================================================
// 14. Suite inputs
// =====================================================================================================
const SCEN = [{ value: 'restart', label: 'Cold restart after a shutdown' }, { value: 'steady', label: 'Steady production' }, { value: 'turndown', label: 'Turndown (reduced rate)' }, { value: 'shutdown', label: 'Shutdown cooldown (shut-in line)' }];
const flowing = (v) => v.scenario !== 'shutdown', shut = (v) => v.scenario === 'restart' || v.scenario === 'shutdown';
const INPUTS = [
  { group: 'Operating scenario', tab: 'inputs', help: 'In steady production the reference line runs far above the hydrate temperature, so hydrates are assessed for the off-design events (turndown, shutdown, restart) while wax, scale, asphaltene and sand are assessed for steady production over the period since the last pig run.', fields: [
    { key: 'scenario', label: 'Scenario', type: 'select', value: 'restart', options: SCEN, help: 'Cold restart: the line is shut in for the stated hours (cooling towards ambient at the settle-out pressure) and then restarted at a reduced rate.' },
    { key: 'shutHours', label: 'Shut-in duration', unit: 'h', value: 48, min: 0.5, max: 2000, typical: [8, 96], showIf: shut, help: 'Hours without flow before the restart, or the length of the cooldown that is simulated.' },
    { key: 'restartPct', label: 'Restart rate', unit: '% of case rate', value: 50, min: 5, max: 100, typical: [20, 60], showIf: (v) => v.scenario === 'restart' },
    { key: 'turndownPct', label: 'Turndown rate', unit: '% of case rate', value: 30, min: 5, max: 100, typical: [20, 70], showIf: (v) => v.scenario === 'turndown' },
    { key: 'simHours', label: 'Simulated flowing period', unit: 'h', value: 24, min: 0.5, max: 720, typical: [6, 72], showIf: flowing, help: 'Length of the hydrate event that is marched in time after the scenario starts.' },
    { key: 'depositDays', label: 'Production period since the last pig run', unit: 'd', value: 30, min: 1, max: 720, typical: [7, 180], help: 'Wax and scale layers are grown over this period of steady production and are present when the hydrate event starts.' },
    { key: 'pShut', label: 'Topsides pressure during shut-in', unit: 'bara', value: 0, min: 0, max: 600, showIf: shut, help: '0 = the line is held at the flowing arrival pressure. A lower value represents a partly depressurised line.' },
    { key: 'pInMax', label: 'Available inlet pressure', unit: 'bara', value: 250, min: 5, max: 1400, typical: [100, 400], help: 'Highest pressure the wells or the pump can supply at the flowline inlet. When a restriction needs more than this the rate falls and the line eventually stalls.' },
  ] },
  { group: 'Fluid system and line', tab: 'inputs', help: 'Leave these at the linked values to stay consistent with the other suites. Any change makes this suite recompute its own steady flow picture with the shared kernel.', fields: [
    { key: 'fluidSystem', label: 'Fluid system', type: 'select', value: 'case', options: [{ value: 'case', label: 'Case fluid and rates' }, { value: 'wetgas', label: 'Lean wet gas (gas-dominated)' }, { value: 'highwc', label: 'Case fluid at a high water cut' }] },
    { key: 'gasRate', label: 'Gas rate', unit: 'million Sm³/d', value: 3, min: 0.1, max: 30, showIf: (v) => v.fluidSystem === 'wetgas' },
    { key: 'gasWater', label: 'Free water rate', unit: 'Sm³/d', value: 20, min: 0, max: 2000, showIf: (v) => v.fluidSystem === 'wetgas' },
    { key: 'highWc', label: 'Water cut', unit: '%', value: 70, min: 1, max: 95, showIf: (v) => v.fluidSystem === 'highwc' },
    { key: 'idMm', label: 'Inner diameter', unit: 'mm', value: BASE.idMm, min: 25, max: 1500 },
    { key: 'roughUm', label: 'Clean wall roughness', unit: 'µm', value: BASE.roughUm, min: 0.5, max: 3000 },
    { key: 'uValue', label: 'Overall heat-transfer coefficient (on ID)', unit: 'W/m²K', value: BASE.U, min: 0.2, max: 200, typical: [1, 25], help: 'Controls the heat flux to the wall (wax), the cooldown time and how fast heat of formation can leave.' },
    { key: 'tSeabed', label: 'Seabed / ambient temperature', unit: '°C', value: BASE.tSeabed, min: -5, max: 40 },
    { key: 'lengthScale', label: 'Route length multiplier', unit: '–', value: 1, min: 0.1, max: 10, help: '1 = the case route. Larger values stretch the horizontal distance (long tie-back study).' },
    { key: 'thermalMass', label: 'Thermal mass of pipe and contents', unit: 'kJ/m/K', value: 0, min: 0, max: 5000, help: '0 = computed from the fluid in place and the steel wall (insulation storage neglected, which is conservative for cooldown).' },
    { key: 'dep0Mm', label: 'Initial hydrate deposit thickness', unit: 'mm', value: 0, min: 0, max: 100 },
    { key: 'inletHydPct', label: 'Hydrate particles entering at the inlet', unit: 'vol % of liquid', value: 0, min: 0, max: 30 },
  ] },
  { group: 'Thermodynamic inhibition', tab: 'inputs', fields: [
    { key: 'inhibitor', label: 'Inhibitor in the aqueous phase', type: 'select', value: 'case', options: [{ value: 'case', label: 'As the case fluid' }, { value: 'none', label: 'None' }, { value: 'MeOH', label: 'Methanol' }, { value: 'MEG', label: 'MEG' }, { value: 'DEG', label: 'DEG' }, { value: 'TEG', label: 'TEG' }, { value: 'EtOH', label: 'Ethanol' }] },
    { key: 'inhWt', label: 'Inhibitor concentration', unit: 'wt % of aqueous phase', value: 0, min: 0, max: 90, showIf: (v) => v.inhibitor !== 'case' && v.inhibitor !== 'none' },
    { key: 'inhEff', label: 'Inhibitor effectiveness', unit: '% of ideal depression', value: 100, min: 10, max: 150, help: 'Scales the hydrate-temperature depression of the inhibitor (calibrate against rocking-cell or autoclave data; below 100 % for poor mixing or lean inhibitor).' },
    { key: 'marginC', label: 'Required margin outside the hydrate region', unit: '°C', value: 3, min: 0, max: 15 },
    { key: 'meohVapK', label: 'Methanol loss to gas', unit: 'kg per million Sm³ per wt %', value: 16, min: 0, max: 200, help: 'Screening partition coefficient at cold high-pressure conditions.' },
    { key: 'meohOilK', label: 'Methanol loss to hydrocarbon liquid', unit: 'kg/kg oil per mass fraction', value: 0.004, min: 0, max: 0.1 },
    { key: 'khiLimit', label: 'Kinetic inhibitor subcooling limit (screening)', unit: '°C', value: 10, min: 3, max: 20 },
    { key: 'aaWcLimit', label: 'Anti-agglomerant water-cut limit (screening)', unit: '%', value: 50, min: 10, max: 90 },
  ] },
  { group: 'Hydrate nucleation and growth', tab: 'setup', help: 'Kinetic constants are system-specific: fit them on the Calibration tab to flow-loop, rocking-cell or autoclave data before using the results for design.', fields: [
    { key: 'regime', label: 'Hydrate system', type: 'select', value: 'auto', options: [{ value: 'auto', label: 'Automatic from the flow picture' }, { value: 'oil', label: 'Oil-dominated (water-in-oil emulsion, shrinking core)' }, { value: 'gas', label: 'Gas-dominated (wall film and entrained water)' }, { value: 'water', label: 'Water-dominated (gas absorption limited)' }] },
    { key: 'nucleation', label: 'Nucleation', type: 'select', value: 'heterogeneous', options: [{ value: 'heterogeneous', label: 'Heterogeneous (contact-angle factor)' }, { value: 'homogeneous', label: 'Homogeneous' }] },
    { key: 'nucA', label: 'Nucleation pre-exponential, log₁₀', unit: 'log₁₀(1/m³/s)', value: 7.5, min: 0, max: 40, help: 'About 7–8 for heterogeneous nucleation with the default angle; about 35 for homogeneous nucleation.' },
    { key: 'contactAngle', label: 'Contact angle of the nucleus on the substrate', unit: '°', value: 40, min: 5, max: 180 },
    { key: 'sigmaHW', label: 'Hydrate–water interfacial energy', unit: 'mJ/m²', value: 20, min: 5, max: 60 },
    { key: 'nucVolume', label: 'Water sample volume for the induction time', unit: 'L', value: 1, min: 0.001, max: 1000, help: 'The induction time is the mean waiting time for the first nucleus in this much water, as measured in an autoclave.' },
    { key: 'kinK', label: 'Intrinsic rate constant at 4 °C', unit: '10⁻¹⁰ mol/m²/Pa/s', value: 1, min: 0, max: 1e4, typical: [0.05, 20], help: 'Kim–Bishnoi / Englezos constant on the fugacity difference. 0 switches kinetics off.' },
    { key: 'kinEa', label: 'Activation temperature E/R', unit: 'K', value: 13600, min: 0, max: 30000 },
    { key: 'shellD', label: 'Diffusivity through the hydrate shell', unit: '10⁻¹³ m²/s', value: 5, min: 0.001, max: 1e5 },
    { key: 'mtMult', label: 'Gas-absorption (mass-transfer) multiplier', unit: '–', value: 1, min: 0.001, max: 1000 },
    { key: 'htMult', label: 'Heat-removal multiplier', unit: '–', value: 1, min: 0.05, max: 20, help: 'Scales the heat that can be removed from a forming slurry before it reaches the hydrate temperature.' },
    { key: 'hydNumber', label: 'Hydration number', unit: 'mol water/mol gas', value: 0, min: 0, max: 9, help: '0 = from the Langmuir cage occupancy of methane structure I at line conditions (van der Waals–Platteeuw).' },
    { key: 'rhoHyd', label: 'Hydrate particle density', unit: 'kg/m³', value: 920, min: 800, max: 1100 },
  ] },
  { group: 'Particles, agglomeration and slurry', tab: 'setup', fields: [
    { key: 'primaryUm', label: 'Primary particle size', unit: 'µm', value: 0, min: 0, max: 2000, help: '0 = droplet size from the Boxall inertial correlation, d/D = 0.063 We^−3/5.' },
    { key: 'cohesion', label: 'Cohesive force per unit particle size', unit: 'mN/m', value: 2, min: 0, max: 200, typical: [0.1, 50], help: 'Micromechanical force between hydrate particles divided by their size. Anti-agglomerants reduce it by one to two orders of magnitude.' },
    { key: 'aggEff', label: 'Collision (sticking) efficiency', unit: '–', value: 0.05, min: 0, max: 1 },
    { key: 'kBreak', label: 'Breakage coefficient', unit: '–', value: 0.04, min: 0, max: 10, help: 'Breakage frequency = coefficient × shear rate × (size / cohesive-limit size)³.' },
    { key: 'fractal', label: 'Fractal dimension of agglomerates', unit: '–', value: 2.5, min: 1.8, max: 2.95 },
    { key: 'phiMax', label: 'Maximum packing fraction', unit: '–', value: 0.571, min: 0.3, max: 0.74 },
    { key: 'viscModel', label: 'Slurry viscosity model', type: 'select', value: 'mills', options: [{ value: 'mills', label: 'Mills' }, { value: 'krieger', label: 'Krieger–Dougherty' }, { value: 'thomas', label: 'Thomas' }] },
  ] },
  { group: 'Wall deposition and plugging', tab: 'setup', fields: [
    { key: 'adhesion', label: 'Wall-capture efficiency', unit: '–', value: 0.02, min: 0, max: 1, help: 'Probability that a particle reaching a sub-cooled wall stays there, before the force-balance reduction.' },
    { key: 'adhForce', label: 'Wall adhesion force per unit particle size', unit: 'mN/m', value: 5, min: 0.001, max: 500 },
    { key: 'filmMult', label: 'Vapour-film deposition multiplier', unit: '–', value: 1, min: 0, max: 100, help: 'Scales hydrate growth from water vapour condensing on the gas-wetted cold wall.' },
    { key: 'tauCrit', label: 'Critical wall shear for removal (fresh deposit)', unit: 'Pa', value: 30, min: 0.1, max: 5000, help: 'Above this shear the deposit is eroded and sloughs back into the stream. It rises as the deposit ages and densifies.' },
    { key: 'kRemove', label: 'Removal / detachment rate constant', unit: '1/h', value: 0.5, min: 0, max: 100 },
    { key: 'porosity0', label: 'Porosity of a fresh deposit', unit: '–', value: 0.6, min: 0.05, max: 0.95 },
    { key: 'porosityInf', label: 'Porosity of an aged deposit', unit: '–', value: 0.2, min: 0.01, max: 0.9 },
    { key: 'ageHours', label: 'Ageing time constant', unit: 'h', value: 12, min: 0.1, max: 2000 },
    { key: 'plugGrainUm', label: 'Grain size for plug permeability', unit: 'µm', value: 0, min: 0, max: 5000, help: '0 = the primary particle size. Used in the Kozeny–Carman permeability.' },
    { key: 'plugBlockPct', label: 'Area blockage that counts as a plug', unit: '%', value: 90, min: 30, max: 99.5 },
    { key: 'plugLength', label: 'Plug length for remediation estimates', unit: 'm', value: 50, min: 1, max: 5000 },
    { key: 'depressP', label: 'Depressurisation pressure for plug melting', unit: 'bara', value: 5, min: 1, max: 100 },
  ] },
  { group: 'Wax', tab: 'setup', fields: [
    { key: 'wat', label: 'Wax appearance temperature', unit: '°C', value: 45, min: -20, max: 90 },
    { key: 'waxContent', label: 'Wax content of the oil', unit: 'wt %', value: 5, min: 0, max: 40 },
    { key: 'waxSlope', label: 'Solubility slope below the WAT', unit: '1/K', value: 0.04, min: 0.005, max: 0.2, help: 'Dissolved wax = total × exp(−slope × (WAT − T)).' },
    { key: 'waxDiff', label: 'Diffusivity correlation', type: 'select', value: 'haydukMinhas', options: [{ value: 'haydukMinhas', label: 'Hayduk–Minhas' }, { value: 'wilkeChang', label: 'Wilke–Chang' }] },
    { key: 'waxMult', label: 'Wax deposition multiplier', unit: '–', value: 1, min: 0, max: 100, help: 'Multiplies the diffusive flux (fit to flow-loop or cold-finger data).' },
    { key: 'waxOil', label: 'Oil trapped in a fresh deposit', unit: 'vol %', value: 80, min: 5, max: 97 },
    { key: 'waxK', label: 'Deposit thermal conductivity', unit: 'W/m/K', value: 0.25, min: 0.05, max: 1 },
    { key: 'waxLimitMm', label: 'Thickness that triggers pigging', unit: 'mm', value: 2, min: 0.1, max: 50 },
    { key: 'pourOffset', label: 'Gelling starts this far below the WAT', unit: '°C', value: 15, min: 0, max: 60 },
    { key: 'gelCoef', label: 'Gel yield-stress coefficient', unit: 'Pa per (wt % solid wax)²', value: 4, min: 0, max: 500 },
  ] },
  { group: 'Scale', tab: 'setup', fields: [
    { key: 'water', label: 'Produced-water analysis', type: 'table', columns: ION_IDS.map((id) => ({ key: id, label: id, unit: 'mg/L' })), value: [{ ...WATER0 }], help: 'First row is used. Bicarbonate is the alkalinity.' },
    { key: 'co2Pct', label: 'CO₂ in the gas phase', unit: 'mol %', value: 3, min: 0.001, max: 60 },
    { key: 'swFrac', label: 'Seawater fraction in the produced water', unit: '%', value: 0, min: 0, max: 100, help: 'Injection-water breakthrough: produced water is mixed with standard seawater before the indices are evaluated.' },
    { key: 'actModel', label: 'Activity-coefficient model', type: 'select', value: 'truesdellJones', options: [{ value: 'truesdellJones', label: 'Truesdell–Jones extended Debye–Hückel' }, { value: 'davies', label: 'Davies' }] },
    { key: 'scaleK', label: 'Scale surface-reaction constant at 25 °C', unit: '10⁻⁸ mol/m²/s', value: 1, min: 0, max: 1e4 },
  ] },
  { group: 'Asphaltene', tab: 'setup', fields: [
    { key: 'saraSat', label: 'Saturates', unit: 'wt %', value: 55, min: 0, max: 100 }, { key: 'saraAro', label: 'Aromatics', unit: 'wt %', value: 28, min: 0, max: 100 }, { key: 'saraRes', label: 'Resins', unit: 'wt %', value: 14, min: 0, max: 100 }, { key: 'saraAsp', label: 'Asphaltenes', unit: 'wt %', value: 3, min: 0, max: 100 },
    { key: 'asphDelta', label: 'Asphaltene solubility parameter at 25 °C', unit: 'MPa^½', value: 20, min: 17, max: 24 }, { key: 'asphMV', label: 'Asphaltene molar volume', unit: 'm³/kmol', value: 2, min: 0.3, max: 10 },
  ] },
  { group: 'Sand and solids', tab: 'setup', fields: [
    { key: 'sandRate', label: 'Sand production', unit: 'kg/d', value: 50, min: 0, max: 1e5 }, { key: 'sandUm', label: 'Sand particle size', unit: 'µm', value: 150, min: 10, max: 3000 }, { key: 'sandRho', label: 'Sand density', unit: 'kg/m³', value: 2650, min: 1100, max: 5000 },
  ] },
  { group: 'Discretisation', tab: 'mesh', help: 'Axial cells, time step of the hydrate march and number of particle-size classes; the Monte Carlo runs use a coarser copy of the same model.', fields: [
    { key: 'nAxial', label: 'Axial cells', value: 40, min: 8, max: 120 },
    { key: 'dtMin', label: 'Time step of the hydrate march', unit: 'min', value: 15, min: 0.5, max: 240 },
    { key: 'nClasses', label: 'Particle-size classes', value: 14, min: 6, max: 40 },
    { key: 'nMC', label: 'Monte Carlo samples for plugging probability', value: 32, min: 0, max: 400 },
    { key: 'seed', label: 'Random seed', value: 7, min: 1, max: 1e6 },
  ] },
];
export const INPUT_FIELDS = INPUTS.flatMap((g) => g.fields);
const FIELDS = INPUTS.flatMap((g) => g.fields), DEF = Object.fromEntries(FIELDS.map((f) => [f.key, f.value]));
/** Fill defaults, coerce and check the inputs; throws a readable Error for impossible input. */
function clean(v0 = {}) {
  const v = {};
  for (const f of FIELDS) {
    const x = v0[f.key];
    if (f.type === 'select') { v[f.key] = f.options.some((o) => o.value === x) ? x : f.value; continue; }
    if (f.type === 'table') { v[f.key] = Array.isArray(x) && x.length ? x : f.value; continue; }
    const y = x === null || x === undefined || x === '' ? f.value : +x;
    need(Number.isFinite(y), `${f.label} must be a number.`);
    need(f.min === undefined || y >= f.min - 1e-12, `${f.label} must be at least ${f.min}${f.unit ? ' ' + f.unit : ''} (got ${y}).`);
    need(f.max === undefined || y <= f.max + 1e-12, `${f.label} must not exceed ${f.max}${f.unit ? ' ' + f.unit : ''} (got ${y}).`);
    v[f.key] = y;
  }
  for (const k of ['nAxial', 'nClasses', 'nMC', 'seed']) v[k] = Math.round(v[k]);
  need(v.porosityInf <= v.porosity0, 'The aged deposit porosity cannot exceed the fresh deposit porosity.');
  need(v.saraSat + v.saraAro + v.saraRes + v.saraAsp > 0, 'The SARA analysis is empty: enter the saturate, aromatic, resin and asphaltene fractions.');
  need(v.dep0Mm < 0.45 * v.idMm, 'The initial hydrate deposit is thicker than the pipe can hold.');
  need(v.pInMax > 1.5, 'The available inlet pressure must exceed the arrival pressure.');
  return v;
}
/** Phases of the hydrate march for a scenario: [{ dur (s), dt (s), frac }] and the rate fraction of the steady production period. */
function scenarioPhases(v, dtScale = 1) {
  const dt = v.dtMin * 60 * dtScale, sim = v.simHours * 3600, sh = v.shutHours * 3600, dtS = Math.min(Math.max(4 * dt, sh / 60), sh);
  if (v.scenario === 'steady') return { phases: [{ dur: sim, dt, frac: 1 }], prod: 1 };
  if (v.scenario === 'turndown') return { phases: [{ dur: sim, dt, frac: v.turndownPct / 100 }], prod: 1 };
  if (v.scenario === 'shutdown') return { phases: [{ dur: sh, dt: dtS, frac: 0 }], prod: 1 };
  return { phases: [{ dur: sh, dt: dtS, frac: 0 }, { dur: sim, dt, frac: v.restartPct / 100 }], prod: 1 };
}

// =====================================================================================================
// 15. The run
// =====================================================================================================
const r3 = (x, s = 3) => (Number.isFinite(x) ? +x.toPrecision(s) : null), tx = (x, s = 3) => (Number.isFinite(x) ? +x.toPrecision(s) : '—');
const argmax = (a) => { let k = 0; for (let i = 1; i < a.length; i++) if (a[i] > a[k]) k = i; return k; };
const SC_LABEL = Object.fromEntries(SCEN.map((s) => [s.value, s.label]));
/** Hydrate-temperature depression (°C) of the aqueous phase of a set-up, with the inhibitor effectiveness applied to the inhibitor part. */
const aqDepression = (S, eff = 1, wt = S.aq.inhWt) => { const salt = hydrateDepression({ S: S.aq.S, inhWt: 0, inh: S.aq.inh }); return wt > 0 ? salt + eff * (hydrateDepression({ S: S.aq.S, inhWt: wt, inh: S.aq.inh }) - salt) : salt; };

/** Monte Carlo over uncertain kinetics, cohesion, adhesion and the stochastic (Poisson) onset on a coarse copy of the model. */
async function plugMonteCarlo(v, ctx, S, p, slow, tick) {
  const N = v.nMC, out = { n: N, plugTimes: [], onsetTimes: [], blk: [], visc: [], prob: 0, lo: 0, hi: 0 };
  if (!(N > 0)) return out;
  const n2 = Math.min(S.n, 12), S2 = buildSetup(v, ctx, n2), g2 = pbeGrid(8, p.mode === 'oil' ? p.dPrim / 4 : 2e-6, 2e-2), { phases } = scenarioPhases(v, 2), rg = rng(v.seed);
  const dW = S2.s.map((s) => interp1(S.s, slow.dWax, s) + interp1(S.s, slow.dScale, s)), Dbase = dW.map((d) => Math.max(S.D0 - 2 * d, 0.2 * S.D0)), ln = (sd) => Math.exp(rg.normal(0, sd));
  for (let k = 0; k < N; k++) {
    const pk = { ...p, kinK: p.kinK * ln(0.7), cohesion: p.cohesion * ln(0.6), adhesion: Math.min(p.adhesion * ln(0.6), 1), adhForce: p.adhForce * ln(0.5), tauCrit: p.tauCrit * ln(0.4), nucA: p.nucA * ln(1.5), shellD: p.shellD * ln(0.7), lamStar: -Math.log(1 - rg.uniform(0, 1) * 0.999999) };
    const r = runHydrateMarch(S2, pk, { phases, grid: g2, Dbase, dep0: v.dep0Mm / 1000, cheap: true, rows: 2 });
    if (r.plug) out.plugTimes.push(r.plug.t / 3600); if (r.onset) out.onsetTimes.push(r.onset.t / 3600);
    out.blk.push(Math.max(...r.ser.blk, 0)); out.visc.push(Math.max(...r.ser.visc, 1));
    if (k % 8 === 7) await tick(0.7 + (0.2 * k) / N);
  }
  const ph = out.plugTimes.length / N, z = 1.645, den = 1 + (z * z) / N, c = (ph + (z * z) / (2 * N)) / den, hw = (z * Math.sqrt((ph * (1 - ph)) / N + (z * z) / (4 * N * N))) / den; // Wilson 90 % interval
  out.prob = ph; out.lo = Math.max(c - hw, 0); out.hi = Math.min(c + hw, 1);
  return out;
}

async function run(v0, ctx = {}) {
  const v = clean(v0), prog = (f, m) => { try { ctx.progress?.(f, m); } catch { /* progress is optional */ } }, tick = async (f, m) => { if (f !== undefined) prog(f, m || 'Working'); if (ctx.tick) await ctx.tick(); };
  prog(0.02, 'Building the flow picture');
  const S = buildSetup(v, ctx, clamp(v.nAxial, 8, 120)), n = S.n, D0 = S.D0, ds = S.ds, xkm = S.x.map((a) => a / 1000), A0 = (PI * D0 * D0) / 4;
  const { phases } = scenarioPhases(v), fProd = v.scenario === 'turndown' ? v.turndownPct / 100 : 1, fEvent = phases[phases.length - 1].frac || 1;
  for (const ph of phases) if (ph.frac > 0) S.prof(ph.frac);
  const base = S.prof(1), pb = S.prof(fProd);
  await tick(0.08, 'Wax and scale in steady production');
  const slow = slowDeposits(S, v, fProd);
  const Dbase = slow.dWax.map((d, i) => Math.max(D0 - 2 * (d + slow.dScale[i]), 0.2 * D0)), roughBase = slow.dWax.map((d, i) => Math.min(S.rough0 + 0.05 * d + 0.3 * slow.dScale[i], 0.05 * D0));
  const p = hydrateParams(v, S), oil = p.mode === 'oil', grid = pbeGrid(v.nClasses, oil ? p.dPrim / 4 : 2e-6, 2e-2);
  // ---- hydrate event marched in time
  const m = hydrateMarch(S, p, { phases, grid, Dbase, roughBase, dep0: v.dep0Mm / 1000, rows: 48 });
  let kk = 0; while (!m.done) { m.step(); if (++kk % 10 === 0) await tick(0.1 + 0.55 * m.progress, 'Marching hydrate formation, transport and deposition'); }
  const h = m.result(), rec = h.rec, ser = h.ser, pk = h.peak, tShut = v.scenario === 'restart' ? v.shutHours : 0;
  await tick(0.68, 'Plugging probability (Monte Carlo)');
  const mc = await plugMonteCarlo(v, ctx, S, p, slow, tick);
  await tick(0.92, 'Wax, scale, asphaltene, sand and remediation');

  // ---- hydrate driving force, exposure, conversion
  const depNow = aqDepression(S, p.inhEff), depSalt = aqDepression(S, 1, 0), maxSub = Math.max(...rec.subMax), iSub = argmax(rec.subMax), stableLen = rec.subMax.filter((s) => s > 0).length * ds, stableEnd = rec.sub.filter((s) => s > 0).length * ds;
  let expoH = 0, degH = 0; for (let k = 0; k < ser.t.length; k++) { const dtH = ser.t[k] - (k ? ser.t[k - 1] : 0); if (ser.sub[k] > 0) { expoH += dtH; degH += ser.sub[k] * dtH; } }
  const teqSteady = base.P.map((P) => S.hT0(P) - depNow), subSteady = base.T.map((T, i) => teqSteady[i] - T);
  const tauC = S.C.map((C) => C / (S.U * PI * D0)), coolT = base.T.map((T, i) => (T <= teqSteady[i] ? 0 : S.tAmb[i] >= teqSteady[i] ? Infinity : tauC[i] * Math.log((T - S.tAmb[i]) / (teqSteady[i] - S.tAmb[i]))) / 3600), cooldown = Math.min(...coolT), iCool = coolT.indexOf(cooldown);
  const peakPhi = pk.phi, peakBlkHyd = pk.blk, peakVisc = Math.max(...ser.visc, 1), waterIn = ser.t.length ? base.P.reduce((a, P, i) => a + S.props(P, base.T[i]).rhoW * S.props(P, base.T[i]).wcut * A0 * base.holdup[i] * ds, 0) : 0;
  const led = h.ledger, waterUsed = (led.formed - led.dissociated) * p.wfH, gasUsed = led.formed * (1 - p.wfH), peakRate = Math.max(...ser.rate, 0), convMax = Math.max(...rec.X, pk.phiX ? Math.max(...pk.phiX.map((f, i) => (f * p.rhoH * p.wfH) / Math.max(S.props(h.P[i], h.T[i]).rhoW * S.props(h.P[i], h.T[i]).wcut, 1e-9))) : 0);
  // ---- combined deposit profile at the worst time of the event
  const dHyd = pk.dHyd || h.dHyd, total = dHyd.map((d, i) => d + slow.dWax[i] + slow.dScale[i]), Deff = total.map((d) => Math.max(D0 - 2 * d, 0.03 * D0)), iMin = argmax(total), blockage = 1 - (Deff[iMin] / D0) ** 2;
  const roughEff = Math.min(S.rough0 + Math.max(...total.map((_, i) => 0.1 * dHyd[i] + 0.05 * slow.dWax[i] + 0.3 * slow.dScale[i])), 0.05 * D0), grain = v.plugGrainUm > 0 ? v.plugGrainUm * 1e-6 : p.dPrim;
  const porPk = pk.por || h.por, perm = porPk.map((e) => kozenyCarman(e, grain));
  // clean versus fouled pressure drop in production (kernel gradient with the restricted bore and the rougher wall)
  let dpClean = 0, dpFoul = 0; for (let i = 0; i < n; i++) { const g0 = S.grad(i, pb.P[i], pb.T[i], D0, S.rough0, 1, fProd), g1 = S.grad(i, pb.P[i], pb.T[i], Math.max(D0 - 2 * (slow.dWax[i] + slow.dScale[i]), 0.2 * D0), roughBase[i], 1, fProd); dpClean += (g0.dpdx * ds) / 1e5; dpFoul += (g1.dpdx * ds) / 1e5; }
  const pIn0 = ser.pIn.find((_, k) => ser.frac[k] > 0) ?? ser.pIn[0] ?? base.pIn, pInPeak = Math.max(...ser.pIn.filter((_, k) => ser.frac[k] > 0), pIn0), dpRise = pInPeak - pIn0;
  // ---- plugging indicator
  const plug = h.plug, pOn = h.onset ? 1 : 1 - Math.exp(-Math.max(...rec.lam, 0)), sev = Math.max(peakBlkHyd / p.plugBlock, Math.log10(peakVisc) / 3, dpRise / Math.max(p.pInMax - pIn0, 1));
  const risk = plug ? 1 : clamp(Math.max(pOn * Math.max(maxSub > 0 ? 0.15 : 0, sev), mc.prob), 0, 1), plugT = mc.plugTimes.slice().sort((a, b) => a - b);
  // ---- inhibition
  const inhId = S.aq.inhId !== 'none' ? S.aq.inhId : 'MEG', inh = INHIBITORS[inhId], needDep = maxSub + (depNow - depSalt) + v.marginC, wReqIdeal = needDep > 0 ? inhibitorFor(depSalt + needDep / Math.max(p.inhEff, 0.1), inhId, S.aq.S) : 0;
  const prE = S.props(base.P[n - 1], base.T[n - 1]), mWater = prE.mW * fEvent, wR = Math.min(wReqIdeal, 94) / 100, mInh = (wR / (1 - wR)) * mWater, inhRate = (mInh / inh.rho) * 86400, qGasStd = (S.fm?.rates?.qGasStd ?? 0) * fEvent;
  const lossVap = inhId === 'MeOH' ? (v.meohVapK * wReqIdeal * qGasStd) / 1e6 : 0, lossOil = inhId === 'MeOH' ? v.meohOilK * wR * prE.mO * fEvent * 86400 : 0, lossPct = mInh > 0 ? (100 * (lossVap + lossOil)) / (mInh * 86400) : 0;
  const khiOk = maxSub <= v.khiLimit, khiHold = maxSub > 0 ? 48 * 2 ** ((v.khiLimit - maxSub) / 1.5) : Infinity, aaOk = prE.wcut * 100 <= v.aaWcLimit && prE.qO > 0 && p.mode !== 'gas';
  // ---- remediation at the plug (or at the coldest point when there is no plug)
  const iR = plug ? plug.i : iSub, Rb = D0 / 2, pLoc = h.P[iR], tA = S.tAmb[iR], tdTwo = S.hT0(v.depressP) - depNow, melt2 = plugMeltTime({ R: Rb, U: S.U, Tamb: tA, Td: tdTwo, eps: v.porosityInf }), tdUp = S.hT0(pLoc) - depNow;
  const alphaH = HYDRATE.k / (p.rhoH * HYDRATE.cp), ste = (HYDRATE.cp * Math.max(tA - tdTwo, 0)) / HYDRATE.latent, lamS = ste > 0 ? stefanLambda(ste) : 0, melt1 = lamS > 0 ? v.plugLength ** 2 / (4 * lamS * lamS * alphaH) : Infinity;
  const projV = Math.sqrt((2 * Math.max(pLoc - v.depressP, 0) * 1e5 * 100) / (p.rhoH * (1 - v.porosityInf) * v.plugLength)), mPlug = p.rhoH * (1 - v.porosityInf) * A0 * v.plugLength;
  const wEq = inhibitorFor(depSalt + Math.max(tdUp + depNow - depSalt - tA, 0) + 1, 'MeOH', S.aq.S) / 100, meohMelt = (mPlug * p.wfH * wEq) / Math.max(1 - wEq, 0.05) / INHIBITORS.MeOH.rho, tMeoh = meohMelt > 0 ? (meohMelt * INHIBITORS.MeOH.rho) / (INHIBITORS.MeOH.rho * A0 * 1e-5 * Math.max(1 - wEq, 0.05)) : 0;
  const kdis = hydrateDissociationRate({ TK: tA + KEL, P: v.depressP, Peq: S.peq(tA + depNow), zG: 0.95 }), tKin = kdis > 0 ? (p.rhoH * grain) / (6 * (p.mwG + p.hydN * MW_W) * kdis) : Infinity;
  const tHold = Math.max(...teqSteady) + v.marginC, dehW = base.T.reduce((a, _, i) => a + Math.max(S.U * PI * D0 * (tHold - S.tAmb[i]) * ds, 0), 0), meltKW = (mPlug * HYDRATE.latent) / (24 * 3600) / 1000;
  // ---- sand
  const sandQ = v.sandRate / 86400 / v.sandRho, sand = pb.P.map((P, i) => { const c = slow.cells[i], vsl = c.vL * pb.holdup[i], C = sandQ / Math.max(c.pr.qL * fProd, 1e-12), cr = sandCriticalVelocity({ d: v.sandUm * 1e-6, D: D0, rhoP: v.sandRho, rhoF: c.pr.rhoL, mu: c.pr.muL, C, vsl, vm: c.vm }), flat = Math.abs(S.theta[i]) < 0.5, st = settlingVelocity(v.sandUm * 1e-6, v.sandRho, c.pr.rhoL, c.pr.muL, { phi: C });
    return { C, cr, vm: c.vm, vL: c.vL, flat, bed: v.sandRate > 0 && (flat ? c.vm < cr.governing : c.vL < 3 * Math.abs(st.v)), st, hold: (C * c.vL) / Math.max(c.vL - Math.abs(st.vHindered), 0.05 * c.vL) }; });
  const flatS = sand.filter((s) => s.flat), sandCrit = flatS.length ? Math.max(...flatS.map((s) => s.cr.governing)) : Math.max(...sand.map((s) => s.cr.governing)), sandBed = sand.some((s) => s.bed), bedLen = sand.filter((s) => s.bed).length * ds;
  const sandMargin = Math.min(...sand.map((s) => (s.flat ? s.vm / s.cr.governing : 9))), sandMinRate = sandMargin > 0 ? (100 * fProd) / sandMargin : 100, sandInv = sand.reduce((a, s, i) => a + s.hold * v.sandRho * A0 * pb.holdup[i] * ds, 0), iV = argmax(sand.map((s) => s.vm)), erosion = sandErosionScreen(v.sandRate, sand[iV].vm, v.sandUm, D0 * 1000, pb.rhoM[iV]);
  const relax = particleRelaxation({ d: v.sandUm * 1e-6, rhoP: v.sandRho, rhoF: slow.cells[0].pr.rhoL, mu: slow.cells[0].pr.muL });
  // ---- asphaltene
  const spec = S.fm?.spec || {}, pRes = num(spec.Pres, 300), tRes = num(spec.Tres, 90), tArr = base.T[n - 1], aPath = linspace(pRes, Math.max(S.pOut, 2), 40).map((P, k) => { const T = tRes + ((tArr - tRes) * k) / 39, o = S.fm ? S.fm.at(P, T) : S.props(P, T); return { P, T, wG: o.wG ?? 0, rho: o.rhoO, mw: o.mwO || 150 }; });
  let pBub = aPath[aPath.length - 1].P; for (const q of aPath) if (q.wG > 1e-4) { pBub = q.P; break; }
  const sara = { sat: v.saraSat, aro: v.saraAro, res: v.saraRes, asp: v.saraAsp }, cii = colloidalInstability(sara), db = deBoer(aPath[0].rho, Math.max(pRes - pBub, 0)), aspPhi = ((v.saraAsp / 100) * aPath[0].rho) / 1200;
  const fh = aPath.map((q) => asphalteneSolubility({ rhoL: q.rho, mwL: q.mw, TK: q.T + KEL, deltaA: v.asphDelta, vA: v.asphMV })), iFh = fh.findIndex((f) => f.phiMax < aspPhi), fhMin = Math.min(...fh.map((f) => f.phiMax));
  const aScore = (db.cls === 'severe problems' ? 2 : db.cls === 'slight problems' ? 1 : 0) + (cii.cls === 'unstable' ? 2 : cii.cls === 'uncertain' ? 1 : 0) + (iFh >= 0 ? 2 : 0), aRisk = aScore >= 4 ? 'high' : aScore >= 2 ? 'medium' : 'low';
  // ---- scale summary and seawater mixing curve at arrival conditions
  const iSI = argmax(slow.si), scaleSI = slow.si[iSI], scaleMineral = scaleSI > 0 ? slow.mineral[iSI] : 'none', mixF = linspace(0, 1, 21), w0 = waterOf(v), tMix = pb.T[0], pMix = pb.P[0];
  const mix = mixF.map((f) => scaleIndices(mixWaters(w0, SEAWATER, f), tMix, pMix, { yCO2: v.co2Pct / 100, model: v.actModel })), sIn = slow.sis[0], mic = scaleSI > 0 ? clamp(2 * 10 ** (0.5 * scaleSI), 1, 200) : 0;
  // ---- sectional versus quadrature method of moments in the worst cell (short batch: aggregation + breakage)
  const qm = (() => {
    const phi = Math.max(peakPhi, 0.02), env = { shear: pk.gdot || 200, eps: 0, nu: (pk.muC || 3e-3) / (pk.rhoC || 800), mu: pk.muC || 3e-3, TK: (pk.T ?? 4) + KEL, dRho: Math.abs(p.rhoH - (pk.rhoC || 800)), alpha: p.aggEff }, g2 = pbeGrid(40, p.dPrim / 2, 60 * p.dPrim), kP = 5, L0 = g2.L[kP], N0 = phi / g2.v[kP];
    const b0 = aggregationKernel(L0, L0, env).total, tB = b0 * N0 > 0 ? 3 / (b0 * N0) : 1, be = new Float64Array(g2.n * g2.n); for (let i = 0; i < g2.n; i++) for (let j = 0; j < g2.n; j++) be[i * g2.n + j] = aggregationKernel(g2.L[i], g2.L[j], env).total;
    const Ni = new Float64Array(g2.n); Ni[kP] = N0; const sct = solvePBE(g2, Ni, tB, { beta: be, maxSub: 4000, frac: 0.05 }), ms = pbeMoments(g2, sct.N);
    const sg = 0.05, mom = [0, 1, 2, 3, 4, 5].map((k) => N0 * L0 ** k * Math.exp((k * k * sg * sg) / 2)), q = solveQMOM(mom, tB, { beta: (a, b) => aggregationKernel(a, b, env).total, steps: 120 });
    return { tB, m0s: ms.m0 / N0, m0q: q.m[0] / N0, m3s: ms.m3 / (N0 * L0 ** 3), m3q: q.m[3] / mom[3], d43s: ms.d43, d43q: q.m[4] / q.m[3], ok: q.ok, nodes: q.nodes };
  })();
  const vdw = vdwpMethane(Math.max(mean(S.tAmb), 0.5)), occ = langmuirOccupancy(Math.max(mean(S.tAmb), 0) + KEL, methaneFugacity(mean(base.P), Math.max(mean(S.tAmb), 0) + KEL));

  // ---- threats ranked by zone
  const nz = Math.min(8, n), zones = Array.from({ length: nz }, (_, zI) => {
    const a = Math.floor((zI * n) / nz), b = Math.max(Math.floor(((zI + 1) * n) / nz), a + 1), ix = Array.from({ length: b - a }, (_, k) => a + k), mx = (f) => Math.max(...ix.map(f));
    const sc = { Hydrate: clamp(Math.max(mx((i) => rec.subMax[i]) / 10, mx((i) => dHyd[i] / (0.15 * D0)), mx((i) => (pk.phiX ? pk.phiX[i] : 0) / 0.2)), 0, 1.5) * (mx((i) => rec.subMax[i]) > 0 ? 1 : 0), Wax: clamp(mx((i) => slow.dWax[i] * 1000) / v.waxLimitMm, 0, 1.5), Scale: clamp(Math.max(mx((i) => slow.si[i]) / 1.5, mx((i) => slow.scaleRate[i]) / 2), 0, 1.5), Sand: clamp(mx((i) => (v.sandRate > 0 ? sand[i].cr.governing / Math.max(sand[i].vm, 1e-6) : 0)) - 0.3, 0, 1.5) };
    const top = Object.entries(sc).sort((p1, p2) => p2[1] - p1[1])[0];
    return { from: (a * ds) / 1000, to: (b * ds) / 1000, sc, top: top[1] > 0.05 ? top[0] : 'None' };
  });

  // ---- outputs for the other suites
  const tOnset = h.onset ? h.onset.t / 3600 : null, piggingInterval = slow.piggingInterval;
  const outputs = {
    hydrateRisk: r3(risk), maxSubcooling: r3(maxSub, 4), onsetX: h.onset ? r3(h.onset.x, 5) : null, onsetTime: tOnset === null ? null : r3(tOnset, 4), hydrateFraction: r3(peakPhi, 4), hydrateRate: r3(peakRate, 4),
    depositProfile: { x: S.x.map((a) => +a.toFixed(1)), hydrate: dHyd.map((d) => +d.toExponential(4)), wax: slow.dWax.map((d) => +d.toExponential(4)), scale: slow.dScale.map((d) => +d.toExponential(4)), total: total.map((d) => +d.toExponential(4)) },
    effectiveId: r3(Deff[iMin], 5), roughnessEff: r3(roughEff, 4), blockage: r3(blockage, 4), plugTime: plug ? r3(plug.t / 3600, 4) : null, plugX: plug ? r3(plug.x, 5) : null, plugProbability: r3(plug && mc.n === 0 ? 1 : mc.prob),
    inhibitorRequired: r3(wReqIdeal, 4), inhibitorRate: r3(inhRate, 4), wat: v.wat, waxRate: r3(slow.waxRate0, 4), waxMass: r3(slow.waxMass, 4), piggingInterval: piggingInterval === null ? null : r3(piggingInterval, 4),
    scaleSI: r3(scaleSI, 4), scaleMineral, asphalteneRisk: aRisk, sandCriticalVelocity: r3(sandCrit, 4), sandBed, slurryViscosityFactor: r3(peakVisc, 4),
    // extras
    scenario: v.scenario, hydrateMode: p.mode, hydrateLength: r3(stableLen, 5), exposureHours: r3(expoH, 4), cooldownTime: Number.isFinite(cooldown) ? r3(cooldown, 4) : null, hydrateMass: r3(Math.max(...ser.susp.map((s, k) => s + ser.dep[k]), 0), 4), waterConversion: r3(convMax, 4), gasConsumed: r3(gasUsed, 4),
    plugProbabilityInterval: [r3(mc.lo), r3(mc.hi)], plugTimeP10: plugT.length ? r3(quantile(plugT, 0.1), 4) : null, plugTimeP50: plugT.length ? r3(quantile(plugT, 0.5), 4) : null, plugTimeP90: plugT.length ? r3(quantile(plugT, 0.9), 4) : null, plugMechanism: plug ? plug.mech : null,
    dpIncrease: r3(dpRise + (dpFoul - dpClean), 4), dpFouling: r3(dpFoul - dpClean, 4), inhibitor: inhId, inhibitorDepression: r3(depNow - depSalt, 4), meltTimeTwoSided: Number.isFinite(melt2.tNumeric) ? r3(melt2.tNumeric / 3600, 4) : null, heatingPower: r3(dehW / 1000, 4),
    waxRestartPressure: r3(slow.restartDp, 4), sandRate: r3(v.sandRate / 86400, 4), sandSize: v.sandUm * 1e-6, sandDensity: v.sandRho, sandMinRateFraction: r3(sandMinRate / 100, 4), sandErosionScreen: r3(erosion, 3), asphalteneCII: r3(cii.cii, 4), particleSize: r3(pk.d43 ? pk.d43[pk.iPhi] : 0, 4),
  };

  // ---- KPIs
  const st = (x, w, b) => (x >= b ? 'bad' : x >= w ? 'warn' : 'ok');
  const kpis = [
    { label: 'Peak subcooling in the event', value: r3(maxSub, 3), unit: '°C', status: st(maxSub, 0.01, 6), help: `Hydrate temperature minus fluid temperature, largest value along the line and through the ${SC_LABEL[v.scenario].toLowerCase()} (positive = inside the hydrate region).` },
    { label: 'Length inside the hydrate region', value: r3(stableLen / 1000, 3), unit: 'km', status: st(stableLen, 1, 0.25 * S.L), help: 'Pipe length that was inside the hydrate region at any time of the event.' },
    { label: 'Hydrate onset', value: tOnset === null ? 'none' : r3(tOnset, 3), unit: tOnset === null ? '' : 'h', status: tOnset === null ? 'ok' : 'bad', help: 'Time from the start of the simulated sequence to the first nucleation (hazard integral reaches one).' },
    { label: 'Peak hydrate fraction of the liquid', value: r3(peakPhi * 100, 3), unit: 'vol %', status: st(peakPhi, 0.01, 0.1) },
    { label: 'Peak slurry viscosity factor', value: r3(peakVisc, 3), unit: '×', status: st(peakVisc, 2, 20), help: 'Relative viscosity of the hydrate slurry from the effective (agglomerate) volume fraction.' },
    { label: 'Peak area blockage (all deposits)', value: r3(blockage * 100, 3), unit: '%', status: st(blockage, 0.1, 0.5) },
    { label: 'Plug', value: plug ? r3(plug.t / 3600, 3) : 'no plug', unit: plug ? `h at ${(plug.x / 1000).toFixed(1)} km` : '', status: plug ? 'bad' : 'ok', help: plug ? plug.mech : 'No plug in the base-case march; see the plugging probability for the effect of uncertain kinetics.' },
    { label: 'Plugging probability', value: r3(outputs.plugProbability * 100, 3), unit: '%', status: st(outputs.plugProbability, 0.05, 0.3), help: `Fraction of ${mc.n} Monte Carlo samples that plug within the event (90 % interval ${(mc.lo * 100).toFixed(0)}–${(mc.hi * 100).toFixed(0)} %).` },
    { label: 'Hydrate risk index', value: r3(risk, 3), unit: '0–1', status: st(risk, 0.2, 0.6) },
    { label: `${inh.name} needed for a ${v.marginC} °C margin`, value: r3(wReqIdeal, 3), unit: 'wt %', status: wReqIdeal > S.aq.inhWt + 0.5 ? (wReqIdeal > 60 ? 'bad' : 'warn') : 'ok', help: `Injection rate ${tx(inhRate)} m³/d at the event rate.` },
    { label: 'Cooldown time to hydrate temperature', value: Number.isFinite(cooldown) ? r3(cooldown, 3) : 'never', unit: Number.isFinite(cooldown) ? 'h' : '', status: Number.isFinite(cooldown) ? st(-cooldown, -12, -4) : 'ok', help: 'Lumped exponential cooldown of the first station to reach its hydrate temperature after a shut-in from steady production.' },
    { label: 'Initial wax build-up rate', value: r3(slow.waxRate0, 3), unit: 'mm/d', status: st(slow.waxRate0, 0.02, 0.2) },
    { label: 'Pigging interval', value: piggingInterval === null ? 'not needed' : r3(piggingInterval, 3), unit: piggingInterval === null ? '' : 'd', status: piggingInterval === null ? 'ok' : st(-piggingInterval, -60, -10) },
    { label: 'Highest scale saturation index', value: r3(scaleSI, 3), unit: scaleSI > 0 ? scaleMineral : '', status: st(scaleSI, 0.01, 1) },
    { label: 'Asphaltene risk', value: aRisk, unit: '', status: aRisk === 'high' ? 'bad' : aRisk === 'medium' ? 'warn' : 'ok' },
    { label: 'Sand transport margin', value: r3(sandMargin, 3), unit: 'v / v_critical', status: v.sandRate > 0 ? st(-sandMargin, -1.3, -1) : 'ok' },
  ];

  // ---- warnings and recommendations
  const warnings = [], recs = [];
  if (plug) warnings.push({ level: 'bad', msg: `The base-case march plugs after ${(plug.t / 3600).toFixed(1)} h at ${(plug.x / 1000).toFixed(1)} km: ${plug.mech}.` });
  if (maxSub > 0 && !plug) warnings.push({ level: maxSub > 6 ? 'bad' : 'warn', msg: `The line is inside the hydrate region over ${(stableLen / 1000).toFixed(1)} km with up to ${maxSub.toFixed(1)} °C of subcooling for ${expoH.toFixed(1)} h.` });
  if (mc.prob >= 0.05) warnings.push({ level: mc.prob > 0.3 ? 'bad' : 'warn', msg: `Plugging probability ${(mc.prob * 100).toFixed(0)} % (90 % interval ${(mc.lo * 100).toFixed(0)}–${(mc.hi * 100).toFixed(0)} %, ${mc.n} samples of uncertain kinetics, cohesion, adhesion and nucleation).` });
  if (h.limited) warnings.push({ level: 'info', msg: 'Agglomeration and breakage are much faster than the transport step in part of the line; those cells were integrated with the implicit conservative scheme, which resolves the equilibrium size but not the sub-second transient.' });
  if (piggingInterval !== null && piggingInterval < 30) warnings.push({ level: 'warn', msg: `Wax reaches ${v.waxLimitMm} mm in ${piggingInterval.toFixed(0)} d at ${(slow.waxRate0).toFixed(3)} mm/d.` });
  if (scaleSI > 0.5) warnings.push({ level: scaleSI > 1 ? 'bad' : 'warn', msg: `${scaleMineral} is supersaturated (SI ${scaleSI.toFixed(2)}) at ${xkm[iSI].toFixed(1)} km.` });
  if (sIn.I > 1 && v.actModel === 'davies') warnings.push({ level: 'warn', msg: `Ionic strength ${sIn.I.toFixed(2)} mol/L is beyond the range of the Davies equation; use the Truesdell–Jones model.` });
  if (Math.abs(ionicStrength(slow.water).balance) > 0.1) warnings.push({ level: 'info', msg: `The water analysis has a charge imbalance of ${(100 * ionicStrength(slow.water).balance).toFixed(0)} %; check the sodium or chloride value.` });
  if (sandBed) warnings.push({ level: 'warn', msg: `Sand settles over ${(bedLen / 1000).toFixed(1)} km: the mixture velocity is below the critical velocity (${sandCrit.toFixed(2)} m/s).` });
  if (aRisk !== 'low') warnings.push({ level: aRisk === 'high' ? 'bad' : 'warn', msg: `Asphaltene screening: de Boer "${db.cls}", colloidal instability index ${cii.cii.toFixed(2)} (${cii.cls})${iFh >= 0 ? `, Flory–Huggins onset near ${aPath[iFh].P.toFixed(0)} bara` : ''}.` });
  if (slow.restartDp > Math.max(p.pInMax - S.pOut, 1)) warnings.push({ level: 'bad', msg: `A gelled line would need ${slow.restartDp.toFixed(0)} bar to break the wax gel, more than the available ${(p.pInMax - S.pOut).toFixed(0)} bar.` });
  if (S.src !== 'flow suite') warnings.push({ level: 'info', msg: S.custom ? 'The flow picture was recomputed here with the shared kernel because the line or fluid inputs differ from the case.' : 'The flow suite has not been run: pressure, temperature and holdup come from the shared kernel estimate.' });
  if (maxSub > 0) {
    recs.push(wReqIdeal >= 94 ? `No practical ${inh.name} dose removes ${maxSub.toFixed(1)} °C of subcooling: shorten the shutdown, depressurise, or displace the line before it cools.` : `Inject ${inh.name} to ${wReqIdeal.toFixed(0)} wt % of the aqueous phase (about ${inhRate.toFixed(1)} m³/d at ${(fEvent * 100).toFixed(0)} % rate) to stay ${v.marginC} °C outside the hydrate region; the present ${S.aq.inhWt.toFixed(0)} wt % leaves ${maxSub.toFixed(1)} °C of subcooling.`);
    if (Number.isFinite(cooldown)) recs.push(`Treat ${cooldown.toFixed(1)} h as the cooldown limit: after a trip, start inhibitor displacement or depressurisation within about ${Math.max(cooldown - 2, 0).toFixed(1)} h (2 h reserved for the operation itself); the coldest point is at ${xkm[iCool].toFixed(1)} km.`);
    recs.push(`Screening for low-dosage inhibitors: a kinetic inhibitor is ${khiOk ? `plausible (subcooling ${maxSub.toFixed(1)} °C ≤ ${v.khiLimit} °C, indicative hold time ${khiHold > 1e4 ? 'very long' : khiHold.toFixed(0) + ' h'})` : `not suitable (subcooling ${maxSub.toFixed(1)} °C exceeds ${v.khiLimit} °C)`}; an anti-agglomerant is ${aaOk ? 'plausible' : 'not suitable'} at ${(prE.wcut * 100).toFixed(0)} % water cut — confirm either with qualification tests.`);
  } else recs.push(`The line stays ${(-maxSub).toFixed(1)} °C outside the hydrate region throughout this scenario; no hydrate inhibitor is needed for it. Check the shutdown and restart scenarios before relaxing the inhibition philosophy.`);
  if (plug || mc.prob > 0.1) recs.push(`If a plug forms, depressurise from both sides to ${v.depressP} bara: a ${v.plugLength} m plug then melts radially in about ${Number.isFinite(melt2.tNumeric) ? (melt2.tNumeric / 86400).toFixed(1) + ' d' : 'no finite time (the seabed is colder than the hydrate temperature at that pressure)'}. Never depressurise from one side only: ${Math.max(pLoc - v.depressP, 0).toFixed(0)} bar across the plug could launch it at the order of ${projV.toFixed(0)} m/s.`);
  if (piggingInterval !== null) recs.push(`Pig for wax every ${Math.max(Math.floor(piggingInterval * 0.8), 1)} d (80 % of the ${piggingInterval.toFixed(0)} d it takes to reach ${v.waxLimitMm} mm; ${slow.waxMass.toFixed(0)} kg of wax after ${v.depositDays} d).`);
  else recs.push(`No wax deposits in this operating mode: the wall stays above the wax appearance temperature of ${v.wat} °C.`);
  if (scaleSI > 0) recs.push(`Dose scale inhibitor against ${scaleMineral} (SI ${scaleSI.toFixed(2)}): a screening minimum inhibitor concentration is about ${mic.toFixed(0)} ppm — confirm with a dynamic tube-blocking test on the real brine.`);
  if (v.sandRate > 0) recs.push(sandBed ? `Raise the rate to at least ${sandMinRate.toFixed(0)} % of the case rate (mixture velocity ≥ ${sandCrit.toFixed(2)} m/s) to keep sand moving, or schedule sand pigging.` : `Sand keeps moving: the lowest velocity is ${sandMargin.toFixed(1)} times the critical velocity; do not run below about ${Math.min(sandMinRate, 100).toFixed(0)} % of the case rate for long periods.`);
  if (aRisk !== 'low') recs.push(`Run an asphaltene onset test on live oil (depressurisation from ${pRes.toFixed(0)} bara through the bubble point near ${pBub.toFixed(0)} bara) before selecting an inhibitor.`);

  // ---- plots
  const tAx = ser.t, pCurve = logspace(Math.max(Math.min(...h.P, S.pOut) * 0.6, 2), Math.max(...base.P, ...ser.pIn) * 1.15, 40), plots = [];
  plots.push({ type: 'line', title: 'Pressure–temperature path over the hydrate curve', xlabel: 'Temperature (°C)', ylabel: 'Pressure (bara)', logy: true, series: [
    { name: 'Hydrate curve (case water and inhibitor)', x: pCurve.map((P) => S.hT0(P) - depNow), y: pCurve }, { name: 'Hydrate curve (fresh water)', x: pCurve.map((P) => S.hT0(P)), y: pCurve, dash: true },
    { name: 'Steady production', x: base.T, y: base.P }, { name: `End of the simulated event (${h.t.toFixed(1)} h)`, x: h.T, y: h.P, mode: 'both' }], note: 'Hydrates are stable to the left of the curve.' });
  plots.push({ type: 'line', title: 'Subcooling and hydrate fraction along the line', xlabel: 'Distance (km)', ylabel: 'Subcooling (°C) · hydrate fraction (vol %)', series: [
    { name: 'Peak subcooling during the event', x: xkm, y: rec.subMax }, { name: 'Subcooling at the end', x: xkm, y: rec.sub, dash: true }, { name: 'Steady production', x: xkm, y: subSteady, dash: true }, { name: 'Peak hydrate fraction (vol % of liquid)', x: xkm, y: rec.phiMaxT.map((f) => f * 100) }], hlines: [{ y: 0, label: 'hydrate curve' }] });
  const dTg = linspace(1, 20, 39), tK = mean(S.tAmb) + KEL;
  plots.push({ type: 'line', title: 'Nucleation: induction time against subcooling', xlabel: 'Subcooling (°C)', ylabel: 'Mean induction time (h)', logy: true, ymin: 1e-3, ymax: 1e6, series: [
    { name: 'Heterogeneous (inputs)', x: dTg, y: dTg.map((d) => clamp(inductionTime(d, { TK: tK, sigma: p.sigma, theta: p.theta, A: p.nucA, het: true, V: p.nucV }) / 3600, 1e-3, 1e6)) },
    { name: 'Homogeneous (A = 10³⁵ 1/m³/s)', x: dTg, y: dTg.map((d) => clamp(inductionTime(d, { TK: tK, sigma: p.sigma, A: 1e35, het: false, V: p.nucV }) / 3600, 1e-3, 1e6)), dash: true }],
    vlines: maxSub > 0 ? [{ x: Math.min(maxSub, 20), label: 'peak subcooling' }] : [], note: `Mean waiting time for the first nucleus in ${v.nucVolume} L of water; the onset itself is exponentially distributed about this mean.` });
  if (pk.N) { const vt = pk.N.reduce((a, N, c) => a + N * grid.v[c], 0) || 1; plots.push({ type: 'line', title: `Particle-size distribution at the peak (${(S.x[pk.iPhi] / 1000).toFixed(1)} km, ${(pk.tPhi / 3600).toFixed(1)} h)`, xlabel: 'Particle / agglomerate size (µm)', ylabel: 'Volume fraction per class', logx: true, series: [{ name: 'Sectional population balance', x: grid.L.map((L) => L * 1e6), y: pk.N.map((N, c) => (N * grid.v[c]) / vt), mode: 'both' }], vlines: [{ x: p.dPrim * 1e6, label: 'primary' }, { x: Math.max(pk.dA[pk.iPhi], p.dPrim) * 1e6, label: 'cohesive limit' }] }); }
  plots.push({ type: 'line', title: 'Deposit thickness by type', xlabel: 'Distance (km)', ylabel: 'Thickness (mm)', zeroY: true, series: [{ name: 'Hydrate (worst time of the event)', x: xkm, y: dHyd.map((d) => d * 1000) }, { name: `Wax after ${v.depositDays} d`, x: xkm, y: slow.dWax.map((d) => d * 1000) }, { name: `Scale after ${v.depositDays} d`, x: xkm, y: slow.dScale.map((d) => d * 1000) }, { name: 'Total', x: xkm, y: total.map((d) => d * 1000), dash: true }] });
  plots.push({ type: 'line', title: 'Effective inner diameter', xlabel: 'Distance (km)', ylabel: 'Diameter (mm)', series: [{ name: 'With deposits', x: xkm, y: Deff.map((d) => d * 1000) }, { name: 'Clean bore', x: [xkm[0], xkm[n - 1]], y: [D0 * 1000, D0 * 1000], dash: true }] });
  plots.push({ type: 'line', title: 'Inlet pressure and blockage through the event', xlabel: 'Time (h)', ylabel: 'Inlet pressure (bara) · blockage (%) · viscosity factor', series: [{ name: 'Inlet pressure (bara)', x: tAx, y: ser.pIn }, { name: 'Peak area blockage (%)', x: tAx, y: ser.blk.map((b) => b * 100) }, { name: 'Peak slurry viscosity factor', x: tAx, y: ser.visc.map((x) => Math.min(x, 100)), dash: true }], hlines: [{ y: p.pInMax, label: 'available inlet pressure' }], vlines: tShut > 0 ? [{ x: tShut, label: 'restart' }] : [] });
  plots.push({ type: 'line', title: 'Hydrate inventory through the event', xlabel: 'Time (h)', ylabel: 'Hydrate mass (t)', zeroY: true, series: [{ name: 'Suspended in the stream', x: tAx, y: ser.susp.map((x) => x / 1000) }, { name: 'Deposited on the wall', x: tAx, y: ser.dep.map((x) => x / 1000) }] });
  if (mc.n > 0) { const data = mc.plugTimes.length >= 3 ? mc.plugTimes : mc.blk.map((b) => b * 100), hg = histogram(data, 10), isT = mc.plugTimes.length >= 3;
    plots.push({ type: 'bar', title: isT ? `Time to plug in ${mc.plugTimes.length} of ${mc.n} Monte Carlo samples` : `Peak hydrate blockage in ${mc.n} Monte Carlo samples`, ylabel: 'Samples', categories: hg.centers.map((c) => (isT ? c.toFixed(1) + ' h' : c.toFixed(1) + ' %')), series: [{ name: 'Samples', values: hg.counts }], note: `Plugging probability ${(mc.prob * 100).toFixed(0)} % (90 % interval ${(mc.lo * 100).toFixed(0)}–${(mc.hi * 100).toFixed(0)} %).` }); }
  plots.push({ type: 'line', title: 'Wax build-up at the worst location', xlabel: 'Time since the last pig run (d)', ylabel: 'Wax thickness (mm)', zeroY: true, series: [{ name: 'Maximum thickness', x: slow.ser.t, y: slow.ser.dMax }], hlines: [{ y: v.waxLimitMm, label: 'pigging limit' }] });
  const minIds = ['calcite', 'barite', 'celestite', 'gypsum', 'anhydrite', 'siderite'], nameOf = Object.fromEntries(MINERALS.map((mn) => [mn.id, mn.name])), siOf = (r, id) => Math.max(r.minerals.find((q) => q.id === id).SI, -6);
  plots.push({ type: 'line', title: 'Scale saturation index along the line', xlabel: 'Distance (km)', ylabel: 'Saturation index', series: minIds.map((id) => ({ name: nameOf[id], x: xkm, y: slow.sis.map((r) => siOf(r, id)) })), hlines: [{ y: 0, label: 'saturated' }] });
  plots.push({ type: 'line', title: 'Seawater mixing: saturation index at inlet conditions', xlabel: 'Seawater fraction (%)', ylabel: 'Saturation index', series: ['calcite', 'barite', 'celestite', 'gypsum'].map((id) => ({ name: nameOf[id], x: mixF.map((f) => f * 100), y: mix.map((r) => siOf(r, id)) })), hlines: [{ y: 0, label: 'saturated' }] });
  const rhoAx = linspace(550, 950, 9);
  plots.push({ type: 'line', title: 'de Boer asphaltene screening', xlabel: 'In-situ oil density (kg/m³)', ylabel: 'Reservoir pressure − saturation pressure (bar)', zeroY: true, series: [{ name: 'Slight problems above', x: rhoAx, y: rhoAx.map((r) => deBoer(r, 0).lower) }, { name: 'Severe problems above', x: rhoAx, y: rhoAx.map((r) => deBoer(r, 0).upper) }, { name: 'This oil', x: [clamp(aPath[0].rho, 550, 950)], y: [Math.max(pRes - pBub, 0)], mode: 'points' }], note: 'Boundaries are straight-line fits to the published screening plot.' });
  plots.push({ type: 'line', title: 'Sand: critical velocity against actual velocity', xlabel: 'Distance (km)', ylabel: 'Velocity (m/s)', zeroY: true, series: [{ name: 'Mixture velocity', x: xkm, y: sand.map((s) => s.vm) }, { name: 'Oroskar–Turian', x: xkm, y: sand.map((s) => s.cr.oroskarTurian), dash: true }, { name: 'Salama', x: xkm, y: sand.map((s) => s.cr.salama), dash: true }, { name: 'Danielson', x: xkm, y: sand.map((s) => s.cr.danielson), dash: true }] });
  if (h.fld.t.length >= 2) {
    plots.push({ type: 'field', title: 'Hydrate fraction in distance and time', xlabel: 'Distance (km)', ylabel: 'Time (h)', zlabel: 'Hydrate fraction', zunit: 'vol %', x: xkm, y: h.fld.t, z: h.fld.phi.map((r) => r.map((f) => f * 100)), cmap: 'viridis', markers: plug ? [{ x: plug.x / 1000, y: plug.t / 3600, label: 'plug' }] : [] });
    plots.push({ type: 'field', title: 'Subcooling in distance and time', xlabel: 'Distance (km)', ylabel: 'Time (h)', zlabel: 'Subcooling', zunit: '°C', x: xkm, y: h.fld.t, z: h.fld.sub, cmap: 'coolwarm' });
    plots.push({ type: 'field', title: 'Hydrate deposit thickness in distance and time', xlabel: 'Distance (km)', ylabel: 'Time (h)', zlabel: 'Thickness', zunit: 'mm', x: xkm, y: h.fld.t, z: h.fld.dep.map((r) => r.map((d) => d * 1000)), cmap: 'thermal' });
  }
  plots.push({ type: 'bar', title: 'Governing solids threat by zone', ylabel: 'Severity index (1 = at the limit)', categories: zones.map((z) => `${z.from.toFixed(1)}–${z.to.toFixed(1)} km`), series: ['Hydrate', 'Wax', 'Scale', 'Sand'].map((k) => ({ name: k, values: zones.map((z) => +z.sc[k].toFixed(3)) })) });

  // ---- tables
  const pick = Array.from({ length: Math.min(12, n) }, (_, k) => Math.round((k * (n - 1)) / Math.max(Math.min(12, n) - 1, 1))), tables = [];
  tables.push({ title: 'Line stations', columns: ['Distance (km)', 'Steady T (°C)', 'End T (°C)', 'End P (bara)', 'Hydrate T (°C)', 'Peak subcooling (°C)', 'Exposure (h)', 'Peak hydrate (vol %)', 'Agglomerate d43 (µm)', 'Hydrate deposit (mm)', 'Wax (mm)', 'Scale (mm)', 'Deposit permeability (m²)', 'Controlling step at the end'],
    rows: pick.map((i) => [tx(xkm[i], 4), tx(base.T[i]), tx(h.T[i]), tx(h.P[i]), tx(rec.teq[i]), tx(rec.subMax[i]), tx(rec.exposure[i]), tx(rec.phiMaxT[i] * 100), tx((pk.d43 ? pk.d43[i] : 0) * 1e6), tx(dHyd[i] * 1000), tx(slow.dWax[i] * 1000), tx(slow.dScale[i] * 1000), dHyd[i] > 0 ? tx(perm[i]) : '—', rec.lim[i] || '—']) });
  tables.push({ title: 'Hydrate formation, transport and plugging', columns: ['Quantity', 'Value', 'Unit'], rows: [
    ['Scenario', SC_LABEL[v.scenario], ''], ['Hydrate system', p.mode === 'oil' ? 'oil-dominated (shrinking-core droplets)' : p.mode === 'gas' ? 'gas-dominated (film and entrained water)' : 'water-dominated (absorption limited)', ''],
    ['Stability margin (most negative = safest)', tx(-maxSub), '°C'], ['Exposure inside the hydrate region', tx(expoH), 'h'], ['Exposure integral', tx(degH), '°C·h'], ['Hydrate-stable length at the end', tx(stableEnd / 1000), 'km'],
    ['Onset location', h.onset ? tx(h.onset.x / 1000) : '—', 'km'], ['Onset time', tOnset === null ? '—' : tx(tOnset), 'h'], ['Peak nucleation rate', tx(Math.max(...rec.J)), '1/m³ water/s'], ['Peak formation rate', tx(peakRate), 'kg/s'],
    ['Hydrate formed', tx(led.formed / 1000), 't'], ['Hydrate dissociated', tx(led.dissociated / 1000), 't'], ['Hydrate carried out of the line', tx(led.exported / 1000), 't'], ['Captured on the wall', tx(led.captured / 1000), 't'], ['Sloughed back into the stream', tx(led.sloughed / 1000), 't'],
    ['Water converted (net)', tx(waterUsed / 1000), 't'], ['Gas consumed', tx(gasUsed / 1000), 't'], ['Peak water conversion in a cell', tx(convMax * 100), '%'], ['Heat of formation released (net)', tx(led.heat / 1e9), 'GJ'], ['Hydration number', tx(p.hydN, 4), 'mol/mol'],
    ['Primary particle size', tx(p.dPrim * 1e6), 'µm'], ['Cohesive-limit agglomerate size at the peak', tx((pk.dA ? pk.dA[pk.iPhi] : 0) * 1e6), 'µm'], ['Agglomerate d43 at the peak', tx((pk.d43 ? pk.d43[pk.iPhi] : 0) * 1e6), 'µm'], ['Effective volume fraction at the peak', tx(pk.phiE ? pk.phiE[pk.iPhi] : 0), '–'], ['Peak slurry viscosity factor', tx(peakVisc), '×'],
    ['Peak wall-capture rate', tx(Math.max(...ser.dep.map((d, k) => (k ? (d - ser.dep[k - 1]) / Math.max((ser.t[k] - ser.t[k - 1]) * 3600, 1) : 0)), 0)), 'kg/s'], ['Inlet-pressure rise in the event', tx(dpRise), 'bar'], ['Pressure-drop increase from wax and scale', tx(dpFoul - dpClean), 'bar'],
    ['Effective roughness with deposits', tx(roughEff * 1e6), 'µm'], ['Plug', plug ? `${(plug.t / 3600).toFixed(2)} h at ${(plug.x / 1000).toFixed(2)} km — ${plug.mech}` : 'none in the base case', ''],
    ['Plugging probability (90 % interval)', `${(mc.prob * 100).toFixed(0)} % (${(mc.lo * 100).toFixed(0)}–${(mc.hi * 100).toFixed(0)} %)`, ''], ['Time to plug P10 / P50 / P90', plugT.length ? `${quantile(plugT, 0.1).toFixed(1)} / ${quantile(plugT, 0.5).toFixed(1)} / ${quantile(plugT, 0.9).toFixed(1)}` : '—', 'h'], ['Risk index', tx(risk), '0–1']],
    note: 'Times are measured from the start of the simulated sequence (the start of the shut-in for a cold restart).' });
  tables.push({ title: 'Population balance cross-check: sectional against quadrature method of moments', columns: ['Quantity', 'Sectional (fixed pivot)', 'QMOM (3 nodes)'], rows: [['Batch time (s)', tx(qm.tB), tx(qm.tB)], ['Number remaining N/N₀', tx(qm.m0s, 4), qm.ok ? tx(qm.m0q, 4) : '—'], ['Third moment m₃/m₃₀ (volume)', tx(qm.m3s, 6), qm.ok ? tx(qm.m3q, 6) : '—'], ['d43 (µm)', tx(qm.d43s * 1e6, 4), qm.ok ? tx(qm.d43q * 1e6, 4) : '—']], note: 'Pure aggregation of a narrow population with the collision kernel of the worst cell over three collision times; both methods must conserve the third moment.' });
  tables.push({ title: 'Inhibition and remediation', columns: ['Item', 'Value', 'Unit', 'Basis'], rows: [
    ['Depression from salt', tx(depSalt), '°C', 'Nielsen–Bucklin on the water mole fraction'], [`Depression from ${S.aq.inhWt > 0 ? S.aq.inh.name : 'inhibitor'} now`, tx(depNow - depSalt), '°C', `${S.aq.inhWt.toFixed(1)} wt %, effectiveness ${v.inhEff} %`], ['Inhibitor effectiveness', tx(maxSub > 0 || depNow > depSalt ? (100 * (depNow - depSalt)) / Math.max(maxSub + depNow - depSalt, 1e-9) : 100), '% of the uninhibited subcooling removed', 'event peak'],
    [`${inh.name} required`, tx(wReqIdeal), 'wt %', `${v.marginC} °C margin on the peak subcooling`], ['Injection rate', tx(inhRate), 'm³/d', `${(mWater * 86.4).toFixed(1)} t/d of water at the event rate`], ['Loss to gas and hydrocarbon liquid', tx(lossPct), '% of injected', inhId === 'MeOH' ? 'partition coefficients (screening)' : 'negligible for glycols (screening)'],
    ['Kinetic inhibitor (screening)', khiOk ? 'plausible' : 'not suitable', '', `subcooling limit ${v.khiLimit} °C`], ['Anti-agglomerant (screening)', aaOk ? 'plausible' : 'not suitable', '', `water-cut limit ${v.aaWcLimit} %`],
    ['Two-sided depressurisation: radial melt time', Number.isFinite(melt2.tNumeric) ? tx(melt2.tNumeric / 3600) : 'no melting', 'h', `Stefan front, hydrate at ${tdTwo.toFixed(1)} °C, ambient ${tA.toFixed(1)} °C`], ['One-sided depressurisation: axial melt time', Number.isFinite(melt1) ? tx(melt1 / 86400) : 'no melting', 'd', `Neumann solution over ${v.plugLength} m — not recommended`],
    ['Plug velocity if released one-sided', tx(projV), 'm/s', 'after 100 m of free travel, friction neglected'], ['Intrinsic dissociation time of a grain', Number.isFinite(tKin) ? tx(tKin) : '—', 's', 'Kim–Bishnoi: far faster than heat supply, so melting is heat-transfer-controlled'],
    ['Methanol to dissolve the plug', tx(meohMelt), 'm³', `${(wEq * 100).toFixed(0)} wt % in the released water`], ['Methanol contact time (mass-transfer-controlled)', tx(tMeoh / 86400), 'd', 'film coefficient 10⁻⁵ m/s on the plug face'], ['Heating to hold the line outside the region', tx(dehW / 1000), 'kW', `${tHold.toFixed(1)} °C along ${(S.L / 1000).toFixed(1)} km`], ['Heat to melt the plug in one day', tx(meltKW), 'kW', 'latent heat only']] });
  tables.push({ title: 'Wax', columns: ['Quantity', 'Value', 'Unit'], rows: [['Wax appearance temperature', v.wat, '°C'], ['Length with a wall below the WAT', tx(slow.last.filter((r) => r && r.Ti < v.wat).length * ds / 1000), 'km'], ['Initial build-up rate', tx(slow.waxRate0), 'mm/d'], [`Maximum thickness after ${v.depositDays} d`, tx(slow.ser.dMax[slow.ser.dMax.length - 1]), 'mm'], ['Wax mass in the line', tx(slow.waxMass), 'kg'], ['Wax fraction of the aged deposit', tx(Math.max(...slow.Fw)), '–'], ['Wax diffusivity at the wall', tx(Math.max(...slow.last.map((r) => r?.Dwo || 0))), 'm²/s'], ['Pigging interval', piggingInterval === null ? 'not needed' : tx(piggingInterval), 'd'], ['Gelled length after a cold shutdown', tx(slow.gelLen / 1000), 'km'], ['Gel yield stress', tx(slow.tauY), 'Pa'], ['Gel-breaking restart pressure', tx(slow.restartDp), 'bar']] });
  tables.push({ title: 'Scale at the inlet and at the worst location', columns: ['Mineral', 'SI at the inlet', 'SI at the worst location', 'Precipitation potential at the inlet (mg/L)', 'SI with 50 % seawater'], rows: minIds.map((id, k) => [nameOf[id], tx(siOf(sIn, id)), tx(Math.max(...slow.sis.map((r) => siOf(r, id)))), tx(sIn.minerals[k].ptb), tx(siOf(mix[10], id))]),
    note: `Ionic strength ${sIn.I.toFixed(2)} mol/L, pH ${sIn.pH.toFixed(2)} at ${pMix.toFixed(0)} bara and ${tMix.toFixed(0)} °C with ${v.co2Pct} mol % CO₂; Oddo–Tomson calcite index ${sIn.oddoTomson.toFixed(2)} for comparison. Maximum deposition rate ${Math.max(...slow.scaleRate).toFixed(3)} mm/y. Ion pairing is neglected, which is conservative for the sulphates.` });
  tables.push({ title: 'Asphaltene and sand', columns: ['Quantity', 'Value', 'Unit'], rows: [['In-situ oil density at reservoir conditions', tx(aPath[0].rho), 'kg/m³'], ['Undersaturation', tx(Math.max(pRes - pBub, 0)), 'bar'], ['de Boer class', db.cls, ''], ['Colloidal instability index', tx(cii.cii), cii.cls], ['Lowest Flory–Huggins solubility on the depressurisation path', tx(fhMin), 'vol fraction'], ['Asphaltene in the oil', tx(aspPhi), 'vol fraction'], ['Flory–Huggins onset pressure', iFh >= 0 ? tx(aPath[iFh].P) : 'none', 'bara'], ['Asphaltene deposition risk', aRisk, ''],
    ['Sand concentration in the liquid', tx(sand[0].C * 1e6), 'ppm by volume'], ['Sand settling velocity (hindered)', tx(Math.abs(sand[0].st.vHindered)), 'm/s'], ['Particle response time', tx(relax.tau), 's'], ['Critical velocity (governing)', tx(sandCrit), 'm/s'], ['Lowest velocity ratio v / v_critical', tx(sandMargin), '–'], ['Length with a sand bed', tx(bedLen / 1000), 'km'], ['Sand hold-up in the line', tx(sandInv), 'kg'], ['Minimum rate to keep sand moving', tx(sandMinRate), '% of case rate'], ['Bend erosion (screening)', tx(erosion), 'mm/y']] });
  tables.push({ title: 'Governing threat by zone', columns: ['From (km)', 'To (km)', 'Hydrate', 'Wax', 'Scale', 'Sand', 'Governing'], rows: zones.map((z) => [tx(z.from), tx(z.to), tx(z.sc.Hydrate), tx(z.sc.Wax), tx(z.sc.Scale), tx(z.sc.Sand), z.top]), note: 'Severity indices: 1 means at the limit (10 °C subcooling or 15 % of the bore for hydrate, the pigging thickness for wax, SI 1.5 or 2 mm/y for scale, velocity at the critical velocity for sand).' });
  tables.push({ title: 'Hydrate thermodynamics check (van der Waals–Platteeuw, methane structure I)', columns: ['Quantity', 'Value', 'Unit'], rows: [['Mean seabed temperature', tx(mean(S.tAmb)), '°C'], ['Methane hydrate equilibrium pressure', vdw.P === null ? '—' : tx(vdw.P), 'bara'], ['Screening curve of the case gas at that temperature', tx(S.peq(mean(S.tAmb))), 'bara'], ['Small-cage occupancy at line pressure', tx(occ.thetaS), '–'], ['Large-cage occupancy at line pressure', tx(occ.thetaL), '–'], ['Hydration number from occupancy', tx(occ.hydrationNumber, 4), 'mol/mol']], note: 'The case gas contains propane and butanes and forms structure II at a lower pressure than pure methane; the statistical model is used for the cage occupancy and as an upper bound on the equilibrium pressure.' });

  const flowArea = Deff.reduce((a, d) => a + (PI / 4) * d * d * ds, 0), depArea = total.reduce((a, d, i) => a + (PI / 4) * (D0 * D0 - Deff[i] ** 2) * ds, 0);
  const balances = [
    { name: 'Hydrate mass (kg): formed + initial + inflow = suspended + deposited + exported + dissociated', in: led.in, out: led.out },
    { name: 'Water bound in hydrate (kg) = hydration-number share of the net hydrate formed', in: (led.formed - led.dissociated) * ((p.hydN * MW_W) / (p.mwG + p.hydN * MW_W)), out: waterUsed },
    { name: 'Pipe volume (m³): flow area + deposit area = clean bore', in: A0 * S.L, out: flowArea + depArea },
  ];
  const summary = `${SC_LABEL[v.scenario]}: ${maxSub > 0 ? `up to ${maxSub.toFixed(1)} °C of subcooling over ${(stableLen / 1000).toFixed(1)} km, hydrate reaches ${(peakPhi * 100).toFixed(1)} vol % of the liquid and ${plug ? `the line plugs after ${(plug.t / 3600).toFixed(1)} h at ${(plug.x / 1000).toFixed(1)} km` : `no plug forms in the base case (plugging probability ${(mc.prob * 100).toFixed(0)} %)`}` : `the line stays ${(-maxSub).toFixed(1)} °C outside the hydrate region`}; wax ${slow.waxRate0 > 1e-4 ? `builds at ${slow.waxRate0.toFixed(3)} mm/d` : 'does not deposit'}, ${scaleSI > 0 ? `${scaleMineral} is supersaturated (SI ${scaleSI.toFixed(2)})` : 'no mineral is supersaturated'}, asphaltene risk is ${aRisk} and sand ${sandBed ? 'settles' : 'keeps moving'}.`;
  prog(1, 'Done');
  return { summary, kpis, warnings, recommendations: recs, plots, tables, balances, outputs };
}

export const runSolids = run;
/*__PART5__*/
