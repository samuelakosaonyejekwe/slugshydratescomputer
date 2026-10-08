// Suite 6 — Integrity, Loads, Risk & Engineering Assessment.
// Turns the pressure, temperature, slug and solids picture of the case into engineering consequences: pipe stress and
// code utilisation, collapse and buckling, slug loads and span dynamics (beam finite elements, Newmark), fatigue and
// crack growth, CO2 corrosion and sand erosion, defect assessment, structural reliability (FORM/SORM/Monte Carlo)
// and risk (fault tree, event tree, Markov, Bayesian network, FMECA, risk matrix). SI inside; MPa, bara, °C, mm at the interfaces.
import { clamp, linspace, logspace, sum, mean, isNum, interp1, brent, rng, lhs, histogram, fmt } from '../core/num.js';
import { fluidModel, waterContent } from '../core/thermo.js';
import { G, slugUnit, frictionFactor } from '../core/pipe.js';
import { flowPicture, caseLine } from '../core/caseflow.js';
import { BASE } from '../data/basecase.js';

const YEAR = 365.25 * 86400, RHO_SW = 1025, PATM = 1.01325e5, RGAS = 8.314462618, FARADAY = 96485.33212, MPA = 1e6, BAR = 1e5;
const pos = (x, d) => (isNum(+x) && x !== null && x !== '' && +x > 0 ? +x : d);
const nz = (x, d = 0) => (isNum(+x) && x !== null && x !== '' ? +x : d);
const fin = (x, d = 0) => (Number.isFinite(x) ? x : d);
const sig = (x, n = 4) => (Number.isFinite(x) ? (x === 0 ? 0 : Number(x.toPrecision(n))) : null);
const cap = (x, hi = 1e6) => (Number.isFinite(x) ? Math.min(x, hi) : hi);
const need = (cond, msg) => { if (!cond) throw new Error(msg); };
const dot = (a, b) => { let s = 0; for (let i = 0; i < a.length; i++) s += a[i] * b[i]; return s; };
const norm2 = (a) => Math.sqrt(dot(a, a));

// ---- probability functions --------------------------------------------------------------------------------------------
/** Standard normal cumulative distribution (Hart's rational approximation, double precision in the tails). */
export function Phi(x) {
  const ax = Math.abs(x);
  let c;
  if (ax > 37) c = 0;
  else {
    const e = Math.exp((-ax * ax) / 2);
    if (ax < 7.07106781186547) {
      let b = 3.52624965998911e-2 * ax + 0.700383064443688; b = b * ax + 6.37396220353165; b = b * ax + 33.912866078383; b = b * ax + 112.079291497871; b = b * ax + 221.213596169931; b = b * ax + 220.206867912376;
      c = e * b;
      b = 8.83883476483184e-2 * ax + 1.75566716318264; b = b * ax + 16.064177579207; b = b * ax + 86.7807322029461; b = b * ax + 296.564248779674; b = b * ax + 637.333633378831; b = b * ax + 793.826512519948; b = b * ax + 440.413735824752;
      c /= b;
    } else { let b = ax + 0.65; b = ax + 4 / b; b = ax + 3 / b; b = ax + 2 / b; b = ax + 1 / b; c = e / b / 2.506628274631; }
  }
  return x > 0 ? 1 - c : c;
}
const phi = (x) => Math.exp(-0.5 * x * x) / 2.5066282746310002;
/** Inverse standard normal distribution (Acklam's algorithm with one Halley refinement). */
export function PhiInv(p) {
  if (!(p > 0)) return -38; if (!(p < 1)) return 38;
  const a = [-39.69683028665376, 220.9460984245205, -275.9285104469687, 138.357751867269, -30.66479806614716, 2.506628277459239], b = [-54.47609879822406, 161.5858368580409, -155.6989798598866, 66.80131188771972, -13.28068155288572];
  const c = [-7.784894002430293e-3, -0.3223964580411365, -2.400758277161838, -2.549732539343734, 4.374664141464968, 2.938163982698783], d = [7.784695709041462e-3, 0.3224671290700398, 2.445134137142996, 3.754408661907416];
  let x;
  if (p < 0.02425) { const q = Math.sqrt(-2 * Math.log(p)); x = (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1); }
  else if (p <= 0.97575) { const q = p - 0.5, r = q * q; x = ((((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q) / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1); }
  else { const q = Math.sqrt(-2 * Math.log(1 - p)); x = -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1); }
  const e = Phi(x) - p, u = e * 2.5066282746310002 * Math.exp((x * x) / 2);
  return Number.isFinite(u) ? x - u / (1 + (x * u) / 2) : x;
}
const LANCZOS = [0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];
/** Gamma function (Lanczos approximation). */
export function gammaFn(z) {
  if (z < 0.5) return Math.PI / (Math.sin(Math.PI * z) * gammaFn(1 - z));
  z -= 1; let x = LANCZOS[0]; for (let i = 1; i < 9; i++) x += LANCZOS[i] / (z + i);
  const t = z + 7.5; return Math.sqrt(2 * Math.PI) * t ** (z + 0.5) * Math.exp(-t) * x;
}
/**
 * Random variable { name, dist: 'normal' | 'lognormal' | 'weibull' | 'gumbel' | 'uniform' | 'det', mean, cov | sd }.
 * Returns { name, dist, mean, sd, x(u) (iso-probabilistic transform from a standard normal u), cdf(x) }.
 */
export function randomVariable(v) {
  const m = +v.mean, sd = v.sd !== undefined ? Math.abs(+v.sd) : Math.abs(m * (+v.cov || 0)), base = { name: v.name || 'X', dist: v.dist || 'normal', mean: m, sd };
  need(Number.isFinite(m) && Number.isFinite(sd), `Random variable "${base.name}" needs a finite mean and scatter.`);
  if (base.dist === 'det' || sd === 0) return { ...base, dist: 'det', sd: 0, x: () => m, cdf: (x) => (x >= m ? 1 : 0) };
  switch (base.dist) {
    case 'lognormal': { need(m > 0, `Lognormal variable "${base.name}" needs a positive mean.`); const z2 = Math.log(1 + (sd / m) ** 2), z = Math.sqrt(z2), l = Math.log(m) - z2 / 2; return { ...base, x: (u) => Math.exp(l + z * u), cdf: (x) => (x > 0 ? Phi((Math.log(x) - l) / z) : 0) }; }
    case 'gumbel': { const b = (sd * Math.sqrt(6)) / Math.PI, a = m - 0.5772156649015329 * b; return { ...base, x: (u) => a - b * Math.log(u > 0 ? -Math.log1p(-Phi(-u)) : -Math.log(Phi(u))), cdf: (x) => Math.exp(-Math.exp(-(x - a) / b)) }; }
    case 'weibull': {
      need(m > 0, `Weibull variable "${base.name}" needs a positive mean.`);
      const cv = sd / m, f = (k) => gammaFn(1 + 2 / k) / gammaFn(1 + 1 / k) ** 2 - 1 - cv * cv, k = brent(f, 0.15, 200, 1e-12), lam = m / gammaFn(1 + 1 / k);
      return { ...base, shape: k, scale: lam, x: (u) => lam * (u < 0 ? -Math.log1p(-Phi(u)) : -Math.log(Phi(-u))) ** (1 / k), cdf: (x) => (x > 0 ? 1 - Math.exp(-((x / lam) ** k)) : 0) };
    }
    case 'uniform': { const h = sd * Math.sqrt(3); return { ...base, x: (u) => m - h + 2 * h * Phi(u), cdf: (x) => clamp((x - m + h) / (2 * h), 0, 1) }; }
    default: return { ...base, dist: 'normal', x: (u) => m + sd * u, cdf: (x) => Phi((x - m) / sd) };
  }
}

// ---- dense linear algebra for the beam model and the reliability methods ---------------------------------------------
const zeros = (n, m = n) => Array.from({ length: n }, () => new Float64Array(m));
/** Cholesky factor (lower) of a symmetric positive-definite matrix, or null when it is not positive definite. */
function chol(A) {
  const n = A.length, L = zeros(n);
  for (let i = 0; i < n; i++) for (let j = 0; j <= i; j++) {
    let s = A[i][j]; const Li = L[i], Lj = L[j];
    for (let k = 0; k < j; k++) s -= Li[k] * Lj[k];
    if (i === j) { if (!(s > 1e-13 * Math.abs(A[i][i]))) return null; Li[i] = Math.sqrt(s); } else Li[j] = s / Lj[j];
  }
  return L;
}
const fwd = (L, b) => { const n = b.length, y = new Float64Array(n); for (let i = 0; i < n; i++) { let s = b[i]; const Li = L[i]; for (let k = 0; k < i; k++) s -= Li[k] * y[k]; y[i] = s / Li[i]; } return y; };
const bwd = (L, y) => { const n = y.length, x = new Float64Array(n); for (let i = n - 1; i >= 0; i--) { let s = y[i]; for (let k = i + 1; k < n; k++) s -= L[k][i] * x[k]; x[i] = s / L[i][i]; } return x; };
const cholSolve = (L, b) => bwd(L, fwd(L, b));
const matVec = (A, x) => { const n = A.length, y = new Float64Array(n); for (let i = 0; i < n; i++) { let s = 0; const Ai = A[i]; for (let j = 0; j < x.length; j++) s += Ai[j] * x[j]; y[i] = s; } return y; };
/** Eigenvalues (ascending) and eigenvectors (columns of V) of a symmetric matrix by cyclic Jacobi rotations. */
export function jacobiEig(A) {
  const n = A.length, a = A.map((r) => Float64Array.from(r)), V = zeros(n);
  for (let i = 0; i < n; i++) V[i][i] = 1;
  let scale = 0; for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) scale += a[i][j] * a[i][j];
  for (let sweep = 0; sweep < 80; sweep++) {
    let off = 0; for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) off += a[i][j] * a[i][j];
    if (off <= 1e-26 * scale) break;
    for (let p = 0; p < n - 1; p++) for (let q = p + 1; q < n; q++) {
      const apq = a[p][q]; if (Math.abs(apq) < 1e-300) continue;
      const th = (a[q][q] - a[p][p]) / (2 * apq), t = (th >= 0 ? 1 : -1) / (Math.abs(th) + Math.sqrt(th * th + 1)), c = 1 / Math.sqrt(t * t + 1), s = t * c;
      for (let k = 0; k < n; k++) { const kp = a[k][p], kq = a[k][q]; a[k][p] = c * kp - s * kq; a[k][q] = s * kp + c * kq; }
      for (let k = 0; k < n; k++) { const pk = a[p][k], qk = a[q][k]; a[p][k] = c * pk - s * qk; a[q][k] = s * pk + c * qk; }
      for (let k = 0; k < n; k++) { const kp = V[k][p], kq = V[k][q]; V[k][p] = c * kp - s * kq; V[k][q] = s * kp + c * kq; }
    }
  }
  const order = Array.from({ length: n }, (_, i) => i).sort((i, j) => a[i][i] - a[j][j]);
  return { values: order.map((i) => a[i][i]), vectors: order.map((i) => V.map((r) => r[i])) }; // vectors[k] = k-th eigenvector
}
/** Generalised symmetric eigenproblem K x = λ M x (M positive definite): all pairs, x M-orthonormal. */
function genEig(K, M) {
  const n = K.length, L = chol(M); need(L, 'The mass matrix of the beam model is not positive definite.');
  const B = zeros(n); // B = L⁻¹ K L⁻ᵀ
  const cols = []; for (let j = 0; j < n; j++) cols.push(fwd(L, K.map((r) => r[j])));
  for (let i = 0; i < n; i++) { const row = fwd(L, cols.map((c) => c[i])); for (let j = 0; j < n; j++) B[i][j] = row[j]; }
  for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) B[i][j] = B[j][i] = 0.5 * (B[i][j] + B[j][i]);
  const e = jacobiEig(B);
  return { values: e.values, vectors: e.vectors.map((y) => bwd(L, y)) };
}
/** Lowest p eigenpairs of K x = λ M x by subspace iteration (K positive definite). */
function subspaceEig(K, M, p) {
  const n = K.length, q = Math.min(n, Math.max(2 * p, p + 8));
  if (n <= 20 || q >= n) { const e = genEig(K, M); return { values: e.values.slice(0, p), vectors: e.vectors.slice(0, p), iterations: 0 }; }
  const Lk = chol(K); need(Lk, 'The beam is a mechanism or has buckled: add supports or reduce the compressive axial force.');
  const r = rng(4711); let X = Array.from({ length: q }, (_, k) => Float64Array.from({ length: n }, (_, i) => (k === 0 ? 1 : r.uniform(-1, 1))));
  let last = null, vals = [], it = 0;
  for (; it < 60; it++) {
    const Y = X.map((x) => matVec(M, x)), Xb = Y.map((y) => cholSolve(Lk, y)), MXb = Xb.map((x) => matVec(M, x));
    const Kr = zeros(q), Mr = zeros(q);
    for (let i = 0; i < q; i++) for (let j = i; j < q; j++) { Kr[i][j] = Kr[j][i] = dot(Xb[i], Y[j]); Mr[i][j] = Mr[j][i] = dot(Xb[i], MXb[j]); }
    const e = genEig(Kr, Mr);
    X = e.vectors.map((c) => { const x = new Float64Array(n); for (let k = 0; k < q; k++) { const ck = c[k], xb = Xb[k]; for (let i = 0; i < n; i++) x[i] += ck * xb[i]; } return x; });
    vals = e.values;
    if (last && vals.slice(0, p).every((v, i) => Math.abs(v - last[i]) <= 1e-11 * Math.abs(v))) break;
    last = vals.slice();
  }
  return { values: vals.slice(0, p), vectors: X.slice(0, p), iterations: it };
}

// ---- pipe stress --------------------------------------------------------------------------------------------------------
/** Lamé thick-cylinder stresses (Pa) at radius r for internal pressure pi and external pressure pe: { radial, hoop, axial (capped ends) }. */
export function lame(pi, pe, ri, ro, r) {
  const d = ro * ro - ri * ri, a = (pi * ri * ri - pe * ro * ro) / d, b = ((pi - pe) * ri * ri * ro * ro) / d;
  return { radial: a - b / (r * r), hoop: a + b / (r * r), axial: a };
}
/** Von Mises equivalent of three normal stresses and one shear stress. */
export const vonMises = (s1, s2, s3 = 0, tau = 0) => Math.sqrt(0.5 * ((s1 - s2) ** 2 + (s2 - s3) ** 2 + (s3 - s1) ** 2) + 3 * tau * tau);
/** Tresca equivalent (largest difference of the principal stresses). */
export const tresca = (s1, s2, s3 = 0) => Math.max(s1, s2, s3) - Math.min(s1, s2, s3);
/**
 * Stress state of a pipe wall (all SI: m, Pa, N, N·m, K).
 * o: { D (steel outer diameter), t, pi, pe, dT (temperature above installation), E, nu, alpha, restrained (bool), axial (true wall force, tension +), moment }
 * Returns { hoopThin (Barlow, outer diameter), hoopMean (mean diameter), hoopInner, hoopOuter, radialInner, radialOuter, thermal (−EαΔT),
 *   endCap, poisson, bending, longitudinal, vonMises, tresca, area, I, Z } with hoop/radial from Lamé at the bore for the equivalents.
 */
export function pipeStress(o) {
  const D = +o.D, t = +o.t, pi = +o.pi || 0, pe = +o.pe || 0, ro = D / 2, ri = ro - t;
  need(D > 0 && t > 0 && t < D / 2, 'Pipe stress needs a positive wall thickness smaller than the pipe radius.');
  const As = Math.PI * (ro * ro - ri * ri), I = (Math.PI / 4) * (ro ** 4 - ri ** 4), Z = I / ro;
  const inner = lame(pi, pe, ri, ro, ri), outer = lame(pi, pe, ri, ro, ro), hoopThin = ((pi - pe) * D) / (2 * t), hoopMean = ((pi - pe) * (D - t)) / (2 * t);
  const thermal = -(+o.E || 0) * (+o.alpha || 0) * (+o.dT || 0), endCap = inner.axial, poisson = (+o.nu || 0) * (inner.hoop + inner.radial); // plane strain: ν(σθ + σr) is uniform through the wall
  const direct = (+o.axial || 0) / As, bending = Math.abs(+o.moment || 0) / Z, base = (o.restrained ? poisson + thermal : endCap) + direct;
  const cand = [base + bending, base - bending].map((sl) => ({ sl, vm: vonMises(inner.hoop, sl, inner.radial), tr: tresca(inner.hoop, sl, inner.radial) })), worst = cand[0].vm >= cand[1].vm ? cand[0] : cand[1];
  return { hoopThin, hoopMean, hoopInner: inner.hoop, hoopOuter: outer.hoop, radialInner: inner.radial, radialOuter: outer.radial, thermal, endCap, poisson, bending, longitudinal: worst.sl, vonMises: worst.vm, tresca: Math.max(cand[0].tr, cand[1].tr), area: As, I, Z };
}
const B318_T = [[121, 1], [149, 0.967], [177, 0.933], [204, 0.9], [232, 0.867]];
const DNV_SC = { low: { pc: 1.046, lb: 1.04, target: 1e-3 }, medium: { pc: 1.138, lb: 1.14, target: 1e-4 }, high: { pc: 1.308, lb: 1.26, target: 1e-5 }, veryHigh: { pc: 1.308, lb: 1.26, target: 1e-6 } };
const DNV_GM = 1.15;
/** Yield derating of C-Mn line pipe with temperature in the DNV-ST-F101 format (MPa): none to 50 °C, 30 MPa at 100 °C, 70 MPa at 200 °C. */
export const dnvDerating = (T) => (T <= 50 ? 0 : T <= 100 ? 0.6 * (T - 50) : Math.min(30 + 0.4 * (T - 100), 110));
/**
 * Pressure-containment design format. o: { code: 'b314' | 'b318' | 'dnv', D, t (m), smys, smts (Pa), T (°C), factor (ASME design factor), safetyClass }
 * Returns { allowDp (Pa, allowable internal minus external pressure for wall t), tReq(dp) (m), incidental (ratio of the checked pressure to
 * the design pressure), hoopAllow (Pa), longAllow, combAllow (Pa), basis (text) }.
 */
export function pressureDesign(o) {
  const D = +o.D, t = Math.max(+o.t, 1e-6), S = +o.smys, T = nz(o.T, 20);
  if (o.code === 'dnv') {
    const sc = DNV_SC[o.safetyClass] || DNV_SC.medium, de = dnvDerating(T) * MPA, fy = (S - de) * 0.96, fu = (+o.smts - de) * 0.96, fcb = Math.min(fy, fu / 1.15), k = (fcb * 2) / Math.sqrt(3) / (DNV_GM * sc.pc);
    return { allowDp: (2 * t * k) / (D - t), tReq: (dp) => (Math.max(dp, 0) * D) / (2 * k + Math.max(dp, 0)), incidental: 1.1, hoopAllow: k, longAllow: fy, combAllow: fy, fy, fu, basis: `DNV-ST-F101 pressure containment: p_li − p_e ≤ 2t/(D−t)·f_cb·(2/√3)/(γ_m·γ_SC), γ_m = 1.15, γ_SC = ${sc.pc}, α_U = 0.96` };
  }
  const F = clamp(nz(o.factor, 0.72), 0.2, 1), Td = o.code === 'b318' ? interp1(B318_T.map((r) => r[0]), B318_T.map((r) => r[1]), T) : 1, k = F * S * Td;
  return { allowDp: (2 * t * k) / D, tReq: (dp) => (Math.max(dp, 0) * D) / (2 * k), incidental: 1, hoopAllow: k, longAllow: 0.8 * S, combAllow: 0.9 * S * Td, fy: S, fu: +o.smts, basis: o.code === 'b318' ? `ASME B31.8 (offshore chapter): (p_i − p_e)·D/(2t) ≤ F·S·T, F = ${F}, T = ${sig(Td, 3)}; longitudinal 0.80 S, combined 0.90 S` : `ASME B31.4 (offshore chapter): (p_i − p_e)·D/(2t) ≤ F·SMYS, F = ${F}; longitudinal 0.80 SMYS, combined 0.90 SMYS` };
}

// ---- collapse and buckling ----------------------------------------------------------------------------------------------
/**
 * External-pressure capacities (Pa). o: { D, t, E, nu, fy, ovality (f0 = (Dmax − Dmin)/D), alphaFab }
 * Returns { pel (elastic 2E(t/D)³/(1−ν²)), pp (plastic), pc (combined collapse with ovality, DNV-ST-F101 cubic), ppr (propagating buckle 35·fy·αfab·(t/D)^2.5) }.
 */
export function collapsePressure(o) {
  const D = +o.D, t = +o.t, fy = +o.fy, af = nz(o.alphaFab, 1), f0 = Math.max(nz(o.ovality, 0.005), 0.005);
  need(D > 0 && t > 0 && fy > 0, 'Collapse needs positive diameter, wall thickness and yield strength.');
  const pel = (2 * +o.E * (t / D) ** 3) / (1 - o.nu * o.nu), pp = (fy * af * 2 * t) / D;
  // root of the cubic (pc − pel)(pc² − pp²) = pc·pel·pp·f0·D/t below both pel and pp, in closed (trigonometric) form
  const b = -pel, c = -(pp * pp + (pp * pel * f0 * D) / t), d = pel * pp * pp, u = (c - (b * b) / 3) / 3, v = 0.5 * ((2 * b ** 3) / 27 - (b * c) / 3 + d), ph = Math.acos(clamp(-v / Math.sqrt(-u * u * u), -1, 1)), pc = -2 * Math.sqrt(-u) * Math.cos(ph / 3 + Math.PI / 3) - b / 3;
  return { pel, pp, pc, ppr: 35 * fy * af * (t / D) ** 2.5 };
}
/**
 * Local buckling under combined loading, load-controlled format of DNV-ST-F101 (design loads in, utilisation ≤ 1 out).
 * o: { D, t, fy, fu, E, nu, pi, pe, pmin (minimum internal pressure), M (design moment), S (design effective axial force), safetyClass, ovality }
 */
export function localBuckling(o) {
  const D = +o.D, t = +o.t, fy = +o.fy, fu = +o.fu, sc = DNV_SC[o.safetyClass] || DNV_SC.medium, g = DNV_GM * sc.lb, dt = D / t;
  const beta = dt < 15 ? 0.5 : dt <= 60 ? (60 - dt) / 90 : 0, ac = 1 - beta + (beta * fu) / fy, Mp = fy * (D - t) ** 2 * t, Sp = fy * Math.PI * (D - t) * t;
  const lin = (g * Math.abs(+o.M || 0)) / (ac * Mp) + ((g * (+o.S || 0)) / (ac * Sp)) ** 2, dp = (+o.pi || 0) - (+o.pe || 0);
  if (dp >= 0) {
    const pb = ((2 * t) / (D - t)) * Math.min(fy, fu / 1.15) * (2 / Math.sqrt(3)), q = dp / pb, ap = q < 2 / 3 ? 1 - beta : 1 - 3 * beta * (1 - q);
    return { util: lin * lin + ((ap * dp) / (ac * pb)) ** 2, mode: 'internal overpressure', Mp, Sp, alphaC: ac, alphaP: ap, pb };
  }
  const pc = collapsePressure({ D, t, E: o.E, nu: o.nu, fy, ovality: o.ovality }).pc, pex = (+o.pe || 0) - nz(o.pmin, 0);
  return { util: lin * lin + ((g * Math.max(pex, 0)) / pc) ** 2, mode: 'external overpressure', Mp, Sp, alphaC: ac, alphaP: 1, pc };
}
const HOBBS = [['Lateral mode 1', 80.76, 6.391e-5, 0.5], ['Lateral mode 2', 4 * Math.PI ** 2, 1.743e-4, 1], ['Lateral mode 3', 34.06, 1.668e-4, 1.294], ['Lateral mode 4', 28.2, 2.144e-4, 1.608]];
/**
 * Hobbs (1984) thermal-buckling forces of a pipe on a rigid seabed. o: { EI, EA, w (submerged weight N/m), muA (axial friction), muL (lateral friction) }
 * Returns { modes: [{ name, force (minimum fully-restrained axial force that sustains the buckle, N), length (m) }], critical (N), governing }.
 */
export function hobbs(o) {
  const EI = +o.EI, EA = +o.EA, w = Math.max(+o.w, 1e-6), fa = Math.max(nz(o.muA, 0.5), 0.01), fl = Math.max(nz(o.muL, 0.5), 0.01), Ls = logspace(Math.max((EI / w) ** (1 / 3) * 0.05, 1), Math.max((EI / w) ** (1 / 3) * 60, 50), 600);
  const best = (f) => { let b = { force: Infinity, length: Ls[0] }; for (const L of Ls) { const v = f(L); if (Number.isFinite(v) && v < b.force) b = { force: v, length: L }; } return b; };
  const modes = [{ name: 'Upheaval (vertical)', ...best((L) => { const r = 1.597e-5 * EA * fa * w * L ** 5 - 0.25 * (fa * EI) ** 2; return r > 0 ? (80.76 * EI) / (L * L) + ((w * L) / EI) * Math.sqrt(r) : NaN; }) }];
  for (const [name, k1, k2, k3] of HOBBS) modes.push({ name, ...best((L) => (k1 * EI) / (L * L) + k3 * fa * w * L * (Math.sqrt(1 + (k2 * EA * fl * fl * w * L ** 5) / (fa * EI * EI)) - 1)) });
  modes.push({ name: 'Lateral mode ∞', ...best((L) => (4 * Math.PI ** 2 * EI) / (L * L) + 4.705e-5 * EA * ((fl * w) / EI) ** 2 * L ** 6) });
  const lat = modes.slice(1).reduce((a, b) => (b.force < a.force ? b : a));
  return { modes, upheaval: modes[0].force, lateral: lat.force, critical: Math.min(modes[0].force, lat.force), governing: modes[0].force < lat.force ? modes[0].name : lat.name };
}
/**
 * Download needed to hold a pipe down on a seabed imperfection (Palmer et al., 1990). o: { EI, P (compressive effective force, N), delta (imperfection height, m), w0 (weight that shaped the imperfection, N/m) }
 * Returns { wReq (N/m), phiL, phiW, L (half-length of the imperfection, m) }.
 */
export function upheavalDownload(o) {
  const EI = +o.EI, P = Math.max(+o.P, 0), d = Math.max(+o.delta, 0), L = ((72 * EI * Math.max(d, 1e-9)) / Math.max(+o.w0, 1e-6)) ** 0.25;
  if (!(P > 0) || !(d > 0)) return { wReq: 0, phiL: 0, phiW: 0, L };
  const pl = L * Math.sqrt(P / EI), pw = pl < 4.49 ? 0.0646 : pl < 8.06 ? 5.68 / pl ** 2 - 88.35 / pl ** 4 : 9.6 / pl ** 2 - 343 / pl ** 4;
  return { wReq: (pw * d * P * P) / EI, phiL: pl, phiW: pw, L };
}
/** Euler buckling load π²EI/(K·L)² (N). */
export const eulerLoad = (EI, L, K = 1) => (Math.PI ** 2 * EI) / (K * L) ** 2;

// ---- beam finite elements -----------------------------------------------------------------------------------------------
const END_K = { 'pinned-pinned': 1, 'fixed-fixed': 0.5, 'fixed-pinned': 0.699, 'fixed-free': 2, springs: 1 };
const endsOf = (kind) => ({ 'pinned-pinned': ['pinned', 'pinned'], 'fixed-fixed': ['fixed', 'fixed'], 'fixed-pinned': ['fixed', 'pinned'], 'fixed-free': ['fixed', 'free'], springs: ['spring', 'spring'] })[kind] || ['pinned', 'pinned'];
/**
 * Beam model with two-node Hermite elements (deflection and rotation at every node).
 * o: { L, EI, m (kg/m incl. content and added mass), n (elements), ends: ['pinned' | 'fixed' | 'free' | 'spring', …], kT: [N/m, N/m], kR: [N·m/rad, N·m/rad] (end springs),
 *      supports: [x | { x, k }] (intermediate simple supports, k = spring stiffness or rigid), nodesAt: [x] (extra nodes), N (axial force, tension +),
 *      kGA (shear stiffness κGA; > 0 switches to the shear-deformable Timoshenko element), rhoI (rotary inertia per length, kg·m) }
 * Returns { L, EI, m, x[] (nodes), ndof, K, M (full matrices), free[] (unconstrained degrees of freedom), fixed[], springs: [{ dof, k }] }.
 */
export function beamModel(o) {
  const L = +o.L, EI = +o.EI, m = +o.m, n = clamp(Math.round(nz(o.n, 24)), 1, 400), ends = o.ends || ['pinned', 'pinned'], N = nz(o.N, 0), kGA = nz(o.kGA, 0), rhoI = nz(o.rhoI, 0);
  need(L > 0 && EI > 0 && m > 0, 'The beam needs a positive length, bending stiffness and mass per length.');
  const sup = (o.supports || []).map((s) => (typeof s === 'number' ? { x: s } : s)).filter((s) => s && s.x > 1e-6 * L && s.x < L * (1 - 1e-6));
  const bp = [0, L, ...sup.map((s) => +s.x), ...(o.nodesAt || []).filter((x) => x > 0 && x < L)].sort((a, b) => a - b).filter((x, i, a) => i === 0 || x - a[i - 1] > 1e-6 * L);
  const x = [0];
  for (let s = 1; s < bp.length; s++) { const a = bp[s - 1], b = bp[s], ne = Math.max(1, Math.round((n * (b - a)) / L)); for (let i = 1; i <= ne; i++) x.push(a + ((b - a) * i) / ne); }
  const nd = x.length, ndof = 2 * nd, K = zeros(ndof), M = zeros(ndof), elems = [];
  for (let e = 0; e < nd - 1; e++) {
    const l = x[e + 1] - x[e], ph = kGA > 0 ? (12 * EI) / (kGA * l * l) : 0, c = EI / ((1 + ph) * l ** 3), g = N / (30 * l), mm = (m * l) / 420, rr = rhoI / (30 * l);
    const ke = [[12, 6 * l, -12, 6 * l], [6 * l, (4 + ph) * l * l, -6 * l, (2 - ph) * l * l], [-12, -6 * l, 12, -6 * l], [6 * l, (2 - ph) * l * l, -6 * l, (4 + ph) * l * l]];
    const kg = [[36, 3 * l, -36, 3 * l], [3 * l, 4 * l * l, -3 * l, -l * l], [-36, -3 * l, 36, -3 * l], [3 * l, -l * l, -3 * l, 4 * l * l]];
    const me = [[156, 22 * l, 54, -13 * l], [22 * l, 4 * l * l, 13 * l, -3 * l * l], [54, 13 * l, 156, -22 * l], [-13 * l, -3 * l * l, -22 * l, 4 * l * l]];
    for (let a = 0; a < 4; a++) for (let b = 0; b < 4; b++) { K[2 * e + a][2 * e + b] += c * ke[a][b] + g * kg[a][b]; M[2 * e + a][2 * e + b] += mm * me[a][b] + rr * kg[a][b]; }
    elems.push(ke.map((r) => r.map((v) => c * v)));
  }
  const fixedSet = new Set(), springs = [];
  [0, nd - 1].forEach((node, i) => {
    const kind = ends[i] || 'pinned';
    if (kind === 'pinned') fixedSet.add(2 * node); else if (kind === 'fixed') { fixedSet.add(2 * node); fixedSet.add(2 * node + 1); }
    else if (kind === 'spring') { const kt = nz(o.kT?.[i], 0), kr = nz(o.kR?.[i], 0); if (kt > 0) springs.push({ dof: 2 * node, k: kt }); if (kr > 0) springs.push({ dof: 2 * node + 1, k: kr }); }
  });
  for (const s of sup) { let node = 0; for (let i = 1; i < nd; i++) if (Math.abs(x[i] - s.x) < Math.abs(x[node] - s.x)) node = i; if (s.k > 0) springs.push({ dof: 2 * node, k: +s.k }); else fixedSet.add(2 * node); }
  for (const s of springs) K[s.dof][s.dof] += s.k;
  const fixed = [...fixedSet].sort((a, b) => a - b), free = []; for (let i = 0; i < ndof; i++) if (!fixedSet.has(i)) free.push(i);
  return { L, EI, m, x, ndof, K, M, free, fixed, springs, elems, timoshenko: kGA > 0 };
}
const sub = (A, idx) => idx.map((i) => { const r = new Float64Array(idx.length); for (let j = 0; j < idx.length; j++) r[j] = A[i][idx[j]]; return r; });
/** Nodal curvature of a degree-of-freedom vector on the Hermite mesh (averaged at interior nodes). */
function curvature(x, u) {
  const nd = x.length, k = new Float64Array(nd), c = new Float64Array(nd);
  for (let e = 0; e < nd - 1; e++) { const l = x[e + 1] - x[e], w1 = u[2 * e], t1 = u[2 * e + 1], w2 = u[2 * e + 2], t2 = u[2 * e + 3]; k[e] += (-6 * w1 - 4 * l * t1 + 6 * w2 - 2 * l * t2) / (l * l); c[e]++; k[e + 1] += (6 * w1 + 2 * l * t1 - 6 * w2 + 4 * l * t2) / (l * l); c[e + 1]++; }
  for (let i = 0; i < nd; i++) k[i] /= c[i];
  return k;
}
/**
 * Natural frequencies and mode shapes of a beam model (generalised eigenproblem by subspace iteration with Jacobi rotations).
 * Returns { f[] (Hz), omega[] (rad/s), shapes[k][node] (deflection, largest value 1), vectors[k] (mass-normalised, all degrees of freedom),
 *   curv[k][node] (curvature of the mass-normalised mode, 1/m per unit modal coordinate), unitCurv[k] (largest curvature per metre of modal amplitude) }.
 */
export function beamModes(model, nModes = 4) {
  const fr = model.free, p = Math.min(Math.max(1, Math.round(nModes)), fr.length), e = subspaceEig(sub(model.K, fr), sub(model.M, fr), p), out = { f: [], omega: [], shapes: [], vectors: [], curv: [], unitCurv: [] };
  for (let k = 0; k < p; k++) {
    need(e.values[k] > (1e-7 * model.EI) / (model.m * model.L ** 4), 'The span is unstable (a natural frequency is zero): it is a mechanism without enough supports, or it has buckled under the compressive axial force.');
    const om = Math.sqrt(e.values[k]), full = new Float64Array(model.ndof); fr.forEach((d, i) => (full[d] = e.vectors[k][i]));
    let big = 0; for (let i = 0; i < model.x.length; i++) if (Math.abs(full[2 * i]) > Math.abs(big)) big = full[2 * i];
    if (big < 0) for (let i = 0; i < full.length; i++) full[i] = -full[i];
    const cv = curvature(model.x, full), amp = Math.abs(big) || 1;
    out.omega.push(om); out.f.push(om / (2 * Math.PI)); out.vectors.push(full); out.curv.push(cv); out.shapes.push(model.x.map((_, i) => full[2 * i] / amp)); out.unitCurv.push(Math.max(...Array.from(cv, Math.abs)) / amp);
  }
  return out;
}
/**
 * Static solution of a beam model. load: { q (uniform N/m in the direction of positive deflection), points: [{ x, F }] (x must be nodes: pass them as nodesAt) }
 * Returns { w[] (m), theta[] (rad), moment[] (N·m at nodes, from element end forces), reactions: [{ x, R (N) }], sumReactions, totalLoad, maxMoment, maxDeflection }.
 */
export function beamStatic(model, load = {}) {
  const { x, ndof, K, free: fr, fixed } = model, nd = x.length, f = new Float64Array(ndof), q = nz(load.q, 0), elem = [];
  for (let e = 0; e < nd - 1; e++) { const l = x[e + 1] - x[e], fe = [(q * l) / 2, (q * l * l) / 12, (q * l) / 2, (-q * l * l) / 12]; elem.push(fe); for (let a = 0; a < 4; a++) f[2 * e + a] += fe[a]; }
  let total = q * model.L;
  for (const p of load.points || []) { let node = 0; for (let i = 1; i < nd; i++) if (Math.abs(x[i] - p.x) < Math.abs(x[node] - p.x)) node = i; f[2 * node] += p.F; total += p.F; }
  const Lc = chol(sub(K, fr)); need(Lc, 'The beam is a mechanism or has buckled: add supports or reduce the compressive axial force.');
  const ur = cholSolve(Lc, fr.map((d) => f[d])), u = new Float64Array(ndof); fr.forEach((d, i) => (u[d] = ur[i]));
  const Ku = matVec(K, u), reactions = [];
  for (const d of fixed) if (d % 2 === 0) reactions.push({ x: x[d / 2], R: Ku[d] - f[d] });
  for (const s of model.springs) if (s.dof % 2 === 0) reactions.push({ x: x[s.dof / 2], R: -s.k * u[s.dof] });
  // nodal bending moments M = EI·w″ from the elastic element end forces (exact at the nodes for consistent loads)
  const moment = new Float64Array(nd), cnt = new Float64Array(nd);
  for (let e = 0; e < nd - 1; e++) {
    const ke = model.elems[e]; let r1 = -elem[e][1], r3 = -elem[e][3];
    for (let a = 0; a < 4; a++) { r1 += ke[1][a] * u[2 * e + a]; r3 += ke[3][a] * u[2 * e + a]; }
    moment[e] -= r1; cnt[e]++; moment[e + 1] += r3; cnt[e + 1]++;
  }
  for (let i = 0; i < nd; i++) moment[i] /= cnt[i];
  const w = x.map((_, i) => u[2 * i]);
  return { w, theta: x.map((_, i) => u[2 * i + 1]), moment: Array.from(moment), reactions, sumReactions: sum(reactions.map((r) => r.R)), totalLoad: total, maxMoment: Math.max(...Array.from(moment, Math.abs)), maxDeflection: Math.max(...w.map(Math.abs)), u };
}
/** LU factorisation with partial pivoting of a small dense matrix; returns a solver b -> x. */
function luSolver(A) {
  const n = A.length, a = A.map((r) => Float64Array.from(r)), piv = Array.from({ length: n }, (_, i) => i);
  for (let k = 0; k < n; k++) {
    let p = k; for (let i = k + 1; i < n; i++) if (Math.abs(a[i][k]) > Math.abs(a[p][k])) p = i;
    need(Math.abs(a[p][k]) > 1e-300, 'Singular matrix in the dynamic solver.');
    if (p !== k) { [a[k], a[p]] = [a[p], a[k]]; [piv[k], piv[p]] = [piv[p], piv[k]]; }
    for (let i = k + 1; i < n; i++) { const f = (a[i][k] /= a[k][k]); if (f !== 0) for (let j = k + 1; j < n; j++) a[i][j] -= f * a[k][j]; }
  }
  return (b) => {
    const x = new Float64Array(n);
    for (let i = 0; i < n; i++) { let s = b[piv[i]]; for (let j = 0; j < i; j++) s -= a[i][j] * x[j]; x[i] = s; }
    for (let i = n - 1; i >= 0; i--) { let s = x[i]; for (let j = i + 1; j < n; j++) s -= a[i][j] * x[j]; x[i] = s / a[i][i]; }
    return x;
  };
}
/**
 * Newmark-β time integration of M ü + C u̇ + K u = f(t) (average acceleration by default: unconditionally stable, no numerical damping).
 * o: { M, C (optional), K (dense), f: (t) => array, u0, v0, dt, steps, beta = 1/4, gamma = 1/2, onStep(i, t, u, v) }
 * Returns { t[], u[] (final), v[] (final), energy[] (kinetic + strain at every step), history[] (u of the tracked degrees of freedom when o.track is given) }.
 */
export function newmark(o) {
  const { M, K, dt } = o, n = M.length, be = o.beta ?? 0.25, ga = o.gamma ?? 0.5, C = o.C || zeros(n), steps = Math.max(1, Math.round(o.steps));
  need(dt > 0 && n > 0, 'The dynamic solver needs a positive time step.');
  const a0 = 1 / (be * dt * dt), a1 = ga / (be * dt), a2 = 1 / (be * dt), a3 = 1 / (2 * be) - 1, a4 = ga / be - 1, a5 = (dt / 2) * (ga / be - 2);
  const Keff = zeros(n); for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) Keff[i][j] = K[i][j] + a0 * M[i][j] + a1 * C[i][j];
  const solveK = luSolver(Keff), solveM = luSolver(M);
  let u = Float64Array.from(o.u0 || new Array(n).fill(0)), v = Float64Array.from(o.v0 || new Array(n).fill(0));
  const f0 = o.f(0), Ku = matVec(K, u), Cv = matVec(C, v);
  let a = solveM(Float64Array.from({ length: n }, (_, i) => f0[i] - Ku[i] - Cv[i]));
  const energy = (uu, vv) => 0.5 * dot(vv, matVec(M, vv)) + 0.5 * dot(uu, matVec(K, uu));
  const t = [0], en = o.energy === false ? null : [energy(u, v)], hist = o.track ? [o.track.map((d) => u[d])] : null, rhs = new Float64Array(n), tm = new Float64Array(n), tc = new Float64Array(n);
  for (let s = 1; s <= steps; s++) {
    const ts = s * dt, f = o.f(ts);
    for (let i = 0; i < n; i++) { tm[i] = a0 * u[i] + a2 * v[i] + a3 * a[i]; tc[i] = a1 * u[i] + a4 * v[i] + a5 * a[i]; }
    const Mt = matVec(M, tm), Ct = matVec(C, tc);
    for (let i = 0; i < n; i++) rhs[i] = f[i] + Mt[i] + Ct[i];
    const un = solveK(rhs), an = new Float64Array(n), vn = new Float64Array(n);
    for (let i = 0; i < n; i++) { an[i] = a0 * (un[i] - u[i]) - a2 * v[i] - a3 * a[i]; vn[i] = v[i] + dt * ((1 - ga) * a[i] + ga * an[i]); }
    u = un; v = vn; a = an; t.push(ts);
    if (en) en.push(energy(u, v));
    if (hist) hist.push(o.track.map((d) => u[d]));
    if (o.onStep) o.onStep(s, ts, u, v);
  }
  return { t, u: Array.from(u), v: Array.from(v), energy: en, history: hist };
}

// ---- slug and flow-induced loads ----------------------------------------------------------------------------------------
/**
 * Momentum force of a fluid stream turned by a bend. o: { rho (kg/m³), A (m²), v (m/s), angle (deg), dlf (dynamic load factor), p (gauge Pa, optional), radius (m, optional) }
 * Returns { force (N, resultant = DLF·ρ·A·v²·√(2(1 − cos θ))), fx, fy (components of the static momentum force), pressureForce (N), centrifugal (N/m along the bend), turning (2 sin(θ/2)) }.
 */
export function bendForce(o) {
  const th = (nz(o.angle, 90) * Math.PI) / 180, k = Math.sqrt(2 * (1 - Math.cos(th))), mom = +o.rho * +o.A * o.v * o.v, dlf = nz(o.dlf, 1);
  return { force: dlf * mom * k, fx: mom * (1 - Math.cos(th)), fy: mom * Math.sin(th), pressureForce: nz(o.p, 0) * o.A * k, centrifugal: o.radius > 0 ? mom / o.radius : 0, turning: k, momentumFlux: mom };
}
/**
 * Energy-Institute-style screening of flow-induced turbulence in piping: kinetic energy ρv² and likelihood of failure.
 * o: { rho, v, Dmm (outer diameter), tmm (wall), support: 'stiff' | 'mediumStiff' | 'medium', fvf (fluid viscosity factor, 1 for multiphase) }
 */
export function fivScreen(o) {
  const ke = +o.rho * o.v * o.v, D = +o.Dmm, lnD = Math.log(D);
  const [al, be] = o.support === 'stiff' ? [446187 + 646 * D + 9.17e-4 * D ** 3, 0.1 * lnD - 1.3739] : o.support === 'medium' ? [150412 + 209 * D, 0.0815 * lnD - 1.3269] : [283921 + 370 * D, 0.1106 * lnD - 1.501];
  const Fv = al * (D / +o.tmm) ** be, lof = (ke / Fv) * nz(o.fvf, 1);
  return { rhoV2: ke, Fv, lof, band: ke < 5000 ? 'low' : ke < 20000 ? 'medium' : 'high', likelihood: lof < 0.3 ? 'low' : lof < 0.5 ? 'medium' : lof < 1 ? 'medium-high' : 'high' };
}
/**
 * Vortex-induced-vibration screening of a free span in a steady current (DNV-RP-F105-type response models, screening level).
 * o: { U (m/s), f1 (Hz, in-line ≈ cross-flow), D (outer incl. coating), me (effective mass kg/m), zeta (total damping ratio), rhoW, gammaOn: { il, cf }, gammaK }
 * Returns { vr, ks, ksd, onsetIL, onsetCF, fShed (Strouhal 0.2), aIL, aCF (amplitude / D), state }.
 */
export function vivScreen(o) {
  const U = Math.max(+o.U, 0), D = +o.D, f1 = +o.f1, vr = U / (f1 * D), ks = (4 * Math.PI * +o.me * +o.zeta) / (nz(o.rhoW, RHO_SW) * D * D), ksd = ks / nz(o.gammaK, 1.15), gIL = o.gammaOn?.il ?? 1.1, gCF = o.gammaOn?.cf ?? 1.2;
  const onIL = (ksd < 0.4 ? 1 : ksd < 1.6 ? 0.6 + ksd : 2.2) / gIL, onCF = 3 / gCF;
  // in-line response model
  const a2 = 0.13 * Math.max(1 - ksd / 1.8, 0), a1 = Math.max(0.18 * Math.max(1 - ksd / 1.2, 0), a2), v1 = 10 * a1 + onIL, vEnd = ksd < 1 ? 4.5 - 0.8 * ksd : 3.7, v2 = vEnd - 2 * a2;
  const aIL = a1 <= 0 || vr <= onIL || vr >= vEnd ? 0 : interp1([onIL, Math.max(v1, onIL + 1e-6), Math.max(v2, v1 + 2e-6, onIL + 2e-6), Math.max(vEnd, v2 + 3e-6)], [0, a1, a2, 0], vr);
  // cross-flow response model, reduced by damping
  const rk = ksd <= 4 ? 1 - 0.15 * ksd : 3.2 * ksd ** -1.5, az = 0.9, c1 = 7 - ((7 - onCF) / 1.15) * (1.3 - az), c2 = 16 - (7 / 1.3) * az;
  const aCF = vr <= 2 || vr >= 16 ? 0 : rk * interp1([2, onCF, c1, c2, 16], [0, 0.15, az, az, 0], vr);
  return { vr, ks, ksd, onsetIL: onIL, onsetCF: onCF, fShed: (0.2 * U) / D, aIL, aCF, state: vr >= onCF ? 'cross-flow VIV' : vr >= onIL ? 'in-line VIV' : 'no VIV' };
}
/**
 * Response of a span to a train of slugs by modal superposition (one-way coupling: the fluid loads the pipe, the pipe does not act back on the flow).
 * Each modal equation is integrated with Newmark-β. The slug is a moving change of distributed weight `dw` (N/m) over its body length and, optionally,
 * a momentum force on a bend at `bend.x` that is `bend.dF` larger while the slug body is in the bend.
 * o: { model, modes, zeta, slugs: [{ t0 (s, front enters), len (m), v (m/s) }], dw, bend: { x, dF } | null, dt, tEnd, ro (outer steel radius), E }
 * Returns { t[], sigma[] (bending stress history at the critical node, Pa), disp[] (deflection at that node, m), node, x, sigmaMax, dispMax, steps }.
 */
export function slugResponse(o) {
  const { model, modes } = o, p = modes.f.length, nd = model.x.length, x = model.x, L = model.L, dt = o.dt, steps = Math.max(2, Math.round(o.tEnd / dt)), zeta = nz(o.zeta, 0.02);
  // cumulative integral of every mode along the span (exact for the Hermite interpolation at the nodes)
  const cum = modes.vectors.map((u) => { const c = new Float64Array(nd); for (let e = 0; e < nd - 1; e++) { const l = x[e + 1] - x[e]; c[e + 1] = c[e] + (l / 2) * (u[2 * e] + u[2 * e + 2]) + ((l * l) / 12) * (u[2 * e + 1] - u[2 * e + 3]); } return c; });
  const at = (arr, s) => { if (s <= 0) return arr[0]; if (s >= L) return arr[nd - 1]; let lo = 0, hi = nd - 1; while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (x[mid] <= s) lo = mid; else hi = mid; } return arr[lo] + ((arr[hi] - arr[lo]) * (s - x[lo])) / (x[hi] - x[lo]); };
  const disp = modes.vectors.map((u) => Float64Array.from({ length: nd }, (_, i) => u[2 * i])), bend = o.bend && o.bend.dF ? o.bend : null, phiB = bend ? disp.map((d) => at(d, bend.x)) : null, ramp = Math.max(nz(o.ramp, 0.3), 1e-3);
  const slugs = o.slugs, dw = nz(o.dw, 0); let first = 0;
  const force = (t) => {
    const Q = new Float64Array(p);
    while (first < slugs.length && (t - slugs[first].t0) * slugs[first].v - slugs[first].len > L + ramp) first++;
    for (let i = first; i < slugs.length; i++) {
      const s = slugs[i]; if (t < s.t0) break;
      const front = (t - s.t0) * s.v, tail = front - s.len, a = Math.max(tail, 0), b = Math.min(front, L);
      if (b > a && dw) for (let k = 0; k < p; k++) Q[k] += dw * (at(cum[k], b) - at(cum[k], a));
      if (bend) { const c = clamp((front - bend.x) / ramp, 0, 1) - clamp((tail - bend.x) / ramp, 0, 1); if (c > 0) for (let k = 0; k < p; k++) Q[k] += bend.dF * c * phiB[k]; }
    }
    return Q;
  };
  // Newmark-β (average acceleration) on the uncoupled modal equations q̈ + 2ζω q̇ + ω² q = Q(t), unit modal mass
  const q = new Float64Array(p * (steps + 1)), qv = new Float64Array(p), qa = Float64Array.from(force(0)), a0 = 4 / (dt * dt), a1 = 2 / dt, a2 = 4 / dt;
  for (let i = 1; i <= steps; i++) {
    const Q = force(i * dt), o0 = (i - 1) * p, o1 = i * p;
    for (let k = 0; k < p; k++) { const w2 = modes.omega[k] ** 2, c = 2 * zeta * modes.omega[k], u = q[o0 + k], un = (Q[k] + a0 * u + a2 * qv[k] + qa[k] + c * (a1 * u + qv[k])) / (w2 + a0 + a1 * c), an = a0 * (un - u) - a2 * qv[k] - qa[k]; qv[k] += 0.5 * dt * (qa[k] + an); qa[k] = an; q[o1 + k] = un; }
  }
  // envelope of the bending curvature at every node, then the history at the worst node
  let node = 0, worst = -1;
  const stride = Math.max(1, Math.floor(steps / 4000));
  for (let j = 0; j < nd; j++) { let lo = 0, hi = 0; for (let i = 0; i <= steps; i += stride) { let c = 0; for (let k = 0; k < p; k++) c += q[i * p + k] * modes.curv[k][j]; if (c < lo) lo = c; if (c > hi) hi = c; } if (hi - lo > worst) { worst = hi - lo; node = j; } }
  const fac = +o.E * +o.ro, t = new Array(steps + 1), sg = new Array(steps + 1), ds = new Array(steps + 1); let sMax = 0, dMax = 0;
  for (let i = 0; i <= steps; i++) { let c = 0, d = 0; for (let k = 0; k < p; k++) { c += q[i * p + k] * modes.curv[k][node]; d += q[i * p + k] * disp[k][node]; } t[i] = i * dt; sg[i] = fac * c; ds[i] = d; if (Math.abs(sg[i]) > sMax) sMax = Math.abs(sg[i]); }
  for (let j = 0; j < nd; j += Math.max(1, Math.floor(nd / 40))) for (let i = 0; i <= steps; i += stride) { let d = 0; for (let k = 0; k < p; k++) d += q[i * p + k] * disp[k][j]; if (Math.abs(d) > dMax) dMax = Math.abs(d); }
  return { t, sigma: sg, disp: ds, node, x: x[node], sigmaMax: sMax, dispMax: dMax, steps };
}

// ---- fatigue and fracture -----------------------------------------------------------------------------------------------
/**
 * Rainflow cycle counting of a load history (ASTM E1049, three-point method with half cycles for the residue).
 * Returns [{ range, mean, count (1 or 0.5) }].
 */
export function rainflow(series) {
  const rev = [];
  for (const v of series) { // turning points
    if (!Number.isFinite(v)) continue;
    const n = rev.length;
    if (n && v === rev[n - 1]) continue;
    if (n >= 2 && (rev[n - 1] - rev[n - 2]) * (v - rev[n - 1]) > 0) rev[n - 1] = v; else rev.push(v);
  }
  const out = [], st = [];
  for (const p of rev) {
    st.push(p);
    while (st.length >= 3) {
      const n = st.length, X = Math.abs(st[n - 1] - st[n - 2]), Y = Math.abs(st[n - 2] - st[n - 3]);
      if (X < Y) break;
      if (n === 3) { out.push({ range: Y, mean: (st[0] + st[1]) / 2, count: 0.5 }); st.shift(); } // the range contains the starting point
      else { out.push({ range: Y, mean: (st[n - 2] + st[n - 3]) / 2, count: 1 }); st.splice(n - 3, 2); }
    }
  }
  for (let i = 0; i + 1 < st.length; i++) out.push({ range: Math.abs(st[i + 1] - st[i]), mean: (st[i + 1] + st[i]) / 2, count: 0.5 });
  return out;
}
// DNV-RP-C203 S–N curves in air: N = 10^(a − m·log S), slope m1 up to 10⁷ cycles, m2 beyond; k = thickness exponent (reference thickness 25 mm).
export const SN_CURVES = Object.freeze({
  B1: { m1: 4, a1: 15.117, m2: 5, a2: 17.146, k: 0 }, B2: { m1: 4, a1: 14.885, m2: 5, a2: 16.856, k: 0 }, C: { m1: 3, a1: 12.592, m2: 5, a2: 16.32, k: 0.05 }, C1: { m1: 3, a1: 12.449, m2: 5, a2: 16.081, k: 0.1 },
  C2: { m1: 3, a1: 12.301, m2: 5, a2: 15.835, k: 0.15 }, D: { m1: 3, a1: 12.164, m2: 5, a2: 15.606, k: 0.2 }, E: { m1: 3, a1: 12.01, m2: 5, a2: 15.35, k: 0.2 }, F: { m1: 3, a1: 11.855, m2: 5, a2: 15.091, k: 0.25 },
  F1: { m1: 3, a1: 11.699, m2: 5, a2: 14.832, k: 0.25 }, F3: { m1: 3, a1: 11.546, m2: 5, a2: 14.576, k: 0.25 }, G: { m1: 3, a1: 11.398, m2: 5, a2: 14.33, k: 0.25 }, W1: { m1: 3, a1: 11.261, m2: 5, a2: 14.101, k: 0.25 },
  W2: { m1: 3, a1: 11.107, m2: 5, a2: 13.845, k: 0.25 }, W3: { m1: 3, a1: 10.97, m2: 5, a2: 13.617, k: 0.25 },
});
/**
 * Cycles to failure for a stress range S (MPa). env: 'air' | 'cp' (seawater with cathodic protection: a1 − 0.4, knee at 10⁶) | 'free' (free corrosion: single slope, a1 − 0.477, m = 3).
 * opt: { t (mm wall, thickness correction above 25 mm), scf }.
 */
export function snCycles(S, cls = 'D', env = 'air', opt = {}) {
  const c = SN_CURVES[cls] || SN_CURVES.D, t = nz(opt.t, 25), s = Math.abs(S) * nz(opt.scf, 1) * (t > 25 ? (t / 25) ** c.k : 1);
  if (!(s > 0)) return Infinity;
  const ls = Math.log10(s);
  if (env === 'free') return 10 ** ((c.m1 === 4 ? c.a1 - 2.681 + (cls === 'B2' ? 0.058 : 0) : c.a1 - 0.477) - 3 * ls);
  const a1 = env === 'cp' ? c.a1 - 0.4 : c.a1, knee = env === 'cp' ? 1e6 : 1e7, n1 = 10 ** (a1 - c.m1 * ls);
  return n1 <= knee ? n1 : 10 ** (c.a2 - c.m2 * ls);
}
/** Palmgren–Miner damage of counted cycles [{ range (MPa), count }]: returns { damage, cycles, sEq (equivalent constant range for slope 3, MPa) }. */
export function snDamage(cycles, cls = 'D', env = 'air', opt = {}) {
  let d = 0, n = 0, s3 = 0;
  for (const c of cycles) { if (!(c.range > 0) || !(c.count > 0)) continue; const N = snCycles(c.range, cls, env, opt); if (Number.isFinite(N)) d += c.count / N; n += c.count; s3 += c.count * c.range ** 3; }
  return { damage: d, cycles: n, sEq: n > 0 ? (s3 / n) ** (1 / 3) : 0 };
}
/** Geometry factor of an edge crack of depth a in a plate of thickness t under tension (Tada), valid to a/t ≈ 0.6. */
export const edgeCrackY = (a, t) => { const r = clamp(a / t, 0, 0.7); return 1.12 - 0.231 * r + 10.55 * r * r - 21.72 * r ** 3 + 30.39 * r ** 4; };
/** Closed-form Paris–Erdogan life for a constant geometry factor: cycles from a0 to ac (m, MPa, MPa√m units for C). */
export function parisClosedForm({ a0, ac, C, m, dS, Y = 1.12 }) {
  const k = C * (Y * dS * Math.sqrt(Math.PI)) ** m;
  return Math.abs(m - 2) < 1e-12 ? Math.log(ac / a0) / k : (ac ** (1 - m / 2) - a0 ** (1 - m / 2)) / (k * (1 - m / 2));
}
/**
 * Paris–Erdogan crack growth da/dN = C·ΔK^m with ΔK = Y·Δσ·√(πa), integrated numerically from a0 to ac.
 * o: { a0, ac (m), C (m/cycle with ΔK in MPa√m), m, dS (MPa), Y (number or function of a), dKth (MPa√m threshold), n (integration points) }
 * Returns { N (cycles; Infinity when the initial flaw is below the threshold), a[] (m), cycles[], dK0 }.
 */
export function parisLife(o) {
  const a0 = +o.a0, ac = +o.ac, Y = typeof o.Y === 'function' ? o.Y : () => nz(o.Y, 1.12), n = Math.max(20, Math.round(nz(o.n, 600)));
  need(a0 > 0 && o.C > 0 && o.m > 0, 'Crack growth needs a positive initial flaw and Paris constants.');
  const dK = (a) => Y(a) * o.dS * Math.sqrt(Math.PI * a), dK0 = dK(a0);
  if (!(ac > a0)) return { N: 0, a: [a0], cycles: [0], dK0 };
  if (!(dK0 > nz(o.dKth, 0)) || !(o.dS > 0)) return { N: Infinity, a: [a0, a0], cycles: [0, 1e12], dK0 };
  const a = logspace(a0, ac, 2 * Math.floor(n / 2) + 1), f = a.map((x) => 1 / (o.C * dK(x) ** o.m)), cyc = [0];
  for (let i = 2; i < a.length; i += 2) { const h1 = a[i - 1] - a[i - 2], h2 = a[i] - a[i - 1], s = ((h1 + h2) / 6) * ((2 - h2 / h1) * f[i - 2] + ((h1 + h2) ** 2 / (h1 * h2)) * f[i - 1] + (2 - h1 / h2) * f[i]); cyc.push(cyc[cyc.length - 1] + s); } // Simpson on the non-uniform grid
  return { N: cyc[cyc.length - 1], a: a.filter((_, i) => i % 2 === 0), cycles: cyc, dK0 };
}

// ---- corrosion ----------------------------------------------------------------------------------------------------------
/** Faraday's law: corrosion current density (A/m²) to penetration rate (mm/y) for iron (M = 55.845 g/mol, n = 2, ρ = 7870 kg/m³). */
export const faradayRate = (i, M = 55.845e-3, n = 2, rho = 7870) => ((i * M) / (n * FARADAY * rho)) * YEAR * 1000;
/** Butler–Volmer current density (A/m²): i0·[exp(αa·F·η/RT) − exp(−αc·F·η/RT)], overpotential η in V, T in K. */
export const butlerVolmer = (eta, i0, T = 298.15, aa = 0.5, ac = 0.5) => i0 * (Math.exp((aa * FARADAY * eta) / (RGAS * T)) - Math.exp((-ac * FARADAY * eta) / (RGAS * T)));
/** Sherwood number of turbulent pipe flow (Berger & Hau): Sh = 0.0165·Re^0.86·Sc^0.33. */
export const sherwood = (Re, Sc) => 0.0165 * Math.max(Re, 1) ** 0.86 * Sc ** 0.33;
/**
 * de Waard–Milliams CO2 corrosion rate (mm/y).
 * o: { T (°C), pCO2 (bar), P (bar total, for the fugacity correction), model: '1991' (nomogram equation) | '1995' (reaction + mass-transfer resistances),
 *      U (liquid velocity m/s), d (hydraulic diameter m), pH (actual; omitted = CO2-saturated water), glycolWt (wt % glycol in the water phase),
 *      inhibEff (0–1), mult, fugacity (false to switch the correction off), scale (false to switch the protective-film factor off) }
 * Returns { rate, base, fCO2, fugacity, pHco2, pH, Fscale, tScale (°C), Fglycol, Vr, Vm }.
 */
export function deWaardMilliams(o) {
  const Tc = nz(o.T, 60), T = Tc + 273.15, p = Math.max(nz(o.pCO2, 0), 0);
  if (!(p > 0)) return { rate: 0, base: 0, fCO2: 0, fugacity: 1, pHco2: 7, pH: nz(o.pH, 7), Fscale: 1, tScale: 0, Fglycol: 1, Vr: 0, Vm: 0 };
  const fug = o.fugacity === false ? 1 : 10 ** (Math.min(Math.max(nz(o.P, p), p), 250) * (0.0031 - 1.4 / T)), f = p * fug, lf = Math.log10(f), pHco2 = 3.71 + 0.00417 * Tc - 0.5 * lf, pH = o.pH > 0 ? +o.pH : pHco2;
  const Fscale = o.scale === false ? 1 : Math.min(1, 10 ** (2400 / T - 0.6 * lf - 6.7)), tScale = 2400 / (6.7 + 0.6 * lf) - 273.15, W = clamp(100 - nz(o.glycolWt, 0), 1, 100), Fglycol = 10 ** (1.6 * (Math.log10(W) - 2));
  let base, Vr = 0, Vm = 0;
  if (o.model === '1991') base = 10 ** (5.8 - 1710 / T + 0.67 * lf);
  else { Vr = 10 ** (4.93 - 1119 / T + 0.58 * lf - 0.34 * (pH - pHco2)); Vm = (2.45 * Math.max(nz(o.U, 1), 0.01) ** 0.8 * f) / pos(o.d, 0.1) ** 0.2; base = 1 / (1 / Vr + 1 / Vm); }
  return { rate: base * Fscale * Fglycol * (1 - clamp(nz(o.inhibEff, 0), 0, 1)) * nz(o.mult, 1), base, fCO2: f, fugacity: fug, pHco2, pH, Fscale, tScale, Fglycol, Vr, Vm };
}
/**
 * Mixed-potential model of iron in CO2-saturated water: anodic iron dissolution against the reduction of H⁺ (mass-transfer limited) and of
 * carbonic acid (limited by the hydration of CO2), Tafel kinetics, solved for the corrosion potential (kinetic constants after Nesic et al., 1996).
 * o: { T (°C), pCO2 (bar, fugacity), pH, U (m/s), d (m), rho, mu (liquid), ionic (mol/L) }
 * Returns { Ecorr (V), icorr (A/m²), rate (mm/y), iLimH, iLimH2CO3 (A/m²), rateLimit (mm/y, mass-transfer / reaction limit), Sh, km (m/s), tafel: { ba, bc } }.
 */
export function mixedPotential(o) {
  const Tc = nz(o.T, 20), T = Tc + 273.15, p = Math.max(nz(o.pCO2, 1), 1e-6), pH = pos(o.pH, 3.71 + 0.00417 * Tc - 0.5 * Math.log10(p)), U = Math.max(nz(o.U, 1), 0.01), d = pos(o.d, 0.1), rho = pos(o.rho, 1000), mu = pos(o.mu, 1e-3);
  const Tf = Tc * 1.8 + 32, ksol = (14.5 / 1.00258) * 10 ** -(2.27 + 5.65e-3 * Tf - 8.06e-6 * Tf * Tf + 0.075 * nz(o.ionic, 0)), cCO2 = ksol * p * 1000, cH = 10 ** -pH * 1000, cH2CO3 = 2.58e-3 * cCO2; // mol/m³
  const dScale = (T / 298.15) * (8.9e-4 / mu), DH = 9.31e-9 * dScale, DC = 1.3e-9 * dScale, Re = (rho * U * d) / mu, Sh = sherwood(Re, mu / (rho * DH)), km = (Sh * DH) / d;
  const iLimH = km * FARADAY * cH, khyd = 10 ** (329.85 - 110.541 * Math.log10(T) - 17265.4 / T), iLimC = FARADAY * cCO2 * Math.sqrt(DC * 2.58e-3 * khyd);
  const arr = (dH, Tref) => Math.exp((-dH / RGAS) * (1 / T - 1 / Tref)), b = (2.303 * RGAS * T) / FARADAY, ba = b / 1.5, bc = b / 0.5, eH = -b * pH, eFe = -0.488;
  const i0H = 0.05 * (cH / 0.1) ** 0.5 * arr(30e3, 298.15), i0C = 0.06 * (cH / 0.01) ** -0.5 * (cH2CO3 / 0.1) * arr(50e3, 293.15), i0Fe = 1 * arr(37.5e3, 298.15);
  const ia = (E) => i0Fe * 10 ** ((E - eFe) / ba), ic = (E) => 1 / (1 / (i0H * 10 ** (-(E - eH) / bc)) + 1 / iLimH) + 1 / (1 / (i0C * 10 ** (-(E - eH) / bc)) + 1 / iLimC);
  const Ecorr = brent((E) => ia(E) - ic(E), -1.5, 0.3, 1e-12), icorr = ia(Ecorr);
  return { Ecorr, icorr, rate: faradayRate(icorr), iLimH, iLimH2CO3: iLimC, rateLimit: faradayRate(iLimH + iLimC), Sh, km, Re, pH, tafel: { ba, bc } };
}
/** Sour-service check in the ISO 15156-2 format: region 0 below 0.3 kPa H2S, regions 1–3 from the in-situ pH (boundaries digitised from the severity diagram). */
export function sourRegion(pH2SkPa, pH) {
  if (!(pH2SkPa >= 0.3)) return { region: 0, sour: false, label: 'Region 0 (below 0.3 kPa H2S: no special requirements)' };
  const up = clamp(Math.log10(pH2SkPa), 0, 2), region = pH >= 4.5 + up ? 1 : pH >= 3.5 + up ? 2 : 3;
  return { region, sour: true, label: `SSC region ${region} (${['mild', 'intermediate', 'severe'][region - 1]} sour service)` };
}
/** Top-of-line condensation rate (g/m²/s) of a water-saturated gas cooling along the pipe: gas rate (Sm³/s) × d(water content)/dT × temperature gradient over the wall area. */
export function condensationRate({ P, T, dTdx, qGasStd, D }) {
  const dW = (waterContent(P, T + 1) - waterContent(P, T - 1)) / 2; // kg/Sm³ per °C
  return Math.max(0, (qGasStd * dW * -nz(dTdx, 0) * 1000) / (Math.PI * D));
}
/**
 * Burst pressure of a pipe with a single metal-loss defect (Pa).
 * o: { D, t, d (depth), L (axial length) (m), smys, smts (Pa), method: 'b31g' | 'modified' | 'dnv' }
 * 'b31g': ASME B31G (flow stress 1.1·SMYS, parabolic area 2/3·d·L, Folias √(1 + 0.8 L²/Dt), rectangular beyond L²/Dt = 20);
 * 'modified': 0.85·d·L area, flow stress SMYS + 69 MPa, two-term Folias factor; 'dnv': DNV-RP-F101 single defect, Q = √(1 + 0.31 (L/√(Dt))²).
 * Returns { pf, intact (Pa), ratio, M (bulging factor), flow (Pa) }.
 */
export function b31g(o) {
  const D = +o.D, t = +o.t, r = clamp(o.d / t, 0, 1), z = (o.L * o.L) / (D * t);
  need(D > 0 && t > 0, 'Defect assessment needs a positive diameter and wall thickness.');
  if (o.method === 'dnv') { const Q = Math.sqrt(1 + 0.31 * z), p0 = (2 * t * o.smts) / (D - t), pf = r >= 1 ? 0 : (p0 * (1 - r)) / (1 - r / Q); return { pf, intact: p0, ratio: pf / p0, M: Q, flow: +o.smts }; }
  if (o.method === 'modified') { const M = z <= 50 ? Math.sqrt(1 + 0.6275 * z - 0.003375 * z * z) : 0.032 * z + 3.3, S = +o.smys + 69e6, p0 = (2 * t * S) / D, pf = (p0 * (1 - 0.85 * r)) / (1 - (0.85 * r) / M); return { pf, intact: p0, ratio: pf / p0, M, flow: S }; }
  const M = Math.sqrt(1 + 0.8 * z), S = 1.1 * o.smys, p0 = (2 * t * S) / D, pf = z <= 20 ? (p0 * (1 - (2 / 3) * r)) / (1 - ((2 / 3) * r) / M) : p0 * (1 - r);
  return { pf, intact: p0, ratio: pf / p0, M, flow: S };
}
/** Statistics of an imported wall map ({ x[], theta[], t[θ][x], kind, unit }) in mm: minimum wall, deepest loss, its position and axial extent. */
export function wallMapStats(map, tNomMm) {
  if (!map || !Array.isArray(map.t) || !Array.isArray(map.x) || !map.t.length) return null;
  const toMM = { m: 1000, cm: 10, in: 25.4, inch: 25.4, um: 1e-3, 'µm': 1e-3 }[map.unit] ?? 1, loss = map.kind === 'corrosion';
  let worst = -Infinity, iw = 0, jw = 0, s = 0, n = 0;
  map.t.forEach((row, i) => (row || []).forEach((v, j) => { if (!Number.isFinite(v)) return; const d = loss ? v * toMM : tNomMm - v * toMM; s += d; n++; if (d > worst) { worst = d; iw = i; jw = j; } }));
  if (!n) return null;
  const depth = clamp(worst, 0, tNomMm), lim = Math.max(0.1 * tNomMm, 0.5 * depth), row = map.t[iw].map((v) => (Number.isFinite(v) ? (loss ? v * toMM : tNomMm - v * toMM) : 0));
  let a = jw, b = jw; while (a > 0 && row[a - 1] >= lim) a--; while (b < row.length - 1 && row[b + 1] >= lim) b++;
  const dx = map.x.length > 1 ? Math.abs(map.x[Math.min(jw + 1, map.x.length - 1)] - map.x[Math.max(jw - 1, 0)]) / (Math.min(jw + 1, map.x.length - 1) - Math.max(jw - 1, 0)) : 0;
  return { minMm: tNomMm - depth, depthMm: depth, meanLossMm: clamp(s / n, 0, tNomMm), x: fin(map.x[jw], 0), theta: fin(map.theta?.[iw], 0), lengthMm: Math.max((Math.abs(map.x[b] - map.x[a]) + dx) * 1000, 10), readings: n };
}

// ---- erosion ------------------------------------------------------------------------------------------------------------
/** API RP 14E erosional velocity (m/s) for a mixture density (kg/m³) and the empirical C-factor in field units (100 continuous, 125 intermittent, up to 150–200 for corrosion-resistant service). */
export const erosionalVelocity = (rho, C = 100) => (1.21951 * C) / Math.sqrt(Math.max(rho, 1e-6));
/**
 * Sand erosion rate (mm/y) and its ingredients.
 * o: { model: 'dnv' | 'salama' | 'finnie' | 'oka', geometry: 'bend' | 'straight' | 'tee', mp (sand kg/s), U (particle ≈ mixture velocity m/s), D (bore m), dp (particle size m),
 *      rhoM, muM (mixture), rhoP (2650), rOverD (bend radius / D), gf (geometry factor), rhoT (target 7800), hv (Vickers hardness GPa), flowStress (Pa), mult }
 * 'dnv': DNV-RP-O501-type bend and straight-pipe equations (K = 2e-9, n = 2.6); 'salama': Salama (2000), S_m = 5.5 for elbows and 68 for plugged tees;
 * 'finnie' and 'oka': single-particle impact models applied over the DNV impact area and angle.
 */
export function erosionRate(o) {
  const D = +o.D, U = Math.max(+o.U, 0), mp = Math.max(nz(o.mp, 0), 0), dp = pos(o.dp, 250e-6), rhoM = pos(o.rhoM, 100), muM = pos(o.muM, 1e-4), rhoP = pos(o.rhoP, 2650), rhoT = pos(o.rhoT, 7800), mult = nz(o.mult, 1), geo = o.geometry || 'bend';
  need(D > 0, 'Erosion needs a positive pipe bore.');
  const alpha = Math.atan(1 / Math.sqrt(2 * pos(o.rOverD, 1.5))), sa = Math.sin(alpha), At = (Math.PI * D * D) / (4 * sa), A = (rhoM * rhoM * Math.tan(alpha) * U * D) / (rhoP * muM), gc0 = rhoM / (rhoP * (1.88 * Math.log(Math.max(A, 1e-12)) - 6.04)), gc = gc0 > 0 && gc0 < 0.1 ? gc0 : 0.1, G1 = dp / D < gc ? dp / D / gc : 1;
  const F = 0.6 * (sa + 7.2 * (sa - sa * sa)) ** 0.6 * (1 - Math.exp(-20 * alpha)), base = { alpha: (alpha * 180) / Math.PI, At, G: G1, F };
  if (!(mp > 0) || !(U > 0)) return { rate: 0, ...base };
  const perYear = 1000 * YEAR; // m/s of wall loss to mm/y
  const straight = () => ({ rate: mult * 2.5e-5 * U ** 2.6 * D ** -2 * mp, ...base, basis: 'DNV-RP-O501 straight pipe' });
  const salama = (Sm) => ({ rate: (mult * (mp * 86400) * U * U * (dp * 1e6)) / (Sm * (D * 1000) ** 2 * rhoM), ...base, Sm, basis: `Salama (2000), S_m = ${Sm}` });
  if (geo === 'straight') return straight(); // every model falls back to the straight-pipe equation away from fittings
  if (o.model === 'salama' || geo === 'tee') return salama(geo === 'tee' ? 68 : 5.5); // plugged tees are only covered by Salama's geometry constant
  if (o.model === 'finnie') { // cutting wear of a ductile target: volume = c·m·V²/(p·ψ·K)·f(α), ψ = K = 2
    const fa = Math.tan(alpha) <= 1 / 3 ? Math.sin(2 * alpha) - 3 * sa * sa : Math.cos(alpha) ** 2 / 3, vol = (nz(o.finnieC, 0.1) * U * U * fa) / (4 * pos(o.flowStress, 1.96e9)); // m³ per kg of sand
    return { rate: (mult * vol * mp * G1 * perYear) / At, ...base, perKg: vol * rhoT, basis: 'Finnie cutting model' };
  }
  if (o.model === 'oka') { // Oka et al. (2005), SiO2 particles: E90 = 65·Hv^-0.12·(v/104)^k2·(d/326 µm)^0.19 mm³/kg
    const hv = pos(o.hv, 1.96), n1 = 0.71 * hv ** 0.14, n2 = 2.4 * hv ** -0.94, g = sa ** n1 * (1 + hv * (1 - sa)) ** n2, e90 = 65 * hv ** -0.12 * (U / 104) ** (2.3 * hv ** 0.038) * (dp / 326e-6) ** 0.19;
    return { rate: (mult * g * e90 * 1e-9 * mp * G1 * perYear) / At, ...base, perKg: g * e90 * 1e-9 * rhoT, basis: 'Oka et al. (2005)' };
  }
  return { rate: ((mult * 2e-9 * F * U ** 2.6) / (rhoT * At)) * G1 * 2.5 * nz(o.gf, 1) * mp * 3.15e10, ...base, basis: 'DNV-RP-O501 bend' };
}

// ---- structural reliability ---------------------------------------------------------------------------------------------
/**
 * First-Order Reliability Method (Hasofer–Lind / Rackwitz–Fiessler iteration with a merit-function line search).
 * g: limit-state function of the physical variables (array) — failure when g < 0; vars: random-variable specifications (independent; each
 * marginal is mapped to a standard normal by F(x) = Φ(u), the Rosenblatt transform of independent variables).
 * Returns { beta, pf, alpha[] (unit vector towards the design point; α² = share of the uncertainty), u[], x[] (design point), names[], iterations, converged, evals, gradNorm }.
 */
export function form(g, vars, opt = {}) {
  const rv = vars.map(randomVariable), n = rv.length, toX = (u) => rv.map((r, i) => r.x(u[i])), h = opt.h || 1e-4; let evals = 0;
  const gu = (u) => { evals++; const v = g(toX(u)); return Number.isFinite(v) ? v : -1e300; };
  const grad = (u) => u.map((_, i) => { if (rv[i].dist === 'det') return 0; const a = u.slice(), b = u.slice(); a[i] += h; b[i] -= h; return (gu(a) - gu(b)) / (2 * h); });
  let u = (opt.start || new Array(n).fill(0)).slice(), gv = gu(u), gr = grad(u), it = 0, converged = false; const scale = Math.abs(gv) || 1;
  for (; it < (opt.maxIter || 100); it++) {
    const nr = norm2(gr); if (!(nr > 0)) break;
    const lam0 = (dot(gr, u) - gv) / (nr * nr), un = gr.map((gi) => lam0 * gi), d = un.map((x, i) => x - u[i]), nu = norm2(u);
    const c = 2 * Math.max(nu / nr, Math.abs(gv) > 1e-9 * scale ? (0.5 * dot(un, un)) / Math.abs(gv) : 0) + 1e-12, merit = (uu, gg) => 0.5 * dot(uu, uu) + c * Math.abs(gg), m0 = merit(u, gv);
    let lam = 1, ut = un, gt = gu(un);
    while (lam > 1 / 64 && merit(ut, gt) > m0) { lam /= 2; ut = u.map((x, i) => x + lam * d[i]); gt = gu(ut); }
    const du = norm2(ut.map((x, i) => x - u[i])); u = ut; gv = gt; gr = grad(u);
    if (du < 1e-7 * (1 + norm2(u)) && Math.abs(gv) < 1e-7 * scale) { converged = true; it++; break; }
  }
  const nr = norm2(gr) || 1, alpha = gr.map((gi) => -gi / nr), beta = dot(alpha, u);
  return { beta, pf: Phi(-beta), alpha, u, x: toX(u), names: rv.map((r) => r.name), iterations: it, converged, evals, gradNorm: nr };
}
/**
 * Second-Order Reliability Method: Breitung's curvature correction of a FORM result, pf ≈ Φ(−β)·Π(1 + β·κᵢ)^(−1/2).
 * The principal curvatures κᵢ come from the Hessian of g in standard-normal space, rotated so that the last axis points at the design point.
 * Returns { beta (generalised), pf, kappa[], valid, form }.
 */
export function sorm(g, vars, f = null, opt = {}) {
  const rv = vars.map(randomVariable), n = rv.length, res = f || form(g, vars), u = res.u, gu = (uu) => g(rv.map((r, i) => r.x(uu[i]))), h = opt.h || 0.05, g0 = gu(u);
  const act = rv.map((r, i) => i).filter((i) => rv[i].dist !== 'det'), m = act.length;
  if (m < 2 || !(res.beta > 0)) return { beta: res.beta, pf: res.pf, kappa: [], valid: m >= 1, form: res };
  const H = zeros(m), sh = (i, a, j, b) => { const x = u.slice(); x[act[i]] += a; if (j >= 0) x[act[j]] += b; return gu(x); };
  for (let i = 0; i < m; i++) { H[i][i] = (sh(i, h, -1, 0) - 2 * g0 + sh(i, -h, -1, 0)) / (h * h); for (let j = i + 1; j < m; j++) H[i][j] = H[j][i] = (sh(i, h, j, h) - sh(i, h, j, -h) - sh(i, -h, j, h) + sh(i, -h, j, -h)) / (4 * h * h); }
  // orthonormal basis with α as the last vector (Gram–Schmidt)
  const al = act.map((i) => res.alpha[i]), na = norm2(al) || 1, basis = [al.map((x) => x / na)];
  for (let k = 0; k < m && basis.length < m; k++) { const e = new Array(m).fill(0); e[k] = 1; for (const b of basis) { const p = dot(e, b); for (let i = 0; i < m; i++) e[i] -= p * b[i]; } const ne = norm2(e); if (ne > 1e-8) basis.push(e.map((x) => x / ne)); }
  const R = basis.slice(1), A = zeros(m - 1); for (let i = 0; i < m - 1; i++) { const Hr = matVec(H, R[i]); for (let j = 0; j < m - 1; j++) A[i][j] = dot(R[j], Hr) / res.gradNorm; }
  const kappa = jacobiEig(A).values; let fac = 1, valid = true;
  for (const k of kappa) { const t = 1 + res.beta * k; if (t <= 0.05) { valid = false; break; } fac /= Math.sqrt(t); }
  const pf = valid ? clamp(res.pf * fac, 0, 1) : res.pf;
  return { beta: -PhiInv(pf), pf, kappa, valid, form: res };
}
/**
 * Sampling estimate of the failure probability. opt: { n, seed, method: 'crude' | 'lhs' (Latin hypercube) | 'is' (importance sampling: unit-variance
 * normal density centred on `center`, normally the FORM design point), center, keep (return the g values) }
 * Returns { pf, cov (coefficient of variation of the estimate), beta (null when no failure was sampled), n, failures, gMean, gStd, g (Float64Array when keep) }.
 */
export function monteCarlo(g, vars, opt = {}) {
  const rv = vars.map(randomVariable), d = rv.length, n = Math.max(10, Math.round(opt.n || 10000)), method = opt.method || 'crude', r = rng(opt.seed ?? 20240), c = opt.center || new Array(d).fill(0), cc = dot(c, c);
  const U = method === 'lhs' ? lhs(n, d, opt.seed ?? 7) : null, keep = opt.keep ? new Float64Array(n) : null, u = new Array(d), x = new Array(d);
  let s1 = 0, s2 = 0, fails = 0, gs = 0, gss = 0;
  for (let k = 0; k < n; k++) {
    let w = 1;
    if (method === 'lhs') for (let i = 0; i < d; i++) u[i] = PhiInv(U[k][i]);
    else if (method === 'is') { let cz = 0; for (let i = 0; i < d; i++) { const z = rv[i].dist === 'det' ? 0 : r.normal(); u[i] = c[i] + z; cz += c[i] * z; } w = Math.exp(-cz - 0.5 * cc); }
    else for (let i = 0; i < d; i++) u[i] = r.normal();
    for (let i = 0; i < d; i++) x[i] = rv[i].x(u[i]);
    const gv = g(x); if (keep) keep[k] = gv; gs += gv; gss += gv * gv;
    if (!(gv >= 0)) { fails++; s1 += w; s2 += w * w; }
  }
  const pf = s1 / n, vr = Math.max(s2 / n - pf * pf, 0) / n;
  return { pf, cov: pf > 0 ? Math.sqrt(vr) / pf : null, beta: pf > 0 && pf < 1 ? -PhiInv(pf) : null, n, failures: fails, gMean: gs / n, gStd: Math.sqrt(Math.max(gss / n - (gs / n) ** 2, 0)), g: keep };
}
/**
 * Response-surface reliability (Bucher–Bourgund): a quadratic polynomial without cross terms is fitted to the limit state in standard-normal space
 * (2n + 1 evaluations), re-centred once on the line towards its design point, then solved by FORM and by importance sampling on the surface.
 * Returns { beta, pf, pfSampling, evals (calls of the true limit state), center[], coeff: { a, b[], c[] } }.
 */
export function responseSurface(g, vars, opt = {}) {
  const rv = vars.map(randomVariable), n = rv.length, gu = (u) => g(rv.map((r, i) => r.x(u[i]))), hh = opt.f || 2, std = rv.map((r) => ({ name: r.name, dist: r.dist === 'det' ? 'det' : 'normal', mean: 0, sd: r.dist === 'det' ? 0 : 1 })); let evals = 0;
  const fit = (c) => { const a = gu(c); evals++; const b = [], q = []; for (let i = 0; i < n; i++) { if (rv[i].dist === 'det') { b.push(0); q.push(0); continue; } const up = c.slice(), um = c.slice(); up[i] += hh; um[i] -= hh; const gp = gu(up), gm = gu(um); evals += 2; b.push((gp - gm) / (2 * hh)); q.push((gp + gm - 2 * a) / (2 * hh * hh)); } return { a, b, c: q, center: c.slice(), f: (u) => { let s = a; for (let i = 0; i < n; i++) { const dlt = u[i] - c[i]; s += b[i] * dlt + q[i] * dlt * dlt; } return s; } }; };
  let s = fit(new Array(n).fill(0)), r = form(s.f, std);
  for (let pass = 0; pass < (opt.passes || 2); pass++) {
    const gD = gu(r.u); evals++;
    const den = s.a - gD, c = Math.abs(den) > 1e-12 * Math.abs(s.a) ? r.u.map((x, i) => s.center[i] + ((x - s.center[i]) * s.a) / den) : r.u.slice();
    s = fit(c.map((x) => clamp(x, -8, 8))); r = form(s.f, std);
  }
  const mc = monteCarlo(s.f, std, { n: opt.n || 4000, method: 'is', center: r.u, seed: opt.seed ?? 99 });
  return { beta: r.beta, pf: r.pf, pfSampling: mc.pf, evals, center: s.center, coeff: { a: s.a, b: s.b, c: s.c }, alpha: r.alpha };
}

// ---- risk models --------------------------------------------------------------------------------------------------------
/**
 * Fault tree. Node: { name, gate: 'OR' | 'AND', children: [...] } or basic event { name, p }. Basic events with the same name are the same event.
 * Returns { top (exact by enumeration of the basic events, independence assumed), cutSets: [{ events[], p }] (minimal, by top-down expansion with absorption),
 *   rareEvent (Σ cut sets), upperBound (1 − Π(1 − P)), importance: [{ name, p, birnbaum, fussellVesely }] }.
 */
export function faultTree(tree) {
  const basics = new Map(); (function walk(nd) { if (nd.children) nd.children.forEach(walk); else basics.set(nd.name, clamp(+nd.p || 0, 0, 1)); })(tree);
  const names = [...basics.keys()], idx = new Map(names.map((nm, i) => [nm, i])), nb = names.length;
  need(nb >= 1 && nb <= 20, 'The fault tree needs between 1 and 20 basic events.');
  const ev = (nd, st) => (nd.children ? (nd.gate === 'AND' ? nd.children.every((c) => ev(c, st)) : nd.children.some((c) => ev(c, st))) : (st >> idx.get(nd.name)) & 1);
  const pv = names.map((nm) => basics.get(nm)); let top = 0; const on = new Float64Array(nb);
  for (let st = 0; st < 1 << nb; st++) { let pr = 1; for (let i = 0; i < nb && pr > 0; i++) pr *= (st >> i) & 1 ? pv[i] : 1 - pv[i]; if (pr > 0 && ev(tree, st)) { top += pr; for (let i = 0; i < nb; i++) if ((st >> i) & 1) on[i] += pr; } }
  const cs = (nd) => { if (!nd.children) return [new Set([nd.name])]; const parts = nd.children.map(cs); if (nd.gate !== 'AND') return parts.flat(); return parts.reduce((acc, p) => acc.flatMap((a) => p.map((b) => new Set([...a, ...b])))); };
  let sets = cs(tree).sort((a, b) => a.size - b.size); const minimal = [];
  for (const s of sets) if (!minimal.some((m) => [...m].every((e) => s.has(e)))) minimal.push(s);
  const cut = minimal.map((s) => ({ events: [...s], p: [...s].reduce((p, e) => p * basics.get(e), 1) })).sort((a, b) => b.p - a.p), rare = sum(cut.map((c) => c.p));
  return { top, cutSets: cut, rareEvent: rare, upperBound: 1 - cut.reduce((p, c) => p * (1 - c.p), 1), importance: names.map((nm, i) => ({ name: nm, p: pv[i], birnbaum: (pv[i] > 0 ? on[i] / pv[i] : 0) - (pv[i] < 1 ? (top - on[i]) / (1 - pv[i]) : 0), fussellVesely: rare > 0 ? sum(cut.filter((c) => c.events.includes(nm)).map((c) => c.p)) / rare : 0 })) };
}
/**
 * Event tree. freq: frequency of the initiating event; branches: [{ name, p (probability of "yes") }] in chronological order.
 * Returns { outcomes: [{ path: [bool…], label, p, freq }], total (Σ freq) }.
 */
export function eventTree(freq, branches) {
  const nb = branches.length, outcomes = [];
  for (let st = 0; st < 1 << nb; st++) { let p = 1; const path = []; for (let i = 0; i < nb; i++) { const yes = !((st >> (nb - 1 - i)) & 1); path.push(yes); p *= yes ? clamp(branches[i].p, 0, 1) : 1 - clamp(branches[i].p, 0, 1); } outcomes.push({ path, label: branches.map((b, i) => `${b.name}: ${path[i] ? 'yes' : 'no'}`).join(' · '), p, freq: freq * p }); }
  return { outcomes, total: sum(outcomes.map((o) => o.freq)) };
}
const matMul = (A, B) => { const n = A.length, C = zeros(n); for (let i = 0; i < n; i++) for (let k = 0; k < n; k++) { const a = A[i][k]; if (a !== 0) for (let j = 0; j < n; j++) C[i][j] += a * B[k][j]; } return C; };
/**
 * Continuous-time Markov chain. rates[i][j] = transition rate from state i to state j (per unit time, off-diagonal); p0 = initial probabilities.
 * The transition matrix over one step is the matrix exponential (scaling and squaring). Returns { t[], p[step][state], P (one-step matrix) }.
 */
export function markov(rates, p0, tEnd, steps = 100) {
  const n = rates.length, dt = tEnd / steps, Q = zeros(n);
  for (let i = 0; i < n; i++) { let s = 0; for (let j = 0; j < n; j++) if (j !== i) { const q = Math.max(+rates[i][j] || 0, 0); Q[i][j] = q * dt; s += q * dt; } Q[i][i] = -s; }
  let nrm = 0; for (let i = 0; i < n; i++) nrm = Math.max(nrm, 2 * Math.abs(Q[i][i]));
  const sq = Math.max(0, Math.ceil(Math.log2(Math.max(nrm, 1e-12) / 0.25))), A = Q.map((r) => r.map((v) => v / 2 ** sq));
  let P = zeros(n), term = zeros(n); for (let i = 0; i < n; i++) P[i][i] = term[i][i] = 1;
  for (let k = 1; k <= 16; k++) { term = matMul(term, A).map((r) => r.map((v) => v / k)); for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) P[i][j] += term[i][j]; }
  for (let s = 0; s < sq; s++) P = matMul(P, P);
  const t = [0], p = [Array.from(p0)];
  for (let k = 1; k <= steps; k++) { const prev = p[k - 1], cur = new Array(n).fill(0); for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) cur[j] += prev[i] * P[i][j]; p.push(cur); t.push(k * dt); }
  return { t, p, P: P.map((r) => Array.from(r)) };
}
/**
 * Discrete Bayesian network by exact enumeration. nodes: [{ name, states: [...], parents: [names], cpt: [[P(state | parent combination)…]…] }] with the
 * parent combinations in row-major order (first parent varies slowest); evidence: { name: state }.
 * Returns { posterior: { name: { state: probability } }, pEvidence }.
 */
export function bayesNet(nodes, evidence = {}) {
  const n = nodes.length, pos = new Map(nodes.map((nd, i) => [nd.name, i])), card = nodes.map((nd) => nd.states.length), total = card.reduce((a, b) => a * b, 1);
  need(total <= 2e5, 'The Bayesian network is too large for exact enumeration.');
  const par = nodes.map((nd) => (nd.parents || []).map((p) => { need(pos.has(p), `Unknown parent "${p}" in the Bayesian network.`); return pos.get(p); })), evi = nodes.map((nd) => (nd.name in evidence ? nd.states.indexOf(evidence[nd.name]) : -1));
  const acc = nodes.map((nd) => new Float64Array(nd.states.length)), st = new Array(n).fill(0); let pe = 0;
  for (let k = 0; k < total; k++) {
    let r = k, okE = true; for (let i = n - 1; i >= 0; i--) { st[i] = r % card[i]; r = Math.floor(r / card[i]); if (evi[i] >= 0 && st[i] !== evi[i]) okE = false; }
    if (!okE) continue;
    let p = 1; for (let i = 0; i < n && p > 0; i++) { let row = 0; for (const j of par[i]) row = row * card[j] + st[j]; p *= nodes[i].cpt[row][st[i]]; }
    if (p > 0) { pe += p; for (let i = 0; i < n; i++) acc[i][st[i]] += p; }
  }
  return { pEvidence: pe, posterior: Object.fromEntries(nodes.map((nd, i) => [nd.name, Object.fromEntries(nd.states.map((s, j) => [s, pe > 0 ? acc[i][j] / pe : 0]))])) };
}
/**
 * Remaining useful life with uncertainty. o: { margin (mm of wall that may still be lost), rate (mm/y model estimate), cov (model uncertainty),
 *   measured: { loss (mm), years, sd (mm measurement error) } (optional: updates the rate, precision-weighted) }
 * Returns { rate, rateSd, p10, p50, p90 (years; p10 = pessimistic), mean, updated }.
 */
export function remainingLife(o) {
  let mu = Math.max(+o.rate, 1e-9), sd = Math.max(mu * nz(o.cov, 0.5), 1e-12), updated = false;
  const me = o.measured;
  if (me && me.years > 0 && Number.isFinite(me.loss)) { const r = Math.max(me.loss, 0) / me.years, so = Math.max(nz(me.sd, 0.3) / me.years, 1e-9), w0 = 1 / (sd * sd), w1 = 1 / (so * so); mu = Math.max((mu * w0 + r * w1) / (w0 + w1), 1e-9); sd = Math.sqrt(1 / (w0 + w1)); updated = true; }
  const rv = randomVariable({ dist: 'lognormal', mean: mu, sd }), m = Math.max(+o.margin, 0), life = (q) => cap(m / rv.x(PhiInv(q)), 500);
  return { rate: mu, rateSd: sd, p10: life(0.9), p50: life(0.5), p90: life(0.1), mean: cap(m / mu, 500), updated };
}
const POF_EDGES = [3e-5, 3e-4, 3e-3, 3e-2], COF_EDGES = [1e5, 1e6, 1e7, 1e8], RISK_NAMES = ['low', 'medium', 'high', 'very high'];
const category = (x, edges) => 1 + edges.filter((e) => x >= e).length;
/** Risk level (0 low … 3 very high) of a probability category and a consequence category (both 1–5) on the 5 × 5 matrix. */
export const riskLevel = (pc, cc) => { const s = pc + cc; return s <= 4 ? 0 : s <= 6 ? 1 : s <= 8 ? 2 : 3; };

// ---- suite declaration --------------------------------------------------------------------------------------------------
const opt = (pairs) => pairs.map(([value, label]) => ({ value, label }));
const SN_OPTIONS = opt(Object.keys(SN_CURVES).map((k) => [k, `${k}${k === 'D' ? ' (girth weld cap, good profile)' : k === 'E' ? ' (girth weld cap)' : k === 'F1' ? ' (single-sided girth weld root)' : k === 'F' ? ' (root with backing)' : k === 'C1' ? ' (ground flush weld / seamless pipe body)' : ''}`]));
const INPUTS = [
  { group: 'Pipe and material', tab: 'inputs', help: 'Steel pipe of the flowline and riser. Linked values come from the network suite.', fields: [
    { key: 'idMm', label: 'Inner diameter', unit: 'mm', value: BASE.idMm, min: 20, max: 1500, typical: [100, 900] },
    { key: 'wtMm', label: 'Nominal wall thickness', unit: 'mm', value: BASE.wtMm, min: 2, max: 80, typical: [6, 40] },
    { key: 'odCoatMm', label: 'Outer diameter including coatings', unit: 'mm', value: BASE.idMm + 2 * BASE.wtMm + 2000 * BASE.insulation.t, min: 0, max: 2500, help: 'Used for buoyancy, added mass and current loads. 0 = bare steel.' },
    { key: 'coatDensity', label: 'Coating density', unit: 'kg/m³', value: 900, min: 100, max: 3500 },
    { key: 'material', label: 'Material grade', type: 'text', value: 'API 5L X65' },
    { key: 'smys', label: 'Specified minimum yield strength', unit: 'MPa', value: BASE.smys, min: 150, max: 900 },
    { key: 'smts', label: 'Specified minimum tensile strength', unit: 'MPa', value: BASE.smts, min: 250, max: 1100 },
    { key: 'eMod', label: 'Young\'s modulus', unit: 'MPa', value: BASE.E, min: 5e4, max: 2.5e5 },
    { key: 'poisson', label: 'Poisson\'s ratio', unit: '–', value: BASE.poisson, min: 0.2, max: 0.4 },
    { key: 'alphaT', label: 'Thermal expansion coefficient', unit: '1/K', value: BASE.alphaT, min: 5e-6, max: 2.5e-5 },
    { key: 'rhoSteel', label: 'Steel density', unit: 'kg/m³', value: BASE.rhoSteel, min: 6000, max: 9000 },
    { key: 'corrAllow', label: 'Corrosion allowance', unit: 'mm', value: BASE.corrosionAllowanceMm, min: 0, max: 15 },
    { key: 'eroAllow', label: 'Erosion allowance', unit: 'mm', value: 0, min: 0, max: 10 },
    { key: 'ovality', label: 'Out-of-roundness (Dmax − Dmin)/D', unit: '%', value: 0.75, min: 0.5, max: 3 },
    { key: 'tFabMm', label: 'Fabrication tolerance on the wall', unit: 'mm', value: 1, min: 0, max: 5, help: 'Deducted in the DNV pressure-containment check.' },
  ] },
  { group: 'Design basis and code', tab: 'inputs', fields: [
    { key: 'code', label: 'Design code format', type: 'select', value: 'dnv', options: opt([['dnv', 'DNV-ST-F101 (limit state, safety class)'], ['b314', 'ASME B31.4 (liquid lines)'], ['b318', 'ASME B31.8 (gas lines)']]) },
    { key: 'safetyClass', label: 'Safety class', type: 'select', value: 'medium', options: opt([['low', 'Low'], ['medium', 'Medium'], ['high', 'High'], ['veryHigh', 'Very high']]), help: 'Sets the DNV resistance factors and the target annual failure probability (1e-3, 1e-4, 1e-5, 1e-6).' },
    { key: 'designFactor', label: 'ASME design factor, pipeline', unit: '–', value: 0.72, min: 0.3, max: 0.8, showIf: (v) => v.code !== 'dnv' },
    { key: 'riserFactor', label: 'ASME design factor, riser', unit: '–', value: 0.6, min: 0.3, max: 0.8, showIf: (v) => v.code !== 'dnv', help: '0.60 in B31.4, 0.50 in B31.8 for risers and platform piping.' },
    { key: 'designPressure', label: 'Design pressure', unit: 'bara', value: BASE.designPressure, min: 2, max: 1500 },
    { key: 'pRefLoc', label: 'Design pressure refers to', type: 'select', value: 'inlet', options: opt([['inlet', 'The inlet (subsea wellhead or manifold)'], ['top', 'The highest point (topsides)'], ['uniform', 'Every point (no static head)']]) },
    { key: 'rhoContent', label: 'Content density for the static head', unit: 'kg/m³', value: 700, min: 1, max: 1300 },
    { key: 'designTemp', label: 'Design temperature', unit: '°C', value: BASE.designTemp, min: -50, max: 250 },
    { key: 'mdmt', label: 'Minimum design metal temperature', unit: '°C', value: -20, min: -120, max: 30 },
    { key: 'blowdownMinT', label: 'Lowest temperature during blowdown', unit: '°C', value: 0, min: -120, max: 60, help: 'From the operations suite when linked.' },
    { key: 'tInstall', label: 'Installation (tie-in) temperature', unit: '°C', value: BASE.tSeabed, min: -10, max: 50 },
    { key: 'layTension', label: 'Residual lay tension', unit: 'kN', value: 0, min: 0, max: 5000 },
    { key: 'restraint', label: 'Axial restraint', type: 'select', value: 'auto', options: opt([['auto', 'Flowline restrained, riser free'], ['restrained', 'Fully restrained everywhere'], ['unrestrained', 'Free to expand (end-cap load)']]) },
    { key: 'hydrotestFactor', label: 'Hydrotest pressure / design pressure', unit: '–', value: 1.25, min: 1, max: 1.5 },
    { key: 'designLife', label: 'Design life', unit: 'y', value: BASE.projectLife, min: 1, max: 60 },
    { key: 'ageY', label: 'Age of the line today', unit: 'y', value: 0, min: 0, max: 60 },
  ] },
  { group: 'Site and seabed', tab: 'inputs', fields: [
    { key: 'waterDepth', label: 'Maximum water depth', unit: 'm', value: BASE.waterDepth, min: 0, max: 4000, help: 'The depths of the route profile are scaled to this value; 0 = onshore.' },
    { key: 'currentSpeed', label: 'Near-bed current speed (long-term mean)', unit: 'm/s', value: BASE.currentSpeed, min: 0, max: 3 },
    { key: 'waveHeight', label: 'Significant wave height', unit: 'm', value: 2, min: 0, max: 20 },
    { key: 'wavePeriod', label: 'Peak wave period', unit: 's', value: 9, min: 2, max: 25 },
    { key: 'muAxial', label: 'Axial pipe–soil friction', unit: '–', value: 0.5, min: 0.05, max: 1.5 },
    { key: 'muLateral', label: 'Lateral pipe–soil friction', unit: '–', value: 0.6, min: 0.05, max: 2 },
    { key: 'coverDepth', label: 'Burial cover above the pipe', unit: 'm', value: 0, min: 0, max: 5 },
    { key: 'soilGamma', label: 'Submerged unit weight of the cover', unit: 'kN/m³', value: 9, min: 3, max: 14 },
    { key: 'imperfection', label: 'Seabed imperfection height', unit: 'm', value: 0.3, min: 0, max: 2 },
  ] },
  { group: 'Slug loads and bends', tab: 'inputs', help: 'Slug data come from the flow suite when it has been run, otherwise from the kernel slug closures. Zero means automatic.', fields: [
    { key: 'slugMode', label: 'Slug loading', type: 'select', value: 'auto', options: opt([['auto', 'Automatic (from the flow regime)'], ['on', 'Always apply'], ['off', 'Ignore']]) },
    { key: 'slugFreq', label: 'Slug frequency', unit: '1/s', value: 0, min: 0, max: 5 },
    { key: 'slugLen', label: 'Slug body length', unit: 'm', value: 0, min: 0, max: 2000 },
    { key: 'slugVel', label: 'Slug velocity', unit: 'm/s', value: 0, min: 0, max: 60 },
    { key: 'slugHoldup', label: 'Liquid holdup in the slug body', unit: '–', value: 0, min: 0, max: 1 },
    { key: 'slugDensity', label: 'Slug body density', unit: 'kg/m³', value: 0, min: 0, max: 1300 },
    { key: 'slugMomentum', label: 'Slug momentum flux ρ·A·v²', unit: 'kN', value: 0, min: 0, max: 5000, help: 'Overrides the density when a measured or simulated momentum flux is available.' },
    { key: 'supportCapacity', label: 'Design capacity of the bend supports / anchors', unit: 'kN', value: 25, min: 0, max: 1e5, help: '0 = not checked.' },
    { key: 'dlf', label: 'Dynamic load factor', unit: '–', value: 2, min: 1, max: 4, help: '2 is the step-load upper bound for an undamped support.' },
    { key: 'bends', label: 'Bends', type: 'table', columns: [{ key: 'x', label: 'Distance', unit: 'm' }, { key: 'angle', label: 'Angle', unit: 'deg' }, { key: 'radius', label: 'Radius', unit: 'm' }], value: [{ x: 0, angle: 90, radius: 1.27 }, { x: 19560, angle: 90, radius: 1.27 }] },
    { key: 'pipingSupport', label: 'Support arrangement (screening)', type: 'select', value: 'mediumStiff', options: opt([['stiff', 'Stiff'], ['mediumStiff', 'Medium stiff'], ['medium', 'Medium']]) },
  ] },
  { group: 'Free span or jumper', tab: 'inputs', fields: [
    { key: 'spanLength', label: 'Span length', unit: 'm', value: 15, min: 1, max: 400 },
    { key: 'spanX', label: 'Span location along the line', unit: 'm', value: 9000, min: 0, max: 1e6 },
    { key: 'spanEnds', label: 'End conditions', type: 'select', value: 'pinned-pinned', options: opt([['pinned-pinned', 'Pinned – pinned'], ['fixed-fixed', 'Fixed – fixed'], ['fixed-pinned', 'Fixed – pinned'], ['fixed-free', 'Fixed – free (cantilever)'], ['springs', 'Elastic supports (springs)']]) },
    { key: 'spanKt', label: 'Support stiffness, translation', unit: 'kN/m', value: 5000, min: 0, max: 1e7, showIf: (v) => v.spanEnds === 'springs' },
    { key: 'spanKr', label: 'Support stiffness, rotation', unit: 'kN·m/rad', value: 20000, min: 0, max: 1e8, showIf: (v) => v.spanEnds === 'springs' },
    { key: 'midSupports', label: 'Intermediate supports', type: 'table', columns: [{ key: 'x', label: 'Position from the left end', unit: 'm' }], value: [] },
    { key: 'spanAxial', label: 'Include the effective axial force', type: 'bool', value: false, help: 'Compression lowers the natural frequency (geometric stiffness).' },
    { key: 'timoshenko', label: 'Shear deformation and rotary inertia', type: 'bool', value: false },
    { key: 'addedMass', label: 'Added-mass coefficient', unit: '–', value: 1, min: 0, max: 3 },
    { key: 'damping', label: 'Damping ratio (structure + soil + fluid)', unit: '%', value: 2, min: 0.1, max: 20 },
    { key: 'bendAt', label: 'Bend position on the span (jumper)', unit: 'fraction', value: 0, min: 0, max: 0.95, help: '0 = straight span: only the moving slug weight loads it.' },
    { key: 'bendAngle', label: 'Angle of that bend', unit: 'deg', value: 90, min: 0, max: 180 },
  ] },
  { group: 'Pressure and temperature history', tab: 'inputs', help: 'Optional operating log. Its pressure and temperature cycles are rainflow-counted into the fatigue spectrum and checked against the limits.', fields: [
    { key: 'opLog', label: 'Operating log', type: 'table', columns: [{ key: 't', label: 'Time', unit: 'h' }, { key: 'p', label: 'Inlet pressure', unit: 'bara' }, { key: 'T', label: 'Inlet temperature', unit: '°C' }], value: [] },
  ] },
  { group: 'Fatigue and fracture', tab: 'inputs', fields: [
    { key: 'snClass', label: 'S–N class', type: 'select', value: 'F1', options: SN_OPTIONS },
    { key: 'snEnv', label: 'Environment', type: 'select', value: 'cp', options: opt([['air', 'Air'], ['cp', 'Seawater with cathodic protection'], ['free', 'Seawater, free corrosion']]) },
    { key: 'scf', label: 'Stress concentration factor', unit: '–', value: 1.3, min: 1, max: 6 },
    { key: 'dff', label: 'Design fatigue factor', unit: '–', value: 6, min: 1, max: 10 },
    { key: 'priorDamage', label: 'Fatigue damage already accumulated', unit: '–', value: 0, min: 0, max: 1 },
    { key: 'eventsShutdown', label: 'Shutdowns', unit: '1/y', value: 6, min: 0, max: 365 },
    { key: 'eventsBlowdown', label: 'Blowdowns', unit: '1/y', value: 1, min: 0, max: 100 },
    { key: 'flawMm', label: 'Initial flaw depth', unit: 'mm', value: 1, min: 0.05, max: 20 },
    { key: 'parisC', label: 'Paris coefficient C', unit: 'm/cycle, MPa√m', value: 1.65e-11, min: 1e-14, max: 1e-8, help: '1.65e-11 with m = 3 is the simplified law for steels in air; 7.3e-11 for a marine environment.' },
    { key: 'parisM', label: 'Paris exponent m', unit: '–', value: 3, min: 2, max: 5 },
    { key: 'dKth', label: 'Threshold ΔK', unit: 'MPa√m', value: 2, min: 0, max: 10 },
    { key: 'kMat', label: 'Fracture toughness', unit: 'MPa√m', value: 150, min: 20, max: 400 },
  ] },
  { group: 'Corrosion', tab: 'inputs', fields: [
    { key: 'co2', label: 'CO2 in the gas phase', unit: 'mol %', value: 2, min: 0, max: 80 },
    { key: 'h2s', label: 'H2S in the gas phase', unit: 'mol %', value: 0, min: 0, max: 40 },
    { key: 'corrModel', label: 'CO2 corrosion model', type: 'select', value: '1995', options: opt([['1995', 'de Waard–Lotz–Milliams 1995 (with flow velocity)'], ['1991', 'de Waard–Milliams 1991 (nomogram equation)']]) },
    { key: 'phAct', label: 'In-situ pH', unit: '–', value: 0, min: 0, max: 8, help: '0 = water saturated with CO2 (no buffering).' },
    { key: 'inhibEff', label: 'Inhibitor efficiency × availability', unit: '%', value: 98.5, min: 0, max: 99.9 },
    { key: 'corrMult', label: 'Corrosion-model multiplier', unit: '–', value: 1, min: 0.05, max: 5 },
    { key: 'glycolWt', label: 'Glycol in the water phase', unit: 'wt %', value: 0, min: 0, max: 95 },
    { key: 'wetModel', label: 'Water wetting', type: 'select', value: 'water', options: opt([['water', 'Always water-wet (conservative)'], ['dewaard', 'Oil-wet when water cut < 30 % and liquid velocity > 1 m/s']]) },
    { key: 'tauCrit', label: 'Wall shear that strips the inhibitor film', unit: 'Pa', value: 150, min: 5, max: 2000 },
    { key: 'minWt', label: 'Measured minimum wall thickness', unit: 'mm', value: 0, min: 0, max: 80, help: '0 = no measurement: the wall is predicted from the age and the rates.' },
    { key: 'wtMap', label: 'Wall-thickness or metal-loss map', type: 'file', value: null },
    { key: 'defects', label: 'Measured metal-loss defects', type: 'table', columns: [{ key: 'x', label: 'Location', unit: 'm' }, { key: 'depth', label: 'Depth', unit: 'mm' }, { key: 'length', label: 'Axial length', unit: 'mm' }], value: [] },
  ] },
  { group: 'Sand and erosion', tab: 'inputs', fields: [
    { key: 'sandKgD', label: 'Sand production', unit: 'kg/d', value: 10, min: 0, max: 50000 },
    { key: 'sandUm', label: 'Particle size', unit: 'µm', value: 250, min: 10, max: 3000 },
    { key: 'sandDensity', label: 'Particle density', unit: 'kg/m³', value: 2650, min: 1000, max: 8000 },
    { key: 'erosionModel', label: 'Erosion model that governs', type: 'select', value: 'governing', options: opt([['governing', 'Larger of DNV-RP-O501 and Salama'], ['dnv', 'DNV-RP-O501'], ['salama', 'Salama (2000)'], ['oka', 'Oka (2005)'], ['finnie', 'Finnie']]) },
    { key: 'erosionMult', label: 'Erosion-model multiplier', unit: '–', value: 1, min: 0.05, max: 20 },
    { key: 'geomFactor', label: 'Geometry factor of the fittings', unit: '–', value: 1, min: 1, max: 4 },
    { key: 'c14e', label: 'API RP 14E C-factor', unit: 'field units', value: 100, min: 50, max: 300 },
  ] },
  { group: 'Blockage and plugs', tab: 'inputs', fields: [
    { key: 'blockage', label: 'Flow area lost to deposits', unit: 'fraction', value: 0, min: 0, max: 1 },
    { key: 'effIdMm', label: 'Smallest effective bore', unit: 'mm', value: 0, min: 0, max: 1500, help: '0 = clean pipe.' },
    { key: 'depositRate', label: 'Deposit growth rate', unit: 'mm/d', value: 0, min: 0, max: 50, help: 'Thickness growth on the wall, for example the wax rate of the solids suite.' },
    { key: 'blockLimit', label: 'Allowed loss of flow area', unit: 'fraction', value: 0.5, min: 0.05, max: 1 },
    { key: 'roughUm', label: 'Effective wall roughness (with deposits)', unit: 'µm', value: BASE.roughUm, min: 0.1, max: 5000 },
    { key: 'plugX', label: 'Plug or worst deposit location', unit: 'm', value: 0, min: 0, max: 1e6 },
    { key: 'plugProb', label: 'Annual probability of a hydrate plug', unit: '–', value: 0.02, min: 0, max: 1 },
    { key: 'plugLen', label: 'Plug length', unit: 'm', value: 30, min: 1, max: 2000 },
    { key: 'plugPorosity', label: 'Plug porosity', unit: '–', value: 0.4, min: 0, max: 0.9 },
    { key: 'pShutIn', label: 'Source shut-in pressure', unit: 'bara', value: 280, min: 2, max: 1500 },
    { key: 'plugDays', label: 'Time to remediate a plug', unit: 'd', value: 30, min: 0, max: 365 },
  ] },
  { group: 'Reliability', tab: 'setup', help: 'Scatter of the random variables in the limit states.', fields: [
    { key: 'covWt', label: 'Wall thickness scatter', unit: '%', value: 3, min: 0, max: 20 },
    { key: 'covYield', label: 'Yield strength scatter', unit: '%', value: 5, min: 0, max: 20 },
    { key: 'yieldBias', label: 'Mean yield / SMYS', unit: '–', value: 1.08, min: 1, max: 1.4 },
    { key: 'covPress', label: 'Annual maximum pressure scatter', unit: '%', value: 5, min: 0, max: 30 },
    { key: 'covCorr', label: 'Wall-loss model uncertainty', unit: '%', value: 50, min: 0, max: 200 },
    { key: 'leakFrac', label: 'Failure criterion: wall loss at a leak / wall thickness', unit: '–', value: 0.8, min: 0.3, max: 1 },
    { key: 'seed', label: 'Random seed', unit: '–', value: 11, min: 1, max: 1e6 },
  ] },
  { group: 'Risk and consequence', tab: 'setup', fields: [
    { key: 'oilPrice', label: 'Oil price', unit: '$/bbl', value: BASE.oilPrice, min: 1, max: 300 },
    { key: 'repairCost', label: 'Repair or intervention cost', unit: 'M$', value: 15, min: 0, max: 2000 },
    { key: 'downtimeDays', label: 'Downtime after a loss of containment', unit: 'd', value: 45, min: 0, max: 730 },
    { key: 'envCost', label: 'Environmental cost of an uncontrolled release', unit: 'M$', value: 40, min: 0, max: 20000 },
    { key: 'safetyCost', label: 'Safety cost of an ignited release', unit: 'M$', value: 100, min: 0, max: 20000 },
    { key: 'pDetect', label: 'Probability that a leak is detected', unit: '–', value: 0.9, min: 0, max: 1 },
    { key: 'pIsolate', label: 'Probability that isolation succeeds', unit: '–', value: 0.95, min: 0, max: 1 },
    { key: 'pIgnite', label: 'Probability of ignition', unit: '–', value: 0.02, min: 0, max: 1 },
    { key: 'pfdProtect', label: 'Overpressure protection, failure on demand', unit: '–', value: 0.01, min: 0, max: 1 },
    { key: 'extFreq', label: 'External damage frequency', unit: '1/y', value: 1e-4, min: 0, max: 1 },
    { key: 'inspInterval', label: 'Planned inspection interval', unit: 'y', value: 5, min: 0.25, max: 30 },
    { key: 'pod', label: 'Probability of detection per inspection', unit: '–', value: 0.8, min: 0, max: 1 },
  ] },
  { group: 'Discretisation', tab: 'mesh', fields: [
    { key: 'nElem', label: 'Beam elements on the span', unit: '–', value: 24, min: 4, max: 200 },
    { key: 'nModes', label: 'Modes kept', unit: '–', value: 6, min: 1, max: 12 },
    { key: 'stepsPerCycle', label: 'Time steps per first-mode period', unit: '–', value: 40, min: 8, max: 400 },
    { key: 'nSlugs', label: 'Slugs simulated', unit: '–', value: 12, min: 3, max: 60 },
    { key: 'nMC', label: 'Monte Carlo samples', unit: '–', value: 20000, min: 500, max: 400000 },
  ] },
];
const FIELDS = INPUTS.flatMap((g) => g.fields), DEFAULTS = Object.fromEntries(FIELDS.map((f) => [f.key, f.value]));

/** Validated engineering parameters (SI) from the raw inputs. */
function params(v0) {
  const v = { ...DEFAULTS, ...v0 }, numKey = (k) => { const x = +v[k]; need(Number.isFinite(x), `"${FIELDS.find((f) => f.key === k)?.label || k}" must be a number.`); return x; };
  const p = {}; for (const f of FIELDS) if (!f.type || f.type === 'number') p[f.key] = numKey(f.key); else p[f.key] = v[f.key];
  need(p.idMm > 0, 'The inner diameter must be positive.'); need(p.wtMm > 0, 'The wall thickness must be positive.'); need(p.wtMm < p.idMm, 'The wall thickness cannot exceed the inner diameter: check the units (both in mm).');
  need(p.smys > 0 && p.smts >= p.smys, 'The tensile strength must be at least the yield strength, and both positive.'); need(p.eMod > 1000, 'Young\'s modulus must be given in MPa (about 207000 for steel).');
  need(p.poisson > 0 && p.poisson < 0.5, 'Poisson\'s ratio must lie between 0 and 0.5.'); need(p.designPressure > 0, 'The design pressure must be positive.'); need(p.spanLength > 0, 'The span length must be positive.');
  need(p.designLife > 0, 'The design life must be positive.'); need(p.rhoSteel > 0, 'The steel density must be positive.');
  const ID = p.idMm / 1000, t = p.wtMm / 1000, D = ID + 2 * t, Dh = Math.max(nz(p.odCoatMm, 0) / 1000, D);
  Object.assign(p, { ID, t, D, Dh, Ai: (Math.PI * ID * ID) / 4, As: (Math.PI * (D * D - ID * ID)) / 4, Isteel: (Math.PI * (D ** 4 - ID ** 4)) / 64, E: p.eMod * MPA, S: p.smys * MPA, Su: p.smts * MPA, nu: p.poisson, f0: Math.max(p.ovality, 0.5) / 100, CA: Math.max(p.corrAllow, 0) + Math.max(p.eroAllow, 0),
    zeta: clamp(p.damping, 0.05, 50) / 100, eta: clamp(p.inhibEff, 0, 100) / 100, age: Math.max(p.ageY, 0), nElem: clamp(Math.round(p.nElem), 2, 200), nModes: clamp(Math.round(p.nModes), 1, 12), spc: clamp(Math.round(p.stepsPerCycle), 6, 2000), nSlugs: clamp(Math.round(p.nSlugs), 2, 200), nMC: clamp(Math.round(p.nMC), 200, 1e6),
    sc: DNV_SC[p.safetyClass] ? p.safetyClass : 'medium', code: ['dnv', 'b314', 'b318'].includes(p.code) ? p.code : 'dnv', snClass: SN_CURVES[p.snClass] ? p.snClass : 'F1', snEnv: ['air', 'cp', 'free'].includes(p.snEnv) ? p.snEnv : 'cp' });
  p.bends = (Array.isArray(p.bends) ? p.bends : []).map((b) => ({ x: nz(b?.x, 0), angle: clamp(nz(b?.angle, 90), 0, 180), radius: Math.max(nz(b?.radius, 0), 0) })).filter((b) => b.angle > 0).slice(0, 40);
  p.midSupports = (Array.isArray(p.midSupports) ? p.midSupports : []).map((s) => nz(s?.x, NaN)).filter((x) => x > 0 && x < p.spanLength).slice(0, 12);
  p.opLog = (Array.isArray(p.opLog) ? p.opLog : []).map((r) => ({ t: nz(r?.t, NaN), p: nz(r?.p, NaN), T: nz(r?.T, NaN) })).filter((r) => isNum(r.t) && isNum(r.p)).sort((a, b) => a.t - b.t).slice(0, 5000);
  p.defects = (Array.isArray(p.defects) ? p.defects : []).map((d) => ({ x: nz(d?.x, 0), depth: nz(d?.depth, 0), length: nz(d?.length, 0) })).filter((d) => d.depth > 0 && d.length > 0).slice(0, 60);
  return p;
}

/** Local flow, fluid and slug data at every node of the flow picture. */
function lineState(p, ctx) {
  const pic = flowPicture(ctx), fm = fluidModel(ctx), line = caseLine(ctx), n = pic.x.length, zMin = Math.min(...pic.z), zTop = Math.max(...pic.z), dMax = Math.max(0, -zMin), scale = p.waterDepth > 0 && dMax > 0 ? p.waterDepth / dMax : 0;
  const nodes = [];
  for (let i = 0; i < n; i++) {
    const a = Math.max(i - 1, 0), b = Math.min(i + 1, n - 1), theta = Math.atan2(pic.z[b] - pic.z[a], Math.max(pic.x[b] - pic.x[a], 1e-9)), P = Math.max(pic.P[i], 1), T = pic.T[i], o = fm.at(P, T);
    const vsl = Math.max(pic.vsl[i], 0), vsg = Math.max(pic.vsg[i], 0), vm = Math.max(pic.vm[i], vsl + vsg, 1e-6), hold = clamp(pic.holdup[i], 0.005, 1), depth = p.waterDepth > 0 ? (dMax > 0 ? Math.max(0, -pic.z[i]) * scale : p.waterDepth) : 0;
    const rhoNS = Math.max((o.rhoL * vsl + o.rhoG * vsg) / vm, 0.5), muNS = Math.max((o.muL * vsl + o.muG * vsg) / vm, 1e-6), su = slugUnit({ vsl: Math.max(vsl, 1e-4), vsg: Math.max(vsg, 1e-4), rhoL: o.rhoL, rhoG: o.rhoG, muL: o.muL, muG: o.muG, D: p.ID, theta });
    const fr = fin(su.freq, 0), len = fin(Math.min(su.length, su.lengthFromFreq), su.length), rhoS = o.rhoL * su.holdupSlug + o.rhoG * (1 - su.holdupSlug), rhoF = o.rhoL * su.holdupFilm + o.rhoG * (1 - su.holdupFilm), regime = String(pic.regime[i] || '');
    nodes.push({ i, x: pic.x[i], z: pic.z[i], theta, P, T, depth, pe: PATM + RHO_SW * G * depth, pi: P * BAR, vsl, vsg, vm, hold, vL: vsl / hold, rhoM: Math.max(pic.rhoM[i], 0.5), rhoNS, muNS, rhoL: o.rhoL, rhoG: o.rhoG, muL: o.muL, wcut: o.wcut, tauW: Math.abs(pic.tauW[i]), regime, slugLike: /slug|churn|intermittent/i.test(regime),
      slug: { freq: fr, len, v: su.vt, hls: su.holdupSlug, hlf: su.holdupFilm, rhoS, rhoF }, riser: pic.x[i] > line.riserBaseX + 1e-6 });
  }
  const at = (x) => { let k = 0; for (let i = 1; i < n; i++) if (Math.abs(nodes[i].x - x) < Math.abs(nodes[k].x - x)) k = i; return nodes[k]; };
  return { pic, fm, line, nodes, at, zTop, zRef: p.pRefLoc === 'top' ? zTop : pic.z[0], source: pic.source };
}

/** Design slug: user inputs first, then the flow suite's slug summary, then the kernel closure at the most loaded slugging node. */
function designSlug(p, st, ctx) {
  const fs = ctx.outputs?.flow?.slug, hasFs = fs && typeof fs === 'object' && fs.type && fs.type !== 'none' && fs.freq > 0 && fs.velocity > 0, cand = st.nodes.filter((nd) => nd.slugLike), pool = cand.length ? cand : st.nodes;
  let nd = pool.reduce((a, b) => (b.slug.rhoS * b.slug.v ** 2 > a.slug.rhoS * a.slug.v ** 2 ? b : a)), source = 'kernel slug closures';
  let { freq, len, v, hls } = nd.slug;
  if (hasFs) { if (isNum(fs.x)) nd = st.at(fs.x); freq = fs.freq; v = fs.velocity; len = pos(fs.length, len); hls = pos(fs.holdupBody, nd.slug.hls); source = `flow suite (${fs.type} slugging)`; }
  if (p.slugFreq > 0) { freq = p.slugFreq; source = 'inputs'; } if (p.slugVel > 0) { v = p.slugVel; source = 'inputs'; } if (p.slugLen > 0) { len = p.slugLen; source = 'inputs'; } if (p.slugHoldup > 0) { hls = clamp(p.slugHoldup, 0.05, 1); source = 'inputs'; }
  const on = p.slugMode === 'on' ? true : p.slugMode === 'off' ? false : hasFs || (fs && fs.type === 'none' ? p.slugFreq > 0 : cand.length > 0 || p.slugFreq > 0);
  freq = Math.max(fin(freq, 0), 0); v = Math.max(fin(v, 1), 0.05); if (freq > 0) len = Math.min(len, (0.95 * v) / freq); len = Math.max(fin(len, 1), 0.1);
  const hlf = Math.min(nd.slug.hlf, 0.95 * hls), rhoS = p.slugMomentum > 0 ? (p.slugMomentum * 1000) / (p.Ai * v * v) : p.slugDensity > 0 ? p.slugDensity : nd.rhoL * hls + nd.rhoG * (1 - hls), rhoF = nd.rhoL * hlf + nd.rhoG * (1 - hlf);
  return { on: on && freq > 0, freq, period: freq > 0 ? 1 / freq : 0, len, v, hls, hlf, rhoS, rhoF: Math.min(rhoF, rhoS), node: nd, source, fixedRho: p.slugMomentum > 0 || p.slugDensity > 0 };
}

// ---- engine sections ----------------------------------------------------------------------------------------------------
/** Corrosion and erosion along the line and at the listed bends. */
function degradation(p, st) {
  const N = st.nodes, yCO2 = clamp(p.co2, 0, 100) / 100, hasWater = (st.fm.rates.qWaterStd || 0) > 0, mp = Math.max(p.sandKgD, 0) / 86400, dp = Math.max(p.sandUm, 1) * 1e-6, hv = (p.smts / 3.2) * 9.80665e-3; // Vickers hardness (GPa) from the tensile strength
  const corr = N.map((nd) => {
    const eta = p.eta * (1 - clamp((nd.tauW - p.tauCrit) / Math.max(p.tauCrit, 1e-6), 0, 1)), dw = deWaardMilliams({ T: nd.T, pCO2: yCO2 * nd.P, P: nd.P, model: p.corrModel, U: nd.vL, d: p.ID, pH: p.phAct, glycolWt: p.glycolWt, inhibEff: eta, mult: p.corrMult });
    const oilWet = p.wetModel === 'dewaard' && (!hasWater || (nd.wcut < 0.3 && nd.vL > 1));
    return { rate: oilWet ? 0 : dw.rate, bare: dw.base * dw.Fscale * dw.Fglycol * p.corrMult, eta, oilWet, dw };
  });
  const eroAt = (nd, rOverD) => {
    const b = { mp, U: nd.vm, D: p.ID, dp, rhoM: nd.rhoNS, muM: nd.muNS, rhoP: p.sandDensity, rOverD, gf: p.geomFactor, mult: p.erosionMult, hv, flowStress: hv * 1e9 }, r = {};
    for (const m of ['dnv', 'salama', 'oka', 'finnie']) r[m] = erosionRate({ ...b, model: m }).rate;
    r.gov = p.erosionModel === 'governing' ? Math.max(r.dnv, r.salama) : r[p.erosionModel] ?? r.dnv; r.model = p.erosionModel === 'governing' ? (r.dnv >= r.salama ? 'DNV-RP-O501' : 'Salama') : { dnv: 'DNV-RP-O501', salama: 'Salama', oka: 'Oka', finnie: 'Finnie' }[p.erosionModel] || 'DNV-RP-O501';
    return r;
  };
  const ero = N.map((nd) => { const ve = erosionalVelocity(nd.rhoNS, p.c14e); return { straight: erosionRate({ model: 'dnv', geometry: 'straight', mp, U: nd.vm, D: p.ID, mult: p.erosionMult }).rate, ve, ratio: nd.vm / ve }; });
  const fast = N.reduce((a, b) => (b.vm > a.vm ? b : a)), list = p.bends.length ? p.bends.map((b, k) => ({ ...b, name: `Bend ${k + 1}`, nd: st.at(b.x), rOverD: b.radius > 0 ? Math.max(b.radius / p.ID, 0.5) : 1.5 })) : [{ x: fast.x, angle: 90, radius: 1.5 * p.ID, name: 'Generic 1.5D elbow at the fastest point', nd: fast, rOverD: 1.5 }];
  const bends = list.map((b) => ({ ...b, ero: eroAt(b.nd, b.rOverD), corr: corr[b.nd.i].rate }));
  let gov = { rate: -1 };
  N.forEach((nd, i) => { const r = corr[i].rate + ero[i].straight; if (r > gov.rate) gov = { rate: r, x: nd.x, corr: corr[i].rate, ero: ero[i].straight, where: 'straight pipe' }; });
  for (const b of bends) { const r = b.corr + b.ero.gov; if (r > gov.rate) gov = { rate: r, x: b.nd.x, corr: b.corr, ero: b.ero.gov, where: b.name }; }
  const iC = corr.reduce((k, c, i) => (c.rate > corr[k].rate ? i : k), 0), eroMax = Math.max(...bends.map((b) => b.ero.gov), ...ero.map((e) => e.straight)), bWorst = bends.reduce((a, b) => (b.ero.gov > a.ero.gov ? b : a));
  const iV = ero.reduce((k, e, i) => (e.ratio > ero[k].ratio ? i : k), 0);
  return { corr, ero, bends, gov, iC, corrMax: corr[iC].rate, eroMax, bWorst, iV, erosionalRatio: ero[iV].ratio, yCO2, hasWater, mp, dp };
}

/** Stress, code utilisation, collapse and global buckling along the line. */
function stressAlong(p, st, wall) {
  const N = st.nodes, dnv = p.code === 'dnv', tNow = wall.tNow, tDes = Math.max(p.t - p.CA / 1000 - (dnv ? p.tFabMm / 1000 : 0), 0.05 * p.t), gLB = DNV_GM * DNV_SC[p.sc].lb, H = p.layTension * 1000;
  const head = (nd) => (p.pRefLoc === 'uniform' ? 0 : p.rhoContent * G * (st.zRef - nd.z)), pLocal = (nd, inc = 1, pRef = p.designPressure) => Math.max(pRef * BAR * inc + head(nd), PATM);
  const fmtOf = (nd, t, T) => pressureDesign({ code: p.code, D: p.D, t, smys: p.S, smts: p.Su, T, factor: nd.riser ? p.riserFactor : p.designFactor, safetyClass: p.sc });
  const fyCold = fmtOf(N[0], tNow, p.tInstall).fy, colOf = (t) => collapsePressure({ D: p.D, t, E: p.E, nu: p.nu, fy: fyCold, ovality: p.f0 }), col = colOf(tNow);
  const rows = N.map((nd) => {
    const fo = fmtOf(nd, tNow, nd.T), fd = fmtOf(nd, tDes, p.designTemp), fn = fmtOf(nd, tNow, p.designTemp), restrained = p.restraint === 'restrained' || (p.restraint === 'auto' && !nd.riser), common = { D: p.D, t: tNow, pe: nd.pe, E: p.E, nu: p.nu, alpha: p.alphaT, restrained, axial: restrained ? H : 0 };
    const sOp = pipeStress({ ...common, pi: nd.pi, dT: nd.T - p.tInstall }), pD = pLocal(nd), sDes = pipeStress({ ...common, pi: pD, dT: p.designTemp - p.tInstall }), pLi = pLocal(nd, fd.incidental);
    const seff = restrained ? H - (nd.pi - PATM) * p.Ai * (1 - 2 * p.nu) - p.As * p.E * p.alphaT * (nd.T - p.tInstall) : 0, pex = Math.max(nd.pe - PATM, 0);
    return { nd, restrained, sOp, sDes, pD, seff, hoopOp: Math.max(nd.pi - nd.pe, 0) / fo.allowDp, hoopDes: Math.max(pLi - nd.pe, 0) / fd.allowDp, longOp: Math.abs(sOp.longitudinal) / fo.longAllow, longDes: Math.abs(sDes.longitudinal) / fn.longAllow, vmOp: sOp.vonMises / fo.combAllow, vmDes: sDes.vonMises / fn.combAllow,
      collapse: (pex * gLB) / col.pc, propagation: (pex * gLB) / col.ppr, tReq: fd.tReq(pLi - nd.pe), allowRef: (t) => (fmtOf(nd, t, p.designTemp).allowDp + nd.pe - head(nd)) / fd.incidental, basis: fd.basis };
  });
  const arg = (f) => rows.reduce((k, r, i) => (f(r) > f(rows[k]) ? i : k), 0), iH = arg((r) => Math.max(r.hoopDes, r.hoopOp)), iV = arg((r) => Math.max(r.vmDes, r.vmOp)), iC = arg((r) => r.collapse), iS = arg((r) => -r.seff);
  const mawpOf = (t) => Math.max(Math.min(...rows.map((r) => r.allowRef(t))) / BAR, 0), tEol = Math.max(tNow - (wall.rate * Math.max(p.designLife - p.age, 0)) / 1000, 0.02 * p.t);
  // global buckling of the restrained flowline: Hobbs forces and hold-down against upheaval
  const rS = rows[iS], sub = rS.nd.depth > 0, wDry = (p.rhoSteel * p.As + p.coatDensity * (Math.PI / 4) * (p.Dh ** 2 - p.D ** 2) + rS.nd.rhoM * p.Ai) * G, wSub = wDry - (sub ? RHO_SW * (Math.PI / 4) * p.Dh ** 2 * G : 0), w = Math.max(wSub, 1), EI = p.E * p.Isteel, comp = Math.max(-rS.seff, 0);
  const hb = hobbs({ EI, EA: p.E * p.As, w, muA: p.muAxial, muL: p.muLateral }), up = upheavalDownload({ EI, P: comp, delta: p.imperfection, w0: w }), soil = (Hc) => p.soilGamma * 1000 * Hc * p.Dh * (1 + (0.5 * Hc) / p.Dh), resist = w + soil(p.coverDepth);
  const coverReq = up.wReq <= w ? 0 : soil(10) + w < up.wReq ? 10 : brent((Hc) => w + soil(Hc) - up.wReq, 0, 10, 1e-9);
  return { rows, iH, iV, iC, iS, col, colOf, tDes, tEol, pLocal, head, fmtOf, gLB, mawp: mawpOf(tNow), mawpEol: mawpOf(tEol), mawpOf, minWall: (Math.max(...rows.map((r) => r.tReq)) + p.CA / 1000 + (dnv ? p.tFabMm / 1000 : 0)) * 1000,
    hoopUtil: Math.max(rows[iH].hoopDes, rows[iH].hoopOp), vmUtil: Math.max(rows[iV].vmDes, rows[iV].vmOp), longUtil: Math.max(...rows.map((r) => Math.max(r.longOp, r.longDes))), collapseUtil: rows[iC].collapse, propUtil: rows[iC].propagation,
    buckling: { hb, up, wSub, comp, resist, coverReq, hobbsUtil: comp / hb.critical, palmerUtil: up.wReq / resist, buried: p.coverDepth > 0, x: rS.nd.x, EI } };
}

/** Train of slugs with log-normal scatter of the period and of the body length (seeded, so repeatable). */
function slugTrain(p, sg) {
  const r = rng(Math.round(p.seed) || 11), train = []; let t0 = 0;
  if (sg.on && sg.freq > 0) for (let i = 0; i < p.nSlugs; i++) { const per = Math.exp(r.normal(-0.02, 0.2)) / sg.freq, len = Math.min(sg.len * Math.exp(r.normal(-0.045, 0.3)), 0.95 * sg.v * per); train.push({ t0, len, v: sg.v, period: per }); t0 += per; }
  return train;
}
/** Slug loads on bends, pressure pulsation and the slug-train force history. */
function slugLoads(p, st, ds, ctx) {
  const A = p.Ai, nd = ds.node, k = ds.v / Math.max(nd.slug.v, 1e-6), vmS = Math.max(nd.vm * k, 0.01), fromKernel = ds.source === 'kernel slug closures';
  const f90 = ds.on ? bendForce({ rho: ds.rhoS, A, v: ds.v, angle: 90, dlf: p.dlf }).force : 0, steady90 = Math.max(...st.nodes.map((n) => bendForce({ rho: n.rhoNS, A, v: n.vm, angle: 90 }).force));
  const dpAcc = ds.on ? nd.rhoL * ds.hls * (ds.v - vmS) ** 2 * Math.max(ds.hls / Math.max(ds.hlf, 1e-3) - 1, 0) : 0, fD = frictionFactor(Math.max((ds.rhoS * vmS * p.ID) / nd.muL, 10), (p.roughUm * 1e-6) / p.ID), dpFric = ds.on ? ((fD * ds.rhoS * vmS * vmS) / (2 * p.ID)) * ds.len : 0;
  const bends = (p.bends.length ? p.bends : [{ x: nd.x, angle: 90, radius: 1.5 * p.ID }]).map((b, i) => {
    const n = st.at(b.x), rho = !ds.on ? n.rhoNS : ds.fixedRho ? ds.rhoS : fromKernel ? n.slug.rhoS : n.rhoL * ds.hls + n.rhoG * (1 - ds.hls), v = !ds.on ? n.vm : fromKernel ? n.slug.v : (ds.v * n.vm) / Math.max(nd.vm, 1e-6);
    const bf = bendForce({ rho, A, v, angle: b.angle, dlf: ds.on ? p.dlf : 1, p: n.pi - n.pe, radius: b.radius }), st0 = bendForce({ rho: n.rhoNS, A, v: n.vm, angle: b.angle }), fiv = fivScreen({ rho: n.rhoNS, v: n.vm, Dmm: p.D * 1000, tmm: p.wtMm, support: p.pipingSupport });
    return { name: `Bend ${i + 1}`, x: n.x, angle: b.angle, radius: b.radius, rho, v, impulse: ds.on ? ((bf.force / p.dlf) * ds.len) / Math.max(v, 1e-6) : 0, force: bf.force, steady: st0.force, pressure: bf.pressureForce, centrifugal: bf.centrifugal, fiv };
  });
  const train = slugTrain(p, ds);
  const Fs = ds.rhoS * A * ds.v ** 2 * Math.SQRT2, Ff = ds.rhoF * A * ds.v ** 2 * Math.SQRT2, tr = Math.max(p.ID / ds.v, 1e-3), ht = [0], hf = [ds.on ? Ff : steady90];
  for (const s of train.slice(0, 12)) { const b = Math.max(s.len / s.v, 1.5 * tr); ht.push(s.t0 + 0.2 * s.period, s.t0 + 0.2 * s.period + tr, s.t0 + 0.2 * s.period + b, s.t0 + 0.2 * s.period + b + tr); hf.push(Ff, Fs, Fs, Ff); }
  if (!ds.on) { ht.push(60); hf.push(steady90); }
  const fPIn = ctx.outputs?.flow?.pInAmplitude;
  return { f90, steady90, peak90: ds.on ? f90 : steady90, dpAcc, dpFric, pulsation: (dpAcc + dpFric) / 2 / BAR, runForce: (dpAcc + dpFric) * A, bends, train, history: { t: ht, f: hf.map((x) => x / 1000) }, flowAmplitude: isNum(fPIn) ? fPIn : null, Fs, Ff };
}

/** Free-span or jumper dynamics: modes, static sag, slug-train response, vortex-induced vibration and the fatigue spectrum. */
function spanDynamics(p, st, ds, sl, sa, ctx) {
  const nd = st.at(p.spanX), sub = nd.depth > 0, L = p.spanLength, areaO = (Math.PI / 4) * p.Dh ** 2, mSteel = p.rhoSteel * p.As, mCoat = p.coatDensity * (Math.PI / 4) * (p.Dh ** 2 - p.D ** 2), mCont = nd.rhoM * p.Ai, mAdd = sub ? p.addedMass * RHO_SW * areaO : 0, me = mSteel + mCoat + mCont + mAdd;
  const wSub = (mSteel + mCoat + mCont) * G - (sub ? RHO_SW * areaO * G : 0), EI = p.E * p.Isteel, seffRaw = sa.rows[nd.i].seff, seff = sa.buckling.buried ? seffRaw : Math.max(seffRaw, -sa.buckling.hb.critical), Naxial = p.spanAxial ? seff : 0, bendX = p.bendAt > 0 ? clamp(p.bendAt, 0.02, 0.98) * L : null, ro = p.D / 2;
  const model = beamModel({ L, EI, m: me, n: p.nElem, ends: endsOf(p.spanEnds), kT: [p.spanKt * 1000, p.spanKt * 1000], kR: [p.spanKr * 1000, p.spanKr * 1000], supports: p.midSupports, nodesAt: bendX ? [bendX] : [], N: Naxial, kGA: p.timoshenko ? (0.5 * p.E * p.As) / (2 * (1 + p.nu)) : 0, rhoI: p.timoshenko ? p.rhoSteel * p.Isteel : 0 });
  const modes = beamModes(model, p.nModes), f1 = modes.f[0], stat = beamStatic(model, { q: Math.abs(wSub) }), sigStatic = (stat.maxMoment * ro) / p.Isteel / MPA, euler = eulerLoad(EI, L, END_K[p.spanEnds] || 1);
  // slug-train response
  const sg = ds.source === 'kernel slug closures' ? { on: ds.on && (nd.slugLike || p.slugMode === 'on') && nd.slug.freq > 0, freq: nd.slug.freq, len: Math.max(nd.slug.len, 0.1), v: Math.max(nd.slug.v, 0.05), hls: nd.slug.hls, hlf: nd.slug.hlf } : { on: ds.on, freq: ds.freq, len: ds.len, v: ds.v, hls: ds.hls, hlf: ds.hlf };
  const train = slugTrain(p, sg), hasSlug = sg.on && train.length > 0, v = sg.v, rhoS = ds.fixedRho ? ds.rhoS : nd.rhoL * sg.hls + nd.rhoG * (1 - sg.hls), rhoF = Math.min(nd.rhoL * sg.hlf + nd.rhoG * (1 - sg.hlf), rhoS), dw = (rhoS - rhoF) * p.Ai * G;
  let resp = null, cycles = [], perYear = 0, note = '', dt = 1 / (f1 * p.spc), simSlugs = 0;
  if (hasSlug) {
    dt = Math.min(dt, L / v / 25);
    const decay = Math.min(4 / (p.zeta * modes.omega[0]), 40 / f1), slugs = []; let t0 = 0.5 / f1, real = 0;
    for (const s of train) { const gap = Math.min(s.period, (L + s.len) / v + decay); if ((t0 + gap) / dt > 60000 && slugs.length >= 2) break; slugs.push({ t0, len: s.len, v }); t0 += gap; real += s.period; }
    simSlugs = slugs.length; if (simSlugs < train.length) note = `Only ${simSlugs} of ${train.length} slugs were simulated to keep the response history below 60,000 steps.`;
    resp = slugResponse({ model, modes, zeta: p.zeta, slugs, dw, bend: bendX ? { x: bendX, dF: (rhoS - rhoF) * p.Ai * v * v * 2 * Math.sin((p.bendAngle * Math.PI) / 360) } : null, dt, tEnd: t0, ro, E: p.E, ramp: Math.max(p.ID, 2 * v * dt) });
    cycles = rainflow(resp.sigma.map((s) => s / MPA)); perYear = YEAR / real;
  }
  const snOpt = { scf: p.scf, t: p.wtMm }, dmgOf = (cyc) => snDamage(cyc, p.snClass, p.snEnv, snOpt).damage, spectrum = [], dSlug = dmgOf(cycles) * perYear;
  if (cycles.length) { // bin the counted ranges for the tables and the S–N plot (damage itself uses every counted cycle)
    const rmax = Math.max(...cycles.map((c) => c.range)), nb = 14, bins = Array.from({ length: nb }, () => ({ n: 0, s3: 0 }));
    for (const c of cycles) { const b = bins[Math.min(nb - 1, Math.floor((c.range / (rmax || 1)) * nb))]; b.n += c.count; b.s3 += c.count * c.range ** 3; }
    bins.forEach((b) => { if (b.n > 0) spectrum.push({ source: 'Slug passage', range: (b.s3 / b.n) ** (1 / 3), perYear: b.n * perYear }); });
  }
  // vortex-induced vibration over a Rayleigh long-term current distribution (mean = the input current), plus the wave-induced velocity at the seabed
  let uw = 0;
  if (sub && p.waveHeight > 0 && p.wavePeriod > 0 && nd.depth < 400) { const om = (2 * Math.PI) / p.wavePeriod; let kk = (om * om) / G; for (let i = 0; i < 40; i++) kk = (om * om) / (G * Math.tanh(kk * nd.depth)); uw = (Math.PI * p.waveHeight) / (p.wavePeriod * Math.sinh(Math.min(kk * nd.depth, 50))); }
  const vivAt = (U) => vivScreen({ U, f1, D: p.Dh, me, zeta: p.zeta }), unit = (p.E * ro * modes.unitCurv[0]) / MPA, viv = sub && p.currentSpeed + uw > 0 ? vivAt(p.currentSpeed + uw) : vivAt(0); let dViv = 0, vivWorst = { aIL: 0, aCF: 0, U: 0 };
  if (sub && p.currentSpeed > 0) {
    const scale = p.currentSpeed / 0.886226925, nbin = 24, umax = 3.2 * p.currentSpeed, cdf = (u) => 1 - Math.exp(-((u / scale) ** 2));
    for (let i = 0; i < nbin; i++) {
      const a = (umax * i) / nbin, b = (umax * (i + 1)) / nbin, pr = cdf(b) - cdf(a), U = 0.5 * (a + b) + uw, r = vivAt(U), n = f1 * YEAR * pr;
      for (const [amp, mult, label] of [[r.aIL, 1, 'In-line VIV'], [r.aCF, 1, 'Cross-flow VIV']]) if (amp > 0) { const range = 2 * amp * p.Dh * unit * mult; dViv += n / snCycles(range, p.snClass, p.snEnv, snOpt); spectrum.push({ source: label, range, perYear: n }); }
      if (r.aCF + r.aIL > vivWorst.aCF + vivWorst.aIL) vivWorst = { ...r, U };
    }
  }
  // start-up / shutdown and blowdown cycles of the restrained line, and the pressure cycles of the flow suite's transient
  const hot = st.nodes.reduce((a, b) => (b.T > a.T ? b : a)), restr = p.restraint !== 'unrestrained', eA = (p.E * p.alphaT) / MPA, hoopOf = (dpPa) => (dpPa * (p.D - p.t)) / (2 * p.t) / MPA;
  const sShut = restr ? eA * Math.max(hot.T - p.tInstall, 0) : 0, sBlow = (restr ? eA * Math.max(hot.T - Math.min(p.blowdownMinT, p.tInstall), 0) : 0) + (restr ? p.nu : 0.5) * hoopOf(Math.max(hot.pi - PATM, 0)); let dOps = 0, dPuls = 0;
  for (const [n, range, source] of [[p.eventsShutdown, sShut, 'Shutdown / restart'], [p.eventsBlowdown, sBlow, 'Blowdown']]) if (n > 0 && range > 0) { dOps += n / snCycles(range, p.snClass, p.snEnv, snOpt); spectrum.push({ source, range, perYear: n }); }
  const ser = ctx.outputs?.flow?.series;
  if (ser && Array.isArray(ser.t) && Array.isArray(ser.pIn) && ser.t.length > 4 && ser.pIn.length === ser.t.length && ser.pIn.every(isNum) && ser.t[ser.t.length - 1] > ser.t[0]) {
    const cyc = rainflow(ser.pIn).map((c) => ({ range: hoopOf(c.range * BAR), count: c.count })).filter((c) => c.range > 0.05), sc = YEAR / (ser.t[ser.t.length - 1] - ser.t[0]);
    if (cyc.length) { const d = snDamage(cyc, p.snClass, p.snEnv, snOpt); dPuls = d.damage * sc; spectrum.push({ source: 'Inlet-pressure cycles (flow transient)', range: d.sEq, perYear: d.cycles * sc }); }
  }
  if (p.opLog.length > 4 && p.opLog[p.opLog.length - 1].t > p.opLog[0].t) { // operating log: pressure cycles as hoop stress, temperature cycles as restrained thermal stress
    const sc = YEAR / ((p.opLog[p.opLog.length - 1].t - p.opLog[0].t) * 3600), sets = [['Operating log: pressure cycles', rainflow(p.opLog.map((r) => r.p)).map((c) => ({ range: hoopOf(c.range * BAR), count: c.count }))], ['Operating log: temperature cycles', restr ? rainflow(p.opLog.map((r) => r.T).filter(isNum)).map((c) => ({ range: eA * c.range, count: c.count })) : []]];
    for (const [source, cyc0] of sets) { const cyc = cyc0.filter((c) => c.range > 0.05); if (cyc.length) { const d = snDamage(cyc, p.snClass, p.snEnv, snOpt); dPuls += d.damage * sc; spectrum.push({ source, range: d.sEq, perYear: d.cycles * sc }); } }
  }
  const dYear = dSlug + dViv + dOps + dPuls, excite = [sg.on ? { f: sg.freq, name: 'slug frequency' } : null, sub && viv.fShed > 0 ? { f: viv.fShed, name: 'vortex shedding' } : null].filter(Boolean);
  const gov = excite.length ? excite.reduce((a, b) => (Math.abs(Math.log(b.f / f1)) < Math.abs(Math.log(a.f / f1)) ? b : a)) : { f: 0, name: 'none' };
  return { nd, sub, L, me, wSub, EI, seff, seffRaw, Naxial, model, modes, f1, stat, sigStatic, euler, eulerUtil: Math.max(-seff, 0) / euler, resp, cycles, perYear, spectrum, dSlug, dViv, dOps, dPuls, dYear, viv, vivWorst, uw, unit, fivRatio: gov.f / f1, excitation: gov, note, dt, simSlugs, dynMax: resp ? resp.sigmaMax / MPA : 0, dispMax: resp ? resp.dispMax : 0, passage: hasSlug ? (L + sg.len) / v : 0, slug: { ...sg, rhoS, rhoF } };
}

/** Limit states, FORM / SORM / sampling comparison and the reliability index over time. */
async function reliability(p, st, wall, sa, deg, span, defectWorst, ctx) {
  const tN = p.wtMm, Dmm = p.D * 1000, ratioU = p.smts / p.smys, cW = p.covWt / 100, cY = p.covYield / 100, cP = p.covPress / 100, cC = p.covCorr / 100, age = p.age, life = Math.max(p.designLife, age + 1);
  const shut = sa.rows.map((r) => (Math.max(p.pShutIn * BAR + sa.head(r.nd), r.nd.pi) - r.nd.pe) / MPA), iB = shut.reduce((k, x, i) => (x > shut[k] ? i : k), 0), pMean = Math.max(shut[iB], 0.01), peMax = Math.max(...st.nodes.map((n) => n.pe - PATM)) / MPA;
  const vT = { name: 'Wall thickness', dist: 'normal', mean: tN, cov: cW }, vY = { name: 'Yield strength', dist: 'lognormal', mean: p.yieldBias * p.smys, cov: cY }, vP = { name: 'Annual extreme pressure', dist: 'gumbel', mean: pMean, cov: cP }, vC = { name: 'Wall-loss model', dist: 'lognormal', mean: 1, cov: cC };
  const lossOf = (rate, loss0) => (T, xc) => loss0 + xc * rate * Math.max(T - age, 0), lossGen = lossOf(wall.rate, wall.lossNow), eroLoss = lossOf(deg.bWorst.ero.gov, Math.min(deg.bWorst.ero.gov * age, 0.9 * tN)), fyCold = sa.fmtOf(st.nodes[0], p.t, p.tInstall).fy / p.S;
  const LS = [
    { key: 'burst', name: 'Burst of the corroding wall at the annual extreme pressure', annual: true, vars: [vT, vY, vP, vC, { name: 'Burst model', dist: 'normal', mean: 1, cov: 0.05 }], g: (T) => (x) => { const tw = Math.max(x[0] - lossGen(T, x[3]), 0.01); return (x[4] * ((2 * tw) / (Dmm - tw)) * Math.min(x[1], (x[1] * ratioU) / 1.15) * 1.1547005 - x[2]) / pMean; } },
    { key: 'leak', name: 'Corrosion: wall loss reaches the leak criterion', annual: false, vars: [vT, vC], g: (T) => (x) => (p.leakFrac * x[0] - lossOf(deg.corrMax, wall.lossNow)(T, x[1])) / tN },
    { key: 'erosion', name: 'Erosion at the worst bend reaches the leak criterion', annual: false, vars: [vT, { ...vC, name: 'Erosion model', cov: Math.max(cC, 1) }], g: (T) => (x) => (p.leakFrac * x[0] - eroLoss(T, x[1])) / tN },
    { key: 'fatigue', name: 'Fatigue: Miner sum reaches its resistance', annual: false, vars: [{ name: 'Miner resistance', dist: 'lognormal', mean: 1, cov: 0.3 }, { name: 'Stress / S–N model', dist: 'lognormal', mean: 1, cov: 0.6 }], g: (T) => (x) => x[0] - (p.priorDamage + span.dYear * Math.max(T - age, 0)) * x[1] },
  ];
  if (peMax > 0.01) LS.push({ key: 'collapse', name: 'External collapse of the depressurised pipe', annual: true, vars: [vT, vY, vC, { name: 'Collapse model', dist: 'normal', mean: 1, cov: 0.08 }, { name: 'Out-of-roundness', dist: 'lognormal', mean: p.f0, cov: 0.25 }], g: (T) => (x) => { const tw = Math.max(x[0] - lossGen(T, x[2]), 0.05); return (x[3] * collapsePressure({ D: Dmm, t: tw, E: p.eMod, nu: p.nu, fy: x[1] * fyCold, ovality: x[4] }).pc) / peMax - 1; } });
  if (defectWorst) LS.push({ key: 'defect', name: `Burst at the worst measured defect (x = ${fmt(defectWorst.x)} m)`, annual: true, vars: [{ name: 'Defect depth sizing', dist: 'normal', mean: 1, cov: 0.1 }, vC, { ...vP, mean: Math.max(defectWorst.dpExt / MPA, 0.01) }, { name: 'Tensile strength', dist: 'lognormal', mean: 1.09 * p.smts, cov: cY }, { name: 'Burst model', dist: 'normal', mean: 1.05, cov: 0.1 }],
    g: (T) => (x) => (x[4] * b31g({ D: Dmm, t: tN, d: Math.min(defectWorst.depth * x[0] + x[1] * defectWorst.rate * Math.max(T - age, 0), 0.999 * tN), L: defectWorst.length, smys: p.smys, smts: x[3], method: 'dnv' }).pf - x[2]) / Math.max(defectWorst.dpExt / MPA, 0.01) });
  const Tg = linspace(age, life, 11), clampB = (b) => clamp(fin(b, 10), -5, 10);
  for (const ls of LS) {
    ls.beta = []; ls.pf = [];
    for (const T of [...Tg, age + 1]) { let r; try { r = form(ls.g(T), ls.vars); } catch { r = { beta: 10 }; } const b = clampB(r.beta); ls.beta.push(b); ls.pf.push(Phi(-b)); }
    const pNext = ls.pf.pop(); ls.beta.pop();
    ls.annualNow = ls.annual ? ls.pf[0] : Math.max(pNext - ls.pf[0], 0) / Math.max(1 - ls.pf[0], 1e-12);
    const dT = Tg[Tg.length - 1] - Tg[Tg.length - 2]; ls.annualEol = ls.annual ? ls.pf[Tg.length - 1] : Math.max(ls.pf[Tg.length - 1] - ls.pf[Tg.length - 2], 0) / dT;
    ls.annualAvg = 0;
    ls.annualT = Tg.map((_, i) => (ls.annual ? ls.pf[i] : i === 0 ? ls.annualNow : Math.max(ls.pf[i] - ls.pf[i - 1], 0) / (Tg[i] - Tg[i - 1])));
  }
  for (const ls of LS) { let a = 0; for (let i = 1; i < Tg.length; i++) a += 0.5 * (ls.annualT[i] + ls.annualT[i - 1]) * (Tg[i] - Tg[i - 1]); ls.annualAvg = clamp(a / (Tg[Tg.length - 1] - Tg[0]), 0, 1); }
  if (ctx.tick) await ctx.tick();
  const series = (f) => 1 - LS.reduce((q, ls) => q * (1 - clamp(f(ls), 0, 1)), 1), pofNow = series((ls) => ls.annualNow), pofEol = series((ls) => ls.annualEol), pofAvg = series((ls) => ls.annualAvg), pofT = Tg.map((_, i) => series((ls) => ls.annualT[i])), target = DNV_SC[p.sc].target;
  let tTarget = null; for (let i = 0; i < Tg.length; i++) if (pofT[i] > target) { tTarget = i === 0 ? Tg[0] : Tg[i - 1] + ((Tg[i] - Tg[i - 1]) * (Math.log(target) - Math.log(Math.max(pofT[i - 1], 1e-300)))) / (Math.log(pofT[i]) - Math.log(Math.max(pofT[i - 1], 1e-300)) || 1); break; }
  // method comparison on the limit state that governs at the end of the design life
  const gov = LS.reduce((a, b) => (b.beta[Tg.length - 1] < a.beta[Tg.length - 1] ? b : a)), g = gov.g(life), f = form(g, gov.vars), so = sorm(g, gov.vars, f), seed = Math.round(p.seed) || 11;
  const mc = monteCarlo(g, gov.vars, { n: p.nMC, seed, keep: true }); if (ctx.tick) await ctx.tick();
  const is = monteCarlo(g, gov.vars, { n: Math.max(1000, Math.round(p.nMC / 4)), method: 'is', center: f.u, seed: seed + 1 }), lh = monteCarlo(g, gov.vars, { n: Math.min(p.nMC, 4000), method: 'lhs', seed: seed + 2 }), rs = responseSurface(g, gov.vars, { seed: seed + 3 });
  return { LS, Tg, pofNow, pofEol, pofAvg, pofT, target, tTarget, gov, f, so, mc, is, lh, rs, life, iB, pMean, peMax };
}

const cell = (x, n = 4) => (typeof x === 'number' ? (Number.isFinite(x) ? sig(x, n) : '—') : x ?? '—');
const stat = (u, warn = 0.8, bad = 1) => (u > bad ? 'bad' : u > warn ? 'warn' : 'ok');
const thin = (a, n = 300) => { if (a.length <= n) return a.slice(); const s = (a.length - 1) / (n - 1); return Array.from({ length: n }, (_, i) => a[Math.round(i * s)]); };

/** Assessment of the measured defects (and of the deepest feature of an imported wall map) by B31G, Modified B31G and DNV-RP-F101. */
function assessDefects(p, st, sa, deg, wall, mapS) {
  const list = p.defects.map((d, i) => ({ ...d, name: `Defect ${i + 1}` })); if (mapS && mapS.depthMm > 0.02 * p.wtMm) list.push({ name: 'Deepest feature of the wall map', x: mapS.x, depth: mapS.depthMm, length: mapS.lengthMm });
  const fOf = (nd) => (p.code === 'dnv' ? 0.9 * 0.72 : nd.riser ? p.riserFactor : p.designFactor), base = { D: p.D, t: p.t, smys: p.S, smts: p.Su };
  const rows = list.map((d) => {
    const nd = st.at(d.x), dp = Math.max(sa.pLocal(nd) - nd.pe, 0), F = fOf(nd), dm = Math.min(d.depth, p.wtMm) / 1000, Lm = d.length / 1000, r = {};
    for (const m of ['b31g', 'modified', 'dnv']) r[m] = b31g({ ...base, d: dm, L: Lm, method: m }).pf;
    const safe = F * r.modified, erf = safe > 0 ? dp / safe : 99, rate = Math.max(deg.corr[nd.i].rate, wall.rate), erfAt = (y) => { const dd = Math.min(dm + (rate * y) / 1000, 0.999 * p.t), s = F * b31g({ ...base, d: dd, L: Lm, method: 'modified' }).pf; return s > 0 ? dp / s : 99; };
    const life = erf >= 1 || d.depth / p.wtMm > 0.8 ? 0 : !(rate > 1e-9) || erfAt(200) < 1 ? 200 : brent((y) => erfAt(y) - 1, 0, 200, 1e-6), lifeDepth = rate > 1e-9 ? Math.max(0.8 * p.wtMm - d.depth, 0) / rate : 200;
    return { ...d, nd, dp, F, ratio: d.depth / p.wtMm, pf: r, safe, erf, ok: erf <= 1 && d.depth / p.wtMm <= 0.8, life: Math.min(life, lifeDepth, 200), rate, dpExt: Math.max(Math.max(p.pShutIn * BAR + sa.head(nd), nd.pi) - nd.pe, 0) };
  });
  // acceptance curves at the largest design differential pressure of the line
  const ndG = sa.rows.reduce((a, r) => (sa.pLocal(r.nd) - r.nd.pe > sa.pLocal(a.nd) - a.nd.pe ? r : a)).nd, dpG = Math.max(sa.pLocal(ndG) - ndG.pe, 0), FG = fOf(ndG), Ls = logspace(10, 3000, 36);
  const curve = (m) => Ls.map((L) => { const f = (r) => FG * b31g({ ...base, d: r * p.t, L: L / 1000, method: m }).pf - dpG; return f(0.8) >= 0 ? 0.8 : f(0) <= 0 ? 0 : brent(f, 0, 0.8, 1e-9); });
  return { rows, Ls, curves: { b31g: curve('b31g'), modified: curve('modified'), dnv: curve('dnv') }, worst: rows.length ? rows.reduce((a, b) => (b.erf > a.erf ? b : a)) : null, dpG };
}

async function run(v0, ctx = {}) {
  const prog = (f, m) => { try { ctx.progress?.(f, m); } catch { /* progress is optional */ } }, tick = async () => { if (ctx.tick) await ctx.tick(); };
  const p = params(v0); prog(0.03, 'Reading the flow picture of the case');
  const st = lineState(p, ctx), N = st.nodes, xkm = N.map((n) => n.x / 1000), warnings = [], recs = [], tables = [], plots = [], viol = [];
  await tick();
  // ---- degradation and wall state
  const deg = degradation(p, st), mapS = wallMapStats(p.wtMap, p.wtMm), meas = [p.minWt > 0 ? p.wtMm - p.minWt : null, mapS ? mapS.depthMm : null].filter(isNum), lossMeas = meas.length ? clamp(Math.max(...meas), 0, 0.95 * p.wtMm) : null;
  const cC = p.covCorr / 100, upd = lossMeas !== null && p.age > 0 ? { loss: lossMeas, years: p.age, sd: 0.3 } : null, r0 = remainingLife({ margin: 1, rate: Math.max(deg.gov.rate, 1e-9), cov: cC, measured: upd }), rate = upd ? r0.rate : deg.gov.rate;
  const lossNow = lossMeas ?? clamp(deg.gov.rate * p.age, 0, 0.95 * p.wtMm), margin = Math.max(p.CA - lossNow, 0), rul = remainingLife({ margin, rate: Math.max(rate, 1e-9), cov: cC, measured: upd }), remLife = cap(margin / Math.max(rate, 1e-9), 200), wall = { tNow: (p.wtMm - lossNow) / 1000, rate, lossNow, margin, measured: lossMeas !== null };
  // ---- stress, slug loads, span dynamics
  prog(0.15, 'Pipe stress and buckling'); const sa = stressAlong(p, st, wall), ds = designSlug(p, st, ctx), sl = slugLoads(p, st, ds, ctx);
  prog(0.3, 'Span modes and slug response'); await tick();
  const span = spanDynamics(p, st, ds, sl, sa, ctx), remYears = Math.max(p.designLife - p.age, 0);
  const fatLife = cap(Math.max(1 - p.priorDamage, 0) / Math.max(span.dYear * p.dff, 1e-12), 1e4), lifeMin = cap(Math.max((wall.tNow - Math.max(...sa.rows.map((r) => r.tReq))) * 1000, 0) / Math.max(rate, 1e-9), 200);
  // local buckling at the span under pressure, effective force and static + dynamic bending (functional load factor 1.2)
  const mSpan = ((span.sigStatic + span.dynMax) * MPA * p.Isteel) / (p.D / 2), lbOp = localBuckling({ D: p.D, t: wall.tNow, fy: sa.fmtOf(span.nd, wall.tNow, span.nd.T).fy, fu: sa.fmtOf(span.nd, wall.tNow, span.nd.T).fu, E: p.E, nu: p.nu, pi: span.nd.pi, pe: span.nd.pe, pmin: PATM, M: 1.2 * mSpan, S: 1.2 * span.seff, safetyClass: p.sc, ovality: p.f0 });
  // ---- fracture mechanics on the fatigue spectrum
  const cycYr = sum(span.spectrum.map((s) => s.perYear)), dSeq = cycYr > 0 ? (sum(span.spectrum.map((s) => s.perYear * s.range ** p.parisM)) / cycYr) ** (1 / p.parisM) * p.scf : 0, tw = wall.tNow, a0 = Math.min(p.flawMm / 1000, 0.5 * tw);
  const sMaxT = Math.max(...sa.rows.map((r) => r.sOp.longitudinal), 0) / MPA + (span.sigStatic + span.dynMax) * p.scf, kOf = (a) => edgeCrackY(a, tw) * sMaxT * Math.sqrt(Math.PI * a) - p.kMat, aC = sMaxT > 0 && kOf(0.8 * tw) > 0 ? (kOf(a0) >= 0 ? a0 : brent(kOf, a0, 0.8 * tw, 1e-9)) : 0.8 * tw;
  const crack = parisLife({ a0, ac: Math.max(aC, a0), C: p.parisC, m: p.parisM, dS: dSeq, Y: (a) => edgeCrackY(a, tw), dKth: p.dKth }), crackLife = cycYr > 0 && Number.isFinite(crack.N) ? cap(crack.N / cycYr, 1e4) : 1e4;
  // ---- defects, blockage
  const def = assessDefects(p, st, sa, deg, wall, mapS), beta = clamp(Math.max(p.blockage, p.effIdMm > 0 ? 1 - (p.effIdMm / p.idMm) ** 2 : 0), 0, 1), ndB = st.at(p.plugX), plugged = beta >= 0.999;
  const dpRestr = plugged ? 0 : 0.5 * ndB.rhoM * ndB.vm ** 2 * (1 / (1 - beta) - 1) ** 2 + ((0.02 * ndB.rhoM * ndB.vm ** 2) / (2 * p.ID)) * p.plugLen * ((1 - beta) ** -2.5 - 1), pUp = plugged ? p.pShutIn : N[0].P + dpRestr / BAR;
  const mawpInlet = sa.mawp + (p.pRefLoc === 'top' ? (p.rhoContent * G * (st.zTop - N[0].z)) / BAR : 0), plugDp = Math.max(p.pShutIn - Math.min(ndB.P, N[N.length - 1].P), 0) * BAR, plugMass = (917 * (1 - p.plugPorosity) + 1000 * p.plugPorosity) * p.Ai * p.plugLen;
  const runLen = Math.max(Math.min(...p.bends.map((b) => b.x - ndB.x).filter((d) => d > 1), N[N.length - 1].x - ndB.x, 1000), 10), plugV = Math.min(Math.sqrt((2 * plugDp * p.Ai * runLen) / Math.max(plugMass, 1)), 300), plugE = 0.5 * plugMass * plugV ** 2;
  // ---- electrochemistry, sour service, top-of-line condensation
  const ndC = N[deg.iC], cw = deg.corr[deg.iC].dw, mix = mixedPotential({ T: ndC.T, pCO2: Math.max(cw.fCO2, 1e-6), pH: p.phAct > 0 ? p.phAct : cw.pHco2, U: ndC.vL, d: p.ID, rho: ndC.rhoL, mu: ndC.muL });
  const pH2S = Math.max(...N.map((n) => (clamp(p.h2s, 0, 100) / 100) * n.P * 100)), sour = sourRegion(pH2S, p.phAct > 0 ? p.phAct : cw.pHco2), qGs = (st.fm.rates.qGasStd || 0) / 86400;
  const tlc = N.map((n, i) => { const j = Math.min(i + 1, N.length - 1), k = Math.max(i - 1, 0), dx = Math.max(N[j].x - N[k].x, 1e-6); return condensationRate({ P: n.P, T: n.T, dTdx: (N[j].T - N[k].T) / dx, qGasStd: qGs, D: p.ID }); }), tlcMax = Math.max(...tlc);
  // ---- reliability
  prog(0.55, 'Reliability: FORM, SORM and sampling');
  const rel = await reliability(p, st, wall, sa, deg, span, def.worst, ctx), LSof = (k) => rel.LS.find((l) => l.key === k), pMiss = 1 - clamp(p.pod, 0, 1), pA = (k) => { const l = LSof(k); return l ? clamp(Math.min(l.annualT[0], l.annualAvg) + Math.max(l.annualAvg - l.annualT[0], 0) * pMiss, 0, 1) : 0; }, pRaw = (k) => LSof(k)?.annualAvg ?? 0;
  prog(0.8, 'Risk models'); await tick();
  // ---- risk: consequence, fault tree, event tree, matrix, Markov, Bayesian network
  const prodDay = (st.fm.rates.qOilStd || 0) * 6.2898 * p.oilPrice, cBase = p.repairCost * 1e6 + p.downtimeDays * prodDay, betaT = -PhiInv(rel.target), utilPof = (u) => (u > 0 ? Phi(-(betaT - Math.log(u) / 0.15)) : 0);
  const pLB = utilPof(lbOp.util), pUph = sa.buckling.buried ? utilPof(sa.buckling.palmerUtil) : 0, pCol = pA('collapse'), pPlugLoc = clamp(p.plugProb, 0, 1) * clamp(p.pfdProtect, 0, 1);
  const tree = { name: 'Loss of containment', gate: 'OR', children: [
    { name: 'Undetected corrosion failure', gate: 'AND', children: [{ name: 'Corrosion wall loss reaches the limit', p: pRaw('leak') }, { name: 'Inspection misses the damage', p: pMiss }] },
    { name: 'Undetected erosion failure', gate: 'AND', children: [{ name: 'Erosion wall loss reaches the limit', p: pRaw('erosion') }, { name: 'Inspection misses the damage', p: pMiss }] },
    { name: 'Fatigue crack through the wall', p: pA('fatigue') }, { name: 'Burst at the annual extreme pressure', p: Math.max(pA('burst'), pA('defect')) }, { name: 'Collapse or local buckling', p: clamp(pCol + pLB + pUph, 0, 1) },
    { name: 'Plug-related failure', gate: 'AND', children: [{ name: 'Hydrate or deposit plug forms', p: clamp(p.plugProb, 0, 1) }, { name: 'Plug remediation barrier fails', p: clamp(p.pfdProtect, 0, 1) }] },
    { name: 'External damage', p: clamp(p.extFreq, 0, 1) }] };
  const ft = faultTree(tree), et = eventTree(ft.top, [{ name: 'Leak detected', p: p.pDetect }, { name: 'Isolation succeeds', p: p.pIsolate }, { name: 'Ignition', p: p.pIgnite }]);
  const etRows = et.outcomes.map((o) => { const [d, iso, ign] = o.path, rel0 = d && iso ? 0.05 : 1, cost = cBase + rel0 * p.envCost * 1e6 + (ign ? p.safetyCost * 1e6 : 0); return { ...o, cost, outcome: `${d ? (iso ? 'Small release, isolated' : 'Release continues until depressurised') : 'Late discovery, full release'}${ign ? ', ignited' : ''}` }; });
  const cLoc = sum(etRows.map((o) => o.p * o.cost)), cPlug = p.plugDays * prodDay + 0.3 * p.repairCost * 1e6;
  const threats = [
    { name: 'CO2 corrosion', pof: pRaw('leak') * pMiss, cof: cLoc, det: 4, cause: `${fmt(deg.corrMax, 3)} mm/y at ${fmt(ndC.x / 1000, 3)} km`, effect: 'Pinhole leak, then rupture', action: 'Inhibition availability, inline inspection' },
    { name: 'Sand erosion', pof: pRaw('erosion') * pMiss, cof: cLoc, det: 5, cause: `${fmt(deg.bWorst.ero.gov, 3)} mm/y at ${deg.bWorst.name}`, effect: 'Wall thinning at fittings', action: 'Sand monitoring, velocity limit, UT at bends' },
    { name: 'Slugging and VIV fatigue', pof: pA('fatigue'), cof: cLoc, det: 7, cause: `damage ${fmt(span.dYear, 3)} per year on the ${fmt(p.spanLength, 3)} m span`, effect: 'Girth-weld crack', action: 'Support the span, suppress slugging' },
    { name: 'Burst at shut-in pressure', pof: Math.max(pA('burst'), pA('defect')), cof: cLoc, det: 6, cause: `annual extreme ${fmt(rel.pMean * 10, 3)} bar differential`, effect: 'Rupture', action: 'Overpressure protection, defect repair' },
    { name: 'External collapse', pof: pCol, cof: cBase, det: 8, cause: `${fmt(sa.collapseUtil, 3)} utilisation when depressurised`, effect: 'Flattened pipe, flooding', action: 'Keep minimum internal pressure, buckle arrestors' },
    { name: 'Local or global buckling', pof: clamp(pLB + pUph, 0, 1), cof: cBase, det: 6, cause: `local buckling utilisation ${fmt(lbOp.util, 3)}`, effect: 'Wrinkle or upheaval', action: 'Span correction, rock cover' },
    { name: 'Hydrate or deposit plug', pof: clamp(p.plugProb, 0, 1), cof: cPlug, det: 3, cause: `${fmt(100 * beta, 3)} % of the bore lost`, effect: `Production stopped for ${fmt(p.plugDays, 3)} d`, action: 'Inhibitor, insulation, no one-sided depressurisation' },
    { name: 'Plug projectile / overpressure', pof: pPlugLoc, cof: cLoc, det: 7, cause: `${fmt(plugDp / BAR, 3)} bar across the plug`, effect: 'Bend rupture', action: 'Two-sided depressurisation procedure' },
    { name: 'External damage', pof: clamp(p.extFreq, 0, 1), cof: cLoc, det: 8, cause: 'Anchors, dropped objects, trawling', effect: 'Dent or rupture', action: 'Protection, exclusion zone' },
  ].map((t) => { const pc = category(t.pof, POF_EDGES), cc = category(t.cof, COF_EDGES), lv = riskLevel(pc, cc); return { ...t, pc, cc, level: lv, risk: t.pof * t.cof, S: Math.min(2 * cc, 10), O: Math.min(2 * pc, 10), rpn: Math.min(2 * cc, 10) * Math.min(2 * pc, 10) * t.det }; });
  const riskCost = sum(threats.map((t) => t.risk)), topLevel = Math.max(...threats.map((t) => t.level)), topThreat = threats.reduce((a, b) => (b.risk > a.risk ? b : a));
  // Markov degradation: as-new → degraded → critical → failed, repair of the critical state after a successful inspection
  const half = Math.max(p.CA / 2, 0.05 * p.wtMm), m3 = Math.max(0.8 * p.wtMm - 2 * half, 0.05 * p.wtMm), l12 = rate / half, l34 = rate / m3, mu = clamp(p.pod, 0, 1) / Math.max(p.inspInterval, 0.05), s0 = lossNow < half ? 0 : lossNow < 2 * half ? 1 : 2, p0 = [0, 0, 0, 0]; p0[s0] = 1;
  const horizon = Math.max(remYears, 1), mk = markov([[0, l12, 0, 0], [0, 0, l12, 0], [mu, 0, 0, l34], [0, 0, 0, 0]], p0, horizon, 60), mk0 = markov([[0, l12, 0, 0], [0, 0, l12, 0], [0, 0, 0, l34], [0, 0, 0, 0]], p0, horizon, 60);
  // Bayesian network: wetting and inhibition → corrosion → wall loss → leak, with an inspection finding as evidence
  const wetFrac = clamp(deg.corr.filter((c) => !c.oilWet).length / N.length, 0.02, 0.98), allowRate = Math.max(p.CA, 0.1) / p.designLife, zc = Math.sqrt(Math.log(1 + cC * cC)) || 0.3, pHigh = (r) => clamp(1 - Phi(Math.log(allowRate / Math.max(r, 1e-9)) / zc + zc / 2), 0.001, 0.999), bare = deg.corr[deg.iC].bare;
  const bnNodes = [{ name: 'Water wetting', states: ['yes', 'no'], cpt: [[wetFrac, 1 - wetFrac]] }, { name: 'Inhibitor working', states: ['yes', 'no'], cpt: [[clamp(p.eta, 0.01, 0.999), 1 - clamp(p.eta, 0.01, 0.999)]] },
    { name: 'Corrosion', states: ['high', 'low'], parents: ['Water wetting', 'Inhibitor working'], cpt: [[pHigh(bare * (1 - p.eta)), 1 - pHigh(bare * (1 - p.eta))], [pHigh(bare), 1 - pHigh(bare)], [0.01, 0.99], [0.02, 0.98]] },
    { name: 'Wall loss', states: ['severe', 'minor'], parents: ['Corrosion'], cpt: [[0.85, 0.15], [0.03, 0.97]] }, { name: 'Leak', states: ['yes', 'no'], parents: ['Wall loss'], cpt: [[0.25, 0.75], [0.001, 0.999]] },
    { name: 'Inspection', states: ['thinning found', 'no thinning'], parents: ['Wall loss'], cpt: [[clamp(p.pod, 0.01, 0.99), 1 - clamp(p.pod, 0.01, 0.99)], [0.05, 0.95]] }];
  const bn = [['No evidence (prior)', {}], ['Inspection finds thinning', { Inspection: 'thinning found' }], ['Inspection finds no thinning', { Inspection: 'no thinning' }]].map(([label, e]) => ({ label, r: bayesNet(bnNodes, e) })), bnNow = wall.measured ? bn[lossNow > 0.5 * Math.max(p.CA, 0.1) ? 1 : 2] : bn[0];
  // risk-based inspection interval
  const toTarget = rel.tTarget === null ? 1e9 : Math.max(rel.tTarget - p.age, 0), inspect = clamp(Math.min(0.5 * remLife, 0.5 * fatLife, toTarget > 0 ? toTarget : 0.5, topLevel >= 3 ? 1 : 10, def.worst ? Math.max(0.5 * def.worst.life, 0.5) : 10), 0.5, 10);

  // ---- warnings and recommendations -----------------------------------------------------------------------------------
  const W = (level, msg) => warnings.push({ level, msg }), V = (name, value, limit, where) => viol.push([name, cell(value), cell(limit), where]);
  const rH = sa.rows[sa.iH], rV = sa.rows[sa.iV], rC = sa.rows[sa.iC], codeName = { dnv: 'DNV-ST-F101', b314: 'ASME B31.4', b318: 'ASME B31.8' }[p.code];
  if (st.source !== 'flow suite') W('info', 'Pressures, temperatures and velocities come from the kernel steady-state estimate; run the flow suite for the detailed profile and the transient slug data.');
  if (sa.hoopUtil > 1) { W('bad', `Pressure containment utilisation ${fmt(sa.hoopUtil, 3)} exceeds 1.0 at ${fmt(rH.nd.x / 1000, 3)} km (${codeName}); the wall needed is ${fmt(sa.minWall, 3)} mm against ${fmt(p.wtMm, 3)} mm.`); V('Pressure containment', sa.hoopUtil, 1, `${fmt(rH.nd.x / 1000, 3)} km`); recs.push(`Increase the wall thickness to at least ${fmt(Math.ceil(sa.minWall * 10) / 10, 3)} mm or lower the design pressure to ${fmt(sa.mawpOf(sa.tDes), 3)} bara: the pressure-containment utilisation is ${fmt(sa.hoopUtil, 3)}.`); }
  if (sa.vmUtil > 1) { W('bad', `Combined (von Mises) stress utilisation ${fmt(sa.vmUtil, 3)} exceeds 1.0 at ${fmt(rV.nd.x / 1000, 3)} km.`); V('Combined stress', sa.vmUtil, 1, `${fmt(rV.nd.x / 1000, 3)} km`); recs.push(`Reduce the restrained thermal stress (expansion loop or lower design temperature than ${fmt(p.designTemp, 3)} °C): the equivalent stress reaches ${fmt(Math.max(rV.sDes.vonMises, rV.sOp.vonMises) / MPA, 3)} MPa.`); }
  if (sa.collapseUtil > 1) { W('bad', `External collapse utilisation ${fmt(sa.collapseUtil, 3)} exceeds 1.0 at ${fmt(rC.nd.depth, 4)} m water depth when the line is depressurised.`); V('External collapse', sa.collapseUtil, 1, `${fmt(rC.nd.depth, 4)} m depth`); recs.push(`Keep at least ${fmt(Math.max((rC.nd.pe - sa.col.pc / sa.gLB) / BAR, 1), 3)} bara inside the line at the deepest point or increase the wall: collapse capacity is ${fmt(sa.col.pc / BAR, 3)} bar.`); }
  else if (sa.propUtil > 1) { W('warn', `A local buckle would propagate: external pressure is ${fmt(sa.propUtil, 3)} times the factored propagation pressure (${fmt(sa.col.ppr / BAR, 3)} bar). Buckle arrestors are required.`); recs.push(`Fit buckle arrestors below about ${fmt(sa.col.ppr / sa.gLB / (RHO_SW * G), 3)} m water depth: the propagation pressure is ${fmt(sa.col.ppr / BAR, 3)} bar against ${fmt((rC.nd.pe - PATM) / BAR, 3)} bar outside.`); }
  if (p.designPressure > sa.mawp * 1.0001) { W('bad', `The design pressure ${fmt(p.designPressure, 4)} bara exceeds the allowable pressure of the present wall, ${fmt(sa.mawp, 4)} bara.`); V('Maximum allowable pressure', p.designPressure, sa.mawp, 'line'); }
  else if (p.designPressure > sa.mawpEol) { W('warn', `With the predicted wall at the end of the design life the allowable pressure drops to ${fmt(sa.mawpEol, 4)} bara, below the design pressure of ${fmt(p.designPressure, 4)} bara.`); recs.push(`Plan a re-rating to ${fmt(sa.mawpEol, 3)} bara or a repair before year ${fmt(p.age + remLife, 3)}: wall loss of ${fmt(rate, 3)} mm/y takes the allowable pressure below the design pressure.`); }
  if (p.blowdownMinT < p.mdmt) { W('bad', `Blowdown cools the metal to ${fmt(p.blowdownMinT, 3)} °C, below the minimum design metal temperature of ${fmt(p.mdmt, 3)} °C: brittle-fracture risk.`); V('Minimum metal temperature', p.blowdownMinT, p.mdmt, 'blowdown'); recs.push(`Slow the blowdown or specify impact-tested steel for ${fmt(p.blowdownMinT - 5, 3)} °C: the predicted minimum is ${fmt(p.mdmt - p.blowdownMinT, 3)} °C below the material limit.`); }
  if (remLife < remYears) { W(remLife < 0.5 * remYears ? 'bad' : 'warn', `The corrosion and erosion allowance (${fmt(p.CA, 3)} mm, ${fmt(margin, 3)} mm left) is consumed in ${fmt(remLife, 3)} y at ${fmt(rate, 3)} mm/y; the remaining design life is ${fmt(remYears, 3)} y.`); V('Remaining life on the allowance (y)', remLife, remYears, `${fmt(deg.gov.x / 1000, 3)} km`);
    const etaNeed = deg.gov.corr > 0 && remYears > 0 ? 1 - ((1 - p.eta) * Math.max(margin / remYears - deg.gov.ero, 0)) / deg.gov.corr : null;
    recs.push(`Remaining life ${fmt(remLife, 3)} y at ${fmt(rate, 3)} mm/y: inspect within ${fmt(inspect, 2)} y${etaNeed !== null && etaNeed < 0.999 && etaNeed > p.eta ? ` and raise the inhibitor efficiency × availability from ${fmt(100 * p.eta, 3)} % to ${fmt(100 * etaNeed, 3)} % to reach the design life` : etaNeed !== null && etaNeed >= 0.999 ? `; inhibition alone cannot reach the design life — a corrosion-resistant liner or ${fmt(rate * remYears, 2)} mm of allowance is needed` : ''}.`); }
  else recs.push(`Remaining life on the allowance is ${fmt(remLife, 3)} y at ${fmt(rate, 3)} mm/y against ${fmt(remYears, 3)} y still required; the risk-based inspection interval is ${fmt(inspect, 2)} y.`);
  if (deg.erosionalRatio > 1) { const nv = N[deg.iV]; W('warn', `Mixture velocity ${fmt(nv.vm, 3)} m/s is ${fmt(deg.erosionalRatio, 3)} times the API RP 14E erosional velocity (C = ${fmt(p.c14e, 3)}) at ${fmt(nv.x / 1000, 3)} km.`); V('Erosional velocity ratio', deg.erosionalRatio, 1, `${fmt(nv.x / 1000, 3)} km`); recs.push(`Cut the rate by ${fmt(100 * (1 - 1 / deg.erosionalRatio), 2)} % or enlarge the bore to ${fmt(p.idMm * Math.sqrt(deg.erosionalRatio), 3)} mm at ${fmt(nv.x / 1000, 3)} km to bring the velocity below the erosional limit of ${fmt(deg.ero[deg.iV].ve, 3)} m/s.`); }
  if (deg.eroMax > 0.1) { W(deg.eroMax > 0.5 ? 'bad' : 'warn', `Sand erosion reaches ${fmt(deg.eroMax, 3)} mm/y at ${deg.bWorst.name} (${deg.bWorst.ero.model}).`); recs.push(`Limit sand to ${fmt((p.sandKgD * 0.1) / deg.eroMax, 2)} kg/d or use a ${fmt(Math.max(deg.bWorst.rOverD * 2, 5), 2)}D bend at ${fmt(deg.bWorst.nd.x / 1000, 3)} km to keep erosion below 0.1 mm/y.`); }
  if (deg.corr.some((c) => c.eta < 0.999 * p.eta)) W('warn', `Wall shear above ${fmt(p.tauCrit, 3)} Pa strips the inhibitor film on part of the line: the efficiency is reduced there (erosion–corrosion interaction).`);
  if (deg.eroMax > 0.05 && deg.corrMax > 0.05) W('info', 'Erosion and corrosion act together: sand removes protective scale and inhibitor films, so the combined rate can exceed the sum used here.');
  if (sour.sour) { W('warn', `H2S partial pressure ${fmt(pH2S, 3)} kPa exceeds 0.3 kPa: sour service, ${sour.label}. Materials must comply with ISO 15156.`); recs.push(`Specify sour-service steel (hardness ≤ 250 HV, ISO 15156-2 ${sour.label.split(' (')[0]}) — H2S partial pressure is ${fmt(pH2S, 3)} kPa.`); if (deg.yCO2 > 0 && deg.yCO2 / Math.max(p.h2s / 100, 1e-12) < 20) W('warn', 'CO2/H2S ratio is below 20: iron-sulphide films control corrosion and the CO2 model is no longer representative.'); }
  if (tlcMax > 0.25) W('warn', `Water condenses at up to ${fmt(tlcMax, 3)} g/m²/s at the top of the line, above the 0.25 g/m²/s screening level for top-of-line corrosion in stratified wet-gas flow.`);
  if (ctx.outputs?.solids?.sandBed) W('warn', 'The solids suite predicts a stationary sand bed: under-deposit corrosion is not covered by the CO2 model.');
  if (ds.on) { if (sl.f90 > 0) W(span.fivRatio > 0.8 && span.fivRatio < 1.25 ? 'bad' : 'info', `Slugs of ${fmt(ds.len, 3)} m at ${fmt(ds.v, 3)} m/s every ${fmt(ds.period, 3)} s load a 90° bend with ${fmt(sl.f90 / 1000, 3)} kN (dynamic load factor ${fmt(p.dlf, 2)}); source: ${ds.source}.`); } else W('info', 'No slug loading is applied: the flow regime at the loaded points is not intermittent.');
  if (span.fivRatio > 0.8 && span.fivRatio < 1.25) { V('Excitation / natural frequency', span.fivRatio, '0.8–1.25 avoided', `span at ${fmt(span.nd.x / 1000, 3)} km`); }
  if (span.slug.on && span.f1 < 3 * span.slug.freq) { const Lneed = p.spanLength * Math.sqrt(span.f1 / (3 * span.slug.freq)); recs.push(`Add a support at ${fmt(p.spanLength / 2, 3)} m or shorten the span to ${fmt(Lneed, 3)} m to move the first natural frequency (${fmt(span.f1, 3)} Hz) above 3 × the slug frequency (${fmt(span.slug.freq, 3)} Hz).`); }
  if (fatLife < remYears) { W(fatLife < 0.5 * remYears ? 'bad' : 'warn', `Fatigue life of the span weld is ${fmt(fatLife, 3)} y with a design fatigue factor of ${fmt(p.dff, 2)} (damage ${fmt(span.dYear, 3)} per year), shorter than the ${fmt(remYears, 3)} y still required.`); V('Fatigue life (y)', fatLife, remYears, `span at ${fmt(span.nd.x / 1000, 3)} km`);
    const share = [['slug passage', span.dSlug], ['vortex-induced vibration', span.dViv], ['shutdown and blowdown cycles', span.dOps], ['inlet-pressure cycles', span.dPuls]].reduce((a, b) => (b[1] > a[1] ? b : a));
    const cut = fmt(100 * (1 - (fatLife / Math.max(remYears, 1)) ** (1 / 3)), 2), pct = fmt((100 * share[1]) / span.dYear, 3);
    recs.push(share[0] === 'slug passage' || share[0] === 'vortex-induced vibration' ? `Fatigue life ${fmt(fatLife, 3)} y is governed by ${share[0]} (${pct} % of the damage): cut the stress range by ${cut} % — a span of ${fmt(p.spanLength * (fatLife / Math.max(remYears, 1)) ** (1 / 6), 3)} m instead of ${fmt(p.spanLength, 3)} m, or a weld detail better than class ${p.snClass}.`
      : `Fatigue life ${fmt(fatLife, 3)} y is governed by ${share[0]} (${pct} % of the damage): cut their stress range by ${cut} % (smaller pressure and temperature swings, slower ramps) or their number by ${fmt(100 * (1 - fatLife / Math.max(remYears, 1)), 2)} %.`); }
  if (span.viv.state !== 'no VIV' && span.sub) W('warn', `Reduced velocity ${fmt(span.viv.vr, 3)} at the span (current ${fmt(p.currentSpeed + span.uw, 3)} m/s, f₁ = ${fmt(span.f1, 3)} Hz) is above the onset of ${span.viv.state} (in-line ${fmt(span.viv.onsetIL, 3)}, cross-flow ${fmt(span.viv.onsetCF, 3)}).`);
  if (span.eulerUtil > 1) { W('bad', `The effective axial compression ${fmt(-span.seff / 1000, 3)} kN exceeds the Euler buckling load of the span, ${fmt(span.euler / 1000, 3)} kN.`); V('Span Euler buckling', span.eulerUtil, 1, `span at ${fmt(span.nd.x / 1000, 3)} km`); }
  if (lbOp.util > 1) { W('bad', `Local buckling utilisation of the span (pressure + bending + axial force, load-controlled) is ${fmt(lbOp.util, 3)}.`); V('Local buckling (combined loading)', lbOp.util, 1, `span at ${fmt(span.nd.x / 1000, 3)} km`); }
  if (span.wSub < 0) W('warn', `The pipe is buoyant at the span (${fmt(span.wSub, 3)} N/m): weight coating or anchoring is needed for on-bottom stability.`);
  const bk = sa.buckling;
  if (bk.buried && bk.palmerUtil > 1) { W('bad', `Upheaval buckling: the download needed on a ${fmt(p.imperfection, 2)} m imperfection is ${fmt(bk.up.wReq / 1000, 3)} kN/m but only ${fmt(bk.resist / 1000, 3)} kN/m is available.`); V('Upheaval buckling (download)', bk.palmerUtil, 1, `${fmt(bk.x / 1000, 3)} km`); recs.push(`Raise the cover from ${fmt(p.coverDepth, 2)} m to ${fmt(bk.coverReq, 2)} m${bk.coverReq >= 10 ? ' or more' : ''} at ${fmt(bk.x / 1000, 3)} km: the effective compression is ${fmt(bk.comp / 1000, 3)} kN.`); }
  else if (!bk.buried && bk.hobbsUtil > 1) { W('warn', `Effective compression ${fmt(bk.comp / 1000, 3)} kN is ${fmt(bk.hobbsUtil, 3)} times the Hobbs buckling force (${bk.hb.governing}): the exposed line will buckle laterally and must be designed for controlled buckling.`); recs.push(`Design buckle initiators (sleepers or snake lay) about every ${fmt(clamp(2 * bk.hb.modes.find((m) => m.name === bk.hb.governing).length * Math.sqrt(bk.hobbsUtil), 500, 5000), 2)} m: the compression of ${fmt(bk.comp / 1000, 3)} kN exceeds the ${fmt(bk.hb.critical / 1000, 3)} kN lateral buckling force.`); }
  for (const d of def.rows) if (!d.ok) { W('bad', `${d.name} at ${fmt(d.x, 4)} m (${fmt(100 * d.ratio, 3)} % deep, ${fmt(d.length, 3)} mm long) fails Modified B31G: safe pressure ${fmt(d.safe / BAR, 3)} bar against ${fmt(d.dp / BAR, 3)} bar required.`); V(`${d.name}: estimated repair factor`, d.erf, 1, `${fmt(d.x, 4)} m`); recs.push(`Repair ${d.name.toLowerCase()} at ${fmt(d.x, 4)} m or de-rate to ${fmt((d.safe + d.nd.pe) / BAR, 3)} bara: its repair factor is ${fmt(d.erf, 3)}.`); }
  for (const d of def.rows) if (d.ok && d.life < remYears) recs.push(`${d.name} at ${fmt(d.x, 4)} m is acceptable today (repair factor ${fmt(d.erf, 3)}) but grows to the limit in ${fmt(d.life, 3)} y at ${fmt(d.rate, 3)} mm/y: re-inspect within ${fmt(Math.max(d.life / 2, 0.5), 2)} y.`);
  if (beta > 0.3 || plugged) W(beta > 0.6 ? 'bad' : 'warn', plugged ? `The line is plugged at ${fmt(ndB.x / 1000, 3)} km: upstream pressure rises to the shut-in pressure of ${fmt(p.pShutIn, 4)} bara.` : `Deposits take ${fmt(100 * beta, 3)} % of the flow area at ${fmt(ndB.x / 1000, 3)} km and add ${fmt(dpRestr / BAR, 3)} bar of pressure drop.`);
  if (p.pShutIn > mawpInlet) { W(plugged ? 'bad' : 'warn', `A full blockage would expose the line to the shut-in pressure of ${fmt(p.pShutIn, 4)} bara, above its allowable pressure of ${fmt(mawpInlet, 4)} bara at the inlet.`); V('Shut-in pressure against allowable pressure (bara)', p.pShutIn, mawpInlet, 'upstream of a plug'); recs.push(`Set the overpressure protection at ${fmt(mawpInlet, 3)} bara or below: the source can deliver ${fmt(p.pShutIn, 3)} bara against a plug.`); }
  if (plugDp > 5 * BAR && (beta > 0.3 || p.plugProb > 0.05)) recs.push(`Never depressurise a plug from one side: ${fmt(plugDp / BAR, 3)} bar across a ${fmt(p.plugLen, 3)} m plug could launch it at up to ${fmt(plugV, 3)} m/s (${fmt(plugE / 1e6, 3)} MJ) into the next bend.`);
  if (rel.pofNow > rel.target) { W(rel.pofNow > 10 * rel.target ? 'bad' : 'warn', `Annual probability of failure ${fmt(rel.pofNow, 2)} exceeds the target ${fmt(rel.target, 1)} of safety class ${p.sc}.`); V('Annual probability of failure', rel.pofNow, rel.target, rel.gov.name); }
  else if (rel.tTarget !== null) recs.push(`The annual failure probability crosses the ${fmt(rel.target, 1)} target in year ${fmt(rel.tTarget, 3)} (${rel.gov.name}): schedule the inspection or repair before then.`);
  if (topLevel >= 2) recs.push(`Highest risk: ${topThreat.name} at ${fmt(topThreat.risk / 1e6, 3)} M$/y (${fmt(topThreat.pof, 2)} per year × ${fmt(topThreat.cof / 1e6, 3)} M$): ${topThreat.action.toLowerCase()}.`);
  if (crackLife < remYears) { W('warn', `A ${fmt(p.flawMm, 2)} mm flaw grows to the critical depth of ${fmt(aC * 1000, 3)} mm in ${fmt(crackLife, 3)} y (Paris law).`); recs.push(`Tighten the weld acceptance flaw size below ${fmt(p.flawMm, 2)} mm or inspect the span welds every ${fmt(Math.max(crackLife / 3, 0.5), 2)} y: crack-growth life is ${fmt(crackLife, 3)} y.`); }
  if (span.note) W('info', span.note);
  if (wall.measured && upd) W('info', `Measured wall loss of ${fmt(lossMeas, 3)} mm in ${fmt(p.age, 3)} y updated the wall-loss rate from ${fmt(deg.gov.rate, 3)} to ${fmt(rate, 3)} mm/y.`);
  const hydro = p.hydrotestFactor * p.designPressure, hydroHoop = ((hydro * BAR - PATM) * (p.D - p.t)) / (2 * p.t) / MPA; if (hydroHoop > 0.96 * p.smys) W('warn', `The hydrotest at ${fmt(hydro, 4)} bara stresses the wall to ${fmt((100 * hydroHoop) / p.smys, 3)} % of yield.`);

  // ---- equipment loading, deposits with time, strains, particle impact and exceedance events ---------------------------
  const bendMax = sl.bends.reduce((a, b) => (b.force > a.force ? b : a)), supUtil = p.supportCapacity > 0 ? bendMax.force / (p.supportCapacity * 1000) : 0;
  if (supUtil > 1) { W('bad', `The slug force on ${bendMax.name.toLowerCase()} (${fmt(bendMax.force / 1000, 3)} kN) exceeds the support capacity of ${fmt(p.supportCapacity, 3)} kN.`); V('Bend support loading (kN)', bendMax.force / 1000, p.supportCapacity, `${fmt(bendMax.x / 1000, 3)} km`); recs.push(`Strengthen the support at ${fmt(bendMax.x / 1000, 3)} km to at least ${fmt(bendMax.force / 1000, 3)} kN or keep the slug velocity below ${fmt(bendMax.v / Math.sqrt(supUtil), 3)} m/s.`); }
  const idEff = p.ID * Math.sqrt(1 - Math.min(beta, 0.9999)), idLim = p.ID * Math.sqrt(1 - clamp(p.blockLimit, 0.01, 0.9999)), daysToLimit = beta >= p.blockLimit ? 0 : p.depositRate > 0 ? ((idEff - idLim) * 1000) / (2 * p.depositRate) : null;
  if (beta >= p.blockLimit) V('Flow area lost to deposits', beta, p.blockLimit, `${fmt(ndB.x / 1000, 3)} km`);
  if (daysToLimit !== null && daysToLimit > 0 && daysToLimit < 365) { W('warn', `Deposits growing at ${fmt(p.depositRate, 3)} mm/d reach the blockage limit of ${fmt(100 * p.blockLimit, 3)} % in ${fmt(daysToLimit, 3)} d.`); recs.push(`Pig or treat the line at least every ${fmt(Math.max(daysToLimit / 2, 1), 2)} d: deposits growing at ${fmt(p.depositRate, 3)} mm/d reach ${fmt(100 * p.blockLimit, 3)} % blockage in ${fmt(daysToLimit, 3)} d.`); }
  const strainOf = (s) => ({ hoop: (s.hoopInner - p.nu * (s.longitudinal + s.radialInner)) / p.E, long: (s.longitudinal - p.nu * (s.hoopInner + s.radialInner)) / p.E }), eDes = strainOf(rV.sDes), dTdes = p.designTemp - p.tInstall, lossEol = Math.min(lossNow + rate * remYears, p.wtMm), tauMaxNode = N.reduce((a, b) => (b.tauW > a.tauW ? b : a));
  const bw = deg.bWorst, eW = erosionRate({ model: 'dnv', mp: deg.mp, U: bw.nd.vm, D: p.ID, dp: deg.dp, rhoM: bw.nd.rhoNS, muM: bw.nd.muNS, rhoP: p.sandDensity, rOverD: bw.rOverD }), mPart = (p.sandDensity * Math.PI * deg.dp ** 3) / 6;
  const log = p.opLog, exceed = [
    ['Operating-log pressure above the allowable pressure', log.filter((r) => r.p > mawpInlet).length, `${log.length} log rows, limit ${fmt(mawpInlet, 4)} bara`], ['Operating-log pressure above the design pressure', log.filter((r) => r.p > p.designPressure).length, `limit ${fmt(p.designPressure, 4)} bara`],
    ['Operating-log temperature above the design temperature', log.filter((r) => isNum(r.T) && r.T > p.designTemp).length, `limit ${fmt(p.designTemp, 3)} °C`], ['Operating-log temperature below the material limit', log.filter((r) => isNum(r.T) && r.T < p.mdmt).length, `limit ${fmt(p.mdmt, 3)} °C`],
    ['Slug force above the support capacity (per year)', supUtil > 1 && ds.on ? cell(ds.freq * YEAR, 3) : 0, p.supportCapacity > 0 ? `${fmt(bendMax.force / 1000, 3)} kN against ${fmt(p.supportCapacity, 3)} kN` : 'support capacity not given'],
    ['Points of the line above the erosional velocity', deg.ero.filter((e) => e.ratio > 1).length, `${N.length} points, API RP 14E C = ${fmt(p.c14e, 3)}`], ['Points where wall shear strips the inhibitor film', N.filter((n) => n.tauW > p.tauCrit).length, `limit ${fmt(p.tauCrit, 3)} Pa, largest ${fmt(tauMaxNode.tauW, 3)} Pa`],
    ['Blowdowns below the material limit (per year)', p.blowdownMinT < p.mdmt ? p.eventsBlowdown : 0, `${fmt(p.blowdownMinT, 3)} °C against ${fmt(p.mdmt, 3)} °C`], ['Shut-ins above the allowable pressure (per year)', p.pShutIn > mawpInlet ? cell(p.eventsShutdown * clamp(p.pfdProtect, 0, 1) + clamp(p.plugProb, 0, 1), 3) : 0, 'shutdowns with failed protection plus plugs'],
    ['Stress cycles above the S–N knee (per year)', cell(sum(span.spectrum.filter((s) => snCycles(s.range, p.snClass, p.snEnv, { scf: p.scf, t: p.wtMm }) < (p.snEnv === 'cp' ? 1e6 : 1e7)).map((s) => s.perYear)), 3), `class ${p.snClass}`],
  ];
  for (const [nm, n, lim] of exceed.slice(0, 4)) if (n > 0) { W('bad', `${nm}: ${n} of ${log.length} log rows (${lim}).`); V(nm, n, 0, 'operating log'); }

  // ---- plots ----------------------------------------------------------------------------------------------------------
  const rows = sa.rows, Tg = rel.Tg, years = linspace(0, Math.max(remYears, 1) * 1.25, 40), wallAt = (y, q) => Math.max(p.wtMm - lossNow - rul.rate * (q === 0.5 ? 1 : randomVariable({ dist: 'lognormal', mean: 1, cov: rul.rateSd / rul.rate }).x(PhiInv(q))) * y, 0);
  plots.push({ type: 'line', title: 'Stress utilisation along the line', xlabel: 'Distance (km)', ylabel: 'Utilisation (–)', zeroY: true, series: [{ name: 'Pressure containment, design case', x: xkm, y: rows.map((r) => r.hoopDes) }, { name: 'Pressure containment, operating', x: xkm, y: rows.map((r) => r.hoopOp), dash: true }, { name: 'Equivalent stress, design case', x: xkm, y: rows.map((r) => r.vmDes) }, { name: 'Equivalent stress, operating', x: xkm, y: rows.map((r) => r.vmOp), dash: true }, { name: 'External collapse (depressurised)', x: xkm, y: rows.map((r) => r.collapse) }], hlines: [{ y: 1, label: 'limit', color: '#dc2626' }], note: `${rows[0].basis}. Equivalent stress against ${p.code === 'dnv' ? 'the derated yield strength (reference check)' : '0.90 SMYS'}.` });
  { const dmax = Math.max(1.3 * Math.max(...N.map((n) => n.depth)), 300), dd = linspace(0, dmax, 40), cE = sa.colOf(Math.max(p.t - p.CA / 1000, 0.05 * p.t));
    plots.push({ type: 'line', title: 'Collapse utilisation against water depth', xlabel: 'Water depth (m)', ylabel: 'Utilisation (–)', zeroY: true, series: [{ name: 'Collapse, present wall', x: dd, y: dd.map((d) => (RHO_SW * G * d * sa.gLB) / sa.col.pc) }, { name: 'Collapse, allowance consumed', x: dd, y: dd.map((d) => (RHO_SW * G * d * sa.gLB) / cE.pc), dash: true }, { name: 'Buckle propagation', x: dd, y: dd.map((d) => (RHO_SW * G * d * sa.gLB) / sa.col.ppr) }], hlines: [{ y: 1, label: 'limit', color: '#dc2626' }], vlines: p.waterDepth > 0 ? [{ x: Math.max(...N.map((n) => n.depth)), label: 'deepest point' }] : [], note: `Elastic ${fmt(sa.col.pel / BAR, 4)} bar, plastic ${fmt(sa.col.pp / BAR, 4)} bar, combined with ${fmt(100 * p.f0, 2)} % out-of-roundness ${fmt(sa.col.pc / BAR, 4)} bar; resistance factor ${fmt(sa.gLB, 4)}.` }); }
  plots.push({ type: 'line', title: 'Corrosion and erosion rate profiles', xlabel: 'Distance (km)', ylabel: 'Wall-loss rate (mm/y)', logy: true, series: [{ name: 'CO2 corrosion, inhibited', x: xkm, y: deg.corr.map((c) => Math.max(c.rate, 1e-6)) }, { name: 'CO2 corrosion, uninhibited', x: xkm, y: deg.corr.map((c) => Math.max(c.bare, 1e-6)), dash: true }, { name: 'Sand erosion, straight pipe', x: xkm, y: deg.ero.map((e) => Math.max(e.straight, 1e-6)) }, { name: 'Sand erosion at the bends', x: deg.bends.map((b) => b.nd.x / 1000), y: deg.bends.map((b) => Math.max(b.ero.gov, 1e-6)), mode: 'points' }], note: `de Waard–Milliams ${p.corrModel} with fugacity, scale and glycol factors; values below 1e-6 mm/y are drawn at 1e-6.` });
  plots.push({ type: 'line', title: 'Wall thickness with time', xlabel: 'Years from today', ylabel: 'Wall thickness (mm)', series: [{ name: 'Median', x: years, y: years.map((y) => wallAt(y, 0.5)) }, { name: '10 % pessimistic', x: years, y: years.map((y) => wallAt(y, 0.9)), dash: true }, { name: '10 % optimistic', x: years, y: years.map((y) => wallAt(y, 0.1)), dash: true }], hlines: [{ y: p.wtMm - p.CA, label: 'allowance consumed', color: '#f59e0b' }, { y: Math.max(...rows.map((r) => r.tReq)) * 1000, label: 'pressure minimum', color: '#dc2626' }], vlines: [{ x: remYears, label: 'design life' }], note: `Rate ${fmt(rul.rate, 3)} ± ${fmt(rul.rateSd, 2)} mm/y${rul.updated ? ' (updated with the measured wall)' : ''}; remaining life P10 / P50 / P90 = ${fmt(rul.p10, 3)} / ${fmt(rul.p50, 3)} / ${fmt(rul.p90, 3)} y.` });
  plots.push({ type: 'line', title: 'Slug force on a 90° bend', xlabel: 'Time (s)', ylabel: 'Resultant force (kN)', zeroY: true, series: [{ name: ds.on ? 'Momentum force of the slug train' : 'Steady momentum force', x: sl.history.t, y: sl.history.f }], hlines: ds.on ? [{ y: sl.f90 / 1000, label: `with dynamic load factor ${fmt(p.dlf, 2)}`, color: '#dc2626' }] : [], note: ds.on ? `F = ρ·A·v²·√2 with slug-body density ${fmt(ds.rhoS, 3)} kg/m³ and film-region density ${fmt(ds.rhoF, 3)} kg/m³ at ${fmt(ds.v, 3)} m/s; lengths and periods scattered log-normally.` : 'No intermittent flow: steady momentum force of the mixture.' });
  plots.push({ type: 'line', title: 'Mode shapes of the span', xlabel: 'Position along the span (m)', ylabel: 'Normalised deflection (–)', series: span.modes.shapes.slice(0, 4).map((s, k) => ({ name: `Mode ${k + 1}: ${fmt(span.modes.f[k], 4)} Hz`, x: span.model.x.slice(), y: s.slice() })), note: `${span.model.x.length - 1} Hermite beam elements, effective mass ${fmt(span.me, 4)} kg/m (content and added mass included)${p.timoshenko ? ', shear deformation and rotary inertia included' : ''}.` });
  { const tt = span.resp ? thin(span.resp.t, 600) : [0, 1], ss = span.resp ? thin(span.resp.sigma.map((s) => s / MPA), 600) : [0, 0];
    plots.push({ type: 'line', title: 'Dynamic bending stress at the critical point of the span', xlabel: 'Time (s)', ylabel: 'Stress (MPa)', series: [{ name: span.resp ? `x = ${fmt(span.resp.x, 3)} m` : 'No slug excitation', x: tt, y: ss }], note: `Modal superposition of ${span.modes.f.length} modes, Newmark-β (average acceleration), Δt = ${fmt(span.dt, 2)} s, damping ${fmt(100 * p.zeta, 2)} %. Static sag stress ${fmt(span.sigStatic, 3)} MPa is not included.` }); }
  { const b = span.cycles.length ? histogram(span.cycles.map((c) => c.range), 12) : { centers: [0], counts: [0] }, w = new Array(b.centers.length).fill(0);
    if (span.cycles.length) for (const c of span.cycles) w[Math.min(w.length - 1, Math.max(0, Math.floor((c.range - (b.centers[0] - b.width / 2)) / b.width)))] += c.count * span.perYear;
    plots.push({ type: 'bar', title: 'Rainflow histogram of the slug-induced stress', ylabel: 'Cycles per year', categories: b.centers.map((c) => `${fmt(c, 2)} MPa`), series: [{ name: 'Cycles per year', values: w }], note: `ASTM E1049 rainflow count of ${span.simSlugs} simulated slug passages, scaled to one year.` }); }
  { const pts = span.spectrum.filter((s) => s.range > 0), smin = Math.max(Math.min(1, ...pts.map((s) => s.range * p.scf)) * 0.5, 0.05), S = logspace(smin, 600, 60);
    plots.push({ type: 'line', title: `S–N curve ${p.snClass} with the stress spectrum`, xlabel: 'Cycles', ylabel: 'Stress range (MPa)', logx: true, logy: true, series: [{ name: `Class ${p.snClass} (${{ air: 'air', cp: 'seawater, cathodic protection', free: 'free corrosion' }[p.snEnv]})`, x: S.map((s) => Math.min(snCycles(s, p.snClass, p.snEnv, { t: p.wtMm }), 1e16)), y: S }, { name: 'Applied ranges × SCF, cycles in the design life', x: pts.map((s) => Math.max(s.perYear * p.designLife, 1e-3)), y: pts.map((s) => s.range * p.scf), mode: 'points' }], note: `Points to the right of the curve fail on their own; Miner sum per year ${fmt(span.dYear, 3)} (slug ${fmt(span.dSlug, 2)}, VIV ${fmt(span.dViv, 2)}, operating cycles ${fmt(span.dOps + span.dPuls, 2)}).` }); }
  { const cy = cycYr > 0 && Number.isFinite(crack.N) ? crack.cycles.map((c) => c / cycYr) : crack.a.map((_, i) => i);
    plots.push({ type: 'line', title: 'Crack growth from the initial flaw', xlabel: cycYr > 0 && Number.isFinite(crack.N) ? 'Years' : 'No growth: the flaw is below the threshold', ylabel: 'Crack depth (mm)', series: [{ name: 'Paris–Erdogan growth', x: cy, y: crack.a.map((a) => a * 1000) }], hlines: [{ y: aC * 1000, label: 'critical depth', color: '#dc2626' }], note: `ΔK₀ = ${fmt(crack.dK0, 3)} MPa√m against a threshold of ${fmt(p.dKth, 2)}; equivalent range ${fmt(dSeq, 3)} MPa, ${fmt(cycYr, 3)} cycles per year, edge-crack geometry factor.` }); }
  plots.push({ type: 'line', title: 'Reliability index with time', xlabel: 'Age of the line (y)', ylabel: 'Reliability index β (–)', series: [...rel.LS.map((l) => ({ name: l.name, x: Tg, y: l.beta })), { name: 'System, annual', x: Tg, y: rel.pofT.map((q) => clamp(-PhiInv(clamp(q, 1e-23, 1 - 1e-12)), -5, 10)), dash: true }], hlines: [{ y: -PhiInv(rel.target), label: `target, safety class ${p.sc}`, color: '#dc2626' }], note: 'First-order reliability index of every limit state (capped at 10); wall-loss and fatigue limit states are cumulative, burst and collapse refer to the annual extreme load.' });
  { const gs = Array.from(rel.mc.g).filter(Number.isFinite), h = gs.length ? histogram(gs, 40) : { centers: [0, 1], counts: [0, 0] };
    plots.push({ type: 'line', title: 'Monte Carlo histogram of the governing limit state', xlabel: 'Limit-state value g (failure when g < 0)', ylabel: 'Samples', series: [{ name: rel.gov.name, x: h.centers, y: h.counts, mode: 'step' }], vlines: [{ x: 0, label: 'failure' }], note: `${rel.mc.n} samples at the end of the design life: ${rel.mc.failures} failures, pf = ${fmt(rel.mc.pf, 3)}; FORM ${fmt(rel.f.pf, 3)}, importance sampling ${fmt(rel.is.pf, 3)}.` }); }
  plots.push({ type: 'field', title: 'Risk matrix', xlabel: 'Consequence category', ylabel: 'Probability category', zlabel: 'Risk level', zunit: '', x: [1, 2, 3, 4, 5], y: [1, 2, 3, 4, 5], z: [1, 2, 3, 4, 5].map((pc) => [1, 2, 3, 4, 5].map((cc) => riskLevel(pc, cc))), zmin: 0, zmax: 3, cmap: 'turbo', markers: threats.map((t, i) => ({ x: t.cc + 0.12 * ((i % 3) - 1), y: t.pc + 0.12 * (Math.floor(i / 3) - 1), label: `${i + 1}` })), note: `Numbers refer to the threat table. Probability categories at 3e-5, 3e-4, 3e-3, 3e-2 per year (average over the remaining design life); consequence categories at 0.1, 1, 10, 100 M$. Levels: 0 low (broadly acceptable), 1 medium and 2 high (ALARP region), 3 very high (intolerable).` });
  plots.push({ type: 'line', title: 'Markov degradation states', xlabel: 'Years from today', ylabel: 'Probability (–)', zeroY: true, ymax: 1, series: [...['As new', 'Degraded', 'Critical', 'Failed'].map((nm, k) => ({ name: nm, x: mk.t, y: mk.p.map((r) => r[k]) })), { name: 'Failed, no inspection', x: mk0.t, y: mk0.p.map((r) => r[3]), dash: true }], note: `Transition rates from the wall-loss rate: ${fmt(l12, 3)} /y between the first states, ${fmt(l34, 3)} /y to failure; repair of the critical state at ${fmt(mu, 3)} /y (inspection every ${fmt(p.inspInterval, 3)} y, detection probability ${fmt(p.pod, 2)}).` });
  plots.push({ type: 'line', title: 'Defect assessment chart', xlabel: 'Axial defect length (mm)', ylabel: 'Allowable depth / wall thickness (–)', logx: true, zeroY: true, ymax: 1, series: [{ name: 'Modified B31G', x: def.Ls, y: def.curves.modified }, { name: 'ASME B31G', x: def.Ls, y: def.curves.b31g, dash: true }, { name: 'DNV-RP-F101', x: def.Ls, y: def.curves.dnv, dash: true }, ...(def.rows.length ? [{ name: 'Measured defects', x: def.rows.map((d) => clamp(d.length, 10, 3000)), y: def.rows.map((d) => Math.min(d.ratio, 1)), mode: 'points' }] : [])], hlines: [{ y: 0.8, label: '80 % depth limit', color: '#dc2626' }], note: `Acceptance at the design differential pressure of ${fmt(def.dpG / BAR, 4)} bar; defects above a curve need repair or de-rating.` });
  plots.push({ type: 'bar', title: 'Importance of the random variables', ylabel: 'α² (share of the uncertainty)', categories: rel.f.names, series: [{ name: rel.gov.name, values: rel.f.alpha.map((a) => a * a) }] });

  // ---- tables ---------------------------------------------------------------------------------------------------------
  const km = (r) => cell(r.nd.x / 1000, 4), mp = (x) => cell(x / MPA, 4);
  tables.push({ title: 'Code checks', columns: ['Check', 'Basis', 'Value', 'Limit', 'Utilisation', 'Location'], rows: [
    ['Pressure containment, design case', rows[0].basis.split(':')[0], cell(Math.max(sa.pLocal(rH.nd, rows[0].basis.startsWith('DNV') ? 1.1 : 1) - rH.nd.pe, 0) / BAR) + ' bar', cell(sa.fmtOf(rH.nd, sa.tDes, p.designTemp).allowDp / BAR) + ' bar', cell(rH.hoopDes), km(rH) + ' km'],
    ['Pressure containment, operating', rows[0].basis.split(':')[0], cell(Math.max(...rows.map((r) => Math.max(r.nd.pi - r.nd.pe, 0))) / BAR) + ' bar', '—', cell(Math.max(...rows.map((r) => r.hoopOp))), '—'],
    ['Hoop stress, Barlow / Lamé (design case)', 'σ = Δp·D/2t; thick cylinder at the bore', `${mp(rH.sDes.hoopThin)} / ${mp(rH.sDes.hoopInner)} MPa`, mp(sa.fmtOf(rH.nd, wall.tNow, p.designTemp).hoopAllow) + ' MPa', '—', km(rH) + ' km'],
    ['Longitudinal stress', rV.restrained ? 'restrained: ν(σθ + σr) − E·α·ΔT' : 'unrestrained: end-cap', mp(Math.abs(rV.sDes.longitudinal) > Math.abs(rV.sOp.longitudinal) ? rV.sDes.longitudinal : rV.sOp.longitudinal) + ' MPa', mp(sa.fmtOf(rV.nd, wall.tNow, p.designTemp).longAllow) + ' MPa', cell(sa.longUtil), km(rV) + ' km'],
    ['Thermal stress, design case', '−E·α·ΔT', mp(rV.sDes.thermal) + ' MPa', '—', '—', km(rV) + ' km'],
    ['Von Mises / Tresca equivalent', p.code === 'dnv' ? 'reference check against the derated yield' : '0.90 SMYS', `${mp(Math.max(rV.sDes.vonMises, rV.sOp.vonMises))} / ${mp(Math.max(rV.sDes.tresca, rV.sOp.tresca))} MPa`, mp(sa.fmtOf(rV.nd, wall.tNow, p.designTemp).combAllow) + ' MPa', cell(sa.vmUtil), km(rV) + ' km'],
    ['External collapse', 'DNV-ST-F101 collapse with out-of-roundness', cell((rC.nd.pe - PATM) / BAR) + ' bar', cell(sa.col.pc / sa.gLB / BAR) + ' bar', cell(sa.collapseUtil), cell(rC.nd.depth) + ' m depth'],
    ['Buckle propagation', '35·f_y·(t/D)^2.5', cell((rC.nd.pe - PATM) / BAR) + ' bar', cell(sa.col.ppr / sa.gLB / BAR) + ' bar', cell(sa.propUtil), cell(rC.nd.depth) + ' m depth'],
    ['Local buckling, combined loading', `DNV-ST-F101 load-controlled, ${lbOp.mode}`, `M = ${cell((1.2 * mSpan) / 1000)} kN·m, S = ${cell((1.2 * span.seff) / 1000)} kN`, `M_p = ${cell(lbOp.Mp / 1000)} kN·m`, cell(lbOp.util), `span at ${cell(span.nd.x / 1000)} km`],
    ['Global buckling (Hobbs)', bk.hb.governing, cell(bk.comp / 1000) + ' kN', cell(bk.hb.critical / 1000) + ' kN', cell(bk.hobbsUtil), cell(bk.x / 1000) + ' km'],
    ['Upheaval hold-down (Palmer)', `imperfection ${cell(p.imperfection)} m`, cell(bk.up.wReq / 1000) + ' kN/m needed', cell(bk.resist / 1000) + ' kN/m available', cell(bk.palmerUtil), `cover needed ${cell(bk.coverReq, 3)} m`],
    ['Euler buckling of the span', `π²EI/(K·L)², K = ${END_K[p.spanEnds] || 1}`, cell(Math.max(-span.seff, 0) / 1000) + ' kN', cell(span.euler / 1000) + ' kN', cell(span.eulerUtil), `span at ${cell(span.nd.x / 1000)} km`],
    ['Minimum required wall', 'pressure containment + allowances', cell(sa.minWall) + ' mm', cell(p.wtMm) + ' mm nominal', cell(sa.minWall / p.wtMm), km(rH) + ' km'],
    ['Allowable pressure, present / end-of-life wall', 'referred to the design-pressure datum', `${cell(sa.mawp)} / ${cell(sa.mawpEol)} bara`, cell(p.designPressure) + ' bara design', cell(p.designPressure / Math.max(sa.mawp, 1e-9)), '—'],
    ['Hydrotest', `${cell(p.hydrotestFactor)} × design pressure`, cell(hydro) + ' bara', `${cell((100 * hydroHoop) / p.smys, 3)} % of yield`, cell(hydroHoop / (0.96 * p.smys)), '—'],
    ['Low-temperature check', 'blowdown minimum against the material limit', cell(p.blowdownMinT) + ' °C', cell(p.mdmt) + ' °C', p.blowdownMinT < p.mdmt ? 'fails' : 'passes', '—'],
  ], note: `Wall used: present ${fmt(wall.tNow * 1000, 4)} mm for stresses, ${fmt(sa.tDes * 1000, 4)} mm (nominal less allowances) for the design case. External pressure from the water depth of every point.` });
  const crit = [
    { x: rH.nd.x, mechanism: 'Pressure containment', utilisation: sa.hoopUtil }, { x: rV.nd.x, mechanism: 'Combined stress', utilisation: sa.vmUtil }, { x: rC.nd.x, mechanism: 'External collapse', utilisation: sa.collapseUtil },
    { x: ndC.x, mechanism: 'CO2 corrosion', utilisation: p.CA > 0 ? (deg.corrMax * p.designLife) / p.CA : deg.corrMax > 0 ? 9.99 : 0 }, { x: deg.bWorst.nd.x, mechanism: 'Sand erosion', utilisation: p.CA > 0 ? (deg.eroMax * p.designLife) / p.CA : deg.eroMax > 0 ? 9.99 : 0 },
    { x: N[deg.iV].x, mechanism: 'Erosional velocity', utilisation: deg.erosionalRatio }, { x: span.nd.x, mechanism: 'Span fatigue', utilisation: (span.dYear * p.dff * p.designLife + p.priorDamage) }, { x: span.nd.x, mechanism: 'Local buckling at the span', utilisation: lbOp.util },
    { x: bk.x, mechanism: bk.buried ? 'Upheaval buckling' : 'Lateral buckling', utilisation: bk.buried ? bk.palmerUtil : bk.hobbsUtil }, ...(def.worst ? [{ x: def.worst.x, mechanism: 'Metal-loss defect', utilisation: def.worst.erf }] : []), ...(beta > 0 ? [{ x: ndB.x, mechanism: 'Blockage', utilisation: beta }] : []),
  ].map((c) => ({ x: sig(c.x, 5) ?? 0, mechanism: c.mechanism, utilisation: sig(cap(c.utilisation, 99), 4) ?? 0 })).sort((a, b) => b.utilisation - a.utilisation);
  tables.push({ title: 'Critical locations', columns: ['Distance (m)', 'Mechanism', 'Utilisation'], rows: crit.map((c) => [c.x, c.mechanism, c.utilisation]), note: 'Utilisation above 1 means the limit is exceeded; corrosion and erosion are measured against the allowance over the design life.' });
  tables.push({ title: 'Design and operating limit violations', columns: ['Limit', 'Value', 'Allowed', 'Where'], rows: viol.length ? viol : [['None', '—', '—', '—']] });
  tables.push({ title: 'Slug and flow-induced loads at the bends', columns: ['Bend', 'Distance (m)', 'Angle (deg)', 'Density (kg/m³)', 'Velocity (m/s)', 'Momentum force (kN)', 'Impulse per slug (kN·s)', 'Support utilisation', 'Steady force (kN)', 'Pressure thrust (kN)', 'Centrifugal load (kN/m)', 'ρv² (Pa)', 'Kinetic-energy band', 'Likelihood of failure'], rows: sl.bends.map((b) => [b.name, cell(b.x, 5), cell(b.angle), cell(b.rho), cell(b.v), cell(b.force / 1000), cell(b.impulse / 1000), p.supportCapacity > 0 ? cell(b.force / (p.supportCapacity * 1000), 3) : '—', cell(b.steady / 1000), cell(b.pressure / 1000), cell(b.centrifugal / 1000), cell(b.fiv.rhoV2), b.fiv.band, `${cell(b.fiv.lof, 3)} (${b.fiv.likelihood})`]),
    note: `Momentum force F = DLF·ρ·A·v²·√(2(1 − cos θ)); pressure thrust p·A·√(2(1 − cos θ)) is carried by the pipe wall unless there is an expansion joint. Pressure pulsation from the slug body: ±${fmt(sl.pulsation, 3)} bar (acceleration ${fmt(sl.dpAcc / BAR, 3)} bar, friction ${fmt(sl.dpFric / BAR, 3)} bar)${sl.flowAmplitude !== null ? `; the flow suite reports ±${fmt(sl.flowAmplitude, 3)} bar at the inlet` : ''}; unbalanced force on a straight run ${fmt(sl.runForce / 1000, 3)} kN. The last two columns are an Energy-Institute-style screening, not a detailed assessment.` });
  tables.push({ title: 'Span dynamics', columns: ['Quantity', 'Value', 'Unit'], rows: [['First natural frequency', cell(span.f1), 'Hz'], ['Higher modes', span.modes.f.slice(1, 4).map((f) => fmt(f, 4)).join(' / ') || '—', 'Hz'], ['Effective mass', cell(span.me), 'kg/m'], ['Submerged weight', cell(span.wSub), 'N/m'], ['Static sag', cell(span.stat.maxDeflection * 1000), 'mm'], ['Static bending stress', cell(span.sigStatic), 'MPa'],
    ['Peak dynamic bending stress', cell(span.dynMax), 'MPa'], ['Peak dynamic deflection', cell(span.dispMax * 1000), 'mm'], ['Slug passage time', cell(span.passage), 's'], ['Governing excitation', `${span.excitation.name}, ${fmt(span.excitation.f, 3)} Hz`, '—'], ['Excitation / first natural frequency', cell(span.fivRatio), '–'],
    ['Reduced velocity (current + waves)', cell(span.viv.vr), '–'], ['Onset in-line / cross-flow', `${fmt(span.viv.onsetIL, 3)} / ${fmt(span.viv.onsetCF, 3)}`, '–'], ['Stability parameter K_sd', cell(span.viv.ksd), '–'], ['Largest VIV amplitude in-line / cross-flow', `${fmt(span.vivWorst.aIL, 3)} / ${fmt(span.vivWorst.aCF, 3)}`, 'diameters'], ['Wave-induced velocity at the seabed', cell(span.uw), 'm/s'],
    ['Effective axial force (fully restrained)', cell(span.seffRaw / 1000), 'kN'], ['Effective axial force used at the span', cell(span.seff / 1000), 'kN'], ['Slug at the span: frequency / length / velocity', span.slug.on ? `${fmt(span.slug.freq, 3)} Hz / ${fmt(span.slug.len, 3)} m / ${fmt(span.slug.v, 3)} m/s` : 'no slugging', '—'], ['Sum of support reactions', cell(-span.stat.sumReactions / 1000), 'kN'], ['Applied static load', cell(span.stat.totalLoad / 1000), 'kN']], note: 'Vortex-induced vibration follows DNV-RP-F105-type response models at screening level (single mode, Rayleigh-distributed current). On an exposed seabed the compression at the span is capped at the lateral buckling force, because the line releases any larger force by buckling.' });
  tables.push({ title: 'Stresses, strains and deformation at the most stressed point', columns: ['Quantity', 'Operating', 'Design case', 'Unit'], rows: [['Hoop stress at the bore (Lamé)', mp(rV.sOp.hoopInner), mp(rV.sDes.hoopInner), 'MPa'], ['Radial stress at the bore', mp(rV.sOp.radialInner), mp(rV.sDes.radialInner), 'MPa'], ['Longitudinal stress', mp(rV.sOp.longitudinal), mp(rV.sDes.longitudinal), 'MPa'], ['Von Mises equivalent', mp(rV.sOp.vonMises), mp(rV.sDes.vonMises), 'MPa'], ['Tresca equivalent', mp(rV.sOp.tresca), mp(rV.sDes.tresca), 'MPa'],
    ['Hoop strain (elastic)', cell(strainOf(rV.sOp).hoop * 1e6), cell(eDes.hoop * 1e6), 'µε'], ['Longitudinal strain (elastic)', cell(strainOf(rV.sOp).long * 1e6), cell(eDes.long * 1e6), 'µε'], ['Free thermal strain α·ΔT', cell(p.alphaT * (rV.nd.T - p.tInstall) * 1e6), cell(p.alphaT * dTdes * 1e6), 'µε'], ['Radial growth of the bore', cell(strainOf(rV.sOp).hoop * p.ID * 500), cell(eDes.hoop * p.ID * 500), 'mm'],
    ['Free end expansion per km if unrestrained', cell(p.alphaT * (rV.nd.T - p.tInstall) * 1e6), cell(p.alphaT * dTdes * 1e6), 'mm/km'], ['Span sag (static) / dynamic deflection', cell(span.stat.maxDeflection * 1000), cell(span.dispMax * 1000), 'mm']], note: `Point at ${fmt(rV.nd.x / 1000, 4)} km, ${rV.restrained ? 'axially restrained' : 'free to expand'}; linear-elastic (Hooke) strains from the three principal stresses.` });
  tables.push({ title: 'Wall shear and particle impact', columns: ['Quantity', 'Value', 'Unit'], rows: [['Largest wall shear stress', cell(tauMaxNode.tauW), 'Pa'], ['Its location', cell(tauMaxNode.x / 1000), 'km'], ['Wall shear / inhibitor-film limit', cell(tauMaxNode.tauW / Math.max(p.tauCrit, 1e-9)), '–'], ['Sand mass rate', cell(deg.mp * 1000), 'g/s'], ['Particle size / density', `${fmt(p.sandUm, 3)} / ${fmt(p.sandDensity, 4)}`, 'µm / kg/m³'], ['Particles per second', cell(deg.mp / Math.max(mPart, 1e-18)), '1/s'],
    ['Impact velocity at the worst fitting', cell(bw.nd.vm), 'm/s'], ['Characteristic impact angle', cell(eW.alpha), 'deg'], ['Kinetic energy of one particle', cell(0.5 * mPart * bw.nd.vm ** 2 * 1e6), 'µJ'], ['Impact power on the fitting', cell(0.5 * deg.mp * bw.nd.vm ** 2 * eW.G), 'W'], ['Share of particles that reach the wall (size correction)', cell(eW.G), '–'], ['Impact area', cell(eW.At * 1e4), 'cm²'], ['Angle function F(α) of the ductile wall', cell(eW.F), '–']], note: `Impact data at ${bw.name} (${fmt(bw.nd.x / 1000, 3)} km); wall shear from the flow picture.` });
  tables.push({ title: 'Exceedance events', columns: ['Event', 'Count', 'Basis'], rows: exceed.map((r) => [r[0], r[1], r[2]]), note: log.length ? 'Counts of the operating log refer to its rows.' : 'No operating log entered: add pressure and temperature history rows to count excursions and to add their cycles to the fatigue spectrum.' });
  tables.push({ title: 'Fatigue spectrum and damage', columns: ['Source', 'Stress range (MPa)', 'Cycles per year', 'Cycles to failure', 'Damage per year'], rows: span.spectrum.length ? span.spectrum.map((s) => { const Nf = snCycles(s.range, p.snClass, p.snEnv, { scf: p.scf, t: p.wtMm }); return [s.source, cell(s.range), cell(s.perYear), cell(Math.min(Nf, 1e30)), cell(s.perYear / Nf)]; }) : [['No cyclic loading', '—', '—', '—', '—']],
    note: `S–N class ${p.snClass}, ${{ air: 'air', cp: 'seawater with cathodic protection', free: 'free corrosion' }[p.snEnv]}, stress concentration ${fmt(p.scf, 3)}; Miner sum ${fmt(span.dYear, 3)} per year, life ${fmt(fatLife, 3)} y with a design fatigue factor of ${fmt(p.dff, 2)}. Crack-growth life from a ${fmt(p.flawMm, 2)} mm flaw: ${fmt(crackLife, 3)} y.` });
  tables.push({ title: 'Corrosion assessment at the governing point', columns: ['Quantity', 'Value', 'Unit'], rows: [['Location', cell(ndC.x / 1000), 'km'], ['Temperature', cell(ndC.T), '°C'], ['CO2 partial pressure / fugacity', `${fmt(deg.yCO2 * ndC.P, 3)} / ${fmt(cw.fCO2, 3)}`, 'bar'], ['pH of CO2-saturated water / used', `${fmt(cw.pHco2, 3)} / ${fmt(cw.pH, 3)}`, '–'], ['Scaling temperature', cell(cw.tScale), '°C'], ['Scale factor', cell(cw.Fscale), '–'], ['Glycol factor', cell(cw.Fglycol), '–'],
    ['Reaction-limited / mass-transfer-limited rate', `${fmt(cw.Vr, 3)} / ${fmt(cw.Vm, 3)}`, 'mm/y'], ['Uninhibited rate', cell(deg.corr[deg.iC].bare), 'mm/y'], ['Inhibitor efficiency applied', cell(100 * deg.corr[deg.iC].eta), '%'], ['Predicted rate', cell(deg.corrMax), 'mm/y'], ['Water wetting', deg.corr[deg.iC].oilWet ? 'oil-wet' : `water-wet (water cut ${fmt(100 * ndC.wcut, 3)} %, liquid ${fmt(ndC.vL, 3)} m/s)`, '—'],
    ['Mixed-potential model: corrosion potential', cell(mix.Ecorr), 'V'], ['Mixed-potential model: current density', cell(mix.icorr), 'A/m²'], ['Mixed-potential model: rate (Faraday)', cell(mix.rate), 'mm/y'], ['Limiting current H⁺ / H2CO3', `${fmt(mix.iLimH, 3)} / ${fmt(mix.iLimH2CO3, 3)}`, 'A/m²'], ['Mass-transfer-limited rate (Sherwood)', cell(mix.rateLimit), 'mm/y'],
    ['H2S partial pressure', cell(pH2S), 'kPa'], ['Sour-service domain', sour.label, '—'], ['Top-of-line condensation (max)', cell(tlcMax), 'g/m²/s'], ['Wall today / lost', `${fmt(wall.tNow * 1000, 4)} / ${fmt(lossNow, 3)}`, 'mm'], ['Cumulative loss at the end of the design life', cell(lossEol), 'mm'], ['Remaining wall at the end of the design life', cell(p.wtMm - lossEol), 'mm'], ['Life to the pressure-minimum wall', cell(lifeMin), 'y']], note: 'The mixed-potential model (Tafel kinetics, uninhibited bare steel) is an independent cross-check of the empirical rate; it does not include protective films.' });
  tables.push({ title: 'Sand erosion at fittings', columns: ['Fitting', 'Distance (m)', 'Velocity (m/s)', 'R/D', 'DNV-RP-O501 (mm/y)', 'Salama (mm/y)', 'Oka (mm/y)', 'Finnie (mm/y)', 'Governing (mm/y)', 'Corrosion there (mm/y)'], rows: deg.bends.map((b) => [b.name, cell(b.nd.x, 5), cell(b.nd.vm), cell(b.rOverD, 3), cell(b.ero.dnv), cell(b.ero.salama), cell(b.ero.oka), cell(b.ero.finnie), cell(b.ero.gov), cell(b.corr)]),
    note: `${fmt(p.sandKgD, 3)} kg/d of ${fmt(p.sandUm, 3)} µm sand. API RP 14E erosional velocity ratio ${fmt(deg.erosionalRatio, 3)} (C = ${fmt(p.c14e, 3)}). Oka and Finnie are single-impact models applied over the DNV impact area; straight-pipe erosion peaks at ${fmt(Math.max(...deg.ero.map((e) => e.straight)), 2)} mm/y.` });
  tables.push({ title: 'Defect assessment', columns: ['Defect', 'Location (m)', 'Depth (mm)', 'Length (mm)', 'Depth / wall', 'B31G burst (bar)', 'Modified B31G burst (bar)', 'DNV-RP-F101 burst (bar)', 'Safe pressure (bar)', 'Required (bar)', 'Repair factor', 'Verdict', 'Life to limit (y)', 'Growth rate (mm/y)', 'Depth at end of life (mm)'], rows: def.rows.length ? def.rows.map((d) => [d.name, cell(d.x, 5), cell(d.depth), cell(d.length), cell(d.ratio, 3), cell(d.pf.b31g / BAR), cell(d.pf.modified / BAR), cell(d.pf.dnv / BAR), cell(d.safe / BAR), cell(d.dp / BAR), cell(d.erf, 3), d.ok ? 'acceptable' : 'repair or de-rate', cell(d.life, 3), cell(d.rate, 3), cell(Math.min(d.depth + d.rate * remYears, p.wtMm), 3)]) : [['No defects entered', '—', '—', '—', '—', '—', '—', '—', '—', '—', '—', '—', '—', '—', '—']],
    note: `Differential (internal minus external) burst pressures of the nominal wall. Safe pressure = factor × Modified B31G burst pressure${mapS ? `; the wall map holds ${mapS.readings} readings, minimum wall ${fmt(mapS.minMm, 4)} mm` : ''}.` });
  tables.push({ title: 'Blockage and plug screening', columns: ['Quantity', 'Value', 'Unit'], rows: [['Flow area lost', cell(100 * beta), '%'], ['Effective bore', cell(idEff * 1000), 'mm'], ['Time to the blockage limit at the deposit growth rate', daysToLimit === null ? 'no growth entered' : cell(daysToLimit), 'd'], ['Effective roughness used for the slug friction', cell(p.roughUm), 'µm'], ['Location', cell(ndB.x / 1000), 'km'], ['Extra pressure drop of the restriction', cell(dpRestr / BAR), 'bar'], ['Upstream pressure', cell(pUp), 'bara'], ['Allowable pressure at the inlet', cell(mawpInlet), 'bara'], ['Shut-in pressure / allowable', cell(p.pShutIn / Math.max(mawpInlet, 1e-9)), '–'], ['Differential pressure across a plug', cell(plugDp / BAR), 'bar'], ['Force on the plug', cell((plugDp * p.Ai) / 1000), 'kN'], ['Plug mass', cell(plugMass), 'kg'], ['Run to the next bend', cell(runLen), 'm'], ['Projectile velocity (upper bound)', cell(plugV), 'm/s'], ['Projectile energy', cell(plugE / 1e6), 'MJ']], note: 'The projectile estimate assumes the full differential pressure acts over the run without gas expansion losses or wall friction: an upper bound for procedures, not a design load.' });
  tables.push({ title: 'Reliability methods on the governing limit state', columns: ['Method', 'Reliability index', 'Probability of failure', 'Evaluations', 'Note'], rows: [['FORM (HL–RF)', cell(rel.f.beta), cell(rel.f.pf, 3), rel.f.evals, rel.f.converged ? `converged in ${rel.f.iterations} iterations` : 'not converged'], ['SORM (Breitung)', cell(rel.so.beta), cell(rel.so.pf, 3), '—', rel.so.valid ? `curvatures ${rel.so.kappa.map((k) => fmt(k, 2)).join(', ') || '—'}` : 'curvature correction not applicable'],
    ['Monte Carlo', cell(rel.mc.beta), cell(rel.mc.pf, 3), rel.mc.n, rel.mc.failures ? `c.o.v. ${fmt(rel.mc.cov, 2)}` : 'no failure sampled'], ['Importance sampling', cell(rel.is.beta), cell(rel.is.pf, 3), rel.is.n, rel.is.cov !== null ? `c.o.v. ${fmt(rel.is.cov, 2)}` : '—'], ['Latin hypercube', cell(rel.lh.beta), cell(rel.lh.pf, 3), rel.lh.n, rel.lh.failures ? `${rel.lh.failures} failures` : 'no failure sampled'], ['Response surface + FORM', cell(rel.rs.beta), cell(rel.rs.pf, 3), rel.rs.evals, `sampling on the surface: ${fmt(rel.rs.pfSampling, 3)}`]],
    note: `${rel.gov.name}, evaluated at year ${fmt(rel.life, 3)}. Target annual failure probability for safety class ${p.sc}: ${fmt(rel.target, 1)} (β = ${fmt(-PhiInv(rel.target), 3)}).` });
  tables.push({ title: 'Random-variable sensitivities', columns: ['Variable', 'Distribution', 'Mean', 'Coefficient of variation', 'Design point', 'α', 'α²'], rows: rel.gov.vars.map((vv, i) => { const rv = randomVariable(vv); return [rv.name, rv.dist, cell(rv.mean), cell(rv.mean ? rv.sd / Math.abs(rv.mean) : 0, 3), cell(rel.f.x[i]), cell(rel.f.alpha[i], 3), cell(rel.f.alpha[i] ** 2, 3)]; }), note: 'α is the unit vector from the mean point to the design point in standard-normal space; a negative α is a resistance variable.' });
  tables.push({ title: 'Limit states: annual probability of failure', columns: ['Limit state', 'β today', 'Annual pf today', 'β at end of life', 'Annual pf at end of life', 'Annual pf, life average'], rows: [...rel.LS.map((l) => [l.name, cell(l.beta[0], 3), cell(l.annualNow, 3), cell(l.beta[Tg.length - 1], 3), cell(l.annualEol, 3), cell(l.annualAvg, 3)]), ['Series system', cell(-PhiInv(clamp(rel.pofNow, 1e-23, 1 - 1e-12)), 3), cell(rel.pofNow, 3), cell(-PhiInv(clamp(rel.pofEol, 1e-23, 1 - 1e-12)), 3), cell(rel.pofEol, 3), cell(rel.pofAvg, 3)]], note: 'The risk models and the published probability of failure use the life average; the reliability target applies to every single year.' });
  tables.push({ title: 'Threats, risk and FMECA ranking', columns: ['No.', 'Threat', 'Cause', 'Effect', 'Annual probability', 'Consequence (M$)', 'Risk (M$/y)', 'Risk level', 'ALARP band', 'Severity', 'Occurrence', 'Detection', 'RPN', 'Action'], rows: threats.map((t, i) => [i + 1, t.name, t.cause, t.effect, cell(t.pof, 3), cell(t.cof / 1e6), cell(t.risk / 1e6, 3), RISK_NAMES[t.level], ['broadly acceptable', 'ALARP: reduce if practicable', 'ALARP: reduce unless grossly disproportionate', 'intolerable'][t.level], t.S, t.O, t.det, t.rpn, t.action]).sort((a, b) => b[12] - a[12]), note: `Probabilities are averages over the remaining design life; the part that grows with wall loss is credited with the inspection barrier (missed with probability ${fmt(pMiss, 2)}). Risk priority number = severity × occurrence × detection (each 1–10); severity and occurrence follow the matrix categories.` });
  tables.push({ title: 'Fault tree: minimal cut sets for loss of containment', columns: ['Cut set', 'Order', 'Probability per year', 'Share of the top event'], rows: ft.cutSets.map((c) => [c.events.join(' AND '), c.events.length, cell(c.p, 3), cell(ft.rareEvent > 0 ? c.p / ft.rareEvent : 0, 3)]), note: `Top event ${fmt(ft.top, 3)} per year (exact, independent basic events); rare-event sum ${fmt(ft.rareEvent, 3)}, upper bound ${fmt(ft.upperBound, 3)}. Most important basic event: ${ft.importance.reduce((a, b) => (b.fussellVesely > a.fussellVesely ? b : a)).name}.` });
  tables.push({ title: 'Event tree: outcomes of a leak', columns: ['Path', 'Outcome', 'Conditional probability', 'Frequency (1/y)', 'Consequence (M$)', 'Risk (M$/y)'], rows: etRows.map((o) => [o.label, o.outcome, cell(o.p, 3), cell(o.freq, 3), cell(o.cost / 1e6), cell((o.freq * o.cost) / 1e6, 3)]), note: `Initiating frequency ${fmt(ft.top, 3)} per year from the fault tree; expected consequence ${fmt(cLoc / 1e6, 4)} M$ per event (repair ${fmt(p.repairCost, 3)} M$, ${fmt(p.downtimeDays, 3)} d of production at ${fmt(prodDay / 1e6, 3)} M$/d, environmental and safety costs by path).` });
  tables.push({ title: 'Bow-tie summary', columns: ['Threat', 'Preventive barriers', 'Top event', 'Mitigating barriers', 'Consequences'], rows: [
    ['Internal corrosion', `Inhibition (${fmt(100 * p.eta, 3)} %), corrosion allowance ${fmt(p.corrAllow, 2)} mm, inspection every ${fmt(p.inspInterval, 2)} y`, 'Loss of containment', `Leak detection (${fmt(p.pDetect, 2)}), isolation (${fmt(p.pIsolate, 2)})`, 'Release, repair, deferred production'], ['Sand erosion', 'Sand control, velocity limit, long-radius bends', 'Loss of containment', 'Wall-thickness monitoring, isolation', 'Release at fittings'],
    ['Fatigue (slugging, VIV)', 'Span correction, slug control, weld quality', 'Loss of containment', 'Leak detection, isolation', 'Crack, rupture'], ['Overpressure / plug', `Protection system (fails ${fmt(p.pfdProtect, 2)} on demand), hydrate management`, 'Loss of containment', 'Emergency shutdown, ignition control', 'Rupture, projectile'], ['Collapse / buckling', 'Wall thickness, minimum internal pressure, buckle arrestors', 'Loss of integrity', 'Arrestors limit the damaged length', 'Flattened section, replacement']] });
  tables.push({ title: 'Bayesian network: corrosion chain with inspection evidence', columns: ['Evidence', 'P(evidence)', 'P(inhibitor working)', 'P(high corrosion)', 'P(severe wall loss)', 'P(leak)'], rows: bn.map((b) => [b.label + (b === bnNow ? ' ← current state of knowledge' : ''), cell(b.r.pEvidence, 3), cell(b.r.posterior['Inhibitor working'].yes, 3), cell(b.r.posterior.Corrosion.high, 3), cell(b.r.posterior['Wall loss'].severe, 3), cell(b.r.posterior.Leak.yes, 3)]), note: 'Water wetting and inhibitor state → corrosion → wall loss → leak, with the inspection result as a child of wall loss. The corrosion probabilities come from the corrosion model and its uncertainty; the wall-loss and leak tables are engineering priors.' });

  // ---- KPIs, outputs --------------------------------------------------------------------------------------------------
  const kp = (label, value, unit, status, help) => ({ label, value: sig(value, 4) ?? 0, unit, status, help }), relBeta = clamp(-PhiInv(clamp(rel.pofNow, 1e-23, 1 - 1e-12)), -5, 10);
  const kpis = [
    kp('Pressure containment utilisation', sa.hoopUtil, '–', stat(sa.hoopUtil, 0.9), `${codeName}, worst of the design and operating cases`), kp('Equivalent stress utilisation', sa.vmUtil, '–', stat(sa.vmUtil, 0.9), 'Von Mises against the combined-stress allowable'), kp('Collapse utilisation', sa.collapseUtil, '–', stat(sa.collapseUtil, 0.9), 'External pressure on the depressurised pipe'),
    kp('Allowable pressure, present wall', sa.mawp, 'bara', p.designPressure > sa.mawp * 1.0001 ? 'bad' : p.designPressure > sa.mawpEol ? 'warn' : 'ok', 'Referred to the design-pressure datum'), kp('Minimum required wall', sa.minWall, 'mm', sa.minWall > p.wtMm ? 'bad' : sa.minWall > 0.95 * p.wtMm ? 'warn' : 'ok', 'Pressure containment plus allowances'),
    kp('Peak force on a 90° bend', sl.peak90 / 1000, 'kN', 'ok', ds.on ? 'Slug momentum force with the dynamic load factor' : 'Steady momentum force (no slugging)'), kp('First natural frequency of the span', span.f1, 'Hz', span.fivRatio > 0.8 && span.fivRatio < 1.25 ? 'bad' : 'ok', 'Beam finite-element model'),
    kp('Excitation / natural frequency', span.fivRatio, '–', span.fivRatio > 0.8 && span.fivRatio < 1.25 ? 'bad' : span.fivRatio > 0.5 && span.fivRatio < 2 ? 'warn' : 'ok', `Governing excitation: ${span.excitation.name}`), kp('Fatigue life', fatLife, 'y', fatLife < 0.5 * remYears ? 'bad' : fatLife < remYears ? 'warn' : 'ok', `With a design fatigue factor of ${fmt(p.dff, 2)}`),
    kp('CO2 corrosion rate', deg.corrMax, 'mm/y', deg.corrMax > 0.5 ? 'bad' : deg.corrMax > 0.1 ? 'warn' : 'ok', 'Largest inhibited rate along the line'), kp('Sand erosion rate', deg.eroMax, 'mm/y', deg.eroMax > 0.5 ? 'bad' : deg.eroMax > 0.1 ? 'warn' : 'ok', `Largest of bends and straight pipe (${deg.bWorst.ero.model})`),
    kp('Erosional velocity ratio', deg.erosionalRatio, '–', stat(deg.erosionalRatio, 0.8), 'API RP 14E'), kp('Remaining life on the allowance', remLife, 'y', remLife < 0.5 * remYears ? 'bad' : remLife < remYears ? 'warn' : 'ok', `At ${fmt(rate, 3)} mm/y combined wall loss`),
    kp('Annual probability of failure, today', rel.pofNow, '1/y', rel.pofNow > 10 * rel.target ? 'bad' : rel.pofNow > rel.target ? 'warn' : 'ok', `Structural limit states in series; target ${fmt(rel.target, 1)} for safety class ${p.sc}`), kp('Annual probability of failure, life average', rel.pofAvg, '1/y', rel.pofEol > rel.target || rel.tTarget !== null ? (rel.tTarget !== null && rel.tTarget - p.age < 0.25 * remYears ? 'bad' : 'warn') : 'ok', rel.tTarget !== null ? `The target is crossed in year ${fmt(rel.tTarget, 3)}; end-of-life value ${fmt(rel.pofEol, 2)}` : 'The target is met over the whole design life'), kp('Reliability index, today', relBeta, '–', rel.pofNow > rel.target ? 'warn' : 'ok', 'Annual, series system of the limit states'),
    kp('Risk cost', riskCost / 1e6, 'M$/y', topLevel >= 3 ? 'bad' : topLevel >= 2 ? 'warn' : 'ok', `Highest level: ${RISK_NAMES[topLevel]} (${topThreat.name})`), kp('Inspection interval', inspect, 'y', inspect < 1 ? 'warn' : 'ok', 'Risk-based: half of the shortest remaining life, capped by the reliability target'),
  ];
  const balances = [{ name: 'Span static equilibrium: applied load against support reactions (N)', in: span.stat.totalLoad, out: -span.stat.sumReactions }, { name: 'Markov state probabilities at the end of the horizon', in: 1, out: sum(mk.p[mk.p.length - 1]) }, { name: 'Event-tree outcome frequencies against the initiating frequency (1/y)', in: ft.top, out: et.total },
    { name: 'Lamé invariant σr + σθ at the bore against the outer wall (MPa)', in: (rH.sDes.hoopInner + rH.sDes.radialInner) / MPA, out: (rH.sDes.hoopOuter + rH.sDes.radialOuter) / MPA }, { name: 'Bayesian network: posterior probabilities of the leak node', in: 1, out: bnNow.r.posterior.Leak.yes + bnNow.r.posterior.Leak.no }];
  const outputs = {
    slugForce: sig(sl.peak90 / 1000, 5) ?? 0, hoopUtil: sig(sa.hoopUtil, 5) ?? 0, vmUtil: sig(sa.vmUtil, 5) ?? 0, collapseUtil: sig(sa.collapseUtil, 5) ?? 0, mawp: sig(sa.mawp, 5) ?? 0, minWallRequired: sig(sa.minWall, 5) ?? 0,
    corrosionRate: sig(deg.corrMax, 5) ?? 0, erosionRate: sig(deg.eroMax, 5) ?? 0, wallLossRate: sig(rate, 5) ?? 0, remainingLife: sig(remLife, 5) ?? 0, fatigueDamagePerYear: sig(span.dYear, 5) ?? 0, fatigueLife: sig(fatLife, 5) ?? 0,
    fivRatio: sig(span.fivRatio, 5) ?? 0, naturalFrequency: sig(span.f1, 9) ?? 0, pof: sig(rel.pofAvg, 5) ?? 0, reliabilityIndex: sig(clamp(-PhiInv(clamp(rel.pofAvg, 1e-23, 1 - 1e-12)), -5, 10), 5) ?? 0, pofToday: sig(rel.pofNow, 5) ?? 0, riskLevel: RISK_NAMES[topLevel], riskCostPerYear: sig(riskCost, 5) ?? 0, consequence: sig(cLoc, 5) ?? 0, inspectionInterval: sig(inspect, 4) ?? 0,
    criticalLocations: crit.slice(0, 10), violations: viol.length, upheavalUtil: sig(bk.buried ? bk.palmerUtil : bk.hobbsUtil, 5) ?? 0, erosionalRatio: sig(deg.erosionalRatio, 5) ?? 0,
    // additional values
    mawpEndOfLife: sig(sa.mawpEol, 5) ?? 0, wallNow: sig(wall.tNow * 1000, 5) ?? 0, pofEndOfLife: sig(rel.pofEol, 5) ?? 0, targetPof: rel.target, lofFrequency: sig(ft.top, 5) ?? 0, slugFrequency: sig(span.slug.on ? span.slug.freq : ds.on ? ds.freq : 0, 5) ?? 0, peakDynamicStress: sig(span.dynMax, 5) ?? 0, staticSpanStress: sig(span.sigStatic, 5) ?? 0,
    pressurePulsation: sig(sl.pulsation, 5) ?? 0, crackGrowthLife: sig(crackLife, 5) ?? 0, remainingLifeP10: sig(rul.p10, 5) ?? 0, remainingLifeP90: sig(rul.p90, 5) ?? 0, localBucklingUtil: sig(lbOp.util, 5) ?? 0, propagationUtil: sig(sa.propUtil, 5) ?? 0, blockage: sig(beta, 5) ?? 0, sourService: sour.sour, hydrotestPressure: sig(hydro, 5) ?? 0,
    slugImpulse: sig(bendMax.impulse / 1000, 5) ?? 0, supportUtil: sig(supUtil, 5) ?? 0, cumulativeLoss: sig(lossNow, 5) ?? 0, cumulativeLossEndOfLife: sig(lossEol, 5) ?? 0, maxWallShear: sig(tauMaxNode.tauW, 5) ?? 0, exceedanceEvents: exceed.map((r) => ({ event: r[0], count: typeof r[1] === 'number' ? r[1] : 0 })), daysToBlockageLimit: daysToLimit === null ? null : sig(daysToLimit, 5),
    mcPof: sig(rel.is.pf, 5) ?? 0, mcBeta: sig(rel.is.beta ?? rel.f.beta, 5) ?? 0, formBeta: sig(rel.f.beta, 6) ?? 0, plugProjectileEnergy: sig(plugE / 1e6, 5) ?? 0, flowSource: st.source,
  };
  prog(1, 'Done');
  const worst = crit[0];
  return {
    summary: `${codeName} pressure containment is used to ${fmt(100 * sa.hoopUtil, 3)} %, collapse to ${fmt(100 * sa.collapseUtil, 3)} %; wall loss of ${fmt(rate, 3)} mm/y leaves ${fmt(remLife, 3)} y on the allowance, the span weld has a fatigue life of ${fmt(fatLife, 3)} y, and the annual failure probability averages ${fmt(rel.pofAvg, 2)} over the remaining life (risk ${RISK_NAMES[topLevel]}, ${fmt(riskCost / 1e6, 3)} M$/y); the most utilised item is “${worst.mechanism}” at ${fmt(worst.x / 1000, 3)} km.`,
    kpis, warnings, recommendations: recs, plots, tables, balances, outputs,
  };
}
// Synthetic "measured" data: generated once from calModel with corrMult = 0.8, inhibEff = 93 %, erosionMult = 1.6, damping = 3.2 % and 4–5 % noise.
const CAL_SAMPLE = [
  {cT: 35, cP: 60, cPco2: 0.8, cVliq: 1.2, cInh: 1, corrRate: 0.1231, cVmix: 8, cRho: 180, cSand: 15, eroRate: 0.002657, cFreq: 1.2, cForce: 400, vibAmp: 1.958},
  {cT: 45, cP: 75, cPco2: 1.5, cVliq: 1.8, cInh: 1, corrRate: 0.2898, cVmix: 10, cRho: 150, cSand: 40, eroRate: 0.0163, cFreq: 1.5, cForce: 400, vibAmp: 3.114},
  {cT: 55, cP: 85, cPco2: 2, cVliq: 2.2, cInh: 0, corrRate: 6.199, cVmix: 12, cRho: 140, cSand: 25, eroRate: 0.01836, cFreq: 1.7, cForce: 300, vibAmp: 5.871},
  {cT: 60, cP: 90, cPco2: 1.2, cVliq: 1.5, cInh: 0, corrRate: 2.921, cVmix: 14, cRho: 120, cSand: 60, eroRate: 0.08105, cFreq: 1.8, cForce: 300, vibAmp: 12.16},
  {cT: 65, cP: 95, cPco2: 2.4, cVliq: 2.6, cInh: 1, corrRate: 0.612, cVmix: 16, cRho: 110, cSand: 35, eroRate: 0.06664, cFreq: 1.85, cForce: 250, vibAmp: 8.953},
  {cT: 70, cP: 100, cPco2: 3, cVliq: 3, cInh: 1, corrRate: 0.8091, cVmix: 18, cRho: 100, cSand: 80, eroRate: 0.2174, cFreq: 1.9, cForce: 250, vibAmp: 6.277},
  {cT: 50, cP: 70, cPco2: 1, cVliq: 0.9, cInh: 0, corrRate: 1.673, cVmix: 9, cRho: 200, cSand: 20, eroRate: 0.004579, cFreq: 2, cForce: 300, vibAmp: 3.96},
  {cT: 75, cP: 110, cPco2: 2.8, cVliq: 2.4, cInh: 1, corrRate: 0.6436, cVmix: 20, cRho: 90, cSand: 120, eroRate: 0.4742, cFreq: 2.2, cForce: 400, vibAmp: 2.426},
  {cT: 40, cP: 65, cPco2: 1.8, cVliq: 2, cInh: 0, corrRate: 5.197, cVmix: 11, cRho: 160, cSand: 50, eroRate: 0.02416, cFreq: 1.75, cForce: 300, vibAmp: 8.133},
  {cT: 80, cP: 120, cPco2: 3.5, cVliq: 3.4, cInh: 1, corrRate: 0.6836, cVmix: 22, cRho: 85, cSand: 150, eroRate: 0.7618, cFreq: 2.5, cForce: 500, vibAmp: 1.562},
];
const CAL_VALIDATION = [
  {cT: 38, cP: 62, cPco2: 1, cVliq: 1.4, cInh: 1, corrRate: 0.1586, cVmix: 13, cRho: 130, cSand: 45, eroRate: 0.03958, cFreq: 1.6, cForce: 350, vibAmp: 3.858},
  {cT: 58, cP: 88, cPco2: 2.2, cVliq: 2.8, cInh: 0, corrRate: 7.563, cVmix: 17, cRho: 105, cSand: 70, eroRate: 0.1478, cFreq: 1.8, cForce: 280, vibAmp: 10.93},
  {cT: 68, cP: 98, cPco2: 2.6, cVliq: 1.9, cInh: 1, corrRate: 0.514, cVmix: 15, cRho: 125, cSand: 30, eroRate: 0.04012, cFreq: 1.95, cForce: 320, vibAmp: 5.6},
  {cT: 48, cP: 78, cPco2: 1.4, cVliq: 2.3, cInh: 0, corrRate: 4.496, cVmix: 19, cRho: 95, cSand: 95, eroRate: 0.2994, cFreq: 2.1, cForce: 380, vibAmp: 3.093},
  {cT: 72, cP: 105, cPco2: 3.2, cVliq: 3.1, cInh: 1, corrRate: 0.9601, cVmix: 21, cRho: 88, cSand: 110, eroRate: 0.5132, cFreq: 1.4, cForce: 420, vibAmp: 2.64},
  {cT: 62, cP: 92, cPco2: 1.6, cVliq: 1.1, cInh: 1, corrRate: 0.2207, cVmix: 10, cRho: 170, cSand: 18, eroRate: 0.006022, cFreq: 2.35, cForce: 450, vibAmp: 1.803},
];

// ---- calibration model: closed-form predictions for one measured operating point ------------------------------------------
function calModel(v0) {
  const v = { ...DEFAULTS, ...v0 }, ID = pos(v.idMm, BASE.idMm) / 1000, t = pos(v.wtMm, BASE.wtMm) / 1000, D = ID + 2 * t, Dh = Math.max(nz(v.odCoatMm, 0) / 1000, D), Ai = (Math.PI * ID * ID) / 4;
  const corr = deWaardMilliams({ T: nz(v.cT, 60), pCO2: Math.max(nz(v.cPco2, 1), 0), P: nz(v.cP, 80), model: v.corrModel, U: Math.max(nz(v.cVliq, 2), 0.01), d: ID, pH: nz(v.phAct, 0), glycolWt: nz(v.glycolWt, 0), inhibEff: (clamp(nz(v.inhibEff, 0), 0, 100) / 100) * clamp(nz(v.cInh, 1), 0, 1), mult: nz(v.corrMult, 1) }).rate;
  const ero = erosionRate({ model: 'dnv', mp: Math.max(nz(v.cSand, 10), 0) / 86400, U: Math.max(nz(v.cVmix, 10), 0), D: ID, dp: pos(v.sandUm, 250) * 1e-6, rhoM: pos(v.cRho, 150), muM: 1e-4, rhoP: pos(v.sandDensity, 2650), rOverD: 1.5, gf: nz(v.geomFactor, 1), mult: nz(v.erosionMult, 1) }).rate;
  // resonance curve of the first span mode (pinned–pinned) under a harmonic mid-span force
  const L = pos(v.spanLength, 15), EI = pos(v.eMod, BASE.E) * MPA * (Math.PI / 64) * (D ** 4 - ID ** 4), me = pos(v.rhoSteel, 7850) * (Math.PI / 4) * (D * D - ID * ID) + nz(v.coatDensity, 900) * (Math.PI / 4) * (Dh * Dh - D * D) + 500 * Ai + nz(v.addedMass, 1) * RHO_SW * (Math.PI / 4) * Dh * Dh;
  const f1 = (Math.PI / (2 * L * L)) * Math.sqrt(EI / me), r = Math.max(nz(v.cFreq, 0.3), 0) / f1, z = clamp(nz(v.damping, 2), 0.01, 100) / 100;
  return { corrRate: corr, eroRate: ero, vibAmp: ((nz(v.cForce, 500) * L ** 3) / (48 * EI) / Math.sqrt((1 - r * r) ** 2 + (2 * z * r) ** 2)) * 1000 };
}

// ---- verification ---------------------------------------------------------------------------------------------------------
function verify() {
  const out = [], chk = (name, got, expected, tol, note) => out.push({ name, expected, got, tol, pass: Number.isFinite(got) && Math.abs(got - expected) <= tol, note });
  const EI = 2.5e7, m = 300, L = 30, w0 = Math.sqrt(EI / m) / (L * L);
  { const s = pipeStress({ D: 0.5, t: 0.01, pi: 10e6, pe: 0 }); chk('Thin-wall hoop stress p·D/2t', s.hoopThin / MPA, 250, 1e-9, 'D = 500 mm, t = 10 mm, 100 bar'); }
  { const a = lame(50e6, 0, 0.1, 0.15, 0.1), b = lame(50e6, 0, 0.1, 0.15, 0.15);
    chk('Lamé hoop stress at the bore', a.hoop / MPA, 130, 1e-9, 'p(ro² + ri²)/(ro² − ri²), ri = 100 mm, ro = 150 mm, 500 bar'); chk('Lamé hoop stress at the outer wall', b.hoop / MPA, 80, 1e-9, '2p·ri²/(ro² − ri²)'); chk('Lamé radial stress at the bore equals −p', a.radial / MPA, -50, 1e-9, 'boundary condition');
    chk('Lamé invariant σr + σθ through the wall', (b.hoop + b.radial) / MPA, (a.hoop + a.radial) / MPA, 1e-9, 'the sum is constant (= 80 MPa)'); }
  chk('Von Mises of (100, 50, 0) MPa with 30 MPa shear', vonMises(100, 50, 0, 30), Math.sqrt(10200), 1e-9, '√(½Σ(σi − σj)² + 3τ²)');
  chk('Tresca of (120, −40, 30) MPa', tresca(120, -40, 30), 160, 1e-12, 'largest minus smallest principal stress');
  chk('Restrained thermal stress −E·α·ΔT', pipeStress({ D: 0.3, t: 0.015, E: 207e9, nu: 0.3, alpha: 1.17e-5, dT: 50, restrained: true }).thermal / MPA, -121.095, 1e-6, '207 GPa, 1.17e-5 /K, 50 K');
  chk('Elastic collapse pressure 2E(t/D)³/(1 − ν²)', collapsePressure({ D: 0.4, t: 0.02, E: 207e9, nu: 0.3, fy: 448e6 }).pel / MPA, (2 * 207e3 * 1.25e-4) / 0.91, 1e-9, 't/D = 0.05');
  { const c = collapsePressure({ D: 0.4, t: 0.02, E: 207e9, nu: 0.3, fy: 448e6, ovality: 0.01 }); chk('Combined collapse satisfies the cubic interaction', ((c.pc - c.pel) * (c.pc ** 2 - c.pp ** 2)) / ((c.pc * c.pel * c.pp * 0.01 * 0.4) / 0.02), 1, 1e-8, '(pc − pel)(pc² − pp²) = pc·pel·pp·f0·D/t'); }
  chk('Euler buckling load π²EI/L²', eulerLoad(1e7, 20), 246740.11005, 1e-3, 'EI = 1e7 N·m², L = 20 m');
  // beam finite elements
  const fSS = (n) => beamModes(beamModel({ L, EI, m, n }), 1).f[0], exact = (Math.PI / 2) * w0;
  chk('Simply-supported beam: first natural frequency', fSS(32), exact, 1e-6 * exact, '(π/2L²)√(EI/m), 32 elements');
  chk('Fixed–fixed beam: frequency constant', beamModes(beamModel({ L, EI, m, n: 40, ends: ['fixed', 'fixed'] }), 1).omega[0] / w0, 22.3733, 1e-3, 'ω₁ = 22.37·√(EI/mL⁴)');
  { const e1 = fSS(2) - exact, e2 = fSS(4) - exact, e3 = fSS(8) - exact; chk('Beam elements: observed order of convergence of the frequency', Math.log2(e1 / e2), 4, 0.25, `errors ${e1.toExponential(2)}, ${e2.toExponential(2)}, ${e3.toExponential(2)} Hz for 2, 4, 8 elements`); }
  { const s = beamStatic(beamModel({ L, EI, m, n: 10, ends: ['fixed', 'free'] }), { points: [{ x: L, F: 1000 }] }); chk('Cantilever tip deflection P·L³/3EI', s.w[s.w.length - 1], (1000 * L ** 3) / (3 * EI), 1e-9, 'exact for Hermite elements'); chk('Cantilever root moment P·L (moment balance)', Math.abs(s.moment[0]), 1000 * L, 1e-5, 'from the element end forces'); }
  { const s = beamStatic(beamModel({ L, EI, m, n: 10, ends: ['fixed', 'free'], kGA: 1e8 }), { points: [{ x: L, F: 1000 }] }); chk('Timoshenko cantilever: P·L³/3EI + P·L/κGA', s.w[s.w.length - 1], (1000 * L ** 3) / (3 * EI) + (1000 * L) / 1e8, 1e-9, 'shear-deformable element'); }
  { const s = beamStatic(beamModel({ L, EI, m, n: 12, ends: ['fixed', 'pinned'] }), { q: 500 }); chk('Static equilibrium: reactions balance the load', -s.sumReactions, 500 * L, 1e-6, 'propped cantilever under uniform load'); chk('Propped cantilever: fixed-end moment q·L²/8', Math.abs(s.moment[0]), (500 * L * L) / 8, 1e-4, 'beam analytical solution'); }
  { const Pe = eulerLoad(EI, L); chk('Compressed beam: frequency falls as √(1 − P/P_Euler)', beamModes(beamModel({ L, EI, m, n: 24, N: -0.5 * Pe }), 1).f[0], exact * Math.SQRT1_2, 1e-5 * exact, 'geometric stiffness, P = 0.5 P_E'); }
  { const md = beamModel({ L, EI, m, n: 8 }), fr = md.free, K = sub(md.K, fr), M = sub(md.M, fr), u0 = fr.map((d) => (d % 2 === 0 ? 0.01 * Math.sin((Math.PI * md.x[d / 2]) / L) : ((0.01 * Math.PI) / L) * Math.cos((Math.PI * md.x[(d - 1) / 2]) / L)));
    const r = newmark({ M, K, f: () => new Float64Array(fr.length), u0, dt: 0.002, steps: 1500 }); chk('Newmark-β: energy conserved in undamped free vibration', Math.max(...r.energy) / Math.min(...r.energy) - 1, 0, 1e-9, 'average acceleration, 1500 steps'); }
  { const sd = (dt) => { const r = newmark({ M: [[1]], K: [[(2 * Math.PI) ** 2]], f: () => [0], u0: [1], dt, steps: Math.round(0.25 / dt), energy: false }); return Math.abs(r.u[0]); }; chk('Newmark-β: second-order convergence in the time step', Math.log2(sd(0.0125) / sd(0.00625)), 2, 0.1, 'error at a quarter period of a 1 Hz oscillator (exact value 0)'); }
  { // fluid-force transfer: a long slug creeping over the span gives the static uniform-load bending stress
    const md = beamModel({ L, EI, m, n: 16 }), mo = beamModes(md, 5), r = slugResponse({ model: md, modes: mo, zeta: 0.3, slugs: [{ t0: 0, len: 500, v: 0.5 }], dw: 100, dt: 0.02, tEnd: 200, ro: 0.15, E: 2e11 }), Ih = EI / 2e11;
    chk('Slug load transfer: quasi-static stress q·L²/8·r/I', Math.abs(r.sigma[r.sigma.length - 1]) / MPA, ((100 * L * L) / 8) * (0.15 / Ih) / MPA, 0.02 * ((100 * L * L) / 8) * (0.15 / Ih) / MPA, 'modal superposition, 5 modes'); }
  // loads
  chk('Bend force on a 90° bend ρ·A·v²·√2', bendForce({ rho: 800, A: 0.05, v: 5, angle: 90 }).force, 1414.2135624, 1e-6, '800 kg/m³, 0.05 m², 5 m/s');
  { const b = bendForce({ rho: 800, A: 0.05, v: 5, angle: 60, radius: 2 }); chk('Momentum-transfer consistency: centrifugal load × chord = bend force', b.centrifugal * 2 * 2 * Math.sin(Math.PI / 6), 1000, 1e-9, '60° bend: ρAv²·2 sin 30° = 1000 N'); }
  // fatigue and fracture
  { const c = rainflow([-2, 1, -3, 5, -1, 3, -4, 4, -2]), cnt = (r) => sum(c.filter((x) => x.range === r).map((x) => x.count)); chk('Rainflow, ASTM E1049 example: cycles of range 4', cnt(4), 1.5, 1e-12, 'sequence −2, 1, −3, 5, −1, 3, −4, 4, −2'); chk('Rainflow, ASTM E1049 example: cycles of ranges 3, 6, 8, 9', cnt(3) + 10 * cnt(6) + 100 * cnt(8) + 1000 * cnt(9), 0.5 + 5 + 100 + 500, 1e-12, '0.5, 0.5, 1.0 and 0.5 cycles'); }
  chk('Miner sum of a constant-amplitude block', snDamage([{ range: 100, count: 1e5 }], 'D', 'air').damage, 1e5 / 10 ** 6.164, 1e-12, '1e5 cycles at 100 MPa, class D in air: N = 10^(12.164 − 3·2)');
  chk('S–N thickness correction (t/25)^k', snCycles(100, 'D', 'air', { t: 50 }) / snCycles(100, 'D', 'air'), 2 ** -0.6, 1e-12, 'class D, k = 0.2, 50 mm wall');
  { const a0 = 1e-3, ac = 1e-2, C = 1.65e-11, dS = 80, Y = 1.12, hand = (ac ** -0.5 - a0 ** -0.5) / (C * (Y * dS * Math.sqrt(Math.PI)) ** 3 * -0.5); chk('Paris-law life against the closed-form integral', parisLife({ a0, ac, C, m: 3, dS, Y }).N, hand, 1e-6 * hand, '1 mm to 10 mm, m = 3, constant geometry factor'); }
  // corrosion and erosion
  chk('Faraday: 1 A/m² on iron in mm/y', faradayRate(1), 1.16, 0.005, 'M = 55.85 g/mol, n = 2, ρ = 7870 kg/m³');
  chk('Butler–Volmer: linear response at small overpotential', butlerVolmer(1e-5, 1) / 1e-5, FARADAY / (RGAS * 298.15), 1e-3, 'i ≈ i0·F·η/RT for αa + αc = 1');
  chk('de Waard–Milliams 1991 at 1 bar CO2, 60 °C', deWaardMilliams({ T: 60, pCO2: 1, model: '1991', fugacity: false, scale: false }).rate, 4.647, 0.002, 'log V = 5.8 − 1710/T + 0.67 log pCO2');
  chk('Corrosion growth: remaining life = margin / rate', remainingLife({ margin: 3, rate: 0.3, cov: 0 }).mean, 10, 1e-9, '3 mm at 0.3 mm/y');
  chk('API RP 14E erosional velocity', erosionalVelocity(100, 100), (100 / Math.sqrt(100 / 16.018463)) * 0.3048, 0.01, 'C = 100, 100 kg/m³ = 6.243 lb/ft³ → 40.0 ft/s');
  chk('DNV-RP-O501 straight-pipe erosion', erosionRate({ model: 'dnv', geometry: 'straight', mp: 0.01, U: 10, D: 0.1 }).rate, 2.5e-5 * 398.1071706 * 100 * 0.01, 1e-9, '2.5e-5·U^2.6·D⁻²·ṁ');
  { const d = { D: 0.762, t: 0.009525, d: 0.00381, L: 0.2, smys: 359e6, smts: 455e6 }; chk('ASME B31G burst pressure of a textbook defect', b31g({ ...d, method: 'b31g' }).pf / MPA, 8.177, 0.003, '30 in × 9.53 mm X52, 40 % deep, 200 mm long: M = 2.326, S = 327.1 MPa'); chk('Modified B31G burst pressure of the same defect', b31g({ ...d, method: 'modified' }).pf / MPA, 8.436, 0.003, 'M = 2.087, flow stress 428 MPa'); }
  // reliability and risk mathematics
  chk('Standard normal distribution Φ(1.96)', Phi(1.96), 0.9750021, 1e-7, 'tabulated value');
  { const vars = [{ name: 'R', mean: 10, sd: 1.5 }, { name: 'S', mean: 5, sd: 1 }], g = (x) => x[0] - x[1], f = form(g, vars), bx = 5 / Math.sqrt(3.25), pfx = Phi(-bx);
    chk('FORM, linear limit state with normal variables', f.beta, bx, 1e-8, 'β = (μR − μS)/√(σR² + σS²)'); chk('FORM sensitivity factor of the resistance', f.alpha[0], -1.5 / Math.sqrt(3.25), 1e-6, 'α_R = −σR/√(σR² + σS²)');
    const mc = monteCarlo(g, vars, { n: 100000, seed: 5 }), se = Math.sqrt((pfx * (1 - pfx)) / 100000); chk('Monte Carlo against the exact failure probability', mc.pf, pfx, 3 * se, '100,000 samples, tolerance = 3 standard errors');
    const is = monteCarlo(g, vars, { n: 10000, seed: 6, method: 'is', center: f.u }); chk('Importance sampling at the design point', is.pf, pfx, 0.05 * pfx, '10,000 samples');
    const se1 = Math.abs(monteCarlo(g, vars, { n: 4000, seed: 9, method: 'is', center: f.u }).cov), se2 = Math.abs(monteCarlo(g, vars, { n: 64000, seed: 9, method: 'is', center: f.u }).cov); chk('Monte Carlo convergence: error falls as 1/√N', Math.log(se1 / se2) / Math.log(16), 0.5, 0.05, 'coefficient of variation for 4,000 and 64,000 samples'); }
  { const b = 3, k = 0.2, s = sorm((x) => b + 0.5 * k * x[0] * x[0] - x[1], [{ name: 'a', mean: 0, sd: 1 }, { name: 'b', mean: 0, sd: 1 }]); let ex = 0; const n = 4000; for (let i = 0; i < n; i++) { const u = -8 + (16 * (i + 0.5)) / n; ex += phi(u) * Phi(-(b + 0.5 * k * u * u)) * (16 / n); }
    chk('SORM (Breitung) on a parabolic limit state', s.pf, ex, 0.04 * ex, 'against numerical integration; FORM alone is 29 % high'); chk('SORM recovers the curvature of the limit state', s.kappa[0], k, 1e-4, 'principal curvature 0.2'); }
  { const U = lhs(50, 3, 3); let okS = 0; for (let j = 0; j < 3; j++) okS += new Set(U.map((r) => Math.floor(r[j] * 50))).size; chk('Latin hypercube: one sample in every stratum', okS, 150, 0, '50 samples, 3 dimensions'); }
  { const rv = randomVariable({ dist: 'weibull', mean: 2, cov: 0.25 }), q = [...Array(2000)].map((_, i) => rv.x(PhiInv((i + 0.5) / 2000))); chk('Weibull variable reproduces its mean through the normal transform', mean(q), 2, 0.004, 'quantile average over 2,000 strata'); }
  chk('Fault tree: A OR (B AND C)', faultTree({ name: 'T', gate: 'OR', children: [{ name: 'A', p: 0.01 }, { name: 'G', gate: 'AND', children: [{ name: 'B', p: 0.1 }, { name: 'C', p: 0.2 }] }] }).top, 1 - 0.99 * 0.98, 1e-12, '1 − (1 − 0.01)(1 − 0.02)');
  { const t = faultTree({ name: 'T', gate: 'AND', children: [{ name: 'G1', gate: 'OR', children: [{ name: 'A', p: 0.01 }, { name: 'B', p: 0.1 }] }, { name: 'G2', gate: 'OR', children: [{ name: 'A', p: 0.01 }, { name: 'C', p: 0.2 }] }] }); chk('Fault tree with a repeated event: (A OR B) AND (A OR C)', t.top, 1 - 0.99 * 0.98, 1e-12, 'Boolean reduction to A OR (B AND C)'); chk('Minimal cut sets after absorption', t.cutSets.length, 2, 0, '{A} and {B, C}'); }
  { const e = eventTree(1e-3, [{ name: 'a', p: 0.9 }, { name: 'b', p: 0.95 }, { name: 'c', p: 0.02 }]); chk('Event tree: outcome frequencies sum to the initiating frequency', e.total, 1e-3, 1e-15, 'probability normalisation'); }
  { const lam = 0.3, mu = 1.2, r = markov([[0, lam], [mu, 0]], [1, 0], 5, 50), last = r.p[r.p.length - 1]; chk('Markov chain: two-state analytical solution', last[0], mu / (lam + mu) + (lam / (lam + mu)) * Math.exp(-(lam + mu) * 5), 1e-10, 'availability at t = 5'); chk('Markov chain: probabilities sum to one', last[0] + last[1], 1, 1e-12, 'probability normalisation'); }
  chk('Bayesian network against Bayes\' theorem', bayesNet([{ name: 'D', states: ['y', 'n'], cpt: [[0.01, 0.99]] }, { name: 'T', states: ['pos', 'neg'], parents: ['D'], cpt: [[0.95, 0.05], [0.1, 0.9]] }], { T: 'pos' }).posterior.D.y, 0.0095 / (0.0095 + 0.099), 1e-12, 'P(D | positive test)');
  return out;
}

// ---- suite object ---------------------------------------------------------------------------------------------------------
const mm = (x) => (isNum(x) && x > 0 ? x * 1000 : undefined), okNum = (x) => (isNum(x) ? x : undefined);
export default {
  id: 'integ', num: 6, title: 'Integrity, Loads, Risk & Engineering Assessment', short: 'Integrity · Risk', icon: '🛡️',
  tagline: 'Turns pressure, temperature, slugs and solids into stress, fatigue, corrosion, erosion, reliability and risk.',
  description: 'Pipe stress and code utilisation along the line, collapse and buckling, slug forces on bends and the dynamic response of a span (beam finite elements, Newmark integration), rainflow fatigue and crack growth, CO2 corrosion and sand erosion with remaining life, defect assessment, structural reliability by FORM, SORM and sampling, and risk through fault tree, event tree, Markov and Bayesian models. Every number comes from the flow picture of the case, so the suite works before and after the upstream suites have been run.',
  guide: ['Check the pipe, material and design basis; pull the linked values from the network, flow, solids and operations suites.', 'Describe the span or jumper to be checked (length, supports, location) and the bends that take slug loads.', 'Set the corrosion and sand inputs; attach inspection data (minimum wall, wall map or a defect list) when available.', 'Run, then read the code checks, the critical locations and the recommendations; use the mesh tab to confirm the beam, time-step and sampling resolution.', 'Calibrate the corrosion multiplier, inhibitor efficiency, erosion multiplier and damping against coupon, probe and vibration data.'],
  implemented: ['conservation of linear momentum', 'conservation of angular momentum', 'stress–strain relations', "hooke's law", 'beam equations', 'euler–bernoulli beam equation', 'timoshenko beam equation', 'finite-element equilibrium equations', 'thin-wall hoop stress', 'lamé thick-cylinder equations', 'longitudinal stress', 'von mises equivalent stress', 'tresca criterion', 'thermal stress', 'combined-loading interaction equations',
    'momentum-flux equation', 'control-volume momentum balance', 'bend-force equation', 'transient pressure-force equation', 'centrifugal force from multiphase bends', 'one-way fsi', 'modal equations', 'vibration equations', 'vortex-induced-vibration models', 's–n curves', 'palmgren–miner cumulative-damage rule', 'rainflow cycle counting', 'fracture-mechanics crack-growth models', 'paris–erdogan law',
    'empirical/semi-mechanistic co₂ corrosion models', 'electrochemical corrosion kinetics', 'butler–volmer equation', "faraday's law", 'mass-transfer-limited corrosion', 'finnie erosion model', 'oka erosion model', 'dnv-type erosional assessment', 'particle-impact erosion correlations', 'erosional-velocity screening', 'euler buckling', 'local buckling equations', 'collapse-pressure models', 'upheaval/lateral buckling formulations',
    'limit-state function', 'first-order reliability method', 'second-order reliability method', 'monte carlo simulation', 'importance sampling', 'latin-hypercube sampling', 'response-surface reliability', 'fault-tree analysis', 'event-tree analysis', 'bow-tie modelling', 'bayesian networks', 'markov models', 'fmea/fmeca', 'consequence-frequency models', 'alarp-type risk assessment',
    'transient slug solver + structural dynamics', 'corrosion + remaining-life model', 'probabilistic reliability + deterministic integrity model', 'digital twin + remaining-useful-life prediction',
    // initial and boundary conditions
    'initial wall thickness', 'corrosion allowance', 'erosion allowance', 'roughness', 'residual stresses where considered', 'pre-existing defects', 'initial corrosion or erosion damage', 'blockage level', 'support condition', 'equipment condition and accumulated fatigue or damage state',
    'pipe supports', 'restraints', 'free spans', 'bends', 'internal pressure and temperature histories', 'external hydrostatic pressure', 'seabed/support interaction where relevant', 'slug impact and momentum loads', 'pressure pulsations', 'thermal loads', 'wall shear', 'solids impact conditions', 'corrosion environment and erosion exposure', 'maximum allowable operating pressure and temperature', 'allowable stress', 'minimum allowable wall thickness',
    'erosion and corrosion limits', 'fatigue limits', 'equipment design capacities', 'blockage limits and applicable safety factors', 'failure criteria', 'consequence categories', 'probability thresholds', 'inspection thresholds and acceptable risk levels',
    // inputs and outputs
    'pressure/temperature histories and transients', 'slug frequency, length, velocity, density, holdup and momentum', 'phase velocities/densities', 'solids/particle loading and impact data', 'evolving deposits/roughness', 'material properties', 'wall thickness and corrosion allowance', 'defects', 'bends, supports, anchors and free spans', 'external hydrostatic/seabed loads', 'corrosion/erosion models and environments', 'fatigue/fracture parameters', 'allowable design limits', 'probability/consequence assumptions and inspection data',
    'slug-induced forces and impulses', 'bend/momentum loads', 'wall shear and particle impact', 'stresses/strains/deformation where structural analysis is enabled', 'vibration/dynamic response', 'erosion/corrosion rates and cumulative loss', 'remaining wall thickness', 'fatigue damage/life and defect growth where modelled', 'buckling/collapse margins where relevant', 'critical locations', 'utilization/integrity margins', 'blockage severity', 'exceedance events', 'reliability/failure probability', 'consequence/risk indicators', 'equipment loading and design/operational-limit violations',
    // calibration and verification
    'damping', 'corrosion coefficients', 'erosion coefficients',
    'static equilibrium', 'force balance', 'moment balance', 'beam analytical solutions', 'pipe-stress analytical solutions', 'pressure-vessel/hoop-stress solutions', 'mesh convergence', 'time-step convergence for dynamics', 'modal-frequency benchmarks', 'natural-frequency verification', 'energy conservation', 'momentum-transfer consistency', 'fluid-force transfer verification', 'fatigue accumulation tests', 'crack-growth benchmark calculations', 'corrosion-growth calculations', 'erosion analytical benchmarks', 'reliability mathematics verification', 'monte-carlo convergence', 'probability normalization', 'limit-state-function verification', 'sensitivity/gradient verification',
    'corrosion coupon data', 'ultrasonic wall-thickness measurements', 'intelligent-pigging data', 'erosion measurements', 'vibration measurements', 'field vibration measurements'],
  referenceOnly: ['cauchy momentum equation', 'elasticity equations', 'shell equations', 'navier–stokes + structural dynamics', 'two-way fsi', 'cfd + fea', 'erosion + particle cfd', 'material condition', 'initial deposit thickness', 'guides', 'connections and equipment interfaces', 'allowable vibration or dynamic loading', 'pipe and equipment geometry',
    'material constitutive parameters', 'elastic/plastic properties', 'soil-pipe interaction', 'erosion threshold', 'wall-shear correlations', 'failure-probability models', 'fatigue parameters', 'support stiffness', 'slug-force coefficients', 'dynamic amplification', 'stress-concentration factors', 'defect-growth parameters', 'reliability distributions', 'consequence-model parameters', 'inspection histories', 'element-order convergence', 'strain-gauge measurements', 'accelerometers', 'pressure transducers', 'load cells', 'full-scale bend-force tests', 'slug-force experiments', 'fatigue tests', 'tensile/material tests', 'burst tests', 'collapse tests', 'failure databases', 'actual damage/failure observations'],
  equationsNote: 'Linear-elastic pipe and beam theory (no plasticity, shells or continuum finite elements); the fluid loads the structure one way. Slug loads use unit-cell closures or the flow suite\'s slug summary, not a slug-tracking solution. Corrosion is the de Waard–Milliams family for sweet service (not valid when H2S controls the film); erosion models assume dilute sand. Code formats (ASME B31.4/B31.8, DNV-ST-F101, DNV-RP-F101/F105/C203/O501, ISO 15156) are implemented from their published equations for screening and concept design; they do not replace a code-compliant design verification. The NORSOK M-506 corrosion model, fluid–structure feedback, shell buckling and computational fluid dynamics are not solved.',
  inputs: INPUTS,
  presets: [
    { name: 'Reference deep-water tie-back', values: {} },
    { name: 'Sour high-CO2 service', values: { co2: 8, h2s: 0.5, phAct: 4.2, inhibEff: 97, corrAllow: 6, wtMm: 23.8, snEnv: 'free', safetyClass: 'high', material: 'API 5L X65 sour service', corrModel: '1995' } },
    { name: 'Sandy late-life well', values: { sandKgD: 1500, sandUm: 350, ageY: 14, minWt: 13.6, erosionModel: 'governing', c14e: 125, inhibEff: 97, bends: [{ x: 0, angle: 90, radius: 0.38 }, { x: 18400, angle: 45, radius: 1.27 }, { x: 19560, angle: 90, radius: 0.38 }], geomFactor: 2, inspInterval: 2 } },
    { name: 'Slugging jumper fatigue', values: { slugMode: 'on', slugFreq: 0.08, slugLen: 25, slugVel: 6, slugHoldup: 0.9, spanLength: 22, spanX: 200, spanEnds: 'fixed-fixed', bendAt: 0.5, bendAngle: 90, odCoatMm: 0, damping: 1.5, snClass: 'F1', scf: 1.5, dlf: 2, supportCapacity: 40, nSlugs: 16 } },
    { name: 'Corroded line: fitness for service', values: { ageY: 12, minWt: 12.4, inhibEff: 97, code: 'b318', designFactor: 0.72, riserFactor: 0.72, designPressure: 250, inspInterval: 3, defects: [{ x: 2400, depth: 5.6, length: 180 }, { x: 7350, depth: 11.6, length: 650 }, { x: 11800, depth: 3.9, length: 900 }, { x: 15200, depth: 11.4, length: 60 }] } },
    { name: 'Shallow-water span with strong current', values: { waterDepth: 60, currentSpeed: 0.45, waveHeight: 2, wavePeriod: 8, spanLength: 24, spanX: 6000, spanEnds: 'fixed-pinned', damping: 1.2, odCoatMm: 365.8, coatDensity: 2400, designPressure: 150, pShutIn: 140, safetyClass: 'medium', snClass: 'F1', spanAxial: false } },
    { name: 'Buried onshore gas line (upheaval check)', values: { waterDepth: 0, coverDepth: 1, soilGamma: 10, imperfection: 0.25, code: 'b318', designFactor: 0.72, riserFactor: 0.72, designPressure: 150, pShutIn: 140, designTemp: 80, odCoatMm: 0, snEnv: 'air', currentSpeed: 0, spanLength: 12, tInstall: 15 } },
  ],
  pull: ({ fluid, site, outputs } = {}) => {
    const n = outputs?.net, f = outputs?.flow, s = outputs?.solids, o = outputs?.ops, sg = f?.slug, has = sg && sg.type && sg.type !== 'none', longest = n?.spanSource === 'terrain grid' && Array.isArray(n?.spans) && n.spans.length ? n.spans.reduce((a, b) => ((b?.length || 0) > (a?.length || 0) ? b : a)) : null, glycol = ['MEG', 'DEG', 'TEG'].includes(fluid?.inhibitor);
    return [
      { key: 'idMm', value: mm(n?.id), from: 'Network: inner diameter' }, { key: 'wtMm', value: mm(n?.wt), from: 'Network: wall thickness' }, { key: 'odCoatMm', value: mm(n?.od), from: 'Network: outer diameter with coatings' },
      { key: 'material', value: n?.material?.grade, from: 'Network: material grade' }, { key: 'smys', value: okNum(n?.material?.smys), from: 'Network: yield strength' }, { key: 'smts', value: okNum(n?.material?.smts), from: 'Network: tensile strength' }, { key: 'eMod', value: okNum(n?.material?.E), from: 'Network: Young\'s modulus' },
      { key: 'poisson', value: okNum(n?.material?.poisson), from: 'Network: Poisson\'s ratio' }, { key: 'alphaT', value: okNum(n?.material?.alphaT), from: 'Network: thermal expansion' }, { key: 'rhoSteel', value: okNum(n?.material?.rho), from: 'Network: steel density' },
      { key: 'designPressure', value: okNum(n?.designPressure), from: 'Network: design pressure' }, { key: 'designTemp', value: okNum(n?.designTemp), from: 'Network: design temperature' }, { key: 'waterDepth', value: okNum(n?.waterDepth), from: 'Network: water depth' }, { key: 'tInstall', value: okNum(n?.tSeabed), from: 'Network: seabed temperature' },
      { key: 'bends', value: Array.isArray(n?.bends) && n.bends.length ? n.bends.filter((b) => b && isNum(b.x)).slice(0, 40).map((b) => ({ x: b.x, angle: nz(b.angle, 90), radius: nz(b.radius, 0) })) : undefined, from: 'Network: bends' },
      { key: 'spanLength', value: longest && longest.length > 0 ? longest.length : undefined, from: 'Network: longest free span on the surveyed terrain' }, { key: 'spanX', value: longest && isNum(longest.x) ? longest.x : undefined, from: 'Network: location of the longest span' },
      { key: 'roughUm', value: isNum(s?.roughnessEff) && s.roughnessEff > 0 ? s.roughnessEff * 1e6 : isNum(n?.roughness) && n.roughness > 0 ? n.roughness * 1e6 : undefined, from: 'Solids / network: effective wall roughness' },
      { key: 'co2', value: okNum(fluid?.comp?.CO2), from: 'Case fluid: CO2 content (well stream, taken for the gas phase)' }, { key: 'h2s', value: okNum(fluid?.comp?.H2S), from: 'Case fluid: H2S content' }, { key: 'glycolWt', value: glycol ? okNum(fluid?.inhWt) : undefined, from: 'Case fluid: glycol in the water phase' },
      { key: 'slugFreq', value: has ? okNum(sg.freq) : undefined, from: 'Flow: slug frequency' }, { key: 'slugLen', value: has ? okNum(sg.length) : undefined, from: 'Flow: slug length' }, { key: 'slugVel', value: has ? okNum(sg.velocity) : undefined, from: 'Flow: slug velocity' }, { key: 'slugHoldup', value: has ? okNum(sg.holdupBody) : undefined, from: 'Flow: slug-body holdup' },
      { key: 'blockage', value: isNum(s?.blockage) ? clamp(s.blockage, 0, 1) : undefined, from: 'Solids: flow area lost' }, { key: 'plugX', value: okNum(s?.plugX ?? s?.onsetX), from: 'Solids: plug or onset location' }, { key: 'effIdMm', value: mm(s?.effectiveId), from: 'Solids: smallest effective bore' }, { key: 'plugProb', value: isNum(s?.plugProbability) ? clamp(s.plugProbability, 0, 1) : undefined, from: 'Solids: plug probability' },
      { key: 'depositRate', value: isNum(s?.waxRate) && s.waxRate >= 0 ? s.waxRate : undefined, from: 'Solids: wax deposition rate' },
      { key: 'blowdownMinT', value: okNum(o?.blowdownMinT), from: 'Operations: lowest blowdown temperature' }, { key: 'eventsShutdown', value: okNum(o?.eventsPerYear?.shutdowns), from: 'Operations: shutdowns per year' }, { key: 'eventsBlowdown', value: okNum(o?.eventsPerYear?.blowdowns), from: 'Operations: blowdowns per year' },
      { key: 'currentSpeed', value: okNum(site?.data?.currentSpeed), from: 'Site: current speed' }, { key: 'oilPrice', value: okNum(site?.data?.oilPrice), from: 'Site: oil price' },
    ].filter((it) => it.value !== undefined && it.value !== null && it.value !== '');
  },
  site: (site) => { const d = site?.data || {}; return [{ key: 'waterDepth', value: okNum(d.depth), from: 'Water depth at the site' }, { key: 'currentSpeed', value: okNum(d.currentSpeed), from: 'Current speed at the site' }, { key: 'waveHeight', value: okNum(d.waveHeight), from: 'Significant wave height at the site' }, { key: 'wavePeriod', value: okNum(d.wavePeriod), from: 'Wave period at the site' }, { key: 'tInstall', value: okNum(d.seabedTemp), from: 'Seabed temperature at the site' }, { key: 'oilPrice', value: okNum(d.oilPrice), from: 'Oil price' }].filter((it) => it.value !== undefined); },
  run,
  mesh: [
    { name: 'Beam elements on the span', keys: ['nElem'], min: 4, note: 'Hermite beam elements converge with the fourth power of the element length in the natural frequency; the dynamic stress converges with the second power.', metrics: [{ label: 'First natural frequency', unit: 'Hz', get: (r) => r.outputs.naturalFrequency }, { label: 'Peak dynamic bending stress', unit: 'MPa', get: (r) => r.outputs.peakDynamicStress }, { label: 'Static span stress', unit: 'MPa', get: (r) => r.outputs.staticSpanStress }] },
    { name: 'Time step of the dynamic response', keys: ['stepsPerCycle'], min: 8, note: 'Newmark-β with average acceleration is second-order accurate in the time step.', metrics: [{ label: 'Peak dynamic bending stress', unit: 'MPa', get: (r) => r.outputs.peakDynamicStress }, { label: 'Fatigue damage per year', unit: '1/y', get: (r) => r.outputs.fatigueDamagePerYear }] },
    { name: 'Monte Carlo sample size', keys: ['nMC'], min: 500, note: 'Sampling error falls with the square root of the sample size, so this study shows statistical scatter rather than a formal order of convergence.', metrics: [{ label: 'Failure probability, importance sampling', unit: '–', get: (r) => r.outputs.mcPof }, { label: 'Reliability index, importance sampling', unit: '–', get: (r) => r.outputs.mcBeta }] },
  ],
  calibration: {
    note: 'Fits the corrosion-model multiplier and the inhibitor efficiency to coupon or ultrasonic wall-loss rates (rows with and without inhibitor separate the two), the erosion multiplier to probe or wall-loss measurements at an elbow, and the damping ratio to measured vibration amplitudes near resonance.',
    params: [{ key: 'corrMult', label: 'Corrosion-model multiplier', lo: 0.1, hi: 4 }, { key: 'inhibEff', label: 'Inhibitor efficiency × availability (%)', lo: 30, hi: 99.9 }, { key: 'erosionMult', label: 'Erosion-model multiplier', lo: 0.05, hi: 20 }, { key: 'damping', label: 'Damping ratio (%)', lo: 0.2, hi: 15 }],
    columns: [{ key: 'cT', label: 'Temperature', unit: '°C' }, { key: 'cP', label: 'Pressure', unit: 'bara' }, { key: 'cPco2', label: 'CO2 partial pressure', unit: 'bar' }, { key: 'cVliq', label: 'Liquid velocity', unit: 'm/s' }, { key: 'cInh', label: 'Inhibitor on (1) or off (0)', unit: '–' }, { key: 'corrRate', label: 'Measured corrosion rate', unit: 'mm/y' },
      { key: 'cVmix', label: 'Mixture velocity at the elbow', unit: 'm/s' }, { key: 'cRho', label: 'Mixture density', unit: 'kg/m³' }, { key: 'cSand', label: 'Sand rate', unit: 'kg/d' }, { key: 'eroRate', label: 'Measured erosion rate', unit: 'mm/y' }, { key: 'cFreq', label: 'Excitation frequency', unit: 'Hz' }, { key: 'cForce', label: 'Force amplitude', unit: 'N' }, { key: 'vibAmp', label: 'Measured mid-span amplitude', unit: 'mm' }],
    targets: [{ key: 'corrRate', label: 'Corrosion rate', unit: 'mm/y' }, { key: 'eroRate', label: 'Erosion rate', unit: 'mm/y' }, { key: 'vibAmp', label: 'Vibration amplitude', unit: 'mm' }],
    model: calModel,
    sample: CAL_SAMPLE,
    validationSample: CAL_VALIDATION,
  },
  verify,
  live: { key: 'opLog', label: 'Operating log (time h, inlet pressure bara, inlet temperature °C)', help: 'Follow an exported historian file: its pressure and temperature cycles are rainflow-counted into the fatigue spectrum and every excursion above the limits is counted.' },
};
