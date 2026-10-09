// Suite 6 — Integrity, Loads, Risk & Engineering Assessment.
// Turns the pressure, temperature, slug and solids picture of the case into engineering consequences: pipe stress and
// code utilisation, collapse, buckle arrestors and managed lateral buckling, slug loads and span dynamics (beam finite
// elements with two-way coupling to the conveyed slug flow), continuum and shell finite elements of the pipe wall, fatigue
// and crack growth, CO2 corrosion and sand erosion (correlations and tracked particles), defect assessment, structural
// reliability (FORM/SORM/Monte Carlo) with an inspection-and-repair plan, and risk (fault tree, event tree, Markov,
// Bayesian network, FMECA, risk matrix). SI inside; MPa, bara, °C, mm at the interfaces.
// Literature constants are listed with their sources in PROVENANCE; measured reference data live in ../data/ref/integ.js.
import { clamp, linspace, logspace, sum, mean, isNum, interp1, brent, rng, lhs, histogram, fmt, nelderMead } from '../core/num.js';
import { fluidModel, waterContent } from '../core/thermo.js';
import { G, slugUnit, frictionFactor } from '../core/pipe.js';
import { flowPicture, caseLine } from '../core/caseflow.js';
import { BASE } from '../data/basecase.js';
import * as PROPS from '../core/props.js';
import * as REF from '../data/ref/integ.js';

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
// Vibration-velocity screening lines of process pipework: log10(v_rms, mm/s) = (log10 f + a)/b — see PROVENANCE.
const VIB_LINES = { concern: [0.48017, 2.127612], problem: [1.871083, 2.084547] };
const B318_T = [[121, 1], [149, 0.967], [177, 0.933], [204, 0.9], [232, 0.867]];
const DNV_SC = { low: { pc: 1.046, lb: 1.04, target: 1e-3 }, medium: { pc: 1.138, lb: 1.14, target: 1e-4 }, high: { pc: 1.308, lb: 1.26, target: 1e-5 }, veryHigh: { pc: 1.308, lb: 1.26, target: 1e-6 } };
const DNV_GM = 1.15;
const ALPHA_FAB = { seamless: 1, uo: 0.93, uoe: 0.85 }; // DNV-ST-F101 fabrication factor
const DNV_STRAIN = { low: 2, medium: 2.5, high: 3.3, veryHigh: 3.3 }; // strain resistance factor γ_ε of the displacement-controlled criterion
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
 * External-pressure capacities (Pa). o: { D, t, E, nu, fy, ovality (f0 = (Dmax − Dmin)/D), alphaFab, ovalityMin (design floor on f0, default 0.005; 0 when a measured out-of-roundness is to be used as it is) }
 * Returns { pel (elastic 2E(t/D)³/(1−ν²)), pp (plastic), pc (combined collapse with ovality, DNV-ST-F101 cubic), ppr (propagating buckle 35·fy·αfab·(t/D)^2.5) }.
 */
export function collapsePressure(o) {
  const D = +o.D, t = +o.t, fy = +o.fy, af = nz(o.alphaFab, 1), f0 = Math.max(nz(o.ovality, 0.005), o.ovalityMin ?? 0.005);
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
  const pc = collapsePressure({ D, t, E: o.E, nu: o.nu, fy, ovality: o.ovality, alphaFab: o.alphaFab }).pc, pex = (+o.pe || 0) - nz(o.pmin, 0);
  return { util: lin * lin + ((g * Math.max(pex, 0)) / pc) ** 2, mode: 'external overpressure', Mp, Sp, alphaC: ac, alphaP: 1, pc };
}
// Hobbs (1984) constants of the lateral modes: [name, k1, k2, k3, k4 (amplitude), k5 (bending moment)]
const HOBBS = [['Lateral mode 1', 80.76, 6.391e-5, 0.5, 2.407e-3, 0.06938], ['Lateral mode 2', 4 * Math.PI ** 2, 1.743e-4, 1, 5.532e-3, 0.1088], ['Lateral mode 3', 34.06, 1.668e-4, 1.294, 1.032e-2, 0.1434], ['Lateral mode 4', 28.2, 2.144e-4, 1.608, 1.047e-2, 0.1483]];
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
/**
 * Post-buckling state of a planned lateral buckle (Hobbs' solution with the feed-in limited by the spacing of the initiators).
 * Compatibility of every mode: geometric shortening k2·k3²·(μ_l·w/EI)²·L⁷ = 2·k3·L·(P0 − P)/EA + feed-in of the two slip zones, P = k1·EI/L²;
 * the slip zones give (P0 − P)²/(EA·μ_a·w) when they are free to develop (Hobbs' closed form) and less when they are cut at half the spacing.
 * o: { EI, EA, w (submerged weight N/m), muA, muL, P0 (fully restrained compressive effective force, N), spacing (m between initiators; Infinity = isolated buckle), ro (outer steel radius) }
 * Returns { modes: [{ name, L (buckle length parameter, m), P (force in the buckle), feedIn (m), slip (m each side), limited (bool), amplitude (m), moment (N·m), strain }], governing (largest strain) | null }.
 */
export function lateralBuckle(o) {
  const EI = +o.EI, EA = +o.EA, w = Math.max(+o.w, 1e-6), fa = Math.max(nz(o.muA, 0.5), 0.01), fl = Math.max(nz(o.muL, 0.5), 0.01), P0 = Math.max(+o.P0, 0), sp = o.spacing > 0 ? +o.spacing : Infinity, modes = [];
  for (const [name, k1, k2, k3, k4, k5] of HOBBS) {
    const Lmin = Math.sqrt((k1 * EI) / Math.max(P0, 1e-9)), half = Math.max(sp / 2, 1e-6);
    const parts = (L) => { const X = P0 - (k1 * EI) / (L * L), ls = X / (fa * w), lim = ls > half, feed = lim ? (2 * (X * half - 0.5 * fa * w * half * half)) / EA : (X * X) / (EA * fa * w); return { X, ls: Math.min(ls, half), lim, feed, geo: k2 * k3 * k3 * ((fl * w) / EI) ** 2 * L ** 7, inner: (2 * k3 * L * X) / EA }; };
    const fn = (L) => { const q = parts(L); return q.geo - q.inner - q.feed; };
    if (!(P0 > 0)) continue;
    const grid = logspace(Lmin * 1.0001, Lmin * 400, 500); let root = null;
    for (let i = grid.length - 1; i > 0; i--) if (fn(grid[i]) > 0 && fn(grid[i - 1]) <= 0) { root = brent(fn, grid[i - 1], grid[i], 1e-10 * grid[i]); break; }
    if (root === null) continue; // the restrained force is below the smallest force that sustains this mode
    const q = parts(root), M = k5 * fl * w * root * root;
    modes.push({ name, L: root, P: P0 - q.X, feedIn: q.geo, slip: q.ls, limited: q.lim, amplitude: (k4 * fl * w * root ** 4) / EI, moment: M, strain: (M * +o.ro) / EI });
  }
  return { modes, governing: modes.length ? modes.reduce((a, b) => (b.strain > a.strain ? b : a)) : null };
}
/**
 * Compressive strain capacity of the DNV-ST-F101 displacement-controlled local-buckling criterion: ε_c = 0.78·(t/D − 0.01)·(1 + 5.75·Δp/p_b)·α_h^−1.5·α_gw.
 * o: { D, t, fy, fu, dp (minimum internal overpressure, Pa), alphaH (yield / tensile ratio, default 0.93), alphaGw (girth-weld factor) }
 */
export function strainCapacity(o) {
  const D = +o.D, t = +o.t, pb = ((2 * t) / (D - t)) * Math.min(+o.fy, +o.fu / 1.15) * (2 / Math.sqrt(3)), agw = o.alphaGw ?? (D / t > 20 ? Math.max(1 - 0.01 * (D / t - 20), 0.6) : 1);
  return 0.78 * Math.max(t / D - 0.01, 0) * (1 + (5.75 * Math.max(nz(o.dp, 0), 0)) / pb) * pos(o.alphaH, 0.93) ** -1.5 * agw;
}
/**
 * Integral buckle arrestor (DNV-ST-F101): crossover pressure p_X = p_pr + (p_pr,BA − p_pr)·[1 − exp(−20·t2·L_BA/D²)], p_pr,BA = 35·f_y·α_fab·(t2/D)^2.5.
 * o: { D, t, t2 (arrestor wall), L (arrestor length), fy, alphaFab }. Returns { ppr, pprBA, pX } (Pa).
 */
export function arrestorCrossover(o) {
  const D = +o.D, af = nz(o.alphaFab, 1), ppr = 35 * o.fy * af * (o.t / D) ** 2.5, pba = 35 * o.fy * af * (o.t2 / D) ** 2.5;
  return { ppr, pprBA: pba, pX: ppr + (pba - ppr) * (1 - Math.exp((-20 * o.t2 * o.L) / (D * D))) };
}
/** Euler buckling load π²EI/(K·L)² (N). */
export const eulerLoad = (EI, L, K = 1) => (Math.PI ** 2 * EI) / (K * L) ** 2;

// ---- beam finite elements -----------------------------------------------------------------------------------------------
const END_K = { 'pinned-pinned': 1, 'fixed-fixed': 0.5, 'fixed-pinned': 0.699, 'fixed-free': 2, springs: 1, 'fixed-guided': 1, connection: 1 };
const endsOf = (kind) => ({ 'pinned-pinned': ['pinned', 'pinned'], 'fixed-fixed': ['fixed', 'fixed'], 'fixed-pinned': ['fixed', 'pinned'], 'fixed-free': ['fixed', 'free'], springs: ['spring', 'spring'], 'fixed-guided': ['fixed', 'guided'], connection: ['flange', 'flange'] })[kind] || ['pinned', 'pinned'];
/**
 * Beam model with two-node Hermite elements (deflection and rotation at every node).
 * o: { L, EI, m (kg/m incl. content and added mass), n (elements), ends: ['pinned' | 'fixed' | 'free' | 'spring' | 'guided' (rotation held, deflection free) | 'flange' (deflection held, rotational spring kR: a connection of finite stiffness), …],
 *      kT: [N/m, N/m], kR: [N·m/rad, N·m/rad] (end springs),
 *      supports: [x | { x, k, type: 'support' (rigid, or a spring when k > 0) | 'guide' (lateral spring k) | 'anchor' (deflection and rotation held) | 'mass' (lumped mass, kg: flange, valve, connector), mass }], nodesAt: [x] (extra nodes), N (axial force, tension +),
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
    else if (kind === 'guided') fixedSet.add(2 * node + 1);
    else if (kind === 'flange') { fixedSet.add(2 * node); const kr = nz(o.kR?.[i], 0); if (kr > 0) springs.push({ dof: 2 * node + 1, k: kr }); }
    else if (kind === 'spring') { const kt = nz(o.kT?.[i], 0), kr = nz(o.kR?.[i], 0); if (kt > 0) springs.push({ dof: 2 * node, k: kt }); if (kr > 0) springs.push({ dof: 2 * node + 1, k: kr }); }
  });
  const lumped = [];
  for (const s of sup) {
    let node = 0; for (let i = 1; i < nd; i++) if (Math.abs(x[i] - s.x) < Math.abs(x[node] - s.x)) node = i;
    if (s.type === 'mass') { if (s.mass > 0) { M[2 * node][2 * node] += +s.mass; lumped.push({ node, mass: +s.mass }); } }
    else if (s.type === 'anchor') { fixedSet.add(2 * node); fixedSet.add(2 * node + 1); }
    else if (s.k > 0) springs.push({ dof: 2 * node, k: +s.k }); else fixedSet.add(2 * node);
  }
  for (const s of springs) K[s.dof][s.dof] += s.k;
  const fixed = [...fixedSet].sort((a, b) => a - b), free = []; for (let i = 0; i < ndof; i++) if (!fixedSet.has(i)) free.push(i);
  return { L, EI, m, x, ndof, K, M, free, fixed, springs, elems, lumped, timoshenko: kGA > 0 };
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
 * o: { rho, v, Dmm (outer diameter), tmm (wall), support: 'stiff' | 'mediumStiff' | 'medium' | 'flexible', fvf (fluid viscosity factor, 1 for multiphase) }
 */
export function fivScreen(o) {
  const ke = +o.rho * o.v * o.v, D = +o.Dmm, lnD = Math.log(D);
  const [al, be] = o.support === 'stiff' ? [446187 + 646 * D + 9.17e-4 * D ** 3, 0.1 * lnD - 1.3739] : o.support === 'medium' ? [150412 + 209 * D, 0.0815 * lnD - 1.3269] : o.support === 'flexible' ? [41.21 * D + 49397, 0.0815 * lnD - 1.3842] : [283921 + 370 * D, 0.1106 * lnD - 1.501];
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
// DNV-RP-C203 (2016 edition) S–N curves in air: N = 10^(a − m·log S), slope m1 up to 10⁷ cycles, m2 beyond; k = thickness exponent (reference thickness 25 mm).
// Seawater with cathodic protection (Table 2-2): log a1 lower by 0.400 (0.200 for B1 and B2) up to 10⁶ cycles, then the air curve. Free corrosion (Table 2-4): one slope m = 3, log a below.
const SN_FREE = Object.freeze({ B1: 12.436, B2: 12.262, C: 12.115, C1: 11.972, C2: 11.824, D: 11.687, E: 11.533, F: 11.378, F1: 11.222, F3: 11.068, G: 10.921, W1: 10.784, W2: 10.63, W3: 10.493 });
/** Standard deviation of log N behind the design curves (mean minus two standard deviations). */
export const SN_SD = 0.2;
export const SN_CURVES = Object.freeze({
  B1: { m1: 4, a1: 15.117, m2: 5, a2: 17.146, k: 0 }, B2: { m1: 4, a1: 14.885, m2: 5, a2: 16.856, k: 0 }, C: { m1: 3, a1: 12.592, m2: 5, a2: 16.32, k: 0.05 }, C1: { m1: 3, a1: 12.449, m2: 5, a2: 16.081, k: 0.1 },
  C2: { m1: 3, a1: 12.301, m2: 5, a2: 15.835, k: 0.15 }, D: { m1: 3, a1: 12.164, m2: 5, a2: 15.606, k: 0.2 }, E: { m1: 3, a1: 12.01, m2: 5, a2: 15.35, k: 0.2 }, F: { m1: 3, a1: 11.855, m2: 5, a2: 15.091, k: 0.25 },
  F1: { m1: 3, a1: 11.699, m2: 5, a2: 14.832, k: 0.25 }, F3: { m1: 3, a1: 11.546, m2: 5, a2: 14.576, k: 0.25 }, G: { m1: 3, a1: 11.398, m2: 5, a2: 14.33, k: 0.25 }, W1: { m1: 3, a1: 11.261, m2: 5, a2: 14.101, k: 0.25 },
  W2: { m1: 3, a1: 11.107, m2: 5, a2: 13.845, k: 0.25 }, W3: { m1: 3, a1: 10.97, m2: 5, a2: 13.617, k: 0.25 },
});
/**
 * Cycles to failure for a stress range S (MPa). env: 'air' | 'cp' (seawater with cathodic protection: a1 − 0.4, or − 0.2 for B1 and B2, knee at 10⁶) | 'free' (free corrosion: single slope m = 3, tabulated log a).
 * opt: { t (mm wall, thickness correction above 25 mm), scf }.
 */
export function snCycles(S, cls = 'D', env = 'air', opt = {}) {
  const c = SN_CURVES[cls] || SN_CURVES.D, t = nz(opt.t, 25), kt = env === 'free' && (cls === 'C' || cls === 'C1') ? 0.15 : c.k, s = Math.abs(S) * nz(opt.scf, 1) * (t > 25 ? (t / 25) ** kt : 1);
  if (!(s > 0)) return Infinity;
  const ls = Math.log10(s);
  const sh = nz(opt.shift, 0);
  if (env === 'free') return 10 ** (sh + (SN_FREE[cls] ?? SN_FREE.D) - 3 * ls);
  const a1 = env === 'cp' ? c.a1 - (c.m1 === 4 ? 0.2 : 0.4) : c.a1, knee = env === 'cp' ? 1e6 : 1e7, n1 = 10 ** (a1 - c.m1 * ls);
  return 10 ** sh * (n1 <= knee ? n1 : 10 ** (c.a2 - c.m2 * ls));
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
  const Fscale = o.scale === false ? 1 : Math.min(1, 10 ** (2400 / T - 0.6 * lf - 6.7)), tScale = 2400 / (6.7 + 0.6 * lf) - 273.15, W = clamp(100 - nz(o.glycolWt, 0), 1, 100), Fglycol = W < 5 ? 0.008 : 10 ** (1.6 * (Math.log10(W) - 2));
  let base, Vr = 0, Vm = 0;
  if (o.model === '1991') base = 10 ** (5.8 - 1710 / T + 0.67 * lf);
  else { Vr = 10 ** (4.93 - 1119 / T + 0.58 * lf - 0.34 * (pH - pHco2)); Vm = (2.45 * Math.max(nz(o.U, 1), 0.01) ** 0.8 * f) / pos(o.d, 0.1) ** 0.2; base = 1 / (1 / Vr + 1 / Vm); }
  return { rate: base * Fscale * Fglycol * (1 - clamp(nz(o.inhibEff, 0), 0, 1)) * nz(o.mult, 1), base, fCO2: f, fugacity: fug, pHco2, pH, Fscale, tScale, Fglycol, Vr, Vm };
}
// NORSOK M-506 (Rev. 2, 2005): K_t and the pH function f(pH)_t at the tabulated temperatures; [pH from, pH to, kind (0 polynomial, 1 exponential), coefficients]
const M506_KT = [[5, 0.42], [15, 1.59], [20, 4.762], [40, 8.927], [60, 10.695], [80, 9.949], [90, 6.25], [120, 7.77], [150, 5.203]];
const M506_FPH = { 5: [[3.5, 4.6, 0, 2.0676, -0.2309, 0, 0], [4.6, 6.5, 0, 4.342, -1.051, 0.0708, 0]], 15: [[3.5, 4.6, 0, 2.0676, -0.2309, 0, 0], [4.6, 6.5, 0, 4.986, -1.191, 0.0708, 0]], 20: [[3.5, 4.6, 0, 2.0676, -0.2309, 0, 0], [4.6, 6.5, 0, 5.1885, -1.2353, 0.0708, 0]], 40: [[3.5, 4.6, 0, 2.0676, -0.2309, 0, 0], [4.6, 6.5, 0, 5.1885, -1.2353, 0.0708, 0]],
  60: [[3.5, 4.6, 0, 1.836, -0.1818, 0, 0], [4.6, 6.5, 0, 15.444, -6.1291, 0.8204, -0.0371]], 80: [[3.5, 4.6, 0, 2.6727, -0.3636, 0, 0], [4.6, 6.5, 1, 331.68, -1.2618]], 90: [[3.5, 4.57, 0, 3.1355, -0.4673, 0, 0], [4.57, 5.62, 1, 21254, -2.1811], [5.62, 6.5, 0, 0.4014, -0.0538, 0, 0]],
  120: [[3.5, 4.3, 0, 1.5375, -0.125, 0, 0], [4.3, 5, 0, 5.9757, -1.157, 0, 0], [5, 6.5, 0, 0.546125, -0.071225, 0, 0]], 150: [[3.5, 3.8, 0, 1, 0, 0, 0], [3.8, 5, 0, 17.634, -7.0945, 0.715, 0], [5, 6.5, 0, 0.037, 0, 0, 0]] };
/** Wall shear stress of NORSOK M-506 (Pa): S = ½·ρ·f·u², f = 0.001375·[1 + (20000·k/D + 10⁶·μ/(ρ·u·D))^0.33]. */
export const norsokShear = (rho, mu, u, D, k = 50e-6) => 0.5 * rho * 0.001375 * (1 + ((20000 * k) / D + (1e6 * mu) / (rho * Math.max(u, 1e-9) * D)) ** 0.33) * u * u;
/**
 * pH of water in equilibrium with CO2 (no iron-carbonate saturation), charge balance in the form of NORSOK M-506: C_H³ + C_bic·C_H² − (K_H·K_1·p + K_W)·C_H − 2·K_H·K_1·K_2·p = 0,
 * with the temperature, pressure and ionic-strength functions of its equations 14–18. Two departures from the printed standard, both needed to obtain physical values: the first
 * dissociation constant uses 1684915/T² (the printed 168491.5 has lost a digit; pK₁ = 6.35 at 25 °C), and because that K₁ already refers to all dissolved CO2 the hydration constant 0.00258 is not applied again.
 * o: { T (°C), pCO2 (bar, fugacity), P (bar total), bicarb (mg/L as HCO3⁻), ionic (mol/L) }.
 */
export function norsokPH(o) {
  const Tc = nz(o.T, 20), T = Tc + 273.15, Tf = Tc * 1.8 + 32, P = nz(o.P, 1) * 14.5038, I = Math.max(nz(o.ionic, 0), 0), p = Math.max(nz(o.pCO2, 1), 1e-9), cb = Math.max(nz(o.bicarb, 0), 0) / 61017;
  const KH = (Tc <= 80 ? 55.5084 * Math.exp(-(4.8 + 3934.4 / T - 941290.2 / (T * T))) : 55.5084 * Math.exp(-((1713.53 * (1 - T / 647) ** (1 / 3)) / T + 3.875 + 3680.09 / T - 1198506.1 / (T * T)))) * 10 ** -(1.79e-4 * P + 0.107 * I);
  const K1 = 10 ** -(356.3094 + 0.06091964 * T - 21834.37 / T - 126.8339 * Math.log10(T) + 1684915 / (T * T) - 2.564e-5 * P - 0.491 * Math.sqrt(I) + 0.379 * I - 0.06506 * I ** 1.5 - 1.458e-3 * I * Tf);
  const K2 = 10 ** -(107.8871 + 0.03252849 * T - 5151.79 / T - 38.92561 * Math.log10(T) + 563713.9 / (T * T) - 2.118e-5 * P - 1.255 * Math.sqrt(I) + 0.867 * I - 0.174 * I ** 1.5 - 1.588e-3 * Tf * I), KW = 10 ** -(29.3868 - 0.0737549 * T + 7.47881e-5 * T * T);
  const a = KH * K1 * p, fn = (lh) => { const h = 10 ** lh; return h ** 3 + cb * h * h - (a + KW) * h - 2 * a * K2; };
  return { pH: -brent(fn, -14, 1, 1e-12), KH, K1, K2, KW };
}
/**
 * NORSOK M-506 (Rev. 2) CO2 corrosion rate of carbon steel (mm/y): CR_t = K_t·f_CO2^0.62·(S/19)^(0.146 + 0.0324·log f_CO2)·f(pH)_t for 20–150 °C,
 * exponent 0.36 at 15 °C, and no shear term at 5 °C; between the tabulated temperatures the rates (not the constants) are interpolated linearly.
 * o: { T (°C), pCO2 (bar), P (bar total, for the fugacity), pH, S (wall shear stress, Pa), glycolWt, inhibEff (0–1), fugacity (false = pCO2 is already a fugacity) }
 * Returns { rate (with the larger of the glycol and inhibitor reductions, as the standard prescribes), base, fCO2, pH, S, inRange (all inputs inside the validity limits), clipped: [names] }.
 */
export function norsokM506(o) {
  const Tc = nz(o.T, 60), T = Tc + 273.15, p = Math.max(nz(o.pCO2, 0), 0), clipped = [];
  if (!(p > 0)) return { rate: 0, base: 0, fCO2: 0, pH: nz(o.pH, 7), S: nz(o.S, 0), inRange: true, clipped };
  const f0 = o.fugacity === false ? p : p * 10 ** (Math.min(Math.max(nz(o.P, p), p), 250) * (0.0031 - 1.4 / T)), lim = (v, lo, hi, name) => { if (v < lo || v > hi) clipped.push(name); return clamp(v, lo, hi); };
  const f = lim(f0, 0.1, 10, 'CO2 fugacity'), S = lim(pos(o.S, 19), 1, 150, 'wall shear stress'), pH = lim(pos(o.pH, 4), 3.5, 6.5, 'pH'), Tt = lim(Tc, 5, 150, 'temperature'), lf = Math.log10(f);
  const fpH = (t) => { const rows = M506_FPH[t], r = rows.find((q) => pH <= q[1]) || rows[rows.length - 1]; return r[2] === 1 ? r[3] * Math.exp(r[4] * pH) : r[3] + r[4] * pH + r[5] * pH * pH + r[6] * pH ** 3; };
  const at = ([t, K]) => (t === 5 ? K * f ** 0.36 : K * f ** (t === 15 ? 0.36 : 0.62) * (S / 19) ** (0.146 + 0.0324 * lf)) * fpH(t);
  let i = 0; while (i < M506_KT.length - 2 && Tt > M506_KT[i + 1][0]) i++;
  const a = M506_KT[i], b = M506_KT[i + 1], base = Math.max(at(a) + ((at(b) - at(a)) * (Tt - a[0])) / (b[0] - a[0]), 0), g = clamp(nz(o.glycolWt, 0), 0, 100), Fg = g > 95 ? 0.008 : 10 ** (1.6 * (Math.log10(100 - g) - 2));
  return { rate: base * Math.min(Fg, 1 - clamp(nz(o.inhibEff, 0), 0, 1)) * nz(o.mult, 1), base, fCO2: f0, pH, S, inRange: !clipped.length, clipped, Fglycol: Fg };
}
// Salt effect on the charge-transfer reactions: reference exchange current densities fitted to sweeps at 1, 3, 10 and 20 wt % NaCl, normalised by the 1 wt % values; activation
// enthalpies (J/mol) that rise below 20 °C — see PROVENANCE. [temperature knots °C], then one row per reaction.
const SALT_I0 = { w: [1, 3, 10, 20], cathodic: [1, 1, 0.75, 0.375], anodic: [1, 0.875, 0.325, 0.225] };
const LOWT_DH = { T: [1, 5, 20], H: [140e3, 100e3, 30e3], C: [140e3, 100e3, 50e3], Fe: [215e3, 185e3, 37.5e3], W: [140e3, 120e3, 30e3] };
/** Retardation of the electrode reactions and activity coefficient of H⁺ in NaCl brine (wt % NaCl): { cathodic, anodic (factors on the exchange current densities), gammaH }. Log-linear between the tabulated concentrations, constant below 1 and above 20 wt %. */
export function saltRetardation(wt) {
  const w = clamp(nz(wt, 0), 0, 26), f = (arr) => (w <= 1 ? 1 : Math.exp(interp1(SALT_I0.w, arr.map(Math.log), Math.min(w, 20))));
  return { cathodic: f(SALT_I0.cathodic), anodic: f(SALT_I0.anodic), gammaH: 3 ** (Math.min(w, 26) / 20) };
}
/** Arrhenius factor i₀(T)/i₀(T_ref) = exp(−ΔH(T)/R·(1/T − 1/T_ref)) with the activation enthalpy tabulated against temperature (piecewise linear between the knots, constant outside them). */
function arrheniusVar(T, Tref, dH) {
  const h = interp1(LOWT_DH.T, dH, clamp(T - 273.15, LOWT_DH.T[0], LOWT_DH.T[LOWT_DH.T.length - 1]));
  return Math.exp((-h / RGAS) * (1 / T - 1 / Tref));
}
/**
 * Mixed-potential model of iron in CO2-saturated water or brine: anodic iron dissolution against the reduction of H⁺ (mass-transfer limited), of carbonic acid (limited by the
 * hydration of CO2) and of water, Tafel kinetics, solved for the corrosion potential (kinetic constants after Nesic et al., 1996, with the later activation enthalpies).
 * Additions for cold and saline service: activation enthalpies that rise below 20 °C (integrated over temperature), retardation of the charge-transfer reactions by chloride,
 * brine density and viscosity in the mass transfer, salting-out of CO2 and the activity coefficient of H⁺ (the pH is the measured one, an activity).
 * o: { T (°C), pCO2 (bar, fugacity), pH, geometry: 'pipe' (U m/s, d m) | 'rce' (rotating cylinder: rpm, d), rho, mu (liquid; default water or brine at T), ionic (mol/L, default from the salt), salt (wt % NaCl),
 *      lowT (false = constant activation enthalpies), water (false = no water reduction), steel: i₀ of iron at 25 °C (A/m², default 1) }
 * Returns { Ecorr (V), icorr (A/m²), rate (mm/y), iLimH, iLimH2CO3 (A/m²), rateLimit (mm/y, mass-transfer / reaction limit), Sh, km (m/s), Re, pH, tafel: { ba, bc }, shear (Pa, rotating cylinder or smooth pipe), salt: { cathodic, anodic, gammaH } }.
 */
export function mixedPotential(o) {
  const Tc = nz(o.T, 20), T = Tc + 273.15, p = Math.max(nz(o.pCO2, 1), 1e-6), pH = pos(o.pH, 3.71 + 0.00417 * Tc - 0.5 * Math.log10(p)), d = pos(o.d, o.geometry === 'rce' ? 0.012 : 0.1), wt = clamp(nz(o.salt, 0), 0, 26), Sg = 10 * wt;
  const rho = pos(o.rho, PROPS.density(clamp(Tc, 0, 150), Math.min(Sg, 260))), mu = pos(o.mu, PROPS.viscosity(clamp(Tc, 0, 150), Math.min(Sg, 260))), sr = saltRetardation(wt), ionic = o.ionic !== undefined && o.ionic !== null ? nz(o.ionic, 0) : ((wt / 100) * rho) / 58.44;
  const U = o.geometry === 'rce' ? (Math.PI * d * pos(o.rpm, 1000)) / 60 : Math.max(nz(o.U, 1), 0.01);
  const Tf = Tc * 1.8 + 32, ksol = (14.5 / 1.00258) * 10 ** -(2.27 + 5.65e-3 * Tf - 8.06e-6 * Tf * Tf + 0.075 * ionic), cCO2 = ksol * p * 1000, cH = (10 ** -pH / sr.gammaH) * 1000, cH2CO3 = 2.58e-3 * cCO2; // mol/m³
  const dScale = (T / 298.15) * (8.9e-4 / mu), DH = 9.31e-9 * dScale, DC = 1.3e-9 * dScale, Re = (rho * U * d) / mu, Sc = mu / (rho * DH), Sh = o.geometry === 'rce' ? 0.0791 * Math.max(Re, 1) ** 0.7 * Sc ** 0.356 : sherwood(Re, Sc), km = (Sh * DH) / d;
  const iLimH = km * FARADAY * cH, khyd = 10 ** (329.85 - 110.541 * Math.log10(T) - 17265.4 / T), iLimC = FARADAY * cCO2 * Math.sqrt(DC * 2.58e-3 * khyd);
  const low = o.lowT !== false, arr = (key, dH, Tref) => (low ? arrheniusVar(T, Tref, LOWT_DH[key]) : Math.exp((-dH / RGAS) * (1 / T - 1 / Tref))), b = (2.303 * RGAS * T) / FARADAY, ba = b / 1.5, bc = b / 0.5, eH = -b * pH, eFe = -0.488;
  const i0H = 0.05 * (cH / 0.1) ** 0.5 * arr('H', 30e3, 298.15) * sr.cathodic, i0C = 0.06 * (cH / 0.01) ** -0.5 * (cH2CO3 / 0.1) * arr('C', 50e3, 293.15) * sr.cathodic, i0W = o.water === false ? 0 : 3e-5 * arr('W', 30e3, 293.15) * sr.cathodic, i0Fe = pos(o.steel, 1) * arr('Fe', 37.5e3, 298.15) * sr.anodic;
  const ia = (E) => i0Fe * 10 ** ((E - eFe) / ba), ic = (E) => 1 / (1 / (i0H * 10 ** (-(E - eH) / bc)) + 1 / iLimH) + 1 / (1 / (i0C * 10 ** (-(E - eH) / bc)) + 1 / iLimC) + i0W * 10 ** (-(E - eH) / bc);
  const Ecorr = brent((E) => ia(E) - ic(E), -1.5, 0.3, 1e-12), icorr = ia(Ecorr), shear = o.geometry === 'rce' ? 0.0791 * Math.max(Re, 1) ** -0.3 * rho * U * U : 0.5 * rho * U * U * 0.046 * Math.max(Re, 1) ** -0.2;
  return { Ecorr, icorr, rate: faradayRate(icorr), iLimH, iLimH2CO3: iLimC, rateLimit: faradayRate(iLimH + iLimC), Sh, km, Re, pH, tafel: { ba, bc }, shear, salt: sr, U, rho, mu };
}
/**
 * CO2 corrosion rate of carbon steel with a model that follows the conditions (mm/y, bare steel before glycol and inhibitor are applied in `rate`).
 * 'auto': NORSOK M-506 inside its range (5–150 °C, pH 3.5–6.5, CO2 fugacity 0.1–10 bar, wall shear 1–150 Pa; its own equations at 5 and 15 °C). Outside the range, and for salinity
 * above 1 wt % NaCl, the NORSOK rate at the nearest valid condition is scaled with the ratio of the mixed-potential model between the actual and that condition — so colder water,
 * strong brine, or a pH or fugacity beyond the limits follow the electrochemistry while the result stays continuous with the standard.
 * o: { model: 'auto' | 'norsok' | '1995' | '1991' | 'mechanistic', T (°C), pCO2 (bar), P (bar), pH (0 or omitted: water in equilibrium with CO2 and the bicarbonate), bicarb (mg/L), salt (wt % NaCl),
 *      S (wall shear, Pa) or rho, mu, U, d, rough (for the NORSOK shear), geometry / rpm (mixed-potential mass transfer), glycolWt, inhibEff (0–1), mult }
 * Returns { rate, base (before glycol, inhibitor and multiplier), model (text), norsok (rate at the clipped condition), factor (mixed-potential ratio), fT, fSalt, fChem (its parts: temperature, salinity, pH and fugacity), pH, S, fCO2, inRange, Fglycol, mech }.
 */
export function co2Corrosion(o) {
  const Tc = nz(o.T, 60), T = Tc + 273.15, p = Math.max(nz(o.pCO2, 0), 0), wt = clamp(nz(o.salt, 0), 0, 26), g = clamp(nz(o.glycolWt, 0), 0, 100), Fg = g > 95 ? 0.008 : 10 ** (1.6 * (Math.log10(100 - g) - 2)), red = Math.min(Fg, 1 - clamp(nz(o.inhibEff, 0), 0, 1)) * nz(o.mult, 1), model = o.model || 'auto';
  if (!(p > 0)) return { rate: 0, base: 0, model, norsok: 0, factor: 1, fT: 1, fSalt: 1, fChem: 1, pH: nz(o.pH, 7), S: nz(o.S, 0), fCO2: 0, inRange: true, Fglycol: Fg };
  const f = p * 10 ** (Math.min(Math.max(nz(o.P, p), p), 250) * (0.0031 - 1.4 / T)), ionic = ((wt / 100) * 1000) / 58.44, pH = o.pH > 0 ? +o.pH : norsokPH({ T: Tc, pCO2: f, P: nz(o.P, p), bicarb: o.bicarb, ionic }).pH;
  const hyd = { geometry: o.geometry, rpm: o.rpm, U: o.U, d: o.d }, S = o.S > 0 ? +o.S : o.geometry === 'rce' ? mixedPotential({ T: Tc, pCO2: f, pH, salt: wt, ...hyd }).shear : o.U > 0 ? norsokShear(pos(o.rho, 1000), pos(o.mu, 1e-3), +o.U, pos(o.d, 0.1), pos(o.rough, 50e-6)) : 19;
  if (model === '1991' || model === '1995') { const dw = deWaardMilliams({ T: Tc, pCO2: p, P: o.P, model, U: o.U, d: o.d, pH: o.pH, glycolWt: 0 }); return { rate: dw.base * dw.Fscale * red, base: dw.base * dw.Fscale, model: `de Waard (${model})`, norsok: null, factor: 1, fT: 1, fSalt: 1, fChem: 1, pH: dw.pH, S, fCO2: dw.fCO2, inRange: true, Fglycol: Fg, dw }; }
  if (model === 'mechanistic') { const m = mixedPotential({ T: Tc, pCO2: f, pH, salt: wt, ...hyd }); return { rate: m.rate * red, base: m.rate, model: 'mixed-potential', norsok: null, factor: 1, fT: 1, fSalt: 1, fChem: 1, pH, S, fCO2: f, inRange: true, Fglycol: Fg, mech: m }; }
  const ns = norsokM506({ T: Tc, pCO2: f, fugacity: false, pH, S });
  if (model === 'norsok') return { rate: ns.base * red, base: ns.base, model: 'NORSOK M-506', norsok: ns.base, factor: 1, fT: 1, fSalt: 1, fChem: 1, pH, S, fCO2: f, inRange: ns.inRange, Fglycol: Fg };
  // automatic: NORSOK at the clipped condition times the mixed-potential ratios, one factor at a time so that each is reported
  const Tr = clamp(Tc, 5, 150), pHr = clamp(pH, 3.5, 6.5), fr = clamp(f, 0.1, 10), wr = Math.min(wt, 1), mp = (t, ph, ff, w) => mixedPotential({ T: t, pCO2: ff, pH: ph, salt: w, ...hyd }).rate;
  const m0 = mp(Tr, pHr, fr, wr), fT = Tr !== Tc ? mp(Tc, pHr, fr, wr) / m0 : 1, fChem = pHr !== pH || fr !== f ? mp(Tr, pH, f, wr) / m0 : 1, fSalt = wt > wr ? mp(Tr, pHr, fr, wt) / m0 : 1, factor = fT * fChem * fSalt, base = ns.base * factor;
  return { rate: base * red, base, model: factor === 1 ? 'NORSOK M-506' : `NORSOK M-506 × mixed-potential ratio (${[fT !== 1 ? 'temperature' : '', fSalt !== 1 ? 'salinity' : '', fChem !== 1 ? 'pH / fugacity' : ''].filter(Boolean).join(', ')})`, norsok: ns.base, factor, fT, fSalt, fChem, pH, S, fCO2: f, inRange: factor === 1 && ns.inRange, Fglycol: Fg };
}
/** Sour-service check in the ISO 15156-2 format: region 0 below 0.3 kPa H2S, regions 1–3 from the in-situ pH (boundaries digitised from the severity diagram). */
export function sourRegion(pH2SkPa, pH) {
  if (!(pH2SkPa >= 0.3)) return { region: 0, sour: false, label: 'Region 0 (below 0.3 kPa H2S: no special requirements)' };
  const lp = Math.log10(pH2SkPa), region = pH >= 4.5 + clamp(lp, -1, 2) ? 1 : pH >= 3.5 + clamp(lp, 0, 2) ? 2 : 3; // corner points of Figure 1 of ISO 15156-2
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
 * 'modified': 0.85·d·L area, flow stress SMYS + 69 MPa, two-term Folias factor; 'dnv': DNV-RP-F101 single defect, Q = √(1 + 0.31 (L/√(Dt))²) (failure pressure of the allowable-stress format);
 * 'dnvCap': the same with the factor 1.05 of the best-estimate burst capacity (use it with the measured tensile strength when comparing with burst tests).
 * Returns { pf, intact (Pa), ratio, M (bulging factor), flow (Pa) }.
 */
export function b31g(o) {
  const D = +o.D, t = +o.t, r = clamp(o.d / t, 0, 1), z = (o.L * o.L) / (D * t);
  need(D > 0 && t > 0, 'Defect assessment needs a positive diameter and wall thickness.');
  if (o.method === 'dnv' || o.method === 'dnvCap') { const Q = Math.sqrt(1 + 0.31 * z), p0 = ((o.method === 'dnvCap' ? 1.05 : 1) * (2 * t * o.smts)) / (D - t), pf = r >= 1 ? 0 : (p0 * (1 - r)) / (1 - r / Q); return { pf, intact: p0, ratio: pf / p0, M: Q, flow: +o.smts }; }
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
// DNV-RP-O501 steel constants, Oka et al. (2005) constants for silica sand (Hv in GPa, mm³/kg) — see PROVENANCE.
const O501 = { K: 2e-9, n: 2.6, rhoT: 7800, C1: 2.5 };
const dnvAngle = (a) => { const sa = Math.sin(a); return 0.6 * (sa + 7.2 * (sa - sa * sa)) ** 0.6 * (1 - Math.exp(-20 * a)); };
const OKA = { s1: 0.71, q1: 0.14, s2: 2.4, q2: -0.94, K: 65, k1: -0.12, k2a: 2.3, k2b: 0.038, k3: 0.19, vRef: 104, dRef: 326e-6 };
/** API RP 14E erosional velocity (m/s) for a mixture density (kg/m³) and the empirical C-factor in field units (100 continuous, 125 intermittent, up to 150–200 for corrosion-resistant service). */
export const erosionalVelocity = (rho, C = 100) => (1.2199 * C) / Math.sqrt(Math.max(rho, 1e-6)); // 0.3048 m/ft × √(16.018463 kg/m³ per lb/ft³)
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
    const fa = Math.tan(alpha) <= 1 / 3 ? Math.sin(2 * alpha) - 3 * sa * sa : Math.cos(alpha) ** 2 / 3, vol = (nz(o.finnieC, 0.5) * U * U * fa) / (4 * pos(o.flowStress, 1.96e9)); // m³ per kg of sand
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

// ---- continuum finite elements: plane stress, plane strain and axisymmetric solids --------------------------------------
// Linear (3-node) and quadratic (6-node, isoparametric) triangles; symmetric banded Cholesky solve. Stresses are
// [σxx, σyy, σtt, τxy]: in the axisymmetric mode x = r, y = z and σtt is the hoop stress; in the plane modes σtt is the out-of-plane stress.
const TRI3 = [[2 / 3, 1 / 6, 1 / 6, 1 / 6], [1 / 6, 2 / 3, 1 / 6, 1 / 6], [1 / 6, 1 / 6, 2 / 3, 1 / 6]]; // [L1, L2, L3, weight] on the reference triangle (area 1/2)
const TRI7 = (() => { const a1 = 0.0597158717897698, b1 = 0.4701420641051151, a2 = 0.7974269853530873, b2 = 0.1012865073234563, w1 = 0.1323941527885062 / 2, w2 = 0.1259391805448271 / 2; return [[1 / 3, 1 / 3, 1 / 3, 0.1125], [a1, b1, b1, w1], [b1, a1, b1, w1], [b1, b1, a1, w1], [a2, b2, b2, w2], [b2, a2, b2, w2], [b2, b2, a2, w2]]; })();
const GAUSS3 = [[-Math.sqrt(0.6), 5 / 9], [0, 8 / 9], [Math.sqrt(0.6), 5 / 9]];
/** Shape functions and their natural derivatives of the 3- or 6-node triangle at area coordinates (L1, L2, L3); ξ = L2, η = L3. */
function triShape(nn, L1, L2, L3) {
  if (nn === 3) return { N: [L1, L2, L3], dxi: [-1, 1, 0], deta: [-1, 0, 1] };
  return { N: [L1 * (2 * L1 - 1), L2 * (2 * L2 - 1), L3 * (2 * L3 - 1), 4 * L1 * L2, 4 * L2 * L3, 4 * L3 * L1],
    dxi: [-(4 * L1 - 1), 4 * L2 - 1, 0, 4 * (L1 - L2), 4 * L3, -4 * L3], deta: [-(4 * L1 - 1), 0, 4 * L3 - 1, -4 * L2, 4 * L2, 4 * (L1 - L3)] };
}
/** Elastic matrix (4 × 4, order xx, yy, tt, xy) of an isotropic material for the mode 'stress' | 'strain' | 'axi'. */
export function elasticMatrix(E, nu, mode) {
  if (mode === 'stress') { const c = E / (1 - nu * nu); return [[c, c * nu, 0, 0], [c * nu, c, 0, 0], [0, 0, 0, 0], [0, 0, 0, (c * (1 - nu)) / 2]]; }
  const lam = (E * nu) / ((1 + nu) * (1 - 2 * nu)), mu = E / (2 * (1 + nu)), a = lam + 2 * mu;
  return [[a, lam, lam, 0], [lam, a, lam, 0], [lam, lam, a, 0], [0, 0, 0, mu]];
}
/**
 * Structured triangular mesh of a mapped quadrilateral region. map(a, b) -> [x, y] for a, b in [0, 1]; na × nb cells, each split into two triangles.
 * order 1: 3-node triangles; order 2: 6-node triangles whose mid-side nodes lie on the mapped (curved) geometry.
 * Returns { X (Float64Array, 2 per node), elems: [[nodes…]], order, NA, NB (lattice size), id(i, j) (node number of lattice point i along a, j along b) }.
 */
export function mapMesh(na, nb, map, order = 2) {
  const s = order === 2 ? 2 : 1, NA = s * na + 1, NB = s * nb + 1, X = new Float64Array(2 * NA * NB), id = (i, j) => i * NB + j, elems = [];
  for (let i = 0; i < NA; i++) for (let j = 0; j < NB; j++) { const p = map(i / (NA - 1), j / (NB - 1)); X[2 * id(i, j)] = p[0]; X[2 * id(i, j) + 1] = p[1]; }
  for (let i = 0; i < na; i++) for (let j = 0; j < nb; j++) {
    const I = s * i, J = s * j, c = (di, dj) => id(I + di * s, J + dj * s), m = (di, dj) => id(I + di, J + dj); // corners and (order 2) mid-points on the fine lattice
    const tris = (i + j) % 2 === 0 ? [[[0, 0], [1, 0], [1, 1]], [[0, 0], [1, 1], [0, 1]]] : [[[0, 0], [1, 0], [0, 1]], [[1, 0], [1, 1], [0, 1]]];
    for (const t of tris) {
      const e = t.map((q) => c(q[0], q[1]));
      if (order === 2) for (let k = 0; k < 3; k++) { const a = t[k], b = t[(k + 1) % 3]; e.push(m(a[0] + b[0], a[1] + b[1])); }
      elems.push(e);
    }
  }
  return { X, elems, order: order === 2 ? 2 : 1, NA, NB, id, na, nb };
}
/**
 * Linear-elastic finite-element solution on a triangular mesh.
 * o: { mesh, E, nu, mode: 'stress' | 'strain' | 'axi', thickness (plane modes, default 1), fix: [{ node, dir: 0 | 1 }] (zero displacement),
 *      pressures: [{ nodes: [edge nodes in order: 2 for linear, 3 (end, middle, end) for quadratic], p (Pa, acting on the solid) , inside: [x, y] (any point of the solid next to the edge) }],
 *      forces: [{ node, dir, F }] }
 * Returns { u (Float64Array, 2 per node), stress (per node: [xx, yy, tt, xy]), vm (per node), gauss: [{ x, y, w (volume weight), s: [4] }] (only with o.keepGauss), hoopArea (axisymmetric: ∫σθ dA over the meridian section), energy (strain energy ½∫σ:ε),
 *   work (½ fᵀu), fext, fint (Float64Arrays), reactions: [{ node, dir, R }], residual (largest out-of-balance nodal force at a free degree of freedom / largest applied force), ndof, bandwidth }.
 */
export function feSolve(o) {
  const { mesh } = o, X = mesh.X, nn = mesh.order === 2 ? 6 : 3, nNode = X.length / 2, ndof = 2 * nNode, axi = o.mode === 'axi', th = axi ? 1 : nz(o.thickness, 1), D = elasticMatrix(+o.E, +o.nu, o.mode || 'strain'), rule = nn === 6 ? TRI7 : TRI3;
  need(o.E > 0 && o.nu > -1 && o.nu < 0.5, 'The finite-element model needs a positive modulus and a Poisson ratio below 0.5.');
  let bw = 0; for (const e of mesh.elems) { let lo = e[0], hi = e[0]; for (const n of e) { if (n < lo) lo = n; if (n > hi) hi = n; } bw = Math.max(bw, 2 * (hi - lo) + 1); }
  const w1 = bw + 1, A = new Float64Array(ndof * w1), f = new Float64Array(ndof), ed = new Array(2 * nn), B = [new Float64Array(2 * nn), new Float64Array(2 * nn), new Float64Array(2 * nn), new Float64Array(2 * nn)], DB = [new Float64Array(2 * nn), new Float64Array(2 * nn), new Float64Array(2 * nn), new Float64Array(2 * nn)];
  // strain–displacement matrix at one point of an element; returns the volume weight per unit quadrature weight
  const bAt = (e, L1, L2, L3) => {
    const sh = triShape(nn, L1, L2, L3); let x = 0, y = 0, xa = 0, xb = 0, ya = 0, yb = 0;
    for (let a = 0; a < nn; a++) { const px = X[2 * e[a]], py = X[2 * e[a] + 1]; x += sh.N[a] * px; y += sh.N[a] * py; xa += sh.dxi[a] * px; xb += sh.deta[a] * px; ya += sh.dxi[a] * py; yb += sh.deta[a] * py; }
    const det = xa * yb - xb * ya; need(Math.abs(det) > 1e-300, 'A finite element is degenerate: check the geometry.');
    for (let a = 0; a < nn; a++) { const dx = (yb * sh.dxi[a] - ya * sh.deta[a]) / det, dy = (-xb * sh.dxi[a] + xa * sh.deta[a]) / det; B[0][2 * a] = dx; B[0][2 * a + 1] = 0; B[1][2 * a] = 0; B[1][2 * a + 1] = dy; B[2][2 * a] = axi ? sh.N[a] / x : 0; B[2][2 * a + 1] = 0; B[3][2 * a] = dy; B[3][2 * a + 1] = dx; }
    return { x, y, vol: Math.abs(det) * (axi ? 2 * Math.PI * x : th) };
  };
  for (const e of mesh.elems) {
    for (let a = 0; a < nn; a++) { ed[2 * a] = 2 * e[a]; ed[2 * a + 1] = 2 * e[a] + 1; }
    for (const q of rule) {
      const g = bAt(e, q[0], q[1], q[2]), wq = g.vol * q[3];
      for (let r = 0; r < 4; r++) for (let c = 0; c < 2 * nn; c++) DB[r][c] = D[r][0] * B[0][c] + D[r][1] * B[1][c] + D[r][2] * B[2][c] + D[r][3] * B[3][c];
      for (let a = 0; a < 2 * nn; a++) { const ia = ed[a]; for (let b = 0; b < 2 * nn; b++) { const ib = ed[b]; if (ib > ia) continue; A[ia * w1 + ib - ia + bw] += wq * (B[0][a] * DB[0][b] + B[1][a] * DB[1][b] + B[2][a] * DB[2][b] + B[3][a] * DB[3][b]); } }
    }
  }
  // consistent nodal forces of the pressure loads (three-point Gauss along each edge, exact normal of the mapped edge)
  for (const pl of o.pressures || []) {
    const en = pl.nodes, k = en.length, p = +pl.p; if (!p) continue;
    const xs = en.map((n) => X[2 * n]), ys = en.map((n) => X[2 * n + 1]), xm = k === 3 ? xs[1] : 0.5 * (xs[0] + xs[1]), ym = k === 3 ? ys[1] : 0.5 * (ys[0] + ys[1]);
    for (const [s, wg] of GAUSS3) {
      const N = k === 3 ? [0.5 * s * (s - 1), 1 - s * s, 0.5 * s * (s + 1)] : [0.5 * (1 - s), 0.5 * (1 + s)], dN = k === 3 ? [s - 0.5, -2 * s, s + 0.5] : [-0.5, 0.5];
      let x = 0, tx = 0, ty = 0; for (let a = 0; a < k; a++) { x += N[a] * xs[a]; tx += dN[a] * xs[a]; ty += dN[a] * ys[a]; }
      let nx = ty, ny = -tx; if (nx * (xm - pl.inside[0]) + ny * (ym - pl.inside[1]) < 0) { nx = -nx; ny = -ny; } // outward normal × edge Jacobian
      const c = -p * wg * (axi ? 2 * Math.PI * x : th);
      for (let a = 0; a < k; a++) { f[2 * en[a]] += c * N[a] * nx; f[2 * en[a] + 1] += c * N[a] * ny; }
    }
  }
  for (const pf of o.forces || []) f[2 * pf.node + pf.dir] += pf.F;
  const fixed = new Uint8Array(ndof); for (const c of o.fix || []) fixed[2 * c.node + c.dir] = 1;
  const rhs = Float64Array.from(f);
  for (let d = 0; d < ndof; d++) if (fixed[d]) { // homogeneous constraints: clear the row and the column, unit diagonal
    for (let j = Math.max(0, d - bw); j < d; j++) A[d * w1 + j - d + bw] = 0;
    for (let i = d + 1; i <= Math.min(ndof - 1, d + bw); i++) A[i * w1 + d - i + bw] = 0;
    A[d * w1 + bw] = 1; rhs[d] = 0;
  }
  // banded Cholesky factorisation A = L·Lᵀ in place, then forward and back substitution
  for (let i = 0; i < ndof; i++) {
    const ri = i * w1 + bw - i, j0 = Math.max(0, i - bw);
    for (let j = j0; j <= i; j++) {
      const rj = j * w1 + bw - j; let s = A[ri + j]; const k0 = Math.max(j0, j - bw);
      for (let k = k0; k < j; k++) s -= A[ri + k] * A[rj + k];
      if (j === i) { need(s > 0, 'The finite-element model is not restrained (singular stiffness matrix): add supports.'); A[ri + i] = Math.sqrt(s); } else A[ri + j] = s / A[rj + j];
    }
  }
  const u = rhs;
  for (let i = 0; i < ndof; i++) { const ri = i * w1 + bw - i; let s = u[i]; for (let k = Math.max(0, i - bw); k < i; k++) s -= A[ri + k] * u[k]; u[i] = s / A[ri + i]; }
  for (let i = ndof - 1; i >= 0; i--) { let s = u[i]; for (let k = i + 1; k <= Math.min(ndof - 1, i + bw); k++) s -= A[k * w1 + bw - k + i] * u[k]; u[i] = s / A[i * w1 + bw]; }
  // stresses at the integration points (internal forces, energy) and at the nodes (averaged over the adjacent elements)
  const fint = new Float64Array(ndof), gauss = o.keepGauss ? [] : null, stress = Array.from({ length: nNode }, () => [0, 0, 0, 0]), cnt = new Uint16Array(nNode), sg = [0, 0, 0, 0], eg = [0, 0, 0, 0]; let energy = 0, hoopArea = 0;
  const sigAt = (e) => { for (let r = 0; r < 4; r++) { let s = 0; for (let a = 0; a < nn; a++) s += B[r][2 * a] * u[2 * e[a]] + B[r][2 * a + 1] * u[2 * e[a] + 1]; eg[r] = s; } for (let r = 0; r < 4; r++) sg[r] = D[r][0] * eg[0] + D[r][1] * eg[1] + D[r][2] * eg[2] + D[r][3] * eg[3]; if (o.mode === 'stress') sg[2] = 0; };
  const NODE_L = nn === 6 ? [[1, 0, 0], [0, 1, 0], [0, 0, 1], [0.5, 0.5, 0], [0, 0.5, 0.5], [0.5, 0, 0.5]] : [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
  for (const e of mesh.elems) {
    for (const q of rule) {
      const g = bAt(e, q[0], q[1], q[2]), wq = g.vol * q[3]; sigAt(e);
      energy += 0.5 * wq * (sg[0] * eg[0] + sg[1] * eg[1] + sg[2] * eg[2] + sg[3] * eg[3]); if (axi) hoopArea += (wq / (2 * Math.PI * g.x)) * sg[2]; if (gauss) gauss.push({ x: g.x, y: g.y, w: wq, s: sg.slice() });
      for (let a = 0; a < nn; a++) { fint[2 * e[a]] += wq * (B[0][2 * a] * sg[0] + B[2][2 * a] * sg[2] + B[3][2 * a] * sg[3]); fint[2 * e[a] + 1] += wq * (B[1][2 * a + 1] * sg[1] + B[3][2 * a + 1] * sg[3]); }
    }
    for (let a = 0; a < nn; a++) { const L = NODE_L[a]; bAt(e, L[0], L[1], L[2]); sigAt(e); const st = stress[e[a]]; for (let r = 0; r < 4; r++) st[r] += sg[r]; cnt[e[a]]++; }
  }
  for (let n = 0; n < nNode; n++) for (let r = 0; r < 4; r++) stress[n][r] /= cnt[n] || 1;
  let fmax = 0, res = 0, work = 0; const reactions = [];
  for (let d = 0; d < ndof; d++) { work += 0.5 * f[d] * u[d]; if (Math.abs(f[d]) > fmax) fmax = Math.abs(f[d]); if (fixed[d]) reactions.push({ node: d >> 1, dir: d & 1, R: fint[d] - f[d] }); else res = Math.max(res, Math.abs(fint[d] - f[d])); }
  return { u, stress, vm: stress.map((s) => vonMises(s[0], s[1], s[2], s[3])), gauss, hoopArea, energy, work, fext: f, fint, reactions, residual: fmax > 0 ? res / fmax : 0, ndof, bandwidth: bw };
}
/** Boundary edges of a mapped mesh along one lattice side ('a0' | 'a1' | 'b0' | 'b1') as node lists ready for a pressure load. */
function meshEdges(mesh, side) {
  const s = mesh.order, out = [], along = side[0] === 'a' ? mesh.nb : mesh.na, fix = side[1] === '0' ? 0 : side[0] === 'a' ? mesh.NA - 1 : mesh.NB - 1, at = (k) => (side[0] === 'a' ? mesh.id(fix, k) : mesh.id(k, fix));
  for (let c = 0; c < along; c++) out.push(s === 2 ? [at(2 * c), at(2 * c + 1), at(2 * c + 2)] : [at(c), at(c + 1)]);
  return out;
}
const smooth01 = (x) => { const c = clamp(x, 0, 1); return c * c * (3 - 2 * c); };
/**
 * Thick-walled cylinder under internal and external pressure, solved with axisymmetric solid elements (r–z slice, plane strain axially) —
 * the numerical counterpart of the Lamé solution. o: { ri, ro, pi, pe, E, nu, nr, nz, order }
 * Returns { uBore, uOuter (radial displacement, m), hoopBore, hoopOuter, radialBore (Pa), energy (J per metre of pipe), hoopForce (∫σθ dr per metre, N/m), fe }.
 */
export function feCylinder(o) {
  const ri = +o.ri, ro = +o.ro, nr = Math.max(1, Math.round(nz(o.nr, 4))), nzc = Math.max(1, Math.round(nz(o.nz, 1))), H = ((ro - ri) * nzc) / nr, mesh = mapMesh(nzc, nr, (a, b) => [ri + (ro - ri) * b, H * a], o.order === 1 ? 1 : 2), fix = [];
  for (let j = 0; j < mesh.NB; j++) { fix.push({ node: mesh.id(0, j), dir: 1 }); fix.push({ node: mesh.id(mesh.NA - 1, j), dir: 1 }); }
  const mid = [(ri + ro) / 2, H / 2], fe = feSolve({ mesh, E: o.E, nu: o.nu, mode: 'axi', fix, pressures: [...meshEdges(mesh, 'b0').map((nodes) => ({ nodes, p: nz(o.pi, 0), inside: mid })), ...meshEdges(mesh, 'b1').map((nodes) => ({ nodes, p: nz(o.pe, 0), inside: mid }))] });
  const iM = Math.floor((mesh.NA - 1) / 2), nb = mesh.id(iM, 0), no = mesh.id(iM, mesh.NB - 1);
  return { uBore: fe.u[2 * nb], uOuter: fe.u[2 * no], hoopBore: fe.stress[nb][2], hoopOuter: fe.stress[no][2], radialBore: fe.stress[nb][0], energy: fe.energy / H, hoopForce: fe.hoopArea / H, fe, mesh };
}
/**
 * Pipe wall with a thinned band (metal loss all around the circumference over an axial length), axisymmetric solid elements on half of the band.
 * o: { ro (outer radius), t, d (depth of the loss), L (axial length of the loss), side: 'in' | 'out', pi, pe, E, nu, nz (cells along the pipe), nr (cells through the wall), order, flow (flow stress, Pa) }
 * Returns { vmMax, vmNominal (Lamé, intact wall at the bore), scf, hoopLigament (mean hoop stress of the remaining ligament at the centre of the band), hoopNominal (mean hoop stress of the intact wall),
 *   ligamentFactor (hoopLigament / hoopNominal), collapse (differential pressure at which the mean ligament hoop stress reaches the flow stress, Pa; null without flow), radialLoad, hoopResultant (N: the two sides of the
 *   radial equilibrium ∫σθ dA = Σ radial nodal forces / 2π), z[], vmSurface[] (von Mises along the thinned surface), fe, mesh }.
 */
export function feGroove(o) {
  const ro = +o.ro, t = +o.t, ri = ro - t, d = clamp(nz(o.d, 0), 0, 0.95 * t), L = Math.max(nz(o.L, 0), 1e-6), inner = o.side !== 'out', ramp = Math.min(Math.max(d, 0.15 * t), 0.45 * L), Rm = ro - t / 2, decay = Math.sqrt(Rm * t);
  need(ro > 0 && t > 0 && t < ro, 'The thinned-wall model needs a wall thickness smaller than the outer radius.');
  const z1 = L / 2 + ramp / 2 + 0.6 * decay, Lm = z1 + 4.5 * decay, nzc = Math.max(4, Math.round(nz(o.nz, 28))), nr = Math.max(1, Math.round(nz(o.nr, 3))), loss = (z) => d * smooth01((L / 2 + ramp / 2 - z) / ramp);
  const mesh = mapMesh(nzc, nr, (a, b) => { const z = a < 0.6 ? (a / 0.6) * z1 : z1 + ((a - 0.6) / 0.4) * (Lm - z1), dl = loss(z), r0 = inner ? ri + dl : ri, r1 = inner ? ro : ro - dl; return [r0 + (r1 - r0) * b, z]; }, o.order === 1 ? 1 : 2), fix = [];
  for (let j = 0; j < mesh.NB; j++) { fix.push({ node: mesh.id(0, j), dir: 1 }); fix.push({ node: mesh.id(mesh.NA - 1, j), dir: 1 }); }
  const pi = nz(o.pi, 0), pe = nz(o.pe, 0), ins = (nodes) => { const n = nodes[nodes.length === 3 ? 1 : 0], a = n - (n % mesh.NB), b = a + mesh.NB - 1; return [0.5 * (mesh.X[2 * a] + mesh.X[2 * b]), mesh.X[2 * n + 1]]; };
  const fe = feSolve({ mesh, E: o.E, nu: o.nu, mode: 'axi', fix, pressures: [...meshEdges(mesh, 'b0').map((nodes) => ({ nodes, p: pi, inside: ins(nodes) })), ...meshEdges(mesh, 'b1').map((nodes) => ({ nodes, p: pe, inside: ins(nodes) }))] });
  let lig = 0; const tl = t - d; // mean hoop stress of the ligament on the symmetry plane from the nodal stresses (Simpson / trapezoid through the wall)
  { const J = mesh.NB, s = (j) => fe.stress[mesh.id(0, j)][2], r = (j) => mesh.X[2 * mesh.id(0, j)]; for (let j = 0; j + 1 < J; j++) lig += 0.5 * (s(j) + s(j + 1)) * (r(j + 1) - r(j)); lig /= tl; }
  const la = lame(pi, pe, ri, ro, ri), vmNom = vonMises(la.hoop, o.nu * (la.hoop + la.radial), la.radial), hoopNom = (pi * ri - pe * ro) / t, jS = inner ? 0 : mesh.NB - 1, zs = [], vs = []; let vmMax = 0, zMax = 0;
  for (let i = 0; i < mesh.NA; i++) { const n = mesh.id(i, jS); zs.push(mesh.X[2 * n + 1]); vs.push(fe.vm[n]); }
  fe.vm.forEach((v, n) => { if (v > vmMax) { vmMax = v; zMax = mesh.X[2 * n + 1]; } });
  let fr = 0; for (let n = 0; n < fe.fext.length / 2; n++) fr += fe.fext[2 * n]; const hr = 2 * Math.PI * fe.hoopArea;
  return { vmMax, zMax, vmNominal: vmNom, scf: vmNom > 0 ? vmMax / vmNom : 1, hoopLigament: lig, hoopNominal: hoopNom, ligamentFactor: hoopNom !== 0 ? lig / hoopNom : 1, collapse: o.flow > 0 && lig > 0 ? ((pi - pe) * o.flow) / lig : null, radialLoad: fr, hoopResultant: hr, z: zs, vmSurface: vs, length: Lm, fe, mesh };
}
/**
 * Cross-section of a pipe as a plane-strain ring (half ring, symmetric about the x axis) with out-of-roundness, a dent and/or a locally thinned wall at θ = 0.
 * o: { ro (nominal outer radius), t, pi, pe, E, nu, nth (cells around the half ring), nr (cells through the wall), order,
 *      ovality (f0 = (Dmax − Dmin)/D; the long axis is the x axis), dent: { depth (m, inward), halfAngle (rad) }, thin: { depth (m), halfAngle (rad), side: 'in' | 'out' } }
 * Returns { hoopMax, hoopMin (Pa, with the angle where they occur), hoopNominal ((pi·ri − pe·ro)/t), scf (largest |hoop| / |nominal|), vmMax, theta[], hoopInner[], hoopOuter[],
 *   crown: { inner, outer } (hoop stress at θ = 0), dDiameter: { x, y } (change of the two diameters, m), balance: { applied, reaction } (N per metre in y on the half ring), fe, mesh }.
 */
export function feRing(o) {
  const ro = +o.ro, t = +o.t, Rm = ro - t / 2, w1 = (nz(o.ovality, 0) * 2 * Rm) / 4, dent = o.dent && o.dent.depth > 0 ? o.dent : null, thin = o.thin && o.thin.depth > 0 ? o.thin : null, local = dent || thin;
  need(ro > 0 && t > 0 && t < ro, 'The ring model needs a wall thickness smaller than the outer radius.');
  const nth = Math.max(4, Math.round(nz(o.nth, 40))), nr = Math.max(1, Math.round(nz(o.nr, 2))), bump = (th, half) => { const x = Math.abs(th) / Math.max(half, 1e-6); return x >= 1 ? 0 : Math.cos((Math.PI * x) / 2) ** 2; }, gam = local ? 1.7 : 1, thOf = (a) => { const u = clamp(a, 0, 1) * nth, k = Math.min(Math.floor(u + 1e-9), nth - 1), t0 = Math.PI * (k / nth) ** gam, t1 = Math.PI * ((k + 1) / nth) ** gam; return t0 + (u - k) * (t1 - t0); }; // corner nodes graded towards θ = 0, mid-side nodes halfway
  const mesh = mapMesh(nth, nr, (a, b) => {
    const th = thOf(a), rm = Rm + w1 * Math.cos(2 * th) - (dent ? dent.depth * bump(th, dent.halfAngle) : 0), dl = thin ? Math.min(thin.depth, 0.95 * t) * bump(th, thin.halfAngle) : 0, r0 = rm - t / 2 + (thin && thin.side !== 'out' ? dl : 0), r1 = rm + t / 2 - (thin && thin.side === 'out' ? dl : 0), r = r0 + (r1 - r0) * b;
    return [r * Math.cos(th), r * Math.sin(th)];
  }, o.order === 1 ? 1 : 2), fix = [];
  for (let j = 0; j < mesh.NB; j++) { fix.push({ node: mesh.id(0, j), dir: 1 }); fix.push({ node: mesh.id(mesh.NA - 1, j), dir: 1 }); }
  fix.push({ node: mesh.id(mesh.NA - 1, 0), dir: 0 }); // removes the rigid translation along the axis of symmetry; its reaction must vanish
  const pi = nz(o.pi, 0), pe = nz(o.pe, 0), ins = (nodes) => { const n = nodes[nodes.length === 3 ? 1 : 0], a = n - (n % mesh.NB), b = a + mesh.NB - 1; return [0.5 * (mesh.X[2 * a] + mesh.X[2 * b]), 0.5 * (mesh.X[2 * a + 1] + mesh.X[2 * b + 1])]; };
  const fe = feSolve({ mesh, E: o.E, nu: o.nu, mode: 'strain', fix, pressures: [...meshEdges(mesh, 'b0').map((nodes) => ({ nodes, p: pi, inside: ins(nodes) })), ...meshEdges(mesh, 'b1').map((nodes) => ({ nodes, p: pe, inside: ins(nodes) }))] });
  const hoopAt = (n) => { const x = mesh.X[2 * n], y = mesh.X[2 * n + 1], r = Math.hypot(x, y) || 1, c = x / r, s = y / r, st = fe.stress[n]; return st[0] * s * s + st[1] * c * c - 2 * st[3] * s * c; };
  const theta = [], hi = [], ho = []; let hMax = -Infinity, hMin = Infinity, thMax = 0, thMin = 0;
  for (let i = 0; i < mesh.NA; i++) { const th = thOf(i / (mesh.NA - 1)), a = hoopAt(mesh.id(i, 0)), b = hoopAt(mesh.id(i, mesh.NB - 1)); theta.push((th * 180) / Math.PI); hi.push(a); ho.push(b); for (const v of [a, b]) { if (v > hMax) { hMax = v; thMax = th; } if (v < hMin) { hMin = v; thMin = th; } } }
  const ri = ro - t, nom = (pi * ri - pe * ro) / t, iq = Array.from({ length: mesh.NA }, (_, i) => i).reduce((b2, i) => (Math.abs(thOf(i / (mesh.NA - 1)) - Math.PI / 2) < Math.abs(thOf(b2 / (mesh.NA - 1)) - Math.PI / 2) ? i : b2), 0); let fy = 0, ry = 0, rx = 0; for (let n = 0; n < fe.fext.length / 2; n++) fy += fe.fext[2 * n + 1]; for (const r of fe.reactions) { if (r.dir === 1) ry += r.R; else rx += r.R; }
  return { hoopMax: hMax, hoopMin: hMin, thetaMax: (thMax * 180) / Math.PI, thetaMin: (thMin * 180) / Math.PI, hoopNominal: nom, scf: nom !== 0 ? Math.max(Math.abs(hMax), Math.abs(hMin)) / Math.abs(nom) : 1, vmMax: Math.max(...fe.vm), theta, hoopInner: hi, hoopOuter: ho, crown: { inner: hi[0], outer: ho[0] },
    dDiameter: { x: fe.u[2 * mesh.id(0, mesh.NB - 1)] - fe.u[2 * mesh.id(mesh.NA - 1, mesh.NB - 1)], y: 2 * fe.u[2 * mesh.id(iq, mesh.NB - 1) + 1] }, balance: { applied: fy, reaction: -ry, axisReaction: rx }, fe, mesh };
}

// ---- axisymmetric thin shell: circular cylinder with wall-thickness steps and ring stiffeners ---------------------------
/**
 * Axisymmetric bending of a thin cylindrical shell (shell of revolution with a straight meridian): d²/dx²(D·w″) + (E·t/R²)·w = p − ν·N_x/R,
 * D = E·t³/12(1 − ν²), solved with Hermite elements along the meridian. Radial deflection w is positive outwards.
 * o: { R (mid-surface radius), L, t (number, or function of x), E, nu, p (net pressure, positive outwards; number or function of x), n (elements),
 *      axial: 'capped' (N_x = p·R/2) | 'open' (N_x = 0) | 'restrained' (ε_x = 0), ends: ['clamped' | 'simple' | 'free' | 'symmetry', …],
 *      rings: [{ x, A (cross-section area of the ring, m²), R (its centroid radius, default the shell radius) }], breaks: [x] (extra nodes, e.g. thickness steps) }
 * Returns { x[], w[], slope[], moment[] (N·m per m of circumference), hoopN[] (N/m), hoopMembrane[], bendingStress[] (meridional, at the surface), hoopSurface[] (largest |hoop| at a surface, Pa),
 *   wMembrane (far-field deflection of the first element's thickness), beta (1/m), radialLoad, radialResistance (N per m of circumference: the two sides of the radial equilibrium) }.
 */
export function shellCylinder(o) {
  const R = +o.R, L = +o.L, E = +o.E, nu = +o.nu, n = clamp(Math.round(nz(o.n, 80)), 2, 400), tOf = typeof o.t === 'function' ? o.t : () => +o.t, pOf = typeof o.p === 'function' ? o.p : () => nz(o.p, 0), axial = o.axial || 'capped', ends = o.ends || ['symmetry', 'symmetry'];
  need(R > 0 && L > 0 && E > 0, 'The shell model needs a positive radius, length and modulus.');
  const bp = [0, L, ...(o.breaks || []), ...(o.rings || []).map((r) => r.x)].filter((x) => x >= 0 && x <= L).sort((a, b) => a - b).filter((x, i, a) => i === 0 || x - a[i - 1] > 1e-9 * L), x = [0];
  for (let s = 1; s < bp.length; s++) { const a = bp[s - 1], b = bp[s], ne = Math.max(1, Math.round((n * (b - a)) / L)); for (let i = 1; i <= ne; i++) x.push(a + ((b - a) * i) / ne); }
  const nd = x.length, ndof = 2 * nd, K = zeros(ndof), f = new Float64Array(ndof), el = []; let load = 0;
  for (let e = 0; e < nd - 1; e++) {
    const l = x[e + 1] - x[e], xm = 0.5 * (x[e] + x[e + 1]), t = tOf(xm), D = (E * t ** 3) / (12 * (1 - nu * nu)), kf = axial === 'restrained' ? (E * t) / ((1 - nu * nu) * R * R) : (E * t) / (R * R), p = pOf(xm), q = axial === 'capped' ? p * (1 - nu / 2) : p, c = D / l ** 3, m = (kf * l) / 420;
    const ke = [[12, 6 * l, -12, 6 * l], [6 * l, 4 * l * l, -6 * l, 2 * l * l], [-12, -6 * l, 12, -6 * l], [6 * l, 2 * l * l, -6 * l, 4 * l * l]], me = [[156, 22 * l, 54, -13 * l], [22 * l, 4 * l * l, 13 * l, -3 * l * l], [54, 13 * l, 156, -22 * l], [-13 * l, -3 * l * l, -22 * l, 4 * l * l]], fe = [(q * l) / 2, (q * l * l) / 12, (q * l) / 2, (-q * l * l) / 12];
    for (let a = 0; a < 4; a++) { f[2 * e + a] += fe[a]; for (let b = 0; b < 4; b++) K[2 * e + a][2 * e + b] += c * ke[a][b] + m * me[a][b]; }
    el.push({ l, t, D, kf, q, p, me: me.map((r) => r.map((v) => m * v)), ke: ke.map((r) => r.map((v) => c * v)), fe }); load += q * l;
  }
  const springs = [];
  for (const r of o.rings || []) { let node = 0; for (let i = 1; i < nd; i++) if (Math.abs(x[i] - r.x) < Math.abs(x[node] - r.x)) node = i; const Rr = pos(r.R, R), k = (E * +r.A) / (Rr * Rr); K[2 * node][2 * node] += k; springs.push({ node, k }); }
  const fixed = new Set(); [0, nd - 1].forEach((node, i) => { const kind = ends[i] || 'symmetry'; if (kind === 'clamped') { fixed.add(2 * node); fixed.add(2 * node + 1); } else if (kind === 'simple') fixed.add(2 * node); else if (kind === 'symmetry') fixed.add(2 * node + 1); });
  const fr = []; for (let i = 0; i < ndof; i++) if (!fixed.has(i)) fr.push(i);
  const Lc = chol(sub(K, fr)); need(Lc, 'The shell model is singular.');
  const ur = cholSolve(Lc, fr.map((d) => f[d])), u = new Float64Array(ndof); fr.forEach((d, i) => (u[d] = ur[i]));
  const mom = new Float64Array(nd), cnt = new Float64Array(nd), tn = new Float64Array(nd); let resist = 0;
  for (let e = 0; e < nd - 1; e++) {
    const { l, D, t, me } = el[e], w1 = u[2 * e], t1 = u[2 * e + 1], w2 = u[2 * e + 2], t2 = u[2 * e + 3];
    mom[e] += (D * (-6 * w1 - 4 * l * t1 + 6 * w2 - 2 * l * t2)) / (l * l); cnt[e]++; mom[e + 1] += (D * (6 * w1 + 2 * l * t1 - 6 * w2 + 4 * l * t2)) / (l * l); cnt[e + 1]++; tn[e] += t; tn[e + 1] += t;
    for (const a of [0, 2]) for (let b = 0; b < 4; b++) resist += me[a][b] * u[2 * e + b]; // foundation (hoop) reaction, weighted with the constant test function
  }
  for (const s of springs) resist += s.k * u[2 * s.node];
  const Ku = matVec(K, u); for (const d of fixed) if (d % 2 === 0) resist += f[d] - Ku[d]; // end shear carried by a support
  const w = [], slope = [], moment = [], hoopN = [], hoopM = [], sb = [], hs = [];
  for (let i = 0; i < nd; i++) {
    const t = tn[i] / cnt[i], M = mom[i] / cnt[i], p = pOf(x[i]), Nx = axial === 'capped' ? (p * R) / 2 : 0, Nt = axial === 'restrained' ? ((E * t) / (1 - nu * nu)) * (u[2 * i] / R) : (E * t * u[2 * i]) / R + nu * Nx;
    w.push(u[2 * i]); slope.push(u[2 * i + 1]); moment.push(M); hoopN.push(Nt); hoopM.push(Nt / t); sb.push((6 * M) / (t * t)); hs.push(Math.abs(Nt / t) + Math.abs((nu * 6 * M) / (t * t)));
  }
  const t0 = el[0].t, k0 = el[0].kf;
  return { x, w, slope, moment, hoopN, hoopMembrane: hoopM, bendingStress: sb, hoopSurface: hs, wMembrane: el[0].q / k0, beta: ((3 * (1 - nu * nu)) / (R * R * t0 * t0)) ** 0.25, radialLoad: load, radialResistance: resist };
}

// ---- two-way fluid–structure interaction: a beam conveying slug flow -----------------------------------------------------
// Pipe conveying fluid (Païdoussis): EI·w⁗ + M_f·U²·w″ + 2·M_f·U·ẇ′ + (M_f + m)·ẅ = q. The conveyed mass per length M_f(x, t) follows the slug
// train of the one-dimensional flow model, so mass, Coriolis (gyroscopic) and centrifugal (curvature) terms all change as a slug crosses the span,
// and the pipe motion acts back on the fluid force. The equations are projected on the span modes and integrated together with them.
const GAUSS4 = [[-0.8611363115940526, 0.3478548451374538], [-0.3399810435848563, 0.6521451548625461], [0.3399810435848563, 0.6521451548625461], [0.8611363115940526, 0.3478548451374538]];
/**
 * Running integrals of the modal products along the span: C = ∫φᵢφⱼ, D = ∫φᵢφⱼ′, H = ∫φᵢφⱼ″ and W = ∫φᵢ from 0 to every node (Hermite interpolation,
 * four-point Gauss). Returns { p, nd, x, C, D, H (Float64Array nd·p·p), W (Float64Array nd·p), full: { C, D, H } (p × p arrays over the whole span) }.
 */
export function modalIntegrals(model, modes) {
  const p = modes.f.length, x = model.x, nd = x.length, pp = p * p, C = new Float64Array(nd * pp), D = new Float64Array(nd * pp), H = new Float64Array(nd * pp), W = new Float64Array(nd * p), f0 = new Float64Array(p), f1 = new Float64Array(p), f2 = new Float64Array(p);
  for (let e = 0; e < nd - 1; e++) {
    const l = x[e + 1] - x[e], o0 = e * pp, o1 = (e + 1) * pp;
    for (let k = 0; k < pp; k++) { C[o1 + k] = C[o0 + k]; D[o1 + k] = D[o0 + k]; H[o1 + k] = H[o0 + k]; }
    for (let k = 0; k < p; k++) W[(e + 1) * p + k] = W[e * p + k];
    for (const [g, wg] of GAUSS4) {
      const s = 0.5 * (1 + g), w = 0.5 * l * wg, N = [1 - 3 * s * s + 2 * s ** 3, l * (s - 2 * s * s + s ** 3), 3 * s * s - 2 * s ** 3, l * (-s * s + s ** 3)], dN = [(-6 * s + 6 * s * s) / l, 1 - 4 * s + 3 * s * s, (6 * s - 6 * s * s) / l, -2 * s + 3 * s * s], ddN = [(-6 + 12 * s) / (l * l), (-4 + 6 * s) / l, (6 - 12 * s) / (l * l), (-2 + 6 * s) / l];
      for (let k = 0; k < p; k++) { const u = modes.vectors[k]; let a = 0, b = 0, c = 0; for (let q = 0; q < 4; q++) { const v = u[2 * e + q]; a += N[q] * v; b += dN[q] * v; c += ddN[q] * v; } f0[k] = a; f1[k] = b; f2[k] = c; W[(e + 1) * p + k] += w * a; }
      for (let i = 0; i < p; i++) for (let j = 0; j < p; j++) { C[o1 + i * p + j] += w * f0[i] * f0[j]; D[o1 + i * p + j] += w * f0[i] * f1[j]; H[o1 + i * p + j] += w * f0[i] * f2[j]; }
    }
  }
  const last = (A) => Array.from({ length: p }, (_, i) => Array.from({ length: p }, (_, j) => A[(nd - 1) * pp + i * p + j]));
  return { p, nd, x, C, D, H, W, full: { C: last(C), D: last(D), H: last(H) } };
}
/**
 * Critical (divergence) velocity of a span conveying fluid of mass per length mf: the lowest U at which K − mf·U²·K_c becomes singular.
 * Analytic value for a pinned–pinned pipe: (π/L)·√(EI/mf). Returns { vc (m/s; null when the supports make the system non-conservative), lambda[] }.
 */
export function fluidCritical(model, modes, mf, mi = null) {
  const I = mi || modalIntegrals(model, modes), p = I.p, Kc = zeros(p), K = zeros(p);
  for (let i = 0; i < p; i++) { K[i][i] = modes.omega[i] ** 2; for (let j = 0; j < p; j++) Kc[i][j] = -0.5 * mf * (I.full.H[i][j] + I.full.H[j][i]); }
  if (!(mf > 0) || !chol(Kc)) return { vc: null, lambda: [] };
  const e = genEig(K, Kc); return { vc: Math.sqrt(Math.max(e.values[0], 0)), lambda: e.values };
}
/**
 * First natural frequency (Hz) of the span with fluid of mass per length mf flowing at velocity v, from the gyroscopic eigenproblem
 * (K − ω²·I + iω·G) q = 0 written as a real symmetric problem of twice the size. Returns 0 at or above the critical velocity.
 */
export function fluidFrequency(model, modes, mf, v, mi = null) {
  const I = mi || modalIntegrals(model, modes), p = Math.min(I.p, 4), K = zeros(p), Gm = zeros(p);
  for (let i = 0; i < p; i++) for (let j = 0; j < p; j++) { K[i][j] = (i === j ? modes.omega[i] ** 2 : 0) + 0.5 * mf * v * v * (I.full.H[i][j] + I.full.H[j][i]); Gm[i][j] = mf * v * (I.full.D[i][j] - I.full.D[j][i]); }
  const lmin = (om) => { const A = zeros(2 * p); for (let i = 0; i < p; i++) for (let j = 0; j < p; j++) { const k = K[i][j] - (i === j ? om * om : 0); A[i][j] = k; A[p + i][p + j] = k; A[i][p + j] = -om * Gm[i][j]; A[p + i][j] = om * Gm[i][j]; } return jacobiEig(A).values[0]; };
  if (!(lmin(0) > 0)) return 0;
  const top = p > 1 ? Math.min(1.05 * modes.omega[0], 0.5 * (modes.omega[0] + modes.omega[1])) : 1.05 * modes.omega[0]; // between the first and the second frequency exactly one eigenvalue is negative
  return lmin(top) < 0 ? brent(lmin, 0, top, 1e-10 * top) / (2 * Math.PI) : modes.f[0];
}
/**
 * Response of a span to a train of slugs with two-way coupling, by modal projection and Newmark-β (average acceleration) with time-dependent matrices.
 * o: { model, modes, zeta, slugs: [{ t0, len, v }], v (convective velocity of the fluid, m/s), mf (conveyed mass per length between slugs, already part of model.m),
 *      dM (extra conveyed mass per length inside a slug body), g (gravity, default 9.80665; the slug weight is dM·g), bend: { x, dF } | null (extra momentum force while the body is in the bend;
 *      it follows the relative velocity of the moving bend), dt, tEnd, ro, E, ramp, coupled (false = the fluid only loads the pipe), q0 (initial modal amplitudes), mi (modalIntegrals) }
 * Returns { t[], sigma[] (Pa at the critical node), disp[] (m), vel[] (m/s at that node), node, x, sigmaMax, dispMax, velMax, velRms (m/s over the history at the node of largest motion), steps, massRatio (largest added modal mass) }.
 */
export function fsiResponse(o) {
  const { model, modes } = o, p = modes.f.length, nd = model.x.length, x = model.x, L = model.L, dt = o.dt, steps = Math.max(2, Math.round(o.tEnd / dt)), zeta = nz(o.zeta, 0.02), I = o.mi || modalIntegrals(model, modes), pp = p * p, coupled = o.coupled !== false;
  const v0 = nz(o.v, 0), mf = coupled ? nz(o.mf, 0) : 0, dM = nz(o.dM, 0), grav = nz(o.g, G), slugs = o.slugs || [], bend = o.bend && o.bend.dF ? o.bend : null, ramp = Math.max(nz(o.ramp, 0.3), 1e-3);
  const disp = modes.vectors.map((u) => Float64Array.from({ length: nd }, (_, i) => u[2 * i]));
  let lf = 0; // fraction inside the element returned by locate()
  const locate = (s) => { if (s <= 0) { lf = 0; return 0; } if (s >= L) { lf = 1; return nd - 2; } let lo = 0, hi = nd - 1; while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (x[mid] <= s) lo = mid; else hi = mid; } lf = (s - x[lo]) / (x[hi] - x[lo]); return lo; };
  const phiB = bend ? (() => { const e = locate(bend.x), f = lf; return disp.map((d) => d[e] + f * (d[e + 1] - d[e])); })() : null;
  // constant parts: unit modal mass, modal damping, modal stiffness, steady conveyed fluid (gyroscopic and centrifugal terms)
  const M0 = new Float64Array(pp), C0 = new Float64Array(pp), K0 = new Float64Array(pp);
  for (let i = 0; i < p; i++) { M0[i * p + i] = 1; C0[i * p + i] = 2 * zeta * modes.omega[i]; K0[i * p + i] = modes.omega[i] ** 2; for (let j = 0; j < p; j++) { C0[i * p + j] += 2 * mf * v0 * I.full.D[i][j]; K0[i * p + j] += mf * v0 * v0 * I.full.H[i][j]; } }
  const Mt = new Float64Array(pp), Ct = new Float64Array(pp), Kt = new Float64Array(pp), S = new Float64Array(pp), Q = new Float64Array(p), rhs = new Float64Array(p), IC = I.C, ID = I.D, IH = I.H, IW = I.W;
  let first = 0, massRatio = 0, dirty = true; const wSlug = o.dw !== undefined ? +o.dw : dM * grav;
  const assemble = (t) => {
    if (dirty) { Mt.set(M0); Ct.set(C0); Kt.set(K0); dirty = false; }
    Q.fill(0);
    while (first < slugs.length && (t - slugs[first].t0) * slugs[first].v - slugs[first].len > L + ramp) first++;
    for (let k = first; k < slugs.length; k++) {
      const s = slugs[k]; if (t < s.t0) break;
      const front = (t - s.t0) * s.v, tail = front - s.len, a = Math.max(tail, 0), b = Math.min(front, L);
      if (b > a && (dM || wSlug)) {
        const ea = locate(a), fa = lf, eb = locate(b), fb = lf, wa = ea * p, wb = eb * p;
        for (let i = 0; i < p; i++) Q[i] += wSlug * (IW[wb + i] + fb * (IW[wb + p + i] - IW[wb + i]) - IW[wa + i] - fa * (IW[wa + p + i] - IW[wa + i]));
        if (coupled) {
          const oa = ea * pp, ob = eb * pp, cG = 2 * dM * s.v, cK = dM * s.v * s.v; dirty = true;
          for (let q = 0; q < pp; q++) {
            Mt[q] += dM * (IC[ob + q] + fb * (IC[ob + pp + q] - IC[ob + q]) - IC[oa + q] - fa * (IC[oa + pp + q] - IC[oa + q]));
            Ct[q] += cG * (ID[ob + q] + fb * (ID[ob + pp + q] - ID[ob + q]) - ID[oa + q] - fa * (ID[oa + pp + q] - ID[oa + q]));
            Kt[q] += cK * (IH[ob + q] + fb * (IH[ob + pp + q] - IH[ob + q]) - IH[oa + q] - fa * (IH[oa + pp + q] - IH[oa + q]));
          }
        }
      }
      if (bend) { const c = clamp((front - bend.x) / ramp, 0, 1) - clamp((tail - bend.x) / ramp, 0, 1); if (c > 0) { for (let i = 0; i < p; i++) Q[i] += bend.dF * c * phiB[i]; if (coupled) { dirty = true; const cc = (2 * bend.dF * c) / s.v; for (let i = 0; i < p; i++) for (let j = 0; j < p; j++) Ct[i * p + j] += cc * phiB[i] * phiB[j]; } } }
    }
    if (coupled && dirty) for (let i = 0; i < p; i++) massRatio = Math.max(massRatio, Mt[i * p + i] - 1);
  };
  const solve = (A, b) => { // Gaussian elimination with partial pivoting, in place (p is small)
    for (let k = 0; k < p; k++) {
      let piv = k; for (let i = k + 1; i < p; i++) if (Math.abs(A[i * p + k]) > Math.abs(A[piv * p + k])) piv = i;
      if (piv !== k) { for (let j = 0; j < p; j++) { const tmp = A[k * p + j]; A[k * p + j] = A[piv * p + j]; A[piv * p + j] = tmp; } const tb = b[k]; b[k] = b[piv]; b[piv] = tb; }
      const d = A[k * p + k]; need(Math.abs(d) > 1e-300, 'The coupled span model is singular: the flow velocity is at the critical velocity of the span.');
      for (let i = k + 1; i < p; i++) { const f = A[i * p + k] / d; if (f !== 0) { for (let j = k; j < p; j++) A[i * p + j] -= f * A[k * p + j]; b[i] -= f * b[k]; } }
    }
    for (let i = p - 1; i >= 0; i--) { let s = b[i]; for (let j = i + 1; j < p; j++) s -= A[i * p + j] * b[j]; b[i] = s / A[i * p + i]; }
  };
  const q = new Float64Array(p * (steps + 1)), qd = new Float64Array(p * (steps + 1)), qa = new Float64Array(p * (steps + 1)), u = new Float64Array(p), vv = new Float64Array(p), acc = new Float64Array(p), up = new Float64Array(p), vp = new Float64Array(p);
  if (o.q0) for (let i = 0; i < p; i++) { u[i] = nz(o.q0[i], 0); q[i] = u[i]; }
  assemble(0); S.set(Mt); for (let i = 0; i < p; i++) { let s = Q[i]; for (let j = 0; j < p; j++) s -= Ct[i * p + j] * vv[j] + Kt[i * p + j] * u[j]; rhs[i] = s; } solve(S, rhs); acc.set(rhs);
  const h2 = dt / 2, h4 = (dt * dt) / 4;
  for (let n = 1; n <= steps; n++) {
    assemble(n * dt);
    for (let i = 0; i < p; i++) { up[i] = u[i] + dt * vv[i] + h4 * acc[i]; vp[i] = vv[i] + h2 * acc[i]; }
    for (let i = 0; i < p; i++) { let s = Q[i]; for (let j = 0; j < p; j++) { const k = i * p + j; s -= Ct[k] * vp[j] + Kt[k] * up[j]; S[k] = Mt[k] + h2 * Ct[k] + h4 * Kt[k]; } rhs[i] = s; }
    solve(S, rhs);
    for (let i = 0; i < p; i++) { u[i] = up[i] + h4 * rhs[i]; vv[i] = vp[i] + h2 * rhs[i]; acc[i] = rhs[i]; q[n * p + i] = u[i]; qd[n * p + i] = vv[i]; qa[n * p + i] = rhs[i]; }
  }
  // envelope of the bending curvature at every node, then the histories at the worst node
  let node = 0, worst = -1, nodeD = 0, dWorst = -1; const stride = Math.max(1, Math.floor(steps / 4000));
  for (let j = 0; j < nd; j++) { let lo = 0, hi = 0, dl = 0, dh = 0; for (let i = 0; i <= steps; i += stride) { let c = 0, d = 0; for (let k = 0; k < p; k++) { c += q[i * p + k] * modes.curv[k][j]; d += q[i * p + k] * disp[k][j]; } if (c < lo) lo = c; if (c > hi) hi = c; if (d < dl) dl = d; if (d > dh) dh = d; } if (hi - lo > worst) { worst = hi - lo; node = j; } if (dh - dl > dWorst) { dWorst = dh - dl; nodeD = j; } }
  const fac = +o.E * +o.ro, t = new Array(steps + 1), sg = new Array(steps + 1), ds = new Array(steps + 1), vs = new Array(steps + 1); let sMax = 0, dMax = 0, vMax = 0, v2 = 0, a2 = 0;
  for (let i = 0; i <= steps; i++) { let c = 0, d = 0, w = 0, dd = 0, aa = 0; for (let k = 0; k < p; k++) { c += q[i * p + k] * modes.curv[k][node]; d += q[i * p + k] * disp[k][node]; w += qd[i * p + k] * disp[k][nodeD]; dd += q[i * p + k] * disp[k][nodeD]; aa += qa[i * p + k] * disp[k][nodeD]; } a2 += aa * aa; t[i] = i * dt; sg[i] = fac * c; ds[i] = d; vs[i] = w; if (Math.abs(sg[i]) > sMax) sMax = Math.abs(sg[i]); if (Math.abs(dd) > dMax) dMax = Math.abs(dd); if (Math.abs(w) > vMax) vMax = Math.abs(w); v2 += w * w; }
  return { t, sigma: sg, disp: ds, vel: vs, node, x: x[node], sigmaMax: sMax, dispMax: dMax, velMax: vMax, velRms: Math.sqrt(v2 / (steps + 1)), accRms: Math.sqrt(a2 / (steps + 1)), steps, massRatio };
}

// ---- sand erosion with tracked particles in a bend ---------------------------------------------------------------------
/** Wear volume per mass of sand (m³/kg) of one impact at speed v (m/s) and angle alpha (rad from the surface). model: 'dnv' | 'oka' | 'finnie'. */
export function impactWear(model, v, alpha, o = {}) {
  const sa = Math.sin(clamp(alpha, 0, Math.PI / 2));
  if (!(v > 0) || !(sa > 0)) return 0;
  if (model === 'finnie') { const fa = Math.tan(alpha) <= 1 / 3 ? Math.sin(2 * alpha) - 3 * sa * sa : Math.cos(alpha) ** 2 / 3; return (nz(o.finnieC, 0.5) * v * v * fa) / (4 * pos(o.flowStress, 1.96e9)); }
  if (model === 'oka') { const hv = pos(o.hv, 1.96), n1 = OKA.s1 * hv ** OKA.q1, n2 = OKA.s2 * hv ** OKA.q2; return sa ** n1 * (1 + hv * (1 - sa)) ** n2 * OKA.K * hv ** OKA.k1 * (v / OKA.vRef) ** (OKA.k2a * hv ** OKA.k2b) * (pos(o.dp, 250e-6) / OKA.dRef) ** OKA.k3 * 1e-9; }
  return (O501.K * dnvAngle(alpha) * v ** O501.n) / pos(o.rhoT, O501.rhoT);
}
/**
 * Plane potential flow through a bend: stream function on a body-fitted grid (s along the centreline, n across; scale factor h = 1 + κ·n),
 * ∂/∂s(h⁻¹·∂ψ/∂s) + ∂/∂n(h·∂ψ/∂n) = 0, finite volumes and successive over-relaxation. The extrados is n = +D/2.
 * o: { D, R (centreline bend radius), angle (deg), U (mean velocity), ns, nn (cells), legs (straight length before and after the bend, in diameters) }
 * Returns { s[], n[], us[i][j], un[i][j] (m/s), psi[i][j], sIn, sOut (ends of the curved part), length, kappa, iterations, residual, toXY(s, n), toSN(x, y), tangent(s) }.
 */
export function bendFlowField(o) {
  const D = +o.D, R = Math.max(+o.R, 0.51 * D), th = (clamp(nz(o.angle, 90), 1, 170) * Math.PI) / 180, U = +o.U, legs = pos(o.legs, 2) * D, arc = R * th, len = 2 * legs + arc, ns = clamp(Math.round(nz(o.ns, 48)), 12, 400), nn = clamp(Math.round(nz(o.nn, 12)), 4, 100), ds = len / ns, dn = D / nn, kap = 1 / R;
  const s = Array.from({ length: ns + 1 }, (_, i) => i * ds), n = Array.from({ length: nn + 1 }, (_, j) => -D / 2 + j * dn), kAt = (x) => (x > legs && x < legs + arc ? kap : 0), psi = Array.from({ length: ns + 1 }, () => new Float64Array(nn + 1));
  for (let i = 0; i <= ns; i++) for (let j = 0; j <= nn; j++) psi[i][j] = U * (n[j] + D / 2);
  const ke = s.map((x, i) => (i < ns ? kAt(x + ds / 2) : 0)), kn = s.map((x) => 0.5 * (kAt(x - 1e-9 * len) + kAt(x + 1e-9 * len))); // curvature at the cell faces and at the nodes
  let it = 0, res = 1; const om = 2 / (1 + Math.sin(Math.PI / Math.max(ns, nn)));
  for (; it < 6000 && res > 1e-11; it++) {
    res = 0;
    for (let i = 1; i <= ns; i++) for (let j = 1; j < nn; j++) {
      const aE = i < ns ? 1 / ((1 + ke[i] * n[j]) * ds * ds) : 0, aW = 1 / ((1 + ke[i - 1] * n[j]) * ds * ds), aN = (1 + kn[i] * (n[j] + dn / 2)) / (dn * dn), aS = (1 + kn[i] * (n[j] - dn / 2)) / (dn * dn);
      const nv = ((i < ns ? aE * psi[i + 1][j] : 0) + aW * psi[i - 1][j] + aN * psi[i][j + 1] + aS * psi[i][j - 1]) / (aE + aW + aN + aS), d = nv - psi[i][j];
      psi[i][j] += om * d; if (Math.abs(d) > res) res = Math.abs(d);
    }
    res /= Math.abs(U * D) || 1;
  }
  const us = psi.map((r) => new Float64Array(nn + 1)), un = psi.map((r) => new Float64Array(nn + 1));
  for (let i = 0; i <= ns; i++) for (let j = 0; j <= nn; j++) {
    us[i][j] = j === 0 ? (psi[i][1] - psi[i][0]) / dn : j === nn ? (psi[i][nn] - psi[i][nn - 1]) / dn : (psi[i][j + 1] - psi[i][j - 1]) / (2 * dn);
    const dps = i === 0 ? (psi[1][j] - psi[0][j]) / ds : i === ns ? 0 : (psi[i + 1][j] - psi[i - 1][j]) / (2 * ds); un[i][j] = -dps / (1 + kn[i] * n[j]);
  }
  const ct = Math.cos(th), sth = Math.sin(th);
  const toXY = (sv, nv) => { if (sv <= legs) return [sv, nv]; if (sv <= legs + arc) { const f = (sv - legs) / R; return [legs + (R + nv) * Math.sin(f), -R + (R + nv) * Math.cos(f)]; } const a = sv - legs - arc; return [legs + (R + nv) * sth + a * ct, -R + (R + nv) * ct - a * sth]; };
  const toSN = (xv, yv) => { const dx = xv - legs, dy = yv + R; if (dx < 0 && dy > 0) return [xv, yv]; const f = Math.atan2(dx, dy); if (dx >= 0 && f <= th) return [legs + R * f, Math.hypot(dx, dy) - R]; return [legs + arc + dx * ct - dy * sth, dx * sth + dy * ct - R]; };
  const tangent = (sv) => { const f = clamp((sv - legs) / R, 0, th); return [Math.cos(f), -Math.sin(f)]; };
  return { s, n, us, un, psi, sIn: legs, sOut: legs + arc, length: len, kappa: kap, iterations: it, residual: res, toXY, toSN, tangent, D, R, U, ds, dn, angle: th };
}
/**
 * Sand erosion of a bend from tracked particles: Lagrangian particles (Schiller–Naumann drag, optional gravity along −y) in the plane potential-flow
 * field of the bend, wall impacts with restitution, and an impact-wear law applied at every impact. The plane result is converted to the pipe with a scar width of π·D/4.
 * o: { D, rOverD, angle (deg), U, rho, mu (carrier fluid), dp, rhoP, mp (kg/s of sand), model: 'dnv' | 'oka' | 'finnie', hv, flowStress, rhoT, nPart, ns, nn, en, et (restitution), bins, field (optional precomputed bendFlowField) }
 * Returns { rateMax, rateMean (mm/y on the extrados scar), sMax (m along the extrados from the bend inlet), angleMax (deg of bend at the peak), s[] (bin centres, m), rate[] (mm/y), impacts, hitFraction (share of the particles that hit the extrados),
 *   meanAngle (deg), meanSpeed (m/s) (wear-weighted), stokes, wearPerKg (m³/kg, all impacts), field, tracks: [{ x[], y[] }] (a few paths for plotting) }.
 */
export function bendErosionTracked(o) {
  const D = +o.D, U = Math.max(+o.U, 0), rho = pos(o.rho, 100), mu = pos(o.mu, 1e-4), dp = pos(o.dp, 250e-6), rhoP = pos(o.rhoP, 2650), mp = Math.max(nz(o.mp, 0), 0), model = o.model || 'dnv', np = clamp(Math.round(nz(o.nPart, 120)), 8, 5000), en = clamp(nz(o.en, 0.8), 0, 1), et = clamp(nz(o.et, 0.9), 0, 1);
  need(D > 0, 'Particle tracking needs a positive bore.');
  const F = o.field || bendFlowField({ D, R: pos(o.rOverD, 1.5) * D, angle: o.angle, U: Math.max(U, 1e-9), ns: o.ns, nn: o.nn }), ns = F.s.length - 1, nn = F.n.length - 1, half = D / 2 - dp / 2, tau = (rhoP * dp * dp) / (18 * mu), nb = clamp(Math.round(nz(o.bins, 24)), 4, 200);
  const Ro = F.R + D / 2, th = F.angle, bw = (F.R * th + 2 * D) / nb, wear = new Float64Array(nb), hits = new Float64Array(nb), tracks = [], gy = o.gravity ? -G * (1 - rho / rhoP) : 0, wo = { hv: o.hv, flowStress: o.flowStress, rhoT: o.rhoT, finnieC: o.finnieC, dp };
  // allocation-free helpers: (x, y) → (s, n) into cs, cn; local tangent into tx0, ty0; fluid velocity into fu, fv
  const legs = F.sIn, Rc = F.R, cth = Math.cos(th), sth = Math.sin(th), arc = Rc * th, usA = F.us, unA = F.un; let cs = 0, cn = 0, tx0 = 1, ty0 = 0, fu = 0, fv = 0;
  const toSN = (xv, yv) => { const dx = xv - legs, dy = yv + Rc; if (dx < 0 && dy > 0) { cs = xv; cn = yv; return; } const f = Math.atan2(dx, dy); if (dx >= 0 && f <= th) { cs = legs + Rc * f; cn = Math.sqrt(dx * dx + dy * dy) - Rc; return; } cs = legs + arc + dx * cth - dy * sth; cn = dx * sth + dy * cth - Rc; };
  const tangent = (sv) => { const f = sv <= legs ? 0 : sv >= legs + arc ? th : (sv - legs) / Rc; tx0 = Math.cos(f); ty0 = -Math.sin(f); };
  const fluid = (sv, nv) => { let a = sv / F.ds, b = (nv + D / 2) / F.dn; if (a < 0) a = 0; else if (a > ns - 1e-9) a = ns - 1e-9; if (b < 0) b = 0; else if (b > nn - 1e-9) b = nn - 1e-9; const i = Math.floor(a), j = Math.floor(b), fa = a - i, fb = b - j, w00 = (1 - fa) * (1 - fb), w01 = (1 - fa) * fb, w10 = fa * (1 - fb), w11 = fa * fb;
    const u1 = w00 * usA[i][j] + w01 * usA[i][j + 1] + w10 * usA[i + 1][j] + w11 * usA[i + 1][j + 1], u2 = w00 * unA[i][j] + w01 * unA[i][j + 1] + w10 * unA[i + 1][j] + w11 * unA[i + 1][j + 1]; tangent(sv); fu = u1 * tx0 - u2 * ty0; fv = u1 * ty0 + u2 * tx0; }; // normal = (−t_y, t_x)
  const dt = (0.3 * Math.min(F.ds, F.dn)) / Math.max(U, 1e-9), maxStep = Math.ceil((6 * F.length) / (Math.max(U, 1e-9) * dt)); let impacts = 0, hitP = 0, wSum = 0, aSum = 0, vSum = 0;
  if (U > 0) for (let k = 0; k < np; k++) {
    const n0 = -half + ((k + 0.5) / np) * 2 * half, share = 1 / np; let x = 0, y = n0, vx = U, vy = 0, hit = false; const keep = k % Math.max(1, Math.floor(np / 12)) === 0, tx = [x], ty = [y];
    for (let step = 0, bounces = 0; step < maxStep && bounces < 30; step++) {
      toSN(x, y); if (cs >= F.length || cs < -1e-9) break;
      toSN(x + 0.5 * dt * vx, y + 0.5 * dt * vy); fluid(cs, cn > half ? half : cn < -half ? -half : cn); // fluid velocity at the mid-point of the step
      const rx = fu - vx, ry = fv - vy, rel = Math.sqrt(rx * rx + ry * ry), fd = 1 + 0.15 * ((rho * rel * dp) / mu) ** 0.687, e = Math.exp((-dt * fd) / tau), nvx = fu - rx * e, nvy = fv - ry * e + gy * dt;
      let xn = x + 0.5 * (vx + nvx) * dt, yn = y + 0.5 * (vy + nvy) * dt; vx = nvx; vy = nvy;
      toSN(xn, yn);
      if (cn > half || cn < -half) { // wall impact: velocity components along the wall and into it
        const s2 = cs, outer = cn > 0; tangent(s2); const vt = vx * tx0 + vy * ty0, vn = -vx * ty0 + vy * tx0, into = outer ? vn : -vn;
        if (into > 0) {
          const speed = Math.sqrt(vt * vt + vn * vn), alpha = Math.atan2(into, Math.abs(vt)), wv = impactWear(model, speed, alpha, wo) * share; bounces++; impacts++;
          if (outer) { hit = true; const b = Math.floor((s2 - F.sIn + D) / bw); if (b >= 0 && b < nb) { wear[b] += wv; hits[b] += share; } wSum += wv; aSum += wv * alpha; vSum += wv * speed; }
          const vn2 = -en * vn, vt2 = et * vt; vx = vt2 * tx0 - vn2 * ty0; vy = vt2 * ty0 + vn2 * tx0;
        }
        const pw = F.toXY(s2, outer ? half : -half); xn = pw[0]; yn = pw[1];
      }
      x = xn; y = yn; if (keep && step % 4 === 0) { tx.push(x); ty.push(y); }
    }
    if (hit) hitP++; if (keep) { tx.push(x); ty.push(y); tracks.push({ x: tx, y: ty }); }
  }
  const width = (Math.PI * D) / 4, sC = Array.from(wear, (_, b) => (b + 0.5) * bw - D), rate = Array.from(wear, (w, b) => (w * mp * 1000 * YEAR) / (bw * (sC[b] > 0 && sC[b] < F.R * th ? Ro / F.R : 1) * width)); let iMax = 0; rate.forEach((r, b) => { if (r > rate[iMax]) iMax = b; });
  const scar = rate.filter((r) => r > 0);
  return { rateMax: rate[iMax] || 0, rateMean: scar.length ? sum(scar) / scar.length : 0, sMax: sC[iMax], angleMax: clamp((sC[iMax] / F.R) * (180 / Math.PI), 0, (th * 180) / Math.PI), s: sC, rate, impacts, hitFraction: hitP / np, meanAngle: wSum > 0 ? ((aSum / wSum) * 180) / Math.PI : 0, meanSpeed: wSum > 0 ? vSum / wSum : 0,
    stokes: (tau * U) / D, wearPerKg: wSum, field: F, tracks, binWidth: bw, scarWidth: width, nPart: np };
}

// ---- suite declaration --------------------------------------------------------------------------------------------------
// ---- resolved fluid–structure interaction: axisymmetric Navier–Stokes in an elastic tube ------------------------------------
/** Symmetric positive-definite band solver. A holds the lower band row-wise (row i, columns i−bw … i at offsets 0 … bw); factorises in place and returns b -> x (b is overwritten). */
export function bandCholesky(A, n, bw) {
  const w1 = bw + 1;
  for (let i = 0; i < n; i++) {
    const ri = i * w1 + bw - i, j0 = Math.max(0, i - bw);
    for (let j = j0; j <= i; j++) {
      const rj = j * w1 + bw - j; let s = A[ri + j];
      for (let k = Math.max(j0, j - bw); k < j; k++) s -= A[ri + k] * A[rj + k];
      if (j === i) { need(s > 0, 'The band matrix is not positive definite (the model is not restrained or the pressure level is not fixed).'); A[ri + i] = Math.sqrt(s); } else A[ri + j] = s / A[rj + j];
    }
  }
  return (u) => {
    for (let i = 0; i < n; i++) { const ri = i * w1 + bw - i; let s = u[i]; for (let k = Math.max(0, i - bw); k < i; k++) s -= A[ri + k] * u[k]; u[i] = s / A[ri + i]; }
    for (let i = n - 1; i >= 0; i--) { let s = u[i]; const hi = Math.min(n - 1, i + bw); for (let k = i + 1; k <= hi; k++) s -= A[k * w1 + bw - k + i] * u[k]; u[i] = s / A[i * w1 + bw]; }
    return u;
  };
}
// complex helpers for the Bessel functions of the oscillatory-flow solution: numbers are [re, im]
const cmul = (a, b) => [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]], cdiv = (a, b) => { const d = b[0] * b[0] + b[1] * b[1]; return [(a[0] * b[0] + a[1] * b[1]) / d, (a[1] * b[0] - a[0] * b[1]) / d]; };
const csqrt = (a) => { const r = Math.hypot(a[0], a[1]), re = Math.sqrt(Math.max((r + a[0]) / 2, 0)), im = Math.sqrt(Math.max((r - a[0]) / 2, 0)); return [re, a[1] < 0 ? -im : im]; };
/** Bessel functions J0 and J1 of a complex argument by their power series (accurate for |z| up to about 20). */
export function besselJ01(z) {
  const q = cmul(z, z), m = [-q[0] / 4, -q[1] / 4]; let t0 = [1, 0], t1 = [1, 0], j0 = [1, 0], j1 = [1, 0];
  for (let k = 1; k < 120; k++) { t0 = cmul(t0, [m[0] / (k * k), m[1] / (k * k)]); t1 = cmul(t1, [m[0] / (k * (k + 1)), m[1] / (k * (k + 1))]); j0 = [j0[0] + t0[0], j0[1] + t0[1]]; j1 = [j1[0] + t1[0], j1[1] + t1[1]]; if (Math.hypot(t0[0], t0[1]) < 1e-17 * Math.hypot(j0[0], j0[1]) && k > 4) break; }
  return { j0, j1: cmul([z[0] / 2, z[1] / 2], j1) };
}
/** Modified Bessel functions I0 and I1 of a real argument (power series). */
export function besselI01(x) { let t0 = 1, t1 = 1, i0 = 1, i1 = 1; const m = (x * x) / 4; for (let k = 1; k < 200; k++) { t0 *= m / (k * k); t1 *= m / (k * (k + 1)); i0 += t0; i1 += t1; if (t0 < 1e-17 * i0) break; } return { i0, i1: (x / 2) * i1 }; }
/**
 * Waves in a liquid-filled elastic tube: the reference solutions used to verify the coupled solver.
 * o: { R (bore radius), h (wall), E, nu, rhoS, rhoF, mu, K (bulk modulus of the liquid; Infinity = incompressible), axial: 'restrained' | 'open', k (wavenumber, 1/m; for the dispersion results) }
 * Returns { c0 (Moens–Korteweg with the mid-surface radius and the axial-restraint factor), cMK (classical √(E·h/(2·ρ·R))), cKorteweg (with the compressibility of the liquid),
 *   kEff (wall stiffness seen by the bore, Pa/m), omegaLamb (rad/s: inviscid standing wave of wavenumber k with the radial inertia of the liquid, the wall inertia and wall bending),
 *   cLamb (omegaLamb / k), omegaWomersley: [re, im] (complex frequency of the free wave in a viscous liquid, tethered thin tube: ω² = c0²k²(1 − F10(α))·…), alpha (Womersley number at that frequency) }.
 */
export function tubeWaves(o) {
  const R = +o.R, h = +o.h, E = +o.E, nu = nz(o.nu, 0.3), rho = +o.rhoF, rhoS = nz(o.rhoS, 7850), Rm = R + h / 2, restr = o.axial !== 'open', ks = (E * h) / (Rm * Rm * (restr ? 1 - nu * nu : 1)), kEff = (ks * Rm) / R, Kb = o.K > 0 ? +o.K : Infinity;
  const c0 = Math.sqrt((R * kEff) / (2 * rho)), cMK = Math.sqrt((E * h) / (2 * rho * R)), cK = Number.isFinite(Kb) ? Math.sqrt(1 / (rho / Kb + 1 / (c0 * c0))) : c0, out = { c0, cMK, cKorteweg: cK, kEff, ks };
  if (o.k > 0) {
    const k = +o.k, D = (E * h ** 3) / (12 * (1 - nu * nu)), B = besselI01(k * R), added = (rho * B.i0) / (k * B.i1), mw = (rhoS * h * Rm) / R, stiff = kEff * (1 + (D * k ** 4) / ks);
    out.omegaLamb = Math.sqrt(stiff / (added + mw)); out.cLamb = out.omegaLamb / k; out.addedMass = added; out.wallMass = mw;
    if (o.mu > 0) { // Womersley: long-wave momentum with the oscillatory velocity profile; wall inertia and bending kept, radial inertia of the liquid neglected
      let w = [out.omegaLamb, -0.01 * out.omegaLamb];
      for (let it = 0; it < 200; it++) {
        const a2 = [0, (-R * R * rho) / o.mu], z = csqrt(cmul(a2, w)), J = besselJ01(z), F = cdiv([2 * J.j1[0], 2 * J.j1[1]], cmul(z, J.j0)), g = [1 - F[0], -F[1]]; // z = α·i^{3/2} with α² = R²ωρ/μ; time factor exp(iωt)
        // ω²·(2ρ/(R k² g) + m_w) = stiff  →  ω = √(stiff / (2ρ/(R k² g) + m_w))
        const den = cdiv([(2 * rho) / (R * k * k), 0], g), wn = csqrt(cdiv([stiff, 0], [den[0] + mw, den[1]]));
        const d = Math.hypot(wn[0] - w[0], wn[1] - w[1]); w = [0.5 * (w[0] + wn[0]), 0.5 * (w[1] + wn[1])]; if (d < 1e-13 * Math.abs(w[0])) break;
      }
      out.omegaWomersley = w; out.alpha = R * Math.sqrt((Math.abs(w[0]) * rho) / o.mu);
    }
  }
  return out;
}
/**
 * Strongly coupled fluid–structure interaction of a liquid in an elastic tube.
 * Fluid: axisymmetric Navier–Stokes equations on a grid that follows the wall (r = ξ·R(x, t); arbitrary Lagrangian–Eulerian form, staggered velocities, finite volumes whose
 * faces move with the wall so that the continuity equation is exact on the moving grid), explicit upwind convection and viscous terms, pressure from the exact discrete
 * pressure equation (optionally with the compressibility of the liquid), time-centred pressure. Metric terms of the momentum equation of the order of the wall slope are neglected.
 * Structure: thin cylindrical shell with radial inertia, hoop stiffness and meridional bending (Hermite elements along the tube, trapezoidal rule in time).
 * Coupling: Dirichlet–Neumann sub-iterations in every time step (the wall velocity is imposed on the fluid, the wall pressure loads the shell). The liquid is nearly incompressible, so its
 * added mass makes the unrelaxed iteration diverge; Aitken's dynamic relaxation converges for short tubes, and the default is a quasi-Newton iteration on the interface residual.
 * o: { R, L, h, E, nu, rhoS, rhoF, mu, K (bulk modulus, Pa; omitted = incompressible), axial: 'restrained' | 'open', nx, nr, dt, tEnd,
 *      inlet, outlet: { type: 'pressure', p: (t) => Pa } | { type: 'wall' }, wallEnds: 'free' | 'clamped', eta0: (x) => m (initial wall displacement, liquid at rest),
 *      coupling: 'newton' (default: interface Newton iteration with a frozen finite-difference Jacobian) | 'aitken' | 'fixed' (constant relaxation o.omega, default 1) | 'one' (rigid wall for the fluid, the pressure loads the shell), tol (relative, default 1e-8), maxIter (default 60),
 *      damping (ratio of critical at the ring frequency, default 0), stations: [x] (histories are kept there), convection (default true), progress(f) }
 * Returns { x[] (nodes), t[], eta[station][] (m), pWall[station][] (Pa), q[station][] (m³/s), etaEnd[], pEnd[] (cell wall pressure), iterations: { mean, max, total }, converged (all steps), diverged (bool),
 *   energy: { fluid[], shell[], interface (work of the liquid on the wall, fluid side, J), shellGain (J, with the damping loss), boundary (work of the end pressures, J), dissipation (J, from the balance) },
 *   volume: { fluidChange, wallSwept, boundary (m³: net inflow) }, divergenceMax (1/s), steps, dt, grid: { nx, nr } }.
 */
export function tubeFsi(o) {
  const R0 = +o.R, L = +o.L, h = +o.h, E = +o.E, nu = nz(o.nu, 0.3), rhoS = nz(o.rhoS, 7850), rho = +o.rhoF, mu = Math.max(nz(o.mu, 0), 0), nuK = mu / rho, Kb = o.K > 0 ? +o.K : Infinity;
  need(R0 > 0 && L > 0 && h > 0 && E > 0 && rho > 0, 'The tube model needs positive radius, length, wall thickness, modulus and liquid density.');
  const nx = clamp(Math.round(nz(o.nx, 60)), 8, 600), nr = clamp(Math.round(nz(o.nr, 8)), 2, 48), dx = L / nx, dxi = 1 / nr, Rm = R0 + h / 2, restr = o.axial !== 'open';
  const ks = (E * h) / (Rm * Rm * (restr ? 1 - nu * nu : 1)), Db = (E * h ** 3) / (12 * (1 - nu * nu)), ms = rhoS * h, lf = R0 / Rm; // per unit area of the mid-surface; lf converts the bore pressure
  let dt = +o.dt; need(dt > 0 && o.tEnd > 0, 'The tube model needs a positive time step and duration.');
  if (nuK > 0) dt = Math.min(dt, (0.2 * (R0 * dxi) ** 2) / nuK, (0.2 * dx * dx) / nuK);
  const steps = Math.max(1, Math.round(o.tEnd / dt)), mode = o.coupling || 'newton', tol = pos(o.tol, 1e-8), maxIter = Math.round(pos(o.maxIter, 60)), conv = o.convection !== false;
  const inl = o.inlet || { type: 'pressure', p: () => 0 }, outl = o.outlet || { type: 'pressure', p: () => 0 };
  // ---- shell: constant matrices, banded (half bandwidth 3)
  const nd = nx + 1, ns = 2 * nd, Ks = zeros(ns), Ms = zeros(ns);
  for (let e = 0; e < nx; e++) {
    const l = dx, c = Db / l ** 3, f = (ks * l) / 420, m = (ms * l) / 420;
    const ke = [[12, 6 * l, -12, 6 * l], [6 * l, 4 * l * l, -6 * l, 2 * l * l], [-12, -6 * l, 12, -6 * l], [6 * l, 2 * l * l, -6 * l, 4 * l * l]], me = [[156, 22 * l, 54, -13 * l], [22 * l, 4 * l * l, 13 * l, -3 * l * l], [54, 13 * l, 156, -22 * l], [-13 * l, -3 * l * l, -22 * l, 4 * l * l]];
    for (let a = 0; a < 4; a++) for (let b = 0; b < 4; b++) { Ks[2 * e + a][2 * e + b] += c * ke[a][b] + f * me[a][b]; Ms[2 * e + a][2 * e + b] += m * me[a][b]; }
  }
  const fixedS = new Uint8Array(ns); if (o.wallEnds === 'clamped') { fixedS[0] = fixedS[1] = fixedS[ns - 2] = fixedS[ns - 1] = 1; }
  const zeta = Math.max(nz(o.damping, 0), 0), cS = 2 * zeta * Math.sqrt(ks / ms); // mass-proportional damping coefficient (1/s)
  const Seff = new Float64Array(ns * 4), sIdx = (i, j) => i * 4 + 3 - i + j;
  for (let i = 0; i < ns; i++) for (let j = Math.max(0, i - 3); j <= i; j++) Seff[sIdx(i, j)] = fixedS[i] || fixedS[j] ? (i === j ? 1 : 0) : (2 / (dt * dt) + cS / dt) * Ms[i][j] + 0.5 * Ks[i][j];
  const solveS = bandCholesky(Seff, ns, 3);
  const bandMul = (A, x) => { const y = new Float64Array(ns); for (let i = 0; i < ns; i++) { let s = 0; const Ai = A[i]; for (let j = Math.max(0, i - 3); j <= Math.min(ns - 1, i + 3); j++) s += Ai[j] * x[j]; y[i] = s; } return y; };
  // ---- fluid arrays
  const nc = nx * nr, U = new Float64Array((nx + 1) * nr), V = new Float64Array(nx * (nr + 1)), Us = new Float64Array(U.length), Vs = new Float64Array(V.length), Un = new Float64Array(U.length), Vn = new Float64Array(V.length), P = new Float64Array(nc), Pold = new Float64Array(nc), rhs = new Float64Array(nc);
  const ui = (i, j) => i * nr + j, vi = (i, j) => i * (nr + 1) + j, pi2 = (i, j) => i * nr + j, xiF = Float64Array.from({ length: nr + 1 }, (_, j) => j * dxi), xiC = Float64Array.from({ length: nr }, (_, j) => (j + 0.5) * dxi), ring = Float64Array.from({ length: nr }, (_, j) => 0.5 * (xiF[j + 1] ** 2 - xiF[j] ** 2));
  let w = new Float64Array(ns), wd = new Float64Array(ns); // shell displacement and velocity (deflection, slope pairs)
  if (o.eta0) { for (let i = 0; i < nd; i++) { const x = i * dx, e = 1e-6 * L; w[2 * i] = o.eta0(x); w[2 * i + 1] = (o.eta0(Math.min(x + e, L)) - o.eta0(Math.max(x - e, 0))) / (Math.min(x + e, L) - Math.max(x - e, 0)); } for (let d = 0; d < ns; d++) if (fixedS[d]) w[d] = 0; }
  const Rn = new Float64Array(nd), Rc = new Float64Array(nx), wallFlux = new Float64Array(nx), A = new Float64Array(nc * (nr + 1)), aIdx = (i, j) => i * (nr + 1) + nr - i + j;
  const cellInt = (vec, e) => (dx / 2) * (vec[2 * e] + vec[2 * e + 2]) + ((dx * dx) / 12) * (vec[2 * e + 1] - vec[2 * e + 3]); // ∫ of the Hermite field over element e
  const geom = (ww) => { for (let i = 0; i < nd; i++) Rn[i] = R0 + ww[2 * i]; for (let e = 0; e < nx; e++) Rc[e] = R0 + cellInt(ww, e) / dx; };
  const volume = () => { let s = 0; for (let e = 0; e < nx; e++) s += Math.PI * dx * Rc[e] * Rc[e]; return s; };
  const kinetic = () => { let s = 0; for (let i = 0; i <= nx; i++) { const wx = i === 0 || i === nx ? 0.5 : 1; for (let j = 0; j < nr; j++) s += wx * ring[j] * Rn[i] * Rn[i] * dx * U[ui(i, j)] ** 2; } for (let i = 0; i < nx; i++) for (let j = 1; j < nr; j++) s += xiF[j] * dxi * Rc[i] * Rc[i] * dx * V[vi(i, j)] ** 2; return Math.PI * rho * s; };
  const shellEnergy = () => { const Kw = bandMul(Ks, w), Mv = bandMul(Ms, wd); return 2 * Math.PI * Rm * 0.5 * (dot(w, Kw) + dot(wd, Mv)); };
  // explicit convective and viscous increments into Us, Vs (velocities at the old time level, grid velocity of the current iterate)
  const explicit = (wdNew) => {
    for (let i = 0; i <= nx; i++) for (let j = 0; j < nr; j++) {
      const k = ui(i, j), u = U[k], R = Rn[i], uW = i > 0 ? U[ui(i - 1, j)] : u, uE = i < nx ? U[ui(i + 1, j)] : u, uS = j > 0 ? U[ui(i, j - 1)] : u, uN = j < nr - 1 ? U[ui(i, j + 1)] : -u; // no slip at the wall, symmetry on the axis
      let acc = nuK * ((uE - 2 * u + uW) / (dx * dx) + (xiF[j + 1] * (uN - u) - xiF[j] * (u - uS)) / (xiC[j] * dxi * dxi * R * R));
      if (conv) {
        const il = Math.max(i - 1, 0), ir = Math.min(i, nx - 1), q = 0.25 * (V[vi(il, j)] + V[vi(il, j + 1)] + V[vi(ir, j)] + V[vi(ir, j + 1)]), W = q - xiC[j] * 0.5 * (wdNew[2 * i] + wd[2 * i]);
        acc -= u * (u > 0 ? (u - uW) / dx : (uE - u) / dx) + (W / R) * (W > 0 ? (u - uS) / dxi : (uN - u) / dxi);
      }
      Us[k] = u + dt * acc;
    }
    for (let i = 0; i < nx; i++) for (let j = 1; j < nr; j++) {
      const k = vi(i, j), v = V[k], R = Rc[i], vW = i > 0 ? V[vi(i - 1, j)] : v, vE = i < nx - 1 ? V[vi(i + 1, j)] : v, vS = V[vi(i, j - 1)], vN = V[vi(i, j + 1)];
      let acc = nuK * ((vE - 2 * v + vW) / (dx * dx) + ((xiF[j + 1] * vN - xiF[j] * v) / xiC[j] - (xiF[j] * v - xiF[j - 1] * vS) / xiC[j - 1]) / (dxi * dxi * R * R));
      if (conv) { const u = 0.25 * (U[ui(i, j - 1)] + U[ui(i, j)] + U[ui(i + 1, j - 1)] + U[ui(i + 1, j)]), W = v - xiF[j] * 0.25 * (wdNew[2 * i] + wdNew[2 * i + 2] + wd[2 * i] + wd[2 * i + 2]); acc -= u * (u > 0 ? (v - vW) / dx : (vE - v) / dx) + (W / R) * (W > 0 ? (v - vS) / dxi : (vN - v) / dxi); }
      Vs[k] = v + dt * acc;
    }
  };
  // one fluid step for a wall iterate: returns the pressure load vector of the shell; leaves the new velocities in Un, Vn and the time-centred pressure in P
  const comp = Number.isFinite(Kb) ? 2 / (Kb * dt) : 0, pw = new Float64Array(nx);
  const fluid = (wNew, wdNew, t) => {
    geom(wNew); explicit(wdNew);
    const pIn = inl.type === 'pressure' ? inl.p(t) : 0, pOut = outl.type === 'pressure' ? outl.p(t) : 0, c = dt / rho; A.fill(0);
    for (let e = 0; e < nx; e++) wallFlux[e] = Rc[e] * cellInt(wdNew, e);
    for (let i = 0; i < nx; i++) for (let j = 0; j < nr; j++) {
      const k = pi2(i, j), aW = ring[j] * Rn[i] * Rn[i], aE = ring[j] * Rn[i + 1] * Rn[i + 1], gS = xiF[j] * Rc[i] * dx, gN = xiF[j + 1] * Rc[i] * dx, vol = ring[j] * Rc[i] * Rc[i] * dx;
      let diag = comp * vol, b = aW * Us[ui(i, j)] - aE * Us[ui(i + 1, j)] + gS * (j > 0 ? Vs[vi(i, j)] : 0) - (j < nr - 1 ? gN * Vs[vi(i, j + 1)] : wallFlux[i]) + comp * vol * Pold[k];
      if (i > 0) { const cf = (c * aW) / dx; diag += cf; A[aIdx(k, k - nr)] -= cf; } else if (inl.type === 'pressure') { const cf = (c * aW) / (dx / 2); diag += cf; b += cf * pIn; } else b -= aW * Us[ui(0, j)] - 0; // closed end: no flux (Us there is ignored below)
      if (i < nx - 1) diag += (c * aE) / dx; else if (outl.type === 'pressure') { const cf = (c * aE) / (dx / 2); diag += cf; b += cf * pOut; } else b += aE * Us[ui(nx, j)];
      if (j > 0) { const cf = (c * gS) / (Rc[i] * dxi); diag += cf; A[aIdx(k, k - 1)] -= cf; }
      if (j < nr - 1) diag += (c * gN) / (Rc[i] * dxi);
      A[aIdx(k, k)] = diag; rhs[k] = b;
    }
    bandCholesky(A, nc, nr)(rhs); P.set(rhs);
    for (let i = 0; i <= nx; i++) for (let j = 0; j < nr; j++) {
      let g;
      if (i === 0) g = inl.type === 'pressure' ? (P[pi2(0, j)] - pIn) / (dx / 2) : null; else if (i === nx) g = outl.type === 'pressure' ? (pOut - P[pi2(nx - 1, j)]) / (dx / 2) : null; else g = (P[pi2(i, j)] - P[pi2(i - 1, j)]) / dx;
      Un[ui(i, j)] = g === null ? 0 : Us[ui(i, j)] - c * g;
    }
    for (let i = 0; i < nx; i++) { Vn[vi(i, 0)] = 0; Vn[vi(i, nr)] = wallFlux[i] / (Rc[i] * dx); for (let j = 1; j < nr; j++) Vn[vi(i, j)] = Vs[vi(i, j)] - (c * (P[pi2(i, j)] - P[pi2(i, j - 1)])) / (Rc[i] * dxi); }
    const f = new Float64Array(ns);
    for (let e = 0; e < nx; e++) { const q = (pw[e] = 1.5 * P[pi2(e, nr - 1)] - 0.5 * P[pi2(e, nr - 2)]) * lf; f[2 * e] += (q * dx) / 2; f[2 * e + 1] += (q * dx * dx) / 12; f[2 * e + 2] += (q * dx) / 2; f[2 * e + 3] -= (q * dx * dx) / 12; }
    return f;
  };
  const structure = (f) => { // trapezoidal rule: (2M/dt² + C/dt + K/2)·w₁ = f + M·(2w₀/dt² + 2ẇ₀/dt) + C·w₀/dt − K·w₀/2
    const a = new Float64Array(ns); for (let d = 0; d < ns; d++) a[d] = (2 / (dt * dt) + cS / dt) * w[d] + (2 / dt) * wd[d];
    const Ma = bandMul(Ms, a), Kw = bandMul(Ks, w), b = new Float64Array(ns); for (let d = 0; d < ns; d++) b[d] = fixedS[d] ? 0 : f[d] + Ma[d] - 0.5 * Kw[d];
    return solveS(b);
  };
  const st = (o.stations && o.stations.length ? o.stations : [L / 2]).map((x) => clamp(Math.round(x / dx), 0, nx)), hist = { t: [0], eta: st.map((i) => [w[2 * i]]), pWall: st.map(() => [0]), q: st.map(() => [0]) };
  geom(w); const vol0 = volume(), eF = [0], eS = [shellEnergy()]; let itTot = 0, itMax = 0, allConv = true, diverged = false, wInt = 0, wBnd = 0, wDamp = 0, swept = 0, inflow = 0, divMax = 0, done = 0;
  let omega = pos(o.omega0, 0.1); let jac = null, rebuilt = false, nJac = 0; const wNew = new Float64Array(ns), wdNew = new Float64Array(ns), res = new Float64Array(ns), resOld = new Float64Array(ns);
  for (let s = 1; s <= steps && !diverged; s++) {
    const t = (s - 0.5) * dt; if (mode === 'fixed') omega = pos(o.omega, 1); let it = 0, ok = false, f = null, wSol = null;
    for (let d = 0; d < ns; d++) wNew[d] = fixedS[d] ? 0 : w[d] + dt * wd[d]; // predictor
    if (mode === 'one') { for (let d = 0; d < ns; d++) { wNew[d] = 0; wdNew[d] = 0; } f = fluid(wNew, wdNew, t); wSol = structure(f); ok = true; it = 1; }
    else for (; it < maxIter; ) {
      for (let d = 0; d < ns; d++) wdNew[d] = (2 * (wNew[d] - w[d])) / dt - wd[d];
      f = fluid(wNew, wdNew, t); wSol = structure(f); it++;
      let rn = 0, sn = 0; for (let d = 0; d < ns; d++) { res[d] = wSol[d] - wNew[d]; if (d % 2 === 0) { rn += res[d] * res[d]; sn += wSol[d] * wSol[d]; } }
      if (!Number.isFinite(rn) || rn > 1e60) { diverged = true; break; }
      if (Math.sqrt(rn) <= tol * Math.max(Math.sqrt(sn), 1e-12 * R0 * Math.sqrt(nd))) { ok = true; break; }
      if (mode === 'newton') { // interface Newton iteration with a frozen Jacobian of the residual, built once by finite differences (one fluid and one shell solve per wall unknown)
        if (!jac || (it > 6 && !rebuilt)) {
          if (jac) rebuilt = true;
          const J = zeros(ns), base = Float64Array.from(res), x0 = Float64Array.from(wNew);
          for (let c = 0; c < ns; c++) {
            if (fixedS[c]) { J[c][c] = -1; continue; }
            const del = (c % 2 === 0 ? 1e-5 * R0 : (1e-5 * R0) / dx); x0[c] += del; for (let d = 0; d < ns; d++) wdNew[d] = (2 * (x0[d] - w[d])) / dt - wd[d];
            const sc = structure(fluid(x0, wdNew, t)); for (let d = 0; d < ns; d++) J[d][c] = (sc[d] - x0[d] - base[d]) / del; x0[c] -= del;
          }
          jac = luSolver(J); nJac++;
        }
        const dxn = jac(res); for (let d = 0; d < ns; d++) if (!fixedS[d]) wNew[d] -= dxn[d];
        continue;
      }
      if (mode === 'aitken' && it > 1) { let num = 0, den = 0; for (let d = 0; d < ns; d++) { const dr = res[d] - resOld[d]; num += resOld[d] * dr; den += dr * dr; } if (den > 0) { const om = (-omega * num) / den; if (Number.isFinite(om) && om !== 0) omega = Math.sign(om) * clamp(Math.abs(om), 1e-4, 2); } }
      for (let d = 0; d < ns; d++) { wNew[d] += omega * res[d]; resOld[d] = res[d]; }
    }
    if (diverged) break;
    if (mode !== 'one' && !ok) allConv = false;
    itTot += it; itMax = Math.max(itMax, it);
    // accept the step: the converged wall iterate is the structural solution
    const wAcc = mode === 'one' ? wSol : wNew, wdAcc = new Float64Array(ns); for (let d = 0; d < ns; d++) wdAcc[d] = (2 * (wAcc[d] - w[d])) / dt - wd[d];
    // energy and volume bookkeeping with the time-centred pressure
    let qIn = 0, qOut = 0, pIn = inl.type === 'pressure' ? inl.p(t) : 0, pOut = outl.type === 'pressure' ? outl.p(t) : 0;
    for (let j = 0; j < nr; j++) { qIn += Math.PI * ring[j] * Rn[0] * Rn[0] * (Un[ui(0, j)] + U[ui(0, j)]); qOut += Math.PI * ring[j] * Rn[nx] * Rn[nx] * (Un[ui(nx, j)] + U[ui(nx, j)]); }
    wBnd += dt * (pIn * qIn - pOut * qOut); inflow += dt * (qIn - qOut);
    if (mode !== 'one') { for (let e = 0; e < nx; e++) { const dv = 2 * Math.PI * R0 * (cellInt(wAcc, e) - cellInt(w, e)); wInt += pw[e] * dv; swept += dv; } const vm = new Float64Array(ns); for (let d = 0; d < ns; d++) vm[d] = 0.5 * (wd[d] + wdAcc[d]); wDamp += 2 * Math.PI * Rm * cS * dt * dot(vm, bandMul(Ms, vm)); }
    else { const vm = new Float64Array(ns); for (let d = 0; d < ns; d++) vm[d] = 0.5 * (wd[d] + wdAcc[d]); wDamp += 2 * Math.PI * Rm * cS * dt * dot(vm, bandMul(Ms, vm)); for (let e = 0; e < nx; e++) wInt += pw[e] * 2 * Math.PI * R0 * (cellInt(wAcc, e) - cellInt(w, e)); }
    if (comp > 0) for (let k = 0; k < nc; k++) Pold[k] = 2 * P[k] - Pold[k];
    U.set(Un); V.set(Vn); w = Float64Array.from(wAcc); wd = wdAcc; geom(w); done = s;
    if (s % Math.max(1, Math.floor(steps / 40)) === 0 || s === steps) { // discrete divergence of the accepted velocity field (incompressible limit)
      if (!(comp > 0)) for (let i = 0; i < nx; i++) for (let j = 0; j < nr; j++) { const dv = ring[j] * (Rn[i + 1] ** 2 * U[ui(i + 1, j)] - Rn[i] ** 2 * U[ui(i, j)]) + Rc[i] * dx * (xiF[j + 1] * V[vi(i, j + 1)] - xiF[j] * V[vi(i, j)]); divMax = Math.max(divMax, Math.abs(dv) / (ring[j] * Rc[i] * Rc[i] * dx)); }
      if (o.progress) o.progress(s / steps);
    }
    hist.t.push(s * dt); eF.push(kinetic()); eS.push(shellEnergy());
    st.forEach((i, k) => { hist.eta[k].push(w[2 * i]); const e = Math.min(i, nx - 1); hist.pWall[k].push(pw[e]); let q = 0; for (let j = 0; j < nr; j++) q += 2 * Math.PI * ring[j] * Rn[i] * Rn[i] * U[ui(i, j)]; hist.q[k].push(q); });
  }
  const eta = Array.from({ length: nd }, (_, i) => w[2 * i]), n1 = eF.length - 1;
  return { x: Array.from({ length: nd }, (_, i) => i * dx), t: hist.t, eta: hist.eta, pWall: hist.pWall, q: hist.q, etaEnd: eta, pEnd: Array.from(pw), iterations: { mean: done ? itTot / done : 0, max: itMax, total: itTot, jacobians: nJac }, converged: allConv && !diverged, diverged,
    energy: { fluid: eF, shell: eS, interface: wInt, shellGain: eS[n1] - eS[0] + wDamp, boundary: wBnd, dissipation: wBnd - (eF[n1] - eF[0]) - wInt, damping: wDamp },
    volume: { fluidChange: volume() - vol0, wallSwept: swept, boundary: inflow }, divergenceMax: divMax, steps: done, dt, grid: { nx, nr } };
}
// ---- resolved flow through a bend (plane, body-fitted grid, variable density) and its loads on the structure --------------------
/**
 * Plane incompressible Navier–Stokes flow through a pipe bend on the body-fitted orthogonal grid of bendFlowField (s along the centreline, n across, scale factor h = 1 + κ·n,
 * extrados at n = +D/2): staggered finite volumes, explicit upwind convection with the centrifugal and Coriolis terms of the curved grid, mixing-length eddy viscosity with a
 * logarithmic wall law (or laminar no-slip walls), variable density carried by a limited (van Leer) conservative scheme — a slug passes as a density wave —, and a projection step
 * with a direct band solution of the variable-density pressure equation. The plane result stands for a pipe through the depth b = π·D/4 (the area of the bore over its diameter).
 * The grid can ride on a moving pipe: `motion` gives the acceleration, the rotation rate and the change of curvature of every section, which enter as body forces and metric changes.
 * o: { D, R (centreline radius), angle (deg), legs (straight length each side, in diameters), U (bulk velocity; number or function of time), rho, mu, rhoSlug, slug: { t0 (s), len (m) } | null,
 *      pOut (gauge pressure at the outlet, Pa), ns, nn (cells along and across), cfl (default 0.4), turbulence (default true), gravity: [gx, gy] (plane of the bend) }
 * Returns a simulation object { grid, step(dt, motion) -> loads (tentative), commit(), time, dtStable(), loads(), fields() }; loads = { s[], pExt[], pInt[], tauExt[], tauInt[] (Pa), fn[], ft[] (N per metre of centreline,
 *   on the pipe: normal towards the extrados, tangential downstream), force: [Fx, Fy] (N, integrated wall tractions), cv: [Fx, Fy] (N, control-volume momentum theorem at the same instant),
 *   pIn (Pa, section mean), momentum: [Mx, My] (N·s), mass (kg), divergenceMax }.
 */
export function bendCfdSim(o) {
  const D = +o.D, R = Math.max(+o.R, 0.6 * D), th = (clamp(nz(o.angle, 90), 0, 175) * Math.PI) / 180, legs = pos(o.legs, 3) * D, arc = R * th, len = 2 * legs + arc, rho0 = pos(o.rho, 1000), mu = pos(o.mu, 1e-3), rhoS = pos(o.rhoSlug, rho0), b = (Math.PI * D) / 4;
  need(D > 0, 'The bend flow model needs a positive bore.');
  const nsT = clamp(Math.round(nz(o.ns, 64)), 12, 600), nn = clamp(Math.round(nz(o.nn, 10)), 4, 64), nA = th > 0 ? Math.max(4, Math.round((nsT * arc) / len)) : 0, nL = Math.max(3, Math.round((nsT - nA) / 2)), ns = 2 * nL + nA, dn = D / nn;
  const sF = new Float64Array(ns + 1), kap0 = new Float64Array(ns); for (let i = 0; i <= ns; i++) sF[i] = i <= nL ? (legs * i) / nL : i <= nL + nA ? legs + (arc * (i - nL)) / nA : legs + arc + (legs * (i - nL - nA)) / nL;
  for (let i = nL; i < nL + nA; i++) kap0[i] = 1 / R;
  const ds = Float64Array.from({ length: ns }, (_, i) => sF[i + 1] - sF[i]), sC = Float64Array.from({ length: ns }, (_, i) => 0.5 * (sF[i] + sF[i + 1])), nC = Float64Array.from({ length: nn }, (_, j) => -D / 2 + (j + 0.5) * dn), nF = Float64Array.from({ length: nn + 1 }, (_, j) => -D / 2 + j * dn);
  const phiOf = (s) => clamp((s - legs) / R, 0, th), tx = (s) => Math.cos(phiOf(s)), ty = (s) => -Math.sin(phiOf(s)); // tangent (tx, ty); normal towards the extrados (−ty, tx)
  const nc = ns * nn, iu = (i, j) => i * nn + j, iv = (i, j) => i * (nn + 1) + j, ic = (i, j) => i * nn + j;
  let us = new Float64Array((ns + 1) * nn), un = new Float64Array(ns * (nn + 1)), rho = new Float64Array(nc).fill(rho0), p = new Float64Array(nc), time = 0;
  const us1 = new Float64Array(us.length), un1 = new Float64Array(un.length), rho1 = new Float64Array(nc), p1 = new Float64Array(nc), A = new Float64Array(nc * (nn + 1)), kap = Float64Array.from(kap0), g = o.gravity || [0, 0], turb = o.turbulence !== false;
  const Uof = typeof o.U === 'function' ? o.U : () => +o.U, slug = o.slug && o.slug.len > 0 ? o.slug : null, rhoIn = (t) => (slug && t >= slug.t0 && (t - slug.t0) * Uof(t) <= slug.len ? rhoS : rho0), pOut = nz(o.pOut, 0);
  us.fill(Uof(0));
  const lim = (a, c) => (a * c <= 0 ? 0 : (2 * a * c) / (a + c)); // van Leer (harmonic) limited slope
  const wallShear = (u, r) => { // wall shear stress from the tangential velocity at the first cell centre
    const y = dn / 2, au = Math.abs(u), nuK = mu / r; if (!(au > 0)) return 0;
    if (!turb) return (mu * u) / y;
    let ut = Math.sqrt((nuK * au) / y); if ((y * ut) / nuK > 11) { for (let k = 0; k < 12; k++) ut = (0.41 * au) / Math.log(Math.max((9 * y * ut) / nuK, 1.5)); } // log law u⁺ = ln(9·y⁺)/0.41
    return Math.sign(u) * r * ut * ut;
  };
  let last = null;
  const evalStep = (dt, motion) => {
    const t1 = time + dt, Uin = Uof(t1), kf = new Float64Array(ns + 1), aS = new Float64Array(ns), aN = new Float64Array(ns), om = new Float64Array(ns);
    for (let i = 0; i < ns; i++) { const m = motion ? motion(sC[i], t1) : null; kap[i] = kap0[i] + (m ? nz(m.dkappa, 0) : 0); aS[i] = m ? nz(m.as, 0) : 0; aN[i] = m ? nz(m.an, 0) : 0; om[i] = m ? nz(m.omega, 0) : 0; }
    for (let i = 0; i <= ns; i++) kf[i] = i === 0 ? kap[0] : i === ns ? kap[ns - 1] : 0.5 * (kap[i - 1] + kap[i]);
    // density transport (conservative, limited upwind)
    for (let i = 0; i < ns; i++) for (let j = 0; j < nn; j++) {
      const fx = (ii) => { const u = us[iu(ii, j)]; if (ii === 0) return u * (u >= 0 ? rhoIn(t1) : rho[ic(0, j)]); if (ii === ns) return u * rho[ic(ns - 1, j)]; const a = u >= 0 ? ii - 1 : ii, s = u >= 0 ? 1 : -1, c = rho[ic(a, j)], up = a - s >= 0 && a - s < ns ? rho[ic(a - s, j)] : c, dnw = a + s >= 0 && a + s < ns ? rho[ic(a + s, j)] : c; return u * (c + 0.5 * lim(c - up, dnw - c)); };
      const fy = (jj) => { if (jj === 0 || jj === nn) return 0; const v = un[iv(i, jj)], a = v >= 0 ? jj - 1 : jj, s = v >= 0 ? 1 : -1, c = rho[ic(i, a)], up = a - s >= 0 && a - s < nn ? rho[ic(i, a - s)] : c, dnw = a + s >= 0 && a + s < nn ? rho[ic(i, a + s)] : c; return (1 + kap[i] * nF[jj]) * v * (c + 0.5 * lim(c - up, dnw - c)); };
      const hc = 1 + kap[i] * nC[j]; rho1[ic(i, j)] = clamp(rho[ic(i, j)] - (dt / (hc * ds[i] * dn)) * ((fx(i + 1) - fx(i)) * dn + (fy(j + 1) - fy(j)) * ds[i]), Math.min(rho0, rhoS), Math.max(rho0, rhoS));
    }
    // explicit momentum: s-faces
    const tauE = new Float64Array(ns), tauI = new Float64Array(ns);
    for (let i = 0; i < ns; i++) { const uE = 0.5 * (us[iu(i, nn - 1)] + us[iu(i + 1, nn - 1)]), uI = 0.5 * (us[iu(i, 0)] + us[iu(i + 1, 0)]); tauE[i] = wallShear(uE, rho1[ic(i, nn - 1)]); tauI[i] = wallShear(uI, rho1[ic(i, 0)]); }
    for (let i = 0; i <= ns; i++) for (let j = 0; j < nn; j++) {
      const k = iu(i, j), u = us[k];
      if (i === 0) { us1[k] = Uin; continue; }
      const il = Math.min(i - 1, ns - 1), ir = Math.min(i, ns - 1), dsf = i === ns ? ds[ns - 1] : 0.5 * (ds[i - 1] + ds[i]), h = 1 + kf[i] * nC[j], r = 0.5 * (rho1[ic(il, j)] + rho1[ic(ir, j)]), nuK = mu / r;
      const uW = us[iu(i - 1, j)], uE = i < ns ? us[iu(i + 1, j)] : u, uS = j > 0 ? us[iu(i, j - 1)] : u, uN = j < nn - 1 ? us[iu(i, j + 1)] : u, v = 0.25 * (un[iv(il, j)] + un[iv(il, j + 1)] + un[iv(ir, j)] + un[iv(ir, j + 1)]);
      const nut = (jf) => { if (!turb || jf <= 0 || jf >= nn) return 0; const y = Math.min(nF[jf] + D / 2, D / 2 - nF[jf]), l = Math.min(0.41 * y, 0.045 * D); return l * l * Math.abs(us[iu(i, jf)] - us[iu(i, jf - 1)]) / dn; };
      let acc = -(u / h) * (u > 0 ? (u - uW) / (i > 1 ? 0.5 * (ds[Math.max(i - 2, 0)] + ds[i - 1]) * 0 + ds[i - 1] : ds[0]) : (uE - u) / ds[ir]) - v * (v > 0 ? (u - uS) / dn : (uN - u) / dn) - (kf[i] * u * v) / h;
      acc += (nuK * (uE - 2 * u + uW)) / (dsf * dsf * h * h) + ((nuK + nut(j + 1)) * (1 + kf[i] * nF[j + 1]) * (uN - u) - (nuK + nut(j)) * (1 + kf[i] * nF[j]) * (u - uS)) / (h * dn * dn);
      if (j === 0) acc -= (0.5 * (tauI[il] + tauI[ir])) / (r * dn); if (j === nn - 1) acc -= (0.5 * (tauE[il] + tauE[ir])) / (r * dn);
      const sm = i === ns ? sF[ns] : sF[i]; acc += g[0] * tx(sm) + g[1] * ty(sm) - 0.5 * (aS[il] + aS[ir]) + (om[il] + om[ir]) * v;
      us1[k] = u + dt * acc;
    }
    for (let i = 0; i < ns; i++) for (let j = 0; j <= nn; j++) {
      const k = iv(i, j); if (j === 0 || j === nn) { un1[k] = 0; continue; }
      const v = un[k], h = 1 + kap[i] * nF[j], r = 0.5 * (rho1[ic(i, j - 1)] + rho1[ic(i, j)]), nuK = mu / r, u = 0.25 * (us[iu(i, j - 1)] + us[iu(i, j)] + us[iu(i + 1, j - 1)] + us[iu(i + 1, j)]);
      const vW = i > 0 ? un[iv(i - 1, j)] : 0, vE = i < ns - 1 ? un[iv(i + 1, j)] : v, vS = un[iv(i, j - 1)], vN = un[iv(i, j + 1)], dsw = i > 0 ? 0.5 * (ds[i - 1] + ds[i]) : ds[0], dse = i < ns - 1 ? 0.5 * (ds[i] + ds[i + 1]) : ds[ns - 1];
      let acc = -(u / h) * (u > 0 ? (v - vW) / dsw : (vE - v) / dse) - v * (v > 0 ? (v - vS) / dn : (vN - v) / dn) + (kap[i] * u * u) / h;
      acc += (nuK * ((vE - v) / dse - (v - vW) / dsw)) / (ds[i] * h * h) + (nuK * (vN - 2 * v + vS)) / (dn * dn);
      acc += -g[0] * ty(sC[i]) + g[1] * tx(sC[i]) - aN[i] - 2 * om[i] * u;
      un1[k] = v + dt * acc;
    }
    // pressure equation with the new density
    A.fill(0); const rhs = new Float64Array(nc), w1 = nn + 1, ai = (r, c) => r * w1 + nn - r + c;
    for (let i = 0; i < ns; i++) for (let j = 0; j < nn; j++) {
      const k = ic(i, j), r = rho1[k]; let diag = 0, q = (us1[iu(i, j)] - us1[iu(i + 1, j)]) * dn + ((1 + kap[i] * nF[j]) * un1[iv(i, j)] - (1 + kap[i] * nF[j + 1]) * un1[iv(i, j + 1)]) * ds[i];
      if (i > 0) { const c = (dt * dn) / (0.5 * (r + rho1[ic(i - 1, j)]) * (1 + kf[i] * nC[j]) * 0.5 * (ds[i - 1] + ds[i])); diag += c; A[ai(k, k - nn)] -= c; }
      if (i < ns - 1) diag += (dt * dn) / (0.5 * (r + rho1[ic(i + 1, j)]) * (1 + kf[i + 1] * nC[j]) * 0.5 * (ds[i] + ds[i + 1])); else { const c = (dt * dn) / (r * (1 + kf[ns] * nC[j]) * 0.5 * ds[i]); diag += c; q += c * pOut; }
      if (j > 0) { const c = (dt * (1 + kap[i] * nF[j]) * ds[i]) / (0.5 * (r + rho1[ic(i, j - 1)]) * dn); diag += c; A[ai(k, k - 1)] -= c; }
      if (j < nn - 1) diag += (dt * (1 + kap[i] * nF[j + 1]) * ds[i]) / (0.5 * (r + rho1[ic(i, j + 1)]) * dn);
      A[ai(k, k)] = diag; rhs[k] = q;
    }
    bandCholesky(A, nc, nn)(rhs); p1.set(rhs);
    for (let i = 1; i <= ns; i++) for (let j = 0; j < nn; j++) {
      if (i < ns) us1[iu(i, j)] -= (dt * (p1[ic(i, j)] - p1[ic(i - 1, j)])) / (0.5 * (rho1[ic(i, j)] + rho1[ic(i - 1, j)]) * (1 + kf[i] * nC[j]) * 0.5 * (ds[i - 1] + ds[i]));
      else us1[iu(i, j)] -= (dt * (pOut - p1[ic(ns - 1, j)])) / (rho1[ic(ns - 1, j)] * (1 + kf[ns] * nC[j]) * 0.5 * ds[ns - 1]);
    }
    for (let i = 0; i < ns; i++) for (let j = 1; j < nn; j++) un1[iv(i, j)] -= (dt * (p1[ic(i, j)] - p1[ic(i, j - 1)])) / (0.5 * (rho1[ic(i, j)] + rho1[ic(i, j - 1)]) * dn);
    // loads on the pipe
    const L = { s: Array.from(sC), pExt: [], pInt: [], tauExt: Array.from(tauE), tauInt: Array.from(tauI), fn: [], ft: [], force: [0, 0], cv: [0, 0], momentum: [0, 0], mass: 0, divergenceMax: 0, t: t1 };
    let mx = 0, my = 0, mass = 0, div = 0;
    for (let i = 0; i < ns; i++) {
      const pe = 1.5 * p1[ic(i, nn - 1)] - 0.5 * p1[ic(i, nn - 2)], pi = 1.5 * p1[ic(i, 0)] - 0.5 * p1[ic(i, 1)], he = 1 + (kap[i] * D) / 2, hi = 1 - (kap[i] * D) / 2, fn = b * (pe * he - pi * hi), ft = b * (tauE[i] * he + tauI[i] * hi), cx = tx(sC[i]), cy = ty(sC[i]);
      L.pExt.push(pe); L.pInt.push(pi); L.fn.push(fn); L.ft.push(ft); L.force[0] += (fn * -cy + ft * cx) * ds[i]; L.force[1] += (fn * cx + ft * cy) * ds[i];
      for (let j = 0; j < nn; j++) { const h = 1 + kap[i] * nC[j], dv = h * ds[i] * dn * b, r = rho1[ic(i, j)], u = 0.5 * (us1[iu(i, j)] + us1[iu(i + 1, j)]), v = 0.5 * (un1[iv(i, j)] + un1[iv(i, j + 1)]); mass += r * dv; mx += r * dv * (u * cx - v * cy); my += r * dv * (u * cy + v * cx);
        div = Math.max(div, Math.abs((us1[iu(i + 1, j)] - us1[iu(i, j)]) * dn + ((1 + kap[i] * nF[j + 1]) * un1[iv(i, j + 1)] - (1 + kap[i] * nF[j]) * un1[iv(i, j)]) * ds[i]) / (h * ds[i] * dn)); }
    }
    let fin = 0, fout = 0, pin = 0;
    for (let j = 0; j < nn; j++) { const pf = 1.5 * p1[ic(0, j)] - 0.5 * p1[ic(1, j)]; pin += pf / nn; fin += (pf + rhoIn(t1) * us1[iu(0, j)] ** 2) * dn * b; fout += (pOut + rho1[ic(ns - 1, j)] * us1[iu(ns, j)] ** 2) * dn * b; }
    const m0 = last ? last.momentum : null, dmx = m0 ? (mx - m0[0]) / dt : 0, dmy = m0 ? (my - m0[1]) / dt : 0;
    L.cv = [fin * tx(0) - fout * tx(len) - dmx + g[0] * mass, fin * ty(0) - fout * ty(len) - dmy + g[1] * mass]; L.momentum = [mx, my]; L.mass = mass; L.divergenceMax = div; L.pIn = pin; L.fluxIn = fin; L.fluxOut = fout;
    return L;
  };
  let pending = null;
  const sim = {
    grid: { ns, nn, s: Array.from(sC), n: Array.from(nC), sFaces: Array.from(sF), ds: Array.from(ds), dn, legs, arc, length: len, R, D, angle: th, depth: b, tangent: (s) => [tx(s), ty(s)],
      toXY: (s, n) => { if (s <= legs) return [s, n]; if (s <= legs + arc) { const f = (s - legs) / R; return [legs + (R + n) * Math.sin(f), -R + (R + n) * Math.cos(f)]; } const a = s - legs - arc; return [legs + (R + n) * Math.sin(th) + a * Math.cos(th), -R + (R + n) * Math.cos(th) - a * Math.sin(th)]; } },
    get time() { return time; },
    dtStable() { let vmax = 1e-9; for (let k = 0; k < us.length; k++) vmax = Math.max(vmax, Math.abs(us[k])); let dmin = dn; for (let i = 0; i < ns; i++) dmin = Math.min(dmin, ds[i] * (1 - D / (2 * R))); return (pos(o.cfl, 0.4) * dmin) / vmax; },
    step(dt, motion) { pending = { L: evalStep(dt, motion), dt }; return pending.L; },
    commit() { if (!pending) return; us.set(us1); un.set(un1); rho.set(rho1); p.set(p1); time += pending.dt; last = pending.L; pending = null; },
    loads() { return last; },
    fields() { const P = [], Rh = [], Us = [], Un = []; for (let j = 0; j < nn; j++) { P.push(Array.from({ length: ns }, (_, i) => p[ic(i, j)])); Rh.push(Array.from({ length: ns }, (_, i) => rho[ic(i, j)])); Us.push(Array.from({ length: ns }, (_, i) => 0.5 * (us[iu(i, j)] + us[iu(i + 1, j)]))); Un.push(Array.from({ length: ns }, (_, i) => 0.5 * (un[iv(i, j)] + un[iv(i, j + 1)]))); } return { p: P, rho: Rh, us: Us, un: Un }; },
  };
  return sim;
}
/**
 * Runs bendCfdSim to a given time (fixed pipe) and returns the history of the wall force with the momentum-theorem values.
 * o: as bendCfdSim plus { tEnd, progress(f) }. Returns { t[], fx[], fy[], f[] (N, wall tractions), cvx[], cvy[] (control volume), pIn[] (Pa), last (loads at the end), fields, grid, steps,
 *   theory: { momentum (ρ·A·v²·√(2(1 − cos θ)) for the carrier fluid), momentumSlug (the same with the slug density) } }.
 */
export function bendCfd(o) {
  const sim = bendCfdSim(o), out = { t: [], fx: [], fy: [], f: [], cvx: [], cvy: [], pIn: [], massErr: 0 }; let steps = 0;
  const tEnd = pos(o.tEnd, (3 * sim.grid.length) / Math.max(typeof o.U === 'function' ? o.U(0) : +o.U, 1e-9));
  while (sim.time < tEnd - 1e-12 && steps < 200000) {
    const dt = Math.min(sim.dtStable(), tEnd - sim.time), L = sim.step(dt, null); sim.commit(); steps++;
    out.t.push(sim.time); out.fx.push(L.force[0]); out.fy.push(L.force[1]); out.f.push(Math.hypot(L.force[0], L.force[1])); out.cvx.push(L.cv[0]); out.cvy.push(L.cv[1]); out.pIn.push(L.pIn);
    if (o.progress && steps % 50 === 0) o.progress(sim.time / tEnd);
  }
  const U0 = typeof o.U === 'function' ? o.U(tEnd) : +o.U, A = (Math.PI * o.D * o.D) / 4, kt = Math.sqrt(2 * (1 - Math.cos(sim.grid.angle)));
  return { ...out, last: sim.loads(), fields: sim.fields(), grid: sim.grid, steps, theory: { momentum: pos(o.rho, 1000) * A * U0 * U0 * kt, momentumSlug: pos(o.rhoSlug, pos(o.rho, 1000)) * A * U0 * U0 * kt, pressure: nz(o.pOut, 0) * A * kt, turning: kt } };
}
// ---- three-dimensional beam (frame) elements and the bend spool coupled to the resolved flow -----------------------------------
/**
 * Space frame of two-node beam elements with six degrees of freedom per node (u, v, w, θx, θy, θz): axial, torsion and Euler–Bernoulli bending about both section axes, consistent mass.
 * o: { nodes: [[x, y, z]], elems: [{ a, b, EA, EIy, EIz, GJ, m (kg/m), rIp (polar radius of gyration², m²) }], fix: [{ node, dofs: [0…5] (default all six) }] }
 * Returns { K, M (dense), ndof, free[], fixed[], nodes, elems, length[] }.
 */
export function frame3d(o) {
  const nodes = o.nodes, N = nodes.length, ndof = 6 * N, K = zeros(ndof), M = zeros(ndof), length = [];
  for (const e of o.elems) {
    const A = nodes[e.a], B = nodes[e.b], dx = B[0] - A[0], dy = B[1] - A[1], dz = B[2] - A[2], l = Math.hypot(dx, dy, dz); need(l > 0, 'A beam element has zero length.'); length.push(l);
    const ex = [dx / l, dy / l, dz / l], ref = Math.abs(ex[2]) < 0.9 ? [0, 0, 1] : [1, 0, 0], ey0 = [ref[1] * ex[2] - ref[2] * ex[1], ref[2] * ex[0] - ref[0] * ex[2], ref[0] * ex[1] - ref[1] * ex[0]], ny = Math.hypot(...ey0), ey = ey0.map((c) => c / ny), ez = [ex[1] * ey[2] - ex[2] * ey[1], ex[2] * ey[0] - ex[0] * ey[2], ex[0] * ey[1] - ex[1] * ey[0]], Rm = [ex, ey, ez];
    const k = zeros(12), m = zeros(12), EA = +e.EA, EIy = +e.EIy, EIz = +e.EIz, GJ = +e.GJ, ml = nz(e.m, 0) * l, set = (A2, i, j, v) => { A2[i][j] += v; if (i !== j) A2[j][i] += v; };
    set(k, 0, 0, EA / l); set(k, 0, 6, -EA / l); set(k, 6, 6, EA / l); set(k, 3, 3, GJ / l); set(k, 3, 9, -GJ / l); set(k, 9, 9, GJ / l);
    set(m, 0, 0, ml / 3); set(m, 0, 6, ml / 6); set(m, 6, 6, ml / 3); const ip = nz(e.rIp, 0) * ml; set(m, 3, 3, ip / 3); set(m, 3, 9, ip / 6); set(m, 9, 9, ip / 3);
    // bending in the local x–y plane (deflection v, rotation θz: indices 1, 5, 7, 11) and in the x–z plane (deflection w, rotation θy with the opposite sign: 2, 4, 8, 10)
    for (const [EI, id, sg] of [[EIz, [1, 5, 7, 11], 1], [EIy, [2, 4, 8, 10], -1]]) {
      const c = EI / l ** 3, kb = [[12, 6 * l, -12, 6 * l], [6 * l, 4 * l * l, -6 * l, 2 * l * l], [-12, -6 * l, 12, -6 * l], [6 * l, 2 * l * l, -6 * l, 4 * l * l]], mb = [[156, 22 * l, 54, -13 * l], [22 * l, 4 * l * l, 13 * l, -3 * l * l], [54, 13 * l, 156, -22 * l], [-13 * l, -3 * l * l, -22 * l, 4 * l * l]], sgn = [1, sg, 1, sg];
      for (let a = 0; a < 4; a++) for (let b2 = 0; b2 < 4; b2++) { k[id[a]][id[b2]] += c * kb[a][b2] * sgn[a] * sgn[b2]; m[id[a]][id[b2]] += (ml / 420) * mb[a][b2] * sgn[a] * sgn[b2]; }
    }
    // rotate to global axes: T = diag(R, R, R, R)
    const T = zeros(12); for (let q = 0; q < 4; q++) for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) T[3 * q + i][3 * q + j] = Rm[i][j];
    const rot = (A2) => { const tmp = zeros(12), out = zeros(12); for (let i = 0; i < 12; i++) for (let j = 0; j < 12; j++) { let s = 0; for (let q = 0; q < 12; q++) s += A2[i][q] * T[q][j]; tmp[i][j] = s; } for (let i = 0; i < 12; i++) for (let j = 0; j < 12; j++) { let s = 0; for (let q = 0; q < 12; q++) s += T[q][i] * tmp[q][j]; out[i][j] = s; } return out; };
    const kg = rot(k), mg = rot(m), map = [...Array(6).keys()].map((d) => 6 * e.a + d).concat([...Array(6).keys()].map((d) => 6 * e.b + d));
    for (let i = 0; i < 12; i++) for (let j = 0; j < 12; j++) { K[map[i]][map[j]] += kg[i][j]; M[map[i]][map[j]] += mg[i][j]; }
  }
  const fx = new Set(); for (const c of o.fix || []) for (const d of c.dofs || [0, 1, 2, 3, 4, 5]) fx.add(6 * c.node + d);
  const free = []; for (let d = 0; d < ndof; d++) if (!fx.has(d)) free.push(d);
  return { K, M, ndof, free, fixed: [...fx].sort((a, b) => a - b), nodes, elems: o.elems, length };
}
/** Static solution of a frame: f (Float64Array of nodal loads, all degrees of freedom). Returns { u, reactions (Float64Array: K·u − f at the fixed degrees of freedom, zero elsewhere), sumReaction: [Fx, Fy, Fz], sumLoad: [Fx, Fy, Fz] }. */
export function frameStatic(fr, f) {
  const L = chol(sub(fr.K, fr.free)); need(L, 'The frame is a mechanism: add restraints.');
  const ur = cholSolve(L, fr.free.map((d) => f[d])), u = new Float64Array(fr.ndof); fr.free.forEach((d, i) => (u[d] = ur[i]));
  const Ku = matVec(fr.K, u), reactions = new Float64Array(fr.ndof), sr = [0, 0, 0], sl = [0, 0, 0];
  for (const d of fr.fixed) { reactions[d] = Ku[d] - f[d]; if (d % 6 < 3) sr[d % 6] += reactions[d]; }
  for (let d = 0; d < fr.ndof; d++) if (d % 6 < 3) sl[d % 6] += f[d];
  return { u, reactions, sumReaction: sr, sumLoad: sl };
}
/** Natural frequencies (Hz) and mode vectors of a frame (lowest p). */
export function frameModes(fr, p = 4) {
  const e = subspaceEig(sub(fr.K, fr.free), sub(fr.M, fr.free), Math.min(p, fr.free.length));
  return { f: e.values.map((v) => Math.sqrt(Math.max(v, 0)) / (2 * Math.PI)), vectors: e.vectors.map((y) => { const u = new Float64Array(fr.ndof); fr.free.forEach((d, i) => (u[d] = y[i])); return u; }) };
}
/** Elbow flexibility and stress-intensification factors in the ASME B31 form: h = t·R/r², k = 1.65/h, i (in-plane) = 0.9/h^⅔, i (out-of-plane) = 0.75/h^⅔ (all at least 1), with von Kármán's first approximation for comparison. */
export function elbowFactors({ t, R, r }) {
  const h = (t * R) / (r * r);
  return { h, k: Math.max(1.65 / h, 1), kKarman: (12 * h * h + 10) / (12 * h * h + 1), iIn: Math.max(0.9 / h ** (2 / 3), 1), iOut: Math.max(0.75 / h ** (2 / 3), 1) };
}
/**
 * A bend spool (two straight legs and an arc, in the x–y plane, both ends anchored) loaded by the resolved flow of bendCfdSim, one-way or two-way.
 * One-way: the wall pressure and shear of the flow on the fixed geometry load the beam model (whose mass includes the contents). Two-way: the beam carries the steel only and its motion is returned
 * to the flow — acceleration, rotation rate and change of curvature of every section enter the fluid equations on the moving grid — so that the inertia of the contents, the Coriolis and the
 * centrifugal forces follow from the flow solution; the two are iterated in every time step with Aitken's relaxation.
 * o: bendCfdSim options plus { t (wall, m), E, nu, rhoSteel, mode: 'one' | 'two', nEl (beam elements, default 24), dt, tEnd, zeta (damping ratio at the first mode), flexibility (apply the elbow factor to the arc; default true),
 *      d0: { mode: k, amplitude } (start from a multiple of a mode shape, liquid at rest), tol, maxIter, relax0, progress }
 * Returns { t[], disp[] (resultant displacement at the bend mid-point, m), ux[], uy[], fFluid[] (N, resultant on the spool), reaction[] (N, resultant of the anchor forces), moment[] (N·m, largest bending moment),
 *   iterations: { mean, max }, converged, energy: { work (of the fluid on the structure, J), gain (strain + kinetic + damping loss, J) }, equilibrium: { fluid: [Fx, Fy], reaction: [Fx, Fy], inertia: [Fx, Fy] } (last step),
 *   modes: { empty[], filled[] (Hz) }, frame, sim, steps }.
 */
export function bendFsi(o) {
  const sim = bendCfdSim(o), g = sim.grid, D = +o.D, t = +o.t, ro = D / 2 + t, E = +o.E, nu = nz(o.nu, 0.3), As = Math.PI * (ro * ro - (D * D) / 4), I = (Math.PI / 4) * (ro ** 4 - (D / 2) ** 4), ms = nz(o.rhoSteel, 7850) * As, mf = pos(o.rho, 1000) * (Math.PI * D * D) / 4, two = o.mode === 'two';
  const nEl = clamp(Math.round(nz(o.nEl, 24)), 6, 120), stride = Math.max(1, Math.round(g.ns / nEl)), idx = []; for (let i = 0; i < g.ns; i += stride) idx.push(i); if (g.ns - idx[idx.length - 1] < Math.max(1, stride / 2)) idx.pop(); idx.push(g.ns);
  const sN = idx.map((i) => g.sFaces[i]), nodes = sN.map((s) => { const q = g.toXY(s, 0); return [q[0], q[1], 0]; }), kEl = o.flexibility === false ? 1 : elbowFactors({ t, R: g.R, r: ro - t / 2 }).k;
  const build = (m) => frame3d({ nodes, elems: sN.slice(1).map((s, e) => { const mid = 0.5 * (s + sN[e]), curved = mid > g.legs && mid < g.legs + g.arc; return { a: e, b: e + 1, EA: E * As, EIy: (E * I) / (curved ? kEl : 1), EIz: (E * I) / (curved ? kEl : 1), GJ: (E * 2 * I) / (2 * (1 + nu)), m, rIp: (2 * I) / As }; }), fix: [{ node: 0 }, { node: nodes.length - 1 }, ...nodes.map((_, i) => ({ node: i, dofs: [2, 3, 4] }))] }); // in-plane motion only: the flow model is plane
  const frS = build(ms), frF = build(ms + mf), fr = two ? frS : frF, fd = fr.free, nf = fd.length, Kf = sub(fr.K, fd), Mf = sub(fr.M, fd), mE = frameModes(frS, 3), mFl = frameModes(frF, 3);
  const om1 = 2 * Math.PI * (two ? mFl.f[0] : mFl.f[0]), zeta = Math.max(nz(o.zeta, 0), 0), Cf = Mf.map((r) => Float64Array.from(r, (v) => 2 * zeta * om1 * v)), dt = +o.dt, steps = Math.max(1, Math.round(o.tEnd / dt)); need(dt > 0, 'The coupled bend model needs a time step.');
  // Newmark (average acceleration) effective matrix
  const Keff = zeros(nf); for (let i = 0; i < nf; i++) for (let j = 0; j < nf; j++) Keff[i][j] = Kf[i][j] + (4 / (dt * dt)) * Mf[i][j] + (2 / dt) * Cf[i][j];
  const Lk = chol(Keff), Lm = chol(Mf); need(Lk && Lm, 'The spool model is singular.');
  let d = new Float64Array(nf), v = new Float64Array(nf), a = new Float64Array(nf);
  if (o.d0) { const md = (two ? mFl : mFl).vectors[o.d0.mode || 0]; let big = 0; for (const q of fd) if (q % 6 < 3) big = Math.max(big, Math.abs(md[q])); fd.forEach((q, i) => (d[i] = (md[q] / big) * o.d0.amplitude)); }
  // section kinematics from a state of the free degrees of freedom: linear interpolation of the nodal values along s
  const full = (x) => { const u = new Float64Array(fr.ndof); fd.forEach((q, i) => (u[q] = x[i])); return u; };
  const motionOf = (dn, vn, an) => { const D2 = full(dn), V2 = full(vn), A2 = full(an); return (s) => { let e = 0; while (e < sN.length - 2 && s > sN[e + 1]) e++; const l = sN[e + 1] - sN[e], w = clamp((s - sN[e]) / l, 0, 1), at = (U2, c) => (1 - w) * U2[6 * e + c] + w * U2[6 * e + 6 + c], tg = g.tangent(s), ax = at(A2, 0), ay = at(A2, 1); return { as: ax * tg[0] + ay * tg[1], an: -ax * tg[1] + ay * tg[0], omega: at(V2, 5), dkappa: -(D2[6 * e + 11] - D2[6 * e + 5]) / l }; }; };
  const loadOf = (L) => { const f = new Float64Array(fr.ndof); for (let i = 0; i < g.ns; i++) { const s = g.s[i], tg = g.tangent(s), fx = (L.fn[i] * -tg[1] + L.ft[i] * tg[0]) * g.ds[i], fy = (L.fn[i] * tg[0] + L.ft[i] * tg[1]) * g.ds[i]; let e = 0; while (e < sN.length - 2 && s > sN[e + 1]) e++; const w = clamp((s - sN[e]) / (sN[e + 1] - sN[e]), 0, 1); f[6 * e] += (1 - w) * fx; f[6 * e + 1] += (1 - w) * fy; f[6 * e + 6] += w * fx; f[6 * e + 7] += w * fy; } return f; };
  const newmarkStep = (fFree) => { const rhs = new Float64Array(nf), tm = new Float64Array(nf), tc = new Float64Array(nf); for (let i = 0; i < nf; i++) { tm[i] = (4 / (dt * dt)) * d[i] + (4 / dt) * v[i] + a[i]; tc[i] = (2 / dt) * d[i] + v[i]; } const Mt = matVec(Mf, tm), Ct = matVec(Cf, tc); for (let i = 0; i < nf; i++) rhs[i] = fFree[i] + Mt[i] + Ct[i]; const dn = cholSolve(Lk, rhs), an = new Float64Array(nf), vn = new Float64Array(nf); for (let i = 0; i < nf; i++) { an[i] = (4 / (dt * dt)) * (dn[i] - d[i]) - (4 / dt) * v[i] - a[i]; vn[i] = v[i] + 0.5 * dt * (a[i] + an[i]); } return { dn, vn, an }; };
  const energy = (dd, vv) => 0.5 * dot(dd, matVec(Kf, dd)) + 0.5 * dot(vv, matVec(Mf, vv));
  // initial acceleration consistent with the initial displacement (and, two-way, with the liquid at rest: effective mass of the filled spool for the first guess)
  { const Kd = matVec(Kf, d); a = cholSolve(two ? chol(sub(frF.M, fd)) : Lm, Float64Array.from(Kd, (x) => -x)); }
  const mid = (() => { const sm = g.legs + g.arc / 2; let k = 0; sN.forEach((s, i) => { if (Math.abs(s - sm) < Math.abs(sN[k] - sm)) k = i; }); return k; })();
  const out = { t: [0], disp: [0], ux: [0], uy: [0], fFluid: [0], reaction: [0], moment: [0] }; { const u0 = full(d); out.ux[0] = u0[6 * mid]; out.uy[0] = u0[6 * mid + 1]; out.disp[0] = Math.hypot(out.ux[0], out.uy[0]); }
  let itTot = 0, itMax = 0, allConv = true, work = 0, damp = 0, omega = pos(o.relax0, 0.5), fPrev = null, eq = null; const e0 = energy(d, v), tol = pos(o.tol, 1e-7), maxIter = Math.round(pos(o.maxIter, 40));
  for (let s = 1; s <= steps; s++) {
    let st = null, L = null, fFull = null;
    if (!two) { L = sim.step(dt, null); fFull = loadOf(L); st = newmarkStep(fd.map((q) => fFull[q])); itTot++; itMax = Math.max(itMax, 1); }
    else {
      let guess = { dn: Float64Array.from(d, (x, i) => x + dt * v[i] + 0.5 * dt * dt * a[i]), vn: Float64Array.from(v, (x, i) => x + dt * a[i]), an: Float64Array.from(a) }, rOld = null, it = 0, ok = false;
      for (; it < maxIter; ) {
        L = sim.step(dt, motionOf(guess.dn, guess.vn, guess.an)); fFull = loadOf(L); st = newmarkStep(fd.map((q) => fFull[q])); it++;
        const r = Float64Array.from(st.dn, (x, i) => x - guess.dn[i]), rn = Math.sqrt(dot(r, r)), sn = Math.sqrt(dot(st.dn, st.dn));
        if (rn <= tol * Math.max(sn, 1e-12)) { ok = true; break; }
        if (rOld) { let num = 0, den = 0; for (let i = 0; i < nf; i++) { const dr = r[i] - rOld[i]; num += rOld[i] * dr; den += dr * dr; } if (den > 0) { const om = (-omega * num) / den; if (Number.isFinite(om) && om !== 0) omega = Math.sign(om) * clamp(Math.abs(om), 0.02, 1.5); } }
        rOld = r; const dn2 = Float64Array.from(guess.dn, (x, i) => x + omega * r[i]); guess = { dn: dn2, an: Float64Array.from(dn2, (x, i) => (4 / (dt * dt)) * (x - d[i]) - (4 / dt) * v[i] - a[i]), vn: null }; guess.vn = Float64Array.from(v, (x, i) => x + 0.5 * dt * (a[i] + guess.an[i]));
      }
      if (!ok) allConv = false; itTot += it; itMax = Math.max(itMax, it);
    }
    sim.commit();
    const fF = fd.map((q) => fFull[q]), fm = fPrev ? fF.map((x, i) => 0.5 * (x + fPrev[i])) : fF, vm = Float64Array.from(v, (x, i) => 0.5 * (x + st.vn[i])); for (let i = 0; i < nf; i++) work += fm[i] * (st.dn[i] - d[i]); damp += dt * dot(vm, matVec(Cf, vm)); fPrev = fF;
    d = st.dn; v = st.vn; a = st.an;
    const u = full(d), Ku = matVec(fr.K, u), Ma = matVec(fr.M, full(a)), Cv = matVec(fr.M, full(v)), rx = [0, 0], inx = [0, 0]; let mom = 0;
    for (const q of fr.fixed) if (q % 6 < 2) rx[q % 6] += Ku[q] + Ma[q] + 2 * zeta * om1 * Cv[q] - fFull[q];
    for (let q = 0; q < fr.ndof; q++) if (q % 6 < 2) inx[q % 6] += Ma[q] + 2 * zeta * om1 * Cv[q];
    for (let e = 0; e < nodes.length - 1; e++) { const l = fr.length[e], ei = fr.elems[e].EIz; mom = Math.max(mom, Math.abs((ei * (u[6 * e + 11] - u[6 * e + 5])) / l)); }
    eq = { fluid: [L.force[0], L.force[1]], reaction: rx, inertia: inx, applied: [0, 1].map((c) => { let q2 = 0; for (let n2 = 0; n2 < nodes.length; n2++) q2 += fFull[6 * n2 + c]; return q2; }) };
    out.t.push(s * dt); out.ux.push(u[6 * mid]); out.uy.push(u[6 * mid + 1]); out.disp.push(Math.hypot(u[6 * mid], u[6 * mid + 1])); out.fFluid.push(Math.hypot(L.force[0], L.force[1])); out.reaction.push(Math.hypot(rx[0], rx[1])); out.moment.push(mom);
    if (o.progress && s % 20 === 0) o.progress(s / steps);
  }
  return { ...out, iterations: { mean: itTot / steps, max: itMax }, converged: allConv, energy: { work, gain: energy(d, v) - e0 + damp, damping: damp }, equilibrium: eq, modes: { empty: mE.f, filled: mFl.f }, frame: fr, grid: g, last: sim.loads(), fields: sim.fields(), steps, mode: two ? 'two' : 'one', flexibilityFactor: kEl, massRatio: mf / ms };
}
// ---- three-dimensional solid finite elements: eight-node hexahedra with incompatible modes, conjugate gradients ---------------
const HEX_S = [[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1], [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]], HEX_G = 1 / Math.sqrt(3);
/** Inverse of a small dense matrix by Gauss–Jordan elimination with partial pivoting (rows are arrays). */
function invSmall(A) {
  const n = A.length, a = A.map((r, i) => { const row = new Float64Array(2 * n); for (let j = 0; j < n; j++) row[j] = r[j]; row[n + i] = 1; return row; });
  for (let c = 0; c < n; c++) { let p = c; for (let i = c + 1; i < n; i++) if (Math.abs(a[i][c]) > Math.abs(a[p][c])) p = i; need(Math.abs(a[p][c]) > 1e-300, 'A solid element is degenerate.'); [a[c], a[p]] = [a[p], a[c]]; const d = a[c][c]; for (let j = 0; j < 2 * n; j++) a[c][j] /= d; for (let i = 0; i < n; i++) if (i !== c) { const f = a[i][c]; if (f !== 0) for (let j = 0; j < 2 * n; j++) a[i][j] -= f * a[c][j]; } }
  return a.map((r) => r.slice(n));
}
/**
 * Linear-elastic solution of a mesh of eight-node hexahedra. With `incompatible` (default) every element carries the nine internal bending modes of Wilson and Taylor
 * (1 − ξ², 1 − η², 1 − ζ² for each displacement component, their gradients taken with the Jacobian of the element centre so that the patch test is passed), condensed out on the
 * element level: one or two elements through a pipe wall then reproduce shell bending. The equations are solved by conjugate gradients preconditioned with one symmetric
 * block Gauss–Seidel sweep over the 3 × 3 node blocks.
 * o: { X (Float64Array, 3 per node), elems: [[8 nodes]], E, nu, fixed (Uint8Array, 3 per node: 1 = zero displacement), f (Float64Array of nodal forces), tol (relative residual, default 1e-8), maxIter, incompatible, progress(f) }
 * Returns { u, stress (per node: [xx, yy, zz, xy, yz, zx], averaged over the adjacent elements), vm (per node), reactions (Float64Array, non-zero at fixed degrees of freedom), energy (½ uᵀK u), work (½ fᵀu),
 *   iterations, residual (relative), converged, ndof, nnz (stored 3 × 3 blocks) }.
 */
export function hexSolve(o) {
  const X = o.X, nN = X.length / 3, ndof = 3 * nN, E = +o.E, nu = +o.nu, inc = o.incompatible !== false, fixed = o.fixed || new Uint8Array(ndof), f = o.f || new Float64Array(ndof);
  need(E > 0 && nu > -1 && nu < 0.5, 'The solid model needs a positive modulus and a Poisson ratio below 0.5.');
  const lam = (E * nu) / ((1 + nu) * (1 - 2 * nu)), G2 = E / (2 * (1 + nu)), C = [[lam + 2 * G2, lam, lam], [lam, lam + 2 * G2, lam], [lam, lam, lam + 2 * G2]];
  // node adjacency → block compressed rows
  const adj = Array.from({ length: nN }, () => new Set()); for (const e of o.elems) for (const a of e) for (const b of e) adj[a].add(b);
  const rowPtr = new Int32Array(nN + 1), cols = []; for (let i = 0; i < nN; i++) { const c = [...adj[i]].sort((p, q) => p - q); rowPtr[i + 1] = rowPtr[i] + c.length; cols.push(c); }
  const col = new Int32Array(rowPtr[nN]), diagPos = new Int32Array(nN); for (let i = 0; i < nN; i++) cols[i].forEach((c, k) => { col[rowPtr[i] + k] = c; if (c === i) diagPos[i] = rowPtr[i] + k; });
  const val = new Float64Array(9 * rowPtr[nN]), find = (i, j) => { let lo = rowPtr[i], hi = rowPtr[i + 1] - 1; while (lo < hi) { const m = (lo + hi) >> 1; if (col[m] < j) lo = m + 1; else hi = m; } return lo; };
  // element kernel: fills B (6 × 24) and Ba (6 × 9) at a natural point, returns the Jacobian determinant (absolute value)
  const B = Array.from({ length: 6 }, () => new Float64Array(24)), Ba = Array.from({ length: 6 }, () => new Float64Array(9)), xe = new Float64Array(24), J = [new Float64Array(3), new Float64Array(3), new Float64Array(3)], Ji = [new Float64Array(3), new Float64Array(3), new Float64Array(3)], J0i = [new Float64Array(3), new Float64Array(3), new Float64Array(3)], dN = new Float64Array(24);
  const jac = (xi, et, ze, out) => {
    for (let a = 0; a < 8; a++) { const s = HEX_S[a]; dN[3 * a] = 0.125 * s[0] * (1 + s[1] * et) * (1 + s[2] * ze); dN[3 * a + 1] = 0.125 * s[1] * (1 + s[0] * xi) * (1 + s[2] * ze); dN[3 * a + 2] = 0.125 * s[2] * (1 + s[0] * xi) * (1 + s[1] * et); }
    for (let r = 0; r < 3; r++) for (let c = 0; c < 3; c++) { let s = 0; for (let a = 0; a < 8; a++) s += dN[3 * a + r] * xe[3 * a + c]; J[r][c] = s; } // J[r][c] = ∂x_c/∂ξ_r
    const d = J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1]) - J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0]) + J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);
    need(Math.abs(d) > 1e-300, 'A solid element is degenerate.');
    out[0][0] = (J[1][1] * J[2][2] - J[1][2] * J[2][1]) / d; out[0][1] = (J[0][2] * J[2][1] - J[0][1] * J[2][2]) / d; out[0][2] = (J[0][1] * J[1][2] - J[0][2] * J[1][1]) / d;
    out[1][0] = (J[1][2] * J[2][0] - J[1][0] * J[2][2]) / d; out[1][1] = (J[0][0] * J[2][2] - J[0][2] * J[2][0]) / d; out[1][2] = (J[0][2] * J[1][0] - J[0][0] * J[1][2]) / d;
    out[2][0] = (J[1][0] * J[2][1] - J[1][1] * J[2][0]) / d; out[2][1] = (J[0][1] * J[2][0] - J[0][0] * J[2][1]) / d; out[2][2] = (J[0][0] * J[1][1] - J[0][1] * J[1][0]) / d; // out[c][r] = ∂ξ_r/∂x_c
    return d;
  };
  const setB = (M, k, gx, gy, gz) => { M[0][3 * k] = gx; M[1][3 * k + 1] = gy; M[2][3 * k + 2] = gz; M[3][3 * k] = gy; M[3][3 * k + 1] = gx; M[4][3 * k + 1] = gz; M[4][3 * k + 2] = gy; M[5][3 * k] = gz; M[5][3 * k + 2] = gx; };
  let det0 = 1;
  const kernel = (xi, et, ze) => {
    const d = jac(xi, et, ze, Ji);
    for (let a = 0; a < 8; a++) setB(B, a, Ji[0][0] * dN[3 * a] + Ji[0][1] * dN[3 * a + 1] + Ji[0][2] * dN[3 * a + 2], Ji[1][0] * dN[3 * a] + Ji[1][1] * dN[3 * a + 1] + Ji[1][2] * dN[3 * a + 2], Ji[2][0] * dN[3 * a] + Ji[2][1] * dN[3 * a + 1] + Ji[2][2] * dN[3 * a + 2]);
    if (inc) { const sc = det0 / d, nat = [-2 * xi, -2 * et, -2 * ze]; for (let m = 0; m < 3; m++) setB(Ba, m, sc * J0i[0][m] * nat[m], sc * J0i[1][m] * nat[m], sc * J0i[2][m] * nat[m]); }
    return Math.abs(d);
  };
  const DB = Array.from({ length: 6 }, () => new Float64Array(24)), DBa = Array.from({ length: 6 }, () => new Float64Array(9)), Ke = Array.from({ length: 24 }, () => new Float64Array(24)), Kua = Array.from({ length: 24 }, () => new Float64Array(9)), Kaa = Array.from({ length: 9 }, () => new Float64Array(9));
  const applyD = (M, out, n) => { for (let c = 0; c < n; c++) { const a = M[0][c], b = M[1][c], d = M[2][c]; out[0][c] = C[0][0] * a + lam * b + lam * d; out[1][c] = lam * a + C[1][1] * b + lam * d; out[2][c] = lam * a + lam * b + C[2][2] * d; out[3][c] = G2 * M[3][c]; out[4][c] = G2 * M[4][c]; out[5][c] = G2 * M[5][c]; } };
  const elementK = (e) => {
    for (let a = 0; a < 8; a++) { xe[3 * a] = X[3 * e[a]]; xe[3 * a + 1] = X[3 * e[a] + 1]; xe[3 * a + 2] = X[3 * e[a] + 2]; }
    for (const r of Ke) r.fill(0); for (const r of Kua) r.fill(0); for (const r of Kaa) r.fill(0);
    if (inc) det0 = jac(0, 0, 0, J0i);
    for (const s of HEX_S) {
      const w = kernel(s[0] * HEX_G, s[1] * HEX_G, s[2] * HEX_G); applyD(B, DB, 24);
      for (let i = 0; i < 24; i++) for (let j = i; j < 24; j++) { let q = 0; for (let r = 0; r < 6; r++) q += B[r][i] * DB[r][j]; Ke[i][j] += w * q; }
      if (inc) { applyD(Ba, DBa, 9); for (let i = 0; i < 24; i++) for (let j = 0; j < 9; j++) { let q = 0; for (let r = 0; r < 6; r++) q += B[r][i] * DBa[r][j]; Kua[i][j] += w * q; } for (let i = 0; i < 9; i++) for (let j = 0; j < 9; j++) { let q = 0; for (let r = 0; r < 6; r++) q += Ba[r][i] * DBa[r][j]; Kaa[i][j] += w * q; } }
    }
    for (let i = 0; i < 24; i++) for (let j = 0; j < i; j++) Ke[i][j] = Ke[j][i];
    if (!inc) return null;
    const Ki = invSmall(Kaa), G = Array.from({ length: 9 }, () => new Float64Array(24)); // G = Kaa⁻¹·Kuaᵀ; the internal modes are α = −G·u
    for (let a = 0; a < 9; a++) for (let j = 0; j < 24; j++) { let q = 0; for (let b = 0; b < 9; b++) q += Ki[a][b] * Kua[j][b]; G[a][j] = q; }
    for (let i = 0; i < 24; i++) for (let j = 0; j < 24; j++) { let q = 0; for (let a = 0; a < 9; a++) q += Kua[i][a] * G[a][j]; Ke[i][j] -= q; }
    return G;
  };
  let ne = 0;
  for (const e of o.elems) {
    elementK(e);
    for (let a = 0; a < 8; a++) for (let b = 0; b < 8; b++) { const pos9 = 9 * find(e[a], e[b]); for (let r = 0; r < 3; r++) for (let c = 0; c < 3; c++) val[pos9 + 3 * r + c] += Ke[3 * a + r][3 * b + c]; }
    if (o.progress && ++ne % 500 === 0) o.progress((0.3 * ne) / o.elems.length);
  }
  const mul = (x, y, skipFixed) => { for (let i = 0; i < nN; i++) { let y0 = 0, y1 = 0, y2 = 0; for (let k = rowPtr[i]; k < rowPtr[i + 1]; k++) { const j = 3 * col[k], v = 9 * k, a = x[j], b = x[j + 1], c = x[j + 2]; y0 += val[v] * a + val[v + 1] * b + val[v + 2] * c; y1 += val[v + 3] * a + val[v + 4] * b + val[v + 5] * c; y2 += val[v + 6] * a + val[v + 7] * b + val[v + 8] * c; } y[3 * i] = y0; y[3 * i + 1] = y1; y[3 * i + 2] = y2; } if (skipFixed) for (let d = 0; d < ndof; d++) if (fixed[d]) y[d] = x[d]; return y; };
  // constrained operator: fixed entries of the iterate stay zero, so the product is taken with them masked
  const Dinv = new Float64Array(9 * nN);
  for (let i = 0; i < nN; i++) { const v = 9 * diagPos[i], m = [[val[v], val[v + 1], val[v + 2]], [val[v + 3], val[v + 4], val[v + 5]], [val[v + 6], val[v + 7], val[v + 8]]]; for (let r = 0; r < 3; r++) if (fixed[3 * i + r]) { for (let c = 0; c < 3; c++) { m[r][c] = 0; m[c][r] = 0; } m[r][r] = 1; } const q = invSmall(m); for (let r = 0; r < 3; r++) for (let c = 0; c < 3; c++) Dinv[9 * i + 3 * r + c] = q[r][c]; }
  const mask = (x) => { for (let d = 0; d < ndof; d++) if (fixed[d]) x[d] = 0; return x; };
  const tmp = new Float64Array(ndof), precond = (r, z) => { // one symmetric block Gauss–Seidel sweep (fixed degrees of freedom masked)
    for (let i = 0; i < nN; i++) { let s0 = r[3 * i], s1 = r[3 * i + 1], s2 = r[3 * i + 2]; for (let k = rowPtr[i]; k < diagPos[i]; k++) { const j = 3 * col[k], v = 9 * k, a = tmp[j], b = tmp[j + 1], c = tmp[j + 2]; s0 -= val[v] * a + val[v + 1] * b + val[v + 2] * c; s1 -= val[v + 3] * a + val[v + 4] * b + val[v + 5] * c; s2 -= val[v + 6] * a + val[v + 7] * b + val[v + 8] * c; }
      if (fixed[3 * i]) s0 = 0; if (fixed[3 * i + 1]) s1 = 0; if (fixed[3 * i + 2]) s2 = 0; const q = 9 * i; tmp[3 * i] = Dinv[q] * s0 + Dinv[q + 1] * s1 + Dinv[q + 2] * s2; tmp[3 * i + 1] = Dinv[q + 3] * s0 + Dinv[q + 4] * s1 + Dinv[q + 5] * s2; tmp[3 * i + 2] = Dinv[q + 6] * s0 + Dinv[q + 7] * s1 + Dinv[q + 8] * s2; if (fixed[3 * i]) tmp[3 * i] = 0; if (fixed[3 * i + 1]) tmp[3 * i + 1] = 0; if (fixed[3 * i + 2]) tmp[3 * i + 2] = 0; }
    for (let i = nN - 1; i >= 0; i--) { let s0 = 0, s1 = 0, s2 = 0; for (let k = diagPos[i] + 1; k < rowPtr[i + 1]; k++) { const j = 3 * col[k], v = 9 * k, a = z[j], b = z[j + 1], c = z[j + 2]; s0 += val[v] * a + val[v + 1] * b + val[v + 2] * c; s1 += val[v + 3] * a + val[v + 4] * b + val[v + 5] * c; s2 += val[v + 6] * a + val[v + 7] * b + val[v + 8] * c; }
      if (fixed[3 * i]) s0 = 0; if (fixed[3 * i + 1]) s1 = 0; if (fixed[3 * i + 2]) s2 = 0; const q = 9 * i; z[3 * i] = tmp[3 * i] - (Dinv[q] * s0 + Dinv[q + 1] * s1 + Dinv[q + 2] * s2); z[3 * i + 1] = tmp[3 * i + 1] - (Dinv[q + 3] * s0 + Dinv[q + 4] * s1 + Dinv[q + 5] * s2); z[3 * i + 2] = tmp[3 * i + 2] - (Dinv[q + 6] * s0 + Dinv[q + 7] * s1 + Dinv[q + 8] * s2); if (fixed[3 * i]) z[3 * i] = 0; if (fixed[3 * i + 1]) z[3 * i + 1] = 0; if (fixed[3 * i + 2]) z[3 * i + 2] = 0; }
    return z;
  };
  const u = new Float64Array(ndof), r = mask(Float64Array.from(f)), z = new Float64Array(ndof), p = new Float64Array(ndof), Ap = new Float64Array(ndof), tol = pos(o.tol, 1e-8), maxIter = Math.round(pos(o.maxIter, 20 * Math.sqrt(ndof) + 2000)), r0 = Math.sqrt(dot(r, r));
  let it = 0, rn = r0, rz = 0;
  if (r0 > 0) { precond(r, z); p.set(z); rz = dot(r, z); for (; it < maxIter && rn > tol * r0; it++) { mask(mul(p, Ap)); const al = rz / dot(p, Ap); for (let d = 0; d < ndof; d++) { u[d] += al * p[d]; r[d] -= al * Ap[d]; } rn = Math.sqrt(dot(r, r)); precond(r, z); const rz1 = dot(r, z), be = rz1 / rz; rz = rz1; for (let d = 0; d < ndof; d++) p[d] = z[d] + be * p[d]; if (o.progress && it % 100 === 0) o.progress(0.3 + 0.6 * Math.min(1, Math.log(r0 / rn) / Math.log(1 / tol))); } }
  const Ku = mul(u, new Float64Array(ndof)), reactions = new Float64Array(ndof); for (let d = 0; d < ndof; d++) if (fixed[d]) reactions[d] = Ku[d] - f[d];
  // stresses at the element nodes (with the internal modes), averaged at the nodes
  const stress = Array.from({ length: nN }, () => new Float64Array(6)), cnt = new Uint16Array(nN), ue = new Float64Array(24), al = new Float64Array(9), eps = new Float64Array(6);
  for (const e of o.elems) {
    const Gm = elementK(e); for (let a = 0; a < 8; a++) { ue[3 * a] = u[3 * e[a]]; ue[3 * a + 1] = u[3 * e[a] + 1]; ue[3 * a + 2] = u[3 * e[a] + 2]; }
    if (Gm) for (let a = 0; a < 9; a++) { let q = 0; for (let j = 0; j < 24; j++) q += Gm[a][j] * ue[j]; al[a] = -q; }
    for (let a = 0; a < 8; a++) { const s = HEX_S[a]; kernel(s[0], s[1], s[2]);
      for (let c = 0; c < 6; c++) { let q = 0; for (let j = 0; j < 24; j++) q += B[c][j] * ue[j]; if (Gm) for (let j = 0; j < 9; j++) q += Ba[c][j] * al[j]; eps[c] = q; }
      const st = stress[e[a]], tr = lam * (eps[0] + eps[1] + eps[2]); st[0] += tr + 2 * G2 * eps[0]; st[1] += tr + 2 * G2 * eps[1]; st[2] += tr + 2 * G2 * eps[2]; st[3] += G2 * eps[3]; st[4] += G2 * eps[4]; st[5] += G2 * eps[5]; cnt[e[a]]++; }
  }
  const vm = new Float64Array(nN); for (let n = 0; n < nN; n++) { const s = stress[n], c = cnt[n] || 1; for (let k = 0; k < 6; k++) s[k] /= c; vm[n] = Math.sqrt(0.5 * ((s[0] - s[1]) ** 2 + (s[1] - s[2]) ** 2 + (s[2] - s[0]) ** 2) + 3 * (s[3] * s[3] + s[4] * s[4] + s[5] * s[5])); }
  return { u, stress, vm, reactions, energy: 0.5 * dot(u, Ku), work: 0.5 * dot(f, u), iterations: it, residual: r0 > 0 ? rn / r0 : 0, converged: !(rn > tol * r0), ndof, nnz: rowPtr[nN] };
}
/**
 * Structured hexahedral mesh of a pipe wall along a centreline in the x–y plane: straight leg, circular arc (bend), straight leg — or a straight pipe when angle = 0.
 * Half model (0 ≤ φ ≤ π, symmetric about the plane of the bend) or full circumference. Local frame of a section: tangent t, n towards the extrados, z out of the plane; a wall point is c(s) + r·(cos φ·n + sin φ·z).
 * o: { ro, t, R (centreline bend radius), angle (deg), legs: [m, m], nc (elements around the modelled circumference), nt (through the wall), nArc, nLeg, half (default true) }
 * Returns { X, elems, id(i, j, k), ns, ncN, nt, s[] (station of every node ring), frame(s) -> { c, t, n }, innerFaces, endFaces(0 | 1) (quads as node lists), ring(i, k) (nodes of a ring), half, length, legs, arc, R, ro, t }.
 */
export function pipeHexMesh(o) {
  const ro = +o.ro, t = +o.t, ri = ro - t, th = (clamp(nz(o.angle, 0), 0, 180) * Math.PI) / 180, R = th > 0 ? Math.max(+o.R, 1.01 * ro) : Infinity, arc = th > 0 ? R * th : 0, legs = o.legs || [2 * ro, 2 * ro], half = o.half !== false;
  const nc = Math.max(half ? 4 : 8, Math.round(nz(o.nc, half ? 12 : 24))), nt = Math.max(1, Math.round(nz(o.nt, 2))), nA = th > 0 ? Math.max(2, Math.round(nz(o.nArc, 12))) : 0, nL = [0, 1].map((k) => (legs[k] > 0 ? Math.max(1, Math.round(nz(Array.isArray(o.nLeg) ? o.nLeg[k] : o.nLeg, 6))) : 0));
  const s = [0]; for (let i = 1; i <= nL[0]; i++) s.push((legs[0] * i) / nL[0]); for (let i = 1; i <= nA; i++) s.push(legs[0] + (arc * i) / nA); for (let i = 1; i <= nL[1]; i++) s.push(legs[0] + arc + (legs[1] * i) / nL[1]);
  const ns = s.length - 1, ncN = half ? nc + 1 : nc, id = (i, j, k) => (i * ncN + (half ? j : ((j % nc) + nc) % nc)) * (nt + 1) + k, X = new Float64Array(3 * (ns + 1) * ncN * (nt + 1));
  const frame = (sv) => { if (!(th > 0) || sv <= legs[0]) return { c: [sv, 0, 0], t: [1, 0, 0], n: [0, 1, 0] }; const f = Math.min((sv - legs[0]) / R, th), cf = Math.cos(f), sf = Math.sin(f), a = Math.max(sv - legs[0] - arc, 0); return { c: [legs[0] + R * sf + a * cf, -R + R * cf - a * sf, 0], t: [cf, -sf, 0], n: [sf, cf, 0] }; };
  for (let i = 0; i <= ns; i++) { const fr = frame(s[i]); for (let j = 0; j < ncN; j++) { const ph = ((half ? Math.PI : 2 * Math.PI) * j) / nc, cp = Math.cos(ph), sp = Math.sin(ph); for (let k = 0; k <= nt; k++) { const r = ri + (t * k) / nt, q = 3 * id(i, j, k); X[q] = fr.c[0] + r * cp * fr.n[0]; X[q + 1] = fr.c[1] + r * cp * fr.n[1]; X[q + 2] = r * sp; } } }
  const elems = [], innerFaces = [], outerFaces = [];
  for (let i = 0; i < ns; i++) for (let j = 0; j < nc; j++) { for (let k = 0; k < nt; k++) elems.push([id(i, j, k), id(i + 1, j, k), id(i + 1, j + 1, k), id(i, j + 1, k), id(i, j, k + 1), id(i + 1, j, k + 1), id(i + 1, j + 1, k + 1), id(i, j + 1, k + 1)]); innerFaces.push({ nodes: [id(i, j, 0), id(i + 1, j, 0), id(i + 1, j + 1, 0), id(i, j + 1, 0)], s: 0.5 * (s[i] + s[i + 1]), phi: ((half ? Math.PI : 2 * Math.PI) * (j + 0.5)) / nc }); outerFaces.push({ nodes: [id(i, j, nt), id(i + 1, j, nt), id(i + 1, j + 1, nt), id(i, j + 1, nt)], s: 0.5 * (s[i] + s[i + 1]), phi: ((half ? Math.PI : 2 * Math.PI) * (j + 0.5)) / nc }); }
  const endFaces = (end) => { const i = end ? ns : 0, out = []; for (let j = 0; j < nc; j++) for (let k = 0; k < nt; k++) out.push({ nodes: [id(i, j, k), id(i, j + 1, k), id(i, j + 1, k + 1), id(i, j, k + 1)] }); return out; };
  const ring = (i, k) => Array.from({ length: ncN }, (_, j) => id(i, j, k));
  return { X, elems, id, ns, nc, ncN, nt, s, iArc: [nL[0], nL[0] + nA], frame, innerFaces, outerFaces, endFaces, ring, half, length: s[ns], legs, arc, R, ro, t, angle: th };
}
/** Consistent nodal forces of a traction field on quadrilateral faces (2 × 2 Gauss). traction(x, normal, face) -> [tx, ty, tz] in Pa; `normal` is the unit normal on the side of `towards` (a point, e.g. the pipe axis side or outside). Adds into f and returns the resultant [Fx, Fy, Fz]. */
export function hexFaceLoad(X, faces, traction, f, towards = null) {
  const sum = [0, 0, 0], g = HEX_G;
  for (const fc of faces) {
    const q = fc.nodes, P = q.map((n) => [X[3 * n], X[3 * n + 1], X[3 * n + 2]]);
    for (const [a, b] of [[-g, -g], [g, -g], [g, g], [-g, g]]) {
      const N = [0.25 * (1 - a) * (1 - b), 0.25 * (1 + a) * (1 - b), 0.25 * (1 + a) * (1 + b), 0.25 * (1 - a) * (1 + b)], da = [-0.25 * (1 - b), 0.25 * (1 - b), 0.25 * (1 + b), -0.25 * (1 + b)], db = [-0.25 * (1 - a), -0.25 * (1 + a), 0.25 * (1 + a), 0.25 * (1 - a)];
      const x = [0, 0, 0], ta = [0, 0, 0], tb = [0, 0, 0]; for (let k = 0; k < 4; k++) for (let c = 0; c < 3; c++) { x[c] += N[k] * P[k][c]; ta[c] += da[k] * P[k][c]; tb[c] += db[k] * P[k][c]; }
      let nx = ta[1] * tb[2] - ta[2] * tb[1], ny = ta[2] * tb[0] - ta[0] * tb[2], nzz = ta[0] * tb[1] - ta[1] * tb[0]; const dA = Math.hypot(nx, ny, nzz); nx /= dA; ny /= dA; nzz /= dA;
      if (towards) { const tw = typeof towards === 'function' ? towards(x, fc) : towards; if (nx * (tw[0] - x[0]) + ny * (tw[1] - x[1]) + nzz * (tw[2] - x[2]) < 0) { nx = -nx; ny = -ny; nzz = -nzz; } }
      const tr = traction(x, [nx, ny, nzz], fc);
      for (let k = 0; k < 4; k++) for (let c = 0; c < 3; c++) { const v = N[k] * tr[c] * dA; f[3 * q[k] + c] += v; sum[c] += v; }
    }
  }
  return sum;
}
/**
 * Pipe segment or elbow of solid elements under internal pressure, an in-plane bending moment or a wall-pressure field (for example mapped from a flow solution).
 * The end s = 0 is held (all displacements, or only the axial one with `ends: 'plane'`); the other end is loaded by the bending stresses of the moment (σ = M·y/I, a pure couple),
 * is free, or is also held axially (`ends: 'plane'`: plane strain along the pipe, for the thick-cylinder check).
 * o: pipeHexMesh options plus { E, nu, p (internal pressure, Pa; number or function (s, φ) -> Pa), M (in-plane moment, N·m, for the full section), shear (N, transverse end force in the plane, full section), ends: 'clamped' | 'plane' | 'both' (both ends held), tol, progress }
 * Returns { mesh, fe, rotation (rad, of the loaded end, least-squares plane through the axial displacements), tip: [ux, uy] (mean displacement of the loaded end), flexibility (rotation·E·I/M per metre of developed length, straight pipe = 1 → use flexFactor),
 *   flexFactor (elbow: (rotation·EI/M − legs)/arc), stressFactor: { hoop, axial, vm } (largest surface stress in the bend over M·ro/I), hoop: { inner, outer } (at mid-length, extrados φ = 0 … intrados, arrays), phi[],
 *   sumLoad: [Fx, Fy, Fz], sumReaction: [Fx, Fy, Fz], I, elements, nodes }.
 */
export function hexPipe(o) {
  const mesh = pipeHexMesh(o), { X, id, ns, ncN, nt, half } = mesh, nN = X.length / 3, ro = mesh.ro, ri = ro - mesh.t, I = (Math.PI / 4) * (ro ** 4 - ri ** 4), fixed = new Uint8Array(3 * nN), f = new Float64Array(3 * nN), ends = o.ends || 'clamped';
  const fr0 = mesh.frame(0), fr1 = mesh.frame(mesh.length);
  for (let j = 0; j < ncN; j++) for (let k = 0; k <= nt; k++) {
    const a = id(0, j, k), b = id(ns, j, k);
    if (ends === 'plane') { fixed[3 * a] = 1; fixed[3 * b] = 1; } else { fixed[3 * a] = fixed[3 * a + 1] = fixed[3 * a + 2] = 1; if (ends === 'both') fixed[3 * b] = fixed[3 * b + 1] = fixed[3 * b + 2] = 1; }
    if (half && (j === 0 || j === ncN - 1)) for (let i = 0; i <= ns; i++) fixed[3 * id(i, j, k) + 2] = 1; // plane of symmetry
  }
  if (ends === 'plane') { const jq = Math.round(mesh.nc / (half ? 2 : 4)); fixed[3 * id(0, jq, 0) + 1] = 1; if (!half) { fixed[3 * id(0, 0, 0) + 2] = 1; fixed[3 * id(0, Math.round(mesh.nc / 2), 0) + 2] = 1; } } // rigid-body motions of the plane-strain slice
  const sum = [0, 0, 0], add = (q) => { sum[0] += q[0]; sum[1] += q[1]; sum[2] += q[2]; };
  if (o.p) { const pf = typeof o.p === 'function' ? o.p : () => +o.p; add(hexFaceLoad(X, mesh.innerFaces, (x, n, fc) => { const p = pf(fc.s, fc.phi, x); return [-p * n[0], -p * n[1], -p * n[2]]; }, f, (x, fc) => mesh.frame(fc.s).c)); } // the normal points to the axis; the pressure pushes away from it
  if (o.M || o.shear) { const A = Math.PI * (ro * ro - ri * ri); add(hexFaceLoad(X, mesh.endFaces(1), (x) => { const y = (x[0] - fr1.c[0]) * fr1.n[0] + (x[1] - fr1.c[1]) * fr1.n[1], sg = (nz(o.M, 0) * y) / I, tv = nz(o.shear, 0) / A; return [sg * fr1.t[0] + tv * fr1.n[0], sg * fr1.t[1] + tv * fr1.n[1], 0]; }, f)); }
  const fe = hexSolve({ X, elems: mesh.elems, E: o.E, nu: o.nu, fixed, f, tol: o.tol, maxIter: o.maxIter, incompatible: o.incompatible, progress: o.progress });
  // rotation and mean displacement of the loaded end (mid-wall ring): axial displacement against the in-plane coordinate
  let sy = 0, syy = 0, su = 0, syu = 0, n = 0, tx = 0, ty = 0; const km = Math.round(nt / 2);
  for (let j = 0; j < ncN; j++) for (let k = 0; k <= nt; k++) { const q = id(ns, j, k), y = (X[3 * q] - fr1.c[0]) * fr1.n[0] + (X[3 * q + 1] - fr1.c[1]) * fr1.n[1], ua = fe.u[3 * q] * fr1.t[0] + fe.u[3 * q + 1] * fr1.t[1]; sy += y; syy += y * y; su += ua; syu += y * ua; n++; tx += fe.u[3 * q]; ty += fe.u[3 * q + 1]; }
  const rotation = (n * syu - sy * su) / (n * syy - sy * sy), Mx = nz(o.M, 0), EI = o.E * I;
  // surface stresses in local axes over the curved part (or the whole pipe when straight)
  const hoopAt = (q, i, j) => { const fr = mesh.frame(mesh.s[i]), ph = ((half ? Math.PI : 2 * Math.PI) * j) / mesh.nc, e = [-Math.sin(ph) * fr.n[0], -Math.sin(ph) * fr.n[1], Math.cos(ph)], s = fe.stress[q]; return { hoop: s[0] * e[0] * e[0] + s[1] * e[1] * e[1] + s[2] * e[2] * e[2] + 2 * (s[3] * e[0] * e[1] + s[4] * e[1] * e[2] + s[5] * e[2] * e[0]), axial: s[0] * fr.t[0] * fr.t[0] + s[1] * fr.t[1] * fr.t[1] + 2 * s[3] * fr.t[0] * fr.t[1] }; };
  let hMax = 0, aMax = 0, vMax = 0; const i0 = mesh.arc > 0 ? mesh.iArc[0] : 1, i1 = mesh.arc > 0 ? mesh.iArc[1] : ns - 1;
  for (let i = Math.max(i0, 0); i <= i1; i++) for (let j = 0; j < ncN; j++) for (const k of [0, nt]) { const q = id(i, j, k), h = hoopAt(q, i, j); hMax = Math.max(hMax, Math.abs(h.hoop)); aMax = Math.max(aMax, Math.abs(h.axial)); vMax = Math.max(vMax, fe.vm[q]); }
  const iM = Math.round((i0 + i1) / 2), phi = [], hi = [], ho = [], ai = [], ao = []; for (let j = 0; j < ncN; j++) { phi.push(((half ? 180 : 360) * j) / mesh.nc); const a = hoopAt(id(iM, j, 0), iM, j), b = hoopAt(id(iM, j, nt), iM, j); hi.push(a.hoop); ho.push(b.hoop); ai.push(a.axial); ao.push(b.axial); }
  const sr = [0, 0, 0]; for (let q = 0; q < nN; q++) for (let c = 0; c < 3; c++) sr[c] += fe.reactions[3 * q + c];
  const nom = Mx ? (Math.abs(Mx) * ro) / I : 1, legsLen = mesh.legs[0] + mesh.legs[1];
  return { mesh, fe, rotation, tip: [tx / n, ty / n], flexibility: Mx ? (rotation * EI) / (Mx * mesh.length) : null, flexFactor: Mx && mesh.arc > 0 ? ((rotation * EI) / Mx - legsLen) / mesh.arc : null, stressFactor: Mx ? { hoop: hMax / nom, axial: aMax / nom, vm: vMax / nom } : null, peak: { hoop: hMax, axial: aMax, vm: vMax },
    hoop: { inner: hi, outer: ho }, axial: { inner: ai, outer: ao }, phi, sumLoad: sum, sumReaction: sr, I, elements: mesh.elems.length, nodes: nN, km };
}
// ---- over-dispersed failure counts: negative-binomial (gamma–Poisson) frequency model with a log-linear trend ------------------
/** Natural logarithm of the gamma function (Lanczos). */
export function lnGamma(z) {
  if (z < 0.5) return Math.log(Math.PI / Math.abs(Math.sin(Math.PI * z))) - lnGamma(1 - z);
  z -= 1; let x = LANCZOS[0]; for (let i = 1; i < 9; i++) x += LANCZOS[i] / (z + i);
  const t = z + 7.5; return 0.5 * Math.log(2 * Math.PI) + (z + 0.5) * Math.log(t) - t + Math.log(x);
}
/** Negative-binomial probability of n events for a mean mu and a shape k (variance mu + mu²/k; k → ∞ is the Poisson distribution). */
export function nbPmf(n, mu, k) {
  if (!(mu > 0)) return n === 0 ? 1 : 0;
  if (!(k < 1e7)) return Math.exp(n * Math.log(mu) - mu - lnGamma(n + 1));
  return Math.exp(lnGamma(n + k) - lnGamma(k) - lnGamma(n + 1) + k * Math.log(k / (k + mu)) + n * Math.log(mu / (k + mu)));
}
// Shape of the gamma–Poisson frequency model: the moment (Pearson) value of the hazardous-liquid incident counts 2010–2024 about their log-linear trend (see nbFit and the reference data sets).
export const FREQ_SHAPE = 183;
/** Regularised lower incomplete gamma function P(a, x) (series for x < a + 1, continued fraction otherwise). */
export function gammaP(a, x) {
  if (!(x > 0)) return 0;
  if (x < a + 1) { let s = 1 / a, d = s; for (let n = 1; n < 500; n++) { d *= x / (a + n); s += d; if (Math.abs(d) < 1e-15 * Math.abs(s)) break; } return s * Math.exp(-x + a * Math.log(x) - lnGamma(a)); }
  let b = x + 1 - a, c = 1e300, d = 1 / b, h = d; for (let i = 1; i < 500; i++) { const an = -i * (i - a); b += 2; d = an * d + b; if (Math.abs(d) < 1e-300) d = 1e-300; c = b + an / c; if (Math.abs(c) < 1e-300) c = 1e-300; d = 1 / d; const dl = d * c; h *= dl; if (Math.abs(dl - 1) < 1e-15) break; }
  return 1 - Math.exp(-x + a * Math.log(x) - lnGamma(a)) * h;
}
/** Quantile of the gamma distribution with shape a and unit scale. */
export const gammaQuantile = (q, a) => brent((x) => gammaP(a, x) - q, 1e-12, a + 40 * Math.sqrt(a) + 50, 1e-10);
/** Quantile of the negative-binomial distribution: the smallest n with P(N ≤ n) ≥ q. */
export function nbQuantile(q, mu, k) { let c = 0, n = 0; const top = Math.ceil(mu + 40 * Math.sqrt(mu + (mu * mu) / Math.min(k, 1e9)) + 50); for (; n <= top; n++) { c += nbPmf(n, mu, k); if (c >= q) return n; } return top; }
/**
 * Fit of yearly counts to a negative-binomial (gamma–Poisson) model with a log-linear trend of the rate: mean_i = exposure_i·exp(a + b·(t_i − t₀)), one shape k for all years.
 * a and b by maximum likelihood; the shape is the smaller of its maximum-likelihood value and the moment value from the Pearson statistic χ²/(n − p) (with few years the likelihood
 * estimate drifts to the Poisson limit). The forecast carries the uncertainty of a and b (inverse Fisher information) as extra dispersion.
 * rows: [{ t (year), n (count), e (exposure) }]; opt: { trend (false = constant rate), t0 }.
 * Returns { a, b, k, kML, t0, rate(t), predict(t, e) -> { mu, k (shape of the predictive distribution), lo, hi (5 % and 95 % quantiles of the count) }, pearson (χ²/(n − p)), logLik, cov: [[..], [..]], n }.
 */
export function nbFit(rows, opt = {}) {
  const t0 = opt.t0 ?? mean(rows.map((r) => r.t)), trend = opt.trend !== false && rows.length >= 3, np = trend ? 2 : 1, r0 = Math.log(sum(rows.map((r) => r.n)) / sum(rows.map((r) => r.e))), KMAX = 1e6;
  const nll = (x) => { const k = Math.exp(clamp(x[2], -3, Math.log(KMAX))); let s = 0; for (const r of rows) s -= Math.log(Math.max(nbPmf(r.n, r.e * Math.exp(x[0] + (trend ? x[1] : 0) * (r.t - t0)), k), 1e-300)); return s; };
  let best = null; for (const lk of [2, 5, 9]) { const f = nelderMead(nll, [r0, 0, lk], { lo: [r0 - 3, -0.5, -3], hi: [r0 + 3, 0.5, Math.log(KMAX)], tol: 1e-10, maxIter: 3000, scale: 0.2 }); if (!best || f.f < best.f) best = f; }
  const [a, b0, lk] = best.x, b = trend ? b0 : 0, kML = Math.exp(clamp(lk, -3, Math.log(KMAX))), muOf = (r) => r.e * Math.exp(a + b * (r.t - t0)), dof = Math.max(rows.length - np, 1);
  const pearson = sum(rows.map((r) => (r.n - muOf(r)) ** 2 / muOf(r))) / dof, mbar = mean(rows.map(muOf)), k = pearson > 1 ? Math.min(kML, mbar / (pearson - 1)) : kML;
  let i00 = 0, i01 = 0, i11 = 0; for (const r of rows) { const mu = muOf(r), w = (mu * k) / (k + mu), x = r.t - t0; i00 += w; i01 += w * x; i11 += w * x * x; }
  const det = trend ? i00 * i11 - i01 * i01 : i00, cov = trend && det > 0 ? [[i11 / det, -i01 / det], [-i01 / det, i00 / det]] : [[1 / i00, 0], [0, 0]];
  const predict = (t, e) => { const x = t - t0, mu = e * Math.exp(a + b * x), vLog = cov[0][0] + 2 * x * cov[0][1] + x * x * cov[1][1], kp = 1 / (1 / k + vLog); return { mu, k: kp, lo: nbQuantile(0.05, mu, kp), hi: nbQuantile(0.95, mu, kp) }; };
  return { a, b, k, kML, t0, rate: (t) => Math.exp(a + b * (t - t0)), predict, pearson, logLik: -best.f, cov, n: rows.length };
}

const opt = (pairs) => pairs.map(([value, label]) => ({ value, label }));
// ---- additional solvers of the set-up tab (run only when selected) -------------------------------------------------------------
const TASKS = opt([['none', 'None: line assessment only'], ['tubeFsi', 'Pressure pulse in the elastic pipe (axisymmetric Navier–Stokes + shell, strongly coupled)'], ['bendCfd', 'Bend: resolved flow → finite elements (one-way)'], ['bendFsi', 'Bend: resolved flow ↔ beam model (two-way)'], ['solid', 'Pipe bend of 3-D solid elements (pressure and bending)']]);
/** Arrival time of a signal at a fraction of its maximum (linear interpolation). */
const arrival = (t, y, frac = 0.5) => { const m = Math.max(...y); if (!(m > 0)) return null; for (let i = 1; i < y.length; i++) if (y[i] >= frac * m) return t[i - 1] + ((t[i] - t[i - 1]) * (frac * m - y[i - 1])) / (y[i] - y[i - 1] || 1); return null; };
/** Wall-pressure field of the plane flow solution on the circumference of the pipe: p(s, φ) = mean + half the extrados–intrados difference × cos φ (φ = 0 at the extrados). The resultant per unit length equals that of the plane solution through the depth π·D/4. */
const mapBendPressure = (L) => (s, phi) => { const i = clamp(L.s.findIndex((x) => x >= s), 0, L.s.length - 1), j = Math.max(i - 1, 0), w = L.s[i] > L.s[j] ? clamp((s - L.s[j]) / (L.s[i] - L.s[j]), 0, 1) : 0, pe = L.pExt[j] + w * (L.pExt[i] - L.pExt[j]), pi = L.pInt[j] + w * (L.pInt[i] - L.pInt[j]); return 0.5 * (pe + pi) + 0.5 * (pe - pi) * Math.cos(phi); };
async function runTask(p, st, ds, sl, span, sa, ctx, prog) {
  const out = { plots: [], tables: [], warnings: [], outputs: {}, text: '' }, tick = async () => { if (ctx.tick) await ctx.tick(); }, nd = ds.node, D = p.ID, t = p.t, cpu0 = Date.now();
  const bend = p.bends.find((b) => b.radius > 0) || { angle: 90, radius: 1.5 * D }, R = Math.max(p.taskRD > 0 ? p.taskRD * D : bend.radius, 1.0 * D), angle = clamp(p.bendAngle > 0 ? p.bendAngle : bend.angle, 10, 170), rhoF = ds.on ? ds.rhoF : nd.rhoNS, rhoS = ds.on ? ds.rhoS : nd.rhoNS, U = ds.on ? ds.v : nd.vm, A = p.Ai;
  if (p.task === 'tubeFsi') {
    prog(0.9, 'Coupled pressure pulse: fluid and shell'); await tick();
    const base = { R: D / 2, h: t, E: p.E, nu: p.nu, rhoS: p.rhoSteel, rhoF: nd.rhoL, mu: nd.muL, K: p.fsiBulk * MPA, axial: 'restrained' }, an = tubeWaves(base), L = p.fsiLen, Tp = L / (4 * an.cKorteweg), nx = p.fsiNx, nr = p.fsiNr, dt = Tp / Math.max(nx / 3, 12), x1 = 0.25 * L, x2 = 0.7 * L, amp = Math.max(p.fsiPulse, 1e-3) * BAR;
    const r = tubeFsi({ ...base, L, nx, nr, dt, tEnd: (0.95 * L) / an.cKorteweg + Tp, inlet: { type: 'pressure', p: (tt) => (tt < Tp ? amp * Math.sin((Math.PI * tt) / Tp) ** 2 : 0) }, stations: [x1, x2], progress: (f) => prog(0.9 + 0.09 * f, 'Coupled pressure pulse: fluid and shell') });
    const tp = [0, 1].map((k) => r.t[r.pWall[k].indexOf(Math.max(...r.pWall[k]))]), cPeak = (x2 - x1) / Math.max(tp[1] - tp[0], 1e-12), t50 = [0, 1].map((k) => arrival(r.t, r.pWall[k])), c50 = t50[0] !== null && t50[1] !== null ? (x2 - x1) / Math.max(t50[1] - t50[0], 1e-12) : cPeak, etaMax = Math.max(...r.eta[0]), hoop = (Math.max(...r.pWall[0]) * (D / 2)) / t;
    out.plots.push({ type: 'line', title: 'Coupled pressure pulse: wall pressure at two stations', xlabel: 'Time (ms)', ylabel: 'Wall pressure (bar)', series: [{ name: `x = ${fmt(x1, 3)} m`, x: r.t.map((q) => q * 1000), y: r.pWall[0].map((q) => q / BAR) }, { name: `x = ${fmt(x2, 3)} m`, x: r.t.map((q) => q * 1000), y: r.pWall[1].map((q) => q / BAR) }], note: `Grid ${nx} × ${nr} cells, ${r.steps} steps of ${fmt(r.dt * 1e6, 3)} µs, ${fmt(r.iterations.mean, 3)} coupling iterations per step.` });
    out.plots.push({ type: 'line', title: 'Coupled pressure pulse: wall displacement along the pipe at the end', xlabel: 'Distance (m)', ylabel: 'Radial displacement (µm)', series: [{ name: 'Shell', x: r.x, y: r.etaEnd.map((q) => q * 1e6) }] });
    out.tables.push({ title: 'Coupled fluid–structure solution: pressure pulse in the elastic pipe', columns: ['Quantity', 'Value', 'Unit'], rows: [['Pipe length modelled / pulse duration', `${fmt(L, 3)} / ${fmt(Tp * 1000, 3)}`, 'm / ms'], ['Grid (axial × radial cells), time steps', `${nx} × ${nr}, ${r.steps}`, '—'], ['Wave speed, peak-to-peak', cell(cPeak), 'm/s'], ['Wave speed, half-rise', cell(c50), 'm/s'], ['Korteweg wave speed (compressible liquid, elastic restrained wall)', cell(an.cKorteweg), 'm/s'], ['Moens–Korteweg wave speed (incompressible liquid)', cell(an.c0), 'm/s'], ['Simulated / Korteweg', cell(cPeak / an.cKorteweg, 4), '–'],
      ['Peak wall displacement at the first station', cell(etaMax * 1e6), 'µm'], ['Peak hoop stress there (p·R/t)', cell(hoop / MPA), 'MPa'], ['Coupling iterations per step (mean / largest)', `${fmt(r.iterations.mean, 3)} / ${r.iterations.max}`, '—'], ['Work of the liquid on the wall / energy gained by the shell', `${fmt(r.energy.interface, 4)} / ${fmt(r.energy.shellGain, 4)}`, 'J'], ['Volume swept by the wall / net inflow', `${fmt(r.volume.wallSwept * 1e6, 4)} / ${fmt(r.volume.boundary * 1e6, 4)}`, 'cm³'], ['Run time', cell(Date.now() - cpu0), 'ms']],
      note: 'Axisymmetric Navier–Stokes equations on a grid that moves with the wall, thin-shell wall, interface Newton iteration in every step. Resolution: the radial grid does not resolve the viscous sub-layer of a turbulent pipe flow; the wave speed, the wall load and the added mass do not depend on it. Refine on the mesh tab (study “Coupled pulse grid”).' });
    out.outputs.fsi = { waveSpeed: sig(cPeak, 5), waveSpeedHalfRise: sig(c50, 5), korteweg: sig(an.cKorteweg, 5), moensKorteweg: sig(an.c0, 5), ratio: sig(cPeak / an.cKorteweg, 5), iterations: sig(r.iterations.mean, 4), converged: r.converged, interfaceWork: sig(r.energy.interface, 5), shellGain: sig(r.energy.shellGain, 5), grid: `${nx} x ${nr}`, steps: r.steps, runMs: Date.now() - cpu0 };
    out.text = ` Coupled pulse: wave speed ${fmt(cPeak, 4)} m/s against ${fmt(an.cKorteweg, 4)} m/s (Korteweg).`;
    if (!r.converged) out.warnings.push({ level: 'warn', msg: 'The coupled pulse solution did not converge in every time step: refine the time step or the grid.' });
  } else if (p.task === 'bendCfd' || p.task === 'bendFsi') {
    const two = p.task === 'bendFsi'; prog(0.9, two ? 'Bend: resolved flow coupled to the beam model' : 'Bend: resolved flow'); await tick();
    const legs = clamp(p.taskLegs, 1, 20), lenTot = 2 * legs * D + (R * angle * Math.PI) / 180, slugLen = ds.on ? Math.min(ds.len, (two ? 1.5 : 4) * lenTot) : 0, tEnd = (lenTot + slugLen) / U + (0.6 * lenTot) / U + (0.3 * lenTot) / U, cfd = { D, R, angle, legs, U, rho: rhoF, rhoSlug: rhoS, slug: slugLen > 0 ? { t0: (0.3 * lenTot) / U, len: slugLen } : null, mu: nd.muL, ns: p.cfdNs, nn: p.cfdNn };
    let r, one = null;
    if (two) { const pre = bendFsi({ ...cfd, U: 0, slug: null, t, E: p.E, nu: p.nu, rhoSteel: p.rhoSteel, mode: 'one', dt: 1, tEnd: 1 }), dt = Math.min(1 / (pre.modes.filled[0] * 12), (0.4 * (lenTot / cfd.ns) * (1 - D / (2 * R))) / U), common = { ...cfd, t, E: p.E, nu: p.nu, rhoSteel: p.rhoSteel, dt, tEnd, zeta: p.zeta };
      one = bendFsi({ ...common, mode: 'one', progress: (f) => prog(0.9 + 0.04 * f, 'Bend: one-way reference') }); await tick(); r = bendFsi({ ...common, mode: 'two', progress: (f) => prog(0.94 + 0.05 * f, 'Bend: two-way coupling') }); }
    else { const c = bendCfd({ ...cfd, tEnd, progress: (f) => prog(0.9 + 0.07 * f, 'Bend: resolved flow') }); r = { t: c.t, fFluid: c.f, cv: c.cvx.map((x, i) => Math.hypot(x, c.cvy[i])), last: c.last, fields: c.fields, grid: c.grid, steps: c.steps, theory: c.theory }; }
    const iPk = r.fFluid.reduce((k, v, i) => (v > r.fFluid[k] ? i : k), 0), fPk = r.fFluid[iPk], kt = Math.sqrt(2 * (1 - Math.cos((angle * Math.PI) / 180))), fFilm = rhoF * A * U * U * kt, fSlug = rhoS * A * U * U * kt, g = r.grid, Lw = r.last;
    out.plots.push({ type: 'line', title: 'Force of the flow on the bend', xlabel: 'Time (s)', ylabel: 'Resultant force (kN)', series: [{ name: two ? 'Wall tractions, two-way' : 'Integrated wall tractions', x: thin(r.t), y: thin(r.fFluid).map((q) => q / 1000) }, ...(two ? [{ name: 'Wall tractions, one-way', x: thin(one.t), y: thin(one.fFluid).map((q) => q / 1000), dash: true }, { name: 'Anchor reactions, two-way', x: thin(r.t), y: thin(r.reaction).map((q) => q / 1000) }] : [{ name: 'Control-volume momentum theorem', x: thin(r.t), y: thin(r.cv).map((q) => q / 1000), dash: true }])], hlines: [{ y: fSlug / 1000, label: 'ρ·A·v²·√(2(1 − cos θ)), slug body' }, { y: fFilm / 1000, label: 'the same, between slugs' }] });
    out.plots.push({ type: 'line', title: 'Wall pressure and shear along the bend at the end of the run', xlabel: 'Distance along the centreline (m)', ylabel: 'Pressure (kPa) / shear (Pa)', series: [{ name: 'Pressure, extrados (kPa)', x: Lw.s, y: Lw.pExt.map((q) => q / 1000) }, { name: 'Pressure, intrados (kPa)', x: Lw.s, y: Lw.pInt.map((q) => q / 1000) }, { name: 'Shear, extrados (Pa)', x: Lw.s, y: Lw.tauExt, dash: true }, { name: 'Shear, intrados (Pa)', x: Lw.s, y: Lw.tauInt, dash: true }], vlines: [{ x: g.legs, label: 'bend starts' }, { x: g.legs + g.arc, label: 'bend ends' }] });
    out.plots.push({ type: 'field', title: 'Pressure in the bend (developed along the centreline)', xlabel: 'Distance along the centreline (m)', ylabel: 'Across the bore (m), extrados up', zlabel: 'Gauge pressure', zunit: 'Pa', x: g.s, y: g.n, z: r.fields.p, cmap: 'coolwarm' });
    // one-way mapping of the peak wall load onto the solid elbow and the beam model: equilibrium of the tractions with the reactions
    prog(0.985, 'Bend: finite elements under the mapped wall load'); await tick();
    const pm = mapBendPressure(Lw), hx = hexPipe({ ro: p.D / 2, t, R, angle, legs: [legs * D, legs * D], nc: p.hexNc, nt: p.hexNt, nArc: p.hexNa, nLeg: Math.max(2, Math.round(p.hexNa / 2)), E: p.E, nu: p.nu, p: (s, phi) => pm(s, phi), ends: 'both', tol: 1e-7 }), fHex = 2 * Math.hypot(hx.sumLoad[0], hx.sumLoad[1]), rHex = 2 * Math.hypot(hx.sumReaction[0], hx.sumReaction[1]), fEnd = Math.hypot(Lw.force[0] - (two ? 0 : 0), Lw.force[1]);
    const pressOnly = b31ThrustCheck(Lw, g);
    out.tables.push({ title: two ? 'Bend: resolved flow and beam model, two-way coupled' : 'Bend: resolved flow mapped onto finite elements (one-way)', columns: ['Quantity', 'Value', 'Unit'], rows: [['Bore, bend radius, angle, straight legs', `${fmt(D * 1000, 4)} mm, ${fmt(R / D, 3)} D, ${fmt(angle, 3)}°, ${fmt(legs, 3)} D`, '—'], ['Flow grid (along × across), time steps', `${g.ns} × ${g.nn}, ${r.steps}`, '—'], ['Velocity, density between slugs / in the slug', `${fmt(U, 3)} m/s, ${fmt(rhoF, 4)} / ${fmt(rhoS, 4)} kg/m³`, '—'],
      ['Peak force from the wall tractions', cell(fPk / 1000), 'kN'], ['Momentum theorem ρ·A·v²·√(2(1 − cos θ)): between slugs / slug body', `${fmt(fFilm / 1000, 4)} / ${fmt(fSlug / 1000, 4)}`, 'kN'], ['Peak / slug-body momentum force (dynamic factor of the passing front)', cell(fPk / Math.max(fSlug, 1e-9), 4), '–'], ['Force at the end of the run: tractions / control volume', `${fmt(Math.hypot(Lw.force[0], Lw.force[1]) / 1000, 4)} / ${fmt(Math.hypot(Lw.cv[0], Lw.cv[1]) / 1000, 4)}`, 'kN'],
      ['Largest extrados − intrados pressure at the end', cell(Math.max(...Lw.pExt.map((q, i) => q - Lw.pInt[i])) / 1000), 'kPa'], ['Centrifugal estimate ρ·v²·D/R', cell((rhoF * U * U * D) / R / 1000), 'kPa'], ['Largest wall shear at the end', cell(Math.max(...Lw.tauExt, ...Lw.tauInt)), 'Pa'],
      ['Solid elbow under the mapped wall pressure: applied resultant / anchor reactions', `${fmt(fHex / 1000, 4)} / ${fmt(rHex / 1000, 4)}`, 'kN'], ['Solid elbow: plane flow resultant of the same pressure field', cell(pressOnly / 1000), 'kN'], ['Solid elbow: largest von Mises stress in the bend', cell(hx.peak.vm / MPA), 'MPa'], ['Solid elements / conjugate-gradient iterations', `${hx.elements} / ${hx.fe.iterations}`, '—'],
      ...(two ? [['First in-plane frequency of the spool: empty / filled (beam model)', `${fmt(r.modes.empty[0], 4)} / ${fmt(r.modes.filled[0], 4)}`, 'Hz'], ['Peak displacement of the bend: two-way / one-way', `${fmt(Math.max(...r.disp) * 1e6, 4)} / ${fmt(Math.max(...one.disp) * 1e6, 4)}`, 'µm'], ['Peak bending moment: two-way / one-way', `${fmt(Math.max(...r.moment) / 1000, 4)} / ${fmt(Math.max(...one.moment) / 1000, 4)}`, 'kN·m'], ['Peak anchor reaction: two-way / one-way', `${fmt(Math.max(...r.reaction) / 1000, 4)} / ${fmt(Math.max(...one.reaction) / 1000, 4)}`, 'kN'], ['Coupling iterations per step (mean / largest)', `${fmt(r.iterations.mean, 3)} / ${r.iterations.max}`, '—'], ['Work of the flow on the beam / energy gained by the beam', `${fmt(r.energy.work, 4)} / ${fmt(r.energy.gain, 4)}`, 'J'], ['Mass of the contents / steel', cell(r.massRatio, 3), '–']] : []), ['Run time', cell(Date.now() - cpu0), 'ms']],
      note: 'Plane (two-dimensional) incompressible flow on a body-fitted grid, mixing-length eddy viscosity with a wall law, the slug as a density wave; the plane solution stands for the pipe through the depth π·D/4 and its wall pressure is spread over the circumference with a cosine. Secondary flow, the gas–liquid interface inside the slug and wall turbulence are not resolved at this resolution: the force and the pressure difference across the bend are converged to a few per cent, the wall shear is an estimate. Production resolution: the OpenFOAM + CalculiX hand-off of the external-solvers page.' });
    out.outputs.cfdFea = { mode: two ? 'two-way' : 'one-way', peakForce: sig(fPk / 1000, 5), momentumForceSlug: sig(fSlug / 1000, 5), momentumForceFilm: sig(fFilm / 1000, 5), dynamicFactor: sig(fPk / Math.max(fSlug, 1e-9), 5), tractionForce: sig(Math.hypot(Lw.force[0], Lw.force[1]) / 1000, 5), controlVolumeForce: sig(Math.hypot(Lw.cv[0], Lw.cv[1]) / 1000, 5), solidApplied: sig(fHex / 1000, 5), solidReaction: sig(rHex / 1000, 5), solidPeakStress: sig(hx.peak.vm / MPA, 5), maxShear: sig(Math.max(...Lw.tauExt, ...Lw.tauInt), 5), grid: `${g.ns} x ${g.nn}`, steps: r.steps,
      peakDisplacement: two ? sig(Math.max(...r.disp), 5) : null, peakDisplacementOneWay: two ? sig(Math.max(...one.disp), 5) : null, peakMoment: two ? sig(Math.max(...r.moment) / 1000, 5) : null, iterations: two ? sig(r.iterations.mean, 4) : null, converged: two ? r.converged : true, work: two ? sig(r.energy.work, 5) : null, energyGain: two ? sig(r.energy.gain, 5) : null, runMs: Date.now() - cpu0 };
    out.text = ` Resolved bend flow (${two ? 'two-way' : 'one-way'}): peak wall force ${fmt(fPk / 1000, 3)} kN against ${fmt(fSlug / 1000, 3)} kN from the momentum theorem.`;
    if (two && !r.converged) out.warnings.push({ level: 'warn', msg: 'The two-way bend solution did not converge in every time step: reduce the time step.' });
  } else if (p.task === 'solid') {
    prog(0.9, 'Solid elements: pressure'); await tick();
    const ro = p.D / 2, legsM = [p.taskLegs * D, p.taskLegs * D], mesh = { ro, t, R, angle, legs: legsM, nc: p.hexNc, nt: p.hexNt, nArc: p.hexNa, nLeg: Math.max(2, Math.round(p.hexNa / 2)), E: p.E, nu: p.nu }, row = sa.rows[sa.iH], dp = Math.max(sa.pLocal(row.nd) - row.nd.pe, 1), Mref = 0.1 * p.S * ((Math.PI / 4) * (ro ** 4 - (ro - t) ** 4)) / ro;
    const hp = hexPipe({ ...mesh, p: dp, ends: 'both', progress: (f) => prog(0.9 + 0.04 * f, 'Solid elements: pressure') }); await tick();
    const hm = hexPipe({ ...mesh, M: Mref, progress: (f) => prog(0.94 + 0.05 * f, 'Solid elements: in-plane bending') }), ef = elbowFactors({ t, R, r: ro - t / 2 }), la = lame(dp, 0, ro - t, ro, ro - t), lo = lame(dp, 0, ro - t, ro, ro), kP = 1 + 6 * (dp / p.E) * ((ro - t / 2) / t) ** (7 / 3) * (R / (ro - t / 2)) ** (1 / 3);
    out.plots.push({ type: 'line', title: 'Solid elbow: hoop stress around the mid-bend section', xlabel: 'Angle from the extrados (°)', ylabel: 'Hoop stress (MPa)', series: [{ name: 'Pressure, inner surface', x: hp.phi, y: hp.hoop.inner.map((q) => q / MPA) }, { name: 'Pressure, outer surface', x: hp.phi, y: hp.hoop.outer.map((q) => q / MPA) }, { name: 'Bending, inner surface', x: hm.phi, y: hm.hoop.inner.map((q) => q / MPA), dash: true }, { name: 'Bending, outer surface', x: hm.phi, y: hm.hoop.outer.map((q) => q / MPA), dash: true }], hlines: [{ y: la.hoop / MPA, label: 'Lamé, straight pipe bore' }], note: `In-plane moment ${fmt(Mref / 1000, 4)} kN·m (10 % of first yield of the straight pipe); ovalisation of the bend shows as hoop bending of the wall.` });
    out.tables.push({ title: 'Pipe bend of three-dimensional solid elements', columns: ['Quantity', 'Value', 'Unit'], rows: [['Outer diameter, wall, bend radius, angle, legs', `${fmt(p.D * 1000, 4)} mm, ${fmt(t * 1000, 3)} mm, ${fmt(R / D, 3)} D, ${fmt(angle, 3)}°, ${fmt(p.taskLegs, 3)} D`, '—'], ['Elements (half model: around × through the wall × along)', `${hp.elements} (${p.hexNc} × ${p.hexNt} × ${hp.mesh.ns})`, '—'], ['Degrees of freedom / conjugate-gradient iterations (pressure, bending)', `${hp.fe.ndof} / ${hp.fe.iterations}, ${hm.fe.iterations}`, '—'],
      ['Pipe characteristic h = t·R/r²', cell(ef.h), '–'], ['Flexibility factor: solid elements', cell(hm.flexFactor), '–'], ['Flexibility factor: ASME B31 (1.65/h)', cell(ef.k), '–'], ['Flexibility factor: von Kármán, first approximation', cell(ef.kKarman), '–'], ['Pressure stiffening of the B31 factor at the design pressure', cell(1 / kP, 4), '–'],
      ['Peak hoop stress in bending / (M·r_o/I): solid elements', cell(hm.stressFactor.hoop), '–'], ['Peak longitudinal stress in bending / (M·r_o/I): solid elements', cell(hm.stressFactor.axial), '–'], ['Stress-intensification factor, in-plane: ASME B31 (0.9/h^⅔)', cell(ef.iIn), '–'], ['Theoretical peak stress factor 1.8/h^⅔ (twice the B31 factor)', cell(2 * ef.iIn), '–'],
      ['Pressure: peak hoop stress in the bend, solid elements', cell(hp.peak.hoop / MPA), 'MPa'], ['Pressure: Lamé hoop stress of the straight pipe at the bore', cell(la.hoop / MPA), 'MPa'], ['Pressure: torus factor at the intrados (2R − r)/(2(R − r))', cell((2 * R - (ro - t / 2)) / (2 * (R - (ro - t / 2))), 4), '–'], ['Pressure: applied resultant against the anchor reactions (half model)', `${fmt(Math.hypot(...hp.sumLoad) / 1000, 4)} / ${fmt(Math.hypot(...hp.sumReaction) / 1000, 4)}`, 'kN'], ['Run time', cell(Date.now() - cpu0), 'ms']],
      note: 'Eight-node hexahedra with incompatible bending modes, half model about the plane of the bend, preconditioned conjugate gradients. The straight legs are held at one end; their length and the end restraint lower the flexibility factor below the value of a bend between long flexible pipes, which is what the code factors assume. Linear elasticity: no plasticity, no pressure stiffening of the bending response. Refine on the mesh tab (study “Solid elements of the bend”).' });
    void lo;
    out.outputs.solid = { flexibilityFactor: sig(hm.flexFactor, 5), flexibilityB31: sig(ef.k, 5), flexibilityKarman: sig(ef.kKarman, 5), stressFactorHoop: sig(hm.stressFactor.hoop, 5), stressFactorAxial: sig(hm.stressFactor.axial, 5), sifB31: sig(ef.iIn, 5), pressureHoopPeak: sig(hp.peak.hoop / MPA, 5), lameHoop: sig(la.hoop / MPA, 5), elements: hp.elements, iterations: hm.fe.iterations, h: sig(ef.h, 5), runMs: Date.now() - cpu0 };
    out.text = ` Solid elbow: flexibility factor ${fmt(hm.flexFactor, 3)} (ASME B31 ${fmt(ef.k, 3)}), peak hoop stress factor ${fmt(hm.stressFactor.hoop, 3)}.`;
  }
  return out;
}
/** Resultant (N) of the mean wall pressure of a plane bend solution on the pipe: Σ b·p̄·(h_e − h_i)·ds along the normal, the pressure-thrust part of the wall force. */
function b31ThrustCheck(L, g) { let fx = 0, fy = 0; for (let i = 0; i < L.s.length; i++) { const tg = g.tangent(L.s[i]), s = L.s[i], inB = s > g.legs && s < g.legs + g.arc, he = inB ? 1 + g.D / (2 * g.R) : 1, hi = inB ? 1 - g.D / (2 * g.R) : 1, fn = g.depth * (L.pExt[i] * he - L.pInt[i] * hi) * g.ds[i]; fx += fn * -tg[1]; fy += fn * tg[0]; } return Math.hypot(fx, fy); }

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
    { key: 'fabrication', label: 'Pipe manufacturing process', type: 'select', value: 'seamless', options: opt([['seamless', 'Seamless (α_fab = 1.00)'], ['uo', 'UO, TRB or ERW (α_fab = 0.93)'], ['uoe', 'UOE (α_fab = 0.85)']]), help: 'Fabrication factor of DNV-ST-F101 on the collapse and propagation resistance.' },
  ] },
  { group: 'Material condition', tab: 'inputs', help: 'State of the steel and of the girth welds today. The defaults describe new pipe.', fields: [
    { key: 'strengthLoss', label: 'Loss of strength from ageing or damage', unit: '%', value: 0, min: 0, max: 40, help: 'Reduces the yield and tensile strength used in every check (thermal ageing, hydrogen damage, fire exposure).' },
    { key: 'residualStress', label: 'Welding residual stress at the girth welds', unit: 'MPa', value: 0, min: 0, max: 700, help: 'Added to the peak stress of the fracture check. As-welded joints can carry residual stress up to yield; 0 = stress-relieved or not considered.' },
    { key: 'hiLoMm', label: 'Girth-weld misalignment (hi-lo)', unit: 'mm', value: 1, min: 0, max: 5, help: 'Gives a stress concentration factor 1 + (3δ/t)·exp(−(D/t)^−0.5) that multiplies the fatigue stress.' },
    { key: 'hardenN', label: 'Strain-hardening exponent (Ramberg–Osgood)', unit: '–', value: 15, min: 3, max: 60, help: 'ε = σ/E + 0.002·(σ/σ_y)^n. Used for the strain at yield-level stress and when fitting tensile tests.' },
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
    { key: 'pMinShut', label: 'Minimum internal pressure kept during a shutdown', unit: 'bara', value: 10, min: 1, max: 500, help: 'Operating rule at the deepest point: the line is not vented below this pressure. Credited in the collapse check of the corroded wall; the installation case (empty pipe, new wall) is always checked as well.' },
  ] },
  { group: 'Buckle arrestors and lateral-buckling management', tab: 'inputs', help: 'Design measures of a deep-water line: arrestors confine a propagating buckle, planned initiators share the thermal expansion between controlled lateral buckles.', fields: [
    { key: 'arrestors', label: 'Buckle arrestors', type: 'select', value: 'auto', options: opt([['auto', 'Installed wherever a buckle could propagate'], ['none', 'None']]) },
    { key: 'arrSpacing', label: 'Arrestor spacing', unit: 'm', value: 300, min: 12, max: 5000 },
    { key: 'arrThick', label: 'Arrestor wall thickness / pipe wall thickness', unit: '–', value: 2, min: 1.1, max: 5 },
    { key: 'arrLength', label: 'Arrestor length', unit: 'm', value: 1.5, min: 0.1, max: 12 },
    { key: 'replaceCost', label: 'Cost of replacing damaged pipe', unit: 'k$/m', value: 20, min: 0, max: 500, help: 'Sets the residual consequence of a confined buckle (one arrestor spacing) against an unconfined one.' },
    { key: 'buckleMgmt', label: 'Lateral buckling', type: 'select', value: 'planned', options: opt([['planned', 'Planned buckle initiators (sleepers or snake lay)'], ['none', 'Not managed (a single buckle takes all the feed-in)']]) },
    { key: 'initSpacing', label: 'Spacing of the buckle initiators', unit: 'm', value: 2000, min: 100, max: 20000 },
    { key: 'strainAllow', label: 'Allowable bending strain in a planned buckle', unit: '%', value: 0, min: 0, max: 5, help: '0 = strain capacity of the DNV-ST-F101 displacement-controlled criterion divided by the strain resistance factor of the safety class.' },
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
    { key: 'slugCf', label: 'Slug-force coefficient', unit: '–', value: 1, min: 0.2, max: 3, help: 'Multiplies ρ·A·v² in the bend force; 1 is the momentum-balance value. Fit it to measured bend forces on the calibration tab.' },
    { key: 'dlf', label: 'Dynamic load factor', unit: '–', value: 2, min: 1, max: 4, help: '2 is the step-load upper bound for an undamped support.' },
    { key: 'bends', label: 'Bends', type: 'table', columns: [{ key: 'x', label: 'Distance', unit: 'm' }, { key: 'angle', label: 'Angle', unit: 'deg' }, { key: 'radius', label: 'Radius', unit: 'm' }], value: [{ x: 0, angle: 90, radius: 1.27 }, { x: 19560, angle: 90, radius: 1.27 }] },
    { key: 'pipingSupport', label: 'Support arrangement (screening)', type: 'select', value: 'mediumStiff', options: opt([['stiff', 'Stiff'], ['mediumStiff', 'Medium stiff'], ['medium', 'Medium'], ['flexible', 'Flexible']]) },
  ] },
  { group: 'Equipment and connection interfaces', tab: 'inputs', help: 'Equipment items along the line with the loads their nozzles or connectors may take. The slug and momentum force at each location is compared with the allowable force, and force × lever arm with the allowable moment.', fields: [
    { key: 'equipment', label: 'Equipment', type: 'table', columns: [{ key: 'name', label: 'Item', type: 'text' }, { key: 'x', label: 'Distance', unit: 'm' }, { key: 'size', label: 'Nozzle or hub bore', unit: 'mm' }, { key: 'arm', label: 'Lever arm to the nearest bend', unit: 'm' }, { key: 'fAllow', label: 'Allowable force', unit: 'kN' }, { key: 'mAllow', label: 'Allowable moment', unit: 'kN·m' }],
      value: [{ name: 'Manifold connector (subsea)', x: 0, size: 254, arm: 1.5, fAllow: 150, mAllow: 200 }, { name: 'Riser hang-off / topsides flange', x: 19560, size: 254, arm: 2, fAllow: 60, mAllow: 80 }, { name: 'Separator inlet nozzle', x: 19570, size: 254, arm: 1, fAllow: 20, mAllow: 15 }] },
  ] },
  { group: 'Free span or jumper', tab: 'inputs', fields: [
    { key: 'spanLength', label: 'Span length', unit: 'm', value: 15, min: 1, max: 400 },
    { key: 'spanX', label: 'Span location along the line', unit: 'm', value: 9000, min: 0, max: 1e6 },
    { key: 'spanEnds', label: 'End conditions', type: 'select', value: 'pinned-pinned', options: opt([['pinned-pinned', 'Pinned – pinned'], ['fixed-fixed', 'Fixed – fixed'], ['fixed-pinned', 'Fixed – pinned'], ['fixed-free', 'Fixed – free (cantilever)'], ['fixed-guided', 'Fixed – guided (sliding end, no rotation)'], ['connection', 'Connections of finite rotational stiffness (flanges, connectors)'], ['springs', 'Elastic supports (springs)']]) },
    { key: 'spanKt', label: 'Support stiffness, translation', unit: 'kN/m', value: 5000, min: 0, max: 1e7, help: 'End springs, and the default stiffness of a guide.' },
    { key: 'spanKr', label: 'Support or connection stiffness, rotation', unit: 'kN·m/rad', value: 20000, min: 0, max: 1e8, showIf: (v) => v.spanEnds === 'springs' || v.spanEnds === 'connection' },
    { key: 'lineSupports', label: 'Supports and restraints along the line', type: 'table', columns: [{ key: 'x', label: 'Distance', unit: 'm' }, { key: 'type', label: 'Type', type: 'text' }, { key: 'k', label: 'Stiffness', unit: 'kN/m' }], value: [{ x: 0, type: 'PLET on a sliding mudmat', k: 2000 }, { x: BASE.riserBaseX, type: 'riser-base anchor', k: 50000 }, { x: BASE.profile[BASE.profile.length - 1].x, type: 'hang-off flex joint', k: 100000 }], help: 'Linked from the network suite. A support inside the span entered above splits it: the longest part is assessed.' },
    { key: 'midSupports', label: 'Intermediate supports, guides, anchors and lumped masses', type: 'table', columns: [{ key: 'x', label: 'Position from the left end', unit: 'm' }, { key: 'type', label: 'Type: support, guide, anchor or mass', type: 'text' }, { key: 'k', label: 'Stiffness (guide or spring support)', unit: 'kN/m' }, { key: 'mass', label: 'Lumped mass', unit: 'kg' }], value: [] },
    { key: 'fsi', label: 'Fluid–structure coupling', type: 'select', value: 'two', options: opt([['two', 'Two-way (pipe conveying slug flow)'], ['one', 'One-way (the fluid only loads the pipe)']]), help: 'Two-way: the conveyed mass, the Coriolis and the centrifugal forces follow the slugs and act back on the pipe motion.' },
    { key: 'vibAllow', label: 'Allowable vibration velocity', unit: 'mm/s rms', value: 0, min: 0, max: 500, help: '0 = frequency-dependent screening lines for process pipework (concern and problem levels).' },
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
    { key: 'weldBasis', label: 'Girth-weld class', type: 'select', value: 'auto', options: opt([['auto', 'From the hi-lo tolerance of the weld procedure (single-sided weld root)'], ['manual', 'As selected below']]), help: 'Automatic: class E for a mean hi-lo up to 1 mm, F up to 2 mm, F1 up to 3 mm, with the misalignment stress concentration added.' },
    { key: 'snClass', label: 'S–N class', type: 'select', value: 'F1', options: SN_OPTIONS, showIf: (v) => v.weldBasis === 'manual' },
    { key: 'spanMgmt', label: 'Span rectification', type: 'select', value: 'rectify', options: opt([['rectify', 'Supports are added where a span would not meet the fatigue criterion'], ['none', 'None: the span is assessed as entered']]) },
    { key: 'fatMargin', label: 'Margin on the fatigue criterion', unit: '–', value: 2, min: 1, max: 20, help: 'The damage with the design fatigue factor over the remaining life must stay below 1 / margin; a span that fails is shortened to the allowable length (support spacing).' },
    { key: 'snEnv', label: 'Environment', type: 'select', value: 'cp', options: opt([['air', 'Air'], ['cp', 'Seawater with cathodic protection'], ['free', 'Seawater, free corrosion']]) },
    { key: 'scf', label: 'Stress concentration factor', unit: '–', value: 1, min: 1, max: 6, help: 'Geometric concentration other than the weld misalignment, which is added from the hi-lo value.' },
    { key: 'dff', label: 'Design fatigue factor', unit: '–', value: 6, min: 1, max: 10 },
    { key: 'priorDamage', label: 'Fatigue damage already accumulated', unit: '–', value: 0, min: 0, max: 1 },
    { key: 'eventsShutdown', label: 'Shutdowns', unit: '1/y', value: 6, min: 0, max: 365 },
    { key: 'eventsBlowdown', label: 'Blowdowns', unit: '1/y', value: 1, min: 0, max: 100 },
    { key: 'flawMm', label: 'Initial flaw depth', unit: 'mm', value: 1, min: 0.05, max: 20 },
    { key: 'parisC', label: 'Paris coefficient C', unit: 'm/cycle, MPa√m', value: 1.65e-11, min: 1e-14, max: 1e-8, help: '1.65e-11 with m = 3 is the simplified law for steels in air; 7.27e-11 for free corrosion in a marine environment (BS 7910 simplified laws).' },
    { key: 'parisM', label: 'Paris exponent m', unit: '–', value: 3, min: 2, max: 5 },
    { key: 'dKth', label: 'Threshold ΔK', unit: 'MPa√m', value: 2, min: 0, max: 10 },
    { key: 'kMat', label: 'Fracture toughness', unit: 'MPa√m', value: 150, min: 20, max: 400 },
    { key: 'snShift', label: 'S–N curve shift Δlog a', unit: '–', value: 0, min: -1.5, max: 1.5, help: 'Added to log a of the selected class (0 = the design curve, mean minus two standard deviations; +0.4 = the mean curve). Fit it to fatigue tests on the calibration tab.' },
  ] },
  { group: 'Corrosion', tab: 'inputs', fields: [
    { key: 'co2', label: 'CO2 in the gas phase', unit: 'mol %', value: 2, min: 0, max: 80 },
    { key: 'h2s', label: 'H2S in the gas phase', unit: 'mol %', value: 0, min: 0, max: 40 },
    { key: 'corrModel', label: 'CO2 corrosion model', type: 'select', value: 'auto', options: opt([['auto', 'Automatic: NORSOK M-506 in its range, mixed-potential scaling outside it and in brine'], ['norsok', 'NORSOK M-506 (with wall shear stress and pH)'], ['mechanistic', 'Mixed-potential (electrochemical) model, bare steel'], ['1995', 'de Waard–Lotz–Dugstad 1995 (with flow velocity)'], ['1991', 'de Waard–Milliams 1991 (nomogram equation)']]), help: 'Automatic follows temperature and water chemistry: the standard where it is valid (5–150 °C, pH 3.5–6.5, CO2 fugacity 0.1–10 bar), scaled by the electrochemical model below 5 °C, outside the pH or fugacity limits and above 1 wt % NaCl.' },
    { key: 'phAct', label: 'In-situ pH', unit: '–', value: 0, min: 0, max: 8, help: '0 = from the CO2 equilibrium of the water with its bicarbonate.' },
    { key: 'bicarb', label: 'Bicarbonate in the produced water', unit: 'mg/L', value: 300, min: 0, max: 20000, help: 'Alkalinity of the formation water as HCO3⁻; it buffers the pH (0 = condensed water). Take it from the water analysis.' },
    { key: 'brineWt', label: 'Salinity of the produced water', unit: 'wt % NaCl', value: 3.5, min: 0, max: 26, help: 'Chloride retards the electrode reactions: above 1 wt % the rate is reduced with the mixed-potential model.' },
    { key: 'corrControl', label: 'Corrosion control basis', type: 'select', value: 'availability', options: opt([['availability', 'Inhibitor availability model (inhibited rate while injecting, bare rate otherwise)'], ['efficiency', 'Inhibitor efficiency × availability (fraction of the bare rate removed)'], ['none', 'No inhibition']]), help: 'Availability model: rate = A × inhibited rate + (1 − A) × uninhibited rate. Design practice limits A to 95 % for an ordinary injection system and allows up to 99 % for a system with redundant pumps, automatic alarms and a shutdown on loss of injection.' },
    { key: 'inhibAvail', label: 'Inhibitor availability', unit: '%', value: 99, min: 0, max: 99.9, showIf: (v) => v.corrControl === 'availability', help: 'Fraction of the life with inhibitor at or above the minimum dose. 99 % assumes the enhanced injection system; field experience with ordinary systems is 85–95 %.' },
    { key: 'inhibRate', label: 'Inhibited corrosion rate', unit: 'mm/y', value: 0, min: 0, max: 1, showIf: (v) => v.corrControl === 'availability', help: '0 = by temperature: 0.05 mm/y up to 70 °C rising to 0.2 mm/y at 120–150 °C. A practical minimum of 0.1 mm/y is used by some operators.' },
    { key: 'inhibEff', label: 'Inhibitor efficiency × availability', unit: '%', value: 98.5, min: 0, max: 99.9, help: 'Used by the efficiency basis and by the calibration tab.' },
    { key: 'phStab', label: 'pH stabilisation', type: 'select', value: 'auto', options: opt([['auto', 'Automatic: on when at least 20 wt % glycol circulates and the water carries no formation-water bicarbonate'], ['on', 'On'], ['off', 'Off']]), help: 'Alkali added to a regenerated glycol loop raises the pH so that a protective iron-carbonate film forms. Not used with formation water (scale).' },
    { key: 'phStabTarget', label: 'Stabilised pH', unit: '–', value: 6.5, min: 5, max: 8, showIf: (v) => v.phStab !== 'off' },
    { key: 'cladMode', label: 'Corrosion-resistant alloy (clad or lined pipe)', type: 'select', value: 'auto', options: opt([['auto', 'Automatic: where the controlled carbon-steel rate exceeds the target'], ['length', 'Over the length given below'], ['none', 'None']]) },
    { key: 'cladFrom', label: 'Clad section starts at', unit: 'm', value: 0, min: 0, max: 1e6, showIf: (v) => v.cladMode === 'length' },
    { key: 'cladLength', label: 'Clad length', unit: 'm', value: 0, min: 0, max: 1e6, showIf: (v) => v.cladMode === 'length' },
    { key: 'corrTarget', label: 'Target corrosion rate of carbon steel', unit: 'mm/y', value: 0.1, min: 0.01, max: 2, help: 'Sections whose controlled rate exceeds it are clad in the automatic mode.' },
    { key: 'corrMult', label: 'Corrosion-model multiplier', unit: '–', value: 1, min: 0.05, max: 5 },
    { key: 'glycolWt', label: 'Glycol in the water phase', unit: 'wt %', value: 0, min: 0, max: 95 },
    { key: 'wetModel', label: 'Water wetting', type: 'select', value: 'water', options: opt([['water', 'Always water-wet (conservative)'], ['dewaard', 'Oil-wet when water cut < 30 % and liquid velocity > 1 m/s']]) },
    { key: 'tauCrit', label: 'Wall shear that strips the inhibitor film', unit: 'Pa', value: 150, min: 5, max: 2000 },
    { key: 'minWt', label: 'Measured minimum wall thickness', unit: 'mm', value: 0, min: 0, max: 80, help: '0 = no measurement: the wall is predicted from the age and the rates.' },
    { key: 'wtMap', label: 'Wall-thickness or metal-loss map', type: 'file', value: null },
    { key: 'defects', label: 'Measured metal-loss defects', type: 'table', columns: [{ key: 'x', label: 'Location', unit: 'm' }, { key: 'depth', label: 'Depth', unit: 'mm' }, { key: 'length', label: 'Axial length', unit: 'mm' }], value: [] },
    { key: 'defectGrowth', label: 'Defect-growth multiplier', unit: '–', value: 1, min: 0.05, max: 10, help: 'Multiplies the growth rate of the measured defects (pits usually grow faster than the general wall loss). Fit it to repeated inspections.' },
    { key: 'dentMm', label: 'Dent depth', unit: 'mm', value: 0, min: 0, max: 100, help: 'A plain dent for the ring model of the cross-section; 0 = none.' },
    { key: 'feDepthMm', label: 'Thinned band for the wall model: depth', unit: 'mm', value: 0, min: 0, max: 60, help: '0 = the deepest measured defect, or the wall loss predicted at the end of the design life.' },
    { key: 'feLenMm', label: 'Thinned band for the wall model: axial length', unit: 'mm', value: 300, min: 5, max: 5000 },
  ] },
  { group: 'Sand and erosion', tab: 'inputs', fields: [
    { key: 'sandKgD', label: 'Sand production', unit: 'kg/d', value: 10, min: 0, max: 50000 },
    { key: 'sandUm', label: 'Particle size', unit: 'µm', value: 250, min: 10, max: 3000 },
    { key: 'sandDensity', label: 'Particle density', unit: 'kg/m³', value: 2650, min: 1000, max: 8000 },
    { key: 'erosionModel', label: 'Erosion model that governs', type: 'select', value: 'governing', options: opt([['governing', 'Larger of DNV-RP-O501 and Salama'], ['dnv', 'DNV-RP-O501'], ['salama', 'Salama (2000)'], ['oka', 'Oka (2005)'], ['finnie', 'Finnie']]) },
    { key: 'erosionMult', label: 'Erosion-model multiplier', unit: '–', value: 1, min: 0.05, max: 20 },
    { key: 'eroVth', label: 'Erosion threshold velocity', unit: 'm/s', value: 0, min: 0, max: 30, help: 'Impacts slower than this cause no wear; the excess velocity enters the erosion models. 0 = no threshold.' },
    { key: 'finnieC', label: 'Finnie model: share of the ideal cutting volume', unit: '–', value: 0.5, min: 0.01, max: 1, help: 'Finnie took one half of the idealised volume as actually removed.' },
    { key: 'restitution', label: 'Particle–wall restitution (normal)', unit: '–', value: 0.8, min: 0, max: 1, help: 'For the tracked-particle model of the worst bend; the tangential value is taken 0.1 higher (at most 1).' },
    { key: 'geomFactor', label: 'Geometry factor of the fittings', unit: '–', value: 1, min: 1, max: 4 },
    { key: 'c14e', label: 'API RP 14E C-factor', unit: 'field units', value: 100, min: 50, max: 300 },
  ] },
  { group: 'Blockage and plugs', tab: 'inputs', fields: [
    { key: 'blockage', label: 'Flow area lost to deposits', unit: 'fraction', value: 0, min: 0, max: 1 },
    { key: 'deposit0', label: 'Initial deposit thickness on the wall', unit: 'mm', value: 0, min: 0, max: 500, help: 'Uniform layer present at the start (wax, scale, hydrate film). It narrows the bore, starts the deposit-growth clock and adds mass to the span.' },
    { key: 'depositDensity', label: 'Deposit density', unit: 'kg/m³', value: 900, min: 100, max: 4000 },
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
    { key: 'rvSamples', label: 'Measured samples for the distributions', type: 'table', columns: [{ key: 'variable', label: 'Variable: wall, yield, pressure or rate', type: 'text' }, { key: 'value', label: 'Value (mm, MPa, bar, mm/y)' }], value: [], help: 'Five or more values of a variable are fitted by maximum likelihood and replace its mean and scatter in the limit states.' },
  ] },
  { group: 'Risk and consequence', tab: 'setup', fields: [
    { key: 'oilPrice', label: 'Oil price', unit: '$/bbl', value: BASE.oilPrice, min: 1, max: 300 },
    { key: 'repairCost', label: 'Repair or intervention cost', unit: 'M$', value: 15, min: 0, max: 2000 },
    { key: 'downtimeDays', label: 'Downtime after a loss of containment', unit: 'd', value: 45, min: 0, max: 730 },
    { key: 'envCost', label: 'Environmental cost of an uncontrolled release', unit: 'M$', value: 40, min: 0, max: 20000 },
    { key: 'safetyCost', label: 'Safety cost of an ignited release', unit: 'M$', value: 100, min: 0, max: 20000 },
    { key: 'pDetect', label: 'Probability that a leak is detected', unit: '–', value: 0.9, min: 0, max: 1 },
    { key: 'pIsolate', label: 'Probability that isolation succeeds', unit: '–', value: 0.95, min: 0, max: 1 },
    { key: 'pIgnite', label: 'Probability of ignition', unit: '–', value: 0.05, min: 0, max: 1 },
    { key: 'pfdProtect', label: 'Overpressure protection, failure on demand', unit: '–', value: 0.01, min: 0, max: 1 },
    { key: 'extFreq', label: 'External damage frequency', unit: '1/y', value: 1e-4, min: 0, max: 1 },
    { key: 'inspInterval', label: 'Planned inspection interval', unit: 'y', value: 0, min: 0, max: 30, help: '0 = risk-based: the longest interval that keeps the annual failure probability below the target of the safety class through the design life.' },
    { key: 'pod', label: 'Probability of detecting a 10 % deep wall loss', unit: '–', value: 0.9, min: 0, max: 0.999, help: 'Defines the detection curve POD(d) = 1 − exp(−d/λ) of the inspection tool.' },
    { key: 'sizingMm', label: 'Depth-sizing error of the inspection (standard deviation)', unit: 'mm', value: 0, min: 0, max: 5, help: '0 = ± 10 % of the wall thickness at 90 % certainty, the typical value quoted for high-resolution magnetic-flux-leakage tools.' },
    { key: 'repairFrac', label: 'Repair threshold: measured wall loss / corrosion allowance', unit: '–', value: 1, min: 0.1, max: 3, help: 'A section whose measured loss exceeds this is repaired or replaced before the next interval.' },
    { key: 'failRate', label: 'Generic loss-of-containment frequency', unit: 'per 1000 km·y', value: 0.5, min: 0, max: 50, help: 'Benchmark from public failure statistics, compared with the predicted failure frequency of this line. 0.5 = subsea well-stream pipelines in open sea (North Sea experience); about 0.05 for processed-fluid subsea lines up to 24 inch and 0.1–0.3 for onshore gas transmission lines.' },
  ] },
  { group: 'Measurements for comparison', tab: 'setup', help: 'Field or test measurements listed here are compared with the model in a table with error metrics. Kind: pressure (bara at x), temperature (°C at x), strain (hoop µε at x), wall (mm at x; condition = age in years), corrosion (mm/y at x), erosion (mm/y at x), force (kN peak at the bend near x), pulsation (bar), frequency (Hz), velocity (mm/s rms on the span), acceleration (m/s² rms), stress (MPa dynamic range on the span), fatigue (cycles to failure; condition = stress range MPa), tensile (MPa; condition = strain %), burst (MPa; x = defect depth mm, condition = defect length mm), collapse (MPa; condition = out-of-roundness %), failure (1 = a failure occurred in that year, 0 = none; condition = age in years).', fields: [
    { key: 'obs', label: 'Measurements', type: 'table', columns: [{ key: 'kind', label: 'Kind', type: 'text' }, { key: 'x', label: 'Location or first condition' }, { key: 'cond', label: 'Condition' }, { key: 'value', label: 'Measured value' }], value: [] },
  ] },
  { group: 'Additional solvers', tab: 'setup', help: 'Resolved coupled solutions that take seconds rather than milliseconds. They run only when selected, on coarse default grids; refine them on the mesh tab. Production resolution (three-dimensional flow, shell or solid models of whole spools) goes through the external-solvers page.', fields: [
    { key: 'task', label: 'Additional solver', type: 'select', value: 'none', options: TASKS },
    { key: 'fsiLen', label: 'Pipe length of the pressure-pulse model', unit: 'm', value: 20, min: 1, max: 500, showIf: (v) => v.task === 'tubeFsi' },
    { key: 'fsiPulse', label: 'Pulse amplitude at the inlet', unit: 'bar', value: 1, min: 0.001, max: 200, showIf: (v) => v.task === 'tubeFsi' },
    { key: 'fsiBulk', label: 'Bulk modulus of the liquid', unit: 'MPa', value: 1500, min: 100, max: 4000, showIf: (v) => v.task === 'tubeFsi', help: 'About 2200 MPa for water, 1000–1700 MPa for crude oil.' },
    { key: 'taskRD', label: 'Bend radius / bore of the modelled bend', unit: '–', value: 1.5, min: 0, max: 20, showIf: (v) => ['bendCfd', 'bendFsi', 'solid'].includes(v.task), help: '0 = the first bend of the table. The bend angle is that of the span group.' },
    { key: 'taskLegs', label: 'Straight leg each side of the bend', unit: 'diameters', value: 3, min: 1, max: 20, showIf: (v) => ['bendCfd', 'bendFsi', 'solid'].includes(v.task) },
  ] },
  { group: 'Discretisation', tab: 'mesh', fields: [
    { key: 'fsiNx', label: 'Coupled pulse: cells along the pipe', unit: '–', value: 60, min: 20, max: 400 },
    { key: 'fsiNr', label: 'Coupled pulse: cells across the radius', unit: '–', value: 6, min: 3, max: 32 },
    { key: 'cfdNs', label: 'Bend flow: cells along the centreline', unit: '–', value: 48, min: 24, max: 320 },
    { key: 'cfdNn', label: 'Bend flow: cells across the bore', unit: '–', value: 8, min: 4, max: 48 },
    { key: 'hexNc', label: 'Solid bend: elements around half the circumference', unit: '–', value: 10, min: 6, max: 40 },
    { key: 'hexNt', label: 'Solid bend: elements through the wall', unit: '–', value: 2, min: 1, max: 4 },
    { key: 'hexNa', label: 'Solid bend: elements along the arc', unit: '–', value: 10, min: 4, max: 60 },
    { key: 'nElem', label: 'Beam elements on the span', unit: '–', value: 24, min: 4, max: 200 },
    { key: 'nModes', label: 'Modes kept', unit: '–', value: 6, min: 1, max: 12 },
    { key: 'stepsPerCycle', label: 'Time steps per first-mode period', unit: '–', value: 40, min: 8, max: 400 },
    { key: 'nSlugs', label: 'Slugs simulated', unit: '–', value: 12, min: 3, max: 60 },
    { key: 'nMC', label: 'Monte Carlo samples', unit: '–', value: 20000, min: 500, max: 400000 },
    { key: 'feOrder', label: 'Order of the wall elements', type: 'select', value: 'quadratic', options: opt([['quadratic', 'Quadratic (6-node triangles)'], ['linear', 'Linear (3-node triangles)']]) },
    { key: 'feNr', label: 'Wall elements through the thickness', unit: '–', value: 2, min: 1, max: 12 },
    { key: 'feNz', label: 'Wall elements along the pipe / around the ring', unit: '–', value: 16, min: 8, max: 160 },
    { key: 'shellN', label: 'Shell elements along the meridian', unit: '–', value: 60, min: 10, max: 400 },
    { key: 'nPart', label: 'Tracked sand particles', unit: '–', value: 40, min: 10, max: 2000 },
    { key: 'bendNs', label: 'Flow-field cells along the bend', unit: '–', value: 32, min: 16, max: 300 },
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
    feOrder: p.feOrder === 'linear' ? 1 : 2, feNr: clamp(Math.round(p.feNr), 1, 12), feNz: clamp(Math.round(p.feNz), 6, 160), shellN: clamp(Math.round(p.shellN), 8, 400), nPart: clamp(Math.round(p.nPart), 8, 2000), bendNs: clamp(Math.round(p.bendNs), 12, 300),
    sc: DNV_SC[p.safetyClass] ? p.safetyClass : 'medium', code: ['dnv', 'b314', 'b318'].includes(p.code) ? p.code : 'dnv', snClass: SN_CURVES[p.snClass] ? p.snClass : 'F1', snEnv: ['air', 'cp', 'free'].includes(p.snEnv) ? p.snEnv : 'cp' });
  p.bends = (Array.isArray(p.bends) ? p.bends : []).map((b) => ({ x: nz(b?.x, 0), angle: clamp(nz(b?.angle, 90), 0, 180), radius: Math.max(nz(b?.radius, 0), 0) })).filter((b) => b.angle > 0).slice(0, 40);
  // material condition: the strength used everywhere, the weld-misalignment stress concentration and the fabrication factor
  const kAge = 1 - clamp(nz(p.strengthLoss, 0), 0, 60) / 100; p.smysNew = p.smys; p.smys *= kAge; p.smts *= kAge; p.S *= kAge; p.Su *= kAge;
  p.scfWeld = 1 + ((3 * Math.max(nz(p.hiLoMm, 0), 0)) / p.wtMm) * Math.exp(-((p.D / p.t) ** -0.5)); p.scfIn = p.scf; p.scf = p.scf * p.scfWeld;
  p.alphaFab = ALPHA_FAB[p.fabrication] ?? 1;
  // girth-weld class from the hi-lo tolerance (single-sided weld root) and the counterpart without the tolerance control (3 mm: class F1)
  const scfOf = (d) => 1 + ((3 * Math.max(d, 0)) / p.wtMm) * Math.exp(-((p.D / p.t) ** -0.5)), autoWeld = p.weldBasis !== 'manual';
  if (autoWeld) p.snClass = p.hiLoMm <= 1 ? 'E' : p.hiLoMm <= 2 ? 'F' : 'F1';
  p.weldAuto = autoWeld; p.unmanagedWeld = autoWeld ? { snClass: 'F1', scf: p.scfIn * scfOf(Math.max(p.hiLoMm, 3)) } : { snClass: p.snClass, scf: p.scf };
  p.fatMargin = clamp(nz(p.fatMargin, 2), 1, 20); p.spanLengthEntered = p.spanLength;
  p.lineSupports = (Array.isArray(p.lineSupports) ? p.lineSupports : []).map((s) => ({ x: nz(s?.x, NaN), type: String(s?.type || 'support').slice(0, 60), k: Math.max(nz(s?.k, 0), 0) })).filter((s) => isNum(s.x)).slice(0, 200);
  { const a = p.spanX - p.spanLength / 2, b = p.spanX + p.spanLength / 2, cuts = [a, ...p.lineSupports.map((s) => s.x).filter((x) => x > a + 0.5 && x < b - 0.5).sort((u, w) => u - w), b]; let best = 0; for (let i = 1; i < cuts.length; i++) if (cuts[i] - cuts[i - 1] > cuts[best + 1] - cuts[best]) best = i - 1; p.spanSplit = cuts.length > 2; if (p.spanSplit) { p.spanLength = cuts[best + 1] - cuts[best]; p.spanX = 0.5 * (cuts[best] + cuts[best + 1]); } }
  // corrosion-control basis
  p.corrControl = ['availability', 'efficiency', 'none'].includes(p.corrControl) ? p.corrControl : 'availability'; p.avail = clamp(p.inhibAvail, 0, 99.9) / 100;
  if (p.corrControl === 'availability') p.eta = p.avail; else if (p.corrControl === 'none') p.eta = 0;
  p.phStabOn = p.phStab === 'on' || (p.phStab !== 'off' && p.glycolWt >= 20 && !(p.bicarb > 0)); p.brineWt = clamp(nz(p.brineWt, 0), 0, 26);
  p.task = ['tubeFsi', 'bendCfd', 'bendFsi', 'solid'].includes(p.task) ? p.task : 'none';
  for (const [k, lo, hi] of [['fsiNx', 20, 400], ['fsiNr', 3, 32], ['cfdNs', 24, 320], ['cfdNn', 4, 48], ['hexNc', 6, 40], ['hexNt', 1, 4], ['hexNa', 4, 60]]) p[k] = clamp(Math.round(p[k]), lo, hi);
  p.midSupports = (Array.isArray(p.midSupports) ? p.midSupports : []).map((s) => { const type = String(s?.type || 'support').trim().toLowerCase(), kind = /^anc/.test(type) ? 'anchor' : /^gui/.test(type) ? 'guide' : /^(mas|fla|val|con)/.test(type) ? 'mass' : 'support', k = Math.max(nz(s?.k, 0), 0) * 1000; return { x: nz(s?.x, NaN), type: kind, k: kind === 'guide' ? (k > 0 ? k : Math.max(p.spanKt, 1) * 1000) : kind === 'support' ? k : 0, mass: Math.max(nz(s?.mass, 0), 0) }; }).filter((s) => s.x > 0 && s.x < p.spanLength).slice(0, 12);
  p.equipment = (Array.isArray(p.equipment) ? p.equipment : []).map((e, i) => ({ name: String(e?.name || `Item ${i + 1}`).slice(0, 60), x: nz(e?.x, 0), size: Math.max(nz(e?.size, 0), 0), arm: Math.max(nz(e?.arm, 0), 0), fAllow: Math.max(nz(e?.fAllow, 0), 0), mAllow: Math.max(nz(e?.mAllow, 0), 0) })).slice(0, 40);
  p.obs = (Array.isArray(p.obs) ? p.obs : []).map((r) => ({ kind: String(r?.kind || '').trim().toLowerCase(), x: nz(r?.x, 0), cond: nz(r?.cond, 0), value: nz(r?.value, NaN) })).filter((r) => r.kind && isNum(r.value)).slice(0, 400);
  p.rvSamples = (Array.isArray(p.rvSamples) ? p.rvSamples : []).map((r) => ({ variable: String(r?.variable || '').trim().toLowerCase(), value: nz(r?.value, NaN) })).filter((r) => r.variable && isNum(r.value) && r.value > 0).slice(0, 2000);
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
  const inhibited = (T) => (p.inhibRate > 0 ? p.inhibRate : T <= 70 ? 0.05 : T >= 120 ? 0.2 : 0.05 + (0.15 * (T - 70)) / 50), oldModel = p.corrModel === '1991' || p.corrModel === '1995';
  const corr = N.map((nd) => {
    const strip = clamp((nd.tauW - p.tauCrit) / Math.max(p.tauCrit, 1e-6), 0, 1), eta = p.eta * (1 - strip), dw = deWaardMilliams({ T: nd.T, pCO2: yCO2 * nd.P, P: nd.P, model: p.corrModel === '1991' ? '1991' : '1995', U: nd.vL, d: p.ID, pH: p.phAct, glycolWt: p.glycolWt, mult: p.corrMult });
    const oilWet = p.wetModel === 'dewaard' && (!hasWater || (nd.wcut < 0.3 && nd.vL > 1));
    const cc = (pH, stab) => co2Corrosion({ model: p.corrModel, T: nd.T, pCO2: yCO2 * nd.P, P: nd.P, pH, bicarb: p.bicarb, salt: p.brineWt, S: norsokShear(nd.rhoNS, nd.muNS, nd.vm, p.ID, p.roughUm * 1e-6), U: Math.max(nd.vL, 0.01), d: p.ID, rho: nd.rhoL, mu: nd.muL, glycolWt: stab ? 0 : p.glycolWt, mult: p.corrMult });
    // unmanaged: no inhibitor, no pH stabilisation, no cladding; managed: the control basis applied to the (stabilised) bare rate
    const un = cc(p.phAct, false), st1 = p.phStabOn ? cc(Math.max(p.phStabTarget, un.pH), true) : un, bare = oilWet ? 0 : oldModel ? dw.rate : un.rate, base = oilWet ? 0 : oldModel ? dw.rate : Math.min(st1.rate, un.rate), cri = Math.min(inhibited(nd.T), base);
    const rate = p.corrControl === 'availability' ? eta * cri + (1 - eta) * base : base * (1 - eta);
    return { rate, steel: rate, bare, base, eta, strip, cri, oilWet, dw, un, clad: false, model: oldModel ? `de Waard (${p.corrModel})` : un.model, pH: oldModel ? dw.pH : (p.phStabOn ? st1 : un).pH };
  });
  // corrosion-resistant alloy: given length, or wherever the controlled carbon-steel rate exceeds the target
  for (let i = 0; i < N.length; i++) { const c = corr[i]; c.clad = p.cladMode === 'length' ? p.cladLength > 0 && N[i].x >= p.cladFrom - 1e-6 && N[i].x <= p.cladFrom + p.cladLength + 1e-6 : p.cladMode === 'auto' ? c.steel > p.corrTarget * 1.0000001 : false; if (c.clad) c.rate = 0; }
  let cladLen = 0; for (let i = 1; i < N.length; i++) if (corr[i].clad && corr[i - 1].clad) cladLen += N[i].x - N[i - 1].x; else if (corr[i].clad !== corr[i - 1].clad) cladLen += 0.5 * (N[i].x - N[i - 1].x);
  const eroAt = (nd, rOverD) => {
    const b = { mp, U: Math.max(nd.vm - p.eroVth, 0), D: p.ID, dp, rhoM: nd.rhoNS, muM: nd.muNS, rhoP: p.sandDensity, rOverD, gf: p.geomFactor, mult: p.erosionMult, hv, flowStress: hv * 1e9, finnieC: p.finnieC }, r = {};
    for (const m of ['dnv', 'salama', 'oka', 'finnie']) r[m] = erosionRate({ ...b, model: m }).rate;
    r.gov = p.erosionModel === 'governing' ? Math.max(r.dnv, r.salama) : r[p.erosionModel] ?? r.dnv; r.model = p.erosionModel === 'governing' ? (r.dnv >= r.salama ? 'DNV-RP-O501' : 'Salama') : { dnv: 'DNV-RP-O501', salama: 'Salama', oka: 'Oka', finnie: 'Finnie' }[p.erosionModel] || 'DNV-RP-O501';
    return r;
  };
  const ero = N.map((nd) => { const ve = erosionalVelocity(nd.rhoNS, p.c14e); return { straight: erosionRate({ model: 'dnv', geometry: 'straight', mp, U: Math.max(nd.vm - p.eroVth, 0), D: p.ID, mult: p.erosionMult }).rate, ve, ratio: nd.vm / ve }; });
  const fast = N.reduce((a, b) => (b.vm > a.vm ? b : a)), list = p.bends.length ? p.bends.map((b, k) => ({ ...b, name: `Bend ${k + 1}`, nd: st.at(b.x), rOverD: b.radius > 0 ? Math.max(b.radius / p.ID, 0.5) : 1.5 })) : [{ x: fast.x, angle: 90, radius: 1.5 * p.ID, name: 'Generic 1.5D elbow at the fastest point', nd: fast, rOverD: 1.5 }];
  const bends = list.map((b) => ({ ...b, ero: eroAt(b.nd, b.rOverD), corr: corr[b.nd.i].rate }));
  let gov = { rate: -1 };
  N.forEach((nd, i) => { const r = corr[i].rate + ero[i].straight; if (r > gov.rate) gov = { rate: r, x: nd.x, corr: corr[i].rate, ero: ero[i].straight, where: 'straight pipe' }; });
  for (const b of bends) { const r = b.corr + b.ero.gov; if (r > gov.rate) gov = { rate: r, x: b.nd.x, corr: b.corr, ero: b.ero.gov, where: b.name }; }
  const iU = corr.reduce((k, c, i) => (c.bare > corr[k].bare ? i : k), 0), iC0 = corr.reduce((k, c, i) => (c.rate > corr[k].rate ? i : k), 0), iC = corr[iC0].rate > 0 ? iC0 : iU, eroMax = Math.max(...bends.map((b) => b.ero.gov), ...ero.map((e) => e.straight)), bWorst = bends.reduce((a, b) => (b.ero.gov > a.ero.gov ? b : a));
  const iV = ero.reduce((k, e, i) => (e.ratio > ero[k].ratio ? i : k), 0);
  const cladX = N.filter((_, i) => corr[i].clad).map((n) => n.x), steelMax = Math.max(...corr.map((c) => c.steel));
  return { corr, ero, bends, gov, iC, iU, corrMax: corr[iC0].rate, corrUnmanaged: corr[iU].bare, steelMax, cladLen, cladFrom: cladX.length ? Math.min(...cladX) : null, cladTo: cladX.length ? Math.max(...cladX) : null, inhibited, eroMax, bWorst, iV, erosionalRatio: ero[iV].ratio, yCO2, hasWater, mp, dp };
}

/** Stress, code utilisation, collapse and global buckling along the line. */
function stressAlong(p, st, wall) {
  const N = st.nodes, dnv = p.code === 'dnv', tNow = wall.tNow, tDes = Math.max(p.t - p.CA / 1000 - (dnv ? p.tFabMm / 1000 : 0), 0.05 * p.t), gLB = DNV_GM * DNV_SC[p.sc].lb, H = p.layTension * 1000;
  const head = (nd) => (p.pRefLoc === 'uniform' ? 0 : p.rhoContent * G * (st.zRef - nd.z)), pLocal = (nd, inc = 1, pRef = p.designPressure) => Math.max(pRef * BAR * inc + head(nd), PATM);
  const fmtOf = (nd, t, T) => pressureDesign({ code: p.code, D: p.D, t, smys: p.S, smts: p.Su, T, factor: nd.riser ? p.riserFactor : p.designFactor, safetyClass: p.sc });
  const fyCold = fmtOf(N[0], tNow, p.tInstall).fy, colOf = (t) => collapsePressure({ D: p.D, t, E: p.E, nu: p.nu, fy: fyCold, ovality: p.f0, alphaFab: p.alphaFab }), col = colOf(tNow), colNew = colOf(p.t), pMin = Math.max(p.pMinShut * BAR, PATM);
  const rows = N.map((nd) => {
    const fo = fmtOf(nd, tNow, nd.T), fd = fmtOf(nd, tDes, p.designTemp), fn = fmtOf(nd, tNow, p.designTemp), restrained = p.restraint === 'restrained' || (p.restraint === 'auto' && !nd.riser), common = { D: p.D, t: tNow, pe: nd.pe, E: p.E, nu: p.nu, alpha: p.alphaT, restrained, axial: restrained ? H : 0 };
    const sOp = pipeStress({ ...common, pi: nd.pi, dT: nd.T - p.tInstall }), pD = pLocal(nd), sDes = pipeStress({ ...common, pi: pD, dT: p.designTemp - p.tInstall }), pLi = pLocal(nd, fd.incidental);
    const seff = restrained ? H - (nd.pi - PATM) * p.Ai * (1 - 2 * p.nu) - p.As * p.E * p.alphaT * (nd.T - p.tInstall) : 0, pex = Math.max(nd.pe - PATM, 0), pexOp = Math.max(nd.pe - pMin, 0);
    return { nd, restrained, sOp, sDes, pD, seff, hoopOp: Math.max(nd.pi - nd.pe, 0) / fo.allowDp, hoopDes: Math.max(pLi - nd.pe, 0) / fd.allowDp, longOp: Math.abs(sOp.longitudinal) / fo.longAllow, longDes: Math.abs(sDes.longitudinal) / fn.longAllow, vmOp: sOp.vonMises / fo.combAllow, vmDes: sDes.vonMises / fn.combAllow,
      collapse: Math.max((pex * gLB) / colNew.pc, (pexOp * gLB) / col.pc), collapseInstall: (pex * gLB) / colNew.pc, collapseOp: (pexOp * gLB) / col.pc, propagation: (pex * gLB) / col.ppr, tReq: fd.tReq(pLi - nd.pe), allowRef: (t) => (fmtOf(nd, t, p.designTemp).allowDp + nd.pe - head(nd)) / fd.incidental, basis: fd.basis };
  });
  const arg = (f) => rows.reduce((k, r, i) => (f(r) > f(rows[k]) ? i : k), 0), iH = arg((r) => Math.max(r.hoopDes, r.hoopOp)), iV = arg((r) => Math.max(r.vmDes, r.vmOp)), iC = arg((r) => r.collapse), iS = arg((r) => -r.seff);
  const mawpOf = (t) => Math.max(Math.min(...rows.map((r) => r.allowRef(t))) / BAR, 0), tEol = Math.max(tNow - (wall.rate * Math.max(p.designLife - p.age, 0)) / 1000, 0.02 * p.t);
  // global buckling of the restrained flowline: Hobbs forces and hold-down against upheaval
  const rS = rows[iS], sub = rS.nd.depth > 0, wDry = (p.rhoSteel * p.As + p.coatDensity * (Math.PI / 4) * (p.Dh ** 2 - p.D ** 2) + rS.nd.rhoM * p.Ai) * G, wSub = wDry - (sub ? RHO_SW * (Math.PI / 4) * p.Dh ** 2 * G : 0), w = Math.max(wSub, 1), EI = p.E * p.Isteel, comp = Math.max(-rS.seff, 0);
  const hb = hobbs({ EI, EA: p.E * p.As, w, muA: p.muAxial, muL: p.muLateral }), up = upheavalDownload({ EI, P: comp, delta: p.imperfection, w0: w }), soil = (Hc) => p.soilGamma * 1000 * Hc * p.Dh * (1 + (0.5 * Hc) / p.Dh), resist = w + soil(p.coverDepth);
  const coverReq = up.wReq <= w ? 0 : soil(10) + w < up.wReq ? 10 : brent((Hc) => w + soil(Hc) - up.wReq, 0, 10, 1e-9);
  return { rows, iH, iV, iC, iS, col, colNew, colOf, fyCold, pMin, tDes, tEol, pLocal, head, fmtOf, gLB, mawp: mawpOf(tNow), mawpEol: mawpOf(tEol), mawpOf, minWall: (Math.max(...rows.map((r) => r.tReq)) + p.CA / 1000 + (dnv ? p.tFabMm / 1000 : 0)) * 1000,
    hoopUtil: Math.max(rows[iH].hoopDes, rows[iH].hoopOp), vmUtil: Math.max(rows[iV].vmDes, rows[iV].vmOp), longUtil: Math.max(...rows.map((r) => Math.max(r.longOp, r.longDes))), collapseUtil: rows[iC].collapse, propUtil: rows[iC].propagation,
    buckling: { hb, up, wSub, w, comp, resist, coverReq, hobbsUtil: comp / hb.critical, palmerUtil: up.wReq / resist, buried: p.coverDepth > 0, x: rS.nd.x, EI, EA: p.E * p.As, nd: rS.nd } };
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
  const cf = ds.on ? p.slugCf : 1, f90 = ds.on ? cf * bendForce({ rho: ds.rhoS, A, v: ds.v, angle: 90, dlf: p.dlf }).force : 0, steady90 = Math.max(...st.nodes.map((n) => bendForce({ rho: n.rhoNS, A, v: n.vm, angle: 90 }).force));
  const dpAcc = ds.on ? nd.rhoL * ds.hls * (ds.v - vmS) ** 2 * Math.max(ds.hls / Math.max(ds.hlf, 1e-3) - 1, 0) : 0, fD = frictionFactor(Math.max((ds.rhoS * vmS * p.ID) / nd.muL, 10), (p.roughUm * 1e-6) / p.ID), dpFric = ds.on ? ((fD * ds.rhoS * vmS * vmS) / (2 * p.ID)) * ds.len : 0;
  const bends = (p.bends.length ? p.bends : [{ x: nd.x, angle: 90, radius: 1.5 * p.ID }]).map((b, i) => {
    const n = st.at(b.x), rho = !ds.on ? n.rhoNS : ds.fixedRho ? ds.rhoS : fromKernel ? n.slug.rhoS : n.rhoL * ds.hls + n.rhoG * (1 - ds.hls), v = !ds.on ? n.vm : fromKernel ? n.slug.v : (ds.v * n.vm) / Math.max(nd.vm, 1e-6);
    const bf = bendForce({ rho, A, v, angle: b.angle, dlf: ds.on ? p.dlf : 1, p: n.pi - n.pe, radius: b.radius }), st0 = bendForce({ rho: n.rhoNS, A, v: n.vm, angle: b.angle }), fiv = fivScreen({ rho: n.rhoNS, v: n.vm, Dmm: p.D * 1000, tmm: p.wtMm, support: p.pipingSupport });
    bf.force *= cf;
    return { name: `Bend ${i + 1}`, x: n.x, angle: b.angle, radius: b.radius, rho, v, impulse: ds.on ? ((bf.force / p.dlf) * ds.len) / Math.max(v, 1e-6) : 0, force: bf.force, steady: st0.force, pressure: bf.pressureForce, centrifugal: bf.centrifugal, fiv };
  });
  const train = slugTrain(p, ds);
  const Fs = cf * ds.rhoS * A * ds.v ** 2 * Math.SQRT2, Ff = cf * ds.rhoF * A * ds.v ** 2 * Math.SQRT2, tr = Math.max(p.ID / ds.v, 1e-3), ht = [0], hf = [ds.on ? Ff : steady90];
  for (const s of train.slice(0, 12)) { const b = Math.max(s.len / s.v, 1.5 * tr); ht.push(s.t0 + 0.2 * s.period, s.t0 + 0.2 * s.period + tr, s.t0 + 0.2 * s.period + b, s.t0 + 0.2 * s.period + b + tr); hf.push(Ff, Fs, Fs, Ff); }
  if (!ds.on) { ht.push(60); hf.push(steady90); }
  const fPIn = ctx.outputs?.flow?.pInAmplitude;
  return { f90, steady90, peak90: ds.on ? f90 : steady90, dpAcc, dpFric, pulsation: (dpAcc + dpFric) / 2 / BAR, runForce: (dpAcc + dpFric) * A, bends, train, history: { t: ht, f: hf.map((x) => x / 1000) }, flowAmplitude: isNum(fPIn) ? fPIn : null, Fs, Ff };
}

/** Free-span or jumper dynamics: modes, static sag, slug-train response, vortex-induced vibration and the fatigue spectrum. */
function spanDynamics(p, st, ds, sl, sa, ctx) {
  const nd = st.at(p.spanX), sub = nd.depth > 0, L = p.spanLength, areaO = (Math.PI / 4) * p.Dh ** 2, mSteel = p.rhoSteel * p.As, mCoat = p.coatDensity * (Math.PI / 4) * (p.Dh ** 2 - p.D ** 2), dep0 = clamp(nz(p.deposit0, 0) / 1000, 0, 0.49 * p.ID), mDep = p.depositDensity * (Math.PI / 4) * (p.ID ** 2 - (p.ID - 2 * dep0) ** 2), mCont = nd.rhoM * (p.Ai - mDep / p.depositDensity) + mDep, mAdd = sub ? p.addedMass * RHO_SW * areaO : 0, me = mSteel + mCoat + mCont + mAdd;
  const wSub = (mSteel + mCoat + mCont) * G - (sub ? RHO_SW * areaO * G : 0), EI = p.E * p.Isteel, seffRaw = sa.rows[nd.i].seff, seff = sa.buckling.buried ? seffRaw : Math.max(seffRaw, -sa.buckling.hb.critical), Naxial = p.spanAxial ? seff : 0, bendX = p.bendAt > 0 ? clamp(p.bendAt, 0.02, 0.98) * L : null, ro = p.D / 2;
  const model = beamModel({ L, EI, m: me, n: p.nElem, ends: endsOf(p.spanEnds), kT: [p.spanKt * 1000, p.spanKt * 1000], kR: [p.spanKr * 1000, p.spanKr * 1000], supports: p.midSupports, nodesAt: bendX ? [bendX] : [], N: Naxial, kGA: p.timoshenko ? (0.5 * p.E * p.As) / (2 * (1 + p.nu)) : 0, rhoI: p.timoshenko ? p.rhoSteel * p.Isteel : 0 });
  const modes = beamModes(model, p.nModes), f1 = modes.f[0], stat = beamStatic(model, { q: Math.abs(wSub) }), sigStatic = (stat.maxMoment * ro) / p.Isteel / MPA, euler = eulerLoad(EI, L, END_K[p.spanEnds] || 1);
  // slug-train response
  const sg = ds.source === 'kernel slug closures' ? { on: ds.on && (nd.slugLike || p.slugMode === 'on') && nd.slug.freq > 0, freq: nd.slug.freq, len: Math.max(nd.slug.len, 0.1), v: Math.max(nd.slug.v, 0.05), hls: nd.slug.hls, hlf: nd.slug.hlf } : { on: ds.on, freq: ds.freq, len: ds.len, v: ds.v, hls: ds.hls, hlf: ds.hlf };
  const train = slugTrain(p, sg), hasSlug = sg.on && train.length > 0, v = sg.v, rhoS = ds.fixedRho ? ds.rhoS : nd.rhoL * sg.hls + nd.rhoG * (1 - sg.hls), rhoF = Math.min(nd.rhoL * sg.hlf + nd.rhoG * (1 - sg.hlf), rhoS), dw = (rhoS - rhoF) * p.Ai * G;
  let resp = null, cmp = null, cycles = [], perYear = 0, note = '', dt = 1 / (f1 * p.spc), simSlugs = 0;
  // conveyed fluid: critical velocity of the span and the first frequency with the flow (gyroscopic eigenproblem)
  const mi = modalIntegrals(model, modes), mConv = nd.rhoM * p.Ai, vConv = hasSlug ? v : nd.vm, crit = fluidCritical(model, modes, mConv, mi), fFlow = fluidFrequency(model, modes, mConv, vConv, mi);
  if (hasSlug) {
    dt = Math.min(dt, L / v / 25);
    const decay = Math.min(4 / (p.zeta * modes.omega[0]), 40 / f1), slugs = []; let t0 = 0.5 / f1, real = 0;
    for (const s of train) { const gap = Math.min(s.period, (L + s.len) / v + decay); if ((t0 + gap) / dt > 4000 && slugs.length >= 2) break; slugs.push({ t0, len: s.len, v }); t0 += gap; real += s.period; }
    simSlugs = slugs.length; if (simSlugs < train.length) note = `Only ${simSlugs} of ${train.length} slugs were simulated to keep the response history below 4,000 steps.`;
    const common = { model, modes, zeta: p.zeta, slugs, dw, v, mf: mConv, dM: Math.max(rhoS - nd.rhoM, 0) * p.Ai, bend: bendX ? { x: bendX, dF: p.slugCf * (rhoS - rhoF) * p.Ai * v * v * 2 * Math.sin((p.bendAngle * Math.PI) / 360) } : null, dt, tEnd: t0, ro, E: p.E, ramp: Math.max(p.ID, 2 * v * dt), mi };
    if (p.fsi === 'one') resp = fsiResponse({ ...common, coupled: false });
    else { // two-way run for the whole train; the one-way comparison covers the first passages only
      resp = fsiResponse({ ...common, coupled: true }); const nWin = Math.min(resp.steps, 1500); let two = 0; for (let i = 0; i <= nWin; i++) two = Math.max(two, Math.abs(resp.sigma[i]));
      const one = fsiResponse({ ...common, coupled: false, tEnd: nWin * dt }); cmp = { two, one: one.sigmaMax, window: nWin * dt };
    }
    cycles = rainflow(resp.sigma.map((s) => s / MPA)); perYear = YEAR / real;
  }
  const snOpt = { scf: p.scf, t: p.wtMm, shift: p.snShift }, dmgOf = (cyc) => snDamage(cyc, p.snClass, p.snEnv, snOpt).damage, spectrum = [], dSlug = dmgOf(cycles) * perYear;
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
  const vivVel = 2 * Math.PI * f1 * Math.max(vivWorst.aCF, vivWorst.aIL) * p.Dh * Math.SQRT1_2, velRms = (resp ? resp.velRms : 0) * 1000, vLim = vibrationLimits(f1), vAllow = p.vibAllow > 0 ? p.vibAllow : vLim.concern;
  return { mi, mConv, vConv, crit, fFlow, cmp, twoWay: !!resp && p.fsi !== 'one', vivVel: vivVel * 1000, velRms, velPeak: resp ? resp.velMax * 1000 : vivVel * Math.SQRT2 * 1000, accRms: resp ? resp.accRms : 0, vLim, vAllow, vibUtil: velRms / vAllow, mDep, snOpt,
    nd, sub, L, me, wSub, EI, seff, seffRaw, Naxial, model, modes, f1, stat, sigStatic, euler, eulerUtil: Math.max(-seff, 0) / euler, resp, cycles, perYear, spectrum, dSlug, dViv, dOps, dPuls, dYear, viv, vivWorst, uw, unit, fivRatio: gov.f / f1, excitation: gov, note, dt, simSlugs, dynMax: resp ? resp.sigmaMax / MPA : 0, dispMax: resp ? resp.dispMax : 0, passage: hasSlug ? (L + sg.len) / v : 0, slug: { ...sg, rhoS, rhoF } };
}

// ---- vibration acceptance, material curve, distribution fitting, inspection planning --------------------------------------
/** Screening lines for the vibration velocity of process pipework (mm/s rms) at frequency f (Hz): below `concern` acceptable, above `problem` not acceptable. */
export function vibrationLimits(f) {
  const lf = Math.log10(clamp(f, 1, 300));
  return { concern: 10 ** ((lf + VIB_LINES.concern[0]) / VIB_LINES.concern[1]), problem: 10 ** ((lf + VIB_LINES.problem[0]) / VIB_LINES.problem[1]) };
}
/** Ramberg–Osgood curve ε = σ/E + 0.002·(σ/σy)^n: strain for a stress, and stress for a strain (same units as E and σy). */
export const rambergStrain = (sigma, E, sy, n) => sigma / E + 0.002 * (Math.abs(sigma) / sy) ** n * Math.sign(sigma);
export function rambergStress(eps, E, sy, n) { const e = Math.abs(eps); if (!(e > 0)) return 0; const hi = Math.min(E * e, sy * (e / 0.002) ** (1 / n) * 1.0000001 + 1e-9); return Math.sign(eps) * brent((s) => rambergStrain(s, E, sy, n) - e, 0, hi, 1e-12 * hi); }
/**
 * Maximum-likelihood fit of a distribution to a sample. dist: 'normal' | 'lognormal' | 'weibull' | 'gumbel' (largest values).
 * Returns { dist, n, mean, sd, cov, logLik, params } — mean and sd are those of the fitted distribution.
 */
export function fitDistribution(values, dist = 'normal') {
  const x = values.filter((v) => Number.isFinite(v)), n = x.length; need(n >= 2, 'Fitting a distribution needs at least two values.');
  const m = mean(x), v = sum(x.map((q) => (q - m) ** 2)) / n, done = (mu, sd, logLik, params) => ({ dist, n, mean: mu, sd, cov: mu !== 0 ? sd / Math.abs(mu) : 0, logLik, params });
  if (dist === 'lognormal') { need(x.every((q) => q > 0), 'A lognormal fit needs positive values.'); const l = x.map(Math.log), ml = mean(l), vl = Math.max(sum(l.map((q) => (q - ml) ** 2)) / n, 1e-300); return done(Math.exp(ml + vl / 2), Math.exp(ml + vl / 2) * Math.sqrt(Math.exp(vl) - 1), -sum(l) - (n / 2) * Math.log(2 * Math.PI * vl) - n / 2, { mu: ml, sigma: Math.sqrt(vl) }); }
  if (dist === 'weibull') { // shape from the likelihood equation Σx^k·ln x / Σx^k − 1/k − mean(ln x) = 0
    need(x.every((q) => q > 0), 'A Weibull fit needs positive values.'); const l = x.map(Math.log), ml = mean(l), xm = Math.max(...x), eq = (k) => { let a = 0, b = 0; for (let i = 0; i < n; i++) { const t = (x[i] / xm) ** k; a += t * l[i]; b += t; } return a / b - 1 / k - ml; };
    const k = brent(eq, 0.05, 500, 1e-12), lam = (sum(x.map((q) => q ** k)) / n) ** (1 / k), mu = lam * gammaFn(1 + 1 / k);
    return done(mu, lam * Math.sqrt(Math.max(gammaFn(1 + 2 / k) - gammaFn(1 + 1 / k) ** 2, 0)), n * Math.log(k / lam) + (k - 1) * sum(x.map((q) => Math.log(q / lam))) - sum(x.map((q) => (q / lam) ** k)), { shape: k, scale: lam });
  }
  if (dist === 'gumbel') { // scale from b = mean(x) − Σx·e^(−x/b)/Σe^(−x/b), then location
    const sd0 = Math.sqrt(Math.max(v, 1e-300)), eq = (b) => { let a = 0, c = 0; for (let i = 0; i < n; i++) { const t = Math.exp(-(x[i] - m) / b); a += x[i] * t; c += t; } return m - a / c - b; };
    const b = brent(eq, 0.02 * sd0, 20 * sd0, 1e-13 * sd0), a = -b * Math.log(sum(x.map((q) => Math.exp(-q / b + m / b))) / n) + m, z = x.map((q) => (q - a) / b);
    return done(a + 0.5772156649015329 * b, (Math.PI * b) / Math.sqrt(6), -n * Math.log(b) - sum(z) - sum(z.map((q) => Math.exp(-q))), { location: a, scale: b });
  }
  const sd = Math.sqrt(Math.max(v, 1e-300)); return done(m, sd, -(n / 2) * Math.log(2 * Math.PI * v) - n / 2, { mu: m, sigma: sd });
}
/**
 * Annual failure probability of wall-loss limit states under an inspection-and-repair plan, with Bayesian updating of the wall-loss factor.
 * The factor X of each limit state (lognormal, mean 1) is discretised into equally likely quantiles; for every quantile the wall loss grows as
 * loss0 + X·rate·(time since the last repair). At an inspection the loss is detected with POD(loss), sized with a normal error, and the section is
 * repaired (loss back to zero) when the measured loss exceeds the threshold. The survivors' weights are the posterior of X given "no repair so far".
 * o: { ls: [{ name, pf(loss) (annual: conditional annual failure probability; else probability that the state has been reached), annual, rate (mm/y), loss0 (mm), cov }],
 *      years, interval (y; 0 or Infinity = no inspection), pod(loss), sizingSd (mm), threshold (mm), nq, cap (largest loss, mm) }
 * Returns { t[] (1 … years), annual[] (series system), perLS: [[…]], inspections: [{ t, pRepair, mean, cov (posterior of X of the first limit state given no repair) }], repairs (expected number over the horizon) }.
 */
export function inspectionPlan(o) {
  const N = Math.max(1, Math.round(o.years)), nq = clamp(Math.round(nz(o.nq, 40)), 4, 400), iv = o.interval > 0 && Number.isFinite(o.interval) ? +o.interval : Infinity, sdm = Math.max(nz(o.sizingSd, 0), 1e-9), thr = +o.threshold, capL = pos(o.cap, 1e9), t = Array.from({ length: N }, (_, i) => i + 1), perLS = [], insp = [];
  const times = []; for (let k = 1; k * iv < N - 1e-9; k++) times.push(k * iv);
  let repairs = 0;
  o.ls.forEach((ls, il) => {
    const rv = randomVariable({ dist: 'lognormal', mean: 1, cov: Math.max(nz(ls.cov, 0), 1e-6) }), out = new Array(N).fill(0), post = times.map(() => ({ w: 0, m1: 0, m2: 0, rep: 0 }));
    for (let q = 0; q < nq; q++) {
      const xq = rv.x(PhiInv((q + 0.5) / nq)), w = 1 / nq, coh = [{ m: 1, t0: 0, l0: Math.max(nz(ls.loss0, 0), 0) }], lossAt = (c, tt) => Math.min(c.l0 + xq * ls.rate * Math.max(tt - c.t0, 0), capL); let ki = 0;
      for (let y = 1; y <= N; y++) {
        while (ki < times.length && times[ki] < y - 1e-9) { // inspections that fall before the end of this year
          const tk = times[ki]; let rep = 0;
          for (const c of coh) { const l = lossAt(c, tk), pr = clamp(o.pod(l), 0, 1) * Phi((l - thr) / sdm); rep += c.m * pr; c.m *= 1 - pr; }
          if (rep > 0) coh.push({ m: rep, t0: tk, l0: 0 });
          const never = coh[0].m; post[ki].w += w * never; post[ki].m1 += w * never * xq; post[ki].m2 += w * never * xq * xq; post[ki].rep += w * rep; ki++;
        }
        let h = 0;
        for (const c of coh) { if (!(c.m > 0)) continue; const a = ls.pf(lossAt(c, y)); h += c.m * (ls.annual ? a : Math.max(a - ls.pf(lossAt(c, Math.max(y - 1, c.t0))), 0)); }
        out[y - 1] += w * h;
      }
    }
    perLS.push(out.map((v) => clamp(v, 0, 1)));
    post.forEach((ps, k) => { if (il === 0) { const mu = ps.w > 0 ? ps.m1 / ps.w : 1, va = ps.w > 0 ? Math.max(ps.m2 / ps.w - mu * mu, 0) : 0; insp.push({ t: times[k], pRepair: ps.rep, mean: mu, cov: mu > 0 ? Math.sqrt(va) / mu : 0 }); } repairs += il === 0 ? ps.rep : 0; });
  });
  return { t, annual: t.map((_, i) => 1 - perLS.reduce((s, a) => s * (1 - a[i]), 1)), perLS, inspections: insp, repairs };
}
/** Brier score and reliability-diagram bins of probability forecasts p[] against outcomes o[] (0 or 1). */
export function brierScore(pr, ob, nb = 5) {
  const n = Math.min(pr.length, ob.length); if (!n) return { n: 0, brier: 0, bins: [] };
  let s = 0; const bins = Array.from({ length: nb }, () => ({ n: 0, p: 0, o: 0 }));
  for (let i = 0; i < n; i++) { s += (pr[i] - ob[i]) ** 2; const b = bins[Math.min(nb - 1, Math.floor(clamp(pr[i], 0, 1) * nb))]; b.n++; b.p += pr[i]; b.o += ob[i]; }
  return { n, brier: s / n, bins: bins.filter((b) => b.n).map((b) => ({ n: b.n, forecast: b.p / b.n, observed: b.o / b.n })) };
}

/** Buckle-arrestor design of the line: where a buckle could propagate, how many arrestors confine it and whether the arrestor holds. */
function arrestorDesign(p, st, sa, ctx) {
  const N = st.nodes, col = sa.col, dLim = col.ppr / sa.gLB / (RHO_SW * G), net = ctx.outputs?.net?.buckleArrestors; let len = 0;
  for (let i = 1; i < N.length; i++) if (0.5 * (N[i].depth + N[i - 1].depth) > dLim) len += Math.hypot(N[i].x - N[i - 1].x, N[i].z - N[i - 1].z);
  const required = sa.propUtil > 1 && len > 0, provided = required && p.arrestors !== 'none', spacing = pos(net?.spacing, p.arrSpacing), count = provided ? (net && net.count > 0 ? Math.round(net.count) : Math.ceil(len / spacing) + 1) : 0;
  const x = arrestorCrossover({ D: p.D, t: sa.rows.length ? p.t : p.t, t2: p.arrThick * p.t, L: p.arrLength, fy: sa.fyCold, alphaFab: p.alphaFab }), pe = Math.max(...N.map((n) => n.pe - PATM), 0);
  return { required, provided, depthLimit: dLim, length: len, spacing, count, pX: x.pX, pprBA: x.pprBA, util: (pe * 1.1 * sa.gLB) / x.pX, pe, confined: Math.min(spacing, len), unconfined: len, source: net && (net.spacing > 0 || net.count > 0) ? 'network suite' : 'inputs', firstDepth: isNum(net?.firstDepth) ? net.firstDepth : dLim };
}
/** Lateral-buckling management: feed-in and post-buckle bending strain of the planned buckles against the strain criterion. */
function lateralDesign(p, st, sa, ctx) {
  const bk = sa.buckling, net = ctx.outputs?.net?.buckleInitiators, managed = p.buckleMgmt !== 'none', spacing = pos(net?.spacing, p.initSpacing), needed = !bk.buried && bk.hobbsUtil > 1, nd = bk.nd;
  const base = { EI: bk.EI, EA: bk.EA, w: bk.w, muA: p.muAxial, muL: p.muLateral, P0: bk.comp, ro: p.D / 2 }, plan = lateralBuckle({ ...base, spacing: managed ? spacing : Infinity }), rogue = managed ? lateralBuckle({ ...base, spacing: Infinity }) : plan;
  const fo = sa.fmtOf(nd, p.t, nd.T), cap = strainCapacity({ D: p.D, t: Math.max(p.t - p.CA / 1000, 0.3 * p.t), fy: fo.fy, fu: fo.fu, dp: 0 }), allow = p.strainAllow > 0 ? p.strainAllow / 100 : cap / (DNV_STRAIN[p.sc] || 2.5), g = plan.governing, restr = st.nodes.filter((n) => !n.riser), Lrestr = restr.length ? restr[restr.length - 1].x - restr[0].x : 0;
  return { needed, managed, spacing, count: needed && managed ? Math.max(1, Math.round(Lrestr / spacing)) : 0, plan, rogue, capacity: cap, allow, strain: g ? g.strain : 0, util: g ? g.strain / allow : 0, rogueStrain: rogue.governing ? rogue.governing.strain : 0, source: net && net.spacing > 0 ? 'network suite' : 'inputs' };
}
/** Loads on the listed equipment items: slug or momentum force at the item, moment = force × lever arm, against the allowable nozzle or hub loads. */
function equipmentLoads(p, st, ds, sl) {
  return p.equipment.map((e) => {
    const n = st.at(e.x), near = sl.bends.reduce((a, b) => (Math.abs(b.x - n.x) < Math.abs(a.x - n.x) ? b : a)), A = e.size > 0 ? (Math.PI / 4) * (e.size / 1000) ** 2 : p.Ai, F = Math.abs(near.x - n.x) <= 50 ? near.force : (ds.on ? p.slugCf * p.dlf : 1) * n.rhoNS * p.Ai * n.vm ** 2 * Math.SQRT2;
    const M = F * e.arm, uF = e.fAllow > 0 ? F / (e.fAllow * 1000) : 0, uM = e.mAllow > 0 ? M / (e.mAllow * 1000) : 0, thrust = Math.max(n.pi - n.pe, 0) * A;
    return { ...e, xNode: n.x, F, M, uF, uM, util: Math.max(uF, uM), thrust, velocity: (n.vm * p.Ai) / A };
  });
}
/** Continuum and shell models of the pipe wall: Lamé check, thinned band (stress concentration, net-section collapse), shell cross-check, ring with out-of-roundness or a dent. */
function continuum(p, st, sa, wall, def, lossEol) {
  const row = sa.rows.reduce((a, r) => (sa.pLocal(r.nd) - r.nd.pe > sa.pLocal(a.nd) - a.nd.pe ? r : a)), nd = row.nd, pD = sa.pLocal(nd), dp = Math.max(pD - nd.pe, 1), ro = p.D / 2, ri = ro - p.t, order = p.feOrder, E = p.E, nu = p.nu;
  const cyl = feCylinder({ ri, ro, pi: pD, pe: nd.pe, E, nu, nr: p.feNr, nz: 1, order }), la = lame(pD, nd.pe, ri, ro, ri);
  // thinned band: user depth, else the deepest measured defect, else the wall loss predicted at the end of the design life
  const w = def.worst, dMm = p.feDepthMm > 0 ? p.feDepthMm : w ? w.depth : Math.max(lossEol, 0.02 * p.wtMm), d = clamp(dMm / 1000, 0.01 * p.t, 0.9 * p.t), Lb = (w && !(p.feDepthMm > 0) ? w.length : p.feLenMm) / 1000, flow = p.S + 69e6;
  const gr = feGroove({ ro, t: p.t, d, L: Lb, side: 'in', pi: dp, pe: 0, E, nu, nz: p.feNz, nr: p.feNr, order, flow }), codes = {}; for (const m of ['b31g', 'modified', 'dnv']) codes[m] = b31g({ D: p.D, t: p.t, d, L: Lb, smys: p.S, smts: p.Su, method: m }).pf;
  // the same band as a thin shell with a thickness step (membrane + edge bending)
  const Rm = ro - p.t / 2, Ls = Lb / 2 + 5 * Math.sqrt(Rm * p.t), sh = shellCylinder({ R: Rm, L: Ls, t: (x) => (x < Lb / 2 ? p.t - d : p.t), E, nu, p: dp, n: p.shellN, axial: 'restrained', ends: ['symmetry', 'symmetry'], breaks: [Lb / 2] }), shellFactor = sh.hoopMembrane[0] / ((dp * Rm) / p.t);
  // cross-section as a ring: out-of-roundness (and a dent when given) under the largest external overpressure
  const deep = st.nodes.reduce((a, b) => (b.pe > a.pe ? b : a)), pex = Math.max(deep.pe - sa.pMin, 0), Rr = wall.tNow, ring = feRing({ ro, t: Rr, pi: 0, pe: Math.max(pex, 1), E, nu, ovality: p.f0, dent: p.dentMm > 0 ? { depth: p.dentMm / 1000, halfAngle: clamp((3 * Math.sqrt((ro * p.dentMm) / 1000)) / ro, 0.15, 1.2) } : null, nth: Math.max(p.feNz, 16), nr: p.feNr, order });
  const w1 = (p.f0 * (2 * ro - Rr)) / 4, ringHand = ((Math.max(pex, 1) * ro) / Rr) * (1 + (6 * w1) / Rr), hoopPeak = Math.max(Math.abs(ring.hoopMax), Math.abs(ring.hoopMin));
  return { nd, pD, dp, cyl, la, gr, codes, d, Lb, flow, sh, shellFactor, ring, ringHand, hoopPeak, pex, deep, ringUtil: pex > 0 ? hoopPeak / sa.fyCold : 0, source: p.feDepthMm > 0 ? 'input' : w ? w.name : 'wall loss at the end of the design life' };
}
/** Tracked-particle erosion of the worst bend, beside the DNV-RP-O501 bend equation. */
function trackedErosion(p, deg) {
  const b = deg.bWorst, nd = b.nd, model = p.erosionModel === 'oka' || p.erosionModel === 'finnie' ? p.erosionModel : 'dnv', hv = (p.smts / 3.2) * 9.80665e-3, U = Math.max(nd.vm - p.eroVth, 0);
  const tr = bendErosionTracked({ D: p.ID, rOverD: b.rOverD, angle: b.angle, U, rho: nd.rhoNS, mu: nd.muNS, dp: deg.dp, rhoP: p.sandDensity, mp: deg.mp, model, hv, flowStress: hv * 1e9, finnieC: p.finnieC, nPart: p.nPart, ns: p.bendNs, nn: Math.max(6, Math.round(p.bendNs / 4)), en: p.restitution, et: Math.min(p.restitution + 0.1, 1) });
  const ref = erosionRate({ model: 'dnv', mp: deg.mp, U, D: p.ID, dp: deg.dp, rhoM: nd.rhoNS, muM: nd.muNS, rhoP: p.sandDensity, rOverD: b.rOverD, gf: 1, mult: 1 });
  return { tr, bend: b, model, rateMax: tr.rateMax * p.erosionMult, rateMean: tr.rateMean * p.erosionMult, dnv: ref.rate * p.erosionMult, dnvNoC1: (ref.rate * p.erosionMult) / O501.C1, G: ref.G, U };
}

// Model uncertainty of the burst capacity (test / prediction), fitted by maximum likelihood to every second test of the burst database (see burstModelFit and the reference data sets):
// defect: DNV-RP-F101 failure pressure with the tensile strength of the test, log-normal (μ_ln = 0.16245, σ_ln = 0.17211); intact: pressure-containment format with the measured strengths, normal.
export const BURST_MODEL = Object.freeze({ defect: Object.freeze({ mu: 0.16245, sigma: 0.17211, mean: 1.1939, cov: 0.17339 }), intact: Object.freeze({ mean: 1.0697, sd: 0.07126 }), old: Object.freeze({ mean: 1.05, sd: 0.105 }) });
/**
 * Burst-model uncertainty from burst tests: ratio measured / predicted of every test, split into a fitting half (rows 1, 3, 5, … of each group in the order of the database) and a checking half.
 * rows: [{ D, t (mm), d, L (mm), yield, uts (MPa), Pb (MPa) }]. Returns { defect: { mu, sigma (of ln ratio), mean, cov, n, check: { n, ks, ksCritical (5 %), coverage90, below5, ratios[] } }, intact: { mean, sd, n, check: { n, mean, sd } }, all: { mean, sd, n } }.
 */
export function burstModelFit(rows) {
  const ratio = (r) => r.Pb / (b31g({ D: r.D / 1000, t: r.t / 1000, d: r.d / 1000, L: r.L / 1000, smys: r.yield * MPA, smts: r.uts * MPA, method: 'dnv' }).pf / MPA), ms = (x) => { const m = mean(x); return [m, Math.sqrt(mean(x.map((v) => (v - m) ** 2)))]; };
  const dr = rows.filter((r) => r.d > 0).map(ratio), cal = dr.filter((_, i) => i % 2 === 0), val = dr.filter((_, i) => i % 2 === 1), [mu, sigma] = ms(cal.map(Math.log)), q = (pr) => Math.exp(mu + sigma * PhiInv(pr)), sv = val.slice().sort((a, b) => a - b); let ks = 0;
  sv.forEach((x, i) => { const F = Phi((Math.log(x) - mu) / sigma); ks = Math.max(ks, Math.abs(F - i / sv.length), Math.abs(F - (i + 1) / sv.length)); });
  const ir = rows.filter((r) => r.d === 0).map((r) => r.Pb / (((2 * r.t) / (r.D - r.t)) * Math.min(r.yield, r.uts / 1.15) * 1.1547005)), ic = ir.filter((_, i) => i % 2 === 0), iv = ir.filter((_, i) => i % 2 === 1), [im, isd] = ms(ic), [vm2, vs2] = iv.length ? ms(iv) : [0, 0], [am, as] = ms(rows.map(ratio));
  return { defect: { mu, sigma, mean: Math.exp(mu + (sigma * sigma) / 2), cov: Math.sqrt(Math.exp(sigma * sigma) - 1), n: cal.length, check: { n: val.length, ks, ksCritical: 1.36 / Math.sqrt(Math.max(val.length, 1)), coverage90: val.filter((x) => x >= q(0.05) && x <= q(0.95)).length / Math.max(val.length, 1), below5: val.filter((x) => x < q(0.05)).length, ratios: sv } }, intact: { mean: im, sd: isd, n: ic.length, check: { n: iv.length, mean: vm2, sd: vs2 } }, all: { mean: am, sd: as, n: rows.length } };
}
/** Limit states, FORM / SORM / sampling comparison and the reliability index over time. */
async function reliability(p, st, wall, sa, deg, span, defectWorst, ctx) {
  const tN = p.wtMm, Dmm = p.D * 1000, ratioU = p.smts / p.smys, cW = p.covWt / 100, cY = p.covYield / 100, cP = p.covPress / 100, cC = p.covCorr / 100, age = p.age, life = Math.max(p.designLife, age + 1);
  const shut = sa.rows.map((r) => (Math.max(p.pShutIn * BAR + sa.head(r.nd), r.nd.pi) - r.nd.pe) / MPA), iB = shut.reduce((k, x, i) => (x > shut[k] ? i : k), 0), pMean = Math.max(shut[iB], 0.01), peMax = Math.max(...st.nodes.map((n) => n.pe - PATM)) / MPA;
  // distributions fitted by maximum likelihood to measured samples replace the assumed means and scatters
  const fits = [], fitOf = (key, label, dist) => { const xs = p.rvSamples.filter((r) => r.variable.startsWith(key)).map((r) => r.value); if (xs.length < 5) return null; try { const ft = fitDistribution(xs, dist); fits.push({ label, ...ft }); return ft; } catch { return null; } };
  const fW = fitOf('wall', 'Wall thickness (mm)', 'normal'), fY = fitOf('yield', 'Yield strength (MPa)', 'lognormal'), fP = fitOf('pres', 'Annual extreme pressure (bar)', 'gumbel'), fR = fitOf('rate', 'Wall-loss rate (mm/y)', 'lognormal'), cCx = fR ? fR.cov : cC, mC = fR && wall.rate > 1e-9 ? fR.mean / wall.rate : 1;
  const vT = fW ? { name: 'Wall thickness', dist: 'normal', mean: fW.mean, sd: fW.sd } : { name: 'Wall thickness', dist: 'normal', mean: tN, cov: cW }, vY = fY ? { name: 'Yield strength', dist: 'lognormal', mean: fY.mean, sd: fY.sd } : { name: 'Yield strength', dist: 'lognormal', mean: p.yieldBias * p.smys, cov: cY }, vP = fP ? { name: 'Annual extreme pressure', dist: 'gumbel', mean: fP.mean / 10, sd: fP.sd / 10 } : { name: 'Annual extreme pressure', dist: 'gumbel', mean: pMean, cov: cP }, vC = { name: 'Wall-loss model', dist: 'lognormal', mean: mC, cov: cCx };
  const peEff = Math.max(Math.max(...st.nodes.map((n) => n.pe - sa.pMin)) / MPA, 1e-3), eroRate = deg.bWorst.ero.gov, ero0 = Math.min(eroRate * age, 0.9 * tN);
  const lossOf = (rate, loss0) => (T, xc) => loss0 + xc * rate * Math.max(T - age, 0), lossGen = lossOf(wall.rate, wall.lossNow), eroLoss = lossOf(deg.bWorst.ero.gov, Math.min(deg.bWorst.ero.gov * age, 0.9 * tN)), fyCold = sa.fmtOf(st.nodes[0], p.t, p.tInstall).fy / p.S;
  const LS = [
    { key: 'burst', name: 'Burst of the corroding wall at the annual extreme pressure', annual: true, wl: { rate: wall.rate, loss0: wall.lossNow, ix: 3 }, vars: [vT, vY, vP, vC, { name: 'Burst model', dist: 'normal', mean: BURST_MODEL.intact.mean, sd: BURST_MODEL.intact.sd }], g: (T) => (x) => { const tw = Math.max(x[0] - lossGen(T, x[3]), 0.01); return (x[4] * ((2 * tw) / (Dmm - tw)) * Math.min(x[1], (x[1] * ratioU) / 1.15) * 1.1547005 - x[2]) / pMean; } },
    { key: 'leak', name: 'Corrosion: wall loss reaches the leak criterion', annual: false, wl: { rate: deg.corrMax, loss0: wall.lossNow, ix: 1 }, vars: [vT, vC], g: (T) => (x) => (p.leakFrac * x[0] - lossOf(deg.corrMax, wall.lossNow)(T, x[1])) / tN },
    { key: 'erosion', name: 'Erosion at the worst bend reaches the leak criterion', annual: false, wl: { rate: eroRate, loss0: ero0, ix: 1 }, vars: [vT, { ...vC, name: 'Erosion model', cov: Math.max(cC, 1) }], g: (T) => (x) => (p.leakFrac * x[0] - eroLoss(T, x[1])) / tN },
    { key: 'fatigue', name: 'Fatigue: Miner sum reaches its resistance', annual: false, vars: [{ name: 'Miner resistance', dist: 'lognormal', mean: 1, cov: 0.3 }, { name: 'Stress / S–N model', dist: 'lognormal', mean: 1, cov: 0.6 }], g: (T) => (x) => x[0] - (p.priorDamage + span.dYear * Math.max(T - age, 0)) * x[1] },
  ];
  if (peMax > 0.01) LS.push({ key: 'collapse', name: 'External collapse of the depressurised pipe', annual: true, wl: { rate: wall.rate, loss0: wall.lossNow, ix: 2 }, vars: [vT, vY, vC, { name: 'Collapse model', dist: 'normal', mean: 1, cov: 0.08 }, { name: 'Out-of-roundness', dist: 'lognormal', mean: p.f0, cov: 0.25 }], g: (T) => (x) => { const tw = Math.max(x[0] - lossGen(T, x[2]), 0.05); return (x[3] * collapsePressure({ D: Dmm, t: tw, E: p.eMod, nu: p.nu, fy: x[1] * fyCold, ovality: x[4], alphaFab: p.alphaFab }).pc) / peEff - 1; } });
  if (defectWorst) LS.push({ key: 'defect', name: `Burst at the worst measured defect (x = ${fmt(defectWorst.x)} m)`, annual: true, vars: [{ name: 'Defect depth sizing', dist: 'normal', mean: 1, cov: 0.1 }, vC, { ...vP, mean: Math.max(defectWorst.dpExt / MPA, 0.01) }, { name: 'Tensile strength', dist: 'lognormal', mean: 1.09 * p.smts, cov: cY }, { name: 'Burst model', dist: 'lognormal', mean: BURST_MODEL.defect.mean, cov: BURST_MODEL.defect.cov }],
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
  // ---- inspection and repair plan: conditional failure probability against wall loss, then the plan with Bayesian updating
  const years = Math.max(1, Math.round(life - age)), podLam = (-0.1 * tN) / Math.log(1 - clamp(p.pod, 1e-6, 0.999)), podOf = (l) => 1 - Math.exp(-Math.max(l, 0) / podLam), sizing = p.sizingMm > 0 ? p.sizingMm : (0.1 * tN) / 1.645, thr = Math.max(p.repairFrac * (p.CA > 0 ? p.CA : 0.2 * tN), 0.02 * tN), capL = 0.95 * tN, lg = linspace(0, capL, 7);
  const wlLS = LS.filter((l) => l.wl && l.wl.rate > 1e-9);
  for (const ls of wlLS) {
    const bt = lg.map((l) => { const vars = ls.vars.map((v, i) => (i === ls.wl.ix ? { name: v.name, dist: 'det', mean: (l - ls.wl.loss0) / ls.wl.rate } : v)); let r; try { r = form(ls.g(age + 1), vars); } catch { r = { beta: 10 }; } return clampB(r.beta); });
    ls.pfCond = (l) => Phi(-interp1(lg, bt, clamp(l, 0, capL))); ls.betaCond = bt; ls.wl.cov = Math.sqrt(Math.log(1 + (ls.vars[ls.wl.ix].cov ?? cCx) ** 2)) > 0 ? ls.vars[ls.wl.ix].cov ?? cCx : cCx; ls.wl.mean = ls.vars[ls.wl.ix].mean;
  }
  const planOf = (iv) => inspectionPlan({ ls: wlLS.map((ls) => ({ name: ls.name, pf: ls.pfCond, annual: ls.annual, rate: ls.wl.rate * ls.wl.mean, loss0: ls.wl.loss0, cov: ls.wl.cov })), years, interval: iv, pod: podOf, sizingSd: sizing, threshold: thr, nq: 40, cap: capL });
  const others = LS.filter((l) => !wlLS.includes(l)), otherAt = (y) => 1 - others.reduce((q, ls) => q * (1 - clamp(interp1(Tg, ls.annualT, Math.min(age + y, life)), 0, 1)), 1), totalOf = (pl) => pl.t.map((y, i) => 1 - (1 - (wlLS.length ? pl.annual[i] : 0)) * (1 - otherAt(y)));
  const cand = p.inspInterval > 0 ? [p.inspInterval] : [10, 8, 6, 5, 4, 3, 2, 1].filter((iv) => iv < years || iv === 1); let chosen = null;
  for (const iv of cand) { const pl = planOf(iv), tot = totalOf(pl), mx = Math.max(...tot); chosen = { interval: iv, plan: pl, annual: tot, max: mx }; if (mx <= target) break; }
  const none = totalOf(planOf(Infinity)), iX = chosen.annual.findIndex((q) => q > target);
  wlLS.forEach((ls, k) => { ls.managedAvg = mean(chosen.plan.perLS[k]); ls.managedMax = Math.max(...chosen.plan.perLS[k]); });
  const managed = { interval: chosen.interval, auto: !(p.inspInterval > 0), t: chosen.plan.t.map((y) => age + y), annual: chosen.annual, max: chosen.max, avg: mean(chosen.annual), met: chosen.max <= target, tCross: iX >= 0 ? age + iX + 1 : null, none, noneMax: Math.max(...none), noneAvg: mean(none), inspections: chosen.plan.inspections, repairs: chosen.plan.repairs, threshold: thr, sizing, podLam, podAt: podOf, years };
  return { LS, Tg, pofNow, pofEol, pofAvg, pofT, target, tTarget, gov, f, so, mc, is, lh, rs, life, iB, pMean, peMax, peEff, managed, fits };
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
    const safe = F * r.modified, erf = safe > 0 ? dp / safe : 99, rate = Math.max(deg.corr[nd.i].rate, wall.rate) * p.defectGrowth, erfAt = (y) => { const dd = Math.min(dm + (rate * y) / 1000, 0.999 * p.t), s = F * b31g({ ...base, d: dd, L: Lm, method: 'modified' }).pf; return s > 0 ? dp / s : 99; };
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
  const remYears = Math.max(p.designLife - p.age, 0), fatUse = (s) => (p.priorDamage + s.dYear * p.dff * remYears) * p.fatMargin;
  let span = spanDynamics(p, st, ds, sl, sa, ctx), rectified = false; const spanEntered = span;
  if (p.spanMgmt === 'rectify' && fatUse(span) > 1 && span.dYear > 0) { // support spacing: the longest span that meets the criterion
    let lo = Math.max(0.2 * p.spanLength, 2), hi = p.spanLength, sLo = null;
    for (let k = 0; k < 6; k++) { const L = k === 0 ? p.spanLength * clamp(fatUse(span) ** (-1 / 8), 0.35, 0.97) : 0.5 * (lo + hi), s = spanDynamics({ ...p, spanLength: L }, st, ds, sl, sa, ctx); if (fatUse(s) <= 1) { lo = L; sLo = s; if (hi - lo < 0.04 * hi) break; } else hi = L; await tick(); }
    if (sLo) { span = sLo; p.spanLength = lo; rectified = true; }
  }
  // the same spectrum of the span as entered with the weld made to the widest tolerance: the fatigue result without the design measures
  const dUnm = sum(spanEntered.spectrum.map((s) => s.perYear / snCycles(s.range, p.unmanagedWeld.snClass, p.snEnv, { scf: p.unmanagedWeld.scf, t: p.wtMm, shift: p.snShift }))), fatCalc = cap(Math.max(1 - p.priorDamage, 0) / Math.max(span.dYear, 1e-12), 1e4), fatCalcUnm = cap(Math.max(1 - p.priorDamage, 0) / Math.max(dUnm, 1e-12), 1e4);
  const spanAllow = span.dYear > 0 ? p.spanLength * clamp(fatUse(span) ** (-1 / 8), 0.2, 5) : null;
  const fatLife = cap(Math.max(1 - p.priorDamage, 0) / Math.max(span.dYear * p.dff, 1e-12), 1e4), lifeMin = cap(Math.max((wall.tNow - Math.max(...sa.rows.map((r) => r.tReq))) * 1000, 0) / Math.max(rate, 1e-9), 200);
  // local buckling at the span under pressure, effective force and static + dynamic bending (functional load factor 1.2)
  const mSpan = ((span.sigStatic + span.dynMax) * MPA * p.Isteel) / (p.D / 2), lbOp = localBuckling({ D: p.D, t: wall.tNow, fy: sa.fmtOf(span.nd, wall.tNow, span.nd.T).fy, fu: sa.fmtOf(span.nd, wall.tNow, span.nd.T).fu, E: p.E, nu: p.nu, pi: span.nd.pi, pe: span.nd.pe, pmin: PATM, M: 1.2 * mSpan, S: 1.2 * span.seff, safetyClass: p.sc, ovality: p.f0, alphaFab: p.alphaFab });
  // ---- fracture mechanics on the fatigue spectrum
  const cycYr = sum(span.spectrum.map((s) => s.perYear)), dSeq = cycYr > 0 ? (sum(span.spectrum.map((s) => s.perYear * s.range ** p.parisM)) / cycYr) ** (1 / p.parisM) * p.scf : 0, tw = wall.tNow, a0 = Math.min(p.flawMm / 1000, 0.5 * tw);
  const sMaxT = Math.max(...sa.rows.map((r) => r.sOp.longitudinal), 0) / MPA + (span.sigStatic + span.dynMax) * p.scf + Math.max(p.residualStress, 0), kOf = (a) => edgeCrackY(a, tw) * sMaxT * Math.sqrt(Math.PI * a) - p.kMat, aC = sMaxT > 0 && kOf(0.8 * tw) > 0 ? (kOf(a0) >= 0 ? a0 : brent(kOf, a0, 0.8 * tw, 1e-9)) : 0.8 * tw;
  const crack = parisLife({ a0, ac: Math.max(aC, a0), C: p.parisC, m: p.parisM, dS: dSeq, Y: (a) => edgeCrackY(a, tw), dKth: p.dKth }), crackLife = cycYr > 0 && Number.isFinite(crack.N) ? cap(crack.N / cycYr, 1e4) : 1e4;
  // ---- defects, blockage
  const def = assessDefects(p, st, sa, deg, wall, mapS), beta = clamp(Math.max(p.blockage, p.effIdMm > 0 ? 1 - (p.effIdMm / p.idMm) ** 2 : 0, p.deposit0 > 0 ? 1 - Math.max(1 - (2 * p.deposit0) / p.idMm, 0) ** 2 : 0), 0, 1), ndB = st.at(p.plugX), plugged = beta >= 0.999;
  const dpRestr = plugged ? 0 : 0.5 * ndB.rhoM * ndB.vm ** 2 * (1 / (1 - beta) - 1) ** 2 + ((0.02 * ndB.rhoM * ndB.vm ** 2) / (2 * p.ID)) * p.plugLen * ((1 - beta) ** -2.5 - 1), pUp = plugged ? p.pShutIn : N[0].P + dpRestr / BAR;
  const mawpInlet = sa.mawp + (p.pRefLoc === 'top' ? (p.rhoContent * G * (st.zTop - N[0].z)) / BAR : 0), plugDp = Math.max(p.pShutIn - Math.min(ndB.P, N[N.length - 1].P), 0) * BAR, plugMass = (917 * (1 - p.plugPorosity) + 1000 * p.plugPorosity) * p.Ai * p.plugLen;
  const runLen = Math.max(Math.min(...p.bends.map((b) => b.x - ndB.x).filter((d) => d > 1), N[N.length - 1].x - ndB.x, 1000), 10), plugV = Math.min(Math.sqrt((2 * plugDp * p.Ai * runLen) / Math.max(plugMass, 1)), 300), plugE = 0.5 * plugMass * plugV ** 2;
  // ---- electrochemistry, sour service, top-of-line condensation
  const ndC = N[deg.iC], cw = deg.corr[deg.iC].dw, cg = deg.corr[deg.iC], mix = mixedPotential({ T: ndC.T, pCO2: Math.max(cw.fCO2, 1e-6), pH: p.phAct > 0 ? p.phAct : cg.pH, U: Math.max(ndC.vL, 0.01), d: p.ID, salt: p.brineWt });
  const pH2S = Math.max(...N.map((n) => (clamp(p.h2s, 0, 100) / 100) * n.P * 100)), sour = sourRegion(pH2S, p.phAct > 0 ? p.phAct : cw.pHco2), qGs = (st.fm.rates.qGasStd || 0) / 86400;
  const tlc = N.map((n, i) => { const j = Math.min(i + 1, N.length - 1), k = Math.max(i - 1, 0), dx = Math.max(N[j].x - N[k].x, 1e-6); return condensationRate({ P: n.P, T: n.T, dTdx: (N[j].T - N[k].T) / dx, qGasStd: qGs, D: p.ID }); }), tlcMax = Math.max(...tlc);
  // ---- reliability
  prog(0.55, 'Reliability: FORM, SORM and sampling');
  const rel = await reliability(p, st, wall, sa, deg, span, def.worst, ctx), LSof = (k) => rel.LS.find((l) => l.key === k), mg = rel.managed, inspInt = mg.interval, pRaw = (k) => LSof(k)?.annualAvg ?? 0, pA = (k) => { const l = LSof(k); return l ? clamp(l.managedAvg ?? l.annualAvg, 0, 1) : 0; }, barrier = (k) => (pRaw(k) > 0 ? clamp(pA(k) / pRaw(k), 0, 1) : 1);
  const arr = arrestorDesign(p, st, sa, ctx), lat = lateralDesign(p, st, sa, ctx), eqp = equipmentLoads(p, st, ds, sl);
  prog(0.8, 'Risk models'); await tick();
  // ---- risk: consequence, fault tree, event tree, matrix, Markov, Bayesian network
  const prodDay = (st.fm.rates.qOilStd || 0) * 6.2898 * p.oilPrice, cBase = p.repairCost * 1e6 + p.downtimeDays * prodDay, betaT = -PhiInv(rel.target), utilPof = (u) => (u > 0 ? Phi(-(betaT - Math.log(u) / 0.15)) : 0);
  const pLB = utilPof(lbOp.util), pUph = sa.buckling.buried ? utilPof(sa.buckling.palmerUtil) : 0, pCol = pA('collapse'), pPlugLoc = clamp(p.plugProb, 0, 1) * clamp(p.pfdProtect, 0, 1);
  const tree = { name: 'Loss of containment', gate: 'OR', children: [
    { name: 'Undetected corrosion failure', gate: 'AND', children: [{ name: 'Corrosion wall loss reaches the limit', p: pRaw('leak') }, { name: 'Inspection and repair do not intervene (corrosion)', p: barrier('leak') }] },
    { name: 'Undetected erosion failure', gate: 'AND', children: [{ name: 'Erosion wall loss reaches the limit', p: pRaw('erosion') }, { name: 'Inspection and repair do not intervene (erosion)', p: barrier('erosion') }] },
    { name: 'Fatigue crack through the wall', p: pA('fatigue') }, { name: 'Burst at the annual extreme pressure', p: Math.max(pA('burst'), pA('defect')) }, { name: 'Collapse or local buckling', p: clamp(pCol + pLB + pUph, 0, 1) },
    { name: 'Plug-related failure', gate: 'AND', children: [{ name: 'Hydrate or deposit plug forms', p: clamp(p.plugProb, 0, 1) }, { name: 'Plug remediation barrier fails', p: clamp(p.pfdProtect, 0, 1) }] },
    { name: 'External damage', p: clamp(p.extFreq, 0, 1) }] };
  const ft = faultTree(tree), et = eventTree(ft.top, [{ name: 'Leak detected', p: p.pDetect }, { name: 'Isolation succeeds', p: p.pIsolate }, { name: 'Ignition', p: p.pIgnite }]);
  const etRows = et.outcomes.map((o) => { const [d, iso, ign] = o.path, rel0 = d && iso ? 0.05 : 1, cost = cBase + rel0 * p.envCost * 1e6 + (ign ? p.safetyCost * 1e6 : 0); return { ...o, cost, outcome: `${d ? (iso ? 'Small release, isolated' : 'Release continues until depressurised') : 'Late discovery, full release'}${ign ? ', ignited' : ''}` }; });
  const cLoc = sum(etRows.map((o) => o.p * o.cost)), cPlug = p.plugDays * prodDay + 0.3 * p.repairCost * 1e6, cBuckle = cBase + p.replaceCost * 1000 * (arr.required ? (arr.provided ? arr.confined : arr.unconfined) : 0), cBuckleOpen = cBase + p.replaceCost * 1000 * (arr.required ? arr.unconfined : 0);
  const threats = [
    { name: 'CO2 corrosion', pof: pA('leak'), cof: cLoc, det: 4, cause: `${fmt(deg.corrMax, 3)} mm/y at ${fmt(ndC.x / 1000, 3)} km`, effect: 'Pinhole leak, then rupture', action: 'Inhibition availability, inline inspection' },
    { name: 'Sand erosion', pof: pA('erosion'), cof: cLoc, det: 5, cause: `${fmt(deg.bWorst.ero.gov, 3)} mm/y at ${deg.bWorst.name}`, effect: 'Wall thinning at fittings', action: 'Sand monitoring, velocity limit, UT at bends' },
    { name: 'Slugging and VIV fatigue', pof: pA('fatigue'), cof: cLoc, det: 7, cause: `damage ${fmt(span.dYear, 3)} per year on the ${fmt(p.spanLength, 3)} m span`, effect: 'Girth-weld crack', action: 'Support the span, suppress slugging' },
    { name: 'Burst at shut-in pressure', pof: Math.max(pA('burst'), pA('defect')), cof: cLoc, det: 6, cause: `annual extreme ${fmt(rel.pMean * 10, 3)} bar differential`, effect: 'Rupture', action: 'Overpressure protection, defect repair' },
    { name: 'External collapse', pof: pCol, cof: cBuckle, det: 8, cause: `${fmt(sa.collapseUtil, 3)} utilisation when depressurised`, effect: arr.provided ? `Flattened pipe, confined to ${fmt(arr.confined, 3)} m by the arrestors` : 'Flattened pipe, flooding', action: 'Keep the minimum internal pressure, inspect the wall, maintain the buckle arrestors' },
    { name: 'Local or global buckling', pof: clamp(pLB + pUph, 0, 1), cof: cBase, det: 6, cause: `local buckling utilisation ${fmt(lbOp.util, 3)}`, effect: 'Wrinkle or upheaval', action: 'Span correction, rock cover' },
    { name: 'Hydrate or deposit plug', pof: clamp(p.plugProb, 0, 1), cof: cPlug, det: 3, cause: `${fmt(100 * beta, 3)} % of the bore lost`, effect: `Production stopped for ${fmt(p.plugDays, 3)} d`, action: 'Inhibitor, insulation, no one-sided depressurisation' },
    { name: 'Plug projectile / overpressure', pof: pPlugLoc, cof: cLoc, det: 7, cause: `${fmt(plugDp / BAR, 3)} bar across the plug`, effect: 'Bend rupture', action: 'Two-sided depressurisation procedure' },
    { name: 'External damage', pof: clamp(p.extFreq, 0, 1), cof: cLoc, det: 8, cause: 'Anchors, dropped objects, trawling', effect: 'Dent or rupture', action: 'Protection, exclusion zone' },
  ].map((t) => { const pc = category(t.pof, POF_EDGES), cc = category(t.cof, COF_EDGES), lv = riskLevel(pc, cc); return { ...t, pc, cc, level: lv, risk: t.pof * t.cof, S: Math.min(2 * cc, 10), O: Math.min(2 * pc, 10), rpn: Math.min(2 * cc, 10) * Math.min(2 * pc, 10) * t.det }; });
  const riskCost = sum(threats.map((t) => t.risk)), topLevel = Math.max(...threats.map((t) => t.level)), topThreat = threats.reduce((a, b) => (b.risk > a.risk ? b : a));
  // the same threats without the design measures (no inspection credit, unconfined buckle): published as the unmanaged counterpart
  const openPof = { 'CO2 corrosion': pRaw('leak'), 'Sand erosion': pRaw('erosion'), 'Burst at shut-in pressure': Math.max(pRaw('burst'), pRaw('defect')), 'External collapse': pRaw('collapse') }, threatsOpen = threats.map((t) => { const pof = openPof[t.name] ?? t.pof, cof = t.name === 'External collapse' ? cBuckleOpen : t.cof; return { pof, cof, level: riskLevel(category(pof, POF_EDGES), category(cof, COF_EDGES)) }; }), riskCostOpen = sum(threatsOpen.map((t) => t.pof * t.cof)), topLevelOpen = Math.max(...threatsOpen.map((t) => t.level));
  // Markov degradation: as-new → degraded → critical → failed, repair of the critical state after a successful inspection
  const half = Math.max(p.CA / 2, 0.05 * p.wtMm), m3 = Math.max(0.8 * p.wtMm - 2 * half, 0.05 * p.wtMm), l12 = rate / half, l34 = rate / m3, mu = clamp(p.pod, 0, 1) / Math.max(inspInt, 0.05), s0 = lossNow < half ? 0 : lossNow < 2 * half ? 1 : 2, p0 = [0, 0, 0, 0]; p0[s0] = 1;
  const horizon = Math.max(remYears, 1), mk = markov([[0, l12, 0, 0], [0, 0, l12, 0], [mu, 0, 0, l34], [0, 0, 0, 0]], p0, horizon, 60), mk0 = markov([[0, l12, 0, 0], [0, 0, l12, 0], [0, 0, 0, l34], [0, 0, 0, 0]], p0, horizon, 60);
  // Bayesian network: wetting and inhibition → corrosion → wall loss → leak, with an inspection finding as evidence
  const wetFrac = clamp(deg.corr.filter((c) => !c.oilWet).length / N.length, 0.02, 0.98), allowRate = Math.max(p.CA, 0.1) / p.designLife, zc = Math.sqrt(Math.log(1 + cC * cC)) || 0.3, pHigh = (r) => clamp(1 - Phi(Math.log(allowRate / Math.max(r, 1e-9)) / zc + zc / 2), 0.001, 0.999), bare = deg.corr[deg.iC].bare;
  const bnNodes = [{ name: 'Water wetting', states: ['yes', 'no'], cpt: [[wetFrac, 1 - wetFrac]] }, { name: 'Inhibitor working', states: ['yes', 'no'], cpt: [[clamp(p.eta, 0.01, 0.999), 1 - clamp(p.eta, 0.01, 0.999)]] },
    { name: 'Corrosion', states: ['high', 'low'], parents: ['Water wetting', 'Inhibitor working'], cpt: [[pHigh(bare * (1 - p.eta)), 1 - pHigh(bare * (1 - p.eta))], [pHigh(bare), 1 - pHigh(bare)], [0.01, 0.99], [0.02, 0.98]] },
    { name: 'Wall loss', states: ['severe', 'minor'], parents: ['Corrosion'], cpt: [[0.85, 0.15], [0.03, 0.97]] }, { name: 'Leak', states: ['yes', 'no'], parents: ['Wall loss'], cpt: [[0.25, 0.75], [0.001, 0.999]] },
    { name: 'Inspection', states: ['thinning found', 'no thinning'], parents: ['Wall loss'], cpt: [[clamp(p.pod, 0.01, 0.99), 1 - clamp(p.pod, 0.01, 0.99)], [0.05, 0.95]] }];
  const bn = [['No evidence (prior)', {}], ['Inspection finds thinning', { Inspection: 'thinning found' }], ['Inspection finds no thinning', { Inspection: 'no thinning' }]].map(([label, e]) => ({ label, r: bayesNet(bnNodes, e) })), bnNow = wall.measured ? bn[lossNow > 0.5 * Math.max(p.CA, 0.1) ? 1 : 2] : bn[0];
  // risk-based inspection interval
  const inspect = clamp(Math.min(inspInt, 0.5 * fatLife, topLevel >= 3 ? 1 : 10, def.worst ? Math.max(0.5 * def.worst.life, 0.5) : 10), 0.5, 10);

  // ---- warnings and recommendations -----------------------------------------------------------------------------------
  const W = (level, msg) => warnings.push({ level, msg }), V = (name, value, limit, where) => viol.push([name, cell(value), cell(limit), where]);
  const rH = sa.rows[sa.iH], rV = sa.rows[sa.iV], rC = sa.rows[sa.iC], codeName = { dnv: 'DNV-ST-F101', b314: 'ASME B31.4', b318: 'ASME B31.8' }[p.code];
  if (st.source !== 'flow suite') W('info', 'Pressures, temperatures and velocities come from the kernel steady-state estimate; run the flow suite for the detailed profile and the transient slug data.');
  if (sa.hoopUtil > 1) { W('bad', `Pressure containment utilisation ${fmt(sa.hoopUtil, 3)} exceeds 1.0 at ${fmt(rH.nd.x / 1000, 3)} km (${codeName}); the wall needed is ${fmt(sa.minWall, 3)} mm against ${fmt(p.wtMm, 3)} mm.`); V('Pressure containment', sa.hoopUtil, 1, `${fmt(rH.nd.x / 1000, 3)} km`); recs.push(`Increase the wall thickness to at least ${fmt(Math.ceil(sa.minWall * 10) / 10, 3)} mm or lower the design pressure to ${fmt(sa.mawpOf(sa.tDes), 3)} bara: the pressure-containment utilisation is ${fmt(sa.hoopUtil, 3)}.`); }
  if (sa.vmUtil > 1) { W('bad', `Combined (von Mises) stress utilisation ${fmt(sa.vmUtil, 3)} exceeds 1.0 at ${fmt(rV.nd.x / 1000, 3)} km.`); V('Combined stress', sa.vmUtil, 1, `${fmt(rV.nd.x / 1000, 3)} km`); recs.push(`Reduce the restrained thermal stress (expansion loop or lower design temperature than ${fmt(p.designTemp, 3)} °C): the equivalent stress reaches ${fmt(Math.max(rV.sDes.vonMises, rV.sOp.vonMises) / MPA, 3)} MPa.`); }
  if (sa.collapseUtil > 1) { W('bad', `External collapse utilisation ${fmt(sa.collapseUtil, 3)} exceeds 1.0 at ${fmt(rC.nd.depth, 4)} m water depth when the line is depressurised.`); V('External collapse', sa.collapseUtil, 1, `${fmt(rC.nd.depth, 4)} m depth`); recs.push(`Keep at least ${fmt(Math.max((rC.nd.pe - sa.col.pc / sa.gLB) / BAR, 1), 3)} bara inside the line at the deepest point or increase the wall: collapse capacity is ${fmt(sa.col.pc / BAR, 3)} bar.`); }
  else if (arr.required && !arr.provided) { W('warn', `A local buckle would propagate: external pressure is ${fmt(sa.propUtil, 3)} times the factored propagation pressure (${fmt(sa.col.ppr / BAR, 3)} bar) and no buckle arrestors are specified; ${fmt(arr.unconfined, 3)} m of line would flatten.`); recs.push(`Fit buckle arrestors below about ${fmt(arr.depthLimit, 3)} m water depth (about ${Math.ceil(arr.length / p.arrSpacing) + 1} at ${fmt(p.arrSpacing, 3)} m spacing): the propagation pressure is ${fmt(sa.col.ppr / BAR, 3)} bar against ${fmt((rC.nd.pe - PATM) / BAR, 3)} bar outside.`); }
  else if (arr.provided && arr.util > 1) { W('warn', `The buckle arrestors are too light: their crossover pressure of ${fmt(arr.pX / BAR, 3)} bar is below the factored external pressure (utilisation ${fmt(arr.util, 3)}), so a propagating buckle would pass them.`); V('Buckle-arrestor crossover', arr.util, 1, `${fmt(rC.nd.depth, 4)} m depth`); recs.push(`Thicken the arrestors beyond ${fmt(p.arrThick, 2)} × the pipe wall or lengthen them beyond ${fmt(p.arrLength, 2)} m: the crossover pressure must exceed ${fmt((arr.pe * 1.1 * sa.gLB) / BAR, 3)} bar.`); }
  if (p.designPressure > sa.mawp * 1.0001) { W('bad', `The design pressure ${fmt(p.designPressure, 4)} bara exceeds the allowable pressure of the present wall, ${fmt(sa.mawp, 4)} bara.`); V('Maximum allowable pressure', p.designPressure, sa.mawp, 'line'); }
  else if (p.designPressure > sa.mawpEol) { W('warn', `With the predicted wall at the end of the design life the allowable pressure drops to ${fmt(sa.mawpEol, 4)} bara, below the design pressure of ${fmt(p.designPressure, 4)} bara.`); recs.push(`Plan a re-rating to ${fmt(sa.mawpEol, 3)} bara or a repair before year ${fmt(p.age + remLife, 3)}: wall loss of ${fmt(rate, 3)} mm/y takes the allowable pressure below the design pressure.`); }
  if (p.blowdownMinT < p.mdmt) { W('bad', `Blowdown cools the metal to ${fmt(p.blowdownMinT, 3)} °C, below the minimum design metal temperature of ${fmt(p.mdmt, 3)} °C: brittle-fracture risk.`); V('Minimum metal temperature', p.blowdownMinT, p.mdmt, 'blowdown'); recs.push(`Slow the blowdown or specify impact-tested steel for ${fmt(p.blowdownMinT - 5, 3)} °C: the predicted minimum is ${fmt(p.mdmt - p.blowdownMinT, 3)} °C below the material limit.`); }
  if (remLife < remYears) { W(remLife < 0.5 * remYears ? 'bad' : 'warn', `The corrosion and erosion allowance (${fmt(p.CA, 3)} mm, ${fmt(margin, 3)} mm left) is consumed in ${fmt(remLife, 3)} y at ${fmt(rate, 3)} mm/y; the remaining design life is ${fmt(remYears, 3)} y.`); V('Remaining life on the allowance (y)', remLife, remYears, `${fmt(deg.gov.x / 1000, 3)} km`);
    const etaNeed = deg.gov.corr > 0 && remYears > 0 ? 1 - ((1 - p.eta) * Math.max(margin / remYears - deg.gov.ero, 0)) / deg.gov.corr : null;
    recs.push(`Remaining life ${fmt(remLife, 3)} y at ${fmt(rate, 3)} mm/y: inspect within ${fmt(inspect, 2)} y${etaNeed !== null && etaNeed < 0.999 && etaNeed > p.eta ? ` and raise the inhibitor efficiency × availability from ${fmt(100 * p.eta, 3)} % to ${fmt(100 * etaNeed, 3)} % to reach the design life` : etaNeed !== null && etaNeed >= 0.999 ? `; inhibition alone cannot reach the design life — a corrosion-resistant liner or ${fmt(rate * remYears, 2)} mm of allowance is needed` : ''}.`); }
  else recs.push(`Remaining life on the allowance is ${fmt(remLife, 3)} y at ${fmt(rate, 3)} mm/y against ${fmt(remYears, 3)} y still required; the risk-based inspection interval is ${fmt(inspect, 2)} y.`);
  if (deg.erosionalRatio > 1) { const nv = N[deg.iV]; W('warn', `Mixture velocity ${fmt(nv.vm, 3)} m/s is ${fmt(deg.erosionalRatio, 3)} times the API RP 14E erosional velocity (C = ${fmt(p.c14e, 3)}) at ${fmt(nv.x / 1000, 3)} km.`); V('Erosional velocity ratio', deg.erosionalRatio, 1, `${fmt(nv.x / 1000, 3)} km`); recs.push(`Cut the rate by ${fmt(100 * (1 - 1 / deg.erosionalRatio), 2)} % or enlarge the bore to ${fmt(p.idMm * Math.sqrt(deg.erosionalRatio), 3)} mm at ${fmt(nv.x / 1000, 3)} km to bring the velocity below the erosional limit of ${fmt(deg.ero[deg.iV].ve, 3)} m/s.`); }
  if (deg.eroMax > 0.1) { W(deg.eroMax > 0.5 ? 'bad' : 'warn', `Sand erosion reaches ${fmt(deg.eroMax, 3)} mm/y at ${deg.bWorst.name} (${deg.bWorst.ero.model}).`); recs.push(`Limit sand to ${fmt((p.sandKgD * 0.1) / deg.eroMax, 2)} kg/d or use a ${fmt(Math.max(deg.bWorst.rOverD * 2, 5), 2)}D bend at ${fmt(deg.bWorst.nd.x / 1000, 3)} km to keep erosion below 0.1 mm/y.`); }
  if (p.corrControl === 'availability' && p.inhibAvail > 95 && deg.corr.some((c) => !c.clad && !c.oilWet && c.bare > c.cri)) W('info', `The corrosion basis assumes ${fmt(p.inhibAvail, 3)} % inhibitor availability, which needs redundant injection pumps, automatic alarms and a shutdown on loss of injection; with an ordinary system (95 %) the rate would be ${fmt(Math.max(...deg.corr.map((c) => (c.clad || c.oilWet ? 0 : 0.95 * (1 - c.strip) * c.cri + (1 - 0.95 * (1 - c.strip)) * c.base))), 3)} mm/y.`);
  if (deg.cladLen > 0) W('info', `Corrosion-resistant alloy is used over ${fmt(deg.cladLen / 1000, 3)} km (${fmt(deg.cladFrom / 1000, 3)}–${fmt(deg.cladTo / 1000, 3)} km): the controlled carbon-steel rate there would be up to ${fmt(Math.max(...deg.corr.filter((c) => c.clad).map((c) => c.steel)), 3)} mm/y against the target of ${fmt(p.corrTarget, 2)} mm/y.`);
  { const hot = Math.max(0, ...deg.corr.filter((c) => !c.clad && !c.oilWet).map((c) => c.base)); if (p.corrControl !== 'none' && hot > 6) { W('warn', `The uninhibited corrosion rate of the carbon-steel sections reaches ${fmt(hot, 3)} mm/y; above about 6 mm/y inhibition is not accepted as the only barrier in operator practice (corrosion-resistant alloy or planned replacement).`); recs.push(`Clad or line the sections above 6 mm/y uninhibited rate (largest ${fmt(hot, 3)} mm/y) or raise the pH: inhibition alone is not relied on at that corrosivity.`); } }
  if (rectified) W('info', `The span of ${fmt(p.spanLengthEntered, 3)} m would not meet the fatigue criterion (design fatigue factor ${fmt(p.dff, 2)}, margin ${fmt(p.fatMargin, 2)}): supports are assumed at ${fmt(p.spanLength, 3)} m spacing.`);
  if (p.spanSplit) W('info', `A support of the line lies inside the span entered: the longest part, ${fmt(p.spanLength, 3)} m centred at ${fmt(p.spanX / 1000, 4)} km, is assessed.`);
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
  else if (lat.needed && !lat.managed) { W(lat.util > 1 ? 'bad' : 'warn', `Effective compression ${fmt(bk.comp / 1000, 3)} kN is ${fmt(bk.hobbsUtil, 3)} times the Hobbs buckling force (${bk.hb.governing}) and lateral buckling is not managed: a single buckle would take all the feed-in and bend the pipe to ${fmt(100 * lat.rogueStrain, 3)} % strain (allowable ${fmt(100 * lat.allow, 3)} %).`); if (lat.util > 1) V('Bending strain in an unplanned lateral buckle', lat.rogueStrain, lat.allow, `${fmt(bk.x / 1000, 3)} km`); recs.push(`Plan buckle initiators (sleepers or snake lay) about every ${fmt(clamp(2000 / Math.max(lat.util, 0.25), 500, 5000), 2)} m so that the expansion is shared: the compression of ${fmt(bk.comp / 1000, 3)} kN exceeds the ${fmt(bk.hb.critical / 1000, 3)} kN lateral buckling force.`); }
  else if (lat.needed && lat.util > 1) { const gm = lat.plan.governing; W('warn', `The planned lateral buckles are overstrained: ${fmt(100 * lat.strain, 3)} % bending strain (${gm.name}, feed-in ${fmt(gm.feedIn, 3)} m at ${fmt(lat.spacing, 4)} m spacing) against ${fmt(100 * lat.allow, 3)} % allowable.`); V('Bending strain in a planned lateral buckle', lat.strain, lat.allow, `${fmt(bk.x / 1000, 3)} km`); recs.push(`Reduce the spacing of the buckle initiators from ${fmt(lat.spacing, 4)} m to about ${fmt(Math.max(lat.spacing / lat.util ** 1.5, 200), 2)} m or lower the operating temperature: each planned buckle must stay below ${fmt(100 * lat.allow, 3)} % strain.`); }
  for (const d of def.rows) if (!d.ok) { W('bad', `${d.name} at ${fmt(d.x, 4)} m (${fmt(100 * d.ratio, 3)} % deep, ${fmt(d.length, 3)} mm long) fails Modified B31G: safe pressure ${fmt(d.safe / BAR, 3)} bar against ${fmt(d.dp / BAR, 3)} bar required.`); V(`${d.name}: estimated repair factor`, d.erf, 1, `${fmt(d.x, 4)} m`); recs.push(`Repair ${d.name.toLowerCase()} at ${fmt(d.x, 4)} m or de-rate to ${fmt((d.safe + d.nd.pe) / BAR, 3)} bara: its repair factor is ${fmt(d.erf, 3)}.`); }
  for (const d of def.rows) if (d.ok && d.life < remYears) recs.push(`${d.name} at ${fmt(d.x, 4)} m is acceptable today (repair factor ${fmt(d.erf, 3)}) but grows to the limit in ${fmt(d.life, 3)} y at ${fmt(d.rate, 3)} mm/y: re-inspect within ${fmt(Math.max(d.life / 2, 0.5), 2)} y.`);
  if (beta > 0.3 || plugged) W(beta > 0.6 ? 'bad' : 'warn', plugged ? `The line is plugged at ${fmt(ndB.x / 1000, 3)} km: upstream pressure rises to the shut-in pressure of ${fmt(p.pShutIn, 4)} bara.` : `Deposits take ${fmt(100 * beta, 3)} % of the flow area at ${fmt(ndB.x / 1000, 3)} km and add ${fmt(dpRestr / BAR, 3)} bar of pressure drop.`);
  if (p.pShutIn > mawpInlet) { W(plugged ? 'bad' : 'warn', `A full blockage would expose the line to the shut-in pressure of ${fmt(p.pShutIn, 4)} bara, above its allowable pressure of ${fmt(mawpInlet, 4)} bara at the inlet.`); V('Shut-in pressure against allowable pressure (bara)', p.pShutIn, mawpInlet, 'upstream of a plug'); recs.push(`Set the overpressure protection at ${fmt(mawpInlet, 3)} bara or below: the source can deliver ${fmt(p.pShutIn, 3)} bara against a plug.`); }
  if (plugDp > 5 * BAR && (beta > 0.3 || p.plugProb > 0.05)) recs.push(`Never depressurise a plug from one side: ${fmt(plugDp / BAR, 3)} bar across a ${fmt(p.plugLen, 3)} m plug could launch it at up to ${fmt(plugV, 3)} m/s (${fmt(plugE / 1e6, 3)} MJ) into the next bend.`);
  if (rel.pofNow > rel.target) { W(rel.pofNow > 10 * rel.target ? 'bad' : 'warn', `Annual probability of failure ${fmt(rel.pofNow, 2)} exceeds the target ${fmt(rel.target, 1)} of safety class ${p.sc}.`); V('Annual probability of failure', rel.pofNow, rel.target, rel.gov.name); }
  else if (!mg.met) { W(mg.max > 10 * rel.target ? 'bad' : 'warn', `With an inspection every ${fmt(mg.interval, 2)} y and repair at ${fmt(mg.threshold, 2)} mm of measured wall loss the annual failure probability still reaches ${fmt(mg.max, 2)} in year ${fmt(mg.tCross ?? rel.life, 3)}, above the target ${fmt(rel.target, 1)} of safety class ${p.sc} (${rel.gov.name}).`); V('Annual probability of failure with the inspection plan', mg.max, rel.target, rel.gov.name); recs.push(mg.auto ? `No inspection interval down to 1 y brings the annual failure probability below ${fmt(rel.target, 1)}: add wall thickness or corrosion allowance, raise the inhibitor availability above ${fmt(100 * p.eta, 3)} % or lower the repair threshold below ${fmt(mg.threshold, 2)} mm.` : `Shorten the inspection interval below ${fmt(mg.interval, 2)} y (or set it to 0 for the risk-based interval): the annual failure probability reaches ${fmt(mg.max, 2)} against the ${fmt(rel.target, 1)} target.`); }
  else recs.push(`Inspect the wall every ${fmt(mg.interval, 2)} y${mg.auto ? ' (risk-based interval)' : ''} and repair at ${fmt(mg.threshold, 2)} mm of measured loss: the annual failure probability then stays at or below ${fmt(mg.max, 2)} against the ${fmt(rel.target, 1)} target${rel.tTarget !== null ? `; without inspection it would cross the target in year ${fmt(rel.tTarget, 3)} (${rel.gov.name})` : ''}.`);
  if (topLevel >= 2) recs.push(`Highest risk: ${topThreat.name} at ${fmt(topThreat.risk / 1e6, 3)} M$/y (${fmt(topThreat.pof, 2)} per year × ${fmt(topThreat.cof / 1e6, 3)} M$): ${topThreat.action.toLowerCase()}.`);
  if (crackLife < remYears) { W('warn', `A ${fmt(p.flawMm, 2)} mm flaw grows to the critical depth of ${fmt(aC * 1000, 3)} mm in ${fmt(crackLife, 3)} y (Paris law).`); recs.push(`Tighten the weld acceptance flaw size below ${fmt(p.flawMm, 2)} mm or inspect the span welds every ${fmt(Math.max(crackLife / 3, 0.5), 2)} y: crack-growth life is ${fmt(crackLife, 3)} y.`); }
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

  // ---- wall models, tracked particles, equipment, vibration acceptance, measurements --------------------------------------
  prog(0.86, 'Wall finite elements, shell and tracked particles'); await tick();
  const cont = continuum(p, st, sa, wall, def, lossEol), tre = deg.mp > 0 && deg.bWorst.nd.vm - p.eroVth > 0 ? trackedErosion(p, deg) : null, eqWorst = eqp.length ? eqp.reduce((a, b) => (b.util > a.util ? b : a)) : null;
  if (cont.gr.collapse !== null && cont.gr.collapse < cont.dp) { W('bad', `The thinned wall (${fmt(cont.d * 1000, 3)} mm deep over ${fmt(cont.Lb * 1000, 3)} mm, ${cont.source}) reaches net-section collapse at ${fmt(cont.gr.collapse / BAR, 3)} bar differential pressure in the finite-element model, below the design value of ${fmt(cont.dp / BAR, 3)} bar.`); V('Net-section collapse of the thinned wall (bar)', cont.dp / BAR, cont.gr.collapse / BAR, `${fmt(cont.nd.x / 1000, 3)} km`); recs.push(`Repair or de-rate before the wall loss reaches ${fmt(cont.d * 1000, 3)} mm over ${fmt(cont.Lb * 1000, 3)} mm: the remaining ligament carries ${fmt(cont.gr.ligamentFactor, 3)} times the nominal hoop stress.`); }
  if (cont.ringUtil > 1) { W('warn', `The out-of-round${p.dentMm > 0 ? ' and dented' : ''} cross-section yields at the surface under the external overpressure of ${fmt(cont.pex / BAR, 3)} bar: peak hoop stress ${fmt(cont.hoopPeak / MPA, 3)} MPa (ring finite elements).`); recs.push(`Limit the out-of-roundness below ${fmt(p.ovality / cont.ringUtil, 2)} %${p.dentMm > 0 ? ' and cut out the dent' : ''}, or keep more than ${fmt(p.pMinShut, 3)} bara inside the line: the ring bending stress is ${fmt(cont.ringUtil, 3)} times yield.`); }
  if (span.vibUtil > 1) { const hard = span.velRms > span.vLim.problem && !(p.vibAllow > 0); W(hard ? 'bad' : 'warn', `The span vibrates at ${fmt(span.velRms, 3)} mm/s rms (${fmt(span.f1, 3)} Hz), above the allowable ${fmt(span.vAllow, 3)} mm/s${p.vibAllow > 0 ? '' : hard ? ' and above the problem level of the screening chart' : ' (concern level of the screening chart)'}.`); V('Vibration velocity (mm/s rms)', span.velRms, span.vAllow, `span at ${fmt(span.nd.x / 1000, 3)} km`); recs.push(`Stiffen or support the span to cut the vibration velocity by ${fmt(100 * (1 - 1 / span.vibUtil), 2)} %: ${fmt(span.velRms, 3)} mm/s rms against ${fmt(span.vAllow, 3)} mm/s allowed.`); }
  for (const e of eqp) if (e.util > 1) { W('bad', `${e.name}: the flow-induced load (${fmt(e.F / 1000, 3)} kN, ${fmt(e.M / 1000, 3)} kN·m) exceeds the allowable ${e.uF >= e.uM ? `force of ${fmt(e.fAllow, 3)} kN` : `moment of ${fmt(e.mAllow, 3)} kN·m`}.`); V(`Equipment loading: ${e.name}`, e.util, 1, `${fmt(e.xNode / 1000, 3)} km`); recs.push(`Add a support or anchor within ${fmt(e.arm / e.util, 2)} m of ${e.name.toLowerCase()} or reduce the slug velocity: its load utilisation is ${fmt(e.util, 3)}.`); }
  if (tre && tre.rateMax > 2 * Math.max(deg.bWorst.ero.gov, 1e-9) && tre.rateMax > 0.1) W('warn', `Particle tracking gives a local erosion peak of ${fmt(tre.rateMax, 3)} mm/y at ${fmt(tre.tr.angleMax, 3)}° of ${tre.bend.name.toLowerCase()}, more than twice the correlation value of ${fmt(deg.bWorst.ero.gov, 3)} mm/y: inspect that spot.`);
  { const gm = /X\s?(\d{2,3})/i.exec(String(p.material || '')), gr = gm ? LINEPIPE['X' + gm[1]] : null; if (gr && Math.abs(p.smysNew - gr[0]) > 0.05 * gr[0]) W('warn', `The yield strength entered (${fmt(p.smysNew, 4)} MPa) differs by more than 5 % from the specified minimum of grade X${gm[1]} (${gr[0]} MPa yield, ${gr[1]} MPa tensile): check the material inputs.`); }
  if (p.strengthLoss > 0) W('info', `Material condition: yield and tensile strength are reduced by ${fmt(p.strengthLoss, 3)} % (${fmt(p.smys, 4)} / ${fmt(p.smts, 4)} MPa) in every check.`);
  // measurements against the model
  const nearest = (list, x, key = 'x') => list.reduce((a, b) => (Math.abs(b[key] - x) < Math.abs(a[key] - x) ? b : a)), rangeMax = span.cycles.length ? Math.max(...span.cycles.map((c) => c.range)) : 2 * span.dynMax;
  const obsPred = (r) => { const n = st.at(r.x); switch (r.kind) {
    case 'pressure': return n.P; case 'temperature': return n.T; case 'strain': return strainOf(sa.rows[n.i].sOp).hoop * 1e6;
    case 'wall': return Math.max(p.wtMm - (deg.corr[n.i].rate + deg.ero[n.i].straight) * (r.cond > 0 ? r.cond : p.age), 0); case 'corrosion': return deg.corr[n.i].rate; case 'erosion': return nearest(deg.bends.map((b) => ({ x: b.nd.x, v: b.ero.gov })), n.x).v;
    case 'force': return nearest(sl.bends, n.x).force / 1000; case 'pulsation': return sl.pulsation; case 'frequency': return span.f1; case 'velocity': return span.velRms; case 'acceleration': return span.accRms; case 'stress': return rangeMax;
    case 'fatigue': return Math.min(snCycles(r.cond, p.snClass, p.snEnv, { ...span.snOpt, shift: p.snShift + 0.4 }), 1e15); case 'tensile': return rambergStress(r.cond / 100, p.eMod, p.yieldBias * p.smys, p.hardenN);
    case 'burst': return b31g({ D: p.D, t: p.t, d: Math.min(r.x / 1000, 0.99 * p.t), L: Math.max(r.cond, 1) / 1000, smys: p.S, smts: p.Su, method: 'dnv' }).pf / MPA; case 'collapse': return collapsePressure({ D: p.D, t: p.t, E: p.E, nu: p.nu, fy: p.yieldBias * p.S, ovality: r.cond > 0 ? r.cond / 100 : p.f0, alphaFab: p.alphaFab }).pc / MPA;
    case 'failure': return interp1(mg.t, mg.annual, clamp(r.cond, mg.t[0], mg.t[mg.t.length - 1])); default: return null; } };
  const obsRows = p.obs.map((r) => ({ ...r, pred: obsPred(r) })).filter((r) => isNum(r.pred)), obsKinds = [...new Set(obsRows.map((r) => r.kind))].map((kind) => { const rs = obsRows.filter((r) => r.kind === kind), n = rs.length, bias = sum(rs.map((r) => r.pred - r.value)) / n, rmse = Math.sqrt(sum(rs.map((r) => (r.pred - r.value) ** 2)) / n), nzr = rs.filter((r) => r.value !== 0), mape = nzr.length ? (100 * sum(nzr.map((r) => Math.abs((r.pred - r.value) / r.value)))) / nzr.length : null; return { kind, n, bias, rmse, mape, brier: kind === 'failure' ? brierScore(rs.map((r) => r.pred), rs.map((r) => (r.value > 0 ? 1 : 0))).brier : null }; });
  const lineKm = (N[N.length - 1].x - N[0].x) / 1000, genericPof = (p.failRate * lineKm) / 1000, genLife = genericPof * Math.max(remYears, 1);

  // ---- plots ----------------------------------------------------------------------------------------------------------
  const rows = sa.rows, Tg = rel.Tg, years = linspace(0, Math.max(remYears, 1) * 1.25, 40), wallAt = (y, q) => Math.max(p.wtMm - lossNow - rul.rate * (q === 0.5 ? 1 : randomVariable({ dist: 'lognormal', mean: 1, cov: rul.rateSd / rul.rate }).x(PhiInv(q))) * y, 0);
  plots.push({ type: 'line', title: 'Stress utilisation along the line', xlabel: 'Distance (km)', ylabel: 'Utilisation (–)', zeroY: true, series: [{ name: 'Pressure containment, design case', x: xkm, y: rows.map((r) => r.hoopDes) }, { name: 'Pressure containment, operating', x: xkm, y: rows.map((r) => r.hoopOp), dash: true }, { name: 'Equivalent stress, design case', x: xkm, y: rows.map((r) => r.vmDes) }, { name: 'Equivalent stress, operating', x: xkm, y: rows.map((r) => r.vmOp), dash: true }, { name: 'External collapse (depressurised)', x: xkm, y: rows.map((r) => r.collapse) }], hlines: [{ y: 1, label: 'limit', color: '#dc2626' }], note: `${rows[0].basis}. Equivalent stress against ${p.code === 'dnv' ? 'the derated yield strength (reference check)' : '0.90 SMYS'}.` });
  { const dmax = Math.max(1.3 * Math.max(...N.map((n) => n.depth)), 300), dd = linspace(0, dmax, 40), cE = sa.colOf(Math.max(p.t - p.CA / 1000, 0.05 * p.t));
    plots.push({ type: 'line', title: 'Collapse utilisation against water depth', xlabel: 'Water depth (m)', ylabel: 'Utilisation (–)', zeroY: true, series: [{ name: 'Collapse, present wall', x: dd, y: dd.map((d) => (RHO_SW * G * d * sa.gLB) / sa.col.pc) }, { name: 'Collapse, allowance consumed', x: dd, y: dd.map((d) => (RHO_SW * G * d * sa.gLB) / cE.pc), dash: true }, { name: 'Buckle propagation', x: dd, y: dd.map((d) => (RHO_SW * G * d * sa.gLB) / sa.col.ppr) }], hlines: [{ y: 1, label: 'limit', color: '#dc2626' }], vlines: p.waterDepth > 0 ? [{ x: Math.max(...N.map((n) => n.depth)), label: 'deepest point' }] : [], note: `Elastic ${fmt(sa.col.pel / BAR, 4)} bar, plastic ${fmt(sa.col.pp / BAR, 4)} bar, combined with ${fmt(100 * p.f0, 2)} % out-of-roundness ${fmt(sa.col.pc / BAR, 4)} bar; resistance factor ${fmt(sa.gLB, 4)}.` }); }
  plots.push({ type: 'line', title: 'Corrosion and erosion rate profiles', xlabel: 'Distance (km)', ylabel: 'Wall-loss rate (mm/y)', logy: true, series: [{ name: 'CO2 corrosion, managed', x: xkm, y: deg.corr.map((c) => Math.max(c.rate, 1e-6)) }, { name: 'CO2 corrosion, unmanaged', x: xkm, y: deg.corr.map((c) => Math.max(c.bare, 1e-6)), dash: true }, { name: 'Sand erosion, straight pipe', x: xkm, y: deg.ero.map((e) => Math.max(e.straight, 1e-6)) }, { name: 'Sand erosion at the bends', x: deg.bends.map((b) => b.nd.x / 1000), y: deg.bends.map((b) => Math.max(b.ero.gov, 1e-6)), mode: 'points' }], note: `Corrosion model: ${cg.model}; control basis: ${{ availability: `inhibitor availability ${fmt(p.inhibAvail, 3)} %`, efficiency: `inhibitor efficiency × availability ${fmt(p.inhibEff, 3)} %`, none: 'no inhibition' }[p.corrControl]}${deg.cladLen > 0 ? `; corrosion-resistant alloy over ${fmt(deg.cladLen / 1000, 3)} km (rate drawn at 1e-6)` : ''}. Values below 1e-6 mm/y are drawn at 1e-6.` });
  plots.push({ type: 'line', title: 'Wall thickness with time', xlabel: 'Years from today', ylabel: 'Wall thickness (mm)', series: [{ name: 'Median', x: years, y: years.map((y) => wallAt(y, 0.5)) }, { name: '10 % pessimistic', x: years, y: years.map((y) => wallAt(y, 0.9)), dash: true }, { name: '10 % optimistic', x: years, y: years.map((y) => wallAt(y, 0.1)), dash: true }], hlines: [{ y: p.wtMm - p.CA, label: 'allowance consumed', color: '#f59e0b' }, { y: Math.max(...rows.map((r) => r.tReq)) * 1000, label: 'pressure minimum', color: '#dc2626' }], vlines: [{ x: remYears, label: 'design life' }], note: `Rate ${fmt(rul.rate, 3)} ± ${fmt(rul.rateSd, 2)} mm/y${rul.updated ? ' (updated with the measured wall)' : ''}; remaining life P10 / P50 / P90 = ${fmt(rul.p10, 3)} / ${fmt(rul.p50, 3)} / ${fmt(rul.p90, 3)} y.` });
  plots.push({ type: 'line', title: 'Slug force on a 90° bend', xlabel: 'Time (s)', ylabel: 'Resultant force (kN)', zeroY: true, series: [{ name: ds.on ? 'Momentum force of the slug train' : 'Steady momentum force', x: sl.history.t, y: sl.history.f }], hlines: ds.on ? [{ y: sl.f90 / 1000, label: `with dynamic load factor ${fmt(p.dlf, 2)}`, color: '#dc2626' }] : [], note: ds.on ? `F = ρ·A·v²·√2 with slug-body density ${fmt(ds.rhoS, 3)} kg/m³ and film-region density ${fmt(ds.rhoF, 3)} kg/m³ at ${fmt(ds.v, 3)} m/s; lengths and periods scattered log-normally.` : 'No intermittent flow: steady momentum force of the mixture.' });
  plots.push({ type: 'line', title: 'Mode shapes of the span', xlabel: 'Position along the span (m)', ylabel: 'Normalised deflection (–)', series: span.modes.shapes.slice(0, 4).map((s, k) => ({ name: `Mode ${k + 1}: ${fmt(span.modes.f[k], 4)} Hz`, x: span.model.x.slice(), y: s.slice() })), note: `${span.model.x.length - 1} Hermite beam elements, effective mass ${fmt(span.me, 4)} kg/m (content and added mass included)${p.timoshenko ? ', shear deformation and rotary inertia included' : ''}.` });
  { const tt = span.resp ? thin(span.resp.t, 600) : [0, 1], ss = span.resp ? thin(span.resp.sigma.map((s) => s / MPA), 600) : [0, 0];
    plots.push({ type: 'line', title: 'Dynamic bending stress at the critical point of the span', xlabel: 'Time (s)', ylabel: 'Stress (MPa)', series: [{ name: span.resp ? `x = ${fmt(span.resp.x, 3)} m` : 'No slug excitation', x: tt, y: ss }], note: `Modal superposition of ${span.modes.f.length} modes, Newmark-β (average acceleration), Δt = ${fmt(span.dt, 2)} s, damping ${fmt(100 * p.zeta, 2)} %. Static sag stress ${fmt(span.sigStatic, 3)} MPa is not included.${span.note ? ' ' + span.note : ''}` }); }
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

  plots.push({ type: 'line', title: 'Annual probability of failure with and without the inspection plan', xlabel: 'Age of the line (y)', ylabel: 'Annual probability of failure (1/y)', logy: true, series: [{ name: 'No inspection', x: mg.t, y: mg.none.map((q) => Math.max(q, 1e-12)), dash: true }, { name: `Inspection every ${fmt(mg.interval, 2)} y with repair`, x: mg.t, y: mg.annual.map((q) => Math.max(q, 1e-12)) }], hlines: [{ y: rel.target, label: `target, safety class ${p.sc}`, color: '#dc2626' }, ...(genericPof > 0 ? [{ y: genericPof, label: 'generic statistics for this length', color: '#64748b' }] : [])], vlines: mg.inspections.map((r) => ({ x: p.age + r.t, label: '' })).slice(0, 12), note: `Series system of all limit states; values below 1e-12 are drawn at 1e-12. Vertical lines mark the inspections; the repair threshold is ${fmt(mg.threshold, 3)} mm of measured wall loss.` });
  plots.push({ type: 'line', title: 'Thinned wall: stress along the corroded surface', xlabel: 'Distance from the centre of the thinned band (mm)', ylabel: 'Stress (MPa)', zeroY: true, series: [{ name: 'Von Mises at the surface, axisymmetric solid elements', x: cont.gr.z.map((z) => z * 1000), y: cont.gr.vmSurface.map((q) => q / MPA) }, { name: 'Hoop membrane stress, thin-shell elements', x: cont.sh.x.map((z) => z * 1000), y: cont.sh.hoopMembrane.map((q) => q / MPA), dash: true }], hlines: [{ y: cont.gr.vmNominal / MPA, label: 'Lamé, intact wall', color: '#64748b' }, { y: cont.flow / MPA, label: 'flow stress', color: '#dc2626' }], note: `${fmt(cont.d * 1000, 3)} mm deep over ${fmt(cont.Lb * 1000, 3)} mm on the bore (${cont.source}) at ${fmt(cont.dp / BAR, 4)} bar differential pressure; ${cont.gr.fe.ndof} degrees of freedom, ${p.feOrder === 2 ? 'quadratic' : 'linear'} triangles.` });
  plots.push({ type: 'line', title: 'Hoop stress around the cross-section (ring model)', xlabel: 'Angle from the long axis (deg)', ylabel: 'Hoop stress (MPa)', series: [{ name: 'Inner surface', x: cont.ring.theta, y: cont.ring.hoopInner.map((q) => q / MPA) }, { name: 'Outer surface', x: cont.ring.theta, y: cont.ring.hoopOuter.map((q) => q / MPA), dash: true }], hlines: [{ y: cont.ring.hoopNominal / MPA, label: 'membrane', color: '#64748b' }], note: `Plane-strain ring with ${fmt(100 * p.f0, 2)} % out-of-roundness${p.dentMm > 0 ? ` and a ${fmt(p.dentMm, 3)} mm dent` : ''} under ${fmt(Math.max(cont.pex, 1) / BAR, 4)} bar external overpressure, present wall ${fmt(wall.tNow * 1000, 4)} mm. First-order bending: no amplification by the deflection.` });
  if (tre) {
    plots.push({ type: 'line', title: 'Erosion scar along the outside of the bend (tracked particles)', xlabel: 'Distance along the bend from its inlet (m)', ylabel: 'Erosion rate (mm/y)', zeroY: true, series: [{ name: `Tracked particles, ${tre.model === 'dnv' ? 'DNV-RP-O501' : tre.model === 'oka' ? 'Oka' : 'Finnie'} impact law`, x: tre.tr.s, y: tre.tr.rate.map((q) => q * p.erosionMult), mode: 'step' }], hlines: [{ y: tre.dnv, label: 'DNV-RP-O501 bend equation', color: '#dc2626' }], vlines: [{ x: 0, label: 'bend inlet' }, { x: tre.tr.field.R * tre.tr.field.angle, label: 'bend outlet' }], note: `${tre.tr.nPart} particles of ${fmt(p.sandUm, 3)} µm in the plane potential flow of ${tre.bend.name.toLowerCase()} (${fmt(tre.U, 3)} m/s, Stokes number ${fmt(tre.tr.stokes, 3)}); ${fmt(100 * tre.tr.hitFraction, 3)} % of them hit the outside wall, at a wear-weighted angle of ${fmt(tre.tr.meanAngle, 3)}°.` });
    const F = tre.tr.field, wallLine = (nv) => { const xs = [], ys = []; for (let i = 0; i <= 60; i++) { const q = F.toXY((F.length * i) / 60, nv); xs.push(q[0]); ys.push(q[1]); } return { x: xs, y: ys }; }, wo = wallLine(F.D / 2), wi = wallLine(-F.D / 2);
    plots.push({ type: 'line', title: 'Particle paths in the bend', xlabel: 'x (m)', ylabel: 'y (m)', series: [{ name: 'Outside wall (extrados)', x: wo.x, y: wo.y, color: '#0f172a' }, { name: 'Inside wall', x: wi.x, y: wi.y, color: '#0f172a' }, ...tre.tr.tracks.slice(0, 10).map((tk, i) => ({ name: i === 0 ? 'Particles' : `Particle ${i + 1}`, x: thin(tk.x, 80), y: thin(tk.y, 80), dash: true, color: '#2563eb' }))], note: 'Mid-plane of the bend; flow enters from the left. Stream function solved on a body-fitted grid, particles with Schiller–Naumann drag and wall restitution.' });
  }

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
    { x: bk.x, mechanism: bk.buried ? 'Upheaval buckling' : 'Lateral buckling', utilisation: bk.buried ? bk.palmerUtil : lat.needed ? (lat.managed ? lat.util : Math.max(lat.util, 1.0001)) : bk.hobbsUtil }, ...(eqWorst && eqWorst.util > 0 ? [{ x: eqWorst.xNode, mechanism: `Equipment loading: ${eqWorst.name}`, utilisation: eqWorst.util }] : []), ...(def.worst ? [{ x: def.worst.x, mechanism: 'Metal-loss defect', utilisation: def.worst.erf }] : []), ...(beta > 0 ? [{ x: ndB.x, mechanism: 'Blockage', utilisation: beta }] : []),
  ].map((c) => ({ x: sig(c.x, 5) ?? 0, mechanism: c.mechanism, utilisation: sig(cap(c.utilisation, 99), 4) ?? 0 })).sort((a, b) => b.utilisation - a.utilisation);
  tables.push({ title: 'Critical locations', columns: ['Distance (m)', 'Mechanism', 'Utilisation'], rows: crit.map((c) => [c.x, c.mechanism, c.utilisation]), note: 'Utilisation above 1 means the limit is exceeded; corrosion and erosion are measured against the allowance over the design life.' });
  tables.push({ title: 'Design and operating limit violations', columns: ['Limit', 'Value', 'Allowed', 'Where'], rows: viol.length ? viol : [['None', '—', '—', '—']] });
  tables.push({ title: 'Slug and flow-induced loads at the bends', columns: ['Bend', 'Distance (m)', 'Angle (deg)', 'Density (kg/m³)', 'Velocity (m/s)', 'Momentum force (kN)', 'Impulse per slug (kN·s)', 'Support utilisation', 'Steady force (kN)', 'Pressure thrust (kN)', 'Centrifugal load (kN/m)', 'ρv² (Pa)', 'Kinetic-energy band', 'Likelihood of failure'], rows: sl.bends.map((b) => [b.name, cell(b.x, 5), cell(b.angle), cell(b.rho), cell(b.v), cell(b.force / 1000), cell(b.impulse / 1000), p.supportCapacity > 0 ? cell(b.force / (p.supportCapacity * 1000), 3) : '—', cell(b.steady / 1000), cell(b.pressure / 1000), cell(b.centrifugal / 1000), cell(b.fiv.rhoV2), b.fiv.band, `${cell(b.fiv.lof, 3)} (${b.fiv.likelihood})`]),
    note: `Momentum force F = DLF·ρ·A·v²·√(2(1 − cos θ)); pressure thrust p·A·√(2(1 − cos θ)) is carried by the pipe wall unless there is an expansion joint. Pressure pulsation from the slug body: ±${fmt(sl.pulsation, 3)} bar (acceleration ${fmt(sl.dpAcc / BAR, 3)} bar, friction ${fmt(sl.dpFric / BAR, 3)} bar)${sl.flowAmplitude !== null ? `; the flow suite reports ±${fmt(sl.flowAmplitude, 3)} bar at the inlet` : ''}; unbalanced force on a straight run ${fmt(sl.runForce / 1000, 3)} kN. The last two columns are an Energy-Institute-style screening, not a detailed assessment.` });
  tables.push({ title: 'Span dynamics', columns: ['Quantity', 'Value', 'Unit'], rows: [['First natural frequency', cell(span.f1), 'Hz'], ['Higher modes', span.modes.f.slice(1, 4).map((f) => fmt(f, 4)).join(' / ') || '—', 'Hz'], ['Effective mass', cell(span.me), 'kg/m'], ['Submerged weight', cell(span.wSub), 'N/m'], ['Static sag', cell(span.stat.maxDeflection * 1000), 'mm'], ['Static bending stress', cell(span.sigStatic), 'MPa'],
    ['Peak dynamic bending stress', cell(span.dynMax), 'MPa'], ['Peak dynamic deflection', cell(span.dispMax * 1000), 'mm'], ['Slug passage time', cell(span.passage), 's'], ['Governing excitation', `${span.excitation.name}, ${fmt(span.excitation.f, 3)} Hz`, '—'], ['Excitation / first natural frequency', cell(span.fivRatio), '–'],
    ['Reduced velocity (current + waves)', cell(span.viv.vr), '–'], ['Onset in-line / cross-flow', `${fmt(span.viv.onsetIL, 3)} / ${fmt(span.viv.onsetCF, 3)}`, '–'], ['Stability parameter K_sd', cell(span.viv.ksd), '–'], ['Largest VIV amplitude in-line / cross-flow', `${fmt(span.vivWorst.aIL, 3)} / ${fmt(span.vivWorst.aCF, 3)}`, 'diameters'], ['Wave-induced velocity at the seabed', cell(span.uw), 'm/s'],
    ['Effective axial force (fully restrained)', cell(span.seffRaw / 1000), 'kN'], ['Effective axial force used at the span', cell(span.seff / 1000), 'kN'], ['Slug at the span: frequency / length / velocity', span.slug.on ? `${fmt(span.slug.freq, 3)} Hz / ${fmt(span.slug.len, 3)} m / ${fmt(span.slug.v, 3)} m/s` : 'no slugging', '—'], ['Sum of support reactions', cell(-span.stat.sumReactions / 1000), 'kN'], ['Applied static load', cell(span.stat.totalLoad / 1000), 'kN'],
    ['Critical (divergence) velocity of the conveyed fluid', span.crit.vc === null ? 'not conservative with these supports' : cell(span.crit.vc), 'm/s'], ['Flow velocity / critical velocity', span.crit.vc ? cell(span.vConv / span.crit.vc) : '—', '–'], ['First frequency with the flow / without', cell(span.fFlow / span.f1, 6), '–'],
    ['Peak dynamic stress: two-way / one-way coupling', span.cmp ? `${fmt(span.cmp.two / MPA, 4)} / ${fmt(span.cmp.one / MPA, 4)} (first ${fmt(span.cmp.window, 3)} s)` : '—', 'MPa'], ['Largest added modal mass in a slug', span.resp && span.twoWay ? cell(100 * span.resp.massRatio, 3) : '—', '%'],
    ['Vibration velocity: slug response / vortex shedding in the strongest current', `${fmt(span.velRms, 3)} / ${fmt(span.vivVel, 3)}`, 'mm/s rms'], ['Allowable vibration velocity: concern / problem', `${fmt(span.vLim.concern, 3)} / ${fmt(span.vLim.problem, 3)}`, 'mm/s rms'], ['Supports, guides, anchors and lumped masses on the span', p.midSupports.length ? p.midSupports.map((q) => `${q.type} at ${fmt(q.x, 3)} m`).join(', ') : 'none', '—']], note: 'Two-way coupling: pipe conveying fluid with the conveyed mass, Coriolis and centrifugal terms following the slug train of the one-dimensional flow model. Vortex-induced vibration follows DNV-RP-F105-type response models at screening level (single mode, Rayleigh-distributed current). On an exposed seabed the compression at the span is capped at the lateral buckling force, because the line releases any larger force by buckling.' });
  tables.push({ title: 'Stresses, strains and deformation at the most stressed point', columns: ['Quantity', 'Operating', 'Design case', 'Unit'], rows: [['Hoop stress at the bore (Lamé)', mp(rV.sOp.hoopInner), mp(rV.sDes.hoopInner), 'MPa'], ['Radial stress at the bore', mp(rV.sOp.radialInner), mp(rV.sDes.radialInner), 'MPa'], ['Longitudinal stress', mp(rV.sOp.longitudinal), mp(rV.sDes.longitudinal), 'MPa'], ['Von Mises equivalent', mp(rV.sOp.vonMises), mp(rV.sDes.vonMises), 'MPa'], ['Tresca equivalent', mp(rV.sOp.tresca), mp(rV.sDes.tresca), 'MPa'],
    ['Hoop strain (elastic)', cell(strainOf(rV.sOp).hoop * 1e6), cell(eDes.hoop * 1e6), 'µε'], ['Longitudinal strain (elastic)', cell(strainOf(rV.sOp).long * 1e6), cell(eDes.long * 1e6), 'µε'], ['Free thermal strain α·ΔT', cell(p.alphaT * (rV.nd.T - p.tInstall) * 1e6), cell(p.alphaT * dTdes * 1e6), 'µε'], ['Radial growth of the bore', cell(strainOf(rV.sOp).hoop * p.ID * 500), cell(eDes.hoop * p.ID * 500), 'mm'],
    ['Free end expansion per km if unrestrained', cell(p.alphaT * (rV.nd.T - p.tInstall) * 1e6), cell(p.alphaT * dTdes * 1e6), 'mm/km'], ['Span sag (static) / dynamic deflection', cell(span.stat.maxDeflection * 1000), cell(span.dispMax * 1000), 'mm']], note: `Point at ${fmt(rV.nd.x / 1000, 4)} km, ${rV.restrained ? 'axially restrained' : 'free to expand'}; linear-elastic (Hooke) strains from the three principal stresses.` });
  tables.push({ title: 'Wall shear and particle impact', columns: ['Quantity', 'Value', 'Unit'], rows: [['Largest wall shear stress', cell(tauMaxNode.tauW), 'Pa'], ['Its location', cell(tauMaxNode.x / 1000), 'km'], ['Wall shear / inhibitor-film limit', cell(tauMaxNode.tauW / Math.max(p.tauCrit, 1e-9)), '–'], ['Sand mass rate', cell(deg.mp * 1000), 'g/s'], ['Particle size / density', `${fmt(p.sandUm, 3)} / ${fmt(p.sandDensity, 4)}`, 'µm / kg/m³'], ['Particles per second', cell(deg.mp / Math.max(mPart, 1e-18)), '1/s'],
    ['Impact velocity at the worst fitting', cell(bw.nd.vm), 'm/s'], ['Characteristic impact angle', cell(eW.alpha), 'deg'], ['Kinetic energy of one particle', cell(0.5 * mPart * bw.nd.vm ** 2 * 1e6), 'µJ'], ['Impact power on the fitting', cell(0.5 * deg.mp * bw.nd.vm ** 2 * eW.G), 'W'], ['Share of particles that reach the wall (size correction)', cell(eW.G), '–'], ['Impact area', cell(eW.At * 1e4), 'cm²'], ['Angle function F(α) of the ductile wall', cell(eW.F), '–']], note: `Impact data at ${bw.name} (${fmt(bw.nd.x / 1000, 3)} km); wall shear from the flow picture.` });
  tables.push({ title: 'Exceedance events', columns: ['Event', 'Count', 'Basis'], rows: exceed.map((r) => [r[0], r[1], r[2]]), note: log.length ? 'Counts of the operating log refer to its rows.' : 'No operating log entered: add pressure and temperature history rows to count excursions and to add their cycles to the fatigue spectrum.' });
  tables.push({ title: 'Fatigue spectrum and damage', columns: ['Source', 'Stress range (MPa)', 'Cycles per year', 'Cycles to failure', 'Damage per year'], rows: span.spectrum.length ? span.spectrum.map((s) => { const Nf = snCycles(s.range, p.snClass, p.snEnv, { scf: p.scf, t: p.wtMm }); return [s.source, cell(s.range), cell(s.perYear), cell(Math.min(Nf, 1e30)), cell(s.perYear / Nf)]; }) : [['No cyclic loading', '—', '—', '—', '—']],
    note: `S–N class ${p.snClass}, ${{ air: 'air', cp: 'seawater with cathodic protection', free: 'free corrosion' }[p.snEnv]}, stress concentration ${fmt(p.scf, 3)}; Miner sum ${fmt(span.dYear, 3)} per year: calculated life ${fmt(fatCalc, 3)} y, ${fmt(fatLife, 3)} y with the design fatigue factor of ${fmt(p.dff, 2)}; without the design measures (class ${p.unmanagedWeld.snClass}, stress concentration ${fmt(p.unmanagedWeld.scf, 3)}, span ${fmt(p.spanLengthEntered, 3)} m) the calculated life is ${fmt(fatCalcUnm, 3)} y.${p.weldAuto ? ` Weld class from the hi-lo tolerance of ${fmt(p.hiLoMm, 2)} mm.` : ''}${spanAllow !== null ? ` Longest span that meets the criterion with a margin of ${fmt(p.fatMargin, 2)}: about ${fmt(spanAllow, 3)} m.` : ''} Crack-growth life from a ${fmt(p.flawMm, 2)} mm flaw: ${fmt(crackLife, 3)} y.` });
  tables.push({ title: 'Corrosion assessment at the governing point', columns: ['Quantity', 'Value', 'Unit'], rows: [['Location', cell(ndC.x / 1000), 'km'], ['Temperature', cell(ndC.T), '°C'], ['CO2 partial pressure / fugacity', `${fmt(deg.yCO2 * ndC.P, 3)} / ${fmt(cw.fCO2, 3)}`, 'bar'], ['pH of CO2-saturated water / used', `${fmt(cw.pHco2, 3)} / ${fmt(cw.pH, 3)}`, '–'], ['Scaling temperature', cell(cw.tScale), '°C'], ['Scale factor', cell(cw.Fscale), '–'], ['Glycol factor', cell(cw.Fglycol), '–'],
    ['Reaction-limited / mass-transfer-limited rate', `${fmt(cw.Vr, 3)} / ${fmt(cw.Vm, 3)}`, 'mm/y'], ['Model used', cg.model, '—'], ['pH used by the model', cell(cg.pH), '–'], ['Bicarbonate / salinity of the water', `${fmt(p.bicarb, 4)} mg/L / ${fmt(p.brineWt, 3)} wt % NaCl`, '—'], ...(cg.un ? [['NORSOK M-506 rate at the nearest valid condition', cell(cg.un.norsok), 'mm/y'], ['Mixed-potential ratio: temperature / salinity / pH and fugacity', `${fmt(cg.un.fT, 3)} / ${fmt(cg.un.fSalt, 3)} / ${fmt(cg.un.fChem, 3)}`, '–']] : []),
    ['Unmanaged rate (no inhibitor, no pH stabilisation, carbon steel)', cell(cg.bare), 'mm/y'], ['Largest unmanaged rate along the line', cell(deg.corrUnmanaged), 'mm/y'], ['Control basis', { availability: 'inhibitor availability', efficiency: 'inhibitor efficiency × availability', none: 'none' }[p.corrControl], '—'], ['Availability or efficiency applied', cell(100 * cg.eta), '%'], ['Inhibited rate assumed', p.corrControl === 'availability' ? cell(cg.cri) : '—', 'mm/y'], ['pH stabilisation', p.phStabOn ? `on, pH ${fmt(p.phStabTarget, 3)}` : 'off', '—'],
    ['Corrosion-resistant alloy', deg.cladLen > 0 ? `${fmt(deg.cladLen / 1000, 3)} km (${fmt(deg.cladFrom / 1000, 3)}–${fmt(deg.cladTo / 1000, 3)} km)` : 'none', '—'], ['Managed rate at this point', cell(cg.rate), 'mm/y'], ['Largest managed rate along the line', cell(deg.corrMax), 'mm/y'], ['The same with 95 % / 90 % availability', p.corrControl === 'availability' ? `${fmt(0.95 * cg.cri + 0.05 * cg.base, 3)} / ${fmt(0.9 * cg.cri + 0.1 * cg.base, 3)}` : '—', 'mm/y'], ['Water wetting', deg.corr[deg.iC].oilWet ? 'oil-wet' : `water-wet (water cut ${fmt(100 * ndC.wcut, 3)} %, liquid ${fmt(ndC.vL, 3)} m/s)`, '—'],
    ['Mixed-potential model: corrosion potential', cell(mix.Ecorr), 'V'], ['Mixed-potential model: current density', cell(mix.icorr), 'A/m²'], ['Mixed-potential model: rate (Faraday)', cell(mix.rate), 'mm/y'], ['Limiting current H⁺ / H2CO3', `${fmt(mix.iLimH, 3)} / ${fmt(mix.iLimH2CO3, 3)}`, 'A/m²'], ['Mass-transfer-limited rate (Sherwood)', cell(mix.rateLimit), 'mm/y'],
    ['H2S partial pressure', cell(pH2S), 'kPa'], ['Sour-service domain', sour.label, '—'], ['Top-of-line condensation (max)', cell(tlcMax), 'g/m²/s'], ['Wall today / lost', `${fmt(wall.tNow * 1000, 4)} / ${fmt(lossNow, 3)}`, 'mm'], ['Cumulative loss at the end of the design life', cell(lossEol), 'mm'], ['Remaining wall at the end of the design life', cell(p.wtMm - lossEol), 'mm'], ['Life to the pressure-minimum wall', cell(lifeMin), 'y']], note: 'The de Waard rows (scale factor, reaction and mass-transfer limits) are shown for comparison whichever model is selected. The mixed-potential model (Tafel kinetics, bare steel in the produced brine) is an independent cross-check; it does not include protective films. Availability basis: rate = A × inhibited rate + (1 − A) × uninhibited rate.' });
  tables.push({ title: 'Sand erosion at fittings', columns: ['Fitting', 'Distance (m)', 'Velocity (m/s)', 'R/D', 'DNV-RP-O501 (mm/y)', 'Salama (mm/y)', 'Oka (mm/y)', 'Finnie (mm/y)', 'Governing (mm/y)', 'Corrosion there (mm/y)'], rows: deg.bends.map((b) => [b.name, cell(b.nd.x, 5), cell(b.nd.vm), cell(b.rOverD, 3), cell(b.ero.dnv), cell(b.ero.salama), cell(b.ero.oka), cell(b.ero.finnie), cell(b.ero.gov), cell(b.corr)]),
    note: `${fmt(p.sandKgD, 3)} kg/d of ${fmt(p.sandUm, 3)} µm sand. API RP 14E erosional velocity ratio ${fmt(deg.erosionalRatio, 3)} (C = ${fmt(p.c14e, 3)}). Oka and Finnie are single-impact models applied over the DNV impact area; straight-pipe erosion peaks at ${fmt(Math.max(...deg.ero.map((e) => e.straight)), 2)} mm/y.` });
  tables.push({ title: 'Defect assessment', columns: ['Defect', 'Location (m)', 'Depth (mm)', 'Length (mm)', 'Depth / wall', 'B31G burst (bar)', 'Modified B31G burst (bar)', 'DNV-RP-F101 burst (bar)', 'Safe pressure (bar)', 'Required (bar)', 'Repair factor', 'Verdict', 'Life to limit (y)', 'Growth rate (mm/y)', 'Depth at end of life (mm)'], rows: def.rows.length ? def.rows.map((d) => [d.name, cell(d.x, 5), cell(d.depth), cell(d.length), cell(d.ratio, 3), cell(d.pf.b31g / BAR), cell(d.pf.modified / BAR), cell(d.pf.dnv / BAR), cell(d.safe / BAR), cell(d.dp / BAR), cell(d.erf, 3), d.ok ? 'acceptable' : 'repair or de-rate', cell(d.life, 3), cell(d.rate, 3), cell(Math.min(d.depth + d.rate * remYears, p.wtMm), 3)]) : [['No defects entered', '—', '—', '—', '—', '—', '—', '—', '—', '—', '—', '—', '—', '—', '—']],
    note: `Differential (internal minus external) burst pressures of the nominal wall. Safe pressure = factor × Modified B31G burst pressure${mapS ? `; the wall map holds ${mapS.readings} readings, minimum wall ${fmt(mapS.minMm, 4)} mm` : ''}.` });
  tables.push({ title: 'Blockage and plug screening', columns: ['Quantity', 'Value', 'Unit'], rows: [['Flow area lost', cell(100 * beta), '%'], ['Effective bore', cell(idEff * 1000), 'mm'], ['Time to the blockage limit at the deposit growth rate', daysToLimit === null ? 'no growth entered' : cell(daysToLimit), 'd'], ['Effective roughness used for the slug friction', cell(p.roughUm), 'µm'], ['Location', cell(ndB.x / 1000), 'km'], ['Extra pressure drop of the restriction', cell(dpRestr / BAR), 'bar'], ['Upstream pressure', cell(pUp), 'bara'], ['Allowable pressure at the inlet', cell(mawpInlet), 'bara'], ['Shut-in pressure / allowable', cell(p.pShutIn / Math.max(mawpInlet, 1e-9)), '–'], ['Differential pressure across a plug', cell(plugDp / BAR), 'bar'], ['Force on the plug', cell((plugDp * p.Ai) / 1000), 'kN'], ['Plug mass', cell(plugMass), 'kg'], ['Run to the next bend', cell(runLen), 'm'], ['Projectile velocity (upper bound)', cell(plugV), 'm/s'], ['Projectile energy', cell(plugE / 1e6), 'MJ']], note: 'The projectile estimate assumes the full differential pressure acts over the run without gas expansion losses or wall friction: an upper bound for procedures, not a design load.' });
  tables.push({ title: 'Reliability methods on the governing limit state', columns: ['Method', 'Reliability index', 'Probability of failure', 'Evaluations', 'Note'], rows: [['FORM (HL–RF)', cell(rel.f.beta), cell(rel.f.pf, 3), rel.f.evals, rel.f.converged ? `converged in ${rel.f.iterations} iterations` : 'not converged'], ['SORM (Breitung)', cell(rel.so.beta), cell(rel.so.pf, 3), '—', rel.so.valid ? `curvatures ${rel.so.kappa.map((k) => fmt(k, 2)).join(', ') || '—'}` : 'curvature correction not applicable'],
    ['Monte Carlo', cell(rel.mc.beta), cell(rel.mc.pf, 3), rel.mc.n, rel.mc.failures ? `c.o.v. ${fmt(rel.mc.cov, 2)}` : 'no failure sampled'], ['Importance sampling', cell(rel.is.beta), cell(rel.is.pf, 3), rel.is.n, rel.is.cov !== null ? `c.o.v. ${fmt(rel.is.cov, 2)}` : '—'], ['Latin hypercube', cell(rel.lh.beta), cell(rel.lh.pf, 3), rel.lh.n, rel.lh.failures ? `${rel.lh.failures} failures` : 'no failure sampled'], ['Response surface + FORM', cell(rel.rs.beta), cell(rel.rs.pf, 3), rel.rs.evals, `sampling on the surface: ${fmt(rel.rs.pfSampling, 3)}`]],
    note: `${rel.gov.name}, evaluated at year ${fmt(rel.life, 3)}. Target annual failure probability for safety class ${p.sc}: ${fmt(rel.target, 1)} (β = ${fmt(-PhiInv(rel.target), 3)}).` });
  tables.push({ title: 'Random-variable sensitivities', columns: ['Variable', 'Distribution', 'Mean', 'Coefficient of variation', 'Design point', 'α', 'α²'], rows: rel.gov.vars.map((vv, i) => { const rv = randomVariable(vv); return [rv.name, rv.dist, cell(rv.mean), cell(rv.mean ? rv.sd / Math.abs(rv.mean) : 0, 3), cell(rel.f.x[i]), cell(rel.f.alpha[i], 3), cell(rel.f.alpha[i] ** 2, 3)]; }), note: 'α is the unit vector from the mean point to the design point in standard-normal space; a negative α is a resistance variable.' });
  tables.push({ title: 'Limit states: annual probability of failure', columns: ['Limit state', 'β today', 'Annual pf today', 'β at end of life', 'Annual pf at end of life', 'Annual pf, life average'], rows: [...rel.LS.map((l) => [l.name, cell(l.beta[0], 3), cell(l.annualNow, 3), cell(l.beta[Tg.length - 1], 3), cell(l.annualEol, 3), cell(l.annualAvg, 3)]), ['Series system', cell(-PhiInv(clamp(rel.pofNow, 1e-23, 1 - 1e-12)), 3), cell(rel.pofNow, 3), cell(-PhiInv(clamp(rel.pofEol, 1e-23, 1 - 1e-12)), 3), cell(rel.pofEol, 3), cell(rel.pofAvg, 3)]], note: 'Without inspection. The reliability target applies to every single year; the risk models and the published probability of failure use the plan below.' });
  tables.push({ title: 'Inspection and repair plan: annual probability of failure', columns: ['Limit state', 'Wall-loss rate (mm/y)', 'Life average, no inspection', 'Life average, with the plan', 'Largest year, with the plan'], rows: [...rel.LS.map((l) => [l.name, l.wl ? cell(l.wl.rate, 3) : '—', cell(l.annualAvg, 3), cell(l.managedAvg ?? l.annualAvg, 3), cell(l.managedMax ?? Math.max(...l.annualT), 3)]), ['Series system', '—', cell(mg.noneAvg, 3), cell(mg.avg, 3), cell(mg.max, 3)]],
    note: `Inspection every ${fmt(mg.interval, 2)} y${mg.auto ? ' (risk-based: the longest of 10, 8, 6, 5, 4, 3, 2, 1 y that meets the target in every year)' : ''}; detection POD(d) = 1 − exp(−d/${fmt(mg.podLam, 3)} mm), depth sizing ± ${fmt(mg.sizing, 3)} mm (one standard deviation), repair when the measured loss exceeds ${fmt(mg.threshold, 3)} mm; ${fmt(mg.repairs, 2)} repairs expected over ${mg.years} y. Target ${fmt(rel.target, 1)} per year: ${mg.met ? 'met in every year' : `exceeded from year ${fmt(mg.tCross ?? 0, 3)}`}. The conditional failure probability is computed by FORM as a function of the wall loss and integrated over the wall-loss factor (40 equally likely quantiles).` });
  if (mg.inspections.length) tables.push({ title: 'Bayesian updating of the wall-loss factor at the inspections', columns: ['Inspection at year', 'Probability of a repair', 'Posterior mean of the factor (no repair so far)', 'Posterior scatter (c.o.v.)'], rows: mg.inspections.map((r) => [cell(p.age + r.t, 3), cell(r.pRepair, 3), cell(r.mean, 4), cell(r.cov, 3)]), note: `Prior: lognormal factor with mean ${fmt(rel.LS[0].wl?.mean ?? 1, 3)} and scatter ${fmt(100 * (rel.LS[0].wl?.cov ?? cC), 3)} % on the wall-loss rate (${rel.LS[0].name}). Every inspection that finds less than the repair threshold removes the fast-corroding tail, so the scatter falls.` });
  if (rel.fits.length) tables.push({ title: 'Distributions fitted by maximum likelihood', columns: ['Variable', 'Distribution', 'Samples', 'Mean', 'Standard deviation', 'Scatter (c.o.v.)', 'Log-likelihood'], rows: rel.fits.map((q) => [q.label, q.dist, q.n, cell(q.mean), cell(q.sd), cell(q.cov, 3), cell(q.logLik)]), note: 'These fits replace the assumed mean and scatter of the corresponding random variable in every limit state.' });
  tables.push({ title: 'Threats, risk and FMECA ranking', columns: ['No.', 'Threat', 'Cause', 'Effect', 'Annual probability', 'Consequence (M$)', 'Risk (M$/y)', 'Risk level', 'ALARP band', 'Severity', 'Occurrence', 'Detection', 'RPN', 'Action'], rows: threats.map((t, i) => [i + 1, t.name, t.cause, t.effect, cell(t.pof, 3), cell(t.cof / 1e6), cell(t.risk / 1e6, 3), RISK_NAMES[t.level], ['broadly acceptable', 'ALARP: reduce if practicable', 'ALARP: reduce unless grossly disproportionate', 'intolerable'][t.level], t.S, t.O, t.det, t.rpn, t.action]).sort((a, b) => b[12] - a[12]), note: `Probabilities are averages over the remaining design life; the wall-loss mechanisms are evaluated with the inspection-and-repair plan (every ${fmt(inspInt, 2)} y). Risk priority number = severity × occurrence × detection (each 1–10); severity and occurrence follow the matrix categories.` });
  tables.push({ title: 'Fault tree: minimal cut sets for loss of containment', columns: ['Cut set', 'Order', 'Probability per year', 'Share of the top event'], rows: ft.cutSets.map((c) => [c.events.join(' AND '), c.events.length, cell(c.p, 3), cell(ft.rareEvent > 0 ? c.p / ft.rareEvent : 0, 3)]), note: `Top event ${fmt(ft.top, 3)} per year (exact, independent basic events); rare-event sum ${fmt(ft.rareEvent, 3)}, upper bound ${fmt(ft.upperBound, 3)}. Most important basic event: ${ft.importance.reduce((a, b) => (b.fussellVesely > a.fussellVesely ? b : a)).name}.` });
  tables.push({ title: 'Event tree: outcomes of a leak', columns: ['Path', 'Outcome', 'Conditional probability', 'Frequency (1/y)', 'Consequence (M$)', 'Risk (M$/y)'], rows: etRows.map((o) => [o.label, o.outcome, cell(o.p, 3), cell(o.freq, 3), cell(o.cost / 1e6), cell((o.freq * o.cost) / 1e6, 3)]), note: `Initiating frequency ${fmt(ft.top, 3)} per year from the fault tree; expected consequence ${fmt(cLoc / 1e6, 4)} M$ per event (repair ${fmt(p.repairCost, 3)} M$, ${fmt(p.downtimeDays, 3)} d of production at ${fmt(prodDay / 1e6, 3)} M$/d, environmental and safety costs by path).` });
  tables.push({ title: 'Bow-tie summary', columns: ['Threat', 'Preventive barriers', 'Top event', 'Mitigating barriers', 'Consequences'], rows: [
    ['Internal corrosion', `Inhibition (${fmt(100 * p.eta, 3)} %), corrosion allowance ${fmt(p.corrAllow, 2)} mm, inspection every ${fmt(inspInt, 2)} y`, 'Loss of containment', `Leak detection (${fmt(p.pDetect, 2)}), isolation (${fmt(p.pIsolate, 2)})`, 'Release, repair, deferred production'], ['Sand erosion', 'Sand control, velocity limit, long-radius bends', 'Loss of containment', 'Wall-thickness monitoring, isolation', 'Release at fittings'],
    ['Fatigue (slugging, VIV)', 'Span correction, slug control, weld quality', 'Loss of containment', 'Leak detection, isolation', 'Crack, rupture'], ['Overpressure / plug', `Protection system (fails ${fmt(p.pfdProtect, 2)} on demand), hydrate management`, 'Loss of containment', 'Emergency shutdown, ignition control', 'Rupture, projectile'], ['Collapse / buckling', 'Wall thickness, minimum internal pressure, buckle arrestors', 'Loss of integrity', 'Arrestors limit the damaged length', 'Flattened section, replacement']] });
  tables.push({ title: 'Bayesian network: corrosion chain with inspection evidence', columns: ['Evidence', 'P(evidence)', 'P(inhibitor working)', 'P(high corrosion)', 'P(severe wall loss)', 'P(leak)'], rows: bn.map((b) => [b.label + (b === bnNow ? ' ← current state of knowledge' : ''), cell(b.r.pEvidence, 3), cell(b.r.posterior['Inhibitor working'].yes, 3), cell(b.r.posterior.Corrosion.high, 3), cell(b.r.posterior['Wall loss'].severe, 3), cell(b.r.posterior.Leak.yes, 3)]), note: 'Water wetting and inhibitor state → corrosion → wall loss → leak, with the inspection result as a child of wall loss. The corrosion probabilities come from the corrosion model and its uncertainty; the wall-loss and leak tables are engineering priors.' });

  { const pg = lat.plan.governing, rg = lat.rogue.governing;
    tables.push({ title: 'Design measures: buckle arrestors and lateral-buckling management', columns: ['Measure', 'Basis', 'Value', 'Limit', 'Utilisation', 'Status'], rows: [
      ['Buckle propagation without arrestors', '35·f_y·α_fab·(t/D)^2.5 with the resistance factor', cell(arr.pe / BAR) + ' bar outside', cell(sa.col.ppr / sa.gLB / BAR) + ' bar', cell(sa.propUtil), arr.required ? `a buckle would run through ${fmt(arr.unconfined, 4)} m of line deeper than ${fmt(arr.depthLimit, 4)} m` : 'a buckle does not propagate'],
      ['Buckle arrestors', `integral, ${fmt(p.arrThick, 3)} × wall, ${fmt(p.arrLength, 3)} m long (${arr.source})`, arr.provided ? `${arr.count} at ${fmt(arr.spacing, 4)} m` : 'none', arr.required ? 'required' : 'not required', arr.provided ? cell(arr.util) : '—', !arr.required ? 'not needed' : arr.provided ? (arr.util <= 1 ? `arrestors provided: a propagating buckle is confined to one spacing (${fmt(arr.confined, 4)} m)` : 'crossover pressure too low') : 'missing'],
      ['Arrestor crossover pressure', 'p_X = p_pr + (p_pr,BA − p_pr)·[1 − exp(−20·t₂·L/D²)]', cell((arr.pe * 1.1 * sa.gLB) / BAR) + ' bar factored (1.1·γ_m·γ_SC)', cell(arr.pX / BAR) + ' bar', cell(arr.util), arr.util <= 1 ? 'holds' : 'insufficient'],
      ['Residual consequence of a collapse', `${fmt(p.replaceCost, 3)} k$/m of pipe replaced plus repair and downtime`, cell(cBuckle / 1e6) + ' M$ confined', cell(cBuckleOpen / 1e6) + ' M$ unconfined', '—', arr.provided ? `the arrestors save ${fmt((cBuckleOpen - cBuckle) / 1e6, 3)} M$ per event` : '—'],
      ['Lateral buckling force (Hobbs)', bk.hb.governing, cell(bk.comp / 1000) + ' kN restrained', cell(bk.hb.critical / 1000) + ' kN', cell(bk.hobbsUtil), lat.needed ? 'the exposed line buckles laterally' : bk.buried ? 'buried line: upheaval check applies' : 'no buckling'],
      ['Planned buckle initiators', `sleepers or snake lay (${lat.source})`, lat.needed && lat.managed ? `${lat.count} at ${fmt(lat.spacing, 4)} m` : lat.needed ? 'none' : 'not required', '—', '—', !lat.needed ? 'not needed' : lat.managed ? 'expansion shared between the planned buckles' : 'unmanaged'],
      ['Feed-in to each planned buckle', pg ? pg.name : '—', pg ? cell(pg.feedIn) + ' m' : '—', rg ? cell(rg.feedIn) + ' m in a single unplanned buckle' : '—', '—', pg ? (pg.limited ? `limited by the spacing (slip length ${fmt(pg.slip, 4)} m each side)` : 'slip zones fully developed') : '—'],
      ['Post-buckle bending strain', 'Hobbs post-buckling moment k5·μ_l·w·L², ε = M·r/EI', pg ? cell(100 * pg.strain) + ' %' : '—', cell(100 * lat.allow) + ' %', pg ? cell(lat.util) : '—', !lat.needed ? '—' : lat.util <= 1 ? 'acceptable' : 'overstrained'],
      ['Strain capacity', p.strainAllow > 0 ? 'allowable strain from the inputs' : `DNV-ST-F101 displacement-controlled, ε_c/γ_ε with γ_ε = ${DNV_STRAIN[p.sc]}`, cell(100 * lat.capacity) + ' % capacity', cell(100 * lat.allow) + ' % allowable', '—', 'end-of-life wall, no internal overpressure credited'],
    ], note: 'Lateral buckles are assessed at the hottest restrained point with the fully restrained effective force; strain-based acceptance does not replace the fatigue and fracture checks of the buckle crown.' });
    if (lat.plan.modes.length) tables.push({ title: 'Planned lateral buckles (Hobbs post-buckling solution)', columns: ['Mode', 'Buckle length L (m)', 'Force in the buckle (kN)', 'Feed-in (m)', 'Slip length each side (m)', 'Amplitude (m)', 'Bending moment (kN·m)', 'Bending strain (%)'], rows: lat.plan.modes.map((m) => [m.name, cell(m.L), cell(m.P / 1000), cell(m.feedIn), cell(m.slip), cell(m.amplitude), cell(m.moment / 1000), cell(100 * m.strain)]), note: `Restrained force ${fmt(bk.comp / 1000, 4)} kN, submerged weight ${fmt(bk.w, 4)} N/m, friction ${fmt(p.muAxial, 2)} axial and ${fmt(p.muLateral, 2)} lateral, initiator spacing ${lat.managed ? fmt(lat.spacing, 4) + ' m' : 'unlimited'}.` }); }
  tables.push({ title: 'Wall models: continuum and shell finite elements', columns: ['Model', 'Quantity', 'Finite elements', 'Reference', 'Ratio'], rows: [
    ['Thick cylinder, axisymmetric solid', 'Hoop stress at the bore (MPa)', mp(cont.cyl.hoopBore), mp(cont.la.hoop) + ' (Lamé)', cell(cont.cyl.hoopBore / cont.la.hoop, 5)],
    ['Thick cylinder, axisymmetric solid', 'Hoop force ∫σθ dr (kN/m)', cell(cont.cyl.hoopForce / 1000, 6), cell((cont.pD * (p.D / 2 - p.t) - cont.nd.pe * (p.D / 2)) / 1000, 6) + ' (equilibrium)', cell(cont.cyl.hoopForce / (cont.pD * (p.D / 2 - p.t) - cont.nd.pe * (p.D / 2)), 7)],
    ['Thinned band, axisymmetric solid', 'Largest von Mises stress (MPa)', mp(cont.gr.vmMax), mp(cont.gr.vmNominal) + ' (intact wall)', cell(cont.gr.scf)],
    ['Thinned band, axisymmetric solid', 'Mean hoop stress of the ligament (MPa)', mp(cont.gr.hoopLigament), mp(cont.gr.hoopNominal) + ' (intact wall)', cell(cont.gr.ligamentFactor)],
    ['Thinned band, thin shell with a thickness step', 'Hoop membrane stress at the centre (MPa)', mp(cont.sh.hoopMembrane[0]), mp(cont.gr.hoopLigament) + ' (solid elements)', cell(cont.sh.hoopMembrane[0] / cont.gr.hoopLigament)],
    ['Thinned band: net-section collapse', 'Differential pressure (bar)', cont.gr.collapse === null ? '—' : cell(cont.gr.collapse / BAR), `${cell(cont.codes.b31g / BAR)} B31G / ${cell(cont.codes.modified / BAR)} Modified B31G / ${cell(cont.codes.dnv / BAR)} DNV-RP-F101`, cont.gr.collapse === null ? '—' : cell(cont.gr.collapse / cont.codes.modified)],
    ['Thinned band: margin on the design case', 'Collapse / design differential pressure', cont.gr.collapse === null ? '—' : cell(cont.gr.collapse / cont.dp), cell(cont.dp / BAR) + ' bar design', cont.gr.collapse !== null && cont.gr.collapse < cont.dp ? 'fails' : 'passes'],
    ['Ring with out-of-roundness' + (p.dentMm > 0 ? ' and a dent' : ''), 'Peak hoop stress (MPa)', mp(cont.hoopPeak), mp(cont.ringHand) + ' (p·r/t·(1 + 6·w₁/t), thin ring)', cell(cont.hoopPeak / cont.ringHand)],
    ['Ring with out-of-roundness' + (p.dentMm > 0 ? ' and a dent' : ''), 'Peak hoop stress / yield', cell(cont.ringUtil), '1', cont.ringUtil > 1 ? 'yields' : 'elastic'],
  ], note: `Linear elasticity (Cauchy equilibrium with Hooke's law) on ${p.feOrder === 2 ? 'quadratic six-node' : 'linear three-node'} triangles, ${p.feNr} cells through the wall. The thinned band runs all around the circumference, which is more severe than a local patch, so its net-section estimate (mean ligament hoop stress = flow stress ${fmt(cont.flow / MPA, 4)} MPa) is a lower bound to the code burst pressures of a patch of the same depth and length.` });
  if (tre) tables.push({ title: 'Sand erosion with tracked particles (worst bend)', columns: ['Quantity', 'Value', 'Unit'], rows: [['Bend', `${tre.bend.name}, ${fmt(tre.bend.angle, 3)}°, r/D = ${fmt(tre.bend.rOverD, 3)}`, '—'], ['Particle Stokes number', cell(tre.tr.stokes), '–'], ['Share of the sand that hits the outside wall', cell(100 * tre.tr.hitFraction), '%'], ['Wall impacts counted', tre.tr.impacts, '—'], ['Wear-weighted impact angle / speed', `${fmt(tre.tr.meanAngle, 3)}° / ${fmt(tre.tr.meanSpeed, 3)} m/s`, '—'],
    ['Peak erosion rate on the scar', cell(tre.rateMax), 'mm/y'], ['Position of the peak', cell(tre.tr.angleMax), '° of bend'], ['Mean erosion rate over the scar', cell(tre.rateMean), 'mm/y'], ['DNV-RP-O501 bend equation', cell(tre.dnv), 'mm/y'], ['DNV-RP-O501 without its model factor C₁ = 2.5', cell(tre.dnvNoC1), 'mm/y'], ['Tracked peak / DNV-RP-O501', tre.dnv > 0 ? cell(tre.rateMax / tre.dnv) : '—', '–'], ['Particle-size correction G of DNV-RP-O501', cell(tre.G), '–']],
    note: 'Plane potential flow on a body-fitted grid (stream function), Lagrangian particles with Schiller–Naumann drag, impact wear law at every wall impact, scar width π·D/4. The tracked value resolves where the sand hits; it carries no model factor, so it is normally below the recommended-practice value.' });
  if (eqp.length) tables.push({ title: 'Equipment and connection loads', columns: ['Item', 'Distance (m)', 'Bore (mm)', 'Velocity in the bore (m/s)', 'Flow-induced force (kN)', 'Allowable force (kN)', 'Moment (kN·m)', 'Allowable moment (kN·m)', 'Utilisation', 'Pressure thrust (kN)'], rows: eqp.map((e) => [e.name, cell(e.xNode, 5), cell(e.size || p.idMm), cell(e.velocity), cell(e.F / 1000), e.fAllow > 0 ? cell(e.fAllow) : '—', cell(e.M / 1000), e.mAllow > 0 ? cell(e.mAllow) : '—', cell(e.util), cell(e.thrust / 1000)]), note: `Force: slug momentum force with the dynamic load factor ${fmt(p.dlf, 2)} and the force coefficient ${fmt(p.slugCf, 3)} when slugging, otherwise the steady momentum force; moment = force × lever arm. The pressure thrust is carried by the pipe wall unless an expansion joint is fitted.` });
  tables.push({ title: 'Material condition and initial state', columns: ['Quantity', 'Value', 'Unit'], rows: [['Yield / tensile strength used', `${fmt(p.smys, 4)} / ${fmt(p.smts, 4)}`, 'MPa'], ['Loss of strength from ageing or damage', cell(p.strengthLoss), '%'], ['Derating at the design temperature (DNV format)', cell(dnvDerating(p.designTemp)), 'MPa'], ['Welding residual stress in the fracture check', cell(p.residualStress), 'MPa'], ['Girth-weld misalignment', cell(p.hiLoMm), 'mm'], ['Misalignment stress concentration factor', cell(p.scfWeld), '–'], ['Total stress concentration factor in fatigue', cell(p.scf), '–'],
    ['Strain at yield (Ramberg–Osgood, n = ' + fmt(p.hardenN, 3) + ')', cell(100 * rambergStrain(p.smys, p.eMod, p.smys, p.hardenN)), '%'], ['Fabrication factor α_fab', cell(p.alphaFab), '–'], ['Initial deposit thickness', cell(p.deposit0), 'mm'], ['Deposit mass on the span', cell(span.mDep), 'kg/m'], ['Wall loss at the start of the assessment', cell(lossNow), 'mm'], ['Fatigue damage already accumulated', cell(p.priorDamage), '–']] });
  tables.push({ title: 'Failure frequency against public statistics', columns: ['Quantity', 'Value', 'Unit'], rows: [['Line length', cell(lineKm), 'km'], ['Generic loss-of-containment frequency', cell(p.failRate), 'per 1000 km·y'], ['Generic frequency for this line', cell(genericPof, 3), '1/y'], ['Predicted loss-of-containment frequency (fault tree)', cell(ft.top, 3), '1/y'], ['Predicted structural failure probability, life average with the plan', cell(mg.avg, 3), '1/y'], ['Predicted / generic', genericPof > 0 ? cell(ft.top / genericPof, 3) : '—', '–'], ['Expected failures over the remaining life at the generic frequency', cell(genLife, 3), '–'], ['Probability of at least one failure: constant rate (Poisson)', cell(1 - Math.exp(-genLife), 3), '–'], ['The same with the rate uncertain between lines (gamma–Poisson, shape ' + fmt(FREQ_SHAPE, 3) + ')', cell(1 - (1 + genLife / FREQ_SHAPE) ** -FREQ_SHAPE, 3), '–'], ['90 % interval of the generic frequency for this line', `${fmt(genericPof * gammaQuantile(0.05, FREQ_SHAPE) / FREQ_SHAPE, 2)} – ${fmt(genericPof * gammaQuantile(0.95, FREQ_SHAPE) / FREQ_SHAPE, 2)}`, '1/y']],
    note: 'The generic value is a population average over many lines, ages and threat mixes; a new, inspected line should lie below it. Yearly incident counts of a pipeline population scatter more than a constant-rate (Poisson) process allows: the frequency model is a negative-binomial (gamma–Poisson) one with a trend, whose forecast intervals are checked against the public incident counts on the calibration tab; the shape used here is the one fitted to those counts.' });
  if (obsRows.length) { tables.push({ title: 'Measurements against the model', columns: ['Kind', 'Location or condition', 'Condition', 'Measured', 'Predicted', 'Difference'], rows: obsRows.map((r) => [r.kind, cell(r.x, 5), cell(r.cond), cell(r.value, 5), cell(r.pred, 5), cell(r.pred - r.value, 4)]) });
    tables.push({ title: 'Measurement comparison: error metrics', columns: ['Kind', 'Points', 'Bias (predicted − measured)', 'Root-mean-square error', 'Mean absolute error (%)', 'Brier score'], rows: obsKinds.map((k) => [k.kind, k.n, cell(k.bias), cell(k.rmse), k.mape === null ? '—' : cell(k.mape, 3), k.brier === null ? '—' : cell(k.brier, 3)]), note: 'Failure observations (1 or 0 per year) are scored against the predicted annual failure probability with the Brier score.' }); }

  // ---- additional solver of the set-up tab ------------------------------------------------------------------------------
  const task = p.task !== 'none' ? await runTask(p, st, ds, sl, span, sa, ctx, prog) : null;
  if (task) { plots.push(...task.plots); tables.push(...task.tables); warnings.push(...task.warnings); }
  // ---- KPIs, outputs --------------------------------------------------------------------------------------------------
  const kp = (label, value, unit, status, help) => ({ label, value: sig(value, 4) ?? 0, unit, status, help }), relBeta = clamp(-PhiInv(clamp(rel.pofNow, 1e-23, 1 - 1e-12)), -5, 10);
  const kpis = [
    kp('Pressure containment utilisation', sa.hoopUtil, '–', stat(sa.hoopUtil, 0.98), `${codeName}, worst of the design and operating cases; the code safety factors are inside the number, 1.0 is the limit`), kp('Equivalent stress utilisation', sa.vmUtil, '–', stat(sa.vmUtil, 0.98), 'Von Mises against the combined-stress allowable'), kp('Collapse utilisation', sa.collapseUtil, '–', stat(sa.collapseUtil, 0.98), 'External pressure on the depressurised pipe'),
    kp('Allowable pressure, present wall', sa.mawp, 'bara', p.designPressure > sa.mawp * 1.0001 ? 'bad' : p.designPressure > sa.mawpEol ? 'warn' : 'ok', 'Referred to the design-pressure datum'), kp('Minimum required wall', sa.minWall, 'mm', sa.minWall > p.wtMm ? 'bad' : 'ok', `Pressure containment plus allowances; ${fmt(p.wtMm - sa.minWall, 2)} mm in hand on the nominal wall`),
    kp('Peak force on a 90° bend', sl.peak90 / 1000, 'kN', 'ok', ds.on ? 'Slug momentum force with the dynamic load factor' : 'Steady momentum force (no slugging)'), kp('First natural frequency of the span', span.f1, 'Hz', span.fivRatio > 0.8 && span.fivRatio < 1.25 ? 'bad' : 'ok', 'Beam finite-element model'),
    kp('Excitation / natural frequency', span.fivRatio, '–', span.fivRatio > 0.8 && span.fivRatio < 1.25 ? 'bad' : span.fivRatio > 0.5 && span.fivRatio < 2 ? 'warn' : 'ok', `Governing excitation: ${span.excitation.name}`), kp('Fatigue life with the design fatigue factor', fatLife, 'y', fatLife < 0.5 * remYears ? 'bad' : fatLife < remYears ? 'warn' : 'ok', `Calculated life ${fmt(fatCalc, 3)} y divided by the design fatigue factor of ${fmt(p.dff, 2)}; weld class ${p.snClass}${rectified ? `, supports at ${fmt(p.spanLength, 3)} m` : ''}. Without the design measures the calculated life is ${fmt(fatCalcUnm, 3)} y`),
    kp('CO2 corrosion rate, managed', deg.corrMax, 'mm/y', deg.corrMax > 0.5 ? 'bad' : deg.corrMax > 0.1000001 ? 'warn' : 'ok', `Largest rate along the line with the corrosion-control basis${deg.cladLen > 0 ? ` (${fmt(deg.cladLen / 1000, 3)} km of corrosion-resistant alloy)` : ''}; unmanaged ${fmt(deg.corrUnmanaged, 3)} mm/y`), kp('Sand erosion rate', deg.eroMax, 'mm/y', deg.eroMax > 0.5 ? 'bad' : deg.eroMax > 0.1 ? 'warn' : 'ok', `Largest of bends and straight pipe (${deg.bWorst.ero.model})`),
    kp('Erosional velocity ratio', deg.erosionalRatio, '–', stat(deg.erosionalRatio, 0.8), 'API RP 14E'), kp('Remaining life on the allowance', remLife, 'y', remLife < 0.5 * remYears ? 'bad' : remLife < remYears ? 'warn' : 'ok', `At ${fmt(rate, 3)} mm/y combined wall loss`),
    kp('Annual probability of failure, today', rel.pofNow, '1/y', rel.pofNow > 10 * rel.target ? 'bad' : rel.pofNow > rel.target ? 'warn' : 'ok', `Structural limit states in series; target ${fmt(rel.target, 1)} for safety class ${p.sc}`), kp('Annual probability of failure, life average', mg.avg, '1/y', mg.met ? 'ok' : mg.max > 10 * rel.target ? 'bad' : 'warn', `With inspection every ${fmt(mg.interval, 2)} y and repair; largest in any year ${fmt(mg.max, 2)}. Without inspection: average ${fmt(mg.noneAvg, 2)}, largest year ${fmt(mg.noneMax, 2)}`), kp('Reliability index, today', relBeta, '–', rel.pofNow > rel.target ? 'warn' : 'ok', 'Annual, series system of the limit states'),
    kp('Risk cost', riskCost / 1e6, 'M$/y', topLevel >= 3 ? 'bad' : topLevel >= 2 ? 'warn' : 'ok', `Highest level: ${RISK_NAMES[topLevel]} (${topThreat.name})`), kp('Inspection interval', inspect, 'y', inspect < 1 ? 'warn' : 'ok', mg.auto ? 'Risk-based: the longest interval that keeps every year below the reliability target, capped by half the fatigue and defect lives' : 'As specified, capped by half the fatigue and defect lives'),
    kp('Lateral-buckle strain utilisation', lat.needed ? lat.util : 0, '–', lat.needed ? stat(lat.util, 0.8) : 'ok', lat.needed ? `${lat.managed ? `Planned buckles every ${fmt(lat.spacing, 4)} m` : 'Unplanned single buckle'}: ${fmt(100 * (lat.managed ? lat.strain : lat.rogueStrain), 3)} % against ${fmt(100 * lat.allow, 3)} % allowable` : 'The line does not buckle laterally'),
    kp('Net-section collapse / design pressure (thinned wall)', cont.gr.collapse !== null ? cont.gr.collapse / cont.dp : 99, '–', cont.gr.collapse !== null && cont.gr.collapse < cont.dp ? 'bad' : 'ok', `Axisymmetric finite elements, ${fmt(cont.d * 1000, 3)} mm deep over ${fmt(cont.Lb * 1000, 3)} mm (${cont.source}); stress concentration ${fmt(cont.gr.scf, 3)}`),
    kp('Span vibration velocity', span.velRms, 'mm/s rms', span.vibUtil > 1 ? (span.velRms > span.vLim.problem && !(p.vibAllow > 0) ? 'bad' : 'warn') : 'ok', `Allowable ${fmt(span.vAllow, 3)} mm/s at ${fmt(span.f1, 3)} Hz`),
  ];
  const balances = [{ name: 'Span static equilibrium: applied load against support reactions (N)', in: span.stat.totalLoad, out: -span.stat.sumReactions }, { name: 'Markov state probabilities at the end of the horizon', in: 1, out: sum(mk.p[mk.p.length - 1]) }, { name: 'Event-tree outcome frequencies against the initiating frequency (1/y)', in: ft.top, out: et.total },
    { name: 'Lamé invariant σr + σθ at the bore against the outer wall (MPa)', in: (rH.sDes.hoopInner + rH.sDes.radialInner) / MPA, out: (rH.sDes.hoopOuter + rH.sDes.radialOuter) / MPA }, { name: 'Bayesian network: posterior probabilities of the leak node', in: 1, out: bnNow.r.posterior.Leak.yes + bnNow.r.posterior.Leak.no },
    { name: 'Thick cylinder finite elements: hoop force ∫σθ dr against pᵢ·rᵢ − pₑ·rₒ (N/m)', in: cont.pD * (p.D / 2 - p.t) - cont.nd.pe * (p.D / 2), out: cont.cyl.hoopForce }, { name: 'Thick cylinder finite elements: external work against strain energy (J)', in: cont.cyl.fe.work, out: cont.cyl.fe.energy },
    { name: 'Thinned wall finite elements: radial pressure load against the hoop resultant (N/rad)', in: cont.gr.radialLoad / (2 * Math.PI), out: cont.gr.hoopResultant / (2 * Math.PI) }, { name: 'Thinned wall finite elements: external work against strain energy (J)', in: cont.gr.fe.work, out: cont.gr.fe.energy },
    { name: 'Ring finite elements: pressure resultant against the support reactions (N/m)', in: cont.ring.balance.applied, out: cont.ring.balance.reaction }, { name: 'Thin shell: pressure load against hoop and ring resistance (N/m)', in: cont.sh.radialLoad, out: cont.sh.radialResistance }];
  const outputs = {
    slugForce: sig(sl.peak90 / 1000, 5) ?? 0, hoopUtil: sig(sa.hoopUtil, 5) ?? 0, vmUtil: sig(sa.vmUtil, 5) ?? 0, collapseUtil: sig(sa.collapseUtil, 5) ?? 0, mawp: sig(sa.mawp, 5) ?? 0, minWallRequired: sig(sa.minWall, 5) ?? 0,
    corrosionRate: sig(deg.corrMax, 5) ?? 0, erosionRate: sig(deg.eroMax, 5) ?? 0, wallLossRate: sig(rate, 5) ?? 0, remainingLife: sig(remLife, 5) ?? 0, fatigueDamagePerYear: sig(span.dYear, 5) ?? 0, fatigueLife: sig(fatCalc, 5) ?? 0,
    fivRatio: sig(span.fivRatio, 5) ?? 0, naturalFrequency: sig(span.f1, 9) ?? 0, pof: sig(mg.avg, 5) ?? 0, reliabilityIndex: sig(clamp(-PhiInv(clamp(mg.avg, 1e-23, 1 - 1e-12)), -5, 10), 5) ?? 0, pofToday: sig(rel.pofNow, 5) ?? 0, riskLevel: RISK_NAMES[topLevel], riskCostPerYear: sig(riskCost, 5) ?? 0, consequence: sig(cLoc, 5) ?? 0, inspectionInterval: sig(inspect, 4) ?? 0,
    criticalLocations: crit.slice(0, 10), violations: viol.length, upheavalUtil: sig(bk.buried ? bk.palmerUtil : lat.needed ? (lat.managed ? lat.util : Math.max(lat.util, 1.0001)) : bk.hobbsUtil, 5) ?? 0, erosionalRatio: sig(deg.erosionalRatio, 5) ?? 0,
    // additional values
    mawpEndOfLife: sig(sa.mawpEol, 5) ?? 0, wallNow: sig(wall.tNow * 1000, 5) ?? 0, pofEndOfLife: sig(rel.pofEol, 5) ?? 0, targetPof: rel.target, lofFrequency: sig(ft.top, 5) ?? 0, slugFrequency: sig(span.slug.on ? span.slug.freq : ds.on ? ds.freq : 0, 5) ?? 0, peakDynamicStress: sig(span.dynMax, 5) ?? 0, staticSpanStress: sig(span.sigStatic, 5) ?? 0,
    pressurePulsation: sig(sl.pulsation, 5) ?? 0, crackGrowthLife: sig(crackLife, 5) ?? 0, remainingLifeP10: sig(rul.p10, 5) ?? 0, remainingLifeP90: sig(rul.p90, 5) ?? 0, localBucklingUtil: sig(lbOp.util, 5) ?? 0, propagationUtil: sig(sa.propUtil, 5) ?? 0, blockage: sig(beta, 5) ?? 0, sourService: sour.sour, hydrotestPressure: sig(hydro, 5) ?? 0,
    slugImpulse: sig(bendMax.impulse / 1000, 5) ?? 0, supportUtil: sig(supUtil, 5) ?? 0, cumulativeLoss: sig(lossNow, 5) ?? 0, cumulativeLossEndOfLife: sig(lossEol, 5) ?? 0, maxWallShear: sig(tauMaxNode.tauW, 5) ?? 0, exceedanceEvents: exceed.map((r) => ({ event: r[0], count: typeof r[1] === 'number' ? r[1] : 0 })), daysToBlockageLimit: daysToLimit === null ? null : sig(daysToLimit, 5),
    // design measures and the unmanaged counterparts
    corrosionRateUnmanaged: sig(deg.corrUnmanaged, 5) ?? 0, corrosionManaged: { basis: p.corrControl, model: cg.model, availability: p.corrControl === 'availability' ? sig(p.inhibAvail, 4) : null, efficiency: p.corrControl === 'efficiency' ? sig(p.inhibEff, 4) : null, inhibitedRate: p.corrControl === 'availability' ? sig(cg.cri, 4) : null, phStabilised: p.phStabOn, pH: sig(cg.pH, 4), bicarbonate: sig(p.bicarb, 5) ?? 0, salinity: sig(p.brineWt, 4) ?? 0, cladLength: sig(deg.cladLen, 6) ?? 0, cladFrom: deg.cladFrom === null ? null : sig(deg.cladFrom, 6), cladTo: deg.cladTo === null ? null : sig(deg.cladTo, 6), carbonSteelRate: sig(deg.steelMax, 5) ?? 0, target: sig(p.corrTarget, 4), rateAt95: p.corrControl === 'availability' ? sig(Math.max(...deg.corr.map((c) => (c.clad || c.oilWet ? 0 : 0.95 * (1 - c.strip) * c.cri + (1 - 0.95 * (1 - c.strip)) * c.base))), 5) : null },
    fatigueLifeUnmanaged: sig(fatCalcUnm, 5) ?? 0, fatigueLifeDesign: sig(fatLife, 5) ?? 0, fatigueUtil: sig(p.priorDamage + span.dYear * p.dff * remYears, 5) ?? 0, fatigueManaged: { weldClass: p.snClass, weldClassUnmanaged: p.unmanagedWeld.snClass, scf: sig(p.scf, 4), scfUnmanaged: sig(p.unmanagedWeld.scf, 4), hiLo: sig(p.hiLoMm, 3) ?? 0, spanLength: sig(p.spanLength, 5), spanEntered: sig(p.spanLengthEntered, 5), rectified, spanAllowable: spanAllow === null ? null : sig(spanAllow, 4), margin: p.fatMargin, designFatigueFactor: p.dff },
    fsi: task ? task.outputs.fsi ?? null : null, cfdFea: task ? task.outputs.cfdFea ?? null : null, solid: task ? task.outputs.solid ?? null : null, task: p.task,
    pofUnmanaged: sig(mg.noneAvg, 5) ?? 0, pofMaxYear: sig(mg.max, 5) ?? 0, pofUnmanagedMaxYear: sig(mg.noneMax, 5) ?? 0, pofTargetMet: mg.met, riskCostUnmanaged: sig(riskCostOpen, 5) ?? 0, riskLevelUnmanaged: RISK_NAMES[topLevelOpen], inspectionPlan: { interval: sig(mg.interval, 4) ?? 0, riskBased: mg.auto, repairThreshold: sig(mg.threshold, 4) ?? 0, expectedRepairs: sig(mg.repairs, 4) ?? 0 },
    buckleArrestors: { required: arr.required, provided: arr.provided, spacing: arr.provided ? sig(arr.spacing, 5) : null, count: arr.count, firstDepth: arr.required ? sig(arr.depthLimit, 5) : null, crossoverUtil: sig(arr.util, 5) ?? 0, confinedLength: arr.required ? sig(arr.provided ? arr.confined : arr.unconfined, 5) : 0 },
    lateralBuckling: { needed: lat.needed, managed: lat.managed, spacing: lat.needed && lat.managed ? sig(lat.spacing, 5) : null, count: lat.count, strain: sig(lat.managed ? lat.strain : lat.rogueStrain, 5) ?? 0, strainAllowable: sig(lat.allow, 5) ?? 0, strainUtil: sig(lat.util, 5) ?? 0, feedIn: sig(lat.plan.governing?.feedIn ?? 0, 5) ?? 0, forceRatio: sig(bk.hobbsUtil, 5) ?? 0 },
    feStressConcentration: sig(cont.gr.scf, 5) ?? 0, feLigamentFactor: sig(cont.gr.ligamentFactor, 5) ?? 0, feCollapsePressure: cont.gr.collapse === null ? null : sig(cont.gr.collapse / BAR, 5), feLameRatio: sig(cont.cyl.hoopBore / cont.la.hoop, 6) ?? 0, shellLigamentRatio: sig(cont.sh.hoopMembrane[0] / cont.gr.hoopLigament, 5) ?? 0, ringPeakHoop: sig(cont.hoopPeak / MPA, 5) ?? 0, ringUtil: sig(cont.ringUtil, 5) ?? 0,
    trackedErosionRate: tre ? sig(tre.rateMax, 5) ?? 0 : 0, trackedErosionAngle: tre ? sig(tre.tr.angleMax, 4) ?? 0 : null, trackedToDnv: tre && tre.dnv > 0 ? sig(tre.rateMax / tre.dnv, 4) : null, criticalFlowVelocity: span.crit.vc === null ? null : sig(span.crit.vc, 5), flowFrequencyRatio: sig(span.fFlow / span.f1, 6) ?? 0, twoWayStressRatio: span.cmp && span.cmp.one > 0 ? sig(span.cmp.two / span.cmp.one, 5) : null,
    vibrationVelocity: sig(span.velRms, 5) ?? 0, vibrationAllowable: sig(span.vAllow, 5) ?? 0, vibrationUtil: sig(span.vibUtil, 5) ?? 0, equipmentUtil: eqWorst ? sig(eqWorst.util, 5) ?? 0 : 0, equipment: eqp.map((e) => ({ name: e.name, x: sig(e.xNode, 6) ?? 0, force: sig(e.F / 1000, 5) ?? 0, moment: sig(e.M / 1000, 5) ?? 0, utilisation: sig(e.util, 5) ?? 0 })), genericFailureFrequency: sig(genericPof, 5) ?? 0,
    mcPof: sig(rel.is.pf, 5) ?? 0, mcBeta: sig(rel.is.beta ?? rel.f.beta, 5) ?? 0, formBeta: sig(rel.f.beta, 6) ?? 0, plugProjectileEnergy: sig(plugE / 1e6, 5) ?? 0, flowSource: st.source,
  };
  prog(1, 'Done');
  const worst = crit[0];
  return {
    summary: `${codeName} pressure containment is used to ${fmt(100 * sa.hoopUtil, 3)} %, collapse to ${fmt(100 * sa.collapseUtil, 3)} %; wall loss of ${fmt(rate, 3)} mm/y leaves ${fmt(remLife, 3)} y on the allowance, the span weld has a calculated fatigue life of ${fmt(fatCalc, 3)} y (${fmt(fatLife, 3)} y with the design fatigue factor), and with inspection every ${fmt(mg.interval, 2)} y the annual failure probability averages ${fmt(mg.avg, 2)} over the remaining life (${mg.met ? 'target met' : 'target exceeded'}; risk ${RISK_NAMES[topLevel]}, ${fmt(riskCost / 1e6, 3)} M$/y); the most utilised item is “${worst.mechanism}” at ${fmt(worst.x / 1000, 3)} km.${task ? task.text : ''}`,
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
  // further test types (one measurement per row), synthesised with modulus 200 GPa, yield 495 MPa, n = 18, S–N shift 0.35, SCF 1.45, C = 2.6e-11, m = 2.9, support stiffness 6000 kN·m/rad,
  // friction 0.72 / 0.41, roughness 80 µm, slug-force coefficient 1.15, defect growth × 1.8, 0.42 failures per 1000 km·y, yield bias 1.10 with 4 % scatter, repair cost 22 M$
  {cStrain: 0.1, cStress: 199.2},
  {cStrain: 0.18, cStress: 354.9},
  {cStrain: 0.22, cStress: 420},
  {cStrain: 0.26, cStress: 451.7},
  {cStrain: 0.4, cStress: 492.1},
  {cStrain: 0.8, cStress: 521.6},
  {cStrain: 1.5, cStress: 551.6},
  {cSrange: 60, logN: 5.806},
  {cSrange: 80, logN: 5.527},
  {cSrange: 100, logN: 5.268},
  {cSrange: 140, logN: 4.718},
  {cSrange: 200, logN: 4.395},
  {cDK: 8, logGrowth: -7.955},
  {cDK: 12, logGrowth: -7.461},
  {cDK: 20, logGrowth: -6.824},
  {cDK: 30, logGrowth: -6.261},
  {cDK: 45, logGrowth: -5.754},
  {cSpanL: 10, natFreq: 5.432},
  {cSpanL: 15, natFreq: 2.598},
  {cSpanL: 20, natFreq: 1.524},
  {cSpanL: 30, natFreq: 0.7327},
  {cSoilW: 300, soilLat: 207.9, soilAx: 124.2},
  {cSoilW: 500, soilLat: 368.9, soilAx: 202.6},
  {cSoilW: 800, soilLat: 587.8, soilAx: 299},
  {cTauV: 1, cTauRho: 800, cTauMu: 0.001, wallShear: 1.79},
  {cTauV: 2, cTauRho: 800, cTauMu: 0.001, wallShear: 6.394},
  {cTauV: 4, cTauRho: 850, cTauMu: 0.002, wallShear: 28.4},
  {cTauV: 6, cTauRho: 1000, cTauMu: 0.001, wallShear: 66.87},
  {cSlugRho: 600, cSlugV: 4, bendForce: 1.675},
  {cSlugRho: 750, cSlugV: 6, bendForce: 4.622},
  {cSlugRho: 850, cSlugV: 9, bendForce: 11.6},
  {cNomStress: 40, hotStress: 60.56},
  {cNomStress: 80, hotStress: 116.7},
  {cNomStress: 120, hotStress: 170.6},
  {cT: 60, cP: 90, cPco2: 1.2, cVliq: 1.5, cInh: 1, cDefD0: 1, cDefYears: 4, defDepth: 2.595},
  {cT: 60, cP: 90, cPco2: 1.2, cVliq: 1.5, cInh: 1, cDefD0: 2, cDefYears: 5, defDepth: 4.171},
  {cT: 60, cP: 90, cPco2: 1.2, cVliq: 1.5, cInh: 1, cDefD0: 1.5, cDefYears: 8, defDepth: 4.308},
  {cExposure: 20000, failures: 8.541},
  {cExposure: 60000, failures: 22.62},
  {cExposure: 150000, failures: 63.85},
  {cQuantile: 0.1, yieldAtP: 520.9},
  {cQuantile: 0.3, yieldAtP: 536.5},
  {cQuantile: 0.5, yieldAtP: 549},
  {cQuantile: 0.7, yieldAtP: 553.1},
  {cQuantile: 0.9, yieldAtP: 572.7},
  {cDownDays: 20, cProd: 8000, incidentCost: 35.65},
  {cDownDays: 45, cProd: 12000, incidentCost: 63.34},
  {cDownDays: 90, cProd: 10000, incidentCost: 88.28},
];
const CAL_VALIDATION = [
  {cT: 38, cP: 62, cPco2: 1, cVliq: 1.4, cInh: 1, corrRate: 0.1586, cVmix: 13, cRho: 130, cSand: 45, eroRate: 0.03958, cFreq: 1.6, cForce: 350, vibAmp: 3.858},
  {cT: 58, cP: 88, cPco2: 2.2, cVliq: 2.8, cInh: 0, corrRate: 7.563, cVmix: 17, cRho: 105, cSand: 70, eroRate: 0.1478, cFreq: 1.8, cForce: 280, vibAmp: 10.93},
  {cT: 68, cP: 98, cPco2: 2.6, cVliq: 1.9, cInh: 1, corrRate: 0.514, cVmix: 15, cRho: 125, cSand: 30, eroRate: 0.04012, cFreq: 1.95, cForce: 320, vibAmp: 5.6},
  {cT: 48, cP: 78, cPco2: 1.4, cVliq: 2.3, cInh: 0, corrRate: 4.496, cVmix: 19, cRho: 95, cSand: 95, eroRate: 0.2994, cFreq: 2.1, cForce: 380, vibAmp: 3.093},
  {cT: 72, cP: 105, cPco2: 3.2, cVliq: 3.1, cInh: 1, corrRate: 0.9601, cVmix: 21, cRho: 88, cSand: 110, eroRate: 0.5132, cFreq: 1.4, cForce: 420, vibAmp: 2.64},
  {cT: 62, cP: 92, cPco2: 1.6, cVliq: 1.1, cInh: 1, corrRate: 0.2207, cVmix: 10, cRho: 170, cSand: 18, eroRate: 0.006022, cFreq: 2.35, cForce: 450, vibAmp: 1.803},
  {cStrain: 0.15, cStress: 305.1},
  {cStrain: 0.3, cStress: 464.4},
  {cStrain: 1, cStress: 537.3},
  {cSrange: 70, logN: 5.587},
  {cSrange: 120, logN: 5.107},
  {cDK: 10, logGrowth: -7.642},
  {cDK: 25, logGrowth: -6.586},
  {cSpanL: 12, natFreq: 3.889},
  {cSpanL: 25, natFreq: 1.036},
  {cSoilW: 400, soilLat: 271.5, soilAx: 162.2},
  {cTauV: 3, cTauRho: 800, cTauMu: 0.001, wallShear: 14.92},
  {cSlugRho: 700, cSlugV: 5, bendForce: 2.744},
  {cNomStress: 60, hotStress: 85.1},
  {cExposure: 80000, failures: 31.59},
  {cQuantile: 0.2, yieldAtP: 527.5},
  {cQuantile: 0.8, yieldAtP: 563.7},
  {cDownDays: 30, cProd: 9000, incidentCost: 42.52},
];

// ---- calibration model: closed-form predictions for one measured operating point ------------------------------------------
function calModel(v0 = {}) {
  // only the test types whose condition columns are present in the row are evaluated (all of them when none is given), so that a fit over many rows stays fast
  const v = (k) => (v0[k] !== undefined && v0[k] !== null && v0[k] !== '' ? v0[k] : DEFAULTS[k]), has = (...ks) => ks.some((k) => isNum(v0[k])), all = !has('cT', 'cPco2', 'cVliq', 'cVmix', 'cSand', 'cFreq', 'cForce', 'cStrain', 'cSrange', 'cDK', 'cSpanL', 'cSoilW', 'cTauV', 'cSlugV', 'cNomStress', 'cDefYears', 'cExposure', 'cQuantile', 'cDownDays'), out = {};
  const ID = pos(v('idMm'), BASE.idMm) / 1000, t = pos(v('wtMm'), BASE.wtMm) / 1000, D = ID + 2 * t, Dh = Math.max(nz(v('odCoatMm'), 0) / 1000, D), Ai = (Math.PI * ID * ID) / 4, E = pos(v('eMod'), BASE.E), sy = pos(v('smys'), BASE.smys), scf = pos(v('scf'), 1);
  const corrOf = () => deWaardMilliams({ T: nz(v('cT'), 60), pCO2: Math.max(nz(v('cPco2'), 1), 0), P: nz(v('cP'), 80), model: v('corrModel') === '1991' ? '1991' : '1995', U: Math.max(nz(v('cVliq'), 2), 0.01), d: ID, pH: nz(v('phAct'), 0), glycolWt: nz(v('glycolWt'), 0), inhibEff: (clamp(nz(v('inhibEff'), 0), 0, 100) / 100) * clamp(nz(v('cInh'), 1), 0, 1), mult: nz(v('corrMult'), 1) }).rate;
  if (all || has('cT', 'cPco2', 'cVliq')) out.corrRate = corrOf();
  if (all || has('cVmix', 'cSand')) out.eroRate = erosionRate({ model: 'dnv', mp: Math.max(nz(v('cSand'), 10), 0) / 86400, U: Math.max(nz(v('cVmix'), 10) - nz(v('eroVth'), 0), 0), D: ID, dp: pos(v('sandUm'), 250) * 1e-6, rhoM: pos(v('cRho'), 150), muM: 1e-4, rhoP: pos(v('sandDensity'), 2650), rOverD: 1.5, gf: nz(v('geomFactor'), 1), mult: nz(v('erosionMult'), 1) }).rate;
  const L = pos(v('spanLength'), 15), EI = E * MPA * (Math.PI / 64) * (D ** 4 - ID ** 4), me = pos(v('rhoSteel'), 7850) * (Math.PI / 4) * (D * D - ID * ID) + nz(v('coatDensity'), 900) * (Math.PI / 4) * (Dh * Dh - D * D) + 500 * Ai + nz(v('addedMass'), 1) * RHO_SW * (Math.PI / 4) * Dh * Dh;
  if (all || has('cFreq', 'cForce')) { // resonance curve of the first span mode (pinned–pinned) under a harmonic mid-span force
    const f1 = (Math.PI / (2 * L * L)) * Math.sqrt(EI / me), r = Math.max(nz(v('cFreq'), 0.3), 0) / f1, z = clamp(nz(v('damping'), 2), 0.01, 100) / 100; out.vibAmp = ((nz(v('cForce'), 500) * L ** 3) / (48 * EI) / Math.sqrt((1 - r * r) ** 2 + (2 * z * r) ** 2)) * 1000; }
  if (all || has('cStrain')) out.cStress = rambergStress(Math.max(nz(v('cStrain'), 0.2), 0) / 100, E, sy, clamp(nz(v('hardenN'), 15), 2, 80));
  if (all || has('cSrange')) out.logN = Math.log10(Math.min(snCycles(pos(v('cSrange'), 100), SN_CURVES[v('snClass')] ? v('snClass') : 'F1', ['air', 'cp', 'free'].includes(v('snEnv')) ? v('snEnv') : 'cp', { scf, t: t * 1000, shift: nz(v('snShift'), 0) }), 1e30));
  if (all || has('cDK')) out.logGrowth = Math.log10(pos(v('parisC'), 1.65e-11)) + pos(v('parisM'), 3) * Math.log10(pos(v('cDK'), 20));
  if (all || has('cSpanL')) { // first frequency of a span whose ends are held by connections of finite rotational stiffness
    const Lf = pos(v('cSpanL'), L), kr = Math.max(nz(v('spanKr'), 0), 0) * 1000, key = `${Lf}|${EI}|${me}|${kr}`; let fK = CAL_F1.get(key);
    if (fK === undefined) { fK = beamModes(beamModel({ L: Lf, EI, m: me, n: 6, ends: ['flange', 'flange'], kR: [kr, kr] }), 1).f[0]; if (CAL_F1.size > 400) CAL_F1.clear(); CAL_F1.set(key, fK); }
    out.natFreq = fK; }
  if (all || has('cSoilW')) { const wS = pos(v('cSoilW'), 500); out.soilLat = nz(v('muLateral'), 0.6) * wS; out.soilAx = nz(v('muAxial'), 0.5) * wS; }
  if (all || has('cTauV')) { const vT = Math.max(nz(v('cTauV'), 2), 1e-6), rT = pos(v('cTauRho'), 800); out.wallShear = (frictionFactor(Math.max((rT * vT * ID) / pos(v('cTauMu'), 1e-3), 10), (pos(v('roughUm'), 45) * 1e-6) / ID) / 8) * rT * vT * vT; }
  if (all || has('cSlugV')) out.bendForce = (nz(v('slugCf'), 1) * nz(v('dlf'), 2) * pos(v('cSlugRho'), 700) * Ai * nz(v('cSlugV'), 5) ** 2 * Math.SQRT2) / 1000;
  if (all || has('cNomStress')) out.hotStress = scf * nz(v('cNomStress'), 50);
  if (all || has('cDefYears')) out.defDepth = nz(v('cDefD0'), 1) + (out.corrRate ?? corrOf()) * nz(v('defectGrowth'), 1) * nz(v('cDefYears'), 5);
  if (all || has('cExposure')) out.failures = (nz(v('failRate'), 0.5) * nz(v('cExposure'), 10000)) / 1000;
  if (all || has('cQuantile')) out.yieldAtP = randomVariable({ dist: 'lognormal', mean: pos(v('yieldBias'), 1.08) * sy, cov: clamp(nz(v('covYield'), 5), 0.1, 50) / 100 }).x(PhiInv(clamp(nz(v('cQuantile'), 0.5), 1e-6, 1 - 1e-6)));
  if (all || has('cDownDays')) out.incidentCost = nz(v('repairCost'), 15) + (nz(v('cDownDays'), 30) * nz(v('cProd'), 10000) * nz(v('oilPrice'), 75)) / 1e6;
  return out;
}
const CAL_F1 = new Map();

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
  // ---- continuum and shell finite elements
  { const ri = 0.1, ro = 0.15, pi = 50e6, pe = 10e6, Ee = 207e9, nu = 0.3, A = (pi * ri * ri - pe * ro * ro) / (ro * ro - ri * ri), Bc = ((pi - pe) * ri * ri * ro * ro) / (ro * ro - ri * ri), uL = (r) => ((1 + nu) / Ee) * ((1 - 2 * nu) * A * r + Bc / r), W = Math.PI * (pi * ri * uL(ri) - pe * ro * uL(ro));
    const c2 = feCylinder({ ri, ro, pi, pe, E: Ee, nu, nr: 4, order: 2 }), err = (order, nr) => Math.abs(feCylinder({ ri, ro, pi, pe, E: Ee, nu, nr, order }).energy - W) / W;
    chk('Axisymmetric solid elements: hoop stress at the bore (Lamé)', c2.hoopBore / MPA, (A + Bc / (ri * ri)) / MPA, 0.005 * (A + Bc / (ri * ri)) / MPA, 'A + B/rᵢ², 4 quadratic elements through the wall');
    chk('Axisymmetric solid elements: radial displacement of the bore (Lamé, plane strain)', c2.uBore * 1e6, uL(ri) * 1e6, 1e-3 * uL(ri) * 1e6, '(1 + ν)/E·[(1 − 2ν)·A·r + B/r] in µm');
    chk('Finite-element equilibrium: ∫σθ dr = pᵢ·rᵢ − pₑ·rₒ', c2.hoopForce, pi * ri - pe * ro, 1e-8 * (pi * ri - pe * ro), 'Cauchy equilibrium of the half pipe, satisfied exactly by the discrete solution');
    chk('Element-order convergence: linear triangles (strain energy)', Math.log2(err(1, 4) / err(1, 8)), 2, 0.2, 'error ∝ h² for 3-node elements');
    chk('Element-order convergence: quadratic triangles (strain energy)', Math.log2(err(2, 2) / err(2, 4)), 4, 0.3, 'error ∝ h⁴ for 6-node elements'); }
  { const m = mapMesh(2, 2, (a, b) => [a, b], 1), fix = [{ node: m.id(0, 0), dir: 0 }, { node: m.id(0, 1), dir: 0 }, { node: m.id(0, 2), dir: 0 }, { node: m.id(0, 0), dir: 1 }], fe = feSolve({ mesh: m, E: 200e9, nu: 0.25, mode: 'stress', fix, pressures: meshEdges(m, 'a1').map((nodes) => ({ nodes, p: -100e6, inside: [0.5, 0.5] })) });
    chk('Plane-stress patch test: uniform tension', fe.u[2 * m.id(2, 1)] * 1e6, (100e6 / 200e9) * 1e6, 1e-6, 'u = σ·L/E on a unit square, in µm'); chk('Plane-stress patch test: Poisson contraction', fe.u[2 * m.id(2, 2) + 1] * 1e6, -0.25 * (100e6 / 200e9) * 1e6, 1e-6, 'v = −ν·σ·L/E, in µm'); }
  { const ro = 0.2, t = 0.01, f0 = 0.01, pe = 2e6, w1 = (f0 * (2 * ro - t)) / 4, r = feRing({ ro, t, pi: 0, pe, E: 207e9, nu: 0.3, ovality: f0, nth: 40, nr: 2 });
    chk('Ring with out-of-roundness: bending stress at the crown', (r.crown.outer - r.crown.inner) / 2 / MPA, (6 * pe * ro * w1) / (t * t) / MPA, 0.03 * ((6 * pe * ro * w1) / (t * t) / MPA), 'M = p·r·w₁ (first order), σ = 6M/t²');
    chk('Ring finite elements: pressure resultant on the half ring', r.balance.reaction, -2 * pe * (ro + w1), 1e-3 * 2 * pe * ro, 'external pressure × projected width 2·(rₒ + w₁), carried by the two cuts'); }
  { const ro = 0.15, t = 0.015, d = 0.006, ri = ro - t, g = feGroove({ ro, t, d, L: 2.5, side: 'in', pi: 10e6, pe: 0, E: 207e9, nu: 0.3, nz: 24, nr: 2 });
    chk('Thinned wall, long band: ligament hoop stress from equilibrium', g.ligamentFactor, ((ri + d) / ri) * (t / (t - d)), 0.01, 'p·(rᵢ + d)/(t − d) against p·rᵢ/t'); }
  { const R = 0.5, t = 0.01, p = 1e6, Ee = 207e9, nu = 0.3, beta = ((3 * (1 - nu * nu)) / (R * R * t * t)) ** 0.25, wm = ((p * R * R) / (Ee * t)) * (1 - nu / 2), xq = Math.PI / (2 * beta), s = shellCylinder({ R, L: 1, t, E: Ee, nu, p, n: 120, axial: 'capped', ends: ['clamped', 'symmetry'], breaks: [xq] }), iq = s.x.findIndex((x) => Math.abs(x - xq) < 1e-9);
    chk('Cylindrical shell: membrane deflection far from the edge', s.w[s.w.length - 1] * 1e6, wm * 1e6, 1e-4 * wm * 1e6, 'p·R²/(E·t)·(1 − ν/2), capped ends, in µm');
    chk('Cylindrical shell: edge bending moment at a clamped end', Math.abs(s.moment[0]), (p * (1 - nu / 2)) / (2 * beta * beta), 0.01 * ((p * (1 - nu / 2)) / (2 * beta * beta)), 'M₀ = p·(1 − ν/2)/(2β²), β⁴ = 3(1 − ν²)/(R²t²)');
    chk('Cylindrical shell: decay of the edge disturbance', s.w[iq] / wm, 1 - Math.exp(-Math.PI / 2), 2e-3, 'w/w_m = 1 − e^(−βx)(cos βx + sin βx) at βx = π/2');
    chk('Cylindrical shell: radial equilibrium', s.radialResistance, s.radialLoad, 1e-8 * s.radialLoad, 'pressure load against hoop resistance and edge shear'); }
  // ---- pipe conveying fluid (two-way coupling)
  { const mf = 100, md = beamModel({ L, EI, m, n: 24 }), mo = beamModes(md, 6), mi = modalIntegrals(md, mo), vc = (Math.PI / L) * Math.sqrt(EI / mf), v = 0.6 * vc;
    chk('Pipe conveying fluid: critical velocity, pinned–pinned', fluidCritical(md, mo, mf, mi).vc, vc, 1e-5 * vc, '(π/L)·√(EI/m_f)');
    const k1 = ((Math.PI / L) ** 4 * EI) / m - (mf * v * v * (Math.PI / L) ** 2) / m, k2 = (16 * (Math.PI / L) ** 4 * EI) / m - (4 * mf * v * v * (Math.PI / L) ** 2) / m, g = ((2 * mf * v) / m) * (8 / 3) / L, bq = k1 + k2 + g * g, f2 = Math.sqrt((bq - Math.sqrt(bq * bq - 4 * k1 * k2)) / 2) / (2 * Math.PI), fF = fluidFrequency(md, mo, mf, v, mi);
    chk('Pipe conveying fluid: frequency at 60 % of the critical velocity', fF, f2, 2e-3 * f2, 'two-mode Galerkin solution with the Coriolis coupling 8/3');
    chk('Pipe conveying fluid: frequency reduction with flow', fF / mo.f[0], Math.sqrt(1 - 0.36), 0.02, 'close to √(1 − (U/U_c)²); the Coriolis coupling lowers it slightly more');
    const dt = 1 / (fF * 400), r = fsiResponse({ model: md, modes: mo, zeta: 0, slugs: [], v, mf, dM: 0, dt, tEnd: 20 / fF, ro: 0.15, E: 207e9, q0: [1e-3, 0, 0, 0, 0, 0], mi }), zc = []; for (let i = 1; i < r.disp.length; i++) if (r.disp[i - 1] < 0 && r.disp[i] >= 0) zc.push(r.t[i - 1] + (dt * -r.disp[i - 1]) / (r.disp[i] - r.disp[i - 1]));
    chk('Two-way coupled time integration reproduces the gyroscopic frequency', (zc.length - 1) / (zc[zc.length - 1] - zc[0]), fF, 2e-3 * fF, 'free vibration with flow, zero crossings over about 20 periods');
    const slugs = [{ t0: 0.2, len: 8, v: 5 }], a = slugResponse({ model: md, modes: mo, zeta: 0.02, slugs, dw: 300, dt: 0.004, tEnd: 10, ro: 0.15, E: 207e9 }), b = fsiResponse({ model: md, modes: mo, zeta: 0.02, slugs, v: 5, mf, dM: 30, dw: 300, dt: 0.004, tEnd: 10, ro: 0.15, E: 207e9, coupled: false, mi });
    chk('One-way limit of the coupled solver', b.sigmaMax / MPA, a.sigmaMax / MPA, 1e-8 * (a.sigmaMax / MPA), 'without feedback the coupled solver equals the modal slug response'); }
  // ---- supports, guides, connections
  chk('Fixed–guided beam: first natural frequency', beamModes(beamModel({ L, EI, m, n: 32, ends: ['fixed', 'guided'] }), 1).f[0], (2.36502 ** 2 / (2 * Math.PI)) * w0, 1e-4 * w0, 'βL = 2.3650 (tan βL = −tanh βL)');
  chk('Span anchored at mid-length', beamModes(beamModel({ L, EI, m, n: 32, supports: [{ x: L / 2, type: 'anchor' }] }), 1).f[0], (15.4182 / (2 * Math.PI)) * Math.sqrt(EI / m) / (L / 2) ** 2, 1e-3 * w0 * 4, 'each half is a fixed–pinned beam: (βL)² = 15.418');
  { const Mc = 2 * m * L, got = beamModes(beamModel({ L, EI, m, n: 32, supports: [{ x: L / 2, type: 'mass', mass: Mc }] }), 1).f[0], ex = Math.sqrt((48 * EI) / (L ** 3 * (Mc + (17 / 35) * m * L))) / (2 * Math.PI); chk('Lumped mass (flange or valve) at mid-span', got, ex, 0.01 * ex, 'Rayleigh: √(48EI/(L³(M + 17/35·m·L)))/2π'); }
  { const kr = 1e12, got = beamModes(beamModel({ L, EI, m, n: 32, ends: ['flange', 'flange'], kR: [kr, kr] }), 1).f[0]; chk('Connection of very high rotational stiffness tends to the fixed end', got, 3.56097 * w0, 2e-3 * w0, 'limit of the flange-type end condition'); }
  // ---- tracked particles
  { const D = 0.2, U = 10, R = 1.5 * D, F = bendFlowField({ D, R, angle: 170, U, ns: 96, nn: 16 }), i = Math.round(F.s.length / 2), C = (U * D) / R / Math.log((R + D / 2) / (R - D / 2));
    chk('Bend flow field: free-vortex velocity inside the bend', F.us[i][4], C / (1 + F.n[4] / R), 2e-3 * U, 'u·r = constant with the same flow rate');
    const Fb = bendFlowField({ D, R, angle: 90, U }), dp = 250e-6, tb = bendErosionTracked({ D, rOverD: 1.5, angle: 90, U, rho: 1e-6, mu: 1e-9, dp, mp: 1, model: 'dnv', field: Fb, nPart: 400, en: 0, et: 0 }); let ex = 0; const half = D / 2 - dp / 2; for (let k = 0; k < 4000; k++) ex += impactWear('dnv', U, Math.acos((R + (-half + ((k + 0.5) / 4000) * 2 * half)) / (R + half))) / 4000;
    chk('Tracked particles, ballistic limit: wear per kg of sand', tb.wearPerKg * 1e12, ex * 1e12, 0.005 * ex * 1e12, 'straight paths hit the outer wall at α = acos((R + n)/R_o); mean of the impact law, in mm³/kg × 10⁻³');
    chk('Tracked particles, ballistic limit: every particle hits the outer wall', tb.hitFraction, 1, 1e-12, 'no fluid drag');
    chk('Tracked particles, tracer limit: no wear', bendErosionTracked({ D, rOverD: 1.5, angle: 90, U, rho: 1000, mu: 1, dp: 5e-6, mp: 1, model: 'dnv', field: Fb, nPart: 40 }).wearPerKg / ex, 0, 1e-4, 'particles of vanishing Stokes number follow the streamlines (wear relative to the ballistic limit)'); }
  // ---- design measures
  { const o = { EI: 5e7, EA: 2.5e9, w: 800, muA: 0.5, muL: 0.6, P0: 3e6, ro: 0.15 }, lb = lateralBuckle({ ...o, spacing: Infinity }).modes.find((q) => q.name === 'Lateral mode 3'), Lr = lb.L, hand = (34.06 * o.EI) / (Lr * Lr) + 1.294 * o.muA * o.w * Lr * (Math.sqrt(1 + (1.668e-4 * o.EA * o.muL ** 2 * o.w * Lr ** 5) / (o.muA * o.EI ** 2)) - 1);
    chk('Planned lateral buckle: unlimited feed-in reproduces Hobbs\' closed form', hand / 1e6, o.P0 / 1e6, 1e-6, 'P₀ = k₁EI/L² + k₃μ_a·w·L·[√(1 + k₂·EA·μ_l²·w·L⁵/(μ_a·EI²)) − 1], mode 3');
    chk('Planned lateral buckle: bending moment k₅·μ_l·w·L²', lb.moment, 0.1434 * 0.6 * 800 * Lr * Lr, 1e-6 * lb.moment, 'Hobbs mode 3');
    const near = lateralBuckle({ ...o, spacing: 1500 }).modes.find((q) => q.name === 'Lateral mode 3'); chk('Planned lateral buckle: closer initiators reduce the strain', near.strain < lb.strain && near.limited ? 1 : 0, 1, 0, 'feed-in limited to half the spacing on each side'); }
  chk('Strain capacity of the displacement-controlled criterion', strainCapacity({ D: 0.3, t: 0.015, fy: 450e6, fu: 535e6, dp: 0, alphaH: 0.93, alphaGw: 1 }), 0.78 * 0.04 * 0.93 ** -1.5, 1e-12, '0.78·(t/D − 0.01)·α_h^−1.5 without internal overpressure');
  { const x = arrestorCrossover({ D: 0.3, t: 0.015, t2: 0.03, L: 1e3, fy: 450e6 }); chk('Buckle arrestor: a long arrestor reaches its own propagation pressure', x.pX / MPA, (35 * 450 * 0.1 ** 2.5), 1e-9, '35·f_y·(t₂/D)^2.5'); chk('Buckle arrestor: crossover pressure of a short arrestor', arrestorCrossover({ D: 0.3, t: 0.015, t2: 0.03, L: 0.3, fy: 450e6 }).pX / MPA, 35 * 450 * 0.05 ** 2.5 + (35 * 450 * 0.1 ** 2.5 - 35 * 450 * 0.05 ** 2.5) * (1 - Math.exp(-2)), 1e-9, 'p_pr + (p_pr,BA − p_pr)·(1 − e^(−20·t₂·L/D²)), exponent 2'); }
  // ---- inspection planning, distribution fitting, scoring
  { const cov = 0.5, z = Math.sqrt(Math.log(1 + cov * cov)), lam = -z * z / 2, tail = (x) => 1 - Phi((Math.log(x) - lam) / z), ls = [{ name: 'step', pf: (l) => (l >= 8 ? 1 : 0), annual: false, rate: 1, loss0: 0, cov }];
    const a = inspectionPlan({ ls, years: 12, interval: 0, pod: () => 0, sizingSd: 0.1, threshold: 3, nq: 400 }); chk('Inspection plan without inspections: cumulative failure probability', sum(a.annual), tail(8 / 12), 0.004, 'P(X·r·T ≥ loss limit), lognormal X');
    const b = inspectionPlan({ ls, years: 12, interval: 5, pod: () => 1, sizingSd: 1e-6, threshold: 3, nq: 400 }); chk('Inspection plan with a perfect inspection and repair', sum(b.annual), 2 * tail(8 / 5), 0.006, 'inspections at years 5 and 10 renew every section above the threshold: a failure needs 8 mm of loss within one 5-year interval, in the first or in the second');
    chk('Bayesian updating: the posterior mean of the wall-loss factor falls when no repair was needed', b.inspections[0].mean < 1 && b.inspections[0].cov < cov ? 1 : 0, 1, 0, 'surviving sections have X < threshold/(r·t)'); }
  { const q = Array.from({ length: 400 }, (_, i) => (i + 0.5) / 400), ln = fitDistribution(q.map((u) => Math.exp(1 + 0.3 * PhiInv(u))), 'lognormal'), wb = fitDistribution(q.map((u) => 5 * (-Math.log(1 - u)) ** (1 / 2.5)), 'weibull'), gb = fitDistribution(q.map((u) => 10 - 2 * Math.log(-Math.log(u))), 'gumbel');
    chk('Maximum likelihood: lognormal parameters', ln.params.sigma, 0.3, 0.004, 'sample of exact quantiles, σ of ln X'); chk('Maximum likelihood: Weibull shape', wb.params.shape, 2.5, 0.05, 'sample of exact quantiles'); chk('Maximum likelihood: Gumbel scale', gb.params.scale, 2, 0.05, 'sample of exact quantiles'); }
  chk('Brier score of two forecasts', brierScore([0.1, 0.9], [0, 1]).brier, 0.01, 1e-12, 'mean of (p − o)²');
  // ---- sourced constants
  chk('NORSOK M-506 at 20 °C, 1 bar, 19 Pa, pH 4', norsokM506({ T: 20, pCO2: 1, pH: 4, S: 19, fugacity: false }).rate, 4.762 * (2.0676 - 0.2309 * 4), 1e-9, 'K_t·f(pH), hand calculation from Tables 1 and 2');
  chk('NORSOK M-506: interpolation between 40 and 60 °C', norsokM506({ T: 50, pCO2: 1, pH: 4, S: 19, fugacity: false }).rate, 0.5 * (8.927 + 10.695 * (1.836 - 0.1818 * 4) / (2.0676 - 0.2309 * 4)) * (2.0676 - 0.2309 * 4), 1e-9, 'mean of the rates at 40 and 60 °C');
  chk('NORSOK M-506 wall shear stress', norsokShear(1000, 1e-3, 2, 0.1, 50e-6), 0.5 * 1000 * 0.001375 * (1 + (10 + 5) ** 0.33) * 4, 1e-9, 'f = 0.001375·[1 + (20000k/D + 10⁶μ/(ρuD))^0.33]');
  chk('pH of CO2-saturated water at 25 °C and 1 bar', norsokPH({ T: 25, pCO2: 1, P: 1 }).pH, 3.71 + 0.00417 * 25, 0.12, 'against the de Waard expression');
  chk('S–N curve B1 with cathodic protection at 200 MPa', snCycles(200, 'B1', 'cp'), 10 ** (14.917 - 4 * Math.log10(200)), 1, 'log a₁ = 14.917 (DNV-RP-C203 Table 2-2)');
  chk('S–N curve F3, free corrosion, at 100 MPa', snCycles(100, 'F3', 'free'), 10 ** (11.068 - 6), 1e-6, 'log a = 11.068, m = 3 (DNV-RP-C203 Table 2-4)');
  chk('Sour-service region at 0.5 kPa H2S and pH 4.3', sourRegion(0.5, 4.3).region, 1, 0, 'above the line pH = 4.5 + log p of ISO 15156-2 Figure 1');
  chk('DNV-RP-F101 best-estimate capacity factor', b31g({ D: 0.5, t: 0.01, d: 0.004, L: 0.2, smys: 450e6, smts: 535e6, method: 'dnvCap' }).pf / b31g({ D: 0.5, t: 0.01, d: 0.004, L: 0.2, smys: 450e6, smts: 535e6, method: 'dnv' }).pf, 1.05, 1e-12, 'capacity equation of section 2.3');
  chk('Ramberg–Osgood strain at the yield strength', rambergStrain(450, 207000, 450, 15), 450 / 207000 + 0.002, 1e-12, 'elastic strain plus 0.2 % offset');
  chk('Girth-weld misalignment factor', params({ hiLoMm: 1.59 }).scfWeld, 1 + ((3 * 1.59) / 15.9) * Math.exp(-Math.sqrt(15.9 / 285.8)), 1e-9, '1 + (3δ/t)·exp(−√(t/D)), δ = 0.1 t on the reference pipe');
  verifyCoupled(chk);
  return out;
}
/** Verification of the coupled fluid–structure solvers, the solid elements, the corrosion extensions and the frequency and model-uncertainty statistics. */
const ELBOW_REF = Object.freeze({ flexFactor: 4.2272, hoopFactor: 3.8893 }); // CalculiX 2.23, 2160 twenty-node elements: bend 285.8 × 15.9 mm, R = 381 mm, legs 0.6 m (reference data set calculix-elbow)
function verifyCoupled(chk) {
  const zero = (t, y) => { const z = []; for (let i = 1; i < y.length; i++) if (y[i - 1] > 0 !== y[i] > 0) z.push(t[i - 1] + ((t[i] - t[i - 1]) * y[i - 1]) / (y[i - 1] - y[i])); return z; };
  // ---- liquid-filled elastic tube: standing wave against the dispersion relation with radial fluid inertia, wall inertia and bending (Lamb)
  const soft = { R: 0.005, L: 0.2, h: 0.001, E: 3e5, nu: 0.3, rhoS: 1200, rhoF: 1000, axial: 'restrained' }, k1 = Math.PI / soft.L, an = tubeWaves({ ...soft, k: k1 }), T1 = (2 * Math.PI) / an.omegaLamb;
  { const r = tubeFsi({ ...soft, mu: 0, nx: 32, nr: 5, dt: T1 / 64, tEnd: 2.3 * T1, eta0: (x) => 1e-5 * Math.sin(k1 * x), stations: [soft.L / 2] }), z = zero(r.t, r.eta[0]), Tn = (2 * (z[z.length - 1] - z[0])) / (z.length - 1);
    chk('Coupled tube: natural frequency of the liquid-filled tube', (2 * Math.PI) / Tn, an.omegaLamb, 0.006 * an.omegaLamb, 'standing wave k = π/L, inviscid liquid; dispersion relation ρ·ω²·I₀(kR)/(k·I₁(kR)) + m_w·ω² = k_wall·(1 + D·k⁴/k_s); 32 × 5 cells');
    chk('Coupled tube: energy of the free oscillation is conserved', (r.energy.fluid[r.energy.fluid.length - 1] + r.energy.shell[r.energy.shell.length - 1]) / r.energy.shell[0], 1, 0.02, 'kinetic energy of the liquid + strain and kinetic energy of the shell after 2.3 periods (upwind convection is the only dissipation)');
    chk('Coupled tube: work of the liquid on the wall equals the energy gained by the shell', r.energy.shellGain / r.energy.interface, 1, 1e-4, 'time-centred pressure, trapezoidal rule in both solvers, consistent load and volume-flux operators');
    chk('Coupled tube: the liquid volume follows the wall (incompressible)', (r.volume.fluidChange - r.volume.boundary) / Math.abs(r.volume.boundary), 0, 2e-3, 'change of the meshed volume against the net inflow through the open ends');
    chk('Coupled tube: discrete divergence of the velocity field', r.divergenceMax * T1, 0, 1e-8, 'largest cell divergence × period');
    const bad = tubeFsi({ ...soft, mu: 0, nx: 32, nr: 5, dt: T1 / 64, tEnd: 0.4 * T1, eta0: (x) => 1e-5 * Math.sin(k1 * x), coupling: 'fixed', omega: 1 });
    chk('Coupled tube: the unrelaxed partitioned iteration diverges (added mass)', bad.diverged ? 1 : 0, 1, 0, `added mass of the liquid / mass of the wall = ${fmt(an.addedMass / an.wallMass, 3)}; the interface Newton iteration needs ${fmt(r.iterations.mean, 2)} iterations per step`); }
  { const b = { ...soft, L: 0.3, mu: 0.0035 }, a2 = tubeWaves(b), Tp = 0.01, r = tubeFsi({ ...b, nx: 60, nr: 5, dt: Tp / 40, tEnd: 0.04, inlet: { type: 'pressure', p: (t) => (t < Tp ? 1333 * Math.sin((Math.PI * t) / Tp) ** 2 : 0) }, stations: [0.06, 0.18] }), tp = [0, 1].map((q) => r.t[r.pWall[q].indexOf(Math.max(...r.pWall[q]))]);
    chk('Coupled tube: pulse speed against Moens–Korteweg', 0.12 / (tp[1] - tp[0]), a2.c0, 0.03 * a2.c0, `√(E·h/(2·ρ·R_m·(1 − ν²))) = ${fmt(a2.c0, 4)} m/s (classical thin-wall value ${fmt(a2.cMK, 4)} m/s); peak-to-peak between two stations, 60 × 5 cells`); }
  { const b = { R: 0.127, L: 40, h: 0.0159, E: 207e9, nu: 0.3, rhoS: 7850, rhoF: 1000, mu: 1e-3, K: 2.2e9, axial: 'restrained' }, a3 = tubeWaves(b), Tp = 0.004, r = tubeFsi({ ...b, nx: 80, nr: 3, dt: Tp / 32, tEnd: 0.027, inlet: { type: 'pressure', p: (t) => (t < Tp ? 1e5 * Math.sin((Math.PI * t) / Tp) ** 2 : 0) }, stations: [8, 28] }), tp = [0, 1].map((q) => r.t[r.pWall[q].indexOf(Math.max(...r.pWall[q]))]);
    chk('Coupled tube: water in a steel pipe against the Korteweg wave speed', 20 / (tp[1] - tp[0]), a3.cKorteweg, 0.02 * a3.cKorteweg, `1/c² = ρ/K + 1/c₀², K = 2.2 GPa: ${fmt(a3.cKorteweg, 5)} m/s (rigid pipe ${fmt(Math.sqrt(2.2e9 / 1000), 5)} m/s); 80 × 3 cells`); }
  { const b = { ...soft, mu: 0.01 }, aw = tubeWaves({ ...b, k: k1 }), Tw = (2 * Math.PI) / aw.omegaWomersley[0], r = tubeFsi({ ...b, nx: 24, nr: 14, dt: Tw / 90, tEnd: 3.3 * Tw, eta0: (x) => 1e-5 * Math.sin(k1 * x), stations: [soft.L / 2] }), e = r.eta[0], z = zero(r.t, e), pk = []; for (let i = 1; i < e.length - 1; i++) if (e[i] > e[i - 1] && e[i] >= e[i + 1] && e[i] > 0) pk.push([r.t[i], e[i]]);
    const amp = (om) => { let c = 0, q = 0; for (let i = 0; i < e.length; i++) { c += e[i] * Math.cos(om * r.t[i]); q += e[i] * Math.sin(om * r.t[i]); } return Math.hypot(c, q); }; let wB = 0, aB = 0; for (let j = 0; j <= 600; j++) { const om = aw.omegaLamb * (0.6 + (0.8 * j) / 600), a = amp(om); if (a > aB) { aB = a; wB = om; } }
    const env = (i0, i1) => { let m = 0; for (let i = i0; i < i1; i++) m = Math.max(m, Math.abs(e[i])); return m; }, nP = Math.round(Tw / r.dt), dec = Math.log(env(Math.round(0.25 * nP), Math.round(1.25 * nP)) / env(Math.round(2.25 * nP), Math.round(3.25 * nP))) / (2 * Tw);
    chk('Coupled tube: frequency of the viscous standing wave (Womersley)', wB, aw.omegaWomersley[0], 0.03 * aw.omegaWomersley[0], `ω² = c₀²k²(1 − F₁₀(α)) with wall inertia, α = ${fmt(aw.alpha, 3)}; inviscid ${fmt(aw.omegaLamb, 4)} rad/s`);
    chk('Coupled tube: decay rate of the viscous standing wave (Womersley)', dec, aw.omegaWomersley[1], 0.15 * aw.omegaWomersley[1], 'imaginary part of the complex frequency; the Stokes layer is resolved by three to four of the 14 radial cells; frequency from the spectral peak of the record, decay from the envelope two periods apart'); void z; void pk; }
  { const J = besselJ01([1, 0]), I = besselI01(1); chk('Bessel J₀(1)', J.j0[0], 0.7651976865579666, 1e-13, 'power series'); chk('Bessel J₁(1)', J.j1[0], 0.44005058574493355, 1e-13, 'power series'); chk('Bessel I₀(1) / I₁(1)', I.i0 / I.i1, 1.2660658777520084 / 0.5651591039924851, 1e-12, 'power series'); }
  // ---- resolved flow through a bend
  const D = 0.2, U = 4, rho = 900, A = (Math.PI * D * D) / 4, cfd = { D, R: 0.3, angle: 90, legs: 2, U, rho, mu: 2e-3, ns: 40, nn: 8 };
  { const r = bendCfd({ ...cfd, tEnd: 1.2 }), L = r.last, F = Math.hypot(L.force[0], L.force[1]), cv = Math.hypot(L.cv[0], L.cv[1]), th = Math.hypot((L.pIn + rho * U * U) * A, rho * U * U * A);
    chk('Bend flow: integrated wall tractions against the control-volume momentum theorem', F / cv, 1, 0.03, 'pressure and shear on both walls against (p + ρu²)·A in − out; 40 × 8 cells');
    chk('Bend flow: wall force against ρ·A·v²·√(2(1 − cos θ)) plus the inlet pressure term', F / th, 1, 0.03, `90° bend, uniform inflow, outlet at zero gauge pressure: |((p_in + ρv²)·A, ρv²·A)| = ${fmt(th, 4)} N`);
    chk('Bend flow: discrete divergence of the velocity field', (L.divergenceMax * D) / U, 0, 1e-9, 'largest cell divergence × D/v');
    chk('Bend flow: pressure difference across the bend against ρ·v²·D/R', (L.pExt[Math.floor(L.s.length / 2)] - L.pInt[Math.floor(L.s.length / 2)]) / ((rho * U * U * D) / 0.3), 1, 0.25, 'centrifugal estimate for a uniform velocity; the resolved velocity is faster at the intrados');
    const r2 = bendCfd({ ...cfd, pOut: 5e5, tEnd: 1.2 }), F2 = [r2.last.force[0] - L.force[0], r2.last.force[1] - L.force[1]];
    chk('Bend flow: pressure term of the bend force p·A·√(2(1 − cos θ))', Math.hypot(F2[0], F2[1]) / (5e5 * A * Math.SQRT2), 1, 1e-3, 'difference between the solutions with 5 bar and with zero outlet pressure; the wall normal is taken at the middle of each of the 14 cells of the arc'); }
  { const r = bendCfd({ ...cfd, rho: 150, rhoSlug: 900, slug: { t0: 0.1, len: 8 }, tEnd: 0.1 + 1.6 / U + 0.25 }), F = Math.hypot(r.last.force[0], r.last.force[1]), th = Math.hypot((r.last.pIn + 900 * U * U) * A, 900 * U * U * A);
    chk('Bend flow: force with the slug body filling the bend', F / th, 1, 0.04, 'the density wave has passed the outlet: momentum theorem with the slug density');
    chk('Bend flow: the passing slug front raises the force above the film value', Math.max(...r.f) / r.f[0] > 4 ? 1 : 0, 1, 0, `density ratio 6: peak ${fmt(Math.max(...r.f), 4)} N against ${fmt(r.f[0], 4)} N before the slug`); }
  // ---- space frame and the two-way coupled spool
  { const L = 3, n = 6, dir = [1 / 3, 2 / 3, 2 / 3], nodes = Array.from({ length: n + 1 }, (_, i) => dir.map((c) => (c * L * i) / n)), EI = 2e6, fr = frame3d({ nodes, elems: Array.from({ length: n }, (_, e) => ({ a: e, b: e + 1, EA: 1e9, EIy: EI, EIz: EI, GJ: 1e6, m: 10, rIp: 0.01 })), fix: [{ node: 0 }] }), f = new Float64Array(fr.ndof); f[6 * n] = 200; f[6 * n + 1] = -100;
    const s = frameStatic(fr, f); chk('Space frame: tip deflection of a skew cantilever', Math.hypot(s.u[6 * n], s.u[6 * n + 1], s.u[6 * n + 2]), (100 * Math.sqrt(5) * L ** 3) / (3 * EI), 1e-9, 'P·L³/(3·EI), load normal to the axis');
    chk('Space frame: reactions balance the load', Math.hypot(s.sumReaction[0] + 200, s.sumReaction[1] - 100, s.sumReaction[2]), 0, 1e-6, 'ΣR + ΣF = 0'); chk('Space frame: first frequency of the cantilever', frameModes(fr, 1).f[0], (1.8751040687 ** 2 / (2 * Math.PI)) * Math.sqrt(EI / (10 * L ** 4)), 2e-3, '(1.8751²/2π)·√(EI/(m·L⁴)), six elements'); }
  { const base = { D: 0.254, R: 0.381, angle: 90, legs: 6, U: 0, rho: 1000, mu: 1e-3, ns: 36, nn: 4, t: 0.0159, E: 207e9, nu: 0.3, rhoSteel: 7850, turbulence: false, nEl: 18 }, pre = bendFsi({ ...base, mode: 'one', dt: 1, tEnd: 1 }), T = 1 / pre.modes.filled[0];
    const r = bendFsi({ ...base, mode: 'two', dt: T / 48, tEnd: 2.3 * T, d0: { mode: 0, amplitude: 1e-3 } }), e = r.ux.map((x, i) => x - r.uy[i]), z = zero(r.t, e);
    chk('Coupled spool: frequency with the liquid resolved by the flow solver (two-way)', (z.length - 1) / (2 * (z[z.length - 1] - z[0])), pre.modes.filled[0], 0.015 * pre.modes.filled[0], `beam model with the mass of the contents: ${fmt(pre.modes.filled[0], 5)} Hz; empty ${fmt(pre.modes.empty[0], 5)} Hz; the steel-only beam is coupled to the flow by Aitken iterations (${fmt(r.iterations.mean, 2)} per step)`);
    chk('Coupled spool: work of the flow on the beam equals its energy gain', r.energy.gain / r.energy.work, 1, 0.02, 'Σ ½(fₙ + fₙ₊₁)·Δd against strain + kinetic energy (average-acceleration rule)');
    const s = bendFsi({ ...base, turbulence: true, U: 4, mode: 'one', dt: 0.004, tEnd: 0.8, zeta: 0.3 }), q = s.equilibrium;
    chk('Coupled spool: anchor reactions balance the fluid tractions (one-way, steady flow)', Math.hypot(q.reaction[0] + q.fluid[0], q.reaction[1] + q.fluid[1]) / Math.hypot(q.fluid[0], q.fluid[1]), 0, 0.01, 'resultant of the mapped wall tractions against the reactions of the beam model once the start-up transient has decayed'); }
  // ---- solid elements
  { const X = Float64Array.from([0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0.1, 0, 1, 1.1, 0, 1, 1.1, 1, 1, 0.1, 1, 1]), fixed = new Uint8Array(24), f = new Float64Array(24); // one brick sheared sideways: the loaded faces are inclined
    fixed[0] = fixed[1] = fixed[2] = 1; fixed[4] = fixed[5] = 1; fixed[11] = 1; const tr = (x, n) => [1e6 * n[0], 0, 0]; hexFaceLoad(X, [{ nodes: [1, 2, 6, 5] }], tr, f, [9, 0.5, 0.5]); hexFaceLoad(X, [{ nodes: [0, 3, 7, 4] }], tr, f, [-9, 0.5, 0.5]);
    const r = hexSolve({ X, elems: [[0, 1, 2, 3, 4, 5, 6, 7]], E: 2e11, nu: 0.3, fixed, f, tol: 1e-13 });
    chk('Solid elements: patch test, uniform tension on a sheared brick', Math.max(...r.stress.map((q) => Math.max(Math.abs(q[0] - 1e6), Math.abs(q[1]), Math.abs(q[3])))) / 1e6, 0, 1e-6, 'σxx = 1 MPa and no other stress at every node, with the incompatible modes'); }
  { const ro = 0.1429, t = 0.0159, ri = ro - t, r = hexPipe({ ro, t, angle: 0, legs: [0.1, 0], nLeg: [2, 0], nc: 12, nt: 3, E: 207e9, nu: 0.3, p: 30e6, ends: 'plane' }), la = lame(30e6, 0, ri, ro, ri);
    chk('Solid elements: hoop stress of the pressurised pipe against Lamé', r.hoop.inner[3] / la.hoop, 1, 0.01, 'plane strain, 12 × 3 elements on half the circumference');
    chk('Solid elements: residual of the conjugate-gradient solution', r.fe.residual, 0, 1e-7, `${r.fe.iterations} iterations, ${r.fe.ndof} degrees of freedom`);
    const b = hexPipe({ ro, t, angle: 0, legs: [2, 0], nLeg: [16, 0], nc: 8, nt: 1, E: 207e9, nu: 0.3, M: 1e5 });
    chk('Solid elements: rotation of a pipe under an end moment against beam theory', b.flexibility, 1, 0.01, 'θ·E·I/(M·L); 8 × 1 elements on half the circumference');
    chk('Solid elements: bending stress M·r_o/I', b.stressFactor.axial, 1, 0.02, 'largest longitudinal stress on the surface'); }
  { const ef = elbowFactors({ t: 0.0159, R: 0.381, r: 0.13495 }), r = hexPipe({ ro: 0.1429, t: 0.0159, R: 0.381, angle: 90, legs: [0.6, 0.6], nc: 12, nt: 2, nArc: 14, nLeg: 6, E: 207e9, nu: 0.3, M: 5e4, tol: 1e-7 });
    chk('Elbow factors: pipe characteristic h = t·R/r²', ef.h, (0.0159 * 0.381) / 0.13495 ** 2, 1e-12, 'ASME B31 form'); chk('Elbow factors: flexibility 1.65/h and stress intensification 0.9/h^⅔', ef.k * ef.h + ef.iIn * ef.h ** (2 / 3), 2.55, 1e-12, 'k·h + i·h^⅔ = 1.65 + 0.9');
    chk('Solid elbow: flexibility factor against the reference solid model', r.flexFactor, ELBOW_REF.flexFactor, 0.08 * ELBOW_REF.flexFactor, `in-plane bending, 90°, R = 1.5 D, legs of 0.6 m held at one end: CalculiX ${fmt(ELBOW_REF.flexFactor, 4)} (20-node elements), ASME B31 ${fmt(ef.k, 4)}, von Kármán ${fmt(ef.kKarman, 4)}; 12 × 2 × 26 elements here`);
    chk('Solid elbow: peak hoop stress factor against the reference solid model', r.stressFactor.hoop, ELBOW_REF.hoopFactor, 0.18 * ELBOW_REF.hoopFactor, `ovalisation bending of the wall over M·r_o/I; the eight-node elements approach it from below (−20 % with 400 elements, −7 % with 1000): CalculiX ${fmt(ELBOW_REF.hoopFactor, 4)}; twice the ASME B31 in-plane factor is ${fmt(2 * ef.iIn, 4)}`); }
  // ---- corrosion extensions
  { const sr = saltRetardation(10); chk('Salt retardation of the electrode reactions at 10 wt % NaCl', sr.cathodic + sr.anodic, 0.75 + 0.325, 1e-12, 'ratios of the fitted exchange current densities to their 1 wt % values');
    const a = mixedPotential({ T: 25, pCO2: 1, pH: 4, geometry: 'rce', rpm: 1000 }), b = mixedPotential({ T: 25, pCO2: 1, pH: 4, geometry: 'rce', rpm: 1000, lowT: false });
    chk('Mixed-potential model: activation enthalpies unchanged at and above 20 °C', a.rate / b.rate, 1, 1e-12, 'the low-temperature table only acts below 20 °C');
    const U = (Math.PI * 0.012 * 1000) / 60, mu = PROPS.viscosity(25, 0), rh = PROPS.density(25, 0), Dh = 9.31e-9 * (8.9e-4 / mu), Re = (rh * U * 0.012) / mu, Sh = 0.0791 * Re ** 0.7 * (mu / (rh * Dh)) ** 0.356;
    chk('Mixed-potential model: H⁺ limiting current on a rotating cylinder', a.iLimH, ((Sh * Dh) / 0.012) * FARADAY * 0.1, 1e-9 * a.iLimH, 'Sh = 0.0791·Re^0.7·Sc^0.356, 12 mm cylinder at 1000 rpm, pH 4');
    const cold = mixedPotential({ T: 5, pCO2: 1, pH: 4, geometry: 'rce', rpm: 1000 }); chk('Mixed-potential model: corrosion falls below the limiting current in cold water', cold.icorr / (cold.iLimH + cold.iLimH2CO3) < 0.5 ? 1 : 0, 1, 0, 'charge-transfer control at 5 °C, as observed');
    const base = { T: 40, pCO2: 1, P: 20, pH: 4.5, S: 20, U: 1, d: 0.2 }, n0 = norsokM506({ T: 40, pCO2: 1, P: 20, pH: 4.5, S: 20 }).base;
    chk('Automatic corrosion model: NORSOK M-506 inside its range', co2Corrosion({ ...base, salt: 0.5 }).base, n0, 1e-12 * n0, '40 °C, pH 4.5, 0.5 wt % NaCl: no scaling');
    chk('Automatic corrosion model: continuous at the 5 °C limit', co2Corrosion({ ...base, T: 4.999, salt: 0 }).base / co2Corrosion({ ...base, T: 5, salt: 0 }).base, 1, 2e-3, 'the mixed-potential ratio tends to one');
    chk('Automatic corrosion model: brine lowers the rate', co2Corrosion({ ...base, salt: 20 }).base < 0.5 * n0 ? 1 : 0, 1, 0, 'salinity factor at 20 wt % NaCl'); }
  // ---- statistics: negative-binomial counts and the burst-model uncertainty
  { let s = 0, m1 = 0; for (let n = 0; n < 400; n++) { const q = nbPmf(n, 30, 8); s += q; m1 += n * q; } chk('Negative-binomial probabilities sum to one', s, 1, 1e-10, 'mean 30, shape 8'); chk('Negative-binomial mean', m1, 30, 1e-8, 'Σ n·p(n)');
    chk('Negative-binomial tends to Poisson for a large shape', nbPmf(25, 30, 1e9), Math.exp(25 * Math.log(30) - 30 - lnGamma(26)), 1e-9, 'p(25; 30)'); chk('Incomplete gamma function P(1, x) = 1 − e^−x', gammaP(1, 0.7), 1 - Math.exp(-0.7), 1e-13, 'x = 0.7');
    const g = rng(5), rows = Array.from({ length: 40 }, (_, i) => { const mu = 200 * Math.exp(-0.03 * (i - 19.5)), lam = mu * Math.exp(g.normal(0, 0.08)); let n = Math.round(lam + Math.sqrt(lam) * g.normal(0, 1)); return { t: i, n: Math.max(n, 0), e: 1 }; }), ft = nbFit(rows);
    chk('Negative-binomial fit recovers the trend of synthetic counts', ft.b, -0.03, 0.006, '40 seeded years, rate falling 3 % per year, 8 % scatter of the rate between years'); chk('Negative-binomial fit finds over-dispersion', ft.pearson > 1.5 ? 1 : 0, 1, 0, `Pearson χ²/(n − 2) = ${fmt(ft.pearson, 3)}, shape ${fmt(ft.k, 3)}`); }
  if (REF.PHMSA_INCIDENTS) { const liq = REF.PHMSA_INCIDENTS.rows.filter((r) => r.system === 'liquid').map((r) => ({ t: r.year, n: r.incidents, e: r.exposure })); chk('Frequency model: shape of the gamma–Poisson model from the incident counts', nbFit(liq).k, FREQ_SHAPE, 0.02 * FREQ_SHAPE, 'hazardous-liquid lines 2010–2024, moment value about the fitted trend'); }
  if (REF.BURST_TESTS) { const bf = burstModelFit(REF.BURST_TESTS.rows); chk('Burst-model uncertainty: log-normal fit of the fitting half (median)', bf.defect.mu, BURST_MODEL.defect.mu, 2e-4, `${bf.defect.n} tests with defects`); chk('Burst-model uncertainty: log-normal fit of the fitting half (scatter)', bf.defect.sigma, BURST_MODEL.defect.sigma, 2e-4, 'standard deviation of ln(test / prediction)');
    chk('Burst-model uncertainty: checking half inside the Kolmogorov–Smirnov limit', bf.defect.check.ks < bf.defect.check.ksCritical ? 1 : 0, 1, 0, `D = ${fmt(bf.defect.check.ks, 3)} against ${fmt(bf.defect.check.ksCritical, 3)} (5 %), ${bf.defect.check.n} tests; ${fmt(100 * bf.defect.check.coverage90, 3)} % inside the 90 % interval`); chk('Burst-model uncertainty: intact pipes, fitting half', bf.intact.mean, BURST_MODEL.intact.mean, 2e-4, `${bf.intact.n} tests; checking half ${fmt(bf.intact.check.mean, 4)} ± ${fmt(bf.intact.check.sd, 3)}`); }
}

// ---- provenance of the literature constants ---------------------------------------------------------------------------------
/** API 5L / ISO 3183 line-pipe grades: specified minimum yield and tensile strength (MPa), metric table. */
export const LINEPIPE = Object.freeze({ X42: [290, 415], X46: [320, 435], X52: [360, 460], X56: [390, 490], X60: [415, 520], X65: [450, 535], X70: [485, 570], X80: [555, 625] });
const D0 = '2026-10-08', U_C203 = 'https://fenix.tecnico.ulisboa.pt/downloadFile/1126518382266636/DNVGL-RP-C203_2016-Fatgiue.pdf', U_F101 = 'https://rules.dnv.com/docs/pdf/dnvpm/codes/docs/2013-10/OS-F101.pdf', S_F101 = 'DNV-OS-F101, Submarine Pipeline Systems, October 2013, Det Norske Veritas (the edition later renamed DNV-ST-F101)',
  U_M506 = 'https://00448349299399495787.googlegroups.com/attach/207f0a05d64cb52e/NORSOK%20CO2%20M-506.pdf', S_M506 = 'NORSOK standard M-506, CO2 corrosion rate calculation model, Rev. 2, June 2005, Standards Norway', U_ICMT = 'https://www.icmt.ohio.edu/documents/Journals2019/Review%20of%20the%20APIRP14%20Eerosional%20velocity%20equation%20Origin,applications,misuses,%20limitations%20and%20alternatives.pdf',
  S_ICMT = 'F. Madani Sani, S. Huizinga, K. A. Esaklul, S. Nesic, Review of the API RP 14E erosional velocity equation: origin, applications, misuses, limitations and alternatives, Wear 426–427 (2019) 620–636', U_O501 = 'https://wiki.pengtools.com/images/9/91/RP_O501_EROSIVE_WEAR_IN_PIPING_SYSTEMS.pdf', S_O501 = 'Det Norske Veritas, Recommended Practice RP O501, Erosive wear in piping systems, Revision 4.2 (2007)',
  U_ORNL = 'https://info.ornl.gov/sites/publications/Files/Pub126720.pdf', S_ORNL = 'B. Oland, M. Lower, S. Rose, Review of methods for determining the strength of corroded natural gas pipelines based on actual remaining wall thickness, Oak Ridge National Laboratory, ORNL/TM-2019/1192', U_NESIC = 'https://www.icmt.ohio.edu/documents/publications/8183.pdf';
const pv = (item, used, source, url, status, note) => ({ item, used, source, url, retrieved: D0, status, note });
export const PROVENANCE = [
  pv('S–N curves in air: m₁, log a₁, m₂, log a₂, thickness exponent k for classes B1–W3', 'SN_CURVES, snCycles()', 'DNVGL-RP-C203, Fatigue design of offshore steel structures, April 2016, Table 2-1', U_C203, 'verified', 'All 14 classes compared number by number: no difference. The exponents of C and C1 (0.05, 0.10) are those of the 2016 edition (0.15 in 2010/2011); later editions were not readable.'),
  pv('S–N curves in seawater with cathodic protection', 'snCycles(…, "cp")', 'DNVGL-RP-C203, April 2016, Table 2-2', U_C203, 'corrected', 'B1 and B2: log a₁ is the air value − 0.200, not − 0.400 (14.717 → 14.917 and 14.485 → 14.685; lives below 10⁶ cycles were a factor 1.58 too short). The twelve m = 3 classes are exactly air − 0.400 with the air curve beyond 10⁶ cycles, as coded.'),
  pv('S–N curves in seawater, free corrosion', 'SN_FREE, snCycles(…, "free")', 'DNVGL-RP-C203, April 2016, Table 2-4', U_C203, 'corrected', 'The rule "air − 0.477" was replaced by the tabulated log a (F3: 11.069 → 11.068; all others equal), and the thickness exponent of C and C1 in this table is 0.15 (was 0.05 and 0.10).'),
  pv('Scatter of S–N data (standard deviation of log N = 0.20) and girth-weld misalignment factor SCF = 1 + (3δ/t)·exp(−√(t/D))', 'SN_SD, params() scfWeld', 'DNVGL-RP-C203, April 2016, section 2.10.1 (eq. 2.10.1) and commentary F.5', U_C203, 'verified', 'Read in the text of the standard; no difference.'),
  pv('Resistance factors γ_m = 1.15, γ_SC (1.046 / 1.138 / 1.308 pressure containment; 1.04 / 1.14 / 1.26 other), α_U = 0.96, α_fab (1.00 / 0.93 / 0.85), strain factors γ_ε (2.0 / 2.5 / 3.3)', 'DNV_SC, DNV_GM, ALPHA_FAB, DNV_STRAIN, pressureDesign()', S_F101 + ', Tables 5-2 to 5-5 and 5-10', U_F101, 'verified', 'No difference. Re-read on 2026-10-09: Table 5-3 of the 2013 edition has three columns only (low, medium, high); the class "very high" appears in Table 2-5 with a footnote to Appendix F, Table F-2 (onshore location classes), which is an image in the mirrored copy and could not be read. The suite reuses the "high" factors for "very high": still an assumption.'),
  pv('Yield-stress derating of C-Mn steel with temperature', 'dnvDerating()', S_F101 + ', Figure 2 of section 5', U_F101, 'verified', 'Read from the plotted curve: 0 at 50 °C, 30 MPa at 100 °C, 50 MPa at 150 °C, 70 MPa at 200 °C — the coded line passes through all four. The figure ends at 200 °C; beyond it the suite extrapolates.'),
  pv('Pressure containment, collapse (elastic, plastic, cubic with out-of-roundness ≥ 0.5 %), propagating-buckle pressure 35·f_y·α_fab·(t/D)^2.5, integral buckle-arrestor crossover pressure', 'pressureDesign(), collapsePressure(), arrestorCrossover()', S_F101 + ', eqs 5.8, 5.10–5.13, 5.16–5.18', U_F101, 'verified', 'Every formula and constant matches. The arrestor design check p_e ≤ p_X/(1.1·γ_m·γ_SC) was added from the same clause.'),
  pv('Local buckling under combined loading (load-controlled): α_c, α_p, β, plastic capacities', 'localBuckling()', S_F101 + ', eqs 5.19–5.28 (cross-read with the October 2007 edition)', U_F101, 'verified', 'Interaction formula, M_p, S_p, α_c and α_p identical. The three-branch β (0.5 below D/t = 15, zero above 60) is the 2007 form; the 2013 edition keeps (60 − D/t)/90 and limits the criterion to 15 ≤ D/t ≤ 45. The external-pressure branch now passes α_fab to the collapse pressure (it was taken as 1).'),
  pv('Compressive strain capacity ε_c = 0.78·(t/D − 0.01)·(1 + 5.75·Δp/p_b)·α_h^−1.5·α_gw', 'strainCapacity(), lateralDesign()', S_F101 + ', eq. 5.30', U_F101, 'verified', 'Constants identical in the 2010 and 2013 editions (the 2007 edition has 5/(p_b·2/√3) instead of 5.75/p_b). α_h = 0.93 is the largest yield-to-tensile ratio allowed for C-Mn grades.'),
  pv('Nominal annual target failure probabilities by safety class (10⁻³, 10⁻⁴, 10⁻⁵, 10⁻⁶)', 'DNV_SC[…].target, reliability()', S_F101 + ', Table 2-5', U_F101, 'verified', 'Equal to the row "ULS / FLS / ALS, all other". The pressure-containment row is one to two orders lower (10⁻⁴–10⁻⁵ low … 10⁻⁷–10⁻⁸ very high) and the serviceability row is 10⁻², 10⁻³, 10⁻³, 10⁻⁴; the suite applies the general row to the series system.'),
  pv('Line-pipe grades: specified minimum yield and tensile strength', 'LINEPIPE; default material X65 = 448 / 531 MPa from the reference case', 'API Specification 5L (45th ed.) / ISO 3183 tensile tables as reproduced by a pipe supplier; DNV-OS-F101 Table 7-5', 'https://tubingchina.com/API-5L-PSL-1-Pipe-Mechanical-Properties-Tensile-Yield-Strength-Elongation.htm', 'corrected', 'The metric table gives X65 (L450) = 450 / 535 MPa; 448 / 531 MPa is the conversion of 65 / 77 ksi used in older editions (0.4 % and 0.8 % lower, on the safe side). The shared reference case keeps 448 / 531; the metric table is exported as LINEPIPE and a note is raised when the grade text and the yield strength disagree by more than 5 %.'),
  pv('Hobbs buckling constants k₁–k₅ (vertical mode, lateral modes 1–4 and ∞)', 'HOBBS, hobbs(), lateralBuckle()', 'R. E. Hobbs, In-service buckling of heated pipelines, J. Transp. Eng. 110 (1984) 175–189, constants as reproduced in M. V. Craveiro, MSc thesis, University of São Paulo (2017) and two further open reproductions', 'https://www.teses.usp.br/teses/disponiveis/3/3144/tde-06122017-082632/publico/MarinaVendlCraveiroOrig17.pdf', 'verified', 'k₁, k₂, k₃ identical to the code; k₄ (2.407e-3, 5.532e-3, 1.032e-2, 1.047e-2) and k₅ (0.06938, 0.1088, 0.1434, 0.1483) added from the same table. One reproduction prints 0.1438 for mode 4 (taken as a transposition). The 1984 paper itself is not openly readable.'),
  pv('Upheaval buckling: dimensionless download curve of Palmer et al. (1990)', 'upheavalDownload()', 'A. C. Palmer et al., Design of submarine pipelines against upheaval buckling, OTC 6335 (1990), as described in an Oil & Gas Journal article', 'https://www.ogj.com/home/article/17230929/method-yields-download-force-to-arrest-upheaval-buckling-in-offshore-lines', 'unverified', 'Only the long-imperfection branch (9.6, 343, Φ_L > 8.06) could be read in an opened text; 0.0646, 4.49, 5.68 and 88.35 were not found in any document that could be opened. The three branches are continuous at 4.49 and 8.06 to 0.0003.'),
  pv('Paris law of BS 7910 (simplified): C = 1.65e-11 m/cycle (MPa√m)^−3, m = 3 in air; 7.27e-11 in a marine environment; threshold 63 N/mm^1.5 = 2.0 MPa√m', 'inputs parisC, parisM, dKth; parisLife()', 'A. Mehmanparast, F. Brennan, I. Tavares, Fatigue crack growth rates for offshore wind monopile weldments in air and seawater, Materials & Design 114 (2017) 494–504 (reproducing BS 7910)', 'https://strathprints.strath.ac.uk/64487/1/Mehmanparast_etal_MD_2017_Fatigue_cracl_growth_rates_for_offshore_wind_monopile_weldments.pdf', 'verified', 'Air constant and exponent identical; the help text quoted 7.3e-11 for the marine constant (now 7.27e-11); the threshold is 1.99 MPa√m against 2 (0.4 %).'),
  pv('Edge-crack geometry factor 1.12 − 0.231r + 10.55r² − 21.72r³ + 30.39r⁴', 'edgeCrackY()', 'M. Surendran, S. Natarajan, S. Bordas, G. S. Palani, Linear smoothed extended finite element method, arXiv:1701.03997, eq. 33 (Brown–Srawley / Tada fit)', 'https://arxiv.org/pdf/1701.03997', 'verified', 'Identical; a second open source prints 21.71 and 30.38 (below 0.1 % on Y). The fit is quoted for a/t ≤ 0.6; the suite clamps at 0.7.'),
  pv('Free-span vortex-induced vibration: onset values, in-line and cross-flow response models, amplitude reduction, safety factors', 'vivScreen()', 'DNV-RP-F105, Free spanning pipelines, February 2006, sections 4.3–4.4 and Table 2-2', 'https://rules.dnv.com/docs/pdf/dnvpm/codes/docs/2006-02/RP-F105.pdf', 'verified', 'Every coded coefficient matches. Not modelled: turbulence, proximity, trench and current-ratio corrections (all taken as 1) and cross-flow amplitudes above 0.9 D; the Strouhal number 0.2 used for the shedding frequency is a textbook value, not a number of the document.'),
  pv('Flow-induced turbulence screening: likelihood-of-failure bands (0.3, 0.5, 1), kinetic-energy level 20 000 kg/(m·s²), flexible-support coefficients', 'fivScreen()', 'Energy Institute, Guidelines for the avoidance of vibration induced fatigue failure in process pipework, 2nd ed. (2008), as reproduced in two open journal papers', 'https://ijcet.evegenis.org/index.php/ijcet/article/download/807/1035/1910', 'verified', 'Bands, class frequencies and the flexible-support formula (α = 41.21·D + 49397, β = 0.0815·ln D − 1.3842, added to the code) read directly. The lower kinetic-energy band of 5000 was not seen in an opened text.'),
  pv('Flow-induced turbulence screening: F_v coefficients for stiff, medium-stiff and medium supports', 'fivScreen()', 'Energy Institute guideline (2008), worked example reproduced in an open journal paper', 'https://ijcet.evegenis.org/index.php/ijcet/article/download/807/1035/1910', 'unverified', 'The formulas were not found printed; evaluated at 355.6 mm they reproduce the published α (717139, 415493, 224732) and β (−0.79, −0.85, −0.85) of a worked example. That is a check at one diameter only.'),
  pv('Allowable vibration velocity of process pipework: "concern" and "problem" lines', 'VIB_LINES, vibrationLimits()', 'Energy Institute guideline (2008) chart as reproduced in S. Salmi, MSc thesis, KTH Royal Institute of Technology (2025), Figure 1', 'https://www.diva-portal.org/smash/get/diva2:1977108/FULLTEXT01.pdf', 'unverified', 'The coefficients (0.48017, 2.127612; 1.871083, 2.084547) were not found printed in an opened document. The lines they give (1.7 → 24 mm/s and 7.9 → 122 mm/s rms over 1–300 Hz) agree with the reproduced chart within reading accuracy. The allowable velocity can be entered directly instead.'),
  pv('de Waard–Milliams nomogram equation log V = 5.8 − 1710/T + 0.67·log f_CO2', 'deWaardMilliams({ model: "1991" })', 'de Waard & Milliams (CORROSION/91 paper 577) as summarised on an engineering reference page', 'https://midstreamcalculator.com/engineering/pipeline-ops/co2-corrosion-fundamentals.html', 'verified', 'Identical, but only a secondary web reference could be opened; the conference paper is not open.'),
  pv('CO2 fugacity coefficient log a = P·(0.0031 − 1.4/T), limited to 250 bar; glycol factor 10^(1.6·(log W − 2))', 'deWaardMilliams(), norsokM506()', S_M506 + ', clause 8.1 (citing de Waard, Lotz & Milliams 1991)', U_M506, 'verified', 'Identical. Added from the same standard: the glycol factor is 0.008 above 95 wt % glycol.'),
  pv('Protective-scale factor log F = 2400/T − 0.6·log f_CO2 − 6.7', 'deWaardMilliams()', 'Three secondary open sources (none the original paper)', 'https://revues.imist.ma/index.php/AJMET/en/article/download/61281/32193/181178', 'unverified', '2400 and 6.7 agree everywhere; the coefficient of log f_CO2 appears as 0.6 (as coded), 0.61 and 0.44 in the three sources. Largest difference between sources: 27 % of that coefficient.'),
  pv('pH of CO2-saturated water 3.71 + 0.00417·t − 0.5·log f_CO2', 'deWaardMilliams()', 'K. S. George, MS thesis, Ohio University (Institute for Corrosion and Multiphase Technology), 2003, chapter 6', 'https://www.icmt.ohio.edu/documents/thesis/Electrochemical%20investigation%20of%20carbon%20dioxide%20corrosion%20of%20mild%20steel%20in%20the%20presence%20of%20acetic%20acid_K.%20S.%20George_2003_MS.pdf', 'verified', 'Identical.'),
  pv('de Waard–Lotz–Dugstad (1995) resistance model: log V_r = 4.93 − 1119/T + 0.58·log f − 0.34·(pH − pH_CO2), V_m = 2.45·U^0.8/d^0.2·f', 'deWaardMilliams({ model: "1995" })', 'Open-source implementation citing CORROSION/95 paper 128 (TNO GEMINI, module dld_model); K. S. George, MS thesis, Ohio University, 2003', 'https://gemini-hvc.westeurope.cloudapp.azure.com/documentation/sphinx/_modules/gemini_model/corrosion/correlation/dld_model.html', 'verified', 'The coded pair 4.93 / 2.45 and all exponents are confirmed by the open-source code. The thesis gives the steel-composition version of the same paper, 4.84 / 2.8 with a carbon factor 1 + 4.5·C %, which is not coded: with it V_r is 19 % lower and V_m 14 % higher before the carbon factor.'),
  pv('NORSOK M-506: K_t table (9 temperatures), pH functions (21 pieces), rate equations for 5, 15 and 20–150 °C, wall shear stress, validity limits', 'M506_KT, M506_FPH, norsokM506(), norsokShear()', S_M506 + ', Tables 1–2, eqs 1–8 and 20–21', U_M506, 'verified', 'Newly implemented from the standard itself (it had been left out). Eight rates evaluated independently from the tables agree with the code to 4 digits. Rates, not constants, are interpolated between temperatures, as the standard prescribes. The current Rev. 3 (2017) was not openly readable.'),
  pv('Carbonate equilibrium for the pH of CO2-saturated water: K_H, K₁, K₂, K_W as functions of temperature, pressure and ionic strength', 'norsokPH()', S_M506 + ', eqs 9–18', U_M506, 'corrected', 'Two departures from the printed text: eq. 16 prints 168491.5/T², which gives pK₁ = −10.7 at 25 °C (1684915 gives 6.35 and is used), and the hydration constant 0.00258 is not multiplied onto that K₁, because doing so gives pH 5.2 instead of 3.9 for pure water under 1 bar CO2. The result agrees with the de Waard expression within 0.1 pH.'),
  pv('Mixed-potential model: exchange current densities, activation enthalpies, Tafel slopes, reversible potential of iron, CO2 hydration, diffusion coefficients, Berger & Hau mass transfer', 'mixedPotential(), sherwood()', 'S. Nesic, J. Postlethwaite, S. Olsen, An electrochemical model for prediction of corrosion of mild steel in aqueous carbon dioxide solutions, Corrosion 52 (1996) 280–294', U_NESIC, 'verified', 'i₀(H⁺) = 0.05 A/m² at pH 4, E_rev(Fe) = −0.488 V, Tafel slopes, K_hyd = 2.58e-3, Sh = 0.0165·Re^0.86·Sc^0.33 and D(H₂CO₃) = 1.3e-9 m²/s identical to the paper. The activation enthalpies (30, 50, 37.5 kJ/mol) and the hydration-rate expression are those of the later parameter set of the same group (2003–2004); the 1996 paper has 30, 30 and 40 kJ/mol. Completed on 2026-10-09 from the same paper: water reduction (i₀ = 3e-5 A/m² at 20 °C, 30 kJ/mol, same Tafel slope and reversible potential as H⁺) and the rotating-cylinder mass transfer Sh = 0.0791·Re^0.7·Sc^0.356 (Eisenberg et al., its eq. 21); its three stated rates (1, 2.5 and 3 mm/y at 20, 50 and 80 °C, pH 4) are in the reference data.'),
  { ...pv('Retardation of the charge-transfer reactions by chloride: reference exchange current densities at 1, 3, 10 and 20 wt % NaCl (H⁺ reduction 0.08, 0.08, 0.06, 0.03 A/m²; iron dissolution 4.0e-3, 3.5e-3, 1.3e-3, 9.0e-4 A/m²), used as ratios to the 1 wt % values; activity coefficient of H⁺ rising to about 3 at 20 wt %', 'SALT_I0, saltRetardation(), mixedPotential(), co2Corrosion()', 'F. Madani Sani, B. Brown, S. Nesic, A mechanistic study on the effect of salt concentration on uniform corrosion rate of pipeline steel in acidic aqueous environments, AMPP Corrosion 2021, paper C2021-16788, Table 5; activity coefficient: the same authors, NACE Corrosion 2019, paper C2019-13026', 'https://www.icmt.ohio.edu/documents/AMPP2021/C2021-16788%20Fazlollah%20Madani%20Sani.pdf', 'verified', 'Eight numbers read in Table 5 (rotating disc, 10 °C, pH 3, 1 bar CO2). Between the tabulated concentrations the factors are interpolated log-linearly; above 20 wt % they are held (a second study reports 25 wt % behaving like 20 wt %). The factors were fitted by their authors at 10 °C and pH 3 and are applied here at all temperatures: checked blind against salt series at 5 °C and 20 °C (reference data set co2-corrosion-lab).'), retrieved: '2026-10-09' },
  { ...pv('Activation enthalpies of the electrode reactions below 20 °C: H⁺ reduction 100 / 140 kJ/mol, H2CO3 reduction 100 / 140, iron dissolution 185 / 215, water reduction 120 / 140 at 5 °C / 1 °C (30, 50, 37.5, 30 kJ/mol from 20 °C up)', 'LOWT_DH, arrheniusVar(), mixedPotential()', 'H. Fang, Low temperature and high salt concentration effects on general CO2 corrosion for carbon steel, MS thesis, Ohio University, 2006, Table 8', 'https://www.icmt.ohio.edu/documents/thesis/Low%20Temperature%20and%20High%20Salt%20Concentration%20Effect%20on%20General%20CO2%20Corrosion%20for%20Carbon%20Steel_2006_MS.pdf', 'verified', 'Twelve numbers read in Table 8. They are used as tabulated (in the Arrhenius factor from the reference temperature, interpolated between 20, 5 and 1 °C), which reproduces the 5 °C rates of the thesis (0.29 against 0.25–0.3 mm/y) and gives 0.055 against 0.1–0.2 mm/y at 1 °C. The author fitted them to her own sweeps, so those rows are a calibration check. The thesis also gives adsorption-isotherm constants for the salt retardation; their concentration unit could not be established from the text and they are not used.'), retrieved: '2026-10-09' },
  { ...pv('Corrosion-control basis: availability model CR = A·CR_inhibited + (1 − A)·CR_uninhibited; inhibited rate 0.05 mm/y up to 70 °C and 0.2 mm/y at 120–150 °C (0.1 mm/y practical minimum of another operator); availability limited to 95 % (99 % with an enhanced injection system); carbon steel with inhibition not accepted above about 6 mm/y uninhibited rate; field experience 85–95 %', 'inputs corrControl, inhibAvail, inhibRate; degradation(); warnings', 'J. Hobbs, Reliable corrosion inhibition in the oil and gas industry, Health and Safety Laboratory for the Health and Safety Executive, Research Report RR1023 (2014), section 5.2, equations 4–5 and Table 4', 'https://web.archive.org/web/2020id_/https://www.hse.gov.uk/research/rrpdf/rr1023.pdf', 'verified', 'Read in the report (Crown copyright, Open Government Licence), which reports the practice of two operators and of NORSOK M-001. Between 70 °C and 120 °C the inhibited rate is interpolated linearly here: the report quotes only the two ends of that table. The default availability of 99 % is the upper end of the quoted practice and is flagged in the results with the 95 % value beside it. NORSOK M-001 itself could not be opened.'), retrieved: '2026-10-09' },
  { ...pv('Classification of single-sided pipeline girth welds by the mean hi-lo: root class E up to 1.0 mm, F up to 2.0 mm, F1 up to 3.0 mm (thickness exponent 0); stress concentration 1 + (3δ/t)·exp(−√(t/D)); improvement factors of Table 7-1 (grinding and TIG dressing 3.5, hammer peening 4.0 on life above 350 MPa yield) not applied to the root', 'params() weld class, input weldBasis', 'DNVGL-RP-C203, Fatigue design of offshore steel structures, April 2016, section 2.10.1, Table 2-5 and Table 7-1', U_C203, 'verified', 'Read in the text of the standard. The suite takes the cap equation (2.10.1) for the root as the standard allows when the weld shape is unknown.'), retrieved: '2026-10-09' },
  { ...pv('Burst-model uncertainty (test / prediction): defects log-normal with μ_ln = 0.16245, σ_ln = 0.17211 (mean 1.194, scatter 17.3 %); intact pipe normal 1.070 ± 0.071', 'BURST_MODEL, burstModelFit(), reliability()', 'Fitted here by maximum likelihood to every second test of the burst database of Amaya-Gómez et al. (Mendeley Data, CC0), 42 tests with defects and 5 intact pipes; the other half checks the fit', 'https://data.mendeley.com/datasets/77t4k43y7g/1', 'verified', 'Replaces the assumed normal(1.05, 10 %) and normal(1.00, 5 %). Checking half: 90.2 % inside the 90 % interval, Kolmogorov–Smirnov distance 0.195 against 0.212 at 5 %; the fit is recomputed in the verification checks. The intact-pipe value rests on ten tests only.'), retrieved: '2026-10-09' },
  { ...pv('Shape of the gamma–Poisson (negative-binomial) frequency model: 183 (variance / mean of yearly counts about 1.8 at 150 incidents per year)', 'FREQ_SHAPE, nbFit(), table “Failure frequency against public statistics”', 'Fitted here to the significant-incident counts of US hazardous-liquid pipelines 2010–2024 (PHMSA flagged incident files, public domain): log-linear trend, dispersion from the Pearson statistic', 'https://web.archive.org/web/20260117080439id_/https://www.phmsa.dot.gov/sites/phmsa.dot.gov/files/data_statistics/pipeline/PHMSA_Pipeline_Safety_Flagged_Incidents.zip', 'verified', 'One-year-ahead forecast intervals hold 16 of 20 yearly counts (constant-rate Poisson: 11 of 20). The gas-transmission counts show no dispersion beyond Poisson about their trend. The shape of a national population is applied to the rate uncertainty of a single line: an assumption.'), retrieved: '2026-10-09' },
  { ...pv('Elbow flexibility and stress-intensification factors: h = t·R/r², k = 1.65/h, i (in-plane) = 0.9/h^⅔, i (out-of-plane) = 0.75/h^⅔; von Kármán (12h² + 10)/(12h² + 1); pressure stiffening 1 + 6·(P/E)·(r/t)^(7/3)·(R/r)^(1/3)', 'elbowFactors(), bendFsi(), task “Pipe bend of 3-D solid elements”', 'ASME B31.3 Appendix D form of the factors, as generally quoted; no open copy was read', 'https://www.calculix.de/', 'unverified', 'Not checked against a printed source. Code-to-code check instead: CalculiX solid models of three 90° bends with 0.6 m legs give flexibility factors 4.23, 7.67 and 2.17 against 4.96, 8.67 and 2.48 from 1.65/h (11–15 % lower, as expected for short restrained tangents), and peak hoop-stress factors 3.89, 5.83 and 2.33 against 3.75, 5.43 and 2.36 from twice 0.9/h^⅔.'), retrieved: '2026-10-09' },
  { ...pv('Wave speed of a pressure pulse in a liquid-filled elastic tube: Moens–Korteweg √(E·h/(2·ρ·R)), with the restraint factor 1/(1 − ν²) and the liquid compressibility 1/c² = ρ/K + 1/c₀² (Korteweg); standing-wave dispersion with the radial inertia of the liquid (I₀/I₁), wall inertia and bending; viscous wave ω² = c₀²k²(1 − F₁₀(α))', 'tubeWaves(), verification of tubeFsi()', 'Derived here from the long-wave equations of the liquid and the shell equation of the wall (classical results of Moens, Korteweg, Lamb and Womersley); no source was read for them', 'https://www.calculix.de/', 'unverified', 'Used only as reference solutions. The coupled solver reproduces them independently: pulse speed within 0.5 %, standing-wave frequency within 0.1 %, viscous frequency within 2 % and decay rate within 10 %.'), retrieved: '2026-10-09' },
  { ...pv('Sand erosion of an elbow in air: six measured rates (6.5–80.3 mm/y at 11–23 m/s, 150 and 300 µm sand, 76.2 mm bore)', 'reference data ELBOW_EROSION (comparison only)', 'N. Khan et al., Numerical study of gas–sand two-phase flow erosion in a standard 90° elbow, Engineering Proceedings 45 (2023) 28, Table 1 (CC BY 4.0), reproducing the tests of Vieira et al. (2016)', 'https://www.mdpi.com/2673-4591/45/1/28', 'verified', 'Only six rows could be found under an open licence, fewer than a data set of the calibration tab needs. Comparison made here: at 23 m/s the DNV-RP-O501 and Oka values are 0.25–0.6 of the measurement and Finnie 1.0–2.4; at 11 and 15 m/s the DNV-RP-O501 bend procedure gives about one hundredth of the measured rate, because in air near atmospheric pressure its critical-particle-size ratio falls outside its range and the small-particle correction is applied at its limit; Salama is 12–19 times too high for dry gas. For low-pressure gas the erosion models of the suite are therefore not validated; dense-gas and liquid conditions remain without open data.'), retrieved: '2026-10-09' },
  pv('API RP 14E erosional velocity V_e = C/√ρ_m and C-factors (100 continuous, 125 intermittent, 150–200 inhibited or corrosion-resistant)', 'erosionalVelocity()', S_ICMT + ' (quoting API RP 14E, 5th ed.)', U_ICMT, 'corrected', 'C-factors agree. The SI conversion constant was 1.21951 and is 0.3048·√16.018463 = 1.2199 (0.03 % higher).'),
  pv('DNV-RP-O501 erosion: K = 2.0e-9, n = 2.6, ρ_t = 7800 kg/m³, angle function, straight-pipe equation, bend procedure (impact angle, critical particle ratio, C₁ = 2.5, unit factor)', 'O501, dnvAngle(), erosionRate(), impactWear()', S_O501, U_O501, 'verified', 'Constants, straight-pipe and bend procedures identical. Revision 4.2 gives the angle function as an 8-term polynomial; the coded trigonometric form is that of the 2015 edition as reproduced in the review by Madani Sani et al. (2019). Of the geometry-factor table only the values 1 and 2 could be read.'),
  pv('Salama (2000) erosion rate ER = W·V²·d/(S_m·D²·ρ_m), S_m = 5.5 for elbows', 'erosionRate({ model: "salama" })', S_ICMT + ', reproducing M. M. Salama, J. Energy Resour. Technol. 122 (2000) 71–77', U_ICMT, 'verified', 'Equation, units and S_m = 5.5 for elbows confirmed in two opened sources. The plugged-tee constant S_m = 68 could not be read in any opened source and remains unverified.'),
  pv('Finnie cutting model: angle function with K = 2, ψ = 2; share of the ideal volume removed', 'erosionRate({ model: "finnie" }), impactWear()', 'L. Del Cid, The Finnie model of erosive wear (2016), summarising I. Finnie, Erosion of surfaces by solid particles, Wear 3 (1960) 87–103', 'https://lizidelcidphd.com/2016/06/09/the-finnie-model-of-erosive-wear/', 'corrected', 'Angle function and the factor 1/(ψ·K) as coded. The default share of the ideal cutting volume was 0.1 and is 0.5 in the source (the coded default predicted one fifth of the Finnie volume). Secondary source only.'),
  pv('Oka et al. (2005) erosion model for silica sand: s₁ = 0.71, q₁ = 0.14, s₂ = 2.4, q₂ = −0.94, K = 65, k₁ = −0.12, k₂ = 2.3·Hv^0.038, k₃ = 0.19, 104 m/s, 326 µm', 'OKA, erosionRate({ model: "oka" }), impactWear()', S_ICMT + ', table of the Oka model constants', U_ICMT, 'verified', 'All constants identical. A conference paper (OTC-27233) lists K = 60 instead of 65 (8 %).'),
  pv('ASME B31G: Folias factor √(1 + 0.8·L²/Dt), flow stress 1.1·SMYS, parabolic area, long-defect limit L²/Dt = 20', 'b31g({ method: "b31g" })', S_ORNL, U_ORNL, 'verified', 'The standard writes A = 0.893·L/√(Dt), that is 0.797·L²/Dt, which the code rounds to 0.8 (below 0.2 % on the factor).'),
  pv('Modified B31G: 0.85·d·L area, flow stress SMYS + 69 MPa, two-term bulging factor, switch at L²/Dt = 50', 'b31g({ method: "modified" })', S_ORNL, U_ORNL, 'verified', 'Identical (the standard has 68.95 MPa).'),
  pv('DNV-RP-F101 single defect: Q = √(1 + 0.31·L²/Dt), capacity equation', 'b31g({ method: "dnv" | "dnvCap" })', 'DNV-RP-F101, Corroded pipelines, October 2010, sections 2.3 and 8.2', 'https://rules.dnv.com/docs/pdf/dnvpm/codes/docs/2010-10/RP-F101.pdf', 'verified', 'Q and the failure pressure of the allowable-stress format are identical. The best-estimate burst capacity carries a further factor 1.05, added as method "dnvCap" and used when comparing with burst tests.'),
  pv('Generic loss-of-containment frequency (default 0.5 per 1000 km·y: subsea well-stream pipelines in open sea)', 'input failRate; table "Failure frequency against public statistics"', 'International Association of Oil & Gas Producers, Risk assessment data directory: Riser & pipeline release frequencies, Report 434-4 (March 2010), Table 2.1 (re-analysis of the PARLOC 2001 North Sea data)', 'https://web.archive.org/web/20140903090935id_/http://www.ogp.org.uk/pubs/434-04.pdf', 'corrected', 'The default was 0.3 from recollection. The report recommends 5.0e-4 per km·y for well-stream and other small pipelines with unprocessed fluid, 5.1e-5 for processed oil or gas up to 24 inch, 1.4e-5 above 24 inch, and 9.1e-4 per riser·y for steel risers up to 16 inch. The 2019 edition was not openly readable. Onshore gas lines: 0.118 per 1000 km·y in 2013–2022 (EGIG 12th report).'),
  pv('External-damage frequency (default 1e-4 per year for the whole line)', 'input extFreq; fault tree', 'Same report, Tables 2.1 and 4.3 (external loads cause 38 % of offshore releases; 7.9e-4 per year inside a platform safety zone for lines up to 16 inch)', 'https://web.archive.org/web/20140903090935id_/http://www.ogp.org.uk/pubs/434-04.pdf', 'unverified', 'The default is an assumption for a deep-water line beyond trawling and anchoring depth and is not taken from a source. The North Sea statistics would give 38 % × 5.0e-4 per km·y = 1.9e-4 per km·y for a well-stream line in open sea, about 40 times the default for a 20 km line: enter a site-specific value where fishing or shipping reaches the route.'),
  pv('Probability of ignition of a release (default 5 %)', 'input pIgnite; event tree', 'EGIG, Gas pipeline incidents: 12th report of the European Gas Pipeline Incident Data Group (1970–2022), Doc. VA 23.0304 (2023), Tables 3 and 7', 'https://www.egig.eu/reports/$60/$178', 'corrected', 'Ignition followed 4.8 % of pinhole or crack releases, 2.2 % of holes and 14.3 % of ruptures; weighted with the 2018–2022 leak-size frequencies this is 5.4 %. The default was 2 %. These are onshore gas-line statistics; a subsea release ignites only if gas reaches a source at the surface.'),
  pv('In-line inspection performance: 10 % of the wall detected with 90 % probability, depth sizing ± 10 % of the wall at 90 % certainty', 'inputs pod, sizingMm; reliability() inspection plan', 'Pipeline Operators Forum, Specifications and requirements for in-line inspection of pipelines, POF 100 (November 2021), guidance to Appendix 4', 'https://pipelineoperators.org/cdn/276341a3-f5e6-49cf-9897-de6ab41bdd5a/POF%20100%20Specifications%20and%20requirements%20for%20ILI%20-%20Nov%202021.pdf', 'corrected', 'Typical values for high-resolution magnetic-flux-leakage tools (general corrosion; 15 % of the wall for pitting; sizing 10–15 % of the wall). The sizing tolerance is quoted at 90 % certainty, not 80 % as first coded: the default standard deviation changed from 0.078 to 0.061 of the wall thickness. The document gives typical values only; the tool specification of the contractor governs.'),
  pv('Sour-service severity regions (0.3 kPa H2S threshold; region boundaries against in-situ pH)', 'sourRegion()', 'ANSI/NACE MR0175 / ISO 15156-2:2015, clause 7.2.1 and Figure 1', 'https://fouladonline.ir/wp-content/uploads/2017/05/NACE-MR-0175-ISO-15156-2015.pdf', 'corrected', 'Between 0.3 and 1 kPa the boundary of regions 1 and 2 keeps falling along pH = 4.5 + log p (3.98 at 0.3 kPa); the code held it at 4.5, which classed up to 0.5 pH units too severely.'),
];

/*VALIDATION_BEGIN*/
// ---- reference data sets with the suite's blind predictions -------------------------------------------------------------------
/** Reference data sets of js/data/ref/integ.js with the engine prediction for every row. Nothing is fitted to these data. */
function buildValidation() {
  const sets = [], add = (d, target, model, tolerance, note, extraCols = []) => { if (d && Array.isArray(d.rows) && d.rows.length >= 8) sets.push({ id: d.id, title: d.title, quantity: d.quantity, unit: d.unit, kind: d.kind, source: d.source, columns: [...d.columns, ...extraCols], rows: d.rows, target, model, tolerance, note }); };
  // burst tests: DNV-RP-F101 best-estimate capacity with the measured (or specified) tensile strength
  const burstOf = (r, method = 'dnvCap') => b31g({ D: r.D / 1000, t: r.t / 1000, d: r.d / 1000, L: r.L / 1000, smys: r.yield * MPA, smts: r.uts * MPA, method }).pf / MPA;
  add(REF.BURST_TESTS, 'Pb', (r) => burstOf(r), { mape: 15, bias: 1.5 }, 'Blind prediction with the DNV-RP-F101 capacity equation (factor 1.05) and the tensile strength listed for each test; D/t from 15 to 95, grades X42 to X100, depth up to 80 % of the wall, real and machined defects. Mean absolute error about 12 %, mean under-prediction about 0.9 MPa (conservative). The largest under-predictions (30–45 %) are long, real corrosion patches whose listed maximum depth is far from their average depth. Modified B31G on the same tests: about 18 %; original B31G: about 28 %, both conservative.');
  if (REF.BURST_TESTS && REF.BURST_TESTS.rows.length >= 8) { // probabilistic check of the burst-model uncertainty on the half of the tests that was not used to fit it
    const bf = burstModelFit(REF.BURST_TESTS.rows), rows = bf.defect.check.ratios.map((x, i, a) => ({ rank: i + 1, ratio: sig(x, 5), pObs: sig((i + 0.5) / a.length, 5) })), old = rows.filter((r) => Math.abs(r.ratio - BURST_MODEL.old.mean) <= 1.645 * BURST_MODEL.old.sd).length;
    sets.push({ id: 'burst-model-uncertainty', title: 'Burst-model uncertainty: measured / predicted burst pressure of the checking half of the tests against the fitted distribution', quantity: 'Cumulative probability', unit: '–', kind: 'experiment', source: REF.BURST_TESTS.source,
      columns: [{ key: 'rank', label: 'Rank in the checking half' }, { key: 'ratio', label: 'Measured / predicted (DNV-RP-F101 failure pressure)', unit: '–' }, { key: 'pObs', label: 'Observed cumulative frequency', unit: '–' }], rows, target: 'pObs',
      model: (r) => Phi((Math.log(r.ratio) - BURST_MODEL.defect.mu) / BURST_MODEL.defect.sigma), tolerance: { maxAbs: 0.2, bias: 0.1 }, note: `Reliability diagram of the random variable “burst model” of the defect limit state. The tests with a defect (${bf.defect.n + bf.defect.check.n}) are split in the order of the database: tests 1, 3, 5, … fit a log-normal distribution by maximum likelihood (median ${fmt(Math.exp(bf.defect.mu), 4)}, mean ${fmt(bf.defect.mean, 4)}, scatter ${fmt(100 * bf.defect.cov, 3)} %), tests 2, 4, 6, … (${bf.defect.check.n}, listed here) check it: ${fmt(100 * bf.defect.check.coverage90, 3)} % of them lie inside its 90 % interval, ${bf.defect.check.below5} below its 5 % quantile (${fmt(0.05 * bf.defect.check.n, 2)} expected), Kolmogorov–Smirnov distance ${fmt(bf.defect.check.ks, 3)} against ${fmt(bf.defect.check.ksCritical, 3)} at 5 %. The distribution assumed before (normal, mean ${BURST_MODEL.old.mean}, scatter 10 %) holds only ${fmt((100 * old) / rows.length, 3)} % of the checking half inside its 90 % interval. Intact pipes (pressure-containment format with the measured strengths): fitting half ${fmt(bf.intact.mean, 4)} ± ${fmt(bf.intact.sd, 3)} (${bf.intact.n} tests), checking half ${fmt(bf.intact.check.mean, 4)} ± ${fmt(bf.intact.check.sd, 3)} (${bf.intact.check.n} tests).` }); }
  // pipe conveying water: frequency ratio from the gyroscopic eigenproblem of the beam model
  if (REF.PIPE_FLOW) { const P = REF.PIPE_FLOW.problem; let cache = null; const ctxOf = () => { if (!cache) { const md = beamModel({ L: P.span, EI: P.E * P.I, m: P.mass, n: 16 }), mo = beamModes(md, 4); cache = { md, mo, mi: modalIntegrals(md, mo) }; } return cache; };
    add(REF.PIPE_FLOW, 'fRatio', (r) => { const c = ctxOf(); return fluidFrequency(c.md, c.mo, P.massFluid, r.vRatio * P.vcTheory, c.mi) / c.mo.f[0]; }, { rmse: 0.09, maxAbs: 0.18 }, `Aluminium pipe ${P.odMm} mm × ${P.wallMm} mm, span ${P.span} m, water. Beam elements with the conveyed-fluid terms give a critical velocity within 0.1 % of the report's theory (${P.vcTheory} m/s; measured divergence ${P.vDivergenceMeasured} m/s) and a first frequency of 28.40 rad/s. The measured frequencies lie up to 0.1 above the prediction at high velocity for pipe 2 and the divergence row (frequency zero at 98.5 % of the theoretical critical velocity) is predicted at 0.17: root-mean-square difference 0.08. The table was read from a scanned page; the legibility of every row is listed.`); }
  // NAFEMS deep beam: flexural modes from the shear-deformable beam elements; the rod and shaft modes from their exact solutions
  if (REF.NAFEMS_BEAM) { const P = REF.NAFEMS_BEAM.problem, a = P.side, A = a * a, I = a ** 4 / 12, Gm = P.E / (2 * (1 + P.poisson)); let fl = null;
    const flex = () => (fl || (fl = beamModes(beamModel({ L: P.length, EI: P.E * I, m: P.density * A, n: 40, kGA: (5 / 6) * Gm * A, rhoI: P.density * I }), 3).f));
    add(REF.NAFEMS_BEAM, 'f', (r) => { const k = REF.NAFEMS_BEAM.rows.filter((q) => q.type === r.type && q.mode < r.mode).length; if (r.type === 'flexural') return flex()[Math.floor(k / 2)]; if (r.type === 'extensional') return ((2 * k + 1) * Math.sqrt(P.E / P.density)) / (4 * P.length); return ((2 * k + 1) * Math.sqrt((Gm * 0.1406 * a ** 4) / (P.density * (a ** 4 / 6)))) / (4 * P.length); }, { mape: 4, maxAbs: 30 },
      'Square beam 2 m × 2 m, 10 m long (length / depth = 5). Flexural modes come from the beam elements of the suite with shear deformation (shear factor 5/6) and rotary inertia; they are 0.4 %, 2.9 % and 9 % below the reference for the first, second and third flexural pair, because the rotary inertia is taken on the slope of the deflection rather than on the cross-section rotation. Pipeline spans (length / diameter above 20) are far from this regime. The beam elements have no axial or torsional degrees of freedom: the extensional and torsional rows are evaluated with the exact rod and shaft formulas (torsion constant 0.1406·a⁴) and only show that the benchmark definition is reproduced.'); }
  // collapse tests: the DNV collapse formula with the measured yield strength and out-of-roundness (no design floor on the out-of-roundness)
  add(REF.COLLAPSE_TESTS, 'Pc', (r) => collapsePressure({ D: r.D, t: r.t, E: 207000, nu: 0.3, fy: r.yield, ovality: r.ovality / 100, ovalityMin: 0 }).pc, { mape: 10, bias: 3 }, 'Seamless and electric-welded tubulars, D/t from 12 to 23, yield 380 to 800 MPa, measured out-of-roundness 0.04 to 0.5 %. Blind prediction with the collapse equation of the suite (elastic and plastic pressure combined through the out-of-roundness), the measured yield strength and the measured out-of-roundness; Young\'s modulus is not tabulated and is taken as 207 GPa. Mean absolute error about 6.5 %, mean over-prediction about 1.2 MPa: residual stress and wall eccentricity, which lower the collapse pressure, are not in the equation. With the 0.5 % design floor on the out-of-roundness the prediction is on the safe side by 1.9 MPa on average. The table was read from a scanned page.');
  // girth-weld fatigue: mean S–N curve (design curve + 2 standard deviations) of the default class for a single-sided weld root
  if (REF.GIRTH_WELD_FATIGUE) { const W = REF.GIRTH_WELD_FATIGUE, t = W.problem.wallMm, below = W.rows.filter((r) => r.N < snCycles(r.sLocal, 'F1', 'air', { t })).length;
    add(W, 'logN', (r) => Math.log10(snCycles(r.sLocal, 'F1', 'air', { t, shift: 2 * SN_SD })), { rmse: 0.9, bias: 0.85 }, `Mean curve of class F1 in air (DNV-RP-C203 design curve shifted by two standard deviations of log N) with the stress range that includes the measured misalignment factor of every weld; the two-slope curve is used as it stands. The pipes last longer than predicted, by 0.75 decades on average: class F1 is conservative for these welds (the authors compare them with the mean of class E of BS 7608), and ${below} of the ${W.rows.length} failures lie below the design curve (2.3 % expected). Loading was at high mean stress (ratio 0.2 to 0.7), failure a through-wall crack from the root. Runouts are not included.`); }
  // failure statistics: one-year-ahead forecasts of a negative-binomial model with trend, each fitted to the earlier years only
  if (REF.PHMSA_INCIDENTS) { const all = REF.PHMSA_INCIDENTS.rows, fc = new Map(); let inside = 0, insideP = 0, n = 0;
    for (const sys of ['liquid', 'gas']) { const rs = all.filter((r) => r.system === sys).map((r) => ({ t: r.year, n: r.incidents, e: r.exposure })); for (const r of rs.filter((q) => q.t > 2014)) { const tr = rs.filter((q) => q.t < r.t), pr = nbFit(tr).predict(r.t, r.e), cst = (sum(tr.map((q) => q.n)) / sum(tr.map((q) => q.e))) * r.e; fc.set(sys + r.t, pr); n++; if (r.n >= pr.lo && r.n <= pr.hi) inside++; if (Math.abs(r.n - cst) <= 1.645 * Math.sqrt(cst)) insideP++; } }
    const test = all.filter((r) => r.year > 2014).map((r) => { const pr = fc.get(r.system + r.year); return { ...r, lo: pr.lo, hi: pr.hi, inside: r.incidents >= pr.lo && r.incidents <= pr.hi ? 'yes' : 'no' }; });
    add({ ...REF.PHMSA_INCIDENTS, rows: test }, 'freq', (r) => fc.get(r.system + r.year).mu / r.exposure, { mape: 16 }, `Probabilistic check of the frequency model behind the risk section: for every year from 2015 the model is fitted to the earlier years only (negative-binomial counts, log-linear trend of the rate, dispersion from the Pearson statistic, uncertainty of the trend carried into the forecast) and the count of that year is compared with the 90 % forecast interval. ${inside} of the ${n} yearly counts fall inside (${fmt(0.9 * n, 3)} expected; the probability of this many or fewer is 0.13, so the intervals are not rejected), against ${insideP} of ${n} for a constant-rate Poisson forecast. All misses are liquid-line years 2017–2020 below the interval: the rate turned from rising to falling in 2016, which no trend fitted to the earlier years can anticipate. Significant incidents include releases at facilities, so these frequencies are not line-pipe loss-of-containment frequencies.`, [{ key: 'lo', label: 'Forecast interval, 5 %' }, { key: 'hi', label: 'Forecast interval, 95 %' }, { key: 'inside', label: 'Inside the interval' }]); }
  // CO2 corrosion on rotating cylinders: the automatic model (blind) and, as a second set, the mixed-potential model alone
  if (REF.CO2_CORROSION_LAB) { const L = REF.CO2_CORROSION_LAB, at = (r, model) => co2Corrosion({ model, T: r.T, pCO2: r.pCO2, P: 1, pH: r.pH, salt: r.nacl, geometry: 'rce', rpm: r.rpm, d: r.d / 1000 }).rate, mp = (model, f = () => true) => { const rs = L.rows.filter(f); return (100 * sum(rs.map((r) => Math.abs(at(r, model) - r.rate) / r.rate))) / rs.length; }, old = (100 * sum(L.rows.map((r) => Math.abs(deWaardMilliams({ T: r.T, pCO2: r.pCO2, P: 1, model: '1991', pH: r.pH }).rate - r.rate) / r.rate))) / L.rows.length;
    add(L, 'rate', (r) => at(r, 'auto'), { mape: 105, bias: 1.1 }, `Bare steel, no inhibitor, no protective film. Default (automatic) model: NORSOK M-506 at the nearest valid condition, scaled with the mixed-potential model below 5 °C and above 1 wt % NaCl; the wall shear of the rotating cylinder comes from its friction law. Nothing is fitted to these rows: the salt retardation comes from sweeps at 10 °C and pH 3 of another study, and only the low-temperature activation enthalpies were fitted (by their author) to the 1 °C and 5 °C rows of the thesis. Mean absolute error ${fmt(mp('auto'), 3)} % against ${fmt(old, 3)} % for the nomogram equation used before, ${fmt(mp('norsok'), 3)} % for NORSOK M-506 alone and ${fmt(mp('mechanistic'), 3)} % for the mixed-potential model alone. The remaining error of the default is the standard itself, which lies a factor of 1.5 to 4 above these film-free glass-cell rates at 5–80 °C (it was fitted to flow-loop data and is on the safe side); the salt series is followed (20 °C: measured 2.0 / 1.46 / 0.46 mm/y at 3 / 10 / 20 wt %, ratios 1 / 0.73 / 0.23; model ratios 1 / ${fmt(at({ T: 20, pCO2: 1, pH: 4, nacl: 10, rpm: 1000, d: 12 }, 'auto') / at({ T: 20, pCO2: 1, pH: 4, nacl: 3, rpm: 1000, d: 12 }, 'auto'), 2)} / ${fmt(at({ T: 20, pCO2: 1, pH: 4, nacl: 20, rpm: 1000, d: 12 }, 'auto') / at({ T: 20, pCO2: 1, pH: 4, nacl: 3, rpm: 1000, d: 12 }, 'auto'), 2)}). Rows in strong brine (≥ 10 wt %): ${fmt(mp('auto', (r) => r.nacl >= 10), 3)} %; rows at 1 °C: ${fmt(mp('auto', (r) => r.T < 5), 3)} %.`);
    add({ ...L, id: 'co2-corrosion-lab-mixed-potential', title: L.title + ' — mixed-potential model' }, 'rate', (r) => at(r, 'mechanistic'), { mape: 45, bias: 0.3 }, `The same measurements against the mixed-potential model alone (selectable as the corrosion model): mean absolute error ${fmt(mp('mechanistic'), 3)} %. Blind for the 20–80 °C rows and for all brine rows; the 1 °C and 5 °C rows in 3 wt % NaCl are the calibration data of the low-temperature activation enthalpies (error on those rows ${fmt(mp('mechanistic', (r) => r.T <= 5 && r.nacl <= 3), 3)} %, on the others ${fmt(mp('mechanistic', (r) => !(r.T <= 5 && r.nacl <= 3)), 3)} %).`); }
  // code-to-code: CalculiX reference models against the in-app finite elements on their default-size meshes
  if (REF.CALCULIX_SECTIONS) { const P = REF.CALCULIX_SECTIONS.problem, cache = {}, get = (k) => { if (cache[k]) return cache[k]; const d = P[k];
      if (k === 'cylinder') { const c = feCylinder({ ...d, nr: 2, nz: 1 }); return (cache[k] = { uBore: c.uBore * 1e6, uOuter: c.uOuter * 1e6, hoopBore: c.hoopBore / MPA, hoopOuter: c.hoopOuter / MPA }); }
      if (k === 'groove') { const g = feGroove({ ...d, pe: 0, nz: 16, nr: 2 }), m = g.mesh; return (cache[k] = { vmMax: g.vmMax / MPA, hoopLigament: g.hoopLigament / MPA, uBoreCentre: g.fe.u[2 * m.id(0, 0)] * 1e6, hoopOuterCentre: g.fe.stress[m.id(0, m.NB - 1)][2] / MPA }); }
      const r = feRing({ ...d, dent: d.dent || undefined, nth: 32, nr: 2 }); return (cache[k] = { crownInner: r.crown.inner / MPA, crownOuter: r.crown.outer / MPA, dDx: r.dDiameter.x * 1e6, dDy: r.dDiameter.y * 1e6 }); };
    add(REF.CALCULIX_SECTIONS, 'value', (r) => get(r.case)[r.quantity], { mape: 3, maxAbs: 20 }, 'Blind comparison of the wall models of the suite (six-node triangles: 2 × 1 cells for the cylinder, 16 × 2 for the thinned band, 32 × 2 for the rings) with finer CalculiX models of the same geometry and loads. Thick cylinder and thinned band agree within 0.3 %, the out-of-round ring within 2 %, the dented ring within 6 % at the crown. The first comparison showed the crown stresses of the dented ring a factor of two to three too large in the suite: the graded mesh placed its mid-side nodes off the element centres; that was corrected and the values here are after the correction.'); }
  if (REF.CALCULIX_ELBOW) { const P = REF.CALCULIX_ELBOW.problem, cache = {}, get = (k) => { if (cache[k]) return cache[k]; const d = P[k], mesh = { ro: d.ro, t: d.t, R: d.R, angle: d.angle, legs: d.legs, nc: 10, nt: 2, nArc: d.R > 0.5 ? 16 : 12, nLeg: 5, E: d.E, nu: d.nu, tol: 1e-7 }, m = hexPipe({ ...mesh, M: d.M }), pr = hexPipe({ ...mesh, p: d.p, ends: 'both' }); return (cache[k] = { flexFactor: m.flexFactor, hoopFactor: m.stressFactor.hoop, axialFactor: m.stressFactor.axial, pressureHoopPeak: pr.peak.hoop / MPA }); };
    add(REF.CALCULIX_ELBOW, 'value', (r) => get(r.elbow)[r.quantity], { mape: 12 }, `Blind comparison of the eight-node solid elements of the suite (incompatible modes, 10 × 2 elements on half the circumference, 22 to 26 along) with CalculiX models of 2160 to 2560 twenty-node elements. Flexibility factors of the three bends: CalculiX ${['elbowA', 'elbowB', 'elbowC'].map((k) => fmt(REF.CALCULIX_ELBOW.rows.find((q) => q.elbow === k && q.quantity === 'flexFactor').value, 3)).join(' / ')}; ASME B31 1.65/h for the same bends ${['elbowA', 'elbowB', 'elbowC'].map((k) => fmt(elbowFactors({ t: P[k].t, R: P[k].R, r: P[k].ro - P[k].t / 2 }).k, 3)).join(' / ')}: with straight legs of only 0.6 m held at one end both solid models give 11 to 15 % less flexibility than the code factor, which assumes long flexible tangents. Mean difference between the coarse in-app mesh and the reference about 9 % (the coarse mesh is the stiffer one; refine it on the mesh tab).`); }
  if (REF.CALCULIX_SPAN) { const P = REF.CALCULIX_SPAN.problem, ri = P.ro - P.t, A = Math.PI * (P.ro ** 2 - ri ** 2), I = (Math.PI / 4) * (P.ro ** 4 - ri ** 4), cache = {}, get = (ends) => cache[ends] || (cache[ends] = beamModes(beamModel({ L: P.L, EI: P.E * I, m: P.rho * A, n: 40, ends: ends === 'pinned' ? ['pinned', 'pinned'] : ['fixed', 'fixed'], kGA: (0.5 * P.E * A) / (2 * (1 + P.nu)), rhoI: P.rho * I }), 4).f);
    add(REF.CALCULIX_SPAN, 'f', (r) => get(r.ends)[r.mode - 1], { mape: 1.5 }, 'Blind comparison of the beam elements of the suite (shear deformation with a shear factor of 0.5 for the thin tube, rotary inertia, 40 elements) with CalculiX beam elements expanded to solids. Empty steel pipe; the added mass of contents and water is a separate input of the span model.'); }
  return sets;
}
/*VALIDATION_END*/
let VALIDATION = null; // built on first use, so that loading the suite stays fast
// ---- suite object ---------------------------------------------------------------------------------------------------------
const mm = (x) => (isNum(x) && x > 0 ? x * 1000 : undefined), okNum = (x) => (isNum(x) ? x : undefined);
export default {
  id: 'integ', num: 6, title: 'Integrity, Loads, Risk & Engineering Assessment', short: 'Integrity · Risk', icon: '🛡️',
  tagline: 'Turns pressure, temperature, slugs and solids into stress, fatigue, corrosion, erosion, reliability and risk.',
  description: 'Pipe stress and code utilisation along the line, collapse, buckle arrestors and managed lateral buckling, slug forces on bends and equipment, the dynamic response of a span with two-way fluid–structure coupling (beam finite elements, pipe conveying slug flow), continuum and shell finite elements of the pipe wall, rainflow fatigue and crack growth, CO2 corrosion and sand erosion (correlations and tracked particles) with remaining life, defect assessment, structural reliability by FORM, SORM and sampling with an inspection-and-repair plan, and risk through fault tree, event tree, Markov and Bayesian models. Every number comes from the flow picture of the case, so the suite works before and after the upstream suites have been run.',
  guide: ['Check the pipe, material and design basis; pull the linked values from the network, flow, solids and operations suites.', 'Describe the span or jumper to be checked (length, supports, location) and the bends that take slug loads.', 'Set the corrosion and sand inputs; attach inspection data (minimum wall, wall map or a defect list) when available.', 'Run, then read the code checks, the critical locations and the recommendations; use the mesh tab to confirm the beam, time-step and sampling resolution.', 'Calibrate the corrosion multiplier, inhibitor efficiency, erosion multiplier and damping against coupon, probe and vibration data.', 'Set the design measures (buckle arrestors, buckle initiators, inspection plan, minimum shutdown pressure) and list equipment with its allowable loads; measurements entered on the set-up tab are compared with the model.'],
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
    'corrosion coupon data', 'ultrasonic wall-thickness measurements', 'intelligent-pigging data', 'erosion measurements', 'vibration measurements', 'field vibration measurements',
    // continuum and shell finite elements, coupled fluid–structure model, tracked particles
    'cauchy momentum equation', 'elasticity equations', 'shell equations', 'two-way fsi', 'erosion + particle cfd', 'element-order convergence',
    // resolved coupled solvers of the set-up tab (axisymmetric and plane Navier–Stokes coupled to shell, beam and solid elements)
    'navier–stokes + structural dynamics', 'cfd + fea',
    // state, interfaces and limits added as inputs
    'material condition', 'initial deposit thickness', 'guides', 'anchors', 'connections and equipment interfaces', 'allowable vibration or dynamic loading', 'pipe and equipment geometry',
    // quantities that can be estimated on the calibration tab (a fitted parameter and a measured column each) or by maximum likelihood in the run
    'material constitutive parameters', 'elastic/plastic properties', 'fatigue parameters', 'fracture parameters', 'support stiffness', 'soil-pipe interaction', 'erosion threshold', 'wall-shear correlations', 'slug-force coefficients', 'dynamic amplification', 'stress-concentration factors', 'defect-growth parameters', 'failure-probability models', 'reliability distributions', 'consequence-model parameters',
    // comparisons supported in the app: sourced reference data sets and the measurement table of the set-up tab (model prediction beside every measurement, error metrics, Brier score for failure records)
    'strain-gauge measurements', 'accelerometers', 'pressure transducers', 'load cells', 'full-scale bend-force tests', 'slug-force experiments', 'fatigue tests', 'tensile/material tests', 'burst tests', 'collapse tests', 'failure databases', 'inspection histories', 'actual damage/failure observations'],
  referenceOnly: [],
  equationsNote: 'Structure: linear-elastic beam elements (span), axisymmetric and plane continuum elements (pipe wall, thinned band, ring with out-of-roundness or a dent) and an axisymmetric thin cylindrical shell with thickness steps and ring stiffeners; no plasticity, so the collapse of a thinned wall is a net-section estimate from the elastic stress field. Fluid–structure interaction is two-way in reduced form: the fluid side is the one-dimensional slug model (conveyed mass, Coriolis and centrifugal terms of a pipe conveying fluid follow the slug train); resolved flow coupled to the structure runs as an additional solver of the set-up tab at browser resolution — axisymmetric Navier–Stokes flow in an elastic pipe (shell wall, strongly coupled) and plane flow through a bend with a slug as a density wave, mapped onto beam and three-dimensional solid elements, one-way or two-way; the plane and axisymmetric solutions do not resolve secondary flow, the gas–liquid interface or wall turbulence, and three-dimensional flow coupled to shell or solid models at production resolution remains with the external open-source solvers. Tracked-particle erosion uses a plane potential-flow field of the bend (no turbulence, no secondary flow), so it locates the scar and gives a rate without model factor; the DNV-RP-O501 equation remains the design value. Slug loads use unit-cell closures or the flow suite\'s slug summary, not a slug-tracking solution. Corrosion: NORSOK M-506 in its range, scaled by a mixed-potential model for cold water, brine and out-of-range chemistry (or the de Waard–Milliams family) for sweet service; the corrosion-control basis is the inhibitor-availability model with optional pH stabilisation and corrosion-resistant alloy (not valid when H2S controls the film); erosion models assume dilute sand. Lateral buckling uses Hobbs\' rigid-seabed solution with the feed-in shared between planned initiators. The inspection plan assumes that a repaired section is as new and that the wall-loss factor of a section stays the same over its life. Code formats (ASME B31.4/B31.8, DNV-ST-F101, DNV-RP-F101/F105/C203/O501, NORSOK M-506, ISO 15156) are implemented from their published equations for screening and concept design; they do not replace a code-compliant design verification. On the validation list a tick means that the comparison is supported in the app — by a sourced reference data set on the calibration tab or by the measurement table on the set-up tab — not that every such test has been reproduced.',
  inputs: INPUTS,
  presets: [
    { name: 'Reference deep-water tie-back', values: {} },
    { name: 'Same line without the design measures', values: { arrestors: 'none', buckleMgmt: 'none', inspInterval: 20, pMinShut: 1.013, corrControl: 'none', cladMode: 'none', phStab: 'off', weldBasis: 'manual', snClass: 'F1', scf: 1.45, spanMgmt: 'none' } },
    { name: 'Reference line with the coupled pressure-pulse solver', values: { task: 'tubeFsi' } },
    { name: 'Corroded and dented section (wall finite elements)', values: { ageY: 10, defects: [{ x: 6000, depth: 7, length: 250 }], dentMm: 12, ovality: 1.5, inspInterval: 3, hiLoMm: 1.5, residualStress: 200, strengthLoss: 5 } },
    { name: 'Sour high-CO2 service', values: { co2: 8, h2s: 0.5, phAct: 4.2, inhibEff: 97, corrAllow: 6, wtMm: 23.8, snEnv: 'free', safetyClass: 'high', material: 'API 5L X65 sour service', corrModel: '1995', corrControl: 'efficiency', cladMode: 'none' } },
    { name: 'Sandy late-life well', values: { sandKgD: 1500, sandUm: 350, ageY: 14, minWt: 13.6, erosionModel: 'governing', c14e: 125, inhibEff: 97, bends: [{ x: 0, angle: 90, radius: 0.38 }, { x: 18400, angle: 45, radius: 1.27 }, { x: 19560, angle: 90, radius: 0.38 }], geomFactor: 2, inspInterval: 2 } },
    { name: 'Slugging jumper fatigue', values: { slugMode: 'on', slugFreq: 0.08, slugLen: 25, slugVel: 6, slugHoldup: 0.9, spanLength: 22, spanX: 200, spanEnds: 'fixed-fixed', bendAt: 0.5, bendAngle: 90, odCoatMm: 0, damping: 1.5, weldBasis: 'manual', snClass: 'F1', scf: 1.5, dlf: 2, supportCapacity: 40, nSlugs: 16 } },
    { name: 'Corroded line: fitness for service', values: { ageY: 12, minWt: 12.4, inhibEff: 97, code: 'b318', designFactor: 0.72, riserFactor: 0.72, designPressure: 250, inspInterval: 3, defects: [{ x: 2400, depth: 5.6, length: 180 }, { x: 7350, depth: 11.6, length: 650 }, { x: 11800, depth: 3.9, length: 900 }, { x: 15200, depth: 11.4, length: 60 }] } },
    { name: 'Shallow-water span with strong current', values: { waterDepth: 60, currentSpeed: 0.45, waveHeight: 2, wavePeriod: 8, spanLength: 24, spanX: 6000, spanEnds: 'fixed-pinned', damping: 1.2, odCoatMm: 365.8, coatDensity: 2400, designPressure: 150, pShutIn: 140, safetyClass: 'medium', weldBasis: 'manual', snClass: 'F1', scf: 1.3, spanAxial: false } },
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
      { key: 'co2', value: okNum(fluid?.comp?.CO2), from: 'Case fluid: CO2 content (well stream, taken for the gas phase)' }, { key: 'h2s', value: okNum(fluid?.comp?.H2S), from: 'Case fluid: H2S content' }, { key: 'glycolWt', value: glycol ? okNum(fluid?.inhWt) : undefined, from: 'Case fluid: glycol in the water phase' }, { key: 'brineWt', value: isNum(fluid?.salinity) ? clamp(fluid.salinity, 0, 26) : undefined, from: 'Case fluid: salinity of the produced water' },
      { key: 'lineSupports', value: Array.isArray(n?.supports) && n.supports.length ? n.supports.filter((q) => q && isNum(q.x)).slice(0, 200).map((q) => ({ x: q.x, type: String(q.type || 'support'), k: nz(q.stiffness ?? q.k, 0) })) : undefined, from: 'Network: supports and restraints' },
      { key: 'slugFreq', value: has ? okNum(sg.freq) : undefined, from: 'Flow: slug frequency' }, { key: 'slugLen', value: has ? okNum(sg.length) : undefined, from: 'Flow: slug length' }, { key: 'slugVel', value: has ? okNum(sg.velocity) : undefined, from: 'Flow: slug velocity' }, { key: 'slugHoldup', value: has ? okNum(sg.holdupBody) : undefined, from: 'Flow: slug-body holdup' },
      { key: 'blockage', value: isNum(s?.blockage) ? clamp(s.blockage, 0, 1) : undefined, from: 'Solids: flow area lost' }, { key: 'plugX', value: okNum(s?.plugX ?? s?.onsetX), from: 'Solids: plug or onset location' }, { key: 'effIdMm', value: mm(s?.effectiveId), from: 'Solids: smallest effective bore' }, { key: 'plugProb', value: isNum(s?.plugProbability) ? clamp(s.plugProbability, 0, 1) : undefined, from: 'Solids: plug probability' },
      { key: 'depositRate', value: isNum(s?.waxRate) && s.waxRate >= 0 ? s.waxRate : undefined, from: 'Solids: wax deposition rate' },
      { key: 'blowdownMinT', value: okNum(o?.blowdownMinT), from: 'Operations: lowest blowdown temperature' }, { key: 'eventsShutdown', value: okNum(o?.eventsPerYear?.shutdowns), from: 'Operations: shutdowns per year' }, { key: 'eventsBlowdown', value: okNum(o?.eventsPerYear?.blowdowns), from: 'Operations: blowdowns per year' },
      { key: 'arrSpacing', value: isNum(n?.buckleArrestors?.spacing) && n.buckleArrestors.spacing > 0 ? n.buckleArrestors.spacing : undefined, from: 'Network: buckle-arrestor spacing' }, { key: 'initSpacing', value: isNum(n?.buckleInitiators?.spacing) && n.buckleInitiators.spacing > 0 ? n.buckleInitiators.spacing : undefined, from: 'Network: spacing of the buckle initiators' },
      { key: 'pMinShut', value: isNum(o?.blowdownEndP) && o.blowdownEndP >= 1 ? o.blowdownEndP : undefined, from: 'Operations: pressure at the end of a blowdown' },
      { key: 'deposit0', value: Array.isArray(s?.depositProfile?.total) && s.depositProfile.total.some((x) => isNum(x) && x > 0) ? Math.min(Math.max(...s.depositProfile.total.filter(isNum)) * 1000, 500) : undefined, from: 'Solids: largest deposit thickness' },
      { key: 'currentSpeed', value: okNum(site?.data?.currentSpeed), from: 'Site: current speed' }, { key: 'oilPrice', value: okNum(site?.data?.oilPrice), from: 'Site: oil price' },
    ].filter((it) => it.value !== undefined && it.value !== null && it.value !== '');
  },
  site: (site) => { const d = site?.data || {}; return [{ key: 'waterDepth', value: okNum(d.depth), from: 'Water depth at the site' }, { key: 'currentSpeed', value: okNum(d.currentSpeed), from: 'Current speed at the site' }, { key: 'waveHeight', value: okNum(d.waveHeight), from: 'Significant wave height at the site' }, { key: 'wavePeriod', value: okNum(d.wavePeriod), from: 'Wave period at the site' }, { key: 'tInstall', value: okNum(d.seabedTemp), from: 'Seabed temperature at the site' }, { key: 'oilPrice', value: okNum(d.oilPrice), from: 'Oil price' }].filter((it) => it.value !== undefined); },
  run,
  mesh: [
    { name: 'Beam elements on the span', keys: ['nElem'], min: 4, note: 'Hermite beam elements converge with the fourth power of the element length in the natural frequency; the dynamic stress converges with the second power.', metrics: [{ label: 'First natural frequency', unit: 'Hz', get: (r) => r.outputs.naturalFrequency }, { label: 'Peak dynamic bending stress', unit: 'MPa', get: (r) => r.outputs.peakDynamicStress }, { label: 'Static span stress', unit: 'MPa', get: (r) => r.outputs.staticSpanStress }] },
    { name: 'Time step of the dynamic response', keys: ['stepsPerCycle'], min: 8, note: 'Newmark-β with average acceleration is second-order accurate in the time step.', metrics: [{ label: 'Peak dynamic bending stress', unit: 'MPa', get: (r) => r.outputs.peakDynamicStress }, { label: 'Fatigue damage per year', unit: '1/y', get: (r) => r.outputs.fatigueDamagePerYear }] },
    { name: 'Wall finite elements (thinned band, ring, thick cylinder)', keys: ['feNr', 'feNz'], min: 1, note: 'Quadratic triangles converge with the fourth power of the element size in strain energy, linear triangles with the second power; switch the element order on the mesh tab to see the difference.', metrics: [{ label: 'Stress concentration of the thinned band', unit: '–', get: (r) => r.outputs.feStressConcentration }, { label: 'Ligament stress factor', unit: '–', get: (r) => r.outputs.feLigamentFactor }, { label: 'Hoop stress at the bore / Lamé', unit: '–', get: (r) => r.outputs.feLameRatio }, { label: 'Peak hoop stress of the ring', unit: 'MPa', get: (r) => r.outputs.ringPeakHoop }] },
    { name: 'Shell elements along the meridian', keys: ['shellN'], min: 10, note: 'Hermite elements on the elastic foundation of the hoop stiffness; the edge-bending length √(R·t) must be resolved.', metrics: [{ label: 'Shell ligament stress / solid elements', unit: '–', get: (r) => r.outputs.shellLigamentRatio }] },
    { name: 'Tracked particles and flow grid of the bend', keys: ['nPart', 'bendNs'], min: 10, note: 'The scar peak depends on the number of particles per bin; the mean over the scar converges faster.', metrics: [{ label: 'Peak tracked erosion rate', unit: 'mm/y', get: (r) => r.outputs.trackedErosionRate }, { label: 'Position of the peak', unit: '° of bend', get: (r) => r.outputs.trackedErosionAngle ?? 0 }] },
    { name: 'Coupled pulse grid', keys: ['fsiNx', 'fsiNr'], min: 20, note: 'Axial and radial cells of the coupled pressure-pulse model; select the task “Pressure pulse in the elastic pipe” first.', metrics: [{ label: 'Wave speed', unit: 'm/s', get: (r) => { if (!r.outputs.fsi) throw new Error('select the task “Pressure pulse in the elastic pipe”'); return r.outputs.fsi.waveSpeed; } }] },
    { name: 'Bend flow grid', keys: ['cfdNs', 'cfdNn'], min: 24, note: 'Cells along and across the bend; select one of the bend tasks first.', metrics: [{ label: 'Peak wall force', unit: 'kN', get: (r) => { if (!r.outputs.cfdFea) throw new Error('select the task “Bend: resolved flow”'); return r.outputs.cfdFea.peakForce; } }, { label: 'Wall force at the end of the run', unit: 'kN', get: (r) => { if (!r.outputs.cfdFea) throw new Error('select the task “Bend: resolved flow”'); return r.outputs.cfdFea.tractionForce; } }] },
    { name: 'Solid elements of the bend', keys: ['hexNc', 'hexNa'], min: 6, note: 'Elements around half the circumference and along the arc; select the task “Pipe bend of 3-D solid elements” first.', metrics: [{ label: 'Flexibility factor', unit: '–', get: (r) => { if (!r.outputs.solid) throw new Error('select the task “Pipe bend of 3-D solid elements”'); return r.outputs.solid.flexibilityFactor; } }, { label: 'Peak hoop stress factor', unit: '–', get: (r) => { if (!r.outputs.solid) throw new Error('select the task “Pipe bend of 3-D solid elements”'); return r.outputs.solid.stressFactorHoop; } }] },
    { name: 'Monte Carlo sample size', keys: ['nMC'], min: 500, note: 'Sampling error falls with the square root of the sample size, so this study shows statistical scatter rather than a formal order of convergence.', metrics: [{ label: 'Failure probability, importance sampling', unit: '–', get: (r) => r.outputs.mcPof }, { label: 'Reliability index, importance sampling', unit: '–', get: (r) => r.outputs.mcBeta }] },
  ],
  calibration: {
    note: 'Every row is one measurement; fill only the columns of its test type. Corrosion multiplier and inhibitor efficiency from coupon or wall-loss rates (rows with and without inhibitor separate the two); erosion multiplier and threshold velocity from elbow measurements at several velocities; damping from vibration amplitudes near resonance; modulus, yield strength and hardening exponent from tensile-test points; the S–N shift from fatigue tests and the Paris constants from crack-growth rates; the rotational support stiffness from a measured natural frequency; pipe–soil friction from pull tests; roughness (wall-shear correlation) from measured wall shear; the slug-force coefficient or the dynamic amplification from measured bend forces (fit one of the two); the stress concentration factor from hot-spot against nominal stress; the defect-growth multiplier from repeated inspections; the generic failure frequency from failure counts and exposure; the yield-strength distribution from ranked mill-test values (probability plot — maximum-likelihood fitting of samples is done in the run from the table “Measured samples for the distributions”); and the repair cost from incident costs.',
    params: [{ key: 'corrMult', label: 'Corrosion-model multiplier', lo: 0.1, hi: 4 }, { key: 'inhibEff', label: 'Inhibitor efficiency × availability (%)', lo: 30, hi: 99.9 }, { key: 'erosionMult', label: 'Erosion-model multiplier', lo: 0.05, hi: 20 }, { key: 'damping', label: 'Damping ratio (%)', lo: 0.2, hi: 15 },
      { key: 'eroVth', label: 'Erosion threshold velocity (m/s)', lo: 0, hi: 6 }, { key: 'eMod', label: 'Young\'s modulus (MPa)', lo: 150000, hi: 230000 }, { key: 'smys', label: 'Yield strength (MPa)', lo: 200, hi: 800 }, { key: 'hardenN', label: 'Strain-hardening exponent', lo: 4, hi: 40 }, { key: 'snShift', label: 'S–N curve shift Δlog a', lo: -1, hi: 1 }, { key: 'scf', label: 'Stress concentration factor', lo: 1, hi: 5 },
      { key: 'parisC', label: 'Paris coefficient C', lo: 1e-13, hi: 1e-9 }, { key: 'parisM', label: 'Paris exponent m', lo: 2, hi: 5 }, { key: 'spanKr', label: 'Rotational support stiffness (kN·m/rad)', lo: 100, hi: 1e7 }, { key: 'muLateral', label: 'Lateral pipe–soil friction', lo: 0.1, hi: 1.5 }, { key: 'muAxial', label: 'Axial pipe–soil friction', lo: 0.1, hi: 1.2 }, { key: 'roughUm', label: 'Wall roughness (µm)', lo: 1, hi: 500 },
      { key: 'slugCf', label: 'Slug-force coefficient', lo: 0.3, hi: 2.5 }, { key: 'dlf', label: 'Dynamic load factor', lo: 1, hi: 3 }, { key: 'defectGrowth', label: 'Defect-growth multiplier', lo: 0.2, hi: 5 }, { key: 'failRate', label: 'Generic failure frequency (per 1000 km·y)', lo: 0.01, hi: 5 }, { key: 'yieldBias', label: 'Mean yield / SMYS', lo: 1, hi: 1.3 }, { key: 'covYield', label: 'Yield strength scatter (%)', lo: 1, hi: 15 }, { key: 'repairCost', label: 'Repair or intervention cost (M$)', lo: 1, hi: 200 }],
    columns: [{ key: 'cT', label: 'Temperature', unit: '°C' }, { key: 'cP', label: 'Pressure', unit: 'bara' }, { key: 'cPco2', label: 'CO2 partial pressure', unit: 'bar' }, { key: 'cVliq', label: 'Liquid velocity', unit: 'm/s' }, { key: 'cInh', label: 'Inhibitor on (1) or off (0)', unit: '–' }, { key: 'corrRate', label: 'Measured corrosion rate', unit: 'mm/y' },
      { key: 'cVmix', label: 'Mixture velocity at the elbow', unit: 'm/s' }, { key: 'cRho', label: 'Mixture density', unit: 'kg/m³' }, { key: 'cSand', label: 'Sand rate', unit: 'kg/d' }, { key: 'eroRate', label: 'Measured erosion rate', unit: 'mm/y' }, { key: 'cFreq', label: 'Excitation frequency', unit: 'Hz' }, { key: 'cForce', label: 'Force amplitude', unit: 'N' }, { key: 'vibAmp', label: 'Measured mid-span amplitude', unit: 'mm' },
      { key: 'cStrain', label: 'Tensile test: strain', unit: '%' }, { key: 'cStress', label: 'Tensile test: measured stress', unit: 'MPa' }, { key: 'cSrange', label: 'Fatigue test: stress range', unit: 'MPa' }, { key: 'logN', label: 'Fatigue test: log10 of the cycles to failure', unit: '–' }, { key: 'cDK', label: 'Crack growth: ΔK', unit: 'MPa√m' }, { key: 'logGrowth', label: 'Crack growth: log10 of da/dN (m/cycle)', unit: '–' },
      { key: 'cSpanL', label: 'Modal test: span length', unit: 'm' }, { key: 'natFreq', label: 'Modal test: measured first frequency', unit: 'Hz' }, { key: 'cSoilW', label: 'Pull test: submerged weight', unit: 'N/m' }, { key: 'soilLat', label: 'Pull test: lateral resistance', unit: 'N/m' }, { key: 'soilAx', label: 'Pull test: axial resistance', unit: 'N/m' },
      { key: 'cTauV', label: 'Shear test: velocity', unit: 'm/s' }, { key: 'cTauRho', label: 'Shear test: density', unit: 'kg/m³' }, { key: 'cTauMu', label: 'Shear test: viscosity', unit: 'Pa·s' }, { key: 'wallShear', label: 'Shear test: measured wall shear', unit: 'Pa' }, { key: 'cSlugRho', label: 'Bend-force test: slug density', unit: 'kg/m³' }, { key: 'cSlugV', label: 'Bend-force test: slug velocity', unit: 'm/s' }, { key: 'bendForce', label: 'Bend-force test: measured peak force (90° bend)', unit: 'kN' },
      { key: 'cNomStress', label: 'Strain gauge: nominal stress', unit: 'MPa' }, { key: 'hotStress', label: 'Strain gauge: hot-spot stress', unit: 'MPa' }, { key: 'cDefD0', label: 'Defect: depth at the first inspection', unit: 'mm' }, { key: 'cDefYears', label: 'Defect: years to the next inspection', unit: 'y' }, { key: 'defDepth', label: 'Defect: depth at the next inspection', unit: 'mm' },
      { key: 'cExposure', label: 'Failure records: exposure', unit: 'km·y' }, { key: 'failures', label: 'Failure records: failures observed', unit: '–' }, { key: 'cQuantile', label: 'Mill tests: rank probability', unit: '–' }, { key: 'yieldAtP', label: 'Mill tests: yield strength at that rank', unit: 'MPa' }, { key: 'cDownDays', label: 'Incident: downtime', unit: 'd' }, { key: 'cProd', label: 'Incident: production lost', unit: 'bbl/d' }, { key: 'incidentCost', label: 'Incident: total cost', unit: 'M$' }],
    targets: [{ key: 'corrRate', label: 'Corrosion rate', unit: 'mm/y' }, { key: 'eroRate', label: 'Erosion rate', unit: 'mm/y' }, { key: 'vibAmp', label: 'Vibration amplitude', unit: 'mm' }, { key: 'cStress', label: 'Tensile stress', unit: 'MPa' }, { key: 'logN', label: 'Fatigue life (log10 N)', unit: '–' }, { key: 'logGrowth', label: 'Crack-growth rate (log10)', unit: '–' }, { key: 'natFreq', label: 'First natural frequency', unit: 'Hz' },
      { key: 'soilLat', label: 'Lateral soil resistance', unit: 'N/m' }, { key: 'soilAx', label: 'Axial soil resistance', unit: 'N/m' }, { key: 'wallShear', label: 'Wall shear stress', unit: 'Pa' }, { key: 'bendForce', label: 'Bend force', unit: 'kN' }, { key: 'hotStress', label: 'Hot-spot stress', unit: 'MPa' }, { key: 'defDepth', label: 'Defect depth', unit: 'mm' }, { key: 'failures', label: 'Failures', unit: '–' }, { key: 'yieldAtP', label: 'Yield strength', unit: 'MPa' }, { key: 'incidentCost', label: 'Incident cost', unit: 'M$' }],
    model: calModel,
    sample: CAL_SAMPLE,
    validationSample: CAL_VALIDATION,
  },
  verify,
  provenance: PROVENANCE,
  get validationData() { return VALIDATION || (VALIDATION = buildValidation()); },
  live: { key: 'opLog', label: 'Operating log (time h, inlet pressure bara, inlet temperature °C)', help: 'Follow an exported historian file: its pressure and temperature cycles are rainflow-counted into the fatigue spectrum and every excursion above the limits is counted.' },
};
