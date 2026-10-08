// Suite 7 — Economics, techno-economics and decision analysis.
// Converts the engineering consequences of suites 1–6 (pressure drop, arrival temperature, hydrate and wax exposure,
// slug surge, corrosion, failure probability, uptime) into cash flows, risk and ranked decisions. The file holds a
// financial-mathematics kernel, parametric CAPEX/OPEX models, a fiscal cash-flow engine, sampling and risk statistics,
// decision-analysis and optimisation routines, and the hybrid physics + economics studies that call the flow kernel.
import { brent, rng, lhs, quantile, histogram, mean, variance, sum, clamp, isNum, linspace, nelderMead, diffEvolution, lstsq, interp1 } from '../core/num.js';
import { fluidModel, inhibitorFor, INHIBITORS } from '../core/thermo.js';
import { uValue, slugUnit } from '../core/pipe.js';
import { caseLine, steadyCase } from '../core/caseflow.js';
import { BASE } from '../data/basecase.js';

export const BBL_PER_M3 = 6.28981077, GJ_PER_MMBTU = 1.05505585, MMBTU_PER_BOE = 5.8;
const MM = 1e6, Z90 = 1.2815515655446004;
const zeros = (n) => new Array(n).fill(0);

// ================================================================================================================
// 1. Time value of money
// ================================================================================================================
/** Future value of a present amount: compounding m times a year for n years at nominal annual rate r (fraction). */
export const futureValue = (pv, r, n, m = 1) => pv * (1 + r / m) ** (n * m);
/** Present value of a future amount (inverse of futureValue). */
export const presentValue = (fv, r, n, m = 1) => fv / (1 + r / m) ** (n * m);
/** Continuous compounding: FV = PV·exp(r·t). */
export const continuousFV = (pv, r, t) => pv * Math.exp(r * t);
/** Discount factor for a flow in year t. mid: flows of year t ≥ 1 fall at t − ½; continuous: exp(−r·t). */
export const discountFactor = (r, t, { mid = false, continuous = false } = {}) => { const tt = mid && t > 0 ? t - 0.5 : t; return continuous ? Math.exp(-r * tt) : (1 + r) ** -tt; };
/** Present value of an ordinary annuity of A per year for n years (optionally growing at g). */
export const annuityPV = (A, r, n, g = 0) => (Math.abs(r - g) < 1e-12 ? (A * n) / (1 + r) : (A * (1 - ((1 + g) / (1 + r)) ** n)) / (r - g));
/** Present value of a perpetuity A/(r − g); requires r > g. */
export const perpetuityPV = (A, r, g = 0) => { if (!(r > g)) throw new Error('A perpetuity needs a discount rate above its growth rate.'); return A / (r - g); };
/** Capital-recovery factor: the level annual amount over n years equivalent to 1 today. */
export const capitalRecovery = (r, n) => (Math.abs(r) < 1e-12 ? 1 / n : r / (1 - (1 + r) ** -n));
/** Fisher relation: real rate from a nominal rate and inflation (fractions). */
export const realRate = (nominal, inflation) => (1 + nominal) / (1 + inflation) - 1;
/** Net present value of cfs[k] (k = 0 is today) at an annual rate (fraction). */
export function npv(rate, cfs, opt = {}) {
  if (!(rate > -1)) throw new Error('The discount rate must be greater than −100 %.');
  let s = 0;
  for (let k = 0; k < cfs.length; k++) s += cfs[k] * discountFactor(rate, k, opt);
  return s;
}
/**
 * Internal rate of return. Scans −99 % … +10,000 % for every sign change of the NPV profile and polishes each bracket.
 * Returns { irr (fraction | null), roots: [], signChanges (of the cash-flow series), multiple (bool) }.
 * With several roots the smallest non-negative one is reported (else the root closest to zero) and `multiple` is set.
 */
export function irr(cfs) {
  const nz = cfs.filter((c) => c !== 0);
  let signChanges = 0;
  for (let i = 1; i < nz.length; i++) if (Math.sign(nz[i]) !== Math.sign(nz[i - 1])) signChanges++;
  if (!signChanges) return { irr: null, roots: [], signChanges, multiple: false };
  const f = (r) => { const q = 1 / (1 + r); let s = 0, d = 1; for (let k = 0; k < cfs.length; k++) { s += cfs[k] * d; d *= q; } return s * (r < 0 ? (1 + r) ** (cfs.length - 1) : 1); }; // scaled for r < 0 to avoid overflow
  const grid = [...linspace(-0.99, -0.01, 50), 0, ...linspace(-4, 2, 220).map((e) => 10 ** e)], roots = [];
  let fa = f(grid[0]);
  for (let i = 1; i < grid.length; i++) {
    const fb = f(grid[i]);
    if (fa === 0) roots.push(grid[i - 1]);
    else if (fa * fb < 0) roots.push(brent(f, grid[i - 1], grid[i], 1e-13));
    fa = fb;
  }
  if (!roots.length) return { irr: null, roots, signChanges, multiple: false };
  const pos = roots.filter((r) => r >= 0), pick = pos.length ? Math.min(...pos) : Math.max(...roots);
  return { irr: pick, roots, signChanges, multiple: roots.length > 1 };
}
/** Modified IRR: negative flows discounted at the finance rate, positive flows compounded at the reinvestment rate. */
export function mirr(cfs, financeRate, reinvestRate) {
  const n = cfs.length - 1;
  let pvNeg = 0, fvPos = 0;
  cfs.forEach((c, k) => { if (c < 0) pvNeg += c / (1 + financeRate) ** k; else fvPos += c * (1 + reinvestRate) ** (n - k); });
  return pvNeg < 0 && fvPos > 0 && n > 0 ? (fvPos / -pvNeg) ** (1 / n) - 1 : null;
}
/** Payback time (years from k = 0, linearly interpolated inside the recovery year); rate > 0 gives the discounted payback. null = never. */
export function payback(cfs, rate = 0, opt = {}) {
  let cum = 0, wasNeg = false;
  for (let k = 0; k < cfs.length; k++) {
    const d = cfs[k] * discountFactor(rate, k, opt), prev = cum;
    cum += d;
    if (cum < 0) wasNeg = true; else if (wasNeg) return k - 1 + -prev / d;
  }
  return wasNeg ? null : 0;
}
/** Equivalent annual value of a present amount over n years. */
export const equivalentAnnual = (pv, r, n) => pv * capitalRecovery(r, n);
/**
 * Depreciation schedule. method: 'sl' straight line | 'db' declining balance (rate per year, default 2/life, switching to
 * straight line when that is larger) | 'uop' units of production (units[] per year). Returns an array of `years` charges
 * that sums to base − salvage.
 */
export function depreciation(method, base, life, { years = life, rate, units, salvage = 0 } = {}) {
  const D = zeros(years), dep = base - salvage, n = Math.min(Math.max(1, Math.round(life)), years);
  if (!(dep > 0)) return D;
  if (method === 'uop' && units && sum(units.slice(0, years)) > 0) { const u = units.slice(0, years), tot = sum(u); u.forEach((x, i) => (D[i] = (dep * x) / tot)); return D; }
  if (method === 'db') { const r = rate ?? 2 / n; let book = dep; for (let i = 0; i < n; i++) { D[i] = i === n - 1 ? book : Math.min(book, Math.max(book * r, book / (n - i))); book -= D[i]; } return D; }
  for (let i = 0; i < n; i++) D[i] = dep / n;
  return D;
}

// ================================================================================================================
// 2. Production: Arps decline and annual profile
// ================================================================================================================
/** Arps decline at time t (years): rate q and cumulative Np (rate × years). b = 0 exponential, 0 < b < 1 hyperbolic, b = 1 harmonic. Di is the initial nominal decline (1/y). */
export function arps(qi, Di, b, t) {
  if (!(Di > 0) || t <= 0) return { q: qi, Np: qi * Math.max(t, 0) };
  if (b < 1e-9) return { q: qi * Math.exp(-Di * t), Np: (qi * (1 - Math.exp(-Di * t))) / Di };
  if (Math.abs(b - 1) < 1e-9) return { q: qi / (1 + Di * t), Np: (qi / Di) * Math.log(1 + Di * t) };
  return { q: qi * (1 + b * Di * t) ** (-1 / b), Np: (qi / ((1 - b) * Di)) * (1 - (1 + b * Di * t) ** ((b - 1) / b)) };
}
/** Annual volumes (rate × years; multiply a per-day rate by 365) for a plateau followed by Arps decline. */
export function productionProfile({ q0, plateau = 0, Di = 0.15, b = 0, life = 20 }) {
  const cum = (t) => (t <= plateau ? q0 * t : q0 * plateau + arps(q0, Di, b, t - plateau).Np);
  return Array.from({ length: life }, (_, j) => cum(j + 1) - cum(j));
}

// ================================================================================================================
// 3. Probability distributions, correlated sampling and risk statistics
// ================================================================================================================
/** Standard normal cumulative distribution (double-precision rational approximation). */
export function normCdf(x) {
  const a = Math.abs(x);
  let c;
  if (a > 37) c = 0;
  else {
    const e = Math.exp((-a * a) / 2);
    if (a < 7.07106781186547) {
      let n = 3.52624965998911e-2 * a + 0.700383064443688; n = n * a + 6.37396220353165; n = n * a + 33.912866078383; n = n * a + 112.079291497871; n = n * a + 221.213596169931; n = n * a + 220.206867912376;
      let d = 8.83883476483184e-2 * a + 1.75566716318264; d = d * a + 16.064177579207; d = d * a + 86.7807322029461; d = d * a + 296.564248779674; d = d * a + 637.333633378831; d = d * a + 793.826512519948; d = d * a + 440.413735824752;
      c = (e * n) / d;
    } else { let b = a + 0.65; b = a + 4 / b; b = a + 3 / b; b = a + 2 / b; b = a + 1 / b; c = e / b / 2.506628274631; }
  }
  return x > 0 ? 1 - c : c;
}
const normPdf = (x) => Math.exp(-0.5 * x * x) / Math.sqrt(2 * Math.PI);
/** Inverse standard normal (Acklam's rational approximation with one Halley refinement). */
export function normInv(p) {
  if (!(p > 0 && p < 1)) return p <= 0 ? -Infinity : Infinity;
  const a = [-3.969683028665376e1, 2.209460984245205e2, -2.759285104469687e2, 1.38357751867269e2, -3.066479806614716e1, 2.506628277459239], b = [-5.447609879822406e1, 1.615858368580409e2, -1.556989798598866e2, 6.680131188771972e1, -1.328068155288572e1];
  const c = [-7.784894002430293e-3, -3.223964580411365e-1, -2.400758277161838, -2.549732539343734, 4.374664141464968, 2.938163982698783], d = [7.784695709041462e-3, 3.224671290700398e-1, 2.445134137142996, 3.754408661907416];
  let x;
  if (p < 0.02425) { const q = Math.sqrt(-2 * Math.log(p)); x = (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1); }
  else if (p > 1 - 0.02425) { const q = Math.sqrt(-2 * Math.log(1 - p)); x = -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1); }
  else { const q = p - 0.5, r = q * q; x = ((((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q) / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1); }
  const e = normCdf(x) - p, u = e / normPdf(x);
  return x - u / (1 + (x * u) / 2);
}
function lgamma(z) {
  if (z < 0.5) return Math.log(Math.PI / Math.sin(Math.PI * z)) - lgamma(1 - z);
  const c = [0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];
  z -= 1;
  let x = c[0];
  for (let i = 1; i < 9; i++) x += c[i] / (z + i);
  const t = z + 7.5;
  return 0.5 * Math.log(2 * Math.PI) + (z + 0.5) * Math.log(t) - t + Math.log(x);
}
/** Regularised incomplete beta function I_x(a, b) (continued fraction). */
export function betaInc(x, a, b) {
  if (x <= 0) return 0;
  if (x >= 1) return 1;
  const cf = (x, a, b) => {
    const tiny = 1e-300, qab = a + b, qap = a + 1, qam = a - 1;
    let c = 1, d = 1 - (qab * x) / qap;
    if (Math.abs(d) < tiny) d = tiny;
    d = 1 / d;
    let h = d;
    for (let m = 1; m <= 300; m++) {
      const m2 = 2 * m;
      let aa = (m * (b - m) * x) / ((qam + m2) * (a + m2));
      d = 1 + aa * d; if (Math.abs(d) < tiny) d = tiny; c = 1 + aa / c; if (Math.abs(c) < tiny) c = tiny; d = 1 / d; h *= d * c;
      aa = (-(a + m) * (qab + m) * x) / ((a + m2) * (qap + m2));
      d = 1 + aa * d; if (Math.abs(d) < tiny) d = tiny; c = 1 + aa / c; if (Math.abs(c) < tiny) c = tiny; d = 1 / d;
      const del = d * c; h *= del;
      if (Math.abs(del - 1) < 3e-16) break;
    }
    return h;
  };
  const bt = Math.exp(lgamma(a + b) - lgamma(a) - lgamma(b) + a * Math.log(x) + b * Math.log(1 - x));
  return x < (a + 1) / (a + b + 2) ? (bt * cf(x, a, b)) / a : 1 - (bt * cf(1 - x, b, a)) / b;
}
/**
 * Build a distribution object from a table row { dist, lo, mode, hi }.
 * triangular / PERT: minimum, most likely, maximum. uniform: minimum, maximum. normal and lognormal: lo and hi are the
 * P10 and P90 values (mode is ignored). lo = hi gives a fixed value. Returns { kind, inv(u), mean, variance }.
 */
export function makeDist(d = {}) {
  const kind = String(d.dist ?? 'triangular').toLowerCase().trim(), a = +d.lo, b = +d.hi, name = d.name || d.id || kind;
  const bad = (why) => { throw new Error(`Distribution "${name}": ${why}.`); };
  if (!Number.isFinite(a) || !Number.isFinite(b)) bad('low and high must be numbers');
  if (b < a) bad('high must not be below low');
  if (b === a) return { kind: 'fixed', inv: () => a, mean: a, variance: 0 };
  const m = Number.isFinite(+d.mode) && d.mode !== null && d.mode !== '' ? +d.mode : (a + b) / 2;
  if (kind === 'uniform') return { kind, inv: (u) => a + (b - a) * u, mean: (a + b) / 2, variance: (b - a) ** 2 / 12 };
  if (kind === 'normal') { const mu = (a + b) / 2, sd = (b - a) / (2 * Z90); return { kind, inv: (u) => mu + sd * normInv(u), mean: mu, variance: sd * sd }; }
  if (kind === 'lognormal') { if (!(a > 0)) bad('a lognormal needs a positive P10'); const mu = (Math.log(a) + Math.log(b)) / 2, s = (Math.log(b) - Math.log(a)) / (2 * Z90); return { kind, inv: (u) => Math.exp(mu + s * normInv(u)), mean: Math.exp(mu + (s * s) / 2), variance: (Math.exp(s * s) - 1) * Math.exp(2 * mu + s * s) }; }
  if (m < a || m > b) bad('the most likely value must lie between low and high');
  if (kind === 'pert') {
    const al = 1 + (4 * (m - a)) / (b - a), be = 1 + (4 * (b - m)) / (b - a), mu = (a + 4 * m + b) / 6;
    let tab = null; // quantile table built on first use (400 intervals, linear in between)
    const nT = 400, inv = (u) => { if (!tab) tab = Array.from({ length: nT + 1 }, (_, i) => (i === 0 ? 0 : i === nT ? 1 : brent((x) => betaInc(x, al, be) - i / nT, 0, 1, 1e-12))); const s = clamp(u, 0, 1) * nT, i = Math.min(nT - 1, Math.floor(s)); return a + (b - a) * (tab[i] + (tab[i + 1] - tab[i]) * (s - i)); };
    return { kind, inv, mean: mu, variance: ((mu - a) * (b - mu)) / 7 };
  }
  if (kind !== 'triangular') bad(`unknown distribution type "${d.dist}" (use triangular, PERT, uniform, normal or lognormal)`);
  const fc = (m - a) / (b - a);
  return { kind, inv: (u) => (u < fc ? a + Math.sqrt(u * (b - a) * (m - a)) : b - Math.sqrt((1 - u) * (b - a) * (b - m))), mean: (a + m + b) / 3, variance: (a * a + b * b + m * m - a * b - a * m - b * m) / 18 };
}
/** Cholesky factor L (lower) of a symmetric positive-definite matrix, or null when it is not positive definite. */
export function cholesky(A) {
  const n = A.length, L = A.map(() => zeros(n));
  for (let i = 0; i < n; i++) for (let j = 0; j <= i; j++) {
    let s = A[i][j];
    for (let k = 0; k < j; k++) s -= L[i][k] * L[j][k];
    if (i === j) { if (!(s > 1e-12)) return null; L[i][i] = Math.sqrt(s); } else L[i][j] = s / L[j][j];
  }
  return L;
}
/** Pearson correlation of two equal-length arrays. */
export function correlation(x, y) {
  const mx = mean(x), my = mean(y);
  let sxy = 0, sxx = 0, syy = 0;
  for (let i = 0; i < x.length; i++) { sxy += (x[i] - mx) * (y[i] - my); sxx += (x[i] - mx) ** 2; syy += (y[i] - my) ** 2; }
  return sxx > 0 && syy > 0 ? sxy / Math.sqrt(sxx * syy) : 0;
}
/** Ranks (0 … n−1) of the values of an array. */
export const ranks = (a) => { const idx = a.map((_, i) => i).sort((p, q) => a[p] - a[q]), r = zeros(a.length); idx.forEach((i, k) => (r[i] = k)); return r; };
// a correlation matrix that is not positive definite is shrunk towards the identity until it is
function usableCorr(C, d) {
  if (!C) return { L: null, shrink: 0 };
  for (let lam = 0; lam <= 1.0001; lam += 0.05) {
    const L = cholesky(C.map((row, i) => row.map((x, j) => (i === j ? 1 : (1 - lam) * x))));
    if (L) return { L, shrink: lam };
  }
  return { L: null, shrink: 1, d };
}
/**
 * Correlated sample of d inputs. dists: distribution objects (makeDist); corr: d × d correlation matrix or null.
 * method 'mc': Gaussian copula (z = L·ε, u = Φ(z), x = F⁻¹(u)). method 'lhs': Latin hypercube with the Iman–Conover
 * rank re-ordering, which keeps every marginal stratified and imposes the target rank correlation.
 * Returns { X: n rows of d values, U (the uniforms), shrink (0 unless the matrix had to be repaired) }.
 */
export function sampleCorrelated(dists, { n = 1000, method = 'lhs', corr = null, seed = 1 } = {}) {
  const d = dists.length, R = rng(seed), plain = !corr || corr.every((row, i) => row.every((x, j) => i === j || x === 0)), { L, shrink } = plain ? { L: null, shrink: 0 } : usableCorr(corr, d);
  const cu = (u) => clamp(u, 1e-12, 1 - 1e-12);
  let U;
  if (method === 'mc') {
    U = Array.from({ length: n }, () => {
      const e = Array.from({ length: d }, () => R.normal());
      return e.map((_, i) => { if (!L) return cu(normCdf(e[i])); let z = 0; for (let k = 0; k <= i; k++) z += L[i][k] * e[k]; return cu(normCdf(z)); });
    });
  } else {
    U = lhs(n, d, seed);
    if (L) {
      const score = Array.from({ length: n }, (_, i) => normInv((i + 1) / (n + 1)));
      const Mc = Array.from({ length: d }, () => { const col = score.slice(); for (let i = n - 1; i > 0; i--) { const k = R.int(i + 1); [col[i], col[k]] = [col[k], col[i]]; } return col; });
      const T = Mc.map((a) => Mc.map((b) => correlation(a, b))), Q = cholesky(T);
      if (Q) {
        // S = L·Q⁻¹ maps the (nearly uncorrelated) scores to scores with the target correlation
        const Qi = Q.map(() => zeros(d));
        for (let j = 0; j < d; j++) for (let i = j; i < d; i++) { let s = i === j ? 1 : 0; for (let k = j; k < i; k++) s -= Q[i][k] * Qi[k][j]; Qi[i][j] = s / Q[i][i]; }
        const S = L.map((row) => zeros(d).map((_, j) => { let s = 0; for (let k = 0; k < d; k++) s += row[k] * Qi[k][j]; return s; }));
        for (let j = 0; j < d; j++) {
          const target = zeros(n);
          for (let i = 0; i < n; i++) { let s = 0; for (let k = 0; k < d; k++) s += S[j][k] * Mc[k][i]; target[i] = s; }
          const rk = ranks(target), sorted = U.map((row) => row[j]).sort((p, q) => p - q);
          for (let i = 0; i < n; i++) U[i][j] = sorted[rk[i]];
        }
      }
    }
  }
  return { X: U.map((row) => row.map((u, j) => dists[j].inv(cu(u)))), U, shrink };
}
/**
 * Risk statistics of a sample of outcomes (higher is better): mean, variance, sd, P10/P50/P90 (10th/50th/90th percentiles,
 * P10 is the low case), probLoss = P(x < 0), var = the (1 − alpha) quantile (value at risk expressed as an outcome),
 * cvar = mean of the outcomes at or below it (conditional value at risk / expected shortfall).
 */
export function riskStats(values, alpha = 0.95) {
  const s = [...values].sort((a, b) => a - b), n = s.length, q = (p) => { const x = (n - 1) * p, i = Math.floor(x), f = x - i; return i + 1 < n ? s[i] * (1 - f) + s[i + 1] * f : s[i]; };
  const v = q(1 - alpha), nTail = Math.max(1, Math.ceil((1 - alpha) * n - 1e-9));
  return { n, mean: mean(s), variance: variance(s), sd: Math.sqrt(variance(s)), p10: q(0.1), p50: q(0.5), p90: q(0.9), min: s[0], max: s[n - 1], probLoss: s.filter((x) => x < 0).length / n, var: v, cvar: mean(s.slice(0, nTail)), sorted: s };
}
/** Monte Carlo / Latin-hypercube propagation of f(x, i) over a set of distributions. Returns { X, y, shrink, ...riskStats }. */
export function monteCarlo(f, dists, opt = {}) {
  const ds = dists.map((d) => (typeof d.inv === 'function' ? d : makeDist(d))), { X, shrink } = sampleCorrelated(ds, opt), y = X.map((x, i) => f(x, i));
  return { X, y, shrink, ...riskStats(y, opt.alpha ?? 0.95) };
}
/** One-at-a-time sensitivity: f(vector) with each input moved to its low and high value. Sorted by swing (largest first). */
export function tornado(f, base, ranges) {
  const f0 = f(base);
  return ranges.map((r, i) => { const at = (x) => f(base.map((b, j) => (j === i ? x : b))), low = at(r.lo), high = at(r.hi); return { name: r.name, lo: r.lo, hi: r.hi, low, high, swing: Math.abs(high - low), base: f0 }; }).sort((p, q) => q.swing - p.swing);
}
/** Sobol sensitivity indices on the unit cube by the Saltelli sampling scheme: { first: [], total: [], variance, evals }. */
export function sobolIndices(f, d, N = 256, seed = 11) {
  const A = lhs(N, d, seed), B = lhs(N, d, seed + 7919), fA = A.map(f), fB = B.map(f), all = fA.concat(fB), V = variance(all) || 1e-300, first = [], total = [];
  for (let i = 0; i < d; i++) {
    let s1 = 0, st = 0;
    for (let k = 0; k < N; k++) { const x = A[k].slice(); x[i] = B[k][i]; const y = f(x); s1 += fB[k] * (y - fA[k]); st += (fA[k] - y) ** 2; }
    first.push(s1 / N / V); total.push(st / (2 * N) / V);
  }
  return { first, total, variance: V, evals: N * (d + 2) };
}
/** Standardised regression coefficients of y on the columns of X (rows are samples): { src: [], r2 }. */
export function standardisedRegression(X, y) {
  const d = X[0].length, cols = Array.from({ length: d }, (_, j) => X.map((r) => r[j])), mu = cols.map(mean), sd = cols.map((c) => Math.sqrt(variance(c)) || 1), my = mean(y), sy = Math.sqrt(variance(y)) || 1;
  const Z = X.map((r) => r.map((x, j) => (x - mu[j]) / sd[j])), yz = y.map((v) => (v - my) / sy), bcoef = lstsq(Z, yz);
  let sse = 0, sst = 0;
  Z.forEach((r, i) => { const p = r.reduce((s, x, j) => s + x * bcoef[j], 0); sse += (yz[i] - p) ** 2; sst += yz[i] ** 2; });
  return { src: bcoef, r2: sst > 0 ? 1 - sse / sst : 0 };
}
/** Annual price-multiplier path around the trend: 'gbm' geometric Brownian motion, 'ou' mean-reverting log price; both have expectation 1. */
export function pricePath(model, K, R, { sigma = 0.25, kappa = 0.3 } = {}) {
  const out = new Array(K).fill(1);
  if (model === 'gbm') { let w = 0; for (let k = 1; k < K; k++) { w += R.normal(); out[k] = Math.exp(sigma * w - 0.5 * sigma * sigma * k); } }
  else if (model === 'ou') { let x = 0; const e = Math.exp(-kappa), sdStep = sigma * Math.sqrt((1 - e * e) / (2 * kappa)); for (let k = 1; k < K; k++) { x = x * e + sdStep * R.normal(); const vr = (sigma * sigma * (1 - Math.exp(-2 * kappa * k))) / (2 * kappa); out[k] = Math.exp(x - 0.5 * vr); } }
  return out;
}

// ================================================================================================================
// 4. Decision analysis
// ================================================================================================================
/**
 * Roll back a decision tree. node: { name, type: 'decision' | 'chance', branches: [{ name, p, cost, value | node }] } or a
 * terminal { name, value }. Chance probabilities are normalised. Returns the evaluated tree: { name, type, emv, choice?, branches }.
 */
export function decisionTree(node) {
  if (!node.branches || !node.branches.length) return { name: node.name ?? '', type: 'end', emv: +node.value || 0 };
  const kids = node.branches.map((b) => { const sub = decisionTree(b.node || { name: b.name, value: b.value }), cost = +b.cost || 0; return { name: b.name, p: +b.p || 0, cost, emv: sub.emv - cost, node: sub }; });
  if (node.type === 'chance') {
    const ps = sum(kids.map((k) => k.p));
    if (!(ps > 0)) throw new Error(`Chance node "${node.name}" has no positive probabilities.`);
    kids.forEach((k) => (k.p /= ps));
    return { name: node.name, type: 'chance', emv: sum(kids.map((k) => k.p * k.emv)), probSum: ps, branches: kids };
  }
  const best = kids.reduce((a, k) => (k.emv > a.emv ? k : a), kids[0]);
  return { name: node.name, type: 'decision', emv: best.emv, choice: best.name, branches: kids };
}
/**
 * Value of information. prior[s], payoff[a][s], optional likelihood[k][s] = P(signal k | state s).
 * Returns { prior, ev: [], emv, best (action index), evWithPI, evpi, signalProb, posterior, signalAction, evWithII, evii }.
 */
export function valueOfInformation({ prior, payoff, likelihood }) {
  const ps = sum(prior);
  if (!(ps > 0)) throw new Error('The prior probabilities must sum to a positive number.');
  const P = prior.map((x) => x / ps), ev = payoff.map((row) => sum(row.map((x, s) => x * P[s]))), emv = Math.max(...ev), best = ev.indexOf(emv);
  const evWithPI = sum(P.map((p, s) => p * Math.max(...payoff.map((row) => row[s])))), out = { prior: P, ev, emv, best, evWithPI, evpi: evWithPI - emv, signalProb: [], posterior: [], signalAction: [], evWithII: emv, evii: 0 };
  if (likelihood) {
    let tot = 0;
    for (const Lk of likelihood) {
      const pk = sum(Lk.map((l, s) => l * P[s])), post = Lk.map((l, s) => (pk > 0 ? (l * P[s]) / pk : P[s])), evs = payoff.map((row) => sum(row.map((x, s) => x * post[s]))), vk = Math.max(...evs);
      out.signalProb.push(pk); out.posterior.push(post); out.signalAction.push(evs.indexOf(vk)); tot += pk * vk;
    }
    out.evWithII = tot; out.evii = tot - emv;
  }
  return out;
}
/** Certainty equivalent of equally likely outcomes under exponential utility U(x) = 1 − exp(−x/R); R is the risk tolerance. */
export function certaintyEquivalent(values, R) {
  if (!(R > 0) || !values.length) return mean(values);
  const m = Math.min(...values);
  return m - R * Math.log(mean(values.map((x) => Math.exp(-(x - m) / R))));
}
/** Exponential utility of an outcome for risk tolerance R. */
export const expUtility = (x, R) => 1 - Math.exp(-x / R);
/** Analytic Hierarchy Process: principal eigenvector of a pairwise comparison matrix. Returns { weights, lambdaMax, ci, cr, n }. */
export function ahp(Mx) {
  const n = Mx.length, RI = [0, 0, 0, 0.58, 0.9, 1.12, 1.24, 1.32, 1.41, 1.45, 1.49];
  if (n < 1 || Mx.some((r) => r.length !== n || r.some((x) => !(x > 0)))) throw new Error('The pairwise comparison matrix must be square with positive entries.');
  let w = zeros(n).map(() => 1 / n), lam = n;
  for (let it = 0; it < 500; it++) {
    const y = Mx.map((r) => r.reduce((s, x, j) => s + x * w[j], 0)), t = sum(y), w2 = y.map((x) => x / t), diff = Math.max(...w2.map((x, i) => Math.abs(x - w[i])));
    w = w2; lam = t;
    if (diff < 1e-14) break;
  }
  lam = mean(Mx.map((r, i) => r.reduce((s, x, j) => s + x * w[j], 0) / w[i]));
  const ci = n > 1 ? (lam - n) / (n - 1) : 0, ri = RI[Math.min(n, 10)];
  return { weights: w, lambdaMax: lam, ci, cr: ri > 0 ? ci / ri : 0, n };
}
/** Weighted-sum multi-criteria score with min–max normalisation. matrix[alt][crit]; benefit[crit] true when more is better. */
export function weightedSum(matrix, weights, benefit) {
  const ws = sum(weights), w = weights.map((x) => x / ws), nc = w.length, lo = zeros(nc).map((_, j) => Math.min(...matrix.map((r) => r[j]))), hi = zeros(nc).map((_, j) => Math.max(...matrix.map((r) => r[j])));
  const norm = matrix.map((r) => r.map((x, j) => (hi[j] - lo[j] < 1e-12 * Math.max(1, Math.abs(hi[j])) ? 1 : benefit[j] ? (x - lo[j]) / (hi[j] - lo[j]) : (hi[j] - x) / (hi[j] - lo[j]))));
  return { scores: norm.map((r) => sum(r.map((x, j) => x * w[j]))), norm, weights: w };
}
/** TOPSIS: closeness to the ideal solution (vector normalisation). Returns { closeness: [], order: [] (best first) }. */
export function topsis(matrix, weights, benefit) {
  const ws = sum(weights), w = weights.map((x) => x / ws), nc = w.length, len = zeros(nc).map((_, j) => Math.sqrt(sum(matrix.map((r) => r[j] * r[j]))) || 1);
  const V = matrix.map((r) => r.map((x, j) => (w[j] * x) / len[j])), col = (j) => V.map((r) => r[j]);
  const best = zeros(nc).map((_, j) => (benefit[j] ? Math.max(...col(j)) : Math.min(...col(j)))), worst = zeros(nc).map((_, j) => (benefit[j] ? Math.min(...col(j)) : Math.max(...col(j))));
  const closeness = V.map((r) => { const dp = Math.sqrt(sum(r.map((x, j) => (x - best[j]) ** 2))), dm = Math.sqrt(sum(r.map((x, j) => (x - worst[j]) ** 2))); return dp + dm > 0 ? dm / (dp + dm) : 0.5; });
  return { closeness, order: closeness.map((_, i) => i).sort((a, b) => closeness[b] - closeness[a]) };
}
/** Black–Scholes–Merton value of a European option on an asset paying a continuous yield q. type: 'call' | 'put'. */
export function blackScholes({ S, K, r, sigma, T, q = 0, type = 'call' }) {
  if (!(S > 0) || !(K > 0)) return type === 'call' ? Math.max(S - K, 0) : Math.max(K - S, 0);
  if (!(sigma > 0) || !(T > 0)) return type === 'call' ? Math.max(S * Math.exp(-q * T) - K * Math.exp(-r * T), 0) : Math.max(K * Math.exp(-r * T) - S * Math.exp(-q * T), 0);
  const d1 = (Math.log(S / K) + (r - q + 0.5 * sigma * sigma) * T) / (sigma * Math.sqrt(T)), d2 = d1 - sigma * Math.sqrt(T);
  return type === 'call' ? S * Math.exp(-q * T) * normCdf(d1) - K * Math.exp(-r * T) * normCdf(d2) : K * Math.exp(-r * T) * normCdf(-d2) - S * Math.exp(-q * T) * normCdf(-d1);
}
/**
 * Cox–Ross–Rubinstein binomial lattice. { S, K, r, sigma, T, q, steps, type: 'call' | 'put', american, payoff(S) }.
 * A custom payoff(S) replaces the vanilla one (used for the option to expand). Returns { value, u, d, p, steps }.
 */
export function binomialOption({ S, K, r, sigma, T, q = 0, steps = 100, type = 'call', american = false, payoff }) {
  const n = Math.max(1, Math.round(steps)), pay = payoff || ((s) => (type === 'call' ? Math.max(s - K, 0) : Math.max(K - s, 0)));
  if (!(sigma > 0) || !(T > 0) || !(S > 0)) return { value: pay(Math.max(S, 0)), u: 1, d: 1, p: 0.5, steps: n };
  const dt = T / n, u = Math.exp(sigma * Math.sqrt(dt)), d = 1 / u, p = clamp((Math.exp((r - q) * dt) - d) / (u - d), 0, 1), disc = Math.exp(-r * dt);
  const V = new Array(n + 1);
  for (let i = 0; i <= n; i++) V[i] = pay(S * u ** (n - i) * d ** i);
  for (let s = n - 1; s >= 0; s--) for (let i = 0; i <= s; i++) { const cont = disc * (p * V[i] + (1 - p) * V[i + 1]); V[i] = american ? Math.max(cont, pay(S * u ** (s - i) * d ** i)) : cont; }
  return { value: V[0], u, d, p, steps: n };
}
/** Minimax-regret choice. payoff[alt][scenario] (higher is better). Returns { regret, maxRegret: [], best (index) }. */
export function minimaxRegret(payoff) {
  const ns = payoff[0].length, top = zeros(ns).map((_, s) => Math.max(...payoff.map((r) => r[s]))), regret = payoff.map((r) => r.map((x, s) => top[s] - x)), maxRegret = regret.map((r) => Math.max(...r));
  return { regret, maxRegret, best: maxRegret.indexOf(Math.min(...maxRegret)) };
}

// ================================================================================================================
// 5. Optimisation
// ================================================================================================================
/**
 * Two-phase simplex method with Bland's rule. Maximise (or minimise) c·x subject to A·x (<= | >= | =) b, x >= 0.
 * { c, A, b, types: ['<=', …], maximize = true }. Returns { status: 'optimal' | 'infeasible' | 'unbounded', x, value, iterations }.
 */
export function simplex({ c, A, b, types, maximize = true, maxIter = 5000 }) {
  const m = A.length, n = c.length, eps = 1e-9;
  const rows = A.map((r, i) => { let row = c.map((_, j) => +r[j] || 0), rhs = +b[i], t = (types && types[i]) || '<='; if (rhs < 0) { row = row.map((x) => -x); rhs = -rhs; t = t === '<=' ? '>=' : t === '>=' ? '<=' : '='; } return { row, rhs, t }; });
  const nSlack = rows.filter((r) => r.t !== '=').length, nArt = rows.filter((r) => r.t !== '<=').length, N = n + nSlack + nArt, T = rows.map(() => zeros(N + 1)), basis = zeros(m);
  let s = n, a = n + nSlack, iter = 0;
  rows.forEach((r, i) => { r.row.forEach((x, j) => (T[i][j] = x)); T[i][N] = r.rhs; if (r.t === '<=') { T[i][s] = 1; basis[i] = s++; } else { if (r.t === '>=') T[i][s++] = -1; T[i][a] = 1; basis[i] = a++; } });
  const isArt = (j) => j >= n + nSlack;
  const pivot = (r, col) => { const p = T[r][col]; for (let j = 0; j <= N; j++) T[r][j] /= p; for (let i = 0; i < m; i++) if (i !== r && T[i][col] !== 0) { const f = T[i][col]; for (let j = 0; j <= N; j++) T[i][j] -= f * T[r][j]; } basis[r] = col; };
  const optimise = (cost, allow) => {
    for (; iter < maxIter; iter++) {
      let enter = -1;
      for (let j = 0; j < N && enter < 0; j++) { if (!allow(j) || basis.includes(j)) continue; let rc = cost[j]; for (let i = 0; i < m; i++) rc -= cost[basis[i]] * T[i][j]; if (rc < -eps) enter = j; }
      if (enter < 0) return 'optimal';
      let leave = -1, best = Infinity;
      for (let i = 0; i < m; i++) if (T[i][enter] > eps) { const ratio = T[i][N] / T[i][enter]; if (ratio < best - 1e-12 || (ratio <= best + 1e-12 && leave >= 0 && basis[i] < basis[leave])) { best = ratio; leave = i; } }
      if (leave < 0) return 'unbounded';
      pivot(leave, enter);
    }
    return 'iteration limit';
  };
  if (nArt) {
    optimise(zeros(N).map((_, j) => (isArt(j) ? 1 : 0)), () => true);
    if (sum(basis.map((bj, i) => (isArt(bj) ? T[i][N] : 0))) > 1e-7) return { status: 'infeasible', x: null, value: null, iterations: iter };
    for (let i = 0; i < m; i++) if (isArt(basis[i])) { const j = T[i].findIndex((x, jj) => jj < n + nSlack && Math.abs(x) > 1e-9); if (j >= 0) pivot(i, j); }
  }
  const status = optimise(zeros(N).map((_, j) => (j < n ? (maximize ? -c[j] : c[j]) : 0)), (j) => !isArt(j));
  if (status !== 'optimal') return { status, x: null, value: null, iterations: iter };
  const x = zeros(n);
  basis.forEach((bj, i) => { if (bj < n) x[bj] = Math.max(0, T[i][N]); });
  return { status, x, value: sum(x.map((v, j) => v * c[j])), iterations: iter };
}
/** Mixed-integer linear programme by depth-first branch-and-bound on the simplex relaxation. integers: indices that must be whole numbers (default all). */
export function branchAndBound({ c, A, b, types, maximize = true, integers, maxNodes = 4000 }) {
  const ints = integers || c.map((_, i) => i), sense = maximize ? 1 : -1, stack = [{ A: A.map((r) => r.slice()), b: b.slice(), types: (types || A.map(() => '<=')).slice() }];
  let best = null, nodes = 0, relaxation = null;
  while (stack.length && nodes < maxNodes) {
    const nd = stack.pop(), r = simplex({ c, A: nd.A, b: nd.b, types: nd.types, maximize });
    nodes++;
    if (nodes === 1 && r.status === 'optimal') relaxation = r.value;
    if (r.status !== 'optimal' || (best && sense * r.value <= sense * best.value + 1e-9)) continue; // infeasible or bounded out
    const k = ints.find((i) => Math.abs(r.x[i] - Math.round(r.x[i])) > 1e-6);
    if (k === undefined) { best = { x: r.x.map((v, i) => (ints.includes(i) ? Math.round(v) : v)), value: r.value }; continue; }
    const row = c.map((_, i) => (i === k ? 1 : 0)), fl = Math.floor(r.x[k]);
    stack.push({ A: [...nd.A, row], b: [...nd.b, fl], types: [...nd.types, '<='] }, { A: [...nd.A, row], b: [...nd.b, fl + 1], types: [...nd.types, '>='] });
  }
  return best ? { status: 'optimal', x: best.x, value: best.value, nodes, relaxation } : { status: 'infeasible', x: null, value: null, nodes, relaxation };
}
/** Projected steepest descent with central-difference gradients and Armijo backtracking (variables scaled by their box). */
export function gradientDescent(f, x0, { lo, hi, maxIter = 80, tol = 1e-10, h = 1e-4 } = {}) {
  const n = x0.length, span = x0.map((x, i) => (lo && hi ? hi[i] - lo[i] : Math.max(1, Math.abs(x)))), proj = (x) => x.map((v, i) => clamp(v, lo ? lo[i] : -Infinity, hi ? hi[i] : Infinity));
  let x = proj(x0), fx = f(x), evals = 1, it = 0, step = 0.25;
  for (; it < maxIter; it++) {
    const g = x.map((_, i) => { const d = h * span[i], xp = x.slice(), xm = x.slice(); xp[i] = Math.min(x[i] + d, hi ? hi[i] : Infinity); xm[i] = Math.max(x[i] - d, lo ? lo[i] : -Infinity); evals += 2; return ((f(xp) - f(xm)) / (xp[i] - xm[i] || 1)) * span[i]; }); // gradient in scaled coordinates
    const gn = Math.sqrt(sum(g.map((v) => v * v)));
    if (!(gn > 0)) break;
    let moved = false;
    for (let k = 0; k < 30; k++, step *= 0.5) {
      const xn = proj(x.map((v, i) => v - (step * span[i] * g[i]) / gn)), fn = f(xn); evals++;
      if (fn < fx - 1e-4 * step * gn * Math.sqrt(sum(xn.map((v, i) => ((v - x[i]) / span[i]) ** 2))) / Math.max(step, 1e-300)) { const df = fx - fn; x = xn; fx = fn; moved = true; step *= 2; if (df < tol * (1 + Math.abs(fx))) it = maxIter; break; }
    }
    if (!moved) break;
  }
  return { x, f: fx, evals, iterations: it };
}
/** Real-coded genetic algorithm (tournament selection, blend crossover, Gaussian mutation, elitism); minimises f inside the box. */
export function geneticAlgorithm(f, lo, hi, { pop = 30, gens = 30, seed = 3, pc = 0.9, pm = 0.2, elite = 2 } = {}) {
  const R = rng(seed);
  let P = Array.from({ length: pop }, () => { const x = lo.map((l, j) => R.uniform(l, hi[j])); return { x, f: f(x) }; }), evals = pop;
  const pick = () => { const a = P[R.int(pop)], b = P[R.int(pop)]; return a.f < b.f ? a : b; };
  for (let g = 0; g < gens; g++) {
    P.sort((a, b) => a.f - b.f);
    const next = P.slice(0, elite);
    while (next.length < pop) {
      const a = pick().x, b = pick().x, cross = R.uniform() < pc;
      const x = a.map((v, j) => { let y = v; if (cross) { const al = R.uniform(-0.25, 1.25); y = al * v + (1 - al) * b[j]; } if (R.uniform() < pm) y += R.normal(0, 0.12 * (hi[j] - lo[j]) * (1 - g / gens) + 1e-12); return clamp(y, lo[j], hi[j]); });
      next.push({ x, f: f(x) }); evals++;
    }
    P = next;
  }
  P.sort((a, b) => a.f - b.f);
  return { x: P[0].x, f: P[0].f, evals };
}
/** Particle-swarm optimisation (inertia–constriction form); minimises f inside the box. */
export function particleSwarm(f, lo, hi, { n = 20, iters = 30, seed = 5, w = 0.72, c1 = 1.49, c2 = 1.49 } = {}) {
  const R = rng(seed), d = lo.length;
  const S = Array.from({ length: n }, () => { const x = lo.map((l, j) => R.uniform(l, hi[j])), fx = f(x); return { x, v: lo.map((l, j) => R.uniform(-0.2, 0.2) * (hi[j] - l)), best: x.slice(), fb: fx }; });
  let g = S.reduce((a, s) => (s.fb < a.fb ? s : a), S[0]), gx = g.best.slice(), gf = g.fb, evals = n;
  for (let it = 0; it < iters; it++) for (const s of S) {
    for (let j = 0; j < d; j++) { s.v[j] = w * s.v[j] + c1 * R.uniform() * (s.best[j] - s.x[j]) + c2 * R.uniform() * (gx[j] - s.x[j]); s.x[j] = clamp(s.x[j] + s.v[j], lo[j], hi[j]); }
    const fx = f(s.x); evals++;
    if (fx < s.fb) { s.fb = fx; s.best = s.x.slice(); if (fx < gf) { gf = fx; gx = s.x.slice(); } }
  }
  return { x: gx, f: gf, evals };
}
const dominates = (a, b) => { let better = false; for (let i = 0; i < a.length; i++) { if (a[i] > b[i]) return false; if (a[i] < b[i]) better = true; } return better; };
/** Indices of the non-dominated rows of F (every objective is minimised). */
export const paretoFront = (F) => F.map((_, i) => i).filter((i) => !F.some((g, j) => j !== i && dominates(g, F[i])));
/** Seeded NSGA-II-style multi-objective search (non-dominated sorting + crowding distance). evalFn(x) → objectives to minimise. Returns the final non-dominated set [{ x, f }]. */
export function nsga2(evalFn, lo, hi, { pop = 40, gens = 20, seed = 9 } = {}) {
  const R = rng(seed), d = lo.length, mk = (x) => ({ x, f: evalFn(x), rank: 0, crowd: 0 });
  const sortFronts = (P) => {
    let rest = P.map((_, i) => i), rank = 0;
    while (rest.length) {
      const front = rest.filter((i) => !rest.some((j) => j !== i && dominates(P[j].f, P[i].f)));
      if (!front.length) { rest.forEach((i) => { P[i].rank = rank; P[i].crowd = 0; }); break; }
      front.forEach((i) => { P[i].rank = rank; P[i].crowd = 0; });
      for (let k = 0; k < P[0].f.length; k++) {
        const s = front.slice().sort((a, b) => P[a].f[k] - P[b].f[k]), span = P[s[s.length - 1]].f[k] - P[s[0]].f[k] || 1;
        P[s[0]].crowd = P[s[s.length - 1]].crowd = 1e9;
        for (let i = 1; i < s.length - 1; i++) P[s[i]].crowd += (P[s[i + 1]].f[k] - P[s[i - 1]].f[k]) / span;
      }
      rest = rest.filter((i) => !front.includes(i)); rank++;
    }
  };
  let P = Array.from({ length: pop }, () => mk(lo.map((l, j) => R.uniform(l, hi[j]))));
  sortFronts(P);
  const better = (a, b) => (a.rank < b.rank || (a.rank === b.rank && a.crowd > b.crowd) ? a : b), pick = () => better(P[R.int(pop)], P[R.int(pop)]);
  for (let g = 0; g < gens; g++) {
    const kids = [];
    while (kids.length < pop) { const a = pick().x, b = pick().x; kids.push(mk(a.map((v, j) => { const al = R.uniform(-0.2, 1.2); let y = al * v + (1 - al) * b[j]; if (R.uniform() < 1 / d) y += R.normal(0, 0.1 * (hi[j] - lo[j])); return clamp(y, lo[j], hi[j]); }))); }
    const all = P.concat(kids);
    sortFronts(all);
    all.sort((a, b) => a.rank - b.rank || b.crowd - a.crowd);
    P = all.slice(0, pop);
  }
  const F = P.map((p) => p.f);
  return paretoFront(F).map((i) => ({ x: P[i].x.slice(), f: P[i].f.slice() }));
}
/**
 * Optimal replacement timing by backward dynamic programming. State = asset age at the start of a year; each year either
 * keep (pay opCost(age)) or replace (pay replaceCost, then opCost(0)). { horizon, maxAge, opCost(age), replaceCost, rate, age0 }.
 * Returns { cost (present value), replaceYears: [], policy: [{ year, age, action }] }.
 */
export function replacementDP({ horizon, maxAge, opCost, replaceCost, rate = 0, age0 = 0 }) {
  const df = 1 / (1 + rate), A = Math.max(1, Math.round(maxAge)), H = Math.max(1, Math.round(horizon)), V = Array.from({ length: H + 1 }, () => zeros(A + 2)), act = Array.from({ length: H }, () => zeros(A + 2));
  for (let t = H - 1; t >= 0; t--) for (let a = 0; a <= A; a++) {
    const keep = a < A ? opCost(a) + df * V[t + 1][a + 1] : Infinity, repl = replaceCost + opCost(0) + df * V[t + 1][1];
    V[t][a] = Math.min(keep, repl); act[t][a] = repl < keep ? 1 : 0;
  }
  const policy = [], replaceYears = [];
  let a = Math.min(Math.max(0, Math.round(age0)), A);
  for (let t = 0; t < H; t++) { const rep = act[t][a] === 1; policy.push({ year: t, age: a, action: rep ? 'replace' : 'keep' }); if (rep) { replaceYears.push(t); a = 1; } else a += 1; }
  return { cost: V[0][Math.min(Math.max(0, Math.round(age0)), A)], replaceYears, policy };
}
/**
 * Two-stage stochastic programme (capacity now, recourse later) solved as its deterministic-equivalent LP:
 * minimise c·x + Σ p_s·q·y_s subject to x + y_s >= d_s. scenarios: [{ p, demand }]. Returns { x, cost, shortfall: [], status }.
 */
export function twoStage({ c, q, scenarios }) {
  const ps = sum(scenarios.map((s) => s.p)), S = scenarios.length;
  if (!(ps > 0)) throw new Error('The scenario probabilities must sum to a positive number.');
  const r = simplex({ c: [c, ...scenarios.map((s) => (s.p / ps) * q)], A: scenarios.map((_, i) => [1, ...zeros(S).map((__, j) => (i === j ? 1 : 0))]), b: scenarios.map((s) => s.demand), types: scenarios.map(() => '>='), maximize: false });
  return r.status === 'optimal' ? { status: r.status, x: r.x[0], cost: r.value, shortfall: r.x.slice(1) } : { status: r.status, x: null, cost: null, shortfall: [] };
}

// ================================================================================================================
// 6. CAPEX: parametric and bottom-up cost models
// ================================================================================================================
/** Default cost basis (2023 US$, cost index 800). ref = purchased cost in M$ at capacity cap; cost = ref·(capacity/cap)^exp; fac = installation (bare-module) factor. */
export const COST_BASIS = Object.freeze([
  { id: 'tree', item: 'Subsea tree, wellhead and controls (per well)', ref: 9, cap: 1, unit: 'well', exp: 1, fac: 1.25 },
  { id: 'manifold', item: 'Production manifold with foundation', ref: 16, cap: 4, unit: 'slots', exp: 0.6, fac: 1.3 },
  { id: 'jumper', item: 'Rigid jumper with connectors (each)', ref: 1.4, cap: 1, unit: 'each', exp: 1, fac: 1.4 },
  { id: 'plet', item: 'Pipeline end termination (each)', ref: 3, cap: 1, unit: 'each', exp: 1, fac: 1.3 },
  { id: 'umbilical', item: 'Control umbilical (electro-hydraulic)', ref: 1.1, cap: 1, unit: 'km', exp: 1, fac: 1.15 },
  { id: 'chemline', item: 'Chemical line in the umbilical (per line)', ref: 0.14, cap: 1, unit: 'km', exp: 1, fac: 1.15 },
  { id: 'slugcatcher', item: 'Slug catcher / inlet separator vessel', ref: 1.6, cap: 60, unit: 'm³', exp: 0.6, fac: 4 },
  { id: 'megregen', item: 'MEG regeneration and reclamation package', ref: 14, cap: 300, unit: 'm³/d lean MEG', exp: 0.65, fac: 2.6 },
  { id: 'cheminj', item: 'Chemical-injection skid (tanks, pumps)', ref: 1.2, cap: 10, unit: 'm³/d', exp: 0.5, fac: 2.8 },
  { id: 'pump', item: 'Booster / export pump with driver', ref: 2.4, cap: 1000, unit: 'kW', exp: 0.7, fac: 2.6 },
  { id: 'compressor', item: 'Gas compressor train with driver', ref: 8.5, cap: 5000, unit: 'kW', exp: 0.75, fac: 2.6 },
  { id: 'dehpower', item: 'Heating power unit, riser cable and feeder', ref: 7, cap: 2000, unit: 'kW', exp: 0.6, fac: 1.8 },
  { id: 'pigtrap', item: 'Pig launcher and receiver (pair)', ref: 1.1, cap: 10, unit: 'in', exp: 0.8, fac: 2.2 },
].map(Object.freeze));
const basisMap = (table) => { const B = Object.fromEntries(COST_BASIS.map((r) => [r.id, r])); for (const r of Array.isArray(table) ? table : []) { const id = String(r?.id ?? '').trim(); if (B[id] && r.ref >= 0 && r.cap > 0 && isNum(+r.exp) && r.fac > 0) B[id] = { ...B[id], ref: +r.ref, cap: +r.cap, exp: clamp(+r.exp, 0, 1.5), fac: +r.fac }; } return B; };
/** Cumulative cost multiplier of n units on a learning curve (unit n costs n^log2(rate) of the first). */
export const learningSum = (n, rate) => { let s = 0; for (let i = 1; i <= Math.round(n); i++) s += i ** Math.log2(rate); return s; };
/** Parametric cost-estimating relationship for an installed pipeline (M$): material cerCoef·(D/10 in)^cerExp per km, lay spread by day rate, mobilisation, learning curve. */
export function pipelineCER({ dIn = 10, lengthKm = 20, depth = 1000, unitNo = 1, cerCoef = 1.35, cerExp = 1.1, layFactor = 1, learnRate = 0.9, vesselRate = 350, layRate = 2.5, mobCost = 6, depthCoef = 0.2 }) {
  const material = cerCoef * (dIn / 10) ** cerExp * lengthKm, days = lengthKm / (layRate * (10 / dIn) ** 0.5), lay = (layFactor * vesselRate * (1 + (depthCoef * depth) / 1000) * days) / 1000;
  return (material + lay + mobCost) * Math.max(unitNo, 1) ** Math.log2(learnRate);
}
/**
 * Capital cost build-up. c: geometry (flowLen, riserLen, id, wt in m, depth), strategy ('none' | 'bare' | 'wet' | 'pip' | 'deh' | 'ldhi' | 'risk'),
 * insT (m), material ('cs' | 'cra'), caExtra (m extra wall), nWells, slugVol (m³), megRate, chemRate (m³/d), pumpKW, compKW, heatKW and the cost basis.
 * Returns { items: [{ group, item, basis, cost }], groups: {…}, direct, contingency, owners, total (US$), steelT, layDays, purchased, moduleCost, langCost, od }.
 */
export function capexEstimate(c) {
  const B = c.basis || basisMap(c.costBasis), items = [], esc = (c.costIndexEval / c.costIndexBase) * c.locFactor, st = c.strategy || 'wet';
  const add = (group, item, cost, basis) => { if (cost > 0) items.push({ group, item, basis: c.brief || !basis ? '' : basis(), cost: cost * esc }); }; // basis text is built only for the reported estimate
  const purchased = (id, cap) => B[id].ref * MM * (Math.max(cap, 0) / B[id].cap) ** B[id].exp;
  // --- line pipe, coatings and thermal system (per metre)
  const wt = c.wt + (c.caExtra || 0), od = c.id + 2 * wt, dIn = c.id / 0.0254, L = c.flowLen + c.riserLen, tCoat = 0.003;
  const insT = st === 'none' || st === 'bare' || st === 'pip' ? 0 : Math.max(c.insT, 0), steelM = (BASE.rhoSteel * Math.PI * (od * od - c.id * c.id)) / 4; // kg/m
  const cer = c.pipeMethod === 'cer', steelPerM = cer ? (c.cerCoef * MM * (dIn / 10) ** c.cerExp) / 1000 : (steelM / 1000) * c.steelPrice * (c.material === 'cra' ? c.craFactor : 1) + Math.PI * od * c.coatPrice + c.fabPerM;
  const oc = od + 2 * tCoat, insVolM = (Math.PI * ((oc + 2 * insT) ** 2 - oc * oc)) / 4, insPerM = insVolM * c.insPrice;
  let pipPerM = 0, steelT = (steelM * L) / 1000;
  if (st === 'pip') { const id2 = od + 0.06, wt2 = Math.max(0.0127, 0.045 * id2), m2 = (BASE.rhoSteel * Math.PI * ((id2 + 2 * wt2) ** 2 - id2 * id2)) / 4; pipPerM = (m2 / 1000) * c.steelPrice + Math.PI * (id2 + 2 * wt2) * c.coatPrice + c.pipPremium; steelT += (m2 * c.flowLen) / 1000; }
  add('Pipeline', cer ? 'Line pipe, coating and fabrication (parametric relationship)' : c.material === 'cra' ? 'Line pipe (CRA-clad), coating and welding' : 'Line pipe, anti-corrosion coating and welding', steelPerM * c.flowLen, () => (cer ? `${c.cerCoef} M$/km × (D/10 in)^${c.cerExp}` : `${((steelM * c.flowLen) / 1000).toFixed(0)} t steel`));
  add('Pipeline', 'Wet insulation', insPerM * c.flowLen, () => `${(insVolM * c.flowLen).toFixed(0)} m³ at ${(insT * 1000).toFixed(0)} mm`);
  add('Pipeline', 'Pipe-in-pipe carrier pipe, annulus insulation, bulkheads', pipPerM * c.flowLen, () => 'outer pipe steel + premium per metre');
  add('Pipeline', 'Direct electrical heating cable and anodes', st === 'deh' ? c.dehCable * c.flowLen : 0, () => `${c.dehCable} $/m`);
  add('Riser', 'Steel catenary riser (pipe, insulation, strakes, fatigue-class welds)', (steelPerM + insPerM + (st === 'pip' ? 0.5 * pipPerM : 0)) * c.riserFactor * c.riserLen, () => `flowline unit cost × ${c.riserFactor}`);
  // --- installation by vessel day rate
  const trees = Math.max(1, Math.round(c.nWells)), jumpers = trees + 2, umbKm = (1.05 * L) / 1000, nLines = { none: 0, bare: 2, wet: 1, pip: 1, deh: 1, ldhi: 2, risk: 0 }[st] ?? 1;
  const layKmD = c.layRate * (10 / dIn) ** 0.5 * (st === 'pip' ? 0.5 : st === 'deh' ? 0.8 : 1) * (1 - Math.min(0.3, insT * 1.5)), depthF = 1 + (c.depthCoef * c.depth) / 1000, dayRate = c.vesselRate * 1000 * depthF * c.layFactor;
  const dPipe = L / 1000 / layKmD, dRiser = c.riserLen > 0 ? 8 : 0, dSub = 1.5 * jumpers + 2 * trees + 3, dUmb = umbKm / 4, layDays = dPipe + dRiser + dSub + dUmb;
  add('Installation', 'Pipelay spread', dPipe * dayRate, () => `${dPipe.toFixed(1)} d at ${layKmD.toFixed(2)} km/d`);
  add('Installation', 'Riser pull-in and hang-off', dRiser * dayRate, () => `${dRiser} d`);
  add('Installation', 'Subsea structures, trees and tie-ins', dSub * dayRate, () => `${dSub.toFixed(1)} d`);
  add('Installation', 'Umbilical lay', dUmb * dayRate * 0.7, () => `${dUmb.toFixed(1)} d`);
  add('Installation', 'Mobilisation, demobilisation, survey and pre-commissioning', c.mobCost * MM * c.layFactor + 0.04 * layDays * dayRate, () => 'lump sum + 4 % of spread cost');
  // --- subsea equipment and wells
  const lc = learningSum(trees, c.learnRate);
  add('Subsea', B.tree.item, B.tree.ref * MM * lc * B.tree.fac, () => `${trees} units on a ${(c.learnRate * 100).toFixed(0)} % learning curve`);
  add('Subsea', B.manifold.item, purchased('manifold', Math.max(trees, 2)) * B.manifold.fac, () => `${Math.max(trees, 2)} slots`);
  add('Subsea', 'Jumpers and connectors', jumpers * B.jumper.ref * MM * B.jumper.fac, () => `${jumpers} jumpers`);
  add('Subsea', 'Pipeline end terminations', 2 * B.plet.ref * MM * B.plet.fac, () => '2 units');
  add('Subsea', B.umbilical.item, purchased('umbilical', umbKm) * B.umbilical.fac, () => `${umbKm.toFixed(1)} km`);
  add('Subsea', 'Chemical lines in the umbilical', nLines * purchased('chemline', umbKm) * B.chemline.fac, () => `${nLines} line(s)`);
  add('Wells', 'Drilling and completion', c.wellCost * MM * lc, () => `${trees} well(s), learning curve`);
  // --- topsides by capacity scaling; bare-module factors or one Lang factor
  const eq = [['slugcatcher', c.slugVol], ['megregen', st === 'bare' ? c.megRate : 0], ['cheminj', c.chemRate], ['pump', c.pumpKW], ['compressor', c.compKW], ['dehpower', st === 'deh' ? c.heatKW : 0], ['pigtrap', dIn]].filter((e) => e[1] > 0).map(([id, cap]) => ({ id, cap, E: purchased(id, cap) }));
  const pE = sum(eq.map((e) => e.E)), moduleCost = sum(eq.map((e) => e.E * B[e.id].fac)), langCost = pE * c.langFactor, useLang = c.costMethod === 'lang';
  for (const e of eq) add('Topsides', B[e.id].item, e.E * (useLang ? c.langFactor : B[e.id].fac), () => `${+e.cap.toPrecision(3)} ${B[e.id].unit}, exponent ${B[e.id].exp}, factor ${useLang ? c.langFactor : B[e.id].fac}`);
  const direct = sum(items.map((i) => i.cost)), contingency = direct * c.contingency, owners = direct * c.owners;
  items.push({ group: 'Indirects', item: 'Contingency', basis: `${(c.contingency * 100).toFixed(0)} % of direct cost`, cost: contingency }, { group: 'Indirects', item: "Owner's costs (project team, insurance, studies)", basis: `${(c.owners * 100).toFixed(0)} % of direct cost`, cost: owners });
  const groups = {};
  for (const it of items) groups[it.group] = (groups[it.group] || 0) + it.cost;
  return { items, groups, direct, contingency, owners, total: direct + contingency + owners, steelT, layDays, purchased: pE * esc, moduleCost: moduleCost * esc, langCost: langCost * esc, od, escalation: esc, pipeInstalled: (groups.Pipeline || 0) + (groups.Riser || 0) + (groups.Installation || 0) };
}

// ================================================================================================================
// 7. Fiscal cash-flow engine
// ================================================================================================================
/**
 * Project cash flow in money of the day. p (US$, real terms of year 0 unless stated):
 *  phase[] (CAPEX fractions by construction year), capex, capexSunk, residual, life, lifeCut, oil[] (bbl/y potential), gas[] (MMBtu/y potential sold),
 *  water[] (m³/y), uptime, deferFrac, oilPrice, gasPrice, infl, costEsc, priceEsc (real escalation), discount (nominal), mid,
 *  opexFixed ($/y), opexVarBoe ($/boe), waterCost ($/m³), opexDown ($/y, scales with downtime), opexBlock ($/y expected blockage cost),
 *  carbonT (t/y), carbonPrice, carbonEsc, includeRisk, consequence ($), haz[] (annual failure probability), royalty, taxRate, regime ('tax' | 'psc'),
 *  costOilCap, profitSplit, deprFrac[], wcDays, abandon, abandonProvision, gearing, loanRate, loanTenor.
 * m: multipliers { price, prod, capex, opex, downtime, failFreq, repair, path[] (price path), events[] (uniforms for failure sampling) }.
 * detail = false returns only the NPV (fast path for sampling). Sign convention: receipts positive, payments negative in `fcf`.
 */
export function cashflow(p, m = {}, detail = true) {
  const mPrice = m.price ?? 1, mProd = m.prod ?? 1, mCapex = m.capex ?? 1, mOpex = m.opex ?? 1, mDown = m.downtime ?? 1, mFail = m.failFreq ?? 1, mRep = m.repair ?? 1;
  const nCon = p.phase.length, life = Math.max(1, Math.min(p.life, p.lifeCut ?? p.life)), K = nCon + life + 1, r = p.discount, psc = p.regime === 'psc';
  const down = clamp((1 - p.uptime) * mDown, 0, 0.95), avail = 1 - down, cg = (1 + p.infl) * (1 + p.costEsc), pg = (1 + p.infl) * (1 + p.priceEsc), kg = (1 + p.infl) * (1 + p.carbonEsc);
  const abNom = p.abandon * cg ** (nCon + life), accr = p.abandonProvision ? abNom / life : 0;
  let deprBase = (p.capexSunk || 0) + (p.residual || 0);
  for (let k = 0; k < nCon; k++) deprBase += p.capex * mCapex * p.phase[k] * cg ** k;
  let book = deprBase, pool = 0, poolL = 0, rec = 0, wcPrev = 0, defOil = 0, defGas = 0, npvSum = 0, cum = 0, cumD = 0, debt = 0, pay = 0;
  const D = detail ? Object.fromEntries(['year', 'oil', 'gas', 'boe', 'potBoe', 'lostBoe', 'defBoe', 'water', 'revenue', 'royalty', 'govShare', 'opex', 'carbon', 'risk', 'ocf', 'depreciation', 'taxable', 'tax', 'atcf', 'capex', 'dwc', 'abandon', 'fcf', 'cum', 'df', 'dcf', 'cumDcf', 'real', 'interest', 'debtService', 'equity', 'price'].map((k) => [k, zeros(K)])) : null;
  for (let k = 0; k < K; k++) {
    const j = k - nCon, cgk = cg ** k, df = (1 + r) ** -(p.mid && k > 0 ? k - 0.5 : k), capex = k < nCon ? p.capex * mCapex * p.phase[k] * cgk : 0;
    let rev = 0, roy = 0, gov = 0, opex = 0, carbon = 0, risk = 0, ocf = 0, dep = 0, taxable = 0, tax = 0, taxL = 0, dwc = 0, ab = 0, oil = 0, gas = 0, boe = 0, pot = 0, lost = 0, defd = 0, water = 0, interest = 0, service = 0, draw = 0, priceF = 0;
    if (psc && capex > 0) rec += capex;
    if (p.gearing > 0 && k < nCon) { draw = p.gearing * capex; debt = debt * (1 + p.loanRate) + draw; }
    if (j >= 0 && j < life) {
      const pOil = p.oil[j] * mProd, pGas = p.gas[j] * mProd;
      oil = pOil * avail; gas = pGas * avail; defOil += pOil * down * p.deferFrac; defGas += pGas * down * p.deferFrac;
      pot = pOil + pGas / MMBTU_PER_BOE; defd = pot * down * p.deferFrac; lost = pot * down * (1 - p.deferFrac);
      if (j === life - 1) { oil += defOil; gas += defGas; defd -= (defOil + defGas / MMBTU_PER_BOE); } // deferred barrels come back in the last year
      boe = oil + gas / MMBTU_PER_BOE; water = p.water[j] * mProd * avail;
      priceF = pg ** k * mPrice * (m.path ? m.path[k] : 1);
      rev = (oil * p.oilPrice + gas * p.gasPrice) * priceF; roy = p.royalty * rev;
      opex = ((p.opexFixed + p.opexVarBoe * boe + p.waterCost * water) * mOpex + p.opexDown * mDown * mOpex + p.opexBlock * mFail * mRep) * cgk;
      carbon = p.carbonT * p.carbonPrice * kg ** k;
      if (p.includeRisk) risk = (m.events ? (m.events[j] < p.haz[j] * mFail ? 1 : 0) : p.haz[j] * mFail) * p.consequence * mRep * cgk;
      dep = j === life - 1 ? book : Math.min(book, deprBase * p.deprFrac[j]); book -= dep;
      if (j === 0 && debt > 0) pay = debt * capitalRecovery(p.loanRate, Math.max(1, Math.min(p.loanTenor, life)));
      if (debt > 1e-6) { interest = debt * p.loanRate; service = Math.min(pay, debt + interest); debt -= service - interest; }
      if (psc) {
        const net = rev - roy; rec += opex + carbon + risk;
        const costOil = Math.min(rec, p.costOilCap * net), profit = net - costOil; rec -= costOil;
        gov = (1 - p.profitSplit) * profit; ocf = rev - roy - gov - opex - carbon - risk;
        taxable = p.profitSplit * profit; tax = taxL = p.taxRate * Math.max(taxable, 0);
      } else {
        ocf = rev - roy - opex - carbon - risk; taxable = ocf - dep - accr;
        if (taxable < 0) pool -= taxable; else { const use = Math.min(pool, taxable); pool -= use; tax = p.taxRate * (taxable - use); }
        const tl = taxable - interest; // levered tax base for the equity view
        if (tl < 0) poolL -= tl; else { const use = Math.min(poolL, tl); poolL -= use; taxL = p.taxRate * (tl - use); }
      }
      const wc = (p.wcDays / 365) * rev; dwc = wc - wcPrev; wcPrev = wc;
    } else if (k === K - 1) { ab = abNom; dwc = -wcPrev; wcPrev = 0; }
    const fcf = ocf - tax - capex - dwc - ab;
    npvSum += fcf * df;
    if (detail) {
      cum += fcf; cumD += fcf * df;
      const set = (key, val) => (D[key][k] = val);
      set('year', k); set('oil', oil); set('gas', gas); set('boe', boe); set('potBoe', pot); set('lostBoe', lost); set('defBoe', defd); set('water', water); set('revenue', rev); set('royalty', roy); set('govShare', gov); set('opex', opex); set('carbon', carbon); set('risk', risk);
      set('ocf', ocf); set('depreciation', dep); set('taxable', taxable); set('tax', tax); set('atcf', ocf - tax); set('capex', capex); set('dwc', dwc); set('abandon', ab); set('fcf', fcf); set('cum', cum); set('df', df); set('dcf', fcf * df); set('cumDcf', cumD);
      set('real', fcf / (1 + p.infl) ** (p.mid && k > 0 ? k - 0.5 : k)); set('interest', interest); set('debtService', service); set('equity', ocf - taxL - capex - dwc - ab + draw - service); set('price', p.oilPrice * priceF);
    }
  }
  if (!detail) return npvSum;
  return { ...D, K, nCon, life, npv: npvSum, deprBase, capexNominal: sum(D.capex), lossPoolEnd: pool, unrecoveredCost: rec, uptime: avail };
}
/**
 * Investment metrics of a detailed cash flow: npv, irr, mirr, payback, discountedPayback, pi, roi, roce, eav, eac, utc ($/boe),
 * liftingCost ($/boe), pvCapex, pvOpex, pvBoe, realNpv (deflated flows at the Fisher real rate), governmentTake, economicLimit (production year or null).
 */
export function investmentMetrics(cf, p) {
  const r = p.discount, opt = { mid: p.mid }, I = irr(cf.fcf), rr = realRate(r, p.infl), pvOf = (a) => npv(r, a, opt), life = cf.life;
  const pvCapex = pvOf(cf.capex), pvOpex = pvOf(cf.opex.map((x, k) => x + cf.carbon[k] + cf.risk[k])), pvAb = pvOf(cf.abandon), pvBoe = npv(rr, cf.boe, opt), totCapex = sum(cf.capex), totBoe = sum(cf.boe);
  const profit = cf.ocf.map((x, k) => x - cf.tax[k] - cf.depreciation[k]);
  let bookV = cf.deprBase, capEmp = 0;
  for (let k = cf.nCon; k < cf.nCon + life; k++) { capEmp += bookV - 0.5 * cf.depreciation[k]; bookV -= cf.depreciation[k]; }
  const take = sum(cf.royalty) + sum(cf.tax) + sum(cf.govShare), preTake = sum(cf.fcf) + take;
  let limit = null;
  for (let k = cf.nCon; k < cf.nCon + life && limit === null; k++) if (cf.revenue[k] - cf.royalty[k] - cf.govShare[k] - cf.opex[k] - cf.carbon[k] < 0) limit = k - cf.nCon;
  return {
    npv: cf.npv, irr: I.irr, irrMultiple: I.multiple, irrRoots: I.roots, mirr: mirr(cf.fcf, r, p.reinvest ?? r), payback: payback(cf.fcf), discountedPayback: payback(cf.fcf, r, opt),
    pi: pvCapex > 0 ? 1 + cf.npv / pvCapex : null, roi: totCapex > 0 ? sum(cf.fcf) / totCapex : null, roce: capEmp > 0 ? sum(profit) / capEmp : null,
    eav: equivalentAnnual(cf.npv, r, cf.K - 1), eac: equivalentAnnual(pvCapex + pvOpex + pvAb, r, cf.K - 1), utc: pvBoe > 0 ? (pvCapex + pvOpex + pvAb) / pvBoe : null, liftingCost: totBoe > 0 ? sum(cf.opex) / totBoe : null,
    pvCapex, pvOpex, pvAb, pvBoe, realNpv: npv(rr, cf.real, opt), realRate: rr, governmentTake: preTake > 0 ? take / preTake : null, economicLimit: limit,
    equityIrr: p.gearing > 0 ? irr(cf.equity).irr : null, equityNpv: p.gearing > 0 ? npv(r, cf.equity, opt) : null,
  };
}
/** Root of the NPV with respect to one scalar driver: returns x in [lo, hi] with npvOf(x) = 0, or null when the NPV does not change sign there. */
export function breakeven(npvOf, lo, hi) {
  const a = npvOf(lo), b = npvOf(hi);
  if (!Number.isFinite(a) || !Number.isFinite(b) || a * b > 0) return null;
  return a === 0 ? lo : b === 0 ? hi : brent(npvOf, lo, hi, 1e-12);
}
/** Weibull wear-out hazard added to a constant base rate: h(t) = h0 + (β/η)(t/η)^(β−1), with η set so that the wear-out failure probability reaches pEnd at the remaining life. */
export function hazardRate(t, { h0 = 0, beta = 3, remLife = 25, pEnd = 0.1 }) {
  const eta = remLife / (-Math.log(1 - clamp(pEnd, 1e-6, 0.999))) ** (1 / beta);
  return h0 + (t > 0 ? (beta / eta) * (t / eta) ** (beta - 1) : beta === 1 ? 1 / eta : 0);
}
/**
 * Lifecycle cost of an inspection interval (risk-based inspection). The wear-out age is reduced to (1 − PoD)·age at every
 * inspection (Kijima virtual age), so more frequent or better inspections cut the expected failure cost.
 * { interval (y), life, inspCost, consequence, h0, beta, remLife, pEnd, pod, rate, age0 }. Returns { total, inspection, failure, pFail (life sum of annual probabilities), hMax (largest annual failure probability) }.
 */
export function inspectionCost({ interval, life, inspCost, consequence, h0, beta, remLife, pEnd, pod, rate, age0 = 0 }) {
  let age = age0, insp = 0, fail = 0, pf = 0, hMax = 0, next = interval;
  for (let y = 0; y < life; y++) {
    const h = Math.min(1, hazardRate(age + 0.5, { h0, beta, remLife, pEnd })), df = (1 + rate) ** -(y + 1);
    fail += h * consequence * df; pf += h; hMax = Math.max(hMax, h); age += 1;
    while (next <= y + 1 + 1e-9 && next < life - 1e-9) { insp += inspCost * (1 + rate) ** -next; age = y + 1 - next + (1 - pod) * (age - (y + 1 - next)); next += interval; }
  }
  return { total: insp + fail, inspection: insp, failure: fail, pFail: pf, hMax };
}

// ================================================================================================================
// 8. Techno-economic model of the case: physics surrogate, flow-assurance strategies, project assembly
// ================================================================================================================
export const STRATEGIES = Object.freeze({
  bare: { name: 'Bare pipe + continuous MEG', operability: 6 },
  wet: { name: 'Wet insulation + methanol at shutdown', operability: 8 },
  pip: { name: 'Pipe-in-pipe', operability: 9 },
  deh: { name: 'Direct electrical heating', operability: 7 },
  ldhi: { name: 'Low-dosage hydrate inhibitor', operability: 6 },
  risk: { name: 'Accept controlled risk', operability: 4 },
});
const EF = Object.freeze({ gasTurbine: 0.6, fuelGJ: 0.0561, flareSm3: 2.3e-3, MeOH: 0.7, MEG: 1.6, LDHI: 3, vesselDay: 96, steel: 1.9 }); // kgCO₂/kWh, t/GJ, t/Sm³, t/t chemicals, t per vessel day, t/t steel
const rd = (x, d = 2) => (x === null || x === undefined || !Number.isFinite(+x) ? '—' : +(+x).toFixed(d));
const mUSD = (x, d = 1) => rd(x / MM, d);
const lin = (xs, ys, x) => { // linear interpolation with linear extrapolation on an ascending table
  const n = xs.length;
  if (n === 1) return ys[0];
  let i = 0;
  while (i < n - 2 && x > xs[i + 1]) i++;
  return ys[i] + ((ys[i + 1] - ys[i]) * (x - xs[i])) / (xs[i + 1] - xs[i] || 1);
};
const argBest = (a, sign = 1) => a.reduce((b, x, i) => (sign * x > sign * a[b] ? i : b), 0);

/** Kernel studies of the case line: insulation, diameter and rate sweeps condensed into fast thermal and hydraulic surrogates. */
async function physics(q, ctx, prog, tick) {
  const fm = fluidModel(ctx), line = caseLine(ctx), id = q.idMm / 1000, wt = q.wtMm / 1000, n = q.nCells, warn = [], Ta = +line.tSeabed, tIn = +line.tIn;
  const uOf = (t, d = id) => uValue({ id: d, wt: (wt * d) / id, layers: [{ name: 'Anti-corrosion coating', t: 0.003, k: 0.3 }, ...(t > 0 ? [{ name: 'Insulation', t, k: q.kIns }] : [])], hIn: 1500, hOut: 500 }).U;
  const id2 = id + 2 * wt + 0.06, uPip = uValue({ id, wt, layers: [{ name: 'Annulus insulation', t: 0.03, k: 0.035 }, { name: 'Carrier pipe', t: Math.max(0.0127, 0.045 * id2), k: 45 }, { name: 'Coating', t: 0.003, k: 0.3 }], hIn: 1500, hOut: 500 }).U;
  let calls = 0;
  const solve = (opt) => { try { calls++; return steadyCase(ctx, { id, n, ...opt }); } catch (e) { return { failed: e.message }; } };
  const pick = (r) => ({ pIn: r.pIn, tArr: r.tOut, sub: Math.max(...r.subcooling), tHyd: mean(r.tHyd), dpFric: r.dpFric, dpGrav: r.dpGrav });
  // --- insulation sweep (always contains the selected thickness)
  const tSel = q.insMm, grid = linspace(0, q.tMaxMm, q.nThick), near = argBest(grid.map((t) => -Math.abs(t - tSel)));
  if (tSel <= q.tMaxMm * 1.2) grid[near] = tSel;
  const ts = [...new Set(grid.map((t) => +t.toFixed(3)))].sort((a, b) => a - b), sweep = [];
  let base = null;
  for (const [i, t] of ts.entries()) {
    const U = uOf(t / 1000), r = solve({ uValue: U });
    if (r.failed) warn.push(`No steady flow solution at ${t} mm of insulation (${r.failed})`); else { sweep.push({ tMm: t, U, ...pick(r) }); if (t === tSel || !base) base = r; }
    prog(0.02 + (0.2 * (i + 1)) / ts.length, 'Insulation sweep on the flow kernel'); await tick();
  }
  if (sweep.length < 2) throw new Error(`The flow kernel found no steady solution for this case, so the techno-economic studies cannot run: ${warn[0] || 'check the fluid rate and line data.'}`);
  const U0 = uOf(tSel / 1000), sw = sweep.slice().sort((a, b) => a.U - b.U), Us = sw.map((s) => s.U), th0 = (tArr) => Math.max((tArr - Ta) / Math.max(tIn - Ta, 1), 1e-3), lnTh = sw.map((s) => Math.log(th0(s.tArr))), tHm = mean(sw.map((s) => s.tHyd)), lnCold = sw.map((s) => Math.log(th0(tHm - s.sub))); // coldest-point temperature = hydrate temperature − subcooling
  if (!sweep.some((s) => s.tMm === tSel)) { const r = solve({ uValue: U0 }); if (!r.failed) base = r; }
  // --- rate sweep at the selected insulation
  const rates = [{ r: 1, ...pick(base) }];
  for (const r of [0.6, 1.4].slice(0, q.nRates)) { const s = solve({ uValue: U0, mScale: r }); if (!s.failed) rates.push({ r, ...pick(s) }); await tick(); }
  rates.sort((a, b) => a.r - b.r);
  const b0 = rates.find((x) => x.r === 1), others = rates.filter((x) => x.r !== 1), l0 = Math.log(th0(b0.tArr));
  const fitExp = (f, dflt) => { const e = others.map(f).filter((x) => Number.isFinite(x) && x > 0.2 && x < 3); return e.length ? mean(e) : dflt; };
  const mExp = fitExp((o) => Math.log(Math.log(th0(o.tArr)) / l0) / Math.log(1 / o.r), 1), nf = fitExp((o) => Math.log(o.dpFric / b0.dpFric) / Math.log(o.r), 1.8);
  const thermal = (U, r = 1) => { const e = 1 / r ** mExp, at = (ys) => Ta + (tIn - Ta) * Math.exp(Math.min(lin(Us, ys, U), -1e-6) * e); return { tArr: at(lnTh), sub: tHm - at(lnCold) }; };
  const pOut = base.pOut, gravOf = (r) => interp1(rates.map((x) => x.r), rates.map((x) => x.dpGrav), r), pInOf = (r, fr = b0.dpFric, gr = null) => pOut + (gr ?? gravOf(r)) + fr * r ** nf;
  prog(0.3, 'Diameter sweep on the flow kernel');
  // --- state near the riser base: thermal mass, slug unit, erosional ratio
  const iRb = Math.max(1, Math.min(base.x.length - 2, base.x.findIndex((x) => x >= line.riserBaseX) - 1)), pr = fm.at(base.P[iRb], base.T[iRb]);
  const erosOf = (r) => { const cSI = q.erosC * 1.22; let m = 0; for (let i = 0; i < r.vm.length; i++) m = Math.max(m, r.vm[i] / (cSI / Math.sqrt(Math.max(r.rhoM[i], 1)))); return m; };
  const slugOf = (r, d) => { try { const p = fm.at(r.P[iRb], r.T[iRb]), u = slugUnit({ vsl: r.vsl[iRb], vsg: r.vsg[iRb], rhoL: r.rhoL[iRb], rhoG: r.rhoG[iRb], muL: r.muL[iRb], muG: p.muG, D: d, theta: 0 }); const v = (u.volume * u.lengthMax) / Math.max(u.length, 1e-9) + r.qL[r.qL.length - 1] * 180; return Number.isFinite(v) && v > 0 ? v : null; } catch { return null; } };
  const slugRaw0 = slugOf(base, id);
  // --- diameter sweep
  const mult = linspace(Math.log(0.7), Math.log(1.4), q.nDiam).map(Math.exp);
  if (q.nDiam > 1) mult[argBest(mult.map((m) => -Math.abs(m - 1)))] = 1; else mult[0] = 1;
  const diam = [];
  for (const m of mult) {
    const d = id * m, U = uOf(tSel / 1000, d), r = m === 1 ? base : solve({ id: d, uValue: U });
    if (r.failed) { diam.push({ d, mult: m, failed: r.failed }); continue; }
    const sv = slugOf(r, d);
    diam.push({ d, mult: m, U, ...pick(r), eros: erosOf(r), slugVol: sv && slugRaw0 ? (q.slugVol * sv) / slugRaw0 : q.slugVol * m * m, vMax: Math.max(...r.vm) });
    await tick();
  }
  // --- cooldown: lumped thermal mass over the heat-loss conductance at the cold end of the flowline
  const A = (Math.PI * id * id) / 4, od = id + 2 * wt, hold = base.holdup[iRb], cFluid = A * (hold * pr.rhoL * pr.cpL + (1 - hold) * pr.rhoG * pr.cpG), cSteel = (BASE.rhoSteel * 470 * Math.PI * (od * od - id * id)) / 4;
  const tHydShut = fm.hydrateT(mean(base.P)), cdModel = (U, tArr, insT, pip) => {
    const cIns = pip ? 0.6 * cSteel * 1.6 : 0.5 * 900 * 1700 * (Math.PI * ((od + 2 * insT) ** 2 - od * od)) / 4, tau = (cFluid + cSteel + cIns) / (U * Math.PI * id);
    if (tArr <= tHydShut) return 0;
    return tHydShut <= Ta + 0.05 ? 1e4 : Math.min(1e4, (tau * Math.log((tArr - Ta) / (tHydShut - Ta))) / 3600);
  };
  const cd0 = cdModel(U0, b0.tArr, tSel / 1000, false), kCd = q.cooldownBase > 0 && cd0 > 0.1 ? q.cooldownBase / cd0 : 1;
  const wcIn = mean(base.wcut), gasInvStd = ((base.volume - base.liquidInventory) * mean(base.P)) / 1.01325 * (288.15 / (mean(base.T) + 273.15)) / 0.85;
  return {
    fm, line, id, wt, od, Ta, tIn, U0, uPip, uOf, sweep, rates, diam, thermal, pInOf, nf, mExp, calls, warn,
    cooldown: (U, tArr, insT, pip) => Math.min(1e4, kCd * cdModel(U, tArr, insT, pip)), kCd,
    base: { pIn: base.pIn, pOut, tArr: base.tOut, dpFric: base.dpFric, dpGrav: base.dpGrav, sub: Math.max(...base.subcooling), eros: q.erosIn > 0 ? q.erosIn : erosOf(base), erosKernel: erosOf(base), volume: base.volume, liqInv: base.liquidInventory, waterInv: base.liquidInventory * wcIn, gasInvStd, holdup: hold, regime: String(base.regime[iRb]), heatLossKW: base.heatLoss / 1000, tHydShut, mdot: base.mdot },
    sal: fm.aq.S, dep: fm.depression ?? 0, shut: null,
  };
}

/** Annual potential volumes for a rate multiplier r (same recoverable volume produced faster or slower): { oil (bbl/y), gas (MMBtu/y sold), gasSm3, water (m³/y), boe }. */
function profileOf(q, r = 1) {
  const b = q.declineType === 'exp' ? 0 : q.declineType === 'har' ? 1 : q.bHyp, days = productionProfile({ q0: 365 * r, plateau: q.plateau / r, Di: (q.Di / 100) * r, b, life: q.life });
  const oil = days.map((d) => d * q.qOil * BBL_PER_M3), gasSm3 = days.map((d) => d * q.qGas), gas = gasSm3.map((g) => (g * (q.gasSalesFrac / 100) * q.gasHV) / 1000 / GJ_PER_MMBTU);
  const wcAt = (t) => clamp((q.wc0 + ((q.wcEnd - q.wc0) * (1 - Math.exp((-3 * t) / q.life))) / (1 - Math.exp(-3))) / 100, 0, 0.98), water = days.map((d, j) => { const w = wcAt(j + 0.5); return (d * q.qOil * w) / (1 - w); });
  return { days, oil, gas, gasSm3, water, boe: oil.map((o, j) => o + gas[j] / MMBTU_PER_BOE), wcAt };
}

/** Shared economic context: discounting, production, unit values of downtime and blockage, energy prices, CAPEX arguments. */
function context(q, S) {
  const life = q.life, nCon = q.phase.length, infl = q.inflation / 100, r = q.discount / 100, cg = (1 + infl) * (1 + q.costEsc / 100), rr = realRate(r, infl), prof = profileOf(q, 1);
  let pvf = 0;
  for (let j = 0; j < life; j++) pvf += cg ** (nCon + j) * discountFactor(r, nCon + j, { mid: q.mid });
  const up = q.uptime / 100, revPot = prof.oil.map((o, j) => o * q.oilPrice + prof.gas[j] * q.gasPrice), marginH = (mean(revPot) * (1 - q.royalty / 100)) / 8760, lossFrac = 1 - (q.deferFrac / 100) * (1 + rr) ** -(life / 2);
  const gasPower = q.powerSource === 'gas', enPrice = gasPower ? (q.gasPrice / GJ_PER_MMBTU) * 0.011 : q.elecPrice, enCarbon = gasPower ? EF.gasTurbine : q.gridCarbon, heatPrice = q.gasPrice / GJ_PER_MMBTU / 0.85;
  const capexArgs = {
    basis: basisMap(q.costBasis), costIndexEval: q.costIndexEval, costIndexBase: q.costIndexBase, locFactor: q.locFactor, flowLen: Math.max(q.lineLen - q.riserLen, 0), riserLen: q.riserLen, id: S.id, wt: S.wt, depth: q.depth, material: 'cs', caExtra: 0,
    nWells: q.nWells, slugVol: q.slugVol, pumpKW: q.pumpKW, compKW: q.compKW, steelPrice: q.steelPrice, coatPrice: q.coatPrice, fabPerM: q.fabPerM, insPrice: q.insPrice, pipPremium: q.pipPremium, dehCable: q.dehCable, riserFactor: q.riserFactor,
    vesselRate: q.vesselRate, layRate: q.layRate, mobCost: q.mobCost, layFactor: q.layFactor, depthCoef: q.depthCoef, wellCost: q.wellCost, learnRate: q.learnRate / 100, langFactor: q.langFactor, costMethod: q.costMethod, pipeMethod: q.pipeMethod,
    cerCoef: q.cerCoef, cerExp: q.cerExp, contingency: q.contingency / 100, owners: q.owners / 100, craFactor: q.craFactor,
  };
  const qLiq0 = q.qOil + (q.qOil * (q.wc0 / 100)) / (1 - Math.min(q.wc0, 98) / 100), pAvail = (x) => S.base.pIn + q.chokeDp + ((1 - x) * qLiq0) / q.pi;
  const X = { life, nCon, infl, r, rr, cg, pvf, up, prof, marginH, lossFrac, enPrice, enCarbon, heatPrice, capexArgs, qLiq0, pAvail, waterMean: mean(prof.water) / 365, waterPeak: Math.max(...prof.water) / 365, wcMean: prof.wcAt(life / 2) };
  const h0 = q.ealOverride > 0 && q.consequence > 0 ? q.ealOverride / q.consequence : q.pof;
  X.h0 = h0; X.haz = prof.oil.map((_, j) => Math.min(1, hazardRate(j + 0.5, { h0, beta: q.weibullBeta, remLife: q.remLife, pEnd: q.pEnd / 100 })));
  X.deprFrac = q.deprMethod === 'uop' ? null : depreciation(q.deprMethod, 1, q.deprLife, { years: life, rate: q.dbRate / 100 });
  X.blockCost = q.remedDays * q.spreadRate * 1000 + (q.blockDays + q.vesselWait) * 24 * marginH * lossFrac; // remediation spread + outage, including the wait for a vessel
  X.capexNone = capexEstimate({ ...capexArgs, brief: true, strategy: 'none', insT: 0, megRate: 0, chemRate: 0, heatKW: 0 }).total;
  X.rMin = q.qMinFrac; X.rMax = Math.max(q.capacityFrac, q.qMinFrac + 0.05);
  // operating constraints g(r) <= 0: erosion, deliverability, turndown (slugging), capacity
  X.constraints = (x) => [{ name: 'erosional velocity', g: S.base.eros * x - 1 }, { name: 'deliverability (inflow and back-pressure)', g: (S.pInOf(x) - pAvail(x)) / Math.max(S.base.pIn, 1) }, { name: 'minimum stable rate (slugging / turndown)', g: X.rMin - x }, { name: 'facility capacity', g: x - X.rMax }, { name: 'design pressure', g: (S.pInOf(x) - q.mawp) / q.mawp }];
  return X;
}

/**
 * One flow-assurance strategy on the case line: thermal state from the kernel surrogate, inhibitor demand from the hydrate
 * depression correlation, cooldown from thermal mass over U, long-shutdown frequency, residual blockage frequency, pigging.
 * over: { insT (m), rate, treat (fraction of long shutdowns preserved), thermal: { tArr, sub }, inhibRate, heatKW }.
 */
function strategyModel(key, q, S, X, over = {}) {
  const st = STRATEGIES[key], pip = key === 'pip', insT = key === 'bare' || pip ? 0 : over.insT ?? q.insMm / 1000, r = over.rate ?? 1, treat = clamp(over.treat ?? 1, 0, 1), up = X.up;
  const U = pip ? S.uPip : S.uOf(insT), th = over.thermal || S.thermal(U, r), need = th.sub + q.hydMargin, water = X.waterMean * r, waterPeak = X.waterPeak * r;
  const wtFor = (dT, inh) => (dT > 0 ? inhibitorFor(S.dep + dT, inh, S.sal) : 0), cdH = S.cooldown(U, th.tArr, insT, pip), pLong = cdH >= 1e4 ? 0 : Math.exp(-Math.max(cdH - q.reactH, 0) / q.shutdownMean), nLong = q.shutdowns * pLong;
  const o = { key, name: st.name, insT, U, tArr: th.tArr, sub: th.sub, cooldownH: cdH, pLong, nLong, inhWt: 0, megRate: 0, chemRate: 0.5, heatKW: 0, chem: 0, energy: 0, shutdown: 0, extraDownH: 0, events: 0, kWh: 0, heatGJ: 0, chemT: 0, contRate: 0, feasible: true, note: '' };
  const shutNeed = S.base.tHydShut - S.Ta + q.hydMargin, perKgW = (w) => w / Math.max(100 - w, 1), av = q.inhibAvail / 100;
  if (!S.shut) { const w = wtFor(shutNeed, 'MeOH'); S.shut = { w, megW: wtFor(shutNeed, 'MEG'), vol: (S.base.waterInv * 1000 * perKgW(w)) / INHIBITORS.MeOH.rho }; } // m³ methanol to protect the line contents
  const volShut = S.shut.vol;
  if (key === 'bare') {
    o.inhWt = need > shutNeed ? wtFor(need, 'MEG') : S.shut.megW;
    if (o.inhWt >= 80) { o.feasible = false; o.note = 'MEG demand exceeds the lean-glycol strength'; }
    const x = o.inhWt / Math.max(90 - o.inhWt, 5), lean = (m3) => (m3 * 1000 * x) / 1100; // m³/d of 90 wt % lean MEG
    o.megRate = lean(waterPeak); o.chemRate = o.megRate;
    const makeup = lean(water) * 0.9 * (q.megLoss / 100) * 365 * up; // m³/y of MEG lost
    o.chem = makeup * q.megPrice; o.chemT = makeup * 1.113 * EF.MEG;
    o.heatGJ = water * 3.0 * 365 * up; o.energy = o.heatGJ * X.heatPrice; // regeneration boils the produced water off: about 3 MJ per kg
    o.events = (q.shutdowns + (need > 0 ? 12 : 0)) * q.plugProb * (1 - av); o.note = o.note || `${o.inhWt.toFixed(0)} wt % MEG in the water phase, ${o.megRate.toFixed(0)} m³/d lean MEG`;
  } else if (key === 'ldhi') {
    const rate = (q.ldhiDose / 100) * water; o.chemRate = (q.ldhiDose / 100) * waterPeak + 0.5;
    o.chem = rate * q.ldhiPrice * 365 * up; o.chemT = rate * 365 * up * 0.95 * EF.LDHI; o.contRate = o.chemRate;
    const eff = clamp(0.9 - 1.5 * Math.max(X.wcMean - 0.5, 0), 0, 0.9) * av;
    o.events = (need > 0 ? 12 : nLong) * q.plugProb * (1 - eff); o.note = `${rate.toFixed(1)} m³/d anti-agglomerant, effectiveness ${(eff * 100).toFixed(0)} % at ${(X.wcMean * 100).toFixed(0)} % water cut`;
  } else if (key === 'risk') {
    o.events = (need > 0 ? 12 : nLong) * q.plugProb; o.note = need > 0 ? 'flowing inside the hydrate region with no inhibition' : `${nLong.toFixed(1)} unprotected long shutdowns a year`;
  } else {
    if (need > 0) { // steady flow inside the hydrate region: continuous once-through methanol
      o.inhWt = wtFor(need, 'MeOH'); const rate = (water * 1000 * perKgW(o.inhWt)) / INHIBITORS.MeOH.rho;
      o.chemRate = (waterPeak * 1000 * perKgW(o.inhWt)) / INHIBITORS.MeOH.rho; o.contRate = o.chemRate; o.chem = rate * q.meohPrice * 365 * up; o.chemT = rate * 365 * up * 0.792 * EF.MeOH;
      o.note = `continuous methanol ${rate.toFixed(1)} m³/d (${o.inhWt.toFixed(0)} wt %) because steady flow is ${need.toFixed(1)} °C inside the hydrate margin`;
    }
    if (key === 'deh') {
      o.heatKW = (U * Math.PI * S.id * X.capexArgs.flowLen * Math.max(S.base.tHydShut + q.hydMargin + 5 - S.Ta, 0)) / 0.75 / 1000;
      o.kWh = nLong * q.shutdownMean * o.heatKW; o.energy = o.kWh * q.elecPrice; o.events = (nLong * 0.03 + (need > 0 ? 12 * (1 - av) : 0)) * q.plugProb; o.note = o.note || `${o.heatKW.toFixed(0)} kW of heating holds the line above the hydrate temperature`;
    } else {
      const vol = nLong * volShut * treat; o.chemRate = Math.max(o.chemRate, Math.min(volShut, 100)); // skid sized to dose the line contents within a day
      o.shutdown = vol * q.meohPrice; o.chemT += vol * 0.792 * EF.MeOH; o.extraDownH = nLong * q.restartH * treat; o.events = (nLong * (1 - 0.97 * av * treat) + (need > 0 ? 12 * (1 - av) : 0)) * q.plugProb;
      o.note = o.note || `cooldown ${cdH >= 1e4 ? 'never reaches hydrate conditions' : cdH.toFixed(0) + ' h'}, ${nLong.toFixed(1)} long shutdowns a year, ${volShut.toFixed(0)} m³ methanol each`;
    }
  }
  if (over.inhibRate > 0) { // continuous inhibitor rate reported by operations replaces the modelled demand
    const meg = q.inhibitor === 'MEG', price = meg ? q.megPrice * (q.megLoss / 100) : q.inhibitor === 'LDHI' ? q.ldhiPrice : q.meohPrice;
    o.chem = over.inhibRate * price * 365 * up; o.chemT = over.inhibRate * 365 * up * (meg ? (1.113 * EF.MEG * q.megLoss) / 100 : q.inhibitor === 'LDHI' ? 0.95 * EF.LDHI : 0.792 * EF.MeOH); o.chemRate = Math.max(o.chemRate, over.inhibRate);
  }
  if (over.heatKW > 0) { o.kWh += over.heatKW * 8760 * up; o.energy += over.heatKW * 8760 * up * X.enPrice; }
  o.pigRuns = q.pigRuns * (1 + Math.max(q.wat - th.tArr, 0) / 10); o.pig = o.pigRuns * q.pigCost * 1000;
  o.block = o.events * X.blockCost; o.deferral = o.extraDownH * X.marginH * X.lossFrac;
  o.carbonT = o.chemT + (o.kWh * X.enCarbon) / 1000 + (o.heatGJ * EF.fuelGJ) / 0.85 + o.events * q.remedDays * EF.vesselDay;
  // engineering safety screening: an option that breaks one of these is reported but never ranked best
  const relies = key === 'wet' || key === 'pip' || key === 'risk', v1 = Math.max(0, o.events / Math.max(q.maxBlockFreq, 1e-9) - 1), v2 = relies && cdH < q.minCooldown ? 1 - cdH / q.minCooldown : 0, v3 = (o.contRate || 0) > q.maxInject ? o.contRate / q.maxInject - 1 : 0;
  o.viol = [];
  if (v1 > 0) o.viol.push(`blockage frequency ${o.events.toFixed(3)} per year above the tolerable ${q.maxBlockFreq}`);
  if (v2 > 0) o.viol.push(`cooldown ${cdH.toFixed(1)} h shorter than the required ${q.minCooldown} h`);
  if (v3 > 0) o.viol.push(`continuous inhibitor demand ${o.contRate.toFixed(0)} m³/d above the ${q.maxInject} m³/d injection capacity`);
  if (key === 'risk' && need > 0) o.viol.push('steady flow inside the hydrate region with no inhibition');
  if (!o.feasible) o.viol.push(o.note);
  o.vio = v1 * v1 + v2 * v2 + v3 * v3 + (o.viol.length > (v1 > 0) + (v2 > 0) + (v3 > 0) ? 1 : 0); o.safe = o.viol.length === 0;
  return o;
}

/** Emission inventory (tCO₂e per year) for a strategy at rate multiplier r. */
function emissions(q, S, X, sm, r = 1) {
  const up = X.up, kWh = (q.pumpKW + q.compKW) * 8760 * up * r, gasY = mean(X.prof.gasSm3) * r;
  const items = [
    { source: 'Power for pumping and compression', t: (kWh * X.enCarbon) / 1000 },
    { source: 'Heating and electrical tracing', t: (sm.kWh * X.enCarbon) / 1000 },
    { source: 'Glycol regeneration heat (fuel gas)', t: (sm.heatGJ * EF.fuelGJ) / 0.85 },
    { source: 'Routine flaring', t: (q.flareFrac / 100) * gasY * EF.flareSm3 },
    { source: 'Blowdown flaring', t: q.blowdowns * S.base.gasInvStd * EF.flareSm3 },
    { source: 'Chemicals (embodied)', t: sm.chemT },
    { source: 'Intervention and inspection vessels', t: (sm.events * q.remedDays + 4 / q.inspInterval) * EF.vesselDay },
  ];
  return { items, total: sum(items.map((i) => i.t)) };
}

/** Assemble the cash-flow model of the project for one strategy result and CAPEX estimate. */
function buildProject(q, S, X, sm, capex, over = {}) {
  const prof = over.profile || X.prof, life = q.life, cap = capex.total, r = over.rate ?? 1;
  const fixedItems = [
    ['Operations, labour and logistics', q.opsFixed * MM], ['Maintenance', (q.maintPct / 100) * cap], ['Insurance', (q.insurPct / 100) * cap], ['Inspection', (q.inspCost * MM) / q.inspInterval], ['Corrosion management', q.corrMgmt * MM],
    ['Hydrate inhibitor and flow-assurance chemicals', sm.chem], ['Energy (pumping, compression, heating)', (q.pumpKW + q.compKW) * 8760 * X.up * r * X.enPrice + sm.energy], ['Pigging', sm.pig],
  ];
  const em = over.carbonT ?? emissions(q, S, X, sm, r).total, dFrac = X.deprFrac || depreciation('uop', 1, q.deprLife, { years: life, units: prof.boe });
  return {
    phase: q.phase, capex: cap, capexSunk: q.capexSunk * MM, residual: q.residual * MM, life, lifeCut: over.lifeCut ?? life, oil: prof.oil, gas: prof.gas, water: prof.water,
    uptime: clamp(X.up - sm.extraDownH / 8760, 0.05, 1), deferFrac: q.deferFrac / 100, oilPrice: q.oilPrice, gasPrice: q.gasPrice, infl: X.infl, costEsc: q.costEsc / 100, priceEsc: q.priceEsc / 100, carbonEsc: q.carbonEsc / 100,
    discount: X.r, mid: q.mid, reinvest: q.reinvest / 100, opexFixed: sum(fixedItems.map((i) => i[1])), opexVarBoe: q.tariff + q.chemOther, waterCost: q.waterCost, opexDown: sm.shutdown, opexBlock: sm.block,
    carbonT: em, carbonPrice: q.carbonPrice, includeRisk: over.includeRisk ?? false, consequence: q.consequence * MM, haz: X.haz,
    royalty: q.royalty / 100, taxRate: q.taxRate / 100, regime: q.regime, costOilCap: q.costOilCap / 100, profitSplit: q.profitSplit / 100, deprFrac: dFrac, wcDays: q.wcDays, abandon: q.abandon * MM, abandonProvision: q.abandonProvision,
    gearing: q.gearing / 100, loanRate: q.loanRate / 100, loanTenor: q.loanTenor, fixedItems,
  };
}
/** Strategy → CAPEX → project in one call (used by every option study). */
function optionProject(key, q, S, X, over = {}) {
  const sm = strategyModel(key, q, S, X, over), capex = capexEstimate({ ...X.capexArgs, brief: !!over.brief, strategy: key, insT: sm.insT, megRate: sm.megRate, chemRate: sm.chemRate, heatKW: sm.heatKW, ...(over.capex || {}) });
  return { sm, capex, p: buildProject(q, S, X, sm, capex, over) };
}

// ================================================================================================================
// 9. The suite engine
// ================================================================================================================
const VARS = ['price', 'prod', 'capex', 'opex', 'downtime', 'failFreq', 'repair'];
const VAR_LABEL = { price: 'Commodity price', prod: 'Production volume', capex: 'CAPEX', opex: 'OPEX', downtime: 'Downtime', failFreq: 'Failure frequency', repair: 'Repair / consequence cost' };
const mOf = (x) => ({ price: x[0], prod: x[1], capex: x[2], opex: x[3], downtime: x[4], failFreq: x[5], repair: x[6] });

async function run(v, ctx = {}) {
  const q = readInputs(v), prog = (f, m) => { try { ctx.progress?.(f, m); } catch { /* progress is optional */ } }, tick = async () => { if (typeof ctx.tick === 'function') await ctx.tick(); };
  const warnings = [], recs = [], tables = [], plots = [], balances = [], warn = (level, msg) => warnings.push({ level, msg });
  prog(0.01, 'Solving the case line on the flow kernel');
  const S = await physics(q, ctx, prog, tick), X = context(q, S), sel = q.strategy, keys = Object.keys(STRATEGIES), nCon = X.nCon, rf = X.r, dOpt = { mid: q.mid };
  S.warn.forEach((m) => warn('warn', m));

  // ---------------------------------------------------------------- base project, cash flow, metrics
  prog(0.42, 'Cash flow and investment metrics');
  let base = optionProject(sel, q, S, X, { inhibRate: q.inhibRate, heatKW: q.heatKW }), p0 = base.p, cf = cashflow(p0), met = investmentMetrics(cf, p0), lifeEff = q.life;
  const limit0 = met.economicLimit;
  if (q.stopAtLimit && limit0 !== null && limit0 >= 1 && limit0 < q.life) { lifeEff = limit0; p0 = { ...p0, lifeCut: lifeEff }; base.p = p0; cf = cashflow(p0); met = investmentMetrics(cf, p0); }
  X.lifeCut = lifeEff;
  const cut = (p) => ({ ...p, lifeCut: lifeEff }), K = cf.K, pR = { ...p0, includeRisk: true }, cfR = cashflow(pR), riskedNpv = cfR.npv, npv0 = cf.npv;
  const capex = base.capex, sm0 = base.sm, em0 = emissions(q, S, X, sm0), y1 = nCon, boeY1 = cf.boe[y1], opexY1 = cf.opex[y1] / X.cg ** y1, meanBoe = sum(cf.boe) / lifeEff;
  const eal = p0.haz[0] * p0.consequence, lifeFail = npv0 - riskedNpv, downH = (1 - cf.uptime) * 8760, deferredCost = downH * X.marginH * X.lossFrac;
  // break-evens (root finding on the fast cash-flow path)
  const up = X.up, inhPrice = q.inhibitor === 'MEG' ? q.megPrice : q.inhibitor === 'LDHI' ? q.ldhiPrice : q.meohPrice;
  const be = {
    price: breakeven((x) => cashflow({ ...p0, oilPrice: x }, {}, false), 0, 3000),
    prod: breakeven((x) => cashflow(p0, { prod: x }, false), 0.005, 50),
    capex: breakeven((x) => cashflow(p0, { capex: x }, false), 0.01, 100),
    inhib: breakeven((x) => cashflow({ ...p0, opexFixed: p0.opexFixed + x * inhPrice * 365 * up }, {}, false), 0, 1e6),
    block: breakeven((x) => cashflow({ ...p0, opexBlock: x * X.blockCost }, {}, false), 0, 1e4),
    uptime: breakeven((x) => cashflow({ ...p0, uptime: x }, {}, false), 0.05, 1),
  };
  const bePriceCheck = be.price === null ? 0 : cashflow({ ...p0, oilPrice: be.price }, {}, false);

  // ---------------------------------------------------------------- hybrid 1: flow-assurance strategies
  prog(0.48, 'Flow-assurance strategy comparison');
  const pvfK = sum(Array.from({ length: lifeEff }, (_, j) => ((1 + X.infl) * (1 + q.carbonEsc / 100)) ** (nCon + j) * discountFactor(rf, nCon + j, dOpt))), pvfL = sum(Array.from({ length: lifeEff }, (_, j) => X.cg ** (nCon + j) * discountFactor(rf, nCon + j, dOpt)));
  const options = keys.map((k) => {
    const o = optionProject(k, q, S, X, { brief: true }), s = o.sm, emis = emissions(q, S, X, s).total, parts = { capex: o.capex.total - X.capexNone, chem: pvfL * (s.chem + s.energy), pig: pvfL * s.pig, shut: pvfL * (s.shutdown + s.deferral), block: pvfL * s.block, carbon: pvfK * s.carbonT * q.carbonPrice };
    return { key: k, name: s.name, sm: s, capex: o.capex, p: cut(o.p), emis, parts, lcc: sum(Object.values(parts)), npv: cashflow(cut(o.p), {}, false), npvRisked: cashflow({ ...cut(o.p), includeRisk: true }, {}, false) };
  });
  // safety screening: options that break an engineering constraint are listed last and never recommended
  const safeOpts = options.filter((o) => o.sm.safe), pool = safeOpts.length ? safeOpts : options, noneSafe = !safeOpts.length, byCost = (a, b) => a.lcc - b.lcc;
  const byLcc = [...safeOpts.slice().sort(byCost), ...options.filter((o) => !o.sm.safe).sort(byCost)], bestOpt = byLcc[0], selOpt = options.find((o) => o.key === sel);

  // ---------------------------------------------------------------- hybrid 2: insulation thickness optimum
  const tGrid = linspace(0, q.tMaxMm, Math.round(q.tMaxMm) + 1), insCurve = tGrid.map((t) => { const o = optionProject('wet', q, S, X, { insT: t / 1000, brief: true }), s = o.sm; return { t, safe: s.safe, capex: o.capex.total - X.capexNone, opex: pvfL * (s.chem + s.energy + s.pig + s.shutdown), defer: pvfL * s.deferral, block: pvfL * s.block, carbon: pvfK * s.carbonT * q.carbonPrice, tArr: s.tArr, sub: s.sub, cd: s.cooldownH, o }; });
  insCurve.forEach((c) => (c.lcc = c.capex + c.opex + c.defer + c.block + c.carbon));
  const insSafe = insCurve.filter((c) => c.safe), insPool = insSafe.length ? insSafe : insCurve, insOpt = insPool[argBest(insPool.map((c) => -c.lcc))], iOpt = insCurve.indexOf(insOpt), insFree = insCurve[argBest(insCurve.map((c) => -c.lcc))], insSel = insCurve[argBest(insCurve.map((c) => -Math.abs(c.t - q.insMm)))], insBare = insCurve[0];

  // ---------------------------------------------------------------- hybrid 3: diameter optimisation
  prog(0.52, 'Diameter and rate optimisation');
  const diam = S.diam.map((d) => {
    if (d.failed) return { ...d, feasible: false, npv: null, why: d.failed };
    const need = (x) => S.pInOf(x, d.dpFric, d.dpGrav) - X.pAvail(x), rate = need(1) <= 0 ? 1 : need(0.05) >= 0 ? 0.05 : brent(need, 0.05, 1, 1e-9), wtD = (S.wt * d.d) / S.id;
    const o = optionProject(sel, q, S, X, { rate, thermal: { tArr: d.tArr, sub: d.sub }, profile: profileOf(q, rate), capex: { id: d.d, wt: wtD, slugVol: d.slugVol }, inhibRate: q.inhibRate, heatKW: q.heatKW, brief: true }), er = (d.eros * rate * S.base.eros) / Math.max(S.base.erosKernel, 1e-9), pInD = S.pInOf(rate, d.dpFric, d.dpGrav), feasible = er <= 1 && pInD <= q.mawp;
    return { ...d, eros: er / Math.max(rate, 1e-9), pInD, rate, feasible, capexTotal: o.capex.total, slugCost: o.capex.items.find((i) => i.item.startsWith('Slug catcher'))?.cost ?? 0, npv: cashflow(cut(o.p), {}, false), why: feasible ? (rate < 0.999 ? 'back-pressure limits the rate' : 'full rate') : er > 1 ? 'excluded: erosional velocity exceeded' : 'excluded: inlet pressure above the design pressure' };
  });
  const dOk = diam.filter((d) => d.feasible && d.npv !== null), dBest = dOk.length ? dOk[argBest(dOk.map((d) => d.npv))] : null, dBase = diam.find((d) => d.mult === 1);

  // ---------------------------------------------------------------- hybrid 6: NPV-optimal operating point (rate, preservation fraction)
  const pen = (x) => sum(X.constraints(x).map((c) => Math.max(0, c.g) ** 2)), opProj = (x, t, key = sel) => optionProject(key, q, S, X, { rate: x, treat: t, profile: profileOf(q, x), inhibRate: key === sel ? q.inhibRate : 0, heatKW: key === sel ? q.heatKW : 0, brief: true });
  const evalOp = (x, t, key) => { const o = opProj(x, t, key); return { npv: cashflow(cut(o.p), {}, false), vio: o.sm.vio }; }, npvOp = (x, t, key) => evalOp(x, t, key).npv, vio0 = selOpt.sm.vio;
  const penAll = (x, t, key) => pen(x) + Math.max(0, evalOp(x, t, key).vio - (key && key !== sel ? 0 : vio0)), fOp = (z) => { const e = evalOp(z[0], z[1]); return -e.npv / MM + 1e5 * (pen(z[0]) + Math.max(0, e.vio - vio0)); }, lo = [Math.max(0.2, X.rMin - 0.2), 0], hi = [X.rMax + 0.2, 1], z0 = [1, 1];
  const solvers = [
    ['Nelder–Mead simplex', () => { const r = nelderMead(fOp, z0, { lo, hi, tol: 1e-9, maxIter: 70, scale: 0.15 }); return { x: r.x, f: r.f, evals: r.evals }; }],
    ['Projected gradient (finite differences)', () => gradientDescent(fOp, z0, { lo, hi, maxIter: 20 })],
    ['Differential evolution', () => { const r = diffEvolution(fOp, lo, hi, { pop: 10, gens: 14, seed: q.seed }); return { x: r.x, f: r.f, evals: r.evals ?? 10 * 15 }; }],
    ['Genetic algorithm', () => geneticAlgorithm(fOp, lo, hi, { pop: 14, gens: 11, seed: q.seed })],
    ['Particle swarm', () => particleSwarm(fOp, lo, hi, { n: 10, iters: 14, seed: q.seed })],
  ].map(([name, fn]) => { const r = fn(); return { name, x: r.x, f: r.f, evals: r.evals ?? null, npv: npvOp(r.x[0], r.x[1]), pen: penAll(r.x[0], r.x[1]) }; });
  const feasS = solvers.filter((s) => s.pen < 1e-6), opBest = (feasS.length ? feasS : solvers).reduce((a, s) => (s.f < a.f ? s : a)), opCons = X.constraints(opBest.x[0]), active = opCons.filter((c) => Math.abs(c.g) < 0.01).map((c) => c.name);
  const rGrid = linspace(lo[0], hi[0], 31), rateCurve = rGrid.map((x) => ({ r: x, npv: npvOp(x, 1), ok: pen(x) < 1e-9 })), envOk = rateCurve.filter((c) => c.ok && c.npv > 0), envelope = envOk.length ? { lo: envOk[0].r, hi: envOk[envOk.length - 1].r } : null;
  // mixed-integer nonlinear: enumerate the discrete strategy, optimise the continuous rate inside each
  const minlp = keys.map((k) => { const f = (z) => { const e = evalOp(z[0], 1, k); return -e.npv / MM + 1e5 * (pen(z[0]) + e.vio); }, r = nelderMead(f, [1], { lo: [lo[0]], hi: [hi[0]], tol: 1e-7, maxIter: 25, scale: 0.1 }), e = evalOp(r.x[0], 1, k); return { key: k, name: STRATEGIES[k].name, rate: r.x[0], npv: e.npv, feasible: pen(r.x[0]) + e.vio < 1e-6 }; }).sort((a, b) => b.feasible - a.feasible || b.npv - a.npv);
  await tick();

  // ---------------------------------------------------------------- reliability economics: inspection interval, spares, replacement
  prog(0.6, 'Reliability economics');
  const relArg = { life: lifeEff, inspCost: q.inspCost * MM, consequence: p0.consequence, h0: p0.haz[0] - (hazardRate(0.5, { h0: 0, beta: q.weibullBeta, remLife: q.remLife, pEnd: q.pEnd / 100 })), beta: q.weibullBeta, remLife: q.remLife, pEnd: q.pEnd / 100, pod: q.pod / 100, rate: X.rr };
  relArg.h0 = Math.max(relArg.h0, 0);
  const iGrid = linspace(0.5, Math.max(1, Math.min(lifeEff, 20)), Math.round(2 * Math.max(1, Math.min(lifeEff, 20)))), rbi = iGrid.map((t) => ({ t, ...inspectionCost({ ...relArg, interval: t }) })), rbiOk = rbi.filter((x) => x.hMax <= q.maxPof), rbiPool = rbiOk.length ? rbiOk : [rbi[argBest(rbi.map((x) => -x.hMax))]], rbiBest = rbiPool[argBest(rbiPool.map((x) => -x.total))], rbiFree = rbi[argBest(rbi.map((x) => -x.total))], rbiNow = inspectionCost({ ...relArg, interval: q.inspInterval }), rbiNone = inspectionCost({ ...relArg, interval: 1e9 });
  const dayValue = 24 * X.marginH * X.lossFrac, spare = { hold: q.spareCost * MM * (X.rr + 0.03), without: q.itemFailRate * q.leadNo * dayValue * (q.spareShare / 100), with: q.itemFailRate * q.leadWith * dayValue * (q.spareShare / 100) };
  spare.saving = spare.without - spare.with - spare.hold; spare.beRate = spare.hold / Math.max((q.leadNo - q.leadWith) * dayValue * (q.spareShare / 100), 1e-9);
  const opAge = (a) => q.maintAsset * MM * (1 + 0.04 * a) + Math.min(1, hazardRate(a + 0.5, { h0: q.pof, beta: q.weibullBeta, remLife: q.remLife + q.assetAge, pEnd: q.pEnd / 100 })) * p0.consequence;
  const dp = replacementDP({ horizon: lifeEff, maxAge: Math.round(q.assetAge + lifeEff + 1), opCost: opAge, replaceCost: q.replCost * MM, rate: X.rr, age0: q.assetAge }), dpNever = sum(Array.from({ length: lifeEff }, (_, t) => opAge(q.assetAge + t) * (1 + X.rr) ** -t));
  // hybrid 5: pigging-interval optimisation, probability × consequence
  const slope = X.qLiq0 / q.pi + S.nf * S.base.dpFric, pigAt = (tau) => { const dEnd = q.waxRate * tau, fac = (1 - Math.min((dEnd / 1000) / S.id, 0.45)) ** -5 - 1, lossFr = Math.max(0, S.base.dpFric * fac - q.chokeDp) / slope, runs = 365 / tau, pStuck = 1 - Math.exp(-((dEnd / q.waxCrit) ** 3)); return { tau, runs, pig: runs * q.pigCost * 1000, loss: lossFr * 8760 * up * X.marginH * X.lossFrac, stuck: runs * pStuck * X.blockCost, pStuck, dEnd }; };
  const pigCurve = linspace(Math.log(2), Math.log(365), 80).map(Math.exp).map(pigAt).map((c) => ({ ...c, total: c.pig + c.loss + c.stuck })), pigOk = pigCurve.filter((c) => c.dEnd <= q.waxCrit), pigPool = pigOk.length ? pigOk : [pigCurve[0]], pigBest = pigPool[argBest(pigPool.map((c) => -c.total))], pigNowTau = 365 / Math.max(sm0.pigRuns, 0.1), pigNow = (() => { const c = pigAt(clamp(pigNowTau, 1, 3650)); return { ...c, total: c.pig + c.loss + c.stuck }; })();
  // hybrid 4: corrosion allowance vs inhibition vs CRA
  const eff = q.inhEff / 100, crAs = q.corrRate, crUn = q.corrInhibited ? crAs / Math.max(1 - eff, 0.02) : crAs, crIn = q.corrInhibited ? crAs : crAs * (1 - eff), capArgs0 = { ...X.capexArgs, strategy: sel, insT: sm0.insT, megRate: sm0.megRate, chemRate: sm0.chemRate, heatKW: sm0.heatKW };
  const matOpt = (name, rate, ca, capOver, opexY, h0f) => {
    const caUse = Math.min(ca, 10), rl = rate > 1e-6 ? caUse / rate : 1e6, cap = capexEstimate({ ...capArgs0, brief: true, ...capOver, caExtra: Math.max(caUse - q.corrAllow, 0) / 1000 }), dCap = cap.total - capex.total;
    let fail = 0;
    for (let j = 0; j < lifeEff; j++) fail += Math.min(1, hazardRate(j + 0.5, { h0: q.pof * h0f, beta: q.weibullBeta, remLife: Math.max(rl, 0.5), pEnd: q.pEnd / 100 })) * p0.consequence * X.cg ** (nCon + j) * discountFactor(rf, nCon + j, dOpt);
    const repl = rl < lifeEff ? cap.pipeInstalled * X.cg ** (nCon + rl) * discountFactor(rf, nCon + rl) : 0;
    return { name, rate, ca: caUse, life: Math.min(rl, 999), dCap, opex: pvfL * opexY, fail, repl, lcc: dCap + pvfL * opexY + fail + repl, practical: ca <= 10 };
  };
  const mats = [matOpt('Carbon steel + corrosion allowance, no inhibitor', crUn, Math.max(crUn * lifeEff, q.corrAllow), {}, 0, 1), matOpt('Carbon steel + corrosion inhibitor', crIn, Math.max(crIn * lifeEff, q.corrAllow), {}, q.corrInhCost * MM, 1), matOpt('CRA-clad pipe', 0, 0, { material: 'cra' }, 0, 0.2)], matPool = mats.filter((m) => m.practical), matBest = matPool[argBest(matPool.map((m) => -m.lcc))];

  // ---------------------------------------------------------------- uncertainty: sampling, scenarios, sensitivities
  prog(0.66, 'Monte Carlo simulation');
  const ds = VARS.map((id) => q.dists[id]), Rm = rng(q.seed + 17), alpha = q.alpha / 100;
  const mc = monteCarlo((x) => { const m = mOf(x); if (q.priceModel !== 'static') m.path = pricePath(q.priceModel, K, Rm, { sigma: q.priceVol / 100, kappa: q.priceKappa }); if (q.failEvents) m.events = Array.from({ length: lifeEff }, () => Rm.uniform()); return cashflow(pR, m, false); }, ds, { n: q.nMC, method: q.sampling, corr: q.corr, seed: q.seed, alpha });
  if (mc.shrink > 0) warn('warn', `The correlation matrix is not positive definite; its off-diagonal terms were shrunk by ${(mc.shrink * 100).toFixed(0)} % before sampling.`);
  const conv = []; { let s = 0; const step = Math.max(1, Math.floor(mc.y.length / 60)); mc.y.forEach((y, i) => { s += y; if ((i + 1) % step === 0 || i === mc.y.length - 1) conv.push({ n: i + 1, mean: s / (i + 1) }); }); }
  const se = mc.sd / Math.sqrt(mc.n), ce = certaintyEquivalent(mc.y, q.riskTol * MM);
  await tick();
  const scen = q.scenarios.map((s) => ({ ...s, npv: cashflow(pR, { price: s.price, prod: s.prod, capex: s.capex, opex: s.opex }, false) })), wSum = sum(scen.map((s) => s.weight)), scenEv = sum(scen.map((s) => (s.weight / wSum) * s.npv));
  const fVec = (x) => cashflow(pR, mOf(x), false), ones = VARS.map(() => 1), torn = tornado(fVec, ones, VARS.map((id, j) => ({ name: VAR_LABEL[id], lo: ds[j].inv(0.1), hi: ds[j].inv(0.9) })));
  const spiderF = linspace(0.7, 1.3, 7), spider = [['Commodity price', (f) => cashflow(pR, { price: f }, false)], ['Production volume', (f) => cashflow(pR, { prod: f }, false)], ['CAPEX', (f) => cashflow(pR, { capex: f }, false)], ['OPEX', (f) => cashflow(pR, { opex: f }, false)], ['Downtime', (f) => cashflow(pR, { downtime: f }, false)], ['Discount rate', (f) => cashflow({ ...pR, discount: X.r * f }, {}, false)]].map(([name, f]) => ({ name, y: spiderF.map((x) => f(x) / MM) }));
  const live = VARS.map((_, j) => j).filter((j) => ds[j].variance > 0), sob = live.length ? sobolIndices((u) => { const x = ones.slice(); live.forEach((j, i) => (x[j] = ds[j].inv(clamp(u[i], 1e-9, 1 - 1e-9)))); return fVec(x); }, live.length, q.nSobol, q.seed + 3) : { first: [], total: [], evals: 0 };
  const srcFit = live.length > 0 && mc.sd > 0 ? standardisedRegression(mc.X.map((r) => live.map((j) => r[j])), mc.y) : { src: [], r2: 0 };

  // ---------------------------------------------------------------- hybrid 7: insulation decision under uncertainty (common random numbers)
  prog(0.74, 'Probabilistic lifecycle optimisation');
  const crn = sampleCorrelated(ds, { n: q.nMCopt, method: 'lhs', corr: q.corr, seed: q.seed + 5 }).X.map(mOf), probIns = insCurve.filter((c, i) => (i % Math.max(1, Math.round(insCurve.length / 15)) === 0 || i === iOpt) && (c.safe || !insSafe.length)).map((c) => { const pp = { ...cut(c.o.p), includeRisk: true }, st = riskStats(crn.map((m) => cashflow(pp, m, false)), alpha); return { t: c.t, mean: st.mean, cvar: st.cvar, p10: st.p10 }; });
  const pMean = probIns[argBest(probIns.map((c) => c.mean))], pCvar = probIns[argBest(probIns.map((c) => c.cvar))];
  await tick();

  // ---------------------------------------------------------------- decision analysis
  prog(0.8, 'Decision analysis');
  const alts = [...pool.map((o) => ({ name: o.name, pay: scen.map((s) => cashflow({ ...o.p, includeRisk: true }, { price: s.price, prod: s.prod, capex: s.capex, opex: s.opex }, false)) })), { name: 'Do not sanction', pay: scen.map(() => 0) }];
  const prior = scen.map((s) => s.weight / wSum), nS = scen.length, rel = q.testRel / 100, like = scen.map((_, k) => scen.map((__, s) => (nS === 1 ? 1 : k === s ? rel : (1 - rel) / (nS - 1))));
  const voi = valueOfInformation({ prior, payoff: alts.map((a) => a.pay), likelihood: like });
  const chance = (name, probs, pay) => ({ name, type: 'chance', branches: scen.map((s, i) => ({ name: s.name, p: probs[i], value: pay[i] })) });
  const tree = decisionTree({ name: 'Sanction decision', type: 'decision', branches: [
    { name: 'Decide now', node: { name: 'Concept selection', type: 'decision', branches: alts.map((a) => ({ name: a.name, node: a.name === 'Do not sanction' ? { name: a.name, value: 0 } : chance(a.name, prior, a.pay) })) } },
    { name: 'Appraise first', cost: q.testCost * MM, node: { name: 'Appraisal result', type: 'chance', branches: scen.map((s, k) => ({ name: `Result indicates ${s.name}`, p: voi.signalProb[k], node: { name: 'Concept selection', type: 'decision', branches: alts.map((a) => ({ name: a.name, node: a.name === 'Do not sanction' ? { name: a.name, value: 0 } : chance(a.name, voi.posterior[k], a.pay) })) } })) } },
  ] });
  const ahpR = ahp(q.ahpM), crit = ['Value (risked NPV)', 'Blockage risk', 'Carbon', 'Operability'], benefit = [true, false, false, true], mcdaM = pool.map((o) => [o.npvRisked / MM, o.sm.block / MM, o.emis, STRATEGIES[o.key].operability]);
  const ws = weightedSum(mcdaM, ahpR.weights, benefit), tp = topsis(mcdaM, ahpR.weights, benefit);
  options.forEach((o) => { o.score = 0; o.closeness = 0; });
  pool.forEach((o, i) => { o.score = ws.scores[i]; o.closeness = tp.closeness[i]; });
  const ranking = byLcc.map((o) => ({ option: o.name, npv: o.npv, capex: o.capex.total, risk: o.sm.block, score: +o.score.toFixed(4), lifecycleCost: o.lcc, feasible: o.sm.safe, violations: o.sm.viol.join('; ') })), mcdaBest = pool[argBest(ws.scores)], topsisBest = pool[tp.order[0]];
  const regret = minimaxRegret(pool.map((o, i) => alts[i].pay)), regretBest = pool[regret.best];
  // real options on the project value
  const invPV = sum(cf.capex.map((c, k) => c * cf.df[k])), Vop = Math.max(npv0 + invPV, 0), sig = q.optVol / 100, rfree = q.riskFree / 100, yld = q.optYield / 100, optArg = { S: Vop, K: Math.max(invPV, 1), r: rfree, sigma: sig, T: q.optYears, q: yld };
  const defer = binomialOption({ ...optArg, steps: q.nLattice, type: 'call', american: true }).value, deferEu = binomialOption({ ...optArg, steps: q.nLattice, type: 'call' }).value, bs = blackScholes({ ...optArg, type: 'call' });
  const expand = binomialOption({ ...optArg, steps: q.nLattice, payoff: (s) => Math.max((q.expandFrac / 100) * s - q.expandCost * MM, 0) }).value, abandonOpt = binomialOption({ S: Vop, K: (q.salvage / 100) * capex.total, r: rfree, sigma: sig, T: Math.min(lifeEff, 10), q: yld, steps: q.nLattice, type: 'put', american: true }).value;
  const volCurve = linspace(5, 60, 12).map((s) => ({ s, lat: binomialOption({ ...optArg, sigma: s / 100, steps: q.nLattice, type: 'call' }).value, am: binomialOption({ ...optArg, sigma: s / 100, steps: q.nLattice, type: 'call', american: true }).value, bs: blackScholes({ ...optArg, sigma: s / 100, type: 'call' }) }));
  await tick();

  // ---------------------------------------------------------------- optimisation: LP, MILP, two-stage, Pareto
  prog(0.86, 'Portfolio and multi-objective optimisation');
  const pf = q.portfolio, lpArg = { c: pf.map((r) => r.npv), A: [pf.map((r) => r.capex), pf.map((r) => r.days), ...pf.map((_, i) => pf.map((__, j) => (i === j ? 1 : 0)))], b: [q.budget, q.vesselDays, ...pf.map(() => 1)] };
  const lp = pf.length ? simplex(lpArg) : { status: 'empty', x: [], value: 0 }, milp = pf.length ? branchAndBound(lpArg) : { status: 'empty', x: [], value: 0, nodes: 0 };
  const annF = annuityPV(1, X.rr, lifeEff), ts2 = twoStage({ c: q.waterCapCost * 1000, q: q.waterPenalty * 365 * annF, scenarios: scen.map((s) => ({ p: s.weight, demand: X.waterPeak * s.prod })) });
  const tMaxM = q.tMaxMm / 1000, pfEval = (z) => { const o = optionProject('wet', q, S, X, { insT: z[0], rate: z[1], treat: z[2], profile: profileOf(q, z[1]), brief: true }), pe = pen(z[1]) + o.sm.vio; return [-cashflow(cut(o.p), {}, false) / MM + 1e5 * pe, (o.sm.block + o.sm.deferral + eal) / MM + 1e3 * pe, emissions(q, S, X, o.sm, z[1]).total + 1e7 * pe]; };
  const front = nsga2(pfEval, [0, lo[0], 0], [tMaxM, hi[0], 1], { pop: q.nPop, gens: q.nGens, seed: q.seed + 9 }).filter((s) => pen(s.x[1]) < 1e-9 && strategyModel('wet', q, S, X, { insT: s.x[0], rate: s.x[1], treat: s.x[2] }).safe).map((s) => ({ t: s.x[0] * 1000, rate: s.x[1], treat: s.x[2], npv: -s.f[0], risk: s.f[1], carbon: s.f[2] })).sort((a, b) => a.npv - b.npv);

  // ---------------------------------------------------------------- sustainability
  const carbonT = em0.total, intensity = meanBoe > 0 ? (carbonT * 1000) / meanBoe : 0, crf = capitalRecovery(X.r, lifeEff), refEm = options[argBest(options.map((o) => o.emis))];
  const mac = options.map((o) => { const dT = refEm.emis - o.emis, dC = (o.lcc - o.parts.carbon - (refEm.lcc - refEm.parts.carbon)) * crf; return { name: o.name, emis: o.emis, abate: dT, cost: dC, mac: dT > 1e-6 ? dC / dT : null }; }).sort((a, b) => (a.mac ?? 1e99) - (b.mac ?? 1e99));

  // ================================================================ results
  prog(0.93, 'Assembling results');
  const fx = q.fx, cur = q.currency || 'USD', pct = (x, d = 1) => (x === null || x === undefined ? '—' : +(100 * x).toFixed(d)), yrs = (x) => (x === null ? '—' : +x.toFixed(2));
  const hurdle = q.hurdle / 100, okNpv = npv0 >= q.minNpv * MM, okIrr = met.irr !== null && met.irr >= hurdle, okPay = met.payback !== null && met.payback <= q.maxPayback;
  const kpis = [
    { label: 'CAPEX', value: mUSD(capex.total), unit: 'M$', status: 'ok', help: `Evaluation-year money, ${sm0.name.toLowerCase()}; includes ${pct(q.contingency / 100, 0)} % contingency and ${pct(q.owners / 100, 0)} % owner's costs` },
    { label: 'OPEX (first production year)', value: mUSD(opexY1), unit: 'M$/y', status: 'ok', help: `${rd(opexY1 / Math.max(boeY1, 1), 2)} $/boe in real terms` },
    { label: 'NPV', value: mUSD(npv0), unit: 'M$', status: okNpv ? 'ok' : 'bad', help: `After tax at ${q.discount} % nominal, ${q.mid ? 'mid-year' : 'end-year'} discounting; threshold ${q.minNpv} M$` },
    { label: 'Risk-adjusted NPV', value: mUSD(riskedNpv), unit: 'M$', status: riskedNpv >= q.minNpv * MM ? 'ok' : 'bad', help: 'NPV less the expected cost of integrity failures (probability × consequence each year)' },
    { label: 'IRR', value: pct(met.irr), unit: '%/y', status: okIrr ? 'ok' : met.irr === null ? 'warn' : 'bad', help: `Hurdle rate ${q.hurdle} %${met.irrMultiple ? '; the cash flow changes sign more than once, so more than one IRR exists' : ''}` },
    { label: 'MIRR', value: pct(met.mirr), unit: '%/y', status: met.mirr !== null && met.mirr >= X.r ? 'ok' : 'warn', help: `Finance rate ${q.discount} %, reinvestment rate ${q.reinvest} %` },
    { label: 'Payback', value: yrs(met.payback), unit: 'y', status: okPay ? 'ok' : met.payback === null ? 'bad' : 'warn', help: `From the evaluation date, interpolated; limit ${q.maxPayback} y. Discounted payback ${yrs(met.discountedPayback)} y` },
    { label: 'Profitability index', value: rd(met.pi, 2), unit: '–', status: met.pi !== null && met.pi >= 1 ? 'ok' : 'bad', help: 'Present value of the net inflows per dollar of discounted CAPEX' },
    { label: 'Unit technical cost', value: rd(met.utc, 2), unit: '$/boe', status: met.utc !== null && met.utc < q.oilPrice * 0.6 ? 'ok' : 'warn', help: 'Discounted CAPEX + OPEX + abandonment over discounted production' },
    { label: 'Break-even oil price', value: rd(be.price, 1), unit: '$/bbl', status: be.price !== null && be.price < q.oilPrice * 0.8 ? 'ok' : be.price !== null && be.price < q.oilPrice ? 'warn' : 'bad', help: 'Flat real oil price at which the after-tax NPV is zero' },
    { label: 'NPV P10 / P50 / P90', value: `${mUSD(mc.p10, 0)} / ${mUSD(mc.p50, 0)} / ${mUSD(mc.p90, 0)}`, unit: 'M$', status: mc.p10 >= 0 ? 'ok' : 'warn', help: `${q.sampling === 'lhs' ? 'Latin-hypercube' : 'Monte Carlo'} sample of ${mc.n}; P10 is the low case` },
    { label: 'Probability of loss', value: pct(mc.probLoss), unit: '%', status: mc.probLoss < 0.1 ? 'ok' : mc.probLoss < 0.3 ? 'warn' : 'bad', help: 'Share of simulated outcomes with a negative NPV' },
    { label: `CVaR (${q.alpha} %)`, value: mUSD(mc.cvar), unit: 'M$', status: mc.cvar >= 0 ? 'ok' : 'warn', help: `Mean NPV of the worst ${rd(100 - q.alpha, 1)} % of outcomes; value at risk ${mUSD(mc.var)} M$` },
    { label: 'Expected annual loss', value: mUSD(eal, 2), unit: 'M$/y', status: eal < 0.01 * capex.total ? 'ok' : 'warn', help: 'Annual failure probability × consequence in the first year' },
    { label: 'Deferred and lost production', value: mUSD(deferredCost, 2), unit: 'M$/y', status: deferredCost < 0.03 * Math.max(cf.revenue[y1], 1) ? 'ok' : 'warn', help: `${rd(downH, 0)} h of downtime a year; ${q.deferFrac} % of the volume is recovered at the end of field life` },
    { label: 'Carbon intensity', value: rd(intensity, 1), unit: 'kgCO₂e/boe', status: intensity < 20 ? 'ok' : intensity < 40 ? 'warn' : 'bad', help: `${rd(carbonT, 0)} tCO₂e a year` },
    { label: 'Best flow-assurance strategy', value: noneSafe ? 'none passes the safety screen' : bestOpt.name, unit: '', status: noneSafe ? 'bad' : bestOpt.key === sel ? 'ok' : 'warn', help: `Lowest lifecycle cost among the options that meet the safety constraints (${mUSD(bestOpt.lcc)} M$); the case uses "${sm0.name}"` },
    { label: 'Optimum insulation', value: rd(insOpt.t, 0), unit: 'mm', status: Math.abs(insOpt.t - q.insMm) <= 15 ? 'ok' : 'warn', help: `Minimum lifecycle cost ${mUSD(insOpt.lcc)} M$; the case has ${q.insMm} mm` },
  ];

  // ---- engineering safety constraints of the case itself (from suites 3–6 when linked, else the kernel estimate)
  const flags = [];
  if (S.base.eros > 1) flags.push(`erosional velocity ratio ${rd(S.base.eros, 2)} exceeds 1.0 at the case rate`);
  if (S.base.pIn > q.mawp) flags.push(`inlet pressure ${rd(S.base.pIn, 0)} bara exceeds the allowable ${q.mawp} bara`);
  if (q.integUtil > 1) flags.push(`structural utilisation ${rd(q.integUtil, 2)} exceeds 1.0`);
  if (q.integViol > 0) flags.push(`${q.integViol} integrity code check(s) fail`);
  if (p0.haz[0] > q.maxPof) flags.push(`annual failure probability ${rd(p0.haz[0], 4)} exceeds the tolerable ${q.maxPof}`);
  if (q.severeSlug) flags.push('severe slugging is predicted and not suppressed');
  if (!selOpt.sm.safe) flags.push(...selOpt.sm.viol.map((x) => `case strategy: ${x}`));
  flags.forEach((f) => warn('bad', `Safety constraint: ${f}.`));
  if (noneSafe) warn('bad', 'No flow-assurance strategy meets the safety constraints; the ranking shows the least-cost option for reference only.');
  kpis.push({ label: 'Safety constraints', value: flags.length ? `${flags.length} violated` : 'all met', unit: '', status: flags.length ? 'bad' : 'ok', help: flags.length ? flags.join('; ') : 'Erosion, pressure, structural utilisation, failure probability, cooldown, blockage frequency and injection capacity are inside their limits' });

  // ---- warnings
  if (!okNpv) warn('bad', `NPV of ${mUSD(npv0)} M$ is below the ${q.minNpv} M$ threshold.`);
  if (met.irr === null) warn('warn', 'No internal rate of return exists: the cash flow never changes sign.');
  else if (!okIrr) warn('warn', `IRR of ${pct(met.irr)} % is below the ${q.hurdle} % hurdle rate.`);
  if (met.irrRoots.filter((x) => x >= 0).length > 1) warn('warn', `The cash flow has ${met.irrRoots.length} internal rates of return (${met.irrRoots.map((x) => pct(x)).join(' %, ')} %); use NPV and MIRR for the decision.`);
  if (met.payback === null) warn('bad', 'The investment is never paid back within the evaluation horizon.');
  if (limit0 !== null && limit0 < q.life) warn(limit0 < 1 ? 'bad' : 'info', limit0 < 1 ? 'Operating costs exceed net revenue from the first production year.' : `Economic limit reached after ${limit0} production years${q.stopAtLimit ? '; the evaluation stops there and abandonment is brought forward' : ' (the evaluation continues to the full project life)'}.`);
  if (mc.probLoss > 0.3) warn('bad', `Probability of a negative NPV is ${pct(mc.probLoss)} %.`);
  if (sm0.sub + q.hydMargin > 0 && (sel === 'risk' || sel === 'ldhi')) warn('warn', `Steady flow is ${rd(sm0.sub + q.hydMargin, 1)} °C inside the hydrate margin with the selected strategy.`);
  if (cf.lossPoolEnd > 1) warn('info', `${mUSD(cf.lossPoolEnd)} M$ of tax losses remain unused at the end of the evaluation.`);
  if (q.regime === 'psc' && cf.unrecoveredCost > 1) warn('warn', `${mUSD(cf.unrecoveredCost)} M$ of costs are never recovered under the cost-oil cap.`);
  if (ahpR.cr > 0.1) warn('warn', `The pairwise comparison matrix is inconsistent (consistency ratio ${rd(ahpR.cr, 2)} > 0.10); review the judgements.`);
  if (q.remLife < lifeEff) warn('warn', `Remaining asset life (${q.remLife} y) is shorter than the evaluation horizon (${lifeEff} y): wear-out failures dominate the late years.`);
  if (S.kCd !== 1) warn('info', `Cooldown model scaled by ${rd(S.kCd, 2)} to match the published cooldown time of ${q.cooldownBase} h.`);

  // ---- recommendations (decision-support voice)
  recs.push(okNpv && okIrr ? `Sanction case holds: NPV ${mUSD(npv0)} M$, IRR ${pct(met.irr)} %, payback ${yrs(met.payback)} y, break-even oil price ${rd(be.price, 1)} $/bbl against ${q.oilPrice} $/bbl assumed.` : `The project does not meet the investment criteria (NPV ${mUSD(npv0)} M$, IRR ${pct(met.irr)} %); it needs an oil price above ${rd(be.price, 1)} $/bbl or ${be.prod === null ? 'a larger resource' : pct(be.prod - 1, 0) + ' % more production'} to break even.`);
  if (flags.length) recs.push(`Resolve the safety findings before using the economics: ${flags.join('; ')}.`);
  recs.push(noneSafe ? `No flow-assurance strategy meets the safety limits (cooldown ≥ ${q.minCooldown} h, blockage ≤ ${q.maxBlockFreq} per year, injection ≤ ${q.maxInject} m³/d); the least-cost one, "${bestOpt.name}", fails on: ${bestOpt.sm.viol.join('; ')}.` : bestOpt.key === sel ? `Keep "${sm0.name}": it has the lowest lifecycle cost (${mUSD(bestOpt.lcc)} M$), ${mUSD(byLcc[1].lcc - bestOpt.lcc)} M$ below "${byLcc[1].name}"${byLcc[1].sm.safe ? '' : ' (which fails the safety screen)'}.` : !selOpt.sm.safe ? `Replace "${sm0.name}", which fails the safety screen (${selOpt.sm.viol.join('; ')}), with "${bestOpt.name}": lifecycle cost ${mUSD(selOpt.lcc)} → ${mUSD(bestOpt.lcc)} M$, NPV ${bestOpt.npv >= selOpt.npv ? '+' : ''}${mUSD(bestOpt.npv - selOpt.npv)} M$.` : `Switch from "${sm0.name}" to "${bestOpt.name}": lifecycle cost falls by ${mUSD(selOpt.lcc - bestOpt.lcc)} M$ (${mUSD(selOpt.lcc)} → ${mUSD(bestOpt.lcc)} M$) and NPV changes by ${mUSD(bestOpt.npv - selOpt.npv)} M$.`);
  { const cheaper = options.filter((o) => !o.sm.safe && o.lcc < bestOpt.lcc && !noneSafe); if (cheaper.length) recs.push(`Not recommended although cheaper: ${cheaper.map((o) => `"${o.name}" (${mUSD(o.lcc)} M$; ${o.sm.viol[0]})`).join(', ')}.`); }
  { const thin = insCurve.find((c) => c.t >= Math.min(20, insOpt.t)) || insBare, dC = insOpt.capex - thin.capex, dO = thin.lcc - thin.capex - (insOpt.lcc - insOpt.capex);
    recs.push(`Wet insulation at ${rd(insOpt.t, 0)} mm minimises lifecycle cost${insSafe.length ? ` among the thicknesses that keep ${q.minCooldown} h of cooldown and a tolerable blockage frequency` : ' (no thickness meets the cooldown and blockage limits)'}: against ${rd(thin.t, 0)} mm, ${dC >= 0 ? '+' : ''}${mUSD(dC)} M$ of CAPEX ${dO >= 0 ? 'saves' : 'costs'} ${mUSD(Math.abs(dO))} M$ of chemicals, deferred production and blockage risk over ${lifeEff} y; the case thickness of ${q.insMm} mm is ${mUSD(insSel.lcc - insOpt.lcc)} M$ from the optimum${insFree !== insOpt ? ` (ignoring the safety limits the cost minimum would be ${rd(insFree.t, 0)} mm)` : ''}.`); }
  if (dBest && dBase) recs.push(dBest.mult === 1 ? `The ${rd(S.id * 1000, 0)} mm bore is the best of the ${diam.length} diameters studied (NPV ${mUSD(dBest.npv)} M$).` : `A ${rd(dBest.d * 1000, 0)} mm bore raises NPV by ${mUSD(dBest.npv - (dBase.npv ?? 0))} M$ against ${rd(S.id * 1000, 0)} mm (rate ${pct(dBest.rate, 0)} % of plan, slug catcher ${rd(dBest.slugVol, 0)} m³, erosional ratio ${rd(dBest.eros * dBest.rate, 2)}).`);
  recs.push(`Operate at ${pct(opBest.x[0], 0)} % of the case rate: NPV ${mUSD(opBest.npv)} M$ (${opBest.npv >= npv0 ? '+' : ''}${mUSD(opBest.npv - npv0)} M$ against the case rate)${active.length ? `, limited by ${active.join(' and ')}` : ''}.`);
  recs.push(!rbiOk.length ? `No inspection interval keeps the annual failure probability below ${q.maxPof}: the lowest reachable is ${rd(rbiBest.hMax, 4)} at ${rd(rbiBest.t, 1)} y — repair or replace rather than inspect.` : rbiBest.t < Math.min(lifeEff, 20) - 0.6 ? `Inspect every ${rd(rbiBest.t, 1)} y: lifecycle inspection plus expected failure cost is ${mUSD(rbiBest.total)} M$ against ${mUSD(rbiNow.total)} M$ at the current ${q.inspInterval} y interval${rbiFree !== rbiBest ? `; the cost minimum alone would be ${rd(rbiFree.t, 1)} y but breaks the failure-probability limit of ${q.maxPof} per year` : ''}.` : `Inspection does not pay for itself on failure-cost grounds with these inputs (expected failure cost ${mUSD(rbiNone.failure)} M$ over the life); keep the regulatory minimum interval.`);
  recs.push(`Material selection: "${matBest.name}" has the lowest lifecycle cost (${mUSD(matBest.lcc)} M$ relative to the case design) at a corrosion rate of ${rd(matBest.rate, 2)} mm/y${mats.some((m) => !m.practical) ? `; ${mats.filter((m) => !m.practical).map((m) => `"${m.name}"`).join(', ')} is excluded because the allowance would exceed 10 mm` : ''}.`);
  if (q.waxRate > 0) recs.push(`Pig every ${rd(pigBest.tau, 0)} d (${rd(pigBest.runs, 1)} runs a year, deposit kept below ${rd(pigBest.dEnd, 1)} mm): ${mUSD(pigBest.total, 2)} M$/y against ${mUSD(pigNow.total, 2)} M$/y at the present ${rd(pigNow.tau, 0)} d interval.`);
  recs.push(voi.evii > q.testCost * MM ? `Appraise before sanction: imperfect information is worth ${mUSD(voi.evii)} M$ against a cost of ${q.testCost} M$ (perfect information would be worth ${mUSD(voi.evpi)} M$).` : `Do not pay ${q.testCost} M$ for further appraisal: it is worth only ${mUSD(voi.evii)} M$ at ${q.testRel} % reliability (perfect information ${mUSD(voi.evpi)} M$).`);
  recs.push(defer - Math.max(npv0, 0) > 0.02 * Math.abs(npv0) + 1e5 ? `Waiting has value: the option to defer for up to ${q.optYears} y is worth ${mUSD(defer)} M$ against ${mUSD(Math.max(npv0, 0))} M$ for investing now.` : `Invest now: deferral adds no value (option ${mUSD(defer)} M$ ≈ static NPV ${mUSD(Math.max(npv0, 0))} M$) because ${q.optYield} %/y of project value is lost while waiting.`);
  recs.push(pMean.t === pCvar.t ? `Under uncertainty the insulation choice is robust: ${rd(pMean.t, 0)} mm maximises both the expected NPV and the ${q.alpha} % CVaR.` : `Under uncertainty ${rd(pMean.t, 0)} mm maximises the expected NPV, while a risk-averse owner would choose ${rd(pCvar.t, 0)} mm (best ${q.alpha} % CVaR).`);
  recs.push(spare.saving > 0 ? `Hold the critical spare: it saves ${mUSD(spare.saving, 2)} M$/y net of ${mUSD(spare.hold, 2)} M$/y holding cost.` : `Do not hold the critical spare: holding cost ${mUSD(spare.hold, 2)} M$/y exceeds the expected downtime saving unless the failure rate rises above ${rd(spare.beRate, 3)} per year.`);
  const macBest = mac.find((x) => x.mac !== null && pool.some((o) => o.name === x.name));
  recs.push(`Carbon: ${rd(carbonT, 0)} tCO₂e/y (${rd(intensity, 1)} kg/boe), costing ${mUSD(carbonT * q.carbonPrice, 2)} M$/y at ${q.carbonPrice} $/t; the largest source is ${em0.items.slice().sort((a, b) => b.t - a.t)[0].source.toLowerCase()}.${macBest ? ` Cheapest abatement among the strategies: "${macBest.name}" at ${rd(macBest.mac, 0)} $/tCO₂e relative to "${refEm.name}"${macBest.mac < q.carbonPrice ? ', below the carbon price' : ''}.` : ''}`);
  if (em0.items[3].t + em0.items[4].t > 0.3 * carbonT) recs.push('Flaring is more than 30 % of the footprint: recover blowdown gas and reduce routine flaring before investing elsewhere.');

  // ---- plots
  const years = cf.year, M = (a) => a.map((x) => x / MM), cdfN = Math.min(200, mc.n), cdfI = Array.from({ length: cdfN }, (_, i) => Math.round((i * (mc.n - 1)) / Math.max(cdfN - 1, 1))), hist = histogram(mc.y.map((y) => y / MM), 30);
  const gKeys = Object.keys(capex.groups), opexCat = [...p0.fixedItems.map(([n, x]) => [n, x]), ['Tariff and production chemicals', p0.opexVarBoe * boeY1], ['Produced-water handling', p0.waterCost * cf.water[y1]], ['Shutdown preservation', p0.opexDown], ['Hydrate / wax remediation (expected)', p0.opexBlock]];
  plots.push(
    { type: 'bar', title: 'Annual cash flow (money of the day)', ylabel: 'M$', categories: years.map(String), stacked: true, series: [{ name: 'Revenue', values: M(cf.revenue) }, { name: 'Royalty, tax and state share', values: M(cf.royalty.map((x, k) => -(x + cf.tax[k] + cf.govShare[k]))) }, { name: 'OPEX, carbon and risk cost', values: M(cf.opex.map((x, k) => -(x + cf.carbon[k] + (cfR.risk[k] || 0)))) }, { name: 'CAPEX, working capital, abandonment', values: M(cf.capex.map((x, k) => -(x + cf.dwc[k] + cf.abandon[k]))) }] },
    { type: 'line', title: 'Cumulative cash flow', xlabel: 'Year from the evaluation date', ylabel: 'M$', zeroY: true, series: [{ name: 'Cumulative', x: years, y: M(cf.cum) }, { name: 'Discounted cumulative', x: years, y: M(cf.cumDcf) }, { name: 'Discounted cumulative, risk-adjusted', x: years, y: M(cfR.cumDcf), dash: true }], hlines: [{ y: 0, label: 'break-even' }], vlines: met.payback !== null && met.payback > 0 ? [{ x: met.payback, label: 'payback' }] : [] },
    { type: 'bar', title: 'CAPEX build-up', ylabel: 'M$', categories: gKeys, series: [{ name: 'CAPEX', values: gKeys.map((g) => capex.groups[g] / MM) }] },
    { type: 'bar', title: 'OPEX build-up (first production year, real terms)', ylabel: 'M$/y', categories: opexCat.map((c) => c[0]), series: [{ name: 'OPEX', values: opexCat.map((c) => c[1] / MM) }] },
    { type: 'line', title: 'Production and revenue profile', xlabel: 'Year from the evaluation date', ylabel: 'Thousand boe/d · M$/y ÷ 10', series: [{ name: 'Produced (kboe/d)', x: years, y: cf.boe.map((b) => b / 365 / 1000), mode: 'step' }, { name: 'Potential (kboe/d)', x: years, y: cf.potBoe.map((b) => b / 365 / 1000), mode: 'step', dash: true }, { name: 'Revenue (M$/y ÷ 10)', x: years, y: cf.revenue.map((x) => x / MM / 10) }] },
    { type: 'line', title: 'NPV distribution (histogram)', xlabel: 'Risk-adjusted NPV (M$)', ylabel: 'Count', series: [{ name: 'Outcomes', x: hist.centers, y: hist.counts, mode: 'step' }], vlines: [{ x: mc.p10 / MM, label: 'P10' }, { x: mc.p50 / MM, label: 'P50' }, { x: mc.p90 / MM, label: 'P90' }] },
    { type: 'line', title: 'NPV cumulative distribution', xlabel: 'Risk-adjusted NPV (M$)', ylabel: 'Cumulative probability', ymin: 0, ymax: 1, series: [{ name: 'CDF', x: cdfI.map((i) => mc.sorted[i] / MM), y: cdfI.map((i) => (i + 0.5) / mc.n) }], vlines: [{ x: mc.p10 / MM, label: 'P10' }, { x: mc.p50 / MM, label: 'P50' }, { x: mc.p90 / MM, label: 'P90' }, { x: 0, label: 'loss' }] },
    { type: 'line', title: 'Convergence of the mean NPV with sample size', xlabel: 'Samples', ylabel: 'Mean NPV (M$)', logx: true, series: [{ name: 'Running mean', x: conv.map((c) => c.n), y: conv.map((c) => c.mean / MM) }, { name: '+2 standard errors', x: conv.map((c) => c.n), y: conv.map((c) => (mc.mean + (2 * mc.sd) / Math.sqrt(c.n)) / MM), dash: true }, { name: '−2 standard errors', x: conv.map((c) => c.n), y: conv.map((c) => (mc.mean - (2 * mc.sd) / Math.sqrt(c.n)) / MM), dash: true }] },
    { type: 'bar', title: 'Tornado: NPV swing for each input at its P10 and P90', ylabel: 'Change in NPV (M$)', categories: torn.map((t) => t.name), series: [{ name: 'Input at P10', values: torn.map((t) => (t.low - t.base) / MM) }, { name: 'Input at P90', values: torn.map((t) => (t.high - t.base) / MM) }] },
    { type: 'line', title: 'Spider plot', xlabel: 'Input as a percentage of its base value', ylabel: 'Risk-adjusted NPV (M$)', series: spider.map((s) => ({ name: s.name, x: spiderF.map((f) => f * 100), y: s.y })) },
    { type: 'bar', title: 'Lifecycle cost of the flow-assurance strategies', ylabel: 'M$ (present value)', stacked: true, categories: byLcc.map((o) => o.name), series: [['CAPEX', 'capex'], ['Chemicals and energy', 'chem'], ['Pigging', 'pig'], ['Shutdown preservation and deferral', 'shut'], ['Expected blockage cost', 'block'], ['Carbon cost', 'carbon']].map(([name, k]) => ({ name, values: byLcc.map((o) => o.parts[k] / MM) })) },
    { type: 'line', title: 'Insulation thickness: lifecycle cost', xlabel: 'Wet insulation thickness (mm)', ylabel: 'M$ (present value)', ymin: 0, ymax: Math.max(3 * insOpt.lcc / MM, 10), series: [{ name: 'Lifecycle cost', x: tGrid, y: insCurve.map((c) => c.lcc / MM) }, { name: 'Insulation CAPEX', x: tGrid, y: insCurve.map((c) => c.capex / MM), dash: true }, { name: 'Chemicals, pigging, shutdowns', x: tGrid, y: insCurve.map((c) => (c.opex + c.defer) / MM), dash: true }, { name: 'Expected blockage cost', x: tGrid, y: insCurve.map((c) => c.block / MM), dash: true }, { name: 'Kernel solutions', x: S.sweep.map((s) => s.tMm), y: S.sweep.map((s) => lin(tGrid, insCurve.map((c) => c.lcc / MM), s.tMm)), mode: 'points' }], vlines: [{ x: insOpt.t, label: 'optimum' }], note: 'Points mark thicknesses solved on the flow kernel; the curve interpolates the arrival temperature and hydrate margin between them.' },
    { type: 'line', title: 'Insulation thickness: arrival temperature and cooldown time', xlabel: 'Wet insulation thickness (mm)', ylabel: '°C · h', series: [{ name: 'Arrival temperature (°C)', x: tGrid, y: insCurve.map((c) => c.tArr) }, { name: 'Cooldown to hydrate temperature (h)', x: tGrid, y: insCurve.map((c) => Math.min(c.cd, 200)) }, { name: 'Subcooling into the hydrate region (°C)', x: tGrid, y: insCurve.map((c) => c.sub), dash: true }], hlines: [{ y: q.wat, label: 'wax appearance' }] },
    { type: 'line', title: 'Pipe diameter optimisation', xlabel: 'Inner diameter (mm)', ylabel: 'M$ · %', series: [{ name: 'NPV (M$)', x: diam.filter((d) => d.npv !== null).map((d) => d.d * 1000), y: diam.filter((d) => d.npv !== null).map((d) => d.npv / MM), mode: 'both' }, { name: 'CAPEX (M$)', x: diam.filter((d) => d.npv !== null).map((d) => d.d * 1000), y: diam.filter((d) => d.npv !== null).map((d) => d.capexTotal / MM), mode: 'both', dash: true }, { name: 'Deliverable rate (% of plan)', x: diam.filter((d) => d.npv !== null).map((d) => d.d * 1000), y: diam.filter((d) => d.npv !== null).map((d) => d.rate * 100), mode: 'both', dash: true }], vlines: dBest ? [{ x: dBest.d * 1000, label: 'best' }] : [] },
    { type: 'line', title: 'Economic operating envelope: NPV against operating rate', xlabel: 'Rate (% of the case rate)', ylabel: 'NPV (M$)', series: [{ name: 'Inside the constraints', x: rateCurve.filter((c) => c.ok).map((c) => c.r * 100), y: rateCurve.filter((c) => c.ok).map((c) => c.npv / MM) }, { name: 'Outside the constraints', x: rateCurve.filter((c) => !c.ok).map((c) => c.r * 100), y: rateCurve.filter((c) => !c.ok).map((c) => c.npv / MM), mode: 'points' }], vlines: [{ x: opBest.x[0] * 100, label: 'optimum' }] },
    { type: 'line', title: 'Pareto front: value against risk and carbon', xlabel: 'NPV (M$)', ylabel: 'M$/y · ktCO₂e/y', series: [{ name: 'Expected annual loss (M$/y)', x: front.map((f) => f.npv), y: front.map((f) => f.risk), mode: 'points' }, { name: 'Emissions (ktCO₂e/y)', x: front.map((f) => f.npv), y: front.map((f) => f.carbon / 1000), mode: 'points' }], note: `${front.length} non-dominated designs (insulation thickness, operating rate, share of long shutdowns preserved).` },
    { type: 'bar', title: 'Decision tree: expected monetary value by branch', ylabel: 'EMV (M$)', categories: [...tree.branches[0].node.branches.map((b) => b.name), 'Appraise first (net of cost)'], series: [{ name: 'EMV', values: [...tree.branches[0].node.branches.map((b) => b.emv / MM), tree.branches[1].emv / MM] }] },
    { type: 'line', title: 'Real option to defer: value against volatility', xlabel: 'Volatility of project value (%/y)', ylabel: 'Option value (M$)', series: [{ name: 'Binomial lattice, American', x: volCurve.map((c) => c.s), y: volCurve.map((c) => c.am / MM) }, { name: 'Binomial lattice, European', x: volCurve.map((c) => c.s), y: volCurve.map((c) => c.lat / MM), mode: 'points' }, { name: 'Black–Scholes', x: volCurve.map((c) => c.s), y: volCurve.map((c) => c.bs / MM), dash: true }], hlines: [{ y: Math.max(npv0, 0) / MM, label: 'invest now' }] },
    { type: 'line', title: 'Risk-based inspection interval', xlabel: 'Inspection interval (y)', ylabel: 'M$ (present value over the life)', series: [{ name: 'Total', x: rbi.map((r) => r.t), y: rbi.map((r) => r.total / MM) }, { name: 'Inspection cost', x: rbi.map((r) => r.t), y: rbi.map((r) => r.inspection / MM), dash: true }, { name: 'Expected failure cost', x: rbi.map((r) => r.t), y: rbi.map((r) => r.failure / MM), dash: true }], vlines: [{ x: rbiBest.t, label: 'optimum' }] },
    { type: 'line', title: 'Pigging interval: cost of intervention against probability × consequence', xlabel: 'Interval between pig runs (d)', ylabel: 'M$/y', logx: true, series: [{ name: 'Total', x: pigCurve.map((c) => c.tau), y: pigCurve.map((c) => Math.min(c.total / MM, 50 * Math.max(pigBest.total / MM, 0.02))) }, { name: 'Pigging', x: pigCurve.map((c) => c.tau), y: pigCurve.map((c) => c.pig / MM), dash: true }, { name: 'Stuck pig / restriction risk', x: pigCurve.map((c) => c.tau), y: pigCurve.map((c) => Math.min((c.stuck + c.loss) / MM, 50 * Math.max(pigBest.total / MM, 0.02))), dash: true }], vlines: [{ x: pigBest.tau, label: 'optimum' }] },
    { type: 'line', title: 'Insulation decision under uncertainty', xlabel: 'Wet insulation thickness (mm)', ylabel: 'M$', series: [{ name: 'Expected NPV', x: probIns.map((c) => c.t), y: probIns.map((c) => c.mean / MM), mode: 'both' }, { name: `CVaR ${q.alpha} %`, x: probIns.map((c) => c.t), y: probIns.map((c) => c.cvar / MM), mode: 'both' }, { name: 'P10', x: probIns.map((c) => c.t), y: probIns.map((c) => c.p10 / MM), dash: true }] },
    { type: 'bar', title: 'Emission inventory', ylabel: 'tCO₂e/y', categories: em0.items.map((i) => i.source), series: [{ name: 'Emissions', values: em0.items.map((i) => i.t) }] },
    { type: 'bar', title: 'Marginal abatement cost of the strategies', ylabel: '$/tCO₂e', categories: mac.filter((x) => x.mac !== null).map((x) => x.name), series: [{ name: `Relative to "${refEm.name}"`, values: mac.filter((x) => x.mac !== null).map((x) => x.mac) }, { name: 'Carbon price', values: mac.filter((x) => x.mac !== null).map(() => q.carbonPrice) }] },
  );

  // ---- tables
  const col = (a, d = 2) => a.map((x) => rd(x / MM, d));
  tables.push({ title: 'Cash-flow statement (M$, money of the day)', columns: ['Year', 'Oil (kbbl)', 'Gas sold (GJ ×10³)', 'Oil price ($/bbl)', 'Revenue', 'Royalty', 'State profit share', 'OPEX', 'Carbon cost', 'Operating cash flow', 'Depreciation', 'Taxable income', 'Tax', 'After-tax cash flow', 'CAPEX', 'Δ working capital', 'Abandonment', 'Free cash flow', 'Cumulative', 'Discount factor', 'Discounted', 'Cumulative discounted', 'Free cash flow (real)'],
    rows: years.map((k) => [k, rd(cf.oil[k] / 1000, 0), rd((cf.gas[k] * GJ_PER_MMBTU) / 1000, 0), rd(cf.price[k], 2), ...[cf.revenue, cf.royalty, cf.govShare, cf.opex, cf.carbon, cf.ocf, cf.depreciation, cf.taxable, cf.tax, cf.atcf, cf.capex, cf.dwc, cf.abandon, cf.fcf, cf.cum].map((a) => rd(a[k] / MM, 2)), rd(cf.df[k], 4), rd(cf.dcf[k] / MM, 2), rd(cf.cumDcf[k] / MM, 2), rd(cf.real[k] / MM, 2)]),
    note: `Year 0 is the evaluation date (${q.evalYear}); CAPEX is phased over ${nCon} year(s); production starts in year ${nCon}. ${q.regime === 'psc' ? `Production-sharing contract: cost oil capped at ${q.costOilCap} % of net revenue, contractor profit share ${q.profitSplit} %.` : `Tax and royalty regime with ${q.deprMethod === 'sl' ? 'straight-line' : q.deprMethod === 'db' ? 'declining-balance' : 'units-of-production'} depreciation and loss carry-forward.`} Receipts are positive and payments negative in the free cash flow.` });
  tables.push({ title: 'Investment metrics', columns: ['Metric', 'Value', 'Unit', 'Basis'], rows: [
    ['Net present value', mUSD(npv0, 2), 'M$', `${q.discount} % nominal, ${q.mid ? 'mid-year' : 'end-year'}`], ['NPV in real terms', mUSD(met.realNpv, 2), 'M$', `deflated at ${q.inflation} %, real rate ${pct(met.realRate, 2)} % (Fisher)`], ['Risk-adjusted NPV', mUSD(riskedNpv, 2), 'M$', 'less expected failure cost'],
    ['Internal rate of return', pct(met.irr, 2), '%/y', met.irrRoots.length > 1 ? `${met.irrRoots.length} roots: ${met.irrRoots.map((x) => pct(x, 1)).join(', ')} %` : 'unique'], ['Modified IRR', pct(met.mirr, 2), '%/y', `finance ${q.discount} %, reinvest ${q.reinvest} %`], ['Return on investment', pct(met.roi, 1), '%', 'undiscounted net cash / CAPEX'], ['Return on capital employed', pct(met.roce, 1), '%/y', 'after-tax profit / average book value'],
    ['Profitability index', rd(met.pi, 3), '–', '1 + NPV / PV(CAPEX)'], ['Simple payback', yrs(met.payback), 'y', 'interpolated, from year 0'], ['Discounted payback', yrs(met.discountedPayback), 'y', 'interpolated, from year 0'], ['Equivalent annual value', mUSD(met.eav, 2), 'M$/y', `${K - 1} y`], ['Equivalent annual cost', mUSD(met.eac, 2), 'M$/y', 'CAPEX + OPEX + abandonment'],
    ['Unit technical cost', rd(met.utc, 2), '$/boe', 'discounted cost / discounted production'], ['Lifting cost', rd(met.liftingCost, 2), '$/boe', 'undiscounted OPEX / production'], ['Government take', pct(met.governmentTake, 1), '%', 'royalty, tax and state share of pre-take cash'], ['Economic limit', limit0 === null ? 'beyond the horizon' : limit0, 'production year', 'first year with negative operating margin'],
    ['Equity IRR (geared)', pct(met.equityIrr, 2), '%/y', q.gearing > 0 ? `${q.gearing} % debt at ${q.loanRate} % over ${q.loanTenor} y` : 'no debt in this case'], ['Equity NPV (geared)', q.gearing > 0 ? mUSD(met.equityNpv, 2) : '—', 'M$', 'cash flow to equity at the discount rate'],
    [`NPV in ${cur}`, rd((npv0 / MM) * fx, 1), `M ${cur}`, `${fx} ${cur} per US$`], [`CAPEX in ${cur}`, rd((capex.total / MM) * fx, 1), `M ${cur}`, `${fx} ${cur} per US$`],
  ] });
  tables.push({ title: 'Break-evens (root of the NPV)', columns: ['Driver', 'Break-even', 'Unit', 'Case value', 'Margin'], rows: [
    ['Oil price', rd(be.price, 2), '$/bbl', q.oilPrice, be.price === null ? '—' : `${pct(1 - be.price / q.oilPrice, 0)} % below the case`], ['Production', be.prod === null ? '—' : rd(be.prod * q.qOil, 0), 'Sm³/d oil (initial)', q.qOil, be.prod === null ? '—' : `${pct(1 - be.prod, 0)} % below the case`], ['CAPEX', be.capex === null ? '—' : mUSD(be.capex * capex.total), 'M$', mUSD(capex.total), be.capex === null ? '—' : `${pct(be.capex - 1, 0)} % overrun tolerated`],
    ['Continuous inhibitor dosage', rd(be.inhib, 1), `m³/d ${q.inhibitor}`, rd(q.inhibRate, 2), be.inhib === null ? 'no break-even in range' : 'additional rate that takes the NPV to zero'], ['Blockage / intervention frequency', rd(be.block, 2), 'events/y', rd(sm0.events, 3), be.block === null ? 'no break-even in range' : `at ${mUSD(X.blockCost)} M$ per event`], ['Uptime', be.uptime === null ? '—' : pct(be.uptime, 1), '%', pct(cf.uptime, 1), be.uptime === null ? 'NPV positive at any uptime above 5 %' : 'minimum availability'],
    ['Discount rate (= IRR)', pct(met.irr, 2), '%/y', q.discount, '—'], ['Economic cut-off', limit0 === null ? `> ${q.life}` : limit0, 'production year', q.life, q.stopAtLimit ? 'evaluation stops at the limit' : 'evaluation runs the full life'],
  ] });
  tables.push({ title: 'CAPEX build-up', columns: ['Group', 'Item', 'Basis', 'Cost (M$)'], rows: [...capex.items.map((i) => [i.group, i.item, i.basis, rd(i.cost / MM, 2)]), ['Total', 'CAPEX', `cost index ${q.costIndexEval}/${q.costIndexBase}, location factor ${q.locFactor}`, rd(capex.total / MM, 2)], ['Check', 'Topsides by bare-module factors', 'purchased cost × item factors', rd(capex.moduleCost / MM, 2)], ['Check', 'Topsides by one Lang factor', `purchased cost ${rd(capex.purchased / MM, 2)} M$ × ${q.langFactor}`, rd(capex.langCost / MM, 2)], ['Check', 'Installed line by the parametric relationship', `${q.cerCoef} M$/km, exponent ${q.cerExp}, learning ${q.learnRate} %`, rd(pipelineCER({ dIn: S.id / 0.0254, lengthKm: q.lineLen / 1000, depth: q.depth, cerCoef: q.cerCoef, cerExp: q.cerExp, layFactor: q.layFactor, learnRate: q.learnRate / 100, vesselRate: q.vesselRate, layRate: q.layRate, mobCost: q.mobCost, depthCoef: q.depthCoef }) * capex.escalation, 2)]], note: `Cost basis: 2023 US$ (index ${q.costIndexBase}), escalated to ${q.evalYear} (index ${q.costIndexEval}). ${rd(capex.steelT, 0)} t of steel, ${rd(capex.layDays, 1)} vessel days. Order-of-magnitude (class 4–5) estimate.` });
  tables.push({ title: 'OPEX build-up (first production year, real terms)', columns: ['Category', 'M$/y', '$/boe', 'Driver'], rows: [...opexCat.map(([n, x]) => [n, rd(x / MM, 3), rd(x / Math.max(boeY1, 1), 2), '']), ['Total', rd(sum(opexCat.map((c) => c[1])) / MM, 3), rd(sum(opexCat.map((c) => c[1])) / Math.max(boeY1, 1), 2), ''], ['Carbon cost', rd((carbonT * q.carbonPrice) / MM, 3), rd((carbonT * q.carbonPrice) / Math.max(boeY1, 1), 2), `${rd(carbonT, 0)} t at ${q.carbonPrice} $/t`]].map((r, i) => { const drv = [`fixed`, `${q.maintPct} % of CAPEX`, `${q.insurPct} % of CAPEX`, `${q.inspCost} M$ every ${q.inspInterval} y`, 'monitoring and chemicals', sm0.note, `${rd(q.pumpKW + q.compKW, 0)} kW at ${rd(X.enPrice, 3)} $/kWh`, `${rd(sm0.pigRuns, 1)} runs at ${q.pigCost} k$`, `${q.tariff} + ${q.chemOther} $/boe`, `${q.waterCost} $/m³`, `${rd(sm0.nLong, 2)} long shutdowns a year`, `${rd(sm0.events, 3)} events/y × ${mUSD(X.blockCost)} M$`][i]; return drv !== undefined ? [r[0], r[1], r[2], drv] : r; }) });
  tables.push({ title: 'Production economics', columns: ['Quantity', 'Value', 'Unit'], rows: [
    ['Recoverable volume produced', rd(sum(cf.boe) / 1e6, 2), 'million boe'], ['Potential volume (no downtime)', rd(sum(cf.potBoe) / 1e6, 2), 'million boe'], ['Production efficiency', pct(sum(cf.boe) / Math.max(sum(cf.potBoe), 1), 2), '%'], ['Uptime used', pct(cf.uptime, 2), '%'],
    ['Deferred volume (recovered in the last year)', rd(sum(cf.potBoe) * (1 - cf.uptime) * p0.deferFrac / 1e6, 3), 'million boe'], ['Lost volume', rd(sum(cf.lostBoe) / 1e6, 3), 'million boe'], ['Value lost by deferral', mUSD(downH * X.marginH * p0.deferFrac * (X.lossFrac - (1 - p0.deferFrac)) / Math.max(p0.deferFrac, 1e-9), 3), 'M$/y'], ['Value of lost production', mUSD(downH * X.marginH * (1 - p0.deferFrac), 3), 'M$/y'],
    ['Net revenue per hour of production', rd(X.marginH / 1000, 1), 'k$/h'], ['Decline model', q.declineType === 'exp' ? 'exponential' : q.declineType === 'har' ? 'harmonic' : `hyperbolic, b = ${q.bHyp}`, `${q.Di} %/y after ${q.plateau} y plateau`], ['Water cut', `${q.wc0} → ${q.wcEnd}`, '%'], ['Total revenue', mUSD(sum(cf.revenue), 0), 'M$'],
  ] });
  tables.push({ title: 'Flow-assurance strategy ranking (lifecycle cost)', columns: ['Rank', 'Strategy', 'U (W/m²K)', 'Arrival T (°C)', 'Cooldown (h)', 'Blockage events/y', 'FA CAPEX (M$)', 'PV OPEX (M$)', 'PV blockage (M$)', 'Lifecycle cost (M$)', 'NPV (M$)', 'tCO₂e/y', 'MCDA score', 'TOPSIS', 'Safety screening', 'Basis'],
    rows: byLcc.map((o, i) => [o.sm.safe ? i + 1 : '—', o.name, rd(o.sm.U, 2), rd(o.sm.tArr, 1), o.sm.cooldownH >= 1e4 ? 'never' : rd(o.sm.cooldownH, 1), rd(o.sm.events, 3), mUSD(o.parts.capex), mUSD(o.parts.chem + o.parts.pig + o.parts.shut + o.parts.carbon), mUSD(o.parts.block), mUSD(o.lcc), mUSD(o.npv), rd(o.emis, 0), rd(o.score, 3), rd(o.closeness, 3), o.sm.safe ? 'passes' : 'EXCLUDED: ' + o.sm.viol.join('; '), o.sm.note]),
    note: `Lifecycle cost = flow-assurance CAPEX + PV(OPEX) + PV(expected blockage cost), before tax. Blockage frequency = long shutdowns × plug probability (${q.plugProb}) × residual exposure; one blockage costs ${mUSD(X.blockCost)} M$. MCDA winner: ${mcdaBest.name}; TOPSIS winner: ${topsisBest.name}; minimax-regret choice: ${regretBest.name}. Options that break a safety constraint are listed last and take no part in the multi-criteria ranking.` });
  tables.push({ title: 'Criteria weights (Analytic Hierarchy Process)', columns: ['Criterion', 'Weight', 'Direction'], rows: [...crit.map((c, i) => [c, rd(ahpR.weights[i], 4), benefit[i] ? 'maximise' : 'minimise']), ['λmax', rd(ahpR.lambdaMax, 4), ''], ['Consistency index', rd(ahpR.ci, 4), ''], ['Consistency ratio', rd(ahpR.cr, 4), ahpR.cr <= 0.1 ? 'acceptable (≤ 0.10)' : 'inconsistent']], note: 'Operability scores (1–10): ' + keys.map((k) => `${STRATEGIES[k].name} ${STRATEGIES[k].operability}`).join(', ') + '.' });
  tables.push({ title: 'Uncertainty: input distributions and sensitivities', columns: ['Input (multiplier)', 'Distribution', 'P10', 'Mean', 'P90', 'NPV at P10 (M$)', 'NPV at P90 (M$)', 'Swing (M$)', 'Sobol first-order', 'Sobol total', 'Std. regression coeff.'],
    rows: VARS.map((id, j) => { const t = torn.find((x) => x.name === VAR_LABEL[id]), li = live.indexOf(j); return [VAR_LABEL[id], ds[j].kind, rd(ds[j].inv(0.1), 3), rd(ds[j].mean, 3), rd(ds[j].inv(0.9), 3), mUSD(t.low), mUSD(t.high), mUSD(t.swing), li >= 0 ? rd(sob.first[li], 3) : 0, li >= 0 ? rd(sob.total[li], 3) : 0, li >= 0 && srcFit.src.length ? rd(srcFit.src[li], 3) : 0]; }),
    note: `Sobol indices by the Saltelli scheme with ${sob.evals} model runs (independent inputs; sampling error about ±${rd(1.5 / Math.sqrt(q.nSobol), 2)}). Regression R² = ${rd(srcFit.r2, 3)}. Sampling: ${q.sampling === 'lhs' ? 'Latin hypercube with Iman–Conover rank correlation' : 'Monte Carlo with a Gaussian copula'}${q.priceModel === 'static' ? '' : `, ${q.priceModel === 'gbm' ? 'geometric Brownian' : 'mean-reverting'} price paths`}${q.failEvents ? ', failures sampled as discrete events' : ''}.` });
  tables.push({ title: 'Uncertainty: NPV statistics and scenarios', columns: ['Quantity', 'Value (M$)', 'Note'], rows: [
    ['Mean', mUSD(mc.mean), `standard error ${mUSD(se, 2)} M$ with ${mc.n} samples`], ['Standard deviation', mUSD(mc.sd), `variance ${rd(mc.variance / MM / MM, 0)} (M$)²`], ['P10 (low)', mUSD(mc.p10), '10th percentile'], ['P50', mUSD(mc.p50), 'median'], ['P90 (high)', mUSD(mc.p90), '90th percentile'], ['Minimum / maximum', `${mUSD(mc.min, 0)} / ${mUSD(mc.max, 0)}`, ''],
    [`Value at risk (${q.alpha} %)`, mUSD(mc.var), `NPV exceeded with ${q.alpha} % probability`], [`Conditional value at risk (${q.alpha} %)`, mUSD(mc.cvar), `mean of the worst ${rd(100 - q.alpha, 1)} %`], ['Probability of loss', `${pct(mc.probLoss, 2)} %`, 'P(NPV < 0)'], ['Certainty equivalent', mUSD(ce), `exponential utility, risk tolerance ${q.riskTol} M$; risk premium ${mUSD(mc.mean - ce)} M$`],
    ...scen.map((s) => [`Scenario: ${s.name}`, mUSD(s.npv), `weight ${rd(s.weight / wSum, 3)}; price ×${s.price}, production ×${s.prod}, CAPEX ×${s.capex}, OPEX ×${s.opex}`]), ['Scenario-weighted NPV', mUSD(scenEv), 'weights normalised to 1'],
  ] });
  tables.push({ title: 'Decision analysis', columns: ['Item', 'Value', 'Unit', 'Note'], rows: [
    ['Best action without further information', alts[voi.best].name, '', `EMV ${mUSD(voi.emv)} M$`], ['Expected value with perfect information', mUSD(voi.evWithPI), 'M$', ''], ['Value of perfect information', mUSD(voi.evpi, 2), 'M$', 'upper bound on any appraisal spend'], ['Value of imperfect information', mUSD(voi.evii, 2), 'M$', `appraisal reliability ${q.testRel} % (Bayesian revision)`], ['Appraisal cost', q.testCost, 'M$', voi.evii > q.testCost * MM ? 'worth buying' : 'not worth buying'],
    ...scen.map((s, k) => [`If the appraisal indicates "${s.name}"`, alts[voi.signalAction[k]].name, '', `probability ${rd(voi.signalProb[k], 3)}; posterior ${voi.posterior[k].map((x) => rd(x, 2)).join(' / ')}`]),
    ['Decision tree root', tree.choice, '', `EMV ${mUSD(tree.emv)} M$`], ['Minimax-regret strategy', regretBest.name, '', `largest regret ${mUSD(regret.maxRegret[regret.best])} M$`],
    ['Option to defer (American, lattice)', mUSD(defer, 2), 'M$', `${q.nLattice} steps, σ ${q.optVol} %, ${q.optYears} y, value leakage ${q.optYield} %/y`], ['Option to defer (European, lattice)', mUSD(deferEu, 2), 'M$', `Black–Scholes ${mUSD(bs, 2)} M$`], ['Value of waiting', mUSD(defer - Math.max(npv0, 0), 2), 'M$', 'option − invest-now NPV'], ['Option to expand', mUSD(expand, 2), 'M$', `+${q.expandFrac} % for ${q.expandCost} M$ at year ${q.optYears}`], ['Option to abandon', mUSD(abandonOpt, 2), 'M$', `salvage ${q.salvage} % of CAPEX`], ['Flexible NPV', mUSD(npv0 + expand + abandonOpt, 2), 'M$', 'static NPV + expansion + abandonment options'],
  ] });
  tables.push({ title: 'Optimisation results', columns: ['Problem', 'Method', 'Solution', 'Objective', 'Note'], rows: [
    ...solvers.map((s) => ['Operating point (NPV-optimal rate and preservation)', s.name, `rate ${pct(s.x[0], 1)} %, preservation ${pct(s.x[1], 0)} %`, `NPV ${mUSD(s.npv, 2)} M$`, `${s.evals ?? '—'} evaluations${s.pen > 1e-6 ? ', constraint violated' : ''}`]),
    ['Operating envelope', 'constraint scan', envelope ? `${pct(envelope.lo, 0)}–${pct(envelope.hi, 0)} % of the case rate` : 'no feasible profitable rate', `optimum ${pct(opBest.x[0], 0)} %`, opCons.map((c) => `${c.name} ${c.g <= 0 ? 'ok' : 'violated'}`).join('; ')],
    ...minlp.slice(0, 3).map((mR, i) => ['Strategy × rate (mixed-integer nonlinear)', 'enumeration + Nelder–Mead', `${mR.name} at ${pct(mR.rate, 0)} %`, `NPV ${mUSD(mR.npv, 2)} M$`, i === 0 ? 'best combination' : `rank ${i + 1}`]),
    ['Project portfolio (linear relaxation)', 'two-phase simplex', lp.status === 'optimal' ? pf.map((r, i) => (lp.x[i] > 1e-6 ? `${r.name} ${pct(lp.x[i], 0)} %` : null)).filter(Boolean).join(', ') || 'none' : lp.status, lp.status === 'optimal' ? `NPV ${rd(lp.value, 2)} M$` : '—', `budget ${q.budget} M$, ${q.vesselDays} vessel days`],
    ['Project portfolio (go / no-go)', 'branch and bound', milp.status === 'optimal' ? pf.filter((_, i) => milp.x[i] > 0.5).map((r) => r.name).join(', ') || 'none' : milp.status, milp.status === 'optimal' ? `NPV ${rd(milp.value, 2)} M$` : '—', `${milp.nodes} nodes`],
    ['Water-handling capacity (two-stage stochastic)', 'deterministic-equivalent LP', ts2.status === 'optimal' ? `${rd(ts2.x, 0)} m³/d installed` : ts2.status, ts2.status === 'optimal' ? `expected cost ${mUSD(ts2.cost, 2)} M$` : '—', `${q.waterCapCost} k$ per m³/d now against ${q.waterPenalty} $/m³ later`],
    ['Asset replacement timing', 'dynamic programming', dp.replaceYears.length ? `replace in year ${dp.replaceYears.join(', ')}` : 'never replace', `PV cost ${mUSD(dp.cost, 2)} M$`, `never replacing costs ${mUSD(dpNever, 2)} M$`],
    ['Inspection interval', 'grid search on the virtual-age model', `${rd(rbiBest.t, 1)} y`, `PV ${mUSD(rbiBest.total, 2)} M$`, `current ${q.inspInterval} y: ${mUSD(rbiNow.total, 2)} M$`],
    ['Pigging interval', 'grid search', `${rd(pigBest.tau, 0)} d`, `${mUSD(pigBest.total, 3)} M$/y`, `stuck-pig probability ${pct(pigBest.pStuck, 2)} % per run`],
    ['Critical spare', 'expected-cost comparison', spare.saving > 0 ? 'hold the spare' : 'no spare', `net ${mUSD(spare.saving, 3)} M$/y`, `break-even failure rate ${rd(spare.beRate, 3)} /y`],
    ['Insulation under uncertainty', 'sample-average approximation', `expected NPV: ${rd(pMean.t, 0)} mm; CVaR: ${rd(pCvar.t, 0)} mm`, `E[NPV] ${mUSD(pMean.mean)} M$`, `${q.nMCopt} common random samples per design`],
    ['Value – risk – carbon', 'non-dominated sorting GA', `${front.length} Pareto designs`, front.length ? `NPV ${rd(front[0].npv, 0)}–${rd(front[front.length - 1].npv, 0)} M$` : '—', `population ${q.nPop}, ${q.nGens} generations`],
  ] });
  tables.push({ title: 'Economically optimised strategies subject to the engineering safety constraints', columns: ['Area', 'Decision', 'Optimum', 'Economic result', 'Constraints applied', 'Status'], rows: [
    ['Mitigation', 'Flow-assurance strategy', bestOpt.name, `lifecycle cost ${mUSD(bestOpt.lcc)} M$, NPV ${mUSD(bestOpt.npv)} M$`, `cooldown ≥ ${q.minCooldown} h; blockage ≤ ${q.maxBlockFreq} /y; injection ≤ ${q.maxInject} m³/d`, noneSafe ? 'no option passes; shown for reference' : `${safeOpts.length} of ${options.length} options pass`],
    ['Design', 'Insulation thickness', `${rd(insOpt.t, 0)} mm`, `lifecycle cost ${mUSD(insOpt.lcc)} M$`, 'same limits as the strategy', insSafe.length ? (insFree === insOpt ? 'constraints not binding' : `binding (cost minimum alone: ${rd(insFree.t, 0)} mm)`) : 'no thickness passes'],
    ['Design', 'Inner diameter', dBest ? `${rd(dBest.d * 1000, 0)} mm` : '—', dBest ? `NPV ${mUSD(dBest.npv)} M$` : '—', `erosional ratio ≤ 1; inlet pressure ≤ ${q.mawp} bara; deliverability`, dBest ? `${dOk.length} of ${diam.length} diameters pass` : 'no diameter passes'],
    ['Design', 'Corrosion control / material', matBest.name, `lifecycle cost ${mUSD(matBest.lcc)} M$`, 'corrosion allowance ≤ 10 mm', `${matPool.length} of ${mats.length} options pass`],
    ['Operation', 'Production rate', `${pct(opBest.x[0], 0)} % of the case rate`, `NPV ${mUSD(opBest.npv)} M$`, opCons.map((c) => c.name).join('; '), opBest.pen < 1e-6 ? (active.length ? 'binding: ' + active.join(', ') : 'inside all limits') : 'limits cannot all be met'],
    ['Operation', 'Long shutdowns preserved', `${pct(opBest.x[1], 0)} %`, `${rd(opProj(opBest.x[0], opBest.x[1]).sm.events, 3)} blockages a year`, `blockage ≤ ${q.maxBlockFreq} /y`, ''],
    ['Inspection', 'In-line inspection interval', `${rd(rbiBest.t, 1)} y`, `PV ${mUSD(rbiBest.total, 2)} M$`, `annual failure probability ≤ ${q.maxPof}`, rbiOk.length ? (rbiFree === rbiBest ? 'constraint not binding' : `binding (cost minimum alone: ${rd(rbiFree.t, 1)} y)`) : 'limit cannot be met by inspection'],
    ['Intervention', 'Pigging interval', `${rd(pigBest.tau, 0)} d`, `${mUSD(pigBest.total, 3)} M$/y`, `deposit ≤ ${q.waxCrit} mm`, pigOk.length ? 'passes' : 'limit cannot be met'],
    ['Intervention', 'Critical spare', spare.saving > 0 ? 'hold' : 'do not hold', `${mUSD(spare.saving, 3)} M$/y net`, 'none', ''],
    ['Intervention', 'Replacement of the ageing asset', dp.replaceYears.length ? `year ${dp.replaceYears.join(', ')}` : 'not within the horizon', `PV ${mUSD(dp.cost, 2)} M$`, 'none', ''],
  ], note: flags.length ? `The case itself breaks ${flags.length} safety constraint(s): ${flags.join('; ')}.` : 'The case design meets every safety constraint that is checked here.' });
  tables.push({ title: 'Insulation and diameter studies', columns: ['Study', 'Design', 'U (W/m²K)', 'Arrival T (°C)', 'Hydrate subcooling (°C)', 'Rate (% plan)', 'Erosional ratio', 'Slug catcher (m³)', 'CAPEX (M$)', 'Result (M$)', 'Status'], rows: [
    ...S.sweep.map((s) => { const c = insCurve[argBest(insCurve.map((x) => -Math.abs(x.t - s.tMm)))]; return ['Insulation', `${rd(s.tMm, 0)} mm`, rd(s.U, 2), rd(s.tArr, 1), rd(s.sub, 1), 100, rd(S.base.eros, 2), rd(q.slugVol, 0), mUSD(c.capex), mUSD(c.lcc), `lifecycle cost; cooldown ${c.cd >= 1e4 ? 'never' : rd(c.cd, 0) + ' h'}`]; }),
    ...diam.map((d) => (d.npv === null ? ['Diameter', `${rd(d.d * 1000, 0)} mm`, '—', '—', '—', '—', '—', '—', '—', '—', d.why] : ['Diameter', `${rd(d.d * 1000, 0)} mm`, rd(d.U, 2), rd(d.tArr, 1), rd(d.sub, 1), pct(d.rate, 0), rd(d.eros * d.rate, 2), rd(d.slugVol, 0), mUSD(d.capexTotal), mUSD(d.npv), `NPV; ${d.why}`])),
  ], note: `${S.calls} steady solutions of the case line on the flow kernel (${q.nCells} cells). Rate response fitted from the kernel: friction ∝ rate^${rd(S.nf, 2)}, thermal exponent ${rd(S.mExp, 2)}.` });
  tables.push({ title: 'Integrity and reliability economics', columns: ['Item', 'Value', 'Unit', 'Note'], rows: [
    ['Annual failure probability (year 1)', rd(p0.haz[0], 5), '1/y', 'base rate + wear-out hazard'], ['Annual failure probability (last year)', rd(p0.haz[lifeEff - 1], 5), '1/y', `Weibull shape ${q.weibullBeta}, ${q.pEnd} % cumulative at ${q.remLife} y`], ['Consequence per failure', q.consequence, 'M$', ''], ['Expected annual loss (year 1)', mUSD(eal, 3), 'M$/y', 'probability × consequence'], ['Lifecycle expected failure cost', mUSD(lifeFail, 2), 'M$', 'present value, = NPV − risk-adjusted NPV'],
    ['Downtime cost of one failure', mUSD(q.repairDays * dayValue, 2), 'M$', `${q.repairDays} d at ${mUSD(dayValue, 2)} M$/d after deferral credit`], ['Optimum inspection interval', rd(rbiBest.t, 1), 'y', `PV ${mUSD(rbiBest.total, 2)} M$ (inspection ${mUSD(rbiBest.inspection, 2)}, failure ${mUSD(rbiBest.failure, 2)})`], ['Current inspection interval', q.inspInterval, 'y', `PV ${mUSD(rbiNow.total, 2)} M$`],
    ...mats.map((mR) => [`Material: ${mR.name}`, mUSD(mR.lcc, 2), 'M$ lifecycle', `rate ${rd(mR.rate, 3)} mm/y, allowance ${rd(mR.ca, 1)} mm${mR.practical ? '' : ' (capped; line replaced at year ' + rd(mR.life, 0) + ')'}, ΔCAPEX ${mUSD(mR.dCap, 2)}, PV OPEX ${mUSD(mR.opex, 2)}, PV failure ${mUSD(mR.fail, 2)}${mR.repl > 0 ? ', PV replacement ' + mUSD(mR.repl, 2) : ''}`]),
  ] });
  tables.push({ title: 'Emissions and abatement', columns: ['Item', 'tCO₂e/y', 'Abatement (t/y)', 'Annual cost difference (M$/y)', 'Abatement cost ($/t)'], rows: [...em0.items.map((i) => [i.source, rd(i.t, 0), '—', '—', '—']), ['Total (case strategy)', rd(carbonT, 0), '—', rd((carbonT * q.carbonPrice) / MM, 3) + ' carbon cost', '—'], ['Embodied in line-pipe steel (one-off)', rd(capex.steelT * EF.steel, 0), '—', '—', '—'], ...mac.map((x) => [`Strategy: ${x.name}`, rd(x.emis, 0), rd(x.abate, 0), rd(x.cost / MM, 3), rd(x.mac, 0)])], note: `Power factor ${rd(X.enCarbon, 2)} kgCO₂/kWh (${q.powerSource === 'gas' ? 'own gas turbines' : 'grid'}); abatement is measured against "${refEm.name}", the highest-emitting strategy.` });
  if (front.length) tables.push({ title: 'Pareto designs (sample)', columns: ['Insulation (mm)', 'Rate (% of case)', 'Shutdowns preserved (%)', 'NPV (M$)', 'Expected annual loss (M$/y)', 'Emissions (tCO₂e/y)'], rows: front.filter((_, i) => i % Math.max(1, Math.floor(front.length / 12)) === 0 || i === front.length - 1).map((f) => [rd(f.t, 0), pct(f.rate, 0), pct(f.treat, 0), rd(f.npv, 1), rd(f.risk, 3), rd(f.carbon, 0)]) });

  // ---- balances (conservation and consistency identities)
  const phaseSum = sum(q.phase), treeP = sum(tree.branches[1].node.branches.map((b) => b.p));
  balances.push(
    { name: 'Σ discounted cash flows = NPV (M$)', in: npv(X.r, cf.fcf, dOpt) / MM, out: npv0 / MM },
    { name: 'Real cash flows at the real rate = nominal NPV (M$)', in: met.realNpv / MM, out: npv0 / MM },
    { name: 'CAPEX phasing sums to the total (M$, real)', in: sum(q.phase.map((f) => f * capex.total)) / MM, out: (capex.total * phaseSum) / MM },
    { name: 'CAPEX items sum to the total (M$)', in: sum(capex.items.map((i) => i.cost)) / MM, out: capex.total / MM },
    { name: 'Depreciation charges sum to the depreciable base (M$)', in: sum(cf.depreciation) / MM, out: cf.deprBase / MM },
    { name: 'Working-capital changes net to zero (M$)', in: 1 + sum(cf.dwc) / MM, out: 1 },
    { name: 'Produced + lost = potential production (million boe)', in: (sum(cf.boe) + sum(cf.lostBoe)) / 1e6, out: sum(cf.potBoe) / 1e6 },
    { name: 'OPEX categories sum to first-year OPEX (M$)', in: sum(opexCat.map((c) => c[1])) / MM, out: opexY1 / MM },
    { name: 'Scenario weights sum to one', in: sum(scen.map((s) => s.weight / wSum)), out: 1 },
    { name: 'Appraisal-signal probabilities sum to one', in: treeP, out: 1 },
    { name: 'Decision tree with appraisal = EMV + value of information − cost (M$)', in: tree.branches[1].emv / MM, out: (voi.emv + voi.evii) / MM - q.testCost },
    { name: 'Histogram counts = sample size', in: sum(hist.counts), out: mc.n },
    { name: 'Currency round trip (M$)', in: ((npv0 / MM) * fx) / fx, out: npv0 / MM },
    { name: 'Break-even price returns zero NPV (M$, offset by 1)', in: 1 + bePriceCheck / MM, out: 1 },
  );

  const summary = `${sm0.name} on the ${rd(q.lineLen / 1000, 1)} km line: CAPEX ${mUSD(capex.total, 0)} M$, NPV ${mUSD(npv0, 0)} M$ (risk-adjusted ${mUSD(riskedNpv, 0)} M$), IRR ${met.irr === null ? 'undefined' : pct(met.irr) + ' %'}, break-even ${be.price === null ? 'not reached' : rd(be.price, 1) + ' $/bbl'}; P10–P90 ${mUSD(mc.p10, 0)} to ${mUSD(mc.p90, 0)} M$ with ${pct(mc.probLoss)} % chance of loss; ${noneSafe ? 'no flow-assurance strategy passes the safety screen' : 'lowest lifecycle-cost strategy that meets the safety constraints: ' + bestOpt.name}${flags.length ? `; ${flags.length} safety constraint(s) violated by the case` : ''}.`;
  prog(1, 'Done');
  return {
    summary, kpis, warnings, recommendations: recs, plots, tables, balances,
    outputs: {
      capex: capex.total, opexPerYear: opexY1, npv: npv0, riskedNpv, irr: met.irr === null ? null : 100 * met.irr, mirr: met.mirr === null ? null : 100 * met.mirr, payback: met.payback, discountedPayback: met.discountedPayback, pi: met.pi, utc: met.utc, breakevenPrice: be.price,
      eal, deferredCost, chemicalCost: sm0.chem + sm0.shutdown + q.chemOther * boeY1, energyCost: p0.fixedItems[6][1], carbon: carbonT, carbonIntensity: intensity, npvP10: mc.p10, npvP50: mc.p50, npvP90: mc.p90, probLoss: mc.probLoss, bestOption: noneSafe ? null : bestOpt.name, ranking,
      npvMean: mc.mean, npvSd: mc.sd, var: mc.var, cvar: mc.cvar, certaintyEquivalent: ce, realNpv: met.realNpv, eav: met.eav, eac: met.eac, roi: met.roi, roce: met.roce, liftingCost: met.liftingCost, governmentTake: met.governmentTake, economicLimit: limit0, lifeEvaluated: lifeEff,
      breakevenProduction: be.prod === null ? null : be.prod * q.qOil, breakevenInhibitorRate: be.inhib, breakevenBlockageFrequency: be.block, lifecycleFailureCost: lifeFail, lifecycleCost: selOpt.lcc, totex: met.pvCapex + met.pvOpex + met.pvAb,
      maintenanceCost: p0.fixedItems[1][1], inspectionCost: p0.fixedItems[3][1], piggingCost: sm0.pig, remediationCost: sm0.block, interventionCost: sm0.block + sm0.shutdown + sm0.pig, downtimeLoss: deferredCost, revenueFirstYear: cf.revenue[y1], totalRevenue: sum(cf.revenue),
      scenarios: scen.map((s) => ({ name: s.name, weight: s.weight / wSum, npv: s.npv })), scenarioNpv: scenEv, safetyFlags: flags, safe: flags.length === 0, bestOptionFeasible: !noneSafe, paretoCount: front.length,
      optimumInsulation: insOpt.t, optimumDiameter: dBest ? dBest.d : null, optimumRate: opBest.x[0], operatingEnvelope: envelope, inspectionIntervalOptimum: rbiBest.t, piggingIntervalOptimum: pigBest.tau, materialChoice: matBest.name,
      evpi: voi.evpi, evii: voi.evii, emv: voi.emv, optionValue: defer, optionExpand: expand, optionAbandon: abandonOpt, blackScholes: bs, mcdaBest: mcdaBest.name, topsisBest: topsisBest.name, ahpWeights: ahpR.weights, consistencyRatio: ahpR.cr,
      strategy: sm0.name, currency: cur, fxPerUSD: fx, npvLocal: npv0 * fx, kernelCalls: S.calls,
    },
  };
}

// ================================================================================================================
// 10. Declarations: inputs, presets, linked data, convergence studies, calibration, verification
// ================================================================================================================
const L0 = (() => { const l = caseLine({}), x = l.profile.x, z = l.profile.z; let fl = 0; for (let i = 1; i < x.length; i++) if (x[i] <= l.riserBaseX + 1e-9) fl += Math.hypot(x[i] - x[i - 1], z[i] - z[i - 1]); return { total: Math.round(l.length), riser: Math.round(l.length - fl) }; })();
const N = (key, label, unit, value, min, max, help, extra = {}) => ({ key, label, unit, value, min, max, help, ...extra });
const sel = (key, label, value, options, help) => ({ key, label, type: 'select', value, options: options.map(([v, l]) => ({ value: v, label: l })), help });
const DIST_DEFAULT = [
  { id: 'price', dist: 'lognormal', lo: 0.65, mode: 1, hi: 1.45 }, { id: 'prod', dist: 'pert', lo: 0.6, mode: 1, hi: 1.3 }, { id: 'capex', dist: 'triangular', lo: 0.9, mode: 1, hi: 1.5 }, { id: 'opex', dist: 'triangular', lo: 0.85, mode: 1, hi: 1.35 },
  { id: 'downtime', dist: 'triangular', lo: 0.5, mode: 1, hi: 2.5 }, { id: 'failFreq', dist: 'lognormal', lo: 0.4, mode: 1, hi: 2.5 }, { id: 'repair', dist: 'uniform', lo: 0.7, mode: 1, hi: 1.6 },
];
const INPUTS = [
  { group: 'Project frame and fiscal terms', tab: 'inputs', help: 'Cash flows are built in money of the day and discounted at the nominal rate; the same result in real terms is reported as a check.', fields: [
    N('evalYear', 'Evaluation year', '', 2026, 1990, 2100, 'Project evaluation date; year 0 of every cash flow.', { int: true }),
    N('life', 'Production life', 'y', BASE.projectLife, 1, 50, 'Years of production after first oil.', { int: true, typical: [10, 30] }),
    { key: 'phasing', label: 'CAPEX phasing', type: 'table', help: 'Percentage of CAPEX spent in each construction year; first production follows the last row. Rows are normalised to 100 %.', columns: [{ key: 'year', label: 'Construction year' }, { key: 'pct', label: 'Share', unit: '%' }], value: [{ year: 1, pct: 40 }, { year: 2, pct: 60 }] },
    N('discount', 'Discount rate (nominal)', '%/y', BASE.discountRate, 0, 40, 'Required rate of return in money-of-the-day terms.', { typical: [8, 15] }),
    { key: 'mid', label: 'Mid-year discounting', type: 'bool', value: false, help: 'Flows of year k ≥ 1 fall at k − ½ instead of the year end.' },
    N('reinvest', 'Reinvestment rate (MIRR)', '%/y', 10, 0, 40, 'Rate earned on positive cash flows for the modified IRR.'),
    N('inflation', 'Inflation', '%/y', 2.5, -2, 30, 'General inflation; converts between real and nominal money.'),
    N('costEsc', 'Cost escalation above inflation', '%/y', 0, -5, 15, 'Real escalation of CAPEX and OPEX.'),
    N('priceEsc', 'Price escalation above inflation', '%/y', 0, -10, 15, 'Real escalation of oil and gas prices (0 = flat real prices).'),
    N('hurdle', 'Hurdle rate (minimum IRR)', '%/y', 12, 0, 60, 'Minimum acceptable internal rate of return.'),
    N('maxPayback', 'Maximum acceptable payback', 'y', 6, 0.5, 40, 'Decision threshold for the simple payback.'),
    N('minNpv', 'Minimum NPV', 'M$', 0, -1e4, 1e5, 'Decision threshold for the net present value.'),
    sel('regime', 'Fiscal regime', 'tax', [['tax', 'Royalty and income tax'], ['psc', 'Production-sharing contract']], 'Royalty/tax concession or a simple production-sharing contract.'),
    N('royalty', 'Royalty', '% of revenue', 10, 0, 60, 'Taken off gross revenue before anything else.'),
    N('taxRate', 'Income tax rate', '%', 30, 0, 90, 'Applied to taxable income after depreciation and losses brought forward (contractor profit share under a PSC).'),
    N('costOilCap', 'Cost-oil cap', '% of net revenue', 60, 5, 100, 'Largest share of net revenue available for cost recovery each year.', { showIf: (v) => v.regime === 'psc' }),
    N('profitSplit', 'Contractor profit share', '%', 40, 1, 100, 'Contractor share of profit oil.', { showIf: (v) => v.regime === 'psc' }),
    sel('deprMethod', 'Depreciation', 'sl', [['sl', 'Straight line'], ['db', 'Declining balance'], ['uop', 'Units of production']], 'Tax depreciation of capitalised cost from first production.'),
    N('deprLife', 'Depreciation life', 'y', 10, 1, 40, 'Years over which CAPEX is written off.', { int: true }),
    N('dbRate', 'Declining-balance rate', '%/y', 20, 1, 100, 'Annual rate; switches to straight line when that is larger.', { showIf: (v) => v.deprMethod === 'db' }),
    N('wcDays', 'Working capital', 'days of revenue', 30, 0, 180, 'Receivables less payables, released at the end.'),
    N('abandon', 'Abandonment cost', 'M$', 60, 0, 5000, 'Decommissioning cost in evaluation-year money, paid the year after production stops.'),
    { key: 'abandonProvision', label: 'Provide for abandonment (tax-deductible accrual)', type: 'bool', value: true, help: 'When off, the cost is only deductible when spent and usually finds no income to shelter.' },
    { key: 'stopAtLimit', label: 'Stop at the economic limit', type: 'bool', value: true, help: 'Cease production in the first year whose operating margin is negative.' },
    N('capexSunk', 'CAPEX already committed', 'M$', 0, 0, 1e5, 'Sunk cost: excluded from the forward NPV but still depreciable.'),
    N('residual', 'Book value of existing equipment', 'M$', 0, 0, 1e5, 'Remaining equipment value carried into the depreciable base.'),
    N('gearing', 'Debt share of CAPEX', '%', 0, 0, 90, 'Financing assumption for the equity view (0 = all equity).'),
    N('loanRate', 'Loan interest rate', '%/y', 7, 0, 30, 'Interest during construction is capitalised.', { showIf: (v) => v.gearing > 0 }),
    N('loanTenor', 'Loan tenor', 'y', 8, 1, 30, 'Level annuity repayment from first production.', { int: true, showIf: (v) => v.gearing > 0 }),
    { key: 'currency', label: 'Reporting currency', type: 'text', value: 'USD', help: 'Headline results are also shown in this currency.' },
    N('fx', 'Exchange rate', 'per US$', 1, 1e-6, 1e7, 'Units of the reporting currency per US dollar.'),
  ] },
  { group: 'Production and prices', tab: 'inputs', fields: [
    N('qOil', 'Initial oil rate', 'Sm³/d', 3000, 0, 2e5, 'Plateau stock-tank oil or condensate rate.'),
    N('qGas', 'Initial gas rate', 'Sm³/d', 542000, 0, 2e8, 'Plateau produced-gas rate.'),
    N('gasSalesFrac', 'Gas sold', '% of produced', 85, 0, 100, 'The rest is fuel, flare and shrinkage.'),
    N('gasHV', 'Gas heating value', 'MJ/Sm³', 39, 20, 60, 'Gross heating value used to convert volume to energy.'),
    N('plateau', 'Plateau length', 'y', 3, 0, 30, 'Years at the initial rate before decline.'),
    sel('declineType', 'Decline model', 'exp', [['exp', 'Arps exponential'], ['hyp', 'Arps hyperbolic'], ['har', 'Arps harmonic']], 'Arps decline after the plateau.'),
    N('Di', 'Initial decline', '%/y', 15, 0, 90, 'Nominal decline rate at the end of the plateau.'),
    N('bHyp', 'Hyperbolic exponent b', '–', 0.5, 0.05, 0.95, 'Arps b for the hyperbolic model.', { showIf: (v) => v.declineType === 'hyp' }),
    N('wc0', 'Initial water cut', '%', 20, 0, 95, 'Water cut of the standard liquid at first oil.'),
    N('wcEnd', 'Final water cut', '%', 70, 0, 98, 'Water cut at the end of the project life.'),
    N('uptime', 'Uptime', '%', 95, 5, 100, 'Production availability over the year (planned and unplanned stops).', { typical: [88, 98] }),
    N('deferFrac', 'Downtime volume recovered later', '%', 60, 0, 100, 'Deferred (recovered in the last year, discounted) rather than lost.'),
    N('oilPrice', 'Oil price', '$/bbl', BASE.oilPrice, 0, 500, 'Evaluation-year price, real terms.'),
    N('gasPrice', 'Gas price', '$/MMBtu', BASE.gasPrice, 0, 100, 'Evaluation-year price, real terms.'),
    N('tariff', 'Transport and processing tariff', '$/boe', 2, 0, 50, 'Existing tariffs and transportation charges.'),
    N('pi', 'Productivity index', 'Sm³/d/bar', BASE.pi, 0.1, 1e4, 'Liquid inflow per bar of drawdown; turns extra back-pressure into lost rate.'),
    N('chokeDp', 'Choke pressure margin', 'bar', 10, 0, 300, 'Pressure drop across the production choke at the case rate; available to absorb extra line losses.'),
    N('capacityFrac', 'Facility capacity', '× case rate', 1.15, 0.3, 3, 'Largest rate the facilities and the operating envelope accept.'),
    N('qMinFrac', 'Minimum stable rate', '× case rate', 0.45, 0.05, 1, 'Turndown limit below which slugging or solids make operation unstable.'),
  ] },
  { group: 'Line and flow-assurance strategy', tab: 'inputs', help: 'The arrival temperature, hydrate margin, pressure drop and slug volume come from the flow kernel on the case line; lengths and sizes here drive the cost.', fields: [
    sel('strategy', 'Flow-assurance strategy of the case', 'wet', Object.entries(STRATEGIES).map(([k, s]) => [k, s.name]), 'The strategy that is costed in the base cash flow; all six are compared.'),
    N('insMm', 'Wet insulation thickness', 'mm', BASE.insulation.t * 1000, 0, 250, 'Thickness on the flowline for the insulated strategies.'),
    N('kIns', 'Insulation conductivity', 'W/m/K', BASE.insulation.k, 0.02, 1, 'Thermal conductivity of the wet insulation.'),
    N('lineLen', 'Line length (flowline + riser)', 'm', L0.total, 100, 5e5, 'Total length that is bought and installed.'),
    N('riserLen', 'Riser length', 'm', L0.riser, 0, 5000, 'Part of the line that is riser.'),
    N('idMm', 'Inner diameter', 'mm', BASE.idMm, 50, 1500, 'Flowline bore.'),
    N('wtMm', 'Wall thickness', 'mm', BASE.wtMm, 3, 80, 'Steel wall.'),
    N('depth', 'Water depth', 'm', BASE.waterDepth, 0, 4000, 'Raises the installation day rate.'),
    N('nWells', 'Wells', '', 2, 1, 40, 'Subsea wells tied back.', { int: true }),
    N('slugVol', 'Slug-catcher volume', 'm³', BASE.slugCatcherVol, 1, 5000, 'Surge volume the receiving vessel must hold.'),
    N('pumpKW', 'Pumping power', 'kW', 0, 0, 1e5, 'Continuous pump power.'),
    N('compKW', 'Compression power', 'kW', 2000, 0, 2e5, 'Continuous compressor power.'),
    N('heatKW', 'Continuous heating power', 'kW', 0, 0, 1e5, 'Heating that runs all year (0 when none).'),
    N('shutdowns', 'Unplanned shutdowns', '1/y', 6, 0, 100, 'Shutdowns a year that can lead to cooldown.'),
    N('shutdownMean', 'Mean shutdown duration', 'h', 12, 0.5, 500, 'Durations are taken as exponentially distributed.'),
    N('reactH', 'Time needed to preserve the line', 'h', 4, 0, 100, 'Part of the cooldown time used up before the line is safe.'),
    N('restartH', 'Restart delay after preservation', 'h', 8, 0, 200, 'Extra production loss when a line has been displaced or dosed.'),
    N('plugProb', 'Hydrate plug probability per unprotected event', '–', 0.05, 0, 1, 'Chance that an unprotected excursion into the hydrate region ends in a blockage.'),
    N('hydMargin', 'Hydrate design margin', '°C', 3, 0, 15, 'Safety margin on the hydrate temperature.'),
    N('wat', 'Wax appearance temperature', '°C', 30, -20, 90, 'Arrival below this temperature raises the pigging frequency.'),
    N('waxRate', 'Wax build-up rate', 'mm/d', 0.01, 0, 5, 'Deposit growth between pig runs.'),
    N('waxCrit', 'Largest deposit thickness allowed', 'mm', 6, 0.2, 50, 'Scale of the stuck-pig probability and upper limit for the pigging interval.'),
    N('pigRuns', 'Routine pig runs', '1/y', 4, 0, 365, 'Baseline pigging frequency.'),
    N('blowdowns', 'Blowdowns', '1/y', 1, 0, 50, 'Depressurisations a year that flare the gas inventory.'),
    sel('inhibitor', 'Inhibitor reported by operations', 'MeOH', [['MeOH', 'Methanol'], ['MEG', 'MEG (regenerated)'], ['LDHI', 'Low-dosage inhibitor']], 'Prices the continuous rate below.'),
    N('inhibRate', 'Continuous inhibitor rate', 'm³/d', 0, 0, 5000, 'When above zero this replaces the modelled inhibitor demand of the case strategy.'),
    N('cooldownBase', 'Cooldown time from operations', 'h', 0, 0, 2000, 'When above zero the lumped cooldown model is scaled to match it.'),
  ] },
  { group: 'Engineering safety constraints and availability', tab: 'inputs', help: 'Limits taken from the flow, solids, operations and integrity studies. An option that breaks one is reported but never recommended; optimisations are restricted to the designs that pass.', fields: [
    N('minCooldown', 'Required cooldown time', 'h', 8, 0, 200, 'Shortest time to hydrate conditions accepted for strategies that rely on insulation (no-touch time plus preservation).'),
    N('maxBlockFreq', 'Tolerable hydrate blockage frequency', '1/y', 0.05, 0.0001, 10, 'Options with a higher expected blockage frequency are excluded.'),
    N('maxInject', 'Chemical injection capacity', 'm³/d', 150, 0.1, 5000, 'Largest continuous once-through inhibitor rate that can be supplied, stored and injected.'),
    N('inhibAvail', 'Availability of the inhibition or heating system', '%', 98, 50, 100, 'Chemical supply and equipment availability; the unavailable share leaves the line unprotected.'),
    N('vesselWait', 'Wait for an intervention vessel', 'd', 15, 0, 365, 'Vessel availability: added to the outage of every blockage.'),
    N('mawp', 'Allowable inlet pressure', 'bara', BASE.designPressure, 5, 2000, 'Design or maximum allowable working pressure of the line.'),
    N('maxPof', 'Tolerable annual probability of failure', '1/y', 0.02, 0.00001, 1, 'Target that the inspection interval must respect.'),
    N('integUtil', 'Largest structural utilisation', '–', 0.7, 0, 10, 'Hoop, combined-stress or collapse utilisation from the integrity study; above 1 the design fails.'),
    N('integViol', 'Integrity code checks failed', '', 0, 0, 1000, 'Count of violated checks from the integrity study.', { int: true }),
    N('erosIn', 'Erosional velocity ratio from the flow study', '–', 0, 0, 20, '0 uses the kernel estimate on the case line.'),
    { key: 'severeSlug', label: 'Severe slugging predicted and not suppressed', type: 'bool', value: false },
  ] },
  { group: 'Integrity and reliability', tab: 'inputs', fields: [
    N('pof', 'Annual probability of failure', '1/y', 0.002, 0, 0.5, 'Random (time-independent) failure rate.'),
    N('consequence', 'Consequence of one failure', 'M$', 150, 0, 1e4, 'Repair, clean-up and production loss.'),
    N('ealOverride', 'Expected annual loss from the integrity study', 'M$/y', 0, 0, 1e4, 'When above zero it sets the base failure rate (= this ÷ consequence).'),
    N('remLife', 'Remaining asset life', 'y', 25, 0.5, 200, 'Time until wear-out (wall loss, fatigue) is expected.'),
    N('pEnd', 'Wear-out failure probability at the remaining life', '%', 50, 1, 99, 'Fixes the scale of the Weibull hazard.'),
    N('weibullBeta', 'Weibull shape', '–', 3, 1, 8, 'Above 1 the hazard rises with age.'),
    N('assetAge', 'Present age of the asset', 'y', 0, 0, 60, 'Starting state of the replacement study.', { int: true }),
    N('inspInterval', 'Inspection interval', 'y', 5, 0.25, 30, 'Current in-line inspection interval.'),
    N('inspCost', 'Cost of one inspection', 'M$', 1.5, 0, 100, 'In-line inspection campaign.'),
    N('pod', 'Probability of detection', '%', 80, 0, 100, 'Share of the accumulated degradation found and repaired at an inspection.'),
    N('corrRate', 'Corrosion rate', 'mm/y', 0.12, 0, 10, 'Wall-loss rate as currently managed.'),
    { key: 'corrInhibited', label: 'That corrosion rate is with inhibitor', type: 'bool', value: true },
    N('inhEff', 'Corrosion-inhibitor effectiveness × availability', '%', 90, 0, 99, 'Reduction of the uninhibited rate.'),
    N('corrAllow', 'Corrosion allowance', 'mm', BASE.corrosionAllowanceMm, 0, 12, 'Allowance in the case wall thickness.'),
    N('corrInhCost', 'Corrosion-inhibitor cost', 'M$/y', 0.6, 0, 50, 'Chemical and injection cost of the inhibited option.'),
    N('repairDays', 'Downtime for a repair', 'd', 45, 0, 720, 'Production outage of one failure.'),
    N('replCost', 'Replacement cost of the ageing asset', 'M$', 30, 0, 5000, 'Replacement study: cost of a new section or unit.'),
    N('maintAsset', 'Its maintenance cost when new', 'M$/y', 0.5, 0, 100, 'Rises 4 % for each year of age.'),
    N('spareCost', 'Cost of the critical spare', 'M$', 3, 0, 500, 'Spare held against a long lead time.'),
    N('itemFailRate', 'Failure rate of that item', '1/y', 0.08, 0, 5, 'Random failures a year.'),
    N('leadNo', 'Outage without a spare', 'd', 120, 0, 1000, 'Lead time to buy and install.'),
    N('leadWith', 'Outage with a spare', 'd', 14, 0, 365, 'Time to mobilise and install.'),
    N('spareShare', 'Production affected', '%', 50, 0, 100, 'Share of production lost while the item is out.'),
  ] },
  { group: 'CAPEX cost basis', tab: 'setup', help: 'Order-of-magnitude industry values in 2023 US dollars (cost index 800). Purchased cost = ref × (capacity ÷ reference capacity)^exponent; installed cost = purchased × factor.', fields: [
    { key: 'costBasis', label: 'Equipment cost basis', type: 'table', help: 'Edit ref (M$), reference capacity, exponent (0.6 is the six-tenths rule) and installation factor. Keep the id column.', columns: [{ key: 'id', label: 'Id', type: 'text' }, { key: 'item', label: 'Item', type: 'text' }, { key: 'ref', label: 'Ref. cost', unit: 'M$' }, { key: 'cap', label: 'Ref. capacity' }, { key: 'unit', label: 'Capacity unit', type: 'text' }, { key: 'exp', label: 'Exponent' }, { key: 'fac', label: 'Installation factor' }], value: COST_BASIS.map((r) => ({ ...r })) },
    sel('costMethod', 'Topsides installed cost', 'module', [['module', 'Bare-module factor for each item'], ['lang', 'One Lang factor on purchased cost']], 'Both are shown in the CAPEX table.'),
    N('langFactor', 'Lang factor', '–', 3.6, 1, 10, 'Installed cost ÷ purchased equipment cost for the whole topsides scope.'),
    sel('pipeMethod', 'Line-pipe cost', 'bottom', [['bottom', 'Bottom-up: steel tonnage, coating, welding'], ['cer', 'Parametric cost-estimating relationship']], 'The parametric relationship can be calibrated to past projects.'),
    N('steelPrice', 'Line-pipe steel', '$/t', 1800, 300, 20000, 'Delivered X65 line pipe.'),
    N('coatPrice', 'Anti-corrosion coating', '$/m²', 60, 0, 1000, 'Three-layer polypropylene or equivalent.'),
    N('fabPerM', 'Welding and field joints', '$/m', 120, 0, 5000, 'Double-jointing, NDT and field-joint coating.'),
    N('insPrice', 'Wet insulation applied', '$/m³', 5000, 200, 50000, 'Syntactic or solid polyurethane / polypropylene.'),
    N('pipPremium', 'Pipe-in-pipe premium', '$/m', 650, 0, 10000, 'Annulus insulation, centralisers, bulkheads and assembly.'),
    N('dehCable', 'Heating cable and anodes', '$/m', 450, 0, 10000, 'Piggy-back cable for direct electrical heating.'),
    N('craFactor', 'CRA-clad pipe cost factor', '× carbon steel', 4.5, 1, 20, 'Line-pipe cost multiplier for corrosion-resistant cladding.'),
    N('riserFactor', 'Riser cost factor', '× flowline per metre', 2.5, 1, 20, 'Fatigue-class welds, strakes, flex joint.'),
    N('vesselRate', 'Installation vessel day rate', 'k$/d', 350, 20, 3000, 'Pipelay or construction vessel spread.'),
    N('layRate', 'Lay rate at 10 in', 'km/d', 2.5, 0.05, 20, 'Scaled with diameter and thermal system.'),
    N('mobCost', 'Mobilisation and demobilisation', 'M$', 6, 0, 500, 'Lump sum.'),
    N('depthCoef', 'Depth factor on the day rate', 'per 1000 m', 0.2, 0, 2, 'Day rate × (1 + this × depth / 1000 m).'),
    N('wellCost', 'Drilling and completion per well', 'M$', 70, 0, 1000, 'First well; later wells follow the learning curve.'),
    N('cerCoef', 'Parametric line cost at 10 in', 'M$/km', 0.4, 0.05, 20, 'Coefficient of the cost-estimating relationship (line pipe, coating, welding).'),
    N('cerExp', 'Parametric diameter exponent', '–', 1.3, 0.2, 3, 'Cost ∝ (diameter ÷ 10 in)^exponent.'),
    N('layFactor', 'Installation day-rate factor', '–', 1, 0.2, 5, 'Multiplies the vessel spread cost (market tightness, weather downtime).'),
    N('learnRate', 'Learning-curve rate', '%', 90, 60, 100, 'Each doubling of repeated units costs this share of the previous.'),
    N('costIndexBase', 'Cost index of the basis year', '', 800, 100, 5000, 'Plant cost index of the reference costs.'),
    N('costIndexEval', 'Cost index of the evaluation year', '', 830, 100, 5000, 'Escalates the basis to the evaluation year.'),
    N('locFactor', 'Location factor', '–', 1, 0.3, 4, 'Regional cost level relative to the basis.'),
    N('contingency', 'Contingency', '% of direct', 15, 0, 100, 'Allowance for undefined scope.'),
    N('owners', "Owner's costs", '% of direct', 8, 0, 60, 'Project team, insurance, studies.'),
  ] },
  { group: 'OPEX and carbon basis', tab: 'setup', fields: [
    N('opsFixed', 'Operations, labour and logistics', 'M$/y', 12, 0, 2000, 'Fixed production-operations cost.'),
    N('maintPct', 'Maintenance', '% of CAPEX per year', 2.5, 0, 20, 'Routine maintenance.'),
    N('insurPct', 'Insurance', '% of CAPEX per year', 0.6, 0, 10, ''),
    N('corrMgmt', 'Corrosion management', 'M$/y', 0.8, 0, 100, 'Monitoring, coupons, cathodic protection.'),
    N('chemOther', 'Production chemicals', '$/boe', 0.6, 0, 20, 'Demulsifier, scale and corrosion inhibitors.'),
    N('waterCost', 'Produced-water handling', '$/m³', 2.5, 0, 100, 'Treatment and disposal.'),
    N('pigCost', 'Cost of a pig run', 'k$', 60, 0, 5000, 'Pigs, labour, deferred production.'),
    N('meohPrice', 'Methanol', '$/m³', 550, 50, 5000, 'Delivered offshore.'),
    N('megPrice', 'MEG', '$/m³', 1100, 100, 8000, 'Delivered offshore.'),
    N('megLoss', 'MEG losses', '% of circulation', 1, 0, 100, 'Make-up needed with regeneration.'),
    N('ldhiPrice', 'Low-dosage inhibitor', '$/m³', 9000, 500, 60000, ''),
    N('ldhiDose', 'Low-dosage inhibitor dose', 'vol % of water', 0.5, 0.05, 5, ''),
    N('spreadRate', 'Intervention vessel spread', 'k$/d', 250, 10, 3000, 'Vessel, coiled tubing or ROV spread for remediation.'),
    N('remedDays', 'Remediation campaign', 'd', 20, 0, 365, 'Vessel days to clear one blockage.'),
    N('blockDays', 'Production outage of a blockage', 'd', 30, 0, 720, ''),
    sel('powerSource', 'Power source', 'grid', [['grid', 'Grid / power from shore'], ['gas', 'Own gas turbines (fuel gas)']], 'Sets the energy price and the emission factor.'),
    N('elecPrice', 'Electricity price', '$/kWh', 0.12, 0, 2, ''),
    N('gridCarbon', 'Grid carbon intensity', 'kgCO₂/kWh', 0.45, 0, 2, ''),
    N('carbonPrice', 'Carbon price', '$/tCO₂e', 50, 0, 1000, 'Charged on the emission inventory each year.'),
    N('carbonEsc', 'Carbon-price escalation above inflation', '%/y', 3, -5, 20, ''),
    N('flareFrac', 'Routine flaring', '% of produced gas', 0.5, 0, 100, ''),
    N('waterCapCost', 'Water-handling capacity', 'k$ per m³/d', 12, 0.1, 500, 'First-stage cost in the capacity study.'),
    N('waterPenalty', 'Cost of water above capacity', '$/m³', 15, 0.1, 500, 'Recourse cost in the capacity study.'),
  ] },
  { group: 'Uncertainty', tab: 'setup', help: 'Each row is a multiplier on the base value. Triangular and PERT use low / most likely / high; normal and lognormal read low and high as P10 and P90.', fields: [
    { key: 'dists', label: 'Distributions', type: 'table', columns: [{ key: 'id', label: 'Input', type: 'text' }, { key: 'dist', label: 'Type', type: 'text' }, { key: 'lo', label: 'Low' }, { key: 'mode', label: 'Most likely' }, { key: 'hi', label: 'High' }], value: DIST_DEFAULT.map((r) => ({ ...r })), help: 'Inputs: price, prod, capex, opex, downtime, failFreq, repair. Types: triangular, pert, uniform, normal, lognormal.' },
    { key: 'corr', label: 'Correlations', type: 'table', columns: [{ key: 'a', label: 'Input A', type: 'text' }, { key: 'b', label: 'Input B', type: 'text' }, { key: 'rho', label: 'Correlation' }], value: [{ a: 'price', b: 'opex', rho: 0.4 }, { a: 'price', b: 'capex', rho: 0.3 }, { a: 'capex', b: 'opex', rho: 0.3 }, { a: 'downtime', b: 'failFreq', rho: 0.5 }], help: 'Pairs not listed are independent.' },
    { key: 'scenarios', label: 'Scenarios', type: 'table', columns: [{ key: 'name', label: 'Scenario', type: 'text' }, { key: 'weight', label: 'Weight' }, { key: 'price', label: 'Price ×' }, { key: 'prod', label: 'Production ×' }, { key: 'capex', label: 'CAPEX ×' }, { key: 'opex', label: 'OPEX ×' }], value: [{ name: 'Low', weight: 0.25, price: 0.65, prod: 0.75, capex: 1.25, opex: 1.15 }, { name: 'Base', weight: 0.5, price: 1, prod: 1, capex: 1, opex: 1 }, { name: 'High', weight: 0.25, price: 1.35, prod: 1.2, capex: 0.95, opex: 0.95 }], help: 'Weights are normalised.' },
    sel('sampling', 'Sampling', 'lhs', [['lhs', 'Latin hypercube (Iman–Conover correlation)'], ['mc', 'Monte Carlo (Gaussian copula)']], ''),
    sel('priceModel', 'Price path in the simulation', 'ou', [['static', 'One multiplier for the whole life'], ['gbm', 'Geometric Brownian motion'], ['ou', 'Mean-reverting']], 'Annual path around the escalated trend.'),
    N('priceVol', 'Price volatility', '%/y', 25, 0, 150, ''),
    N('priceKappa', 'Mean-reversion speed', '1/y', 0.3, 0.01, 5, '', { showIf: (v) => v.priceModel === 'ou' }),
    { key: 'failEvents', label: 'Sample failures as discrete events', type: 'bool', value: true, help: 'Otherwise the expected failure cost is charged every year.' },
    N('alpha', 'Confidence level for VaR and CVaR', '%', 95, 50, 99.9, ''),
    N('riskTol', 'Risk tolerance', 'M$', 300, 1, 1e6, 'Parameter of the exponential utility.'),
    N('seed', 'Random seed', '', 2026, 0, 1e9, 'Same seed, same result.', { int: true }),
  ] },
  { group: 'Decision analysis', tab: 'setup', fields: [
    { key: 'ahp', label: 'Pairwise comparison of criteria', type: 'table', help: 'Saaty scale 1–9: how much more important the row is than the column. Only the upper triangle is read; the lower one is its reciprocal.', columns: [{ key: 'name', label: 'Criterion', type: 'text' }, { key: 'c1', label: 'vs Value' }, { key: 'c2', label: 'vs Risk' }, { key: 'c3', label: 'vs Carbon' }, { key: 'c4', label: 'vs Operability' }], value: [{ name: 'Value (risked NPV)', c1: 1, c2: 2, c3: 4, c4: 3 }, { name: 'Blockage risk', c1: 0.5, c2: 1, c3: 3, c4: 2 }, { name: 'Carbon', c1: 0.25, c2: 0.333, c3: 1, c4: 0.5 }, { name: 'Operability', c1: 0.333, c2: 0.5, c3: 2, c4: 1 }] },
    N('testRel', 'Reliability of the appraisal', '%', 80, 34, 100, 'Probability that it points to the true scenario.'),
    N('testCost', 'Cost of the appraisal', 'M$', 15, 0, 1000, ''),
    N('optVol', 'Volatility of project value', '%/y', 30, 1, 150, ''),
    N('optYears', 'Deferral window', 'y', 3, 0.25, 20, ''),
    N('riskFree', 'Risk-free rate', '%/y', 4, 0, 25, ''),
    N('optYield', 'Value lost while waiting', '%/y', 5, 0, 40, 'Cash-flow leakage (dividend yield) of the undeveloped project.'),
    N('expandFrac', 'Expansion size', '% of project value', 30, 0, 300, ''),
    N('expandCost', 'Expansion cost', 'M$', 80, 0, 1e4, ''),
    N('salvage', 'Salvage on abandonment', '% of CAPEX', 25, 0, 100, ''),
    { key: 'portfolio', label: 'Candidate projects', type: 'table', help: 'Chosen to maximise NPV within the budget and the vessel days available.', columns: [{ key: 'name', label: 'Project', type: 'text' }, { key: 'capex', label: 'CAPEX', unit: 'M$' }, { key: 'npv', label: 'NPV', unit: 'M$' }, { key: 'days', label: 'Vessel days' }], value: [{ name: 'Slug-catcher upgrade', capex: 14, npv: 9, days: 0 }, { name: 'Subsea chemical-injection upgrade', capex: 22, npv: 17, days: 12 }, { name: 'Infill well', capex: 75, npv: 46, days: 30 }, { name: 'Water-injection debottleneck', capex: 30, npv: 21, days: 6 }, { name: 'Riser-base gas lift', capex: 26, npv: 20, days: 16 }, { name: 'Online wax monitoring', capex: 4, npv: 3.5, days: 0 }] },
    N('budget', 'Capital budget', 'M$', 100, 0, 1e5, ''),
    N('vesselDays', 'Vessel days available', 'd', 40, 0, 5000, ''),
  ] },
  { group: 'Resolution', tab: 'mesh', help: 'Sample sizes and grid counts. The kernel studies cost about 0.1 s for each steady solution.', fields: [
    N('nMC', 'Simulation samples', '', 2000, 100, 400000, 'Monte Carlo / Latin-hypercube sample size.', { int: true }),
    N('nLattice', 'Lattice steps', '', 120, 5, 4000, 'Time steps of the binomial lattice.', { int: true }),
    N('nSobol', 'Base sample of the Sobol study', '', 256, 16, 20000, 'Model runs = this × (inputs + 2).', { int: true }),
    N('nMCopt', 'Samples per design under uncertainty', '', 300, 50, 20000, '', { int: true }),
    N('nPop', 'Population of the Pareto search', '', 36, 8, 400, '', { int: true }),
    N('nGens', 'Generations of the Pareto search', '', 16, 2, 400, '', { int: true }),
    N('nCells', 'Cells along the line (flow kernel)', '', 40, 20, 400, '', { int: true }),
    N('nThick', 'Insulation thicknesses solved on the kernel', '', 5, 3, 8, '', { int: true }),
    N('tMaxMm', 'Largest insulation thickness studied', 'mm', 140, 20, 300, ''),
    N('nDiam', 'Diameters solved on the kernel', '', 5, 1, 7, '', { int: true }),
    N('nRates', 'Extra rates solved on the kernel', '', 1, 0, 2, '0 uses standard exponents for the rate response.', { int: true }),
    N('erosC', 'Erosional constant C', '(lb/ft³)^½·ft/s', 100, 50, 400, 'API RP 14E.'),
  ] },
];
const FIELDS = INPUTS.flatMap((g) => g.fields), FIELD = Object.fromEntries(FIELDS.map((f) => [f.key, f]));

/** Validate and normalise the raw input values; throws a plain-language Error for impossible input. */
function readInputs(v = {}) {
  const q = {};
  for (const f of FIELDS) {
    const raw = v[f.key];
    if (f.type === 'select') { q[f.key] = f.options.some((o) => o.value === raw) ? raw : f.value; continue; }
    if (f.type === 'bool') { q[f.key] = raw === undefined || raw === null ? f.value : !!raw; continue; }
    if (f.type === 'text') { q[f.key] = typeof raw === 'string' && raw.trim() ? raw.trim().slice(0, 12) : f.value; continue; }
    if (f.type === 'table') { q[f.key] = Array.isArray(raw) ? raw.filter((r) => r && typeof r === 'object') : f.value; continue; }
    const x = raw === undefined || raw === null || raw === '' ? f.value : +raw;
    if (!Number.isFinite(x)) throw new Error(`${f.label} must be a number.`);
    if (x < f.min || x > f.max) throw new Error(`${f.label} must lie between ${f.min} and ${f.max}${f.unit ? ' ' + f.unit : ''} (got ${x}).`);
    q[f.key] = f.int ? Math.round(x) : x;
  }
  if (!(q.qOil > 0) && !(q.qGas > 0)) throw new Error('Enter an oil rate or a gas rate: there is nothing to sell.');
  if (q.riserLen >= q.lineLen) throw new Error('The riser cannot be longer than the whole line.');
  if (q.discount <= -100) throw new Error('The discount rate must be greater than −100 %.');
  const pos = (x, d = 1) => (Number.isFinite(+x) && +x > 0 ? +x : d), ph = q.phasing.map((r) => +r.pct).filter((x) => Number.isFinite(x) && x >= 0).slice(0, 10), phS = sum(ph);
  if (!(phS > 0)) throw new Error('CAPEX phasing needs at least one construction year with a positive share.');
  q.phase = ph.map((x) => x / phS);
  const row = (id) => q.dists.find((r) => String(r.id ?? '').trim().toLowerCase() === id.toLowerCase());
  q.dists = Object.fromEntries(VARS.map((id) => { const r = row(id); return [id, r ? makeDist({ ...r, name: VAR_LABEL[id] }) : makeDist({ lo: 1, hi: 1 })]; }));
  for (const id of VARS) if (q.dists[id].inv(1e-6) < 0) throw new Error(`Distribution "${VAR_LABEL[id]}" reaches negative multipliers; narrow it or use a lognormal.`);
  const C = VARS.map((_, i) => VARS.map((__, j) => (i === j ? 1 : 0))), ix = (s) => VARS.findIndex((k) => k.toLowerCase() === String(s ?? '').trim().toLowerCase());
  for (const r of q.corr) { const a = ix(r.a), b = ix(r.b), rho = +r.rho; if (a >= 0 && b >= 0 && a !== b && Number.isFinite(rho)) C[a][b] = C[b][a] = clamp(rho, -0.99, 0.99); }
  q.corr = C;
  q.scenarios = q.scenarios.filter((s) => +s.weight > 0).slice(0, 7).map((s, i) => ({ name: String(s.name ?? `Scenario ${i + 1}`).slice(0, 30), weight: +s.weight, price: pos(s.price), prod: pos(s.prod), capex: pos(s.capex), opex: pos(s.opex) }));
  if (!q.scenarios.length) throw new Error('The scenario table needs at least one scenario with a positive weight.');
  q.ahpM = [0, 1, 2, 3].map((i) => [0, 1, 2, 3].map((j) => (i === j ? 1 : i < j ? pos(q.ahp[i]?.[`c${j + 1}`]) : 1 / pos(q.ahp[j]?.[`c${i + 1}`]))));
  q.portfolio = q.portfolio.filter((r) => Number.isFinite(+r.capex) && +r.capex >= 0 && Number.isFinite(+r.npv)).slice(0, 12).map((r, i) => ({ name: String(r.name ?? `Project ${i + 1}`).slice(0, 40), capex: +r.capex, npv: +r.npv, days: Math.max(+r.days || 0, 0) }));
  return q;
}

const PRESETS = [
  { name: 'Deep-water oil tie-back (reference)', values: { strategy: 'wet', oilPrice: 75, qOil: 3000 } },
  { name: 'Marginal field — low-price stress test', values: { oilPrice: 48, gasPrice: 2.5, qOil: 1400, qGas: 250000, plateau: 1, Di: 22, nWells: 1, wellCost: 85, discount: 12, abandon: 45, uptime: 92, hurdle: 15, priceModel: 'gbm' } },
  { name: 'Gas-condensate export with MEG loop', values: { strategy: 'bare', inhibitor: 'MEG', qOil: 700, qGas: 4500000, gasSalesFrac: 96, gasPrice: 6.5, wc0: 3, wcEnd: 15, plateau: 6, Di: 10, declineType: 'hyp', bHyp: 0.4, compKW: 9000, powerSource: 'gas', tariff: 1.2, life: 25 } },
  { name: 'Brownfield late life with integrity spend', values: { life: 10, phasing: [{ year: 1, pct: 100 }], qOil: 1500, qGas: 200000, wc0: 55, wcEnd: 88, plateau: 0, Di: 12, wellCost: 0, capexSunk: 300, residual: 40, pof: 0.01, remLife: 8, assetAge: 18, corrRate: 0.35, inspInterval: 3, maintPct: 4, uptime: 90, abandon: 80, consequence: 220, corrMgmt: 2.5, deprLife: 5, pEnd: 60 } },
  { name: 'Production-sharing contract', values: { regime: 'psc', royalty: 5, costOilCap: 60, profitSplit: 35, taxRate: 30, deprMethod: 'uop', gearing: 50, loanRate: 8, loanTenor: 7 } },
  { name: 'High carbon price', values: { carbonPrice: 150, carbonEsc: 5, powerSource: 'gas', flareFrac: 2, compKW: 4000, blowdowns: 3, deprMethod: 'db', mid: true } },
];

function pull({ fluid, outputs } = {}) {
  const o = outputs || {}, items = [];
  const add = (key, value, from) => {
    const f = FIELD[key];
    if (!f || items.some((i) => i.key === key) || value === undefined || value === null) return;
    if (f.type === 'select') { if (f.options.some((x) => x.value === value)) items.push({ key, value, from }); return; }
    if (typeof value !== 'number' || !Number.isFinite(value)) return;
    const x = clamp(value, f.min, f.max);
    items.push({ key, value: f.int ? Math.round(x) : x, from });
  };
  const oil = !fluid?.rateBasis || fluid.rateBasis === 'oil';
  add('qOil', o.pvt?.rates?.qOilStd, 'Fluid suite: stock-tank oil rate'); add('qOil', oil ? fluid?.qOil : undefined, 'Case fluid: oil rate');
  add('qGas', o.pvt?.rates?.qGasStd, 'Fluid suite: produced gas rate'); add('qGas', fluid?.rateBasis === 'gas' && isNum(fluid?.qGas) ? fluid.qGas * 1e6 : undefined, 'Case fluid: gas rate');
  add('wc0', oil ? fluid?.wc : undefined, 'Case fluid: water cut');
  add('inhibitor', fluid?.inhibitor === 'MEG' ? 'MEG' : fluid?.inhibitor === 'MeOH' ? 'MeOH' : undefined, 'Case fluid: inhibitor');
  add('lineLen', o.net?.length, 'Network suite: line length'); add('riserLen', isNum(o.net?.riserHeight) && isNum(o.net?.length) ? Math.min(o.net.riserHeight * 1.08, 0.9 * o.net.length) : undefined, 'Network suite: riser height × 1.08');
  add('idMm', isNum(o.net?.id) ? o.net.id * 1000 : undefined, 'Network suite: inner diameter'); add('wtMm', isNum(o.net?.wt) ? o.net.wt * 1000 : undefined, 'Network suite: wall thickness');
  add('depth', o.net?.waterDepth, 'Network suite: water depth'); add('slugVol', o.net?.slugCatcherVol, 'Network suite: slug-catcher volume'); add('slugVol', isNum(o.flow?.slug?.surge) && o.flow.slug.surge > 0 ? o.flow.slug.surge * 1.5 : undefined, 'Flow suite: slug surge volume × 1.5');
  add('pumpKW', o.net?.pumpPower, 'Network suite: pump power'); add('compKW', o.net?.compressorPower, 'Network suite: compressor power'); add('pi', o.net?.ipr?.pi, 'Network suite: productivity index');
  add('chokeDp', isNum(o.net?.whp) && isNum(o.flow?.pIn) && o.net.whp > o.flow.pIn ? o.net.whp - o.flow.pIn : undefined, 'Wellhead pressure − flowline inlet pressure');
  add('uptime', isNum(o.ops?.uptime) ? o.ops.uptime * 100 : undefined, 'Operations suite: uptime'); add('inhibRate', o.ops?.inhibitorRate, 'Operations suite: inhibitor rate'); add('inhibRate', o.solids?.inhibitorRate, 'Solids suite: inhibitor rate');
  add('shutdowns', o.ops?.eventsPerYear?.shutdowns, 'Operations suite: shutdowns a year'); add('blowdowns', o.ops?.eventsPerYear?.blowdowns, 'Operations suite: blowdowns a year');
  add('pigRuns', o.ops?.eventsPerYear?.pigRuns, 'Operations suite: pig runs a year'); add('pigRuns', isNum(o.solids?.piggingInterval) && o.solids.piggingInterval > 0 ? 365 / o.solids.piggingInterval : undefined, 'Solids suite: pigging interval');
  add('cooldownBase', o.ops?.cooldownTime, 'Operations suite: cooldown time'); add('heatKW', o.ops?.heatingPower, 'Operations suite: heating power');
  add('capacityFrac', o.ops?.envelope?.qMax, 'Operations suite: upper rate limit'); add('qMinFrac', o.ops?.envelope?.qMin, 'Operations suite: lower rate limit'); add('qMinFrac', o.flow?.turndownRate, 'Flow suite: minimum stable rate');
  add('plugProb', o.solids?.plugProbability, 'Solids suite: plug probability'); add('wat', o.solids?.wat, 'Solids suite: wax appearance temperature'); add('wat', o.pvt?.wat, 'Fluid suite: wax appearance temperature'); add('waxRate', o.solids?.waxRate, 'Solids suite: wax build-up rate');
  add('pof', o.integ?.pof, 'Integrity suite: annual probability of failure'); add('consequence', isNum(o.integ?.consequence) ? o.integ.consequence / MM : undefined, 'Integrity suite: consequence of failure');
  add('ealOverride', isNum(o.integ?.riskCostPerYear) ? o.integ.riskCostPerYear / MM : undefined, 'Integrity suite: risk cost per year'); add('remLife', o.integ?.remainingLife, 'Integrity suite: remaining life');
  add('inspInterval', o.integ?.inspectionInterval, 'Integrity suite: inspection interval'); add('corrRate', o.integ?.corrosionRate, 'Integrity suite: corrosion rate');
  add('mawp', o.integ?.mawp, 'Integrity suite: maximum allowable working pressure'); add('mawp', o.net?.designPressure, 'Network suite: design pressure');
  { const u = [o.integ?.hoopUtil, o.integ?.vmUtil, o.integ?.collapseUtil, o.integ?.upheavalUtil].filter(isNum); add('integUtil', u.length ? Math.max(...u) : undefined, 'Integrity suite: largest utilisation'); }
  add('integViol', o.integ?.violations, 'Integrity suite: failed code checks');
  { const e = [o.flow?.erosionalRatio, o.integ?.erosionalRatio].filter(isNum); add('erosIn', e.length ? Math.max(...e) : undefined, 'Flow / integrity suite: erosional velocity ratio'); }
  if (typeof o.flow?.severeSlugging === 'boolean') items.push({ key: 'severeSlug', value: o.flow.severeSlugging && !o.ops?.slugSuppressed, from: 'Flow and operations suites: severe slugging not suppressed' });
  return items;
}
function siteHook(site) {
  const d = site?.data || {}, items = [], add = (key, value, from) => { const f = FIELD[key]; if (f.type === 'text') { if (typeof value === 'string' && value.trim()) items.push({ key, value: value.trim().slice(0, 12), from }); return; } if (isNum(value)) items.push({ key, value: clamp(value, f.min, f.max), from }); };
  add('inflation', d.inflation, 'Inflation at the site'); add('discount', isNum(d.lendingRate) ? d.lendingRate + 2 : undefined, 'Lending rate + 2 points as a discount-rate suggestion'); add('loanRate', d.lendingRate, 'Lending rate');
  add('fx', d.fxPerUSD, 'Exchange rate per US$'); add('currency', d.currency, 'Local currency'); add('elecPrice', d.electricityPrice, 'Electricity price'); add('gridCarbon', d.gridCarbon, 'Grid carbon intensity');
  add('oilPrice', d.oilPrice, 'Oil price'); add('gasPrice', d.gasPrice, 'Gas price'); add('carbonPrice', d.carbonPrice, 'Carbon price'); add('taxRate', d.taxRate, 'Corporate tax rate'); add('depth', d.depth, 'Water depth at the site');
  return items;
}

const CAL_TRUE = { cerCoef: 0.52, cerExp: 1.22, layFactor: 1.18, learnRate: 93 };
const calModel = (v) => ({ cost: pipelineCER({ dIn: v.calD ?? 10, lengthKm: v.calL ?? 20, depth: v.calDepth ?? 1000, unitNo: v.calUnit ?? 1, cerCoef: v.cerCoef ?? 0.4, cerExp: v.cerExp ?? 1.3, layFactor: v.layFactor ?? 1, learnRate: (v.learnRate ?? 90) / 100, vesselRate: v.vesselRate ?? 350, layRate: v.layRate ?? 2.5, mobCost: v.mobCost ?? 6, depthCoef: v.depthCoef ?? 0.2 }) });
// catalogue items that the engine computes (full item names, lower case) and those it does not
const IMPLEMENTED = ["present-value equation","future-value equation","compound-interest equation","continuous-compounding equation","discount-factor equation","annuity equation","perpetuity equation","discounted cash flow","cumulative cash flow","free cash flow","operating cash flow","after-tax cash flow","net present value","internal rate of return","modified internal rate of return","return on investment","return on capital employed","profitability index","discounted payback period","simple payback period","equivalent annual value","equivalent annual cost","equipment-cost scaling equations","capacity-factor/scaling-law model","six-tenths-rule-type scaling","installation factors","lang-factor methodology","bare-module costing","pipeline cost-per-length models","subsea installation cost models","vessel/day-rate calculations","compressor/pump costing","insulation costing","chemical-injection-system capex","slug-catcher sizing/cost relations","energy-consumption cost","pumping cost","compression cost","chemical/inhibitor cost","meg/methanol consumption cost","heating cost","pigging cost","inspection cost","maintenance cost","corrosion-management cost","hydrate-remediation cost","vessel/intervention cost","production-operations cost","production-revenue equation","oil/gas price models","production-decline models","cumulative production","uptime/availability","production-efficiency equation","deferred-production calculation","lost-production calculation","expected failure cost","expected annual loss","probability × consequence formulation","expected monetary value","lifecycle failure cost","intervention-cost model","downtime-cost model","monte carlo simulation","latin-hypercube sampling","probability distributions","expected-value analysis","variance","value at risk","conditional value at risk","stochastic cash-flow modelling","scenario analysis","sensitivity analysis","tornado analysis","decision trees","bayesian decision analysis","utility theory","multi-criteria decision analysis","analytic hierarchy process","topsis where appropriate","real-options analysis","linear programming","nonlinear programming","mixed-integer linear programming","mixed-integer nonlinear programming","dynamic programming","stochastic programming","robust optimization","multi-objective optimization","pareto-front optimization","genetic algorithms","particle-swarm optimization","hydrate risk + economics","slugging + economics","thermal hydraulics + economics","integrity + economics","reliability + economics","production + flow assurance + economics","physics + economics + uncertainty","historical capex","installation campaign costs","vessel/day rates","cost-estimating relationships","parametric cost regression","learning curves","historical back-fitting","npv analytical benchmarks","irr benchmark cases","mirr verification","discount-factor verification","nominal-vs-real cash-flow consistency","inflation calculations","tax calculations","depreciation schedules","royalty calculations","working-capital calculations","escalation calculations","capex phasing","opex aggregation","revenue calculations","production-decline integration","unit conversion","currency conversion","cash-flow sign conventions","payback calculation","breakeven root finding","probability-weighted cash flows","monte-carlo convergence","latin-hypercube sampling convergence","correlation-matrix handling","probability-distribution sampling tests","sensitivity calculations","tornado-chart calculations","scenario-weight normalization","decision-tree arithmetic","predicted capex vs actual capex","project evaluation date","project life","remaining asset life","base currency","exchange rates where multiple currencies are involved","initial commodity prices","initial production rates","initial water and gas handling rates","remaining equipment value","existing hydrate/wax/scale management strategy","initial capex already committed","baseline electricity or fuel prices","inhibitor and chemical prices","pigging costs","inspection and maintenance costs","logistics and offshore-vessel costs","disposal and treatment costs","existing tariffs","transportation charges","taxes","royalties and other fiscal assumptions relevant to the project","project evaluation horizon","discount rate","inflation assumptions","escalation rates","commodity-price scenarios","exchange-rate assumptions","tax and royalty structure","financing assumptions where considered","production limits","equipment capacity constraints","chemical availability and maximum injection capacity","emissions or carbon costs where applicable","abandonment/decommissioning obligations","minimum economic return requirements","decision thresholds such as minimum npv","maximum acceptable payback period or required rate of return","base and high cases","probability distributions or scenario ranges for production","commodity prices","capex","opex","downtime","failure frequency","repair costs and other uncertain economic drivers","engineering outputs and constraints from modules 1-6","project/evaluation life","base currency and exchange-rate assumptions","production/revenue forecasts","fixed/variable opex","energy prices","inhibitor/chemical costs","heating, compression and pumping costs","pigging, inspection, maintenance, intervention and repair costs","downtime and deferred-production assumptions","taxes/royalties","discount and inflation/escalation rates","equipment/chemical availability","failure frequencies and consequences","uncertainty distributions/scenarios","decision and return thresholds","capex/opex breakdowns","chemical and energy expenditure","maintenance/intervention and remediation costs","production revenue, deferment and downtime losses","expected failure/risk cost","lifecycle cost and unit production cost","cash-flow profiles","npv","irr","payback and discounted payback","profitability/break-even metrics","sensitivity and uncertainty results","scenario comparisons","pareto/decision metrics where applicable","economically optimized design, operating, mitigation, inspection and intervention strategies subject to engineering safety constraints"];
const REFERENCE_ONLY = ["bayesian optimization","actual procurement costs","epc cost data","actual chemical consumption","historical chemical prices","electricity/fuel consumption","energy tariffs","actual maintenance expenditures","historical downtime","production deferment","equipment availability","actual production profiles","commodity-price histories","inflation indices","escalation indices","tax/royalty histories","decommissioning costs","econometric regression","maximum-likelihood estimation","bayesian calibration","time-series calibration","cost-index normalization","location-factor calibration","predicted opex vs actual opex","forecast vs actual production","predicted vs actual chemical consumption","predicted vs actual energy consumption","predicted vs actual maintenance","predicted vs actual intervention frequency","predicted vs actual downtime","predicted vs actual production losses","forecast vs actual project cash flow","forecast vs actual abandonment costs","existing installed equipment and infrastructure","existing maintenance state","available chemical inventory","initial working capital and the starting economic scenario","labour costs","electricity and fuel-price scenarios","contractual sales limits","maintenance and vessel availability","planned turnaround periods","regulatory constraints"];

/** Verification: textbook cases with hand or closed-form answers, sampling tests against analytic moments, and solver benchmarks. */
function verify() {
  const out = [], chk = (name, expected, got, tol, note = '') => out.push({ name, expected, got, tol, pass: typeof got === 'number' && Number.isFinite(got) && Math.abs(got - expected) <= tol, note });
  const flag = (name, cond, got, expected, note = '') => out.push({ name, expected, got, tol: 0, pass: !!cond, note });
  // a minimal project: every fiscal term off unless overridden
  const mini = (o = {}) => ({ phase: [1], capex: 100, capexSunk: 0, residual: 0, life: 3, oil: [1, 1, 1], gas: [0, 0, 0], water: [0, 0, 0], uptime: 1, deferFrac: 0, oilPrice: 60, gasPrice: 0, infl: 0, costEsc: 0, priceEsc: 0, carbonEsc: 0, discount: 0.1, mid: false, opexFixed: 0, opexVarBoe: 0, waterCost: 0, opexDown: 0, opexBlock: 0, carbonT: 0, carbonPrice: 0, includeRisk: false, consequence: 0, haz: [0, 0, 0], royalty: 0, taxRate: 0, regime: 'tax', costOilCap: 0.5, profitSplit: 0.4, deprFrac: [1, 0, 0], wcDays: 0, abandon: 0, abandonProvision: false, gearing: 0, loanRate: 0, loanTenor: 1, ...o });
  // --- time value of money
  chk('NPV of a 10-year annuity vs closed form', 671.0081398941447, npv(0.08, [0, ...zeros(10).map(() => 100)]), 1e-9, 'A·(1 − (1+r)^−n)/r with A = 100, r = 8 %');
  chk('Annuity formula', 671.0081398941447, annuityPV(100, 0.08, 10), 1e-9, 'closed form');
  chk('Growing perpetuity vs a 3000-year sum', 1666.6666666666667, npv(0.08, [0, ...zeros(3000).map((_, k) => 100 * 1.02 ** k)]), 1e-6, 'A/(r − g) = 100/0.06');
  chk('Perpetuity formula', 1666.6666666666667, perpetuityPV(100, 0.08, 0.02), 1e-9, 'A/(r − g)');
  chk('Future value, annual compounding', 1628.8946267774422, futureValue(1000, 0.05, 10), 1e-9, '1000 × 1.05^10');
  chk('Compound interest, monthly', 1126.8250301319697, futureValue(1000, 0.12, 1, 12), 1e-9, '1000 × 1.01^12');
  chk('Present value inverts future value', 1000, presentValue(1628.8946267774422, 0.05, 10), 1e-9, '');
  chk('Continuous compounding', 1648.7212707001281, continuousFV(1000, 0.1, 5), 1e-9, '1000·e^0.5');
  chk('Discount factor, end of year 5 at 10 %', 0.6209213230591549, discountFactor(0.1, 5), 1e-12, '1.1^−5');
  chk('Discount factor, mid-year convention', 0.6512277776419588, discountFactor(0.1, 5, { mid: true }), 1e-12, '1.1^−4.5');
  chk('Capital-recovery factor', 0.26379748079474524, capitalRecovery(0.1, 5), 1e-12, 'r/(1 − (1+r)^−n)');
  // --- investment metrics
  chk('IRR of −100 now, +121 in two years', 0.1, irr([-100, 0, 121]).irr, 1e-10, 'analytic: √1.21 − 1');
  const multi = irr([-100, 230, -132]);
  chk('Multiple IRR: both roots found (10 % and 20 %)', 0.3, sum(multi.roots), 1e-9, 'NPV = −100 + 230/(1+r) − 132/(1+r)² has roots 0.10 and 0.20');
  flag('Multiple sign changes are flagged', multi.multiple && multi.signChanges === 2, multi.signChanges, 2, 'two sign changes in the cash-flow series');
  flag('No IRR when the cash flow never changes sign', irr([10, 20, 30]).irr === null, 'null', 'null');
  chk('MIRR hand calculation', 0.09815669244631531, mirr([-1000, 300, 400, 500], 0.1, 0.12), 1e-12, '(1324.32/1000)^(1/3) − 1');
  chk('Simple payback, interpolated', 2.5, payback([-100, 40, 40, 40]), 1e-12, '2 y + 20/40');
  chk('Discounted payback, interpolated', 1.9166666666666667, payback([-100, 60, 60], 0.1), 1e-10, '1 y + 45.4545/49.5868');
  // --- nominal vs real (Fisher)
  { const infl = 0.03, rr = 0.05, rn = (1 + rr) * (1 + infl) - 1; chk('Nominal flows at the nominal rate = real flows at the real rate', 432.9476670630823, npv(rn, [0, ...zeros(5).map((_, k) => 100 * (1 + infl) ** (k + 1))]), 1e-9, 'Fisher: (1+r_n) = (1+r_r)(1+i); closed-form real annuity'); chk('Real rate from the Fisher relation', 0.05, realRate(rn, infl), 1e-12, ''); }
  { const p = mini({ infl: 0.03, opexFixed: 10, life: 5, oil: zeros(5).map(() => 1), gas: zeros(5), water: zeros(5), haz: zeros(5), deprFrac: [1, 0, 0, 0, 0] }), cf = cashflow(p), m = investmentMetrics(cf, p); chk('Cash-flow engine: deflated flows at the real rate reproduce the NPV', cf.npv, m.realNpv, 1e-9 * Math.abs(cf.npv), 'engine identity on an inflated case'); }
  chk('Cost escalation and inflation compound', 10.93363663608, cashflow(mini({ infl: 0.02, costEsc: 0.01, opexFixed: 10 })).opex[3], 1e-9, '10 × (1.02 × 1.01)³');
  // --- depreciation, tax, royalty, phasing, working capital
  chk('Straight-line schedule sums to the base', 100, sum(depreciation('sl', 100, 4, { years: 6 })), 1e-12, '4 × 25');
  { const d = depreciation('db', 100, 5); chk('Declining-balance schedule sums to the base', 100, sum(d), 1e-12, '40, 24, 14.4, then straight line 10.8, 10.8'); chk('Declining-balance third-year charge', 14.4, d[2], 1e-12, '100 × 0.6² × 0.4'); }
  chk('Units-of-production depreciation follows production', 50, depreciation('uop', 100, 3, { units: [5, 3, 2] })[0], 1e-12, '100 × 5/10');
  { const cf = cashflow(mini({ taxRate: 0.3 })); chk('Tax with loss carry-forward (hand case)', 24, sum(cf.tax), 1e-9, 'taxable −40, 60, 60 → tax 0, 6, 18'); chk('Loss carried forward shelters year-2 income', 6, cf.tax[2], 1e-9, '0.3 × (60 − 40)'); chk('Sign convention: year-0 free cash flow is −CAPEX', -100, cf.fcf[0], 1e-12, 'payments negative, receipts positive'); chk('Σ free cash flow = revenue − tax − CAPEX', 56, sum(cf.fcf), 1e-9, '180 − 24 − 100'); }
  chk('Royalty on gross revenue', 22.5, sum(cashflow(mini({ royalty: 0.125 })).royalty), 1e-9, '12.5 % × 180');
  { const cf = cashflow(mini({ phase: [0.3, 0.5, 0.2], capex: 200 })); chk('CAPEX phasing sums to the total', 200, sum(cf.capex), 1e-9, '60 + 100 + 40'); chk('CAPEX phasing: second-year tranche', 100, cf.capex[1], 1e-9, '50 % of 200'); }
  { const cf = cashflow(mini({ wcDays: 36.5 })); chk('Working capital builds in the first production year', 6, cf.dwc[1], 1e-9, '36.5/365 × 60'); chk('Working capital is fully released', 0, sum(cf.dwc), 1e-9, 'Σ Δ = 0'); }
  chk('Production-sharing contract (hand case)', 3, cashflow(mini({ regime: 'psc', royalty: 0.1, capex: 0, opexFixed: 60, oilPrice: 100, life: 1, oil: [1], gas: [0], water: [0], haz: [0], deprFrac: [1] })).ocf[1], 1e-9, 'net 90, cost oil 45 (cap), profit oil 45, contractor 18: 45 + 18 − 60');
  chk('OPEX aggregation: fixed + variable + water', 10 + 2 * 3 + 4 * 5, cashflow(mini({ opexFixed: 10, opexVarBoe: 2, oil: [3, 3, 3], waterCost: 4, water: [5, 5, 5] })).opex[1], 1e-9, '10 + 2 × 3 + 4 × 5');
  chk('Revenue = volume × price × uptime', 0.9 * 7 * 60 + 0.9 * 100 * 4, cashflow(mini({ uptime: 0.9, oil: [7, 7, 7], gas: [100, 100, 100], gasPrice: 4 })).revenue[1], 1e-9, 'oil and gas streams');
  chk('Expected annual loss = probability × consequence', 1.5, sum(cashflow(mini({ includeRisk: true, consequence: 50, haz: [0.01, 0.01, 0.01] })).risk), 1e-12, '3 × 0.01 × 50');
  { const cf = cashflow(mini({ uptime: 0.9, deferFrac: 1, oil: [10, 10, 10] })); chk('Deferred production is recovered in the last year', 30, sum(cf.oil), 1e-9, 'nothing is lost when 100 % is deferred'); }
  // --- production decline
  chk('Arps exponential: cumulative vs analytic integral', 4323.323583816937, sum(productionProfile({ q0: 1000, plateau: 0, Di: 0.2, b: 0, life: 10 })), 1e-8, 'q_i (1 − e^(−Dt))/D');
  chk('Arps hyperbolic (b = 0.5): cumulative', 5000, arps(1000, 0.2, 0.5, 10).Np, 1e-9, 'q_i/((1−b)D) (1 − (1+bDt)^((b−1)/b))');
  chk('Arps harmonic: cumulative', 5493.061443340548, arps(1000, 0.2, 1, 10).Np, 1e-9, '(q_i/D) ln(1 + Dt)');
  { const n = 20000, h = 10 / n; let s = 0; for (let i = 0; i < n; i++) s += 0.5 * h * (arps(1000, 0.2, 0.5, i * h).q + arps(1000, 0.2, 0.5, (i + 1) * h).q); chk('Arps hyperbolic: cumulative equals the numerical integral of the rate', s, arps(1000, 0.2, 0.5, 10).Np, 1e-4, 'trapezoidal rule, 20,000 steps'); }
  // --- break-even root finding
  { const p = mini({ capex: 1000, life: 5, oil: zeros(5).map(() => 100), gas: zeros(5), water: zeros(5), haz: zeros(5), deprFrac: [1, 0, 0, 0, 0], opexFixed: 20 }), x = breakeven((pr) => cashflow({ ...p, oilPrice: pr }, {}, false), 0, 100); chk('Break-even price vs closed form', 2.8379748079474525, x, 1e-9, '(I × capital-recovery factor + C)/Q'); chk('Break-even root gives NPV = 0', 0, cashflow({ ...p, oilPrice: x }, {}, false), 1e-7, 'residual of the root'); }
  // --- decision analysis
  { const t = decisionTree({ name: 'root', type: 'decision', branches: [{ name: 'A', node: { name: 'A', type: 'chance', branches: [{ name: 'up', p: 0.3, value: 100 }, { name: 'down', p: 0.7, value: -20 }] } }, { name: 'B', value: 10 }] }); chk('Decision-tree roll-back (EMV)', 16, t.emv, 1e-12, '0.3 × 100 − 0.7 × 20 = 16 beats 10'); flag('Decision node picks the larger EMV', t.choice === 'A', t.choice, 'A'); }
  chk('Chance-node weights are normalised', 0.5 * 40 + 0.25 * 0 + 0.25 * -20, decisionTree({ name: 's', type: 'chance', branches: [{ name: 'a', p: 2, value: 40 }, { name: 'b', p: 1, value: 0 }, { name: 'c', p: 1, value: -20 }] }).emv, 1e-12, 'weights 2 : 1 : 1 → 0.5, 0.25, 0.25');
  chk('Branch cost is charged in the roll-back', 4, decisionTree({ name: 'r', type: 'decision', branches: [{ name: 'test', cost: 6, value: 10 }, { name: 'none', value: 3 }] }).emv, 1e-12, '10 − 6 = 4 > 3');
  { const r = valueOfInformation({ prior: [0.4, 0.6], payoff: [[100, -50], [20, 20]], likelihood: [[0.8, 0.2], [0.2, 0.8]] }); chk('Value of perfect information (hand case)', 32, r.evpi, 1e-12, '0.4 × 100 + 0.6 × 20 − 20'); chk('Value of imperfect information (Bayesian revision)', 17.2, r.evii, 1e-12, 'signals 0.44 / 0.56, posteriors 0.727 / 0.143'); flag('0 ≤ EVII ≤ EVPI', r.evii >= 0 && r.evii <= r.evpi, r.evii, '≤ 32'); }
  chk('Certainty equivalent under exponential utility', 37.988549304172246, certaintyEquivalent([0, 100], 100), 1e-10, '−R ln(½(1 + e^−1))');
  { const r = ahp([[1, 2, 4], [0.5, 1, 2], [0.25, 0.5, 1]]); chk('AHP: consistency ratio of a consistent matrix', 0, r.cr, 1e-10, 'λmax = n'); chk('AHP: principal eigenvector', 4 / 7, r.weights[0], 1e-10, 'weights 4/7, 2/7, 1/7'); }
  { const r = topsis([[10, 1], [5, 2], [7, 1.5]], [0.5, 0.5], [true, false]); chk('TOPSIS: dominating alternative has closeness 1', 1, r.closeness[0], 1e-12, ''); chk('TOPSIS: dominated alternative has closeness 0', 0, r.closeness[1], 1e-12, ''); }
  chk('Weighted sum with min–max normalisation', 0.5, weightedSum([[10, 1], [5, 2], [7.5, 1.5]], [1, 1], [true, false]).scores[2], 1e-12, 'mid-point alternative scores 0.5');
  flag('Minimax regret (hand case)', minimaxRegret([[10, 0], [6, 5]]).best === 1, minimaxRegret([[10, 0], [6, 5]]).best, 1, 'largest regrets 5 and 4');
  // --- options
  { const o = { S: 100, K: 100, r: 0.05, sigma: 0.2, T: 1 }, bs = blackScholes(o); chk('Black–Scholes call (textbook case)', 10.450583572185565, bs, 1e-9, 'S = K = 100, r = 5 %, σ = 20 %, T = 1'); chk('Binomial lattice converges to Black–Scholes', bs, binomialOption({ ...o, steps: 2000 }).value, 3e-3, '2000 steps'); chk('Put–call parity', 100 - 100 * Math.exp(-0.05), bs - blackScholes({ ...o, type: 'put' }), 1e-10, 'C − P = S − K e^(−rT)'); chk('American put (lattice benchmark)', 6.09, binomialOption({ ...o, steps: 1000, type: 'put', american: true }).value, 0.01, 'reference value 6.090'); }
  // --- sampling
  { const d = [{ dist: 'triangular', lo: 0, mode: 1, hi: 3 }, { dist: 'pert', lo: 0, mode: 1, hi: 3 }, { dist: 'lognormal', lo: 0.5, hi: 2 }].map(makeDist), X = sampleCorrelated(d, { n: 3000, method: 'lhs', seed: 5 }).X, c = (j) => X.map((r) => r[j]);
    chk('Triangular sample mean', 4 / 3, mean(c(0)), 2e-3, '(a + m + b)/3'); chk('Triangular sample variance', 7 / 18, variance(c(0)), 4e-3, '(a² + b² + m² − ab − am − bm)/18'); chk('PERT sample mean', 7 / 6, mean(c(1)), 2e-3, '(a + 4m + b)/6'); chk('PERT sample variance', (7 / 6) * (11 / 6) / 7, variance(c(1)), 4e-3, '(μ − a)(b − μ)/7');
    chk('Lognormal sample mean', Math.exp(0.5 * (Math.log(2) / Z90) ** 2), mean(c(2)), 4e-3, 'exp(μ + σ²/2) with P10 = 0.5, P90 = 2'); chk('Lognormal P90', 2, quantile(c(2), 0.9), 0.02, 'defined by its P10 and P90'); }
  { const n = 200, u = lhs(n, 1, 3).map((r) => Math.floor(r[0] * n)); chk('Latin hypercube: one sample in every stratum', n, new Set(u).size, 0, '200 strata'); }
  { const C = [[1, 0.6, 0.2], [0.6, 1, -0.3], [0.2, -0.3, 1]], L = cholesky(C); let e = 0; for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) { let s = 0; for (let k = 0; k < 3; k++) s += L[i][k] * L[j][k]; e = Math.max(e, Math.abs(s - C[i][j])); } chk('Cholesky factor reproduces the correlation matrix', 0, e, 1e-14, 'max |L·Lᵀ − C|');
    const d = [0, 1, 2].map(() => makeDist({ dist: 'uniform', lo: 0, hi: 1 })), target = (6 / Math.PI) * Math.asin(0.3);
    for (const method of ['mc', 'lhs']) { const X = sampleCorrelated(d, { n: 4000, method, corr: C, seed: 8 }).X; chk(`Correlated sampling (${method === 'mc' ? 'Gaussian copula' : 'Latin hypercube, Iman–Conover'}): rank correlation`, target, correlation(ranks(X.map((r) => r[0])), ranks(X.map((r) => r[1]))), 0.03, 'Spearman = (6/π) asin(ρ/2) for ρ = 0.6'); } }
  { const sdOf = (n, lh) => { const m = []; for (let k = 0; k < 60; k++) { const R = rng(1000 + 17 * k + n); m.push(lh ? mean(lhs(n, 1, 1000 + 17 * k + n).map((r) => r[0])) : mean(Array.from({ length: n }, () => R.uniform()))); } return Math.sqrt(variance(m)); }, s1 = sdOf(100), s2 = sdOf(1600);
    chk('Monte Carlo error falls as 1/√N', -0.5, Math.log(s2 / s1) / Math.log(16), 0.12, 'slope of log(standard error) between N = 100 and 1600, 60 replicates'); chk('Standard error of the mean of U(0,1), N = 1600', Math.sqrt(1 / 12 / 1600), s2, 0.0025, 'σ/√N');
    flag('Latin hypercube converges faster than Monte Carlo', sdOf(100, true) < 0.2 * s1, sdOf(100, true), `< ${(0.2 * s1).toExponential(2)}`, 'standard error of the mean at N = 100'); }
  { const y = lhs(20000, 1, 4).map((r) => 100 + 50 * normInv(r[0])), s = riskStats(y, 0.95); chk('Value at risk of a normal outcome', 17.75731865242639, s.var, 0.2, 'μ − 1.645σ'); chk('Conditional value at risk of a normal outcome', -3.1356403753713806, s.cvar, 0.3, 'μ − σ φ(z)/(1 − α)'); chk('Sample variance of a normal outcome', 2500, s.variance, 15, 'σ²'); }
  chk('Normal distribution function', 0.9750021048517795, normCdf(1.96), 1e-12, 'Φ(1.96)');
  chk('Inverse normal', 1.959963984540054, normInv(0.975), 1e-9, 'Φ⁻¹(0.975)');
  { const R = rng(9); let s = 0; const n = 20000; for (let i = 0; i < n; i++) s += pricePath('gbm', 6, R, { sigma: 0.25 })[5]; chk('Geometric Brownian price path keeps its expected value', 1, s / n, 0.02, 'E[exp(σW − σ²t/2)] = 1 at t = 5'); }
  // --- sensitivities
  { const t = tornado((x) => 5 * x[0] + x[1] + 3 * x[2], [0, 0, 0], [{ name: 'a', lo: -1, hi: 1 }, { name: 'b', lo: -1, hi: 1 }, { name: 'c', lo: -1, hi: 1 }]); flag('Tornado bars are ordered by swing', t.map((x) => x.name).join('') === 'acb', t.map((x) => x.name).join(' > '), 'a > c > b'); chk('Tornado swing of the largest driver', 10, t[0].swing, 1e-12, '5 × (1 − (−1))'); }
  { const s = sobolIndices((u) => u[0] + 2 * u[1] + 3 * u[2], 3, 4096, 2); chk('Sobol first-order index of an additive model', 9 / 14, s.first[2], 0.04, 'a_i²/Σa² = 9/14'); chk('Sobol indices of an additive model sum to one', 1, sum(s.first), 0.06, ''); }
  { const R = rng(4), X = Array.from({ length: 400 }, () => [R.normal(), R.normal()]), r = standardisedRegression(X, X.map((x) => 3 * x[0] + 4 * x[1])); chk('Standardised regression coefficients', 0.8, r.src[1], 0.03, '4/√(3² + 4²) for independent unit-variance inputs'); }
  // --- optimisation
  { const r = simplex({ c: [3, 5], A: [[1, 0], [0, 2], [3, 2]], b: [4, 12, 18] }); chk('Simplex: textbook LP optimum', 36, r.value, 1e-9, 'max 3x + 5y; x ≤ 4, 2y ≤ 12, 3x + 2y ≤ 18 → (2, 6)'); chk('Simplex: optimal vertex', 6, r.x[1], 1e-9, 'y = 6'); }
  chk('Simplex with ≥ constraints (two-phase)', 13, simplex({ c: [2, 3], A: [[1, 1], [1, 0], [0, 1]], b: [10, 2, 3], types: ['<=', '>=', '>='], maximize: false }).value, 1e-9, 'min 2x + 3y; x ≥ 2, y ≥ 3');
  flag('Simplex detects an infeasible problem', simplex({ c: [1, 1], A: [[1, 1], [1, 1]], b: [1, 3], types: ['<=', '>='] }).status === 'infeasible', simplex({ c: [1, 1], A: [[1, 1], [1, 1]], b: [1, 3], types: ['<=', '>='] }).status, 'infeasible');
  { const r = branchAndBound({ c: [60, 100, 120], A: [[10, 20, 30], [1, 0, 0], [0, 1, 0], [0, 0, 1]], b: [50, 1, 1, 1] }); chk('Branch and bound: 0/1 knapsack', 220, r.value, 1e-9, 'items 2 and 3; the LP relaxation gives 240'); chk('Knapsack LP relaxation bound', 240, r.relaxation, 1e-9, ''); }
  chk('Dynamic programming: replacement timing (hand case)', 55, replacementDP({ horizon: 3, maxAge: 3, opCost: (a) => [10, 20, 40][a], replaceCost: 15 }).cost, 1e-12, 'keep, keep, replace: 10 + 20 + 25');
  { const r = twoStage({ c: 1, q: 2, scenarios: [{ p: 0.3, demand: 100 }, { p: 0.4, demand: 200 }, { p: 0.3, demand: 300 }] }); chk('Two-stage stochastic programme: first-stage capacity', 200, r.x, 1e-9, 'newsvendor critical fractile 1 − c/q = 0.5'); chk('Two-stage stochastic programme: expected cost', 260, r.cost, 1e-9, '200 + 2 × 0.3 × 100'); }
  { const f = (x) => (x[0] - 1) ** 2 + 2 * (x[1] + 0.5) ** 2, lo = [-2, -2], hi = [2, 2]; chk('Projected gradient reaches the minimum', 1, gradientDescent(f, [0, 0], { lo, hi }).x[0], 1e-4, 'quadratic bowl at (1, −0.5)'); chk('Genetic algorithm reaches the minimum', 0, geneticAlgorithm(f, lo, hi, { pop: 30, gens: 40 }).f, 1e-3, ''); chk('Particle swarm reaches the minimum', 0, particleSwarm(f, lo, hi, { n: 20, iters: 60 }).f, 1e-4, ''); }
  flag('Non-dominated filter', paretoFront([[1, 2], [2, 1], [2, 2]]).join(',') === '0,1', paretoFront([[1, 2], [2, 1], [2, 2]]).join(','), '0,1', '(2, 2) is dominated');
  { const fr = nsga2((x) => [x[0] ** 2, (x[0] - 2) ** 2], [-5], [5], { pop: 30, gens: 25, seed: 1 }); flag('Multi-objective search finds the Pareto set', fr.length > 5 && fr.every((s) => s.x[0] > -0.05 && s.x[0] < 2.05), fr.length, '> 5 points in [0, 2]', 'Schaffer problem'); }
  // --- cost models, reliability, units
  { const a = { costIndexEval: 800, costIndexBase: 800, locFactor: 1, flowLen: 1000, riserLen: 0, id: 0.254, wt: 0.0159, depth: 0, insT: 0, material: 'cs', nWells: 1, megRate: 0, chemRate: 0, pumpKW: 0, compKW: 0, heatKW: 0, steelPrice: 1800, coatPrice: 60, fabPerM: 120, insPrice: 5000, pipPremium: 650, dehCable: 450, riserFactor: 2.5, vesselRate: 350, layRate: 2.5, mobCost: 6, layFactor: 1, depthCoef: 0.2, wellCost: 70, learnRate: 0.9, langFactor: 4, costMethod: 'module', contingency: 0.15, owners: 0.08, craFactor: 4.5, strategy: 'none' }, sc = (vol) => capexEstimate({ ...a, slugVol: vol }).items.find((i) => i.item.startsWith('Slug catcher')).cost, ce = capexEstimate({ ...a, slugVol: 60 });
    chk('Six-tenths rule: doubling capacity', 1.515716566510398, sc(120) / sc(60), 1e-12, '2^0.6'); chk('Contingency and owner\'s costs on direct cost', 1.23, ce.total / ce.direct, 1e-12, '1 + 15 % + 8 %'); chk('Line-pipe steel tonnage', (7850 * Math.PI * (0.2858 ** 2 - 0.254 ** 2)) / 4, ce.steelT, 1e-9, 'ρ π (D_o² − D_i²)/4 × 1 km, in tonnes'); chk('Cost-index escalation', 1.05, capexEstimate({ ...a, slugVol: 60, costIndexEval: 840 }).total / ce.total, 1e-12, '840/800'); }
  chk('Learning curve: four units at 90 %', 3.556205986312342, learningSum(4, 0.9), 1e-12, '1 + 0.9 + 0.8462 + 0.81');
  { const n = 4000, h = 25 / n; let s = 0; for (let i = 0; i < n; i++) s += h * hazardRate((i + 0.5) * h, { h0: 0, beta: 3, remLife: 25, pEnd: 0.5 }); chk('Weibull hazard integrates to the stated end-of-life probability', Math.log(2), s, 1e-6, '∫h dt = −ln(1 − 0.5)'); }
  { const a = { life: 20, inspCost: 1, consequence: 100, h0: 0.001, beta: 3, remLife: 20, pEnd: 0.5, rate: 0 }; flag('Inspection lowers the expected failure cost', inspectionCost({ ...a, interval: 5, pod: 0.9 }).failure < inspectionCost({ ...a, interval: 1e9, pod: 0.9 }).failure, inspectionCost({ ...a, interval: 5, pod: 0.9 }).failure, `< ${inspectionCost({ ...a, interval: 1e9, pod: 0.9 }).failure.toFixed(2)}`); chk('Inspection count over the life', 3, inspectionCost({ ...a, interval: 5, pod: 0.9 }).inspection, 1e-12, 'at years 5, 10, 15'); }
  chk('Unit conversion: barrels per cubic metre', 6.289810770432105, BBL_PER_M3, 1e-8, '1 bbl = 0.158987294928 m³');
  chk('Unit conversion: one boe in GJ', 6.11932393, MMBTU_PER_BOE * GJ_PER_MMBTU, 1e-8, '5.8 MMBtu × 1.05505585 GJ/MMBtu');
  { const usd = 123.45, fx = 1500; chk('Currency conversion round trip', usd, (usd * fx) / fx, 1e-12, 'US$ → local → US$'); }
  chk('Equivalent annual value × annuity factor = NPV', 1000, annuityPV(equivalentAnnual(1000, 0.1, 7), 0.1, 7), 1e-9, '');
  return out;
}

export default {
  id: 'econ',
  num: 7,
  title: 'Economics, Techno-Economics & Decision Analysis',
  short: 'Economics',
  icon: '💲',
  tagline: 'Turns pressure drop, temperature, hydrate and wax exposure, slugging, corrosion and downtime into cash flow, risk and ranked decisions.',
  description: 'Builds the CAPEX and OPEX of the case line from parametric cost models, runs a fiscal cash flow (royalty and tax or production sharing) and reports NPV, IRR, payback, unit cost and break-evens. The flow kernel is solved for a set of insulation thicknesses, diameters and rates so that six flow-assurance strategies, the insulation thickness, the bore and the operating rate are optimised on lifecycle cost and NPV. Uncertainty is propagated by correlated Monte Carlo or Latin-hypercube sampling, and the decision is supported by decision trees, value of information, multi-criteria ranking, real options and mathematical programming.',
  guide: [
    'Run suites 1–6 first if you can: rates, line data, uptime, inhibitor demand, plug probability, failure probability and remaining life are then offered as linked values. The suite also works alone on the reference case.',
    'Set the fiscal frame (discount rate, inflation, royalty, tax or production sharing, depreciation) and the production profile (plateau, decline, water cut, uptime).',
    'Choose the flow-assurance strategy of the case. All six strategies are always compared on lifecycle cost, and the insulation thickness, bore and operating rate are optimised with the flow kernel.',
    'Review the cost basis on the setup tab: reference costs, exponents and installation factors are editable; the parametric line-cost relationship can be fitted to past projects on the calibration tab.',
    'Edit the distributions, correlations and scenarios, then read P10/P50/P90, probability of loss, CVaR, the tornado and the Sobol indices.',
    'Use the convergence tab to confirm that the sample size and the lattice steps are large enough for the decision.',
  ],
  equationsNote: 'Screening-level (class 4–5) cost models in 2023 US dollars with editable reference costs; a single-field, single-line project; one price multiplier drives oil and gas together. The thermal and hydraulic response between the kernel solutions is interpolated (arrival temperature as an exponential in U and in 1/rate, friction as a power of rate). Cooldown uses a lumped thermal mass. Downtime volume marked as deferred is recovered in the last production year. The production-sharing option has one cost-oil cap and one profit split (no R-factor or sliding scale). Hazard is a constant rate plus a Weibull wear-out term; inspection acts through a virtual-age reduction. Real options assume a lognormal project value. Bayesian optimisation, time-series and maximum-likelihood calibration of price and cost histories are not implemented.',
  implemented: IMPLEMENTED,
  referenceOnly: REFERENCE_ONLY,
  inputs: INPUTS,
  presets: PRESETS,
  pull,
  site: siteHook,
  run,
  mesh: [
    { name: 'Sample size / convergence', keys: ['nMC'], min: 200, note: 'The simulation is repeated with 1×, 2× and 4× the sample size. Sampling error falls as 1/√N for Monte Carlo and faster for Latin hypercube, so the observed order is statistical rather than a discretisation order.', metrics: [{ label: 'Mean NPV', unit: 'M$', get: (res) => res.outputs.npvMean / MM }, { label: 'NPV P10', unit: 'M$', get: (res) => res.outputs.npvP10 / MM }, { label: 'CVaR', unit: 'M$', get: (res) => res.outputs.cvar / MM }] },
    { name: 'Lattice steps', keys: ['nLattice'], min: 10, note: 'Binomial-lattice values converge to the continuous-time limit at first order in the time step, with the usual odd–even oscillation.', metrics: [{ label: 'Option to defer', unit: 'M$', get: (res) => res.outputs.optionValue / MM }, { label: 'Option to expand', unit: 'M$', get: (res) => res.outputs.optionExpand / MM }] },
  ],
  calibration: {
    note: 'Benchmark and back-cast: fit the cost-estimating relationship for an installed pipeline — cost = [coefficient × (D/10 in)^exponent × length + day-rate factor × vessel spread × lay days + mobilisation] × unit number^log₂(learning rate) — to historical project costs normalised to the basis year and location, then check it on projects that were not used in the fit. The sample is synthetic: generated from the model with different parameters and 4 % noise.',
    params: [{ key: 'cerCoef', label: 'Line cost at 10 in (M$/km)', lo: 0.1, hi: 3 }, { key: 'cerExp', label: 'Diameter exponent', lo: 0.5, hi: 2.5 }, { key: 'layFactor', label: 'Installation day-rate factor', lo: 0.4, hi: 3 }, { key: 'learnRate', label: 'Learning-curve rate (%)', lo: 70, hi: 100 }],
    columns: [{ key: 'calD', label: 'Diameter', unit: 'in' }, { key: 'calL', label: 'Length', unit: 'km' }, { key: 'calDepth', label: 'Water depth', unit: 'm' }, { key: 'calUnit', label: 'Project sequence number', unit: '' }, { key: 'cost', label: 'Installed cost', unit: 'M$' }],
    targets: [{ key: 'cost', label: 'Installed cost', unit: 'M$' }],
    model: calModel,
    sample: [{ calD: 8, calL: 12, calDepth: 300, calUnit: 1, cost: 13.15 }, { calD: 10, calL: 18, calDepth: 1350, calUnit: 1, cost: 17.96 }, { calD: 12, calL: 35, calDepth: 900, calUnit: 2, cost: 33.98 }, { calD: 16, calL: 60, calDepth: 150, calUnit: 1, cost: 75.58 }, { calD: 10, calL: 8, calDepth: 2000, calUnit: 3, cost: 10.77 }, { calD: 14, calL: 45, calDepth: 1100, calUnit: 2, cost: 44.71 }, { calD: 20, calL: 95, calDepth: 400, calUnit: 4, cost: 123.27 }, { calD: 6, calL: 6, calDepth: 600, calUnit: 1, cost: 8.87 }, { calD: 12, calL: 22, calDepth: 1800, calUnit: 5, cost: 22.13 }, { calD: 24, calL: 140, calDepth: 120, calUnit: 3, cost: 217.45 }, { calD: 18, calL: 70, calDepth: 750, calUnit: 6, cost: 89.94 }, { calD: 8, calL: 28, calDepth: 1500, calUnit: 2, cost: 20.9 }],
    validationSample: [{ calD: 10, calL: 25, calDepth: 1000, calUnit: 2, cost: 22.41 }, { calD: 14, calL: 30, calDepth: 500, calUnit: 1, cost: 36.98 }, { calD: 16, calL: 80, calDepth: 1300, calUnit: 4, cost: 80.88 }, { calD: 8, calL: 15, calDepth: 250, calUnit: 3, cost: 12.26 }, { calD: 22, calL: 110, calDepth: 200, calUnit: 2, cost: 170.52 }, { calD: 12, calL: 50, calDepth: 1600, calUnit: 7, cost: 40.77 }],
  },
  verify,
};
