// Suite 1 — Fluid, PVT & Phase Behaviour.
// Builds on the shared cubic-EOS kernel (core/thermo.js) and adds: alternative C7+ critical-property correlations and tuning,
// PH / PS / TV flashes, saturation points and a traced phase envelope with quality lines, the standard PVT-laboratory
// experiments (CCE, DLE, CVD, separator test, swelling, multiple-contact miscibility), derived properties (speed of sound,
// compressibility, black-oil correlations), residual-Helmholtz equations of state beside the cubic one (GERG-2008, PC-SAFT,
// cubic-plus-association, Lee–Kesler modified BWR), Lohrenz–Bray–Clark and Pedersen corresponding-states viscosity, aqueous-phase
// activity models (NRTL, UNIQUAC, Wilson, Debye–Hückel family, Pitzer), Henry's-law gas solubility, the van der Waals–Platteeuw
// hydrate model (Munck constants or Kihara cell potential; structures I, II and H) with the water-side Gibbs minimum, ideal-solution
// wax precipitation, Flory–Huggins asphaltene equilibrium and asphaltene screening. Every literature constant set is listed in
// PROVENANCE with the source it was checked against; the sourced reference data sets live in ../data/ref/pvt.js.
// Units at the interface: bara, °C, mol %; SI inside.
import { brent, clamp, linspace, logspace, interp1, isNum, solveLinear } from '../core/num.js';
import {
  R, P_STD, T_STD, VM_STD, MW_AIR, COMPONENTS, COMP_IDS, COMP_LABELS, INHIBITORS, EOS, DEFAULT_FLUID, makeFluid, eosPhase, rachfordRice, stability, flashPT,
  phaseProps, props, interfacialTension, saturationP, stdFlash, streams, aqueous, waterProps, waterContent, hydrateDepression, inhibitorFor, hydrateT0, hydrateScreening, pseudoProps, buildTable, lookup,
} from '../core/thermo.js';
import { psat as psatWater, density as rhoBrine } from '../core/props.js';
import { BASE } from '../data/basecase.js';
import { SOURCES, GERG, REF_ISO, REF_SAT, BENCH, BENCH_META, HYDRATE_DATA, HYDRATE_REFS, NACL_25C, MEOH_FREEZING, ACTIVITY_CALIBRATION, C1_WATER_VLE, C1_MEOH_VLE, C1_MEOH_WATER_VLE, C1_SOLUBILITY } from '../data/ref/pvt.js';

const KEL = 273.15, ATM = 1.01325, MW_W = 18.015, MW_NACL = 58.443, G0 = 9.80665;
const sumA = (a) => { let s = 0; for (let i = 0; i < a.length; i++) s += a[i]; return s; };
const num = (x, d) => (Number.isFinite(+x) ? +x : d);
const fin = (x, d = null) => (Number.isFinite(x) ? x : d);
const rd = (x, n = 4) => (Number.isFinite(x) ? +(+x).toPrecision(n) : null);
const cell = (x, n = 4) => (Number.isFinite(x) ? +(+x).toPrecision(n) : '—');

// ---- C7+ characterisation: alternative critical-property correlations ----------------------------------------------
const acentricLK = (Tb, Tc, PcBar, SG) => {
  const Tbr = Tb / Tc, Kw = (Tb * 1.8) ** (1 / 3) / SG, p = PcBar / ATM;
  return Tbr < 0.8
    ? (-Math.log(p) - 5.92714 + 6.09648 / Tbr + 1.28862 * Math.log(Tbr) - 0.169347 * Tbr ** 6) / (15.2518 - 15.6875 / Tbr - 13.4721 * Math.log(Tbr) + 0.43577 * Tbr ** 6)
    : -7.904 + 0.1352 * Kw - 0.007465 * Kw * Kw + 8.359 * Tbr + (1.408 - 0.01063 * Kw) / Tbr;
};
/** Riazi–Daubert (1980) critical temperature (K) and pressure (bar) of a petroleum fraction from Tb (K) and specific gravity. */
export function critRiaziDaubert(Tb, SG) {
  const t = Tb * 1.8;
  return { Tc: (24.2787 * t ** 0.58848 * SG ** 0.3596) / 1.8, Pc: 3.12281e9 * t ** -2.3125 * SG ** 2.3201 * 0.0689476 };
}
/** Twu (1984) critical temperature (K) and pressure (bar): n-alkane reference plus a specific-gravity perturbation. */
export function critTwu(Tb, SG) {
  const t = Tb * 1.8, Tc0 = t / (0.533272 + 0.191017e-3 * t + 0.779681e-7 * t * t - 0.284376e-10 * t ** 3 + 0.959468e28 / t ** 13), a = 1 - t / Tc0;
  const Pc0 = (3.83354 + 1.19629 * Math.sqrt(a) + 34.8888 * a + 36.1952 * a * a + 104.193 * a ** 4) ** 2;
  const Vc0 = (1 - (0.419869 - 0.505839 * a - 1.56436 * a ** 3 - 9481.7 * a ** 14)) ** -8, SG0 = 0.843593 - 0.128624 * a - 3.36159 * a ** 3 - 13749.5 * a ** 12;
  const sq = Math.sqrt(t), dT = Math.exp(5 * (SG0 - SG)) - 1, fT = dT * (-0.362456 / sq + (0.0398285 - 0.948125 / sq) * dT), Tc = Tc0 * ((1 + 2 * fT) / (1 - 2 * fT)) ** 2;
  const dV = Math.exp(4 * (SG0 * SG0 - SG * SG)) - 1, fV = dV * (0.46659 / sq + (-0.182421 + 3.01721 / sq) * dV), Vc = Vc0 * ((1 + 2 * fV) / (1 - 2 * fV)) ** 2;
  const dP = Math.exp(0.5 * (SG0 - SG)) - 1, fP = dP * (2.53262 - 46.1955 / sq - 0.00127885 * t + (-11.4277 + 252.14 / sq + 0.00230535 * t) * dP);
  return { Tc: Tc / 1.8, Pc: Pc0 * (Tc / Tc0) * (Vc0 / Vc) * ((1 + 2 * fP) / (1 - 2 * fP)) ** 2 * 0.0689476 };
}
/** Critical properties and acentric factor of a pseudo-component by 'KL' (Kesler–Lee), 'RD' (Riazi–Daubert) or 'Twu'. */
export function critProps(M, SG, method = 'KL') {
  const kl = pseudoProps(M, SG);
  if (method === 'KL') return { Tb: kl.Tb, Tc: kl.Tc, Pc: kl.Pc, w: kl.w };
  const c = method === 'Twu' ? critTwu(kl.Tb, SG) : critRiaziDaubert(kl.Tb, SG), Tc = clamp(c.Tc, 500, 1150), Pc = clamp(c.Pc, 5, 40);
  return { Tb: kl.Tb, Tc, Pc, w: clamp(acentricLK(kl.Tb, Tc, Pc, SG), 0.2, 1.7) };
}
/**
 * Tuned characterisation of a case fluid. `tune`: { eos, nPseudo, method, tcMult, pcMult, wMult, kijScale, vcMult, shift }.
 * The alternative correlations enter as plus-fraction-average ratios to Kesler–Lee so that the same model can be handed to
 * the kernel (`makeFluid` / `buildTable`) through its multipliers. Returns { opts, fluid, ratios, pseudo[] }.
 */
export function characterise(spec = DEFAULT_FLUID, tune = {}) {
  const t = { method: 'KL', tcMult: 1, pcMult: 1, wMult: 1, kijScale: 1, vcMult: 1, shift: true, ...tune };
  const eos = EOS[t.eos] ? t.eos : EOS[spec.eos] ? spec.eos : 'PR', nPseudo = clamp(Math.round(num(t.nPseudo, num(spec.nPseudo, 3))), 1, 3);
  const base = makeFluid(spec, { eos, nPseudo }), ps = base.comps.filter((c) => c.pseudo);
  let rT = 1, rP = 1, rW = 1;
  if (t.method !== 'KL' && ps.length) {
    let zt = 0, a = 0, b = 0, c = 0;
    for (const p of ps) { const q = critProps(p.MW, p.SG, t.method); zt += p.z; a += (p.z * q.Tc) / p.Tc; b += (p.z * q.Pc) / p.Pc; c += (p.z * q.w) / p.w; }
    rT = a / zt; rP = b / zt; rW = c / zt;
  }
  const opts = { eos, nPseudo, kijScale: clamp(num(t.kijScale, 1), 0, 3), vcMult: clamp(num(t.vcMult, 1), 0.5, 2), tcMult: clamp(num(t.tcMult, 1), 0.7, 1.3) * rT, pcMult: clamp(num(t.pcMult, 1), 0.6, 1.5) * rP, wMult: clamp(num(t.wMult, 1), 0.5, 1.6) * rW, shift: t.shift !== false };
  const fluid = makeFluid(spec, opts);
  return { opts, fluid, ratios: { Tc: rT, Pc: rP, w: rW }, pseudo: ps.map((p) => ({ id: p.id, z: p.z, MW: p.MW, SG: p.SG, Tb: p.Tb, KL: critProps(p.MW, p.SG, 'KL'), RD: critProps(p.MW, p.SG, 'RD'), Twu: critProps(p.MW, p.SG, 'Twu') })) };
}

function lnGamma(x) { // Lanczos
  const g = [76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5];
  let y = x, t = x + 5.5, s = 1.000000000190015; t -= (x + 0.5) * Math.log(t);
  for (let j = 0; j < 6; j++) s += g[j] / ++y;
  return -t + Math.log((2.5066282746310005 * s) / x);
}
function gammaP(a, x) { // regularised lower incomplete gamma function by its series
  if (x <= 0) return 0;
  let term = 1 / a, s = term;
  for (let k = 1; k < 600; k++) { term *= x / (a + k); s += term; if (term < s * 1e-14) break; }
  return clamp(s * Math.exp(-x + a * Math.log(x) - lnGamma(a)), 0, 1);
}
/**
 * Single-carbon-number distribution of the plus fraction: 'exp' (Pedersen: ln z linear in carbon number) or 'gamma' (Whitson,
 * shape alpha). Returns [{ n, z (mole fraction of the whole fluid), M, SG }] for C7 … C(nMax), honouring the C7+ amount, molar
 * mass (exponential model exactly, gamma model through its mean) and specific gravity (Søreide factor).
 */
export function scnDistribution(z7, M7, SG7, { model = 'exp', alpha = 1, nMax = 80 } = {}) {
  if (!(z7 > 0)) return [];
  M7 = clamp(M7, 100, 600);
  const ns = []; for (let n = 7; n <= nMax; n++) ns.push(n);
  const Ms = ns.map((n) => 14 * n - 4);
  let w;
  if (model === 'gamma') {
    const a = clamp(alpha, 0.5, 3), eta = 89, beta = (M7 - eta) / a, cdf = (M) => gammaP(a, Math.max(0, M - eta) / beta);
    w = ns.map((n, i) => (i === ns.length - 1 ? 1 : cdf(14 * n + 3)) - (i === 0 ? 0 : cdf(14 * n - 11)));
  } else {
    const mean = (B) => { let a = 0, b = 0; ns.forEach((n, i) => { const e = Math.exp(B * (n - 7)); a += e * Ms[i]; b += e; }); return a / b; };
    let B = -0.1; try { B = brent((q) => mean(q) - M7, -3, -1e-4, 1e-12); } catch { B = mean(-1e-4) < M7 ? -1e-4 : -3; }
    w = ns.map((n) => Math.exp(B * (n - 7)));
  }
  const tot = sumA(w), sgOf = (Cf) => Ms.map((M) => clamp(0.2855 + Cf * (M - 66) ** 0.13, 0.6, 1.25));
  const mixSG = (Cf) => { const s = sgOf(Cf); let a = 0, b = 0; w.forEach((wi, i) => { a += wi * Ms[i]; b += (wi * Ms[i]) / s[i]; }); return a / b; };
  let Cf = 0.29; try { Cf = brent((c) => mixSG(c) - SG7, 0.18, 0.5, 1e-10); } catch { /* keep the Søreide default */ }
  const SGs = sgOf(Cf);
  return ns.map((n, i) => ({ n, z: (z7 * w[i]) / tot, M: Ms[i], SG: SGs[i] }));
}

// ---- EOS thermodynamics with analytic temperature derivative -----------------------------------------------------
function aDeriv(f, x, T) { // mixture attraction parameter a and da/dT
  const n = f.n, s = new Array(n), ds = new Array(n), id = f.eosId;
  for (let i = 0; i < n; i++) {
    const c = f.comps[i], ra = Math.sqrt(c.ac);
    if (id === 'vdW') { s[i] = ra; ds[i] = 0; }
    else if (id === 'RK') { s[i] = ra * (T / c.Tc) ** -0.25; ds[i] = (-0.25 * s[i]) / T; }
    else { s[i] = ra * (1 + c.m * (1 - Math.sqrt(T / c.Tc))); ds[i] = (-ra * c.m) / (2 * Math.sqrt(T * c.Tc)); }
  }
  let a = 0, da = 0;
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) { const k = x[i] * x[j] * (1 - f.kij[i][j]); a += k * s[i] * s[j]; da += k * (s[i] * ds[j] + s[j] * ds[i]); }
  return { a, da };
}
const cpIg = (f, x, T) => { let s = 0; for (let i = 0; i < f.n; i++) { const c = f.comps[i].cp; s += x[i] * (c[0] + c[1] * T + c[2] * T * T + c[3] * T ** 3); } return s; };
const hIg = (f, x, T) => { let s = 0; for (let i = 0; i < f.n; i++) { const c = f.comps[i].cp; s += x[i] * (c[0] * (T - KEL) + (c[1] / 2) * (T * T - KEL * KEL) + (c[2] / 3) * (T ** 3 - KEL ** 3) + (c[3] / 4) * (T ** 4 - KEL ** 4)); } return s; };
const sIg = (f, x, T, Pbar) => { // ideal-gas entropy relative to 0 °C, 1 bar, pure components
  let s = 0;
  for (let i = 0; i < f.n; i++) { const c = f.comps[i].cp; s += x[i] * (c[0] * Math.log(T / KEL) + c[1] * (T - KEL) + (c[2] / 2) * (T * T - KEL * KEL) + (c[3] / 3) * (T ** 3 - KEL ** 3)); if (x[i] > 0) s -= R * x[i] * Math.log(x[i]); }
  return s - R * Math.log(Pbar);
};
/**
 * Molar enthalpy (J/mol, ideal gas at 0 °C = 0), entropy (J/mol/K, ideal gas at 0 °C and 1 bar = 0) and residual parts of
 * one phase of composition x at P (bara), T (K) with the analytic da/dT. kind: 'liquid' | 'vapour' | 'stable'.
 */
export function phaseHS(f, x, Pbar, TK, kind = 'stable') {
  const ph = eosPhase(f, x, Pbar, TK, kind), { a, da } = aDeriv(f, x, TK);
  const hRes = R * TK * (ph.Z - 1) + ((TK * da - a) * ph.L * Pbar * 1e5) / (R * TK);
  let gRes = 0; for (let i = 0; i < f.n; i++) gRes += x[i] * ph.lnphi[i];
  gRes *= R * TK;
  const sRes = (hRes - gRes) / TK;
  return { h: hIg(f, x, TK) + hRes, s: sIg(f, x, TK, Pbar) + sRes, hRes, sRes, gRes, Z: ph.Z, lnphi: ph.lnphi, cpIg: cpIg(f, x, TK) };
}
/** A copy of an EOS fluid with another overall composition (mole fractions, normalised here). */
export function withZ(f, z) {
  const t = sumA(z), zz = z.map((v) => Math.max(v, 0) / t);
  let MW = 0; zz.forEach((v, i) => (MW += v * f.comps[i].MW));
  return { ...f, z: zz, comps: f.comps.map((c, i) => ({ ...c, z: zz[i] })), MW };
}
/**
 * Equilibrium state of the feed at P (bara), T (°C): flash plus total molar enthalpy H (J/mol), entropy S (J/mol/K) and
 * volume V (m³/mol, volume-translated). Returns { phase, beta, x, y, K, H, S, V, vL, vV, liq, vap }.
 */
export function stateHSV(f, Pbar, Tc) {
  const TK = Tc + KEL, fl = flashPT(f, Pbar, Tc);
  const vol = (x, kind) => phaseProps(f, x, Pbar, Tc, kind, { thermal: false }).vm;
  if (fl.phase !== 'two') {
    const kind = fl.phase === 'gas' ? 'vapour' : 'liquid', q = phaseHS(f, f.z, Pbar, TK, 'stable'), v = vol(f.z, kind);
    return { ...fl, H: q.h, S: q.s, V: v, vL: v, vV: v };
  }
  const l = phaseHS(f, fl.x, Pbar, TK, 'liquid'), g = phaseHS(f, fl.y, Pbar, TK, 'vapour'), vL = vol(fl.x, 'liquid'), vV = vol(fl.y, 'vapour'), b = fl.beta;
  return { ...fl, H: b * g.h + (1 - b) * l.h, S: b * g.s + (1 - b) * l.s, V: b * vV + (1 - b) * vL, vL, vV };
}
const solveMono = (g, lo, hi, tol) => { const a = g(lo), b = g(hi); if (!(a * b <= 0)) return null; return brent(g, lo, hi, tol, 80); };
/** Isenthalpic flash: temperature (°C) at which the feed has total molar enthalpy H (J/mol) at P (bara). Returns the state plus T, or null. */
export function flashPH(f, Pbar, H, { Tmin = -150, Tmax = 700 } = {}) {
  const T = solveMono((t) => stateHSV(f, Pbar, t).H - H, Tmin, Tmax, 1e-7);
  return T === null ? null : { T, ...stateHSV(f, Pbar, T) };
}
/** Isentropic flash: temperature (°C) at which the feed has total molar entropy S (J/mol/K) at P (bara). */
export function flashPS(f, Pbar, S, { Tmin = -150, Tmax = 700 } = {}) {
  const T = solveMono((t) => stateHSV(f, Pbar, t).S - S, Tmin, Tmax, 1e-7);
  return T === null ? null : { T, ...stateHSV(f, Pbar, T) };
}
/** Isochoric–isothermal (TV) flash: pressure (bara) at which the feed occupies molar volume V (m³/mol) at T (°C). */
export function flashTV(f, Tc, V, { Pmin = 0.05, Pmax = 3000 } = {}) {
  const lp = solveMono((q) => Math.log(stateHSV(f, Math.exp(q), Tc).V / V), Math.log(Pmin), Math.log(Pmax), 1e-10);
  return lp === null ? null : { P: Math.exp(lp), ...stateHSV(f, Math.exp(lp), Tc) };
}

// ---- saturation points and phase envelope ---------------------------------------------------------------------------
const isTwo = (f, P, Tc) => flashPT(f, P, Tc).phase === 'two';
// Equations of a line of constant phase fraction: the feed splits into phase a (fraction 1 − beta) and phase b (fraction beta);
// unknowns X = [ln K_1 … ln K_n, ln T, ln P] with K = b/a. beta = 0 is the saturation line (b incipient: dew and bubble branches).
function lineSystem(f, beta) {
  const n = f.n, z = f.z;
  return (X) => {
    const T = Math.exp(X[n]), P = Math.exp(X[n + 1]), a = new Array(n), b = new Array(n); let sa = 0, sb = 0;
    for (let i = 0; i < n; i++) { const K = Math.exp(X[i]); a[i] = z[i] / (1 - beta + beta * K); b[i] = K * a[i]; sa += a[i]; sb += b[i]; }
    for (let i = 0; i < n; i++) { a[i] /= sa; b[i] /= sb; }
    const pa = eosPhase(f, a, P, T), pb = eosPhase(f, b, P, T), F = new Array(n + 2);
    for (let i = 0; i < n; i++) F[i] = X[i] + pb.lnphi[i] - pa.lnphi[i];
    F[n] = sb - sa;
    return F;
  };
}
// Newton iteration on the line equations with the specification X[s] = S. Returns { ok, X, J, it }.
function lineNewton(sys, X0, s, S, { maxIter = 16, tol = 1e-9 } = {}) {
  const N = X0.length, n = N - 2; let X = X0.slice(), J = null;
  const jac = (X, F) => { const M = Array.from({ length: N }, () => new Array(N).fill(0)); for (let j = 0; j < N; j++) { const h = 1e-6, Xj = X.slice(); Xj[j] += h; const Fj = sys(Xj); for (let i = 0; i <= n; i++) M[i][j] = (Fj[i] - F[i]) / h; } M[N - 1][s] = 1; return M; };
  for (let it = 0; it < maxIter; it++) {
    const F = sys(X); F[N - 1] = X[s] - S;
    let err = 0; for (let i = 0; i < N; i++) { if (!Number.isFinite(F[i])) return { ok: false }; err = Math.max(err, Math.abs(F[i])); }
    J = jac(X, F);
    if (err < tol) { let mk = 0; for (let i = 0; i < n; i++) mk = Math.max(mk, Math.abs(X[i])); return { ok: mk > 1e-5, X, J, it }; }
    let dX; try { dX = solveLinear(J, F.map((v) => -v)); } catch { return { ok: false }; }
    let mx = 0; for (const d of dX) mx = Math.max(mx, Math.abs(d));
    if (!Number.isFinite(mx)) return { ok: false };
    const lam = mx > 0.8 ? 0.8 / mx : 1;
    for (let i = 0; i < N; i++) X[i] += lam * dX[i];
  }
  return { ok: false };
}
// Starting point of a line at a low pressure on the dew side: phase b is the liquid of mole fraction beta.
function dewStart(f, beta, Pfloor, td = dewT(f, Pfloor, { Tmax: 900 })) {
  if (td === null) return null;
  let t = td - 0.5;
  if (beta > 0) {
    const g = (T) => { const fl = flashPT(f, Pfloor, T); return (fl.phase === 'two' ? 1 - fl.beta : fl.phase === 'gas' ? 0 : 1) - beta; };
    let lo = td - 0.5; for (let k = 0; k < 60 && g(lo) < 0; k++) lo -= 15;
    if (g(lo) < 0) return null;
    try { t = brent(g, lo, td - 0.5, 1e-6, 60); } catch { return null; }
  }
  // the flash can miss a vanishing liquid fraction just inside the dew point: step further in until it resolves two phases
  let fl = flashPT(f, Pfloor, t);
  for (const dt of [2, 5, 10, 20]) { if (fl.phase === 'two' || beta > 0) break; t = td - dt; fl = flashPT(f, Pfloor, t); }
  if (fl.phase !== 'two') return null;
  return [...fl.K.map((k) => -Math.log(k)), Math.log(t + KEL), Math.log(Pfloor)];
}
// Bubble point at T (°C) from Wilson K-values refined by Newton and confirmed by two flashes; null when that fails.
function bubbleWilson(f, Tc) {
  const TK = Tc + KEL, n = f.n, kw = f.comps.map((c) => c.Pc * Math.exp(5.373 * (1 + c.w) * (1 - c.Tc / TK))), P0 = f.z.reduce((s, z, i) => s + z * kw[i], 0);
  if (!(P0 > 0.01 && P0 < 5000)) return null;
  const lt = Math.log(TK), r = lineNewton(lineSystem(f, 0), [...kw.map((k) => Math.log(k / P0)), lt, Math.log(P0)], n, lt, { maxIter: 40 });
  if (!r.ok) return null;
  const P = Math.exp(r.X[n + 1]), lo = flashPT(f, P * 0.97, Tc);
  return P > 0.05 && P < 3000 && lo.phase === 'two' && lo.beta < 0.5 && !isTwo(f, P * 1.03, Tc) ? r.X : null;
}
// Starting point on the bubble side at a low temperature: phase b is the vapour of mole fraction q.
function bubbleStart(f, q, Tmin, Pfloor, Pmax) {
  for (let k = 0; k < 3; k++) {
    const t = Tmin + 25 * k; let X = bubbleWilson(f, t);
    if (!X) { const s = saturationPoint(f, t, { Pmax }); if (s.P === null || s.capped || !s.X || s.type !== 'bubble') continue; X = s.X; }
    const Ps = Math.exp(X[f.n + 1]);
    if (Ps > Pmax || Ps < Pfloor) continue;
    if (!(q > 0)) return X;
    const g = (lp) => { const fl = flashPT(f, Math.exp(lp), t); return (fl.phase === 'two' ? fl.beta : fl.phase === 'gas' ? 1 : 0) - q; };
    const lo = Math.log(Pfloor), hi = Math.log(Ps * 0.9995);
    if (!(g(lo) > 0) || !(g(hi) < 0)) return null;
    let lp; try { lp = brent(g, lo, hi, 1e-7, 60); } catch { return null; }
    const fl = flashPT(f, Math.exp(lp), t);
    return fl.phase === 'two' ? [...fl.K.map((v) => Math.log(v)), Math.log(t + KEL), lp] : null;
  }
  return null;
}
/**
 * Continuation trace of a line of constant phase fraction (Michelsen-type): Newton on [ln K, ln T, ln P] with the fastest-moving
 * variable as the specification. from: 'dew' starts at Pfloor on the dew side (phase b = liquid of fraction beta), 'bubble' at
 * Tmin on the bubble side (phase b = vapour of fraction beta); the trace stops at the critical point, at Tmin, Pfloor or Pmax.
 * Returns { pts: [{ T (°C), P }], critical: { T, P } | null }.
 */
export function traceLine(f, beta = 0, { from = 'dew', Tmin = -40, Pmax = 1200, Pfloor = 1, n = 40, maxPts = 0, tDew } = {}) {
  const N = f.n + 2, nc = f.n, sys = lineSystem(f, beta), out = { pts: [], critical: null };
  if (nc < 2) return out;
  const X0 = from === 'dew' ? dewStart(f, beta, Pfloor, tDew) : bubbleStart(f, beta, Tmin, Pfloor, Pmax);
  if (!X0) return out;
  const sc = clamp(40 / n, 0.1, 2); maxPts = maxPts || Math.round(6 * n + 60);
  let s = from === 'dew' ? N - 1 : N - 2, r = lineNewton(sys, X0, s, X0[s], { maxIter: 30 });
  if (!r.ok) return out;
  let X = r.X, dS = (from === 'dew' ? 0.15 : 0.03) * sc, still = 0;
  out.pts.push({ T: Math.exp(X[nc]) - KEL, P: Math.exp(X[nc + 1]), X });
  while (out.pts.length < maxPts) {
    // sensitivities along the line and choice of the next specification variable
    const e = new Array(N).fill(0); e[N - 1] = 1; let v; try { v = solveLinear(r.J, e); } catch { break; }
    let s2 = 0; for (let i = 0; i < N; i++) if (Math.abs(v[i]) > Math.abs(v[s2])) s2 = i;
    if (!(Math.abs(v[s2]) > 0)) break;
    const dir = Math.sign(v[s2] * dS) || 1, t = v.map((q) => q / v[s2]);
    let iK = 0, dk = 0; for (let i = 0; i < nc; i++) { if (Math.abs(X[i]) > Math.abs(X[iK])) iK = i; dk = Math.max(dk, Math.abs(t[i])); }
    const mk = Math.abs(X[iK]);
    let step = Math.abs(v[s2] * dS) * (r.it <= 2 ? 1.6 : r.it <= 4 ? 1.1 : 0.6);
    step = Math.min(step, (0.035 * sc) / Math.max(Math.abs(t[nc]), 1e-12), (0.22 * sc) / Math.max(Math.abs(t[nc + 1]), 1e-12), ((mk < 1.5 ? 0.25 : 0.7) * sc) / Math.max(dk, 1e-12));
    // approaching the critical point (all ln K -> 0): close in geometrically, then extrapolate along the tangent and stop
    const toZero = -X[iK] / (t[iK] * dir); // step that would bring the largest ln K to zero
    if (toZero > 0 && step > 0.6 * toZero) {
      if (mk < 0.05) { out.critical = { T: Math.exp(X[nc] + t[nc] * dir * toZero) - KEL, P: Math.exp(X[nc + 1] + t[nc + 1] * dir * toZero) }; break; }
      step = 0.6 * toZero;
    }
    let ok = false, r2 = null;
    for (let tr = 0; tr < 8; tr++) {
      const d = dir * step;
      r2 = lineNewton(sys, X.map((x, i) => x + t[i] * d), s2, X[s2] + d);
      if (r2.ok && Math.abs(r2.X[nc] - X[nc]) < 0.12 && Math.abs(r2.X[nc + 1] - X[nc + 1]) < 0.8 && Math.sign(r2.X[iK]) === Math.sign(X[iK])) { ok = true; dS = d; break; }
      step *= 0.4;
    }
    if (!ok) break;
    const Xn = r2.X, T = Math.exp(Xn[nc]) - KEL, P = Math.exp(Xn[nc + 1]);
    if (!(P >= Pfloor * 0.999) || P > Pmax || T < Tmin - 1e-6 || T > 1200) break;
    still = Math.abs(Xn[nc] - X[nc]) < 2e-5 && Math.abs(Xn[nc + 1] - X[nc + 1]) < 2e-5 ? still + 1 : 0;
    if (still > 3) break;
    out.pts.push({ T, P, X: Xn });
    X = Xn; r = r2; s = s2;
  }
  return out;
}
/**
 * Saturation pressure (bara) at T (°C). A converged neighbouring point (`seed`: { X } from an earlier call) is refined by Newton;
 * otherwise the kernel's robust pressure scan brackets the upper boundary, which is then polished by Newton on the
 * incipient-phase equations. Returns { P, type: 'bubble' | 'dew', X } or { P: null }.
 */
export function saturationPoint(f, Tc, { seed = null, Pmax = 1200 } = {}) {
  if (f.n < 2) return { P: null, type: null };
  const n = f.n, sys = lineSystem(f, 0), lt = Math.log(Tc + KEL), iL = f.comps.reduce((m, c, i) => (c.MW < f.comps[m].MW ? i : m), 0);
  const finish = (r) => ({ P: Math.exp(r.X[n + 1]), type: r.X[iL] > 0 ? 'bubble' : 'dew', X: r.X });
  if (seed && seed.length === n + 2) {
    const X0 = seed.slice(); X0[n] = lt; const r = lineNewton(sys, X0, n, lt, { maxIter: 25 });
    if (r.ok && Math.exp(r.X[n + 1]) < Pmax * 1.5 && Math.abs(r.X[n + 1] - seed[n + 1]) < 1.5) return finish(r);
  }
  const s = saturationP(f, Tc, { Pmax });
  if (s.P === null || s.capped) return { P: s.P, type: s.type, capped: !!s.capped, X: null };
  for (const q of [0.999, 0.99, 0.96]) {
    const fl = flashPT(f, s.P * q, Tc);
    if (fl.phase !== 'two') continue;
    const bub = fl.beta < 0.5, r = lineNewton(sys, [...fl.K.map((k) => (bub ? Math.log(k) : -Math.log(k))), lt, Math.log(s.P)], n, lt, { maxIter: 25 });
    if (r.ok && Math.abs(Math.exp(r.X[n + 1]) / s.P - 1) < 0.05) return finish(r);
  }
  return { P: s.P, type: s.type, X: null };
}
/** Bubble-point pressure (bara) at T (°C), or null when the upper saturation point is a dew point or does not exist. */
export const bubbleP = (f, Tc, o) => { const s = saturationPoint(f, Tc, o); return s.P !== null && s.type === 'bubble' ? s.P : null; };
/** Upper dew-point pressure (bara) at T (°C), or null. */
export const dewP = (f, Tc, o) => { const s = saturationPoint(f, Tc, o); return s.P !== null && s.type === 'dew' ? s.P : null; };
/** Highest temperature (°C) at which the feed is two-phase at P (bara) — the dew-point temperature — or null (bisection on the flash). */
export function dewT(f, Pbar, { Tmin = -80, Tmax = 750, tol = 0.01 } = {}) {
  let a = null, b = null, top = Tmax;
  // start the downward search a little above the Wilson dew temperature (Σ z/K = 1) instead of at Tmax: far fewer single-phase flashes
  { const g = (t) => { let q = 0; for (let i = 0; i < f.n; i++) { const c = f.comps[i]; q += (f.z[i] * Pbar) / (c.Pc * Math.exp(5.373 * (1 + c.w) * (1 - c.Tc / (t + KEL)))); } return q - 1; };
    if (g(Tmin) > 0 && g(Tmax) < 0) { let t0 = Math.min(Tmax, Tmin + 25 * Math.ceil((brent(g, Tmin, Tmax, 1e-3, 60) + 60 - Tmin) / 25)); while (t0 < Tmax && isTwo(f, Pbar, t0)) t0 = Math.min(Tmax, t0 + 50); top = t0; } }
  for (let t = top; t >= Tmin; t -= 25) { if (isTwo(f, Pbar, t)) { a = t; break; } b = t; }
  if (a === null || b === null) return null;
  while (b - a > tol) { const m = 0.5 * (a + b); if (isTwo(f, Pbar, m)) a = m; else b = m; }
  return 0.5 * (a + b);
}
const parabolaMax = (x, y) => { // vertex of the parabola through three points (falls back to the middle point)
  const [x0, x1, x2] = x, [y0, y1, y2] = y, d = (x0 - x1) * (x0 - x2) * (x1 - x2);
  if (!(Math.abs(d) > 1e-300)) return { x: x1, y: y1 };
  const A = (x2 * (y1 - y0) + x1 * (y0 - y2) + x0 * (y2 - y1)) / d, B = (x2 * x2 * (y0 - y1) + x1 * x1 * (y2 - y0) + x0 * x0 * (y1 - y2)) / d;
  if (!(A < 0)) return { x: x1, y: y1 };
  const xv = clamp(-B / (2 * A), Math.min(x0, x2), Math.max(x0, x2)), C = y1 - A * x1 * x1 - B * x1;
  return { x: xv, y: A * xv * xv + B * xv + C };
};
/**
 * Phase envelope by continuation. The saturation line is traced from the low-pressure dew point (through the cricondentherm and
 * cricondenbar) and from the cold end of the bubble line, both up to the critical point where they meet; lines of constant vapour
 * fraction are traced the same way. Returns { T, P, type } (upper line from the cricondentherm, T ascending), { lowT, lowP } (lower dew
 * line from the cricondentherm downwards), full { T, P } (whole line, dew side first), critical { T, P } | null,
 * cricondenbar { T, P } | null, cricondentherm { T, P } | null, quality: [{ q (vapour mole fraction), T: [], P: [] }].
 */
export function traceEnvelope(f, { n = 40, Tmin = -40, Pmax = 1200, Pfloor = 1, qualities = [0.1, 0.5, 0.9] } = {}) {
  const out = { T: [], P: [], type: [], lowT: [], lowP: [], full: { T: [], P: [] }, dew: { T: [], P: [] }, bubble: { T: [], P: [] }, seeds: [], closed: false, critical: null, cricondenbar: null, cricondentherm: null, quality: [] };
  if (f.n < 2) return out;
  n = clamp(Math.round(n), 8, 400);
  const o = { Tmin, Pmax, Pfloor, n, tDew: dewT(f, Pfloor, { Tmax: 900 }) }, d = traceLine(f, 0, { ...o, from: 'dew' }), dEnd = d.pts.length ? d.pts[d.pts.length - 1] : null, b = d.critical || (dEnd && dEnd.T > Tmin + 5 && dEnd.P > Pfloor * 1.5 && dEnd.P < Pmax * 0.9) ? traceLine(f, 0, { ...o, from: 'bubble' }) : { pts: [], critical: null };
  let crit = d.critical || b.critical;
  if (d.critical && b.critical && Math.abs(d.critical.T - b.critical.T) < 5) crit = { T: 0.5 * (d.critical.T + b.critical.T), P: 0.5 * (d.critical.P + b.critical.P) };
  const p = [...d.pts.map((q) => ({ ...q, type: 'dew' })), ...(crit ? [{ ...crit, type: 'dew', X: null }] : []), ...b.pts.slice().reverse().map((q) => ({ ...q, type: 'bubble' }))];
  out.closed = !!(d.critical && b.critical);
  out.dew = { T: [...d.pts.map((q) => q.T), ...(crit ? [crit.T] : [])], P: [...d.pts.map((q) => q.P), ...(crit ? [crit.P] : [])] };
  out.bubble = b.pts.length ? { T: [...(crit ? [crit.T] : []), ...b.pts.slice().reverse().map((q) => q.T)], P: [...(crit ? [crit.P] : []), ...b.pts.slice().reverse().map((q) => q.P)] } : { T: [], P: [] };
  if (p.length < 3) return out;
  out.full = { T: p.map((q) => q.T), P: p.map((q) => q.P) }; out.critical = crit;
  let iT = 0, iP = 0; p.forEach((q, i) => { if (q.T > p[iT].T) iT = i; if (q.P > p[iP].P) iP = i; });
  out.cricondentherm = { T: p[iT].T, P: p[iT].P }; out.cricondenbar = { T: p[iP].T, P: p[iP].P };
  if (iT > 0 && iT < p.length - 1) { const v = parabolaMax([Math.log(p[iT - 1].P), Math.log(p[iT].P), Math.log(p[iT + 1].P)], [p[iT - 1].T, p[iT].T, p[iT + 1].T]); out.cricondentherm = { T: v.y, P: Math.exp(v.x) }; }
  if (iP > 0 && iP < p.length - 1) { const v = parabolaMax([p[iP - 1].T, p[iP].T, p[iP + 1].T], [p[iP - 1].P, p[iP].P, p[iP + 1].P]); out.cricondenbar = { T: v.x, P: Math.max(v.y, p[iP].P) }; }
  for (let i = iT; i >= 0; i--) { out.lowT.push(p[i].T); out.lowP.push(p[i].P); }
  let last = Infinity; const up = [];
  for (let i = iT; i < p.length; i++) if (p[i].T < last - 1e-6) { up.push(p[i]); last = p[i].T; }
  up.reverse(); out.T = up.map((q) => q.T); out.P = up.map((q) => q.P); out.type = up.map((q) => q.type); out.seeds = up.filter((q) => q.X).map((q) => ({ T: q.T, X: q.X }));
  // constant-vapour-fraction lines: started on the dew side (phase b = liquid, fraction 1 − q) or, failing that, on the bubble side
  for (const q of qualities) {
    if (!(q > 0 && q < 1)) continue;
    const nn = Math.max(12, Math.round(n * 0.6)); let l = traceLine(f, 1 - q, { ...o, n: nn, from: 'dew' });
    if (l.pts.length < 2) l = traceLine(f, q, { ...o, n: nn, from: 'bubble' });
    const pts = l.critical ? [...l.pts, l.critical] : l.pts;
    if (pts.length > 1) out.quality.push({ q, T: pts.map((x) => x.T), P: pts.map((x) => x.P) });
  }
  return out;
}

// ---- derived single-phase properties ---------------------------------------------------------------------------------
/**
 * Derivative properties of one phase at P (bara), T (°C): isothermal compressibility kT (1/bar), expansivity (1/K), cv and cp
 * (J/kg/K), speed of sound (m/s) and Joule–Thomson coefficient (K/bar), from numerical derivatives of the translated volume.
 */
export function derivedProps(f, x, Pbar, Tc, kind) {
  const p = phaseProps(f, x, Pbar, Tc, kind), dP = Math.max(1e-3, Pbar * 2e-4), dT = 0.05, v = (P, T) => phaseProps(f, x, P, T, kind, { thermal: false }).vm;
  const dvdP = (v(Pbar + dP, Tc) - v(Pbar - dP, Tc)) / (2 * dP * 1e5), dvdT = (v(Pbar, Tc + dT) - v(Pbar, Tc - dT)) / (2 * dT), TK = Tc + KEL, M = p.MW * 1e-3;
  const cpm = p.cp * M, cvm = Math.max(cpm + (TK * dvdT * dvdT) / dvdP, 0.3 * cpm), gamma = cpm / cvm;
  return { ...p, kT: (-dvdP / p.vm) * 1e5, beta: dvdT / p.vm, cv: cvm / M, gamma, sound: p.vm * Math.sqrt(Math.max(-gamma / (M * dvdP), 0)), jtBar: p.jt * 1e5 };
}
const LK0 = [0.1181193, 0.265728, 0.15479, 0.030323, 0.0236744, 0.0186984, 0, 0.042724, 0.155488e-4, 0.623689e-4, 0.65392, 0.060167];
const LKR = [0.2026579, 0.331511, 0.027655, 0.203488, 0.0313385, 0.0503618, 0.016901, 0.041577, 0.48736e-4, 0.0740336e-4, 1.226, 0.03754];
function lkVr(c, Tr, Pr, liquid = false) { // reduced volume Vr = Pc·V/(R·Tc) of the Lee–Kesler modified BWR equation (vapour-like root unless liquid)
  const B = c[0] - c[1] / Tr - c[2] / Tr ** 2 - c[3] / Tr ** 3, C = c[4] - c[5] / Tr + c[6] / Tr ** 3, D = c[8] + c[9] / Tr;
  const Z = (V) => 1 + B / V + C / V ** 2 + D / V ** 5 + (c[7] / (Tr ** 3 * V * V)) * (c[10] + c[11] / (V * V)) * Math.exp(-c[11] / (V * V));
  const g = (lv) => { const V = Math.exp(lv); return Z(V) - (Pr * V) / Tr; };
  if (liquid) { let lo = Math.log(0.04), glo = g(lo); for (let k = 0; k < 400; k++) { const hi = lo + 0.02, ghi = g(hi); if (glo * ghi <= 0 && glo > 0) return Math.exp(brent(g, lo, hi, 1e-13)); lo = hi; glo = ghi; } return null; }
  // scan from the ideal-gas side towards small volumes for the first sign change
  let hi = Math.log((Tr / Pr) * 3), ghi = g(hi);
  for (let k = 0; k < 400; k++) { const lo = hi - 0.02, glo = g(lo); if (glo * ghi <= 0) return Math.exp(brent(g, lo, hi, 1e-13)); hi = lo; ghi = glo; }
  return null;
}
/** Lee–Kesler (1975) corresponding-states compressibility factor (vapour-like root) at reduced T, P and acentric factor w. */
export function leeKeslerZ(Tr, Pr, w = 0) { const r = leeKesler(Tr, Pr, w); return r ? r.Z : null; }
const kay = (f, x) => { let Tc = 0, Pc = 0, w = 0, M = 0, Vc = 0, Zc = 0; for (let i = 0; i < f.n; i++) { const c = f.comps[i]; Tc += x[i] * c.Tc; Pc += x[i] * c.Pc; w += x[i] * c.w; M += x[i] * c.MW; Vc += x[i] * c.Vc; Zc += x[i] * (0.2905 - 0.085 * c.w); } return { Tc, Pc, w, M, Vc, Zc }; };
/** Lee–Gonzalez–Eakin gas viscosity (Pa·s) from molar mass (g/mol), density (kg/m³) and T (K). */
export function viscosityLGE(M, rho, TK) {
  const T = TK * 1.8, K = ((9.4 + 0.02 * M) * T ** 1.5) / (209 + 19 * M + T), X = 3.5 + 986 / T + 0.01 * M, Y = 2.4 - 0.2 * X;
  return 1e-7 * K * Math.exp(X * (rho / 1000) ** Y);
}
/** Stiel–Thodos corresponding-states thermal conductivity (W/m/K) of a dense fluid from its dilute-gas value and reduced density. */
export function conductivityStielThodos(k0, rhoR, { Tc, Pc, M, Zc }) {
  const G = 210 * ((Tc * M ** 3) / Pc ** 4) ** (1 / 6), z5 = Zc ** 5, r = clamp(rhoR, 0, 2.8);
  const ex = r < 0.5 ? 1.22e-2 * (Math.exp(0.535 * r) - 1) : r < 2 ? 1.14e-2 * (Math.exp(0.67 * r) - 1.069) : 2.6e-3 * (Math.exp(1.155 * r) + 2.016);
  return k0 + ex / (G * z5);
}

// ---- black-oil correlations (field units inside: psia, °F, scf/STB) ---------------------------------------------------
const PSI = 14.5038, SCF_STB = 5.61458, degF = (c) => c * 1.8 + 32;
/**
 * Black-oil correlations as an independent cross-check of the EOS. Inputs: Rs (Sm³/Sm³), gas gravity, API, T (°C), P (bara).
 * Returns bubble points (bara), Bo (m³/Sm³) and viscosities (Pa·s) by Standing, Vasquez–Beggs, Glasø and Beggs–Robinson.
 */
export function blackOil({ rs, gasSG, api, T, P = null }) {
  const Rs = Math.max(rs * SCF_STB, 1), g = gasSG, F = degF(T), so = 141.5 / (131.5 + api);
  const pbStanding = (18.2 * ((Rs / g) ** 0.83 * 10 ** (0.00091 * F - 0.0125 * api) - 1.4)) / PSI;
  const vb = api <= 30 ? [0.0362, 1.0937, 25.724, 4.677e-4, 1.751e-5, -1.811e-8] : [0.0178, 1.187, 23.931, 4.67e-4, 1.1e-5, 1.337e-9];
  const pbVB = (Rs / (vb[0] * g * Math.exp((vb[2] * api) / (F + 460)))) ** (1 / vb[1]) / PSI;
  const lg = Math.log10((Rs / g) ** 0.816 * F ** 0.172 / api ** 0.989), pbGlaso = 10 ** (1.7669 + 1.7447 * lg - 0.30218 * lg * lg) / PSI;
  const boStanding = 0.9759 + 0.00012 * (Rs * Math.sqrt(g / so) + 1.25 * F) ** 1.2;
  const boVB = 1 + vb[3] * Rs + (F - 60) * (api / g) * (vb[4] + vb[5] * Rs);
  const lb = Math.log10(Rs * (g / so) ** 0.526 + 0.968 * F), boGlaso = 1 + 10 ** (-6.58511 + 2.91329 * lb - 0.27683 * lb * lb);
  const muDead = 10 ** (10 ** (3.0324 - 0.02023 * api) * F ** -1.163) - 1, muSat = 10.715 * (Rs + 100) ** -0.515 * muDead ** (5.44 * (Rs + 150) ** -0.338);
  let muP = muSat;
  if (P !== null && P > pbStanding) { const pp = P * PSI; muP = muSat * (P / pbStanding) ** (2.6 * pp ** 1.187 * Math.exp(-11.513 - 8.98e-5 * pp)); }
  return { pbStanding, pbVB, pbGlaso, boStanding, boVB, boGlaso, muDead: muDead * 1e-3, muSat: muSat * 1e-3, muP: muP * 1e-3 };
}
/** Standing solution gas–oil ratio (Sm³/Sm³) at P (bara) below the bubble point. */
export const rsStanding = (P, gasSG, api, T) => (gasSG * (((P * PSI) / 18.2 + 1.4) * 10 ** (0.0125 * api - 0.00091 * degF(T))) ** (1 / 0.83)) / SCF_STB;

// ---- laboratory experiments ---------------------------------------------------------------------------------------------
const pv = (f, P, Tc) => props(f, P, Tc, { thermal: false });
const molarV = (s) => (s.phase === 'two' ? s.beta * s.gas.vm + (1 - s.beta) * s.oil.vm : s.gas.vm);
/**
 * Constant composition (mass) expansion at T (°C). Returns rows [{ P, vRel, liqPct (% of the saturation volume), Y, Z, rho, co (1/bar) }].
 */
export function simulateCCE(f, Tc, Psat, pressures) {
  const vs = molarV(pv(f, Psat * 1.00002, Tc)), TK = Tc + KEL;
  return pressures.map((P) => {
    const s = pv(f, P, Tc), V = molarV(s), single = s.phase !== 'two';
    let co = null;
    if (single && P > Psat) { const d = Math.max(0.05, 1e-3 * P); co = -(molarV(pv(f, P + d, Tc)) - molarV(pv(f, Math.max(P - d, Psat * 1.00002), Tc))) / ((P + d - Math.max(P - d, Psat * 1.00002)) * V); }
    return { P, vRel: V / vs, liqPct: single ? (s.phase === 'oil' ? 100 * (V / vs) : 0) : (100 * (1 - s.beta) * s.oil.vm) / vs, Y: P < Psat * 0.999 && V > vs ? (Psat - P) / (P * (V / vs - 1)) : null, Z: (P * 1e5 * V) / (R * TK), rho: (f.MW * 1e-3) / V, co, beta: s.beta };
  });
}
/**
 * Differential liberation at T (°C) from the bubble point down to atmospheric pressure. Returns { rows: [{ P, rs, bo, bt, rhoO, gasSG, Z, bg }], residual: { rho, api } }
 * with volumes referred to the residual oil at standard conditions.
 */
export function simulateDLE(f, Tc, Psat, pressures) {
  const TK = Tc + KEL, steps = [];
  let z = f.z.slice(), nL = 1;
  const first = pv(f, Psat * 1.00002, Tc); steps.push({ P: Psat, vL: nL * molarV(first), nG: 0, rhoO: first.oil.rho, gasSG: null, Z: null });
  for (const P of pressures) {
    const fz = withZ(f, z), s = pv(fz, P, Tc);
    if (s.phase === 'two') { const nG = nL * s.beta; nL *= 1 - s.beta; z = s.x.slice(); steps.push({ P, vL: nL * s.oil.vm, nG, rhoO: s.oil.rho, gasSG: s.gas.MW / MW_AIR, Z: s.gas.Z }); }
    else if (s.phase === 'oil') steps.push({ P, vL: nL * s.oil.vm, nG: 0, rhoO: s.oil.rho, gasSG: null, Z: null });
    else { steps.push({ P, vL: 0, nG: nL, rhoO: null, gasSG: s.gas.MW / MW_AIR, Z: s.gas.Z }); nL = 0; break; }
  }
  if (!(nL > 0)) return { rows: [], residual: null };
  const res = phaseProps(withZ(f, z), z, P_STD, T_STD, 'liquid', { thermal: false }), vRes = nL * res.vm, sg = res.rho / 999.016;
  let after = 0; const rows = [];
  for (let i = steps.length - 1; i >= 0; i--) { // gas still in solution at step i = everything liberated at lower pressures
    const st = steps[i], bg = st.Z ? (st.Z * TK * P_STD) / (st.P * (T_STD + KEL)) : null, rs = (after * VM_STD) / vRes;
    rows.unshift({ P: st.P, rs, bo: st.vL / vRes, rhoO: st.rhoO, gasSG: st.gasSG, Z: st.Z, bg });
    after += st.nG;
  }
  const rsb = rows[0].rs; for (const r of rows) r.bt = r.bg ? r.bo + (rsb - r.rs) * r.bg : r.bo;
  return { rows, residual: { rho: res.rho, api: 141.5 / sg - 131.5 } };
}
/**
 * Constant volume depletion at T (°C) from the saturation pressure. Returns rows [{ P, liqPct (% of cell volume), produced (cumulative mol %), Zgas, Z2 (two-phase), heavy (C7+ mol % of produced gas) }].
 */
export function simulateCVD(f, Tc, Psat, pressures) {
  const TK = Tc + KEL, Vcell = molarV(pv(f, Psat * 1.00002, Tc)), rows = [{ P: Psat, liqPct: 0, produced: 0, Zgas: (Psat * 1e5 * Vcell) / (R * TK), Z2: (Psat * 1e5 * Vcell) / (R * TK), heavy: null }];
  let z = f.z.slice(), n = 1;
  const heavyOf = (y) => 100 * f.comps.reduce((s, c, i) => s + (c.pseudo ? y[i] : 0), 0);
  for (const P of pressures) {
    const s = pv(withZ(f, z), P, Tc), V = n * molarV(s), excess = V - Vcell;
    if (s.phase === 'oil' || excess <= 0) { rows.push({ P, liqPct: s.phase === 'oil' ? 100 : (100 * n * (1 - s.beta) * s.oil.vm) / Vcell, produced: 100 * (1 - n), Zgas: null, Z2: (P * 1e5 * Vcell) / (n * R * TK), heavy: null }); continue; }
    const beta = s.phase === 'gas' ? 1 : s.beta, dn = Math.min(excess / s.gas.vm, n * beta), nV = n * beta - dn, nLq = n * (1 - beta);
    z = s.x.map((xi, i) => nV * s.y[i] + nLq * xi); n = nV + nLq; z = z.map((q) => q / n);
    rows.push({ P, liqPct: (100 * nLq * (s.phase === 'gas' ? 0 : s.oil.vm)) / Vcell, produced: 100 * (1 - n), Zgas: s.gas.Z, Z2: (P * 1e5 * Vcell) / (n * R * TK), heavy: heavyOf(s.y) });
  }
  return rows;
}
/**
 * Multi-stage separator test of one mole of feed. stages: [{ p (bara), t (°C) }]; a stock-tank stage at standard conditions is appended.
 * Returns { stages: [{ p, t, gor, gasSG, liqMol }], gor (total, Sm³/Sm³), vSto (m³ per mol feed), rhoSto, api, gasSG (rate-weighted), mwGas }.
 */
export function separatorTest(f, stages = []) {
  const list = [...stages.filter((s) => isNum(+s.p) && isNum(+s.t) && +s.p > P_STD * 1.01).map((s) => ({ p: +s.p, t: +s.t })).sort((a, b) => b.p - a.p), { p: P_STD, t: T_STD }];
  let z = f.z.slice(), nL = 1, last = null; const raw = [];
  for (const st of list) {
    if (!(nL > 0)) { raw.push({ ...st, nG: 0, mw: null, liqMol: 0 }); continue; }
    const s = pv(withZ(f, z), st.p, st.t);
    if (s.phase === 'two') { raw.push({ ...st, nG: nL * s.beta, mw: s.gas.MW, liqMol: nL * (1 - s.beta) }); nL *= 1 - s.beta; z = s.x.slice(); last = s.oil; }
    else if (s.phase === 'oil') { raw.push({ ...st, nG: 0, mw: null, liqMol: nL }); last = s.oil; }
    else { raw.push({ ...st, nG: nL, mw: s.gas.MW, liqMol: 0 }); nL = 0; last = null; }
  }
  const vSto = nL > 0 && last ? nL * last.vm : 0, nG = raw.reduce((s, r) => s + r.nG, 0), mG = raw.reduce((s, r) => s + r.nG * (r.mw || 0), 0);
  return { stages: raw.map((r) => ({ p: r.p, t: r.t, gor: vSto > 0 ? (r.nG * VM_STD) / vSto : null, gasSG: r.mw ? r.mw / MW_AIR : null, liqMol: r.liqMol })), gor: vSto > 0 ? (nG * VM_STD) / vSto : null, vSto, rhoSto: last ? last.rho : null, api: last ? 141.5 / (last.rho / 999.016) - 131.5 : null, mwGas: nG > 0 ? mG / nG : null, gasSG: nG > 0 ? mG / nG / MW_AIR : null };
}
export const INJECTION_GASES = Object.freeze({
  lean: { name: 'Lean hydrocarbon gas', comp: { N2: 1, CO2: 1, C1: 90, C2: 6, C3: 2 } },
  CO2: { name: 'Carbon dioxide', comp: { CO2: 100 } },
  N2: { name: 'Nitrogen', comp: { N2: 100 } },
  rich: { name: 'Rich (enriched) gas', comp: { CO2: 1, C1: 65, C2: 15, C3: 12, nC4: 7 } },
});
/**
 * Swelling test: saturation pressure and swollen volume against moles of injection gas added per mole of original fluid.
 * Returns rows [{ inj (mol/mol), psat, type, swell (saturated volume / original saturated volume) }].
 */
export function swellingTest(spec, opts, gasComp, Tc, fractions = [0, 0.1, 0.2, 0.4, 0.6]) {
  const tot = COMP_IDS.reduce((s, k) => s + (+spec.comp[k] || 0), 0), gt = Object.values(gasComp).reduce((s, v) => s + v, 0), rows = [];
  let v0 = null, guess = null;
  for (const r of fractions) {
    const comp = Object.fromEntries(COMP_IDS.map((k) => [k, (+spec.comp[k] || 0) / tot + (r * (+gasComp[k] || 0)) / gt]));
    const f = makeFluid({ ...spec, comp }, opts), s = saturationPoint(f, Tc, { seed: guess });
    if (s.P === null || s.capped) { rows.push({ inj: r, psat: null, type: null, swell: null }); guess = null; continue; }
    const v = (1 + r) * molarV(pv(f, s.P * 1.00002, Tc)); if (v0 === null) v0 = v;
    rows.push({ inj: r, psat: s.P, type: s.type, swell: v / v0 }); guess = s.X;
  }
  return rows;
}

// ---- aqueous phase: non-electrolyte activity models -------------------------------------------------------------------
// Component 1 = inhibitor, 2 = water. Every model carries a temperature-dependent pair of dimensionless interaction terms
// e12 = a12 + b12/T, e21 = a21 + b21/T (T in K): NRTL τ_ij = e_ij (non-randomness α), Wilson Λ_12 = (V_w/V_1)·exp(−e12), Λ_21 = (V_1/V_w)·exp(−e21),
// UNIQUAC τ_ij = exp(−e_ij). They were regressed in this work to measurements that are NOT hydrate data (see PROVENANCE and ACTIVITY_CALIBRATION):
// MEG — isothermal P–x at 333.15 and 353.15 K (Horstmann et al. 2004) and isobaric boiling points at 50–101 kPa (Chouireb et al. 2018; Kamihama et al. 2012);
// methanol — P–x at 323.15 K (Bernatová et al. 2006), limiting activity coefficients of methanol in water at 273–323 K (Vrbka et al. 2005) and the
// ice line of aqueous methanol (4–34 wt %). `range` is the temperature span of those data plus the extrapolation allowed; T is clamped to it.
export const ACTIVITY_PARAMS = Object.freeze({
  MeOH: { name: 'Methanol', V: 40.73, r: 1.4311, q: 1.432, alpha: 0.3, range: [233, 340], nrtl: { a: [4.92084, -2.98763], b: [-1605.765, 1194.933] }, wilson: { a: [-1.60948, 4.14077], b: [471.836, -1069.907] }, uniquac: { a: [2.10383, -1.6903], b: [-776.333, 672.898] } },
  MEG: { name: 'Mono-ethylene glycol', V: 55.92, r: 2.4088, q: 2.248, alpha: 0.3, range: [240, 475], nrtl: { a: [2.43276, -3.42053], b: [-1069.57, 1426.039] }, wilson: { a: [-4.01808, 2.55337], b: [1280.248, -716.111] }, uniquac: { a: [-2.64535, 3.50825], b: [731.46, -1132.993] } },
});
const V_W = 18.07, R_W = 0.92, Q_W = 1.4;
const actPair = (set, TK, range) => { const T = clamp(TK, range[0], range[1]); return [set.a[0] + set.b[0] / T, set.a[1] + set.b[1] / T]; };
/**
 * Activity coefficients in a binary inhibitor(1)–water(2) solution. model: 'NRTL' | 'UNIQUAC' | 'Wilson' | 'ideal';
 * x1 = inhibitor mole fraction. Returns { g1, g2, aw (water activity), ideal }. Inhibitors without parameters are treated as ideal.
 */
export function activityBinary(model, inhId, x1, TK, params = ACTIVITY_PARAMS) {
  const p = params[inhId], x2 = 1 - x1;
  if (!p || model === 'ideal' || !(x1 > 0)) return { g1: 1, g2: 1, aw: x2, ideal: true };
  if (!(x2 > 0)) return { g1: 1, g2: 1, aw: 0, ideal: false };
  let l1, l2;
  if (model === 'Wilson') {
    const [a12, a21] = actPair(p.wilson, TK, p.range), L12 = (V_W / p.V) * Math.exp(-a12), L21 = (p.V / V_W) * Math.exp(-a21), br = L12 / (x1 + L12 * x2) - L21 / (x2 + L21 * x1);
    l1 = -Math.log(x1 + L12 * x2) + x2 * br; l2 = -Math.log(x2 + L21 * x1) - x1 * br;
  } else if (model === 'UNIQUAC') {
    const [a12, a21] = actPair(p.uniquac, TK, p.range), z = 10, r1 = p.r, q1 = p.q, t12 = Math.exp(-a12), t21 = Math.exp(-a21);
    const phi1 = (x1 * r1) / (x1 * r1 + x2 * R_W), phi2 = 1 - phi1, th1 = (x1 * q1) / (x1 * q1 + x2 * Q_W), th2 = 1 - th1, la = (z / 2) * (r1 - q1) - (r1 - 1), lb = (z / 2) * (R_W - Q_W) - (R_W - 1);
    l1 = Math.log(phi1 / x1) + (z / 2) * q1 * Math.log(th1 / phi1) + phi2 * (la - (r1 / R_W) * lb) - q1 * Math.log(th1 + th2 * t21) + th2 * q1 * (t21 / (th1 + th2 * t21) - t12 / (th2 + th1 * t12));
    l2 = Math.log(phi2 / x2) + (z / 2) * Q_W * Math.log(th2 / phi2) + phi1 * (lb - (R_W / r1) * la) - Q_W * Math.log(th2 + th1 * t12) + th1 * Q_W * (t12 / (th2 + th1 * t12) - t21 / (th1 + th2 * t21));
  } else {
    const [t12, t21] = actPair(p.nrtl, TK, p.range), G12 = Math.exp(-p.alpha * t12), G21 = Math.exp(-p.alpha * t21);
    l1 = x2 * x2 * (t21 * (G21 / (x1 + x2 * G21)) ** 2 + (t12 * G12) / (x2 + x1 * G12) ** 2);
    l2 = x1 * x1 * (t12 * (G12 / (x2 + x1 * G12)) ** 2 + (t21 * G21) / (x1 + x2 * G21) ** 2);
  }
  return { g1: Math.exp(l1), g2: Math.exp(l2), aw: x2 * Math.exp(l2), ideal: false };
}

// ---- aqueous phase: electrolyte models (NaCl) ---------------------------------------------------------------------------
/** Debye–Hückel constants of water at T (K): { Aphi (osmotic, natural log), A (log10, activity), B (1/Å per sqrt(mol/kg)), eps }. */
export function debyeHuckel(TK) {
  const t = clamp(TK - KEL, -30, 150), eps = 87.74 - 0.40008 * t + 9.398e-4 * t * t - 1.41e-6 * t ** 3, rho = rhoBrine(clamp(t, 0, 150), 0);
  const e = 1.602176634e-19, kB = 1.380649e-23, NA = 6.02214076e23, e0 = 8.8541878128e-12, lB = (e * e) / (4 * Math.PI * e0 * eps * kB * TK);
  const Aphi = (Math.sqrt(2 * Math.PI * NA * rho) * lB ** 1.5) / 3;
  return { Aphi, A: (3 * Aphi) / Math.LN10, B: Math.sqrt((2 * e * e * NA * rho) / (e0 * eps * kB * TK)) * 1e-10, eps };
}
// Pitzer parameters of Na+–Cl-: β0, β1 and Cφ with the temperature function of the USGS PHREEQC database pitzer.dat
// (P = a0 + a1(1/T − 1/Tr) + a2 ln(T/Tr) + a3(T − Tr) + a4(T² − Tr²) + a5(1/T² − 1/Tr²), Tr = 298.15 K), read 2026-10-08.
export const PITZER_NACL = Object.freeze({ b0: [7.534e-2, 9598.4, 35.48, -5.8731e-2, 1.798e-5, -5e5], b1: [0.2769, 1.377e4, 46.8, -6.9512e-2, 2e-5, -7.4823e5], c: [1.48e-3, -120.5, -0.2081, 0, 1.166e-7, 11121] });
const pitzerNaCl = (TK) => { const Tr = 298.15, f = (a) => a[0] + a[1] * (1 / TK - 1 / Tr) + a[2] * Math.log(TK / Tr) + a[3] * (TK - Tr) + a[4] * (TK * TK - Tr * Tr) + a[5] * (1 / (TK * TK) - 1 / (Tr * Tr)); return { b0: f(PITZER_NACL.b0), b1: f(PITZER_NACL.b1), c: f(PITZER_NACL.c) }; };
/** Mean ionic activity coefficient of aqueous NaCl at molality m. model: 'dh' | 'edh' | 'davies' | 'pitzer'. */
export function gammaNaCl(model, m, TK = 298.15) {
  if (!(m > 0)) return 1;
  const dh = debyeHuckel(TK), I = m, s = Math.sqrt(I);
  if (model === 'dh') return 10 ** (-dh.A * s);
  if (model === 'edh') return 10 ** ((-dh.A * s) / (1 + dh.B * 4 * s)); // ion-size parameter 4 Å
  if (model === 'davies') return 10 ** (-dh.A * (s / (1 + s) - 0.3 * I));
  const p = pitzerNaCl(TK), b = 1.2, al = 2, fg = -dh.Aphi * (s / (1 + b * s) + (2 / b) * Math.log(1 + b * s));
  const Bg = 2 * p.b0 + ((2 * p.b1) / (al * al * I)) * (1 - (1 + al * s - (al * al * I) / 2) * Math.exp(-al * s));
  return Math.exp(fg + m * Bg + 1.5 * m * m * p.c);
}
/** Osmotic coefficient of aqueous NaCl: analytic for Pitzer, by Gibbs–Duhem integration of ln γ± for the Debye–Hückel family. */
export function osmoticNaCl(model, m, TK = 298.15) {
  if (!(m > 0)) return 1;
  if (model === 'pitzer') { const dh = debyeHuckel(TK), p = pitzerNaCl(TK), s = Math.sqrt(m); return 1 - (dh.Aphi * s) / (1 + 1.2 * s) + m * (p.b0 + p.b1 * Math.exp(-2 * s)) + m * m * p.c; }
  // φ = 1 + ln γ(m) − (1/m) ∫0^m ln γ dm', integrated in s = sqrt(m') by Simpson's rule
  const n = 200, sm = Math.sqrt(m), h = sm / n, g = (s) => (s > 0 ? Math.log(gammaNaCl(model, s * s, TK)) * 2 * s : 0);
  let q = g(0) + g(sm); for (let i = 1; i < n; i++) q += (i % 2 ? 4 : 2) * g(i * h);
  return 1 + Math.log(gammaNaCl(model, m, TK)) - ((q * h) / 3) / m;
}
/** Water activity of an NaCl brine of molality m (mol per kg water). */
export const waterActivityNaCl = (model, m, TK = 298.15) => (m > 0 ? Math.exp(-2 * m * MW_W * 1e-3 * osmoticNaCl(model, m, TK)) : 1);
/**
 * Water activity of the aqueous phase: salinity S (g NaCl-equivalent per kg brine), inhibitor mass fraction inhWt (wt % of the
 * aqueous phase). The electrolyte and inhibitor contributions are multiplied (ln a_w additive), the salt molality being referred
 * to the water of the brine. Returns { aw, awSalt, awInh, xw, x1, molality, gammaW }.
 */
export function waterActivity({ S = 0, inhId = 'none', inhWt = 0 }, TK, { act = 'NRTL', elec = 'pitzer' } = {}, tune = null) {
  const s = clamp(S, 0, 260) / 1000, m = s > 0 ? (s / (MW_NACL * 1e-3)) / (1 - s) : 0, awSalt = waterActivityNaCl(elec, m, TK) ** (tune?.salt ?? 1);
  const w = inhId && inhId !== 'none' ? clamp(inhWt, 0, 95) / 100 : 0, mwI = (INHIBITORS[inhId] || INHIBITORS.none).MW;
  // salt-free mole fraction of inhibitor: w kg inhibitor with (1 - w)(1 - s) kg water
  const nI = w / mwI, nW = ((1 - w) * (1 - s)) / MW_W, x1 = nI + nW > 0 ? nI / (nI + nW) : 0, a0 = activityBinary(act, inhId, x1, TK);
  const a = tune && tune.act !== undefined && tune.act !== 1 && !a0.ideal ? { ...a0, g2: a0.g2 ** tune.act, aw: (1 - x1) * a0.g2 ** tune.act } : a0; // calibration multiplier on ln γ of water
  return { aw: awSalt * a.aw, awSalt, awInh: a.aw, xw: 1 - x1, x1, molality: m, gammaW: a.g2, idealInh: w > 0 && a.ideal && act !== 'ideal' };
}

// ---- Henry's law, Raoult's law and Poynting correction ------------------------------------------------------------------
// Henry solubility constants at 25 °C (mol/kg/bar), van 't Hoff temperature coefficient (K), partial molar volume at infinite
// dilution (cm³/mol) and Setschenow salting-out constant for NaCl (kg/mol): Sander compilation and Wilhelm et al.
const HENRY = { C1: [1.4e-3, 1600, 37, 0.127], C2: [1.9e-3, 2300, 51, 0.162], C3: [1.5e-3, 2700, 67, 0.194], N2: [6.4e-4, 1300, 33, 0.121], CO2: [3.4e-2, 2400, 33, 0.1], H2S: [1.0e-1, 2100, 35, 0.064], nC4: [1.2e-3, 3100, 77, 0.217], iC4: [8.4e-4, 2700, 77, 0.217] };
/**
 * Gas solubility in the aqueous phase by Henry's law with Poynting correction and Setschenow salting-out.
 * fug: { id: fugacity (bar) }; returns { x: { id: mole fraction }, total (mole fraction), sm3 (Sm³ gas per m³ of water) }.
 */
export function gasSolubility(fug, Pbar, Tc, molality = 0) {
  const TK = Tc + KEL, x = {}; let tot = 0;
  for (const [id, f] of Object.entries(fug)) {
    const h = HENRY[id]; if (!h || !(f > 0)) continue;
    const kH = h[0] * Math.exp(h[1] * (1 / TK - 1 / 298.15)), poy = Math.exp((h[2] * 1e-6 * (Pbar - 1) * 1e5) / (R * TK)), mol = (kH * f) / poy / 10 ** (h[3] * molality);
    x[id] = mol / (1000 / MW_W + mol); tot += x[id];
  }
  return { x, total: tot, sm3: (tot / (1 - tot)) * (1000 / MW_W) * 1000 * VM_STD };
}
/** Water content of gas (kg/Sm³) by Raoult's law with the Poynting correction and water activity (ideal-gas vapour phase). */
export function waterContentRaoult(Pbar, Tc, aw = 1) {
  const TK = Tc + KEL, ps = psatWater(clamp(Tc, 0.01, 370)) / 1e5, y = clamp((aw * ps * Math.exp((MW_W * 1e-6 * (Pbar - ps) * 1e5) / (R * TK))) / Pbar, 0, 0.999);
  return ((y / (1 - y)) * MW_W * 1e-3) / VM_STD;
}

// ---- gas hydrates: van der Waals–Platteeuw ---------------------------------------------------------------------------------
// Langmuir constants C = (A/T) exp(B/T) in 1/atm, [A (K/atm), B (K)] for the small and large cavities of structures I and II,
// and reference properties of the empty lattice relative to water — the parameter set of Munck, Skjold-Jørgensen & Rasmussen (1988).
const T0H = 273.15, SH_NU = [3 / 34, 2 / 34], SH_FIT = [528.66, -5623.6];
export const HYDRATE_STRUCTURES = Object.freeze({
  sI: { nu: [1 / 23, 3 / 23], dmu0: 1264, dhL: -4858, dhI: 1151, dvL: 4.6e-6, dvI: 3.0e-6, waters: 23 },
  sII: { nu: [2 / 17, 1 / 17], dmu0: 883, dhL: -5201, dhI: 808, dvL: 5.0e-6, dvI: 3.4e-6, waters: 17 },
  sH: { nu: SH_NU, dmu0: SH_FIT[0], dhL: SH_FIT[1], dhI: SH_FIT[1] + 6009.5, dvL: 3.77e-6, dvI: 2.14e-6, waters: 34 },
});
const DCP_L = -39.16; // J/mol/K, empty lattice minus liquid water
export const LANGMUIR = Object.freeze({
  C1: { sI: [[0.7228e-3, 3187], [23.35e-3, 2653]], sII: [[0.2207e-3, 3453], [100e-3, 1916]] },
  C2: { sI: [[0, 0], [3.039e-3, 3861]], sII: [[0, 0], [240e-3, 2967]] },
  C3: { sI: [[0, 0], [0, 0]], sII: [[0, 0], [5.455e-3, 4638]] },
  iC4: { sI: [[0, 0], [0, 0]], sII: [[0, 0], [189.3e-3, 3800]] },
  nC4: { sI: [[0, 0], [0, 0]], sII: [[0, 0], [30.51e-3, 3699]] },
  N2: { sI: [[1.617e-3, 2905], [6.078e-3, 2431]], sII: [[0.1742e-3, 3082], [18e-3, 1728]] },
  CO2: { sI: [[0.2474e-3, 3410], [42.46e-3, 2813]], sII: [[0.0845e-3, 3615], [851e-3, 2025]] },
  H2S: { sI: [[0.025e-3, 4568], [16.34e-3, 3737]], sII: [[0.0298e-3, 4878], [87.2e-3, 2633]] },
});
/** Langmuir constant (1/bar) of guest `id` in cavity m (0 small, 1 large) of structure s at T (K). */
export const langmuirC = (id, s, m, TK, table = LANGMUIR) => {
  if (table.kihara || s === 'sH') { const k = KIHARA.guests[id], cv = KIHARA.cav[s]?.[m]; if (k && cv && !KIHARA.munck[id]?.includes(s)) return s !== 'sH' || SH_HELP.includes(id) ? kiharaLangmuir(id, s, m, TK) : 0; if (s === 'sH') return 0; table = LANGMUIR; } // guests without a Kihara set keep the Munck constants
  const c = table[id]?.[s]?.[m]; return c && c[0] > 0 ? ((c[0] / TK) * Math.exp(c[1] / TK)) / ATM : 0;
};
let HYD_ACTIVE = { table: null, sH: false }; // constants set and structure-H switch of the run in progress (null: the default Kihara set)
const SH_HELP = ['C1', 'N2']; // small help gases admitted to the 5¹² and 4³5⁶6³ cavities of structure H
/**
 * Chemical-potential differences of water (divided by RT) at T (K), P (bara) for hydrate structure s:
 * hydrate side Σ ν ln(1 + Σ C f) and water side (liquid with activity aw, or ice when that is the stable phase).
 * fug: { id: fugacity in bar }. Returns { dmuH, dmuW, ice, theta: [{ id: occupancy } small, large], drive = dmuH − dmuW (> 0: hydrate stable) }.
 */
export function hydrateState(s, TK, Pbar, fug, aw = 1, table = LANGMUIR) {
  const S = table.ref?.[s] ? { ...HYDRATE_STRUCTURES[s], ...table.ref[s] } : HYDRATE_STRUCTURES[s], theta = [{}, {}], dcp = (s !== 'sH' && table.dcp) || [DCP_L, 0]; let dmuH = 0;
  for (let m = 0; m < 2; m++) {
    let sm = 0; const cf = {};
    for (const id in fug) { const v = langmuirC(id, s, m, TK, table) * fug[id]; if (v > 0) { cf[id] = v; sm += v; } }
    for (const id in cf) theta[m][id] = cf[id] / (1 + sm);
    dmuH += S.nu[m] * Math.log(1 + sm);
  }
  const P = Pbar * 1e5, Tm = 0.5 * (TK + T0H), base = S.dmu0 / (R * T0H), inv = 1 / T0H - 1 / TK;
  // ∫ Δh/(RT²) dT with Δcp = a + b(T − T0): Δh = Δh0 + a(T − T0) + b(T − T0)²/2
  const liq = base - (((S.dhL - dcp[0] * T0H + (dcp[1] * T0H * T0H) / 2) / R) * inv + ((dcp[0] - dcp[1] * T0H) / R) * Math.log(TK / T0H) + (dcp[1] / (2 * R)) * (TK - T0H)) + (S.dvL * P) / (R * Tm) - Math.log(Math.max(aw, 1e-6));
  const ice = base - (S.dhI / R) * inv + (S.dvI * P) / (R * Tm);
  const dmuW = Math.max(liq, ice);
  return { dmuH, dmuW, ice: ice > liq, theta, drive: dmuH - dmuW };
}
/** Fugacities (bar) of the hydrate formers in the hydrocarbon system at P (bara), T (°C) — vapour phase when two phases coexist. */
export function formerFugacities(f, Pbar, Tc) {
  const fl = flashPT(f, Pbar, Tc), x = fl.phase === 'two' ? fl.y : f.z, ph = fl.phase === 'two' ? fl.vap : fl.liq, fug = {};
  f.comps.forEach((c, i) => { if (LANGMUIR[c.id] && x[i] > 0) fug[c.id] = x[i] * Math.exp(ph.lnphi[i]) * Pbar; });
  return fug;
}
const occupancy = (s, st) => {
  const S = HYDRATE_STRUCTURES[s], tot = st.theta.map((t) => Object.values(t).reduce((a, b) => a + b, 0)), guests = S.nu[0] * tot[0] + S.nu[1] * tot[1];
  return { small: tot[0], large: tot[1], byGuest: st.theta, hydrationNumber: guests > 0 ? 1 / guests : null };
};
/**
 * Hydrate dissociation temperature (°C) at P (bara) for an EOS fluid: both structures are solved and the stable one (highest
 * temperature, i.e. lowest pressure) is returned. aw: water activity of the gas-free aqueous phase, a number or a function of
 * T (K); the dissolved gas (Henry's law, salting-out at `molality`) lowers it further. The guest fugacities are re-flashed at the
 * solution until the temperature settles (`passes` = 1 keeps the fugacities of the guess temperature). Returns { T, structure, TsI, TsII, ice, occupancy, xGas, fug } or null.
 */
export function hydrateTofP(f, Pbar, aw = 1, { Tmin = -45, Tmax = 50, table = HYD_ACTIVE.table || KIHARA, molality = 0, solubility = true, guess = 10, passes = 10, sH = HYD_ACTIVE.sH } = {}) {
  const structs = sH ? ['sI', 'sII', 'sH'] : ['sI', 'sII'];
  const awAt = typeof aw === 'function' ? aw : () => aw;
  const awTot = (fug, t) => awAt(t + KEL) * (1 - (solubility ? gasSolubility(fug, Pbar, t, molality).total : 0));
  const solve = (fug) => {
    const res = {};
    for (const s of structs) {
      const g = (t) => hydrateState(s, t + KEL, Pbar, fug, awTot(fug, t), table).drive;
      if (!(g(Tmin) > 0)) { res[s] = null; continue; }
      if (g(Tmax) > 0) { res[s] = Tmax; continue; }
      res[s] = brent(g, Tmin, Tmax, 2e-4, 60); // the driving force falls monotonically with temperature
    }
    return res;
  };
  let T = clamp(guess, Tmin, Tmax), fug = formerFugacities(f, Pbar, T), res = null;
  for (let it = 0; it < passes; it++) {
    if (!Object.keys(fug).length) return null;
    res = solve(fug);
    if (structs.every((k) => res[k] === null)) return null;
    const Tn = Math.max(...structs.map((k) => res[k] ?? -1e9)), done = Math.abs(Tn - T) < 5e-3;
    T = it > 5 ? 0.5 * (T + Tn) : Tn;
    if (done || it === passes - 1) break;
    fug = formerFugacities(f, Pbar, T);
  }
  const structure = structs.reduce((b, k) => ((res[k] ?? -1e9) > (res[b] ?? -1e9) ? k : b), 'sI'), st = hydrateState(structure, T + KEL, Pbar, fug, awTot(fug, T), table);
  return { T, structure, TsI: res.sI, TsII: res.sII, TsH: res.sH ?? null, ice: st.ice, occupancy: occupancy(structure, st), xGas: solubility ? gasSolubility(fug, Pbar, T, molality).total : 0, fug };
}
/** Lowest hydrate dissociation pressure (bara) at T (°C); null when no hydrate forms below Pmax. Returns { P, structure, occupancy }. */
export function hydratePofT(f, Tc, aw = 1, { Pmin = 0.2, Pmax = 1500, table = HYD_ACTIVE.table || KIHARA, molality = 0, solubility = true, sH = HYD_ACTIVE.sH } = {}) {
  const a = typeof aw === 'function' ? aw(Tc + KEL) : aw, TK = Tc + KEL, structs = sH ? ['sI', 'sII', 'sH'] : ['sI', 'sII'];
  const both = (P) => { const fug = formerFugacities(f, P, Tc), w = a * (1 - (solubility ? gasSolubility(fug, P, Tc, molality).total : 0)); return structs.map((s) => hydrateState(s, TK, P, fug, w, table)); };
  const drive = (P) => Math.max(...both(P).map((q) => q.drive));
  const grid = logspace(Pmin, Pmax, 36);
  if (drive(grid[0]) > 0) return { P: Pmin, structure: null, occupancy: null };
  let k = 1; while (k < grid.length && !(drive(grid[k]) > 0)) k++;
  if (k === grid.length) return null;
  let lo = Math.log(grid[k - 1]), hi = Math.log(grid[k]);
  for (let i = 0; i < 50 && hi - lo > 1e-6; i++) { const m = 0.5 * (lo + hi); if (drive(Math.exp(m)) > 0) hi = m; else lo = m; }
  const P = Math.exp(0.5 * (lo + hi)), all = both(P), ib = all.reduce((b, q, k) => (q.drive > all[b].drive ? k : b), 0), structure = structs[ib];
  return { P, structure, occupancy: occupancy(structure, all[ib]) };
}
/**
 * Freezing point (°C) of an aqueous phase { S, inhId, inhWt }: temperature at which the water activity of the solution equals that of ice
 * (enthalpy of fusion and heat-capacity difference as implied by the liquid and ice reference properties of the hydrate model). Null below −60 °C.
 */
export function freezingPoint(aq, models = {}, tune = null) {
  const S = HYDRATE_STRUCTURES.sI, lnIce = (T) => ((S.dhI - S.dhL + DCP_L * T0H) / R) * (1 / T0H - 1 / T) - (DCP_L / R) * Math.log(T / T0H);
  const g = (T) => lnIce(T) - Math.log(waterActivity(aq, T, models, tune).aw);
  if (!(g(T0H) > 0)) return 0; if (g(213.15) > 0) return null;
  return brent(g, 213.15, T0H, 1e-8, 80) - KEL;
}
/**
 * Conservative allowance on the hydrate depression credited to a thermodynamic inhibitor: the fraction of the model depression that is
 * NOT credited in the curve handed to the other suites. Zero up to 30 wt %, rising linearly to 0.12 at 50 wt % and held above: the largest
 * amount by which the model has been shown to over-predict the depression of any measured set (methane–propane with 50 wt % MEG, Song & Kobayashi 1989).
 */
export const inhibitorAllowance = (inhWt) => 0.12 * clamp(((+inhWt || 0) - 30) / 20, 0, 1);
/**
 * Conservative allowance (K) added to every published hydrate temperature at high pressure: above 250 bar the model drifts cold against the measured
 * methane, nitrogen and natural-gas points (about −0.4 K at 350 bar, −1.4 K at 500 bar, −1.9 K at 600–730 bar): 0.004 K per bar above 250 bar, at most 3 K.
 */
export const pressureAllowance = (Pbar) => clamp(0.004 * ((+Pbar || 0) - 250), 0, 3);
/**
 * Published (conservative) hydrate temperature from the model value Tm, the model fresh-water value T0 at the same pressure, the inhibitor content
 * (wt % of the aqueous phase) and the pressure (bara): model value plus the inhibitor allowance on the depression plus the high-pressure allowance.
 */
export const publishedHydrateT = (Tm, T0, inhWt = 0, Pbar = 0) => Tm + (Number.isFinite(T0) && T0 > Tm ? inhibitorAllowance(inhWt) * (T0 - Tm) : 0) + pressureAllowance(Pbar);
/** Hammerschmidt hydrate depression (°C) of w wt % inhibitor. */
export const hammerschmidt = (w, inhId) => { const i = INHIBITORS[inhId]; return i && w > 0 && w < 100 ? (i.K * w) / (i.MW * (100 - w)) : 0; };

// ---- wax (ideal-solution solid–liquid equilibrium) and asphaltene screening ------------------------------------------------
/** Won (1986) melting temperature (K) and enthalpy of fusion (J/mol) of a paraffin of molar mass M. */
export const wonFusion = (M) => { const Tf = 374.5 + 0.02617 * M - 20172 / M; return { Tf, dH: 0.1426 * M * Tf * 4.184 }; };
/** Pedersen (1995) wax-forming fraction of a carbon-number cut of molar mass M and density SG. */
export const waxFormingFraction = (M, SG) => { const rp = 0.3915 + 0.0675 * Math.log(M); return SG <= rp ? 1 : clamp(1 - (1.074 + 6.584e-4 * M) * ((SG - rp) / rp) ** 0.1915, 0, 1); };
/**
 * Wax model of an oil. dist: SCN distribution (z summing to the C7+ mole fraction of the whole fluid); x7: C7+ mole fraction of
 * the liquid considered; mwLiq: its molar mass. Returns { wat(detect wt %) -> °C | null, solidWt(T °C) -> wt % of the liquid, waxContent (wt % wax-forming C18+) }.
 */
export function waxModel(dist, x7, mwLiq, { hfMult = 1, nMin = 15 } = {}) {
  const z7 = sumA(dist.map((d) => d.z));
  if (!(z7 > 0) || !(x7 > 0)) return { wat: () => null, solidWt: () => 0, waxContent: 0, comps: [] };
  const comps = dist.filter((d) => d.n >= nMin).map((d) => { const fu = wonFusion(d.M); return { n: d.n, M: d.M, z: ((x7 * d.z) / z7) * waxFormingFraction(d.M, d.SG), Tf: fu.Tf, dH: fu.dH * hfMult }; }).filter((c) => c.z > 1e-14);
  const zw = sumA(comps.map((c) => c.z)), mw = sumA(comps.map((c) => c.z * c.M)), zs = 1 - zw, M = Math.max(mwLiq, mw + 1e-9);
  const Ks = (TK) => comps.map((c) => Math.exp((c.dH / (R * TK)) * (1 - TK / c.Tf)));
  const solidWt = (Tc) => {
    const K = Ks(Tc + KEL), g = (S) => { let s = -zs / (1 - S); comps.forEach((c, i) => (s += (c.z * (K[i] - 1)) / (1 + S * (K[i] - 1)))); return s; };
    if (!(g(0) > 0)) return 0;
    const S = brent(g, 0, 1 - 1e-12, 1e-14); let ms = 0;
    comps.forEach((c, i) => (ms += ((c.z * K[i]) / (1 + S * (K[i] - 1))) * S * c.M));
    return clamp((100 * ms) / M, 0, 100);
  };
  const wat = (detect = 0.02) => { const g = (t) => solidWt(t) - detect; if (!(g(-60) > 0)) return null; if (g(150) > 0) return 150; let lo = -60, hi = 150; for (let k = 0; k < 50 && hi - lo > 1e-3; k++) { const m = 0.5 * (lo + hi); if (g(m) > 0) lo = m; else hi = m; } return 0.5 * (lo + hi); };
  return { wat, solidWt, waxContent: (100 * sumA(comps.filter((c) => c.n >= 18).map((c) => c.z * c.M))) / M, comps };
}
/**
 * Asphaltene screening. de Boer: undersaturation (Pres − Pbub, bar) against in-situ density (kg/m³) with approximate
 * digitised boundaries of the published plot; colloidal instability index CII = (saturates + asphaltenes)/(aromatics + resins).
 * Returns { deBoer: 'none' | 'slight' | 'severe', dP, limitSlight, limitSevere, cii, ciiClass, risk: 'low' | 'medium' | 'high' }.
 */
export function asphalteneScreen({ rhoRes, pRes, pBub, sara }) {
  const dP = Math.max(0, pRes - (pBub ?? pRes)), limitSlight = 95 * Math.exp((clamp(rhoRes, 450, 1000) - 600) / 125), limitSevere = 1.7 * limitSlight;
  const deBoer = dP > limitSevere ? 'severe' : dP > limitSlight ? 'slight' : 'none';
  const s = sara || {}, den = (+s.aro || 0) + (+s.res || 0), cii = den > 0 ? ((+s.sat || 0) + (+s.asp || 0)) / den : null, ciiClass = cii === null ? 'unknown' : cii >= 0.9 ? 'unstable' : cii >= 0.7 ? 'uncertain' : 'stable';
  const score = (deBoer === 'severe' ? 2 : deBoer === 'slight' ? 1 : 0) + (ciiClass === 'unstable' ? 2 : ciiClass === 'uncertain' ? 1 : 0);
  return { deBoer, dP, limitSlight, limitSevere, cii, ciiClass, risk: score >= 3 ? 'high' : score >= 1 ? 'medium' : 'low' };
}

// ---- residual-Helmholtz equations of state: GERG-2008, PC-SAFT, cubic-plus-association ------------------------------------
// A Helmholtz model is { ids, n, M[] (g/mol), crit[] ({ Tc, Pc, w }), rhoMax(x), ar(T, rho, x) } where ar is the residual Helmholtz
// energy per mole divided by RT at temperature T (K), molar density rho (mol/m³) and mole fractions x. Every property follows
// from ar and its derivatives (taken numerically here): Z = 1 + ρ ∂ar/∂ρ, ln φ_i = ∂(n ar)/∂n_i − ln Z, h_res = RT(Z − 1 − T ∂ar/∂T).
const N_AV = 6.02214076e23;
const d1 = (g, x, h) => (g(x - 2 * h) - 8 * g(x - h) + 8 * g(x + h) - g(x + 2 * h)) / (12 * h); // fourth-order central difference
/** Compressibility factor of a Helmholtz model at T (K), rho (mol/m³), composition x. */
export const hZ = (m, T, rho, x) => 1 + rho * d1((r) => m.ar(T, r, x), rho, rho * 2e-4);
const hP = (m, T, rho, x) => rho * R * T * hZ(m, T, rho, x); // Pa
/**
 * Molar density (mol/m³) of a Helmholtz model at T (K), P (bara). kind 'vapour' walks up from the ideal gas, 'liquid' walks down
 * from the close-packed side; each returns null when that branch does not reach the pressure with (∂P/∂ρ)T > 0. 'stable' takes the
 * root of lower Gibbs energy.
 */
function hRoot(m, T, Pbar, x, branch) { // strict branch root: 'vapour' from the ideal-gas side, 'liquid' from the dense side
  const P = Pbar * 1e5, top = m.rhoMax(x), g = (r) => hP(m, T, r, x) - P;
  const walk = (r0, fac) => {
    let a = r0, ga = g(a);
    for (let k = 0; k < 500; k++) {
      const b = a * fac; if (b > top || b < 1e-9) return null;
      const gb = g(b);
      if (!Number.isFinite(gb) || (gb - ga) * (b - a) <= 0) return null; // mechanical stability lost: this branch ends here
      if (ga * gb <= 0) return brent(g, Math.min(a, b), Math.max(a, b), 1e-13 * b, 80);
      a = b; ga = gb;
    }
    return null;
  };
  if (branch === 'vapour') return walk(Math.min((P / (R * T)) * 0.25, top * 0.02), 1.12);
  let r0 = top, g0 = g(r0); for (let k = 0; k < 80 && !(g0 > 0 && Number.isFinite(g0)); k++) { r0 *= 0.97; g0 = g(r0); }
  return g0 > 0 ? walk(r0, 0.975) : null;
}
export function hRho(m, T, Pbar, x, kind = 'stable') {
  if (kind === 'vapour') return hRoot(m, T, Pbar, x, 'vapour') ?? hRoot(m, T, Pbar, x, 'liquid');
  if (kind === 'liquid') return hRoot(m, T, Pbar, x, 'liquid') ?? hRoot(m, T, Pbar, x, 'vapour');
  const v = hRoot(m, T, Pbar, x, 'vapour'), l = hRoot(m, T, Pbar, x, 'liquid');
  if (v === null || l === null || Math.abs(v / l - 1) < 1e-6) return v ?? l;
  const gib = (r) => { const Z = (Pbar * 1e5) / (r * R * T); return m.ar(T, r, x) + Z - 1 - Math.log(Z); };
  return gib(l) < gib(v) ? l : v;
}
/** Fugacity coefficients (natural log) of a Helmholtz model at T, rho, x. */
export function hLnPhi(m, T, rho, x, only = null) { // only: optional list of component indices (the others are returned as NaN)
  const lnZ = Math.log(hZ(m, T, rho, x)), out = new Array(m.n).fill(NaN), e = 1e-4;
  const F = (i, d) => { const N = 1 + d, xx = x.map((v, k) => (v + (k === i ? d : 0)) / N); return N * m.ar(T, rho * N, xx); }; // n·ar at constant T and V
  for (let i = 0; i < m.n; i++) if (!only || only.includes(i)) out[i] = (x[i] > 2 * e ? d1((d) => F(i, d), 0, e) : (F(i, 2 * e) * -1 + 4 * F(i, e) - 3 * F(i, 0)) / (2 * e)) - lnZ;
  return out;
}
/**
 * Thermal and volumetric properties of a Helmholtz model at T, rho, x: { Z, P (bara), hRes (J/mol), cvRes (J/mol/K), dPdrho, dPdT }
 * and, when the ideal-gas heat capacity cp0 (J/mol/K) is given, cv, cp (J/mol/K), speed of sound w (m/s) and Joule–Thomson coefficient jt (K/bar).
 */
export function hProps(m, T, rho, x, cp0 = null) {
  const hT = T * 1e-3, a = (t) => m.ar(t, rho, x), a0 = a(T), aT = d1(a, T, hT), aTT = (-a(T - 2 * hT) + 16 * a(T - hT) - 30 * a0 + 16 * a(T + hT) - a(T + 2 * hT)) / (12 * hT * hT);
  const Z = hZ(m, T, rho, x), dPdrho = d1((r) => hP(m, T, r, x), rho, rho * 2e-4), dPdT = d1((t) => hP(m, t, rho, x), T, hT);
  let M = 0; for (let i = 0; i < m.n; i++) M += x[i] * m.M[i];
  const o = { Z, P: (rho * R * T * Z) / 1e5, rho, rhoMass: rho * M * 1e-3, M, hRes: R * T * (Z - 1 - T * aT), cvRes: -R * (2 * T * aT + T * T * aTT), dPdrho, dPdT };
  if (cp0 !== null) { o.cv = cp0 - R + o.cvRes; o.cp = o.cv + (T * dPdT * dPdT) / (rho * rho * dPdrho); o.w = Math.sqrt(Math.max(((o.cp / o.cv) * dPdrho) / (M * 1e-3), 0)); o.jt = (((T * dPdT) / (rho * rho * dPdrho) - 1 / rho) / o.cp) * 1e5; }
  return o;
}
/** Vapour pressure (bara) of pure component i of a Helmholtz model at T (K) by equal fugacity; { P, rhoL, rhoV } or null. */
export function hPsat(m, T, i = 0) {
  const c = m.crit[i], x = m.ids.map((_, k) => (k === i ? 1 : 0));
  if (!(T < c.Tc)) return null;
  let P = c.Pc * Math.exp(5.373 * (1 + c.w) * (1 - c.Tc / T)), rl = null, rv = null;
  for (let it = 0; it < 120; it++) {
    rl = hRoot(m, T, P, x, 'liquid'); rv = hRoot(m, T, P, x, 'vapour');
    if (rl === null) { P *= 1.06; continue; } // below the liquid spinodal pressure
    if (rv === null) { P *= 0.94; continue; } // above the vapour spinodal pressure
    if (rl / rv < 1.0005) return null;
    const k = Math.exp(hLnPhi(m, T, rl, x)[i] - hLnPhi(m, T, rv, x)[i]); // φL/φV → 1 at saturation
    P *= k; if (Math.abs(k - 1) < 1e-9) return { P, rhoL: rl, rhoV: rv };
  }
  return null;
}
/** Two-phase PT flash of a Helmholtz model by successive substitution from Wilson K-values: { phase, beta, x, y, rhoL, rhoV, K }. */
export function hFlash(m, T, Pbar, z, K0 = null) {
  let K = K0 || m.crit.map((c) => (c.Pc / Pbar) * Math.exp(5.373 * (1 + c.w) * (1 - c.Tc / T))), beta = 0.5, x = z, y = z, rl = null, rv = null;
  const single = () => { const r = hRho(m, T, Pbar, z, 'stable'), liq = r / m.rhoMax(z) > 0.3; return { phase: liq ? 'liquid' : 'gas', beta: liq ? 0 : 1, x: z.slice(), y: z.slice(), rhoL: r, rhoV: r, K: z.map(() => 1) }; };
  for (let it = 0; it < 300; it++) {
    beta = rachfordRice(z, K);
    if (beta <= 0 || beta >= 1) { const sK = z.reduce((s, v, i) => s + v * K[i], 0), sI = z.reduce((s, v, i) => s + v / K[i], 0); if (sK <= 1 || sI <= 1) return single(); }
    x = z.map((v, i) => v / (1 + beta * (K[i] - 1))); y = x.map((v, i) => v * K[i]);
    const sx = sumA(x), sy = sumA(y); x = x.map((v) => v / sx); y = y.map((v) => v / sy);
    rl = hRoot(m, T, Pbar, x, 'liquid'); rv = hRoot(m, T, Pbar, y, 'vapour');
    if (rl === null || rv === null) return single();
    const pl = hLnPhi(m, T, rl, x), pvp = hLnPhi(m, T, rv, y); let err = 0, tr = 0;
    for (let i = 0; i < m.n; i++) { const Kn = Math.exp(pl[i] - pvp[i]); err += (Kn / K[i] - 1) ** 2; K[i] = Kn; tr += Math.log(Kn) ** 2; }
    if (tr < 1e-8) return single();
    if (err < 1e-14) break;
  }
  return beta > 1e-10 && beta < 1 - 1e-10 ? { phase: 'two', beta, x, y, rhoL: rl, rhoV: rv, K } : single();
}

// GERG-2008 (Kunz & Wagner 2012) for the natural-gas components of the kernel; coefficients in ../data/ref/pvt.js (from NIST teqp).
const GERG_KEY = { C7: 'C6' };
/** GERG-2008 Helmholtz model for kernel component ids (subset of N2, CO2, H2S, C1 … C6); null when a component is not covered. */
export function gergModel(ids) {
  const gi = ids.map((id) => GERG.ids.indexOf(id)); if (gi.some((k) => k < 0)) return null;
  const n = ids.length, P = gi.map((k) => GERG.pure[k]), Tc = P.map((p) => p.Tc), vc = P.map((p) => 1 / p.rhoc), pairs = [];
  for (let a = 0; a < n; a++) for (let b = a + 1; b < n; b++) {
    const [i, j] = gi[a] < gi[b] ? [a, b] : [b, a], q = GERG.pairs[`${ids[i]}-${ids[j]}`]; // i is the component that comes first in the GERG order
    pairs.push({ i, j, bT2: q[2] * q[2], YT: q[2] * q[3] * Math.sqrt(Tc[i] * Tc[j]), bV2: q[0] * q[0], Yv: (q[0] * q[1] * (Math.cbrt(vc[i]) + Math.cbrt(vc[j])) ** 3) / 8, F: q[4] || 0, dep: q[4] ? GERG.dep[q[5]] : null }); // q = [βv, γv, βT, γT, F, departure function]
  }
  const pure = (p, tau, delta) => { const f = GERG.forms[p.form], lt = Math.log(tau), ld = Math.log(delta); let s = 0; for (let k = 0; k < p.n.length; k++) s += p.n[k] * Math.exp(f.t[k] * lt + f.d[k] * ld - (f.c[k] ? delta ** f.l[k] : 0)); return s; };
  const depf = (d, tau, delta) => { const lt = Math.log(tau), ld = Math.log(delta); let s = 0; for (let k = 0; k < d.n.length; k++) s += d.n[k] * Math.exp(d.t[k] * lt + d.d[k] * ld - d.eta[k] * (delta - d.epsilon[k]) ** 2 - d.beta[k] * (delta - d.gamma[k])); return s; };
  const reduce = (x) => { let Tr = 0, vr = 0; for (let i = 0; i < n; i++) { Tr += x[i] * x[i] * Tc[i]; vr += x[i] * x[i] * vc[i]; } for (const p of pairs) { const xi = x[p.i], xj = x[p.j]; if (!(xi > 0 && xj > 0)) continue; const c = 2 * xi * xj * (xi + xj); Tr += (c / (p.bT2 * xi + xj)) * p.YT; vr += (c / (p.bV2 * xi + xj)) * p.Yv; } return { Tr, vr }; };
  return {
    name: 'GERG-2008', ids, n, M: P.map((p) => p.M), crit: ids.map((id) => COMPONENTS[id]), reduce,
    rhoMax: (x) => 3.9 / reduce(x).vr,
    ar(T, rho, x) { const { Tr, vr } = reduce(x), tau = Tr / T, delta = rho * vr; let s = 0; for (let i = 0; i < n; i++) if (x[i] > 0) s += x[i] * pure(P[i], tau, delta); for (const p of pairs) if (p.F && x[p.i] > 0 && x[p.j] > 0) s += x[p.i] * x[p.j] * p.F * depf(p.dep, tau, delta); return s; },
  };
}

// PC-SAFT (Gross & Sadowski 2001): hard-chain reference plus dispersion with the universal constants a, b of the original paper.
// Pure-component segment number m, segment diameter σ (Å) and energy ε/k (K): Gross & Sadowski (2001), H2S from Tihic et al. (2006),
// as tabulated in the open Clapeyron.jl database (PCSAFT_like.csv); binary k_ij from the same database (PCSAFT_unlike.csv).
const PCS_A = [[0.9105631445, 0.6361281449, 2.6861347891, -26.547362491, 97.759208784, -159.59154087, 91.297774084], [-0.3084016918, 0.1860531159, -2.5030047259, 21.419793629, -65.25588533, 83.318680481, -33.74692293], [-0.0906148351, 0.4527842806, 0.5962700728, -1.7241829131, -4.1302112531, 13.77663187, -8.6728470368]];
const PCS_B = [[0.7240946941, 2.2382791861, -4.0025849485, -21.003576815, 26.855641363, 206.55133841, -355.60235612], [-0.5755498075, 0.6995095521, 3.892567339, -17.215471648, 192.67226447, -161.82646165, -165.20769346], [0.0976883116, -0.2557574982, -9.155856153, 20.642075974, -38.804430052, 93.626774077, -29.666905585]];
export const PCSAFT_PARAMS = Object.freeze({
  C1: [1, 3.7039, 150.03], C2: [1.6069, 3.5206, 191.42], C3: [2.002, 3.6184, 208.11], nC4: [2.3316, 3.7086, 222.88], iC4: [2.2616, 3.7574, 216.53], nC5: [2.6896, 3.7729, 231.2], iC5: [2.562, 3.8296, 230.75],
  C6: [3.0576, 3.7983, 236.77], N2: [1.2053, 3.313, 90.96], CO2: [2.0729, 2.7852, 169.21], H2S: [1.6941, 3.0214, 226.79], nC7: [3.4831, 3.8049, 238.4], nC10: [4.6627, 3.8384, 243.87], nC16: [6.6485, 3.9552, 254.7], nC20: [7.9849, 3.9869, 257.75],
});
const PCS_KIJ = { 'C1-nC4': 0.022, 'C1-nC5': 0.024, 'C1-C6': 0.021, 'C1-iC4': 0.028, 'CO2-C1': 0.065, 'CO2-C3': 0.109, 'CO2-nC4': 0.12, 'CO2-nC5': 0.143, 'N2-C6': 0.119 };
// Pseudo-components: linear interpolation in molar mass between the n-alkane sets (n-C7, n-C10, n-C16, n-C20) of m, m·σ³ and m·ε/k.
const PCS_ALK = [[100.2, 'nC7'], [142.29, 'nC10'], [226.45, 'nC16'], [282.55, 'nC20']];
function pcsaftOf(c) {
  if (c.pcs) return c.pcs;
  if (PCSAFT_PARAMS[c.id]) return PCSAFT_PARAMS[c.id];
  if (!c.pseudo) return null;
  let k = 0; while (k < PCS_ALK.length - 2 && c.MW > PCS_ALK[k + 1][0]) k++;
  const [Ma, ia] = PCS_ALK[k], [Mb, ib] = PCS_ALK[k + 1], u = (c.MW - Ma) / (Mb - Ma), A = PCSAFT_PARAMS[ia], B = PCSAFT_PARAMS[ib], li = (f) => f(A) + u * (f(B) - f(A));
  const m = li((p) => p[0]), ms3 = li((p) => p[0] * p[1] ** 3), me = li((p) => p[0] * p[2]), base = [m, Math.cbrt(ms3 / m), me / m];
  if (!(c.SG > 0)) return base;
  // the n-alkane line is lighter than a real cut: the segment diameter is scaled once so that the liquid density at 15.6 °C, 1 atm equals the specific gravity
  const key = `${c.MW.toFixed(3)}|${c.SG.toFixed(5)}`; let sg = pcsCache.get(key);
  if (sg === undefined) { sg = base[1]; for (let it = 0; it < 3; it++) { const one = pcsaftModel([{ id: 'x', MW: c.MW, Tc: c.Tc, Pc: c.Pc, w: c.w, pcs: [base[0], sg, base[2]] }]), r = hRho(one, 288.71, 1.01325, [1], 'liquid'); if (r === null) break; sg *= Math.cbrt((r * c.MW * 1e-3) / (c.SG * 999.016)); } if (pcsCache.size > 500) pcsCache.clear(); pcsCache.set(key, sg); }
  return [base[0], sg, base[2]];
}
const pcsCache = new Map();
/** PC-SAFT Helmholtz model for a list of components ({ id, MW, Tc, Pc, w, pseudo }); null when a component has no parameters. */
export function pcsaftModel(comps) {
  const pr = comps.map(pcsaftOf); if (pr.some((p) => !p)) return null;
  const n = comps.length, ms = pr.map((p) => p[0]), sg = pr.map((p) => p[1]), ek = pr.map((p) => p[2]);
  const kij = comps.map((a) => comps.map((b) => PCS_KIJ[`${a.id}-${b.id}`] ?? PCS_KIJ[`${b.id}-${a.id}`] ?? 0));
  const s3 = comps.map((_, i) => comps.map((__, j) => ((sg[i] + sg[j]) / 2) ** 3)), eij = comps.map((_, i) => comps.map((__, j) => Math.sqrt(ek[i] * ek[j]) * (1 - kij[i][j])));
  const dOf = (T) => sg.map((s, i) => s * (1 - 0.12 * Math.exp((-3 * ek[i]) / T)));
  return {
    name: 'PC-SAFT', ids: comps.map((c) => c.id), n, M: comps.map((c) => c.MW), crit: comps.map((c) => ({ Tc: c.Tc, Pc: c.Pc, w: c.w })), params: pr,
    rhoMax(x) { let s = 0; for (let i = 0; i < n; i++) s += x[i] * ms[i] * sg[i] ** 3; return (0.72 * 6) / (Math.PI * s) / (N_AV * 1e-30); }, // packing fraction 0.72 (close packing is 0.74)
    ar(T, rho, x) {
      const d = dOf(T), rn = rho * N_AV * 1e-30, z = [0, 0, 0, 0]; let mb = 0;
      for (let i = 0; i < n; i++) { mb += x[i] * ms[i]; let p = 1; for (let k = 0; k < 4; k++) { z[k] += x[i] * ms[i] * p; p *= d[i]; } }
      for (let k = 0; k < 4; k++) z[k] *= (Math.PI / 6) * rn;
      const eta = z[3], om = 1 - eta, ahs = ((3 * z[1] * z[2]) / om + z[2] ** 3 / (eta * om * om) + (z[2] ** 3 / (eta * eta) - z[0]) * Math.log(om)) / z[0];
      let ch = 0; for (let i = 0; i < n; i++) if (x[i] > 0 && ms[i] !== 1) { const h = d[i] / 2, g = 1 / om + (h * 3 * z[2]) / (om * om) + (h * h * 2 * z[2] * z[2]) / om ** 3; ch += x[i] * (ms[i] - 1) * Math.log(g); }
      let m2e = 0, m2e2 = 0; for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) { const q = x[i] * x[j] * ms[i] * ms[j] * s3[i][j], e = eij[i][j] / T; m2e += q * e; m2e2 += q * e * e; }
      const c1 = (mb - 1) / mb, c2 = c1 * ((mb - 2) / mb); let I1 = 0, I2 = 0, p = 1;
      for (let k = 0; k < 7; k++) { I1 += (PCS_A[0][k] + c1 * PCS_A[1][k] + c2 * PCS_A[2][k]) * p; I2 += (PCS_B[0][k] + c1 * PCS_B[1][k] + c2 * PCS_B[2][k]) * p; p *= eta; }
      const C1 = 1 / (1 + (mb * (8 * eta - 2 * eta * eta)) / om ** 4 + ((1 - mb) * (20 * eta - 27 * eta * eta + 12 * eta ** 3 - 2 * eta ** 4)) / (om * (2 - eta)) ** 2);
      return mb * ahs - ch - 2 * Math.PI * rn * I1 * m2e - Math.PI * rn * mb * C1 * I2 * m2e2;
    },
  };
}

// Cubic-plus-association (Kontogeorgis et al. 1996; simplified CPA of Kontogeorgis et al. 1999): SRK plus the Wertheim association
// term. Associating compounds: a0 (Pa·m⁶/mol²), b (m³/mol), c1, association energy ε (J/mol), volume β, sites of each kind
// (4C: two donors + two acceptors, 2B: one + one) — the published sets as held in the open NeqSim component database (COMP.csv);
// non-associating compounds use SRK from Tc, Pc and ω. Binary k_ij of the cubic part: NeqSim interaction table (INTER.csv).
export const CPA_PARAMS = Object.freeze({
  W: { name: 'Water', MW: 18.015, Tc: 647.3, Pc: 220.89, w: 0.344, a0: 0.12277, b: 1.4515e-5, c1: 0.67359, eps: 16655, beta: 0.0692, sites: 2 },
  MeOH: { name: 'Methanol', MW: 32.042, Tc: 512.64, Pc: 80.96, w: 0.559, a0: 0.40531, b: 3.0978e-5, c1: 0.43102, eps: 24591, beta: 0.0161, sites: 1 },
  MEG: { name: 'Mono-ethylene glycol', MW: 62.069, Tc: 720, Pc: 90, w: 0.5347, a0: 1.0819, b: 5.14e-5, c1: 0.6744, eps: 19752, beta: 0.0141, sites: 2 },
});
const CPA_KIJ = { 'W-C1': (T) => -0.827413423 + 0.0026055 * T, 'W-N2': (T) => -1.389072518 + 0.003879719 * T, 'W-CO2': (T) => -0.27686 + 0.001121 * T, 'W-C3': () => 0.11, 'W-nC4': () => 0.0875, 'W-H2S': () => 0.1913, 'W-MeOH': () => -0.153, 'W-MEG': () => -0.115, 'MeOH-C1': () => 0.0134, 'MEG-C1': () => 0.134 };
/** CPA Helmholtz model. comps: [{ id, MW, Tc, Pc, w }] with ids 'W', 'MeOH', 'MEG' taken from CPA_PARAMS (associating). */
export function cpaModel(comps) {
  const n = comps.length, cs = comps.map((c) => { const p = CPA_PARAMS[c.id]; return p ? { ...p, id: c.id } : { id: c.id, MW: c.MW, Tc: c.Tc, Pc: c.Pc, w: c.w, a0: (0.42748023 * (R * c.Tc) ** 2) / (c.Pc * 1e5), b: (0.08664035 * R * c.Tc) / (c.Pc * 1e5), c1: 0.48 + 1.574 * c.w - 0.176 * c.w * c.w, sites: 0 }; });
  const kf = cs.map((a) => cs.map((b) => CPA_KIJ[`${a.id}-${b.id}`] || CPA_KIJ[`${b.id}-${a.id}`] || null)), as = cs.map((c, i) => (c.sites ? i : -1)).filter((i) => i >= 0);
  return {
    name: 'CPA', ids: cs.map((c) => c.id), n, M: cs.map((c) => c.MW), crit: cs.map((c) => ({ Tc: c.Tc, Pc: c.Pc, w: c.w })),
    rhoMax(x) { let b = 0; for (let i = 0; i < n; i++) b += x[i] * cs[i].b; return 0.985 / b; },
    ar(T, rho, x) {
      const sa = cs.map((c) => Math.sqrt(c.a0) * (1 + c.c1 * (1 - Math.sqrt(T / c.Tc)))); let a = 0, b = 0;
      for (let i = 0; i < n; i++) { b += x[i] * cs[i].b; for (let j = 0; j < n; j++) a += x[i] * x[j] * sa[i] * sa[j] * (1 - (i !== j && kf[i][j] ? kf[i][j](T) : 0)); }
      let out = -Math.log(1 - b * rho) - (a / (b * R * T)) * Math.log(1 + b * rho);
      const act = as.filter((i) => x[i] > 0); if (!act.length) return out;
      const g = 1 / (1 - 0.475 * b * rho), D = act.map((i) => act.map((j) => g * (Math.exp((cs[i].eps + cs[j].eps) / (2 * R * T)) - 1) * ((cs[i].b + cs[j].b) / 2) * Math.sqrt(cs[i].beta * cs[j].beta)));
      const na = act.length, X = new Array(na).fill(0.2), w = act.map((j) => rho * x[j] * cs[j].sites);
      if (na === 1) { const q = w[0] * D[0][0]; X[0] = (Math.sqrt(1 + 4 * q) - 1) / (2 * q); } // one associating compound: closed form
      else for (let it = 0; it < 500; it++) { // fraction of sites not bonded: damped successive substitution
        let err = 0;
        for (let p = 0; p < na; p++) { let sm = 0; for (let q = 0; q < na; q++) sm += w[q] * X[q] * D[p][q]; const xn = 0.5 * (X[p] + 1 / (1 + sm)); err = Math.max(err, Math.abs(xn - X[p])); X[p] = xn; }
        if (err < 1e-13) break;
      }
      act.forEach((i, p) => (out += x[i] * 2 * cs[i].sites * (Math.log(X[p]) - X[p] / 2 + 0.5)));
      return out;
    },
  };
}
/**
 * Water (and inhibitor) in the gas phase and gas dissolved in the aqueous phase by CPA: the dry gas (EOS fluid f, mole fractions yDry)
 * is equilibrated at P (bara), T (°C) with an aqueous phase of inhibitor mass fraction (wt %) whose water activity is further lowered
 * by salt (awSalt). Returns { yW, yInh (mole fractions in the gas), wc (kg water per Sm³ dry gas), inhLoss (kg per Sm³), xGas (dissolved gas mole fraction), ok }.
 */
export function cpaWater(f, yDry, Pbar, Tc, { inhId = 'none', inhWt = 0, awSalt = 1 } = {}) {
  const TK = Tc + KEL, hasI = !!CPA_PARAMS[inhId] && inhWt > 0, idx = []; yDry.forEach((v, i) => { if (v > 1e-9) idx.push(i); });
  const comps = [{ id: 'W' }, ...(hasI ? [{ id: inhId }] : []), ...idx.map((i) => f.comps[i])], m = cpaModel(comps), n0 = hasI ? 2 : 1, n = comps.length;
  const wI = hasI ? inhWt / 100 : 0, nI = wI / CPA_PARAMS[hasI ? inhId : 'W'].MW, nW = (1 - wI) / 18.015, xI = hasI ? nI / (nI + nW) : 0, tot = idx.reduce((s, i) => s + yDry[i], 0);
  let x = [1 - xI, ...(hasI ? [xI] : []), ...idx.map(() => 0)], y = [1e-3, ...(hasI ? [1e-4] : []), ...idx.map((i) => (yDry[i] / tot) * 0.999)], ok = true, rl = null, rv = null;
  const aqi = hasI ? [0, 1] : [0];
  for (let it = 0; it < 6; it++) { // water and inhibitor between the two phases (the dissolved gas barely changes the aqueous fugacities)
    rl = hRho(m, TK, Pbar, x, 'liquid'); rv = hRho(m, TK, Pbar, y, 'vapour'); if (rl === null || rv === null) { ok = false; break; }
    const pl = hLnPhi(m, TK, rl, x, aqi), pg = hLnPhi(m, TK, rv, y, aqi), yW = x[0] * awSalt * Math.exp(pl[0] - pg[0]), yI = hasI ? x[1] * Math.exp(pl[1] - pg[1]) : 0, dry = Math.max(1 - yW - yI, 1e-6);
    const yn = [yW, ...(hasI ? [yI] : []), ...idx.map((i) => (yDry[i] / tot) * dry)], err = Math.abs(yn[0] - y[0]); y = yn; if (err < 1e-7 * Math.max(y[0], 1e-9)) break;
  }
  if (ok) { // gas dissolved in the aqueous phase from equal fugacities, one pass
    const pl = hLnPhi(m, TK, rl, x), pg = hLnPhi(m, TK, rv, y), xg = idx.map((_, k) => y[n0 + k] * Math.exp(pg[n0 + k] - pl[n0 + k])), sg = sumA(xg);
    x = [(1 - xI) * (1 - sg), ...(hasI ? [xI * (1 - sg)] : []), ...xg];
  }
  const dry = 1 - y[0] - (hasI ? y[1] : 0);
  return { yW: y[0], yInh: hasI ? y[1] : 0, wc: ((y[0] / dry) * 18.015e-3) / VM_STD, inhLoss: hasI ? ((y[1] / dry) * CPA_PARAMS[inhId].MW * 1e-3) / VM_STD : 0, xGas: sumA(x.slice(n0)), xGasBy: Object.fromEntries(idx.map((i, k) => [f.comps[i].id, x[n0 + k]])), ok };
}
// ---- Lee–Kesler (modified Benedict–Webb–Rubin) as a residual-Helmholtz model -----------------------------------------------
// ar(Tr, Vr) of one Lee–Kesler reference fluid (constants c = [b1..b4, c1..c4, d1, d2, β, γ]); Z = 1 − Vr ∂ar/∂Vr reproduces the BWR form.
const lkAr = (c, Tr, V) => { const B = c[0] - c[1] / Tr - c[2] / Tr ** 2 - c[3] / Tr ** 3, C = c[4] - c[5] / Tr + c[6] / Tr ** 3, D = c[8] + c[9] / Tr, u = c[11] / (V * V); return B / V + C / (2 * V * V) + D / (5 * V ** 5) + (c[7] / (2 * Tr ** 3 * c[11])) * (c[10] + 1 - (c[10] + 1 + u) * Math.exp(-u)); };
const lkState = (c, Tr, Pr, liquid) => { const V = lkVr(c, Tr, Pr, liquid); if (V === null) return null; const Z = (Pr * V) / Tr; return { Z, V, hDep: Tr * (Z - 1 - Tr * d1((t) => lkAr(c, t, V), Tr, Tr * 1e-3)), lnPhi: lkAr(c, Tr, V) + Z - 1 - Math.log(Z) }; };
/**
 * Lee–Kesler (1975) three-parameter corresponding states from the modified BWR equation of the simple and reference (n-octane)
 * fluids: { Z, hDep ((H − H_ideal)/(R·Tc)), lnPhi (ln of the fugacity coefficient) } at reduced T, P and acentric factor w.
 * liquid = true takes the dense root of both fluids.
 */
export function leeKesler(Tr, Pr, w = 0, liquid = false) {
  const a = lkState(LK0, Tr, Pr, liquid), b = lkState(LKR, Tr, Pr, liquid); if (!a || !b) return null;
  const q = w / 0.3978, mix = (k) => a[k] + q * (b[k] - a[k]);
  return { Z: mix('Z'), hDep: mix('hDep'), lnPhi: mix('lnPhi'), Z0: a.Z, Zr: b.Z };
}

// ---- Pedersen corresponding-states viscosity (Pedersen et al. 1984; Pedersen & Fredenslund 1987) --------------------------
// Methane is the reference fluid: viscosity correlation of Hanley, McCarty & Haynes (1975) with the low-temperature branch of
// Pedersen & Fredenslund (1987) — coefficients as in the open NeqSim implementation (PFCTViscosityMethodMod86.java) — and the
// methane density from GERG-2008 instead of the 1974 BWR equation of McCarty.
const HMH = { GV: [-2.090975e5, 2.647269e5, -1.472818e5, 4.71674e4, -9.491872e3, 1.219979e3, -9.627993e1, 4.274152, -8.141531e-2], A: 1.696985927, B: -0.133372346, C: 1.4, F: 168, j: [-10.35060586, 17.571599671, -3019.3918656, 188.73011594, 0.042903609488, 145.29023444, 6127.6818706], k: [-9.74602, 18.0834, -4126.66, 44.6055, 0.976544, 81.8134, 15649.9], rhoc: 0.16266 };
/** Viscosity of methane (Pa·s) at T (K) and density rho (kg/m³): Hanley et al. (1975), dense-liquid branch of Pedersen & Fredenslund (1987) below 91 K. */
export function methaneViscosity(T, rho) {
  const r = rho / 1000, th = (r - HMH.rhoc) / HMH.rhoc; let e0 = 0; for (let i = 0; i < 9; i++) e0 += HMH.GV[i] * T ** ((i - 3) / 3);
  const e1 = (HMH.A + HMH.B * (HMH.C - Math.log(T / HMH.F)) ** 2) * r, br = (c) => Math.exp(c[0] + c[3] / T) * (Math.exp(r ** 0.1 * (c[1] + c[2] / T ** 1.5) + th * Math.sqrt(r) * (c[4] + c[5] / T + c[6] / (T * T))) - 1);
  const ht = Math.tanh(T - 91);
  return (e0 + e1 + ((ht + 1) / 2) * br(HMH.j) + ((1 - ht) / 2) * br(HMH.k)) * 1e-7;
}
let gergC1 = null;
const methaneRho = (T, Pbar) => { gergC1 ||= gergModel(['C1']); const Tq = Math.max(T, 40), r = hRho(gergC1, Tq, Pbar, [1], 'stable'); return r === null ? null : r * 16.04246e-3; }; // kg/m³
/**
 * Pedersen corresponding-states viscosity (Pa·s) of a mixture of composition x at P (bara), T (K). comps: [{ Tc, Pc, MW }].
 * Returns { mu, T0, P0 (reference-state methane conditions), Tcm, Pcm, Mmix, alpha } or null outside the reference equation.
 */
export function viscosityPedersen(comps, x, Pbar, TK) {
  const n = comps.length, Tc0 = 190.564, Pc0 = 45.992, M0 = 16.043; let t1 = 0, t2 = 0, mw = 0, mn = 0;
  for (let i = 0; i < n; i++) { if (!(x[i] > 0)) continue; mn += x[i] * comps[i].MW; mw += x[i] * comps[i].MW ** 2; const ci = Math.cbrt(comps[i].Tc / comps[i].Pc); for (let j = 0; j < n; j++) { if (!(x[j] > 0)) continue; const v = x[i] * x[j] * (ci + Math.cbrt(comps[j].Tc / comps[j].Pc)) ** 3; t1 += v * Math.sqrt(comps[i].Tc * comps[j].Tc); t2 += v; } }
  if (!(t2 > 0)) return null;
  const Tcm = t1 / t2, Pcm = (8 * t1) / (t2 * t2), Mmix = mn + 1.304e-4 * ((mw / mn) ** 2.303 - mn ** 2.303);
  const r0 = methaneRho((TK * Tc0) / Tcm, (Pbar * Pc0) / Pcm); if (r0 === null) return null;
  const rr = r0 / (10.15 * 16.043), al = 1 + 7.378e-3 * rr ** 1.847 * Mmix ** 0.5173, al0 = 1 + 7.378e-3 * rr ** 1.847 * M0 ** 0.5173; // reduced density with ρc = 10.15 mol/L
  const T0 = ((TK * Tc0) / Tcm) * (al0 / al), P0 = ((Pbar * Pc0) / Pcm) * (al0 / al), rho0 = methaneRho(T0, P0); if (rho0 === null) return null;
  const mu = (Tcm / Tc0) ** (-1 / 6) * (Pcm / Pc0) ** (2 / 3) * Math.sqrt(Mmix / M0) * (al / al0) * methaneViscosity(Math.max(T0, 40), rho0);
  return Number.isFinite(mu) && mu > 0 ? { mu, T0, P0, Tcm, Pcm, Mmix, alpha: al } : null;
}

// ---- gas hydrates: Kihara cell potential (Lennard-Jones–Devonshire smoothed cell) and structure H --------------------------
// Kihara parameters of the guest–water interaction [hard-core radius a (Å), collision diameter σ (Å), ε/k (K)]. C1–iC4: the set
// optimised by Avaji et al. (2023, Fluid Phase Equilib. 567, 113716; open manuscript, Table 8) together with the Sloan & Koh reference
// properties listed there; N2, CO2, nC4 and H2S: hydrate columns of the open NeqSim component database (not regressed with the
// same reference properties — less accurate). Cavity radii (Å) and coordination numbers: Sloan & Koh as tabulated by Herri et al.
export const KIHARA = Object.freeze({
  kihara: true,
  guests: { C1: [0.383, 3.1436, 155.8], C2: [0.59, 3.2998, 178.709], C3: [0.647, 3.419, 191.855], iC4: [0.8921, 3.20691, 198.332], N2: [0.359188902, 3.132506748, 126.578386713], CO2: [0.68463388, 3.03720716, 170.162382832], nC4: [0.9379, 2.9125, 209], H2S: [0.36, 3.1688, 206.61] },
  cav: { sI: [[3.95, 20], [4.33, 24]], sII: [[3.91, 20], [4.73, 28]], sH: [[3.91, 20], [4.06, 20]] },
  ref: { sI: { dmu0: 1263.6, dhL: -4858.9, dhI: 1151, dvL: 4.6e-6, dvI: 3.0e-6 }, sII: { dmu0: 882.8, dhL: -5202.2, dhI: 808, dvL: 5.0e-6, dvI: 3.4e-6 } },
  dcp: [-38.12, 0.141],
  munck: { nC4: ['sI', 'sII'], CO2: ['sII'] }, // guest → structures for which the Munck Langmuir constants are kept inside the Kihara set
});
const kihCache = new Map();
/** Forget the tabulated Kihara Langmuir constants (needed only after KIHARA.guests has been changed, e.g. by a regression). */
export const kiharaReset = () => kihCache.clear();
/**
 * Langmuir constant (1/bar) from the Kihara potential averaged over a spherical cell of radius Rc (Å) with z water molecules
 * (McKoy & Sinanoğlu form of the Lennard-Jones–Devonshire theory): C = (4π/kT) ∫₀^{Rc−a} exp(−w(r)/kT) r² dr.
 */
export function kiharaC(a, sig, epsK, Rc, z, TK) {
  const dl = (N, r) => ((1 - r / Rc - a / Rc) ** -N - (1 + r / Rc - a / Rc) ** -N) / N;
  const w = (r) => 2 * z * epsK * ((sig ** 12 / (Rc ** 11 * r)) * (dl(10, r) + (a / Rc) * dl(11, r)) - (sig ** 6 / (Rc ** 5 * r)) * (dl(4, r) + (a / Rc) * dl(5, r))); // in K
  const n = 100, top = (Rc - a) * 0.999, h = top / n; let s = 0;
  for (let i = 1; i <= n; i++) { const r = i * h, f = Math.exp(-w(r) / TK) * r * r; s += (i === n ? 1 : i % 2 ? 4 : 2) * (Number.isFinite(f) ? f : 0); }
  return ((4 * Math.PI) / (1.380649e-23 * TK)) * ((s * h) / 3) * 1e-30 * 1e5; // Å³ → m³, 1/Pa → 1/bar
}
// ln C is tabulated every 5 K (220–340 K, nodes filled on demand) for each guest and cavity and interpolated linearly in 1/T.
function kiharaLangmuir(id, s, m, TK) {
  const key = `${id}|${s}|${m}`; let t = kihCache.get(key); if (!t) { t = new Array(25).fill(null); kihCache.set(key, t); }
  const T = clamp(TK, 220, 339.999), i = Math.floor((T - 220) / 5), Ta = 220 + 5 * i, u = (1 / T - 1 / Ta) / (1 / (Ta + 5) - 1 / Ta), k = KIHARA.guests[id], cv = KIHARA.cav[s][m];
  for (const j of [i, i + 1]) if (t[j] === null) t[j] = Math.log(Math.max(kiharaC(k[0], k[1], k[2], cv[0], cv[1], 220 + 5 * j), 1e-300));
  return Math.exp(t[i] + u * (t[i + 1] - t[i]));
}
// ---- solid phases: asphaltene (Flory–Huggins / Hirschberg), water-side Gibbs minimum, phase inventory ------------------------
/** Hildebrand solubility parameter (MPa^0.5) of a liquid from the EOS: square root of the cohesive energy density −U_res/v. */
export function solubilityParameter(f, x, Pbar, Tc) {
  const TK = Tc + KEL, q = phaseHS(f, x, Pbar, TK, 'liquid'), p = phaseProps(f, x, Pbar, Tc, 'liquid', { thermal: false });
  return Math.sqrt(Math.max(-(q.hRes - R * TK * (q.Z - 1)) / p.vm, 0)) / 1000;
}
/**
 * Asphaltene solid–liquid equilibrium by the Flory–Huggins regular-solution model of Hirschberg et al. (1984): the largest
 * asphaltene volume fraction the liquid can hold solves ln φ + (1 − Va/VL)(1 − φ) + χ(1 − φ)² = 0 with χ = Va(δa − δL)²/RT.
 * Inputs: deltaL (MPa^0.5), vL (m³/mol), TK, wAsp (mass fraction of asphaltene in the liquid), rhoL (kg/m³), aspMW (g/mol), aspRho (kg/m³),
 * aspDelta (MPa^0.5 at 25 °C, falling by the fraction aspDeltaT per K). Returns { chi, phiMax, phiA, precipWt (wt % of the liquid precipitated), stable }.
 */
export function asphalteneFH({ deltaL, vL, TK, wAsp, rhoL, aspMW = 750, aspRho = 1100, aspDelta = 21, aspDeltaT = 1.07e-3 }) {
  const Va = (aspMW * 1e-3) / aspRho, dA = aspDelta * (1 - aspDeltaT * (TK - 298.15)), chi = (Va * ((dA - deltaL) * 1e3) ** 2) / (R * TK), r = Va / vL, g = (lp) => { const p = Math.exp(lp); return lp + (1 - r) * (1 - p) + chi * (1 - p) ** 2; };
  let phiMax = 1; const grid = linspace(Math.log(1e-14), Math.log(0.999), 60);
  for (let k = 1; k < grid.length; k++) if (g(grid[k - 1]) < 0 && g(grid[k]) >= 0) { phiMax = Math.exp(brent(g, grid[k - 1], grid[k], 1e-12, 80)); break; }
  const phiA = clamp((wAsp * rhoL) / aspRho, 0, 1), ex = Math.max(0, phiA - phiMax);
  return { chi, phiMax, phiA, precipWt: (100 * ex * aspRho) / rhoL, stable: !(ex > 0) };
}
/**
 * Asphaltene precipitation along an isothermal depletion: at each pressure the EOS flash gives the liquid, its solubility parameter and
 * molar volume, and the Flory–Huggins equilibrium gives the precipitated amount. wAspSto: asphaltene mass fraction of the stock-tank oil.
 * Returns { rows: [{ P, deltaL, phiMax, phiA, precipWt }], upperOnset, lowerOnset (bara | null), maxPrecip (wt % of the liquid), pAtMax }.
 */
export function asphalteneCurve(f, Tc, pressures, { wAspSto = 0.025, aspMW = 750, aspRho = 1100, aspDelta = 21, aspDeltaT = 1.07e-3 } = {}) {
  const TK = Tc + KEL, at = (P) => {
    const s = props(f, P, Tc, { thermal: false }); if (s.phase === 'gas') return null;
    const x = s.x, sf = stdFlash(withZ(f, x)), mSto = sf.vOil > 0 ? ((1 - sf.beta) * (sf.mwOil || s.oil.MW)) / s.oil.MW : 1; // mass of stock-tank oil per mass of live liquid
    const dL = solubilityParameter(f, x, P, Tc), r = asphalteneFH({ deltaL: dL, vL: s.oil.vm, TK, wAsp: wAspSto * clamp(mSto, 0, 1), rhoL: s.oil.rho, aspMW, aspRho, aspDelta, aspDeltaT });
    return { P, deltaL: dL, phiMax: r.phiMax, phiA: r.phiA, precipWt: r.precipWt, margin: r.phiA - r.phiMax };
  };
  const rows = pressures.map(at).filter(Boolean).sort((a, b) => b.P - a.P); let upperOnset = null, lowerOnset = null, maxPrecip = 0, pAtMax = null;
  const cross = (a, b) => { try { return brent((P) => at(P)?.margin ?? -1, Math.min(a, b), Math.max(a, b), 1e-3, 40); } catch { return 0.5 * (a + b); } };
  for (let k = 0; k < rows.length; k++) {
    if (rows[k].precipWt > maxPrecip) { maxPrecip = rows[k].precipWt; pAtMax = rows[k].P; }
    if (k > 0 && rows[k - 1].margin <= 0 && rows[k].margin > 0 && upperOnset === null) upperOnset = cross(rows[k - 1].P, rows[k].P);
    if (k > 0 && rows[k - 1].margin > 0 && rows[k].margin <= 0) lowerOnset = cross(rows[k - 1].P, rows[k].P);
  }
  if (rows.length && rows[0].margin > 0) upperOnset = rows[0].P;
  return { rows, upperOnset, lowerOnset, maxPrecip, pAtMax };
}
/**
 * Water-side Gibbs-energy minimum at P (bara), T (°C): free water of an aqueous phase aq ({ S g/kg, inhId, inhWt }) converts to hydrate
 * (structure of highest driving force) until the chemical potential of water in the concentrating brine / inhibitor solution equals that
 * in the hydrate, the water is used up, or the gas runs out (gasPerWater: mol of hydrate formers per mol of water). Returns
 * { phases: [names], structure, conversion (fraction of the water in hydrate), drive0 (Δμ/RT before conversion; > 0 hydrate stable), ice,
 *   dG (J per mol of water, ≤ 0), awFinal, saltFinal (g/kg), inhFinal (wt %), hydrateKgPerKgWater, hydrationNumber, limitedBy }.
 */
export function waterPhaseEquilibrium(f, Pbar, Tc, aq = {}, models = {}, { table = HYD_ACTIVE.table || KIHARA, sH = HYD_ACTIVE.sH, gasPerWater = Infinity, allowHydrate = true, tune = null } = {}) {
  const TK = Tc + KEL, fug = formerFugacities(f, Pbar, Tc), S0 = clamp(aq.S || 0, 0, 260) / 1000, w0 = aq.inhId && aq.inhId !== 'none' ? clamp(aq.inhWt || 0, 0, 95) / 100 : 0;
  const structs = ['sI', 'sII', ...(sH ? ['sH'] : [])], mW0 = (1 - w0) * (1 - S0), mS = (1 - w0) * S0; // per kg of aqueous phase
  const conc = (xi) => { const mW = mW0 * (1 - xi), tot = mW + mS + w0; return { S: clamp((1000 * mS) / Math.max(mW + mS, 1e-12), 0, 260), inhId: aq.inhId || 'none', inhWt: clamp((100 * w0) / Math.max(tot, 1e-12), 0, 95) }; };
  const best = (xi) => { const c = conc(xi), wa = waterActivity(c, TK, models, tune), a = wa.aw * (1 - gasSolubility(fug, Pbar, Tc, wa.molality).total); let b = null; for (const s of structs) { const st = hydrateState(s, TK, Pbar, fug, a, table); if (!b || st.drive > b.st.drive) b = { s, st }; } return { ...b, c, aw: a }; };
  const b0 = best(0), base = { structure: null, conversion: 0, drive0: b0.st.drive, ice: b0.st.ice, dG: 0, awFinal: b0.aw, saltFinal: b0.c.S, inhFinal: b0.c.inhWt, hydrateKgPerKgWater: 0, hydrationNumber: null, limitedBy: null };
  if (!Object.keys(fug).length || !(b0.st.drive > 0) || !allowHydrate) return { ...base, phases: [b0.st.ice ? 'ice' : 'aqueous liquid'], limitedBy: !allowHydrate && b0.st.drive > 0 ? 'hydrate not permitted (metastable aqueous phase)' : null };
  const occ = occupancy(b0.s, b0.st), nH = occ.hydrationNumber || HYDRATE_STRUCTURES[b0.s].waters / 4, xiGas = Math.min(1, gasPerWater * nH);
  let xi = Math.min(0.9999, xiGas), limitedBy = xiGas < 1 ? 'hydrate formers' : 'water';
  if ((mS > 0 || w0 > 0) && best(xi).st.drive > 0 === false) { xi = brent((q) => best(q).st.drive, 0, xi, 1e-10, 80); limitedBy = 'water activity (salt / inhibitor concentrating)'; }
  else if (mS <= 0 && w0 <= 0 && xiGas >= 1) xi = 1;
  const n = 12; let integ = 0; for (let k = 0; k <= n; k++) { const q = (xi * k) / n; integ += (k === 0 || k === n ? 0.5 : 1) * best(Math.min(q, 0.9999)).st.drive; } // ∫ drive dξ
  const e = best(Math.min(xi, 0.9999)), mwGuest = Object.entries(occ.byGuest.reduce((o, t, m) => { for (const id in t) o[id] = (o[id] || 0) + HYDRATE_STRUCTURES[b0.s].nu[m] * t[id]; return o; }, {})).reduce((s, [id, v]) => s + v * (COMPONENTS[id]?.MW || 16), 0);
  return { ...base, phases: xi >= 1 ? ['hydrate ' + b0.s] : ['hydrate ' + b0.s, 'aqueous liquid'], structure: b0.s, conversion: xi, dG: -R * TK * ((integ * xi) / n), awFinal: e.aw, saltFinal: e.c.S, inhFinal: e.c.inhWt, hydrateKgPerKgWater: xi * (1 + mwGuest / MW_W), hydrationNumber: nH, limitedBy };
}
/** Reduced molar Gibbs energies g/RT of the feed as one phase and as the equilibrium split at P (bara), T (°C): { g1, g2, dg = g2 − g1 (≤ 0), tpd }. */
export function gibbsSplit(f, Pbar, Tc) {
  const TK = Tc + KEL, one = eosPhase(f, f.z, Pbar, TK), fl = flashPT(f, Pbar, Tc), st = f.n > 1 ? stability(f, Pbar, TK) : { tpd: 0 };
  const g2 = fl.phase === 'two' ? fl.beta * eosPhase(f, fl.y, Pbar, TK, 'vapour').g + (1 - fl.beta) * eosPhase(f, fl.x, Pbar, TK, 'liquid').g : one.g;
  return { g1: one.g, g2, dg: g2 - one.g, tpd: st.tpd };
}

// ---- suite: fluid library, model assembly and reporting ------------------------------------------------------------------
const LIB = (name, kind, comp, c7MW, c7SG, extra = {}) => Object.freeze({ name, kind, spec: Object.freeze({ name, comp: Object.freeze(comp), c7MW, c7SG, ...extra }) });
/** Built-in fluids that the suite can analyse instead of the case fluid (mol %, water-free). */
export const FLUID_LIBRARY = Object.freeze({
  blackOil: LIB('Black oil (GOR ≈ 100)', 'oil', { N2: 0.3, CO2: 0.9, H2S: 0, C1: 32, C2: 5.5, C3: 4.5, iC4: 1, nC4: 2.6, iC5: 1.3, nC5: 1.8, C6: 3.1, C7p: 47 }, 245, 0.87),
  volatileOil: LIB('Volatile oil', 'oil', { N2: 0.5, CO2: 1.8, H2S: 0, C1: 60, C2: 8.5, C3: 5, iC4: 1.2, nC4: 2.3, iC5: 1, nC5: 1.2, C6: 1.8, C7p: 16.7 }, 185, 0.82),
  gasCondensate: LIB('Gas condensate', 'gas', { N2: 0.6, CO2: 2.4, H2S: 0, C1: 73, C2: 7.8, C3: 3.6, iC4: 0.7, nC4: 1.4, iC5: 0.6, nC5: 0.6, C6: 1.1, C7p: 8.2 }, 150, 0.79, { rateBasis: 'gas', qGas: 3, qWater: 40 }),
  leanGas: LIB('Lean gas', 'gas', { N2: 1, CO2: 1.5, H2S: 0, C1: 90, C2: 4.5, C3: 1.5, iC4: 0.3, nC4: 0.4, iC5: 0.15, nC5: 0.15, C6: 0.2, C7p: 0.3 }, 110, 0.75, { rateBasis: 'gas', qGas: 5, qWater: 15 }),
  sourGas: LIB('Sour CO₂-rich gas', 'gas', { N2: 1, CO2: 12, H2S: 6, C1: 70, C2: 5, C3: 2.5, iC4: 0.5, nC4: 0.8, iC5: 0.3, nC5: 0.3, C6: 0.4, C7p: 1.2 }, 120, 0.76, { rateBasis: 'gas', qGas: 4, qWater: 25 }),
});
const mergeSpec = (fluid) => ({ ...DEFAULT_FLUID, ...(fluid || {}), comp: { ...DEFAULT_FLUID.comp, ...((fluid || {}).comp || {}) } });
function resolveSpec(fluid, source) {
  const base = mergeSpec(fluid), lib = FLUID_LIBRARY[source];
  const spec = lib ? { ...base, ...lib.spec, comp: { ...lib.spec.comp } } : base;
  for (const k of COMP_IDS) { const x = +spec.comp[k]; if (!Number.isFinite(x) || x < 0) throw new Error(`The mole percentage of ${COMP_LABELS[k]} must be a number of zero or more.`); }
  if (!(COMP_IDS.reduce((s, k) => s + +spec.comp[k], 0) > 0)) throw new Error('The fluid composition is empty: enter at least one component with a positive mole percentage on the case page.');
  return { spec, library: lib || null };
}
const tuneOf = (v, spec) => ({
  eos: EOS[v.eosSel] ? v.eosSel : spec.eos, nPseudo: ['1', '2', '3'].includes(String(v.nPseudoSel)) ? +v.nPseudoSel : spec.nPseudo, method: ['KL', 'RD', 'Twu'].includes(v.critCorr) ? v.critCorr : 'KL',
  tcMult: num(v.tcMult, 1), pcMult: num(v.pcMult, 1), wMult: num(v.wMult, 1), kijScale: num(v.kijScale, 1), vcMult: num(v.vcMult, 1), shift: v.shift !== false,
});
const ION_DATA = { Na: [22.99, 1], K: [39.098, 1], Ca: [40.078, 2], Mg: [24.305, 2], Ba: [137.33, 2], Sr: [87.62, 2], Fe: [55.845, 2], Cl: [35.453, 1], SO4: [96.06, 2], HCO3: [61.017, 1], Br: [79.904, 1] };
/**
 * Ionic strength of a brine analysis and its NaCl equivalent. ions: [{ ion, mgL }]. Returns { I (mol/kg water), tds (g/L), S (g NaCl-equivalent
 * per kg brine at the same ionic strength), balance (cation − anion equivalents over their sum) } or null when no ion is given.
 */
export function brineFromIons(ions) {
  let tds = 0, I = 0, cat = 0, an = 0;
  for (const r of Array.isArray(ions) ? ions : []) { const d = ION_DATA[String(r.ion || '').trim()], c = +r.mgL; if (!d || !(c > 0)) continue; const mol = c / 1000 / d[0]; tds += c / 1000; I += 0.5 * mol * d[1] ** 2; if (['Cl', 'SO4', 'HCO3', 'Br'].includes(String(r.ion).trim())) an += mol * d[1]; else cat += mol * d[1]; }
  if (!(tds > 0)) return null;
  const kgW = Math.max(0.5, 1 + 0.0007 * tds - tds / 1000), m = I / kgW; // kg of water per litre of brine
  return { I: m, tds, S: (1000 * m * MW_NACL) / (1000 + m * MW_NACL), balance: cat + an > 0 ? (cat - an) / (cat + an) : 0 };
}
// User binary-interaction overrides: rows { a, b, kij } with component ids ('C7+' addresses every pseudo-component).
function applyKij(f, rows) {
  let n = 0; const match = (c, id) => c.id === id || (c.pseudo && (id === 'C7+' || id === 'C7p'));
  for (const r of Array.isArray(rows) ? rows : []) {
    const a = String(r.a || '').trim(), b = String(r.b || '').trim(), k = +r.kij;
    if (!a || !b || !Number.isFinite(k)) continue;
    if (Math.abs(k) > 0.5) throw new Error(`The binary interaction parameter ${a}–${b} = ${k} is outside the plausible range of −0.5 to 0.5.`);
    f.comps.forEach((ci, i) => f.comps.forEach((cj, j) => { if (i !== j && ((match(ci, a) && match(cj, b)) || (match(ci, b) && match(cj, a)))) { f.kij[i][j] = k; n++; } }));
  }
  return n / 2;
}
const TABLE_FIELDS = ['wG', 'rhoG', 'rhoO', 'muG', 'muO', 'cpG', 'cpO', 'kG', 'kO', 'hG', 'hO', 'jtG', 'jtO', 'sigma', 'zG', 'mwG', 'mwO'];
// Refill a kernel property table with a fluid whose interaction parameters differ from the kernel defaults (same layout and stand-in rules).
function refillTable(t, spec, f) {
  t.P.forEach((p, i) => t.T.forEach((tc, j) => {
    const s = props(f, p, tc), o = { wG: s.wG, rhoG: s.gas.rho, rhoO: s.oil.rho, muG: s.gas.mu, muO: s.oil.mu, cpG: s.gas.cp, cpO: s.oil.cp, kG: s.gas.k, kO: s.oil.k, hG: s.gas.h, hO: s.oil.h, jtG: s.gas.jt, jtO: s.oil.jt, sigma: s.sigma, zG: s.gas.Z, mwG: s.gas.MW, mwO: s.oil.MW };
    if (o.wG >= 1) { o.rhoO = Math.max(o.rhoO, 500); o.muO = Math.max(o.muO, 2e-4); }
    if (o.wG <= 0) { o.rhoG = Math.min(o.rhoG, (p * 1e5 * 0.02) / (0.9 * R * (tc + KEL))); o.muG = Math.min(o.muG, 1.5e-5); o.zG = 0.9; o.mwG = 20; o.cpG = 2300; o.kG = 0.035; o.jtG = 3e-6; }
    for (const k of TABLE_FIELDS) t[k][i][j] = o[k];
  }));
  const st = streams(spec, f);
  t.gasSG = st.std.gasSG ?? f.MW / MW_AIR; t.rates = { nHC: st.nHC, mHC: st.mHC, mW: st.mW, qOilStd: st.qOilStd, qGasStd: st.qGasStd, qWaterStd: st.qWaterStd, gor: Number.isFinite(st.gor) ? st.gor : null, wc: st.wc, api: st.std.api, rhoOilStd: st.std.rhoOilStd, mwGas: st.std.mwGas };
  return t;
}
/** Tuned EOS fluid of a specification: characterisation plus optional binary-interaction overrides. Returns { ch, f, nKij }. */
export function tunedFluid(spec, tune = {}, kijRows = [], tx = TX0) {
  const ch = characterise(spec, tune), nKij = applyKij(ch.fluid, kijRows), f = ch.fluid;
  if (tx.shift !== 1 || tx.par !== 1 || tx.cp !== 1) for (const c of f.comps) { c.c *= tx.shift; c.par *= tx.par; c.cp = c.cp.map((q) => q * tx.cp); } // calibration multipliers: volume shift, parachor, ideal-gas heat capacity
  return { ch, f, nKij };
}
// Calibration multipliers beyond the C7+ characterisation: volume shift, parachor, ideal-gas cp, gas viscosity, water density, water content of gas,
// hydrate reference chemical potential, inhibitor activity (ln γ of water) and salt (ln a_w of the brine).
const TX0 = Object.freeze({ shift: 1, par: 1, cp: 1, muG: 1, rhoW: 1, wc: 1, hyd: 1, act: 1, salt: 1 });
const tuneX = (v) => ({ shift: clamp(num(v.shiftMult, 1), 0, 3), par: clamp(num(v.parMult, 1), 0.5, 1.5), cp: clamp(num(v.cpMult, 1), 0.7, 1.3), muG: clamp(num(v.muGMult, 1), 0.5, 2), rhoW: clamp(num(v.rhoWMult, 1), 0.9, 1.1), wc: clamp(num(v.wcMult, 1), 0.3, 3), hyd: clamp(num(v.hydMult, 1), 0.8, 1.2), act: clamp(num(v.actMult, 1), 0, 3), salt: clamp(num(v.saltMult, 1), 0, 3) });
const txFluid = (tx) => tx.shift !== 1 || tx.par !== 1 || tx.cp !== 1;
/** Hydrate constants set of a model selection ('vdwp' Munck, 'kihara') with the reference chemical potential Δμ⁰ scaled by `mult`. */
export function hydrateTable(model = 'kihara', mult = 1) {
  const base = model === 'vdwp' ? LANGMUIR : KIHARA; if (mult === 1) return base;
  const ref = {}; for (const st of ['sI', 'sII', 'sH']) { const S = { ...HYDRATE_STRUCTURES[st], ...(base.ref?.[st] || {}) }; ref[st] = { ...S, dmu0: S.dmu0 * mult }; }
  return { ...base, ref };
}
/**
 * Multiple-contact miscibility pressure (bara) of an injection gas with the oil at T (°C): the lowest pressure at which repeated
 * contacts develop a single phase — forward (equilibrium gas against fresh oil: vaporising drive) or backward (equilibrium liquid against
 * fresh gas: condensing drive); the first contact is the first-contact miscibility test. A single-cell estimate of the slim-tube
 * minimum miscibility pressure (combined condensing/vaporising drives can be miscible at a somewhat lower pressure).
 * Returns { mmp, mechanism: 'first contact' | 'vaporising' | 'condensing' | null, capped }.
 */
export function miscibilityPressure(f, gasComp, Tc, { Pmin = 20, Pmax = 1200, contacts = 6 } = {}) {
  const gt = Object.values(gasComp).reduce((a, b) => a + b, 0), zg = f.comps.map((c) => (c.pseudo ? 0 : (+gasComp[c.id] || 0) / gt)), zo = f.z;
  if (!(sumA(zg) > 0.5)) return { mmp: null, mechanism: null, capped: false };
  const mixf = (a, b) => { const z = a.map((q, i) => 0.5 * (q + b[i])), fz = withZ(f, z), fl = flashPT(fz, P_, Tc); return fl.phase === 'two' ? fl : null; };
  let P_ = Pmin;
  const test = (P) => { P_ = P; let g = zg, l = zo;
    for (let k = 0; k < contacts; k++) { const a = mixf(g, zo); if (!a) return k === 0 ? 'first contact' : 'vaporising'; g = a.y; }
    for (let k = 0; k < contacts; k++) { const a = mixf(l, zg); if (!a) return k === 0 ? 'first contact' : 'condensing'; l = a.x; }
    return null; };
  let mech = test(Pmax); if (mech === null) return { mmp: Pmax, mechanism: null, capped: true };
  const low = test(Pmin); if (low !== null) return { mmp: Pmin, mechanism: low, capped: true };
  let lo = Pmin, hi = Pmax; for (let k = 0; k < 7; k++) { const m = Math.sqrt(lo * hi), q = test(m); if (q === null) lo = m; else { hi = m; mech = q; } } // resolves the pressure to about 3 %
  return { mmp: hi, mechanism: mech, capped: false };
}
function gergFor(f, x) { // GERG-2008 model of a phase: C7+ up to 1 mol % is counted as n-hexane; null when a larger heavy end or an uncovered component is present
  let ps = 0, tot = 0; const ids = [], xs = [];
  f.comps.forEach((c, i) => { if (!(x[i] > 0)) return; tot += x[i]; if (c.pseudo) ps += x[i]; const id = c.pseudo ? 'C6' : c.id, k = ids.indexOf(id); if (k < 0) { ids.push(id); xs.push(x[i]); } else xs[k] += x[i]; });
  if (ps > 0.01 * tot) return null;
  const m = gergModel(ids); return m ? { m, x: xs.map((q) => q / tot), lumped: ps / tot } : null;
}
/**
 * One phase of composition x at P (bara), T (°C) by the alternative models: PC-SAFT and GERG-2008 (density kg/m³, Z, cp J/kg/K, speed of
 * sound m/s, Joule–Thomson K/bar, residual enthalpy J/mol), Lee–Kesler with Kay's rule (vapour only) and the Pedersen viscosity (Pa·s).
 */
export function altPhase(f, x, Pbar, Tc, kind, { visc = true, thermal = true } = {}) {
  const TK = Tc + KEL, cp0 = cpIg(f, x, TK), out = { pcsaft: null, gerg: null, lk: null, muPedersen: null };
  const run = (m, xx) => { const r = hRho(m, TK, Pbar, xx, kind === 'vapour' ? 'vapour' : 'liquid'); if (r === null) return null; if (!thermal) { let M = 0; xx.forEach((q, i) => (M += q * m.M[i])); return { rho: r * M * 1e-3, Z: (Pbar * 1e5) / (r * R * TK) }; } const q = hProps(m, TK, r, xx, cp0); return Number.isFinite(q.w) ? { rho: q.rhoMass, Z: q.Z, cp: q.cp / (q.M * 1e-3), w: q.w, jt: q.jt, hRes: q.hRes } : null; };
  try { const pm = pcsaftModel(f.comps); if (pm) out.pcsaft = run(pm, x); } catch { /* outside the model */ }
  try { const g = gergFor(f, x); if (g) { const r = run(g.m, g.x); if (r) out.gerg = { ...r, lumped: g.lumped }; } } catch { /* outside the model */ }
  if (kind === 'vapour') { const k = kay(f, x), l = leeKesler(TK / k.Tc, Pbar / k.Pc, k.w); if (l && l.Z > 0) out.lk = { Z: l.Z, rho: (Pbar * 1e5 * k.M * 1e-3) / (l.Z * R * TK), hRes: l.hDep * R * k.Tc }; }
  if (visc) { const q = viscosityPedersen(f.comps, x, Pbar, TK); out.muPedersen = q ? q.mu : null; }
  return out;
}
let lastModel = null; // spec and kernel options of the most recent run: used by the calibration model and the flash calculator
const seedFor = (env, T) => { let best = null, d = Infinity; for (const s of env?.seeds || []) { const q = Math.abs(s.T - T); if (q < d) { d = q; best = s.X; } } return d < 40 ? best : null; };
/** Complete state at P (bara), T (°C) with derivative properties of each phase present: { phase, beta, wG, sigma, x, y, gas | null, oil | null }. */
export function fluidState(f, Pbar, Tc) {
  const s = props(f, Pbar, Tc, { thermal: false });
  return { phase: s.phase, beta: s.beta, wG: s.wG, sigma: s.phase === 'two' ? s.sigma : null, x: s.x, y: s.y, gas: s.phase !== 'oil' ? derivedProps(f, s.y, Pbar, Tc, 'vapour') : null, oil: s.phase !== 'gas' ? derivedProps(f, s.x, Pbar, Tc, 'liquid') : null };
}
const phaseName = { two: 'gas + oil', gas: 'single-phase gas', oil: 'single-phase liquid' };
const aad = (meas, pred) => { let s = 0, n = 0; meas.forEach((m, i) => { if (isNum(m) && isNum(pred[i]) && m !== 0) { s += Math.abs(pred[i] / m - 1); n++; } }); return n ? (100 * s) / n : null; };

// Synthetic "laboratory" report of the reference fluid: generated once from the model with C7+ Tc × 1.025, Pc × 0.97, kij × 1.15 and
// viscosity Vc × 1.06, plus 0.3–3 % noise.
const LAB = {
  psat: [{ t: 40, p: 193 }, { t: 60, p: 207.6 }, { t: 90, p: 224.5 }, { t: 120, p: 238.1 }],
  cce: [{ p: 300, v: 0.98159 }, { p: 260, v: 0.988 }, { p: 224.5, v: 1 }, { p: 180, v: 1.0915 }, { p: 140, v: 1.2433 }, { p: 100, v: 1.5594 }, { p: 60, v: 2.3802 }],
  dle: [{ p: 224.5, rs: 179.3, bo: 1.566, rho: 645.2 }, { p: 180, rs: 140.3, bo: 1.457, rho: 665.6 }, { p: 140, rs: 107.9, bo: 1.381, rho: 686.3 }, { p: 100, rs: 80.64, bo: 1.301, rho: 704.9 }, { p: 60, rs: 54.19, bo: 1.237, rho: 728.6 }, { p: 30, rs: 34.73, bo: 1.183, rho: 744.3 }, { p: 10, rs: 19.41, bo: 1.131, rho: 760.5 }],
  visc: [{ p: 300, mu: 0.806 }, { p: 224.5, mu: 0.669 }, { p: 140, mu: 1.2 }, { p: 60, mu: 2.24 }, { p: 10, mu: 2.34 }],
};
const opt = (o) => Object.entries(o).map(([value, label]) => ({ value, label }));

function engine(v, ctx = {}) { const keep = HYD_ACTIVE; try { return engineCore(v, ctx); } finally { HYD_ACTIVE = keep; } }
function engineCore(v, ctx = {}) {
  const prog = (x, m) => ctx.progress?.(x, m), warnings = [], recommendations = [], plots = [], tables = [], balances = [];
  const { spec, library } = resolveSpec(ctx.fluid, v.fluidSource), tune = tuneOf(v, spec);
  const tRes = clamp(num(v.tRes, 90), -20, 250), pRes = clamp(num(v.pRes, 300), 2, 1400), pRef = clamp(num(v.pRef, 100), 2, 1000), tSea = clamp(num(v.tSeabed, 4), -5, 60), pArr = clamp(num(v.pArr, 25), 1.05, 600), tArr = clamp(num(v.tArr, 30), -20, 150), tIn = clamp(num(spec.Tin, 70), -20, 250);
  if (!(num(v.pRes, 300) > 0)) throw new Error('The reservoir pressure must be a positive absolute pressure in bara.');
  prog(0.02, 'Characterising the fluid');
  const tx = tuneX(v), { ch, f, nKij } = tunedFluid(spec, tune, v.kijTable, tx), opts = ch.opts;
  HYD_ACTIVE = { table: hydrateTable(v.hydModel === 'vdwp' ? 'vdwp' : 'kihara', tx.hyd), sH: v.shFormer === true };
  lastModel = { spec: JSON.parse(JSON.stringify(spec)), kij: JSON.parse(JSON.stringify(Array.isArray(v.kijTable) ? v.kijTable : [])), aq: aqueous(spec) };
  const z7 = f.comps.reduce((s, c) => s + (c.pseudo ? c.z : 0), 0), isPseudo = f.comps.map((c) => !!c.pseudo);

  // -- property table for the flow solvers (kernel grid, tuned model)
  const nP = clamp(Math.round(num(v.nP, 22)), 6, 80), nT = clamp(Math.round(num(v.nT, 17)), 5, 70);
  const tPmin = clamp(num(v.tblPmin, 1), 0.5, 50), tPmax = clamp(num(v.tblPmax, 600), 100, 1500), tTmin = clamp(num(v.tblTmin, -30), -60, 20), tTmax = clamp(num(v.tblTmax, 170), 60, 300);
  if (num(v.tblPmax, 600) <= num(v.tblPmin, 1) || num(v.tblTmax, 170) <= num(v.tblTmin, -30)) throw new Error('The upper limit of the property-table range must be above its lower limit.');
  const table = buildTable(spec, { nP, nT, Pmin: tPmin, Pmax: tPmax, Tmin: tTmin, Tmax: tTmax, opts, onProgress: (x) => prog(0.03 + 0.25 * x, 'Building the pressure–temperature property table') });
  if (txFluid(tx) && !(nKij > 0)) refillTable(table, spec, f);
  if (tx.muG !== 1) for (const row of table.muG) for (let j = 0; j < row.length; j++) row[j] *= tx.muG;
  if (nKij > 0) { refillTable(table, spec, f); if (tx.muG !== 1) for (const row of table.muG) for (let j = 0; j < row.length; j++) row[j] *= tx.muG; warnings.push({ level: 'info', msg: `${nKij} binary interaction parameter(s) are overridden by the user table; the property table was recomputed with them.` }); }
  for (const c of [['reservoir', pRes, tRes], ['cold reference', pRef, tSea], ['arrival', pArr, tArr]]) if (c[1] > tPmax || c[1] < tPmin || c[2] > tTmax || c[2] < tTmin) warnings.push({ level: 'warn', msg: `The ${c[0]} condition (${c[1]} bara, ${c[2]} °C) lies outside the property-table range (${tPmin}–${tPmax} bara, ${tTmin}–${tTmax} °C): the flow solvers will clamp to the edge of the table.` });
  const std = stdFlash(f), rates = table.rates, hasLiq = std.vOil > 0 && std.beta < 1, hasOil = hasLiq && std.gor < 20000, gasSG = std.gasSG ?? f.MW / MW_AIR;

  // -- phase envelope
  prog(0.3, 'Tracing the phase envelope');
  const first = (() => { const X = f.n > 1 ? bubbleWilson(f, tRes) : null, q = X ? saturationPoint(f, tRes, { seed: X }) : null; return q && q.P !== null && q.X ? q : saturationPoint(f, tRes); })(), dewLike = first.P === null || first.type === 'dew';
  const env = traceEnvelope(f, { n: clamp(Math.round(num(v.nEnv, 40)), 10, 300), qualities: dewLike ? [0.5, 0.9, 0.99] : [0.1, 0.5, 0.9] });
  let sat = first;
  if (first.P !== null && !first.X) { const s2 = saturationPoint(f, tRes, { seed: seedFor(env, tRes) }); if (s2.P !== null && s2.X) sat = s2; }
  const psat = sat.P !== null && !sat.capped ? sat.P : null, psatType = psat !== null ? sat.type : null;
  if (sat.capped) warnings.push({ level: 'warn', msg: `The fluid is still two-phase at 1,200 bara and ${tRes} °C: no saturation pressure is reported. Check the heavy-end description (C7+ molar mass ${spec.c7MW} g/mol).` });

  // -- states at the key conditions
  prog(0.42, 'Flashing the key conditions');
  const conds = [{ name: 'Reservoir', P: pRes, T: tRes }, { name: 'Flowline inlet', P: pRef, T: tIn }, { name: 'Cold reference (seabed)', P: pRef, T: tSea }, { name: 'Arrival', P: pArr, T: tArr }, { name: 'Standard', P: P_STD, T: T_STD }];
  const muScale = (st) => { if (tx.muG !== 1 && st.gas) st.gas.mu *= tx.muG; return st; };
  for (const c of conds) { c.s = muScale(fluidState(f, c.P, c.T)); c.st = f.n > 1 ? stability(f, c.P, c.T + KEL) : { stable: true, tpd: 0 }; }
  const res = conds[0].s, resOil = res.oil || res.gas, vRes = res.phase === 'two' ? res.beta * res.gas.vm + (1 - res.beta) * res.oil.vm : (res.oil || res.gas).vm;
  const bo = hasOil ? vRes / std.vOil : null, bg = !hasOil || res.phase !== 'oil' ? ((res.gas || res.oil).Z * (tRes + KEL) * P_STD) / (pRes * (T_STD + KEL)) : null;

  // -- laboratory experiments
  prog(0.5, 'Simulating the laboratory experiments');
  let cce = [], dle = { rows: [], residual: null }, cvd = [], swell = [];
  const sep = separatorTest(f, Array.isArray(v.sepStages) ? v.sepStages : []);
  if (psat !== null) {
    const top = Math.max(pRes, 1.25 * psat), above = linspace(top, psat, 5).slice(0, -1), below = logspace(psat * 0.9, Math.max(0.12 * psat, 3), 9);
    cce = simulateCCE(f, tRes, psat, [...above, psat, ...below]);
    if (psatType === 'bubble') dle = simulateDLE(f, tRes, psat, [...linspace(psat, P_STD, 9).slice(1, -1), 5, P_STD].filter((p, i, a) => i === 0 || p < a[i - 1] - 0.5 || p === P_STD).sort((a, b) => b - a));
    cvd = simulateCVD(f, tRes, psat, linspace(psat, Math.max(0.08 * psat, 10), 9).slice(1));
    const gas = INJECTION_GASES[v.swellGas];
    if (gas) swell = swellingTest(spec, opts, gas.comp, tRes);
  }

  // -- property curves
  prog(0.58, 'Evaluating property curves');
  const nIso = clamp(Math.round(num(v.nIso, 24)), 8, 120), pTop = Math.max(pRes * 1.1, (psat || 0) * 1.15, 50), isoP = logspace(1, pTop, nIso);
  const iso = isoP.map((P) => ({ P, s: muScale(fluidState(f, P, tRes)) })), isoCold = isoP.map((P) => ({ P, s: props(f, P, tSea, { thermal: false }) }));
  const ser = (name, arr, get, extra = {}) => { const x = [], y = []; for (const r of arr) { const q = get(r); if (Number.isFinite(q)) { x.push(r.P ?? r.T); y.push(q); } } return { name, x, y, ...extra }; };
  const isoT = linspace(Math.min(tSea, 0) - 5, Math.max(tRes, tIn) + 10, Math.max(8, Math.round(nIso * 0.7))).map((T) => ({ T, s: props(f, pRef, T) }));

  // -- aqueous phase
  prog(0.64, 'Aqueous phase and hydrate equilibrium');
  const models = { act: ['NRTL', 'UNIQUAC', 'Wilson', 'ideal'].includes(v.actModel) ? v.actModel : 'NRTL', elec: ['pitzer', 'davies', 'edh', 'dh'].includes(v.elecModel) ? v.elecModel : 'pitzer' };
  const aqK = aqueous(spec), aqCase = { S: aqK.S, inhId: aqK.inhId, inhWt: aqK.inhWt };
  const custom = v.aqSource === 'custom', brine = custom ? brineFromIons(v.aqIons) : null, aqCus = custom ? { S: brine ? clamp(brine.S, 0, 260) : clamp(num(v.salinityIn, 3.5), 0, 26) * 10, inhId: INHIBITORS[v.inhIn] && v.inhIn !== 'none' && num(v.inhWtIn, 0) > 0 ? v.inhIn : 'none', inhWt: clamp(num(v.inhWtIn, 0), 0, 90) } : null;
  const awOf = (aq) => (TK) => waterActivity(aq, TK, models, tx).aw, molOf = (aq) => waterActivity(aq, 277, models).molality;
  const awCase = waterActivity(aqCase, tSea + KEL, models, tx);
  if (awCase.idealInh) warnings.push({ level: 'warn', msg: `No activity-coefficient parameters are held for ${INHIBITORS[aqCase.inhId].name}: its solution with water is treated as ideal (Raoult's law), which overstates the inhibition of alcohols and understates that of glycols by up to about 20 %.` });
  if (models.elec !== 'pitzer' && molOf(aqCase) > 0.7) warnings.push({ level: 'warn', msg: `The ${models.elec === 'davies' ? 'Davies' : 'Debye–Hückel'} model is used at ${molOf(aqCase).toFixed(2)} mol/kg, beyond its range of validity (about 0.5 mol/kg): select the Pitzer model for this brine.` });

  // -- hydrate curve
  const nHyd = clamp(Math.round(num(v.nHyd, 24)), 20, 120), hP = logspace(5, 500, nHyd), useVdw = v.hydModel !== 'corr';
  const dep = hydrateDepression(aqK), corr = hydrateScreening(table, gasSG); // what the other suites use before this suite has run
  const curveRaw = (aq) => { let g = 8; const T = [], st = []; for (const P of hP) { const r = hydrateTofP(f, P, aq ? awOf(aq) : 1, { molality: aq ? molOf(aq) : 0, guess: g }); if (!r || !Number.isFinite(r.T)) return null; T.push(r.T); st.push(r.structure); g = r.T; } return { T, st }; };
  let hyd0 = useVdw ? curveRaw(null) : null; if (hyd0) hyd0 = { T: hyd0.T.map((q, i) => publishedHydrateT(q, q, 0, hP[i])), st: hyd0.st, raw: hyd0.T };
  // inhibited curves: the model value plus the conservative allowance on the inhibitor depression (publishedHydrateT); `raw` keeps the model value
  const curve = (aq) => { const c = curveRaw(aq); if (!c || !hyd0) return c; const w = aq && aq.inhId !== 'none' ? aq.inhWt : 0; return { T: c.T.map((q, i) => publishedHydrateT(q, hyd0.raw[i], w, hP[i])), st: c.st, raw: c.T }; };
  let hydC = hyd0 ? curve(aqCase) : null, hydModelUsed = 'van der Waals–Platteeuw';
  if (!hyd0 || !hydC) {
    if (useVdw) warnings.push({ level: 'warn', msg: 'The van der Waals–Platteeuw solution did not exist over the whole 5–500 bara range (no hydrate former in the fluid or no root): the published curve falls back to the gas-gravity screening correlation with the corrected Nielsen–Bucklin depression.' });
    hyd0 = { T: hP.map(corr), st: hP.map(() => 'sII'), raw: hP.map(corr) }; hydC = { T: hP.map((P) => corr(P) - dep), st: hyd0.st, raw: hP.map((P) => corr(P) - dep) }; hydModelUsed = 'gas-gravity screening correlation with the corrected Nielsen–Bucklin depression';
  }
  const hydX = custom && hydModelUsed.startsWith('van') ? curve(aqCus) : null;
  const lnHP = hP.map(Math.log), hydAt = (arr, P) => interp1(lnHP, arr, Math.log(clamp(P, 5, 500)));
  const vdw = hydModelUsed.startsWith('van'), refFresh = vdw ? hydrateTofP(f, pRef, 1) : null, refCase = vdw ? hydrateTofP(f, pRef, awOf(aqCase), { molality: molOf(aqCase) }) : null;
  const wCase = aqCase.inhId !== 'none' ? aqCase.inhWt : 0, tHyd0Model = refFresh ? refFresh.T : corr(pRef), tHyd0 = vdw ? publishedHydrateT(tHyd0Model, tHyd0Model, 0, pRef) : tHyd0Model, tHydModel = refCase ? refCase.T : corr(pRef) - dep, tHyd = vdw ? publishedHydrateT(tHydModel, tHyd0Model, wCase, pRef) : tHydModel, allowance = tHyd - tHydModel;
  const subcool = tHyd - tSea, depression = tHyd0 - tHyd, structure = refCase?.structure || refFresh?.structure || 'sII';
  const margin = clamp(num(v.margin, 3), 0, 15), inhDesign = aqCase.inhId !== 'none' ? aqCase.inhId : INHIBITORS[v.inhDesign] && v.inhDesign !== 'none' ? v.inhDesign : 'MEG';
  const tWith = (w, id = inhDesign) => { if (!vdw) return corr(pRef) - hydrateDepression({ S: aqCase.S, inhWt: w, inh: INHIBITORS[id] }); const aq = { S: aqCase.S, inhId: w > 0 ? id : 'none', inhWt: w }, r = hydrateTofP(f, pRef, awOf(aq), { molality: molOf(aq), guess: tHyd }); return r ? publishedHydrateT(r.T, tHyd0Model, w, pRef) : -45; };
  let inhReq = 0; const target = tSea - margin;
  const tNear = (w) => { if (!vdw) return tWith(w); const aq = { S: aqCase.S, inhId: w > 0 ? inhDesign : 'none', inhWt: w }, r = hydrateTofP(f, pRef, awOf(aq), { molality: molOf(aq), guess: target, passes: 1 }); return r ? publishedHydrateT(r.T, tHyd0Model, w, pRef) : -45; }; // guest fugacities taken at the target temperature
  if (tWith(0) > target) { const g = (w) => tNear(w) - target; inhReq = g(85) > 0 ? 85 : brent(g, 0, 85, 1e-3, 40); }
  const inhNB = Math.max(0, inhibitorFor(tHyd0 - target, inhDesign, aqCase.S)), doseRows = [0, 10, 20, 30, 40, 50, 60].map((w) => ({ w, model: tHyd0 - tWith(w), hamm: hydrateDepression({ S: aqCase.S }) + hammerschmidt(w, inhDesign), nb: hydrateDepression({ S: aqCase.S, inhWt: w, inh: INHIBITORS[inhDesign] }) }));
  const pDiss = vdw ? hydratePofT(f, tSea, awOf(aqCase), { molality: molOf(aqCase) }) : null, pDissSea = pDiss ? pDiss.P : null;
  const occ = refCase?.occupancy || refFresh?.occupancy || null;
  // sensitivity of the inhibitor requirement to pressure, seabed temperature, margin and salinity
  const reqAt = (P, tS, S, mg) => {
    const tg = tS - mg, tw = (w) => { if (!vdw) return corr(P) - hydrateDepression({ S, inhWt: w, inh: INHIBITORS[inhDesign] }); const aq = { S, inhId: w > 0 ? inhDesign : 'none', inhWt: w }, r = hydrateTofP(f, P, awOf(aq), { molality: molOf(aq), guess: tg, passes: 1 }); return r ? publishedHydrateT(r.T, hydAt(hyd0.raw, P), w, P) : -45; };
    if (!(tw(0) > tg)) return 0; const g = (w) => tw(w) - tg; return g(85) > 0 ? 85 : brent(g, 0, 85, 1e-2, 30);
  };
  const sens = [['Base case', pRef, tSea, aqCase.S, margin], ['Pressure +50 %', Math.min(pRef * 1.5, 1000), tSea, aqCase.S, margin], ['Pressure −50 %', Math.max(pRef * 0.5, 2), tSea, aqCase.S, margin], ['Seabed 2 °C colder', pRef, tSea - 2, aqCase.S, margin], ['Seabed 2 °C warmer', pRef, tSea + 2, aqCase.S, margin], ['No safety margin', pRef, tSea, aqCase.S, 0], ['Fresh (condensed) water', pRef, tSea, 0, margin], ['Salinity doubled', pRef, tSea, Math.min(2 * aqCase.S, 260), margin]].map((r, i) => ({ name: r[0], P: r[1], tS: r[2], S: r[3], mg: r[4], w: i === 0 ? inhReq : reqAt(r[1], r[2], r[3], r[4]) }));

  // -- water content, gas solubility
  const wcBuk = waterContent(pArr, tArr) * awOf(aqCase)(tArr + KEL) * tx.wc, wcRaoult = waterContentRaoult(pArr, tArr, awOf(aqCase)(tArr + KEL));
  const fugRef = formerFugacities(f, pRef, tSea), sol = gasSolubility(fugRef, pRef, tSea, molOf(aqCase));
  // water distribution between the gas and the aqueous phase (kg/d) for the case rates
  const waterDist = [conds[1], conds[2], conds[3]].map((c) => { const wc = waterContent(c.P, c.T) * awOf(aqCase)(c.T + KEL) * tx.wc, q = lookup(table, c.P, c.T), sat = wc * ((q.wG * rates.mHC) / (q.mwG * 1e-3)) * VM_STD * 86400, tot = rates.mW * 86400 * (1 - aqCase.inhWt / 100) * (1 - aqCase.S / 1000), inGas = Math.min(sat, tot); return { name: c.name, P: c.P, T: c.T, wc, tot, inGas, free: tot - inGas, satCap: sat }; });

  // -- wax and asphaltenes
  prog(0.86, 'Wax and asphaltene screening');
  const dist = scnDistribution(z7, clamp(num(spec.c7MW, 210), 96, 600), clamp(num(spec.c7SG, 0.84), 0.7, 1.05), { model: v.scnModel === 'gamma' ? 'gamma' : 'exp', alpha: num(v.gammaAlpha, 1) });
  const x7Of = (x) => x.reduce((s, xi, i) => s + (isPseudo[i] ? xi : 0), 0), detect = clamp(num(v.waxDetect, 0.02), 0.001, 1), hfMult = clamp(num(v.waxHfMult, 1), 0.5, 1.5);
  const waxDead = hasLiq ? waxModel(dist, x7Of(std.x), std.mwOil, { hfMult }) : waxModel([], 0, 1), wat = waxDead.wat(detect);
  const cold = conds[2].s, waxLive = cold.oil && z7 > 0 ? waxModel(dist, x7Of(cold.x), cold.oil.MW, { hfMult }) : null, watLive = waxLive ? waxLive.wat(detect) : null;
  const waxT = linspace(-10, Math.max(wat ?? 20, 20) + 6, 30), waxAtSea = waxDead.solidWt(tSea);
  const asph = hasOil ? asphalteneScreen({ rhoRes: resOil.rho, pRes, pBub: psatType === 'bubble' ? psat : pRes, sara: { sat: num(v.saraSat, 0), aro: num(v.saraAro, 0), res: num(v.saraRes, 0), asp: num(v.saraAsp, 0) } }) : null;

  // -- isenthalpic / isentropic expansion from the reservoir to arrival pressure
  const s0 = stateHSV(f, pRes, tRes), jtOut = flashPH(f, pArr, s0.H), isOut = flashPS(f, pArr, s0.S);

  // ---- cross-checks -----------------------------------------------------------------------------------------------------
  const sepGor = sep.gor, bo0 = hasOil && sep.api !== null && sep.gasSG ? blackOil({ rs: sepGor ?? 0, gasSG: sep.gasSG, api: sep.api, T: tRes, P: pRes }) : null;
  const satState = psat !== null ? fluidState(f, psat * 1.0005, tRes) : null, boSat = psat !== null && sep.vSto > 0 && psatType === 'bubble' ? satState.oil.vm / sep.vSto : null;
  const refGas = conds.find((c) => c.s.gas && c.name !== 'Standard') || conds[4], gq = refGas.s.gas || refGas.s.oil, kg = kay(f, refGas.s.gas ? refGas.s.y : refGas.s.x);
  const zLK = leeKeslerZ((refGas.T + KEL) / kg.Tc, refGas.P / kg.Pc, kg.w), muLGE = viscosityLGE(gq.MW, gq.rho, refGas.T + KEL), k0 = phaseProps(f, refGas.s.gas ? refGas.s.y : refGas.s.x, 1, refGas.T, 'vapour').k, kST = conductivityStielThodos(k0, (kg.Vc * 1e-6) / gq.vm, kg);

  // ---- laboratory comparison ------------------------------------------------------------------------------------------------
  const lab = { rows: [], aad: {} };
  if (v.useLab !== false) {
    const tb = (a) => (Array.isArray(a) ? a : []);
    const ps = tb(v.labPsat).filter((r) => isNum(+r.t) && isNum(+r.p) && +r.p > 0).map((r) => ({ t: +r.t, m: +r.p, c: saturationPoint(f, +r.t, { seed: seedFor(env, +r.t) }).P }));
    lab.psat = ps; if (ps.length) lab.aad.psat = aad(ps.map((r) => r.m), ps.map((r) => r.c));
    if (psat !== null) {
      const cr = tb(v.labCCE).filter((r) => isNum(+r.p) && isNum(+r.v) && +r.p > 1), cm = cr.length ? simulateCCE(f, tRes, psat, cr.map((r) => +r.p)) : [];
      lab.cce = cr.map((r, i) => ({ p: +r.p, m: +r.v, c: cm[i].vRel })); if (cr.length) lab.aad.cce = aad(lab.cce.map((r) => r.m), lab.cce.map((r) => r.c));
      const vr = tb(v.labVisc).filter((r) => isNum(+r.p) && isNum(+r.mu) && +r.p > 0);
      lab.visc = vr.map((r) => { const s = props(f, +r.p, tRes, { thermal: false }); return { p: +r.p, m: +r.mu, c: s.oil.mu * 1e3 }; }); if (vr.length) lab.aad.visc = aad(lab.visc.map((r) => r.m), lab.visc.map((r) => r.c));
      if (psatType === 'bubble') {
        const dr = tb(v.labDLE).filter((r) => isNum(+r.p) && +r.p > 0 && +r.p < psat * 0.995).sort((a, b) => b.p - a.p), dm = dr.length ? simulateDLE(f, tRes, psat, [...dr.map((r) => +r.p), P_STD]).rows : [];
        lab.dle = dr.map((r, i) => ({ p: +r.p, rs: +r.rs, bo: +r.bo, rho: +r.rho, c: dm[i + 1] })).filter((r) => r.c);
        if (lab.dle.length) { lab.aad.rs = aad(lab.dle.map((r) => r.rs), lab.dle.map((r) => r.c.rs)); lab.aad.bo = aad(lab.dle.map((r) => r.bo), lab.dle.map((r) => r.c.bo)); lab.aad.rho = aad(lab.dle.map((r) => r.rho), lab.dle.map((r) => r.c.rhoO)); }
      }
    }
    const hr = tb(v.labHyd).filter((r) => isNum(+r.p) && isNum(+r.t) && +r.p >= 5 && +r.p <= 500);
    lab.hyd = hr.map((r) => ({ p: +r.p, m: +r.t, c: hydAt(hyd0.T, +r.p) })); if (hr.length) lab.hydDev = lab.hyd.reduce((a, r) => a + Math.abs(r.c - r.m), 0) / hr.length;
    const lw = num(v.labWat, 0); if (lw !== 0 && wat !== null) lab.wat = { m: lw, c: wat };
  }

  // ---- balances ---------------------------------------------------------------------------------------------------------
  { const c = conds[2].s; let out = 0; f.z.forEach((_, i) => (out += c.phase === 'two' ? c.beta * c.y[i] + (1 - c.beta) * c.x[i] : c.x[i])); balances.push({ name: 'Moles over the flash at the cold reference point (per mole of feed)', in: 1, out }); }
  { const i = f.comps.findIndex((c) => c.id === 'C1'); if (i >= 0) { const c = conds[3].s; balances.push({ name: 'Methane over the arrival flash (mole fraction of feed)', in: f.z[i], out: c.phase === 'two' ? c.beta * c.y[i] + (1 - c.beta) * c.x[i] : c.x[i] }); } }
  balances.push({ name: 'Hydrocarbon mass at standard conditions (kg/s)', in: rates.mHC, out: (rates.qOilStd * (rates.rhoOilStd || 0) + (rates.qGasStd * (rates.mwGas || 0) * 1e-3) / VM_STD) / 86400 });
  { const last = sep.stages[sep.stages.length - 1], gasMol = sep.vSto > 0 ? sep.stages.reduce((q, r) => q + (r.gor * sep.vSto) / VM_STD, 0) : 1; balances.push({ name: 'Moles over the separator train (per mole of feed)', in: 1, out: gasMol + (sep.vSto > 0 ? last.liqMol : 0) }); }
  if (jtOut) balances.push({ name: 'Enthalpy across the isenthalpic expansion (J/mol)', in: s0.H, out: jtOut.H });

  // ---- warnings and recommendations ----------------------------------------------------------------------------------------
  const F = (x, d = 1) => (Number.isFinite(x) ? (+x).toFixed(d) : '—'), inhName = INHIBITORS[inhDesign].name;
  if (library) warnings.push({ level: 'info', msg: `The built-in fluid “${library.name}” is analysed instead of the case fluid. The published property table belongs to this library fluid, so the other suites keep using their own table of the case fluid.` });
  if (tune.method !== 'KL') warnings.push({ level: 'info', msg: `${tune.method === 'Twu' ? 'Twu' : 'Riazi–Daubert'} critical properties are applied as plus-fraction-average ratios to Kesler–Lee (Tc × ${F(ch.ratios.Tc, 4)}, Pc × ${F(ch.ratios.Pc, 4)}, ω × ${F(ch.ratios.w, 4)}).` });
  if (psat !== null && psatType === 'bubble' && pRes < psat) warnings.push({ level: 'warn', msg: `The reservoir pressure (${F(pRes)} bara) is below the bubble point (${F(psat)} bara): the fluid is already two-phase in the reservoir (${F(100 * res.beta)} mol % vapour).` });
  if (psat !== null && psatType === 'dew' && pRes < psat) warnings.push({ level: 'warn', msg: `The reservoir pressure (${F(pRes)} bara) is below the dew point (${F(psat)} bara): retrograde liquid drops out in the reservoir.` });
  if (subcool > 0) warnings.push({ level: subcool > 3 ? 'bad' : 'warn', msg: `At ${F(pRef)} bara and the seabed temperature of ${F(tSea)} °C the fluid is ${F(subcool)} °C inside the hydrate region (equilibrium ${F(tHyd)} °C, ${structure}).` });
  if (wat !== null && wat > tSea) warnings.push({ level: 'warn', msg: `The wax appearance temperature of the stock-tank oil (${F(wat)} °C) is above the seabed temperature (${F(tSea)} °C): ${F(waxAtSea, 2)} wt % of the oil can precipitate at the wall.` });
  if (asph && asph.risk !== 'low') warnings.push({ level: asph.risk === 'high' ? 'bad' : 'warn', msg: `Asphaltene screening: de Boer class “${asph.deBoer}” (undersaturation ${F(asph.dP)} bar at ${F(resOil.rho, 0)} kg/m³)${asph.cii !== null ? `, colloidal instability index ${F(asph.cii, 2)} (${asph.ciiClass})` : ''}.` });
  if (env.T.length && !env.critical) warnings.push({ level: 'info', msg: 'No critical point lies on the traced part of the saturation line (above −40 °C and between 1 and 1,200 bara).' });
  if (f.eosId === 'RK' || f.eosId === 'vdW') warnings.push({ level: 'warn', msg: `The ${EOS[f.eosId].name} equation is kept for teaching and benchmarking: its liquid densities and saturation pressures are not of engineering quality. Use Peng–Robinson or Soave–Redlich–Kwong for design.` });
  for (const [k, lbl] of [['psat', 'saturation pressure'], ['cce', 'CCE relative volume'], ['rs', 'solution gas–oil ratio'], ['bo', 'oil formation-volume factor'], ['rho', 'oil density'], ['visc', 'oil viscosity']]) if (lab.aad[k] > (k === 'visc' ? 15 : 4)) warnings.push({ level: 'warn', msg: `Laboratory ${lbl}: average absolute deviation ${F(lab.aad[k])} % — the model is not yet tuned to this measurement.` });

  if (subcool > -margin) {
    const wNeed = inhReq, vol = rates.mW > 0 && wNeed < 85 ? ((rates.mW * wNeed) / (100 - wNeed) / INHIBITORS[inhDesign].rho) * 86400 : null;
    recommendations.push(`Hydrate control: the dissociation temperature at ${F(pRef)} bara is ${F(tHyd)} °C against ${F(tSea)} °C at the seabed. To hold a ${F(margin)} °C margin, dose ${F(wNeed)} wt % ${inhName} in the aqueous phase${vol !== null ? ` (about ${F(vol, 2)} m³/d for the case water rate of ${F(rates.qWaterStd)} Sm³/d)` : ''}; the Nielsen–Bucklin estimate is ${F(inhNB)} wt %. Alternatively keep the fluid above ${F(tHyd + margin)} °C with insulation or heating — pass this curve to the flow and operations suites to size the cooldown time.`);
  } else recommendations.push(`Hydrates: the cold reference point (${F(pRef)} bara, ${F(tSea)} °C) is ${F(-subcool)} °C outside the hydrate region with the present aqueous phase (${aqCase.inhId !== 'none' ? `${F(aqCase.inhWt)} wt % ${INHIBITORS[aqCase.inhId].name}, ` : ''}${F(aqCase.S / 10)} wt % salt): no additional inhibitor is needed for a ${F(margin)} °C margin. Re-check at shut-in pressure.`);
  if (inhReq > 0) { const hi = sens.slice(1).reduce((a, r) => (r.w > a.w ? r : a), sens[1]), lo = sens.slice(1).reduce((a, r) => (r.w < a.w ? r : a), sens[1]); recommendations.push(`Inhibitor sensitivity: the ${inhName} requirement ranges from ${F(lo.w)} wt % (${lo.name.toLowerCase()}) to ${F(hi.w)} wt % (${hi.name.toLowerCase()}) around the base value of ${F(inhReq)} wt %; size the injection system for the upper value and for the start-up water cut.`); }
  if (pDissSea !== null && subcool > 0) { const gasLine = conds[2].s.wG > 0.5, head = (num(v.waterDepth, 0) * G0 * (gasLine ? (conds[3].s.gas || conds[2].s.gas).rho : conds[2].s.oil.rho)) / 1e5; recommendations.push(`Depressurisation: hydrate dissociates below ${F(pDissSea)} bara at ${F(tSea)} °C. A ${gasLine ? 'gas' : 'liquid'}-filled riser of ${F(num(v.waterDepth, 0), 0)} m exerts about ${F(head)} bar at its base, so blowdown from the host ${head + 2 > pDissSea ? 'cannot bring the flowline out of the hydrate region on its own — plan displacement, dead-oil circulation or subsea depressurisation' : 'can bring the flowline out of the hydrate region'}.`); }
  if (wat !== null) recommendations.push(wat > tSea ? `Wax: keep the wall above ${F(wat)} °C (stock-tank oil; ${watLive !== null ? `${F(watLive)} °C for the live oil at ${F(pRef)} bara` : 'no live-oil value'}) or plan pigging — about ${F(waxAtSea, 2)} wt % precipitates at ${F(tSea)} °C and the wax-forming C18+ content is ${F(waxDead.waxContent)} wt %. Confirm the appearance temperature by cross-polar microscopy or DSC and enter it as laboratory data.` : `Wax: the predicted appearance temperature (${F(wat)} °C) is below the seabed temperature (${F(tSea)} °C); wax deposition is not expected in steady operation.`);
  if (psat !== null) recommendations.push(psatType === 'bubble' ? `Phase behaviour: the bubble point at ${F(tRes)} °C is ${F(psat)} bara, ${pRes > psat ? `${F(pRes - psat)} bar below reservoir pressure — gas first breaks out where the flowing pressure falls to that value, so expect two-phase flow from ${psat > pRef ? 'the well tubing onwards' : 'the flowline onwards'}` : 'above the reservoir pressure'}. Design the separator train on the multi-stage result (${sepGor !== null ? F(sepGor) : '—'} Sm³/Sm³, ${sep.api !== null ? F(sep.api) : '—'} °API) rather than on the single-stage flash (${rates.gor !== null ? F(rates.gor) : '—'} Sm³/Sm³).` : `Phase behaviour: the dew point at ${F(tRes)} °C is ${F(psat)} bara and the cricondentherm is ${F(env.cricondentherm?.T)} °C; the maximum retrograde liquid drop-out in constant-volume depletion is ${F(Math.max(0, ...cvd.map((r) => r.liqPct)), 2)} % of the pore volume. Keep the flowing bottom-hole pressure above the dew point as long as practical.`);
  const worst = Object.entries(lab.aad).filter(([, x]) => x !== null).sort((a, b) => b[1] - a[1])[0];
  if (worst && worst[1] > 2) recommendations.push(`Tuning: the largest deviation from the laboratory report is ${F(worst[1])} % (${worst[0]}). Open the calibration tab and regress the C7+ critical temperature and pressure multipliers and the methane–C7+ interaction scale to the saturation pressure and density data, then the critical-volume multiplier to viscosity; validate on the points held back.`);
  else if (v.useLab === false || !worst) recommendations.push('Tuning: no laboratory data are compared in this run. Enter the saturation pressure, CCE, differential-liberation and viscosity points of the PVT report on the inputs tab — an untuned cubic EOS is typically within 5–10 % on saturation pressure and 20–50 % on oil viscosity.');
  if (jtOut) recommendations.push(`Expansion: throttling the reservoir fluid from ${F(pRes)} bara, ${F(tRes)} °C to ${F(pArr)} bara at constant enthalpy gives ${F(jtOut.T)} °C (isentropic limit ${isOut ? F(isOut.T) : '—'} °C)${jtOut.T < tHyd0 ? ', which is inside the hydrate region — check choke and restart scenarios' : ''}.`);
  if (lab.hydDev > 1.5) warnings.push({ level: 'warn', msg: `The laboratory hydrate points deviate by ${F(lab.hydDev, 2)} °C on average from the model: check the gas analysis used in the test, or switch the hydrate model.` });
  if (brine) warnings.push({ level: Math.abs(brine.balance) > 0.1 ? 'warn' : 'info', msg: `Custom brine analysis: total dissolved solids ${F(brine.tds)} g/L, ionic strength ${F(brine.I, 3)} mol/kg, treated as ${F(brine.S / 10, 2)} wt % NaCl of equal ionic strength${Math.abs(brine.balance) > 0.1 ? `; the charge balance is off by ${F(100 * brine.balance, 0)} % — check the analysis` : ''}.` });
  if (asph && asph.risk !== 'low') recommendations.push(`Asphaltenes: screening risk is ${asph.risk}. Commission a depressurisation test (onset pressure) on a live sample before finalising the completion and chemical-injection design.`);

  // ---- plots -----------------------------------------------------------------------------------------------------------------
  prog(0.93, 'Assembling the results');
  const ops = conds.slice(0, 4), yTop = Math.max(pRes, env.cricondenbar?.P || 0, 60) * 1.12, watLine = wat !== null ? [{ name: `Wax appearance (${F(wat)} °C)`, x: [wat, wat], y: [1, yTop], dash: true, color: '#b45309' }] : [];
  const envSeries = [
    ...(env.bubble.T.length ? [{ name: 'Bubble-point line', x: env.bubble.T, y: env.bubble.P, color: '#1d4ed8' }] : []),
    ...(env.dew.T.length ? [{ name: 'Dew-point line', x: env.dew.T, y: env.dew.P, color: '#dc2626' }] : []),
    ...env.quality.map((q) => ({ name: `${Math.round(q.q * 100)} mol % vapour`, x: q.T, y: q.P, dash: true })),
    ...(env.critical ? [{ name: 'Critical point', x: [env.critical.T], y: [env.critical.P], mode: 'points', color: '#111827' }] : []),
    { name: 'Hydrate curve, fresh water', x: hyd0.T, y: hP, color: '#0891b2', dash: true }, { name: 'Hydrate curve, case aqueous phase', x: hydC.T, y: hP, color: '#0891b2' },
    ...(hydX ? [{ name: 'Hydrate curve, custom aqueous phase', x: hydX.T, y: hP, color: '#7c3aed' }] : []),
    ...watLine, { name: 'Operating points', x: ops.map((c) => c.T), y: ops.map((c) => c.P), mode: 'points', color: '#16a34a' },
  ];
  plots.push({ type: 'line', title: 'Phase envelope with hydrate curve and wax appearance temperature', xlabel: 'Temperature (°C)', ylabel: 'Pressure (bara)', ymin: 0, ymax: yTop, xmin: -45, series: envSeries, note: `Operating points: ${ops.map((c) => `${c.name} ${F(c.P)} bara / ${F(c.T)} °C`).join('; ')}. ${f.eosId} equation of state, ${f.n} components.` });
  plots.push({ type: 'line', title: 'Flow-assurance window (zoom)', xlabel: 'Temperature (°C)', ylabel: 'Pressure (bara)', ymin: 0, ymax: Math.max(pRes, pRef) * 1.15, xmin: Math.min(-10, tSea - 10), xmax: Math.max(tRes, tIn) + 15, series: envSeries.map((s) => ({ ...s })) });
  plots.push({ type: 'field', title: 'Gas mass fraction over the property table', xlabel: 'Temperature (°C)', ylabel: 'Pressure (bara)', zlabel: 'Gas mass fraction', zunit: '–', x: table.T, y: table.P, z: table.wG, zmin: 0, zmax: 1, cmap: 'viridis', contours: 8,
    shapes: [...(env.T.length ? [{ x: env.T, y: env.P, color: '#ffffff' }] : []), { x: hydC.T, y: hP, color: '#38bdf8', dash: true }], markers: ops.map((c) => ({ x: c.T, y: c.P, label: c.name })) });
  plots.push({ type: 'field', title: 'Mixture density over the property table (no-slip, hydrocarbon phases)', xlabel: 'Temperature (°C)', ylabel: 'Pressure (bara)', zlabel: 'Density', zunit: 'kg/m³', x: table.T, y: table.P, z: table.P.map((_, i) => table.T.map((__, j) => { const w = table.wG[i][j]; return 1 / (w / table.rhoG[i][j] + (1 - w) / table.rhoO[i][j]); })), cmap: 'turbo', contours: 8 });
  const tl = `at ${F(tRes)} °C`;
  plots.push({ type: 'line', title: `Phase densities ${tl}`, xlabel: 'Pressure (bara)', ylabel: 'Density (kg/m³)', series: [ser('Oil', iso, (r) => r.s.oil?.rho), ser('Gas', iso, (r) => r.s.gas?.rho), ser(`Oil at ${F(tSea)} °C`, isoCold, (r) => (r.s.phase !== 'gas' ? r.s.oil.rho : NaN), { dash: true })], vlines: psat !== null ? [{ x: psat, label: 'Psat' }] : [] });
  plots.push({ type: 'line', title: `Phase viscosities ${tl} (Lohrenz–Bray–Clark)`, xlabel: 'Pressure (bara)', ylabel: 'Viscosity (mPa·s)', logy: true, series: [ser('Oil', iso, (r) => (r.s.oil ? r.s.oil.mu * 1e3 : NaN)), ser('Gas', iso, (r) => (r.s.gas ? r.s.gas.mu * 1e3 : NaN)), ser('Gas, Lee–Gonzalez–Eakin', iso, (r) => (r.s.gas ? viscosityLGE(r.s.gas.MW, r.s.gas.rho, tRes + KEL) * 1e3 : NaN), { dash: true }), ser(`Oil at ${F(tSea)} °C`, isoCold, (r) => (r.s.phase !== 'gas' ? r.s.oil.mu * 1e3 : NaN), { dash: true })] });
  plots.push({ type: 'line', title: `Compressibility factor and vapour fraction ${tl}`, xlabel: 'Pressure (bara)', ylabel: 'Z, vapour mole fraction (–)', series: [ser('Z gas', iso, (r) => r.s.gas?.Z), ser('Z oil', iso, (r) => r.s.oil?.Z), ser('Vapour mole fraction', iso, (r) => r.s.beta, { dash: true })] });
  plots.push({ type: 'line', title: `Heat capacity and Joule–Thomson coefficient ${tl}`, xlabel: 'Pressure (bara)', ylabel: 'cp (kJ/kg/K), JT (K/bar)', series: [ser('cp oil', iso, (r) => (r.s.oil ? r.s.oil.cp / 1e3 : NaN)), ser('cp gas', iso, (r) => (r.s.gas ? r.s.gas.cp / 1e3 : NaN)), ser('JT gas', iso, (r) => r.s.gas?.jtBar, { dash: true }), ser('JT oil', iso, (r) => r.s.oil?.jtBar, { dash: true })] });
  plots.push({ type: 'line', title: `Speed of sound and isothermal compressibility ${tl}`, xlabel: 'Pressure (bara)', ylabel: 'c (m/s), 10⁵·κT (1/bar)', logy: true, series: [ser('Sound speed, oil', iso, (r) => r.s.oil?.sound), ser('Sound speed, gas', iso, (r) => r.s.gas?.sound), ser('κT oil × 10⁵', iso, (r) => (r.s.oil ? r.s.oil.kT * 1e5 : NaN), { dash: true }), ser('κT gas × 10⁵', iso, (r) => (r.s.gas ? r.s.gas.kT * 1e5 : NaN), { dash: true })] });
  plots.push({ type: 'line', title: `Thermal conductivity and interfacial tension ${tl}`, xlabel: 'Pressure (bara)', ylabel: 'k (W/m/K), σ (mN/m ÷ 100)', series: [ser('k oil', iso, (r) => r.s.oil?.k), ser('k gas', iso, (r) => r.s.gas?.k), ser('Gas–oil tension ÷ 100 (parachor)', iso, (r) => (r.s.sigma !== null ? r.s.sigma * 10 : NaN), { dash: true })] });
  plots.push({ type: 'line', title: `Properties against temperature at ${F(pRef)} bara`, xlabel: 'Temperature (°C)', ylabel: 'Enthalpy (kJ/kg ÷ 100), gas mass fraction, oil viscosity (mPa·s)', series: [ser('Mixture enthalpy ÷ 100', isoT, (r) => (r.s.wG * r.s.gas.h + (1 - r.s.wG) * r.s.oil.h) / 1e5), ser('Gas mass fraction', isoT, (r) => r.s.wG), ser('Oil viscosity', isoT, (r) => (r.s.phase !== 'gas' ? r.s.oil.mu * 1e3 : NaN), { dash: true })], vlines: [{ x: tHyd, label: 'Hydrate' }, ...(wat !== null ? [{ x: wat, label: 'WAT' }] : [])] });
  if (cce.length) {
    plots.push({ type: 'line', title: `Constant composition expansion ${tl}`, xlabel: 'Pressure (bara)', ylabel: 'Relative volume V/Vsat, liquid (fraction of Vsat)', series: [ser('Relative volume', cce, (r) => r.vRel, { mode: 'both' }), ser('Liquid volume / Vsat', cce, (r) => r.liqPct / 100, { mode: 'both', dash: true }), ...(lab.cce?.length ? [{ name: 'Laboratory relative volume', x: lab.cce.map((r) => r.p), y: lab.cce.map((r) => r.m), mode: 'points' }] : [])], vlines: [{ x: psat, label: psatType === 'bubble' ? 'Bubble point' : 'Dew point' }] });
    if (cce.some((r) => r.Y !== null)) plots.push({ type: 'line', title: 'CCE Y-function below the saturation pressure', xlabel: 'Pressure (bara)', ylabel: 'Y = (Psat − P) / (P (Vrel − 1))', series: [ser('Y-function', cce, (r) => r.Y ?? NaN, { mode: 'both' })], note: 'A straight Y-function is the usual consistency check of the measured two-phase volumes.' });
  }
  if (dle.rows.length) {
    plots.push({ type: 'line', title: `Differential liberation ${tl}: solution gas and formation-volume factor`, xlabel: 'Pressure (bara)', ylabel: 'Rs (Sm³/Sm³ ÷ 100), Bo (m³/Sm³)', series: [ser('Rs ÷ 100', dle.rows, (r) => r.rs / 100, { mode: 'both' }), ser('Bo', dle.rows, (r) => r.bo, { mode: 'both' }), ser('Bt', dle.rows.filter((r) => r.P > 0.15 * psat), (r) => r.bt, { dash: true }), ...(lab.dle?.length ? [{ name: 'Laboratory Rs ÷ 100', x: lab.dle.map((r) => r.p), y: lab.dle.map((r) => r.rs / 100), mode: 'points' }, { name: 'Laboratory Bo', x: lab.dle.map((r) => r.p), y: lab.dle.map((r) => r.bo), mode: 'points' }] : [])] });
    plots.push({ type: 'line', title: `Differential liberation ${tl}: oil density and gas gravity`, xlabel: 'Pressure (bara)', ylabel: 'Oil density (kg/m³ ÷ 1000), gas gravity (air = 1)', series: [ser('Oil density ÷ 1000', dle.rows, (r) => (r.rhoO ?? NaN) / 1000, { mode: 'both' }), ser('Liberated-gas gravity', dle.rows, (r) => r.gasSG ?? NaN, { mode: 'both' }), ser('Gas Z', dle.rows, (r) => r.Z ?? NaN, { dash: true })] });
  }
  if (cvd.length > 1) plots.push({ type: 'line', title: `Constant volume depletion ${tl}`, xlabel: 'Pressure (bara)', ylabel: 'Liquid (% of cell ÷ 100), produced (mol fraction), Z', series: [ser('Liquid saturation', cvd, (r) => r.liqPct / 100, { mode: 'both' }), ser('Cumulative produced', cvd, (r) => r.produced / 100, { mode: 'both' }), ser('Two-phase Z', cvd, (r) => r.Z2, { dash: true })], note: psatType === 'bubble' ? 'For an oil the constant-volume depletion describes solution-gas drive below the bubble point; for a gas condensate it gives the retrograde liquid drop-out.' : '' });
  plots.push({ type: 'bar', title: 'Separator test: gas released per stage', ylabel: 'GOR (Sm³ gas per Sm³ stock-tank oil)', categories: sep.stages.map((s) => `${F(s.p, s.p < 10 ? 2 : 0)} bara / ${F(s.t, 0)} °C`), series: [{ name: 'Stage GOR', values: sep.stages.map((s) => s.gor ?? 0) }] });
  if (swell.length) plots.push({ type: 'line', title: `Swelling test with ${INJECTION_GASES[v.swellGas].name.toLowerCase()} ${tl}`, xlabel: 'Injected gas (mol per mol of reservoir fluid)', ylabel: 'Psat (bara ÷ 100), swelling factor', series: [{ name: 'Saturation pressure ÷ 100', x: swell.filter((r) => r.psat !== null).map((r) => r.inj), y: swell.filter((r) => r.psat !== null).map((r) => r.psat / 100), mode: 'both' }, { name: 'Swelling factor', x: swell.filter((r) => r.swell !== null).map((r) => r.inj), y: swell.filter((r) => r.swell !== null).map((r) => r.swell), mode: 'both' }] });
  plots.push({ type: 'line', title: `Hydrate dissociation curve (${hydModelUsed})`, xlabel: 'Temperature (°C)', ylabel: 'Pressure (bara)', logy: true, series: [{ name: 'Fresh water', x: hyd0.T, y: hP }, { name: `Case aqueous phase (${F(aqCase.S / 10)} wt % salt${aqCase.inhId !== 'none' ? `, ${F(aqCase.inhWt)} wt % ${aqCase.inhId}` : ''})`, x: hydC.T, y: hP }, ...(hydX ? [{ name: `Custom (${F(aqCus.S / 10)} wt % salt${aqCus.inhId !== 'none' ? `, ${F(aqCus.inhWt)} wt % ${aqCus.inhId}` : ''})`, x: hydX.T, y: hP }] : []), { name: 'Gas-gravity correlation (Motiee)', x: hP.map(corr), y: hP, dash: true }, { name: 'Cold reference point', x: [tSea], y: [pRef], mode: 'points' }], note: `Structure ${structure} at ${F(pRef)} bara. Stable hydrate lies to the left of the curve.` });
  plots.push({ type: 'line', title: `Hydrate depression by ${inhName} at ${F(pRef)} bara`, xlabel: 'Inhibitor in the aqueous phase (wt %)', ylabel: 'Depression of the hydrate temperature (°C)', series: [{ name: `${vdw ? 'van der Waals–Platteeuw + ' + models.act : 'Nielsen–Bucklin'}`, x: doseRows.map((r) => r.w), y: doseRows.map((r) => r.model), mode: 'both' }, { name: 'Hammerschmidt (+ salt)', x: doseRows.map((r) => r.w), y: doseRows.map((r) => r.hamm), dash: true }, { name: 'Nielsen–Bucklin', x: doseRows.map((r) => r.w), y: doseRows.map((r) => r.nb), dash: true }], hlines: [{ y: tHyd0 - target, label: `Required (${F(margin)} °C margin)` }] });
  { const gT = linspace(Math.min(-5, tSea - 5), 35, 41); plots.push({ type: 'field', title: 'Hydrate stability region: subcooling with the case aqueous phase', xlabel: 'Temperature (°C)', ylabel: 'Pressure (bara)', zlabel: 'Subcooling (positive = hydrate stable)', zunit: '°C', x: gT, y: hP, z: hP.map((_, i) => gT.map((t) => hydC.T[i] - t)), cmap: 'coolwarm', contours: 8, shapes: [{ x: hydC.T, y: hP, color: '#111827' }], markers: [{ x: tSea, y: clamp(pRef, 5, 500), label: 'Cold reference' }] }); }
  plots.push({ type: 'bar', title: `Sensitivity of the ${inhName} requirement`, ylabel: 'Required inhibitor (wt % of the aqueous phase)', categories: sens.map((r) => r.name), series: [{ name: 'Required', values: sens.map((r) => r.w) }] });
  if (occ) { const ids = [...new Set([...Object.keys(occ.byGuest[0]), ...Object.keys(occ.byGuest[1])])]; plots.push({ type: 'bar', title: `Cage occupancy of the ${structure} hydrate at ${F(pRef)} bara`, ylabel: 'Fractional occupancy', categories: ['Small cavity', 'Large cavity'], stacked: true, series: ids.map((id) => ({ name: COMP_LABELS[id] || id, values: [occ.byGuest[0][id] || 0, occ.byGuest[1][id] || 0] })) }); }
  if (hasLiq && z7 > 0 && wat !== null) plots.push({ type: 'line', title: 'Wax precipitation curve (ideal-solution solid–liquid equilibrium)', xlabel: 'Temperature (°C)', ylabel: 'Solid wax (wt % of the oil)', zeroY: true, series: [{ name: 'Stock-tank oil', x: waxT, y: waxT.map((t) => waxDead.solidWt(t)) }, ...(waxLive ? [{ name: `Live oil at ${F(pRef)} bara`, x: waxT, y: waxT.map((t) => waxLive.solidWt(t)), dash: true }] : [])], vlines: [{ x: tSea, label: 'Seabed' }, ...(wat !== null ? [{ x: wat, label: 'WAT' }] : [])] });
  if (dist.length) plots.push({ type: 'line', title: 'Carbon-number distribution of the plus fraction', xlabel: 'Carbon number', ylabel: 'Mole fraction of the fluid', logy: true, series: [{ name: v.scnModel === 'gamma' ? 'Gamma (Whitson)' : 'Exponential (Pedersen)', x: dist.filter((d) => d.z > 1e-9).map((d) => d.n), y: dist.filter((d) => d.z > 1e-9).map((d) => d.z) }, { name: 'Wax-forming part', x: waxDead.comps.filter((c) => c.z > 1e-9).map((c) => c.n), y: waxDead.comps.filter((c) => c.z > 1e-9).map((c) => (c.z * z7) / Math.max(x7Of(std.x), 1e-12)), dash: true }] });
  { const ms = linspace(0.05, 6, 40), ref = [[0.1, 0.9324], [0.5, 0.9209], [1, 0.9355], [2, 0.9833], [3, 1.0453], [4, 1.1158], [5, 1.1916], [6, 1.2706]];
    plots.push({ type: 'line', title: 'Osmotic coefficient of NaCl brine at 25 °C: electrolyte models against tabulated data', xlabel: 'Molality (mol/kg water)', ylabel: 'Osmotic coefficient φ', ymin: 0.6, ymax: 1.5, series: [{ name: 'Pitzer', x: ms, y: ms.map((m) => osmoticNaCl('pitzer', m)) }, { name: 'Davies', x: ms, y: ms.map((m) => osmoticNaCl('davies', m)), dash: true }, { name: 'Extended Debye–Hückel', x: ms, y: ms.map((m) => osmoticNaCl('edh', m)), dash: true }, { name: 'Debye–Hückel limiting law', x: ms, y: ms.map((m) => osmoticNaCl('dh', m)), dash: true }, { name: 'Robinson & Stokes', x: ref.map((r) => r[0]), y: ref.map((r) => r[1]), mode: 'points' }], vlines: molOf(aqCase) > 0 ? [{ x: molOf(aqCase), label: 'Case brine' }] : [] }); }
  if (lab.psat?.length) plots.push({ type: 'line', title: 'Saturation pressure: model against laboratory', xlabel: 'Temperature (°C)', ylabel: 'Saturation pressure (bara)', series: [{ name: 'Model', x: env.T.filter((t) => t > -20 && t < 200), y: env.P.filter((_, i) => env.T[i] > -20 && env.T[i] < 200) }, { name: 'Laboratory', x: lab.psat.map((r) => r.t), y: lab.psat.map((r) => r.m), mode: 'points' }] });

  // ---- tables ----------------------------------------------------------------------------------------------------------------
  const cs = conds[2].s, ca = conds[3].s;
  tables.push({ title: 'Fluid composition and characterisation', columns: ['Component', 'Feed (mol %)', 'MW (g/mol)', 'Tc (K)', 'Pc (bar)', 'ω', 'Vc (cm³/mol)', 'Volume shift c/b', `Liquid at arrival (mol %)`, 'Vapour at arrival (mol %)', 'K at arrival'], rows: f.comps.map((c, i) => [c.name, cell(100 * c.z, 5), cell(c.MW, 5), cell(c.Tc, 5), cell(c.Pc, 4), cell(c.w, 4), cell(c.Vc, 4), cell(c.b > 0 ? c.c / c.b : 0, 3), cell(100 * ca.x[i], 4), cell(100 * ca.y[i], 4), ca.phase === 'two' && ca.x[i] > 1e-300 ? cell(ca.y[i] / ca.x[i], 4) : '—']), note: `${spec.name || 'Case fluid'}: ${f.eosId} equation of state${opts.shift ? ' with volume translation' : ''}, ${opts.nPseudo} pseudo-component(s) for C7+ (M = ${spec.c7MW} g/mol, SG = ${spec.c7SG}), mixture molar mass ${F(f.MW, 2)} g/mol. Arrival: ${F(pArr)} bara, ${F(tArr)} °C.` });
  if (ch.pseudo.length) tables.push({ title: 'C7+ pseudo-components: critical-property correlations', columns: ['Pseudo', 'mol %', 'MW', 'SG', 'Tb (K)', 'Tc Kesler–Lee (K)', 'Tc Riazi–Daubert', 'Tc Twu', 'Pc Kesler–Lee (bar)', 'Pc Riazi–Daubert', 'Pc Twu', 'ω Kesler–Lee', 'ω Riazi–Daubert', 'ω Twu'], rows: ch.pseudo.map((p) => [p.id, cell(100 * p.z), cell(p.MW), cell(p.SG), cell(p.Tb), cell(p.KL.Tc), cell(p.RD.Tc), cell(p.Twu.Tc), cell(p.KL.Pc), cell(p.RD.Pc), cell(p.Twu.Pc), cell(p.KL.w), cell(p.RD.w), cell(p.Twu.w)]), note: `Applied multipliers on the plus fraction: Tc × ${F(opts.tcMult, 4)}, Pc × ${F(opts.pcMult, 4)}, ω × ${F(opts.wMult, 4)}, kij scale ${F(opts.kijScale, 3)}, viscosity Vc × ${F(opts.vcMult, 3)}. Riazi–Daubert is extrapolated beyond its fitted range for cuts heavier than about C25.` });
  { const ids = f.comps.map((c) => c.id), rows = []; f.comps.forEach((a, i) => f.comps.forEach((b, j) => { if (j > i && f.kij[i][j] !== 0) rows.push([a.id, b.id, cell(f.kij[i][j], 4)]); })); if (rows.length) tables.push({ title: 'Binary interaction parameters in use (non-zero pairs)', columns: ['Component i', 'Component j', 'kij'], rows, note: `Kernel defaults × scale ${F(opts.kijScale, 3)}${nKij ? `, with ${nKij} user override(s)` : ''}. Component ids: ${ids.join(', ')}.` }); }
  tables.push({ title: 'Equilibrium states at the key conditions', columns: ['Condition', 'P (bara)', 'T (°C)', 'State', 'Vapour (mol %)', 'Gas (mass %)', 'ρ gas (kg/m³)', 'ρ oil (kg/m³)', 'μ gas (mPa·s)', 'μ oil (mPa·s)', 'Z gas', 'cp gas (J/kg/K)', 'cp oil (J/kg/K)', 'k gas (W/m/K)', 'k oil (W/m/K)', 'JT gas (K/bar)', 'c gas (m/s)', 'c oil (m/s)', 'κT oil (1/bar)', 'σ (mN/m)', 'Tangent-plane distance', 'Feed stable'],
    rows: conds.map((c) => { const g = c.s.gas, o = c.s.oil; return [c.name, cell(c.P), cell(c.T), phaseName[c.s.phase], cell(100 * c.s.beta), cell(100 * c.s.wG), g ? cell(g.rho) : '—', o ? cell(o.rho) : '—', g ? cell(g.mu * 1e3) : '—', o ? cell(o.mu * 1e3) : '—', g ? cell(g.Z) : '—', g ? cell(g.cp) : '—', o ? cell(o.cp) : '—', g ? cell(g.k) : '—', o ? cell(o.k) : '—', g ? cell(g.jtBar) : '—', g ? cell(g.sound) : '—', o ? cell(o.sound) : '—', o ? cell(o.kT) : '—', c.s.sigma !== null ? cell(c.s.sigma * 1e3) : '—', cell(c.st.tpd, 3), c.st.stable ? 'yes' : 'no (splits)']; }), note: 'Michelsen tangent-plane test: a negative distance means the single-phase feed is unstable and the Gibbs energy is lowered by splitting into two phases.' });
  tables.push({ title: 'Saturation, critical point and standard-condition yields', columns: ['Quantity', 'Value', 'Unit'], rows: [
    [`Saturation pressure at ${F(tRes)} °C (${psatType || 'none'})`, cell(psat), 'bara'], ['Critical temperature', cell(env.critical?.T), '°C'], ['Critical pressure', cell(env.critical?.P), 'bara'], ['Cricondenbar', cell(env.cricondenbar?.P), 'bara'], ['Temperature at the cricondenbar', cell(env.cricondenbar?.T), '°C'], ['Cricondentherm', cell(env.cricondentherm?.T), '°C'], ['Pressure at the cricondentherm', cell(env.cricondentherm?.P), 'bara'],
    ['Single-stage GOR', cell(rates.gor), 'Sm³/Sm³'], ['Multi-stage separator GOR', cell(sepGor), 'Sm³/Sm³'], ['Stock-tank gravity, single stage', cell(std.api), '°API'], ['Stock-tank gravity, separator train', cell(sep.api), '°API'], ['Stock-tank oil density', cell(hasOil ? std.rhoOilStd : null), 'kg/m³'], ['Gas gravity (air = 1)', cell(std.gasSG), '–'],
    ['Formation-volume factor at reservoir conditions (single-stage basis)', cell(bo), 'm³/Sm³'], ['Formation-volume factor at Psat (separator basis)', cell(boSat), 'm³/Sm³'], ['Gas formation-volume factor at reservoir conditions', cell(bg), 'm³/Sm³'],
    ['Temperature after isenthalpic expansion to arrival pressure', cell(jtOut?.T), '°C'], ['Temperature after isentropic expansion to arrival pressure', cell(isOut?.T), '°C']] });
  if (cce.length) tables.push({ title: `Constant composition expansion at ${F(tRes)} °C`, columns: ['P (bara)', 'Relative volume', 'Liquid (% of Vsat)', 'Y-function', 'Z', 'Density (kg/m³)', 'Compressibility (1/bar)', 'Vapour (mol %)'], rows: cce.map((r) => [cell(r.P), cell(r.vRel, 5), cell(r.liqPct), cell(r.Y), cell(r.Z), cell(r.rho), cell(r.co, 3), cell(100 * r.beta)]) });
  if (dle.rows.length) tables.push({ title: `Differential liberation at ${F(tRes)} °C`, columns: ['P (bara)', 'Rs (Sm³/Sm³)', 'Bo (m³/Sm³)', 'Bt (m³/Sm³)', 'Oil density (kg/m³)', 'Gas gravity', 'Gas Z', 'Bg (m³/Sm³)'], rows: dle.rows.map((r) => [cell(r.P), cell(r.rs), cell(r.bo, 5), cell(r.bt, 5), cell(r.rhoO), cell(r.gasSG), cell(r.Z), cell(r.bg)]), note: `Residual oil at standard conditions: ${F(dle.residual.rho)} kg/m³ (${F(dle.residual.api)} °API). Volumes are per volume of residual oil.` });
  if (cvd.length > 1) tables.push({ title: `Constant volume depletion at ${F(tRes)} °C`, columns: ['P (bara)', 'Liquid (% of cell volume)', 'Cumulative produced (mol %)', 'Z of produced gas', 'Two-phase Z', 'C7+ in produced gas (mol %)'], rows: cvd.map((r) => [cell(r.P), cell(r.liqPct), cell(r.produced), cell(r.Zgas), cell(r.Z2), cell(r.heavy)]) });
  tables.push({ title: 'Multi-stage separator test', columns: ['Stage', 'P (bara)', 'T (°C)', 'GOR (Sm³/Sm³ STO)', 'Gas gravity', 'Liquid leaving (mol per mol feed)'], rows: [...sep.stages.map((s, i) => [i === sep.stages.length - 1 ? 'Stock tank' : `Separator ${i + 1}`, cell(s.p), cell(s.t), cell(s.gor), cell(s.gasSG), cell(s.liqMol)]), ['Total', '—', '—', cell(sep.gor), cell(sep.gasSG), '—']], note: sep.vSto > 0 ? `Stock-tank oil ${F(sep.rhoSto)} kg/m³ (${F(sep.api)} °API).` : 'No stock-tank liquid is formed: the fluid is a dry gas at these conditions.' });
  if (swell.length) tables.push({ title: `Swelling test with ${INJECTION_GASES[v.swellGas].name.toLowerCase()}`, columns: ['Gas added (mol/mol)', 'Saturation pressure (bara)', 'Type', 'Swelling factor'], rows: swell.map((r) => [cell(r.inj), cell(r.psat), r.type || '—', cell(r.swell, 5)]) });
  if (bo0) tables.push({ title: 'Black-oil correlations against the equation of state', columns: ['Quantity', 'EOS', 'Standing', 'Vasquez–Beggs', 'Glasø', 'Unit'], rows: [
    ['Bubble-point pressure', cell(psatType === 'bubble' ? psat : null), cell(bo0.pbStanding), cell(bo0.pbVB), cell(bo0.pbGlaso), 'bara'], ['Oil formation-volume factor at the bubble point', cell(boSat), cell(bo0.boStanding), cell(bo0.boVB), cell(bo0.boGlaso), 'm³/Sm³'],
    ['Saturated-oil viscosity (Beggs–Robinson)', cell(satState?.oil ? satState.oil.mu * 1e3 : null), cell(bo0.muSat * 1e3), '—', '—', 'mPa·s'], ['Oil viscosity at reservoir pressure (Vasquez–Beggs)', cell(res.oil ? res.oil.mu * 1e3 : null), '—', cell(bo0.muP * 1e3), '—', 'mPa·s'], ['Dead-oil viscosity at reservoir temperature (Beggs–Robinson)', '—', cell(bo0.muDead * 1e3), '—', '—', 'mPa·s'],
    [`Solution GOR at ${F(pRef)} bara (Standing)`, cell(dle.rows.length ? interp1(dle.rows.map((r) => r.P).reverse(), dle.rows.map((r) => r.rs).reverse(), pRef) : null), cell(Math.min(rsStanding(pRef, sep.gasSG, sep.api, tRes), sepGor ?? Infinity)), '—', '—', 'Sm³/Sm³']], note: `Inputs to the correlations: separator GOR ${F(sepGor)} Sm³/Sm³, gas gravity ${F(sep.gasSG, 3)}, ${F(sep.api)} °API, ${F(tRes)} °C. The correlations are an independent sanity check, not a replacement for the compositional model.` });
  tables.push({ title: `Alternative property models (gas at ${refGas.name.toLowerCase()} conditions)`, columns: ['Property', 'Kernel model', 'Value', 'Alternative model', 'Value', 'Difference (%)'], rows: [
    ['Compressibility factor', `${f.eosId} with volume translation`, cell(gq.Z), 'Lee–Kesler corresponding states (Kay mixing)', cell(zLK), zLK ? cell(100 * (zLK / gq.Z - 1), 3) : '—'],
    ['Viscosity (mPa·s)', 'Lohrenz–Bray–Clark', cell(gq.mu * 1e3), 'Lee–Gonzalez–Eakin', cell(muLGE * 1e3), cell(100 * (muLGE / gq.mu - 1), 3)],
    ['Thermal conductivity (W/m/K)', 'Modified Eucken', cell(gq.k), 'Stiel–Thodos corresponding states', cell(kST), cell(100 * (kST / gq.k - 1), 3)]] });
  const setName = (m) => (m === 'vdwp' ? 'Munck constants' : 'Kihara cell potential'), hydSel = v.hydModel === 'vdwp' ? 'vdwp' : 'kihara', hydOther = hydSel === 'vdwp' ? 'kihara' : 'vdwp';
  tables.push({ title: 'Hydrate dissociation curve', columns: ['P (bara)', 'T fresh water (°C)', 'T case aqueous phase, published (°C)', 'T case aqueous phase, model before the allowance (°C)', 'Structure', 'Gas-gravity screening correlation (°C)', ...(hydX ? ['T custom aqueous phase (°C)'] : [])], rows: hP.map((P, i) => [cell(P), cell(hyd0.T[i]), cell(hydC.T[i]), cell((hydC.raw || hydC.T)[i]), hydC.st[i], cell(corr(P)), ...(hydX ? [cell(hydX.T[i])] : [])]), note: `${hydModelUsed}. ${vdw ? `Langmuir constants: ${setName(hydSel)}; guest fugacities from the equation of state; water activity from the selected activity and electrolyte models and from dissolved gas. The published curve credits ${F(100 * (1 - inhibitorAllowance(wCase)))} % of the inhibitor depression of the model (allowance ${F(allowance, 2)} °C at ${F(pRef)} bara): nothing is withheld up to 30 wt % inhibitor, 12 % at 50 wt % and above — the largest over-prediction found against measured data — and 0.004 K per bar is added above 250 bara, where the model drifts cold against measured methane, nitrogen and natural-gas points.` : ''}` });
  if (vdw) { const cmp = [25, 50, 100, 200, 300, 500].map((P) => { const a = hydAt(hyd0.raw, P), r = hydrateTofP(f, P, 1, { table: hydrateTable(hydOther, tx.hyd), guess: a }); return [P, cell(a), r ? cell(r.T) : '—', r ? cell(a - r.T, 3) : '—', cell(corr(P)), cell(corr(P) - a, 3)]; });
    tables.push({ title: 'Hydrate constants sets compared (fresh water)', columns: ['P (bara)', `${setName(hydSel)} — selected (°C)`, `${setName(hydOther)} (°C)`, 'Selected − other (°C)', 'Kernel screening correlation (°C)', 'Screening − selected (°C)'], rows: cmp, note: 'Against 309 measured dissociation points (1–960 bar; pure guests, mixtures, natural gases) the Kihara set has a bias of +0.7 K below 100 bar and +0.3 K above (warm), the Munck set −0.3 K below 100 bar, −1.2 K above and −2.1 K above 300 bar (cold: hydrate forms warmer than predicted); the Kihara set is therefore the default. The kernel screening correlation is what the other suites use before this suite has been run: it is fitted to the Kihara model for natural gases and is meant to lie on the warm side.' }); }
  tables.push({ title: `Inhibitor response at ${F(pRef)} bara (${inhName})`, columns: ['Inhibitor (wt %)', 'Depression, this model (°C)', 'Hammerschmidt + salt (°C)', 'Nielsen–Bucklin (°C)', 'Hydrate temperature (°C)'], rows: doseRows.map((r) => [r.w, cell(r.model), cell(r.hamm), cell(r.nb), cell(tHyd0 - r.model)]) });
  tables.push({ title: `Sensitivity of the ${inhName} requirement`, columns: ['Case', 'P (bara)', 'Seabed T (°C)', 'Salinity (wt %)', 'Margin (°C)', 'Required (wt %)', 'Change (wt %)'], rows: sens.map((r) => [r.name, cell(r.P), cell(r.tS), cell(r.S / 10), cell(r.mg), cell(r.w), cell(r.w - inhReq, 3)]) });
  tables.push({ title: 'Water distribution between gas and aqueous phase (case rates)', columns: ['Condition', 'P (bara)', 'T (°C)', 'Saturated water content of gas (mg/Sm³)', 'Produced water (kg/d)', 'Carried in the gas (kg/d)', 'Free aqueous water (kg/d)', 'Gas capacity (kg/d)'], rows: waterDist.map((r) => [r.name, cell(r.P), cell(r.T), cell(r.wc * 1e6), cell(r.tot), cell(r.inGas), cell(r.free), cell(r.satCap)]), note: 'Water dissolved in the liquid hydrocarbon is neglected. When the gas capacity exceeds the produced water, the line carries no free water at that point and hydrates can only form from condensed water further downstream.' });
  tables.push({ title: 'Aqueous phase', columns: ['Quantity', 'Value', 'Unit'], rows: [['Salinity (NaCl equivalent)', cell(aqCase.S / 10), 'wt %'], ['Salt molality', cell(awCase.molality), 'mol/kg water'], ['Inhibitor', aqCase.inhId === 'none' ? 'none' : INHIBITORS[aqCase.inhId].name, ''], ['Inhibitor concentration', cell(aqCase.inhWt), 'wt %'], ['Water activity from the electrolyte model', cell(awCase.awSalt, 5), '–'], ['Water activity from the inhibitor model', cell(awCase.awInh, 5), '–'], ['Water activity coefficient in the inhibitor solution', cell(awCase.gammaW, 5), '–'], ['Water activity of the aqueous phase', cell(awCase.aw, 5), '–'],
    [`Dissolved gas at ${F(pRef)} bara, ${F(tSea)} °C (Henry's law)`, cell(sol.sm3), 'Sm³/m³ water'], ...Object.entries(sol.x).map(([id, x]) => [`  mole fraction of ${COMP_LABELS[id]} in water`, cell(x, 3), '–']), [`Water content of gas at arrival (Bukacek × activity)`, cell(wcBuk * 1e6), 'mg/Sm³'], ['Water content of gas at arrival (Raoult + Poynting)', cell(wcRaoult * 1e6), 'mg/Sm³']], note: `Activity model ${models.act}; electrolyte model ${models.elec}. The salt molality refers to the water of the brine; the electrolyte and inhibitor contributions to ln a_w are added.` });
  if (hasLiq) tables.push({ title: 'Wax and asphaltene screening', columns: ['Quantity', 'Value', 'Unit'], rows: [['Wax appearance temperature, stock-tank oil', cell(wat), '°C'], [`Wax appearance temperature, live oil at ${F(pRef)} bara`, cell(watLive), '°C'], ['Wax-forming C18+ content', cell(waxDead.waxContent), 'wt %'], [`Solid wax at ${F(tSea)} °C`, cell(waxAtSea), 'wt %'], ['Laboratory wax appearance temperature', lab.wat ? cell(lab.wat.m) : '—', '°C'],
    ['In-situ oil density', asph ? cell(resOil.rho) : '—', 'kg/m³'], ['Undersaturation Pres − Pbub', cell(asph?.dP), 'bar'], ['de Boer limit: slight problems above', cell(asph?.limitSlight), 'bar'], ['de Boer limit: severe problems above', cell(asph?.limitSevere), 'bar'], ['de Boer class', asph?.deBoer || '—', ''], ['Colloidal instability index', cell(asph?.cii), '–'], ['Asphaltene screening risk', asph?.risk || '—', '']], note: `Wax: Won melting properties, Pedersen wax-forming fraction, ideal solid and liquid solutions, detection limit ${detect} wt %. de Boer boundaries are an approximate digitisation of the published plot — a screening, not a prediction of onset pressure.` });
  if (v.useLab !== false) {
    const rows = [];
    for (const r of lab.psat || []) rows.push(['Saturation pressure (bara)', `${F(r.t)} °C`, cell(r.m), cell(r.c), r.c ? cell(100 * (r.c / r.m - 1), 3) : '—']);
    for (const r of lab.cce || []) rows.push(['CCE relative volume', `${F(r.p)} bara`, cell(r.m, 5), cell(r.c, 5), cell(100 * (r.c / r.m - 1), 3)]);
    for (const r of lab.dle || []) rows.push(['DLE Rs (Sm³/Sm³)', `${F(r.p)} bara`, cell(r.rs), cell(r.c.rs), cell(100 * (r.c.rs / r.rs - 1), 3)], ['DLE Bo (m³/Sm³)', `${F(r.p)} bara`, cell(r.bo, 5), cell(r.c.bo, 5), cell(100 * (r.c.bo / r.bo - 1), 3)], ['DLE oil density (kg/m³)', `${F(r.p)} bara`, cell(r.rho), cell(r.c.rhoO), cell(100 * (r.c.rhoO / r.rho - 1), 3)]);
    for (const r of lab.visc || []) rows.push(['Oil viscosity (mPa·s)', `${F(r.p)} bara`, cell(r.m), cell(r.c), cell(100 * (r.c / r.m - 1), 3)]);
    for (const r of lab.hyd || []) rows.push(['Hydrate dissociation temperature, fresh water (°C)', `${F(r.p)} bara`, cell(r.m), cell(r.c), cell(r.c - r.m, 3) + ' °C']);
    if (lab.wat) rows.push(['Wax appearance temperature (°C)', 'stock-tank oil', cell(lab.wat.m), cell(lab.wat.c), cell(lab.wat.c - lab.wat.m, 3) + ' °C']);
    if (rows.length) tables.push({ title: 'Model against the laboratory report', columns: ['Measurement', 'Condition', 'Laboratory', 'Model', 'Deviation (%)'], rows, note: `Average absolute deviations: ${Object.entries(lab.aad).filter(([, x]) => x !== null).map(([k, x]) => `${k} ${F(x, 2)} %`).join(', ') || '—'}.` });
  }

  // ---- KPIs --------------------------------------------------------------------------------------------------------------------
  const kpis = [
    { label: psatType === 'dew' ? 'Dew point' : 'Bubble point', value: psat !== null ? rd(psat) : '—', unit: 'bara', status: psat !== null && pRes < psat ? 'warn' : 'ok', help: `Saturation pressure at ${F(tRes)} °C` },
    hasOil || rates.gor === null ? { label: 'Gas–oil ratio', value: rates.gor !== null ? rd(rates.gor) : '—', unit: 'Sm³/Sm³', help: 'Single-stage flash to standard conditions' } : { label: 'Condensate yield', value: rd(1e6 / rates.gor), unit: 'Sm³/MSm³', help: 'Stock-tank liquid per million Sm³ of gas, single-stage flash' },
    hasLiq ? { label: 'Stock-tank gravity', value: std.api !== null ? rd(std.api) : '—', unit: '°API' } : { label: 'Gas gravity', value: rd(gasSG), unit: 'air = 1' },
    bo !== null ? { label: 'Formation-volume factor', value: rd(bo), unit: 'm³/Sm³', help: 'Reservoir volume per standard volume of stock-tank oil' } : { label: 'Gas formation-volume factor', value: bg !== null ? rd(bg) : '—', unit: 'm³/Sm³', help: 'Reservoir volume per standard volume of gas' },
    { label: 'Critical point', value: env.critical ? `${F(env.critical.T, 0)} °C / ${F(env.critical.P, 0)} bara` : '—' },
    { label: 'Cricondenbar', value: env.cricondenbar ? rd(env.cricondenbar.P) : '—', unit: 'bara' },
    { label: 'Cricondentherm', value: env.cricondentherm ? rd(env.cricondentherm.T) : '—', unit: '°C' },
    { label: 'Hydrate temperature', value: rd(tHyd), unit: '°C', status: subcool > 3 ? 'bad' : subcool > 0 ? 'warn' : 'ok', help: `At ${F(pRef)} bara with the case aqueous phase (${structure})` },
    { label: 'Subcooling at seabed', value: rd(subcool, 3), unit: '°C', status: subcool > 3 ? 'bad' : subcool > 0 ? 'warn' : 'ok', help: 'Hydrate temperature minus seabed temperature at the reference pressure' },
    { label: `${inhDesign} required`, value: rd(inhReq, 3), unit: 'wt %', status: inhReq > 45 ? 'bad' : inhReq > 0 ? 'warn' : 'ok', help: `For a ${F(margin)} °C margin at the cold reference point` },
    { label: 'Wax appearance', value: wat !== null ? rd(wat, 3) : '—', unit: '°C', status: wat !== null && wat > tSea ? 'warn' : 'ok' },
    { label: 'Oil viscosity, reservoir', value: res.oil ? rd(res.oil.mu * 1e3, 3) : '—', unit: 'mPa·s' },
    { label: 'Water in gas at arrival', value: rd(wcBuk * 1e6, 3), unit: 'mg/Sm³' },
    { label: 'Temperature after choke', value: jtOut ? rd(jtOut.T, 3) : '—', unit: '°C', help: `Isenthalpic expansion from reservoir conditions to ${F(pArr)} bara` },
  ];

  // ---- outputs -------------------------------------------------------------------------------------------------------------------
  const outputs = {
    table, eos: f.eosId, gor: fin(rates.gor), api: fin(std.api), gasSG: fin(gasSG), mwGas: fin(std.mwGas), rhoOilStd: hasOil ? fin(std.rhoOilStd) : null, psat, psatType, bo: fin(bo),
    envelope: { T: env.T, P: env.P, type: env.type, lowT: env.lowT, lowP: env.lowP }, critical: env.critical ? { T: env.critical.T, P: env.critical.P } : null, cricondenbar: fin(env.cricondenbar?.P), cricondentherm: fin(env.cricondentherm?.T),
    hydrateCurve: { P: hP, T0: hyd0.T, T: hydC.T }, hydrateStructure: structure, hydrateDepression: fin(depression, 0), waterContent: fin(wcBuk, 0), wat, waxContent: fin(waxDead.waxContent, 0), inhibitorWt: fin(inhReq, 0),
    components: f.comps.map((c) => ({ id: c.id, z: c.z, MW: c.MW, Tc: c.Tc, Pc: c.Pc, w: c.w })), rates: { mHC: rates.mHC, mW: rates.mW, qOilStd: rates.qOilStd, qGasStd: rates.qGasStd, qWaterStd: rates.qWaterStd },
    // extras
    fluidName: spec.name || 'Case fluid', isCaseFluid: !library, modelOptions: { ...opts }, hydrateModel: hydModelUsed, hydrateT: fin(tHyd), hydrateT0: fin(tHyd0), subcooling: fin(subcool), hydratePressureAtSeabed: pDissSea, inhibitorId: inhDesign, waterActivity: fin(awCase.aw),
    inhibitorSensitivity: sens.map((r) => ({ name: r.name, wt: r.w })), waterDistribution: waterDist.map((r) => ({ name: r.name, inGas: r.inGas, free: r.free })), margin, watLive, waxAtSeabed: fin(waxAtSea, 0), waxCurve: { T: waxT, wt: waxT.map((t) => waxDead.solidWt(t)) }, asphalteneRisk: asph?.risk || null, sepGor: fin(sepGor), sepApi: fin(sep.api), bg: fin(bg), jtOutletT: jtOut ? jtOut.T : null, dissolvedGas: fin(sol.sm3, 0),
    qualityLines: env.quality, labAad: Object.fromEntries(Object.entries(lab.aad).map(([k, x]) => [k, fin(x)])),
  };
  // ---- alternative equations of state, association, solid phases, miscibility, initial solid inventory -------------------------
  prog(0.93, 'Alternative models and solid-phase equilibria');
  const altOn = v.altEos !== 'none', altFull = v.altEos === 'all', pct = (a, b) => (Number.isFinite(a) && Number.isFinite(b) && b !== 0 ? cell(100 * (a / b - 1), 3) : '—');
  const alt = { rows: [], zLine: null };
  if (altOn) {
    for (const c of altFull ? [conds[0], conds[2], conds[3]] : [conds[0], conds[3]]) for (const [kindP, st, xx] of [['vapour', c.s.gas, c.s.y], ['liquid', c.s.oil, c.s.x]]) {
      if (!st) continue;
      const a = altPhase(f, xx, c.P, c.T, kindP);
      alt.rows.push({ name: c.name, phase: kindP === 'vapour' ? 'gas' : 'oil', P: c.P, T: c.T, rho: st.rho, sound: st.sound, cp: st.cp, mu: st.mu, a });
    }
    // gas compressibility factor of the separator gas along pressure by every model (plot)
    if (altFull && std.beta > 0 && std.y) { const fg = withZ(f, std.y), Ps = logspace(5, Math.max(pRes, 200), 9), zs = { pr: [], pc: [], ge: [], lk: [] }, tz = Math.max(tArr, 15);
      for (const P of Ps) { const k = phaseProps(fg, fg.z, P, tz, 'vapour', { thermal: false }), a = altPhase(fg, fg.z, P, tz, 'vapour', { visc: false, thermal: false }), zOf = (r) => (r ? (P * 1e5 * k.MW * 1e-3) / (r.rho * R * (tz + KEL)) : NaN); zs.pr.push(k.Z); zs.pc.push(zOf(a.pcsaft)); zs.ge.push(zOf(a.gerg)); zs.lk.push(a.lk ? a.lk.Z : NaN); }
      alt.zLine = { P: Ps, T: tz, ...zs }; }
  }
  const altOil = alt.rows.find((r) => r.phase === 'oil' && r.name === 'Reservoir') || alt.rows.find((r) => r.phase === 'oil') || null, altGas = alt.rows.find((r) => r.phase === 'gas') || null;
  // water in gas and inhibitor partitioning by cubic-plus-association
  const cpaRows = [];
  if (altOn) for (const c of [conds[3]]) { if (!c.s.gas) continue; try { const r = cpaWater(f, c.s.y, c.P, c.T, { inhId: aqCase.inhId, inhWt: aqCase.inhWt, awSalt: waterActivity({ S: aqCase.S }, c.T + KEL, models, tx).awSalt }); if (r.ok && Number.isFinite(r.wc)) cpaRows.push({ name: c.name, P: c.P, T: c.T, cpa: r.wc * tx.wc, buk: waterContent(c.P, c.T) * awOf(aqCase)(c.T + KEL) * tx.wc, raoult: waterContentRaoult(c.P, c.T, awOf(aqCase)(c.T + KEL)), inhLoss: r.inhLoss, xGas: r.xGas }); } catch { /* CPA has no root at this state */ } }
  // hydrate: second constants set, structure H and the water-side Gibbs minimum at the cold reference point
  const solids = ['all', 'hydrate', 'waxAsph', 'none'].includes(v.solidsAllowed) ? v.solidsAllowed : 'all', allowHyd = solids === 'all' || solids === 'hydrate', allowWax = solids === 'all' || solids === 'waxAsph';
  const hydAlt = vdw ? (() => { const other = v.hydModel === 'vdwp' ? 'kihara' : 'vdwp', r = hydrateTofP(f, pRef, 1, { table: hydrateTable(other, tx.hyd), guess: tHyd0 }); return r ? { model: other === 'kihara' ? 'Kihara cell potential' : 'Munck constants', T: r.T, structure: r.structure } : null; })() : null;
  const formersZ = f.comps.reduce((q, c, i) => q + (LANGMUIR[c.id] ? f.z[i] : 0), 0), nWater = rates.mW / (MW_W * 1e-3), gasPerWater = nWater > 0 ? (table.rates.nHC * formersZ) / nWater : Infinity;
  const wEq = vdw ? waterPhaseEquilibrium(f, pRef, tSea, aqCase, models, { gasPerWater, allowHydrate: allowHyd, tune: tx }) : null;
  const hyd0In = clamp(num(v.hyd0, 0), 0, 100) / 100, solid0In = clamp(num(v.solid0, 0), 0, 60);
  const gib = [conds[0], conds[2], conds[3]].map((c) => ({ name: c.name, ...gibbsSplit(f, c.P, c.T) }));
  // asphaltene solid–liquid equilibrium along the depletion at reservoir temperature
  const aspOpt = { wAspSto: clamp(num(v.saraAsp, 0), 0, 60) / 100, aspMW: clamp(num(v.aspMW, 750), 300, 5000), aspRho: clamp(num(v.aspRho, 1100), 900, 1300), aspDelta: clamp(num(v.aspDelta, 21), 15, 30), aspDeltaT: clamp(num(v.aspDeltaT, 1.07e-3), 0, 3e-3) };
  const aspP = hasOil ? [...linspace(Math.max(pRes, (psat || pRes) * 1.2), Math.max((psat || 50) * 0.25, 5), 10)] : [], asp2 = hasOil && aspOpt.wAspSto > 0 && allowWax ? asphalteneCurve(f, tRes, aspP, aspOpt) : null;
  // multiple-contact miscibility with the swelling gas
  const gasInj = INJECTION_GASES[v.swellGas], mmp = gasInj && hasOil && psatType === 'bubble' && v.mmpCalc === true ? miscibilityPressure(f, gasInj.comp, tRes) : null;
  // initial solid inventory against the equilibrium at the cold reference point
  const waxEq = allowWax ? waxAtSea : 0, initRows = [
    ['Hydrate (fraction of the water inventory)', cell(100 * hyd0In), wEq ? cell(100 * wEq.conversion) : '—', '%', !allowHyd ? 'hydrate phase not permitted: aqueous phase kept metastable' : wEq ? (hyd0In > wEq.conversion + 1e-9 ? `${F(100 * (hyd0In - wEq.conversion))} % of the water dissociates from hydrate` : hyd0In < wEq.conversion - 1e-9 ? `${F(100 * (wEq.conversion - hyd0In))} % of the water can still convert` : 'at equilibrium') : 'gas-gravity correlation: no phase amounts'],
    ['Wax solids in the stock-tank oil', cell(solid0In), cell(waxEq), 'wt %', !allowWax ? 'solid hydrocarbon phases not permitted' : solid0In > waxEq + 1e-9 ? `${F(solid0In - waxEq, 2)} wt % re-dissolves at ${F(tSea)} °C` : `${F(waxEq - solid0In, 2)} wt % can still precipitate at ${F(tSea)} °C`],
    ['Asphaltene precipitated at reservoir temperature (maximum over pressure)', '0', asp2 ? cell(asp2.maxPrecip) : '—', 'wt % of liquid', asp2 ? (asp2.maxPrecip > 0 ? `onset between ${asp2.lowerOnset !== null ? F(asp2.lowerOnset) : '—'} and ${asp2.upperOnset !== null ? F(asp2.upperOnset) : '—'} bara` : 'no precipitation predicted') : allowWax ? 'no asphaltene content or no liquid' : 'not permitted'],
  ];
  // additional measurements (comparison path for data types without a bundled reference set)
  const xl = [];
  if (v.useLab !== false) {
    const gasAt = (P, T) => { const s = props(f, P, T, { thermal: true }); return { s, g: s.phase !== 'oil' ? derivedProps(f, s.y, P, T, 'vapour') : null, o: s.phase !== 'gas' ? derivedProps(f, s.x, P, T, 'liquid') : null }; };
    const Q = {
      rhoG: ['Gas density', 'kg/m³', (P, T) => gasAt(P, T).g?.rho], zG: ['Gas compressibility factor', '–', (P, T) => gasAt(P, T).g?.Z], muG: ['Gas viscosity', 'mPa·s', (P, T) => { const g = gasAt(P, T).g; return g ? g.mu * 1e3 * tx.muG : NaN; }], cpG: ['Gas heat capacity', 'kJ/kg/K', (P, T) => gasAt(P, T).g?.cp / 1000],
      cpO: ['Oil heat capacity', 'kJ/kg/K', (P, T) => gasAt(P, T).o?.cp / 1000], soundG: ['Speed of sound in the gas', 'm/s', (P, T) => gasAt(P, T).g?.sound], soundO: ['Speed of sound in the oil', 'm/s', (P, T) => gasAt(P, T).o?.sound], co: ['Oil isothermal compressibility', '1/bar', (P, T) => gasAt(P, T).o?.kT],
      sigma: ['Gas–oil interfacial tension', 'mN/m', (P, T) => { const s = gasAt(P, T).s; return s.phase === 'two' ? s.sigma * 1e3 : NaN; }], wc: ['Water content of the gas', 'mg/Sm³', (P, T) => waterContent(P, T) * awOf(aqCase)(T + KEL) * tx.wc * 1e6], rhoW: ['Aqueous-phase density', 'kg/m³', (P, T) => waterProps(P, T, aqK).rho * tx.rhoW],
      psat: ['Saturation pressure', 'bara', (P, T) => saturationPoint(f, T, { seed: seedFor(env, T) }).P], hydT: ['Hydrate temperature (case aqueous phase)', '°C', (P) => { const r = vdw ? hydrateTofP(f, P, awOf(aqCase), { molality: molOf(aqCase) }) : null; return r ? r.T : corr(P) - dep; }], jt: ['Joule–Thomson coefficient of the gas', 'K/bar', (P, T) => gasAt(P, T).g?.jtBar],
    };
    for (const r of Array.isArray(v.labProps) ? v.labProps : []) { const q = Q[String(r.q || '').trim()], P = +r.p, T = +r.t, m = +r.value; if (!q || !isNum(P) || !isNum(T) || !isNum(m) || !(P > 0)) continue; let c = NaN; try { c = +q[2](clamp(P, 1.02, 1400), clamp(T, -40, 250)); } catch { c = NaN; } xl.push({ name: q[0], unit: q[1], P, T, m, c }); }
    const lm = num(v.labMmp, 0); if (lm > 0 && mmp && mmp.mmp !== null) xl.push({ name: 'Slim-tube minimum miscibility pressure', unit: 'bara', P: lm, T: tRes, m: lm, c: mmp.mmp });
    const la = num(v.labAop, 0); if (la > 0 && asp2) xl.push({ name: 'Upper asphaltene onset pressure', unit: 'bara', P: la, T: tRes, m: la, c: asp2.upperOnset ?? NaN });
  }
  if (alt.rows.length) {
    tables.push({ title: 'Alternative equations of state against the cubic model', columns: ['Condition', 'Phase', 'P (bara)', 'T (°C)', `ρ ${f.eosId} (kg/m³)`, 'ρ PC-SAFT', 'Δ (%)', 'ρ GERG-2008', 'Δ (%)', 'ρ Lee–Kesler', 'Δ (%)', `Sound ${f.eosId} (m/s)`, 'Sound PC-SAFT', 'Sound GERG-2008'], rows: alt.rows.map((r) => [r.name, r.phase, cell(r.P), cell(r.T), cell(r.rho), cell(r.a.pcsaft?.rho), pct(r.a.pcsaft?.rho, r.rho), cell(r.a.gerg?.rho), pct(r.a.gerg?.rho, r.rho), cell(r.a.lk?.rho), pct(r.a.lk?.rho, r.rho), cell(r.sound), cell(r.a.pcsaft?.w), cell(r.a.gerg?.w)]), note: 'Residual-Helmholtz models evaluated on the phase compositions of the cubic flash. PC-SAFT (Gross & Sadowski 2001): pure-component parameters of the light components, C7+ pseudo-components interpolated between the n-alkane sets by molar mass. GERG-2008 (reduced to the eleven kernel components): only for phases with at most 1 mol % C7+, counted as n-hexane; “—” otherwise. Lee–Kesler: modified Benedict–Webb–Rubin equation with Kay mixing, gas phase only. The published property table remains the cubic model.' });
    tables.push({ title: 'Viscosity models', columns: ['Condition', 'Phase', 'P (bara)', 'T (°C)', 'Lohrenz–Bray–Clark (mPa·s)', 'Pedersen corresponding states (mPa·s)', 'Δ (%)'], rows: alt.rows.map((r) => [r.name, r.phase, cell(r.P), cell(r.T), cell(r.mu * 1e3), cell(r.a.muPedersen !== null ? r.a.muPedersen * 1e3 : NaN), pct(r.a.muPedersen, r.mu)]), note: 'Pedersen et al. corresponding-states model with methane as the reference fluid (Hanley et al. viscosity correlation, GERG-2008 methane density). It needs no critical volumes; for oils heavier than about 30 °API it usually reads higher than an untuned Lohrenz–Bray–Clark model.' });
  }
  if (alt.zLine) plots.push({ type: 'line', title: `Compressibility factor of the separator gas at ${F(alt.zLine.T)} °C by four equations of state`, xlabel: 'Pressure (bara)', ylabel: 'Z (–)', logx: true, series: [{ name: f.eosId + ' (volume-translated)', x: alt.zLine.P, y: alt.zLine.pr }, ...[['PC-SAFT', alt.zLine.pc], ['GERG-2008', alt.zLine.ge], ['Lee–Kesler', alt.zLine.lk]].map(([name, y]) => { const xs = [], ys = []; y.forEach((q, i) => { if (Number.isFinite(q)) { xs.push(alt.zLine.P[i]); ys.push(q); } }); return { name, x: xs, y: ys, dash: true }; }).filter((q) => q.x.length > 1)], note: 'GERG-2008 is the reference-quality model for natural gas; the spread of the other curves around it is their model error.' });
  if (cpaRows.length) tables.push({ title: 'Water in gas and inhibitor partitioning by cubic-plus-association', columns: ['Condition', 'P (bara)', 'T (°C)', 'Water in gas, CPA (mg/Sm³)', 'Bukacek (mg/Sm³)', 'Raoult + Poynting (mg/Sm³)', 'Inhibitor lost to gas (mg/Sm³)', 'Gas dissolved in the aqueous phase (mol %)'], rows: cpaRows.map((r) => [r.name, cell(r.P), cell(r.T), cell(r.cpa * 1e6), cell(r.buk * 1e6), cell(r.raoult * 1e6), cell(r.inhLoss * 1e6), cell(100 * r.xGas)]), note: 'CPA (SRK plus Wertheim association; water and MEG four-site, methanol two-site) solves the gas–aqueous equilibrium rigorously, including the inhibitor that leaves with the gas. The salt enters through the water activity of the brine. The flow suites keep the Bukacek value.' });
  if (vdw) tables.push({ title: 'Hydrate models and stable water-side phases at the cold reference point', columns: ['Quantity', 'Value', 'Unit'], rows: [
    [`Hydrate temperature, fresh water (${v.hydModel === 'vdwp' ? 'Munck constants' : 'Kihara cell potential'})`, cell(tHyd0), '°C'], ...(hydAlt ? [[`Hydrate temperature, fresh water (${hydAlt.model})`, cell(hydAlt.T), '°C'], ['Difference between the two constants sets', cell(hydAlt.T - tHyd0, 3), '°C']] : []),
    ...(v.shFormer === true && refFresh ? [['Structure-H temperature with the heavy former present', refFresh.TsH !== null && refFresh.TsH !== undefined ? cell(refFresh.TsH) : '—', '°C']] : []),
    ['Driving force Δμ/RT of water before conversion (> 0: hydrate stable)', wEq ? cell(wEq.drive0) : '—', '–'], ['Stable water-side phases', wEq ? wEq.phases.join(' + ') : '—', ''], ['Water converted to hydrate at equilibrium', wEq ? cell(100 * wEq.conversion) : '—', '%'], ['Conversion limited by', wEq?.limitedBy || '—', ''],
    ['Gibbs energy released', wEq ? cell(wEq.dG) : '—', 'J per mol of water'], ['Salinity of the remaining brine', wEq ? cell(wEq.saltFinal / 10) : '—', 'wt %'], ['Inhibitor in the remaining aqueous phase', wEq ? cell(wEq.inhFinal) : '—', 'wt %'], ['Hydrate formed', wEq ? cell(wEq.hydrateKgPerKgWater) : '—', 'kg per kg of water'],
  ], note: 'The water-side Gibbs energy is minimised over aqueous liquid, ice and the hydrate structures: water converts until its chemical potential in the concentrating brine or inhibitor solution equals that in the hydrate, or the water or the hydrate formers run out. Formation kinetics are the solids suite\'s subject.' });
  tables.push({ title: 'Gibbs energy of the hydrocarbon split', columns: ['Condition', 'g/RT as one phase', 'g/RT at equilibrium', 'Δg/RT (≤ 0)', 'Tangent-plane distance'], rows: gib.map((r) => [r.name, cell(r.g1, 6), cell(r.g2, 6), cell(r.dg, 3), cell(r.tpd, 3)]), note: 'The flash is the minimum of the Gibbs energy: the split never has a higher Gibbs energy than the single phase, and a negative tangent-plane distance marks an unstable feed.' });
  if (asp2 && asp2.rows.length) {
    tables.push({ title: `Asphaltene solid–liquid equilibrium at ${F(tRes)} °C (Flory–Huggins)`, columns: ['P (bara)', 'Liquid solubility parameter (MPa^0.5)', 'Soluble limit (vol %)', 'Asphaltene in liquid (vol %)', 'Precipitated (wt % of liquid)'], rows: asp2.rows.map((r) => [cell(r.P), cell(r.deltaL), cell(100 * r.phiMax), cell(100 * r.phiA), cell(r.precipWt)]), note: `Hirschberg-type regular-solution model: liquid solubility parameter from the cohesive energy of the equation of state, asphaltene molar mass ${F(aspOpt.aspMW, 0)} g/mol, density ${F(aspOpt.aspRho, 0)} kg/m³ and solubility parameter ${F(aspOpt.aspDelta, 2)} MPa^0.5 at 25 °C (inputs; tune the solubility parameter to a measured onset pressure). ${asp2.maxPrecip > 0 ? `Precipitation is predicted between ${asp2.lowerOnset !== null ? F(asp2.lowerOnset) : 'the lowest pressure examined'} and ${asp2.upperOnset !== null ? F(asp2.upperOnset) : '—'} bara, at most ${F(asp2.maxPrecip, 2)} wt % near ${F(asp2.pAtMax)} bara.` : 'No precipitation is predicted over the pressure range.'}` });
    plots.push({ type: 'line', title: 'Asphaltene stability along the depletion', xlabel: 'Pressure (bara)', ylabel: 'Volume fraction (%)', series: [{ name: 'Soluble limit', x: asp2.rows.map((r) => r.P), y: asp2.rows.map((r) => Math.min(100 * r.phiMax, 100)) }, { name: 'Asphaltene in the liquid', x: asp2.rows.map((r) => r.P), y: asp2.rows.map((r) => 100 * r.phiA), dash: true }], vlines: psat !== null ? [{ x: psat, label: 'saturation' }] : [], note: 'Asphaltene precipitates where its content exceeds the soluble limit; the limit is lowest near the bubble point, where the liquid is lightest.' });
  }
  if (mmp && mmp.mmp !== null) tables.push({ title: `Miscibility of ${gasInj.name.toLowerCase()} with the reservoir fluid at ${F(tRes)} °C`, columns: ['Quantity', 'Value', 'Unit'], rows: [['Multiple-contact miscibility pressure', mmp.capped && !mmp.mechanism ? `> ${F(mmp.mmp, 0)}` : cell(mmp.mmp), 'bara'], ['Mechanism', mmp.mechanism || 'not miscible below the limit', ''], ['Reservoir pressure', cell(pRes), 'bara'], ['Displacement at reservoir pressure', mmp.mechanism && pRes >= mmp.mmp ? 'miscible' : 'immiscible', '']], note: 'Single-cell multiple-contact test with the equation of state (six forward and six backward contacts, pressure resolved to about 3 %): an estimate of the slim-tube minimum miscibility pressure. A slim-tube displacement with a combined condensing/vaporising drive can be miscible at a somewhat lower pressure.' });
  tables.push({ title: `Initial solid inventory against equilibrium at ${F(pRef)} bara, ${F(tSea)} °C`, columns: ['Solid phase', 'Initial', 'Equilibrium', 'Unit', 'Consequence'], rows: initRows, note: `Solid phases permitted: ${{ all: 'hydrate, ice, wax and asphaltene', hydrate: 'hydrate and ice only', waxAsph: 'wax and asphaltene only', none: 'none (fluid phases only)' }[solids]}. The initial amounts are inputs (for example the state after a shut-in); the equilibrium amounts are what the thermodynamic model allows at the cold reference point.` });
  if (xl.length) tables.push({ title: 'Additional measurements against the model', columns: ['Measurement', 'P (bara)', 'T (°C)', 'Measured', 'Model', 'Unit', 'Deviation (%)'], rows: xl.map((r) => [r.name, cell(r.P), cell(r.T), cell(r.m), cell(r.c), r.unit, pct(r.c, r.m)]), note: 'Comparison path for measurement types that have no bundled reference set (interfacial tension, water content, calorimetry, speed of sound, compressibility, slim tube, asphaltene onset): enter the points in the “Other measurements” table.' });
  if (v.shFormer === true && refFresh && refFresh.structure === 'sH') warnings.push({ level: 'warn', msg: `With a structure-H former in the liquid the hydrate is stable up to ${F(refFresh.T)} °C at ${F(pRef)} bara as structure H — warmer than structures I and II would allow.` });
  if (asp2 && asp2.maxPrecip > 0) warnings.push({ level: 'warn', msg: `The asphaltene solid-phase model predicts up to ${F(asp2.maxPrecip, 2)} wt % precipitation near ${F(asp2.pAtMax)} bara at ${F(tRes)} °C (upper onset ${asp2.upperOnset !== null ? F(asp2.upperOnset) + ' bara' : 'above the range'}); confirm with a measured onset pressure and tune the asphaltene solubility parameter.` });
  if (mmp && mmp.mechanism && pRes < mmp.mmp) recommendations.push(`${gasInj.name} is not multiple-contact miscible with the oil at ${F(pRes)} bara: the estimated miscibility pressure is ${F(mmp.mmp)} bara (${mmp.mechanism}); a slim-tube test would confirm it.`);
  if (wEq && wEq.conversion > 0 && wEq.conversion < 1) recommendations.push(`At the cold reference point ${F(100 * wEq.conversion)} % of the water can convert to hydrate before the remaining aqueous phase (${F(wEq.saltFinal / 10)} wt % salt, ${F(wEq.inhFinal)} wt % inhibitor) is self-inhibited: the plugging potential is bounded by that amount.`);
  if (altGas && altGas.a.gerg && Math.abs(altGas.a.gerg.rho / altGas.rho - 1) > 0.03) warnings.push({ level: 'info', msg: `The ${f.eosId} gas density at ${altGas.name.toLowerCase()} conditions differs from GERG-2008 by ${F(100 * (altGas.rho / altGas.a.gerg.rho - 1))} %: use the GERG value for fiscal or line-pack calculations.` });
  kpis.push({ label: 'Oil viscosity, Pedersen', value: altOil && altOil.a.muPedersen !== null ? rd(altOil.a.muPedersen * 1e3, 3) : '—', unit: 'mPa·s', help: 'Corresponding-states cross-check of the Lohrenz–Bray–Clark value' },
    { label: 'Asphaltene onset', value: asp2 && asp2.upperOnset !== null ? rd(asp2.upperOnset) : '—', unit: 'bara', status: asp2 && asp2.maxPrecip > 0 ? 'warn' : 'ok', help: 'Upper onset pressure of the Flory–Huggins solid-phase model at reservoir temperature' });
  Object.assign(outputs, {
    altEos: alt.rows.map((r) => ({ name: r.name, phase: r.phase, rho: fin(r.rho), rhoPcSaft: fin(r.a.pcsaft?.rho), rhoGerg: fin(r.a.gerg?.rho), rhoLeeKesler: fin(r.a.lk?.rho), muPedersen: fin(r.a.muPedersen) })),
    waterContentCpa: cpaRows.length ? fin(cpaRows[0].cpa) : null, inhibitorLossToGas: cpaRows.length ? fin(cpaRows[0].inhLoss) : null, hydrateConstants: v.hydModel === 'vdwp' ? 'munck' : 'kihara', hydrateTAlt: hydAlt ? fin(hydAlt.T) : null,
    hydrateConversion: wEq ? fin(wEq.conversion) : null, hydratePhases: wEq ? wEq.phases : null, hydrateGibbs: wEq ? fin(wEq.dG) : null, solidsPermitted: solids, initialHydrate: hyd0In, initialSolids: solid0In,
    asphalteneOnset: asp2 ? { upper: fin(asp2.upperOnset), lower: fin(asp2.lowerOnset), maxPrecip: fin(asp2.maxPrecip, 0) } : null, mmp: mmp ? fin(mmp.mmp) : null, mmpMechanism: mmp ? mmp.mechanism : null, calibrationMultipliers: { ...tx },
  });
  const kind = !hasOil || rates.gor > 50000 ? 'dry gas' : psatType === 'dew' || (psat === null && rates.gor > 900) ? (rates.gor > 10000 ? 'wet gas' : 'gas condensate') : rates.gor > 350 ? 'volatile oil' : 'black oil';
  const summary = `${spec.name || 'The case fluid'} behaves as a ${kind} (${f.eosId}, ${f.n} components): ${psat !== null ? `${psatType} point ${F(psat)} bara at ${F(tRes)} °C` : `no saturation pressure at ${F(tRes)} °C`}${rates.gor !== null && hasOil ? `, GOR ${F(rates.gor)} Sm³/Sm³, ${F(std.api)} °API` : `, gas gravity ${F(gasSG, 3)}`}; hydrate equilibrium ${F(tHyd)} °C at ${F(pRef)} bara (${subcool > 0 ? `${F(subcool)} °C subcooling at the seabed` : 'outside the hydrate region at the seabed'})${wat !== null ? `, wax appearance ${F(wat)} °C` : ''}.`;
  prog(1, 'Done');
  return { summary, kpis, warnings, recommendations, plots, tables, balances, outputs };
}

// ---- calibration model ----------------------------------------------------------------------------------------------------------
const calCache = new Map();
/**
 * Predictions for one laboratory point at (calP, calT): saturation pressure at calT, liquid density, viscosity and solution GOR, and —
 * evaluated only when a row asks for them — formation-volume factor, properties of the produced (stock-tank flash) gas, liquid heat
 * capacity and compressibility, interfacial tension, water content, aqueous density, water activity and hydrate temperature for the
 * aqueous phase of the row (calSal wt % NaCl, calInhWt wt % of the design inhibitor) and the multiple-contact miscibility pressure.
 */
export function calibrationModel(v) {
  const src = FLUID_LIBRARY[v.fluidSource] ? resolveSpec(null, v.fluidSource).spec : lastModel?.spec || mergeSpec(null), tune = tuneOf(v, src), tx = tuneX(v);
  const T = clamp(num(v.calT, num(v.tRes, 90)), -20, 250), P = clamp(num(v.calP, num(v.pRef, 100)), 1.02, 1400);
  const key = JSON.stringify([v.fluidSource, src.comp, src.c7MW, src.c7SG, tune, v.kijTable || null, tx.shift, tx.par, tx.cp]);
  let m = calCache.get(key);
  if (!m) { m = { f: tunedFluid(src, tune, v.kijTable, tx).f, sat: new Map() }; if (calCache.size > 400) calCache.clear(); calCache.set(key, m); }
  const f = m.f;
  let sat = m.sat.get(T);
  if (!sat) { const seedKey = `${v.fluidSource}|${T}`, seed = calCache.get(seedKey); sat = saturationPoint(f, T, { seed: seed || null }); if (sat.X) calCache.set(seedKey, sat.X); m.sat.set(T, sat); }
  const s = pv(f, P, T), liq = s.phase === 'gas' ? null : s.oil;
  let rs = null;
  if (liq) { const fl = stdFlash(withZ(f, s.x)); rs = fl.vOil > 0 ? fl.vGas / fl.vOil : null; }
  const out = { psat: sat.P ?? NaN, rho: (liq || s.gas).rho, mu: (liq || s.gas).mu * 1e3, rs: rs ?? NaN };
  const lazy = (k, fn) => Object.defineProperty(out, k, { enumerable: true, configurable: true, get() { let q; try { q = +fn(); } catch { q = NaN; } if (!Number.isFinite(q)) q = NaN; Object.defineProperty(out, k, { value: q, enumerable: true, configurable: true }); return q; } });
  const inhId = INHIBITORS[v.inhDesign] && v.inhDesign !== 'none' ? v.inhDesign : 'MEG', inhWt = clamp(num(v.calInhWt, 0), 0, 90), S = clamp(num(v.calSal, (lastModel?.aq?.S ?? 35) / 10), 0, 26) * 10, aq = { S, inhId: inhWt > 0 ? inhId : 'none', inhWt };
  const models = { act: ['NRTL', 'UNIQUAC', 'Wilson', 'ideal'].includes(v.actModel) ? v.actModel : 'NRTL', elec: ['pitzer', 'davies', 'edh', 'dh'].includes(v.elecModel) ? v.elecModel : 'pitzer' }, awT = (TK) => waterActivity(aq, TK, models, tx).aw;
  let gm = null, lm = null;
  const gas = () => (gm ||= (() => { const st = stdFlash(f), fg = withZ(f, st.beta > 0 ? st.y : f.z); return derivedProps(fg, fg.z, P, T, 'vapour'); })()), liqD = () => (lm ||= liq ? derivedProps(f, s.x, P, T, 'liquid') : gas());
  lazy('bo', () => { if (!liq) return NaN; const fl = stdFlash(withZ(f, s.x)); return fl.vOil > 0 ? liq.vm / fl.vOil : NaN; });
  lazy('rhoG', () => gas().rho); lazy('zG', () => gas().Z); lazy('muG', () => gas().mu * 1e3 * tx.muG); lazy('sound', () => gas().sound);
  lazy('cp', () => liqD().cp / 1000); lazy('co', () => liqD().kT * 1e4); lazy('sigma', () => (s.phase === 'two' ? s.sigma * 1e3 : NaN));
  lazy('rhoW', () => waterProps(P, T, { S, inhWt, inh: INHIBITORS[inhId] }).rho * tx.rhoW); lazy('wc', () => waterContent(P, T) * awT(T + KEL) * tx.wc * 1e6); lazy('aw', () => awT(T + KEL));
  lazy('hydT', () => { const r = hydrateTofP(f, P, awT, { table: hydrateTable(v.hydModel === 'vdwp' ? 'vdwp' : 'kihara', tx.hyd), sH: v.shFormer === true, molality: waterActivity(aq, 277, models).molality, guess: 10 }); return r ? r.T : NaN; });
  lazy('mmp', () => miscibilityPressure(f, (INJECTION_GASES[v.swellGas] || INJECTION_GASES.lean).comp, T).mmp);
  return out;
}

// ---- verification -----------------------------------------------------------------------------------------------------------------
function purePsat(id, TK) { // vapour pressure of a pure component from equal fugacity of the liquid and vapour roots
  const f = makeFluid({ comp: { [id]: 100 } }), c = f.comps[0], g = (lp) => { const P = Math.exp(lp); return eosPhase(f, [1], P, TK, 'liquid').lnphi[0] - eosPhase(f, [1], P, TK, 'vapour').lnphi[0]; };
  const guess = c.Pc * Math.exp(5.373 * (1 + c.w) * (1 - c.Tc / TK));
  let lo = Math.log(guess * 0.85), hi = Math.log(Math.min(guess * 1.15, c.Pc * 0.999)); // liquid root has the higher fugacity below the vapour pressure
  for (let k = 0; k < 60; k++) { const m = 0.5 * (lo + hi); if (g(m) > 0) lo = m; else hi = m; }
  return Math.exp(0.5 * (lo + hi));
}
function verifyChecks() {
  const out = [], chk = (name, expected, got, tol, note) => out.push({ name, expected, got: Number.isFinite(got) ? +(+got).toPrecision(7) : got, tol, pass: Number.isFinite(got) && Math.abs(got - expected) <= tol, note });
  const f = makeFluid(DEFAULT_FLUID), c1 = makeFluid({ comp: { C1: 100 } });
  chk('Pure-component EOS: methane vapour pressure at 150 K', 10.4, purePsat('C1', 150), 0.25, 'Peng–Robinson against the reference equation of state (1.040 MPa), bar');
  chk('Pure-component EOS: propane vapour pressure at 300 K', 9.98, purePsat('C3', 300), 0.2, 'Reference value 0.998 MPa, bar');
  chk('Pure-component EOS: carbon dioxide vapour pressure at 273.15 K', 34.85, purePsat('CO2', 273.15), 0.6, 'Reference value 3.485 MPa, bar');
  chk('Pure-component EOS: n-pentane at its normal boiling point (309.2 K)', 1.013, purePsat('nC5', 309.22), 0.03, 'Vapour pressure must equal one atmosphere, bar');
  chk('Critical-point verification: Peng–Robinson critical compressibility', 0.3074, eosPhase(c1, [1], c1.comps[0].Pc, c1.comps[0].Tc).Z, 2e-3, 'Analytic Zc of the Peng–Robinson equation at (Tc, Pc)');
  chk('Limiting single-phase case: ideal-gas limit of Z', 1, eosPhase(c1, [1], 0.01, 300).Z, 1e-4, 'Methane at 0.01 bara, 300 K');
  chk('Speed of sound of methane at 300 K, 1 bara', 450, derivedProps(c1, [1], 1, 26.85, 'vapour').sound, 4, 'Reference value about 450 m/s (ideal-gas sqrt(γRT/M) = 450.1 m/s)');
  chk('Compressibility factor of methane at 300 K, 100 bara (Lee–Kesler)', 0.849, leeKeslerZ(300 / 190.56, 100 / 45.99, 0.0115), 0.012, 'Reference equation of state: Z ≈ 0.849');
  chk('Rachford–Rice equation: analytic binary case', 0.5, rachfordRice([0.5, 0.5], [2, 0.5]), 1e-10, 'z = (0.5, 0.5), K = (2, 0.5) gives a vapour fraction of exactly one half');
  const fl = flashPT(f, 100, 60);
  { let mb = 0, fe = 0; f.z.forEach((z, i) => { mb = Math.max(mb, Math.abs(fl.beta * fl.y[i] + (1 - fl.beta) * fl.x[i] - z)); if (fl.x[i] > 1e-200 && fl.y[i] > 1e-200) fe = Math.max(fe, Math.abs(Math.log(fl.x[i]) + fl.liq.lnphi[i] - Math.log(fl.y[i]) - fl.vap.lnphi[i])); });
    chk('Component material-balance closure of the multicomponent PT flash', 0, mb, 1e-12, 'Largest |βy + (1 − β)x − z| at 100 bara, 60 °C (machine precision)');
    chk('Fugacity equality at equilibrium', 0, fe, 1e-6, 'Largest |ln f(liquid) − ln f(vapour)| over the components'); }
  { // Gibbs–Duhem: Σ x d ln φ / dP = (Z − 1)/P at constant T and composition
    const x = fl.y, P = 100, T = 333.15, d = 1e-3, a = eosPhase(f, x, P + d, T, 'vapour'), b = eosPhase(f, x, P - d, T, 'vapour'), m = eosPhase(f, x, P, T, 'vapour'); let sm = 0; x.forEach((xi, i) => (sm += (xi * (a.lnphi[i] - b.lnphi[i])) / (2 * d)));
    chk('Gibbs–Duhem consistency of the fugacity coefficients', 1, sm / ((m.Z - 1) / P), 1e-5, 'Σ x ∂ln φ/∂P divided by (Z − 1)/P for the equilibrium vapour'); }
  { // analytic versus numerical temperature derivative: residual enthalpy from the Gibbs–Helmholtz relation
    const x = fl.x, P = 100, T = 333.15, d = 0.01, g = (t) => { const p = eosPhase(f, x, P, t, 'liquid'); let s = 0; x.forEach((xi, i) => (s += xi * p.lnphi[i])); return s; };
    chk('Analytic-vs-numerical derivative: residual enthalpy', 1, phaseHS(f, x, P, T, 'liquid').hRes / (-R * T * T * (g(T + d) - g(T - d)) / (2 * d)), 1e-5, 'Analytic da/dT against −RT² ∂(Σ x ln φ)/∂T by central differences');
    const s1 = phaseHS(f, x, P + 0.01, T, 'liquid'), s2 = phaseHS(f, x, P - 0.01, T, 'liquid'), v = (t) => (eosPhase(f, x, P, t, 'liquid').Z * R * t) / (P * 1e5);
    chk('Maxwell relation (∂S/∂P)T = −(∂V/∂T)P', 1, ((s1.s - s2.s) / (0.02 * 1e5)) / (-(v(T + d) - v(T - d)) / (2 * d)), 1e-4, 'Liquid phase of the reference oil at 100 bara, 60 °C'); }
  const st = stateHSV(f, 100, 60), ph = flashPH(f, 30, st.H), back = ph ? flashPH(f, 100, stateHSV(f, 30, ph.T).H) : null;
  chk('PH-flash round trip', 60, back ? back.T : NaN, 1e-4, `Isenthalpic expansion 100 → 30 bara (${ph ? ph.T.toFixed(2) : '—'} °C) and recompression at the same enthalpy must return 60 °C`);
  const pS = flashPS(f, 30, st.S), backS = pS ? flashPS(f, 100, stateHSV(f, 30, pS.T).S) : null;
  chk('PS-flash round trip', 60, backS ? backS.T : NaN, 1e-4, 'Isentropic expansion 100 → 30 bara and back');
  chk('Energy-balance closure of the PH flash', st.H, ph ? ph.H : NaN, 1e-3, 'Total molar enthalpy before and after the isenthalpic flash, J/mol');
  const tv = flashTV(f, 60, st.V); chk('TV-flash test', 100, tv ? tv.P : NaN, 1e-5, 'Pressure recovered from temperature and molar volume, bara');
  { const sp = saturationPoint(f, 90), ref = saturationP(f, 90).P; let sk = 0; f.z.forEach((z, i) => (sk += z * Math.exp(sp.X[i])));
    chk('Bubble-point solver: Σ z K = 1 at the bubble point', 1, sk, 1e-8, 'Incipient-vapour mole fractions sum to one');
    chk('Cross-implementation benchmark: Newton bubble point against flash bisection', ref, sp.P, 0.05, 'Two independent algorithms for the saturation pressure of the reference oil at 90 °C, bara');
    chk('Bubble-point limit of the flash', 0, flashPT(f, sp.P * 0.9995, 90).beta, 2e-3, 'Vapour fraction just below the bubble point tends to zero');
    chk('Gibbs-energy stability test above the bubble point', 1, stability(f, sp.P * 1.05, 363.15).stable ? 1 : 0, 0, 'Tangent-plane test reports a stable single phase 5 % above the bubble point');
    chk('Gibbs-energy stability test below the bubble point', 0, stability(f, sp.P * 0.9, 363.15).stable ? 1 : 0, 0, 'Tangent-plane distance is negative 10 % below the bubble point'); }
  { const b = makeFluid({ comp: { C3: 50, nC4: 50 } }), sp = saturationPoint(b, 26.85), dp = (() => { const lg = makeFluid({ comp: FLUID_LIBRARY.leanGas.spec.comp, c7MW: 110, c7SG: 0.75 }), t = dewT(lg, 30), s = t !== null ? flashPT(lg, 30, t - 0.05) : null; return s ? s.beta : NaN; })();
    chk("Binary-mixture verification: Raoult's-law limit for propane–n-butane at 300 K", 0.5 * 9.98 + 0.5 * 2.58, sp.P, 0.32, 'Bubble pressure of a near-ideal equimolar mixture against the mole-fraction average of the pure vapour pressures (within 5 %), bar');
    chk('Dew-point solver: vapour fraction at the dew temperature', 1, dp, 5e-3, 'Lean gas at 30 bara just inside its dew point'); }
  { const env = traceEnvelope(f, { n: 30, qualities: [] }), s = env.critical ? flashPT(f, env.critical.P * 0.97, env.critical.T - 4) : null;
    chk('Phase-envelope tracing: K-values approach one at the traced critical point', 0, s && s.phase === 'two' ? Math.max(...s.K.slice(0, 6).map((k) => Math.abs(Math.log(k)))) : s ? 0 : NaN, 0.35, 'Largest |ln K| of the light components 3 % in pressure and 4 K below the critical point'); }
  const hyd = (comp, P) => { const r = hydrateTofP(makeFluid({ comp }), P); return r ? r.T + KEL : NaN; };
  { const pick = (id, lo, hi) => HYDRATE_DATA.pure.find((r) => r.y[id] === 100 && r.P >= lo && r.P <= hi), src = (r) => `${HYDRATE_REFS[r.ref].s}: ${r.T} K at ${r.P} bar`;
    for (const [nm, id, lo, hi, tol] of [['methane', 'C1', 30, 100, 1.5], ['ethane', 'C2', 5, 30, 2.3], ['propane', 'C3', 2, 6, 1], ['carbon dioxide', 'CO2', 15, 45, 1.5], ['hydrogen sulphide', 'H2S', 2, 20, 1]]) { const r = pick(id, lo, hi); if (r) chk(`Hydrate dissociation of ${nm} (default Kihara set) against a measured point`, r.T, hyd({ [id]: 100 }, r.P), tol, `${src(r)}; the default set errs warm (ethane by about 1.6 K)`); }
    const m = makeFluid({ comp: { C1: 100 } }), h7 = hydrateTofP(m, 700, 1).T;
    chk('Published hydrate curve: high-pressure allowance at 700 bara', 1.8, publishedHydrateT(h7, h7, 0, 700) - h7, 1e-9, '0.004 K per bar above 250 bara');
    chk('Published hydrate curve: share of the inhibitor depression withheld at 20, 40 and 60 wt %', 0 + 0.06 + 0.12, inhibitorAllowance(20) + inhibitorAllowance(40) + inhibitorAllowance(60), 1e-12, 'Zero to 30 wt %, 12 % from 50 wt %');
    const ng = makeFluid({ comp: { C1: 88, C2: 6, C3: 3, nC4: 1, N2: 1, CO2: 1 } }); let worst = -Infinity;
    for (const [inhId, w, S] of [['MEG', 30, 0], ['MEG', 50, 0], ['MEG', 60, 0], ['MeOH', 20, 0], ['MeOH', 40, 0], ['MeOH', 50, 0], ['none', 0, 35], ['none', 0, 150], ['none', 0, 250]]) for (const P of [30, 100]) {
      const aq = { S, inhId, inhWt: w }, t0 = hydrateTofP(ng, P, 1).T, tm = hydrateTofP(ng, P, (TK) => waterActivity(aq, TK).aw, { molality: waterActivity(aq, 277).molality, guess: t0 - 8, Tmin: -70 }).T;
      worst = Math.max(worst, hydrateDepression({ S, inhWt: w, inh: INHIBITORS[inhId] }) - (t0 - publishedHydrateT(tm, t0, w, 0))); }
    chk('Kernel screening depression never exceeds the published model depression', 0, Math.max(0, worst), 1e-9, `Natural gas at 30 and 100 bara, MEG to 60 wt %, methanol to 50 wt %, NaCl to 25 wt %: largest (kernel − model) ${worst.toFixed(2)} K`);
    const tb = buildTable(DEFAULT_FLUID, { nP: 14, nT: 9 }), scr = hydrateScreening(tb), fr = makeFluid(DEFAULT_FLUID), res = [25, 100, 300].map((P) => { const q = hydrateTofP(fr, P, 1).T; return scr(P) - publishedHydrateT(q, q, 0, P); });
    chk('Kernel screening hydrate curve of the reference fluid lies on the warm side of the model, within 3 K', 1, res.every((q) => q > -0.3 && q < 3) ? 1 : 0, 0, `Screening − model at 25, 100 and 300 bara: ${res.map((q) => q.toFixed(2)).join(', ')} K`);
    let mono = 0; for (const g of [0.554, 0.65, 0.8, 1]) for (let P = 2; P < 700; P *= 1.3) if (hydrateT0(P * 1.3, g) < hydrateT0(P, g) || hydrateT0(P, g + 0.05) < hydrateT0(P, g) - 1e-9) mono++;
    chk('Kernel gas-gravity correlation rises with pressure and never falls with gravity', 0, mono, 0, '1–700 bara, gravity 0.554–1.0');
    const lg = makeFluid({ comp: { C1: 92.7, C2: 4.56, C3: 1.73, N2: 0.46, CO2: 0.55 } }), q70 = hydrateTofP(lg, 70, 1).T;
    chk('Kernel gas-gravity correlation against the model for a 0.60-gravity natural gas at 70 bara', q70, hydrateT0(70, lg.MW / MW_AIR), 2, 'Correlation fitted to the model with a 1 K warm shift; root-mean-square residual 1.5 K over 91 gases');
    chk('Activity model: limiting activity coefficient of methanol in water at 273.35 K', 1.23, activityBinary('NRTL', 'MeOH', 1e-9, 273.35).g1, 0.12, 'Vrbka et al. (2005), a calibration point');
    chk('Activity model: ice point of 50 wt % MEG from parameters fitted to vapour–liquid data only', -33.6, freezingPoint({ inhId: 'MEG', inhWt: 50 }), 0.6, 'Extrapolation from 333–470 K to the ice line; handbook values lie between −33 and −37 °C');
    const b = BENCH.state.find((r) => r.mix === 'methane-propane 70/30' && r.T === 300 && r.P === 200), fb = makeFluid({ comp: b.y }), ga = altPhase(fb, fb.z, b.P, b.T - KEL, 'vapour', { visc: false }).gerg;
    chk('GERG-2008 mixture model against independent reference software', b.rho, ga ? ga.rho : NaN, 2e-3 * b.rho, `${BENCH_META.software} ${BENCH_META.version}: methane–propane 70/30 at 300 K, 200 bar, kg/m³`); }
  chk('Structure selection: 1 % propane turns methane hydrate into structure II', 2, hydrateTofP(makeFluid({ comp: { C1: 99, C3: 1 } }), 43.6).structure === 'sII' ? 2 : 1, 0, 'Lowest dissociation pressure decides the structure');
  chk('Inhibitor response: 50 wt % MEG on methane hydrate at 99 bara', 22.6, (() => { const m = makeFluid({ comp: { C1: 100 } }); return hydrateTofP(m, 98.9).T - hydrateTofP(m, 98.9, (TK) => waterActivity({ inhId: 'MEG', inhWt: 50 }, TK).aw).T; })(), 2.6, 'Measured depression of about 22–23 K (Robinson & Ng 1986); 16–17 K for methane–propane (Song & Kobayashi 1989): the model (20.6 K) lies between the two laboratories');
  chk('Debye–Hückel osmotic constant of water at 25 °C', 0.3915, debyeHuckel(298.15).Aphi, 0.002, 'From the density and dielectric constant of water');
  chk('Pitzer osmotic coefficient of 1 mol/kg NaCl at 25 °C', 0.936, osmoticNaCl('pitzer', 1), 0.003, 'Hamer & Wu (1972) Table 16');
  chk('Pitzer osmotic coefficient of 6 mol/kg NaCl at 25 °C', 1.27, osmoticNaCl('pitzer', 6), 0.01, 'Hamer & Wu (1972) Table 16');
  chk('Pitzer mean activity coefficient of 1 mol/kg NaCl at 25 °C', 0.657, gammaNaCl('pitzer', 1), 0.005, 'Hamer & Wu (1972) Table 16');
  chk('Davies mean activity coefficient of 0.1 mol/kg NaCl at 25 °C', 0.778, gammaNaCl('davies', 0.1), 0.01, 'Measured value 0.778');
  { const m = 2, d = 1e-4, lhs = (osmoticNaCl('pitzer', m + d) * (m + d) - osmoticNaCl('pitzer', m - d) * (m - d)) / (2 * d) - 1, rhs = (m * (Math.log(gammaNaCl('pitzer', m + d)) - Math.log(gammaNaCl('pitzer', m - d)))) / (2 * d);
    chk('Gibbs–Duhem consistency of the Pitzer model', 0, lhs - rhs, 1e-6, 'd[m(φ − 1)]/dm = m d ln γ±/dm at 2 mol/kg'); }
  { let worst = 0; for (const md of ['NRTL', 'UNIQUAC', 'Wilson']) for (const inh of ['MeOH', 'MEG']) { const x = 0.3, d = 1e-5, a = activityBinary(md, inh, x + d, 290), b = activityBinary(md, inh, x - d, 290); worst = Math.max(worst, Math.abs((x * Math.log(a.g1 / b.g1) + (1 - x) * Math.log(a.g2 / b.g2)) / (2 * d))); }
    chk('Gibbs–Duhem consistency of NRTL, UNIQUAC and Wilson', 0, worst, 1e-6, 'x₁ d ln γ₁ + x₂ d ln γ₂ = 0 for methanol–water and MEG–water at 290 K'); }
  chk("Henry's law: methane solubility in water at 25 °C, 1 atm", 2.55e-5, gasSolubility({ C1: ATM }, ATM, 25).x.C1, 1.5e-6, 'Mole fraction, IUPAC solubility data');
  chk("Henry's law: carbon dioxide solubility in water at 25 °C, 1 atm", 6.15e-4, gasSolubility({ CO2: ATM }, ATM, 25).x.CO2, 3e-5, 'Mole fraction');
  chk('Twu critical temperature of n-heptane', 540.2, critTwu(371.6, 0.6882).Tc, 1, 'From Tb = 371.6 K and SG = 0.6882; measured 540.2 K');
  chk('Twu critical pressure of n-heptane', 27.4, critTwu(371.6, 0.6882).Pc, 0.3, 'Measured 27.4 bar');
  chk('Riazi–Daubert critical temperature of n-decane', 617.7, critRiaziDaubert(447.3, 0.7342).Tc, 3, 'From Tb = 447.3 K and SG = 0.7342; measured 617.7 K');
  chk('Lee–Gonzalez–Eakin viscosity of methane at 300 K, 1 bara', 11.2, viscosityLGE(16.043, 0.6443, 300) * 1e6, 0.6, 'Reference value 11.2 µPa·s');
  chk('Unit consistency: molar volume of an ideal gas at standard conditions', 23.6445, VM_STD * 1e3, 1e-3, '15 °C and 1.01325 bara, litres per mole');
  chk('Unit consistency: API gravity of water', 10, 141.5 / 1 - 131.5, 1e-12, '°API = 141.5/SG − 131.5');
  { const t = buildTable(DEFAULT_FLUID, { nP: 9, nT: 7, Pmin: 20, Pmax: 200, Tmin: 20, Tmax: 100 }), q = lookup(t, 77, 43), d = props(f, 77, 43, { thermal: false });
    chk('Property-table interpolation error (oil density at an off-grid point)', d.oil.rho, q.rhoO, 0.01 * d.oil.rho, 'Bilinear look-up in ln P and T against a direct flash at 77 bara, 43 °C, kg/m³'); }
  { const a = calibrationModel({ tcMult: 1.02, pcMult: 0.97, calT: 90, calP: 150 }), b = calibrationModel({ tcMult: 1.02, pcMult: 0.97, calT: 90, calP: 150 });
    chk('Regression reproducibility', a.psat, b.psat, 0, 'Identical parameters give bit-identical predictions (deterministic model)'); }
  // ---- sourced constants, Helmholtz models, association, corresponding states, solids --------------------------------------------
  { const ref = (fl, T, P) => REF_ISO.find((r) => r.f === fl && r.T === T && r.P === P), a = ref('C1', 300, 101), g = gergModel(['C1']), pc = pcsaftModel([{ ...COMPONENTS.C1 }]), sat = REF_SAT.find((r) => r.f === 'C1' && r.T === 155);
    chk('Kernel constants: critical temperature of methane equals the source value', 190.564, COMPONENTS.C1.Tc, 1e-9, 'CoolProp fluid file (Setzmann & Wagner reference equation), K');
    chk('GERG-2008: methane density at 300 K, 101 bar', a.rho, hRho(g, 300, 101, [1]) * g.M[0] * 1e-3, 0.002 * a.rho, 'NIST WebBook value, kg/m³ (within 0.2 %)');
    chk('GERG-2008: vapour pressure of methane at 155 K', sat.P, hPsat(g, 155).P, 0.003 * sat.P, 'NIST WebBook saturation line, bar');
    chk('PC-SAFT: vapour pressure of methane at 155 K', sat.P, hPsat(pc, 155).P, 0.02 * sat.P, 'NIST WebBook saturation line, bar (Gross & Sadowski parameters, within 2 %)');
    chk('PC-SAFT: ideal-gas limit of Z', 1, hZ(pc, 300, 0.4, [1]), 1e-4, 'Methane at 0.4 mol/m³ (about 0.01 bar), 300 K');
    const b = ref('C1', 300, 1), q = hProps(g, 300, hRho(g, 300, 1, [1]), [1], cpIg(c1, [1], 300));
    chk('Helmholtz formulation: speed of sound of methane at 300 K, 1 bar (GERG-2008)', b.w, q.w, 0.004 * b.w, 'From the residual Helmholtz energy and its derivatives plus the kernel ideal-gas heat capacity; NIST value, m/s');
    const mu = ref('C1', 300, 101); chk('Pedersen viscosity: methane at 300 K, 101 bar', mu.mu, viscosityPedersen([COMPONENTS.C1], [1], 101, 300).mu * 1e6, 0.03 * mu.mu, 'Reference fluid of the corresponding-states model; NIST value, µPa·s');
    chk('Methane reference viscosity (Hanley et al.) at 300 K, 1 bar', b.mu, methaneViscosity(300, b.rho) * 1e6, 0.02 * b.mu, 'Dilute-gas limit of the reference correlation; NIST value, µPa·s');
    const lk = leeKesler(300 / 190.564, 101 / 45.992, 0.01142), zN = (101e5 * 16.0428e-3) / (a.rho * R * 300); chk('Lee–Kesler (modified BWR): Z of methane at 300 K, 101 bar', zN, lk.Z, 0.012, 'NIST value from the reference density');
    const e = ref('C1', 300, 1), hN = (a.h - e.h) * 16.0428, hi = leeKesler(300 / 190.564, 1 / 45.992, 0.01142); chk('Lee–Kesler enthalpy departure: methane 1 → 101 bar at 300 K', hN, (lk.hDep - hi.hDep) * R * 190.564, 0.06 * Math.abs(hN), 'Isothermal enthalpy change from the NIST WebBook, J/mol (within 6 %)'); }
  { const m = gergModel(['C1', 'C2', 'CO2', 'N2']), x = [0.85, 0.08, 0.04, 0.03], T = 280, r = hRho(m, T, 80, x), P = hProps(m, T, r, x).P, d = 0.02, lp = (p) => hLnPhi(m, T, hRho(m, T, p, x), x), a = lp(P + d), b = lp(P - d); let sm = 0; x.forEach((xi, i) => (sm += (xi * (a[i] - b[i])) / (2 * d)));
    chk('Helmholtz formulation: Gibbs–Duhem consistency of the GERG-2008 fugacity coefficients', 1, sm / ((hZ(m, T, r, x) - 1) / P), 2e-4, 'Σ x ∂ln φ/∂P divided by (Z − 1)/P for a four-component natural gas at 280 K, 80 bar');
    const gT = (t) => { const rr = hRho(m, t, 80, x), l = hLnPhi(m, t, rr, x); let s = 0; x.forEach((xi, i) => (s += xi * l[i])); return s; };
    chk('Helmholtz formulation: residual enthalpy against the Gibbs–Helmholtz relation', 1, hProps(m, T, r, x).hRes / ((-R * T * T * (gT(T + 0.05) - gT(T - 0.05))) / 0.1), 2e-4, 'RT(Z − 1 − T ∂ar/∂T) against −RT² ∂(Σ x ln φ)/∂T'); }
  { const w = C1_WATER_VLE.find((r) => r.T > 298 && r.T < 299) || C1_WATER_VLE[0], y = cpaWater(c1, [1], w.P, w.T - KEL).yW; chk('Cubic-plus-association: water in methane', w.yw, y, 0.35 * w.yw, `Frost et al. (2014) at ${w.T} K, ${w.P} bar, mole fraction (single measurement, within 35 %)`);
    const m = cpaModel([{ id: 'W' }]), rw = hRho(m, 298.15, 1, [1], 'liquid') * 18.015e-3; chk('Cubic-plus-association: density of liquid water at 25 °C', 997, rw, 30, 'Published CPA water parameters reproduce the liquid density within about 3 %, kg/m³');
    const ps = hPsat(m, 373.15); chk('Cubic-plus-association: vapour pressure of water at 100 °C', 1.013, ps ? ps.P : NaN, 0.04, 'One atmosphere at the normal boiling point, bar'); }
  { const r1 = HYDRATE_DATA.pure.find((r) => r.y.C1 === 100), r2 = HYDRATE_DATA.sH.find((r) => r.ds === 986), m1 = makeFluid({ comp: { C1: 100 } });
    chk('Kihara cell-potential hydrate model (default): a measured methane point', r1.T, hydrateTofP(m1, r1.P, 1, { table: KIHARA, sH: false, guess: r1.T - KEL }).T + KEL, 1.2, `${HYDRATE_REFS[r1.ref].s}: ${r1.T} K at ${r1.P} bar`);
    chk('Kihara Langmuir constant: tabulated interpolation against the direct integral', 1, langmuirC('C1', 'sI', 1, 281.3, KIHARA) / kiharaC(...KIHARA.guests.C1, ...KIHARA.cav.sI[1], 281.3), 5e-4, 'ln C interpolated in 1/T on a 5 K grid');
    chk('Structure-H hydrate: methane + methylcyclohexane point not used in the regression', r2.T, hydrateTofP(m1, r2.P, 1, { table: KIHARA, sH: true, guess: r2.T - KEL }).T + KEL, 0.6, `${HYDRATE_REFS[r2.ref].s}: ${r2.T} K at ${r2.P} bar`);
    const eq = waterPhaseEquilibrium(m1, 100, 4, { S: 35 }, {}, { table: KIHARA, sH: false }), back = hydrateTofP(m1, 100, (TK) => waterActivity({ S: eq.saltFinal }, TK).aw, { table: KIHARA, sH: false, molality: waterActivity({ S: eq.saltFinal }, 277).molality, guess: 4 });
    chk('Water-side Gibbs minimum: the residual brine is at hydrate equilibrium', 4, back ? back.T : NaN, 0.35, 'Methane, 100 bara, 4 °C, 3.5 wt % brine: the hydrate temperature of the brine left after conversion, solved independently, equals the system temperature');
    chk('Water-side Gibbs minimum: fresh water converts completely', 1, waterPhaseEquilibrium(m1, 100, 4, { S: 0 }, {}, { table: KIHARA, sH: false }).conversion, 1e-9, 'No solute to stop the conversion and excess gas');
    chk('Gibbs energy decreases on hydrate formation', 1, eq.dG < 0 ? 1 : 0, 0, `ΔG = ${eq.dG.toFixed(0)} J per mol of water`); }
  { const g = gibbsSplit(f, 100, 60); chk('Gibbs-energy minimisation: the two-phase split lowers the Gibbs energy', 1, g.dg < 0 && g.tpd < 0 ? 1 : 0, 0, `Δg/RT = ${g.dg.toFixed(4)}, tangent-plane distance ${g.tpd.toFixed(3)} at 100 bara, 60 °C`);
    const h = NACL_25C.find((r) => r.m === 1), h6 = NACL_25C.find((r) => r.m === 6); chk('Pitzer (PHREEQC parameters): osmotic coefficient at 1 and 6 mol/kg', h.phi + h6.phi, osmoticNaCl('pitzer', 1) + osmoticNaCl('pitzer', 6), 0.006, 'Hamer & Wu (1972) Table 16');
    const fz = MEOH_FREEZING.find((r) => r.w === 20.6); chk('Freezing point of 20.6 wt % methanol (calibration point of the activity model)', fz.tf, freezingPoint({ inhId: 'MeOH', inhWt: 20.6 }), 1.5, 'Tabulated value −15.6 °C; the NRTL fit lies 1.3 K warm here (it errs towards less inhibition)');
    const fh = asphalteneFH({ deltaL: 16, vL: 2e-4, TK: 350, wAsp: 0.03, rhoL: 800, aspMW: 750, aspRho: 1100, aspDelta: 21, aspDeltaT: 0 }), Va = 0.75 / 1100, rr = Va / 2e-4;
    chk('Asphaltene Flory–Huggins equilibrium: residual of the solubility equation', 0, Math.log(fh.phiMax) + (1 - rr) * (1 - fh.phiMax) + fh.chi * (1 - fh.phiMax) ** 2, 1e-8, 'ln φ + (1 − Va/VL)(1 − φ) + χ(1 − φ)² at the returned soluble limit');
    chk('Asphaltene Flory–Huggins equilibrium: complete miscibility when the solubility parameters match', 1, asphalteneFH({ deltaL: 21, vL: 2e-4, TK: 350, wAsp: 0.03, rhoL: 800, aspDelta: 21, aspDeltaT: 0 }).phiMax, 1e-12, 'χ = 0 and Va > VL: no solid phase');
    const lean = miscibilityPressure(f, INJECTION_GASES.lean.comp, 90).mmp, rich = miscibilityPressure(f, INJECTION_GASES.rich.comp, 90).mmp; chk('Multiple-contact miscibility: enrichment lowers the miscibility pressure', 1, rich < lean ? 1 : 0, 0, `Lean gas ${lean.toFixed(0)} bara, enriched gas ${rich.toFixed(0)} bara at 90 °C`);
    const a = calibrationModel({ calT: 60, calP: 80 }), b = calibrationModel({ calT: 60, calP: 80, shiftMult: 1.3, parMult: 1.1 }); chk('Calibration parameters act on their targets', 1, Math.abs(b.rho - a.rho) > 0.5 && Math.abs(b.sigma / a.sigma - 1) > 0.05 ? 1 : 0, 0, `Volume-shift multiplier 1.3 moves the liquid density from ${a.rho.toFixed(1)} to ${b.rho.toFixed(1)} kg/m³; parachor multiplier 1.1 the interfacial tension from ${a.sigma.toFixed(2)} to ${b.sigma.toFixed(2)} mN/m`); }
  return out;
}

// ---- sourced reference data sets (rows in ../data/ref/pvt.js) with the engine's blind predictions ------------------------------
const pureCache = new Map();
const pureFluid = (id) => { let f = pureCache.get(id); if (!f) { f = makeFluid({ comp: { [id]: 100 } }); pureCache.set(id, f); } return f; };
const pureKind = (r) => (r.ph === 'l' ? 'liquid' : 'vapour');
const purePhase = (r) => derivedProps(pureFluid(r.f), [1], r.P, r.T - KEL, pureKind(r));
const mixCache = new Map();
const mixFluid = (y) => { const k = JSON.stringify(y); let f = mixCache.get(k); if (!f) { f = makeFluid({ comp: y }); if (mixCache.size > 300) mixCache.clear(); mixCache.set(k, f); } return f; };
const saltOf = (r) => (r.ws !== undefined ? r.ws * 10 : r.xs ? (1000 * r.xs * MW_NACL) / (r.xs * MW_NACL + (1 - r.xs) * MW_W) : 0); // g NaCl per kg brine
// Blind hydrate prediction of a data row. published = true returns the curve handed to the other suites (model value plus the inhibitor and high-pressure allowances).
const hydModelT = (table, sH = false, published = false) => (r) => {
  const S = saltOf(r), aq = { S, inhId: r.inh || 'none', inhWt: r.w || 0 }, plain = !(S > 0) && !(r.w > 0), f = mixFluid(r.y), o = { table, sH, guess: r.T - KEL, Tmax: 60, Tmin: -60 };
  const h = hydrateTofP(f, r.P, plain ? 1 : (TK) => waterActivity(aq, TK).aw, { ...o, molality: plain ? 0 : waterActivity(aq, 277).molality }); if (!h) return NaN;
  if (!published) return h.T + KEL;
  const h0 = plain ? h : hydrateTofP(f, r.P, 1, { ...o, guess: h.T + 5 }); return publishedHydrateT(h.T, h0 ? h0.T : h.T, r.w || 0, r.P) + KEL;
};
const gasLabel = (y) => Object.entries(y).map(([k, q]) => `${k} ${q}`).join(', ');
const withGas = (rows) => rows.map((r) => ({ ...r, gas: gasLabel(r.y), aq: r.w ? `${r.w} wt % ${r.inh}` : r.ws !== undefined ? `${r.ws} wt % NaCl` : r.xs ? `${(saltOf(r) / 10).toFixed(2)} wt % NaCl` : 'water', src: HYDRATE_REFS[r.ref]?.s || '' }));
const C_FL = { key: 'f', label: 'Fluid', type: 'text' }, C_T = { key: 'T', label: 'Temperature', unit: 'K' }, C_P = { key: 'P', label: 'Pressure', unit: 'bar' };
const isoZ = REF_ISO.filter((r) => r.ph !== 'l').map((r) => ({ ...r, Z: +((r.P * 1e5 * COMPONENTS[r.f].MW * 1e-3) / (r.rho * R * r.T)).toPrecision(5) }));
const isoLiq = [...REF_ISO.filter((r) => r.ph === 'l'), ...REF_SAT.filter((r) => r.T / COMPONENTS[r.f].Tc < 0.9).map((r) => ({ f: r.f, T: r.T, P: r.P, rho: r.rhoL, ph: 'l', sat: true }))];
const satRows = REF_SAT.map((r) => ({ ...r, dh: +(r.hV - r.hL).toPrecision(5), sigma: r.sig === null ? null : +(r.sig * 1e3).toPrecision(4) }));
const pureSat = (r) => { const f = pureFluid(r.f), P = purePsat(r.f, r.T); return { f, P, l: phaseProps(f, [1], P, r.T - KEL, 'liquid', { thermal: false }), v: phaseProps(f, [1], P, r.T - KEL, 'vapour', { thermal: false }) }; };
const cpSrc = (what) => ({ ...SOURCES.coolprop, citation: SOURCES.coolprop.citation + ` — ${what}` }), mixSrc = (what) => ({ ...SOURCES.coolpropMix, citation: SOURCES.coolpropMix.citation + ` — ${what}` });
const hydCite = (keys) => { const ids = [...new Set(keys.flatMap((k) => HYDRATE_DATA[k].map((r) => r.ref)))]; return { ...SOURCES.nistHyd, citation: `${ids.map((id) => HYDRATE_REFS[id].c).join(' | ')}. ${SOURCES.nistHyd.citation}` }; };
const hydCols = [{ key: 'src', label: 'Source', type: 'text' }, { key: 'gas', label: 'Dry gas (mol %)', type: 'text' }, { key: 'aq', label: 'Aqueous phase', type: 'text' }, C_P, { key: 'T', label: 'Measured dissociation temperature', unit: 'K' }];
const gergPure = new Map(), gergOf = (id) => { let m = gergPure.get(id); if (!m) { m = gergModel([id]); gergPure.set(id, m); } return m; };
const pcsPure = new Map(), pcsOf = (id) => { let m = pcsPure.get(id); if (!m) { m = pcsaftModel([{ ...COMPONENTS[id] }]); pcsPure.set(id, m); } return m; };
const hDens = (mOf) => (r) => { const m = mOf(r.f), q = hRho(m, r.T, r.P, [1], r.ph === 'l' ? 'liquid' : r.ph === 'v' ? 'vapour' : 'stable'); return q === null ? NaN : q * m.M[0] * 1e-3; };
const c1 = () => pureFluid('C1'), wcOfY = (y) => ((y / (1 - y)) * MW_W * 1e-3) / VM_STD * 1e6; // mg/Sm³
// Code-to-code benchmark helpers: the cubic model and the suite's own Helmholtz models at a state of a defined mixture.
const C_MIX = { key: 'mix', label: 'Mixture', type: 'text' }, C_COMP = { key: 'comp', label: 'Composition (mol %)', type: 'text' };
const bRows = (k, every = 1) => BENCH[k].filter((_, i) => i % every === 0).map((r) => ({ ...r, comp: gasLabel(r.y) }));
const bKind = (r) => (r.liq ? 'liquid' : 'vapour'), bCub = (r) => { const f = mixFluid(r.y); return derivedProps(f, f.z, r.P, r.T - KEL, bKind(r)); };
const altCache = new Map(), bAlt = (r, which) => { const k = `${which}|${r.mix}|${r.T}|${r.P}`; if (!altCache.has(k)) { const f = mixFluid(r.y), a = altPhase(f, f.z, r.P, r.T - KEL, bKind(r), { visc: false }); if (altCache.size > 2000) altCache.clear(); altCache.set(`gerg|${r.mix}|${r.T}|${r.P}`, a.gerg); altCache.set(`pcsaft|${r.mix}|${r.T}|${r.P}`, a.pcsaft); } return altCache.get(k) || {}; };
const bLiq = bRows('liquid').map((r) => ({ ...r, liq: true }));
const hOf = (y, T, P, which) => { const f = mixFluid(y), TK = T; if (which === 'pr') return phaseHS(f, f.z, P, TK, 'vapour').h / f.MW; const a = altPhase(f, f.z, P, T - KEL, 'vapour', { visc: false }).gerg; return a ? (hIg(f, f.z, TK) + a.hRes) / f.MW : NaN; }; // kJ/kg
const benchNote = `Reference values computed with ${BENCH_META.software} ${BENCH_META.version} (${BENCH_META.backend}) on ${BENCH_META.date}: a code-to-code comparison with independent reference software, reproducible from the stated compositions.`;
const hydAll = ['pure', 'highP', 'mix'];
const VALIDATION = [
  { id: 'ref-gas-z', title: 'Compressibility factor of gases and supercritical fluids, 250–450 K, 1–1000 bar', quantity: 'Compressibility factor Z', unit: '–', kind: 'reference-fluid', source: cpSrc('isotherms of methane, ethane, propane, n-butane, carbon dioxide, nitrogen and hydrogen sulphide'),
    columns: [C_FL, C_T, C_P, { key: 'Z', label: 'Z from the reference density' }], rows: isoZ, target: 'Z', model: (r) => purePhase(r).Z, tolerance: { mape: 3 }, note: 'Volume-translated Peng–Robinson with the kernel constants, no fitting. Largest deviations sit at the highest pressures and near the critical point; the GERG-2008 option reproduces the same data within 0.1 %.' },
  { id: 'ref-liquid-density', title: 'Liquid density of light hydrocarbons, CO₂ and H₂S (compressed and saturated)', quantity: 'Density', unit: 'kg/m³', kind: 'reference-fluid', source: cpSrc('compressed-liquid isotherms and saturated-liquid lines'),
    columns: [C_FL, C_T, C_P, { key: 'rho', label: 'Reference density', unit: 'kg/m³' }], rows: isoLiq, target: 'rho', model: (r) => phaseProps(pureFluid(r.f), [1], r.P, r.T - KEL, 'liquid', { thermal: false }).rho, tolerance: { mape: 5 }, note: 'A cubic equation of state with a constant volume shift: a few per cent is the expected accuracy, worst towards the critical temperature.' },
  { id: 'ref-vapour-pressure', title: 'Vapour pressure (bubble = dew pressure of the pure components)', quantity: 'Saturation pressure', unit: 'bar', kind: 'reference-fluid', source: cpSrc('saturation lines'),
    columns: [C_FL, C_T, { key: 'P', label: 'Reference vapour pressure', unit: 'bar' }], rows: REF_SAT, target: 'P', model: (r) => purePsat(r.f, r.T), tolerance: { mape: 2.5 }, note: 'Equal-fugacity solution of the Peng–Robinson equation with the acentric factors of the kernel; the phase boundary of each pure component from 0.55 to 0.98 of its critical temperature.' },
  { id: 'ref-heat-capacity', title: 'Isobaric heat capacity of gases, supercritical fluids and liquids', quantity: 'Heat capacity cp', unit: 'J/g/K', kind: 'reference-fluid', source: cpSrc('isotherms (calorimetric properties of the reference equations of state)'),
    columns: [C_FL, C_T, C_P, { key: 'cp', label: 'Reference cp', unit: 'J/g/K' }], rows: REF_ISO, target: 'cp', model: (r) => purePhase(r).cp / 1000, tolerance: { mape: 8 }, note: 'Ideal-gas polynomial plus the residual heat capacity of the cubic equation. Liquid heat capacities of a cubic equation are its weakest calorimetric property.' },
  { id: 'ref-enthalpy-vaporisation', title: 'Enthalpy of vaporisation along the saturation line', quantity: 'Enthalpy difference vapour − liquid', unit: 'kJ/kg', kind: 'reference-fluid', source: cpSrc('saturation lines (enthalpy of the coexisting phases)'),
    columns: [C_FL, C_T, { key: 'dh', label: 'Reference enthalpy of vaporisation', unit: 'kJ/kg' }], rows: satRows, target: 'dh', model: (r) => { const s = pureSat(r); return (phaseHS(s.f, [1], s.P, r.T, 'vapour').h - phaseHS(s.f, [1], s.P, r.T, 'liquid').h) / COMPONENTS[r.f].MW; }, tolerance: { mape: 7 }, note: 'Residual enthalpies of both phases with the analytic temperature derivative of the attraction term, at the model\'s own vapour pressure.' },
  { id: 'ref-speed-of-sound', title: 'Speed of sound in gases, supercritical fluids and liquids', quantity: 'Speed of sound', unit: 'm/s', kind: 'reference-fluid', source: cpSrc('isotherms'),
    columns: [C_FL, C_T, C_P, { key: 'w', label: 'Reference speed of sound', unit: 'm/s' }], rows: REF_ISO, target: 'w', model: (r) => purePhase(r).sound, tolerance: { mape: 12 }, note: 'Gas-phase values agree within a few per cent; in compressed liquids the cubic equation is known to be poor (the translated volume does not correct (∂P/∂v)), which dominates the average. Use the GERG-2008 or PC-SAFT option where the liquid speed of sound matters.' },
  { id: 'ref-joule-thomson', title: 'Joule–Thomson coefficient of gases and supercritical fluids', quantity: 'Joule–Thomson coefficient', unit: 'K/bar', kind: 'reference-fluid', source: cpSrc('isotherms'),
    columns: [C_FL, C_T, C_P, { key: 'jt', label: 'Reference coefficient', unit: 'K/bar' }], rows: REF_ISO.filter((r) => r.ph !== 'l'), target: 'jt', model: (r) => purePhase(r).jtBar, tolerance: { rmse: 0.08 }, note: 'Includes the inversion region at high pressure, where the coefficient changes sign; judged by the absolute error.' },
  { id: 'ref-compressibility', title: 'Isothermal compressibility from the reference densities of methane, ethane and carbon dioxide', quantity: 'Density change over a 100-bar (50-bar) step', unit: 'kg/m³', kind: 'reference-fluid', source: cpSrc('isotherms; differences between neighbouring pressures'),
    columns: [C_FL, C_T, { key: 'P', label: 'Lower pressure', unit: 'bar' }, { key: 'P2', label: 'Upper pressure', unit: 'bar' }, { key: 'drho', label: 'Reference density increase', unit: 'kg/m³' }],
    rows: (() => { const o = []; for (let k = 1; k < REF_ISO.length; k++) { const a = REF_ISO[k - 1], b = REF_ISO[k]; if (a.f === b.f && a.T === b.T && a.ph === b.ph && ['C1', 'CO2', 'C2'].includes(a.f)) o.push({ f: a.f, T: a.T, P: a.P, P2: b.P, ph: a.ph, drho: +(b.rho - a.rho).toPrecision(5) }); } return o; })(), target: 'drho',
    model: (r) => { const f = pureFluid(r.f), k = pureKind(r); return phaseProps(f, [1], r.P2, r.T - KEL, k, { thermal: false }).rho - phaseProps(f, [1], r.P, r.T - KEL, k, { thermal: false }).rho; }, tolerance: { mape: 15 }, note: 'The compressibility is the slope of density with pressure: the density increase over each pressure step is compared.' },
  { id: 'ref-viscosity-lbc', title: 'Viscosity by Lohrenz–Bray–Clark', quantity: 'Viscosity', unit: 'µPa·s', kind: 'reference-fluid', source: cpSrc('isotherms (reference viscosity correlations)'),
    columns: [C_FL, C_T, C_P, { key: 'mu', label: 'Reference viscosity', unit: 'µPa·s' }], rows: REF_ISO.filter((r) => r.mu !== null), target: 'mu', model: (r) => purePhase(r).mu * 1e6, tolerance: { mape: 12 }, note: 'Untuned; dense-fluid values depend on the fourth power of a polynomial in reduced density and so on the density error of the cubic equation.' },
  { id: 'ref-viscosity-pedersen', title: 'Viscosity by the Pedersen corresponding-states model', quantity: 'Viscosity', unit: 'µPa·s', kind: 'reference-fluid', source: cpSrc('isotherms (reference viscosity correlations)'),
    columns: [C_FL, C_T, C_P, { key: 'mu', label: 'Reference viscosity', unit: 'µPa·s' }], rows: REF_ISO.filter((r) => r.mu !== null), target: 'mu', model: (r) => { const q = viscosityPedersen([COMPONENTS[r.f]], [1], r.P, r.T); return q ? q.mu * 1e6 : NaN; }, tolerance: { mape: 8 }, note: 'Methane itself is reproduced within about 1 % (reference fluid); carbon dioxide near its critical point is the worst case.' },
  { id: 'ref-surface-tension', title: 'Surface tension of pure hydrocarbons, CO₂ and H₂S (parachor method)', quantity: 'Surface tension', unit: 'mN/m', kind: 'reference-fluid', source: cpSrc('saturation lines (surface tension)'),
    columns: [C_FL, C_T, { key: 'sigma', label: 'Reference surface tension', unit: 'mN/m' }], rows: satRows.filter((r) => r.sigma > 0.5 && r.f !== 'N2'), target: 'sigma', model: (r) => { const s = pureSat(r); return interfacialTension(s.f, [1], [1], s.l.vm, s.v.vm) * 1e3; }, tolerance: { mape: 22 }, note: 'Macleod–Sugden with the kernel parachors and the saturated densities of the equation of state: the fourth power amplifies any density error, so the deviation grows from a few per cent at low reduced temperature to −40 % and more above 0.9 Tc. Nitrogen is left out: its tabulated parachor (41) is the value effective in hydrocarbon mixtures and under-predicts the surface tension of pure nitrogen by about 80 %. No openly licensed gas–oil interfacial-tension measurements of mixtures were found: for mixtures the parachor method is checked only through these pure-component limits.' },
  { id: 'ref-hpht-density', title: 'High-pressure, high-temperature density (400–450 K, 300–1000 bar)', quantity: 'Density', unit: 'kg/m³', kind: 'reference-fluid', source: cpSrc('isotherms of methane (450 K), propane and carbon dioxide (400 K)'),
    columns: [C_FL, C_T, C_P, { key: 'rho', label: 'Reference density', unit: 'kg/m³' }], rows: REF_ISO.filter((r) => r.T >= 400 && r.P >= 300), target: 'rho', model: (r) => purePhase(r).rho, tolerance: { mape: 5 }, note: 'High-pressure/high-temperature range of the property table.' },
  { id: 'ref-deepwater-density', title: 'High-pressure, low-temperature density (≤ 280 K, 50–600 bar: deep-water flowline conditions)', quantity: 'Density', unit: 'kg/m³', kind: 'reference-fluid', source: cpSrc('isotherms of methane (250 and 277.15 K), nitrogen (250 K), ethane and carbon dioxide (280 K)'),
    columns: [C_FL, C_T, C_P, { key: 'rho', label: 'Reference density', unit: 'kg/m³' }], rows: REF_ISO.filter((r) => r.T <= 280 && r.P >= 50), target: 'rho', model: (r) => purePhase(r).rho, tolerance: { mape: 5 }, note: 'Seabed-temperature range at flowline and shut-in pressures; the hydrate sets below cover the same range for the solid boundary.' },
  { id: 'gerg-density', title: 'GERG-2008 option: density of all seven fluids over the whole range', quantity: 'Density', unit: 'kg/m³', kind: 'reference-fluid', source: cpSrc('isotherms'),
    columns: [C_FL, C_T, C_P, { key: 'rho', label: 'Reference density', unit: 'kg/m³' }], rows: REF_ISO, target: 'rho', model: hDens(gergOf), tolerance: { mape: 0.2 }, note: 'Checks the bundled GERG-2008 coefficients and the Helmholtz machinery against the reference equations of state of the pure fluids (GERG-2008 uses shorter pure-fluid equations, hence the small residual for propane and n-butane).' },
  { id: 'pcsaft-density', title: 'PC-SAFT option: density of all seven fluids over the whole range', quantity: 'Density', unit: 'kg/m³', kind: 'reference-fluid', source: cpSrc('isotherms'),
    columns: [C_FL, C_T, C_P, { key: 'rho', label: 'Reference density', unit: 'kg/m³' }], rows: REF_ISO, target: 'rho', model: hDens(pcsOf), tolerance: { mape: 2 }, note: 'Gross & Sadowski parameters without polar terms: carbon dioxide (a quadrupolar molecule) near its critical point carries the largest errors.' },
  // ---- code-to-code benchmarks against CoolProp / GERG-2008 for defined mixtures
  { id: 'bench-dew-temperature', title: 'Phase envelope of natural-gas and light-oil mixtures: dew temperatures and cricondentherm', quantity: 'Dew temperature', unit: 'K', kind: 'benchmark', source: mixSrc('traced phase envelopes of five mixtures'),
    columns: [C_MIX, C_COMP, { key: 'kind', label: 'Point', type: 'text' }, C_P, { key: 'T', label: 'Reference dew temperature', unit: 'K' }], rows: bRows('envelope'), target: 'T', model: (r) => { const t = dewT(mixFluid(r.y), r.P, { Tmin: -150, Tmax: 400 }); return t === null ? NaN : t + KEL; }, tolerance: { rmse: 6, bias: 4 }, note: `${benchNote} Peng–Robinson with the kernel interaction parameters; the dew line of a gas with heavy ends is the most sensitive part of an envelope, and a cubic equation differs from GERG-2008 by a few kelvin there.` },
  { id: 'bench-bubble-pressure', title: 'Bubble-point pressure of a six-component light oil (260–380 K) and of methane–propane (180–260 K)', quantity: 'Bubble pressure', unit: 'bar', kind: 'benchmark', source: mixSrc('bubble branch of the traced envelopes'),
    columns: [C_MIX, C_COMP, C_T, { key: 'P', label: 'Reference bubble pressure', unit: 'bar' }], rows: bRows('bubble'), target: 'P', model: (r) => { const s = saturationPoint(mixFluid(r.y), r.T - KEL); return s.P === null ? NaN : s.P; }, tolerance: { mape: 6 }, note: `${benchNote} The saturation-pressure solver of the suite with Peng–Robinson; a defined multicomponent reservoir-type fluid without a plus fraction.` },
  { id: 'bench-saturated-liquid-density', title: 'Saturated-liquid density at the bubble point of the same mixtures', quantity: 'Density', unit: 'kg/m³', kind: 'benchmark', source: mixSrc('incipient-phase densities of the traced envelopes'),
    columns: [C_MIX, C_COMP, C_T, { key: 'rhoL', label: 'Reference saturated-liquid density', unit: 'kg/m³' }], rows: bRows('bubble'), target: 'rhoL', model: (r) => { const f = mixFluid(r.y), s = saturationPoint(f, r.T - KEL); return s.P === null ? NaN : phaseProps(f, f.z, s.P * 1.0001, r.T - KEL, 'liquid', { thermal: false }).rho; }, tolerance: { mape: 6 }, note: `${benchNote} Volume-translated Peng–Robinson at its own bubble pressure (oil density of a live fluid).` },
  { id: 'bench-flash-split', title: 'Constant-composition expansion of the light oil: vapour mole fraction below the bubble point (300 and 340 K)', quantity: 'Vapour mole fraction', unit: '–', kind: 'benchmark', source: mixSrc('two-phase PT flash'),
    columns: [C_MIX, C_COMP, C_T, C_P, { key: 'beta', label: 'Reference vapour fraction' }], rows: bRows('split'), target: 'beta', model: (r) => { const fl = flashPT(mixFluid(r.y), r.P, r.T - KEL); return fl.phase === 'two' ? fl.beta : fl.beta; }, tolerance: { rmse: 0.04 }, note: `${benchNote} The PVT-cell experiment on a defined fluid: phase split of the two-phase flash at 20–90 % of the bubble pressure.` },
  { id: 'bench-mixture-density', title: 'Density of natural-gas mixtures, 250–400 K, 10–600 bar (Peng–Robinson)', quantity: 'Density', unit: 'kg/m³', kind: 'benchmark', source: mixSrc('single-phase states of four gas mixtures'),
    columns: [C_MIX, C_COMP, C_T, C_P, { key: 'rho', label: 'Reference density', unit: 'kg/m³' }], rows: bRows('state'), target: 'rho', model: (r) => bCub(r).rho, tolerance: { mape: 3 }, note: `${benchNote} Lean, rich and CO₂-rich natural gas and methane–propane 70/30, gas and dense phase.` },
  { id: 'bench-mixture-density-gerg', title: 'Density of the same mixtures by the suite\'s own GERG-2008 implementation', quantity: 'Density', unit: 'kg/m³', kind: 'benchmark', source: mixSrc('single-phase states of four gas mixtures'),
    columns: [C_MIX, C_COMP, C_T, C_P, { key: 'rho', label: 'Reference density', unit: 'kg/m³' }], rows: bRows('state', 2), target: 'rho', model: (r) => bAlt(r, 'gerg').rho ?? NaN, tolerance: { mape: 0.05 }, note: `${benchNote} Two independent implementations of the same equation of state (coefficients here from NIST teqp): mixing rules, reducing functions and departure terms must agree to round-off.` },
  { id: 'bench-mixture-density-pcsaft', title: 'Density of the same mixtures by PC-SAFT', quantity: 'Density', unit: 'kg/m³', kind: 'benchmark', source: mixSrc('single-phase states of four gas mixtures'),
    columns: [C_MIX, C_COMP, C_T, C_P, { key: 'rho', label: 'Reference density', unit: 'kg/m³' }], rows: bRows('state', 2), target: 'rho', model: (r) => bAlt(r, 'pcsaft').rho ?? NaN, tolerance: { mape: 2.5 }, note: `${benchNote} PC-SAFT with the Gross & Sadowski parameters and the few binary parameters held.` },
  { id: 'bench-mixture-sound', title: 'Speed of sound of natural-gas mixtures (Peng–Robinson)', quantity: 'Speed of sound', unit: 'm/s', kind: 'benchmark', source: mixSrc('single-phase states of four gas mixtures'),
    columns: [C_MIX, C_COMP, C_T, C_P, { key: 'w', label: 'Reference speed of sound', unit: 'm/s' }], rows: bRows('state'), target: 'w', model: (r) => bCub(r).sound, tolerance: { mape: 6 }, note: `${benchNote} Derivative property of the cubic equation; errors grow in the dense phase.` },
  { id: 'bench-mixture-sound-gerg', title: 'Speed of sound of the same mixtures by the suite\'s own GERG-2008 implementation', quantity: 'Speed of sound', unit: 'm/s', kind: 'benchmark', source: mixSrc('single-phase states of four gas mixtures'),
    columns: [C_MIX, C_COMP, C_T, C_P, { key: 'w', label: 'Reference speed of sound', unit: 'm/s' }], rows: bRows('state', 2), target: 'w', model: (r) => bAlt(r, 'gerg').w ?? NaN, tolerance: { mape: 0.6 }, note: `${benchNote} The residual part is GERG-2008 in both codes; the ideal-gas heat capacity here is the kernel polynomial, not the GERG ideal-gas function, which accounts for the small difference.` },
  { id: 'bench-mixture-joule-thomson', title: 'Joule–Thomson coefficient of natural-gas mixtures (Peng–Robinson)', quantity: 'Joule–Thomson coefficient', unit: 'K/bar', kind: 'benchmark', source: mixSrc('single-phase states of four gas mixtures'),
    columns: [C_MIX, C_COMP, C_T, C_P, { key: 'jt', label: 'Reference coefficient', unit: 'K/bar' }], rows: bRows('state'), target: 'jt', model: (r) => bCub(r).jtBar, tolerance: { rmse: 0.05 }, note: `${benchNote} The coefficient that sets the temperature drop across chokes and along a gas flowline; includes the inversion at high pressure.` },
  { id: 'bench-mixture-joule-thomson-gerg', title: 'Joule–Thomson coefficient by the suite\'s own GERG-2008 implementation', quantity: 'Joule–Thomson coefficient', unit: 'K/bar', kind: 'benchmark', source: mixSrc('single-phase states of four gas mixtures'),
    columns: [C_MIX, C_COMP, C_T, C_P, { key: 'jt', label: 'Reference coefficient', unit: 'K/bar' }], rows: bRows('state', 2), target: 'jt', model: (r) => bAlt(r, 'gerg').jt ?? NaN, tolerance: { rmse: 0.01 }, note: benchNote },
  { id: 'bench-mixture-heat-capacity', title: 'Isobaric heat capacity of natural-gas mixtures (Peng–Robinson)', quantity: 'Heat capacity cp', unit: 'J/g/K', kind: 'benchmark', source: mixSrc('single-phase states of four gas mixtures'),
    columns: [C_MIX, C_COMP, C_T, C_P, { key: 'cp', label: 'Reference cp', unit: 'J/g/K' }], rows: bRows('state'), target: 'cp', model: (r) => bCub(r).cp / 1000, tolerance: { mape: 5 }, note: `${benchNote} Calorimetric property of the gas phase.` },
  { id: 'bench-mixture-enthalpy', title: 'Enthalpy differences of natural-gas mixtures: isobaric heating and isothermal expansion (Peng–Robinson)', quantity: 'Enthalpy difference', unit: 'kJ/kg', kind: 'benchmark', source: mixSrc('enthalpy of single-phase states of three gas mixtures'),
    columns: [C_MIX, { key: 'path', label: 'Path', type: 'text' }, { key: 'T1', label: 'T₁', unit: 'K' }, { key: 'P1', label: 'P₁', unit: 'bar' }, { key: 'T2', label: 'T₂', unit: 'K' }, { key: 'P2', label: 'P₂', unit: 'bar' }, { key: 'dh', label: 'Reference h₂ − h₁', unit: 'kJ/kg' }], rows: bRows('dh'), target: 'dh', model: (r) => hOf(r.y, r.T2, r.P2, 'pr') - hOf(r.y, r.T1, r.P1, 'pr'), tolerance: { rmse: 6 }, note: `${benchNote} The quantities of a heat balance: sensible heat at constant pressure and the isothermal enthalpy change that drives Joule–Thomson cooling.` },
  { id: 'bench-mixture-enthalpy-gerg', title: 'Enthalpy differences by the suite\'s own GERG-2008 implementation', quantity: 'Enthalpy difference', unit: 'kJ/kg', kind: 'benchmark', source: mixSrc('enthalpy of single-phase states of three gas mixtures'),
    columns: [C_MIX, { key: 'path', label: 'Path', type: 'text' }, { key: 'T1', label: 'T₁', unit: 'K' }, { key: 'P1', label: 'P₁', unit: 'bar' }, { key: 'T2', label: 'T₂', unit: 'K' }, { key: 'P2', label: 'P₂', unit: 'bar' }, { key: 'dh', label: 'Reference h₂ − h₁', unit: 'kJ/kg' }], rows: bRows('dh'), target: 'dh', model: (r) => hOf(r.y, r.T2, r.P2, 'gerg') - hOf(r.y, r.T1, r.P1, 'gerg'), tolerance: { rmse: 1.5 }, note: `${benchNote} Residual enthalpy from GERG-2008, ideal-gas part from the kernel heat-capacity polynomials.` },
  { id: 'bench-mixture-compressibility', title: 'Isothermal compressibility of natural-gas mixtures and of the compressed light oil (Peng–Robinson)', quantity: 'Isothermal compressibility', unit: '1/bar', kind: 'benchmark', source: mixSrc('single-phase states; compressed liquid of the six-component light oil at 280–360 K, 200–600 bar'),
    columns: [C_MIX, C_COMP, C_T, C_P, { key: 'kT', label: 'Reference compressibility', unit: '1/bar' }], rows: [...bRows('state', 2), ...bLiq], target: 'kT', model: (r) => bCub(r).kT, tolerance: { mape: 12 }, note: `${benchNote} Gas compressibility is reproduced within a few per cent; the compressibility of the under-saturated oil is the weak point of a volume-translated cubic equation.` },
  { id: 'bench-oil-density', title: 'Compressed-liquid density of the six-component light oil, 280–360 K, 200–600 bar (Peng–Robinson)', quantity: 'Density', unit: 'kg/m³', kind: 'benchmark', source: mixSrc('compressed liquid above the bubble point'),
    columns: [C_MIX, C_COMP, C_T, C_P, { key: 'rho', label: 'Reference density', unit: 'kg/m³' }], rows: bLiq, target: 'rho', model: (r) => bCub(r).rho, tolerance: { mape: 5 }, note: `${benchNote} Under-saturated oil at reservoir-type pressures.` },
  { id: 'bench-oil-density-pcsaft', title: 'Compressed-liquid density of the light oil by PC-SAFT and speed of sound cross-check', quantity: 'Density', unit: 'kg/m³', kind: 'benchmark', source: mixSrc('compressed liquid above the bubble point'),
    columns: [C_MIX, C_COMP, C_T, C_P, { key: 'rho', label: 'Reference density', unit: 'kg/m³' }], rows: bLiq, target: 'rho', model: (r) => bAlt(r, 'pcsaft').rho ?? NaN, tolerance: { mape: 3 }, note: benchNote },
  // ---- gas hydrates
  { id: 'hydrate-pure', title: 'Hydrate dissociation of single guests in pure water below 100 bar (CH₄, C₂H₆, C₃H₈, i-C₄H₁₀, CO₂, H₂S)', quantity: 'Dissociation temperature', unit: 'K', kind: 'experiment', source: hydCite(['pure']),
    columns: hydCols, rows: withGas(HYDRATE_DATA.pure), target: 'T', model: hydModelT(KIHARA), tolerance: { rmse: 1.3, bias: 1.1 }, note: 'Default model: van der Waals–Platteeuw with Langmuir constants from the Kihara cell potential, guest fugacities from Peng–Robinson and Henry\'s-law gas solubility. The model errs warm (hydrate predicted slightly before it is measured). The H₂S Kihara parameters were regressed in this work to Selleck et al. (1952): rows of that publication show the fit, not a validation.' },
  { id: 'hydrate-high-pressure', title: 'Hydrate dissociation at 100–960 bar: methane, nitrogen, CO₂, mixtures and natural gases', quantity: 'Dissociation temperature', unit: 'K', kind: 'experiment', source: hydCite(['highP']),
    columns: hydCols, rows: withGas(HYDRATE_DATA.highP), target: 'T', model: hydModelT(KIHARA), tolerance: { rmse: 1.6, bias: 0.8 }, note: 'Deep-water shut-in pressures, model value without allowance. Unbiased up to about 300 bar; above, the model drifts cold (about −0.4 K at 350 bar, −1.4 K at 500 bar, −1.9 K at 600–730 bar for methane and nitrogen): the published curve therefore carries 0.004 K per bar above 250 bar.' },
  { id: 'hydrate-mixtures', title: 'Hydrate dissociation of gas mixtures and natural gases below 100 bar', quantity: 'Dissociation temperature', unit: 'K', kind: 'experiment', source: hydCite(['mix']),
    columns: hydCols, rows: withGas(HYDRATE_DATA.mix), target: 'T', model: hydModelT(KIHARA), tolerance: { rmse: 2.8, bias: 1.5 }, note: 'Structure selection (I or II) is part of the prediction. Six points of propane-rich methane–H₂S–propane ternaries (Schroeter et al. 1983) are predicted about 10 K too warm by both constants sets and dominate the root-mean-square error; without them it is 1.2 K. Inside the Kihara set, n-butane (both structures) and CO₂ in structure II keep the Munck constants because the Kihara parameters taken from the open NeqSim database mis-predicted methane–n-butane (+5 to +8 K) and CO₂–propane (down to −12 K): that choice was made on points of this set, which is therefore not blind for those two binaries.' },
  { id: 'hydrate-munck', title: 'Alternative Munck constants on the same pure, high-pressure and mixture points', quantity: 'Dissociation temperature', unit: 'K', kind: 'experiment', source: hydCite(hydAll),
    columns: hydCols, rows: withGas(hydAll.flatMap((k) => HYDRATE_DATA[k].filter((_, i) => i % 3 === 0))), target: 'T', model: hydModelT(LANGMUIR), tolerance: { rmse: 2, bias: 1 }, note: 'The selectable alternative (the default of earlier versions). On all points the Munck set is biased cold (−0.4 K below 100 bar, −1.2 K above, −2.1 K above 300 bar: hydrate forms warmer than predicted), which is why the Kihara set is the default.' },
  { id: 'hydrate-methanol', title: 'Methanol-inhibited hydrates (0.6–50 wt % methanol): methane, CO₂, H₂S, propane, mixtures and natural gases', quantity: 'Dissociation temperature', unit: 'K', kind: 'experiment', source: hydCite(['meoh']),
    columns: hydCols, rows: withGas(HYDRATE_DATA.meoh), target: 'T', model: hydModelT(KIHARA), tolerance: { rmse: 2, bias: 1.6 }, note: 'Model value without allowance. Water activity from NRTL regressed in this work to vapour–liquid, limiting-activity-coefficient and ice-line data (ACTIVITY_CALIBRATION) — not to hydrate data. The model predicts less depression than most sets measure (it errs warm by 1–3 K at 10–35 wt %); the largest cold deviations are −1.9 K (7 wt %, methane–propane) and, in the wider database, −2.7 K at 50 wt % and 233–255 K (Robinson & Ng 1986), an extrapolation of the activity model below its calibration range.' },
  { id: 'hydrate-meg', title: 'MEG-inhibited hydrates (5–50 wt % MEG): methane, ethane, propane, CO₂, H₂S and mixtures', quantity: 'Dissociation temperature', unit: 'K', kind: 'experiment', source: hydCite(['meg']),
    columns: hydCols, rows: withGas(HYDRATE_DATA.meg), target: 'T', model: hydModelT(KIHARA), tolerance: { rmse: 1.6, bias: 1 }, note: 'Model value without allowance. Water activity from NRTL regressed in this work to isothermal and isobaric vapour–liquid data of water–MEG (ACTIVITY_CALIBRATION) — not to hydrate data. From 5 to 40 wt % every set is matched within about +1.4 K (warm). At 50 wt % the two laboratories that measured it disagree by 4.6 K: the model lies between them, 2.3 K warm against methane (Robinson & Ng 1986) and ethane (Ng & Robinson 1985) and 2.3 K cold against methane–propane (Song & Kobayashi 1989). The allowance of the published curve covers that cold deviation (next set).' },
  { id: 'hydrate-inhibited-published', title: 'Published (conservative) curve on all methanol and MEG points: model value plus allowances', quantity: 'Dissociation temperature', unit: 'K', kind: 'experiment', source: hydCite(['meoh', 'meg']),
    columns: hydCols, rows: withGas([...HYDRATE_DATA.meoh, ...HYDRATE_DATA.meg]), target: 'T', model: hydModelT(KIHARA, false, true), tolerance: { rmse: 2.4, bias: 1.9 }, note: 'What the other suites receive: the model value plus 0–12 % of the inhibitor depression (zero up to 30 wt %, 12 % from 50 wt %) plus 0.004 K per bar above 250 bar. The bias is positive by design (the curve must not be optimistic); the tolerance bounds how conservative it is.' },
  { id: 'hydrate-brine', title: 'Hydrates in sodium-chloride brines (3–24 wt % NaCl): methane, ethane, propane, i-butane, CO₂, H₂S and mixtures', quantity: 'Dissociation temperature', unit: 'K', kind: 'experiment', source: hydCite(['nacl']),
    columns: hydCols, rows: withGas(HYDRATE_DATA.nacl), target: 'T', model: hydModelT(KIHARA), tolerance: { rmse: 1.4, bias: 1 }, note: 'Water activity of the brine from the Pitzer model with the PHREEQC parameters. The bias equals that of the same gases in pure water: the salt response itself is within about 0.3 K up to 24 wt %.' },
  { id: 'hydrate-structure-h', title: 'Structure-H hydrate of methane with methylcyclohexane', quantity: 'Dissociation temperature', unit: 'K', kind: 'experiment', source: hydCite(['sH']),
    columns: [{ key: 'src', label: 'Source', type: 'text' }, C_P, { key: 'T', label: 'Measured dissociation temperature', unit: 'K' }], rows: withGas(HYDRATE_DATA.sH), target: 'T', model: hydModelT(KIHARA, true), tolerance: { rmse: 0.5 }, note: 'The two lumped structure-H constants were regressed to Nakamura et al. (2003, 19 points, not shown here); the rows come from three other publications, which were not used.' },
  { id: 'nacl-osmotic', title: 'Osmotic coefficient of NaCl(aq) at 25 °C, 0.01–6 mol/kg', quantity: 'Osmotic coefficient', unit: '–', kind: 'experiment', source: SOURCES.hamerWu,
    columns: [{ key: 'm', label: 'Molality', unit: 'mol/kg' }, { key: 'phi', label: 'Evaluated osmotic coefficient' }], rows: NACL_25C, target: 'phi', model: (r) => osmoticNaCl('pitzer', r.m), tolerance: { mape: 0.3 }, note: 'Pitzer model with the Na–Cl parameters of the PHREEQC database: the water activity of brines follows directly from this coefficient.' },
  { id: 'water-content-cpa', title: 'Water content of methane in equilibrium with liquid water, 284–323 K, 50–200 bar (CPA)', quantity: 'Water in gas', unit: 'mg/Sm³', kind: 'experiment', source: SOURCES.frost,
    columns: [C_T, C_P, { key: 'wc', label: 'Measured water content', unit: 'mg/Sm³' }], rows: C1_WATER_VLE.map((r) => ({ ...r, wc: +wcOfY(r.yw).toPrecision(4) })), target: 'wc', model: (r) => cpaWater(c1(), [1], r.P, r.T - KEL).wc * 1e6, tolerance: { mape: 20 }, note: 'Cubic-plus-association with the published water parameters and the water–methane interaction of the NeqSim table. The measurements themselves scatter by 10–20 % between neighbouring points.' },
  { id: 'water-content-bukacek', title: 'Water content of methane: Bukacek correlation used by the flow suites', quantity: 'Water in gas', unit: 'mg/Sm³', kind: 'experiment', source: SOURCES.frost,
    columns: [C_T, C_P, { key: 'wc', label: 'Measured water content', unit: 'mg/Sm³' }], rows: C1_WATER_VLE.map((r) => ({ ...r, wc: +wcOfY(r.yw).toPrecision(4) })), target: 'wc', model: (r) => waterContent(r.P, r.T - KEL) * 1e6, tolerance: { mape: 20 }, note: 'The kernel correlation (sweet gas, fresh water) on the same points.' },
  { id: 'methane-solubility', title: 'Solubility of methane in water at 283 and 298 K, 12–100 bar', quantity: 'Dissolved methane', unit: 'mol/kg', kind: 'experiment', source: SOURCES.bottger,
    columns: [C_T, C_P, { key: 'm', label: 'Measured molality', unit: 'mol/kg' }], rows: C1_SOLUBILITY, target: 'm', model: (r) => { const f = c1(), x = gasSolubility({ C1: Math.exp(eosPhase(f, [1], r.P, r.T, 'vapour').lnphi[0]) * r.P }, r.P, r.T - KEL).x.C1; return (x / (1 - x)) * (1000 / MW_W); }, tolerance: { mape: 10 }, note: 'Henry\'s law with the Peng–Robinson gas fugacity and the Poynting correction — the dissolved-gas term of the hydrate model and of the water balance.' },
  { id: 'methanol-in-gas', title: 'Methanol and water in methane over aqueous methanol (280–323 K, 50–180 bar)', quantity: 'Methanol in gas', unit: 'mol ppm', kind: 'experiment', source: SOURCES.frost,
    columns: [C_T, C_P, { key: 'xw', label: 'Water mole fraction in the liquid' }, { key: 'ppm', label: 'Measured methanol in gas', unit: 'mol ppm' }], rows: C1_MEOH_WATER_VLE.map((r) => ({ ...r, ppm: +(r.ym * 1e6).toPrecision(4) })), target: 'ppm',
    model: (r) => { const wI = ((1 - r.xw) * 32.042) / ((1 - r.xw) * 32.042 + r.xw * MW_W); return cpaWater(c1(), [1], r.P, r.T - KEL, { inhId: 'MeOH', inhWt: 100 * wI }).yInh * 1e6; }, tolerance: { mape: 35 }, note: 'Inhibitor partitioning into the gas by cubic-plus-association (methanol two-site, water four-site, cross-association by the CR-1 rule). Only nine points, with a stated experimental uncertainty of 10–20 %; the comparison is an order-of-magnitude check of the methanol loss.' },
];
// ---- provenance of every constant set ----------------------------------------------------------------------------------------------
const RET = '2026-10-09', GH = 'https://raw.githubusercontent.com/', TML = 'https://trc.nist.gov/ThermoML/';
const PV = (item, used, source, url, status, note) => ({ item, used, source, url, retrieved: RET, status, note });
/** Where each constant set of this suite and of the thermodynamic kernel was checked ('verified' only when the source was opened and the numbers compared). */
export const PROVENANCE = [
  PV('Critical temperature, critical pressure, acentric factor, molar mass and critical volume of N2, CO2, H2S, C1–nC6', 'kernel COMPONENTS: every equation-of-state and viscosity calculation', 'CoolProp fluid files (reference equation of state of each fluid)', GH + 'CoolProp/CoolProp/master/dev/fluids/Methane.json (and Ethane, n-Propane, IsoButane, n-Butane, Isopentane, n-Pentane, n-Hexane, Nitrogen, CarbonDioxide, HydrogenSulfide)', 'corrected', 'All eleven components replaced by the source values. Largest changes: H2S Tc 373.4 → 373.1 K, Pc 89.63 → 90.0 bar, ω 0.0942 → 0.1005; n-hexane Pc 30.25 → 30.441 bar; i-butane Vc 262.7 → 257.75 cm³/mol; nitrogen ω 0.0377 → 0.0372; the others moved by less than 0.2 %.'),
  PV('Ideal-gas heat-capacity polynomials', 'kernel COMPONENTS.cp: enthalpy, entropy, heat capacity, speed of sound', 'Poling, Prausnitz & O\'Connell, The Properties of Gases and Liquids, 5th ed., as tabulated in the open `chemicals` library (PolingDatabank.tsv)', GH + 'CalebBell/chemicals/master/chemicals/Heat%20Capacity/PolingDatabank.tsv', 'corrected', 'The cubic polynomials held before deviated from the source quartics by up to 3.7 % at 250 K and 9 % at 200 K (methane, ethane, n-pentane). Replaced by cubic least-squares fits of the source over 180–560 K (largest residual 0.10 %). Checked further against the NIST heat capacities (data set nist-heat-capacity).'),
  PV('Peng–Robinson binary interaction parameters of the light components', 'kernel KIJ', 'ChemSep interaction-parameter table (DECHEMA-regressed) distributed with the open `thermo` library', GH + 'CalebBell/thermo/master/thermo/Interaction%20Parameters/ChemSep/pr.json', 'corrected', 'Replaced by the source: N2–C2 0.010 → 0.0533, N2–H2S 0.13 → 0.1652, CO2–C1 0.105 → 0.0978, CO2–C3 0.125 → 0.1315, N2–C1 0.025 → 0.0289, H2S–C2 0.085 → 0.0952, H2S–C3 0.08 → 0.0878; 21 hydrocarbon–hydrocarbon pairs that were zero now carry the source values (|kij| ≤ 0.04). H2S–methane (0.07) is not in the table and is unchanged.'),
  PV('Default interaction parameters with C7+ (N2 0.10, CO2 0.115, H2S 0.06; methane–C7+ 0.14·SG − 0.0668)', 'kernel kijOf', '—', '', 'unverified', 'Customary characterisation defaults; no openly readable source was found. They are tuning parameters (binary-interaction scale and override table).'),
  PV('Peng–Robinson volume-shift parameters of the light components, C7+ shift 1 − 2.258/M^0.1823, Péneloux SRK shift', 'kernel makeFluid: liquid density', '—', '', 'unverified', 'Jhaveri–Youngren-type values written from the literature; not found in an openly readable source. Their effect is validated instead: liquid densities of seven fluids against NIST within 1.6 % on average (data set nist-liquid-density). A volume-shift multiplier is now a calibration parameter.'),
  PV('Parachors of the light components', 'kernel COMPONENTS.par: gas–oil interfacial tension', 'NeqSim component database', GH + 'equinor/neqsim/master/src/main/resources/data/COMP.csv', 'verified', 'Within 1.1 % of the source for all eleven components (largest: propane 150.3 here, 151.9 in the source; n-pentane 231.5 / 233.9). Pure-component surface tension against NIST: data set nist-surface-tension. The C7+ relation 59.3 + 2.34·M is unverified.'),
  PV('Equation-of-state constants (Ωa, Ωb, m(ω) of Peng–Robinson 1976/1978, Soave)', 'kernel EOS', 'open `thermo` library, eos.py (classes PR, PR78, SRK, RK)', GH + 'CalebBell/thermo/master/thermo/eos.py', 'verified', 'Identical: Peng–Robinson Ωa 0.45723552892, Ωb 0.07779607390, κ = 0.37464 + 1.54226ω − 0.26992ω² up to ω = 0.491 and 0.379642 + 1.48503ω − 0.164423ω² + 0.016666ω³ above; Soave/Redlich–Kwong Ωa 0.42748023354, Ωb 0.08664034996, m = 0.480 + 1.574ω − 0.176ω².'),
  PV('Lohrenz–Bray–Clark viscosity polynomial and Stiel–Thodos dilute-gas viscosity', 'kernel viscosityLBC', 'open `chemicals` library (viscosity.py: Lorentz_Bray_Clarke, Stiel_Thodos) and NeqSim LBCViscosityMethod.java', GH + 'equinor/neqsim/master/src/main/java/neqsim/physicalproperties/methods/commonphasephysicalproperties/viscosity/LBCViscosityMethod.java', 'verified', 'All coefficients agree with NeqSim (0.1023, 0.023364, 0.058533, −0.040758, 0.0093324; 34e-5·Tr^0.94 and 17.78e-5·(4.58Tr − 1.67)^0.625). The `chemicals` file carries 0.0093724 for the last coefficient — a difference between the two sources, not changed here.'),
  PV('Stiel–Thodos dense-fluid thermal conductivity', 'conductivityStielThodos()', 'open `chemicals` library (thermal_conductivity.py: Stiel_Thodos_dense)', GH + 'CalebBell/chemicals/master/chemicals/thermal_conductivity.py', 'verified', 'All constants identical (210, 1.22e-2/0.535, 1.14e-2/0.67/1.069, 2.60e-3/1.155/2.016).'),
  PV('Lee–Kesler modified BWR constants of the simple and reference fluids', 'leeKesler(), leeKeslerZ()', 'thermopack (SINTEF/NTNU) leekesler.f90', GH + 'thermotools/thermopack/main/src/leekesler.f90', 'verified', 'All 24 constants and ω_ref = 0.3978 identical.'),
  PV('Lee–Kesler / Kesler–Lee acentric-factor correlations', 'kernel pseudoProps, acentricLK()', 'open `chemicals` library (acentric.py: LK_omega) and NeqSim TBPfractionModel.java', GH + 'CalebBell/chemicals/master/chemicals/acentric.py', 'verified', 'Identical (5.92714, 6.09648, 1.28862, 0.169347; 15.2518, 15.6875, 13.4721, 0.43577; heavy branch −7.904, 0.1352, 0.007465, 8.359, 1.408, 0.01063). NeqSim writes 6.09649 for the second constant.'),
  PV('Kesler–Lee critical temperature and pressure of petroleum fractions', 'kernel pseudoProps', 'NeqSim TBPfractionModel.java (class LeeKesler, kelvin / MPa form)', GH + 'equinor/neqsim/master/src/main/java/neqsim/thermo/characterization/TBPfractionModel.java', 'verified', 'After conversion from °R/psia to K/MPa every coefficient agrees within 0.01 % except one (9.9099 in the source against 9.9010 from the constant held here, 0.09 %).'),
  PV('Twu (1984) critical-property correlations', 'critTwu()', 'NeqSim TBPfractionModel.java (class TwuModel)', GH + 'equinor/neqsim/master/src/main/java/neqsim/thermo/characterization/TBPfractionModel.java', 'verified', 'All coefficients of the n-alkane reference (Tc, Pc, Vc, SG) and of the three perturbation functions agree after unit conversion (largest difference 0.004 %).'),
  PV('Riazi–Daubert (1980) critical-property correlations', 'critRiaziDaubert()', '—', '', 'unverified', 'The open implementations found use the 1987 molar-mass form, not the 1980 boiling-point form held here. Checked only against the critical temperature of n-decane (within 3 K).'),
  PV('Søreide boiling-point correlation', 'kernel tbSoreide', 'NeqSim TBPfractionModel.java', GH + 'equinor/neqsim/master/src/main/java/neqsim/thermo/characterization/TBPfractionModel.java', 'verified', 'Identical. The Søreide specific-gravity relation 0.2855 + Cf(M − 66)^0.13 is unverified.'),
  PV('Won (1986) melting temperature and enthalpy of fusion of paraffins', 'wonFusion()', 'NeqSim ComponentWonWax.java', GH + 'equinor/neqsim/master/src/main/java/neqsim/thermo/component/ComponentWonWax.java', 'verified', 'Identical (374.5 + 0.02617·M − 20172/M; 0.1426·M·Tf cal/mol). The Pedersen wax-forming fraction is unverified.'),
  PV('Munck et al. (1988) Langmuir constants and reference properties of the empty hydrate lattice', 'LANGMUIR, HYDRATE_STRUCTURES: selectable Munck hydrate set; n-butane and CO2 (structure II) inside the default set', 'U.S. Patent 6,871,118 B2 (hydrate module after Munck et al. 1988: tables of A, B and of the reference properties; patent text, public domain)', 'https://patents.google.com/patent/US6871118B2/en', 'corrected', 'All eight guests compared with the table of the source. Nitrogen, H2S, methane, ethane, propane, i-butane and n-butane and the reference properties (Δμ⁰ 1264 / 883 J/mol, ΔH⁰ −4858 / −5201 J/mol, ΔCp 39.16 J/mol/K) are identical. Carbon dioxide differed and was replaced by the source: structure I 0.00588e-3 / 5410 and 3.36e-3 / 3202 → 0.2474e-3 / 3410 and 42.46e-3 / 2813; structure II 0.0846e-3 / 3602 and 846e-3 / 2030 → 0.0845e-3 / 3615 and 851e-3 / 2025 (K/atm, K). Effect on the Munck option: methane–CO2 mixtures bias +1.2 → 0.0 K, CO2 with methanol or MEG +2.0 / +2.4 → +1.6 / +0.9 K.'),
  PV('Kihara parameters, cavity radii and coordination numbers; Sloan & Koh reference properties', 'KIHARA, kiharaC(): default hydrate model and structure H', 'Avaji et al. (2023) Fluid Phase Equilib. 567, 113716 (open manuscript, Tables 1, 2 and 8); NeqSim COMP.csv and ComponentHydrate.java', 'https://bradscholars.brad.ac.uk/server/api/core/bitstreams/d47ed2ca-b934-4bc8-ae3b-a8631c4e0d8d/content', 'verified', 'C1–iC4 from Avaji et al.; N2 and CO2 from NeqSim; radii 3.95/4.33 and 3.91/4.73 Å, coordination 20/24 and 20/28; structure-H cavity data from Herri et al., https://hal-emse.ccsd.cnrs.fr/emse-00724388v1/document, Table 1. Now the default set, chosen on 476 measured points: bias +0.7 K below 100 bar, +0.3 K at 100–960 bar, +1.1 K for mixtures (Munck: −0.3, −1.2, +0.1 K; −2.1 K above 300 bar). Three departures from the published parameters, all made in this work: (1) H2S (a 0.36 Å, σ 3.1688 Å, ε/k 206.61 K) regressed to the 15 points of Selleck et al. (1952) and checked on three other publications (19 points, bias +0.7 K); (2) n-butane keeps the Munck constants — the NeqSim Kihara values predicted methane–n-butane 5–8 K too warm; (3) CO2 in structure II keeps the Munck constants — the NeqSim values predicted CO2–propane up to 12 K too cold. (2) and (3) were decided on the mixture data set itself.'),
  PV('Structure-H effective reference constants (Δμ⁰ 528.7 J/mol, Δh⁰ −5624 J/mol)', 'HYDRATE_STRUCTURES.sH', 'regressed in this work to NIST hydrate database set 275', 'https://gashydrates.nist.gov/hydrate-browser/dataset_jsons.json', 'unverified', 'No published structure-H parameter set could be read from an open source. Two lumped constants (large-cavity occupation by methylcyclohexane included) were fitted to 19 points of Nakamura et al. (2003) and reproduce three independent sets within 0.3 K on average. Valid for methane with a methylcyclohexane-type former only.'),
  PV('NRTL / Wilson / UNIQUAC parameters, methanol–water and MEG–water (temperature-dependent, e = a + b/T)', 'ACTIVITY_PARAMS: water activity of inhibited aqueous phases', 'regressed in this work to non-hydrate measurements in the NIST ThermoML archive and to the ice line of aqueous methanol (record: ACTIVITY_CALIBRATION in data/ref/pvt.js)', TML + '10.1021/je0342522.json', 'corrected', 'Calibration (never hydrate data). MEG: P–x at 333.15 and 353.15 K (Horstmann et al. 2004, 80 points, NRTL bias −0.2 %, r.m.s. 0.6 % in pressure), boiling points at 50–101 kPa (Chouireb et al. 2018, 49 points, 4.6 %; Kamihama et al. 2012, 19 points, 2.6 %). Methanol: P–x at 323.15 K (Bernatová et al. 2006, 7 points, 1.2 %), limiting activity coefficient at 273–323 K (Vrbka et al. 2005, 8 points, 4.2 %), ice line 4–34 wt % (8 points, within 1.5 K, model warm). Validation: the hydrate sets. The earlier MEG set had no temperature dependence (fitted to ice points of unstated origin and used at all temperatures) and gave 3–6 K too much depression at 40–50 wt % on methane–propane; the new set extrapolates from 333 K to ice points of −3.3, −7.5, −13.1, −20.9, −33.4 °C at 10–50 wt % without having been fitted to them. Remaining uncertainty: the three model forms fit the vapour–liquid data equally but differ at 265 K and 50 wt % MEG by ±0.02 in ln a_w (about ±1.5 K of hydrate temperature); NRTL, in the middle, is the default. Methanol below 240 K and above 40 wt % is an extrapolation.'),
  PV('Pitzer parameters of NaCl (β0, β1, Cφ and their temperature functions)', 'PITZER_NACL: water activity of brines', 'USGS PHREEQC database pitzer.dat', GH + 'usgs-coupled/phreeqc3/master/database/pitzer.dat', 'corrected', 'β0 0.0765 → 0.07534, β1 0.2664 → 0.2769, Cφ 0.00127 → 0.00148, and the linear temperature slopes written from memory replaced by the six-term PHREEQC temperature function. Against Hamer & Wu (1972): osmotic coefficient within 0.06 % on average, activity coefficient 0.10 % (0.23 % before); methane-hydrate temperatures in NaCl brine (De Roo et al.) bias −0.33 K (−0.62 K before).'),
  PV('Debye–Hückel constants from the density and permittivity of water', 'debyeHuckel()', 'Malmberg, C. G. & Maryott, A. A. (1956). Dielectric constant of water from 0° to 100 °C. J. Res. Natl. Bur. Stand. 56, 1–8 (U.S. Government work), Table 2', 'https://nvlpubs.nist.gov/nistpubs/jres/56/jresv56n1p1_A1b.pdf', 'verified', 'The polynomial 87.740 − 0.40008 t + 9.398e-4 t² − 1.410e-6 t³ reproduces the tabulated values of the source at 0, 10, 20 and 25 °C (87.740, 83.832, 80.103, 78.304) to the last digit. Aφ = 0.3920 at 25 °C against the customary 0.3915.'),
  PV('Henry constants of gases in water at 25 °C and their temperature coefficients', 'HENRY: gas solubility, dissolved-gas correction of the hydrate model', 'Sander compilation as held in the open `thermo` library (Sander_henry_T_dep.json)', GH + 'CalebBell/thermo/master/thermo/Interaction%20Parameters/Sander_henry_T_dep.json', 'verified', 'Within 11 % of the source at 25 °C for CH4 (1.4e-3 here, 1.27e-3 mol/kg/bar in the file), C2H6, C3H8, n-C4H10, N2, CO2 and H2S; i-butane differs by a factor of two (8.4e-4 here, 1.6e-3 in the file) and is left unchanged pending a second source. Methane solubility against Böttger et al. (2016): 7 % (data set methane-solubility). Partial molar volumes and Setschenow constants are unverified.'),
  PV('PC-SAFT universal constants and pure-component parameters', 'pcsaftModel()', 'NIST teqp (PCSAFT.cpp, constants of Gross & Sadowski 2001) and Clapeyron.jl database (PCSAFT_like.csv, PCSAFT_unlike.csv)', GH + 'ClapeyronThermo/Clapeyron.jl/master/database/SAFT/PCSAFT/PCSAFT_like.csv', 'verified', 'Transcribed from the sources in this pass; reproduces NIST densities of seven fluids within 1.1 % on average and vapour pressures within 0.7 % except CO2 (up to 8 %: no quadrupole term). C7+ pseudo-components: n-alkane sets interpolated in molar mass with the segment diameter matched to the specific gravity (this work).'),
  PV('Cubic-plus-association parameters of water, methanol and MEG; water–gas interaction parameters', 'cpaModel(), cpaWater()', 'NeqSim component and interaction databases (published sets of Kontogeorgis and co-workers)', GH + 'equinor/neqsim/master/src/main/resources/data/INTER.csv', 'verified', 'Transcribed from the sources in this pass (water 4C: a0 0.12277 Pa·m⁶/mol², b 14.515 cm³/mol, c1 0.67359, ε 16655 J/mol, β 0.0692). Water content of methane against Frost et al. (2014): 17 % (data set water-content-cpa).'),
  PV('GERG-2008 coefficients (eleven components, 55 binary pairs, 7 departure functions)', 'gergModel()', 'NIST teqp, GERG.hpp', SOURCES.gerg.url, 'verified', 'Extracted by script from the source file into data/ref/pvt.js; reproduces NIST densities within 0.03 % on average (0.3 % at most) and vapour pressures within 0.2 % (data set gerg-density).'),
  PV('Pedersen corresponding-states viscosity: mixing rules and methane reference correlation (Hanley et al. 1975)', 'viscosityPedersen(), methaneViscosity()', 'NeqSim PFCTViscosityMethodMod86.java', GH + 'equinor/neqsim/master/src/main/java/neqsim/physicalproperties/methods/commonphasephysicalproperties/viscosity/PFCTViscosityMethodMod86.java', 'verified', 'All constants transcribed from the source; methane viscosity against NIST within 1.1 % on average. The tanh(T − 91 K) switch between the two dense-fluid branches is written from the published description and unverified.'),
  PV('Flory–Huggins asphaltene model defaults (molar mass 750 g/mol, density 1100 kg/m³, solubility parameter 21 MPa^0.5)', 'asphalteneFH(), asphalteneCurve()', 'NeqSim FloryHugginsAsphalteneModel.java', GH + 'equinor/neqsim/master/src/main/java/neqsim/pvtsimulation/flowassurance/FloryHugginsAsphalteneModel.java', 'verified', 'Equation and default values as in the source. The temperature coefficient of the solubility parameter (1.07e-3 per K, Hirschberg et al. 1984) is unverified. All four are inputs.'),
  PV('de Boer screening boundaries and colloidal-instability thresholds', 'asphalteneScreen()', '—', '', 'unverified', 'Approximate digitisation of the published plot; thresholds 0.7 / 0.9 written from the literature. The only open implementation found (NeqSim DeBoerAsphalteneScreening.java, its own quadratic digitisation) puts the boundaries about a factor of two lower in undersaturation (severe above 165 bar at 700 kg/m³ against 359 bar here), so the two digitisations do not confirm each other: screening only, and the Flory–Huggins equilibrium should be preferred.'),
  PV('Black-oil correlations (Standing, Vasquez–Beggs, Glasø, Beggs–Robinson), Lee–Gonzalez–Eakin gas viscosity', 'blackOil(), rsStanding(), viscosityLGE()', '—', '', 'unverified', 'Not compared with an open source in this pass; they are cross-checks beside the equation of state, not part of the published properties.'),
  PV('Kernel screening hydrate curve: gas-gravity correlation (8 coefficients) and gravity of the gas present at the pressure', 'kernel hydrateT0, hydrateScreening, fluidModel().hydrateT0 — used by the other suites until this suite has run', 'fitted in this work to the default hydrate model of this suite (not a literature correlation)', '', 'corrected', 'Replaces the Motiee (1991) coefficients. Least-squares fit of T = c0 + c1 L + c2 L² + c3 L³ + (c4 + c5 L + c7 L²) g + c6 g², L = ln P, g = ln(gravity/0.554), to the published model curve of 110 fluids (91 natural gases of gravity 0.57–0.96, 3 condensates, 16 oils) at 5–700 bara, shifted 1.0 K warm. Residual (correlation − model): gases +0.8 K (r.m.s. 1.5 K, −2.1 to +6.2 K, 74 % warm); condensates +0.8 K; oils +1.9 K (92 % warm) when the gravity of the gas present at the pressure is used, as fluidModel does. Motiee on the same gases: −2.0 K (−5.8 K at worst, 9 % warm), and −1.3 to −2.5 K against measured lean natural gases of gravity 0.58–0.65. A fixed gravity cannot represent composition: at one gravity the model spreads by about ±3 K (propane and butanes warm, CO2 and N2 cold), so individual gases can lie up to 2 K on the cold side. Reference case, fresh water at 25 / 100 / 300 bara: screening 12.6 / 20.4 / 24.7 °C, model 10.8 / 19.7 / 24.4 °C.'),
  PV('Kernel screening depression ΔT = −K ln x_w, K = 60 K (Nielsen–Bucklin form; original constant 72 K)', 'kernel hydrateDepression, inhibitorFor', 'constant set in this work from the hydrate model and the measured inhibitor data behind it', '', 'corrected', 'With 72 K the screening depression exceeded the model (structure-II gas, 30–100 bara) by up to 2.0 K at 50 wt % MEG, 4.4 K at 50 wt % methanol and 0.2 K for sea water — non-conservative. With 60 K it stays at or below the published model depression for 0–60 wt % MEG, 0–50 wt % methanol and 0–25 wt % NaCl (largest difference −0.05 K, verification check). The original Nielsen–Bucklin constant could not be read in an open source.'),
  PV('Bukacek water content; Hammerschmidt depression constants (1297 and 2222 K·g/mol)', 'kernel waterContent; hammerschmidt() (comparison column only)', '—', '', 'unverified', 'Coefficients not found in an openly readable source. The Bukacek correlation is demonstrated by validation data instead (measured water content of methane, 17 %, data set water-content-bukacek); Hammerschmidt is shown beside the model as a comparison and is not used in any published result.'),
  PV('Molar masses of the inhibitors', 'kernel INHIBITORS', 'NeqSim component database', GH + 'equinor/neqsim/master/src/main/resources/data/COMP.csv', 'verified', 'Methanol 32.042, MEG 62.068 (62.069), DEG 106.12 (106.122), TEG 150.17 (150.175), ethanol 46.069 g/mol. Densities are unverified.'),
  PV('Conservative allowances of the published hydrate curve (0–12 % of the inhibitor depression above 30 wt %; 0.004 K per bar above 250 bara)', 'inhibitorAllowance(), pressureAllowance(), publishedHydrateT(): outputs.hydrateCurve, inhibitor requirement', 'set in this work from the residuals of the hydrate data sets', '', 'verified', 'Demonstrated uncertainty of the model, applied towards the safe side. Inhibitor: the largest cold deviation of any inhibited set is −2.3 K at 50 wt % MEG on methane–propane (Song & Kobayashi 1989), 12 % of the model depression; the two laboratories that measured 50 wt % MEG differ by 4.6 K. Pressure: against methane and nitrogen points the model is unbiased to 300 bar and −0.4, −1.4, −1.9 K at 350, 500 and 600–730 bar. With the allowances the published curve has a bias of +1.6 K on the 68 inhibited points (coldest point −1.9 K, a 7 wt % methanol set whose uninhibited gas shows the same offset).'),
  PV('GERG-2008 binary reducing parameters: order of β and γ', 'gergModel()', 'CoolProp 6.7.0 (HEOS mixture backend) and NIST teqp GERG.hpp', SOURCES.gerg.url, 'corrected', 'The bundled arrays follow the teqp structure [βv, γv, βT, γT] but were read as [βT, γT, βv, γv]: pure fluids were unaffected, mixtures were wrong by up to 5.6 % in density (methane–propane 70/30 at 300 K, 200 bar). Found by the code-to-code benchmark and corrected; the suite now agrees with CoolProp within 0.01 % in density, 0.05 % in speed of sound and 0.2 kJ/kg in enthalpy differences on 53 mixture states.'),
];

// ---- interactive flash calculator (custom tab) ---------------------------------------------------------------------------------
function flashView(el, api) {
  const h = api.h, m = lastModel || { spec: mergeSpec(null), kij: [] }, v = api.values();
  let f; try { f = tunedFluid(m.spec, tuneOf(v, m.spec), m.kij).f; } catch (e) { el.append(h('p', { class: 'bad' }, e.message)); return; }
  const out = h('div'), st = { P: clamp(num(v.pRef, 100), 1, 1000), T: clamp(num(v.tSeabed, 4), -40, 250), P2: clamp(num(v.pArr, 25), 1, 1000) };
  const show = () => {
    try {
      const s = fluidState(f, st.P, st.T), q = stateHSV(f, st.P, st.T), jt = flashPH(f, st.P2, q.H), is = flashPS(f, st.P2, q.S), hy = hydrateTofP(f, st.P), tp = f.n > 1 ? stability(f, st.P, st.T + KEL) : { tpd: 0, stable: true };
      const g = s.gas, o = s.oil, n = (x, d = 4) => (Number.isFinite(x) ? +x.toPrecision(d) : '—');
      out.replaceChildren(
        api.kpiGrid([{ label: 'State', value: phaseName[s.phase] }, { label: 'Vapour fraction', value: n(100 * s.beta), unit: 'mol %' }, { label: 'Gas mass fraction', value: n(100 * s.wG), unit: '%' }, { label: 'Tangent-plane distance', value: n(tp.tpd, 3), status: tp.stable ? 'ok' : 'warn' },
          { label: 'Molar enthalpy', value: n(q.H / 1000), unit: 'kJ/mol' }, { label: 'Molar entropy', value: n(q.S), unit: 'J/mol/K' }, { label: 'Molar volume', value: n(q.V * 1e6), unit: 'cm³/mol' }, { label: 'Hydrate temperature (fresh water)', value: hy ? n(hy.T) : '—', unit: '°C', status: hy && hy.T > st.T ? 'warn' : 'ok' },
          { label: `Isenthalpic to ${n(st.P2)} bara`, value: jt ? n(jt.T) : '—', unit: '°C' }, { label: `Isentropic to ${n(st.P2)} bara`, value: is ? n(is.T) : '—', unit: '°C' }]),
        api.dataTable({ title: 'Phase properties', columns: ['Property', 'Gas', 'Oil', 'Unit'], rows: [['Density', 'rho', 1, 'kg/m³'], ['Viscosity', 'mu', 1e3, 'mPa·s'], ['Compressibility factor', 'Z', 1, '–'], ['Molar mass', 'MW', 1, 'g/mol'], ['Heat capacity cp', 'cp', 1, 'J/kg/K'], ['Heat capacity cv', 'cv', 1, 'J/kg/K'], ['Thermal conductivity', 'k', 1, 'W/m/K'], ['Joule–Thomson coefficient', 'jtBar', 1, 'K/bar'], ['Speed of sound', 'sound', 1, 'm/s'], ['Isothermal compressibility', 'kT', 1, '1/bar'], ['Enthalpy', 'h', 1e-3, 'kJ/kg']].map(([lbl, k, c, u]) => [lbl, g ? n(g[k] * c) : '—', o ? n(o[k] * c) : '—', u]).concat([['Interfacial tension', s.sigma !== null ? n(s.sigma * 1e3) : '—', '', 'mN/m']]) }),
        api.dataTable({ title: 'Phase compositions', columns: ['Component', 'Feed (mol %)', 'Liquid (mol %)', 'Vapour (mol %)', 'K = y/x'], rows: f.comps.map((c, i) => [c.name, n(100 * c.z), s.phase !== 'gas' ? n(100 * s.x[i]) : '—', s.phase !== 'oil' ? n(100 * s.y[i]) : '—', s.phase === 'two' && s.x[i] > 1e-300 ? n(s.y[i] / s.x[i]) : '—']) }));
    } catch (e) { out.replaceChildren(h('p', { class: 'bad' }, 'Flash failed: ' + e.message)); }
  };
  const slider = (label, key, min, max, step, unit) => { const val = h('strong', null, `${st[key]} ${unit}`), inp = h('input', { type: 'range', min, max, step, value: st[key], 'aria-label': label, oninput: (e) => { st[key] = +e.target.value; val.textContent = `${st[key]} ${unit}`; show(); } }); return h('label', { class: 'field' }, h('span', null, label + ' '), inp, ' ', val); };
  el.append(h('p', { class: 'summary' }, `Interactive PT, PH and PS flash of ${m.spec.name || 'the case fluid'} with the tuning of the model-setup tab (${f.eosId}, ${f.n} components). Run the suite first to pick up the current case fluid.`),
    h('fieldset', { class: 'group' }, h('legend', null, 'Conditions'), slider('Pressure', 'P', 1, 600, 1, 'bara'), slider('Temperature', 'T', -40, 250, 1, '°C'), slider('Expansion outlet pressure', 'P2', 1, 600, 1, 'bara')), out);
  show();
}

// ---- suite declaration -------------------------------------------------------------------------------------------------------------
const CAL_SAMPLE = [
  { calT: 40, calP: 250, psat: 192.4, rho: 688.1, rs: 172.3, mu: 1.32 },
  { calT: 40, calP: 120, psat: 191.1, rho: 714.5, rs: 109.1, mu: 1.9 },
  { calT: 60, calP: 300, psat: 209.2, rho: 678.8, rs: 170.4, mu: 1.16 },
  { calT: 60, calP: 150, psat: 206.9, rho: 689.7, rs: 123.3, mu: 1.35 },
  { calT: 60, calP: 60, psat: 205.7, rho: 742.6, rs: 49.77, mu: 2.57 },
  { calT: 90, calP: 300, psat: 226.1, rho: 655.9, rs: 174, mu: 0.828 },
  { calT: 90, calP: 200, psat: 228.2, rho: 656, rs: 147.8, mu: 0.785 },
  { calT: 90, calP: 120, psat: 225.2, rho: 699.4, rs: 83.79, mu: 1.41 },
  { calT: 90, calP: 40, psat: 226, rho: 745.5, rs: 26.37, mu: 2.46 },
  { calT: 120, calP: 280, psat: 241.1, rho: 631.2, rs: 169.7, mu: 0.547 },
  { calT: 120, calP: 160, psat: 238.5, rho: 662.2, rs: 102.8, mu: 0.814 },
  { calT: 120, calP: 80, psat: 241.5, rho: 704.2, rs: 47.34, mu: 1.48 },
];
const CAL_VALID = [
  { calT: 50, calP: 200, psat: 201.5, rho: 670.3, rs: 174.7, mu: 1.05 },
  { calT: 50, calP: 80, psat: 199.9, rho: 728.4, rs: 69.92, mu: 2.38 },
  { calT: 75, calP: 260, psat: 217.5, rho: 661.6, rs: 170.2, mu: 0.879 },
  { calT: 75, calP: 100, psat: 216.6, rho: 711.2, rs: 74.27, mu: 1.77 },
  { calT: 105, calP: 220, psat: 231.3, rho: 640.6, rs: 158.3, mu: 0.616 },
  { calT: 105, calP: 60, psat: 232.7, rho: 721.1, rs: 37.81, mu: 2.04 },
  { calT: 130, calP: 180, psat: 241.1, rho: 645, rs: 118, mu: 0.65 },
];
const CAL_SAMPLE2 = [
  {calT:40,calP:120,calSal:3.5,calInhWt:0,bo:1.309,rhoG:168.5,zG:0.6619,muG:0.0207,co:1.856,cp:2.232,sound:361.2,sigma:4.433,wc:681,rhoW:1029,aw:0.9814,hydT:16.95},
  {calT:60,calP:150,calSal:3.5,calInhWt:0,bo:1.374,rhoG:178.1,zG:0.7291,muG:0.02247,co:2.22,cp:2.355,sound:397.8,sigma:3.161,wc:1425,rhoW:1019,aw:0.98,hydT:18.04},
  {calT:60,calP:60,calSal:0,calInhWt:0,bo:1.173,rhoG:62.31,zG:0.8263,muG:0.01428,co:1.761,cp:2.276,sound:345.9,sigma:8.959,wc:2775,rhoW:990.3,aw:0.9996,hydT:14.63},
  {calT:90,calP:200,calSal:3.5,calInhWt:0,bo:1.471,rhoG:192,zG:0.8262,muG:0.02433,co:3.008,cp:2.454,sound:448.1,sigma:1.844,wc:3722,rhoW:1004,aw:0.9792,mmp:386.1,hydT:18.73},
  {calT:90,calP:120,calSal:10,calInhWt:0,bo:1.305,rhoG:116.6,zG:0.816,muG:0.01812,co:2.598,cp:2.429,sound:384.3,sigma:4.672,wc:5303,rhoW:1048,aw:0.9369,hydT:13.8},
  {calT:90,calP:40,calSal:0,calInhWt:20,bo:1.139,rhoG:34.7,zG:0.9191,muG:0.01394,co:2.01,cp:2.409,sound:372.1,sigma:10.17,wc:12610,rhoW:988.4,aw:0.9286,hydT:6.84},
  {calT:120,calP:160,calSal:3.5,calInhWt:30,bo:1.397,rhoG:134.6,zG:0.8676,muG:0.02025,co:3.356,cp:2.603,sound:421.3,sigma:2.812,wc:10250,rhoW:997.2,aw:0.8591,hydT:8.82},
  {calT:120,calP:80,calSal:0,calInhWt:40,bo:1.234,rhoG:65.58,zG:0.8972,muG:0.01609,co:2.896,cp:2.546,sound:388.2,sigma:6.405,wc:17220,rhoW:983.2,aw:0.822,hydT:1.74},
  {calT:20,calP:100,calSal:3.5,calInhWt:20,bo:1.275,rhoG:167.4,zG:0.5831,muG:0.02021,co:1.595,cp:2.18,sound:338,sigma:5.504,wc:249.3,rhoW:1050,aw:0.9074,hydT:10.87},
  {calT:10,calP:60,calSal:8,calInhWt:0,bo:1.2,rhoG:90.63,zG:0.6675,muG:0.01368,co:1.273,cp:2.07,sound:295.7,sigma:9.2,wc:196.1,rhoW:1070,aw:0.9522,hydT:11.43},
];
const CAL_VALID2 = [
  {calT:50,calP:80,calSal:3.5,calInhWt:0,bo:1.219,rhoG:93.84,zG:0.7575,muG:0.01541,co:1.769,cp:2.263,sound:335,sigma:7.491,wc:1400,rhoW:1022,aw:0.981,hydT:14.83},
  {calT:75,calP:100,calSal:0,calInhWt:10,bo:1.264,rhoG:104.4,zG:0.7947,muG:0.0166,co:2.229,cp:2.376,sound:363.1,sigma:5.854,wc:3564,rhoW:992.4,aw:0.9688,hydT:14.95},
  {calT:105,calP:60,calSal:5,calInhWt:0,bo:1.184,rhoG:51.26,zG:0.8945,muG:0.01468,co:2.42,cp:2.444,sound:375.7,sigma:8.431,wc:15590,rhoW:998.5,aw:0.9709,hydT:13.03},
  {calT:30,calP:150,calSal:3.5,calInhWt:30,bo:1.378,rhoG:226.1,zG:0.6274,muG:0.02745,co:1.837,cp:2.23,sound:417.3,sigma:2.774,wc:306.6,rhoW:1055,aw:0.8564,hydT:8.51},
  {calT:15,calP:80,calSal:0,calInhWt:0,bo:1.238,rhoG:132.9,zG:0.6014,muG:0.01657,co:1.38,cp:2.111,sound:309.5,sigma:7.163,wc:229.8,rhoW:1006,aw:1.001,hydT:16.14},
];
export default {
  id: 'pvt', num: 1, title: 'Fluid, PVT & Phase Behaviour', short: 'Fluid · PVT', icon: '🧪',
  tagline: 'Equation-of-state characterisation, phase envelope, laboratory experiments, properties, hydrate, wax and aqueous-phase thermodynamics of the case fluid.',
  description: 'The case fluid is characterised with a cubic equation of state (C7+ split, alternative critical-property correlations, tuning multipliers and interaction parameters) and flashed at constant PT, PH, PS and TV. The phase envelope is traced by continuation through the critical point, the standard PVT-cell experiments are simulated for comparison with a laboratory report, and every property the flow solvers need is tabulated. Hydrate equilibrium follows van der Waals–Platteeuw (Kihara cell potential by default, Munck constants as the alternative; structures I, II and H) and is published on the conservative side of its demonstrated uncertainty, with the water activity from activity-coefficient and electrolyte models; wax and asphaltene solid phases, the water-side Gibbs minimum and the miscibility pressure come from the same characterisation. GERG-2008, PC-SAFT, cubic-plus-association, Lee–Kesler and the Pedersen viscosity model run beside the cubic equation as cross-checks, and every constant set and reference data set carries its source.',
  guide: [
    'Define the fluid composition, rates, water cut, salinity and inhibitor on the case page; choose here whether to analyse that fluid or one of the built-in library fluids.',
    'Enter reservoir and reference conditions and, if a PVT report exists, its saturation pressures, CCE, differential-liberation, viscosity and hydrate points.',
    'On the model-setup tab select the equation of state, C7+ description, activity, electrolyte and hydrate models; leave the multipliers at 1 for an untuned model.',
    'Run, then read the phase envelope with the hydrate curve and wax appearance temperature against the operating points, and the comparison with the laboratory data.',
    'Regress the C7+ multipliers on the calibration tab, apply them, and re-run: the tuned property table and hydrate curve are then used by every other suite.',
  ],
  implemented: [
    'conservation of mass', 'component material-balance', 'gibbs phase-equilibrium criterion', 'equality of component fugacities', 'chemical-potential equilibrium', 'gibbs–duhem relation', 'first and second laws', 'fundamental thermodynamic relation', 'gibbs free-energy minimization',
    'phase-fraction/material-balance', 'rachford–rice', 'isothermal flash', 'isenthalpic flash', 'isentropic flash', 'multiphase stability',
    'peng–robinson', 'soave–redlich–kwong eos', 'redlich–kwong eos', 'van der waals eos', 'lee–kesler corresponding-states model', 'activity-coefficient models including nrtl, uniquac and wilson', "henry's-law", "raoult's-law", 'poynting correction',
    'electrolyte activity models', 'debye–hückel', 'davies equation', 'pitzer electrolyte model', 'corresponding-states viscosity', 'lohrenz–bray–clark', 'parachor', 'corresponding-states thermal-property', 'mixing rules', 'volume-translation',
    'van der waals–platteeuw', 'langmuir adsorption', 'fugacity-based water–hydrate equilibrium', 'electrolyte/inhibitor activity corrections', 'hydrate phase-stability',
    'eos + activity-coefficient model', 'eos + electrolyte model', 'eos + hydrate statistical thermodynamics', 'eos + empirical pvt regression', 'eos + corresponding-states transport', 'compositional flash + hydrate equilibrium', 'compositional flash + wax/asphaltene',
    // initial and boundary conditions
    'initial fluid pressure and temperature', 'overall hydrocarbon composition', 'oil composition', 'gas composition', 'aqueous composition', 'water content', 'salinity and ionic', 'initial oil', 'gas and water phase fractions', 'initial phase state', 'initial density and viscosity', 'initial enthalpy', 'hydrate inhibitor concentration', 'initial dissolved-gas',
    'pressure and temperature ranges', 'permitted phase states', 'composition entering through each fluid source', 'water and electrolyte composition entering', 'inhibitor composition and concentration', 'reservoir-fluid composition at well', 'separator or process conditions', 'specified thermodynamic constraints', 'saturation', 'hydrate-equilibrium or other phase-behaviour',
    // inputs and outputs
    'hydrocarbon and non-hydrocarbon composition', 'heavy fractions/pseudocomponents', 'molecular weights', 'critical properties', 'acentric factors', 'binary interaction parameters', 'laboratory pvt and phase-equilibrium measurements', 'equation-of-state/activity/electrolyte/hydrate model selections', 'hydrate inhibitor type and concentration',
    'consistent thermophysical-property package', 'phase equilibrium and phase envelope', 'bubble/dew points', 'phase fractions and compositions', 'density', 'viscosity', 'compressibility and z-factor', 'enthalpy', 'heat capacity', 'thermal conductivity', 'surface/interfacial tension', 'joule-thomson coefficient', 'water distribution/content', 'hydrate equilibrium/dissociation boundary', 'hydrate structure/stability region', 'subcooling/safety margin', 'inhibitor requirement and sensitivity',
    // calibration, verification, validation supported by the engine
    'constant composition expansion', 'constant volume depletion', 'differential liberation', 'separator tests', 'swelling tests', 'saturation-pressure measurements', 'bubble-point calibration', 'dew-point calibration', 'gas-oil-ratio calibration', 'oil-density calibration', 'oil-viscosity calibration',
    'heavy-end/pseudocomponent characterization', 'plus-fraction splitting/lumping', 'critical-property regression', 'acentric-factor adjustment', 'binary-interaction-parameter regression', 'eos parameter regression',
    'pure-component eos verification', 'binary-mixture verification', 'multicomponent flash verification', 'pt-flash', 'ph-flash', 'ps-flash', 'tv-flash', 'bubble-point solver', 'dew-point solver', 'critical-point calculation', 'phase-envelope tracing', 'material-balance closure', 'component mass-balance closure', 'energy-balance closure', 'thermodynamic consistency',
    'gibbs-energy minimum/stability', 'fugacity equality', 'gibbs–duhem consistency', 'maxwell-relation', 'eos derivative verification', 'analytic-vs-numerical derivative', 'unit/dimensional consistency', 'limiting/single-phase', 'regression reproducibility', 'cross-implementation benchmark', 'machine-precision conservation',
    // added with the residual-Helmholtz models, association, corresponding states, solid phases and the extended calibration
    'helmholtz free-energy formulation', 'benedict–webb–rubin family', 'cubic-plus-association eos', 'saft', 'pc-saft', 'gerg-type multiparameter eos', 'pedersen viscosity model', 'kihara potential', 'lennard–jones/devonshire-type intermolecular potentials',
    'hydrate structure-i, structure-ii and structure-h occupancy models', 'hydrate phase-stability/gibbs-energy minimization', 'compositional flash + wax/asphaltene solid-phase equilibrium', 'where solids are permitted', 'initial hydrate', 'asphaltene or other solid phase fractions',
    'slim-tube measurements', 'formation-volume-factor calibration', 'gas-density calibration', 'water-density calibration', 'gas-viscosity calibration', 'compressibility calibration', 'z-factor calibration', 'enthalpy/heat-capacity calibration', 'speed-of-sound calibration',
    'interfacial-tension calibration', 'surface-tension calibration', 'water-content calibration', 'volume-shift regression', 'hydrate-equilibrium calibration', 'hydrate dissociation p–t calibration', 'hydrate inhibitor-response calibration', 'brine/salinity correction calibration', 'water-activity calibration',
    'phase-envelope measurements', 'compressibility measurements', 'calorimetric data', 'interfacial-tension measurements', 'water-content measurements', 'high-pressure/high-temperature measurements', 'high-pressure/low-temperature deepwater-condition measurements', 'saline-water hydrate experiments',
    'pvt-cell experiments', 'bubble-point measurements', 'dew-point measurements', 'density measurements', 'viscosity measurements', 'hydrate equilibrium experiments', 'hydrate dissociation experiments', 'meg/methanol inhibitor', 'multicomponent reservoir-fluid',
  ],
  referenceOnly: [],
  equationsNote: "Cubic equations of state (Peng–Robinson 1978, Soave–Redlich–Kwong, Redlich–Kwong, van der Waals) with van der Waals one-fluid mixing rules and Péneloux-type volume translation give the two-phase vapour–liquid equilibrium and the published property table; the hydrocarbon flash carries no water. Residual-Helmholtz models are solved beside it on the phases of the cubic flash: GERG-2008 (the eleven kernel components only — phases with more than 1 mol % C7+ are not evaluated), PC-SAFT (non-associating, no polar terms; C7+ by interpolated n-alkane parameters with the segment diameter matched to the specific gravity), cubic-plus-association for water, methanol and MEG with the gas (water content, inhibitor loss, dissolved gas) and the Lee–Kesler modified Benedict–Webb–Rubin equation with Kay mixing for the gas; their derivatives are numerical. They are cross-checks and options, not the source of the table used by the other suites. Riazi–Daubert and Twu critical properties act as plus-fraction-average ratios to Kesler–Lee because the shared kernel accepts uniform multipliers only. The critical point is located on the traced saturation line, not by the Heidemann–Khalil criterion. Viscosity: Lohrenz–Bray–Clark in the table; Pedersen corresponding states (methane reference, GERG-2008 methane density) and Lee–Gonzalez–Eakin as alternatives. Hydrates: van der Waals–Platteeuw for structures I and II with Langmuir constants integrated from the Kihara potential in the Lennard-Jones–Devonshire cell (default; n-butane and CO₂ in structure II keep the Munck constants, H₂S parameters regressed in this work) or the Munck Langmuir constants, guest fugacities from the cubic equation, water activity from NRTL / UNIQUAC / Wilson (methanol, MEG), Pitzer / Davies / Debye–Hückel (NaCl equivalent) and Henry's-law gas solubility. On 309 measured dissociation points (at most six per publication) the default Kihara set has a bias of +0.7 K below 100 bar, +0.3 K at 100–960 bar and +1.1 K for mixtures; the Munck set is 1.2 K cold above 100 bar and 2.1 K cold above 300 bar. Above 300 bar the default set drifts cold as well (−1.4 K at 500 bar). Structure H is modelled for methane with a methylcyclohexane-type former only, with two constants regressed in this work. The water-side Gibbs minimum decides between aqueous liquid, ice and the hydrate structures and gives the equilibrium conversion; it is not a general multiphase flash of all components, and hydrate formation from a water-free gas is not modelled. The hydrate curve handed to the other suites is the model value made conservative by its demonstrated uncertainty: none of the inhibitor depression is withheld up to 30 wt %, 12 % from 50 wt % (the two laboratories that measured 50 wt % MEG differ by 4.6 K and the model lies between them), and 0.004 K per bar is added above 250 bara; above 50 wt % MEG or 40 wt % methanol, and below 240 K, the activity models are extrapolated and the result should be confirmed by measurement. Ethanol, DEG and TEG solutions are treated as ideal. The methanol–water and MEG–water activity parameters (NRTL, Wilson, UNIQUAC, each with a + b/T) were regressed in this work to vapour–liquid, limiting-activity-coefficient and ice-line measurements, never to hydrate data, which serve as validation only. Before this suite has run, the other suites use the kernel screening curve: a gas-gravity correlation fitted to this model (1 K warm on average, individual gases −2 to +6 K) with the depression −60 K·ln x_w. Multi-ion brines are reduced to NaCl of equal ionic strength. Wax: ideal-solution solid–liquid equilibrium (Won melting properties, Pedersen wax-forming fraction) on a carbon-number distribution. Asphaltenes: Flory–Huggins regular-solution solid–liquid equilibrium (Hirschberg) with the liquid solubility parameter from the cubic equation — strongly dependent on the asphaltene solubility parameter, which must be tuned to a measured onset — beside the de Boer and colloidal-instability screens. The miscibility pressure is a single-cell multiple-contact estimate of the slim-tube value. Validation: bundled data sets compare the engine with measurements (hydrate dissociation, water content, gas solubility — ThermoML and cited publications) and, code-to-code, with independent reference software (CoolProp: pure-fluid reference equations; GERG-2008 for mixtures — phase envelope, bubble points, flash split, density, speed of sound, Joule–Thomson coefficient, heat capacity, enthalpy differences, compressibility of defined natural-gas and light-oil mixtures). No openly licensed measurements were found for gas–oil interfacial tension of mixtures, slim-tube miscibility pressure, asphaltene onset or reservoir fluids with a plus fraction: for those the app only supports the comparison with the user's own laboratory values (table inputs with deviations), and the corresponding ticks mean no more than that.",

  inputs: [
    { group: 'Fluid and conditions', tab: 'inputs', help: 'The composition, rates, water cut, salinity and inhibitor come from the case fluid; the results tables show them back.', fields: [
      { key: 'fluidSource', label: 'Fluid analysed', type: 'select', value: 'case', options: [{ value: 'case', label: 'Case fluid (from the case page)' }, ...Object.entries(FLUID_LIBRARY).map(([value, l]) => ({ value, label: 'Library: ' + l.name }))], help: 'Library fluids demonstrate other fluid types without changing the case; only the case fluid feeds the other suites.' },
      { key: 'tRes', label: 'Reservoir temperature', unit: '°C', value: BASE.tRes, min: -20, max: 250, typical: [50, 150], help: 'Temperature of the saturation pressure and of the simulated laboratory experiments.' },
      { key: 'pRes', label: 'Reservoir pressure', unit: 'bara', value: BASE.pRes, min: 2, max: 1400, typical: [100, 700] },
      { key: 'pRef', label: 'Flowline reference pressure', unit: 'bara', value: 100, min: 2, max: 1000, typical: [30, 300], help: 'Pressure at the coldest point of the flowline (operating or shut-in) used for subcooling, inhibitor demand and live-oil wax.' },
      { key: 'tSeabed', label: 'Seabed (coldest wall) temperature', unit: '°C', value: BASE.tSeabed, min: -5, max: 60, typical: [2, 20] },
      { key: 'pArr', label: 'Arrival pressure', unit: 'bara', value: BASE.pOut, min: 1.05, max: 600, typical: [10, 80] },
      { key: 'tArr', label: 'Arrival temperature', unit: '°C', value: 30, min: -20, max: 150, typical: [10, 60], help: 'Used for the water content of the gas and the arrival flash.' },
      { key: 'waterDepth', label: 'Water depth (riser height)', unit: 'm', value: BASE.waterDepth, min: 0, max: 4000, typical: [100, 2500], help: 'Liquid head that limits depressurisation from the host.' },
      { key: 'margin', label: 'Hydrate safety margin', unit: '°C', value: 3, min: 0, max: 15, typical: [2, 5], help: 'The inhibitor requirement keeps the hydrate temperature this far below the seabed temperature.' },
    ] },
    { group: 'Laboratory PVT report', tab: 'inputs', help: 'Measured points are compared with the model and reported as deviations; clear a table to skip it. The defaults are a synthetic report of the reference fluid.', fields: [
      { key: 'useLab', label: 'Compare with laboratory data', type: 'bool', value: true },
      { key: 'labPsat', label: 'Saturation pressures', type: 'table', columns: [{ key: 't', label: 'Temperature', unit: '°C' }, { key: 'p', label: 'Saturation pressure', unit: 'bara' }], value: LAB.psat, showIf: (v) => v.useLab !== false },
      { key: 'labCCE', label: 'Constant composition expansion at reservoir temperature', type: 'table', columns: [{ key: 'p', label: 'Pressure', unit: 'bara' }, { key: 'v', label: 'Relative volume V/Vsat', unit: '–' }], value: LAB.cce, showIf: (v) => v.useLab !== false },
      { key: 'labDLE', label: 'Differential liberation at reservoir temperature', type: 'table', columns: [{ key: 'p', label: 'Pressure', unit: 'bara' }, { key: 'rs', label: 'Rs', unit: 'Sm³/Sm³' }, { key: 'bo', label: 'Bo', unit: 'm³/Sm³' }, { key: 'rho', label: 'Oil density', unit: 'kg/m³' }], value: LAB.dle, showIf: (v) => v.useLab !== false },
      { key: 'labVisc', label: 'Oil viscosity at reservoir temperature', type: 'table', columns: [{ key: 'p', label: 'Pressure', unit: 'bara' }, { key: 'mu', label: 'Viscosity', unit: 'mPa·s' }], value: LAB.visc, showIf: (v) => v.useLab !== false },
      { key: 'labHyd', label: 'Hydrate dissociation points (fresh water)', type: 'table', columns: [{ key: 'p', label: 'Pressure', unit: 'bara' }, { key: 't', label: 'Temperature', unit: '°C' }], value: [{ p: 20, t: 8.4 }, { p: 50, t: 13.6 }, { p: 100, t: 17.9 }], showIf: (v) => v.useLab !== false },
      { key: 'labProps', label: 'Other measurements (quantity codes: rhoG, zG, muG, cpG, cpO, soundG, soundO, co, sigma, wc, rhoW, psat, hydT, jt)', type: 'table', columns: [{ key: 'q', label: 'Quantity code', type: 'text' }, { key: 'p', label: 'Pressure', unit: 'bara' }, { key: 't', label: 'Temperature', unit: '°C' }, { key: 'value', label: 'Measured value' }], value: [], showIf: (v) => v.useLab !== false, help: 'Comparison path for gas density (kg/m³), Z, gas viscosity (mPa·s), heat capacities (kJ/kg/K), speed of sound (m/s), oil compressibility (1/bar), gas–oil interfacial tension (mN/m), water content of gas (mg/Sm³), aqueous density (kg/m³), saturation pressure (bara; the pressure column is ignored), hydrate temperature (°C) and Joule–Thomson coefficient (K/bar).' },
      { key: 'labMmp', label: 'Slim-tube minimum miscibility pressure (0 = not measured)', unit: 'bara', value: 0, min: 0, max: 1500, showIf: (v) => v.useLab !== false },
      { key: 'labAop', label: 'Upper asphaltene onset pressure at reservoir temperature (0 = not measured)', unit: 'bara', value: 0, min: 0, max: 1500, showIf: (v) => v.useLab !== false },
      { key: 'labWat', label: 'Wax appearance temperature (0 = not measured)', unit: '°C', value: 0, min: -30, max: 90, showIf: (v) => v.useLab !== false },
    ] },
    { group: 'Separator train, swelling gas and SARA', tab: 'inputs', fields: [
      { key: 'sepStages', label: 'Separator stages (a stock tank at 1.01325 bara, 15 °C is added)', type: 'table', columns: [{ key: 'p', label: 'Pressure', unit: 'bara' }, { key: 't', label: 'Temperature', unit: '°C' }], value: [{ p: BASE.separatorP, t: 40 }, { p: 5, t: 30 }] },
      { key: 'swellGas', label: 'Injection gas for the swelling test', type: 'select', value: 'lean', options: [{ value: 'none', label: 'No swelling test' }, ...Object.entries(INJECTION_GASES).map(([value, g]) => ({ value, label: g.name }))] },
      { key: 'saraSat', label: 'Saturates', unit: 'wt %', value: 44, min: 0, max: 100, help: 'SARA analysis of the stock-tank oil for the colloidal instability index.' },
      { key: 'saraAro', label: 'Aromatics', unit: 'wt %', value: 33, min: 0, max: 100 },
      { key: 'saraRes', label: 'Resins', unit: 'wt %', value: 20.5, min: 0, max: 100 },
      { key: 'saraAsp', label: 'Asphaltenes', unit: 'wt %', value: 2.5, min: 0, max: 100 },
    ] },
    { group: 'Initial state and permitted solid phases', tab: 'inputs', help: 'Which solid phases the equilibrium may form, and the solids already present at the start (for example after a shut-in). The equilibrium at the cold reference point tells whether they grow or dissolve.', fields: [
      { key: 'solidsAllowed', label: 'Solid phases permitted', type: 'select', value: 'all', options: opt({ all: 'Hydrate, ice, wax and asphaltene', hydrate: 'Hydrate and ice only', waxAsph: 'Wax and asphaltene only', none: 'None (fluid phases only)' }) },
      { key: 'hyd0', label: 'Initial hydrate (share of the water inventory)', unit: '%', value: 0, min: 0, max: 100 },
      { key: 'solid0', label: 'Initial wax solids in the oil', unit: 'wt %', value: 0, min: 0, max: 60 },
      { key: 'shFormer', label: 'Structure-H former present in the liquid (methylcyclohexane-type naphthene)', type: 'bool', value: false, help: 'Adds structure H (methane in the small and medium cavities, the heavy former in the large one) to the hydrate structures considered.' },
      { key: 'aspMW', label: 'Asphaltene molar mass', unit: 'g/mol', value: 750, min: 300, max: 5000, typical: [500, 2500] },
      { key: 'aspRho', label: 'Asphaltene density', unit: 'kg/m³', value: 1100, min: 900, max: 1300 },
      { key: 'aspDelta', label: 'Asphaltene solubility parameter at 25 °C', unit: 'MPa^0.5', value: 21, min: 15, max: 30, typical: [19, 24], help: 'Tune to a measured onset pressure; the amount comes from the SARA asphaltene content.' },
      { key: 'aspDeltaT', label: 'Relative decrease of that parameter per K', unit: '1/K', value: 1.07e-3, min: 0, max: 3e-3 },
    ] },
    { group: 'Equation of state and characterisation', tab: 'setup', fields: [
      { key: 'eosSel', label: 'Equation of state', type: 'select', value: 'case', options: [{ value: 'case', label: 'As set in the case' }, ...Object.entries(EOS).map(([value, e]) => ({ value, label: e.name }))] },
      { key: 'nPseudoSel', label: 'C7+ pseudo-components', type: 'select', value: 'case', options: opt({ case: 'As set in the case', 1: '1 (lumped)', 2: '2 (gamma quadrature)', 3: '3 (gamma quadrature)' }) },
      { key: 'critCorr', label: 'Critical-property correlation', type: 'select', value: 'KL', options: opt({ KL: 'Kesler–Lee', RD: 'Riazi–Daubert', Twu: 'Twu' }) },
      { key: 'shift', label: 'Volume translation of liquid density', type: 'bool', value: true },
      { key: 'tcMult', label: 'C7+ critical-temperature multiplier', unit: '–', value: 1, min: 0.7, max: 1.3, typical: [0.95, 1.05] },
      { key: 'pcMult', label: 'C7+ critical-pressure multiplier', unit: '–', value: 1, min: 0.6, max: 1.5, typical: [0.9, 1.1] },
      { key: 'wMult', label: 'C7+ acentric-factor multiplier', unit: '–', value: 1, min: 0.5, max: 1.6, typical: [0.9, 1.1] },
      { key: 'kijScale', label: 'Binary-interaction scale', unit: '–', value: 1, min: 0, max: 3, typical: [0.5, 1.5], help: 'Scales every default kij (including methane–C7+).' },
      { key: 'vcMult', label: 'C7+ critical-volume multiplier (viscosity)', unit: '–', value: 1, min: 0.5, max: 2, typical: [0.9, 1.2], help: 'Tunes the Lohrenz–Bray–Clark liquid viscosity.' },
      { key: 'kijTable', label: 'Binary interaction overrides (ids: N2, CO2, H2S, C1 … C6, C7+)', type: 'table', columns: [{ key: 'a', label: 'Component i', type: 'text' }, { key: 'b', label: 'Component j', type: 'text' }, { key: 'kij', label: 'kij', unit: '–' }], value: [] },
      { key: 'altEos', label: 'Alternative equations of state and viscosity model', type: 'select', value: 'key', options: opt({ key: 'Cross-check at the key conditions: PC-SAFT, GERG-2008, Lee–Kesler, CPA, Pedersen', all: 'As above plus the Z-factor curves and the second hydrate constants set', none: 'Cubic model only (fastest)' }), help: 'The alternative models are evaluated on the phases of the cubic flash at the key conditions; the published property table stays the cubic model.' },
      { key: 'mmpCalc', label: 'Estimate the miscibility pressure of the swelling gas (adds a few tenths of a second)', type: 'bool', value: false },
      { key: 'scnModel', label: 'Carbon-number distribution (wax, plot)', type: 'select', value: 'exp', options: opt({ exp: 'Exponential (Pedersen)', gamma: 'Gamma (Whitson)' }) },
      { key: 'gammaAlpha', label: 'Gamma shape α', unit: '–', value: 1, min: 0.5, max: 3, showIf: (v) => v.scnModel === 'gamma' },
    ] },
    { group: 'Aqueous phase, hydrate and wax models', tab: 'setup', fields: [
      { key: 'hydModel', label: 'Hydrate model', type: 'select', value: 'kihara', options: opt({ kihara: 'van der Waals–Platteeuw, Kihara cell potential (default: best validated, errs warm)', vdwp: 'van der Waals–Platteeuw, Munck Langmuir constants (1–2.4 K cold above 100 bara)', corr: 'Gas-gravity screening correlation + Nielsen–Bucklin' }), help: 'Both constants sets are compared in the results. On 309 measured dissociation points from 1 to 960 bar the Kihara set has a bias of +0.3 to +1.1 K (warm, i.e. conservative); the Munck set runs 1.2 K cold above 100 bar and 2.1 K cold above 300 bar.' },
      { key: 'actModel', label: 'Activity-coefficient model (inhibitor–water)', type: 'select', value: 'NRTL', options: opt({ NRTL: 'NRTL', UNIQUAC: 'UNIQUAC', Wilson: 'Wilson', ideal: "Ideal solution (Raoult's law)" }) },
      { key: 'elecModel', label: 'Electrolyte model', type: 'select', value: 'pitzer', options: opt({ pitzer: 'Pitzer', davies: 'Davies', edh: 'Extended Debye–Hückel', dh: 'Debye–Hückel limiting law' }) },
      { key: 'inhDesign', label: 'Inhibitor for the requirement (when the case has none)', type: 'select', value: 'MEG', options: Object.entries(INHIBITORS).filter(([k]) => k !== 'none').map(([value, i]) => ({ value, label: i.name })) },
      { key: 'aqSource', label: 'Additional what-if aqueous phase', type: 'select', value: 'case', options: opt({ case: 'None (case aqueous phase only)', custom: 'Custom salinity / ions / inhibitor' }) },
      { key: 'salinityIn', label: 'Custom salinity (NaCl equivalent)', unit: 'wt %', value: 3.5, min: 0, max: 26, showIf: (v) => v.aqSource === 'custom' },
      { key: 'aqIons', label: 'Custom brine analysis (replaces the salinity when filled; ions: Na, K, Ca, Mg, Ba, Sr, Fe, Cl, SO4, HCO3, Br)', type: 'table', columns: [{ key: 'ion', label: 'Ion', type: 'text' }, { key: 'mgL', label: 'Concentration', unit: 'mg/L' }], value: [], showIf: (v) => v.aqSource === 'custom' },
      { key: 'inhIn', label: 'Custom inhibitor', type: 'select', value: 'MEG', options: Object.entries(INHIBITORS).map(([value, i]) => ({ value, label: i.name })), showIf: (v) => v.aqSource === 'custom' },
      { key: 'inhWtIn', label: 'Custom inhibitor concentration', unit: 'wt %', value: 30, min: 0, max: 90, showIf: (v) => v.aqSource === 'custom' },
      { key: 'waxDetect', label: 'Wax detection limit for the appearance temperature', unit: 'wt %', value: 0.02, min: 0.001, max: 1 },
      { key: 'waxHfMult', label: 'Wax enthalpy-of-fusion multiplier', unit: '–', value: 1, min: 0.5, max: 1.5, help: 'Tune to a measured wax appearance temperature.' },
    ] },
    { group: 'Calibration multipliers (estimated on the calibration tab)', tab: 'setup', help: 'Each multiplier moves one measured quantity; leave at 1 for the untuned model. They act on the run, on the published property table (volume shift, parachor, heat capacity, gas viscosity) and on the hydrate curve.', fields: [
      { key: 'shiftMult', label: 'Volume-shift multiplier (densities, Z, formation-volume factor)', unit: '–', value: 1, min: 0, max: 3, typical: [0.7, 1.3] },
      { key: 'parMult', label: 'Parachor multiplier (interfacial and surface tension)', unit: '–', value: 1, min: 0.5, max: 1.5, typical: [0.9, 1.1] },
      { key: 'cpMult', label: 'Ideal-gas heat-capacity multiplier (enthalpy, heat capacity, speed of sound)', unit: '–', value: 1, min: 0.7, max: 1.3, typical: [0.95, 1.05] },
      { key: 'muGMult', label: 'Gas-viscosity multiplier', unit: '–', value: 1, min: 0.5, max: 2, typical: [0.9, 1.1] },
      { key: 'rhoWMult', label: 'Aqueous-density multiplier', unit: '–', value: 1, min: 0.9, max: 1.1, typical: [0.99, 1.01] },
      { key: 'wcMult', label: 'Water-content multiplier', unit: '–', value: 1, min: 0.3, max: 3, typical: [0.8, 1.2] },
      { key: 'hydMult', label: 'Hydrate reference chemical-potential multiplier (shifts the dissociation curve)', unit: '–', value: 1, min: 0.8, max: 1.2, typical: [0.97, 1.03] },
      { key: 'actMult', label: 'Inhibitor activity multiplier (ln γ of water: inhibitor response)', unit: '–', value: 1, min: 0, max: 3, typical: [0.7, 1.3] },
      { key: 'saltMult', label: 'Salt activity multiplier (ln a_w of the brine: salinity correction)', unit: '–', value: 1, min: 0, max: 3, typical: [0.8, 1.2] },
    ] },
    { group: 'Property-table range and resolution', tab: 'mesh', help: 'The table is interpolated by the flow, solids and operations suites.', fields: [
      { key: 'nP', label: 'Pressure points (logarithmic)', value: 22, min: 6, max: 80 },
      { key: 'nT', label: 'Temperature points', value: 17, min: 5, max: 70 },
      { key: 'tblPmin', label: 'Lowest pressure', unit: 'bara', value: 1, min: 0.5, max: 50 },
      { key: 'tblPmax', label: 'Highest pressure', unit: 'bara', value: 600, min: 100, max: 1500 },
      { key: 'tblTmin', label: 'Lowest temperature', unit: '°C', value: -30, min: -60, max: 20 },
      { key: 'tblTmax', label: 'Highest temperature', unit: '°C', value: 170, min: 60, max: 300 },
      { key: 'nEnv', label: 'Phase-envelope resolution (points per 40 standard steps)', value: 40, min: 10, max: 300 },
      { key: 'nHyd', label: 'Hydrate-curve points (5–500 bara)', value: 24, min: 20, max: 120 },
      { key: 'nIso', label: 'Points on the property curves', value: 24, min: 8, max: 120 },
    ] },
  ],
  presets: [
    { name: 'Reference light oil with laboratory report', values: { fluidSource: 'case' } },
    { name: 'Black oil, low GOR', values: { fluidSource: 'blackOil', tRes: 80, pRes: 250, pRef: 80, useLab: false, swellGas: 'CO2' } },
    { name: 'Volatile oil, near-critical', values: { fluidSource: 'volatileOil', tRes: 120, pRes: 420, pRef: 150, useLab: false, swellGas: 'lean' } },
    { name: 'Gas condensate with constant-volume depletion', values: { fluidSource: 'gasCondensate', tRes: 110, pRes: 420, pRef: 180, tArr: 20, useLab: false, swellGas: 'none', sepStages: [{ p: 60, t: 25 }, { p: 10, t: 20 }] } },
    { name: 'Lean gas tie-back with MEG injection', values: { fluidSource: 'leanGas', altEos: 'all', tRes: 95, pRes: 280, pRef: 150, pArr: 70, tArr: 8, useLab: false, swellGas: 'none', aqSource: 'custom', salinityIn: 0.5, inhIn: 'MEG', inhWtIn: 45, inhDesign: 'MEG' } },
    { name: 'Sour CO₂-rich gas, SRK and methanol', values: { fluidSource: 'sourGas', eosSel: 'SRK', tRes: 100, pRes: 300, pRef: 120, pArr: 60, tArr: 10, useLab: false, swellGas: 'none', inhDesign: 'MeOH', aqSource: 'custom', salinityIn: 2, inhIn: 'MeOH', inhWtIn: 25 } },
    { name: 'Reference oil tuned to the laboratory report', values: { fluidSource: 'case', tcMult: 1.025, pcMult: 0.97, kijScale: 1.15, vcMult: 1.06 } },
    { name: 'Munck hydrate constants, structure-H former, miscibility and solids inventory', values: { fluidSource: 'case', altEos: 'all', hydModel: 'vdwp', shFormer: true, mmpCalc: true, swellGas: 'rich', hyd0: 20, solid0: 0.5, labProps: [{ q: 'sigma', p: 100, t: 60, value: 6.1 }, { q: 'wc', p: 25, t: 30, value: 1500 }] } },
    { name: 'Heavy-end sensitivity: SRK, Twu, two pseudo-components', values: { fluidSource: 'case', eosSel: 'SRK', critCorr: 'Twu', nPseudoSel: '2', useLab: false } },
  ],
  pull: ({ fluid, outputs } = {}) => [
    { key: 'tRes', value: fluid?.Tres, from: 'Case fluid: reservoir temperature' }, { key: 'pRes', value: fluid?.Pres, from: 'Case fluid: reservoir pressure' }, { key: 'pArr', value: fluid?.Pout, from: 'Case fluid: arrival pressure' },
    { key: 'tSeabed', value: outputs?.net?.tSeabed, from: 'Network suite: seabed temperature' }, { key: 'waterDepth', value: outputs?.net?.waterDepth, from: 'Network suite: water depth' },
    { key: 'pRef', value: outputs?.flow?.pIn, from: 'Flow suite: flowline inlet pressure' }, { key: 'tArr', value: outputs?.flow?.tOut, from: 'Flow suite: arrival temperature' },
  ].filter((it) => Number.isFinite(it.value)),
  site: (site) => [
    { key: 'tSeabed', value: site?.data?.seabedTemp, from: 'Seabed temperature at the site' }, { key: 'waterDepth', value: site?.data?.depth, from: 'Water depth at the site' },
  ].filter((it) => Number.isFinite(it.value)),

  run: (v, ctx) => engine(v, ctx),

  mesh: [
    { name: 'Property table and hydrate curve resolution', keys: ['nP', 'nT', 'nHyd'], min: 8, note: 'Bilinear interpolation in ln P and T: the metrics are values interpolated at an off-grid point, as the flow solvers see them.', metrics: [
      { label: 'Interpolated oil density at 77 bara, 43 °C', unit: 'kg/m³', get: (r) => lookup(r.outputs.table, 77, 43).rhoO },
      { label: 'Interpolated gas mass fraction at 77 bara, 43 °C', unit: '–', get: (r) => lookup(r.outputs.table, 77, 43).wG },
      { label: 'Interpolated oil viscosity at 77 bara, 43 °C', unit: 'mPa·s', get: (r) => lookup(r.outputs.table, 77, 43).muO * 1e3 },
      { label: 'Hydrate temperature interpolated at 77 bara', unit: '°C', get: (r) => interp1(r.outputs.hydrateCurve.P.map(Math.log), r.outputs.hydrateCurve.T, Math.log(77)) },
    ] },
    { name: 'Phase-envelope step size', keys: ['nEnv'], min: 12, note: 'Continuation step length of the envelope trace; cricondenbar and cricondentherm are refined by parabolic interpolation.', metrics: [
      { label: 'Cricondenbar', unit: 'bara', get: (r) => r.outputs.cricondenbar ?? 0 }, { label: 'Cricondentherm', unit: '°C', get: (r) => r.outputs.cricondentherm ?? 0 }, { label: 'Critical temperature', unit: '°C', get: (r) => r.outputs.critical?.T ?? 0 }, { label: 'Critical pressure', unit: 'bara', get: (r) => r.outputs.critical?.P ?? 0 },
    ] },
  ],

  calibration: {
    note: 'Regress the C7+ critical-property multipliers, the interaction scale and the viscosity critical-volume multiplier to saturation pressure, liquid density, solution gas–oil ratio and viscosity, and the single-purpose multipliers to their own measurement: volume shift → densities, Z and formation-volume factor; parachor → interfacial tension; ideal-gas heat capacity → heat capacity and speed of sound; gas viscosity; aqueous density; water content; hydrate reference chemical potential → fresh-water hydrate points; inhibitor activity → inhibited hydrate points and water activity; salt activity → brine hydrate points. A row only needs the columns that were measured; give the salinity and inhibitor of the aqueous phase with hydrate, water-activity, water-content and aqueous-density points. The model uses the fluid of the most recent run (the reference fluid before any run). Fit one group of parameters at a time (tick only those the data can identify).',
    params: [
      { key: 'tcMult', label: 'C7+ Tc multiplier', lo: 0.9, hi: 1.1 }, { key: 'pcMult', label: 'C7+ Pc multiplier', lo: 0.8, hi: 1.2 }, { key: 'wMult', label: 'C7+ acentric-factor multiplier', lo: 0.8, hi: 1.2 },
      { key: 'kijScale', label: 'Binary-interaction scale', lo: 0, hi: 2.5 }, { key: 'vcMult', label: 'C7+ Vc multiplier (viscosity)', lo: 0.8, hi: 1.4 },
      { key: 'shiftMult', label: 'Volume-shift multiplier', lo: 0.3, hi: 2 }, { key: 'parMult', label: 'Parachor multiplier', lo: 0.7, hi: 1.3 }, { key: 'cpMult', label: 'Ideal-gas heat-capacity multiplier', lo: 0.85, hi: 1.15 },
      { key: 'muGMult', label: 'Gas-viscosity multiplier', lo: 0.7, hi: 1.4 }, { key: 'rhoWMult', label: 'Aqueous-density multiplier', lo: 0.95, hi: 1.05 }, { key: 'wcMult', label: 'Water-content multiplier', lo: 0.5, hi: 2 },
      { key: 'hydMult', label: 'Hydrate reference chemical-potential multiplier', lo: 0.9, hi: 1.1 }, { key: 'actMult', label: 'Inhibitor activity multiplier', lo: 0.3, hi: 2 }, { key: 'saltMult', label: 'Salt activity multiplier', lo: 0.5, hi: 1.6 },
    ],
    columns: [{ key: 'calT', label: 'Temperature', unit: '°C' }, { key: 'calP', label: 'Pressure', unit: 'bara' }, { key: 'calSal', label: 'Salinity of the aqueous phase', unit: 'wt %' }, { key: 'calInhWt', label: 'Inhibitor in the aqueous phase', unit: 'wt %' },
      { key: 'psat', label: 'Saturation pressure at T', unit: 'bara' }, { key: 'rho', label: 'Liquid density at P, T', unit: 'kg/m³' }, { key: 'rs', label: 'Solution GOR at P, T', unit: 'Sm³/Sm³' }, { key: 'mu', label: 'Liquid viscosity at P, T', unit: 'mPa·s' },
      { key: 'bo', label: 'Formation-volume factor', unit: 'm³/Sm³' }, { key: 'rhoG', label: 'Gas density', unit: 'kg/m³' }, { key: 'zG', label: 'Gas Z-factor', unit: '–' }, { key: 'muG', label: 'Gas viscosity', unit: 'mPa·s' }, { key: 'co', label: 'Liquid compressibility', unit: '10⁻⁴/bar' },
      { key: 'cp', label: 'Liquid heat capacity', unit: 'kJ/kg/K' }, { key: 'sound', label: 'Speed of sound in the gas', unit: 'm/s' }, { key: 'sigma', label: 'Interfacial tension', unit: 'mN/m' }, { key: 'wc', label: 'Water content of gas', unit: 'mg/Sm³' }, { key: 'rhoW', label: 'Aqueous density', unit: 'kg/m³' },
      { key: 'aw', label: 'Water activity', unit: '–' }, { key: 'hydT', label: 'Hydrate temperature', unit: '°C' }, { key: 'mmp', label: 'Minimum miscibility pressure', unit: 'bara' }],
    targets: [{ key: 'psat', label: 'Saturation pressure', unit: 'bara' }, { key: 'rho', label: 'Liquid density', unit: 'kg/m³' }, { key: 'rs', label: 'Solution GOR', unit: 'Sm³/Sm³' }, { key: 'mu', label: 'Liquid viscosity', unit: 'mPa·s' },
      { key: 'bo', label: 'Formation-volume factor', unit: 'm³/Sm³' }, { key: 'rhoG', label: 'Gas density (stock-tank flash gas at P, T)', unit: 'kg/m³' }, { key: 'zG', label: 'Gas Z-factor', unit: '–' }, { key: 'muG', label: 'Gas viscosity', unit: 'mPa·s' }, { key: 'co', label: 'Liquid isothermal compressibility', unit: '10⁻⁴/bar' },
      { key: 'cp', label: 'Liquid heat capacity', unit: 'kJ/kg/K' }, { key: 'sound', label: 'Speed of sound in the gas', unit: 'm/s' }, { key: 'sigma', label: 'Gas–oil interfacial (surface) tension', unit: 'mN/m' }, { key: 'wc', label: 'Water content of gas', unit: 'mg/Sm³' }, { key: 'rhoW', label: 'Aqueous-phase density', unit: 'kg/m³' },
      { key: 'aw', label: 'Water activity of the aqueous phase', unit: '–' }, { key: 'hydT', label: 'Hydrate dissociation temperature at P', unit: '°C' }, { key: 'mmp', label: 'Slim-tube minimum miscibility pressure at T', unit: 'bara' }],
    model: calibrationModel,
    sample: [...CAL_SAMPLE, ...CAL_SAMPLE2],
    validationSample: [...CAL_VALID, ...CAL_VALID2],
  },

  verify: verifyChecks,
  validationData: VALIDATION,

  views: [{ id: 'flash', label: 'Flash calculator', tip: 'Interactive PT / PH / PS flash with phase compositions and properties', render: flashView }],
};
