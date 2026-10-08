// Suite 1 — Fluid, PVT & Phase Behaviour.
// Builds on the shared cubic-EOS kernel (core/thermo.js) and adds: alternative C7+ critical-property correlations and tuning,
// PH / PS / TV flashes, saturation points and a traced phase envelope with quality lines, the standard PVT-laboratory
// experiments (CCE, DLE, CVD, separator test, swelling), derived properties (speed of sound, compressibility, alternative
// viscosity / conductivity / Z-factor models, black-oil correlations), aqueous-phase activity models (NRTL, UNIQUAC, Wilson,
// Debye–Hückel family, Pitzer), Henry's-law gas solubility, the van der Waals–Platteeuw hydrate model, ideal-solution wax
// precipitation and asphaltene screening. Units at the interface: bara, °C, mol %; SI inside.
import { brent, clamp, linspace, logspace, interp1, isNum, solveLinear } from '../core/num.js';
import {
  R, P_STD, T_STD, VM_STD, MW_AIR, COMP_IDS, COMP_LABELS, INHIBITORS, EOS, DEFAULT_FLUID, makeFluid, eosPhase, rachfordRice, stability, flashPT,
  phaseProps, props, saturationP, stdFlash, streams, aqueous, waterContent, hydrateDepression, inhibitorFor, hydrateT0, pseudoProps, buildTable, lookup,
} from '../core/thermo.js';
import { psat as psatWater, density as rhoBrine } from '../core/props.js';
import { BASE } from '../data/basecase.js';

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
  const fl = flashPT(f, Pfloor, t);
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
  let a = null, b = null;
  for (let t = Tmax; t >= Tmin; t -= 25) { if (isTwo(f, Pbar, t)) { a = t; break; } b = t; }
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
function lkZ(c, Tr, Pr) { // vapour-like root of the Lee–Kesler modified BWR equation
  const B = c[0] - c[1] / Tr - c[2] / Tr ** 2 - c[3] / Tr ** 3, C = c[4] - c[5] / Tr + c[6] / Tr ** 3, D = c[8] + c[9] / Tr;
  const Z = (V) => 1 + B / V + C / V ** 2 + D / V ** 5 + (c[7] / (Tr ** 3 * V * V)) * (c[10] + c[11] / (V * V)) * Math.exp(-c[11] / (V * V));
  const g = (lv) => { const V = Math.exp(lv); return Z(V) - (Pr * V) / Tr; };
  // scan from the ideal-gas side towards small volumes for the first sign change
  let hi = Math.log((Tr / Pr) * 3), ghi = g(hi);
  for (let k = 0; k < 400; k++) { const lo = hi - 0.02, glo = g(lo); if (glo * ghi <= 0) { const lv = brent(g, lo, hi, 1e-13); return Z(Math.exp(lv)); } hi = lo; ghi = glo; }
  return null;
}
/** Lee–Kesler (1975) corresponding-states compressibility factor (vapour-like root) at reduced T, P and acentric factor w. */
export function leeKeslerZ(Tr, Pr, w = 0) {
  const z0 = lkZ(LK0, Tr, Pr), zr = lkZ(LKR, Tr, Pr);
  return z0 === null || zr === null ? null : z0 + (w / 0.3978) * (zr - z0);
}
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
const RCAL = 1.9872;
// Component 1 = inhibitor, 2 = water. Each model carries two parameter pairs: `hi` for vapour–liquid temperatures (≥ 350 K) and
// `lo` for the cold end (≤ 255 K); the energies are interpolated linearly in between. NRTL / Wilson energies in cal/mol, UNIQUAC in K.
// Methanol `hi`: NRTL and Wilson from the Gmehling–Onken VLE compilation as tabulated by Smith, Van Ness & Abbott; UNIQUAC matched to
// that NRTL set. `lo` sets (methanol and MEG): regressed in this work to the ice-point depression of the aqueous solutions
// (methanol 10–40 wt %: −6.5, −15.0, −25.9, −38.6 °C; MEG 10–50 wt %: −3.4, −7.9, −14.0, −22.3, −33.8 °C); MEG uses them at all temperatures.
export const ACTIVITY_PARAMS = Object.freeze({
  MeOH: { name: 'Methanol', V: 40.73, r: 1.4311, q: 1.432, alpha: 0.2994, nrtl: { lo: [-463.53, 284.31], hi: [-253.88, 845.21] }, wilson: { lo: [44.66, -63.04], hi: [107.38, 469.55] }, uniquac: { lo: [-224.88, 172.67], hi: [-199.18, 319.94] } },
  MEG: { name: 'Mono-ethylene glycol', V: 55.92, r: 2.4088, q: 2.248, alpha: 0.3, nrtl: { lo: [-639.84, 387.04], hi: [-639.84, 387.04] }, wilson: { lo: [-582.74, 260.72], hi: [-582.74, 260.72] }, uniquac: { lo: [-147.71, -44.66], hi: [-147.71, -44.66] } },
});
const V_W = 18.07, R_W = 0.92, Q_W = 1.4, T_LO = 255, T_HI = 350;
const actPair = (set, TK) => { const u = clamp((TK - T_LO) / (T_HI - T_LO), 0, 1); return [set.lo[0] + u * (set.hi[0] - set.lo[0]), set.lo[1] + u * (set.hi[1] - set.lo[1])]; };
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
    const [a12, a21] = actPair(p.wilson, TK), L12 = (V_W / p.V) * Math.exp(-a12 / (RCAL * TK)), L21 = (p.V / V_W) * Math.exp(-a21 / (RCAL * TK)), br = L12 / (x1 + L12 * x2) - L21 / (x2 + L21 * x1);
    l1 = -Math.log(x1 + L12 * x2) + x2 * br; l2 = -Math.log(x2 + L21 * x1) - x1 * br;
  } else if (model === 'UNIQUAC') {
    const [a12, a21] = actPair(p.uniquac, TK), z = 10, r1 = p.r, q1 = p.q, t12 = Math.exp(-a12 / TK), t21 = Math.exp(-a21 / TK);
    const phi1 = (x1 * r1) / (x1 * r1 + x2 * R_W), phi2 = 1 - phi1, th1 = (x1 * q1) / (x1 * q1 + x2 * Q_W), th2 = 1 - th1, la = (z / 2) * (r1 - q1) - (r1 - 1), lb = (z / 2) * (R_W - Q_W) - (R_W - 1);
    l1 = Math.log(phi1 / x1) + (z / 2) * q1 * Math.log(th1 / phi1) + phi2 * (la - (r1 / R_W) * lb) - q1 * Math.log(th1 + th2 * t21) + th2 * q1 * (t21 / (th1 + th2 * t21) - t12 / (th2 + th1 * t12));
    l2 = Math.log(phi2 / x2) + (z / 2) * Q_W * Math.log(th2 / phi2) + phi1 * (lb - (R_W / r1) * la) - Q_W * Math.log(th2 + th1 * t12) + th1 * Q_W * (t12 / (th2 + th1 * t12) - t21 / (th1 + th2 * t21));
  } else {
    const [b12, b21] = actPair(p.nrtl, TK), t12 = b12 / (RCAL * TK), t21 = b21 / (RCAL * TK), G12 = Math.exp(-p.alpha * t12), G21 = Math.exp(-p.alpha * t21);
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
const pitzerNaCl = (TK) => { const d = TK - 298.15; return { b0: 0.0765 + 7.159e-4 * d, b1: 0.2664 + 7.005e-4 * d, c: 0.00127 - 1.054e-4 * d }; };
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
export function waterActivity({ S = 0, inhId = 'none', inhWt = 0 }, TK, { act = 'NRTL', elec = 'pitzer' } = {}) {
  const s = clamp(S, 0, 260) / 1000, m = s > 0 ? (s / (MW_NACL * 1e-3)) / (1 - s) : 0, awSalt = waterActivityNaCl(elec, m, TK);
  const w = inhId && inhId !== 'none' ? clamp(inhWt, 0, 95) / 100 : 0, mwI = (INHIBITORS[inhId] || INHIBITORS.none).MW;
  // salt-free mole fraction of inhibitor: w kg inhibitor with (1 - w)(1 - s) kg water
  const nI = w / mwI, nW = ((1 - w) * (1 - s)) / MW_W, x1 = nI + nW > 0 ? nI / (nI + nW) : 0, a = activityBinary(act, inhId, x1, TK);
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
const T0H = 273.15;
export const HYDRATE_STRUCTURES = Object.freeze({
  sI: { nu: [1 / 23, 3 / 23], dmu0: 1264, dhL: -4858, dhI: 1151, dvL: 4.6e-6, dvI: 3.0e-6, waters: 23 },
  sII: { nu: [2 / 17, 1 / 17], dmu0: 883, dhL: -5201, dhI: 808, dvL: 5.0e-6, dvI: 3.4e-6, waters: 17 },
});
const DCP_L = -39.16; // J/mol/K, empty lattice minus liquid water
export const LANGMUIR = Object.freeze({
  C1: { sI: [[0.7228e-3, 3187], [23.35e-3, 2653]], sII: [[0.2207e-3, 3453], [100e-3, 1916]] },
  C2: { sI: [[0, 0], [3.039e-3, 3861]], sII: [[0, 0], [240e-3, 2967]] },
  C3: { sI: [[0, 0], [0, 0]], sII: [[0, 0], [5.455e-3, 4638]] },
  iC4: { sI: [[0, 0], [0, 0]], sII: [[0, 0], [189.3e-3, 3800]] },
  nC4: { sI: [[0, 0], [0, 0]], sII: [[0, 0], [30.51e-3, 3699]] },
  N2: { sI: [[1.617e-3, 2905], [6.078e-3, 2431]], sII: [[0.1742e-3, 3082], [18e-3, 1728]] },
  CO2: { sI: [[0.00588e-3, 5410], [3.36e-3, 3202]], sII: [[0.0846e-3, 3602], [846e-3, 2030]] },
  H2S: { sI: [[0.025e-3, 4568], [16.34e-3, 3737]], sII: [[0.0298e-3, 4878], [87.2e-3, 2633]] },
});
/** Langmuir constant (1/bar) of guest `id` in cavity m (0 small, 1 large) of structure s at T (K). */
export const langmuirC = (id, s, m, TK, table = LANGMUIR) => { const c = table[id]?.[s]?.[m]; return c && c[0] > 0 ? ((c[0] / TK) * Math.exp(c[1] / TK)) / ATM : 0; };
/**
 * Chemical-potential differences of water (divided by RT) at T (K), P (bara) for hydrate structure s:
 * hydrate side Σ ν ln(1 + Σ C f) and water side (liquid with activity aw, or ice when that is the stable phase).
 * fug: { id: fugacity in bar }. Returns { dmuH, dmuW, ice, theta: [{ id: occupancy } small, large], drive = dmuH − dmuW (> 0: hydrate stable) }.
 */
export function hydrateState(s, TK, Pbar, fug, aw = 1, table = LANGMUIR) {
  const S = HYDRATE_STRUCTURES[s], theta = [{}, {}]; let dmuH = 0;
  for (let m = 0; m < 2; m++) {
    let sm = 0; const cf = {};
    for (const id in fug) { const v = langmuirC(id, s, m, TK, table) * fug[id]; if (v > 0) { cf[id] = v; sm += v; } }
    for (const id in cf) theta[m][id] = cf[id] / (1 + sm);
    dmuH += S.nu[m] * Math.log(1 + sm);
  }
  const P = Pbar * 1e5, Tm = 0.5 * (TK + T0H), base = S.dmu0 / (R * T0H), inv = 1 / T0H - 1 / TK;
  const liq = base - (((S.dhL - DCP_L * T0H) / R) * inv + (DCP_L / R) * Math.log(TK / T0H)) + (S.dvL * P) / (R * Tm) - Math.log(Math.max(aw, 1e-6));
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
export function hydrateTofP(f, Pbar, aw = 1, { Tmin = -45, Tmax = 50, table = LANGMUIR, molality = 0, solubility = true, guess = 10, passes = 10 } = {}) {
  const awAt = typeof aw === 'function' ? aw : () => aw;
  const awTot = (fug, t) => awAt(t + KEL) * (1 - (solubility ? gasSolubility(fug, Pbar, t, molality).total : 0));
  const solve = (fug) => {
    const res = {};
    for (const s of ['sI', 'sII']) {
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
    if (res.sI === null && res.sII === null) return null;
    const Tn = Math.max(res.sI ?? -1e9, res.sII ?? -1e9), done = Math.abs(Tn - T) < 5e-3;
    T = it > 5 ? 0.5 * (T + Tn) : Tn;
    if (done || it === passes - 1) break;
    fug = formerFugacities(f, Pbar, T);
  }
  const structure = (res.sII ?? -1e9) > (res.sI ?? -1e9) ? 'sII' : 'sI', st = hydrateState(structure, T + KEL, Pbar, fug, awTot(fug, T), table);
  return { T, structure, TsI: res.sI, TsII: res.sII, ice: st.ice, occupancy: occupancy(structure, st), xGas: solubility ? gasSolubility(fug, Pbar, T, molality).total : 0, fug };
}
/** Lowest hydrate dissociation pressure (bara) at T (°C); null when no hydrate forms below Pmax. Returns { P, structure, occupancy }. */
export function hydratePofT(f, Tc, aw = 1, { Pmin = 0.2, Pmax = 1500, table = LANGMUIR, molality = 0, solubility = true } = {}) {
  const a = typeof aw === 'function' ? aw(Tc + KEL) : aw, TK = Tc + KEL;
  const both = (P) => { const fug = formerFugacities(f, P, Tc), w = a * (1 - (solubility ? gasSolubility(fug, P, Tc, molality).total : 0)); return [hydrateState('sI', TK, P, fug, w, table), hydrateState('sII', TK, P, fug, w, table)]; };
  const drive = (P) => { const [x, y] = both(P); return Math.max(x.drive, y.drive); };
  const grid = logspace(Pmin, Pmax, 36);
  if (drive(grid[0]) > 0) return { P: Pmin, structure: null, occupancy: null };
  let k = 1; while (k < grid.length && !(drive(grid[k]) > 0)) k++;
  if (k === grid.length) return null;
  let lo = Math.log(grid[k - 1]), hi = Math.log(grid[k]);
  for (let i = 0; i < 50 && hi - lo > 1e-6; i++) { const m = 0.5 * (lo + hi); if (drive(Math.exp(m)) > 0) hi = m; else lo = m; }
  const P = Math.exp(0.5 * (lo + hi)), [sI, sII] = both(P), structure = sII.drive > sI.drive ? 'sII' : 'sI';
  return { P, structure, occupancy: occupancy(structure, structure === 'sII' ? sII : sI) };
}
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
export function tunedFluid(spec, tune = {}, kijRows = []) { const ch = characterise(spec, tune), nKij = applyKij(ch.fluid, kijRows); return { ch, f: ch.fluid, nKij }; }
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

function engine(v, ctx = {}) {
  const prog = (x, m) => ctx.progress?.(x, m), warnings = [], recommendations = [], plots = [], tables = [], balances = [];
  const { spec, library } = resolveSpec(ctx.fluid, v.fluidSource), tune = tuneOf(v, spec);
  const tRes = clamp(num(v.tRes, 90), -20, 250), pRes = clamp(num(v.pRes, 300), 2, 1400), pRef = clamp(num(v.pRef, 100), 2, 1000), tSea = clamp(num(v.tSeabed, 4), -5, 60), pArr = clamp(num(v.pArr, 25), 1.05, 600), tArr = clamp(num(v.tArr, 30), -20, 150), tIn = clamp(num(spec.Tin, 70), -20, 250);
  if (!(num(v.pRes, 300) > 0)) throw new Error('The reservoir pressure must be a positive absolute pressure in bara.');
  prog(0.02, 'Characterising the fluid');
  const { ch, f, nKij } = tunedFluid(spec, tune, v.kijTable), opts = ch.opts;
  lastModel = { spec: JSON.parse(JSON.stringify(spec)), kij: JSON.parse(JSON.stringify(Array.isArray(v.kijTable) ? v.kijTable : [])) };
  const z7 = f.comps.reduce((s, c) => s + (c.pseudo ? c.z : 0), 0), isPseudo = f.comps.map((c) => !!c.pseudo);

  // -- property table for the flow solvers (kernel grid, tuned model)
  const nP = clamp(Math.round(num(v.nP, 22)), 6, 80), nT = clamp(Math.round(num(v.nT, 17)), 5, 70);
  const tPmin = clamp(num(v.tblPmin, 1), 0.5, 50), tPmax = clamp(num(v.tblPmax, 600), 100, 1500), tTmin = clamp(num(v.tblTmin, -30), -60, 20), tTmax = clamp(num(v.tblTmax, 170), 60, 300);
  if (num(v.tblPmax, 600) <= num(v.tblPmin, 1) || num(v.tblTmax, 170) <= num(v.tblTmin, -30)) throw new Error('The upper limit of the property-table range must be above its lower limit.');
  const table = buildTable(spec, { nP, nT, Pmin: tPmin, Pmax: tPmax, Tmin: tTmin, Tmax: tTmax, opts, onProgress: (x) => prog(0.03 + 0.25 * x, 'Building the pressure–temperature property table') });
  if (nKij > 0) { refillTable(table, spec, f); warnings.push({ level: 'info', msg: `${nKij} binary interaction parameter(s) are overridden by the user table; the property table was recomputed with them.` }); }
  for (const c of [['reservoir', pRes, tRes], ['cold reference', pRef, tSea], ['arrival', pArr, tArr]]) if (c[1] > tPmax || c[1] < tPmin || c[2] > tTmax || c[2] < tTmin) warnings.push({ level: 'warn', msg: `The ${c[0]} condition (${c[1]} bara, ${c[2]} °C) lies outside the property-table range (${tPmin}–${tPmax} bara, ${tTmin}–${tTmax} °C): the flow solvers will clamp to the edge of the table.` });
  const std = stdFlash(f), rates = table.rates, hasLiq = std.vOil > 0 && std.beta < 1, hasOil = hasLiq && std.gor < 20000, gasSG = std.gasSG ?? f.MW / MW_AIR;

  // -- phase envelope
  prog(0.3, 'Tracing the phase envelope');
  const first = saturationPoint(f, tRes), dewLike = first.P === null || first.type === 'dew';
  const env = traceEnvelope(f, { n: clamp(Math.round(num(v.nEnv, 40)), 10, 300), qualities: dewLike ? [0.5, 0.9, 0.99] : [0.1, 0.5, 0.9] });
  let sat = first;
  if (first.P !== null && !first.X) { const s2 = saturationPoint(f, tRes, { seed: seedFor(env, tRes) }); if (s2.P !== null && s2.X) sat = s2; }
  const psat = sat.P !== null && !sat.capped ? sat.P : null, psatType = psat !== null ? sat.type : null;
  if (sat.capped) warnings.push({ level: 'warn', msg: `The fluid is still two-phase at 1,200 bara and ${tRes} °C: no saturation pressure is reported. Check the heavy-end description (C7+ molar mass ${spec.c7MW} g/mol).` });

  // -- states at the key conditions
  prog(0.42, 'Flashing the key conditions');
  const conds = [{ name: 'Reservoir', P: pRes, T: tRes }, { name: 'Flowline inlet', P: pRef, T: tIn }, { name: 'Cold reference (seabed)', P: pRef, T: tSea }, { name: 'Arrival', P: pArr, T: tArr }, { name: 'Standard', P: P_STD, T: T_STD }];
  for (const c of conds) { c.s = fluidState(f, c.P, c.T); c.st = f.n > 1 ? stability(f, c.P, c.T + KEL) : { stable: true, tpd: 0 }; }
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
  const iso = isoP.map((P) => ({ P, s: fluidState(f, P, tRes) })), isoCold = isoP.map((P) => ({ P, s: props(f, P, tSea, { thermal: false }) }));
  const ser = (name, arr, get, extra = {}) => { const x = [], y = []; for (const r of arr) { const q = get(r); if (Number.isFinite(q)) { x.push(r.P ?? r.T); y.push(q); } } return { name, x, y, ...extra }; };
  const isoT = linspace(Math.min(tSea, 0) - 5, Math.max(tRes, tIn) + 10, Math.max(8, Math.round(nIso * 0.7))).map((T) => ({ T, s: props(f, pRef, T) }));

  // -- aqueous phase
  prog(0.64, 'Aqueous phase and hydrate equilibrium');
  const models = { act: ['NRTL', 'UNIQUAC', 'Wilson', 'ideal'].includes(v.actModel) ? v.actModel : 'NRTL', elec: ['pitzer', 'davies', 'edh', 'dh'].includes(v.elecModel) ? v.elecModel : 'pitzer' };
  const aqK = aqueous(spec), aqCase = { S: aqK.S, inhId: aqK.inhId, inhWt: aqK.inhWt };
  const custom = v.aqSource === 'custom', brine = custom ? brineFromIons(v.aqIons) : null, aqCus = custom ? { S: brine ? clamp(brine.S, 0, 260) : clamp(num(v.salinityIn, 3.5), 0, 26) * 10, inhId: INHIBITORS[v.inhIn] && v.inhIn !== 'none' && num(v.inhWtIn, 0) > 0 ? v.inhIn : 'none', inhWt: clamp(num(v.inhWtIn, 0), 0, 90) } : null;
  const awOf = (aq) => (TK) => waterActivity(aq, TK, models).aw, molOf = (aq) => waterActivity(aq, 277, models).molality;
  const awCase = waterActivity(aqCase, tSea + KEL, models);
  if (awCase.idealInh) warnings.push({ level: 'warn', msg: `No activity-coefficient parameters are held for ${INHIBITORS[aqCase.inhId].name}: its solution with water is treated as ideal (Raoult's law), which overstates the inhibition of alcohols and understates that of glycols by up to about 20 %.` });
  if (models.elec !== 'pitzer' && molOf(aqCase) > 0.7) warnings.push({ level: 'warn', msg: `The ${models.elec === 'davies' ? 'Davies' : 'Debye–Hückel'} model is used at ${molOf(aqCase).toFixed(2)} mol/kg, beyond its range of validity (about 0.5 mol/kg): select the Pitzer model for this brine.` });

  // -- hydrate curve
  const nHyd = clamp(Math.round(num(v.nHyd, 24)), 20, 120), hP = logspace(5, 500, nHyd), useVdw = v.hydModel !== 'corr';
  const dep = hydrateDepression(aqK), corr = (P) => hydrateT0(P, gasSG);
  const curve = (aq) => { let g = 8; const T = [], st = []; for (const P of hP) { const r = hydrateTofP(f, P, aq ? awOf(aq) : 1, { molality: aq ? molOf(aq) : 0, guess: g }); if (!r || !Number.isFinite(r.T)) return null; T.push(r.T); st.push(r.structure); g = r.T; } return { T, st }; };
  let hyd0 = useVdw ? curve(null) : null, hydC = hyd0 ? curve(aqCase) : null, hydModelUsed = 'van der Waals–Platteeuw';
  if (!hyd0 || !hydC) {
    if (useVdw) warnings.push({ level: 'warn', msg: 'The van der Waals–Platteeuw solution did not exist over the whole 5–500 bara range (no hydrate former in the fluid or no root): the published curve falls back to the gas-gravity correlation with the Nielsen–Bucklin depression.' });
    hyd0 = { T: hP.map(corr), st: hP.map(() => 'sII') }; hydC = { T: hP.map((P) => corr(P) - dep), st: hyd0.st }; hydModelUsed = 'gas-gravity correlation (Motiee) with Nielsen–Bucklin depression';
  }
  const hydX = custom && hydModelUsed.startsWith('van') ? curve(aqCus) : null;
  const lnHP = hP.map(Math.log), hydAt = (arr, P) => interp1(lnHP, arr, Math.log(clamp(P, 5, 500)));
  const vdw = hydModelUsed.startsWith('van'), refFresh = vdw ? hydrateTofP(f, pRef, 1) : null, refCase = vdw ? hydrateTofP(f, pRef, awOf(aqCase), { molality: molOf(aqCase) }) : null;
  const tHyd0 = refFresh ? refFresh.T : corr(pRef), tHyd = refCase ? refCase.T : corr(pRef) - dep, subcool = tHyd - tSea, depression = tHyd0 - tHyd, structure = refCase?.structure || refFresh?.structure || 'sII';
  const margin = clamp(num(v.margin, 3), 0, 15), inhDesign = aqCase.inhId !== 'none' ? aqCase.inhId : INHIBITORS[v.inhDesign] && v.inhDesign !== 'none' ? v.inhDesign : 'MEG';
  const tWith = (w, id = inhDesign) => { if (!vdw) return corr(pRef) - hydrateDepression({ S: aqCase.S, inhWt: w, inh: INHIBITORS[id] }); const aq = { S: aqCase.S, inhId: w > 0 ? id : 'none', inhWt: w }, r = hydrateTofP(f, pRef, awOf(aq), { molality: molOf(aq), guess: tHyd }); return r ? r.T : -45; };
  let inhReq = 0; const target = tSea - margin;
  const tNear = (w) => { if (!vdw) return tWith(w); const aq = { S: aqCase.S, inhId: w > 0 ? inhDesign : 'none', inhWt: w }, r = hydrateTofP(f, pRef, awOf(aq), { molality: molOf(aq), guess: target, passes: 1 }); return r ? r.T : -45; }; // guest fugacities taken at the target temperature
  if (tWith(0) > target) { const g = (w) => tNear(w) - target; inhReq = g(85) > 0 ? 85 : brent(g, 0, 85, 1e-3, 40); }
  const inhNB = Math.max(0, inhibitorFor(tHyd0 - target, inhDesign, aqCase.S)), doseRows = [0, 10, 20, 30, 40, 50, 60].map((w) => ({ w, model: tHyd0 - tWith(w), hamm: hydrateDepression({ S: aqCase.S }) + hammerschmidt(w, inhDesign), nb: hydrateDepression({ S: aqCase.S, inhWt: w, inh: INHIBITORS[inhDesign] }) }));
  const pDiss = vdw ? hydratePofT(f, tSea, awOf(aqCase), { molality: molOf(aqCase) }) : null, pDissSea = pDiss ? pDiss.P : null;
  const occ = refCase?.occupancy || refFresh?.occupancy || null;
  // sensitivity of the inhibitor requirement to pressure, seabed temperature, margin and salinity
  const reqAt = (P, tS, S, mg) => {
    const tg = tS - mg, tw = (w) => { if (!vdw) return corr(P) - hydrateDepression({ S, inhWt: w, inh: INHIBITORS[inhDesign] }); const aq = { S, inhId: w > 0 ? inhDesign : 'none', inhWt: w }, r = hydrateTofP(f, P, awOf(aq), { molality: molOf(aq), guess: tg, passes: 1 }); return r ? r.T : -45; };
    if (!(tw(0) > tg)) return 0; const g = (w) => tw(w) - tg; return g(85) > 0 ? 85 : brent(g, 0, 85, 1e-2, 30);
  };
  const sens = [['Base case', pRef, tSea, aqCase.S, margin], ['Pressure +50 %', Math.min(pRef * 1.5, 1000), tSea, aqCase.S, margin], ['Pressure −50 %', Math.max(pRef * 0.5, 2), tSea, aqCase.S, margin], ['Seabed 2 °C colder', pRef, tSea - 2, aqCase.S, margin], ['Seabed 2 °C warmer', pRef, tSea + 2, aqCase.S, margin], ['No safety margin', pRef, tSea, aqCase.S, 0], ['Fresh (condensed) water', pRef, tSea, 0, margin], ['Salinity doubled', pRef, tSea, Math.min(2 * aqCase.S, 260), margin]].map((r, i) => ({ name: r[0], P: r[1], tS: r[2], S: r[3], mg: r[4], w: i === 0 ? inhReq : reqAt(r[1], r[2], r[3], r[4]) }));

  // -- water content, gas solubility
  const wcBuk = waterContent(pArr, tArr) * awOf(aqCase)(tArr + KEL), wcRaoult = waterContentRaoult(pArr, tArr, awOf(aqCase)(tArr + KEL));
  const fugRef = formerFugacities(f, pRef, tSea), sol = gasSolubility(fugRef, pRef, tSea, molOf(aqCase));
  // water distribution between the gas and the aqueous phase (kg/d) for the case rates
  const waterDist = [conds[1], conds[2], conds[3]].map((c) => { const wc = waterContent(c.P, c.T) * awOf(aqCase)(c.T + KEL), q = lookup(table, c.P, c.T), sat = wc * ((q.wG * rates.mHC) / (q.mwG * 1e-3)) * VM_STD * 86400, tot = rates.mW * 86400 * (1 - aqCase.inhWt / 100) * (1 - aqCase.S / 1000), inGas = Math.min(sat, tot); return { name: c.name, P: c.P, T: c.T, wc, tot, inGas, free: tot - inGas, satCap: sat }; });

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
  tables.push({ title: 'Hydrate dissociation curve', columns: ['P (bara)', 'T fresh water (°C)', 'T case aqueous phase (°C)', 'Structure', 'Gas-gravity correlation (°C)', ...(hydX ? ['T custom aqueous phase (°C)'] : [])], rows: hP.map((P, i) => [cell(P), cell(hyd0.T[i]), cell(hydC.T[i]), hydC.st[i], cell(corr(P)), ...(hydX ? [cell(hydX.T[i])] : [])]), note: `${hydModelUsed}. ${vdw ? 'Langmuir constants and reference properties of Munck et al. (1988); guest fugacities from the equation of state; water activity from the selected activity and electrolyte models and from dissolved gas.' : ''}` });
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
  const kind = !hasOil || rates.gor > 50000 ? 'dry gas' : psatType === 'dew' || (psat === null && rates.gor > 900) ? (rates.gor > 10000 ? 'wet gas' : 'gas condensate') : rates.gor > 350 ? 'volatile oil' : 'black oil';
  const summary = `${spec.name || 'The case fluid'} behaves as a ${kind} (${f.eosId}, ${f.n} components): ${psat !== null ? `${psatType} point ${F(psat)} bara at ${F(tRes)} °C` : `no saturation pressure at ${F(tRes)} °C`}${rates.gor !== null && hasOil ? `, GOR ${F(rates.gor)} Sm³/Sm³, ${F(std.api)} °API` : `, gas gravity ${F(gasSG, 3)}`}; hydrate equilibrium ${F(tHyd)} °C at ${F(pRef)} bara (${subcool > 0 ? `${F(subcool)} °C subcooling at the seabed` : 'outside the hydrate region at the seabed'})${wat !== null ? `, wax appearance ${F(wat)} °C` : ''}.`;
  prog(1, 'Done');
  return { summary, kpis, warnings, recommendations, plots, tables, balances, outputs };
}

// ---- calibration model ----------------------------------------------------------------------------------------------------------
const calCache = new Map();
/** Predictions for one laboratory point: saturation pressure at calT, liquid density, viscosity and solution GOR at (calP, calT). */
export function calibrationModel(v) {
  const src = FLUID_LIBRARY[v.fluidSource] ? resolveSpec(null, v.fluidSource).spec : lastModel?.spec || mergeSpec(null), tune = tuneOf(v, src);
  const T = clamp(num(v.calT, num(v.tRes, 90)), -20, 250), P = clamp(num(v.calP, num(v.pRes, 300)), 1.02, 1400);
  const key = JSON.stringify([v.fluidSource, src.comp, src.c7MW, src.c7SG, tune, v.kijTable || null]);
  let m = calCache.get(key);
  if (!m) { m = { f: tunedFluid(src, tune, v.kijTable).f, sat: new Map() }; if (calCache.size > 400) calCache.clear(); calCache.set(key, m); }
  const f = m.f;
  let sat = m.sat.get(T);
  if (!sat) { const seedKey = `${v.fluidSource}|${T}`, seed = calCache.get(seedKey); sat = saturationPoint(f, T, { seed: seed || null }); if (sat.X) calCache.set(seedKey, sat.X); m.sat.set(T, sat); }
  const s = pv(f, P, T), liq = s.phase === 'gas' ? null : s.oil;
  let rs = null;
  if (liq) { const fl = stdFlash(withZ(f, s.x)); rs = fl.vOil > 0 ? fl.vGas / fl.vOil : null; }
  return { psat: sat.P ?? NaN, rho: (liq || s.gas).rho, mu: (liq || s.gas).mu * 1e3, rs: rs ?? NaN };
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
  chk('Hydrate dissociation: methane at 53 bara', 280, hyd({ C1: 100 }, 53), 1, 'Deaton & Frost type data: 280 K at 5.3 MPa (structure I)');
  chk('Hydrate dissociation: methane at 26.5 bara', 273.2, hyd({ C1: 100 }, 26.5), 1, 'About 2.6–2.7 MPa at 273.2 K');
  chk('Hydrate dissociation: ethane at 5.3 bara', 273.7, hyd({ C2: 100 }, 5.3), 1, '0.53 MPa at 273.7 K (structure I)');
  chk('Hydrate dissociation: propane at 5.5 bara', 278.3, hyd({ C3: 100 }, 5.5), 1, '0.55 MPa at 278.3 K (structure II)');
  chk('Hydrate dissociation: carbon dioxide at 45 bara', 283, hyd({ CO2: 100 }, 45), 1, '4.5 MPa at 283 K, with the dissolved CO₂ lowering the water activity');
  chk('Structure selection: 1 % propane turns methane hydrate into structure II', 2, hydrateTofP(makeFluid({ comp: { C1: 99, C3: 1 } }), 43.6).structure === 'sII' ? 2 : 1, 0, 'Lowest dissociation pressure decides the structure');
  chk('Inhibitor response: 50 wt % MEG on methane hydrate at 99 bara', 22.6, (() => { const m = makeFluid({ comp: { C1: 100 } }); return hydrateTofP(m, 98.9).T - hydrateTofP(m, 98.9, (TK) => waterActivity({ inhId: 'MEG', inhWt: 50 }, TK).aw).T; })(), 2, 'Measured depression of about 22–23 K (Robinson & Ng type data)');
  chk('Debye–Hückel osmotic constant of water at 25 °C', 0.3915, debyeHuckel(298.15).Aphi, 0.002, 'From the density and dielectric constant of water');
  chk('Pitzer osmotic coefficient of 1 mol/kg NaCl at 25 °C', 0.9355, osmoticNaCl('pitzer', 1), 0.003, 'Robinson & Stokes tabulation');
  chk('Pitzer osmotic coefficient of 6 mol/kg NaCl at 25 °C', 1.2706, osmoticNaCl('pitzer', 6), 0.01, 'Robinson & Stokes tabulation');
  chk('Pitzer mean activity coefficient of 1 mol/kg NaCl at 25 °C', 0.657, gammaNaCl('pitzer', 1), 0.005, 'Robinson & Stokes tabulation');
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
  return out;
}

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
export default {
  id: 'pvt', num: 1, title: 'Fluid, PVT & Phase Behaviour', short: 'Fluid · PVT', icon: '🧪',
  tagline: 'Equation-of-state characterisation, phase envelope, laboratory experiments, properties, hydrate, wax and aqueous-phase thermodynamics of the case fluid.',
  description: 'The case fluid is characterised with a cubic equation of state (C7+ split, alternative critical-property correlations, tuning multipliers and interaction parameters) and flashed at constant PT, PH, PS and TV. The phase envelope is traced by continuation through the critical point, the standard PVT-cell experiments are simulated for comparison with a laboratory report, and every property the flow solvers need is tabulated. Hydrate equilibrium follows van der Waals–Platteeuw with the water activity from activity-coefficient and electrolyte models; wax and asphaltenes are screened from the same characterisation.',
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
    'pvt-cell experiments', 'bubble-point measurements', 'dew-point measurements', 'density measurements', 'viscosity measurements', 'hydrate equilibrium experiments', 'hydrate dissociation experiments', 'meg/methanol inhibitor', 'multicomponent reservoir-fluid',
  ],
  referenceOnly: ['wax/asphaltene solid-phase equilibrium', 'gibbs-energy minimization', 'benedict', 'cubic-plus-association', 'saft', 'gerg', 'pedersen viscosity', 'kihara', 'lennard', 'structure-h', 'helmholtz', 'slim-tube', 'where solids are permitted', 'asphaltene or other solid phase fractions', 'volume-shift regression', 'speed-of-sound calibration', 'gas-density calibration', 'water-density calibration', 'gas-viscosity calibration', 'enthalpy/heat-capacity calibration', 'water-content calibration', 'compressibility calibration', 'z-factor calibration', 'tension calibration', 'formation-volume-factor calibration', 'hydrate-equilibrium calibration', 'hydrate dissociation p-t calibration', 'inhibitor-response calibration', 'salinity correction calibration', 'water-activity calibration', 'saline-water hydrate', 'calorimetric', 'interfacial-tension measurements', 'water-content measurements', 'high-pressure', 'phase-envelope measurements', 'compressibility measurements'],
  equationsNote: 'Cubic equations of state (Peng–Robinson 1978, Soave–Redlich–Kwong, Redlich–Kwong, van der Waals) with van der Waals one-fluid mixing rules, Péneloux-type volume translation and two-phase vapour–liquid equilibrium; the hydrocarbon flash carries no water and no third (aqueous or solid) phase. Riazi–Daubert and Twu critical properties act as plus-fraction-average ratios to Kesler–Lee because the shared kernel accepts uniform multipliers only; volume shift is on/off, not a regressed parameter. The critical point is located on the traced saturation line, not by the Heidemann–Khalil criterion. Viscosity is Lohrenz–Bray–Clark (a corresponding-states residual-viscosity correlation) with Lee–Gonzalez–Eakin and Beggs–Robinson as cross-checks; the Pedersen corresponding-states model is not implemented. Lee–Kesler is used for the gas compressibility cross-check only. Hydrates: van der Waals–Platteeuw for structures I and II with the Langmuir constants and reference properties of Munck et al. (1988), guest fugacities from the EOS and water activity from NRTL / UNIQUAC / Wilson (methanol, MEG), Pitzer / Davies / Debye–Hückel (NaCl equivalent) and Henry\'s-law gas solubility; checked here against pure-gas dissociation data to within about 1.2 K up to 200 bara and 2 K at 450 bara; the stable structure is the one with the lower water chemical potential (no general multiphase Gibbs minimisation); structure H, Kihara cell-potential integration and hydrate formation from a liquid-water-free gas are not modelled. Water dissolved in the liquid hydrocarbon is neglected in the water distribution. Ethanol, DEG and TEG solutions are treated as ideal. Methanol and MEG low-temperature activity parameters were regressed to ice-point data in this work. Multi-ion brines are reduced to NaCl of equal ionic strength. Wax: ideal-solution solid–liquid equilibrium (Won melting properties, Pedersen wax-forming fraction) on a carbon-number distribution; asphaltenes: de Boer and colloidal-instability screening only. SAFT, CPA, GERG and BWR-type equations are listed for reference and not solved.',

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
      { key: 'scnModel', label: 'Carbon-number distribution (wax, plot)', type: 'select', value: 'exp', options: opt({ exp: 'Exponential (Pedersen)', gamma: 'Gamma (Whitson)' }) },
      { key: 'gammaAlpha', label: 'Gamma shape α', unit: '–', value: 1, min: 0.5, max: 3, showIf: (v) => v.scnModel === 'gamma' },
    ] },
    { group: 'Aqueous phase, hydrate and wax models', tab: 'setup', fields: [
      { key: 'hydModel', label: 'Hydrate model', type: 'select', value: 'vdwp', options: opt({ vdwp: 'van der Waals–Platteeuw (Munck parameters)', corr: 'Gas-gravity correlation + Nielsen–Bucklin' }) },
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
    { name: 'Lean gas tie-back with MEG injection', values: { fluidSource: 'leanGas', tRes: 95, pRes: 280, pRef: 150, pArr: 70, tArr: 8, useLab: false, swellGas: 'none', aqSource: 'custom', salinityIn: 0.5, inhIn: 'MEG', inhWtIn: 45, inhDesign: 'MEG' } },
    { name: 'Sour CO₂-rich gas, SRK and methanol', values: { fluidSource: 'sourGas', eosSel: 'SRK', tRes: 100, pRes: 300, pRef: 120, pArr: 60, tArr: 10, useLab: false, swellGas: 'none', inhDesign: 'MeOH', aqSource: 'custom', salinityIn: 2, inhIn: 'MeOH', inhWtIn: 25 } },
    { name: 'Reference oil tuned to the laboratory report', values: { fluidSource: 'case', tcMult: 1.025, pcMult: 0.97, kijScale: 1.15, vcMult: 1.06 } },
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
    note: 'Regress the C7+ critical-property multipliers, the interaction scale and the viscosity critical-volume multiplier to saturation pressure, liquid density, solution gas–oil ratio and viscosity. The model uses the fluid of the most recent run (the reference fluid before any run). Fit Tc, Pc and kij to saturation pressure, density and GOR first; the critical-volume multiplier only moves viscosity.',
    params: [
      { key: 'tcMult', label: 'C7+ Tc multiplier', lo: 0.9, hi: 1.1 }, { key: 'pcMult', label: 'C7+ Pc multiplier', lo: 0.8, hi: 1.2 }, { key: 'wMult', label: 'C7+ acentric-factor multiplier', lo: 0.8, hi: 1.2 },
      { key: 'kijScale', label: 'Binary-interaction scale', lo: 0, hi: 2.5 }, { key: 'vcMult', label: 'C7+ Vc multiplier (viscosity)', lo: 0.8, hi: 1.4 },
    ],
    columns: [{ key: 'calT', label: 'Temperature', unit: '°C' }, { key: 'calP', label: 'Pressure', unit: 'bara' }, { key: 'psat', label: 'Saturation pressure at T', unit: 'bara' }, { key: 'rho', label: 'Liquid density at P, T', unit: 'kg/m³' }, { key: 'rs', label: 'Solution GOR at P, T', unit: 'Sm³/Sm³' }, { key: 'mu', label: 'Liquid viscosity at P, T', unit: 'mPa·s' }],
    targets: [{ key: 'psat', label: 'Saturation pressure', unit: 'bara' }, { key: 'rho', label: 'Liquid density', unit: 'kg/m³' }, { key: 'rs', label: 'Solution GOR', unit: 'Sm³/Sm³' }, { key: 'mu', label: 'Liquid viscosity', unit: 'mPa·s' }],
    model: calibrationModel,
    sample: CAL_SAMPLE,
    validationSample: CAL_VALID,
  },

  verify: verifyChecks,

  views: [{ id: 'flash', label: 'Flash calculator', tip: 'Interactive PT / PH / PS flash with phase compositions and properties', render: flashView }],
};
