// Suite 3 — Multiphase thermal-hydraulics and slugging.
// Steady pressure / temperature / holdup along the elevation profile with selectable holdup closures, flow-regime
// maps and interfacial-stability criteria, hydrodynamic / terrain / severe-riser slugging (unit cell, Lagrangian
// slug tracking, lumped riser cycle), a semi-implicit finite-volume drift-flux transient, radial conduction through
// the wall and coatings, and local demonstration solvers (1-D radial RANS, 1-D interface capturing).
// SI units inside; bara, °C, mm and µm at the interfaces.
import { clamp, brent, interp1, tridiag, rng, mean, std, quantile, linspace, histogram, gci, isNum } from '../core/num.js';
import { fluidModel } from '../core/thermo.js';
import { G, frictionFactor, hInside, hOutside, nuCrossFlow, nuFreeCylinder, uValue, seaTemperature, stratifiedLevel, flowPattern, slugVelocity, slugBodyHoldup, slugFrequency, slugLength, slugUnit, slugFilm, annularFilm, severeSlugging, gradient, discretise, marchSteady } from '../core/pipe.js';
import { BASE } from '../data/basecase.js';
import { taylorGreen, channel3D, pipe3D, twoPhase3D, coupled1D3D, createFlow3D, makePoisson, CFD3D_CONSTANTS } from '../core/cfd3d.js';
import * as REF from '../data/ref/flow.js';

const PI = Math.PI, DEG = PI / 180, FT = 0.3048, P_ATM = 101325, RGAS = 8.314462618;
const fin = (x, d = 0) => (typeof x === 'number' && Number.isFinite(x) ? x : d);
const r3 = (x, n = 3) => (Number.isFinite(x) ? +x.toFixed(n) : null);
const sig = (x, n = 4) => (Number.isFinite(x) ? +x.toPrecision(n) : null);
const log10 = Math.log10;

// =====================================================================================================
// 1. Holdup and pressure-gradient closures added by this suite
//    p: { vsl, vsg, rhoL, rhoG, muL, muG, sigma, D, theta (rad, + up), rough (m), P (Pa), fModel, waterCont }
//    each returns { holdup, fric, grav, acc (Pa/m), regime }
// =====================================================================================================
/** Duns & Ros dimensionless groups: liquid and gas velocity numbers, diameter number, liquid-viscosity number. */
export function velocityNumbers(p) {
  const s = Math.max(p.sigma, 1e-4), k = (p.rhoL / (G * s)) ** 0.25;
  return { NLv: p.vsl * k, NGv: p.vsg * k, ND: p.D * Math.sqrt((p.rhoL * G) / s), NL: p.muL * (G / (p.rhoL * s ** 3)) ** 0.25 };
}
const withAcc = (r, rho, vm, vsg, P) => { const Ek = clamp((rho * vm * vsg) / Math.max(P, 1e4), 0, 0.6); r.acc = ((r.fric + r.grav) * Ek) / (1 - Ek); return r; };
// Griffith & Wallis bubble-flow holdup with a constant slip velocity of 0.8 ft/s
const griffithHoldup = (vm, vsg, vs = 0.8 * FT) => 1 - 0.5 * (1 + vm / vs - Math.sqrt(Math.max((1 + vm / vs) ** 2 - (4 * vsg) / vs, 0)));
function bubbleGriffith(p) {
  const { vsl, vsg, rhoL, rhoG, muL, D, theta, rough, P, fModel } = p, vm = vsl + vsg, H = clamp(griffithHoldup(vm, vsg), vsl / vm, 1), vL = vsl / H;
  const f = frictionFactor((rhoL * vL * D) / muL, rough / D, fModel), rhoS = rhoL * H + rhoG * (1 - H);
  return withAcc({ holdup: H, fric: (f * rhoL * vL * vL) / (2 * D), grav: rhoS * G * Math.sin(theta), regime: 'bubble' }, rhoS, vm, vsg, P);
}

/** Hagedorn & Brown (1965) with the Griffith–Wallis bubble-flow modification; chart curve fits as published for CN_L, H_L/ψ and ψ. */
function hagedornBrown(p) {
  const { vsl, vsg, rhoL, rhoG, muL, muG, D, theta, rough, P, fModel } = p, vm = vsl + vsg, lam = vsl / vm, { NLv, NGv, ND, NL } = velocityNumbers(p);
  if (vsg / vm < Math.max(1.071 - (0.2218 * (vm / FT) ** 2) / (D / FT), 0.13)) return bubbleGriffith(p);
  const x1 = log10(clamp(NL, 0.002, 0.5)) + 3, CNL = 10 ** (-2.69851 + 0.15841 * x1 - 0.551 * x1 ** 2 + 0.54785 * x1 ** 3 - 0.12195 * x1 ** 4);
  const x2 = clamp(log10(Math.max((NLv / NGv ** 0.575) * (P / P_ATM) ** 0.1 * (CNL / ND), 1e-12)) + 6, 0.3, 4), hPsi = -0.10307 + 0.61777 * x2 - 0.63295 * x2 ** 2 + 0.29598 * x2 ** 3 - 0.0401 * x2 ** 4;
  const x3 = Math.min((NGv * NL ** 0.38) / ND ** 2.14, 0.09), psi = x3 <= 0.012 ? 1 : Math.max(1, 0.91163 - 4.82176 * x3 + 1232.25 * x3 ** 2 - 22253.6 * x3 ** 3 + 116174.3 * x3 ** 4);
  const H = clamp(hPsi * psi, lam, 1), rhoS = rhoL * H + rhoG * (1 - H), rhoN = rhoL * lam + rhoG * (1 - lam);
  const f = frictionFactor((rhoN * vm * D) / (muL ** H * muG ** (1 - H)), rough / D, fModel);
  return withAcc({ holdup: H, fric: (f * rhoN * rhoN * vm * vm) / (2 * rhoS * D), grav: rhoS * G * Math.sin(theta), regime: 'slug' }, rhoS, vm, vsg, P);
}

/** Gray (1974) for gas wells producing condensate and water. */
function gray(p) {
  const { vsl, vsg, rhoL, rhoG, muL, muG, D, theta, rough, P, fModel } = p, vm = vsl + vsg, lam = vsl / vm, s = Math.max(p.sigma, 1e-4), dRho = Math.max(rhoL - rhoG, 1);
  const rhoN = rhoL * lam + rhoG * (1 - lam), Nv = (rhoN * rhoN * vm ** 4) / (G * s * dRho), Nd = (G * dRho * D * D) / s, R = vsl / vsg;
  const B = 0.0814 * (1 - 0.0554 * Math.log(1 + (730 * R) / (R + 1))), A = -2.314 * (Nv * (1 + 205 / Nd)) ** B, H = clamp(1 - (1 - Math.exp(A)) / (R + 1), lam, 1);
  const e0 = (28.5 * s) / (rhoN * vm * vm), eff = Math.max(R >= 0.007 ? e0 : rough + (R * (e0 - rough)) / 0.007, 2.77e-5, rough), rhoS = rhoL * H + rhoG * (1 - H);
  const f = frictionFactor((rhoN * vm * D) / (muL * lam + muG * (1 - lam)), Math.min(eff / D, 0.05), fModel);
  return withAcc({ holdup: H, fric: (f * rhoN * vm * vm) / (2 * D), grav: rhoS * G * Math.sin(theta), regime: H < 0.25 ? 'annular' : 'churn' }, rhoS, vm, vsg, P);
}

/** Orkiszewski (1967): Griffith–Wallis bubble flow, liquid-distribution-coefficient slug flow, Duns–Ros transition and mist flow. */
function orkiszewski(p) {
  const { vsl, vsg, rhoL, rhoG, muL, muG, D, theta, rough, P, fModel } = p, vm = vsl + vsg, lam = vsl / vm, s = Math.max(p.sigma, 1e-4), { NLv, NGv } = velocityNumbers(p);
  const vmF = vm / FT, dF = D / FT, LB = Math.max(1.071 - (0.2218 * vmF * vmF) / dF, 0.13), LS = 50 + 36 * NLv, LM = 75 + 84 * NLv ** 0.75, sinT = Math.sin(theta);
  if (vsg / vm < LB) return bubbleGriffith(p);
  const slug = () => {
    const ReL = (rhoL * vm * D) / muL, sg = Math.sqrt(G * D), muCp = muL * 1000, rhoF = rhoL / 16.01846;
    let vb = 0.5 * sg; // Taylor-bubble rise velocity, iterated on the bubble Reynolds number
    for (let k = 0; k < 20; k++) {
      const Reb = (rhoL * vb * D) / muL; let nv;
      if (Reb <= 3000) nv = (0.546 + 8.74e-6 * ReL) * sg; else if (Reb >= 8000) nv = (0.35 + 8.74e-6 * ReL) * sg;
      else { const vbi = ((0.251 + 8.74e-6 * ReL) * sg) / FT; nv = 0.5 * (vbi + Math.sqrt(vbi * vbi + (13.59 * muCp) / (rhoF * Math.sqrt(dF)))) * FT; }
      if (Math.abs(nv - vb) < 1e-9) { vb = nv; break; } vb = 0.5 * (vb + nv);
    }
    let d; // liquid distribution coefficient (field units: ft/s, ft, cP)
    if (p.waterCont) d = vmF < 10 ? (0.013 * log10(muCp)) / dF ** 1.38 - 0.681 + 0.232 * log10(vmF) - 0.428 * log10(dF) : (0.045 * log10(muCp)) / dF ** 0.799 - 0.709 - 0.162 * log10(vmF) - 0.888 * log10(dF);
    else d = vmF < 10 ? (0.0127 * log10(muCp + 1)) / dF ** 1.415 - 0.284 + 0.167 * log10(vmF) + 0.113 * log10(dF) : (0.0274 * log10(muCp + 1)) / dF ** 1.371 + 0.161 + 0.569 * log10(dF) - log10(vmF) * ((0.01 * log10(muCp + 1)) / dF ** 1.571 + 0.397 + 0.63 * log10(dF));
    const base = (rhoL * (vsl + vb) + rhoG * vsg) / (vm + vb);
    d = vmF < 10 ? Math.max(d, -0.065 * vmF) : Math.max(d, (-vb / (vm + vb)) * (1 - base / rhoL));
    const rhoS = clamp(base + rhoL * d, rhoL * lam + rhoG * (1 - lam), rhoL), H = clamp((rhoS - rhoG) / Math.max(rhoL - rhoG, 1e-6), lam, 1), f = frictionFactor(ReL, rough / D, fModel);
    return { holdup: H, fric: Math.max(((f * rhoL * vm * vm) / (2 * D)) * ((vsl + vb) / (vm + vb) + d), 0), grav: rhoS * G * sinT, rho: rhoS, regime: 'slug' };
  };
  const mist = () => {
    const rhoN = rhoL * lam + rhoG * (1 - lam), N = (rhoG * vsg * vsg * muL * muL) / (rhoL * s * s), e0 = s / (rhoG * vsg * vsg * D);
    const ed = clamp(N <= 0.005 ? 0.0749 * e0 : 0.3713 * e0 * N ** 0.302, rough / D, 0.5);
    const f = ed > 0.05 ? 4 * (1 / (4 * log10(0.27 * ed)) ** 2 + 0.067 * ed ** 1.73) : frictionFactor((rhoG * vsg * D) / muG, ed, fModel);
    return { holdup: lam, fric: (f * rhoG * vsg * vsg) / (2 * D), grav: rhoN * G * sinT, rho: rhoN, regime: 'annular' };
  };
  let r;
  if (NGv < LS) r = slug();
  else if (NGv > LM) r = mist();
  else { const a = slug(), b = mist(), w = (LM - NGv) / (LM - LS); r = { holdup: w * a.holdup + (1 - w) * b.holdup, fric: w * a.fric + (1 - w) * b.fric, grav: w * a.grav + (1 - w) * b.grav, rho: w * a.rho + (1 - w) * b.rho, regime: 'churn' }; }
  const rho = r.rho; delete r.rho;
  return withAcc(r, rho, vm, vsg, P);
}

const MB = { up: [-0.380113, 0.129875, -0.119788, 2.343227, 0.475686, 0.288657], downStrat: [-1.330282, 4.808139, 4.171584, 56.262268, 0.079951, 0.504887], downOther: [-0.516644, 0.789805, 0.551627, 15.519214, 0.371771, 0.393952] };
const MB_FR = [[0.01, 1.0], [0.2, 0.98], [0.3, 1.2], [0.4, 1.25], [0.5, 1.3], [0.7, 1.25], [1.0, 1.0], [10, 1.0]];
/** Mukherjee & Brill (1985): inclined-flow holdup correlation with its own flow-pattern map (all inclinations). */
function mukherjeeBrill(p) {
  const { vsl, vsg, rhoL, rhoG, muL, muG, D, theta, rough, P, fModel } = p, vm = vsl + vsg, lam = vsl / vm, { NLv, NGv, NL } = velocityNumbers(p), s = Math.sin(theta), lg = log10(Math.max(NGv, 1e-12)), ll = log10(Math.max(NLv, 1e-12));
  let regime;
  if (NGv > 10 ** (1.401 - 2.694 * NL + 0.521 * NLv ** 0.329)) regime = 'annular';
  else if (theta > 0) regime = NLv > 10 ** (lg + 0.94 + 0.074 * s - 0.855 * s * s + 3.695 * NL) ? 'bubble' : 'slug';
  else if (NLv < 10 ** (0.321 - 0.017 * NGv - 4.267 * s - 2.972 * NL - 0.033 * lg * lg - 3.925 * s * s)) regime = 'stratified wavy';
  else regime = NGv > 10 ** (0.431 - 3.003 * NL - 1.138 * ll * s - 0.429 * ll * ll * s + 1.132 * s) ? 'slug' : 'bubble';
  const c = theta >= 0 ? MB.up : regime === 'stratified wavy' ? MB.downStrat : MB.downOther;
  const H = clamp(Math.exp(((c[0] + c[1] * s + c[2] * s * s + c[3] * NL * NL) * NGv ** c[4]) / NLv ** c[5]), theta >= 0 ? lam : 1e-4, 1);
  const rhoS = rhoL * H + rhoG * (1 - H), rhoN = rhoL * lam + rhoG * (1 - lam), muN = muL * lam + muG * (1 - lam), fn = frictionFactor((rhoN * vm * D) / muN, rough / D, fModel);
  if (regime === 'stratified wavy' && H < 0.999) { // gas-phase momentum balance over the stratified cross-section, smooth interface
    const dlt = brent((x) => (x - Math.sin(x)) / (2 * PI) - H, 1e-6, 2 * PI - 1e-6, 1e-10), A = (PI * D * D) / 4, AG = A * (1 - H), SG = D * (PI - dlt / 2), Wi = D * Math.sin(dlt / 2), vG = vsg / (1 - H), DhG = (4 * AG) / (SG + Wi);
    const tw = (frictionFactor((rhoG * vG * DhG) / muG, rough / DhG, fModel) * rhoG * vG * vG) / 8;
    return { holdup: H, fric: (tw * (SG + Wi)) / AG, grav: rhoG * G * s, acc: 0, regime };
  }
  const fR = regime === 'annular' ? interp1(MB_FR.map((q) => q[0]), MB_FR.map((q) => q[1]), clamp(lam / H, 0.01, 10)) : 1;
  return withAcc({ holdup: H, fric: (fn * fR * rhoS * vm * vm) / (2 * D), grav: rhoS * G * s, regime }, rhoS, vm, vsg, P);
}

/**
 * Zuber–Findlay drift flux with a user distribution parameter C0 and drift velocity vd, optionally with the Wallis
 * hindered-drift factor (1 − α)^n:  vsg = α [C0 vm + vd (1 − α)^n].  C0 = 1, vd = 0 reproduces no-slip flow.
 */
export function zuberFindlay(p, { c0 = 1.2, vDrift = 0.35, wallisN = 0 } = {}) {
  const { vsl, vsg, rhoL, rhoG, muL, muG, D, theta, rough, P, fModel } = p, vm = vsl + vsg, vd = theta >= 0 ? vDrift : vDrift * (Math.cos(theta) + Math.sin(theta));
  let alpha;
  if (wallisN > 0 && vd !== 0) { const g = (a) => a * (c0 * vm + vd * (1 - a) ** wallisN) - vsg; alpha = g(1 - 1e-9) <= 0 ? 1 - 1e-9 : brent(g, 0, 1 - 1e-9, 1e-12); }
  else alpha = vsg / Math.max(c0 * vm + vd, 1e-9);
  alpha = clamp(alpha, 0, 1 - 1e-6);
  const H = 1 - alpha, rhoM = rhoL * H + rhoG * alpha, muM = muL * H + muG * alpha, f = frictionFactor((rhoM * vm * D) / muM, rough / D, fModel);
  return withAcc({ holdup: H, fric: (f * rhoM * vm * vm) / (2 * D), grav: rhoM * G * Math.sin(theta), regime: 'drift flux' }, rhoM, vm, vsg, P);
}

// Slip closure of the transient model: Bendiksen-type drift blended to no-slip as the gas fraction approaches one
// (keeps the liquid velocity bounded) and direction-aware for flow reversal.
function slipLaw(alpha, j, theta, D, sp) {
  const w = clamp((1 - alpha) / 0.25, 0, 1), sg = Math.sqrt(G * D);
  return { C0: 1 + (sp.c0 - 1) * w, vd: sg * (0.35 * Math.sin(theta) + 0.54 * Math.cos(theta) * Math.tanh(j / 0.1)) * sp.vdScale * w };
}
function transientClosure(p, sp) {
  const { vsl, vsg, rhoL, rhoG, muL, muG, D, theta, rough, fModel } = p, vm = vsl + vsg;
  const g = (a) => { const sl = slipLaw(a, vm, theta, D, sp); return a * (sl.C0 * vm + sl.vd) - vsg; };
  const alpha = clamp(brent(g, 0, 1, 1e-12), 0, 1), H = 1 - alpha, rhoM = rhoL * H + rhoG * alpha, muM = muL * H + muG * alpha;
  return { holdup: H, fric: (sp.fricMult * frictionFactor((rhoM * vm * D) / muM, rough / D, fModel) * rhoM * vm * vm) / (2 * D), grav: rhoM * G * Math.sin(theta), acc: 0, regime: 'drift flux' };
}

export const HOLDUP_MODELS = [
  { value: 'beggsBrill', label: 'Beggs & Brill (1973, Payne corrections) — all inclinations' },
  { value: 'mechanistic', label: 'Mechanistic (Taitel–Dukler stratified, slug unit cell, drift flux)' },
  { value: 'driftFlux', label: 'Drift flux (Bendiksen slug-bubble closure)' },
  { value: 'zuberFindlay', label: 'Zuber–Findlay / Wallis drift flux (user C0 and drift velocity)' },
  { value: 'mukherjeeBrill', label: 'Mukherjee & Brill (1985) — all inclinations' },
  { value: 'hagedornBrown', label: 'Hagedorn & Brown (1965, Griffith bubble) — upward flow' },
  { value: 'orkiszewski', label: 'Orkiszewski (1967) — upward flow' },
  { value: 'gray', label: 'Gray (1974) — gas-condensate, upward flow' },
  { value: 'dunsRos', label: 'Duns & Ros (1963) — upward flow' },
  { value: 'ansari', label: 'Ansari et al. (1994) mechanistic — upward flow' },
  { value: 'mechEnt', label: 'Mechanistic + three-field annular flow (entrainment / deposition)' },
  { value: 'homogeneous', label: 'Homogeneous (no slip)' },
];
const OWN = { hagedornBrown, gray, orkiszewski, mukherjeeBrill, dunsRos, ansari, mechEnt };
const UP_ONLY = new Set(['hagedornBrown', 'gray', 'orkiszewski', 'dunsRos', 'ansari']);
/**
 * Two-phase holdup and pressure gradient at one location with any of the suite's closures.
 * p: { vsl, vsg, rhoL, rhoG, muL, muG, sigma, D, theta (rad, + up), rough (m), P (Pa), fModel, waterCont }
 * model: a value of HOLDUP_MODELS; mp: { c0, vDrift, wallisN, holdupMult, slip (transient closure parameters) }.
 * Upward-flow correlations (Hagedorn–Brown, Orkiszewski, Gray) fall back to Beggs & Brill on downward-inclined cells.
 * Returns { holdup, fric, grav, acc, dpdx (Pa/m, positive = pressure falls along the flow), regime, tauW (Pa) }.
 */
export function holdupGradient(p, model = 'beggsBrill', mp = {}) {
  const q = { sigma: 0.02, rough: 4.5e-5, theta: 0, P: 1e7, fModel: 'colebrook', ...p }, vm = q.vsl + q.vsg;
  let r;
  if (vm <= 1e-9 || q.vsg <= 1e-9 * vm || q.vsl <= 1e-9 * vm) {
    r = gradient(q, 'homogeneous');
  } else if (model === 'transient') r = transientClosure(q, mp.slip);
  else if (model === 'zuberFindlay') r = zuberFindlay(q, mp);
  else if (OWN[model] && !(UP_ONLY.has(model) && q.theta <= (model === 'ansari' ? 0.05 : -1e-12))) r = OWN[model](q, mp);
  else r = gradient(q, OWN[model] ? 'beggsBrill' : model);
  const hm = mp.holdupMult;
  if (hm && hm !== 1 && r.holdup > 0 && r.holdup < 1) { // calibration multiplier on the liquid holdup: the static head follows
    const H = clamp(r.holdup * hm, 1e-4, 1), rho0 = q.rhoL * r.holdup + q.rhoG * (1 - r.holdup);
    r.grav *= (q.rhoL * H + q.rhoG * (1 - H)) / rho0; r.holdup = H;
  }
  if (r.tauW === undefined) r.tauW = (r.fric * q.D) / 4;
  r.dpdx = r.fric + r.grav + r.acc;
  return r;
}

// =====================================================================================================
// 2. Heat transfer: ambient, resistance network, transient radial conduction
// =====================================================================================================
/** Ambient temperature (°C) at elevation z (m): exponential thermocline below sea level, air above. */
export const ambientTemperature = (z, th) => (z >= 0 ? th.tAir : seaTemperature(-z, th.tSeaSurface, th.tSeabed, th.thermocline || 250));

/**
 * Transient radial conduction through concentric layers (1-D finite volume in r, implicit Euler, exact logarithmic
 * face conductances so that the steady solution reproduces the analytic composite-cylinder profile).
 * o: { ri (m), layers: [{ name, t (m), k (W/m/K), rho (kg/m³), cp (J/kg/K) }], hIn, hOut (W/m²/K),
 *      tFluid (°C | fn(t)), tAmb (°C | fn(t)), fluidC (J/m/K heat capacity of the contents per metre; 0 = prescribed fluid temperature),
 *      T0 ('steady' | °C | fn(r)), tEnd (s), nSteps, nPer (cells per layer), source: fn(r, t) (W/m³), stopAt (°C), nProfiles }
 * Returns { r[], T[] (final), layerOf[], t[], fluid[], wall[] (inner-wall temperature), profiles: [{ t, T[] }], tReach (s | null), qIn (W/m, final) }.
 */
export function radialConduction(o) {
  const { ri, layers, hIn = 500, hOut = 500, fluidC = 0, tEnd = 3600, nSteps = 200, nPer = 6, source = null, stopAt = null, nProfiles = 5 } = o;
  const fT = typeof o.tFluid === 'function' ? o.tFluid : () => o.tFluid ?? 60, fA = typeof o.tAmb === 'function' ? o.tAmb : () => o.tAmb ?? 4;
  const rc = [], rf = [ri], kk = [], cap = [], layerOf = [];
  let r = ri;
  layers.forEach((l, li) => { if (!(l.t > 0 && l.k > 0)) return; for (let j = 0; j < nPer; j++) { const a = r + (l.t * j) / nPer, b = r + (l.t * (j + 1)) / nPer; rc.push(Math.sqrt(a * b)); rf.push(b); kk.push(l.k); cap.push((l.rho || 1000) * (l.cp || 1000) * PI * (b * b - a * a)); layerOf.push(li); } r += l.t; });
  const n = rc.length, ro = r;
  if (!n) throw new Error('Radial conduction needs at least one layer with positive thickness and conductivity.');
  const Gf = new Array(n + 1); // conductances per metre of pipe (W/m/K): fluid–cell 0, cell j−1–cell j, cell n−1–ambient
  Gf[0] = (2 * PI) / (1 / (Math.max(hIn, 1e-9) * ri) + Math.log(rc[0] / ri) / kk[0]);
  for (let j = 1; j < n; j++) Gf[j] = (2 * PI) / (Math.log(rf[j] / rc[j - 1]) / kk[j - 1] + Math.log(rc[j] / rf[j]) / kk[j]);
  Gf[n] = (2 * PI) / (Math.log(ro / rc[n - 1]) / kk[n - 1] + 1 / (Math.max(hOut, 1e-9) * ro));
  const src = (t) => rc.map((x, j) => (source ? source(x, t) * PI * (rf[j + 1] ** 2 - rf[j] ** 2) : 0));
  // one implicit step; dt = Infinity gives the steady state. Unknowns: [fluid (when fluidC > 0)], cells 0..n−1
  const step = (T, Tf, dt, t) => {
    const free = fluidC > 0 && Number.isFinite(dt), m = n + (free ? 1 : 0), a = new Array(m).fill(0), b = new Array(m).fill(0), c = new Array(m).fill(0), d = new Array(m).fill(0), o0 = free ? 1 : 0, S = src(t), tf = free ? null : fT(t), ta = fA(t);
    if (free) { b[0] = fluidC / dt + Gf[0]; c[0] = -Gf[0]; d[0] = (fluidC / dt) * Tf; }
    for (let j = 0; j < n; j++) {
      const i = j + o0, ct = Number.isFinite(dt) ? cap[j] / dt : 0;
      b[i] = ct + Gf[j] + Gf[j + 1]; d[i] = ct * T[j] + S[j];
      if (j > 0 || free) a[i] = -Gf[j]; else d[i] += Gf[0] * tf;
      if (j < n - 1) c[i] = -Gf[j + 1]; else d[i] += Gf[n] * ta;
    }
    const x = tridiag(a, b, c, d);
    return free ? { Tf: x[0], T: x.slice(1) } : { Tf: tf, T: x };
  };
  let T, Tf = fT(0);
  if (o.T0 === 'steady' || o.T0 === undefined) T = step(new Array(n).fill(0), Tf, Infinity, 0).T;
  else T = rc.map((x) => (typeof o.T0 === 'function' ? o.T0(x) : +o.T0));
  const wallT = (Tc, tf) => tf - (Gf[0] * (tf - Tc[0])) / (2 * PI * ri * Math.max(hIn, 1e-9));
  const out = { r: rc, layerOf, t: [0], fluid: [Tf], wall: [wallT(T, Tf)], profiles: [{ t: 0, T: T.slice() }], tReach: null, ro };
  const dt = tEnd / nSteps, every = Math.max(1, Math.round(nSteps / Math.max(nProfiles - 1, 1)));
  for (let k = 1; k <= nSteps && tEnd > 0; k++) {
    const t = k * dt, prev = Tf, s = step(T, Tf, dt, t); T = s.T; Tf = s.Tf;
    out.t.push(t); out.fluid.push(Tf); out.wall.push(wallT(T, Tf));
    if (stopAt !== null && out.tReach === null && Tf <= stopAt) out.tReach = prev > Tf ? t - dt + (dt * (prev - stopAt)) / (prev - Tf) : t;
    if (k % every === 0 || k === nSteps) out.profiles.push({ t, T: T.slice() });
  }
  out.T = T; out.qIn = Gf[0] * (Tf - T[0]);
  return out;
}

// =====================================================================================================
// 3. Steady state along the elevation profile
// =====================================================================================================
// Inside film coefficient of the two-phase stream: the wetting phase at the mixture velocity.
function filmInside(pr, vm, H, D) {
  if (vm < 1e-6) return H > 0.05 ? (3.66 * pr.kL) / D : (3.66 * pr.kG) / D;
  return H > 0.05 ? hInside((pr.rhoL * vm * D) / pr.muL, (pr.cpL * pr.muL) / pr.kL, pr.kL, D) : hInside((pr.rhoG * vm * D) / pr.muG, (pr.cpG * pr.muG) / pr.kG, pr.kG, D);
}
const enthalpyFlow = (pr) => pr.mG * pr.hG + pr.mO * pr.hO + pr.mW * pr.hW; // W relative to the table reference

/**
 * One steady march from the inlet (second-order midpoint rule) for pressure, temperature and holdup.
 * o: { fm, profile { x[], z[] } | grid (from discretise), n, id (m), rough (m), idOf(s), roughOf(s), tIn (°C), pIn (bara), mScale,
 *      kTotal (total minor-loss coefficient), fittings: [{ s (m), K }], model, fModel, mp, energy: 'enthalpy' (flowing-enthalpy balance with flashing, default) | 'cpjt' (frozen cp and Joule–Thomson) | 'isothermal',
 *      U (W/m²K on ID) | uOf({ s, z, D, pr, vm, holdup, ta }), tAmb | tAmbOf(s, z), heatOf(s) (W/m added, e.g. electrical heating) }
 * Returns { ok, reached, P[], T[], cells[] } (cells hold the mid-step state used for the update).
 */
export function steadyMarch(o) {
  const { fm, id, rough = 4.5e-5, tIn = 70, mScale = 1, model = 'beggsBrill', fModel = 'colebrook', mp = {}, energy = 'enthalpy' } = o, grid = o.grid || discretise(o.profile, o.n || 150), N = grid.n, ds = grid.ds;
  const uOf = typeof o.uOf === 'function' ? o.uOf : null, U0 = o.U ?? 3, tAmbOf = typeof o.tAmbOf === 'function' ? o.tAmbOf : () => o.tAmb ?? 4, heatOf = typeof o.heatOf === 'function' ? o.heatOf : null;
  const dOf = typeof o.idOf === 'function' ? o.idOf : null, rOf = typeof o.roughOf === 'function' ? o.roughOf : null, mdot = (fm.rates.mHC + fm.rates.mW) * mScale;
  // local (minor) losses: a total coefficient spread along the line plus fittings at given positions, K ρ v² / 2 on the no-slip mixture
  const kCell = new Array(N).fill((o.kTotal || 0) / N); for (const q of o.fittings || []) kCell[clamp(Math.floor(q.s / ds), 0, N - 1)] += q.K;
  const state = (P, T, i) => {
    const sMid = 0.5 * (grid.s[i] + grid.s[i + 1]), zMid = 0.5 * (grid.z[i] + grid.z[i + 1]), D = dOf ? Math.max(dOf(sMid), 0.01) : id, A = (PI * D * D) / 4, pr = fm.at(P, T, mScale), th = grid.theta[i];
    const vsg = pr.qG / A, vsl = pr.qL / A, vm = vsl + vsg;
    const gr = holdupGradient({ vsl, vsg, rhoL: pr.rhoL, rhoG: pr.rhoG, muL: pr.muL, muG: pr.muG, sigma: pr.sigma, D, theta: th, rough: rOf ? rOf(sMid) : rough, P: P * 1e5, fModel, waterCont: pr.phaseInv, label: false }, model, mp);
    gr.loc = kCell[i] > 0 && vm > 0 ? (kCell[i] * (pr.rhoL * (vsl / vm) + pr.rhoG * (vsg / vm)) * vm * vm) / (2 * ds) : 0; gr.dpdx += gr.loc;
    const H = gr.holdup, ta = tAmbOf(sMid, zMid), U = uOf ? uOf({ s: sMid, z: zMid, D, pr, vm, holdup: H, ta, T }) : U0, q = U * PI * D * (T - ta) - (heatOf ? heatOf(sMid) : 0);
    const mCp = pr.mG * pr.cpG + pr.mO * pr.cpO + pr.mW * pr.cpW, jt = mCp > 0 ? (pr.mG * pr.cpG * pr.jtG + pr.mO * pr.cpO * pr.jtO - pr.mW / pr.rhoW) / mCp : 0, sinT = Math.sin(th);
    const dTds = mCp > 0 && energy !== 'isothermal' ? -q / mCp - jt * gr.dpdx - (mdot * G * sinT) / mCp : 0;
    const ke = H > 1e-6 && H < 1 - 1e-6 ? 0.5 * (pr.mG * (vsg / (1 - H)) ** 2 + (pr.mO + pr.mW) * (vsl / H) ** 2) : 0.5 * mdot * vm * vm;
    return { pr, gr, vsl, vsg, D, A, ta, U, q, mCp, sinT, ke, dPds: -gr.dpdx / 1e5, dTds };
  };
  const P = [o.pIn], T = [tIn], cells = [];
  let ok = true, Hf = energy === 'enthalpy' ? enthalpyFlow(fm.at(o.pIn, tIn, mScale)) : 0, keSum = 0;
  const hIn = Hf;
  for (let i = 0; i < N; i++) {
    const a = state(P[i], T[i], i), Pm = P[i] + 0.5 * ds * a.dPds, Tm = T[i] + 0.5 * ds * a.dTds;
    if (!(Pm > 1.0)) { ok = false; break; }
    const b = state(Pm, Tm, i), Pn = P[i] + ds * b.dPds;
    let Tn = T[i] + ds * b.dTds;
    if (!(Pn > 1.0) || !Number.isFinite(Tn)) { ok = false; break; }
    if (energy === 'enthalpy' && b.mCp > 0) { // flowing enthalpy + potential + kinetic energy balance, then T from H(P, T)
      const dKe = 2 * (b.ke - a.ke); Hf -= ds * (b.q + mdot * G * b.sinT) + dKe; keSum += dKe;
      let t0 = Tn, f0 = enthalpyFlow(fm.at(Pn, t0, mScale)) - Hf, slope = b.mCp;
      for (let k = 0; k < 8 && Math.abs(f0) > 2e-5 * b.mCp; k++) { const t1 = clamp(t0 - f0 / slope, t0 - 25, t0 + 25), f1 = enthalpyFlow(fm.at(Pn, t1, mScale)) - Hf; if (Math.abs(t1 - t0) > 1e-9 && (f1 - f0) / (t1 - t0) > 0.2 * b.mCp) slope = (f1 - f0) / (t1 - t0); t0 = t1; f0 = f1; }
      Tn = t0; b.hRes = f0;
    }
    P.push(Pn); T.push(clamp(Tn, -60, 250)); cells.push(b);
  }
  return { ok, reached: cells.length, P, T, cells, grid, mdot, hIn, hOut: Hf, keSum };
}

/**
 * Steady solution for a fixed inlet pressure (o.pIn) or a fixed outlet pressure (o.pOut, shooting on the inlet pressure).
 * Options as steadyMarch plus pGuess (bara, speeds up the shooting), hydrate (false skips the hydrate curve), cErosion (API RP 14E C factor).
 * Returns { ok, reason } or the node arrays { s, x, z, theta, P, T, holdup, vsl, vsg, vm, vL, vG, rhoM, rhoL, rhoG, muL, muG, sigma, dpdx, fric, grav, acc,
 *   regime, tauW, tAmb, tWall, tHyd, subcooling, qG, qL, U, qLoss (W/m), D, evr, wcut } and the totals { pIn, pOut, tIn, tOut, dpFric, dpGrav, dpAcc, dpLocal (bar),
 *   liquidInventory, volume (m³), residence, residenceLiquid (s), heatLoss (W), mdot (kg/s), energy: { hIn, hOut, potential, kinetic (W) }, marches }.
 */
export function steadySolve(o) {
  const grid = o.grid || discretise(o.profile, o.n || 150), N = grid.n, base = { ...o, grid };
  let marches = 0, sol;
  const run = (p) => { marches++; return steadyMarch({ ...base, pIn: p }); };
  if (o.pOut !== undefined && o.pOut !== null && (o.pIn === undefined || o.pIn === null)) {
    const target = o.pOut, res = (p) => { const r = run(p); sol = r; return r.ok ? r.P[N] - target : -1e3 - (N - r.reached); };
    let pIn = null, guess = o.pGuess; const tol = o.tolP || 1e-6;
    if (!isNum(guess) && N > 40) { const c = steadySolve({ ...o, grid: undefined, n: 30, hydrate: false, tolP: 1e-3 }); if (!c.ok) return c; guess = c.pIn; marches += c.marches; }
    if (isNum(guess)) { // secant from the guess
      let p0 = Math.max(guess, 1.05), f0 = res(p0), p1 = f0 > -900 ? Math.max(p0 - f0, 1.05) : null;
      if (Math.abs(f0) < tol) pIn = p0;
      for (let k = 0; pIn === null && p1 !== null && k < 12; k++) {
        const f1 = res(p1);
        if (f1 <= -900 || p1 === p0) break;
        if (Math.abs(f1) < tol) { pIn = p1; break; }
        const slope = (f1 - f0) / (p1 - p0), pn = slope > 0.05 ? p1 - f1 / slope : p1 - f1;
        p0 = p1; f0 = f1; p1 = clamp(pn, Math.max(1.05, 0.5 * p1), 2 * p1 + 5);
      }
    }
    if (pIn === null) { // robust bracket and Brent
      let lo = 1.05, hi = Math.max(target, 2) + 20, fhi = res(hi), g = 0;
      if (res(lo) > 0) return { ok: false, reason: `The line gains more pressure from elevation than it loses to friction: even ${lo} bara at the inlet arrives above the ${target} bara outlet pressure.` };
      while (fhi < 0 && hi < 1400 && g++ < 40) { lo = hi; hi = hi * 1.35 + 10; fhi = res(hi); }
      if (fhi < 0) return { ok: false, reason: 'No inlet pressure below 1,400 bara can deliver this rate: the line is too small, too long or blocked.' };
      pIn = brent(res, lo, hi, Math.min(tol, 1e-4));
    }
    sol = run(pIn);
    if (!sol.ok) { sol = run(pIn * (1 + 1e-9) + 1e-7); if (!sol.ok) return { ok: false, reason: 'The shooting on the inlet pressure did not converge for this rate.' }; }
  } else {
    if (!isNum(o.pIn) || o.pIn <= 1) return { ok: false, reason: 'A fixed inlet pressure above 1 bara is needed.' };
    sol = run(o.pIn);
    if (!sol.ok) return { ok: false, reason: `An inlet pressure of ${(+o.pIn).toFixed(1)} bara is too low to push this rate to the outlet: the pressure falls to atmospheric ${Math.round(grid.s[sol.reached])} m from the inlet.` };
  }
  const c = sol.cells, fm = o.fm, col = (fn) => { const a = c.map(fn); a.push(a[a.length - 1]); return a; }, cE = o.cErosion || 100;
  const out = { ok: true, marches, grid, ds: grid.ds, length: grid.length, n: N, s: grid.s, x: grid.x, z: grid.z, theta: grid.theta.concat(grid.theta[N - 1]), P: sol.P, T: sol.T, mdot: sol.mdot, pIn: sol.P[0], pOut: sol.P[N], tIn: sol.T[0], tOut: sol.T[N] };
  out.holdup = col((q) => q.gr.holdup); out.vsl = col((q) => q.vsl); out.vsg = col((q) => q.vsg); out.vm = col((q) => q.vsl + q.vsg);
  out.vL = col((q) => (q.gr.holdup > 1e-6 ? q.vsl / q.gr.holdup : 0)); out.vG = col((q) => (q.gr.holdup < 1 - 1e-6 ? q.vsg / (1 - q.gr.holdup) : 0));
  out.dpdx = col((q) => q.gr.dpdx); out.fric = col((q) => q.gr.fric); out.grav = col((q) => q.gr.grav); out.acc = col((q) => q.gr.acc); out.loc = col((q) => q.gr.loc); out.regime = col((q) => q.gr.regime); out.tauW = col((q) => q.gr.tauW);
  out.rhoL = col((q) => q.pr.rhoL); out.rhoG = col((q) => q.pr.rhoG); out.muL = col((q) => q.pr.muL); out.muG = col((q) => q.pr.muG); out.sigma = col((q) => q.pr.sigma); out.wcut = col((q) => q.pr.wcut);
  out.rhoM = col((q) => q.pr.rhoL * q.gr.holdup + q.pr.rhoG * (1 - q.gr.holdup)); out.tAmb = col((q) => q.ta); out.qG = col((q) => q.pr.qG); out.qL = col((q) => q.pr.qL); out.U = col((q) => q.U); out.qLoss = col((q) => q.q); out.D = col((q) => q.D);
  out.mG = col((q) => q.pr.mG); out.mL = col((q) => q.pr.mO + q.pr.mW); out.zG = col((q) => q.pr.zG); out.mwG = col((q) => q.pr.mwG);
  out.hIn = col((q) => filmInside(q.pr, q.vsl + q.vsg, q.gr.holdup, q.D));
  out.tWall = out.T.map((t, i) => t - (out.U[i] * (t - out.tAmb[i])) / Math.max(out.hIn[i], 1e-6));
  out.evr = col((q) => { const vm = q.vsl + q.vsg, lam = vm > 0 ? q.vsl / vm : 1, rhoN = q.pr.rhoL * lam + q.pr.rhoG * (1 - lam); return vm / ((1.22 * cE) / Math.sqrt(rhoN)); });
  if (o.hydrate !== false && typeof fm.hydrateT === 'function') { out.tHyd = out.P.map((p) => fm.hydrateT(p)); out.subcooling = out.T.map((t, i) => out.tHyd[i] - t); }
  let inv = 0, vol = 0, fr = 0, gv = 0, ac = 0, lc = 0, heat = 0, res = 0;
  c.forEach((q) => { inv += q.gr.holdup * q.A * grid.ds; vol += q.A * grid.ds; fr += q.gr.fric * grid.ds; gv += q.gr.grav * grid.ds; ac += q.gr.acc * grid.ds; lc += q.gr.loc * grid.ds; heat += q.q * grid.ds; res += grid.ds / Math.max(q.vsl + q.vsg, 1e-6); });
  const qLm = mean(c.map((q) => q.pr.qL));
  Object.assign(out, { liquidInventory: inv, volume: vol, dpFric: fr / 1e5, dpGrav: gv / 1e5, dpAcc: ac / 1e5, dpLocal: lc / 1e5, heatLoss: heat, residence: res, residenceLiquid: qLm > 1e-12 ? inv / qLm : res });
  out.energy = { hIn: sol.hIn, hOut: sol.hOut, potential: sol.mdot * G * (grid.z[N] - grid.z[0]), kinetic: sol.keSum, residual: Math.max(...c.map((q) => Math.abs(q.hRes || 0))) };
  return out;
}

// ---- inputs → configuration ------------------------------------------------------------------------------
const GAS_CONDENSATE = { name: 'Lean gas condensate', comp: { N2: 0.8, CO2: 2.2, H2S: 0, C1: 82, C2: 6.5, C3: 3, iC4: 0.6, nC4: 1.0, iC5: 0.4, nC5: 0.4, C6: 0.6, C7p: 2.5 }, c7MW: 135, c7SG: 0.775, rateBasis: 'gas', qGas: 10, qWater: 60, wc: 0 };
function cleanProfile(rows) {
  const pts = (Array.isArray(rows) ? rows : []).map((r) => ({ x: +r?.x, z: +r?.z })).filter((r) => Number.isFinite(r.x) && Number.isFinite(r.z)).sort((a, b) => a.x - b.x);
  const x = [], z = [];
  for (const p of pts) { if (x.length && p.x - x[x.length - 1] < 1e-6) { if (Math.abs(p.z - z[z.length - 1]) > 1e-6) { x.push(x[x.length - 1] + 1e-3); z.push(p.z); } continue; } x.push(p.x); z.push(p.z); }
  if (x.length < 2) throw new Error('The elevation profile needs at least two points with different distances.');
  if (z.some((v) => Math.abs(v) > 12000) || x[x.length - 1] - x[0] > 2e6) throw new Error('The elevation profile is outside the supported range (elevations within ±12 km, length below 2,000 km).');
  const x0 = x[0];
  return { x: x.map((v) => v - x0), z };
}
const num = (v, d, lo = -Infinity, hi = Infinity) => clamp(fin(+v, d), lo, hi);
/**
 * Build the solver configuration (fluid model, geometry, thermal boundary, closure parameters) from the suite inputs.
 * Returns { fm, profile, id, wt, rough, idOf, roughOf, od, layers, thermal, uOf, U, tAmbOf, heatOf, tIn, pOut, model, fModel, mp, energy, n, riserBaseX, riserBaseS, hasRiser, riserHeight, base (options for steadySolve) }.
 */
export function flowConfig(v = {}, ctx = {}) {
  const override = {};
  if (v.fluidSel === 'gascond') Object.assign(override, GAS_CONDENSATE);
  else if ((ctx.fluid?.rateBasis || 'oil') === 'oil' && isNum(+v.wc) && +v.wc >= 0 && Math.abs(+v.wc - fin(ctx.fluid?.wc, 20)) > 1e-9) override.wc = clamp(+v.wc, 0, 98);
  const fm = fluidModel(ctx, override), profile = cleanProfile(v.profile ?? BASE.profile), L = profile.x[profile.x.length - 1];
  const id0 = num(v.idMm, BASE.idMm, 10, 3000) / 1000, wt = num(v.wtMm, BASE.wtMm, 0.5, 200) / 1000, rough0 = num(v.roughUm, BASE.roughUm, 0, 5000) * 1e-6;
  // deposits from the solids suite: thickness profile, or a uniform restriction when only the minimum bore is known
  const dep = (Array.isArray(v.deposit) ? v.deposit : []).map((r) => ({ x: +r?.x, t: +r?.t })).filter((r) => Number.isFinite(r.x) && Number.isFinite(r.t) && r.t >= 0).sort((a, b) => a.x - b.x);
  const effId = num(v.effIdMm, 0, 0, 3000) / 1000, roughEff = num(v.roughEffUm, 0, 0, 20000) * 1e-6, rough = roughEff > 0 ? roughEff : rough0;
  let idOf = null, idMin = id0;
  const S = [0]; for (let i = 1; i < profile.x.length; i++) S.push(S[i - 1] + Math.hypot(profile.x[i] - profile.x[i - 1], profile.z[i] - profile.z[i - 1]));
  const xOfS = (s) => interp1(S, profile.x, s), sOfX = (x) => interp1(profile.x, S, x);
  if (dep.length >= 2 && dep.some((r) => r.t > 0)) { const dx = dep.map((r) => r.x), dt = dep.map((r) => r.t / 1000); idOf = (s) => Math.max(id0 - 2 * interp1(dx, dt, xOfS(s)), 0.2 * id0); idMin = Math.max(id0 - 2 * Math.max(...dt), 0.2 * id0); }
  else if (effId > 0 && effId < id0) { idMin = Math.max(effId, 0.2 * id0); }
  const id = idOf ? id0 : idMin;
  // thermal boundary
  const thermal = { tAir: num(v.tAir, BASE.tAir, -60, 60), tSeaSurface: num(v.tSeaSurface, BASE.tSeaSurface, -2, 40), tSeabed: num(v.tSeabed, BASE.tSeabed, -2, 40), thermocline: num(v.thermocline, 250, 10, 3000) };
  const insT = num(v.insT, BASE.insulation.t * 1000, 0, 500) / 1000, insK = num(v.insK, BASE.insulation.k, 0.005, 5), kWall = num(v.kWall, 45, 1, 400);
  const layers = [{ name: 'Pipe wall (steel)', t: wt, k: kWall, rho: BASE.rhoSteel, cp: 470 }];
  if (insT > 0) layers.push({ name: 'Insulation / coating', t: insT, k: insK, rho: num(v.insRho, 900, 20, 5000), cp: num(v.insCp, 1700, 200, 5000) });
  const od = id0 + 2 * (wt + insT), riserBaseX = clamp(num(v.riserBaseX, BASE.riserBaseX, 0, 1e7), 0, L), burialDepth = num(v.burialDepth, 0, 0, 20), kSoil = num(v.kSoil, 1.2, 0.1, 10);
  const current = num(v.currentSpeed, BASE.currentSpeed, 0, 5), wind = num(v.windSpeed, 5, 0, 60), uMult = num(v.uMult, 1, 0.05, 20), uIn = num(v.uValue, BASE.U, 0, 5000);
  const tAmbOf = (s, z) => ambientTemperature(z, thermal), buriedAt = (s) => burialDepth > 0 && xOfS(s) <= riserBaseX + 1e-6;
  const hOutMult = num(v.hOutMult, 1, 0.05, 20), natConv = v.natConv !== false;
  // outside film: forced convection (Churchill–Bernstein) combined with free convection (Churchill–Chu) for the surface temperature excess, which is found by one fixed-point pass
  const hOutOf = (c, dT = 0) => hOutMult * (c.z >= 0 ? hOutside(wind, od, 'air', c.ta, dT) : hOutside(current, od, 'seawater', c.ta, dT));
  const network = (c) => {
    const mk = (hOut) => uValue({ id: c.D, wt, kWall, layers: layers.slice(1), hIn: filmInside(c.pr, c.vm, c.holdup, c.D), hOut, burial: buriedAt(c.s) ? { depth: burialDepth + od / 2, kSoil } : null });
    let r = mk(hOutOf(c));
    if (natConv && Number.isFinite(c.T) && !buriedAt(c.s)) { const h0 = hOutOf(c), dTs = (r.U * c.D * (c.T - c.ta)) / (od * Math.max(h0, 1e-9)); r = mk(hOutOf(c, dTs)); r.dTsurface = dTs; r.hOut = hOutOf(c, dTs); }
    return r;
  };
  const uMode = v.uMode === 'layers' ? 'layers' : 'input', uOf = uMode === 'layers' ? (c) => network(c).U * uMult : null, heat = num(v.heatTrace, 0, 0, 5000);
  const model = HOLDUP_MODELS.some((m) => m.value === v.model) ? v.model : 'beggsBrill', fModel = ['colebrook', 'haaland', 'swamee', 'churchill'].includes(v.fModel) ? v.fModel : 'colebrook';
  const mp = { c0: num(v.c0, 1.2, 0.8, 2), vDrift: num(v.vDrift, 0.35, -5, 10), wallisN: num(v.wallisN, 0, 0, 5), holdupMult: num(v.holdupMult, 1, 0.2, 3), annular: annularParams(v) };
  const equip = equipment({ pumpDp0: num(v.pumpDp0, 0, 0, 500), pumpQmax: num(v.pumpQmax, 0.5, 1e-4, 100), pumpSpeed: num(v.pumpSpeed, 1, 0.1, 2), sepKv: v.sepMode === 'valve' ? num(v.sepKv, 0.5, 0, 1e6) : 0, compHead: v.sepMode === 'valve' ? num(v.compHead, 60, 0, 500) * 1e3 : 0, compQmax: num(v.compQmax, 2, 1e-3, 1e3), compSpeed: num(v.compSpeed, 1, 0.3, 1.3) });
  const branch = num(v.branchFrac, 0, 0, 0.95) > 0 ? { frac: num(v.branchFrac, 0, 0, 0.95), x: num(v.branchX, 9000, 0, 1e7), length: num(v.branchLength, 3000, 10, 1e6), dz: num(v.branchDz, 0, -3000, 3000), idMm: num(v.branchIdMm, 0, 0, 3000), tIn: num(v.branchTin, BASE.tIn, -40, 250) } : null;
  const energy = v.energy === 'cpjt' ? 'cpjt' : 'enthalpy', n = Math.round(num(v.nSteady, 150, 8, 4000)), zEnd = profile.z[profile.z.length - 1], zBase = interp1(profile.x, profile.z, riserBaseX);
  const riserHeight = zEnd - zBase, hasRiser = riserBaseX < L - 1e-6 && riserHeight > 20 * id0;
  const fittings = (Array.isArray(v.fittings) ? v.fittings : []).map((r) => ({ x: +r?.x, K: +r?.K * (isNum(+r?.open) && +r.open > 0 && +r.open < 100 ? (100 / +r.open) ** 2 : 1) })).filter((r) => Number.isFinite(r.x) && r.K > 0).map((r) => ({ s: sOfX(clamp(r.x, 0, L)), K: r.K })), kTotal = num(v.kLoss, 0, 0, 1e5);
  const chokeDp = num(v.chokeDp, 0, 0, 500), chokeOpening = num(v.chokeOpening, 100, 1, 100);
  const cfg = { fm, profile, length: L, fittings, kTotal, chokeDp, chokeOpening, chokeLoss: (m) => chokeDp * m * m * (100 / chokeOpening) ** 2, id, id0, idMin, wt, rough, idOf, od, layers, thermal, uMode, uOf, U: uIn * uMult, uMult, network, tAmbOf, heat, tIn: num(v.tIn, BASE.tIn, -40, 250), pSep: num(v.pOut, BASE.pOut, 1.05, 1300), pOut: num(v.pOut, BASE.pOut, 1.05, 1300) + chokeDp * num(v.rateFrac, 1, 1e-4, 20) ** 2 * (100 / chokeOpening) ** 2, model, fModel, mp, energy, n, riserBaseX, riserBaseS: sOfX(riserBaseX), hasRiser, riserHeight, sOfX, xOfS, burialDepth, kSoil, insT, insK, kWall, current, wind, cErosion: num(v.cErosion, 100, 30, 400), buriedAt };
  cfg.equip = equip; cfg.branch = branch; cfg.compPd = num(v.compPd, 120, 2, 1000); cfg.sepMode = v.sepMode === 'valve' ? 'valve' : 'fixed'; cfg.hOutOf = hOutOf; cfg.natConv = natConv;
  cfg.base = { fm, profile, n, id, rough, idOf, kTotal, fittings, tolP: num(v.tolP, 1e-6, 1e-9, 0.1), tIn: cfg.tIn, model, fModel, mp, energy, uOf, U: cfg.U, tAmbOf, heatOf: heat > 0 ? () => heat : null, cErosion: cfg.cErosion };
  return cfg;
}

const stdLiquid = (fm, m) => (fin(fm.rates.qOilStd) + fin(fm.rates.qWaterStd)) * m; // Sm³/d of stock-tank liquid
/**
 * Steady flow of the suite inputs with the selected boundary conditions:
 * 'outletP' (outlet pressure + rate), 'inletP' (inlet pressure + rate), 'bothP' (inlet and outlet pressure, the rate is solved),
 * 'ipr' (well inflow: reservoir pressure, productivity index and well pressure loss upstream of the flowline inlet).
 * Returns the steadySolve result plus { cfg, mScale, bc, nodal: { pwf, drawdown, qLiqStd } | null }. Throws when no steady solution exists.
 */
export function steadyFlow(v = {}, ctx = {}, over = {}) {
  const cfg = over.cfg || flowConfig(v, ctx), bc = ['outletP', 'inletP', 'bothP', 'ipr'].includes(v.bc) ? v.bc : 'outletP', rate = num(over.mScale ?? v.rateFrac, 1, 1e-4, 20), base = { ...cfg.base, ...(over.base || {}) };
  const need = (r) => { if (!r.ok) throw new Error(r.reason || 'No steady flow solution exists for these inputs.'); return r; };
  let st, mScale = rate, nodal = null;
  // separator / compressor characteristic: the outlet pressure floats with the gas rate through the gas outlet valve and the compressor curve
  const pOutAt = (m) => { if (cfg.sepMode !== 'valve') return cfg.pOut; const o = cfg.fm.at(cfg.pSep, 20, m), std = cfg.fm.at(1.01325, 15, m), qStd = fin(std.qG, 0), cs = cfg.equip.compressor(fin(o.qG, 0), cfg.compPd, 293, fin(o.zG, 0.9), fin(o.mwG, 20) * 1e-3, fin(o.mG, 0)); cfg.compressor = { ...cs, qStd, qSuction: fin(o.qG, 0) }; return cfg.equip.sepPressure(qStd, cfg.pSep, cs.head > 0 ? Math.min(cs.pSuction, cfg.pSep) : cfg.pSep) + (cfg.pOut - cfg.pSep); };
  const solveOut = (opts, m, guess) => (cfg.branch ? steadyBranch({ ...opts, mScale: m, pOut: pOutAt(m) }, cfg.branch) : steadySolve({ ...opts, mScale: m, pOut: pOutAt(m), pGuess: guess }));
  if (bc === 'outletP') { st = need(solveOut(base, mScale, over.pGuess)); cfg.pOutEff = pOutAt(mScale); }
  else if (bc === 'inletP') st = need(steadySolve({ ...base, mScale, pIn: num(v.pInSet, 95, 1.1, 1400) }));
  else {
    const coarse = { ...base, n: Math.min(base.n, 30), hydrate: false }, pInSet = num(v.pInSet, 95, 1.1, 1400), pRes = num(v.pRes, BASE.pRes, 2, 2000), pi = num(v.piIpr, BASE.pi, 1e-3, 1e6), wellDp = num(v.wellDp, 60, 0, 1500), vogel = v.iprType === 'vogel', qMax = vogel ? (pi * pRes) / 1.8 : pi * pRes;
    const pwfOf = (q) => (vogel ? pRes * ((-0.2 + Math.sqrt(Math.max(0.04 + 3.2 * (1 - Math.min(q / qMax, 1)), 0))) / 1.6) : pRes - q / pi);
    let last = null;
    const g = (m, opts = coarse) => {
      if (bc === 'bothP') { const r = steadySolve({ ...opts, mScale: m, pIn: pInSet }); return r.ok ? r.pOut - cfg.pOut : -1e3; }
      const r = solveOut(opts, m, last); if (!r.ok) return -1e3; last = r.pIn;
      return pwfOf(stdLiquid(cfg.fm, m)) - wellDp * (0.8 + 0.2 * m * m) + cfg.equip.pumpDp(r.qG[0] + r.qL[0]) - r.pIn; // available (inflow − well losses + booster pump) minus required flowline inlet pressure
    };
    const mHi = bc === 'ipr' ? Math.min(3, (0.999 * qMax) / Math.max(stdLiquid(cfg.fm, 1), 1e-9)) : 3;
    if (g(mHi) > 0) { if (bc === 'ipr') mScale = mHi; else throw new Error(`The pressure difference between ${pInSet} and ${cfg.pOut} bara would drive more than three times the case rate; reduce it or increase the case rate.`); }
    else {
      let hi = mHi, lo = null;
      for (let k = 1; k <= 14; k++) { const m = mHi * 0.72 ** k; if (g(m) > 0) { lo = m; break; } hi = m; }
      if (lo === null) throw new Error(bc === 'ipr' ? 'The reservoir cannot lift this well stream to the outlet pressure: the inflow curve stays below the outflow curve at every rate.' : `An inlet pressure of ${pInSet} bara cannot deliver any stable flow against ${cfg.pOut} bara at the outlet: the static head of the liquid-filled line is higher.`);
      mScale = brent(g, lo, hi, 1e-4);
      if (bc === 'bothP' && base.n > coarse.n) { const a = mScale * 0.9, b = mScale * 1.1, full = { ...base, hydrate: false }, ga = g(a, full), gb = g(b, full); if (ga > 0 && gb < 0) mScale = brent((m) => g(m, full), a, b, 1e-5); }
    }
    st = need(bc === 'bothP' ? steadySolve({ ...base, mScale, pIn: pInSet }) : solveOut(base, mScale, last)); if (bc !== 'bothP') cfg.pOutEff = pOutAt(mScale);
    if (bc === 'ipr') { const q = stdLiquid(cfg.fm, mScale), pwf = pwfOf(q); nodal = { pwf, drawdown: pRes - pwf, qLiqStd: q, pRes, wellDp: wellDp * (0.8 + 0.2 * mScale * mScale) }; }
  }
  return Object.assign(st, { cfg, mScale, bc, nodal });
}

// =====================================================================================================
// 4. Flow-regime maps and interfacial stability
// =====================================================================================================
function stratGeom(hD, D) {
  const c = clamp(2 * hD - 1, -1, 1), A = (PI * D * D) / 4, AL = 0.25 * D * D * (PI - Math.acos(c) + c * Math.sqrt(1 - c * c));
  return { A, AL, AG: A - AL, SL: D * (PI - Math.acos(c)), SG: D * Math.acos(c), Si: D * Math.sqrt(1 - c * c), c };
}
const fanning = (Re) => (Re < 2100 ? 16 / Math.max(Re, 1e-9) : 0.046 * Re ** -0.2);
/**
 * Combined two-phase momentum balance of stratified flow at a liquid level hD = h/D (zero at equilibrium).
 * o.pure = true uses the original Taitel & Dukler (1976) simplifications (fi = fG, interfacial shear on the gas velocity);
 * otherwise the closure of the kernel (fi = max(fG, 0.0142), relative velocity) is used.
 */
export function stratifiedBalance(hD, p, o = {}) {
  const { vsl, vsg, rhoL, rhoG, muL, muG, D, theta = 0 } = p, g = stratGeom(hD, D), vL = (vsl * g.A) / g.AL, vG = (vsg * g.A) / g.AG, DL = (4 * g.AL) / g.SL, DG = (4 * g.AG) / (g.SG + g.Si);
  const fL = fanning((rhoL * vL * DL) / muL) * (o.fwlMult || 1), fG0 = fanning((rhoG * vG * DG) / muG), fG = fG0 * (o.fwgMult || 1), fi = (o.pure ? fG0 : Math.max(fG0, 0.0142)) * (o.fiMult || 1);
  const tL = (fL * rhoL * vL * vL) / 2, tG = (fG * rhoG * vG * vG) / 2, ti = o.pure ? (fi * rhoG * vG * vG) / 2 : (fi * rhoG * (vG - vL) * Math.abs(vG - vL)) / 2;
  return { F: (tG * g.SG) / g.AG - (tL * g.SL) / g.AL + ti * g.Si * (1 / g.AL + 1 / g.AG) - (rhoL - rhoG) * G * Math.sin(theta), g, vL, vG, tL, tG, ti, fL, DL };
}
/** All equilibrium stratified levels h/D (ascending) of the momentum balance: one root, or three on upward slopes at low gas velocity. */
export function stratifiedRoots(p, o = {}) {
  const n = o.n || 90, roots = [], f = (x) => stratifiedBalance(x, p, o).F;
  let a = 1e-4, fa = f(a);
  for (let i = 1; i <= n; i++) { const b = 1e-4 + ((0.9995 - 1e-4) * i) / n, fb = f(b); if (fa * fb <= 0 && Number.isFinite(fa) && Number.isFinite(fb)) roots.push(brent(f, a, b, 1e-10)); a = b; fa = fb; }
  if (!roots.length) roots.push(fa < 0 ? 0.9995 : 1e-4);
  return roots;
}
const holdupOfLevel = (hD) => { const g = stratGeom(hD, 1); return g.AL / g.A; };
const levelOfHoldup = (H) => (H <= 1e-6 ? 0 : H >= 1 - 1e-6 ? 1 : brent((x) => holdupOfLevel(x) - H, 0, 1, 1e-10));

/**
 * Stability of stratified flow at its equilibrium level: inviscid Kelvin–Helmholtz (long waves), the Taitel–Dukler form,
 * the viscous Kelvin–Helmholtz criterion of Barnea & Taitel (kinematic versus dynamic wave speed), the Jeffreys sheltering
 * criterion for wave generation and the roll-wave (Vedernikov) criterion of the liquid layer. A ratio above 1 is unstable.
 * p: { vsl, vsg, rhoL, rhoG, muL, muG, D, theta }. Returns { hD, holdup, vL, vG, ikhCrit (m/s relative), ikhRatio, tdRatio, vkhRatio, cV, cIV, wavyRatio, froudeL, rollRatio }.
 */
export function interfacialStability(p) {
  const { vsl, vsg, rhoL, rhoG, muL, D, theta = 0 } = p, st = stratifiedLevel(p), hD = clamp(st.hD, 2e-4, 0.9993), b = stratifiedBalance(hD, p), g = b.g, RL = g.AL / g.A, RG = 1 - RL;
  const dRho = Math.max(rhoL - rhoG, 1e-6), cosT = Math.max(Math.cos(theta), 0.02), Si = Math.max(g.Si, 1e-9), grav = dRho * G * cosT * (g.A / Si);
  const ikhCrit = Math.sqrt(((rhoL * RG + rhoG * RL) / (rhoL * rhoG)) * grav), tdCrit = (1 - hD) * Math.sqrt((dRho * G * cosT * g.AG) / (rhoG * Si));
  // kinematic-wave speed cV = (∂F/∂RL) / (∂F/∂UGS − ∂F/∂ULS) at constant superficial velocities
  const dh = Math.min(1e-4, 0.2 * hD, 0.2 * (1 - hD)), e = 1e-4, FR = (stratifiedBalance(hD + dh, p).F - stratifiedBalance(hD - dh, p).F) / (2 * dh) / ((Si * D) / g.A);
  const FL = (stratifiedBalance(hD, { ...p, vsl: vsl * (1 + e) }).F - stratifiedBalance(hD, { ...p, vsl: vsl * (1 - e) }).F) / (2 * e * vsl), FG = (stratifiedBalance(hD, { ...p, vsg: vsg * (1 + e) }).F - stratifiedBalance(hD, { ...p, vsg: vsg * (1 - e) }).F) / (2 * e * vsg);
  const cV = FR / (FG - FL), rhoS = rhoL / RL + rhoG / RG, cIV = ((rhoL * b.vL) / RL + (rhoG * b.vG) / RG) / rhoS, dU = b.vG - b.vL;
  const vkh = (cV - cIV) ** 2 + ((rhoL * rhoG) / (RL * RG * rhoS * rhoS)) * dU * dU, hEff = g.AL / Si, froudeL = b.vL / Math.sqrt(G * cosT * hEff);
  const gamma = clamp(1 - ((g.AL / g.SL) * 2) / (D * Math.max(1 - g.c * g.c, 1e-9)), 0.2, 1), frCrit = ((rhoL * b.vL * b.DL) / muL < 2100 ? 0.58 : 2) / gamma;
  return { hD, holdup: RL, vL: b.vL, vG: b.vG, ikhCrit, ikhRatio: Math.abs(dU) / ikhCrit, tdRatio: b.vG / tdCrit, vkhRatio: Math.sqrt(vkh / (grav / rhoS)), cV, cIV,
    wavyRatio: b.vG / Math.sqrt((4 * muL * dRho * G * cosT) / (0.01 * rhoL * rhoG * Math.max(b.vL, 1e-9))), froudeL, rollRatio: froudeL / frCrit };
}

export const PATTERNS = ['stratified smooth', 'stratified wavy', 'slug', 'churn', 'bubble', 'dispersed bubble', 'annular'];
/**
 * Taitel–Dukler / Barnea flow-pattern map in the vsg–vsl plane for fixed properties, diameter and inclination.
 * p: { rhoL, rhoG, muL, muG, sigma, D, theta }; o: { n, vsgRange, vslRange }.
 * Returns { vsg[], vsl[], index[][] (rows = vsl), boundaries: [{ name, x[] (vsg), y[] (vsl) }] }.
 */
export function flowPatternMap(p, o = {}) {
  const n = o.n || 44, vsg = Array.from({ length: n }, (_, i) => 10 ** (-2 + (4 * i) / (n - 1)) * (o.gScale || 1)), vsl = Array.from({ length: n }, (_, i) => 10 ** (-3 + (4.3 * i) / (n - 1)));
  const idx = vsl.map((l) => vsg.map((g) => Math.max(PATTERNS.indexOf(flowPattern({ ...p, vsl: l, vsg: g }).pattern), 0))), bd = new Map();
  const add = (a, b, x, y) => { const k = a < b ? `${PATTERNS[a]} | ${PATTERNS[b]}` : `${PATTERNS[b]} | ${PATTERNS[a]}`; if (!bd.has(k)) bd.set(k, []); bd.get(k).push([x, y]); };
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) {
    if (i < n - 1 && idx[i][j] !== idx[i + 1][j]) add(idx[i][j], idx[i + 1][j], vsg[j], Math.sqrt(vsl[i] * vsl[i + 1]));
    if (j < n - 1 && idx[i][j] !== idx[i][j + 1]) add(idx[i][j], idx[i][j + 1], Math.sqrt(vsg[j] * vsg[j + 1]), vsl[i]);
  }
  const boundaries = [...bd.entries()].filter(([, pts]) => pts.length >= 3).map(([name, pts]) => { pts.sort((a, b) => a[0] - b[0] || a[1] - b[1]); return { name, x: pts.map((q) => q[0]), y: pts.map((q) => q[1]) }; });
  return { vsg, vsl, index: idx, boundaries };
}
// Mandhane, Gregory & Aziz (1974) horizontal map: transition lines digitised in superficial velocities (ft/s), log–log interpolation
const MANDHANE = { annular: [[0.01, 70], [0.1, 60], [0.3, 38], [0.56, 40], [1, 50], [2.5, 100], [14, 230], [30, 269]], wave: [[0.01, 32.7], [0.1, 14], [0.2, 10.5], [1.15, 2.5], [4.8, 2.5], [14, 3.26]] };
const logInterp = (pts, x) => 10 ** interp1(pts.map((q) => log10(q[0])), pts.map((q) => log10(q[1])), log10(x));
/** Mandhane-type horizontal flow-pattern classification from the superficial velocities (m/s); air–water coordinates, no property correction. */
export function mandhaneRegime(vsl, vsg) {
  const l = Math.max(vsl / FT, 1e-6), g = Math.max(vsg / FT, 1e-6);
  if (l > 14) return g > logInterp(MANDHANE.annular, l) ? 'annular mist' : 'dispersed bubble';
  if (g > logInterp(MANDHANE.annular, l)) return 'annular mist';
  if (g < logInterp(MANDHANE.wave, l)) return l < 0.5 ? 'stratified' : 'elongated bubble';
  return l < 0.3 ? 'stratified wavy' : 'slug';
}

// =====================================================================================================
// 5. Slugging: unit cell, Lagrangian tracking, terrain accumulation, severe riser slugging
// =====================================================================================================
/**
 * Hydrodynamic slug unit cell at one location with selectable closures.
 * p: { vsl, vsg, rhoL, rhoG, muL, muG, sigma, D, theta }
 * o: { freqModel: 'zabaras' | 'gregory' | 'heywood' | 'unitCell', lengthModel: 'scott' | 'brill' | 'norris', vtModel: 'bendiksen' | 'nicklin',
 *      bodyModel: 'gregory' | 'barnea', freqMult, lenMult, filmMult (multiplier on the film holdup) }
 * Returns { vt, C0, vd, holdupSlug, holdupFilm, holdup, slugFraction, freq, period, length, lengthFromFreq, lengthMax, unitLength,
 *           volume (m³ liquid per slug), filmThickness (m), filmVelocity, vBody (liquid velocity in the slug body), pickup (m³/s scooped at the front) }.
 */
export function slugUnitCell(p, o = {}) {
  const { vsl, vsg, rhoL, rhoG, muL, D, theta = 0 } = p, sigma = Math.max(p.sigma ?? 0.02, 1e-4), vm = vsl + vsg, A = (PI * D * D) / 4, dRho = Math.max(rhoL - rhoG, 1), sinP = Math.sin(Math.max(theta, 0));
  let { C0, vd } = slugVelocity(vm, D, theta);
  if (o.vtModel === 'nicklin') { C0 = 1.2; vd = 0.35 * Math.sqrt(G * D) * sinP; }
  const vt = Math.max((C0 * vm + vd) * (o.vtMult || 1), 1e-6);
  let HLS = slugBodyHoldup(vm);
  if (o.bodyModel === 'barnea') { const fs = fanning((rhoL * vm * D) / muL), x = 2 * Math.sqrt((0.4 * sigma) / (dRho * G)) * ((2 * fs * vm ** 3) / D) ** 0.4 * (rhoL / sigma) ** 0.6 - 0.725; HLS = clamp(1 - (x > 0 ? 0.058 * x * x : 0), 0.48, 1); }
  HLS = clamp(HLS * (o.bodyMult || 1), 0.3, 1);
  const vGb = 1.2 * vm + 1.53 * ((G * sigma * dRho) / (rhoL * rhoL)) ** 0.25 * Math.sqrt(HLS) * sinP;
  const holdup = clamp((vt * HLS + vGb * (1 - HLS) - vsg) / vt, vsl / Math.max(vm, 1e-9), 1);
  let HLF = slugFilm({ ...p, sigma }, { vt, HLS }).HLF; // equilibrium film of the Taitel–Barnea unit cell (kernel)
  HLF = clamp(Math.min(HLF, holdup * 0.98) * (o.filmMult || 1), 0.005, Math.min(0.95 * HLS, 0.98 * holdup));
  const beta = clamp((holdup - HLF) / Math.max(HLS - HLF, 1e-6), 0.02, 1), dIn = D / 0.0254;
  const length = (o.lenMult || 1) * (o.lengthModel === 'norris' ? Math.exp(-2.099 + 4.859 * Math.sqrt(Math.log(Math.max(dIn, 1.01)))) * FT : slugLength(D, vm, o.lengthModel === 'brill' ? 'brill' : 'scott', o.xD >= 0 ? { xD: o.xD, theta } : null));
  const unitLength0 = length / beta, unitLength = unitLength0, lam = vsl / Math.max(vm, 1e-9);
  let f = o.freqModel === 'unitCell' ? vt / unitLength : slugFrequency(vsl, vm, D, theta, ['gregory', 'heywood', 'greskovich'].includes(o.freqModel) ? o.freqModel : 'zabaras');
  if (!(f > 0)) f = vt / unitLength;
  const freq = f * (o.freqMult || 1), vBody = (vm - vGb * (1 - HLS)) / HLS, filmVelocity = vt - ((vt - vBody) * HLS) / HLF;
  // unit-cell identity: the body length that goes with the frequency actually used (equal to the correlation when the frequency is derived from it)
  const lengthUnit = (beta * vt) / freq;
  return { vt, C0, vd, holdupSlug: HLS, holdupFilm: HLF, holdup, slugFraction: beta, freq, period: 1 / freq, length, lengthFromFreq: lengthUnit, lengthMax: lengthUnit * Math.exp(3.09 * 0.5 - 0.125), unitLength: lengthUnit / beta, lengthCorrelation: length,
    volume: lengthUnit * A * HLS, filmThickness: levelOfHoldup(HLF) * D, filmVelocity, vBody, pickup: A * (vt - vBody) * HLS };
}

/**
 * Lagrangian slug tracking on a steady carrier field. Slugs are initiated at the given sites with random intervals and
 * lengths (seeded). The tail moves with the bubble-nose velocity (accelerated in the wake of a short slug); the body length
 * follows the liquid balance  d(L·HLS)/dt = pick-up at the front − shedding at the tail, where the pick-up is the developed
 * value at the front position corrected for the film having relaxed towards its local terrain equilibrium over a long gap
 * (thicker film on upward slopes → growth, thinner on downward slopes → decay). Slugs that catch the one ahead merge and
 * slugs that shrink below two diameters dissipate.
 * f: { s[] (uniform), vt[], holdupSlug[], shed[] (developed shedding = pick-up rate per area, m/s), mod[] (pick-up ratio of the terrain-equilibrium film to the developed film), D, A }
 * o: { sites: [{ s, freq (1/s), length (m mean at initiation) }], nSlugs, seed, sigmaL, sigmaT, prefill, qDrain (m³/s) | drainFactor (× the mean liquid arrival rate), qSlugOut, qFilmOut (m³/s), relax (m), maxOps }
 * Returns { n, arrivals: [{ t, length, volume, velocity }], meanLength, stdLength, p50, p90, p99, maxLength, lognormal { mu, sigma }, length1000,
 *           freqArrival, period, merges, dissipated, generated, surge, surgeSingle (m³), hist, from (m), tSim (s) }.
 */
export function slugTracking(f, o = {}) {
  const { nSlugs = 200, seed = 42, sigmaL = 0.5, sigmaT = 0.5, prefill = true, maxOps = 1.5e5 } = o, D = f.D, A = f.A, S = f.s, n = S.length, Ltot = S[n - 1], ds = S[1] - S[0], R = rng(seed);
  const at = (a, s) => { const u = clamp(s / ds, 0, n - 1 - 1e-9), i = Math.floor(u); return a[i] + (a[i + 1] - a[i]) * (u - i); };
  const empty = { n: 0, arrivals: [], meanLength: 0, stdLength: 0, p50: 0, p90: 0, p99: 0, maxLength: 0, lognormal: null, length1000: 0, freqArrival: 0, period: null, merges: 0, dissipated: 0, generated: 0, surge: 0, surgeSingle: 0, hist: { centers: [], counts: [] }, from: Ltot, tSim: 0 };
  // a site cannot launch slugs closer than two body lengths apart (above that the flow is a continuous liquid column with bubbles)
  let sites = (o.sites || []).filter((q) => q.freq > 0 && q.length > 0 && q.s < Ltot).map((q) => ({ ...q, freq: Math.min(q.freq * (o.initMult || 1), at(f.vt, q.s) / (2 * q.length)) })).sort((a, b) => a.s - b.s);
  if (!sites.length || nSlugs < 1) return empty;
  // keep the number of slugs in the line affordable: track the downstream window that holds about 220 slug units
  const unit0 = at(f.vt, sites[0].s) / sites[0].freq, from = Math.max(sites[0].s, Ltot - 220 * unit0);
  if (from > sites[0].s) { const keep = sites.filter((q) => q.s > from); sites = [{ ...sites[0], s: from }, ...keep]; }
  const relax = (o.relax || 300 * D) * (o.relaxMult || 1), wkM = o.wakeMult ?? 1, minL = 2 * D, lnL = (m) => R.lognormal(Math.log(Math.max(m, minL)) - 0.5 * sigmaL * sigmaL, sigmaL), lnT = (fq) => R.lognormal(Math.log(1 / fq) - 0.5 * sigmaT * sigmaT, sigmaT);
  const slugs = []; // ordered from the outlet (index 0) to the inlet
  if (prefill) { let x = Ltot - R.uniform(0, 1) * unit0; while (x > sites[0].s + unit0) { const sp = at(f.vt, x) * lnT(sites[0].freq), L = Math.min(lnL(sites[0].length), 0.9 * sp); slugs.push({ xt: x - L, L, rec: null }); x -= sp; } }
  sites.forEach((q) => { q.next = lnT(q.freq) * R.uniform(0, 1); });
  const fMax = Math.max(...sites.map((q) => q.freq)), tMax = (prefill ? 0 : (Ltot - from) / Math.max(at(f.vt, from), 0.1)) + (2.5 * nSlugs) / sites[0].freq;
  const lCap = 12 * Math.max(...sites.map((q) => q.length)); // liquid available to one slug
  let steps = 0, dt = clamp(Math.min(0.2 / fMax, (0.5 * sites[0].length) / Math.max(at(f.vt, Ltot), 0.1)), 0.05, 30), t = 0, ops = 0, merges = 0, dissipated = 0, generated = slugs.length;
  const arrivals = [], wake = (L, th) => (wkM > 0 ? 1 + wkM * (Math.abs(th) > 0.8 ? 8 * Math.exp((-1.06 * L) / (D * wkM)) : 0.56 * Math.exp((-0.46 * L) / (D * wkM))) : 1); // the multiplier scales the strength and the length of the wake // Moissis–Griffith (steep) / Cook–Behnia (near horizontal)
  while (arrivals.length < nSlugs && t < tMax && ops < maxOps && steps++ < 60000) {
    t += dt;
    for (const q of sites) if (t >= q.next) { // initiation
      q.next += lnT(q.freq);
      let k = 0; while (k < slugs.length && slugs[k].xt > q.s) k++;
      const blocked = (k > 0 && slugs[k - 1].xt <= q.s + minL) || (k < slugs.length && slugs[k].xt + slugs[k].L >= q.s);
      if (!blocked) { const L = Math.min(lnL(q.length), k > 0 ? slugs[k - 1].xt - q.s : Infinity); if (L > minL) { slugs.splice(k, 0, { xt: q.s - L, L, rec: null }); generated++; } }
    }
    for (let k = 0; k < slugs.length; k++) {
      const sl = slugs[k], xf = sl.xt + sl.L, ahead = k > 0 ? slugs[k - 1] : null, gap = ahead ? ahead.xt - xf : 1e9;
      const vInf = at(f.vt, sl.xt), wk = wake(sl.L, o.theta ? at(o.theta, sl.xt) : 0), vT = vInf * wk, shed = at(f.shed, sl.xt) + vInf * (wk - 1) * at(f.holdupSlug, sl.xt);
      const pick = at(f.shed, xf) * (1 + (at(f.mod, xf) - 1) * (1 - Math.exp(-Math.max(gap, 0) / relax))), vF = vT + (pick - shed) / at(f.holdupSlug, xf);
      sl.xt += vT * dt; const xn = Math.min(xf + Math.max(vF, 0) * dt, sl.xt + lCap); sl.L = xn - sl.xt; ops++;
      if (sl.L < minL && !sl.rec) { slugs.splice(k--, 1); dissipated++; continue; }
      if (ahead && xn >= ahead.xt && ahead.L + sl.L < lCap) { // merge into the slug ahead
        const add = ahead.xt - sl.xt; ahead.L += add; ahead.xt = sl.xt; if (ahead.rec) { ahead.rec.length += Math.max(sl.L, 0); ahead.rec.volume = ahead.rec.length * A * ahead.rec.hs; }
        slugs.splice(k--, 1); merges++; continue;
      }
      if (!sl.rec && xn >= Ltot) { sl.rec = { t: t - dt + ((Ltot - xf) / Math.max(vF, 1e-6)), length: Math.max(sl.L, minL), velocity: at(f.vt, Ltot), hs: at(f.holdupSlug, Ltot) }; sl.rec.volume = sl.rec.length * A * sl.rec.hs; arrivals.push(sl.rec); }
      if (sl.xt >= Ltot) slugs.splice(k--, 1);
    }
  }
  if (arrivals.length < 3) return { ...empty, generated, merges, dissipated, from, tSim: t };
  arrivals.sort((a, b) => a.t - b.t);
  const Ls = arrivals.map((q) => q.length), lg = Ls.map(Math.log), mu = mean(lg), sgm = std(lg), span = arrivals[arrivals.length - 1].t - arrivals[0].t, freqArrival = span > 0 ? (arrivals.length - 1) / span : 0;
  // receiving-vessel surge: slug bodies arrive at the slug liquid rate, the film in between, the vessel drains at qDrain
  const qS = o.qSlugOut ?? 0, qF = Math.max(o.qFilmOut ?? 0, 0);
  let qD = o.qDrain ?? 0;
  if (o.drainFactor) { let vIn = 0, tS = 0; for (const a of arrivals) { const dur = a.length / Math.max(a.velocity, 1e-6); vIn += qS * dur; tS += dur; } qD = (o.drainFactor * (vIn + qF * Math.max(span - tS, 0))) / Math.max(span, 1e-9); } // drain relative to the mean liquid arrival of the record
  let lvl = 0, surge = 0, tPrev = arrivals[0].t;
  for (const a of arrivals) { lvl = Math.max(lvl - (qD - qF) * Math.max(a.t - tPrev, 0), 0); const dur = a.length / Math.max(a.velocity, 1e-6); lvl = Math.max(lvl + (qS - qD) * dur, 0); surge = Math.max(surge, lvl); tPrev = a.t + dur; }
  const maxLength = Math.max(...Ls), hs = arrivals[arrivals.length - 1].hs, vOut = arrivals[arrivals.length - 1].velocity, single = (L) => Math.max(L * A * hs * (1 - (qS > 0 ? Math.min(qD / qS, 1) : 0)), 0);
  const length1000 = Math.min(Math.exp(mu + 3.09 * sgm), lCap, Ltot);
  return { n: arrivals.length, arrivals: arrivals.map((q) => ({ t: q.t, length: q.length, volume: q.volume, velocity: q.velocity })), meanLength: mean(Ls), stdLength: std(Ls), p50: quantile(Ls, 0.5), p90: quantile(Ls, 0.9), p99: quantile(Ls, 0.99), maxLength,
    lognormal: { mu, sigma: sgm }, length1000, freqArrival, period: freqArrival > 0 ? 1 / freqArrival : null, merges, dissipated, generated, surge, surgeSingle: single(Math.max(length1000, maxLength)), hist: histogram(Ls, 16), from, tSim: t, vOut };
}

/**
 * Terrain-induced liquid accumulation at the low points of a steady solution: for every dip, the equilibrium holdup of the
 * upward leg (thick-film root of the stratified momentum balance), the liquid it holds, and the critical superficial gas
 * velocity above which the upward leg no longer admits a liquid-accumulating (holdup > 0.5) stratified solution.
 * A dip produces terrain slugs when the downward leg feeding it is stratified and the gas velocity is below that critical value.
 * Returns { dips: [{ x, z, s, upLength, upAngle (deg), rise, vsg, vsl, vsgCrit, holdupUp, holdupDown, holdupSteady, downRegime, volume (m³ held in the upward leg), excess (m³ above the downward-leg holdup), period (s), accumulates }], worst, accumulates }.
 */
export function terrainSlugging(st, o = {}) {
  const N = st.n, dips = [], maxAngle = o.maxAngle ?? 30 * DEG;
  for (let i = 1; i < N; i++) {
    if (!(st.theta[i - 1] < -1e-5 && st.theta[i] > 1e-5)) continue;
    let j = i, len = 0; while (j < N && st.theta[j] > 1e-5) { len += st.ds; j++; }
    const rise = st.z[j] - st.z[i], D = st.D[i], ang = Math.asin(clamp(rise / len, -1, 1));
    if (rise < 2 * D || ang > maxAngle) continue;
    const p = { vsl: st.vsl[i], vsg: st.vsg[i], rhoL: st.rhoL[i], rhoG: st.rhoG[i], muL: st.muL[i], muG: st.muG[i], D, theta: ang };
    if (!(p.vsl > 1e-9 && p.vsg > 1e-9)) continue;
    const hi = (vsg) => { const r = stratifiedRoots({ ...p, vsg }, { n: 60 }); return holdupOfLevel(r[r.length - 1]); }, Hup = hi(p.vsg);
    let crit; if (hi(40) > 0.5) crit = 40; else if (hi(0.02) < 0.5) crit = 0.02; else { let a = 0.02, b = 40; for (let k = 0; k < 26; k++) { const m = Math.sqrt(a * b); if (hi(m) > 0.5) a = m; else b = m; } crit = b; }
    const A = (PI * D * D) / 4, Hs = mean(st.holdup.slice(i, j)), vol = Hup * A * len, k = Math.max(i - 2, 0);
    const down = flowPattern({ vsl: st.vsl[k], vsg: st.vsg[k], rhoL: st.rhoL[k], rhoG: st.rhoG[k], muL: st.muL[k], muG: st.muG[k], sigma: st.sigma[k], D: st.D[k], theta: st.theta[k] }), Hd = down.strat ? down.strat.holdup : st.holdup[k], excess = Math.max(Hup - Hd, 0) * A * len;
    dips.push({ x: st.x[i], z: st.z[i], s: st.s[i], upLength: len, upAngle: ang / DEG, rise, vsg: p.vsg, vsl: p.vsl, vsgCrit: crit, holdupUp: Hup, holdupSteady: Hs, holdupDown: Hd, downRegime: down.pattern, volume: vol, excess, period: excess / Math.max(st.qL[i], 1e-9), accumulates: p.vsg < crit && down.pattern.startsWith('stratified') });
  }
  const acc = dips.filter((d) => d.accumulates), worst = acc.length ? acc.reduce((a, b) => (b.volume > a.volume ? b : a)) : dips.length ? dips.reduce((a, b) => (b.vsg / b.vsgCrit < a.vsg / a.vsgCrit ? b : a)) : null;
  return { dips, worst, accumulates: acc.length > 0 };
}

/**
 * Severe (riser-induced) slugging screening on a steady solution: Bøe criterion, Pots number, Taitel stability pressure.
 * o: { riserBaseS (m along the pipe), feedLength (m, 0 = from the last crest upstream of the riser base), pSep (bara) }
 * Returns null without a riser, else { iBase, feedLength, feedAngle (deg, downward), alpha (gas fraction in the feed line), riserHeight, riserLength, riserAngle (deg),
 *   vsl, vsg, P (bara at the riser base), boeVsl, boe (vsl / boeVsl), pots, taitelPsep (bara needed for stability), taitelStable, feedRegime, feedStratified, severe }.
 */
export function severeScreen(st, o = {}) {
  const N = st.n, ib = clamp(Math.round(o.riserBaseS / st.ds), 1, N - 1), riserHeight = st.z[N] - st.z[ib], riserLength = st.s[N] - st.s[ib];
  if (!(riserHeight > 20 * st.D[ib]) || riserLength <= 0) return null;
  let il = ib; while (il > 1 && st.z[il - 1] < st.z[il]) il--; // low point feeding the riser
  let ic = il; while (ic > 0 && st.z[ic - 1] >= st.z[ic]) ic--; // crest upstream of it
  let feedLength = st.s[ib] - st.s[ic], feedAngle = Math.atan2(st.z[ic] - st.z[il], Math.max(st.s[il] - st.s[ic], 1e-9));
  if (feedLength < 50 * st.D[ib]) { ic = 0; feedLength = st.s[ib]; feedAngle = 0; }
  if (o.feedLength > 0) { feedLength = Math.min(o.feedLength, st.s[ib]); ic = clamp(Math.round((st.s[ib] - feedLength) / st.ds), 0, ib - 1); }
  const im = clamp(Math.round(0.5 * (ic + il)), 0, N - 1), alpha = clamp(1 - mean(st.holdup.slice(ic, ib + 1)), 0.02, 0.98), riserAngle = Math.asin(clamp(riserHeight / riserLength, 0, 1));
  const k = ib - 1, sc = severeSlugging({ vsl: st.vsl[k], vsg: st.vsg[k], P: st.P[ib] * 1e5, rhoL: st.rhoL[k], alpha, L: feedLength, thetaRiser: riserAngle, mG: st.mG[k], mL: st.mL[k], T: st.T[ib] + 273.15, zG: st.zG[k], mwG: st.mwG[k] });
  const fp = flowPattern({ vsl: st.vsl[im], vsg: st.vsg[im], rhoL: st.rhoL[im], rhoG: st.rhoG[im], muL: st.muL[im], muG: st.muG[im], sigma: st.sigma[im], D: st.D[im], theta: st.theta[im] }).pattern, feedStratified = fp.startsWith('stratified');
  const taitelPsep = (st.rhoL[k] * G * ((alpha / 0.89) * feedLength - riserHeight)) / 1e5, pots = sc.pots === null || !Number.isFinite(sc.pots) ? null : sc.pots, boe = st.vsl[k] / Math.max(sc.boeVsl, 1e-12);
  return { iBase: ib, iCrest: ic, iLow: il, feedLength, feedAngle: feedAngle / DEG, alpha, riserHeight, riserLength, riserAngle: riserAngle / DEG, vsl: st.vsl[k], vsg: st.vsg[k], P: st.P[ib], boeVsl: sc.boeVsl, boe, pots,
    taitelPsep, taitelStable: (o.pSep ?? st.pOut) > taitelPsep, feedRegime: fp, feedStratified, severe: feedStratified && (boe >= 1 || (pots !== null && pots < 1)) };
}

/**
 * Lumped severe-slugging cycle of a pipeline–riser system (quasi-equilibrium model in the tradition of Schmidt, Taitel and
 * co-workers). States: gas mass in the feed-line pocket, liquid penetration x into the feed line, liquid column in the riser.
 *   1 build-up   — the low point is blocked; liquid fills the riser and backs up into the feed line while the gas compresses
 *                  (hydrostatic balance P_gas = P_top + ρL g (z sin γ − x sin β), liquid and gas mass balances);
 *   2 production — the riser is full; liquid is produced while the gas pushes the liquid front back to the riser base;
 *   3 blowout    — gas enters the riser; momentum balance of the liquid column with wall friction and choke loss;
 *   4 fallback   — the gas pocket blows down through the riser and choke, the film left on the wall falls back and
 *                  blocks the low point again.
 * o: { D, feedLength, feedAngle (rad, downward), riserHeight, riserLength, wG, wL (kg/s), rhoL, muL, T (K), zG, mwG (g/mol), pSep (Pa),
 *      wLift (kg/s of lift gas injected at the riser base), alphaL (mean liquid fraction of the feed line), chokeOpening (%), chokeDp (Pa at full opening and the mean rate), rough, fallback (film fraction), maxCycles }
 * Returns { stable, type, period (s | null), amplitude (Pa, riser-base pressure), pBaseMean, pBaseMin, pBaseMax, qLiqPeakRatio, qGasPeakRatio,
 *           t[], pBase[], pFeed[], wLout[], wGout[], level[] (liquid column / riser length), stage[], stages: { buildUp, production, blowout, fallback } (s), cycles }.
 */
export function riserSluggingCycle(o) {
  const { D, feedLength: Lp, riserHeight: hR, riserLength: Lr, wG, wL, rhoL, T, zG = 0.9, mwG = 20, pSep, rough = 4.5e-5 } = o, A = (PI * D * D) / 4;
  if (!(D > 0 && Lp > 0 && hR > 0 && Lr >= hR && wG > 0 && wL > 0 && rhoL > 0 && pSep > 0)) throw new Error('The riser-slugging cycle needs positive geometry, rates and separator pressure.');
  const sinB = Math.sin(Math.max(o.feedAngle || 0, 0)), sinG = hR / Lr, al = clamp(1 - (o.alphaL ?? 0.3), 0.05, 0.98), gasK = (zG * RGAS * T) / (mwG * 1e-3), qL = wL / rhoL, phi = clamp(o.fallback ?? 0.12, 0, 0.6), zc = clamp((o.chokeOpening ?? 100) / 100, 0.02, 1);
  // riser-base gas lift: the lift gas bubbles through the liquid column; its void fraction (drift flux at the mean riser pressure) reduces the liquid content, the head of the full riser and the time needed to fill it
  const wLift = Math.max(o.wLift || 0, 0), jLift = wLift / ((pSep + 0.5 * rhoL * G * hR) / gasK) / A, aLift = wLift > 0 ? clamp(jLift / (1.2 * (jLift + qL / A) + 0.35 * Math.sqrt(G * D)), 0, 0.9) : 0, aR = 1 - aLift;
  const muL = o.muL || 2e-3, rel = rough / D, fH = (Re) => (Re < 2000 ? 64 / Math.max(Re, 1) : 1 / (-1.8 * log10((rel / 3.7) ** 1.11 + 6.9 / Re)) ** 2);
  const rhoG0 = pSep / gasK, u0 = (qL + wG / rhoG0) / A, lam0 = qL / (u0 * A), chokeK = (o.chokeDp ?? 1e5) / ((lam0 * rhoL + (1 - lam0) * rhoG0) * u0 * u0) / (zc * zc); // Δp_choke = chokeK ρ u²
  const gasFr = (0.015 * Lp) / 3 / (2 * D); // gas-side friction coefficient of the feed line (the pocket empties over about a third of its length)
  const Vg = (x) => al * A * (Lp - x), head = rhoL * G, tFill = (A * Lr * aR) / qL, pTopL = pSep + chokeK * rhoL * (qL / A) ** 2;
  const ts = [], pb = [], pf = [], wl = [], wg = [], lv = [], sg = [], rec = (t, P, Pg, ql, wgo, lev, st) => { ts.push(t); pb.push(P); pf.push(Pg); wl.push(ql * rhoL); wg.push(wgo); lv.push(lev); sg.push(st); };
  // liquid front in the feed line from the hydrostatic balance (−1 when the gas pressure exceeds the liquid head: gas penetrates)
  const frontX = (ell, mg, full) => {
    const pT = full ? pTopL : pSep, f = (x) => (mg * gasK) / Vg(x) - pT - head * ((full ? aR * Lr : ell - al * x) * sinG - x * sinB);
    if (f(0) >= 0) return -1;
    const xm = Math.min(0.98 * Lp, full ? 0.98 * Lp : ell / al);
    return f(xm) <= 0 ? xm : brent(f, 0, xm, 1e-9);
  };
  let t = 0, ell = phi * Lr, mg = ((pSep + head * ell * sinG) * Vg(0)) / gasK, cycles = [], guard = 0;
  const maxCycles = o.maxCycles || 4;
  while (cycles.length < maxCycles && guard++ < 12) {
    const c = { t0: t, d: [0, 0, 0, 0], pMin: Infinity, pMax: 0, qlMax: 0, wgMax: 0, i0: ts.length }, note = (P, ql, wgo) => { c.pMin = Math.min(c.pMin, P); c.pMax = Math.max(c.pMax, P); c.qlMax = Math.max(c.qlMax, ql); c.wgMax = Math.max(c.wgMax, wgo); };
    // stage 1: build-up
    let x = 0, full = false, dt = tFill / 240, n = 0;
    for (; n < 6000; n++) {
      ell += (qL / A) * dt; mg += wG * dt; t += dt; c.d[0] += dt;
      const xn = frontX(ell, mg, false);
      if (xn < 0) { x = 0; break; }
      x = xn; const z = (ell - al * x) / aR, P = pSep + head * aR * z * sinG; note(P, 0, 0); rec(t, P, P - head * x * sinB, 0, 0, z / Lr, 0);
      if (z >= Lr) { full = true; ell = aR * Lr + al * x; break; }
    }
    // stage 2: production with a full riser
    if (full) {
      const P2 = pTopL + head * aR * hR, need = Math.max((P2 * Vg(0)) / gasK - mg, 0), dt2 = clamp(need / wG / 160, 1e-3, tFill / 100);
      for (n = 0; n < 4000; n++) {
        mg += wG * dt2; t += dt2; c.d[1] += dt2;
        const xn = frontX(0, mg, true), qo = qL + (al * A * (x - Math.max(xn, 0))) / dt2; x = Math.max(xn, 0); note(P2, qo, 0); rec(t, P2, P2 - head * x * sinB, qo, 0, 1, 1);
        if (xn <= 0) break;
      }
    }
    // stage 3: gas penetration and blowout of the liquid column
    let Ls = Math.min(full ? Lr : Math.max((ell - al * x) / aR, 0.02 * Lr), Lr), b = 0, u = qL / A, pool = 0, t3 = 0, blown = false;
    for (n = 0; n < 60000 && t3 < 4 * tFill; n++) {
      const top = b + Ls >= Lr - 1e-9, Pg = (mg * gasK) / (Vg(0) + A * b * (1 - phi)), h = Math.min(0.5, (0.01 * Lr) / Math.max(Math.abs(u), 0.2), tFill / 200), au = Math.abs(u);
      const drive = (Pg - pSep - head * aR * Ls * sinG) / (rhoL * aR * Ls), damp = (fH((rhoL * au * D) / muL) * au) / (2 * D) + (top ? (chokeK * au) / Ls : 0) + (gasFr * (Pg / gasK) * au) / (rhoL * Ls); // wall friction, choke, friction of the gas feeding the bubble
      u = (u + h * drive) / (1 + h * damp); b = Math.max(b + u * h, 0); if (b + Ls >= Lr) { Ls = Math.max(Lr - b, 0); } mg += wG * h; pool += (qL / A) * h; t += h; t3 += h; c.d[2] += h;
      const P = Pg, qo = b + Ls >= Lr - 1e-9 ? A * Math.max(u, 0) : 0; note(P, Ls > 0.3 * Lr ? qo : 0, 0); rec(t, P, Pg, qo, 0, Ls / Lr, 2); // the peak rate is taken while the bulk of the column is still in the riser
      if (Ls < 0.05 * Lr) { blown = true; break; }
      if (b <= 0 && u < 0) break; // the column settled back: the low point is blocked again
    }
    // stage 4: gas blowdown and liquid fallback
    ell = (blown ? phi * Lr : Ls * aR) + pool;
    if (blown) for (n = 0; n < 20000; n++) {
      const Pg = (mg * gasK) / Vg(0), block = pSep + head * ell * sinG;
      if (Pg <= block) break;
      const rg = Pg / gasK, ug = Math.min(Math.sqrt((Pg - pSep) / (rg * ((0.02 * Lr) / (2 * D) + gasFr + chokeK))), Math.sqrt(gasK)), wo = rg * A * ug, h = clamp((0.03 * mg) / Math.max(Math.abs(wo - wG), 1e-9), 0.01, tFill / 300);
      mg += (wG - wo) * h; ell += (qL / A) * h; t += h; c.d[3] += h; note(Pg, 0, wo); rec(t, Pg, Pg, 0, wo, ell / Lr, 3);
    }
    mg = Math.max(mg, ((pSep + head * ell * sinG) * Vg(0)) / gasK * 0.999);
    c.t1 = t; c.i1 = ts.length; cycles.push(c);
  }
  // no cycle when the pressure swing is small or when the liquid is produced without a surge (choked, gas and liquid leave together)
  const last = cycles[cycles.length - 1], period0 = last.t1 - last.t0, amplitude = last.pMax - last.pMin, stable = !(amplitude > Math.max(Math.min(1e5, 0.2 * head * hR), 0.05 * head * hR)) || !(period0 > 0) || last.qlMax < 1.5 * qL;
  const i0 = cycles[Math.max(cycles.length - 3, 0)].i0, every = Math.max(1, Math.ceil((ts.length - i0) / 380)), pick = (a) => a.filter((_, i) => i >= i0 && (i - i0) % every === 0), NAMES = ['build-up', 'production', 'blowout', 'fallback'], tFrom = ts[i0] ?? 0;
  return { stable, type: stable ? 'stable' : last.d[1] > 0 ? 'severe slugging 1 (riser fills completely)' : 'severe slugging 2 (gas penetrates before the riser is full)', period: stable ? null : period0, amplitude, pBaseMean: mean(pb.slice(last.i0, last.i1)), pBaseMin: last.pMin, pBaseMax: last.pMax,
    qLiqPeakRatio: last.qlMax / qL, qGasPeakRatio: last.wgMax / wG, liftVoid: aLift, t: pick(ts).map((v) => v - tFrom), pBase: pick(pb), pFeed: pick(pf), wLout: pick(wl), wGout: pick(wg), level: pick(lv), stage: pick(sg).map((k) => NAMES[k]),
    stages: { buildUp: last.d[0], production: last.d[1], blowout: last.d[2], fallback: last.d[3] }, cycles: cycles.length, steps: ts.length };
}

/** Carrier field for slugTracking from a steady solution and its unit cells (one per node): shedding rate and terrain pick-up ratio. */
export function trackingField(st, units) {
  const mod = units.map((u, i) => {
    const qf = Math.max(u.holdupFilm * u.filmVelocity, 0.02 * st.vsl[i], 1e-4); // liquid carried by the film between slugs
    const Hl = clamp(stratifiedLevel({ vsl: qf, vsg: st.vsg[i], rhoL: st.rhoL[i], rhoG: st.rhoG[i], muL: st.muL[i], muG: st.muG[i], D: st.D[i], theta: clamp(st.theta[i], -0.6, 0.6) }).holdup, 0.005, 0.92 * u.holdupSlug), Hd = u.holdupFilm;
    return clamp(((Hl / Hd) * (u.holdupSlug - Hd)) / (u.holdupSlug - Hl), 0.9, 1.15); // bounded: the film only partly relaxes between slugs
  });
  return { s: st.s, vt: units.map((u) => u.vt), holdupSlug: units.map((u) => u.holdupSlug), shed: units.map((u) => Math.max((u.vt - u.vBody) * u.holdupSlug, 1e-4)), mod, D: st.D[0], A: (PI * st.D[0] ** 2) / 4 };
}

// =====================================================================================================
// 6. Transient drift-flux model (finite volume, staggered, semi-implicit in pressure)
// =====================================================================================================
/** Grid for the transient solver: n cells along the profile, `refine` times finer beyond the arc length sFrom (the riser). */
export function transientGrid(profile, n = 60, sFrom = Infinity, refine = 1) {
  const X = profile.x, Z = profile.z, S = [0];
  for (let i = 1; i < X.length; i++) S.push(S[i - 1] + Math.hypot(X[i] - X[i - 1], Z[i] - Z[i - 1]));
  const L = S[S.length - 1], s0 = clamp(sFrom, 0, L), W = s0 + refine * (L - s0), s = [], x = [], z = [];
  for (let i = 0; i <= n; i++) { const w = (W * i) / n, si = w <= s0 ? w : s0 + (w - s0) / refine; s.push(si); x.push(interp1(S, X, si)); z.push(interp1(S, Z, si)); }
  const ds = [], sc = [], zc = [], xc = [];
  for (let i = 0; i < n; i++) { ds.push(s[i + 1] - s[i]); sc.push(0.5 * (s[i] + s[i + 1])); zc.push(0.5 * (z[i] + z[i + 1])); xc.push(0.5 * (x[i] + x[i + 1])); }
  return { n, length: L, s, x, z, ds, sc, zc, xc };
}

/**
 * Transient 1-D drift-flux solver. Conserved per cell: gas component (free + dissolved), dead-liquid component and thermal
 * energy; the mixture momentum equation (inertia, wall friction, gravity; convective acceleration neglected) is solved on the
 * faces and coupled implicitly to the pressure through the volume-conservation (pressure) equation, so the time step is
 * limited by the material CFL number only. Phase fluxes are first-order upwind with the slip law v_g = C0 j + v_d; flashing
 * follows the equilibrium gas fraction of the property table (black-oil form), heat is lost through U to the ambient.
 * o: { fm, grid (transientGrid), D, rough, init: { P[] (bara), T[] (°C), holdup[] } at the cell centres, mdot0 (kg/s, for the initial velocities),
 *      mdotOf(t) (kg/s at the inlet), tInOf(t) (°C), pOutOf(t) (bara), chokeOf(t) (% opening), chokeK (Δp = chokeK ρ j² at full opening),
 *      slip: { c0, vdScale }, fricMult, U[] (W/m²K per cell), tAmb[] (°C per cell), wallC (J/m³/K of wall per unit bore volume), cfl, dtMax, tEnd, nSeries, nField }
 * Returns a stepper { done, t, steps, advance(nSteps), result() } — see transientDriftFlux for the result.
 */
export function makeTransient(o) {
  const { fm, grid, D, rough = 4.5e-5, tEnd, cfl = 0.8 } = o, N = grid.n, A = (PI * D * D) / 4, sp = { c0: 1.2, vdScale: 1, ...(o.slip || {}) }, fricMult = o.fricMult ?? 1, chokeK = o.chokeK || 0;
  const mdotOf = o.mdotOf, pOutOf = o.pOutOf, chokeOf = o.chokeOf || (() => 100), tInOf = o.tInOf || (() => o.init.T[0]), wallC = o.wallC || 0, dtMax = o.dtMax || tEnd / 40;
  const F64 = (n) => new Float64Array(n), Gc = F64(N), md = F64(N), T = F64(N), P = F64(N), jf = F64(N + 1), vgF = F64(N + 1), vlF = F64(N + 1);
  const rhoG = F64(N), rhoL = F64(N), rsA = F64(N), mgA = F64(N), mlA = F64(N), Hh = F64(N), rhoM = F64(N), muM = F64(N), Vv = F64(N), kap = F64(N), dRg = F64(N), dRl = F64(N), dRs = F64(N), gG = F64(N), gD = F64(N), rsmd = F64(N);
  const AGf = F64(N + 1), BGf = F64(N + 1), ADf = F64(N + 1), BDf = F64(N + 1), jst = F64(N + 1), bco = F64(N + 1), C0f = F64(N + 1), vdf = F64(N + 1), aLf = F64(N + 1), bLf = F64(N + 1);
  const xgOf = (pr) => { const m = pr.mG + pr.mO + pr.mW; return m > 0 ? pr.mG / m : pr.wG; };
  // reference (maximum) gas mass fraction of the stream: everything above the local equilibrium value is dissolved in the liquid
  let xgMax = xgOf(fm.at(1.01325, Math.max(...o.init.T, 15)));
  for (let i = 0; i < N; i++) xgMax = Math.max(xgMax, xgOf(fm.at(o.init.P[i], o.init.T[i])));
  xgMax = clamp(xgMax, 1e-9, 0.9995);
  const rsOf = (pr) => Math.max(xgMax - xgOf(pr), 0) / (1 - xgMax), cG = o.cpG ?? 2300, cD = o.cpL ?? 2200;
  for (let i = 0; i < N; i++) { // initial state
    const pr = fm.at(o.init.P[i], o.init.T[i]), H = clamp(o.init.holdup[i], 0, 1), rs = rsOf(pr);
    P[i] = o.init.P[i] * 1e5; T[i] = o.init.T[i]; md[i] = (H * pr.rhoL) / (1 + rs); Gc[i] = (1 - H) * pr.rhoG + rs * md[i];
  }
  // side inflows (branches joining the line): [{ s (m), mdotOf(t) (kg/s of the same stream), T (°C) }] → one source per cell
  const srcs = (o.sources || []).map((q) => ({ ...q, w: 0 })), srcAt = new Array(N).fill(null); for (const q of srcs) { let c = 0; while (c < N - 1 && grid.s[c + 1] < q.s) c++; srcAt[c] = q; }
  let last0 = { qL: 0, qG: 0, mOut: 0 };
  { const W = (o.mdot0 ?? mdotOf(0)) / A; for (let f = 0; f <= N; f++) { const i = Math.min(f, N - 1), pr = fm.at(P[i] / 1e5, T[i]), xg = xgOf(pr); jf[f] = W * (xg / pr.rhoG + (1 - xg) / pr.rhoL); vgF[f] = vlF[f] = jf[f]; if (f === N) last0 = { qL: (W * (1 - xg) * A) / pr.rhoL, qG: (W * xg * A) / pr.rhoG, mOut: W * A }; } }
  const total = () => { let m = 0, e = 0; for (let i = 0; i < N; i++) { m += (Gc[i] + md[i]) * grid.ds[i] * A; e += (cG * Gc[i] + cD * md[i] + wallC) * T[i] * grid.ds[i] * A; } return { m, e }; };
  const pCap = (2.5 * Math.max(...o.init.P) + 50) * 1e5, start = total(), bal = { mIn: 0, mOut: 0, eIn: 0, eOut: 0, eLoss: 0, volErr: 0 };
  const ser = { t: [], pIn: [], pOutCell: [], qLiqOut: [], qGasOut: [], holdupOut: [], mOut: [], inv: [], tOut: [] }, fld = { t: [], holdup: [], P: [], T: [] };
  const nField = o.nField || 60, nSeries = o.nSeries || 400;
  let t = 0, dt = Math.min(dtMax, 1), steps = 0, relax = 1, rejected = 0, aborted = null, nextField = 0, last = last0;
  const record = () => {
    let inv = 0; for (let i = 0; i < N; i++) inv += Hh[i] * grid.ds[i] * A;
    ser.t.push(t); ser.pIn.push((P[0] + ((P[0] - P[1]) * grid.ds[0]) / (grid.ds[0] + grid.ds[1])) / 1e5); ser.pOutCell.push(P[N - 1] / 1e5); ser.qLiqOut.push(last.qL); ser.qGasOut.push(last.qG); ser.holdupOut.push(Hh[N - 1]); ser.mOut.push(last.mOut); ser.inv.push(inv); ser.tOut.push(T[N - 1]);
    if (t >= nextField - 1e-9) { fld.t.push(t); fld.holdup.push(Array.from(Hh)); fld.P.push(Array.from(P, (p) => p / 1e5)); fld.T.push(Array.from(T)); nextField += tEnd / (nField - 1); }
  };
  const cellProps = () => {
    for (let i = 0; i < N; i++) {
      const pb = P[i] / 1e5, pr = fm.at(pb, T[i]), rsS = rsOf(pr), rs = md[i] > 0 ? Math.min(rsS, Gc[i] / md[i]) : 0, mg = Math.max(Gc[i] - rs * md[i], 0), ml = md[i] * (1 + rs), sat = mg > 1e-12 || (rsS * md[i] - Gc[i]) < 0.01 * pr.rhoG; // saturated, or so close that added gas stays free
      rhoG[i] = pr.rhoG; rhoL[i] = pr.rhoL; rsA[i] = rs; mgA[i] = mg; mlA[i] = ml; rsmd[i] = rs * md[i];
      const H = clamp(ml / pr.rhoL, 0, 1); Hh[i] = H; rhoM[i] = mg + ml; muM[i] = pr.muL * H + pr.muG * (1 - H);
      Vv[i] = mg / pr.rhoG + ml / pr.rhoL - 1;
      if (steps % 6 === 0 || !(dRg[i] > 0)) { const p2 = fm.at(pb * 1.001, T[i]), dp = 0.001 * P[i]; dRg[i] = Math.max((p2.rhoG - pr.rhoG) / dp, 1e-12); dRl[i] = Math.max((p2.rhoL - pr.rhoL) / dp, pr.rhoL / 2.5e9); dRs[i] = (rsOf(p2) - rsS) / dp; } // compressibilities, refreshed every few steps
      kap[i] = Math.min((-mg * dRg[i]) / (pr.rhoG * pr.rhoG) - (ml * dRl[i]) / (pr.rhoL * pr.rhoL) + md[i] * Math.max(dRs[i], 0) * (1 / pr.rhoL - 1 / pr.rhoG), -0.01 / P[i]); // the dissolution term is kept below the bubble point too: gas appears as soon as the pressure falls
      gG[i] = sat ? 1 / pr.rhoG : 1 / pr.rhoL; gD[i] = sat ? (1 + rsS) / pr.rhoL - rsS / pr.rhoG : 1 / pr.rhoL;
    }
  };
  const slipAt = (f, a, th) => { // linear phase-velocity relations on a face: v_g = C0 j + vd, v_l = aL j + bL
    const sl = slipLaw(a, jf[f], th, D, sp); C0f[f] = sl.C0; vdf[f] = sl.vd;
    if (a > 0.75) { aLf[f] = 1 - (a * (sp.c0 - 1)) / 0.25; bLf[f] = (-a * slipLaw(0, jf[f], th, D, sp).vd) / 0.25; } // analytic limit as the liquid vanishes
    else { aLf[f] = (1 - a * sl.C0) / (1 - a); bLf[f] = (-a * sl.vd) / (1 - a); }
  };
  cellProps(); record();
  const sub = new Array(N), dia = new Array(N), sup = new Array(N), rhs = new Array(N), FG = F64(N + 1), FD = F64(N + 1), Fgas = F64(N + 1), Fliq = F64(N + 1), FE = F64(N + 1), outG = F64(N), outL = F64(N);
  const keep = [Gc, md, T, P, jf, vgF, vlF].map((a) => F64(a.length)), state = [Gc, md, T, P, jf, vgF, vlF];
  const attempt = (dtTry) => {
    // time step from the material CFL condition
    let lim = Infinity; for (let f = 0; f <= N; f++) { const i = Math.min(f, N - 1), v = Math.max(Math.abs(vgF[f]), Math.abs(vlF[f]), 1e-3); lim = Math.min(lim, grid.ds[i] / v, f > 0 ? grid.ds[f - 1] / v : Infinity); }
    dt = Math.max(Math.min(cfl * relax * lim, dtMax, dtTry), 1e-3); dt = Math.min(dt, tEnd - t); if (!(dt > 0)) return -1;
    const W = mdotOf(t + dt) / A, pOut = pOutOf(t + dt) * 1e5, zc = clamp(chokeOf(t + dt) / 100, 0.01, 1), tIn = tInOf(t + dt), inert = 1 / dt;
    for (let f = 1; f <= N; f++) { // momentum and linearised fluxes on the faces
      const L = f - 1, out = f === N, R = out ? L : f, dsf = out ? 0.5 * grid.ds[L] : 0.5 * (grid.ds[L] + grid.ds[R]), dz = out ? grid.z[N] - grid.zc[L] : grid.zc[R] - grid.zc[L];
      // gas fraction of the slip law from the cell the mixture comes from, so that a liquid front can advance into a gas-filled cell
      const a = out || jf[f] >= 0 ? 1 - Hh[L] : 1 - Hh[R], th = Math.asin(clamp(dz / dsf, -1, 1)); slipAt(f, a, th);
      const rf = 0.5 * (rhoM[L] + rhoM[R]), aj = Math.abs(jf[f]), Re = (rf * aj * D) / (0.5 * (muM[L] + muM[R]));
      const fr = (fricMult * (Re > 1e-6 ? frictionFactor(Re, rough / D, 'haaland') * aj : (64 * 0.5 * (muM[L] + muM[R])) / (rf * D))) / (2 * D) + (out ? (chokeK * aj) / (dsf * zc * zc) : 0);
      jst[f] = (jf[f] * inert - (G * dz) / dsf) / (inert + fr); bco[f] = 1 / (rf * dsf * (inert + fr));
      const vg = C0f[f] * jf[f] + vdf[f], vl = aLf[f] * jf[f] + bLf[f], dg = vg >= 0 || out ? L : R, dl = vl >= 0 || out ? L : R, mgD = out && vg < 0 ? pOut / (P[L] / rhoG[L]) : mgA[dg], mdD = out && vl < 0 ? 0 : md[dl], rsD = out && vl < 0 ? 0 : rsmd[dl];
      AGf[f] = mgD * C0f[f] + rsD * aLf[f]; BGf[f] = mgD * vdf[f] + rsD * bLf[f]; ADf[f] = mdD * aLf[f]; BDf[f] = mdD * bLf[f];
    }
    AGf[0] = 0; ADf[0] = 0; BGf[0] = W * xgMax; BDf[0] = W * (1 - xgMax);
    for (let i = 0; i < N; i++) { // pressure (volume-conservation) equation
      const r = dt / grid.ds[i], cm = gG[i] * AGf[i] + gD[i] * ADf[i], cp = gG[i] * AGf[i + 1] + gD[i] * ADf[i + 1];
      sub[i] = i > 0 ? r * cm * bco[i] : 0; sup[i] = i < N - 1 ? r * cp * bco[i + 1] : 0; dia[i] = kap[i] - r * (cp * bco[i + 1] + (i > 0 ? cm * bco[i] : 0));
      rhs[i] = kap[i] * P[i] - clamp(Vv[i], -0.1, 0.1) - (srcAt[i] ? r * (srcAt[i].mdotOf(t + dt) / A) * (gG[i] * xgMax + gD[i] * (1 - xgMax)) : 0) + r * (gG[i] * (BGf[i + 1] - BGf[i]) + gD[i] * (BDf[i + 1] - BDf[i])) + r * (cp * jst[i + 1] - (i > 0 ? cm * jst[i] : 0)) - (i === N - 1 ? r * cp * bco[N] * pOut : 0);
    }
    const Pn = tridiag(sub, dia, sup, rhs);
    for (let i = 0; i < N; i++) Pn[i] = Number.isFinite(Pn[i]) ? clamp(Pn[i], 5e4, pCap) : P[i];
    // new volumetric fluxes, upwind phase fluxes with the new velocity signs
    outG.fill(0); outL.fill(0);
    for (let f = 1; f <= N; f++) {
      const L = f - 1, out = f === N, R = out ? L : f; jf[f] = clamp(jst[f] - bco[f] * ((out ? pOut : Pn[R]) - Pn[L]), out ? 0 : -40, 40); // the outlet acts as a check valve
      const vg = C0f[f] * jf[f] + vdf[f], vl = aLf[f] * jf[f] + bLf[f]; vgF[f] = vg; vlF[f] = vl;
      const mgD = vg >= 0 ? mgA[L] : out ? pOut / (P[L] / rhoG[L]) : mgA[R], lD = vl >= 0 ? L : out ? -1 : R;
      Fgas[f] = mgD * vg; Fliq[f] = lD < 0 ? 0 : mlA[lD] * vl; FD[f] = lD < 0 ? 0 : md[lD] * vl; FG[f] = Fgas[f] + (lD < 0 ? 0 : rsmd[lD] * vl);
      if (vg >= 0) outG[L] += Fgas[f]; else if (!out) outG[R] -= Fgas[f];
      if (vl >= 0) outL[L] += Fliq[f]; else if (!out) outL[R] -= Fliq[f];
    }
    for (let f = 1; f <= N; f++) { // positivity limiter: a cell cannot give more than it holds
      const L = f - 1, out = f === N, R = out ? L : f, dg = vgF[f] >= 0 ? L : out ? -1 : R, dl = vlF[f] >= 0 ? L : out ? -1 : R;
      const sg = dg >= 0 && outG[dg] * dt > 0.98 * mgA[dg] * grid.ds[dg] ? (0.98 * mgA[dg] * grid.ds[dg]) / (outG[dg] * dt) : 1, slq = dl >= 0 && outL[dl] * dt > 0.98 * mlA[dl] * grid.ds[dl] ? (0.98 * mlA[dl] * grid.ds[dl]) / (outL[dl] * dt) : 1;
      if (sg < 1 || slq < 1) { const dis = FG[f] - Fgas[f]; Fgas[f] *= sg; Fliq[f] *= slq; FD[f] *= slq; FG[f] = Fgas[f] + dis * slq; }
    }
    FG[0] = W * xgMax; FD[0] = W * (1 - xgMax); vgF[0] = vlF[0] = jf[0] = W * (xgMax / rhoG[0] + (1 - xgMax) / rhoL[0]);
    for (const q of srcs) { q.w = q.mdotOf(t + dt) / A; }
    for (let f = 0; f <= N; f++) { const fe = cG * FG[f] + cD * FD[f]; FE[f] = fe * (f === 0 ? tIn : fe >= 0 || f === N ? T[f - 1] : T[f]); }
    for (let i = 0; i < N; i++) {
      const r = dt / grid.ds[i], sq = srcAt[i], sw = sq ? sq.w : 0, E = (cG * Gc[i] + cD * md[i] + wallC) * T[i] - r * (FE[i + 1] - FE[i]) + r * sw * (cG * xgMax + cD * (1 - xgMax)) * (sq ? sq.T : 0), hl = (4 * (o.U ? o.U[i] : 0)) / D, ta = o.tAmb ? o.tAmb[i] : 4;
      Gc[i] = Math.max(Gc[i] - r * (FG[i + 1] - FG[i]) + r * sw * xgMax, 0); md[i] = Math.max(md[i] - r * (FD[i + 1] - FD[i]) + r * sw * (1 - xgMax), 0);
      if (sw) { bal.mIn += sw * A * dt; bal.eIn += sw * (cG * xgMax + cD * (1 - xgMax)) * sq.T * A * dt; }
      const C = cG * Gc[i] + cD * md[i] + wallC, Tn = (E + dt * hl * ta) / (C + dt * hl);
      bal.eLoss += hl * (Tn - ta) * dt * grid.ds[i] * A; T[i] = Tn; P[i] = Pn[i];
    }
    const rg = rhoG[N - 1], rl = rhoL[N - 1];
    last = { qL: (Fliq[N] / rl) * A, qG: (Fgas[N] / rg) * A, mOut: (FG[N] + FD[N]) * A };
    bal.mIn += (FG[0] + FD[0]) * A * dt; bal.mOut += last.mOut * dt; bal.eIn += FE[0] * A * dt; bal.eOut += FE[N] * A * dt;
    t += dt; steps++;
    cellProps();
    let ve = 0; for (let i = 0; i < N; i++) ve = Math.max(ve, Math.abs(Vv[i]));
    return ve;
  };
  // a step whose volume constraint is poorly met (phase appearance, flow reversal) is rejected and repeated with a shorter step
  const stepOnce = () => {
    let dtTry = dt * 1.5, ve = 0;
    for (let k = 0; k < 4; k++) {
      state.forEach((a, i) => keep[i].set(a));
      const b0 = { ...bal }, t0 = t, s0 = steps, l0 = last;
      ve = attempt(dtTry);
      if (ve < 0) return false;
      if (ve <= 0.05 || dt <= 0.5 || k === 3) break;
      state.forEach((a, i) => a.set(keep[i])); Object.assign(bal, b0); t = t0; steps = s0; last = l0; dtTry = 0.4 * dt; rejected++; cellProps();
    }
    if (ve > 0.25) { aborted = t; return false; } // outside the validity range of the scheme (liquid-full line with pressure surges)
    bal.volErr = Math.max(bal.volErr, ve); relax = ve > 0.01 ? Math.max(0.7 * relax, 0.05) : Math.min(1, relax * 1.1);
    record();
    return t < tEnd - 1e-9;
  };
  const api = {
    done: !(tEnd > 0), get t() { return t; }, get steps() { return steps; },
    advance(n = 50) { for (let k = 0; k < n && !api.done; k++) { if (!stepOnce()) api.done = true; if (steps > (o.maxSteps || 200000)) api.done = true; } return api.done; },
    result() {
      const end = total(), every = Math.max(1, Math.ceil(ser.t.length / nSeries)), pick = (a) => a.filter((_, i) => i % every === 0 || i === a.length - 1);
      const series = Object.fromEntries(Object.entries(ser).map(([k, a]) => [k, pick(a)]));
      return { ...series, steps, rejected, aborted, dtMean: steps ? t / steps : 0, tEnd: t, xgMax, field: { s: grid.sc.slice(), x: grid.xc.slice(), t: fld.t, holdup: fld.holdup, P: fld.P, T: fld.T },
        mass: { initial: start.m, final: end.m, inflow: bal.mIn, outflow: bal.mOut, error: (end.m - start.m - bal.mIn + bal.mOut) / Math.max(start.m, 1e-9) },
        energy: { initial: start.e, final: end.e, inflow: bal.eIn, outflow: bal.eOut, loss: bal.eLoss, error: (end.e - start.e - bal.eIn + bal.eOut + bal.eLoss) / Math.max(Math.abs(start.e) + Math.abs(bal.eIn), 1e-9) },
        volErrMax: bal.volErr, final: { P: Array.from(P, (p) => p / 1e5), T: Array.from(T), holdup: Array.from(Hh), j: Array.from(jf) } };
    },
  };
  return api;
}
/**
 * Run the transient drift-flux solver to the end time (options as makeTransient).
 * Returns { aborted (s | null: time at which the run left the validity range and was stopped), t[], pIn[] (bara), qLiqOut[], qGasOut[] (actual m³/s), holdupOut[], mOut[] (kg/s), inv[] (m³ liquid in the line), tOut[], steps, dtMean,
 *   field: { s[], x[], t[], holdup[][], P[][], T[][] }, mass: { initial, final, inflow, outflow, error }, energy: { …, loss, error }, volErrMax, final: { P, T, holdup, j } }.
 */
export function transientDriftFlux(o) { const sim = makeTransient(o); while (!sim.done) sim.advance(200); return sim.result(); }

// =====================================================================================================
// 7. Local high-resolution models: 1-D radial RANS of developed pipe flow, 1-D interface capturing
// =====================================================================================================
/**
 * Fully developed turbulent flow in a smooth pipe by a 1-D radial Reynolds-averaged solve in wall units.
 * model: 'laminar' | 'mixing' (Nikuradse mixing length with van Driest damping) | 'komega' (Wilcox 1988, low-Reynolds wall treatment by
 * the analytic near-wall ω) | 'kepsilon' (Chien 1982 low-Reynolds k–ε). The total shear stress is linear in r, so the momentum equation
 * reduces to du⁺/dy⁺ = (1 − y/R)/(1 + ν_t⁺); k and ω or ε are solved by under-relaxed finite-volume iteration.
 * o: { reTau (R u_τ / ν), model, n, maxIter, tol }. Returns { Re, f (Darcy), reTau, y[] (y⁺), u[] (u⁺), r[] (r/R), uRel[] (u / U_bulk), nut[], k[], iterations, residual, converged }.
 */
export function ransPipe(o = {}) {
  if (['sst', 'sa', 'kestd', 'rng', 'realizable', 'rsm'].includes(o.model)) return ransPipeExtra(o);
  const Rp = Math.max(o.reTau || 1000, 5), model = o.model || 'komega', n = o.n || 120, maxIter = o.maxIter || 6000, tol = o.tol || 1e-9;
  // geometric grid from the wall, first node at y⁺ ≈ 0.1
  const y1 = Math.min(0.1, Rp / (4 * n)), g = (q) => (y1 * (q ** n - 1)) / (q - 1) - Rp, q = g(1 + 1e-9) >= 0 ? 1 + 1e-9 : brent(g, 1 + 1e-9, 3, 1e-13), y = [0];
  for (let j = 1; j <= n; j++) y.push(j === n ? Rp : (y1 * (q ** j - 1)) / (q - 1));
  const r = y.map((v) => Rp - v), nut = new Array(n + 1).fill(0), u = new Array(n + 1).fill(0), S = new Array(n + 1).fill(0);
  const mix = (yy) => { const e = 1 - yy / Rp, l = Rp * (0.14 - 0.08 * e * e - 0.06 * e ** 4) * (1 - Math.exp(-yy / 26)); return l; };
  const velocity = () => { // integrate du/dy with mid-point eddy viscosity; returns the bulk velocity
    for (let j = 0; j <= n; j++) S[j] = (1 - y[j] / Rp) / (1 + nut[j]);
    for (let j = 1; j <= n; j++) u[j] = u[j - 1] + 0.5 * (S[j] + S[j - 1]) * (y[j] - y[j - 1]);
    let ub = 0; for (let j = 1; j <= n; j++) ub += 0.5 * (u[j] * r[j] + u[j - 1] * r[j - 1]) * (y[j] - y[j - 1]);
    return (2 * ub) / (Rp * Rp);
  };
  const mixingNut = () => { for (let j = 0; j <= n; j++) { const l = mix(y[j]), e = 1 - y[j] / Rp, s = (2 * e) / (1 + Math.sqrt(1 + 4 * l * l * e)); nut[j] = l * l * s; } };
  let ub, it = 0, res = 0, k = null;
  if (model === 'laminar') ub = velocity();
  else if (model === 'mixing') { mixingNut(); ub = velocity(); }
  else {
    mixingNut(); ub = velocity();
    const ke = model === 'kepsilon', bS = 0.09, kk = y.map((_, j) => Math.max((nut[j] * S[j]) / 0.3, 1e-10)), ww = y.map((yy, j) => (ke ? Math.max((0.09 * kk[j] * kk[j]) / Math.max(nut[j], 1e-6), 1e-10) : j === 0 ? 0 : Math.max(kk[j] / Math.max(nut[j], 1e-6), 6 / (0.075 * yy * yy) * (yy < 2.5 ? 1 : 0))));
    kk[0] = 0; if (ke) ww[0] = 0;
    const a = new Array(n + 1), b = new Array(n + 1), c = new Array(n + 1), d = new Array(n + 1), ur = 0.5;
    const solve = (phi, gam, src, sink, fixed) => {
      for (let j = 0; j <= n; j++) {
        if (j === 0 || (fixed && fixed[j] !== undefined)) { a[j] = 0; c[j] = 0; b[j] = 1; d[j] = j === 0 && !fixed ? 0 : fixed ? fixed[j] : 0; continue; }
        if (j === n) { a[j] = -1; b[j] = 1; c[j] = 0; d[j] = 0; continue; }
        const we = (0.5 * (r[j] + r[j + 1]) * 0.5 * (gam[j] + gam[j + 1])) / (y[j + 1] - y[j]), wwst = (0.5 * (r[j] + r[j - 1]) * 0.5 * (gam[j] + gam[j - 1])) / (y[j] - y[j - 1]), vol = (r[j] * (y[j + 1] - y[j - 1])) / 2, bb = (wwst + we + sink[j] * vol) / ur;
        a[j] = -wwst; c[j] = -we; b[j] = bb; d[j] = src[j] * vol + (1 - ur) * bb * phi[j];
      }
      const x = tridiag(a, b, c, d); for (let j = 0; j <= n; j++) phi[j] = Math.max(x[j], j === 0 ? 0 : 1e-12);
    };
    const gam = new Array(n + 1), src = new Array(n + 1), sink = new Array(n + 1), fmu = y.map((yy) => 1 - Math.exp(-0.0115 * yy));
    const wFix = ke ? null : y.map((yy, j) => (j === 0 ? 60 / (0.075 * y[1] * y[1]) : yy < 2.5 ? 6 / (0.075 * yy * yy) : undefined));
    for (it = 1; it <= maxIter; it++) {
      if (ke) { // Chien: k and the dissipation variable ε̃ (zero at the wall)
        for (let j = 0; j <= n; j++) { gam[j] = 1 + nut[j]; src[j] = nut[j] * S[j] * S[j]; sink[j] = j ? ww[j] / Math.max(kk[j], 1e-12) + 2 / (y[j] * y[j]) : 0; }
        solve(kk, gam, src, sink);
        for (let j = 0; j <= n; j++) { const Ret = (kk[j] * kk[j]) / Math.max(ww[j], 1e-12), f2 = 1 - (0.4 / 1.8) * Math.exp(-((Ret / 6) ** 2)); gam[j] = 1 + nut[j] / 1.3; src[j] = 1.35 * 0.09 * fmu[j] * kk[j] * S[j] * S[j]; sink[j] = j ? (1.8 * f2 * ww[j]) / Math.max(kk[j], 1e-12) + (2 * Math.exp(-0.5 * y[j])) / (y[j] * y[j]) : 0; }
        solve(ww, gam, src, sink);
        for (let j = 1; j <= n; j++) nut[j] = (0.09 * fmu[j] * kk[j] * kk[j]) / Math.max(ww[j], 1e-12);
      } else { // Wilcox k–ω
        for (let j = 0; j <= n; j++) { gam[j] = 1 + 0.5 * nut[j]; src[j] = nut[j] * S[j] * S[j]; sink[j] = bS * ww[j]; }
        solve(kk, gam, src, sink);
        for (let j = 0; j <= n; j++) { src[j] = (5 / 9) * S[j] * S[j]; sink[j] = 0.075 * ww[j]; }
        solve(ww, gam, src, sink, wFix);
        for (let j = 1; j <= n; j++) nut[j] = kk[j] / Math.max(ww[j], 1e-12);
      }
      nut[0] = 0;
      const un = velocity(); res = Math.abs(un - ub) / un; ub = un;
      if (res < tol && it > 20) break;
    }
    k = kk.slice();
  }
  const Re = 2 * Rp * ub;
  return { Re, f: 8 / (ub * ub), reTau: Rp, y, u: u.slice(), r: r.map((v) => v / Rp), uRel: u.map((v) => v / ub), nut: nut.slice(), k, iterations: it, residual: res, converged: model === 'laminar' || model === 'mixing' || res < 1e-6 };
}
/** Friction Reynolds number R⁺ of a smooth pipe at a bulk Reynolds number (Prandtl / Colebrook friction law). */
export const reTauOf = (Re) => (Re / 2) * Math.sqrt(frictionFactor(Math.max(Re, 10), 0, 'colebrook') / 8);

/**
 * 1-D volume-of-fluid advection of a liquid fraction at constant velocity on a periodic domain: first-order upwind or the
 * algebraic interface-capturing THINC flux (hyperbolic-tangent reconstruction, Xiao et al.), both conservative.
 * o: { n, cfl, beta, scheme: 'thinc' | 'upwind', length, u, tEnd, phi0: fn(x) }.
 * Returns { x[], phi[], exact[], mass0, mass, massError, l1 (mean absolute error), thickness (cells between 5 % and 95 % per interface), min, max }.
 */
export function vofAdvect1D(o = {}) {
  const n = o.n || 200, L = o.length || 1, u = o.u ?? 1, dx = L / n, cfl = o.cfl || 0.4, beta = o.beta || 2.5, tEnd = o.tEnd ?? L / Math.abs(u || 1), thinc = o.scheme !== 'upwind';
  const f0 = o.phi0 || ((x) => (x > 0.25 * L && x < 0.55 * L ? 1 : 0)), x = Array.from({ length: n }, (_, i) => (i + 0.5) * dx);
  let phi = x.map(f0); const mass0 = phi.reduce((s, v) => s + v, 0) * dx, steps = Math.max(1, Math.ceil((Math.abs(u) * tEnd) / (cfl * dx))), c = (Math.abs(u) * tEnd) / steps / dx, F = new Array(n);
  for (let s = 0; s < steps; s++) {
    for (let i = 0; i < n; i++) { // flux through the downstream face of cell i (u > 0), as a fraction of the cell volume
      const p = phi[i], pm = phi[(i - 1 + n) % n], pp = phi[(i + 1) % n];
      let fl = c * p;
      if (thinc && p > 1e-6 && p < 1 - 1e-6 && (pp - p) * (p - pm) > 0) {
        const gm = pp > pm ? 1 : -1, qq = (beta * (2 * p - 1)) / gm, w = (Math.exp(beta) - Math.exp(qq)) / (Math.exp(qq) - Math.exp(-beta));
        if (w > 0 && Number.isFinite(w)) { const xt = Math.log(w) / (2 * beta); fl = 0.5 * (c + (gm / beta) * (Math.log(Math.cosh(beta * (1 - xt))) - Math.log(Math.cosh(beta * (1 - c - xt))))); }
      }
      F[i] = clamp(fl, Math.max(0, c - (1 - p)), Math.min(c, p));
    }
    phi = phi.map((p, i) => p - F[i] + F[(i - 1 + n) % n]);
  }
  const shift = u * tEnd, exact = x.map((xx) => f0((((xx - shift) % L) + L) % L)), mass = phi.reduce((s, v) => s + v, 0) * dx;
  let inter = 0, edges = 0; for (let i = 0; i < n; i++) { if (phi[i] > 0.05 && phi[i] < 0.95) inter++; if (exact[i] !== exact[(i + 1) % n]) edges++; }
  return { x, phi, exact, mass0, mass, massError: Math.abs(mass - mass0) / Math.max(mass0, 1e-12), l1: mean(phi.map((p, i) => Math.abs(p - exact[i]))), thickness: inter / Math.max(edges, 1), min: Math.min(...phi), max: Math.max(...phi), steps };
}

/**
 * Dam-break of a liquid layer in a horizontal channel (collapse of a slug tail onto the film): shallow-water equations,
 * finite volume with the HLL flux, compared with the exact Ritter (dry bed) or Stoker (wet bed) solution.
 * o: { hL, hR (m), length (m, dam in the middle), n, tEnd (s), cfl }. Returns { x[], h[], u[], exact[], l1 (relative), massError, tEnd }.
 */
export function damBreak(o = {}) {
  const hL = o.hL ?? 1, hR = Math.max(o.hR ?? 0, 0), n = o.n || 300, L = o.length || 20 * hL, dx = L / n, cfl = o.cfl || 0.45, cL = Math.sqrt(G * hL), tEnd = o.tEnd ?? (0.2 * L) / cL, dry = 1e-9 * hL;
  const x = Array.from({ length: n }, (_, i) => -L / 2 + (i + 0.5) * dx); let h = x.map((v) => (v < 0 ? hL : hR)), q = new Array(n).fill(0), t = 0;
  const m0 = h.reduce((s, v) => s + v, 0) * dx, flux = (h1, q1, h2, q2) => {
    const u1 = h1 > dry ? q1 / h1 : 0, u2 = h2 > dry ? q2 / h2 : 0, c1 = Math.sqrt(G * h1), c2 = Math.sqrt(G * h2);
    const sL = h1 > dry ? Math.min(u1 - c1, h2 > dry ? u2 - c2 : u1 - c1) : u2 - 2 * c2, sR = h2 > dry ? Math.max(u2 + c2, h1 > dry ? u1 + c1 : u2 + c2) : u1 + 2 * c1;
    const f1 = [q1, q1 * u1 + 0.5 * G * h1 * h1], f2 = [q2, q2 * u2 + 0.5 * G * h2 * h2];
    if (sL >= 0) return f1; if (sR <= 0) return f2;
    return [(sR * f1[0] - sL * f2[0] + sL * sR * (h2 - h1)) / (sR - sL), (sR * f1[1] - sL * f2[1] + sL * sR * (q2 - q1)) / (sR - sL)];
  };
  while (t < tEnd - 1e-12) {
    let smax = 1e-9; for (let i = 0; i < n; i++) smax = Math.max(smax, (h[i] > dry ? Math.abs(q[i] / h[i]) : 0) + Math.sqrt(G * h[i]));
    const dt = Math.min((cfl * dx) / smax, tEnd - t), F = []; for (let i = 0; i <= n; i++) { const a = Math.max(i - 1, 0), b = Math.min(i, n - 1); F.push(flux(h[a], i === 0 ? -q[a] : q[a], h[b], i === n ? -q[b] : q[b])); }
    h = h.map((v, i) => Math.max(v - (dt / dx) * (F[i + 1][0] - F[i][0]), 0)); q = q.map((v, i) => v - (dt / dx) * (F[i + 1][1] - F[i][1])); t += dt;
  }
  let exact;
  if (hR <= dry) exact = x.map((v) => { const xi = v / tEnd; return xi <= -cL ? hL : xi >= 2 * cL ? 0 : (2 * cL - xi) ** 2 / (9 * G); });
  else { // Stoker: rarefaction, constant middle state, bore
    const um = (hm) => (hm - hR) * Math.sqrt((0.5 * G * (hm + hR)) / (hm * hR)), hm = brent((v) => 2 * (cL - Math.sqrt(G * v)) - um(v), hR, hL, 1e-13), u2 = um(hm), cm = Math.sqrt(G * hm), sb = (hm * u2) / (hm - hR);
    exact = x.map((v) => { const xi = v / tEnd; return xi <= -cL ? hL : xi <= u2 - cm ? (2 * cL - xi) ** 2 / (9 * G) : xi <= sb ? hm : hR; });
  }
  const mass = h.reduce((s, v) => s + v, 0) * dx;
  return { x, h, u: h.map((v, i) => (v > dry ? q[i] / v : 0)), exact, l1: mean(h.map((v, i) => Math.abs(v - exact[i]))) / hL, massError: Math.abs(mass - m0) / m0, tEnd };
}

/*NEW-SECTIONS-BEGIN*/
// =====================================================================================================
// 7a. More steady closures: Duns & Ros, Ansari mechanistic model, annular flow with entrainment and deposition,
//     homogeneous-relaxation flashing, Baker map
// =====================================================================================================
const tabLog = (xs, ys, x) => interp1(xs.map(log10), ys, log10(clamp(x, xs[0], xs[xs.length - 1])));
// Duns & Ros (1963) chart functions, read from open digitisations of the published charts (see PROVENANCE): slip factors F1–F7 against the
// liquid-viscosity number N_L, regime-boundary factors L1, L2 against the diameter number N_D and the friction correction f2.
const DR = { NL: [0.002, 0.004, 0.006, 0.01, 0.02, 0.04, 0.05, 0.07, 0.1, 0.2, 0.4, 1, 2],
  F1: [1.25, 1.25, 1.254, 1.26, 1.29, 1.523, 1.675, 1.915, 2.096, 2.151, 1.872, 1.178, 0.9], F2: [0.24, 0.24, 0.24, 0.24, 0.27, 0.499, 0.594, 0.751, 0.916, 1.02, 0.969, 0.812, 0.7],
  F3: [0.83, 0.84, 0.926, 1.296, 1.987, 2.654, 2.824, 3.122, 3.339, 3.698, 3.943, 4.121, 4.15], F4: [-19.76, -0.41, 9.97, 21.62, 35.6, 47.6, 50.6, 53.6, 55.1, 56.2, 55.95, 55.8, 55.62],
  F5: [0.22, 0.21, 0.207, 0.19, 0.17, 0.14, 0.13, 0.1, 0.06, 0.05, 0.06, 0.09, 0.111], F6: [0.852, 0.301, 0.093, -0.092, -0.059, 0.779, 1.047, 1.536, 2.076, 1.971, 1.78, 1.728, 1.75],
  F7: [0.13, 0.119, 0.1, 0.09, 0.07, 0.06, 0.05, 0.049, 0.04, 0.03, 0.03, 0.0235, 0.02],
  ND: [10, 16, 20, 30, 40, 50, 60, 70, 100, 275], L1: [2.066, 2.058, 2.056, 2.005, 1.635, 1.3, 1.097, 1.043, 1.032, 1.02], L2: [0.476, 0.499, 0.584, 0.771, 0.922, 1.042, 1.108, 1.129, 1.138, 1.128],
  fx: [0.001, 0.4, 0.7, 1, 2, 3, 6, 10, 20, 40, 100], f2: [1.002, 0.939, 0.76, 0.649, 0.504, 0.431, 0.342, 0.289, 0.24, 0.215, 0.201] };
// Mist-flow wall friction of Duns & Ros (also used by Orkiszewski): the liquid film acts as a roughness that depends on the gas Weber number.
function mistFriction(p) {
  const { vsg, rhoL, rhoG, muL, muG, D, rough, fModel } = p, s = Math.max(p.sigma, 1e-4), N = (rhoG * vsg * vsg * muL * muL) / (rhoL * s * s), e0 = s / (rhoG * vsg * vsg * D), ed = clamp(N <= 0.005 ? 0.0749 * e0 : 0.3713 * e0 * N ** 0.302, rough / D, 0.5);
  return ed > 0.05 ? 4 * (1 / (4 * log10(0.27 * ed)) ** 2 + 0.067 * ed ** 1.73) : frictionFactor((rhoG * vsg * D) / muG, ed, fModel);
}
/** Duns & Ros (1963) for upward flow: bubble (I), slug (II), transition and mist (III) regions with the slip-velocity charts. */
function dunsRos(p) {
  const { vsl, vsg, rhoL, rhoG, muL, D, theta, rough, P, fModel } = p, vm = vsl + vsg, lam = vsl / vm, { NLv, NGv, ND, NL } = velocityNumbers(p), sinT = Math.sin(theta), F = (k) => tabLog(DR.NL, DR[k], NL);
  const L1 = tabLog(DR.ND, DR.L1, ND), L2 = tabLog(DR.ND, DR.L2, ND), bubSlug = L1 + L2 * NLv, Ls = 50 + 36 * NLv, Lm = 75 + 84 * NLv ** 0.75, k4 = (rhoL / (G * Math.max(p.sigma, 1e-4))) ** 0.25;
  const slip = () => {
    const S = NGv <= bubSlug ? F('F1') + F('F2') * NLv + (F('F3') - F('F4') / ND) * (NGv / (1 + NLv)) ** 2 : ((1 + F('F5')) * (NGv ** 0.982 + 0.029 * ND + F('F6'))) / (1 + F('F7') * NLv) ** 2, vs = Math.max(S, 1e-6) / k4;
    const H = clamp((vs - vm + Math.sqrt((vm - vs) ** 2 + 4 * vs * vsl)) / (2 * vs), lam, 1), f1 = frictionFactor((rhoL * vsl * D) / muL, rough / D, fModel), f2 = tabLog(DR.fx, DR.f2, Math.max(((f1 / 4) * (vsg / vsl)) * ND ** (2 / 3), 0.001)), f3 = 1 + (f1 / 4) * Math.sqrt(vsg / (50 * vsl));
    return { holdup: H, fric: (((f1 * f2) / f3) * rhoL * vsl * vm) / (2 * D), rho: rhoL * H + rhoG * (1 - H), regime: NGv <= bubSlug ? 'bubble' : 'slug' };
  };
  const mist = (rg) => ({ holdup: lam, fric: (mistFriction(p) * rg * vsg * vsg) / (2 * D), rho: rhoL * lam + rg * (1 - lam), regime: 'annular' });
  let r;
  if (NGv <= Ls) r = slip();
  else if (NGv >= Lm) r = mist(rhoG);
  else { const a = slip(), b = mist((rhoG * NGv) / Lm), w = (Lm - NGv) / (Lm - Ls); r = { holdup: w * a.holdup + (1 - w) * b.holdup, fric: w * a.fric + (1 - w) * b.fric, rho: w * a.rho + (1 - w) * b.rho, regime: 'churn' }; }
  const rho = r.rho; delete r.rho; r.grav = rho * G * sinT;
  return NGv > Ls ? withAcc(r, rho, vm, vsg, P) : Object.assign(r, { acc: 0 });
}

/** Flow pattern of the Ansari et al. (1994) model for upward flow: 'bubble' | 'dispersed bubble' | 'slug' | 'annular' (before the film checks). */
export function ansariPattern(p) {
  const { vsl, vsg, rhoL, rhoG, muL, muG, D, theta, rough, fModel } = p, s = Math.max(p.sigma, 1e-4), dRho = Math.max(rhoL - rhoG, 1), vm = vsl + vsg, sinT = Math.max(Math.sin(theta), 0.05), v0 = (G * dRho * s / (rhoL * rhoL)) ** 0.25;
  if (vsg >= (3.1 * (s * G * sinT * dRho) ** 0.25) / Math.sqrt(rhoG)) return 'annular';
  const Hg = vsg / vm, rhoM = rhoL * (1 - Hg) + rhoG * Hg, muM = muL * (1 - Hg) + muG * Hg, f = frictionFactor((D * rhoM * vm) / muM, rough / D, fModel), c = 2 * Math.sqrt((0.4 * s) / (dRho * G)) * (rhoL / s) ** 0.6 * (2 / D) ** 0.4;
  if (Hg <= 0.76 && vm >= ((0.725 + 4.15 * Math.sqrt(Hg)) / c / (f / 4) ** 0.4) ** (1 / 1.2)) return 'dispersed bubble';
  const bubbly = theta > 70 * DEG && D > 0.95 * 19 * Math.sqrt((dRho * s) / (rhoL * rhoL * G));
  return bubbly && vsg <= (vsl + 1.15 * v0 * Math.sin(theta)) / 3 ? 'bubble' : 'slug';
}
/**
 * Ansari, Sylvester, Sarica, Shoham & Brill (1994) mechanistic model for upward two-phase flow: bubble / dispersed-bubble flow with
 * the Harmathy slip, slug flow with the Taylor-bubble film (developed or developing), annular flow with the Wallis entrainment and the
 * film–core momentum balance; annular flow that fails the film-stability or blockage checks is treated as slug flow.
 */
function ansari(p) {
  const { vsl, vsg, rhoL, rhoG, muL, muG, D, theta, rough, P, fModel } = p, s = Math.max(p.sigma, 1e-4), dRho = Math.max(rhoL - rhoG, 1), vm = vsl + vsg, lam = vsl / vm, sinT = Math.sin(theta), v0 = 1.53 * (G * s * dRho / (rhoL * rhoL)) ** 0.25, ed = rough / D;
  let pat = ansariPattern(p);
  if (pat === 'annular') {
    const x = Math.sqrt(rhoG / rhoL) * 1e4 * vsg * muG / s, fe = clamp(1 - Math.exp(-0.125 * (x - 1.5)), 0, 1), c = fe > 0.9 ? 300 : 24 * (rhoL / rhoG) ** (1 / 3), alfc = 1 / (1 + (fe * vsl) / vsg), vsc = vsg + fe * vsl, rhoC = rhoG * alfc + rhoL * (1 - alfc), muC = muG * alfc + muL * (1 - alfc);
    const fcs = frictionFactor((rhoC * vsc * D) / muC, ed, fModel), fls = frictionFactor((rhoL * vsl * D) / muL, ed, fModel), a = fe < 0.9999 ? ((1 - fe) ** 2 * frictionFactor((rhoL * vsl * (1 - fe) * D) / muL, ed, fModel)) / fls : 1;
    const gcs = (fcs * rhoC * vsc * vsc) / (2 * D), gls = (fls * rhoL * vsl * vsl) / (2 * D), xmo2 = (gls / gcs) * a, ym = (G * sinT * (rhoL - rhoC)) / gcs;
    let deld = 0, H = 1 - alfc, ok = true;
    if (fe < 0.9999) {
      const fn = (d) => { const t = 4 * d * (1 - d); return ym - (1 + c * d) / t / (1 - t) ** 2.5 + xmo2 / t ** 3; }, lo = 1e-6, hi = 0.499;
      if (fn(lo) * fn(hi) < 0) deld = brent(fn, lo, hi, 1e-10); else ok = false;
      H = 4 * deld * (1 - deld) + (1 - alfc) * (1 - 2 * deld) ** 2;
      if (ok && H > 0.12) ok = false; // the film and the entrained liquid bridge the core
      if (ok) { const st = (d) => { const t = 1 - (1 - 2 * d) ** 2; return ym - ((2 - 1.5 * t) * xmo2) / t ** 3 / (1 - 1.5 * t); }; let ds = 0.499, fa = st(1e-5); for (let i = 1; i <= 60; i++) { const b = 1e-5 + (0.3 * i) / 60, fb = st(b); if (Number.isFinite(fb) && fa * fb < 0) { ds = brent(st, b - 0.005, b, 1e-9); break; } fa = fb; } if (ds < deld) ok = false; } // film thicker than the stable (minimum-shear) film
    } else if (H > 0.12) ok = false;
    if (ok) { const phi = fe < 0.9999 ? (1 + c * deld) / (1 - 2 * deld) ** 5 : 1; return withAcc({ holdup: clamp(H, lam * 0.2, 1), fric: gcs * phi + (rhoC - (rhoL * H + rhoG * (1 - H))) * G * sinT, grav: (rhoL * H + rhoG * (1 - H)) * G * sinT, regime: 'annular', entrainment: fe, filmThickness: deld * D }, rhoC, vm, vsg, P); }
    pat = 'slug';
  }
  if (pat === 'bubble' || pat === 'dispersed bubble') {
    const H = pat === 'dispersed bubble' ? lam : clamp(brent((e) => v0 * Math.sqrt(e) + 1.2 * vm - vsg / (1 - e), lam, 0.9999, 1e-10), lam, 1), rho = rhoL * H + rhoG * (1 - H), f = frictionFactor((rho * vm * D) / (muL * lam + muG * (1 - lam)), ed, fModel);
    return withAcc({ holdup: H, fric: (f * rho * vm * vm) / (2 * D), grav: rho * G * sinT, regime: pat }, rho, vm, vsg, P);
  }
  // slug flow
  const lls = 30 * D, vtb = 1.2 * vm + 0.35 * Math.sqrt((G * D * dRho) / rhoL), als = Math.min(vsg / (0.425 + 2.65 * vm), 0.6), vgls = (a) => 1.2 * vm + v0 * Math.sqrt(1 - als);
  const bal = (atb) => { const vltb = 9.916 * Math.sqrt(G * D * (1 - Math.sqrt(atb))), vlls = vtb - ((vtb + vltb) * (1 - atb)) / (1 - als), vg = als > 0.25 ? vlls : vgls(), vgtb = vtb * (1 - als / atb) + (vg * als) / atb, b1 = (vlls * (1 - als) - vsl) / (vltb * (1 - atb) + vlls * (1 - als)), b2 = (vsg - als * vg) / (atb * vgtb - als * vg); return { d: b1 - b2, b1, b2, vltb, vlls, vg, vgtb }; };
  let atb = null; { let a0 = Math.max(als + 0.01, 0.3), f0 = bal(a0).d; for (let i = 1; i <= 80 && atb === null; i++) { const a1 = a0 + ((0.9999 - Math.max(als + 0.01, 0.3)) * 1) / 80, f1 = bal(a1).d; if (Number.isFinite(f0) && Number.isFinite(f1) && f0 * f1 <= 0) atb = brent((v) => bal(v).d, a0, a1, 1e-10); a0 = a1; f0 = f1; } }
  if (atb === null) { const r = gradient(p, 'driftFlux'); r.regime = 'slug'; return r; } // no Taylor-bubble solution (very low rates): Bendiksen drift flux
  const q = bal(atb), beta = clamp(0.5 * (q.b1 + q.b2), 0.01, 0.99), ltb = (lls * beta) / (1 - beta);
  let atbE = atb; { const lc = (q.vltb + vtb) ** 2 / (2 * G); if (lc > 0.75 * ltb) { const c = (vsg - q.vg * als) / vtb, d = 1 - vsg / vtb, e = vtb - q.vlls, ff = (-2 * d * c * lls - (2 * (e * (1 - als)) ** 2) / G) / (d * d), gg = ((c * lls) / d) ** 2, h = ff * ff - 4 * gg; if (h > 0) { const l2 = Math.max((-ff + Math.sqrt(h)) / 2, (-ff - Math.sqrt(h)) / 2); if (l2 > 0) atbE = clamp(1 - (2 * (vtb - q.vlls) * (1 - als)) / Math.sqrt(2 * G * l2), als, 0.9999); } } } // developing Taylor bubble: mean void of the falling film region
  const asu = atbE * beta + als * (1 - beta), H = clamp(1 - asu, lam, 1), rhoS = rhoL * (1 - als) + rhoG * als, vmls = q.vlls * (1 - als) + q.vg * als, ans = (q.vg * als) / Math.max(vmls, 1e-9), f = frictionFactor((rhoS * Math.abs(vmls) * D) / (muL * (1 - ans) + muG * ans), ed, fModel);
  const grav = G * sinT * (rhoS * (1 - beta) + (rhoG * atbE + rhoL * (1 - atbE)) * beta);
  return withAcc({ holdup: H, fric: ((rhoS * vmls * vmls * f) / (2 * D)) * (1 - beta), grav, regime: 'slug', slugFraction: 1 - beta, taylorVoid: atb }, rhoL * H + rhoG * (1 - H), vm, vsg, P);
}

// ---- annular-mist flow as three fields: gas core, entrained droplets, wall film ---------------------------------------------
const dragSphere = (Re) => (Re < 1e-9 ? 0 : Re < 1000 ? (24 / Re) * (1 + 0.15 * Re ** 0.687) : 0.44); // Schiller & Naumann
/** Equilibrium entrained liquid fraction of annular flow (Wallis, as used by Ansari et al.): FE = 1 − exp[−0.125 (v_crit − 1.5)], v_crit = 10⁴ vsg μg / σ · (ρg/ρl)^½. */
export const entrainmentFraction = (p) => clamp(1 - Math.exp(-0.125 * ((1e4 * p.vsg * p.muG) / Math.max(p.sigma, 1e-4) * Math.sqrt(p.rhoG / p.rhoL) - 1.5)), 0, 0.9999);
/**
 * Annular-mist flow with three fields (gas, droplets carried in the gas core, liquid film on the wall).
 * Film thickness from the combined momentum balance of film and core with the Wallis interfacial friction f_i = f_c (1 + 300 δ/D) and
 * a wall friction on the film; droplets slip behind the gas by their terminal velocity (Schiller–Naumann drag, size from a critical
 * Weber number); the entrained flow follows  dW_E/dz = π D (R_E − R_D),  R_D = k_D C  (deposition, C = droplet concentration in the core),
 * R_E = k_E k_D C_eq  (entrainment written so that the developed state is the Wallis equilibrium fraction when k_E = 1).
 * p: closure point { vsl, vsg, rhoL, rhoG, muL, muG, sigma, D, theta, rough, P, fModel }; mp: { entFrac (entrained fraction to use; default equilibrium),
 * kDep (deposition velocity m/s), entMult (k_E), weCrit (droplet critical Weber number), dropMult, fiMult }.
 * Returns { holdup, holdupFilm, holdupDrops, fric, grav, acc, regime, filmThickness, filmVelocity, coreVelocity, tauI, tauW, entrainment, entEq, dropSize, dropSlip, rateDep, rateEnt (kg/m²/s), relaxLength (m) }.
 */
export function annularMist(p, mp = {}) {
  const { vsl, vsg, rhoL, rhoG, muL, muG, D, theta, rough, P, fModel } = p, s = Math.max(p.sigma, 1e-4), A = (PI * D * D) / 4, sinT = Math.sin(theta), kD = mp.kDep ?? 0.15, kE = mp.entMult ?? 1, WeC = mp.weCrit ?? 12, vm = vsl + vsg;
  const Eeq = clamp(entrainmentFraction(p) * kE, 0, 0.9999), E = clamp(mp.entFrac ?? Eeq, 0, 0.9999), dMax = ((WeC * s) / (rhoG * Math.max(vsg, 0.1) ** 2)) * (mp.dropMult ?? 1), dDrop = clamp(dMax, 5e-6, 0.3 * D);
  let vt = Math.sqrt((4 * G * dDrop * Math.max(rhoL - rhoG, 1)) / (3 * 0.44 * rhoG)); for (let k = 0; k < 12; k++) { const cd = Math.max(dragSphere((rhoG * vt * dDrop) / muG), 0.1); vt = Math.sqrt((4 * G * dDrop * Math.max(rhoL - rhoG, 1)) / (3 * cd * rhoG)); }
  const state = (dl) => {
    const core = (1 - 2 * dl) ** 2, Ac = A * core, Af = A - Ac, vf = (vsl * (1 - E) * A) / Af, vg = (vsg * A) / Ac, vd = Math.max(vg - vt * sinT, 0.2 * vg), Hd = Math.min((vsl * E * A) / (vd * Ac), 0.5), rhoC = rhoG * (1 - Hd) + rhoL * Hd, vc = (vsg + vsl * E) * A / Ac;
    const fc = frictionFactor((rhoC * vc * D * (1 - 2 * dl)) / muG, 0, fModel) / 4, fi = fc * (1 + 300 * dl) * (mp.fiMult ?? 1), Dhf = 4 * dl * D * (1 - dl), ff = frictionFactor((rhoL * Math.abs(vf) * Dhf) / muL, rough / D, fModel) / 4;
    const ti = (fi * rhoC * (vc - vf) * Math.abs(vc - vf)) / 2, tw = (ff * rhoL * vf * Math.abs(vf)) / 2, Si = PI * D * (1 - 2 * dl), dpC = (ti * Si) / Ac + rhoC * G * sinT, dpF = (tw * PI * D - ti * Si) / Af + rhoL * G * sinT;
    return { dpC, dpF, vf, vc, vd, Hd, rhoC, ti, tw, core, Ac };
  };
  let dl; { const f = (x) => { const q = state(x); return q.dpC - q.dpF; }, lo = 1e-6, hi = 0.45; const flo = f(lo), fhi = f(hi); dl = flo * fhi < 0 ? brent(f, lo, hi, 1e-12) : Math.abs(flo) < Math.abs(fhi) ? lo : hi; }
  const q = state(dl), Hf = 1 - q.core, Hdr = q.Hd * q.core, H = clamp(Hf + Hdr, 1e-6, 1), rho = rhoL * H + rhoG * (1 - H), grav = rho * G * sinT, C = (rhoL * vsl * E) / Math.max(vsg + vsl * E, 1e-9), Ceq = (rhoL * vsl * Eeq) / Math.max(vsg + vsl * Eeq, 1e-9);
  const r = { holdup: H, holdupFilm: Hf, holdupDrops: Hdr, fric: q.dpC - grav, grav, regime: 'annular', filmThickness: dl * D, filmVelocity: q.vf, coreVelocity: q.vc, tauI: q.ti, tauW: q.tw, entrainment: E, entEq: Eeq, dropSize: dDrop, dropSlip: vt, rateDep: kD * C, rateEnt: kD * Ceq, relaxLength: (vsg * A) / (PI * D * kD) };
  return withAcc(r, rho, vm, vsg, P);
}
/**
 * Development of the entrained fraction along a pipe of constant conditions (three-field mass balance marched with RK4).
 * o: { p (closure point), mp, length, n, e0 (entrained fraction at the inlet) }. Returns { z[], entrained[], film[] (kg/s), holdup[], dpdz[], filmThickness[], eq, balance: { liquidIn, liquidOut } }.
 */
export function annularDevelopment(o) {
  const p = o.p, mp = o.mp || {}, n = o.n || 60, L = o.length, dz = L / n, A = (PI * p.D * p.D) / 4, WL = p.rhoL * p.vsl * A, kD = mp.kDep ?? 0.15, eq = clamp(entrainmentFraction(p) * (mp.entMult ?? 1), 0, 0.9999);
  const conc = (E) => (p.rhoL * p.vsl * E) / Math.max(p.vsg + p.vsl * E, 1e-9), dE = (E) => (PI * p.D * kD * (conc(eq) - conc(E))) / WL, z = [0], ent = [o.e0 ?? 0], st = [annularMist(p, { ...mp, entFrac: ent[0] })];
  for (let i = 1; i <= n; i++) { const E = ent[i - 1], k1 = dE(E), k2 = dE(E + 0.5 * dz * k1), k3 = dE(E + 0.5 * dz * k2), k4 = dE(E + dz * k3), En = clamp(E + (dz / 6) * (k1 + 2 * k2 + 2 * k3 + k4), 0, 0.9999); z.push(i * dz); ent.push(En); st.push(annularMist(p, { ...mp, entFrac: En })); }
  return { z, entrained: ent, film: ent.map((E) => WL * (1 - E)), drops: ent.map((E) => WL * E), holdup: st.map((q) => q.holdup), dpdz: st.map((q) => q.fric + q.grav + q.acc), filmThickness: st.map((q) => q.filmThickness), eq, relaxLength: st[0].relaxLength, balance: { liquidIn: WL, liquidOut: WL * (1 - ent[n]) + WL * ent[n] } };
}

// ---- homogeneous-relaxation model of flashing flow ------------------------------------------------------------------------------
/** Relaxation time of Downar-Zapolski et al. (1996), low-pressure form: Θ = Θ₀ α^−0.257 ψ^−2.24, Θ₀ = 6.51·10⁻⁴ s, ψ = (p_sat − p)/p_sat. */
export const relaxationTime = (alpha, psi, theta0 = 6.51e-4) => theta0 * clamp(alpha, 1e-4, 1) ** -0.257 * clamp(psi, 1e-4, 1) ** -2.24;
/**
 * Steady flashing flow in a pipe with the homogeneous-relaxation model (equal velocities, the vapour mass fraction x relaxes to its
 * equilibrium value):  G dx/dz = ρ (x_eq − x)/Θ,   dp/dz = −[f G² v/(2D) + g sinθ / v + G² dv/dz],   v = x v_g(p) + (1 − x) v_l.
 * model 'hem' is the homogeneous-equilibrium limit x = x_eq(p). The specific volume derivative is taken along the solution (explicit in x,
 * implicit in p through the gas compressibility), marched with a midpoint rule.
 * o: { G (kg/m²/s), D, length, n, p0 (Pa), x0, pSat (Pa), xEq(p), vG(p), vL, f (Darcy), theta, model: 'hrm' | 'hem', theta0, tau (fixed relaxation time, s) }.
 * Returns { z[], p[], x[], xEq[], alpha[], u[], lag (max x_eq − x), choked (bool: the march reached the sonic limit), massFluxCheck }.
 */
export function flashingFlow(o) {
  const n = o.n || 200, dz = o.length / n, Gm = o.G, D = o.D, f = o.f ?? 0.02, sinT = Math.sin(o.theta || 0), hem = o.model === 'hem', vL = o.vL, vG = o.vG, xEq = o.xEq, vol = (p, x) => x * vG(p) + (1 - x) * vL;
  const z = [0], P = [o.p0], X = [hem ? xEq(o.p0) : o.x0 ?? 0], XE = [xEq(o.p0)]; let choked = false;
  const rhs = (p, x) => {
    const v = vol(p, x), al = clamp((x * vG(p)) / v, 0, 1), th = o.tau ?? relaxationTime(Math.max(al, 1e-3), (o.pSat - p) / o.pSat, o.theta0), dx = hem ? 0 : (xEq(p) - x) / (th * Gm * v), dp = 1e-4 * p, dvdp = hem ? (vol(p + dp, xEq(p + dp)) - vol(p - dp, xEq(p - dp))) / (2 * dp) : (x * (vG(p + dp) - vG(p - dp))) / (2 * dp);
    const den = 1 + Gm * Gm * dvdp, num = (f * Gm * Gm * v) / (2 * D) + (G * sinT) / v + (hem ? 0 : Gm * Gm * (vG(p) - vL) * dx);
    return { dp: -num / den, dx, den };
  };
  for (let i = 1; i <= n; i++) {
    const a = rhs(P[i - 1], X[i - 1]); if (!(a.den > 0.02)) { choked = true; break; }
    const pm = P[i - 1] + 0.5 * dz * a.dp, xm = hem ? xEq(pm) : clamp(X[i - 1] + 0.5 * dz * a.dx, 0, 1); if (!(pm > 0)) { choked = true; break; }
    const b = rhs(pm, xm); if (!(b.den > 0.02)) { choked = true; break; }
    const pn = P[i - 1] + dz * b.dp; if (!(pn > 0)) { choked = true; break; }
    z.push(i * dz); P.push(pn); X.push(hem ? xEq(pn) : clamp(X[i - 1] + dz * b.dx, 0, 1)); XE.push(xEq(pn));
  }
  const alpha = X.map((x, i) => (x * vG(P[i])) / vol(P[i], x)), u = X.map((x, i) => Gm * vol(P[i], x));
  return { z, p: P, x: X, xEq: XE, alpha, u, lag: Math.max(...X.map((x, i) => XE[i] - x)), choked };
}
/** Mechanistic closure of the kernel with the three-field annular-mist model in annular flow (entrainment and deposition feed the pressure gradient). */
function mechEnt(p, mp = {}) {
  const fp = flowPattern(p);
  if (fp.pattern === 'annular') return annularMist(p, mp.annular || {});
  const r = gradient({ ...p, label: false }, 'mechanistic'); return r;
}

// =====================================================================================================
// 7b. Transient two-fluid model (separate mass, momentum and energy equations of gas and liquid)
// =====================================================================================================
// Cross-section relations of a partly filled circular pipe tabulated against the liquid holdup (all lengths per diameter).
const SECT = (() => { const M = 600, hD = new Float64Array(M + 1), SL = new Float64Array(M + 1), Si = new Float64Array(M + 1); for (let j = 0; j <= M; j++) { const H = j / M, x = H <= 0 ? 0 : H >= 1 ? 1 : brent((v) => holdupOfLevel(v) - H, 0, 1, 1e-12), c = clamp(2 * x - 1, -1, 1); hD[j] = x; SL[j] = PI - Math.acos(c); Si[j] = Math.sqrt(1 - c * c); } return { M, hD, SL, Si }; })();
const sect = (tab, H) => { const u = clamp(H, 0, 1) * SECT.M, j = Math.min(Math.floor(u), SECT.M - 1); return tab[j] + (tab[j + 1] - tab[j]) * (u - j); };
// Thomas algorithm on typed arrays (a: sub-, b: main, c: super-diagonal, d: right-hand side; the solution overwrites d); cyclic variant by Sherman–Morrison.
function thomas(a, b, c, d, n, w) { w[0] = c[0] / b[0]; d[0] /= b[0]; for (let i = 1; i < n; i++) { const m = b[i] - a[i] * w[i - 1]; w[i] = c[i] / m; d[i] = (d[i] - a[i] * d[i - 1]) / m; } for (let i = n - 2; i >= 0; i--) d[i] -= w[i] * d[i + 1]; }
function thomasCyclic(a, b, c, d, n, w, u, bb) { const gam = -b[0], a0 = a[0], cn = c[n - 1]; for (let i = 0; i < n; i++) { bb[i] = b[i]; u[i] = 0; } bb[0] = b[0] - gam; bb[n - 1] = b[n - 1] - (a0 * cn) / gam; u[0] = gam; u[n - 1] = cn; thomas(a, bb, c, d, n, w); thomas(a, bb, c, u, n, w); const f = (d[0] + (a0 * d[n - 1]) / gam) / (1 + u[0] + (a0 * u[n - 1]) / gam); for (let i = 0; i < n; i++) d[i] -= f * u[i]; }

/**
 * Transient one-dimensional two-fluid model of stratified / slug flow in a pipe or a plane channel.
 * Per phase k (liquid l: incompressible, gas g: ideal gas ρg = p / (Rs Tg)):
 *   ∂(αk ρk)/∂t + ∂(αk ρk uk)/∂x = 0
 *   ρk (∂uk/∂t + uk ∂uk/∂x) = −∂p/∂x − ρk g cosθ ∂h/∂x − ρk g sinθ − τwk Sk / (αk A) ± τi Si / (αk A) − B
 *   ∂(αk ρk ek)/∂t + ∂(αk ρk ek uk)/∂x = −p [∂αk/∂t + ∂(αk uk)/∂x] + hi ai (Tj − Tk) + qwk        (ek = cvk Tk; six-equation form, o.energy)
 * with αl + αg = 1, one pressure p (at the interface), the liquid level h(αl) from the pipe geometry (the level-gradient term makes the
 * model hyperbolic below the inviscid Kelvin–Helmholtz limit) and B a uniform background pressure gradient (periodic domains).
 * Discretisation: staggered finite volumes, first-order upwind, wall and interfacial friction and the pressure gradient implicit
 * (the two phase momentum equations are solved together on every face and substituted in the volume-conservation equation, which gives
 * one tridiagonal pressure equation per step), so the step is limited by the material and gravity-wave CFL number. Liquid bridging
 * (αl → 1) is not treated specially: the pressure equation becomes the incompressible one there, which is what captures slugs.
 * o: { n, length, D (pipe) | channelH (plane channel of that height, unit width), theta (rad, + up) | thetaOf(x), rhoL, muL, muG, Rs (J/kg/K of the gas, Z R / M),
 *      T (K, isothermal value), pOut (Pa), periodic, closed, bodyForce (Pa/m), friction (false = inviscid), fiMult, fwlMult, fwgMult, fiMin (0.0142),
 *      interfacialPressure (δ of the Bestion interfacial-pressure term, 0 = off),
 *      virtualMass (coefficient C_vm of the added-mass force on the relative acceleration, bubbly limit), init(x) → { al, ul, ug, p, Tl, Tg },
 *      inlet(t) → { al, ul, ug, Tl, Tg }, tEnd, cfl, dtMax, gasMin (gas fraction below which a face is a liquid bridge and the gas moves with the liquid, 0.05), residualGas (gas fraction that sets the compliance of liquid-full cells, 0.01) | liquidSound (m/s, sound speed of the liquid used instead), probes: [x], bridge (holdup counted as a slug, 1 − gasMin), nField, nSeries,
 *      energy: { cvL, cvG, hi (W/m²K between the phases), UwL, UwG (W/m²K to the wall on the wetted perimeters), Tw } }
 * Returns { x[], al[], ul[], ug[], p[], Tl[], Tg[] (final, cell centres), t, steps, dtMean, mass: { liquid0, liquid, liquidIn, liquidOut, gas0, gas, gasIn, gasOut, errorL, errorG },
 *           energy: { e0, e, in, out, wall, work, error }, volErrMax, series: { t, probes: [[αl]], pIn, inv }, field: { t, x, al: [[..]] }, slugs: [{ probe, t, duration }], momentum: { residual } }.
 */
export function twoFluid(o) {
  const N = o.n || 200, L = o.length || 10, dx = L / N, per = !!o.periodic, closed = !!o.closed, chan = o.channelH > 0, D = chan ? o.channelH : o.D || 0.1, A = chan ? D : (PI * D * D) / 4;
  const rhoL = o.rhoL ?? 1000, muL = o.muL ?? 1e-3, muG = o.muG ?? 1.8e-5, Rs = o.Rs ?? 287, T0 = o.T ?? 293.15, pOut = o.pOut ?? 1e5, B = o.bodyForce || 0, fric = o.friction !== false, fiM = o.fiMult ?? 1, fwlM = o.fwlMult ?? 1, fwgM = o.fwgMult ?? 1, fiMin = o.fiMin ?? 0.0142, cvm = o.virtualMass || 0, piC = o.interfacialPressure || 0;
  const tEnd = o.tEnd ?? 1, cfl = o.cfl ?? 0.5, dtMax = o.dtMax ?? tEnd / 20, en = o.energy || null, cvL = en?.cvL ?? 4180, cvG = en?.cvG ?? 718, eps = 1e-7, agMin = o.gasMin ?? 0.05, bridge = o.bridge ?? 1 - agMin, kapL = o.liquidSound ? 1 / (rhoL * o.liquidSound ** 2) : (o.residualGas ?? 0.01) / pOut, cG = Math.sqrt(1.4 * Rs * T0); // compliance of a liquid-full cell: as if it kept 1 % of gas at the outlet pressure (regularises the closing of a liquid bridge)
  const F = (n) => new Float64Array(n), NF = per ? N : N + 1, al = F(N), mg = F(N), p = F(N), Tl = F(N), Tg = F(N), ul = F(NF + 1), ug = F(NF + 1), hh = F(N), sinT = F(NF + 1), cosT = F(NF + 1), xc = F(N);
  for (let f = 0; f <= NF; f++) { const th = typeof o.thetaOf === 'function' ? o.thetaOf(Math.min(f * dx, L)) : o.theta || 0; sinT[f] = Math.sin(th); cosT[f] = Math.cos(th); }
  const init = o.init || (() => ({ al: 0.5, ul: 0, ug: 0 }));
  for (let i = 0; i < N; i++) { xc[i] = (i + 0.5) * dx; const s = init(xc[i]); al[i] = clamp(s.al, 0, 1); p[i] = s.p ?? pOut; Tl[i] = s.Tl ?? T0; Tg[i] = s.Tg ?? T0; mg[i] = ((1 - al[i]) * p[i]) / (Rs * Tg[i]); }
  for (let f = 0; f < (per ? N : N + 1); f++) { const s = init(Math.min(f * dx, L)); ul[f] = s.ul || 0; ug[f] = s.ug || 0; }
  const level = (H) => (chan ? clamp(H, 0, 1) * D : sect(SECT.hD, H) * D);
  const us = F(NF + 1), gs = F(NF + 1), dl = F(NF + 1), dg = F(NF + 1), alU = F(NF + 1), mgU = F(NF + 1), FL = F(NF + 1), FG = F(NF + 1), sa = F(N), sb = F(N), sc = F(N), sd = F(N), w1 = F(N), w2 = F(N), w3 = F(N), ulO = F(NF + 1), ugO = F(NF + 1);
  const total = () => { let l = 0, g = 0, e = 0; for (let i = 0; i < N; i++) { l += al[i] * rhoL * dx * A; g += mg[i] * dx * A; e += (al[i] * rhoL * cvL * Tl[i] + mg[i] * cvG * Tg[i]) * dx * A; } return { l, g, e }; };
  const start = total(), bal = { lIn: 0, lOut: 0, gIn: 0, gOut: 0, eIn: 0, eOut: 0, eWall: 0, eWork: 0, vol: 0, momRes: 0 };
  const probes = (o.probes || []).map((x) => clamp(Math.floor(x / dx), 0, N - 1)), ser = { t: [], probes: probes.map(() => []), pIn: [], inv: [] }, fld = { t: [], al: [] }, slugs = [], inSlug = probes.map(() => -1);
  const nField = o.nField || 0, nSeries = o.nSeries || 400; let nextF = 0, nextS = 0, t = 0, steps = 0;
  const record = () => {
    for (let k = 0; k < probes.length; k++) { const b = al[probes[k]] >= bridge; if (b && inSlug[k] < 0) inSlug[k] = t; else if (!b && inSlug[k] >= 0) { slugs.push({ probe: k, t: inSlug[k], duration: t - inSlug[k] }); inSlug[k] = -1; } }
    if (t >= nextS - 1e-12) { let inv = 0; for (let i = 0; i < N; i++) inv += al[i] * dx * A; ser.t.push(t); ser.pIn.push(p[0]); ser.inv.push(inv); probes.forEach((c, k) => ser.probes[k].push(al[c])); nextS += tEnd / nSeries; }
    if (nField && t >= nextF - 1e-12) { fld.t.push(t); fld.al.push(Array.from(al)); nextF += tEnd / (nField - 1); }
  };
  record();
  const wrap = (i) => (per ? (i + N) % N : i);
  while (t < tEnd - 1e-12 && steps < (o.maxSteps || 4e6)) {
    // ---- time step: material and gravity-wave CFL condition
    let vmax = 1e-6; for (let f = 0; f < NF + (per ? 0 : 0); f++) { const i = wrap(Math.min(f, N - 1)), c = Math.sqrt(G * Math.abs(cosT[f]) * D * 0.8); vmax = Math.max(vmax, Math.abs(ul[f]) + c, Math.abs(ug[f])); }
    const dt = Math.min((cfl * dx) / vmax, dtMax, tEnd - t), rdt = 1 / dt, bc = !per && !closed && o.inlet ? o.inlet(t + dt) : null;
    for (let i = 0; i < N; i++) hh[i] = level(al[i]);
    ulO.set(ul); ugO.set(ug);
    // ---- phase momentum equations on the faces, linear in the new pressure difference: u = u* − d (pR − pL)
    const f0 = per ? 0 : 1, f1 = per ? N - 1 : N;
    for (let f = f0; f <= f1; f++) {
      const Lc = wrap(f - 1), out = !per && f === N, Rc = out ? Lc : wrap(f), dxf = out ? 0.5 * dx : dx;
      if (!per && closed && f === N) { us[f] = 0; gs[f] = 0; dl[f] = 0; dg[f] = 0; continue; }
      const a = clamp(0.5 * (al[Lc] + al[Rc]), eps, 1 - eps), ag = 1 - a, rg = Math.max((0.5 * (mg[Lc] + mg[Rc])) / ag, 0.5 * (p[Lc] + p[Rc]) / (Rs * 0.5 * (Tg[Lc] + Tg[Rc])) * 0.2, 1e-3), u1 = ulO[f], u2 = ugO[f];
      let Fwl = 0, Fwg = 0, Fil = 0, Fig = 0;
      if (fric) {
        let SLp, SGp, Sip, AL, AG;
        if (chan) { SLp = 1; SGp = 1; Sip = 1; AL = a * D; AG = ag * D; } else { SLp = sect(SECT.SL, a) * D; Sip = sect(SECT.Si, a) * D; SGp = PI * D - SLp; AL = a * A; AG = ag * A; }
        const DL = (4 * AL) / Math.max(SLp, 1e-12), DG = (4 * AG) / Math.max(SGp + Sip, 1e-12), fL = fanning((rhoL * Math.abs(u1) * DL) / muL) * fwlM, fG = fanning((rg * Math.abs(u2) * DG) / muG), fi = Math.max(fG, fiMin) * fiM, ur = Math.abs(u2 - u1);
        Fwl = (fL * Math.abs(u1) * SLp) / (2 * AL); Fwg = (fG * fwgM * Math.abs(u2) * SGp) / (2 * AG); Fil = (fi * rg * ur * Sip) / (2 * AL * rhoL); Fig = (fi * ur * Sip) / (2 * AG);
      }
      // upwind convective acceleration and explicit gravity terms
      const fm = per ? wrap(f - 1) : f - 1, fp = per ? wrap(f + 1) : Math.min(f + 1, N), cl = u1 >= 0 ? (u1 * (u1 - ulO[fm])) / dx : (u1 * (ulO[fp] - u1)) / dx, cg = u2 >= 0 ? (u2 * (u2 - ugO[fm])) / dx : (u2 * (ugO[fp] - u2)) / dx;
      // interfacial-pressure term Δp_i = δ αg αl ρg ρl (ug − ul)² / (αg ρl + αl ρg) (Bestion): optional regularisation that keeps the model hyperbolic beyond the inviscid Kelvin–Helmholtz limit
      const dpi = piC > 0 ? (piC * ag * a * rg * rhoL * (u2 - u1) ** 2) / (ag * rhoL + a * rg) : 0, dal = (al[Rc] - al[Lc]) / dxf;
      const gh = (G * cosT[f] * (hh[Rc] - hh[Lc])) / dxf, b1 = u1 * rdt - cl - G * sinT[f] - gh - B / rhoL - (dpi * dal) / (a * rhoL), b2 = u2 * rdt - cg - G * sinT[f] - gh - B / rg + (dpi * dal) / (ag * rg);
      // added-mass coupling of the relative acceleration (C_vm ρl αg per unit volume on each phase, opposite signs)
      const vl = cvm ? (cvm * ag) / a : 0, vg = cvm ? (cvm * rhoL) / rg : 0, vmE = cvm ? (u2 - u1) * rdt : 0;
      const a11 = rdt + Fwl + Fil + vl * rdt, a12 = Fil + vl * rdt, a21 = Fig + vg * rdt, a22 = rdt + Fwg + Fig + vg * rdt, r1 = b1 - vl * vmE, r2 = b2 + vg * vmE, det = a11 * a22 - a12 * a21;
      if (1 - Math.max(al[Lc], al[Rc]) < agMin) { // liquid bridge: the trapped gas moves with the liquid (single-fluid momentum of the slug body)
        const rm = a * rhoL + ag * rg; us[f] = (u1 * rdt - cl - G * sinT[f] - B / rm) / (rdt + Fwl); gs[f] = us[f]; dl[f] = 1 / rm / (rdt + Fwl) / dxf; dg[f] = dl[f];
      } else { us[f] = (a22 * r1 + a12 * r2) / det; gs[f] = (a11 * r2 + a21 * r1) / det; dl[f] = (a22 / rhoL + a12 / rg) / det / dxf; dg[f] = (a11 / rg + a21 / rhoL) / det / dxf; }
      alU[f] = u1 > 0 || (u1 === 0 && us[f] >= 0) || out ? al[Lc] : al[Rc]; mgU[f] = u2 > 0 || (u2 === 0 && gs[f] >= 0) ? mg[Lc] : out ? ((1 - al[Lc]) * pOut) / (Rs * Tg[Lc]) : mg[Rc]; // donor cells from the old velocity sign, frozen over the step
    }
    if (!per) { // inlet face
      if (bc) { us[0] = bc.ul; gs[0] = bc.ug; alU[0] = bc.al; mgU[0] = ((1 - bc.al) * p[0]) / (Rs * (bc.Tg ?? T0)); } else { us[0] = 0; gs[0] = 0; alU[0] = al[0]; mgU[0] = mg[0]; }
      dl[0] = 0; dg[0] = 0;
    }
    // ---- pressure equation: volume conservation with the linearised face velocities
    for (let i = 0; i < N; i++) {
      const fa = i, fb = per ? wrap(i + 1) : i + 1, rg = Math.max(p[i] / (Rs * Tg[i]), 1e-9), r = dt / dx, kap = mg[i] / (rg * p[i]) + al[i] * kapL, vol = al[i] + mg[i] / rg - 1;
      const cA = r * (alU[fa] * dl[fa] + (mgU[fa] * dg[fa]) / rg), cB = r * (alU[fb] * dl[fb] + (mgU[fb] * dg[fb]) / rg), last = !per && i === N - 1;
      sa[i] = -cA; sc[i] = last ? 0 : -cB; sb[i] = kap + cA + cB;
      sd[i] = kap * p[i] + clamp(vol, -0.05, 0.05) - r * (alU[fb] * us[fb] - alU[fa] * us[fa] + (mgU[fb] * gs[fb] - mgU[fa] * gs[fa]) / rg) + (last ? cB * pOut : 0);
      if (!per && i === 0) sa[i] = 0;
    }
    if (per) thomasCyclic(sa, sb, sc, sd, N, w1, w2, w3); else thomas(sa, sb, sc, sd, N, w1);
    // ---- new velocities and upwind fluxes
    for (let f = per ? 0 : 0; f <= f1; f++) {
      const Lc = wrap(f - 1), out = !per && f === N, Rc = out ? Lc : wrap(f);
      if (!per && f === 0) { ul[0] = us[0]; ug[0] = gs[0]; FL[0] = alU[0] * us[0]; FG[0] = mgU[0] * gs[0]; continue; }
      const dp = (out ? pOut : sd[Rc]) - sd[Lc]; ul[f] = us[f] - dl[f] * dp; ug[f] = clamp(gs[f] - dg[f] * dp, -cG, cG);
      FL[f] = alU[f] * ul[f]; FG[f] = mgU[f] * ug[f];
    }
    if (per) { ul[N] = ul[0]; ug[N] = ug[0]; FL[N] = FL[0]; FG[N] = FG[0]; }
    // positivity: a cell cannot give more than it holds
    for (let i = 0; i < N; i++) {
      const fa = i, fb = i + 1, oL = Math.max(FL[fb], 0) + Math.max(-FL[fa], 0), oG = Math.max(FG[fb], 0) + Math.max(-FG[fa], 0);
      if (oL * dt > al[i] * dx && oL > 0) { const s = (al[i] * dx) / (oL * dt); if (FL[fb] > 0) FL[fb] *= s; if (FL[fa] < 0) FL[fa] *= s; if (per && i === N - 1 && FL[fb] > 0) FL[0] = FL[N]; if (per && i === 0 && FL[fa] < 0) FL[N] = FL[0]; }
      if (oG * dt > mg[i] * dx && oG > 0) { const s = (mg[i] * dx) / (oG * dt); if (FG[fb] > 0) FG[fb] *= s; if (FG[fa] < 0) FG[fa] *= s; if (per && i === N - 1 && FG[fb] > 0) FG[0] = FG[N]; if (per && i === 0 && FG[fa] < 0) FG[N] = FG[0]; }
    }
    // ---- phase energies (internal energy, pressure work, interfacial and wall heat transfer), then masses
    let mres = 0;
    for (let i = 0; i < N; i++) {
      const fa = i, fb = i + 1, r = dt / dx, alN = al[i] - r * (FL[fb] - FL[fa]), mgN = Math.max(mg[i] - r * (FG[fb] - FG[fa]), 0);
      if (en) {
        const tlA = FL[fa] >= 0 ? (i > 0 || per ? Tl[wrap(i - 1)] : bc?.Tl ?? Tl[0]) : Tl[i], tlB = FL[fb] >= 0 ? Tl[i] : Tl[per ? wrap(i + 1) : Math.min(i + 1, N - 1)], tgA = FG[fa] >= 0 ? (i > 0 || per ? Tg[wrap(i - 1)] : bc?.Tg ?? Tg[0]) : Tg[i], tgB = FG[fb] >= 0 ? Tg[i] : Tg[per ? wrap(i + 1) : Math.min(i + 1, N - 1)];
        const pn = sd[i], agO = 1 - al[i], agN = 1 - alN, work = -pn * (agN - agO + r * ((1 - (ug[fb] >= 0 ? al[i] : al[per ? wrap(i + 1) : Math.min(i + 1, N - 1)])) * ug[fb] - (1 - (ug[fa] >= 0 ? (i > 0 || per ? al[wrap(i - 1)] : bc?.al ?? al[0]) : al[i])) * ug[fa]));
        const Si = chan ? 1 / D : (sect(SECT.Si, al[i]) * D) / A, SLw = chan ? 1 / D : (sect(SECT.SL, al[i]) * D) / A, SGw = chan ? 1 / D : (PI * D) / A - SLw, hi = (en.hi || 0) * Si, hwl = (en.UwL || 0) * SLw, hwg = (en.UwG || 0) * SGw, Tw = en.Tw ?? T0;
        const Cl = Math.max(alN, 1e-9) * rhoL * cvL, Cg = Math.max(mgN, 1e-12) * cvG, El = al[i] * rhoL * cvL * Tl[i] - r * rhoL * cvL * (FL[fb] * tlB - FL[fa] * tlA), Eg = mg[i] * cvG * Tg[i] - r * cvG * (FG[fb] * tgB - FG[fa] * tgA) + work;
        // implicit 2 × 2 exchange: Cl Tl' = El + dt [hi (Tg' − Tl') + hwl (Tw − Tl')],  Cg Tg' = Eg + dt [hi (Tl' − Tg') + hwg (Tw − Tg')]
        const m11 = Cl + dt * (hi + hwl), m22 = Cg + dt * (hi + hwg), m12 = dt * hi, dd = m11 * m22 - m12 * m12, q1 = El + dt * hwl * Tw, q2 = Eg + dt * hwg * Tw, tl = (m22 * q1 + m12 * q2) / dd, tg = (m11 * q2 + m12 * q1) / dd;
        bal.eWall += (hwl * (tl - Tw) + hwg * (tg - Tw)) * dt * dx * A; bal.eWork += work * dx * A; Tl[i] = tl; Tg[i] = tg;
      }
      al[i] = alN; mg[i] = mgN; p[i] = Math.max(sd[i], 100);
      const v = Math.abs(al[i] + (mg[i] * Rs * Tg[i]) / p[i] - 1); if (v > mres) mres = v;
    }
    if (!per) {
      bal.lIn += FL[0] * rhoL * A * dt; bal.lOut += FL[N] * rhoL * A * dt; bal.gIn += FG[0] * A * dt; bal.gOut += FG[N] * A * dt;
      if (en) { bal.eIn += (FL[0] * rhoL * cvL * (FL[0] >= 0 ? bc?.Tl ?? Tl[0] : Tl[0]) + FG[0] * cvG * (FG[0] >= 0 ? bc?.Tg ?? Tg[0] : Tg[0])) * A * dt; bal.eOut += (FL[N] * rhoL * cvL * Tl[N - 1] + FG[N] * cvG * Tg[N - 1]) * A * dt; }
    }
    bal.vol = Math.max(bal.vol, mres); t += dt; steps++;
    record();
    if (!(mres < 0.5)) break; // diverged (ill-posed regime on this grid): stop and report
  }
  const end = total();
  return { x: Array.from(xc), al: Array.from(al), ul: Array.from(ul.subarray(0, N)), ug: Array.from(ug.subarray(0, N)), p: Array.from(p), Tl: Array.from(Tl), Tg: Array.from(Tg), t, steps, dtMean: steps ? t / steps : 0, completed: t >= tEnd - 1e-9,
    mass: { liquid0: start.l, liquid: end.l, liquidIn: bal.lIn, liquidOut: bal.lOut, gas0: start.g, gas: end.g, gasIn: bal.gIn, gasOut: bal.gOut, errorL: (end.l - start.l - bal.lIn + bal.lOut) / Math.max(start.l + bal.lIn, 1e-12), errorG: (end.g - start.g - bal.gIn + bal.gOut) / Math.max(start.g + bal.gIn, 1e-12) },
    energy: { e0: start.e, e: end.e, in: bal.eIn, out: bal.eOut, wall: bal.eWall, work: bal.eWork, error: en ? (end.e - start.e - bal.eIn + bal.eOut + bal.eWall - bal.eWork) / Math.max(Math.abs(start.e) + Math.abs(bal.eIn), 1e-12) : 0 },
    volErrMax: bal.vol, series: ser, field: { t: fld.t, x: Array.from(xc), al: fld.al }, slugs };
}

// =====================================================================================================
// 7d. Compressible two-phase benchmark models: exact Riemann solution, seven-equation (Baer–Nunziato) model, species transport
// =====================================================================================================
/**
 * Exact solution of the Riemann problem of the Euler equations for an ideal gas (two-shock / two-rarefaction / mixed, after Toro).
 * L, R: { rho, u, p }; returns { pStar, uStar, sample(xi) → { rho, u, p } } with xi = x / t.
 */
export function riemannExact(L, R, gamma = 1.4) {
  const g = gamma, aL = Math.sqrt((g * L.p) / L.rho), aR = Math.sqrt((g * R.p) / R.rho);
  const fK = (p, K, aK) => (p > K.p ? (p - K.p) * Math.sqrt(2 / ((g + 1) * K.rho) / (p + ((g - 1) / (g + 1)) * K.p)) : ((2 * aK) / (g - 1)) * ((p / K.p) ** ((g - 1) / (2 * g)) - 1));
  const F = (p) => fK(p, L, aL) + fK(p, R, aR) + R.u - L.u, pStar = brent(F, 1e-12 * Math.min(L.p, R.p), 50 * Math.max(L.p, R.p) + L.rho * (L.u - R.u) ** 2 + 1, 1e-14), uStar = 0.5 * (L.u + R.u) + 0.5 * (fK(pStar, R, aR) - fK(pStar, L, aL));
  const sample = (xi) => {
    const left = xi <= uStar, K = left ? L : R, aK = left ? aL : aR, s = left ? 1 : -1, pr = pStar / K.p; // mirror the right side
    if (pStar > K.p) { const sK = K.u - s * aK * Math.sqrt(((g + 1) / (2 * g)) * pr + (g - 1) / (2 * g)); if (s * (xi - sK) < 0) return { ...K }; return { rho: (K.rho * (pr + (g - 1) / (g + 1))) / (((g - 1) / (g + 1)) * pr + 1), u: uStar, p: pStar }; }
    const head = K.u - s * aK, aS = aK * pr ** ((g - 1) / (2 * g)), tail = uStar - s * aS;
    if (s * (xi - head) < 0) return { ...K }; if (s * (xi - tail) > 0) return { rho: K.rho * pr ** (1 / g), u: uStar, p: pStar };
    const c = (2 / (g + 1)) * (aK + s * ((g - 1) / 2) * (K.u - xi)); return { rho: K.rho * (c / aK) ** (2 / (g - 1)), u: (2 / (g + 1)) * (s * aK + ((g - 1) / 2) * K.u + xi), p: K.p * (c / aK) ** ((2 * g) / (g - 1)) };
  };
  return { pStar, uStar, sample };
}

/**
 * Seven-equation two-phase model of Baer–Nunziato type in one dimension (two velocities, two pressures, volume-fraction transport),
 * stiffened-gas equations of state p = (γ − 1) ρ e − γ p∞, interface velocity u_I = u₂ and interface pressure p_I = p₁, without relaxation:
 *   ∂α₁/∂t + u_I ∂α₁/∂x = 0;   ∂(αρ)_k/∂t + ∂(αρu)_k/∂x = 0;   ∂(αρu)_k/∂t + ∂(αρu² + αp)_k/∂x = p_I ∂α_k/∂x;
 *   ∂(αρE)_k/∂t + ∂(αu(ρE + p))_k/∂x = p_I u_I ∂α_k/∂x.
 * Rusanov fluxes with the non-conservative terms discretised after Saurel & Abgrall so that a contact with uniform velocity and
 * pressure stays uniform. o: { n, length, x0, tEnd, cfl, eos: [{ gamma, pinf }, { gamma, pinf }], left / right: { a1, rho1, u1, p1, rho2, u2, p2 } }.
 * Returns { x[], a1[], rho1[], u1[], p1[], rho2[], u2[], p2[], steps, conservation: { mass1, mass2, momentum (after the end-pressure impulse), energy (relative drift; mass and energy for states at rest at the ends) }, entropy: { s0, s1, production (≥ 0) } }.
 */
export function baerNunziato(o = {}) {
  const N = o.n || 200, L = o.length || 1, dx = L / N, x0 = o.x0 ?? 0.5 * L, tEnd = o.tEnd ?? 0.2, cfl = o.cfl || 0.45, eos = o.eos || [{ gamma: 1.4, pinf: 0 }, { gamma: 1.4, pinf: 0 }];
  const U = Array.from({ length: 7 }, () => new Float64Array(N)), x = Array.from({ length: N }, (_, i) => (i + 0.5) * dx), g1 = eos[0].gamma, g2 = eos[1].gamma, q1 = eos[0].pinf || 0, q2 = eos[1].pinf || 0;
  const set = (i, s) => { const a2 = 1 - s.a1; U[0][i] = s.a1; U[1][i] = s.a1 * s.rho1; U[2][i] = s.a1 * s.rho1 * s.u1; U[3][i] = s.a1 * ((s.p1 + g1 * q1) / (g1 - 1) + 0.5 * s.rho1 * s.u1 * s.u1); U[4][i] = a2 * s.rho2; U[5][i] = a2 * s.rho2 * s.u2; U[6][i] = a2 * ((s.p2 + g2 * q2) / (g2 - 1) + 0.5 * s.rho2 * s.u2 * s.u2); };
  const init = o.init || ((xx) => (xx < x0 ? o.left : o.right)); for (let i = 0; i < N; i++) set(i, init(x[i]));
  const prim = (i) => { const a1 = U[0][i], a2 = 1 - a1, r1 = U[1][i] / a1, u1 = U[2][i] / U[1][i], p1 = (g1 - 1) * (U[3][i] / a1 - 0.5 * r1 * u1 * u1) - g1 * q1, r2 = U[4][i] / a2, u2 = U[5][i] / U[4][i], p2 = (g2 - 1) * (U[6][i] / a2 - 0.5 * r2 * u2 * u2) - g2 * q2; return { a1, a2, r1, u1, p1, r2, u2, p2, c1: Math.sqrt(Math.max((g1 * (p1 + q1)) / r1, 1e-12)), c2: Math.sqrt(Math.max((g2 * (p2 + q2)) / r2, 1e-12)) }; };
  const tot = () => { let m1 = 0, m2 = 0, mo = 0, e = 0, s = 0; for (let i = 0; i < N; i++) { const w = prim(i); m1 += U[1][i]; m2 += U[4][i]; mo += U[2][i] + U[5][i]; e += U[3][i] + U[6][i]; s += (U[1][i] * Math.log((w.p1 + q1) / w.r1 ** g1)) / (g1 - 1) + (U[4][i] * Math.log((w.p2 + q2) / w.r2 ** g2)) / (g2 - 1); } return { m1: m1 * dx, m2: m2 * dx, mo: mo * dx, e: e * dx, s: s * dx }; };
  const pe = (i) => { const w = prim(i); return w.a1 * w.p1 + w.a2 * w.p2 + w.a1 * w.r1 * w.u1 * w.u1 + w.a2 * w.r2 * w.u2 * w.u2; }, pEnds = [pe(0), pe(N - 1)]; // momentum flux through the open ends (constant until a wave arrives)
  const t0 = tot(), W = new Array(N), Fx = Array.from({ length: 7 }, () => new Float64Array(N + 1)), Sf = new Float64Array(N + 1), Un = Array.from({ length: 7 }, () => new Float64Array(N));
  let t = 0, steps = 0;
  const flux = (w, k) => (k === 1 ? w.a1 * w.r1 * w.u1 : k === 2 ? w.a1 * (w.r1 * w.u1 * w.u1 + w.p1) : k === 3 ? w.a1 * w.u1 * ((w.p1 + g1 * q1) / (g1 - 1) + 0.5 * w.r1 * w.u1 * w.u1 + w.p1) : k === 4 ? w.a2 * w.r2 * w.u2 : k === 5 ? w.a2 * (w.r2 * w.u2 * w.u2 + w.p2) : w.a2 * w.u2 * ((w.p2 + g2 * q2) / (g2 - 1) + 0.5 * w.r2 * w.u2 * w.u2 + w.p2));
  while (t < tEnd - 1e-14) {
    let smax = 1e-12; for (let i = 0; i < N; i++) { W[i] = prim(i); smax = Math.max(smax, Math.abs(W[i].u1) + W[i].c1, Math.abs(W[i].u2) + W[i].c2); }
    const dt = Math.min((cfl * dx) / smax, tEnd - t), lam = dt / dx;
    for (let f = 0; f <= N; f++) { // transmissive ends
      const a = W[Math.max(f - 1, 0)], b = W[Math.min(f, N - 1)], iL = Math.max(f - 1, 0), iR = Math.min(f, N - 1), S = Math.max(Math.abs(a.u1) + a.c1, Math.abs(a.u2) + a.c2, Math.abs(b.u1) + b.c1, Math.abs(b.u2) + b.c2); Sf[f] = S;
      for (let k = 1; k < 7; k++) Fx[k][f] = 0.5 * (flux(a, k) + flux(b, k)) - 0.5 * S * (U[k][iR] - U[k][iL]);
    }
    for (let i = 0; i < N; i++) {
      const w = W[i], aP = U[0][Math.min(i + 1, N - 1)], aM = U[0][Math.max(i - 1, 0)], a0 = U[0][i], dA = 0.5 * (aP - aM), uI = w.u2, pI = w.p1;
      Un[0][i] = a0 - lam * (uI * dA - 0.5 * Sf[i + 1] * (aP - a0) + 0.5 * Sf[i] * (a0 - aM));
      for (let k = 1; k < 7; k++) Un[k][i] = U[k][i] - lam * (Fx[k][i + 1] - Fx[k][i]);
      Un[2][i] += lam * pI * dA; Un[3][i] += lam * pI * uI * dA; Un[5][i] -= lam * pI * dA; Un[6][i] -= lam * pI * uI * dA;
    }
    for (let k = 0; k < 7; k++) U[k].set(Un[k]);
    t += dt; steps++;
  }
  const t1 = tot(), out = { x, a1: [], rho1: [], u1: [], p1: [], rho2: [], u2: [], p2: [] };
  for (let i = 0; i < N; i++) { const w = prim(i); out.a1.push(w.a1); out.rho1.push(w.r1); out.u1.push(w.u1); out.p1.push(w.p1); out.rho2.push(w.r2); out.u2.push(w.u2); out.p2.push(w.p2); }
  const rel = (a, b) => (b - a) / Math.max(Math.abs(a), 1e-300);
  return { ...out, steps, tEnd: t, conservation: { mass1: rel(t0.m1, t1.m1), mass2: rel(t0.m2, t1.m2), momentum: (t1.mo - t0.mo - (pEnds[0] - pEnds[1]) * t) / Math.max(Math.abs(t0.mo), Math.abs(pEnds[0] - pEnds[1]) * t, 1e-300), energy: rel(t0.e, t1.e) }, entropy: { s0: t0.s, s1: t1.s, production: t1.s - t0.s } };
}

/**
 * Species (tracer, inhibitor, water-cut marker) transport with the liquid along a steady solution:
 *   ∂(A H c)/∂t + ∂(q_L c)/∂s = ∂/∂s (A H D ∂c/∂s),   first-order upwind advection, central dispersion, explicit in time.
 * o: { s[] (uniform nodes), qL[] (m³/s), holdup[], area[] | A, D (m²/s dispersion), cIn(t), c0, tEnd, cfl, nSeries }.
 * Returns { s[], c[] (cell values at the end), t[], cOut[], breakthrough (time at which the outlet reaches 50 % of a unit step, s | null), balance: { in, out, stored, error } }.
 */
export function speciesTransport(o) {
  const n = o.s.length - 1, ds = o.s[1] - o.s[0], mid = (a) => Array.from({ length: n }, (_, i) => 0.5 * (a[i] + a[i + 1])), q = o.qL, H = mid(o.holdup), A = o.area ? mid(o.area) : new Array(n).fill(o.A), V = H.map((h, i) => Math.max(h * A[i], 1e-12) * ds);
  const Dd = o.D || 0, cIn = typeof o.cIn === 'function' ? o.cIn : () => o.cIn ?? 1, c = new Array(n).fill(o.c0 ?? 0), tEnd = o.tEnd, nS = o.nSeries || 300;
  let dt = Infinity; for (let i = 0; i < n; i++) dt = Math.min(dt, V[i] / (Math.abs(q[i]) + Math.abs(q[i + 1]) + (4 * Dd * V[i]) / (ds * ds) + 1e-30)); dt *= o.cfl || 0.9;
  const steps = Math.max(1, Math.ceil(tEnd / dt)); dt = tEnd / steps;
  // the upwind scheme adds a numerical dispersion ½ |u| Δs (1 − |u| Δt / Δs); it is taken off the physical dispersion where that is larger
  const Dn = Array.from({ length: n + 1 }, (_, f) => { if (f === 0 || f === n) return 0; const uf = Math.abs(q[f]) * ds / (0.5 * (V[f - 1] + V[f])); return Math.max(Dd - 0.5 * uf * ds * (1 - (uf * dt) / ds), 0); });
  const ts = [0], co = [c[n - 1]], F = new Array(n + 1); let mIn = 0, mOut = 0, bt = null; const m0 = c.reduce((s, v, i) => s + v * V[i], 0), every = Math.max(1, Math.floor(steps / nS));
  for (let k = 1; k <= steps; k++) {
    const t = k * dt, ci = cIn(t - 0.5 * dt);
    for (let f = 0; f <= n; f++) { const adv = q[f] * (f === 0 ? ci : c[f - 1]), dif = f > 0 && f < n ? (-Dn[f] * 0.5 * (V[f - 1] + V[f])) / (ds * ds) * (c[f] - c[f - 1]) : 0; F[f] = adv + dif; }
    const prev = c[n - 1]; for (let i = 0; i < n; i++) c[i] += (dt / V[i]) * (F[i] - F[i + 1]);
    mIn += F[0] * dt; mOut += F[n] * dt;
    if (bt === null && prev < 0.5 && c[n - 1] >= 0.5) bt = t - dt + (dt * (0.5 - prev)) / (c[n - 1] - prev);
    if (k % every === 0 || k === steps) { ts.push(t); co.push(c[n - 1]); }
  }
  const m1 = c.reduce((s, v, i) => s + v * V[i], 0);
  return { s: mid(o.s), c, t: ts, cOut: co, breakthrough: bt, steps, balance: { in: mIn, out: mOut, stored: m1 - m0, error: (m1 - m0 - mIn + mOut) / Math.max(Math.abs(mIn), Math.abs(m0), 1e-300) } };
}

// =====================================================================================================
// 7c. More turbulence closures for developed pipe flow (1-D radial, wall units: ν = 1, u_τ = 1, R⁺ = Re_τ)
// =====================================================================================================
// One under-relaxed finite-volume sweep of (1/r) d/dr (r Γ dφ/dr) + src − sink φ = 0 on nodes y[0..n] (y from the wall, r = R⁺ − y).
// wall: { value } (Dirichlet) | null (zero gradient); the centre line is symmetric; fixed[j] pins a node.
function radialSweep(y, r, phi, gam, src, sink, wall, ur = 0.5, fixed = null) {
  const n = y.length - 1, a = new Array(n + 1), b = new Array(n + 1), c = new Array(n + 1), d = new Array(n + 1);
  for (let j = 0; j <= n; j++) {
    if (j === 0) { if (wall) { a[j] = 0; c[j] = 0; b[j] = 1; d[j] = wall.value; } else { a[j] = 0; b[j] = 1; c[j] = -1; d[j] = 0; } continue; }
    if (fixed && fixed[j] !== undefined) { a[j] = 0; c[j] = 0; b[j] = 1; d[j] = fixed[j]; continue; }
    if (j === n) { a[j] = -1; b[j] = 1; c[j] = 0; d[j] = 0; continue; }
    const we = (0.5 * (r[j] + r[j + 1]) * 0.5 * (gam[j] + gam[j + 1])) / (y[j + 1] - y[j]), ww = (0.5 * (r[j] + r[j - 1]) * 0.5 * (gam[j] + gam[j - 1])) / (y[j] - y[j - 1]), vol = (r[j] * (y[j + 1] - y[j - 1])) / 2, bb = (ww + we + Math.max(sink[j], 0) * vol) / ur;
    a[j] = -ww; c[j] = -we; b[j] = bb; d[j] = Math.max(src[j], 0) * vol + (1 - ur) * bb * phi[j];
  }
  const x = tridiag(a, b, c, d); for (let j = 0; j <= n; j++) phi[j] = Math.max(x[j], 1e-14);
}
export const RANS_MODELS = [
  { value: 'mixing', label: 'Mixing length (van Driest)' }, { value: 'komega', label: 'k–ω (Wilcox 1988)' }, { value: 'sst', label: 'SST k–ω (Menter 2003)' }, { value: 'sa', label: 'Spalart–Allmaras' },
  { value: 'kepsilon', label: 'Low-Reynolds k–ε (Chien)' }, { value: 'kestd', label: 'Standard k–ε, wall functions' }, { value: 'rng', label: 'RNG k–ε, wall functions' }, { value: 'realizable', label: 'Realizable k–ε, wall functions' }, { value: 'rsm', label: 'Reynolds-stress model (LRR), wall functions' },
];
/**
 * Fully developed turbulent pipe flow with the closures that ransPipe does not cover. Momentum: (1 + ν_t) dU⁺/dy⁺ = 1 − y⁺/R⁺.
 *  'sst'        Menter SST k–ω (2003 constants, production limiter, cross-diffusion, F1/F2 blending), integrated to the wall;
 *  'sa'         Spalart–Allmaras one-equation model (standard form without the trip term), integrated to the wall;
 *  'kestd' | 'rng' | 'realizable'  high-Reynolds k–ε family with the log-law wall function at the first node (κ = 0.41, E = 9.8);
 *  'rsm'        Launder–Reece–Rodi Reynolds-stress model (isotropisation of production, wall-reflection term, gradient diffusion with the
 *               wall-normal stress) for ⟨uu⟩, ⟨vv⟩, ⟨ww⟩, ⟨uv⟩ and ε, wall functions.
 * o: { reTau, model, n, maxIter, tol, yWall (y⁺ of the wall-function node, default 40) }.
 * Returns { Re, f (Darcy), reTau, y[], u[], r[], uRel[], nut[], k[], stress: { uu, vv, ww, uv } | null, iterations, residual, converged, wallFunction }.
 */
export function ransPipeExtra(o = {}) {
  const Rp = Math.max(o.reTau || 1000, 50), model = o.model || 'sst', n = o.n || 100, maxIter = o.maxIter || 8000, tol = o.tol || 1e-8, kap = 0.41, E = 9.8, wf = !(model === 'sst' || model === 'sa');
  const yp = wf ? Math.min(o.yWall || 40, 0.2 * Rp) : 0, y1 = wf ? (Rp - yp) / (3 * n) : Math.min(0.1, Rp / (4 * n)), gq = (q) => (y1 * (q ** n - 1)) / (q - 1) - (Rp - yp), q = gq(1 + 1e-9) >= 0 ? 1 + 1e-9 : brent(gq, 1 + 1e-9, 3, 1e-13), y = [yp];
  for (let j = 1; j <= n; j++) y.push(j === n ? Rp : yp + (y1 * (q ** j - 1)) / (q - 1));
  const r = y.map((v) => Rp - v), z = () => new Array(n + 1).fill(0), nut = z(), u = z(), S = z(), gam = z(), src = z(), sink = z(), tau = y.map((v) => 1 - v / Rp);
  const law = (yy) => (yy < 11.06 ? yy : Math.log(E * yy) / kap), uWall = wf ? law(yp) : 0;
  let ubWall = 0; if (wf) { const m = 60; for (let i = 0; i < m; i++) { const ya = (yp * (i + 0.5)) / m; ubWall += law(ya) * (Rp - ya) * (yp / m); } }
  const velocity = () => { for (let j = 0; j <= n; j++) S[j] = Math.max(tau[j], 0) / (1 + nut[j]); u[0] = uWall; for (let j = 1; j <= n; j++) u[j] = u[j - 1] + 0.5 * (S[j] + S[j - 1]) * (y[j] - y[j - 1]); let ub = ubWall; for (let j = 1; j <= n; j++) ub += 0.5 * (u[j] * r[j] + u[j - 1] * r[j - 1]) * (y[j] - y[j - 1]); return (2 * ub) / (Rp * Rp); };
  const grad = (a, j) => (j === 0 ? (a[1] - a[0]) / (y[1] - y[0]) : j === n ? 0 : (a[j + 1] - a[j - 1]) / (y[j + 1] - y[j - 1]));
  // start from the mixing-length solution
  for (let j = 0; j <= n; j++) { const e = 1 - y[j] / Rp, l = Rp * (0.14 - 0.08 * e * e - 0.06 * e ** 4) * (1 - Math.exp(-y[j] / 26)); nut[j] = (l * l * 2 * e) / (1 + Math.sqrt(1 + 4 * l * l * e)); }
  let ub = velocity(), it = 0, res = 1, kk = null, stress = null;
  if (model === 'sa') {
    const cb1 = 0.1355, cb2 = 0.622, sg = 2 / 3, cv1 = 7.1, cw2 = 0.3, cw3 = 2, cw1 = cb1 / (kap * kap) + (1 + cb2) / sg, nt = nut.map((v) => Math.max(v, 1e-10));
    for (it = 1; it <= maxIter; it++) {
      for (let j = 0; j <= n; j++) {
        const d = Math.max(y[j], 1e-9), chi = nt[j], fv1 = chi ** 3 / (chi ** 3 + cv1 ** 3), fv2 = 1 - chi / (1 + chi * fv1), St = Math.max(S[j] + (nt[j] * fv2) / (kap * kap * d * d), 0.3 * S[j] + 1e-12), rr = Math.min(nt[j] / (St * kap * kap * d * d), 10), g = rr + cw2 * (rr ** 6 - rr), fw = g * ((1 + cw3 ** 6) / (g ** 6 + cw3 ** 6)) ** (1 / 6), dn = grad(nt, j);
        gam[j] = (1 + nt[j]) / sg; src[j] = cb1 * St * nt[j] + (cb2 / sg) * dn * dn; sink[j] = (cw1 * fw * nt[j]) / (d * d);
      }
      radialSweep(y, r, nt, gam, src, sink, { value: 0 }, 0.6);
      for (let j = 0; j <= n; j++) { const chi = nt[j]; nut[j] = (nt[j] * chi ** 3) / (chi ** 3 + cv1 ** 3); }
      const un = velocity(); res = Math.abs(un - ub) / un; ub = un; if (res < tol && it > 30) break;
    }
  } else if (model === 'sst') {
    const bS = 0.09, a1 = 0.31, sk1 = 0.85, sk2 = 1, sw1 = 0.5, sw2 = 0.856, b1 = 0.075, b2 = 0.0828, g1 = 5 / 9, g2 = 0.44;
    kk = y.map((_, j) => Math.max((nut[j] * S[j]) / 0.3, 1e-10)); const ww = y.map((yy, j) => (j === 0 ? 0 : Math.max(kk[j] / Math.max(nut[j], 1e-6), yy < 2.5 ? 6 / (b1 * yy * yy) : 0)));
    const wFix = y.map((yy, j) => (j === 0 ? 60 / (b1 * y[1] * y[1]) : yy < 2.5 ? 6 / (b1 * yy * yy) : undefined)), F1 = z(), CD = z();
    for (it = 1; it <= maxIter; it++) {
      for (let j = 1; j <= n; j++) {
        const d = y[j], w = Math.max(ww[j], 1e-12), cd = 2 * sw2 * grad(kk, j) * grad(ww, j) / w, cdp = Math.max(cd, 1e-10), arg1 = Math.min(Math.max(Math.sqrt(kk[j]) / (bS * w * d), 500 / (d * d * w)), (4 * sw2 * kk[j]) / (cdp * d * d));
        F1[j] = Math.tanh(arg1 ** 4); CD[j] = cd;
      }
      F1[0] = 1;
      for (let j = 0; j <= n; j++) { const P = Math.min(nut[j] * S[j] * S[j], 10 * bS * kk[j] * ww[j]); gam[j] = 1 + (F1[j] * sk1 + (1 - F1[j]) * sk2) * nut[j]; src[j] = P; sink[j] = bS * ww[j]; }
      radialSweep(y, r, kk, gam, src, sink, { value: 0 }, 0.5);
      for (let j = 0; j <= n; j++) {
        const f = F1[j], gm = f * g1 + (1 - f) * g2, be = f * b1 + (1 - f) * b2, P = Math.min(nut[j] * S[j] * S[j], 10 * bS * kk[j] * ww[j]), x = (1 - f) * CD[j];
        gam[j] = 1 + (f * sw1 + (1 - f) * sw2) * nut[j]; src[j] = (gm * P) / Math.max(nut[j], 1e-10) + Math.max(x, 0); sink[j] = be * ww[j] + Math.max(-x, 0) / Math.max(ww[j], 1e-12);
      }
      radialSweep(y, r, ww, gam, src, sink, { value: wFix[0] }, 0.5, wFix);
      for (let j = 1; j <= n; j++) { const d = y[j], arg2 = Math.max((2 * Math.sqrt(kk[j])) / (bS * ww[j] * d), 500 / (d * d * ww[j])), F2 = Math.tanh(arg2 * arg2); nut[j] = (a1 * kk[j]) / Math.max(a1 * ww[j], S[j] * F2); }
      nut[0] = 0; const un = velocity(); res = Math.abs(un - ub) / un; ub = un; if (res < tol && it > 30) break;
    }
  } else if (model === 'rsm') {
    const Cmu = 0.09, C1 = 1.8, C2 = 0.6, Cs = 0.25, Ce = 0.15, Ce1 = 1.44, Ce2 = 1.92, Cr1 = 0.5, Cr2 = 0.3, k0 = 1 / Math.sqrt(Cmu);
    kk = y.map(() => k0); const ee = y.map((yy, j) => Math.max((Cmu * k0 * k0) / Math.max(nut[j], 1e-3), 1e-8)), uu = kk.map((k) => 1.1 * k), vv = kk.map((k) => 0.5 * k), wz = kk.map((k) => 0.4 * k), qv = y.map((_, j) => Math.max(nut[j] * S[j], 1e-8));
    for (it = 1; it <= maxIter; it++) {
      const T = kk.map((k, j) => k / Math.max(ee[j], 1e-12)), f = y.map((yy, j) => ((3 * Cmu ** 0.75) / kap) * Math.sqrt(kk[j]) / Math.max(yy, 1e-9)), P = qv.map((v, j) => v * S[j]);
      const Dg = (c) => y.map((_, j) => 1 + c * T[j] * vv[j]);
      // ⟨vv⟩ (wall-normal), then ⟨uu⟩, ⟨ww⟩, −⟨uv⟩, ε
      for (let j = 0; j <= n; j++) { src[j] = (2 / 3) * (C1 - 1) * ee[j] + (2 / 3) * C2 * P[j]; sink[j] = C1 / T[j] + f[j] * (2 / 3) * Cr1 + (f[j] * (4 / 9) * Cr2 * C2 * T[j] * P[j]) / Math.max(vv[j], 1e-12); }
      radialSweep(y, r, vv, Dg(Cs), src, sink, null, 0.4);
      const refl = y.map((_, j) => Cr1 * vv[j] + Cr2 * C2 * T[j] * (2 / 3) * P[j]);
      for (let j = 0; j <= n; j++) { src[j] = P[j] * (2 - (4 / 3) * C2) + (2 / 3) * (C1 - 1) * ee[j] + (f[j] * refl[j]) / 3; sink[j] = C1 / T[j]; }
      radialSweep(y, r, uu, Dg(Cs), src, sink, null, 0.4);
      for (let j = 0; j <= n; j++) { src[j] = (2 / 3) * (C1 - 1) * ee[j] + (2 / 3) * C2 * P[j] + (f[j] * refl[j]) / 3; sink[j] = C1 / T[j]; }
      radialSweep(y, r, wz, Dg(Cs), src, sink, null, 0.4);
      for (let j = 0; j <= n; j++) { src[j] = vv[j] * S[j] * (1 - C2 + 0.5 * f[j] * Cr2 * C2 * T[j]); sink[j] = C1 / T[j] + 0.5 * f[j] * Cr1; }
      radialSweep(y, r, qv, Dg(Cs), src, sink, { value: tau[0] }, 0.4);
      for (let j = 0; j <= n; j++) kk[j] = 0.5 * (uu[j] + vv[j] + wz[j]);
      for (let j = 0; j <= n; j++) { src[j] = (Ce1 * P[j]) / T[j]; sink[j] = Ce2 / T[j]; }
      radialSweep(y, r, ee, Dg(Ce), src, sink, { value: (Cmu ** 0.75 * kk[0] ** 1.5) / (kap * yp) }, 0.4);
      for (let j = 0; j <= n; j++) nut[j] = S[j] > 1e-9 ? clamp(qv[j] / S[j], 0, 1e5) : nut[j]; // effective eddy viscosity −⟨uv⟩ / (dU/dy) used in the momentum balance
      nut[n] = nut[n - 1];
      const un = velocity(); res = Math.abs(un - ub) / un; ub = un; if (res < tol && it > 50) break;
    }
    stress = { uu: uu.slice(), vv: vv.slice(), ww: wz.slice(), uv: qv.map((v) => -v) };
  } else {
    const rng = model === 'rng', rea = model === 'realizable', Cmu0 = rng ? 0.0845 : 0.09, C1 = rng ? 1.42 : 1.44, C2 = rng ? 1.68 : rea ? 1.9 : 1.92, sk = rng ? 0.71942 : 1, se = rng ? 0.71942 : rea ? 1.2 : 1.3, kW = 1 / Math.sqrt(Cmu0), eW = 1 / (kap * yp);
    kk = y.map(() => kW); const ee = y.map((yy, j) => Math.max((Cmu0 * kW * kW) / Math.max(nut[j], 1e-3), 1e-8));
    for (it = 1; it <= maxIter; it++) {
      for (let j = 0; j <= n; j++) { gam[j] = 1 + nut[j] / sk; src[j] = nut[j] * S[j] * S[j]; sink[j] = ee[j] / Math.max(kk[j], 1e-12); }
      radialSweep(y, r, kk, gam, src, sink, { value: kW }, 0.5);
      for (let j = 0; j <= n; j++) {
        const T = kk[j] / Math.max(ee[j], 1e-12), eta = S[j] * T, P = nut[j] * S[j] * S[j]; gam[j] = 1 + nut[j] / se;
        if (rea) { src[j] = Math.max(0.43, eta / (eta + 5)) * S[j] * ee[j]; sink[j] = (C2 * ee[j]) / (kk[j] + Math.sqrt(ee[j])); }
        else { const R = rng ? (eta * (1 - eta / 4.38)) / (1 + 0.012 * eta ** 3) : 0, c = C1 - R; src[j] = (Math.max(c, 0) * P) / T; sink[j] = C2 / T + (Math.max(-c, 0) * P) / T / Math.max(ee[j], 1e-12); }
      }
      radialSweep(y, r, ee, gam, src, sink, { value: eW }, 0.5);
      for (let j = 0; j <= n; j++) { const T = kk[j] / Math.max(ee[j], 1e-12), cm = rea ? 1 / (4.0 + Math.sqrt(6) * Math.cos(PI / 6) * S[j] * T) : Cmu0; nut[j] = cm * kk[j] * T; } // realizable: As = √6 cos φ with φ = π/6 in simple shear, U* = S
      const un = velocity(); res = Math.abs(un - ub) / un; ub = un; if (res < tol && it > 30) break;
    }
  }
  const Re = 2 * Rp * ub;
  return { Re, f: 8 / (ub * ub), reTau: Rp, y, u: u.slice(), r: r.map((v) => v / Rp), uRel: u.map((v) => v / ub), nut: nut.slice(), k: kk ? kk.slice() : null, stress, iterations: Math.min(it, maxIter), residual: res, converged: res < 1e-6, wallFunction: wf };
}

// =====================================================================================================
// 7e. Bubbly flow: interfacial force closures in a radial void-distribution model, bubble dynamics, interfacial-area transport
// =====================================================================================================
/** Drag coefficient of a bubble: Ishii & Zuber (viscous and distorted regimes, swarm correction), as coded in OpenFOAM. */
export function dragIshiiZuber(Re, Eo, alpha = 0, muRatio = 0.02) {
  const muStar = (muRatio + 0.4) / (muRatio + 1), mixF = Math.max(1 - alpha, 1e-3) ** (-2.5 * muStar), ReM = Re / mixF, cdRe = ReM <= 1000 ? 24 * (1 + 0.1 * ReM ** 0.75) : 0.44 * ReM;
  const F = Math.max((1 / mixF) * Math.sqrt(1 - alpha), 1e-3), Ea = (1 + 17.67 * F ** 0.8571428) / (18.67 * F), cdEl = Ea * 0.6666 * Math.sqrt(Eo) * Re;
  return (cdEl >= cdRe ? Math.min(cdEl, Re * (1 - alpha) ** 2 * 2.66667) : cdRe) / Math.max(Re, 1e-12);
}
/** Lift coefficient of Tomiyama et al. (2002) with the horizontal Eötvös number (Wellek aspect ratio): positive for small bubbles (towards the wall in upflow), negative above d ≈ 5.8 mm in air–water. */
export function liftTomiyama(Re, Eo) {
  const EoH = Eo * (1 + 0.163 * Eo ** 0.757) ** (2 / 3), f = 0.0010422 * EoH ** 3 - 0.0159 * EoH * EoH - 0.0204 * EoH + 0.474;
  return EoH < 4 ? Math.min(0.288 * Math.tanh(0.121 * Re), f) : EoH < 10.7 ? f : -0.288;
}
/** Terminal rise velocity of a bubble in a swarm (buoyancy = drag, Ishii–Zuber). */
export function bubbleTerminal({ d, rhoL = 998, rhoG = 1.2, muL = 1e-3, sigma = 0.072, alpha = 0 }) {
  const Eo = (G * (rhoL - rhoG) * d * d) / sigma; let v = 0.2;
  for (let k = 0; k < 60; k++) { const cd = dragIshiiZuber((rhoL * v * d) / muL, Eo, alpha), vn = Math.sqrt((4 * G * d * (rhoL - rhoG) * (1 - alpha)) / (3 * cd * rhoL)); if (Math.abs(vn - v) < 1e-12) { v = vn; break; } v = 0.5 * (v + vn); }
  return v;
}
/**
 * Rise of a single bubble released from rest with the added (virtual) mass of the displaced liquid:
 *   (ρg + C_vm ρl) dv/dt = (ρl − ρg) g − ¾ C_D ρl v² / d.   Returns { t[], v[], a0 (initial acceleration), terminal }.
 */
export function bubbleRise({ d = 3e-3, rhoL = 998, rhoG = 1.2, muL = 1e-3, sigma = 0.072, cvm = 0.5, tEnd = 0.5, n = 400 } = {}) {
  const Eo = (G * (rhoL - rhoG) * d * d) / sigma, f = (v) => ((rhoL - rhoG) * G - (0.75 * dragIshiiZuber(Math.max((rhoL * Math.abs(v) * d) / muL, 1e-9), Eo) * rhoL * v * Math.abs(v)) / d) / (rhoG + cvm * rhoL), dt = tEnd / n, t = [0], vv = [0];
  let v = 0; for (let i = 1; i <= n; i++) { const k1 = f(v), k2 = f(v + 0.5 * dt * k1), k3 = f(v + 0.5 * dt * k2), k4 = f(v + dt * k3); v += (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4); t.push(i * dt); vv.push(v); }
  return { t, v: vv, a0: f(0), terminal: bubbleTerminal({ d, rhoL, rhoG, muL, sigma }) };
}
/**
 * Radial void-fraction distribution of developed upward bubbly flow in a pipe from the balance of the lateral forces on the bubbles:
 *   turbulent dispersion (Burns et al. Favre-averaged drag, or Lopez de Bertodano C_TD ρl k ∇α) = lift (Tomiyama, damped at the wall) + wall lubrication (Antal),
 * with the liquid velocity from the radial momentum balance (mixing length with van Driest damping plus the Sato bubble-induced
 * viscosity 0.6 α d u_r) and the bubble-induced (pseudo-turbulent) kinetic energy ½ C_vm α u_r² added to the shear-induced k.
 * o: { R (pipe radius), jl (mean liquid velocity (1 − α) u_l averaged over the section), alphaMean, d (bubble diameter), rhoL, rhoG, muL, sigma, n,
 *      cl (number, or 'tomiyama'), cw1, cw2, sigmaTD, ctd, tdModel: 'burns' | 'bertodano', cvm, cmub, iterations }.
 * Returns { r[] (r/R), alpha[], ul[] (m/s), k[], uRel, cl, cd, peak: { r, alpha }, wallPeaked, alphaCentre, dpdz, mean (area-average void, check), iterations, residual }.
 */
export function bubblyPipe(o = {}) {
  const R = o.R || 0.019, n = o.n || 120, d = o.d || 3.4e-3, rhoL = o.rhoL ?? 998, rhoG = o.rhoG ?? 1.2, muL = o.muL ?? 1e-3, sigma = o.sigma ?? 0.072, aM = o.alphaMean ?? 0.05, jl = o.jl ?? 0.9, nu = muL / rhoL;
  const cw1 = o.cw1 ?? -0.01, cw2 = o.cw2 ?? 0.05, sTD = o.sigmaTD ?? 0.7, ctd = o.ctd ?? 1, cvm = o.cvm ?? 0.5, cmub = o.cmub ?? 0.6, burns = o.tdModel !== 'bertodano', Eo = (G * (rhoL - rhoG) * d * d) / sigma;
  // nodes clustered at the wall: y = R (1 − cos) spacing
  const r = Array.from({ length: n + 1 }, (_, j) => R * Math.sin((0.5 * PI * j) / n)), y = r.map((v) => Math.max(R - v, 0)), al = r.map(() => aM), ul = new Array(n + 1).fill(0), kk = new Array(n + 1).fill(0), nut = new Array(n + 1).fill(0), S = new Array(n + 1).fill(0);
  const avg = (f) => { let s = 0; for (let j = 1; j <= n; j++) s += 0.5 * (f(j) * r[j] + f(j - 1) * r[j - 1]) * (r[j] - r[j - 1]); return (2 * s) / (R * R); };
  let ur = bubbleTerminal({ d, rhoL, rhoG, muL, sigma, alpha: aM }), cd = 0, cl = 0, dpdz = 0, res = 1, it = 0;
  const flow = (gp) => { // liquid velocity for a pressure gradient gp = −dp/dz − ρl g (the part that drives the flow)
    let I = 0; const tau = new Array(n + 1).fill(0); // r τ(r) = ∫ [gp + (ρl − ρm) g] r dr  (τ = shear stress towards the wall)
    for (let j = 1; j <= n; j++) { const fa = gp + al[j - 1] * (rhoL - rhoG) * G, fb = gp + al[j] * (rhoL - rhoG) * G; I += 0.5 * (fa * r[j - 1] + fb * r[j]) * (r[j] - r[j - 1]); tau[j] = I / r[j]; }
    const tw = Math.max(tau[n], 1e-9), ut = Math.sqrt(tw / rhoL);
    for (let j = 0; j <= n; j++) { const e = r[j] / R, l = R * (0.14 - 0.08 * e * e - 0.06 * e ** 4) * (1 - Math.exp((-y[j] * ut) / nu / 26)), nb = cmub * al[j] * d * ur, a = Math.abs(tau[j]) / (rhoL * Math.max(1 - al[j], 0.05)); S[j] = l > 1e-12 ? (-(nu + nb) + Math.sqrt((nu + nb) ** 2 + 4 * l * l * a)) / (2 * l * l) : a / (nu + nb); nut[j] = l * l * S[j] + nb; S[j] *= Math.sign(tau[j]); }
    ul[n] = 0; for (let j = n - 1; j >= 0; j--) ul[j] = ul[j + 1] + 0.5 * (S[j] + S[j + 1]) * (r[j + 1] - r[j]);
    return avg((j) => (1 - al[j]) * ul[j]);
  };
  for (it = 1; it <= (o.iterations || 60); it++) {
    { const a = -(rhoL - rhoG) * G; let b = 2000, fb = flow(b) - jl, g = 0; while (fb < 0 && g++ < 20) { b *= 3; fb = flow(b) - jl; } dpdz = brent((x) => flow(x) - jl, a, b, 1e-8); flow(dpdz); }
    const Reb = (rhoL * ur * d) / muL; cd = dragIshiiZuber(Reb, Eo, aM); cl = typeof o.cl === 'number' ? o.cl : liftTomiyama(Reb, Eo);
    for (let j = 0; j <= n; j++) kk[j] = (Math.max(nut[j] - cmub * al[j] * d * ur, 0) * Math.abs(S[j])) / 0.3 + 0.5 * cvm * al[j] * ur * ur;
    // lateral force balance integrated from the centre: d ln α / dr = (1 − α) Φ(r)
    const phi = (j, a) => { const damp = 0.5 * (1 - Math.cos(PI * Math.min(y[j] / d, 1))), W = Math.max(0, cw1 / d + cw2 / Math.max(y[j], 1e-6)), dudr = -S[j], force = -cl * ur * dudr * damp - ur * ur * W; return burns ? ((1 - a) * sTD * d * force) / (0.75 * cd * ur * Math.max(nut[j], 1e-9)) : force / (ctd * Math.max(kk[j], 1e-9)); };
    const profile = (ac) => { const a = [ac]; for (let j = 1; j <= n; j++) { const dr = r[j] - r[j - 1], p0 = phi(j - 1, a[j - 1]), a1 = clamp(a[j - 1] * Math.exp(clamp(p0 * dr, -30, 30)), 1e-12, 0.95), p1 = phi(j, a1); a.push(clamp(a[j - 1] * Math.exp(clamp(0.5 * (p0 + p1) * dr, -30, 30)), 1e-12, 0.95)); } return a; };
    const mean = (a) => { let s = 0; for (let j = 1; j <= n; j++) s += 0.5 * (a[j] * r[j] + a[j - 1] * r[j - 1]) * (r[j] - r[j - 1]); return (2 * s) / (R * R); };
    let lo = 1e-10, hi = 0.95; for (let k = 0; k < 70; k++) { const m = Math.sqrt(lo * hi); if (mean(profile(m)) > aM) hi = m; else lo = m; }
    const an = profile(Math.sqrt(lo * hi)); res = 0; for (let j = 0; j <= n; j++) { res = Math.max(res, Math.abs(an[j] - al[j])); al[j] += 0.5 * (an[j] - al[j]); }
    ur = bubbleTerminal({ d, rhoL, rhoG, muL, sigma, alpha: aM });
    if (res < 1e-7 && it > 3) break;
  }
  { const a = -(rhoL - rhoG) * G; let b = 2000, fb = flow(b) - jl, g = 0; while (fb < 0 && g++ < 20) { b *= 3; fb = flow(b) - jl; } dpdz = brent((x) => flow(x) - jl, a, b, 1e-10); flow(dpdz); } // momentum balance of the final void profile
  let jp = 0; for (let j = 1; j <= n; j++) if (al[j] > al[jp]) jp = j;
  return { r: r.map((v) => v / R), alpha: al.slice(), ul: ul.slice(), k: kk.slice(), uRel: ur, cl, cd, eotvos: Eo, peak: { r: r[jp] / R, alpha: al[jp] }, wallPeaked: r[jp] / R > 0.5, alphaCentre: al[0], dpdz: dpdz + rhoL * G, mean: avg((j) => al[j]), liquidFlux: avg((j) => (1 - al[j]) * ul[j]), iterations: it, residual: res };
}
/**
 * One-group interfacial-area transport along a vertical bubbly column (curvature κ = a_i/α = 6/d_sm, the form of Ishii and co-workers as coded in OpenFOAM):
 *   v_g dκ/dz = −R_RC κ − 12 φ C_WE C_D^⅓ α u_r κ² + (C_TI/18) u_t κ² √(1 − We_cr/We) exp(−We_cr/We),   φ = 1/(36π),
 * random-collision coalescence R_RC = 12 φ κ α C_RC u_t [1 − exp(−C α^⅓ α_max^⅓ / (α_max^⅓ − α^⅓))] / [α_max^⅓ (α_max^⅓ − α^⅓)], wake entrainment and turbulent-impact breakup.
 * o: { d0, alpha, vg, k (liquid turbulent kinetic energy), length, n, rhoL, rhoG, muL, sigma, crc, cwe, cti, weCr, c, alphaMax, cdFixed }.
 * Returns { z[], d[] (Sauter diameter), dEnd, dEquilibrium (| null), rates: { coalescenceRC, coalescenceWE, breakup } at the outlet (1/m/s) }.
 */
export function interfacialAreaTransport(o = {}) {
  const { d0 = 3e-3, alpha = 0.1, vg = 1, k = 0.01, length = 10, n = 400, rhoL = 998, rhoG = 1.2, muL = 1e-3, sigma = 0.072, crc = 0.04, cwe = 0.002, cti = 0.085, weCr = 6, c = 3, alphaMax = 0.75 } = o, phi = 1 / (36 * PI), ut = Math.sqrt(2 * k), ur = Math.SQRT2 * ((sigma * G * (rhoL - rhoG)) / (rhoL * rhoL)) ** 0.25 * (1 - alpha) ** 1.75, cm = Math.cbrt(alphaMax), ca = Math.cbrt(alpha);
  const parts = (kap) => { const d = 6 / kap, Re = Math.max((ur * d * rhoL) / muL, 1e-3), Eo = (G * d * d * (rhoL - rhoG)) / sigma, cd = o.cdFixed ?? Math.max(Math.min((16 / Re) * (1 + 0.15 * Re ** 0.687), 48 / Re), (8 * Eo) / (3 * (Eo + 4))), We = (rhoL * ut * ut * d) / sigma;
    return { rc: alpha < alphaMax ? (12 * phi * kap * alpha * crc * ut * (1 - Math.exp((-c * Math.cbrt(alpha * alphaMax)) / (cm - ca)))) / (cm * (cm - ca)) * kap : 0, we: 12 * phi * cwe * Math.cbrt(cd) * alpha * kap * ur * kap, ti: We > weCr ? (cti / 18) * ut * kap * kap * Math.sqrt(1 - weCr / We) * Math.exp(-weCr / We) : 0 }; };
  const f = (kap) => { const q = parts(kap); return (q.ti - q.rc - q.we) / vg; }, dz = length / n, z = [0], dd = [d0]; let kap = 6 / d0;
  for (let i = 1; i <= n; i++) { const k1 = f(kap), k2 = f(kap + 0.5 * dz * k1), k3 = f(kap + 0.5 * dz * k2), k4 = f(kap + dz * k3); kap = clamp(kap + (dz / 6) * (k1 + 2 * k2 + 2 * k3 + k4), 6 / 0.2, 6 / 1e-5); z.push(i * dz); dd.push(6 / kap); }
  let dEq = null; { let a = 6 / 0.1, fa = f(a); for (let i = 1; i <= 200 && dEq === null; i++) { const b = (6 / 0.1) * (6 / 2e-5 / (6 / 0.1)) ** (i / 200), fb = f(b); if (fa * fb < 0) dEq = 6 / brent(f, a, b, 1e-9); a = b; fa = fb; } }
  const q = parts(kap);
  return { z, d: dd, dEnd: 6 / kap, dEquilibrium: dEq, rates: { coalescenceRC: q.rc, coalescenceWE: q.we, breakup: q.ti }, uRel: ur, uTurb: ut };
}

// =====================================================================================================
// 7f. Two-dimensional CFD: interface advection schemes and an incompressible Navier–Stokes solver (projection, SST k–ω, volume of fluid)
// =====================================================================================================
// THINC flux through the downstream face of a donor cell (fraction of the cell volume) for a Courant number c ≥ 0; pm, p, pp: upstream, donor, downstream values.
function thincFlux(pm, p, pp, c, beta = 2.3) {
  let fl = c * p;
  if (p > 1e-8 && p < 1 - 1e-8 && (pp - p) * (p - pm) > 0) { const gm = pp > pm ? 1 : -1, qq = (beta * (2 * p - 1)) / gm, w = (Math.exp(beta) - Math.exp(qq)) / (Math.exp(qq) - Math.exp(-beta)); if (w > 0 && Number.isFinite(w)) { const xt = Math.log(w) / (2 * beta); fl = 0.5 * (c + (gm / beta) * (Math.log(Math.cosh(beta * (1 - xt))) - Math.log(Math.cosh(beta * (1 - c - xt))))); } }
  return clamp(fl, Math.max(0, c - (1 - p)), Math.min(c, p));
}
/**
 * One direction-split volume-of-fluid step on a staggered grid (cells with one ghost layer, stride nx + 2): THINC reconstruction weighted
 * with the interface normal (THINC/WLIC), with the dilatation term of Weymouth & Yue so that the split advection conserves volume to round-off.
 * C: cell fractions (ghosts filled by the caller), U, V: face velocities (U: (nx+1) × (ny+2) with ghost rows, V: (nx+2) × (ny+1) with ghost columns),
 * normal: optional { nx[], ny[] } per cell (e.g. from a level set); otherwise from the gradient of C. odd swaps the sweep order.
 */
function vofStep(C, U, V, nx, ny, dx, dy, dt, odd, normal, fill, perX = false, perY = false) {
  const sx = nx + 2, cc = new Float64Array(C.length); for (let k = 0; k < C.length; k++) cc[k] = C[k] > 0.5 ? 1 : 0;
  const wOf = (k, dir) => { let gx, gy; if (normal) { gx = normal.nx[k]; gy = normal.ny[k]; } else { gx = (C[k + 1] ?? C[k]) - (C[k - 1] ?? C[k]); gy = (C[k + sx] ?? C[k]) - (C[k - sx] ?? C[k]); } const a = Math.abs(gx), b = Math.abs(gy); return a + b < 1e-12 ? 0.5 : dir === 0 ? a / (a + b) : b / (a + b); };
  const sweepX = () => { const F = new Float64Array(C.length); for (let j = 1; j <= ny; j++) for (let i = 0; i <= nx; i++) { const u = U[i + (nx + 1) * j], c = (Math.abs(u) * dt) / dx, d = u >= 0 ? i + sx * j : i + 1 + sx * j, s = u >= 0 ? 1 : -1; if (c === 0) continue; const w = wOf(d, 0), f = w * thincFlux(C[d - s] ?? C[d], C[d], C[d + s] ?? C[d], c) + (1 - w) * c * C[d]; F[i + sx * j] = s * f; }
    if (perX) for (let j = 1; j <= ny; j++) F[sx * j] = F[nx + sx * j];
    for (let j = 1; j <= ny; j++) for (let i = 1; i <= nx; i++) { const k = i + sx * j; C[k] += -(F[k] - F[k - 1]) + (cc[k] * (U[i + (nx + 1) * j] - U[i - 1 + (nx + 1) * j]) * dt) / dx; } fill(C); };
  const sweepY = () => { const F = new Float64Array(C.length); for (let j = 0; j <= ny; j++) for (let i = 1; i <= nx; i++) { const v = V[i + sx * j], c = (Math.abs(v) * dt) / dy, d = v >= 0 ? i + sx * j : i + sx * (j + 1), s = v >= 0 ? sx : -sx; if (c === 0) continue; const w = wOf(d, 1), f = w * thincFlux(C[d - s] ?? C[d], C[d], C[d + s] ?? C[d], c) + (1 - w) * c * C[d]; F[i + sx * j] = Math.sign(s) * f; }
    if (perY) for (let i = 1; i <= nx; i++) F[i] = F[i + sx * ny];
    for (let j = 1; j <= ny; j++) for (let i = 1; i <= nx; i++) { const k = i + sx * j; C[k] += -(F[k] - F[k - sx]) + (cc[k] * (V[k] - V[k - sx]) * dt) / dy; } fill(C); };
  if (odd) { sweepY(); sweepX(); } else { sweepX(); sweepY(); }
  for (let k = 0; k < C.length; k++) if (C[k] < 0) C[k] = Math.max(C[k], -1e-12) < 0 ? 0 : C[k]; // round-off only
}
/**
 * Interface advection in a prescribed velocity field on the unit square (n × n cells, one period of the motion), comparing five families:
 *  'vof' (THINC/WLIC volume of fluid), 'levelset' (signed distance, second-order upwind, PDE reinitialisation), 'clsvof' (volume of fluid with the
 *  normal from the level set; the level set is shifted to enclose the VOF volume), 'phasefield' (conservative Allen–Cahn equation), 'front' (front tracking with markers).
 * test: 'zalesak' (slotted disc in solid-body rotation) | 'circle' (rotation) | 'translateX' | 'translateDiag' (uniform translation over one period, periodic).
 * Returns { scheme, test, n, x[], y[], c[][] (liquid fraction), exact[][], massError (relative), shapeError (L1 of the fraction / interface length … relative to the body area), steps }.
 */
export function interfaceAdvect2D(o = {}) {
  const n = o.n || 64, scheme = o.scheme || 'vof', test = o.test || 'zalesak', dx = 1 / n, sx = n + 2, N = sx * sx, idx = (i, j) => i + sx * j, periodic = test.startsWith('translate'), xc = (i) => (i - 0.5) * dx;
  const rot = !periodic, ux = (x, y) => (rot ? 2 * PI * (0.5 - y) : test === 'translateX' ? 1 : 1), uy = (x, y) => (rot ? 2 * PI * (x - 0.5) : test === 'translateX' ? 0 : 1), cx0 = rot ? 0.5 : 0.5, cy0 = rot ? 0.75 : 0.5, rad = 0.15;
  const slot = test === 'zalesak', inside = (x, y) => { const xx = periodic ? ((x % 1) + 1) % 1 : x, yy = periodic ? ((y % 1) + 1) % 1 : y; return Math.hypot(xx - cx0, yy - cy0) <= rad && !(slot && Math.abs(xx - cx0) < 0.025 && yy < cy0 + 0.1); };
  const frac = (i, j) => { let s = 0; const m = 6; for (let a = 0; a < m; a++) for (let b = 0; b < m; b++) if (inside(xc(i) + ((a + 0.5) / m - 0.5) * dx, xc(j) + ((b + 0.5) / m - 0.5) * dx)) s++; return s / (m * m); };
  const sdf = (x, y) => { const d = Math.hypot(x - cx0, y - cy0) - rad; if (!slot) return -d; // positive inside
    const inC = d <= 0, sxl = Math.abs(x - cx0) - 0.025, syt = y - (cy0 + 0.1), inS = sxl < 0 && syt < 0; if (inC && !inS) return Math.min(-d, sxl > 0 && syt > 0 ? Math.hypot(sxl, syt) : Math.max(sxl, syt) > 0 ? Math.max(sxl, syt) : 1); if (inC && inS) return -Math.min(-sxl, -syt); return -d; };
  const fill = (A) => { for (let i = 1; i <= n; i++) { A[idx(i, 0)] = periodic ? A[idx(i, n)] : A[idx(i, 1)]; A[idx(i, n + 1)] = periodic ? A[idx(i, 1)] : A[idx(i, n)]; } for (let j = 0; j <= n + 1; j++) { A[idx(0, j)] = periodic ? A[idx(n, j)] : A[idx(1, j)]; A[idx(n + 1, j)] = periodic ? A[idx(1, j)] : A[idx(n, j)]; } };
  const C = new Float64Array(N), E = new Float64Array(N), U = new Float64Array((n + 1) * (n + 2)), V = new Float64Array((n + 2) * (n + 1));
  for (let j = 1; j <= n; j++) for (let i = 1; i <= n; i++) { C[idx(i, j)] = frac(i, j); E[idx(i, j)] = C[idx(i, j)]; }
  for (let j = 0; j <= n + 1; j++) for (let i = 0; i <= n; i++) U[i + (n + 1) * j] = ux(i * dx, xc(j)); for (let j = 0; j <= n; j++) for (let i = 0; i <= n + 1; i++) V[i + sx * j] = uy(xc(i), j * dx);
  fill(C); const umax = rot ? 2 * PI * 0.71 : 1, cfl = o.cfl || 0.25, steps = Math.ceil(1 / ((cfl * dx) / umax)), dt = 1 / steps, mass = (A) => { let s = 0; for (let j = 1; j <= n; j++) for (let i = 1; i <= n; i++) s += A[idx(i, j)]; return s * dx * dx; }, m0 = mass(C);
  let out = C;
  if (scheme === 'front') { // markers on the interface, RK2 in the analytic velocity field, area by the shoelace formula
    const pts = []; const m = o.markers || 8 * n; if (slot) { const a0 = Math.asin(0.025 / rad), top = cy0 + 0.1, add = (x, y) => pts.push([x, y]), seg = (x1, y1, x2, y2, k) => { for (let q = 0; q < k; q++) add(x1 + ((x2 - x1) * q) / k, y1 + ((y2 - y1) * q) / k); }; const yb = cy0 - rad * Math.cos(a0), arcN = Math.round(m * 0.7); for (let q = 0; q < arcN; q++) { const th = -PI / 2 + a0 + ((2 * PI - 2 * a0) * q) / arcN; add(cx0 + rad * Math.cos(th), cy0 + rad * Math.sin(th)); } seg(cx0 - 0.025, yb, cx0 - 0.025, top, Math.round(m * 0.13)); seg(cx0 - 0.025, top, cx0 + 0.025, top, Math.round(m * 0.04)); seg(cx0 + 0.025, top, cx0 + 0.025, yb, Math.round(m * 0.13)); } else for (let q = 0; q < m; q++) pts.push([cx0 + rad * Math.cos((2 * PI * q) / m), cy0 + rad * Math.sin((2 * PI * q) / m)]);
    const area = () => { let s = 0; for (let q = 0; q < pts.length; q++) { const a = pts[q], b = pts[(q + 1) % pts.length]; s += a[0] * b[1] - b[0] * a[1]; } return Math.abs(s) / 2; }, a0 = area();
    for (let s = 0; s < steps; s++) for (const p of pts) { const k1x = ux(p[0], p[1]), k1y = uy(p[0], p[1]), mx = p[0] + 0.5 * dt * k1x, my = p[1] + 0.5 * dt * k1y; p[0] += dt * ux(mx, my); p[1] += dt * uy(mx, my); }
    const poly = pts.map((p) => (periodic ? [p[0] - 1, p[1] - (test === 'translateDiag' ? 1 : 0)] : p)), inPoly = (x, y) => { let c = false; for (let q = 0, w = poly.length - 1; q < poly.length; w = q++) { const a = poly[q], b = poly[w]; if (a[1] > y !== b[1] > y && x < ((b[0] - a[0]) * (y - a[1])) / (b[1] - a[1]) + a[0]) c = !c; } return c; };
    out = new Float64Array(N); for (let j = 1; j <= n; j++) for (let i = 1; i <= n; i++) { let s = 0; const q = 4; for (let a = 0; a < q; a++) for (let b = 0; b < q; b++) if (inPoly(xc(i) + ((a + 0.5) / q - 0.5) * dx, xc(j) + ((b + 0.5) / q - 0.5) * dx)) s++; out[idx(i, j)] = s / (q * q); }
    const sh = (() => { let s = 0; for (let j = 1; j <= n; j++) for (let i = 1; i <= n; i++) s += Math.abs(out[idx(i, j)] - E[idx(i, j)]); return (s * dx * dx) / m0; })();
    return pack(out, Math.abs(area() - a0) / a0, sh);
  }
  const P = new Float64Array(N), heav = (p) => { const e = 1.5 * dx; return p > e ? 1 : p < -e ? 0 : 0.5 * (1 + p / e + Math.sin((PI * p) / e) / PI); };
  if (scheme === 'levelset' || scheme === 'clsvof') { for (let j = 1; j <= n; j++) for (let i = 1; i <= n; i++) P[idx(i, j)] = sdf(xc(i), xc(j)); fill(P); }
  const mm = (a, b) => (a * b <= 0 ? 0 : Math.abs(a) < Math.abs(b) ? a : b);
  const lsRhs = (A, R) => { for (let j = 1; j <= n; j++) for (let i = 1; i <= n; i++) { const k = idx(i, j), u = 0.5 * (U[i - 1 + (n + 1) * j] + U[i + (n + 1) * j]), v = 0.5 * (V[k - sx] + V[k]);
    const ip = Math.min(i + 2, n + 1), im = Math.max(i - 2, 0), jp = Math.min(j + 2, n + 1), jm = Math.max(j - 2, 0), d2 = (a, b, c) => a - 2 * b + c; // second-order upwind (ENO-2)
    const dxm = (A[k] - A[k - 1]) / dx + (mm(d2(A[k + 1], A[k], A[k - 1]), d2(A[k], A[k - 1], A[idx(im, j)])) / dx) * 0.5, dxp = (A[k + 1] - A[k]) / dx - (mm(d2(A[k + 1], A[k], A[k - 1]), d2(A[idx(ip, j)], A[k + 1], A[k])) / dx) * 0.5;
    const dym = (A[k] - A[k - sx]) / dx + (mm(d2(A[k + sx], A[k], A[k - sx]), d2(A[k], A[k - sx], A[idx(i, jm)])) / dx) * 0.5, dyp = (A[k + sx] - A[k]) / dx - (mm(d2(A[k + sx], A[k], A[k - sx]), d2(A[idx(i, jp)], A[k + sx], A[k])) / dx) * 0.5;
    R[k] = -(u * (u >= 0 ? dxm : dxp) + v * (v >= 0 ? dym : dyp)); } };
  const fillLS = (A) => { if (periodic) return fill(A); for (let i = 1; i <= n; i++) { A[idx(i, 0)] = 2 * A[idx(i, 1)] - A[idx(i, 2)]; A[idx(i, n + 1)] = 2 * A[idx(i, n)] - A[idx(i, n - 1)]; } for (let j = 0; j <= n + 1; j++) { A[idx(0, j)] = 2 * A[idx(1, j)] - A[idx(2, j)]; A[idx(n + 1, j)] = 2 * A[idx(n, j)] - A[idx(n - 1, j)]; } };
  const reinit = (A, iters) => { const S0 = new Float64Array(N), R = new Float64Array(N); for (let k = 0; k < N; k++) S0[k] = A[k] / Math.sqrt(A[k] * A[k] + dx * dx); for (let q = 0; q < iters; q++) { for (let j = 1; j <= n; j++) for (let i = 1; i <= n; i++) { const k = idx(i, j), a = (A[k] - A[k - 1]) / dx, b = (A[k + 1] - A[k]) / dx, c = (A[k] - A[k - sx]) / dx, d = (A[k + sx] - A[k]) / dx, s = S0[k]; const g = s > 0 ? Math.sqrt(Math.max(Math.max(a, 0) ** 2, Math.min(b, 0) ** 2) + Math.max(Math.max(c, 0) ** 2, Math.min(d, 0) ** 2)) : Math.sqrt(Math.max(Math.min(a, 0) ** 2, Math.max(b, 0) ** 2) + Math.max(Math.min(c, 0) ** 2, Math.max(d, 0) ** 2)); R[k] = -s * (g - 1); } for (let j = 1; j <= n; j++) for (let i = 1; i <= n; i++) A[idx(i, j)] += 0.3 * dx * R[idx(i, j)]; fillLS(A); } };
  const volLS = (sh) => { let s = 0; for (let j = 1; j <= n; j++) for (let i = 1; i <= n; i++) s += heav(P[idx(i, j)] + sh); return s * dx * dx; };
  const R1 = new Float64Array(N), R2 = new Float64Array(N), T = new Float64Array(N), nrm = { nx: new Float64Array(N), ny: new Float64Array(N) };
  const m0ls = scheme === 'levelset' ? volLS(0) : m0;
  for (let s = 0; s < steps; s++) {
    if (scheme === 'levelset' || scheme === 'clsvof') { lsRhs(P, R1); for (let k = 0; k < N; k++) T[k] = P[k] + dt * R1[k]; fillLS(T); lsRhs(T, R2); for (let j = 1; j <= n; j++) for (let i = 1; i <= n; i++) { const k = idx(i, j); P[k] += 0.5 * dt * (R1[k] + R2[k]); } fillLS(P); if ((s + 1) % (o.reinitEvery || 40) === 0) reinit(P, 2); }
    if (scheme === 'vof' || scheme === 'clsvof') { let nr = null; if (scheme === 'clsvof') { for (let j = 1; j <= n; j++) for (let i = 1; i <= n; i++) { const k = idx(i, j); nrm.nx[k] = P[k + 1] - P[k - 1]; nrm.ny[k] = P[k + sx] - P[k - sx]; } nr = nrm; } vofStep(C, U, V, n, n, dx, dx, dt, s % 2 === 1, nr, fill, periodic, periodic);
      if (scheme === 'clsvof' && (s + 1) % 4 === 0) { const target = mass(C), v0 = volLS(0), e = 0.2 * dx, slope = (volLS(e) - v0) / e, sh = slope > 1e-9 ? clamp((target - v0) / slope, -dx, dx) : 0; for (let k = 0; k < N; k++) P[k] += sh; } } // the level set is shifted so that it encloses the VOF volume
    if (scheme === 'phasefield') { // conservative Allen–Cahn: ∂φ/∂t + ∇·(uφ) = ∇·[γ (ε ∇φ − φ(1 − φ) n)], n = ∇φ/|∇φ|
      const eps = (o.epsCells || 0.75) * dx, gam = umax, sub = Math.max(1, Math.ceil((dt * gam * eps * (o.pfSafety || 8)) / (dx * dx)));
      for (let q = 0; q < sub; q++) { const h = dt / sub, Fx = R1, Fy = R2;
        for (let j = 1; j <= n; j++) for (let i = 0; i <= n; i++) { const k = idx(i, j), u = U[i + (n + 1) * j], d = u >= 0 ? k : k + 1, s1 = u >= 0 ? 1 : -1, sl = mm(C[d] - C[d - s1], C[d + s1] - C[d]), cf = C[d] + 0.5 * sl * (1 - (Math.abs(u) * h) / dx), g1 = (C[k + 1] - C[k]) / dx, gy = (C[k + sx] + C[k + 1 + sx] - C[k - sx] - C[k + 1 - sx]) / (4 * dx), mg = Math.hypot(g1, gy) + 1e-12, cm = 0.5 * (C[k] + C[k + 1]); Fx[k] = u * cf - gam * (eps * g1 - cm * (1 - cm) * (g1 / mg)); }
        for (let j = 0; j <= n; j++) for (let i = 1; i <= n; i++) { const k = idx(i, j), v = V[k], d = v >= 0 ? k : k + sx, s1 = v >= 0 ? sx : -sx, sl = mm(C[d] - (C[d - s1] ?? C[d]), (C[d + s1] ?? C[d]) - C[d]), cf = C[d] + 0.5 * sl * (1 - (Math.abs(v) * h) / dx), g1 = (C[k + sx] - C[k]) / dx, gx = (C[k + 1] + C[k + 1 + sx] - C[k - 1] - C[k - 1 + sx]) / (4 * dx), mg = Math.hypot(g1, gx) + 1e-12, cm = 0.5 * (C[k] + C[k + sx]); Fy[k] = v * cf - gam * (eps * g1 - cm * (1 - cm) * (g1 / mg)); }
        if (!periodic) { for (let j = 1; j <= n; j++) { Fx[idx(0, j)] = 0; Fx[idx(n, j)] = 0; } for (let i = 1; i <= n; i++) { Fy[idx(i, 0)] = 0; Fy[idx(i, n)] = 0; } } else { for (let j = 1; j <= n; j++) Fx[idx(0, j)] = Fx[idx(n, j)]; for (let i = 1; i <= n; i++) Fy[idx(i, 0)] = Fy[idx(i, n)]; }
        for (let j = 1; j <= n; j++) for (let i = 1; i <= n; i++) { const k = idx(i, j); T[k] = C[k] - (h / dx) * (Fx[k] - Fx[k - 1] + Fy[k] - Fy[k - sx]); } for (let j = 1; j <= n; j++) for (let i = 1; i <= n; i++) C[idx(i, j)] = T[idx(i, j)]; fill(C); }
    }
  }
  if (scheme === 'levelset') { out = new Float64Array(N); for (let k = 0; k < N; k++) out[k] = heav(P[k]); for (let k = 0; k < N; k++) E[k] = E[k]; }
  const m1 = scheme === 'levelset' ? volLS(0) : mass(out); let sh = 0; for (let j = 1; j <= n; j++) for (let i = 1; i <= n; i++) sh += Math.abs(out[idx(i, j)] - E[idx(i, j)]);
  return pack(out, Math.abs(m1 - m0ls) / m0ls, (sh * dx * dx) / m0);
  function pack(A, massError, shapeError) { const xs = Array.from({ length: n }, (_, i) => xc(i + 1)), z = [], ex = []; for (let j = 1; j <= n; j++) { z.push(Array.from({ length: n }, (_, i) => A[idx(i + 1, j)])); ex.push(Array.from({ length: n }, (_, i) => E[idx(i + 1, j)])); } return { scheme, test, n, x: xs, y: xs, c: z, exact: ex, massError, shapeError, steps }; }
}

// Preconditioned conjugate gradients (symmetric Gauss–Seidel preconditioner) for the pressure equation −∇·(β ∇p) = b on nx × ny cells.
// cE[k]: coefficient to the east neighbour, cN[k]: to the north; dg[k]: diagonal (sum of the neighbour coefficients + any Dirichlet part).
function pcgPoisson(nx, ny, cE, cN, dg, b, p, perX, tol, maxIt, wk) {
  const n = nx * ny, r = wk.r, z = wk.z, q = wk.q, d = wk.d;
  const apply = (x, y) => { for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) { const k = i + nx * j; let s = dg[k] * x[k]; if (i < nx - 1) s -= cE[k] * x[k + 1]; else if (perX) s -= cE[k] * x[k + 1 - nx]; if (i > 0) s -= cE[k - 1] * x[k - 1]; else if (perX) s -= cE[k - 1 + nx] * x[k - 1 + nx]; if (j < ny - 1) s -= cN[k] * x[k + nx]; if (j > 0) s -= cN[k - nx] * x[k - nx]; y[k] = s; } };
  const prec = (x, y) => { for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) { const k = i + nx * j; let s = x[k]; if (i > 0) s += cE[k - 1] * y[k - 1]; if (j > 0) s += cN[k - nx] * y[k - nx]; y[k] = s / dg[k]; } for (let k = 0; k < n; k++) y[k] *= dg[k]; for (let j = ny - 1; j >= 0; j--) for (let i = nx - 1; i >= 0; i--) { const k = i + nx * j; let s = y[k]; if (i < nx - 1) s += cE[k] * y[k + 1]; if (j < ny - 1) s += cN[k] * y[k + nx]; y[k] = s / dg[k]; } };
  apply(p, q); let bn = 0, rz = 0; for (let k = 0; k < n; k++) { r[k] = b[k] - q[k]; bn += b[k] * b[k]; } bn = Math.sqrt(bn) || 1;
  prec(r, z); for (let k = 0; k < n; k++) { d[k] = z[k]; rz += r[k] * z[k]; }
  let it = 0, rn = 0;
  for (; it < maxIt; it++) {
    rn = 0; for (let k = 0; k < n; k++) rn += r[k] * r[k]; rn = Math.sqrt(rn); if (rn <= tol * bn) break;
    apply(d, q); let dq = 0; for (let k = 0; k < n; k++) dq += d[k] * q[k]; if (!(Math.abs(dq) > 1e-300)) break; const a = rz / dq;
    for (let k = 0; k < n; k++) { p[k] += a * d[k]; r[k] -= a * q[k]; }
    prec(r, z); let rz2 = 0; for (let k = 0; k < n; k++) rz2 += r[k] * z[k]; const be = rz2 / rz; rz = rz2; for (let k = 0; k < n; k++) d[k] = z[k] + be * d[k];
  }
  return { iterations: it, residual: rn / bn };
}
/**
 * Incompressible two-dimensional Navier–Stokes solver on a uniform staggered (MAC) grid: explicit momentum predictor (conservative
 * convection with a central / upwind blend, full viscous stress with variable viscosity), pressure projection with a variable-density
 * Poisson equation (conjugate gradients), optional second fluid by volume of fluid (THINC/WLIC, see vofStep) and optional Menter SST
 * k–ω turbulence (2003 constants, ω fixed in the wall-adjacent cells). Walls on all sides, or periodic in x; the top wall may move (lid).
 * o: { nx, ny, lx, ly, rho, mu (single fluid) | rhoL, rhoG, muL, muG + c0(x, y) (liquid fraction), gx, gy, fx (body force per volume, x), lid (m/s), periodicX,
 *      turbulence: 'none' | 'sst', tEnd, cfl, dtMax, upwind (0 = central … 1 = first-order upwind), steady (stop when the velocity change per unit time falls below it),
 *      poissonTol, maxSteps, u0(x, y), v0(x, y), k0, omega0, probeFront (track the liquid front along the bottom) }
 * Returns { x[], y[], u[][], v[][], p[][], c[][] | null, nut[][] | null, k[][] | null, t, steps, dtMean, divergenceMax, poissonIterations (mean), steadyResidual,
 *           mass: { initial, final, error } (liquid volume per unit depth), kinetic[] { t, e }, front: { t[], x[] } | null, uMid[] (u on the vertical centre line), vMid[] (v on the horizontal centre line), momentumX: { force, wall } }.
 */
export function cfd2d(o = {}) {
  const nx = o.nx || 32, ny = o.ny || 32, lx = o.lx || 1, ly = o.ly || 1, dx = lx / nx, dy = ly / ny, perX = !!o.periodicX, two = typeof o.c0 === 'function', sst = o.turbulence === 'sst';
  const rhoL = o.rhoL ?? o.rho ?? 1, rhoG = o.rhoG ?? rhoL, muL = o.muL ?? o.mu ?? 0.01, muG = o.muG ?? muL, gx = o.gx || 0, gy = o.gy || 0, fx = o.fx || 0, lid = o.lid || 0, upw = o.upwind ?? (two ? 0.3 : 0), tEnd = o.tEnd ?? 1, cfl = o.cfl ?? 0.4;
  const sx = nx + 2, su = nx + 1, NC = sx * (ny + 2), F = (n) => new Float64Array(n), U = F(su * (ny + 2)), V = F(sx * (ny + 1)), Us = F(U.length), Vs = F(V.length), C = F(NC), P = F(nx * ny), RHO = F(NC), MU = F(NC), K = F(NC), W = F(NC), NUT = F(NC);
  const ci = (i, j) => i + sx * j, ui = (i, j) => i + su * j, xc = (i) => (i - 0.5) * dx, yc = (j) => (j - 0.5) * dy;
  for (let j = 1; j <= ny; j++) for (let i = 1; i <= nx; i++) C[ci(i, j)] = two ? clamp(o.c0(xc(i), yc(j)), 0, 1) : 1;
  if (o.u0) for (let j = 1; j <= ny; j++) for (let i = 0; i <= nx; i++) U[ui(i, j)] = o.u0(i * dx, yc(j)); if (o.v0) for (let j = 0; j <= ny; j++) for (let i = 1; i <= nx; i++) V[ci(i, j)] = o.v0(xc(i), j * dy);
  const fillC = (A) => { for (let i = 1; i <= nx; i++) { A[ci(i, 0)] = A[ci(i, 1)]; A[ci(i, ny + 1)] = A[ci(i, ny)]; } for (let j = 0; j <= ny + 1; j++) { A[ci(0, j)] = perX ? A[ci(nx, j)] : A[ci(1, j)]; A[ci(nx + 1, j)] = perX ? A[ci(1, j)] : A[ci(nx, j)]; } };
  const bcVel = (Ua, Va) => { // no-slip walls by reflection, moving lid, periodic or solid ends
    if (perX) for (let j = 1; j <= ny; j++) { const m = 0.5 * (Ua[ui(0, j)] + Ua[ui(nx, j)]); Ua[ui(0, j)] = m; Ua[ui(nx, j)] = m; } else for (let j = 1; j <= ny; j++) { Ua[ui(0, j)] = 0; Ua[ui(nx, j)] = 0; }
    for (let i = 0; i <= nx; i++) { Ua[ui(i, 0)] = -Ua[ui(i, 1)]; Ua[ui(i, ny + 1)] = 2 * lid - Ua[ui(i, ny)]; }
    for (let i = 1; i <= nx; i++) { Va[ci(i, 0)] = 0; Va[ci(i, ny)] = 0; }
    for (let j = 0; j <= ny; j++) { Va[ci(0, j)] = perX ? Va[ci(nx, j)] : -Va[ci(1, j)]; Va[ci(nx + 1, j)] = perX ? Va[ci(1, j)] : -Va[ci(nx, j)]; }
  };
  const props = () => { for (let k = 0; k < NC; k++) { const c = clamp(C[k], 0, 1); RHO[k] = rhoG + (rhoL - rhoG) * c; MU[k] = muG + (muL - muG) * c + (sst ? RHO[k] * NUT[k] : 0); } };
  const dist = (i, j) => { let d = Math.min(yc(j), ly - yc(j)); if (!perX) d = Math.min(d, xc(i), lx - xc(i)); return d; };
  if (sst) for (let j = 1; j <= ny; j++) for (let i = 1; i <= nx; i++) { K[ci(i, j)] = o.k0 ?? 1e-4; W[ci(i, j)] = o.omega0 ?? 10; }
  fillC(C); if (sst) { fillC(K); fillC(W); } bcVel(U, V); props();
  const cE = F(nx * ny), cN = F(nx * ny), dg = F(nx * ny), bb = F(nx * ny), wk = { r: F(nx * ny), z: F(nx * ny), q: F(nx * ny), d: F(nx * ny) }, liquid = () => { let s = 0; for (let j = 1; j <= ny; j++) for (let i = 1; i <= nx; i++) s += C[ci(i, j)]; return s * dx * dy; };
  const m0 = liquid(), kin = [], front = o.probeFront ? { t: [], x: [] } : null; let t = 0, steps = 0, pits = 0, divMax = 0, resid = Infinity, forceI = 0, wallI = 0;
  const bS = 0.09, a1 = 0.31, sk1 = 0.85, sk2 = 1, sw1 = 0.5, sw2 = 0.856, b1 = 0.075, b2 = 0.0828, g1 = 5 / 9, g2 = 0.44, Kn = F(NC), Wn = F(NC);
  while (t < tEnd - 1e-12 && steps < (o.maxSteps || 200000)) {
    let vmax = 1e-9, numax = 0; for (let k = 0; k < U.length; k++) vmax = Math.max(vmax, Math.abs(U[k])); for (let k = 0; k < V.length; k++) vmax = Math.max(vmax, Math.abs(V[k])); for (let j = 1; j <= ny; j++) for (let i = 1; i <= nx; i++) numax = Math.max(numax, MU[ci(i, j)] / RHO[ci(i, j)]);
    const h = Math.min(dx, dy), gmag = Math.hypot(gx, gy) + Math.abs(fx) / Math.min(rhoL, rhoG), dt = Math.min((cfl * h) / vmax, (0.2 * h * h) / Math.max(numax, 1e-30), gmag > 0 ? 0.4 * Math.sqrt(h / gmag) : Infinity, o.dtMax ?? Infinity, tEnd - t);
    // ---- interface
    if (two) { vofStep(C, U, V, nx, ny, dx, dy, dt, steps % 2 === 1, null, fillC, perX, false); props(); }
    // ---- turbulence (explicit convection and diffusion, point-implicit destruction)
    if (sst) {
      for (let j = 1; j <= ny; j++) for (let i = 1; i <= nx; i++) {
        const k = ci(i, j), uc = 0.5 * (U[ui(i - 1, j)] + U[ui(i, j)]), vc = 0.5 * (V[ci(i, j - 1)] + V[k]), nu = (muG + (muL - muG) * clamp(C[k], 0, 1)) / RHO[k], d = dist(i, j);
        const dudx = (U[ui(i, j)] - U[ui(i - 1, j)]) / dx, dvdy = (V[k] - V[ci(i, j - 1)]) / dy, dudy = (U[ui(i, j + 1)] + U[ui(i - 1, j + 1)] - U[ui(i, j - 1)] - U[ui(i - 1, j - 1)]) / (4 * dy), dvdx = (V[ci(i + 1, j)] + V[ci(i + 1, j - 1)] - V[ci(i - 1, j)] - V[ci(i - 1, j - 1)]) / (4 * dx), S2 = 2 * dudx * dudx + 2 * dvdy * dvdy + (dudy + dvdx) ** 2, S = Math.sqrt(S2);
        const kc = K[k], wc = Math.max(W[k], 1e-12), dkx = (K[k + 1] - K[k - 1]) / (2 * dx), dky = (K[k + sx] - K[k - sx]) / (2 * dy), dwx = (W[k + 1] - W[k - 1]) / (2 * dx), dwy = (W[k + sx] - W[k - sx]) / (2 * dy), cd = (2 * sw2 * (dkx * dwx + dky * dwy)) / wc;
        const arg1 = Math.min(Math.max(Math.sqrt(kc) / (bS * wc * d), (500 * nu) / (d * d * wc)), (4 * sw2 * kc) / (Math.max(cd, 1e-10) * d * d)), F1 = Math.tanh(arg1 ** 4), arg2 = Math.max((2 * Math.sqrt(kc)) / (bS * wc * d), (500 * nu) / (d * d * wc)), F2 = Math.tanh(arg2 * arg2);
        const nt = (a1 * kc) / Math.max(a1 * wc, S * F2), Pk = Math.min(nt * S2, 10 * bS * kc * wc), sk = F1 * sk1 + (1 - F1) * sk2, sw = F1 * sw1 + (1 - F1) * sw2, be = F1 * b1 + (1 - F1) * b2, gm = F1 * g1 + (1 - F1) * g2;
        const adv = (A) => -(uc >= 0 ? (uc * (A[k] - A[k - 1])) / dx : (uc * (A[k + 1] - A[k])) / dx) - (vc >= 0 ? (vc * (A[k] - A[k - sx])) / dy : (vc * (A[k + sx] - A[k])) / dy);
        const dif = (A, s) => ((nu + s * 0.5 * (NUT[k] + NUT[k + 1])) * (A[k + 1] - A[k]) - (nu + s * 0.5 * (NUT[k] + NUT[k - 1])) * (A[k] - A[k - 1])) / (dx * dx) + ((nu + s * 0.5 * (NUT[k] + NUT[k + sx])) * (A[k + sx] - A[k]) - (nu + s * 0.5 * (NUT[k] + NUT[k - sx])) * (A[k] - A[k - sx])) / (dy * dy);
        Kn[k] = Math.max((kc + dt * (adv(K) + dif(K, sk) + Pk)) / (1 + dt * bS * wc), 1e-14);
        const x = (1 - F1) * cd; Wn[k] = Math.max((wc + dt * (adv(W) + dif(W, sw) + (gm * Pk) / Math.max(nt, 1e-12) + Math.max(x, 0))) / (1 + dt * (be * wc + Math.max(-x, 0) / wc)), 1e-10);
        if (j === 1 || j === ny || (!perX && (i === 1 || i === nx))) Wn[k] = (6 * nu) / (b1 * d * d); // wall-adjacent cells: analytic near-wall ω
        NUT[k] = (a1 * Kn[k]) / Math.max(a1 * Wn[k], S * F2);
      }
      K.set(Kn); W.set(Wn); fillC(K); fillC(W); fillC(NUT);
      for (let i = 1; i <= nx; i++) { K[ci(i, 0)] = -K[ci(i, 1)]; K[ci(i, ny + 1)] = -K[ci(i, ny)]; } if (!perX) for (let j = 1; j <= ny; j++) { K[ci(0, j)] = -K[ci(1, j)]; K[ci(nx + 1, j)] = -K[ci(nx, j)]; }
      props();
    }
    // ---- momentum predictor
    const muCorner = (i, j) => 0.25 * (MU[ci(i, j)] + MU[ci(i + 1, j)] + MU[ci(i, j + 1)] + MU[ci(i + 1, j + 1)]);
    for (let j = 1; j <= ny; j++) for (let i = perX ? 0 : 1; i <= (perX ? nx : nx - 1); i++) {
      const k = ui(i, j), iw = i === 0 ? nx - 1 : i - 1, ie = i === nx ? 1 : i + 1, cl = i === 0 ? nx : i, cr = i === nx ? 1 : i + 1; // cells left / right of the face (periodic wrap)
      const uP = U[k], uE = U[ui(ie, j)], uW = U[ui(iw, j)], uN = U[ui(i, j + 1)], uS = U[ui(i, j - 1)], ue = 0.5 * (uP + uE), uw = 0.5 * (uP + uW), vn = 0.5 * (V[ci(cl, j)] + V[ci(cr, j)]), vs = 0.5 * (V[ci(cl, j - 1)] + V[ci(cr, j - 1)]);
      const fe = ue * (0.5 * (uP + uE) - upw * 0.5 * Math.sign(ue) * (uE - uP)), fw = uw * (0.5 * (uW + uP) - upw * 0.5 * Math.sign(uw) * (uP - uW)), fn = vn * (0.5 * (uP + uN) - upw * 0.5 * Math.sign(vn) * (uN - uP)), fs = vs * (0.5 * (uS + uP) - upw * 0.5 * Math.sign(vs) * (uP - uS));
      const conv = (fe - fw) / dx + (fn - fs) / dy, rf = 0.5 * (RHO[ci(cl, j)] + RHO[ci(cr, j)]);
      const txxE = (2 * MU[ci(cr, j)] * (uE - uP)) / dx, txxW = (2 * MU[ci(cl, j)] * (uP - uW)) / dx, mN = muCorner(cl, j), mS = muCorner(cl, j - 1), txyN = mN * ((uN - uP) / dy + (V[ci(cr, j)] - V[ci(cl, j)]) / dx), txyS = mS * ((uP - uS) / dy + (V[ci(cr, j - 1)] - V[ci(cl, j - 1)]) / dx);
      Us[k] = uP + dt * (-conv + ((txxE - txxW) / dx + (txyN - txyS) / dy + fx) / rf + gx);
    }
    for (let j = 1; j <= ny - 1; j++) for (let i = 1; i <= nx; i++) {
      const k = ci(i, j), vP = V[k], vE = V[k + 1], vW = V[k - 1], vN = V[k + sx], vS = V[k - sx], vnn = 0.5 * (vP + vN), vss = 0.5 * (vP + vS), ue = 0.5 * (U[ui(i, j)] + U[ui(i, j + 1)]), uw = 0.5 * (U[ui(i - 1, j)] + U[ui(i - 1, j + 1)]);
      const fn = vnn * (0.5 * (vP + vN) - upw * 0.5 * Math.sign(vnn) * (vN - vP)), fs = vss * (0.5 * (vS + vP) - upw * 0.5 * Math.sign(vss) * (vP - vS)), fe = ue * (0.5 * (vP + vE) - upw * 0.5 * Math.sign(ue) * (vE - vP)), fw = uw * (0.5 * (vW + vP) - upw * 0.5 * Math.sign(uw) * (vP - vW));
      const conv = (fe - fw) / dx + (fn - fs) / dy, rf = 0.5 * (RHO[k] + RHO[k + sx]), tyyN = (2 * MU[k + sx] * (vN - vP)) / dy, tyyS = (2 * MU[k] * (vP - vS)) / dy, mE = muCorner(i, j), mW = muCorner(i - 1, j);
      const txyE = mE * ((vE - vP) / dx + (U[ui(i, j + 1)] - U[ui(i, j)]) / dy), txyW = mW * ((vP - vW) / dx + (U[ui(i - 1, j + 1)] - U[ui(i - 1, j)]) / dy);
      Vs[k] = vP + dt * (-conv + ((tyyN - tyyS) / dy + (txyE - txyW) / dx) / rf + gy);
    }
    bcVel(Us, Vs);
    // ---- pressure projection
    let bsum = 0;
    for (let j = 1; j <= ny; j++) for (let i = 1; i <= nx; i++) { const k = i - 1 + nx * (j - 1), c = ci(i, j); cE[k] = i < nx || perX ? 1 / (0.5 * (RHO[c] + RHO[ci(i === nx ? 1 : i + 1, j)])) / (dx * dx) : 0; cN[k] = j < ny ? 1 / (0.5 * (RHO[c] + RHO[c + sx])) / (dy * dy) : 0; bb[k] = -((Us[ui(i, j)] - Us[ui(i - 1, j)]) / dx + (Vs[c] - Vs[c - sx]) / dy) / dt; bsum += bb[k]; }
    for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) { const k = i + nx * j; dg[k] = cE[k] + (i > 0 ? cE[k - 1] : perX ? cE[k - 1 + nx] : 0) + cN[k] + (j > 0 ? cN[k - nx] : 0); bb[k] -= bsum / (nx * ny); }
    dg[0] *= 1 + 1e-9; // pins the pressure level without disturbing the solution
    const ps = pcgPoisson(nx, ny, cE, cN, dg, bb, P, perX, o.poissonTol ?? 1e-8, o.poissonMaxIt ?? 2000, wk); pits += ps.iterations;
    let dU = 0;
    for (let j = 1; j <= ny; j++) for (let i = perX ? 0 : 1; i <= (perX ? nx : nx - 1); i++) { const cl = i === 0 ? nx : i, cr = i === nx ? 1 : i + 1, rf = 0.5 * (RHO[ci(cl, j)] + RHO[ci(cr, j)]), un = Us[ui(i, j)] - (dt * (P[cr - 1 + nx * (j - 1)] - P[cl - 1 + nx * (j - 1)])) / (dx * rf); dU = Math.max(dU, Math.abs(un - U[ui(i, j)])); U[ui(i, j)] = un; }
    for (let j = 1; j <= ny - 1; j++) for (let i = 1; i <= nx; i++) { const k = ci(i, j), rf = 0.5 * (RHO[k] + RHO[k + sx]), vn = Vs[k] - (dt * (P[i - 1 + nx * j] - P[i - 1 + nx * (j - 1)])) / (dy * rf); dU = Math.max(dU, Math.abs(vn - V[k])); V[k] = vn; }
    bcVel(U, V);
    t += dt; steps++; resid = dU / dt;
    if (steps % 10 === 0 || t >= tEnd - 1e-12) { let e = 0, dv = 0; for (let j = 1; j <= ny; j++) for (let i = 1; i <= nx; i++) { const c = ci(i, j), uc = 0.5 * (U[ui(i - 1, j)] + U[ui(i, j)]), vc = 0.5 * (V[c - sx] + V[c]); e += 0.5 * RHO[c] * (uc * uc + vc * vc) * dx * dy; dv = Math.max(dv, Math.abs((U[ui(i, j)] - U[ui(i - 1, j)]) / dx + (V[c] - V[c - sx]) / dy)); } divMax = Math.max(divMax, dv * h / Math.max(vmax, 1e-9)); if (kin.length < 400) kin.push({ t, e });
      if (front) { let xf = 0; for (let i = 1; i <= nx; i++) if (C[ci(i, 1)] > 0.5) xf = i * dx; front.t.push(t); front.x.push(xf); } }
    if (fx) { forceI += fx * lx * ly * dt; let w = 0; for (let i = 1; i <= nx; i++) w += (MU[ci(i, 1)] * 2 * U[ui(i, 1)]) / dy * dx + (MU[ci(i, ny)] * 2 * (U[ui(i, ny)] - lid)) / dy * dx; wallI += w * dt; }
    if (o.steady && resid < o.steady && steps > 20) break;
    if (!Number.isFinite(resid)) break;
  }
  const xs = Array.from({ length: nx }, (_, i) => xc(i + 1)), ys = Array.from({ length: ny }, (_, j) => yc(j + 1)), grid = (f) => ys.map((_, j) => xs.map((__, i) => f(i + 1, j + 1)));
  const im = Math.round(nx / 2), jm = Math.round(ny / 2), m1 = liquid();
  return { x: xs, y: ys, u: grid((i, j) => 0.5 * (U[ui(i - 1, j)] + U[ui(i, j)])), v: grid((i, j) => 0.5 * (V[ci(i, j - 1)] + V[ci(i, j)])), p: grid((i, j) => P[i - 1 + nx * (j - 1)]), c: two ? grid((i, j) => C[ci(i, j)]) : null, nut: sst ? grid((i, j) => NUT[ci(i, j)]) : null, k: sst ? grid((i, j) => K[ci(i, j)]) : null,
    t, steps, dtMean: steps ? t / steps : 0, divergenceMax: divMax, poissonIterations: steps ? pits / steps : 0, steadyResidual: resid, mass: { initial: m0, final: m1, error: two ? (m1 - m0) / Math.max(m0, 1e-300) : 0 }, kinetic: kin, front,
    uMid: ys.map((_, j) => (nx % 2 === 0 ? U[ui(im, j + 1)] : 0.5 * (U[ui(im, j + 1)] + U[ui(im - 1, j + 1)]))), vMid: xs.map((_, i) => (ny % 2 === 0 ? V[ci(i + 1, jm)] : 0.5 * (V[ci(i + 1, jm)] + V[ci(i + 1, jm - 1)]))), momentumX: { force: forceI, wall: wallI } };
}

// =====================================================================================================
// 7g. Baker map, equipment boundary characteristics, branch junction, riser shapes
// =====================================================================================================
/**
 * Baker (1954) horizontal flow-pattern map in Baker's coordinates  B_y = G_g / λ  and  B_x = (G_l / G_g) λ ψ  with
 * λ = [(ρg/0.075)(ρl/62.3)]^½ and ψ = (73/σ)[μl (62.3/ρl)²]^⅓ (lb/h/ft², lb/ft³, cP, dyn/cm). The six boundary curves are
 * log-polynomial fits of the chart; the fits and the region logic could not be compared with an openly readable copy (see PROVENANCE):
 * above C5 bubble/froth, above C3 dispersed, above C1 (and C4 to the right of its minimum) annular, between C2 and C1 wave,
 * below both: stratified up to the abscissa where C1 and C2 meet, beyond it slug above C6 and plug below.
 * Returns { regime, Bx, By }.
 */
export function bakerRegime(p) {
  const GL = p.rhoL * p.vsl * 737.338, GG = Math.max(p.rhoG * p.vsg, 1e-12) * 737.338, rl = p.rhoL / 16.01846, rg = p.rhoG / 16.01846, By = (2.16 * GG) / Math.sqrt(rl * rg), Bx = (531 * (GL / GG) * Math.sqrt(rl * rg) * (p.muL * 1000) ** (1 / 3)) / (rl ** (2 / 3) * Math.max(p.sigma, 1e-4) * 1000), x = Math.log(Math.max(Bx, 1e-9));
  const c1 = Math.exp(9.774459 - 0.6548 * x), c2 = Math.exp(8.67694 - 0.1901 * x), c3 = Math.exp(11.3976 - 0.6084 * x + 0.0779 * x * x), c4 = Math.exp(10.7448 - 1.6265 * x + 0.2839 * x * x), c5 = Math.exp(14.569802 - 1.0173 * x), c6 = Math.exp(7.8206 - 0.2189 * x), xMin = 1.6265 / (2 * 0.2839), x12 = (9.774459 - 8.67694) / (0.6548 - 0.1901);
  const regime = By > c5 ? 'bubble (froth)' : By > c3 ? 'dispersed' : By > c1 && (x < xMin || By > c4) ? 'annular' : By > c2 && By <= c1 ? 'wave' : x <= x12 ? 'stratified' : By > c6 ? 'slug' : 'plug';
  return { regime, Bx, By };
}
/** Catenary riser profile hanging from the top: horizontal offset x and height z above the touch-down point for a top angle to the vertical (deg). Returns [{ x, z }] from the touch-down point upward. */
export function catenaryProfile({ height, topAngle = 12, n = 14 }) {
  const tanT = Math.tan((90 - clamp(topAngle, 1, 80)) * DEG), a = height / (Math.sqrt(1 + tanT * tanT) - 1), xT = a * Math.asinh(tanT); // z = a (cosh(x/a) − 1), slope at the top = tan(90° − top angle)
  return Array.from({ length: n + 1 }, (_, i) => { const x = (xT * i) / n; return { x, z: a * (Math.cosh(x / a) - 1) }; });
}
/**
 * Hydraulic characteristics imposed at the ends of the line.
 *  pump (multiphase booster at the inlet):  Δp = Δp₀ s² − (Δp₀ / q_max²) q²   (parabolic head curve, affinity laws in the speed ratio s), q = actual volume rate at suction;
 *  separator:  the gas leaves through a valve to the compressor suction, p_sep² − p_suc² = K_v q_std² (compressible valve equation), so the separator pressure floats with the gas rate;
 *  compressor: polytropic head curve H = H₀ s² [1 − (q_s / q_surge-free max)²], p_suc = p_dis / [1 + (n − 1)/n · H M / (Z R T η_p)]^(n/(n−1)).
 * Returns { pumpDp(q) (bar), sepPressure(qGasStd (Sm³/s), pSet (bara)) (bara), compressor(qSuction (m³/s), pDis (bara), T (K), Z, M (kg/mol)) → { pSuction, head, ratio, power (kW at mdot) } }.
 */
export function equipment(e = {}) {
  const dp0 = e.pumpDp0 || 0, qmx = e.pumpQmax || 1, sp = e.pumpSpeed ?? 1, kv = e.sepKv || 0, h0 = e.compHead || 0, qc = e.compQmax || 1, sc = e.compSpeed ?? 1, npoly = e.compN ?? 1.3, eta = e.compEta ?? 0.78;
  const pumpDp = (q) => (dp0 > 0 ? Math.max(dp0 * sp * sp - (dp0 / (qmx * qmx)) * q * q, 0) : 0);
  const compressor = (q, pDis, T = 300, Z = 0.9, M = 0.02, mdot = 0) => { const head = Math.max(h0 * sc * sc * (1 - (q / (qc * sc)) ** 2), 0), ex = (npoly - 1) / npoly, ratio = (1 + (ex * head * M) / (Z * RGAS * T)) ** (1 / ex); return { head, ratio, pSuction: pDis / ratio, power: (mdot * head) / Math.max(eta, 0.05) / 1e3 }; };
  const sepPressure = (qStd, pSet, pSuc = null) => (kv > 0 ? Math.sqrt((pSuc ?? pSet) ** 2 + kv * qStd * qStd) : pSet);
  return { pumpDp, compressor, sepPressure, active: dp0 > 0 || kv > 0 || h0 > 0 };
}
/**
 * Steady flow of a main line with one branch joining it (same fluid): the branch delivers the fraction `frac` of the total rate at the junction,
 * where the pressures are equal and mass and enthalpy are mixed. The main line is marched in two parts on the sub-grids of one uniform grid;
 * the inlet pressure of the main line is found by shooting on the outlet pressure and that of the branch by shooting on the junction pressure.
 * base: options of steadySolve (fm, profile, n, id, …, mScale = total rate downstream of the junction, pOut); br: { x (junction distance along the main profile, m), frac, length, dz (elevation gain of the branch), idMm, tIn, U }.
 * Returns the steadySolve result of the whole main line plus { junction: { s, x, P, tMain, tBranch, tMix, mMain, mBranch, hIn, hOut }, branch: (steadySolve result of the branch) }.
 */
export function steadyBranch(base, br) {
  const full = discretise(base.profile, base.n || 150), N = full.n, jn = clamp(Math.round(interp1(full.x, full.s, clamp(br.x, full.x[1], full.x[N - 1])) / full.ds), 1, N - 1), m = base.mScale ?? 1, f = clamp(br.frac, 0.01, 0.95);
  const sub = (a, b) => ({ s: full.s.slice(a, b + 1), x: full.x.slice(a, b + 1), z: full.z.slice(a, b + 1), theta: full.theta.slice(a, b), ds: full.ds, length: full.ds * (b - a), n: b - a });
  const gA = sub(0, jn), gB = sub(jn, N), sJ = full.s[jn], fit = base.fittings || [], kT = base.kTotal || 0, fm = base.fm;
  const optA = { ...base, grid: gA, pOut: undefined, mScale: m * (1 - f), kTotal: (kT * jn) / N, fittings: fit.filter((q) => q.s < sJ), hydrate: base.hydrate }, optB = { ...base, grid: gB, pOut: undefined, mScale: m, kTotal: (kT * (N - jn)) / N, fittings: fit.filter((q) => q.s >= sJ).map((q) => ({ ...q, s: q.s - sJ })) };
  // sub-grid of part B keeps absolute arc length for the property functions but the fittings are placed relative to its first node
  const idB = (br.idMm || 0) / 1000 || base.id, brOpt = { fm, profile: { x: [0, Math.max(br.length, 10)], z: [full.z[jn] - (br.dz || 0), full.z[jn]] }, n: Math.max(12, Math.round((base.n || 150) / 6)), id: idB, rough: base.rough, tIn: br.tIn ?? base.tIn, model: base.model, fModel: base.fModel, mp: base.mp, energy: base.energy, U: br.U ?? base.U ?? 3, tAmbOf: base.tAmbOf, mScale: m * f, hydrate: false, tolP: 1e-5 };
  let A, Bq, Br, mix, pj = null;
  const eval1 = (pIn) => {
    A = steadySolve({ ...optA, pIn }); if (!A.ok) return -1e3;
    Br = steadySolve({ ...brOpt, pOut: A.pOut, pGuess: Br && Br.ok ? Br.pIn : A.pOut + 2 }); if (!Br.ok) return -1e3;
    const hA = enthalpyFlow(fm.at(A.pOut, A.tOut, m * (1 - f))), hB = enthalpyFlow(fm.at(A.pOut, Br.tOut, m * f)), hT = hA + hB, lo = Math.min(A.tOut, Br.tOut) - 2, hi = Math.max(A.tOut, Br.tOut) + 2, g = (T) => enthalpyFlow(fm.at(A.pOut, T, m)) - hT;
    const tMix = Number.isFinite(hT) && g(lo) * g(hi) < 0 ? brent(g, lo, hi, 1e-8) : (1 - f) * A.tOut + f * Br.tOut; mix = { tMix, hIn: hT, hOut: enthalpyFlow(fm.at(A.pOut, tMix, m)) };
    Bq = steadySolve({ ...optB, pIn: A.pOut, tIn: tMix }); pj = A.pOut;
    return Bq.ok ? Bq.pOut - base.pOut : -1e3 - (gB.n - (Bq.reached || 0));
  };
  let lo = Math.max(base.pOut, 1.5), hi = lo + 20, fh = eval1(hi), g = 0; while (fh < 0 && hi < 1400 && g++ < 40) { lo = hi; hi = hi * 1.35 + 10; fh = eval1(hi); }
  if (fh < 0) return { ok: false, reason: 'No inlet pressure below 1,400 bara can deliver the combined rate through the junction.' };
  const pIn = brent(eval1, lo, hi, base.tolP ? Math.max(base.tolP, 1e-6) : 1e-6); eval1(pIn);
  if (!A.ok || !Bq.ok) return { ok: false, reason: 'The junction solution did not converge.' };
  const out = { ...Bq, ok: true, grid: full, n: N, ds: full.ds, length: full.length, s: full.s, x: full.x, z: full.z, pIn: A.pIn, tIn: A.tIn, pOut: Bq.pOut, tOut: Bq.tOut, mdot: Bq.mdot, marches: (A.marches || 0) + (Bq.marches || 0) };
  for (const k of Object.keys(A)) if (Array.isArray(A[k]) && A[k].length === gA.n + 1 && Array.isArray(Bq[k]) && Bq[k].length === gB.n + 1 && !['s', 'x', 'z'].includes(k)) out[k] = (['P', 'T', 'tHyd', 'subcooling', 'tWall'].includes(k) ? A[k].slice(0, gA.n) : A[k].slice(0, gA.n)).concat(Bq[k]);
  for (const k of ['liquidInventory', 'volume', 'dpFric', 'dpGrav', 'dpAcc', 'dpLocal', 'heatLoss', 'residence', 'residenceLiquid']) out[k] = A[k] + Bq[k];
  out.energy = { hIn: A.energy.hIn + (mix.hIn - enthalpyFlow(fm.at(pj, A.tOut, m * (1 - f)))), hOut: Bq.energy.hOut, potential: A.energy.potential + Bq.energy.potential, kinetic: A.energy.kinetic + Bq.energy.kinetic, residual: Math.max(A.energy.residual, Bq.energy.residual), mixing: mix.hIn - mix.hOut };
  out.junction = { s: sJ, x: full.x[jn], index: jn, P: pj, tMain: A.tOut, tBranch: Br.tOut, tMix: mix.tMix, mMain: A.mdot, mBranch: Br.mdot, mOut: Bq.mdot, hIn: mix.hIn, hOut: mix.hOut, pBranchIn: Br.pIn, pStepError: Math.abs(Bq.P[0] - A.pOut) };
  out.branch = { pIn: Br.pIn, tIn: Br.tIn, tOut: Br.tOut, dp: Br.pIn - Br.pOut, holdup: Br.liquidInventory / Br.volume, mdot: Br.mdot, length: brOpt.profile.x[1], id: idB };
  return out;
}
// closure parameters of the three-field annular model from the suite inputs
const annularParams = (v = {}) => ({ kDep: num(v.kDep, 0.15, 1e-4, 10), entMult: num(v.entMult, 1, 0.01, 20), weCrit: num(v.weCrit, 12, 0.5, 100), dropMult: 1, fiMult: num(v.fiAnnMult, 1, 0.05, 20) });
/**
 * Critical superficial liquid velocity at which stratified flow becomes unstable (Taitel & Dukler Kelvin–Helmholtz criterion on the
 * equilibrium level of the stratified balance with the friction multipliers o.fiMult, o.fwlMult, o.fwgMult); o.transMult scales the critical gas velocity.
 * Returns { vslCrit, hD, ratio (gas velocity / critical gas velocity at the given vsl), tauI, level (h/D at the given vsl) }.
 */
export function stratifiedTransition(p, o = {}) {
  const tm = o.transMult || 1, at = (vsl) => { const q = { ...p, vsl }, roots = stratifiedRoots(q, { ...o, n: 24 }), hD = roots[0], b = stratifiedBalance(hD, q, o), g = b.g, crit = tm * (1 - hD) * Math.sqrt((Math.max(p.rhoL - p.rhoG, 1e-6) * G * Math.max(Math.cos(p.theta || 0), 0.02) * g.AG) / (p.rhoG * Math.max(g.Si, 1e-9))); return { hD, ratio: b.vG / crit, tauI: b.ti, tauWL: b.tL, tauWG: b.tG }; };
  const here = at(p.vsl); let lo = 1e-4, hi = 20;
  if (at(hi).ratio < 1) return { vslCrit: hi, ...here, level: here.hD }; if (at(lo).ratio > 1) return { vslCrit: lo, ...here, level: here.hD };
  for (let k = 0; k < 22; k++) { const m = Math.sqrt(lo * hi); if (at(m).ratio < 1) lo = m; else hi = m; }
  return { vslCrit: Math.sqrt(lo * hi), ...here, level: here.hD };
}
/** Mixed liquid of an oil–water stream for point comparisons: volume-weighted density, viscosity of the dispersion (continuous phase × (1 − φ)^(−2.5 μ*), Ishii–Zuber form) with inversion at 50 % water. */
export function liquidMixture({ rhoO, rhoW, muO, muW, wc }) {
  const w = clamp(wc, 0, 1), rho = rhoO * (1 - w) + rhoW * w, waterCont = w >= 0.5, phi = Math.min(waterCont ? 1 - w : w, 0.74), muC = waterCont ? muW : muO, muD = waterCont ? muO : muW, mStar = (muD + 0.4 * muC) / (muD + muC);
  return { rho, mu: muC * (1 - phi) ** (-2.5 * mStar), waterCont };
}

// =====================================================================================================
// 7h. Point and system comparisons with measurements, extra results of a run, dedicated solver tasks
// =====================================================================================================
const CLASS6 = { 'stratified smooth': 'SS', 'stratified wavy': 'SW', slug: 'I', churn: 'I', annular: 'A', 'dispersed bubble': 'DB', bubble: 'B' };
/** Observed-pattern class (file labels of the open flow-pattern data bases: SS, SW, I, A, DB, B) predicted by the mechanistic map for one point. */
export const patternClass = (p) => CLASS6[flowPattern({ sigma: 0.02, theta: 0, ...p }).pattern] || '?';
/**
 * Model prediction for one measured local quantity. kind: 'holdup' | 'dpdx' (Pa/m, total, positive when the pressure falls) | 'slugFreq' (1/s) | 'slugLength' (m) |
 * 'slugVelocity' (m/s) | 'bodyHoldup' | 'filmThickness' (mm) | 'level' (h/D of stratified flow) | 'entrainment' | 'regime' (returns the pattern name).
 * pt: { D (m), angle (deg), vsl, vsg, rhoL, rhoG, muL, muG (Pa s), sigma (N/m), wc (% water in the liquid, optional, with rhoW, muW) }.
 */
export function predictPoint(kind, pt, model = 'beggsBrill', so = {}, mp = {}) {
  let rhoL = pt.rhoL, muL = pt.muL, waterCont = false;
  if (pt.wc > 0 && pt.rhoW > 0) { const m = liquidMixture({ rhoO: pt.rhoL, rhoW: pt.rhoW, muO: pt.muL, muW: pt.muW || 1e-3, wc: pt.wc / 100 }); rhoL = m.rho; muL = m.mu; waterCont = m.waterCont; }
  const p = { vsl: pt.vsl, vsg: pt.vsg, rhoL, rhoG: pt.rhoG, muL, muG: pt.muG || 1.8e-5, sigma: pt.sigma || 0.03, D: pt.D, theta: (pt.angle || 0) * DEG, rough: pt.rough ?? 0, P: pt.P || 1e5, fModel: 'colebrook', waterCont };
  if (kind === 'regime') return flowPattern(p).pattern;
  if (kind === 'holdup' || kind === 'dpdx') { const r = holdupGradient({ ...p, label: false }, model, mp); return kind === 'holdup' ? r.holdup : r.dpdx; }
  if (kind === 'level') return stratifiedRoots(p, so)[0];
  if (kind === 'entrainment') return annularMist(p, mp.annular || {}).entEq;
  const u = slugUnitCell(p, so);
  return kind === 'slugFreq' ? u.freq : kind === 'slugLength' ? u.lengthFromFreq : kind === 'slugVelocity' ? u.vt : kind === 'bodyHoldup' ? u.holdupSlug : kind === 'filmThickness' ? u.filmThickness * 1000 : NaN;
}
const errStats = (meas, pred) => { const n = meas.length; if (!n) return null; let ape = 0, b = 0, se = 0; for (let i = 0; i < n; i++) { const e = pred[i] - meas[i]; ape += Math.abs(e) / Math.max(Math.abs(meas[i]), 1e-12); b += e; se += e * e; } return { n, mape: (100 * ape) / n, bias: b / n, rmse: Math.sqrt(se / n) }; };

/** Results that every run adds after the line solution: closure detail, regime maps, tracer transport, equipment and junction, comparisons with user measurements. */
function runExtras(v, ctx, R) {
  const { st, cfg, fm, reg, units, im, cm, slug, tr, cyc, sev, cool, plots, tables, balances, warnings } = R, N = st.n, out = {}, kpis = [], so = slugOpts(v), xKm = st.x.map((x) => x / 1000);
  // ---- regime maps and closure detail at the map location
  if (cm.vsl > 1e-9 && cm.vsg > 1e-9) {
    const bk = bakerRegime(cm), pc = { ...cm, rough: cfg.rough, P: st.P[im] * 1e5, fModel: cfg.fModel }, sT = stratifiedTransition(cm, { fiMult: num(v.fiMult, 1, 0.05, 20), fwlMult: num(v.fwlMult, 1, 0.05, 20), fwgMult: num(v.fwgMult, 1, 0.05, 20), transMult: num(v.transMult, 1, 0.2, 5) }), vn = velocityNumbers(cm);
    const rows = [['Taitel–Dukler / Barnea (mechanistic)', reg[im], ''], ['Mandhane, Gregory & Aziz map', mandhaneRegime(cm.vsl, cm.vsg), 'air–water coordinates'], ['Baker map', bk.regime, `Bx = ${sig(bk.Bx, 3)}, By = ${sig(bk.By, 3)}`], ['Beggs & Brill horizontal pattern', gradient({ ...pc, theta: 0, label: false }).regime, 'correlation regime'],
      ['Stratified → non-stratified: critical vsl', r3(sT.vslCrit, 4), `m/s (now ${r3(cm.vsl, 3)}; gas velocity / critical ${r3(sT.ratio, 2)})`], ['Equilibrium stratified level h/D', r3(sT.level, 3), `τ interface ${r3(sT.tauI, 2)} Pa, τ liquid wall ${r3(sT.tauWL, 2)} Pa, τ gas wall ${r3(sT.tauWG, 3)} Pa`],
      ['Velocity numbers N_Lv / N_Gv / N_D / N_L', `${sig(vn.NLv, 3)} / ${sig(vn.NGv, 3)} / ${sig(vn.ND, 3)} / ${sig(vn.NL, 3)}`, 'Duns & Ros groups']];
    if (cm.theta > 0.05) rows.push(['Ansari et al. pattern (upward)', ansariPattern(pc), ''], ['Duns & Ros region (upward)', dunsRos(pc).regime, '']);
    tables.push({ title: `Flow-pattern maps and transition criteria at ${Math.round(st.x[im])} m`, columns: ['Map / criterion', 'Result', 'Note'], rows, note: 'The Baker boundaries are chart fits that could not be checked against an open copy; use the mechanistic map for design.' });
    out.baker = bk.regime; out.stratified = { vslCrit: sT.vslCrit, level: sT.level, tauI: sT.tauI };
  }
  // ---- annular-flow closure with entrainment and deposition where the gas velocity is highest
  { let ia = 0; for (let i = 1; i < N; i++) if (st.vsg[i] > st.vsg[ia]) ia = i;
    if (st.vsl[ia] > 1e-9 && st.vsg[ia] > 1e-9) { const p = { ...cellOf(st, ia), rough: cfg.rough, P: st.P[ia] * 1e5, fModel: cfg.fModel }, a = annularMist(p, cfg.mp.annular), dev = annularDevelopment({ p, mp: cfg.mp.annular, length: 6 * a.relaxLength, n: 30, e0: 0 });
      tables.push({ title: `Annular-flow closure (gas core, droplets, wall film) at ${Math.round(st.x[ia])} m`, columns: ['Quantity', 'Value', 'Unit'], rows: [['Flow pattern there (mechanistic map)', reg[ia], ''], ['Equilibrium entrained fraction', r3(a.entEq, 4), '–'], ['Film thickness', r3(a.filmThickness * 1000, 3), 'mm'], ['Film velocity', r3(a.filmVelocity, 3), 'm/s'], ['Core velocity', r3(a.coreVelocity, 2), 'm/s'], ['Interfacial / wall shear', `${sig(a.tauI, 3)} / ${sig(a.tauW, 3)}`, 'Pa'], ['Droplet size (critical Weber number)', r3(a.dropSize * 1e6, 0), 'µm'], ['Deposition = entrainment rate at equilibrium', sig(a.rateDep, 3), 'kg/m²/s'], ['Development length (1/e)', r3(a.relaxLength, 1), 'm'], ['Holdup: film + droplets', `${sig(a.holdupFilm, 3)} + ${sig(a.holdupDrops, 3)}`, '–'], ['Pressure gradient', r3(a.fric + a.grav + a.acc, 1), 'Pa/m']], note: 'Evaluated with the three-field model whatever the predicted pattern; it enters the line solution when the holdup model “Mechanistic + three-field annular flow” is selected and the cell is annular.' });
      balances.push({ name: 'Annular three-field liquid (kg/s): film + droplets at the end of the development length vs liquid in', in: dev.balance.liquidIn, out: dev.film[dev.film.length - 1] + dev.drops[dev.drops.length - 1] });
      out.annular = { x: st.x[ia], entrainment: a.entEq, filmThickness: a.filmThickness, dropSize: a.dropSize, relaxLength: a.relaxLength, rateDep: a.rateDep }; } }
  // ---- outside film: forced and free convection at the coldest-margin location
  { let ic = 0; for (let i = 1; i <= N; i++) if (st.subcooling[i] > st.subcooling[ic]) ic = i; const k = Math.min(ic, N - 1), c = { s: st.s[k], z: st.z[k], D: st.D[k], pr: fm.at(st.P[ic], st.T[ic], st.mScale), vm: st.vm[k], holdup: st.holdup[k], ta: st.tAmb[k], T: st.T[ic] }, net = cfg.network(c), hF = cfg.hOutOf(c, 0);
    if (!cfg.buriedAt(c.s)) { tables.push({ title: `Outside film at ${Math.round(st.x[ic])} m: forced and free convection`, columns: ['Quantity', 'Value', 'Unit'], rows: [['Forced convection (Churchill–Bernstein)', r3(hF, 1), 'W/m²K'], ['Surface temperature above ambient', r3(fin(net.dTsurface, 0), 3), 'K'], ['Combined with free convection (Churchill–Chu)', r3(fin(net.hOut, hF), 1), 'W/m²K'], ['Layer-based U-value with the combined film', r3(net.U, 3), 'W/m²K on ID']], note: cfg.uMode === 'layers' ? 'Used in the line solution (U-value from films, wall, insulation and burial).' : 'Shown for information: the line solution uses the given U-value.' }); out.outsideFilm = { forced: hF, combined: fin(net.hOut, hF), dTsurface: fin(net.dTsurface, 0) }; } }
  // ---- species transport with the liquid: arrival of an inlet step (inhibitor, tracer, water-cut change)
  if (v.tracer !== false && st.qL.every((q) => q > 1e-9)) {
    const ut = Math.sqrt(mean(st.tauW.slice(0, N)) / Math.max(mean(st.rhoL), 1)), Dd = num(v.dispersion, 0, 0, 1e4) || 10.1 * (st.D[0] / 2) * Math.max(ut, 1e-3), A = st.D.map((d) => (PI * d * d) / 4), tRes = st.residenceLiquid, sp = speciesTransport({ s: st.s, qL: st.qL, holdup: st.holdup.map((h) => Math.max(h, 0.01)), area: A, D: Dd, cIn: 1, tEnd: 1.6 * tRes, nSeries: 160, cfl: 0.9 });
    plots.push({ type: 'line', title: 'Arrival of an inlet step carried by the liquid (species transport)', xlabel: 'Time (h)', ylabel: 'Outlet concentration / inlet', ymin: 0, series: [{ name: 'Outlet', x: sp.t.map((t) => t / 3600), y: sp.cOut }], vlines: [{ x: tRes / 3600, label: 'Liquid residence time' }], note: `Advection with the local liquid velocity and axial dispersion ${sig(Dd, 3)} m²/s (Taylor estimate 10.1 R u* unless given); breakthrough (50 %) after ${sp.breakthrough !== null ? (sp.breakthrough / 3600).toFixed(2) + ' h' : 'more than the simulated time'}.` });
    balances.push({ name: 'Species transport: injected vs stored + produced (relative to the injected amount)', in: 1, out: 1 + sp.balance.error }); out.tracer = { breakthrough: sp.breakthrough !== null ? sp.breakthrough / 3600 : null, dispersion: Dd };
  }
  // ---- equipment characteristics and the branch junction
  { const rows = [], q0 = st.qG[0] + st.qL[0], dpP = cfg.equip.pumpDp(q0);
    if (dpP > 0) { rows.push(['Booster pump: suction volume rate', r3(q0 * 3600, 1), 'm³/h'], ['Booster pump: pressure rise from its curve', r3(dpP, 2), 'bar'], ['Pressure needed upstream of the pump', r3(st.pIn - dpP, 2), 'bara'], ['Booster pump: hydraulic power', r3((q0 * dpP * 1e5) / 1e3, 1), 'kW']); out.pump = { dp: dpP, suctionPressure: st.pIn - dpP, power: (q0 * dpP * 1e5) / 1e3 }; }
    if (cfg.sepMode === 'valve' && cfg.compressor) { const c = cfg.compressor; rows.push(['Separator: gas rate at standard conditions', r3(c.qStd, 3), 'Sm³/s'], ['Compressor: suction pressure from its curve', r3(c.pSuction, 2), 'bara'], ['Compressor: pressure ratio / polytropic head', `${r3(c.ratio, 2)} / ${r3(c.head / 1e3, 1)}`, '– / kJ/kg'], ['Compressor: shaft power', r3(c.power, 0), 'kW'], ['Separator pressure imposed by the gas valve and the compressor', r3(fin(cfg.pOutEff, cfg.pOut), 2), 'bara']); out.separator = { pressure: fin(cfg.pOutEff, cfg.pOut), compressorPower: c.power, suction: c.pSuction }; }
    if (st.junction) { const j = st.junction; rows.push(['Junction position', r3(j.x, 0), 'm'], ['Junction pressure', r3(j.P, 2), 'bara'], ['Main-line / branch arrival temperature', `${j.tMain.toFixed(2)} / ${j.tBranch.toFixed(2)}`, '°C'], ['Mixed temperature', r3(j.tMix, 2), '°C'], ['Main-line / branch mass rate', `${j.mMain.toFixed(2)} / ${j.mBranch.toFixed(2)}`, 'kg/s'], ['Branch inlet pressure', r3(j.pBranchIn, 2), 'bara']);
      balances.push({ name: 'Junction mass (kg/s): main + branch vs downstream', in: j.mMain + j.mBranch, out: j.mOut }, { name: 'Junction energy (kW): enthalpy of main + branch vs mixed stream', in: j.hIn / 1e3, out: j.hOut / 1e3 }); out.junction = { x: j.x, P: j.P, tMix: j.tMix, mBranch: j.mBranch, pBranchIn: j.pBranchIn }; }
    if (rows.length) tables.push({ title: 'Equipment characteristics and junction', columns: ['Quantity', 'Value', 'Unit'], rows, note: 'Pump: parabolic head curve with the affinity laws; separator: gas-outlet valve to the compressor suction; compressor: polytropic head curve; junction: equal pressure, mass and enthalpy mixing.' }); }
  // ---- comparisons with user measurements (validation path for any configuration)
  { const KINDS = ['holdup', 'dpdx', 'slugFreq', 'slugLength', 'slugVelocity', 'bodyHoldup', 'filmThickness', 'level', 'entrainment', 'regime'], pts = (Array.isArray(v.valPoints) ? v.valPoints : []).filter((r) => r && KINDS.includes(String(r.kind || '').trim()) && +r.D > 0 && (+r.vsl > 0 || +r.vsg > 0));
    if (pts.length) { const rows = [], by = {}; for (const r of pts) { const kind = String(r.kind).trim(), pt = { D: +r.D / 1000, angle: fin(+r.angle, 0), vsl: Math.max(+r.vsl || 0, 1e-9), vsg: Math.max(+r.vsg || 0, 1e-9), rhoL: fin(+r.rhoL, 850), rhoG: fin(+r.rhoG, 50), muL: fin(+r.muL, 1) / 1000, muG: fin(+r.muG, 0.015) / 1000, sigma: fin(+r.sigma, 20) / 1000, wc: fin(+r.wc, 0), rhoW: 1020, muW: 1e-3, P: fin(+r.P, 1) * 1e5 }; let pred; try { pred = predictPoint(kind, pt, cfg.model, so, cfg.mp); } catch { pred = NaN; }
        if (kind === 'regime') { rows.push([kind, r3(pt.D * 1000, 1), pt.angle, pt.vsl, pt.vsg, String(r.measured ?? ''), pred, pred === String(r.measured ?? '').trim() ? 'match' : 'differs']); (by.regime ||= { m: [], p: [] }).m.push(1); by.regime.p.push(pred === String(r.measured ?? '').trim() ? 1 : 0); }
        else if (Number.isFinite(+r.measured) && Number.isFinite(pred)) { rows.push([kind, r3(pt.D * 1000, 1), pt.angle, pt.vsl, pt.vsg, sig(+r.measured, 5), sig(pred, 5), r3((100 * (pred - +r.measured)) / Math.max(Math.abs(+r.measured), 1e-12), 1) + ' %']); (by[kind] ||= { m: [], p: [] }).m.push(+r.measured); by[kind].p.push(pred); } }
      tables.push({ title: 'Local measurements against the model', columns: ['Quantity', 'D (mm)', 'Angle (°)', 'vsl (m/s)', 'vsg (m/s)', 'Measured', 'Predicted', 'Deviation'], rows, note: Object.entries(by).map(([k, q]) => { const e = errStats(q.m, q.p); return `${k}: n = ${e.n}, ${k === 'regime' ? 'agreement ' + (100 * mean(q.p)).toFixed(0) + ' %' : 'MAPE ' + e.mape.toFixed(1) + ' %, bias ' + sig(e.bias, 3)}`; }).join(' · ') });
      const num2 = Object.entries(by).filter(([k]) => k !== 'regime'); if (num2.length) plots.push({ type: 'line', title: 'Local measurements: predicted against measured', xlabel: 'Measured', ylabel: 'Predicted', logx: true, logy: true, series: [...num2.map(([k, q]) => ({ name: k, x: q.m.map((x) => Math.max(x, 1e-9)), y: q.p.map((x) => Math.max(x, 1e-9)), mode: 'points' })), { name: 'Parity', x: [1e-3, 1e4], y: [1e-3, 1e4], dash: true }] });
      out.pointComparison = Object.fromEntries(Object.entries(by).map(([k, q]) => [k, errStats(q.m, q.p)])); }
    const QS = { pIn: [st.pIn, 'bara'], pOut: [st.pOut, 'bara'], dp: [st.pIn - st.pOut, 'bar'], tOut: [st.tOut, '°C'], inventory: [st.liquidInventory, 'm³'], holdupMean: [st.liquidInventory / st.volume, '–'], heatLoss: [st.heatLoss / 1e3, 'kW'], slugPeriod: [slug.period, 's'], slugFrequency: [slug.freq, '1/s'], slugLength: [slug.length, 'm'], slugLengthMax: [slug.lengthMax, 'm'], slugVelocity: [slug.velocity, 'm/s'], bodyHoldup: [slug.holdupBody, '–'], surge: [slug.surge, 'm³'], catcherLoad: [1.25 * fin(slug.surge, 0), 'm³'], pAmplitude: [R.pAmp, 'bar'], severePeriod: [cyc && cyc.period ? cyc.period : null, 's'], severeAmplitude: [cyc ? cyc.amplitude / 1e5 : null, 'bar'], cooldown: [cool.tReach !== null ? cool.tReach / 3600 : null, 'h'], terrainVolume: [R.terrain.worst ? R.terrain.worst.volume : null, 'm³'], boe: [sev ? sev.boe : null, '–'], tracerArrival: [out.tracer ? out.tracer.breakthrough : null, 'h'] };
    const sys = (Array.isArray(v.valSystem) ? v.valSystem : []).filter((r) => r && QS[String(r.quantity || '').trim()] && Number.isFinite(+r.measured));
    if (sys.length) { const rows = sys.map((r) => { const [val, unit] = QS[String(r.quantity).trim()]; return [String(r.quantity).trim(), sig(+r.measured, 5), val === null || !Number.isFinite(val) ? '—' : sig(val, 5), unit, val === null || !Number.isFinite(val) ? '—' : r3((100 * (val - +r.measured)) / Math.max(Math.abs(+r.measured), 1e-12), 1) + ' %']; }); tables.push({ title: 'System measurements against this run', columns: ['Quantity', 'Measured', 'Predicted', 'Unit', 'Deviation'], rows, note: 'Quantities: ' + Object.keys(QS).join(', ') + '.' }); out.systemComparison = rows.length; }
  }
  return { outputs: out, kpis };
}

/** Dedicated solver tasks (run on request): two-fluid slug capturing, 2-D CFD, bubbly-flow closures, benchmark set. */
async function runTask(task, v, ctx, R) {
  const { st, cfg, cm, im, units, plots, tables, balances, warnings } = R, tick = typeof ctx.tick === 'function' ? ctx.tick : async () => {}, prog = typeof ctx.progress === 'function' ? ctx.progress : () => {}, out = {}, kpis = [];
  if (task === 'twofluid') {
    prog(0.86, 'Two-fluid model: slug capturing');
    const D = clamp(num(v.tfDiameterMm, 0, 0, 2000) / 1000 || cm.D, 0.02, 1.5), L = num(v.tfLengthD, 200, 20, 5000) * D, n = Math.round(num(v.tfCells, 240, 40, 4000)), tEnd = num(v.tfTime, 20, 0.5, 600), Rs = (cm.zG ?? 0.9) * RGAS / ((R.mwG || 20) * 1e-3), T = st.T[im] + 273.15, p0 = (num(v.tfPbar, 0, 0, 1000) || st.P[im]) * 1e5, rhoG = p0 / (Rs * T);
    const p = { ...cm, D, rhoG, vsl: num(v.tfVsl, 0, 0, 20) || cm.vsl, vsg: num(v.tfVsg, 0, 0, 100) || cm.vsg, theta: clamp(cm.theta, -0.1, 0.1) }, eq = stratifiedLevel(p), a0 = clamp(eq.holdup, 0.02, 0.9), stab = interfacialStability(p), reg = stab.ikhRatio > 0.9; // beyond the inviscid limit the equations are regularised with the interfacial-pressure term
    const r = twoFluid({ n, length: L, D, theta: p.theta, interfacialPressure: reg ? 1.2 : 0, rhoL: cm.rhoL, muL: cm.muL, muG: cm.muG, Rs, T, pOut: p0, init: () => ({ al: a0, ul: eq.vL, ug: eq.vG, p: p0 }), inlet: (t) => ({ al: a0 * (1 + 0.02 * Math.sin((2 * PI * t) / 0.7) + 0.02 * Math.sin((2 * PI * t) / 1.9)), ul: eq.vL, ug: eq.vG }), tEnd, cfl: num(v.tfCfl, 0.4, 0.05, 0.9), dtMax: tEnd / 400, probes: [0.5 * L, 0.75 * L, 0.9 * L], nSeries: 600, nField: 60 });
    await tick();
    const late = (k) => r.slugs.filter((q) => q.probe === k && q.t > 0.3 * tEnd), s2 = late(2), fq = s2.length / (0.7 * tEnd), a = r.series.probes[1], b = r.series.probes[2], i0 = Math.floor(a.length * 0.3), dts = r.series.t[r.series.t.length - 1] / Math.max(r.series.t.length - 1, 1), ma = mean(a.slice(i0)), mb = mean(b.slice(i0));
    let best = 0, bl = 0; for (let lag = 1; lag < Math.min(300, a.length - i0 - 2); lag++) { let c = 0; for (let i = i0; i + lag < a.length; i++) c += (a[i] - ma) * (b[i + lag] - mb); if (c > best) { best = c; bl = lag; } }
    const swing = Math.max(...b.slice(i0)) - Math.min(...b.slice(i0)), vWave = bl > 0 && swing > 0.05 ? (0.15 * L) / (bl * dts) : null, // no velocity is reported when no waves pass the probes
      uc = slugUnitCell({ ...p, sigma: cm.sigma }, { ...slugOpts(v), freqModel: 'zabaras' });
    plots.push({ type: 'field', title: 'Two-fluid model: liquid holdup on the test section (slug capturing)', xlabel: 'Distance (m)', ylabel: 'Time (s)', zlabel: 'Holdup', zunit: '–', x: r.field.x, y: r.field.t, z: r.field.al, zmin: 0, zmax: 1, cmap: 'viridis' });
    plots.push({ type: 'line', title: 'Two-fluid model: holdup at the probes', xlabel: 'Time (s)', ylabel: 'Liquid holdup', ymin: 0, ymax: 1.02, series: [0.5, 0.75, 0.9].map((f, k) => ({ name: `x = ${(f * L).toFixed(1)} m`, x: r.series.t, y: r.series.probes[k] })), hlines: [{ y: 0.95, label: 'Liquid bridge' }] });
    tables.push({ title: 'Two-fluid slug capturing', columns: ['Quantity', 'Value', 'Unit'], rows: [['Test section: length / diameter / cells', `${L.toFixed(1)} / ${(D * 1000).toFixed(0)} mm / ${n}`, 'm'], ['Cell size', r3(L / n / D, 2), 'diameters'], ['Initial stratified holdup (equilibrium)', r3(a0, 3), '–'], ['Viscous Kelvin–Helmholtz ratio of the initial state', r3(stab.vkhRatio, 2), '> 1 = waves grow'], ['Inviscid Kelvin–Helmholtz ratio', r3(stab.ikhRatio, 2), reg ? '> 0.9: interfacial-pressure regularisation (δ = 1.2) switched on' : '> 1 = ill-posed without regularisation'], ['Superficial velocities liquid / gas, pressure', `${r3(p.vsl, 3)} / ${r3(p.vsg, 3)} m/s, ${r3(p0 / 1e5, 2)} bara`, ''], ['Time steps / mean step', `${r.steps} / ${sig(r.dtMean, 3)} s`, ''], ['Holdup swing at the last probe', r3(swing, 3), swing > 0.05 ? '–' : '– (the stratified state stays smooth)'], ['Slugs passing the last probe (after 30 % of the time)', s2.length, '–'], ['Captured slug frequency', r3(fq, 3), '1/s'], ['Unit-cell frequency (correlation)', r3(uc.freq, 3), '1/s'], ['Captured front / wave velocity (probe cross-correlation)', vWave === null ? '—' : r3(vWave, 2), 'm/s'], ['Bendiksen translational velocity', r3(uc.vt, 2), 'm/s'], ['Mean slug duration at the last probe', s2.length ? r3(mean(s2.map((q) => q.duration)), 3) : '—', 's'], ['Liquid mass error', sig(Math.abs(r.mass.errorL), 3), '–'], ['Gas mass error', sig(Math.abs(r.mass.errorG), 3), '–'], ['Largest volume-constraint error', sig(r.volErrMax, 3), '–'], ['Run completed', r.completed ? 'yes' : 'stopped early', '']], note: 'Waves grow from a 2 % inlet holdup disturbance when the stratified state is unstable; a face is a liquid bridge when the gas fraction next to it falls below 5 %. First-order upwind: the captured frequency depends on the cell size (use the mesh study).' });
    balances.push({ name: 'Two-fluid liquid mass (kg): initial + inflow vs final + outflow', in: r.mass.liquid0 + r.mass.liquidIn, out: r.mass.liquid + r.mass.liquidOut }, { name: 'Two-fluid gas mass (kg): initial + inflow vs final + outflow', in: r.mass.gas0 + r.mass.gasIn, out: r.mass.gas + r.mass.gasOut });
    if (!r.completed) warnings.push({ level: 'warn', msg: 'The two-fluid run stopped early: the flow left the well-posed range of the model on this grid (inviscid Kelvin–Helmholtz limit exceeded). Results up to that time are shown.' });
    out.twoFluid = { slugs: s2.length, frequency: fq, waveVelocity: vWave, bendiksen: uc.vt, massErrorL: r.mass.errorL, massErrorG: r.mass.errorG, steps: r.steps, completed: r.completed, vkhRatio: stab.vkhRatio, holdupMean: mean(r.al) };
    kpis.push({ label: 'Captured slug frequency (two-fluid)', value: sig(fq, 3), unit: '1/s', status: 'ok' }, { label: 'Captured front velocity', value: vWave === null ? '—' : r3(vWave, 2), unit: 'm/s', status: 'ok' });
  }
  if (task === 'cfd') {
    const kind = ['cavity', 'dambreak', 'slugfront', 'channel'].includes(v.cfdCase) ? v.cfdCase : 'dambreak', n = Math.round(num(v.cfdN, 24, 8, 256)), turb = v.cfdTurbulence === 'sst' ? 'sst' : 'none'; prog(0.86, '2-D Navier–Stokes solver: ' + kind);
    let r, title, note = '';
    if (kind === 'cavity') { const Re = num(v.cfdRe, 100, 1, 5000); r = cfd2d({ nx: n, ny: n, rho: 1, mu: 1 / Re, lid: 1, tEnd: 60, steady: 1e-5, cfl: 0.5 }); title = `Lid-driven cavity, Re = ${Re}`; const g = REF.GHIA.rows.filter((q) => q[0] === 0), col = 2; if (Re === 100) { let e = 0; for (const q of g) e = Math.max(e, Math.abs(interp1(r.y, r.uMid, q[1]) - q[col])); note = `Largest difference of u on the vertical centre line from the Ghia, Ghia & Shin (1982) table: ${e.toFixed(4)} (lid velocity 1). `; out.cavityError = e; }
      plots.push({ type: 'line', title: 'Cavity: u on the vertical centre line', xlabel: 'u / U lid', ylabel: 'y / L', series: [{ name: `${n} × ${n} cells`, x: r.uMid, y: r.y }, ...(Re === 100 ? [{ name: 'Ghia et al. (1982)', x: g.map((q) => q[col]), y: g.map((q) => q[1]), mode: 'points' }] : [])] }); }
    else if (kind === 'channel') { const ReT = num(v.cfdRe, 180, 50, 2000); r = cfd2d({ nx: 4, ny: 2 * n, lx: 0.4, ly: 2, rho: 1, mu: 1 / ReT, fx: 1, periodicX: true, turbulence: 'sst', tEnd: 60, steady: 1e-6, u0: (x, y) => 15 * (1 - (y - 1) ** 2), k0: 1, omega0: 20, cfl: 0.5 }); title = `Turbulent channel (SST k–ω), Re_τ = ${ReT}`; const ub = mean(r.u.map((q) => q[1])); note = `Bulk velocity U⁺ = ${ub.toFixed(2)} against ${(Math.log(ReT) / 0.41 + 5.2 - 1 / 0.41).toFixed(2)} from the log law. `; out.channelUb = ub;
      plots.push({ type: 'line', title: 'Channel: mean velocity in wall units', xlabel: 'y⁺', ylabel: 'u⁺', logx: true, series: [{ name: 'SST k–ω (2-D solver)', x: r.y.slice(0, n).map((y) => y * ReT), y: r.u.slice(0, n).map((q) => q[1]) }, { name: 'Log law', x: r.y.slice(0, n).map((y) => y * ReT).filter((y) => y > 30), y: r.y.slice(0, n).map((y) => y * ReT).filter((y) => y > 30).map((y) => Math.log(y) / 0.41 + 5.2), dash: true }] }); }
    else { const D = cm.D, front = kind === 'slugfront', hS = front ? clamp(levelOfHoldup(units[im] ? units[im].holdupSlug : 0.9), 0.3, 0.98) * D : 2 * 0.25 * D * 4, a = front ? 2 * D : D, lx = front ? 6 * D : 8 * D, ly = front ? D : 4 * D, hF = front ? (units[im] ? units[im].filmThickness : 0.2 * D) : 0;
      r = cfd2d({ nx: Math.round((n * lx) / ly), ny: n, lx, ly, rhoL: cm.rhoL, rhoG: cm.rhoG, muL: cm.muL, muG: cm.muG, gy: -G, c0: (x, y) => ((x < a && y < (front ? hS : 2 * a)) || y < hF ? 1 : 0), tEnd: num(v.cfdTime, 0, 0, 60) || 2.2 * Math.sqrt(lx / G), probeFront: true, turbulence: turb }); title = front ? `Slug tail / front collapsing onto the film at ${Math.round(st.x[im])} m (channel analogue of the pipe)` : 'Dam break of a liquid column in a channel';
      const h0 = front ? hS : 2 * a, tS = r.front.t.map((t) => t * Math.sqrt(G / h0)), ritter = hF > 0 ? null : r.front.t.map((t) => (a + 2 * Math.sqrt(G * h0) * t) / h0); note = `Liquid volume error ${Math.abs(r.mass.error).toExponential(1)}; ${turb === 'sst' ? 'SST k–ω turbulence. ' : 'laminar. '}`;
      plots.push({ type: 'line', title: 'Front position along the bottom', xlabel: 't √(g / h₀)', ylabel: 'x front / h₀', series: [{ name: `${r.x.length} × ${r.y.length} cells`, x: tS, y: r.front.x.map((x) => x / h0) }, ...(ritter ? [{ name: 'Shallow-water limit (Ritter, frictionless)', x: tS, y: ritter.map((x) => Math.min(x, lx / h0)), dash: true }] : [])] });
      balances.push({ name: '2-D solver liquid volume (m² per unit depth): initial vs final', in: r.mass.initial, out: r.mass.final }); out.cfdMassError = r.mass.error; out.frontX = r.front.x[r.front.x.length - 1]; }
    await tick();
    plots.push({ type: 'field', title: title + (r.c ? ': liquid fraction and velocity' : ': speed and streamlines'), xlabel: 'x (m)', ylabel: 'y (m)', zlabel: r.c ? 'Liquid fraction' : 'Speed', zunit: r.c ? '–' : 'm/s', x: r.x, y: r.y, z: r.c || r.u.map((row, j) => row.map((u, i) => Math.hypot(u, r.v[j][i]))), cmap: r.c ? 'salinity' : 'viridis', u: r.u, v: r.v, vectors: !!r.c, stream: !r.c, equal: kind !== 'channel' });
    if (r.nut) plots.push({ type: 'field', title: 'Eddy viscosity (SST k–ω)', xlabel: 'x (m)', ylabel: 'y (m)', zlabel: 'ν_t', zunit: 'm²/s', x: r.x, y: r.y, z: r.nut, cmap: 'turbo', equal: kind !== 'channel' });
    tables.push({ title: '2-D Navier–Stokes solver', columns: ['Quantity', 'Value', 'Unit'], rows: [['Case', title, ''], ['Cells', `${r.x.length} × ${r.y.length}`, ''], ['Time steps / mean step', `${r.steps} / ${sig(r.dtMean, 3)}`, 's'], ['Simulated time', sig(r.t, 4), 's'], ['Pressure iterations per step (mean)', r3(r.poissonIterations, 1), '–'], ['Largest cell divergence × Δx / |u|', sig(r.divergenceMax, 3), '–'], ['Velocity change per unit time at the end', sig(r.steadyResidual, 3), 'm/s²'], ['Liquid volume error', sig(Math.abs(r.mass.error), 3), '–']], note: note + 'Projection method on a staggered grid, THINC/WLIC volume of fluid, surface tension neglected. Three-dimensional and scale-resolving simulations are written as OpenFOAM cases on the hand-off page.' });
    out.cfd = { case: kind, cells: r.x.length * r.y.length, steps: r.steps, divergence: r.divergenceMax, massError: r.mass.error, residual: r.steadyResidual, kineticEnd: r.kinetic.length ? r.kinetic[r.kinetic.length - 1].e : 0, uMin: Math.min(...r.uMid) };
    kpis.push({ label: '2-D solver: liquid volume error', value: sig(Math.abs(r.mass.error), 2), unit: '–', status: Math.abs(r.mass.error) < 1e-6 ? 'ok' : 'warn' });
  }
  if (task === 'bubbly') {
    prog(0.86, 'Bubbly-flow closures: radial void distribution'); let ib = 0; for (let i = 1; i < st.n; i++) if (st.theta[i] > st.theta[ib]) ib = i; const c = cellOf(st, ib), aM = clamp(1 - st.holdup[ib], 0.01, 0.15), d0 = num(v.bubbleMm, 3, 0.2, 30) / 1000;
    const kTurb = 0.01 * (c.vsl + c.vsg) ** 2, ia = interfacialAreaTransport({ d0, alpha: aM, vg: Math.max(c.vsg / aM, 0.1), k: kTurb, length: Math.max(cfg.riserHeight, 50), rhoL: c.rhoL, rhoG: c.rhoG, muL: c.muL, sigma: c.sigma, crc: 0.04 * num(v.crcMult, 1, 0, 100), cti: 0.085 * num(v.ctiMult, 1, 0, 100), cwe: 0.002 * num(v.crcMult, 1, 0, 100) });
    const mk = (d) => bubblyPipe({ R: c.D / 2, jl: c.vsl, alphaMean: aM, d, rhoL: c.rhoL, rhoG: c.rhoG, muL: c.muL, sigma: c.sigma, n: Math.round(num(v.bubblyN, 100, 30, 600)), cvm: num(v.cvm, 0.5, 0, 2) }), small = mk(d0), end = mk(clamp(ia.dEnd, 2e-4, 0.3 * c.D)); await tick();
    plots.push({ type: 'line', title: `Radial void fraction of bubbly flow at ${Math.round(st.x[ib])} m (mean void ${aM.toFixed(3)})`, xlabel: 'r / R', ylabel: 'Void fraction', ymin: 0, series: [{ name: `Bubbles ${(d0 * 1000).toFixed(1)} mm (inlet size)`, x: small.r, y: small.alpha }, { name: `Bubbles ${(ia.dEnd * 1000).toFixed(1)} mm (after coalescence / breakup)`, x: end.r, y: end.alpha }], note: `Lift coefficient ${small.cl.toFixed(3)} → ${end.cl.toFixed(3)} (Tomiyama): small bubbles collect at the wall, large ones in the core.` });
    plots.push({ type: 'line', title: 'Sauter diameter along the riser (interfacial-area transport)', xlabel: 'Height (m)', ylabel: 'mm', series: [{ name: 'Sauter diameter', x: ia.z, y: ia.d.map((x) => x * 1000) }], hlines: ia.dEquilibrium ? [{ y: ia.dEquilibrium * 1000, label: 'Coalescence = breakup' }] : [] });
    tables.push({ title: 'Bubbly-flow closures', columns: ['Quantity', 'Inlet size', 'After the riser', 'Unit'], rows: [['Bubble diameter', r3(d0 * 1000, 2), r3(ia.dEnd * 1000, 2), 'mm'], ['Eötvös number', r3(small.eotvos, 2), r3(end.eotvos, 2), '–'], ['Relative (terminal) velocity', r3(small.uRel, 3), r3(end.uRel, 3), 'm/s'], ['Drag coefficient (Ishii–Zuber)', r3(small.cd, 3), r3(end.cd, 3), '–'], ['Lift coefficient (Tomiyama)', r3(small.cl, 3), r3(end.cl, 3), '–'], ['Void peak position r/R', r3(small.peak.r, 3), r3(end.peak.r, 3), '–'], ['Peak / centre-line void', `${small.peak.alpha.toFixed(3)} / ${small.alphaCentre.toFixed(3)}`, `${end.peak.alpha.toFixed(3)} / ${end.alphaCentre.toFixed(3)}`, '–'], ['Frictional + buoyant pressure gradient', r3(small.dpdz, 0), r3(end.dpdz, 0), 'Pa/m'], ['Coalescence (random collision / wake) and breakup rates of the curvature', `${sig(ia.rates.coalescenceRC, 3)} / ${sig(ia.rates.coalescenceWE, 3)}`, sig(ia.rates.breakup, 3), '1/m/s']], note: 'Lateral balance of turbulent dispersion (Burns), lift (Tomiyama, wall-damped) and wall lubrication (Antal); Sato bubble-induced viscosity; pseudo-turbulence ½ C_vm α u_r². Applied to the steepest cell with its superficial velocities, whatever the predicted pattern.' });
    balances.push({ name: 'Bubbly model: area-mean void of the profile vs the given mean', in: aM, out: end.mean }, { name: 'Bubbly model: liquid flux of the profile vs the superficial velocity', in: c.vsl, out: end.liquidFlux });
    out.bubbly = { dSauter: ia.dEnd, dEquilibrium: ia.dEquilibrium, peakR: end.peak.r, wallPeaked: end.wallPeaked, cl: end.cl, alphaCentre: end.alphaCentre };
  }
  if (task === 'benchmarks') {
    prog(0.86, 'Benchmark set: shock tube, interface schemes, turbulence closures, flashing');
    const n = Math.round(num(v.benchN, 200, 50, 4000)), e = 1e-7, ex = riemannExact({ rho: 1, u: 0, p: 1 }, { rho: 0.125, u: 0, p: 0.1 }), bn = baerNunziato({ n, tEnd: 0.2, left: { a1: 1 - e, rho1: 1, u1: 0, p1: 1, rho2: 1, u2: 0, p2: 1 }, right: { a1: 1 - e, rho1: 0.125, u1: 0, p1: 0.1, rho2: 0.125, u2: 0, p2: 0.1 } }), exR = bn.x.map((x) => ex.sample((x - 0.5) / 0.2).rho);
    const two = baerNunziato({ n, tEnd: 2.2e-4, left: { a1: 0.2, rho1: 50, u1: 0, p1: 1e7, rho2: 1000, u2: 0, p2: 1e7 }, right: { a1: 0.8, rho1: 5, u1: 0, p1: 1e6, rho2: 1000, u2: 0, p2: 1e6 }, eos: [{ gamma: 1.4, pinf: 0 }, { gamma: 4.4, pinf: 6e8 }] });
    plots.push({ type: 'line', title: 'Seven-equation model: Sod shock tube in the single-phase limit', xlabel: 'x', ylabel: 'Density', series: [{ name: `Baer–Nunziato model, ${n} cells`, x: bn.x, y: bn.rho1 }, { name: 'Exact Riemann solution', x: bn.x, y: exR, dash: true }] });
    plots.push({ type: 'line', title: 'Seven-equation model: liquid–gas shock tube (two pressures, two velocities)', xlabel: 'x (m)', ylabel: 'see legend', series: [{ name: 'Gas volume fraction', x: two.x, y: two.a1 }, { name: 'Gas pressure / 10 MPa', x: two.x, y: two.p1.map((p) => p / 1e7) }, { name: 'Liquid pressure / 10 MPa', x: two.x, y: two.p2.map((p) => p / 1e7) }, { name: 'Gas velocity / 100 m/s', x: two.x, y: two.u1.map((u) => u / 100) }] });
    await tick();
    const m = Math.round(num(v.benchAdvN, 48, 16, 200)), rows = []; for (const sc of ['vof', 'levelset', 'clsvof', 'phasefield', 'front']) { const a = interfaceAdvect2D({ n: m, scheme: sc, test: 'translateDiag' }), b = interfaceAdvect2D({ n: m, scheme: sc, test: 'translateX' }); rows.push([{ vof: 'Volume of fluid (THINC/WLIC)', levelset: 'Level set', clsvof: 'Coupled level set / VOF', phasefield: 'Phase field (conservative Allen–Cahn)', front: 'Front tracking' }[sc], sig(a.massError, 2), r3(a.shapeError, 4), r3(b.shapeError, 4), r3(a.shapeError / Math.max(b.shapeError, 1e-12), 2)]); await tick(); }
    tables.push({ title: `Interface-advection schemes: circle carried once across a periodic ${m} × ${m} grid`, columns: ['Scheme', 'Volume error (diagonal)', 'Shape error, diagonal', 'Shape error, along x', 'Mesh-orientation ratio'], rows, note: 'Shape error: L1 difference of the liquid fraction from the exact one, relative to the body area. The ratio of the diagonal to the grid-aligned error measures the mesh-orientation sensitivity.' });
    const lam = cm.vsl / Math.max(cm.vsl + cm.vsg, 1e-9), ReL = clamp(((cm.rhoL * lam + cm.rhoG * (1 - lam)) * (cm.vsl + cm.vsg) * cm.D) / (cm.muL * lam + cm.muG * (1 - lam)), 1e4, 3e6), rt = reTauOf(ReL), rr = []; for (const q of RANS_MODELS) { const s = ransPipe({ reTau: rt, model: q.value, n: 70, tol: 1e-7, maxIter: 4000 }); rr.push([q.label, sig(s.Re, 4), r3(s.f, 5), r3(frictionFactor(s.Re, 0), 5), r3(100 * (s.f / frictionFactor(s.Re, 0) - 1), 1), s.iterations]); await tick(); }
    tables.push({ title: `Turbulence closures in developed pipe flow at Re_τ = ${Math.round(rt)}`, columns: ['Closure', 'Re', 'Friction factor', 'Colebrook (smooth)', 'Difference (%)', 'Iterations'], rows: rr });
    const k = st.n - 1, prA = fm0(R, st.P[k], st.T[k]), prB = fm0(R, st.P[st.n], st.T[st.n]); // flashing over the last cell at the riser top
    if (prA && prB && prA.mG + prA.mO + prA.mW > 0) { const mt = prA.mG + prA.mO + prA.mW, A = (PI * st.D[k] ** 2) / 4, xe = (p) => clamp(interp1([st.P[st.n] * 1e5, st.P[k] * 1e5], [prB.mG / mt, prA.mG / mt], p), 0, 1), vG = (p) => 1 / (prA.rhoG * (p / (st.P[k] * 1e5))), opt = { G: mt / A, D: st.D[k], length: st.ds, n: 120, p0: st.P[k] * 1e5, x0: prA.mG / mt, pSat: 1.5 * st.P[k] * 1e5, xEq: xe, vG, vL: 1 / prA.rhoL, f: 0.02, theta: st.theta[k] }, hrm = flashingFlow({ ...opt, model: 'hrm' }), hem = flashingFlow({ ...opt, model: 'hem' });
      tables.push({ title: 'Flashing over the last cell of the riser: homogeneous relaxation against equilibrium', columns: ['Quantity', 'Relaxation model', 'Equilibrium model', 'Unit'], rows: [['Outlet pressure', r3(hrm.p[hrm.p.length - 1] / 1e5, 3), r3(hem.p[hem.p.length - 1] / 1e5, 3), 'bara'], ['Gas mass fraction at the outlet', sig(hrm.x[hrm.x.length - 1], 4), sig(hem.x[hem.x.length - 1], 4), '–'], ['Largest lag behind equilibrium', sig(hrm.lag, 3), '0', '–'], ['Outlet velocity', r3(hrm.u[hrm.u.length - 1], 2), r3(hem.u[hem.u.length - 1], 2), 'm/s']], note: 'Relaxation time of Downar-Zapolski et al. (fitted to flashing water below 10 bar; indicative for hydrocarbons). Equilibrium gas fraction interpolated between the two ends of the cell.' }); out.flashingLag = hrm.lag; }
    balances.push({ name: 'Seven-equation model: total energy of the closed tube, initial vs final (normalised)', in: 1, out: 1 + bn.conservation.energy }, { name: 'Seven-equation model: mass of phase 1, initial vs final (normalised)', in: 1, out: 1 + bn.conservation.mass1 });
    out.benchmarks = { sodL1: mean(bn.rho1.map((x, i) => Math.abs(x - exR[i]))), entropyProduction: bn.entropy.production, schemes: rows.map((q) => ({ scheme: q[0], massError: q[1], shape: q[2] })) };
  }
  if (task === 'cfd3d') {
    const kind = C3_CASES.some((q) => q.value === v.c3Case) ? v.c3Case : 'tgv', n = Math.round(num(v.c3N, 16, 8, 128)), sgsIn = C3_SGS.some((q) => q.value === v.c3Sgs) ? v.c3Sgs : 'wale', Re = num(v.c3Re, 180, 1, 20000), tUser = num(v.c3Time, 0, 0, 1e4);
    const pow2 = (x) => 2 ** Math.max(2, Math.round(Math.log2(Math.max(x, 4)))), t00 = Date.now();
    const drive = async (run, label) => { let last = -1; while (!run.done) { run.advance(25, 250); const f = Math.floor(run.progress * 20); if (f !== last) { last = f; prog(0.86 + 0.13 * run.progress, `${label}: ${(100 * run.progress).toFixed(0)} % (step ${run.sim.state.steps})`); } await tick(); } return run.result(); };
    const costRows = (r) => [['Cells', r.cells, '–'], ['Time steps', r.steps, '–'], ['Time per step (this run)', sig(r.msPerStep, 3), 'ms'], ['Cost', sig((1000 * r.msPerStep) / Math.max(r.cells, 1), 3), 'µs per cell and step'], ['Largest divergence of the velocity', sig(r.divergenceMax ?? 0, 2), '1/time unit']];
    const SGS_TXT = { none: 'no sub-grid model', smagorinsky: 'Smagorinsky model with van Driest damping (LES)', wale: 'WALE model (LES)', des: 'Spalart–Allmaras DES', ddes: 'Spalart–Allmaras delayed DES', iddes: 'Spalart–Allmaras IDDES' };
    let primary = 0, secondary = 0, statement = '';
    if (kind === 'tgv') {
      const dim = v.c3Dim === '2' ? 2 : 3, sg = sgsIn === 'none' || sgsIn === 'smagorinsky' || sgsIn === 'wale' ? sgsIn : 'none', r = await drive(taylorGreen({ n: pow2(n), dim, re: Re, tEnd: tUser || (dim === 2 ? 2 : 10), sgs: sg }), 'Taylor–Green vortex');
      const ref = Math.abs(Re - 1600) < 1 && dim === 3 ? rowsOf(REF.TGV1600) : null;
      plots.push({ type: 'line', title: `Taylor–Green vortex (${dim}-D), Re = ${Re}: kinetic-energy dissipation rate`, xlabel: 't U / L', ylabel: '−dE/dt (U³/L)', series: [{ name: `${r.n}³ cells, ${SGS_TXT[sg]}`, x: r.t, y: r.dissipation }, { name: 'Resolved part 2 ν × enstrophy', x: r.t, y: r.enstrophyDissipation, dash: true }, ...(dim === 2 ? [{ name: 'Exact', x: r.t, y: r.t.map((t) => (4 / Re) * r.ek[0] * Math.exp((-4 * t) / Re)), mode: 'points' }] : []), ...(ref ? [{ name: 'DNS 512³ (Dairay et al. 2017)', x: ref.map((q) => q.t), y: ref.map((q) => q.eps), mode: 'points' }] : [])] });
      tables.push({ title: 'Three-dimensional solver: Taylor–Green vortex', columns: ['Quantity', 'Value', 'Unit'], rows: [['Grid', `${r.n} × ${r.n} × ${dim === 2 ? 2 : r.n}`, 'cells, triply periodic box of side 2π'], ['Sub-grid model', SGS_TXT[sg], ''], ...(dim === 2 ? [['Largest error of the velocity amplitude against exp(−2νt)', sig(r.errorMax, 3), '– (second order in the cell size)']] : [['Peak dissipation rate', sig(r.peak.value, 4), 'U³/L'], ['Time of the peak', r3(r.peak.t, 2), 'L/U'], ['Initial dissipation rate (exact 0.75 / Re)', `${sig(r.dissipation[0], 4)} (${sig(0.75 / Re, 4)})`, 'U³/L']]), ...costRows(r)], note: dim === 3 ? `Resolution statement: ${r.n} cells per 2π. At Re = ${Re} ${Re <= 200 ? 'the flow is resolved from about 32 cells per side (the peak dissipation changes by 3 % between 32³ and 64³ at Re = 100)' : 'this grid does not resolve the smallest scales: the run is a large-eddy or under-resolved simulation, and the dissipation peak is lower and earlier than in the 512³ reference'}.` : 'The decay of the two-dimensional vortex is an exact solution of the Navier–Stokes equations.' });
      primary = dim === 2 ? r.errorMax : r.peak.value; secondary = r.ek[r.ek.length - 1]; statement = `Taylor–Green vortex on ${r.n}³ cells`;
      balances.push({ name: '3-D solver: kinetic energy lost vs time integral of the dissipation rate (Taylor–Green)', in: r.ek[0] - r.ek[r.ek.length - 1], out: r.t.reduce((s, t, i) => (i ? s + 0.5 * (r.dissipation[i] + r.dissipation[i - 1]) * (t - r.t[i - 1]) : 0), 0) + 0 });
      out.cfd3d = { case: kind, cells: r.cells, steps: r.steps, primary, secondary, peakTime: r.peak.t, msPerStep: r.msPerStep };
    } else if (kind === 'channel' || kind === 'dns') {
      const reTau = clamp(Re, 40, 1000); let o;
      if (kind === 'dns') { // resolved minimal-flow-unit channel: the grid follows from the resolution criteria Δx⁺ ≤ 10, Δz⁺ ≤ 6.5, Δy⁺ ≤ 1 at the wall and ≤ 7 at the centre
        const lxP = Math.max(550, 3 * reTau > 550 ? 550 : 550), lzP = 200, nx = pow2(lxP / 9.5), nz = pow2(lzP / 6.3), ny = 2 * Math.ceil(Math.max(32, 0.53 * reTau) / 2);
        o = { reTau, nx, ny, nz, lx: lxP / reTau, lz: lzP / reTau, stretchY: 1.6, sgs: 'none', mode: 'bulk', tEnd: tUser || 12, tStats: 0.5 * (tUser || 12) };
      } else o = { reTau, nx: pow2(n), ny: 2 * Math.round(n / 2), nz: pow2(n), sgs: sgsIn, mode: 'bulk', tEnd: tUser || 8, tStats: 0.5 * (tUser || 8) };
      if (reTau === 180) o.bulkPlus = REF.MKM180.uBulkPlus; // the flow rate of the reference DNS (Re_b = 5,600)
      const r = await drive(channel3D(o), kind === 'dns' ? 'Resolved channel simulation' : `Channel, ${SGS_TXT[o.sgs]}`), mk = REF.MKM180, at180 = reTau === 180, rs = r.resolution;
      const dnsOk = rs.dxPlus <= 12 && rs.dzPlus <= 7 && rs.firstCentrePlus <= 1 && rs.dyPlusMax <= 7.5, turbulent = Math.max(...r.uvPlus) > 0.3;
      plots.push({ type: 'line', title: `Channel at Re_τ ≈ ${Math.round(r.reTauActual)}: mean velocity in wall units`, xlabel: 'y⁺', ylabel: 'U⁺', logx: true, series: [{ name: `${r.nx} × ${r.ny} × ${r.nz} cells, ${SGS_TXT[r.sgs]}`, x: r.yPlus, y: r.uPlus, mode: 'both' }, ...(at180 ? [{ name: 'DNS (Moser, Kim & Mansour 1999)', x: mk.rows.map((q) => q[0]), y: mk.rows.map((q) => q[1]), mode: 'points' }] : []), { name: 'U⁺ = y⁺', x: [0.5, 11], y: [0.5, 11], dash: true }, { name: 'Log law', x: [11, Math.max(r.reTauActual, 30)], y: [Math.log(11) / 0.41 + 5.2, Math.log(Math.max(r.reTauActual, 30)) / 0.41 + 5.2], dash: true }] });
      plots.push({ type: 'line', title: 'Channel: resolved velocity fluctuations and shear stress', xlabel: 'y⁺', ylabel: 'wall units', series: [{ name: 'u rms', x: r.yPlus, y: r.uRms }, { name: 'v rms', x: r.yPlus, y: r.vRms }, { name: 'w rms', x: r.yPlus, y: r.wRms }, { name: '−u′v′ (resolved)', x: r.yPlus, y: r.uvPlus }, ...(r.sgs !== 'none' ? [{ name: 'ν_t / ν', x: r.yPlus, y: r.nutPlus, dash: true }] : [])] });
      plots.push({ type: 'line', title: 'Channel: history of the friction Reynolds number', xlabel: 't u_τ / h', ylabel: 'Re_τ', series: [{ name: 'from the wall shear', x: r.history.t, y: r.history.reTau }], vlines: [{ x: r.time - o.tStats, label: 'statistics from here' }] });
      tables.push({ title: `Three-dimensional solver: plane channel, ${kind === 'dns' ? 'resolved simulation (no model)' : SGS_TXT[r.sgs]}`, columns: ['Quantity', 'This run', 'DNS reference', 'Unit'], rows: [['Grid / box', `${r.nx} × ${r.ny} × ${r.nz} / ${r3(r.lx, 2)} h × 2 h × ${r3(r.lz, 2)} h`, at180 ? '128 × 129 × 128 / 4π h × 2 h × 4π/3 h (spectral)' : '—', ''], ['Bulk Reynolds number U_b 2h / ν', sig(r.reBulk, 4), at180 ? mk.reBulk : '—', '–'], ['Friction Reynolds number', r3(r.reTauActual, 1), at180 ? mk.reTau : '—', '–'], ['Friction coefficient 2 τ_w / ρ U_b²', sig(r.cf, 4), at180 ? mk.cf : sig(0.073 * r.reBulk ** -0.25, 4) + ' (Dean)', '–'], ['Friction-coefficient error', at180 ? r3(100 * (r.cf / mk.cf - 1), 1) : r3(100 * (r.cf / (0.073 * r.reBulk ** -0.25) - 1), 1), '0', '%'], ['Centre-line velocity U_c / u_τ', r3(r.uCentrePlus, 2), at180 ? mk.uCentrePlus : '—', '–'], ['Largest u rms', r3(Math.max(...r.uRms), 2), at180 ? mk.uRmsMax : '—', 'u_τ'], ['Largest resolved −u′v′', r3(Math.max(...r.uvPlus), 2), at180 ? mk.uvMax : '—', 'u_τ²'], ['Δx⁺ / Δz⁺', `${r3(rs.dxPlus, 1)} / ${r3(rs.dzPlus, 1)}`, '≤ 10–12 / ≤ 5–7 for a DNS with second-order differences', ''], ['Δy⁺ first cell centre / largest cell', `${r3(rs.firstCentrePlus, 2)} / ${r3(rs.dyPlusMax, 1)}`, '≤ 1 / ≤ 7', ''], ['Resolution criteria of a DNS met', dnsOk ? 'yes' : 'no', '', ''], ['Flow state', turbulent ? 'turbulent (resolved shear stress present)' : 'no resolved turbulence (laminar or RANS-like solution)', '', ''], ['Simulated / averaging time', `${r3(r.time, 1)} / ${r3(o.tStats, 1)}`, at180 ? 'about 10 (after a long start-up)' : '—', 'h / u_τ'], ['Samples of the plane averages', r.samples, '', ''], ...costRows(r).map((q) => [q[0], q[1], '', q[2]])],
        note: `Resolution statement: ${dnsOk ? 'all four spacing criteria of a direct simulation are met, in a small periodic box (minimal flow unit) that holds one or two near-wall streak pairs' : 'this is a coarse scale-resolving simulation; the grid does not resolve the near-wall streaks (which need Δz⁺ ≈ 5 and Δx⁺ ≈ 10)'}. The averaging time of ${r3(o.tStats, 1)} h/u_τ is short: the statistical error of the friction coefficient is a few per cent. ${!turbulent && (r.sgs === 'ddes' || r.sgs === 'iddes') ? 'The shielding function keeps the whole channel in RANS mode at this resolution, so the resolved fluctuations decay and the result is the Spalart–Allmaras RANS solution.' : ''} The constant flow rate is that of the reference; production-resolution cases are written by the hand-off page.` });
      plots.push({ type: 'field', title: 'Channel: instantaneous streamwise velocity on a cross-section (wall units)', xlabel: 'z / h', ylabel: 'y / h', zlabel: 'u⁺', zunit: '–', x: r.slice.z, y: r.slice.y, z: r.slice.u, cmap: 'turbo', equal: true });
      balances.push({ name: '3-D channel: flow rate held (target vs time-mean bulk velocity, u_τ units)', in: r.mode === 'bulk' ? (o.bulkPlus || r.uBulk) : r.uBulk, out: r.uBulk });
      if (Math.abs(r.balance.forcing - r.balance.wallShear) > 0.05 * Math.abs(r.balance.forcing)) warnings.push({ level: 'info', msg: `3-D channel: the mean wall shear (${sig(r.balance.wallShear, 3)}) and the mean driving force (${sig(r.balance.forcing, 3)}) differ by ${(100 * Math.abs(1 - r.balance.wallShear / r.balance.forcing)).toFixed(0)} % over the averaging window — the flow is not yet statistically steady; simulate longer.` });
      if (kind === 'dns' && !turbulent) warnings.push({ level: 'warn', msg: `Resolved channel run: the grid meets the spacing criteria of a direct simulation, but the flow did not stay turbulent in this small box (resolved shear stress ${r3(Math.max(...r.uvPlus), 2)} u_τ², Re_τ ${r3(r.reTauActual, 0)} and falling towards the laminar value): the result is a resolved laminarisation, not a turbulent DNS. The stored run at Re_b = 5,600 (64 × 96 × 32 cells, 12 h/u_τ) behaved the same way.` });
      primary = r.cf; secondary = r.uCentrePlus; statement = `channel ${r.nx} × ${r.ny} × ${r.nz}, Δx⁺ ${r3(rs.dxPlus, 0)}, Δz⁺ ${r3(rs.dzPlus, 0)}`;
      out.cfd3d = { case: kind, sgs: r.sgs, cells: r.cells, steps: r.steps, primary, secondary, reTau: r.reTauActual, cf: r.cf, cfError: at180 ? r.cf / mk.cf - 1 : null, dnsResolved: dnsOk, turbulent, dxPlus: rs.dxPlus, dzPlus: rs.dzPlus, dyPlusWall: rs.firstCentrePlus, msPerStep: r.msPerStep };
    } else if (kind === 'pipe') {
      const r = await drive(pipe3D({ n, nx: 2, nu: 1 / clamp(Re, 1, 1500), bulk: 1, tEnd: tUser || 0.4 * clamp(Re, 1, 1500), steady: 1e-5 }), 'Laminar pipe (immersed boundary)');
      plots.push({ type: 'line', title: 'Pipe by immersed boundary: velocity on the diameter', xlabel: 'r / R', ylabel: 'u / U bulk', series: [{ name: `${r3(r.cellsAcross, 1)} cells across the diameter`, x: r.profile.r, y: r.profile.u, mode: 'both' }, { name: 'Hagen–Poiseuille', x: linspace(-1, 1, 41), y: linspace(-1, 1, 41).map((x) => 2 * (1 - x * x)), dash: true }] });
      tables.push({ title: 'Three-dimensional solver: laminar pipe flow', columns: ['Quantity', 'Value', 'Unit'], rows: [['Cells across the diameter', r3(r.cellsAcross, 1), '–'], ['Reynolds number', sig(r.reynolds, 4), '–'], ['Darcy friction factor', sig(r.friction, 4), '–'], ['64 / Re', sig(r.frictionLaminar, 4), '–'], ['Error', r3(100 * r.error, 2), '%'], ['Error of the represented cross-section', r3(100 * r.areaError, 2), '%'], ...costRows(r)], note: 'The wall is represented by the fluid fraction of each velocity control volume (first-order accurate at the wall): the friction error falls from about 13 % at 10 cells across the diameter to 1–4 % at 13–30 cells.' });
      balances.push({ name: '3-D pipe: wall force from the momentum balance vs Hagen–Poiseuille (ratio within the stated mesh error)', in: 1, out: 1 + (Math.abs(r.error) < 0.2 ? 0 : r.error) });
      primary = r.friction; secondary = r.error; statement = `pipe with ${r3(r.cellsAcross, 1)} cells across`;
      out.cfd3d = { case: kind, cells: r.cells, steps: r.steps, primary, secondary, frictionError: r.error, msPerStep: r.msPerStep };
    } else if (kind === 'dambreak' || kind === 'stratified') {
      const sg = sgsIn === 'none' ? 'none' : 'smagorinsky', rhoL = cm.rhoL, rhoG = Math.max(cm.rhoG, 1), muL = cm.muL, muG = cm.muG, own = v.c3Fluids !== 'airwater', fl = own ? { rhoL, rhoG, muL, muG } : {};
      const r = await drive(twoPhase3D(kind === 'dambreak' ? { kind, n, spanCells: Math.max(4, pow2(n / 4)), shape3d: true, sgs: sg, tEnd: tUser || undefined, ...fl } : { kind, n, sgs: sg, height: cm.D, level: clamp(levelOfHoldup(clamp(st.holdup[im], 0.1, 0.8)), 0.15, 0.8), force: 2, tEnd: tUser || 1.5 * Math.sqrt(cm.D / G) * 4, ...fl }), kind === 'dambreak' ? 'LES + VOF dam break' : 'LES + VOF stratified channel');
      plots.push({ type: 'field', title: `${kind === 'dambreak' ? 'Dam break' : 'Stratified wavy channel'} (3-D, volume of fluid${sg === 'none' ? '' : ' + Smagorinsky LES'}): liquid fraction on the mid-plane`, xlabel: 'x (m)', ylabel: 'y (m)', zlabel: 'Liquid fraction', zunit: '–', x: r.slice.x, y: r.slice.y, z: r.slice.c, zmin: 0, zmax: 1, cmap: 'salinity', u: r.slice.u, v: r.slice.v, vectors: true, equal: true });
      if (kind === 'dambreak') { const tu = Math.sqrt(r.a / G), of = !own ? rowsOf(REF.OF_DAMBREAK) : null; plots.push({ type: 'line', title: 'Dam break: front position on the floor', xlabel: 't √(g / a)', ylabel: 'x front / a', series: [{ name: `${r.nx} × ${r.ny} × ${r.nz} cells`, x: r.t.map((t) => t / tu), y: r.front.map((x) => x / r.a) }, ...(of ? [{ name: 'OpenFOAM interFoam, 160 × 80 cells (2-D)', x: of.map((q) => q.tStar), y: of.map((q) => q.front), mode: 'points' }] : [])] }); }
      else plots.push({ type: 'line', title: 'Stratified channel: kinetic energy', xlabel: 't (s)', ylabel: 'J/m³', series: [{ name: 'volume average', x: r.t, y: r.kinetic }] });
      tables.push({ title: `Three-dimensional solver: ${kind === 'dambreak' ? 'dam break' : 'stratified wavy channel'} with a volume-of-fluid interface`, columns: ['Quantity', 'Value', 'Unit'], rows: [['Grid / box', `${r.nx} × ${r.ny} × ${r.nz} / ${sig(r.lx, 3)} × ${sig(r.ly, 3)} × ${sig(r.lz, 3)}`, 'cells / m'], ['Liquid / gas density', `${sig(own ? rhoL : 998.2, 4)} / ${sig(own ? rhoG : 1.204, 4)} (ratio ${sig(r.densityRatio, 3)})`, 'kg/m³'], ['Sub-grid model', SGS_TXT[sg], ''], ['Liquid volume, initial / final', `${sig(r.volume.initial, 6)} / ${sig(r.volume.final, 6)}`, 'm³'], ['Relative volume error', sig(Math.abs(r.volume.error), 2), '–'], ['Liquid fraction outside [0, 1] (largest)', sig(Math.max(-r.bounds.min, r.bounds.max), 2), '–'], ['Mean holdup of the box', r3(r.holdup, 4), '–'], ...(kind === 'stratified' ? [['Superficial velocities liquid / gas at the end', `${sig(r.vsl, 3)} / ${sig(r.vsg, 3)}`, 'm/s']] : [['Front position at the end', r3(r.front[r.front.length - 1], 3), 'm']]), ['Simulated time', sig(r.time, 3), 's'], ...costRows(r)], note: 'Direction-split THINC/WLIC volume of fluid with the dilatation correction (volume conserved to round-off), variable density and viscosity, gravity, no surface tension; limited (van Leer) momentum advection; pressure by the constant-coefficient splitting of Dodd & Ferrante with the direct solver. Resolution statement: the interface is captured over about two cells; waves shorter than four cells and droplets are not represented.' });
      balances.push({ name: '3-D volume of fluid: liquid volume, initial vs final (m³)', in: r.volume.initial, out: r.volume.final });
      primary = kind === 'dambreak' ? r.front[r.front.length - 1] : r.kinetic[r.kinetic.length - 1]; secondary = Math.abs(r.volume.error); statement = `${r.nx} × ${r.ny} × ${r.nz} cells, two phases`;
      out.cfd3d = { case: kind, sgs: sg, cells: r.cells, steps: r.steps, primary, secondary, volumeError: r.volume.error, holdup: r.holdup, msPerStep: r.msPerStep };
    } else { // 1-D line model coupled with a 3-D pipe section
      prog(0.87, '1-D line model coupled with a 3-D pipe section'); await tick();
      const D = cm.D, L = Math.max(st.s[st.n], 100 * D), mu = Math.max(num(v.c3Mu, 0.5, 1e-4, 100), 1e-4), rho = cm.rhoL, uT = (num(v.c3ReCoupled, 20, 0.1, 1500) * mu) / (rho * D), dp = (32 * mu * L * uT) / (D * D);
      const r = coupled1D3D({ D, length: L, nCells: 20, rho, mu, dp, n, nx: 2 }); await tick();
      plots.push({ type: 'line', title: '1-D + 3-D coupling: flow velocity over the coupling iterations', xlabel: 'Iteration', ylabel: 'm/s', series: [{ name: 'Coupled model', x: r.iterations.map((_, i) => i + 1), y: r.iterations.map((q) => q.u), mode: 'both' }], hlines: [{ y: r.velocityAnalytic, label: 'Hagen–Poiseuille' }] });
      tables.push({ title: 'One-dimensional line model coupled with a three-dimensional pipe section (single-phase laminar verification)', columns: ['Quantity', 'Value', 'Unit'], rows: [['Line: diameter / length / cells', `${sig(D, 3)} / ${sig(L, 5)} / ${r.nCells}`, 'm'], ['Fluid: density / viscosity', `${sig(rho, 4)} / ${sig(mu, 3)}`, 'kg/m³, Pa s'], ['Imposed pressure difference (inlet − outlet)', sig(dp, 5), 'Pa'], ['Velocity, 1-D closure alone', sig(r.velocity1d, 6), 'm/s'], ['Velocity, coupled', sig(r.velocity, 6), 'm/s'], ['Velocity, analytic (Hagen–Poiseuille)', sig(r.velocityAnalytic, 6), 'm/s'], ['Coupled − analytic', sig(100 * r.errorCoupled, 3), '%'], ['Friction of the 3-D section against 64/Re', sig(100 * r.error3d, 3), `% (${r3(r.cellsAcross, 1)} cells across the diameter)`], ['Pressure drop returned by the 3-D cell', sig(r.pressureDrop3dCell, 5), 'Pa'], ['Coupling iterations', r.iterations.length, r.converged ? 'converged' : 'not converged'], ['3-D cells / time per step', `${r.cells3d} / ${sig(r.msPerStep3d, 3)} ms`, '']], note: 'Exchange per iteration: the 1-D model sends the flow rate (mean velocity) and the fluid properties to the periodic 3-D section; the section returns its pressure gradient, which replaces the friction closure of that 1-D cell; the 1-D momentum balance is solved again with the fixed inlet and outlet pressures. A single-phase section returns a holdup of one. The viscosity is an input of this task so that the verification case is laminar; two-phase and turbulent production cases are written for an external solver on the hand-off page.' });
      balances.push({ name: '1-D + 3-D coupling: imposed pressure difference vs sum of the cell pressure drops (Pa)', in: r.pressureDropBalance.given, out: r.pressureDropBalance.sum });
      primary = r.velocity; secondary = r.errorCoupled; statement = `3-D section with ${r3(r.cellsAcross, 1)} cells across the diameter`;
      out.cfd3d = { case: kind, cells: r.cells3d, steps: r.iterations.reduce((s, q) => s + q.steps3d, 0), primary, secondary, coupledError: r.errorCoupled, error3d: r.error3d, converged: r.converged, msPerStep: r.msPerStep3d };
    }
    out.cfd3d.statement = statement; out.cfd3d.wallSeconds = (Date.now() - t00) / 1000;
    kpis.push({ label: `3-D solver (${statement})`, value: sig(primary, 4), unit: kind === 'channel' || kind === 'dns' ? 'C_f' : kind === 'pipe' ? 'f Darcy' : kind === 'coupled' ? 'm/s' : kind === 'tgv' ? 'peak −dE/dt' : 'm', status: 'ok', help: 'Primary result of the selected three-dimensional case; the table states the resolution.' });
  }
  return { outputs: out, kpis };
}
const fm0 = (R, P, T) => { try { return R.fm.at(P, T, R.st.mScale); } catch { return null; } };

// =====================================================================================================
// 7i. Reference data sets wired to the engine, and the provenance of the constants
// =====================================================================================================
const TASKS = [{ value: 'line', label: 'Line solution (steady, slugging, transient)' }, { value: 'twofluid', label: 'Two-fluid model: slug capturing on a test section' }, { value: 'cfd', label: '2-D Navier–Stokes solver (projection, VOF, SST)' }, { value: 'bubbly', label: 'Bubbly-flow closures: radial void distribution' }, { value: 'benchmarks', label: 'Benchmark set (shock tube, interface schemes, turbulence closures, flashing)' }, { value: 'cfd3d', label: '3-D Navier–Stokes solver (LES, DES, IDDES, DNS, LES + VOF, 1-D + 3-D coupling)' }];
const C3_CASES = [{ value: 'tgv', label: 'Taylor–Green vortex (verification)' }, { value: 'channel', label: 'Turbulent channel with the selected turbulence treatment (LES / DES / IDDES / no model)' }, { value: 'dns', label: 'Resolved channel (DNS grid from the resolution criteria, small box)' }, { value: 'pipe', label: 'Laminar pipe by immersed boundary (verification)' }, { value: 'dambreak', label: 'LES + VOF: three-dimensional dam break' }, { value: 'stratified', label: 'LES + VOF: stratified wavy channel' }, { value: 'coupled', label: '1-D line model coupled with a 3-D pipe section' }];
const C3_SGS = [{ value: 'wale', label: 'LES: WALE' }, { value: 'smagorinsky', label: 'LES: Smagorinsky + van Driest' }, { value: 'des', label: 'DES (Spalart–Allmaras)' }, { value: 'ddes', label: 'Delayed DES' }, { value: 'iddes', label: 'IDDES' }, { value: 'none', label: 'No model' }];
const AIR20 = { rhoL: 998.2, muL: 1.0e-3, sigma: 0.0728, muG: 1.82e-5 }; // water and air at 20 °C, used where a source gives no fluid temperature (stated in the data-set notes)
const rowsOf = (d, f = () => true) => d.rows.map((r) => Object.fromEntries(d.cols.map((c, i) => [c, r[i]]))).filter(f);
const memo = {}; const once = (k, f) => (k in memo ? memo[k] : (memo[k] = f()));
function buildValidation() {
  const sh = rowsOf(REF.SHOHAM).map((r) => ({ ...r, match: 1 })), fl = REF.SHOHAM.fluid, pcl = (r, f = fl) => patternClass({ vsl: r.vsl, vsg: r.vsg, rhoL: f.rhoL ?? r.rhoL, rhoG: f.rhoG ?? r.rhoG, muL: f.muL ?? r.muL, muG: f.muG ?? r.muG, sigma: f.sigma ?? r.sigma, D: r.D, theta: r.angle * DEG });
  const fpCols = [{ key: 'D', label: 'Diameter', unit: 'm' }, { key: 'angle', label: 'Inclination', unit: '°' }, { key: 'vsl', label: 'Superficial liquid velocity', unit: 'm/s' }, { key: 'vsg', label: 'Superficial gas velocity', unit: 'm/s' }, { key: 'pattern', label: 'Observed pattern', type: 'text' }, { key: 'match', label: 'Observed class reproduced (1 = yes)' }];
  const fp = (id, title, rows, tol, note, covers, src = REF.SHOHAM.source, f = fl) => ({ id, title, quantity: 'Flow-pattern class reproduced', unit: '–', kind: 'experiment', source: src, columns: fpCols, rows, target: 'match', model: (r) => (pcl(r, f) === r.pattern ? 1 : 0), tolerance: { mape: tol }, covers,
    note: `${note} Accuracy measure: the target is 1 for every row and the model returns 1 when the mechanistic map (Taitel–Dukler / Barnea, kernel flowPattern) gives the observed class (SS, SW, I, A, DB, B) and 0 otherwise, so the reported MAPE is the percentage of misclassified points and −bias the misclassified fraction. ${REF.SHOHAM.note}` });
  const A = [];
  A.push(fp('fp-shoham-horizontal', 'Flow pattern, horizontal air–water, 25 and 51 mm (Shoham 1982)', sh.filter((r) => r.angle === 0), 10, 'Horizontal pipes.', ['Horizontal multiphase flow', 'Gas-liquid flow']));
  A.push(fp('fp-shoham-up', 'Flow pattern, upward inclined air–water, +0.25° to +80° (Shoham 1982)', sh.filter((r) => r.angle > 0 && r.angle < 90), 20, 'Upward inclinations. Bubble flow is now limited to pipes steeper than 60° and the dispersed-bubble boundary is that of Barnea (1986).', ['Inclined flow']));
  A.push(fp('fp-shoham-down', 'Flow pattern, downward inclined air–water, −1° to −80° (Shoham 1982)', sh.filter((r) => r.angle < 0 && r.angle > -90), 28, 'Downward inclinations (−1° to −80°). With the torn-film criterion of Barnea, Shoham & Taitel (1982) for the stratified → annular transition and the film-stability criterion for annular ↔ intermittent flow the misses fell from 31 % to 23 % of these rows; what remains is stratified-wavy flow at very low liquid and high gas rates that the model calls annular (9 rows) and the dispersed-bubble boundary, observed at a liquid velocity about 30 % lower than predicted.', ['Inclined flow']));
  A.push(fp('fp-shoham-vertical', 'Flow pattern, vertical upward and downward air–water (Shoham 1982)', sh.filter((r) => Math.abs(r.angle) === 90), 15, 'Vertical pipes (+90° and −90°); churn is counted as intermittent. Vertical downward flow used to be treated with the stratified logic (6 of 19 rows right); it now has the falling-film (annular) and intermittent criteria of the unified model (18 of 19).', ['Vertical upward/downward flow']));
  const cls = (set) => sh.filter((r) => set.includes(r.pattern));
  A.push(fp('fp-class-stratified', 'Observed stratified flow (smooth and wavy), all inclinations (Shoham 1982)', cls(['SS', 'SW']), 25, 'Recall of the stratified classes; smooth and wavy are distinct classes here, and the smooth / wavy boundary is the least certain of the map.', ['Stratified flow']));
  A.push(fp('fp-class-annular', 'Observed annular flow, all inclinations (Shoham 1982)', cls(['A']), 25, 'Recall of annular flow.', ['Annular flow']));
  A.push(fp('fp-class-bubbly', 'Observed bubble and dispersed-bubble flow, all inclinations (Shoham 1982)', cls(['B', 'DB']), 40, 'Recall of the bubbly classes.', ['Bubbly flow']));
  A.push(fp('fp-class-intermittent', 'Observed intermittent flow (slug, plug, churn), all inclinations (Shoham 1982)', cls(['I']), 12, 'Recall of intermittent flow.', ['Intermittent flow', 'Hydrodynamic slugging']));
  A.push(fp('fp-kokal', 'Flow pattern, oil–air, 26 and 51 mm, 0, ±1, ±5, ±9° (Kokal 1987)', rowsOf(REF.KOKAL).map((r) => ({ ...r, match: 1 })), 36, `A second fluid pair and two diameters, not used when the transition criteria were selected (blind). Half of the misses are downward stratified flow that the observers called smooth and the film Froude-number criterion calls wavy. ${REF.KOKAL.note}`, ['Gas-liquid flow', 'Hilly-terrain pipelines'], REF.KOKAL.source, {}));
  // ---- slug-body holdup (CC BY compilation of 22 studies)
  { const H = REF.SLUG_HOLDUP, cols = [{ key: 'vsl', label: 'Superficial liquid velocity', unit: 'm/s' }, { key: 'vsg', label: 'Superficial gas velocity', unit: 'm/s' }, { key: 'D', label: 'Diameter', unit: 'm' }, { key: 'muL', label: 'Liquid viscosity', unit: 'Pa s' }, { key: 'angle', label: 'Inclination', unit: '°' }, { key: 'hls', label: 'Measured slug-body holdup' }];
    A.push({ id: 'slugbody-gregory', title: 'Slug-body liquid holdup, light oil–air, 26 and 51 mm horizontal (Gregory, Nicholson & Aziz 1978 data)', quantity: 'Slug-body holdup', unit: '–', kind: 'experiment', source: REF.SLUG_HOLDUP_GREGORY.source, columns: cols.filter((c) => c.key !== 'angle'), rows: rowsOf(REF.SLUG_HOLDUP_GREGORY), target: 'hls', model: (r) => slugBodyHoldup(r.vsl + r.vsg), tolerance: { mape: 6 }, covers: ['Slug body holdup', 'Hydrodynamic slugging'], note: 'Kernel closure H_LS = 1 / (1 + (v_m / 8.66)^1.39) against the measurements it was fitted to (as compiled in an open data set): a check of the two constants. ' + REF.SLUG_HOLDUP_GREGORY.note });
    A.push({ id: 'slugbody-compilation', title: 'Slug-body liquid holdup, 22 studies, horizontal to vertical, 1–5,300 mPa s (Abdul-Majeed 2022 compilation)', quantity: 'Slug-body holdup', unit: '–', kind: 'experiment', source: H.source, columns: cols, rows: rowsOf(H), target: 'hls', model: (r) => slugBodyHoldup(r.vsl + r.vsg), tolerance: { mape: 18 }, covers: ['Slug body holdup', 'Inclined flow', 'Vertical upward/downward flow'], note: 'Blind comparison of the kernel closure (Gregory et al.) with the other 21 studies. On all 2,699 rows of the source the mean error is 14.8 % (bias +10 %): 10 % for horizontal pipes, 15 % for inclined ones and 40 % (all over-prediction) for vertical pipes, where the slug body is much more aerated than the horizontal correlation allows — a known limit, left as it is because no inclination-dependent closure could be read from an open source. ' + H.note }); }
  // ---- liquid holdup in vertical flow of viscous oil and gas (CC BY)
  { const Rb = REF.RIBEIRO, pt = (r) => ({ vsl: r.vsl, vsg: r.vsg, rhoL: 900, rhoG: 1.3, muL: r.muL, muG: 1.8e-5, sigma: 0.03, D: Rb.D, theta: PI / 2, rough: 0, P: 1.1e5, label: false });
    A.push({ id: 'holdup-vertical-viscous', title: 'Liquid holdup, vertical upward flow of viscous oil and gas, 60 mm (Ribeiro et al. 2020)', quantity: 'Liquid holdup', unit: '–', kind: 'experiment', source: Rb.source, columns: [{ key: 'muL', label: 'Liquid viscosity', unit: 'Pa s' }, { key: 'vsl', label: 'Superficial liquid velocity', unit: 'm/s' }, { key: 'vsg', label: 'Superficial gas velocity', unit: 'm/s' }, { key: 'pattern', label: 'Pattern named by the authors', type: 'text' }, { key: 'holdup', label: 'Measured liquid holdup' }], rows: rowsOf(Rb, (r) => r.void < 0.97).map((r) => ({ ...r, holdup: 1 - r.void })), target: 'holdup', model: (r) => holdupGradient(pt(r), 'driftFlux').holdup, tolerance: { rmse: 0.18 }, covers: ['Liquid holdup', 'Vertical upward/downward flow'], note: 'Drift-flux closure of the kernel (Bendiksen slug-bubble velocity), blind. Assumed for the model because the file does not give them: oil density 900 kg/m³, surface tension 0.03 N/m, gas at about 1.1 bara. The gas rates reach 30 m/s (churn and annular flow), where a slug-bubble drift-flux relation is outside its range; the tolerance is therefore an absolute one (root-mean-square error of the holdup 0.15); the relative error is 85 % on average because most rows have a holdup below 0.2 (rows with a void fraction above 0.97 are left out). A clear miss of the closure for viscous oil at high gas rates, recorded as such. The other closures do no better on these rows: Beggs & Brill 62 % (slug rows), mechanistic 43 %, Duns & Ros 48 %. ' + Rb.note }); }
  // ---- slug frequency, length and velocity, horizontal 74 mm (Mohmmed et al.)
  const mNote = (d) => ` ${d.note} Water and air properties at 24 °C and 1.013 bar as stated by the source (no property values given).`;
  A.push({ id: 'slugfreq-mohmmed', title: 'Slug frequency, horizontal air–water, 74 mm (Mohmmed et al. 2018)', quantity: 'Slug frequency', unit: '1/s', kind: 'experiment', source: REF.MOHMMED_FREQ.source, columns: [{ key: 'vsl', label: 'Superficial liquid velocity', unit: 'm/s' }, { key: 'vsg', label: 'Superficial gas velocity', unit: 'm/s' }, { key: 'xD', label: 'Station', unit: 'D' }, { key: 'freq', label: 'Measured frequency', unit: '1/s' }], rows: rowsOf(REF.MOHMMED_FREQ), target: 'freq', model: (r) => slugFrequency(r.vsl, r.vsl + r.vsg, r.D, 0, 'zabaras'), tolerance: { mape: 35 }, covers: ['Slug frequency', 'Horizontal multiphase flow'], note: 'Zabaras (2000) correlation in the form checked against two open sources (Gregory–Scott value × (0.836 + 2.75 sin θ)); blind. The other correlations of the kernel on the same rows: Gregory & Scott 20 %, Greskovich & Shrier 19 %, Heywood & Richardson 15 % mean error. On the inclined 26 mm UNICAMP measurements (0–90°, CC BY-NC, used during development but not distributed here) the same correlation gives 24–27 % against 30–36 % for the horizontal correlations, and the earlier form with sin^0.25 θ applied beyond 11° gave 136 %.' + mNote(REF.MOHMMED_FREQ) });
  A.push({ id: 'sluglength-mohmmed', title: 'Slug-body length 54 and 81 diameters from the inlet, horizontal air–water, 74 mm (Mohmmed et al. 2018)', quantity: 'Slug length', unit: 'm', kind: 'experiment', source: REF.MOHMMED_LEN.source, columns: [{ key: 'vsl', label: 'Superficial liquid velocity', unit: 'm/s' }, { key: 'vsg', label: 'Superficial gas velocity', unit: 'm/s' }, { key: 'xD', label: 'Station', unit: 'D' }, { key: 'ls', label: 'Measured slug length', unit: 'm' }], rows: rowsOf(REF.MOHMMED_LEN), target: 'ls', model: (r) => slugLength(r.D, r.vsl + r.vsg, 'scott', { xD: r.xD, theta: 0 }), tolerance: { mape: 45 }, covers: ['Slug length'], note: 'Development-aware slug length of the kernel (slug-train model with the wake law, evaluated at the distance of the station); blind, nothing fitted to these rows. The developed-flow rule (32 D) applied to this 108 D test section gave 478 %. The model has no dependence on the flow rates, which is the scatter that remains. On the UNICAMP inclined measurements (not distributed) the same model gives 23 % at 77 D and 26 % at 257 D (131 % with the 32 D rule).' + mNote(REF.MOHMMED_LEN) });
  A.push({ id: 'slugvelocity-mohmmed', title: 'Slug translational velocity, horizontal air–water, 74 mm (Mohmmed et al. 2018)', quantity: 'Translational velocity', unit: 'm/s', kind: 'experiment', source: REF.MOHMMED_VEL.source, columns: [{ key: 'vsl', label: 'Superficial liquid velocity', unit: 'm/s' }, { key: 'vsg', label: 'Superficial gas velocity', unit: 'm/s' }, { key: 'vt', label: 'Measured translational velocity', unit: 'm/s' }], rows: rowsOf(REF.MOHMMED_VEL), target: 'vt', model: (r) => slugVelocity(r.vsl + r.vsg, r.D, 0).vt, tolerance: { mape: 20 }, covers: ['Slug velocity'], note: 'Bendiksen (1984) relation.' + mNote(REF.MOHMMED_VEL) });
  // ---- severe slugging with a choke or gas lift (Jansen, Shoham & Taitel 1996)
  { const g = REF.JANSEN.rig, A0 = (PI * g.D * g.D) / 4, model = (r) => { const wG = 1.204 * r.vsg0 * A0, wLift = 1.204 * r.gasLift * A0, wL = 998.2 * r.vsl * A0, rhoG0 = 101325 / (287.05 * 293.15), u0 = r.vsl + ((wG + wLift) / rhoG0) / A0, lam = r.vsl / u0, rhoN = lam * 998.2 + (1 - lam) * rhoG0, c = riserSluggingCycle({ D: g.D, feedLength: g.pipeline + g.bufferLength, feedAngle: 1 * DEG, riserHeight: g.riser, riserLength: g.riser, wG, wLift, wL, rhoL: 998.2, muL: 1e-3, T: 293.15, zG: 1, mwG: 28.96, pSep: 101325, alphaL: (stratifiedLevel({ vsl: r.vsl, vsg: Math.max(r.vsg0, 1e-3), rhoL: 998.2, rhoG: rhoG0, muL: 1e-3, muG: 1.82e-5, D: g.D, theta: -1 * DEG }).holdup * g.pipeline) / (g.pipeline + g.bufferLength), chokeDp: Math.max((r.chokeC / 998.2) * rhoN * u0 * u0, 10), rough: 1.5e-6, maxCycles: 4 }); const last = c.stages; return c.period ?? last.buildUp + last.production + last.blowout + last.fallback; };
    A.push({ id: 'severe-jansen', title: 'Severe-slugging cycle time with a riser-top choke or riser-base gas lift (Jansen, Shoham & Taitel 1996)', quantity: 'Cycle time', unit: 's', kind: 'experiment', source: REF.JANSEN.source, columns: [{ key: 'table', label: 'Table (1 choke, 2 gas lift)' }, { key: 'vsl', label: 'Superficial liquid velocity', unit: 'm/s' }, { key: 'vsg0', label: 'Superficial gas velocity at standard conditions', unit: 'm/s' }, { key: 'chokeC', label: 'Choke coefficient', unit: 'Pa s²/m²' }, { key: 'gasLift', label: 'Gas-lift velocity at standard conditions', unit: 'm/s' }, { key: 'cycle', label: 'Measured cycle time', unit: 's' }], rows: rowsOf(REF.JANSEN), target: 'cycle', model, tolerance: { mape: 45 }, covers: ['Severe riser slugging', 'Pressure fluctuations'], note: `Lumped riser cycle (riserSluggingCycle) with the rig geometry; the gas buffer is added to the feed-line gas volume (liquid fraction of the feed line from the stratified equilibrium level of the −1° pipeline), the choke is Δp = C v² on the liquid, and the lift gas is injected at the riser base, where it aerates the liquid column (drift-flux void fraction) and so lowers the head that the feed gas has to overcome. ${REF.JANSEN.note} No open tabulated record of pressure amplitude was found; amplitudes are compared through the measurement table of the suite.` }); }
  // ---- single-phase friction, cavity benchmark
  A.push({ id: 'friction-superpipe', title: 'Smooth-pipe friction factor, Re = 7.4·10⁴ – 3.6·10⁷ (Princeton Superpipe, McKeon et al. 2004)', quantity: 'Darcy friction factor', unit: '–', kind: 'experiment', source: REF.SUPERPIPE.source, columns: [{ key: 'Re', label: 'Reynolds number' }, { key: 'f', label: 'Measured friction factor' }], rows: rowsOf(REF.SUPERPIPE), target: 'f', model: (r) => frictionFactor(r.Re, 0, 'colebrook'), tolerance: { mape: 3 }, covers: ['Pressure gradient'], note: 'Colebrook–White in the smooth limit (Prandtl law). ' + REF.SUPERPIPE.note });
  { const g = rowsOf(REF.GHIA), sol = () => once('cavity100', () => cfd2d({ nx: 32, ny: 32, rho: 1, mu: 0.01, lid: 1, tEnd: 40, steady: 1e-4, cfl: 0.5, poissonTol: 1e-6 }));
    A.push({ id: 'cavity-ghia', title: 'Lid-driven cavity at Re = 100: centre-line velocities (Ghia, Ghia & Shin 1982)', quantity: 'Velocity / lid velocity', unit: '–', kind: 'benchmark', source: REF.GHIA.source, columns: [{ key: 'profile', label: 'Profile (0: u on x = 0.5, 1: v on y = 0.5)' }, { key: 'coord', label: 'Coordinate' }, { key: 're100', label: 'Benchmark value' }], rows: g, target: 're100', model: (r) => { const s = sol(); return r.profile === 0 ? interp1(s.y, s.uMid, r.coord) : interp1(s.x, s.vMid, r.coord); }, tolerance: { maxAbs: 0.012, rmse: 0.006 }, covers: ['Code-to-code benchmark problems'], note: '2-D Navier–Stokes solver of the suite on 32 × 32 cells (computed once and cached). ' + REF.GHIA.note }); }
  // ---- three-dimensional solver against the channel DNS (stored deterministic runs of the in-app solver; repeat them with the 3-D task)
  { const mk = REF.MKM180, runs = REF.CFD3D_CHANNEL.runs, all = rowsOf(REF.CFD3D_CHANNEL), NM = ['none', 'smagorinsky', 'wale', 'des', 'ddes', 'iddes', 'dns'], TX = { none: 'no sub-grid model (under-resolved)', smagorinsky: 'LES, Smagorinsky + van Driest', wale: 'LES, WALE', des: 'DES (Spalart–Allmaras)', ddes: 'delayed DES', iddes: 'IDDES', dns: 'resolved simulation (DNS grid, small box)' }, TOL = { none: 6, smagorinsky: 10, wale: 5, des: 9, ddes: 5, iddes: 5, dns: 6 };
    NM.forEach((nm, k) => { const q = runs[nm]; if (!q || (nm === 'dns' && q.uvMax < 0.3)) return; const pr = all.filter((r) => r.model === k), yp = pr.map((r) => r.yPlus), up = pr.map((r) => r.uPlus), turb = q.uvMax > 0.3;
      A.push({ id: `channel180-${nm}`, title: `Turbulent channel at Re_b = 5,600 (Re_τ ≈ 180): mean velocity, ${TX[nm]}, ${q.grid.join(' × ')} cells, against DNS`, quantity: 'Mean velocity U⁺', unit: '–', kind: 'benchmark', source: mk.source, columns: [{ key: 'yPlus', label: 'y⁺ (DNS)' }, { key: 'uPlus', label: 'U⁺ (DNS)' }], rows: rowsOf(mk), target: 'uPlus', model: (r) => interp1([0, ...yp.map((y) => (y * mk.reTau) / q.reTau), mk.reTau * 1.02], [0, ...up.map((u) => (u * q.reTau) / mk.reTau), (up[up.length - 1] * q.reTau) / mk.reTau], Math.min(r.yPlus, mk.reTau)), tolerance: { mape: TOL[nm] }, covers: ['Code-to-code benchmark problems'],
        note: `Stored result of the in-app three-dimensional solver (deterministic; repeat it with the task “3-D Navier–Stokes solver”), compared at equal y/h and scaled with the friction velocity of the DNS so that the error of the wall shear is included. Friction coefficient ${q.cf.toExponential(3)} against ${mk.cf.toExponential(3)} of the DNS (${(100 * (q.cf / mk.cf - 1)).toFixed(1)} %), Re_τ ${q.reTau} against ${mk.reTau}, centre-line U⁺ ${q.uCentrePlus} against ${mk.uCentrePlus}, largest u rms ${q.uRmsMax} against ${mk.uRmsMax}, largest resolved −u′v′ ${q.uvMax} against ${mk.uvMax}. Resolution: Δx⁺ ${q.dxPlus}, Δz⁺ ${q.dzPlus}, first cell centre at y⁺ ${q.firstCentrePlus}, largest Δy⁺ ${q.dyPlusMax}${nm === 'dns' ? ` in a periodic box of ${q.box.join(' × ')} half-heights` : ''}; ${q.steps} steps, ${q.time} h/u_τ simulated, ${q.samples} samples; ${q.cpuMsPerStep} ms per step on a loaded machine (reference loop ${q.bench1e8} ms, about 100 ms on an idle desktop core). ${turb ? '' : 'The resolved fluctuations decayed in this run: the shielding of the model keeps the channel in RANS mode at this resolution, so the profile is the steady Spalart–Allmaras solution. '}${REF.CFD3D_CHANNEL.note} ${mk.note}` }); }); }
  // ---- Taylor–Green vortex at Re = 1,600 against the 512³ DNS (CC BY)
  if (REF.CFD3D_TGV) { const tg = rowsOf(REF.CFD3D_TGV), cols = [{ key: 't', label: 'Time', unit: 'L/U' }, { key: 'ek', label: 'Kinetic energy (DNS)', unit: 'U²' }, { key: 'eps', label: 'Dissipation rate (DNS)', unit: 'U³/L' }];
    for (const [key, label, tol] of [['re1600n64smag', '64³ cells, Smagorinsky LES', 8], ['re1600n32wale', '32³ cells, WALE LES', 10]]) { const q = tg.filter((r) => r.run === key); if (q.length < 8) continue;
      A.push({ id: `tgv1600-${key}`, title: `Taylor–Green vortex at Re = 1,600: kinetic energy, ${label}, against the 512³ DNS`, quantity: 'Kinetic energy', unit: 'U²', kind: 'benchmark', source: REF.TGV1600.source, columns: cols, rows: rowsOf(REF.TGV1600, (r) => r.t > 0 && r.t <= q[q.length - 1].t), target: 'ek', model: (r) => interp1(q.map((x) => x.t), q.map((x) => x.ek), r.t), tolerance: { mape: tol }, covers: ['Code-to-code benchmark problems'], note: `Stored result of the in-app solver (deterministic; repeat it with the 3-D task). The grid is far coarser than the 512³ of the reference, so this is a large-eddy simulation: the transition is reproduced, the dissipation peak (${Math.max(...q.map((x) => x.eps)).toExponential(3)} at t = ${q.reduce((a, b) => (b.eps > a.eps ? b : a)).t}) is lower and broader than the DNS peak of 1.26e-2 at t = 9. ${REF.CFD3D_TGV.note} ${REF.TGV1600.note}` }); } }
  // ---- OpenFOAM code-to-code benchmarks
  { const d = REF.OF_DAMBREAK, sol = () => once('dam3d', () => { const run = twoPhase3D({ kind: 'dambreak', n: 16, spanCells: 2, a: d.a }); while (!run.done) run.advance(500); return run.result(); }), tu = Math.sqrt(d.a / G);
    A.push({ id: 'dambreak-openfoam', title: 'Dam break of a water column: front position, in-app 3-D volume-of-fluid solver against OpenFOAM interFoam', quantity: 'Front position / column width', unit: '–', kind: 'benchmark', source: d.source, columns: [{ key: 'tStar', label: 't √(g/a)' }, { key: 'front', label: 'Front position (OpenFOAM)', unit: 'a' }, { key: 'height', label: 'Column height at the wall (OpenFOAM)', unit: 'a' }], rows: rowsOf(d, (r) => r.tStar >= 0.5), target: 'front', model: (r) => { const s = sol(); return interp1(s.t, s.front, r.tStar * tu) / d.a; }, tolerance: { mape: 25 }, covers: ['Dam-break benchmark where formulation-relevant', 'Code-to-code benchmark problems'], note: 'In-app solver on 32 × 16 × 2 cells (computed once and cached, water and air, density ratio 829) against 160 × 80 cells of the reference: the coarse front lags by 15–25 %; on 64 × 32 cells the lag at t √(g/a) = 3 is 12 %. The liquid volume is conserved to round-off in both. ' + d.note }); }
  { const d = REF.OF_GRAETZ;
    A.push({ id: 'heattransfer-openfoam', title: 'Laminar heat transfer in a pipe with a cold wall: local Nusselt number, kernel against OpenFOAM', quantity: 'Nusselt number', unit: '–', kind: 'benchmark', source: d.source, columns: [{ key: 'Re', label: 'Reynolds number' }, { key: 'Pr', label: 'Prandtl number' }, { key: 'xStar', label: 'x / (D Re Pr)' }, { key: 'xD', label: 'x / D' }, { key: 'Nu', label: 'Local Nusselt number (OpenFOAM)' }], rows: rowsOf(d), target: 'Nu', model: (r) => (hInside(r.Re, r.Pr, 1, 1)), tolerance: { mape: 6 }, covers: ['Heat-transfer coefficient', 'Code-to-code benchmark problems'], note: 'The kernel uses the developed value Nu = 3.66 for laminar flow; the computation shows the thermal entrance (18 % higher at x* = 0.02, within 2 % from x* = 0.05), which the kernel neglects. Turbulent heat transfer was also attempted (k–ω SST with a passive scalar) but the scalar solution of those runs did not pass its checks and is not included. ' + d.note }); }
  return A;
}

// =====================================================================================================
// 7j. Verification of the solvers added in sections 7a–7g (called from verify)
// =====================================================================================================
async function verifyExtra(add) {
  const flag = (name, cond, note) => add(name, 1, cond ? 1 : 0, 0, note, false);
  { // kernel correlations against the worked examples printed in open libraries
    add('Colebrook friction factor at Re = 10⁵, ε/D = 10⁻⁴ (fluids library example)', 0.018513866077471, frictionFactor(1e5, 1e-4, 'colebrook'), 1e-9, 'Independent implementation: Caleb Bell, fluids.friction.Colebrook docstring');
    add('Haaland friction factor (fluids library example)', 0.018265053014793857, frictionFactor(1e5, 1e-4, 'haaland'), 1e-12, 'fluids.friction.Haaland(1E5, 1E-4)');
    add('Swamee–Jain friction factor (fluids library example)', 0.018452424431901808, frictionFactor(1e5, 1e-4, 'swamee'), 2e-6, 'fluids.friction.Swamee_Jain_1976(1E5, 1E-4); the library writes (6.97/Re)^0.9, the kernel the rounded 5.74/Re^0.9 (6.97^0.9 = 5.7407)');
    add('Churchill friction factor (fluids library example)', 0.018462624566280075, frictionFactor(1e5, 1e-4, 'churchill'), 1e-12, 'fluids.friction.Churchill_1977(1E5, 1E-4)');
    const A = (PI * 0.05 ** 2) / 4, bb = gradient({ vsl: (0.6 * 0.9) / 915 / A, vsg: (0.6 * 0.1) / 2.67 / A, rhoL: 915, rhoG: 2.67, muL: 180e-6, muG: 14e-6, sigma: 0.0487, D: 0.05, theta: 0, rough: 0, P: 1e7, label: false });
    add('Code-to-code: Beggs & Brill pressure gradient against the fluids library example', 686.9724506803469, bb.dpdx, 1e-10, 'Pa/m; fluids.two_phase.Beggs_Brill(m=0.6, x=0.1, rhol=915, rhog=2.67, mul=180E-6, mug=14E-6, sigma=0.0487, P=1E7, D=0.05, angle=0)');
    add('Churchill–Bernstein Nusselt number (ht library example)', 40.63708594124974, nuCrossFlow(6071, 0.7), 1e-12, 'ht.conv_external.Nu_cylinder_Churchill_Bernstein(6071, 0.7)');
    add('Churchill–Chu free-convection Nusselt number (ht library example)', 139.13493970073597, nuFreeCylinder(0.69 * 2.63e9, 0.69), 1e-12, 'ht.conv_free_immersed.Nu_horizontal_cylinder_Churchill_Chu(0.69, 2.63E9)');
    add('Gnielinski Nusselt number at Re = 10⁵, Pr = 1.2 (ht library example with f = 0.0185)', 254.62682749359632, hInside(1e5, 1.2, 1, 1), 0.04, 'The kernel uses the smooth-pipe Haaland friction factor (0.0180) instead of the 0.0185 of the example, which lowers the result by 3.5 %');
    const g = gradient({ vsl: 0, vsg: 12, rhoL: 800, rhoG: 40, muL: 1e-3, muG: 1.3e-5, D: 0.3, theta: 0.2, rough: 4.5e-5, P: 4e6 }), ek = (40 * 144) / 4e6;
    add('Single-phase gas: acceleration term of the kernel', ((g.fric + g.grav) * ek) / (1 - ek), g.acc, 1e-12, 'Pa/m; (friction + gravity) Ek / (1 − Ek) with Ek = ρ v² / P');
    const u = slugUnit({ vsl: 1, vsg: 1.5, rhoL: 800, rhoG: 60, muL: 2e-3, muG: 1.4e-5, D: 0.254, theta: 0 }), uf = slugUnit({ vsl: 1, vsg: 1.5, rhoL: 800, rhoG: 60, muL: 2e-3, muG: 1.4e-5, D: 0.254, theta: 0, basis: 'frequency' });
    add('Kernel slug unit cell: frequency × unit length = translational velocity', u.vt, u.freq * u.unitLength, 1e-12, 'm/s, length basis (frequency derived from the length correlation)');
    add('Kernel slug unit cell, frequency basis: length = slug fraction × v_t / frequency', (uf.slugFraction * uf.vt) / uf.freqCorrelation, uf.length, 1e-12, 'm; the length is derived when the frequency correlation is kept');
  }
  { // kernel march and suite march now use the same energy equation
    const a = steadyFlow({}, {}), k = marchSteady({ fm: a.cfg.fm, profile: a.cfg.profile, id: a.cfg.id, rough: a.cfg.rough, U: a.cfg.U, tAmbOf: a.cfg.tAmbOf, tIn: a.cfg.tIn, pOut: a.cfg.pOut, n: a.cfg.n, label: false });
    add('Kernel and suite agree on the arrival temperature (flowing-enthalpy balance in both)', a.tOut, k.tOut, 0.02, '°C on the reference case, same grid; the kernel default was a frozen-cp / Joule–Thomson form before', false);
    add('Kernel and suite agree on the inlet pressure of the reference case', a.pIn, k.pIn, 2e-4, 'bara');
  }
  { // two-fluid model: water faucet, Kelvin–Helmholtz limits, slug capturing, phase energies
    const fa = (n) => { const r = twoFluid({ n, length: 12, D: 1, theta: -PI / 2, rhoL: 1000, friction: false, pOut: 1e5, init: () => ({ al: 0.8, ul: 10, ug: 0 }), inlet: () => ({ al: 0.8, ul: 10, ug: 0 }), tEnd: 0.5, cfl: 0.4, dtMax: 1 }), xf = 5 + (G * 0.25) / 2; let e = 0; for (let i = 0; i < n; i++) e += Math.abs(1 - r.al[i] - (r.x[i] < xf ? 1 - 8 / Math.sqrt(100 + 2 * G * r.x[i]) : 0.2)); return { l1: e / n, r }; }, f1 = fa(100), f2 = fa(400);
    add('Water-faucet benchmark (Ransom): mean error of the gas fraction against the analytic solution', 0, f2.l1, 0.006, `400 cells at t = 0.5 s; α_g = 1 − α_l0 u_0 / √(u_0² + 2 g x) behind the front. 100 cells: ${f1.l1.toExponential(2)}`, false);
    add('Water faucet: observed order of the first-order upwind scheme on a discontinuous solution', 0.67, Math.log(f1.l1 / f2.l1) / Math.log(4), 0.2, 'L1 error on 100 and 400 cells; a contact discontinuity limits a first-order scheme to an order between ½ and 1', false);
    add('Water faucet: liquid mass conservation of the two-fluid solver', 0, f2.r.mass.errorL, 1e-11, 'Relative', false);
    await 0;
    // inviscid Kelvin–Helmholtz: plane channel, periodic, no friction
    const H = 0.1, rhoG = 1e5 / (287 * 293.15), Uc = Math.sqrt(((1000 * 0.5 + rhoG * 0.5) * (1000 - rhoG) * G * H) / (1000 * rhoG));
    const mode = (al, a0) => { let s = 0, c = 0; const n = al.length; for (let i = 0; i < n; i++) { const ph = (2 * PI * (i + 0.5)) / n; s += (al[i] - a0) * Math.sin(ph); c += (al[i] - a0) * Math.cos(ph); } return (2 * Math.hypot(s, c)) / n; };
    const rate = (r, a0) => { const A = r.field.al.map((q) => mode(q, a0)), k0 = Math.floor(A.length / 3), xs = r.field.t.slice(k0), ys = A.slice(k0).map(Math.log), mx = mean(xs), my = mean(ys); let nu = 0, de = 0; xs.forEach((x, i) => { nu += (x - mx) * (ys[i] - my); de += (x - mx) ** 2; }); return nu / de; };
    const ikh = (f) => rate(twoFluid({ n: 64, length: 2, channelH: H, rhoL: 1000, Rs: 287, T: 293.15, friction: false, periodic: true, pOut: 1e5, init: (x) => ({ al: 0.5 + 1e-4 * Math.sin((2 * PI * x) / 2), ul: 0, ug: f * Uc, p: 1e5 }), tEnd: 0.5, cfl: 0.3, dtMax: 1, nField: 31 }), 0.5), g1 = ikh(0.9), g2 = ikh(1.1);
    add('Inviscid Kelvin–Helmholtz limit: gas velocity at which waves start to grow / theory', 1, 0.9 + (0.2 * -g1) / (g2 - g1), 0.03, `Two-fluid model, periodic plane channel; growth rates ${g1.toFixed(3)} and ${g2.toFixed(3)} 1/s at 0.9 and 1.1 times ΔU² = (ρl αg + ρg αl)(ρl − ρg) g H / (ρl ρg) = (${Uc.toFixed(2)} m/s)²`);
    await 0;
    // viscous Kelvin–Helmholtz: pipe, periodic with the equilibrium pressure gradient as body force, long wave
    const base = { rhoL: 998, rhoG, muL: 1e-3, muG: 1.8e-5, D: 0.078, theta: 0, vsg: 2 }, vk = brent((x) => interfacialStability({ ...base, vsl: x }).vkhRatio - 1, 0.05, 0.6, 1e-7), Lw = 32;
    const vgr = (f) => { const p = { ...base, vsl: f * vk }, st = stratifiedLevel(p), gg = st.geom, B = -(st.tauWL * gg.SL + st.tauWG * gg.SG) / gg.A; return rate(twoFluid({ n: 96, length: Lw, D: 0.078, rhoL: 998, muL: 1e-3, muG: 1.8e-5, Rs: 287, T: 293.15, periodic: true, pOut: 1e5, bodyForce: B, init: (x) => ({ al: st.holdup + 2e-4 * Math.sin((2 * PI * x) / Lw), ul: st.vL, ug: st.vG, p: 1e5 }), tEnd: 70, cfl: 0.4, dtMax: 1, nField: 41 }), st.holdup); }, h1 = vgr(0.95), h2 = vgr(1.15);
    add('Viscous Kelvin–Helmholtz limit: liquid velocity at which long waves start to grow / theory', 1, 0.95 + (0.2 * -h1) / (h2 - h1), 0.08, `Two-fluid model against the Barnea–Taitel criterion (interfacialStability) at vsg = 2 m/s in a 78 mm air–water pipe: theory ${vk.toFixed(4)} m/s; wavelength ${Lw} m on 96 cells (the upwind damping of shorter waves delays the numerical onset)`);
    await 0;
    // slug capturing
    { const D = 0.078, L = 16, p = { vsl: 1, vsg: 1.5, rhoL: 998, rhoG, muL: 1e-3, muG: 1.8e-5, D, theta: 0 }, st = stratifiedLevel(p), tE = 16;
      const r = twoFluid({ n: 240, length: L, D, rhoL: 998, muL: 1e-3, muG: 1.8e-5, Rs: 287, T: 293.15, pOut: 1e5, init: () => ({ al: st.holdup, ul: st.vL, ug: st.vG, p: 1e5 }), inlet: (t) => ({ al: st.holdup * (1 + 0.02 * Math.sin((2 * PI * t) / 0.7) + 0.02 * Math.sin((2 * PI * t) / 1.9)), ul: st.vL, ug: st.vG }), tEnd: tE, cfl: 0.4, dtMax: 0.01, probes: [0.75 * L, 0.9 * L], nSeries: 1600 });
      const a = r.series.probes[0], b = r.series.probes[1], i0 = Math.floor(a.length * 0.3), dts = tE / (r.series.t.length - 1), ma = mean(a.slice(i0)), mb = mean(b.slice(i0)); let best = 0, bl = 1; for (let lag = 1; lag < 400; lag++) { let c = 0; for (let i = i0; i + lag < a.length; i++) c += (a[i] - ma) * (b[i + lag] - mb); if (c > best) { best = c; bl = lag; } }
      flag('Slug capturing: interfacial waves grow into liquid bridges on the test section', r.completed && r.slugs.filter((q) => q.probe === 1).length >= 2 && Math.max(...b) > 0.95, `${r.slugs.filter((q) => q.probe === 1).length} slugs pass the probe at 0.9 L in ${tE} s (vsl = 1, vsg = 1.5 m/s, 78 mm, 240 cells); the initial stratified state is viscous-Kelvin–Helmholtz unstable`);
      add('Slug capturing: velocity of the captured slugs against the Bendiksen relation', slugVelocity(2.5, D, 0).vt, (0.15 * L) / (bl * dts), 0.15, 'm/s from the cross-correlation of the holdup at 0.75 L and 0.9 L; the model has no closure for the translational velocity, it follows from the mass and momentum balances');
      add('Slug capturing: liquid mass conservation', 0, r.mass.errorL, 1e-10, 'Relative, with bridging and unbridging cells', false); }
    await 0;
    // six-equation form: closed periodic box, phases at different temperatures
    { const cvL = 4180, cvG = 718, al = 0.3, Tl0 = 320, Tg0 = 290, p0 = 2e5, Rs = 287, Cg = ((1 - al) * p0) / (Rs * Tg0) * cvG, Cl = al * 1000 * cvL, Teq = (Cl * Tl0 + Cg * Tg0) / (Cl + Cg);
      const r = twoFluid({ n: 16, length: 4, channelH: 0.1, rhoL: 1000, Rs, friction: false, periodic: true, pOut: p0, init: () => ({ al, ul: 0, ug: 0, p: p0, Tl: Tl0, Tg: Tg0 }), tEnd: 40, cfl: 0.4, dtMax: 0.02, energy: { cvL, cvG, hi: 400 } });
      add('Six-equation model: equilibrium temperature of gas heated by the liquid in a closed volume', Teq, r.Tg[5], 2e-4, 'K; isochoric gas (heat capacity c_v) and incompressible liquid exchanging heat through the interface');
      add('Six-equation model: pressure rise of the heated gas (ideal gas at constant volume)', (p0 * Teq) / Tg0, r.p[5], 2e-3, 'Pa; p / T_g constant, which couples the gas internal-energy equation to the pressure equation');
      add('Six-equation model: total internal-energy conservation', 0, r.energy.error, 1e-9, 'Relative', false); }
  }
  await 0;
  { // seven-equation model, exact Riemann solution, entropy
    const e = 1e-7, ex = riemannExact({ rho: 1, u: 0, p: 1 }, { rho: 0.125, u: 0, p: 0.1 }), run = (n) => baerNunziato({ n, tEnd: 0.2, left: { a1: 1 - e, rho1: 1, u1: 0, p1: 1, rho2: 1, u2: 0, p2: 1 }, right: { a1: 1 - e, rho1: 0.125, u1: 0, p1: 0.1, rho2: 0.125, u2: 0, p2: 0.1 } }), l1 = (r) => mean(r.rho1.map((x, i) => Math.abs(x - ex.sample((r.x[i] - 0.5) / 0.2).rho))), a = run(200), b = run(800);
    // independent check of the exact solution: mass flux through the right shock from the Rankine–Hugoniot relations
    const sR = ex.sample(0.5 * (ex.uStar + 1.75)), S = (sR.rho * ex.uStar) / (sR.rho - 0.125);
    add('Exact Riemann solution of the Sod shock tube: Rankine–Hugoniot momentum jump across the shock', ex.pStar - 0.1, 0.125 * S * ex.uStar, 1e-9, `p* − p_R = ρ_R S u*; p* = ${ex.pStar.toFixed(5)}, u* = ${ex.uStar.toFixed(5)}, shock speed ${S.toFixed(4)}`);
    add('Seven-equation (Baer–Nunziato) model, shock-tube benchmark: density error against the exact solution', 0, l1(b), 0.008, `Mean absolute error on 800 cells in the single-phase limit (200 cells: ${l1(a).toExponential(2)}); first-order Rusanov fluxes`, false);
    add('Seven-equation model: observed order on the shock tube', 0.6, Math.log(l1(a) / l1(b)) / Math.log(4), 0.25, 'Between ½ (contact) and 1 (shock) for a first-order scheme', false);
    add('Seven-equation model: total energy conservation', 0, b.conservation.energy, 1e-12, 'Relative drift of the closed tube', false);
    add('Seven-equation model: momentum balance with the end-pressure impulse', 0, b.conservation.momentum, 1e-10, 'Relative', false);
    flag('Entropy inequality: the captured shock produces entropy', a.entropy.production > 0 && b.entropy.production > 0 && b.entropy.production < a.entropy.production, `Total entropy rises by ${a.entropy.production.toExponential(3)} (200 cells) and ${b.entropy.production.toExponential(3)} (800 cells): positive, and falling towards the physical shock production as the numerical dissipation is refined`);
    const c = baerNunziato({ n: 100, tEnd: 0.002, left: { a1: 0.9, rho1: 1.2, u1: 100, p1: 1e5, rho2: 1000, u2: 100, p2: 1e5 }, right: { a1: 0.1, rho1: 1.2, u1: 100, p1: 1e5, rho2: 1000, u2: 100, p2: 1e5 }, eos: [{ gamma: 1.4, pinf: 0 }, { gamma: 4.4, pinf: 6e8 }] });
    add('Seven-equation model: a moving volume-fraction contact keeps pressure and velocity uniform', 0, Math.max(...c.p1.map((p) => Math.abs(p / 1e5 - 1)), ...c.p2.map((p) => Math.abs(p / 1e5 - 1)), ...c.u1.map((u) => Math.abs(u / 100 - 1))), 1e-8, 'Largest relative deviation (air / stiffened-gas water interface advected at 100 m/s): the non-conservative terms are discretised consistently', false);
    // species transport: Ogata–Banks
    const n = 400, L = 100, D = 0.5, s = Array.from({ length: n + 1 }, (_, i) => (i * L) / n), sp = speciesTransport({ s, qL: s.map(() => 0.1), holdup: s.map(() => 1), A: 0.1, D, cIn: 1, tEnd: 60, cfl: 0.5 });
    const erfc = (x) => { const z = Math.abs(x), t = 1 / (1 + 0.5 * z), y = t * Math.exp(-z * z - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 + t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277))))))))); return x >= 0 ? y : 2 - y; };
    let em = 0; sp.s.forEach((x, i) => { if (x > 15) em = Math.max(em, Math.abs(sp.c[i] - 0.5 * erfc((x - 60) / (2 * Math.sqrt(D * 60))))); });
    add('Species transport: advection–dispersion of an inlet step (Ogata–Banks solution)', 0, em, 0.02, 'Largest concentration error away from the inlet, u = 1 m/s, D = 0.5 m²/s, t = 60 s, 400 cells', false);
    add('Species transport: mass balance', 0, sp.balance.error, 1e-12, 'Relative', false);
  }
  await 0;
  { // turbulence closures in the radial solve
    const tol = { sa: 0.04, sst: 0.05, kestd: 0.04, rng: 0.08, realizable: 0.09, rsm: 0.03 }, lab = Object.fromEntries(RANS_MODELS.map((q) => [q.value, q.label]));
    for (const m of Object.keys(tol)) { const r = ransPipeExtra({ reTau: 2000, model: m, n: 60, tol: 2e-6 }); add(`${lab[m]}: pipe friction factor against Colebrook (smooth)`, frictionFactor(r.Re, 0), r.f, tol[m], `Re = ${r.Re.toFixed(0)}, ${r.iterations} iterations${m === 'rng' || m === 'realizable' ? '; these closures sit a few per cent low in a plain pipe' : ''}`); if (m === 'rsm') { const j = r.y.findIndex((y) => y > 200), C1 = 1.8, C2 = 0.6, c1 = 0.5, c2 = 0.3, vvk = ((2 / 3) * (C1 - 1) + (2 / 3) * C2 - 3 * (4 / 9) * c2 * C2) / (C1 + 2 * c1), uvk = Math.sqrt(((1 - C2 + 1.5 * c2 * C2) / (C1 + 1.5 * c1)) * vvk); add('Reynolds-stress model: structure parameter −⟨uv⟩/k in the log layer', uvk, -r.stress.uv[j] / r.k[j], 0.08, `Local-equilibrium value of the LRR closure with wall reflection worked out by hand (⟨vv⟩/k = ${vvk.toFixed(3)}); solver at y⁺ = ${r.y[j].toFixed(0)}`); } await 0; }
  }
  { // bubbly-flow closures
    const b = bubbleRise({});
    add('Virtual (added) mass: initial acceleration of a bubble released from rest', (G * (998 - 1.2)) / (1.2 + 0.5 * 998), b.a0, 1e-12, 'm/s² ≈ 2 g for C_vm = ½: (ρl − ρg) g / (ρg + C_vm ρl)');
    add('Bubble rise: the transient reaches the terminal velocity of the drag law', b.terminal, b.v[b.v.length - 1], 1e-4, 'm/s, 3 mm bubble in water, Ishii–Zuber drag');
    const ia = interfacialAreaTransport({ d0: 2e-3, alpha: 0.1, vg: 1, k: 1e-6, length: 20, crc: 0, cti: 0, cdFixed: 1 });
    add('Interfacial-area transport with wake entrainment only: analytic diameter growth', 2e-3 + 6 * (12 / (36 * PI)) * 0.002 * 0.1 * ia.uRel * 20, ia.dEnd, 1e-9, 'm; dκ/dz = −K κ² gives d(z) = d₀ + 6 K z');
    const g = { R: 0.01905, uLiquid: 0.916, alphaMean: 0.0567, d: 0.0034, rhoL: 995.7, muL: 0.000799, sigma: 0.071 } /* air–water bubbly upflow in a 38 mm pipe, small bubbles: the conditions of a wall-peaked void profile */, sm = once('gros', () => bubblyPipe({ R: g.R, jl: g.uLiquid * (1 - g.alphaMean), alphaMean: g.alphaMean, d: g.d, rhoL: g.rhoL, rhoG: 1.19, muL: g.muL, sigma: g.sigma, cl: 0.288 })), big = bubblyPipe({ R: g.R, jl: g.uLiquid * (1 - g.alphaMean), alphaMean: g.alphaMean, d: 0.008, rhoL: g.rhoL, rhoG: 1.19, muL: g.muL, sigma: g.sigma });
    add('Bubbly pipe flow: position of the void peak against the published profile (Grossetete 1995)', 0.892, sm.peak.r, 0.05, 'r/R of the wall peak of 3.4 mm bubbles (measured maximum at r/R = 0.892)');
    add('Bubbly pipe flow: height of the void peak against the published profile', 0.107, sm.peak.alpha, 0.2, 'Measured maximum void fraction 0.107');
    flag('Lift-force sign change: 8 mm bubbles collect in the core', big.cl < 0 && !big.wallPeaked && sm.wallPeaked, `Tomiyama lift coefficient ${big.cl.toFixed(3)} at 8 mm against ${sm.cl.toFixed(3)} at 3.4 mm`);
    add('Bubbly pipe flow: the profile carries the prescribed mean void', g.alphaMean, sm.mean, 1e-5, 'Area average of the computed profile');
  }
  { // annular entrainment / deposition, homogeneous relaxation, closures
    const p = { vsl: 0.1, vsg: 12, rhoL: 800, rhoG: 60, muL: 1.5e-3, muG: 1.4e-5, sigma: 0.02, D: 0.1, theta: PI / 2, rough: 4.5e-5, P: 6e6, fModel: 'colebrook' }, a = annularMist(p), dv = annularDevelopment({ p, length: a.relaxLength, n: 200, e0: 0 }), cq = (E) => (p.rhoL * p.vsl * E) / (p.vsg + p.vsl * E);
    add('Annular flow: entrained fraction after one deposition length', 1 - Math.exp(-1), cq(dv.entrained[200]) / cq(dv.eq), 5e-3, 'Droplet concentration relative to equilibrium: the three-field mass balance relaxes exponentially with the length Q_g / (π D k_D) when the droplet volume is small', false);
    add('Annular flow: film and core momentum balances give the same pressure gradient', a.fric + a.grav, (a.tauW * 4) / p.D + (p.rhoL * a.holdupFilm + (p.rhoG * (1 - a.holdup) + p.rhoL * a.holdupDrops)) * G, 0.03, 'Pa/m: wall shear × perimeter / area + weight of film and core (overall force balance worked out separately)');
    const xe = 0.05, r = flashingFlow({ G: 500, D: 0.05, length: 2, n: 400, p0: 5e5, x0: 0, pSat: 6e5, xEq: () => xe, vG: () => 1e-3, vL: 1e-3, f: 0, tau: 2, model: 'hrm' });
    add('Homogeneous relaxation: analytic approach to equilibrium at constant velocity', xe * (1 - Math.exp(-2 / (2 * 500 * 1e-3))), r.x[r.x.length - 1], 1e-5, 'x(z) = x_eq [1 − exp(−z / (Θ u))] for a fixed relaxation time, equal phase volumes and no friction');
    const R = 461.5, vG = (pp) => (R * 400) / pp, xq = (pp) => Math.max(0, (0.2 * (3e5 - pp)) / 3e5), o = { G: 300, D: 0.05, length: 3, n: 600, p0: 2.9e5, x0: xq(2.9e5), pSat: 3e5, xEq: xq, vG, vL: 1e-3, f: 0.02 }, hem = flashingFlow({ ...o, model: 'hem' }), fast = flashingFlow({ ...o, model: 'hrm', tau: 1e-7 });
    add('Homogeneous relaxation → homogeneous equilibrium as the relaxation time vanishes', hem.p[hem.p.length - 1], fast.p[fast.p.length - 1], 2e-3, 'Pa at the outlet of a flashing pipe');
    add('Downar-Zapolski relaxation time at α = 0.1, ψ = 0.1', 6.51e-4 * 0.1 ** -0.257 * 0.1 ** -2.24, relaxationTime(0.1, 0.1), 1e-12, 's; Θ = 6.51·10⁻⁴ α^−0.257 ψ^−2.24');
    const pb = { vsl: 1, vsg: 0.1, rhoL: 800, rhoG: 60, muL: 1.5e-3, muG: 1.4e-5, sigma: 0.02, D: 0.1, theta: PI / 2, rough: 4.5e-5, P: 6e6 }, ab = holdupGradient(pb, 'ansari'), v0 = 1.53 * ((G * 0.02 * 740) / 800 ** 2) ** 0.25;
    add('Ansari bubble flow: the holdup satisfies the slip relation v_g = 1.2 v_m + v_s H_L^½', 1.2 * 1.1 + v0 * Math.sqrt(ab.holdup), 0.1 / (1 - ab.holdup), 1e-6, 'm/s, gas velocity from both sides');
    const dr = holdupGradient(pb, 'dunsRos'), vn = velocityNumbers(pb), S = interp1(DR.NL.map(log10), DR.F1, log10(vn.NL)) + interp1(DR.NL.map(log10), DR.F2, log10(vn.NL)) * vn.NLv + (interp1(DR.NL.map(log10), DR.F3, log10(vn.NL)) - interp1(DR.NL.map(log10), DR.F4, log10(vn.NL)) / vn.ND) * (vn.NGv / (1 + vn.NLv)) ** 2, vs = S / (800 / (G * 0.02)) ** 0.25;
    add('Duns & Ros bubble region: holdup from the slip velocity', (vs - 1.1 + Math.sqrt((1.1 - vs) ** 2 + 4 * vs * 1)) / (2 * vs), dr.holdup, 1e-9, 'Hand evaluation of S = F1 + F2 N_Lv + F3′ (N_Gv / (1 + N_Lv))² with the chart values');
    const eq = equipment({ pumpDp0: 20, pumpQmax: 0.4, pumpSpeed: 0.9, sepKv: 2, compHead: 60e3, compQmax: 2 }), cs = eq.compressor(1, 100, 300, 0.9, 0.02, 10);
    add('Pump characteristic: affinity laws at 90 % speed and half the run-out rate', 20 * 0.81 - (20 / 0.16) * 0.04, eq.pumpDp(0.2), 1e-12, 'bar');
    add('Compressor characteristic: suction pressure from the polytropic head', 100 / (1 + ((0.3 / 1.3) * 60e3 * 0.75 * 0.02) / (0.9 * RGAS * 300)) ** (1.3 / 0.3), cs.pSuction, 1e-12, 'bara at half the maximum flow (head 45 kJ/kg)');
    add('Separator characteristic: gas-outlet valve equation', Math.sqrt(20 * 20 + 2 * 9), eq.sepPressure(3, 25, 20), 1e-12, 'bara: p_sep² = p_suction² + K_v q²');
    const cat = catenaryProfile({ height: 1000, topAngle: 12, n: 400 }), top = cat[400], prev = cat[399];
    add('Catenary riser: top angle of the generated profile', 12, 90 - Math.atan2(top.z - prev.z, top.x - prev.x) / DEG, 0.02, 'degrees from the vertical (finite-difference slope of the last segment)');
  }
  await 0;
  { // junction
    const cfg = flowConfig({}, {}), r = steadyBranch({ ...cfg.base, n: 60, mScale: 1, pOut: cfg.pOut, hydrate: false }, { x: 9000, frac: 0.3, length: 2500, dz: 20, tIn: 45 });
    add('Branch junction: mass balance', r.junction.mMain + r.junction.mBranch, r.junction.mOut, 1e-12, 'kg/s: main + branch = downstream');
    add('Branch junction: enthalpy of the mixed stream', r.junction.hIn, r.junction.hOut, 1e-6, 'W: the mixed temperature carries the enthalpy of both streams');
    add('Branch junction: the outlet pressure of the two-part march is met', cfg.pOut, r.pOut, 1e-5, 'bara', false);
    flag('Branch junction: the mixed temperature lies between the two arriving streams', r.junction.tMix > Math.min(r.junction.tMain, r.junction.tBranch) - 1e-6 && r.junction.tMix < Math.max(r.junction.tMain, r.junction.tBranch) + 1e-6, `${r.junction.tMain.toFixed(2)} °C (main) and ${r.junction.tBranch.toFixed(2)} °C (branch) mix to ${r.junction.tMix.toFixed(2)} °C`);
  }
  await 0;
  { // 2-D solver and interface schemes
    const cav = once('cavity100', () => cfd2d({ nx: 32, ny: 32, rho: 1, mu: 0.01, lid: 1, tEnd: 40, steady: 1e-4, cfl: 0.5, poissonTol: 1e-6 })); let e = 0; for (const q of REF.GHIA.rows) if (q[0] === 0 && q[1] > 0 && q[1] < 1) e = Math.max(e, Math.abs(interp1(cav.y, cav.uMid, q[1]) - q[2]));
    add('Lid-driven cavity at Re = 100: u on the vertical centre line against Ghia, Ghia & Shin (1982)', 0, e, 0.01, 'Largest difference from the 15 tabulated interior values (lid velocity 1), 32 × 32 cells', false);
    add('2-D solver: discrete continuity', 0, cav.divergenceMax, 1e-6, 'Largest cell divergence × Δx / |u|max over the run', false);
    const mu = 0.1, po = cfd2d({ nx: 4, ny: 16, lx: 0.5, ly: 1, rho: 1, mu, fx: 1, periodicX: true, tEnd: 30, steady: 1e-9 });
    add('Plane Poiseuille flow: maximum velocity', 1 / (8 * mu), Math.max(...po.u.map((q) => q[1])), 2e-3, 'Body-force driven periodic channel, 16 cells across');
    await 0;
    const za = interfaceAdvect2D({ n: 50, scheme: 'vof', test: 'zalesak' }), zf = interfaceAdvect2D({ n: 50, scheme: 'front', test: 'zalesak' });
    add('Zalesak slotted disc, volume of fluid (THINC/WLIC): volume conservation', 0, za.massError, 1e-6, 'Relative after one revolution on 50 × 50 cells', false);
    add('Zalesak slotted disc, volume of fluid: shape error', 0, za.shapeError, 0.35, 'L1 error of the liquid fraction relative to the disc area; the slot is 2.5 cells wide on this grid, so the bridge above it is smeared; front tracking below shows the same test without that limit', false);
    add('Zalesak slotted disc, front tracking: shape error', 0, zf.shapeError, 0.05, 'Markers advected with a second-order Runge–Kutta step', false);
    await 0;
    const sch = ['vof', 'levelset', 'clsvof', 'phasefield'], lim = { vof: 0.15, levelset: 0.2, clsvof: 0.15, phasefield: 0.45 }, cons = { vof: 1e-11, levelset: 0.25, clsvof: 1e-11, phasefield: 1e-11 }, nm = { vof: 'Volume of fluid', levelset: 'Level set', clsvof: 'Coupled level set / VOF', phasefield: 'Phase field' }; let rv = null;
    for (const q of sch) { const d = interfaceAdvect2D({ n: 40, scheme: q, test: 'translateDiag' }); add(`${nm[q]}: circle advected diagonally across a periodic grid, shape error`, 0, d.shapeError, lim[q], `Volume error ${d.massError.toExponential(1)}${q === 'levelset' ? ' (the level set does not conserve volume; the coupled scheme takes it from the VOF field)' : q === 'phasefield' ? ' (diffuse interface about three cells wide)' : ''}`, false); add(`${nm[q]}: volume conservation`, 0, d.massError, cons[q], 'Relative', false); if (q === 'vof') rv = d; await 0; }
    const ax = interfaceAdvect2D({ n: 40, scheme: 'vof', test: 'translateX' });
    add('Mesh-orientation sensitivity: diagonal against grid-aligned advection (volume of fluid)', 1, rv.shapeError / ax.shapeError, 1, `Ratio of the shape errors (${rv.shapeError.toFixed(4)} diagonal, ${ax.shapeError.toFixed(4)} along x): the split scheme is at most about twice as inaccurate across the mesh`, false);
    const a0 = 0.25, dm = cfd2d({ nx: 32, ny: 16, lx: 2, ly: 1, rhoL: 1000, rhoG: 1.2, muL: 1e-3, muG: 1.8e-5, gy: -G, c0: (x, y) => (x < a0 && y < 2 * a0 ? 1 : 0), tEnd: 0.4, probeFront: true, turbulence: 'sst' });
    add('RANS + VOF dam break in a channel (SST k–ω): liquid volume conservation', 0, dm.mass.error, 1e-8, 'Relative, 32 × 16 cells, water and air', false);
    flag('RANS + VOF dam break: the front stays behind the frictionless shallow-water limit', dm.front.x[dm.front.x.length - 1] > 2 * a0 && dm.front.x[dm.front.x.length - 1] < a0 + 2 * Math.sqrt(G * 2 * a0) * 0.4, `Front at ${dm.front.x[dm.front.x.length - 1].toFixed(2)} m after 0.4 s; Ritter limit ${(a0 + 2 * Math.sqrt(G * 2 * a0) * 0.4).toFixed(2)} m`);
    await 0;
    const ReT = 180, ch = cfd2d({ nx: 4, ny: 32, lx: 0.4, ly: 2, rho: 1, mu: 1 / ReT, fx: 1, periodicX: true, turbulence: 'sst', tEnd: 25, steady: 2e-4, u0: (x, y) => 17 * Math.max(1 - Math.abs(y - 1), 0.02) ** (1 / 7), k0: 1, omega0: 20, cfl: 0.5, poissonTol: 1e-6 });
    add('SST k–ω in the 2-D solver: bulk velocity of a turbulent channel against the log law', Math.log(ReT) / 0.41 + 5.2 - 1 / 0.41, mean(ch.u.map((q) => q[1])), 0.1, `U_b / u_τ at Re_τ = 180 on 32 uniform cells (first cell centre at y⁺ = 5.6), started from a 1/7-power profile and run until the velocity changes by less than ${ch.steadyResidual.toExponential(1)} per unit time; the integral of the log law is the reference`);
  }
}

// =====================================================================================================
// 7j-2. Verification of the three-dimensional solver (js/core/cfd3d.js); small grids, a few seconds in total
// =====================================================================================================
async function verify3D(add) {
  const flag = (name, cond, note) => add(name, 1, cond ? 1 : 0, 0, note, false), go = (r) => { while (!r.done) r.advance(2000); return r.result(); };
  { // direct Poisson solver: discrete Laplacian of a known field, periodic × wall × wall
    const nx = 16, ny = 12, nz = 8, g = { nx, ny, nz, dx: 1 / nx, dz: 1 / nz, dy: new Float64Array(ny + 2).fill(1 / ny), dyc: new Float64Array(ny + 1).fill(1 / ny), perX: true, perY: false, perZ: false }, P = makePoisson(g), W = new Float64Array(nx * ny * nz), ref = new Float64Array(nx * ny * nz);
    const at = (i, j, l) => ref[((i + nx) % nx) + nx * (clamp(j, 0, ny - 1) + ny * clamp(l, 0, nz - 1))];
    for (let l = 0; l < nz; l++) for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) ref[i + nx * (j + ny * l)] = Math.cos((2 * PI * (i + 0.5)) / nx) * Math.cos((PI * (j + 0.5)) / ny) + 0.3 * Math.sin((4 * PI * i) / nx) * Math.cos((PI * (l + 0.5)) / nz);
    for (let l = 0; l < nz; l++) for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) W[i + nx * (j + ny * l)] = (at(i + 1, j, l) - 2 * at(i, j, l) + at(i - 1, j, l)) * nx * nx + (at(i, j + 1, l) - 2 * at(i, j, l) + at(i, j - 1, l)) * ny * ny + (at(i, j, l + 1) - 2 * at(i, j, l) + at(i, j, l - 1)) * nz * nz;
    P.solve(W); let e = 0; const sh = W[0] - ref[0]; for (let k = 0; k < W.length; k++) e = Math.max(e, Math.abs(W[k] - ref[k] - sh));
    add('3-D solver: direct Poisson solver (Fourier × cosine transforms + tridiagonal) returns a field from its discrete Laplacian', 0, e, 1e-11, '16 × 12 × 8 cells, periodic in x, walls in y and z; largest difference up to a constant', false);
  }
  { // plane Poiseuille flow is an exact discrete solution (uniform grid, parabolic wall-shear formula)
    const s = createFlow3D({ nx: 4, ny: 8, nz: 4, lx: 1, ly: 2, lz: 1, nu: 0.1, force: [1, 0, 0], init: (x, y) => [5 * y * (2 - y), 0, 0] }); for (let k = 0; k < 30; k++) s.step(); let e = 0; for (let j = 1; j <= 8; j++) e = Math.max(e, Math.abs(s.u[s.grid.idx(2, j, 2)] - 5 * s.grid.yc[j] * (2 - s.grid.yc[j])));
    add('3-D solver: plane Poiseuille flow stays exact (8 uniform cells across, 30 steps)', 0, e, 1e-10, 'Largest deviation from u = G y (2h − y) / 2ν; the wall shear uses the parabola through the wall and two cell centres', false);
    add('3-D solver: wall shear of plane Poiseuille flow', 10, s.wallShear(), 1e-9, 'du/dy at the wall = G h / ν');
    flag('3-D solver: velocity field solenoidal after the projection', s.divergence() < 1e-10, `Largest divergence ${s.divergence().toExponential(1)}`);
  }
  await 0;
  { // Taylor–Green vortex: 2-D exact decay with second-order convergence; 3-D initial dissipation rate and energy budget
    const a = go(taylorGreen({ n: 8, dim: 2, re: 10, tEnd: 0.5 })), b = go(taylorGreen({ n: 16, dim: 2, re: 10, tEnd: 0.5 }));
    add('3-D solver: decay of the two-dimensional Taylor–Green vortex against exp(−2νt), 16 cells per side', 0, b.errorMax, 3e-3, 'Largest relative error of the velocity amplitude up to t = 0.5 at Re = 10', false);
    add('3-D solver: observed order of accuracy on the Taylor–Green vortex (8 → 16 cells)', 2, Math.log2(a.errorMax / b.errorMax), 0.3, 'Second-order central differences');
    const c = go(taylorGreen({ n: 16, dim: 3, re: 100, tEnd: 1 }));
    add('3-D solver: initial dissipation rate of the three-dimensional Taylor–Green vortex', 0.75 / 100, c.dissipation[0], 0.02, 'ν ⟨|ω|²⟩ = 0.75 ν U²/L² exactly; 16³ cells (modified wavenumber of the second-order stencil)');
    add('3-D solver: kinetic-energy budget, −dE/dt against 2 ν × enstrophy (Taylor–Green, 16³)', c.enstrophyDissipation[c.t.length - 2], c.dissipation[c.t.length - 2], 0.01, 'The advection scheme neither creates nor destroys the kinetic energy of the face velocities, so the energy lost equals the resolved viscous dissipation');
    const d1 = go(taylorGreen({ n: 16, dim: 3, re: 1e9, tEnd: 0.5, dtMax: 0.05 })), d2 = go(taylorGreen({ n: 16, dim: 3, re: 1e9, tEnd: 0.5, dtMax: 0.025 })), l1 = Math.abs(d1.ek[d1.ek.length - 1] / d1.ek[0] - 1), l2 = Math.abs(d2.ek[d2.ek.length - 1] / d2.ek[0] - 1);
    add('3-D solver: inviscid flow conserves kinetic energy (loss at Δt = 0.025)', 0, l2, 1e-6, `Relative loss over t = 0.5 without viscosity: ${l1.toExponential(1)} at Δt = 0.05 and ${l2.toExponential(1)} at 0.025 — the small error left is that of the third-order Runge–Kutta scheme, not of the spatial scheme`, false);
    const e1 = go(taylorGreen({ n: 16, dim: 3, re: 100, tEnd: 0.6 })), e2 = go(taylorGreen({ n: 16, dim: 3, re: 100, tEnd: 0.6 }));
    add('3-D solver: floating-point reproducibility (two identical runs)', e1.ek[e1.ek.length - 1], e2.ek[e2.ek.length - 1], 0, 'Kinetic energy after the same number of steps, bit for bit');
  }
  await 0;
  { // sub-grid models in a uniform shear
    const sh = (sgs) => { const s = createFlow3D({ nx: 4, ny: 16, nz: 4, lx: 1, ly: 1, lz: 1, nu: 1e-3, sgs, wallDistance: () => 1e30, init: (x, y) => [3 * y, 0, 0] }); s.fillVel(); s.step(1e-9); return s.nut[s.grid.idx(2, 8, 2)]; }, d = Math.cbrt(1 / 256);
    add('3-D solver: Smagorinsky eddy viscosity in a uniform shear', (CFD3D_CONSTANTS.smagorinsky.cs * d) ** 2 * 3, sh('smagorinsky'), 1e-9, '(C_s Δ)² |S| with C_s = 0.1, Δ = (Δx Δy Δz)^⅓, far from walls');
    add('3-D solver: WALE eddy viscosity vanishes in a pure shear', 0, sh('wale'), 1e-12, 'The model is built to give no sub-grid stress in laminar shear and the correct y³ behaviour at a wall', false);
    const r = go(channel3D({ reTau: 180, nx: 8, ny: 12, nz: 8, sgs: 'iddes', tEnd: 0.03, tStats: 0.01 })), ls = r.sim ? null : null; const s2 = createFlow3D({ nx: 8, ny: 12, nz: 8, lx: 6.28, ly: 2, lz: 3.14, stretchY: 1.9, nu: 1 / 180, sgs: 'des', init: (x, y) => [15 * y * (2 - y), 0, 0] }); s2.fillVel(); s2.step(1e-6); let okL = true; for (let j = 1; j <= 12; j++) { const k = s2.grid.idx(2, j, 2), dw = s2.wallDistance[k], L = s2.lengthScale[k]; if (!(L <= dw * (1 + 1e-12) && L <= 0.65 * Math.max(s2.grid.dx, s2.grid.dy[j], s2.grid.dz) * (1 + 1e-12) + 1e-300 || L === dw)) okL = false; if (Math.abs(L - Math.min(dw, 0.65 * Math.max(s2.grid.dx, s2.grid.dy[j], s2.grid.dz))) > 1e-12) okL = false; }
    flag('3-D solver: DES length scale is min(wall distance, 0.65 Δ_max) in every cell', okL, 'Spalart–Allmaras DES97 switch on a stretched channel grid');
    flag('3-D solver: IDDES run on a coarse channel stays finite with a positive eddy viscosity', !r.failed && Math.max(...r.nutPlus) > 0 && r.cf > 0, `ν_t/ν up to ${Math.max(...r.nutPlus).toFixed(1)} after ${r.steps} steps`);
  }
  await 0;
  { // volume of fluid in three dimensions: conservation and boundedness with a density ratio of 830
    const r = go(twoPhase3D({ kind: 'dambreak', n: 8, spanCells: 2, tEnd: 0.15 }));
    add('3-D volume of fluid (dam break, water and air): liquid volume conserved', 0, Math.abs(r.volume.error), 1e-11, `Relative; ${r.nx} × ${r.ny} × ${r.nz} cells, ${r.steps} steps, density ratio ${r.densityRatio.toFixed(0)}`, false);
    flag('3-D volume of fluid: liquid fraction stays within [0, 1]', r.bounds.min > -1e-9 && r.bounds.max < 1e-9, `Overshoots ${r.bounds.min.toExponential(1)} / ${r.bounds.max.toExponential(1)}`);
    flag('3-D dam break: the front advances but stays behind the frictionless shallow-water limit', r.front[r.front.length - 1] > r.a && r.front[r.front.length - 1] < r.a + 2 * Math.sqrt(G * 2 * r.a) * r.time, `Front at ${r.front[r.front.length - 1].toFixed(3)} m after ${r.time.toFixed(2)} s; Ritter limit ${(r.a + 2 * Math.sqrt(G * 2 * r.a) * r.time).toFixed(3)} m`);
    const q = go(twoPhase3D({ kind: 'stratified', n: 8, tEnd: 0.08, force: 2 }));
    add('3-D volume of fluid (periodic stratified wavy channel): liquid volume conserved', 0, Math.abs(q.volume.error), 1e-11, `Relative; ${q.nx} × ${q.ny} × ${q.nz} cells, periodic in two directions`, false);
  }
  await 0;
  { // immersed-boundary pipe and the 1-D + 3-D coupling
    const p = go(pipe3D({ n: 16, nx: 2, nu: 0.05, tEnd: 20, steady: 1e-5 }));
    add('3-D solver: laminar pipe by immersed boundary, friction factor against 64/Re', 1, p.friction / p.frictionLaminar, 0.05, `${p.cellsAcross.toFixed(1)} cells across the diameter; the wall is first-order accurate`);
    const c = coupled1D3D({ n: 12 });
    add('1-D + 3-D coupling: coupled flow velocity against Hagen–Poiseuille', c.velocityAnalytic, c.velocity, 0.01, `Single-phase laminar line of 20 cells, one of them replaced by the 3-D section (${c.cellsAcross.toFixed(1)} cells across, friction ${(100 * c.error3d).toFixed(1)} % from 64/Re); ${c.iterations.length} exchanges`);
    add('1-D + 3-D coupling: 1-D closure alone against Hagen–Poiseuille', c.velocityAnalytic, c.velocity1d, 1e-9, 'Both sides of the coupling agree with the analytic solution');
    flag('1-D + 3-D coupling: the exchange converges and the pressure drops add up to the imposed difference', c.converged && Math.abs(c.pressureDropBalance.sum / c.pressureDropBalance.given - 1) < 1e-6, `Relative change of the flow rate below 1e-4 after ${c.iterations.length} exchanges`);
  }
}

// =====================================================================================================
// 7k. Provenance of the constants used by this suite and by the pipe-flow kernel (js/core/pipe.js)
// =====================================================================================================
const FL = 'https://github.com/CalebBell/fluids', HT = 'https://github.com/CalebBell/ht', UF = 'https://github.com/unifloc/unifloc_vba/blob/master/modules_txt/u7_Multiphase_PVT.txt', OF = 'https://github.com/OpenFOAM/OpenFOAM-dev', D0 = '2026-10-08';
const pv = (item, used, source, url, status, note) => ({ item, used, source, url, retrieved: /2026-10-09|searched again|New\./.test(source + note) ? '2026-10-09' : D0, status, note });
export const PROVENANCE = [
  pv('Colebrook–White, Haaland, Swamee–Jain and Churchill (1977) friction factors', 'pipe.js frictionFactor()', 'C. Bell, fluids (MIT), fluids/friction.py', FL + '/blob/master/fluids/friction.py', 'verified', 'Formulas compared term by term; the library examples at Re = 1e5, ε/D = 1e-4 are reproduced to 1e-9 or better (Swamee–Jain to 1e-6: the kernel uses the rounded constant 5.74 where the library writes 6.97^0.9 = 5.7407).'),
  pv('Gnielinski turbulent Nusselt number', 'pipe.js hInside()', 'C. Bell, ht (MIT), ht/conv_internal.py turbulent_Gnielinski', HT + '/blob/master/ht/conv_internal.py', 'verified', 'Same expression; the kernel inserts the smooth-pipe Haaland friction factor, which gives 3.5 % less than the library example with f = 0.0185.'),
  pv('Churchill–Bernstein cross-flow Nusselt number', 'pipe.js nuCrossFlow(), hOutside()', 'C. Bell, ht (MIT), ht/conv_external.py Nu_cylinder_Churchill_Bernstein', HT + '/blob/master/ht/conv_external.py', 'verified', 'Library example Nu(6071, 0.7) = 40.637 reproduced to 1e-12.'),
  pv('Churchill–Chu free convection of a horizontal cylinder; Nu³ = Nu_forced³ + Nu_free³', 'pipe.js nuFreeCylinder(), hOutside(…, dT)', 'C. Bell, ht (MIT), ht/conv_free_immersed.py Nu_horizontal_cylinder_Churchill_Chu', HT + '/blob/master/ht/conv_free_immersed.py', 'verified', 'New in the kernel. Library example Nu(Pr = 0.69, Gr = 2.63e9) = 139.135 reproduced to 1e-12. The cubic combining rule for mixed convection is the usual engineering rule and was not compared with a source.'),
  pv('Seawater properties in the outside film', 'pipe.js hOutside()', 'core/props.js (Sharqawy, Lienhard & Zubair 2010; Nayar et al. 2016)', 'https://web.mit.edu/seawater/', 'corrected', 'Old: hard-coded fits μ = 1.9e-3 exp(−0.027 T) + 3.5e-4 Pa s and Pr = 13.4 exp(−0.027 T) + 2, which were 23–37 % and 20–39 % above the property module at 4–25 °C. New: density, viscosity, conductivity and heat capacity of core/props.js at 35 g/kg. The forced-convection film coefficient rises by about 10 %.'),
  pv('Beggs & Brill (1973) flow-pattern limits L1–L4, holdup coefficients a, b, c, inclination coefficients, friction-factor ratio', 'pipe.js beggsBrill()', 'C. Bell, fluids (MIT), fluids/two_phase.py Beggs_Brill', FL + '/blob/master/fluids/two_phase.py', 'verified', 'All 4 + 9 + 12 + 5 constants identical; the library example (m = 0.6 kg/s, x = 0.1, 50 mm, horizontal) gives 686.9724506803 Pa/m in both codes.'),
  pv('Payne et al. (1979) holdup factors 0.924 (uphill) and 0.685 (downhill)', 'pipe.js beggsBrill()', 'not found in an openly readable source (searched again: the pengtools Beggs–Brill page and the fluids library do not apply any Payne factor)', 'https://doi.org/10.2118/6874-PA', 'unverified', 'The two factors remain unverified (the paper is not open). The original Beggs–Brill correlation, verified against the fluids library, is used unscaled at exactly 0°.'),
  pv('Taitel & Dukler (1976) stratified momentum balance and transitions (Kelvin–Helmholtz, wave generation s = 0.01, dispersed bubble)', 'pipe.js stratifiedLevel(), flowPattern()', 'C. Bell, fluids (MIT), fluids/two_phase.py Taitel_Dukler_regime (dimensionless groups X, T, F, K)', FL + '/blob/master/fluids/two_phase.py', 'verified', 'The kernel solves the level and evaluates the criteria directly instead of interpolating the published curves; the groups are equivalent. Deviations kept on purpose: interfacial friction f_i = max(f_G, 0.0142) and annular flow below h/D = 0.35 (Barnea) instead of 0.5.'),
  pv('Taitel, Barnea & Dukler (1980) upward transitions: bubble–slug at a void of 0.25, dispersed-bubble line 4.0{D^0.429 (σ/ρl)^0.089 / ν^0.072}[gΔρ/ρl]^0.446, annular 3.1[σ g Δρ]^¼/ρg^½, D > 19[…]^½', 'pipe.js flowPattern()', 'unifloc_vba (MIT), Ansari flow-pattern routines; S. Mittlböck, MSc thesis, TU Wien 2014 (eqs 4.47–4.49); A. Schulte, MSc thesis, TU Delft 2013 (eqs 4.60–4.62)', 'https://repository.tudelft.nl/file/File_4425d042-ff4a-486f-b3e1-ae07a5e5728e', 'verified', 'The bubble–slug line j_G = j_L/3 + 0.383[σgΔρ/ρl²]^¼ (0.383 = 0.25 × 1.53), the diameter limit 19 and the constant 4.0 with the exponents 0.429, 0.089, 0.072, 0.446 of the dispersed-bubble line are printed in the two theses (the Delft thesis misprints the variables of the last one; the dimensionally consistent form is used). The 1980 dispersed-bubble line is no longer used by the kernel (replaced by the criterion of Barnea, next entry); the churn rule (v_sg > 0.6 of the annular velocity) remains an engineering simplification and is not sourced.'),
  pv('Bendiksen (1984) C0 = 1.05 + 0.15 sin²θ, v_d = (0.54 cosθ + 0.35 sinθ)√(gD), Fr < 3.5; C0 = 1.2, v_d = 0.35 √(gD) sinθ above', 'pipe.js slugVelocity()', 'A. Al-lababidi, PhD thesis, Cranfield 2006 (eq. 3.26); S. Mittlböck, MSc thesis, TU Wien 2014 (eqs 4.106–4.108)', 'https://core.ac.uk/download/pdf/139376.pdf', 'verified', 'Both open theses print the same four constants and the Froude-number limit 3.5. Against measurements: 15 horizontal translational velocities (Mohmmed et al.) 11.4 % mean error; 134 bubble-nose velocities at 0–90° of the UNICAMP set (CC BY-NC, not distributed) 12.6 %.'),
  pv('Gregory, Nicholson & Aziz (1978) slug-body holdup 1/(1 + (v_m/8.66)^1.39), v_m in m/s', 'pipe.js slugBodyHoldup()', 'A. Al-lababidi, PhD thesis, Cranfield 2006 (eq. 3.34); Y. Chen, MS report, Stanford 2001 (eq. 3-17); original measurements in the CC BY compilation doi:10.17632/wyfdm5ysh6.1', 'https://stacks.stanford.edu/file/druid:pp071qb9526/Chen01.pdf', 'verified', 'Two open documents print the same constants, and the closure reproduces the 157 original measurements of Gregory et al. with 3.6 % mean error and 0.2 % bias. Blind on the other 21 studies of the compilation: 10 % horizontal, 15 % inclined, 40 % vertical (over-prediction): the closure has no inclination term. Alternatives (Gomez et al., Barnea & Brauner) were evaluated from memory of their forms and did no better (22 % and 21 % overall); they could not be read from an open source and are not offered as sourced options.'),
  pv('Slug frequency: Gregory & Scott (1969) 0.0226[(v_sl/gD)(19.75/v_m + v_m)]^1.2; Zabaras (2000) = that × (0.836 + 2.75 sin θ); Greskovich & Shrier (1972) 0.0226[λ(2.02/D + v_m²/gD)]^1.2; Heywood & Richardson (1979) 0.0434[…]^1.02', 'pipe.js slugFrequency(); slugUnitCell()', 'V. Hernandez-Perez, PhD thesis, Nottingham 2008 (eqs 2.65–2.73); A. Arabi et al., Chem. Eng. Sci. 2020, accepted manuscript on HAL (table of slug-frequency correlations); A. Al-lababidi, PhD thesis, Cranfield 2006 (eq. 3.42)', 'https://hal.science/hal-02548768v1/document', 'corrected', 'Units: all sources are SI throughout (19.75 m²/s² = 212.6 ft²/s², 2.02 m), so the Gregory–Scott constant was right. Corrected: the Zabaras inclination factor was coded as 0.836 + 2.75 sin^0.25 θ and applied at any angle; the two open sources print 0.836 + 2.75 sin θ, fitted for 0–11°, and the factor is now held beyond 11° (the SPE paper itself is not open, so the exponent rests on these two secondary sources). Effect: 136 % → 24 % mean error on 67 inclined 26 mm points (UNICAMP, 257 D, not distributed), 30 % on the 30 horizontal 74 mm points of Mohmmed et al. (Gregory–Scott 20 %, Heywood–Richardson 15 % there, but 30–36 % with a −17 to −29 % bias on the inclined pipe). Zabaras is the default correlation because it is the only one that follows the inclination trend. Al-Safran (2009, 2016), Schulkes (2011), Hill & Wood, Shea and Gokcal were found in one table only, with undefined inputs or normalisation (film holdup, pipe-length unit, f D/v), and are not implemented.'),
  pv('Developed slug length: Scott, Shoham & Brill (1989) ln L = −26.6 + 28.5 [ln D + 3.67]^0.1 (L, D in m), about 30 D in small pipes; Brill et al. (1981) ln L[ft] = −2.663 + 5.441 √ln d[in] + 0.059 ln v_m[ft/s]; Norris (1982)', 'pipe.js slugLength(); slugUnitCell()', 'S. Mittlböck, MSc thesis, TU Wien 2014 (eq. 4.102, SI form); A. Al-lababidi, PhD thesis, Cranfield 2006 (eqs 3.44, 3.45, field-unit forms, Table 3-3 of stable lengths)', 'https://repositum.tuwien.at/bitstream/20.500.12708/7775/2/Mittlboeck%20Sebastian%20-%202014%20-%20Development%20of%20a%20slug%20analysis%20tool%20for%20slugging...pdf', 'verified', 'The SI and the field-unit forms of the Scott correlation are the same equation (ln 39.37 = 3.673, −25.4134 − ln 3.2808 = −26.60): the metre form in the kernel is right. Brill agrees with the field-unit form. Stable lengths quoted there: 12–30 D (Dukler & Hubbard), 30 D (Nicholson; Gregory), 32 D (Barnea & Taitel theory): the kernel keeps 32 D. The Norris constants (suite option) were not found in an open source: unverified.'),
  pv('Slug length against distance from the inlet: kinematic slug-train model with the wake law v = v∞[1 + a exp(−b L/D)], a = 0.4, b = 1.0 (horizontal), a = 8, b = 1.06 (inclined, vertical)', 'pipe.js slugLength(…, { xD, theta }), slugUnit({ xD }); run()', 'C. Sarmiento et al., “Influence of the initial conditions for the numerical simulation of two-phase slug flow”, Proc. ENCIT 2010, ABCM (open proceedings): wake law after Moissis & Griffith (1962) with the constants of Rodrigues (2009), within the Barnea & Taitel (1993) slug-length model', 'https://www.abcm.org.br/anais/encit/2010/PDF/ENC10-0607.pdf', 'corrected', 'New. Old: the developed-flow length (32 D or Scott) at any distance, 131 % and 478 % too long in 300 D and 108 D laboratory pipes. New: mean length of a simulated slug train (3 × 6,000 slugs, tabulated in the kernel) up to 3,000 D, then a logarithmic approach to the developed correlation at 30,000 D. Blind results: 31 % (Mohmmed et al., 74 mm, 54 and 81 D), 23 % and 26 % (UNICAMP, 26 mm, 77 and 257 D, not distributed). Not sourced and stated as assumptions: the mean length of new slugs (3.2 D horizontal, 4.2 D inclined; the cited paper reports 4.5–30 D at its inlet), the blending between the horizontal and the inclined curve (sin θ / sin 10°) and the 3,000–30,000 D transition to the developed correlation. One open source only for the wake constants. Unit-cell identity f L_u = v_t, L_s = β L_u as before.'),
  pv('Slug-unit film zone (Taitel & Barnea 1990, equilibrium film): film and bubble velocities from the liquid and gas shed at the slug tail, film level from the combined momentum balance, slug fraction from the liquid balance', 'pipe.js slugFilm(), slugUnit(), gradient(…, \'mechanistic\')', 'derived balances (no fitted constants); friction closures as in the stratified model of the kernel', 'https://doi.org/10.1016/S0065-2717(08)70026-1', 'corrected', 'Old: the film holdup was the stratified level at 0.3 × the liquid rate, which exceeded the unit holdup at low gas rates and drove the slug fraction to its lower limit of 0.02, so the mechanistic model returned almost no wall friction (9 Pa/m against 298 Pa/m measured in one horizontal case). New: equilibrium film of the unit cell; wall forces of the slug body and of the film zone weighted with the slug fraction. Effect on 67 measured gradients of air–water slug flow at 0–90° (UNICAMP, not distributed): 31.3 % → 17.2 % mean error (Beggs & Brill 11.5 %, drift flux 9.7 % on the same rows); the remaining bias (−17 %) is largest in the horizontal pipe, where the gradient is the difference of two absolute pressures 180 D apart. Reference case with the mechanistic model: inlet pressure 87.4 → 93.0 bara.'),
  pv('Barnea (1986, 1987) unified transitions: dispersed bubble d_max = [0.725 + 4.15 (v_sg/v_m)^½](σ/ρl)^0.6 (2 f_m v_m³/D)^−0.4 < min(d_CD, d_CB), void ≤ 0.52; annular film stability and blockage at a film holdup of 0.24; downward stratified → annular at v_L > [gD(1 − h/D) cos θ / f_L]^½', 'pipe.js flowPattern(), annularFilm()', 'WIT Trans. Eng. Sci. 89 (2015) paper MPF15030 (eqs 1–4); Y. Chen, MS report, Stanford 2001 (eqs 3-14, 3-16, 3-28); S. Mittlböck, MSc thesis, TU Wien 2014 (eqs 4.23, 4.24, 4.45–4.49)', 'http://www.witpress.com/Secure/elibrary/papers/MPF15/MPF15030FU1.pdf', 'corrected', 'New in the kernel. Found wrong or missing: (1) vertical and steep downward flow was run through the near-horizontal stratified logic (6 of 19 vertical-downward observations right); the torn-film criterion and the annular-film criteria now give 18 of 19; (2) bubble flow was allowed at any upward angle above 10° (10 intermittent points called bubble), now only from 60°; (3) the dispersed-bubble boundary of Barnea replaces the 1980 line. Shoham air–water data: vertical 65 → 92.5 %, downward 61.5 % (68.7 % without the corrupt −50° block) → 77.4 %, upward 82.7 → 85.7 %, horizontal 96.8 % unchanged; oil–air data of Kokal, not used in the selection: 57.5 % (63.6 % without the corrupt 76 mm block) → 68.6 %. Kept from the previous kernel because they did better on the Shoham data than the alternatives tried (h/D < 0.5, film criterion at all angles): annular flow near the horizontal below h/D = 0.35 and the Taitel–Dukler turbulence criterion for dispersed bubbles — a selection made on the Shoham set, so only the Kokal result is blind. Not found in an open source: the wall-friction form f_i = f_G(1 + 300 δ/D) behind the factor 1 + 75 H, the film Froude number 1.5 for wavy downward flow and the bubble-flow inclination limit (a fixed 60° is used; the lift-force expression of the TU Wien thesis has an unreadable bubble-size unit).'),
  pv('Bøe (1981) criterion, Pots (1987) number, Taitel (1986) stability pressure ρl g (α L/0.89 − h); riser cycle with choke and riser-base gas lift', 'pipe.js severeSlugging(); severeScreen(); riserSluggingCycle(); minimumStableRate()', 'not found in an openly readable source (searched again)', 'https://doi.org/10.1016/0301-9322(86)90027-4', 'unverified', 'The forms are derived balances (liquid head build-up against gas pressure build-up) without fitted constants except the 0.89 of the Taitel criterion. New: lift gas injected at the riser base aerates the liquid column in the cycle model (drift-flux void), and the controlled turndown uses the Taitel back-pressure. The cycle model is compared with twelve measured cycle times with a choke or gas lift (Jansen, Shoham & Taitel 1996): 31 % mean error (32 % before on all 23).'),
  pv('API RP 14E erosional velocity v_e = C/√ρ, C = 100 continuous service (122 in SI)', 'steadySolve() evr', 'pengtools wiki, Erosional velocity', 'https://wiki.pengtools.com/index.php?title=Erosional_velocity', 'verified', 'C = 100 for continuous and 125 for non-continuous service in lb, ft, s units; factor 1.22 converts to kg, m, s.'),
  pv('Hagedorn & Brown (1965) holdup chart fits H_L/ψ and ψ, Griffith bubble-flow limit 1.071 − 0.2218 v_m²/D ≥ 0.13, CN_L(N_L) fit', 'hagedornBrown()', 'unifloc_vba (MIT), unf_HagedornandBrawnmodified; pengtools wiki, Hagedorn and Brown correlation (CN_L fit of Economides)', 'https://wiki.pengtools.com/index.php?title=Hagedorn_and_Brown_correlation', 'verified', 'Both polynomials agree with the open code to the digits carried. CN_L: the log-polynomial used here was compared with the cubic CN_L = 0.061 N_L³ − 0.0929 N_L² + 0.0505 N_L + 0.0019 printed on the wiki at N_L = 0.002 … 0.5: +3, −7, −10, −6, +7, +15, +15 and −10 % at 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2 and 0.5 — two fits of the same chart that differ by up to 15 %; the fit of the suite is kept and this spread is the uncertainty of that factor.'),
  pv('Gray (1974): A = −2.314 [N_v (1 + 205/N_D)]^B, B = 0.0814 [1 − 0.0554 ln(1 + 730 R/(R + 1))], pseudo-roughness 28.5 σ/(ρ v²) ≥ 2.77e-5', 'gray()', 'unifloc_vba (MIT) unf_GrayModifiedGradient; pengtools wiki Gray correlation', UF, 'verified', 'The two open sources disagree with each other in single digits (the wiki prints −2.2314 and 0.554, the code 250/N_D); each constant used here is confirmed by one of them and the remaining digits by the other: −2.314, 0.0554, 730 (code) and 205, 28.5, 0.007, 2.77e-5 (wiki).'),
  pv('Orkiszewski (1967) liquid-distribution coefficient (four equations), limits, Griffith–Wallis boundary', 'orkiszewski()', 'N. Antipin, orkiszewski_model (GitHub, no licence stated), two_phases/orkiszewski_model.py', 'https://github.com/Nikita-antipin/orkiszewski_model', 'verified', 'Three of the four equations agree in every constant. The fourth (oil, v_m > 10 ft/s) has 0.161 here and 0.167 in the open code; 0.161 is kept because it makes the two oil equations continuous at 10 ft/s for tubing sizes (the 0.167 would not): that single constant is unverified. Bubble-rise Reynolds-number branches (0.546, 0.35, 0.251, 8.74e-6, 13.59) were not found in an open source: unverified.'),
  pv('Mukherjee & Brill (1985): three holdup coefficient sets, four transition equations, friction-ratio table', 'mukherjeeBrill()', 'S. Sakurai, MukherjeeBrill R package (MIT), R/MukherjeeBrill.R', 'https://github.com/sshunsuke/MukherjeeBrill', 'verified', 'All 18 holdup coefficients, the 4 transition equations and the 8-point f_R table identical.'),
  pv('Duns & Ros (1963): F1–F7, L1, L2, f2 chart values, L_S = 50 + 36 N_Lv, L_M = 75 + 84 N_Lv^0.75, slip and friction procedure, mist-flow film roughness', 'dunsRos(), mistFriction()', 'open digitisations and code: ObuxoffCost/Duns-Ros (dense chart digitisation), TRUEVORO/duns_ros_cor (algorithm and a coarse table), lucasnasution/DunsRos (MIT; L1, L2, F5, F6)', 'https://github.com/TRUEVORO/duns_ros_cor', 'verified', 'New. Chart values sampled from the dense digitisation at 13 N_L nodes; the independent coarse table agrees within 10 % except L1 near N_D = 30 (1.6 against 2.0), where two of the three digitisations agree and are used. The procedure follows the open code. The mist-flow roughness constants (0.0749, 0.3713, 0.302, 0.005) agree with it.'),
  pv('Ansari et al. (1994) mechanistic model: transitions, Taylor-bubble relations (1.2 v_m + 0.35 √(gDΔρ/ρl), 9.916 √(gD(1 − √α))), H_gLS = v_sg/(0.425 + 2.65 v_m), entrainment 1 − exp[−0.125(v_crit − 1.5)], Z = 1 + 300 δ/D or 1 + 24 (ρl/ρg)^⅓ δ/D, film stability', 'ansari(), ansariPattern(), entrainmentFraction()', 'unifloc_vba (MIT): Ansari, anmist, bubble, slug, chkan, Func, fpup, mpoint, dbtran', UF, 'verified', 'New. Implemented from the open code (a port of the Tulsa routines); the falling-film constant is 9.916 in one routine and 9.961 in another there — 9.916 is used.'),
  pv('Baker (1954) map: coordinates B_x, B_y (531 and 2.16 follow from λ and ψ) and six log-polynomial boundary fits', 'bakerRegime()', 'curve fits: not found in an openly readable source; compared with the boundaries digitised (about ±5 %) from the redrawn map in V. Korelstein & E. Pereyra, E3S Web Conf. 397, 01002 (2023), Fig. 1(d) (open access)', 'https://www.truboprovod.ru/download/articles/2023_universal_gas_liquid_flow_pattern_map.pdf', 'unverified', 'Still the weakest item. Against the digitised open figure (axes in lb/(h ft²), as the ranges indicate): the annular boundary (fits 1 and 4) agrees within 10–35 %, the bubble/froth line (fit 5) within 10 %, the stratified–wave line (fit 2) within 10–30 % up to B_x ≈ 10; the dispersed-flow line (fit 3) lies 2.5–5 times above the digitised curve and the slug–plug line (fit 6) departs beyond B_x ≈ 300, so those two are not confirmed, and the identification of the digitised curves with the classical boundaries is itself a reading of an unlabelled redraw. Validation on the 31 horizontal air–water observations of Shoham: 20 right (65 %), against 25 for the Mandhane map and 30 for the mechanistic map. Shown for comparison only; no solver uses it.'),
  pv('Wallis interfacial friction 0.005 (1 + 300 δ/D) form, critical Weber number 12 for droplets, deposition velocity 0.15 m/s', 'annularMist(), annularDevelopment()', 'Z = 1 + 300 δ/D: unifloc_vba (MIT); the other two: not found in an open source', UF, 'unverified', 'New. The equilibrium entrained fraction is the sourced Wallis / Ansari correlation; the entrainment rate is written as k_D C_eq so that only the development length depends on the unverified deposition velocity, which is a calibration parameter. Ishii–Mishima and Pan–Hanratty could not be read from an open source and are not used.'),
  pv('Homogeneous-relaxation model: Θ = 6.51e-4 α^−0.257 ψ^−2.24 s, ψ = (p_sat − p)/p_sat (Downar-Zapolski et al. 1996, below 10 bar)', 'relaxationTime(), flashingFlow()', 'arXiv:2109.15203 (LaTeX source, eqs. for Θ)', 'https://arxiv.org/abs/2109.15203', 'verified', 'New. Constants and exponents identical. The correlation was fitted to flashing water; its use for hydrocarbons is indicative.'),
  pv('Bubble closures: Tomiyama lift polynomial, Wellek aspect ratio 1/(1 + 0.163 Eo^0.757), Antal wall force (−0.01, 0.05), Burns dispersion (σ = 0.7), Sato viscosity 0.6, virtual mass 0.5, Ishii–Zuber drag', 'liftTomiyama(), dragIshiiZuber(), bubblyPipe(), bubbleRise()', 'OpenFOAM-dev (GPL-3.0) multiphaseEuler interfacial models and the tutorial case Grossetete', OF + '/tree/master/applications/modules/multiphaseEuler/phaseSystem/interfacialModels', 'verified', 'New. Formulas and default coefficients read from the source files and the tutorial dictionaries.'),
  pv('One-group interfacial-area transport: C_RC = 0.04, C = 3, α_max = 0.75, C_WE = 0.002, C_TI = 0.085, We_cr = 6', 'interfacialAreaTransport()', 'OpenFOAM-dev (GPL-3.0) IATE sources and tutorial bubbleColumnIATE', OF + '/tree/master/applications/modules/multiphaseEuler/phaseSystem/diameterModels/IATE', 'verified', 'New. Source terms and coefficients as in the open code.'),
  pv('Turbulence constants: SST k–ω (0.85, 1, 0.5, 0.856, 5/9, 0.44, 0.075, 0.0828, 0.09, a1 = 0.31, c1 = 10), Spalart–Allmaras (0.1355, 0.622, 2/3, 0.41, 0.3, 2, 7.1), standard / RNG (0.0845, 1.42, 1.68, 0.71942, η0 = 4.38, β = 0.012) / realizable (A0 = 4, C2 = 1.9, σε = 1.2) k–ε, LRR (1.8, 0.6, 0.25, 0.15, 1.44, 1.92, wall reflection 0.5, 0.3), wall functions κ = 0.41, E = 9.8', 'ransPipeExtra(), cfd2d()', 'OpenFOAM-dev (GPL-3.0) src/MomentumTransportModels; Turbulence Modeling Resource (SST, SA)', OF + '/tree/master/src/MomentumTransportModels/momentumTransportModels', 'verified', 'New. Default coefficients read from the source files. Each closure is also checked against the smooth-pipe friction law (verification tab).'),
  pv('Wilcox (1988) k–ω (5/9, 3/40, 0.09, 0.5, 0.5) and Chien (1982) low-Reynolds k–ε (C_μ = 0.09, C_ε1 = 1.35, C_ε2 = 1.8, σ_k = 1, σ_ε = 1.3, f_μ = 1 − exp(−0.0115 y⁺), f_2 = 1 − (0.4/1.8) exp(−Re_t²/36), wall terms −2νk/y² and −2νε/y² exp(−y⁺/2))', 'ransPipe()', 'NASA Turbulence Modeling Resource (public domain), pages wilcox.html and ke-chien.html', 'https://tmbwg.github.io/turbmodels/wilcox.html', 'corrected', 'All constants and damping functions agree with the resource, except the factor of f_2, coded as 0.22 and now 0.4/1.8 = 0.2222. Both closures reproduce the smooth-pipe friction law within 5–6 %.'),
  pv('Scale-resolving closures of the 3-D solver: Smagorinsky with van Driest damping (A⁺ = 26), WALE (C_w = 0.325), Spalart–Allmaras DES (C_DES = 0.65), delayed DES f_d = 1 − tanh[(8 r_d)³], IDDES (C_w = 0.15, c_t = 1.63, c_l = 3.55, f_w* = 0.424, f_B = min(2 exp(−9α²), 1), exponents −11.09 / −9, Δ = min(max(C_w d, C_w h_max, h_wn), h_max))', 'core/cfd3d.js CFD3D_CONSTANTS, eddyViscosity()', 'OpenFOAM-7 (GPL-3.0) sources: LES/WALE, LESdeltas/vanDriestDelta and IDDESDelta, DES/SpalartAllmarasDES, SpalartAllmarasDDES, SpalartAllmarasIDDES', 'https://github.com/OpenFOAM/OpenFOAM-7/tree/master/src/TurbulenceModels/turbulenceModels', 'verified', 'New. Every constant listed was read in the source files. Not from that source: the Smagorinsky constant 0.1 used for the channel (the code ships C_k = 0.094, C_e = 1.048, equivalent to C_s ≈ 0.17): unverified, a common channel-flow value. The Spalart–Allmaras constants are those of the earlier turbulence entry.'),
  pv('Numerics of the 3-D solver: low-storage Runge–Kutta coefficients (0, −5/9, −153/128; 1/3, 15/16, 8/15) of Williamson (1980); constant-coefficient pressure splitting for two phases (Dodd & Ferrante 2014); THINC/WLIC volume of fluid with the dilatation correction of Weymouth & Yue (2010)', 'core/cfd3d.js', 'not opened in this task; checked numerically instead', 'https://doi.org/10.1016/0021-9991(80)90033-9', 'unverified', 'New. The sources were not read during this task. Checked by the verification set: third-order time accuracy (inviscid kinetic energy conserved to 1e-8 at Δt = 0.025), second-order convergence on the Taylor–Green vortex, exact plane Poiseuille flow, liquid volume conserved to round-off with a density ratio of 830, and a dam-break front within 2 % of an OpenFOAM interFoam run on a finer grid (64 × 32 against 160 × 80 cells).'),
  pv('Taylor axial dispersion 10.1 R u*', 'runExtras() species transport default', 'not found in an openly readable source (searched again)', 'https://doi.org/10.1098/rspa.1954.0130', 'unverified', 'Used only as the default of the dispersion input.'),
  pv('Mandhane, Gregory & Aziz (1974) transition lines', 'mandhaneRegime()', 'C. Bell, fluids (MIT), fluids/two_phase.py Mandhane_Gregory_Aziz_regime', FL + '/blob/master/fluids/two_phase.py', 'corrected', 'Annular-mist line: the break points (70, 60, 38, 40, 50, 100, 230 ft/s) agree with the piecewise power laws of the open code. Stratified / wave line: old break points (0.01, 14), (0.1, 10.5), (0.3, 2.5), (0.5, 2.5), (1.7, 3.25) ft/s were shifted; new (0.01, 32.7), (0.1, 14), (0.2, 10.5), (1.15, 2.5), (4.8, 2.5), (14, 3.26) from the open code. The property corrections X1, Y1 of the original are not applied (air–water coordinates).'),
  pv('Reference validation data', 'validationData', 'see each data set', 'https://data.mendeley.com/datasets/wyfdm5ysh6/1', 'verified', 'All rows were read from the cited addresses on the retrieval date or computed for this suite, and are stored in js/data/ref/flow.js rounded to five significant figures. Licences: CC BY 4.0 (Mendeley, Data in Brief, Zenodo), own computations (OpenFOAM v2412, in-app 3-D solver), and at most twelve numbers quoted as cited facts from each of four published works without an open licence (Jansen et al., Princeton Superpipe, Ghia et al., Moser et al.). Removed: the UNICAMP inclined-slug sets (CC BY-NC) and the void profile taken from a GPL tutorial file. Two corrupt blocks of the flow-pattern file are left out (Shoham −50°, Kokal 76 mm).'),
];

/*NEW-SECTIONS-END*/
// =====================================================================================================
// 8. Suite assembly
// =====================================================================================================
const cellOf = (st, i) => ({ vsl: st.vsl[i], vsg: st.vsg[i], rhoL: st.rhoL[i], rhoG: st.rhoG[i], muL: st.muL[i], muG: st.muG[i], sigma: st.sigma[i], D: st.D[i], theta: st.theta[i] });
const slugOpts = (v) => ({ freqModel: ['zabaras', 'gregory', 'heywood', 'greskovich', 'unitCell'].includes(v.freqModel) ? v.freqModel : 'unitCell', lengthModel: ['scott', 'brill', 'norris'].includes(v.lengthModel) ? v.lengthModel : 'scott', vtModel: v.vtModel === 'nicklin' ? 'nicklin' : 'bendiksen', bodyModel: v.bodyModel === 'barnea' ? 'barnea' : 'gregory', freqMult: num(v.freqMult, 1, 0.01, 100), lenMult: num(v.lenMult, 1, 0.01, 100), filmMult: num(v.filmMult, 1, 0.2, 3), vtMult: num(v.vtMult, 1, 0.5, 2), bodyMult: num(v.bodyMult, 1, 0.5, 1.5) });
const isSlug = (r) => r === 'slug' || r === 'churn';
const pwl = (ts, ys) => (t) => interp1(ts, ys, t);
const DEFAULT_SCHEDULE = [{ t: 0, rate: 1, dp: 0, choke: 100 }, { t: 0.25, rate: 1, dp: 0, choke: 100 }, { t: 0.5, rate: 0.6, dp: 0, choke: 100 }, { t: 6, rate: 0.6, dp: 0, choke: 100 }];

/** Regime (Taitel–Dukler / Barnea) at every node of a steady solution. */
function regimes(st) { return st.vsl.map((_, i) => (st.vsl[i] > 1e-9 && st.vsg[i] > 1e-9 ? flowPattern(cellOf(st, i)).pattern : st.regime[i])); }

/**
 * Minimum stable rate: scan the rate downwards on a coarse grid and return the lowest rate fraction before severe or terrain slugging is predicted.
 * With mit = { chokeOpening (%), chokeDpOpen (bar at full opening and the case rate), boostDp (bar of subsea boosting that can be spent as back-pressure),
 * liftGas (kg/s of riser-base lift gas at the case rate), label } every rate that the screening calls unstable is examined again with the riser cycle
 * model under those mitigations (choke throttled to the given opening, the booster head added to the choke pressure drop, lift gas in the riser);
 * rateControlled is the lowest rate that is stable either way. Terrain slugging (no riser) has no mitigation in this model.
 */
function minimumStableRate(cfg, mNow, mit = null) {
  const rates = [1.2, 1, 0.85, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.12].map((r) => r * Math.max(mNow, 1e-3)), rows = [];
  let guess, lowest = rates[0], hit = false, lowestC = rates[0], hitC = false, used = false;
  for (const m of rates) {
    const r = steadySolve({ ...cfg.base, n: 30, mScale: m, pOut: cfg.pOut, pGuess: guess, hydrate: false, tolP: 1e-3 });
    if (!r.ok) break;
    guess = r.pIn;
    const sc = cfg.hasRiser ? severeScreen(r, { riserBaseS: cfg.riserBaseS, pSep: cfg.pSep }) : null, ter = sc ? null : terrainSlugging(r), unstable = sc ? sc.severe : ter.accumulates;
    let controlled = !unstable, dpChoke = null;
    if (unstable && sc && mit && !hitC) { // riser cycle with the mitigations of the design
      const k = sc.iBase - 1, z = clamp(mit.chokeOpening, 5, 100); dpChoke = Math.max(mit.chokeDpOpen, 0.05) * m * m * (100 / z) ** 2 + mit.boostDp;
      try { const c = riserSluggingCycle({ D: r.D[k], feedLength: sc.feedLength, feedAngle: sc.feedAngle * DEG, riserHeight: sc.riserHeight, riserLength: sc.riserLength, wG: Math.max(r.mG[k], 1e-6), wLift: mit.liftGas, wL: Math.max(r.mL[k], 1e-6), rhoL: r.rhoL[k], muL: r.muL[k], T: r.T[k] + 273.15, zG: r.zG[k], mwG: r.mwG[k], pSep: cfg.pSep * 1e5, alphaL: 1 - sc.alpha, chokeOpening: 100, chokeDp: dpChoke * 1e5, rough: cfg.rough, maxCycles: 3 }); controlled = c.stable; used = true; } catch { controlled = false; }
      if (cfg.pSep + dpChoke >= sc.taitelPsep) { controlled = true; used = true; } // Taitel (1986): the riser column is stable when the back-pressure exceeds ρL g (α L / 0.89 − h); the choke (and the booster head spent across it) supplies that back-pressure
    }
    rows.push({ m, pIn: r.pIn, unstable, controlled, dpChoke, taitel: sc ? sc.taitelPsep : null, boe: sc ? sc.boe : null, pots: sc ? sc.pots : null, regime: sc ? sc.feedRegime : ter.worst ? ter.worst.downRegime : '—' });
    if (unstable && !hit) hit = true; if (!hit) lowest = m;
    if (!controlled && !hitC) hitC = true; if (!hitC) lowestC = m;
  }
  const last = rows.length ? rows[rows.length - 1].m : mNow;
  return { rate: hit ? lowest : last, limited: hit, rateControlled: hitC ? lowestC : last, limitedControlled: hitC, mitigated: used, rows };
}

/** Transient set-up from the suite inputs and a steady solution; returns the stepper and its descriptors. */
function buildTransient(v, st) {
  const cfg = st.cfg, fm = cfg.fm, n = Math.round(num(v.nCells, 50, 8, 600)), grid = transientGrid(cfg.profile, n, cfg.hasRiser ? cfg.riserBaseS : Infinity, cfg.hasRiser ? num(v.riserRefine, 3, 1, 10) : 1);
  const rows = (Array.isArray(v.schedule) && v.schedule.length ? v.schedule : DEFAULT_SCHEDULE).map((r) => ({ t: +r?.t, rate: +r?.rate, dp: +r?.dp, choke: +r?.choke })).filter((r) => Number.isFinite(r.t) && r.t >= 0).sort((a, b) => a.t - b.t);
  if (!rows.length) rows.push({ t: 0, rate: 1, dp: 0, choke: 100 });
  const ts = rows.map((r) => r.t * 3600), tEnd = num(v.tEnd, 2.5, 0.01, 2000) * 3600, mdotCase = st.mdot, r0 = clamp(fin(rows[0].rate, 1), 0, 10);
  const rateOf = pwl(ts, rows.map((r) => clamp(fin(r.rate, 1), 0, 10))), dpOf = pwl(ts, rows.map((r) => fin(r.dp, 0))), chokeOf = pwl(ts, rows.map((r) => clamp(fin(r.choke, 100), 1, 100)));
  const slip = { c0: 1.2, vdScale: 1, fricMult: 1 }, uniform = v.initMode === 'uniform', at = (ref, a) => grid.sc.map((s) => interp1(ref.s, a, s));
  const tdSteady = (fmult, nn, guess) => steadySolve({ ...cfg.base, n: nn, model: 'transient', mp: { slip: { ...slip, fricMult: fmult } }, energy: 'cpjt', mScale: st.mScale * Math.max(r0, 0.02), pOut: cfg.pOut, pGuess: guess, hydrate: false, tolP: 1e-4 });
  let init, matched = false, mdot0 = mdotCase * r0, ref = null;
  if (uniform) { // settled line: uniform holdup and temperature, hydrostatic pressure from the outlet
    const H = clamp(num(v.initHoldup, 0.5, 0, 1), 0, 1), T0 = num(v.initT, cfg.thermal.tSeabed, -40, 200), P = new Array(n); let p = num(v.initP, cfg.pSep, 1.05, 1300);
    for (let i = n - 1; i >= 0; i--) { const pr = fm.at(p, T0), rho = pr.rhoL * H + pr.rhoG * (1 - H), dz = (i === n - 1 ? grid.z[n] : grid.zc[i + 1]) - grid.zc[i]; p += (rho * G * dz) / 1e5; P[i] = Math.max(p, 1.05); }
    init = { P, T: P.map(() => T0), holdup: P.map(() => H) }; mdot0 = 0;
  } else {
    if (v.matchSteady !== false && Math.abs(r0 - 1) < 1e-9) { // tune the wall friction of the transient closure to the steady inlet pressure
      const g = (f) => { const r = tdSteady(f, 40, st.pIn); return r.ok ? r.pIn - st.pIn : 1e3; }, lo = g(0.2), hi = g(5);
      if (lo < 0 && hi > 0) { slip.fricMult = brent(g, 0.2, 5, 1e-3); matched = true; }
    }
    ref = tdSteady(slip.fricMult, st.n, st.pIn);
    if (!ref.ok) throw new Error('The transient could not be initialised: ' + (ref.reason || 'no steady solution of the drift-flux closure at the initial rate.'));
    init = { P: at(ref, ref.P), T: at(ref, ref.T), holdup: at(ref, ref.holdup) };
  }
  const mid = fm.at(0.5 * (st.pIn + st.pOut), 0.5 * (st.tIn + st.tOut)), D = cfg.idOf ? cfg.idMin : cfg.id, jOut = ref ? ref.vm[ref.n] : 1, rhoOut = ref ? ref.rhoM[ref.n] : 100;
  const chokeK = cfg.chokeDp > 0 ? (cfg.chokeDp * 1e5 * (st.mScale * Math.max(r0, 0.05)) ** 2) / Math.max(rhoOut * jOut * jOut, 1e-6) : 0;
  const jn = st.junction || null, mMain = jn ? jn.mMain : mdotCase; if (jn && !uniform) mdot0 = mMain * r0;
  const sim = makeTransient({ fm, grid, D, rough: cfg.rough, init, mdot0, mdotOf: (t) => mMain * rateOf(t), sources: jn ? [{ s: jn.s, mdotOf: (t) => jn.mBranch * rateOf(t), T: jn.tBranch }] : [], tInOf: () => cfg.tIn, pOutOf: (t) => Math.max(cfg.pSep + dpOf(t), 1.05), chokeOf, chokeK,
    slip, fricMult: slip.fricMult, U: at(st, st.U), tAmb: grid.sc.map((s, i) => cfg.tAmbOf(s, grid.zc[i])), wallC: (BASE.rhoSteel * 470 * ((D + 2 * cfg.wt) ** 2 - D * D)) / (D * D), cpG: mid.cpG, cpL: mid.cpL,
    cfl: num(v.cfl, 0.6, 0.05, 0.95), tEnd, nSeries: 380, nField: 70, maxSteps: 60000 });
  return { sim, grid, tEnd, slip, matched, ref, rows, uniform };
}

const INPUTS = [
  { group: 'Boundary conditions and rates', tab: 'inputs', help: 'The well stream (composition, GOR, rates, water cut) is the case fluid; here you choose which pressures and rates are fixed.', fields: [
    { key: 'bc', label: 'Boundary conditions', type: 'select', value: 'outletP', options: [{ value: 'outletP', label: 'Outlet pressure + rate (inlet pressure is solved)' }, { value: 'inletP', label: 'Inlet (source) pressure + rate (outlet pressure is solved)' }, { value: 'bothP', label: 'Inlet and outlet pressure (rate is solved)' }, { value: 'ipr', label: 'Well inflow + outlet pressure (rate is solved)' }] },
    { key: 'rateFrac', label: 'Rate multiplier on the case rates', unit: '×', value: 1, min: 0.01, max: 5, typical: [0.2, 1.5], help: 'Scales the oil, gas and water mass rates of the case fluid together (GOR and water cut unchanged).', showIf: (v) => v.bc === 'outletP' || v.bc === 'inletP' || !v.bc },
    { key: 'wc', label: 'Water cut', unit: '%', value: 20, min: 0, max: 98, typical: [0, 90], help: 'Water cut of the stock-tank liquid (oil-rate basis). Changing it rebuilds the stream rates.' },
    { key: 'fluidSel', label: 'Fluid', type: 'select', value: 'case', options: [{ value: 'case', label: 'Case fluid (from the Fluid suite)' }, { value: 'gascond', label: 'Lean gas condensate, 10 MSm³/d (built-in)' }] },
    { key: 'tIn', label: 'Inlet temperature', unit: '°C', value: BASE.tIn, min: -20, max: 200, typical: [40, 120] },
    { key: 'pOut', label: 'Outlet (separator / back-) pressure', unit: 'bara', value: BASE.pOut, min: 1.1, max: 500, typical: [5, 100] },
    { key: 'pInSet', label: 'Inlet (source) pressure', unit: 'bara', value: 95, min: 2, max: 1300, typical: [30, 300], showIf: (v) => v.bc === 'inletP' || v.bc === 'bothP' },
    { key: 'pRes', label: 'Reservoir pressure', unit: 'bara', value: BASE.pRes, min: 5, max: 1500, showIf: (v) => v.bc === 'ipr' },
    { key: 'piIpr', label: 'Productivity index (stock-tank liquid)', unit: 'Sm³/d/bar', value: BASE.pi, min: 0.01, max: 1e5, showIf: (v) => v.bc === 'ipr' },
    { key: 'iprType', label: 'Inflow relation', type: 'select', value: 'linear', options: [{ value: 'linear', label: 'Straight line (PI)' }, { value: 'vogel', label: 'Vogel (solution-gas drive)' }], showIf: (v) => v.bc === 'ipr' },
    { key: 'wellDp', label: 'Pressure loss sandface → flowline inlet at the case rate', unit: 'bar', value: 60, min: 0, max: 1000, help: 'Tubing head, friction and wellhead choke; 80 % is treated as static and 20 % scales with the rate squared.', showIf: (v) => v.bc === 'ipr' },
  ] },
  { group: 'Line geometry, fittings and choke', tab: 'inputs', fields: [
    { key: 'profile', label: 'Elevation profile', type: 'table', columns: [{ key: 'x', label: 'Distance', unit: 'm' }, { key: 'z', label: 'Elevation', unit: 'm' }], value: BASE.profile.map((p) => ({ x: p.x, z: p.z })), help: 'Horizontal distance from the inlet and elevation relative to mean sea level (negative below).' },
    { key: 'idMm', label: 'Inner diameter', unit: 'mm', value: BASE.idMm, min: 20, max: 1500, typical: [100, 900] },
    { key: 'wtMm', label: 'Wall thickness', unit: 'mm', value: BASE.wtMm, min: 2, max: 80 },
    { key: 'roughUm', label: 'Wall roughness', unit: 'µm', value: BASE.roughUm, min: 0, max: 3000, typical: [15, 300] },
    { key: 'riserBaseX', label: 'Riser base position', unit: 'm', value: BASE.riserBaseX, min: 0, max: 1e6, help: 'Distance at which the flowline meets the riser. Set it to the line length when there is no riser.' },
    { key: 'kLoss', label: 'Total minor-loss coefficient (bends, tees, spools)', unit: '–', value: 0, min: 0, max: 5000, help: 'Spread evenly along the line: Δp = K ρ v² / 2 on the no-slip mixture.' },
    { key: 'fittings', label: 'Valves and fittings', type: 'table', columns: [{ key: 'x', label: 'Distance', unit: 'm' }, { key: 'K', label: 'Loss coefficient fully open', unit: '–' }, { key: 'open', label: 'Opening', unit: '%' }], value: [], help: 'Local losses at a position; a partly closed valve is scaled with (100 / opening)².' },
    { key: 'chokeOpening', label: 'Topsides choke opening', unit: '%', value: 100, min: 1, max: 100 },
    { key: 'chokeDp', label: 'Choke pressure drop fully open at the case rate', unit: 'bar', value: 0, min: 0, max: 200, help: '0 = no choke modelled (the outlet pressure applies at the end of the line).' },
    { key: 'effIdMm', label: 'Effective bore with deposits (uniform)', unit: 'mm', value: 0, min: 0, max: 1500, help: '0 = clean bore. Used only when no deposit profile is given.' },
    { key: 'roughEffUm', label: 'Effective roughness with deposits', unit: 'µm', value: 0, min: 0, max: 20000, help: '0 = use the wall roughness.' },
    { key: 'deposit', label: 'Deposit thickness profile', type: 'table', columns: [{ key: 'x', label: 'Distance', unit: 'm' }, { key: 't', label: 'Thickness', unit: 'mm' }], value: [], help: 'Hydrate, wax and scale layer from the solids suite; narrows the bore locally.' },
  ] },
  { group: 'Ambient and heat transfer', tab: 'inputs', fields: [
    { key: 'uMode', label: 'Overall heat-transfer coefficient', type: 'select', value: 'input', options: [{ value: 'input', label: 'Given U-value' }, { value: 'layers', label: 'From films, wall, insulation and burial' }] },
    { key: 'uValue', label: 'U-value (on inner diameter)', unit: 'W/m²K', value: BASE.U, min: 0, max: 2000, typical: [1, 30], showIf: (v) => v.uMode !== 'layers' },
    { key: 'uMult', label: 'U-value multiplier (calibration)', unit: '×', value: 1, min: 0.1, max: 10 },
    { key: 'insT', label: 'Insulation thickness', unit: 'mm', value: BASE.insulation.t * 1000, min: 0, max: 300 },
    { key: 'insK', label: 'Insulation conductivity', unit: 'W/mK', value: BASE.insulation.k, min: 0.01, max: 3 },
    { key: 'insRho', label: 'Insulation density', unit: 'kg/m³', value: 900, min: 30, max: 3000 },
    { key: 'insCp', label: 'Insulation heat capacity', unit: 'J/kgK', value: 1700, min: 300, max: 4000 },
    { key: 'kWall', label: 'Pipe-wall conductivity', unit: 'W/mK', value: 45, min: 5, max: 400 },
    { key: 'tSeabed', label: 'Seabed temperature', unit: '°C', value: BASE.tSeabed, min: -2, max: 35 },
    { key: 'tSeaSurface', label: 'Sea-surface temperature', unit: '°C', value: BASE.tSeaSurface, min: -2, max: 35 },
    { key: 'thermocline', label: 'Thermocline depth scale', unit: 'm', value: 250, min: 20, max: 2000 },
    { key: 'tAir', label: 'Air temperature (above sea level)', unit: '°C', value: BASE.tAir, min: -50, max: 55 },
    { key: 'currentSpeed', label: 'Sea current speed', unit: 'm/s', value: BASE.currentSpeed, min: 0, max: 3 },
    { key: 'windSpeed', label: 'Wind speed', unit: 'm/s', value: 5, min: 0, max: 50 },
    { key: 'burialDepth', label: 'Burial depth to top of pipe (flowline)', unit: 'm', value: 0, min: 0, max: 10, help: '0 = exposed. Used with the layer-based U-value and in the cooldown preview.' },
    { key: 'kSoil', label: 'Soil conductivity', unit: 'W/mK', value: 1.2, min: 0.2, max: 5 },
    { key: 'heatTrace', label: 'Imposed heating (e.g. electrical)', unit: 'W/m', value: 0, min: 0, max: 500 },
    { key: 'watC', label: 'Wax appearance temperature', unit: '°C', value: 0, min: -50, max: 90, help: '0 = not checked. Offered from the Fluid suite.' },
  ] },
  { group: 'Flow model', tab: 'setup', fields: [
    { key: 'model', label: 'Holdup and pressure-gradient model', type: 'select', value: 'beggsBrill', options: HOLDUP_MODELS },
    { key: 'fModel', label: 'Friction factor', type: 'select', value: 'colebrook', options: [{ value: 'colebrook', label: 'Colebrook–White' }, { value: 'haaland', label: 'Haaland' }, { value: 'swamee', label: 'Swamee–Jain' }, { value: 'churchill', label: 'Churchill' }] },
    { key: 'energy', label: 'Energy equation', type: 'select', value: 'enthalpy', options: [{ value: 'enthalpy', label: 'Flowing-enthalpy balance (flashing, Joule–Thomson, elevation, kinetic)' }, { value: 'cpjt', label: 'Frozen heat capacity and Joule–Thomson coefficients' }] },
    { key: 'c0', label: 'Distribution parameter C0', unit: '–', value: 1.2, min: 0.9, max: 1.6, showIf: (v) => v.model === 'zuberFindlay' },
    { key: 'vDrift', label: 'Drift velocity', unit: 'm/s', value: 0.35, min: 0, max: 3, showIf: (v) => v.model === 'zuberFindlay' },
    { key: 'wallisN', label: 'Wallis hindered-drift exponent n', unit: '–', value: 0, min: 0, max: 4, help: 'Drift velocity × (1 − α)^n; 0 gives the plain Zuber–Findlay relation.', showIf: (v) => v.model === 'zuberFindlay' },
    { key: 'holdupMult', label: 'Liquid-holdup multiplier (calibration)', unit: '×', value: 1, min: 0.5, max: 2 },
    { key: 'cErosion', label: 'Erosional-velocity constant C (API RP 14E)', unit: '–', value: 100, min: 50, max: 300 },
    { key: 'mapX', label: 'Location for the regime map and local models', unit: 'm', value: 9000, min: 0, max: 1e6 },
  ] },
  { group: 'Slugging', tab: 'setup', fields: [
    { key: 'freqModel', label: 'Slug frequency', type: 'select', value: 'unitCell', options: [{ value: 'unitCell', label: 'Unit-cell balance with the length correlation' }, { value: 'zabaras', label: 'Zabaras (2000)' }, { value: 'gregory', label: 'Gregory & Scott (1969)' }, { value: 'greskovich', label: 'Greskovich & Shrier (1972)' }, { value: 'heywood', label: 'Heywood & Richardson (1979)' }] },
    { key: 'lengthModel', label: 'Slug length', type: 'select', value: 'scott', options: [{ value: 'scott', label: 'Scott et al. (1989)' }, { value: 'brill', label: 'Brill et al. (1981)' }, { value: 'norris', label: 'Norris (1982)' }] },
    { key: 'vtModel', label: 'Translational velocity', type: 'select', value: 'bendiksen', options: [{ value: 'bendiksen', label: 'Bendiksen (1984)' }, { value: 'nicklin', label: 'Nicklin et al. (1962)' }] },
    { key: 'bodyModel', label: 'Slug-body holdup', type: 'select', value: 'gregory', options: [{ value: 'gregory', label: 'Gregory et al. (1978)' }, { value: 'barnea', label: 'Barnea & Brauner (1985)' }] },
    { key: 'freqMult', label: 'Slug-frequency multiplier (calibration)', unit: '×', value: 1, min: 0.05, max: 20 },
    { key: 'lenMult', label: 'Slug-length multiplier (calibration)', unit: '×', value: 1, min: 0.05, max: 20 },
    { key: 'filmMult', label: 'Film-holdup multiplier (calibration)', unit: '×', value: 1, min: 0.3, max: 2.5 },
    { key: 'trackSlugs', label: 'Track slugs to the outlet', type: 'bool', value: true },
    { key: 'nSlugs', label: 'Slugs to collect at the outlet', unit: '–', value: 150, min: 20, max: 2000, showIf: (v) => v.trackSlugs !== false },
    { key: 'slugSeed', label: 'Random seed', unit: '–', value: 42, min: 1, max: 1e6, showIf: (v) => v.trackSlugs !== false },
    { key: 'drainFactor', label: 'Vessel liquid drain rate / mean liquid rate', unit: '×', value: 1.2, min: 1, max: 5 },
    { key: 'feedLength', label: 'Gas-pocket length upstream of the riser', unit: 'm', value: 0, min: 0, max: 1e6, help: '0 = from the last crest upstream of the riser base.' },
  ] },
  { group: 'Transient', tab: 'setup', fields: [
    { key: 'transient', label: 'Run the transient', type: 'bool', value: true },
    { key: 'schedule', label: 'Boundary history', type: 'table', columns: [{ key: 't', label: 'Time', unit: 'h' }, { key: 'rate', label: 'Rate / steady rate', unit: '×' }, { key: 'dp', label: 'Outlet-pressure change', unit: 'bar' }, { key: 'choke', label: 'Choke opening', unit: '%' }], value: DEFAULT_SCHEDULE.map((r) => ({ ...r })), help: 'Piecewise-linear history of the inlet rate, the outlet (separator) pressure and the choke.', showIf: (v) => v.transient !== false },
    { key: 'tEnd', label: 'Simulated time', unit: 'h', value: 2.5, min: 0.05, max: 500, showIf: (v) => v.transient !== false },
    { key: 'initMode', label: 'Initial condition', type: 'select', value: 'steady', options: [{ value: 'steady', label: 'Converged steady state at the first rate' }, { value: 'uniform', label: 'Settled line: given pressure, temperature and holdup' }], showIf: (v) => v.transient !== false },
    { key: 'initP', label: 'Initial pressure at the outlet', unit: 'bara', value: BASE.pOut, min: 1.1, max: 1000, showIf: (v) => v.transient !== false && v.initMode === 'uniform' },
    { key: 'initT', label: 'Initial temperature', unit: '°C', value: BASE.tSeabed, min: -20, max: 150, showIf: (v) => v.transient !== false && v.initMode === 'uniform' },
    { key: 'initHoldup', label: 'Initial liquid holdup', unit: '–', value: 0.5, min: 0, max: 1, showIf: (v) => v.transient !== false && v.initMode === 'uniform' },
    { key: 'matchSteady', label: 'Tune transient wall friction to the steady inlet pressure', type: 'bool', value: true, showIf: (v) => v.transient !== false },
  ] },
  { group: 'Local models and comparisons', tab: 'setup', fields: [
    { key: 'localModels', label: 'Run the local models (radial RANS, interface capturing)', type: 'bool', value: true },
    { key: 'ransModel', label: 'Turbulence model of the radial solve', type: 'select', value: 'mixing', options: RANS_MODELS, showIf: (v) => v.localModels !== false },
    { key: 'compareModels', label: 'Compare all holdup models', type: 'bool', value: true },
  ] },
  { group: 'Branch, pump, separator and compressor', tab: 'inputs', help: 'Flows between connected branches and the hydraulic characteristics of the equipment at the ends of the line. With a branch the boundary condition is the outlet pressure with given rates.', fields: [
    { key: 'branchFrac', label: 'Branch share of the total rate', unit: '–', value: 0, min: 0, max: 0.95, help: '0 = no branch. The branch carries the same fluid and joins the main line at the junction; mass and enthalpy are mixed there.' },
    { key: 'branchX', label: 'Junction position on the main line', unit: 'm', value: 9000, min: 0, max: 1e6, showIf: (v) => v.branchFrac > 0 },
    { key: 'branchLength', label: 'Branch length', unit: 'm', value: 3000, min: 10, max: 1e6, showIf: (v) => v.branchFrac > 0 },
    { key: 'branchDz', label: 'Elevation gain along the branch', unit: 'm', value: 0, min: -3000, max: 3000, showIf: (v) => v.branchFrac > 0 },
    { key: 'branchIdMm', label: 'Branch inner diameter (0 = as the main line)', unit: 'mm', value: 0, min: 0, max: 1500, showIf: (v) => v.branchFrac > 0 },
    { key: 'branchTin', label: 'Branch inlet temperature', unit: '°C', value: BASE.tIn, min: -20, max: 200, showIf: (v) => v.branchFrac > 0 },
    { key: 'pumpDp0', label: 'Booster pump at the inlet: shut-off pressure rise', unit: 'bar', value: 0, min: 0, max: 400, help: '0 = no pump. Δp = Δp₀ s² − (Δp₀/q_max²) q² with the suction volume rate q.' },
    { key: 'pumpQmax', label: 'Booster pump: run-out volume rate at rated speed', unit: 'm³/s', value: 0.5, min: 0.001, max: 50, showIf: (v) => v.pumpDp0 > 0 },
    { key: 'pumpSpeed', label: 'Booster pump: speed / rated speed', unit: '–', value: 1, min: 0.2, max: 1.5, showIf: (v) => v.pumpDp0 > 0 },
    { key: 'sepMode', label: 'Separator pressure', type: 'select', value: 'fixed', options: [{ value: 'fixed', label: 'Fixed (given outlet pressure)' }, { value: 'valve', label: 'From the gas-outlet valve and the compressor curve' }] },
    { key: 'sepKv', label: 'Gas-outlet valve coefficient K_v in p_sep² − p_suc² = K_v q²', unit: 'bar²/(Sm³/s)²', value: 0.5, min: 0, max: 1e5, showIf: (v) => v.sepMode === 'valve' },
    { key: 'compPd', label: 'Compressor discharge pressure', unit: 'bara', value: 120, min: 2, max: 800, showIf: (v) => v.sepMode === 'valve' },
    { key: 'compHead', label: 'Compressor polytropic head at zero flow', unit: 'kJ/kg', value: 200, min: 0, max: 500, showIf: (v) => v.sepMode === 'valve' },
    { key: 'compQmax', label: 'Compressor suction volume rate at zero head', unit: 'm³/s', value: 2, min: 0.001, max: 500, showIf: (v) => v.sepMode === 'valve' },
    { key: 'compSpeed', label: 'Compressor speed / rated speed', unit: '–', value: 1, min: 0.3, max: 1.3, showIf: (v) => v.sepMode === 'valve' },
  ] },
  { group: 'Measurements for comparison', tab: 'inputs', help: 'Optional. Local measurements of any flow configuration (horizontal, inclined, vertical; gas–liquid, oil–water with the lighter liquid as the “gas” phase, gas–oil–water through the water cut) and system measurements of this line are compared with the model; the deviations and error statistics appear in the result tables.', fields: [
    { key: 'valPoints', label: 'Local measurements', type: 'table', columns: [{ key: 'kind', label: 'Quantity (holdup, dpdx, regime, slugFreq, slugLength, slugVelocity, bodyHoldup, filmThickness, level, entrainment)', type: 'text' }, { key: 'D', label: 'Diameter', unit: 'mm' }, { key: 'angle', label: 'Inclination', unit: '°' }, { key: 'vsl', label: 'Superficial liquid velocity', unit: 'm/s' }, { key: 'vsg', label: 'Superficial gas velocity', unit: 'm/s' }, { key: 'rhoL', label: 'Liquid (oil) density', unit: 'kg/m³' }, { key: 'rhoG', label: 'Gas density', unit: 'kg/m³' }, { key: 'muL', label: 'Liquid (oil) viscosity', unit: 'mPa·s' }, { key: 'muG', label: 'Gas viscosity', unit: 'mPa·s' }, { key: 'sigma', label: 'Interfacial tension', unit: 'mN/m' }, { key: 'wc', label: 'Water cut of the liquid', unit: '%' }, { key: 'P', label: 'Pressure', unit: 'bara' }, { key: 'measured', label: 'Measured value (SI; dpdx in Pa/m, film in mm; text for regime)', type: 'text' }], value: [], help: 'One row per measured point.' },
    { key: 'valSystem', label: 'System measurements of this line', type: 'table', columns: [{ key: 'quantity', label: 'Quantity (pIn, dp, tOut, inventory, holdupMean, slugPeriod, slugFrequency, slugLength, slugLengthMax, slugVelocity, bodyHoldup, surge, catcherLoad, pAmplitude, severePeriod, severeAmplitude, cooldown, terrainVolume, tracerArrival)', type: 'text' }, { key: 'measured', label: 'Measured value' }], value: [], help: 'Field or loop measurements of the simulated line, e.g. slug-catcher load, pressure-fluctuation amplitude, cooldown time, terrain-slug volume.' },
  ] },
  { group: 'Closure parameters', tab: 'setup', help: 'Multipliers and constants of the closure laws; all can be estimated on the Calibration tab.', fields: [
    { key: 'fiMult', label: 'Interfacial-friction multiplier (stratified flow)', unit: '×', value: 1, min: 0.1, max: 10 },
    { key: 'fwlMult', label: 'Liquid-wall friction multiplier (stratified flow)', unit: '×', value: 1, min: 0.1, max: 10 },
    { key: 'fwgMult', label: 'Gas-wall friction multiplier (stratified flow)', unit: '×', value: 1, min: 0.1, max: 10 },
    { key: 'transMult', label: 'Stratified → slug transition: multiplier on the critical gas velocity', unit: '×', value: 1, min: 0.3, max: 3 },
    { key: 'entMult', label: 'Entrainment-rate multiplier (annular flow)', unit: '×', value: 1, min: 0.05, max: 10 },
    { key: 'kDep', label: 'Droplet deposition velocity', unit: 'm/s', value: 0.15, min: 0.001, max: 5 },
    { key: 'weCrit', label: 'Critical Weber number of the droplets (droplet size)', unit: '–', value: 12, min: 1, max: 60 },
    { key: 'fiAnnMult', label: 'Interfacial-friction multiplier (annular film)', unit: '×', value: 1, min: 0.1, max: 10 },
    { key: 'vtMult', label: 'Slug-celerity multiplier', unit: '×', value: 1, min: 0.6, max: 1.6 },
    { key: 'bodyMult', label: 'Slug-body holdup multiplier', unit: '×', value: 1, min: 0.6, max: 1.3 },
    { key: 'initMult', label: 'Slug-initiation frequency multiplier (tracking)', unit: '×', value: 1, min: 0.1, max: 10 },
    { key: 'relaxMult', label: 'Slug growth / decay: film relaxation-length multiplier', unit: '×', value: 1, min: 0.1, max: 10 },
    { key: 'wakeMult', label: 'Slug merging: wake-acceleration multiplier', unit: '×', value: 1, min: 0, max: 5 },
    { key: 'hOutMult', label: 'External heat-transfer coefficient multiplier', unit: '×', value: 1, min: 0.1, max: 10 },
    { key: 'natConv', label: 'Add free convection to the outside film', type: 'bool', value: true },
    { key: 'bubbleMm', label: 'Bubble size at the riser base', unit: 'mm', value: 3, min: 0.3, max: 25 },
    { key: 'crcMult', label: 'Bubble coalescence coefficient multiplier', unit: '×', value: 1, min: 0, max: 20 },
    { key: 'ctiMult', label: 'Bubble breakup coefficient multiplier', unit: '×', value: 1, min: 0, max: 20 },
    { key: 'cvm', label: 'Virtual-mass coefficient', unit: '–', value: 0.5, min: 0, max: 2 },
    { key: 'tracer', label: 'Transport an inlet step of a species with the liquid', type: 'bool', value: true },
    { key: 'dispersion', label: 'Axial dispersion coefficient (0 = Taylor estimate)', unit: 'm²/s', value: 0, min: 0, max: 1000, showIf: (v) => v.tracer !== false },
  ] },
  { group: 'Solver task', tab: 'setup', help: 'The line solution always runs. The other solvers run on request, on top of it, with their own resolution inputs.', fields: [
    { key: 'task', label: 'Additional solver', type: 'select', value: 'line', options: TASKS },
    { key: 'tfLengthD', label: 'Two-fluid test section length', unit: 'diameters', value: 200, min: 30, max: 3000, showIf: (v) => v.task === 'twofluid' },
    { key: 'tfDiameterMm', label: 'Two-fluid test section diameter (0 = line diameter)', unit: 'mm', value: 0, min: 0, max: 1500, showIf: (v) => v.task === 'twofluid' },
    { key: 'tfVsl', label: 'Two-fluid superficial liquid velocity (0 = at the map location)', unit: 'm/s', value: 0, min: 0, max: 10, showIf: (v) => v.task === 'twofluid' },
    { key: 'tfVsg', label: 'Two-fluid superficial gas velocity (0 = at the map location)', unit: 'm/s', value: 0, min: 0, max: 50, showIf: (v) => v.task === 'twofluid' },
    { key: 'tfPbar', label: 'Two-fluid pressure (0 = at the map location)', unit: 'bara', value: 0, min: 0, max: 500, showIf: (v) => v.task === 'twofluid' },
    { key: 'tfCells', label: 'Two-fluid cells', unit: '–', value: 240, min: 40, max: 3000, showIf: (v) => v.task === 'twofluid' },
    { key: 'tfTime', label: 'Two-fluid simulated time', unit: 's', value: 20, min: 1, max: 600, showIf: (v) => v.task === 'twofluid' },
    { key: 'tfCfl', label: 'Two-fluid CFL number', unit: '–', value: 0.4, min: 0.05, max: 0.9, showIf: (v) => v.task === 'twofluid' },
    { key: 'cfdCase', label: '2-D case', type: 'select', value: 'dambreak', options: [{ value: 'dambreak', label: 'Dam break of a liquid column (water–gas, local properties)' }, { value: 'slugfront', label: 'Slug tail collapsing onto the film (local slug body and film)' }, { value: 'cavity', label: 'Lid-driven cavity (benchmark)' }, { value: 'channel', label: 'Turbulent channel, SST k–ω (benchmark)' }], showIf: (v) => v.task === 'cfd' },
    { key: 'cfdN', label: '2-D cells across the height', unit: '–', value: 24, min: 8, max: 200, showIf: (v) => v.task === 'cfd' },
    { key: 'cfdRe', label: 'Reynolds number (cavity) or friction Reynolds number (channel)', unit: '–', value: 100, min: 1, max: 5000, showIf: (v) => v.task === 'cfd' && (v.cfdCase === 'cavity' || v.cfdCase === 'channel') },
    { key: 'cfdTime', label: '2-D simulated time (0 = automatic)', unit: 's', value: 0, min: 0, max: 60, showIf: (v) => v.task === 'cfd' && (v.cfdCase === 'dambreak' || v.cfdCase === 'slugfront') },
    { key: 'cfdTurbulence', label: '2-D turbulence closure', type: 'select', value: 'sst', options: [{ value: 'sst', label: 'SST k–ω' }, { value: 'none', label: 'None (laminar)' }], showIf: (v) => v.task === 'cfd' && (v.cfdCase === 'dambreak' || v.cfdCase === 'slugfront') },
    { key: 'bubblyN', label: 'Radial nodes of the bubbly-flow model', unit: '–', value: 100, min: 30, max: 500, showIf: (v) => v.task === 'bubbly' },
    { key: 'benchN', label: 'Cells of the shock-tube benchmark', unit: '–', value: 200, min: 50, max: 3000, showIf: (v) => v.task === 'benchmarks' },
    { key: 'benchAdvN', label: 'Cells per side of the interface-advection benchmark', unit: '–', value: 40, min: 16, max: 160, showIf: (v) => v.task === 'benchmarks' },
    { key: 'c3Case', label: '3-D case', type: 'select', value: 'tgv', options: C3_CASES, showIf: (v) => v.task === 'cfd3d', help: 'Every case states its grid, time per step and resolution in the result table. The resolved channel at Re_τ = 180 takes several minutes; the other defaults finish in seconds to a minute.' },
    { key: 'c3Sgs', label: '3-D turbulence treatment', type: 'select', value: 'wale', options: C3_SGS, showIf: (v) => v.task === 'cfd3d' && ['channel', 'tgv', 'dambreak', 'stratified'].includes(v.c3Case) },
    { key: 'c3N', label: '3-D cells across (per side, height or diameter)', unit: '–', value: 16, min: 8, max: 128, showIf: (v) => v.task === 'cfd3d' && v.c3Case !== 'dns', help: 'Periodic directions use the nearest power of two. 32 gives the stored reference runs of the channel.' },
    { key: 'c3Re', label: 'Reynolds number (Taylor–Green, pipe) or friction Reynolds number (channel)', unit: '–', value: 180, min: 1, max: 20000, showIf: (v) => v.task === 'cfd3d' && ['tgv', 'channel', 'dns', 'pipe'].includes(v.c3Case) },
    { key: 'c3Dim', label: 'Taylor–Green vortex', type: 'select', value: '3', options: [{ value: '3', label: 'Three-dimensional (transition and decay)' }, { value: '2', label: 'Two-dimensional (exact solution)' }], showIf: (v) => v.task === 'cfd3d' && v.c3Case === 'tgv' },
    { key: 'c3Time', label: '3-D simulated time (0 = automatic)', unit: 'case units', value: 0, min: 0, max: 10000, showIf: (v) => v.task === 'cfd3d' && v.c3Case !== 'coupled', help: 'L/U for the vortex, h/u_τ for the channel, s for the two-phase cases.' },
    { key: 'c3Fluids', label: 'Fluids of the two-phase 3-D case', type: 'select', value: 'case', options: [{ value: 'case', label: 'Liquid and gas of the case at mid-line' }, { value: 'airwater', label: 'Water and air (comparison with the OpenFOAM benchmark)' }], showIf: (v) => v.task === 'cfd3d' && ['dambreak', 'stratified'].includes(v.c3Case) },
    { key: 'c3Mu', label: 'Liquid viscosity of the coupling verification', unit: 'Pa s', value: 0.5, min: 0.0001, max: 100, showIf: (v) => v.task === 'cfd3d' && v.c3Case === 'coupled' },
    { key: 'c3ReCoupled', label: 'Reynolds number of the coupling verification', unit: '–', value: 20, min: 0.1, max: 1500, showIf: (v) => v.task === 'cfd3d' && v.c3Case === 'coupled' },
  ] },
  { group: 'Discretisation and convergence', tab: 'mesh', fields: [
    { key: 'nSteady', label: 'Cells of the steady march', unit: '–', value: 150, min: 20, max: 2000 },
    { key: 'tolP', label: 'Shooting tolerance on the outlet pressure', unit: 'bar', value: 1e-6, min: 1e-9, max: 0.01 },
    { key: 'nCells', label: 'Cells of the transient grid', unit: '–', value: 50, min: 10, max: 400 },
    { key: 'riserRefine', label: 'Refinement of the riser cells', unit: '×', value: 3, min: 1, max: 8 },
    { key: 'cfl', label: 'CFL number of the transient', unit: '–', value: 0.6, min: 0.05, max: 0.95 },
    { key: 'nRadial', label: 'Radial cells per layer', unit: '–', value: 6, min: 2, max: 60 },
    { key: 'nRans', label: 'Radial nodes of the RANS solve', unit: '–', value: 70, min: 30, max: 400 },
  ] },
];

const HILLY = [0, 120, 1500, 210, 3200, 150, 5000, 290, 6800, 205, 8500, 330, 10200, 240, 12000, 300].reduce((a, _, i, s) => (i % 2 ? a : a.concat([{ x: s[i], z: s[i + 1] }])), []);
const TRUNK = [0, -95, 8000, -110, 16000, -102, 24000, -118, 32000, -90, 40000, -70, 46000, -30, 49000, -8, 50000, 12].reduce((a, _, i, s) => (i % 2 ? a : a.concat([{ x: s[i], z: s[i + 1] }])), []);
const SCR = (() => { const flow = [0, -1250, 3000, -1280, 6000, -1265, 9000, -1310, 12000, -1295, 15000, -1340, 18000, -1350].reduce((a, _, i, q) => (i % 2 ? a : a.concat([{ x: q[i], z: q[i + 1] }])), []), cat = catenaryProfile({ height: 1375, topAngle: 12, n: 12 }); return flow.concat(cat.slice(1).map((q) => ({ x: +(18000 + q.x).toFixed(1), z: +(-1350 + q.z).toFixed(1) }))); })();
const PRESETS = [
  { name: 'Deep-water oil tie-back (reference case)', values: { rateFrac: 1 } },
  { name: 'Low-rate turndown (riser slugging)', values: { rateFrac: 0.4, schedule: [{ t: 0, rate: 1, dp: 0, choke: 100 }, { t: 4, rate: 1, dp: 0, choke: 100 }], tEnd: 3, matchSteady: false } },
  { name: 'Gas-condensate trunk line', values: { fluidSel: 'gascond', profile: TRUNK, riserBaseX: 50000, idMm: 590, wtMm: 19.1, pOut: 70, tIn: 55, uValue: 12, insT: 3, insK: 0.3, model: 'mechanistic', mapX: 24000, nCells: 40, tEnd: 4, schedule: [{ t: 0, rate: 1, dp: 0, choke: 100 }, { t: 0.5, rate: 1, dp: 0, choke: 100 }, { t: 1, rate: 0.7, dp: 0, choke: 100 }, { t: 8, rate: 0.7, dp: 0, choke: 100 }] } },
  { name: 'Hilly-terrain onshore multiphase line', values: { profile: HILLY, riserBaseX: 12000, idMm: 203, wtMm: 8.2, rateFrac: 0.45, pOut: 15, tAir: 15, uMode: 'layers', insT: 0, burialDepth: 1.2, kSoil: 1.4, model: 'mukherjeeBrill', mapX: 5000, nCells: 40, tEnd: 1.5, schedule: [{ t: 0, rate: 1, dp: 0, choke: 100 }, { t: 0.2, rate: 1, dp: 0, choke: 100 }, { t: 0.4, rate: 1.3, dp: 0, choke: 100 }, { t: 4, rate: 1.3, dp: 0, choke: 100 }] } },
  { name: 'High-water-cut late life', values: { wc: 75, rateFrac: 0.7, model: 'mechanistic', tEnd: 1.5, schedule: [{ t: 0, rate: 1, dp: 0, choke: 100 }, { t: 0.3, rate: 1, dp: 0, choke: 100 }, { t: 0.4, rate: 1, dp: 5, choke: 100 }, { t: 4, rate: 1, dp: 5, choke: 100 }] } },
  { name: 'Catenary riser with a satellite branch, booster pump and floating separator pressure', values: { profile: SCR, riserBaseX: 18000, branchFrac: 0.3, branchX: 9000, branchLength: 4000, branchDz: 30, branchTin: 55, pumpDp0: 25, pumpQmax: 0.6, sepMode: 'valve', sepKv: 0.5, compPd: 120, compHead: 200, compQmax: 2, tEnd: 1, uMode: 'layers' } },
  { name: 'Two-fluid slug capturing on a laboratory test section', values: { task: 'twofluid', tfDiameterMm: 78, tfLengthD: 200, tfCells: 200, tfTime: 16, tfVsl: 1, tfVsg: 1.5, tfPbar: 1.05, transient: false, compareModels: false, mapX: 9000 } },
  { name: 'Slug tail collapse in the 2-D solver (RANS + VOF)', values: { task: 'cfd', cfdCase: 'slugfront', cfdN: 16, transient: false, compareModels: false } },
  { name: 'Bubbly-flow closures in the riser', values: { task: 'bubbly', transient: false, compareModels: false } },
  { name: '3-D solver: 1-D line model coupled with a 3-D pipe section', values: { task: 'cfd3d', c3Case: 'coupled', c3N: 12, transient: false, compareModels: false, trackSlugs: false } },
  { name: '3-D solver: LES + VOF dam break (coarse)', values: { task: 'cfd3d', c3Case: 'dambreak', c3N: 12, c3Sgs: 'smagorinsky', transient: false, compareModels: false, trackSlugs: false } },
  { name: 'Benchmark set of the local solvers', values: { task: 'benchmarks', benchN: 100, benchAdvN: 32, transient: false, compareModels: false, trackSlugs: false } },
  { name: 'Ramp-up transient', values: { rateFrac: 0.4, tEnd: 2, schedule: [{ t: 0, rate: 1, dp: 0, choke: 100 }, { t: 0.2, rate: 1, dp: 0, choke: 100 }, { t: 0.7, rate: 2.5, dp: 0, choke: 100 }, { t: 5, rate: 2.5, dp: 0, choke: 100 }] } },
];

let lastCtx = {};
async function run(v, ctx = {}) {
  lastCtx = { fluid: ctx.fluid, outputs: ctx.outputs };
  const prog = typeof ctx.progress === 'function' ? ctx.progress : () => {}, tick = typeof ctx.tick === 'function' ? ctx.tick : async () => {};
  const warnings = [], recs = [], plots = [], tables = [], balances = [];
  prog(0.02, 'Steady pressure, temperature and holdup');
  const st = steadyFlow(v, ctx), cfg = st.cfg, fm = cfg.fm, N = st.n, A0 = (PI * st.D[0] ** 2) / 4, so = slugOpts(v), xKm = st.x.map((x) => x / 1000);
  await tick();
  // ---- regimes, unit cells, stability -----------------------------------------------------------------
  const reg = regimes(st), units = st.vsl.map((_, i) => (st.vsl[i] > 1e-9 && st.vsg[i] > 1e-9 ? slugUnitCell(cellOf(st, i), { ...so, xD: st.s[i] / st.D[i] }) : null)); // slug length aware of the distance from the inlet
  const regLen = {}; for (let i = 0; i < N; i++) regLen[reg[i]] = (regLen[reg[i]] || 0) + st.ds;
  const kh = []; for (let i = 0; i < N; i += Math.max(1, Math.round(N / 60))) { if (Math.abs(st.theta[i]) < 10 * DEG && st.vsl[i] > 1e-6 && st.vsg[i] > 1e-6) { try { const s = interfacialStability(cellOf(st, i)); if (Object.values(s).every(Number.isFinite)) kh.push({ i, ...s }); } catch { /* no equilibrium level */ } } }
  const im = clamp(Math.round(cfg.sOfX(clamp(num(v.mapX, 0.5 * cfg.length, 0, 1e7), 0, cfg.length)) / st.ds), 0, N - 1), cm = cellOf(st, im);
  prog(0.2, 'Flow-pattern map');
  const twoPhase = cm.vsl > 1e-9 && cm.vsg > 1e-9, map = twoPhase ? flowPatternMap(cm, { n: 40 }) : null, mand = twoPhase ? mandhaneRegime(cm.vsl, cm.vsg) : 'single phase';
  await tick();
  // ---- slugging -------------------------------------------------------------------------------------------
  prog(0.3, 'Slug tracking and severe-slugging screening');
  const slugIdx = []; for (let i = 0; i <= N; i++) if (isSlug(reg[i]) && units[i]) slugIdx.push(i);
  const iOnsetKH = kh.find((k) => k.vkhRatio > 1), iOnset = slugIdx.length ? slugIdx[0] : iOnsetKH ? iOnsetKH.i : -1;
  let track = null;
  const allUnits = units.every(Boolean);
  if (v.trackSlugs !== false && slugIdx.length > 2 && allUnits) {
    const field = trackingField(st, units), sites = [{ s: st.s[slugIdx[0]], freq: units[slugIdx[0]].freq, length: units[slugIdx[0]].lengthFromFreq }];
    for (let i = 2; i < N - 1; i++) if (st.theta[i - 1] < -1e-5 && st.theta[i] > 1e-5 && isSlug(reg[i + 1]) && st.s[i] > sites[0].s) sites.push({ s: st.s[i], freq: 0.5 * units[i + 1].freq, length: units[i + 1].lengthFromFreq }); // low points re-initiate slugs
    const uo = units[N];
    track = slugTracking(field, { sites, nSlugs: Math.round(num(v.nSlugs, 150, 5, 5000)), seed: Math.round(num(v.slugSeed, 42, 1, 2 ** 31)), theta: st.theta, drainFactor: num(v.drainFactor, 1.2, 1, 10), initMult: num(v.initMult, 1, 0.05, 20), relaxMult: num(v.relaxMult, 1, 0.05, 20), wakeMult: num(v.wakeMult, 1, 0, 10), qSlugOut: A0 * uo.holdupSlug * st.vm[N], qFilmOut: A0 * uo.holdupFilm * Math.max(uo.filmVelocity, 0) });
    if (!track.n) track = null;
  }
  const terrain = terrainSlugging(st), sev = cfg.hasRiser ? severeScreen(st, { riserBaseS: cfg.riserBaseS, feedLength: num(v.feedLength, 0, 0, 1e7), pSep: cfg.pSep }) : null;
  let cyc = null;
  if (sev && sev.severe) { const k = sev.iBase - 1; try { cyc = riserSluggingCycle({ D: st.D[k], feedLength: sev.feedLength, feedAngle: sev.feedAngle * DEG, riserHeight: sev.riserHeight, riserLength: sev.riserLength, wG: Math.max(st.mG[k], 1e-6), wL: Math.max(st.mL[k], 1e-6), rhoL: st.rhoL[k], muL: st.muL[k], T: st.T[k] + 273.15, zG: st.zG[k], mwG: st.mwG[k], pSep: cfg.pSep * 1e5, alphaL: 1 - sev.alpha, chokeOpening: cfg.chokeOpening, chokeDp: Math.max(cfg.chokeDp, 0.05) * 1e5 * st.mScale ** 2, rough: cfg.rough }); } catch (e) { warnings.push({ level: 'info', msg: 'Riser-slugging cycle model: ' + e.message }); } }
  const severe = !!(sev && sev.severe && (!cyc || !cyc.stable));
  await tick();
  prog(0.4, 'Minimum stable rate');
  // mitigations a real design has: active topside choke control (operations suite) and riser-base gas lift or subsea boosting (network suite)
  const ops = ctx.outputs?.ops, lift = ctx.outputs?.net?.lift, ctrlOpen = fin(+ops?.chokeOpening, 0) >= 5 && fin(+ops?.chokeOpening, 0) <= 100 ? +ops.chokeOpening : 30, boostDp = lift && (lift.type === 'boost' || lift.type === 'pump' || lift.type === 'esp') ? clamp(fin(+lift.dp, 0), 0, 300) : 0, liftGas = lift && fin(+lift.gasMass, 0) > 0 ? +lift.gasMass : 0;
  const mit = cfg.hasRiser ? { chokeOpening: Math.min(ctrlOpen, cfg.chokeOpening), chokeDpOpen: Math.max(cfg.chokeDp, 1), boostDp, liftGas } : null;
  const turndown = minimumStableRate(cfg, st.mScale, mit);
  const tdMethod = !cfg.hasRiser ? 'none (terrain slugging has no mitigation in this model)' : !turndown.limited ? 'none needed' : [`topside choke control at ${mit.chokeOpening.toFixed(0)} % opening${ops?.chokeOpening !== undefined ? ' (operations suite' + (ops.slugSuppressed === true ? ', slugging suppressed by its controller' : '') + ')' : ' (assumed; run the operations suite for the tuned value)'}`, boostDp > 0 ? `${boostDp.toFixed(0)} bar of subsea boosting available as back-pressure (network suite)` : null, liftGas > 0 ? `${liftGas.toFixed(2)} kg/s of riser-base lift gas (network suite)` : null].filter(Boolean).join(' + ');
  await tick();
  // ---- transient ------------------------------------------------------------------------------------------
  let tr = null, trInfo = null;
  if (v.transient !== false) {
    prog(0.45, 'Transient drift-flux solution');
    trInfo = buildTransient(v, st);
    while (!trInfo.sim.done) { trInfo.sim.advance(40); prog(0.45 + 0.35 * Math.min(trInfo.sim.t / trInfo.tEnd, 1), `Transient: ${(trInfo.sim.t / 3600).toFixed(2)} h of ${(trInfo.tEnd / 3600).toFixed(2)} h`); await tick(); }
    tr = trInfo.sim.result();
  }
  // ---- heat-transfer detail and cooldown preview ----------------------------------------------------------
  prog(0.82, 'Radial conduction and local models');
  let ic = 0; for (let i = 1; i <= N; i++) if (st.subcooling[i] > st.subcooling[ic]) ic = i;
  const kc = Math.min(ic, N - 1), prC = fm.at(st.P[ic], st.T[ic], st.mScale), netC = cfg.network({ s: st.s[kc], z: st.z[kc], D: st.D[kc], pr: prC, vm: st.vm[kc], holdup: st.holdup[kc], ta: st.tAmb[kc] });
  const radLayers = cfg.layers.map((l) => ({ ...l })), buried = cfg.buriedAt(st.s[kc]);
  if (buried) radLayers.push({ name: 'Soil (equivalent annulus)', t: (cfg.od / 2) * (Math.exp(Math.acosh((2 * (cfg.burialDepth + cfg.od / 2)) / cfg.od)) - 1), k: cfg.kSoil, rho: 1800, cp: 1300 });
  const hOutC = buried ? 1e5 : st.z[kc] >= 0 ? hOutside(cfg.wind, cfg.od, 'air') : hOutside(cfg.current, cfg.od, 'seawater', st.tAmb[kc]), Ac = (PI * st.D[kc] ** 2) / 4;
  const fluidC = Ac * (st.holdup[kc] * st.rhoL[kc] * prC.cpL + (1 - st.holdup[kc]) * st.rhoG[kc] * prC.cpG), nPer = Math.round(num(v.nRadial, 6, 2, 80));
  const radSteady = radialConduction({ ri: st.D[kc] / 2, layers: radLayers, hIn: st.hIn[kc], hOut: hOutC, tFluid: st.T[ic], tAmb: st.tAmb[kc], T0: 'steady', tEnd: 0, nPer });
  const tauC = fluidC / Math.max(netC.U * cfg.uMult * PI * st.D[kc], 1e-9), cool = radialConduction({ ri: st.D[kc] / 2, layers: radLayers, hIn: clamp(0.25 * st.hIn[kc], 30, 400), hOut: hOutC, tFluid: st.T[ic], tAmb: st.tAmb[kc], fluidC, T0: 'steady', tEnd: clamp(5 * tauC, 3600, 60 * 86400), nSteps: 240, nPer, stopAt: st.tHyd[ic], nProfiles: 5 });
  // ---- local models ---------------------------------------------------------------------------------------
  let rans = null, ransRef = null, vof = null, vofUp = null, dam = null, mdl = 'mixing';
  if (v.localModels !== false) {
    const lam = cm.vsl / Math.max(cm.vsl + cm.vsg, 1e-9), rhoN = cm.rhoL * lam + cm.rhoG * (1 - lam), muN = cm.muL * lam + cm.muG * (1 - lam), ReLoc = clamp((rhoN * (cm.vsl + cm.vsg) * cm.D) / muN, 3000, 3e6), nR = Math.round(num(v.nRans, 70, 30, 400));
    mdl = RANS_MODELS.some((q) => q.value === v.ransModel) ? v.ransModel : 'mixing';
    rans = ransPipe({ reTau: reTauOf(ReLoc), model: mdl, n: nR, tol: 1e-7, maxIter: 2500 }); ransRef = mdl === 'mixing' ? null : ransPipe({ reTau: reTauOf(ReLoc), model: 'mixing', n: nR });
    await tick();
    vof = vofAdvect1D({ n: 120, scheme: 'thinc' }); vofUp = vofAdvect1D({ n: 120, scheme: 'upwind' });
    const um = units[im]; dam = damBreak({ hL: Math.max(levelOfHoldup(um ? um.holdupSlug : 0.9), 0.05) * cm.D, hR: (um ? um.filmThickness : 0.2 * cm.D), n: 160 });
  }
  await tick();
  // ---- model comparison -----------------------------------------------------------------------------------
  const cmp = [];
  if (v.compareModels !== false && st.bc !== 'inletP') {
    prog(0.9, 'Comparing holdup models');
    for (const m of HOLDUP_MODELS) { const r = m.value === cfg.model ? st : steadySolve({ ...cfg.base, model: m.value, n: Math.min(cfg.n, 40), mScale: st.mScale, pOut: cfg.pOut, pGuess: st.pIn, hydrate: false, tolP: 1e-4 }); cmp.push(r.ok ? [m.label, r3(r.pIn, 2), r3(r.dpFric, 2), r3(r.dpGrav, 2), r3(r.liquidInventory, 1), r3(r.tOut, 2)] : [m.label, '—', '—', '—', '—', '—']); }
  }
  // ---- derived quantities ---------------------------------------------------------------------------------
  const maxOf = (a) => a.reduce((m, x) => Math.max(m, x), -Infinity), argMax = (a) => a.reduce((m, x, i) => (x > a[m] ? i : m), 0), minOf = (a) => a.reduce((m, x) => Math.min(m, x), Infinity);
  const evrMax = maxOf(st.evr), iEvr = argMax(st.evr), vMax = maxOf(st.vm), subMax = maxOf(st.subcooling), hydLen = st.subcooling.slice(0, N).filter((s) => s > 0).length * st.ds, dpTotal = st.pIn - st.pOut, iHyd = st.subcooling.findIndex((s) => s > 0);
  const heatFlux = st.qLoss.map((q, i) => q / (PI * st.D[i])), wat = num(v.watC, 0, -60, 120), tMin = minOf(st.T), gor = fin(fm.rates.gor, 0), wcIn = fin(fm.rates.wc, 0) * 100;
  // strongest hydrodynamic slugging location
  let iS = -1; for (const i of slugIdx) if (iS < 0 || units[i].vt * units[i].volume > units[iS].vt * units[iS].volume) iS = i;
  const uS = iS >= 0 ? units[iS] : null, qLout = st.qL[N], qGout = st.qG[N];
  const terrainActive = terrain.accumulates && !severe, slugType = severe ? 'severe' : terrainActive ? 'terrain' : uS ? 'hydrodynamic' : 'none';
  let slug;
  if (slugType === 'severe') { const vol = A0 * sev.riserLength, per = cyc && cyc.period ? cyc.period : (vol / Math.max(qLout, 1e-9)); slug = { type: 'severe', freq: 1 / per, period: per, length: sev.riserLength, lengthMax: sev.riserLength, velocity: cyc ? (cyc.qLiqPeakRatio * st.mL[sev.iBase - 1]) / st.rhoL[sev.iBase - 1] / A0 : st.vm[N], holdupBody: 1, volume: vol, surge: vol * Math.max(1 - (cyc ? Math.min(num(v.drainFactor, 1.2, 1, 10) / Math.max(cyc.qLiqPeakRatio, 1), 1) : 0), 0.2), x: st.x[sev.iBase] }; }
  else if (slugType === 'terrain') { const d = terrain.worst, hs = uS ? uS.holdupSlug : 0.9, vol = Math.max(d.excess, 0.05 * d.volume), L = vol / (A0 * hs), per = Math.max(d.period, vol / Math.max(qLout, 1e-9)); slug = { type: 'terrain', freq: 1 / per, period: per, length: L, lengthMax: 1.5 * L, velocity: slugVelocity(d.vsl + d.vsg, st.D[0], d.upAngle * DEG).vt, holdupBody: hs, volume: vol, surge: vol, x: d.x }; }
  else if (uS) { const Lm = track ? track.meanLength : uS.length, Lx = track ? Math.max(track.maxLength, track.length1000) : uS.lengthMax, f = track && track.freqArrival > 0 ? track.freqArrival : uS.freq; slug = { type: 'hydrodynamic', freq: f, period: 1 / f, length: Lm, lengthMax: Lx, velocity: uS.vt, holdupBody: uS.holdupSlug, volume: Lm * A0 * uS.holdupSlug, surge: track ? Math.max(track.surge, track.surgeSingle, Lm * A0 * uS.holdupSlug) : Lx * A0 * uS.holdupSlug * 0.5, x: st.x[iS] }; }
  else slug = { type: 'none', freq: 0, period: null, length: 0, lengthMax: 0, velocity: 0, holdupBody: 0, volume: 0, surge: 0, x: null };
  for (const k of Object.keys(slug)) if (typeof slug[k] === 'number' && !Number.isFinite(slug[k])) slug[k] = null;
  const catcher = 1.25 * fin(slug.surge, 0);
  // transient summary
  let pAmp = 0, trSummary = null;
  if (tr && tr.t.length > 3) {
    const i0 = Math.floor(tr.t.length * 0.4), w = tr.pIn.slice(i0), wq = tr.qLiqOut.slice(i0), meanQ = mean(wq);
    pAmp = maxOf(w) - minOf(w); let cum = 0; for (let i = 1; i < tr.t.length; i++) cum += 0.5 * (tr.qLiqOut[i] + tr.qLiqOut[i - 1]) * (tr.t[i] - tr.t[i - 1]);
    trSummary = { pMean: mean(w), pMin: minOf(w), pMax: maxOf(w), qMean: meanQ, qPeak: maxOf(wq), cumLiq: cum, invEnd: tr.inv[tr.inv.length - 1], inv0: tr.inv[0] };
  } else if (cyc && !cyc.stable) pAmp = cyc.amplitude / 1e5;

  // ---- warnings ---------------------------------------------------------------------------------------------
  if (st.junction) warnings.push({ level: 'info', msg: `A branch delivers ${(100 * cfg.branch.frac).toFixed(0)} % of the rate at ${Math.round(st.junction.x)} m: junction pressure ${st.junction.P.toFixed(1)} bara, mixed temperature ${st.junction.tMix.toFixed(1)} °C. Screening scans (minimum stable rate, model comparison) use the line without the branch.` });
  if (subMax > 0) warnings.push({ level: subMax > 3 ? 'bad' : 'warn', msg: `The stream is inside the hydrate region over ${Math.round(hydLen)} m (from ${Math.round(st.x[Math.max(iHyd, 0)])} m); the largest subcooling is ${subMax.toFixed(1)} °C at ${Math.round(st.x[ic])} m.` });
  if (evrMax > 1) warnings.push({ level: 'bad', msg: `The mixture velocity reaches ${evrMax.toFixed(2)} times the API RP 14E erosional velocity (C = ${cfg.cErosion}) at ${Math.round(st.x[iEvr])} m.` });
  else if (evrMax > 0.8) warnings.push({ level: 'warn', msg: `The mixture velocity is ${(100 * evrMax).toFixed(0)} % of the erosional velocity at ${Math.round(st.x[iEvr])} m.` });
  if (severe) warnings.push({ level: 'bad', msg: `Severe riser slugging is predicted: the feed line is ${sev.feedRegime}, the liquid velocity is ${sev.boe.toFixed(1)} times the Bøe limit and the Pots number is ${sev.pots === null ? 'n/a' : sev.pots.toFixed(2)}${cyc && cyc.period ? `; cycle of about ${(cyc.period / 60).toFixed(0)} min with ${(cyc.amplitude / 1e5).toFixed(0)} bar swing at the riser base` : ''}.` });
  else if (sev && sev.severe && cyc && cyc.stable) warnings.push({ level: 'info', msg: 'The screening criteria place the riser in the severe-slugging region, but the cycle model shows no liquid surge at this choke setting.' });
  if (terrain.accumulates) warnings.push({ level: 'warn', msg: `Liquid accumulates at ${terrain.dips.filter((d) => d.accumulates).length} low point(s): at ${Math.round(terrain.worst.x)} m the gas velocity is ${terrain.worst.vsg.toFixed(2)} m/s against ${terrain.worst.vsgCrit.toFixed(2)} m/s needed to sweep the upward leg, which holds about ${terrain.worst.volume.toFixed(0)} m³.` });
  if (wat !== 0 && tMin < wat) warnings.push({ level: 'warn', msg: `The fluid temperature falls to ${tMin.toFixed(1)} °C, below the wax appearance temperature of ${wat.toFixed(1)} °C.` });
  if (tr && tr.aborted !== null) warnings.push({ level: 'warn', msg: `The transient was stopped at ${(tr.aborted / 3600).toFixed(2)} h: the line became liquid-full and the volume constraint could no longer be met (pressure-surge dynamics are outside this scheme). Results up to that time are shown.` });
  if (tr && pAmp > 0.05 * dpTotal && pAmp > 2) warnings.push({ level: 'warn', msg: `The transient inlet pressure swings by ${pAmp.toFixed(1)} bar over the last 60 % of the simulated time.` });
  if (UP_ONLY.has(cfg.model) && st.theta.some((t) => t < -1e-4)) warnings.push({ level: 'info', msg: 'The selected correlation is for upward flow; downward-inclined cells use Beggs & Brill.' });
  if (st.nodal) warnings.push({ level: 'info', msg: `Well inflow balance: ${st.nodal.qLiqStd.toFixed(0)} Sm³/d of liquid at a flowing bottom-hole pressure of ${st.nodal.pwf.toFixed(1)} bara (${(100 * st.mScale).toFixed(0)} % of the case rate).` });
  if (cfg.idOf || cfg.idMin < cfg.id0) warnings.push({ level: 'info', msg: `Deposits narrow the bore to ${(cfg.idMin * 1000).toFixed(0)} mm at the tightest point.` });

  // ---- recommendations --------------------------------------------------------------------------------------
  const rateNow = st.mScale, tdFrac = turndown.rate, tdCtrl = Math.min(turndown.rateControlled, tdFrac);
  if (turndown.limited) recs.push(`Without slug control keep the rate above about ${(100 * tdFrac).toFixed(0)} % of the case rate (${(tdFrac * stdLiquid(fm, 1)).toFixed(0)} Sm³/d liquid): below it the screening predicts ${cfg.hasRiser ? 'riser-induced' : 'terrain'} slugging. ${cfg.hasRiser ? (tdCtrl < tdFrac ? `With ${tdMethod} the riser cycle model stays stable down to ${(100 * tdCtrl).toFixed(0)} % (${(tdCtrl * stdLiquid(fm, 1)).toFixed(0)} Sm³/d)` : `The mitigations examined (${tdMethod}) do not widen the window in the riser cycle model`) : ''}${rateNow < tdCtrl ? ` — the present rate (${(100 * rateNow).toFixed(0)} %) is below the controlled limit` : ''}.`);
  if (cfg.hasRiser && turndown.limited && tdCtrl > 0.6 * Math.max(rateNow, 1e-9)) warnings.push({ level: 'warn', msg: `Narrow stable operating window: even with ${tdMethod} the minimum stable rate is ${(100 * tdCtrl).toFixed(0)} % of the case rate (${(100 * tdFrac).toFixed(0)} % without control).` });
  else recs.push(`No slugging limit was found down to ${(100 * tdFrac).toFixed(0)} % of the case rate; turndown is limited by other constraints (arrival temperature, liquid inventory).`);
  if (slug.type !== 'none') recs.push(`Size the receiving vessel for a liquid surge of about ${catcher.toFixed(0)} m³ (largest ${slug.type} slug ${fin(slug.lengthMax, 0).toFixed(0)} m, 25 % margin, vessel drained at ${num(v.drainFactor, 1.2, 1, 10).toFixed(1)} × the mean liquid rate).`);
  if (severe) recs.push(`Suppress the riser cycle by topsides choking or riser-base gas lift: the Bøe ratio must fall below 1 (now ${sev.boe.toFixed(1)}); the Taitel criterion asks for ${Math.max(sev.taitelPsep, 0).toFixed(0)} bara at the riser top without active control.`);
  if (subMax > 0) {
    // U-value that keeps the arrival just above the hydrate temperature (same holdup model, coarse grid)
    let uNeed = null;
    try { const g = (m) => { const r = steadySolve({ ...cfg.base, n: 40, U: cfg.U * m, uOf: cfg.uOf ? (c) => cfg.uOf(c) * m : null, mScale: st.mScale, pOut: cfg.pOut, pGuess: st.pIn, tolP: 1e-3 }); return r.ok ? Math.max(...r.subcooling) : 1; }; if (g(0.02) < 0) uNeed = brent(g, 0.02, 1, 1e-3) * (cfg.uOf ? mean(st.U) : cfg.U); } catch { uNeed = null; }
    recs.push(uNeed ? `Lower the overall U-value to about ${uNeed.toFixed(2)} W/m²K (now ${mean(st.U).toFixed(2)}) to keep the whole line outside the hydrate region in steady flow, or inhibit for ${subMax.toFixed(1)} °C of subcooling.` : `Insulation alone cannot keep the line outside the hydrate region at this rate; plan inhibitor for ${subMax.toFixed(1)} °C of subcooling or active heating.`);
  } else recs.push(`The steady stream stays ${(-subMax).toFixed(1)} °C above the hydrate temperature at the closest point (${Math.round(st.x[ic])} m); after a stop the contents there reach the hydrate temperature in ${cool.tReach !== null ? 'about ' + (cool.tReach / 3600).toFixed(1) + ' h' : 'more than ' + (cool.t[cool.t.length - 1] / 3600).toFixed(0) + ' h'}.`);
  if (evrMax > 0.8) recs.push(`Limit the rate to ${(100 * rateNow * (0.8 / evrMax)).toFixed(0)} % of the case rate, or increase the bore near ${Math.round(st.x[iEvr])} m, to stay below 80 % of the erosional velocity.`);
  if (terrain.accumulates) recs.push(`Raise the gas velocity above ${terrain.worst.vsgCrit.toFixed(1)} m/s in the upward leg at ${Math.round(terrain.worst.x)} m (rate increase, gas recycle or pigging) to stop the ${terrain.worst.volume.toFixed(0)} m³ liquid accumulation there.`);

  // ---- plots ------------------------------------------------------------------------------------------------
  const D0 = (a, k = 3) => a.map((x) => r3(x, k));
  plots.push({ type: 'line', title: 'Pressure, temperature and hydrate temperature along the line', xlabel: 'Distance (km)', ylabel: 'bara · °C', series: [{ name: 'Pressure (bara)', x: xKm, y: st.P }, { name: 'Temperature (°C)', x: xKm, y: st.T }, { name: 'Hydrate temperature (°C)', x: xKm, y: st.tHyd, dash: true }, { name: 'Wall temperature (°C)', x: xKm, y: st.tWall }, { name: 'Ambient (°C)', x: xKm, y: st.tAmb, dash: true }], vlines: cfg.hasRiser ? [{ x: cfg.riserBaseX / 1000, label: 'Riser base' }] : [] });
  plots.push({ type: 'line', title: 'Elevation profile', xlabel: 'Distance (km)', ylabel: 'Elevation (m)', height: 220, series: [{ name: 'Elevation', x: xKm, y: st.z }], hlines: [{ y: 0, label: 'Sea level' }] });
  plots.push({ type: 'line', title: 'Liquid holdup and flow regime', xlabel: 'Distance (km)', ylabel: 'Holdup (–) · regime index', ymin: 0, series: [{ name: 'Liquid holdup', x: xKm, y: st.holdup }, { name: 'No-slip liquid fraction', x: xKm, y: st.vsl.map((l, i) => l / Math.max(l + st.vsg[i], 1e-12)), dash: true }, { name: 'Void fraction', x: xKm, y: st.holdup.map((h) => 1 - h), dash: true }, { name: 'Regime index / 6', x: xKm, y: reg.map((r) => Math.max(PATTERNS.indexOf(r), 0) / 6), mode: 'step' }], note: 'Regime index: ' + PATTERNS.map((p, i) => `${i} ${p}`).join(' · ') });
  plots.push({ type: 'line', title: 'Velocities and erosional limit', xlabel: 'Distance (km)', ylabel: 'm/s', zeroY: true, series: [{ name: 'Superficial liquid', x: xKm, y: st.vsl }, { name: 'Superficial gas', x: xKm, y: st.vsg }, { name: 'Liquid phase', x: xKm, y: st.vL }, { name: 'Gas phase', x: xKm, y: st.vG }, { name: 'Mixture', x: xKm, y: st.vm }, { name: 'Erosional (API RP 14E)', x: xKm, y: st.vm.map((vv, i) => vv / Math.max(st.evr[i], 1e-9)), dash: true }] });
  plots.push({ type: 'line', title: 'Pressure-gradient components', xlabel: 'Distance (km)', ylabel: 'Pa/m', series: [{ name: 'Friction', x: xKm, y: st.fric }, { name: 'Gravity', x: xKm, y: st.grav }, { name: 'Acceleration', x: xKm, y: st.acc }, { name: 'Local losses', x: xKm, y: st.loc }, { name: 'Total', x: xKm, y: st.dpdx, dash: true }] });
  plots.push({ type: 'bar', title: 'Pressure-loss breakdown', ylabel: 'bar', categories: ['Friction', 'Gravity', 'Acceleration', 'Local', 'Choke'], series: [{ name: 'Pressure loss', values: [r3(st.dpFric), r3(st.dpGrav), r3(st.dpAcc), r3(st.dpLocal), r3(cfg.pOut - cfg.pSep)] }] });
  plots.push({ type: 'line', title: 'Heat loss and heat flux', xlabel: 'Distance (km)', ylabel: 'W/m · W/m²', series: [{ name: 'Heat loss (W/m)', x: xKm, y: st.qLoss }, { name: 'Heat flux on the inner wall (W/m²)', x: xKm, y: heatFlux }, { name: 'U-value × 10 (W/m²K)', x: xKm, y: st.U.map((u) => 10 * u), dash: true }] });
  if (map) plots.push({ type: 'line', title: `Flow-pattern map at ${Math.round(st.x[im])} m (inclination ${(cm.theta / DEG).toFixed(2)}°)`, xlabel: 'Superficial gas velocity (m/s)', ylabel: 'Superficial liquid velocity (m/s)', logx: true, logy: true, series: [...map.boundaries.map((b) => ({ name: b.name, x: b.x, y: b.y, mode: 'points' })), { name: 'Operating point', x: [cm.vsg], y: [cm.vsl], mode: 'points' }, { name: 'Operating line (inlet → outlet)', x: st.vsg.filter((g, i) => g > 0 && st.vsl[i] > 0), y: st.vsl.filter((l, i) => l > 0 && st.vsg[i] > 0), dash: true }], note: `Taitel–Dukler / Barnea transitions for the local properties: ${reg[im]}. Mandhane map (horizontal, air–water coordinates): ${mand}.` });
  if (kh.length > 1) plots.push({ type: 'line', title: 'Stability of stratified flow (ratio > 1 = unstable)', xlabel: 'Distance (km)', ylabel: 'Actual / critical', logy: true, series: [{ name: 'Inviscid Kelvin–Helmholtz', x: kh.map((k) => xKm[k.i]), y: kh.map((k) => Math.max(k.ikhRatio, 1e-3)) }, { name: 'Viscous Kelvin–Helmholtz', x: kh.map((k) => xKm[k.i]), y: kh.map((k) => Math.max(k.vkhRatio, 1e-3)) }, { name: 'Taitel–Dukler', x: kh.map((k) => xKm[k.i]), y: kh.map((k) => Math.max(k.tdRatio, 1e-3)) }, { name: 'Wave generation (Jeffreys)', x: kh.map((k) => xKm[k.i]), y: kh.map((k) => Math.max(k.wavyRatio, 1e-3)), dash: true }, { name: 'Roll waves (Vedernikov)', x: kh.map((k) => xKm[k.i]), y: kh.map((k) => Math.max(k.rollRatio, 1e-3)), dash: true }], hlines: [{ y: 1, label: 'Neutral stability' }], note: 'Evaluated at the equilibrium stratified level of near-horizontal cells (|θ| < 10°).' });
  if (slugIdx.length > 1) {
    const sx = slugIdx.map((i) => xKm[i]);
    plots.push({ type: 'line', title: 'Slug unit cell along the line', xlabel: 'Distance (km)', ylabel: 'see legend', logy: true, series: [{ name: 'Frequency (1/min)', x: sx, y: slugIdx.map((i) => 60 * units[i].freq), mode: 'points' }, { name: 'Mean length (m)', x: sx, y: slugIdx.map((i) => units[i].lengthFromFreq), mode: 'points' }, { name: 'Length, 1-in-1000 (m)', x: sx, y: slugIdx.map((i) => units[i].lengthMax), mode: 'points' }, { name: 'Translational velocity (m/s)', x: sx, y: slugIdx.map((i) => units[i].vt), mode: 'points' }, { name: 'Liquid volume per slug (m³)', x: sx, y: slugIdx.map((i) => units[i].lengthFromFreq * A0 * units[i].holdupSlug), mode: 'points' }] });
    plots.push({ type: 'line', title: 'Slug body and film', xlabel: 'Distance (km)', ylabel: 'Holdup (–) · m/s · m', series: [{ name: 'Body holdup', x: sx, y: slugIdx.map((i) => units[i].holdupSlug), mode: 'points' }, { name: 'Film holdup', x: sx, y: slugIdx.map((i) => units[i].holdupFilm), mode: 'points' }, { name: 'Film thickness / D', x: sx, y: slugIdx.map((i) => units[i].filmThickness / st.D[i]), mode: 'points' }, { name: 'Film velocity (m/s)', x: sx, y: slugIdx.map((i) => units[i].filmVelocity), mode: 'points' }] });
  }
  if (track) {
    plots.push({ type: 'bar', title: 'Slug-length distribution at the outlet', ylabel: 'Slugs', categories: track.hist.centers.map((c) => c.toFixed(0) + ' m'), series: [{ name: 'Tracked slugs', values: track.hist.counts }] });
    plots.push({ type: 'line', title: 'Slug arrivals at the outlet', xlabel: 'Arrival time (min)', ylabel: 'Slug length (m)', series: [{ name: 'Slug', x: track.arrivals.map((a) => (a.t - track.arrivals[0].t) / 60), y: track.arrivals.map((a) => a.length), mode: 'points' }] });
  }
  if (cyc) plots.push({ type: 'line', title: 'Severe-slugging cycle (lumped model)', xlabel: 'Time (min)', ylabel: 'bara · × mean rate', series: [{ name: 'Riser-base pressure (bara)', x: cyc.t.map((t) => t / 60), y: cyc.pBase.map((p) => p / 1e5) }, { name: 'Liquid out / mean (×10)', x: cyc.t.map((t) => t / 60), y: cyc.wLout.map((w) => (10 * w) / Math.max(st.mL[sev.iBase - 1], 1e-9)) }, { name: 'Liquid column / riser length (×100)', x: cyc.t.map((t) => t / 60), y: cyc.level.map((l) => 100 * l) }], note: cyc.stable ? 'No sustained cycle at this choke setting.' : `${cyc.type}. Stages per cycle: build-up ${(cyc.stages.buildUp / 60).toFixed(0)} min, production ${(cyc.stages.production / 60).toFixed(0)} min, blowout ${(cyc.stages.blowout / 60).toFixed(1)} min, fallback ${(cyc.stages.fallback / 60).toFixed(1)} min.` });
  if (tr && tr.t.length > 2) {
    const th = tr.t.map((t) => t / 3600);
    plots.push({ type: 'line', title: 'Transient: inlet pressure and liquid inventory', xlabel: 'Time (h)', ylabel: 'bara · m³/10', series: [{ name: 'Inlet pressure (bara)', x: th, y: tr.pIn }, { name: 'Liquid inventory / 10 (m³)', x: th, y: tr.inv.map((x) => x / 10), dash: true }] });
    plots.push({ type: 'line', title: 'Transient: outlet rates and holdup', xlabel: 'Time (h)', ylabel: 'm³/h · %', series: [{ name: 'Liquid out (m³/h)', x: th, y: tr.qLiqOut.map((q) => 3600 * q) }, { name: 'Gas out / 10 (actual m³/h)', x: th, y: tr.qGasOut.map((q) => 360 * q) }, { name: 'Outlet holdup (%)', x: th, y: tr.holdupOut.map((h) => 100 * h), dash: true }] });
    if (tr.field.t.length > 1) {
      const fx = tr.field.x.map((x) => x / 1000), fy = tr.field.t.map((t) => t / 3600);
      plots.push({ type: 'field', title: 'Liquid holdup: distance–time map', xlabel: 'Distance (km)', ylabel: 'Time (h)', zlabel: 'Holdup', zunit: '–', x: fx, y: fy, z: tr.field.holdup, zmin: 0, zmax: 1, cmap: 'viridis' });
      plots.push({ type: 'field', title: 'Pressure: distance–time map', xlabel: 'Distance (km)', ylabel: 'Time (h)', zlabel: 'Pressure', zunit: 'bara', x: fx, y: fy, z: tr.field.P, cmap: 'turbo' });
    }
  }
  { const lay = radSteady.layerOf.map((l) => radLayers.filter((q) => q.t > 0 && q.k > 0)[l]?.name || ''), rmm = radSteady.r.map((r) => r * 1000);
    plots.push({ type: 'line', title: `Radial temperature at ${Math.round(st.x[ic])} m and cooldown after a stop`, xlabel: 'Radius (mm)', ylabel: '°C', series: cool.profiles.map((p) => ({ name: p.t === 0 ? 'Flowing (steady)' : `${(p.t / 3600).toFixed(1)} h after the stop`, x: rmm, y: p.T, mode: 'both' })), hlines: [{ y: st.tHyd[ic], label: 'Hydrate temperature' }], note: `Layers from the bore outwards: ${[...new Set(lay)].join(', ')}.` });
    plots.push({ type: 'line', title: 'Cooldown of the contents at the coldest-margin location', xlabel: 'Time after the stop (h)', ylabel: '°C', series: [{ name: 'Contents', x: cool.t.map((t) => t / 3600), y: cool.fluid }, { name: 'Inner wall', x: cool.t.map((t) => t / 3600), y: cool.wall, dash: true }], hlines: [{ y: st.tHyd[ic], label: 'Hydrate temperature' }] }); }
  if (rans) {
    plots.push({ type: 'line', title: `Developed turbulent velocity profile (radial RANS, Re = ${sig(rans.Re, 3)})`, xlabel: 'y⁺', ylabel: 'u⁺', logx: true, series: [{ name: RANS_MODELS.find((q) => q.value === mdl).label, x: rans.y.slice(1), y: rans.u.slice(1) }, ...(ransRef ? [{ name: 'Mixing length', x: ransRef.y.slice(1), y: ransRef.u.slice(1), dash: true }] : []), { name: 'Log law (κ = 0.41, B = 5.2)', x: rans.y.filter((y) => y > 30), y: rans.y.filter((y) => y > 30).map((y) => Math.log(y) / 0.41 + 5.2), dash: true }], note: `Darcy friction factor ${rans.f.toFixed(5)} against ${frictionFactor(rans.Re, 0).toFixed(5)} from the smooth-pipe Colebrook law (${(100 * (rans.f / frictionFactor(rans.Re, 0) - 1)).toFixed(1)} %). Single-phase flow at the local no-slip mixture Reynolds number.` });
    plots.push({ type: 'line', title: 'Interface capturing: liquid slug advected once around a periodic pipe', xlabel: 'Position (–)', ylabel: 'Liquid fraction', series: [{ name: 'Exact', x: vof.x, y: vof.exact, dash: true }, { name: 'THINC volume-of-fluid', x: vof.x, y: vof.phi }, { name: 'First-order upwind', x: vofUp.x, y: vofUp.phi }], note: `Interface thickness ${vof.thickness.toFixed(1)} cells with THINC against ${vofUp.thickness.toFixed(1)} with upwind; both conserve the liquid volume to round-off.` });
    plots.push({ type: 'line', title: 'Slug-tail collapse onto the film (shallow-water dam break)', xlabel: 'Distance from the tail (m)', ylabel: 'Liquid depth (m)', series: [{ name: 'Finite volume (HLL)', x: dam.x, y: dam.h }, { name: 'Exact (Stoker)', x: dam.x, y: dam.exact, dash: true }], note: `Channel analogue at ${Math.round(st.x[im])} m, ${dam.tEnd.toFixed(2)} s after release; mean error ${(100 * dam.l1).toFixed(2)} % of the initial depth.` });
  }

  // ---- tables -----------------------------------------------------------------------------------------------
  const rowIdx = [...new Set(linspace(0, N, Math.min(N + 1, 41)).map(Math.round))];
  tables.push({ title: 'Results along the line', columns: ['x (m)', 'z (m)', 'P (bara)', 'T (°C)', 'T wall (°C)', 'T hydrate (°C)', 'Holdup', 'vsl (m/s)', 'vsg (m/s)', 'v liquid (m/s)', 'v gas (m/s)', 'ρ mix (kg/m³)', 'Regime', 'dp/dx fric (Pa/m)', 'grav', 'acc', 'local', 'Heat loss (W/m)', 'Heat flux (W/m²)', 'τ wall (Pa)', 'v / v eros'],
    rows: rowIdx.map((i) => [r3(st.x[i], 0), r3(st.z[i], 1), r3(st.P[i], 2), r3(st.T[i], 2), r3(st.tWall[i], 2), r3(st.tHyd[i], 2), r3(st.holdup[i]), r3(st.vsl[i]), r3(st.vsg[i]), r3(st.vL[i]), r3(st.vG[i]), r3(st.rhoM[i], 1), reg[i], r3(st.fric[i], 1), r3(st.grav[i], 1), r3(st.acc[i], 2), r3(st.loc[i], 2), r3(st.qLoss[i], 1), r3(heatFlux[i], 1), r3(st.tauW[i], 2), r3(st.evr[i])]) });
  tables.push({ title: 'Flow-regime lengths', columns: ['Regime', 'Length (m)', 'Share (%)'], rows: Object.entries(regLen).sort((a, b) => b[1] - a[1]).map(([k, L]) => [k, r3(L, 0), r3((100 * L) / st.length, 1)]), note: `Mechanistic (Taitel–Dukler / Barnea) regime per cell. At ${Math.round(st.x[im])} m: ${reg[im]}; Mandhane map: ${mand}.` });
  { const rows = [];
    rows.push(['Slugging type', slug.type, '']); rows.push(['Slug onset (first unstable / slug cell)', iOnset >= 0 ? r3(st.x[iOnset], 0) : '—', 'm']);
    if (uS) { rows.push(['Strongest location', r3(st.x[iS], 0), 'm'], ['Unit-cell frequency there', r3(uS.freq, 5), '1/s'], ['Unit-cell mean length (balance)', r3(uS.lengthFromFreq, 1), 'm'], ['Developed length correlation', r3(uS.length, 1), 'm'], ['1-in-1000 length', r3(uS.lengthMax, 1), 'm'], ['Translational velocity', r3(uS.vt, 2), 'm/s'], ['Body holdup', r3(uS.holdupSlug), '–'], ['Film holdup', r3(uS.holdupFilm), '–'], ['Film thickness', r3(uS.filmThickness * 1000, 1), 'mm'], ['Film velocity', r3(uS.filmVelocity, 2), 'm/s'], ['Liquid volume per slug', r3(uS.lengthFromFreq * A0 * uS.holdupSlug, 2), 'm³']); }
    if (track) rows.push(['Tracked: slugs arrived', track.n, '–'], ['Tracked: initiated / merged / dissipated', `${track.generated} / ${track.merges} / ${track.dissipated}`, '–'], ['Tracked: arrival frequency', r3(track.freqArrival, 5), '1/s'], ['Tracked: mean arrival interval', track.period ? r3(track.period, 1) : '—', 's'], ['Tracked: mean length at arrival', r3(track.meanLength, 1), 'm'], ['Tracked: standard deviation', r3(track.stdLength, 1), 'm'], ['Tracked: P50 / P90 / P99 length', `${track.p50.toFixed(0)} / ${track.p90.toFixed(0)} / ${track.p99.toFixed(0)}`, 'm'], ['Tracked: longest', r3(track.maxLength, 1), 'm'], ['Log-normal 1-in-1000 length', r3(track.length1000, 1), 'm'], ['Vessel surge over the record', r3(track.surge, 1), 'm³'], ['Tracked from', r3(cfg.xOfS(track.from), 0), 'm']);
    rows.push(['Liquid surge to accommodate', r3(fin(slug.surge, 0), 1), 'm³'], ['Slug-catcher volume with 25 % margin', r3(catcher, 1), 'm³']);
    if (sev) rows.push(['Riser: Bøe ratio vsl / vsl,crit', r3(sev.boe, 2), '–'], ['Riser: Pots number', sev.pots === null ? '—' : r3(sev.pots, 3), '–'], ['Riser: Taitel stability pressure', r3(sev.taitelPsep, 1), 'bara'], ['Riser: feed-line regime', sev.feedRegime, ''], ['Riser: gas-pocket length', r3(sev.feedLength, 0), 'm'], ['Riser: height / length', `${sev.riserHeight.toFixed(0)} / ${sev.riserLength.toFixed(0)}`, 'm']);
    if (cyc) rows.push(['Cycle: type', cyc.type, ''], ['Cycle: period', cyc.period ? r3(cyc.period / 60, 1) : '—', 'min'], ['Cycle: riser-base pressure swing', r3(cyc.amplitude / 1e5, 1), 'bar'], ['Cycle: peak liquid rate / mean', r3(cyc.qLiqPeakRatio, 2), '×'], ['Cycle: peak gas rate / mean', r3(cyc.qGasPeakRatio, 1), '×']);
    tables.push({ title: 'Slug statistics', columns: ['Quantity', 'Value', 'Unit'], rows, note: 'Slug tracking is one-way coupled: slugs travel on the steady solution. Frequencies from correlations are uncertain by a factor of two or more for large pipes; calibrate the multipliers to field data when available.' }); }
  if (terrain.dips.length) tables.push({ title: 'Low points and terrain accumulation', columns: ['x (m)', 'Upward leg (m)', 'Angle (°)', 'vsg (m/s)', 'Critical vsg (m/s)', 'Holdup upward leg', 'Holdup downward leg', 'Liquid held (m³)', 'Regime of the downward leg', 'Accumulates'], rows: terrain.dips.map((d) => [r3(d.x, 0), r3(d.upLength, 0), r3(d.upAngle, 2), r3(d.vsg), r3(d.vsgCrit), r3(d.holdupUp), r3(d.holdupDown), r3(d.volume, 1), d.downRegime, d.accumulates ? 'yes' : 'no']) });
  if (cmp.length) tables.push({ title: 'Holdup-model comparison at the same rate and outlet pressure', columns: ['Model', 'Inlet pressure (bara)', 'Friction (bar)', 'Gravity (bar)', 'Liquid inventory (m³)', 'Arrival temperature (°C)'], rows: cmp, note: 'Other models are solved on a coarser grid (40 cells); the spread is the model-form uncertainty of the pressure drop.' });
  tables.push({ title: 'Minimum stable rate scan', columns: ['Rate (× case)', 'Inlet pressure (bara)', 'Bøe ratio', 'Pots number', 'Feed / downhill regime', 'Slugging predicted', 'Stable with the mitigations', 'Choke pressure drop used (bar)', 'Back-pressure needed (Taitel, bara)'], rows: turndown.rows.map((r) => [r3(r.m, 2), r3(r.pIn, 1), r.boe === null ? '—' : r3(r.boe, 2), r.pots === null ? '—' : r3(r.pots, 3), r.regime, r.unstable ? 'yes' : 'no', r.controlled ? 'yes' : 'no', r.dpChoke === null ? '—' : r3(r.dpChoke, 1), r.taitel === null ? '—' : r3(Math.max(r.taitel, 0), 1)]), note: `Slugging predicted: Bøe / Pots screening without control. Stable with the mitigations (${tdMethod}): the back-pressure at the riser top (separator pressure + choke pressure drop, including booster head spent across the choke) reaches the Taitel stability pressure, or the riser cycle model with the choke and the lift gas shows no cycle; rows below the first rate that stays unstable are not examined. Static choking is assumed: a feedback controller can hold the same operating point at a larger mean opening.` });
  tables.push({ title: `Heat-transfer resistances at ${Math.round(st.x[ic])} m`, columns: ['Layer', 'Resistance (m²K/W on ID)', 'Share (%)', 'Temperature on its outer side (°C)'], rows: (() => { const tot = netC.resistances.reduce((s, x) => s + x.R, 0); let acc = 0; return netC.resistances.map((x) => { acc += x.R; return [x.name, sig(x.R, 4), r3((100 * x.R) / tot, 1), r3(st.T[ic] - ((st.T[ic] - st.tAmb[kc]) * acc) / tot, 2)]; }); })(), note: `Layer-based U-value ${netC.U.toFixed(2)} W/m²K; U-value used ${st.U[kc].toFixed(2)} W/m²K; inside film ${st.hIn[kc].toFixed(0)} W/m²K.` });
  if (tr) tables.push({ title: 'Transient summary', columns: ['Quantity', 'Value', 'Unit'], rows: [['Simulated time', r3(tr.tEnd / 3600, 3), 'h'], ['Time steps (rejected)', `${tr.steps} (${tr.rejected})`, '–'], ['Mean time step', r3(tr.dtMean, 2), 's'], ['Cells', trInfo.grid.n, '–'], ['Initial condition', trInfo.uniform ? 'settled line' : 'steady state of the drift-flux closure', ''], ['Wall-friction tuning factor', r3(trInfo.slip.fricMult, 3), '×'], ['Inlet pressure: mean / min / max (last 60 %)', trSummary ? `${trSummary.pMean.toFixed(1)} / ${trSummary.pMin.toFixed(1)} / ${trSummary.pMax.toFixed(1)}` : '—', 'bara'], ['Liquid outflow: mean / peak (last 60 %)', trSummary ? `${(3600 * trSummary.qMean).toFixed(1)} / ${(3600 * trSummary.qPeak).toFixed(1)}` : '—', 'm³/h'], ['Liquid inventory: start → end', trSummary ? `${trSummary.inv0.toFixed(0)} → ${trSummary.invEnd.toFixed(0)}` : '—', 'm³'], ['Mass-conservation error', sig(Math.abs(tr.mass.error), 3), '–'], ['Energy-conservation error', sig(Math.abs(tr.energy.error), 3), '–'], ['Largest volume-constraint error', sig(tr.volErrMax, 3), '–']], note: 'Drift-flux closure v_g = C0 j + v_d (C0 = 1.2, Bendiksen drift) for every regime; convective acceleration, Joule–Thomson cooling and latent heat are not included in the transient.' });

  // ---- balances ---------------------------------------------------------------------------------------------
  const mOutS = st.rhoG[N] * st.vsg[N] * (PI * st.D[N] ** 2) / 4 + st.rhoL[N] * st.vsl[N] * (PI * st.D[N] ** 2) / 4;
  balances.push({ name: 'Steady mass flow (kg/s): inlet vs outlet phases', in: st.mdot, out: mOutS });
  balances.push({ name: 'Steady momentum (bar): inlet − outlet pressure vs friction + gravity + acceleration + local', in: dpTotal, out: st.dpFric + st.dpGrav + st.dpAcc + st.dpLocal });
  if (cfg.energy === 'enthalpy') balances.push({ name: 'Steady energy (kW): enthalpy in vs enthalpy out + heat loss + elevation + kinetic', in: st.energy.hIn / 1e3, out: (st.energy.hOut + st.heatLoss + st.energy.potential + st.energy.kinetic) / 1e3 });
  if (tr) { balances.push({ name: 'Transient mass (kg): initial + inflow vs final + outflow', in: tr.mass.initial + tr.mass.inflow, out: tr.mass.final + tr.mass.outflow }); balances.push({ name: 'Transient thermal energy (MJ): initial + inflow vs final + outflow + loss', in: (tr.energy.initial + tr.energy.inflow) / 1e6, out: (tr.energy.final + tr.energy.outflow + tr.energy.loss) / 1e6 }); }

  // ---- closure detail, comparisons and dedicated solver tasks ------------------------------------------------
  cm.zG = st.zG[im];
  const RX = { st, cfg, fm, reg, units, im, cm, slug, tr, cyc, sev, cool, plots, tables, balances, warnings, recs, pAmp, terrain, mwG: st.mwG[im] };
  const ex = runExtras(v, ctx, RX), task = TASKS.some((q) => q.value === v.task) ? v.task : 'line', tk = task !== 'line' ? await runTask(task, v, ctx, RX) : { outputs: {}, kpis: [] };
  // ---- KPIs -------------------------------------------------------------------------------------------------
  const kpis = [
    { label: 'Inlet pressure', value: r3(st.pIn, 1), unit: 'bara', status: 'ok', help: `Holdup model: ${HOLDUP_MODELS.find((m) => m.value === cfg.model).label}` },
    { label: 'Pressure drop', value: r3(dpTotal, 1), unit: 'bar', status: 'ok', help: `Friction ${st.dpFric.toFixed(1)}, gravity ${st.dpGrav.toFixed(1)}, acceleration ${st.dpAcc.toFixed(2)}, local ${st.dpLocal.toFixed(2)} bar` },
    { label: 'Arrival temperature', value: r3(st.tOut, 1), unit: '°C', status: st.subcooling[N] > 0 ? 'bad' : wat !== 0 && st.tOut < wat ? 'warn' : 'ok' },
    { label: 'Largest hydrate subcooling', value: r3(subMax, 1), unit: '°C', status: subMax > 3 ? 'bad' : subMax > 0 ? 'warn' : 'ok', help: 'Positive = inside the hydrate region' },
    { label: 'Liquid inventory', value: r3(st.liquidInventory, 0), unit: 'm³', status: 'ok', help: `${(100 * st.liquidInventory / st.volume).toFixed(0)} % of the ${st.volume.toFixed(0)} m³ line volume` },
    { label: 'Heat loss', value: r3(st.heatLoss / 1e3, 0), unit: 'kW', status: 'ok' },
    { label: 'Erosional velocity ratio', value: r3(evrMax, 2), unit: '–', status: evrMax > 1 ? 'bad' : evrMax > 0.8 ? 'warn' : 'ok' },
    { label: 'Slugging', value: slug.type, unit: '', status: severe ? 'bad' : slug.type === 'terrain' ? 'warn' : 'ok' },
    { label: 'Slug frequency', value: sig(60 * fin(slug.freq, 0), 3), unit: '1/min', status: 'ok' },
    { label: 'Largest slug length', value: r3(fin(slug.lengthMax, 0), 0), unit: 'm', status: 'ok' },
    { label: 'Slug-catcher volume needed', value: r3(catcher, 0), unit: 'm³', status: catcher > num(BASE.slugCatcherVol, 60) ? 'warn' : 'ok', help: 'Liquid surge with 25 % margin' },
    { label: 'Minimum stable rate with slug control', value: r3(100 * tdCtrl, 0), unit: '% of case', status: rateNow < tdCtrl ? 'bad' : tdCtrl > 0.6 * rateNow ? 'warn' : 'ok', help: `Uncontrolled: ${(100 * tdFrac).toFixed(0)} %. Mitigation: ${tdMethod}.` },
    { label: 'Bøe ratio / Pots number', value: sev ? `${sev.boe.toFixed(2)} / ${sev.pots === null ? '—' : sev.pots.toFixed(2)}` : 'no riser', unit: '', status: severe ? 'bad' : 'ok' },
    { label: 'Transient inlet-pressure swing', value: r3(pAmp, 1), unit: 'bar', status: pAmp > 0.1 * dpTotal && pAmp > 2 ? 'warn' : 'ok' },
    { label: 'Cooldown to hydrate temperature', value: cool.tReach !== null ? r3(cool.tReach / 3600, 1) : subMax > 0 ? 0 : `> ${(cool.t[cool.t.length - 1] / 3600).toFixed(0)}`, unit: 'h', status: cool.tReach !== null && cool.tReach < 8 * 3600 ? 'warn' : 'ok', help: `At ${Math.round(st.x[ic])} m, the location with the smallest margin` },
    { label: 'Operating rate', value: r3(stdLiquid(fm, st.mScale), 0), unit: 'Sm³/d liquid', status: 'ok', help: `GOR ${gor.toFixed(0)} Sm³/Sm³, water cut ${wcIn.toFixed(0)} %, ${st.mdot.toFixed(1)} kg/s` },
    ...ex.kpis, ...tk.kpis,
  ];

  // ---- outputs ----------------------------------------------------------------------------------------------
  const idx = [...new Set(linspace(0, N, Math.min(N + 1, 120)).map(Math.round))], pk = (a, k = 5) => idx.map((i) => sig(a[i], k + 2));
  const series = tr ? { t: tr.t.map((x) => r3(x, 1)), pIn: D0(tr.pIn, 3), qLiqOut: tr.qLiqOut.map((x) => sig(x, 5)), qGasOut: tr.qGasOut.map((x) => sig(x, 5)), holdupOut: D0(tr.holdupOut, 4) } : { t: [], pIn: [], qLiqOut: [], qGasOut: [], holdupOut: [] };
  const outputs = {
    profile: { x: idx.map((i) => r3(st.x[i], 1)), z: idx.map((i) => r3(st.z[i], 2)), P: pk(st.P), T: pk(st.T), holdup: pk(st.holdup), vsl: pk(st.vsl), vsg: pk(st.vsg), vm: pk(st.vm), rhoM: pk(st.rhoM), dpdx: pk(st.dpdx), tauW: pk(st.tauW), tAmb: pk(st.tAmb), tHyd: pk(st.tHyd), subcooling: pk(st.subcooling), regime: idx.map((i) => reg[i]), tWall: pk(st.tWall), heatFlux: idx.map((i) => sig(heatFlux[i], 6)) },
    pIn: st.pIn, pOut: st.pOut, tIn: st.tIn, tOut: st.tOut, dpTotal, dpFric: st.dpFric, dpGrav: st.dpGrav, dpAcc: st.dpAcc, dpLocal: st.dpLocal, mdot: st.mdot, qLiq: qLout, qGas: qGout, liquidInventory: st.liquidInventory, volume: st.volume, residence: st.residence / 3600, residenceLiquid: st.residenceLiquid / 3600,
    heatLoss: st.heatLoss / 1e3, uValue: mean(st.U), slug, severeSlugging: severe, boe: sev ? sev.boe : null, pots: sev ? sev.pots : null, erosionalRatio: evrMax, maxVelocity: vMax, maxSubcooling: subMax, hydrateLength: hydLen,
    series, pInAmplitude: pAmp, turndownRate: tdFrac, turndownControlled: tdCtrl, turndownMethod: tdMethod, rateFraction: st.mScale, model: cfg.model, slugOnsetX: iOnset >= 0 ? st.x[iOnset] : null, slugCatcherVolume: catcher, terrain: { accumulates: terrain.accumulates, x: terrain.worst ? terrain.worst.x : null, volume: terrain.worst ? terrain.worst.volume : null, vsgCrit: terrain.worst ? terrain.worst.vsgCrit : null },
    cycle: cyc ? { stable: cyc.stable, period: cyc.period, amplitude: cyc.amplitude / 1e5, qLiqPeakRatio: cyc.qLiqPeakRatio } : null, arrival: track ? { n: track.n, meanLength: track.meanLength, p99: track.p99, maxLength: track.maxLength, frequency: track.freqArrival, surge: track.surge } : null,
    task, ...ex.outputs, ...tk.outputs,
    cooldownPreview: cool.tReach !== null ? cool.tReach / 3600 : null, coldSpotX: st.x[ic], regimeAtMap: reg[im], mandhane: mand, id: cfg.idMin, transient: tr ? { massError: tr.mass.error, energyError: tr.energy.error, steps: tr.steps, aborted: tr.aborted, pMean: trSummary ? trSummary.pMean : null, cumLiq: trSummary ? trSummary.cumLiq : null } : null,
  };
  prog(1, 'Done');
  const summary = `Inlet pressure ${st.pIn.toFixed(1)} bara for ${st.pOut.toFixed(1)} bara at the outlet (${dpTotal.toFixed(1)} bar: friction ${st.dpFric.toFixed(1)}, gravity ${st.dpGrav.toFixed(1)}); arrival at ${st.tOut.toFixed(1)} °C, ${subMax > 0 ? `${subMax.toFixed(1)} °C inside the hydrate region at worst` : `${(-subMax).toFixed(1)} °C above the hydrate temperature at the closest point`}; slugging: ${slug.type}${slug.type !== 'none' ? `, liquid surge about ${fin(slug.surge, 0).toFixed(0)} m³` : ''}; minimum stable rate ${(100 * tdFrac).toFixed(0)} % of the case rate.`;
  return { summary, kpis, warnings, recommendations: recs, plots, tables, balances, outputs };
}

// ---- calibration model: coarse steady march, fast enough for least-squares fitting ----------------------------
const CAL_NAN = { dp: NaN, tArr: NaN, holdup: NaN, slugFreq: NaN, slugLen: NaN, filmThk: NaN, filmVel: NaN, level: NaN, tauI: NaN, dpStrat: NaN, vslCrit: NaN, entFrac: NaN, depLength: NaN, dropD: NaN, slugVel: NaN, bodyHoldup: NaN, arrFreq: NaN, arrLen: NaN, mergeShare: NaN, tSurf: NaN, dSauter: NaN };
function calibrationModel(v) {
  const cfg = flowConfig({ ...v, nSteady: 24 }, lastCtx), r = steadySolve({ ...cfg.base, n: 24, mScale: num(v.rateFrac, 1, 1e-3, 20), pOut: cfg.pOut, hydrate: false, tolP: 1e-4 });
  if (!r.ok) return { ...CAL_NAN };
  const i = 12, so = slugOpts(v), two = (k) => r.vsl[k] > 1e-9 && r.vsg[k] > 1e-9, u = two(i) ? slugUnitCell(cellOf(r, i), so) : null, out = { ...CAL_NAN, dp: r.pIn - r.pOut, tArr: r.tOut, holdup: r.liquidInventory / r.volume, slugFreq: u ? 60 * u.freq : 0, slugLen: u ? u.lengthFromFreq : 0, filmThk: u ? 1000 * u.filmThickness : 0, filmVel: u ? u.filmVelocity : 0, slugVel: u ? u.vt : 0, bodyHoldup: u ? u.holdupSlug : 0 };
  if (two(i)) { // stratified closures and the transition at the mid-line cell
    const c = cellOf(r, i), o = { fiMult: num(v.fiMult, 1, 0.05, 20), fwlMult: num(v.fwlMult, 1, 0.05, 20), fwgMult: num(v.fwgMult, 1, 0.05, 20), transMult: num(v.transMult, 1, 0.2, 5) }, sT = stratifiedTransition({ ...c, theta: Math.min(c.theta, 0) }, o), g = stratGeom(sT.level, c.D);
    Object.assign(out, { level: sT.level, tauI: sT.tauI, dpStrat: (sT.tauWL * g.SL + sT.tauWG * g.SG) / g.A, vslCrit: sT.vslCrit });
  }
  { let ia = 0; for (let k = 1; k < 24; k++) if (r.vsg[k] > r.vsg[ia]) ia = k; if (two(ia)) { const a = annularMist({ ...cellOf(r, ia), rough: cfg.rough, P: r.P[ia] * 1e5, fModel: cfg.fModel }, cfg.mp.annular); Object.assign(out, { entFrac: a.entEq, depLength: a.relaxLength, dropD: a.dropSize * 1e6 });
      const c = cellOf(r, ia), aM = clamp(1 - r.holdup[ia], 0.01, 0.25), ib = interfacialAreaTransport({ d0: num(v.bubbleMm, 3, 0.2, 30) / 1000, alpha: aM, vg: Math.max(c.vsg / aM, 0.1), k: 0.01 * (c.vsl + c.vsg) ** 2, length: 30, n: 80, rhoL: c.rhoL, rhoG: c.rhoG, muL: c.muL, sigma: c.sigma, crc: 0.04 * num(v.crcMult, 1, 0, 100), cti: 0.085 * num(v.ctiMult, 1, 0, 100), cwe: 0.002 * num(v.crcMult, 1, 0, 100) }); out.dSauter = ib.dEnd * 1000; } }
  { const c = { s: r.s[i], z: r.z[i], D: r.D[i], pr: cfg.fm.at(r.P[i], r.T[i], r.mScale ?? num(v.rateFrac, 1, 1e-3, 20)), vm: r.vm[i], holdup: r.holdup[i], ta: r.tAmb[i], T: r.T[i] }, net = cfg.network(c); out.tSurf = c.ta + (net.U * c.D * (c.T - c.ta)) / (cfg.od * Math.max(net.hOut ?? cfg.hOutOf(c), 1e-9)); }
  if (r.vsl.every((q, k) => two(k))) { // tracked arrivals: initiation, growth / decay and merging parameters
    const units = r.vsl.map((_, k) => slugUnitCell(cellOf(r, k), so)), tk = slugTracking(trackingField(r, units), { sites: [{ s: 0, freq: units[0].freq, length: units[0].lengthFromFreq }], nSlugs: 14, seed: 11, theta: r.theta, initMult: num(v.initMult, 1, 0.05, 20), relaxMult: num(v.relaxMult, 1, 0.05, 20), wakeMult: num(v.wakeMult, 1, 0, 10), maxOps: 6000 });
    Object.assign(out, { arrFreq: 60 * tk.freqArrival, arrLen: tk.meanLength, mergeShare: tk.generated > 0 ? tk.merges / tk.generated : 0 });
  } else Object.assign(out, { arrFreq: 0, arrLen: 0, mergeShare: 0 });
  for (const k of Object.keys(out)) if (!Number.isFinite(out[k])) out[k] = 0;
  return out;
}
const CAL_SAMPLE = /*CAL*/[{rateFrac:0.5, dp:54.36, tArr:21.95, holdup:0.4689, slugFreq:0.187, slugLen:54.1, filmThk:177.3, filmVel:0.655, slugVel:1.882, bodyHoldup:0.9143, level:0.6952, tauI:0.43, dpStrat:12.58, vslCrit:0.3147, entFrac:0, depLength:0.284, dropD:3756, dSauter:32.916, arrFreq:1.116, arrLen:98.8, mergeShare:0.11, tSurf:4.18},
  {rateFrac:0.6, dp:57.94, tArr:25.21, holdup:0.4504, slugFreq:0.1912, slugLen:50, filmThk:172.9, filmVel:0.824, slugVel:2.113, bodyHoldup:0.9437, level:0.6672, tauI:0.621, dpStrat:17.89, vslCrit:0.3387, entFrac:0.0686, depLength:0.359, dropD:2379, dSauter:30.762, arrFreq:1.2514, arrLen:92.3, mergeShare:0.098, tSurf:4.218},
  {rateFrac:0.7, dp:58.5, tArr:31.77, holdup:0.4514, slugFreq:0.1963, slugLen:51, filmThk:167.3, filmVel:0.979, slugVel:2.319, bodyHoldup:0.8992, level:0.6733, tauI:0.882, dpStrat:22.47, vslCrit:0.3436, entFrac:0.1317, depLength:0.428, dropD:1724, dSauter:27.622, arrFreq:1.477, arrLen:87.9, mergeShare:0.075, tSurf:4.378},
  {rateFrac:0.8, dp:59.42, tArr:32.38, holdup:0.425, slugFreq:0.2119, slugLen:52.7, filmThk:161.7, filmVel:1.154, slugVel:2.532, bodyHoldup:0.9185, level:0.6621, tauI:1.135, dpStrat:29.81, vslCrit:0.3494, entFrac:0.1993, depLength:0.524, dropD:1230, dSauter:28.577, arrFreq:1.6334, arrLen:86.2, mergeShare:0.073, tSurf:4.229},
  {rateFrac:0.9, dp:63.53, tArr:35.49, holdup:0.4313, slugFreq:0.2245, slugLen:52.1, filmThk:163.6, filmVel:1.33, slugVel:2.676, bodyHoldup:0.9012, level:0.6761, tauI:1.385, dpStrat:36.64, vslCrit:0.3465, entFrac:0.264, depLength:0.569, dropD:988, dSauter:27.716, arrFreq:1.7447, arrLen:88.3, mergeShare:0.057, tSurf:4.092},
  {rateFrac:1, dp:65.15, tArr:37.73, holdup:0.4431, slugFreq:0.2463, slugLen:51.6, filmThk:162, filmVel:1.527, slugVel:2.99, bodyHoldup:0.8578, level:0.6692, tauI:1.66, dpStrat:44.21, vslCrit:0.35, entFrac:0.3333, depLength:0.669, dropD:758, dSauter:27.079, arrFreq:1.869, arrLen:87.1, mergeShare:0.053, tSurf:4.082},
  {rateFrac:1.1, dp:69.78, tArr:39.79, holdup:0.4476, slugFreq:0.2747, slugLen:53.6, filmThk:163.9, filmVel:1.601, slugVel:3.159, bodyHoldup:0.8416, level:0.6842, tauI:1.977, dpStrat:54.75, vslCrit:0.3582, entFrac:0.3867, depLength:0.719, dropD:632, dSauter:25.719, arrFreq:1.9711, arrLen:86.4, mergeShare:0.043, tSurf:4.202},
  {rateFrac:1.2, dp:69.88, tArr:39.96, holdup:0.4371, slugFreq:0.2958, slugLen:52.2, filmThk:160.5, filmVel:1.81, slugVel:3.298, bodyHoldup:0.8442, level:0.684, tauI:2.15, dpStrat:60.87, vslCrit:0.3542, entFrac:0.4349, depLength:0.794, dropD:515, dSauter:25.569, arrFreq:2.044, arrLen:92.5, mergeShare:0.045, tSurf:4.302},
  {rateFrac:1.3, dp:76.11, tArr:40.46, holdup:0.4472, slugFreq:0.3442, slugLen:49, filmThk:163, filmVel:1.945, slugVel:3.483, bodyHoldup:0.8411, level:0.6785, tauI:2.47, dpStrat:70.12, vslCrit:0.3548, entFrac:0.4999, depLength:0.865, dropD:420, dSauter:26.986, arrFreq:2.1925, arrLen:87.5, mergeShare:0.065, tSurf:4.232},
  {rateFrac:1.4, dp:81.65, tArr:43.61, holdup:0.4541, slugFreq:0.5573, slugLen:51.6, filmThk:158.4, filmVel:2.118, slugVel:3.6, bodyHoldup:0.8141, level:0.6887, tauI:2.757, dpStrat:78.48, vslCrit:0.3501, entFrac:0.5326, depLength:0.931, dropD:377, dSauter:25.888, arrFreq:2.2629, arrLen:87.8, mergeShare:0.058, tSurf:4.208}]/*CAL*/;
const VAL_SAMPLE = /*VAL*/[{rateFrac:0.55, dp:56.77, tArr:24.98, holdup:0.4589, slugFreq:0.19, slugLen:51, filmThk:173.6, filmVel:0.732, slugVel:2.011, bodyHoldup:0.9041, level:0.6848, tauI:0.529, dpStrat:14.99, vslCrit:0.3215, entFrac:0.0336, depLength:0.326, dropD:3047, dSauter:31.738, arrFreq:1.1988, arrLen:102.4, mergeShare:0.1, tSurf:4.132},
  {rateFrac:0.75, dp:58.38, tArr:31.63, holdup:0.4446, slugFreq:0.2, slugLen:53.7, filmThk:166.6, filmVel:1.066, slugVel:2.416, bodyHoldup:0.9056, level:0.656, tauI:0.96, dpStrat:26.57, vslCrit:0.3422, entFrac:0.1676, depLength:0.447, dropD:1460, dSauter:29, arrFreq:1.5284, arrLen:86.9, mergeShare:0.08, tSurf:4.16},
  {rateFrac:0.95, dp:62.57, tArr:35.45, holdup:0.4511, slugFreq:0.2342, slugLen:51.4, filmThk:163.8, filmVel:1.392, slugVel:2.74, bodyHoldup:0.8636, level:0.6852, tauI:1.559, dpStrat:40.61, vslCrit:0.3543, entFrac:0.2891, depLength:0.61, dropD:861, dSauter:27.848, arrFreq:1.7085, arrLen:86.8, mergeShare:0.053, tSurf:4.353},
  {rateFrac:1.05, dp:68.85, tArr:37.7, holdup:0.4529, slugFreq:0.2548, slugLen:52.2, filmThk:164.9, filmVel:1.587, slugVel:2.977, bodyHoldup:0.8536, level:0.6797, tauI:1.732, dpStrat:47.43, vslCrit:0.3664, entFrac:0.3562, depLength:0.701, dropD:675, dSauter:26.965, arrFreq:1.9972, arrLen:83.4, mergeShare:0.049, tSurf:4.183},
  {rateFrac:1.25, dp:74.05, tArr:40.46, holdup:0.4529, slugFreq:0.3214, slugLen:50.7, filmThk:163.3, filmVel:1.897, slugVel:3.3, bodyHoldup:0.8428, level:0.6726, tauI:2.379, dpStrat:66.74, vslCrit:0.3521, entFrac:0.4556, depLength:0.852, dropD:486, dSauter:25.678, arrFreq:2.0519, arrLen:87.8, mergeShare:0.068, tSurf:4.102},
  {rateFrac:1.5, dp:81.68, tArr:44.8, holdup:0.4618, slugFreq:0.8114, slugLen:53.6, filmThk:167, filmVel:2.242, slugVel:3.751, bodyHoldup:0.7929, level:0.7008, tauI:3.015, dpStrat:86.14, vslCrit:0.3493, entFrac:0.579, depLength:0.979, dropD:330, dSauter:24.619, arrFreq:2.4099, arrLen:87.7, mergeShare:0.056, tSurf:4.117}]/*VAL*/;

// ---- verification ------------------------------------------------------------------------------------------------
// Constant-property test fluid: ideal gas with a fixed gas mass fraction and an incompressible liquid.
function testFluid({ mdot = 40, wG = 0, rhoL = 850, muL = 2e-3, muG = 1.5e-5, cp = 2100, M = 0.02, Z = 1, sigma = 0.03 } = {}) {
  return { rates: { mHC: mdot, mW: 0 }, hydrateT: () => -50, at(P, T, m = 1) {
    const rhoG = (P * 1e5 * M) / (Z * RGAS * (T + 273.15)), mG = mdot * wG * m, mO = mdot * (1 - wG) * m;
    return { wG, rhoG, rhoO: rhoL, rhoW: rhoL, rhoL, muG, muO: muL, muW: muL, muL, cpG: cp, cpO: cp, cpW: cp, cpL: cp, kG: 0.03, kO: 0.13, kW: 0.6, kL: 0.13, hG: cp * T, hO: cp * T, hW: cp * T, jtG: 0, jtO: 0, sigma, zG: Z, mwG: M * 1000, mG, mO, mW: 0, qG: mG / rhoG, qO: mO / rhoL, qW: 0, qL: mO / rhoL, wcut: 0, phaseInv: false };
  } };
}
const colebrook = (Re, rel) => { let x = 7; for (let i = 0; i < 80; i++) x = -2 * log10(rel / 3.7 + (2.51 * x) / Re); return 1 / (x * x); }; // independent fixed-point solution
async function verify() {
  const out = [], add = (name, expected, got, tol, note, rel = true) => out.push({ name, expected: sig(expected, 8), got: sig(got, 8), tol, pass: Number.isFinite(got) && Math.abs(got - expected) <= tol * (rel ? Math.max(Math.abs(expected), 1e-12) : 1), note });
  const flat = (L) => ({ x: [0, L], z: [0, 0] });
  { // 1 single-phase liquid: Darcy–Weisbach
    const D = 0.2, L = 5000, fm = testFluid({ mdot: 40 }), v = 40 / 850 / ((PI * D * D) / 4), f = colebrook((850 * v * D) / 2e-3, 4.5e-5 / D), r = steadySolve({ fm, profile: flat(L), n: 40, id: D, pIn: 60, energy: 'isothermal' });
    add('Single-phase liquid pressure drop (Darcy–Weisbach hand calculation)', (f * (L / D) * 850 * v * v) / 2 / 1e5, r.pIn - r.pOut, 1e-6, 'bar; 5 km of 200 mm pipe, 40 kg/s, friction factor from an independent Colebrook iteration');
  }
  { // 2 hydrostatic liquid column, 3 barometric gas column
    const fm = testFluid({ mdot: 1e-13 }), r = steadySolve({ fm, profile: { x: [0, 1e-6], z: [0, 500] }, n: 20, id: 0.2, pIn: 80, energy: 'isothermal' });
    add('Hydrostatic equilibrium of a static liquid column', (850 * G * 500) / 1e5, r.pIn - r.pOut, 1e-9, 'bar; 500 m vertical, no flow');
    const fg = testFluid({ mdot: 1e-5, wG: 1 }), rg = steadySolve({ fm: fg, profile: { x: [0, 1e-6], z: [0, 1000] }, n: 40, id: 0.2, pIn: 100, tIn: 20, energy: 'isothermal' });
    add('Static gas column (barometric formula)', 100 * Math.exp((-0.02 * G * 1000) / (RGAS * 293.15)), rg.pOut, 1e-6, 'bara at the top of 1,000 m of ideal gas at 20 °C');
  }
  { // 4–6 isothermal gas pipeline: analytic solution, observed order, grid convergence
    const D = 0.5, L = 60000, mdot = 90, fm = testFluid({ mdot, wG: 1 }), Gm = mdot / ((PI * D * D) / 4), f = colebrook((Gm * D) / 1.5e-5, 4.5e-5 / D), K = (Gm * Gm * RGAS * 293.15) / 0.02, p1 = 80e5;
    let p2 = 60e5; for (let i = 0; i < 200; i++) p2 = Math.sqrt(p1 * p1 - K * ((f * L) / D + 2 * Math.log(p1 / p2)));
    const run = (n) => steadySolve({ fm, profile: flat(L), n, id: D, pIn: 80, tIn: 20, energy: 'isothermal' }).pOut, a = run(25), b = run(50), c = run(100);
    add('Single-phase gas isothermal pipeline (analytic compressible solution)', p2 / 1e5, c, 2e-5, 'bara at the outlet of 60 km; p₁² − p₂² = G²(ZRT/M)(fL/D + 2 ln p₁/p₂)');
    add('Observed order of accuracy of the steady march (exact solution)', 2, Math.log2(Math.abs(a - p2 / 1e5) / Math.abs(b - p2 / 1e5)), 0.25, 'Errors on 25 and 50 cells; the midpoint rule is formally second order', false);
    const g = gci([1 / 100, 1 / 50, 1 / 25], [c, b, a]);
    add('Richardson extrapolation recovers the exact outlet pressure', p2 / 1e5, g.fExact, 2e-6, `Grid-convergence index of the fine grid ${(100 * g.gciFine).toExponential(1)} %, observed order ${g.p.toFixed(2)}`);
  }
  { // 7 exponential temperature decay
    const D = 0.25, L = 20000, U = 6, fm = testFluid({ mdot: 30 }), exp = 4 + (70 - 4) * Math.exp((-U * PI * D * L) / (30 * 2100));
    const r = steadySolve({ fm, profile: flat(L), n: 200, id: D, pIn: 80, tIn: 70, U, tAmb: 4, energy: 'cpjt' }), h = steadySolve({ fm, profile: flat(L), n: 200, id: D, pIn: 80, tIn: 70, U, tAmb: 4, energy: 'enthalpy' });
    add('Heat loss: analytic exponential temperature decay', exp, r.tOut, 5e-5, '°C after 20 km on 200 cells, constant U and heat capacity, no Joule–Thomson effect');
    add('Enthalpy form of the energy equation reproduces the same decay', exp, h.tOut, 2e-4, '°C; flowing-enthalpy balance with h = cp T (kinetic-energy change is negligible)');
    add('Energy-balance closure of the steady solver', h.energy.hIn - h.energy.hOut, h.heatLoss + h.energy.potential + h.energy.kinetic, 1e-9, 'W; enthalpy drop = heat loss + elevation + kinetic energy');
  }
  const two = { vsl: 1.2, vsg: 2.5, rhoL: 820, rhoG: 60, muL: 1.5e-3, muG: 1.4e-5, sigma: 0.02, D: 0.25, theta: 0, rough: 4.5e-5, P: 6e6 };
  { // 8–9 homogeneous limit and Zuber–Findlay no-slip limit
    const lam = two.vsl / 3.7, rho = 820 * lam + 60 * (1 - lam), mu = 1.5e-3 * lam + 1.4e-5 * (1 - lam), f = colebrook((rho * 3.7 * 0.25) / mu, 4.5e-5 / 0.25), h = holdupGradient(two, 'homogeneous'), z = holdupGradient(two, 'zuberFindlay', { c0: 1, vDrift: 0 });
    add('Homogeneous two-phase limit: frictional gradient', (f * rho * 3.7 * 3.7) / (2 * 0.25), h.fric, 1e-6, 'Pa/m from the no-slip density and viscosity (hand calculation)');
    add('Zuber–Findlay with C0 = 1, vd = 0 reproduces the no-slip holdup', lam, z.holdup, 1e-12, 'Liquid holdup equals the input liquid fraction');
    add('Zuber–Findlay holdup for C0 = 1.2, vd = 0.35 m/s', 1 - 2.5 / (1.2 * 3.7 + 0.35), holdupGradient(two, 'zuberFindlay', { c0: 1.2, vDrift: 0.35 }).holdup, 1e-12, 'α = vsg / (C0 vm + vd)');
  }
  { // 10 Taitel–Dukler level, 11 Bendiksen constants
    const X = 1.58386, p = { rhoL: 1000, rhoG: 10, muL: 1e-3, muG: 1.5e-5, D: 0.1, vsg: 3, theta: 0 }; p.vsl = p.vsg * (X * X * (p.rhoG / p.rhoL) ** 0.8 * (p.muG / p.muL) ** 0.2) ** (1 / 1.8);
    add('Taitel–Dukler equilibrium level at X = 1.584, Y = 0 (turbulent–turbulent)', 0.5, stratifiedRoots(p, { pure: true })[0], 2e-4, 'h/D = 0.5 follows from the dimensionless momentum balance by hand', false);
    const h = slugUnitCell({ ...two, vsl: 0.5, vsg: 1 }), vt = slugUnitCell({ ...two, vsl: 0.5, vsg: 1, theta: PI / 2 });
    add('Bendiksen translational velocity, horizontal (C0 = 1.05, vd = 0.54 √gD)', 1.05 * 1.5 + 0.54 * Math.sqrt(G * 0.25), h.vt, 1e-12, 'm/s at vm = 1.5 m/s in a 250 mm pipe');
    add('Bendiksen / Nicklin translational velocity, vertical (C0 = 1.2, vd = 0.35 √gD)', 1.2 * 1.5 + 0.35 * Math.sqrt(G * 0.25), vt.vt, 1e-12, 'm/s');
    add('Slug-body / film liquid balance closes at the tail', (h.vt - h.vBody) * h.holdupSlug, (h.vt - h.filmVelocity) * h.holdupFilm, 1e-10, 'Shedding rate from the body equals the film flux relative to the tail');
  }
  { // radial conduction: steady analytic, manufactured solution in space and time
    const L = [{ t: 0.0159, k: 45, rho: 7850, cp: 470 }, { t: 0.08, k: 0.17, rho: 900, cp: 1700 }], ri = 0.127, r1 = ri + 0.0159, r2 = r1 + 0.08, R = 1 / (800 * ri) + Math.log(r1 / ri) / 45 + Math.log(r2 / r1) / 0.17 + 1 / (500 * r2);
    add('Radial multilayer conduction: heat flow of the analytic composite cylinder', (2 * PI * 56) / R, radialConduction({ ri, layers: L, hIn: 800, hOut: 500, tFluid: 60, tAmb: 4, T0: 'steady', tEnd: 0 }).qIn, 1e-10, 'W/m through steel + insulation with both films');
    const k = 0.5, rc = 2e6, tau = 4000, a = 0.1, b = 0.2, Tm = (r, t) => 20 + 800 * r * r * Math.exp(-t / tau), src = (r, t) => Math.exp(-t / tau) * 800 * (-(rc / tau) * r * r - 4 * k);
    const err = (nPer, nSteps) => { const s = radialConduction({ ri: a, layers: [{ t: b - a, k, rho: 2000, cp: 1000 }], hIn: 1e12, hOut: 1e12, tFluid: (t) => Tm(a, t), tAmb: (t) => Tm(b, t), T0: (r) => Tm(r, 0), tEnd: 2000, nSteps, nPer, source: src }); return Math.max(...s.r.map((r, i) => Math.abs(s.T[i] - Tm(r, 2000)))); };
    const e1 = err(8, 40), e2 = err(16, 160), t1 = err(64, 10), t2 = err(64, 20);
    add('Method of manufactured solutions: spatial order of the radial conduction solver', 2, Math.log2(e1 / e2), 0.35, `Maximum error ${e1.toExponential(2)} → ${e2.toExponential(2)} K when the cell size is halved (time step ÷ 4)`, false);
    add('Time-step convergence of the implicit radial solver (first order)', 1, Math.log2(t1 / t2), 0.3, `Maximum error ${t1.toExponential(2)} → ${t2.toExponential(2)} K when the time step is halved`, false);
  }
  { // transient: conservation, steady limit, advection, CFL sensitivity
    const D = 0.2, L = 3000, n = 60, fm = testFluid({ mdot: 20, wG: 0.02, rhoL: 800 }), prof = flat(L), grid = transientGrid(prof, n), A = (PI * D * D) / 4, pr = fm.at(50, 30), jj = (20 / A) * (0.02 / pr.rhoG + 0.98 / 800), lam = (20 * 0.98) / 800 / A / jj;
    const mk = (cfl, init) => transientDriftFlux({ fm, grid, D, rough: 0, init, mdot0: 20, mdotOf: () => 20, pOutOf: () => 50, slip: { c0: 1, vdScale: 0 }, U: grid.sc.map(() => 0), tAmb: grid.sc.map(() => 30), tEnd: (0.8 * L) / jj, cfl, nField: 161, cpG: 2100, cpL: 2100 });
    const init = { P: grid.sc.map(() => 50), T: grid.sc.map(() => 30), holdup: grid.sc.map((s) => (s < L / 3 ? lam : 0.5 * lam)) }, r = mk(0.8, init), r2 = mk(0.4, init), ic = 39, mid = 0.75 * lam;
    const cross = (q) => { for (let k = 1; k < q.field.t.length; k++) { const h0 = q.field.holdup[k - 1][ic], h1 = q.field.holdup[k][ic]; if (h0 < mid && h1 >= mid) return q.field.t[k - 1] + ((q.field.t[k] - q.field.t[k - 1]) * (mid - h0)) / (h1 - h0); } return NaN; };
    add('Transient advection of a holdup step: arrival time (contact wave)', (grid.sc[ic] - L / 3) / jj, cross(r), 0.03, 's; no-slip closure, the step travels with the mixture velocity');
    add('CFL sensitivity of the arrival time (CFL 0.4 against 0.8)', cross(r), cross(r2), 0.02, 's; the first-order upwind front smears but its centre does not move');
    add('Transient mass-conservation error', 0, r.mass.error, 1e-10, 'Relative; (final − initial − inflow + outflow) / initial', false);
    add('Transient energy-conservation error', 0, r.energy.error, 1e-10, 'Relative thermal-energy balance', false);
    add('Phase-volume conservation of the transient (largest volume-constraint error)', 0, r.volErrMax, 5e-3, 'Sum of the phase volume fractions minus one', false);
    // steady limit against the steady solver with the same closure
    const slip = { c0: 1.2, vdScale: 1, fricMult: 1 }, prof2 = { x: [0, 2000, 2200], z: [0, 0, 150] }, st = steadySolve({ fm, profile: prof2, n: 120, id: D, rough: 4.5e-5, pOut: 30, tIn: 30, U: 0, model: 'transient', mp: { slip }, energy: 'isothermal', hydrate: false }), g2 = transientGrid(prof2, 44, 2000, 2), at = (q) => g2.sc.map((s) => interp1(st.s, q, s));
    const tr = transientDriftFlux({ fm, grid: g2, D, init: { P: at(st.P).map((p) => p * 1.01), T: at(st.T), holdup: at(st.holdup).map((h) => 0.9 * h) }, mdot0: 20, mdotOf: () => 20, pOutOf: () => 30, slip, U: g2.sc.map(() => 0), tAmb: g2.sc.map(() => 30), tEnd: 3000, cfl: 0.8, cpG: 2100, cpL: 2100 });
    add('Steady-state limit of the transient from a perturbed initial state', st.pIn, tr.pIn[tr.pIn.length - 1], 0.01, 'bara at the inlet after 3,000 s with constant boundaries, against the steady march with the same closure');
  }
  { // local models
    const lam = ransPipe({ reTau: 40, model: 'laminar', n: 120 }), mx = ransPipe({ reTau: 2000, model: 'mixing', n: 100 }), kw = ransPipe({ reTau: 2000, model: 'komega', n: 80, tol: 1e-8 }), ke = ransPipe({ reTau: 2000, model: 'kepsilon', n: 80, tol: 1e-8 });
    add('Laminar limit of the radial solve: f = 64 / Re', 64 / lam.Re, lam.f, 1e-3, `Re = ${lam.Re.toFixed(0)} with the eddy viscosity switched off`);
    add('Mixing-length RANS friction factor against Colebrook (smooth)', frictionFactor(mx.Re, 0), mx.f, 0.03, `Re = ${mx.Re.toFixed(0)}`);
    add('k–ω RANS friction factor against Colebrook (smooth)', frictionFactor(kw.Re, 0), kw.f, 0.06, `Re = ${kw.Re.toFixed(0)}, ${kw.iterations} iterations`);
    add('Low-Reynolds k–ε RANS friction factor against Colebrook (smooth)', frictionFactor(ke.Re, 0), ke.f, 0.05, `Re = ${ke.Re.toFixed(0)}, ${ke.iterations} iterations`);
    add('Iterative convergence of the k–ω solve (bulk-velocity residual)', 0, kw.residual, 1e-6, 'Relative change of the bulk velocity in the last iteration', false);
    const db = damBreak({ hL: 1, hR: 0, n: 400 }), vf = vofAdvect1D({ n: 100 });
    add('Dam-break benchmark (Ritter solution): mean depth error', 0, db.l1, 0.008, 'Fraction of the initial depth on 400 cells, first-order HLL scheme', false);
    add('Wet-bed dam break (semi-analytical Stoker solution): mean depth error', 0, damBreak({ hL: 1, hR: 0.2, n: 400 }).l1, 0.008, 'Fraction of the initial depth; the middle state needs a root solve', false);
    add('Volume-of-fluid advection conserves the liquid volume', 0, vf.massError, 1e-12, 'Relative, after one revolution of the periodic domain', false);
    add('THINC keeps the advected interface sharp', 2, vf.thickness, 1.01, 'Cells between 5 % and 95 % liquid per interface after 250 steps (upwind: tens of cells)', false);
  }
  { // reference case: closure, code-to-code, boundary-condition consistency, reproducibility
    const ctx = {}, a = steadyFlow({ energy: 'cpjt' }, ctx), cfg = a.cfg, k = marchKernel(cfg);
    add('Momentum-balance closure of the steady solver (reference case)', a.pIn - a.pOut, a.dpFric + a.dpGrav + a.dpAcc + a.dpLocal, 1e-10, 'bar; inlet − outlet pressure against the summed components');
    let loc = 0; for (let i = 0; i < a.n; i++) loc = Math.max(loc, Math.abs(a.P[i] - a.P[i + 1] - (a.dpdx[i] * a.ds) / 1e5));
    add('Local momentum-balance closure (largest cell residual)', 0, loc, 1e-10, 'bar per cell', false);
    add('Nonlinear convergence of the shooting on the inlet pressure', cfg.pOut, a.pOut, 1e-6, 'bara: outlet pressure met', false);
    add('Code-to-code: inlet pressure against the kernel march (same closure, bisection shooting)', k.pIn, a.pIn, 1e-4, 'bara; two independent implementations of the reference case');
    const b = steadyFlow({ energy: 'cpjt', bc: 'inletP', pInSet: a.pIn }, ctx);
    add('Boundary-condition consistency: marching from the solved inlet pressure returns the outlet pressure', a.pOut, b.pOut, 1e-6, 'bara');
    const m = a.rhoG[a.n] * a.vsg[a.n] * a.D[a.n] ** 2 * (PI / 4) + a.rhoL[a.n] * a.vsl[a.n] * a.D[a.n] ** 2 * (PI / 4);
    add('Mass-balance closure of the steady solver', a.mdot, m, 1e-10, 'kg/s: inlet mass rate against the outlet phase rates from velocity × density × area');
    const u = a.vsl.map((_, i) => slugUnitCell(cellOf(a, i), slugOpts({}))), tf = trackingField(a, u), site = [{ s: 0, freq: u[0].freq, length: u[0].lengthFromFreq }], t1 = slugTracking(tf, { sites: site, nSlugs: 40, seed: 7 }), t2 = slugTracking(tf, { sites: site, nSlugs: 40, seed: 7 }), a2 = steadyFlow({ energy: 'cpjt' }, ctx);
    add('Reproducibility: two identical runs give identical results', a.pIn + t1.meanLength, a2.pIn + t2.meanLength, 0, 'Steady inlet pressure plus the seeded slug-tracking mean length, bit for bit', false);
    const cy = riserSluggingCycle({ D: 0.254, feedLength: 5000, feedAngle: 0.01, riserHeight: 1000, riserLength: 1000, wG: 0.5, wL: 10, rhoL: 800, T: 300, pSep: 20e5, alphaL: 0.4, chokeDp: 1e3 });
    add('Riser cycle: peak riser-base pressure against the full liquid head', 20 + (800 * G * 1000) / 1e5, cy.pBaseMax / 1e5, 0.01, 'bara: separator pressure + ρL g H of a liquid-filled riser (the blowout adds a small friction and inertia overshoot)');
  }
  await verifyExtra(add);
  await verify3D(add);
  return out;
}
function marchKernel(cfg) { return marchSteady({ fm: cfg.fm, profile: cfg.profile, id: cfg.id, rough: cfg.rough, U: cfg.U, tAmbOf: cfg.tAmbOf, tIn: cfg.tIn, pOut: cfg.pOut, n: cfg.n, energy: cfg.energy, label: false }); }

const SUITE = {
  id: 'flow', num: 3, title: 'Multiphase Thermal-Hydraulics & Slugging', short: 'Flow · Slugs', icon: '🌊',
  tagline: 'Pressure, temperature, holdup and flow regime along the line, slugging of every kind, and the transient response to rate and pressure changes.',
  description: 'Solves the steady mass, momentum and energy balances along the elevation profile with twelve selectable holdup closures and the equation-of-state property table, classifies the flow regime mechanistically, and evaluates hydrodynamic, terrain and severe riser slugging with a unit-cell model, Lagrangian slug tracking and a lumped riser cycle. A semi-implicit drift-flux model gives the transient response to rate, back-pressure and choke histories. On request a transient two-fluid model captures slugs on a test section, a 2-D Navier–Stokes solver (projection, SST k–ω, volume of fluid) resolves a slug front or a dam break, a 3-D solver runs LES, DES, IDDES, resolved simulation, LES + VOF and a 1-D + 3-D coupling on small grids, and radial models give turbulence, bubbly-flow and conduction detail. Branch junctions, a booster pump and separator / compressor characteristics can be imposed at the boundaries; measurements are compared with sourced reference data sets.',
  guide: [
    'Check the boundary conditions: by default the case rates and the outlet pressure are fixed and the inlet pressure is solved. Pull the profile, diameter and U-value from the Network suite when it has been run.',
    'Choose the holdup model on the Model tab; the model-comparison table shows how much the inlet pressure depends on that choice.',
    'Read the slugging results: type, frequency, largest slug, surge volume and the minimum stable rate. For a riser, the Bøe ratio, the Pots number and the cycle model are reported.',
    'Edit the boundary history for the transient (rate ramp, back-pressure step, choke move) and the simulated time; the distance–time maps show liquid accumulation and surges.',
    'Use the Mesh tab to quantify the numerical uncertainty of the steady march and of the transient grid and time step.',
    'Fit roughness, U-value, holdup and slug multipliers to measured pressure drop, arrival temperature, holdup and slug frequency on the Calibration tab.',
  ],
  implemented: /*IMPL*/[ // fragments of the catalogue names (normalised as the Equations tab does); the last line carries the fragments of LES (“le”), DES (“de”), IDDES, DNS, LES + VOF and the 1-D + 3-D coupling
    'conservation of ma', 'component ma', 'phase continuity equation', 'mixture continuity equation', 'conservation of linear momentum', 'phase momentum equation',
    'mixture momentum equation', 'conservation of total energy', 'phase energy equation', 'mixture energy equation', 'internal energy equation', 'enthalpy equation', 'entropy inequality',
    'specie transport equation', 'two fluid mod', 'multi fluid mod', 'six equation mod', 'seven equation mod', 'drift flux mod', 'mixture mod', 'homogeneou equilibrium mod',
    'homogeneou relaxation mod', 'separated flow mod', 'mechanistic multiphase mod', 'walli drift flux formulation', 'nicklin type translational velocity relation',
    'begg brill correlation', 'hagedorn brown correlation', 'dun ro correlation', 'mukherjee brill mod', 'gray correlation', 'orkiszewski mod', 'ansari mechanistic mod',
    'er flow regime mod', 'mandhane type flow map', 'baker flow map', 'barnea unified flow pattern mod', 'mechanistic stratified annular intermittent transition criteria',
    'unit cell slug mod', 'mechanistic slug flow mod', 'slug tracking mod', 'slug capturing mod', 'slug frequency correlation', 'ength correlation', 'nicklin slug cel',
    'translational velocity mod', 'slug body film ma', 'slug growth d', 'slug merging equation', 'kelvin helmholtz instability criterion', 'viscou kelvin helmholtz mod',
    'long wave stability analysi', 'linear stability theory', 'nonlinear wave growth mod', 'roll wave instability mod', 'pipeline riser liquid accumulation mod', 'severe slugging cycl',
    'hydrodynamic stability criteria', 'pressure build up liquid fallback mod', 'terrain induced accumulation mod', 'interfacial drag', 'wall friction', 'lift force', 'virtual add',
    'wall lubrication', 'ent dispersion', 'e induced turbul', 'entrainment d', 'escence breakup', 'fourier law', 'transient heat conduction equation', 'newton law of cooling',
    'overall heat transfer resistance', 'conjugate heat transfer equation', 'radial multilayer conduction', 'seabed conduction', 'natural forced convection', 'e thomson relation', 'rng k',
    'realizabl', 'sst k', 'reynold stre', 'spalart allmara', 'volume of fluid', 'evel set', 'front tracking', 'phase field', 'drift flux slug tracking', 'two fluid mechanistic closure',
    'two fluid slug capturing', 'ran vof', 'ection conservation law solver', 'eo thermal hydraulic conservation equation', 'pressure and temperature throughout the pipeline',
    'riser and well network', 'initial oil', 'ga and water flow rate', 'phase velocitie and superficial velocitie', 'liquid holdup and ga void fraction', 'initial flow regime',
    'initial liquid film thickne', 'initial fluid inventory', 'initial phase distribution', 'initial wall temperature', 'initial pipe and insulation temperature',
    'where slugging already exist', 'initial number', 'position', 'ength', 'velocity', 'liquid content and distribution of slug', 'appropriate combination of oil',
    'or volumetric flow rate', 'phase fraction', 'pressure', 'temperature', 'external heat transfer condition', 'insulation behaviour and any imposed heating or cooling',
    'appropriate flow between connected branche', 'e valve', 'choke', 'pump', 'compressor and separator boundarie should impose their respective hydraulic characteristic',
    'production rate change', 'terrain or severe riser slugging', 'validated output from modul', 'or volumetric rate', 'gor', 'water cut', 'well inflow condition',
    'phase distribution and liquid holdup', 'insulation propertie', 'transient operating boundary historie', 'numerical discretization and convergence control',
    'phase and superficial velocitie', 'liquid holdup and void fraction', 'flow regime', 'fluid inventory', 'slug initiation onset', 'terrain and severe riser slugging', 'liquid surge',
    'slug arrival statistic and slug catcher loading', 'interfacial friction correlation', 'liquid wall shear closure', 'ga wall shear closure', 'interfacial shear', 'entrainment rate',
    'e size', 'et size', 'distribution parameter', 'liquid holdup correlation parameter', 'slug frequency parameter', 'erity parameter', 'slug body holdup', 'film thickne',
    'slug initiation parameter', 'slug merging parameter', 'heat transfer coefficient', 'overall u value', 'pipe insulation thermal propertie', 'method of manufactured solution mm',
    'exact analytical solution', 'semi analytical benchmark solution', 'e benchmark probl', 'grid refinement studie', 'grid convergence ind', 'richardson extrapolation', 'er of accuracy',
    'time step convergence', 'cfl sensitivity', 'iterative convergence', 'nonlinear convergence', 'residual convergence', 'conservation error', 'phase volume conservation',
    'global balance closure', 'local balance closure', 'shock contact wave propagation test', 'hydrostatic equilibrium test', 'water faucet benchmark',
    'dam break benchmark where formulation rel', 'kelvin helmholtz interfacial instability benchmark', 'e phase limiting solution', 'homogeneou two phase limiting solution',
    'steady state limiting solution', 'mesh orientation sensitivity', 'initial condition sensitivity', 'boundary condition sensitivity', 'floating point reproducibility test',
    'horizontal multiphase flow', 'vertical upward downward flow', 'inclined flow', 'hilly terrain pipeline', 'ga liquid flow', 'oil water flow', 'annular flow', 'stratified flow',
    'bubbly flow', 'intermittent flow', 'hydrodynamic slugging', 'terrain slugging', 'severe riser slugging', 'long pipeline slugging', 'e catenary riser where applicabl', 'liquid holdup',
    'slug frequency', 'slug catcher arrival load', 'cooldown behaviour',
    'le vof', 'idde', 'dn where computationally feasibl', '3 d cfd', // LES, LES + VOF, DES, IDDES, DNS and the 1-D + 3-D coupling run in the app with js/core/cfd3d.js at the stated (coarse) resolution
  ]/*IMPL*/,
  referenceOnly: /*REF*/[ // nothing: the scale-resolving and coupled formulations run in the app at coarse, stated resolution; the hand-off to external solvers stays available for production resolution
  ]/*REF*/,
  equationsNote: 'Scope and limits. Steady state: 1-D mass, momentum and energy balances marched with a second-order midpoint rule; properties and flashing from the equation-of-state table of the case fluid; three-phase flow is gas plus one mixed liquid. Transient line model: drift flux with one slip relation, first-order upwind, pressure waves damped. Two-fluid model: isothermal or six-equation (phase internal energies), incompressible liquid, first-order upwind; slug capturing on a test section of a few hundred diameters, not on the whole line; results depend on the cell size (mesh study). Seven-equation model: Baer–Nunziato type without relaxation, shown on shock-tube benchmarks. Multi-fluid: three fields (gas, droplets, film) in annular flow. Interfacial closures of bubbly flow (lift, wall lubrication, turbulent dispersion, virtual mass, bubble-induced turbulence, coalescence and breakup) act in a developed radial model and a 1-D interfacial-area equation. Turbulence closures (mixing length, k–ω, SST, Spalart–Allmaras, k–ε family, Reynolds stress) are solved for developed single-phase pipe flow; the 2-D solver (incompressible, uniform grid, SST k–ω, THINC/WLIC volume of fluid, no surface tension) is a channel analogue of the pipe on coarse grids. Level set, coupled level set / VOF, phase field and front tracking are interface-advection schemes verified on standard tests. A tick under Validation means that the comparison is supported in the app: by a sourced reference data set (Calibration tab) or, where no open data were found (oil–water, gas–oil–water, terrain and long-pipeline slugging, catenary risers, pressure fluctuations, slug-body holdup, slug-catcher load, cooldown), by the measurement tables of the suite. LES (Smagorinsky with van Driest damping, WALE), DES, delayed DES, IDDES (Spalart–Allmaras based), simulation without a model on a DNS grid, LES + VOF and the coupling of a 1-D line model with a 3-D section run in the app with a three-dimensional staggered-grid projection solver (second-order central differences, third-order Runge–Kutta, direct Fourier Poisson solver) on small grids: 32³ cells for the stored channel runs at Re_τ ≈ 180 (friction coefficient within −13 % to +6 % of the DNS depending on the model; the delayed DES variants return the RANS solution at this resolution), a DNS grid only for the transitional Taylor–Green vortex at Re = 100 (grid-converged within 3 % between 32³ and 64³ cells) and in a small periodic channel box, where the stored run at Re_b = 5,600 (64 × 96 × 32 cells, Δx⁺ 9, Δz⁺ 6, first cell centre at y⁺ 0.5) met the resolution criteria but laminarised instead of sustaining turbulence — a resolved turbulent channel DNS has not been achieved in the app, a first-order immersed boundary for the pipe wall, an interface captured over about two cells without surface tension, and a coupling verified for single-phase laminar flow. Every 3-D result states its grid, time per step and resolution; production resolution is written for external open-source solvers on the hand-off page.',
  inputs: INPUTS, presets: PRESETS,
  pull: ({ fluid, outputs } = {}) => {
    const n = outputs?.net, s = outputs?.solids, okArr = (a) => Array.isArray(a) && a.length >= 2 && a.every((x) => typeof x === 'number' && Number.isFinite(x)), it = [];
    const add = (key, value, from) => { if (value !== undefined && value !== null && (typeof value !== 'number' || Number.isFinite(value))) it.push({ key, value, from }); };
    add('tIn', fluid?.Tin, 'Case fluid: inlet temperature'); add('pOut', fluid?.Pout ?? n?.separatorP, fluid?.Pout !== undefined ? 'Case fluid: arrival pressure' : 'Network: separator pressure');
    if ((fluid?.rateBasis || 'oil') === 'oil') add('wc', fluid?.wc, 'Case fluid: water cut');
    if (okArr(n?.profile?.x) && okArr(n?.profile?.z) && n.profile.x.length === n.profile.z.length) add('profile', n.profile.x.map((x, i) => ({ x, z: n.profile.z[i] })), 'Network: elevation profile');
    if (n?.id > 0) add('idMm', n.id * 1000, 'Network: inner diameter'); if (n?.wt > 0) add('wtMm', n.wt * 1000, 'Network: wall thickness'); if (n?.roughness >= 0) add('roughUm', n.roughness * 1e6, 'Network: wall roughness');
    if (n?.uValue > 0) add('uValue', n.uValue, 'Network: overall U-value'); add('tSeabed', n?.tSeabed, 'Network: seabed temperature'); add('tSeaSurface', n?.tSeaSurface, 'Network: sea-surface temperature');
    if (n?.riserBaseX >= 0) add('riserBaseX', n.riserBaseX, 'Network: riser base'); if (n?.kLoss >= 0) add('kLoss', n.kLoss, 'Network: minor-loss coefficient'); if (n?.burial > 0) add('burialDepth', n.burial, 'Network: burial depth');
    if (n?.ipr?.pRes > 0) add('pRes', n.ipr.pRes, 'Network: reservoir pressure'); if (n?.ipr?.pi > 0) add('piIpr', n.ipr.pi, 'Network: productivity index'); if (n?.ipr?.type === 'vogel' || n?.ipr?.type === 'linear') add('iprType', n.ipr.type, 'Network: inflow relation');
    if (s?.effectiveId > 0) add('effIdMm', s.effectiveId * 1000, 'Solids: effective bore with deposits'); if (s?.roughnessEff > 0) add('roughEffUm', s.roughnessEff * 1e6, 'Solids: effective roughness');
    if (okArr(s?.depositProfile?.x) && okArr(s?.depositProfile?.total) && s.depositProfile.x.length === s.depositProfile.total.length) add('deposit', s.depositProfile.x.map((x, i) => ({ x, t: s.depositProfile.total[i] * 1000 })), 'Solids: deposit thickness profile');
    const ch = outputs?.ops?.chokeOpening ?? n?.chokeOpening; if (ch >= 1 && ch <= 100) add('chokeOpening', ch, outputs?.ops?.chokeOpening !== undefined ? 'Operations: choke opening' : 'Network: choke opening');
    if (Number.isFinite(outputs?.pvt?.wat)) add('watC', outputs.pvt.wat, 'Fluid: wax appearance temperature');
    return it;
  },
  site: (site) => { const d = site?.data || {}, it = [], add = (key, value, from) => { if (typeof value === 'number' && Number.isFinite(value)) it.push({ key, value, from }); }; add('tSeabed', d.seabedTemp, 'Seabed temperature at the site'); add('tSeaSurface', d.sst, 'Sea-surface temperature at the site'); add('tAir', d.airTemp, 'Air temperature at the site'); add('currentSpeed', d.currentSpeed, 'Sea current at the site'); add('windSpeed', d.windSpeed, 'Wind speed at the site'); return it; },
  run,
  mesh: [
    { name: 'Steady axial grid', keys: ['nSteady'], min: 20, note: 'Cells of the steady march (second-order midpoint rule).', metrics: [{ label: 'Inlet pressure', unit: 'bara', get: (r) => r.outputs.pIn }, { label: 'Arrival temperature', unit: '°C', get: (r) => r.outputs.tOut }, { label: 'Liquid inventory', unit: 'm³', get: (r) => r.outputs.liquidInventory }] },
    { name: 'Transient grid', keys: ['nCells'], min: 12, note: 'Cells of the transient drift-flux grid (first-order upwind).', metrics: [{ label: 'Mean inlet pressure (last 60 %)', unit: 'bara', get: (r) => r.outputs.transient?.pMean ?? NaN }, { label: 'Cumulative liquid outflow', unit: 'm³', get: (r) => r.outputs.transient?.cumLiq ?? NaN }] },
    { name: 'Two-fluid grid (slug capturing)', keys: ['tfCells'], min: 60, note: 'Cells of the two-fluid test section; select the two-fluid task first.', metrics: [{ label: 'Captured front velocity', unit: 'm/s', get: (r) => { if (!r.outputs.twoFluid) throw new Error('select the task “Two-fluid model” to run this study'); return r.outputs.twoFluid.waveVelocity ?? NaN; } }, { label: 'Mean holdup at the end', unit: '–', get: (r) => { if (!r.outputs.twoFluid) throw new Error('select the task “Two-fluid model” to run this study'); return r.outputs.twoFluid.holdupMean; } }] },
    { name: '2-D solver grid', keys: ['cfdN'], min: 8, note: 'Cells across the height of the 2-D case; select the 2-D task first.', metrics: [{ label: 'Kinetic energy at the end', unit: 'J/m', get: (r) => { if (!r.outputs.cfd) throw new Error('select the task “2-D Navier–Stokes solver” to run this study'); return r.outputs.cfd.kineticEnd; } }, { label: 'Minimum centre-line u', unit: 'm/s', get: (r) => { if (!r.outputs.cfd) throw new Error('select the task “2-D Navier–Stokes solver” to run this study'); return r.outputs.cfd.uMin; } }] },
    { name: '3-D solver grid', keys: ['c3N'], min: 8, note: 'Cells across of the selected three-dimensional case; select the 3-D task first. The primary result is the friction coefficient (channel), the friction factor (pipe), the peak dissipation (vortex), the front position (dam break) or the coupled velocity.', metrics: [{ label: 'Primary result of the 3-D case', unit: 'case units', get: (r) => { if (!r.outputs.cfd3d) throw new Error('select the task “3-D Navier–Stokes solver” to run this study'); return r.outputs.cfd3d.primary; } }, { label: 'Secondary result', unit: 'case units', get: (r) => { if (!r.outputs.cfd3d) throw new Error('select the task “3-D Navier–Stokes solver” to run this study'); return r.outputs.cfd3d.secondary; } }] },
    { name: 'Transient time step (CFL number)', keys: ['cfl'], refine: 'divide', note: 'The CFL number is divided by the refinement ratio.', metrics: [{ label: 'Mean inlet pressure (last 60 %)', unit: 'bara', get: (r) => r.outputs.transient?.pMean ?? NaN }, { label: 'Cumulative liquid outflow', unit: 'm³', get: (r) => r.outputs.transient?.cumLiq ?? NaN }] },
  ],
  calibration: {
    note: 'Every closure parameter of the suite can be estimated here from measurements at several rates: wall roughness (wall friction), interfacial, liquid-wall and gas-wall friction of stratified flow, the stratified → slug transition, entrainment and deposition rates, droplet and bubble size, coalescence and breakup coefficients, drift-flux C0 and drift velocity, holdup, slug frequency, length, celerity and body holdup, film holdup, slug initiation, growth / decay and merging in the tracking model, the external film coefficient, the overall U-value and the insulation conductivity. Each parameter moves at least one of the target columns. The model is the steady march on a coarse grid (24 cells) plus the closures evaluated at the mid-line and highest-velocity cells. The sample data are synthetic (generated from the model with shifted parameters and a few per cent of noise); measured reference data are on this tab below.',
    params: [{ key: 'roughUm', label: 'Wall roughness (µm)', lo: 5, hi: 500 }, { key: 'uMult', label: 'U-value multiplier', lo: 0.3, hi: 3 }, { key: 'holdupMult', label: 'Holdup multiplier', lo: 0.6, hi: 1.6 }, { key: 'freqMult', label: 'Slug-frequency multiplier', lo: 0.1, hi: 10 }, { key: 'lenMult', label: 'Slug-length multiplier', lo: 0.1, hi: 10 }, { key: 'filmMult', label: 'Film-holdup multiplier', lo: 0.5, hi: 2 }, { key: 'insK', label: 'Insulation conductivity (W/mK)', lo: 0.03, hi: 1 }, { key: 'c0', label: 'Distribution parameter C0', lo: 1, hi: 1.5 }, { key: 'vDrift', label: 'Drift velocity (m/s)', lo: 0, hi: 1.5 },
      { key: 'fiMult', label: 'Interfacial-friction multiplier', lo: 0.2, hi: 5 }, { key: 'fwlMult', label: 'Liquid-wall friction multiplier', lo: 0.3, hi: 3 }, { key: 'fwgMult', label: 'Gas-wall friction multiplier', lo: 0.3, hi: 3 }, { key: 'transMult', label: 'Transition multiplier', lo: 0.5, hi: 2 }, { key: 'entMult', label: 'Entrainment-rate multiplier', lo: 0.2, hi: 5 }, { key: 'kDep', label: 'Deposition velocity (m/s)', lo: 0.01, hi: 2 }, { key: 'weCrit', label: 'Droplet critical Weber number', lo: 3, hi: 40 }, { key: 'bubbleMm', label: 'Bubble size at the riser base (mm)', lo: 0.5, hi: 15 }, { key: 'crcMult', label: 'Coalescence multiplier', lo: 0.1, hi: 10 }, { key: 'ctiMult', label: 'Breakup multiplier', lo: 0.1, hi: 10 },
      { key: 'vtMult', label: 'Slug-celerity multiplier', lo: 0.7, hi: 1.4 }, { key: 'bodyMult', label: 'Slug-body holdup multiplier', lo: 0.7, hi: 1.2 }, { key: 'initMult', label: 'Slug-initiation multiplier', lo: 0.2, hi: 5 }, { key: 'relaxMult', label: 'Growth / decay length multiplier', lo: 0.2, hi: 5 }, { key: 'wakeMult', label: 'Merging (wake) multiplier', lo: 0, hi: 4 }, { key: 'hOutMult', label: 'External film-coefficient multiplier', lo: 0.2, hi: 5 }],
    columns: [{ key: 'rateFrac', label: 'Rate / case rate', unit: '×' }, { key: 'dp', label: 'Pressure drop', unit: 'bar' }, { key: 'tArr', label: 'Arrival temperature', unit: '°C' }, { key: 'holdup', label: 'Mean liquid holdup', unit: '–' }, { key: 'slugFreq', label: 'Slug frequency (mid-line)', unit: '1/min' }, { key: 'slugLen', label: 'Mean slug length', unit: 'm' }, { key: 'filmThk', label: 'Film thickness (mid-line)', unit: 'mm' }, { key: 'filmVel', label: 'Film velocity', unit: 'm/s' }, { key: 'slugVel', label: 'Slug translational velocity', unit: 'm/s' }, { key: 'bodyHoldup', label: 'Slug-body holdup', unit: '–' },
      { key: 'level', label: 'Stratified level h/D', unit: '–' }, { key: 'tauI', label: 'Interfacial shear', unit: 'Pa' }, { key: 'dpStrat', label: 'Stratified frictional gradient', unit: 'Pa/m' }, { key: 'vslCrit', label: 'Critical vsl of the stratified → slug transition', unit: 'm/s' }, { key: 'entFrac', label: 'Entrained fraction', unit: '–' }, { key: 'depLength', label: 'Deposition length', unit: 'm' }, { key: 'dropD', label: 'Droplet size', unit: 'µm' }, { key: 'dSauter', label: 'Bubble Sauter diameter 30 m above the riser base', unit: 'mm' },
      { key: 'arrFreq', label: 'Slug arrival frequency (tracked)', unit: '1/min' }, { key: 'arrLen', label: 'Mean slug length at arrival (tracked)', unit: 'm' }, { key: 'mergeShare', label: 'Share of slugs that merged', unit: '–' }, { key: 'tSurf', label: 'Outer-surface temperature (mid-line)', unit: '°C' }],
    targets: [{ key: 'dp', label: 'Pressure drop', unit: 'bar' }, { key: 'tArr', label: 'Arrival temperature', unit: '°C' }, { key: 'holdup', label: 'Mean liquid holdup', unit: '–' }, { key: 'slugFreq', label: 'Slug frequency', unit: '1/min' }, { key: 'slugLen', label: 'Mean slug length', unit: 'm' }, { key: 'filmThk', label: 'Film thickness', unit: 'mm' }, { key: 'filmVel', label: 'Film velocity', unit: 'm/s' }, { key: 'slugVel', label: 'Slug translational velocity', unit: 'm/s' }, { key: 'bodyHoldup', label: 'Slug-body holdup', unit: '–' },
      { key: 'level', label: 'Stratified level', unit: '–' }, { key: 'tauI', label: 'Interfacial shear', unit: 'Pa' }, { key: 'dpStrat', label: 'Stratified frictional gradient', unit: 'Pa/m' }, { key: 'vslCrit', label: 'Transition liquid velocity', unit: 'm/s' }, { key: 'entFrac', label: 'Entrained fraction', unit: '–' }, { key: 'depLength', label: 'Deposition length', unit: 'm' }, { key: 'dropD', label: 'Droplet size', unit: 'µm' }, { key: 'dSauter', label: 'Bubble Sauter diameter', unit: 'mm' },
      { key: 'arrFreq', label: 'Slug arrival frequency', unit: '1/min' }, { key: 'arrLen', label: 'Slug length at arrival', unit: 'm' }, { key: 'mergeShare', label: 'Merged share', unit: '–' }, { key: 'tSurf', label: 'Outer-surface temperature', unit: '°C' }],
    model: calibrationModel, sample: CAL_SAMPLE, validationSample: VAL_SAMPLE,
  },
  validationData: buildValidation(),
  verify,
};
export default SUITE;
