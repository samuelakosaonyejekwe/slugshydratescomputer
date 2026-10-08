// Suite 3 — Multiphase thermal-hydraulics and slugging.
// Steady pressure / temperature / holdup along the elevation profile with selectable holdup closures, flow-regime
// maps and interfacial-stability criteria, hydrodynamic / terrain / severe-riser slugging (unit cell, Lagrangian
// slug tracking, lumped riser cycle), a semi-implicit finite-volume drift-flux transient, radial conduction through
// the wall and coatings, and local demonstration solvers (1-D radial RANS, 1-D interface capturing).
// SI units inside; bara, °C, mm and µm at the interfaces.
import { clamp, brent, interp1, tridiag, rng, mean, std, quantile, linspace, histogram, gci, isNum } from '../core/num.js';
import { fluidModel } from '../core/thermo.js';
import { G, frictionFactor, hInside, hOutside, uValue, seaTemperature, stratifiedLevel, flowPattern, slugVelocity, slugBodyHoldup, slugFrequency, slugLength, severeSlugging, gradient, discretise } from '../core/pipe.js';
import { BASE } from '../data/basecase.js';

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
  { value: 'homogeneous', label: 'Homogeneous (no slip)' },
];
const OWN = { hagedornBrown, gray, orkiszewski, mukherjeeBrill };
const UP_ONLY = new Set(['hagedornBrown', 'gray', 'orkiszewski']);
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
    if (vm > 1e-9 && q.vsl <= 1e-9 * vm) withAcc(r, q.rhoG, vm, vm, q.P); // single-phase gas: expansion acceleration
  } else if (model === 'transient') r = transientClosure(q, mp.slip);
  else if (model === 'zuberFindlay') r = zuberFindlay(q, mp);
  else if (OWN[model] && !(UP_ONLY.has(model) && q.theta < 0)) r = OWN[model](q);
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
    const gr = holdupGradient({ vsl, vsg, rhoL: pr.rhoL, rhoG: pr.rhoG, muL: pr.muL, muG: pr.muG, sigma: pr.sigma, D, theta: th, rough: rOf ? rOf(sMid) : rough, P: P * 1e5, fModel, waterCont: pr.phaseInv }, model, mp);
    gr.loc = kCell[i] > 0 && vm > 0 ? (kCell[i] * (pr.rhoL * (vsl / vm) + pr.rhoG * (vsg / vm)) * vm * vm) / (2 * ds) : 0; gr.dpdx += gr.loc;
    const H = gr.holdup, ta = tAmbOf(sMid, zMid), U = uOf ? uOf({ s: sMid, z: zMid, D, pr, vm, holdup: H, ta }) : U0, q = U * PI * D * (T - ta) - (heatOf ? heatOf(sMid) : 0);
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
    let pIn = null, guess = o.pGuess;
    if (!isNum(guess) && N > 40) { const c = steadySolve({ ...o, grid: undefined, n: 30, hydrate: false }); if (!c.ok) return c; guess = c.pIn; marches += c.marches; }
    if (isNum(guess)) { // secant from the guess
      let p0 = Math.max(guess, 1.05), f0 = res(p0), p1 = f0 > -900 ? Math.max(p0 - f0, 1.05) : null;
      if (Math.abs(f0) < 1e-6) pIn = p0;
      for (let k = 0; pIn === null && p1 !== null && k < 12; k++) {
        const f1 = res(p1);
        if (f1 <= -900 || p1 === p0) break;
        if (Math.abs(f1) < 1e-6) { pIn = p1; break; }
        const slope = (f1 - f0) / (p1 - p0), pn = slope > 0.05 ? p1 - f1 / slope : p1 - f1;
        p0 = p1; f0 = f1; p1 = clamp(pn, Math.max(1.05, 0.5 * p1), 2 * p1 + 5);
      }
    }
    if (pIn === null) { // robust bracket and Brent
      let lo = 1.05, hi = Math.max(target, 2) + 20, fhi = res(hi), g = 0;
      if (res(lo) > 0) return { ok: false, reason: `The line gains more pressure from elevation than it loses to friction: even ${lo} bara at the inlet arrives above the ${target} bara outlet pressure.` };
      while (fhi < 0 && hi < 1400 && g++ < 40) { lo = hi; hi = hi * 1.35 + 10; fhi = res(hi); }
      if (fhi < 0) return { ok: false, reason: 'No inlet pressure below 1,400 bara can deliver this rate: the line is too small, too long or blocked.' };
      pIn = brent(res, lo, hi, 1e-7);
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
  const network = (c) => uValue({ id: c.D, wt, kWall, layers: layers.slice(1), hIn: filmInside(c.pr, c.vm, c.holdup, c.D), hOut: c.z >= 0 ? hOutside(wind, od, 'air') : hOutside(current, od, 'seawater', c.ta), burial: buriedAt(c.s) ? { depth: burialDepth + od / 2, kSoil } : null });
  const uMode = v.uMode === 'layers' ? 'layers' : 'input', uOf = uMode === 'layers' ? (c) => network(c).U * uMult : null, heat = num(v.heatTrace, 0, 0, 5000);
  const model = HOLDUP_MODELS.some((m) => m.value === v.model) ? v.model : 'beggsBrill', fModel = ['colebrook', 'haaland', 'swamee', 'churchill'].includes(v.fModel) ? v.fModel : 'colebrook';
  const mp = { c0: num(v.c0, 1.2, 0.8, 2), vDrift: num(v.vDrift, 0.35, -5, 10), wallisN: num(v.wallisN, 0, 0, 5), holdupMult: num(v.holdupMult, 1, 0.2, 3) };
  const energy = v.energy === 'cpjt' ? 'cpjt' : 'enthalpy', n = Math.round(num(v.nSteady, 150, 8, 4000)), zEnd = profile.z[profile.z.length - 1], zBase = interp1(profile.x, profile.z, riserBaseX);
  const riserHeight = zEnd - zBase, hasRiser = riserBaseX < L - 1e-6 && riserHeight > 20 * id0;
  const fittings = (Array.isArray(v.fittings) ? v.fittings : []).map((r) => ({ x: +r?.x, K: +r?.K * (isNum(+r?.open) && +r.open > 0 && +r.open < 100 ? (100 / +r.open) ** 2 : 1) })).filter((r) => Number.isFinite(r.x) && r.K > 0).map((r) => ({ s: sOfX(clamp(r.x, 0, L)), K: r.K })), kTotal = num(v.kLoss, 0, 0, 1e5);
  const cfg = { fm, profile, length: L, fittings, kTotal, id, id0, idMin, wt, rough, idOf, od, layers, thermal, uMode, uOf, U: uIn * uMult, uMult, network, tAmbOf, heat, tIn: num(v.tIn, BASE.tIn, -40, 250), pOut: num(v.pOut, BASE.pOut, 1.05, 1300), model, fModel, mp, energy, n, riserBaseX, riserBaseS: sOfX(riserBaseX), hasRiser, riserHeight, sOfX, xOfS, burialDepth, kSoil, insT, insK, kWall, current, wind, cErosion: num(v.cErosion, 100, 30, 400), buriedAt };
  cfg.base = { fm, profile, n, id, rough, idOf, kTotal, fittings, tIn: cfg.tIn, model, fModel, mp, energy, uOf, U: cfg.U, tAmbOf, heatOf: heat > 0 ? () => heat : null, cErosion: cfg.cErosion };
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
  if (bc === 'outletP') st = need(steadySolve({ ...base, mScale, pOut: cfg.pOut, pGuess: over.pGuess }));
  else if (bc === 'inletP') st = need(steadySolve({ ...base, mScale, pIn: num(v.pInSet, 95, 1.1, 1400) }));
  else {
    const coarse = { ...base, n: Math.min(base.n, 30), hydrate: false }, pInSet = num(v.pInSet, 95, 1.1, 1400), pRes = num(v.pRes, BASE.pRes, 2, 2000), pi = num(v.piIpr, BASE.pi, 1e-3, 1e6), wellDp = num(v.wellDp, 60, 0, 1500), vogel = v.iprType === 'vogel', qMax = vogel ? (pi * pRes) / 1.8 : pi * pRes;
    const pwfOf = (q) => (vogel ? pRes * ((-0.2 + Math.sqrt(Math.max(0.04 + 3.2 * (1 - Math.min(q / qMax, 1)), 0))) / 1.6) : pRes - q / pi);
    let last = null;
    const g = (m, opts = coarse) => {
      if (bc === 'bothP') { const r = steadySolve({ ...opts, mScale: m, pIn: pInSet }); return r.ok ? r.pOut - cfg.pOut : -1e3; }
      const r = steadySolve({ ...opts, mScale: m, pOut: cfg.pOut, pGuess: last }); if (!r.ok) return -1e3; last = r.pIn;
      return pwfOf(stdLiquid(cfg.fm, m)) - wellDp * (0.8 + 0.2 * m * m) - r.pIn; // available minus required flowline inlet pressure
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
    st = need(steadySolve(bc === 'bothP' ? { ...base, mScale, pIn: pInSet } : { ...base, mScale, pOut: cfg.pOut }));
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
  const fL = fanning((rhoL * vL * DL) / muL), fG = fanning((rhoG * vG * DG) / muG), fi = o.pure ? fG : Math.max(fG, 0.0142);
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
const MANDHANE = { annular: [[0.01, 70], [0.1, 60], [0.3, 38], [0.56, 40], [1, 50], [2.5, 100], [14, 230], [30, 269]], wave: [[0.01, 14], [0.1, 10.5], [0.3, 2.5], [0.5, 2.5], [1.7, 3.25], [14, 3.25]] };
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
 *      bodyModel: 'gregory' | 'barnea', freqMult, lenMult }
 * Returns { vt, C0, vd, holdupSlug, holdupFilm, holdup, slugFraction, freq, period, length, lengthFromFreq, lengthMax, unitLength,
 *           volume (m³ liquid per slug), filmThickness (m), filmVelocity, vBody (liquid velocity in the slug body), pickup (m³/s scooped at the front) }.
 */
export function slugUnitCell(p, o = {}) {
  const { vsl, vsg, rhoL, rhoG, muL, D, theta = 0 } = p, sigma = Math.max(p.sigma ?? 0.02, 1e-4), vm = vsl + vsg, A = (PI * D * D) / 4, dRho = Math.max(rhoL - rhoG, 1), sinP = Math.sin(Math.max(theta, 0));
  let { C0, vd } = slugVelocity(vm, D, theta);
  if (o.vtModel === 'nicklin') { C0 = 1.2; vd = 0.35 * Math.sqrt(G * D) * sinP; }
  const vt = Math.max(C0 * vm + vd, 1e-6);
  let HLS = slugBodyHoldup(vm);
  if (o.bodyModel === 'barnea') { const fs = fanning((rhoL * vm * D) / muL), x = 2 * Math.sqrt((0.4 * sigma) / (dRho * G)) * ((2 * fs * vm ** 3) / D) ** 0.4 * (rhoL / sigma) ** 0.6 - 0.725; HLS = clamp(1 - (x > 0 ? 0.058 * x * x : 0), 0.48, 1); }
  const vGb = 1.2 * vm + 1.53 * ((G * sigma * dRho) / (rhoL * rhoL)) ** 0.25 * Math.sqrt(HLS) * sinP;
  const holdup = clamp((vt * HLS + vGb * (1 - HLS) - vsg) / vt, vsl / Math.max(vm, 1e-9), 1);
  let HLF = clamp(stratifiedLevel({ ...p, vsl: Math.max(vsl * 0.3, 1e-4), theta: Math.min(theta, 0.15) }).holdup, 0.01, 0.9 * HLS);
  HLF = Math.min(HLF, holdup * 0.98);
  const beta = clamp((holdup - HLF) / Math.max(HLS - HLF, 1e-6), 0.02, 1), dIn = D / 0.0254;
  const length = (o.lenMult || 1) * (o.lengthModel === 'norris' ? Math.exp(-2.099 + 4.859 * Math.sqrt(Math.log(Math.max(dIn, 1.01)))) * FT : slugLength(D, vm, o.lengthModel === 'brill' ? 'brill' : 'scott'));
  const unitLength = length / beta, lam = vsl / Math.max(vm, 1e-9);
  let f = o.freqModel === 'unitCell' ? vt / unitLength : o.freqModel === 'heywood' ? 0.0434 * (lam * (2.02 / D + (vm * vm) / (G * D))) ** 1.02 : slugFrequency(vsl, vm, D, theta, o.freqModel === 'gregory' ? 'gregory' : 'zabaras');
  if (!(f > 0)) f = vt / unitLength;
  const freq = f * (o.freqMult || 1), vBody = (vm - vGb * (1 - HLS)) / HLS, filmVelocity = vt - ((vt - vBody) * HLS) / HLF;
  return { vt, C0, vd, holdupSlug: HLS, holdupFilm: HLF, holdup, slugFraction: beta, freq, period: 1 / freq, length, lengthFromFreq: (beta * vt) / freq, lengthMax: length * Math.exp(3.09 * 0.5 - 0.125), unitLength,
    volume: length * A * HLS, filmThickness: levelOfHoldup(HLF) * D, filmVelocity, vBody, pickup: A * (vt - vBody) * HLS };
}

/**
 * Lagrangian slug tracking on a steady carrier field. Slugs are initiated at the given sites with random intervals and
 * lengths (seeded). The tail moves with the bubble-nose velocity (accelerated in the wake of a short slug); the body length
 * follows the liquid balance  d(L·HLS)/dt = pick-up at the front − shedding at the tail, where the pick-up is the developed
 * value at the front position corrected for the film having relaxed towards its local terrain equilibrium over a long gap
 * (thicker film on upward slopes → growth, thinner on downward slopes → decay). Slugs that catch the one ahead merge and
 * slugs that shrink below two diameters dissipate.
 * f: { s[] (uniform), vt[], holdupSlug[], shed[] (developed shedding = pick-up rate per area, m/s), mod[] (pick-up ratio of the terrain-equilibrium film to the developed film), D, A }
 * o: { sites: [{ s, freq (1/s), length (m mean at initiation) }], nSlugs, seed, sigmaL, sigmaT, prefill, qDrain (m³/s), qSlugOut, qFilmOut (m³/s), relax (m), maxOps }
 * Returns { n, arrivals: [{ t, length, volume, velocity }], meanLength, stdLength, p50, p90, p99, maxLength, lognormal { mu, sigma }, length1000,
 *           freqArrival, period, merges, dissipated, generated, surge, surgeSingle (m³), hist, from (m), tSim (s) }.
 */
export function slugTracking(f, o = {}) {
  const { nSlugs = 200, seed = 42, sigmaL = 0.5, sigmaT = 0.5, prefill = true, maxOps = 1.5e5 } = o, D = f.D, A = f.A, S = f.s, n = S.length, Ltot = S[n - 1], ds = S[1] - S[0], R = rng(seed);
  const at = (a, s) => { const u = clamp(s / ds, 0, n - 1 - 1e-9), i = Math.floor(u); return a[i] + (a[i + 1] - a[i]) * (u - i); };
  const empty = { n: 0, arrivals: [], meanLength: 0, stdLength: 0, p50: 0, p90: 0, p99: 0, maxLength: 0, lognormal: null, length1000: 0, freqArrival: 0, period: null, merges: 0, dissipated: 0, generated: 0, surge: 0, surgeSingle: 0, hist: { centers: [], counts: [] }, from: Ltot, tSim: 0 };
  // a site cannot launch slugs closer than two body lengths apart (above that the flow is a continuous liquid column with bubbles)
  let sites = (o.sites || []).filter((q) => q.freq > 0 && q.length > 0 && q.s < Ltot).map((q) => ({ ...q, freq: Math.min(q.freq, at(f.vt, q.s) / (2 * q.length)) })).sort((a, b) => a.s - b.s);
  if (!sites.length || nSlugs < 1) return empty;
  // keep the number of slugs in the line affordable: track the downstream window that holds about 220 slug units
  const unit0 = at(f.vt, sites[0].s) / sites[0].freq, from = Math.max(sites[0].s, Ltot - 220 * unit0);
  if (from > sites[0].s) { const keep = sites.filter((q) => q.s > from); sites = [{ ...sites[0], s: from }, ...keep]; }
  const relax = o.relax || 300 * D, minL = 2 * D, lnL = (m) => R.lognormal(Math.log(Math.max(m, minL)) - 0.5 * sigmaL * sigmaL, sigmaL), lnT = (fq) => R.lognormal(Math.log(1 / fq) - 0.5 * sigmaT * sigmaT, sigmaT);
  const slugs = []; // ordered from the outlet (index 0) to the inlet
  if (prefill) { let x = Ltot - R.uniform(0, 1) * unit0; while (x > sites[0].s + unit0) { const sp = at(f.vt, x) * lnT(sites[0].freq), L = Math.min(lnL(sites[0].length), 0.9 * sp); slugs.push({ xt: x - L, L, rec: null }); x -= sp; } }
  sites.forEach((q) => { q.next = lnT(q.freq) * R.uniform(0, 1); });
  const fMax = Math.max(...sites.map((q) => q.freq)), tMax = (prefill ? 0 : (Ltot - from) / Math.max(at(f.vt, from), 0.1)) + (2.5 * nSlugs) / sites[0].freq;
  const lCap = 40 * Math.max(...sites.map((q) => q.length));
  let steps = 0, dt = clamp(Math.min(0.2 / fMax, (0.5 * sites[0].length) / Math.max(at(f.vt, Ltot), 0.1)), 0.05, 30), t = 0, ops = 0, merges = 0, dissipated = 0, generated = slugs.length;
  const arrivals = [], wake = (L, th) => (Math.abs(th) > 0.8 ? 1 + 8 * Math.exp((-1.06 * L) / D) : 1 + 0.56 * Math.exp((-0.46 * L) / D)); // Moissis–Griffith (steep) / Cook–Behnia (near horizontal)
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
  const qS = o.qSlugOut ?? 0, qF = Math.max(o.qFilmOut ?? 0, 0), qD = o.qDrain ?? 0;
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
 *      alphaL (mean liquid fraction of the feed line), chokeOpening (%), chokeDp (Pa at full opening and the mean rate), rough, fallback (film fraction), maxCycles }
 * Returns { stable, type, period (s | null), amplitude (Pa, riser-base pressure), pBaseMean, pBaseMin, pBaseMax, qLiqPeakRatio, qGasPeakRatio,
 *           t[], pBase[], pFeed[], wLout[], wGout[], level[] (liquid column / riser length), stage[], stages: { buildUp, production, blowout, fallback } (s), cycles }.
 */
export function riserSluggingCycle(o) {
  const { D, feedLength: Lp, riserHeight: hR, riserLength: Lr, wG, wL, rhoL, T, zG = 0.9, mwG = 20, pSep, rough = 4.5e-5 } = o, A = (PI * D * D) / 4;
  if (!(D > 0 && Lp > 0 && hR > 0 && Lr >= hR && wG > 0 && wL > 0 && rhoL > 0 && pSep > 0)) throw new Error('The riser-slugging cycle needs positive geometry, rates and separator pressure.');
  const sinB = Math.sin(Math.max(o.feedAngle || 0, 0)), sinG = hR / Lr, al = clamp(1 - (o.alphaL ?? 0.3), 0.05, 0.98), gasK = (zG * RGAS * T) / (mwG * 1e-3), qL = wL / rhoL, phi = clamp(o.fallback ?? 0.12, 0, 0.6), zc = clamp((o.chokeOpening ?? 100) / 100, 0.02, 1);
  const muL = o.muL || 2e-3, rel = rough / D, fH = (Re) => (Re < 2000 ? 64 / Math.max(Re, 1) : 1 / (-1.8 * log10((rel / 3.7) ** 1.11 + 6.9 / Re)) ** 2);
  const rhoG0 = pSep / gasK, u0 = (qL + wG / rhoG0) / A, lam0 = qL / (u0 * A), chokeK = (o.chokeDp ?? 1e5) / ((lam0 * rhoL + (1 - lam0) * rhoG0) * u0 * u0) / (zc * zc); // Δp_choke = chokeK ρ u²
  const Vg = (x) => al * A * (Lp - x), head = rhoL * G, tFill = (A * Lr) / qL, pTopL = pSep + chokeK * rhoL * (qL / A) ** 2;
  const ts = [], pb = [], pf = [], wl = [], wg = [], lv = [], sg = [], rec = (t, P, Pg, ql, wgo, lev, st) => { ts.push(t); pb.push(P); pf.push(Pg); wl.push(ql * rhoL); wg.push(wgo); lv.push(lev); sg.push(st); };
  // liquid front in the feed line from the hydrostatic balance (−1 when the gas pressure exceeds the liquid head: gas penetrates)
  const frontX = (ell, mg, full) => {
    const pT = full ? pTopL : pSep, f = (x) => (mg * gasK) / Vg(x) - pT - head * ((full ? Lr : ell - al * x) * sinG - x * sinB);
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
      x = xn; const z = ell - al * x, P = pSep + head * z * sinG; note(P, 0, 0); rec(t, P, P - head * x * sinB, 0, 0, z / Lr, 0);
      if (z >= Lr) { full = true; ell = Lr + al * x; break; }
    }
    // stage 2: production with a full riser
    if (full) {
      const P2 = pTopL + head * hR, need = Math.max((P2 * Vg(0)) / gasK - mg, 0), dt2 = clamp(need / wG / 160, 1e-3, tFill / 100);
      for (n = 0; n < 4000; n++) {
        mg += wG * dt2; t += dt2; c.d[1] += dt2;
        const xn = frontX(0, mg, true), qo = qL + (al * A * (x - Math.max(xn, 0))) / dt2; x = Math.max(xn, 0); note(P2, qo, 0); rec(t, P2, P2 - head * x * sinB, qo, 0, 1, 1);
        if (xn <= 0) break;
      }
    }
    // stage 3: gas penetration and blowout of the liquid column
    let Ls = Math.min(full ? Lr : Math.max(ell - al * x, 0.02 * Lr), Lr), b = 0, u = qL / A, pool = 0, t3 = 0, blown = false;
    for (n = 0; n < 60000 && t3 < 4 * tFill; n++) {
      const top = b + Ls >= Lr - 1e-9, Pg = (mg * gasK) / (Vg(0) + A * b * (1 - phi)), h = Math.min(0.5, (0.01 * Lr) / Math.max(Math.abs(u), 0.2), tFill / 200), au = Math.abs(u);
      const drive = (Pg - pSep - head * Ls * sinG) / (rhoL * Ls), damp = (fH((rhoL * au * D) / muL) * au) / (2 * D) + (top ? (chokeK * au) / Ls : 0);
      u = (u + h * drive) / (1 + h * damp); b = Math.max(b + u * h, 0); if (b + Ls >= Lr) { Ls = Math.max(Lr - b, 0); } mg += wG * h; pool += (qL / A) * h; t += h; t3 += h; c.d[2] += h;
      const P = Pg, qo = b + Ls >= Lr - 1e-9 ? A * Math.max(u, 0) : 0; note(P, qo, 0); rec(t, P, Pg, qo, 0, Ls / Lr, 2);
      if (Ls < 0.05 * Lr) { blown = true; break; }
      if (b <= 0 && u < 0) break; // the column settled back: the low point is blocked again
    }
    // stage 4: gas blowdown and liquid fallback
    ell = (blown ? phi * Lr : Ls) + pool;
    if (blown) for (n = 0; n < 20000; n++) {
      const Pg = (mg * gasK) / Vg(0), block = pSep + head * ell * sinG;
      if (Pg <= block) break;
      const rg = Pg / gasK, ug = Math.min(Math.sqrt((Pg - pSep) / (rg * ((0.02 * Lr) / (2 * D) + chokeK))), Math.sqrt(gasK)), wo = rg * A * ug, h = clamp((0.03 * mg) / Math.max(Math.abs(wo - wG), 1e-9), 0.01, tFill / 300);
      mg += (wG - wo) * h; ell += (qL / A) * h; t += h; c.d[3] += h; note(Pg, 0, wo); rec(t, Pg, Pg, 0, wo, ell / Lr, 3);
    }
    mg = Math.max(mg, ((pSep + head * ell * sinG) * Vg(0)) / gasK * 0.999);
    c.t1 = t; c.i1 = ts.length; cycles.push(c);
  }
  // no cycle when the pressure swing is small or when the liquid is produced without a surge (choked, gas and liquid leave together)
  const last = cycles[cycles.length - 1], period0 = last.t1 - last.t0, amplitude = last.pMax - last.pMin, stable = !(amplitude > Math.max(1e5, 0.05 * head * hR)) || !(period0 > 0) || last.qlMax < 1.5 * qL;
  const i0 = cycles[Math.max(cycles.length - 3, 0)].i0, every = Math.max(1, Math.ceil((ts.length - i0) / 380)), pick = (a) => a.filter((_, i) => i >= i0 && (i - i0) % every === 0), NAMES = ['build-up', 'production', 'blowout', 'fallback'], tFrom = ts[i0] ?? 0;
  return { stable, type: stable ? 'stable' : last.d[1] > 0 ? 'severe slugging 1 (riser fills completely)' : 'severe slugging 2 (gas penetrates before the riser is full)', period: stable ? null : period0, amplitude, pBaseMean: mean(pb.slice(last.i0, last.i1)), pBaseMin: last.pMin, pBaseMax: last.pMax,
    qLiqPeakRatio: last.qlMax / qL, qGasPeakRatio: last.wgMax / wG, t: pick(ts).map((v) => v - tFrom), pBase: pick(pb), pFeed: pick(pf), wLout: pick(wl), wGout: pick(wg), level: pick(lv), stage: pick(sg).map((k) => NAMES[k]),
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
  { const W = (o.mdot0 ?? mdotOf(0)) / A; for (let f = 0; f <= N; f++) { const i = Math.min(f, N - 1), pr = fm.at(P[i] / 1e5, T[i]), xg = xgOf(pr); jf[f] = W * (xg / pr.rhoG + (1 - xg) / pr.rhoL); vgF[f] = vlF[f] = jf[f]; } }
  const total = () => { let m = 0, e = 0; for (let i = 0; i < N; i++) { m += (Gc[i] + md[i]) * grid.ds[i] * A; e += (cG * Gc[i] + cD * md[i] + wallC) * T[i] * grid.ds[i] * A; } return { m, e }; };
  const pCap = (2.5 * Math.max(...o.init.P) + 50) * 1e5, start = total(), bal = { mIn: 0, mOut: 0, eIn: 0, eOut: 0, eLoss: 0, volErr: 0 };
  const ser = { t: [], pIn: [], pOutCell: [], qLiqOut: [], qGasOut: [], holdupOut: [], mOut: [], inv: [], tOut: [] }, fld = { t: [], holdup: [], P: [], T: [] };
  const nField = o.nField || 60, nSeries = o.nSeries || 400;
  let t = 0, dt = Math.min(dtMax, 1), steps = 0, relax = 1, rejected = 0, nextField = 0, last = { qL: 0, qG: 0, mOut: 0 };
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
      rhs[i] = kap[i] * P[i] - 0.5 * clamp(Vv[i], -0.03, 0.03) + r * (gG[i] * (BGf[i + 1] - BGf[i]) + gD[i] * (BDf[i + 1] - BDf[i])) + r * (cp * jst[i + 1] - (i > 0 ? cm * jst[i] : 0)) - (i === N - 1 ? r * cp * bco[N] * pOut : 0);
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
    for (let f = 0; f <= N; f++) { const fe = cG * FG[f] + cD * FD[f]; FE[f] = fe * (f === 0 ? tIn : fe >= 0 || f === N ? T[f - 1] : T[f]); }
    for (let i = 0; i < N; i++) {
      const r = dt / grid.ds[i], E = (cG * Gc[i] + cD * md[i] + wallC) * T[i] - r * (FE[i + 1] - FE[i]), hl = (4 * (o.U ? o.U[i] : 0)) / D, ta = o.tAmb ? o.tAmb[i] : 4;
      Gc[i] = Math.max(Gc[i] - r * (FG[i + 1] - FG[i]), 0); md[i] = Math.max(md[i] - r * (FD[i + 1] - FD[i]), 0);
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
      return { ...series, steps, rejected, dtMean: steps ? t / steps : 0, tEnd: t, xgMax, field: { s: grid.sc.slice(), x: grid.xc.slice(), t: fld.t, holdup: fld.holdup, P: fld.P, T: fld.T },
        mass: { initial: start.m, final: end.m, inflow: bal.mIn, outflow: bal.mOut, error: (end.m - start.m - bal.mIn + bal.mOut) / Math.max(start.m, 1e-9) },
        energy: { initial: start.e, final: end.e, inflow: bal.eIn, outflow: bal.eOut, loss: bal.eLoss, error: (end.e - start.e - bal.eIn + bal.eOut + bal.eLoss) / Math.max(Math.abs(start.e) + Math.abs(bal.eIn), 1e-9) },
        volErrMax: bal.volErr, final: { P: Array.from(P, (p) => p / 1e5), T: Array.from(T), holdup: Array.from(Hh), j: Array.from(jf) } };
    },
  };
  return api;
}
/**
 * Run the transient drift-flux solver to the end time (options as makeTransient).
 * Returns { t[], pIn[] (bara), qLiqOut[], qGasOut[] (actual m³/s), holdupOut[], mOut[] (kg/s), inv[] (m³ liquid in the line), tOut[], steps, dtMean,
 *   field: { s[], x[], t[], holdup[][], P[][], T[][] }, mass: { initial, final, inflow, outflow, error }, energy: { …, loss, error }, volErrMax, final: { P, T, holdup, j } }.
 */
export function transientDriftFlux(o) { const sim = makeTransient(o); while (!sim.done) sim.advance(200); return sim.result(); }
