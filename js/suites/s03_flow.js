// Suite 3 — Multiphase thermal-hydraulics and slugging.
// Steady pressure / temperature / holdup along the elevation profile with selectable holdup closures, flow-regime
// maps and interfacial-stability criteria, hydrodynamic / terrain / severe-riser slugging (unit cell, Lagrangian
// slug tracking, lumped riser cycle), a semi-implicit finite-volume drift-flux transient, radial conduction through
// the wall and coatings, and local demonstration solvers (1-D radial RANS, 1-D interface capturing).
// SI units inside; bara, °C, mm and µm at the interfaces.
import { clamp, brent, interp1, tridiag, rng, mean, std, quantile, linspace, histogram, gci, isNum } from '../core/num.js';
import { fluidModel } from '../core/thermo.js';
import { G, frictionFactor, hInside, hOutside, uValue, seaTemperature, stratifiedLevel, flowPattern, slugVelocity, slugBodyHoldup, slugFrequency, slugLength, severeSlugging, gradient, discretise, marchSteady } from '../core/pipe.js';
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
  const network = (c) => uValue({ id: c.D, wt, kWall, layers: layers.slice(1), hIn: filmInside(c.pr, c.vm, c.holdup, c.D), hOut: c.z >= 0 ? hOutside(wind, od, 'air') : hOutside(current, od, 'seawater', c.ta), burial: buriedAt(c.s) ? { depth: burialDepth + od / 2, kSoil } : null });
  const uMode = v.uMode === 'layers' ? 'layers' : 'input', uOf = uMode === 'layers' ? (c) => network(c).U * uMult : null, heat = num(v.heatTrace, 0, 0, 5000);
  const model = HOLDUP_MODELS.some((m) => m.value === v.model) ? v.model : 'beggsBrill', fModel = ['colebrook', 'haaland', 'swamee', 'churchill'].includes(v.fModel) ? v.fModel : 'colebrook';
  const mp = { c0: num(v.c0, 1.2, 0.8, 2), vDrift: num(v.vDrift, 0.35, -5, 10), wallisN: num(v.wallisN, 0, 0, 5), holdupMult: num(v.holdupMult, 1, 0.2, 3) };
  const energy = v.energy === 'cpjt' ? 'cpjt' : 'enthalpy', n = Math.round(num(v.nSteady, 150, 8, 4000)), zEnd = profile.z[profile.z.length - 1], zBase = interp1(profile.x, profile.z, riserBaseX);
  const riserHeight = zEnd - zBase, hasRiser = riserBaseX < L - 1e-6 && riserHeight > 20 * id0;
  const fittings = (Array.isArray(v.fittings) ? v.fittings : []).map((r) => ({ x: +r?.x, K: +r?.K * (isNum(+r?.open) && +r.open > 0 && +r.open < 100 ? (100 / +r.open) ** 2 : 1) })).filter((r) => Number.isFinite(r.x) && r.K > 0).map((r) => ({ s: sOfX(clamp(r.x, 0, L)), K: r.K })), kTotal = num(v.kLoss, 0, 0, 1e5);
  const chokeDp = num(v.chokeDp, 0, 0, 500), chokeOpening = num(v.chokeOpening, 100, 1, 100);
  const cfg = { fm, profile, length: L, fittings, kTotal, chokeDp, chokeOpening, chokeLoss: (m) => chokeDp * m * m * (100 / chokeOpening) ** 2, id, id0, idMin, wt, rough, idOf, od, layers, thermal, uMode, uOf, U: uIn * uMult, uMult, network, tAmbOf, heat, tIn: num(v.tIn, BASE.tIn, -40, 250), pSep: num(v.pOut, BASE.pOut, 1.05, 1300), pOut: num(v.pOut, BASE.pOut, 1.05, 1300) + chokeDp * num(v.rateFrac, 1, 1e-4, 20) ** 2 * (100 / chokeOpening) ** 2, model, fModel, mp, energy, n, riserBaseX, riserBaseS: sOfX(riserBaseX), hasRiser, riserHeight, sOfX, xOfS, burialDepth, kSoil, insT, insK, kWall, current, wind, cErosion: num(v.cErosion, 100, 30, 400), buriedAt };
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
 *      bodyModel: 'gregory' | 'barnea', freqMult, lenMult, filmMult (multiplier on the film holdup) }
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
  HLF = clamp(Math.min(HLF, holdup * 0.98) * (o.filmMult || 1), 0.005, Math.min(0.95 * HLS, 0.98 * holdup));
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
 * o: { sites: [{ s, freq (1/s), length (m mean at initiation) }], nSlugs, seed, sigmaL, sigmaT, prefill, qDrain (m³/s) | drainFactor (× the mean liquid arrival rate), qSlugOut, qFilmOut (m³/s), relax (m), maxOps }
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
  const lCap = 12 * Math.max(...sites.map((q) => q.length)); // liquid available to one slug
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
  const gasFr = (0.015 * Lp) / 3 / (2 * D); // gas-side friction coefficient of the feed line (the pocket empties over about a third of its length)
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
      const drive = (Pg - pSep - head * Ls * sinG) / (rhoL * Ls), damp = (fH((rhoL * au * D) / muL) * au) / (2 * D) + (top ? (chokeK * au) / Ls : 0) + (gasFr * (Pg / gasK) * au) / (rhoL * Ls); // wall friction, choke, friction of the gas feeding the bubble
      u = (u + h * drive) / (1 + h * damp); b = Math.max(b + u * h, 0); if (b + Ls >= Lr) { Ls = Math.max(Lr - b, 0); } mg += wG * h; pool += (qL / A) * h; t += h; t3 += h; c.d[2] += h;
      const P = Pg, qo = b + Ls >= Lr - 1e-9 ? A * Math.max(u, 0) : 0; note(P, Ls > 0.3 * Lr ? qo : 0, 0); rec(t, P, Pg, qo, 0, Ls / Lr, 2); // the peak rate is taken while the bulk of the column is still in the riser
      if (Ls < 0.05 * Lr) { blown = true; break; }
      if (b <= 0 && u < 0) break; // the column settled back: the low point is blocked again
    }
    // stage 4: gas blowdown and liquid fallback
    ell = (blown ? phi * Lr : Ls) + pool;
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
      rhs[i] = kap[i] * P[i] - clamp(Vv[i], -0.1, 0.1) + r * (gG[i] * (BGf[i + 1] - BGf[i]) + gD[i] * (BDf[i + 1] - BDf[i])) + r * (cp * jst[i + 1] - (i > 0 ? cm * jst[i] : 0)) - (i === N - 1 ? r * cp * bco[N] * pOut : 0);
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
        for (let j = 0; j <= n; j++) { const Ret = (kk[j] * kk[j]) / Math.max(ww[j], 1e-12), f2 = 1 - 0.22 * Math.exp(-((Ret / 6) ** 2)); gam[j] = 1 + nut[j] / 1.3; src[j] = 1.35 * 0.09 * fmu[j] * kk[j] * S[j] * S[j]; sink[j] = j ? (1.8 * f2 * ww[j]) / Math.max(kk[j], 1e-12) + (2 * Math.exp(-0.5 * y[j])) / (y[j] * y[j]) : 0; }
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

// =====================================================================================================
// 8. Suite assembly
// =====================================================================================================
const cellOf = (st, i) => ({ vsl: st.vsl[i], vsg: st.vsg[i], rhoL: st.rhoL[i], rhoG: st.rhoG[i], muL: st.muL[i], muG: st.muG[i], sigma: st.sigma[i], D: st.D[i], theta: st.theta[i] });
const slugOpts = (v) => ({ freqModel: ['zabaras', 'gregory', 'heywood', 'unitCell'].includes(v.freqModel) ? v.freqModel : 'unitCell', lengthModel: ['scott', 'brill', 'norris'].includes(v.lengthModel) ? v.lengthModel : 'scott', vtModel: v.vtModel === 'nicklin' ? 'nicklin' : 'bendiksen', bodyModel: v.bodyModel === 'barnea' ? 'barnea' : 'gregory', freqMult: num(v.freqMult, 1, 0.01, 100), lenMult: num(v.lenMult, 1, 0.01, 100), filmMult: num(v.filmMult, 1, 0.2, 3) });
const isSlug = (r) => r === 'slug' || r === 'churn';
const pwl = (ts, ys) => (t) => interp1(ts, ys, t);
const DEFAULT_SCHEDULE = [{ t: 0, rate: 1, dp: 0, choke: 100 }, { t: 0.25, rate: 1, dp: 0, choke: 100 }, { t: 0.5, rate: 0.6, dp: 0, choke: 100 }, { t: 6, rate: 0.6, dp: 0, choke: 100 }];

/** Regime (Taitel–Dukler / Barnea) at every node of a steady solution. */
function regimes(st) { return st.vsl.map((_, i) => (st.vsl[i] > 1e-9 && st.vsg[i] > 1e-9 ? flowPattern(cellOf(st, i)).pattern : st.regime[i])); }

/** Minimum stable rate: scan the rate downwards on a coarse grid and return the lowest rate fraction before severe or terrain slugging is predicted. */
function minimumStableRate(cfg, mNow) {
  const rates = [1.2, 1, 0.85, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.12].map((r) => r * Math.max(mNow, 1e-3)), rows = [];
  let guess, lowest = rates[0], hit = false;
  for (const m of rates) {
    const r = steadySolve({ ...cfg.base, n: 30, mScale: m, pOut: cfg.pOut, pGuess: guess, hydrate: false, tolP: 1e-3 });
    if (!r.ok) break;
    guess = r.pIn;
    const sc = cfg.hasRiser ? severeScreen(r, { riserBaseS: cfg.riserBaseS, pSep: cfg.pSep }) : null, ter = sc ? null : terrainSlugging(r), unstable = sc ? sc.severe : ter.accumulates;
    rows.push({ m, pIn: r.pIn, unstable, boe: sc ? sc.boe : null, pots: sc ? sc.pots : null, regime: sc ? sc.feedRegime : ter.worst ? ter.worst.downRegime : '—' });
    if (unstable && !hit) hit = true;
    if (!hit) lowest = m;
  }
  return { rate: hit ? lowest : rows.length ? rows[rows.length - 1].m : mNow, limited: hit, rows };
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
  const sim = makeTransient({ fm, grid, D, rough: cfg.rough, init, mdot0, mdotOf: (t) => mdotCase * rateOf(t), tInOf: () => cfg.tIn, pOutOf: (t) => Math.max(cfg.pSep + dpOf(t), 1.05), chokeOf, chokeK,
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
    { key: 'freqModel', label: 'Slug frequency', type: 'select', value: 'unitCell', options: [{ value: 'unitCell', label: 'Unit-cell balance with the length correlation' }, { value: 'zabaras', label: 'Zabaras (2000)' }, { value: 'gregory', label: 'Gregory & Scott (1969)' }, { value: 'heywood', label: 'Heywood & Richardson (1979)' }] },
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
    { key: 'ransModel', label: 'Turbulence model of the radial solve', type: 'select', value: 'mixing', options: [{ value: 'mixing', label: 'Mixing length (van Driest)' }, { value: 'komega', label: 'k–ω (Wilcox)' }, { value: 'kepsilon', label: 'Low-Reynolds k–ε (Chien)' }], showIf: (v) => v.localModels !== false },
    { key: 'compareModels', label: 'Compare all holdup models', type: 'bool', value: true },
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
const PRESETS = [
  { name: 'Deep-water oil tie-back (reference case)', values: { rateFrac: 1 } },
  { name: 'Low-rate turndown (riser slugging)', values: { rateFrac: 0.4, schedule: [{ t: 0, rate: 1, dp: 0, choke: 100 }, { t: 4, rate: 1, dp: 0, choke: 100 }], tEnd: 3, matchSteady: false } },
  { name: 'Gas-condensate trunk line', values: { fluidSel: 'gascond', profile: TRUNK, riserBaseX: 50000, idMm: 590, wtMm: 19.1, pOut: 70, tIn: 55, uValue: 12, insT: 3, insK: 0.3, model: 'mechanistic', mapX: 24000, nCells: 40, tEnd: 4, schedule: [{ t: 0, rate: 1, dp: 0, choke: 100 }, { t: 0.5, rate: 1, dp: 0, choke: 100 }, { t: 1, rate: 0.7, dp: 0, choke: 100 }, { t: 8, rate: 0.7, dp: 0, choke: 100 }] } },
  { name: 'Hilly-terrain onshore multiphase line', values: { profile: HILLY, riserBaseX: 12000, idMm: 203, wtMm: 8.2, rateFrac: 0.45, pOut: 15, tAir: 15, uMode: 'layers', insT: 0, burialDepth: 1.2, kSoil: 1.4, model: 'mukherjeeBrill', mapX: 5000, nCells: 40, tEnd: 1.5, schedule: [{ t: 0, rate: 1, dp: 0, choke: 100 }, { t: 0.2, rate: 1, dp: 0, choke: 100 }, { t: 0.4, rate: 1.3, dp: 0, choke: 100 }, { t: 4, rate: 1.3, dp: 0, choke: 100 }] } },
  { name: 'High-water-cut late life', values: { wc: 75, rateFrac: 0.7, model: 'mechanistic', tEnd: 1.5, schedule: [{ t: 0, rate: 1, dp: 0, choke: 100 }, { t: 0.3, rate: 1, dp: 0, choke: 100 }, { t: 0.4, rate: 1, dp: 5, choke: 100 }, { t: 4, rate: 1, dp: 5, choke: 100 }] } },
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
  const reg = regimes(st), units = st.vsl.map((_, i) => (st.vsl[i] > 1e-9 && st.vsg[i] > 1e-9 ? slugUnitCell(cellOf(st, i), so) : null));
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
    track = slugTracking(field, { sites, nSlugs: Math.round(num(v.nSlugs, 150, 5, 5000)), seed: Math.round(num(v.slugSeed, 42, 1, 2 ** 31)), theta: st.theta, drainFactor: num(v.drainFactor, 1.2, 1, 10), qSlugOut: A0 * uo.holdupSlug * st.vm[N], qFilmOut: A0 * uo.holdupFilm * Math.max(uo.filmVelocity, 0) });
    if (!track.n) track = null;
  }
  const terrain = terrainSlugging(st), sev = cfg.hasRiser ? severeScreen(st, { riserBaseS: cfg.riserBaseS, feedLength: num(v.feedLength, 0, 0, 1e7), pSep: cfg.pSep }) : null;
  let cyc = null;
  if (sev && sev.severe) { const k = sev.iBase - 1; try { cyc = riserSluggingCycle({ D: st.D[k], feedLength: sev.feedLength, feedAngle: sev.feedAngle * DEG, riserHeight: sev.riserHeight, riserLength: sev.riserLength, wG: Math.max(st.mG[k], 1e-6), wL: Math.max(st.mL[k], 1e-6), rhoL: st.rhoL[k], muL: st.muL[k], T: st.T[k] + 273.15, zG: st.zG[k], mwG: st.mwG[k], pSep: cfg.pSep * 1e5, alphaL: 1 - sev.alpha, chokeOpening: cfg.chokeOpening, chokeDp: Math.max(cfg.chokeDp, 0.05) * 1e5 * st.mScale ** 2, rough: cfg.rough }); } catch (e) { warnings.push({ level: 'info', msg: 'Riser-slugging cycle model: ' + e.message }); } }
  const severe = !!(sev && sev.severe && (!cyc || !cyc.stable));
  await tick();
  prog(0.4, 'Minimum stable rate');
  const turndown = minimumStableRate(cfg, st.mScale);
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
  let rans = null, ransRef = null, vof = null, vofUp = null, dam = null;
  if (v.localModels !== false) {
    const lam = cm.vsl / Math.max(cm.vsl + cm.vsg, 1e-9), rhoN = cm.rhoL * lam + cm.rhoG * (1 - lam), muN = cm.muL * lam + cm.muG * (1 - lam), ReLoc = clamp((rhoN * (cm.vsl + cm.vsg) * cm.D) / muN, 3000, 3e6), nR = Math.round(num(v.nRans, 70, 30, 400));
    const mdl = ['komega', 'kepsilon', 'mixing'].includes(v.ransModel) ? v.ransModel : 'mixing';
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
  const rateNow = st.mScale, tdFrac = turndown.rate;
  if (turndown.limited) recs.push(`Keep the rate above about ${(100 * tdFrac).toFixed(0)} % of the case rate (${(tdFrac * stdLiquid(fm, 1)).toFixed(0)} Sm³/d liquid): below it the screening predicts ${cfg.hasRiser ? 'riser-induced' : 'terrain'} slugging${rateNow < tdFrac ? ` — the present rate (${(100 * rateNow).toFixed(0)} %) is already below it` : ''}.`);
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
    plots.push({ type: 'line', title: `Developed turbulent velocity profile (radial RANS, Re = ${sig(rans.Re, 3)})`, xlabel: 'y⁺', ylabel: 'u⁺', logx: true, series: [{ name: `${v.ransModel === 'komega' ? 'k–ω' : v.ransModel === 'kepsilon' ? 'k–ε (Chien)' : 'Mixing length'}`, x: rans.y.slice(1), y: rans.u.slice(1) }, ...(ransRef ? [{ name: 'Mixing length', x: ransRef.y.slice(1), y: ransRef.u.slice(1), dash: true }] : []), { name: 'Log law (κ = 0.41, B = 5.2)', x: rans.y.filter((y) => y > 30), y: rans.y.filter((y) => y > 30).map((y) => Math.log(y) / 0.41 + 5.2), dash: true }], note: `Darcy friction factor ${rans.f.toFixed(5)} against ${frictionFactor(rans.Re, 0).toFixed(5)} from the smooth-pipe Colebrook law (${(100 * (rans.f / frictionFactor(rans.Re, 0) - 1)).toFixed(1)} %). Single-phase flow at the local no-slip mixture Reynolds number.` });
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
  tables.push({ title: 'Minimum stable rate scan', columns: ['Rate (× case)', 'Inlet pressure (bara)', 'Bøe ratio', 'Pots number', 'Feed / downhill regime', 'Slugging predicted'], rows: turndown.rows.map((r) => [r3(r.m, 2), r3(r.pIn, 1), r.boe === null ? '—' : r3(r.boe, 2), r.pots === null ? '—' : r3(r.pots, 3), r.regime, r.unstable ? 'yes' : 'no']) });
  tables.push({ title: `Heat-transfer resistances at ${Math.round(st.x[ic])} m`, columns: ['Layer', 'Resistance (m²K/W on ID)', 'Share (%)', 'Temperature on its outer side (°C)'], rows: (() => { const tot = netC.resistances.reduce((s, x) => s + x.R, 0); let acc = 0; return netC.resistances.map((x) => { acc += x.R; return [x.name, sig(x.R, 4), r3((100 * x.R) / tot, 1), r3(st.T[ic] - ((st.T[ic] - st.tAmb[kc]) * acc) / tot, 2)]; }); })(), note: `Layer-based U-value ${netC.U.toFixed(2)} W/m²K; U-value used ${st.U[kc].toFixed(2)} W/m²K; inside film ${st.hIn[kc].toFixed(0)} W/m²K.` });
  if (tr) tables.push({ title: 'Transient summary', columns: ['Quantity', 'Value', 'Unit'], rows: [['Simulated time', r3(tr.tEnd / 3600, 3), 'h'], ['Time steps (rejected)', `${tr.steps} (${tr.rejected})`, '–'], ['Mean time step', r3(tr.dtMean, 2), 's'], ['Cells', trInfo.grid.n, '–'], ['Initial condition', trInfo.uniform ? 'settled line' : 'steady state of the drift-flux closure', ''], ['Wall-friction tuning factor', r3(trInfo.slip.fricMult, 3), '×'], ['Inlet pressure: mean / min / max (last 60 %)', trSummary ? `${trSummary.pMean.toFixed(1)} / ${trSummary.pMin.toFixed(1)} / ${trSummary.pMax.toFixed(1)}` : '—', 'bara'], ['Liquid outflow: mean / peak (last 60 %)', trSummary ? `${(3600 * trSummary.qMean).toFixed(1)} / ${(3600 * trSummary.qPeak).toFixed(1)}` : '—', 'm³/h'], ['Liquid inventory: start → end', trSummary ? `${trSummary.inv0.toFixed(0)} → ${trSummary.invEnd.toFixed(0)}` : '—', 'm³'], ['Mass-conservation error', sig(Math.abs(tr.mass.error), 3), '–'], ['Energy-conservation error', sig(Math.abs(tr.energy.error), 3), '–'], ['Largest volume-constraint error', sig(tr.volErrMax, 3), '–']], note: 'Drift-flux closure v_g = C0 j + v_d (C0 = 1.2, Bendiksen drift) for every regime; convective acceleration, Joule–Thomson cooling and latent heat are not included in the transient.' });

  // ---- balances ---------------------------------------------------------------------------------------------
  const mOutS = st.rhoG[N] * st.vsg[N] * (PI * st.D[N] ** 2) / 4 + st.rhoL[N] * st.vsl[N] * (PI * st.D[N] ** 2) / 4;
  balances.push({ name: 'Steady mass flow (kg/s): inlet vs outlet phases', in: st.mdot, out: mOutS });
  balances.push({ name: 'Steady momentum (bar): inlet − outlet pressure vs friction + gravity + acceleration + local', in: dpTotal, out: st.dpFric + st.dpGrav + st.dpAcc + st.dpLocal });
  if (cfg.energy === 'enthalpy') balances.push({ name: 'Steady energy (kW): enthalpy in vs enthalpy out + heat loss + elevation + kinetic', in: st.energy.hIn / 1e3, out: (st.energy.hOut + st.heatLoss + st.energy.potential + st.energy.kinetic) / 1e3 });
  if (tr) { balances.push({ name: 'Transient mass (kg): initial + inflow vs final + outflow', in: tr.mass.initial + tr.mass.inflow, out: tr.mass.final + tr.mass.outflow }); balances.push({ name: 'Transient thermal energy (MJ): initial + inflow vs final + outflow + loss', in: (tr.energy.initial + tr.energy.inflow) / 1e6, out: (tr.energy.final + tr.energy.outflow + tr.energy.loss) / 1e6 }); }

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
    { label: 'Minimum stable rate', value: r3(100 * tdFrac, 0), unit: '% of case', status: rateNow < tdFrac ? 'bad' : 'ok' },
    { label: 'Bøe ratio / Pots number', value: sev ? `${sev.boe.toFixed(2)} / ${sev.pots === null ? '—' : sev.pots.toFixed(2)}` : 'no riser', unit: '', status: severe ? 'bad' : 'ok' },
    { label: 'Transient inlet-pressure swing', value: r3(pAmp, 1), unit: 'bar', status: pAmp > 0.1 * dpTotal && pAmp > 2 ? 'warn' : 'ok' },
    { label: 'Cooldown to hydrate temperature', value: cool.tReach !== null ? r3(cool.tReach / 3600, 1) : subMax > 0 ? 0 : `> ${(cool.t[cool.t.length - 1] / 3600).toFixed(0)}`, unit: 'h', status: cool.tReach !== null && cool.tReach < 8 * 3600 ? 'warn' : 'ok', help: `At ${Math.round(st.x[ic])} m, the location with the smallest margin` },
    { label: 'Operating rate', value: r3(stdLiquid(fm, st.mScale), 0), unit: 'Sm³/d liquid', status: 'ok', help: `GOR ${gor.toFixed(0)} Sm³/Sm³, water cut ${wcIn.toFixed(0)} %, ${st.mdot.toFixed(1)} kg/s` },
  ];

  // ---- outputs ----------------------------------------------------------------------------------------------
  const idx = [...new Set(linspace(0, N, Math.min(N + 1, 120)).map(Math.round))], pk = (a, k = 5) => idx.map((i) => sig(a[i], k + 2));
  const series = tr ? { t: tr.t.map((x) => r3(x, 1)), pIn: D0(tr.pIn, 3), qLiqOut: tr.qLiqOut.map((x) => sig(x, 5)), qGasOut: tr.qGasOut.map((x) => sig(x, 5)), holdupOut: D0(tr.holdupOut, 4) } : { t: [], pIn: [], qLiqOut: [], qGasOut: [], holdupOut: [] };
  const outputs = {
    profile: { x: idx.map((i) => r3(st.x[i], 1)), z: idx.map((i) => r3(st.z[i], 2)), P: pk(st.P), T: pk(st.T), holdup: pk(st.holdup), vsl: pk(st.vsl), vsg: pk(st.vsg), vm: pk(st.vm), rhoM: pk(st.rhoM), dpdx: pk(st.dpdx), tauW: pk(st.tauW), tAmb: pk(st.tAmb), tHyd: pk(st.tHyd), subcooling: pk(st.subcooling), regime: idx.map((i) => reg[i]), tWall: pk(st.tWall), heatFlux: idx.map((i) => sig(heatFlux[i], 6)) },
    pIn: st.pIn, pOut: st.pOut, tIn: st.tIn, tOut: st.tOut, dpTotal, dpFric: st.dpFric, dpGrav: st.dpGrav, dpAcc: st.dpAcc, dpLocal: st.dpLocal, mdot: st.mdot, qLiq: qLout, qGas: qGout, liquidInventory: st.liquidInventory, volume: st.volume, residence: st.residence / 3600, residenceLiquid: st.residenceLiquid / 3600,
    heatLoss: st.heatLoss / 1e3, uValue: mean(st.U), slug, severeSlugging: severe, boe: sev ? sev.boe : null, pots: sev ? sev.pots : null, erosionalRatio: evrMax, maxVelocity: vMax, maxSubcooling: subMax, hydrateLength: hydLen,
    series, pInAmplitude: pAmp, turndownRate: tdFrac, rateFraction: st.mScale, model: cfg.model, slugOnsetX: iOnset >= 0 ? st.x[iOnset] : null, slugCatcherVolume: catcher, terrain: { accumulates: terrain.accumulates, x: terrain.worst ? terrain.worst.x : null, volume: terrain.worst ? terrain.worst.volume : null, vsgCrit: terrain.worst ? terrain.worst.vsgCrit : null },
    cycle: cyc ? { stable: cyc.stable, period: cyc.period, amplitude: cyc.amplitude / 1e5, qLiqPeakRatio: cyc.qLiqPeakRatio } : null, arrival: track ? { n: track.n, meanLength: track.meanLength, p99: track.p99, maxLength: track.maxLength, frequency: track.freqArrival, surge: track.surge } : null,
    cooldownPreview: cool.tReach !== null ? cool.tReach / 3600 : null, coldSpotX: st.x[ic], regimeAtMap: reg[im], mandhane: mand, id: cfg.idMin, transient: tr ? { massError: tr.mass.error, energyError: tr.energy.error, steps: tr.steps, aborted: tr.aborted, pMean: trSummary ? trSummary.pMean : null, cumLiq: trSummary ? trSummary.cumLiq : null } : null,
  };
  prog(1, 'Done');
  const summary = `Inlet pressure ${st.pIn.toFixed(1)} bara for ${st.pOut.toFixed(1)} bara at the outlet (${dpTotal.toFixed(1)} bar: friction ${st.dpFric.toFixed(1)}, gravity ${st.dpGrav.toFixed(1)}); arrival at ${st.tOut.toFixed(1)} °C, ${subMax > 0 ? `${subMax.toFixed(1)} °C inside the hydrate region at worst` : `${(-subMax).toFixed(1)} °C above the hydrate temperature at the closest point`}; slugging: ${slug.type}${slug.type !== 'none' ? `, liquid surge about ${fin(slug.surge, 0).toFixed(0)} m³` : ''}; minimum stable rate ${(100 * tdFrac).toFixed(0)} % of the case rate.`;
  return { summary, kpis, warnings, recommendations: recs, plots, tables, balances, outputs };
}

// ---- calibration model: coarse steady march, fast enough for least-squares fitting ----------------------------
function calibrationModel(v) {
  const cfg = flowConfig({ ...v, nSteady: 24 }, lastCtx), r = steadySolve({ ...cfg.base, n: 24, mScale: num(v.rateFrac, 1, 1e-3, 20), pOut: cfg.pOut, hydrate: false, tolP: 1e-4 });
  if (!r.ok) return { dp: NaN, tArr: NaN, holdup: NaN, slugFreq: NaN, slugLen: NaN, filmThk: NaN };
  const i = 12, u = r.vsl[i] > 1e-9 && r.vsg[i] > 1e-9 ? slugUnitCell(cellOf(r, i), slugOpts(v)) : null;
  return { dp: r.pIn - r.pOut, tArr: r.tOut, holdup: r.liquidInventory / r.volume, slugFreq: u ? 60 * u.freq : 0, slugLen: u ? u.lengthFromFreq : 0, filmThk: u ? 1000 * u.filmThickness : 0 };
}
const CAL_SAMPLE = /*CAL*/[{rateFrac:0.5, dp:64.5, tArr:19.41, holdup:0.5397, slugFreq:0.2481, slugLen:39.2, filmThk:189.2},
  {rateFrac:0.6, dp:65.53, tArr:23.71, holdup:0.5469, slugFreq:0.2361, slugLen:45, filmThk:183.9},
  {rateFrac:0.7, dp:64.97, tArr:27.44, holdup:0.5243, slugFreq:0.2378, slugLen:39.6, filmThk:178.3},
  {rateFrac:0.8, dp:66.25, tArr:30.32, holdup:0.5188, slugFreq:0.2452, slugLen:40.4, filmThk:167.7},
  {rateFrac:0.9, dp:71.13, tArr:33.98, holdup:0.5145, slugFreq:0.299, slugLen:37.2, filmThk:177.6},
  {rateFrac:1, dp:70.95, tArr:35.36, holdup:0.5187, slugFreq:0.2942, slugLen:38.6, filmThk:171.2},
  {rateFrac:1.1, dp:74.75, tArr:37.25, holdup:0.5362, slugFreq:0.3136, slugLen:41.8, filmThk:175.6},
  {rateFrac:1.2, dp:78.33, tArr:39.31, holdup:0.511, slugFreq:0.33, slugLen:41.5, filmThk:163.6},
  {rateFrac:1.3, dp:84.7, tArr:40.54, holdup:0.5394, slugFreq:0.3793, slugLen:38.4, filmThk:166.8},
  {rateFrac:1.4, dp:84.2, tArr:41.83, holdup:0.5497, slugFreq:0.4491, slugLen:38.9, filmThk:168.7}]/*CAL*/;
const VAL_SAMPLE = /*VAL*/[{rateFrac:0.55, dp:64.92, tArr:22.51, holdup:0.539, slugFreq:0.2281, slugLen:39.6, filmThk:179.3},
  {rateFrac:0.75, dp:65.67, tArr:29.42, holdup:0.53, slugFreq:0.2464, slugLen:38.8, filmThk:165.4},
  {rateFrac:0.95, dp:70.73, tArr:35.29, holdup:0.5178, slugFreq:0.2908, slugLen:40.4, filmThk:170.1},
  {rateFrac:1.05, dp:75.01, tArr:36.48, holdup:0.5145, slugFreq:0.3092, slugLen:41.8, filmThk:168.7},
  {rateFrac:1.25, dp:80.17, tArr:39.22, holdup:0.5405, slugFreq:0.3703, slugLen:43.5, filmThk:170.4},
  {rateFrac:1.5, dp:90.86, tArr:42.79, holdup:0.5375, slugFreq:0.6039, slugLen:39.7, filmThk:163.7}]/*VAL*/;

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
  return out;
}
function marchKernel(cfg) { return marchSteady({ fm: cfg.fm, profile: cfg.profile, id: cfg.id, rough: cfg.rough, U: cfg.U, tAmbOf: cfg.tAmbOf, tIn: cfg.tIn, pOut: cfg.pOut, n: cfg.n }); }

const SUITE = {
  id: 'flow', num: 3, title: 'Multiphase Thermal-Hydraulics & Slugging', short: 'Flow · Slugs', icon: '🌊',
  tagline: 'Pressure, temperature, holdup and flow regime along the line, slugging of every kind, and the transient response to rate and pressure changes.',
  description: 'Solves the steady mass, momentum and energy balances along the elevation profile with nine selectable holdup closures and the equation-of-state property table, classifies the flow regime mechanistically, and evaluates hydrodynamic, terrain and severe riser slugging with a unit-cell model, Lagrangian slug tracking and a lumped riser cycle. A semi-implicit finite-volume drift-flux model gives the transient response to rate, back-pressure and choke histories; radial conduction, a radial RANS solve and 1-D interface capturing resolve local detail.',
  guide: [
    'Check the boundary conditions: by default the case rates and the outlet pressure are fixed and the inlet pressure is solved. Pull the profile, diameter and U-value from the Network suite when it has been run.',
    'Choose the holdup model on the Model tab; the model-comparison table shows how much the inlet pressure depends on that choice.',
    'Read the slugging results: type, frequency, largest slug, surge volume and the minimum stable rate. For a riser, the Bøe ratio, the Pots number and the cycle model are reported.',
    'Edit the boundary history for the transient (rate ramp, back-pressure step, choke move) and the simulated time; the distance–time maps show liquid accumulation and surges.',
    'Use the Mesh tab to quantify the numerical uncertainty of the steady march and of the transient grid and time step.',
    'Fit roughness, U-value, holdup and slug multipliers to measured pressure drop, arrival temperature, holdup and slug frequency on the Calibration tab.',
  ],
  implemented: /*IMPL*/[
    'conservation of mass', 'component mass conservation', 'phase continuity equations', 'mixture continuity equation', 'conservation of linear momentum', 'mixture mo',
    'conservation of total energy', 'mixture-energy equation', 'enthalpy equation', 'drift-flux', 'homogeneous-equilibrium', 'separated-flow', 'mechanistic multiphase', 'nicklin',
    'beggs–brill', 'hagedorn–brown', 'mukherjee–brill', 'gray correlation', 'orkiszewski', 'taitel', 'mandhane', 'barnea unified', 'mechanistic stratified/annular/intermittent',
    'unit-cell slug', 'mechanistic slug-flow', 'slug-tracking', 'slug-frequency', 'translational-velocity', 'slug-body/film mass', 'slug growth', 'slug merging',
    'kelvin–helmholtz instability criterion', 'viscous kelvin–helmholtz', 'long-wave stability', 'linear stability', 'roll-wave', 'pipeline–riser liquid accumulation', 'severe-slugging cyc',
    'hydrodynamic stability criteria', 'pressure-build-up/liquid-fallback', 'terrain-induced accumulation', 'interfacial drag', 'wall friction', "fourier's law", 'transient heat-conduction',
    "newton's law of cooling", 'overall heat-transfer resistance', 'conjugate heat-transfer', 'radial multilayer conduction', 'seabed conduction', 'thomson relation', 'volume of fluid',
    'mechanistic flow-regime se', 'eos + thermal-hydraulic', 'pressure', 'temperature', 'initial oil', 'gas and water flow rates', 'phase velocities and superficial velocities',
    'liquid holdup', 'flow regime', 'initial liquid-film', 'fluid inventory', 'initial phase distribution', 'where slugging already exists', 'initial number', 'position', 'ngth', 'velocity',
    'liquid content and distribution of slugs', 'appropriate combination of oil', 'gas and water mass or volumetric flow rates', 'phase fractions', 'water cut',
    'external heat-transfer conditions', 'insulation', 'valve', 'choke', 'production-rate changes', 'terrain or severe', 'validated outputs from mo', 'oil/gas/water mass', 'gor',
    'well inflow conditions', 'transient operating boundary histories', 'numerical discretization and convergence controls', 'phase and superficial velocities', 'slug initiation/onset',
    'terrain and severe', 'liquid surge', 'slug arrival statistics', 'film thickness/velocity', 'slip/drift velocity', 'distribution parameter', 'heat-transfer coefficients',
    'overall u-value', 'method of manufactured solutions', 'exact analytical solutions', 'semi-analytical benchmark', 'benchmark prob', 'grid-refinement', 'grid convergence in',
    'richardson extrapolation', 'observed or', 'formal or', 'time-step convergence', 'cfl sensitivity', 'iterative convergence', 'nonlinear convergence', 'residual convergence',
    'mass-conservation error', 'momentum-conservation error', 'energy-conservation error', 'phase-volume conservation', 'global balance closure', 'local balance closure',
    'shock/contact-wave propagation', 'hydrostatic-equilibrium test', 'dam-break', 'phase limiting solution', 'steady-state limiting', 'initial-condition sensitivity',
    'boundary-condition sensitivity', 'floating-point/reproducibility',
  ]/*IMPL*/,
  referenceOnly: /*REF*/[
    'phase momentum equations', 'phase-energy equations', 'internal-energy equation', 'entropy inequality', 'species-transport', 'two-fluid', 'multi-fluid', 'six-equation', 'seven-equation',
    'homogeneous-relaxation', 'duns–ros', 'ansari', 'baker flow map', 'slug-capturing', 'nonlinear wave-growth', 'lift force', 'virtual/added mass', 'wall lubrication',
    'turbulent dispersion', 'bubble-induced turbulence', 'entrainment', 'deposition rate', 'coalescence/breakup', 'natural/forced convection', 'rng k–ε', 'realizable k–ε', 'sst k–ω',
    'reynolds-stress', 'spalart–allmaras', 'iddes', 'dns where computationally feasible', 'level set', 'front tracking', 'phase field', '1-d transient flow + 3-d cfd', 'rans + vof',
    'les + vof', 'riser and well network', 'appropriate flows between connected branches', 'pump', 'compressor and separator boundaries', 'interfacial-friction correlations',
    'liquid-wall shear closure', 'gas-wall shear closure', 'interfacial shear', 'bubble size', 'droplet size', 'flow-regime transition parameters', 'slug-celerity parameters',
    'slug-body holdup', 'slug initiation parameters', 'slug growth/decay parameters', 'slug merging parameters', 'external heat-transfer coefficient', 'water-faucet',
    'kelvin–helmholtz/interfacial-instability benchmark', 'mesh-orientation sensitivity', 'horizontal multiphase flow', 'vertical upward/downward flow', 'inclined flow',
    'hilly-terrain pipelines', 'gas-liquid flow', 'oil-water flow', 'gas-oil-water flow', 'annular flow', 'stratified flow', 'bubbly flow', 'intermittent flow', 'hydrodynamic slugging',
    'terrain slugging', 'long-pipeline slugging', 'flexible/catenary risers', 'pressure fluctuations', 'slug velocity', 'slug-catcher arrival/load', 'cooldown behaviour',
  ]/*REF*/,
  equationsNote: 'Scope and limits. Steady state: 1-D mass, momentum and energy balances marched with a second-order midpoint rule; properties and flashing from the equation-of-state table of the case fluid; three-phase flow is treated as gas plus one mixed liquid. Transient: drift-flux model with one slip relation for all regimes, first-order upwind, convective acceleration, Joule–Thomson cooling and latent heat neglected; pressure waves are damped rather than resolved and the run stops if the line becomes liquid-full. Slug tracking is one-way coupled to the steady solution; the riser cycle is a lumped quasi-equilibrium model. The RANS solve is single-phase, fully developed and smooth-walled; the interface-capturing solvers are 1-D demonstrations. Two-fluid, slug-capturing, 3-D CFD, LES/DES/DNS and Reynolds-stress models are not solved here.',
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
    { name: 'Transient time step (CFL number)', keys: ['cfl'], refine: 'divide', note: 'The CFL number is divided by the refinement ratio.', metrics: [{ label: 'Mean inlet pressure (last 60 %)', unit: 'bara', get: (r) => r.outputs.transient?.pMean ?? NaN }, { label: 'Cumulative liquid outflow', unit: 'm³', get: (r) => r.outputs.transient?.cumLiq ?? NaN }] },
  ],
  calibration: {
    note: 'Fit wall roughness, the U-value multiplier, the holdup multiplier (or C0 and drift velocity when the Zuber–Findlay model is selected) and the slug-frequency, slug-length and film-holdup multipliers and the insulation conductivity to measured pressure drop, arrival temperature, mean holdup, slug frequency, slug length and film thickness at several rates. The model is the steady march on a coarse grid (24 cells). The sample data are synthetic.',
    params: [{ key: 'roughUm', label: 'Wall roughness (µm)', lo: 5, hi: 500 }, { key: 'uMult', label: 'U-value multiplier', lo: 0.3, hi: 3 }, { key: 'holdupMult', label: 'Holdup multiplier', lo: 0.6, hi: 1.6 }, { key: 'freqMult', label: 'Slug-frequency multiplier', lo: 0.1, hi: 10 }, { key: 'lenMult', label: 'Slug-length multiplier', lo: 0.1, hi: 10 }, { key: 'filmMult', label: 'Film-holdup multiplier', lo: 0.5, hi: 2 }, { key: 'insK', label: 'Insulation conductivity (W/mK)', lo: 0.03, hi: 1 }, { key: 'c0', label: 'Distribution parameter C0', lo: 1, hi: 1.5 }, { key: 'vDrift', label: 'Drift velocity (m/s)', lo: 0, hi: 1.5 }],
    columns: [{ key: 'rateFrac', label: 'Rate / case rate', unit: '×' }, { key: 'dp', label: 'Pressure drop', unit: 'bar' }, { key: 'tArr', label: 'Arrival temperature', unit: '°C' }, { key: 'holdup', label: 'Mean liquid holdup', unit: '–' }, { key: 'slugFreq', label: 'Slug frequency (mid-line)', unit: '1/min' }, { key: 'slugLen', label: 'Mean slug length', unit: 'm' }, { key: 'filmThk', label: 'Film thickness (mid-line)', unit: 'mm' }],
    targets: [{ key: 'dp', label: 'Pressure drop', unit: 'bar' }, { key: 'tArr', label: 'Arrival temperature', unit: '°C' }, { key: 'holdup', label: 'Mean liquid holdup', unit: '–' }, { key: 'slugFreq', label: 'Slug frequency', unit: '1/min' }, { key: 'slugLen', label: 'Mean slug length', unit: 'm' }, { key: 'filmThk', label: 'Film thickness', unit: 'mm' }],
    model: calibrationModel, sample: CAL_SAMPLE, validationSample: VAL_SAMPLE,
  },
  verify,
};
export default SUITE;
