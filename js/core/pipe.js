// Pipe-flow kernel shared by the suites: friction factors, heat-transfer resistances, two-phase flow-pattern
// prediction (Taitel–Dukler / Taitel–Barnea–Dukler), liquid holdup and pressure-gradient closures (Beggs–Brill,
// drift flux, mechanistic stratified / slug unit cell), slug closures, severe-slugging criteria, sea-temperature
// profile and a steady-state marching solver for pressure and temperature along an elevation profile.
// SI units throughout except pressure in bara and temperature in °C at the interface of marchSteady.
import { clamp, brent, interp1 } from './num.js';

export const G = 9.80665;

// ---- single-phase friction (Darcy friction factor) -----------------------------------------------------
const fTurb = {
  haaland: (Re, e) => 1 / (-1.8 * Math.log10((e / 3.7) ** 1.11 + 6.9 / Re)) ** 2,
  swamee: (Re, e) => 0.25 / Math.log10(e / 3.7 + 5.74 / Re ** 0.9) ** 2,
  colebrook: (Re, e) => { let x = -1.8 * Math.log10((e / 3.7) ** 1.11 + 6.9 / Re); for (let i = 0; i < 30; i++) { const xn = -2 * Math.log10(e / 3.7 + (2.51 * x) / Re); if (Math.abs(xn - x) < 1e-12) { x = xn; break; } x = xn; } return 1 / (x * x); },
};
/** Darcy friction factor. model: 'colebrook' | 'haaland' | 'swamee' | 'churchill'. rel = roughness / diameter. */
export function frictionFactor(Re, rel = 0, model = 'colebrook') {
  Re = Math.max(Re, 1e-9);
  if (model === 'churchill') { const A = (2.457 * Math.log(1 / ((7 / Re) ** 0.9 + 0.27 * rel))) ** 16, B = (37530 / Re) ** 16; return 8 * ((8 / Re) ** 12 + (A + B) ** -1.5) ** (1 / 12); }
  const lam = 64 / Re;
  if (Re <= 2000) return lam;
  const turb = (fTurb[model] || fTurb.colebrook)(Math.max(Re, 4000), rel);
  if (Re >= 4000) return turb;
  const w = (Re - 2000) / 2000; // smooth bridge across the critical zone
  return (1 - w) * (64 / 2000) + w * turb;
}

// ---- heat transfer -------------------------------------------------------------------------------------
/** Inside film coefficient (W/m²/K): Nu = 3.66 laminar, Gnielinski turbulent. */
export function hInside(Re, Pr, k, D) {
  if (Re < 2300) return (3.66 * k) / D;
  const f = frictionFactor(Re, 0, 'haaland'), Nu = ((f / 8) * (Re - 1000) * Pr) / (1 + 12.7 * Math.sqrt(f / 8) * (Pr ** (2 / 3) - 1));
  return (Math.max(Nu, 3.66) * k) / D;
}
/** Outside film coefficient for cross-flow over a cylinder (Churchill–Bernstein). medium: 'seawater' | 'air'. */
export function hOutside(v, Do, medium = 'seawater', T = 4) {
  const p = medium === 'air' ? { rho: 1.2, mu: 1.8e-5, k: 0.026, Pr: 0.71 } : { rho: 1027, mu: 1.9e-3 * Math.exp(-0.027 * T) + 3.5e-4, k: 0.57, Pr: 13.4 * Math.exp(-0.027 * T) + 2 };
  const Re = Math.max((p.rho * Math.max(v, 0.01) * Do) / p.mu, 1), Nu = 0.3 + ((0.62 * Math.sqrt(Re) * p.Pr ** (1 / 3)) / (1 + (0.4 / p.Pr) ** (2 / 3)) ** 0.25) * (1 + (Re / 282000) ** 0.625) ** 0.8;
  return (Nu * p.k) / Do;
}
/**
 * Overall heat-transfer coefficient referred to the inner diameter (W/m²/K).
 * { id, wt, kWall, layers: [{ t (m), k (W/m/K) }], hIn, hOut, burial: { depth (m to pipe centre), kSoil } | null }
 * Returns { U, resistances: [{ name, R (m²K/W referred to ID) }], od }.
 */
export function uValue({ id, wt = 0.0159, kWall = 45, layers = [], hIn = 1500, hOut = 500, burial = null }) {
  const ri = id / 2, res = [{ name: 'Inside film', R: 1 / Math.max(hIn, 1e-6) }];
  let r = ri;
  const add = (name, t, k) => { if (t > 0 && k > 0) { res.push({ name, R: (ri * Math.log((r + t) / r)) / k }); r += t; } };
  add('Pipe wall', wt, kWall);
  layers.forEach((l, i) => add(l.name || `Layer ${i + 1}`, +l.t || 0, +l.k || 0));
  const od = 2 * r;
  if (burial && burial.depth > od / 2 && burial.kSoil > 0) res.push({ name: 'Soil (buried)', R: (ri * Math.acosh((2 * burial.depth) / od)) / burial.kSoil });
  else res.push({ name: 'Outside film', R: ri / (r * Math.max(hOut, 1e-6)) });
  return { U: 1 / res.reduce((s, x) => s + x.R, 0), resistances: res, od };
}
/** Sea temperature (°C) at a depth (m, positive down): exponential thermocline between surface and seabed values. */
export const seaTemperature = (depth, tSurface = 24, tBed = 4, scale = 250) => tBed + (tSurface - tBed) * Math.exp(-Math.max(depth, 0) / scale);

// ---- stratified-flow geometry and equilibrium level ----------------------------------------------------
function stratGeom(hD, D) {
  const c = clamp(2 * hD - 1, -1, 1), A = (Math.PI * D * D) / 4, AL = 0.25 * D * D * (Math.PI - Math.acos(c) + c * Math.sqrt(1 - c * c));
  return { A, AL, AG: A - AL, SL: D * (Math.PI - Math.acos(c)), SG: D * Math.acos(c), Si: D * Math.sqrt(1 - c * c) };
}
const fanning = (Re) => (Re < 2100 ? 16 / Math.max(Re, 1e-9) : 0.046 * Re ** -0.2);
/** Equilibrium stratified liquid level (Taitel & Dukler momentum balance). theta in radians, positive upward. */
export function stratifiedLevel({ vsl, vsg, rhoL, rhoG, muL, muG, D, theta = 0 }) {
  const bal = (hD) => {
    const g = stratGeom(hD, D), vL = (vsl * g.A) / g.AL, vG = (vsg * g.A) / g.AG, DL = (4 * g.AL) / g.SL, DG = (4 * g.AG) / (g.SG + g.Si);
    const fL = fanning((rhoL * vL * DL) / muL), fG = fanning((rhoG * vG * DG) / muG), fi = Math.max(fG, 0.0142);
    const tL = (fL * rhoL * vL * vL) / 2, tG = (fG * rhoG * vG * vG) / 2, ti = (fi * rhoG * (vG - vL) * Math.abs(vG - vL)) / 2;
    return { r: (tG * g.SG) / g.AG - (tL * g.SL) / g.AL + ti * g.Si * (1 / g.AL + 1 / g.AG) - (rhoL - rhoG) * G * Math.sin(theta), g, vL, vG, tL, tG, ti, fL, fG };
  };
  // scan from the bottom for the first sign change (the thin-film root is the stable one when several exist)
  let a = 1e-4, fa = bal(a).r, hD = null;
  for (let i = 1; i <= 80; i++) { const b = 1e-4 + ((0.9995 - 1e-4) * i) / 80, fb = bal(b).r; if (fa * fb <= 0) { hD = brent((x) => bal(x).r, a, b, 1e-10); break; } a = b; fa = fb; }
  if (hD === null) hD = fa < 0 ? 0.9995 : 1e-4;
  const s = bal(hD);
  return { hD, holdup: s.g.AL / s.g.A, vL: s.vL, vG: s.vG, tauWL: s.tL, tauWG: s.tG, tauI: s.ti, geom: s.g };
}

/**
 * Flow pattern. Returns { pattern, strat } with pattern one of 'stratified smooth', 'stratified wavy', 'slug',
 * 'annular', 'dispersed bubble', 'bubble', 'churn', 'single-phase liquid', 'single-phase gas'.
 */
export function flowPattern({ vsl, vsg, rhoL, rhoG, muL, muG, sigma = 0.02, D, theta = 0 }) {
  if (vsg <= 1e-9) return { pattern: 'single-phase liquid' };
  if (vsl <= 1e-9) return { pattern: 'single-phase gas' };
  const vm = vsl + vsg, dRho = Math.max(rhoL - rhoG, 1e-6);
  if (theta > (10 * Math.PI) / 180) { // steeply upward: Taitel–Barnea–Dukler
    const vAnn = (3.1 * (sigma * G * dRho) ** 0.25) / Math.sqrt(rhoG);
    if (vsg >= vAnn && vsl / vm < 0.24) return { pattern: 'annular' }; // a thick film bridges the pipe: no stable annular flow above ~24 % liquid
    const vDB = 4 * ((D ** 0.429 * (sigma / rhoL) ** 0.089) / (muL / rhoL) ** 0.072) * ((G * dRho) / rhoL) ** 0.446;
    if (vm >= vDB && vsg / vm < 0.52) return { pattern: 'dispersed bubble' };
    const vBub = 1.53 * ((G * sigma * dRho) / (rhoL * rhoL)) ** 0.25;
    if (D > 19 * Math.sqrt((sigma * dRho) / (rhoL * rhoL * G)) && vsg < (vsl + 1.15 * vBub * Math.sin(theta)) / 3) return { pattern: 'bubble' };
    return { pattern: vsg > 0.6 * vAnn ? 'churn' : 'slug' };
  }
  const st = stratifiedLevel({ vsl, vsg, rhoL, rhoG, muL, muG, D, theta }), g = st.geom, cosT = Math.max(Math.cos(theta), 0.02);
  const stable = st.vG < (1 - st.hD) * Math.sqrt((dRho * G * cosT * g.AG) / (rhoG * Math.max(g.Si, 1e-9)));
  if (stable && st.hD < 0.999) {
    const wavy = st.vG >= Math.sqrt((4 * muL * dRho * G * cosT) / (0.01 * rhoL * rhoG * Math.max(st.vL, 1e-9)));
    return { pattern: wavy || theta < -0.005 ? 'stratified wavy' : 'stratified smooth', strat: st };
  }
  if (st.hD < 0.35) return { pattern: 'annular', strat: st };
  const DL = (4 * g.AL) / g.SL, fL = fanning((rhoL * st.vL * DL) / muL);
  if (st.vL >= Math.sqrt((4 * g.AG * G * cosT * (1 - rhoG / rhoL)) / (Math.max(g.Si, 1e-9) * fL))) return { pattern: 'dispersed bubble', strat: st };
  return { pattern: 'slug', strat: st };
}

// ---- slug closures --------------------------------------------------------------------------------------
/** Bendiksen (1984) slug translational velocity: { C0, vd, vt }. */
export function slugVelocity(vm, D, theta = 0) {
  const Fr = vm / Math.sqrt(G * D), s = Math.sin(theta), c = Math.max(Math.cos(theta), 0);
  const C0 = Fr < 3.5 ? 1.05 + 0.15 * s * s : 1.2, vd = Fr < 3.5 ? Math.sqrt(G * D) * (0.35 * s + 0.54 * c) : 0.35 * Math.sqrt(G * D) * s;
  return { C0, vd, vt: Math.max(C0 * vm + vd, 1e-6) };
}
/** Liquid holdup in the slug body (Gregory et al., 1978). */
export const slugBodyHoldup = (vm) => 1 / (1 + (vm / 8.66) ** 1.39);
/** Slug frequency (1/s). model: 'zabaras' (inclination-corrected) | 'gregory' (Gregory & Scott). */
export function slugFrequency(vsl, vm, D, theta = 0, model = 'zabaras') {
  if (vm <= 0 || vsl <= 0) return 0;
  if (model === 'gregory') return 0.0226 * ((vsl / (G * D)) * (19.75 / vm + vm)) ** 1.2;
  const ft = 3.28084, vslF = vsl * ft, vmF = vm * ft, dF = D * ft;
  return 0.0226 * (vslF / (32.174 * dF)) ** 1.2 * (212.6 / vmF + vmF) ** 1.2 * (0.836 + 2.75 * Math.sin(Math.max(theta, 0)) ** 0.25);
}
/** Mean slug-body length (m): Scott et al. (1989) for large pipes, 32 diameters for small ones. model 'brill' uses Brill et al. (1981). */
export function slugLength(D, vm = 3, model = 'scott') {
  const dIn = D / 0.0254;
  if (model === 'brill') return Math.exp(-2.663 + 5.441 * Math.sqrt(Math.log(dIn)) + 0.059 * Math.log(Math.max(vm * 3.28084, 0.1))) * 0.3048;
  return D < 0.1 ? 32 * D : Math.max(32 * D, Math.exp(-26.6 + 28.5 * (Math.log(D) + 3.67) ** 0.1)); // D and length in metres
}
/**
 * Hydrodynamic slug unit-cell summary at one location.
 * Returns { vt, C0, vd, holdupSlug, holdupFilm, holdup (unit average), freq (1/s), length (m, developed-slug correlation), lengthFromFreq (m, consistent with freq), lengthMax, unitLength, slugFraction, volume (m³ liquid per slug), period (s) }.
 */
export function slugUnit({ vsl, vsg, rhoL, rhoG, muL, muG, D, theta = 0, freqModel = 'zabaras', lengthModel = 'scott' }) {
  const vm = vsl + vsg, A = (Math.PI * D * D) / 4, { vt, C0, vd } = slugVelocity(vm, D, theta), HLS = slugBodyHoldup(vm);
  const vGb = 1.2 * vm + 1.53 * ((G * 0.02 * Math.max(rhoL - rhoG, 1)) / (rhoL * rhoL)) ** 0.25 * Math.sqrt(HLS) * Math.sin(Math.max(theta, 0));
  const holdup = clamp((vt * HLS + vGb * (1 - HLS) - vsg) / vt, vsl / Math.max(vm, 1e-9), 1);
  let HLF = clamp(stratifiedLevel({ vsl: Math.max(vsl * 0.3, 1e-4), vsg, rhoL, rhoG, muL, muG, D, theta: Math.min(theta, 0.15) }).holdup, 0.01, 0.9 * HLS);
  HLF = Math.min(HLF, holdup * 0.98);
  const beta = clamp((holdup - HLF) / Math.max(HLS - HLF, 1e-6), 0.02, 1), length = slugLength(D, vm, lengthModel);
  const fModel = slugFrequency(vsl, vm, D, theta, freqModel), unitLength = length / beta, freq = fModel > 0 ? fModel : vt / unitLength;
  return { vt, C0, vd, holdupSlug: HLS, holdupFilm: HLF, holdup, freq, length, lengthFromFreq: freq > 0 ? (beta * vt) / freq : length, lengthMax: length * Math.exp(3.09 * 0.5 - 0.125), unitLength, slugFraction: beta, volume: length * A * HLS, period: freq > 0 ? 1 / freq : Infinity };
}
/**
 * Severe (riser-induced) slugging screening. Bøe criterion and the Pots number.
 * { vsl, vsg (at riser-base conditions), P (Pa at riser base), rhoL, alpha (gas fraction in the flowline), L (flowline length feeding the riser), thetaRiser (rad), mG, mL, T (K), zG, mwG (g/mol) }
 * Returns { boeVsl (liquid superficial velocity above which severe slugging is possible), severe, pots }.
 */
export function severeSlugging({ vsl, vsg, P, rhoL, alpha = 0.5, L, thetaRiser = Math.PI / 2, mG, mL, T, zG = 0.9, mwG = 20 }) {
  const boeVsl = (P * vsg) / (rhoL * G * Math.max(alpha, 0.02) * L * Math.max(Math.sin(thetaRiser), 0.05));
  const pots = mG !== undefined && mL > 0 ? ((zG * 8.314462618 * T) / (mwG * 1e-3)) * mG / (G * Math.max(alpha, 0.02) * L * mL) : null;
  return { boeVsl, severe: vsl >= boeVsl, pots };
}

// ---- holdup and pressure-gradient closures --------------------------------------------------------------
function beggsBrill(p) {
  const { vsl, vsg, rhoL, rhoG, muL, muG, sigma, D, theta, rough, P, fModel } = p, vm = vsl + vsg, lam = clamp(vsl / vm, 1e-6, 1), Fr = (vm * vm) / (G * D);
  const L1 = 316 * lam ** 0.302, L2 = 0.0009252 * lam ** -2.4684, L3 = 0.1 * lam ** -1.4516, L4 = 0.5 * lam ** -6.738;
  let reg = 'distributed';
  if ((lam < 0.01 && Fr < L1) || (lam >= 0.01 && Fr < L2)) reg = 'segregated';
  else if (lam >= 0.01 && Fr >= L2 && Fr <= L3) reg = 'transition';
  else if ((lam >= 0.01 && lam < 0.4 && Fr > L3 && Fr <= L1) || (lam >= 0.4 && Fr > L3 && Fr <= L4)) reg = 'intermittent';
  const NLv = vsl * (rhoL / (G * Math.max(sigma, 1e-4))) ** 0.25, up = theta >= 0;
  const HL = (r) => {
    const [a, b, c] = { segregated: [0.98, 0.4846, 0.0868], intermittent: [0.845, 0.5351, 0.0173], distributed: [1.065, 0.5824, 0.0609] }[r];
    const H0 = Math.max((a * lam ** b) / Fr ** c, lam);
    let C = 0;
    if (!up) C = (1 - lam) * Math.log(4.7 * lam ** -0.3692 * NLv ** 0.1244 * Fr ** -0.5056);
    else if (r === 'segregated') C = (1 - lam) * Math.log(0.011 * lam ** -3.768 * NLv ** 3.539 * Fr ** -1.614);
    else if (r === 'intermittent') C = (1 - lam) * Math.log(2.96 * lam ** 0.305 * NLv ** -0.4473 * Fr ** 0.0978);
    C = Math.max(C, 0);
    const a18 = 1.8 * theta, psi = 1 + C * (Math.sin(a18) - Math.sin(a18) ** 3 / 3);
    return clamp(H0 * psi * (up ? 0.924 : 0.685), up ? lam : 1e-4, 1); // Payne et al. corrections
  };
  let holdup;
  if (reg === 'transition') { const Aw = (L3 - Fr) / (L3 - L2); holdup = Aw * HL('segregated') + (1 - Aw) * HL('intermittent'); } else holdup = HL(reg);
  const rhoN = rhoL * lam + rhoG * (1 - lam), muN = muL * lam + muG * (1 - lam), fn = frictionFactor((rhoN * vm * D) / muN, rough / D, fModel);
  const y = lam / (holdup * holdup), ly = Math.log(y), s = y > 1 && y < 1.2 ? Math.log(2.2 * y - 1.2) : ly / (-0.0523 + 3.182 * ly - 0.8725 * ly * ly + 0.01853 * ly ** 4);
  const rhoS = rhoL * holdup + rhoG * (1 - holdup), fric = (fn * Math.exp(Number.isFinite(s) ? clamp(s, -1, 2) : 0) * rhoN * vm * vm) / (2 * D), grav = rhoS * G * Math.sin(theta);
  const Ek = clamp((rhoS * vm * vsg) / Math.max(P, 1e4), 0, 0.6);
  return { holdup, fric, grav, acc: ((fric + grav) * Ek) / (1 - Ek), regime: reg === 'segregated' ? 'stratified wavy' : reg === 'distributed' ? (lam > 0.5 ? 'dispersed bubble' : 'annular') : 'slug', tauW: (fric * D) / 4 };
}
function driftFlux(p, forced) {
  const { vsl, vsg, rhoL, rhoG, muL, muG, sigma, D, theta, rough, P, fModel } = p, vm = vsl + vsg;
  let { C0, vd } = slugVelocity(vm, D, theta);
  if (forced === 'annular') { C0 = 1; vd = 0; }
  else if (forced === 'dispersed bubble' || forced === 'bubble') { C0 = 1.2; vd = 1.53 * ((G * Math.max(sigma, 1e-4) * Math.max(rhoL - rhoG, 1)) / (rhoL * rhoL)) ** 0.25 * Math.sin(theta); }
  const alpha = clamp(vsg / Math.max(C0 * vm + vd, 1e-9), 0, 1 - 1e-4), holdup = clamp(1 - alpha, Math.min(vsl / vm, 0.9999), 1);
  const rhoM = rhoL * holdup + rhoG * (1 - holdup), muM = muL * holdup + muG * (1 - holdup), f = frictionFactor((rhoM * vm * D) / muM, rough / D, fModel);
  const fric = (f * rhoM * vm * vm) / (2 * D), grav = rhoM * G * Math.sin(theta), Ek = clamp((rhoM * vm * vsg) / Math.max(P, 1e4), 0, 0.6);
  return { holdup, fric, grav, acc: ((fric + grav) * Ek) / (1 - Ek), regime: forced || 'slug', tauW: (fric * D) / 4 };
}
function mechanistic(p) {
  const { vsl, vsg, rhoL, rhoG, muL, muG, D, theta, rough, P, fModel } = p, fp = flowPattern(p), vm = vsl + vsg;
  if (fp.pattern.startsWith('stratified')) {
    const st = fp.strat, g = st.geom, holdup = st.holdup, fricG = (st.tauWG * g.SG + st.tauI * g.Si) / g.AG, fricTot = (st.tauWL * g.SL + st.tauWG * g.SG) / g.A;
    const rhoM = rhoL * holdup + rhoG * (1 - holdup);
    return { holdup, fric: fricTot, grav: rhoM * G * Math.sin(theta), acc: 0, regime: fp.pattern, tauW: st.tauWL, fricGas: fricG };
  }
  if (fp.pattern === 'slug' || fp.pattern === 'churn') {
    const u = slugUnit(p), rhoS = rhoL * u.holdupSlug + rhoG * (1 - u.holdupSlug), muS = muL * u.holdupSlug + muG * (1 - u.holdupSlug), f = frictionFactor((rhoS * vm * D) / muS, rough / D, fModel);
    const rhoU = rhoL * u.holdup + rhoG * (1 - u.holdup), fric = ((f * rhoS * vm * vm) / (2 * D)) * u.slugFraction + ((0.02 * rhoG * vm * vm) / (2 * D)) * (1 - u.slugFraction);
    const grav = rhoU * G * Math.sin(theta), Ek = clamp((rhoU * vm * vsg) / Math.max(P, 1e4), 0, 0.6);
    return { holdup: u.holdup, fric, grav, acc: ((fric + grav) * Ek) / (1 - Ek), regime: fp.pattern, tauW: (f * rhoS * vm * vm) / 8, slug: u };
  }
  return driftFlux(p, fp.pattern);
}
/**
 * Two-phase holdup and pressure gradient at one location.
 * p: { vsl, vsg, rhoL, rhoG, muL, muG, sigma, D, theta (rad, + up), rough (m), P (Pa), fModel }
 * model: 'beggsBrill' | 'driftFlux' | 'mechanistic' | 'homogeneous'
 * Returns { holdup, fric, grav, acc (Pa/m, positive = pressure falls in the flow direction), dpdx (total), regime, tauW (Pa) }.
 */
export function gradient(p, model = 'beggsBrill') {
  const q = { sigma: 0.02, rough: 4.5e-5, theta: 0, P: 1e7, fModel: 'colebrook', ...p }, vm = q.vsl + q.vsg;
  let r;
  if (vm <= 1e-9) { const liquid = q.vsg <= 1e-12; r = { holdup: liquid ? 1 : 0, fric: 0, grav: (liquid ? q.rhoL : q.rhoG) * G * Math.sin(q.theta), acc: 0, regime: 'static', tauW: 0 }; }
  else if (q.vsg <= 1e-9 * vm || q.vsl <= 1e-9 * vm) {
    const liquid = q.vsg <= 1e-9 * vm, rho = liquid ? q.rhoL : q.rhoG, mu = liquid ? q.muL : q.muG, f = frictionFactor((rho * vm * q.D) / mu, q.rough / q.D, q.fModel), fric = (f * rho * vm * vm) / (2 * q.D);
    r = { holdup: liquid ? 1 : 0, fric, grav: rho * G * Math.sin(q.theta), acc: 0, regime: liquid ? 'single-phase liquid' : 'single-phase gas', tauW: (fric * q.D) / 4 };
  } else if (model === 'homogeneous') {
    const lam = q.vsl / vm, rho = q.rhoL * lam + q.rhoG * (1 - lam), mu = q.muL * lam + q.muG * (1 - lam), f = frictionFactor((rho * vm * q.D) / mu, q.rough / q.D, q.fModel), fric = (f * rho * vm * vm) / (2 * q.D);
    r = { holdup: lam, fric, grav: rho * G * Math.sin(q.theta), acc: 0, regime: 'homogeneous', tauW: (fric * q.D) / 4 };
  } else r = model === 'driftFlux' ? driftFlux(q) : model === 'mechanistic' ? mechanistic(q) : beggsBrill(q);
  if (model === 'beggsBrill' || model === 'driftFlux') { try { const fp = flowPattern(q).pattern; if (!fp.startsWith('single')) r.regime = fp; } catch { /* keep the correlation's own regime label */ } }
  r.dpdx = r.fric + r.grav + r.acc;
  return r;
}

// ---- elevation profile helpers ---------------------------------------------------------------------------
/** Resample an elevation profile { x[], z[] } (horizontal distance, elevation) to n cells of equal arc length. Returns node arrays { s, x, z, theta (per cell, rad), ds }. */
export function discretise(profile, n = 200) {
  const X = profile.x, Z = profile.z, S = [0];
  for (let i = 1; i < X.length; i++) S.push(S[i - 1] + Math.hypot(X[i] - X[i - 1], Z[i] - Z[i - 1]));
  const L = S[S.length - 1], s = [], x = [], z = [], theta = [];
  for (let i = 0; i <= n; i++) { const si = (L * i) / n; s.push(si); x.push(interp1(S, X, si)); z.push(interp1(S, Z, si)); }
  for (let i = 0; i < n; i++) theta.push(Math.asin(clamp((z[i + 1] - z[i]) / (s[i + 1] - s[i] || 1), -1, 1)));
  return { s, x, z, theta, ds: L / n, length: L, n };
}

/**
 * Steady-state pressure and temperature along a pipe by marching from the inlet (shooting on the inlet pressure
 * when the outlet pressure is the boundary condition).
 * o: { fm (fluidModel), profile: { x[], z[] }, id (m), rough (m), U (W/m²/K on ID) | uOf(s), tAmb (°C) | tAmbOf(s, z), tIn (°C),
 *      pOut (bara) | pIn (bara), mScale (rate multiplier), model, fModel, n, idOf(s) (effective inner diameter), roughOf(s) }
 * Returns { ok, s, x, z, theta, P, T, holdup, vsl, vsg, vm, rhoM, dpdx, regime, tauW, tAmb, tHyd, subcooling, qG, qL,
 *           pIn, pOut, tOut, dpFric, dpGrav, liquidInventory (m³), volume (m³), residence (s), heatLoss (W), mdot }.
 */
export function marchSteady(o) {
  const { fm, profile, id, rough = 4.5e-5, tIn = 70, mScale = 1, model = 'beggsBrill', fModel = 'colebrook', n = 200 } = o, grid = discretise(profile, n);
  const uOf = typeof o.uOf === 'function' ? o.uOf : () => o.U ?? 3, tAmbOf = typeof o.tAmbOf === 'function' ? o.tAmbOf : () => o.tAmb ?? 4;
  const dOf = typeof o.idOf === 'function' ? o.idOf : () => id, rOf = typeof o.roughOf === 'function' ? o.roughOf : () => rough;
  const mdot = (fm.rates.mHC + fm.rates.mW) * mScale;
  const state = (P, T, i) => {
    const sMid = 0.5 * (grid.s[i] + grid.s[i + 1]), D = Math.max(dOf(sMid), 0.01), A = (Math.PI * D * D) / 4, pr = fm.at(P, T, mScale);
    const vsg = pr.qG / A, vsl = pr.qL / A, gr = gradient({ vsl, vsg, rhoL: pr.rhoL, rhoG: pr.rhoG, muL: pr.muL, muG: pr.muG, sigma: pr.sigma, D, theta: grid.theta[i], rough: rOf(sMid), P: P * 1e5, fModel }, model);
    const mCp = pr.mG * pr.cpG + pr.mO * pr.cpO + pr.mW * pr.cpW, jt = mCp > 0 ? (pr.mG * pr.cpG * pr.jtG + pr.mO * pr.cpO * pr.jtO - (pr.mW / pr.rhoW)) / mCp : 0;
    const ta = tAmbOf(sMid, 0.5 * (grid.z[i] + grid.z[i + 1])), U = uOf(sMid), q = U * Math.PI * D * (T - ta); // W/m
    const dTds = mCp > 0 ? -q / mCp - jt * gr.dpdx - (mdot * G * Math.sin(grid.theta[i])) / mCp : 0;
    return { pr, gr, vsg, vsl, D, A, ta, q, dPds: -gr.dpdx / 1e5, dTds };
  };
  const run = (pIn) => {
    const P = [pIn], T = [tIn], cells = [];
    let ok = true;
    for (let i = 0; i < grid.n; i++) {
      const a = state(P[i], T[i], i), Pm = P[i] + 0.5 * grid.ds * a.dPds, Tm = T[i] + 0.5 * grid.ds * a.dTds;
      if (!(Pm > 1.0)) { ok = false; break; }
      const b = state(Pm, Tm, i), Pn = P[i] + grid.ds * b.dPds, Tn = T[i] + grid.ds * b.dTds;
      if (!(Pn > 1.0) || !Number.isFinite(Tn)) { ok = false; break; }
      P.push(Pn); T.push(clamp(Tn, -60, 250)); cells.push(b);
    }
    return { ok, P, T, cells };
  };
  let sol, pIn = o.pIn;
  if (o.pOut !== undefined && o.pOut !== null && o.pIn === undefined) {
    const res = (p) => { const r = run(p); return r.ok ? r.P[r.P.length - 1] - o.pOut : -1e3 - (grid.n - r.cells.length); };
    let lo = Math.max(o.pOut, 1.5), hi = lo + 20, fhi = res(hi), guard = 0;
    while (fhi < 0 && hi < 1400 && guard++ < 40) { lo = hi; hi = hi * 1.35 + 10; fhi = res(hi); }
    if (fhi < 0) return { ok: false, reason: 'No inlet pressure below 1,400 bara can deliver this rate: the line is too small, too long or blocked.' };
    for (let k = 0; k < 60 && hi - lo > 1e-6 * hi; k++) { const m = 0.5 * (lo + hi); if (res(m) < 0) lo = m; else hi = m; }
    pIn = hi;
  }
  sol = run(pIn);
  if (!sol.ok) return { ok: false, reason: 'The inlet pressure is too low to push this rate to the outlet.' };
  const N = grid.n, out = { ok: true, s: grid.s, x: grid.x, z: grid.z, theta: grid.theta.concat(grid.theta[N - 1]), P: sol.P, T: sol.T, mdot, pIn: sol.P[0], pOut: sol.P[N], tOut: sol.T[N], ds: grid.ds, length: grid.length };
  const col = (fn) => { const a = sol.cells.map(fn); a.push(a[a.length - 1]); return a; };
  out.holdup = col((c) => c.gr.holdup); out.vsl = col((c) => c.vsl); out.vsg = col((c) => c.vsg); out.vm = col((c) => c.vsl + c.vsg); out.dpdx = col((c) => c.gr.dpdx); out.regime = col((c) => c.gr.regime); out.tauW = col((c) => c.gr.tauW);
  out.rhoM = col((c) => c.pr.rhoL * c.gr.holdup + c.pr.rhoG * (1 - c.gr.holdup)); out.tAmb = col((c) => c.ta); out.qG = col((c) => c.pr.qG); out.qL = col((c) => c.pr.qL); out.rhoL = col((c) => c.pr.rhoL); out.rhoG = col((c) => c.pr.rhoG); out.muL = col((c) => c.pr.muL); out.wcut = col((c) => c.pr.wcut);
  out.tHyd = out.P.map((p) => fm.hydrateT(p)); out.subcooling = out.T.map((t, i) => out.tHyd[i] - t);
  let inv = 0, vol = 0, fr = 0, gv = 0, heat = 0, res = 0;
  sol.cells.forEach((c) => { inv += c.gr.holdup * c.A * grid.ds; vol += c.A * grid.ds; fr += (c.gr.fric + c.gr.acc) * grid.ds; gv += c.gr.grav * grid.ds; heat += c.q * grid.ds; res += grid.ds / Math.max(c.vsl + c.vsg, 1e-6); });
  Object.assign(out, { liquidInventory: inv, volume: vol, dpFric: fr / 1e5, dpGrav: gv / 1e5, heatLoss: heat, residence: res });
  return out;
}
