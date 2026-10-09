// Pipe-flow kernel shared by the suites: friction factors, heat-transfer resistances, two-phase flow-pattern
// prediction at any inclination (Taitel–Dukler / Barnea unified model), liquid holdup and pressure-gradient closures (Beggs–Brill,
// drift flux, mechanistic stratified / slug unit cell with the Taitel–Barnea film), slug closures (frequency, development-aware
// length), severe-slugging criteria, sea-temperature
// profile and a steady-state marching solver for pressure and temperature along an elevation profile.
// SI units throughout except pressure in bara and temperature in °C at the interface of marchSteady.
// Where each constant set was checked against an openly readable source is listed in PROVENANCE of js/suites/s03_flow.js.
import { clamp, brent, interp1 } from './num.js';
import { density as swDensity, viscosity as swViscosity, cp as swCp, conductivityThermal as swConductivity } from './props.js';

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
/**
 * Outside film coefficient (W/m²/K) of a cylinder in cross-flow (Churchill & Bernstein, 1977). medium: 'seawater' | 'air';
 * T = ambient temperature (°C). Seawater properties from core/props.js at 35 g/kg (Sharqawy et al.); air at 1 atm, 20 °C.
 * With dT (surface minus ambient temperature, K) the free convection of a horizontal cylinder (Churchill & Chu, 1975) is
 * combined with the forced convection as Nu³ = Nu_forced³ + Nu_free³; dT = 0 (default) gives forced convection alone.
 */
export function hOutside(v, Do, medium = 'seawater', T = 4, dT = 0) {
  let p;
  if (medium === 'air') p = { rho: 1.2, mu: 1.8e-5, k: 0.026, Pr: 0.71, beta: 1 / (T + 273.15) };
  else { const rho = swDensity(T, 35), mu = swViscosity(T, 35), k = swConductivity(T, 35); p = { rho, mu, k, Pr: (swCp(T, 35) * mu) / k, beta: Math.max((swDensity(T - 0.5, 35) - swDensity(T + 0.5, 35)) / rho, 1e-6) }; }
  const NuF = nuCrossFlow(Math.max((p.rho * Math.max(v, 0.01) * Do) / p.mu, 1), p.Pr);
  if (!(Math.abs(dT) > 0)) return (NuF * p.k) / Do;
  const NuN = nuFreeCylinder(((G * p.beta * Math.abs(dT) * Do ** 3 * p.rho * p.rho) / (p.mu * p.mu)) * p.Pr, p.Pr);
  return (Math.cbrt(NuF ** 3 + NuN ** 3) * p.k) / Do;
}
/** Nusselt number of a cylinder in cross-flow (Churchill & Bernstein, 1977). */
export const nuCrossFlow = (Re, Pr) => 0.3 + ((0.62 * Math.sqrt(Re) * Pr ** (1 / 3)) / (1 + (0.4 / Pr) ** (2 / 3)) ** 0.25) * (1 + (Re / 282000) ** 0.625) ** 0.8;
/** Nusselt number of free convection around a horizontal cylinder (Churchill & Chu, 1975); Ra = Gr Pr on the diameter. */
export const nuFreeCylinder = (Ra, Pr) => (0.6 + (0.387 * Ra ** (1 / 6)) / (1 + (0.559 / Pr) ** (9 / 16)) ** (8 / 27)) ** 2;
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

// section geometry of a stratified layer as a function of the liquid holdup (wetted wall arc and interface chord, both per diameter)
const SECT = (() => { const M = 400, SL = new Float64Array(M + 1), Si = new Float64Array(M + 1); for (let j = 0; j <= M; j++) { const H = j / M, hD = j === 0 ? 0 : j === M ? 1 : brent((x) => { const c = 2 * x - 1; return (Math.PI - Math.acos(c) + c * Math.sqrt(1 - c * c)) / Math.PI - H; }, 0, 1, 1e-13), c = 2 * hD - 1; SL[j] = Math.PI - Math.acos(c); Si[j] = Math.sqrt(Math.max(1 - c * c, 0)); } return { M, SL, Si }; })();
const sect = (t, H) => { const u = clamp(H, 0, 1) * SECT.M, j = Math.min(Math.floor(u), SECT.M - 1); return t[j] + (t[j + 1] - t[j]) * (u - j); };

/**
 * Annular-flow film of the Barnea (1986, 1987) annular ↔ intermittent transition: film holdup H from the combined momentum balance
 * Y = (1 + 75 H) / [(1 − H)^2.5 H] − X² / H³ (interfacial friction f_i = f_G (1 + 300 δ/D), i.e. 1 + 75 H), X² = (dp/dx)_SL / (dp/dx)_SG,
 * Y = (ρL − ρG) g sinθ / (dp/dx)_SG. The film is unstable when Y ≥ (2 − 1.5 H) X² / [H³ (1 − 1.5 H)]; the gas core is blocked when H ≥ 0.24
 * (half of the maximum packing 0.48 of the liquid in a slug). Returns { H, stable, blocked, X2, Y }.
 */
export function annularFilm({ vsl, vsg, rhoL, rhoG, muL, muG, D, theta = 0 }) {
  const dpL = (2 * fanning((rhoL * vsl * D) / muL) * rhoL * vsl * vsl) / D, dpG = (2 * fanning((rhoG * vsg * D) / muG) * rhoG * vsg * vsg) / D, X2 = dpL / dpG, Y = ((rhoL - rhoG) * G * Math.sin(theta)) / dpG;
  const F = (h) => (1 + 75 * h) / ((1 - h) ** 2.5 * h) - X2 / h ** 3 - Y;
  let a = 1e-6, fa = F(a), H = null;
  for (let i = 1; i <= 96; i++) { const b = 10 ** (-6 + (6 * i) / 96) * 0.999, fb = F(b); if (fa * fb <= 0) { H = brent(F, a, b, 1e-12); break; } a = b; fa = fb; } // the thinnest film is the stable solution
  if (H === null) return { H: 1, stable: false, blocked: true, X2, Y };
  return { H, stable: H < 2 / 3 && Y < ((2 - 1.5 * H) / (H ** 3 * (1 - 1.5 * H))) * X2, blocked: H >= 0.24, X2, Y };
}

/**
 * Flow pattern at any inclination (unified model after Taitel & Dukler 1976, Taitel, Barnea & Dukler 1980 and Barnea 1986, 1987).
 * Order of the tests: (1) stratified flow exists (not steeper upward than +10°) when the equilibrium level is stable against the
 * Kelvin–Helmholtz criterion of Taitel & Dukler and, in downward flow, the film is slower than the velocity at which liquid is torn off
 * and carried to the upper wall, v_L < [g D (1 − h/D) cosθ / f_L]^½ (Barnea, Shoham & Taitel 1982); it is wavy above the wind criterion
 * (sheltering coefficient 0.01) or, downward, above a film Froude number of 1.5; (2) dispersed bubbles when turbulence keeps the bubbles
 * smaller than both the deformation and the creaming size (Barnea 1986) at a gas fraction below 0.52, or, near the horizontal, when the
 * turbulent fluctuations of the liquid overcome buoyancy (Taitel & Dukler); (3) annular flow: near the horizontal (|θ| ≤ 10°) when the
 * equilibrium level is below 0.35 D, otherwise when the annular film is stable and does not block the core (annularFilm);
 * (4) bubble flow in steep upward pipes (≥ 60° from the horizontal) wide enough for small bubbles to rise slower than a Taylor bubble, at a
 * void fraction below 0.25; otherwise (5) intermittent: slug, or churn in upward flow above 0.6 of the annular (droplet-lifting) velocity.
 * Returns { pattern, strat, film } with pattern one of 'stratified smooth', 'stratified wavy', 'slug', 'annular', 'dispersed bubble',
 * 'bubble', 'churn', 'single-phase liquid', 'single-phase gas'.
 */
export function flowPattern({ vsl, vsg, rhoL, rhoG, muL, muG, sigma = 0.02, D, theta = 0 }) {
  if (vsg <= 1e-9) return { pattern: 'single-phase liquid' };
  if (vsl <= 1e-9) return { pattern: 'single-phase gas' };
  const vm = vsl + vsg, dRho = Math.max(rhoL - rhoG, 1e-6), sinT = Math.sin(theta), cosT = Math.cos(theta), deg = (theta * 180) / Math.PI;
  let st = null;
  if (deg <= 10 && deg > -89.99) { // (1) stratified
    st = stratifiedLevel({ vsl, vsg, rhoL, rhoG, muL, muG, D, theta }); const g = st.geom, c = Math.max(cosT, 1e-6), fL = fanning((rhoL * st.vL * 4 * g.AL) / (g.SL * muL));
    const stable = st.vG < (1 - st.hD) * Math.sqrt((dRho * G * c * g.AG) / (rhoG * Math.max(g.Si, 1e-9))), torn = theta < 0 && st.vL > Math.sqrt((G * D * (1 - st.hD) * c) / fL);
    if (stable && !torn && st.hD < 0.999) {
      const wind = st.vG >= Math.sqrt((4 * muL * dRho * G * c) / (0.01 * rhoL * rhoG * Math.max(st.vL, 1e-9))), froude = st.vL / Math.sqrt(G * Math.max(st.hD * D, 1e-9));
      return { pattern: wind || (theta < 0 && froude > 1.5) ? 'stratified wavy' : 'stratified smooth', strat: st };
    }
  }
  // (2) dispersed bubble
  const fm = fanning((rhoL * vm * D) / muL), dMax = (0.725 + 4.15 * Math.sqrt(vsg / vm)) * (sigma / rhoL) ** 0.6 * ((2 * fm * vm ** 3) / D) ** -0.4, dCD = 2 * Math.sqrt((0.4 * sigma) / (dRho * G)), dCB = (0.375 * rhoL * fm * vm * vm) / (dRho * G * Math.max(Math.abs(cosT), 1e-9));
  if (dMax < Math.min(dCD, dCB) && vsg / vm <= 0.52) return { pattern: 'dispersed bubble', strat: st };
  if (st && Math.abs(deg) <= 10 && st.hD >= 0.35) { const g = st.geom, fL = fanning((rhoL * st.vL * 4 * g.AL) / (g.SL * muL)); if (st.vL >= Math.sqrt((4 * g.AG * G * Math.max(cosT, 0.02) * (1 - rhoG / rhoL)) / (Math.max(g.Si, 1e-9) * fL))) return { pattern: 'dispersed bubble', strat: st }; }
  // (3) annular
  if (st && Math.abs(deg) <= 10) { if (st.hD < 0.35) return { pattern: 'annular', strat: st }; }
  else { const film = annularFilm({ vsl, vsg, rhoL, rhoG, muL, muG, D, theta }); if (film.stable && !film.blocked) return { pattern: 'annular', strat: st, film }; }
  // (4) bubble, (5) intermittent
  if (deg > 10) {
    const vBub = 1.53 * ((G * sigma * dRho) / (rhoL * rhoL)) ** 0.25;
    if (deg >= 60 && D > 19 * Math.sqrt((sigma * dRho) / (rhoL * rhoL * G)) && vsg < (vsl + 0.75 * vBub * sinT) / 3) return { pattern: 'bubble' }; // α = 0.25: vsl = 3 vsg − 0.75 v∞ sinθ (0.75 × 1.53 = 1.15)
    return { pattern: vsg > (0.6 * 3.1 * (sigma * G * dRho) ** 0.25) / Math.sqrt(rhoG) ? 'churn' : 'slug' };
  }
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
/**
 * Slug frequency (1/s) from a correlation of developed flow. model:
 *  'gregory'    Gregory & Scott (1969): 0.0226 [(vsl / gD)(19.75 / vm + vm)]^1.2 (SI; 19.75 m²/s²), horizontal, 19 and 35 mm;
 *  'zabaras'    Zabaras (2000): the Gregory–Scott value times (0.836 + 2.75 sin θ), fitted for 0 ≤ θ ≤ 11° (the factor is held beyond 11°);
 *  'greskovich' Greskovich & Shrier (1972): 0.0226 [λ (2.02 / D + vm² / gD)]^1.2;
 *  'heywood'    Heywood & Richardson (1979): 0.0434 [λ (2.02 / D + vm² / gD)]^1.02, λ = vsl / vm (2.02 in metres).
 */
export function slugFrequency(vsl, vm, D, theta = 0, model = 'zabaras') {
  if (vm <= 0 || vsl <= 0) return 0;
  if (model === 'heywood' || model === 'greskovich') { const x = (vsl / vm) * (2.02 / D + (vm * vm) / (G * D)); return model === 'heywood' ? 0.0434 * x ** 1.02 : 0.0226 * x ** 1.2; }
  const gs = 0.0226 * ((vsl / (G * D)) * (19.75 / vm + vm)) ** 1.2;
  return model === 'gregory' ? gs : gs * (0.836 + 2.75 * Math.sin(clamp(theta, 0, (11 * Math.PI) / 180)));
}
// Mean slug-body length (diameters) against the distance from the point of slug formation (diameters), computed once with the kinematic
// slug-train model of Barnea & Taitel (1993): slugs enter with random short lengths (uniform, mean 3.2 D horizontal / 4.2 D inclined), every bubble
// nose moves at v∞ [1 + a exp(−b L/D)] with L the slug ahead of it (Moissis & Griffith wake law; a = 0.4, b = 1.0 horizontal and a = 8, b = 1.06
// inclined and vertical, as used with this model in the open literature), a slug overtaken by its bubble disappears. Mean of 3 × 6,000 slugs.
const TRAIN = { x: [0, 10, 20, 30, 50, 75, 100, 150, 200, 300, 500, 750, 1000, 1500, 2000, 3000], hor: [3.2, 3.8, 4.2, 4.48, 4.9, 5.3, 5.6, 6.11, 6.47, 7.05, 7.93, 8.73, 9.34, 10.34, 10.79, 11.36], inc: [4.2, 5.95, 6.75, 7.25, 7.98, 8.61, 9.1, 9.9, 10.5, 11.49, 12.93, 14.34, 14.8, 15.2, 15.51, 15.88] };
/**
 * Mean slug-body length (m).
 * Developed flow (default): Scott, Shoham & Brill (1989) ln L = −26.6 + 28.5 [ln D + 3.67]^0.1 (L and D in metres) for large pipes, 32 diameters
 * for small ones; model 'brill' uses Brill et al. (1981), ln L[ft] = −2.663 + 5.441 √ln d[in] + 0.059 ln vm[ft/s].
 * With o.xD (distance from the inlet or from the point of slug formation, in diameters) the length is development-aware: up to 3,000 D it is the
 * mean length of the slug-train model (table above; horizontal and inclined curves blended with sin|θ| / sin 10°), i.e. 4–8 D at 50–100 D and
 * 11–16 D at 3,000 D; from there it approaches the developed-flow value on a logarithmic distance scale and reaches it at 30,000 D (this
 * blending range is an assumption: no open measurements of slug length between 10³ and 10⁵ diameters were found). o: { xD, theta }.
 */
export function slugLength(D, vm = 3, model = 'scott', o = null) {
  const dIn = D / 0.0254;
  const dev = model === 'brill' ? Math.exp(-2.663 + 5.441 * Math.sqrt(Math.log(Math.max(dIn, 1.01))) + 0.059 * Math.log(Math.max(vm * 3.28084, 0.1))) * 0.3048 : D < 0.1 ? 32 * D : Math.max(32 * D, Math.exp(-26.6 + 28.5 * (Math.log(D) + 3.67) ** 0.1));
  if (!o || !(o.xD >= 0)) return dev;
  const w = clamp(Math.abs(Math.sin(o.theta || 0)) / Math.sin((10 * Math.PI) / 180), 0, 1), x = Math.min(o.xD, 3000), lt = ((1 - w) * interp1(TRAIN.x, TRAIN.hor, x) + w * interp1(TRAIN.x, TRAIN.inc, x)) * D;
  return o.xD <= 3000 ? lt : lt + (dev - lt) * clamp(Math.log10(o.xD / 3000), 0, 1);
}
/**
 * Film zone of a slug unit (Taitel & Barnea 1990, equilibrium-film form). In the frame moving with the slug the liquid shed at the tail,
 * x = (v_t − v_LLS) H_LS, flows back through the film: v_LF = v_t − x / H_LF, and the gas of the bubble v_GF = v_t − (v_t − v_GLS)(1 − H_LS) / (1 − H_LF).
 * The film holdup is the level at which the film and the gas above it are in momentum equilibrium,
 *   τ_F S_F / A_F − τ_G S_G / A_G − τ_i S_i (1/A_F + 1/A_G) + (ρL − ρG) g sinθ = 0
 * (stratified section; a symmetric falling film steeper than 75°), the root reached first when the film thins from the slug tail.
 * The liquid balance then gives the slug fraction β = (H_U − H_LF) / (H_LS − H_LF) with the unit holdup H_U = (v_sl + x) / v_t.
 * p: { vsl, vsg, rhoL, rhoG, muL, muG, sigma, D, theta }; o: { vt, HLS } to override the closures.
 * Returns { vt, HLS, HLF, beta, HU, vLLS, vGLS, vLF, vGF, x, tauF, tauG (Pa, signed), SF, SG (m), wall (film-zone wall force per pipe volume, Pa/m) }.
 */
export function slugFilm(p, o = {}) {
  const { vsl, vsg, rhoL, rhoG, muL, muG, D, theta = 0, sigma = 0.02 } = p, vm = vsl + vsg, A = (Math.PI * D * D) / 4, vt = o.vt ?? slugVelocity(vm, D, theta).vt, HLS = clamp(o.HLS ?? slugBodyHoldup(vm), 0.05, 1), dRho = Math.max(rhoL - rhoG, 1e-6), sinT = Math.sin(theta);
  const vGLS = Math.min(1.2 * vm + 1.53 * ((G * sigma * Math.max(dRho, 1)) / (rhoL * rhoL)) ** 0.25 * Math.sqrt(HLS) * Math.sin(Math.max(theta, 0)), vt), vLLS = (vm - vGLS * (1 - HLS)) / HLS, x = Math.max((vt - vLLS) * HLS, 1e-9 * vt), y = (vt - vGLS) * (1 - HLS), ring = theta > (75 * Math.PI) / 180;
  const at = (H) => {
    const AF = H * A, AG = (1 - H) * A, SF = ring ? Math.PI * D : D * sect(SECT.SL, H), SG = ring ? 0 : Math.PI * D - SF, Si = ring ? Math.PI * D * Math.sqrt(1 - H) : D * sect(SECT.Si, H);
    const vLF = vt - x / H, vGF = vt - y / (1 - H), fF = fanning((rhoL * Math.abs(vLF) * 4 * AF) / (SF * muL)), fG = fanning((rhoG * Math.abs(vGF) * 4 * AG) / (Math.max(SG + Si, 1e-12) * muG)), fi = Math.max(fG, 0.0142);
    const tF = 0.5 * fF * rhoL * vLF * Math.abs(vLF), tG = 0.5 * fG * rhoG * vGF * Math.abs(vGF), ti = 0.5 * fi * rhoG * (vGF - vLF) * Math.abs(vGF - vLF);
    return { F: (tF * SF) / AF - (tG * SG) / AG - ti * Si * (1 / AF + 1 / AG) + dRho * G * sinT, tF, tG, SF, SG, vLF, vGF };
  };
  const hi = Math.min(HLS, 0.999) * 0.9999, lo = 1e-4, n = 40; let hb = hi, fb = at(hb).F, HLF = null;
  for (let i = 1; i <= n; i++) { const h = hi * (lo / hi) ** (i / n), fh = at(h).F; if (fb * fh <= 0) { HLF = brent((q) => at(q).F, h, hb, 1e-10); break; } hb = h; fb = fh; }
  if (HLF === null) HLF = Math.abs(at(lo).F) < Math.abs(at(hi).F) ? lo : hi;
  const e = at(HLF), HU = clamp((vsl + x) / vt, Math.min(vsl / Math.max(vm, 1e-9), 1), 1), beta = clamp((HU - HLF) / Math.max(HLS - HLF, 1e-9), 0, 1);
  return { vt, HLS, HLF, beta, HU, vLLS, vGLS, vLF: e.vLF, vGF: e.vGF, x, tauF: e.tF, tauG: e.tG, SF: e.SF, SG: e.SG, wall: (e.tF * e.SF + e.tG * e.SG) / A };
}
/**
 * Hydrodynamic slug unit-cell summary at one location. The unit cell obeys freq × unitLength = vt and length = slugFraction × unitLength,
 * so only one of slug length and slug frequency can come from a correlation; the other is derived:
 *   basis 'length' (default): the slug length (lengthModel; development-aware when xD, the distance from the inlet in diameters, is given) is kept and freq = slugFraction · vt / length;
 *   basis 'frequency': the frequency correlation (freqModel) is kept and length = slugFraction · vt / freq.
 * The value of the correlation that was not used is still reported (freqCorrelation, lengthCorrelation) for comparison.
 * Returns { vt, C0, vd, holdupSlug, holdupFilm, holdup (unit average), freq (1/s), period (s), length (m), lengthFromFreq (= length, kept for older callers),
 *           lengthMax (1-in-1000 slug of a log-normal distribution, σ = 0.5), unitLength, slugFraction, volume (m³ liquid per slug), freqCorrelation, lengthCorrelation, basis }.
 */
export function slugUnit({ vsl, vsg, rhoL, rhoG, muL, muG, D, theta = 0, freqModel = 'zabaras', lengthModel = 'scott', basis = 'length', sigma = 0.02, xD = null }) {
  const vm = vsl + vsg, A = (Math.PI * D * D) / 4, { vt, C0, vd } = slugVelocity(vm, D, theta), HLS = slugBodyHoldup(vm);
  const film = slugFilm({ vsl, vsg, rhoL, rhoG, muL, muG, D, theta, sigma }, { vt, HLS }), holdup = film.HU, HLF = clamp(film.HLF, 1e-4, Math.min(0.999 * HLS, holdup)); // film zone from the Taitel–Barnea equilibrium film
  const beta = clamp((holdup - HLF) / Math.max(HLS - HLF, 1e-6), 0.02, 1), lengthCorrelation = slugLength(D, vm, lengthModel, xD === null || xD === undefined ? null : { xD, theta }), freqCorrelation = slugFrequency(vsl, vm, D, theta, freqModel);
  const useF = basis === 'frequency' && freqCorrelation > 0, length = useF ? (beta * vt) / freqCorrelation : lengthCorrelation, freq = useF ? freqCorrelation : (beta * vt) / length;
  return { vt, C0, vd, holdupSlug: HLS, holdupFilm: HLF, holdup, freq, length, lengthFromFreq: length, lengthMax: length * Math.exp(3.09 * 0.5 - 0.125), unitLength: length / beta, slugFraction: beta, volume: length * A * HLS, period: freq > 0 ? 1 / freq : Infinity, freqCorrelation, lengthCorrelation, basis: useF ? 'frequency' : 'length', film };
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
    return clamp(H0 * psi * (theta > 0 ? 0.924 : theta < 0 ? 0.685 : 1), up ? lam : 1e-4, 1); // Payne et al. (1979) factors for uphill and downhill flow; the horizontal correlation is not rescaled
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
    // slug unit: wall friction of the slug body over the slug fraction β plus the wall forces of the film zone (liquid film and gas bubble) over 1 − β;
    // the static head carries the unit holdup. The film is the equilibrium film of slugFilm(), so a film falling back in upward flow lowers the gradient.
    const u = slugUnit(p), rhoS = rhoL * u.holdupSlug + rhoG * (1 - u.holdupSlug), muS = muL * u.holdupSlug + muG * (1 - u.holdupSlug), f = frictionFactor((rhoS * vm * D) / muS, rough / D, fModel);
    const rhoU = rhoL * u.holdup + rhoG * (1 - u.holdup), b = u.film.beta, fric = ((f * rhoS * vm * vm) / (2 * D)) * b + u.film.wall * (1 - b);
    const grav = rhoU * G * Math.sin(theta), Ek = clamp((rhoU * vm * vsg) / Math.max(P, 1e4), 0, 0.6);
    return { holdup: u.holdup, fric, grav, acc: ((fric + grav) * Ek) / (1 - Ek), regime: fp.pattern, tauW: (f * rhoS * vm * vm) / 8, slug: u };
  }
  return driftFlux(p, fp.pattern);
}
/**
 * Two-phase holdup and pressure gradient at one location.
 * p: { vsl, vsg, rhoL, rhoG, muL, muG, sigma, D, theta (rad, + up), rough (m), P (Pa), fModel }
 * model: 'beggsBrill' | 'driftFlux' | 'mechanistic' | 'homogeneous'
 * p.label = false skips the mechanistic flow-pattern label of the 'beggsBrill' and 'driftFlux' models (the costly part of a call:
 * the stratified equilibrium level has to be solved); the correlation's own regime name is then returned. Numbers are unaffected.
 * Returns { holdup, fric, grav, acc (Pa/m, positive = pressure falls in the flow direction), dpdx (total), regime, tauW (Pa) }.
 */
export function gradient(p, model = 'beggsBrill') {
  const q = { sigma: 0.02, rough: 4.5e-5, theta: 0, P: 1e7, fModel: 'colebrook', ...p }, vm = q.vsl + q.vsg;
  let r;
  if (vm <= 1e-9) { const liquid = q.vsg <= 1e-12; r = { holdup: liquid ? 1 : 0, fric: 0, grav: (liquid ? q.rhoL : q.rhoG) * G * Math.sin(q.theta), acc: 0, regime: 'static', tauW: 0 }; }
  else if (q.vsg <= 1e-9 * vm || q.vsl <= 1e-9 * vm) {
    const liquid = q.vsg <= 1e-9 * vm, rho = liquid ? q.rhoL : q.rhoG, mu = liquid ? q.muL : q.muG, f = frictionFactor((rho * vm * q.D) / mu, q.rough / q.D, q.fModel), fric = (f * rho * vm * vm) / (2 * q.D);
    const grav = rho * G * Math.sin(q.theta), Ek = liquid ? 0 : clamp((rho * vm * vm) / Math.max(q.P, 1e4), 0, 0.6); // gas expansion: dv/v = −dP/P (isothermal ideal-gas estimate)
    r = { holdup: liquid ? 1 : 0, fric, grav, acc: ((fric + grav) * Ek) / (1 - Ek), regime: liquid ? 'single-phase liquid' : 'single-phase gas', tauW: (fric * q.D) / 4 };
  } else if (model === 'homogeneous') {
    const lam = q.vsl / vm, rho = q.rhoL * lam + q.rhoG * (1 - lam), mu = q.muL * lam + q.muG * (1 - lam), f = frictionFactor((rho * vm * q.D) / mu, q.rough / q.D, q.fModel), fric = (f * rho * vm * vm) / (2 * q.D);
    const grav = rho * G * Math.sin(q.theta), Ek = clamp((rho * vm * q.vsg) / Math.max(q.P, 1e4), 0, 0.6);
    r = { holdup: lam, fric, grav, acc: ((fric + grav) * Ek) / (1 - Ek), regime: 'homogeneous', tauW: (fric * q.D) / 4 };
  } else r = model === 'driftFlux' ? driftFlux(q) : model === 'mechanistic' ? mechanistic(q) : beggsBrill(q);
  if (q.label !== false && (model === 'beggsBrill' || model === 'driftFlux')) { try { const fp = flowPattern(q).pattern; if (!fp.startsWith('single')) r.regime = fp; } catch { /* keep the correlation's own regime label */ } }
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
 *      pOut (bara) | pIn (bara), mScale (rate multiplier), model, fModel, n, idOf(s) (effective inner diameter), roughOf(s),
 *      energy: 'enthalpy' (default: flowing enthalpy + potential + kinetic energy balance, temperature from H(P, T) of the fluid model, so
 *      flashing, latent heat and the Joule–Thomson effect are included) | 'cpjt' (frozen heat capacities and Joule–Thomson coefficients) }
 * Returns { ok, s, x, z, theta, P, T, holdup, vsl, vsg, vm, rhoM, dpdx, regime, tauW, tAmb, tHyd, subcooling, qG, qL,
 *           pIn, pOut, tOut, dpFric, dpGrav, liquidInventory (m³), volume (m³), residence (s), heatLoss (W), mdot, energy (the form used) }.
 */
export function marchSteady(o) {
  const { fm, profile, id, rough = 4.5e-5, tIn = 70, mScale = 1, model = 'beggsBrill', fModel = 'colebrook', n = 200 } = o, grid = discretise(profile, n);
  const uOf = typeof o.uOf === 'function' ? o.uOf : () => o.U ?? 3, tAmbOf = typeof o.tAmbOf === 'function' ? o.tAmbOf : () => o.tAmb ?? 4;
  const dOf = typeof o.idOf === 'function' ? o.idOf : () => id, rOf = typeof o.roughOf === 'function' ? o.roughOf : () => rough;
  const mdot = (fm.rates.mHC + fm.rates.mW) * mScale, hFlow = (pr) => pr.mG * pr.hG + pr.mO * pr.hO + pr.mW * pr.hW;
  const state = (P, T, i, label) => {
    const sMid = 0.5 * (grid.s[i] + grid.s[i + 1]), D = Math.max(dOf(sMid), 0.01), A = (Math.PI * D * D) / 4, pr = fm.at(P, T, mScale), sinT = Math.sin(grid.theta[i]);
    const vsg = pr.qG / A, vsl = pr.qL / A, vm = vsl + vsg, gr = gradient({ vsl, vsg, rhoL: pr.rhoL, rhoG: pr.rhoG, muL: pr.muL, muG: pr.muG, sigma: pr.sigma, D, theta: grid.theta[i], rough: rOf(sMid), P: P * 1e5, fModel, label }, model);
    const mCp = pr.mG * pr.cpG + pr.mO * pr.cpO + pr.mW * pr.cpW, jt = mCp > 0 ? (pr.mG * pr.cpG * pr.jtG + pr.mO * pr.cpO * pr.jtO - (pr.mW / pr.rhoW)) / mCp : 0;
    const ta = tAmbOf(sMid, 0.5 * (grid.z[i] + grid.z[i + 1])), U = uOf(sMid), q = U * Math.PI * D * (T - ta); // W/m
    const dTds = mCp > 0 ? -q / mCp - jt * gr.dpdx - (mdot * G * sinT) / mCp : 0, H = gr.holdup;
    const ke = H > 1e-6 && H < 1 - 1e-6 ? 0.5 * (pr.mG * (vsg / (1 - H)) ** 2 + (pr.mO + pr.mW) * (vsl / H) ** 2) : 0.5 * mdot * vm * vm; // kinetic-energy flow (W)
    return { pr, gr, vsg, vsl, D, A, ta, q, mCp, sinT, ke, dPds: -gr.dpdx / 1e5, dTds };
  };
  const run = (pIn, label = false) => {
    const P = [pIn], T = [tIn], cells = [], pr0 = fm.at(pIn, tIn, mScale), enth = o.energy !== 'cpjt' && Number.isFinite(hFlow(pr0));
    let ok = true, Hf = enth ? hFlow(pr0) : 0;
    for (let i = 0; i < grid.n; i++) {
      const a = state(P[i], T[i], i, false), Pm = P[i] + 0.5 * grid.ds * a.dPds, Tm = T[i] + 0.5 * grid.ds * a.dTds;
      if (!(Pm > 1.0)) { ok = false; break; }
      const b = state(Pm, Tm, i, label), Pn = P[i] + grid.ds * b.dPds;
      let Tn = T[i] + grid.ds * b.dTds;
      if (!(Pn > 1.0) || !Number.isFinite(Tn)) { ok = false; break; }
      if (enth && b.mCp > 0) { // energy balance over the cell, then the temperature that carries the remaining enthalpy flow at the new pressure
        Hf -= grid.ds * (b.q + mdot * G * b.sinT) + 2 * (b.ke - a.ke);
        let t0 = Tn, f0 = hFlow(fm.at(Pn, t0, mScale)) - Hf, slope = b.mCp;
        for (let k = 0; k < 8 && Math.abs(f0) > 2e-5 * b.mCp; k++) { const t1 = clamp(t0 - f0 / slope, t0 - 25, t0 + 25), f1 = hFlow(fm.at(Pn, t1, mScale)) - Hf; if (Math.abs(t1 - t0) > 1e-9 && (f1 - f0) / (t1 - t0) > 0.2 * b.mCp) slope = (f1 - f0) / (t1 - t0); t0 = t1; f0 = f1; }
        Tn = t0;
      }
      P.push(Pn); T.push(clamp(Tn, -60, 250)); cells.push(b);
    }
    return { ok, P, T, cells, enth };
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
  sol = run(pIn, o.label !== false); // the flow-pattern label is evaluated once, on the converged march
  if (!sol.ok) return { ok: false, reason: 'The inlet pressure is too low to push this rate to the outlet.' };
  const N = grid.n, out = { ok: true, s: grid.s, x: grid.x, z: grid.z, theta: grid.theta.concat(grid.theta[N - 1]), P: sol.P, T: sol.T, mdot, pIn: sol.P[0], pOut: sol.P[N], tOut: sol.T[N], ds: grid.ds, length: grid.length, energy: sol.enth ? 'enthalpy' : 'cpjt' };
  const col = (fn) => { const a = sol.cells.map(fn); a.push(a[a.length - 1]); return a; };
  out.holdup = col((c) => c.gr.holdup); out.vsl = col((c) => c.vsl); out.vsg = col((c) => c.vsg); out.vm = col((c) => c.vsl + c.vsg); out.dpdx = col((c) => c.gr.dpdx); out.regime = col((c) => c.gr.regime); out.tauW = col((c) => c.gr.tauW);
  out.rhoM = col((c) => c.pr.rhoL * c.gr.holdup + c.pr.rhoG * (1 - c.gr.holdup)); out.tAmb = col((c) => c.ta); out.qG = col((c) => c.pr.qG); out.qL = col((c) => c.pr.qL); out.rhoL = col((c) => c.pr.rhoL); out.rhoG = col((c) => c.pr.rhoG); out.muL = col((c) => c.pr.muL); out.wcut = col((c) => c.pr.wcut);
  out.tHyd = out.P.map((p) => fm.hydrateT(p)); out.subcooling = out.T.map((t, i) => out.tHyd[i] - t);
  let inv = 0, vol = 0, fr = 0, gv = 0, heat = 0, res = 0;
  sol.cells.forEach((c) => { inv += c.gr.holdup * c.A * grid.ds; vol += c.A * grid.ds; fr += (c.gr.fric + c.gr.acc) * grid.ds; gv += c.gr.grav * grid.ds; heat += c.q * grid.ds; res += grid.ds / Math.max(c.vsl + c.vsg, 1e-6); });
  Object.assign(out, { liquidInventory: inv, volume: vol, dpFric: fr / 1e5, dpGrav: gv / 1e5, heatLoss: heat, residence: res });
  return out;
}
