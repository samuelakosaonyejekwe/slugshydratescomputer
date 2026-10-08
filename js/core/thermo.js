// Hydrocarbon thermodynamics kernel shared by every suite: component data, C7+ characterisation, cubic
// equations of state (Peng–Robinson, Soave–Redlich–Kwong, Redlich–Kwong, van der Waals) with volume
// translation, two-phase PT flash with tangent-plane stability, saturation points, phase envelope,
// thermal and transport properties, standard-condition rates, tabulated properties for the flow solvers
// and a screening hydrate curve with inhibitor/salt depression.
// Units at the interface: pressure bara, temperature °C, composition mol %, SI everywhere else.
import { brent, clamp, linspace, logspace } from './num.js';
import { density as rhoBrine, viscosity as muBrine, cp as cpBrine, conductivityThermal as kBrine, psat as psatWater } from './props.js';

export const R = 8.314462618;
export const P_STD = 1.01325, T_STD = 15; // bara, °C
export const VM_STD = (R * (T_STD + 273.15)) / (P_STD * 1e5); // m³/mol of ideal gas at standard conditions
export const MW_AIR = 28.9647;
const KEL = 273.15;

// id, name, Tc (K), Pc (bar), acentric, MW (g/mol), Vc (cm³/mol), parachor, ideal-gas Cp = a + bT + cT² + dT³ (J/mol/K), PR volume shift s = c/b
// Tc, Pc, acentric factor, molar mass and Vc (= 1/critical molar density): CoolProp fluid files (dev/fluids/*.json, the reference
// equations of state of each fluid), read 2026-10-08. Ideal-gas Cp: cubic least-squares fit (180–560 K, within 0.11 %) of the
// quartic Cp/R polynomials of Poling, Prausnitz & O'Connell (5th ed.) as tabulated in the open `chemicals` library
// (PolingDatabank.tsv). Parachors: Weinaug–Katz set (within 1.1 % of the open NeqSim component database). 'Hexanes' is n-hexane.
const C = (id, name, Tc, Pc, w, MW, Vc, par, cp, s, hyd) => ({ id, name, Tc, Pc, w, MW, Vc, par, cp, s, hyd });
export const COMPONENTS = Object.freeze({
  N2: C('N2', 'Nitrogen', 126.192, 33.958, 0.0372, 28.01348, 89.41, 41, [29.544, -3.6408e-3, 7.0708e-6, 8.8174e-10], -0.1927),
  CO2: C('CO2', 'Carbon dioxide', 304.1282, 73.773, 0.22394, 44.0098, 94.12, 78, [25.891, 2.6372e-2, 5.7351e-5, -6.9062e-8], -0.0817),
  H2S: C('H2S', 'Hydrogen sulphide', 373.1, 90, 0.1005, 34.08088, 98.14, 80, [34.893, -2.1429e-2, 7.7953e-5, -5.0913e-8], -0.1288),
  C1: C('C1', 'Methane', 190.564, 45.992, 0.01142, 16.0428, 98.63, 77, [36.743, -5.9109e-2, 2.3241e-4, -1.5099e-7], -0.1595),
  C2: C('C2', 'Ethane', 305.322, 48.722, 0.099, 30.06904, 145.84, 108, [32.051, -2.6839e-3, 3.1575e-4, -2.5465e-7], -0.1134),
  C3: C('C3', 'Propane', 369.89, 42.512, 0.1521, 44.09562, 200, 150.3, [28.711, 8.4428e-2, 3.0944e-4, -2.8812e-7], -0.0863),
  iC4: C('iC4', 'i-Butane', 407.817, 36.29, 0.18353, 58.1222, 257.75, 181.5, [24.439, 1.9245e-1, 2.5556e-4, -2.8631e-7], -0.0844),
  nC4: C('nC4', 'n-Butane', 425.125, 37.96, 0.20081, 58.1222, 254.92, 189.9, [41.701, 1.0232e-1, 4.1372e-4, -3.8408e-7], -0.0675),
  iC5: C('iC5', 'i-Pentane', 460.35, 33.78, 0.2274, 72.14878, 305.72, 225, [14.016, 3.4664e-1, 6.9297e-5, -1.7209e-7], -0.0608),
  nC5: C('nC5', 'n-Pentane', 469.7, 33.675, 0.25103, 72.14878, 310.99, 231.5, [56.664, 7.5215e-2, 6.2863e-4, -5.5366e-7], -0.039),
  C6: C('C6', 'Hexanes', 507.82, 30.441, 0.30032, 86.17536, 369.55, 271, [65.823, 9.5497e-2, 7.4806e-4, -6.7033e-7], -0.008),
});
export const LIGHT_IDS = Object.freeze(Object.keys(COMPONENTS));
export const COMP_IDS = Object.freeze([...LIGHT_IDS, 'C7p']);
export const COMP_LABELS = Object.freeze({ ...Object.fromEntries(LIGHT_IDS.map((k) => [k, COMPONENTS[k].name])), C7p: 'Heptanes plus (C7+)' });
export const INHIBITORS = Object.freeze({
  none: { name: 'None', MW: 18.015, rho: 1000, K: 0 },
  MeOH: { name: 'Methanol', MW: 32.042, rho: 792, K: 1297 },
  MEG: { name: 'Mono-ethylene glycol (MEG)', MW: 62.068, rho: 1113, K: 1297 },
  DEG: { name: 'Di-ethylene glycol (DEG)', MW: 106.12, rho: 1118, K: 2222 },
  TEG: { name: 'Tri-ethylene glycol (TEG)', MW: 150.17, rho: 1125, K: 2222 },
  EtOH: { name: 'Ethanol', MW: 46.069, rho: 789, K: 1297 },
});

export const EOS = Object.freeze({
  PR: { name: 'Peng–Robinson (1978)', d1: 1 + Math.SQRT2, d2: 1 - Math.SQRT2, oa: 0.457235529, ob: 0.0777960739, m: (w) => (w <= 0.491 ? 0.37464 + 1.54226 * w - 0.26992 * w * w : 0.379642 + 1.48503 * w - 0.164423 * w * w + 0.016666 * w ** 3), alpha: (m, Tr) => (1 + m * (1 - Math.sqrt(Tr))) ** 2 },
  SRK: { name: 'Soave–Redlich–Kwong', d1: 1, d2: 0, oa: 0.42748023, ob: 0.08664035, m: (w) => 0.48 + 1.574 * w - 0.176 * w * w, alpha: (m, Tr) => (1 + m * (1 - Math.sqrt(Tr))) ** 2 },
  RK: { name: 'Redlich–Kwong', d1: 1, d2: 0, oa: 0.42748023, ob: 0.08664035, m: () => 0, alpha: (m, Tr) => 1 / Math.sqrt(Tr) },
  vdW: { name: 'van der Waals', d1: 0, d2: 0, oa: 27 / 64, ob: 1 / 8, m: () => 0, alpha: () => 1 },
});

// ---- C7+ characterisation ---------------------------------------------------------------------------
const LAGUERRE = { 1: [[1], [1]], 2: [[0.585786, 3.414214], [0.853553, 0.146447]], 3: [[0.415775, 2.29428, 6.289945], [0.711093, 0.278518, 0.0103893]] };
/** Søreide boiling point (K) from molar mass and specific gravity. */
export const tbSoreide = (M, SG) => (1928.3 - 1.695e5 * M ** -0.03522 * SG ** 3.266 * Math.exp(-4.922e-3 * M - 4.7685 * SG + 3.462e-3 * M * SG)) / 1.8;
/** Kesler–Lee critical properties and Lee–Kesler / Kesler–Lee acentric factor of a petroleum fraction. */
export function pseudoProps(M, SG) {
  const Tb = tbSoreide(M, SG) * 1.8; // °R
  const Tc = 341.7 + 811 * SG + (0.4244 + 0.1174 * SG) * Tb + ((0.4669 - 3.2623 * SG) * 1e5) / Tb;
  const lnPc = 8.3634 - 0.0566 / SG - (0.24244 + 2.2898 / SG + 0.11857 / SG ** 2) * 1e-3 * Tb + (1.4685 + 3.648 / SG + 0.47227 / SG ** 2) * 1e-7 * Tb ** 2 - (0.42019 + 1.6977 / SG ** 2) * 1e-10 * Tb ** 3;
  const Pc = Math.exp(lnPc), Tbr = Tb / Tc, Kw = Tb ** (1 / 3) / SG;
  const w = Tbr < 0.8
    ? (-Math.log(Pc / 14.696) - 5.92714 + 6.09648 / Tbr + 1.28862 * Math.log(Tbr) - 0.169347 * Tbr ** 6) / (15.2518 - 15.6875 / Tbr - 13.4721 * Math.log(Tbr) + 0.43577 * Tbr ** 6)
    : -7.904 + 0.1352 * Kw - 0.007465 * Kw * Kw + 8.359 * Tbr + (1.408 - 0.01063 * Kw) / Tbr;
  const Vc = (21.573 + 0.015122 * M - 27.656 * SG + 0.070615 * M * SG) * 62.428; // ft³/lbmol -> cm³/mol
  return { Tb: Tb / 1.8, Tc: clamp(Tc / 1.8, 500, 1150), Pc: clamp(Pc * 0.0689476, 5, 40), w: clamp(w, 0.2, 1.7), Vc, Kw };
}
function splitC7(z7, M7, SG7, n) {
  n = clamp(Math.round(n) || 3, 1, 3);
  const eta = 90, [xs, ws] = LAGUERRE[n], Ms = xs.map((x) => (n === 1 ? M7 : eta + Math.max(5, M7 - eta) * x));
  // Søreide specific gravities with the factor chosen so that the mixture reproduces the measured C7+ gravity
  const sgOf = (Cf) => Ms.map((M) => 0.2855 + Cf * Math.max(1, M - 66) ** 0.13), mixSG = (Cf) => { const s = sgOf(Cf); let a = 0, b = 0; ws.forEach((w, i) => { a += w * Ms[i]; b += (w * Ms[i]) / s[i]; }); return a / b; };
  let SGs = [SG7];
  if (n > 1) { let Cf = 0.29; try { Cf = brent((c) => mixSG(c) - SG7, 0.2, 0.42, 1e-10); } catch { /* keep default */ } SGs = sgOf(Cf).map((s) => clamp(s, 0.6, 1.2)); }
  return Ms.map((M, i) => {
    const SG = SGs[i], p = pseudoProps(M, SG);
    return { id: n === 1 ? 'C7+' : `C7+_${i + 1}`, name: n === 1 ? 'C7+' : `C7+ pseudo ${i + 1}`, pseudo: true, z: z7 * ws[i], MW: M, SG, Tb: p.Tb, Tc: p.Tc, Pc: p.Pc, w: p.w, Vc: p.Vc, par: 59.3 + 2.34 * M, cp: [0.39 * M, 0.0042 * M, 0, 0], s: clamp(1 - 2.258 / M ** 0.1823, -0.05, 0.35) };
  });
}
// Peng–Robinson binary interaction parameters: ChemSep table distributed with the open `thermo` library (Interaction Parameters/
// ChemSep/pr.json, regressed to the DECHEMA vapour–liquid data collection), read 2026-10-08. H2S–methane is absent from that table
// (0.07 kept); the 'X-*' entries are the defaults for pairs it does not list and for the C7+ pseudo-components.
const KIJ = {
  'N2-CO2': -0.0122, 'N2-H2S': 0.1652, 'N2-C1': 0.0289, 'N2-C2': 0.0533, 'N2-C3': 0.0878, 'N2-iC4': 0.1033, 'N2-nC4': 0.0711, 'N2-iC5': 0.0922, 'N2-nC5': 0.1, 'N2-C6': 0.1496, 'N2-*': 0.1,
  'CO2-H2S': 0.0967, 'CO2-C1': 0.0978, 'CO2-C2': 0.13, 'CO2-C3': 0.1315, 'CO2-iC4': 0.13, 'CO2-nC4': 0.1352, 'CO2-iC5': 0.1219, 'CO2-nC5': 0.1252, 'CO2-C6': 0.11, 'CO2-*': 0.115,
  'H2S-C1': 0.07, 'H2S-C2': 0.0952, 'H2S-C3': 0.0878, 'H2S-iC4': 0.0474, 'H2S-nC5': 0.063, 'H2S-*': 0.06,
  'C1-C2': -0.0059, 'C1-C3': 0.0119, 'C1-iC4': 0.0256, 'C1-nC4': 0.0185, 'C1-iC5': -0.0056, 'C1-nC5': 0.023, 'C1-C6': 0.04,
  'C2-C3': 0.0011, 'C2-iC4': -0.0067, 'C2-nC4': 0.0089, 'C2-nC5': 0.0078, 'C2-C6': -0.04, 'C3-iC4': -0.0078, 'C3-nC4': 0.0033, 'C3-iC5': 0.0111, 'C3-nC5': 0.0267, 'C3-C6': 0.0007,
  'iC4-nC4': -0.0004, 'nC4-nC5': 0.0174, 'nC4-C6': -0.0056,
};
function kijOf(a, b) {
  if (a.id === b.id) return 0;
  const k = KIJ[`${a.id}-${b.id}`] ?? KIJ[`${b.id}-${a.id}`];
  if (k !== undefined) return k;
  for (const [p, q] of [[a, b], [b, a]]) if (['N2', 'CO2', 'H2S'].includes(p.id)) return KIJ[`${p.id}-*`];
  for (const [p, q] of [[a, b], [b, a]]) if (p.id === 'C1' && q.pseudo) return clamp(0.14 * q.SG - 0.0668, 0.02, 0.09);
  return 0;
}

export const DEFAULT_FLUID = Object.freeze({
  name: 'Light oil with associated gas', comp: Object.freeze({ N2: 0.5, CO2: 2, H2S: 0, C1: 45, C2: 7, C3: 5, iC4: 1, nC4: 2.5, iC5: 1, nC5: 1.5, C6: 2.5, C7p: 32 }), c7MW: 210, c7SG: 0.84,
  rateBasis: 'oil', qOil: 3000, qGas: 1, mdot: 30, wc: 20, qWater: 0, salinity: 3.5, inhibitor: 'none', inhWt: 0, eos: 'PR', nPseudo: 3, Tin: 70, Pout: 25, Tres: 90, Pres: 300,
});

/**
 * Build the EOS model of a hydrocarbon mixture from a case fluid specification.
 * spec: { comp: { N2..C6, C7p } in mol %, c7MW, c7SG }; opts: { eos, nPseudo, kijScale, vcMult, tcMult, pcMult, wMult, shift }.
 */
export function makeFluid(spec = DEFAULT_FLUID, opts = {}) {
  const o = { eos: spec.eos || 'PR', nPseudo: spec.nPseudo ?? 3, kijScale: 1, vcMult: 1, tcMult: 1, pcMult: 1, wMult: 1, shift: true, ...opts };
  const eos = EOS[o.eos] || EOS.PR, comp = spec.comp || DEFAULT_FLUID.comp;
  let list = LIGHT_IDS.filter((k) => +comp[k] > 0).map((k) => ({ ...COMPONENTS[k], z: +comp[k] }));
  if (+comp.C7p > 0) list.push(...splitC7(+comp.C7p, clamp(+spec.c7MW || 210, 96, 600), clamp(+spec.c7SG || 0.84, 0.7, 1.05), o.nPseudo).map((p) => ({ ...p, Tc: p.Tc * o.tcMult, Pc: p.Pc * o.pcMult, w: p.w * o.wMult, Vc: p.Vc * o.vcMult })));
  if (!list.length) list = [{ ...COMPONENTS.C1, z: 100 }];
  const tot = list.reduce((s, c) => s + c.z, 0);
  const comps = list.map((c) => ({ ...c, z: c.z / tot, ac: (eos.oa * (R * c.Tc) ** 2) / (c.Pc * 1e5), b: (eos.ob * R * c.Tc) / (c.Pc * 1e5), m: eos.m(c.w), c: o.shift && o.eos === 'PR' ? c.s * ((eos.ob * R * c.Tc) / (c.Pc * 1e5)) : o.shift && o.eos === 'SRK' ? 0.40768 * ((R * c.Tc) / (c.Pc * 1e5)) * (0.29441 - (0.29056 - 0.08775 * c.w)) : 0 }));
  const n = comps.length, kij = comps.map((a) => comps.map((b) => o.kijScale * kijOf(a, b)));
  return { comps, n, z: comps.map((c) => c.z), kij, eos, eosId: EOS[o.eos] ? o.eos : 'PR', opts: o, MW: comps.reduce((s, c) => s + c.z * c.MW, 0) };
}

// ---- cubic EOS --------------------------------------------------------------------------------------
function cubicRoots(c2, c1, c0) { // real roots of Z³ + c2 Z² + c1 Z + c0
  const q = (3 * c1 - c2 * c2) / 9, r = (9 * c2 * c1 - 27 * c0 - 2 * c2 ** 3) / 54, d = q ** 3 + r * r;
  if (d >= 0) { const s = Math.cbrt(r + Math.sqrt(d)), t = Math.cbrt(r - Math.sqrt(d)); return [s + t - c2 / 3]; }
  const th = Math.acos(clamp(r / Math.sqrt(-(q ** 3)), -1, 1)), k = 2 * Math.sqrt(-q);
  return [k * Math.cos(th / 3) - c2 / 3, k * Math.cos((th + 2 * Math.PI) / 3) - c2 / 3, k * Math.cos((th + 4 * Math.PI) / 3) - c2 / 3];
}
/** Mixture EOS terms at T (K) for composition x. */
function mixTerms(f, x, T) {
  const n = f.n, ai = new Array(n), sa = new Array(n).fill(0);
  for (let i = 0; i < n; i++) ai[i] = Math.sqrt(f.comps[i].ac * f.eos.alpha(f.comps[i].m, T / f.comps[i].Tc));
  let a = 0, b = 0;
  for (let i = 0; i < n; i++) { let s = 0; for (let j = 0; j < n; j++) s += x[j] * ai[j] * (1 - f.kij[i][j]); sa[i] = ai[i] * s; a += x[i] * sa[i]; b += x[i] * f.comps[i].b; }
  return { a, b, sa };
}
/**
 * Compressibility factor and fugacity coefficients of one phase. kind: 'liquid' | 'vapour' | 'stable' (lowest Gibbs energy).
 * Returns { Z, lnphi[], a, b, A, B, g } where g is the reduced molar Gibbs energy of the phase.
 */
export function eosPhase(f, x, Pbar, TK, kind = 'stable') {
  const P = Pbar * 1e5, { a, b, sa } = mixTerms(f, x, TK), A = (a * P) / (R * TK) ** 2, B = (b * P) / (R * TK), { d1, d2 } = f.eos, u = d1 + d2, w = d1 * d2;
  const roots = cubicRoots(-(1 + B - u * B), A + w * B * B - u * B - u * B * B, -(A * B + w * B * B + w * B ** 3)).filter((z) => z > B * (1 + 1e-12));
  const evalZ = (Z) => {
    const L = Math.abs(d1 - d2) > 1e-12 ? Math.log((Z + d1 * B) / (Z + d2 * B)) / ((d1 - d2) * B) : 1 / (Z + d1 * B);
    const lnphi = new Array(f.n);
    let g = 0;
    for (let i = 0; i < f.n; i++) { const bi = f.comps[i].b / b; lnphi[i] = bi * (Z - 1) - Math.log(Z - B) - A * L * ((2 * sa[i]) / a - bi); if (x[i] > 0) g += x[i] * (Math.log(x[i]) + lnphi[i]); }
    return { Z, lnphi, a, b, A, B, g, L };
  };
  if (!roots.length) return evalZ(Math.max(B * 1.0001, 1));
  if (roots.length === 1 || kind === 'liquid' || kind === 'vapour') return evalZ(kind === 'liquid' ? Math.min(...roots) : kind === 'vapour' ? Math.max(...roots) : roots[0]);
  const lo = evalZ(Math.min(...roots)), hi = evalZ(Math.max(...roots));
  return lo.g < hi.g ? lo : hi;
}

const wilsonK = (f, Pbar, TK) => f.comps.map((c) => (c.Pc / Pbar) * Math.exp(5.373 * (1 + c.w) * (1 - c.Tc / TK)));
/** Rachford–Rice vapour fraction for K-values (bracketed, always converges when a root exists). */
export function rachfordRice(z, K) {
  const g = (b) => { let s = 0; for (let i = 0; i < z.length; i++) s += (z[i] * (K[i] - 1)) / (1 + b * (K[i] - 1)); return s; };
  const g0 = g(0), g1 = g(1);
  if (g0 <= 0) return 0;
  if (g1 >= 0) return 1;
  return brent(g, 0, 1, 1e-13, 100);
}
/** Michelsen tangent-plane stability test of a feed at P, T. Returns { stable, K } (K-values for flash initiation when unstable). */
export function stability(f, Pbar, TK) {
  const z = f.z, base = eosPhase(f, z, Pbar, TK), d = z.map((zi, i) => Math.log(Math.max(zi, 1e-300)) + base.lnphi[i]), Kw = wilsonK(f, Pbar, TK);
  let best = { tpd: 0, K: null };
  for (const trial of ['vapour', 'liquid']) {
    let W = z.map((zi, i) => (trial === 'vapour' ? zi * Kw[i] : zi / Kw[i])), tpd = 0, trivial = false;
    for (let it = 0; it < 300; it++) {
      const S = W.reduce((s, v) => s + v, 0), x = W.map((v) => v / S), ph = eosPhase(f, x, Pbar, TK, trial);
      let err = 0;
      const Wn = W.map((v, i) => { const nv = Math.exp(d[i] - ph.lnphi[i]); err += (Math.log(Math.max(nv, 1e-300) / Math.max(v, 1e-300))) ** 2; return nv; });
      W = Wn;
      if (err < 1e-13) break;
      let tr = 0; for (let i = 0; i < f.n; i++) tr += (Math.log(Math.max(W[i], 1e-300) / Math.max(z[i], 1e-300))) ** 2;
      if (tr < 1e-8 && it > 2) { trivial = true; break; }
    }
    const S = W.reduce((s, v) => s + v, 0);
    tpd = 1 - S; // modified tangent-plane distance at the stationary point: negative means unstable
    if (!trivial && tpd < best.tpd - 1e-10) best = { tpd, K: W.map((v, i) => (trial === 'vapour' ? v / S / z[i] : z[i] / (v / S))) };
  }
  return { stable: !best.K, K: best.K, tpd: best.tpd };
}
/**
 * Two-phase PT flash. Returns { phase: 'gas' | 'oil' | 'two', beta (vapour mole fraction), x[], y[], K[], liq, vap (eosPhase results), iterations }.
 */
export function flashPT(f, Pbar, Tc) {
  const TK = Tc + KEL, z = f.z, n = f.n;
  const single = () => {
    const ph = eosPhase(f, z, Pbar, TK), vm = (ph.Z * R * TK) / (Pbar * 1e5), gas = vm / ph.b > 1.75;
    return { phase: gas ? 'gas' : 'oil', beta: gas ? 1 : 0, x: z.slice(), y: z.slice(), K: z.map(() => 1), liq: ph, vap: ph, iterations: 0 };
  };
  if (n === 1) return single();
  const attempt = (K0) => {
    let K = K0.slice(), beta = 0.5, x = z, y = z, liq = null, vap = null, it = 0;
    for (; it < 400; it++) {
      beta = rachfordRice(z, K);
      if (beta <= 0 || beta >= 1) { // outside the two-phase window of these K-values: nudge inside once, else give up
        const sK = z.reduce((s, zi, i) => s + zi * K[i], 0), sInv = z.reduce((s, zi, i) => s + zi / K[i], 0);
        if (sK <= 1 || sInv <= 1) return null;
      }
      x = z.map((zi, i) => zi / (1 + beta * (K[i] - 1))); y = x.map((xi, i) => xi * K[i]);
      const sx = x.reduce((s, v) => s + v, 0), sy = y.reduce((s, v) => s + v, 0);
      x = x.map((v) => v / sx); y = y.map((v) => v / sy);
      liq = eosPhase(f, x, Pbar, TK, 'liquid'); vap = eosPhase(f, y, Pbar, TK, 'vapour');
      let err = 0;
      for (let i = 0; i < n; i++) { const Kn = Math.exp(liq.lnphi[i] - vap.lnphi[i]); err += (Kn / K[i] - 1) ** 2; K[i] = Kn; }
      if (err < 1e-14) break;
      let tr = 0; for (let i = 0; i < n; i++) tr += Math.log(K[i]) ** 2;
      if (tr < 1e-8) return null; // collapsed to the trivial solution
    }
    if (!(beta > 1e-12 && beta < 1 - 1e-12)) return null;
    return { phase: 'two', beta, x, y, K, liq, vap, iterations: it };
  };
  let res = attempt(wilsonK(f, Pbar, TK));
  if (!res) { const st = stability(f, Pbar, TK); if (!st.stable) res = attempt(st.K); }
  return res || single();
}

// ---- properties ---------------------------------------------------------------------------------------
const cpIdeal = (f, x, TK) => { let s = 0; for (let i = 0; i < f.n; i++) { const c = f.comps[i].cp; s += x[i] * (c[0] + c[1] * TK + c[2] * TK * TK + c[3] * TK ** 3); } return s; };
const hIdeal = (f, x, TK) => { let s = 0; const T0 = KEL; for (let i = 0; i < f.n; i++) { const c = f.comps[i].cp; s += x[i] * (c[0] * (TK - T0) + (c[1] / 2) * (TK ** 2 - T0 ** 2) + (c[2] / 3) * (TK ** 3 - T0 ** 3) + (c[3] / 4) * (TK ** 4 - T0 ** 4)); } return s; };
function hResidual(f, x, Pbar, TK, kind) { // J/mol
  const ph = eosPhase(f, x, Pbar, TK, kind), dT = 0.05, a1 = mixTerms(f, x, TK + dT).a, a0 = mixTerms(f, x, TK - dT).a, da = (a1 - a0) / (2 * dT);
  return R * TK * (ph.Z - 1) + ((TK * da - ph.a) * ph.L * (Pbar * 1e5)) / (R * TK);
}
/** Molar enthalpy (J/mol) relative to the ideal gas at 0 °C. */
export const enthalpyMolar = (f, x, Pbar, TK, kind) => hIdeal(f, x, TK) + hResidual(f, x, Pbar, TK, kind);

function viscosityLBC(f, x, TK, rhoMolar) { // rhoMolar mol/m³ -> Pa·s
  let num = 0, den = 0, Tc = 0, Pc = 0, M = 0, Vc = 0;
  for (let i = 0; i < f.n; i++) {
    const c = f.comps[i], Tr = TK / c.Tc, xi = c.Tc ** (1 / 6) / (Math.sqrt(c.MW) * (c.Pc / 1.01325) ** (2 / 3));
    const mu = (Tr < 1.5 ? 34e-5 * Tr ** 0.94 : 17.78e-5 * (4.58 * Tr - 1.67) ** 0.625) / xi;
    num += x[i] * mu * Math.sqrt(c.MW); den += x[i] * Math.sqrt(c.MW); Tc += x[i] * c.Tc; Pc += x[i] * c.Pc; M += x[i] * c.MW; Vc += x[i] * c.Vc;
  }
  const mu0 = num / den, xi = Tc ** (1 / 6) / (Math.sqrt(M) * (Pc / 1.01325) ** (2 / 3)), rr = clamp(rhoMolar * Vc * 1e-6, 0, 3.4);
  const poly = 0.1023 + 0.023364 * rr + 0.058533 * rr * rr - 0.040758 * rr ** 3 + 0.0093324 * rr ** 4;
  return Math.max(mu0, mu0 + (poly ** 4 - 1e-4) / xi) * 1e-3;
}
/**
 * Properties of one hydrocarbon phase of composition x at P (bara), T (°C). kind 'liquid' | 'vapour'.
 * Returns { Z, MW (g/mol), vm (m³/mol), rho (kg/m³), mu (Pa·s), h (J/kg), cp (J/kg/K), k (W/m/K), jt (K/Pa) }.
 */
export function phaseProps(f, x, Pbar, Tc, kind, { thermal = true } = {}) {
  const TK = Tc + KEL, ph = eosPhase(f, x, Pbar, TK, kind);
  let MW = 0, cs = 0; for (let i = 0; i < f.n; i++) { MW += x[i] * f.comps[i].MW; cs += x[i] * f.comps[i].c; }
  const vm = Math.max((ph.Z * R * TK) / (Pbar * 1e5) - cs, 0.6 * ph.b), rho = (MW * 1e-3) / vm, mu = viscosityLBC(f, x, TK, 1 / vm);
  const out = { Z: (Pbar * 1e5 * vm) / (R * TK), MW, vm, rho, mu };
  if (thermal) {
    const dT = 0.5, dP = Math.max(0.02, Pbar * 1e-3), h = enthalpyMolar(f, x, Pbar, TK, kind), hT = enthalpyMolar(f, x, Pbar, TK + dT, kind), hP = enthalpyMolar(f, x, Pbar + dP, TK, kind);
    const cpm = Math.max((hT - h) / dT, 0.6 * cpIdeal(f, x, TK));
    out.h = h / (MW * 1e-3); out.cp = cpm / (MW * 1e-3); out.jt = -((hP - h) / (dP * 1e5)) / cpm;
    out.k = kind === 'vapour' || rho < 250 ? mu * (out.cp + (1.25 * R) / (MW * 1e-3)) : clamp((0.1172 * (1 - 0.00054 * Tc)) / clamp(rho / 1000, 0.5, 1.1), 0.07, 0.2);
  }
  return out;
}
/** Gas–oil interfacial tension (N/m) by the Macleod–Sugden parachor method. */
export function interfacialTension(f, x, y, vmL, vmV) {
  let s = 0; for (let i = 0; i < f.n; i++) s += f.comps[i].par * (x[i] / (vmL * 1e6) - y[i] / (vmV * 1e6));
  return clamp(Math.max(s, 0) ** 4 * 1e-3, 1e-5, 0.04);
}
/**
 * Full hydrocarbon state at P (bara), T (°C): flash plus phase properties.
 * Returns { phase, beta, wG (gas mass fraction), gas, oil, sigma (N/m), x, y } — gas/oil are phaseProps objects (the absent phase mirrors the present one).
 */
export function props(f, Pbar, Tc, opt) {
  const fl = flashPT(f, Pbar, Tc);
  if (fl.phase !== 'two') {
    const kind = fl.phase === 'gas' ? 'vapour' : 'liquid', p = phaseProps(f, f.z, Pbar, Tc, kind, opt);
    return { phase: fl.phase, beta: fl.beta, wG: fl.beta, gas: p, oil: p, sigma: fl.phase === 'gas' ? 0.02 : 0.005, x: fl.x, y: fl.y };
  }
  const gas = phaseProps(f, fl.y, Pbar, Tc, 'vapour', opt), oil = phaseProps(f, fl.x, Pbar, Tc, 'liquid', opt);
  return { phase: 'two', beta: fl.beta, wG: (fl.beta * gas.MW) / (fl.beta * gas.MW + (1 - fl.beta) * oil.MW), gas, oil, sigma: interfacialTension(f, fl.x, fl.y, oil.vm, gas.vm), x: fl.x, y: fl.y };
}

// ---- saturation and phase envelope ------------------------------------------------------------------
/** Upper saturation pressure (bara) at T (°C): { P, type: 'bubble' | 'dew' } or { P: null } when single-phase at every pressure. */
export function saturationP(f, Tc, { Pmax = 1500 } = {}) {
  const two = (P) => flashPT(f, P, Tc).phase === 'two', grid = logspace(1, Pmax, 46);
  let hi = -1;
  for (let i = grid.length - 1; i >= 0; i--) if (two(grid[i])) { hi = i; break; }
  if (hi < 0) return { P: null, type: null };
  if (hi === grid.length - 1) return { P: Pmax, type: 'bubble', capped: true };
  let a = grid[hi], b = grid[hi + 1];
  for (let k = 0; k < 44 && b - a > 1e-4 * a; k++) { const m = 0.5 * (a + b); if (two(m)) a = m; else b = m; }
  const beta = flashPT(f, a, Tc).beta;
  return { P: 0.5 * (a + b), type: beta < 0.5 ? 'bubble' : 'dew', beta };
}
/** Lower dew-point pressure (bara) at T (°C), or null. */
export function lowerDewP(f, Tc) {
  const two = (P) => flashPT(f, P, Tc).phase === 'two', grid = logspace(0.01, 1000, 40);
  let lo = -1;
  for (let i = 0; i < grid.length; i++) if (two(grid[i])) { lo = i; break; }
  if (lo <= 0) return null;
  let a = grid[lo - 1], b = grid[lo];
  for (let k = 0; k < 40 && b - a > 1e-4 * b; k++) { const m = 0.5 * (a + b); if (two(m)) b = m; else a = m; }
  return 0.5 * (a + b);
}
/** Phase envelope: saturation pressure against temperature with the cricondenbar, cricondentherm and an estimate of the critical point. */
export function phaseEnvelope(f, { Tmin = -100, Tmax = 600, n = 48 } = {}) {
  const T = [], P = [], type = [], lowT = [], lowP = [];
  for (const t of linspace(Tmin, Tmax, n)) {
    const s = saturationP(f, t);
    if (s.P === null) { if (T.length) break; continue; }
    T.push(t); P.push(s.P); type.push(s.type);
    const l = lowerDewP(f, t); if (l) { lowT.push(t); lowP.push(l); }
  }
  let crit = null;
  for (let i = 1; i < T.length; i++) if (type[i] !== type[i - 1]) { crit = { T: 0.5 * (T[i] + T[i - 1]), P: 0.5 * (P[i] + P[i - 1]) }; break; }
  const iB = P.indexOf(Math.max(...P));
  return { T, P, type, lowT, lowP, critical: crit, cricondenbar: T.length ? { T: T[iB], P: P[iB] } : null, cricondentherm: T.length ? { T: T[T.length - 1], P: P[P.length - 1] } : null };
}

// ---- standard conditions and stream rates -------------------------------------------------------------
/** Single-stage flash to standard conditions of one mole of feed. */
export function stdFlash(f) {
  const p = props(f, P_STD, T_STD, { thermal: false });
  const vOil = p.phase === 'gas' ? 0 : (1 - p.beta) * p.oil.vm, vGas = p.beta * VM_STD;
  const rhoOil = p.phase === 'gas' ? 0 : p.oil.rho, sgOil = rhoOil / 999.016;
  return { beta: p.beta, vOil, vGas, gor: vOil > 0 ? vGas / vOil : Infinity, rhoOilStd: rhoOil, api: sgOil > 0 ? 141.5 / sgOil - 131.5 : null, mwGas: p.beta > 0 ? p.gas.MW : null, gasSG: p.beta > 0 ? p.gas.MW / MW_AIR : null, mwOil: p.beta < 1 ? p.oil.MW : null, x: p.x, y: p.y };
}
/** Aqueous-phase description of a case fluid: salinity (g/kg), density at standard conditions, inhibitor data. */
export function aqueous(spec = DEFAULT_FLUID) {
  const inh = INHIBITORS[spec.inhibitor] || INHIBITORS.none, w = spec.inhibitor && spec.inhibitor !== 'none' ? clamp(+spec.inhWt || 0, 0, 95) : 0;
  return { S: clamp((+spec.salinity || 0) * 10, 0, 260), inh, inhWt: w, inhId: w > 0 ? spec.inhibitor : 'none' };
}
/**
 * Mass and molar rates of the hydrocarbon and water streams from the rate basis of the case.
 * Returns { nHC (mol/s), mHC (kg/s), mW (kg/s), qOilStd, qGasStd, qWaterStd (Sm³/d), gor (Sm³/Sm³), wc (fraction of std liquid), std }.
 */
export function streams(spec = DEFAULT_FLUID, f = makeFluid(spec)) {
  const std = stdFlash(f), day = 86400, basis = spec.rateBasis || 'oil';
  let nHC;
  if (basis === 'gas' || (basis === 'oil' && std.vOil <= 0)) nHC = ((+spec.qGas || 0) * 1e6) / day / Math.max(std.vGas, 1e-12);
  else if (basis === 'mass') nHC = (+spec.mdot || 0) / (f.MW * 1e-3);
  else nHC = (+spec.qOil || 0) / day / std.vOil;
  const qOilStd = nHC * std.vOil * day, qGasStd = nHC * std.vGas * day;
  const wc = clamp((+spec.wc || 0) / 100, 0, 0.98), qWaterStd = basis === 'oil' && std.vOil > 0 ? (qOilStd * wc) / (1 - wc) : +spec.qWater || 0;
  const aq = aqueous(spec), rhoW = rhoBrine(T_STD, aq.S);
  return { nHC, mHC: nHC * f.MW * 1e-3, mW: (qWaterStd * rhoW) / day, qOilStd, qGasStd, qWaterStd, gor: std.gor, wc: qOilStd + qWaterStd > 0 ? qWaterStd / (qOilStd + qWaterStd) : 0, std, rhoWaterStd: rhoW };
}

// ---- water and aqueous phase ---------------------------------------------------------------------------
/** Aqueous-phase properties at P (bara), T (°C) for salinity S (g/kg) and an inhibitor mass fraction (wt %). */
export function waterProps(Pbar, Tc, aq = { S: 35, inhWt: 0, inh: INHIBITORS.none }) {
  const T = clamp(Tc, -5, 180), w = (aq.inhWt || 0) / 100;
  let rho = rhoBrine(T, aq.S) * (1 + 4.6e-5 * (Pbar - 1)), mu = muBrine(clamp(T, 0, 180), aq.S), cp = cpBrine(T, aq.S), k = kBrine(T, aq.S);
  if (w > 0) { // ideal-volume mixing for density; Arrhenius (log) mixing for viscosity with the inhibitor's own temperature trend
    const rI = aq.inh.rho * (1 - 7e-4 * (T - 20)), muI = aq.inh.rho > 1000 ? 0.0161 * Math.exp(-0.035 * (T - 25)) * (aq.inh.MW / 62.07) ** 1.2 : 0.00054 * Math.exp(-0.012 * (T - 25));
    rho = 1 / (w / rI + (1 - w) / rho); mu = Math.exp(w * Math.log(muI) + (1 - w) * Math.log(mu)); cp = w * (aq.inh.rho > 1000 ? 2400 : 2500) + (1 - w) * cp; k = w * (aq.inh.rho > 1000 ? 0.25 : 0.2) + (1 - w) * k;
  }
  return { rho, mu, cp, k, h: cp * Tc, jt: -1 / (rho * cp) };
}
/** Equilibrium water content of gas (kg water per Sm³ gas), Bukacek correlation. */
export function waterContent(Pbar, Tc) {
  const TR = (Tc + KEL) * 1.8, Ppsia = Pbar * 14.5038, pv = (psatWater(clamp(Tc, 0.01, 370)) / 1e5) * 14.5038, A = 47484 * pv, B = 10 ** (-3083.87 / TR + 6.69449);
  return ((A / Ppsia + B) * 0.45359237) / (1e6 * 0.0283168 * (288.15 / 288.706));
}

// ---- hydrate screening curve ----------------------------------------------------------------------------
/** Hydrate temperature depression (°C) by salt and thermodynamic inhibitor in the aqueous phase (Nielsen–Bucklin on the water mole fraction). */
export function hydrateDepression(aq = {}) {
  const S = clamp(aq.S || 0, 0, 260) / 1000, w = clamp((aq.inhWt || 0) / 100, 0, 0.95), mwI = (aq.inh || INHIBITORS.none).MW;
  // 1 kg of aqueous phase: inhibitor w, the remainder is brine of salinity S
  const mBr = 1 - w, nW = (mBr * (1 - S)) / 0.018015, nS = (2 * mBr * S) / 0.05844, nI = w / (mwI * 1e-3), xw = nW / (nW + nS + nI);
  return -72 * Math.log(xw);
}
/** Inhibitor mass fraction (wt % of the aqueous phase) that gives a required depression (°C) on top of the salt already present. */
export function inhibitorFor(dT, inhId = 'MEG', S = 0) {
  const inh = INHIBITORS[inhId] || INHIBITORS.MEG, base = hydrateDepression({ S, inhWt: 0, inh });
  if (dT <= base) return 0;
  const g = (w) => hydrateDepression({ S, inhWt: w, inh }) - dT;
  return g(94) < 0 ? 94 : brent(g, 0, 94, 1e-8);
}
/** Hydrate equilibrium temperature (°C) of a natural gas of specific gravity sg at P (bara) with fresh water (Motiee correlation). */
export function hydrateT0(Pbar, sg = 0.65) {
  const lp = Math.log10(clamp(Pbar, 1, 700) * 14.5038), g = clamp(sg, 0.554, 1);
  return (-238.24469 + 78.99667 * lp - 5.352544 * lp * lp + 349.473877 * g - 150.854675 * g * g - 27.604065 * lp * g - 32) / 1.8;
}
/** Hydrate equilibrium temperature (°C) including salt and inhibitor. */
export const hydrateT = (Pbar, sg, aq) => hydrateT0(Pbar, sg) - (aq ? hydrateDepression(aq) : 0);
/** Hydrate equilibrium pressure (bara) at T (°C); null above the highest temperature of the curve. */
export function hydrateP(Tc, sg, aq) {
  const g = (P) => hydrateT(P, sg, aq) - Tc;
  if (g(1) > 0) return 1;
  if (g(700) < 0) return null;
  return brent(g, 1, 700, 1e-9);
}

// ---- tabulated properties for the flow solvers ------------------------------------------------------
const FIELDS = ['wG', 'rhoG', 'rhoO', 'muG', 'muO', 'cpG', 'cpO', 'kG', 'kO', 'hG', 'hO', 'jtG', 'jtO', 'sigma', 'zG', 'mwG', 'mwO'];
/**
 * Property table on a pressure–temperature grid. Returns { P[], T[], lnP[], <field>[iP][iT], spec, eosId }.
 * Missing-phase cells carry the values of the present phase so that interpolation across the phase boundary stays smooth.
 */
export function buildTable(spec = DEFAULT_FLUID, { nP = 22, nT = 17, Pmin = 1, Pmax = 600, Tmin = -30, Tmax = 170, opts = {}, onProgress } = {}) {
  const f = makeFluid(spec, opts), P = logspace(Pmin, Pmax, nP), T = linspace(Tmin, Tmax, nT), t = { P, T, lnP: P.map(Math.log), eosId: f.eosId };
  for (const k of FIELDS) t[k] = P.map(() => new Array(nT).fill(0));
  P.forEach((p, i) => {
    T.forEach((tc, j) => {
      const s = props(f, p, tc);
      t.wG[i][j] = s.wG; t.rhoG[i][j] = s.gas.rho; t.rhoO[i][j] = s.oil.rho; t.muG[i][j] = s.gas.mu; t.muO[i][j] = s.oil.mu; t.cpG[i][j] = s.gas.cp; t.cpO[i][j] = s.oil.cp; t.kG[i][j] = s.gas.k; t.kO[i][j] = s.oil.k;
      t.hG[i][j] = s.gas.h; t.hO[i][j] = s.oil.h; t.jtG[i][j] = s.gas.jt; t.jtO[i][j] = s.oil.jt; t.sigma[i][j] = s.sigma; t.zG[i][j] = s.gas.Z; t.mwG[i][j] = s.gas.MW; t.mwO[i][j] = s.oil.MW;
    });
    onProgress?.((i + 1) / nP);
  });
  // Where one phase is absent, give it physically sensible stand-in properties so that averaged mixture values never blow up.
  P.forEach((p, i) => T.forEach((tc, j) => {
    if (t.wG[i][j] >= 1) { t.rhoO[i][j] = Math.max(t.rhoO[i][j], 500); t.muO[i][j] = Math.max(t.muO[i][j], 2e-4); }
    if (t.wG[i][j] <= 0) { const z = 0.9, mw = 20; t.rhoG[i][j] = Math.min(t.rhoG[i][j], (p * 1e5 * mw * 1e-3) / (z * R * (tc + KEL))); t.muG[i][j] = Math.min(t.muG[i][j], 1.5e-5); t.zG[i][j] = z; t.mwG[i][j] = mw; t.cpG[i][j] = 2300; t.kG[i][j] = 0.035; t.jtG[i][j] = 3e-6; }
  }));
  const st = streams(spec, f), aq = aqueous(spec);
  t.spec = JSON.parse(JSON.stringify(spec)); t.gasSG = st.std.gasSG ?? f.MW / MW_AIR; t.rates = { nHC: st.nHC, mHC: st.mHC, mW: st.mW, qOilStd: st.qOilStd, qGasStd: st.qGasStd, qWaterStd: st.qWaterStd, gor: Number.isFinite(st.gor) ? st.gor : null, wc: st.wc, api: st.std.api, rhoOilStd: st.std.rhoOilStd, mwGas: st.std.mwGas };
  t.aq = { S: aq.S, inhWt: aq.inhWt, inhId: aq.inhId };
  return t;
}
/** Bilinear interpolation (in ln P and T) of every tabulated field: returns a plain object of properties at P (bara), T (°C). */
export function lookup(t, Pbar, Tc) {
  const lp = Math.log(clamp(Pbar, t.P[0], t.P[t.P.length - 1])), T = clamp(Tc, t.T[0], t.T[t.T.length - 1]), nP = t.P.length, nT = t.T.length;
  let i = Math.min(nP - 2, Math.max(0, Math.floor(((lp - t.lnP[0]) / (t.lnP[nP - 1] - t.lnP[0])) * (nP - 1)))), j = Math.min(nT - 2, Math.max(0, Math.floor(((T - t.T[0]) / (t.T[nT - 1] - t.T[0])) * (nT - 1))));
  const u = (lp - t.lnP[i]) / (t.lnP[i + 1] - t.lnP[i]), v = (T - t.T[j]) / (t.T[j + 1] - t.T[j]), o = {};
  for (const k of FIELDS) { const a = t[k]; o[k] = (1 - u) * ((1 - v) * a[i][j] + v * a[i][j + 1]) + u * ((1 - v) * a[i + 1][j] + v * a[i + 1][j + 1]); }
  o.wG = clamp(o.wG, 0, 1);
  return o;
}

const tableCache = new Map();
const specKey = (s) => JSON.stringify([s.comp, s.c7MW, s.c7SG, s.rateBasis, s.qOil, s.qGas, s.mdot, s.wc, s.qWater, s.salinity, s.inhibitor, s.inhWt, s.eos, s.nPseudo]);
/**
 * Fluid model used by every downstream suite. `ctx` is the suite run context ({ fluid, outputs }); the table published by
 * the PVT suite is used when it belongs to the same fluid, otherwise a table is built (and cached) from the case fluid.
 * Returns { spec, table, rates, aq, gasSG,
 *   at(P, T) -> { wG, rhoG, rhoO, rhoW, rhoL, muG, muO, muW, muL, cpG, cpO, cpW, cpL, kG, kO, kW, kL, hG, hO, hW, jtG, jtO, sigma (gas–liquid), zG, mwG, mwO,
 *                 mG, mO, mW (kg/s), qG, qO, qW, qL (actual m³/s), wcut (in-situ water fraction of the liquid), phaseInv },
 *   hydrateT(P) (°C, with the salt and inhibitor of the case), hydrateT0(P) (fresh water), hydrateP(T) }.
 */
export function fluidModel(ctx = {}, override = {}) {
  const spec = { ...DEFAULT_FLUID, ...(ctx.fluid || {}), ...override, comp: { ...DEFAULT_FLUID.comp, ...((ctx.fluid || {}).comp || {}), ...(override.comp || {}) } }, key = specKey(spec);
  let table = null;
  const pub = ctx.outputs?.pvt?.table;
  if (pub && pub.spec && specKey({ ...DEFAULT_FLUID, ...pub.spec }) === key && Array.isArray(pub.P) && pub.rhoG) { table = pub; if (!table.lnP) table.lnP = table.P.map(Math.log); }
  if (!table) { table = tableCache.get(key); if (!table) { table = buildTable(spec); if (tableCache.size > 6) tableCache.clear(); tableCache.set(key, table); } }
  const aq = aqueous(spec), rates = table.rates, gasSG = table.gasSG || 0.7, curve = ctx.outputs?.pvt?.hydrateCurve, useCurve = pub === table && curve && Array.isArray(curve.P) && curve.P.length > 3;
  const curveT = (P, arr) => { const lp = Math.log(clamp(P, curve.P[0], curve.P[curve.P.length - 1])); let i = 0; while (i < curve.P.length - 2 && Math.log(curve.P[i + 1]) < lp) i++; const a = Math.log(curve.P[i]), b = Math.log(curve.P[i + 1]); return arr[i] + ((arr[i + 1] - arr[i]) * (lp - a)) / (b - a || 1); };
  const hT0 = (P) => (useCurve && curve.T0 ? curveT(P, curve.T0) : hydrateT0(P, gasSG)), dep = hydrateDepression(aq);
  const model = {
    spec, table, rates, aq, gasSG,
    at(P, T, mScale = 1) {
      const o = lookup(table, P, T), w = waterProps(P, T, aq), mHC = rates.mHC * mScale, mW = rates.mW * mScale;
      o.rhoW = w.rho; o.muW = w.mu; o.cpW = w.cp; o.kW = w.k; o.hW = w.h;
      o.mG = mHC * o.wG; o.mO = mHC * (1 - o.wG); o.mW = mW;
      o.qG = o.mG / o.rhoG; o.qO = o.mO / o.rhoO; o.qW = mW / w.rho; o.qL = o.qO + o.qW;
      const wcut = o.qL > 0 ? o.qW / o.qL : 0, mL = o.mO + mW;
      o.wcut = wcut; o.rhoL = o.qL > 0 ? mL / o.qL : o.rhoO;
      // liquid viscosity: continuous phase with the Brinkman emulsion factor, inversion at 60 % water
      o.phaseInv = wcut > 0.6; o.muL = wcut <= 0 ? o.muO : wcut >= 1 ? o.muW : o.phaseInv ? o.muW * (1 - Math.min(1 - wcut, 0.7)) ** -2.5 : o.muO * (1 - Math.min(wcut, 0.7)) ** -2.5;
      o.cpL = mL > 0 ? (o.mO * o.cpO + mW * o.cpW) / mL : o.cpO; o.kL = (1 - wcut) * o.kO + wcut * o.kW;
      o.sigmaGW = clamp(0.0756 - 1.5e-4 * T - 2e-5 * P, 0.03, 0.076); o.sigmaOW = clamp(0.032 - 1e-4 * (T - 20), 0.012, 0.035);
      o.sigma = wcut > 0.6 ? o.sigmaGW : Math.max(o.sigma, 1e-4);
      return o;
    },
    hydrateT0: hT0,
    hydrateT: (P) => (useCurve && curve.T ? curveT(P, curve.T) : hT0(P) - dep),
    hydrateP: (T) => { const g = (P) => model.hydrateT(P) - T; if (g(1) > 0) return 1; if (g(700) < 0) return null; return brent(g, 1, 700, 1e-8); },
    depression: dep,
  };
  return model;
}
