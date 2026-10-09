// Suite 4 — Hydrate & multiphase solids flow assurance.
// Hydrate thermodynamics (van der Waals–Platteeuw with Parrish–Prausnitz or Kihara constants, Gibbs minimisation),
// heterogeneous nucleation at the measured rate (stochastic onset, ramp integration; classical theory for comparison),
// intrinsic / transfer-limited growth with an exact shrinking-core advance, coupled to the slug unit cell, population balances
// (sectional, quadrature moments, Monte Carlo, and sectional on a resolved axisymmetric flow field), slurry rheology, wall deposition and plugging marched in time on the
// case line with injection, heating, depressurisation and equipment boundaries; Lagrangian particle tracking
// (Maxey–Riley) and Eulerian–Eulerian solids transport; wax solid–liquid equilibrium and deposition; mineral scale
// with the PHREEQC Pitzer and ion-association models, induction time and threshold inhibition, crystallisation moments;
// asphaltene onset; sand transport by flow regime; preservation of the shut-in line with the unprotected event as comparison;
// one combined deposit profile for the backward coupling to the network and flow suites.
// Literature constants are listed with their sources in PROVENANCE; measured reference data live in ../data/ref/solids.js.
// SI units inside; bara and °C at the interfaces.
import { brent, clamp, interp1, linspace, logspace, rng, histogram, mean, quantile, isNum, rk45, tridiag, metrics } from '../core/num.js';
import { R, INHIBITORS, VM_STD, DEFAULT_FLUID, makeFluid, eosPhase, flashPT, fluidModel, waterContent, hydrateDepression, inhibitorFor, hydrateT0 } from '../core/thermo.js';
import { G, gradient, frictionFactor, hInside, slugUnit } from '../core/pipe.js';
import { waxModel, scnDistribution } from './s01_pvt.js';
import { caseLine, steadyCase, flowPicture } from '../core/caseflow.js';
import { BASE } from '../data/basecase.js';
import { density as waterDensity } from '../core/props.js';
import { PHREEQC, REF_SETS, PHREEQC_WATERS } from '../data/ref/solids.js';

const KEL = 273.15, KB = 1.380649e-23, MW_W = 0.018015, PI = Math.PI;
/** Hydrate solid properties used throughout: density of structure-II hydrate, latent heat, thermal conductivity and heat capacity of structure I (sources in PROVENANCE). */
export const HYDRATE = Object.freeze({ rho: 914, latent: 4.6e5, k: 0.49, cp: 2080 }); // kg/m³, J/kg, W/m/K, J/kg/K
const need = (cond, msg) => { if (!cond) throw new Error(msg); };
const num = (x, d) => (isNum(+x) ? +x : d);

// =====================================================================================================
// 1. Hydrate thermodynamics: van der Waals–Platteeuw statistical model for methane structure I,
//    Langmuir constants from the Parrish–Prausnitz fit or from the Kihara cell potential, Gibbs-energy minimisation
// =====================================================================================================
let c1Fluid = null;
/** Fugacity (bar) of pure methane from the kernel Peng–Robinson model. */
export function methaneFugacity(Pbar, TK) { c1Fluid ||= makeFluid({ comp: { C1: 100 } }); return Pbar * Math.exp(eosPhase(c1Fluid, c1Fluid.z, Pbar, TK, 'vapour').lnphi[0]); }
/** Constants of the statistical hydrate model (structure I, methane). Sources and checks are listed in PROVENANCE. */
export const HYD_REF = Object.freeze({
  // Parrish & Prausnitz (1972) Langmuir fit C = (A/T) exp(B/T), 1/atm
  ppSmall: [3.7237e-3, 2708.8], ppLarge: [1.8372e-2, 2737.9],
  // empty lattice − liquid water at 273.15 K: Δμ⁰ (J/mol), Δh⁰ (J/mol), Δcp = a + b (T − 273.15) (J/mol/K), Δv (m³/mol)
  dmu0: 1264, dh0: -4858, dcpA: -39.16, dcpB: 0, dv: 4.6e-6,
  // Kihara spherical-core parameters of methane in the hydrate lattice and the structure-I cages (radius m, coordination number, cages per water molecule)
  kihara: { a: 0.3834e-10, sigma: 3.1650e-10, epsK: 154.54 },
  // structure-I cages as single water shells: mean radius (m), coordination number, cages per water molecule
  cages: { small: { R: 3.95e-10, z: 20, nu: 2 / 46 }, large: { R: 4.33e-10, z: 24, nu: 6 / 46 } },
});
/**
 * Langmuir constant (1/atm) of a guest in a hydrate cage from the Kihara spherical-core cell potential
 * (Lennard-Jones–Devonshire smearing, McKoy–Sinanoğlu), summed over the water shells of the cage:
 * C = 4π/(kT) ∫₀^{R₁−a} exp(−Σ w_k(r)/kT) r² dr with
 * w_k(r) = 2 z_k ε [σ¹²/(R_k¹¹ r)(δ¹⁰ + a/R_k δ¹¹) − σ⁶/(R_k⁵ r)(δ⁴ + a/R_k δ⁵)], δᴺ = [(1 − r/R_k − a/R_k)⁻ᴺ − (1 + r/R_k − a/R_k)⁻ᴺ]/N.
 * cage: { shells: [[R (m), z], …] } (or { R, z } for one shell); guest: { a (m), sigma (m), epsK (K) }.
 */
export function kiharaLangmuir(TK, cage, guest = HYD_REF.kihara, n = 240) {
  const shells = cage.shells || [[cage.R, cage.z]], { a, sigma, epsK } = guest, rMax = shells[0][0] - a;
  const wk = (r) => { let w = 0; for (const [Rc, z] of shells) { const x = r / Rc, aR = a / Rc, dl = (N) => ((1 - x - aR) ** -N - (1 + x - aR) ** -N) / N; w += 2 * z * epsK * ((sigma ** 12 / (Rc ** 11 * r)) * (dl(10) + aR * dl(11)) - (sigma ** 6 / (Rc ** 5 * r)) * (dl(4) + aR * dl(5))); } return w; }; // in kelvin
  let sum = 0; const h = rMax / n;
  for (let k = 0; k < n; k++) { const r = (k + 0.5) * h, e = -wk(r) / TK; if (e > -700) sum += Math.exp(Math.min(e, 700)) * r * r * h; } // midpoint rule: the integrand vanishes at the cage wall
  return ((4 * PI * sum) / (KB * TK)) * 101325;
}
/**
 * Langmuir cage occupancies of methane in structure I. model: 'parrish' (Parrish–Prausnitz constants, C = A/T·exp(B/T) in 1/atm)
 * or 'kihara' (cell-potential integral). Returns { Cs, Cl (1/atm), thetaS, thetaL, hydrationNumber = 46 / (2θs + 6θl) }.
 */
export function langmuirOccupancy(TK, fBar, model = 'parrish') {
  const f = fBar / 1.01325, kh = model === 'kihara', Cs = kh ? kiharaLangmuir(TK, HYD_REF.cages.small) : (HYD_REF.ppSmall[0] / TK) * Math.exp(HYD_REF.ppSmall[1] / TK), Cl = kh ? kiharaLangmuir(TK, HYD_REF.cages.large) : (HYD_REF.ppLarge[0] / TK) * Math.exp(HYD_REF.ppLarge[1] / TK);
  const thetaS = (Cs * f) / (1 + Cs * f), thetaL = (Cl * f) / (1 + Cl * f);
  return { Cs, Cl, thetaS, thetaL, hydrationNumber: 46 / Math.max(2 * thetaS + 6 * thetaL, 1e-9) };
}
/** (μ_w^β − μ_w^L,pure)/RT of water between the empty hydrate lattice and pure liquid water at T (K), P (bara). */
export function emptyLatticeDmu(T, Pbar, H = HYD_REF) {
  const T0 = KEL, n = 60; let I = 0; // ∫ Δh/(R T²) dT from T0 to T
  for (let k = 0; k < n; k++) { const t = T0 + ((T - T0) * (k + 0.5)) / n, dh = H.dh0 + H.dcpA * (t - T0) + 0.5 * H.dcpB * (t - T0) ** 2; I += ((dh / (R * t * t)) * (T - T0)) / n; }
  return H.dmu0 / (R * T0) - I + (H.dv * Pbar * 1e5) / (R * T);
}
/**
 * Three-phase (liquid water – hydrate – vapour) equilibrium pressure of methane hydrate from the equality of the
 * chemical potential of water in the hydrate lattice and in the aqueous phase (van der Waals–Platteeuw, sI, T ≥ 0 °C).
 * o: { langmuir: 'parrish' | 'kihara', aw (activity of water in the aqueous phase, 1 = fresh water) }.
 * Returns { P (bara) | null, thetaS, thetaL, hydrationNumber }.
 */
export function vdwpMethane(Tc, { langmuir = 'parrish', aw = 1 } = {}) {
  const T = Tc + KEL;
  const g = (P) => {
    const f = methaneFugacity(P, T), o = langmuirOccupancy(T, f, langmuir), xg = (f * 1e5) / (4.0e9 * Math.exp(-1700 * (1 / T - 1 / 298.15))); // dissolved methane lowers the water activity slightly
    return -((2 / 46) * Math.log(1 - o.thetaS) + (6 / 46) * Math.log(1 - o.thetaL)) - (emptyLatticeDmu(T, P) - Math.log(1 - xg) - Math.log(aw));
  };
  if (!(T >= 272.5) || g(1) > 0 || g(900) < 0) return { P: null, thetaS: null, thetaL: null, hydrationNumber: null };
  const P = brent(g, 1, 900, 1e-8), o = langmuirOccupancy(T, methaneFugacity(P, T), langmuir);
  return { P, thetaS: o.thetaS, thetaL: o.thetaL, hydrationNumber: o.hydrationNumber };
}
/**
 * Hydrate phase amounts by direct minimisation of the Gibbs energy of a closed water–methane–solute system at fixed T, P
 * (gas in excess at its fugacity; the aqueous phase is an ideal solution of water and dissolved particles, so water that
 * goes into hydrate concentrates the salt or inhibitor and raises its own resistance to further conversion).
 *   G(ξ)/RT = ξ·Δg + (n_w − ξ) ln x_w + n_s ln x_s,  Δg = (μ_w^H − μ_w^L,pure)/RT,  0 ≤ ξ ≤ min(n_w, n·n_gas)
 * o: { Tc, P (bara), nW, nGas (mol), nSolute (mol of dissolved particles: 2 per NaCl, 1 per inhibitor molecule), langmuir }.
 * Returns { xi (mol water in hydrate at the minimum), conversion, xw (final water mole fraction), dg (Δg), G0, Gmin, hydrationNumber, limiting: 'none' | 'equilibrium' | 'water' | 'gas', xiAnalytic }.
 */
export function hydrateGibbsMin({ Tc, P, nW = 1, nGas = 1, nSolute = 0, langmuir = 'parrish' }) {
  const T = Tc + KEL, f = methaneFugacity(P, T), o = langmuirOccupancy(T, f, langmuir), dg = (2 / 46) * Math.log(1 - o.thetaS) + (6 / 46) * Math.log(1 - o.thetaL) + emptyLatticeDmu(T, P);
  const hi = Math.min(nW * (1 - 1e-12), nGas * o.hydrationNumber), Gf = (xi) => { const w = nW - xi, tot = w + nSolute; return xi * dg + (w > 0 ? w * Math.log(w / tot) : 0) + (nSolute > 0 ? nSolute * Math.log(nSolute / tot) : 0); };
  let a = 0, b = hi; const gr = (Math.sqrt(5) - 1) / 2; let c = b - gr * (b - a), d = a + gr * (b - a), fc = Gf(c), fd = Gf(d);
  for (let k = 0; k < 200 && b - a > 1e-13 * Math.max(nW, 1e-12); k++) { if (fc < fd) { b = d; d = c; fd = fc; c = b - gr * (b - a); fc = Gf(c); } else { a = c; c = d; fc = fd; d = a + gr * (b - a); fd = Gf(d); } }
  let xi = 0.5 * (a + b); if (Gf(0) <= Gf(xi)) xi = 0; if (Gf(hi) < Gf(xi)) xi = hi;
  const xs = Math.exp(dg), xiA = dg >= 0 ? 0 : nSolute > 0 ? clamp(nW - (nSolute * xs) / (1 - xs), 0, hi) : hi; // stationarity: ln x_w = Δg
  return { xi, conversion: xi / nW, xw: (nW - xi) / (nW - xi + nSolute), dg, G0: Gf(0), Gmin: Gf(xi), hydrationNumber: o.hydrationNumber, limiting: xi <= 1e-12 * nW ? 'none' : xi >= hi * (1 - 1e-9) ? (hi < nW * (1 - 1e-9) ? 'gas' : 'water') : 'equilibrium', xiAnalytic: xiA };
}

// =====================================================================================================
// 2. Nucleation (classical nucleation theory) and induction time
// =====================================================================================================
/**
 * Classical nucleation rate. o: { TK, dT (K subcooling), TeqK, sigma (J/m², hydrate–water), theta (deg contact angle),
 * A (pre-exponential, 1/m³/s), het (true = heterogeneous) }. The volumetric driving force is ρ·L·ΔT/Teq.
 * Returns { J (1/m³ water/s), dG (J, barrier), rc (m, critical radius), f (contact-angle factor), exponent }.
 */
export function nucleationRate({ TK, dT, TeqK, sigma = 0.032, theta = 40, A = 3e7, het = true }) {
  const c = Math.cos((clamp(theta, 0, 180) * PI) / 180), f = het ? ((2 + c) * (1 - c) ** 2) / 4 : 1;
  if (!(dT > 0)) return { J: 0, dG: Infinity, rc: Infinity, f, exponent: Infinity };
  const dgv = (HYDRATE.rho * HYDRATE.latent * dT) / TeqK, dG = ((16 * PI * sigma ** 3) / (3 * dgv * dgv)) * f, ex = dG / (KB * TK);
  return { J: A * Math.exp(-Math.min(ex, 700)), dG, rc: (2 * sigma) / dgv, f, exponent: ex };
}
/** Mean induction time (s) of a water sample of volume V (m³) at a subcooling dT: 1 / (J·V), the mean of the exponential (Poisson) waiting time. */
export function inductionTime(dT, { TK = 277.15, TeqK, sigma, theta, A, het, V = 1e-3 } = {}) {
  const J = nucleationRate({ TK, dT, TeqK: TeqK ?? TK + dT, sigma, theta, A, het }).J;
  return J * V > 1e-30 ? 1 / (J * V) : 1e30;
}

/**
 * Measured heterogeneous nucleation of gas hydrate at a gas–water interface (stirred high-pressure cells, thousands of
 * formation events): two populations of nucleation sites,
 *   J/α = Σ A_i · exp(Δs_e ΔT / (k T)) · exp(−B′_i / (T ΔT²)),   ΔT = T_eq − T,
 * per unit of gas–water interfacial area α. `low` governs the first few kelvin of subcooling (isothermal induction times,
 * ΔT 2–4 K), `high` the fast formation seen in rapid ramps (ΔT above about 5 K). A per cell (s⁻¹) is divided by the
 * interface area of the 10.9 cm³ cell it was measured in (alpha0). dse = entropy of dissociation per hydrate building unit (k).
 */
export const NUCLEATION = Object.freeze({ low: Object.freeze({ A: 0.58e-3, B: 2.66e3 }), high: Object.freeze({ A: 0.10, B: 3.1e4 }), dse: 22.2, alpha0: 11e-4 });
/**
 * Heterogeneous nucleation rate per unit interfacial area (1/m²/s) at a subcooling dT (K) below the equilibrium temperature TeqK.
 * o: { dT, TeqK, mult (multiplier on both pre-factors: site density, agitation), A1, B1, A2, B2 (override the published sets; A per m²) }
 */
export function hydrateNucleationRate({ dT, TeqK = 283.15, mult = 1, A1 = NUCLEATION.low.A / NUCLEATION.alpha0, B1 = NUCLEATION.low.B, A2 = NUCLEATION.high.A / NUCLEATION.alpha0, B2 = NUCLEATION.high.B, dse = NUCLEATION.dse }) {
  if (!(dT > 0)) return 0;
  const T = Math.max(TeqK - dT, 150), e = Math.exp(Math.min((dse * dT) / T, 50)), q = T * dT * dT;
  return mult * e * (A1 * Math.exp(-Math.min(B1 / q, 700)) + A2 * Math.exp(-Math.min(B2 / q, 700)));
}
/** Mean induction time (s) at a constant subcooling for an interface of `area` m²: 1/(J·α), the mean of the exponential waiting time. */
export function hydrateInductionTime(dT, { TeqK = 283.15, area = NUCLEATION.alpha0, ...o } = {}) {
  const J = hydrateNucleationRate({ dT, TeqK, ...o }) * area;
  return J > 1e-30 ? 1 / J : 1e30;
}
/**
 * Onset statistics of a constant-cooling ramp started at the equilibrium temperature: the survival probability is
 * exp(−(α/β) ∫₀^ΔT J dΔT′) for a cooling rate β. o: { rate (K/s), area (m²), TeqK, …rate options, dTmax, n }.
 * Returns { median, mean, sd, p10, p90 (K of subcooling at onset), quantile(p), cdf(dT), hazard: { dT[], H[] } }.
 */
export function onsetRamp({ rate, area = NUCLEATION.alpha0, TeqK = 283.15, dTmax = 30, n = 3000, ...o }) {
  need(rate > 0 && area > 0, 'onsetRamp needs a positive cooling rate and interface area.');
  const h = dTmax / n, dTs = new Float64Array(n + 1), H = new Float64Array(n + 1); let mean = 0, m2 = 0;
  for (let k = 1; k <= n; k++) { const x = k * h, J = hydrateNucleationRate({ dT: x - 0.5 * h, TeqK, ...o }); dTs[k] = x; H[k] = H[k - 1] + (J * h * area) / rate; // midpoint rule
    const s0 = Math.exp(-H[k - 1]), s1 = Math.exp(-H[k]), xm = x - 0.5 * h; mean += (s0 - s1) * xm; m2 += (s0 - s1) * xm * xm; }
  const left = Math.exp(-H[n]); mean += left * dTmax; m2 += left * dTmax * dTmax; // runs that have not formed by dTmax are counted there
  const quantile = (p) => { const t = -Math.log(1 - clamp(p, 1e-12, 1 - 1e-12)); if (H[n] < t) return dTmax; let lo = 0, hi = n; while (hi - lo > 1) { const m = (lo + hi) >> 1; if (H[m] < t) lo = m; else hi = m; } return dTs[lo] + ((t - H[lo]) / Math.max(H[hi] - H[lo], 1e-300)) * h; };
  const cdf = (x) => { const q = clamp(x / h, 0, n), k = Math.min(Math.floor(q), n - 1); return 1 - Math.exp(-(H[k] + (H[k + 1] - H[k]) * (q - k))); };
  return { median: quantile(0.5), mean, sd: Math.sqrt(Math.max(m2 - mean * mean, 0)), p10: quantile(0.1), p90: quantile(0.9), quantile, cdf, hazard: { dT: Array.from(dTs), H: Array.from(H) } };
}
// =====================================================================================================
// 3. Growth and dissociation kinetics
// =====================================================================================================
const fug = (Pbar, z) => Pbar * 1e5 * Math.exp(clamp(z - 1, -1.2, 0.3)); // Pa, first-order virial fugacity coefficient ln φ = Z − 1
const fugAt = (Pq, P, z) => Pq * 1e5 * Math.exp(clamp(((z - 1) * Pq) / Math.max(P, 1e-6), -1.2, 0.3));
/**
 * Hydrate growth flux on a particle surface with resistances in series (Kim–Bishnoi / Englezos intrinsic step on the
 * fugacity difference, liquid-film mass transfer, diffusion through a hydrate shell, heat removal from the particle).
 * o: { TK, P, Peq (bara), zG, kRef (mol/m²/Pa/s at 277.15 K), EaR (K), H (Pa·m³/mol, Henry constant of the gas in the
 *      continuous liquid), kFilm (m/s), kShell (m/s, Infinity = no shell), hPart (W/m²/K), dT (K), dHmol (J/mol gas) }
 * Returns { j (mol gas/m²/s, ≥ 0), jKin, jFilm, jShell, jHeat, kp (m/s, kinetic+film+shell conductance), dc (mol/m³), df (Pa), limiting }.
 */
export function hydrateGrowthRate({ TK, P, Peq, zG = 0.85, kRef = 1e-10, EaR = 13600, H = 2500, kFilm = Infinity, kShell = Infinity, hPart = Infinity, dT = 0, dHmol = 6e4 }) {
  const df = fug(P, zG) - fugAt(Peq, P, zG), z = { j: 0, jKin: 0, jFilm: 0, jShell: 0, jHeat: 0, kp: 0, dc: 0, df, limiting: 'none' };
  if (!(df > 0) || !(kRef > 0)) return z;
  const kStar = kRef * Math.exp(-EaR * (1 / TK - 1 / 277.15)), dc = df / H, kKin = kStar * H;
  const jKin = kKin * dc, jFilm = kFilm * dc, jShell = kShell * dc, jHeat = Number.isFinite(hPart) ? (hPart * Math.max(dT, 0)) / dHmol : Infinity;
  const kp = 1 / (1 / kKin + 1 / kFilm + 1 / kShell), j = 1 / (1 / (kp * dc) + 1 / jHeat);
  const m = Math.min(jKin, jFilm, jShell, jHeat);
  return { j, jKin, jFilm, jShell, jHeat, kp, dc, df, limiting: m === jKin ? 'intrinsic kinetics' : m === jFilm ? 'mass transfer' : m === jShell ? 'shell diffusion' : 'heat transfer' };
}
/**
 * Kim–Bishnoi dissociation flux (mol gas/m²/s) with an Arrhenius constant: K0·exp(−E/RT)·(f_eq − f).
 * Defaults are the methane values of Clarke & Bishnoi (K0 = 3.6e4 mol/m²/Pa/s, E = 81 kJ/mol).
 */
export function hydrateDissociationRate({ TK, P, Peq, zG = 0.85, K0 = 3.6e4, E = 81e3 }) {
  return K0 * Math.exp(-E / (R * TK)) * Math.max(fugAt(Peq, P, zG) - fug(P, zG), 0);
}
/** Shrinking-core shell conductance (m/s on the outer surface) for a converted fraction X of a droplet of radius Rd. */
export const shellConductance = (X, Rd, Dshell) => { const rc = Rd * Math.cbrt(clamp(1 - X, 0, 1)); return rc <= 0 ? 0 : Rd - rc < 1e-12 * Rd ? Infinity : (Dshell * rc) / (Rd * (Rd - rc)); };
/**
 * Exact advance of a shrinking-core conversion over a time step. With a constant resistance Ra (s/m: surface kinetics, liquid
 * film, heat removal) in series with the growing shell, (Ra + Rd(1 − y)/(D y)) dX = c dt with y = (1 − X)^⅓ integrates to
 *   G(X) = Ra·X + (Rd/2D)(1 − 3y² + 2y³),   G(X₁) = G(X₀) + a,   a = c·Δt  (m·s/m, i.e. driving concentration × volume factor × time).
 * An explicit step with the rate at X₀ overshoots by orders of magnitude while the shell is thin; this form does not.
 * Returns X₁ (≤ 1).
 */
export function shrinkingCoreAdvance(X0, a, Ra, Rd, D) {
  if (!(a > 0) || !Number.isFinite(Ra)) return X0; // no driving force, or a blocked step (zero rate constant): nothing converts
  const kd = D > 0 && Number.isFinite(D) ? Rd / (2 * D) : 0, Gf = (X) => { const y = Math.cbrt(Math.max(1 - X, 0)); return Ra * X + kd * (1 - 3 * y * y + 2 * y * y * y); }, target = Gf(clamp(X0, 0, 1)) + a;
  if (Gf(1) <= target) return 1;
  let lo = clamp(X0, 0, 1), hi = 1; for (let it = 0; it < 60; it++) { const m = 0.5 * (lo + hi); if (Gf(m) < target) lo = m; else hi = m; if (hi - lo < 1e-13) break; }
  return 0.5 * (lo + hi);
}
/** Exact advance of a flat film of thickness d0 (m) growing through its own diffusion resistance: (Ra + δ/D) dδ = a/Δt·dt → δ₁. */
export const filmAdvance = (d0, a, Ra, D) => (!(a > 0) || !Number.isFinite(Ra) ? d0 : D > 0 && Number.isFinite(D) ? -Ra * D + Math.sqrt((Ra * D + d0) ** 2 + 2 * D * a) : d0 + a / Math.max(Ra, 1e-300));

// =====================================================================================================
// 4. Population balance: sectional (fixed pivot) and quadrature method of moments
// =====================================================================================================
/**
 * Geometric size grid with the fixed-pivot allocation tables for binary aggregation and binary breakage.
 * dCol: collision diameter floor (m) — particles smaller than this collide as if they had this size.
 * Returns { n, L[], v[], pk, pa (pair → lower pivot and number fraction assigned to it), bk, ba (breakage), gS, gD, gB (kernel geometry) }.
 */
export function pbeGrid(n = 14, Lmin = 2e-6, Lmax = 2e-2, dCol = 0) {
  n = Math.max(3, Math.round(n));
  const L = logspace(Lmin, Lmax, n), v = L.map((x) => (PI / 6) * x ** 3), find = (vs) => { let k = 0; while (k < n - 2 && v[k + 1] <= vs) k++; return k; };
  const pk = new Int32Array(n * n), pa = new Float64Array(n * n), pb = new Float64Array(n * n), gS = new Float64Array(n * n), gD = new Float64Array(n * n), gB = new Float64Array(n * n), bk = new Int32Array(n), ba = new Float64Array(n), bb = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      const vs = v[i] + v[j], q = i * n + j, a = Math.max(L[i], dCol), b = Math.max(L[j], dCol);
      if (vs >= v[n - 1]) { pk[q] = n - 2; pa[q] = 0; pb[q] = vs / v[n - 1]; } // overflow: volume kept in the top class
      else { const k = find(vs); pk[q] = k; pa[q] = (v[k + 1] - vs) / (v[k + 1] - v[k]); pb[q] = 1 - pa[q]; }
      gS[q] = (a + b) ** 3; gD[q] = (a + b) ** 2 * Math.abs(a * a - b * b); gB[q] = (a + b) ** 2 / (a * b);
    }
    const vd = v[i] / 2; // two equal daughters shared between the neighbouring pivots: conserves number (2) and volume
    if (vd < v[0]) { bk[i] = -1; } else { const k = find(vd); bk[i] = k; ba[i] = (2 * (v[k + 1] - vd)) / (v[k + 1] - v[k]); bb[i] = 2 - ba[i]; }
  }
  let ub = 1; for (let i = 0; i < n; i++) if (bk[i] >= 0) ub = Math.max(ub, i - bk[i]); // classes a daughter can fall below its parent
  return { n, L, v, pk, pa, pb, bk, ba, bb, gS, gD, gB, ub };
}
/**
 * Collision-frequency kernel (m³/s) between spheres of diameter Li and Lj and its parts.
 * e: { shear (1/s, laminar velocity gradient), eps (W/kg, turbulent dissipation), nu (m²/s), mu (Pa·s), TK, dRho (kg/m³), alpha (collision efficiency) }
 * Returns { total, shear (Smoluchowski orthokinetic), turbulent (Saffman–Turner), settling (differential Stokes), brownian }.
 */
export function aggregationKernel(Li, Lj, { shear = 0, eps = 0, nu = 1e-6, mu = 1e-3, TK = 277, dRho = 100, alpha = 1 } = {}) {
  const s = Li + Lj, sh = (shear / 6) * s ** 3, tu = 0.1618 * Math.sqrt(eps / nu) * s ** 3, se = (PI / 4) * s * s * ((Math.abs(dRho) * G) / (18 * mu)) * Math.abs(Li * Li - Lj * Lj), br = ((2 * KB * TK) / (3 * mu)) * (s * s) / (Li * Lj);
  return { total: alpha * (sh + tu + se + br), shear: sh, turbulent: tu, settling: se, brownian: br };
}
/**
 * Sectional population balance over a time t (s) in one well-mixed volume: aggregation (fixed pivot, explicit sub-steps),
 * binary breakage, growth/shrinkage in volume space (pivot shift with the exact volume rate and, for growth, exact number) and a source in one class.
 * N: number per class (any consistent basis). o: { beta: Float64Array n×n (m³/s on the same basis) | number (constant kernel) |
 *   [cS, cD, cB] (coefficients of the grid's shear, differential-settling and Brownian geometry tables),
 *   gv: [dv/dt per class, m³/s], S: [breakage frequency 1/s], src: { k, rate (1/s) }, frac (max fractional loss per sub-step), maxSub }.
 * Returns { N, sub, limited } — `limited` is true when the interval was too stiff for maxSub explicit sub-steps and
 * linearly implicit (modified Patankar–Euler) steps on the class volumes were used instead: unconditionally positive,
 * volume-conserving and convergent to the aggregation–breakage equilibrium.
 */
export function solvePBE(grid, N0, t, { beta = null, gv = null, S = null, src = null, frac = 0.25, maxSub = 400 } = {}) {
  const { n, v, pk, pa, pb, bk, ba, bb } = grid, ws = (grid.ws ||= { dN: new Float64Array(n), tmp: new Float64Array(n), dth: new Float64Array(n), lim: new Float64Array(n), V: new Float64Array(n), M: new Float64Array(n * n) });
  const N = new Float64Array(n), { dN, tmp, dth, lim, V, M } = ws; for (let i = 0; i < n; i++) N[i] = N0[i];
  const cst = typeof beta === 'number', kc = beta !== null && !cst && beta.length === 3 ? beta : null, { gS, gD, gB } = grid, k0 = kc ? kc[0] : 0, k1 = kc ? kc[1] : 0, k2 = kc ? kc[2] : 0;
  const B = (i, j) => { if (cst) return beta; const q = i * n + j; return kc ? k0 * gS[q] + k1 * gD[q] + k2 * gB[q] : beta[q]; };
  // death frequency of every class and the largest one among the classes that carry a noticeable share of the volume
  const rates = () => { let r = 0, vt = 0; for (let i = 0; i < n; i++) vt += N[i] * v[i]; for (let i = 0; i < n; i++) { let d = 0; if (beta !== null && N[i] > 0) for (let j = 0; j < n; j++) if (N[j] > 0) d += B(i, j) * N[j]; dth[i] = d; const tot = d + (S ? S[i] : 0); if (N[i] * v[i] > 1e-5 * vt && tot > r) r = tot; } return r; };
  let rate = rates(), sub = Math.max(1, Math.ceil((t * rate) / frac)), limited = false, dt = t / sub;
  if (sub > maxSub) { // too stiff for explicit sub-steps: linearly implicit, positivity-preserving and volume-conserving steps
    limited = true; sub = 0;
    const nImp = 3, h = t / nImp;
    const ub = grid.ub, flow = (src, dst, fl) => { if (src === dst) return; const c = (h * fl) / V[src]; if (!(c > 1e-18)) return; M[dst * n + src] -= c; M[src * n + src] += c; }; // negligible transfers are skipped (they would only cost denormal arithmetic)
    for (let st = 0; st < nImp; st++) {
      M.fill(0); for (let i = 0; i < n; i++) { V[i] = N[i] * v[i]; M[i * n + i] = 1; }
      if (beta !== null) for (let i = 0; i < n; i++) {
        if (!(N[i] > 0)) continue;
        for (let j = i; j < n; j++) {
          if (!(N[j] > 0)) continue;
          const q = i * n + j, r = (i === j ? 0.5 : 1) * B(i, j) * N[i] * N[j], k = pk[q], vs = v[i] + v[j], wk = (pa[q] * v[k]) / vs, wk1 = (pb[q] * v[k + 1]) / vs;
          flow(i, k, r * v[i] * wk); flow(i, k + 1, r * v[i] * wk1); flow(j, k, r * v[j] * wk); flow(j, k + 1, r * v[j] * wk1);
        }
      }
      if (S) for (let i = 0; i < n; i++) { const k = bk[i]; if (k < 0 || !(S[i] > 0) || !(N[i] > 0)) continue; const r = S[i] * N[i]; flow(i, k, r * ba[i] * v[k]); flow(i, k + 1, r * bb[i] * v[k + 1]); }
      // (I − h·M) V' = V by elimination. The matrix is lower triangular plus `ub` super-diagonals (aggregation only moves volume
      // up, a breakage daughter falls at most `ub` classes) and column-diagonally dominant, so no pivoting and O(n²·ub) work.
      for (let c = 0; c < n; c++) { const pv = M[c * n + c], ce = Math.min(c + ub, n - 1); for (let r2 = c + 1; r2 < n; r2++) { const m = M[r2 * n + c] / pv; if (m === 0) continue; for (let cc = c + 1; cc <= ce; cc++) M[r2 * n + cc] -= m * M[c * n + cc]; V[r2] -= m * V[c]; } }
      for (let r2 = n - 1; r2 >= 0; r2--) { let x = V[r2]; const ce = Math.min(r2 + ub, n - 1); for (let cc = r2 + 1; cc <= ce; cc++) x -= M[r2 * n + cc] * V[cc]; V[r2] = x / M[r2 * n + r2]; }
      for (let i = 0; i < n; i++) N[i] = Math.max(V[i], 0) / v[i];
    }
  }
  for (let s = 0; s < sub; s++) {
    if (s) rate = rates();
    if (rate > 0 && dt > 0) {
      dN.fill(0);
      // sparse classes with a very high death frequency are depleted exponentially, never below zero (pair-consistent, so volume is conserved)
      for (let i = 0; i < n; i++) { const x = dth[i] * dt; lim[i] = x > 1e-6 ? (1 - Math.exp(-x)) / x : 1; }
      if (beta !== null) for (let i = 0; i < n; i++) {
        if (!(N[i] > 0)) continue;
        for (let j = i; j < n; j++) {
          if (!(N[j] > 0)) continue;
          const q = i * n + j, r = (i === j ? 0.5 : 1) * B(i, j) * N[i] * N[j] * lim[i] * lim[j], k = pk[q];
          dN[i] -= r; dN[j] -= r; dN[k] += r * pa[q]; dN[k + 1] += r * pb[q];
        }
      }
      for (let i = 0; i < n; i++) N[i] = Math.max(N[i] + dt * dN[i], 0);
      if (S) { dN.fill(0); for (let i = 0; i < n; i++) { const k = bk[i]; if (k < 0 || !(S[i] > 0) || !(N[i] > 0)) continue; const r = N[i] * (1 - Math.exp(-S[i] * dt)); dN[i] -= r; dN[k] += r * ba[i]; dN[k + 1] += r * bb[i]; } for (let i = 0; i < n; i++) N[i] += dN[i]; }
    }
    if (src && src.rate > 0) N[src.k] += src.rate * dt;
  }
  if (gv) { // growth / shrinkage over the whole interval: every class is moved by its volume increment and shared between the neighbouring pivots
    tmp.fill(0);
    for (let i = 0; i < n; i++) {
      if (!(N[i] > 0)) continue;
      const dv = gv[i] * t, vn = v[i] + dv;
      if (dv === 0) { tmp[i] += N[i]; continue; }
      if (vn <= 0) continue; // dissolved completely
      if (vn <= v[0]) { tmp[0] += (N[i] * vn) / v[0]; continue; }
      if (vn >= v[n - 1]) { tmp[n - 1] += (N[i] * vn) / v[n - 1]; continue; }
      let k = dv > 0 ? i : 0; while (k < n - 2 && v[k + 1] <= vn) k++;
      const a = (v[k + 1] - vn) / (v[k + 1] - v[k]); tmp[k] += N[i] * a; tmp[k + 1] += N[i] * (1 - a);
    }
    N.set(tmp);
  }
  // drop numerically empty classes (below 1e-24 of the largest volume share): they cost denormal arithmetic and carry nothing
  let vmx = 0; for (let i = 0; i < n; i++) if (N[i] * v[i] > vmx) vmx = N[i] * v[i];
  for (let i = 0; i < n; i++) if (N[i] * v[i] < 1e-24 * vmx) N[i] = 0;
  return { N, sub, limited };
}
/** Moments of a sectional distribution: { m0, m1, m2, m3, vol (Σ N v), d10, d32, d43 }. */
export function pbeMoments(grid, N) {
  let m0 = 0, m1 = 0, m2 = 0, m3 = 0, m4 = 0, vol = 0;
  for (let i = 0; i < grid.n; i++) { const L = grid.L[i], w = N[i]; m0 += w; m1 += w * L; m2 += w * L * L; m3 += w * L ** 3; m4 += w * L ** 4; vol += w * grid.v[i]; }
  return { m0, m1, m2, m3, vol, d10: m0 > 0 ? m1 / m0 : 0, d32: m2 > 0 ? m3 / m2 : 0, d43: m3 > 0 ? m4 / m3 : 0 };
}
function symEig(A) { // cyclic Jacobi for a small symmetric matrix: { val[], vec[][] (columns) }
  const n = A.length, a = A.map((r) => r.slice()), V = a.map((_, i) => a.map((__, j) => +(i === j)));
  for (let sweep = 0; sweep < 60; sweep++) {
    let off = 0; for (let p = 0; p < n; p++) for (let q = p + 1; q < n; q++) off += a[p][q] ** 2;
    if (off < 1e-26) break;
    for (let p = 0; p < n; p++) for (let q = p + 1; q < n; q++) {
      if (Math.abs(a[p][q]) < 1e-300) continue;
      const th = (a[q][q] - a[p][p]) / (2 * a[p][q]), t = Math.sign(th || 1) / (Math.abs(th) + Math.sqrt(th * th + 1)), c = 1 / Math.sqrt(t * t + 1), s = t * c;
      for (let k = 0; k < n; k++) { const x = a[k][p], y = a[k][q]; a[k][p] = c * x - s * y; a[k][q] = s * x + c * y; }
      for (let k = 0; k < n; k++) { const x = a[p][k], y = a[q][k]; a[p][k] = c * x - s * y; a[q][k] = s * x + c * y; }
      for (let k = 0; k < n; k++) { const x = V[k][p], y = V[k][q]; V[k][p] = c * x - s * y; V[k][q] = s * x + c * y; }
    }
  }
  return { val: a.map((r, i) => r[i]), vec: V };
}
/**
 * Quadrature nodes and weights from the first 2N moments of a size distribution (Wheeler algorithm).
 * Returns { L[], w[], ok } with Σ w L^k = m_k for k = 0 … 2N−1.
 */
export function qmomNodes(m) {
  const N = m.length >> 1, sc = m[1] / m[0], mm = m.map((x, k) => x / (m[0] * sc ** k)), a = new Array(N).fill(0), b = new Array(N).fill(0);
  let prev = new Array(2 * N).fill(0), cur = mm.slice();
  a[0] = mm[1];
  for (let k = 1; k < N; k++) {
    const nx = new Array(2 * N).fill(0);
    for (let l = k; l < 2 * N - k; l++) nx[l] = cur[l + 1] - a[k - 1] * cur[l] - b[k - 1] * prev[l];
    a[k] = nx[k + 1] / nx[k] - cur[k] / cur[k - 1]; b[k] = nx[k] / cur[k - 1];
    prev = cur; cur = nx;
  }
  if (b.slice(1).some((x) => !(x > 0)) || !a.every(Number.isFinite)) return { L: [m[1] / m[0]], w: [m[0]], ok: false };
  const J = a.map((_, i) => a.map((__, j) => (i === j ? a[i] : Math.abs(i - j) === 1 ? Math.sqrt(b[Math.max(i, j)]) : 0))), e = symEig(J);
  const L = e.val.map((x) => x * sc), w = e.val.map((_, j) => m[0] * e.vec[0][j] ** 2), ok = L.every((x) => x > 0 && Number.isFinite(x)) && w.every((x) => x >= 0);
  return { L, w, ok };
}
/**
 * Quadrature method of moments for the same processes as the sectional solver (length-based moments m0 … m(2N−1)).
 * o: { G(L) → dL/dt, beta(Li, Lj) → m³/s, S(L) → 1/s (binary equal-volume breakage), J (1/s source), L0 (m), steps }.
 * Returns { m[], nodes: { L[], w[] }, ok }.
 */
export function solveQMOM(m0, t, { G: Gf = null, beta = null, S = null, J = 0, L0 = 1e-6, steps = 200 } = {}) {
  const K = m0.length;
  let nodes = qmomNodes(m0), ok = nodes.ok;
  const rhs = (m) => {
    const q = qmomNodes(m); if (q.ok) nodes = q; else ok = false;
    const { L, w } = nodes, d = new Array(K).fill(0);
    for (let k = 0; k < K; k++) {
      let s = J * L0 ** k;
      for (let i = 0; i < L.length; i++) {
        if (Gf && k > 0) s += k * w[i] * Gf(L[i]) * L[i] ** (k - 1);
        if (S) s += w[i] * S(L[i]) * L[i] ** k * (2 ** (1 - k / 3) - 1);
        if (beta) for (let j = 0; j < L.length; j++) s += 0.5 * w[i] * w[j] * beta(L[i], L[j]) * ((L[i] ** 3 + L[j] ** 3) ** (k / 3) - L[i] ** k - L[j] ** k);
      }
      d[k] = s;
    }
    return d;
  };
  let m = m0.slice(); const dt = t / steps;
  for (let s = 0; s < steps; s++) { const k1 = rhs(m), mp = m.map((x, i) => x + dt * k1[i]), k2 = rhs(mp); m = m.map((x, i) => x + 0.5 * dt * (k1[i] + k2[i])); }
  const q = qmomNodes(m);
  return { m, nodes: q.ok ? { L: q.L, w: q.w } : { L: nodes.L, w: nodes.w }, ok: ok && q.ok };
}

// =====================================================================================================
// 4a. Population balance on a resolved flow field: axisymmetric pipe flow (Reynolds-averaged momentum equation
//     with an algebraic eddy viscosity) carrying the size classes by advection and turbulent diffusion, with
//     aggregation and breakage at the local shear and wall deposition as a boundary flux
// =====================================================================================================
/** Constants of the eddy-viscosity closure (Cess composite of the van Driest wall damping and the Reichardt core profile). */
export const RANS = Object.freeze({ kappa: 0.42, Aplus: 27, ReLam: 2300, ScT: 0.9, Cmu: 0.09 }); // κ and A⁺: the pipe-flow fit of the Cess formula; k = νt |du/dr| / √Cμ
/**
 * Fully developed axisymmetric pipe flow. The Reynolds-averaged axial momentum balance τ(r) = τw·r/R = (μ + μt)(−du/dr) is
 * integrated from the wall on a fine wall-resolved grid with the total viscosity
 *   (ν + νt)/ν = ½ + ½ {1 + (κ² R⁺²/9)(1 − η²)²(1 + 2η²)² [1 − exp(−y⁺/A⁺)]²}^½,  η = r/R,
 * and the friction velocity is iterated until the bulk velocity is met. Below Re = 2300 the flow is laminar (νt = 0).
 * A radially varying molecular viscosity (slurry) may be given as muRel(η) → μ/μ0; wall units then use the local viscosity.
 * o: { D (m), U (m/s bulk), rho, mu (Pa·s), nr (transport cells, centre → wall), muRel, nFine }
 * Returns { R, uStar, tauW, f (Darcy), Re, laminar, rf[nr+1] (faces), rc[nr], dA[nr] (cell areas), u[nr], gamma[nr] (mean |du/dr|),
 *           epsT[nr] (turbulent dissipation = production, W/kg), nut[nr], k[nr], nutF[nr−1] (harmonic eddy viscosity between
 *           cell centres), fine: { y[], u[], nut[] } }.
 */
export function pipeFlowField({ D, U, rho, mu, nr = 16, muRel = null, nFine = 300, kappa = RANS.kappa, Aplus = RANS.Aplus }) {
  need(D > 0 && U > 0 && rho > 0 && mu > 0, 'pipeFlowField needs a positive diameter, velocity, density and viscosity.');
  const Rp = D / 2, Re = (rho * U * D) / mu, laminar = Re < RANS.ReLam, yl = [0];
  for (let dy = 2e-6 * Rp, dmax = Rp / Math.max(nFine / 2, 40); yl[yl.length - 1] < Rp * (1 - 1e-9); dy = Math.min(dy * 1.045, dmax)) yl.push(Math.min(yl[yl.length - 1] + dy, Rp)); // wall distance: geometric from 2e-6 R, then uniform
  const n = yl.length - 1, y = Float64Array.from(yl);
  const u = new Float64Array(n + 1), nut = new Float64Array(n + 1), muL = new Float64Array(n + 1);
  for (let j = 0; j <= n; j++) muL[j] = mu * (muRel ? Math.max(muRel(1 - y[j] / Rp), 1) : 1);
  const profile = (us) => { // returns the bulk velocity for a friction velocity us
    const tw = rho * us * us; let q = 0, prevS = 0;
    for (let j = 0; j <= n; j++) {
      const eta = 1 - y[j] / Rp, nu = muL[j] / rho, yp = (y[j] * us) / nu, Rpl = (Rp * us) / nu;
      const tot = laminar ? 1 : 0.5 + 0.5 * Math.sqrt(1 + ((kappa * kappa * Rpl * Rpl) / 9) * (1 - eta * eta) ** 2 * (1 + 2 * eta * eta) ** 2 * (1 - Math.exp(-yp / Aplus)) ** 2);
      nut[j] = nu * (tot - 1);
      const s = (tw * eta) / (muL[j] * tot); // du/dy
      if (j) { u[j] = u[j - 1] + 0.5 * (s + prevS) * (y[j] - y[j - 1]); const r0 = Rp - y[j - 1], r1 = Rp - y[j]; q += 0.5 * (u[j - 1] * r0 + u[j] * r1) * (y[j] - y[j - 1]); } else u[0] = 0;
      prevS = s;
    }
    return (2 * q) / (Rp * Rp);
  };
  let us;
  if (laminar && !muRel) us = Math.sqrt((8 * mu * U) / (rho * D)); // Hagen–Poiseuille
  else { const u0 = U * Math.sqrt(frictionFactor(Math.max(Re, 10), 0) / 8); us = brent((x) => profile(x) - U, 0.2 * u0, 5 * u0, 1e-10 * u0, 80); }
  profile(us);
  // transport grid: nr annular cells from the axis to the wall, refined toward the wall
  const rf = new Array(nr + 1), rc = new Array(nr), dA = new Array(nr), uc = new Array(nr).fill(0), gam = new Array(nr).fill(0), eps = new Array(nr).fill(0), nuc = new Array(nr).fill(0);
  // The wall cell is a wall-function cell: in turbulent flow it reaches at least y⁺ = 60 (centre at 30), because the deposition
  // velocity applied at the wall already contains the transport resistance of the viscous and buffer layers.
  const dW = laminar ? 0 : Math.min((60 * muL[0]) / (rho * us), 0.25 * Rp);
  if (dW > Rp * (1 / nr) ** 1.5 && nr > 2) { for (let j = 0; j < nr; j++) rf[j] = (Rp - dW) * (1 - (1 - j / (nr - 1)) ** 1.2); rf[nr] = Rp; }
  else for (let j = 0; j <= nr; j++) rf[j] = Rp * (1 - (1 - j / nr) ** 1.5);
  for (let j = 0; j < nr; j++) { rc[j] = 0.5 * (rf[j] + rf[j + 1]); dA[j] = PI * (rf[j + 1] ** 2 - rf[j] ** 2); }
  // cell averages by area-weighted integration of the fine profile (r = R − y)
  const tw = rho * us * us; let jc = nr - 1;
  for (let i = 1; i <= n; i++) {
    const r1 = Rp - y[i - 1], r0 = Rp - y[i], rm = 0.5 * (r0 + r1), um = 0.5 * (u[i] + u[i - 1]), nm = 0.5 * (nut[i] + nut[i - 1]), mm = 0.5 * (muL[i] + muL[i - 1]);
    const s = (tw * (rm / Rp)) / (mm + rho * nm);
    while (jc > 0 && r1 <= rf[jc]) jc--;
    for (let j = jc; j >= 0 && rf[j + 1] > r0; j--) { const a = Math.max(r0, rf[j]), b = Math.min(r1, rf[j + 1]); if (!(b > a)) continue; const w = PI * (b * b - a * a); uc[j] += um * w; gam[j] += s * w; eps[j] += nm * s * s * w; nuc[j] += nm * w; } // a fine annulus is shared between the cells it overlaps
  }
  for (let j = 0; j < nr; j++) { uc[j] /= dA[j]; gam[j] /= dA[j]; eps[j] /= dA[j]; nuc[j] /= dA[j]; }
  { let q = 0; for (let j = 0; j < nr; j++) q += uc[j] * dA[j]; const c = (U * PI * Rp * Rp) / q; for (let j = 0; j < nr; j++) uc[j] *= c; } // the cell velocities carry exactly the bulk flow
  const nutF = new Array(Math.max(nr - 1, 0)).fill(0);
  for (let j = 0; j < nr - 1; j++) { // harmonic mean between the cell centres: the resistances to radial mixing add up
    let res = 0, len = 0;
    for (let i = 1; i <= n; i++) { const r1 = Rp - y[i - 1], r0 = Rp - y[i]; if (r1 <= rc[j] || r0 >= rc[j + 1]) continue; const a = Math.max(r0, rc[j]), b = Math.min(r1, rc[j + 1]), nm = 0.5 * (nut[i] + nut[i - 1]); if (b > a) { res += (b - a) / Math.max(nm, 1e-30); len += b - a; } }
    nutF[j] = len > 0 && res > 0 ? len / res : 0;
  }
  return { R: Rp, D, U, rho, mu, uStar: us, tauW: tw, f: 8 * (us / U) ** 2, Re, laminar, nr, rf, rc, dA, u: uc, gamma: gam, epsT: eps, nut: nuc, k: nuc.map((v, j) => (v * gam[j]) / Math.sqrt(RANS.Cmu)), nutF,
    fine: { y: Array.from(y), u: Array.from(u), nut: Array.from(nut) } };
}
/**
 * Sectional population balance transported through the resolved pipe flow (steady, marching downstream):
 *   u(r) ∂N_c/∂x = (1/r) ∂/∂r [ r (D_B,c + νt/Sc_t) ∂N_c/∂r ] + aggregation and breakage at the local shear,
 *   wall:  −Γ ∂N_c/∂r = V_d,c · N_c  (deposition velocity of the class times the sticking probability).
 * Each step is an implicit radial diffusion solve per class (conservative finite volumes, upwind in x) followed by the
 * local kinetics over the time the fluid of that annulus needs to cross the step.
 * o: { field (pipeFlowField), grid (pbeGrid), Nin: number per m³ per class at the inlet (uniform) | [nr][nC], L (m), nx,
 *      kin: { alpha (collision efficiency), mu, rhoF, rhoP, TK, kBreak, dA (m, cohesive size limit) | dAof(shear) → m, dPrim, brownian (default true) },
 *      wall: { stick (0…1), adhForce (N/m; 0 = no shear limit on sticking), vd(L) → m/s (optional override) },
 *      ScT, mixing (multiplier on the radial diffusivity; a large value gives the well-mixed limit), maxSub }
 * Returns { x[nx+1], r[nr], N (outlet [nr][nC]), flux: { m0[], m3[] (flow-weighted moment fluxes along x) }, dep[nx] (volume
 *   deposition flux m³/m²/s), phi[nx+1][nr], d43[nx+1][nr], mixed: { N[nC] (flow-weighted outlet), m0, m3, d43, d32 },
 *   ledger: { in, out, deposited (m³ of particles per second), numIn, numOut, numDeposited }, vdBulk (m/s, volume-weighted) }.
 */
export function solvePBEField({ field: F, grid: g, Nin, L, nx = 20, kin = {}, wall = {}, ScT = RANS.ScT, mixing = 1, maxSub = 40 }) {
  const nr = F.nr, nC = g.n, dx = L / nx, nu = (kin.mu ?? F.mu) / (kin.rhoF ?? F.rho), mu = kin.mu ?? F.mu, TK = kin.TK ?? 277, rhoF = kin.rhoF ?? F.rho, rhoP = kin.rhoP ?? HYDRATE.rho, alpha = kin.alpha ?? 0;
  let N = Array.from({ length: nr }, (_, j) => Float64Array.from(Array.isArray(Nin[0]) || ArrayBuffer.isView(Nin[0]) ? Nin[j] : Nin));
  const stick = clamp(wall.stick ?? 0, 0, 1), tw = F.tauW, vd = new Float64Array(nC), dB = new Float64Array(nC);
  for (let c = 0; c < nC; c++) { const Lc = g.L[c]; dB[c] = kin.brownian === false ? 0 : (KB * TK) / (3 * PI * mu * Lc); vd[c] = stick > 0 ? (wall.vd ? wall.vd(Lc) : depositionVelocity(Lc, rhoP, F.uStar, nu, rhoF, TK)) * stick * (wall.adhForce > 0 ? Math.min(1, wall.adhForce / (8 * Math.max(tw, 1e-9) * Lc)) : 1) : 0; }
  const a = new Array(nr), b = new Array(nr), cc = new Array(nr), d = new Array(nr), Sb = new Float64Array(nC), beta = new Float64Array(3), L3 = g.L.map((x) => x ** 3), L4 = g.L.map((x) => x ** 4);
  const cD = (alpha * (PI / 4) * Math.abs(rhoP - rhoF) * G) / (18 * mu), cB = kin.brownian === false ? 0 : (alpha * 2 * KB * TK) / (3 * mu), wallP = 2 * PI * F.R;
  const kernels = F.u.map((_, j) => { const gd = Math.max(F.gamma[j], 1e-9), gT = Math.sqrt(Math.max(F.epsT[j], 0) / nu), dAj = kin.dAof ? kin.dAof(gd + gT) : kin.dA; return { cS: alpha * (gd / 6 + 0.1618 * gT), kb: kin.kBreak > 0 && dAj > 0 ? (kin.kBreak * (gd + gT)) / dAj ** 3 : 0, dA: dAj || 0, shear: gd + gT }; });
  const fluxOf = (w) => { let s = 0; for (let j = 0; j < nr; j++) { let q = 0; for (let c = 0; c < nC; c++) q += N[j][c] * w[c]; s += q * F.u[j] * F.dA[j]; } return s; };
  const one = new Array(nC).fill(1), x = [0], fm0 = [fluxOf(one)], fm3 = [fluxOf(g.v)], dep = [], phi = [], d43 = [];
  const snap = () => { const p = new Array(nr), q = new Array(nr); for (let j = 0; j < nr; j++) { let v = 0, m3 = 0, m4 = 0; for (let c = 0; c < nC; c++) { v += N[j][c] * g.v[c]; m3 += N[j][c] * L3[c]; m4 += N[j][c] * L4[c]; } p[j] = v; q[j] = m3 > 0 ? m4 / m3 : 0; } phi.push(p); d43.push(q); };
  snap();
  let depV = 0, depN = 0, limited = false;
  for (let s = 0; s < nx; s++) {
    let depStep = 0;
    for (let c = 0; c < nC; c++) { // radial transport of one class, implicit
      let any = false; for (let j = 0; j < nr; j++) if (N[j][c] > 0) { any = true; break; }
      if (!any) continue;
      for (let j = 0; j < nr; j++) {
        const cap = (F.u[j] * F.dA[j]) / dx, gm = j > 0 ? (2 * PI * F.rf[j] * ((mixing * F.nutF[j - 1]) / ScT + dB[c])) / (F.rc[j] - F.rc[j - 1]) : 0, gp = j < nr - 1 ? (2 * PI * F.rf[j + 1] * ((mixing * F.nutF[j]) / ScT + dB[c])) / (F.rc[j + 1] - F.rc[j]) : 0;
        a[j] = -gm; cc[j] = -gp; b[j] = cap + gm + gp + (j === nr - 1 ? wallP * vd[c] : 0); d[j] = cap * N[j][c];
      }
      const sol = tridiag(a, b, cc, d);
      for (let j = 0; j < nr; j++) N[j][c] = Math.max(sol[j], 0);
      const lost = wallP * vd[c] * N[nr - 1][c]; depStep += lost * g.v[c]; depN += lost * dx;
    }
    depV += depStep * dx; dep.push(depStep / wallP);
    if (alpha > 0 || kin.kBreak > 0) for (let j = 0; j < nr; j++) { // local kinetics along the streamline of the annulus
      const K = kernels[j]; let any = false; for (let c = 0; c < nC; c++) if (N[j][c] > 0) { any = true; break; }
      if (!any) continue;
      beta[0] = K.cS; beta[1] = cD; beta[2] = cB;
      const lim = Math.max(kin.dPrim || 0, 0) * 1.01; for (let c = 0; c < nC; c++) Sb[c] = K.kb > 0 && g.L[c] > lim ? K.kb * L3[c] : 0;
      const r = solvePBE(g, N[j], dx / Math.max(F.u[j], 1e-9), { beta: alpha > 0 ? beta : null, S: K.kb > 0 ? Sb : null, maxSub, frac: 0.25 });
      if (r.limited) limited = true; N[j] = r.N;
    }
    x.push((s + 1) * dx); fm0.push(fluxOf(one)); fm3.push(fluxOf(g.v)); snap();
  }
  const Q = F.U * PI * F.R * F.R, Nm = new Array(nC).fill(0);
  for (let j = 0; j < nr; j++) for (let c = 0; c < nC; c++) Nm[c] += (N[j][c] * F.u[j] * F.dA[j]) / Q;
  const mom = pbeMoments(g, Nm); let vdB = 0, vt = 0; for (let c = 0; c < nC; c++) { vdB += vd[c] * Nm[c] * g.v[c]; vt += Nm[c] * g.v[c]; }
  return { x, r: F.rc.slice(), N: N.map((q) => Array.from(q)), flux: { m0: fm0, m3: fm3 }, dep, phi, d43, kernels, limited, mixed: { N: Nm, m0: mom.m0, m3: mom.m3, vol: mom.vol, d43: mom.d43, d32: mom.d32 },
    ledger: { in: fm3[0], out: fm3[nx], deposited: depV, numIn: fm0[0], numOut: fm0[nx], numDeposited: depN }, vdBulk: vt > 0 ? vdB / vt : 0 };
}
// =====================================================================================================
// 5. Cohesion, slurry rheology, settling, particle momentum, wall deposition, porous plug
// =====================================================================================================
/**
 * Largest stable agglomerate from the Camargo–Palermo balance between the cohesive force and the shear stress.
 * { dp (m, primary particle), Fa (N, cohesive force between two primaries), mu0 (Pa·s), shear (1/s), phi (hydrate volume fraction), phiMax, fr (fractal dimension) }
 * Returns { ratio (dA/dp ≥ 1), dA (m), phiEff (effective volume fraction of the porous agglomerates) }.
 */
export function maxAgglomerateSize({ dp, Fa, mu0, shear, phi, phiMax = 4 / 7, fr = 2.5 }) {
  const p = clamp(phi, 0, 0.99 * phiMax), e = 3 - fr, cap = p > 1e-9 ? (phiMax / p) ** (1 / e) : 1e6;
  const g = (x) => x ** (4 - fr) - (Fa * (1 - (p / phiMax) * x ** e) ** 2) / (dp * dp * mu0 * Math.max(shear, 1e-9) * (1 - p * x ** e));
  let ratio = 1;
  if (g(1) < 0) { const hi = Math.min(cap * (1 - 1e-9), 1e6); ratio = g(hi) <= 0 ? hi : brent(g, 1, hi, 1e-6, 60); }
  return { ratio, dA: ratio * dp, phiEff: Math.min(p * ratio ** e, phiMax) };
}
/**
 * Relative viscosity of a suspension. model: 'mills' | 'krieger' (Krieger–Dougherty) | 'thomas' | 'einstein'.
 * phi is the effective volume fraction (agglomerates count with their trapped liquid). Capped at 1e4.
 */
export function slurryViscosity(phi, model = 'mills', { phiMax = 4 / 7, intrinsic = 2.5 } = {}) {
  const p = Math.max(phi, 0);
  if (model === 'einstein') return 1 + intrinsic * p;
  if (model === 'thomas') return Math.min(1 + 2.5 * p + 10.05 * p * p + 0.00273 * Math.exp(16.6 * p), 1e4);
  if (p >= phiMax * 0.9999) return 1e4;
  return Math.min(model === 'krieger' ? (1 - p / phiMax) ** (-intrinsic * phiMax) : (1 - p) / (1 - p / phiMax) ** 2, 1e4);
}
/**
 * Terminal settling velocity of a sphere. model: 'stokes' | 'schiller' (Schiller–Naumann drag, iterated).
 * Returns { v (m/s, positive = sinks, negative = rises), Re, Cd, n (Richardson–Zaki exponent), vHindered (at volume fraction phi) }.
 */
export function settlingVelocity(d, rhoP, rhoF, mu, { model = 'schiller', phi = 0, shape = 1 } = {}) {
  const dr = rhoP - rhoF, sgn = Math.sign(dr), vSt = (Math.abs(dr) * G * d * d) / (18 * mu);
  let v = vSt, Re = (rhoF * v * d) / mu, Cd = Re > 0 ? 24 / Re : Infinity;
  if (model !== 'stokes' && vSt > 0) for (let it = 0; it < 80; it++) {
    Re = (rhoF * v * d) / mu; Cd = shape * (Re < 1000 ? (24 / Re) * (1 + 0.15 * Re ** 0.687) : 0.44);
    const vn = Math.sqrt((4 * Math.abs(dr) * G * d) / (3 * Cd * rhoF));
    if (Math.abs(vn - v) < 1e-12 * v) { v = vn; break; } v = 0.5 * (v + vn);
  }
  Re = (rhoF * v * d) / mu;
  const n = Re < 0.2 ? 4.65 : Re < 1 ? 4.4 * Re ** -0.03 : Re < 500 ? 4.4 * Re ** -0.1 : 2.4;
  return { v: sgn * v, Re, Cd, n, vHindered: sgn * v * (1 - clamp(phi, 0, 0.99)) ** n };
}
/**
 * Particle momentum equation in a quiescent fluid: (ρp + ½ρf) dv/dt = (ρp − ρf) g − ¾ Cd ρf |v| v / d (drag by Schiller–Naumann,
 * buoyancy and added mass; history force neglected). Returns { t[], v[], vTerminal, tau (s, 63 % response time) }.
 */
export function particleRelaxation({ d, rhoP, rhoF, mu, tEnd = null }) {
  const vt = settlingVelocity(d, rhoP, rhoF, mu).v, tauS = ((rhoP + 0.5 * rhoF) * d * d) / (18 * mu), T = tEnd ?? 8 * tauS;
  const f = (t, y) => { const v = y[0], Re = Math.max((rhoF * Math.abs(v) * d) / mu, 1e-12), Cd = Re < 1000 ? (24 / Re) * (1 + 0.15 * Re ** 0.687) : 0.44; return [((rhoP - rhoF) * G - (0.75 * Cd * rhoF * Math.abs(v) * v) / d) / (rhoP + 0.5 * rhoF)]; };
  const r = rk45(f, [0], 0, T, { rtol: 1e-8, atol: 1e-14 }), v = r.y.map((y) => y[0]);
  let tau = T; for (let i = 1; i < v.length; i++) if (Math.abs(v[i]) >= 0.6321 * Math.abs(vt)) { const a = Math.abs(v[i - 1]), b = Math.abs(v[i]); tau = r.t[i - 1] + ((0.6321 * Math.abs(vt) - a) / (b - a || 1)) * (r.t[i] - r.t[i - 1]); break; }
  return { t: r.t, v, vTerminal: vt, tau };
}
/**
 * Turbulent deposition velocity of a particle onto a pipe wall (m/s): diffusion + eddy-impaction regimes capped at the
 * inertia-moderated plateau. V⁺ = min(0.13, 0.057 Sc^−2/3 + 4.5e-4 τ⁺²).
 */
export function depositionVelocity(d, rhoP, uStar, nu, rhoF, TK = 277) {
  if (!(uStar > 0)) return 0;
  const mu = nu * rhoF, tauP = (rhoP * d * d * uStar * uStar) / (18 * mu * nu), Sc = nu / (KB * TK / (3 * PI * mu * d));
  return uStar * Math.min(0.13, 0.057 * Sc ** (-2 / 3) + 4.5e-4 * tauP * tauP);
}
/** Kozeny–Carman permeability (m²) of a packed bed of porosity eps and grain size dp. */
export const kozenyCarman = (eps, dp) => (eps ** 3 * dp * dp) / (180 * (1 - eps) ** 2);
/**
 * Pressure gradient (Pa/m) for a superficial velocity v through a porous plug: Darcy term μv/k plus the Forchheimer
 * (Ergun inertial) term β ρ v². Returns { k (m²), beta (1/m), darcy, forchheimer, total }.
 */
export function porousGradient(v, mu, rho, eps, dp) {
  const k = kozenyCarman(eps, dp), beta = (1.75 * (1 - eps)) / (eps ** 3 * dp), darcy = (mu * v) / k, forchheimer = beta * rho * v * v;
  return { k, beta, darcy, forchheimer, total: darcy + forchheimer };
}
/** Superficial velocity (m/s) through a porous plug under a pressure gradient (Pa/m): root of the Darcy–Forchheimer quadratic. */
export function porousVelocity(dpdx, mu, rho, eps, dp) { const k = kozenyCarman(eps, dp), b = (1.75 * (1 - eps)) / (eps ** 3 * dp) * rho, a = mu / k; return b > 0 ? (-a + Math.sqrt(a * a + 4 * b * dpdx)) / (2 * b) : dpdx / a; }

// =====================================================================================================
// 6. Dissociation of a plug: Stefan moving boundary
// =====================================================================================================
/** Similarity constant λ of the one-phase Stefan problem: λ·exp(λ²)·erf(λ) = Ste/√π. */
export function stefanLambda(Ste) {
  const erf = (x) => { const t = 1 / (1 + 0.3275911 * x), y = 1 - ((((1.061405429 * t - 1.453152027) * t + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * Math.exp(-x * x); return y; };
  return brent((l) => l * Math.exp(l * l) * erf(l) - Ste / Math.sqrt(PI), 1e-9, 5, 1e-12);
}
/**
 * Radial melting of a hydrate plug that fills the bore after depressurisation (heat-transfer-controlled, quasi-steady
 * Stefan problem): heat flows from the surroundings through the overall coefficient U (on the bore radius) and the
 * annulus of released water to the dissociation front at Td.
 * { R (m, bore radius), U (W/m²/K), kW (W/m/K of the melt annulus), Tamb, Td (°C), eps (plug porosity), latent, rho, steps }
 * Returns { tAnalytic, tNumeric (s; Infinity when Tamb ≤ Td), r: [], t: [] (front radius against time) }.
 */
export function plugMeltTime({ R: Rb, U, kW = 0.58, Tamb, Td, eps = 0.4, latent = HYDRATE.latent, rho = HYDRATE.rho, steps = 400 }) {
  const dT = Tamb - Td, q = rho * (1 - eps) * latent;
  if (!(dT > 0)) return { tAnalytic: Infinity, tNumeric: Infinity, r: [Rb], t: [0] };
  const tAnalytic = (q * (Rb / (2 * U) + (Rb * Rb) / (4 * kW))) / dT, r = [Rb], t = [0];
  let tt = 0;
  for (let i = 0; i < steps; i++) { // front tracking: time for the front to move one radial increment at the local heat flux
    const s1 = Rb * (1 - i / steps), s2 = Rb * (1 - (i + 1) / steps), s = 0.5 * (s1 + s2), flux = dT / (1 / (U * Rb) + Math.log(Rb / s) / kW); // W/m per 2π
    tt += (q * s * (s1 - s2)) / flux; r.push(s2); t.push(tt);
  }
  return { tAnalytic, tNumeric: tt, r, t };
}

// =====================================================================================================
// 7. Wax
// =====================================================================================================
/**
 * Diffusivity of wax molecules in oil (m²/s). model: 'haydukMinhas' | 'wilkeChang'. mu in Pa·s, VA molar volume of the
 * wax (cm³/mol), MB solvent molar mass (g/mol).
 */
export function waxDiffusivity(TK, mu, { model = 'haydukMinhas', VA = 430, MB = 200, assoc = 1 } = {}) {
  const cP = mu * 1000;
  return model === 'wilkeChang' ? (7.4e-12 * Math.sqrt(assoc * MB) * TK) / (cP * VA ** 0.6) : (13.3e-12 * TK ** 1.47 * cP ** (10.2 / VA - 0.791)) / VA ** 0.71;
}
/**
 * Wax solubility curve: dissolved mass fraction falls exponentially below the wax appearance temperature.
 * Returns { dissolved, solid (mass fractions of the oil), dCdT (1/K) }.
 */
export function waxSolubility(T, wat, wTot, slope = 0.04, curve = null) {
  if (curve && curve.wTot > 0) { const q = curve.at(T), sc = wTot / curve.wTot; return { dissolved: q.dissolved * sc, solid: q.solid * sc, dCdT: q.dCdT * sc }; } // tabulated solid–liquid equilibrium, scaled with the wax still in the oil
  if (T >= wat) return { dissolved: wTot, solid: 0, dCdT: 0 };
  const c = wTot * Math.exp(-slope * (wat - T));
  return { dissolved: c, solid: wTot - c, dCdT: slope * c };
}
/**
 * Local wax deposition rates on a cold wall.
 * o: { Tb, Tamb (°C), U (W/m²/K, clean overall coefficient on the bore), hIn (W/m²/K), kOil, rhoOil, muOil, wat, wTot (mass fraction),
 *      slope (1/K), delta (m, present thickness), Fw (wax fraction of the deposit), kDep (W/m/K), D (m), vL (m/s), gammaW (1/s wall shear rate),
 *      rhoMix, regime, mult (deposition multiplier), diffModel, wetFrac }
 *      wmodel: 'matzain' (default: Fick flux with the empirical enhancement Π1 = C1/(1 − C_oil), C1 = 15, the deposit oil fraction
 *      C_oil from the Reynolds-number closure unless `coil` is given, and shear stripping) | 'mechanistic' (Fick + shear dispersion +
 *      Brownian crystal flux with ageing by hindered diffusion into the gel; needs a fitted multiplier), muBulk (Pa·s, for N_Re) }
 * Returns { Ti (deposit surface °C), q (W/m²), dTdr (K/m), Dwo, jMol, jShear, jBrown (kg wax/m²/s), strip (shear-stripping factor 0–1),
 *           dDelta (m/s), dFw (1/s), Fset (wax fraction of the deposit set by the closure, Matzain model only), Ueff }.
 */
export function waxDeposition({ Tb, Tamb, U, hIn, kOil = 0.13, rhoOil = 800, muOil = 3e-3, wat, wTot, slope = 0.04, delta = 0, Fw = 0.2, kDep = 0.25, D = 0.25, vL = 1, gammaW = 100, rhoMix = null, regime = '', mult = 1, diffModel = 'haydukMinhas', wetFrac = 1, rhoWax = 900, curve = null, stripC = 0.055, stripN = 1.4, aspect = 8, wmodel = 'matzain', c1 = 15, coil = null, muBulk = null }) {
  const Ueff = 1 / (1 / U + delta / kDep), q = Ueff * (Tb - Tamb), Ti = Tb - q / Math.max(hIn, 1e-6), dTdr = q / kOil;
  const z = { Ti, q, dTdr, Dwo: 0, jMol: 0, jShear: 0, jBrown: 0, strip: 1, dDelta: 0, dFw: 0, Ueff };
  if (!(Ti < wat) || !(q > 0) || !(wTot > 0)) return z;
  const TK = Ti + KEL, Dwo = waxDiffusivity(TK, muOil, { model: diffModel }), sol = waxSolubility(Ti, wat, wTot, slope, curve), bulk = waxSolubility(Tb, wat, wTot, slope, curve);
  const jMol = rhoOil * Dwo * sol.dCdT * dTdr; // Fick's law with the solubility slope and the radial temperature gradient
  const dCr = 10e-6, sub = (5 * muOil) / (rhoOil * Math.max(Math.sqrt((gammaW * muOil) / rhoOil), 1e-6)); // crystal size, viscous sub-layer thickness
  const jShear = bulk.solid > 0 ? (rhoOil * 0.1 * (dCr / 2) ** 2 * gammaW * bulk.solid * bulk.solid) / (D / 2) : 0; // shear dispersion of precipitated crystals (Eckstein diffusivity, gradient over the radius)
  const jBrown = bulk.solid > 0 ? (rhoOil * ((KB * TK) / (3 * PI * muOil * dCr)) * bulk.solid) / sub : 0; // Brownian diffusion of crystals to the wall
  const nsr = ((/slug|bubble|churn/.test(regime) ? rhoMix ?? rhoOil : /annular/.test(regime) ? Math.sqrt((rhoMix ?? rhoOil) * rhoOil) : rhoOil) * vL * Math.max(delta, 1e-5)) / muOil;
  const strip = 1 / (1 + stripC * (nsr / 1000) ** stripN); // shear stripping (Matzain form on a film Reynolds number in thousands)
  if (wmodel === 'matzain') { // dδ/dt = Π1/(1 + Π2) · D_wo · (dw/dT · dT/dr), Π1 = C1/(1 − C_oil), Π2 = C2 N_SR^C3, C_oil = 1 − N_Re^0.15/8
    const NRe = (rhoOil * vL * D) / (muBulk ?? muOil), Coil = coil !== null && coil > 0 ? clamp(coil, 0.05, 0.97) : clamp(1 - NRe ** 0.15 / 8, 0.05, 0.97), pi1 = c1 / (1 - Coil);
    return { Ti, q, dTdr, Dwo, jMol, jShear: 0, jBrown: 0, strip, dDelta: mult * wetFrac * pi1 * strip * Dwo * sol.dCdT * dTdr, dFw: 0, Fset: 1 - Coil, Coil, pi1, NRe, Ueff };
  }
  const j = mult * wetFrac * (jMol + jShear + jBrown), F = clamp(Fw, 0.02, 0.98), De = 1 / (1 + (aspect * aspect * F * F) / (1 - F)), psi = De * (1 - F); // ageing: Cussler hindered diffusion into the gel (crystal aspect ratio α)
  const dDelta = ((j * (1 - psi)) / (rhoWax * F)) * strip, dFw = (j * psi) / (rhoWax * Math.max(delta, 2e-5));
  return { Ti, q, dTdr, Dwo, jMol, jShear, jBrown, strip, dDelta, dFw, Ueff };
}

// =====================================================================================================
// 8. Mineral scale: aqueous speciation with Pitzer or ion-association activity models (USGS PHREEQC data)
// =====================================================================================================
const IONS = { Na: { z: 1, M: 22.9898 }, K: { z: 1, M: 39.0983 }, Ca: { z: 2, M: 40.08 }, Mg: { z: 2, M: 24.305 }, Ba: { z: 2, M: 137.33 }, Sr: { z: 2, M: 87.62 }, Fe: { z: 2, M: 55.847 }, Cl: { z: -1, M: 35.453 }, SO4: { z: -2, M: 96.064 }, HCO3: { z: -1, M: 61.017 } };
export const ION_IDS = Object.freeze(Object.keys(IONS));
/** Standard seawater (mg/L) used for the injection-water mixing curve. */
export const SEAWATER = Object.freeze({ Na: 10781, K: 399, Ca: 412, Mg: 1284, Ba: 0.02, Sr: 7.9, Fe: 0.003, Cl: 19353, SO4: 2712, HCO3: 142 });
// Minerals: molar mass (g/mol) and density (kg/m³) for the deposit volume; solubility products come from the PHREEQC tables.
const MINERALS = [
  { id: 'calcite', key: 'Calcite', name: 'Calcite (CaCO₃)', M: 100.09, rho: 2710, cat: 'Ca', an: 'CO3' },
  { id: 'barite', key: 'Barite', name: 'Barite (BaSO₄)', M: 233.39, rho: 4480, cat: 'Ba', an: 'SO4' },
  { id: 'celestite', key: 'Celestite', name: 'Celestite (SrSO₄)', M: 183.68, rho: 3960, cat: 'Sr', an: 'SO4' },
  { id: 'gypsum', key: 'Gypsum', name: 'Gypsum (CaSO₄·2H₂O)', M: 172.17, rho: 2320, cat: 'Ca', an: 'SO4', nW: 2 },
  { id: 'anhydrite', key: 'Anhydrite', name: 'Anhydrite (CaSO₄)', M: 136.14, rho: 2960, cat: 'Ca', an: 'SO4' },
  { id: 'siderite', key: 'Siderite', name: 'Siderite (FeCO₃)', M: 115.85, rho: 3870, cat: 'Fe', an: 'CO3' },
];
const LN10 = Math.LN10, R_ATM = 82.0597; // cm³·atm/(mol·K), the value used by PHREEQC

/**
 * Pure-water density, compressibility, relative dielectric constant and the Debye–Hückel slopes at Tc (°C), Pbar (bara),
 * ported from PHREEQC (utilities.cpp: Wagner–Pruss saturation density with a pressure polynomial; Bradley–Pitzer dielectric
 * constant). Returns { rho (g/cm³), kappa (1/atm), eps, A (log10 slope), B (1/Å), Aphi (Pitzer), Av (cm³/mol per √molal), QBrn }.
 */
const wMemo = new Map();
export function waterDH(Tc, Pbar = 1.01325) {
  const key = Math.round(Tc * 100) * 1e6 + Math.round(Pbar * 10); let o = wMemo.get(key); if (o) return o;
  const tc = Math.min(Tc, 350), T = tc + KEL, th = 1 - T / 647.096;
  const rs = 322 * (1 + 1.99274064 * th ** (1 / 3) + 1.09965342 * th ** (2 / 3) - 0.510839303 * th ** (5 / 3) - 1.75493479 * th ** (16 / 3) - 45.5170352 * th ** (43 / 3) - 6.7469445e5 * th ** (110 / 3));
  const p0 = 5.188e-2 + tc * (-4.1885519e-4 + tc * (6.6780748e-6 + tc * (-3.6648699e-8 + tc * 8.3501912e-11))), p1 = -6.0251348e-6 + tc * (3.6696407e-7 + tc * (-9.2056269e-9 + tc * (6.7024182e-11 + tc * -1.5947241e-13)));
  const p2 = -2.2983596e-9 + tc * (-4.0133819e-10 + tc * (1.2619821e-11 + tc * (-9.8952363e-14 + tc * 2.3363281e-16))), p3 = 7.0517647e-11 + tc * (6.8566831e-12 + tc * (-2.282975e-13 + tc * (1.8113313e-15 + tc * -4.2475324e-18)));
  const psat = Math.exp(11.6702 - 3816.44 / (T - 46.13)), pa = Math.max(Pbar / 1.01325, psat), pe = pa - (psat - 1e-6);
  const rho0 = Math.max(rs + pe * (p0 + pe * (p1 + pe * (p2 + Math.sqrt(pe) * p3))), 0.01), kappa = (p0 + pe * (2 * p1 + pe * (3 * p2 + Math.sqrt(pe) * 3.5 * p3))) / rho0, rho = rho0 / 1e3;
  const d1000 = 3.4279e2 * Math.exp(T * (-5.0866e-3 + T * 9.469e-7)), c = -2.0525 + 3.1159e3 / (-1.8289e2 + T), b = -8.0325e3 + 4.2142e6 / T + 2.1417 * T, pb = pa * 1.01325;
  const eps = Math.max(d1000 + c * Math.log((b + pb) / (b + 1e3)), 10), e2 = 1.671008e-3 / (eps * T), Bcm = Math.sqrt((8 * PI * 6.02252e23 * e2 * rho) / 1e3);
  o = { rho, kappa, eps, A: (Bcm * e2) / (2 * LN10), B: Bcm / 1e8, Aphi: (Bcm * e2) / 6, Av: Bcm * e2 * 0.0820597 * 1e3 * T * (((c / (b + pb)) * 1.01325) / eps - kappa / 3), QBrn: (c / (b + pb) / eps / eps) * 41.84004, T, pb };
  if (wMemo.size > 4000) wMemo.clear(); wMemo.set(key, o); return o;
}
/** log10 K(T) from a PHREEQC entry: the analytic expression when present, else van 't Hoff with ΔH, else the 25 °C value. */
const logKT = (e, T) => (e.an ? e.an[0] + (e.an[1] || 0) * T + (e.an[2] || 0) / T + (e.an[3] || 0) * Math.log10(T) + (e.an[4] || 0) / (T * T) + (e.an[5] || 0) * T * T : e.logK - ((e.dH || 0) / (LN10 * R)) * (1 / T - 1 / 298.15));
/** Molar volume (cm³/mol) of an aqueous species from its PHREEQC -Vm parameters at the state W (waterDH) and ionic strength I. */
function speciesVm(vm, z, W, I) {
  if (!vm || !vm.length) return 0;
  const pbs = 2600 + W.pb, Ts = W.T - 228, sq = Math.sqrt(I), a0 = vm[5] || 0;
  let v = 41.84004 * (0.1 * vm[0] + (100 * (vm[1] || 0)) / pbs + (vm[2] || 0) / Ts + (1e4 * (vm[3] || 0)) / (pbs * Ts)) - (vm[4] || 0) * 1e5 * W.QBrn;
  if (z) { v += (z * z * 0.5 * W.Av * sq) / (a0 < 1e-5 ? 1 : 1 + a0 * W.B * sq); const bi = (vm[6] || 0) + (vm[7] || 0) / Ts + (vm[8] || 0) * Ts; if (bi && I > 0) v += bi * I ** (vm[9] ?? 1); }
  return v;
}
const spId = (name) => name.replace(/[-+]\d*$/, ''), spZ = (name) => { const m = /([-+])(\d*)$/.exec(name); return m ? (m[1] === '+' ? 1 : -1) * (m[2] ? +m[2] : 1) : 0; };
/** Parse 'a A + b B = C' into { product, z, nu: { id: coef } } (reactant side only; H2O is kept under the id H2O). */
function parseReaction(eq) {
  const [lhs, rhs] = eq.split(' = '), nu = {};
  for (const term of lhs.split(' + ')) { const m = /^(\d+(?:\.\d+)?)\s+(.+)$/.exec(term.trim()), nm = m ? m[2] : term.trim(); nu[spId(nm)] = (nu[spId(nm)] || 0) + (m ? +m[1] : 1); }
  return { product: spId(rhs.trim()), z: spZ(rhs.trim()), nu };
}
const MASTERS = ['Na', 'K', 'Ca', 'Mg', 'Ba', 'Sr', 'Fe', 'Cl', 'SO4'], CARB = ['H', 'OH', 'HCO3', 'CO3', 'CO2'];
const chemCache = {};
/** Static structure of an aqueous model ('pitzer' | 'truesdellJones' | 'davies'): species, charges, volume and activity parameters, complexes. */
function chemModel(model) {
  if (chemCache[model]) return chemCache[model];
  const pz = model === 'pitzer', db = pz ? PHREEQC.pitzer : PHREEQC.dh, names = [...MASTERS, ...CARB], z = [1, 1, 2, 2, 2, 2, 2, -1, -2, 1, -1, -1, -2, 0], vm = [], gam = [];
  const full = { Na: 'Na+', K: 'K+', Ca: 'Ca+2', Mg: 'Mg+2', Ba: 'Ba+2', Sr: 'Sr+2', Fe: 'Fe+2', Cl: 'Cl-', SO4: 'SO4-2', H: 'H+', CO3: 'CO3-2' }, rx = db.reactions;
  const kW = rx['H2O = OH- + H+'], kB = rx['CO3-2 + H+ = HCO3-'], kA = rx['CO3-2 + 2 H+ = CO2 + H2O'];
  for (const n of names) { const e = n === 'OH' ? kW : n === 'HCO3' ? kB : n === 'CO2' ? kA : db.species[full[n]] || {}; vm.push(e.vm || null); gam.push(model === 'davies' ? null : e.gamma || null); }
  const cx = [];
  for (const [eq, e] of Object.entries(rx)) {
    if (/^H2O =|= HCO3-$|= CO2 \+ H2O$|^2 CO2 =/.test(eq)) continue;
    const r = parseReaction(eq), idx = [], coef = []; let alk = 0;
    for (const [id, c] of Object.entries(r.nu)) { idx.push(names.indexOf(id)); coef.push(c); alk += id === 'CO3' ? 2 * c : id === 'HCO3' ? c : id === 'H' ? -c : 0; }
    cx.push({ name: r.product, z: r.z, idx, coef, alk, e, gamma: model === 'davies' ? null : e.gamma || null });
  }
  const o = { model, pz, names, z, vm, gam, cx, n: names.length, kW, kB, kA, kH: db.phases['CO2(g)'], phases: { ...PHREEQC.dh.phases, ...(pz ? PHREEQC.pitzer.phases : {}) } };
  if (pz) { // Pitzer parameter lists with species indices; alpha values by charge type as in PHREEQC (pitzer.cpp)
    const ix = (s) => names.indexOf(spId(s)), P = db.params, pair = (list, type) => list.map((q) => { const i = ix(q[0]), j = ix(q[1]), zi = Math.abs(z[i]), zj = Math.abs(z[j]), ord = zi === 1 || zj === 1 ? 1 : zi === 2 && zj === 2 ? 2 : 3; return { i, j, a: q[2], alpha: type === 'B1' ? (ord === 2 ? 1.4 : 2) : type === 'B2' ? (ord === 3 ? 50 : 12) : 0 }; });
    o.B0 = pair(P.B0, 'B0'); o.B1 = pair(P.B1, 'B1'); o.B2 = pair(P.B2, 'B2'); o.C0 = pair(P.C0, 'C0'); o.TH = pair(P.THETA, 'TH'); o.LA = pair(P.LAMBDA, 'LA');
    o.PSI = [...P.PSI, ...P.ZETA].map((q) => ({ i: ix(q[0]), j: ix(q[1]), k: ix(q[2]), a: q[3] }));
    o.like = []; for (let i = 0; i < names.length; i++) for (let j = i + 1; j < names.length; j++) if (z[i] * z[j] > 0 && z[i] !== z[j]) o.like.push([i, j]); // unsymmetrical like-charge pairs (higher-order electrostatic terms)
    const kcl = (L) => L.find((q) => (names[q.i] === 'K' && names[q.j] === 'Cl') || (names[q.i] === 'Cl' && names[q.j] === 'K'));
    o.mac = [kcl(o.B0), kcl(o.B1), kcl(o.C0)]; o.iCl = names.indexOf('Cl');
  }
  return (chemCache[model] = o);
}
const pzT = (a, T) => { const Tr = 298.15; return Math.abs(T - Tr) < 1e-3 ? a[0] : a[0] + (a[1] || 0) * (1 / T - 1 / Tr) + (a[2] || 0) * Math.log(T / Tr) + (a[3] || 0) * (T - Tr) + (a[4] || 0) * (T * T - Tr * Tr) + (a[5] || 0) * (1 / (T * T) - 1 / (Tr * Tr)); };
// Chebyshev coefficients of Pitzer's J integrals for the higher-order electrostatic terms (Harvie–Weare, as coded in PHREEQC pitzer.cpp)
const AKX = [1.925154014814667, -0.060076477753119, -0.029779077456514, -0.007299499690937, 0.000388260636404, 0.000636874599598, 0.000036583601823, -0.000045036975204, -0.00000453789571, 0.000002937706971, 0.000000396566462, -0.000000202099617, -0.000000025267769, 0.00000001352261, 0.000000001229405, -0.000000000821969, -0.000000000050847, 0.000000000046333, 0.000000000001943, -0.000000000002563, -0.000000000010991,
  0.628023320520852, 0.462762985338493, 0.150044637187895, -0.028796057604906, -0.036552745910311, -0.001668087945272, 0.006519840398744, 0.001130378079086, -0.000887171310131, -0.000242107641309, 0.000087294451594, 0.000034682122751, -0.000004583768938, -0.000003548684306, -0.00000025045388, 0.000000216991779, 0.00000008077957, 0.000000004558555, -0.000000006944757, -0.000000002849257, 0.000000000237816];
function jFun(X) {
  const lo = X <= 1, p = lo ? X ** 0.2 : X ** -0.1, Z = lo ? 4 * p - 2 : (40 * p - 22) / 9, DZ = lo ? (0.8 * p) / 2 : (-4 * p) / 18, o = lo ? 0 : 21, BK = new Array(23).fill(0), DK = new Array(23).fill(0);
  BK[20] = AKX[o + 20]; BK[19] = Z * AKX[o + 20] + AKX[o + 19]; DK[19] = AKX[o + 20];
  for (let i = 18; i >= 0; i--) { BK[i] = Z * BK[i + 1] - BK[i + 2] + AKX[o + i]; DK[i] = BK[i + 1] + Z * DK[i + 1] - DK[i + 2]; }
  return [X / 4 - 1 + 0.5 * (BK[0] - BK[2]), X * 0.25 + DZ * (DK[0] - DK[2])];
}
function eTheta(zj, zk, I, Aphi) {
  const xc = 6 * Aphi * Math.sqrt(I), zz = zj * zk, [j1, p1] = jFun(xc * zz), [j2, p2] = jFun(xc * zj * zj), [j3, p3] = jFun(xc * zk * zk), et = (zz * (j1 - j2 / 2 - j3 / 2)) / (4 * I);
  return [et, (zz * (p1 - p2 / 2 - p3 / 2)) / (8 * I * I) - et / I];
}
const gPz = (y) => (y ? (2 * (1 - (1 + y) * Math.exp(-y))) / (y * y) : 0), gpPz = (y) => (y ? (-2 * (1 - (1 + y + (y * y) / 2) * Math.exp(-y))) / (y * y) : 0);
/**
 * Pitzer ion-interaction model (Harvie–Møller–Weare form with the PHREEQC pitzer.dat parameter set, MacInnes scaling).
 * m: molalities in the species order of the model; TK; W = waterDH state; lg: output array of ln γ.
 * Returns { aw (water activity), phi (osmotic coefficient), I }.
 */
export function pitzerGamma(C, m, TK, W, lg) {
  const n = C.n, z = C.z, patm = W.pb / 1.01325; let I = 0, bigZ = 0, osum = 0;
  for (let i = 0; i < n; i++) { lg[i] = 0; I += 0.5 * m[i] * z[i] * z[i]; bigZ += m[i] * Math.abs(z[i]); osum += m[i]; }
  if (!(I > 1e-14)) return { aw: 1, phi: 1, I: 0 };
  const DI = Math.sqrt(I), A0 = W.Aphi, B = 1.2, fOf = (b) => -A0 * (DI / (1 + b * DI) + (2 * Math.log(1 + b * DI)) / b);
  let F = fOf(B), F1 = F, F2 = F;
  if (patm > 1) { const pa1 = (7e-5 + 1.93e-9 * (TK - 250) ** 2) * patm, pa2 = TK > 263 ? 9.65e-10 * (TK - 263) ** 2.773 * patm ** 0.623 : pa1; F1 = fOf(B - Math.min(pa1, 0.2)); F2 = fOf(B - Math.min(pa2, 0.2)); }
  const x2 = 2 * DI, xxx = (1 - (1 + x2 - x2 * x2 * 0.5) * Math.exp(-x2)) / (x2 * x2);
  let gamclm = F1, csum = 0, osmot = (-A0 * I ** 1.5) / (1 + B * DI), Fv = 0;
  if (C.mac[0]) gamclm += I * 2 * pzT(C.mac[0].a, TK); if (C.mac[1]) gamclm += I * 2 * pzT(C.mac[1].a, TK) * xxx; if (C.mac[2]) gamclm += 1.5 * pzT(C.mac[2].a, TK) * I * I;
  for (const q of C.B0) { const mi = m[q.i], mj = m[q.j]; if (!(mi > 0) && !(mj > 0)) continue; const p = pzT(q.a, TK); lg[q.i] += mj * 2 * p; lg[q.j] += mi * 2 * p; osmot += mi * mj * p; }
  for (const L of [C.B1, C.B2]) for (const q of L) { const mi = m[q.i], mj = m[q.j]; if (!(mi > 0) && !(mj > 0)) continue; const p = pzT(q.a, TK); if (!p) continue; const y = q.alpha * DI, g = gPz(y); Fv += (mi * mj * p * gpPz(y)) / I; lg[q.i] += mj * 2 * p * g; lg[q.j] += mi * 2 * p * g; osmot += mi * mj * p * Math.exp(-y); }
  for (const q of C.C0) { const mi = m[q.i], mj = m[q.j]; if (!(mi > 0) && !(mj > 0)) continue; const p = pzT(q.a, TK), s = 2 * Math.sqrt(Math.abs(z[q.i] * z[q.j])); csum += (mi * mj * p) / s; lg[q.i] += (mj * bigZ * p) / s; lg[q.j] += (mi * bigZ * p) / s; osmot += (mi * mj * bigZ * p) / s; }
  for (const q of C.TH) { const mi = m[q.i], mj = m[q.j]; if (!(mi > 0) && !(mj > 0)) continue; const p = pzT(q.a, TK); lg[q.i] += 2 * mj * p; lg[q.j] += 2 * mi * p; osmot += mi * mj * p; }
  for (const [i, j] of C.like) { const mi = m[i], mj = m[j]; if (!(mi > 0) && !(mj > 0)) continue; const [et, etp] = eTheta(z[i], z[j], I, A0); Fv += mi * mj * etp; lg[i] += 2 * mj * et; lg[j] += 2 * mi * et; osmot += mi * mj * (et + I * etp); }
  for (const q of C.PSI) { const mi = m[q.i], mj = m[q.j], mk = m[q.k]; if (!(mk > 0)) continue; const p = pzT(q.a, TK); lg[q.i] += mj * mk * p; lg[q.j] += mi * mk * p; lg[q.k] += mi * mj * p; osmot += mi * mj * mk * p; }
  for (const q of C.LA) { const mi = m[q.i], mj = m[q.j]; if (!(mi > 0) && !(mj > 0)) continue; const p = pzT(q.a, TK), same = q.i === q.j, lc = same ? 1 : 2; lg[q.i] += mj * p * lc; lg[q.j] += mi * p * lc; osmot += mi * mj * p * (same ? 0.5 : 1); }
  F += Fv; F1 += Fv; F2 += Fv;
  for (let i = 0; i < n; i++) { const za = Math.abs(z[i]); if (za) lg[i] += za * za * (za === 1 ? F1 : za === 2 ? F2 : F) + za * csum; }
  const phimac = lg[C.iCl] - gamclm; for (let i = 0; i < n; i++) lg[i] += z[i] * phimac;
  const phi = 1 + (2 * osmot) / osum;
  return { aw: Math.exp((-osum * phi) / 55.50837), phi, I };
}
let co2Fluids = null;
/** Fugacity (bar) of CO₂ at mole fraction y in a methane-rich gas from the kernel Peng–Robinson model. */
export function co2Fugacity(y, Pbar, TK) {
  if (!(y > 0)) return 0;
  co2Fluids ||= new Map(); const key = Math.round(clamp(y, 1e-6, 1) * 1e5); let f = co2Fluids.get(key);
  if (!f) { const yy = key / 1e5; f = makeFluid({ comp: yy >= 0.99999 ? { CO2: 100 } : { C1: 100 * (1 - yy), CO2: 100 * yy } }); f._i = f.comps.findIndex((c) => c.id === 'CO2'); if (co2Fluids.size > 200) co2Fluids.clear(); co2Fluids.set(key, f); }
  return y * Pbar * Math.exp(eosPhase(f, f.z, Pbar, TK, 'vapour').lnphi[f._i]);
}
/**
 * Aqueous speciation at Tc (°C), Pbar (bara).
 * tot: { Na, K, Ca, Mg, Ba, Sr, Fe, Cl, SO4 } total molalities (mol/kg water); o: { alk (eq/kg: carbonate alkalinity), fCO2 (bar;
 *   open system) | cT (mol/kg total inorganic carbon; closed system), model: 'pitzer' | 'truesdellJones' | 'davies' }.
 * With neither fCO2 nor cT the water is carbonate-free and neutral.
 * Returns { m, g (ln γ), la (ln activity) by species index, names, I, aw, pH, fCO2 (bar, equilibrium value), lnK(mineral key) → ln Ksp(T, P, I), si(mineral) }.
 */
export function speciate(tot, Tc, Pbar, { alk = 0, fCO2 = 0, cT = 0, model = 'pitzer' } = {}) {
  const C = chemModel(model), n = C.n, T = Tc + KEL, W = waterDH(Tc, Pbar), nm = C.names, m = new Float64Array(n), lg = new Float64Array(n), la = new Float64Array(n), mc = new Float64Array(C.cx.length), gc = new Float64Array(C.cx.length), lk = new Float64Array(C.cx.length);
  const tt = MASTERS.map((id) => Math.max(num(tot?.[id], 0), 0)), carb = fCO2 > 0 || cT > 0, iH = 9, iOH = 10, iHCO3 = 11, iCO3 = 12, iCO2 = 13, pTerm = (dv) => (-dv * (W.pb / 1.01325 - 1)) / (R_ATM * T);
  for (let i = 0; i < 9; i++) m[i] = tt[i];
  let I = 0, aw = 1, pH = 7, lnKw = 0, lnKb = 0, lnKa = 0, lnKh = 0, vmS = new Array(n).fill(0), fOut = fCO2;
  const dh = (zz, g, Im) => { const s = Math.sqrt(Im); if (!zz) return LN10 * (g ? g[1] : 0.1) * Im; return g && C.model !== 'davies' ? LN10 * ((-W.A * zz * zz * s) / (1 + W.B * g[0] * s) + g[1] * Im) : -LN10 * W.A * zz * zz * (s / (1 + s) - 0.3 * Im); };
  const gammas = () => {
    I = 0; for (let i = 0; i < n; i++) I += 0.5 * m[i] * C.z[i] ** 2; for (let k = 0; k < mc.length; k++) I += 0.5 * mc[k] * C.cx[k].z ** 2;
    if (C.pz) { aw = pitzerGamma(C, m, T, W, lg).aw; for (let k = 0; k < mc.length; k++) gc[k] = 0; }
    else { let sm = 0; for (let i = 0; i < n; i++) { lg[i] = dh(C.z[i], C.gam[i], I); sm += m[i]; } for (let k = 0; k < mc.length; k++) { gc[k] = dh(C.cx[k].z, C.cx[k].gamma, I); sm += mc[k]; } aw = 1 - 0.017 * sm; }
    for (let i = 0; i < n; i++) vmS[i] = speciesVm(C.vm[i], C.z[i], W, I);
    const vW = 18.016 / W.rho;
    lnKw = LN10 * logKT(C.kW, T) + pTerm(vmS[iOH] - vW); lnKb = LN10 * logKT(C.kB, T) + pTerm(vmS[iHCO3] - vmS[iCO3]); lnKa = LN10 * logKT(C.kA, T) + pTerm(vmS[iCO2] + vW - vmS[iCO3]); lnKh = LN10 * logKT(C.kH, T) + pTerm(vmS[iCO2]);
    for (let k = 0; k < mc.length; k++) { const c = C.cx[k]; let dv = speciesVm(c.e.vm, c.z, W, I); for (let q = 0; q < c.idx.length; q++) dv -= c.coef[q] * vmS[c.idx[q]]; lk[k] = LN10 * logKT(c.e, T) + pTerm(dv); }
  };
  const bound = new Float64Array(9);
  const inner = (lnH) => { // free molalities at a given proton activity; returns the alkalinity
    const lnAw = Math.log(aw); la[iH] = lnH; m[iH] = Math.exp(lnH - lg[iH]); la[iOH] = lnKw + lnAw - lnH; m[iOH] = Math.exp(la[iOH] - lg[iOH]);
    let lnCO2 = carb && fCO2 > 0 ? lnKh + Math.log(fCO2) : -700;
    for (let it = 0; it < 60; it++) {
      if (carb) { la[iCO2] = lnCO2; la[iCO3] = lnCO2 + lnAw - lnKa - 2 * lnH; la[iHCO3] = lnKb + la[iCO3] + lnH; m[iCO2] = Math.exp(lnCO2 - lg[iCO2]); m[iCO3] = Math.exp(la[iCO3] - lg[iCO3]); m[iHCO3] = Math.exp(la[iHCO3] - lg[iHCO3]); }
      else { la[iCO2] = la[iCO3] = la[iHCO3] = -700; m[iCO2] = m[iCO3] = m[iHCO3] = 0; }
      for (let i = 0; i < 9; i++) la[i] = m[i] > 0 ? Math.log(m[i]) + lg[i] : -700;
      let ch = 0, cCarb = 0; bound.fill(0);
      for (let k = 0; k < mc.length; k++) { const c = C.cx[k]; let s = lk[k]; for (let q = 0; q < c.idx.length; q++) s += c.coef[q] * la[c.idx[q]]; mc[k] = s > -600 ? Math.exp(s - gc[k]) : 0; for (let q = 0; q < c.idx.length; q++) { const j = c.idx[q]; if (j < 9) bound[j] += c.coef[q] * mc[k]; else if (j === iCO3 || j === iHCO3) cCarb += c.coef[q] * mc[k]; } }
      for (let i = 0; i < 9; i++) { if (!(tt[i] > 0)) continue; const x = (tt[i] * m[i]) / (m[i] + bound[i]); ch = Math.max(ch, Math.abs(x - m[i]) / tt[i]); m[i] = x; }
      if (cT > 0 && !(fCO2 > 0)) { const sumC = m[iCO2] + m[iHCO3] + m[iCO3] + cCarb, nl = lnCO2 + Math.log(cT / Math.max(sumC, 1e-300)); ch = Math.max(ch, Math.abs(nl - lnCO2)); lnCO2 = nl; }
      if (ch < 1e-10) break;
    }
    let a = m[iHCO3] + 2 * m[iCO3] + m[iOH] - m[iH]; for (let k = 0; k < mc.length; k++) a += C.cx[k].alk * mc[k];
    if (carb) fOut = Math.exp(lnCO2 - lnKh);
    return a;
  };
  let lnH = -7 * LN10;
  for (let outer = 0; outer < 40; outer++) {
    const I0 = I; gammas();
    if (carb) { const f = (x) => inner(x) - alk; const wd = outer === 0 ? 1.5 : 0.1; let lo = lnH - wd, hi = lnH + wd, flo = f(lo), fhi = f(hi), guard = 0; while (flo < 0 && guard++ < 30) { lo -= 2; flo = f(lo); } guard = 0; while (fhi > 0 && guard++ < 30) { hi += 2; fhi = f(hi); } lnH = flo >= 0 && fhi <= 0 ? brent(f, lo, hi, 1e-9, 80) : flo < 0 ? lo : hi; inner(lnH); } // alkalinity rises with pH, i.e. falls with ln a(H+)
    else { lnH = 0.5 * (lnKw + Math.log(aw) + lg[iH] - lg[iOH]); inner(lnH); }
    if (outer > 0 && Math.abs(I - I0) < 2e-8 * Math.max(I, 1e-6)) { gammas(); inner(lnH); break; }
  }
  pH = -lnH / LN10;
  const lnK = (key) => { const e = C.phases[key], mn = MINERALS.find((q) => q.key === key); let dv = -(e.vm ? e.vm[0] : 0); if (mn) dv += vmS[nm.indexOf(mn.cat)] + vmS[nm.indexOf(mn.an)] + (mn.nW || 0) * (18.016 / W.rho); return LN10 * logKT(e, T) + pTerm(dv); };
  const si = (mn) => { const a = la[nm.indexOf(mn.cat)] + la[nm.indexOf(mn.an)] + (mn.nW || 0) * Math.log(aw); return a > -600 ? (a - lnK(mn.key)) / LN10 : -99; };
  return { m, g: lg, la, names: nm, I, aw, pH, fCO2: fOut, lnK, si, complexes: C.cx.map((c, k) => ({ name: c.name, m: mc[k] })) };
}
/**
 * Solubility (mol/kg water) of a mineral in a background electrolyte.
 * o: { Tc, Pbar, bg: { Na, Cl, … } background totals (mol/kg), pCO2 (bar; calcite/siderite in an open system — omit for a closed system), model }.
 */
export function mineralSolubility(id, { Tc = 25, Pbar = 1.01325, bg = {}, pCO2 = 0, model = 'pitzer' } = {}) {
  const mn = MINERALS.find((q) => q.id === id); need(mn, `Unknown mineral ${id}.`);
  const carbM = mn.an === 'CO3', f = (lx) => { const x = Math.exp(lx), t = { ...bg, [mn.cat]: num(bg[mn.cat], 0) + x }; if (!carbM) t.SO4 = num(bg.SO4, 0) + x; return speciate(t, Tc, Pbar, carbM ? (pCO2 > 0 ? { alk: 2 * x, fCO2: pCO2, model } : { alk: 2 * x, cT: x, model }) : { model }).si(mn); };
  let lo = Math.log(1e-9), hi = Math.log(carbM ? 0.5 : 2); if (f(lo) > 0) return 1e-9; if (f(hi) < 0) return Math.exp(hi);
  return Math.exp(brent(f, lo, hi, 1e-9, 80));
}
/** Ionic strength (mol/L) and molar concentrations of a water analysis given in mg/L: { I, m: { ion: mol/L }, tds (mg/L), balance (charge-balance error, fraction) }. */
export function ionicStrength(water) {
  const m = {}; let I = 0, tds = 0, cat = 0, an = 0;
  for (const id of ION_IDS) { const c = Math.max(num(water?.[id], 0), 0), ion = IONS[id]; m[id] = c / (ion.M * 1000); I += 0.5 * m[id] * ion.z * ion.z; tds += c; if (ion.z > 0) cat += m[id] * ion.z; else an -= m[id] * ion.z; }
  return { I, m, tds, balance: cat + an > 0 ? (cat - an) / (cat + an) : 0 };
}
/**
 * Single-ion activity coefficient. model: 'davies' | 'truesdellJones' (extended Debye–Hückel with the ion-size and
 * salting-out parameters a, b of phreeqc.dat). ion: optional { a (Å), b }. The slopes A and B follow the dielectric constant of water.
 */
export function activityCoefficient(z, I, Tc = 25, model = 'davies', ion = null) {
  const W = waterDH(Tc, 1.01325), s = Math.sqrt(I);
  if (model === 'truesdellJones' && ion) return 10 ** ((-W.A * z * z * s) / (1 + W.B * ion.a * s) + ion.b * I);
  return 10 ** (-W.A * z * z * (s / (1 + s) - 0.3 * I));
}
/** Mix two water analyses (mg/L) with a volume fraction f of the second one. */
export const mixWaters = (a, b, f) => Object.fromEntries(ION_IDS.map((id) => [id, (1 - f) * num(a?.[id], 0) + f * num(b?.[id], 0)]));
/** Molalities (mol/kg water) of a water analysis in mg/L, with the kilograms of water per litre of the analysed sample. */
export function waterMolality(water) {
  const st = ionicStrength(water), kgw = Math.max(waterDensity(25, Math.min(st.tds / 1000, 260)) / 1000 - st.tds / 1e6, 0.5), tot = {};
  for (const id of ION_IDS) tot[id] = st.m[id] / kgw;
  return { tot, kgw, st };
}
/**
 * Saturation indices of the common oilfield scales at T (°C), P (bara) for a water analysis in mg/L (bicarbonate = alkalinity).
 * opt: { yCO2 (mole fraction of CO2 in the gas), model: 'pitzer' (default; PHREEQC pitzer.dat) | 'truesdellJones' (PHREEQC
 *        phreeqc.dat ion-association model with ion pairs) | 'davies' }.
 * The carbonate system follows the CO2 fugacity (Peng–Robinson) and the alkalinity; solubility products carry their
 * temperature dependence (analytic expressions) and pressure dependence (molar volumes of the ions and of the solid).
 * Returns { I (mol/kg), tds, pH, fCO2 (bar), minerals: [{ id, name, SI, ptb (mg/L that can precipitate), logK }], max: { id, name, SI }, aw, gamma (total activity coefficients) }.
 */
export function scaleIndices(water, Tc, Pbar, { yCO2 = 0.03, model = 'pitzer' } = {}) {
  const { tot, kgw, st } = waterMolality(water), T = Tc + KEL, lg = Math.log10, fCO2 = Math.max(co2Fugacity(yCO2, Pbar, T), 1e-9), mdl = model === 'davies' || model === 'truesdellJones' ? model : 'pitzer';
  const sp = speciate(tot, Tc, Pbar, { alk: tot.HCO3, fCO2, model: mdl }), ix = (id) => sp.names.indexOf(id), gT = {};
  for (const id of ['Na', 'K', 'Ca', 'Mg', 'Ba', 'Sr', 'Fe', 'Cl', 'SO4']) gT[id] = tot[id] > 0 ? Math.exp(sp.la[ix(id)]) / tot[id] : Math.exp(sp.g[ix(id)]); // total (stoichiometric) activity coefficients
  gT.HCO3 = tot.HCO3 > 0 ? Math.exp(sp.la[ix('HCO3')]) / tot.HCO3 : 1;
  const lnAw = Math.log(sp.aw), lnKb = sp.la[ix('HCO3')] - sp.la[ix('CO3')] - sp.la[ix('H')], kc = Math.exp(2 * sp.la[ix('HCO3')] - sp.la[ix('CO3')]); // a(HCO3)²/a(CO3) is fixed by the CO2 fugacity
  const minerals = MINERALS.map((mn) => {
    const SI = tot[mn.cat] > 0 ? sp.si(mn) : -99, lnK = sp.lnK(mn.key), K = Math.exp(lnK - (mn.nW || 0) * lnAw), carb = mn.an === 'CO3'; let x = 0;
    if (SI > 0) {
      if (carb) { const fn = (y) => (tot[mn.cat] - y) * (tot.HCO3 - 2 * y) ** 2 * gT[mn.cat] * gT.HCO3 ** 2 - K * kc, hi = Math.min(tot[mn.cat], tot.HCO3 / 2); x = fn(hi) >= 0 ? hi : brent(fn, 0, hi, 1e-14); } // M²⁺ + 2 HCO3⁻ → MCO3 + CO2 + H2O at constant CO2 fugacity
      else { const a = tot[mn.cat], b = tot.SO4; x = 0.5 * (a + b - Math.sqrt((a - b) ** 2 + (4 * K) / (gT[mn.cat] * gT.SO4))); }
    }
    return { id: mn.id, name: mn.name, SI, ptb: Math.max(x, 0) * mn.M * 1000 * kgw, logK: lnK / LN10, M: mn.M, rho: mn.rho };
  });
  const max = minerals.reduce((a, b) => (b.SI > a.SI ? b : a));
  return { I: sp.I, tds: st.tds, pH: sp.pH, fCO2, aw: sp.aw, gamma: gT, kgw, tot, minerals, max: { id: max.id, name: max.name, SI: max.SI }, lnKb };
}

/**
 * Threshold scale inhibitors: class coefficients [a8 (per SI), a9 (K), a10 (per pH unit), a11 (per log R)] and the constant a7 of
 * each product in  log10 b = a7 + a8·SI + a9/T + a10·pH + a11·log10 R,  log10(t_inhibited / t0) = b·C  (C in mg/L of active inhibitor).
 * Fitted by the source to barite induction times at 4–175 °C, pH 3–8, SI 0.56–2.74 (± 0.5 log units).
 */
export const SCALE_INHIBITORS = Object.freeze({
  classes: Object.freeze({ phosphonate: [-1.41, 1329.29, 0.15, 0.18], polycarboxylate: [-1.28, 1007.30, 0.03, 0.13], sulfonated: [-1.72, 1229.65, -0.01, 0.20] }),
  DTPMP: { cls: 'phosphonate', a7: -2.19, name: 'DTPMP (phosphonate)' }, BHPMP: { cls: 'phosphonate', a7: -2.19, name: 'BHPMP (phosphonate)' }, NTMP: { cls: 'phosphonate', a7: -2.47, name: 'NTMP (phosphonate)' }, HEDP: { cls: 'phosphonate', a7: -2.19, name: 'HEDP (phosphonate)' },
  PPCA: { cls: 'polycarboxylate', a7: -0.69, name: 'PPCA (phosphino-polycarboxylate)' }, PASP: { cls: 'polycarboxylate', a7: -0.95, name: 'PASP (polyaspartate)' }, PVS: { cls: 'sulfonated', a7: -0.28, name: 'PVS (polyvinyl sulfonate)' }, SPCA: { cls: 'sulfonated', a7: -0.19, name: 'SPCA (sulfonated polycarboxylate)' },
  // celestite: √(t_inh/t0) − 1 = b·C with b = β1 + β2/T + β3/SI²
  celestite: Object.freeze({ DTPMP: [-13.41, 3526.7, 7.36], PPCA: [-8.24, 1946.5, 6.67], PVS: [-6.82, 1898.0, 3.30] }),
});
/**
 * Induction time (s) of barite or celestite nucleation without inhibitor.
 * Barite: log10 t0 = 1.5232 − 10.8784/SI − 895.6683/T + 5476.992/(SI·T) + 0.8286 [Ca²⁺] + 0.225 log10 R  (4–90 °C; above 90 °C the
 * extended fit −2.11 − 4.29/SI + 279.29/T + 3332.26/(SI·T) + 0.8286 [Ca] + 0.99 log R), R = [Ba]/[SO4], [Ca] in mol/L.
 * Celestite: log10 t0 = −1.324 − 3.301/SI + 2462.5/(SI·T). Other minerals have no induction model here: 0 (immediate).
 * Returns Infinity at or below saturation; capped at 1e15 s.
 */
export function scaleInductionTime(id, { SI, TK, Ca = 0, R: rr = 1 }) {
  if (id !== 'barite' && id !== 'celestite') return SI > 0 ? 0 : Infinity;
  if (!(SI > 0)) return Infinity;
  const s = Math.max(SI, 0.02), lr = Math.log10(clamp(rr, 1e-3, 1e3));
  const lt = id === 'celestite' ? -1.324 - 3.301 / s + 2462.5 / (s * TK) : TK <= 363.15 ? 1.5232 - 10.8784 / s - 895.6683 / TK + 5476.992 / (s * TK) + 0.8286 * Ca + 0.225 * lr : -2.11 - 4.29 / s + 279.29 / TK + 3332.26 / (s * TK) + 0.8286 * Ca + 0.99 * lr;
  return 10 ** clamp(lt, -3, 15);
}
/**
 * Effect of a threshold inhibitor on a supersaturated sulphate scale.
 * o: { SI, TK, pH, R ([cation]/[SO4]), Ca (mol/L), inh (key of SCALE_INHIBITORS), dose (mg/L active), tProtect (s) }
 * Returns { t0, tInh (s), b (L/mg), mic (mg/L for tInh = tProtect; 0 when the brine is safe untreated), eff (1 − t0/tInh: share of
 *           the nucleation and growth sites blocked), risk (tProtect/tInh, capped at 1), model (false: no inhibition model for this mineral) }.
 */
export function scaleInhibition(id, { SI, TK, pH = 6, R: rr = 1, Ca = 0, inh = 'DTPMP', dose = 0, tProtect = 86400 }) {
  const t0 = scaleInductionTime(id, { SI, TK, Ca, R: rr }), P = SCALE_INHIBITORS[inh] || SCALE_INHIBITORS.DTPMP;
  if (!(SI > 0)) return { t0, tInh: Infinity, b: 0, mic: 0, eff: 1, risk: 0, model: true };
  if (id === 'barite') {
    const c = SCALE_INHIBITORS.classes[P.cls], b = 10 ** (P.a7 + c[0] * SI + c[1] / TK + c[2] * pH + c[3] * Math.log10(clamp(rr, 1e-3, 1e3))), tInh = t0 * 10 ** Math.min(b * dose, 30);
    return { t0, tInh, b, mic: t0 >= tProtect ? 0 : Math.log10(tProtect / t0) / b, eff: 1 - t0 / tInh, risk: Math.min(tProtect / tInh, 1), model: true };
  }
  if (id === 'celestite') {
    const q = SCALE_INHIBITORS.celestite[inh] || SCALE_INHIBITORS.celestite[P.cls === 'phosphonate' ? 'DTPMP' : P.cls === 'sulfonated' ? 'PVS' : 'PPCA'], b = Math.max(q[0] + q[1] / TK + q[2] / (SI * SI), 0), tInh = t0 * (1 + b * dose) ** 2;
    return { t0, tInh, b, mic: t0 >= tProtect ? 0 : b > 0 ? (Math.sqrt(tProtect / t0) - 1) / b : Infinity, eff: 1 - t0 / tInh, risk: Math.min(tProtect / tInh, 1), model: true };
  }
  return { t0: 0, tInh: 0, b: 0, mic: Infinity, eff: 0, risk: clamp(SI, 0, 1), model: false }; // carbonates and calcium sulphates: no inhibition credit is taken
}
/** Saturation index at which the induction time equals tProtect (the kinetic "critical" SI), with or without inhibitor. */
export function scaleCriticalSI(id, o) {
  const f = (SI) => Math.log10(Math.max(scaleInhibition(id, { ...o, SI }).tInh, 1e-9)) - Math.log10(o.tProtect);
  if (id !== 'barite' && id !== 'celestite') return 0;
  if (f(0.03) < 0) return 0.03; if (f(6) > 0) return 6;
  return brent(f, 0.03, 6, 1e-6, 80);
}
// =====================================================================================================
// 9. Asphaltene
// =====================================================================================================
/** de Boer plot: boundary points (bar of undersaturation) against in-situ oil density (g/cm³), read from a reproduction of the plot. */
const DE_BOER = Object.freeze({ rho: [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.775, 0.8, 0.85], lower: [64.8, 84.1, 93.1, 114.5, 156.5, 245.5, 311.6, 395.1, 620.5], upper: [120.0, 163.4, 181.3, 214.4, 299.2, 467.5, 597.1] });
/**
 * de Boer screening: undersaturation (reservoir pressure − saturation pressure, bar) against the in-situ oil density (kg/m³).
 * The two boundaries are interpolated between points read from a reproduction of the published plot (reading error about ± 10 bar and
 * ± 5 kg/m³; the severe-problem line is extrapolated above 775 kg/m³); for screening only.
 * Returns { cls: 'no problems' | 'slight problems' | 'severe problems', lower, upper (bar at this density) }.
 */
export function deBoer(rho, dP) {
  const x = rho / 1000, lower = interp1(DE_BOER.rho, DE_BOER.lower, clamp(x, 0.5, 0.85)), upper = x <= 0.775 ? interp1(DE_BOER.rho.slice(0, 7), DE_BOER.upper, Math.max(x, 0.5)) : DE_BOER.upper[6] + ((DE_BOER.upper[6] - DE_BOER.upper[5]) / 0.025) * (x - 0.775);
  return { cls: dP <= lower ? 'no problems' : dP <= upper ? 'slight problems' : 'severe problems', lower, upper };
}
/** Colloidal instability index from a SARA analysis (wt %): (saturates + asphaltenes) / (aromatics + resins). */
export function colloidalInstability({ sat, aro, res, asp }) {
  const cii = (sat + asp) / Math.max(aro + res, 1e-9);
  return { cii, cls: cii < 0.7 ? 'stable' : cii <= 0.9 ? 'uncertain' : 'unstable' };
}
/**
 * Flory–Huggins (Hirschberg) maximum volume fraction of asphaltene soluble in a live oil.
 * { rhoL (kg/m³), mwL (g/mol), TK, deltaA (MPa^0.5 at 25 °C), vA (m³/kmol) }. The oil solubility parameter follows
 * δ = 17.347 ρ(g/cm³) + 2.904 MPa^0.5. Returns { phiMax, deltaL, deltaA }.
 */
export function asphalteneSolubility({ rhoL, mwL, TK, deltaA = 20.0, vA = 2.0 }) {
  const deltaL = 17.347 * (rhoL / 1000) + 2.904, dA = deltaA * (1 - 1.07e-3 * (TK - 298.15)), vL = mwL / rhoL; // m³/kmol
  const ex = (vA / vL) * (1 - vL / vA) - (vA * 1e-3 * ((dA - deltaL) * 1e3) ** 2) / (R * TK); // δ in Pa^0.5, v in m³/mol
  return { phiMax: Math.min(Math.exp(Math.min(ex, 0)), 1), deltaL, deltaA: dA };
}

// =====================================================================================================
// 10. Sand
// =====================================================================================================
/**
 * Threshold Shields number of a grain bed (Soulsby–Whitehouse): θcr = 0.30/(1 + 1.2 D*) + 0.055 [1 − exp(−0.020 D*)],
 * D* = d [g (s − 1)/ν²]^⅓. It rises from 0.055 for coarse grains in water to 0.30 in the viscous limit.
 */
export function shieldsCritical(d, rhoP, rhoF, mu) {
  const nu = mu / rhoF, Ds = d * Math.cbrt((G * Math.max(rhoP / rhoF - 1, 1e-9)) / (nu * nu));
  return 0.30 / (1 + 1.2 * Ds) + 0.055 * (1 - Math.exp(-0.020 * Ds));
}
/** Shields number of a fully mobile grain layer in laminar flow: upper end of the range of the viscous bed-load theory (Ouriemi et al.), where the moving layer is θ/(2θc) ≈ 6 grains thick. */
export const SHIELDS_MOBILE = 1.5;
/**
 * Minimum transport (critical) velocity of sand in a mostly horizontal pipe, selected by flow regime.
 * { d (m), D (m), rhoP, rhoF (carrier liquid), mu (Pa·s), C (sand volume fraction), vsl, vm (m/s, for Salama's liquid fraction), thetaM }
 *  · turbulent carrier, grain larger than the viscous sub-layer (d⁺ ≥ 5): the largest of Oroskar–Turian, Salama and Danielson
 *    (turbulent suspension);
 *  · turbulent carrier, grain inside the viscous sub-layer (d⁺ < 5): the Thomas limit of a bed sliding in the sub-layer,
 *    V = 9.0 [g ν (s − 1)]^0.37 (D/ν)^0.11;
 *  · laminar carrier (Re < 2300 at that velocity): turbulence cannot suspend the grains; the deposit is cleared as a fully
 *    mobile bed-load layer when the wall Shields number 8 μ V/(D Δρ g d) reaches thetaM.
 * Returns { oroskarTurian, salama, danielson, thomas, laminar (m/s), governing, regime, dPlus, Re, shieldsCr, settling }.
 */
export function sandCriticalVelocity({ d, D, rhoP = 2650, rhoF = 800, mu = 2e-3, C = 1e-4, vsl = 1, vm = 2, coef = 1, shape = 1, thetaM = SHIELDS_MOBILE }) {
  const s = rhoP / rhoF, nu = mu / rhoF, c = clamp(C, 1e-7, 0.5), sm = Math.max(s - 1, 1e-6), root = Math.sqrt(G * d * sm);
  const oroskarTurian = root * 1.85 * c ** 0.1536 * (1 - c) ** 0.3564 * (d / D) ** -0.378 * ((D * rhoF * root) / mu) ** 0.09 * 0.95 ** 0.3;
  const salama = clamp(vsl / Math.max(vm, 1e-9), 0.01, 1) ** 0.53 * d ** 0.17 * nu ** -0.09 * sm ** 0.55 * D ** 0.47;
  const danielson = 0.23 * nu ** (-1 / 9) * d ** (1 / 9) * (G * D * sm) ** (5 / 9);
  const thomas = 9.0 * (G * nu * sm) ** 0.37 * (D / nu) ** 0.11, laminar = (thetaM * sm * rhoF * G * d * D) / (8 * mu);
  const turb = Math.max(oroskarTurian, salama, danielson), dPl = (V) => { const Re = (V * D) / nu; return (d * V * Math.sqrt(frictionFactor(Math.max(Re, 2300), 0) / 8)) / nu; };
  let governing = turb, regime = 'turbulent suspension', dPlus = dPl(turb);
  if (dPlus < 5) { governing = thomas; regime = 'viscous sub-layer (sliding bed)'; dPlus = dPl(thomas); } // grains hidden in the sub-layer are not reached by the eddies the suspension correlations rely on
  if ((governing * D) / nu < 2300) { // the carrier would be laminar at the turbulent requirement: bed-load limit, capped where the flow turns turbulent
    const vT = (2300 * nu) / D; if (laminar <= vT) { governing = laminar; regime = 'laminar (mobile bed load)'; dPlus = (d * Math.sqrt((8 * nu * laminar) / D)) / nu; } else { governing = Math.max(vT, Math.min(governing, laminar)); regime = 'transitional'; }
  }
  return { oroskarTurian: coef * oroskarTurian, salama: coef * salama, danielson: coef * danielson, thomas: coef * thomas, laminar: coef * laminar, governing: coef * governing, regime, dPlus, Re: (governing * D) / nu, shieldsCr: shieldsCritical(d, rhoP, rhoF, mu), settling: settlingVelocity(d, rhoP, rhoF, mu, { shape }).v };
}
/** Screening erosion rate (mm/y) in a bend from sand: Salama (2000), E = W·V²·d / (Sm·D²·ρm) with W kg/d, d µm, D mm, Sm = 5.5. */
export const sandErosionScreen = (Wkgd, vm, dUm, Dmm, rhoM, Sm = 5.5) => (Wkgd * vm * vm * dUm) / (Sm * Dmm * Dmm * Math.max(rhoM, 1));

// =====================================================================================================
// 10a. Stochastic and field solvers: Monte Carlo population balance, Maxey–Riley particle dynamics,
//      Eulerian–Lagrangian tracking, Eulerian–Eulerian solids transport, reaction–diffusion
// =====================================================================================================
/** Error function (Abramowitz–Stegun 7.1.26, absolute error below 1.5×10⁻⁷) and its complement. */
export const erf = (x) => { const s = Math.sign(x), a = Math.abs(x), t = 1 / (1 + 0.3275911 * a); return s * (1 - ((((1.061405429 * t - 1.453152027) * t + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * Math.exp(-a * a)); };
export const erfc = (x) => 1 - erf(x);
/**
 * Monte Carlo (direct simulation, majorant-kernel acceptance–rejection) solution of the aggregation–breakage population
 * balance in a well-mixed volume. L0: sizes (m) of the simulated particles; conc: number concentration they represent (1/m³).
 * o: { beta: number (constant kernel, m³/s) | (Li, Lj) → m³/s, S: (L) → 1/s (binary equal-volume breakage; must not decrease with size), seed }.
 * The sample is doubled when it has halved (and thinned with an exact volume correction when it has doubled), so the
 * statistical resolution stays constant. Returns { L[], n, m0 (1/m³), vol (m³/m³), d43, d10, events, accepted }.
 */
export function solvePBEMonteCarlo(L0, t, { beta = 0, S = null, conc = 1, seed = 1, maxEvents = 4e6 } = {}) {
  const n0 = L0.length, cap = 2 * n0, v = new Float64Array(cap + 2), rg = rng(seed), cst = typeof beta === 'number', lOf = (x) => Math.cbrt((6 * x) / PI);
  let N = n0, Vs = n0 / conc, time = 0, vmax = 0, vmin = Infinity, events = 0, accepted = 0;
  for (let i = 0; i < n0; i++) { v[i] = (PI / 6) * L0[i] ** 3; if (v[i] > vmax) vmax = v[i]; if (v[i] < vmin) vmin = v[i]; }
  const maj = () => { if (cst) return beta; const a = lOf(vmax), b = lOf(vmin); return beta(a, a) + beta(a, b) + beta(b, b); };
  let bm = maj(), sm = S ? S(lOf(vmax)) : 0;
  while (events < maxEvents) {
    const Ra = (bm * N * (N - 1)) / (2 * Vs), Rb = sm * N, Rt = Ra + Rb; if (!(Rt > 0)) break;
    const tau = -Math.log(1 - rg.uniform()) / Rt; if (time + tau > t) break; time += tau; events++;
    if (rg.uniform() * Rt < Ra) {
      const i = rg.int(N); let j = rg.int(N - 1); if (j >= i) j++;
      if (cst || rg.uniform() * bm < beta(lOf(v[i]), lOf(v[j]))) { accepted++; v[i] += v[j]; v[j] = v[N - 1]; N--; if (v[i] > vmax) { vmax = v[i]; bm = maj(); if (S) sm = S(lOf(vmax)); } }
    } else {
      const i = rg.int(N);
      if (rg.uniform() * sm < S(lOf(v[i]))) { accepted++; v[i] /= 2; v[N++] = v[i]; if (v[i] < vmin) { vmin = v[i]; bm = maj(); } }
    }
    if (N <= n0 / 2) { for (let i = 0; i < N; i++) v[N + i] = v[i]; N *= 2; Vs *= 2; }
    else if (N >= cap) { let tot = 0, kept = 0; for (let i = 0; i < N; i++) tot += v[i]; let k = 0; for (let i = 0; i < N; i++) if (rg.uniform() < 0.5) { v[k++] = v[i]; kept += v[i]; } N = k; Vs *= kept / tot; }
  }
  const L = new Array(N); let m3 = 0, m4 = 0, m1 = 0, vol = 0; for (let i = 0; i < N; i++) { const l = lOf(v[i]); L[i] = l; m1 += l; m3 += l ** 3; m4 += l ** 4; vol += v[i]; }
  return { L, n: N, m0: N / Vs, vol: vol / Vs, d43: m3 > 0 ? m4 / m3 : 0, d10: N ? m1 / N : 0, events, accepted };
}
/**
 * One-dimensional Maxey–Riley stepper for a sphere: Stokes/Schiller–Naumann drag, weight minus buoyancy, added mass (½),
 * fluid-acceleration (pressure-gradient) force and the Basset history force with a windowed kernel.
 *   (m_p + ½ m_f) dv/dt = (m_p − m_f) g + (3/2) m_f Du/Dt + 3πμ d φ(Re)(u − v) + (3/2) d² √(π ρ_f μ) ∫ (d(u − v)/dτ) / √(t − τ) dτ
 * The drag and the newest history increment are implicit, so the step is stable for dt far above the particle response time.
 * p: { d, rhoP, rhoF, mu, g (m/s² along the axis, positive = direction of positive velocity), drag: 'schiller' | 'stokes', basset (bool), window (steps kept in the history) }.
 * Returns { step(dt, u1, dudt) → v, v (getter), reset(v, u) }.
 */
export function maxeyRileyStepper({ d, rhoP, rhoF, mu, g = G, drag = 'schiller', basset = true, window = 64 }) {
  const vol = (PI / 6) * d ** 3, mp = rhoP * vol, mf = rhoF * vol, Bc = basset ? 1.5 * d * d * Math.sqrt(PI * rhoF * mu) : 0, hist = new Float64Array(3 * window); // per kept step: d(u − v)/dτ, start time, end time
  let v = 0, u0 = 0, t = 0, nH = 0, head = 0;
  return {
    get v() { return v; }, get t() { return t; },
    reset(v0 = 0, uf = 0) { v = v0; u0 = uf; t = 0; nH = 0; head = 0; },
    step(dt, u1, dudt = (u1 - u0) / dt) {
      const w0 = u0 - v, Re = (rhoF * Math.abs(w0) * d) / mu, kd = 3 * PI * mu * d * (drag === 'stokes' || Re <= 0 ? 1 : Re < 1000 ? 1 + 0.15 * Re ** 0.687 : (0.44 * Re) / 24), t1 = t + dt;
      let H = 0; // history of the earlier steps inside the window: piecewise-constant d(u − v)/dτ integrated exactly against 1/√(t − τ)
      if (Bc) for (let k = 0; k < nH; k++) { const q = 3 * k; H += hist[q] * 2 * (Math.sqrt(t1 - hist[q + 1]) - Math.sqrt(t1 - hist[q + 2])); }
      const cB = Bc ? (Bc * 2 * Math.sqrt(dt)) / dt : 0, M = mp + 0.5 * mf;
      // M (v1 − v)/dt = (mp − mf) g + 1.5 mf Du/Dt + kd (u1 − v1) + Bc H + cB ((u1 − v1) − (u0 − v))
      const v1 = ((M * v) / dt + (mp - mf) * g + 1.5 * mf * dudt + kd * u1 + Bc * H + cB * (u1 - w0)) / (M / dt + kd + cB);
      if (Bc) { const q = 3 * head; hist[q] = (u1 - v1 - w0) / dt; hist[q + 1] = t; hist[q + 2] = t1; head = (head + 1) % window; if (nH < window) nH++; }
      v = v1; u0 = u1; t = t1; return v;
    },
  };
}
/**
 * Maxey–Riley response of a sphere released from rest in a fluid whose velocity is uf(t) (default: still fluid, settling).
 * Returns { t[], v[], vTerminal (Schiller–Naumann or Stokes) }.
 */
export function maxeyRiley({ d, rhoP, rhoF, mu, g = G, tEnd, n = 400, basset = true, drag = 'stokes', uf = null, v0 = 0 }) {
  const st = maxeyRileyStepper({ d, rhoP, rhoF, mu, g, drag, basset, window: Math.min(n, 4096) }), dt = tEnd / n, t = [0], v = [v0]; st.reset(v0, uf ? uf(0) : 0);
  for (let k = 1; k <= n; k++) { st.step(dt, uf ? uf(k * dt) : 0); t.push(k * dt); v.push(st.v); }
  return { t, v, vTerminal: settlingVelocity(d, rhoP, rhoF, mu, { model: drag === 'stokes' ? 'stokes' : 'schiller' }).v * Math.sign(g || 1) };
}
/**
 * Eulerian–Lagrangian tracking of solid parcels through the one-dimensional flow field of the line.
 * Each parcel carries an axial position and a wall-normal position y on a pipe diameter. Axially it moves with the local
 * mean velocity profile (1/7 power law, normalised so that the diameter average is the liquid velocity) minus its axial
 * settling slip. Across the pipe it obeys the Maxey–Riley equation (drag, weight and buoyancy, added mass, fluid
 * acceleration and windowed Basset history) in a fluctuating fluid velocity from a discrete random walk whose intensity
 * follows the eddy diffusivity of pipe flow D_t(y) (Reichardt profile), including the drift dD_t/dy that keeps a passive
 * tracer uniformly mixed. The step is `stepFrac` of the cross-pipe mixing time in the core and shorter near the walls.
 * Walls: a parcel that reaches a wall is reflected with the restitution coefficient unless the shear velocity is below
 * its settling velocity (Rouse criterion), in which case it joins the bed; while it is inside the wall layer of a cold
 * wall (below the hydrate temperature) it is captured at the rate adhesion × V_d / δ, V_d being the turbulent deposition
 * velocity. Equipment (choke, low point, strainer …) traps a passing parcel with its efficiency; a parcel reaching the outlet escapes.
 * field: { ds, D[], vL[], uStar[], rhoF[], mu[], theta[], cold[] }; o: { n, sizes[] (m, sampled per parcel) | d, rhoP, adhesion, restitution, equip: [{ x (m), eff, name }], seed, x0 (m), basset, stepFrac, maxSteps }.
 * Returns { n, wall, bed, equipment, escaped, flying, depX: [] (m, where parcels were captured), hist: [] (captured per cell), eq: [{ name, x, n }], meanTransit (s), meanTravel (m), steps, reflections }.
 */
export function trackParticles(field, { n = 60, sizes = null, d = 1e-4, rhoP = 920, adhesion = 0.01, restitution = 0.8, equip = [], seed = 3, x0 = 0, basset = true, maxSteps = 400000, window = 6, stepFrac = 0.05, TK = 277 } = {}) {
  const nc = field.D.length, ds = field.ds, Ltot = nc * ds, rg = rng(seed), hist = new Array(nc).fill(0), depX = [], eq = equip.map((e) => ({ name: e.name || 'equipment', x: e.x, eff: clamp(num(e.eff, 0), 0, 1), n: 0 })).sort((a, b) => a.x - b.x);
  let wall = 0, bed = 0, equipment = 0, escaped = 0, flying = 0, steps = 0, reflections = 0, sumT = 0, sumX = 0;
  for (let k = 0; k < n; k++) {
    const dp = sizes ? sizes[rg.int(sizes.length)] : d; let i = clamp(Math.floor(x0 / ds), 0, nc - 1), s = x0, y = rg.uniform(0.02, 0.98) * field.D[i], t = 0, fate = 'flying', nEq = eq.findIndex((e) => e.x > s);
    let mr = null, iMr = -1, vf = 0, wsl = 0, ct = 1, st = 0, us = 0, D = 0, R2 = 0, dtCore = 0, kCap = 0, dWall = 0;
    for (let it = 0; it < maxSteps; it++) {
      if (i !== iMr) { // new cell: local fluid, inclination, settling and deposition velocity (y is measured upwards from the bottom, so gravity is negative)
        const vy = mr ? mr.v : 0; ct = Math.cos(field.theta[i]); st = Math.sin(field.theta[i]); D = field.D[i]; R2 = D / 2; us = Math.max(field.uStar[i], 1e-5); dWall = 0.1 * R2;
        mr = maxeyRileyStepper({ d: dp, rhoP, rhoF: field.rhoF[i], mu: field.mu[i], g: -G * ct, basset, window }); mr.reset(vy, vf); iMr = i;
        wsl = settlingVelocity(dp, rhoP, field.rhoF[i], field.mu[i]).v; dtCore = Math.min((stepFrac * D) / (0.07 * us), (0.25 * ds) / Math.max(Math.abs(field.vL[i]), 1e-3));
        kCap = field.cold[i] && adhesion > 0 ? (adhesion * depositionVelocity(dp, rhoP, us, field.mu[i] / field.rhoF[i], field.rhoF[i], TK)) / dWall : 0;
      }
      const yw = Math.max(Math.min(y, D - y), 0), q = Math.min(Math.abs(y - R2) / R2, 1), q2 = q * q, cD = (0.4 * us * R2) / 6, Dt = cD * (1 - q2) * (1 + 2 * q2), drift = ((cD * (2 * q - 8 * q2 * q)) / R2) * (y < R2 ? 1 : -1); // Reichardt eddy diffusivity and its gradient in y
      const ywc = Math.max(yw, dWall), dt = Math.min(dtCore, (0.5 * ywc * ywc) / Math.max(Dt, (0.4 * us * ywc) / 6)), ua = (8 / 7) * field.vL[i] * (1 - Math.min(q, 0.999)) ** (1 / 7) - wsl * st;
      const vNew = drift + Math.sqrt((2 * Dt) / dt) * rg.normal(), vy = mr.step(dt, vNew, (vNew - vf) / dt); vf = vNew;
      if (kCap > 0 && yw < dWall && rg.uniform() < 1 - Math.exp(-kCap * dt)) { fate = 'wall'; break; }
      y += vy * dt; const sOld = s; s += ua * dt; t += dt; steps++;
      if (nEq >= 0 && nEq < eq.length && sOld < eq[nEq].x && s >= eq[nEq].x) { if (rg.uniform() < eq[nEq].eff) { fate = 'equipment'; eq[nEq].n++; s = eq[nEq].x; break; } nEq++; }
      if (s >= Ltot) { fate = 'escaped'; break; } if (s < 0) s = 0;
      if (y <= dp / 2 || y >= D - dp / 2) {
        const bottom = y <= dp / 2;
        if ((bottom ? wsl : -wsl) * ct > us) { fate = 'bed'; break; } // the turbulence cannot lift the grain again (u* < w_s); a buoyant particle is held at the top in the same way
        y = bottom ? dp - y : 2 * (D - dp / 2) - y; if (!(y > dp / 2 && y < D - dp / 2)) y = bottom ? dp : D - dp; mr.reset(-restitution * mr.v, vf); reflections++;
      }
      i = clamp(Math.floor(s / ds), 0, nc - 1);
    }
    if (fate === 'wall') wall++; else if (fate === 'bed') bed++; else if (fate === 'equipment') equipment++; else if (fate === 'escaped') escaped++; else flying++;
    if (fate !== 'escaped' && fate !== 'flying') { depX.push(s); hist[clamp(Math.floor(s / ds), 0, nc - 1)]++; }
    if (fate === 'escaped') sumT += t; sumX += Math.min(s, Ltot) - x0;
  }
  return { n, wall, bed, equipment, escaped, flying, depX, hist, eq: eq.map((e) => ({ name: e.name, x: e.x, n: e.n })), meanTransit: escaped ? sumT / escaped : null, meanTravel: n ? sumX / n : 0, steps, reflections };
}
/**
 * Eulerian–Eulerian transport of a dispersed solid phase on the line: the solids concentration c (volume fraction of the
 * carrier liquid) is advected with its own phase velocity (liquid velocity minus the axial slip of the settling solids),
 * spreads by axial dispersion, and exchanges mass with a stationary bed (deposition by hindered settling below the
 * deposition shear, entrainment above the critical Shields stress):
 *   ∂(A_L c)/∂t + ∂(A_L u_s c)/∂x = ∂/∂x (A_L D_ax ∂c/∂x) − W (w_h c p_d − E),   ∂V_b/∂t = W (w_h c p_d − E)
 * Advection is explicit (flux-limited, van Leer) in sub-steps below the Courant limit, dispersion and bed exchange are implicit.
 * f: { ds, A[] (liquid flow area, m²), u[] (solids phase velocity), W[] (bed chord width), Dax[], ws[] (hindered settling normal to the wall),
 *      pd[] (deposition probability 0–1), E[] (entrainment capacity m/s as solids volume flux when a bed exists) };
 * o: { c0[] | 0, bed0[] | 0 (m³ solids per m), cIn (inlet concentration, number or (t) → number), tEnd, dt, trap: [{ i, eff }] (equipment that removes a fraction of the passing solids), limiter (bool) }.
 * Returns { c[], bed[] (m³/m), t, ledger: { in, out, suspended, bed, trapped, initial }, steps }.
 */
export function solidsTransport(f, { c0 = 0, bed0 = 0, cIn = 0, tEnd, dt, trap = [], limiter = true } = {}) {
  const n = f.A.length, ds = f.ds, c = Array.from({ length: n }, (_, i) => (Array.isArray(c0) ? c0[i] : c0)), bed = Array.from({ length: n }, (_, i) => (Array.isArray(bed0) ? bed0[i] : bed0));
  const cin = typeof cIn === 'function' ? cIn : () => cIn, led = { in: 0, out: 0, trapped: 0, initial: 0 }, trapEff = new Array(n).fill(0), trapped = new Array(n).fill(0); for (const tr of trap) if (tr.i >= 0 && tr.i < n) trapEff[tr.i] = clamp(tr.eff, 0, 1);
  for (let i = 0; i < n; i++) led.initial += c[i] * f.A[i] * ds + bed[i] * ds;
  let cfl = 0; for (let i = 0; i < n; i++) cfl = Math.max(cfl, Math.abs(f.u[i]) / ds);
  const nSteps = Math.max(1, Math.ceil(tEnd / dt)), h = tEnd / nSteps, nSub = Math.max(1, Math.ceil(h * cfl / 0.8)), hs = h / nSub, F = new Float64Array(n + 1), Fo = new Float64Array(n), a = new Array(n), b = new Array(n), cc = new Array(n), r = new Array(n);
  const phi = (th) => (limiter ? (th + Math.abs(th)) / (1 + Math.abs(th)) : 0); let t = 0;
  for (let st = 0; st < nSteps; st++) {
    for (let sb = 0; sb < nSub; sb++) { // explicit conservative advection of A·c (second order in smooth regions)
      F[0] = f.A[0] * f.u[0] * cin(t); led.in += F[0] * hs;
      for (let i = 0; i < n; i++) { const q = f.A[i] * f.u[i], nu = (Math.abs(f.u[i]) * hs) / ds; let face = c[i]; if (i > 0 && i < n - 1) { const dC = c[i + 1] - c[i], th = Math.abs(dC) > 1e-300 ? (c[i] - c[i - 1]) / dC : 0; face += 0.5 * phi(th) * (1 - nu) * dC; } Fo[i] = q * Math.max(face, 0); const lost = Fo[i] * trapEff[i]; F[i + 1] = Fo[i] - lost; trapped[i] += lost * hs; led.trapped += lost * hs; }
      for (let i = 0; i < n; i++) c[i] += (hs * (F[i] - Fo[i])) / (f.A[i] * ds);
      led.out += F[n] * hs; t += hs;
    }
    // implicit dispersion and deposition: (A/h + W ws pd) c − ∂/∂x(A D ∂c/∂x) = A c*/h + W E_eff
    for (let i = 0; i < n; i++) {
      const dl = i > 0 ? (0.5 * (f.A[i] * f.Dax[i] + f.A[i - 1] * f.Dax[i - 1])) / (ds * ds) : 0, dr = i < n - 1 ? (0.5 * (f.A[i] * f.Dax[i] + f.A[i + 1] * f.Dax[i + 1])) / (ds * ds) : 0;
      const ent = Math.min(f.W[i] * f.E[i], bed[i] / h); // the bed cannot give more than it holds
      a[i] = -dl; cc[i] = -dr; b[i] = f.A[i] / h + dl + dr + f.W[i] * f.ws[i] * f.pd[i]; r[i] = (f.A[i] * c[i]) / h + ent; F[i] = ent;
    }
    const x = tridiag(a, b, cc, r);
    for (let i = 0; i < n; i++) { c[i] = Math.max(x[i], 0); bed[i] = Math.max(bed[i] + h * (f.W[i] * f.ws[i] * f.pd[i] * c[i] - F[i]), 0); }
  }
  let susp = 0, bd = 0; for (let i = 0; i < n; i++) { susp += c[i] * f.A[i] * ds; bd += bed[i] * ds; }
  return { c, bed, trapped, t, steps: nSteps * nSub, ledger: { ...led, suspended: susp, bed: bd } };
}
/** Ogata–Banks solution of advection–dispersion of a step of concentration c0 entering a semi-infinite line at x = 0: c(x, t)/c0. */
export const ogataBanks = (x, t, u, D) => 0.5 * (erfc((x - u * t) / (2 * Math.sqrt(D * t))) + Math.exp(Math.min((u * x) / D, 700)) * erfc((x + u * t) / (2 * Math.sqrt(D * t))));
/**
 * Transient reaction–diffusion in a slab: ∂c/∂t = D ∂²c/∂x² − k c with c = c0 at x = 0 and no flux at x = L (implicit Euler,
 * second-order in space). Returns { x[], c[], flux (D ∂c/∂x into the slab at x = 0), thiele (L√(k/D)), effectiveness (tanh φ/φ, analytic) }.
 */
export function reactionDiffusion({ L, D, k, c0 = 1, n = 60, tEnd = null, steps = 200 }) {
  const dx = L / n, x = Array.from({ length: n + 1 }, (_, i) => i * dx), tt = tEnd ?? (20 * L * L) / D, dt = tt / steps, lam = (D * dt) / (dx * dx); let c = new Array(n + 1).fill(0); c[0] = c0;
  const a = new Array(n).fill(-lam), b = new Array(n).fill(1 + 2 * lam + k * dt), cu = new Array(n).fill(-lam); a[0] = 0; cu[n - 1] = 0; a[n - 1] = -2 * lam; // unknowns c1 … cn; mirror node at the closed end
  for (let s = 0; s < steps; s++) { const r = c.slice(1); r[0] += lam * c0; const y = tridiag(a, b, cu, r); for (let i = 0; i < n; i++) c[i + 1] = y[i]; }
  const th = L * Math.sqrt(k / D);
  return { x, c, flux: (D * (3 * c[0] - 4 * c[1] + c[2])) / (2 * dx), thiele: th, effectiveness: th > 1e-9 ? Math.tanh(th) / th : 1 };
}
/**
 * Mass-transfer-controlled dissociation of a hydrate plug face in contact with a thermodynamic inhibitor.
 * The inhibitor must reach the face faster than the released water dilutes it: the interface stays at the equilibrium
 * concentration wEq (the mass fraction at which hydrate is just stable at the local P and T), the bulk is at wBulk.
 *  - mode 'film': convective film, flux = kFilm ρ (wBulk − wEq): constant recession speed;
 *  - mode 'stagnant': one-dimensional diffusion into a quiescent liquid (solved on a grid): recession ∝ √t;
 *  - mode 'porous': the inhibitor also diffuses into the pores of the plug and is consumed on the grains (reaction–diffusion
 *    with rate k = kGrain·a_s): flux = ρ (wBulk − wEq) √(D_eff k) tanh φ.
 * o: { wBulk, wEq (mass fractions), D (m²/s), rhoL, rhoH, eps (plug porosity), wfH (water mass fraction of hydrate), kFilm (m/s), dGrain (m), tort, Lp (plug length m), tEnd (s), n, steps }
 * Returns { s (m recession of the face after tEnd), rate (m/s at tEnd), flux (kg inhibitor/m²/s), tMelt (s to consume Lp/2 from each side | null), sAnalytic (stagnant: similarity solution), thiele }.
 */
export function inhibitorDissociation({ mode = 'film', wBulk, wEq, D = 1.3e-9, rhoL = 950, rhoH = 920, eps = 0.4, wfH = 0.86, kFilm = 1e-5, dGrain = 1e-4, tort = 2, Lp = 50, tEnd = 86400, n = 80, steps = 400 }) {
  const dw = wBulk - wEq, z = { s: 0, rate: 0, flux: 0, tMelt: null, sAnalytic: 0, thiele: 0 };
  if (!(dw > 0) || !(wEq > 0)) return z;
  const perFlux = (1 - wEq) / (wEq * rhoH * (1 - eps) * wfH); // metres of plug dissolved per kg/m² of inhibitor delivered: the melt water must end at wEq
  if (mode === 'stagnant') {
    const Lx = 6 * Math.sqrt(D * tEnd), dx = Lx / n, dt = tEnd / steps, lam = (D * dt) / (dx * dx), c = new Array(n + 1).fill(dw); c[0] = 0; let s = 0, fl = 0;
    const a = new Array(n - 1).fill(-lam), b = new Array(n - 1).fill(1 + 2 * lam), cu = new Array(n - 1).fill(-lam); a[0] = 0; cu[n - 2] = 0;
    for (let k = 0; k < steps; k++) { const r = c.slice(1, n); r[n - 2] += lam * dw; const y = tridiag(a, b, cu, r); for (let i = 0; i < n - 1; i++) c[i + 1] = y[i]; fl = (rhoL * D * (-3 * c[0] + 4 * c[1] - c[2])) / (2 * dx); s += fl * perFlux * dt; }
    const sA = 2 * rhoL * dw * Math.sqrt((D * tEnd) / PI) * perFlux, kk = 2 * rhoL * dw * Math.sqrt(D / PI) * perFlux;
    return { s, rate: fl * perFlux, flux: fl, tMelt: (Lp / 2 / kk) ** 2, sAnalytic: sA, thiele: 0 };
  }
  let flux = kFilm * rhoL * dw, th = 0;
  if (mode === 'porous') { const De = (D * eps) / tort, k = (kFilm * 6 * (1 - eps)) / (dGrain * eps); th = (Lp / 2) * Math.sqrt(k / De); flux += rhoL * dw * eps * Math.sqrt(De * k) * Math.tanh(th); }
  const rate = flux * perFlux;
  return { s: rate * tEnd, rate, flux, tMelt: Lp / 2 / rate, sAnalytic: rate * tEnd, thiele: th };
}
/** Effective wall roughness (m) of a fouled wall: clean roughness plus a fraction of each deposit layer, capped at 5 % of the bore. */
export const roughnessUpdate = (rough0, dHyd, dWax, dScale, D, k = ROUGH_K) => Math.min(rough0 + k.hyd * dHyd + k.wax * dWax + k.scale * dScale, 0.05 * D);
const ROUGH_K = Object.freeze({ hyd: 0.1, wax: 0.05, scale: 0.3 });

// =====================================================================================================
// 10b. Coupled sub-models: slug-unit mass transfer, wax solid–liquid equilibrium, scale crystallisation, asphaltene onset
// =====================================================================================================
/**
 * Gas-to-liquid mass transfer in hydrodynamic slug flow from the kernel slug unit cell. The slug body carries entrained
 * bubbles (Hinze maximum stable size in the body's turbulence, Sauter size 0.6 of it) and is strongly mixed; the film
 * zone behind it has a stratified interface. Conductances k_L·a (1/s, per unit liquid volume) follow the small-eddy
 * (Lamont–Scott) coefficient k_L = 0.4 (ε ν)^¼ Sc^−½ and are weighted with the liquid volume of the two zones.
 * p: { vsl, vsg, rhoL, rhoG, muL, muG, sigma, D, theta, Dg (m²/s diffusivity of the gas in the liquid), rough }.
 * Returns { kLa, kLaSlug, kLaFilm, kLaStratified (same rates without slugging), enhancement, aSlug, aFilm (1/m), dBubble (m), slugFraction, liquidInSlug (fraction of the liquid inside slug bodies), freq, length, holdupSlug, holdupFilm, epsSlug (W/kg), Dg }.
 */
export function slugMassTransfer({ vsl, vsg, rhoL, rhoG, muL, muG, sigma = 0.02, D, theta = 0, Dg = 2e-9, rough = 4.5e-5 }) {
  const u = slugUnit({ vsl, vsg, rhoL, rhoG, muL, muG, D, theta }), vm = vsl + vsg, nu = muL / rhoL, A = (PI * D * D) / 4, kL = (eps) => 0.4 * (Math.max(eps, 1e-9) * nu) ** 0.25 * Math.sqrt(Dg / nu);
  const chord = (H) => D * Math.sqrt(Math.max(1 - (2 * H - 1) ** 2, 0.05)), fS = frictionFactor((rhoL * vm * D) / muL, rough / D), epsSlug = (fS * vm ** 3) / (2 * D);
  const dMax = 0.725 * (sigma / rhoL) ** 0.6 * Math.max(epsSlug, 1e-9) ** -0.4, dB = clamp(0.6 * dMax, 5e-5, 0.3 * D), HLS = u.holdupSlug, HLF = u.holdupFilm, aSlug = (6 * (1 - HLS)) / (dB * HLS);
  const vF = Math.max(vsl * 0.3, 1e-4) / HLF, fF = frictionFactor((rhoL * vF * D) / muL, rough / D), epsFilm = (fF * vF ** 3) / (2 * D), aFilm = chord(HLF) / (A * HLF);
  const b = u.slugFraction, wS = (b * HLS) / (b * HLS + (1 - b) * HLF), kLaSlug = kL(epsSlug) * aSlug, kLaFilm = kL(epsFilm) * aFilm, kLa = wS * kLaSlug + (1 - wS) * kLaFilm;
  const H0 = clamp(u.holdup, 0.02, 0.98), v0 = vsl / H0, f0 = frictionFactor((rhoL * v0 * D) / muL, rough / D), kLaStratified = (kL((f0 * v0 ** 3) / (2 * D)) * chord(H0)) / (A * H0);
  return { kLa, kLaSlug, kLaFilm, kLaStratified, enhancement: kLa / Math.max(kLaStratified, 1e-300), aSlug, aFilm, dBubble: dB, slugFraction: b, liquidInSlug: wS, freq: u.freq, length: u.length, holdupSlug: HLS, holdupFilm: HLF, epsSlug, Dg };
}
/**
 * Wax solubility curve from the solid–liquid equilibrium of the n-paraffin distribution of the oil (ideal solution, the
 * fluid suite's wax model on the single-carbon-number split of the C7+ fraction): tabulated solid fraction against
 * temperature, the wax appearance temperature at a detection threshold and the precipitation slope.
 * spec: case fluid { comp, c7MW, c7SG }; o: { detect (wt % solid that counts as the cloud point), hfMult (multiplier on the heats of fusion), n }.
 * Returns { ok, wat (°C) | null, T[] (°C ascending), solid[] (mass fraction of the oil), wTot (precipitable wax, mass fraction at the lowest temperature), waxContent (wt % C18+ wax formers), at(T) → { dissolved, solid, dCdT } }.
 */
export function waxEquilibrium(spec, { detect = 0.02, hfMult = 1, n = 48, Tmin = -30 } = {}) {
  const none = { ok: false, wat: null, T: [], solid: [], wTot: 0, waxContent: 0, at: () => ({ dissolved: 0, solid: 0, dCdT: 0 }) };
  try {
    const comp = spec?.comp || {}, tot = Object.values(comp).reduce((a, b) => a + (+b || 0), 0), z7 = tot > 0 ? (+comp.C7p || 0) / tot : 0; if (!(z7 > 0)) return none;
    const f = makeFluid(spec), fl = flashPT(f, 1.01325, 15), x = fl.phase === 'gas' ? null : fl.x; if (!x) return none;
    let x7 = 0, mw = 0; f.comps.forEach((c, i) => { if (c.pseudo) x7 += x[i]; mw += x[i] * c.MW; });
    const dist = scnDistribution(z7, clamp(num(spec.c7MW, 210), 96, 600), clamp(num(spec.c7SG, 0.84), 0.7, 1.05)), wm = waxModel(dist, x7, mw, { hfMult }), wat = wm.wat(detect);
    if (wat === null) return { ...none, ok: true, waxContent: wm.waxContent };
    const T = linspace(Tmin, Math.max(wat + 2, Tmin + 5), n), solid = T.map((t) => wm.solidWt(t) / 100), wTot = solid[0];
    const at = (t) => { if (t >= T[n - 1]) return { dissolved: wTot, solid: 0, dCdT: 0 }; const tt = Math.max(t, T[0]), s = interp1(T, solid, tt), h = 0.25, sl = (interp1(T, solid, Math.max(tt - h, T[0])) - interp1(T, solid, Math.min(tt + h, T[n - 1]))) / (Math.min(tt + h, T[n - 1]) - Math.max(tt - h, T[0])); return { dissolved: Math.max(wTot - s, 0), solid: s, dCdT: Math.max(sl, 0) }; };
    return { ok: true, wat, T, solid, wTot, waxContent: wm.waxContent, at };
  } catch { return none; }
}
/**
 * Population-balance crystallisation of a sparingly soluble salt in plug flow (method of moments for nucleation and
 * size-independent growth, with depletion of the supersaturation):
 *   dμ₀/dt = J,  dμ₁/dt = G μ₀ + J L*,  dμ₂/dt = 2 G μ₁,  dμ₃/dt = 3 G μ₂,  dx/dt = (k_v ρ_c / M)(dμ₃/dt)
 *   J = A exp(−16π σ³ v_m² / (3 (kT)³ (ν ln S)²))  (classical nucleation, ν ions per formula unit),  G = k_g (S − 1)^g
 * S is the saturation ratio (IAP/Ksp)^(1/ν); it falls as the precipitated amount x (mol/m³) removes the lattice ions.
 * o: { S0, cA, cB (mol/m³ of the two lattice ions before precipitation), nuB (stoichiometric demand on B per mole: 1 sulphates, 2 for bicarbonate),
 *      M (kg/mol), rho (kg/m³), TK, t (s), sigma (J/m²), A (1/m³/s), kg (m/s), g, seeds: { n (1/m³), L (m) }, steps }.
 * Returns { mu: [μ₀…μ₃], x (mol/m³ precipitated), S (final), d10, d32 (m), number (1/m³), massConc (kg/m³), J0, G0 (initial rates), t[], Sser[], Lser[] }.
 */
export function scaleCrystallisation({ S0, cA, cB, nuB = 1, M, rho, TK = 330, t, sigma = 0.09, A = 1e30, kg = 1e-10, g = 2, seeds = null, init = null, steps = 200, deplete = true }) {
  const vm = M / rho / 6.02214076e23, nu = 2, kv = PI / 6, Jof = (S) => { if (!(S > 1.0001)) return 0; const ex = (16 * PI * sigma ** 3 * vm * vm) / (3 * (KB * TK) ** 3 * (nu * Math.log(S)) ** 2); return ex > 700 ? 0 : A * Math.exp(-ex); }, Gof = (S) => (S > 1 ? kg * (S - 1) ** g : 0);
  const Lc = (S) => (S > 1.0001 ? (4 * sigma * vm) / (KB * TK * nu * Math.log(S)) : 0), Sof = (x) => (deplete ? S0 * (Math.max((cA - x) * Math.max(cB - nuB * x, 0) ** nuB, 0) / (cA * cB ** nuB)) ** 0.5 : S0);
  const y = init ? init.slice(0, 5) : seeds ? [seeds.n, seeds.n * seeds.L, seeds.n * seeds.L ** 2, seeds.n * seeds.L ** 3, 0] : [0, 0, 0, 0, 0], xMax = Math.min(cA, cB / nuB) * 0.999999, f = (yy) => { const S = Sof(yy[4]), J = Jof(S), Gr = Gof(S), ls = Lc(S), d3 = 3 * Gr * yy[2] + J * ls ** 3; return [J, Gr * yy[0] + J * ls, 2 * Gr * yy[1] + J * ls * ls, d3, (kv * rho * d3) / M]; };
  const xSat = !deplete ? Infinity : S0 <= 1 ? 0 : Sof(xMax) >= 1 ? xMax : brent((x) => Sof(x) - 1, 0, xMax, 1e-14 * xMax, 80);
  const xLim = Math.max(xSat, y[4]), ts = [0], Ss = [S0], Ls = [seeds ? seeds.L : 0]; let h = t / steps;
  for (let k = 0; k < steps; k++) { // classical fourth-order Runge–Kutta with a guard on the precipitable amount
    const k1 = f(y), y2 = y.map((v, i) => v + 0.5 * h * k1[i]), k2 = f(y2), y3 = y.map((v, i) => v + 0.5 * h * k2[i]), k3 = f(y3), y4 = y.map((v, i) => v + h * k3[i]), k4 = f(y4);
    for (let i = 0; i < 5; i++) y[i] += (h * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i])) / 6;
    if (deplete && y[4] > xLim) { y[3] = Math.max(y[3] - ((y[4] - xLim) * M) / (kv * rho), 0); y[4] = xLim; } // growth cannot carry the solution below saturation
    ts.push((k + 1) * h); Ss.push(Sof(y[4])); Ls.push(y[0] > 0 ? y[1] / y[0] : 0);
  }
  return { mu: y.slice(0, 4), state: y.slice(), x: y[4], S: Sof(y[4]), d10: y[0] > 0 ? y[1] / y[0] : 0, d32: y[2] > 0 ? y[3] / y[2] : 0, number: y[0], massConc: kv * rho * y[3], J0: Jof(S0), G0: Gof(S0), t: ts, Sser: Ss, Lser: Ls };
}
/**
 * Solubility parameter (MPa^½) and molar volume (m³/mol) of a phase from the cubic equation of state:
 * δ² = −U_res / V with U_res = [(a − T da/dT) / (b (d₁ − d₂))] ln[(V + d₂ b)/(V + d₁ b)] (cohesive energy density).
 */
export function eosSolubilityParameter(f, x, Pbar, TK, kind = 'liquid') {
  const e0 = eosPhase(f, x, Pbar, TK, kind), h = 0.5, dadT = (eosPhase(f, x, Pbar, TK + h, kind).a - eosPhase(f, x, Pbar, TK - h, kind).a) / (2 * h), V = (e0.Z * R * TK) / (Pbar * 1e5), { d1, d2 } = f.eos;
  const U = ((e0.a - TK * dadT) / (e0.b * (d1 - d2))) * Math.log((V + d2 * e0.b) / (V + d1 * e0.b));
  return { delta: Math.sqrt(Math.max(-U / V, 0)) / 1000, V, U };
}
/**
 * Asphaltene onset along a depressurisation path from the regular-solution / Flory–Huggins (Hirschberg) model, with the
 * live-oil solubility parameter and molar volume taken either from the equation of state (flash at each pressure) or
 * from the density correlation. Precipitation is predicted where the soluble volume fraction falls below the amount present.
 * A measured onset pressure can be honoured by adjusting the asphaltene solubility parameter (solid-phase equilibrium with a
 * reference fugacity fitted at that point); dvS adds the Poynting term of the volume change on precipitation.
 * o: { spec (case fluid), pRes, tRes, pEnd, tEnd, n, phiA (asphaltene volume fraction of the stock-tank oil), deltaA (MPa^½ at 25 °C), vA (m³/kmol),
 *      method: 'eos' | 'density', refOnset (bara, measured upper onset at reservoir temperature; 0 = none), dvS (m³/kmol volume change on precipitation), path: optional [{ P, T, rho, mw }] for the density method }.
 * Returns { P[], T[], deltaL[], vL[] (m³/kmol), phiMax[], phiA, pBubble, upperOnset, lowerOnset (bara | null), minSolubility, pMin, deltaA (used), method, precipitates }.
 */
export function asphalteneOnset({ spec, pRes, tRes, pEnd = 1.5, tEnd = null, n = 36, phiA, deltaA = 20, vA = 2, method = 'eos', refOnset = 0, dvS = 0, path = null, pre = null, anchor = true }) {
  const st = pre || asphaltenePath({ spec, pRes, tRes, pEnd, tEnd, n, method, path }), { P, T, vL } = st, N = P.length, p0 = P[0], pE = P[N - 1];
  // the cubic equation of state gives the change of the oil solubility parameter along the path; its level is anchored to the density correlation at the first point
  const shift = anchor && st.rho0 > 0 && Number.isFinite(st.dL[0]) && /^eos/.test(st.method) ? 17.347 * (st.rho0 / 1000) + 2.904 - st.dL[0] : 0, dL = st.dL.map((x) => x + shift);
  const curve = (dA25) => P.map((p, k) => { if (!Number.isFinite(dL[k])) return NaN; const TK = T[k] + KEL, dA = dA25 * (1 - 1.07e-3 * (TK - 298.15)), ex = vA / vL[k] - 1 - (vA * 1e-3 * ((dA - dL[k]) * 1e3) ** 2) / (R * TK) - (dvS * 1e-3 * (p - p0) * 1e5) / (R * TK); return Math.min(Math.exp(Math.min(ex, 0)), 1); });
  const onsets = (ph) => { let up = null, lo = null; for (let k = 1; k < N; k++) { const a = ph[k - 1] - phiA, b = ph[k] - phiA; if (!Number.isFinite(a) || !Number.isFinite(b)) continue; if (a >= 0 && b < 0 && up === null) up = P[k - 1] + ((P[k] - P[k - 1]) * a) / (a - b); if (a < 0 && b >= 0) lo = P[k - 1] + ((P[k] - P[k - 1]) * a) / (a - b); } if (up === null && ph[0] < phiA) up = P[0]; return { up, lo }; };
  let dUse = deltaA;
  if (refOnset > pE && refOnset < p0) { const g = (d) => (onsets(curve(d)).up ?? pE) - refOnset; try { if (g(12) * g(30) < 0) dUse = brent(g, 12, 30, 1e-6, 80); } catch { dUse = deltaA; } } // a larger δ_a makes the asphaltene less soluble and raises the onset pressure
  const ph = curve(dUse), o = onsets(ph), fin = ph.filter(Number.isFinite), mn = fin.length ? Math.min(...fin) : 1, kMin = ph.indexOf(mn);
  return { P, T, deltaL: dL, vL, phiMax: ph, phiA, pBubble: st.pBubble, upperOnset: o.up, lowerOnset: o.lo, minSolubility: mn, pMin: kMin >= 0 ? P[kMin] : null, deltaA: dUse, method: st.method, precipitates: mn < phiA, shift };
}
const aspMemo = new Map();
/**
 * Solubility parameter and molar volume of the live oil along a depressurisation path (the expensive part of the onset
 * model; cached). method 'eos': flash and cohesive energy from the equation of state; 'density': δ = 17.347 ρ + 2.904 on the
 * supplied path [{ rho (kg/m³), mw (g/mol), wG }]. Returns { P[], T[], dL[] (MPa^½), vL[] (m³/kmol), pBubble, method }.
 */
export function asphaltenePath({ spec, pRes, tRes, pEnd = 1.5, tEnd = null, n = 36, method = 'eos', path = null }) {
  const key = method === 'eos' && !path ? JSON.stringify([spec?.comp, spec?.c7MW, spec?.c7SG, spec?.eos, spec?.nPseudo, pRes, tRes, pEnd, tEnd, n]) : null; if (key && aspMemo.has(key)) return aspMemo.get(key);
  const P = linspace(pRes, pEnd, n), T = P.map((_, k) => tRes + (((tEnd ?? tRes) - tRes) * k) / (n - 1)), dL = [], vL = []; let pBub = null, used = method, f = null, rho0 = 0;
  if (method === 'eos') { try { f = makeFluid(spec); } catch { f = null; } }
  for (let k = 0; k < n; k++) {
    let ok = false;
    if (f) { try { const fl = flashPT(f, P[k], T[k]); if (fl.phase !== 'gas') { const x = fl.phase === 'two' ? fl.x : f.z, e = eosSolubilityParameter(f, x, P[k], T[k] + KEL, 'liquid'); if (e.delta > 5 && e.delta < 30) { dL.push(e.delta); vL.push(e.V * 1e3); ok = true; if (k === 0) rho0 = (f.comps.reduce((s, c, q) => s + x[q] * c.MW, 0) * 1e-3) / e.V; if (fl.phase === 'two' && pBub === null) pBub = k > 0 ? 0.5 * (P[k] + P[k - 1]) : P[k]; } } } catch { ok = false; } }
    if (!ok) { const q = path ? path[Math.min(k, path.length - 1)] : null; if (!q) { dL.push(NaN); vL.push(NaN); continue; } if (f) used = 'eos with density fallback'; else used = 'density'; dL.push(17.347 * (q.rho / 1000) + 2.904); vL.push(q.mw / q.rho); if (pBub === null && q.wG > 1e-4) pBub = P[k]; }
  }
  const out = { P, T, dL, vL, pBubble: pBub, method: used, rho0 };
  if (key) { if (aspMemo.size > 8) aspMemo.clear(); aspMemo.set(key, out); }
  return out;
}

// =====================================================================================================
// 11. Case set-up on the solids grid
// =====================================================================================================
const WETGAS = { name: 'Lean wet gas', comp: { N2: 1, CO2: 2, H2S: 0, C1: 86, C2: 6, C3: 3, iC4: 0.5, nC4: 0.8, iC5: 0.2, nC5: 0.2, C6: 0.2, C7p: 0.1 }, c7MW: 120, c7SG: 0.76, rateBasis: 'gas', wc: 0 };
const C_STEEL = 470; // J/kg/K

/**
 * Build the line, fluid and reference flow pictures on n equal cells. Returns the set-up object used by the marching
 * solvers: geometry arrays, props(P, T) (memoised kernel properties at the case rate), prof(frac) (steady picture at a
 * rate fraction), grad(...) (kernel pressure gradient), the hydrate curve and its inverse, thermal mass per metre.
 */
export function buildSetup(v, ctx, n) {
  const c = { fluid: ctx?.fluid, site: ctx?.site, outputs: ctx?.outputs };
  const fo = v.fluidSystem === 'wetgas' ? { ...WETGAS, qGas: v.gasRate, qWater: v.gasWater } : v.fluidSystem === 'highwc' ? { rateBasis: 'oil', wc: v.highWc } : null;
  const line0 = caseLine(c), over = { id: v.idMm / 1000, roughness: v.roughUm * 1e-6, uValue: v.uValue, tSeabed: v.tSeabed };
  const differs = (a, b) => Math.abs(a - b) > 1e-6 * Math.max(Math.abs(b), 1e-9);
  const custom = !!fo || differs(over.id, line0.id) || differs(over.roughness, line0.roughness) || differs(over.uValue, line0.uValue) || differs(over.tSeabed, line0.tSeabed) || differs(v.lengthScale, 1);
  if (differs(v.lengthScale, 1)) over.profile = { x: line0.profile.x.map((x) => x * v.lengthScale), z: line0.profile.z.slice() };
  if (fo) over.override = fo;
  const fm = fluidModel(c, fo || {}), line = caseLine(c, over), D0 = line.id;
  const toCells = (pic) => {
    const sp = [0]; for (let i = 1; i < pic.x.length; i++) sp.push(sp[i - 1] + Math.max(Math.hypot(pic.x[i] - pic.x[i - 1], pic.z[i] - pic.z[i - 1]), 1e-9));
    const L = sp[sp.length - 1], ds = L / n, sc = Array.from({ length: n }, (_, i) => (i + 0.5) * ds), f = (a) => sc.map((s) => interp1(sp, a, s)), ze = Array.from({ length: n + 1 }, (_, i) => interp1(sp, pic.z, i * ds));
    return { L, ds, s: sc, x: f(pic.x), z: f(pic.z), theta: ze.slice(1).map((z, i) => Math.asin(clamp((z - ze[i]) / ds, -1, 1))), P: f(pic.P), T: f(pic.T), holdup: f(pic.holdup).map((h) => clamp(h, 0.01, 1)), dpdx: f(pic.dpdx), tauW: f(pic.tauW).map((t) => Math.max(t, 0)), rhoM: f(pic.rhoM), tAmb: f(pic.tAmb),
      regime: sc.map((s) => { let k = 0; while (k < sp.length - 1 && sp[k + 1] < s) k++; return String(pic.regime[Math.min(k, pic.regime.length - 1)] ?? ''); }), pOut: pic.P[pic.P.length - 1], pIn: pic.P[0], src: pic.source || 'kernel estimate' };
  };
  const memo = new Map();
  const prof = (frac) => {
    const key = Math.round(frac * 1e4);
    if (!memo.has(key)) {
      let pic;
      try { pic = key === 10000 ? flowPicture(c, custom ? { ...over, force: true } : {}) : steadyCase(c, { ...over, mScale: frac }); }
      catch (e) { throw new Error(`No steady flow solution at ${(frac * 100).toFixed(0)} % of the case rate: ${e.message}`); }
      memo.set(key, toCells(pic));
    }
    return memo.get(key);
  };
  const base = prof(1), pc = new Map();
  const props = (P, T) => {
    const kp = Math.round(48 * Math.log(clamp(P, 1, 900))), kt = Math.round(clamp(T, -30, 170)), key = kt * 4096 + kp;
    let o = pc.get(key); if (!o) { o = fm.at(Math.exp(kp / 48), kt); pc.set(key, o); } return o;
  };
  const inj = (Array.isArray(v.injections) ? v.injections : []).map((r) => ({ x: num(r?.x, 0) * 1000, q: Math.max(num(r?.rate, 0), 0) })).filter((r) => r.q > 0);
  const aqS = fm.aq.S, inhId0 = v.inhibitor === 'case' ? fm.aq.inhId : v.inhibitor, inhId = inhId0 === 'none' && inj.length ? v.injChem || 'MEG' : inhId0, inh = INHIBITORS[inhId] || INHIBITORS.none, inhWt = inhId0 === 'none' ? 0 : v.inhibitor === 'case' ? fm.aq.inhWt : v.inhWt;
  const eqOff = num(v.hydEqOffset, 0), Pg = logspace(1, 700, 90), Tg = Pg.map((p) => fm.hydrateT0(p) + eqOff); for (let i = 1; i < Tg.length; i++) if (Tg[i] <= Tg[i - 1]) Tg[i] = Tg[i - 1] + 1e-9;
  const lnPg = Pg.map(Math.log), wt = line.wt;
  const S = { n, ds: base.ds, L: base.L, s: base.s, x: base.x, z: base.z, theta: base.theta, tAmb: base.tAmb, D0, rough0: line.roughness, U: line.uValue, wt, fm, line, props, prof, custom, src: base.src,
    hT0: (P) => fm.hydrateT0(P) + eqOff, peq: (Tfresh) => Math.exp(interp1(Tg, lnPg, Tfresh)), aq: { S: aqS, inh, inhId, inhWt }, pOut: base.pOut };
  { // inhibitor concentration along the line: the dose in the produced water plus every injection point upstream of the cell (at the rate of the event)
    const ph = scenarioPhases(v).phases, fEv = ph[ph.length - 1].frac || 1, mW = Math.max((fm.rates?.mW ?? 0) * fEv, 1e-9), w0 = clamp(inhWt / 100, 0, 0.97), m0 = (w0 / (1 - w0)) * mW;
    S.aq.inhWtCell = S.x.map((x) => { let mi = m0; for (const q of inj) if (q.x <= x + 1e-6) mi += (inh.rho * q.q) / 86400; return (100 * mi) / (mi + mW); }); S.aq.injections = inj; S.aq.injected = inj.reduce((a, q) => a + q.q, 0);
  }
  const sm = new Map();
  S.slugMT = (i, frac) => { const key = i * 20000 + Math.round(frac * 1e4); if (!sm.has(key)) { const r = prof(frac), pr = props(r.P[i], r.T[i]), A = (PI * D0 * D0) / 4, oilC = pr.wcut <= 0.6, muC = oilC ? pr.muO : pr.muW; sm.set(key, slugMassTransfer({ vsl: (pr.qL * frac) / A, vsg: (pr.qG * frac) / A, rhoL: pr.rhoL, rhoG: pr.rhoG, muL: pr.muL, muG: pr.muG, sigma: pr.sigma, D: D0, theta: S.theta[i], rough: S.rough0, Dg: (7.4e-12 * Math.sqrt(oilC ? pr.mwO || 150 : 46.8) * (r.T[i] + KEL)) / (muC * 1000 * 37.7 ** 0.6) })); } return sm.get(key); };
  S.grad = (i, P, T, D, rough, muFac, frac) => { const pr = props(P, T), A = (PI * D * D) / 4; return gradient({ vsl: (pr.qL * frac) / A, vsg: (pr.qG * frac) / A, rhoL: pr.rhoL, rhoG: pr.rhoG, muL: pr.muL * muFac, muG: pr.muG, sigma: pr.sigma, D, theta: S.theta[i], rough, P: P * 1e5 }); };
  const gm = new Map();
  S.gref = (frac) => { const key = Math.round(frac * 1e4); if (!gm.has(key)) { const r = prof(frac); gm.set(key, r.P.map((p, i) => S.grad(i, p, r.T[i], D0, S.rough0, 1, frac))); } return gm.get(key); };
  S.C = base.P.map((p, i) => { if (v.thermalMass > 0) return v.thermalMass * 1000; const pr = props(p, base.T[i]), A = (PI * D0 * D0) / 4, H = base.holdup[i]; return A * (H * pr.rhoL * pr.cpL + (1 - H) * pr.rhoG * pr.cpG) + BASE.rhoSteel * C_STEEL * PI * (D0 + wt) * wt; });
  return S;
}
/** Uniform laboratory-style set-up with constant fluid properties (used by the verification cases and the calibration model). */
export function labSetup({ n = 20, L = 2000, D = 0.1, T = 4, P = 80, U = 0, tAmb = 4, vsl = 1, vsg = 1, wcut = 0.2, sg = 0.7, C = 6e4 } = {}) {
  const A = (PI * D * D) / 4, qL = vsl * A, qG = vsg * A, rhoO = 800, rhoW = 1000, rhoG = 70, qW = qL * wcut, qO = qL - qW, rhoL = (qO * rhoO + qW * rhoW) / qL;
  const pr = { rhoG, rhoO, rhoW, rhoL, muG: 1.3e-5, muO: 3e-3, muW: 1.5e-3, muL: 3e-3 * (1 - Math.min(wcut, 0.7)) ** -2.5, cpG: 2500, cpO: 2100, cpW: 4100, cpL: (qO * rhoO * 2100 + qW * rhoW * 4100) / (qL * rhoL), kG: 0.04, kO: 0.14, kW: 0.58, kL: 0.2, sigma: 0.02, sigmaOW: 0.03, zG: 0.85, mwG: 18, mwO: 150, mG: qG * rhoG, mO: qO * rhoO, mW: qW * rhoW, qG, qO, qW, qL, wcut, phaseInv: wcut > 0.6 };
  const ds = L / n, arr = (x) => new Array(n).fill(x), Pg = logspace(1, 700, 90), Tg = Pg.map((p) => hydrateT0(p, sg)), lnPg = Pg.map(Math.log);
  const S = { n, ds, L, s: arr(0).map((_, i) => (i + 0.5) * ds), x: arr(0).map((_, i) => (i + 0.5) * ds), z: arr(0), theta: arr(0), tAmb: arr(tAmb), D0: D, rough0: 4.5e-5, U, wt: 0.008, props: () => pr, custom: true, src: 'laboratory set-up',
    hT0: (p) => hydrateT0(p, sg), peq: (Tf) => Math.exp(interp1(Tg, lnPg, Tf)), aq: { S: 0, inh: INHIBITORS.none, inhId: 'none', inhWt: 0 }, pOut: P, C: arr(C) };
  S.grad = (i, Pp, Tt, Dd, rough, muFac, frac) => { const a = (PI * Dd * Dd) / 4; return gradient({ vsl: (qL * frac) / a, vsg: (qG * frac) / a, rhoL, rhoG, muL: pr.muL * muFac, muG: pr.muG, sigma: pr.sigma, D: Dd, theta: 0, rough, P: Pp * 1e5 }); };
  const g1 = S.grad(0, P, T, D, 4.5e-5, 1, 1), pm = new Map();
  S.prof = (frac) => { const k = Math.round(frac * 1e4); if (!pm.has(k)) { const g = S.grad(0, P, T, D, 4.5e-5, 1, frac); pm.set(k, { L, ds, P: arr(P), T: arr(T), holdup: arr(g.holdup), dpdx: arr(0), tauW: arr(g.tauW), rhoM: arr(rhoL * g.holdup + rhoG * (1 - g.holdup)), tAmb: arr(tAmb), regime: arr(g.regime || 'slug'), pOut: P, pIn: P }); } return pm.get(k); };
  S.gref = (frac) => arr(S.grad(0, P, T, D, 4.5e-5, 1, frac)); S.g1 = g1; S.aq.inhWtCell = arr(0);
  S.slugMT = (i, frac) => slugMassTransfer({ vsl: vsl * frac, vsg: vsg * frac, rhoL, rhoG, muL: pr.muL, muG: pr.muG, sigma: pr.sigma, D, theta: 0, Dg: (7.4e-12 * Math.sqrt(150) * (T + KEL)) / (pr.muO * 1000 * 37.7 ** 0.6) });
  return S;
}

// =====================================================================================================
// 12. Transient hydrate march: transport, kinetics, population balance, deposition, bore feedback
// =====================================================================================================
/** Kinetic and deposition parameters of the hydrate march from the input values (SI). */
export function hydrateParams(v, S, over = {}) {
  const b = S.prof(1), mid = Math.floor(S.n / 2), pr = S.props(b.P[mid], b.T[mid]), lamL = pr.qL / Math.max(pr.qL + pr.qG, 1e-12);
  const mode = v.regime && v.regime !== 'auto' ? v.regime : lamL < 0.1 ? 'gas' : pr.wcut > 0.6 ? 'water' : 'oil';
  // primary particle: Boxall inertial droplet correlation d/D = 0.063 We^-3/5 on the mixture velocity
  const vm = (pr.qL + pr.qG) / ((PI * S.D0 * S.D0) / 4), We = ((mode === 'water' ? pr.rhoW : pr.rhoO) * vm * vm * S.D0) / Math.max(mode === 'oil' ? pr.sigmaOW : pr.sigma, 1e-4);
  const dAuto = clamp(0.063 * We ** -0.6 * S.D0, 10e-6, 400e-6), dPrim = v.primaryUm > 0 ? v.primaryUm * 1e-6 : dAuto;
  let hydN = v.hydNumber;
  if (!(hydN > 0)) { const Tm = mean(S.tAmb) + KEL, Pm = mean(b.P); hydN = clamp(langmuirOccupancy(Math.max(Tm, 273.2), methaneFugacity(Pm, Math.max(Tm, 273.2)), v.langmuirModel === 'kihara' ? 'kihara' : 'parrish').hydrationNumber, 5.75, 7.5); }
  const mwG = clamp(pr.mwG || 18, 16, 30) * 1e-3, mHyd = mwG + hydN * MW_W;
  const rhoH = v.rhoHyd > 0 ? v.rhoHyd : HYDRATE.rho;
  return { mode, dPrim, dAuto, hydN, mwG, rhoH, wfH: (hydN * MW_W) / mHyd, vmh: mHyd / rhoH, dHmol: HYDRATE.latent * mHyd,
    kinK: v.kinK * 1e-10, EaR: v.kinEa, shellD: v.shellD * 1e-15, H: mode === 'oil' ? 2500 : 7e4, mtMult: v.mtMult,
    nucModel: v.nucleation === 'heterogeneous' ? 'measured' : 'classical', nucMult: 10 ** v.nucA, nucB: v.nucB, nucA: 10 ** (v.nucA + (v.nucleation === 'homogeneous' ? 35 : 7.5)), theta: v.contactAngle, sigma: v.sigmaHW * 1e-3, nucV: v.nucVolume * 1e-3, het: v.nucleation !== 'homogeneous', lamStar: 1,
    cohesion: v.cohesion * 1e-3, fr: v.fractal, phiMax: v.phiMax, viscModel: v.viscModel, aggEff: v.aggEff, kBreak: v.kBreak,
    inhEff: v.inhEff / 100, htMult: v.htMult, inletPhi: v.inletHydPct / 100, disK0: v.disK0 * 1e4, disE: v.disE * 1000, sigG: Math.max(v.primarySigma, 1),
    adhesion: v.adhesion, adhForce: v.adhForce * 1e-3, tauCrit: v.tauCrit, kRemove: v.kRemove / 3600, por0: v.porosity0, porInf: Math.min(v.porosityInf, v.porosity0), tAge: v.ageHours * 3600, filmMult: v.filmMult,
    plugBlock: v.plugBlockPct / 100, pInMax: v.pInMax, pShut: v.pShut, maxSub: 4, slugCouple: !!v.slugCouple, gasSat: clamp(num(v.gasWaterSatPct, 100) / 100, 0, 1.5), ...over };
}
/**
 * Transient hydrate solver on the line (a stepper: call step() until done, then result()).
 * Fluid-borne quantities (particle numbers per size class, primary-particle count, water removed to deposits, the
 * nucleation hazard integral) are carried as cell inventories with implicit upwind transport; in each cell the kinetics
 * act on the mixed state: nucleation hazard → onset → growth (resistances in series, limited by water, gas and the heat
 * that can be removed) → aggregation and breakage (sectional population balance) → wall capture, vapour-film growth,
 * shear removal and ageing of the deposit. The bore, the slurry viscosity and the wall shear feed back on the hydraulics.
 * o: { phases: [{ dur (s), dt (s), frac }], grid (pbeGrid), Dbase[] (bore before hydrate, m), roughBase[], dep0 (m), fracInit, cheap, rows }
 */
export function hydrateMarch(S, p, o) {
  const n = S.n, ds = S.ds, g = o.grid, nC = g.n, D0 = S.D0, RHO = p.rhoH, LAT = HYDRATE.latent, vp = (PI / 6) * p.dPrim ** 3, eexp = (S.props(S.pOut, 4).rhoW / RHO) / p.wfH; // hydrate volume per unit water volume
  const NP = Array.from({ length: n }, () => new Float64Array(nC)), PS = new Float64Array(n), WD = new Float64Array(n), LV = new Float64Array(n), HV = new Float64Array(n), film = new Float64Array(n), mDep = new Float64Array(n), por = new Float64Array(n).fill(p.por0);
  const Dbase = o.Dbase || new Array(n).fill(D0), roughBase = o.roughBase || new Array(n).fill(S.rough0), init = S.prof(o.fracInit ?? 1);
  for (let i = 0; i < n; i++) if (o.dep0 > 0) mDep[i] = RHO * (1 - p.por0) * (PI / 4) * (Dbase[i] ** 2 - Math.max(Dbase[i] - 2 * o.dep0, 0.05 * D0) ** 2) * ds;
  const mDep0 = mDep.reduce((a, b) => a + b, 0);
  // optional boundary and initial data: inhibitor per cell (mid-line injection), heated section, depressurisation, equipment traps, water distribution at rest
  const inhC = o.inhWt || null, heat = o.heat || null, heatT0 = o.heatStart ?? 0, depr = o.depress || null, trapEff = new Float64Array(n), wMult = o.wMult || null, gasSat = p.gasSat ?? 1, slugK = new Array(n).fill(null);
  for (const tr of o.traps || []) if (tr.i >= 0 && tr.i < n) trapEff[tr.i] = clamp(tr.eff, 0, 1);
  const T = init.T.slice(), P = init.P.slice(), hold = init.holdup.slice(), dpdx = init.dpdx.slice(), tauW = init.tauW.slice(), muFac = new Array(n).fill(1), Dn = Dbase.slice();
  const cD = new Array(n).fill(-1), cMu = new Array(n).fill(1), cT = new Array(n).fill(0), cP = new Array(n).fill(0), cF = new Array(n).fill(0), cG = new Array(n).fill(null);
  const rec = { sub: new Array(n).fill(0), teq: new Array(n).fill(0), phi: new Array(n).fill(0), phiE: new Array(n).fill(0), d43: new Array(n).fill(0), dA: new Array(n).fill(0), J: new Array(n).fill(0), lam: new Array(n).fill(0), rate: new Array(n).fill(0), vL: new Array(n).fill(0), X: new Array(n).fill(0), tw: new Array(n).fill(0), lim: new Array(n).fill(''), depRate: new Array(n).fill(0), freeW: new Array(n).fill(1), subMax: new Array(n).fill(-99), phiMaxT: new Array(n).fill(0), Jmax: 0, removalMax: 0, captureMax: 0, VL: new Array(n).fill(1), gdot: new Array(n).fill(0), muC: new Array(n).fill(1e-3), rhoC: new Array(n).fill(800), exposure: new Array(n).fill(0), removal: new Array(n).fill(0), capture: new Array(n).fill(0) };
  const ser = { exp: [], t: [], pIn: [], dp: [], susp: [], dep: [], blk: [], phi: [], rate: [], sub: [], frac: [], visc: [] }, fld = { t: [], phi: [], dep: [], sub: [] };
  const cpC = new Array(n).fill(null), L3 = g.L.map((L) => L ** 3), L4 = g.L.map((L) => L ** 4);
  const beta = new Float64Array(3), gv = new Float64Array(nC), Sb = new Float64Array(nC), FinN = new Float64Array(nC);
  const peak = { phi: 0, blk: 0, tPhi: 0, tBlk: 0, iPhi: 0, iBlk: 0, N: null, phiX: null, dHyd: null };
  const led = { inflow: 0, formed: 0, dissociated: 0, exported: 0, sloughed: 0, captured: 0, filmWall: 0, heat: 0, trapped: 0, heatIn: 0, formedSlug: 0, susp0: 0 };
  const total = o.phases.reduce((a, b) => a + b.dur, 0), nSteps = o.phases.reduce((a, b) => a + Math.ceil(b.dur / b.dt - 1e-9), 0), stride = Math.max(1, Math.ceil(nSteps / (o.rows || 48)));
  let ph = 0, tPh = 0, t = 0, k = 0, done = nSteps === 0, ref = null, gRef = null, fe = 0, pst = null, tm0 = 1, onset = null, plug = null, limited = false, pInNow = P[0], released = true, newPhase = true;
  const dep0Of = (cf, i = -1) => { const Sx = Math.min(S.aq.S * cf, 260), salt = hydrateDepression({ S: Sx, inhWt: 0, inh: S.aq.inh }), wt = inhC && i >= 0 ? inhC[i] : S.aq.inhWt; if (!(wt > 0)) return salt; const w = wt / 100; return salt + (p.inhEff ?? 1) * (hydrateDepression({ S: Sx, inhWt: (100 * w) / (w + (1 - w) / cf), inh: S.aq.inh }) - salt); };
  const depBase = Array.from({ length: n }, (_, i) => dep0Of(1, i)), kStarAt = (TK) => p.kinK * Math.exp(-p.EaR * (1 / TK - 1 / 277.15));
  const kP = (() => { let b = 0; for (let j = 1; j < nC; j++) if (Math.abs(Math.log(g.L[j] / p.dPrim)) < Math.abs(Math.log(g.L[b] / p.dPrim))) b = j; return b; })();
  const seedW = (() => { const w = new Float64Array(nC), lg = Math.log(p.sigG || 1); if (!(lg > 1e-6)) { w[kP] = 1; return w; } let sm = 0; for (let c = 0; c < nC; c++) { w[c] = Math.exp(-0.5 * (Math.log(g.L[c] / p.dPrim) / lg) ** 2); sm += w[c]; } for (let c = 0; c < nC; c++) w[c] /= sm; return w; })(); // share of the seeded particle volume per class
  const kSl = (L) => { let b = 0; for (let j = 1; j < nC; j++) if (Math.abs(Math.log(g.L[j] / L)) < Math.abs(Math.log(g.L[b] / L))) b = j; return b; };
  if (o.init && o.init.phi > 0) { // initial suspended hydrate: a log-normal population of particles or agglomerates of size d in every cell
    const dI = o.init.d > 0 ? o.init.d : p.dPrim, lg = Math.log(Math.max(o.init.sig || 1, 1)), w = new Float64Array(nC); let sm = 0;
    if (lg > 1e-6) for (let c = 0; c < nC; c++) { w[c] = Math.exp(-0.5 * (Math.log(g.L[c] / dI) / lg) ** 2); sm += w[c]; } else { w[kSl(dI)] = 1; sm = 1; }
    for (let i = 0; i < n; i++) { const hv = o.init.phi * (PI / 4) * Dbase[i] ** 2 * hold[i] * ds; for (let c = 0; c < nC; c++) NP[i][c] = (hv * w[c]) / sm / g.v[c]; HV[i] = hv; PS[i] = hv / (vp * eexp); led.susp0 += hv * RHO; }
  }

  function step() {
    if (done) return;
    const phz = o.phases[ph], dt = Math.min(phz.dt, phz.dur - tPh), flowing = phz.frac > 0;
    if (newPhase) {
      newPhase = false;
      if (flowing) { fe = phz.frac; ref = S.prof(phz.frac); gRef = S.gref(phz.frac); released = false; cD.fill(-1); for (let i = 0; i < n; i++) { hold[i] = ref.holdup[i]; dpdx[i] = ref.dpdx[i]; tauW[i] = ref.tauW[i]; slugK[i] = p.slugCouple && S.slugMT && /slug/i.test(ref.regime[i]) ? S.slugMT(i, phz.frac) : null; } }
      else { // static pressure: outlet pressure plus the head of the settled column, later scaled with the absolute gas temperature
        const b = S.prof(1), rm = mean(b.rhoM); pst = new Array(n); pst[n - 1] = p.pShut > 0 ? p.pShut : P[n - 1];
        for (let i = n - 2; i >= 0; i--) pst[i] = pst[i + 1] + (rm * G * (S.z[i + 1] - S.z[i])) / 1e5;
        tm0 = mean(T) + KEL;
      }
    }
    // ---- pressure along the line
    if (flowing) {
      P[n - 1] = S.pOut + (dpdx[n - 1] * ds) / 2e5; for (let i = n - 2; i >= 0; i--) P[i] = P[i + 1] + ((dpdx[i] + dpdx[i + 1]) * ds) / 2e5;
      pInNow = P[0] + (dpdx[0] * ds) / 2e5;
      if (pInNow > p.pInMax) { // the source cannot push harder: the rate falls and the profile is capped at the available pressure
        const r = (p.pInMax - S.pOut) / (pInNow - S.pOut); fe *= clamp(Math.sqrt(Math.max(r, 0)), 0.3, 1);
        for (let i = 0; i < n; i++) P[i] = S.pOut + (P[i] - S.pOut) * r; pInNow = p.pInMax;
        if (fe < 0.05 * phz.frac && !plug) { let im = 0; for (let i = 1; i < n; i++) if (dpdx[i] > dpdx[im]) im = i; plug = { t: t + dt, x: S.x[im], i: im, mech: 'flow stalled at the available inlet pressure' }; }
      }
      for (let i = 0; i < n; i++) P[i] = clamp(P[i], 1.05, 1500);
    } else { const f = (mean(T) + KEL) / tm0, dp = depr && t >= depr.t ? depr.P / Math.max(pst[n - 1] * f, 1e-9) : 1; for (let i = 0; i < n; i++) P[i] = Math.max(pst[i] * f * Math.min(dp, 1), 1.05); pInNow = P[0]; } // a depressurisation boundary lowers the whole settled column in proportion
    // ---- march downstream through the cells
    let thUp = 0, FinPS = 0, FinWD = 0, FinLV = 0, FinHV = 0, hydIn = 0, sumRate = 0, sumSusp = 0, maxBlk = 0, iBlk = 0, iPlugDep = -1, maxPhi = 0, iPhi = 0, maxVisc = 1, iVisc = 0, maxSub = -99;
    FinN.fill(0);
    if (flowing && p.inletPhi > 0) { const qin = S.props(P[0], T[0]).qL * fe * p.inletPhi; for (let c = 0; c < nC; c++) FinN[c] = (qin * seedW[c]) / g.v[c]; FinHV = qin; FinPS = qin / (vp * eexp); led.inflow += qin * RHO * dt; }
    for (let i = 0; i < n; i++) {
      const pr0 = S.props(P[i], T[i]), mcp = flowing ? fe * (pr0.mG * pr0.cpG + pr0.mO * pr0.cpO + pr0.mW * pr0.cpW) : 0, Cds = (S.C[i] * ds) / dt, UA = S.U * PI * D0 * ds, den = Cds + mcp + UA;
      const qh = heat && t >= heatT0 ? heat[i] * ds : 0; if (qh) led.heatIn += qh * dt; // direct heating of the section (W)
      let Ti = flowing ? ref.T[i] + ((T[i] - ref.T[i]) * Cds + mcp * thUp + qh) / den : (T[i] * Cds + UA * S.tAmb[i] + qh) / den;
      const pr = S.props(P[i], Ti), Adep = mDep[i] / (RHO * (1 - por[i]) * ds), D = Math.sqrt(Math.max(Dbase[i] ** 2 - (4 * Adep) / PI, (0.03 * D0) ** 2)), A = (PI * D * D) / 4, dHyd = (Dbase[i] - D) / 2;
      Dn[i] = D;
      if (flowing) {
        const mu = muFac[i] * pr.muL;
        if (o.cheap) { const sc = (D0 / D) ** 4.8 * (mu / S.props(ref.P[i], ref.T[i]).muL) ** 0.2 * (fe / phz.frac) ** 1.8, fr0 = Math.max(gRef[i].fric, 1e-6); dpdx[i] = ref.dpdx[i] + fr0 * (sc - 1); tauW[i] = Math.max(ref.tauW[i], gRef[i].tauW) * sc * (D / D0); }
        else {
          // kernel gradient at anchor states; between anchors the friction follows the Blasius scaling in bore, viscosity and rate
          if (cD[i] < 0 || Math.abs(D / cD[i] - 1) > 0.06 || mu / cMu[i] > 2.5 || mu / cMu[i] < 0.4 || Math.abs(Ti - cT[i]) > 8 || Math.abs(P[i] / cP[i] - 1) > 0.15 || Math.abs(fe / cF[i] - 1) > 0.15) {
            const gk = S.grad(i, P[i], Ti, D, Math.min(roughBase[i] + 0.1 * dHyd, 0.05 * D), muFac[i], fe), g0 = gRef[i];
            cD[i] = D; cMu[i] = mu; cT[i] = Ti; cP[i] = P[i]; cF[i] = fe; cG[i] = { dp: ref.dpdx[i] + gk.dpdx - g0.dpdx, fric: Math.max(gk.fric, 0), tau: g0.tauW > 1e-9 ? (Math.max(ref.tauW[i], 1e-9) * gk.tauW) / g0.tauW : gk.tauW };
            hold[i] = clamp(ref.holdup[i] + gk.holdup - g0.holdup, 0.02, 1);
          }
          const sc = (cD[i] / D) ** 4.8 * (mu / cMu[i]) ** 0.2 * (fe / cF[i]) ** 1.8, an = cG[i];
          dpdx[i] = an.dp + an.fric * (sc - 1); tauW[i] = an.tau * sc * (D / cD[i]);
        }
      } else tauW[i] = 0;
      const HL = hold[i], VL = Math.max(A * HL * ds, 1e-12), qL = pr.qL * fe, tau = flowing ? VL / Math.max(qL, 1e-15) : Infinity, a = flowing ? dt / tau : 0, inv = 1 / (1 + a), h = dt * inv, q = 1 + a;
      const vsl = flowing ? qL / A : 0, vsg = flowing ? (pr.qG * fe) / A : 0, vm = vsl + vsg, vL = vsl / HL, wcut = pr.wcut, mW0 = pr.rhoW * wcut * VL * (!flowing && wMult ? wMult[i] : 1);
      const Nc = NP[i], oil = p.mode === 'oil';
      let pVol = 0; // particle volume carried in the cell
      for (let c = 0; c < nC; c++) { Nc[c] = (Nc[c] + dt * FinN[c]) * inv; pVol += Nc[c] * g.v[c]; }
      // hydrate volume: a transported scalar for shelled droplets (oil-continuous), the particle volume itself for crystals
      let hydVol = oil ? (HV[i] + dt * FinHV) * inv : pVol, ps = (PS[i] + dt * FinPS) * inv, wd = (WD[i] + dt * FinWD) * inv, lv = (LV[i] + dt * FinLV) * inv;
      if (flowing && !released && film[i] > 0) { const fv = film[i] / q, c = kSl(Math.max(film[i] / (D * ds), p.dPrim)); Nc[c] += fv / g.v[c]; pVol += fv; ps += fv / (vp * eexp); hydVol += fv; film[i] = 0; } // the shut-in interface film breaks up into the stream
      const volBefore = hydVol, mWfree = Math.max(mW0 - (hydVol + film[i]) * RHO * p.wfH - wd, 0), cf = mW0 > 0 ? mW0 / Math.max(mWfree, 0.03 * mW0) : 1;
      const dep = cf > 1.0005 ? dep0Of(cf, i) : depBase[i], Teq = S.hT0(P[i]) - dep, dTs = Teq - Ti, TK = Ti + KEL, zG = pr.zG || 0.85;
      // nucleation hazard: measured rate per unit of gas–liquid interface of the cell, or classical theory per unit of water volume
      const meas = p.nucModel === 'measured', Jn = dTs > 0 ? (meas ? hydrateNucleationRate({ dT: dTs, TeqK: Teq + KEL, mult: p.nucMult, B1: p.nucB }) : nucleationRate({ TK, dT: dTs, TeqK: Teq + KEL, sigma: p.sigma, theta: p.theta, A: p.nucA, het: p.het }).J) : 0;
      if (mWfree > 0 && dTs > 0) lv += h * Jn * (meas ? D * Math.sqrt(Math.max(1 - (2 * HL - 1) ** 2, 0.05)) * ds : p.nucV) * VL;
      if (!(dTs > 0) && pVol <= 0 && hydVol <= 0 && film[i] <= 0 && mDep[i] <= 0) { // nothing can happen here: warm, clean and particle-free
        PS[i] = ps; WD[i] = wd; LV[i] = lv; HV[i] = 0; T[i] = Ti; muFac[i] = 1;
        if (flowing) { for (let c = 0; c < nC; c++) FinN[c] = 0; FinPS = ps / tau; FinWD = wd / tau; FinLV = lv / tau; FinHV = 0; hydIn = 0; thUp = Ti - ref.T[i]; }
        rec.sub[i] = dTs; rec.teq[i] = Teq; rec.phi[i] = 0; rec.phiE[i] = 0; rec.d43[i] = 0; rec.J[i] = 0; rec.lam[i] = lv / VL; rec.rate[i] = 0; rec.vL[i] = vL; rec.lim[i] = 'outside the hydrate region'; rec.freeW[i] = mW0 > 0 ? mWfree / mW0 : 0; rec.X[i] = 0; rec.capture[i] = 0; rec.removal[i] = 0; rec.depRate[i] = 0;
        if (dTs > rec.subMax[i]) rec.subMax[i] = dTs; if (dTs > maxSub) maxSub = dTs;
        { const blk = 1 - (D * D) / (D0 * D0); if (blk > maxBlk) { maxBlk = blk; iBlk = i; } }
        continue;
      }
      const lam = lv / VL, rhoC = oil ? pr.rhoO : p.mode === 'water' ? pr.rhoW : pr.rhoL, muC = oil ? pr.muO : p.mode === 'water' ? pr.muW : pr.muL, nuC = muC / rhoC, kC = oil ? pr.kO : p.mode === 'water' ? pr.kW : pr.kL;
      const rhoM = pr.rhoL * HL + pr.rhoG * (1 - HL), epsT = (4 * tauW[i] * vm) / (rhoM * D), gdot = Math.max(Math.sqrt(epsT / nuC), (8 * vL) / D, 0.05);
      const Peq = S.peq(Ti + dep), Dg = (7.4e-12 * Math.sqrt(oil ? pr.mwO || 150 : 46.8) * TK) / (muC * 1000 * 37.7 ** 0.6); // methane in the continuous liquid (Wilke–Chang)
      let cpR = 1e9, dVol = 0, limTxt = dTs > 0 ? (lam >= p.lamStar ? 'no free water' : 'induction (no nuclei yet)') : 'outside the hydrate region', depVol = 0, mFilmWall = 0;
      if (!flowing) {
        // shut-in: phases are segregated, hydrate grows as a film at the water interface, fed by diffusion
        const Ai = D * Math.sqrt(Math.max(1 - (2 * HL - 1) ** 2, 0.05)) * ds;
        if (lam >= p.lamStar && dTs > 0 && mWfree > 0) {
          if (!onset) onset = { t: t + dt, x: S.x[i], i };
          const dc = (fug(P[i], zG) - fugAt(Peq, P[i], zG)) / p.H, d0 = film[i] / Ai, d1 = dc > 0 ? filmAdvance(d0, dc * p.vmh * dt, 1 / (kStarAt(TK) * p.H) + (0.25 * D * HL + 1e-3) / Dg, p.shellD) : d0; // film thickness advanced exactly through its own resistance
          dVol = Math.min((d1 - d0) * Ai, (0.98 * mWfree) / (RHO * p.wfH), (0.9 * dTs * den * dt) / (RHO * LAT)); limTxt = 'diffusion through the interface film (shut-in)';
        } else if (dTs < 0 && film[i] > 0) { dVol = -Math.min(film[i], hydrateDissociationRate({ TK, P: P[i], Peq, zG, K0: p.disK0, E: p.disE }) * Ai * p.vmh * dt, (0.9 * -dTs * den * dt) / (RHO * LAT)); limTxt = 'dissociating'; }
        film[i] += dVol;
      } else {
        if (lam >= p.lamStar && dTs > 0 && mWfree > 0.02 * mW0) { // onset: the water phase is seeded as primary particles
          const target = (mWfree / pr.rhoW + hydVol / eexp) / vp;
          if (ps < target * 0.999) {
            const add = target - ps;
            if (oil) { for (let c = 0; c < nC; c++) if (seedW[c] > 0) Nc[c] += (add * vp * seedW[c]) / g.v[c]; pVol += add * vp; hydVol += add * vp * eexp * 1e-4; } else { Nc[0] += add; hydVol += add * g.v[0]; pVol = hydVol; }
            ps = target; if (!onset) onset = { t: t + dt, x: S.x[i], i };
          }
        }
        let any = false, dV = 0;
        if (hydVol > 0 && dTs > 0 && mWfree > 0) {
          const dc = (fug(P[i], zG) - fugAt(Peq, P[i], zG)) / p.H, kKin = kStarAt(TK) * p.H;
          let rPart = 0, lim = 'intrinsic kinetics';
          if (dc > 0) {
            if (oil) {
              const X = clamp(hydVol / Math.max(ps * vp * eexp, 1e-300), 0, 1), Sh = 2 + 0.6 * Math.sqrt((gdot * p.dPrim * p.dPrim) / nuC) * Math.cbrt(nuC / Dg);
              const gr = hydrateGrowthRate({ TK, P: P[i], Peq, zG, kRef: p.kinK, EaR: p.EaR, H: p.H, kFilm: (Sh * Dg) / p.dPrim, kShell: shellConductance(X, p.dPrim / 2, p.shellD), hPart: (2 * kC) / p.dPrim, dT: dTs, dHmol: p.dHmol });
              // conversion advanced exactly over the step (the shell resistance grows inside it); rPart is the step-mean rate
              const Ra = 1 / (kStarAt(TK) * p.H) + p.dPrim / (Sh * Dg) + (Number.isFinite(gr.jHeat) && gr.jHeat > 0 ? dc / gr.jHeat : 0), X1 = X < 1 ? shrinkingCoreAdvance(X, (dc * h * 6 * p.vmh) / (p.dPrim * eexp), Ra, p.dPrim / 2, p.shellD) : 1;
              rPart = ((X1 - X) * ps * vp * eexp) / (p.vmh * h); lim = gr.limiting;
            } else {
              for (let c = 0; c < nC; c++) { const L = g.L[c], j = 1 / (1 / (dc / (1 / kKin + L / (2 * Dg))) + p.dHmol / (((2 * kC) / L) * dTs)); gv[c] = PI * L * L * j; rPart += Nc[c] * gv[c]; }
              lim = 'surface growth (intrinsic + film)';
            }
            let s = 1;
            if (p.mode !== 'gas' && rPart > 0) { // gas must first dissolve in the liquid: absorption in series (Skovborg–Rasmussen)
              const kL = 0.4 * (Math.max(epsT, 1e-9) * nuC) ** 0.25 * Math.sqrt(Dg / nuC), aGL = ((D * Math.sqrt(Math.max(1 - (2 * HL - 1) ** 2, 0.05))) / (A * HL)) * (1 + (vm * vm) / (G * D)), rSup = (slugK[i] ? slugK[i].kLa * Math.sqrt(Dg / slugK[i].Dg) : kL * aGL) * VL * dc * p.mtMult; // slug flow: unit-cell conductance (dispersed bubbles in the slug body, stratified film behind it)
              s = 1 / (1 + rPart / Math.max(rSup, 1e-300)); if (s < 0.5) lim = 'gas absorption (mass transfer)';
            }
            const req = rPart * s * p.vmh * h, capW = (0.98 * mWfree) / (RHO * p.wfH), capG = Math.max(((0.9 * fe * pr.mG - hydIn * (1 - p.wfH)) * dt) / (q * RHO * (1 - p.wfH)), 0), capH = ((p.htMult ?? 1) * 0.9 * dTs * den * dt) / (RHO * LAT * q);
            const cap = Math.min(capW, capG, capH), sc = req > cap ? cap / req : 1;
            if (sc < 1) lim = cap === capH ? 'heat removal (fluid held at the hydrate temperature)' : cap === capW ? 'water availability' : 'gas availability';
            if (oil) dV = req * sc; else { const f = s * sc * p.vmh; for (let c = 0; c < nC; c++) gv[c] *= f; any = req > 0; }
            limTxt = lim;
          }
        } else if (hydVol > 0 && dTs < 0) {
          const jd = hydrateDissociationRate({ TK, P: P[i], Peq, zG, K0: p.disK0, E: p.disE }) * p.vmh, capH = (0.9 * -dTs * den * dt) / (RHO * LAT * q);
          if (oil) dV = -Math.min(hydVol, jd * ps * PI * p.dPrim * p.dPrim * h, capH);
          else { let tot = 0; for (let c = 0; c < nC; c++) { gv[c] = -PI * g.L[c] ** 2 * jd; tot -= Nc[c] * gv[c]; } const sc = tot * h > capH ? capH / (tot * h) : 1; for (let c = 0; c < nC; c++) gv[c] *= sc; any = tot > 0; }
          limTxt = 'dissociating';
        }
        if (pVol > 0) {
          const phiNow = hydVol / VL, cc = cpC[i]; let cp = cc && Math.abs(phiNow / cc.phi - 1) < 0.03 && Math.abs(gdot / cc.gdot - 1) < 0.03 && Math.abs(muC / cc.mu - 1) < 0.03 ? cc.cp : null;
          if (!cp) { cp = maxAgglomerateSize({ dp: p.dPrim, Fa: p.cohesion * p.dPrim, mu0: muC, shear: gdot, phi: phiNow, phiMax: p.phiMax, fr: p.fr }); cpC[i] = { phi: phiNow || 1e-300, gdot, mu: muC, cp }; }
          const agg = p.aggEff > 0 && dTs > 0;
          rec.dA[i] = cp.dA; cpR = cp.ratio;
          if (agg) {
            const cS = (p.aggEff * ((8 * vL) / D / 6 + 0.1618 * Math.sqrt(epsT / nuC))) / VL, cDf = (p.aggEff * (PI / 4) * Math.abs(RHO - rhoC) * G) / (18 * muC) / VL, cB = (p.aggEff * 2 * KB * TK) / (3 * muC) / VL;
            beta[0] = cS; beta[1] = cDf; beta[2] = cB;
          }
          { const kb = (p.kBreak * gdot) / cp.dA ** 3; for (let c = 0; c < nC; c++) Sb[c] = g.L[c] > 1.01 * p.dPrim ? kb * L3[c] : 0; } // agglomerates above the cohesive limit break even when nothing sticks any more
          const r = solvePBE(g, Nc, h, { beta: agg ? beta : null, gv: any ? gv : null, S: Sb, maxSub: p.maxSub, frac: 0.3 });
          if (r.limited) limited = true;
          pVol = 0; for (let c = 0; c < nC; c++) { Nc[c] = r.N[c]; pVol += Nc[c] * g.v[c]; }
        }
        hydVol = oil ? hydVol + dV : pVol;
        dVol = hydVol - volBefore;
        if (!oil && dTs < 0 && volBefore > 0) ps *= clamp(hydVol / volBefore, 0, 1);
        if (hydVol < 1e-30 || (dTs < 0 && hydVol < 1e-9 * pVol)) { Nc.fill(0); ps = 0; hydVol = 0; pVol = 0; dVol = -volBefore; }
        // ---- wall: particle capture on a sub-cooled wall, hydrate film from water vapour on the gas-wetted wall
        const hIn = hInside((pr.rhoL * vm * D) / pr.muL, (pr.cpL * pr.muL) / pr.kL, pr.kL, D), Tw = Ti - (S.U * (Ti - S.tAmb[i])) / Math.max(hIn, 1), cold = Tw < Teq, wet = /strat/.test(ref.regime[i]) ? Math.acos(clamp(1 - 2 * HL, -1, 1)) / PI : 1;
        rec.tw[i] = Tw;
        if (cold && pVol > 0 && p.adhesion > 0) {
          const uS = Math.sqrt(tauW[i] / pr.rhoL), geo = (PI * D * wet) / (A * HL); let lost = 0;
          for (let c = 0; c < nC; c++) {
            if (!(Nc[c] > 0)) continue;
            const L = g.L[c], lamK = depositionVelocity(L, RHO, uS, nuC, rhoC, TK) * p.adhesion * Math.min(1, p.adhForce / (8 * Math.max(tauW[i], 1e-6) * L)) * geo, loss = Nc[c] * (1 - Math.exp(-lamK * h));
            Nc[c] -= loss; lost += loss * g.v[c];
          }
          if (lost > 0) { const fP = lost / pVol; depVol = oil ? hydVol * fP : lost; ps *= 1 - fP; wd += depVol * RHO * p.wfH; hydVol -= depVol; pVol -= lost; }
        }
        if (cold && wet < 1 && vsg > 0 && mWfree > 0 && p.filmMult > 0) {
          const conv = pr.rhoG / ((pr.mwG * 1e-3) / VM_STD), dC = (gasSat * waterContent(P[i], Ti) - waterContent(P[i], Math.max(Tw, -20))) * conv, Dwg = 2.2e-5 * (1.013 / P[i]) * (TK / 273.15) ** 1.75, Sh = 0.023 * ((pr.rhoG * vsg * D) / pr.muG) ** 0.83 * (pr.muG / (pr.rhoG * Dwg)) ** 0.33;
          mFilmWall = Math.min((p.filmMult * Math.max(dC, 0) * ((Sh * Dwg) / D) * PI * D * (1 - wet) * ds * dt) / p.wfH, (0.5 * mWfree * q) / p.wfH); wd += (mFilmWall * p.wfH) / q;
        }
      }
      // ---- ledger, heat of formation, deposit
      const mNew = dVol * RHO * (flowing ? q : 1);
      if (mNew >= 0) { led.formed += mNew; if (slugK[i]) led.formedSlug += mNew; } else led.dissociated -= mNew;
      led.formed += mFilmWall; led.filmWall += mFilmWall; led.heat += mNew * LAT;
      Ti += (mNew * LAT) / dt / den;
      let m = mDep[i]; rec.removal[i] = 0;
      const tauC = p.tauCrit * ((1 - por[i]) / (1 - p.por0)) ** 2;
      if (m > 0 && flowing && tauW[i] > tauC && p.kRemove > 0) { // shear removal / sloughing back into the stream
        const dm = m * (1 - Math.exp(-p.kRemove * (tauW[i] / tauC - 1) * dt)), c = kSl(clamp(dHyd, p.dPrim, g.L[nC - 1])), fv = dm / (q * RHO);
        rec.removal[i] = dm / dt; m -= dm; Nc[c] += fv / g.v[c]; pVol += fv; ps += fv / (vp * eexp); wd -= (dm * p.wfH) / q; hydVol += fv; led.sloughed += dm;
      }
      if (m > 0 && dTs < 0) { // a deposit outside the hydrate region melts at the Kim–Bishnoi rate on its exposed surface, as fast as heat arrives
        const dm = Math.min(m, hydrateDissociationRate({ TK, P: P[i], Peq, zG, K0: p.disK0, E: p.disE }) * p.vmh * RHO * PI * D * ds * dt, (0.9 * Math.max(Ti - Teq, 0) * den * dt) / LAT);
        if (dm > 0) { m -= dm; led.dissociated += dm; led.heat -= dm * LAT; Ti -= (dm * LAT) / dt / den; if (flowing) wd -= (dm * p.wfH) / q; }
      }
      // equipment boundary (choke, low point, strainer): a fraction of the solids leaving the cell is retained there
      const eta = flowing ? trapEff[i] : 0, mTrap = eta > 0 ? eta * a * hydVol * RHO : 0;
      if (mTrap > 0) { led.trapped += mTrap; wd += eta * hydVol * RHO * p.wfH * (a / q); }
      const gain = depVol * RHO * q + mFilmWall + mTrap;
      if (gain > 0) { const Vo = m / (RHO * (1 - por[i])), Vn = gain / (RHO * (1 - p.por0)); por[i] = (Vo * por[i] + Vn * p.por0) / (Vo + Vn); m += gain; led.captured += depVol * RHO * q + mTrap; }
      por[i] = p.porInf + (por[i] - p.porInf) * Math.exp(-dt / p.tAge);
      rec.depRate[i] = (m - mDep[i]) / dt; mDep[i] = m;
      // ---- slurry state for the hydraulics of the next step
      const phiH = hydVol / VL; let m3 = 0, m4 = 0; for (let c = 0; c < nC; c++) { m3 += Nc[c] * L3[c]; m4 += Nc[c] * L4[c]; }
      const d43 = m3 > 0 ? m4 / m3 : 0, ratio = Math.max(1, Math.min(d43 / p.dPrim, 1.5 * cpR)), // the cohesive limit caps the size that enters the rheology
        phiE = Math.min(phiH * ratio ** (3 - p.fr), 1);
      muFac[i] = slurryViscosity(phiE, p.viscModel, { phiMax: p.phiMax });
      HV[i] = hydVol; PS[i] = ps; WD[i] = wd; LV[i] = lv; T[i] = Ti;
      if (flowing) { const pass = 1 - eta; for (let c = 0; c < nC; c++) FinN[c] = (pass * Nc[c]) / tau; FinPS = (pass * ps) / tau; FinWD = wd / tau; FinLV = lv / tau; FinHV = (pass * hydVol) / tau; hydIn = (pass * hydVol * RHO) / tau; thUp = Ti - ref.T[i]; }
      const blk = 1 - (D * D) / (D0 * D0), rate = mNew / dt + mFilmWall / dt;
      rec.gdot[i] = gdot; rec.muC[i] = muC; rec.rhoC[i] = rhoC; rec.VL[i] = VL; if (rec.removal[i] > rec.removalMax) rec.removalMax = rec.removal[i]; if ((depVol * RHO * q + mFilmWall) / dt > rec.captureMax) rec.captureMax = (depVol * RHO * q + mFilmWall) / dt; if (dTs > 0) rec.exposure[i] += dt / 3600; rec.capture[i] = (depVol * RHO * q + mFilmWall) / dt;
      rec.sub[i] = dTs; rec.teq[i] = Teq; rec.phi[i] = phiH + film[i] / VL; rec.phiE[i] = phiE; rec.d43[i] = d43; rec.J[i] = Jn; if (Jn > rec.Jmax) rec.Jmax = Jn; rec.lam[i] = lam; rec.rate[i] = rate; rec.vL[i] = vL; rec.lim[i] = limTxt; rec.freeW[i] = mW0 > 0 ? mWfree / mW0 : 0;
      rec.X[i] = mW0 > 0 ? clamp(((hydVol + film[i]) * RHO * p.wfH) / mW0, 0, 1) : 0;
      if (dTs > rec.subMax[i]) rec.subMax[i] = dTs; if (rec.phi[i] > rec.phiMaxT[i]) rec.phiMaxT[i] = rec.phi[i];
      sumRate += rate; sumSusp += (hydVol + film[i]) * RHO;
      if (blk > maxBlk) { maxBlk = blk; iBlk = i; } if (blk >= p.plugBlock && mDep[i] > 0 && iPlugDep < 0) iPlugDep = i; if (rec.phi[i] > maxPhi) { maxPhi = rec.phi[i]; iPhi = i; } if (muFac[i] > maxVisc) { maxVisc = muFac[i]; iVisc = i; } if (dTs > maxSub) maxSub = dTs;
    }
    if (flowing) { led.exported += dt * hydIn; released = true; }
    t += dt; tPh += dt; k++;
    if (!plug && iPlugDep >= 0) plug = { t, x: S.x[iPlugDep], i: iPlugDep, mech: 'hydrate wall deposit closes the bore' };
    if (!plug && flowing && maxVisc >= 1000) plug = { t, x: S.x[iVisc], i: iVisc, mech: 'slurry jams (effective solids fraction at the packing limit)' };
    let depTot = 0; for (let i = 0; i < n; i++) depTot += mDep[i];
    if (maxPhi > peak.phi) Object.assign(peak, { phi: maxPhi, tPhi: t, iPhi, N: Array.from(NP[iPhi]), phiX: rec.phi.slice(), phiE: rec.phiE.slice(), d43: rec.d43.slice(), dA: rec.dA.slice(), gdot: rec.gdot[iPhi], muC: rec.muC[iPhi], rhoC: rec.rhoC[iPhi], VL: rec.VL[iPhi], vL: rec.vL[iPhi], T: T[iPhi], visc: muFac.slice() });
    if (maxBlk > peak.blk) Object.assign(peak, { blk: maxBlk, tBlk: t, iBlk, dHyd: Dn.map((d, i) => Math.max((Dbase[i] - d) / 2, 0)), por: Array.from(por) });
    ser.exp.push(led.exported); ser.t.push(t / 3600); ser.pIn.push(pInNow); ser.dp.push(pInNow - P[n - 1]); ser.susp.push(sumSusp); ser.dep.push(depTot); ser.blk.push(maxBlk); ser.phi.push(maxPhi); ser.rate.push(sumRate); ser.sub.push(maxSub); ser.frac.push(flowing ? fe : 0); ser.visc.push(maxVisc);
    if (k % stride === 0 || plug || tPh >= phz.dur - 1e-9) { fld.t.push(t / 3600); fld.phi.push(rec.phi.slice()); fld.dep.push(Dn.map((d, i) => (Dbase[i] - d) / 2)); fld.sub.push(rec.sub.slice()); }
    if (tPh >= phz.dur - 1e-9) { ph++; tPh = 0; newPhase = true; }
    if (plug || ph >= o.phases.length) done = true;
    if (o.stopWhenClear && flowing && onset && maxSub < 0 && sumSusp <= 0 && depTot <= 0) done = true; // nothing left that could still plug
  }
  return {
    step, get done() { return done; }, get progress() { return total > 0 ? t / total : 1; },
    result() {
      let susp = 0, depTot = 0; for (let i = 0; i < n; i++) { let hv = film[i] + HV[i]; susp += hv * RHO; depTot += mDep[i]; }
      const Dend = Dbase.map((d, i) => Math.sqrt(Math.max(d * d - (4 * mDep[i]) / (PI * RHO * (1 - por[i]) * ds), (0.03 * D0) ** 2))), dHyd = Dend.map((d, i) => Math.max((Dbase[i] - d) / 2, 0));
      return { t: t / 3600, T: T.slice(), P: P.slice(), D: Dend, dHyd, mDep: Array.from(mDep), por: Array.from(por), hold: hold.slice(), tauW: tauW.slice(), dpdx: dpdx.slice(), muFac: muFac.slice(), N: NP.map((a) => Array.from(a)), rec, ser, fld, peak, onset, plug, limited, pIn: pInNow,
        ledger: { ...led, suspended: susp, deposited: depTot, deposit0: mDep0, in: led.formed + mDep0 + led.inflow + led.susp0, out: susp + depTot + led.exported + led.dissociated } };
    },
  };
}
/** Run a hydrate march to completion (synchronous helper for the Monte Carlo loop, the verification cases and scripts). */
export function runHydrateMarch(S, p, o) { const m = hydrateMarch(S, p, o); while (!m.done) m.step(); return m.result(); }
/**
 * Implicit upwind transport of cell inventories E through cells of residence time tau (s) over a step dt with an inlet
 * flux fin (per second): E_i ← (E_i + dt·F_in)/(1 + dt/τ_i), F_out = E_i/τ_i. Returns the outlet flux. Conservative and
 * unconditionally stable; this is the transport operator of hydrateMarch.
 */
export function advectImplicit(E, tau, dt, fin = 0) { let F = fin; for (let i = 0; i < E.length; i++) { E[i] = (E[i] + dt * F) / (1 + dt / tau[i]); F = E[i] / tau[i]; } return F; }

// =====================================================================================================
// 13. Slow deposits in steady production: wax and mineral scale along the line
// =====================================================================================================
const WATER0 = Object.freeze({ Na: 12500, K: 250, Ca: 900, Mg: 150, Ba: 40, Sr: 90, Fe: 5, Cl: 21500, SO4: 15, HCO3: 650 });
/**
 * Boundary and initial data of the hydrate march for a set-up S: inhibitor per cell (mid-line injection points), heated
 * section, depressurisation of the shut-in line, equipment that retains solids, free-water distribution at rest and the
 * initial suspended hydrate population. Returns the option fields understood by hydrateMarch.
 */
export function marchExtras(S, v) {
  const n = S.n, heat = new Array(n).fill(0), a = Math.min(v.heatFromKm, v.heatToKm) * 1000, b = Math.max(v.heatFromKm, v.heatToKm) * 1000; let heated = 0;
  if (v.heatWm > 0) for (let i = 0; i < n; i++) if (S.s[i] >= a && S.s[i] <= b) { heat[i] = v.heatWm; heated += S.ds; }
  const traps = (Array.isArray(v.equip) ? v.equip : []).map((r) => ({ x: num(r?.x, 0) * 1000, eff: clamp(num(r?.eff, 0) / 100, 0, 1), name: String(r?.name || 'equipment') })).filter((r) => r.eff > 0).map((r) => ({ ...r, i: clamp(Math.floor(r.x / S.ds), 0, n - 1) }));
  let wMult = null;
  if (v.waterSettle === 'lowpoints') { // free water drains down-slope while the line is at rest and collects at the low points (up to a full bore)
    const b0 = S.prof(1), A = (PI * S.D0 * S.D0) / 4, V0 = b0.P.map((P, i) => A * b0.holdup[i] * S.props(P, b0.T[i]).wcut * S.ds), V = V0.slice(), cap = A * S.ds;
    for (let pass = 0; pass < 4 * n; pass++) { let moved = 0; for (let i = 0; i < n; i++) { if (!(V[i] > 0)) continue; const zl = i > 0 ? S.z[i - 1] : Infinity, zr = i < n - 1 ? S.z[i + 1] : Infinity, j = zl < zr ? i - 1 : i + 1; if (!(Math.min(zl, zr) < S.z[i] - 1e-3)) continue; const mv = Math.min(V[i], cap - V[j]); if (mv > 1e-12) { V[i] -= mv; V[j] += mv; moved += mv; } } if (moved < 1e-9) break; }
    wMult = V.map((x, i) => (V0[i] > 0 ? x / V0[i] : 1));
  }
  return { inhWt: S.aq.inhWtCell, heat: v.heatWm > 0 ? heat : null, heatStart: v.heatStartH * 3600, heatedLength: heated, depress: v.depressAtH > 0 ? { t: v.depressAtH * 3600, P: v.depressP } : null, traps, wMult,
    init: v.initHydPct > 0 ? { phi: v.initHydPct / 100, d: v.initAggUm * 1e-6, sig: v.initAggSigma } : null };
}
const waterOf = (v) => { const row = Array.isArray(v.water) && v.water[0] ? v.water[0] : WATER0; return Object.fromEntries(ION_IDS.map((id) => [id, Math.max(num(row[id], 0), 0)])); };
/**
 * Wax and scale layers after `days` of steady production at a rate fraction.
 * Wax: Fick diffusion on the solubility slope and the wall heat flux, shear dispersion, Brownian flux, shear stripping,
 * ageing and the insulating feedback of the layer on the heat flux and on the bulk temperature downstream.
 * Scale: saturation indices at the local P, T, surface-reaction kinetics in series with ion mass transfer to the wall.
 * Returns { dWax[], Fw[], dScale[], si[] (max SI), mineral[], scaleRate[] (mm/y), ser: { t (d), dMax (mm), mass (kg) }, waxMass, waxRate0 (mm/d), piggingInterval (d) | null, ... }.
 */
export function slowDeposits(S, v, frac = 1, nt = 40) {
  const wc = v.waxThermo === 'sle' && S.waxCurve && S.waxCurve.ok && S.waxCurve.wat !== null && S.waxCurve.wTot > 0 ? S.waxCurve : null, wat = wc ? wc.wat : v.wat; // computed solid–liquid equilibrium or the entered WAT and slope
  const b = S.prof(frac), n = S.n, D0 = S.D0, ds = S.ds, days = v.depositDays, dt = (days * 86400) / nt, wTot = wc ? wc.wTot : v.waxContent / 100, F0 = clamp(1 - (v.waxOil > 0 ? v.waxOil : 80) / 100, 0.03, 0.95), inhW = v.waxInhOn ? 1 - v.waxInhEff / 100 : 1, gelC = v.gelCoef * (v.waxInhOn ? 1 - v.waxInhGel / 100 : 1);
  const dWax = new Array(n).fill(0), Fw = new Array(n).fill(F0), sea = v.srpOn ? { ...SEAWATER, SO4: Math.min(v.srpSO4, SEAWATER.SO4) } : SEAWATER, water = mixWaters(waterOf(v), sea, clamp(v.swFrac / 100, 0, 1));
  const cells = b.P.map((P, i) => { const pr = S.props(P, b.T[i]), A = (PI * D0 * D0) / 4, vsl = (pr.qL * frac) / A, vsg = (pr.qG * frac) / A, vm = vsl + vsg, H = b.holdup[i];
    return { pr, vL: vsl / H, vm, hIn: hInside((pr.rhoL * vm * D0) / pr.muL, (pr.cpL * pr.muL) / pr.kL, pr.kL, D0), mcp: frac * (pr.mG * pr.cpG + pr.mO * pr.cpO + pr.mW * pr.cpW), gammaW: b.tauW[i] / pr.muL, oilWet: pr.phaseInv ? 0.15 : 1 - 0.5 * pr.wcut, waterWet: pr.phaseInv ? 1 : clamp(pr.wcut, 0, 1) }; });
  const ser = { t: [0], dMax: [0], mass: [0] }, last = new Array(n).fill(null), Tb = b.T.slice();
  let waxRate0 = 0;
  for (let k = 0; k < nt; k++) {
    let cum = 0, dMax = 0, mass = 0, used = 0;
    for (let i = 0; i < n; i++) {
      const c = cells[i], left = clamp(1 - used / Math.max(c.pr.mO * frac * wTot, 1e-12), 0, 1); Tb[i] = S.tAmb[i] + (b.T[i] - S.tAmb[i]) * Math.exp(Math.min(cum, 3));
      const r = waxDeposition({ Tb: Tb[i], Tamb: S.tAmb[i], U: S.U, hIn: c.hIn, kOil: c.pr.kO, rhoOil: c.pr.rhoO, muOil: c.pr.muO, wat, wTot: wTot * left, slope: v.waxSlope, delta: dWax[i], Fw: Fw[i], kDep: v.waxK, D: D0 - 2 * dWax[i], vL: c.vL, gammaW: c.gammaW, rhoMix: b.rhoM[i], regime: b.regime[i], mult: v.waxMult * inhW, diffModel: v.waxDiff, wetFrac: c.oilWet, curve: wc, stripC: v.waxStripC, aspect: v.waxAspect, wmodel: v.waxModel, c1: v.waxC1, coil: v.waxOil > 0 ? v.waxOil / 100 : null }); if (r.Fset !== undefined) Fw[i] = r.Fset;
      used += 900 * (r.dDelta * Fw[i] + r.dFw * dWax[i]) * PI * D0 * ds; // wax leaving the oil: the dissolved wax available downstream is depleted
      if (k === 0 && r.dDelta * 86400e3 > waxRate0) waxRate0 = r.dDelta * 86400e3;
      cum += ((S.U - r.Ueff) * PI * D0 * ds) / Math.max(c.mcp, 1e-9); last[i] = r;
      dWax[i] = Math.min(dWax[i] + r.dDelta * dt, 0.4 * D0); Fw[i] = Math.min(Fw[i] + r.dFw * dt, 0.95);
      if (dWax[i] > dMax) dMax = dWax[i]; mass += 900 * Fw[i] * (PI / 4) * (D0 * D0 - (D0 - 2 * dWax[i]) ** 2) * ds;
    }
    ser.t.push(((k + 1) * dt) / 86400); ser.dMax.push(dMax * 1000); ser.mass.push(mass);
  }
  const lim = v.waxLimitMm, end = ser.dMax[nt]; let pig = null;
  if (end >= lim) { for (let k = 1; k <= nt; k++) if (ser.dMax[k] >= lim) { pig = ser.t[k - 1] + ((lim - ser.dMax[k - 1]) / (ser.dMax[k] - ser.dMax[k - 1] || 1)) * (ser.t[k] - ser.t[k - 1]); break; } }
  else if (end > 1e-6) { const rate = Math.max((end - ser.dMax[Math.floor(nt / 2)]) / (ser.t[nt] - ser.t[Math.floor(nt / 2)]), 1e-9); pig = Math.min(days + (lim - end) / rate, 3650); }
  // gel strength after a cold shutdown: yield stress from the wax precipitated at ambient temperature
  let restartDp = 0, gelLen = 0, tauY = 0;
  for (let i = 0; i < n; i++) { const sol = waxSolubility(S.tAmb[i], wat, wTot, v.waxSlope, wc).solid * 100, ty = sol > 0 && S.tAmb[i] < wat - v.pourOffset ? gelC * sol * sol : 0; if (ty > 0) { restartDp += (4 * ty * ds) / D0; gelLen += ds; tauY = Math.max(tauY, ty); } }
  // scale
  const si = [], mineral = [], scaleRate = [], dScale = [], sis = [], ctlAll = [], opt = { yCO2: v.co2Pct / 100, model: v.actModel };
  const nSt = Math.min(n, 13), stI = Array.from({ length: nSt }, (_, k) => Math.round((k * (n - 1)) / Math.max(nSt - 1, 1))), stR = stI.map((i) => scaleIndices(water, b.T[i], b.P[i], opt));
  const siAt = (i) => { // full speciation at up to 13 stations; saturation indices, precipitation potentials and pH interpolated linearly in between
    let k = 0; while (k < nSt - 2 && stI[k + 1] < i) k++; const i0 = stI[k], i1 = stI[Math.min(k + 1, nSt - 1)], f = i1 > i0 ? clamp((i - i0) / (i1 - i0), 0, 1) : 0, A = stR[k], Bq = stR[Math.min(k + 1, nSt - 1)]; if (f === 0) return A; if (f === 1) return Bq;
    const minerals = A.minerals.map((m, q) => ({ ...m, SI: m.SI + f * (Bq.minerals[q].SI - m.SI), ptb: m.ptb + f * (Bq.minerals[q].ptb - m.ptb), logK: m.logK + f * (Bq.minerals[q].logK - m.logK) })), max = minerals.reduce((x, y) => (y.SI > x.SI ? y : x));
    return { ...A, pH: A.pH + f * (Bq.pH - A.pH), I: A.I + f * (Bq.I - A.I), aw: A.aw + f * (Bq.aw - A.aw), fCO2: A.fCO2 + f * (Bq.fCO2 - A.fCO2), minerals, max: { id: max.id, name: max.name, SI: max.SI } }; };
  for (let i = 0; i < n; i++) {
    const r = siAt(i), c = cells[i], T = b.T[i] + KEL, ReW = (c.pr.rhoL * c.vm * D0) / c.pr.muL, Sc = c.pr.muW / (c.pr.rhoW * 1e-9), km = (0.023 * ReW ** 0.83 * Sc ** 0.33 * 1e-9) / D0;
    let rate = 0; const ctl = {};
    for (const m of r.minerals) { if (!(m.SI > 0) || !(c.pr.mW > 0)) continue;
      // nucleation kinetics: nothing grows on the wall before the induction time has passed; a threshold inhibitor lengthens it and blocks growth sites
      const cat = MINERALS.find((q) => q.id === m.id).cat, q = scaleInhibition(m.id, { SI: m.SI, TK: T, pH: r.pH, R: r.tot.SO4 > 0 ? r.tot[cat] / r.tot.SO4 : 1, Ca: r.tot.Ca * r.kgw, inh: v.scaleInh, dose: v.scaleInhOn ? v.scaleInhDose : 0, tProtect: v.scaleProtectH * 3600 }), gate = q.model ? clamp(1 - q.tInh / (days * 86400), 0, 1) * (1 - q.eff) : 1;
      ctl[m.id] = q; if (!(gate > 0)) continue;
      const kin = gate * v.scaleK * 1e-8 * Math.exp((-30800 / R) * (1 / T - 1 / 298.15)) * (10 ** (m.SI / 2) - 1) ** 2, mt = km * (m.ptb / m.M), fl = 1 / (1 / Math.max(kin, 1e-300) + 1 / Math.max(mt, 1e-300)); rate += (c.waterWet * fl * m.M * 1e-3) / (m.rho * 0.8); } // m/s of a 20 % porous layer
    ctlAll.push(ctl); si.push(r.max.SI); mineral.push(r.max.SI > 0 ? r.max.name : 'none'); scaleRate.push(rate * 3.156e10); dScale.push(Math.min(rate * days * 86400, 0.2 * D0)); sis.push(r);
  }
  // population-balance crystallisation in the bulk water (plug flow along the line): nucleation and growth moments of the minerals that are supersaturated
  const cryst = [], A0 = (PI * D0 * D0) / 4;
  for (const mn of MINERALS) {
    const k = MINERALS.indexOf(mn); if (!sis.some((r) => r.minerals[k].SI > 0.02)) continue;
    const nSeed = v.scaleSeedLog > 0 ? 10 ** v.scaleSeedLog : 0, Ls = 1e-6, mu3s = nSeed * Ls ** 3; let st = [nSeed, nSeed * Ls, nSeed * Ls * Ls, mu3s, 0], peakJ = 0; const carb = mn.an === 'CO3', L10 = [], num = []; // suspended fines act as seed crystals
    for (let i = 0; i < n; i++) {
      const c = cells[i], r = sis[i], SI = r.minerals[k].SI, T = b.T[i] + KEL, tau = (A0 * b.holdup[i] * ds) / Math.max(c.pr.qL * frac, 1e-12), conv = r.kgw * 1000, cA = Math.max(r.tot[mn.cat] * conv, 1e-9), cB = Math.max((carb ? r.tot.HCO3 : r.tot.SO4) * conv, 1e-9);
      const kq = ctlAll[i][mn.id]; if (kq && kq.model && kq.tInh > v.scaleProtectH * 3600) { L10.push(st[0] > 0 ? st[1] / st[0] : 0); num.push(st[0]); continue; } // inside the induction time: no bulk crystallisation in this cell
      const q = scaleCrystallisation({ S0: 10 ** (Math.max(SI, -3) / 2), cA, cB, nuB: carb ? 2 : 1, M: mn.M / 1000, rho: mn.rho, TK: T, t: tau, sigma: v.scaleSigma * 1e-3, A: 10 ** v.scaleNucA, kg: v.scaleKg * 1e-10 * Math.exp((-30800 / R) * (1 / T - 1 / 298.15)), init: st, steps: 12 });
      st = q.state; if (q.J0 > peakJ) peakJ = q.J0; L10.push(q.d10); num.push(q.number);
      if (q.massConc > 0 && c.pr.mW > 0) { // crystals carried to the wall by turbulence and held there
        const uS = Math.sqrt(b.tauW[i] / c.pr.rhoL), add = (c.waterWet * v.scaleDepEff * depositionVelocity(Math.max(q.d32, 1e-8), mn.rho, uS, c.pr.muW / c.pr.rhoW, c.pr.rhoW, T) * q.massConc) / (mn.rho * 0.8);
        scaleRate[i] += add * 3.156e10; dScale[i] = Math.min(dScale[i] + add * days * 86400, 0.2 * D0);
      }
    }
    cryst.push({ id: mn.id, name: mn.name, d10: st[0] > 0 ? st[1] / st[0] : 0, d32: st[2] > 0 ? st[3] / st[2] : 0, number: st[0], x: st[4], massMgL: (PI / 6) * mn.rho * (st[3] - mu3s) * 1000, molCheck: ((PI / 6) * mn.rho * (st[3] - mu3s)) / (mn.M / 1000), peakJ, L10, num });
  }
  const control = MINERALS.map((mn, k) => { let iM = 0; for (let i = 1; i < n; i++) if (sis[i].minerals[k].SI > sis[iM].minerals[k].SI) iM = i; let worst = null, iW = iM; for (let i = 0; i < n; i++) { const q = ctlAll[i][mn.id]; if (q && (!worst || q.risk > worst.risk)) { worst = q; iW = i; } }
    const SI = sis[iW].minerals[k].SI, o = { TK: b.T[iW] + KEL, pH: sis[iW].pH, R: sis[iW].tot.SO4 > 0 ? sis[iW].tot[mn.cat] / sis[iW].tot.SO4 : 1, Ca: sis[iW].tot.Ca * sis[iW].kgw, inh: v.scaleInh, tProtect: v.scaleProtectH * 3600 };
    return { id: mn.id, name: mn.name, SImax: sis[iM].minerals[k].SI, i: iW, SI, q: worst, raw: worst ? scaleInhibition(mn.id, { ...o, SI, dose: 0 }) : null, siCrit0: worst && worst.model ? scaleCriticalSI(mn.id, { ...o, dose: 0 }) : null, siCrit: worst && worst.model ? scaleCriticalSI(mn.id, { ...o, dose: v.scaleInhOn ? v.scaleInhDose : 0 }) : null }; });
  return { control, sea, dWax, Fw, dScale, si, mineral, scaleRate, sis, ser, waxMass: ser.mass[nt], waxRate0, piggingInterval: pig, last, Tb, restartDp: restartDp / 1e5, gelLen, tauY, water, prof: b, cells, wat, waxTotal: wTot, waxCurve: wc, cryst };
}

// =====================================================================================================
// 14. Suite inputs
// =====================================================================================================
const SCEN = [{ value: 'restart', label: 'Cold restart after a shutdown' }, { value: 'steady', label: 'Steady production' }, { value: 'turndown', label: 'Turndown (reduced rate)' }, { value: 'shutdown', label: 'Shutdown cooldown (shut-in line)' }];
const flowing = (v) => v.scenario !== 'shutdown', shut = (v) => v.scenario === 'restart' || v.scenario === 'shutdown';
const INPUTS = [
  { group: 'Operating scenario', tab: 'inputs', help: 'In steady production the reference line runs far above the hydrate temperature, so hydrates are assessed for the off-design events (turndown, shutdown, restart) while wax, scale, asphaltene and sand are assessed for steady production over the period since the last pig run.', fields: [
    { key: 'scenario', label: 'Scenario', type: 'select', value: 'restart', options: SCEN, help: 'Cold restart: the line is shut in for the stated hours (cooling towards ambient at the settle-out pressure) and then restarted at a reduced rate.' },
    { key: 'shutHours', label: 'Shut-in duration', unit: 'h', value: 48, min: 0.5, max: 2000, typical: [8, 96], showIf: shut, help: 'Hours without flow before the restart, or the length of the cooldown that is simulated.' },
    { key: 'restartPct', label: 'Restart rate', unit: '% of case rate', value: 50, min: 5, max: 100, typical: [20, 60], showIf: (v) => v.scenario === 'restart' },
    { key: 'turndownPct', label: 'Turndown rate', unit: '% of case rate', value: 30, min: 5, max: 100, typical: [20, 70], showIf: (v) => v.scenario === 'turndown' },
    { key: 'simHours', label: 'Simulated flowing period', unit: 'h', value: 24, min: 0.5, max: 720, typical: [6, 72], showIf: flowing, help: 'Length of the hydrate event that is marched in time after the scenario starts.' },
    { key: 'depositDays', label: 'Production period since the last pig run', unit: 'd', value: 30, min: 1, max: 720, typical: [7, 180], help: 'Wax and scale layers are grown over this period of steady production and are present when the hydrate event starts.' },
    { key: 'pShut', label: 'Topsides pressure during shut-in', unit: 'bara', value: 0, min: 0, max: 600, showIf: shut, help: '0 = the line is held at the flowing arrival pressure. A lower value represents a partly depressurised line.' },
    { key: 'pInMax', label: 'Available inlet pressure', unit: 'bara', value: 250, min: 5, max: 1400, typical: [100, 400], help: 'Highest pressure the wells or the pump can supply at the flowline inlet. When a restriction needs more than this the rate falls and the line eventually stalls.' },
    { key: 'depressAtH', label: 'Depressurise the shut-in line after', unit: 'h', value: 0, min: 0, max: 2000, showIf: shut, help: '0 = no depressurisation. After this time the topsides pressure of the shut-in line is lowered to the depressurisation pressure entered with the plugging data (a depressurisation boundary inside the time march).' },
  ] },
  { group: 'Initial state of the line', tab: 'inputs', help: 'What is in the line when the simulated sequence starts, in addition to the deposits of the production period.', fields: [
    { key: 'initHydPct', label: 'Hydrate already suspended in the liquid', unit: 'vol % of liquid', value: 0, min: 0, max: 30, help: 'Initial hydrate concentration (saturation) in every cell.' },
    { key: 'initAggUm', label: 'Size of the initial particles or agglomerates', unit: 'µm', value: 0, min: 0, max: 20000, help: '0 = the primary particle size. A larger value starts the population balance from an agglomerate population.' },
    { key: 'initAggSigma', label: 'Geometric standard deviation of the initial population', unit: '–', value: 1.5, min: 1, max: 4 },
    { key: 'waterSettle', label: 'Free water while the line is at rest', type: 'select', value: 'uniform', options: [{ value: 'uniform', label: 'Stays where the flowing holdup left it' }, { value: 'lowpoints', label: 'Drains to the low points of the profile' }] },
    { key: 'gasWaterSatPct', label: 'Water dissolved in the gas (saturation at the inlet)', unit: '%', value: 100, min: 0, max: 100, help: '100 % = the gas is saturated with water vapour; a dehydrated gas carries less and builds a thinner hydrate film on the cold gas-wetted wall.' },
  ] },
  { group: 'Fluid system and line', tab: 'inputs', help: 'Leave these at the linked values to stay consistent with the other suites. Any change makes this suite recompute its own steady flow picture with the shared kernel.', fields: [
    { key: 'fluidSystem', label: 'Fluid system', type: 'select', value: 'case', options: [{ value: 'case', label: 'Case fluid and rates' }, { value: 'wetgas', label: 'Lean wet gas (gas-dominated)' }, { value: 'highwc', label: 'Case fluid at a high water cut' }] },
    { key: 'gasRate', label: 'Gas rate', unit: 'million Sm³/d', value: 3, min: 0.1, max: 30, showIf: (v) => v.fluidSystem === 'wetgas' },
    { key: 'gasWater', label: 'Free water rate', unit: 'Sm³/d', value: 20, min: 0, max: 2000, showIf: (v) => v.fluidSystem === 'wetgas' },
    { key: 'highWc', label: 'Water cut', unit: '%', value: 70, min: 1, max: 95, showIf: (v) => v.fluidSystem === 'highwc' },
    { key: 'idMm', label: 'Inner diameter', unit: 'mm', value: BASE.idMm, min: 25, max: 1500 },
    { key: 'roughUm', label: 'Clean wall roughness', unit: 'µm', value: BASE.roughUm, min: 0.5, max: 3000 },
    { key: 'uValue', label: 'Overall heat-transfer coefficient (on ID)', unit: 'W/m²K', value: BASE.U, min: 0.2, max: 200, typical: [1, 25], help: 'Controls the heat flux to the wall (wax), the cooldown time and how fast heat of formation can leave.' },
    { key: 'tSeabed', label: 'Seabed / ambient temperature', unit: '°C', value: BASE.tSeabed, min: -5, max: 40 },
    { key: 'lengthScale', label: 'Route length multiplier', unit: '–', value: 1, min: 0.1, max: 10, help: '1 = the case route. Larger values stretch the horizontal distance (long tie-back study).' },
    { key: 'thermalMass', label: 'Thermal mass of pipe and contents', unit: 'kJ/m/K', value: 0, min: 0, max: 5000, help: '0 = computed from the fluid in place and the steel wall (insulation storage neglected, which is conservative for cooldown).' },
    { key: 'dep0Mm', label: 'Initial hydrate deposit thickness', unit: 'mm', value: 0, min: 0, max: 100 },
    { key: 'inletHydPct', label: 'Hydrate particles entering at the inlet', unit: 'vol % of liquid', value: 0, min: 0, max: 30 },
  ] },
  { group: 'Thermodynamic inhibition', tab: 'inputs', fields: [
    { key: 'inhibitor', label: 'Inhibitor in the aqueous phase', type: 'select', value: 'case', options: [{ value: 'case', label: 'As the case fluid' }, { value: 'none', label: 'None' }, { value: 'MeOH', label: 'Methanol' }, { value: 'MEG', label: 'MEG' }, { value: 'DEG', label: 'DEG' }, { value: 'TEG', label: 'TEG' }, { value: 'EtOH', label: 'Ethanol' }] },
    { key: 'inhWt', label: 'Inhibitor concentration', unit: 'wt % of aqueous phase', value: 0, min: 0, max: 90, showIf: (v) => v.inhibitor !== 'case' && v.inhibitor !== 'none' },
    { key: 'inhEff', label: 'Inhibitor effectiveness', unit: '% of ideal depression', value: 100, min: 10, max: 150, help: 'Scales the hydrate-temperature depression of the inhibitor (calibrate against rocking-cell or autoclave data; below 100 % for poor mixing or lean inhibitor).' },
    { key: 'marginC', label: 'Required margin outside the hydrate region', unit: '°C', value: 3, min: 0, max: 15 },
    { key: 'meohVapK', label: 'Methanol loss to gas', unit: 'kg per million Sm³ per wt %', value: 16, min: 0, max: 200, help: 'Screening partition coefficient at cold high-pressure conditions.' },
    { key: 'meohOilK', label: 'Methanol loss to hydrocarbon liquid', unit: 'kg/kg oil per mass fraction', value: 0.004, min: 0, max: 0.1 },
    { key: 'khiLimit', label: 'Kinetic inhibitor subcooling limit (screening)', unit: '°C', value: 10, min: 3, max: 20 },
    { key: 'aaWcLimit', label: 'Anti-agglomerant water-cut limit (screening)', unit: '%', value: 50, min: 10, max: 90 },
    { key: 'injections', label: 'Chemical-injection points along the line', type: 'table', columns: [{ key: 'x', label: 'Distance from the inlet', unit: 'km' }, { key: 'rate', label: 'Inhibitor rate', unit: 'm³/d' }], value: [{ x: 0, rate: 0 }], help: 'Each row adds inhibitor to the water from that point downstream (rows with a zero rate are ignored).' },
    { key: 'injChem', label: 'Injected chemical when the case fluid carries none', type: 'select', value: 'MEG', options: [{ value: 'MeOH', label: 'Methanol' }, { value: 'MEG', label: 'MEG' }, { value: 'DEG', label: 'DEG' }, { value: 'TEG', label: 'TEG' }, { value: 'EtOH', label: 'Ethanol' }] },
    { key: 'hydEqOffset', label: 'Hydrate equilibrium temperature offset', unit: '°C', value: 0, min: -5, max: 5, help: 'Shifts the hydrate curve of the fluid (calibration against measured dissociation points).' },
  ] },
  { group: 'Preservation of the shut-in line', tab: 'inputs', showIf: shut, help: 'The design case of a shutdown or cold restart is the preserved line. The same event without preservation is marched as well and reported beside it, so the benefit and the residual risk are both visible.', fields: [
    { key: 'preserve', label: 'Preservation strategy', type: 'select', value: 'inhibit', options: [{ value: 'inhibit', label: 'Inhibitor placed in the line before it cools' }, { value: 'depressurise', label: 'Depressurise, restart with inhibitor' }, { value: 'heat', label: 'Hold the line warm (heating)' }, { value: 'none', label: 'None (unprotected)' }], help: 'Offered from the operations suite when it has been run.' },
    { key: 'preserveChem', label: 'Preservation chemical', type: 'select', value: 'MEG', options: [{ value: 'MEG', label: 'MEG' }, { value: 'MeOH', label: 'Methanol' }], showIf: (v) => v.preserve === 'inhibit' || v.preserve === 'depressurise', help: 'Used when the aqueous phase of the case carries no inhibitor of its own.' },
    { key: 'preserveDose', label: 'Preservation dose in the water', unit: 'wt % of aqueous phase', value: 0, min: 0, max: 80, showIf: (v) => v.preserve === 'inhibit' || v.preserve === 'depressurise', help: '0 = sized here: the dose that removes the peak subcooling of the unprotected event plus the required margin. The operations suite publishes its own dose, which can be linked.' },
    { key: 'preserveLeadH', label: 'Time needed to put the preservation in place', unit: 'h', value: 4, min: 0, max: 200, help: 'Decision time plus the displacement, depressurisation or heat-up time. It is compared with the cooldown time of the line.' },
  ] },
  { group: 'Heating and equipment boundaries', tab: 'inputs', help: 'A directly heated section acts inside the time march; equipment rows retain a fraction of the solids that pass (choke, low point, strainer, dead leg).', fields: [
    { key: 'heatWm', label: 'Heating power per metre', unit: 'W/m', value: 0, min: 0, max: 2000, typical: [20, 150], help: '0 = no heating. Direct electrical heating or heat tracing of the section below.' },
    { key: 'heatFromKm', label: 'Heated section from', unit: 'km', value: 0, min: 0, max: 1000 },
    { key: 'heatToKm', label: 'Heated section to', unit: 'km', value: 1000, min: 0, max: 1000 },
    { key: 'heatStartH', label: 'Heating starts after', unit: 'h', value: 0, min: 0, max: 2000 },
    { key: 'equip', label: 'Equipment that retains solids', type: 'table', columns: [{ key: 'x', label: 'Distance from the inlet', unit: 'km' }, { key: 'eff', label: 'Solids retained', unit: '% of the passing solids' }, { key: 'name', label: 'Item', type: 'text' }], value: [{ x: 0, eff: 0, name: 'none' }], help: 'Rows with zero retention are ignored.' },
  ] },
  { group: 'Hydrate nucleation and growth', tab: 'setup', help: 'Kinetic constants are system-specific: fit them on the Calibration tab to flow-loop, rocking-cell or autoclave data before using the results for design.', fields: [
    { key: 'regime', label: 'Hydrate system', type: 'select', value: 'auto', options: [{ value: 'auto', label: 'Automatic from the flow picture' }, { value: 'oil', label: 'Oil-dominated (water-in-oil emulsion, shrinking core)' }, { value: 'gas', label: 'Gas-dominated (wall film and entrained water)' }, { value: 'water', label: 'Water-dominated (gas absorption limited)' }] },
    { key: 'nucleation', label: 'Nucleation', type: 'select', value: 'heterogeneous', options: [{ value: 'heterogeneous', label: 'Heterogeneous: measured rate at the gas–water interface' }, { value: 'classical', label: 'Classical theory with a contact-angle factor' }, { value: 'homogeneous', label: 'Classical theory, homogeneous' }], help: 'The measured-rate model reproduces the onset statistics of stirred-cell and rocking-cell tests (two populations of nucleation sites, rate per unit of gas–water interface). Classical theory with bulk properties is kept for comparison: it nucleates far too late.' },
    { key: 'nucB', label: 'Nucleation barrier constant B′ of the first site population', unit: 'K³', value: 2660, min: 100, max: 1e6, showIf: (v) => v.nucleation === 'heterogeneous', help: 'J ∝ exp(−B′/(T ΔT²)). 2 660 K³ from isothermal induction times of methane hydrate at 2–4 K subcooling; the second population (3.1 × 10⁴ K³) takes over above about 5 K.' },
    { key: 'nucArea', label: 'Gas–water interface of a test cell', unit: 'cm²', value: 11, min: 0.01, max: 1e6, help: 'Used for the induction-time curve and for onset and induction measurements entered for comparison (11 cm² is the cell in which the rate constants were measured). In the line the interface of each cell is used.' },
    { key: 'nucA', label: 'Nucleation rate multiplier, log₁₀', unit: 'log₁₀', value: 0, min: -99, max: 40, help: 'Multiplies the nucleation rate (site density, agitation, memory water). 0 = the measured rate constants. With the classical options it multiplies the pre-exponential 10^7.5 (contact-angle) or 10³⁵ (homogeneous) 1/m³/s.' },
    { key: 'contactAngle', label: 'Contact angle of the nucleus on the substrate', unit: '°', value: 40, min: 5, max: 180, showIf: (v) => v.nucleation === 'classical' },
    { key: 'sigmaHW', label: 'Hydrate–water interfacial energy', unit: 'mJ/m²', value: 32, min: 5, max: 60, showIf: (v) => v.nucleation !== 'heterogeneous' },
    { key: 'nucVolume', label: 'Water sample volume for the induction time', unit: 'L', value: 1, min: 0.001, max: 1000, showIf: (v) => v.nucleation !== 'heterogeneous', help: 'Classical theory only: the induction time is the mean waiting time for the first nucleus in this much water.' },
    { key: 'kinK', label: 'Intrinsic rate constant at 4 °C', unit: '10⁻¹⁰ mol/m²/Pa/s', value: 0.25, min: 0, max: 1e4, typical: [0.05, 20], help: 'Englezos–Bishnoi constant on the fugacity difference: 0.055–0.065 in the original methane work, 0.21–0.31 after the later correction for the particle surface (default). 0 switches kinetics off.' },
    { key: 'kinEa', label: 'Activation temperature E/R', unit: 'K', value: 13600, min: 0, max: 30000 },
    { key: 'shellD', label: 'Diffusivity through the hydrate shell', unit: '10⁻¹⁵ m²/s', value: 0.8, min: 0.001, max: 1e9, help: 'Effective diffusivity of the hydrate former through the shell on a converting droplet. Line simulators use annealed-film values of 0.01–0.1; fresh films measured by calorimetry give 300–800. The default reproduces the day-long conversion of a stirred water-in-oil autoclave with 40 µm droplets.' },
    { key: 'mtMult', label: 'Gas-absorption (mass-transfer) multiplier', unit: '–', value: 1, min: 0.001, max: 1000 },
    { key: 'htMult', label: 'Heat-removal multiplier', unit: '–', value: 1, min: 0.05, max: 20, help: 'Scales the heat that can be removed from a forming slurry before it reaches the hydrate temperature.' },
    { key: 'hydNumber', label: 'Hydration number', unit: 'mol water/mol gas', value: 0, min: 0, max: 9, help: '0 = from the Langmuir cage occupancy of methane structure I at line conditions (van der Waals–Platteeuw).' },
    { key: 'disK0', label: 'Dissociation constant K₀ (Kim–Bishnoi)', unit: '10⁴ mol/m²/Pa/s', value: 3.6, min: 0, max: 1e4, help: 'Arrhenius pre-exponential of the dissociation flux K₀ exp(−E/RT)(f_eq − f); 3.6 and 81 kJ/mol are the methane values of Clarke and Bishnoi.' },
    { key: 'disE', label: 'Dissociation activation energy', unit: 'kJ/mol', value: 81, min: 10, max: 200 },
    { key: 'rhoHyd', label: 'Hydrate particle density', unit: 'kg/m³', value: 914, min: 800, max: 1100 },
    { key: 'langmuirModel', label: 'Langmuir constants of the statistical model', type: 'select', value: 'parrish', options: [{ value: 'parrish', label: 'Parrish–Prausnitz fit' }, { value: 'kihara', label: 'Kihara cell potential' }] },
    { key: 'slugCouple', label: 'Couple gas absorption to the slug unit cell', type: 'bool', value: true, help: 'In slug flow the gas-to-liquid conductance is taken from the slug body (entrained bubbles, strong mixing) and the film zone instead of the stratified interface.' },
    { key: 'inhDiff', label: 'Inhibitor diffusivity in water', unit: '10⁻⁹ m²/s', value: 1.3, min: 0.1, max: 10, help: 'Used for the mass-transfer-controlled dissociation of a plug in contact with inhibitor.' },
    { key: 'inhFilmK', label: 'Inhibitor film coefficient at the plug face', unit: '10⁻⁵ m/s', value: 1, min: 0.001, max: 1000 },
  ] },
  { group: 'Particles, agglomeration and slurry', tab: 'setup', fields: [
    { key: 'primaryUm', label: 'Primary particle size', unit: 'µm', value: 0, min: 0, max: 2000, help: '0 = droplet size from the Boxall inertial correlation, d/D = 0.063 We^−3/5.' },
    { key: 'primarySigma', label: 'Geometric standard deviation of the primary size', unit: '–', value: 1, min: 1, max: 4, help: '1 = all primary particles in one size class; larger values seed a log-normal size distribution around the primary size.' },
    { key: 'cohesion', label: 'Cohesive force per unit particle size', unit: 'mN/m', value: 2, min: 0, max: 200, typical: [0.1, 50], help: 'Micromechanical force between hydrate particles divided by their size. Anti-agglomerants reduce it by one to two orders of magnitude.' },
    { key: 'aggEff', label: 'Collision (sticking) efficiency', unit: '–', value: 0.05, min: 0, max: 1 },
    { key: 'kBreak', label: 'Breakage coefficient', unit: '–', value: 0.04, min: 0, max: 10, help: 'Breakage frequency = coefficient × shear rate × (size / cohesive-limit size)³.' },
    { key: 'fractal', label: 'Fractal dimension of agglomerates', unit: '–', value: 2.5, min: 1.8, max: 2.95 },
    { key: 'phiMax', label: 'Maximum packing fraction', unit: '–', value: 0.571, min: 0.3, max: 0.74 },
    { key: 'viscModel', label: 'Slurry viscosity model', type: 'select', value: 'mills', options: [{ value: 'mills', label: 'Mills' }, { value: 'krieger', label: 'Krieger–Dougherty' }, { value: 'thomas', label: 'Thomas' }] },
  ] },
  { group: 'Wall deposition and plugging', tab: 'setup', fields: [
    { key: 'adhesion', label: 'Wall-capture efficiency', unit: '–', value: 0.005, min: 0, max: 1, typical: [0.0005, 0.05], help: 'Probability that a particle reaching a sub-cooled wall stays there, before the force-balance reduction.' },
    { key: 'adhForce', label: 'Wall adhesion force per unit particle size', unit: 'mN/m', value: 5, min: 0.001, max: 500 },
    { key: 'filmMult', label: 'Vapour-film deposition multiplier', unit: '–', value: 1, min: 0, max: 100, help: 'Scales hydrate growth from water vapour condensing on the gas-wetted cold wall.' },
    { key: 'tauCrit', label: 'Critical wall shear for removal (fresh deposit)', unit: 'Pa', value: 30, min: 0.1, max: 5000, help: 'Above this shear the deposit is eroded and sloughs back into the stream. It rises as the deposit ages and densifies.' },
    { key: 'kRemove', label: 'Removal / detachment rate constant', unit: '1/h', value: 0.5, min: 0, max: 100 },
    { key: 'porosity0', label: 'Porosity of a fresh deposit', unit: '–', value: 0.6, min: 0.05, max: 0.95 },
    { key: 'porosityInf', label: 'Porosity of an aged deposit', unit: '–', value: 0.2, min: 0.01, max: 0.9 },
    { key: 'ageHours', label: 'Ageing time constant', unit: 'h', value: 12, min: 0.1, max: 2000 },
    { key: 'plugGrainUm', label: 'Grain size for plug permeability', unit: 'µm', value: 0, min: 0, max: 5000, help: '0 = the primary particle size. Used in the Kozeny–Carman permeability.' },
    { key: 'plugBlockPct', label: 'Area blockage that counts as a plug', unit: '%', value: 90, min: 30, max: 99.5 },
    { key: 'plugLength', label: 'Plug length for remediation estimates', unit: 'm', value: 50, min: 1, max: 5000 },
    { key: 'depressP', label: 'Depressurisation pressure for plug melting', unit: 'bara', value: 5, min: 1, max: 100 },
  ] },
  { group: 'Wax', tab: 'setup', fields: [
    { key: 'wat', label: 'Wax appearance temperature', unit: '°C', value: 45, min: -20, max: 90 },
    { key: 'waxContent', label: 'Wax content of the oil', unit: 'wt %', value: 5, min: 0, max: 40 },
    { key: 'waxThermo', label: 'Wax thermodynamics', type: 'select', value: 'sle', options: [{ value: 'sle', label: 'Solid–liquid equilibrium of the n-paraffin distribution (computed WAT and solubility curve)' }, { value: 'input', label: 'Entered WAT, wax content and solubility slope' }], help: 'The computed option splits the C7+ fraction of the case fluid into single carbon numbers and solves the ideal-solution equilibrium; it falls back to the entered values when the fluid has no wax-forming fraction.' },
    { key: 'waxDetect', label: 'Solid wax that defines the cloud point', unit: 'wt %', value: 0.02, min: 0.001, max: 1, showIf: (v) => v.waxThermo === 'sle' },
    { key: 'waxHf', label: 'Heat-of-fusion multiplier', unit: '–', value: 1, min: 0.5, max: 2, showIf: (v) => v.waxThermo === 'sle' },
    { key: 'waxSlope', label: 'Solubility slope below the WAT', unit: '1/K', value: 0.04, min: 0.005, max: 0.2, help: 'Dissolved wax = total × exp(−slope × (WAT − T)).' },
    { key: 'waxStripC', label: 'Shear-stripping coefficient', unit: '–', value: 0.055, min: 0, max: 5, help: 'Coefficient C₂ of the Matzain shear-stripping factor 1/(1 + C₂ N_SR^1.4).' },
    { key: 'waxAspect', label: 'Wax crystal aspect ratio in the gel', unit: '–', value: 8, min: 1, max: 60, help: 'α of the hindered diffusivity D/(1 + α²F²/(1 − F)) that drives ageing.' },
    { key: 'waxDiff', label: 'Diffusivity correlation', type: 'select', value: 'haydukMinhas', options: [{ value: 'haydukMinhas', label: 'Hayduk–Minhas' }, { value: 'wilkeChang', label: 'Wilke–Chang' }] },
    { key: 'waxMult', label: 'Wax deposition multiplier', unit: '–', value: 1, min: 0, max: 100, help: 'Multiplies the diffusive flux (fit to flow-loop or cold-finger data).' },
    { key: 'waxModel', label: 'Wax deposition model', type: 'select', value: 'matzain', options: [{ value: 'matzain', label: 'Matzain (enhanced molecular diffusion with shear stripping)' }, { value: 'mechanistic', label: 'Mechanistic (diffusion, shear dispersion, Brownian flux, ageing)' }], help: 'The Matzain model multiplies the Fick flux by the empirical constant C1 = 15 divided by the wax fraction of the deposit. The mechanistic option has no such constant and must be fitted with the multiplier.' },
    { key: 'waxC1', label: 'Matzain enhancement constant C1', unit: '–', value: 15, min: 0.1, max: 200, showIf: (v) => v.waxModel === 'matzain' },
    { key: 'waxOil', label: 'Oil trapped in the deposit', unit: 'vol %', value: 0, min: 0, max: 97, help: '0 = from the Reynolds-number closure of the Matzain model, 100 (1 − Re^0.15/8), or 80 % for the mechanistic model. Measured deposits hold 40–95 % oil.' },
    { key: 'waxInhOn', label: 'Wax inhibitor / pour-point depressant injected', type: 'bool', value: false, help: 'Scales the deposition flux and the gel strength by the efficiencies below (from cold-finger and gel-strength tests of the chemical on the real oil).' },
    { key: 'waxInhEff', label: 'Reduction of the deposition rate by the inhibitor', unit: '%', value: 50, min: 0, max: 99, showIf: (v) => v.waxInhOn },
    { key: 'waxInhGel', label: 'Reduction of the gel yield stress by the inhibitor', unit: '%', value: 60, min: 0, max: 99, showIf: (v) => v.waxInhOn },
    { key: 'waxK', label: 'Deposit thermal conductivity', unit: 'W/m/K', value: 0.25, min: 0.05, max: 1 },
    { key: 'waxLimitMm', label: 'Thickness that triggers pigging', unit: 'mm', value: 2, min: 0.1, max: 50 },
    { key: 'pourOffset', label: 'Gelling starts this far below the WAT', unit: '°C', value: 15, min: 0, max: 60 },
    { key: 'gelCoef', label: 'Gel yield-stress coefficient', unit: 'Pa per (wt % solid wax)²', value: 4, min: 0, max: 500 },
  ] },
  { group: 'Scale', tab: 'setup', fields: [
    { key: 'water', label: 'Produced-water analysis', type: 'table', columns: ION_IDS.map((id) => ({ key: id, label: id, unit: 'mg/L' })), value: [{ ...WATER0 }], help: 'First row is used. Bicarbonate is the alkalinity.' },
    { key: 'co2Pct', label: 'CO₂ in the gas phase', unit: 'mol %', value: 3, min: 0.001, max: 60 },
    { key: 'swFrac', label: 'Seawater fraction in the produced water', unit: '%', value: 0, min: 0, max: 100, help: 'Injection-water breakthrough: produced water is mixed with standard seawater before the indices are evaluated.' },
    { key: 'srpOn', label: 'Sulphate removal on the injected seawater', type: 'bool', value: true, help: 'The seawater that breaks through carries the sulphate left by a sulphate-removal (nanofiltration) plant instead of the 2 712 mg/L of raw seawater.' },
    { key: 'srpSO4', label: 'Sulphate in the treated injection water', unit: 'mg/L', value: 40, min: 0, max: 2712, showIf: (v) => v.srpOn },
    { key: 'scaleInhOn', label: 'Scale inhibitor injected', type: 'bool', value: true, help: 'Threshold inhibitor injected upstream of the first supersaturated point (downhole or at the tree). It lengthens the nucleation induction time; the dose is compared with the minimum inhibitor concentration.' },
    { key: 'scaleInh', label: 'Scale inhibitor', type: 'select', value: 'DTPMP', options: ['DTPMP', 'BHPMP', 'NTMP', 'HEDP', 'PPCA', 'PASP', 'PVS', 'SPCA'].map((k) => ({ value: k, label: SCALE_INHIBITORS[k].name })), showIf: (v) => v.scaleInhOn },
    { key: 'scaleInhDose', label: 'Scale-inhibitor dose (active)', unit: 'mg/L of water', value: 5, min: 0, max: 500, showIf: (v) => v.scaleInhOn, help: 'Typical continuous doses are 1–20 mg/L of active inhibitor.' },
    { key: 'scaleProtectH', label: 'Time the brine must stay scale-free', unit: 'h', value: 24, min: 0.1, max: 2000, help: 'Residence time from the first supersaturated point to disposal (flowline, separators, produced-water system). 2 h and 24 h are the usual design values.' },
    { key: 'actModel', label: 'Activity-coefficient model', type: 'select', value: 'pitzer', options: [{ value: 'pitzer', label: 'Pitzer ion interaction (PHREEQC pitzer.dat)' }, { value: 'truesdellJones', label: 'Ion association, Truesdell–Jones (PHREEQC phreeqc.dat)' }, { value: 'davies', label: 'Ion association, Davies' }], help: 'Pitzer is valid to halite saturation; the ion-association models to about 1–2 mol/kg.' },
    { key: 'scaleK', label: 'Scale surface-reaction constant at 25 °C', unit: '10⁻⁸ mol/m²/s', value: 1.26, min: 0, max: 1e4, help: 'Neutral-mechanism rate constant of barite (10^−7.90 mol/m²/s at 25 °C, activation energy 30.8 kJ/mol) from the USGS compilation of mineral rate parameters; other minerals react faster and are then limited by ion transport to the wall.' },
    { key: 'scaleSigma', label: 'Crystal–water interfacial energy', unit: 'mJ/m²', value: 79, min: 10, max: 200, help: 'Controls the nucleation rate in the bulk water (classical nucleation theory).' },
    { key: 'scaleNucA', label: 'Nucleation pre-exponential, log₁₀', unit: 'log₁₀(1/m³/s)', value: 30, min: 5, max: 40 },
    { key: 'scaleKg', label: 'Crystal growth constant at 25 °C', unit: '10⁻¹⁰ m/s', value: 1, min: 0, max: 1e5, help: 'Linear growth rate G = k (S − 1)² (parabolic law).' },
    { key: 'scaleDepEff', label: 'Crystal sticking efficiency at the wall', unit: '–', value: 0.01, min: 0, max: 1 },
    { key: 'scaleSeedLog', label: 'Seed particles in the water, log₁₀', unit: 'log₁₀(1/m³)', value: 9, min: 0, max: 14, help: 'Suspended fines of 1 µm on which the scale mineral grows (0 = none: only homogeneous nucleation, which is negligible below a saturation index of about 2).' },
  ] },
  { group: 'Asphaltene', tab: 'setup', fields: [
    { key: 'saraSat', label: 'Saturates', unit: 'wt %', value: 55, min: 0, max: 100 }, { key: 'saraAro', label: 'Aromatics', unit: 'wt %', value: 28, min: 0, max: 100 }, { key: 'saraRes', label: 'Resins', unit: 'wt %', value: 14, min: 0, max: 100 }, { key: 'saraAsp', label: 'Asphaltenes', unit: 'wt %', value: 3, min: 0, max: 100 },
    { key: 'asphDelta', label: 'Asphaltene solubility parameter at 25 °C', unit: 'MPa^½', value: 20, min: 17, max: 24 }, { key: 'asphMV', label: 'Asphaltene molar volume', unit: 'm³/kmol', value: 2, min: 0.3, max: 10 },
    { key: 'asphMethod', label: 'Oil solubility parameter', type: 'select', value: 'density', options: [{ value: 'density', label: 'Density correlation' }, { value: 'eos', label: 'Equation of state (cohesive energy of the flashed liquid)' }], help: 'The equation-of-state option flashes the case fluid along the depletion path and takes the change of the oil solubility parameter and molar volume from the cohesive energy of the liquid, anchored to the density correlation at reservoir conditions.' },
    { key: 'aspInhOn', label: 'Asphaltene inhibitor / dispersant injected', type: 'bool', value: true, help: 'Reduces the deposition tendency by the efficiency below; the screening indices of the oil are reported unchanged beside the managed risk.' },
    { key: 'aspInhEff', label: 'Deposit reduction by the inhibitor', unit: '%', value: 70, min: 0, max: 99, showIf: (v) => v.aspInhOn, help: 'From a dispersant or deposition test of the chemical on the real oil at the planned dose.' },
    { key: 'asphOnsetRef', label: 'Measured upper onset pressure at reservoir temperature', unit: 'bara', value: 0, min: 0, max: 1500, help: '0 = none. When given, the asphaltene solubility parameter is adjusted so that the model reproduces it (solid-phase reference state).' },
  ] },
  { group: 'Sand and solids', tab: 'setup', fields: [
    { key: 'sandRate', label: 'Sand production', unit: 'kg/d', value: 50, min: 0, max: 1e5 }, { key: 'sandUm', label: 'Sand particle size', unit: 'µm', value: 150, min: 10, max: 3000 }, { key: 'sandRho', label: 'Sand density', unit: 'kg/m³', value: 2650, min: 1100, max: 5000 },
    { key: 'sandShape', label: 'Drag multiplier for particle shape', unit: '–', value: 1, min: 0.5, max: 5, help: '1 = spheres. Angular grains settle more slowly (larger drag coefficient).' },
    { key: 'sandCoef', label: 'Critical-velocity multiplier', unit: '–', value: 1, min: 0.2, max: 5, help: 'Multiplies the minimum transport velocity of the correlations (fit to transport tests).' },
    { key: 'sandShields', label: 'Critical Shields number for resuspension', unit: '–', value: 0, min: 0, max: 0.5, help: '0 = from the grain size and the liquid viscosity (Soulsby–Whitehouse threshold curve: 0.055 for coarse grains in water, up to 0.30 in viscous liquids).' },
    { key: 'sandEntrain', label: 'Bed entrainment coefficient', unit: '10⁻⁴ m/s', value: 1, min: 0, max: 1e4, help: 'Solids volume flux picked up from a bed per unit of excess Shields stress.' },
    { key: 'erosionSm', label: 'Erosion geometry constant', unit: '–', value: 5.5, min: 0.5, max: 100, help: 'Constant of the Salama screening relation (5.5 for elbows).' },
    { key: 'sandHours', label: 'Transient sand transport period', unit: 'h', value: 24, min: 0.5, max: 2000 },
  ] },
  { group: 'Measurements for comparison', tab: 'setup', help: 'Optional. Enter laboratory or field measurements; the run adds the model prediction, a parity plot and error metrics. Types: equilibrium (a = pressure bara → hydrate temperature °C) · inhibitor (a = wt %, b = 1 methanol / 2 MEG / 3 ethanol → depression °C) · onset (a = cooling rate °C/h, b = gas–water interface cm², c = hydrate temperature °C → median subcooling at onset °C) · induction (a = subcooling °C, b = interface cm², c = hydrate temperature °C → mean induction time h) · autoclave, rockingcell, flowloop, growth (a = time h, b = subcooling °C, c = velocity m/s → water conversion %) · loopdp (same → pressure-drop ratio) · plugging (a = subcooling °C, b = velocity m/s → time to plug h) · psd, agglomeration (a = time h, b = subcooling °C, c = velocity m/s → d43 µm) · rheology (a = hydrate vol %, b = shear rate 1/s → relative viscosity) · deposition (a = time h, b = subcooling °C, c = velocity m/s → deposit mm) · adhesion (a = particle size µm → force µN) · dissociation (a = time min, b = °C above equilibrium, c = particle size µm → % dissociated) · restart (a = shut-in h, b = subcooling °C, c = restart velocity m/s → peak pressure-drop ratio) · field-plugtime (→ h), field-plugx (→ km), field-dp (→ bar; compared with this run) · waxloop (a = time h, b = oil − coolant °C, c = velocity m/s → mm) · scale-barite, scale-celestite, scale-gypsum, scale-anhydrite, scale-calcite (a = NaCl mol/kg, b = °C, c = bara → µmol/kg, calcite at the CO₂ fraction of the scale inputs) · asphaltene (a = temperature °C → upper onset pressure bara) · sand (a = particle µm, b = pipe diameter mm, c = liquid viscosity mPa·s → critical velocity m/s) · settling (a = particle µm, b = particle density kg/m³, c = liquid viscosity mPa·s → mm/s) · erosion (a = sand kg/d, b = velocity m/s → mm/y).', fields: [
    { key: 'expData', label: 'Measurements', type: 'table', columns: [{ key: 'type', label: 'Type', type: 'text' }, { key: 'a', label: 'a' }, { key: 'b', label: 'b' }, { key: 'c', label: 'c' }, { key: 'meas', label: 'Measured' }], value: [{ type: '', a: 0, b: 0, c: 0, meas: 0 }] },
  ] },
  { group: 'Discretisation', tab: 'mesh', help: 'Axial cells, time step of the hydrate march and number of particle-size classes; the Monte Carlo runs use a coarser copy of the same model.', fields: [
    { key: 'nAxial', label: 'Axial cells', value: 36, min: 8, max: 120 },
    { key: 'dtMin', label: 'Time step of the hydrate march', unit: 'min', value: 20, min: 0.5, max: 240 },
    { key: 'nClasses', label: 'Particle-size classes', value: 14, min: 6, max: 40 },
    { key: 'nMC', label: 'Monte Carlo samples for plugging probability', value: 20, min: 0, max: 400 },
    { key: 'mcParticles', label: 'Simulated particles of the Monte Carlo population balance', value: 500, min: 0, max: 50000, help: '0 = skip. Stochastic cross-check of the sectional solver in the worst cell.' },
    { key: 'trackOn', label: 'Run the Lagrangian particle tracker', type: 'bool', value: false, help: 'Tracks stochastic parcels (Maxey–Riley equation with turbulent dispersion and wall rules) through the line and reports where they are captured.' },
    { key: 'trackWhat', label: 'Tracked solid', type: 'select', value: 'hydrate', options: [{ value: 'hydrate', label: 'Hydrate particles of the peak size distribution' }, { value: 'sand', label: 'Sand grains' }], showIf: (v) => v.trackOn },
    { key: 'nParcels', label: 'Parcels', value: 40, min: 4, max: 5000, showIf: (v) => v.trackOn },
    { key: 'trackStepPct', label: 'Tracker step', unit: '% of the cross-pipe mixing time', value: 5, min: 0.2, max: 25, showIf: (v) => v.trackOn },
    { key: 'eeSteps', label: 'Time steps of the Eulerian solids transport', value: 60, min: 10, max: 4000 },
    { key: 'seed', label: 'Random seed', value: 7, min: 1, max: 1e6 },
    { key: 'cfdOn', label: 'Population balance on the resolved flow field', type: 'bool', value: true, help: 'Solves the size classes on the axisymmetric velocity, shear and eddy-viscosity field of the pipe section with the most hydrate: radial transport by turbulent diffusion, aggregation and breakage at the local shear, deposition as a wall flux.' },
    { key: 'cfdNr', label: 'Radial cells of the flow-field population balance', value: 8, min: 4, max: 80, showIf: (v) => v.cfdOn },
    { key: 'cfdNx', label: 'Axial steps of the flow-field population balance', value: 12, min: 4, max: 400, showIf: (v) => v.cfdOn },
    { key: 'cfdLenD', label: 'Resolved pipe length', unit: 'diameters', value: 200, min: 10, max: 5000, showIf: (v) => v.cfdOn },
  ] },
];
/** Flat list of the input field declarations (key, label, unit, default, limits). */
export const INPUT_FIELDS = INPUTS.flatMap((g) => g.fields);
const FIELDS = INPUT_FIELDS, DEF = Object.fromEntries(FIELDS.map((f) => [f.key, f.value]));
/** Fill defaults, coerce and check the inputs; throws a readable Error for impossible input. */
function clean(v0 = {}) {
  const v = {};
  for (const f of FIELDS) {
    const x = v0[f.key];
    if (f.type === 'select') { v[f.key] = f.options.some((o) => o.value === x) ? x : f.value; continue; }
    if (f.type === 'table') { v[f.key] = Array.isArray(x) && x.length ? x : f.value; continue; }
    const y = x === null || x === undefined || x === '' ? f.value : +x;
    need(Number.isFinite(y), `${f.label} must be a number.`);
    need(f.min === undefined || y >= f.min - 1e-12, `${f.label} must be at least ${f.min}${f.unit && f.unit !== '–' ? ' ' + f.unit : ''} (got ${y}).`);
    need(f.max === undefined || y <= f.max + 1e-12, `${f.label} must not exceed ${f.max}${f.unit && f.unit !== '–' ? ' ' + f.unit : ''} (got ${y}).`);
    v[f.key] = y;
  }
  for (const k of ['nAxial', 'nClasses', 'nMC', 'seed', 'mcParticles', 'nParcels', 'eeSteps', 'cfdNr', 'cfdNx']) v[k] = Math.round(v[k]);
  for (const f of FIELDS) if (f.type === 'bool') v[f.key] = v0[f.key] === undefined || v0[f.key] === null ? f.value : !!v0[f.key];
  need(v.porosityInf <= v.porosity0, 'The aged deposit porosity cannot exceed the fresh deposit porosity.');
  need(v.saraSat + v.saraAro + v.saraRes + v.saraAsp > 0, 'The SARA analysis is empty: enter the saturate, aromatic, resin and asphaltene fractions.');
  need(v.dep0Mm < 0.45 * v.idMm, 'The initial hydrate deposit is thicker than the pipe can hold.');
  need(v.pInMax > 1.5, 'The available inlet pressure must exceed the arrival pressure.');
  return v;
}
/** Phases of the hydrate march for a scenario: [{ dur (s), dt (s), frac }] and the rate fraction of the steady production period. */
function scenarioPhases(v, dtScale = 1) {
  const dt = v.dtMin * 60 * dtScale, sim = v.simHours * 3600, sh = v.shutHours * 3600, dtS = Math.min(Math.max(4 * dt, sh / 60), sh);
  if (v.scenario === 'steady') return { phases: [{ dur: sim, dt, frac: 1 }], prod: 1 };
  if (v.scenario === 'turndown') return { phases: [{ dur: sim, dt, frac: v.turndownPct / 100 }], prod: 1 };
  if (v.scenario === 'shutdown') return { phases: [{ dur: sh, dt: dtS, frac: 0 }], prod: 1 };
  return { phases: [{ dur: sh, dt: dtS, frac: 0 }, { dur: sim, dt, frac: v.restartPct / 100 }], prod: 1 };
}

// =====================================================================================================
// 15. The run
// =====================================================================================================
const r3 = (x, s = 3) => (Number.isFinite(x) ? +x.toPrecision(s) : null), tx = (x, s = 3) => (Number.isFinite(x) ? +x.toPrecision(s) : '—');
const argmax = (a) => { let k = 0; for (let i = 1; i < a.length; i++) if (a[i] > a[k]) k = i; return k; };
const SC_LABEL = Object.fromEntries(SCEN.map((s) => [s.value, s.label]));
/** Hydrate-temperature depression (°C) of the aqueous phase of a set-up, with the inhibitor effectiveness applied to the inhibitor part. */
const aqDepression = (S, eff = 1, wt = S.aq.inhWt) => { const salt = hydrateDepression({ S: S.aq.S, inhWt: 0, inh: S.aq.inh }); return wt > 0 ? salt + eff * (hydrateDepression({ S: S.aq.S, inhWt: wt, inh: S.aq.inh }) - salt) : salt; };

const waxMemo = new Map();
/** Wax solid–liquid equilibrium curve of a fluid specification (cached: the split and the flash do not change between runs of the same fluid). */
function waxCurveFor(spec, v) {
  if (!spec) return null;
  const key = JSON.stringify([spec.comp, spec.c7MW, spec.c7SG, spec.eos, spec.nPseudo, v.waxDetect, v.waxHf]);
  if (!waxMemo.has(key)) { if (waxMemo.size > 8) waxMemo.clear(); waxMemo.set(key, waxEquilibrium(spec, { detect: v.waxDetect, hfMult: v.waxHf, n: 28 })); }
  return waxMemo.get(key);
}
const oilLike = (spec) => { const c = spec?.comp || {}, tot = Object.values(c).reduce((a, b) => a + (+b || 0), 0); return tot > 0 && (+c.C7p || 0) / tot > 0.02; };

/** Preservation plan of a scenario with a shut-in: { id, active, chem, lead (h), dose (wt %, 0 = size it) }. Flowing scenarios have none. */
export function preservationPlan(v) {
  const id = v.scenario === 'restart' || v.scenario === 'shutdown' ? v.preserve : 'none';
  return { id, active: id !== 'none', chem: v.preserveChem, lead: v.preserveLeadH, dose: v.preserveDose };
}
/**
 * Apply a preservation plan to a set-up and its march options (copies; the unprotected objects stay untouched).
 * inhibit: the water of the whole line carries at least `dose` wt % of inhibitor from the start of the shut-in;
 * depressurise: the shut-in line is lowered to the depressurisation pressure once the lead time has passed and is restarted with the dose;
 * heat: every metre receives the power that holds it `margin` above its hydrate temperature.
 */
export function applyPreservation(S, mx, v, plan, dose = plan.dose) {
  if (!plan.active) return { S, mx };
  if (plan.id === 'heat') {
    const b = S.prof(1), salt = hydrateDepression({ S: S.aq.S, inhWt: 0, inh: S.aq.inh }), tHold = S.hT0(Math.max(...b.P, v.pShut || 0)) - salt + v.marginC;
    const heat = S.tAmb.map((ta, i) => (mx.heat ? mx.heat[i] : 0) + Math.max(S.U * PI * S.D0 * (tHold - ta), 0));
    return { S, mx: { ...mx, heat, heatStart: 0, heatedLength: S.L, holdT: tHold, holdPower: heat.reduce((a, q) => a + q * S.ds, 0) } };
  }
  const inhId = S.aq.inhId !== 'none' ? S.aq.inhId : plan.chem, cell = S.aq.inhWtCell.map((w) => Math.max(w, dose));
  const S1 = { ...S, aq: { ...S.aq, inh: INHIBITORS[inhId], inhId, inhWt: Math.max(S.aq.inhWt, dose), inhWtCell: cell } };
  return { S: S1, mx: { ...mx, inhWt: cell, ...(plan.id === 'depressurise' ? { depress: { t: plan.lead * 3600, P: v.depressP } } : {}) } };
}
/** Headline numbers of one marched hydrate event (used for the unprotected comparison case): subcooling, fraction, plug, risk index. */
export function eventSummary(h, mc, p, pIn0) {
  const maxSub = Math.max(...h.rec.subMax), peakVisc = Math.max(...h.ser.visc, 1), pInPeak = Math.max(...h.ser.pIn.filter((_, k) => h.ser.frac[k] > 0), pIn0), dpRise = pInPeak - pIn0;
  const pOn = h.onset ? 1 : 1 - Math.exp(-Math.max(...h.rec.lam, 0)), sev = Math.max(h.peak.blk / p.plugBlock, Math.log10(peakVisc) / 3, dpRise / Math.max(p.pInMax - pIn0, 1));
  let expo = 0; for (let k = 0; k < h.ser.t.length; k++) if (h.ser.sub[k] > 0) expo += h.ser.t[k] - (k ? h.ser.t[k - 1] : 0);
  return { maxSub, peakPhi: h.peak.phi, peakVisc, peakBlk: h.peak.blk, dpRise, plug: h.plug, onset: h.onset, prob: mc ? (h.plug && mc.n === 0 ? 1 : mc.prob) : h.plug ? 1 : 0, lo: mc ? mc.lo : 0, hi: mc ? mc.hi : 0, n: mc ? mc.n : 0, expo,
    pk: h.peak, stableLen: h.rec.subMax.filter((x) => x > 0).length, risk: h.plug ? 1 : clamp(Math.max(pOn * Math.max(maxSub > 0 ? 0.15 : 0, sev), mc ? mc.prob : 0), 0, 1), subMax: h.rec.subMax.slice() };
}
/** Monte Carlo over uncertain kinetics, cohesion, adhesion and the stochastic (Poisson) onset on a coarse copy of the model. */
async function plugMonteCarlo(v, ctx, S, p, slow, tick, skip = false, prep = null, N = v.nMC) {
  const out = { n: N, plugX: [], plugTimes: [], onsetTimes: [], blk: [], visc: [], prob: 0, lo: 0, hi: 0 };
  if (!(N > 0) || skip) return out;
  const n2 = Math.min(S.n, 10), g2 = pbeGrid(8, p.mode === 'oil' ? p.dPrim / 4 : 2e-6, 2e-2), { phases } = scenarioPhases(v, 3), rg = rng(v.seed);
  let S2 = buildSetup(v, ctx, n2), mx2 = marchExtras(S2, v); if (prep) ({ S: S2, mx: mx2 } = prep(S2, mx2));
  const dW = S2.s.map((s) => interp1(S.s, slow.dWax, s) + interp1(S.s, slow.dScale, s)), Dbase = dW.map((d) => Math.max(S.D0 - 2 * d, 0.2 * S.D0)), ln = (sd) => Math.exp(rg.normal(0, sd));
  for (let k = 0; k < N; k++) {
    const pk = { ...p, kinK: p.kinK * ln(0.7), cohesion: p.cohesion * ln(0.6), adhesion: Math.min(p.adhesion * ln(0.6), 1), adhForce: p.adhForce * ln(0.5), tauCrit: p.tauCrit * ln(0.4), nucA: p.nucA * ln(1.5), nucMult: p.nucMult * ln(1.5), shellD: p.shellD * ln(0.7), lamStar: -Math.log(1 - rg.uniform(0, 1) * 0.999999) };
    const r = runHydrateMarch(S2, pk, { phases, grid: g2, Dbase, dep0: v.dep0Mm / 1000, cheap: true, rows: 2, stopWhenClear: true, ...mx2 });
    if (r.plug) { out.plugTimes.push(r.plug.t / 3600); out.plugX.push(r.plug.x); } if (r.onset) out.onsetTimes.push(r.onset.t / 3600);
    out.blk.push(Math.max(...r.ser.blk, 0)); out.visc.push(Math.max(...r.ser.visc, 1));
    if (k % 8 === 7) await tick(0.7 + (0.2 * k) / N);
  }
  const ph = out.plugTimes.length / N, z = 1.645, den = 1 + (z * z) / N, c = (ph + (z * z) / (2 * N)) / den, hw = (z * Math.sqrt((ph * (1 - ph)) / N + (z * z) / (4 * N * N))) / den; // Wilson 90 % interval
  out.prob = ph; out.lo = Math.max(c - hw, 0); out.hi = Math.min(c + hw, 1);
  return out;
}

async function run(v0, ctx = {}) {
  const v = clean(v0), prog = (f, m) => { try { ctx.progress?.(f, m); } catch { /* progress is optional */ } }, tick = async (f, m) => { if (f !== undefined) prog(f, m || 'Working'); if (ctx.tick) await ctx.tick(); };
  prog(0.02, 'Building the flow picture');
  let S = buildSetup(v, ctx, clamp(v.nAxial, 8, 120)); const n = S.n, D0 = S.D0, ds = S.ds, xkm = S.x.map((a) => a / 1000), A0 = (PI * D0 * D0) / 4;
  if (v.waxThermo === 'sle') S.waxCurve = waxCurveFor(S.fm?.spec, v);
  let mx = marchExtras(S, v);
  const { phases } = scenarioPhases(v), fProd = v.scenario === 'turndown' ? v.turndownPct / 100 : 1, fEvent = phases[phases.length - 1].frac || 1;
  for (const ph of phases) if (ph.frac > 0) S.prof(ph.frac);
  const base = S.prof(1), pb = S.prof(fProd);
  await tick(0.08, 'Wax and scale in steady production');
  const slow = slowDeposits(S, v, fProd);
  const Dbase = slow.dWax.map((d, i) => Math.max(D0 - 2 * (d + slow.dScale[i]), 0.2 * D0)), roughBase = slow.dWax.map((d, i) => Math.min(S.rough0 + 0.05 * d + 0.3 * slow.dScale[i], 0.05 * D0));
  const p = hydrateParams(v, S), oil = p.mode === 'oil', grid = pbeGrid(v.nClasses, oil ? p.dPrim / 4 : 2e-6, 2e-2);
  // ---- hydrate event marched in time
  const phCoarse = scenarioPhases(v, 2).phases, marchOf = async (Sx, mxx, f0, f1, rows, msg, coarse = false) => { const m = hydrateMarch(Sx, p, { phases: coarse === true ? phCoarse : phases, cheap: !!coarse, grid, Dbase, roughBase, dep0: v.dep0Mm / 1000, rows, ...mxx }); let kk = 0; while (!m.done) { m.step(); if (++kk % 10 === 0) await tick(f0 + (f1 - f0) * m.progress, msg); } return m.result(); };
  const idle = (hh) => Math.max(...hh.rec.subMax) <= 0 && !hh.plug && !(v.dep0Mm > 0) && !(v.initHydPct > 0); // without any subcooling no sample can form hydrate
  // The design case of a shut-in is the preserved line. The unprotected event is marched first: it is the comparison case and it sizes an automatic dose.
  const plan = preservationPlan(v), Sun = S, pIn0u = S.prof(fEvent).pIn; let unprot = null, doseUsed = 0, doseCapped = false;
  if (plan.active) {
    const hU = await marchOf(S, mx, 0.1, 0.3, 12, 'Marching the unprotected event (comparison case)', true), mcU = await plugMonteCarlo(v, ctx, S, p, slow, async (f, m) => tick(f === undefined ? f : 0.3 + (f - 0.7), m), idle(hU), null, Math.ceil(v.nMC / 3));
    unprot = eventSummary(hU, mcU, p, pIn0u);
    if (plan.id !== 'heat') {
      const chem = S.aq.inhId !== 'none' ? S.aq.inhId : plan.chem, salt = aqDepression(S, 1, 0), iU = argmax(hU.rec.subMax), now = aqDepression(S, p.inhEff, S.aq.inhWtCell[iU]) - salt, cap = chem === 'MeOH' ? 60 : 70;
      const need = unprot.maxSub + now + v.marginC, auto = need > 0 ? inhibitorFor(salt + need / Math.max(p.inhEff, 0.1), chem, S.aq.S) : 0;
      doseUsed = plan.dose > 0 ? plan.dose : Math.min(auto, cap); doseCapped = !(plan.dose > 0) && auto > cap;
    }
    const plan1 = { ...plan, dose: doseUsed }; ({ S, mx } = applyPreservation(S, mx, v, plan1)); plan.prep = (S2, mx2) => applyPreservation(S2, mx2, v, plan1);
  }
  let h = await marchOf(S, mx, plan.active ? 0.4 : 0.1, 0.6, 48, 'Marching hydrate formation, transport and deposition');
  // An automatic dose is sized on the subcooling of the unprotected event, where the heat of formation holds the fluid near the hydrate
  // temperature; the preserved line releases no such heat and gets colder, so the dose is corrected until the margin is met.
  for (let it = 0; it < 3 && plan.active && plan.id !== 'heat' && !(plan.dose > 0) && !doseCapped; it++) {
    const short = Math.max(...h.rec.subMax) + v.marginC; if (!(short > 0.05)) break; // once no hydrate forms the temperatures no longer depend on the dose, so one correction is normally exact
    const chem = S.aq.inhId, salt = aqDepression(Sun, 1, 0), cap = chem === 'MeOH' ? 60 : 70, have = hydrateDepression({ S: S.aq.S, inhWt: doseUsed, inh: S.aq.inh }) - salt, want = inhibitorFor(salt + have + short / Math.max(p.inhEff, 0.1), chem, S.aq.S);
    doseCapped = want > cap; doseUsed = Math.min(Math.max(want, doseUsed + 0.2), cap);
    const plan1 = { ...plan, dose: doseUsed }; ({ S, mx } = applyPreservation(Sun, marchExtras(Sun, v), v, plan1)); plan.prep = (S2, mx2) => applyPreservation(S2, mx2, v, plan1);
    h = await marchOf(S, mx, 0.6, 0.65, 48, 'Marching the preserved line with the corrected dose');
  }
  const rec = h.rec, ser = h.ser, pk = h.peak, tShut = v.scenario === 'restart' ? v.shutHours : 0;
  await tick(0.68, 'Plugging probability (Monte Carlo)');
  const mc = await plugMonteCarlo(v, ctx, S, p, slow, tick, idle(h), plan.prep || null);
  await tick(0.92, 'Wax, scale, asphaltene, sand and remediation');

  // ---- hydrate driving force, exposure, conversion
  const depNowC = S.aq.inhWtCell.map((w) => aqDepression(S, p.inhEff, w)), depNow = Math.min(...depNowC), depSalt = aqDepression(S, 1, 0), maxSub = Math.max(...rec.subMax), iSub = argmax(rec.subMax), stableLen = rec.subMax.filter((s) => s > 0).length * ds, stableEnd = rec.sub.filter((s) => s > 0).length * ds;
  let expoH = 0, degH = 0; for (let k = 0; k < ser.t.length; k++) { const dtH = ser.t[k] - (k ? ser.t[k - 1] : 0); if (ser.sub[k] > 0) { expoH += dtH; degH += ser.sub[k] * dtH; } }
  const depProd = Sun.aq.inhWtCell.map((w) => aqDepression(Sun, p.inhEff, w)); // steady production and the cooldown clock run with the inhibition of the produced stream, not with the preservation dose
  const teqSteady = base.P.map((P, i) => S.hT0(P) - depProd[i]), subSteady = base.T.map((T, i) => teqSteady[i] - T);
  const tauC = S.C.map((C) => C / (S.U * PI * D0)), coolT = base.T.map((T, i) => (T <= teqSteady[i] ? 0 : S.tAmb[i] >= teqSteady[i] ? Infinity : tauC[i] * Math.log((T - S.tAmb[i]) / (teqSteady[i] - S.tAmb[i]))) / 3600), cooldown = Math.min(...coolT), iCool = coolT.indexOf(cooldown);
  const peakPhi = pk.phi, peakBlkHyd = pk.blk, peakVisc = Math.max(...ser.visc, 1);
  const led = h.ledger, waterUsed = led.formed * p.wfH, gasUsed = led.formed * (1 - p.wfH), peakRate = Math.max(...ser.rate, 0), convMax = Math.max(...rec.X, pk.phiX ? Math.max(...pk.phiX.map((f, i) => (f * p.rhoH * p.wfH) / Math.max(S.props(h.P[i], h.T[i]).rhoW * S.props(h.P[i], h.T[i]).wcut, 1e-9))) : 0);
  // ---- combined deposit profile at the worst time of the event
  const dHyd = pk.dHyd || h.dHyd, total = dHyd.map((d, i) => d + slow.dWax[i] + slow.dScale[i]), Deff = total.map((d) => Math.max(D0 - 2 * d, 0.03 * D0)), iMin = argmax(total), blockage = 1 - (Deff[iMin] / D0) ** 2;
  const roughEff = Math.max(...total.map((_, i) => roughnessUpdate(S.rough0, dHyd[i], slow.dWax[i], slow.dScale[i], D0))), grain = v.plugGrainUm > 0 ? v.plugGrainUm * 1e-6 : p.dPrim;
  const porPk = pk.por || h.por, perm = porPk.map((e) => kozenyCarman(e, grain));
  // clean versus fouled pressure drop in production (kernel gradient with the restricted bore and the rougher wall)
  let dpClean = 0, dpFoul = 0; for (let i = 0; i < n; i++) { const g0 = S.grad(i, pb.P[i], pb.T[i], D0, S.rough0, 1, fProd), g1 = S.grad(i, pb.P[i], pb.T[i], Math.max(D0 - 2 * (slow.dWax[i] + slow.dScale[i]), 0.2 * D0), roughBase[i], 1, fProd); dpClean += (g0.dpdx * ds) / 1e5; dpFoul += (g1.dpdx * ds) / 1e5; }
  const pIn0 = S.prof(fEvent).pIn, pInPeak = Math.max(...ser.pIn.filter((_, k) => ser.frac[k] > 0), pIn0), dpRise = pInPeak - pIn0;
  // slurry transport at the peak: settling of the agglomerates (hindered) and the velocity needed to keep them suspended
  const dAg = Math.max(pk.d43 ? pk.d43[pk.iPhi] : 0, p.dPrim), setH = settlingVelocity(dAg, p.rhoH, pk.rhoC || 800, pk.muC || 3e-3, { phi: pk.phiE ? Math.min(pk.phiE[pk.iPhi], 0.6) : 0 }), critH = sandCriticalVelocity({ d: dAg, D: D0, rhoP: Math.max(p.rhoH, (pk.rhoC || 800) + 1), rhoF: pk.rhoC || 800, mu: pk.muC || 3e-3, C: Math.max(peakPhi, 1e-6) }).oroskarTurian;
  const numPk = pk.N ? pk.N.reduce((a, b) => a + b, 0) / Math.max(pk.VL || 1, 1e-12) : 0;
  // flow that still leaks through a fully formed porous plug at the available pressure difference (Darcy–Forchheimer, Kozeny–Carman)
  const iPl = h.plug ? h.plug.i : iSub, prPlug = S.props(h.P[iPl], S.tAmb[iPl]), dpAvail = Math.max(p.pInMax - S.pOut, 0) * 1e5, vLeak = porousVelocity(dpAvail / v.plugLength, prPlug.muL, prPlug.rhoL, v.porosityInf, grain), pgr = porousGradient(vLeak, prPlug.muL, prPlug.rhoL, v.porosityInf, grain);
  // ---- plugging indicator
  const plug = h.plug, pOn = h.onset ? 1 : 1 - Math.exp(-Math.max(...rec.lam, 0)), sev = Math.max(peakBlkHyd / p.plugBlock, Math.log10(peakVisc) / 3, dpRise / Math.max(p.pInMax - pIn0, 1));
  const risk = plug ? 1 : clamp(Math.max(pOn * Math.max(maxSub > 0 ? 0.15 : 0, sev), mc.prob), 0, 1), plugT = mc.plugTimes.slice().sort((a, b) => a - b);
  // ---- inhibition
  const inhId = S.aq.inhId !== 'none' ? S.aq.inhId : 'MEG', inh = INHIBITORS[inhId], needDep = maxSub + (depNowC[iSub] - depSalt) + v.marginC, wReqIdeal = needDep > 0 ? inhibitorFor(depSalt + needDep / Math.max(p.inhEff, 0.1), inhId, S.aq.S) : 0;
  const prE = S.props(base.P[n - 1], base.T[n - 1]), mWater = prE.mW * fEvent, wR = Math.min(wReqIdeal, 94) / 100, mInh = (wR / (1 - wR)) * mWater, inhRate = (mInh / inh.rho) * 86400, qGasStd = (S.fm?.rates?.qGasStd ?? 0) * fEvent;
  const lossVap = inhId === 'MeOH' ? (v.meohVapK * wReqIdeal * qGasStd) / 1e6 : 0, lossOil = inhId === 'MeOH' ? v.meohOilK * wR * prE.mO * fEvent * 86400 : 0, lossPct = mInh > 0 ? (100 * (lossVap + lossOil)) / (mInh * 86400) : 0;
  const khiOk = maxSub <= v.khiLimit, khiHold = maxSub > 0 ? 48 * 2 ** ((v.khiLimit - maxSub) / 1.5) : Infinity, aaOk = prE.wcut * 100 <= v.aaWcLimit && prE.qO > 0 && p.mode !== 'gas';
  // ---- preservation: inhibitor needed to dose the water held in the line, power of the heating option
  const waterVol = base.P.reduce((a, P, i) => a + A0 * base.holdup[i] * S.props(P, base.T[i]).wcut * ds, 0), wPz = doseUsed / 100, presVol = plan.active && plan.id !== 'heat' ? (waterVol * 1000 * wPz) / Math.max(1 - wPz, 0.05) / S.aq.inh.rho : 0;
  const presName = { inhibit: `${S.aq.inh.name} placed in the line before it cools`, depressurise: `depressurisation to ${v.depressP} bara and restart with ${S.aq.inh.name}`, heat: 'heating that holds the line above its hydrate temperature', none: 'none' }[plan.id];
  // ---- remediation at the plug (or at the coldest point when there is no plug)
  const iR = plug ? plug.i : iSub, Rb = D0 / 2, pLoc = h.P[iR], tA = S.tAmb[iR], tdTwo = S.hT0(v.depressP) - depNow, melt2 = plugMeltTime({ R: Rb, U: S.U, Tamb: tA, Td: tdTwo, eps: v.porosityInf }), tdUp = S.hT0(pLoc) - depNow;
  const alphaH = HYDRATE.k / (p.rhoH * HYDRATE.cp), ste = (HYDRATE.cp * Math.max(tA - tdTwo, 0)) / HYDRATE.latent, lamS = ste > 0 ? stefanLambda(ste) : 0, melt1 = lamS > 0 ? v.plugLength ** 2 / (4 * lamS * lamS * alphaH) : Infinity;
  const projV = Math.sqrt((2 * Math.max(pLoc - v.depressP, 0) * 1e5 * 100) / (p.rhoH * (1 - v.porosityInf) * v.plugLength)), mPlug = p.rhoH * (1 - v.porosityInf) * A0 * v.plugLength;
  const wEq = inhibitorFor(depSalt + Math.max(tdUp + depNow - depSalt - tA, 0) + 1, 'MeOH', S.aq.S) / 100, meohMelt = (mPlug * p.wfH * wEq) / Math.max(1 - wEq, 0.05) / INHIBITORS.MeOH.rho, tMeoh = meohMelt > 0 ? (meohMelt * INHIBITORS.MeOH.rho) / (INHIBITORS.MeOH.rho * A0 * 1e-5 * Math.max(1 - wEq, 0.05)) : 0;
  const kdis = hydrateDissociationRate({ TK: tA + KEL, P: v.depressP, Peq: S.peq(tA + depNow), zG: 0.95 }), tKin = kdis > 0 ? (p.rhoH * grain) / (6 * (p.mwG + p.hydN * MW_W) * kdis) : Infinity;
  const tHold = Math.max(...teqSteady) + v.marginC, dehW = base.T.reduce((a, _, i) => a + Math.max(S.U * PI * D0 * (tHold - S.tAmb[i]) * ds, 0), 0), meltKW = (mPlug * HYDRATE.latent) / (24 * 3600) / 1000;
  // ---- sand
  const sandQ = v.sandRate / 86400 / v.sandRho, sand = pb.P.map((P, i) => { const c = slow.cells[i], vsl = c.vL * pb.holdup[i], C = sandQ / Math.max(c.pr.qL * fProd, 1e-12), cr = sandCriticalVelocity({ d: v.sandUm * 1e-6, D: D0, rhoP: v.sandRho, rhoF: c.pr.rhoL, mu: c.pr.muL, C, vsl, vm: c.vm, coef: v.sandCoef, shape: v.sandShape }), flat = Math.abs(S.theta[i]) < 0.5, st = settlingVelocity(v.sandUm * 1e-6, v.sandRho, c.pr.rhoL, c.pr.muL, { phi: C, shape: v.sandShape });
    return { C, cr, vm: c.vm, vL: c.vL, flat, bed: v.sandRate > 0 && (flat ? c.vm < cr.governing : c.vL < 3 * Math.abs(st.v)), st, hold: (C * c.vL) / Math.max(c.vL - Math.abs(st.vHindered), 0.05 * c.vL) }; });
  const flatS = sand.filter((s) => s.flat), sandCrit = flatS.length ? Math.max(...flatS.map((s) => s.cr.governing)) : Math.max(...sand.map((s) => s.cr.governing)), sandBed = sand.some((s) => s.bed), bedLen = sand.filter((s) => s.bed).length * ds;
  const sandMargin = Math.min(...sand.map((s) => (s.flat ? s.vm / s.cr.governing : 9))), sandMinRate = sandMargin > 0 ? (100 * fProd) / sandMargin : 100, sandInv = sand.reduce((a, s, i) => a + s.hold * v.sandRho * A0 * pb.holdup[i] * ds, 0), iV = argmax(sand.map((s) => s.vm)), erosion = sandErosionScreen(v.sandRate, sand[iV].vm, v.sandUm, D0 * 1000, pb.rhoM[iV], v.erosionSm);
  const relax = particleRelaxation({ d: v.sandUm * 1e-6, rhoP: v.sandRho, rhoF: slow.cells[0].pr.rhoL, mu: slow.cells[0].pr.muL });
  // ---- asphaltene
  const spec = S.fm?.spec || {}, pRes = num(spec.Pres, 300), tRes = num(spec.Tres, 90), tArr = base.T[n - 1], aPath = linspace(pRes, Math.max(S.pOut, 2), 40).map((P, k) => { const T = tRes + ((tArr - tRes) * k) / 39, o = S.fm ? S.fm.at(P, T) : S.props(P, T); return { P, T, wG: o.wG ?? 0, rho: o.rhoO, mw: o.mwO || 150 }; });
  let pBub = aPath[aPath.length - 1].P; for (const q of aPath) if (q.wG > 1e-4) { pBub = q.P; break; }
  const sara = { sat: v.saraSat, aro: v.saraAro, res: v.saraRes, asp: v.saraAsp }, cii = colloidalInstability(sara), db = deBoer(aPath[0].rho, Math.max(pRes - pBub, 0)), aspPhi = ((v.saraAsp / 100) * aPath[0].rho) / 1200;
  // Flory–Huggins onset on an isothermal depletion at reservoir temperature: density correlation and equation-of-state versions
  const isoPath = linspace(pRes, Math.max(S.pOut, 2), 30).map((P) => { const o = S.fm ? S.fm.at(P, tRes) : S.props(P, tRes); return { P, rho: o.rhoO, mw: o.mwO || 150, wG: o.wG ?? 0 }; });
  const aDen = asphalteneOnset({ pRes, tRes, pEnd: Math.max(S.pOut, 2), n: 30, phiA: aspPhi, deltaA: v.asphDelta, vA: v.asphMV, method: 'density', refOnset: v.asphOnsetRef, path: isoPath });
  const aEos = S.fm && oilLike(spec) && (v.asphMethod === 'eos' || v.asphOnsetRef > 0) ? asphalteneOnset({ spec, pRes, tRes, pEnd: Math.max(S.pOut, 2), n: 14, phiA: aspPhi, deltaA: v.asphDelta, vA: v.asphMV, method: 'eos', refOnset: v.asphOnsetRef, path: null }) : null, aUse = v.asphMethod === 'eos' && aEos ? aEos : aDen;
  const iFh = aUse.precipitates ? 0 : -1, fhMin = aUse.minSolubility, fhOnset = aUse.upperOnset;
  const aScore = (db.cls === 'severe problems' ? 2 : db.cls === 'slight problems' ? 1 : 0) + (cii.cls === 'unstable' ? 2 : cii.cls === 'uncertain' ? 1 : 0) + (iFh >= 0 ? 2 : 0), oilSys = p.mode !== 'gas', cls = (x) => (x >= 4 ? 'high' : x >= 2 ? 'medium' : 'low'), aRiskRaw = !oilSys ? 'low' : cls(aScore), aRisk = !oilSys ? 'low' : cls(v.aspInhOn ? aScore * (1 - v.aspInhEff / 100) : aScore);
  // ---- scale summary and seawater mixing curve at arrival conditions
  const iSI = argmax(slow.si), scaleSI = slow.si[iSI], scaleMineral = scaleSI > 0 ? slow.mineral[iSI] : 'none', mixF = linspace(0, 1, 11), w0 = waterOf(v), tMix = pb.T[0], pMix = pb.P[0];
  const mix = mixF.map((f) => scaleIndices(mixWaters(w0, slow.sea, f), tMix, pMix, { yCO2: v.co2Pct / 100, model: v.actModel })), mixRaw = v.srpOn ? mixF.map((f) => scaleIndices(mixWaters(w0, SEAWATER, f), tMix, pMix, { yCO2: v.co2Pct / 100, model: v.actModel })) : mix, sIn = slow.sis[0];
  // scale management: every supersaturated mineral against its induction time with the inhibitor dose; carbonates and calcium sulphates carry no inhibition model and are judged on the index alone
  const scaleCtl = slow.control.filter((c) => c.q), tProt = v.scaleProtectH * 3600, doseS = v.scaleInhOn ? v.scaleInhDose : 0, riskOf = (c) => (c.q.model ? c.q.risk : c.SI < 0.3 ? 0.3 * (c.SI / 0.3) : clamp(c.SI, 0.3, 1));
  const scaleRisk = scaleCtl.length ? Math.max(...scaleCtl.map(riskOf)) : 0, scaleRiskRaw = scaleCtl.length ? Math.max(...scaleCtl.map((c) => (c.q.model ? c.raw.risk : riskOf(c)))) : 0, gov = scaleCtl.length ? scaleCtl.reduce((a, c) => (riskOf(c) > riskOf(a) ? c : a)) : null;
  const micS = scaleCtl.length ? Math.max(...scaleCtl.map((c) => (c.q.model && Number.isFinite(c.q.mic) ? c.q.mic : 0))) : 0, scaleManaged = scaleRisk < 1 && scaleCtl.every((c) => (c.q.model ? c.q.tInh >= tProt : c.SI < 0.3)), siRawSea = v.srpOn && v.swFrac > 0 ? Math.max(...scaleIndices(mixWaters(w0, SEAWATER, clamp(v.swFrac / 100, 0, 1)), pb.T[iSI], pb.P[iSI], { yCO2: v.co2Pct / 100, model: v.actModel }).minerals.map((q) => q.SI)) : scaleSI;
  const hrs = (t) => (t >= 3.6e12 ? 'more than 100 000 years' : t >= 3.156e8 ? `${(t / 3.156e7).toPrecision(2)} years` : t >= 172800 ? `${(t / 86400).toPrecision(3)} d` : `${(t / 3600).toPrecision(3)} h`);
  // ---- sectional versus quadrature method of moments in the worst cell (short batch: aggregation + breakage)
  const qm = (() => {
    const phi = Math.max(peakPhi, 0.02), env = { shear: pk.gdot || 200, eps: 0, nu: (pk.muC || 3e-3) / (pk.rhoC || 800), mu: pk.muC || 3e-3, TK: (pk.T ?? 4) + KEL, dRho: Math.abs(p.rhoH - (pk.rhoC || 800)), alpha: p.aggEff }, g2 = pbeGrid(40, p.dPrim / 2, 60 * p.dPrim), kP = 5, L0 = g2.L[kP], N0 = phi / g2.v[kP];
    const b0 = aggregationKernel(L0, L0, env).total, tB = b0 * N0 > 0 ? 3 / (b0 * N0) : 1, be = new Float64Array(g2.n * g2.n); for (let i = 0; i < g2.n; i++) for (let j = 0; j < g2.n; j++) be[i * g2.n + j] = aggregationKernel(g2.L[i], g2.L[j], env).total;
    const Ni = new Float64Array(g2.n); Ni[kP] = N0; const sct = solvePBE(g2, Ni, tB, { beta: be, maxSub: 4000, frac: 0.05 }), ms = pbeMoments(g2, sct.N);
    const sg = 0.05, mom = [0, 1, 2, 3, 4, 5].map((k) => N0 * L0 ** k * Math.exp((k * k * sg * sg) / 2)), q = solveQMOM(mom, tB, { beta: (a, b) => aggregationKernel(a, b, env).total, steps: 50 });
    const mcr = v.mcParticles > 0 ? solvePBEMonteCarlo(new Array(v.mcParticles).fill(L0), tB, { beta: (a, b) => aggregationKernel(a, b, env).total, conc: N0, seed: v.seed }) : null;
    return { tB, m0s: ms.m0 / N0, m0q: q.m[0] / N0, m3s: ms.m3 / (N0 * L0 ** 3), m3q: q.m[3] / mom[3], d43s: ms.d43, d43q: q.m[4] / q.m[3], ok: q.ok, nodes: q.nodes, mc: mcr ? { m0: mcr.m0 / N0, vol: mcr.vol / (N0 * g2.v[kP]), d43: mcr.d43, events: mcr.events, n: mcr.n } : null };
  })();
  // ---- population balance on the resolved flow field of the pipe section with the most hydrate (of the design case; of the unprotected event when the preserved line forms none)
  const cfd = v.cfdOn ? (() => {
    const src = pk.N && peakPhi > 1e-6 ? { k: pk, what: 'design case' } : unprot && unprot.pk.N && unprot.peakPhi > 1e-6 ? { k: unprot.pk, what: 'unprotected event' } : { k: null, what: 'nominal population (2 vol % of primary particles)' }, k = src.k;
    const cI = slow.cells[k ? k.iPhi : iSub], muC0 = k?.muC || (oil ? cI.pr.muO : cI.pr.muL), rhoC = k?.rhoC || (oil ? cI.pr.rhoO : cI.pr.rhoL), TK = (k?.T ?? S.tAmb[iSub]) + KEL, U = Math.max(k?.vL ?? cI.vL * fEvent, 0.05), Dp = D0;
    const Nin = new Array(grid.n).fill(0); let phi = 0;
    if (k) { for (let c = 0; c < grid.n; c++) { Nin[c] = k.N[c] / Math.max(k.VL, 1e-12); phi += Nin[c] * grid.v[c]; } } else { const kP = grid.L.reduce((b, L, j) => (Math.abs(Math.log(L / p.dPrim)) < Math.abs(Math.log(grid.L[b] / p.dPrim)) ? j : b), 0); phi = 0.02; Nin[kP] = phi / grid.v[kP]; }
    const muS = muC0 * clamp(k?.visc ? k.visc[k.iPhi] : slurryViscosity(phi, p.viscModel, { phiMax: p.phiMax }), 1, 50), F = pipeFlowField({ D: Dp, U, rho: rhoC, mu: muS, nr: v.cfdNr });
    const dAof = (shear) => maxAgglomerateSize({ dp: p.dPrim, Fa: p.cohesion * p.dPrim, mu0: muC0, shear, phi, phiMax: p.phiMax, fr: p.fr }).dA;
    const r = solvePBEField({ field: F, grid, Nin, L: v.cfdLenD * Dp, nx: v.cfdNx, kin: { alpha: p.aggEff, kBreak: p.kBreak, dAof, dPrim: p.dPrim, mu: muC0, rhoF: rhoC, rhoP: p.rhoH, TK }, wall: { stick: p.adhesion, adhForce: p.adhForce }, maxSub: 8 });
    const last = r.d43.length - 1, dep1d = r.vdBulk * r.mixed.vol, gk = S.grad(k ? k.iPhi : iSub, base.P[k ? k.iPhi : iSub], base.T[k ? k.iPhi : iSub], D0, 0, 1, fEvent);
    return { F, r, what: src.what, phi, U, mu: muS, rho: rhoC, L: v.cfdLenD * Dp, dCore: r.d43[last][0], dWall: r.d43[last][F.nr - 1], dMean: r.mixed.d43, d1d: k?.d43 ? k.d43[k.iPhi] : 0, depRatio: dep1d > 0 ? r.dep[r.dep.length - 1] / dep1d : null, phiCore: r.phi[last][0], phiWall: r.phi[last][F.nr - 1], tauKernel: gk.tauW };
  })() : null;
  const vdw = vdwpMethane(Math.max(mean(S.tAmb), 0.5)), occ = langmuirOccupancy(Math.max(mean(S.tAmb), 0) + KEL, methaneFugacity(mean(base.P), Math.max(mean(S.tAmb), 0) + KEL));

  // ---- threats ranked by zone
  const nz = Math.min(8, n), zones = Array.from({ length: nz }, (_, zI) => {
    const a = Math.floor((zI * n) / nz), b = Math.max(Math.floor(((zI + 1) * n) / nz), a + 1), ix = Array.from({ length: b - a }, (_, k) => a + k), mx = (f) => Math.max(...ix.map(f));
    const sc = { Hydrate: clamp(Math.max(mx((i) => rec.subMax[i]) / 10, mx((i) => dHyd[i] / (0.15 * D0)), mx((i) => (pk.phiX ? pk.phiX[i] : 0) / 0.2)), 0, 1.5) * (mx((i) => rec.subMax[i]) > 0 ? 1 : 0), Wax: clamp(mx((i) => slow.dWax[i] * 1000) / v.waxLimitMm, 0, 1.5), Scale: clamp(Math.max(mx((i) => slow.si[i]) / 1.5, mx((i) => slow.scaleRate[i]) / 2), 0, 1.5), Sand: clamp(mx((i) => (v.sandRate > 0 ? sand[i].cr.governing / Math.max(sand[i].vm, 1e-6) : 0)) - 0.3, 0, 1.5) };
    const top = Object.entries(sc).sort((p1, p2) => p2[1] - p1[1])[0];
    return { from: (a * ds) / 1000, to: (b * ds) / 1000, sc, top: top[1] > 0.05 ? top[0] : 'None' };
  });

  // ---- additional solvers: statistical thermodynamics, slug coupling, stochastic population balance, particle tracking, Eulerian solids transport, inhibitor contact
  await tick(0.95, 'Particle tracking, solids transport and cross-checks');
  const X = {};
  { // Kihara-potential Langmuir constants and Gibbs-energy minimisation at the coldest point of the event
    const tC = Math.max(mean(S.tAmb), 0.5), TKc = tC + KEL, pC = clamp(base.P[iSub], 5, 800), prC = S.props(pC, tC), wI = clamp(S.aq.inhWtCell[iSub] / 100, 0, 0.95), nW = Math.max(prC.mW / MW_W, 1e-9), nG = Math.max(prC.mG / (clamp(prC.mwG || 18, 16, 60) * 1e-3), 1e-9);
    const solPerKg = (2 * Math.min(S.aq.S, 260)) / 58.44 + ((wI / (1 - wI)) * 1000) / (S.aq.inh.MW || 62.07), nS = nW * MW_W * solPerKg;
    X.th = { tC, pC, pp: vdw, kh: vdwpMethane(tC, { langmuir: 'kihara' }), cK: { s: kiharaLangmuir(TKc, HYD_REF.cages.small), l: kiharaLangmuir(TKc, HYD_REF.cages.large) }, cP: langmuirOccupancy(TKc, 1), gm: hydrateGibbsMin({ Tc: tC, P: pC, nW, nGas: nG, nSolute: nS, langmuir: v.langmuirModel === 'kihara' ? 'kihara' : 'parrish' }), nW, nG, nS };
  }
  { // slug unit cell: where the event flows in slug flow, how much the gas-to-liquid conductance exceeds the stratified value
    const pe = S.prof(fEvent), ix = []; for (let i = 0; i < n; i++) if (/slug/i.test(pe.regime[i])) ix.push(i);
    const iS = ix.length ? ix.reduce((a, b) => (rec.subMax[b] > rec.subMax[a] ? b : a)) : -1, sm = iS >= 0 ? S.slugMT(iS, fEvent) : null, fs = ctx.outputs?.flow?.slug;
    X.slug = { n: ix.length, len: ix.length * ds, i: iS, sm, formedSlug: led.formedSlug, flowFreq: isNum(fs?.freq) ? fs.freq : null, flowType: typeof fs?.type === 'string' ? fs.type : null };
  }
  if (v.trackOn) { // Eulerian–Lagrangian: stochastic parcels in the flow field of the event
    const pe = S.prof(fEvent), sandT = v.trackWhat === 'sand', fld = { ds, D: Deff.slice(), vL: [], uStar: [], rhoF: [], mu: [], theta: S.theta.slice(), cold: rec.subMax.map((s) => !sandT && s > 0) };
    for (let i = 0; i < n; i++) { const pr = S.props(pe.P[i], pe.T[i]), A = (PI * Deff[i] ** 2) / 4; fld.vL.push(Math.max((pr.qL * fEvent) / A / pe.holdup[i], 1e-3)); fld.uStar.push(Math.sqrt(Math.max(pe.tauW[i], 1e-6) / pr.rhoL)); fld.rhoF.push(sandT ? pr.rhoL : oil ? pr.rhoO : pr.rhoL); fld.mu.push(sandT ? pr.muL : (oil ? pr.muO : pr.muL)); }
    let sizes = null; if (!sandT && pk.N) { const vt = pk.N.reduce((a, N, c) => a + N * grid.v[c], 0); if (vt > 0) { sizes = []; pk.N.forEach((N, c) => { const k = Math.round((200 * N * grid.v[c]) / vt); for (let q = 0; q < k; q++) sizes.push(grid.L[c]); }); if (!sizes.length) sizes = null; } }
    X.track = trackParticles(fld, { n: v.nParcels, sizes, d: sandT ? v.sandUm * 1e-6 : p.dPrim, rhoP: sandT ? v.sandRho : p.rhoH, adhesion: sandT ? 0 : p.adhesion, equip: mx.traps.map((t) => ({ x: t.x, eff: t.eff, name: t.name })), seed: v.seed, stepFrac: v.trackStepPct / 100, TK: mean(pe.T) + KEL });
    X.track.what = sandT ? 'sand' : 'hydrate';
  }
  { // Eulerian–Eulerian sand transport: advection, axial dispersion, settling into a bed and re-entrainment over the production period
    const dS = v.sandUm * 1e-6, f = { ds, A: [], u: [], W: [], Dax: [], ws: [], pd: [], E: [] }, tauC = [];
    for (let i = 0; i < n; i++) { const c = slow.cells[i], H = pb.holdup[i], ct = Math.cos(S.theta[i]), st = Math.sin(S.theta[i]), wh = Math.abs(sand[i].st.vHindered), tc = (v.sandShields > 0 ? v.sandShields : shieldsCritical(dS, v.sandRho, c.pr.rhoL, c.pr.muL)) * Math.max(v.sandRho - c.pr.rhoL, 1) * G * dS, tw = Math.max(pb.tauW[i], 0);
      f.A.push(A0 * H); f.u.push(Math.max(c.vL - wh * Math.max(st, 0), 0.02 * c.vL, 1e-4)); f.W.push(D0 * Math.sqrt(Math.max(1 - (2 * Math.min(H, 0.5) - 1) ** 2, 0.05)) * 0.5); f.Dax.push(10.1 * (D0 / 2) * Math.sqrt(tw / c.pr.rhoL)); f.ws.push(wh * Math.abs(ct)); f.pd.push(clamp(1 - tw / tc, 0, 1)); f.E.push(v.sandEntrain * 1e-4 * Math.max(tw / tc - 1, 0)); tauC.push(tc); }
    const tEnd = v.sandHours * 3600, ee = v.sandRate > 0 ? solidsTransport(f, { cIn: sand[0].C, tEnd, dt: tEnd / v.eeSteps, trap: mx.traps.map((t) => ({ i: t.i, eff: t.eff })) }) : null;
    const bedH = ee ? ee.bed.map((b) => Math.min(Math.cbrt((9 * (b / 0.6) ** 2) / (16 * D0)), D0)) : new Array(n).fill(0); // height of a circular segment of area a ≈ (4/3) h √(D h), bed at 40 % porosity
    X.ee = { ee, bedH, tEnd, tauC, bedMass: ee ? ee.ledger.bed * v.sandRho : 0, trapMass: ee ? ee.ledger.trapped * v.sandRho : 0, iMax: argmax(bedH) };
  }
  { // mass-transfer-controlled dissociation: inhibitor reaching the plug face (film, stagnant diffusion, diffusion with reaction in the pores)
    const o = { wBulk: 0.98, wEq: clamp(wEq, 0.01, 0.9), D: v.inhDiff * 1e-9, rhoL: INHIBITORS.MeOH.rho, rhoH: p.rhoH, eps: v.porosityInf, wfH: p.wfH, kFilm: v.inhFilmK * 1e-5, dGrain: grain, Lp: v.plugLength };
    X.mt = { film: inhibitorDissociation({ ...o, mode: 'film' }), stag: inhibitorDissociation({ ...o, mode: 'stagnant', tEnd: 30 * 86400 }), por: inhibitorDissociation({ ...o, mode: 'porous' }), wEq: o.wEq };
  }
  { // measurements entered for comparison
    const rows = (Array.isArray(v.expData) ? v.expData : []).filter((r) => r && String(r.type || '').trim() && isNum(+r.meas)), env = { run: null, spec, pCO2: null };
    X.exp = rows.map((r) => ({ type: String(r.type).trim().toLowerCase(), a: num(r.a, 0), b: num(r.b, 0), c: num(r.c, 0), meas: +r.meas, env }));
  }

  // ---- outputs for the other suites
  const tOnset = h.onset ? h.onset.t / 3600 : null, piggingInterval = slow.piggingInterval;
  const outputs = {
    hydrateRisk: r3(risk), maxSubcooling: r3(maxSub, 4), onsetX: h.onset ? r3(h.onset.x, 5) : null, onsetTime: tOnset === null ? null : r3(tOnset, 4), hydrateFraction: r3(peakPhi, 4), hydrateRate: r3(peakRate, 4),
    depositProfile: { x: S.x.map((a) => +a.toFixed(1)), hydrate: dHyd.map((d) => +d.toExponential(4)), wax: slow.dWax.map((d) => +d.toExponential(4)), scale: slow.dScale.map((d) => +d.toExponential(4)), total: total.map((d) => +d.toExponential(4)) },
    effectiveId: r3(Deff[iMin], 5), roughnessEff: r3(roughEff, 4), blockage: r3(blockage, 4), plugTime: plug ? r3(plug.t / 3600, 4) : null, plugX: plug ? r3(plug.x, 5) : null, plugProbability: r3(plug && mc.n === 0 ? 1 : mc.prob),
    inhibitorRequired: r3(wReqIdeal, 4), inhibitorRate: r3(inhRate, 4), wat: r3(slow.wat, 4), waxRate: r3(slow.waxRate0, 4), waxMass: r3(slow.waxMass, 4), piggingInterval: piggingInterval === null ? null : r3(piggingInterval, 4),
    scaleSI: r3(scaleSI, 4), scaleMineral, scaleManaged, scaleInhibitorDose: r3(doseS, 4), residualScaleRisk: scaleRisk < 1e-6 ? 0 : r3(scaleRisk, 4), scaleRiskUntreated: scaleRiskRaw < 1e-6 ? 0 : r3(scaleRiskRaw, 4), scaleMIC: r3(micS, 4), scaleInductionHours: gov && gov.q.model ? r3(Math.min(gov.q.tInh / 3600, 1e12), 4) : null, scaleSIUnmanaged: r3(siRawSea, 4), sulphateRemoval: !!v.srpOn, asphalteneRisk: aRisk, asphalteneRiskUnmitigated: aRiskRaw, asphalteneInhibitor: !!v.aspInhOn, waxInhibitor: !!v.waxInhOn, sandCriticalVelocity: r3(sandCrit, 4), sandBed, slurryViscosityFactor: r3(peakVisc, 4),
    hydrateRiskUnprotected: r3(unprot ? unprot.risk : risk), plugProbabilityUnprotected: r3(unprot ? unprot.prob : plug && mc.n === 0 ? 1 : mc.prob), maxSubcoolingUnprotected: r3(unprot ? unprot.maxSub : maxSub, 4), hydrateFractionUnprotected: r3(unprot ? unprot.peakPhi : peakPhi, 4),
    preservation: plan.id, preservationDose: r3(doseUsed, 4), preservationVolume: r3(presVol, 4), preservationPower: r3((mx.holdPower || 0) / 1000, 4), preservationLeadTime: plan.active ? plan.lead : null,
    // extras
    scenario: v.scenario, hydrateMode: p.mode, hydrateLength: r3(stableLen, 5), exposureHours: r3(expoH, 4), cooldownTime: Number.isFinite(cooldown) ? r3(cooldown, 4) : null, hydrateMass: r3(Math.max(...ser.susp.map((s, k) => s + ser.dep[k]), 0), 4), waterConversion: r3(convMax, 4), gasConsumed: r3(gasUsed, 4),
    plugProbabilityInterval: [r3(mc.lo), r3(mc.hi)], plugTimeP10: plugT.length ? r3(quantile(plugT, 0.1), 4) : null, plugTimeP50: plugT.length ? r3(quantile(plugT, 0.5), 4) : null, plugTimeP90: plugT.length ? r3(quantile(plugT, 0.9), 4) : null, plugMechanism: plug ? plug.mech : null,
    dpIncrease: r3(dpRise + (dpFoul - dpClean), 4), dpFouling: r3(dpFoul - dpClean, 4), inhibitor: inhId, inhibitorDepression: r3(depNow - depSalt, 4), meltTimeTwoSided: Number.isFinite(melt2.tNumeric) ? r3(melt2.tNumeric / 3600, 4) : null, heatingPower: r3(dehW / 1000, 4),
    waxRestartPressure: r3(slow.restartDp, 4), sandRate: r3(v.sandRate / 86400, 4), sandSize: v.sandUm * 1e-6, sandDensity: v.sandRho, sandMinRateFraction: r3(sandMinRate / 100, 4), sandErosionScreen: r3(erosion, 3), plugXP50: mc.plugX.length ? r3(quantile(mc.plugX, 0.5), 5) : null, particleNumber: r3(numPk, 4), slurryCriticalVelocity: r3(critH, 4), depositPermeability: r3(kozenyCarman(v.porosityInf, grain), 3), captureRate: r3(rec.captureMax, 4), removalRate: r3(rec.removalMax, 4), asphalteneCII: r3(cii.cii, 4), particleSize: r3(pk.d43 ? pk.d43[pk.iPhi] : 0, 4),
  };

  // ---- KPIs
  const st = (x, w, b) => (x >= b ? 'bad' : x >= w ? 'warn' : 'ok');
  const kpis = [
    { label: 'Peak subcooling in the event', value: r3(maxSub, 3), unit: '°C', status: st(maxSub, 0.01, 6), help: `Hydrate temperature minus fluid temperature, largest value along the line and through the ${SC_LABEL[v.scenario].toLowerCase()} (positive = inside the hydrate region).` },
    { label: 'Length inside the hydrate region', value: r3(stableLen / 1000, 3), unit: 'km', status: st(stableLen, 1, 0.25 * S.L), help: 'Pipe length that was inside the hydrate region at any time of the event.' },
    { label: 'Hydrate onset', value: tOnset === null ? 'none' : r3(tOnset, 3), unit: tOnset === null ? '' : 'h', status: tOnset === null ? 'ok' : 'bad', help: 'Time from the start of the simulated sequence to the first nucleation (hazard integral reaches one).' },
    { label: 'Peak hydrate fraction of the liquid', value: r3(peakPhi * 100, 3), unit: 'vol %', status: st(peakPhi, 0.01, 0.1) },
    { label: 'Peak slurry viscosity factor', value: r3(peakVisc, 3), unit: '×', status: st(peakVisc, 2, 20), help: 'Relative viscosity of the hydrate slurry from the effective (agglomerate) volume fraction.' },
    { label: 'Peak area blockage (all deposits)', value: r3(blockage * 100, 3), unit: '%', status: st(blockage, 0.1, 0.5) },
    { label: 'Plug', value: plug ? r3(plug.t / 3600, 3) : 'no plug', unit: plug ? `h at ${(plug.x / 1000).toFixed(1)} km` : '', status: plug ? 'bad' : 'ok', help: plug ? plug.mech : 'No plug in the base-case march; see the plugging probability for the effect of uncertain kinetics.' },
    { label: 'Plugging probability', value: r3(outputs.plugProbability * 100, 3), unit: '%', status: st(outputs.plugProbability, 0.05, 0.3), help: `Fraction of ${mc.n} Monte Carlo samples that plug within the event (90 % interval ${(mc.lo * 100).toFixed(0)}–${(mc.hi * 100).toFixed(0)} %).` },
    { label: 'Hydrate risk index', value: r3(risk, 3), unit: '0–1', status: st(risk, 0.2, 0.6) },
    ...(unprot ? [{ label: 'Hydrate risk index without preservation', value: r3(unprot.risk, 3), unit: '0–1', status: 'ok', help: `Comparison case: the same event with no preservation reaches ${unprot.maxSub.toFixed(1)} °C of subcooling and ${(unprot.peakPhi * 100).toFixed(1)} vol % hydrate${unprot.plug ? ` and plugs after ${(unprot.plug.t / 3600).toFixed(1)} h` : ''}; plugging probability ${(unprot.prob * 100).toFixed(0)} %. Shown for comparison; the design case is the preserved line.` }] : []),
    { label: `${inh.name} needed for a ${v.marginC} °C margin`, value: r3(wReqIdeal, 3), unit: 'wt %', status: wReqIdeal > S.aq.inhWt + 0.5 ? (wReqIdeal > 60 ? 'bad' : 'warn') : 'ok', help: `Injection rate ${tx(inhRate)} m³/d at the event rate.` },
    { label: 'Cooldown time to hydrate temperature', value: Number.isFinite(cooldown) ? r3(cooldown, 3) : 'never', unit: Number.isFinite(cooldown) ? 'h' : '', status: Number.isFinite(cooldown) ? st(-cooldown, -12, -4) : 'ok', help: 'Lumped exponential cooldown of the first station to reach its hydrate temperature after a shut-in from steady production.' },
    { label: 'Initial wax build-up rate', value: r3(slow.waxRate0, 3), unit: 'mm/d', status: st(slow.waxRate0, 0.02, 0.2) },
    { label: 'Pigging interval', value: piggingInterval === null ? 'not needed' : r3(piggingInterval, 3), unit: piggingInterval === null ? '' : 'd', status: piggingInterval === null ? 'ok' : st(-piggingInterval, -60, -10) },
    { label: 'Highest scale saturation index', value: r3(scaleSI, 3), unit: scaleSI > 0 ? scaleMineral : '', status: scaleSI <= 0 || scaleManaged ? 'ok' : scaleRisk >= 1 && scaleSI > 1 ? 'bad' : 'warn', help: scaleSI > 0 && gov ? (gov.q.model ? `Supersaturated, but nucleation needs ${hrs(gov.q.tInh)} ${doseS > 0 ? `with ${doseS} mg/L of ${v.scaleInh}` : 'untreated'} (${hrs(gov.raw.t0)} untreated) against ${v.scaleProtectH} h of required protection.` : 'No inhibition model for this mineral: judged on the index alone.') : 'No mineral is supersaturated.' },
    { label: 'Residual scale risk', value: scaleRisk < 1e-6 ? 0 : r3(scaleRisk, 3), unit: '0–1', status: st(scaleRisk, 0.3, 1), help: `Required protection time divided by the nucleation induction time of the governing mineral with the inhibitor in place (1 = scale forms inside the system). Untreated: ${scaleRiskRaw.toPrecision(2)}. Minimum inhibitor concentration ${micS.toPrecision(2)} mg/L.` },
    { label: 'Asphaltene risk', value: aRisk, unit: v.aspInhOn && oilSys ? `managed (screening: ${aRiskRaw})` : '', status: aRisk === 'high' ? 'bad' : aRisk === 'medium' ? 'warn' : 'ok', help: `de Boer “${db.cls}”, colloidal instability index ${cii.cii.toFixed(2)} (${cii.cls})${iFh >= 0 ? ', thermodynamic onset on the depletion path' : ', no thermodynamic onset on the depletion path'}${v.aspInhOn ? `; deposition tendency reduced by ${v.aspInhEff} % by the injected dispersant` : ''}.` },
    { label: 'Sand transport margin', value: r3(sandMargin, 3), unit: 'v / v_critical', status: v.sandRate > 0 ? st(-sandMargin, -1.3, -1) : 'ok' },
  ];

  // ---- warnings and recommendations
  const warnings = [], recs = [];
  if (plug) warnings.push({ level: 'bad', msg: `The base-case march plugs after ${(plug.t / 3600).toFixed(1)} h at ${(plug.x / 1000).toFixed(1)} km: ${plug.mech}.` });
  if (maxSub > 0 && !plug) warnings.push({ level: maxSub > 6 ? 'bad' : 'warn', msg: `The line is inside the hydrate region over ${(stableLen / 1000).toFixed(1)} km with up to ${maxSub.toFixed(1)} °C of subcooling for ${expoH.toFixed(1)} h.` });
  if (mc.prob >= 0.05) warnings.push({ level: mc.prob > 0.3 ? 'bad' : 'warn', msg: `Plugging probability ${(mc.prob * 100).toFixed(0)} % (90 % interval ${(mc.lo * 100).toFixed(0)}–${(mc.hi * 100).toFixed(0)} %, ${mc.n} samples of uncertain kinetics, cohesion, adhesion and nucleation).` });
  if (peakPhi > 0.005 && (pk.vL ?? 9) < critH) warnings.push({ level: 'warn', msg: `At the peak the liquid moves at ${(pk.vL ?? 0).toFixed(2)} m/s, below the ${critH.toFixed(2)} m/s needed to keep ${(dAg * 1e6).toFixed(0)} µm agglomerates suspended: a moving or stationary hydrate bed is likely.` });
  if (h.limited) warnings.push({ level: 'info', msg: 'Agglomeration and breakage are much faster than the transport step in part of the line; those cells were integrated with the implicit conservative scheme, which resolves the equilibrium size but not the sub-second transient.' });
  if (Math.max(...slow.dWax.map((d, i) => 1 - (Dbase[i] / D0) ** 2)) > 0.5) warnings.push({ level: 'bad', msg: `Wax and scale alone close more than half of the bore after ${v.depositDays} d without pigging.` });
  if (piggingInterval !== null && piggingInterval < 30) warnings.push({ level: 'warn', msg: `Wax reaches ${v.waxLimitMm} mm in ${piggingInterval.toFixed(0)} d at ${(slow.waxRate0).toFixed(3)} mm/d.` });
  if (scaleSI > 0 && !scaleManaged && gov) warnings.push({ level: scaleRisk >= 1 && gov.SI > 1 ? 'bad' : 'warn', msg: gov.q.model
    ? `${gov.name} (SI ${gov.SI.toFixed(2)} at ${xkm[gov.i].toFixed(1)} km) nucleates after ${hrs(gov.q.tInh)}${doseS > 0 ? ` with ${doseS} mg/L of ${v.scaleInh}` : ' without inhibitor'}, inside the ${v.scaleProtectH} h the brine spends in the system: the minimum inhibitor concentration is ${Number.isFinite(gov.q.mic) ? gov.q.mic.toPrecision(2) + ' mg/L' : 'not reachable with this product'}.`
    : `${gov.name} is supersaturated (SI ${gov.SI.toFixed(2)} at ${xkm[gov.i].toFixed(1)} km); no inhibition model is available for this mineral, so no credit is taken for the inhibitor.` });
  if (v.srpOn && v.swFrac > 0 && siRawSea > scaleSI + 0.05) warnings.push({ level: 'info', msg: `Sulphate removal (${v.srpSO4} mg/L in the injection water) holds the highest saturation index at ${scaleSI.toFixed(2)}; with raw seawater at ${v.swFrac} % breakthrough it would be ${siRawSea.toFixed(2)}.` });
  if (sIn.I > 1 && v.actModel !== 'pitzer') warnings.push({ level: 'warn', msg: `Ionic strength ${sIn.I.toFixed(2)} mol/kg is beyond the range of the ion-association models; use the Pitzer model.` });
  if (Math.abs(ionicStrength(slow.water).balance) > 0.1) warnings.push({ level: 'info', msg: `The water analysis has a charge imbalance of ${(100 * ionicStrength(slow.water).balance).toFixed(0)} %; check the sodium or chloride value.` });
  if (sandBed) warnings.push({ level: 'warn', msg: `Sand settles over ${(bedLen / 1000).toFixed(1)} km: the mixture velocity is below the critical velocity (${sandCrit.toFixed(2)} m/s).` });
  if (aRisk !== 'low') warnings.push({ level: aRisk === 'high' ? 'bad' : 'warn', msg: `Asphaltene screening: de Boer "${db.cls}", colloidal instability index ${cii.cii.toFixed(2)} (${cii.cls})${iFh >= 0 && fhOnset ? `, Flory–Huggins onset near ${fhOnset.toFixed(0)} bara` : ''}${v.aspInhOn ? `; a ${v.aspInhEff} % effective dispersant does not bring the risk down to low` : ''}.` });
  if (slow.restartDp > Math.max(p.pInMax - S.pOut, 1)) warnings.push({ level: 'bad', msg: `A gelled line would need ${slow.restartDp.toFixed(0)} bar to break the wax gel, more than the available ${(p.pInMax - S.pOut).toFixed(0)} bar.` });
  if (S.src !== 'flow suite') warnings.push({ level: 'info', msg: S.custom ? 'The flow picture was recomputed here with the shared kernel because the line or fluid inputs differ from the case.' : 'The flow suite has not been run: pressure, temperature and holdup come from the shared kernel estimate.' });
  if (plan.active && Number.isFinite(cooldown) && plan.lead > cooldown) warnings.push({ level: 'bad', msg: `The preservation needs ${plan.lead} h to put in place but the line reaches its hydrate temperature ${cooldown.toFixed(1)} h after the trip at ${xkm[iCool].toFixed(1)} km: hydrate can form before the line is protected.` });
  if (plan.active && doseCapped) warnings.push({ level: 'bad', msg: `The dose that would keep the shut-in line ${v.marginC} °C outside the hydrate region exceeds the practical limit of ${doseUsed.toFixed(0)} wt % ${S.aq.inh.name}; ${maxSub.toFixed(1)} °C of subcooling remains with the preservation in place.` });
  if (unprot && unprot.maxSub > 0) warnings.push({ level: 'info', msg: `Comparison case without preservation: ${unprot.maxSub.toFixed(1)} °C of subcooling for ${unprot.expo.toFixed(0)} h, hydrate up to ${(unprot.peakPhi * 100).toFixed(1)} vol % of the liquid, ${unprot.plug ? `a plug after ${(unprot.plug.t / 3600).toFixed(1)} h` : 'no plug in the base march'}, plugging probability ${(unprot.prob * 100).toFixed(0)} % (risk index ${unprot.risk.toFixed(2)}). The design case above is the line preserved by ${presName}.` });
  if (plan.active) recs.push(plan.id === 'heat' ? `Preservation by heating: ${((mx.holdPower || 0) / 1000).toFixed(0)} kW holds the line at ${(mx.holdT ?? 0).toFixed(1)} °C or warmer; the heating must be on within ${Number.isFinite(cooldown) ? cooldown.toFixed(1) : 'the cooldown'} h of the trip.`
    : `Preservation: bring the water in the line to ${doseUsed.toFixed(0)} wt % ${S.aq.inh.name} (about ${presVol.toFixed(1)} m³ of inhibitor for ${waterVol.toFixed(0)} m³ of water held in the line${plan.id === 'depressurise' ? `, and depressurise to ${v.depressP} bara` : ''}) within ${Number.isFinite(cooldown) ? cooldown.toFixed(1) : 'the cooldown'} h of the trip; the planned ${plan.lead} h leaves ${Number.isFinite(cooldown) ? Math.max(cooldown - plan.lead, 0).toFixed(1) + ' h of no-touch time' : 'ample time'}. Without it the event reaches ${unprot.maxSub.toFixed(1)} °C of subcooling (risk index ${unprot.risk.toFixed(2)}).`);
  if (maxSub > 0) {
    recs.push(wReqIdeal >= 94 ? `No practical ${inh.name} dose removes ${maxSub.toFixed(1)} °C of subcooling: shorten the shutdown, depressurise, or displace the line before it cools.` : `Inject ${inh.name} to ${wReqIdeal.toFixed(0)} wt % of the aqueous phase (about ${inhRate.toFixed(1)} m³/d at ${(fEvent * 100).toFixed(0)} % rate) to stay ${v.marginC} °C outside the hydrate region; the present ${S.aq.inhWt.toFixed(0)} wt % leaves ${maxSub.toFixed(1)} °C of subcooling.`);
    if (cooldown <= 0) recs.push(`The line is already inside the hydrate region in steady production at ${xkm[iCool].toFixed(1)} km: continuous inhibition (or more insulation) is needed, not only shutdown procedures.`);
    else if (Number.isFinite(cooldown)) recs.push(`Treat ${cooldown.toFixed(1)} h as the cooldown limit: after a trip, start inhibitor displacement or depressurisation within about ${Math.max(cooldown - 2, 0).toFixed(1)} h (2 h reserved for the operation itself); the coldest point is at ${xkm[iCool].toFixed(1)} km.`);
    recs.push(`Screening for low-dosage inhibitors: a kinetic inhibitor is ${khiOk ? `plausible (subcooling ${maxSub.toFixed(1)} °C ≤ ${v.khiLimit} °C, indicative hold time ${khiHold > 1e4 ? 'very long' : khiHold.toFixed(0) + ' h'})` : `not suitable (subcooling ${maxSub.toFixed(1)} °C exceeds ${v.khiLimit} °C)`}; an anti-agglomerant is ${aaOk ? 'plausible' : 'not suitable'} at ${(prE.wcut * 100).toFixed(0)} % water cut — confirm either with qualification tests.`);
  } else if (!plan.active) recs.push(`The line stays ${(-maxSub).toFixed(1)} °C outside the hydrate region throughout this scenario; no hydrate inhibitor is needed for it. Check the shutdown and restart scenarios before relaxing the inhibition philosophy.`);
  if (plug || mc.prob > 0.1) recs.push(`If a plug forms, depressurise from both sides to ${v.depressP} bara: a ${v.plugLength} m plug then melts radially in about ${Number.isFinite(melt2.tNumeric) ? (melt2.tNumeric / 86400).toFixed(1) + ' d' : 'no finite time (the seabed is colder than the hydrate temperature at that pressure)'}. Never depressurise from one side only: ${Math.max(pLoc - v.depressP, 0).toFixed(0)} bar across the plug could launch it at the order of ${projV.toFixed(0)} m/s.`);
  if (piggingInterval !== null) recs.push(`Pig for wax every ${Math.max(Math.floor(piggingInterval * 0.8), 1)} d (80 % of the ${piggingInterval.toFixed(0)} d it takes to reach ${v.waxLimitMm} mm; ${slow.waxMass.toFixed(0)} kg of wax after ${v.depositDays} d).`);
  else recs.push(`No wax deposits in this operating mode: the wall stays above the wax appearance temperature of ${(+slow.wat).toFixed(1)} °C.`);
  if (scaleSI > 0 && gov) recs.push(gov.q.model
    ? `${gov.name} is supersaturated (SI ${gov.SI.toFixed(2)}) but ${scaleManaged ? 'kinetically held' : 'NOT held'}: nucleation takes ${hrs(gov.raw.t0)} untreated and ${hrs(gov.q.tInh)} with ${doseS} mg/L of ${v.scaleInh}, against ${v.scaleProtectH} h of residence. Minimum inhibitor concentration ${gov.q.mic.toPrecision(2)} mg/L; the index at which the untreated brine would scale within that time is ${gov.siCrit0.toFixed(2)} (${gov.siCrit.toFixed(2)} with the dose). Keep a safety factor of about 3 on the dose and confirm it with a dynamic tube-blocking test on the real brine.`
    : `${gov.name} is supersaturated (SI ${gov.SI.toFixed(2)}): select the inhibitor and its dose by a dynamic tube-blocking test; no induction-time model is available for this mineral.`);
  if (v.sandRate > 0) recs.push(sandBed ? `Raise the rate to at least ${sandMinRate.toFixed(0)} % of the case rate (mixture velocity ≥ ${sandCrit.toFixed(2)} m/s) to keep sand moving, or schedule sand pigging.` : `Sand keeps moving: the lowest velocity is ${sandMargin.toFixed(1)} times the critical velocity; do not run below about ${Math.min(sandMinRate, 100).toFixed(0)} % of the case rate for long periods.`);
  if (!oilSys) warnings.push({ level: 'info', msg: 'Asphaltene screening is not meaningful for a gas-dominated system with little stock-tank oil; the class is reported as low.' });
  if (aRiskRaw !== 'low' && aRisk === 'low') recs.push(`Asphaltene: the screening class of the oil is ${aRiskRaw} (colloidal instability index ${cii.cii.toFixed(2)}, de Boer "${db.cls}"${iFh >= 0 ? '' : ', no thermodynamic onset on the depletion path'}); it is managed here by a dispersant assumed ${v.aspInhEff} % effective. Qualify the chemical and its dose with a dispersant test on the real oil.`);
  if (aRisk !== 'low') recs.push(`Run an asphaltene onset test on live oil (depressurisation from ${pRes.toFixed(0)} bara through the bubble point near ${pBub.toFixed(0)} bara) before selecting an inhibitor.`);

  // ---- plots
  const tAx = ser.t, pCurve = logspace(Math.max(Math.min(...h.P, S.pOut) * 0.6, 2), Math.max(...base.P, ...ser.pIn) * 1.15, 40), plots = [];
  plots.push({ type: 'line', title: 'Pressure–temperature path over the hydrate curve', xlabel: 'Temperature (°C)', ylabel: 'Pressure (bara)', logy: true, series: [
    { name: 'Hydrate curve (case water and inhibitor)', x: pCurve.map((P) => S.hT0(P) - depNow), y: pCurve }, { name: 'Hydrate curve (fresh water)', x: pCurve.map((P) => S.hT0(P)), y: pCurve, dash: true },
    { name: 'Steady production', x: base.T, y: base.P }, { name: `End of the simulated event (${h.t.toFixed(1)} h)`, x: h.T, y: h.P, mode: 'both' }], note: 'Hydrates are stable to the left of the curve.' });
  plots.push({ type: 'line', title: 'Subcooling and hydrate fraction along the line', xlabel: 'Distance (km)', ylabel: 'Subcooling (°C) · hydrate fraction (vol %)', series: [
    { name: 'Peak subcooling during the event', x: xkm, y: rec.subMax }, { name: 'Subcooling at the end', x: xkm, y: rec.sub, dash: true }, { name: 'Steady production', x: xkm, y: subSteady, dash: true }, { name: 'Peak hydrate fraction (vol % of liquid)', x: xkm, y: rec.phiMaxT.map((f) => f * 100) }], hlines: [{ y: 0, label: 'hydrate curve' }] });
  if (unprot) plots.push({ type: 'line', title: 'Preserved line against the unprotected event: peak subcooling', xlabel: 'Distance (km)', ylabel: 'Peak subcooling (°C)', series: [{ name: `Preserved (${presName})`, x: xkm, y: rec.subMax }, { name: 'No preservation', x: xkm, y: unprot.subMax, dash: true }], hlines: [{ y: 0, label: 'hydrate curve' }, { y: -v.marginC, label: 'required margin' }], note: 'Positive values are inside the hydrate region.' });
  const dTg = linspace(1, 20, 39), tK = mean(S.tAmb) + KEL;
  { const aCell = D0 * ds, teqK = Math.max(...teqSteady) + KEL, o = { TeqK: teqK, mult: p.nucMult, B1: p.nucB }, ramp = onsetRamp({ rate: 1 / 3600, area: v.nucArea * 1e-4, ...o });
    plots.push({ type: 'line', title: 'Nucleation: induction time against subcooling', xlabel: 'Subcooling (°C)', ylabel: 'Mean induction time (h)', logy: true, ymin: 1e-3, ymax: 1e6, series: [
      { name: `Measured-rate model, one cell of the line (${aCell.toFixed(0)} m² of interface)`, x: dTg, y: dTg.map((d) => clamp(hydrateInductionTime(d, { ...o, area: aCell }) / 3600, 1e-3, 1e6)) },
      { name: `Measured-rate model, test cell (${v.nucArea} cm²)`, x: dTg, y: dTg.map((d) => clamp(hydrateInductionTime(d, { ...o, area: v.nucArea * 1e-4 }) / 3600, 1e-3, 1e6)), dash: true },
      { name: `Classical theory, contact angle ${v.contactAngle}° (${v.nucVolume} L of water)`, x: dTg, y: dTg.map((d) => clamp(inductionTime(d, { TK: tK, sigma: v.sigmaHW * 1e-3, theta: v.contactAngle, A: 3e7, het: true, V: p.nucV }) / 3600, 1e-3, 1e6)), dash: true }],
      vlines: maxSub > 0 ? [{ x: Math.min(maxSub, 20), label: 'peak subcooling' }] : [], note: `Mean waiting time for the first nucleus at constant subcooling; the onset itself is exponentially distributed about this mean. In a test cell cooled at 1 °C/h the measured-rate model gives onset at ${ramp.median.toFixed(1)} °C of subcooling (10–90 %: ${ramp.p10.toFixed(1)}–${ramp.p90.toFixed(1)} °C). Classical theory with bulk interfacial energy is shown for comparison only.` }); }
  if (pk.N) { const vt = pk.N.reduce((a, N, c) => a + N * grid.v[c], 0) || 1; plots.push({ type: 'line', title: `Particle-size distribution at the peak (${(S.x[pk.iPhi] / 1000).toFixed(1)} km, ${(pk.tPhi / 3600).toFixed(1)} h)`, xlabel: 'Particle / agglomerate size (µm)', ylabel: 'Volume fraction per class', logx: true, series: [{ name: 'Sectional population balance', x: grid.L.map((L) => L * 1e6), y: pk.N.map((N, c) => (N * grid.v[c]) / vt), mode: 'both' }], vlines: [{ x: p.dPrim * 1e6, label: 'primary' }, { x: Math.max(pk.dA[pk.iPhi], p.dPrim) * 1e6, label: 'cohesive limit' }] }); }
  plots.push({ type: 'line', title: 'Deposit thickness by type', xlabel: 'Distance (km)', ylabel: 'Thickness (mm)', zeroY: true, series: [{ name: 'Hydrate (worst time of the event)', x: xkm, y: dHyd.map((d) => d * 1000) }, { name: `Wax after ${v.depositDays} d`, x: xkm, y: slow.dWax.map((d) => d * 1000) }, { name: `Scale after ${v.depositDays} d`, x: xkm, y: slow.dScale.map((d) => d * 1000) }, { name: 'Total', x: xkm, y: total.map((d) => d * 1000), dash: true }] });
  plots.push({ type: 'line', title: 'Effective inner diameter', xlabel: 'Distance (km)', ylabel: 'Diameter (mm)', series: [{ name: 'With deposits', x: xkm, y: Deff.map((d) => d * 1000) }, { name: 'Clean bore', x: [xkm[0], xkm[n - 1]], y: [D0 * 1000, D0 * 1000], dash: true }] });
  plots.push({ type: 'line', title: 'Inlet pressure and blockage through the event', xlabel: 'Time (h)', ylabel: 'Inlet pressure (bara) · blockage (%) · viscosity factor', series: [{ name: 'Inlet pressure (bara)', x: tAx, y: ser.pIn }, { name: 'Peak area blockage (%)', x: tAx, y: ser.blk.map((b) => b * 100) }, { name: 'Peak slurry viscosity factor', x: tAx, y: ser.visc.map((x) => Math.min(x, 100)), dash: true }], hlines: [{ y: p.pInMax, label: 'available inlet pressure' }], vlines: tShut > 0 ? [{ x: tShut, label: 'restart' }] : [] });
  plots.push({ type: 'line', title: 'Hydrate inventory through the event', xlabel: 'Time (h)', ylabel: 'Hydrate mass (t)', zeroY: true, series: [{ name: 'Suspended in the stream', x: tAx, y: ser.susp.map((x) => x / 1000) }, { name: 'Deposited on the wall', x: tAx, y: ser.dep.map((x) => x / 1000) }] });
  if (mc.n > 0) { const data = mc.plugTimes.length >= 3 ? mc.plugTimes : mc.blk.map((b) => b * 100), hg = histogram(data, 10), isT = mc.plugTimes.length >= 3;
    plots.push({ type: 'bar', title: isT ? `Time to plug in ${mc.plugTimes.length} of ${mc.n} Monte Carlo samples` : `Peak hydrate blockage in ${mc.n} Monte Carlo samples`, ylabel: 'Samples', categories: hg.centers.map((c) => (isT ? c.toFixed(1) + ' h' : c.toFixed(1) + ' %')), series: [{ name: 'Samples', values: hg.counts }], note: `Plugging probability ${(mc.prob * 100).toFixed(0)} % (90 % interval ${(mc.lo * 100).toFixed(0)}–${(mc.hi * 100).toFixed(0)} %).` }); }
  plots.push({ type: 'line', title: 'Wax build-up at the worst location', xlabel: 'Time since the last pig run (d)', ylabel: 'Wax thickness (mm)', zeroY: true, series: [{ name: 'Maximum thickness', x: slow.ser.t, y: slow.ser.dMax }], hlines: [{ y: v.waxLimitMm, label: 'pigging limit' }] });
  const minIds = ['calcite', 'barite', 'celestite', 'gypsum', 'anhydrite', 'siderite'], nameOf = Object.fromEntries(MINERALS.map((mn) => [mn.id, mn.name])), siOf = (r, id) => Math.max(r.minerals.find((q) => q.id === id).SI, -6);
  plots.push({ type: 'line', title: 'Scale saturation index along the line', xlabel: 'Distance (km)', ylabel: 'Saturation index', series: minIds.map((id) => ({ name: nameOf[id], x: xkm, y: slow.sis.map((r) => siOf(r, id)) })), hlines: [{ y: 0, label: 'saturated' }] });
  plots.push({ type: 'line', title: 'Seawater mixing: saturation index at inlet conditions', xlabel: 'Seawater fraction (%)', ylabel: 'Saturation index', series: [...['calcite', 'barite', 'celestite', 'gypsum'].map((id) => ({ name: nameOf[id], x: mixF.map((f) => f * 100), y: mix.map((r) => siOf(r, id)) })), ...(v.srpOn ? [{ name: 'Barite with raw seawater (no sulphate removal)', x: mixF.map((f) => f * 100), y: mixRaw.map((r) => siOf(r, 'barite')), dash: true }] : [])], hlines: [{ y: 0, label: 'saturated' }], note: v.srpOn ? `Injection water after sulphate removal (${v.srpSO4} mg/L sulphate).` : 'Raw seawater.' });
  const rhoAx = linspace(550, 950, 9);
  plots.push({ type: 'line', title: 'de Boer asphaltene screening', xlabel: 'In-situ oil density (kg/m³)', ylabel: 'Reservoir pressure − saturation pressure (bar)', zeroY: true, series: [{ name: 'Slight problems above', x: rhoAx, y: rhoAx.map((r) => deBoer(r, 0).lower) }, { name: 'Severe problems above', x: rhoAx, y: rhoAx.map((r) => deBoer(r, 0).upper) }, { name: 'This oil', x: [clamp(aPath[0].rho, 550, 950)], y: [Math.max(pRes - pBub, 0)], mode: 'points' }], note: 'Boundaries interpolated between points read from a reproduction of the published screening plot (about ± 10 bar).' });
  plots.push({ type: 'line', title: 'Sand: critical velocity against actual velocity', xlabel: 'Distance (km)', ylabel: 'Velocity (m/s)', zeroY: true, series: [{ name: 'Mixture velocity', x: xkm, y: sand.map((s) => s.vm) }, { name: 'Oroskar–Turian', x: xkm, y: sand.map((s) => s.cr.oroskarTurian), dash: true }, { name: 'Salama', x: xkm, y: sand.map((s) => s.cr.salama), dash: true }, { name: 'Danielson', x: xkm, y: sand.map((s) => s.cr.danielson), dash: true }] });
  if (h.fld.t.length >= 2) {
    plots.push({ type: 'field', title: 'Hydrate fraction in distance and time', xlabel: 'Distance (km)', ylabel: 'Time (h)', zlabel: 'Hydrate fraction', zunit: 'vol %', x: xkm, y: h.fld.t, z: h.fld.phi.map((r) => r.map((f) => f * 100)), cmap: 'viridis', markers: plug ? [{ x: plug.x / 1000, y: plug.t / 3600, label: 'plug' }] : [] });
    plots.push({ type: 'field', title: 'Subcooling in distance and time', xlabel: 'Distance (km)', ylabel: 'Time (h)', zlabel: 'Subcooling', zunit: '°C', x: xkm, y: h.fld.t, z: h.fld.sub, cmap: 'coolwarm' });
    plots.push({ type: 'field', title: 'Hydrate deposit thickness in distance and time', xlabel: 'Distance (km)', ylabel: 'Time (h)', zlabel: 'Thickness', zunit: 'mm', x: xkm, y: h.fld.t, z: h.fld.dep.map((r) => r.map((d) => d * 1000)), cmap: 'thermal' });
  }
  plots.push({ type: 'bar', title: 'Governing solids threat by zone', ylabel: 'Severity index (1 = at the limit)', categories: zones.map((z) => `${z.from.toFixed(1)}–${z.to.toFixed(1)} km`), series: ['Hydrate', 'Wax', 'Scale', 'Sand'].map((k) => ({ name: k, values: zones.map((z) => +z.sc[k].toFixed(3)) })) });

  // ---- tables
  const pick = Array.from({ length: Math.min(12, n) }, (_, k) => Math.round((k * (n - 1)) / Math.max(Math.min(12, n) - 1, 1))), tables = [];
  tables.push({ title: 'Line stations', columns: ['Distance (km)', 'Steady T (°C)', 'End T (°C)', 'End P (bara)', 'Hydrate T (°C)', 'Peak subcooling (°C)', 'Exposure (h)', 'Peak hydrate (vol %)', 'Agglomerate d43 (µm)', 'Hydrate deposit (mm)', 'Wax (mm)', 'Scale (mm)', 'Deposit permeability (m²)', 'Controlling step at the end'],
    rows: pick.map((i) => [tx(xkm[i], 4), tx(base.T[i]), tx(h.T[i]), tx(h.P[i]), tx(rec.teq[i]), tx(rec.subMax[i]), tx(rec.exposure[i]), tx(rec.phiMaxT[i] * 100), tx((pk.d43 ? pk.d43[i] : 0) * 1e6), tx(dHyd[i] * 1000), tx(slow.dWax[i] * 1000), tx(slow.dScale[i] * 1000), dHyd[i] > 0 ? tx(perm[i]) : '—', rec.lim[i] || '—']) });
  tables.push({ title: 'Hydrate formation, transport and plugging', columns: ['Quantity', 'Value', 'Unit'], rows: [
    ['Scenario', SC_LABEL[v.scenario], ''], ['Hydrate system', p.mode === 'oil' ? 'oil-dominated (shrinking-core droplets)' : p.mode === 'gas' ? 'gas-dominated (film and entrained water)' : 'water-dominated (absorption limited)', ''],
    ['Stability margin (most negative = safest)', tx(-maxSub), '°C'], ['Exposure inside the hydrate region', tx(expoH), 'h'], ['Exposure integral', tx(degH), '°C·h'], ['Hydrate-stable length at the end', tx(stableEnd / 1000), 'km'],
    ['Onset location', h.onset ? tx(h.onset.x / 1000) : '—', 'km'], ['Onset time', tOnset === null ? '—' : tx(tOnset), 'h'], ['Peak nucleation rate', tx(rec.Jmax), p.nucModel === 'measured' ? '1/m² of gas–liquid interface/s' : '1/m³ water/s'], ['Peak formation rate', tx(peakRate), 'kg/s'],
    ['Hydrate formed', tx(led.formed / 1000), 't'], ['Hydrate dissociated', tx(led.dissociated / 1000), 't'], ['Hydrate carried out of the line', tx(led.exported / 1000), 't'], ['Captured on the wall', tx(led.captured / 1000), 't'], ['Sloughed back into the stream', tx(led.sloughed / 1000), 't'],
    ['Water converted to hydrate', tx(waterUsed / 1000), 't'], ['Gas consumed', tx(gasUsed / 1000), 't'], ['Peak water conversion in a cell', tx(convMax * 100), '%'], ['Heat of formation released (net of melting)', tx(led.heat / 1e9), 'GJ'], ['Hydration number', tx(p.hydN, 4), 'mol/mol'],
    ['Primary particle size', tx(p.dPrim * 1e6), 'µm'], ['Cohesive-limit agglomerate size at the peak', tx((pk.dA ? pk.dA[pk.iPhi] : 0) * 1e6), 'µm'], ['Agglomerate d43 at the peak', tx((pk.d43 ? pk.d43[pk.iPhi] : 0) * 1e6), 'µm'], ['Particle number concentration at the peak', tx(numPk), '1/m³ liquid'], ['Agglomerate settling velocity (hindered; negative = rises)', tx(setH.vHindered), 'm/s'], ['Velocity needed to keep the slurry suspended', tx(critH), 'm/s'], ['Liquid velocity at the peak', tx(pk.vL ?? 0), 'm/s'], ['Effective volume fraction at the peak', tx(pk.phiE ? pk.phiE[pk.iPhi] : 0), '–'], ['Peak slurry viscosity factor', tx(peakVisc), '×'],
    ['Peak wall-capture rate in a cell', tx(rec.captureMax), 'kg/s'], ['Peak shear-removal rate in a cell', tx(rec.removalMax), 'kg/s'], ['Deposit porosity (aged → fresh)', `${Math.min(...porPk).toFixed(2)} – ${Math.max(...porPk).toFixed(2)}`, '–'], ['Deposit permeability (Kozeny–Carman, range)', `${Math.min(...perm).toExponential(1)} – ${Math.max(...perm).toExponential(1)}`, 'm²'], ['Inlet-pressure rise in the event', tx(dpRise), 'bar'], ['Pressure-drop increase from wax and scale', tx(dpFoul - dpClean), 'bar'],
    ['Effective roughness with deposits', tx(roughEff * 1e6), 'µm'], ['Plug', plug ? `${(plug.t / 3600).toFixed(2)} h at ${(plug.x / 1000).toFixed(2)} km — ${plug.mech}` : 'none in the base case', ''],
    ['Plugging probability (90 % interval)', `${(mc.prob * 100).toFixed(0)} % (${(mc.lo * 100).toFixed(0)}–${(mc.hi * 100).toFixed(0)} %)`, ''], ['Plug location P10 / P50 / P90', mc.plugX.length ? `${(quantile(mc.plugX, 0.1) / 1000).toFixed(1)} / ${(quantile(mc.plugX, 0.5) / 1000).toFixed(1)} / ${(quantile(mc.plugX, 0.9) / 1000).toFixed(1)}` : '—', 'km'], ['Time to plug P10 / P50 / P90', plugT.length ? `${quantile(plugT, 0.1).toFixed(1)} / ${quantile(plugT, 0.5).toFixed(1)} / ${quantile(plugT, 0.9).toFixed(1)}` : '—', 'h'], ['Risk index', tx(risk), '0–1']],
    note: 'Times are measured from the start of the simulated sequence (the start of the shut-in for a cold restart).' });
  if (unprot) tables.push({ title: 'Preservation: design case against the unprotected event', columns: ['Quantity', 'Preserved line', 'No preservation', 'Unit'], rows: [
    ['Strategy', presName, 'none', ''], ['Inhibitor in the water of the line', tx(doseUsed), tx(Sun.aq.inhWt), 'wt %'], ['Inhibitor volume to dose the line', tx(presVol), '0', 'm³'], ['Heating power', tx((mx.holdPower || 0) / 1000), '0', 'kW'],
    ['Peak subcooling', tx(maxSub), tx(unprot.maxSub), '°C'], ['Exposure inside the hydrate region', tx(expoH), tx(unprot.expo), 'h'], ['Peak hydrate fraction of the liquid', tx(peakPhi * 100), tx(unprot.peakPhi * 100), 'vol %'], ['Peak slurry viscosity factor', tx(peakVisc), tx(unprot.peakVisc), '×'], ['Peak hydrate blockage', tx(peakBlkHyd * 100), tx(unprot.peakBlk * 100), '% of area'],
    ['Onset', tOnset === null ? 'none' : tx(tOnset), unprot.onset ? tx(unprot.onset.t / 3600) : 'none', 'h'], ['Plug in the base march', plug ? `${(plug.t / 3600).toFixed(1)} h` : 'none', unprot.plug ? `${(unprot.plug.t / 3600).toFixed(1)} h at ${(unprot.plug.x / 1000).toFixed(1)} km` : 'none', ''],
    ['Plugging probability', tx(outputs.plugProbability * 100), tx(unprot.prob * 100), '%'], ['Hydrate risk index', tx(risk), tx(unprot.risk), '0–1'], ['Time to put the preservation in place / cooldown time', `${plan.lead} / ${Number.isFinite(cooldown) ? cooldown.toFixed(1) : 'never'}`, '—', 'h']],
    note: `The unprotected event is marched on the same axial grid with twice the time step and ${unprot.n} Monte Carlo samples. An automatic dose is the one that leaves the required margin of ${v.marginC} °C at the coldest point and time of the preserved event.` });
  tables.push({ title: 'Population balance cross-check: sectional against quadrature method of moments', columns: ['Quantity', 'Sectional (fixed pivot)', 'QMOM (3 nodes)', 'Monte Carlo (direct simulation)'], rows: [['Batch time (s)', tx(qm.tB), tx(qm.tB), tx(qm.tB)], ['Number remaining N/N₀', tx(qm.m0s, 4), qm.ok ? tx(qm.m0q, 4) : '—', qm.mc ? tx(qm.mc.m0, 4) : '—'], ['Third moment m₃/m₃₀ (volume)', tx(qm.m3s, 6), qm.ok ? tx(qm.m3q, 6) : '—', qm.mc ? tx(qm.mc.vol, 6) : '—'], ['d43 (µm)', tx(qm.d43s * 1e6, 4), qm.ok ? tx(qm.d43q * 1e6, 4) : '—', qm.mc ? tx(qm.mc.d43 * 1e6, 4) : '—']], note: `Pure aggregation of a narrow population with the collision kernel of the worst cell over three collision times; all methods must conserve the third moment.${qm.mc ? ` The stochastic solution used ${v.mcParticles} simulated particles and ${qm.mc.events} collision trials (seed ${v.seed}).` : ''}` });
  if (cfd) { const F = cfd.F, r = cfd.r, um = (x) => tx(x * 1e6);
    tables.push({ title: 'Population balance on the resolved flow field (axisymmetric pipe section)', columns: ['Quantity', 'Value', 'Unit'], rows: [
      ['Population taken from', cfd.what, ''], ['Liquid velocity / slurry viscosity / density', `${cfd.U.toFixed(2)} / ${(cfd.mu * 1000).toFixed(2)} / ${cfd.rho.toFixed(0)}`, 'm/s / mPa·s / kg/m³'], ['Reynolds number', tx(F.Re), F.laminar ? 'laminar' : 'turbulent'], ['Friction velocity', tx(F.uStar), 'm/s'], ['Darcy friction factor of the resolved profile', tx(F.f, 4), `smooth-pipe law ${frictionFactor(F.Re, 0).toPrecision(4)}`],
      ['Centre-line velocity / bulk velocity', tx(F.fine.u[F.fine.u.length - 1] / F.U, 4), '–'], ['Shear rate: core cell / wall cell', `${F.gamma[0].toPrecision(3)} / ${F.gamma[F.nr - 1].toPrecision(3)}`, '1/s'], ['Turbulent dissipation: core cell / wall cell', `${F.epsT[0].toPrecision(3)} / ${F.epsT[F.nr - 1].toPrecision(3)}`, 'W/kg'], ['Largest eddy viscosity / molecular viscosity', tx(Math.max(...F.nut) / (cfd.mu / cfd.rho)), '–'],
      ['Resolved length', tx(cfd.L), `m (${v.cfdNx} steps × ${v.cfdNr} radial cells × ${grid.n} size classes)`], ['Agglomerate d43 at the outlet: core / wall / flow-weighted mean', `${um(cfd.dCore)} / ${um(cfd.dWall)} / ${um(cfd.dMean)}`, 'µm'], ['Agglomerate d43 of the line model in that cell', cfd.d1d > 0 ? um(cfd.d1d) : '—', 'µm'],
      ['Cohesive size limit: core / wall', `${um(r.kernels[0].dA)} / ${um(r.kernels[F.nr - 1].dA)}`, 'µm'], ['Particle volume fraction at the outlet: core / wall', `${cfd.phiCore.toPrecision(3)} / ${cfd.phiWall.toPrecision(3)}`, '–'], ['Particle number flux: inlet → outlet', `${r.ledger.numIn.toExponential(2)} → ${r.ledger.numOut.toExponential(2)}`, '1/s'],
      ['Wall deposition flux at the outlet', tx(r.dep[r.dep.length - 1] * 1e9), 'mm³/m²/s'], ['Deposition flux: resolved field / section-mean closure', cfd.depRatio === null ? '—' : tx(cfd.depRatio), '–'], ['Share of the particle volume deposited over the length', tx((100 * r.ledger.deposited) / Math.max(r.ledger.in, 1e-300)), '%']],
      note: `Reynolds-averaged momentum balance with the algebraic eddy viscosity of Cess (κ = ${RANS.kappa}, A⁺ = ${RANS.Aplus}); size classes carried by advection and turbulent diffusion (turbulent Schmidt number ${RANS.ScT}) with the collision and breakage kernels of the line model evaluated at the local mean shear and dissipation, and the deposition velocity of each class as the wall flux. The liquid is treated as filling the pipe at the liquid velocity of the cell. Three-dimensional and gas–liquid resolved cases are written by the external-solver hand-off.${r.limited ? ' Aggregation and breakage are stiff on this grid in part of the section; those cells used the implicit conservative scheme.' : ''}` });
    plots.push({ type: 'field', title: 'Flow-field population balance: agglomerate size d43', xlabel: 'Distance along the resolved section (m)', ylabel: 'Radius (mm)', zlabel: 'd43', zunit: 'µm', x: r.x, y: r.r.map((q) => q * 1000), z: r.r.map((_, j) => r.d43.map((row) => row[j] * 1e6)), cmap: 'viridis' });
    plots.push({ type: 'line', title: 'Flow-field population balance: radial profiles at the outlet', xlabel: 'Radius / pipe radius', ylabel: 'Normalised value', series: [{ name: 'Velocity / bulk velocity', x: F.rc.map((q) => q / F.R), y: F.u.map((q) => q / F.U) }, { name: 'Shear rate / wall-cell shear rate', x: F.rc.map((q) => q / F.R), y: F.gamma.map((q) => q / Math.max(F.gamma[F.nr - 1], 1e-30)) }, { name: 'd43 / flow-mean d43', x: F.rc.map((q) => q / F.R), y: r.d43[r.d43.length - 1].map((q) => q / Math.max(cfd.dMean, 1e-30)), mode: 'both' }, { name: 'Particle fraction / inlet fraction', x: F.rc.map((q) => q / F.R), y: r.phi[r.phi.length - 1].map((q) => q / Math.max(cfd.phi, 1e-30)), dash: true }] }); }
  tables.push({ title: 'Inhibition and remediation', columns: ['Item', 'Value', 'Unit', 'Basis'], rows: [
    ['Depression from salt', tx(depSalt), '°C', 'Nielsen–Bucklin on the water mole fraction'], [`Depression from ${S.aq.inhWt > 0 ? S.aq.inh.name : 'inhibitor'} now`, tx(depNow - depSalt), '°C', `${S.aq.inhWt.toFixed(1)} wt %, effectiveness ${v.inhEff} %`], ['Inhibitor effectiveness', tx(maxSub > 0 || depNow > depSalt ? (100 * (depNow - depSalt)) / Math.max(maxSub + depNow - depSalt, 1e-9) : 100), '% of the uninhibited subcooling removed', 'event peak'],
    [`${inh.name} required`, tx(wReqIdeal), 'wt %', `${v.marginC} °C margin on the peak subcooling`], ['Injection rate', tx(inhRate), 'm³/d', `${(mWater * 86.4).toFixed(1)} t/d of water at the event rate`], ['Loss to gas and hydrocarbon liquid', tx(lossPct), '% of injected', inhId === 'MeOH' ? 'partition coefficients (screening)' : 'negligible for glycols (screening)'],
    ['Kinetic inhibitor (screening)', khiOk ? 'plausible' : 'not suitable', '', `subcooling limit ${v.khiLimit} °C`], ['Anti-agglomerant (screening)', aaOk ? 'plausible' : 'not suitable', '', `water-cut limit ${v.aaWcLimit} %`],
    ['Two-sided depressurisation: radial melt time', Number.isFinite(melt2.tNumeric) ? tx(melt2.tNumeric / 3600) : 'no melting', 'h', `Stefan front, hydrate at ${tdTwo.toFixed(1)} °C, ambient ${tA.toFixed(1)} °C`], ['One-sided depressurisation: axial melt time', !Number.isFinite(melt1) ? 'no melting' : melt1 > 3.15e8 ? 'more than 10 years' : tx(melt1 / 86400), 'd', `Neumann solution over ${v.plugLength} m — not recommended`],
    ['Leak through a formed plug at the available pressure', tx(vLeak * A0 * 86400), 'm³/d', `Darcy–Forchheimer over ${v.plugLength} m, permeability ${kozenyCarman(v.porosityInf, grain).toExponential(1)} m², inertial share ${(100 * pgr.forchheimer / Math.max(pgr.total, 1e-300)).toFixed(0)} %`],
    ['Plug velocity if released one-sided', tx(projV), 'm/s', 'after 100 m of free travel, friction neglected'], ['Intrinsic dissociation time of a grain', Number.isFinite(tKin) ? tx(tKin) : '—', 's', 'Kim–Bishnoi: far faster than heat supply, so melting is heat-transfer-controlled'],
    ['Methanol to dissolve the plug', tx(meohMelt), 'm³', `${(wEq * 100).toFixed(0)} wt % in the released water`], ['Methanol contact time (mass-transfer-controlled)', tx(X.mt.film.tMelt / 86400), 'd', `film coefficient ${v.inhFilmK}×10⁻⁵ m/s on both plug faces (see the dissociation table)`], ['Heating to hold the line outside the region', tx(dehW / 1000), 'kW', `${tHold.toFixed(1)} °C along ${(S.L / 1000).toFixed(1)} km`], ['Heat to melt the plug in one day', tx(meltKW), 'kW', 'latent heat only']] });
  tables.push({ title: 'Wax', columns: ['Quantity', 'Value', 'Unit'], rows: [['Wax appearance temperature', tx(slow.wat), slow.waxCurve ? '°C (solid–liquid equilibrium of the n-paraffin distribution)' : '°C (entered)'], ['Precipitable wax', tx(slow.waxTotal * 100), 'wt % of the oil'], ['Solid wax at the seabed temperature', tx(waxSolubility(Math.min(...S.tAmb), slow.wat, slow.waxTotal, v.waxSlope, slow.waxCurve).solid * 100), 'wt % of the oil'], ['Length with a wall below the WAT', tx(slow.last.filter((r) => r && r.Ti < slow.wat).length * ds / 1000), 'km'], ['Initial build-up rate', tx(slow.waxRate0), 'mm/d'], [`Maximum thickness after ${v.depositDays} d`, tx(slow.ser.dMax[slow.ser.dMax.length - 1]), 'mm'], ['Wax mass in the line', tx(slow.waxMass), 'kg'], ['Wax fraction of the aged deposit', tx(Math.max(...slow.Fw)), '–'], ['Wax diffusivity at the wall', tx(Math.max(...slow.last.map((r) => r?.Dwo || 0))), 'm²/s'], ['Pigging interval', piggingInterval === null ? 'not needed' : tx(piggingInterval), 'd'], ['Gelled length after a cold shutdown', tx(slow.gelLen / 1000), 'km'], ['Gel yield stress', tx(slow.tauY), 'Pa'], ['Gel-breaking restart pressure', tx(slow.restartDp), 'bar']] });
  tables.push({ title: 'Scale at the inlet and at the worst location', columns: ['Mineral', 'SI at the inlet', 'SI at the worst location', 'Precipitation potential at the inlet (mg/L)', 'SI with 50 % seawater'], rows: minIds.map((id, k) => [nameOf[id], tx(siOf(sIn, id)), tx(Math.max(...slow.sis.map((r) => siOf(r, id)))), tx(sIn.minerals[k].ptb), tx(siOf(mix[5], id))]),
    note: `Ionic strength ${sIn.I.toFixed(2)} mol/kg, pH ${sIn.pH.toFixed(2)} at ${pMix.toFixed(0)} bara and ${tMix.toFixed(0)} °C with ${v.co2Pct} mol % CO₂. Maximum deposition rate ${Math.max(...slow.scaleRate).toFixed(3)} mm/y. Activity model: ${v.actModel === 'pitzer' ? 'Pitzer ion interaction (PHREEQC pitzer.dat)' : v.actModel === 'davies' ? 'ion association with Davies coefficients (PHREEQC phreeqc.dat pairs)' : 'ion association with Truesdell–Jones coefficients (PHREEQC phreeqc.dat)'}; water activity ${sIn.aw.toFixed(3)}.` });
  if (scaleCtl.length) tables.push({ title: 'Scale management: induction time and inhibitor', columns: ['Mineral', 'SI', 'At (km)', 'Induction time untreated', 'Induction time with the dose', 'Minimum inhibitor concentration (mg/L)', 'Critical SI untreated / with the dose', 'Residual risk (0–1)'],
    rows: scaleCtl.map((c) => [c.name, tx(c.SI), tx(xkm[c.i]), c.q.model ? hrs(c.raw.t0) : 'no model', c.q.model ? hrs(c.q.tInh) : 'no model', c.q.model ? (Number.isFinite(c.q.mic) ? tx(c.q.mic) : 'not reachable') : '—', c.q.model ? `${c.siCrit0.toFixed(2)} / ${c.siCrit.toFixed(2)}` : '—', tx(riskOf(c))]),
    note: `${v.scaleInhOn ? `${doseS} mg/L of ${SCALE_INHIBITORS[v.scaleInh].name}` : 'No inhibitor'}; the brine must stay scale-free for ${v.scaleProtectH} h. Barite and celestite: semi-empirical induction-time and threshold-inhibition relations (4–175 °C, SI 0.56–2.74, ± 0.5 log units); the wall deposit grows only after the induction time and at the uninhibited share of the growth sites. Other minerals are judged on the saturation index alone.${v.swFrac > 0 ? ` Seawater breakthrough ${v.swFrac} % with ${v.srpOn ? `${v.srpSO4} mg/L sulphate after sulphate removal` : 'raw seawater sulphate'}.` : ''}` });
  tables.push({ title: 'Asphaltene and sand', columns: ['Quantity', 'Value', 'Unit'], rows: [['In-situ oil density at reservoir conditions', tx(aPath[0].rho), 'kg/m³'], ['Undersaturation', tx(Math.max(pRes - pBub, 0)), 'bar'], ['de Boer class', db.cls, ''], ['Colloidal instability index', tx(cii.cii), cii.cls], ['Lowest Flory–Huggins solubility on the depressurisation path', tx(fhMin), 'vol fraction'], ['Asphaltene in the oil', tx(aspPhi), 'vol fraction'], ['Flory–Huggins upper onset pressure (density correlation)', aDen.upperOnset ? tx(aDen.upperOnset) : 'none', 'bara'], ['Flory–Huggins upper onset pressure (equation of state)', aEos ? (aEos.upperOnset ? tx(aEos.upperOnset) : 'none') : 'not evaluated (select the equation-of-state option or enter a measured onset)', aEos ? 'bara' : ''], ['Lower onset pressure (re-dissolution below the bubble point)', aUse.lowerOnset ? tx(aUse.lowerOnset) : 'none', 'bara'], ['Bubble point on the depletion path', aUse.pBubble ? tx(aUse.pBubble) : '—', 'bara'], ['Asphaltene solubility parameter used (25 °C)', tx(aUse.deltaA, 4), v.asphOnsetRef > 0 ? 'MPa^½, fitted to the measured onset' : 'MPa^½'], ['Asphaltene screening class of the oil', aRiskRaw, ''], ['Asphaltene deposition risk with the measures in place', aRisk, v.aspInhOn ? `dispersant, ${v.aspInhEff} % effective` : 'no inhibitor'],
    ['Sand concentration in the liquid', tx(sand[0].C * 1e6), 'ppm by volume'], ['Sand settling velocity (hindered)', tx(Math.abs(sand[0].st.vHindered)), 'm/s'], ['Particle response time', tx(relax.tau), 's'], ['Critical velocity (governing)', tx(sandCrit), 'm/s'], ['Lowest velocity ratio v / v_critical', tx(sandMargin), '–'], ['Length with a sand bed', tx(bedLen / 1000), 'km'], ['Sand hold-up in the line', tx(sandInv), 'kg'], ['Minimum rate to keep sand moving', tx(sandMinRate), '% of case rate'], ['Bend erosion (screening)', tx(erosion), 'mm/y']] });
  tables.push({ title: 'Governing threat by zone', columns: ['From (km)', 'To (km)', 'Hydrate', 'Wax', 'Scale', 'Sand', 'Governing'], rows: zones.map((z) => [tx(z.from), tx(z.to), tx(z.sc.Hydrate), tx(z.sc.Wax), tx(z.sc.Scale), tx(z.sc.Sand), z.top]), note: 'Severity indices: 1 means at the limit (10 °C subcooling or 15 % of the bore for hydrate, the pigging thickness for wax, SI 1.5 or 2 mm/y for scale, velocity at the critical velocity for sand).' });
  tables.push({ title: 'Hydrate thermodynamics check (van der Waals–Platteeuw, methane structure I)', columns: ['Quantity', 'Value', 'Unit'], rows: [['Mean seabed temperature', tx(mean(S.tAmb)), '°C'], ['Methane hydrate equilibrium pressure', vdw.P === null ? '—' : tx(vdw.P), 'bara'], ['Screening curve of the case gas at that temperature', tx(S.peq(mean(S.tAmb))), 'bara'], ['Small-cage occupancy at line pressure', tx(occ.thetaS), '–'], ['Large-cage occupancy at line pressure', tx(occ.thetaL), '–'], ['Hydration number from occupancy', tx(occ.hydrationNumber, 4), 'mol/mol']], note: 'The case gas contains propane and butanes and forms structure II at a lower pressure than pure methane; the statistical model is used for the cage occupancy and as an upper bound on the equilibrium pressure.' });

  // ---- results of the additional solvers
  for (const e of X.exp) { e.env.run = outputs; let q = NaN; try { q = experimentModel(e.type, e, v, e.env); } catch { q = NaN; } e.pred = q; }
  { const t = X.th, g = t.gm;
    tables.push({ title: 'Statistical hydrate model: Kihara cell potential and Gibbs-energy minimisation', columns: ['Quantity', 'Value', 'Unit'], rows: [
      ['Langmuir constant, small cage: Kihara potential / Parrish–Prausnitz fit', `${t.cK.s.toPrecision(3)} / ${t.cP.Cs.toPrecision(3)}`, '1/atm'], ['Langmuir constant, large cage: Kihara potential / Parrish–Prausnitz fit', `${t.cK.l.toPrecision(3)} / ${t.cP.Cl.toPrecision(3)}`, '1/atm'],
      ['Methane hydrate pressure with Kihara constants', t.kh.P === null ? '—' : tx(t.kh.P), `bara at ${t.tC.toFixed(1)} °C`], ['Methane hydrate pressure with Parrish–Prausnitz constants', t.pp.P === null ? '—' : tx(t.pp.P), 'bara'],
      ['Gibbs minimisation: water that can convert at the coldest point', tx(g.conversion * 100), `% at ${t.pC.toFixed(0)} bara`], ['Limiting factor', g.limiting === 'none' ? 'hydrate not stable (methane sI)' : g.limiting === 'equilibrium' ? 'salt and inhibitor concentrate until the water activity stops the conversion' : g.limiting + ' runs out', ''],
      ['Water mole fraction of the remaining aqueous phase', tx(g.xw, 4), '–'], ['Driving force Δg = (μ_w,hydrate − μ_w,liquid)/RT', tx(g.dg, 4), '–'], ['Gibbs energy change at the minimum', tx((g.Gmin - g.G0) / t.nW, 4), 'RT per mole of water']],
      note: 'Pure methane structure I: an upper bound on the pressure and a lower bound on the conversion for a gas that also forms structure II. The minimisation is closed for water, gas and dissolved particles at fixed pressure and temperature.' }); }
  { const s = X.slug, m = s.sm;
    tables.push({ title: 'Hydrate kinetics coupled to the slug unit cell', columns: ['Quantity', 'Value', 'Unit'], rows: m ? [
      ['Length in slug flow during the event', tx(s.len / 1000), 'km'], ['Evaluated at', tx(xkm[s.i]), 'km'], ['Slug frequency (kernel unit cell)', tx(m.freq), '1/s'], ['Slug frequency published by the flow suite', s.flowFreq === null ? '—' : tx(s.flowFreq), '1/s'], ['Slug-body length', tx(m.length), 'm'], ['Liquid holdup: slug body / film zone', `${m.holdupSlug.toFixed(2)} / ${m.holdupFilm.toFixed(2)}`, '–'],
      ['Share of the liquid inside slug bodies', tx(m.liquidInSlug * 100), '%'], ['Entrained bubble size in the slug body (Sauter)', tx(m.dBubble * 1e6), 'µm'], ['Interfacial area: slug body / film zone', `${m.aSlug.toPrecision(3)} / ${m.aFilm.toPrecision(3)}`, 'm²/m³ liquid'], ['Conductance k_L·a: slug body / film zone', `${m.kLaSlug.toPrecision(3)} / ${m.kLaFilm.toPrecision(3)}`, '1/s'],
      ['Unit-cell conductance against a stratified interface', tx(m.enhancement), '×'], ['Hydrate formed in slug-flow cells', tx(s.formedSlug / 1000), `t of ${tx(led.formed / 1000)} t`], ['Coupling', v.slugCouple ? 'on' : 'off (stratified interface everywhere)', '']] : [['Slug flow during the event', 'none', '']],
      note: 'Gas must dissolve before it can reach the hydrate particles; in slug flow the absorption step uses the liquid-volume-weighted conductance of slug body and film zone.' }); }
  { const inj = S.aq.injections || [], wC = S.aq.inhWtCell;
    tables.push({ title: 'Initial state and boundaries inside the time march', columns: ['Item', 'Value', 'Unit'], rows: [
      ['Chemical-injection points', inj.length ? inj.map((q) => `${q.q} m³/d at ${(q.x / 1000).toFixed(1)} km`).join('; ') : 'none', inj.length ? S.aq.inh.name : ''], ['Inhibitor in the water: inlet / outlet', `${wC[0].toFixed(1)} / ${wC[n - 1].toFixed(1)}`, 'wt %'], ['Hydrate temperature depression: least / most protected cell', `${Math.min(...depNowC).toFixed(1)} / ${Math.max(...depNowC).toFixed(1)}`, '°C'],
      ['Heated length', tx(mx.heatedLength / 1000), 'km'], ['Heating power', tx((v.heatWm * mx.heatedLength) / 1000), 'kW'], ['Heat supplied during the event', tx(led.heatIn / 3.6e9), 'MWh'], ['Depressurisation of the shut-in line', mx.depress ? `to ${mx.depress.P} bara after ${v.depressAtH} h` : 'none', ''],
      ['Equipment that retains solids', mx.traps.length ? mx.traps.map((t) => `${t.name} at ${(t.x / 1000).toFixed(1)} km (${(t.eff * 100).toFixed(0)} %)`).join('; ') : 'none', ''], ['Hydrate retained by equipment', tx(led.trapped), 'kg'], ['Initial suspended hydrate', tx(led.susp0), 'kg'], ['Free water at rest', mx.wMult ? `drained to the low points (largest local enrichment ${Math.max(...mx.wMult).toFixed(1)} ×)` : 'as left by the flowing holdup', ''], ['Water vapour in the gas', v.gasWaterSatPct, '% of saturation']] }); }
  if (X.track) { const t = X.track;
    tables.push({ title: `Lagrangian particle tracking (${t.what}, Maxey–Riley equation with turbulent dispersion)`, columns: ['Fate of the parcels', 'Parcels', 'Share (%)'], rows: [['Captured on a cold wall', t.wall, tx((100 * t.wall) / t.n)], ['Settled into a bed', t.bed, tx((100 * t.bed) / t.n)], ['Retained by equipment', t.equipment, tx((100 * t.equipment) / t.n)], ['Escaped through the outlet', t.escaped, tx((100 * t.escaped) / t.n)], ['Still in flight at the step limit', t.flying, tx((100 * t.flying) / t.n)], ...t.eq.map((e) => [`  retained at ${e.name} (${(e.x / 1000).toFixed(1)} km)`, e.n, tx((100 * e.n) / t.n)])],
      note: `${t.n} parcels, ${t.steps} steps, ${t.reflections} wall reflections; mean transit time of the escaped parcels ${t.meanTransit === null ? '—' : (t.meanTransit / 3600).toFixed(2) + ' h'}; median capture location ${t.depX.length ? (quantile(t.depX, 0.5) / 1000).toFixed(2) + ' km' : '—'}.` });
    plots.push({ type: 'bar', title: 'Where tracked parcels are captured', ylabel: 'Parcels', categories: xkm.map((x) => x.toFixed(1)), series: [{ name: 'Captured (wall, bed or equipment)', values: t.hist }], note: 'Distance in km; stochastic parcels released at the inlet.' }); }
  { const e = X.ee;
    if (e.ee) { plots.push({ type: 'line', title: `Sand bed after ${v.sandHours} h (Eulerian–Eulerian transport)`, xlabel: 'Distance (km)', ylabel: 'Bed height (mm) · suspended sand (ppm by volume)', zeroY: true, series: [{ name: 'Bed height (mm)', x: xkm, y: e.bedH.map((b) => b * 1000) }, { name: 'Suspended concentration (ppmv)', x: xkm, y: e.ee.c.map((c) => c * 1e6), dash: true }] });
      tables.push({ title: 'Eulerian–Eulerian sand transport', columns: ['Quantity', 'Value', 'Unit'], rows: [['Simulated period', v.sandHours, 'h'], ['Sand fed', tx(e.ee.ledger.in * v.sandRho), 'kg'], ['Sand carried out', tx(e.ee.ledger.out * v.sandRho), 'kg'], ['Sand in the bed', tx(e.bedMass), 'kg'], ['Sand in suspension', tx(e.ee.ledger.suspended * v.sandRho), 'kg'], ['Sand retained by equipment', tx(e.trapMass), 'kg'], ['Largest bed height', tx(e.bedH[e.iMax] * 1000), `mm at ${xkm[e.iMax].toFixed(1)} km`], ['Critical shear stress for resuspension', tx(Math.max(...e.tauC)), 'Pa'], ['Lowest wall shear stress', tx(Math.min(...pb.tauW)), 'Pa']], note: 'Solids phase velocity = liquid velocity minus the axial slip; deposition by hindered settling below the critical Shields stress, entrainment above it.' }); } }
  { const mt = X.mt; tables.push({ title: 'Mass-transfer-controlled dissociation by inhibitor contact', columns: ['Regime', 'Face recession rate (mm/h)', 'Time to dissolve the plug from both sides (d)', 'Basis'], rows: [
      ['Convective film at the face', tx(mt.film.rate * 3.6e6), tx(mt.film.tMelt / 86400), `film coefficient ${v.inhFilmK}×10⁻⁵ m/s`], ['Stagnant inhibitor column (diffusion only)', tx(mt.stag.rate * 3.6e6), mt.stag.tMelt / 86400 > 36500 ? 'more than 100 years' : tx(mt.stag.tMelt / 86400), `recession ${(mt.stag.s * 1000).toFixed(1)} mm in 30 d, similarity solution ${(mt.stag.sAnalytic * 1000).toFixed(1)} mm`], ['Film plus diffusion and reaction in the pores', tx(mt.por.rate * 3.6e6), tx(mt.por.tMelt / 86400), `Thiele modulus ${mt.por.thiele.toExponential(1)}`]],
      note: `Methanol at 98 wt % in the bulk; the interface stays at the equilibrium concentration ${(mt.wEq * 100).toFixed(0)} wt % at which hydrate is just stable at the plug; every kilogram of inhibitor that arrives must be diluted to that concentration by melt water.` }); }
  if (slow.cryst.length) tables.push({ title: 'Scale crystallisation in the bulk water (population balance, method of moments)', columns: ['Mineral', 'Crystals at the outlet (1/m³)', 'Number-mean size (µm)', 'Sauter size d32 (µm)', 'Precipitated in the bulk (mg/L)', 'Peak nucleation rate (1/m³/s)'], rows: slow.cryst.map((c) => [c.name, tx(c.number), tx(c.d10 * 1e6), tx(c.d32 * 1e6), tx(c.massMgL), tx(c.peakJ)]), note: `Classical nucleation with ${v.scaleSigma} mJ/m² interfacial energy and parabolic growth; crystals reach the wall with the turbulent deposition velocity and stick with efficiency ${v.scaleDepEff}.` });
  if (X.exp.length) { const ok = X.exp.filter((e) => Number.isFinite(e.pred)), byType = [...new Set(ok.map((e) => e.type))];
    tables.push({ title: 'Comparison with the entered measurements', columns: ['Type', 'a', 'b', 'c', 'Measured', 'Model', 'Model − measured'], rows: X.exp.map((e) => [e.type, tx(e.a), tx(e.b), tx(e.c), tx(e.meas, 4), tx(e.pred, 4), Number.isFinite(e.pred) ? tx(e.pred - e.meas) : 'unknown type']),
      note: byType.map((t) => { const q = ok.filter((e) => e.type === t), mt = metrics(q.map((e) => e.meas), q.map((e) => e.pred)); return `${t}: n = ${q.length}, bias ${mt.bias.toPrecision(3)}, RMSE ${mt.rmse.toPrecision(3)}${Number.isFinite(mt.mape) ? `, MAPE ${mt.mape.toFixed(1)} %` : ''}`; }).join(' · ') || 'No row has a known type.' });
    if (ok.length) plots.push({ type: 'line', title: 'Entered measurements: predicted against measured', xlabel: 'Measured', ylabel: 'Predicted', logx: ok.every((e) => e.meas > 0 && e.pred > 0), logy: ok.every((e) => e.meas > 0 && e.pred > 0), series: [...byType.map((t) => ({ name: t, x: ok.filter((e) => e.type === t).map((e) => e.meas), y: ok.filter((e) => e.type === t).map((e) => e.pred), mode: 'points' })), { name: 'Perfect agreement', x: [Math.min(...ok.map((e) => e.meas)), Math.max(...ok.map((e) => e.meas))], y: [Math.min(...ok.map((e) => e.meas)), Math.max(...ok.map((e) => e.meas))], dash: true }] }); }
  Object.assign(outputs, { watComputed: !!slow.waxCurve, inhibitorAlongLine: S.aq.inhWtCell.map((w) => +w.toFixed(3)), heatedLength: r3(mx.heatedLength, 5), heatInput: r3(led.heatIn / 3.6e9, 4), trappedAtEquipment: r3(led.trapped, 4), slugEnhancement: X.slug.sm ? r3(X.slug.sm.enhancement, 4) : null, slugLength: r3(X.slug.len, 5),
    methaneHydratePKihara: X.th.kh.P === null ? null : r3(X.th.kh.P, 4), gibbsConversion: r3(X.th.gm.conversion, 4), sandBedMass: r3(X.ee.bedMass, 4), sandBedHeight: r3(X.ee.bedH[X.ee.iMax], 4), sandBedX: X.ee.bedMass > 0 ? r3(S.x[X.ee.iMax], 5) : null,
    asphalteneOnsetPressure: aUse.upperOnset ? r3(aUse.upperOnset, 4) : null, asphalteneOnsetEos: aEos && aEos.upperOnset ? r3(aEos.upperOnset, 4) : null, inhibitorMeltTime: r3(X.mt.film.tMelt / 3600, 4), scaleCrystalSize: slow.cryst.length ? r3(Math.max(...slow.cryst.map((c) => c.d32)), 3) : null,
    pbeFieldD43: cfd ? r3(cfd.dMean, 4) : null, pbeFieldD43Core: cfd ? r3(cfd.dCore, 4) : null, pbeFieldD43Wall: cfd ? r3(cfd.dWall, 4) : null, pbeFieldDeposition: cfd ? r3(cfd.r.dep[cfd.r.dep.length - 1], 4) : null, pbeFieldDepositionRatio: cfd && cfd.depRatio !== null ? r3(cfd.depRatio, 4) : null, pbeFieldFriction: cfd ? r3(cfd.F.f, 4) : null,
    trackedCaptured: X.track ? r3((X.track.wall + X.track.bed + X.track.equipment) / X.track.n, 4) : null, pbeMonteCarlo: qm.mc ? r3(qm.mc.m0, 4) : null });
  const flowArea = Deff.reduce((a, d) => a + (PI / 4) * d * d * ds, 0), depArea = total.reduce((a, d, i) => a + (PI / 4) * (D0 * D0 - Deff[i] ** 2) * ds, 0);
  const balances = [
    { name: 'Hydrate mass (kg): formed + initial + inflow = suspended + deposited + exported + dissociated', in: led.in, out: led.out },
    { name: 'Water + gas consumed (kg) = hydrate formed (hydration number ' + p.hydN.toFixed(2) + ')', in: led.formed, out: waterUsed + gasUsed },
    { name: 'Latent heat (J): heat put into the fluid = (formed − dissociated) × latent heat', in: (led.formed - led.filmWall - led.dissociated) * HYDRATE.latent, out: led.heat },
    { name: 'Pipe volume (m³): flow area + deposit area = clean bore', in: A0 * S.L, out: flowArea + depArea },
  ];
  if (X.ee.ee) { const l = X.ee.ee.ledger; balances.push({ name: 'Sand volume, Eulerian–Eulerian transport (m³): fed + initial = out + suspended + bed + retained', in: l.in + l.initial, out: l.out + l.suspended + l.bed + l.trapped }); }
  if (X.track) balances.push({ name: 'Tracked parcels: released = wall + bed + equipment + escaped + in flight', in: X.track.n, out: X.track.wall + X.track.bed + X.track.equipment + X.track.escaped + X.track.flying });
  if (cfd) balances.push({ name: 'Flow-field population balance (m³/s of particles): inlet flux = outlet flux + wall deposition', in: cfd.r.ledger.in, out: cfd.r.ledger.out + cfd.r.ledger.deposited });
  if (qm.mc) balances.push({ name: 'Monte Carlo population balance: particle volume after aggregation / before', in: 1, out: qm.mc.vol });
  for (const c of slow.cryst) balances.push({ name: `${c.name} crystallisation (mol/m³): ions removed from solution = crystal volume × density / molar mass`, in: c.x, out: c.molCheck });
  { const g = X.th.gm; balances.push({ name: 'Gibbs minimisation: converted water at the numerical minimum = stationary point ln x_w = Δg (mol per mol water fed)', in: g.xiAnalytic / X.th.nW, out: g.xi / X.th.nW }); }
  const summary = `${SC_LABEL[v.scenario]}: ${maxSub > 0 ? `up to ${maxSub.toFixed(1)} °C of subcooling over ${(stableLen / 1000).toFixed(1)} km, hydrate reaches ${(peakPhi * 100).toFixed(1)} vol % of the liquid and ${plug ? `the line plugs after ${(plug.t / 3600).toFixed(1)} h at ${(plug.x / 1000).toFixed(1)} km` : `no plug forms in the base case (plugging probability ${(mc.prob * 100).toFixed(0)} %)`}` : `the line stays ${(-maxSub).toFixed(1)} °C outside the hydrate region${plan.active ? ` with the preservation in place (${presName}${doseUsed > 0 ? `, ${doseUsed.toFixed(0)} wt %` : ''})` : ''}`}${unprot && unprot.maxSub > 0 ? ` — without preservation the same event reaches ${unprot.maxSub.toFixed(1)} °C of subcooling and a risk index of ${unprot.risk.toFixed(2)}` : ''}; wax ${slow.waxRate0 > 1e-4 ? `builds at ${slow.waxRate0.toFixed(3)} mm/d` : 'does not deposit'}, ${scaleSI > 0 ? `${scaleMineral} is supersaturated (SI ${scaleSI.toFixed(2)}) ${scaleManaged ? 'but kinetically held' : 'and not held'}${doseS > 0 ? ` with ${doseS} mg/L of inhibitor` : ''}` : 'no mineral is supersaturated'}, asphaltene risk is ${aRisk}${aRisk !== aRiskRaw ? ` with the dispersant (screening class ${aRiskRaw})` : ''} and sand ${sandBed ? 'settles' : 'keeps moving'}.`;
  prog(1, 'Done');
  return { summary, kpis, warnings, recommendations: recs, plots, tables, balances, outputs };
}

export const runSolids = run;
// =====================================================================================================
// 16. Calibration model: a jacketed flow loop at constant subcooling (hydrate) and a cold-wall wax test
// =====================================================================================================
/**
 * Fast loop model for parameter estimation. v: input values plus the operating point of a test:
 * calT (h of test), calDT (°C subcooling without inhibitor), calInh (wt % MEG in the water), calV (m/s loop velocity), calDTw (°C oil-to-coolant
 * difference of the wax test).
 * Returns { phi (hydrate vol % of the liquid), dpRatio (pressure drop relative to the hydrate-free loop), waxMm (wax thickness) }.
 */
export function loopModel(v0) {
  const v = { ...DEF, ...v0 }, tH = num(v0.calT, 4), dT0 = num(v0.calDT, 8), dT = Math.max(dT0 - (v.inhEff / 100) * hydrateDepression({ S: 0, inhWt: clamp(num(v0.calInh, 0), 0, 90), inh: INHIBITORS.MEG }), 0), vel = num(v0.calV, 1.5), dTw = num(v0.calDTw, 15), D = 0.0508, wc = 0.2, P = 80;
  const rho = 820, mu = 4e-3, nu = mu / rho, sg = 0.7, Teq = hydrateT0(P, sg) + num(v.hydEqOffset, 0), TK = Teq - dT + KEL, Peq = brent((p) => hydrateT0(p, sg) + num(v.hydEqOffset, 0) - (Teq - dT), 1, 700, 1e-6), dp = v.primaryUm > 0 ? v.primaryUm * 1e-6 : 40e-6;
  const hydN = v.hydNumber > 0 ? v.hydNumber : 6, mHyd = 0.018 + hydN * MW_W, rhoH = v.rhoHyd > 0 ? v.rhoHyd : HYDRATE.rho, vmh = mHyd / rhoH, eexp = (1000 / rhoH) / ((hydN * MW_W) / mHyd);
  const Re = (rho * vel * D) / mu, f = frictionFactor(Re, 4.5e-5 / D), tauW = (f / 8) * rho * vel * vel, eps = (4 * tauW * vel) / (rho * D), gdot = Math.sqrt(eps / nu), Dg = (7.4e-12 * Math.sqrt(150) * TK) / (mu * 1000 * 37.7 ** 0.6), Sh = 2 + 0.6 * Math.sqrt((gdot * dp * dp) / nu) * Math.cbrt(nu / Dg);
  const kLa = 0.4 * (eps * nu) ** 0.25 * Math.sqrt(Dg / nu) * 40 * v.mtMult, nd = wc / ((PI / 6) * dp ** 3), n = 80, dt = (tH * 3600) / n, hJ = hInside(Re, (2000 * mu) / 0.13, 0.13, D) * (v.htMult ?? 1);
  let X = 1e-4, dep = 0, phi = 0, mur = 1, dEq = dp, rMax = 0, por = v.porosity0;
  for (let k = 0; k < n; k++) {
    const gr = hydrateGrowthRate({ TK, P, Peq, zG: 0.85, kRef: v.kinK * 1e-10, EaR: v.kinEa, H: 2500, kFilm: (Sh * Dg) / dp, kShell: shellConductance(X, dp / 2, v.shellD * 1e-15), hPart: (2 * 0.14) / dp, dT, dHmol: HYDRATE.latent * mHyd });
    const Ra = gr.dc > 0 ? 1 / (v.kinK * 1e-10 * Math.exp(-v.kinEa * (1 / TK - 1 / 277.15)) * 2500) + dp / (Sh * Dg) + gr.dc / gr.jHeat : 0, X1 = gr.dc > 0 ? shrinkingCoreAdvance(X, (gr.dc * dt * 6 * vmh) / (dp * eexp), Ra, dp / 2, v.shellD * 1e-15) : X;
    const rP = ((X1 - X) * wc * eexp) / (vmh * dt), r = gr.dc > 0 ? rP / (1 + rP / (kLa * gr.dc)) : 0; if (r > rMax) rMax = r; // step-mean particle rate (exact shrinking-core advance) in series with gas absorption
    X = Math.min(X + (r * vmh * dt) / (wc * eexp), 1); phi = X * wc * eexp;
    const cp = maxAgglomerateSize({ dp, Fa: v.cohesion * 1e-3 * dp, mu0: mu, shear: gdot, phi, phiMax: v.phiMax, fr: v.fractal });
    // agglomerate size where the collision and breakage frequencies of the kernels balance: (2.5 α φ / k_b)^(1/3) of the cohesive limit
    dEq = clamp(cp.dA * Math.cbrt((2.5 * v.aggEff * phi) / Math.max(v.kBreak, 1e-9)), dp, 1.5 * cp.dA); const phiE = Math.min(phi * (dEq / dp) ** (3 - v.fractal), 0.98 * v.phiMax);
    mur = slurryViscosity(phiE, v.viscModel, { phiMax: v.phiMax });
    const uS = Math.sqrt(tauW / rho), Vd = depositionVelocity(dEq, rhoH, uS, nu, rho, TK) * v.adhesion * Math.min(1, (v.adhForce * 1e-3) / (8 * tauW * dEq)), tauC = v.tauCrit * ((1 - por) / (1 - v.porosity0)) ** 2;
    dep = Math.min(dep + (Vd * phi * dt) / (1 - v.porosity0), 0.45 * D); if (tauW > tauC && v.kRemove > 0) dep *= Math.exp((-v.kRemove / 3600) * (tauW / tauC - 1) * dt); // shear removal above the critical wall shear
    por = v.porosityInf + (por - v.porosityInf) * Math.exp(-dt / (v.ageHours * 3600));
  }
  // wax: oil 3 °C above the WAT in the bulk, coolant dTw below the oil, film coefficient of the loop
  const hIn = hInside(Re, (2000 * mu) / 0.13, 0.13, D), Tb = v.wat + 3; let dW = 0, Fw = clamp(1 - (v.waxOil > 0 ? v.waxOil : 80) / 100, 0.03, 0.95);
  for (let k = 0; k < 40; k++) { const r = waxDeposition({ wmodel: v.waxModel, c1: v.waxC1, coil: v.waxOil > 0 ? v.waxOil / 100 : null, Tb, Tamb: Tb - dTw, U: 1 / (1 / hIn + 1 / 400), hIn, kOil: 0.13, rhoOil: rho, muOil: mu, wat: v.wat, wTot: v.waxContent / 100, slope: v.waxSlope, delta: dW, Fw, kDep: v.waxK, D, vL: vel, gammaW: tauW / mu, mult: v.waxMult, diffModel: v.waxDiff, stripC: v.waxStripC, aspect: v.waxAspect }); dW += (r.dDelta * tH * 3600) / 40; Fw = Math.min(Fw + (r.dFw * tH * 3600) / 40, 0.95); }
  // other observables of the same loop and of bench tests at the conditions of the row
  const grain = v.plugGrainUm > 0 ? v.plugGrainUm * 1e-6 : dp, sup = num(v0.calSuper, 2), TKd = Teq + sup + KEL, PeqD = brent((p) => hydrateT0(p, sg) + num(v.hydEqOffset, 0) - (Teq + sup), 1, 900, 1e-6), jd = hydrateDissociationRate({ TK: TKd, P, Peq: PeqD, zG: 0.85, K0: v.disK0 * 1e4, E: v.disE * 1000 });
  const SI = num(v0.calSI, 1), Ssat = 10 ** (SI / 2), TKs = 343.15, cr = scaleCrystallisation({ S0: Ssat, cA: 1, cB: 1, M: 0.23339, rho: 4480, TK: TKs, t: tH * 3600, sigma: v.scaleSigma * 1e-3, A: 10 ** v.scaleNucA, kg: v.scaleKg * 1e-10 * Math.exp((-30800 / R) * (1 / TKs - 1 / 298.15)), steps: 40, deplete: false });
  const kinS = v.scaleK * 1e-8 * Math.exp((-30800 / R) * (1 / TKs - 1 / 298.15)) * (Ssat - 1) ** 2, dS = v.sandUm * 1e-6, rhoWl = 1000, muWl = 1e-3, cv = sandCriticalVelocity({ d: dS, D, rhoP: v.sandRho, rhoF: rhoWl, mu: muWl, C: 1e-3, vsl: 1, vm: 1, coef: v.sandCoef, shape: v.sandShape });
  const tauCs = (v.sandShields > 0 ? v.sandShields : shieldsCritical(dS, v.sandRho, rhoWl, muWl)) * (v.sandRho - rhoWl) * G * dS, fW = frictionFactor((rhoWl * vel * D) / muWl, 4.5e-5 / D), tauWs = (fW / 8) * rhoWl * vel * vel, vRes = Math.sqrt((8 * tauCs) / (frictionFactor((rhoWl * Math.max(cv.governing, 0.1) * D) / muWl, 4.5e-5 / D) * rhoWl));
  const ap = loopAsphPath(), ao = asphalteneOnset({ pre: ap, phiA: ((v.saraAsp / 100) * 800) / 1200, deltaA: v.asphDelta, vA: v.asphMV });
  return { phi: phi * 100, conv: X * 100, dpRatio: mur ** 0.25 * (D / (D - 2 * dep)) ** 4.75, waxMm: dW * 1000,
    teqC: Teq - (v.inhEff / 100) * hydrateDepression({ S: 0, inhWt: clamp(num(v0.calInh, 0), 0, 90), inh: INHIBITORS.MEG }), tIndH: Math.min((v.nucleation === 'heterogeneous' ? hydrateInductionTime(Math.max(dT, 1e-6), { TeqK: Teq + KEL, area: v.nucArea * 1e-4, mult: 10 ** v.nucA, B1: v.nucB }) : inductionTime(Math.max(dT, 1e-6), { TK, sigma: v.sigmaHW * 1e-3, theta: v.contactAngle, A: 10 ** (v.nucA + (v.nucleation === 'homogeneous' ? 35 : 7.5)), het: v.nucleation !== 'homogeneous', V: v.nucVolume * 1e-3 })) / 3600, 1e9),
    kLa, dTex: (rMax * HYDRATE.latent * mHyd * D) / (4 * hJ), d43Um: dEq * 1e6, muRel: mur, depMm: dep * 1000, permM2: kozenyCarman(por, grain), tDisS: jd > 0 ? (rhoH * dp) / (6 * mHyd * jd) : 1e9,
    waxSolidPct: waxSolubility(v.wat - dTw, v.wat, v.waxContent / 100, v.waxSlope).solid * 100, scaleUm: (cr.G0 * tH * 3600) * 1e6, scaleMmY: (kinS * 0.23339 * 3.156e10) / (4480 * 0.8) + v.scaleDepEff * cr.massConc * 3.156e10 * depositionVelocity(Math.max(cr.d32, 1e-8), 4480, Math.sqrt(tauWs / rhoWl), 1e-6, rhoWl, TKs) / (4480 * 0.8),
    aopBar: ao.upperOnset ?? 0, vSetMm: settlingVelocity(dS, v.sandRho, rhoWl, muWl, { shape: v.sandShape }).v * 1000, vCrit: cv.governing, vResus: vRes, entrRate: v.sandEntrain * 1e-4 * Math.max(tauWs / tauCs - 1, 0) * v.sandRho, erosMmY: sandErosionScreen(v.sandRate, vel, v.sandUm, D * 1000, rhoWl, v.erosionSm) };
}
/**
 * Water conversion (%) in a stirred autoclave holding a water-in-oil emulsion under pure methane at constant pressure and
 * temperature: every droplet converts as a shrinking core (surface kinetics, liquid film and the growing hydrate shell in
 * series) under the fugacity difference between the cell pressure and the methane hydrate pressure of the statistical model.
 * o: { tH (h since nucleation), Tc (°C), P (bara), dUm (droplet size, µm; default 40) }.
 */
export function autoclaveConversion({ tH, Tc, P, dUm = 40 }, v0 = {}) {
  const v = { ...DEF, ...v0 }, TK = Tc + KEL, Peq = vdwpMethane(Tc).P, dp = dUm * 1e-6, hydN = v.hydNumber > 0 ? v.hydNumber : 6, mHyd = 0.016 + hydN * MW_W, rhoH = v.rhoHyd > 0 ? v.rhoHyd : HYDRATE.rho, vmh = mHyd / rhoH, eexp = (1000 / rhoH) / ((hydN * MW_W) / mHyd);
  if (!(Peq > 0) || !(P > Peq)) return 0;
  const dc = (methaneFugacity(P, TK) - methaneFugacity(Peq, TK)) * 1e5 / 2500, kK = v.kinK * 1e-10 * Math.exp(-v.kinEa * (1 / TK - 1 / 277.15)) * 2500, Dg = (7.4e-12 * Math.sqrt(150) * TK) / (4 * 37.7 ** 0.6);
  return 100 * shrinkingCoreAdvance(0, (dc * tH * 3600 * 6 * vmh) / (dp * eexp), 1 / kK + dp / (2 * Dg), dp / 2, v.shellD * 1e-15);
}
let loopAsph = null;
/** Live-oil solubility-parameter path of the reference fluid used by the calibration model (computed once). */
function loopAsphPath() { return (loopAsph ||= asphaltenePath({ spec: DEFAULT_FLUID, pRes: num(DEFAULT_FLUID.Pres, 300), tRes: num(DEFAULT_FLUID.Tres, 90), pEnd: 5, n: 24, method: 'eos' })); }
/**
 * Model prediction for one laboratory or field measurement (the comparison path of the Measurements table and of the
 * reference data sets). type: see the Measurements input; r: { a, b, c } conditions; v: input values; env: { run (outputs of this run), spec (case fluid), yCO2 }.
 * Returns a number (NaN when the type is unknown).
 */
export function experimentModel(type, r, v0 = {}, env = {}) {
  const v = { ...DEF, ...v0 }, a = num(r?.a, 0), b = num(r?.b, 0), c = num(r?.c, 0), t = String(type || '').trim().toLowerCase(), lm = (o) => loopModel({ ...v, ...o });
  switch (t) {
    case 'equilibrium': return hydrateT0(Math.max(a, 1.5), b > 0 ? b : 0.7) + num(v.hydEqOffset, 0);
    case 'equilibrium-methane': { const q = vdwpMethane(a, { langmuir: v.langmuirModel === 'kihara' || b === 2 ? 'kihara' : 'parrish' }); return q.P ?? NaN; } // a = temperature °C → pressure bara
    case 'inhibitor': return (v.inhEff / 100) * hydrateDepression({ S: 0, inhWt: clamp(a, 0, 90), inh: b === 1 ? INHIBITORS.MeOH : b === 3 ? INHIBITORS.EtOH : INHIBITORS.MEG });
    case 'onset': { // linear cooling at a °C/h from the hydrate temperature: median subcooling at onset (the hazard ∫ J α dt reaches ln 2); b = gas–water interface (cm²), c = hydrate temperature (°C)
      if (v.nucleation === 'heterogeneous') return onsetRamp({ rate: Math.max(a, 1e-6) / 3600, area: (b > 0 ? b : v.nucArea) * 1e-4, TeqK: (c > 0 ? c : 10) + KEL, mult: 10 ** v.nucA, B1: v.nucB }).median;
      const rate = Math.max(a, 1e-6) / 3600, V = v.nucVolume * 1e-3; let hz = 0, dT = 0; const st = 0.02;
      while (dT < 40 && hz < Math.LN2) { dT += st; hz += (nucleationRate({ TK: 277.15, dT, TeqK: 277.15 + dT, sigma: v.sigmaHW * 1e-3, theta: v.contactAngle, A: 10 ** (v.nucA + (v.nucleation === 'homogeneous' ? 35 : 7.5)), het: v.nucleation !== 'homogeneous' }).J * V * st) / rate; }
      return dT; }
    case 'induction': return Math.min((v.nucleation === 'heterogeneous' ? hydrateInductionTime(Math.max(a, 1e-6), { TeqK: (c > 0 ? c : 10) + KEL, area: (b > 0 ? b : v.nucArea) * 1e-4, mult: 10 ** v.nucA, B1: v.nucB }) : inductionTime(Math.max(a, 1e-6), { TK: 277.15, sigma: v.sigmaHW * 1e-3, theta: v.contactAngle, A: 10 ** (v.nucA + (v.nucleation === 'homogeneous' ? 35 : 7.5)), het: v.nucleation !== 'homogeneous', V: v.nucVolume * 1e-3 })) / 3600, 1e9);
    case 'autoclave': case 'rockingcell': case 'flowloop': case 'growth': return lm({ calT: a, calDT: b, calV: c > 0 ? c : t === 'rockingcell' ? 0.3 : 1.5 }).conv; // % of the water converted
    case 'loopdp': return lm({ calT: a, calDT: b, calV: c > 0 ? c : 1.5 }).dpRatio;
    case 'plugging': { const f = (th) => lm({ calT: th, calDT: a, calV: b > 0 ? b : 1.5 }).dpRatio - 5; if (f(500) < 0) return 500; if (f(0.01) > 0) return 0.01; return brent(f, 0.01, 500, 1e-3, 60); } // time for the pressure drop to reach five times the clean value
    case 'psd': case 'agglomeration': return lm({ calT: a, calDT: b, calV: c > 0 ? c : 1.5 }).d43Um;
    case 'rheology': { const phi = a / 100, dp = v.primaryUm > 0 ? v.primaryUm * 1e-6 : 40e-6, cp = b > 0 ? maxAgglomerateSize({ dp, Fa: v.cohesion * 1e-3 * dp, mu0: c > 0 ? c * 1e-3 : 4e-3, shear: b, phi, phiMax: v.phiMax, fr: v.fractal }) : { phiEff: phi }; return slurryViscosity(cp.phiEff, v.viscModel, { phiMax: v.phiMax }); }
    case 'suspension': return slurryViscosity(a / 100, b === 2 ? 'thomas' : b === 3 ? 'mills' : 'krieger', { phiMax: c > 0 ? c : 0.64 }); // hard spheres: a = vol %, b = model code, c = packing fraction
    case 'deposition': return lm({ calT: a, calDT: b, calV: c > 0 ? c : 1.5 }).depMm;
    case 'adhesion': return v.adhForce * 1e-3 * a * 1e-6 * 1e6; // µN for a particle of a µm
    case 'dissociation': { const q = lm({ calSuper: Math.max(b, 0.01), primaryUm: c > 0 ? c : v.primaryUm }); return 100 * (1 - Math.max(1 - (a * 60) / q.tDisS, 0) ** 3); } // shrinking sphere at the intrinsic rate
    case 'restart': { const S = labSetup({ n: 8, L: 400, D: 0.0508, T: 4, P: 80, U: 20, tAmb: 4, vsl: (c > 0 ? c : 1) * 0.6, vsg: (c > 0 ? c : 1) * 0.4, sg: 0.7 }), teq = hydrateT0(80, 0.7), S2 = labSetup({ n: 8, L: 400, D: 0.0508, T: teq - Math.max(b, 0.1), P: 80, U: 20, tAmb: teq - Math.max(b, 0.1), vsl: (c > 0 ? c : 1) * 0.6, vsg: (c > 0 ? c : 1) * 0.4, sg: 0.7 }); void S;
      const p = hydrateParams({ ...v, regime: 'oil', nucA: Math.max(v.nucA, 12) }, S2), q = runHydrateMarch(S2, p, { phases: [{ dur: Math.max(a, 0.05) * 3600, dt: Math.max(a * 3600 / 12, 60), frac: 0 }, { dur: 2 * 3600, dt: 120, frac: 1 }], grid: pbeGrid(10, 5e-6, 5e-3), cheap: true, rows: 2 }), fl = q.ser.dp.filter((_, k) => q.ser.frac[k] > 0);
      void fl; let pk = 1; q.ser.visc.forEach((m, k) => { if (q.ser.frac[k] > 0) pk = Math.max(pk, m ** 0.25 / (1 - Math.min(q.ser.blk[k], 0.97)) ** 2.375); }); return pk; } // friction scaling with the slurry viscosity and the restricted bore
    case 'field-plugtime': return env.run?.plugTime ?? env.run?.plugTimeP50 ?? NaN;
    case 'field-plugx': return (env.run?.plugX ?? env.run?.plugXP50 ?? NaN) / 1000;
    case 'field-dp': return env.run?.dpIncrease ?? NaN;
    case 'waxloop': return lm({ calT: a, calDTw: b, calV: c > 0 ? c : 1.5 }).waxMm;
    case 'scale-barite': case 'scale-celestite': case 'scale-gypsum': case 'scale-anhydrite': return 1e6 * mineralSolubility(t.slice(6), { Tc: b, Pbar: c > 0 ? c : 1.01325, bg: { Na: a, Cl: a }, model: v.actModel });
    case 'scale-calcite': return 1e6 * mineralSolubility('calcite', { Tc: b, Pbar: Math.max(c, 1.01325), bg: { Na: a, Cl: a }, pCO2: env.pCO2 ?? (v.co2Pct / 100) * Math.max(c, 1.01325), model: v.actModel });
    case 'asphaltene': { const sp = env.spec || DEFAULT_FLUID, pR = num(sp.Pres, 300); return asphalteneOnset({ pre: asphaltenePath({ spec: sp, pRes: pR, tRes: a, pEnd: 5, n: 24, method: 'eos' }), phiA: ((v.saraAsp / 100) * 800) / 1200, deltaA: v.asphDelta, vA: v.asphMV }).upperOnset ?? 0; }
    case 'sand': return sandCriticalVelocity({ d: a * 1e-6, D: b / 1000, rhoP: v.sandRho, rhoF: 1000, mu: (c > 0 ? c : 1) * 1e-3, C: 1e-3, vsl: 1, vm: 1, coef: v.sandCoef, shape: v.sandShape }).governing;
    case 'settling': return 1000 * settlingVelocity(a * 1e-6, b > 0 ? b : v.sandRho, 1000, (c > 0 ? c : 1) * 1e-3, { shape: v.sandShape }).v;
    case 'erosion': return sandErosionScreen(a, b, v.sandUm, v.idMm, 1000, v.erosionSm);
    default: return NaN;
  }
}
const CAL_SAMPLE = [
  {calT: 1, calDT: 4, calInh: 0, calV: 1.5, calDTw: 15, calSI: 1, calSuper: 2, phi: 3.027, dpRatio: 1.054, kLa: 0.006705, dTex: 0.4317, d43Um: 255.1, muRel: 1.17, depMm: 0.04275},
  {calT: 2.5, calDT: 4, calInh: 0, calV: 1.5, calDTw: 15, calSI: 1, calSuper: 2, phi: 5.036, dpRatio: 1.088, kLa: 0.006456, dTex: 0.2867, d43Um: 253.3, muRel: 1.31, depMm: 0.2001},
  {calT: 5, calDT: 6, calInh: 0, calV: 1, calDTw: 15, calSI: 1, calSuper: 2, phi: 7.728, dpRatio: 1.251, kLa: 0.00476, dTex: 0.3434, d43Um: 332.9, muRel: 1.793, depMm: 0.6043},
  {calT: 10, calDT: 8, calInh: 0, calV: 1.5, calDTw: 15, calSI: 1, calSuper: 2, phi: 11.83, dpRatio: 1.646, kLa: 0.006472, dTex: 0.1841, d43Um: 220.5, muRel: 2.102, depMm: 1.99},
  {calT: 2.5, calDT: 10, calInh: 10, calV: 1.5, calDTw: 15, calSI: 1, calSuper: 2, phi: 6.217, dpRatio: 1.201, kLa: 0.006136, dTex: 0.3474, d43Um: 259.4, muRel: 1.462, depMm: 0.2716},
  {calT: 5, calDT: 12, calInh: 20, calV: 2, calDTw: 15, calSI: 1, calSuper: 2, phi: 8.476, dpRatio: 1.243, kLa: 0.00781, dTex: 0.1921, d43Um: 199, muRel: 1.502, depMm: 0.8642},
  {calT: 1, calDT: 6, calInh: 0, calV: 1.5, calDTw: 15, calSI: 1, calSuper: 2, phi: 3.745, dpRatio: 1.156, kLa: 0.006656, dTex: 0.5138, d43Um: 238.3, muRel: 1.206, depMm: 0.05161},
  {calT: 4, calDT: 6, calInh: 0, calV: 2.5, calDTw: 15, calSI: 1, calSuper: 2, phi: 7.059, dpRatio: 1.194, kLa: 0.00905, dTex: 0.1689, d43Um: 159.6, muRel: 1.387, depMm: 0.504},
  {calT: 8, calDT: 5, calInh: 0, calV: 3.2, calDTw: 15, calSI: 1, calSuper: 2, phi: 9.011, dpRatio: 1.288, kLa: 0.01012, dTex: 0.08685, d43Um: 143.8, muRel: 1.489, depMm: 0.8125},
  {calT: 12, calDT: 8, calInh: 0, calV: 1.5, calDTw: 15, calSI: 1, calSuper: 2, phi: 12.32, dpRatio: 1.935, kLa: 0.006101, dTex: 0.171, d43Um: 233.1, muRel: 2.178, depMm: 2.468},
  {calT: 4, calDT: 3, calInh: 0, calV: 1.5, calDTw: 15, calSI: 1, calSuper: 2, teqC: 21.18, tIndH: 0.5268},
  {calT: 4, calDT: 5, calInh: 0, calV: 1.5, calDTw: 15, calSI: 1, calSuper: 2, teqC: 20.81, tIndH: 0.04372},
  {calT: 4, calDT: 8, calInh: 0, calV: 1.5, calDTw: 15, calSI: 1, calSuper: 2, teqC: 21.32, tIndH: 0.003116},
  {calT: 4, calDT: 12, calInh: 15, calV: 1.5, calDTw: 15, calSI: 1, calSuper: 2, teqC: 18.49, tIndH: 0.001809},
  {calT: 4, calDT: 6, calInh: 0, calV: 1.5, calDTw: 15, calSI: 1, calSuper: 1, permM2: 1.029e-11, tDisS: 74.88},
  {calT: 4, calDT: 6, calInh: 0, calV: 1.5, calDTw: 15, calSI: 1, calSuper: 3, permM2: 1.113e-11, tDisS: 17.02},
  {calT: 4, calDT: 6, calInh: 0, calV: 1.5, calDTw: 15, calSI: 1, calSuper: 6, permM2: 1.089e-11, tDisS: 5.909},
  {calT: 6, calDT: 6, calInh: 0, calV: 1.5, calDTw: 10, calSI: 1, calSuper: 2, waxMm: 0.338, waxSolidPct: 1.979},
  {calT: 12, calDT: 6, calInh: 0, calV: 1.5, calDTw: 15, calSI: 1, calSuper: 2, waxMm: 0.8763, waxSolidPct: 2.652},
  {calT: 24, calDT: 6, calInh: 0, calV: 1, calDTw: 15, calSI: 1, calSuper: 2, waxMm: 2.319, waxSolidPct: 2.634},
  {calT: 24, calDT: 6, calInh: 0, calV: 2.5, calDTw: 25, calSI: 1, calSuper: 2, waxMm: 3.087, waxSolidPct: 3.419},
  {calT: 48, calDT: 6, calInh: 0, calV: 3, calDTw: 20, calSI: 1, calSuper: 2, waxMm: 5.682, waxSolidPct: 3.138},
  {calT: 2, calDT: 6, calInh: 0, calV: 1.5, calDTw: 15, calSI: 0.6, calSuper: 2, scaleUm: 9.177, scaleMmY: 0.3099},
  {calT: 6, calDT: 6, calInh: 0, calV: 1.5, calDTw: 15, calSI: 1, calSuper: 2, scaleUm: 123.5, scaleMmY: 1.453},
  {calT: 12, calDT: 6, calInh: 0, calV: 1.5, calDTw: 15, calSI: 1.5, calSuper: 2, scaleUm: 1138, scaleMmY: 501300},
  {calT: 24, calDT: 6, calInh: 0, calV: 1.5, calDTw: 15, calSI: 2, calSuper: 2, scaleUm: 9121, scaleMmY: 146600000000000000000},
  {calT: 4, calDT: 6, calInh: 0, calV: 0.3, calDTw: 15, calSI: 1, calSuper: 2, aopBar: 0, vSetMm: 12.29, vCrit: 0.4079, vResus: 0.166, entrRate: 1.416, erosMmY: 0.00003777},
  {calT: 4, calDT: 6, calInh: 0, calV: 0.8, calDTw: 15, calSI: 1, calSuper: 2, aopBar: 0, vSetMm: 12.49, vCrit: 0.4228, vResus: 0.1774, entrRate: 11.42, erosMmY: 0.0002611},
  {calT: 4, calDT: 6, calInh: 0, calV: 1.5, calDTw: 15, calSI: 1, calSuper: 2, aopBar: 0, vSetMm: 12.76, vCrit: 0.4252, vResus: 0.163, entrRate: 38.44, erosMmY: 0.0009668},
  {calT: 4, calDT: 6, calInh: 0, calV: 3, calDTw: 15, calSI: 1, calSuper: 2, aopBar: 0, vSetMm: 12.86, vCrit: 0.4316, vResus: 0.1702, entrRate: 134.9, erosMmY: 0.003567},
];
const CAL_VALID = [
  {calT: 1.5, calDT: 5, calInh: 0, calV: 2, calDTw: 15, calSI: 1, calSuper: 2, phi: 4.326, dpRatio: 1.099, kLa: 0.008045, dTex: 0.311, d43Um: 205.2, muRel: 1.248, depMm: 0.09689},
  {calT: 4, calDT: 12, calInh: 15, calV: 1.5, calDTw: 15, calSI: 1, calSuper: 2, phi: 8.066, dpRatio: 1.223, kLa: 0.006414, dTex: 0.3045, d43Um: 270, muRel: 1.615, depMm: 0.5742},
  {calT: 2, calDT: 7, calInh: 0, calV: 1.5, calDTw: 15, calSI: 1, calSuper: 2, phi: 5.319, dpRatio: 1.135, kLa: 0.006204, dTex: 0.3831, d43Um: 255.4, muRel: 1.308, depMm: 0.1777},
  {calT: 6, calDT: 4, calInh: 0, calV: 2.8, calDTw: 15, calSI: 1, calSuper: 2, phi: 7.192, dpRatio: 1.217, kLa: 0.009751, dTex: 0.1033, d43Um: 157.3, muRel: 1.395, depMm: 0.7215},
  {calT: 4, calDT: 6, calInh: 0, calV: 1.5, calDTw: 15, calSI: 1, calSuper: 2, teqC: 21.3, tIndH: 0.0127},
  {calT: 4, calDT: 6, calInh: 0, calV: 1.5, calDTw: 15, calSI: 1, calSuper: 4, permM2: 1.076e-11, tDisS: 10.96},
  {calT: 18, calDT: 6, calInh: 0, calV: 1.2, calDTw: 12, calSI: 1, calSuper: 2, waxMm: 1.068, waxSolidPct: 2.223},
  {calT: 36, calDT: 6, calInh: 0, calV: 2, calDTw: 18, calSI: 1, calSuper: 2, waxMm: 3.336, waxSolidPct: 2.964},
  {calT: 8, calDT: 6, calInh: 0, calV: 1.5, calDTw: 15, calSI: 1.2, calSuper: 2, scaleUm: 324.1, scaleMmY: 2.726},
  {calT: 4, calDT: 6, calInh: 0, calV: 2, calDTw: 15, calSI: 1, calSuper: 2, aopBar: 0, vSetMm: 13.63, vCrit: 0.4123, vResus: 0.174, entrRate: 63.96, erosMmY: 0.001572},
];

// =====================================================================================================
// 17a. Reference data sets: the suite's blind prediction for every measured row
// =====================================================================================================
/**
 * Wax deposit thickness (mm) in a cooled single-phase flow loop with the suite's deposition model and default wax
 * parameters. o: { D (m), vel (m/s), Toil, Tcool (°C), tH (h), wat (°C), wTot (mass fraction), mu (Pa·s), rho (kg/m³) }.
 */
export function waxLoopCase({ D, vel, Toil, Tcool, tH, wat, wTot, mu, rho, muOf = null, curve = null, kOil = 0.13, coil = null, n = 120 }, v0 = {}) {
  const v = { ...DEF, ...v0 }, muB = muOf ? muOf(Toil) : mu, Re = (rho * vel * D) / muB, Pr = (2000 * muB) / kOil, hL = (3.66 * kOil) / D,
    hIn = Re < 2300 ? hL : Re > 4000 ? hInside(Re, Pr, kOil, D) : hL + ((hInside(4000, Pr, kOil, D) - hL) * (Re - 2300)) / 1700, // laminar–turbulent transition: interpolated between the two branches
    f = frictionFactor(Re, 4.5e-5 / D), tauW = (f / 8) * rho * vel * vel, U = 1 / (1 / hIn + 1 / 2000);
  let dW = 0, Fw = clamp(1 - (v.waxOil > 0 ? v.waxOil : 80) / 100, 0.03, 0.95), Ti = Tcool;
  for (let k = 0; k < n; k++) { const muI = muOf ? muOf(Ti) : mu; // the diffusivity is evaluated at the deposit surface, the film Reynolds number in the bulk
    const r = waxDeposition({ Tb: Toil, Tamb: Tcool, U, hIn, kOil, rhoOil: rho, muOil: muI, muBulk: muB, wat, wTot, slope: v.waxSlope, delta: dW, Fw, kDep: v.waxK, D, vL: vel, gammaW: tauW / muI, mult: v.waxMult, diffModel: v.waxDiff, stripC: v.waxStripC, aspect: v.waxAspect, curve, wmodel: v.waxModel, c1: v.waxC1, coil: coil ?? (v.waxOil > 0 ? v.waxOil / 100 : null) });
    Ti = r.Ti; dW += (r.dDelta * tH * 3600) / n; Fw = r.Fset !== undefined ? r.Fset : Math.min(Fw + (r.dFw * tH * 3600) / n, 0.95); }
  return dW * 1000;
}
const refMemo = new Map(), memo = (key, f) => { if (!refMemo.has(key)) refMemo.set(key, f()); return refMemo.get(key); };
const methaneT = (Pbar) => memo('mT' + Pbar, () => brent((t) => (vdwpMethane(t).P ?? 1e9) - Pbar, 0, 40, 1e-6, 80)); // °C at which methane hydrate is stable at Pbar (statistical model)
// Tulsa loop crude (Cote Blanche Island): viscosity correlation and wax solubility curve printed in the report of the tests
const CBI = { muOf: (Tc) => { const TF = Tc * 1.8 + 32; return 1e-3 * 10 ** (-1571.5 / (TF * TF) + 130.37 / TF + 0.1977); }, T: [-0.02, 2.47, 4.97, 7.48, 9.98, 12.48, 14.97, 17.47, 19.97, 22.48, 24.98, 27.47, 29.97, 32.47, 34.96, 37.48, 39.97, 42.47, 64.99], C: [0, 0.41, 0.818, 1.21, 1.578, 1.918, 2.228, 2.507, 2.752, 2.943, 3.091, 3.197, 3.27, 3.33, 3.385, 3.426, 3.436, 3.441, 3.441].map((x) => x / 100) };
CBI.curve = { wTot: 0.03441, at: (T) => { const c = interp1(CBI.T, CBI.C, T); return { dissolved: c, solid: 0.03441 - c, dCdT: Math.max(interp1(CBI.T, CBI.C, T + 0.5) - interp1(CBI.T, CBI.C, T - 0.5), 0) }; } };
const phreeqcRow = (r) => memo('ph' + r.w + '|' + r.T + '|' + r.P + '|' + r.cT, () => { const w = PHREEQC_WATERS[r.w - 1], sp = speciate(w, r.T, r.P, { alk: w.HCO3, cT: r.cT }), si = (id) => sp.si(MINERALS.find((q) => q.id === id)); return { pH: sp.pH, cal: si('calcite'), bar: si('barite'), cel: si('celestite'), gyp: si('gypsum'), anh: si('anhydrite') }; });
const rampDefault = () => memo('ramp1', () => onsetRamp({ rate: 1 / 3600, TeqK: 283.5 }));
const REF_MODEL = {
  'barite-nacl-50c': (r) => 1e6 * mineralSolubility('barite', { Tc: 50, bg: { Na: r.mNaCl, Cl: r.mNaCl } }),
  'barite-t-500bar': (r) => 1e6 * mineralSolubility('barite', { Tc: r.T, Pbar: 499.5 }),
  'gypsum-nacl': (r) => 1e3 * mineralSolubility('gypsum', { Tc: r.T, bg: { Na: r.mNaCl, Cl: r.mNaCl } }),
  'celestite-t-p': (r) => 1e6 * mineralSolubility('celestite', { Tc: r.T, Pbar: r.P }),
  'calcite-co2-1bar': (r) => 1e3 * mineralSolubility('calcite', { Tc: r.T, pCO2: 0.9993 }),
  'calcite-pressure': (r) => 40.08e3 * mineralSolubility('calcite', { Tc: 25, Pbar: r.P * 1.01325 }),
  'seawater-calcite-omega': (r) => { const k = r.S / 35.17, w = Object.fromEntries(ION_IDS.map((id) => [id, SEAWATER[id] * k])); w.HCO3 = r.TA * 1e-6 * 61.017 * 1000 * 1.022; return 10 ** scaleIndices(w, r.T, 1.01325, { yCO2: (r.pCO2 * 1e-6) }).minerals[0].SI; },
  'methane-hydrate-equilibrium': (r) => { const m = r.wNaCl > 0 ? ((r.wNaCl / (100 - r.wNaCl)) * 1000) / 58.443 : 0, aw = m > 0 ? speciate({ Na: m, Cl: m }, r.T - KEL, 80).aw : 1; return (vdwpMethane(r.T - KEL, { aw }).P ?? NaN) / 10; },
  'autoclave-conversion': (r) => autoclaveConversion({ tH: r.t, Tc: 1.0, P: 66.5 }),
  'sphere-drag': (r) => (r.Re < 1000 ? (24 / r.Re) * (1 + 0.15 * r.Re ** 0.687) : 0.44),
  'hindered-settling': (r) => 1000 * r.vs * (1 - r.phi) ** settlingVelocity(r.d * 1e-6, 2530, 1240, 0.22).n,
  'sand-minimum-transport': (r) => sandCriticalVelocity({ d: 189e-6, D: r.D, rhoP: 2650, rhoF: r.rhoL, mu: r.mu * 1e-3, C: r.Cv, vsl: 1, vm: 1 }).governing,
  'slurry-critical-velocity': (r) => sandCriticalVelocity({ d: r.d50 * 1e-6, D: 0.078, rhoP: r.rhoP, rhoF: r.rhoF, mu: r.mu * 1e-3, C: r.C, vsl: 1, vm: 1 }).governing,
  'wax-flow-loop': (r) => waxLoopCase({ D: 0.04356, vel: 1.852, Toil: r.Toil, Tcool: r.Tcool, tH: r.t, wat: 35, wTot: 0.03441, rho: 910, muOf: CBI.muOf, curve: CBI.curve }),
  'hydrate-onset-methane-ramp': (r) => rampDefault().quantile(r.p),
  'hydrate-onset-natural-gas-ramp': () => rampDefault().median,
  'barite-induction-time': (r) => Math.log10(scaleInductionTime('barite', { SI: r.SI, TK: r.T + KEL, R: 1 })),
  'phreeqc-barite-si': (r) => phreeqcRow(r).bar, 'phreeqc-celestite-si': (r) => phreeqcRow(r).cel, 'phreeqc-gypsum-si': (r) => phreeqcRow(r).gyp, 'phreeqc-anhydrite-si': (r) => phreeqcRow(r).anh, 'phreeqc-calcite-si': (r) => phreeqcRow(r).cal, 'phreeqc-ph': (r) => phreeqcRow(r).pH,
};
// What each comparison shows. Where the model misses, the miss is stated, not hidden.
const REF_NOTE = {
  'barite-nacl-50c': 'Pitzer model with the pitzer.dat parameters, nothing fitted here. The salting-in by a factor of fourteen between pure water and 5 mol/kg NaCl is reproduced within about 5 %; the same residual is seen in the PHREEQC test plot (digitised points).',
  'barite-t-500bar': 'Tests the temperature function of the solubility product and the pressure term (molar volumes of the ions and of the solid) together, 493 atm.',
  'gypsum-nacl': 'Six isotherms from 0.5 to 95 °C up to halite saturation: the solubility maximum near 3 mol/kg and its decline at high salinity come from the Pitzer mixing terms.',
  'celestite-t-p': 'Pure water, pressure up to 600 bar. The 250 °C point of the source is outside the range used here and is left out.',
  'calcite-co2-1bar': 'Open system at a CO2 fugacity of 0.999 bar. The scatter between the laboratories in the source is itself about 10 %.',
  'calcite-pressure': 'Closed system (no gas phase): calcite dissolving in pure water under pressure; tests the pressure dependence of the carbonate equilibria.',
  'seawater-calcite-omega': 'The model is given salinity, total alkalinity and CO2 partial pressure and returns Ω = 10^SI. It is biased high by about 1 unit of Ω because the whole measured alkalinity is treated as carbonate alkalinity (seawater holds roughly 4 % of it as borate) and because the reference values were themselves calculated with seawater-scale constants; this is the accuracy to expect when a water analysis gives only a total alkalinity.',
  'methane-hydrate-equilibrium': 'van der Waals–Platteeuw model with the Parrish–Prausnitz constants; for the NaCl series the water activity comes from the Pitzer model. Pure methane, structure I only.',
  'autoclave-conversion': 'CALIBRATION SET for the shell diffusivity. Shrinking-core conversion of 40 µm droplets at the driving force of the test, advanced exactly through the growing shell; the default shell diffusivity (0.8 × 10⁻¹⁵ m²/s) was chosen on this test, inside the range between annealed films (10⁻¹⁷–10⁻¹⁶) and fresh films (3–8 × 10⁻¹³ m²/s). The two-hour lag before bulk conversion starts (0.6 % at 2 h) is not modelled, which is most of the remaining error. Earlier versions stepped the shell growth explicitly and converted all the water within about two hours. No independent open conversion–time series was found for a blind check; flow-loop conversions (about 15 % in one hour, 60 % in a day with the same parameters) are in the range reported for loops.',
  'sphere-drag': 'Schiller–Naumann drag law with the constant 0.44 above Re = 1000; the measured coefficients between Re = 10⁴ and 10⁵ scatter from 0.39 to 0.64. All 44 rows were checked against the page image of the report.',
  'hindered-settling': 'Richardson–Zaki exponent 4.65 at low Reynolds number applied to the Stokes speed printed in the source (the authors fit 4.48 ± 0.04). Eight of the published points.',
  'sand-minimum-transport': 'Regime-selected criterion, nothing fitted: turbulent suspension in water (grain at the edge of the viscous sub-layer), the Thomas sliding-bed limit in the 7 and 20 mPa·s polymer solutions (grain inside the sub-layer; over-predicted by 20–55 %, the solutions are shear-thinning), and the mobile bed-load limit at a wall Shields number of 1.5 in the laminar oil tests (105–340 mPa·s), where the turbulent correlations used before do not apply. The Shields value is the upper end of the range of the viscous bed-load theory, not a fit, but these oil rows are the only check of it.',
  'slurry-critical-velocity': 'BLIND. The same regime selection applied to an independent public data set: fine dense particles in water and in kaolin carriers (grain inside the viscous sub-layer for all but two rows, so the Thomas limit governs). The kaolin carriers are Bingham fluids; only their consistency viscosity is used and the yield stress is ignored.',
  'wax-flow-loop': 'Matzain model with its published constants, the viscosity correlation and the wax solubility curve printed in the report, nothing fitted. The comparison is with the thickness after solvent wash (the wax deposit); the thickness as recovered includes a gelled-oil layer of about 0.25 mm that is present after three hours and holds about 93 % oil, which a deposition model does not describe. Earlier versions used a viscosity five times too low and compared with the unwashed thickness. The tests are laminar to transitional (Re about 2 200–2 900): the fully developed laminar film coefficient under-predicts the 12–24 h points by a factor of two to three, and the 35 °C test, which the kernel correlation treats as transitional, is over-predicted.',
  'hydrate-onset-methane-ramp': 'BLIND. Measured-rate nucleation model with the rate constants published for a stirred 10.9 cm³ cell (isothermal induction times at 2–4 K and fast ramps; those experiments are the calibration data, fitted by their authors) applied without change to rocking cells of another laboratory. Rows are the 18 measured onsets in ascending order against the predicted quantiles of the onset distribution at 1 K/h (predicted mean 3.1 K, standard deviation 0.7 K; measured 2.3 ± 0.7 K with the hydrate curve of this suite, 2.65 ± 0.68 K as evaluated by the authors). The interface of the rocking cells is taken equal to that of the calibration cell. Classical nucleation theory with bulk properties, used before, gave onset near 14 K.',
  'hydrate-onset-natural-gas-ramp': 'BLIND. Same model and constants (measured for methane) applied to a structure-II natural gas and to methane in steel rocking cells of a third laboratory; each row is the mean of 5–10 cells. The subcooling uses the equilibrium temperature at the filling pressure, so the true values are a few tenths of a kelvin smaller. Stirred-cell measurements with a structure-II gas at the source laboratory of the constants nucleate 2–3 K later than these rocking cells: the rate multiplier is apparatus-dependent.',
  'barite-induction-time': 'Semi-empirical induction-time relation of the Rice University Brine Chemistry Consortium as published (open access); the relation was fitted by its authors to a larger data base that includes measurements of this kind, so this is a check of the implementation and of the scatter (± 0.4 log units), not an independent validation.',
  'phreeqc-barite-si': 'Blind code-to-code comparison: same molalities, alkalinity, total carbon, temperature and pressure in both codes. Agreement to about 0.002 units over the whole matrix.',
  'phreeqc-celestite-si': 'Blind code-to-code comparison, as for barite.', 'phreeqc-gypsum-si': 'Blind code-to-code comparison, as for barite.', 'phreeqc-anhydrite-si': 'Blind code-to-code comparison, as for barite.',
  'phreeqc-calcite-si': 'Blind code-to-code comparison. Agreement within 0.005 units up to 60 °C. At 100–150 °C and low CO2 fugacity (pH above about 7) the engine is up to 0.6 units high because it does not carry the MgOH⁺ species of pitzer.dat, which takes up hydroxide in hot magnesium-bearing brines; produced waters under CO2 pressure are below that pH.',
  'phreeqc-ph': 'Blind code-to-code comparison; same remark as for calcite: the differences are confined to hot (≥ 100 °C), high-pH states where MgOH⁺ matters.',
};
const REF_TOL = { 'barite-nacl-50c': { mape: 7 }, 'barite-t-500bar': { mape: 6 }, 'gypsum-nacl': { mape: 3 }, 'celestite-t-p': { mape: 7 }, 'calcite-co2-1bar': { mape: 15 }, 'calcite-pressure': { mape: 5 }, 'seawater-calcite-omega': { mape: 30 }, 'methane-hydrate-equilibrium': { mape: 8 }, 'sphere-drag': { mape: 12 }, 'hindered-settling': { mape: 12 },
  'autoclave-conversion': { rmse: 8 }, 'sand-minimum-transport': { mape: 30 }, 'slurry-critical-velocity': { mape: 30 }, 'wax-flow-loop': { bias: 0.15, rmse: 0.35 }, 'hydrate-onset-methane-ramp': { mape: 45, bias: 1 }, 'hydrate-onset-natural-gas-ramp': { mape: 25 }, 'barite-induction-time': { rmse: 0.6 },
  'phreeqc-barite-si': { maxAbs: 0.01 }, 'phreeqc-celestite-si': { maxAbs: 0.01 }, 'phreeqc-gypsum-si': { maxAbs: 0.01 }, 'phreeqc-anhydrite-si': { maxAbs: 0.01 }, 'phreeqc-calcite-si': { rmse: 0.12 }, 'phreeqc-ph': { rmse: 0.12 } };
const validationData = REF_SETS.map((s) => ({ ...s, model: REF_MODEL[s.id], ...(REF_TOL[s.id] ? { tolerance: REF_TOL[s.id] } : {}), note: REF_NOTE[s.id] || s.note }));
export const VALIDATION_DATA = validationData;

// =====================================================================================================
// 17b. Provenance of the literature constants
// =====================================================================================================
const PH_DB = 'https://raw.githubusercontent.com/usgs-coupled/phreeqc3/master/database/', PH_SRC = 'https://raw.githubusercontent.com/usgs-coupled/phreeqc3/master/src/', RD = '2026-10-08', RD2 = '2026-10-09';
const prov = (item, used, source, url, status, note) => ({ item, used, source, url, retrieved: RD, status, note });
/** Every literature constant set the engine relies on, with the address it was checked against on the retrieval date. */
export const PROVENANCE = [
  prov('Mineral solubility products with temperature functions: calcite, barite, celestite, gypsum, anhydrite (analytic log K) and solid molar volumes', 'speciate(), scaleIndices(), mineralSolubility()', 'USGS PHREEQC version 3 database pitzer.dat (PHASES)', PH_DB + 'pitzer.dat', 'corrected', 'The earlier 25 °C constants with van \'t Hoff slopes are replaced by the database expressions, extracted by script. Largest change at 25 °C: barite log K −9.97 → −9.91 (analytic expression), anhydrite −4.36 → −4.25 at 25 °C and −4.52 → −4.90 at 70 °C; reaction volumes −42 … −58 cm³/mol are replaced by pressure- and ionic-strength-dependent ion volumes.'),
  prov('Siderite solubility product (log K −10.89, ΔH −2.48 kcal/mol)', 'speciate()', 'USGS PHREEQC version 3 database phreeqc.dat (pitzer.dat has no siderite)', PH_DB + 'phreeqc.dat', 'verified', 'log K identical; ΔH −10.38 kJ/mol identical.'),
  prov('Pitzer interaction parameters β⁰, β¹, β², C⁰, θ, ψ, λ, ζ with temperature coefficients for Na-K-Ca-Mg-Ba-Sr-Fe-H-Cl-SO4-HCO3-CO3-OH-CO2', 'pitzerGamma()', 'USGS PHREEQC version 3 database pitzer.dat (PITZER block; Appelo 2015 parameter set)', PH_DB + 'pitzer.dat', 'verified', 'New in this version (extracted by script, 125 parameter lines). Check: barite solubility at 50 °C in 0–4 mol/kg NaCl reproduces the PHREEQC output file mytest/Barite_NaCl.out to four significant figures (largest difference 0.02 %).'),
  prov('Pitzer model equations: α values by charge type, temperature function, higher-order electrostatic terms (Chebyshev coefficients), MacInnes scaling, pressure correction of the Debye–Hückel term', 'pitzerGamma()', 'USGS PHREEQC version 3 source, pitzer.cpp', PH_SRC + 'pitzer.cpp', 'verified', 'Ported line by line; γ± of 1 mol/kg NaCl = 0.657 and water activity 0.9669 against tabulated values.'),
  prov('Water density, dielectric constant (Bradley–Pitzer), Debye–Hückel slopes, Born function; aqueous molar-volume equation', 'waterDH(), speciesVm()', 'USGS PHREEQC version 3 source, utilities.cpp, prep.cpp and read.cpp', PH_SRC + 'utilities.cpp', 'corrected', 'The earlier polynomial dielectric constant and density are replaced. Debye–Hückel A at 25 °C: 0.509 → 0.510.'),
  prov('Ion-size and salting-out parameters (Truesdell–Jones a, b) of the ion-association model', 'speciate() with model "truesdellJones"', 'USGS PHREEQC version 3 database phreeqc.dat (-gamma)', PH_DB + 'phreeqc.dat', 'corrected', 'Ba²⁺: a 5.0 → 7.1 Å, b 0 → 0.0553; Na⁺: 4.0/0.075 → 4.08/0.082; Cl⁻: 3.5/0.015 → 3.63/0.017; Ca, Mg, Sr, K, SO4 identical.'),
  prov('Ion-pair association constants and their temperature functions (CaSO4°, MgSO4°, NaSO4⁻, Na2SO4°, KSO4⁻, SrSO4°, FeSO4°, CaHCO3⁺, CaCO3°, MgHCO3⁺, MgCO3°, NaHCO3°, KHCO3°, SrHCO3⁺, SrCO3°, BaHCO3⁺, BaCO3°, FeHCO3⁺, FeCO3°, FeCl⁺)', 'speciate() with model "truesdellJones" or "davies"', 'USGS PHREEQC version 3 database phreeqc.dat (SOLUTION_SPECIES)', PH_DB + 'phreeqc.dat', 'corrected', 'log K at 25 °C, old → database: CaSO4° 2.30 → 2.14, MgSO4° 2.37 → 2.42, NaSO4⁻ 0.70 → 0.94, KSO4⁻ 0.85 → 1.18, NaHCO3° −0.25 → −0.43; BaSO4° 2.70 removed (not in the database). Carbonate pairs and the temperature dependence are new.'),
  prov('Carbonate system: CO2 solubility, first and second dissociation of carbonic acid, water ionisation (analytic log K)', 'speciate()', 'USGS PHREEQC version 3 databases pitzer.dat and phreeqc.dat (Plummer–Busenberg expressions)', PH_DB + 'pitzer.dat', 'verified', 'The Plummer–Busenberg coefficients used before agree with the database to the printed digits; the CO2 fugacity coefficient now comes from the kernel Peng–Robinson model instead of the Oddo–Tomson exponential.'),
  prov('Calcite dissolution rate constants (Plummer–Wigley–Parkhurst k1, k2, k3)', 'not used by the engine (listed for reference)', 'USGS PHREEQC version 3 database phreeqc.dat (RATES)', PH_DB + 'phreeqc.dat', 'verified', 'Read but not adopted: the scale surface-reaction constant remains a calibration input.'),
  { ...prov('Code-to-code check of the aqueous model (replaces the Oddo–Tomson cross-check, whose gas-phase coefficients could not be verified and which has been removed)', 'speciate(), scaleIndices()', 'USGS PHREEQC 3.8.6 with pitzer.dat, run for this suite (216 states, 4–150 °C, 1–600 bar)', 'https://water.usgs.gov/water-resources/software/PHREEQC/iphreeqc-3.8.6-17100.tar.gz', 'verified', 'Sulphate saturation indices agree within 0.002 units over the whole matrix; pH and calcite within 0.005 up to 60 °C and up to 0.6 at 100–150 °C and pH above 7, where the MgOH⁺ species of pitzer.dat (not carried here) matters.'), retrieved: RD2 },
  { ...prov('Barite induction time (1.5232, −10.8784, −895.6683, 5476.992, 0.8286, 0.225; above 90 °C −2.11, −4.29, 279.29, 3332.26, 0.99) and threshold-inhibitor coefficients (phosphonates −1.41, 1329.29, 0.15, 0.18; polycarboxylates −1.28, 1007.30, 0.03, 0.13; sulfonated −1.72, 1229.65, −0.01, 0.20; product constants)', 'scaleInductionTime(), scaleInhibition(), SCALE_INHIBITORS', 'Dai, C. et al. (2021) Sustainability 13, 8533, Eq. 11–12 and Table 2 (CC BY 4.0); minimum-inhibitor-concentration form from the Rice University Brine Chemistry Consortium course notes', 'https://mdpi-res.com/d_attachment/sustainability/sustainability-13-08533/article_deploy/sustainability-13-08533.pdf', 'verified', 'All coefficients identical. Against 104 measured induction times of a Rice University thesis the relation gives a mean error of +0.12 and a scatter of 0.39 in log10 t; the source quotes minimum inhibitor concentrations of 0.15–0.9 mg/L for barite at SI 1.1–1.4 and 100 °C, the implementation gives 0.35 mg/L at SI 1.08 for 24 h.'), retrieved: RD2 },
  { ...prov('Celestite induction time (−1.324, −3.301, 2462.5) and inhibition constants (DTPMP −13.41, 3526.7, 7.36; PPCA −8.24, 1946.5, 6.67; PVS −6.82, 1898.0, 3.30)', 'scaleInductionTime(), scaleInhibition()', 'Zhao, Y. (2021) PhD thesis, Rice University', 'https://repository.rice.edu/server/api/core/bitstreams/214c4ca7-7d77-4096-8491-e0c9019d3f90/content', 'verified', 'Identical (thesis, constants quoted with citation).'), retrieved: RD2 },
  { ...prov('Scale surface-reaction constant and activation energy (barite, neutral mechanism: log k = −7.90 mol/m²/s at 25 °C, 30.8 kJ/mol); crystal–water interfacial energy', 'slowDeposits(), scaleCrystallisation(), inputs scaleK, scaleSigma', 'Palandri, J.L. and Kharaka, Y.K. (2004) U.S. Geological Survey Open-File Report 2004-1068, Table 34; interfacial-energy fits quoted in Dai et al. (2021)', 'https://pubs.usgs.gov/of/2004/1068/pdf/OFR_2004_1068.pdf', 'corrected', 'Rate constant 1.0 → 1.26 × 10⁻⁸ mol/m²/s and activation energy 45 → 30.8 kJ/mol; interfacial energy 90 → 79 mJ/m² (published fits 75–79 for barite, literature range 38–150).'), retrieved: RD2 },
  { ...prov('Sulphate left by a sulphate-removal plant, 40 mg/L', 'input srpSO4', 'Offshore magazine, Heidrun sulphate-removal plant ("guaranteed quality of 40 mg SO4/liter")', 'https://www.offshore-mag.com/business-briefs/equipment-engineering/article/16756906/new-system-at-heidrun-removes-sulfate-from-seawater', 'verified', 'Single quoted value; raw seawater 2 712 mg/L here against 2 790 mg/L quoted in Dai et al. (2021).'), retrieved: RD2 },
  prov('Bulk crystallisation: nucleation pre-exponential 10³⁰ 1/m³/s, parabolic growth constant 10⁻¹⁰ m/s, crystal sticking efficiency, seed concentration', 'scaleCrystallisation(), slowDeposits()', 'none', 'https://pubs.usgs.gov/of/2004/1068/pdf/OFR_2004_1068.pdf', 'unverified', 'Order-of-magnitude model inputs and calibration parameters; with the induction-time gate they act only once a brine has passed its induction time.'),
  prov('Kim–Bishnoi dissociation: K0 = 3.6e4 mol/m²/Pa/s, E = 81 kJ/mol (Clarke and Bishnoi 2001)', 'hydrateDissociationRate()', 'OSTI abstract of Clarke and Bishnoi (2001), Can. J. Chem. Eng. 79, 143', 'https://www.osti.gov/etdeweb/biblio/20155246', 'verified', 'Both numbers identical. Flow-loop fits on natural-gas hydrate slurry (Lv et al. 2021, RSC Advances 11) give apparent values of 56–65 kJ/mol and 4.8–6.8e4: a lumped slurry constant, not the intrinsic one, so the two are not compared as a data set any more.'),
  prov('Intrinsic formation rate constant of methane hydrate', 'hydrateGrowthRate(), input kinK', 'Englezos (1986) and Malegaonkar (1996) theses, University of Calgary', 'https://ucalgary.scholaris.ca/server/api/core/bitstreams/0f91465c-e6d9-41fb-b506-40daccd19cf5/content', 'corrected', 'Default 1.0e-10 → 0.25e-10 mol/m²/Pa/s: Englezos et al. give 0.55–0.65e-5 mol/m²/MPa/s (274–282 K), the later corrected values are 0.21–0.31e-4. The reference case is shell-diffusion limited and does not change.'),
  prov('Temperature coefficient of the formation constant, E/R = 13 600 K; rate scaling 1/500', 'input kinEa', 'Qu et al. (2024) Energies 17, 6101; Zerpa (2013) and Boxall (2009) theses, Colorado School of Mines', 'https://mdpi-res.com/d_attachment/energies/energies-17-06101/article_deploy/energies-17-06101.pdf', 'verified', 'k2 = −13 600 K identical (it belongs to the subcooling-driven CSMHyK law and is used here as the Arrhenius temperature). The 1/500 scaling of that model is not used; fitted values range from 0.0006 to 0.02.'),
  prov('Hydrate density, heat capacity, thermal conductivity', 'HYDRATE, input rhoHyd', 'Zerpa (2013) thesis quoting Sloan and Koh (2008)', 'https://repository.mines.edu/server/api/core/bitstreams/3108bae4-8282-456b-a0be-cd19dab3b929/content', 'corrected', 'Density 920 → 914 kg/m³ (structure II; structure I 910), heat capacity 2100 → 2080 J/kg/K, conductivity 0.6 → 0.49 W/m/K (structure I).'),
  prov('Latent heat of hydrate formation', 'HYDRATE.latent', 'Davies (2009) thesis quoting Handa (1986); Qu et al. (2024) for simulator defaults', 'https://repository.mines.edu/server/api/core/bitstreams/b3c44b83-0c2e-449c-8ff7-3fea1c9359fe/content', 'corrected', '4.4e5 → 4.6e5 J/kg (460 J/g for structure I; simulator defaults 469–477 kJ/kg). The 54.2 kJ/mol of Handa was not seen in an open source.'),
  prov('Hydration number (5.75 ideal, about 6 measured)', 'langmuirOccupancy(), input hydNumber', 'Qu et al. (2024) Energies 17, 6101; Boxall (2009) Table D.2', 'https://mdpi-res.com/d_attachment/energies/energies-17-06101/article_deploy/energies-17-06101.pdf', 'verified', 'Ideal structure I 5.75; 6.089 used for methane in the autoclave evaluation; the engine computes 5.8–6.1 from the occupancies.'),
  prov('Parrish–Prausnitz Langmuir constants of methane in structure I (A, B small cage 3.7237e-3, 2708.8; large cage 1.8372e-2, 2737.9)', 'langmuirOccupancy()', 'Zhang et al. (2022) RSC Advances 12 (table of the Parrish–Prausnitz constants)', 'https://pmc.ncbi.nlm.nih.gov/articles/PMC9133728/', 'verified', 'All four numbers identical. The pressure unit (1/atm) is not printed in the open table; with it the model gives 26.6 bara at 0 °C against 25.6–26.3 measured.'),
  prov('Empty-lattice reference properties of structure I: Δμ⁰ = 1264 J/mol, Δh⁰ = −4858 J/mol, Δcp, Δv = 4.6 cm³/mol', 'emptyLatticeDmu()', 'open-source teaching notebook (PyTherm, methane hydrates)', 'https://raw.githubusercontent.com/iurisegtovich/PyTherm-applied-thermodynamics/master/contents/main-lectures/HYD1-methane-hydrates.ipynb', 'corrected', 'Δμ⁰, Δh⁰ and Δv identical; Δcp −38.12 + 0.141 (T − 273.15) → constant −39.16 J/mol/K (the only value seen in an open source). Effect on the equilibrium pressure below 0.3 %.'),
  prov('Kihara parameters of methane in hydrate: a = 0.3834 Å, σ = 3.1650 Å, ε/k = 154.54 K', 'kiharaLangmuir()', 'Zhang et al. (2022) RSC Advances 12', 'https://pmc.ncbi.nlm.nih.gov/articles/PMC9133728/', 'verified', 'Identical.'),
  { ...prov('Structure-I cage radii and coordination numbers for the single-shell cell potential (3.95 Å, 20; 4.33 Å, 24)', 'kiharaLangmuir()', 'Scientific Reports 7 (2017), PMC5526936, and npj Clean Water (2025), PMC12158787, Table 1 (radii, CC BY); ACS Omega 7 (2022), PMC8991894, Table 2 (radii with coordination numbers)', 'https://pmc.ncbi.nlm.nih.gov/articles/PMC5526936/', 'verified', 'Radii identical in two open sources; coordination numbers identical in the third.'), retrieved: RD2 },
  prov('Hydrate–water interfacial energy', 'nucleationRate(), input sigmaHW', 'Chemical Reviews 125 (2025) 5003, quoting the measurement of Anderson et al.', 'https://pmc.ncbi.nlm.nih.gov/articles/PMC12123632/', 'corrected', '20 → 32 mJ/m² (experimental value; simulations give 38–44). The 20 mJ/m² attributed to Kashchiev and Firoozabadi was not found in an open source.'),
  { ...prov('Heterogeneous nucleation rate of methane hydrate, J = A exp(Δs ΔT/kT) exp(−B′/(T ΔT²)): A = 0.58 × 10⁻³ 1/s, B′ = 2.66 × 10³ K³ (isothermal, 2–4 K) and A = 0.10 1/s, B′ = 3.1 × 10⁴ K³ (ramps, above 5 K) per 10.9 cm³ stirred cell; Δs = 22.2 k; gas–water interface of the cell 11 cm²', 'hydrateNucleationRate(), onsetRamp(), NUCLEATION', 'Lim et al. (2021) Chem. Eng. J. 411, 128478, Table S4; Barwood et al. (2022) Chem. Eng. J.; Li (2025) PhD thesis, University of Western Australia (accepted manuscripts and thesis in the university repository)', 'https://api.research-repository.uwa.edu.au/ws/files/189542301/2022_Barwood_Extracting_J_Constant_Cooling.pdf', 'verified', 'Constants as published (calibration by their authors on about 2 000 isothermal and 860 ramped formation events). Check of the implementation: mean onset 6.4 K at 1 K/min and 7.6 K at 3 K/min against 6.60 ± 1.3 K and 7.28 ± 0.96 K measured in the source. Replaces the classical-theory defaults (pre-exponential 10^7.5 1/m³/s, contact angle 40°), which remain as a comparison option and nucleate near 14 K.'), retrieved: RD2 },
  prov('Droplet size in water-in-oil flow: d/D = 0.063 We^−3/5', 'hydrateParams()', 'Boxall (2009) thesis, Colorado School of Mines', 'https://repository.mines.edu/server/api/core/bitstreams/443d7dac-3295-4035-be9b-7392451cd688/content', 'verified', 'Fitted prefactor 0.0628. The viscous branch (0.016 Re^½ We⁻¹) is not implemented.'),
  prov('Cohesive force between hydrate particles, 2 mN/m', 'input cohesion', 'Zerpa (2013) and Dieker (2009) theses; Qu et al. (2024)', 'https://repository.mines.edu/server/api/core/bitstreams/3108bae4-8282-456b-a0be-cd19dab3b929/content', 'verified', 'Inside the measured range: 4.3 mN/m for cyclopentane hydrate, about 0.5 with crude oil, 1.6 back-calculated by Camargo and Palermo; simulator defaults are 27–50 mN/m.'),
  prov('Camargo–Palermo agglomerate balance (fractal dimension 2.5, packing 4/7) and Mills viscosity (1 − φ)/(1 − φ/φmax)²', 'maxAgglomerateSize(), slurryViscosity()', 'Zerpa (2013) thesis; Qu et al. (2024) Energies 17, 6101', 'https://mdpi-res.com/d_attachment/energies/energies-17-06101/article_deploy/energies-17-06101.pdf', 'verified', 'Equation form and both constants identical.'),
  prov('Hayduk–Minhas (13.3e-12, 1.47, 10.2/V − 0.791, 0.71) and Wilke–Chang (7.4e-12, 0.6) diffusivities', 'waxDiffusivity()', 'Sarica and Volk (2004) Tulsa University Paraffin Deposition Projects, U.S. DOE final report', 'https://www.netl.doe.gov/sites/default/files/2018-05/BC15150_Final.pdf', 'verified', 'All coefficients identical in SI form.'),
  prov('Hindered diffusivity in the wax gel, D/(1 + α²F²/(1 − F)) (Singh et al. 2000, Cussler)', 'waxDeposition()', 'Giraldo et al. (2015) CT&F 6(1), reproducing Singh et al.', 'http://www.scielo.org.co/scielo.php?pid=S0122-53832015000100003&script=sci_arttext', 'verified', 'Form identical. The aspect ratio (8) is a model input, not a sourced constant.'),
  prov('Matzain shear-stripping constants C2 = 0.055, C3 = 1.4 and flow-pattern forms of the stripping number', 'waxDeposition()', 'Aalborg University thesis quoting Matzain; equation form in the Tulsa/DOE report', 'https://projekter.aau.dk/projekter/files/335444794/K10_OG_2_F20.pdf', 'verified', 'Identical (secondary source).'),
  { ...prov('Won melting temperature and heat of fusion of n-paraffins (374.5, 0.02617, 20172; 0.1426)', 'waxEquilibrium() through the fluid suite’s wax model', 'PetroWiki, Thermodynamic models for wax precipitation (archived copy; equations printed as images)', 'https://web.archive.org/web/2id_/https://petrowiki.spe.org/Thermodynamic_models_for_wax_precipitation', 'verified', 'All four constants identical.'), retrieved: RD2 },
  { ...prov('Matzain wax model: enhancement constant C1 = 15, deposit oil fraction 100 (1 − Re^0.15/8); wax solubility slope 0.04 1/K', 'waxDeposition(), inputs waxC1, waxSlope', 'Aalborg University thesis quoting Matzain; solubility curve of the Cote Blanche Island crude in Sarica and Volk (2004), Fig. 59', 'https://projekter.aau.dk/projekter/files/335444794/K10_OG_2_F20.pdf', 'verified', 'C1 and the closure identical. The measured solubility curve has d ln C/dT = 0.068 1/K at 12 °C and 0.016 1/K at 25 °C, bracketing the default 0.04. The sources disagree on the viscosity unit of the stripping number; mPa·s is used. The measured deposit of that crude held about 93 % oil against 60 % from the closure.'), retrieved: RD2 },
  prov('Mechanistic wax option: shear-dispersion coefficient 0.1, crystal size 10 µm, crystal aspect ratio 8; gel yield-stress coefficient', 'waxDeposition() with model "mechanistic", slowDeposits()', 'none', 'https://www.netl.doe.gov/sites/default/files/2018-05/BC15150_Final.pdf', 'unverified', 'Not the default model; inputs to be fitted with the deposition multiplier. No open gel-restart data set was found for the yield-stress coefficient.'),
  prov('Schiller–Naumann drag (0.15, 0.687, 0.44 above Re 1000)', 'settlingVelocity(), maxeyRileyStepper()', 'OpenFOAM 7 source, SchillerNaumann.C', 'https://raw.githubusercontent.com/OpenFOAM/OpenFOAM-7/master/applications/solvers/multiphase/reactingEulerFoam/interfacialModels/dragModels/SchillerNaumann/SchillerNaumann.C', 'verified', 'Identical.'),
  prov('Richardson–Zaki exponent (4.65; 4.4 Re^−0.03; 4.4 Re^−0.1; 2.4)', 'settlingVelocity()', 'arXiv:2008.07137; alternative set in arXiv:1711.00336', 'https://arxiv.org/pdf/2008.07137', 'verified', 'Identical to the first source; the more common variant is 4.35 Re^−0.03, 4.45 Re^−0.1, 2.39 (largest difference 1.1 %).'),
  prov('Suspension viscosity: Einstein 2.5, Krieger–Dougherty exponent −2.5 φmax, Thomas 1 + 2.5φ + 10.05φ² + 0.00273 exp(16.6φ)', 'slurryViscosity()', 'arXiv:1207.3774; OpenFOAM 7 source slurry.C', 'https://raw.githubusercontent.com/OpenFOAM/OpenFOAM-7/master/applications/solvers/multiphase/driftFluxFoam/mixtureViscosityModels/slurry/slurry.C', 'verified', 'Identical.'),
  prov('Oroskar–Turian critical velocity (1.85, 0.1536, 0.3564, 0.378, 0.09, x^0.30)', 'sandCriticalVelocity()', 'Miedema, Slurry Transport (LibreTexts); PNNL-17639 for the factor x', 'https://eng.libretexts.org/Bookshelves/Civil_Engineering/Slurry_Transport_(Miedema)/06%3A_Slurry_Transport_a_Historical_Overview/6.27%3A_The_Limit_Deposit_Velocity_(LDV)', 'corrected', 'Coefficient and exponents identical; x 0.96 → 0.95 (the value of the design guide read in PNNL-17639; effect 0.3 %).'),
  prov('Salama minimum velocity (exponents 0.53, 0.17, −0.09, 0.55, 0.47) and Danielson (K = 0.23, exponents −1/9, 1/9, 5/9)', 'sandCriticalVelocity()', 'Yan (2010) thesis, Cranfield University; Bello (2013) thesis, Robert Gordon University', 'https://core.ac.uk/reader/19209312', 'verified', 'Exponents and K identical; the unit coefficient of the Salama form was checked by reproducing the predictions tabulated by Yan (0.36 m/s for 200 µm sand in a 0.1 m pipe).'),
  prov('Salama erosion relation, geometry constant 5.5', 'sandErosionScreen(), input erosionSm', 'Salama (2000)', 'https://www.slideshare.net/slideshow/an-alternative-to-api-14-e-erosional-velocity-limits-for-sand-laden-fluids-mamdouh-m-salama-conoco-inc-2000-7-pgs/222775017', 'unverified', 'Form and units confirmed from a transcript; the constant 5.5 was not seen. Screening only.'),
  prov('Kozeny–Carman constant 180 and Ergun constants 150, 1.75', 'kozenyCarman(), porousGradient()', 'arXiv:1710.09314; OpenFOAM 7 source Ergun.C', 'https://raw.githubusercontent.com/OpenFOAM/OpenFOAM-7/master/applications/solvers/multiphase/reactingEulerFoam/interfacialModels/dragModels/Ergun/Ergun.C', 'verified', 'Identical.'),
  prov('Turbulent deposition velocity: 0.057 Sc^−2/3 + 4.5e-4 τ⁺², plateau', 'depositionVelocity()', 'Clarkson University ME637 lecture notes on Wood (1981)', 'https://webspace.clarkson.edu/projects/crcd/public_html/me637/downloads/P_TurbDep.pdf', 'corrected', 'Both regime constants identical; plateau 0.14 → 0.13.'),
  prov('Collision kernels: shear (G/6)(di + dj)³, Saffman–Turner 0.1618 √(ε/ν)(di + dj)³, Brownian 2kT/(3μ)', 'aggregationKernel()', 'arXiv:2510.05319, arXiv:2506.03881, arXiv:1909.10608', 'https://arxiv.org/pdf/2510.05319', 'verified', 'Identical after conversion between radius and diameter forms (√(8π/15) = 1.2944 on radii).'),
  prov('Maxey–Riley equation: added-mass coefficient ½, Basset kernel (3/2) d² √(πρμ), settling benchmark', 'maxeyRileyStepper()', 'arXiv:2006.16577; Prasath, Vasan and Govindarajan, arXiv:1808.08769', 'https://arxiv.org/pdf/1808.08769', 'verified', 'Form identical; the long-time approach to the terminal velocity ∝ t^−½ is reproduced within 0.1 %.'),
  prov('Axial dispersion in turbulent pipe flow, 10.1 a u*', 'run(): Eulerian solids transport', 'Hart (2013) thesis, University of Warwick (Taylor 1954)', 'https://wrap.warwick.ac.uk/id/eprint/57725/1/WRAP_THESIS_Hart_2013.pdf', 'verified', 'Identical; valid above Re = 20 000.'),
  { ...prov('Eddy-viscosity profile of pipe flow (κ/6)(1 − r²)(1 + 2r²) u* R (Reichardt) and the composite with wall damping (Cess), κ = 0.42, A⁺ = 27; turbulent Schmidt number 0.9, Cμ = 0.09', 'trackParticles(), pipeFlowField(), RANS', 'arXiv:2604.09454, Eq. 7 and 12; arXiv:2008.13486, Eq. 35, and arXiv:2302.14408, Eq. 19 (CC BY 4.0), constants fitted to the Princeton pipe data', 'https://arxiv.org/abs/2302.14408', 'verified', 'Profile and constants identical. The resolved profile gives a friction factor within 1.5 % of the smooth-pipe law from Re = 10⁵ to 3 × 10⁶ (8 % high at Re = 5 000). The particle tracker keeps κ = 0.4.'), retrieved: RD2 },
  { ...prov('Hinze bubble size 0.725 (σ/ρ)^0.6 ε^−0.4; small-eddy mass-transfer coefficient 0.4 (εν)^¼ Sc^−½', 'slugMassTransfer(), hydrateMarch()', 'arXiv:2011.00963 (Hinze constant); Farsoiya et al. (2023) J. Fluid Mech. (CC BY) and Katul et al. (2024) for the surface-renewal coefficient', 'https://arxiv.org/abs/2011.00963', 'verified', 'Hinze constant identical; the mass-transfer scaling is confirmed and its coefficient is reported as 0.39–0.46 (0.4 used).'), retrieved: RD2 },
  prov('Sauter-mean to maximum bubble size ratio 0.6', 'slugMassTransfer()', 'none', 'https://arxiv.org/abs/2011.00963', 'unverified', 'Not found in an opened source; it enters the interfacial area of the slug body linearly and is covered by the absorption multiplier.'),
  prov('Asphaltene: δ = 17.347 ρ + 2.904, temperature coefficient 1.07e-3 1/K, colloidal-instability thresholds 0.7 and 0.9', 'asphalteneOnset(), colloidalInstability()', 'patent US 8271248 B2; Scientific Reports 2022 (PMC9643540); ACS Omega (PMC7469123)', 'https://patents.google.com/patent/US8271248B2/en', 'verified', 'Coefficients identical. The density correlation is stated for stock-tank oil and maltenes; applying it to the live oil is this suite\'s simplification. The source prints the temperature law with 273.15 K as reference.'),
  { ...prov('Asphaltene solubility parameter 20 MPa^½ and molar volume 2 m³/kmol', 'inputs asphDelta, asphMV', 'PetroWiki, Thermodynamic models for asphaltene precipitation (archived copy): molar volume 1–4 m³/kmol used by Hirschberg et al.; Panuganti (2013) thesis, Rice University: solubility parameter 19–24 MPa^½', 'https://repository.rice.edu/server/api/core/bitstreams/b56b8d35-a0c9-4c73-9995-8451678abb6a/content', 'verified', 'Both defaults lie inside the quoted ranges; they are fitted when a measured onset pressure is entered.'), retrieved: RD2 },
  prov('Flory–Huggins (Hirschberg) solubility expression', 'asphalteneOnset()', 'Scientific Reports 2022 (PMC9643540)', 'https://pmc.ncbi.nlm.nih.gov/articles/PMC9643540/', 'verified', 'The implicit form ln φ + (1 − φ)(1 − Va/Vs) + (1 − φ)² Va (δa − δs)²/RT = 0 was read; the explicit expression used here is its limit for a small asphaltene fraction.'),
  { ...prov('de Boer screening boundaries: 15 points of undersaturation against in-situ density', 'deBoer(), DE_BOER', 'Ahmed et al. (2023) Iraqi J. Chem. Pet. Eng. 24(1), Fig. 5 (reproduction of the de Boer plot with readable axes)', 'https://ijcpe.uobaghdad.edu.iq/index.php/ijcpe/article/download/933/822', 'corrected', 'The straight lines used before (70 + 1.3 (ρ − 600) bar and 130 bar above) are replaced by points read from the figure, with a reading error of about ± 10 bar and ± 5 kg/m³: at 700 kg/m³ the no-problem limit is 157 bar (200 bar before) and the severe limit 299 bar (330 bar before). A third-hand redraw: systematic error unknown.'), retrieved: RD2 },
  { ...prov('Threshold Shields number 0.30/(1 + 1.2 D*) + 0.055 (1 − exp(−0.020 D*)); viscous bed load: critical Shields number 0.12, mobile-layer thickness θ d/(2 θc), validity of the cubic flux law up to θ ≈ 1.5', 'shieldsCritical(), SHIELDS_MOBILE, sandCriticalVelocity()', 'U.S. Geological Survey Open-File Report 2012-1234, Appendix 3 (Soulsby–Whitehouse); Ouriemi, Aussillous and Guazzelli (2009) J. Fluid Mech. 636, 295, Eq. 3.13–3.16 (HAL)', 'https://hal.science/hal-01432016', 'corrected', 'The constant Shields number 0.05 is replaced by the threshold curve (0.05 in water, 0.24–0.27 for 189 µm sand in 100–340 mPa·s oil).'), retrieved: RD2 },
  { ...prov('Thomas sliding-bed limit for grains inside the viscous sub-layer, V = 9.0 [g ν (s − 1)]^0.37 (D/ν)^0.11', 'sandCriticalVelocity()', 'Poloski et al. (2009) PNNL-17639, Eq. 3.2', 'https://www.pnnl.gov/main/publications/external/technical_reports/PNNL-17639.pdf', 'verified', 'Identical. Blind on the 16 critical velocities of the same report: mean absolute error 23 %.'), retrieved: RD2 },
  { ...prov('Shell diffusivity of a converting droplet, 0.8 × 10⁻¹⁵ m²/s', 'hydrateMarch(), autoclaveConversion(), input shellD', 'Davies (2009) thesis, Colorado School of Mines: 10⁻¹⁷–10⁻¹⁶ m²/s for annealed films (simulator default 10⁻¹⁶), a sensitivity base case of 5 × 10⁻¹⁵ m²/s with 40 µm droplets, 3.4–7.6 × 10⁻¹³ m²/s for fresh films by calorimetry', 'https://repository.mines.edu/server/api/core/bitstreams/b3c44b83-0c2e-449c-8ff7-3fea1c9359fe/content', 'corrected', '5 × 10⁻¹³ → 0.8 × 10⁻¹⁵ m²/s, chosen inside the published range on the autoclave conversion series (calibration set).'), retrieved: RD2 },
  prov('Roughness contribution of deposits (0.1, 0.05, 0.3 of the layer thickness), bed entrainment coefficient, inhibitor film coefficient and diffusivity', 'roughnessUpdate(), run()', 'none', 'https://www.netl.doe.gov/sites/default/files/2018-05/BC15150_Final.pdf', 'unverified', 'Modelling assumptions exposed as inputs or calibration parameters; no relation between deposit thickness and roughness was found (the Tulsa report assumes the bare-pipe roughness).'),
];

// =====================================================================================================
// 17. Verification
// =====================================================================================================
function verify() {
  const out = [], chk = (name, expected, got, tol, note) => out.push({ name, expected, got: Number.isFinite(got) ? got : null, tol, pass: Number.isFinite(got) && Math.abs(got - expected) <= tol, note });
  // --- particles and rheology
  chk('Stokes settling velocity (hand value)', 8.9894e-3, settlingVelocity(1e-4, 2650, 1000, 1e-3, { model: 'stokes' }).v, 2e-6, '100 µm quartz in water: Δρ g d²/(18 μ) = 1650 × 9.80665 × 10⁻⁸ / 0.018');
  { const s = settlingVelocity(1e-3, 2650, 1000, 1e-3), cdBal = (4 * 1650 * G * 1e-3) / (3 * 1000 * s.v * s.v), cdSN = (24 / s.Re) * (1 + 0.15 * s.Re ** 0.687); chk('Schiller–Naumann drag closes the force balance', 1, cdBal / cdSN, 1e-6, `1 mm sand grain: Re = ${s.Re.toFixed(0)}, Cd from weight/drag balance against 24/Re (1 + 0.15 Re^0.687)`); }
  { const d = 5e-5, rp = 2650, rf = 1000, mu = 1e-3, r = particleRelaxation({ d, rhoP: rp, rhoF: rf, mu }); chk('Particle momentum equation: Stokes response time', ((rp + 0.5 * rf) * d * d) / (18 * mu), r.tau, 0.06 * (((rp + 0.5 * rf) * d * d) / (18 * mu)), 'time to 63 % of the terminal velocity equals (ρp + ½ρf) d²/(18 μ) in the Stokes regime'); }
  chk('Krieger–Dougherty reduces to Einstein at small fraction', 1.0025, slurryViscosity(0.001, 'krieger', { phiMax: 4 / 7 }), 1e-5, '1 + 2.5 φ at φ = 0.001');
  chk('Mills relative viscosity (hand value)', 0.8 / (1 - 0.35) ** 2, slurryViscosity(0.2, 'mills', { phiMax: 4 / 7 }), 1e-9, '(1 − φ)/(1 − φ/φmax)² at φ = 0.2, φmax = 4/7 → 1.8935');
  chk('No-cohesion limit of the agglomerate force balance', 1, maxAgglomerateSize({ dp: 4e-5, Fa: 0, mu0: 3e-3, shear: 200, phi: 0.1 }).ratio, 1e-12, 'without cohesive force the stable agglomerate is the primary particle');
  chk('Kozeny–Carman permeability (hand value)', 9.8765e-12, kozenyCarman(0.4, 1e-4), 1e-15, 'ε³ d²/(180 (1 − ε)²) = 0.064 × 10⁻⁸ / 64.8');
  { const g = porousGradient(0.01, 1e-3, 800, 0.4, 1e-4); chk('Darcy–Forchheimer: velocity recovered from the gradient', 0.01, porousVelocity(g.total, 1e-3, 800, 0.4, 1e-4), 1e-9, 'root of μv/k + βρv² for the gradient computed at 0.01 m/s'); }
  chk('Danielson sand critical velocity (hand value)', 0.7958, sandCriticalVelocity({ d: 2e-4, D: 0.2, rhoP: 2650, rhoF: 1000, mu: 1e-3 }).danielson, 2e-3, '0.23 ν^−1/9 d^1/9 (g D (s − 1))^5/9 with ν = 10⁻⁶ m²/s, d = 200 µm, D = 0.2 m');
  // --- population balance
  { const g = pbeGrid(46, 1e-6, 1e-6 * 2 ** 15), n = g.n; let eN = 0, eV = 0; for (let q = 0; q < n * n; q++) { const i = Math.floor(q / n), j = q % n, k = g.pk[q], vs = g.v[i] + g.v[j]; if (vs >= g.v[n - 1]) continue; eN = Math.max(eN, Math.abs(g.pa[q] + g.pb[q] - 1)); eV = Math.max(eV, Math.abs((g.pa[q] * g.v[k] + g.pb[q] * g.v[k + 1]) / vs - 1)); }
    chk('Fixed-pivot allocation: one particle per collision, volume conserved', 0, Math.max(eN, eV), 1e-12, 'largest error of the number and volume assigned to the two neighbouring pivots over all class pairs');
    const N0 = new Float64Array(n); N0[0] = 1e12; const K = 1e-15, t = 4000, r = solvePBE(g, N0, t, { beta: K, frac: 0.02, maxSub: 1e6 }), m = pbeMoments(g, r.N);
    chk('Smoluchowski constant kernel: N(t) = N₀/(1 + K N₀ t/2)', 1 / 3, m.m0 / 1e12, 0.01, 'sectional solution after two characteristic times (K N₀ t/2 = 2)');
    chk('Pure aggregation conserves particle volume', 1, m.vol / (1e12 * g.v[0]), 1e-9, 'Σ N v after aggregation over Σ N v before');
    const m2 = (gg, N) => { let x = 0; for (let i = 0; i < gg.n; i++) x += N[i] * gg.v[i] ** 2; return x / (1e12 * gg.v[0] ** 2); }, g2 = pbeGrid(16, 1e-6, 1e-6 * 2 ** 15), Nc = new Float64Array(16); Nc[0] = 1e12;
    const e1 = Math.abs(m2(g2, solvePBE(g2, Nc, t, { beta: K, frac: 0.02, maxSub: 1e6 }).N) / 5 - 1), e2 = Math.abs(m2(g, r.N) / 5 - 1);
    chk('Size-class refinement reduces the error of the second volume moment', 1, e2 < 0.5 * e1 ? 1 : 0, 0, `Σ N v² = N₀ v₀² (1 + K N₀ t) analytically; relative error ${e1.toFixed(3)} with 16 classes, ${e2.toFixed(3)} with 46 classes over the same size range`);
    const rc = solvePBE(g, N0, t, { beta: K, frac: 0.2, maxSub: 1e6 }), ec = Math.abs(pbeMoments(g, rc.N).m0 / 1e12 - 1 / 3), ef = Math.abs(m.m0 / 1e12 - 1 / 3);
    chk('Sub-step refinement reduces the aggregation error', 0.1, ef / ec, 0.03, 'first-order explicit sub-steps: ten times smaller steps give a ten times smaller error in N(t)');
    const ri = solvePBE(g, N0, t, { beta: K, maxSub: 2 }), mi = pbeMoments(g, ri.N);
    chk('Implicit (stiff) aggregation step conserves volume', 1, mi.vol / (1e12 * g.v[0]), 1e-9, 'modified Patankar–Euler steps on the class volumes');
    const Ng = new Float64Array(n); Ng[10] = 1e10; const gr = 2e-20, rg = solvePBE(g, Ng, 100, { gv: new Float64Array(n).fill(gr) }), mg = pbeMoments(g, rg.N);
    chk('Pure growth conserves particle number', 1e10, mg.m0, 1, 'number after a volume growth step');
    chk('Pure growth adds exactly N·(dv/dt)·t of volume', 1e10 * (g.v[10] + gr * 100), mg.vol, 1e-9 * 1e10 * g.v[10], 'Σ N v after growth');
    const Gl = 2e-9, L0 = g.L[12], Nl = new Float64Array(n); Nl[12] = 1e9; const gvl = Float64Array.from(g.L, (L) => ((PI / 6) * ((L + Gl * 500) ** 3 - L ** 3)) / 500), rl = solvePBE(g, Nl, 500, { gv: gvl }), ml = pbeMoments(g, rl.N);
    chk('Size-independent growth shifts the size by G·t', L0 + Gl * 500, Math.cbrt((6 * ml.vol) / (PI * ml.m0)), 1e-9 * L0, 'volume-mean size after linear growth at G = 2 nm/s for 500 s');
    const Nb = new Float64Array(n), Sb = new Float64Array(n); Nb[20] = 1e8; Sb[20] = 1e-3; const rb = solvePBE(g, Nb, 500, { S: Sb }), mb = pbeMoments(g, rb.N);
    chk('Binary breakage: two daughters per event', 2 - Math.exp(-0.5), mb.m0 / 1e8, 1e-9, 'N/N₀ = 2 − exp(−S t) when only the parent class breaks');
    chk('Breakage conserves particle volume', 1, mb.vol / (1e8 * g.v[20]), 1e-12, 'Σ N v after breakage over Σ N v before');
    const q = qmomNodes([0, 1, 2, 3].map((k) => 2 + 3 ** k)); chk('Wheeler quadrature recovers a two-point distribution', 3, Math.max(...q.L), 1e-9, 'moments of 2 δ(L − 1) + δ(L − 3) give the nodes 1 and 3');
    const qm = solveQMOM([0, 1, 2, 3, 4, 5].map((k) => 1e12 * 1e-5 ** k * Math.exp((k * k * 0.01) / 2)), t, { beta: () => K, steps: 400 });
    chk('Method of moments, constant kernel: m₀(t)', 1 / 3, qm.m[0] / 1e12, 1e-4, 'quadrature method of moments against the Smoluchowski solution');
    chk('Method of moments conserves the third moment', 1, qm.m[3] / (1e12 * 1e-15 * Math.exp(0.045)), 1e-9, 'm₃ after aggregation over m₃ before'); }
  // --- nucleation, kinetics, thermodynamics
  chk('Critical nucleus radius (hand value)', (0.04 * 290) / (914 * 4.6e5 * 10), nucleationRate({ TK: 280, dT: 10, TeqK: 290, sigma: 0.02 }).rc, 2e-15, '2σ Teq/(ρ L ΔT) = 0.04 × 290 / (914 × 4.6×10⁵ × 10) = 2.759 nm');
  chk('Contact-angle factor at 90°', 0.5, nucleationRate({ TK: 280, dT: 10, TeqK: 290, theta: 90 }).f, 1e-12, '(2 + cos θ)(1 − cos θ)²/4');
  { const r = rng(11), ti = inductionTime(7, { TK: 277 }); let s = 0; const N = 6000; for (let k = 0; k < N; k++) s += -Math.log(1 - r.uniform(0, 1)) * ti; chk('Poisson onset: sampled mean equals the induction time', 1, s / N / ti, 0.04, 'mean of 6,000 exponentially distributed onset times over 1/(J V)'); }
  chk('Equilibrium limit: no growth at the hydrate curve', 0, hydrateGrowthRate({ TK: 285, P: 60, Peq: 60 }).j + hydrateDissociationRate({ TK: 285, P: 60, Peq: 60 }), 1e-30, 'growth and dissociation fluxes vanish when the pressure equals the equilibrium pressure');
  { const a = hydrateGrowthRate({ TK: 277.15, P: 80, Peq: 40, zG: 1, kRef: 1e-10, H: 2500 }); chk('Intrinsic kinetics (hand value)', 1e-10 * 40e5, a.j, 1e-9, 'K*(f − f_eq) with ideal gas: 10⁻¹⁰ × 4.0 MPa = 4×10⁻⁴ mol/m²/s'); }
  { const a = hydrateGrowthRate({ TK: 277.15, P: 80, Peq: 40, zG: 1, kRef: 1, H: 2500, kFilm: 1e-5 }); chk('Fast kinetics: growth becomes mass-transfer-limited', 1e-5 * (40e5 / 2500), a.j, 1e-4 * 16, 'flux tends to k_film × Δc when the intrinsic constant is very large'); }
  chk('Kim–Bishnoi dissociation constant (Arrhenius, hand value)', 3.6e4 * Math.exp(-81e3 / (8.314462618 * 280)) * 20e5, hydrateDissociationRate({ TK: 280, P: 30, Peq: 50, zG: 1 }), 1e-9, 'K₀ exp(−E/RT) (f_eq − f) at 280 K with 20 bar of driving force');
  chk('Methane hydrate equilibrium at 10 °C (van der Waals–Platteeuw)', 72.5, vdwpMethane(10).P, 4, 'measured three-phase pressure of methane hydrate at 283.15 K is 71–73 bara');
  chk('Methane hydrate equilibrium at 0 °C', 26.0, vdwpMethane(0).P, 1.5, 'measured quadruple-point pressure is about 25.6–26.3 bara');
  // --- moving boundary
  { const m = plugMeltTime({ R: 0.127, U: 3, Tamb: 4, Td: -1, eps: 0.4 }); chk('Stefan radial melting: front tracking against the analytic time', 1, m.tNumeric / m.tAnalytic, 5e-3, 'ρ(1 − ε)L (R/2U + R²/4k)/ΔT'); }
  chk('Stefan similarity constant at small Stefan number', Math.sqrt(0.005), stefanLambda(0.01), 1e-3, 'λ → √(Ste/2); exact 0.07059 at Ste = 0.01');
  // --- transport operator
  { const n = 60, tau = new Array(n).fill(10), E = new Array(n).fill(0); let t = 0, f = 0, t50 = 0; while (t < 2000) { const fo = advectImplicit(E, tau, 1, 1); t += 1; if (!t50 && fo >= 0.5) { t50 = t - (fo - 0.5) / Math.max(fo - f, 1e-12); } f = fo; } chk('Advection of a tracer front: arrival time', 600, t50, 18, 'half-height breakthrough after 60 cells of 10 s residence time'); }
  // --- wax and scale
  chk('Hayduk–Minhas wax diffusivity (hand value)', 2.287e-10, waxDiffusivity(300, 5e-3), 2e-12, '13.3×10⁻¹² T^1.47 μ^(10.2/V − 0.791)/V^0.71 at 300 K, 5 cP, V = 430 cm³/mol');
  chk('Wax diffusion flux (hand value)', 3.832e-7, waxDeposition({ Tb: 30, Tamb: 10, U: 10, hIn: 1e12, kOil: 0.13, rhoOil: 800, muOil: 5e-3, wat: 40, wTot: 0.05, slope: 0.04 }).jMol, 4e-9, 'ρ D (dC/dT)(q/k): 800 × 2.322×10⁻¹⁰ × 1.3406×10⁻³ × 1538.5 kg/m²/s');
  chk('Davies activity coefficient (hand value)', 0.3733, activityCoefficient(2, 0.1, 25, 'davies'), 2e-3, 'divalent ion at I = 0.1 mol/L, 25 °C: log γ = −0.509 × 4 × (0.2402 − 0.03)');
  { const r = scaleIndices({ ...SEAWATER, HCO3: 110 }, 25, 1.01325, { yCO2: 4.0e-4, model: 'truesdellJones' }); chk('Calcite saturation of surface seawater', 0.7, r.minerals[0].SI, 0.3, `surface seawater (1.8 mmol/L bicarbonate, 400 µatm CO₂) is 4–6 times supersaturated with calcite, SI 0.6–0.8; computed pH ${r.pH.toFixed(2)}`); }
  { const c = mineralSolubility('barite', { Tc: 25 }), w = { Ba: c * 137330 * 0.99705, SO4: c * 96064 * 0.99705 }; chk('Barite: a water at its own solubility is exactly saturated', 0, scaleIndices(w, 25, 1.01325, { yCO2: 1e-9 }).minerals[1].SI, 2e-3, 'round trip mol/kg → mg/L → saturation index through the solubility solver and the index routine'); }
  chk('Colloidal instability index (hand value)', 58 / 42, colloidalInstability({ sat: 55, aro: 28, res: 14, asp: 3 }).cii, 1e-12, '(55 + 3)/(28 + 14)');
  // --- the marching solver on a uniform cold loop
  { const S = labSetup({ n: 16, L: 1600, D: 0.1, T: 4, P: 80, U: 0, tAmb: 4 }), v = { ...DEF, regime: 'oil', nucA: 30 }, g = pbeGrid(12, 5e-6, 5e-3), ph = [{ dur: 1800, dt: 60, frac: 1 }];
    const base = runHydrateMarch(S, hydrateParams(v, S), { phases: ph, grid: g, cheap: true }), L = base.ledger;
    chk('Hydrate mass balance of the march', 1, L.out / L.in, 1e-9, 'suspended + deposited + exported + dissociated over formed');
    let fa = 0; for (let i = 0; i < S.n; i++) fa += (PI / 4) * base.D[i] ** 2 * S.ds + base.mDep[i] / (HYDRATE.rho * (1 - base.por[i])); chk('Flow area + deposit area = clean bore', (PI / 4) * 0.01 * 1600, fa, 1e-9 + (base.ledger.deposited > 0 ? 0 : 1e9), 'pipe volume is conserved as the deposit grows');
    const zero = runHydrateMarch(S, hydrateParams({ ...v, kinK: 0 }, S), { phases: ph, grid: g, cheap: true });
    chk('Zero-kinetics limit: no hydrate', 0, Math.max(...zero.ser.phi), 1e-4, 'only the seed shells (0.01 % of the water) exist when the rate constant is zero');
    const S2 = labSetup({ n: 16, L: 1600, D: 0.1, T: 4, P: 80, U: 5e4, tAmb: 4, wcut: 0.05 }), fast = runHydrateMarch(S2, hydrateParams({ ...v, kinK: 1e4, shellD: 1e7, mtMult: 1000, htMult: 20, adhesion: 0 }, S2), { phases: ph, grid: g, cheap: true });
    chk('Fast-kinetics limit: the limiting reactant (water) is consumed', 1, Math.max(...fast.rec.X), 0.05, 'water conversion with very fast kinetics, mass transfer and heat removal');
    const nd = runHydrateMarch(S, hydrateParams({ ...v, adhesion: 0, filmMult: 0 }, S), { phases: ph, grid: g, cheap: true });
    chk('No-deposition limit', 0, nd.ledger.deposited, 1e-12, 'zero capture efficiency leaves the wall clean');
    const na = runHydrateMarch(S, hydrateParams({ ...v, aggEff: 0, kBreak: 0 }, S), { phases: ph, grid: g, cheap: true }), pP = hydrateParams(v, S); let kP = 0; for (let j = 1; j < g.n; j++) if (Math.abs(Math.log(g.L[j] / pP.dPrim)) < Math.abs(Math.log(g.L[kP] / pP.dPrim))) kP = j;
    chk('No-agglomeration limit: particles stay at the primary size', g.L[kP], Math.max(...na.rec.d43), 1e-9 * g.L[kP], 'd43 with zero collision efficiency equals the seeded size class');
    // steady adiabatic loop: water conservation through the hydration number, and the latent heat carried by the stream
    const Sa = labSetup({ n: 12, L: 1200, D: 0.1, T: 4, P: 80, U: 0, tAmb: 4 }), pa = hydrateParams({ ...v, hydNumber: 6, adhesion: 0, filmMult: 0, aggEff: 0, cohesion: 0 }, Sa), ad = runHydrateMarch(Sa, pa, { phases: [{ dur: 6 * 3600, dt: 120, frac: 1 }], grid: g, cheap: true }), prA = Sa.props(80, 4), K = ad.ser.t.length - 1;
    const expRate = (ad.ser.exp[K] - ad.ser.exp[K - 5]) / (5 * 120), mcpA = prA.mG * prA.cpG + prA.mO * prA.cpO + prA.mW * prA.cpW, freeOut = ad.rec.freeW[Sa.n - 1] * prA.mW;
    chk('Water balance through the hydration number', prA.mW, freeOut + expRate * ((6 * 18.015) / (18 + 6 * 18.015)), 0.01 * prA.mW, 'steady loop: water fed = free water leaving + 0.8572 × hydrate leaving (n = 6, gas 18 g/mol)');
    chk('Energy balance with the heat of formation', expRate * HYDRATE.latent, mcpA * (ad.T[Sa.n - 1] - 4), 0.02 * expRate * HYDRATE.latent, 'adiabatic steady loop: ṁ cp (T_out − T_in) = hydrate formation rate × latent heat');
    // lumped cooldown against the exponential and its time-step order
    const Sc = labSetup({ n: 4, L: 400, D: 0.25, T: 60, P: 80, U: 3, tAmb: 4, C: 1e5 }), pc = hydrateParams({ ...v, kinK: 0, nucA: -99 }, Sc), tEnd = 6 * 3600, exact = 4 + 56 * Math.exp((-3 * PI * 0.25 * tEnd) / 1e5), Tn = (dt) => runHydrateMarch(Sc, pc, { phases: [{ dur: tEnd, dt, frac: 0 }], grid: g, cheap: true }).T[0];
    const T1 = Tn(150), T2 = Tn(300), T4 = Tn(600);
    chk('Shut-in cooldown against the exponential solution', exact, T1, 0.05, 'T = Tamb + (T₀ − Tamb) exp(−U π D t/C) after 6 h');
    chk('Time-step convergence is first order', 1, Math.log2((T4 - T2) / (T2 - T1)), 0.1, 'observed order from three step sizes (150, 300, 600 s)'); }
  // --- statistical thermodynamics: Kihara potential and Gibbs-energy minimisation
  { const cs = kiharaLangmuir(273.15, HYD_REF.cages.small), cl = kiharaLangmuir(273.15, HYD_REF.cages.large), pp = langmuirOccupancy(273.15, 1);
    chk('Kihara cell potential: small-cage Langmuir constant against the Parrish–Prausnitz fit', pp.Cs, cs, 0.15 * pp.Cs, 'independent routes to the methane constant at 273.15 K (potential integral against the empirical A/T·exp(B/T)), 1/atm');
    chk('Kihara cell potential: large-cage Langmuir constant against the Parrish–Prausnitz fit', pp.Cl, cl, 0.25 * pp.Cl, 'same comparison for the 5¹²6² cage');
    chk('Kihara potential: quadrature converged', 1, kiharaLangmuir(280, HYD_REF.cages.small, HYD_REF.kihara, 960) / kiharaLangmuir(280, HYD_REF.cages.small, HYD_REF.kihara, 240), 1e-4, 'Langmuir constant with 960 against 240 integration points');
    const g = hydrateGibbsMin({ Tc: 4, P: 100, nW: 1, nGas: 1, nSolute: 0.02 }), xw = Math.exp(g.dg);
    chk('Gibbs-energy minimisation: equilibrium conversion with a dissolved solute', 1 - (0.02 * xw) / (1 - xw), g.conversion, 1e-6, 'minimum of G(ξ) against the stationarity condition ln x_w = Δg, i.e. ξ = n_w − n_s x_w/(1 − x_w)');
    chk('Gibbs-energy minimisation: gas-limited conversion', 0.05 * g.hydrationNumber, hydrateGibbsMin({ Tc: 4, P: 100, nW: 1, nGas: 0.05 }).conversion, 1e-9, 'without solute G falls linearly, so the minimum sits where the gas is exhausted: ξ = n·n_gas');
    chk('Gibbs-energy minimisation: no hydrate above the equilibrium temperature', 0, hydrateGibbsMin({ Tc: 20, P: 50 }).conversion, 1e-12, 'Δg > 0 at 20 °C and 50 bara, the minimum is at zero conversion'); }
  // --- Pitzer model and mineral solubility against PHREEQC's own published results and limiting laws
  chk('Debye–Hückel slope at 25 °C', 0.5092, waterDH(25).A, 2e-3, 'A = 0.509 (mol/kg)^−½ from the dielectric constant of water (78.4)');
  { const C = { Na: 1, Cl: 1 }, sp = speciate(C, 25, 1.01325, { model: 'pitzer' }), gNa = Math.exp(sp.g[0]), gCl = Math.exp(sp.g[7]);
    chk('Pitzer: mean activity coefficient of 1 mol/kg NaCl at 25 °C', 0.657, Math.sqrt(gNa * gCl), 0.004, 'tabulated γ± of sodium chloride (Robinson and Stokes: 0.657)');
    chk('Pitzer: water activity of 1 mol/kg NaCl at 25 °C', 0.9669, sp.aw, 6e-4, 'osmotic coefficient 0.936: a_w = exp(−2 × 0.936/55.51)');
    chk('Barite solubility in water at 50 °C (PHREEQC, pitzer.dat)', 1.502e-5, mineralSolubility('barite', { Tc: 50 }), 2e-8, 'mol/kg water printed in the USGS test output Barite_NaCl.out');
    chk('Barite solubility in 1 mol/kg NaCl at 50 °C (PHREEQC, pitzer.dat)', 1.294e-4, mineralSolubility('barite', { Tc: 50, bg: { Na: 1, Cl: 1 } }), 2e-7, 'same test output: the salting-in by a factor 8.6 comes from the Pitzer interaction terms');
    chk('Barite solubility in 4 mol/kg NaCl at 50 °C (PHREEQC, pitzer.dat)', 2.002e-4, mineralSolubility('barite', { Tc: 50, bg: { Na: 4, Cl: 4 } }), 3e-7, 'same test output');
    chk('Ion-association model at infinite dilution', 1, mineralSolubility('celestite', { Tc: 25, model: 'truesdellJones' }) / mineralSolubility('celestite', { Tc: 25, model: 'davies' }), 0.03, 'Truesdell–Jones and Davies activity models agree for a dilute SrSO₄ solution (ionic strength 0.002)'); }
  // --- Monte Carlo population balance
  { const N0 = 1e12, K = 1e-15, mc = solvePBEMonteCarlo(new Array(4000).fill(1e-6), 4000, { beta: K, conc: N0, seed: 5 });
    chk('Monte Carlo population balance, constant kernel: N(t)', 1 / 3, mc.m0 / N0, 0.012, 'direct simulation with 4,000 particles against N₀/(1 + K N₀ t/2)');
    chk('Monte Carlo population balance conserves particle volume', 1, mc.vol / ((N0 * PI) / 6 * 1e-18), 1e-9, 'Σ v per unit volume after aggregation over the initial value');
    const env = { shear: 200, mu: 3e-3, nu: 3e-3 / 800, TK: 277, dRho: 120, alpha: 0.05 }, g = pbeGrid(40, 2e-5, 2.4e-3), kP = 5, L0 = g.L[kP], c0 = 0.05 / g.v[kP], tB = 3 / (aggregationKernel(L0, L0, env).total * c0), Ni = new Float64Array(g.n); Ni[kP] = c0;
    const be = new Float64Array(g.n * g.n); for (let i = 0; i < g.n; i++) for (let j = 0; j < g.n; j++) be[i * g.n + j] = aggregationKernel(g.L[i], g.L[j], env).total;
    const sec = pbeMoments(g, solvePBE(g, Ni, tB, { beta: be, maxSub: 4000, frac: 0.05 }).N), m2 = solvePBEMonteCarlo(new Array(1500).fill(L0), tB, { beta: (a, b) => aggregationKernel(a, b, env).total, conc: c0, seed: 5 });
    chk('Monte Carlo against the sectional solver (shear, settling and Brownian kernel)', sec.m0 / c0, m2.m0 / c0, 0.04, 'number remaining after three collision times; the sectional value carries the numerical diffusion of the geometric grid');
    const mb = solvePBEMonteCarlo(new Array(6000).fill(1e-4), 500, { S: () => 1e-3, conc: 1e8, seed: 2 });
    chk('Monte Carlo breakage: N(t) = N₀ exp(S t)', Math.exp(0.5), mb.m0 / 1e8, 0.05, 'every particle breaks in two at the constant frequency S = 10⁻³ 1/s for 500 s'); }
  // --- Maxey–Riley particle equation and Lagrangian tracking
  { const d = 1e-4, rp = 2650, rf = 1000, mu = 1e-3, tau = ((rp + 0.5 * rf) * d * d) / (18 * mu), vt = ((rp - rf) * G * d * d) / (18 * mu);
    chk('Maxey–Riley without history force: exponential approach to the Stokes velocity', 1 - Math.exp(-3), maxeyRiley({ d, rhoP: rp, rhoF: rf, mu, tEnd: 3 * tau, n: 3000, basset: false }).v.at(-1) / vt, 2e-4, 'v/v_t after three response times (ρp + ½ρf) d²/(18 μ)');
    const tE = 400 * tau, b = maxeyRiley({ d, rhoP: rp, rhoF: rf, mu, tEnd: tE, n: 1000, basset: true });
    chk('Maxey–Riley with the Basset history force: algebraic approach to the terminal velocity', d / 2 / Math.sqrt((PI * mu * tE) / rf), 1 - b.v.at(-1) / vt, 2e-3, 'long-time deficit 1 − v/v_t = a/√(π ν t) (a = radius), far slower than the exponential without history');
    const n = 20, f = { ds: 250, D: new Array(n).fill(0.2), vL: new Array(n).fill(1.5), uStar: new Array(n).fill(0.07), rhoF: new Array(n).fill(800), mu: new Array(n).fill(3e-3), theta: new Array(n).fill(0), cold: new Array(n).fill(false) };
    const tr = trackParticles(f, { n: 40, d: 2e-5, rhoP: 800, seed: 4, equip: [{ x: 2500, eff: 1, name: 'filter' }] });
    chk('Lagrangian tracking: a fully retaining equipment boundary traps every parcel', 40, tr.equipment, 0, 'trap rule: efficiency 1 at 2.5 km, no parcel escapes');
    const tn = trackParticles(f, { n: 60, d: 2e-5, rhoP: 800, seed: 4 });
    chk('Lagrangian tracking: mean transit time of neutrally buoyant parcels', (n * 250) / 1.5, tn.meanTransit, 0.06 * ((n * 250) / 1.5), 'tracer parcels stay uniformly mixed (random walk with the drift of the eddy diffusivity), so the mean transit time is L/v');
    chk('Lagrangian tracking: parcels are conserved', 60, tn.wall + tn.bed + tn.equipment + tn.escaped + tn.flying, 0, 'every released parcel is trapped, has escaped or is still in flight');
    const ts = trackParticles({ ...f, uStar: new Array(n).fill(0.008), vL: new Array(n).fill(0.15) }, { n: 30, d: 3e-4, rhoP: 2650, seed: 4 });
    chk('Lagrangian tracking: grains settle where the shear velocity is below the settling velocity', 30, ts.bed, 0, '300 µm sand, u* = 8 mm/s against a settling velocity of about 30 mm/s: all parcels join the bed (reflect/trap rule)'); }
  // --- Eulerian–Eulerian solids transport: advection–diffusion benchmark, conservation, deposition-layer grid convergence
  { const mk = (n, L, u, D, W, ws) => ({ ds: L / n, A: new Array(n).fill(0.05), u: new Array(n).fill(u), W: new Array(n).fill(W), Dax: new Array(n).fill(D), ws: new Array(n).fill(ws), pd: new Array(n).fill(1), E: new Array(n).fill(0) });
    const err = (n) => { const r = solidsTransport(mk(n, 1000, 1, 0.5, 0, 0), { cIn: 1, tEnd: 500, dt: 1000 / n / 2 }); let e = 0; for (let i = 0; i < n; i++) e = Math.max(e, Math.abs(r.c[i] - ogataBanks((i + 0.5) * (1000 / n), 500, 1, 0.5))); return { e, r }; }, e2 = err(200), e4 = err(400);
    chk('Advection–diffusion benchmark (Ogata–Banks front)', 0, e4.e, 0.004, 'largest error of the concentration front after 500 s on 400 cells (u = 1 m/s, D = 0.5 m²/s)');
    chk('Advection–diffusion benchmark: error falls with grid refinement', 1, e4.e < 0.5 * e2.e ? 1 : 0, 0, `largest error ${e2.e.toFixed(4)} on 200 cells and ${e4.e.toFixed(4)} on 400 cells`);
    chk('Eulerian solids transport conserves the solids volume', e4.r.ledger.in, e4.r.ledger.out + e4.r.ledger.suspended + e4.r.ledger.bed + e4.r.ledger.trapped, 1e-9, 'fed = left + suspended + bed + retained');
    const lam = (2e-4 * 0.2) / 0.05, bedAt = (n) => { const r = solidsTransport(mk(n, 2000, 1, 0, 0.2, 2e-4), { cIn: 1e-3, tEnd: 40000, dt: 20, limiter: false }); return interp1(r.bed.map((_, i) => (i + 0.5) * (2000 / n)), r.bed, 1000) / 40000; }, exact = 2e-4 * 0.2 * 1e-3 * Math.exp(-lam * 1000);
    const b1 = bedAt(20), b2 = bedAt(40), b3 = bedAt(80);
    chk('Deposition-layer growth against the analytic profile', 1, b3 / exact, 0.04, 'steady deposition from plug flow: layer growth rate w W c₀ exp(−w W x/(A u)) at mid-length, 80 cells (start-up of 1,000 s included in the 40,000 s average)');
    chk('Deposition-layer grid convergence is first order', 1, Math.log2((b1 - b2) / (b2 - b3)), 0.25, 'observed order of the layer growth rate on 20, 40 and 80 cells (upwind transport)'); }
  // --- reaction–diffusion and mass-transfer-controlled dissociation
  { const r = reactionDiffusion({ L: 1, D: 1e-3, k: 4e-3, n: 80 });
    chk('Reaction–diffusion benchmark: flux into a slab (Thiele modulus 2)', 1e-3 * 2 * Math.tanh(2), r.flux, 2e-6, 'steady solution c = cosh(φ(1 − x))/cosh φ gives the flux D c₀ φ tanh φ / L');
    chk('Reaction–diffusion benchmark: concentration at the closed end', 1 / Math.cosh(2), r.c.at(-1), 2e-4, '1/cosh φ');
    const a = inhibitorDissociation({ mode: 'stagnant', wBulk: 0.9, wEq: 0.3, tEnd: 86400 });
    chk('Mass-transfer-controlled dissociation: diffusion-limited recession', a.sAnalytic, a.s, 0.005 * a.sAnalytic, 'numerical diffusion of inhibitor to the plug face against the similarity solution 2 ρ Δw √(D t/π) × (1 − w_eq)/(w_eq ρ_h (1 − ε) w_water)');
    chk('Mass-transfer-controlled dissociation: no attack below the equilibrium concentration', 0, inhibitorDissociation({ mode: 'film', wBulk: 0.2, wEq: 0.3 }).rate, 0, 'an inhibitor weaker than the equilibrium concentration cannot dissolve the plug'); }
  // --- slug coupling, wax equilibrium, crystallisation, asphaltene, roughness
  { const q = { vsl: 1, vsg: 2, rhoL: 800, rhoG: 70, muL: 3e-3, muG: 1.3e-5, sigma: 0.02, D: 0.25 }, s = slugMassTransfer(q);
    chk('Slug unit cell: Hinze bubble size in the slug body (hand value)', 0.6 * 0.725 * (0.02 / 800) ** 0.6 * s.epsSlug ** -0.4, s.dBubble, 1e-12, 'd32 = 0.6 × 0.725 (σ/ρ)^0.6 ε^−0.4 at the dissipation rate of the slug body');
    chk('Slug unit cell: liquid-volume weighting of the conductance', s.liquidInSlug * s.kLaSlug + (1 - s.liquidInSlug) * s.kLaFilm, s.kLa, 1e-12, 'k_L a of the unit cell is bounded by the film-zone and slug-body values');
    const S1 = labSetup({ n: 8, L: 800, D: 0.1, T: 4, P: 80, U: 5e4, tAmb: 4, wcut: 0.9 }), vv = { ...DEF, regime: 'water', nucA: 30, adhesion: 0, filmMult: 0, kinK: 1e4, mtMult: 0.05 }, gg = pbeGrid(10, 2e-6, 5e-3), ph = [{ dur: 600, dt: 60, frac: 1 }];
    const on = runHydrateMarch(S1, hydrateParams({ ...vv, slugCouple: true }, S1), { phases: ph, grid: gg, cheap: true }), off = runHydrateMarch(S1, hydrateParams({ ...vv, slugCouple: false }, S1), { phases: ph, grid: gg, cheap: true });
    chk('Slug coupling raises absorption-limited hydrate formation', 1, on.ledger.formed > 1.05 * off.ledger.formed ? 1 : 0, 0, `water-dominated loop in slug flow: ${on.ledger.formed.toFixed(1)} kg formed with the unit-cell conductance against ${off.ledger.formed.toFixed(1)} kg with the stratified interface`); }
  { const w = waxEquilibrium(DEFAULT_FLUID), i = w.T.findIndex((t) => t > w.wat - 8);
    chk('Wax solid–liquid equilibrium: no solid above the computed WAT', 0, w.at(w.wat + 1).solid, 1e-12, 'the solubility curve returns the full wax in solution above the cloud point');
    chk('Wax solid–liquid equilibrium: solid fraction rises monotonically on cooling', 1, w.solid.every((s, k) => k === 0 || s <= w.solid[k - 1] + 1e-12) ? 1 : 0, 0, `tabulated curve below the cloud point of ${w.wat === null ? '—' : w.wat.toFixed(1)} °C`);
    chk('Wax solid–liquid equilibrium: precipitation slope is the derivative of the curve', (w.solid[i - 1] - w.solid[i + 1]) / (w.T[i + 1] - w.T[i - 1]), w.at(w.T[i]).dCdT, 0.25 * w.at(w.T[i]).dCdT, 'dC/dT used in the diffusion flux against a central difference of the tabulated solid fraction'); }
  { const t = 1000, r = scaleCrystallisation({ S0: 5, cA: 1, cB: 1, M: 0.2334, rho: 4480, TK: 330, t, sigma: 0.05, A: 1e20, kg: 1e-9 / 16, deplete: false });
    chk('Crystallisation moments: number of crystals at constant supersaturation', r.J0 * t, r.mu[0], 1e-9 * r.J0 * t, 'μ₀ = J t');
    chk('Crystallisation moments: third moment at constant supersaturation', (r.J0 * r.G0 ** 3 * t ** 4) / 4, r.mu[3], 1e-2 * ((r.J0 * r.G0 ** 3 * t ** 4) / 4), 'μ₃ = J G³ t⁴/4 for size-independent growth; the nuclei enter at the critical size, which adds less than 1 %');
    const q = scaleCrystallisation({ S0: 12, cA: 0.3, cB: 0.15, M: 0.2334, rho: 4480, TK: 343, t: 6 * 3600, sigma: 0.09, A: 1e30, kg: 1e-10 });
    chk('Crystallisation conserves the lattice ions', q.x, ((PI / 6) * 4480 * q.mu[3]) / 0.2334, 1e-6 * q.x, 'moles removed from solution = crystal volume × density / molar mass');
    chk('Crystallisation relaxes to saturation', 1, q.S, 0.05, 'a strongly supersaturated barite solution (S = 12) ends at S ≈ 1 after 6 h'); }
  { const st = { P: [300, 200, 100], T: [90, 90, 90], dL: [16, 15.2, 16.5], vL: [0.38, 0.4, 0.45], pBubble: 200, method: 'density', rho0: 0 }, vA = 2, phiA = 0.02, TK = 363.15, dA = 20 * (1 - 1.07e-3 * (TK - 298.15)), hand = Math.exp(vA / 0.4 - 1 - (vA * 1e-3 * ((dA - 15.2) * 1e3) ** 2) / (R * TK));
    chk('Flory–Huggins asphaltene solubility (hand value)', hand, asphalteneOnset({ pre: st, phiA, deltaA: 20, vA }).phiMax[1], 1e-12 * hand + 1e-300, 'φ_max = exp[v_a/v_L − 1 − v_a (δ_a − δ_L)²/RT] at the bubble point of a three-point path');
    const o = asphalteneOnset({ pre: st, phiA: 1e-9, deltaA: 20, vA, refOnset: 250 });
    chk('Asphaltene onset fitted to a measured pressure', 250, o.upperOnset, 0.01, 'the asphaltene solubility parameter is adjusted until the upper onset equals the measurement (solid-phase reference state)');
    const f = makeFluid({ comp: { nC5: 100 } }), e = eosSolubilityParameter(f, f.z, 1.01325, 298.15, 'liquid');
    chk('Equation-of-state solubility parameter of n-pentane', 14.4, e.delta, 1.0, 'Hildebrand value 14.4–14.5 MPa^½ at 25 °C; cohesive energy density from the Peng–Robinson residual internal energy'); }
  { chk('Roughness update: clean wall', 4.5e-5, roughnessUpdate(4.5e-5, 0, 0, 0, 0.25), 0, 'no deposit leaves the clean roughness');
    chk('Roughness update: layers add in proportion and are capped', 0.05 * 0.25, roughnessUpdate(4.5e-5, 0.2, 0.01, 0.01, 0.25), 1e-15, 'a thick deposit cannot make the wall rougher than 5 % of the bore');
    const r1 = roughnessUpdate(4.5e-5, 2e-3, 1e-3, 5e-4, 0.25), D1 = 0.25 - 2 * 3.5e-3, S0 = labSetup({ n: 2, L: 100, D: 0.25, vsl: 2, vsg: 1e-9 }), g0 = S0.grad(0, 80, 4, 0.25, 4.5e-5, 1, 1), g1 = S0.grad(0, 80, 4, D1, r1, 1, 1), Re0 = (g0.rhoM ?? 840) * 0, pr = S0.props();
    const fr = (D, e) => { const v = (2 * 0.25 ** 2) / D ** 2; return (frictionFactor((pr.rhoL * v * D) / pr.muL, e / D) * v * v) / D; }; void Re0;
    chk('Roughness-update consistency with the hydraulics', fr(D1, r1) / fr(0.25, 4.5e-5), g1.fric / g0.fric, 0.03 * (fr(D1, r1) / fr(0.25, 4.5e-5)), 'single-phase liquid: the kernel friction gradient with the updated bore and roughness equals the Colebrook ratio f(ε/D) v²/D computed by hand'); }
  // --- boundaries inside the time march: heating, injection, equipment, initial population
  { const Sc = labSetup({ n: 4, L: 400, D: 0.25, T: 4, P: 80, U: 3, tAmb: 4, C: 1e5 }), g = pbeGrid(10, 5e-6, 5e-3), pc = hydrateParams({ ...DEF, regime: 'oil', kinK: 0, nucA: -99 }, Sc);
    const hot = runHydrateMarch(Sc, pc, { phases: [{ dur: 60 * 3600, dt: 600, frac: 0 }], grid: g, cheap: true, heat: new Array(4).fill(50) });
    chk('Heating boundary: temperature of a heated shut-in section', 4 + (50 / (3 * PI * 0.25)) * (1 - Math.exp((-3 * PI * 0.25 * 60 * 3600) / 1e5)), hot.T[0], 0.02, 'T = T_amb + q′/(U π D) (1 − exp(−U π D t/C)) for 50 W/m after 60 h');
    chk('Heating boundary: energy supplied', 50 * 400 * 60 * 3600, hot.ledger.heatIn, 1, 'q′ × length × time');
    const Sl = labSetup({ n: 12, L: 1200, D: 0.1, T: 4, P: 80, U: 0, tAmb: 4 }), pl = hydrateParams({ ...DEF, regime: 'oil', nucA: 30, kRemove: 0 }, Sl), ph = [{ dur: 3600, dt: 60, frac: 1 }];
    const tr = runHydrateMarch(Sl, pl, { phases: ph, grid: g, cheap: true, traps: [{ i: 8, eff: 0.5 }] }), Lg = tr.ledger;
    chk('Equipment boundary: hydrate mass balance with a trap', 1, Lg.out / Lg.in, 1e-9, 'solids retained at the equipment are part of the deposited mass');
    chk('Equipment boundary retains solids', 1, Lg.trapped > 0 && tr.mDep[8] >= Lg.trapped * 0.5 ? 1 : 0, 0, `${Lg.trapped.toFixed(2)} kg retained in the trap cell`);
    const ini = runHydrateMarch(Sl, hydrateParams({ ...DEF, regime: 'oil', kinK: 0, nucA: -99, adhesion: 0, filmMult: 0, aggEff: 0, kBreak: 0 }, Sl), { phases: [{ dur: 60, dt: 60, frac: 0 }], grid: g, cheap: true, init: { phi: 0.05, d: 2e-4, sig: 1 } });
    chk('Initial hydrate population: suspended mass', 0.05 * (PI / 4) * 0.01 * 1200 * Sl.prof(1).holdup[0] * HYDRATE.rho, ini.ledger.susp0, 1e-6, 'φ × liquid volume × hydrate density');
    const Sj = labSetup({ n: 6 }); Sj.aq = { ...Sj.aq, inh: INHIBITORS.MEG, inhId: 'MEG', inhWtCell: [0, 0, 0, 30, 30, 30] };
    const pj = hydrateParams({ ...DEF, regime: 'oil', kinK: 0, nucA: -99 }, Sj), rj = runHydrateMarch(Sj, pj, { phases: [{ dur: 600, dt: 300, frac: 1 }], grid: g, cheap: true, inhWt: Sj.aq.inhWtCell });
    chk('Mid-line injection: hydrate temperature falls by the inhibitor depression downstream of the point', hydrateDepression({ S: 0, inhWt: 30, inh: INHIBITORS.MEG }), rj.rec.teq[1] - rj.rec.teq[4], 1e-9, '30 wt % MEG from the fourth cell on'); }
  // --- population balance on the resolved flow field
  { const Fl = pipeFlowField({ D: 0.25, U: 1, rho: 800, mu: 0.2, nr: 12 }), Ft = pipeFlowField({ D: 0.25, U: 0.4, rho: 1000, mu: 1e-3, nr: 12 }), nuT = 1e-6;
    chk('Flow field, laminar limit: Darcy friction factor 64/Re', 64 / 1000, Fl.f, 1e-9, 'Re = 1000, eddy viscosity switched off below Re = 2300');
    chk('Flow field, laminar limit: centre-line velocity is twice the bulk velocity', 2, Fl.fine.u[Fl.fine.u.length - 1] / Fl.U, 2e-3, 'Hagen–Poiseuille parabola recovered by the wall-to-axis integration');
    chk('Flow field: the cell velocities carry the bulk flow', 0.4 * PI * 0.125 ** 2, Ft.u.reduce((q, u, j) => q + u * Ft.dA[j], 0), 1e-12, 'Σ u_j A_j = U A');
    chk('Flow field, turbulent: friction factor against the smooth-pipe law', frictionFactor(1e5, 0), Ft.f, 0.02 * frictionFactor(1e5, 0), 'Re = 10⁵: Colebrook–Prandtl 0.0180; the Cess eddy viscosity with the pipe constants κ = 0.42, A⁺ = 27 is within a few per cent');
    { const yq = (100 * nuT) / Ft.uStar, j = Ft.fine.y.findIndex((y) => y >= yq), up = (Ft.fine.u[j - 1] + ((Ft.fine.u[j] - Ft.fine.u[j - 1]) * (yq - Ft.fine.y[j - 1])) / (Ft.fine.y[j] - Ft.fine.y[j - 1])) / Ft.uStar;
      chk('Flow field, turbulent: velocity at y⁺ = 100 against the logarithmic law', Math.log(100) / 0.41 + 5.2, up, 0.5, 'u⁺ = ln(y⁺)/0.41 + 5.2 = 16.4');
      const y5 = (5 * nuT) / Ft.uStar, j5 = Ft.fine.y.findIndex((y) => y >= y5); chk('Flow field, turbulent: viscous sub-layer u⁺ = y⁺', 5, Ft.fine.u[j5] / Ft.uStar * (y5 / Ft.fine.y[j5]), 0.3, 'at y⁺ = 5 the damped eddy viscosity is still small'); }
    const gF = pbeGrid(12, 1e-5, 2e-2), N0 = new Array(12).fill(0); N0[3] = 0.05 / gF.v[3];
    const Fo = pipeFlowField({ D: 0.25, U: 1.2, rho: 800, mu: 3e-3, nr: 10 }), kinF = { alpha: 0.05, kBreak: 0.04, dA: 5e-4, dPrim: gF.L[3], rhoP: 914 };
    const r1 = solvePBEField({ field: Fo, grid: gF, Nin: N0, L: 50, nx: 16, kin: kinF, wall: { stick: 0.5, adhForce: 5e-3 } });
    chk('Flow-field population balance: third moment (particle volume) flux in = out + deposited', 0, (r1.ledger.in - r1.ledger.out - r1.ledger.deposited) / r1.ledger.in, 1e-10, 'aggregation, breakage, turbulent diffusion and wall deposition together');
    const r2 = solvePBEField({ field: Fo, grid: gF, Nin: N0, L: 50, nx: 16, kin: { alpha: 0 }, wall: { stick: 1 } });
    chk('Flow-field population balance: zeroth moment (number) flux conserved without aggregation or breakage', 0, (r2.ledger.numIn - r2.ledger.numOut - r2.ledger.numDeposited) / r2.ledger.numIn, 1e-10, 'number in = number out + number deposited');
    const r3n = solvePBEField({ field: Fo, grid: gF, Nin: N0, L: 50, nx: 16, kin: { alpha: 0 }, wall: { stick: 0 } });
    chk('Flow-field population balance: no kinetics and no sticking leaves the population unchanged', gF.L[3], r3n.d43[16][5], 1e-12 * gF.L[3], 'uniform inlet stays uniform');
    { // well-mixed limit: with very fast radial mixing the section behaves as one stirred volume with the area-mean kernel coefficients
      const rr = solvePBEField({ field: Fo, grid: gF, Nin: N0, L: 5, nx: 80, kin: kinF, wall: { stick: 0 }, mixing: 1e6, maxSub: 4000 }), A = PI * 0.125 ** 2, mA = (a) => a.reduce((q, x, j) => q + x * Fo.dA[j], 0) / A, cS = mA(rr.kernels.map((k) => k.cS)), kb = mA(rr.kernels.map((k) => k.kb));
      const z = solvePBE(gF, N0, 5 / 1.2, { beta: [cS, (0.05 * (PI / 4) * Math.abs(914 - 800) * G) / (18 * 3e-3), (0.05 * 2 * KB * 277) / (3 * 3e-3)], S: gF.L.map((L) => (L > 1.01 * gF.L[3] ? kb * L ** 3 : 0)), maxSub: 100000, frac: 0.02 }), mz = pbeMoments(gF, z.N);
      chk('Flow-field population balance: well-mixed limit reproduces the 0-D solution (number)', 1, rr.mixed.m0 / mz.m0, 0.03, 'radial diffusivity × 10⁶; 0-D sectional solver with area-mean shear over the residence time L/U');
      chk('Flow-field population balance: well-mixed limit reproduces the 0-D solution (d43)', 1, rr.mixed.d43 / mz.d43, 0.06, 'same case'); }
    { // grid convergence of the wall deposition (radial cells and axial steps refined together)
      const dep = [6, 12, 24].map((nr) => solvePBEField({ field: pipeFlowField({ D: 0.25, U: 1.2, rho: 800, mu: 3e-3, nr }), grid: gF, Nin: N0, L: 50, nx: 2 * nr, kin: kinF, wall: { stick: 0.5, adhForce: 5e-3 }, maxSub: 400 }).ledger.deposited), q = (dep[1] - dep[0]) / (dep[2] - dep[1]);
      chk('Flow-field population balance: grid convergence of the deposited volume (ratio of successive differences)', 1.6, q, 0.6, `6, 12 and 24 radial cells: ${dep.map((d) => d.toExponential(3)).join(', ')} m³/s; a ratio above one means convergence (first-order scheme: about 2)`);
      chk('Flow-field population balance: change on the last refinement', 0, (dep[2] - dep[1]) / dep[2], 0.06, 'deposited volume changes by less than 6 % from 12 to 24 radial cells'); }
    { // analytical deposition case: plug flow, constant diffusivity, partly absorbing wall (Biot number 1): first eigenvalue of λ J1(λ) = Bi J0(λ)
      const nr = 40, Rr = 0.1, Uu = 1, Gm = 1e-3, Vd = Gm / Rr, rf = Array.from({ length: nr + 1 }, (_, j) => (Rr * j) / nr), rc = rf.slice(1).map((r, j) => 0.5 * (r + rf[j]));
      const Fs = { nr, R: Rr, U: Uu, rho: 1000, mu: 1e-3, uStar: 0.05, tauW: 2.5, rf, rc, dA: rf.slice(1).map((r, j) => PI * (r * r - rf[j] ** 2)), u: new Array(nr).fill(Uu), gamma: new Array(nr).fill(0), epsT: new Array(nr).fill(0), nut: new Array(nr).fill(Gm), nutF: new Array(nr - 1).fill(Gm) };
      const ra = solvePBEField({ field: Fs, grid: gF, Nin: N0, L: 40, nx: 400, kin: { alpha: 0, brownian: false }, wall: { stick: 1, vd: () => Vd }, ScT: 1 });
      chk('Analytical deposition case: decay rate of the particle flux (first Bessel eigenvalue, Biot number 1)', (1.2558 ** 2 * Gm) / (Uu * Rr * Rr), Math.log(ra.flux.m0[200] / ra.flux.m0[400]) / 20, 0.02 * (1.2558 ** 2 * Gm) / (Uu * Rr * Rr), 'N ∝ J0(λ1 r/R) exp(−λ1² Γ x / (U R²)), λ1 = 1.2558 for V_d R/Γ = 1'); }
  }
  // --- measured-rate nucleation, ramp statistics, exact shell growth
  { const T = 285 - 3, hand = (0.58e-3 / 11e-4) * Math.exp((22.2 * 3) / T) * Math.exp(-2660 / (T * 9)) + (0.10 / 11e-4) * Math.exp((22.2 * 3) / T) * Math.exp(-3.1e4 / (T * 9));
    chk('Measured-rate nucleation: rate per unit interface at 3 K subcooling (hand value)', hand, hydrateNucleationRate({ dT: 3, TeqK: 285 }), 1e-9 * hand, 'two site populations, A exp(Δs ΔT/kT) exp(−B′/(T ΔT²)) with the published constants');
    chk('Measured-rate nucleation: no nucleation at or above the hydrate temperature', 0, hydrateNucleationRate({ dT: 0, TeqK: 285 }) + hydrateNucleationRate({ dT: -2, TeqK: 285 }), 1e-300, 'J = 0 for ΔT ≤ 0');
    chk('Induction time × nucleation rate × interface = 1', 1, hydrateInductionTime(4, { TeqK: 285, area: 2 }) * hydrateNucleationRate({ dT: 4, TeqK: 285 }) * 2, 1e-12, 'mean of the exponential waiting time');
    const flat = onsetRamp({ rate: 1e-3, area: 0.01, A1: 5, B1: 0, A2: 0, dse: 0, dTmax: 0.6, n: 6000 });
    chk('Ramp integration, constant-rate limit: median onset = ln 2 · β/(J α)', (Math.LN2 * 1e-3) / (5 * 0.01), flat.median, 2e-5, 'with a subcooling-independent rate the hazard is linear in ΔT');
    chk('Ramp integration, constant-rate limit: mean onset = β/(J α)', 1e-3 / (5 * 0.01), flat.mean, 2e-5, 'exponential distribution');
    const r1 = onsetRamp({ rate: 1 / 60, TeqK: 287.9 }), r3 = onsetRamp({ rate: 3 / 60, TeqK: 287.9 });
    chk('Ramp at 1 K/min in the calibration cell: mean onset subcooling against the measured value', 6.60, r1.mean, 0.6, 'published statistics of 408 methane-hydrate formation events: 6.60 ± 1.3 K');
    chk('Ramp at 3 K/min in the calibration cell: mean onset subcooling against the measured value', 7.28, r3.mean, 0.6, 'published statistics of 456 events: 7.28 ± 0.96 K');
    chk('Ramp at 1 K/min: spread of the onset against the measured standard deviation', 1.3, r1.sd, 0.4, 'the stochastic model also reproduces the width of the distribution');
    chk('Faster cooling delays the onset to a larger subcooling', 1, +(r3.median > r1.median && r1.median > onsetRamp({ rate: 1 / 3600, TeqK: 287.9 }).median), 0, '1 K/h < 1 K/min < 3 K/min');
    const kd = 2e-5 / (2 * 1e-15), g = (X) => 1 - 3 * (1 - X) ** (2 / 3) + 2 * (1 - X);
    chk('Shrinking core, shell-diffusion limit: X after a = (Rd/2D) g(0.5) (analytical)', 0.5, shrinkingCoreAdvance(0, kd * g(0.5), 0, 2e-5, 1e-15), 1e-9, 't/τ = 1 − 3(1 − X)^⅔ + 2(1 − X)');
    chk('Shrinking core, surface-kinetics limit: X = a/Ra', 0.3, shrinkingCoreAdvance(0, 0.3 * 7e9, 7e9, 2e-5, Infinity), 1e-9, 'constant resistance: linear conversion');
    chk('Shrinking core: two half steps equal one step (exact advance, no step-size error)', shrinkingCoreAdvance(0, 4e9, 1e9, 2e-5, 1e-15), shrinkingCoreAdvance(shrinkingCoreAdvance(0, 2e9, 1e9, 2e-5, 1e-15), 2e9, 1e9, 2e-5, 1e-15), 1e-9, 'the explicit step used before over-predicted the conversion of a thin shell by orders of magnitude');
    chk('Film growth through its own resistance: parabolic law δ = √(2 D a)', Math.sqrt(2 * 1e-15 * 3e-5), filmAdvance(0, 3e-5, 0, 1e-15), 1e-15, 'no other resistance');
    chk('Autoclave conversion: none below the methane hydrate pressure', 0, autoclaveConversion({ tH: 24, Tc: 10, P: 50 }), 1e-12, 'methane hydrate needs about 72 bara at 10 °C'); }
  // --- scale kinetics and inhibition, sand regimes, screening curves
  { const lt = 1.5232 - 10.8784 / 2 - 895.6683 / 298.15 + 5476.992 / (2 * 298.15);
    chk('Barite induction time at SI 2, 25 °C (hand value)', 10 ** lt, scaleInductionTime('barite', { SI: 2, TK: 298.15 }), 1e-6 * 10 ** lt, 'log10 t0 = 1.5232 − 10.8784/SI − 895.6683/T + 5476.992/(SI T)');
    const q = scaleInhibition('barite', { SI: 1.3, TK: 353.15, pH: 6, R: 1, dose: 0, tProtect: 86400 }), q2 = scaleInhibition('barite', { SI: 1.3, TK: 353.15, pH: 6, R: 1, dose: q.mic, tProtect: 86400 });
    chk('Minimum inhibitor concentration: the induction time at that dose equals the protection time', 86400, q2.tInh, 1e-6 * 86400, 'MIC = log10(t_protect/t0)/b');
    chk('Threshold inhibition raises the critical saturation index', 1, +(scaleCriticalSI('barite', { TK: 298.15, tProtect: 7200, dose: 5 }) > scaleCriticalSI('barite', { TK: 298.15, tProtect: 7200, dose: 0 }) + 0.5), 0, '5 mg/L of phosphonate at 25 °C: 1.40 → 2.48');
    chk('Critical barite saturation index for a 2 h residence at 25 °C', 1.40, scaleCriticalSI('barite', { TK: 298.15, tProtect: 7200, dose: 0 }), 0.02, 'the index at which the published induction-time relation gives 2 h');
    chk('No inhibition credit for minerals without an induction model', 0, scaleInhibition('calcite', { SI: 0.8, TK: 330, dose: 50 }).eff, 0, 'carbonates and calcium sulphates are judged on the index alone');
    const sv = sandCriticalVelocity({ d: 189e-6, D: 0.0776, rhoP: 2650, rhoF: 880, mu: 0.2, C: 5e-5 });
    chk('Sand in laminar oil: mobile bed-load limit (hand value)', (1.5 * (2650 - 880) * G * 189e-6 * 0.0776) / (8 * 0.2), sv.governing, 1e-9, 'wall Shields number 8 μ V/(D Δρ g d) = 1.5; regime "' + sv.regime + '"');
    const st = sandCriticalVelocity({ d: 15e-6, D: 0.078, rhoP: 7950, rhoF: 1000, mu: 1e-3, C: 0.09 });
    chk('Fine particles inside the viscous sub-layer: Thomas limit (hand value)', 9.0 * (G * 1e-6 * 6.95) ** 0.37 * (0.078 / 1e-6) ** 0.11, st.governing, 1e-9, 'regime "' + st.regime + '", d⁺ = ' + st.dPlus.toFixed(2));
    chk('Coarse sand in water stays on the turbulent-suspension correlations', 1, +(sandCriticalVelocity({ d: 4e-4, D: 0.2, rhoP: 2650, rhoF: 1000, mu: 1e-3 }).regime === 'turbulent suspension'), 0, 'grain larger than the viscous sub-layer');
    chk('Threshold Shields number, viscous limit', 0.30, shieldsCritical(1e-7, 2650, 900, 1), 2e-3, 'D* → 0');
    chk('Threshold Shields number, coarse-grain limit', 0.055, shieldsCritical(0.05, 2650, 1000, 1e-3), 1e-3, 'D* → ∞');
    chk('de Boer boundary at 700 kg/m³ (tabulated point)', 156.5, deBoer(700, 0).lower, 1e-9, '2 270 psi read from the plot');
    chk('de Boer classification above the severe boundary', 1, +(deBoer(650, 300).cls === 'severe problems' && deBoer(800, 100).cls === 'no problems'), 0, '650 kg/m³ with 300 bar of undersaturation; 800 kg/m³ with 100 bar');
    const wm = waxDeposition({ Tb: 30, Tamb: 10, U: 10, hIn: 1e12, kOil: 0.13, rhoOil: 800, muOil: 5e-3, wat: 35, wTot: 0.05, slope: 0.04, delta: 0, D: 0.25, vL: 1, gammaW: 0, wmodel: 'matzain' });
    chk('Matzain deposit oil fraction from the Reynolds-number closure (hand value)', 1 - (800 * 1 * 0.25 / 5e-3) ** 0.15 / 8, wm.Coil, 1e-12, '1 − Re^0.15/8 at Re = 40 000');
    chk('Matzain enhancement Π1 = C1/(1 − C_oil)', 15 / (1 - wm.Coil), wm.pi1, 1e-9, 'C1 = 15'); }
  // --- preservation of the shut-in line
  { const Sp = labSetup({ n: 4, L: 1000, D: 0.2, U: 5, tAmb: 4 }), mxp = { inhWt: Sp.aq.inhWtCell, heat: null, heatStart: 0 }, vP = { ...DEF, marginC: 3, pShut: 0 };
    const a = applyPreservation(Sp, mxp, vP, { id: 'inhibit', active: true, chem: 'MEG', dose: 30, lead: 4 });
    chk('Preservation by inhibitor: every cell carries the dose', 30, Math.min(...a.S.aq.inhWtCell), 1e-12, 'and the unprotected set-up keeps its own concentration');
    chk('Preservation by inhibitor: the unprotected set-up is not modified', 0, Math.max(...Sp.aq.inhWtCell), 1e-12, 'copies, not references');
    const hh = applyPreservation(Sp, mxp, vP, { id: 'heat', active: true, chem: 'MEG', dose: 0, lead: 4 });
    chk('Preservation by heating: holding power (hand value)', 5 * PI * 0.2 * (hydrateT0(80, 0.7) + 3 - 4) * 1000, hh.mx.holdPower, 1e-6, 'U π D (T_hold − T_ambient) L with T_hold = hydrate temperature + margin'); }
  return out;
}

// =====================================================================================================
// 18. Suite declaration
// =====================================================================================================
const CAL_TARGETS = [{ key: 'phi', label: 'Hydrate fraction', unit: 'vol %' }, { key: 'dpRatio', label: 'Pressure-drop ratio', unit: '–' }, { key: 'teqC', label: 'Hydrate equilibrium temperature at 80 bara', unit: '°C' }, { key: 'tIndH', label: 'Induction time', unit: 'h' }, { key: 'kLa', label: 'Gas-uptake (mass-transfer) coefficient', unit: '1/s' }, { key: 'dTex', label: 'Temperature excess during formation (heat transfer)', unit: '°C' }, { key: 'd43Um', label: 'Agglomerate size d43', unit: 'µm' }, { key: 'muRel', label: 'Relative (effective) viscosity', unit: '–' }, { key: 'depMm', label: 'Hydrate deposit thickness', unit: 'mm' }, { key: 'permM2', label: 'Deposit permeability', unit: 'm²' }, { key: 'tDisS', label: 'Time to dissociate a particle', unit: 's' },
  { key: 'waxMm', label: 'Wax thickness', unit: 'mm' }, { key: 'waxSolidPct', label: 'Precipitated wax at the wall temperature', unit: 'wt %' }, { key: 'scaleUm', label: 'Scale crystal size', unit: 'µm' }, { key: 'scaleMmY', label: 'Scale deposition rate', unit: 'mm/y' }, { key: 'aopBar', label: 'Asphaltene onset pressure', unit: 'bara' }, { key: 'vSetMm', label: 'Particle settling velocity', unit: 'mm/s' }, { key: 'vCrit', label: 'Sand deposition (critical) velocity', unit: 'm/s' }, { key: 'vResus', label: 'Resuspension velocity', unit: 'm/s' }, { key: 'entrRate', label: 'Sand entrainment rate', unit: 'kg/m²/s' }, { key: 'erosMmY', label: 'Erosion rate', unit: 'mm/y' }];
const okNum = (x, lo, hi) => (typeof x === 'number' && Number.isFinite(x) && x >= lo && x <= hi ? x : null);
export default {
  id: 'solids', num: 4, title: 'Hydrate & Multiphase Solids Flow Assurance', short: 'Hydrate · Solids', icon: '❄️',
  tagline: 'Where, when and how fast hydrates form, agglomerate, deposit and plug — with wax, scale, asphaltene and sand on the same line.',
  validationData,
  description: 'The pressure–temperature picture of the line is carried through an operating scenario (steady production, turndown, shut-in cooldown, cold restart). A shutdown or cold restart is solved as the preserved line (inhibitor placed before cooldown, depressurisation or heating) with the unprotected event marched beside it. Hydrate onset follows the nucleation rate measured at gas–water interfaces (two populations of sites, stochastic waiting time; classical theory is kept for comparison); growth is intrinsic kinetics in series with mass transfer, diffusion through the hydrate shell (advanced exactly over each step) and heat removal; particles agglomerate and break in a sectional population balance and travel with the liquid; the wall captures them, the deposit ages, is sheared off and narrows the bore, which feeds back on velocity, shear and pressure. Wax (solid–liquid equilibrium of the paraffin distribution, diffusion-driven deposition and ageing), mineral scale (Pitzer speciation, saturation indices, nucleation induction time with threshold inhibitor and sulphate removal, population-balance crystallisation), asphaltene onset (regular-solution model on the depletion path, with a dispersant option) and sand transport (critical velocity by flow regime, Eulerian–Eulerian bed model, optional Lagrangian tracking) are evaluated for steady production and all layers are combined into one deposit profile for the network and flow suites.',
  guide: [
    'Choose the operating scenario. Steady production of the reference case is far outside the hydrate region, so use turndown, shutdown cooldown or cold restart to test the hydrate strategy.',
    'Link the line, fluid and flow picture from the other suites (or leave the reference case) and enter the present inhibitor dose.',
    'For a shutdown or restart choose the preservation strategy (linked from the operations suite when it has been run); the design case is the preserved line and the unprotected event is reported beside it.',
    'Enter the design measures that are in place: scale inhibitor and dose, sulphate removal on the injection water, asphaltene dispersant, wax inhibitor.',
    'Review the kinetic, cohesion and wall-capture parameters; fit them to flow-loop or rocking-cell data on the Calibration tab when you have measurements.',
    'Enter the produced-water analysis, wax data, SARA analysis and sand rate for the steady-production threats.',
    'Run, then read the plugging probability with its interval, the inhibitor dose, the cooldown limit, the pigging interval and the governing threat by zone.',
    'Use the mesh and step study to quantify the numerical uncertainty before relying on a time to plug.',
    'Optional: switch on the Lagrangian particle tracker, enter chemical-injection points, a heated section or equipment that retains solids, and paste laboratory or field measurements into the Measurements table to compare them with the model.',
  ],
  implemented: [
    // equations
    'van der waals–platteeuw theory', 'langmuir occupancy equations', 'kihara potential', 'chemical-potential equality', 'fugacity equilibrium', 'gibbs-energy minimization', 'hydrate phase-stability equations', 'classical nucleation theory', 'homogeneous nucleation', 'heterogeneous nucleation', 'stochastic/poisson nucleation models', 'induction-time distributions', 'kim–bishnoi kinetic model', 'englezos–bishnoi-type kinetic models', 'intrinsic kinetic growth equations', 'mass-transfer-limited growth', 'heat-transfer-limited growth', 'combined heat/mass-transfer growth', 'shrinking-core-type formulations', 'population balance equation', 'number-density transport equation', 'method of moments', 'sectional methods', 'monte carlo population models', 'smoluchowski coagulation equation', 'collision-frequency kernels', 'aggregation kernels', 'breakage kernels', 'shear-induced collision', 'turbulent collision', 'differential-settling collision', 'cohesive-force/agglomeration models', 'eulerian solids transport', 'eulerian–eulerian multiphase equations', 'eulerian–lagrangian particle dynamics', 'maxey–riley-type particle equations', 'drag-force equations', 'settling equations', 'stokes settling', 'hindered settling', 'mass-transfer deposition equations', 'wall-capture models', 'adhesion probability models', 'shear-removal models', 'deposition–erosion competition', 'moving-boundary equations', 'stefan-type moving-interface formulation', 'effective-diameter reduction', 'effective-area reduction', 'permeability models', 'darcy flow through porous deposits', 'darcy–forchheimer equation', 'kozeny–carman permeability model', 'evolving-roughness models', 'arrhenius kinetics', 'kim–bishnoi-type dissociation', 'heat-transfer-controlled dissociation', 'mass-transfer-controlled dissociation', 'stefan moving-boundary models', 'solid–liquid equilibrium', 'wax appearance thermodynamics', 'molecular-diffusion deposition', 'shear dispersion', 'brownian diffusion', 'aging models', 'ionic activity equations', 'saturation-index equations', 'solubility-product equilibrium', 'precipitation kinetics', 'population-balance crystallization', 'solubility models', 'flory–huggins theory', 'eos-based precipitation', 'solid-phase equilibrium models', 'particle momentum equation', 'stokes/schiller–naumann drag', 'erosion/deposition models',
    // hybrid models
    'thermodynamic hydrate stability + kinetic formation', 'eulerian fluid + lagrangian hydrate particles', 'hydrate kinetics + slug-flow model', 'population balance + cfd', 'deposition + evolving hydraulic diameter', 'hydrate + wax + scale + sand competing-solids framework',
    // initial and boundary conditions
    'initial hydrate-free or hydrate-containing state of the system', 'initial dissolved water and free-water distribution', 'hydrate saturation or concentration', 'hydrate particle number and size distribution', 'initial agglomerate population', 'initial deposited-hydrate thickness', 'initial effective pipe diameter and roughness', 'initial inhibitor concentration throughout the fluid', 'initial wax', 'sand and other solid concentrations where these phenomena are enabled',
    'hydrate-forming gas', 'liquid and water entering the domain', 'inlet particle or solid concentrations where applicable', 'inhibitor concentration and injection rate at chemical-injection points', 'wall temperature and heat-transfer conditions controlling hydrate formation', 'wall adhesion', 'detachment and resuspension behaviour', 'particle behaviour at inlets', 'outlets and equipment', 'reflection', 'trapping', 'entrainment or escape as appropriate', 'local thermodynamic environment for hydrate formation and dissociation', 'while shutdown', 'restart', 'depressurization and heating boundaries should permit hydrate formation', 'growth and melting to evolve with operations',
    // inputs and outputs
    'coupled pressure/temperature/velocity/holdup fields from module 3', 'thermodynamic hydrate stability properties from module 1', 'water availability', 'hydrate nucleation, induction, growth and dissociation parameters', 'heat/mass-transfer parameters', 'particle density, size distribution and population parameters', 'agglomeration/breakage kernels', 'slurry rheology', 'wall adhesion/capture, deposition, erosion, detachment and resuspension parameters', 'inhibitor concentration/effectiveness', 'deposit porosity/permeability', 'wax/asphaltene/scale thermodynamic and kinetic properties where enabled', 'sand/solids properties', 'hydrate stability margin and exposure', 'onset/induction location and time', 'nucleation and growth rates', 'hydrate mass/volume fraction and water/gas conversion', 'particle number/concentration/size distribution', 'slurry properties and transport', 'wall-capture/deposition/removal rates', 'deposit thickness and distribution', 'dissociation/melting', 'pressure-drop increase', 'effective diameter/area reduction', 'evolving roughness/permeability', 'blockage/plugging indicator, location and time with stated uncertainty', 'inhibitor effectiveness', 'wax/scale/asphaltene/sand precipitation, transport, deposition, resuspension and erosion outputs where enabled',
    // calibration: every quantity has a fitted parameter and a target column in the calibration model
    'hydrate equilibrium parameters', 'nucleation-rate parameters', 'induction-time parameters', 'growth-rate constants', 'mass-transfer coefficients', 'heat-transfer coefficients', 'particle-growth parameters', 'particle-size distribution', 'agglomeration kernels', 'breakage kernels', 'cohesion/adhesion parameters', 'wall-capture efficiency', 'deposition coefficient', 'erosion/removal coefficient', 'critical wall shear', 'deposit porosity', 'deposit permeability', 'effective viscosity', 'dissociation kinetics', 'wax precipitation parameters', 'wax deposition coefficient', 'wax diffusion coefficients', 'wax-removal parameters', 'scale nucleation/growth constants', 'scale deposition parameters', 'asphaltene precipitation parameters', 'particle settling velocity', 'sand entrainment', 'sand deposition', 'resuspension threshold', 'erosion coefficients',
    // verification
    'species mass conservation', 'solid-phase mass conservation', 'energy conservation including latent/reaction heat', 'population-balance conservation', 'particle-number conservation where applicable', 'analytical nucleation/growth cases', 'analytical deposition cases', 'zero-kinetics limit', 'infinite/very-fast kinetics limiting behaviour', 'no-deposition limit', 'no-agglomeration limit', 'equilibrium limit', 'dissociation-limit tests', 'particle transport benchmarks', 'advection-diffusion benchmarks', 'reaction-diffusion benchmarks', 'population-balance benchmark solutions', 'mesh convergence', 'time-step convergence', 'deposition-layer grid convergence', 'moving-boundary verification', 'flow-area conservation as deposits evolve', 'roughness-update consistency', 'thermodynamic/kinetic coupling consistency',
    // validation: a tick means that the comparison is supported in the app (a sourced reference data set and/or the Measurements table with prediction and error metrics)
    'high-pressure hydrate loop experiments', 'rocking-cell experiments', 'autoclave experiments', 'flow-loop experiments', 'hydrate onset measurements', 'induction-time measurements', 'hydrate growth measurements', 'particle-size measurements', 'agglomeration measurements', 'slurry-rheology measurements', 'deposition experiments', 'wall-adhesion measurements', 'plugging experiments', 'dissociation experiments', 'inhibitor experiments', 'restart experiments', 'deepwater p–t condition experiments', 'field hydrate incidents/data where reliable', 'wax flow-loop data', 'scale experiments', 'asphaltene experiments', 'sand-transport/erosion experiments',
  ],
  referenceOnly: [],
  equationsNote: 'Scope and limits. The hydrate curve is the kernel screening curve (gas-gravity correlation with salt and inhibitor depression) so that all suites agree; the van der Waals–Platteeuw model (Parrish–Prausnitz or Kihara-potential Langmuir constants) and the Gibbs-energy minimisation are solved for methane structure I and supply the hydration number, a cross-check and the equilibrium conversion in a concentrating brine. The line is one-dimensional. In the hydrate march particles travel with the liquid (Eulerian inventories, implicit upwind transport); slip, settling, wall rules and equipment traps are resolved by the Lagrangian tracker (optional mode: Maxey–Riley equation in a random-walk turbulence field built from the one-dimensional flow picture) and, for sand, by the Eulerian–Eulerian transport with bed exchange. The population balance is also solved on a resolved flow field: the size classes are carried through the axisymmetric velocity, shear and eddy-viscosity field of the pipe section with the most hydrate (Reynolds-averaged momentum balance with an algebraic eddy viscosity; the liquid is taken to fill the pipe; wall deposition by a deposition-velocity wall function); three-dimensional and gas–liquid resolved cases are handed to an external solver. Scale chemistry uses the PHREEQC pitzer.dat model (valid to halite saturation, checked to 200 °C and 600 bar on the reference sets) or the phreeqc.dat ion-association model; iron chemistry is limited to Fe²⁺ without sulphide. Wax thermodynamics is an ideal-solution solid–liquid equilibrium on the single-carbon-number split. The asphaltene onset is a regular-solution model whose equation-of-state variant is anchored to the density correlation at reservoir conditions; de Boer and the colloidal instability index remain cross-checks. Nucleation uses rate constants measured for methane in a stirred cell and referred to the gas–liquid interface; rocking-cell onsets of methane and of a natural gas are reproduced blind, but the rate is apparatus-dependent and extrapolating it to the interface of a pipeline cell is an assumption (multiplier input). The shell diffusivity is calibrated on one autoclave series. Scale induction times and inhibition are fitted relations for barite and celestite (± 0.5 log units); carbonates and calcium sulphates take no inhibitor credit. Sand in laminar oil uses a mobile-bed Shields limit that has only one data set behind it. The wax-loop comparison is against solvent-washed thickness and is within a factor of two to three. Cohesion, wall capture, agglomeration and crystallisation constants remain effective values to be fitted on the Calibration tab. A validation tick means that the comparison is supported in the app (sourced data set and/or the Measurements table); slurry viscosity against hydrate fraction, particle size, plugging time and location, deposit thickness, restart and field incidents have no openly licensed point data behind them (the search is recorded in the notes of the data sets); the status of every literature constant is listed in PROVENANCE. The plugging probability reflects only the parameter uncertainty that is sampled.',
  inputs: INPUTS,
  presets: [
    { name: 'Cold restart after a 48 h shutdown, uninhibited', values: { scenario: 'restart', shutHours: 48, restartPct: 50, inhibitor: 'none', preserve: 'none' } },
    { name: 'MEG-inhibited steady production, under-insulated line', values: { scenario: 'steady', uValue: 16, inhibitor: 'MEG', inhWt: 40, simHours: 24 } },
    { name: 'Gas-dominated wet-gas line', values: { scenario: 'steady', fluidSystem: 'wetgas', gasRate: 3, gasWater: 20, uValue: 12, simHours: 48, sandRate: 5, waxContent: 0 } },
    { name: 'Waxy crude on a long tie-back', values: { scenario: 'steady', lengthScale: 2, wat: 52, waxContent: 9, depositDays: 60, waxLimitMm: 3 } },
    { name: 'High-barium water with seawater breakthrough (scale)', values: { scenario: 'steady', swFrac: 40, water: [{ Na: 28000, K: 600, Ca: 2400, Mg: 300, Ba: 260, Sr: 420, Fe: 12, Cl: 49000, SO4: 5, HCO3: 480 }], co2Pct: 4 } },
    { name: 'Sand-producing late-life well at turndown', values: { scenario: 'turndown', turndownPct: 25, fluidSystem: 'highwc', highWc: 70, sandRate: 900, sandUm: 300, simHours: 36 } },
    { name: 'Shutdown cooldown for 24 h', values: { scenario: 'shutdown', shutHours: 24 } },
    { name: 'Cold restart with mid-line MEG injection, heated riser base and a choke that retains solids', values: { scenario: 'restart', shutHours: 36, injections: [{ x: 0, rate: 20 }, { x: 9, rate: 15 }], injChem: 'MEG', heatWm: 80, heatFromKm: 15, heatToKm: 19, heatStartH: 12, equip: [{ x: 19, eff: 20, name: 'riser-base choke' }], waterSettle: 'lowpoints' } },
    { name: 'Particle tracking of sand at turndown (Lagrangian parcels and bed model)', values: { scenario: 'turndown', turndownPct: 30, sandRate: 600, sandUm: 350, trackOn: true, trackWhat: 'sand', nParcels: 40, simHours: 12, sandHours: 48, nMC: 0 } },
    { name: 'Shut-in with hydrate already present, depressurised after 12 h', values: { scenario: 'shutdown', shutHours: 48, initHydPct: 3, initAggUm: 300, depressAtH: 12, depressP: 15 } },
  ],
  pull: ({ fluid, outputs } = {}) => {
    const pvt = outputs?.pvt, net = outputs?.net, ops = outputs?.ops, tm = okNum(net?.thermalMass, 1e-3, 1e9), inhOk = fluid?.inhibitor && fluid.inhibitor !== 'none' && INHIBITORS[fluid.inhibitor];
    return [
      okNum(pvt?.wat, -20, 90) !== null && { key: 'wat', value: pvt.wat, from: 'Fluid suite: wax appearance temperature' },
      okNum(pvt?.waxContent, 0, 40) !== null && { key: 'waxContent', value: pvt.waxContent, from: 'Fluid suite: wax content' },
      okNum(net?.id, 0.025, 1.5) !== null && { key: 'idMm', value: net.id * 1000, from: 'Network suite: inner diameter' },
      okNum(net?.roughness, 5e-7, 3e-3) !== null && { key: 'roughUm', value: net.roughness * 1e6, from: 'Network suite: wall roughness' },
      okNum(net?.uValue, 0.2, 200) !== null && { key: 'uValue', value: net.uValue, from: 'Network suite: overall heat-transfer coefficient' },
      okNum(net?.tSeabed, -5, 40) !== null && { key: 'tSeabed', value: net.tSeabed, from: 'Network suite: seabed temperature' },
      tm !== null && okNum(tm > 5000 ? tm / 1000 : tm, 0, 5000) !== null && { key: 'thermalMass', value: tm > 5000 ? tm / 1000 : tm, from: 'Network suite: thermal mass of pipe and contents' },
      inhOk && { key: 'inhibitor', value: fluid.inhibitor, from: 'Case fluid: inhibitor' },
      inhOk && okNum(fluid?.inhWt, 0, 90) !== null && { key: 'inhWt', value: fluid.inhWt, from: 'Case fluid: inhibitor concentration' },
      inhOk && okNum(ops?.inhibitorDose, 0, 90) !== null && { key: 'inhWt', value: ops.inhibitorDose, from: 'Operations suite: inhibitor dose' },
      ['inhibit', 'depressurise', 'heat', 'none'].includes({ inhibit: 'inhibit', inhibitor: 'inhibit', blowdown: 'depressurise', depressurise: 'depressurise', heat: 'heat', hotoil: 'heat', none: 'none' }[String(ops?.preservation ?? ops?.preserve ?? '')]) && { key: 'preserve', value: { inhibit: 'inhibit', inhibitor: 'inhibit', blowdown: 'depressurise', depressurise: 'depressurise', heat: 'heat', hotoil: 'heat', none: 'none' }[String(ops.preservation ?? ops.preserve)], from: 'Operations suite: preservation strategy' },
      okNum(ops?.inhibitorDose, 1, 80) !== null && { key: 'preserveDose', value: ops.inhibitorDose, from: 'Operations suite: inhibitor dose for the shut-in line' },
      okNum(ops?.maxShutdown, 0.5, 2000) !== null && { key: 'shutHours', value: ops.maxShutdown, from: 'Operations suite: longest planned shutdown' },
      okNum(fluid?.comp?.CO2, 0.001, 40) !== null && { key: 'co2Pct', value: Math.min(fluid.comp.CO2 * 1.5, 60), from: 'Case fluid: CO₂ (enriched in the gas phase)' },
    ].filter(Boolean);
  },
  site: (site) => [okNum(site?.data?.seabedTemp, -5, 40) !== null && (site.data.depth === undefined || site.data.depth > 0) && { key: 'tSeabed', value: site.data.seabedTemp, from: `Seabed temperature at site${isNum(site.data.depth) ? ` (${site.data.depth} m water depth)` : ''}` }].filter(Boolean),
  run,
  mesh: [
    { name: 'Axial cells', keys: ['nAxial'], min: 8, note: 'First-order upwind transport of particles and heat along the line.', metrics: [
      { label: 'Peak hydrate fraction', unit: 'vol fraction', get: (r) => r.outputs.hydrateFraction }, { label: 'Peak hydrate inventory', unit: 'kg', get: (r) => r.outputs.hydrateMass }, { label: 'Peak blockage', unit: 'area fraction', get: (r) => r.outputs.blockage }, { label: 'Wax mass', unit: 'kg', get: (r) => r.outputs.waxMass }] },
    { name: 'Time step of the hydrate march', keys: ['dtMin'], refine: 'divide', note: 'Implicit first-order time integration; the shut-in steps scale with the same input.', metrics: [
      { label: 'Peak hydrate inventory', unit: 'kg', get: (r) => r.outputs.hydrateMass }, { label: 'Peak blockage', unit: 'area fraction', get: (r) => r.outputs.blockage }, { label: 'Peak slurry viscosity factor', unit: '×', get: (r) => r.outputs.slurryViscosityFactor }, { label: 'Exposure inside the hydrate region', unit: 'h', get: (r) => r.outputs.exposureHours }] },
    { name: 'Monte Carlo particles of the stochastic population balance', keys: ['mcParticles'], min: 100, note: 'Statistical error falls with the square root of the number of simulated particles.', metrics: [{ label: 'Number remaining after three collision times (Monte Carlo)', unit: 'N/N₀', get: (r) => r.outputs.pbeMonteCarlo }] },
    { name: 'Time steps of the Eulerian solids transport', keys: ['eeSteps'], min: 10, note: 'Flux-limited advection with implicit dispersion and bed exchange.', metrics: [{ label: 'Sand in the bed', unit: 'kg', get: (r) => r.outputs.sandBedMass }, { label: 'Largest bed height', unit: 'm', get: (r) => r.outputs.sandBedHeight }] },
    { name: 'Flow-field population balance: radial cells and axial steps', keys: ['cfdNr', 'cfdNx'], min: 4, note: 'Conservative finite volumes in the radius, implicit upwind marching along the pipe; first order in both.', metrics: [
      { label: 'Flow-weighted d43 at the outlet', unit: 'm', get: (r) => r.outputs.pbeFieldD43 }, { label: 'Wall deposition flux at the outlet', unit: 'm³/m²/s', get: (r) => r.outputs.pbeFieldDeposition }, { label: 'Friction factor of the resolved profile', unit: '–', get: (r) => r.outputs.pbeFieldFriction }] },
    { name: 'Particle-size classes', keys: ['nClasses'], min: 6, note: 'Geometric size grid of the sectional population balance.', metrics: [
      { label: 'Agglomerate d43 at the peak', unit: 'm', get: (r) => r.outputs.particleSize }, { label: 'Peak slurry viscosity factor', unit: '×', get: (r) => r.outputs.slurryViscosityFactor }, { label: 'Peak blockage', unit: 'area fraction', get: (r) => r.outputs.blockage }] },
  ],
  calibration: {
    note: 'Every quantity of the calibration matrix has a parameter that moves it and a measured column. A row is one test: fill only the columns that were measured (a hydrate loop test gives hydrate fraction, pressure-drop ratio, agglomerate size, relative viscosity, deposit thickness, gas-uptake coefficient and temperature excess after a time at constant subcooling and velocity; a cell test gives induction time and equilibrium temperature; a dissociation test the time to dissociate a particle at a superheat; a cold-wall test wax thickness and precipitated wax; a tube test scale crystal size and deposition rate at a saturation index; a sand test settling, critical and resuspension velocities, pick-up rate and erosion). Untick the parameters your data cannot separate — the fit reports weak identifiability. The hydrate model is a 2-inch jacketed loop at 80 bara with 20 % water cut; scale tests are barite at 70 °C; the asphaltene onset is that of the reference fluid at reservoir temperature; sand tests are in water. The sample rows are synthetic (generated with other parameter values plus noise).',
    params: [{ key: 'hydEqOffset', label: 'Hydrate equilibrium offset (°C)', lo: -5, hi: 5 }, { key: 'nucA', label: 'Nucleation rate multiplier, log₁₀', lo: -3, hi: 3 }, { key: 'nucB', label: 'Nucleation barrier constant B′ (K³)', lo: 1000, hi: 10000 }, { key: 'kinK', label: 'Intrinsic rate constant (10⁻¹⁰ mol/m²/Pa/s)', lo: 0.01, hi: 100 }, { key: 'shellD', label: 'Shell diffusivity (10⁻¹⁵ m²/s)', lo: 0.01, hi: 1000 }, { key: 'mtMult', label: 'Gas-absorption (mass-transfer) multiplier', lo: 0.01, hi: 100 }, { key: 'htMult', label: 'Heat-removal multiplier', lo: 0.05, hi: 20 }, { key: 'inhEff', label: 'Inhibitor effectiveness (%)', lo: 20, hi: 150 },
      { key: 'primaryUm', label: 'Primary particle size (µm)', lo: 5, hi: 500 }, { key: 'aggEff', label: 'Collision efficiency (agglomeration kernel)', lo: 0.001, hi: 1 }, { key: 'kBreak', label: 'Breakage coefficient', lo: 0.001, hi: 5 }, { key: 'cohesion', label: 'Cohesive force per particle size (mN/m)', lo: 0.01, hi: 100 }, { key: 'phiMax', label: 'Maximum packing fraction (slurry rheology)', lo: 0.45, hi: 0.74 },
      { key: 'adhesion', label: 'Wall-capture efficiency (deposition coefficient)', lo: 0.0005, hi: 1 }, { key: 'adhForce', label: 'Wall adhesion force per particle size (mN/m)', lo: 0.01, hi: 500 }, { key: 'kRemove', label: 'Removal rate constant (1/h)', lo: 0, hi: 100 }, { key: 'tauCrit', label: 'Critical wall shear (Pa)', lo: 0.5, hi: 500 }, { key: 'porosity0', label: 'Fresh deposit porosity', lo: 0.2, hi: 0.9 }, { key: 'plugGrainUm', label: 'Grain size for the deposit permeability (µm)', lo: 1, hi: 2000 }, { key: 'disK0', label: 'Dissociation constant K₀ (10⁴ mol/m²/Pa/s)', lo: 0.01, hi: 1000 }, { key: 'disE', label: 'Dissociation activation energy (kJ/mol)', lo: 40, hi: 120 },
      { key: 'waxMult', label: 'Wax deposition multiplier (deposition and diffusion coefficient)', lo: 0.05, hi: 20 }, { key: 'waxSlope', label: 'Wax solubility slope (1/K)', lo: 0.01, hi: 0.15 }, { key: 'waxStripC', label: 'Wax shear-stripping coefficient (removal)', lo: 0, hi: 5 },
      { key: 'scaleKg', label: 'Scale crystal growth constant (10⁻¹⁰ m/s)', lo: 0.001, hi: 1e4 }, { key: 'scaleSigma', label: 'Scale crystal interfacial energy (mJ/m²)', lo: 20, hi: 150 }, { key: 'scaleK', label: 'Scale surface-reaction constant (10⁻⁸ mol/m²/s)', lo: 0.001, hi: 1e4 }, { key: 'scaleDepEff', label: 'Scale crystal sticking efficiency', lo: 0, hi: 1 }, { key: 'asphDelta', label: 'Asphaltene solubility parameter (MPa^½)', lo: 17, hi: 24 },
      { key: 'sandShape', label: 'Sand drag multiplier (settling velocity)', lo: 0.5, hi: 5 }, { key: 'sandCoef', label: 'Sand critical-velocity multiplier (deposition)', lo: 0.2, hi: 5 }, { key: 'sandShields', label: 'Critical Shields number (resuspension)', lo: 0.005, hi: 0.5 }, { key: 'sandEntrain', label: 'Bed entrainment coefficient (10⁻⁴ m/s)', lo: 0.01, hi: 1e4 }, { key: 'erosionSm', label: 'Erosion geometry constant', lo: 0.5, hi: 100 }],
    columns: [{ key: 'calT', label: 'Test time', unit: 'h' }, { key: 'calDT', label: 'Subcooling without inhibitor', unit: '°C' }, { key: 'calInh', label: 'MEG in the water', unit: 'wt %' }, { key: 'calV', label: 'Velocity', unit: 'm/s' }, { key: 'calDTw', label: 'Oil − coolant (wax test)', unit: '°C' }, { key: 'calSI', label: 'Saturation index (scale test)', unit: '–' }, { key: 'calSuper', label: 'Temperature above equilibrium (dissociation test)', unit: '°C' }, ...CAL_TARGETS],
    targets: CAL_TARGETS,
    model: loopModel, sample: CAL_SAMPLE, validationSample: CAL_VALID,
  },
  verify,
};
