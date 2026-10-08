// Bridge to external open-source solvers. Formulations that cannot run in a browser at useful resolution —
// three-dimensional LES/DES/DNS, interface-capturing CFD, shell and solid finite elements, two-way fluid–structure
// interaction, particle CFD, multi-parameter reference equations of state, dynamic process simulation — are handed
// off: this module writes a complete, ready-to-run case for an established open solver from the data the case already
// holds, says how to run it, and reads the results back so that they join the rest of the study.
// Generators are pure functions of plain data and return { files: [{ path, text }], readme, summary, ... };
// importers are size-capped parsers of untrusted text. No DOM, no network, no storage.
import { clamp, interp1, isNum, linspace, mean } from './num.js';
import { G, slugUnit, stratifiedLevel } from './pipe.js';
import { COMP_IDS, DEFAULT_FLUID, R, fluidModel, makeFluid, flashPT } from './thermo.js';
import { caseLine, flowPicture } from './caseflow.js';
import { BASE } from '../data/basecase.js';

export const BRIDGE_LIMITS = Object.freeze({ chars: 40e6, rows: 400000, series: 64, files: 400, nodes: 2e6, zipBytes: 200e6 });
export const OPENFOAM_FLAVOUR = 'OpenFOAM v2412 (ESI, openfoam.com); dictionaries are valid for v2212 and later';
const KEL = 273.15;

// ---- small helpers ------------------------------------------------------------------------------------
const num = (v, d) => (isNum(+v) && v !== null && v !== '' ? +v : d);
const pos = (v, d) => (isNum(+v) && +v > 0 ? +v : d);
/** Number formatting for solver input files: fixed number of significant digits, never NaN or Infinity. */
export const ff = (x, sig = 8) => (isNum(x) ? String(+x.toPrecision(sig)) : '0');
const vec = (v, sig = 8) => `(${v.map((c) => ff(Math.abs(c) < 1e-14 ? 0 : c, sig)).join(' ')})`;
const add = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const mul = (a, s) => [a[0] * s, a[1] * s, a[2] * s];
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const hyp = (a) => Math.hypot(a[0], a[1], a[2]);
const safeName = (s) => String(s || 'case').normalize('NFKD').replace(/[^\w.-]+/g, '_').replace(/^_+|_+$/g, '').slice(0, 60) || 'case';
const text = (lines) => lines.flat(Infinity).filter((l) => l !== null && l !== undefined && l !== false).join('\n') + '\n';

// ---- catalogue hand-off map ------------------------------------------------------------------------------
/** Normalisation used by the workspace when it matches catalogue item names. */
export const normItem = (s) => String(s).toLowerCase().normalize('NFKD').replace(/[‐-―]/g, '-').replace(/[’'`]/g, '').replace(/s\b/g, '').replace(/[^a-z0-9]+/g, ' ').trim();

// [exact catalogue item, solver, generator, generator options, what the hand-off resolves]
const VOF = 'Three-dimensional two-phase pipe section (volume-of-fluid interface) ';
const ROWS = {
  pvt: [
    ['GERG-type multiparameter EOS', 'CoolProp (GERG-2008 Helmholtz mixture model)', 'coolpropScript', {}, 'Reference-quality gas density, compressibility, heat capacity and Joule–Thomson coefficient from the GERG-2008 multi-fluid Helmholtz model on the grid of the property table; the result is re-imported to replace or cross-check the cubic-EOS table.'],
    ['Helmholtz free-energy formulation', 'CoolProp (multi-fluid Helmholtz energy model)', 'coolpropScript', {}, 'Properties derived from an explicit Helmholtz energy a(T, ρ, x) with departure functions, evaluated by CoolProp on the grid of the property table.'],
    ['PC-SAFT', 'teqp (PC-SAFT)', 'teqpScript', { model: 'PCSAFT' }, 'Perturbed-chain SAFT densities and derivative properties of both phases at the equilibrium compositions of the case, on the grid of the property table.'],
    ['SAFT', 'teqp (PC-SAFT, the SAFT-family equation for chain molecules)', 'teqpScript', { model: 'PCSAFT' }, 'Statistical associating fluid theory in its perturbed-chain form, evaluated at the phase compositions of the case.'],
  ],
  net: [],
  flow: [
    ['SST k–ω', 'OpenFOAM interFoam (RAS kOmegaSST)', 'openfoamPipeCase', { turbulence: 'kOmegaSST' }, VOF + 'with the k–ω SST closure.'],
    ['RNG k–ε', 'OpenFOAM interFoam (RAS RNGkEpsilon)', 'openfoamPipeCase', { turbulence: 'RNGkEpsilon' }, VOF + 'with the RNG k–ε closure.'],
    ['Realizable k–ε', 'OpenFOAM interFoam (RAS realizableKE)', 'openfoamPipeCase', { turbulence: 'realizableKE' }, VOF + 'with the realizable k–ε closure.'],
    ['Reynolds-stress model', 'OpenFOAM interFoam (RAS LRR)', 'openfoamPipeCase', { turbulence: 'LRR' }, 'Transport of all six Reynolds stresses (Launder–Reece–Rodi): secondary flows and anisotropy near the gas–liquid interface that eddy-viscosity models cannot give.'],
    ['Spalart–Allmaras', 'OpenFOAM interFoam (RAS SpalartAllmaras)', 'openfoamPipeCase', { turbulence: 'SpalartAllmaras' }, VOF + 'with the one-equation Spalart–Allmaras closure.'],
    ['LES', 'OpenFOAM interFoam (LES, WALE sub-grid model)', 'openfoamPipeCase', { turbulence: 'WALE' }, 'Large-eddy simulation of the pipe section: resolved turbulent structures, interface waves and slug initiation.'],
    ['DES', 'OpenFOAM interFoam (kOmegaSSTDES)', 'openfoamPipeCase', { turbulence: 'kOmegaSSTDES' }, 'Detached-eddy simulation: RANS in the wall layer, LES in the core.'],
    ['IDDES', 'OpenFOAM interFoam (kOmegaSSTIDDES)', 'openfoamPipeCase', { turbulence: 'kOmegaSSTIDDES' }, 'Improved delayed detached-eddy simulation with a wall-modelled LES branch.'],
    ['DNS where computationally feasible', 'OpenFOAM interFoam (no turbulence model, resolved)', 'openfoamPipeCase', { turbulence: 'DNS' }, 'Direct numerical simulation set-up with the cell count that full resolution needs — feasible only for short sections at low Reynolds number; the read-me gives the estimate.'],
    ['RANS + VOF', 'OpenFOAM interFoam (RAS kOmegaSST + volume of fluid)', 'openfoamPipeCase', { turbulence: 'kOmegaSST', interface: 'vof' }, 'Reynolds-averaged turbulence with a resolved gas–liquid interface in three dimensions.'],
    ['LES + VOF', 'OpenFOAM interFoam (LES WALE + volume of fluid)', 'openfoamPipeCase', { turbulence: 'WALE', interface: 'vof' }, 'Large-eddy simulation with the volume-of-fluid interface: slug front, wave growth and gas entrainment resolved in three dimensions.'],
    ['Coupled Level-Set/VOF', 'OpenFOAM interIsoFoam (geometric VOF with a reconstructed distance function, plicRDF)', 'openfoamPipeCase', { interface: 'plicRDF', turbulence: 'kOmegaSST' }, 'Mass-conserving geometric VOF advection whose interface normal and curvature come from a signed-distance function reconstructed around the interface — the coupling of a level-set description with volume-of-fluid transport.'],
    ['1-D transient flow + 3-D CFD', 'OpenFOAM interFoam driven by the 1-D transient results', 'openfoamCoupledCase', {}, 'The three-dimensional section is driven by the time histories of the one-dimensional transient model (phase flow rates at its inlet) and returns holdup, pressure drop and forces to compare with the 1-D closures.'],
    ['Two-fluid model', 'OpenFOAM reactingTwoPhaseEulerFoam (three-dimensional two-fluid model)', 'openfoamEulerCase', { phases: 2 }, 'Three-dimensional Euler–Euler two-fluid model with interfacial drag, lift, virtual mass, wall lubrication and turbulent dispersion.'],
    ['Two-fluid + mechanistic closures', 'OpenFOAM reactingTwoPhaseEulerFoam (two-fluid model with mechanistic interfacial closures)', 'openfoamEulerCase', { phases: 2 }, 'Two-fluid transport in three dimensions closed by mechanistic models of each interfacial force.'],
    ['Multi-fluid model', 'OpenFOAM multiphaseEulerFoam (gas, oil and water as separate fluids)', 'openfoamEulerCase', { phases: 3 }, 'Three-dimensional multi-fluid model with one momentum equation per phase (gas, oil, water).'],
    ['Lift force', 'OpenFOAM reactingTwoPhaseEulerFoam (Tomiyama lift)', 'openfoamEulerCase', { phases: 2 }, 'Lateral lift on bubbles in the resolved shear field of the pipe section.'],
    ['Virtual/added mass', 'OpenFOAM reactingTwoPhaseEulerFoam (virtual-mass force)', 'openfoamEulerCase', { phases: 2 }, 'Added-mass force of accelerating bubbles in the three-dimensional two-fluid model.'],
    ['Wall lubrication', 'OpenFOAM reactingTwoPhaseEulerFoam (Antal wall lubrication)', 'openfoamEulerCase', { phases: 2 }, 'Wall-lubrication force that keeps bubbles off the wall, resolved across the pipe section.'],
    ['Turbulent dispersion', 'OpenFOAM reactingTwoPhaseEulerFoam (Burns turbulent dispersion)', 'openfoamEulerCase', { phases: 2 }, 'Dispersion of the gas phase by the resolved turbulence field.'],
    ['Bubble-induced turbulence', 'OpenFOAM reactingTwoPhaseEulerFoam (Lahey k–ε)', 'openfoamEulerCase', { phases: 2 }, 'Turbulence produced by the relative motion of bubbles, as a source in the liquid k–ε equations.'],
    ['Coalescence/breakup', 'OpenFOAM reactingTwoPhaseEulerFoam (population balance of bubble size classes)', 'openfoamPopulationCase', { dispersed: 'gas' }, 'Bubble coalescence and breakup between size classes, resolved on the three-dimensional flow field.'],
  ],
  solids: [
    ['Population balance + CFD', 'OpenFOAM reactingTwoPhaseEulerFoam (size-group population balance)', 'openfoamPopulationCase', {}, 'Size-resolved hydrate particles (agglomeration and breakage between size classes) transported by a three-dimensional Euler–Euler flow.'],
    ['Eulerian–Lagrangian particle dynamics', 'OpenFOAM DPMFoam (Lagrangian parcels in a resolved carrier flow)', 'openfoamErosionCase', { particles: 'sand' }, 'Individual particle parcels tracked through the three-dimensional flow with drag, gravity, turbulent dispersion and wall rebound.'],
    ['Eulerian fluid + Lagrangian hydrate particles', 'OpenFOAM DPMFoam (Lagrangian hydrate parcels)', 'openfoamErosionCase', { particles: 'hydrate' }, 'Hydrate particles as Lagrangian parcels in the resolved carrier flow: where they travel and where they hit the wall.'],
    ['Eulerian–Eulerian multiphase equations', 'OpenFOAM reactingTwoPhaseEulerFoam (solids as a continuum, kinetic theory)', 'openfoamEulerCase', { phases: 2, dispersed: 'solids' }, 'Solids as an interpenetrating continuum with kinetic-theory closures in three dimensions.'],
  ],
  ops: [],
  integ: [
    ['Shell equations', 'CalculiX (S8R shell elements)', 'calculixPipeCase', { element: 'shell' }, 'Shell finite-element model of the span, bend or jumper: local bending at supports and bends, ovalisation and the stress field that beam theory averages out.'],
    ['CFD + FEA', 'OpenFOAM + CalculiX coupled through preCICE', 'fsiCase', { coupling: 'explicit' }, 'Wall loads from the flow solution applied to a finite-element model of the pipe wall.'],
    ['Navier–Stokes + structural dynamics', 'OpenFOAM + CalculiX coupled through preCICE', 'fsiCase', { coupling: 'implicit' }, 'The fluid equations and the structural dynamics solved together, exchanging force and displacement every time step.'],
    ['Two-way FSI', 'OpenFOAM + CalculiX coupled through preCICE (implicit, quasi-Newton)', 'fsiCase', { coupling: 'implicit' }, 'Strongly coupled two-way fluid–structure interaction of a jumper or span under slug loading.'],
    ['Erosion + particle CFD', 'OpenFOAM DPMFoam with a particle-erosion cloud function', 'openfoamErosionCase', { particles: 'sand' }, 'Sand parcels tracked through a bend; impact angle, velocity and the eroded-volume field on the wall.'],
  ],
  econ: [],
};
/** Suite id → hand-off entries. `match` equals the normalised whole catalogue item name; `item` is the exact catalogue string. */
export const HANDOFF = Object.freeze(Object.fromEntries(Object.entries(ROWS).map(([id, rows]) => [id, Object.freeze(rows.map(([item, solver, generator, options, what]) => Object.freeze({ item, match: normItem(item), solver, generator, options: Object.freeze({ ...options }), what })))])));
/** Hand-off entry of one catalogue item of a suite (exact match on the normalised whole name), or null. */
export function handoffFor(suiteId, itemName) {
  const n = normItem(itemName), list = HANDOFF[suiteId];
  return (list && list.find((e) => e.match === n)) || null;
}
/** Exact catalogue strings of the items of a suite that are handed off. */
export const HANDOFF_ITEMS = (suiteId) => (HANDOFF[suiteId] || []).map((e) => e.item);

// ---- zip writer (stored, no compression) and reader -----------------------------------------------------
const CRC_TABLE = (() => { const t = new Uint32Array(256); for (let n = 0; n < 256; n++) { let c = n; for (let k = 0; k < 8; k++) c = c & 1 ? 0xedb88320 ^ (c >>> 1) : c >>> 1; t[n] = c >>> 0; } return t; })();
/** CRC-32 (IEEE 802.3) of a byte array. */
export function crc32(bytes) {
  let c = 0xffffffff;
  for (let i = 0; i < bytes.length; i++) c = CRC_TABLE[(c ^ bytes[i]) & 0xff] ^ (c >>> 8);
  return (c ^ 0xffffffff) >>> 0;
}
/**
 * Zip archive with stored (uncompressed) entries. files: [{ path, text | bytes, mode }]; shell scripts (Allrun, Allclean,
 * *.sh, *.py) get the executable bit. `date` fixes the time stamp so that the archive is reproducible. Returns a Uint8Array.
 */
export function zipStore(files, { date = new Date(Date.UTC(2026, 0, 1)), root = '' } = {}) {
  const enc = new TextEncoder(), parts = [], central = [];
  const dosTime = ((date.getUTCHours() & 31) << 11) | ((date.getUTCMinutes() & 63) << 5) | ((date.getUTCSeconds() >> 1) & 31), dosDate = ((Math.max(0, date.getUTCFullYear() - 1980) & 127) << 9) | (((date.getUTCMonth() + 1) & 15) << 5) | (date.getUTCDate() & 31);
  const seen = new Set();
  let offset = 0;
  if (files.length > 65000) throw new Error('Too many files for one archive.');
  for (const f of files) {
    const path = (root ? root.replace(/\/+$/, '') + '/' : '') + String(f.path).replace(/\\/g, '/').replace(/^\/+/, '');
    if (!path || path.split('/').some((p) => p === '..' || p === '')) throw new Error('Unsafe path in archive: ' + f.path);
    if (seen.has(path)) throw new Error('Duplicate path in archive: ' + path);
    seen.add(path);
    const name = enc.encode(path), data = f.bytes instanceof Uint8Array ? f.bytes : enc.encode(String(f.text ?? '')), crc = crc32(data);
    const exec = f.mode ? f.mode === 0o755 : /(^|\/)(Allrun|Allclean|Allrun\.[\w-]+|run|clean)$|\.(sh|py)$/.test(path), attr = ((exec ? 0o100755 : 0o100644) << 16) >>> 0;
    const lh = new Uint8Array(30 + name.length), dv = new DataView(lh.buffer);
    dv.setUint32(0, 0x04034b50, true); dv.setUint16(4, 20, true); dv.setUint16(6, 0x0800, true); dv.setUint16(8, 0, true); dv.setUint16(10, dosTime, true); dv.setUint16(12, dosDate, true);
    dv.setUint32(14, crc, true); dv.setUint32(18, data.length, true); dv.setUint32(22, data.length, true); dv.setUint16(26, name.length, true); dv.setUint16(28, 0, true); lh.set(name, 30);
    const ch = new Uint8Array(46 + name.length), cv = new DataView(ch.buffer);
    cv.setUint32(0, 0x02014b50, true); cv.setUint16(4, (3 << 8) | 20, true); cv.setUint16(6, 20, true); cv.setUint16(8, 0x0800, true); cv.setUint16(10, 0, true); cv.setUint16(12, dosTime, true); cv.setUint16(14, dosDate, true);
    cv.setUint32(16, crc, true); cv.setUint32(20, data.length, true); cv.setUint32(24, data.length, true); cv.setUint16(28, name.length, true); cv.setUint32(38, attr, true); cv.setUint32(42, offset, true); ch.set(name, 46);
    parts.push(lh, data); central.push(ch); offset += lh.length + data.length;
    if (offset > BRIDGE_LIMITS.zipBytes) throw new Error('Archive would be larger than ' + BRIDGE_LIMITS.zipBytes / 1e6 + ' MB.');
  }
  const cdSize = central.reduce((s, c) => s + c.length, 0), end = new Uint8Array(22), ev = new DataView(end.buffer);
  ev.setUint32(0, 0x06054b50, true); ev.setUint16(8, files.length, true); ev.setUint16(10, files.length, true); ev.setUint32(12, cdSize, true); ev.setUint32(16, offset, true);
  const out = new Uint8Array(offset + cdSize + 22);
  let p = 0;
  for (const part of [...parts, ...central, end]) { out.set(part, p); p += part.length; }
  return out;
}
/** Read the stored entries of a zip archive (central directory walk, CRC verified). Returns [{ path, bytes, text, crcOk, method }]. */
export function unzipStored(u8) {
  if (!(u8 instanceof Uint8Array) || u8.length < 22 || u8.length > BRIDGE_LIMITS.zipBytes) throw new Error('Not a zip archive, or too large.');
  const dv = new DataView(u8.buffer, u8.byteOffset, u8.byteLength), dec = new TextDecoder(), out = [];
  let eocd = -1;
  for (let i = u8.length - 22; i >= Math.max(0, u8.length - 66000); i--) if (dv.getUint32(i, true) === 0x06054b50) { eocd = i; break; }
  if (eocd < 0) throw new Error('Zip end record not found.');
  const n = dv.getUint16(eocd + 10, true);
  let p = dv.getUint32(eocd + 16, true);
  for (let k = 0; k < n && k < BRIDGE_LIMITS.files * 10; k++) {
    if (p + 46 > u8.length || dv.getUint32(p, true) !== 0x02014b50) throw new Error('Corrupt zip directory.');
    const method = dv.getUint16(p + 10, true), crc = dv.getUint32(p + 16, true), csize = dv.getUint32(p + 20, true), nl = dv.getUint16(p + 28, true), el = dv.getUint16(p + 30, true), cl = dv.getUint16(p + 32, true), lho = dv.getUint32(p + 42, true);
    const path = dec.decode(u8.subarray(p + 46, p + 46 + nl));
    if (lho + 30 > u8.length || dv.getUint32(lho, true) !== 0x04034b50) throw new Error('Corrupt zip entry: ' + path);
    const start = lho + 30 + dv.getUint16(lho + 26, true) + dv.getUint16(lho + 28, true);
    if (start + csize > u8.length) throw new Error('Truncated zip entry: ' + path);
    const bytes = u8.subarray(start, start + csize);
    out.push({ path, bytes, method, crcOk: method === 0 ? crc32(bytes) === crc : null, text: method === 0 ? dec.decode(bytes) : null, mode: dv.getUint32(p + 38, true) >>> 16 });
    p += 46 + nl + el + cl;
  }
  return out;
}

// ---- the case as plain data --------------------------------------------------------------------------------
/** Reference numbers used when a generator is called without data (the reference tie-back near mid-line). */
export const DEFAULT_BRIDGE_CASE = Object.freeze({
  name: BASE.name, x: 9000, z: -1310, angleDeg: 0.29, D: BASE.idMm / 1000, wt: BASE.wtMm / 1000, roughness: BASE.roughUm * 1e-6, lineLength: 19900,
  P: 85.7, T: 59.1, holdup: 0.42, vsl: 1.04, vsg: 0.79, regime: 'slug', dpdx: 38, tauW: 5.9,
  rhoL: 795.5, rhoG: 72.1, muL: 2.99e-3, muG: 1.41e-5, sigma: 0.01, rhoO: 760, rhoW: 1010, muO: 2.2e-3, muW: 5e-4, wcut: 0.2, cpL: 2500, cpG: 2600, kL: 0.2, kG: 0.045, mwG: 21,
  tAmb: 4, pExt: 133, tInstall: 4,
  steel: Object.freeze({ grade: 'API 5L X65', E: BASE.E * 1e6, nu: BASE.poisson, rho: BASE.rhoSteel, alphaT: BASE.alphaT, smys: BASE.smys * 1e6 }),
  sand: Object.freeze({ rateKgD: 50, dUm: 150, rho: 2650 }),
  slug: Object.freeze({ freq: 0.03, length: 60, velocity: 2.9, holdupBody: 0.9, forceN: 4500 }),
  riser: Object.freeze({ baseX: BASE.riserBaseX, height: BASE.riserHeight, flowlineLength: 400, flowlineAngleDeg: -1.5 }),
  series: null, source: 'reference case',
});

/**
 * Plain-data picture of the case at one location x (m along the line) — what every generator takes as input.
 * ctx = { fluid, site, outputs, inputs } (any part may be missing). Uses the flow suite's profile when it has been
 * run, otherwise the kernel estimate. Returns the fields of DEFAULT_BRIDGE_CASE filled from the case.
 */
export function bridgeCase(ctx = {}, opt = {}) {
  const line = caseLine(ctx), pic = flowPicture(ctx), fm = fluidModel(ctx), o = ctx.outputs || {}, n = pic.x.length;
  const x = clamp(num(opt.x, 0.5 * pic.x[n - 1]), pic.x[0], pic.x[n - 1]), at = (k) => interp1(pic.x, pic[k], x);
  let i = 0; while (i < n - 2 && pic.x[i + 1] <= x) i++;
  const dx = pic.x[i + 1] - pic.x[i], angle = Math.atan2(pic.z[i + 1] - pic.z[i], Math.max(dx, 1e-9)), P = at('P'), T = at('T'), f = fm.at(P, T), z = at('z');
  const vsl = Math.max(at('vsl'), 1e-4), vsg = Math.max(at('vsg'), 1e-4), D = line.id, A = (Math.PI * D * D) / 4;
  const su = slugUnit({ vsl, vsg, rhoL: f.rhoL, rhoG: f.rhoG, muL: f.muL, muG: f.muG, D, theta: angle }), fs = o.flow?.slug, mat = o.net?.material || {};
  const slug = { freq: pos(fs?.freq, su.freq), length: pos(fs?.length, su.length), velocity: pos(fs?.velocity, su.vt), holdupBody: pos(fs?.holdupBody, su.holdupSlug) };
  slug.forceN = isNum(o.integ?.slugForce) && o.integ.slugForce > 0 ? o.integ.slugForce * 1000 : Math.SQRT2 * (slug.holdupBody * f.rhoL + (1 - slug.holdupBody) * f.rhoG) * A * slug.velocity ** 2;
  const si = ctx.inputs?.solids || {}, ii = ctx.inputs?.integ || {};
  const rb = line.riserBaseX, zb = interp1(line.profile.x, line.profile.z, rb), zTop = line.profile.z[line.profile.z.length - 1], xa = Math.max(line.profile.x[0], rb - 400), za = interp1(line.profile.x, line.profile.z, xa);
  const s = o.flow?.series, okS = s && Array.isArray(s.t) && s.t.length > 3 && Array.isArray(s.qLiqOut) && Array.isArray(s.qGasOut) && s.qLiqOut.length === s.t.length && s.qGasOut.length === s.t.length;
  return {
    name: String(ctx.name || 'case'), x, z, angleDeg: (angle * 180) / Math.PI, D, wt: line.wt, roughness: line.roughness, lineLength: pic.x[n - 1] - pic.x[0],
    P, T, holdup: clamp(at('holdup'), 0.01, 0.99), vsl, vsg, regime: String(pic.regime[Math.round(interp1(pic.x, pic.x.map((_, k) => k), x))] || ''), dpdx: at('dpdx'), tauW: Math.max(at('tauW'), 1e-3),
    rhoL: f.rhoL, rhoG: f.rhoG, muL: f.muL, muG: f.muG, sigma: f.sigma, rhoO: f.rhoO, rhoW: f.rhoW, muO: f.muO, muW: f.muW, wcut: f.wcut, cpL: f.cpL, cpG: f.cpG, kL: f.kL, kG: f.kG, mwG: f.mwG,
    tAmb: at('tAmb'), pExt: z < 0 ? 1.01325 + (1025 * G * -z) / 1e5 : 1.01325, tInstall: num(line.tSeabed, 4),
    steel: { grade: String(mat.grade || 'API 5L X65'), E: pos(mat.E, BASE.E) * 1e6, nu: pos(mat.poisson, BASE.poisson), rho: pos(mat.rho, BASE.rhoSteel), alphaT: pos(mat.alphaT, BASE.alphaT), smys: pos(mat.smys, BASE.smys) * 1e6 },
    sand: { rateKgD: num(ii.sandKgD ?? si.sandRate, 50), dUm: pos(ii.sandUm ?? si.sandUm, 150), rho: pos(ii.sandDensity ?? si.sandRho, 2650) },
    slug,
    riser: { baseX: rb, height: Math.max(zTop - zb, 1), flowlineLength: Math.max(rb - xa, 1), flowlineAngleDeg: (Math.atan2(zb - za, Math.max(rb - xa, 1e-9)) * 180) / Math.PI },
    series: okS ? { t: s.t.slice(0, 400), vsl: s.qLiqOut.slice(0, 400).map((q) => Math.max(q, 0) / A), vsg: s.qGasOut.slice(0, 400).map((q) => Math.max(q, 0) / A), pIn: Array.isArray(s.pIn) && s.pIn.length === s.t.length ? s.pIn.slice(0, 400) : null } : null,
    source: pic.source,
  };
}
const withCase = (opts = {}) => {
  const c = { ...DEFAULT_BRIDGE_CASE, ...(opts.case || {}) };
  for (const k of ['steel', 'sand', 'slug', 'riser']) c[k] = { ...DEFAULT_BRIDGE_CASE[k], ...((opts.case || {})[k] || {}) };
  for (const k of ['D', 'wt', 'P', 'vsl', 'vsg', 'rhoL', 'rhoG', 'muL', 'muG', 'sigma', 'tauW']) c[k] = pos(c[k], DEFAULT_BRIDGE_CASE[k]);
  c.holdup = clamp(num(c.holdup, 0.4), 0.01, 0.99); c.angleDeg = clamp(num(c.angleDeg, 0), -90, 90);
  return c;
};
const BANNER = (title) => ['/*--------------------------------*- C++ -*----------------------------------*\\', '| HydraSlug hand-off case: ' + title, '| Written for ' + OPENFOAM_FLAVOUR, '\\*---------------------------------------------------------------------------*/'];
const foamHeader = (cls, object, title, location) => [...BANNER(title), 'FoamFile', '{', '    version     2.0;', '    format      ascii;', `    class       ${cls};`, location ? `    location    "${location}";` : null, `    object      ${object};`, '}', '// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //', ''];

// ---- swept block-structured meshes (blockMesh) --------------------------------------------------------------
/** Liquid level h/D of a stratified layer that fills the area fraction `holdup` of a circular section. */
export function levelForHoldup(holdup) {
  const H = clamp(holdup, 1e-6, 1 - 1e-6), area = (hD) => { const c = 2 * hD - 1; return (Math.PI - Math.acos(c) + c * Math.sqrt(1 - c * c)) / Math.PI; };
  let lo = 0, hi = 1;
  for (let i = 0; i < 60; i++) { const m = 0.5 * (lo + hi); if (area(m) < H) lo = m; else hi = m; }
  return 0.5 * (lo + hi);
}
/**
 * Centreline stations of a path made of straight legs and circular bends, with a transported frame.
 * legs: [{ type: 'straight', L } | { type: 'bend', R, angleDeg, toward: 'up' | 'down' | 'left' | 'right' }]; `div` on a leg fixes its number of subdivisions.
 * Returns [{ c, t, e1, e2, s (arc length), kind ('straight' | 'bend' for the segment that ENDS here), mid: frame at the middle of that segment }].
 * The frame is right-handed with e1 × e2 = t; it starts as t = +x, e1 = +y, e2 = +z ("up" before gravity is resolved).
 */
export function sweepStations(legs, { maxBendStep = 45, split = 1 } = {}) {
  let f = { c: [0, 0, 0], t: [1, 0, 0], e1: [0, 1, 0], e2: [0, 0, 1], s: 0 };
  const out = [{ ...f, kind: 'start', mid: null }];
  const turn = (fr, n, Rb, psi) => { // rotate the frame about the bend axis; n is the unit vector toward the centre of curvature
    const rot = (v) => { const a = dot(v, fr.t), b = dot(v, n), rest = sub(sub(v, mul(fr.t, a)), mul(n, b)); return add(rest, add(mul(fr.t, a * Math.cos(psi) - b * Math.sin(psi)), mul(n, a * Math.sin(psi) + b * Math.cos(psi)))); };
    return { c: add(fr.c, add(mul(fr.t, Rb * Math.sin(psi)), mul(n, Rb * (1 - Math.cos(psi))))), t: rot(fr.t), e1: rot(fr.e1), e2: rot(fr.e2), s: fr.s + Rb * psi };
  };
  for (const leg of legs) {
    if (leg.type === 'bend') {
      const ang = (clamp(Math.abs(num(leg.angleDeg, 90)), 1, 180) * Math.PI) / 180, Rb = pos(leg.R, 1), nSub = leg.div > 0 ? Math.round(leg.div) : Math.max(1, Math.ceil((ang * 180) / Math.PI / maxBendStep - 1e-9)) * split;
      const n = leg.toward === 'down' ? mul(f.e2, -1) : leg.toward === 'left' ? f.e1 : leg.toward === 'right' ? mul(f.e1, -1) : f.e2, f0 = f;
      for (let k = 1; k <= nSub; k++) { f = turn(f0, n, Rb, (ang * k) / nSub); out.push({ ...f, kind: 'bend', R: Rb, mid: turn(f0, n, Rb, (ang * (k - 0.5)) / nSub), dPsi: ang / nSub, n }); }
    } else {
      const L = pos(leg.L, 1), nSub = leg.div > 0 ? Math.round(leg.div) : Math.max(1, Math.round(split * (leg.parts || 1)));
      for (let k = 1; k <= nSub; k++) { const prev = f; f = { ...f, c: add(prev.c, mul(prev.t, L / nSub)), s: prev.s + L / nSub }; out.push({ ...f, kind: 'straight', mid: { ...f, c: add(prev.c, mul(prev.t, L / nSub / 2)), s: prev.s + L / nSub / 2 } }); }
    }
  }
  return out;
}
const sectionPoint = (fr, r, phi) => add(fr.c, add(mul(fr.e1, r * Math.cos(phi)), mul(fr.e2, r * Math.sin(phi))));

/**
 * Block topology of a pipe swept along `legs` with an O-grid (butterfly) cross-section: one square core block and four
 * wall blocks per segment. spec: { D, legs, nc (cells along one side of the core = cells per quarter of the circumference),
 * nr (radial cells in the wall blocks), wallRatio (wall cell / core-side cell of the wall blocks), dz (axial cell length),
 * core (corner radius of the core square / pipe radius) }.
 * Returns { vertices, blocks: [{ v[8], n[3], grading[3], kind }], edges: [{ a, b, mid }], patches: { inlet, outlet, wall },
 *   stations, nCells, length, volume (exact, m³), area }.
 */
export function ogridMesh(spec = {}) {
  const D = pos(spec.D, 0.254), Rp = D / 2, rc = clamp(num(spec.core, 0.55), 0.3, 0.8) * Rp, nc = Math.max(2, Math.round(num(spec.nc, 8))), nr = Math.max(1, Math.round(num(spec.nr, 6))), dz = pos(spec.dz, D / 8);
  const legs = spec.legs && spec.legs.length ? spec.legs : [{ type: 'straight', L: 20 * D }], st = sweepStations(legs), phi = [0, 1, 2, 3].map((k) => Math.PI / 4 + (k * Math.PI) / 2);
  const vertices = [], blocks = [], edges = [], patches = { inlet: [], outlet: [], wall: [] }, bulge = 0.884; // mid-point radius of the curved core sides relative to the corner radius
  const I = (s, k) => 8 * s + k, O = (s, k) => 8 * s + 4 + k;
  st.forEach((fr) => { for (const r of [rc, Rp]) for (let k = 0; k < 4; k++) vertices.push(sectionPoint(fr, r, phi[k])); });
  const faces = (s) => [[I(s, 2), I(s, 3), I(s, 0), I(s, 1)], ...[0, 1, 2, 3].map((k) => [I(s, k), O(s, k), O(s, (k + 1) % 4), I(s, (k + 1) % 4)])];
  let nCells = 0;
  st.forEach((fr, s) => {
    for (let k = 0; k < 4; k++) { // curved sides of the section at every station
      const pm = phi[k] + Math.PI / 4;
      edges.push({ a: O(s, k), b: O(s, (k + 1) % 4), mid: sectionPoint(fr, Rp, pm) }, { a: I(s, k), b: I(s, (k + 1) % 4), mid: sectionPoint(fr, rc * bulge, pm) });
    }
    if (s === 0) return;
    const prev = st[s - 1], len = fr.s - prev.s, nz = Math.max(1, Math.round(len / dz)), f0 = faces(s - 1), f1 = faces(s);
    blocks.push({ v: [...f0[0], ...f1[0]], n: [nc, nc, nz], grading: [1, 1, 1], kind: 'core' });
    for (let k = 0; k < 4; k++) blocks.push({ v: [...f0[k + 1], ...f1[k + 1]], n: [nr, nc, nz], grading: [pos(spec.wallRatio, 1), 1, 1], kind: 'wall' });
    nCells += (nc * nc + 4 * nr * nc) * nz;
    if (fr.kind === 'bend') for (const r of [rc, Rp]) for (let k = 0; k < 4; k++) edges.push({ a: 8 * (s - 1) + (r === Rp ? 4 : 0) + k, b: 8 * s + (r === Rp ? 4 : 0) + k, mid: sectionPoint(fr.mid, r, phi[k]) });
    for (let k = 0; k < 4; k++) patches.wall.push([O(s - 1, k), O(s - 1, (k + 1) % 4), O(s, (k + 1) % 4), O(s, k)]);
  });
  patches.inlet = faces(0); patches.outlet = faces(st.length - 1);
  const length = st[st.length - 1].s, area = Math.PI * Rp * Rp;
  return { type: 'ogrid', D, vertices, blocks, edges, patches, stations: st, nCells, length, area, volume: area * length, nc, nr, sectionCells: nc * nc + 4 * nr * nc };
}
/**
 * Planar (one cell thick) channel of height D swept along `legs` in the vertical plane: a 2-D model of a pipeline–riser.
 * spec: { D, legs, nD (cells across the height), dz, width }. Patches: inlet, outlet, wall, frontAndBack (empty).
 */
export function slotMesh(spec = {}) {
  const D = pos(spec.D, 0.254), w = pos(spec.width, D), nD = Math.max(2, Math.round(num(spec.nD, 16))), dz = pos(spec.dz, D / 8), legs = spec.legs && spec.legs.length ? spec.legs : [{ type: 'straight', L: 20 * D }], st = sweepStations(legs, { maxBendStep: 30 });
  const vertices = [], blocks = [], edges = [], patches = { inlet: [], outlet: [], wall: [], frontAndBack: [] }, corner = [[-1, -1], [1, -1], [1, 1], [-1, 1]];
  st.forEach((fr) => corner.forEach(([a, b]) => vertices.push(add(fr.c, add(mul(fr.e1, (a * w) / 2), mul(fr.e2, (b * D) / 2))))));
  let nCells = 0;
  st.forEach((fr, s) => {
    if (s === 0) return;
    const a = 4 * (s - 1), b = 4 * s, nz = Math.max(1, Math.round((fr.s - st[s - 1].s) / dz));
    blocks.push({ v: [a, a + 1, a + 2, a + 3, b, b + 1, b + 2, b + 3], n: [1, nD, nz], grading: [1, 1, 1], kind: 'slot' }); nCells += nD * nz;
    if (fr.kind === 'bend') corner.forEach(([p, q], k) => edges.push({ a: a + k, b: b + k, mid: add(fr.mid.c, add(mul(fr.mid.e1, (p * w) / 2), mul(fr.mid.e2, (q * D) / 2))) }));
    patches.wall.push([a, a + 1, b + 1, b], [a + 3, a + 2, b + 2, b + 3]); patches.frontAndBack.push([a, a + 3, b + 3, b], [a + 1, a + 2, b + 2, b + 1]);
  });
  patches.inlet = [[0, 1, 2, 3]]; const e = 4 * (st.length - 1); patches.outlet = [[e, e + 1, e + 2, e + 3]];
  const length = st[st.length - 1].s;
  return { type: 'slot', D, width: w, vertices, blocks, edges, patches, stations: st, nCells, length, area: D * w, volume: D * w * length, nD, sectionCells: nD };
}
/** Geometric and topological checks of a block mesh description: positive corner Jacobians, faces shared by two blocks or on exactly one patch, consistent cell counts across shared faces. */
export function checkBlockMesh(m) {
  const issues = [], key = (f) => f.slice().sort((a, b) => a - b).join(','), count = new Map(), F = [[0, 3, 2, 1], [4, 5, 6, 7], [0, 1, 5, 4], [3, 7, 6, 2], [0, 4, 7, 3], [1, 2, 6, 5]];
  let minJ = Infinity;
  const edgeN = new Map(), EDG = [[0, 1, 0], [3, 2, 0], [7, 6, 0], [4, 5, 0], [0, 3, 1], [1, 2, 1], [5, 6, 1], [4, 7, 1], [0, 4, 2], [1, 5, 2], [2, 6, 2], [3, 7, 2]];
  m.blocks.forEach((b, bi) => {
    if (b.v.some((v) => !(v >= 0 && v < m.vertices.length))) issues.push(`block ${bi} references a missing vertex`);
    const P = b.v.map((v) => m.vertices[v]), XYZ = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], at = (i, j, l) => P[XYZ.findIndex((q) => q[0] === i && q[1] === j && q[2] === l)];
    XYZ.forEach(([i, j, l], k) => { const j3 = dot(cross(sub(at(1, j, l), at(0, j, l)), sub(at(i, 1, l), at(i, 0, l))), sub(at(i, j, 1), at(i, j, 0))); minJ = Math.min(minJ, j3); if (!(j3 > 0)) issues.push(`block ${bi} corner ${k} has a non-positive Jacobian`); });
    F.forEach((f) => { const kf = key(f.map((i) => b.v[i])); count.set(kf, (count.get(kf) || 0) + 1); });
    EDG.forEach(([a, c, dir]) => { const ke = [b.v[a], b.v[c]].sort((x, y) => x - y).join(','); if (edgeN.has(ke) && edgeN.get(ke) !== b.n[dir]) issues.push(`cell count mismatch along edge ${ke}`); edgeN.set(ke, b.n[dir]); });
  });
  const patchFaces = new Map();
  for (const [name, faces] of Object.entries(m.patches)) for (const f of faces) { const kf = key(f); if (patchFaces.has(kf)) issues.push(`face ${kf} is on two patches`); patchFaces.set(kf, name); if (!count.has(kf)) issues.push(`patch ${name} face ${kf} is not a block face`); }
  let internal = 0, boundary = 0;
  for (const [kf, c] of count) { if (c === 2) { internal++; if (patchFaces.has(kf)) issues.push(`internal face ${kf} is listed on a patch`); } else if (c === 1) { boundary++; if (!patchFaces.has(kf)) issues.push(`boundary face ${kf} is on no patch`); } else issues.push(`face ${kf} is shared by ${c} blocks`); }
  for (const e of m.edges) if (!(e.a >= 0 && e.a < m.vertices.length && e.b >= 0 && e.b < m.vertices.length) || !e.mid.every(isNum)) issues.push('edge with a missing vertex or mid-point');
  return { ok: !issues.length, issues, minJacobian: minJ, internalFaces: internal, boundaryFaces: boundary };
}
/** Volume of the meshed region from the block description itself: every curved side of a section is integrated as the circular arc through its three points; bends by Pappus' theorem (both section shapes have their centroid on the centreline). */
export function blockMeshVolume(m, nArc = 96) {
  const st = m.stations, outer = m.type === 'ogrid' ? [4, 5, 6, 7] : [0, 1, 2, 3], per = m.type === 'ogrid' ? 8 : 4;
  const arcPts = (A, M, B) => { // points of the circular arc A–M–B, B excluded
    const u = sub(M, A), v = sub(B, A), w = cross(u, v), ww = dot(w, w);
    if (ww < 1e-24) return [A];
    const c = add(A, mul(add(mul(cross(v, w), dot(u, u)), mul(cross(w, u), dot(v, v))), 1 / (2 * ww))), ra = sub(A, c), r = hyp(ra), ex = mul(ra, 1 / r), ey = cross(mul(w, 1 / Math.sqrt(ww)), ex);
    let ang = Math.atan2(dot(sub(B, c), ey), dot(sub(B, c), ex)); if (ang <= 0) ang += 2 * Math.PI;
    return Array.from({ length: nArc }, (_, i) => add(c, add(mul(ex, r * Math.cos((ang * i) / nArc)), mul(ey, r * Math.sin((ang * i) / nArc)))));
  };
  const sectionArea = (s) => {
    const fr = st[s], ring = [];
    outer.forEach((k, i) => { const a = per * s + k, b = per * s + outer[(i + 1) % 4], e = m.edges.find((x) => (x.a === a && x.b === b) || (x.a === b && x.b === a)); ring.push(...(e ? arcPts(m.vertices[a], e.mid, m.vertices[b]) : [m.vertices[a]])); });
    const q = ring.map((p) => [dot(sub(p, fr.c), fr.e1), dot(sub(p, fr.c), fr.e2)]);
    let A = 0; q.forEach((p, i) => { const n = q[(i + 1) % q.length]; A += (p[0] * n[1] - n[0] * p[1]) / 2; });
    return Math.abs(A);
  };
  let V = 0;
  for (let s = 1; s < st.length; s++) { const A = 0.5 * (sectionArea(s - 1) + sectionArea(s)), fr = st[s]; V += fr.kind === 'bend' ? A * fr.dPsi * fr.R : A * (fr.s - st[s - 1].s); }
  return V;
}

/** blockMeshDict text of a mesh description from ogridMesh / slotMesh. */
export function blockMeshDict(m, title = 'pipe section', inletName = 'inletGas') {
  const patchType = { wall: 'wall', frontAndBack: 'empty' };
  return text([foamHeader('dictionary', 'blockMeshDict', title, 'system'), 'scale   1;', '', 'vertices', '(',
    m.vertices.map((v, i) => `    ${vec(v, 10)} // ${i}`), ');', '', 'blocks', '(',
    m.blocks.map((b) => `    hex (${b.v.join(' ')}) (${b.n.join(' ')}) simpleGrading (${b.grading.map((g) => ff(g, 6)).join(' ')})`), ');', '', 'edges', '(',
    m.edges.map((e) => `    arc ${e.a} ${e.b} ${vec(e.mid, 10)}`), ');', '', 'boundary', '(',
    Object.entries(m.patches).map(([name, faces]) => [`    ${name === 'inlet' ? inletName : name}`, '    {', `        type ${patchType[name] || 'patch'};`, '        faces', '        (', faces.map((f) => `            (${f.join(' ')})`), '        );', '    }']),
    ');', '', 'mergePatchPairs', '(', ');']);
}

// ---- mesh sizing and cost ------------------------------------------------------------------------------------
export const TURBULENCE = Object.freeze({
  kOmegaSST: { label: 'k–ω SST (RANS)', type: 'RAS', fields: ['k', 'omega', 'nut'], yPlus: 50, aspect: 2 },
  kEpsilon: { label: 'Standard k–ε (RANS)', type: 'RAS', fields: ['k', 'epsilon', 'nut'], yPlus: 50, aspect: 2 },
  RNGkEpsilon: { label: 'RNG k–ε (RANS)', type: 'RAS', fields: ['k', 'epsilon', 'nut'], yPlus: 50, aspect: 2 },
  realizableKE: { label: 'Realizable k–ε (RANS)', type: 'RAS', fields: ['k', 'epsilon', 'nut'], yPlus: 50, aspect: 2 },
  LRR: { label: 'Reynolds-stress transport, LRR (RANS)', type: 'RAS', fields: ['R', 'k', 'epsilon', 'nut'], yPlus: 50, aspect: 2 },
  SpalartAllmaras: { label: 'Spalart–Allmaras (RANS)', type: 'RAS', fields: ['nuTilda', 'nut'], yPlus: 30, aspect: 2 },
  WALE: { label: 'LES, WALE sub-grid model', type: 'LES', fields: ['nut'], yPlus: 30, aspect: 1 },
  Smagorinsky: { label: 'LES, Smagorinsky sub-grid model', type: 'LES', fields: ['nut'], yPlus: 30, aspect: 1 },
  kEqn: { label: 'LES, one-equation k sub-grid model', type: 'LES', fields: ['k', 'nut'], yPlus: 30, aspect: 1 },
  SpalartAllmarasDES: { label: 'DES (Spalart–Allmaras)', type: 'LES', fields: ['nuTilda', 'nut'], yPlus: 30, aspect: 1.5 },
  SpalartAllmarasDDES: { label: 'Delayed DES (Spalart–Allmaras)', type: 'LES', fields: ['nuTilda', 'nut'], yPlus: 30, aspect: 1.5 },
  SpalartAllmarasIDDES: { label: 'IDDES (Spalart–Allmaras)', type: 'LES', fields: ['nuTilda', 'nut'], yPlus: 30, aspect: 1.5 },
  kOmegaSSTDES: { label: 'DES (k–ω SST)', type: 'LES', fields: ['k', 'omega', 'nut'], yPlus: 30, aspect: 1.5 },
  kOmegaSSTIDDES: { label: 'IDDES (k–ω SST)', type: 'LES', fields: ['k', 'omega', 'nut'], yPlus: 30, aspect: 1.5 },
  DNS: { label: 'No model — DNS resolution', type: 'laminar', fields: [], yPlus: 1, aspect: 1 },
  laminar: { label: 'Laminar', type: 'laminar', fields: [], yPlus: 1, aspect: 2 },
});
/** Core-seconds per cell and time step used for the run-cost estimates (measured with interFoam on a workstation core; the read-me repeats it). */
export const COST_PER_CELL_STEP = 5e-5;

/**
 * Mesh plan of a pipe section from a target y+ and a number of cells per diameter.
 * Returns { nc, nr, wallRatio, y1 (first cell height, m), yPlus (achieved), dz, cells, dt, steps, coreHours, uTau, reTau, dns: { cells, coreHours } , notes[] }.
 */
export function pipeMeshPlan(c, o = {}) {
  const tm = TURBULENCE[o.turbulence] || TURBULENCE.kOmegaSST, D = c.D, Rp = D / 2, notes = [];
  const nPerD = Math.round(clamp(num(o.cellsPerDiameter, tm.type === 'RAS' ? 24 : 40), 8, 400)), yPlus = clamp(num(o.yPlus, tm.yPlus), 0.2, 300), length = pos(o.length, 20 * D);
  const vm = c.vsl + c.vsg, uL = c.vsl / c.holdup, uG = c.vsg / (1 - c.holdup), uMax = 1.5 * Math.max(uL, uG, vm);
  // viscous length scale in each phase; the thinner one sets the first cell
  // viscous length scale μ/√(τw ρ) of each layer (single-phase friction of the layer, the liquid at least at the 1-D wall shear); the thinner one sets the first cell
  const tauOf = (rho, mu, u) => 0.5 * rho * u * u * 0.046 * Math.max((rho * u * D) / mu, 2100) ** -0.2, tauL = Math.max(tauOf(c.rhoL, c.muL, uL), c.tauW), tauG = tauOf(c.rhoG, c.muG, uG);
  const lv = [c.muL / Math.sqrt(tauL * c.rhoL), c.muG / Math.sqrt(tauG * c.rhoG)], lMin = Math.min(...lv), uTau = Math.sqrt(c.tauW / (c.holdup * c.rhoL + (1 - c.holdup) * c.rhoG));
  let nc = Math.max(4, Math.round(0.4 * nPerD)), nr = Math.max(2, Math.round((nPerD - nc) / 2));
  const Lr = Rp * (1 - 0.55 * 0.9), core = (Math.SQRT2 * 0.55 * Rp) / nc; // radial length of the wall blocks, cell size of the core
  const dz = pos(o.dz, (tm.aspect * D) / nPerD), side = Math.max(dz, (Math.PI * D) / (4 * nc)), maxAspect = 1000; // checkMesh fails cells above this aspect ratio
  let y1 = Math.min(yPlus * lMin, Lr / nr), ratio = 1, growth = 1;
  if (side / y1 > maxAspect) { y1 = side / maxAspect; notes.push(`The first cell was thickened to ${ff(y1 * 1e6, 3)} µm (y+ ≈ ${ff(y1 / lMin, 3)}) to keep the cell aspect ratio at ${maxAspect}; use more cells per diameter to reach the target y+ of ${ff(yPlus, 3)}.`); }
  const solve = (n) => { // geometric growth rate g with y1 (g^n - 1)/(g - 1) = Lr
    if (y1 * n >= Lr * 0.999) return 1;
    let lo = 1, hi = 5; for (let i = 0; i < 80; i++) { const g = 0.5 * (lo + hi); if ((y1 * (g ** n - 1)) / (g - 1) < Lr) lo = g; else hi = g; } return 0.5 * (lo + hi);
  };
  growth = solve(nr);
  while (growth > 1.25 && nr < 400) { nr++; growth = solve(nr); }
  if (growth > 1.25) { notes.push('The wall layer cannot be reached with a growth rate below 1.25 within 400 radial cells; the first cell was thickened.'); while (growth > 1.25) { y1 *= 1.1; growth = solve(nr); } }
  ratio = growth ** (nr - 1);
  const nz = Math.max(1, Math.round(length / dz)), cells = (nc * nc + 4 * nr * nc) * nz, maxCo = pos(o.maxCo, tm.type === 'RAS' ? 1 : 0.5);
  const dt = (maxCo * Math.min(dz, core)) / Math.max(uMax, 1e-6), endTime = pos(o.endTime, (pos(o.flowThroughs, 5) * length) / Math.max(vm, 1e-6)), steps = Math.ceil(endTime / dt), coreHours = (cells * steps * COST_PER_CELL_STEP) / 3600;
  // DNS estimate for the same length: Δr+ ≈ 1 at the wall growing to 5, RΔθ+ ≈ 5, Δz+ ≈ 10 (pipe DNS practice)
  const reTau = Rp / lMin, dnsCells = Math.max(1, ((2 * Math.PI * reTau) / 5) * (reTau / 3) * (length / lMin / 10)) * 0.6, dnsDt = (0.3 * 10 * lMin) / Math.max(uMax, 1e-6), dnsHours = (dnsCells * (endTime / dnsDt) * COST_PER_CELL_STEP) / 3600;
  const aspect = side / y1;
  if (aspect > 400) notes.push(`Wall cells have an aspect ratio of about ${Math.round(aspect)}: acceptable with wall functions, but the pressure solver converges more slowly.`);
  if (tm.type === 'LES' && y1 / lMin > 5) notes.push(`First-cell y+ ≈ ${ff(y1 / lMin, 3)}: this is a wall-modelled LES (wall functions carry the wall layer). A wall-resolved LES needs y+ ≈ 1.`);
  if (o.turbulence === 'DNS') notes.push(`Full DNS resolution of this section needs about ${dnsCells.toExponential(2)} cells and ${dnsHours.toExponential(2)} core-hours; the mesh written here uses the resolution you asked for and is a coarse DNS (under-resolved) unless it reaches that count.`);
  return { nc, nr, nPerD, growth, aspect, wallRatio: 1 / ratio, y1, yPlus: y1 / lMin, dz, nz, length, cells, sectionCells: nc * nc + 4 * nr * nc, maxCo, dt, endTime, steps, coreHours, uTau, uMax, reTau, viscousLength: lMin, dns: { cells: dnsCells, coreHours: dnsHours }, notes };
}

// ---- OpenFOAM: shared writers ---------------------------------------------------------------------------------
const dictFile = (object, location, title, body) => text([foamHeader('dictionary', object, title, location), body]);
const fieldFile = (cls, name, dims, internal, bf, title) => text([foamHeader(cls, name, title, '0'), `dimensions      ${dims};`, '', `internalField   ${internal};`, '', 'boundaryField', '{', Object.entries(bf).map(([patch, lines]) => [`    ${patch}`, '    {', lines.map((l) => '        ' + l), '    }']), '}']);
const INLETS = '"inlet.*"';
/** Turbulence set-up: constant/turbulenceProperties plus the 0/ fields the chosen model needs. */
function turbulenceSetup(key, { U, D, nu, title, empty = false, suffix = '' }) {
  const tm = TURBULENCE[key] || TURBULENCE.kOmegaSST, name = key === 'DNS' ? 'laminar' : TURBULENCE[key] ? key : 'kOmegaSST', k = Math.max(1.5 * (0.05 * U) ** 2, 1e-8), l = 0.07 * D, eps = (0.09 ** 0.75 * k ** 1.5) / l, om = Math.sqrt(k) / (0.09 ** 0.25 * l), files = [];
  const props = tm.type === 'laminar' ? ['simulationType  laminar;'] : tm.type === 'RAS' ? ['simulationType  RAS;', '', 'RAS', '{', `    RASModel        ${name};`, '    turbulence      on;', '    printCoeffs     on;', '}'] : ['simulationType  LES;', '', 'LES', '{', `    LESModel        ${name};`, '    turbulence      on;', '    printCoeffs     on;', ...(/IDDES$/.test(name) ? ['    delta           IDDESDelta;', '    IDDESDeltaCoeffs', '    {', '        hmax            maxDeltaxyzCubeRoot;', '        maxDeltaxyzCubeRootCoeffs', '        {', '        }', '    }'] : ['    delta           cubeRootVol;', '    cubeRootVolCoeffs', '    {', '        deltaCoeff      1;', '    }']), '}'];
  const e = empty ? { frontAndBack: ['type            empty;'] } : {}, fv = (v) => [`type            fixedValue;`, `value           uniform ${v};`], io = (v) => ['type            inletOutlet;', `inletValue      uniform ${v};`, `value           uniform ${v};`], wf = (type, v) => [`type            ${type};`, `value           uniform ${v};`];
  const kWall = tm.type === 'RAS' ? 'nutkWallFunction' : 'nutUSpaldingWallFunction';
  for (const f of tm.fields) {
    if (f === 'k') files.push({ path: `0.orig/k${suffix}`, text: fieldFile('volScalarField', 'k' + suffix, '[0 2 -2 0 0 0 0]', `uniform ${ff(k, 5)}`, { [INLETS]: fv(ff(k, 5)), outlet: io(ff(k, 5)), wall: wf('kqRWallFunction', ff(k, 5)), ...e }, title) });
    if (f === 'omega') files.push({ path: `0.orig/omega${suffix}`, text: fieldFile('volScalarField', 'omega' + suffix, '[0 0 -1 0 0 0 0]', `uniform ${ff(om, 5)}`, { [INLETS]: fv(ff(om, 5)), outlet: io(ff(om, 5)), wall: wf('omegaWallFunction', ff(om, 5)), ...e }, title) });
    if (f === 'epsilon') files.push({ path: `0.orig/epsilon${suffix}`, text: fieldFile('volScalarField', 'epsilon' + suffix, '[0 2 -3 0 0 0 0]', `uniform ${ff(eps, 5)}`, { [INLETS]: fv(ff(eps, 5)), outlet: io(ff(eps, 5)), wall: wf('epsilonWallFunction', ff(eps, 5)), ...e }, title) });
    if (f === 'nuTilda') files.push({ path: `0.orig/nuTilda${suffix}`, text: fieldFile('volScalarField', 'nuTilda' + suffix, '[0 2 -1 0 0 0 0]', `uniform ${ff(4 * nu, 5)}`, { [INLETS]: fv(ff(4 * nu, 5)), outlet: io(ff(4 * nu, 5)), wall: fv('0'), ...e }, title) });
    if (f === 'R') { const r = `(${ff((2 * k) / 3, 5)} 0 0 ${ff((2 * k) / 3, 5)} 0 ${ff((2 * k) / 3, 5)})`; files.push({ path: `0.orig/R${suffix}`, text: fieldFile('volSymmTensorField', 'R' + suffix, '[0 2 -2 0 0 0 0]', `uniform ${r}`, { [INLETS]: fv(r), outlet: io(r), wall: wf('kqRWallFunction', r), ...e }, title) }); }
    if (f === 'nut') files.push({ path: `0.orig/nut${suffix}`, text: fieldFile('volScalarField', 'nut' + suffix, '[0 2 -1 0 0 0 0]', 'uniform 0', { [INLETS]: ['type            calculated;', 'value           uniform 0;'], outlet: ['type            calculated;', 'value           uniform 0;'], wall: wf(key === 'SpalartAllmaras' ? 'nutUSpaldingWallFunction' : kWall, '0'), ...e }, title) });
  }
  return { tm, name, props, files, k, eps, om };
}
const table1 = (t, y, sig = 6) => `table (${t.map((ti, i) => `(${ff(ti, 7)} ${ff(y[i], sig)})`).join(' ')})`;
/** Function objects written by every VOF case: holdup (volume and sections), pressure at sections and probes, wall forces. */
function vofFunctions({ mesh, stationsAt, every, forces = true }) {
  const out = ['functions', '{',
    '    holdupVolume', '    {', '        type            volFieldValue;', '        libs            (fieldFunctionObjects);', '        operation       volAverage;', '        fields          (alpha.liquid);', '        regionType      all;', '        writeFields     false;', '        log             false;', '        writeControl    timeStep;', `        writeInterval   ${every};`, '    }'];
  stationsAt.forEach((p, i) => out.push(
    `    section${i + 1}`, '    {', '        type            surfaceFieldValue;', '        libs            (fieldFunctionObjects);', '        regionType      sampledSurface;', `        name            section${i + 1};`,
    '        sampledSurfaceDict', '        {', '            type        plane;', `            point       ${vec(p.c)};`, `            normal      ${vec(p.t)};`, `            bounds      ${vec(sub(p.c, [mesh.D, mesh.D, mesh.D].map((v) => 0.8 * v)))} ${vec(add(p.c, [mesh.D, mesh.D, mesh.D].map((v) => 0.8 * v)))};`, '            triangulate false;', '        }',
    '        operation       areaAverage;', '        fields          (alpha.liquid p);', '        surfaceFormat   none;', '        writeFields     false;', '        log             false;', '        writeControl    timeStep;', `        writeInterval   ${every};`, '    }'));
  out.push('    outletHoldup', '    {', '        type            surfaceFieldValue;', '        libs            (fieldFunctionObjects);', '        regionType      patch;', '        name            outlet;', '        operation       areaAverage;', '        fields          (alpha.liquid);', '        writeFields     false;', '        log             false;', '        writeControl    timeStep;', `        writeInterval   ${every};`, '    }',
    '    probes', '    {', '        type            probes;', '        libs            (sampling);', '        fields          (p p_rgh alpha.liquid);', '        probeLocations', '        (', stationsAt.map((p) => `            ${vec(add(p.c, add(mul(p.e1, 0.07 * mesh.D), mul(p.e2, 0.07 * mesh.D))))}`), '        );', '        writeControl    timeStep;', `        writeInterval   ${every};`, '    }');
  if (forces) out.push('    forces', '    {', '        type            forces;', '        libs            (forces);', '        patches         (wall);', '        rho             rho;', '        CofR            (0 0 0);', '        log             false;', '        writeControl    timeStep;', `        writeInterval   ${every};`, '    }');
  out.push('}');
  return out;
}
const COLLECT = text(['#!/bin/sh', '# Gathers the time series written by the function objects into one text file that HydraSlug imports.', 'cd "${0%/*}" || exit', 'out=hydraslug_results.txt', ': > "$out"', 'for f in $(find postProcessing -type f | sort)', 'do', "    printf '### file: %s\\n' \"$f\" >> \"$out\"", '    cat "$f" >> "$out"', "    printf '\\n' >> \"$out\"", 'done', 'echo "Wrote $out — import this file on the External solvers page."']);
const allrun = (steps, nProcs, extra = []) => text(['#!/bin/sh', 'cd "${0%/*}" || exit                                # run from this directory', '. ${WM_PROJECT_DIR:?}/bin/tools/RunFunctions        # OpenFOAM run functions', '#------------------------------------------------------------------------------', '', ...steps.map((s) => (s.startsWith('!') ? s.slice(1) : 'runApplication ' + s)), '',
  nProcs > 1 ? ['runApplication decomposePar', 'runParallel $(getApplication)', 'runApplication reconstructPar'] : 'runApplication $(getApplication)', ...extra, '', './collect.sh', '', '#------------------------------------------------------------------------------']);
const ALLCLEAN = text(['#!/bin/sh', 'cd "${0%/*}" || exit                                # run from this directory', '. ${WM_PROJECT_DIR:?}/bin/tools/CleanFunctions      # OpenFOAM clean functions', '#------------------------------------------------------------------------------', '', 'cleanCase0', 'rm -f hydraslug_results.txt', '', '#------------------------------------------------------------------------------']);
const decomposeDict = (n, title) => dictFile('decomposeParDict', 'system', title, ['numberOfSubdomains ' + n + ';', '', 'method          hierarchical;', '', 'coeffs', '{', `    n           (${n} 1 1);`, '}']);
const FLAVOUR_NOTE = ['## Other OpenFOAM versions', '', `The dictionaries are written for ${OPENFOAM_FLAVOUR}. With the OpenFOAM Foundation releases (openfoam.org, v10 and later) the same`, 'physics is set up with these differences, which you must apply by hand:', '', '- `constant/turbulenceProperties` is called `constant/momentumTransport`;', '- `constant/transportProperties` is split into `constant/physicalProperties.liquid`, `constant/physicalProperties.gas` and `constant/phaseProperties` (phases and sigma);', '- the solver is started as `foamRun -solver incompressibleVoF` (v11 and later) instead of `interFoam`;', '- in `surfaceFieldValue` / `volFieldValue` the key `regionType` is called `select`, and `name` is called `patch` / `cellZone`;', '- `createPatch` and `topoSet` read the same dictionaries.', ''];
const kv = (rows) => ['| Quantity | Value |', '|---|---|', ...rows.filter(Boolean).map(([a, b]) => `| ${a} | ${b} |`), ''];
const caseTable = (c) => kv([['Location along the line', `${ff(c.x, 6)} m (elevation ${ff(c.z, 5)} m), inclination ${ff(c.angleDeg, 4)}° upward in the flow direction`], ['Internal diameter', `${ff(c.D, 5)} m`], ['Pressure / temperature', `${ff(c.P, 5)} bara / ${ff(c.T, 4)} °C`], ['Superficial velocities', `liquid ${ff(c.vsl, 4)} m/s, gas ${ff(c.vsg, 4)} m/s`], ['Liquid holdup (1-D model)', ff(c.holdup, 4)], ['Flow regime (1-D model)', c.regime || '–'], ['Liquid density / viscosity', `${ff(c.rhoL, 5)} kg/m³ / ${ff(c.muL, 4)} Pa·s`], ['Gas density / viscosity', `${ff(c.rhoG, 5)} kg/m³ / ${ff(c.muG, 4)} Pa·s`], ['Interfacial tension', `${ff(c.sigma, 4)} N/m`], ['Wall shear stress / pressure gradient (1-D model)', `${ff(c.tauW, 4)} Pa / ${ff(c.dpdx, 4)} Pa/m`], ['Data source', c.source]]);
const hours = (h) => (h < 1 ? `${ff(h * 60, 2)} core-minutes` : h < 1e4 ? `${ff(h, 3)} core-hours` : `${h.toExponential(1)} core-hours`);

/**
 * Volume-of-fluid case on a swept mesh. cfg: { title, name, c (case data), mesh, gAngleDeg, solver, turbulence, isoAdvector,
 * levelZ (liquid below this height at the inlet and initially), qL, qG (m³/s) or uL, uG (m/s) or tables { t, qL, qG }, pOut (Pa, compressible only),
 * endTime, dt0, maxCo, writeInterval, nProcs, sectionsAt (fractions of the length), forces, intro[], planRows[] }.
 */
function vofCase(cfg) {
  const { c, mesh, title } = cfg, solver = cfg.solver || 'interFoam', comp = solver === 'compressibleInterFoam', empty = mesh.type === 'slot', th = (num(cfg.gAngleDeg, 0) * Math.PI) / 180;
  const gvec = [-G * Math.sin(th), 0, -G * Math.cos(th)], nuL = c.muL / c.rhoL, nuG = c.muG / c.rhoG, vm = c.vsl + c.vsg, nProcs = Math.max(1, Math.round(num(cfg.nProcs, 4)));
  const ts = turbulenceSetup(cfg.turbulence, { U: vm, D: c.D, nu: nuL, title, empty }), E = empty ? { frontAndBack: ['type            empty;'] } : {};
  const st = mesh.stations, sAt = (f) => { // frame at a fraction of the length
    const s = f * mesh.length; let i = 1; while (i < st.length - 1 && st[i].s < s) i++;
    const a = st[i - 1], b = st[i], w = (s - a.s) / (b.s - a.s || 1); return b.kind === 'bend' ? (w < 0.5 ? a : b) : { ...b, c: add(a.c, mul(a.t, s - a.s)), s };
  };
  const sections = (cfg.sectionsAt || [0.25, 0.5, 0.75]).map(sAt), every = Math.max(1, Math.round(num(cfg.sampleEvery, 5)));
  const big = 1e4, eps = Math.max(1e-6, 1e-4 * c.D), files = [];
  // inlet: the faces below the liquid level become the liquid inlet
  const tab = cfg.tables, inletU = (q, u, key) => (tab ? ['type            flowRateInletVelocity;', `volumetricFlowRate ${table1(tab.t, tab[key])};`, 'value           uniform (0 0 0);'] : isNum(q) ? ['type            flowRateInletVelocity;', `volumetricFlowRate ${ff(q, 7)};`, 'value           uniform (0 0 0);'] : ['type            fixedValue;', `value           uniform (${ff(u, 6)} 0 0);`]);
  const outU = ['type            inletOutlet;', 'inletValue      uniform (0 0 0);', 'value           uniform (0 0 0);'];
  files.push({ path: '0.orig/U', text: fieldFile('volVectorField', 'U', '[0 1 -1 0 0 0 0]', 'uniform (0 0 0)', { inletLiquid: inletU(cfg.qL, cfg.uL, 'qL'), inletGas: inletU(cfg.qG, cfg.uG, 'qG'), outlet: outU, wall: ['type            noSlip;'], ...E }, title) });
  files.push({ path: '0.orig/alpha.liquid', text: fieldFile('volScalarField', 'alpha.liquid', '[0 0 0 0 0 0 0]', 'uniform 0', { inletLiquid: ['type            fixedValue;', 'value           uniform 1;'], inletGas: ['type            fixedValue;', 'value           uniform 0;'], outlet: ['type            inletOutlet;', 'inletValue      uniform 0;', 'value           uniform 0;'], wall: ['type            zeroGradient;'], ...E }, title) });
  const p0 = comp ? pos(cfg.pOut, c.P * 1e5) : 0, ffp = ['type            fixedFluxPressure;', `value           uniform ${ff(p0)};`];
  files.push({ path: '0.orig/p_rgh', text: fieldFile('volScalarField', 'p_rgh', '[1 -1 -2 0 0 0 0]', `uniform ${ff(p0)}`, { [INLETS]: ffp, outlet: cfg.pOutTable ? ['type            uniformFixedValue;', `uniformValue    ${table1(cfg.pOutTable.t, cfg.pOutTable.p, 8)};`, `value           uniform ${ff(p0)};`] : ['type            fixedValue;', `value           uniform ${ff(p0)};`], wall: ffp, ...E }, title) });
  if (comp) {
    const calc = ['type            calculated;', `value           uniform ${ff(p0)};`], TK = c.T + KEL;
    files.push({ path: '0.orig/p', text: fieldFile('volScalarField', 'p', '[1 -1 -2 0 0 0 0]', `uniform ${ff(p0)}`, { [INLETS]: calc, outlet: calc, wall: calc, ...E }, title) });
    files.push({ path: '0.orig/T', text: fieldFile('volScalarField', 'T', '[0 0 0 1 0 0 0]', `uniform ${ff(TK, 6)}`, { [INLETS]: ['type            fixedValue;', `value           uniform ${ff(TK, 6)};`], outlet: ['type            inletOutlet;', `inletValue      uniform ${ff(TK, 6)};`, `value           uniform ${ff(TK, 6)};`], wall: ['type            zeroGradient;'], ...E }, title) });
    if (ts.tm.type !== 'laminar') files.push({ path: '0.orig/alphat', text: fieldFile('volScalarField', 'alphat', '[1 -1 -1 0 0 0 0]', 'uniform 0', { [INLETS]: ['type            calculated;', 'value           uniform 0;'], outlet: ['type            calculated;', 'value           uniform 0;'], wall: ['type            compressible::alphatWallFunction;', 'Prt             0.85;', 'value           uniform 0;'], ...E }, title) });
    const thermo = (eos, mw, cp, mu, Pr, extra) => ['thermoType', '{', '    type            heRhoThermo;', '    mixture         pureMixture;', '    transport       const;', '    thermo          hConst;', `    equationOfState ${eos};`, '    specie          specie;', '    energy          sensibleInternalEnergy;', '}', '', 'mixture', '{', '    specie', '    {', `        molWeight   ${ff(mw, 6)};`, '    }', ...extra, '    thermodynamics', '    {', `        Cp          ${ff(cp, 5)};`, '        Hf          0;', '    }', '    transport', '    {', `        mu          ${ff(mu, 5)};`, `        Pr          ${ff(Pr, 4)};`, '    }', '}'];
    // molar mass chosen so that the perfect-gas law reproduces the gas density of the property table at the operating point
    const mwEff = (c.rhoG * R * TK) / (c.P * 1e5) * 1000;
    files.push({ path: 'constant/thermophysicalProperties', text: dictFile('thermophysicalProperties', 'constant', title, ['phases          (liquid gas);', '', 'pMin            10000;', '', 'sigma', '{', '    type        constant;', `    sigma       ${ff(c.sigma, 5)};`, '}']) });
    files.push({ path: 'constant/thermophysicalProperties.gas', text: dictFile('thermophysicalProperties.gas', 'constant', title, ['// molWeight is an effective value: with it the perfect-gas law gives the gas density of the case', `// (${ff(c.rhoG, 5)} kg/m3 at ${ff(c.P, 5)} bara, ${ff(c.T, 4)} degC), i.e. it contains the compressibility factor.`, ...thermo('perfectGas', mwEff, pos(c.cpG, 2400), c.muG, (pos(c.cpG, 2400) * c.muG) / pos(c.kG, 0.04), [])]) });
    files.push({ path: 'constant/thermophysicalProperties.liquid', text: dictFile('thermophysicalProperties.liquid', 'constant', title, thermo('rhoConst', 100, pos(c.cpL, 2100), c.muL, (pos(c.cpL, 2100) * c.muL) / pos(c.kL, 0.13), ['    equationOfState', '    {', `        rho         ${ff(c.rhoL, 6)};`, '    }'])) });
  } else {
    files.push({ path: 'constant/transportProperties', text: dictFile('transportProperties', 'constant', title, ['phases          (liquid gas);', '', 'liquid', '{', '    transportModel  Newtonian;', `    nu              ${ff(nuL, 6)};`, `    rho             ${ff(c.rhoL, 6)};`, '}', '', 'gas', '{', '    transportModel  Newtonian;', `    nu              ${ff(nuG, 6)};`, `    rho             ${ff(c.rhoG, 6)};`, '}', '', `sigma           ${ff(c.sigma, 5)};`]) });
  }
  files.push(...ts.files);
  files.push({ path: 'constant/turbulenceProperties', text: dictFile('turbulenceProperties', 'constant', title, ts.props) });
  files.push({ path: 'constant/g', text: text([foamHeader('uniformDimensionedVectorField', 'g', title, 'constant'), 'dimensions      [0 1 -2 0 0 0 0];', `value           ${vec(gvec, 7)};`]) });
  files.push({ path: 'system/blockMeshDict', text: blockMeshDict(mesh, title) });
  files.push({ path: 'system/topoSetDict', text: dictFile('topoSetDict', 'system', title, ['// the inlet faces below the liquid level become the liquid inlet', 'actions', '(', '    {', '        name    inletLiquidFaces;', '        type    faceSet;', '        action  new;', '        source  patchToFace;', '        patch   inletGas;', '    }', '    {', '        name    inletLiquidFaces;', '        type    faceSet;', '        action  subset;', '        source  boxToFace;', `        box     (${ff(-eps)} ${-big} ${-big}) (${ff(eps)} ${big} ${ff(cfg.levelZ, 7)});`, '    }', ');']) });
  files.push({ path: 'system/createPatchDict', text: dictFile('createPatchDict', 'system', title, ['pointSync false;', '', 'patches', '(', '    {', '        name            inletLiquid;', '        patchInfo', '        {', '            type        patch;', '        }', '        constructFrom   set;', '        set             inletLiquidFaces;', '    }', ');']) });
  files.push({ path: 'system/setFieldsDict', text: dictFile('setFieldsDict', 'system', title, ['defaultFieldValues', '(', '    volScalarFieldValue alpha.liquid 0', ');', '', 'regions', '(', '    boxToCell', '    {', `        box (${-big} ${-big} ${-big}) (${big} ${big} ${ff(cfg.levelZ, 7)});`, '        fieldValues', '        (', '            volScalarFieldValue alpha.liquid 1', '        );', '    }', ');']) });
  const les = ts.tm.type !== 'RAS', conv = les ? 'Gauss limitedLinear 1' : 'Gauss upwind';
  files.push({ path: 'system/fvSchemes', text: dictFile('fvSchemes', 'system', title, ['ddtSchemes', '{', `    default         ${les && ts.tm.type !== 'laminar' && !cfg.isoAdvector ? 'CrankNicolson 0.9' : 'Euler'};`, '}', '', 'gradSchemes', '{', '    default         Gauss linear;', '}', '', 'divSchemes', '{', `    div(rhoPhi,U)   Gauss ${ts.tm.type === 'RAS' ? 'linearUpwind grad(U)' : ts.tm.type === 'laminar' ? 'linear' : 'LUST grad(U)'};`, '    div(phi,alpha)  Gauss vanLeer;', '    div(phirb,alpha) Gauss linear;',
    comp ? ['    div(rhoPhi,T)   Gauss upwind;', '    div(rhoPhi,K)   Gauss upwind;', '    div(phi,p)      Gauss upwind;'] : null,
    `    "div\\((phi|rhoPhi),(k|omega|epsilon|nuTilda|R)\\)" ${conv};`, '    div(R)          Gauss linear;', '    div((nuEff*dev2(T(grad(U))))) Gauss linear;', '    div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;', '    div(((rho*nu)*dev2(T(grad(U))))) Gauss linear;', '    div((rho*R))    Gauss linear;', '}', '', 'laplacianSchemes', '{', '    default         Gauss linear corrected;', '}', '', 'interpolationSchemes', '{', '    default         linear;', '}', '', 'snGradSchemes', '{', '    default         corrected;', '}', '', 'wallDist', '{', '    method          meshWave;', '}']) });
  const alphaCtl = cfg.isoAdvector ? ['        isoFaceTol      1e-8;', '        surfCellTol     1e-8;', '        nAlphaBounds    3;', '        snapTol         1e-12;', '        clip            true;', `        reconstructionScheme ${cfg.isoAdvector === 'plicRDF' ? 'plicRDF' : 'isoAlpha'};`, '        nAlphaSubCycles 1;', '        cAlpha          1;'] : ['        nAlphaCorr      2;', '        nAlphaSubCycles 1;', '        cAlpha          1;', '        MULESCorr       yes;', '        nLimiterIter    5;', '        solver          smoothSolver;', '        smoother        symGaussSeidel;', '        tolerance       1e-8;', '        relTol          0;'];
  files.push({ path: 'system/fvSolution', text: dictFile('fvSolution', 'system', title, ['solvers', '{', '    "alpha.liquid.*"', '    {', alphaCtl, '    }', '', comp ? ['    ".*(rho|rhoFinal)"', '    {', '        solver          diagonal;', '    }', ''] : null,
    '    "pcorr.*"', '    {', '        solver          PCG;', '        preconditioner  DIC;', '        tolerance       1e-5;', '        relTol          0;', '    }', '', '    p_rgh', '    {', '        solver          GAMG;', '        smoother        DIC;', '        tolerance       1e-7;', '        relTol          0.05;', '    }', '', '    p_rghFinal', '    {', '        solver          PCG;', '        preconditioner', '        {', '            preconditioner  GAMG;', '            smoother        DICGaussSeidel;', '            tolerance       1e-7;', '            relTol          0;', '            nVcycles        2;', '        }', '        tolerance       1e-8;', '        relTol          0;', '        maxIter         50;', '    }', '',
    '    "(U|T|k|epsilon|omega|nuTilda|R).*"', '    {', '        solver          smoothSolver;', '        smoother        symGaussSeidel;', '        tolerance       1e-7;', '        relTol          0;', '        minIter         1;', '    }', '}', '', 'PIMPLE', '{', '    momentumPredictor   no;', comp ? '    transonic           no;' : null, '    nOuterCorrectors    1;', '    nCorrectors         3;', `    nNonOrthogonalCorrectors ${mesh.stations.some((s) => s.kind === 'bend') || mesh.type === 'ogrid' ? 1 : 0};`, '}', '', 'relaxationFactors', '{', '    equations', '    {', '        ".*"            1;', '    }', '}']) });
  files.push({ path: 'system/decomposeParDict', text: decomposeDict(nProcs, title) });
  files.push({ path: 'system/controlDict', text: dictFile('controlDict', 'system', title, [`application     ${solver};`, '', 'startFrom       latestTime;', 'startTime       0;', 'stopAt          endTime;', `endTime         ${ff(cfg.endTime, 6)};`, `deltaT          ${ff(cfg.dt0, 3)};`, '', 'writeControl    adjustable;', `writeInterval   ${ff(cfg.writeInterval, 4)};`, 'purgeWrite      0;', 'writeFormat     binary;', 'writePrecision  8;', 'writeCompression off;', 'timeFormat      general;', 'timePrecision   8;', 'runTimeModifiable yes;', '', 'adjustTimeStep  on;', `maxCo           ${ff(cfg.maxCo, 3)};`, `maxAlphaCo      ${ff(Math.min(cfg.maxCo, cfg.isoAdvector ? 0.5 : 1), 3)};`, `maxDeltaT       ${ff(cfg.writeInterval, 3)};`, '', ...vofFunctions({ mesh, stationsAt: sections, every, forces: cfg.forces !== false })]) });
  files.push({ path: 'Allrun', text: allrun(['blockMesh', 'topoSet', 'createPatch -overwrite', '!restore0Dir', 'setFields', 'checkMesh'], nProcs) }, { path: 'Allclean', text: ALLCLEAN }, { path: 'collect.sh', text: COLLECT }, { path: 'case.foam', text: '' });
  return { files, sections, gvec, ts, nProcs };
}

/**
 * OpenFOAM VOF case of a pipe section at one location of the line. opts: { case (from bridgeCase), turbulence (key of TURBULENCE),
 * interface: 'vof' | 'isoAdvector' | 'plicRDF', geometry: 'straight' | 'bend', lengthD (section length in diameters), bendRadiusD, bendAngle (deg),
 * bendToward: 'up' | 'down' | 'left' | 'right', cellsPerDiameter, yPlus, flowThroughs, endTime, maxCo, nProcs, inletSeries: { t, vsl, vsg } }.
 * Returns { files, readme, summary, plan, mesh, name }.
 */
export function openfoamPipeCase(opts = {}) {
  const c = withCase(opts), tmKey = TURBULENCE[opts.turbulence] ? opts.turbulence : 'kOmegaSST', tm = TURBULENCE[tmKey], D = c.D, rdf = opts.interface === 'plicRDF', iso = opts.interface === 'isoAdvector' || rdf, bend = opts.geometry === 'bend';
  const lengthD = clamp(num(opts.lengthD, 20), 3, 2000), Rb = clamp(num(opts.bendRadiusD, 1.5), 1, 50) * D, ang = clamp(num(opts.bendAngle, 90), 5, 180);
  const legs = bend ? [{ type: 'straight', L: Math.max(2 * D, (lengthD * D - (Rb * ang * Math.PI) / 180) / 2) }, { type: 'bend', R: Rb, angleDeg: ang, toward: opts.bendToward || 'up' }, { type: 'straight', L: Math.max(2 * D, (lengthD * D - (Rb * ang * Math.PI) / 180) / 2) }] : [{ type: 'straight', L: lengthD * D }];
  const L = legs.reduce((s, l) => s + (l.type === 'bend' ? (l.R * l.angleDeg * Math.PI) / 180 : l.L), 0), plan = pipeMeshPlan(c, { ...opts, turbulence: tmKey, length: L });
  const mesh = ogridMesh({ D, legs, nc: plan.nc, nr: plan.nr, wallRatio: plan.wallRatio, dz: plan.dz }), A = mesh.area, hD = levelForHoldup(c.holdup), levelZ = -D / 2 + hD * D;
  const ser = opts.inletSeries && Array.isArray(opts.inletSeries.t) && opts.inletSeries.t.length > 1 ? opts.inletSeries : null, name = safeName(opts.name || `pipe_${tmKey}${rdf ? '_plicRDF' : iso ? '_iso' : ''}${bend ? '_bend' : ''}`);
  const title = `${tm.label}, ${rdf ? 'geometric VOF with reconstructed distance function (plicRDF)' : iso ? 'geometric VOF (isoAdvector)' : 'algebraic VOF (MULES)'}, ${bend ? `${ff(ang, 3)} deg bend` : 'straight section'} at x = ${ff(c.x, 6)} m`;
  const endTime = ser ? ser.t[ser.t.length - 1] - ser.t[0] : plan.endTime;
  const v = vofCase({ title, c, mesh, gAngleDeg: c.angleDeg, solver: iso ? 'interIsoFoam' : 'interFoam', turbulence: tmKey, isoAdvector: rdf ? 'plicRDF' : iso, levelZ, qL: c.vsl * A, qG: c.vsg * A,
    tables: ser ? { t: ser.t.map((t) => t - ser.t[0]), qL: ser.vsl.map((u) => u * A), qG: ser.vsg.map((u) => u * A) } : null,
    endTime, dt0: plan.dt / 10, maxCo: plan.maxCo, writeInterval: endTime / 50, nProcs: num(opts.nProcs, Math.min(64, Math.max(1, Math.round(mesh.nCells / 50000)))), forces: true });
  const steps = Math.ceil(endTime / plan.dt), coreHours = (mesh.nCells * steps * COST_PER_CELL_STEP) / 3600, slugUnitLen = c.slug.freq > 0 ? c.slug.velocity / c.slug.freq : 0;
  const readme = text([`# ${name} — OpenFOAM hand-off`, '', `Three-dimensional two-phase simulation of ${ff(L / D, 4)} diameters of pipe (${ff(L, 4)} m) at ${ff(c.x, 6)} m along the line: ${title}.`, '',
    'The one-dimensional models of the application give section-averaged holdup, pressure gradient and wall shear from closure laws. This case resolves the', 'velocity field, the gas–liquid interface and the turbulence in the section, so that those closures can be checked where it matters.', '',
    '## Case data (taken from the HydraSlug case)', '', ...caseTable(c),
    '## Mesh and cost', '', ...kv([['Block structure', `O-grid: square core + 4 wall blocks, ${mesh.blocks.length} blocks`], ['Cells per diameter (across) / per section', `${plan.nPerD} / ${plan.sectionCells}`], ['Radial cells in the wall blocks, growth rate', `${plan.nr}, ${ff(plan.growth, 4)}`], ['First cell height / y+ (thinner viscous layer of the two phases)', `${ff(plan.y1 * 1000, 3)} mm / ${ff(plan.yPlus, 3)}`], ['Axial cell length', `${ff(plan.dz * 1000, 4)} mm`], ['Total cells', mesh.nCells.toLocaleString('en-US')], ['Simulated time', `${ff(endTime, 4)} s${ser ? ' (length of the 1-D time series)' : ` (${ff(pos(opts.flowThroughs, 5), 3)} flow-through times)`}`], ['Time step at the Courant limit (estimate)', `${ff(plan.dt, 3)} s, about ${steps.toLocaleString('en-US')} steps`], ['Estimated cost', `${hours(coreHours)} (${ff(COST_PER_CELL_STEP * 1e6, 3)} core-µs per cell and step, measured with this set-up on one workstation core — time the first steps and rescale)`], ['Friction Reynolds number Re_τ', ff(plan.reTau, 4)], ['Cells for a full DNS of the same section', `${plan.dns.cells.toExponential(2)} (${plan.dns.coreHours.toExponential(1)} core-hours)`]]),
    plan.notes.map((n) => `- ${n}`), slugUnitLen > L ? `- The slug unit (slug + film) is about ${ff(slugUnitLen, 3)} m long, longer than this section: the case shows slug initiation and the film/front structure, not a train of developed slugs. Increase the length to several slug units to capture slug statistics.` : null, tm.type === 'LES' ? '- The inlet is a smooth two-layer flow without synthetic turbulence; allow two flow-through times before averaging, or map a fluctuating inlet.' : null, Math.abs(c.angleDeg) > 45 ? '- The line is steep here: the two-layer inlet is only a starting condition; the flow reorganises within a few diameters.' : null, '',
    '## Run', '', '```sh', 'source /path/to/OpenFOAM-v2412/etc/bashrc   # or: openfoam2412', `unzip ${name}.zip && cd ${name}`, './Allrun            # mesh, patches, initial fields, solver, result collection', '# ./Allclean       # removes everything that Allrun created', '```', '',
    `\`Allrun\` runs blockMesh → topoSet → createPatch (splits the inlet into a liquid and a gas inlet at the stratified level) → setFields → checkMesh → ${v.nProcs > 1 ? `decomposePar → mpirun -np ${v.nProcs} ` : ''}${iso ? 'interIsoFoam' : 'interFoam'}${v.nProcs > 1 ? ' -parallel → reconstructPar' : ''} → collect.sh.`, `Change \`numberOfSubdomains\` and \`n\` in system/decomposeParDict to use another number of cores.`, '',
    '## Results', '', 'The function objects in system/controlDict write, every ' + 5 + ' time steps:', '', '- `postProcessing/holdupVolume/0/volFieldValue.dat` — volume-averaged liquid holdup;', '- `postProcessing/section1..3/0/surfaceFieldValue.dat` — area-averaged holdup and pressure on three cross-sections (25 %, 50 %, 75 % of the length);', '- `postProcessing/outletHoldup/0/surfaceFieldValue.dat` — holdup at the outlet;', '- `postProcessing/probes/0/p`, `p_rgh`, `alpha.liquid` — point values near the centreline at the same sections;', '- `postProcessing/forces/0/force.dat` — pressure and viscous force of the fluid on the pipe wall (N).', '',
    '`collect.sh` gathers them in `hydraslug_results.txt`. On the External solvers page press “Import results” and choose that file (or the individual files): the application', 'computes mean holdup, slug frequency, pressure gradient, wall shear and peak force and sets them beside its own one-dimensional values.', '',
    ...FLAVOUR_NOTE]);
  return { name, files: [...v.files, { path: 'README.md', text: readme }], readme, plan: { ...plan, cells: mesh.nCells, steps, coreHours, endTime }, mesh, sections: v.sections,
    summary: `${iso ? 'interIsoFoam' : 'interFoam'} ${tm.label}: ${mesh.nCells.toLocaleString('en-US')} cells, y+ ≈ ${ff(plan.yPlus, 3)}, ${ff(endTime, 3)} s simulated, about ${hours(coreHours)}.`,
    commands: ['source /path/to/OpenFOAM-v2412/etc/bashrc', `unzip ${name}.zip && cd ${name}`, './Allrun'], reference: { holdup: c.holdup, dpdx: c.dpdx, tauW: c.tauW, slugFreq: c.slug.freq, forceN: bend ? c.slug.forceN : null, wallArea: Math.PI * D * L, sectionDistance: 0.5 * L, rhoM: c.holdup * c.rhoL + (1 - c.holdup) * c.rhoG } };
}

/**
 * Coupled 1-D / 3-D hand-off: the same three-dimensional section driven by the time histories of the one-dimensional
 * transient model (liquid and gas flow rates at the inlet). opts as openfoamPipeCase plus series: { t (s), vsl, vsg (m/s) };
 * without a series a slug-train history is synthesised from the slug frequency, length and body holdup of the case.
 */
export function openfoamCoupledCase(opts = {}) {
  const c = withCase(opts), s = opts.series || c.series;
  let ser = s && Array.isArray(s.t) && s.t.length > 3 && Array.isArray(s.vsl) && Array.isArray(s.vsg) ? { t: s.t.slice(), vsl: s.vsl.slice(), vsg: s.vsg.slice() } : null, synthetic = false;
  if (ser) { const mL = mean(ser.vsl), mG = mean(ser.vsg); if (!(mL > 0) || !(mG > 0)) ser = null; else { ser.vsl = ser.vsl.map((v) => (v * c.vsl) / mL); ser.vsg = ser.vsg.map((v) => (v * c.vsg) / mG); } } // the 1-D series is at the outlet: rescale its mean to the local superficial velocities
  if (!ser) { // slug train: liquid-rich slug bodies separated by gas-rich films, same time-averaged rates
    synthetic = true;
    const f = pos(c.slug.freq, 0.05), T = 1 / f, vm = c.vsl + c.vsg, tS = clamp(c.slug.length / Math.max(c.slug.velocity, 0.1), 0.05 * T, 0.6 * T), HS = clamp(c.slug.holdupBody, 0.5, 1), n = 5, ramp = 0.05 * T;
    const vslS = Math.min(HS * vm, (c.vsl * T) / tS), vslF = Math.max((c.vsl * T - vslS * tS) / (T - tS), 0.02 * c.vsl), t = [], vsl = [];
    for (let k = 0; k < n; k++) { t.push(k * T, k * T + ramp, k * T + tS, k * T + tS + ramp); vsl.push(vslF, vslS, vslS, vslF); }
    t.push(n * T); vsl.push(vslF);
    ser = { t, vsl, vsg: vsl.map((v) => Math.max(vm - v, 0.02 * c.vsg)) };
  }
  const out = openfoamPipeCase({ ...opts, case: c, inletSeries: ser, name: opts.name || `coupled_1d3d_${TURBULENCE[opts.turbulence] ? opts.turbulence : 'kOmegaSST'}` });
  const csv = ['t_s,vsl_m_s,vsg_m_s', ...ser.t.map((t, i) => `${ff(t - ser.t[0], 7)},${ff(ser.vsl[i], 6)},${ff(ser.vsg[i], 6)}`)].join('\n') + '\n';
  const extra = text(['', '## Coupling with the one-dimensional model', '', synthetic ? 'No transient 1-D result is stored in the case, so the inlet history is a slug train built from the slug frequency, slug length and slug-body holdup of the case' : 'The inlet history is the transient result of the one-dimensional flow model (liquid and gas rates), rescaled so that its mean equals the local superficial velocities', `(${ser.t.length} points over ${ff(ser.t[ser.t.length - 1] - ser.t[0], 4)} s; file \`constant/inletSeries.csv\`, written into 0.orig/U as flow-rate tables).`, '',
    'This is a one-way (1-D → 3-D) coupling per run. To close the loop, import the results: the application compares the three-dimensional holdup, pressure gradient and', 'wall shear with its closure values; adjust the calibration factors of the flow suite with those ratios, rerun its transient and write this case again. Two or three', 'such passes usually settle the exchange.']);
  const readme = out.readme + extra;
  return { ...out, files: [...out.files.filter((f) => f.path !== 'README.md'), { path: 'constant/inletSeries.csv', text: csv }, { path: 'README.md', text: readme }], readme, series: ser, synthetic, summary: 'Coupled 1-D → 3-D: ' + out.summary };
}

/**
 * Pipeline–riser case for severe slugging: a downward-inclined flowline, a riser-base bend and a vertical riser, as a planar 2-D
 * channel (default) or a coarse 3-D O-grid, with a compressible gas phase (compressibleInterFoam) because the cycle is driven by gas compression
 * in the flowline; `solver: 'interFoam'` writes the incompressible variant. opts: { case, dimension: '2d' | '3d', solver,
 * flowlineLengthD, riserHeightD, bendRadiusD, flowlineAngleDeg, cellsPerDiameter, turbulence, cycles, endTime, nProcs }.
 */
export function openfoamRiserSlugCase(opts = {}) {
  const c = withCase(opts), D = c.D, two = opts.dimension !== '3d', solver = opts.solver === 'interFoam' ? 'interFoam' : 'compressibleInterFoam', tmKey = TURBULENCE[opts.turbulence] ? opts.turbulence : 'kOmegaSST';
  const beta = clamp(num(opts.flowlineAngleDeg, Math.min(c.riser.flowlineAngleDeg, -0.5)), -30, -0.1), Lf = clamp(num(opts.flowlineLengthD, 120), 10, 5000) * D, H = clamp(num(opts.riserHeightD, 60), 5, 5000) * D, Rb = clamp(num(opts.bendRadiusD, 3), 1.5, 50) * D, nPerD = Math.round(clamp(num(opts.cellsPerDiameter, 16), 6, 200));
  const legs = [{ type: 'straight', L: Lf, parts: 2 }, { type: 'bend', R: Rb, angleDeg: 90 - beta, toward: 'up' }, { type: 'straight', L: H, parts: 2 }], dz = (2 * D) / nPerD;
  const nc = Math.max(3, Math.round(0.4 * nPerD)), nr = Math.max(2, Math.round((nPerD - nc) / 2)), mesh = two ? slotMesh({ D, legs, nD: nPerD, dz, width: D }) : ogridMesh({ D, legs, nc, nr, wallRatio: 0.5, dz });
  // the gas volume of the real flowline upstream sets the cycle period; the short model keeps the superficial velocities and the riser head
  const nLiq = clamp(Math.round(c.holdup * nPerD), 1, nPerD - 1), levelZ = two ? -D / 2 + (nLiq * D) / nPerD : -D / 2 + levelForHoldup(c.holdup) * D, Hin = two ? nLiq / nPerD : c.holdup;
  const pTop = Math.max(1.5e5, c.P * 1e5 - c.rhoL * G * H * 0.5), vm = c.vsl + c.vsg, period = Math.max((2 * H) / Math.max(c.vsl, 0.05), 10), endTime = pos(opts.endTime, pos(opts.cycles, 4) * period), dt = (0.5 * dz) / Math.max(4 * vm, 1);
  const name = safeName(opts.name || `riser_slugging_${two ? '2d' : '3d'}`), title = `Pipeline-riser severe slugging, ${two ? 'planar 2-D' : 'coarse 3-D'}, ${solver}`;
  const v = vofCase({ title, c: { ...c, P: pTop / 1e5 }, mesh, gAngleDeg: beta, solver, turbulence: tmKey, levelZ, uL: c.vsl / Hin, uG: c.vsg / (1 - Hin), qL: two ? null : c.vsl * mesh.area, qG: two ? null : c.vsg * mesh.area, pOut: pTop,
    endTime, dt0: dt / 10, maxCo: 0.5, writeInterval: endTime / 100, nProcs: num(opts.nProcs, Math.min(32, Math.max(1, Math.round(mesh.nCells / 40000)))), sectionsAt: [Lf / mesh.length * 0.5, (Lf + 0.5 * Rb * ((90 - beta) * Math.PI) / 180) / mesh.length, 0.97], forces: false });
  const steps = Math.ceil(endTime / dt), coreHours = (mesh.nCells * steps * COST_PER_CELL_STEP * (solver === 'interFoam' ? 1 : 1.6)) / 3600;
  const readme = text([`# ${name} — OpenFOAM hand-off`, '', `${title}: ${ff(Lf, 4)} m of flowline inclined ${ff(beta, 3)}°, a riser-base bend of radius ${ff(Rb, 3)} m and ${ff(H, 4)} m of vertical riser.`, '',
    'Severe (riser-induced) slugging is a cycle of liquid blocking the riser base, gas compressing in the flowline, blow-out and fall-back. The application screens it with', 'the Bøe and Pots criteria and a lumped cycle model; this case resolves the cycle itself. The gas is compressible (perfect gas with an effective molar mass that gives', 'the density of the case), the liquid is incompressible.', '',
    '## Case data', '', ...caseTable(c), ...kv([['Real riser height / modelled', `${ff(c.riser.height, 5)} m / ${ff(H, 4)} m`], ['Pressure at the riser top (outlet)', `${ff(pTop / 1e5, 4)} bara`], ['Inlet', `liquid in the lower ${ff(100 * Hin, 3)} % of the height at ${ff(c.vsl / Hin, 4)} m/s, gas above at ${ff(c.vsg / (1 - Hin), 4)} m/s`], ['Cells', `${mesh.nCells.toLocaleString('en-US')} (${nPerD} across the diameter)`], ['Simulated time', `${ff(endTime, 4)} s (about ${ff(endTime / period, 2)} riser fill times)`], ['Estimated cost', `${hours(coreHours)}`]]),
    '- The model is shorter than the real system. The cycle period scales with the compressible gas volume upstream of the riser base: lengthen the flowline (`flowlineLengthD`) or add', '  the real upstream gas volume to compare periods one-to-one; the Bøe criterion (liquid velocity below which the riser base blocks) is reproduced at any length.', two ? '- A planar channel has no pipe curvature: slug fronts are sharper and wall friction differs from a round pipe. Use the 3-D variant for quantitative amplitudes.' : null, '',
    '## Run', '', '```sh', 'source /path/to/OpenFOAM-v2412/etc/bashrc', `unzip ${name}.zip && cd ${name}`, './Allrun', '```', '', '## Results', '', '`collect.sh` writes `hydraslug_results.txt` with the holdup and pressure on three sections (mid-flowline, riser base, riser top) and at probes, plus the outlet holdup.', 'Import it on the External solvers page: the riser-base pressure history gives the cycle period and amplitude, compared with the severe-slugging result of the flow suite.', '', ...FLAVOUR_NOTE]);
  return { name, files: [...v.files, { path: 'README.md', text: readme }], readme, mesh, plan: { cells: mesh.nCells, steps, coreHours, endTime, dt }, sections: v.sections,
    summary: `${solver} ${two ? '2-D' : '3-D'} pipeline–riser: ${mesh.nCells.toLocaleString('en-US')} cells, ${ff(endTime, 3)} s simulated, about ${hours(coreHours)}.`, commands: ['source /path/to/OpenFOAM-v2412/etc/bashrc', `unzip ${name}.zip && cd ${name}`, './Allrun'], reference: { holdup: c.holdup, dpdx: null, tauW: null, slugFreq: 1 / period, forceN: null } };
}

// ---- CalculiX: shell and solid finite-element models of a pipe span, bend or jumper ----------------------------
const FE_GEOMETRY = Object.freeze({
  span: (o, D) => [{ type: 'straight', L: pos(o.spanLength, 15) }],
  bend: (o, D) => [{ type: 'straight', L: pos(o.legLength, 8 * D) }, { type: 'bend', R: clamp(num(o.bendRadiusD, 5), 1, 100) * D, angleDeg: clamp(num(o.bendAngle, 90), 5, 180), toward: o.bendToward || 'up' }, { type: 'straight', L: pos(o.legLength, 8 * D) }],
  // rigid M-shaped jumper: up, across, down (the first leg points along the path axis, which gravity is resolved against)
  jumper: (o, D) => [{ type: 'straight', L: pos(o.jumperHeight, 4) }, { type: 'bend', R: clamp(num(o.bendRadiusD, 3), 1, 100) * D, angleDeg: 90, toward: 'down' }, { type: 'straight', L: pos(o.jumperLength, 12) }, { type: 'bend', R: clamp(num(o.bendRadiusD, 3), 1, 100) * D, angleDeg: 90, toward: 'down' }, { type: 'straight', L: pos(o.jumperHeight, 4) }],
});
/**
 * Structured finite-element mesh of a pipe wall swept along `legs`.
 * spec: { D (inner diameter), wt, legs, element: 'shell' (S8R on the mid-surface) | 'solid' (C3D20R) | 'solid8' (C3D8I), nC (elements around), elemLen, nR (solid: elements through the wall) }.
 * Returns { element, type (CalculiX element type), nodes: [[x, y, z]] (index + 1 = node number), elements: [[node numbers]],
 *   sets: { END0, END1 (end rings), INNER, OUTER (solid: wall faces as element lists), WETTED (nodes of the inner wall), MID (ring nearest the middle), BEND (ring at the middle of the first bend) },
 *   stations, length, ri, ro, rm, nA, nC, nR, midFrame, bendFrame, forceDir }.
 */
export function pipeFeMesh(spec = {}) {
  const D = pos(spec.D, 0.254), wt = pos(spec.wt, 0.0159), ri = D / 2, ro = ri + wt, rm = ri + wt / 2, element = spec.element === 'solid' || spec.element === 'solid8' ? spec.element : 'shell', quad = element !== 'solid8';
  const nC = Math.max(8, 2 * Math.round(num(spec.nC, 16) / 2)), elemLen = pos(spec.elemLen, (Math.PI * (2 * rm)) / nC * 2), nR = element === 'shell' ? 0 : Math.max(1, Math.round(num(spec.nR, element === 'solid8' ? 2 : 1))), q = quad ? 2 : 1;
  const legsIn = spec.legs && spec.legs.length ? spec.legs : FE_GEOMETRY.span({}, D), legs = legsIn.map((l) => ({ ...l, div: q * Math.max(l.type === 'bend' ? 2 : 1, Math.round((l.type === 'bend' ? (l.R * Math.abs(l.angleDeg) * Math.PI) / 180 : l.L) / elemLen)) }));
  const st = sweepStations(legs), nS = st.length, nA = (nS - 1) / q, nJ = q * nC, nK = element === 'shell' ? 1 : q * nR + 1;
  const id = new Int32Array(nS * nJ * nK).fill(0), nodes = [], at = (i, j, k = 0) => id[(i * nJ + ((j % nJ) + nJ) % nJ) * nK + k];
  for (let i = 0; i < nS; i++) for (let j = 0; j < nJ; j++) for (let k = 0; k < nK; k++) {
    if (quad && (i % 2) + (j % 2) + (k % 2) > 1) continue; // quadratic serendipity elements: no face-centre or body-centre nodes
    const r = element === 'shell' ? rm : ri + (wt * k) / (nK - 1);
    nodes.push(sectionPoint(st[i], r, (2 * Math.PI * j) / nJ)); id[(i * nJ + j) * nK + k] = nodes.length;
  }
  const elements = [], inner = [], outer = [];
  for (let a = 0; a < nA; a++) for (let c = 0; c < nC; c++) {
    const i0 = q * a, j0 = q * c, i1 = i0 + q, j1 = j0 + q;
    if (element === 'shell') elements.push([at(i0, j0), at(i0, j1), at(i1, j1), at(i1, j0), at(i0, j0 + 1), at(i0 + 1, j1), at(i1, j0 + 1), at(i0 + 1, j0)]); // normal points outward
    else for (let b = 0; b < nR; b++) {
      const k0 = q * b, k1 = k0 + q, corner = [at(i0, j0, k0), at(i0, j0, k1), at(i0, j1, k1), at(i0, j1, k0), at(i1, j0, k0), at(i1, j0, k1), at(i1, j1, k1), at(i1, j1, k0)];
      elements.push(quad ? [...corner, at(i0, j0, k0 + 1), at(i0, j0 + 1, k1), at(i0, j1, k0 + 1), at(i0, j0 + 1, k0), at(i1, j0, k0 + 1), at(i1, j0 + 1, k1), at(i1, j1, k0 + 1), at(i1, j0 + 1, k0), at(i0 + 1, j0, k0), at(i0 + 1, j0, k1), at(i0 + 1, j1, k1), at(i0 + 1, j1, k0)] : corner);
      if (b === 0) inner.push(elements.length); if (b === nR - 1) outer.push(elements.length);
    }
  }
  const innerNodes = []; for (let i = 0; i < nS; i++) for (let j = 0; j < nJ; j++) { const n = at(i, j, 0); if (n) innerNodes.push(n); }
  const ring = (i) => { const out = []; for (let j = 0; j < nJ; j++) for (let k = 0; k < nK; k++) { const n = at(i, j, k); if (n) out.push(n); } return out; };
  const even = (i) => Math.min(nS - 1, Math.max(0, q * Math.round(i / q))), iMid = even((nS - 1) / 2);
  let iBend = -1; { let i = 1; while (i < nS && st[i].kind !== 'bend') i++; if (i < nS) { let e = i; while (e < nS - 1 && st[e + 1].kind === 'bend') e++; iBend = even((i - 1 + e) / 2); const tin = st[i - 1].t, tout = st[e].t, d = sub(tin, tout), L = hyp(d); spec._dir = L > 1e-9 ? mul(d, 1 / L) : null; spec._turn = L; } }
  return { element, type: element === 'shell' ? 'S8R' : element === 'solid' ? 'C3D20R' : 'C3D8I', nodes, elements, sets: { END0: ring(0), END1: ring(nS - 1), MID: ring(iMid), BEND: iBend >= 0 ? ring(iBend) : [], INNER: inner, OUTER: outer, WETTED: innerNodes },
    stations: st, length: st[nS - 1].s, ri, ro, rm, wt, nA, nC, nR, midFrame: st[iMid], bendFrame: iBend >= 0 ? st[iBend] : null, forceDir: spec._dir || null, turn: spec._turn || 0 };
}
/** Consistency checks of a finite-element mesh: node references, element orientation (positive corner Jacobians / outward shell normals), wall volume. */
export function checkFeMesh(m) {
  const issues = [], N = m.nodes.length, P = (n) => m.nodes[n - 1];
  let volume = 0, minJ = Infinity;
  const used = new Uint8Array(N + 1);
  m.elements.forEach((e, ei) => {
    if (e.some((n) => !(Number.isInteger(n) && n >= 1 && n <= N))) { issues.push(`element ${ei + 1} references a missing node`); return; }
    if (new Set(e).size !== e.length) issues.push(`element ${ei + 1} repeats a node`);
    e.forEach((n) => (used[n] = 1));
    if (m.element === 'shell') {
      const nrm = cross(sub(P(e[1]), P(e[0])), sub(P(e[3]), P(e[0]))), cen = mul(add(add(P(e[0]), P(e[1])), add(P(e[2]), P(e[3]))), 0.25);
      let best = null, bd = Infinity; for (const s of m.stations) { const d = hyp(sub(cen, s.c)); if (d < bd) { bd = d; best = s; } }
      const j = dot(nrm, sub(cen, best.c)); minJ = Math.min(minJ, j); if (!(j > 0)) issues.push(`shell element ${ei + 1} normal does not point outward`);
      // mid-surface area × thickness (two triangles per corner quad is too coarse on a curved wall: use the 8-node quadrature points instead)
      volume += shellArea(e.map(P)) * m.wt;
    } else {
      const c = e.slice(0, 8).map(P);
      const XYZ = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], at = (i, j, l) => c[XYZ.findIndex((q) => q[0] === i && q[1] === j && q[2] === l)];
      XYZ.forEach(([i, j, l]) => { const j3 = dot(cross(sub(at(1, j, l), at(0, j, l)), sub(at(i, 1, l), at(i, 0, l))), sub(at(i, j, 1), at(i, j, 0))); minJ = Math.min(minJ, j3); if (!(j3 > 0)) issues.push(`element ${ei + 1} has a non-positive corner Jacobian`); });
      volume += hexVolume(e.map(P));
    }
  });
  for (let n = 1; n <= N; n++) if (!used[n]) { issues.push(`node ${n} belongs to no element`); break; }
  for (const [name, set] of Object.entries(m.sets)) if (set.some((n) => !(n >= 1 && n <= (name === 'INNER' || name === 'OUTER' ? m.elements.length : N)))) issues.push(`set ${name} has an invalid entry`);
  const exact = Math.PI * (m.ro * m.ro - m.ri * m.ri) * m.length;
  return { ok: !issues.length, issues, volume, exactVolume: exact, minJacobian: minJ, nodes: N, elements: m.elements.length };
}
// 2 × 2 × 2 Gauss volume of a hexahedron (8-node trilinear, or 20-node serendipity when 20 points are given)
function hexVolume(p) {
  const g = 1 / Math.sqrt(3), C = [[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1], [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]], M = [[0, -1, -1], [1, 0, -1], [0, 1, -1], [-1, 0, -1], [0, -1, 1], [1, 0, 1], [0, 1, 1], [-1, 0, 1], [-1, -1, 0], [1, -1, 0], [1, 1, 0], [-1, 1, 0]];
  let V = 0;
  for (const x of [-g, g]) for (const y of [-g, g]) for (const z of [-g, g]) {
    const J = [[0, 0, 0], [0, 0, 0], [0, 0, 0]], addN = (pt, d) => { for (let a = 0; a < 3; a++) for (let b = 0; b < 3; b++) J[a][b] += d[a] * pt[b]; };
    if (p.length === 8) C.forEach(([a, b, c], n) => addN(p[n], [(a * (1 + b * y) * (1 + c * z)) / 8, (b * (1 + a * x) * (1 + c * z)) / 8, (c * (1 + a * x) * (1 + b * y)) / 8]));
    else {
      C.forEach(([a, b, c], n) => { const s = a * x + b * y + c * z - 2, fx = 1 + a * x, fy = 1 + b * y, fz = 1 + c * z; addN(p[n], [(a * fy * fz * (s + fx)) / 8, (b * fx * fz * (s + fy)) / 8, (c * fx * fy * (s + fz)) / 8]); });
      M.forEach(([a, b, c], n) => { let d; if (a === 0) d = [(-2 * x * (1 + b * y) * (1 + c * z)) / 4, (b * (1 - x * x) * (1 + c * z)) / 4, (c * (1 - x * x) * (1 + b * y)) / 4]; else if (b === 0) d = [(a * (1 - y * y) * (1 + c * z)) / 4, (-2 * y * (1 + a * x) * (1 + c * z)) / 4, (c * (1 + a * x) * (1 - y * y)) / 4]; else d = [(a * (1 + b * y) * (1 - z * z)) / 4, (b * (1 + a * x) * (1 - z * z)) / 4, (-2 * z * (1 + a * x) * (1 + b * y)) / 4]; addN(p[8 + n], d); });
    }
    V += dot(J[0], cross(J[1], J[2]));
  }
  return V;
}
// area of an 8-node serendipity quadrilateral (corners 0-3, mid-sides 4-7) by 2 × 2 Gauss quadrature
function shellArea(p) {
  const g = 1 / Math.sqrt(3), C = [[-1, -1], [1, -1], [1, 1], [-1, 1]], M = [[0, -1], [1, 0], [0, 1], [-1, 0]];
  let A = 0;
  for (const x of [-g, g]) for (const y of [-g, g]) {
    let dx = [0, 0, 0], dy = [0, 0, 0];
    C.forEach(([a, b], n) => { dx = add(dx, mul(p[n], (a * (1 + b * y) * (2 * a * x + b * y)) / 4)); dy = add(dy, mul(p[n], (b * (1 + a * x) * (a * x + 2 * b * y)) / 4)); });
    M.forEach(([a, b], n) => { if (a === 0) { dx = add(dx, mul(p[4 + n], -x * (1 + b * y))); dy = add(dy, mul(p[4 + n], (b * (1 - x * x)) / 2)); } else { dx = add(dx, mul(p[4 + n], (a * (1 - y * y)) / 2)); dy = add(dy, mul(p[4 + n], -y * (1 + a * x))); } });
    A += hyp(cross(dx, dy));
  }
  return A;
}
const inpList = (items, per = 8, fmtItem = String) => { const out = []; for (let i = 0; i < items.length; i += per) out.push(items.slice(i, i + per).map(fmtItem).join(', ') + (i + per < items.length ? ',' : '')); return out; };
const inpSet = (kind, name, items) => (items.length ? [`*${kind}, ${kind === 'NSET' ? 'NSET' : 'ELSET'}=${name}`, ...inpList(items, 12).map((l) => l.replace(/,$/, ','))] : []);
/** CalculiX mesh include file: nodes, elements and sets of a mesh from pipeFeMesh. */
export function calculixMeshText(m, title = 'pipe') {
  return text([`** HydraSlug hand-off: ${title}`, `** ${m.nodes.length} nodes, ${m.elements.length} ${m.type} elements; units m, kg, s, N, Pa, K`, '*NODE, NSET=Nall', m.nodes.map((p, i) => `${i + 1}, ${ff(p[0], 10)}, ${ff(p[1], 10)}, ${ff(p[2], 10)}`),
    `*ELEMENT, TYPE=${m.type}, ELSET=Eall`, m.elements.flatMap((e, i) => { const all = [i + 1, ...e], out = []; for (let k = 0; k < all.length; k += 10) out.push(all.slice(k, k + 10).join(', ') + (k + 10 < all.length ? ',' : '')); return out; }),
    inpSet('NSET', 'END0', m.sets.END0), inpSet('NSET', 'END1', m.sets.END1), inpSet('NSET', 'MID', m.sets.MID), inpSet('NSET', 'BEND', m.sets.BEND), inpSet('ELSET', 'INNER', m.sets.INNER), inpSet('ELSET', 'OUTER', m.sets.OUTER), inpSet('NSET', 'Ninterface', m.sets.WETTED)]);
}
/** First natural frequency (Hz) of a straight span as a beam: ends 'pinned' or 'fixed'. */
export const beamFrequency = (L, EI, massPerLength, ends = 'fixed') => (((ends === 'pinned' ? Math.PI : 4.730040745) ** 2) / (2 * Math.PI * L * L)) * Math.sqrt(EI / massPerLength);

/**
 * CalculiX model of a pipe span, bend or jumper with three decks sharing one mesh: static (pressure, temperature, weight),
 * modal (natural frequencies) and dynamic (slug force history as an amplitude table).
 * opts: { case, geometry: 'span' | 'bend' | 'jumper', element: 'shell' | 'solid', ends: 'fixed' | 'pinned', spanLength, legLength,
 *   bendRadiusD, bendAngle, jumperHeight, jumperLength, nC, nR, elemLen, modes, addedMass (bool: external added mass of sea water),
 *   damping (ratio of critical, default 0.01), forceN (peak slug force), maxIncrements (dynamic deck; default 4000 for solids, which use mode superposition, 400 for shells, which use direct integration), name }.
 * Returns { files, readme, summary, mesh, name, commands, reference }.
 */
export function calculixPipeCase(opts = {}) {
  const c = withCase(opts), D = c.D, wt = c.wt, geometry = FE_GEOMETRY[opts.geometry] ? opts.geometry : 'span', element = opts.element === 'solid' ? 'solid' : 'shell', ends = opts.ends === 'pinned' ? 'pinned' : 'fixed';
  const mesh = pipeFeMesh({ D, wt, legs: FE_GEOMETRY[geometry](opts, D), element, nC: opts.nC, nR: opts.nR, elemLen: opts.elemLen }), st = c.steel;
  const ri = mesh.ri, ro = mesh.ro, As = Math.PI * (ro * ro - ri * ri), Ai = Math.PI * ri * ri, Ao = Math.PI * ro * ro, I = (Math.PI / 4) * (ro ** 4 - ri ** 4), rhoM = c.holdup * c.rhoL + (1 - c.holdup) * c.rhoG, sub1 = c.z < 0, rhoSea = 1025;
  // equivalent densities of the wall material: weight (contents, buoyancy when submerged) and inertia (contents, added mass of the surrounding water)
  const rhoWeight = st.rho + (rhoM * Ai - (sub1 ? rhoSea * Ao : 0)) / As, rhoInertia = st.rho + (rhoM * Ai + (sub1 && opts.addedMass !== false ? rhoSea * Ao : 0)) / As, mLin = rhoInertia * As;
  const theta = geometry === 'jumper' ? Math.PI / 2 : (c.angleDeg * Math.PI) / 180, gdir = [-Math.sin(theta), 0, -Math.cos(theta)], pi = c.P * 1e5, pe = c.pExt * 1e5, name = safeName(opts.name || `fea_${geometry}_${element}`);
  const f1 = beamFrequency(geometry === 'span' ? mesh.length : geometry === 'jumper' ? pos(opts.jumperLength, 12) : mesh.length, st.E * I, mLin, ends), zeta = clamp(num(opts.damping, 0.01), 0, 0.2), modes = Math.round(clamp(num(opts.modes, 10), 1, 200));
  // slug load: on a bend the momentum change of the slug body, on a straight span the extra weight of the slug body
  const vt = pos(c.slug.velocity, 3), bend = mesh.sets.BEND.length > 0, HS = clamp(c.slug.holdupBody, 0.3, 1), rhoS = HS * c.rhoL + (1 - HS) * c.rhoG;
  const F = pos(opts.forceN, bend ? (c.slug.forceN * mesh.turn) / Math.SQRT2 : (rhoS - rhoM) * Ai * G * Math.min(c.slug.length, mesh.length)), dirF = bend ? mesh.forceDir : gdir, loadSet = bend ? 'BEND' : 'MID', nLoad = mesh.sets[loadSet].length;
  const ramp = Math.max((bend ? mesh.turn * D : D) / vt, 0.25 / Math.max(f1, 1e-3)), hold = Math.max(c.slug.length / vt, 2 * ramp), period = c.slug.freq > 0 ? 1 / c.slug.freq : 4 * (hold + ramp);
  const dt = Math.min(ramp / 8, 1 / (20 * Math.max(f1, 1e-3))), nPulse = 2, tEndFull = nPulse * period, tEnd = Math.min(tEndFull, Math.round(clamp(num(opts.maxIncrements, opts.element === 'solid' ? 4000 : 400), 20, 50000)) * dt), amp = [];
  for (let k = 0; k < nPulse; k++) { const t0 = k * Math.min(period, tEnd / nPulse); amp.push([t0, 0], [t0 + ramp, 1], [t0 + ramp + Math.min(hold, 0.3 * (tEnd / nPulse)), 1], [t0 + 2 * ramp + Math.min(hold, 0.3 * (tEnd / nPulse)), 0]); }
  amp.push([tEnd, 0]);
  // mode superposition for solid elements; direct implicit integration with Rayleigh damping (fitted at f1 and 3 f1) for shells, whose expanded nodes CalculiX cannot use in a modal dynamic step
  const direct = element === 'shell', w1 = 2 * Math.PI * f1;
  const material = (rho, damp) => ['*MATERIAL, NAME=STEEL', '*ELASTIC', `${ff(st.E, 7)}, ${ff(st.nu, 4)}`, '*DENSITY', ff(rho, 7), '*EXPANSION, ZERO=' + ff(c.tInstall, 5), ff(st.alphaT, 6), damp ? `*DAMPING, ALPHA=${ff(1.5 * zeta * w1, 5)}, BETA=${ff(zeta / (2 * w1), 5)}` : null,
    element === 'shell' ? ['*SHELL SECTION, ELSET=Eall, MATERIAL=STEEL', ff(wt, 6)] : ['*SOLID SECTION, ELSET=Eall, MATERIAL=STEEL']];
  // end conditions on the end rings. fixed: every degree of freedom; pinned: no transverse movement of the ring, zero mean axial movement, the section free to rotate.
  const e1f = mesh.stations[mesh.stations.length - 1], eq = (set) => { const n = mesh.sets[set]; return ['*EQUATION', String(n.length), ...inpList(n, 4, (k) => `${k}, 1, 1.`).map((l) => l.replace(/,$/, ','))]; };
  const constraints = ends === 'fixed' ? ['** fixed ends: all degrees of freedom of the end rings', '*BOUNDARY', element === 'shell' ? ['END0, 1, 6', 'END1, 1, 6'] : ['END0, 1, 3', 'END1, 1, 3']]
    : ['** pinned ends: local axis 1 along the pipe at each end; transverse displacements fixed, mean axial displacement zero', '*TRANSFORM, NSET=END1, TYPE=R', `${e1f.t.map((v) => ff(v, 8)).join(', ')}, ${e1f.e1.map((v) => ff(v, 8)).join(', ')}`, eq('END0'), eq('END1'), '*BOUNDARY', 'END0, 2, 3', 'END1, 2, 3'];
  const pressure = element === 'shell' ? ['** net pressure referred to the mid-surface, (p_i r_i - p_e r_o) / r_m; the element normals point outward and a positive value acts along the normal', `Eall, P, ${ff((pi * ri - pe * ro) / mesh.rm, 8)}`] : [`INNER, P6, ${ff(pi, 8)}`, `OUTER, P4, ${ff(pe, 8)}`];
  const head = (what) => [`** HydraSlug hand-off — CalculiX ${what}`, `** ${geometry}, ${mesh.type} elements, ends ${ends}; units m, kg, s, N, Pa, K`, '*INCLUDE, INPUT=mesh.inc'];
  const out = ['*NODE FILE', 'U', '*EL FILE', 'S', `*NODE PRINT, NSET=${loadSet}`, 'U', '*NODE PRINT, NSET=END0, TOTALS=ONLY', 'RF'];
  const staticDeck = text([head('static deck: internal and external pressure, temperature, submerged weight'), material(rhoWeight), constraints, '*INITIAL CONDITIONS, TYPE=TEMPERATURE', `Nall, ${ff(c.tInstall, 5)}`,
    '** step 1: pressure only', '*STEP', '*STATIC', '*DLOAD', pressure, out, '*END STEP',
    '** step 2: pressure + operating temperature (restrained thermal expansion) + weight', '*STEP', '*STATIC', '*TEMPERATURE', `Nall, ${ff(c.T, 5)}`, '*DLOAD', `Eall, GRAV, ${ff(G, 7)}, ${gdir.map((v) => ff(v, 7)).join(', ')}`, out, '*END STEP']);
  const modalDeck = text([head('modal deck: natural frequencies and mode shapes (contents and added mass in the density)'), material(rhoInertia), constraints, '*STEP', '*FREQUENCY', String(modes), '*NODE FILE', 'U', '*END STEP']);
  const nInc = Math.ceil(tEnd / dt), every = Math.max(1, Math.round(nInc / 100));
  const dynDeck = text([head(direct ? 'dynamic deck: slug force history by implicit direct integration' : 'dynamic deck: eigenmodes (step 1), then the slug force history by mode superposition (step 2)'), material(rhoInertia, direct), constraints, '*AMPLITUDE, NAME=SLUG', inpList(amp, 4, ([t, a]) => `${ff(t, 7)}, ${ff(a, 5)}`),
    direct ? ['*STEP, INC=' + (nInc + 10), '*DYNAMIC, DIRECT, ALPHA=-0.05', `${ff(dt, 5)}, ${ff(tEnd, 7)}`] : ['*STEP', '*FREQUENCY, STORAGE=YES', String(modes), '*END STEP', '*STEP, INC=' + (nInc + 10), '*MODAL DYNAMIC', `${ff(dt, 5)}, ${ff(tEnd, 7)}`, '*MODAL DAMPING', `1, ${modes}, ${ff(zeta, 4)}`], '*CLOAD, AMPLITUDE=SLUG', [1, 2, 3].filter((k) => Math.abs(dirF[k - 1]) > 1e-9).map((k) => `${loadSet}, ${k}, ${ff((F * dirF[k - 1]) / nLoad, 8)}`),
    '*NODE FILE, FREQUENCY=' + every, 'U', '*EL FILE, FREQUENCY=' + every, 'S', `*NODE PRINT, NSET=${loadSet}, FREQUENCY=1`, 'U', '*END STEP']);
  const hoop = (pi * ri - pe * ro) / wt, lame = (pi * (ri * ri + ro * ro) - 2 * pe * ro * ro) / (ro * ro - ri * ri), thermal = -st.E * st.alphaT * (c.T - c.tInstall);
  const run = text(['#!/bin/sh', '# Runs the three decks one after another (CalculiX ccx 2.20 or later).', 'cd "${0%/*}" || exit', 'for job in pipe_static pipe_modal pipe_dynamic', 'do', '    ccx -i "$job" > "$job.log" 2>&1 || echo "ccx failed for $job (see $job.log)"', 'done', 'echo "Import pipe_static.frd, pipe_modal.dat and pipe_dynamic.dat (or .frd) on the External solvers page."']);
  const readme = text([`# ${name} — CalculiX hand-off`, '', `${element === 'shell' ? 'Shell (S8R, eight-node reduced-integration shell on the mid-surface)' : 'Solid (C3D20R, twenty-node reduced-integration brick)'} finite-element model of a pipe ${geometry} at ${ff(c.x, 6)} m along the line: ${mesh.nodes.length.toLocaleString('en-US')} nodes, ${mesh.elements.length.toLocaleString('en-US')} elements.`, '',
    'The integrity suite of the application uses beam and thick-cylinder formulae. This model resolves what they average out: local bending at the supports, ovalisation and', 'stress concentration in bends, the full stress field through the wall, the mode shapes, and the response to the slug force in time.', '',
    '## Model', '', ...kv([['Pipe', `ID ${ff(D * 1000, 5)} mm, wall ${ff(wt * 1000, 4)} mm, developed length ${ff(mesh.length, 5)} m`], ['Steel', `${st.grade}: E = ${ff(st.E / 1e9, 4)} GPa, ν = ${ff(st.nu, 3)}, ρ = ${ff(st.rho, 5)} kg/m³, α = ${ff(st.alphaT, 4)} 1/K, SMYS = ${ff(st.smys / 1e6, 4)} MPa`], ['Internal / external pressure', `${ff(c.P, 5)} / ${ff(c.pExt, 5)} bara`], ['Temperature: installation → operation', `${ff(c.tInstall, 4)} → ${ff(c.T, 4)} °C`], ['Contents', `mixture density ${ff(rhoM, 5)} kg/m³ (holdup ${ff(c.holdup, 3)})`], ['Equivalent wall density: weight / inertia', `${ff(rhoWeight, 5)} / ${ff(rhoInertia, 5)} kg/m³ (${sub1 ? 'submerged: buoyancy and added mass of sea water included' : 'in air'})`], ['End conditions', ends === 'fixed' ? 'fixed: every degree of freedom of both end rings' : 'pinned: end rings held transversally, zero mean axial displacement, sections free to rotate'], ['Slug load', `${ff(F / 1000, 4)} kN peak on node set ${loadSet} (${nLoad} nodes), ramp ${ff(ramp, 3)} s, ${nPulse} pulses, time step ${ff(dt, 3)} s over ${ff(tEnd, 4)} s${tEnd < tEndFull ? ` (slug period ${ff(period, 4)} s shortened to keep the run at ${nInc} increments)` : ''}`], ['Damping', direct ? `${ff(100 * zeta, 3)} % of critical (Rayleigh, fitted at f₁ and 3 f₁; implicit direct integration)` : `${ff(100 * zeta, 3)} % of critical on every mode (mode superposition with ${modes} modes)`]]),
    '## Hand values to compare with', '', ...kv([['Hoop stress, thin wall (p_i r_i − p_e r_o)/t', `${ff(hoop / 1e6, 5)} MPa`], ['Hoop stress at the bore, Lamé', `${ff(lame / 1e6, 5)} MPa`], ['Axial stress from restrained thermal expansion −E α ΔT', `${ff(thermal / 1e6, 5)} MPa`], [`First natural frequency of a ${ends} beam of the ${geometry === 'jumper' ? 'jumper length' : 'developed length'}`, `${ff(f1, 4)} Hz`]]),
    '## Run', '', '```sh', `unzip ${name}.zip && cd ${name}`, './run.sh                 # = ccx -i pipe_static; ccx -i pipe_modal; ccx -i pipe_dynamic', 'cgx pipe_static.frd      # optional: view the results', '```', '',
    '## Results', '', '- `pipe_static.frd` — displacements and stresses of the two static steps (pressure; pressure + temperature + weight);', '- `pipe_modal.dat` — eigenvalue table with the natural frequencies; `pipe_modal.frd` — mode shapes;', `- \`pipe_dynamic.dat\` — displacement history of node set ${loadSet}; \`pipe_dynamic.frd\` — fields at about 100 instants.`, '',
    'On the External solvers page press “Import results” and choose these files: the application extracts the maximum von Mises stress, the maximum displacement, the natural', 'frequencies and the dynamic amplification, and sets them beside the hoop, thermal, von Mises utilisation and natural frequency of the integrity suite.', '',
    '## Code_Aster', '', 'The mesh file `mesh.inc` is plain Abaqus-style node and element lists: Code_Aster users can read it through Salome-Meca (Mesh → Import → Abaqus .inp) and apply the same', 'loads, which are all listed in the table above.', '']);
  return { name, mesh, files: [{ path: 'mesh.inc', text: calculixMeshText(mesh, `${geometry} ${mesh.type}`) }, { path: 'pipe_static.inp', text: staticDeck }, { path: 'pipe_modal.inp', text: modalDeck }, { path: 'pipe_dynamic.inp', text: dynDeck }, { path: 'run.sh', text: run }, { path: 'README.md', text: readme }], readme,
    summary: `CalculiX ${mesh.type} ${geometry}: ${mesh.nodes.length.toLocaleString('en-US')} nodes, ${mesh.elements.length.toLocaleString('en-US')} elements; static, modal (${modes} modes) and dynamic (${nInc} increments, ${direct ? 'direct integration' : 'mode superposition'}) decks.`,
    commands: [`unzip ${name}.zip && cd ${name}`, './run.sh'], plan: { cells: mesh.elements.length, coreHours: (mesh.nodes.length * 0.02 + nInc * (direct ? mesh.nodes.length * 4e-3 : 0.2)) / 3600, steps: nInc },
    reference: { hoop, lame, thermal, f1, forceN: F, smys: st.smys, loadSet, staticDisp: null } };
}

// ---- OpenFOAM: Lagrangian particles through a bend (DPMFoam) with the eroded-volume field ---------------------
const SECONDS_PER_YEAR = 31557600;
/**
 * OpenFOAM DPMFoam case of sand or hydrate particles carried through a pipe bend: Lagrangian parcels in a resolved carrier flow with
 * drag, gravity, turbulent dispersion and wall rebound, and the Finnie eroded-volume field on the wall. opts: { case (from bridgeCase),
 * particles: 'sand' | 'hydrate', particleUm, rateKgD (hydrate; sand comes from the case), carrier: 'mixture' | 'gas' | 'liquid',
 * bendRadiusD, bendAngle (deg), bendToward: 'up' | 'down' | 'left' | 'right', legD (straight legs, diameters each), cellsPerDiameter,
 * turbulence: 'kOmegaSST' | 'kEpsilon', endTime, flowThroughs, nProcs, parcelsPerSecond, name }.
 * Returns { name, files, readme, summary, mesh, plan, commands, reference }.
 */
export function openfoamErosionCase(opts = {}) {
  const c = withCase(opts), D = c.D, Rp = D / 2, hydrate = opts.particles === 'hydrate', carrier = opts.carrier === 'gas' || opts.carrier === 'liquid' ? opts.carrier : 'mixture', tmKey = opts.turbulence === 'kEpsilon' ? 'kEpsilon' : 'kOmegaSST', tm = TURBULENCE[tmKey];
  // carrier: one incompressible fluid
  const H = c.holdup, rhoC = carrier === 'gas' ? c.rhoG : carrier === 'liquid' ? c.rhoL : H * c.rhoL + (1 - H) * c.rhoG, muC = carrier === 'gas' ? c.muG : carrier === 'liquid' ? c.muL : H * c.muL + (1 - H) * c.muG;
  const U = Math.max(carrier === 'gas' ? c.vsg / (1 - H) : carrier === 'liquid' ? c.vsl / H : c.vsl + c.vsg, 1e-3), nuC = muC / rhoC, Re = (U * D) / nuC;
  const carrierText = carrier === 'gas' ? 'the gas alone at its in-situ velocity vsg / (1 − holdup)' : carrier === 'liquid' ? 'the liquid alone at its in-situ velocity vsl / holdup' : 'a homogeneous gas–liquid mixture (holdup-weighted density and viscosity) at the mixture velocity vsl + vsg';
  // particles
  const dP = clamp(hydrate ? pos(opts.particleUm, 500) : pos(c.sand.dUm, 150), 1, 20000) * 1e-6, rhoP = hydrate ? pos(opts.particleDensity, 920) : pos(c.sand.rho, 2650), rateKgD = hydrate ? pos(opts.rateKgD, 1000) : pos(c.sand.rateKgD, 50), rateKgS = rateKgD / 86400;
  const mP = (rhoP * Math.PI * dP ** 3) / 6, particlesPerSecond = rateKgS / mP, eWall = hydrate ? 0.5 : 0.9, muWall = hydrate ? 0.3 : 0.1, kind = hydrate ? 'hydrate' : 'sand';
  const tauP = (rhoP * dP * dP) / (18 * muC), stokes = (tauP * U) / D;
  // geometry and mesh
  const Rb = clamp(num(opts.bendRadiusD, 1.5), 1, 50) * D, ang = clamp(num(opts.bendAngle, 90), 5, 180), legL = clamp(num(opts.legD, 5), 1, 200) * D, toward = ['up', 'down', 'left', 'right'].includes(opts.bendToward) ? opts.bendToward : 'up';
  const legs = [{ type: 'straight', L: legL }, { type: 'bend', R: Rb, angleDeg: ang, toward }, { type: 'straight', L: legL }], arc = (Rb * ang * Math.PI) / 180;
  const nPerD = Math.round(clamp(num(opts.cellsPerDiameter, 20), 6, 200)), nc = Math.max(3, Math.round(0.4 * nPerD)), nr0 = Math.max(2, Math.round((nPerD - nc) / 2)), dz = (1.5 * D) / nPerD, notes = [];
  // wall-function sizing: first cell at y+ of about 100 when the uniform wall cell would be thicker, never thinner than y+ = 30
  const tauW = 0.5 * rhoC * U * U * 0.046 * Math.max(Re, 2100) ** -0.2, uTau = Math.sqrt(tauW / rhoC), lv = nuC / uTau, Lr = Rp * (1 - 0.55 * 0.9), nr = nr0;
  const solveG = (y) => { if (y * nr >= Lr * 0.999) return 1; let lo = 1, hi = 5; for (let i = 0; i < 80; i++) { const g = 0.5 * (lo + hi); if ((y * (g ** nr - 1)) / (g - 1) < Lr) lo = g; else hi = g; } return 0.5 * (lo + hi); };
  let y1 = Lr / nr, growth = 1;
  if (y1 / lv > 150) { y1 = 100 * lv; growth = solveG(y1); while (growth > 1.3) { y1 *= 1.05; growth = solveG(y1); } }
  const wallRatio = 1 / growth ** (nr - 1), yPlus = y1 / lv;
  if (yPlus < 30) notes.push(`The first wall cell is at y+ ≈ ${ff(yPlus, 3)}, below the logarithmic layer: the wall functions blend toward the viscous sub-layer (k–ω SST) — use fewer cells per diameter or accept the blended wall treatment.`);
  if (yPlus > 300) notes.push(`The first wall cell is at y+ ≈ ${ff(yPlus, 3)}, above the usual wall-function range (30–300): use more cells per diameter.`);
  if (Re < 4000) notes.push(`The carrier Reynolds number is ${ff(Re, 4)}: the flow is laminar or transitional and a RANS closure over-predicts mixing.`);
  const mesh = ogridMesh({ D, legs, nc, nr, wallRatio, dz }), L = mesh.length, st = mesh.stations;
  // typical wall face: a quarter of the circumference carries nc faces; axial size of the bend cells on the centreline and on the outer side
  const bendSt = st.filter((s) => s.kind === 'bend'), segArc = Rb * bendSt[0].dPsi, nzBend = Math.max(1, Math.round(segArc / dz)), dTheta = (Math.PI * D) / (4 * nc), axialMid = segArc / nzBend, axialOuter = (axialMid * (Rb + Rp)) / Rb;
  const wallFaceArea = dTheta * axialOuter, wallFaceAreaCentre = dTheta * axialMid;
  // time
  const core = (Math.SQRT2 * 0.55 * Rp) / nc, uMax = 1.5 * U, maxCo = 1, dt = (maxCo * Math.min(dz, core, axialMid)) / uMax, flowThroughs = pos(opts.flowThroughs, 6), endTime = pos(opts.endTime, (flowThroughs * L) / U), steps = Math.ceil(endTime / dt);
  const pps = Math.max(1, Math.round(Math.min(pos(opts.parcelsPerSecond, clamp(2e5 / endTime, 200, 2e6)), particlesPerSecond))), nParticle = particlesPerSecond / pps, parcelsInSystem = (pps * L) / U;
  const coreHours = ((mesh.nCells + 0.2 * parcelsInSystem) * steps * COST_PER_CELL_STEP) / 3600, nProcs = Math.max(1, Math.round(num(opts.nProcs, Math.min(64, Math.max(1, Math.round(mesh.nCells / 50000))))));
  if (stokes < 0.05) notes.push(`The particle Stokes number (response time × velocity / diameter) is ${ff(stokes, 3)}: the particles follow the carrier closely and few reach the wall by inertia — impacts come mainly from turbulent dispersion, which a wall-function mesh resolves poorly.`);
  if (carrier === 'mixture' && H > 0.05 && H < 0.95) notes.push('The carrier is a homogeneous pseudo-fluid. In slug or stratified flow the particles really travel in the liquid and meet the wall at the liquid velocity and through a liquid film; run the `liquid` and `gas` carriers as bounds.');
  const name = safeName(opts.name || `particles_${kind}_${carrier}_bend`), title = `${kind} parcels in a ${ff(ang, 3)} deg bend, DPMFoam, ${tm.label}, carrier: ${carrier}`, ph = 'carrier', Q = 'kinematicCloudQ';
  const th = (c.angleDeg * Math.PI) / 180, gvec = [-G * Math.sin(th), 0, -G * Math.cos(th)], addForces = rhoC / rhoP > 0.05;
  // Finnie: plastic flow stress (indentation pressure) of the wall, about three times a flow stress 10 % above the specified minimum yield strength
  const smys = pos(c.steel.smys, 448e6), pFlow = 3.3 * smys, psi = 2, K = 2, every = Math.max(1, Math.round(num(opts.sampleEvery, 5)));

  const files = [], ts = turbulenceSetup(tmKey, { U, D, nu: nuC, title, suffix: '.' + ph }), phiLine = `phi             phi.${ph};`;
  files.push({ path: `0.orig/U.${ph}`, text: fieldFile('volVectorField', 'U.' + ph, '[0 1 -1 0 0 0 0]', 'uniform (0 0 0)', { inlet: ['type            fixedValue;', `value           uniform (${ff(U, 6)} 0 0);`], outlet: ['type            inletOutlet;', phiLine, 'inletValue      uniform (0 0 0);', 'value           uniform (0 0 0);'], wall: ['type            noSlip;'] }, title) });
  const ffp = ['type            fixedFluxPressure;', phiLine, 'value           uniform 0;'];
  files.push({ path: '0.orig/p', text: fieldFile('volScalarField', 'p', '[0 2 -2 0 0 0 0]', 'uniform 0', { inlet: ffp, outlet: ['type            fixedValue;', 'value           uniform 0;'], wall: ffp }, title) });
  // the flux of the carrier is called phi.carrier: name it in the inletOutlet conditions of the turbulence fields
  files.push(...ts.files.map((f) => ({ path: f.path, text: f.text.replace(/( *)type            inletOutlet;/g, `$1type            inletOutlet;\n$1${phiLine}`) })));
  files.push({ path: 'constant/transportProperties', text: dictFile('transportProperties', 'constant', title, [`// carrier: ${carrierText.replace(/−/g, '-').replace(/–/g, '-')}`, `continuousPhase ${ph};`, '', `rho.${ph}     ${ff(rhoC, 6)};`, '', 'transportModel  Newtonian;', `nu              ${ff(nuC, 6)};`]) });
  files.push({ path: `constant/turbulenceProperties.${ph}`, text: dictFile('turbulenceProperties.' + ph, 'constant', title, ts.props) });
  files.push({ path: 'constant/g', text: text([foamHeader('uniformDimensionedVectorField', 'g', title, 'constant'), 'dimensions      [0 1 -2 0 0 0 0];', `value           ${vec(gvec, 7)};`]) });
  files.push({ path: 'constant/kinematicCloudProperties', text: dictFile('kinematicCloudProperties', 'constant', title, ['solution', '{', '    active          true;', '    coupled         true;', '    transient       yes;', '    cellValueSourceCorrection off;', '    maxCo           0.3;', '',
    '    interpolationSchemes', '    {', `        rho.${ph}     cell;`, `        U.${ph}       cellPoint;`, `        mu.${ph}      cell;`, addForces ? '        DUcDt           cell;' : null, '    }', '', '    integrationSchemes', '    {', '        U               Euler;', '    }', '', '    sourceTerms', '    {', '        schemes', '        {', '            U               semiImplicit 1;', '        }', '    }', '}', '',
    'constantProperties', '{', `    rho0            ${ff(rhoP, 6)};`, '    youngsModulus   1e9;', '    poissonsRatio   0.3;', '    alphaMax        0.99;', '}', '',
    'subModels', '{', '    particleForces', '    {', '        sphereDrag;', '        gravity;', addForces ? ['        pressureGradient', '        {', `            U               U.${ph};`, '        }', '        virtualMass', '        {', '            Cvm             0.5;', '        }'] : null, '    }', '',
    '    injectionModels', '    {', '        inlet', '        {', '            type            patchInjection;', '            patch           inlet;', `            // ${ff(rateKgD, 6)} kg/d of ${kind}: ${ff(particlesPerSecond, 5)} particles/s, ${ff(nParticle, 5)} particles per parcel`, '            parcelBasisType mass;', `            massTotal       ${ff(rateKgS * endTime, 7)}; // kg injected over the duration`, '            SOI             0;', `            duration        ${ff(endTime, 6)};`, `            parcelsPerSecond ${pps};`, '            flowRateProfile constant 1;', `            U0              (${ff(U, 6)} 0 0);`, '            minParticlesPerParcel 0;',
    '            sizeDistribution', '            {', '                type            fixedValue;', '                fixedValueDistribution', '                {', `                    value           ${ff(dP, 6)};`, '                }', '            }', '        }', '    }', '',
    `    dispersionModel ${tm.type === 'RAS' ? 'stochasticDispersionRAS' : 'none'};`, '', '    patchInteractionModel localInteraction;', '', '    localInteractionCoeffs', '    {', '        patches', '        (', '            wall', '            {', '                type            rebound;', `                e               ${ff(eWall, 3)}; // normal restitution`, `                mu              ${ff(muWall, 3)}; // tangential loss`, '            }', '            "(inlet|outlet)"', '            {', '                type            escape;', '            }', '        );', '    }', '',
    '    surfaceFilmModel none;', '    collisionModel  none;', '    stochasticCollisionModel none;', '}', '',
    'cloudFunctions', '{', '    particleErosion1', '    {', '        type            particleErosion;', '        patches         (wall);', `        p               ${ff(pFlow, 5)}; // plastic flow stress of the wall, Pa (3.3 x SMYS of ${String(c.steel.grade).replace(/[^\w .-]/g, '')})`, `        psi             ${psi}; // depth of contact / depth of cut`, `        K               ${K}; // normal / tangential force on the particle`, '    }', '}']) });
  files.push({ path: 'system/blockMeshDict', text: blockMeshDict(mesh, title, 'inlet') });
  const conv = 'Gauss limitedLinear 1';
  files.push({ path: 'system/fvSchemes', text: dictFile('fvSchemes', 'system', title, ['ddtSchemes', '{', '    default         Euler;', '}', '', 'gradSchemes', '{', '    default         Gauss linear;', '}', '', 'divSchemes', '{', '    default         none;', `    div(alphaPhi.${ph},U.${ph}) Gauss linearUpwindV grad(U.${ph});`, `    "div\\(alphaPhi.${ph},(k|omega|epsilon).${ph}\\)" ${conv};`, `    div(((alpha.${ph}*nuEff.${ph})*dev2(T(grad(U.${ph}))))) Gauss linear;`, '}', '', 'laplacianSchemes', '{', '    default         Gauss linear corrected;', '}', '', 'interpolationSchemes', '{', '    default         linear;', '}', '', 'snGradSchemes', '{', '    default         corrected;', '}', '', 'wallDist', '{', '    method          meshWave;', '}']) });
  files.push({ path: 'system/fvSolution', text: dictFile('fvSolution', 'system', title, ['solvers', '{', '    p', '    {', '        solver          GAMG;', '        smoother        GaussSeidel;', '        tolerance       1e-7;', '        relTol          0.01;', '    }', '', '    pFinal', '    {', '        solver          GAMG;', '        smoother        GaussSeidel;', '        tolerance       1e-7;', '        relTol          0;', '    }', '',
    `    "(U|k|omega|epsilon).${ph}"`, '    {', '        solver          smoothSolver;', '        smoother        symGaussSeidel;', '        tolerance       1e-7;', '        relTol          0.1;', '    }', '', `    "(U|k|omega|epsilon).${ph}Final"`, '    {', '        solver          smoothSolver;', '        smoother        symGaussSeidel;', '        tolerance       1e-7;', '        relTol          0;', '    }', '}', '', 'PIMPLE', '{', '    nOuterCorrectors    1;', '    nCorrectors         2;', '    momentumPredictor   yes;', '    nNonOrthogonalCorrectors 1;', '}', '', 'relaxationFactors', '{', '    equations', '    {', '        ".*"            1;', '    }', '}']) });
  files.push({ path: 'system/decomposeParDict', text: decomposeDict(nProcs, title) });
  const sfv = (fo, patch, op, field) => [`    ${fo}`, '    {', '        type            surfaceFieldValue;', '        libs            (fieldFunctionObjects);', '        regionType      patch;', `        name            ${patch};`, `        operation       ${op};`, `        fields          (${field});`, '        writeFields     false;', '        log             false;', '        writeControl    timeStep;', `        writeInterval   ${every};`, '    }'];
  files.push({ path: 'system/controlDict', text: dictFile('controlDict', 'system', title, ['application     DPMFoam;', '', 'startFrom       latestTime;', 'startTime       0;', 'stopAt          endTime;', `endTime         ${ff(endTime, 6)};`, `deltaT          ${ff(dt / 5, 3)};`, '', 'writeControl    adjustable;', `writeInterval   ${ff(endTime / 20, 4)};`, 'purgeWrite      0;', 'writeFormat     binary;', 'writePrecision  8;', 'writeCompression off;', 'timeFormat      general;', 'timePrecision   8;', 'runTimeModifiable yes;', '', 'adjustTimeStep  on;', `maxCo           ${ff(maxCo, 3)};`, `maxDeltaT       ${ff(Math.min(endTime / 20, 5 * dt), 3)};`, '',
    'functions', '{',
    `    // eroded volume of wall material per wall face (m3), accumulated by the cloud function particleErosion in the field ${Q}`,
    ...sfv('erosionTotal', 'wall', 'sum', Q), ...sfv('erosionMax', 'wall', 'max', Q),
    `    // static pressure in Pa from the kinematic pressure of the solver (p x rho.${ph})`,
    '    staticPressure', '    {', '        type            pressure;', '        libs            (fieldFunctionObjects);', '        field           p;', '        mode            static;', '        rho             rhoInf;', `        rhoInf          ${ff(rhoC, 6)};`, '        result          pStatic;', '        log             false;', '        executeControl  timeStep;', `        executeInterval ${every};`, '        writeControl    writeTime;', '    }',
    ...sfv('inletPressure', 'inlet', 'areaAverage', 'pStatic'), ...sfv('outletPressure', 'outlet', 'areaAverage', 'pStatic'),
    '    cloudInfo', '    {', '        type            cloudInfo;', '        libs            (lagrangianFunctionObjects);', '        clouds          (kinematicCloud);', '        log             false;', '        writeControl    timeStep;', `        writeInterval   ${every};`, '    }', '}']) });
  files.push({ path: 'Allrun', text: allrun(['blockMesh', '!restore0Dir', 'checkMesh'], nProcs) }, { path: 'Allclean', text: ALLCLEAN }, { path: 'collect.sh', text: COLLECT }, { path: 'case.foam', text: '' });

  const readme = text([`# ${name} — OpenFOAM hand-off`, '', `Three-dimensional simulation of ${kind} particles carried through a ${ff(ang, 3)}° bend of radius ${ff(Rb / D, 3)} D with ${ff(legL / D, 3)} D of straight pipe on each side, at ${ff(c.x, 6)} m along the line.`, `Solver: DPMFoam (Lagrangian parcels coupled to a resolved incompressible carrier flow), ${tm.label}.`, '',
    'The erosion models of the application are correlations: one impact velocity (the mixture velocity), one impact angle and a geometry factor per component give one', 'wall-loss rate. This case tracks the particles themselves through the resolved flow — drag, gravity, turbulent dispersion and wall rebound — so it shows where on the', 'bend the particles hit, at which angle and speed, how the secondary flow and repeated impacts spread the wear, and it integrates the Finnie cutting-wear model', 'impact by impact into a map of removed wall volume.', '',
    '## Case data (taken from the HydraSlug case)', '', ...caseTable(c),
    '## Carrier and particles', '', ...kv([['Carrier fluid', carrierText], ['Carrier density / viscosity / velocity', `${ff(rhoC, 5)} kg/m³ / ${ff(muC, 4)} Pa·s / ${ff(U, 4)} m/s (Reynolds number ${ff(Re, 4)})`], ['Particles', `${kind}, diameter ${ff(dP * 1e6, 4)} µm, density ${ff(rhoP, 5)} kg/m³`], ['Particle rate', `${ff(rateKgD, 5)} kg/d = ${ff(rateKgS, 4)} kg/s = ${ff(particlesPerSecond, 4)} particles/s`], ['Parcels', `${pps.toLocaleString('en-US')} per second, ${ff(nParticle, 4)} particles each, about ${Math.round(parcelsInSystem).toLocaleString('en-US')} in the section at any time`], ['Particle response time / Stokes number', `${ff(tauP, 3)} s / ${ff(stokes, 3)}`], ['Forces on the particles', `sphere drag, gravity with buoyancy${addForces ? ', pressure gradient, virtual mass' : ''}; stochastic turbulent dispersion; no particle–particle collisions (dilute)`], ['Wall', `rebound with normal restitution ${ff(eWall, 3)} and tangential loss ${ff(muWall, 3)}; particles leave through the outlet`], ['Erosion model', `Finnie cutting wear (cloud function particleErosion): plastic flow stress p = ${ff(pFlow / 1e9, 4)} GPa (3.3 × SMYS of ${c.steel.grade}), ψ = ${psi}, K = ${K}`], ['Gravity in the mesh frame', `${vec(gvec, 5)} m/s² (x is the flow direction at the inlet, inclined ${ff(c.angleDeg, 4)}°)`]]),
    '## Mesh and cost', '', ...kv([['Block structure', `O-grid: square core + 4 wall blocks, ${mesh.blocks.length} blocks`], ['Cells per diameter (across) / per section', `${nPerD} / ${mesh.sectionCells}`], ['Radial cells in the wall blocks, growth rate', `${nr}, ${ff(growth, 4)}`], ['First cell height / y+', `${ff(y1 * 1000, 3)} mm / ${ff(yPlus, 3)}`], ['Wall face in the bend (circumferential × axial)', `${ff(dTheta * 1000, 4)} mm × ${ff(axialOuter * 1000, 4)} mm on the outer side = ${ff(wallFaceArea, 5)} m² (${ff(wallFaceAreaCentre, 5)} m² at the centreline radius)`], ['Total cells', mesh.nCells.toLocaleString('en-US')], ['Simulated time', `${ff(endTime, 4)} s (${ff((endTime * U) / L, 3)} flow-through times of ${ff(L / U, 3)} s)`], ['Time step at the Courant limit (estimate)', `${ff(dt, 3)} s, about ${steps.toLocaleString('en-US')} steps`], ['Estimated cost', `${hours(coreHours)} (cells × steps × ${ff(COST_PER_CELL_STEP * 1e6, 3)} core-µs, a parcel in the section counted as 0.2 cells — time the first steps and rescale)`]]),
    notes.map((n) => `- ${n}`), '- The flow and the particle field start from rest: the first flow-through time is a start-up transient. Read the erosion rate from the slope of the later part of the time series rather than from the total.', '- The erosion statistics on a single wall face converge slowly: the maximum is noisy until every face of the worn area has received many parcels. Raise `parcelsPerSecond` or the simulated time until the maximum grows linearly with time.', '',
    '## Run', '', '```sh', 'source /path/to/OpenFOAM-v2412/etc/bashrc   # or: openfoam2412', `unzip ${name}.zip && cd ${name}`, './Allrun            # mesh, initial fields, solver, result collection', '# ./Allclean       # removes everything that Allrun created', '```', '',
    `\`Allrun\` runs blockMesh → restore0Dir → checkMesh → ${nProcs > 1 ? `decomposePar → mpirun -np ${nProcs} ` : ''}DPMFoam${nProcs > 1 ? ' -parallel → reconstructPar' : ''} → collect.sh.`, 'Change `numberOfSubdomains` and `n` in system/decomposeParDict to use another number of cores.', '',
    '## Results', '', `The function objects in system/controlDict write, every ${every} time steps:`, '',
    `- \`postProcessing/erosionTotal/0/surfaceFieldValue.dat\` — \`sum(${Q})\`: volume of wall material removed on the whole wall since the start (m³);`, `- \`postProcessing/erosionMax/0/surfaceFieldValue.dat\` — \`max(${Q})\`: the largest volume removed on one wall face since the start (m³);`, '- `postProcessing/inletPressure/0/surfaceFieldValue.dat`, `postProcessing/outletPressure/0/surfaceFieldValue.dat` — `areaAverage(pStatic)`: area-averaged static pressure (Pa, relative to the outlet); their difference is the pressure drop of the section, static head included;', '- `postProcessing/cloudInfo/0/kinematicCloud.dat` — number of parcels and particle mass (kg) in the section, and the diameter statistics of the cloud.', '',
    `The field \`${Q}\` itself (eroded volume per wall face, m³) and the parcels (\`lagrangian/kinematicCloud\`) are stored in every time directory; open \`case.foam\` in ParaView to see the wear map on the bend.`, '',
    '### Erosion rate', '', `The wall-loss rate is the removed volume on the worst wall face divided by the area of that face and by the time over which it was removed:`, '', '    rate [m/s] = max(Q) [m³] / A_face [m²] / t [s];   rate [mm/year] = rate [m/s] × 1000 × 31 557 600', '',
    `with A_face = ${ff(wallFaceArea, 5)} m² (wall face on the outer side of the bend, where the maximum normally lies: ${ff(dTheta * 1000, 4)} mm around the circumference × ${ff(axialOuter * 1000, 4)} mm along the pipe) and t the simulated time`, `(${ff(endTime, 4)} s for the complete run; better: the slope of max(Q) against time after the start-up). The particle rate is the real one (${ff(rateKgS, 4)} kg/s), so no scaling to the field rate is needed.`, `For example max(Q) = 1e-15 m³ after ${ff(endTime, 4)} s would be ${ff((1e-15 / wallFaceArea / endTime) * 1000 * SECONDS_PER_YEAR, 3)} mm/year. The specific wear is sum(Q) × steel density (${ff(c.steel.rho, 5)} kg/m³) / mass of ${kind} passed.`, hydrate ? 'Hydrate particles are far softer than steel: for them the field is a map of impact intensity (where and how hard the particles hit — the places where deposits start), not a metal-loss prediction.' : null, '',
    '`collect.sh` gathers the files in `hydraslug_results.txt`. On the External solvers page press “Import results” and choose that file: the application converts the', 'erosion series into a wall-loss rate and sets it, with the pressure drop, beside its own correlation values.', '',
    ...FLAVOUR_NOTE, `For this case: the Foundation equivalent of DPMFoam is \`denseParticleFoam\` (or \`particleFoam\` for one-way coupled particles in a frozen flow; \`foamRun -solver incompressibleDenseParticleFluid\` from v11); the cloud dictionary is \`constant/cloudProperties\` there and the eroded-volume field is called \`Q\`.`, '']);
  return { name, files: [...files, { path: 'README.md', text: readme }], readme, mesh,
    summary: `DPMFoam ${kind} particles in a ${ff(ang, 3)}° bend (${tm.label}, carrier ${carrier}): ${mesh.nCells.toLocaleString('en-US')} cells, ${pps.toLocaleString('en-US')} parcels/s, ${ff(endTime, 3)} s simulated, about ${hours(coreHours)}.`,
    plan: { cells: mesh.nCells, steps, coreHours, endTime, dt, nc, nr, nPerD, wallRatio, growth, y1, yPlus, dz, parcelsPerSecond: pps, particlesPerParcel: nParticle, parcelsInSystem, notes },
    commands: ['source /path/to/OpenFOAM-v2412/etc/bashrc', `unzip ${name}.zip && cd ${name}`, './Allrun'],
    reference: { erosionRateMmY: null, wallFaceArea, wallFaceAreaCentre, simulatedTime: endTime, sandRateKgS: rateKgS, steelDensity: c.steel.rho, carrierDensity: rhoC, carrierVelocity: U, particleDiameter: dP, particleDensity: rhoP, erosionField: Q, flowStress: pFlow } };
}

// ---- OpenFOAM: Euler–Euler models (two-fluid, multi-fluid, size-class population balance) ------------------------
// Cost of one cell and time step relative to COST_PER_CELL_STEP (interFoam), measured on short runs of each variant (5,040 cells, one core,
// three outer correctors): 86, 84, 72 and 88 core-µs for gas–liquid, liquid–solids, three fluids and the population balance with 12 size classes.
const eulerCost = Object.freeze({ gas: 1.7, solids: 1.7, three: 1.4, population: (n) => 1.6 + 0.015 * n });
const eulerFv = (v) => ['type            fixedValue;', `value           uniform ${v};`];
const eulerIo = (v, ph) => ['type            inletOutlet;', `phi             phi.${ph};`, `inletValue      uniform ${v};`, `value           uniform ${v};`];
const eulerCalc = (v) => ['type            calculated;', `value           uniform ${v};`];
const eulerWf = (type, v, extra = []) => [`type            ${type};`, ...extra, `value           uniform ${v};`];
const eulerZg = ['type            zeroGradient;'];
const eulerDim = { one: '[0 0 0 0 0 0 0]', U: '[0 1 -1 0 0 0 0]', p: '[1 -1 -2 0 0 0 0]', T: '[0 0 0 1 0 0 0]', k: '[0 2 -2 0 0 0 0]', eps: '[0 2 -3 0 0 0 0]', nu: '[0 2 -1 0 0 0 0]', alphat: '[1 -1 -1 0 0 0 0]' };
const eulerFld = (cls, name, dims, internal, inlet, outlet, wall, title) => ({ path: '0.orig/' + name, text: fieldFile(cls, name, dims, 'uniform ' + internal, { inlet, outlet, wall }, title) });
const eulerPair = (pair, lines) => [`    (${pair})`, '    {', ...lines.map((l) => '        ' + l), '    }'];
const eulerList = (name, pairs) => [name, '(', ...pairs.filter(Boolean).flat(), ');', ''];
const eulerDrag = (type) => [`type            ${type};`, 'residualRe      1e-3;', 'swarmCorrection', '{', '    type        none;', '}'];
const eulerRas = (model, coeffs = []) => (model === 'laminar' ? ['simulationType  laminar;'] : ['simulationType  RAS;', '', 'RAS', '{', `    RASModel        ${model};`, '    turbulence      on;', '    printCoeffs     on;', ...coeffs.map((l) => '    ' + l), '}']);
const eulerThermo = (eos, mw, cp, mu, Pr, rho) => ['thermoType', '{', '    type            heRhoThermo;', '    mixture         pureMixture;', '    transport       const;', '    thermo          hConst;', `    equationOfState ${eos};`, '    specie          specie;', '    energy          sensibleInternalEnergy;', '}', '', 'mixture', '{', '    specie', '    {', `        molWeight   ${ff(mw, 6)};`, '    }', ...(eos === 'rhoConst' ? ['    equationOfState', '    {', `        rho         ${ff(rho, 6)};`, '    }'] : []), '    thermodynamics', '    {', `        Cp          ${ff(cp, 5)};`, '        Hf          0;', '    }', '    transport', '    {', `        mu          ${ff(mu, 5)};`, `        Pr          ${ff(Pr, 4)};`, '    }', '}'];

/** Mesh and time plan shared by the Euler–Euler cases: O-grid of nPerD cells across, axial cells twice as long, mild wall grading. U is the transport velocity (m/s), peak the ratio of the largest phase velocity to it (sets the time step). */
function eulerPlan(c, opts, U, factor, peak = 1.3) {
  const D = c.D, nPerD = Math.round(clamp(num(opts.cellsPerDiameter, 16), 6, 200)), L = clamp(num(opts.lengthD, 20), 3, 2000) * D, maxCo = 0.5;
  const nc = Math.max(3, Math.round(0.4 * nPerD)), nr = Math.max(2, Math.round((nPerD - nc) / 2)), dz = (2 * D) / nPerD, mesh = ogridMesh({ D, legs: [{ type: 'straight', L }], nc, nr, wallRatio: 0.5, dz });
  const Ut = Math.max(U, 1e-3), endTime = pos(opts.endTime, (pos(opts.flowThroughs, 5) * L) / Ut), dt = (maxCo * dz) / (peak * Ut), steps = Math.max(1, Math.ceil(endTime / dt)), coreHours = (mesh.nCells * steps * COST_PER_CELL_STEP * factor) / 3600;
  const sections = [0.25, 0.5, 0.75].map((f) => ({ c: [f * L, 0, 0], t: [1, 0, 0], e1: [0, 1, 0], e2: [0, 0, 1], s: f * L }));
  return { mesh, nPerD, nc, nr, dz, L, maxCo, endTime, dt, steps, coreHours, factor, sections, nProcs: Math.max(1, Math.round(num(opts.nProcs, Math.min(64, Math.max(1, Math.round(mesh.nCells / 50000)))))), every: 5 };
}
/** Function objects of the Euler–Euler cases; same names as the VOF cases, `first` (a liquid or gas fraction) leads every field list. */
function eulerFunctions({ mesh, sections, every, first, more = [], pre = [], outletMore = [], volumeMore = [] }) {
  const ctl = ['        writeFields     false;', '        log             false;', '        writeControl    timeStep;', `        writeInterval   ${every};`], fl = (a) => `        fields          (${a.join(' ')});`, b = [mesh.D, mesh.D, mesh.D].map((v) => 0.8 * v);
  const out = ['functions', '{', ...pre, '    holdupVolume', '    {', '        type            volFieldValue;', '        libs            (fieldFunctionObjects);', '        operation       volAverage;', fl([first, ...more, ...volumeMore]), '        regionType      all;', ...ctl, '    }'];
  sections.forEach((p, i) => out.push(`    section${i + 1}`, '    {', '        type            surfaceFieldValue;', '        libs            (fieldFunctionObjects);', '        regionType      sampledSurface;', `        name            section${i + 1};`,
    '        sampledSurfaceDict', '        {', '            type        plane;', `            point       ${vec(p.c)};`, `            normal      ${vec(p.t)};`, `            bounds      ${vec(sub(p.c, b))} ${vec(add(p.c, b))};`, '            triangulate false;', '        }',
    '        operation       areaAverage;', fl([first, 'p', ...more]), '        surfaceFormat   none;', ...ctl, '    }'));
  out.push('    outletHoldup', '    {', '        type            surfaceFieldValue;', '        libs            (fieldFunctionObjects);', '        regionType      patch;', '        name            outlet;', '        operation       areaAverage;', fl([first, ...more, ...outletMore]), ...ctl, '    }',
    '    probes', '    {', '        type            probes;', '        libs            (sampling);', fl(['p', 'p_rgh', first]), '        probeLocations', '        (', sections.map((p) => `            ${vec(add(p.c, add(mul(p.e1, 0.07 * mesh.D), mul(p.e2, 0.07 * mesh.D))))}`), '        );', '        writeControl    timeStep;', `        writeInterval   ${every};`, '    }', '}');
  return out;
}
const eulerControl = (solver, P, title, functions) => dictFile('controlDict', 'system', title, [`application     ${solver};`, '', 'startFrom       latestTime;', 'startTime       0;', 'stopAt          endTime;', `endTime         ${ff(P.endTime, 6)};`, `deltaT          ${ff(P.dt / 5, 3)};`, '', 'writeControl    adjustable;', `writeInterval   ${ff(P.endTime / 20, 4)};`, 'purgeWrite      0;', 'writeFormat     binary;', 'writePrecision  8;', 'writeCompression off;', 'timeFormat      general;', 'timePrecision   8;', 'runTimeModifiable yes;', '', 'adjustTimeStep  on;', `maxCo           ${ff(P.maxCo, 3)};`, `maxDeltaT       ${ff(Math.min(P.endTime / 20, 4 * P.dt), 3)};`, '', ...functions]);
const eulerGravity = (c, title) => { const th = (c.angleDeg * Math.PI) / 180; return { path: 'constant/g', text: text([foamHeader('uniformDimensionedVectorField', 'g', title, 'constant'), 'dimensions      [0 1 -2 0 0 0 0];', `value           ${vec([-G * Math.sin(th), 0, -G * Math.cos(th)], 7)};`]) }; };
const eulerShell = (P, title) => [{ path: 'system/decomposeParDict', text: decomposeDict(P.nProcs, title) }, { path: 'Allrun', text: allrun(['blockMesh', '!restore0Dir', 'checkMesh'], P.nProcs) }, { path: 'Allclean', text: ALLCLEAN }, { path: 'collect.sh', text: COLLECT }, { path: 'case.foam', text: '' }];
const eulerFlavour = [...FLAVOUR_NOTE.slice(0, -1).filter((l) => !/incompressibleVoF|transportProperties|createPatch/.test(l)), '- the Euler–Euler solvers of this release (`reactingTwoPhaseEulerFoam`, `multiphaseEulerFoam`) correspond to `multiphaseEulerFoam` (Foundation v8–v10) and `foamRun -solver multiphaseEuler` (v11 and later); there `constant/phaseProperties` keeps the same sub-models under partly different keywords (`blending`, `drag`, `virtualMass`, `lift`, `wallLubrication`, `turbulentDispersion`, `populationBalanceCoeffs`), and the per-phase files are `physicalProperties.<phase>` and `momentumTransport.<phase>`.', ''];
const eulerRunBlock = (name, solver, P) => ['## Run', '', '```sh', 'source /path/to/OpenFOAM-v2412/etc/bashrc   # or: openfoam2412', `unzip ${name}.zip && cd ${name}`, './Allrun            # mesh, initial fields, mesh check, solver, result collection', '# ./Allclean       # removes everything that Allrun created', '```', '', `\`Allrun\` runs blockMesh → restore0Dir (copies 0.orig to 0) → checkMesh → ${P.nProcs > 1 ? `decomposePar → mpirun -np ${P.nProcs} ` : ''}${solver}${P.nProcs > 1 ? ' -parallel → reconstructPar' : ''} → collect.sh.`, 'Change `numberOfSubdomains` and `n` in system/decomposeParDict to use another number of cores.', ''];
const eulerMeshRows = (P, U) => [['Block structure', `O-grid: square core + 4 wall blocks, ${P.mesh.sectionCells} cells per section (${P.nPerD} across the diameter), wall cells half as thick as the inner ones`], ['Axial cell length', `${ff(P.dz * 1000, 4)} mm`], ['Section length', `${ff(P.L, 4)} m (${ff(P.L / P.mesh.D, 4)} diameters)`], ['Total cells', P.mesh.nCells.toLocaleString('en-US')], ['Simulated time', `${ff(P.endTime, 4)} s (${ff((P.endTime * U) / P.L, 3)} flow-through times at ${ff(U, 4)} m/s)`], ['Time step at Courant number ' + ff(P.maxCo, 2) + ' (estimate)', `${ff(P.dt, 3)} s, about ${P.steps.toLocaleString('en-US')} steps`], ['Estimated cost', `${hours(P.coreHours)} — ${ff(COST_PER_CELL_STEP * P.factor * 1e6, 3)} core-µs per cell and step, i.e. ${ff(P.factor, 3)} × the interFoam figure (${ff(COST_PER_CELL_STEP * 1e6, 3)} core-µs), measured on short runs of this set-up (5,040 cells, one workstation core); time the first steps and rescale`]];
const eulerImportNote = ['`collect.sh` gathers them in `hydraslug_results.txt`. On the External solvers page press “Import results” and choose that file (or the individual files): the application', 'computes the mean phase fraction and the pressure gradient between the sections and sets them beside its own one-dimensional values.', ''];

/** Fields, properties and numerics common to the reactingTwoPhaseEulerFoam cases. cfg: { title, c, P, d (dispersed phase), aD, U, TK, turb: { phase: { model, coeffs, ke, granular } }, thermo: { phase: lines }, phaseProps, solids, sizeGroups, popName, wallU, wallTheta, functions }. */
function eulerTwoFluidFiles(cfg) {
  const { title, c, P, d, aD, U, TK } = cfg, p0 = ff(c.P * 1e5), Ui = `(${ff(U, 6)} 0 0)`, files = [], k = Math.max(1.5 * (0.05 * U) ** 2, 1e-8), eps = (0.09 ** 0.75 * k ** 1.5) / (0.07 * c.D), T0 = ff(TK, 6);
  files.push(eulerFld('volScalarField', 'alpha.' + d, eulerDim.one, ff(aD, 6), eulerFv(ff(aD, 6)), eulerIo(ff(aD, 6), d), eulerZg, title));
  files.push(eulerFld('volScalarField', 'p', eulerDim.p, p0, eulerCalc(p0), eulerCalc(p0), eulerCalc(p0), title));
  files.push(eulerFld('volScalarField', 'p_rgh', eulerDim.p, p0, eulerWf('fixedFluxPressure', p0), eulerFv(p0), eulerWf('fixedFluxPressure', p0), title));
  for (const ph of [d, 'liquid']) {
    const t = cfg.turb[ph], wallU = ph === d && cfg.wallU ? cfg.wallU : ['type            noSlip;'];
    files.push(eulerFld('volVectorField', 'U.' + ph, eulerDim.U, Ui, eulerFv(Ui), ['type            pressureInletOutletVelocity;', `phi             phi.${ph};`, `value           uniform ${Ui};`], wallU, title));
    files.push(eulerFld('volScalarField', 'T.' + ph, eulerDim.T, T0, eulerFv(T0), eulerIo(T0, ph), eulerZg, title));
    files.push(eulerFld('volScalarField', 'alphat.' + ph, eulerDim.alphat, '0', eulerCalc('0'), eulerCalc('0'), t.ke ? eulerWf('compressible::alphatWallFunction', '0', ['Prt             0.85;']) : eulerCalc('0'), title));
    if (t.ke) {
      files.push(eulerFld('volScalarField', 'k.' + ph, eulerDim.k, ff(k, 5), eulerFv(ff(k, 5)), eulerIo(ff(k, 5), ph), eulerWf('kqRWallFunction', ff(k, 5)), title));
      files.push(eulerFld('volScalarField', 'epsilon.' + ph, eulerDim.eps, ff(eps, 5), eulerFv(ff(eps, 5)), eulerIo(ff(eps, 5), ph), eulerWf('epsilonWallFunction', ff(eps, 5)), title));
      files.push(eulerFld('volScalarField', 'nut.' + ph, eulerDim.nu, '1e-08', eulerCalc('1e-08'), eulerCalc('1e-08'), eulerWf('nutkWallFunction', '1e-08'), title));
    } else if (t.granular) {
      files.push(eulerFld('volScalarField', 'nut.' + ph, eulerDim.nu, '0', eulerCalc('0'), eulerCalc('0'), eulerCalc('0'), title));
      files.push(eulerFld('volScalarField', 'Theta.' + ph, eulerDim.k, '1e-04', eulerFv('1e-04'), eulerZg, cfg.wallTheta || eulerZg, title));
    }
    files.push({ path: 'constant/turbulenceProperties.' + ph, text: dictFile('turbulenceProperties.' + ph, 'constant', title, eulerRas(t.model, t.coeffs)) });
    files.push({ path: 'constant/thermophysicalProperties.' + ph, text: dictFile('thermophysicalProperties.' + ph, 'constant', title, cfg.thermo[ph]) });
  }
  if (cfg.sizeGroups) files.push(eulerFld('volScalarField', `f.${d}.${cfg.popName}`, eulerDim.one, '0', eulerFv('0'), eulerIo('0', d), eulerZg, title));
  files.push({ path: 'constant/phaseProperties', text: dictFile('phaseProperties', 'constant', title, cfg.phaseProps) }, eulerGravity(c, title), { path: 'system/blockMeshDict', text: blockMeshDict(P.mesh, title, 'inlet') });
  files.push({ path: 'system/fvSchemes', text: dictFile('fvSchemes', 'system', title, ['ddtSchemes', '{', '    default         Euler;', '}', '', 'gradSchemes', '{', '    default         Gauss linear;', '}', '', 'divSchemes', '{', '    default         none;', '    "div\\(phi,alpha.*\\)" Gauss vanLeer;', '    "div\\(phir,alpha.*\\)" Gauss vanLeer;', '    "div\\(alphaRhoPhi.*,U.*\\)" Gauss limitedLinearV 1;', '    "div\\(phi.*,U.*\\)" Gauss limitedLinearV 1;', '    "div\\(alphaRhoPhi.*,(h|e|f).*\\)" Gauss limitedLinear 1;', '    "div\\(alphaRhoPhi.*,K.*\\)" Gauss limitedLinear 1;', '    "div\\(alphaPhi.*,p\\)" Gauss limitedLinear 1;', '    "div\\(alphaRhoPhi.*,(k|epsilon|Theta).*\\)" Gauss upwind;', '    "div\\(phim,(k|epsilon)m\\)" Gauss upwind;',
    '    "div\\(\\(\\(\\(alpha.*\\*thermo:rho.*\\)\\*nuEff.*\\)\\*dev2\\(T\\(grad\\(U.*\\)\\)\\)\\)\\)" Gauss linear;', '    "div\\(\\(\\(\\(thermo:rho.*\\*nut.*\\)\\*dev2\\(T\\(grad\\(U.*\\)\\)\\)\\)\\+\\(\\(\\(thermo:rho.*\\*lambda.*\\)\\*div\\(phi.*\\)\\)\\*I\\)\\)\\)" Gauss linear;', '}', '', 'laplacianSchemes', '{', '    default         Gauss linear corrected;', '    bounded         Gauss linear corrected;', '}', '', 'interpolationSchemes', '{', '    default         linear;', '}', '', 'snGradSchemes', '{', '    default         corrected;', '    bounded         corrected;', '}', '', 'wallDist', '{', '    method          meshWave;', '    nRequired       yes;', '}']) });
  const lin = (name, tol, extra = []) => [`    "${name}"`, '    {', '        solver          smoothSolver;', '        smoother        symGaussSeidel;', `        tolerance       ${tol};`, '        relTol          0;', '        minIter         1;', ...extra, '    }', ''];
  files.push({ path: 'system/fvSolution', text: dictFile('fvSolution', 'system', title, ['solvers', '{', '    "alpha.*"', '    {', '        nAlphaCorr      1;', '        nAlphaSubCycles 2;', ...(cfg.solids ? ['        implicitPhasePressure yes;', '        solver          smoothSolver;', '        smoother        symGaussSeidel;', '        tolerance       1e-9;', '        relTol          0;', '        minIter         1;'] : []), '    }', '',
    cfg.sizeGroups ? [`    ${cfg.popName}`, '    {', '        nCorr           1;', '        tolerance       1e-4;', '        renormalize     false;', '        solveOnFinalIterOnly true;', '    }', '', '    "f.*"', '    {', '        solver          PBiCGStab;', '        preconditioner  DILU;', '        tolerance       1e-6;', '        relTol          0;', '    }', ''] : null,
    '    p_rgh', '    {', '        solver          GAMG;', '        smoother        DIC;', '        tolerance       1e-8;', '        relTol          0.01;', '    }', '', '    p_rghFinal', '    {', '        $p_rgh;', '        relTol          0;', '    }', '', lin('U.*', '1e-7'), lin('(h|e).*', '1e-7', ['        maxIter         20;']), lin('(k|epsilon|Theta).*', '1e-7'), '}', '', 'PIMPLE', '{', '    nOuterCorrectors 3;', '    nCorrectors      1;', '    nNonOrthogonalCorrectors 1;', cfg.solids ? '    faceMomentum     yes;' : null, '}', '', 'relaxationFactors', '{', '    equations', '    {', '        ".*"            1;', '    }', '}']) });
  files.push({ path: 'system/controlDict', text: eulerControl('reactingTwoPhaseEulerFoam', P, title, cfg.functions) }, ...eulerShell(P, title));
  return files;
}
const eulerPhase = (name, diameter, extra = []) => [name, '{', '    type            purePhaseModel;', ...diameter.map((l) => '    ' + l), ...extra.map((l) => '    ' + l), '}', ''];
const eulerConstD = (dm) => ['diameterModel   constant;', 'constantCoeffs', '{', `    d               ${ff(dm, 5)};`, '}'];
const eulerKt = ['kineticTheoryCoeffs', '{', '    equilibrium             off;', '    e                       0.8;', '    alphaMax                0.62;', '    alphaMinFriction        0.5;', '    residualAlpha           1e-4;', '    viscosityModel          Gidaspow;', '    conductivityModel       Gidaspow;', '    granularPressureModel   Lun;', '    frictionalStressModel   JohnsonJackson;', '    radialModel             SinclairJackson;', '    JohnsonJacksonCoeffs', '    {', '        Fr                      0.05;', '        eta                     2;', '        p                       5;', '        phi                     28.5;', '        alphaDeltaMin           0.05;', '    }', '}'];

/** Gas–liquid two-fluid case (reactingTwoPhaseEulerFoam, basicTwoPhaseSystem). */
function eulerGasLiquid(c, opts) {
  const vm = c.vsl + c.vsg, aG = clamp(c.vsg / vm, 1e-4, 1 - 1e-4), P = eulerPlan(c, opts, vm, eulerCost.gas, 2), TK = c.T + KEL, dB = clamp(num(opts.bubbleMm, 3), 0.05, 100) / 1000, dD = clamp(num(opts.dropMm, 0.5), 0.01, 50) / 1000;
  const cl = { lift: true, virtualMass: true, wallLubrication: true, turbulentDispersion: true, bubbleInducedTurbulence: true, ...(opts.closures || {}) }, name = safeName(opts.name || 'euler_two_fluid'), title = `Euler-Euler two-fluid model (gas-liquid), straight section at x = ${ff(c.x, 6)} m`;
  const mwEff = ((c.rhoG * R * TK) / (c.P * 1e5)) * 1000, cpG = pos(c.cpG, 2400), cpL = pos(c.cpL, 2100), gl = 'gas in liquid', lg = 'liquid in gas';
  const phaseProps = ['type            basicTwoPhaseSystem;', '', 'phases          (gas liquid);', '',
    ...eulerPhase('gas', ['diameterModel   isothermal;', 'isothermalCoeffs', '{', `    d0              ${ff(dB, 5)};`, `    p0              ${ff(c.P * 1e5)};`, '}'], ['residualAlpha   1e-6;']), ...eulerPhase('liquid', eulerConstD(dD), ['residualAlpha   1e-6;']),
    'blending', '{', '    default', '    {', '        type            linear;', '        minFullyContinuousAlpha.gas 0.7;', '        minPartlyContinuousAlpha.gas 0.3;', '        minFullyContinuousAlpha.liquid 0.7;', '        minPartlyContinuousAlpha.liquid 0.3;', '    }', '    drag', '    {', '        type            linear;', '        minFullyContinuousAlpha.gas 0.7;', '        minPartlyContinuousAlpha.gas 0.5;', '        minFullyContinuousAlpha.liquid 0.7;', '        minPartlyContinuousAlpha.liquid 0.5;', '    }', '}', '',
    ...eulerList('surfaceTension', [eulerPair('gas and liquid', ['type            constant;', `sigma           ${ff(c.sigma, 5)};`])]),
    ...eulerList('aspectRatio', [eulerPair(gl, ['type            constant;', 'E0              1;']), eulerPair(lg, ['type            constant;', 'E0              1;'])]),
    ...eulerList('drag', [eulerPair(gl, eulerDrag('IshiiZuber')), eulerPair(lg, eulerDrag('SchillerNaumann')), eulerPair('gas and liquid', ['type            segregated;', 'm               0.5;', 'n               8;', 'swarmCorrection', '{', '    type        none;', '}'])]),
    // the entry stays when the force is switched off (coefficient 0): the gas turbulence model looks the virtual-mass model up
    ...eulerList('virtualMass', [eulerPair(gl, ['type            constantCoefficient;', `Cvm             ${cl.virtualMass ? 0.5 : 0};`]), eulerPair(lg, ['type            constantCoefficient;', `Cvm             ${cl.virtualMass ? 0.5 : 0};`])]),
    ...eulerList('heatTransfer', [eulerPair(gl, ['type            RanzMarshall;', 'residualAlpha   1e-4;']), eulerPair(lg, ['type            RanzMarshall;', 'residualAlpha   1e-4;'])]),
    ...eulerList('lift', cl.lift ? [eulerPair(gl, ['type            Tomiyama;'])] : []),
    ...eulerList('phaseTransfer', []),
    ...eulerList('wallLubrication', cl.wallLubrication ? [eulerPair(gl, ['type            Antal;', 'Cw1             -0.01;', 'Cw2             0.05;'])] : []),
    ...eulerList('turbulentDispersion', cl.turbulentDispersion ? [eulerPair(gl, ['type            Burns;', 'sigma           0.7;', 'Ctd             1;', 'residualAlpha   1e-3;'])] : []),
    '// Minimum allowable pressure', 'pMin            10000;'];
  const first = 'alpha.liquid', bit = !!cl.bubbleInducedTurbulence;
  const files = eulerTwoFluidFiles({ title, c, P, d: 'gas', aD: aG, U: vm, TK, phaseProps,
    turb: { gas: { model: 'continuousGasKEpsilon', ke: true }, liquid: { model: bit ? 'LaheyKEpsilon' : 'kEpsilon', ke: true } },
    thermo: { gas: ['// molWeight is an effective value: with it the perfect-gas law gives the gas density of the case', `// (${ff(c.rhoG, 5)} kg/m3 at ${ff(c.P, 5)} bara, ${ff(c.T, 4)} degC), i.e. it contains the compressibility factor.`, ...eulerThermo('perfectGas', mwEff, cpG, c.muG, (cpG * c.muG) / pos(c.kG, 0.04))], liquid: eulerThermo('rhoConst', 100, cpL, c.muL, (cpL * c.muL) / pos(c.kL, 0.13), c.rhoL) },
    functions: eulerFunctions({ mesh: P.mesh, sections: P.sections, every: P.every, first, more: ['alpha.gas'] }) });
  const closures = [['Interfacial drag', 'Ishii–Zuber (`IshiiZuber`) for bubbles in liquid, Schiller–Naumann (`SchillerNaumann`) for drops in gas, `segregated` where both phases are continuous; linear blending between 30/50 % and 70 %'], ['Virtual mass', cl.virtualMass ? '`constantCoefficient`, Cvm = 0.5' : 'off (`constantCoefficient` with Cvm = 0)'], ['Lift', cl.lift ? 'Tomiyama (`Tomiyama`), bubbles in liquid' : 'off'], ['Wall lubrication', cl.wallLubrication ? 'Antal (`Antal`, Cw1 = −0.01, Cw2 = 0.05)' : 'off'], ['Turbulent dispersion', cl.turbulentDispersion ? 'Burns, Favre-averaged drag (`Burns`, σ = 0.7)' : 'off'], ['Turbulence, liquid', bit ? 'k–ε with bubble-induced production (`LaheyKEpsilon`)' : 'standard k–ε (`kEpsilon`)'], ['Turbulence, gas', 'k–ε response to the liquid turbulence (`continuousGasKEpsilon`)'], ['Bubble / drop diameter', `${ff(dB * 1000, 4)} mm at ${ff(c.P, 5)} bara (\`isothermal\`: scales with pressure) / ${ff(dD * 1000, 4)} mm`], ['Surface tension', `${ff(c.sigma, 4)} N/m`], ['Thermophysics', `gas: perfect gas with effective molar mass ${ff(mwEff, 5)} g/mol (density of the case at P, T); liquid: constant density`]];
  return { name, title, files, P, U: vm, solver: 'reactingTwoPhaseEulerFoam', label: 'two-fluid gas–liquid', first, invert: false, closures,
    inlet: `both phases at the mixture velocity ${ff(vm, 4)} m/s, gas fraction ${ff(aG, 4)} (no-slip split vsg/vm); the slip and the holdup develop along the section`, more: 'alpha.gas',
    intro: ['The one-dimensional models of the application give section-averaged holdup and pressure gradient from a drift-flux or two-fluid closure with one velocity per phase and section.', 'This case solves the two-fluid equations in three dimensions: each phase has its own velocity field, and the interfacial forces (drag, lift, virtual mass, wall lubrication,', 'turbulent dispersion) decide how the gas distributes over the cross-section — phase segregation, wall peaking or coring of the void, and the resulting slip.'],
    notes: ['The dispersed-phase diameters are inputs (no interfacial-area transport): results for separated regimes (stratified, annular) depend on the `segregated` drag and on the blending limits.', aG > 0.7 || aG < 0.3 ? null : 'At this gas fraction both phases are partly continuous: the blending between bubble and drop closures is active across the section.'] };
}

/** Liquid + solid particles (kinetic theory of granular flow) with reactingTwoPhaseEulerFoam. */
function eulerSolids(c, opts) {
  const U = c.vsl / c.holdup, aS = clamp(num(opts.solidsFraction, 0.05), 1e-4, 0.5), dP = clamp(num(opts.particleUm, 300), 1, 20000) * 1e-6, rhoS = clamp(num(opts.solidsDensity, 920), 100, 20000), P = eulerPlan(c, opts, U, eulerCost.solids), TK = c.T + KEL, cpL = pos(c.cpL, 2100);
  const name = safeName(opts.name || 'euler_solids'), title = `Euler-Euler liquid-solids model (kinetic theory of granular flow), straight section at x = ${ff(c.x, 6)} m`, sl = 'solids in liquid';
  const phaseProps = ['type            basicTwoPhaseSystem;', '', 'phases          (solids liquid);', '', ...eulerPhase('solids', eulerConstD(dP), ['alphaMax        0.62;', 'residualAlpha   1e-6;']), ...eulerPhase('liquid', eulerConstD(1), ['residualAlpha   0;']),
    'blending', '{', '    default', '    {', '        type            none;', '        continuousPhase liquid;', '    }', '}', '',
    ...eulerList('surfaceTension', [eulerPair('liquid and solids', ['type            constant;', 'sigma           0;'])]), ...eulerList('aspectRatio', []),
    ...eulerList('drag', [eulerPair(sl, eulerDrag('GidaspowErgunWenYu'))]), ...eulerList('virtualMass', [eulerPair(sl, ['type            constantCoefficient;', 'Cvm             0.5;'])]),
    ...eulerList('heatTransfer', [eulerPair(sl, ['type            RanzMarshall;', 'residualAlpha   1e-4;'])]), ...eulerList('lift', []), ...eulerList('phaseTransfer', []), ...eulerList('wallLubrication', []),
    ...eulerList('turbulentDispersion', [eulerPair(sl, ['type            Burns;', 'sigma           0.7;', 'Ctd             1;', 'residualAlpha   1e-3;'])]), '// Minimum allowable pressure', 'pMin            10000;'];
  const first = 'alpha.liquid';
  const jj = ['restitutionCoefficient 0.2;', 'specularityCoefficient 0.1;', 'muF             0.25;', 'sigma           2;'];
  const files = eulerTwoFluidFiles({ title, c, P, d: 'solids', aD: aS, U, TK, phaseProps, solids: true, wallU: eulerWf('JohnsonJacksonParticleSlip', '(0 0 0)', jj), wallTheta: eulerWf('JohnsonJacksonParticleTheta', '1e-04', jj),
    turb: { solids: { model: 'kineticTheory', granular: true, coeffs: [...eulerKt, 'phasePressureCoeffs', '{', '    preAlphaExp     500;', '    expMax          1000;', '    alphaMax        0.62;', '    g0              1000;', '}'] }, liquid: { model: 'kEpsilon', ke: true } },
    thermo: { solids: eulerThermo('rhoConst', 100, 2100, 0, 1, rhoS), liquid: eulerThermo('rhoConst', 100, cpL, c.muL, (cpL * c.muL) / pos(c.kL, 0.13), c.rhoL) },
    functions: eulerFunctions({ mesh: P.mesh, sections: P.sections, every: P.every, first, more: ['alpha.solids'] }) });
  const closures = [['Interfacial drag', 'Gidaspow blend of Ergun (dense) and Wen–Yu (dilute) (`GidaspowErgunWenYu`)'], ['Virtual mass', '`constantCoefficient`, Cvm = 0.5'], ['Turbulent dispersion', 'Burns (`Burns`, σ = 0.7)'], ['Solids stress', 'kinetic theory of granular flow (`kineticTheory`), granular temperature transported (`equilibrium off`), restitution coefficient 0.8'], ['Granular viscosity / conductivity', 'Gidaspow (`Gidaspow`)'], ['Granular pressure', 'Lun (`Lun`)'], ['Radial distribution', 'Sinclair–Jackson (`SinclairJackson`), packing limit 0.62'], ['Frictional stress', 'Johnson–Jackson (`JohnsonJackson`) above a solids fraction of 0.5'], ['Turbulence, liquid', 'standard k–ε (`kEpsilon`)'], ['Particles', `${ff(dP * 1e6, 4)} µm, ${ff(rhoS, 5)} kg/m³, ${ff(100 * aS, 3)} % by volume at the inlet`]];
  return { name, title, files, P, U, solver: 'reactingTwoPhaseEulerFoam', label: 'liquid–solids (kinetic theory)', first, invert: false, closures, extraRef: { solidsFraction: aS },
    inlet: `liquid and particles at the liquid velocity of the 1-D model, vsl/holdup = ${ff(U, 4)} m/s, solids fraction ${ff(aS, 4)}; the pipe runs full of the slurry (the gas of the case is not part of this model)`, more: 'alpha.solids',
    intro: ['The solids models of the application transport hydrate or sand as a section-averaged concentration with correlations for deposition velocity and bed formation.', 'This case treats the particles as an interpenetrating continuum with its own momentum equation and kinetic-theory closures (granular temperature, collisional pressure and', 'viscosity, frictional stress near packing): it shows where the particles travel in the cross-section, whether a moving or stationary bed forms and how it loads the wall.'],
    notes: [`Particle density ${ff(rhoS, 5)} kg/m³ against ${ff(c.rhoL, 5)} kg/m³ of liquid: the particles ${rhoS > c.rhoL ? 'settle toward the bottom' : 'rise toward the top'} of the pipe.`, 'Cohesion between hydrate particles is not part of the kinetic theory; use the population-balance case for agglomeration.'] };
}

/** Gas, oil and water as three fluids with multiphaseEulerFoam. */
function eulerThreeFluid(c, opts) {
  const vm = c.vsl + c.vsg, wc = clamp(num(c.wcut, 0), 0, 1), aG = +ff(clamp(c.vsg / vm, 1e-4, 1 - 1e-4), 6), aW = +ff((wc * c.vsl) / vm, 6), aO = +ff(Math.max(1 - aG - aW, 0), 6), P = eulerPlan(c, opts, vm, eulerCost.three, 2);
  const name = safeName(opts.name || 'euler_three_fluid'), title = `Euler-Euler multi-fluid model (gas, oil, water), straight section at x = ${ff(c.x, 6)} m`, dB = clamp(num(opts.bubbleMm, 3), 0.05, 100) / 1000, dD = clamp(num(opts.dropMm, 0.5), 0.01, 50) / 1000;
  const rhoO = pos(c.rhoO, c.rhoL), rhoW = pos(c.rhoW, 1000), muO = pos(c.muO, c.muL), muW = pos(c.muW, 5e-4), sGW = pos(opts.sigmaGasWater, Math.max(c.sigma, 0.06)), sOW = pos(opts.sigmaOilWater, 0.025), cA = clamp(num(opts.interfaceCompression, 0), 0, 1.5), p0 = ff(c.P * 1e5), Ui = `(${ff(vm, 6)} 0 0)`;
  const ph = [['gas', aG, c.rhoG, c.muG, pos(c.kG, 0.04), pos(c.cpG, 2400), dB], ['oil', aO, rhoO, muO, pos(c.kL, 0.13), pos(c.cpL, 2100), dD], ['water', aW, rhoW, muW, 0.6, 4180, dD]], pairs = [['gas', 'oil'], ['gas', 'water'], ['oil', 'water']], files = [];
  for (const [n, a] of ph) {
    files.push(eulerFld('volScalarField', 'alpha.' + n, eulerDim.one, ff(a, 6), eulerFv(ff(a, 6)), eulerIo(ff(a, 6), n), eulerZg, title));
    files.push(eulerFld('volVectorField', 'U.' + n, eulerDim.U, Ui, eulerFv(Ui), ['type            pressureInletOutletVelocity;', `phi             phi.${n};`, `value           uniform ${Ui};`], ['type            noSlip;'], title));
  }
  files.push(eulerFld('volScalarField', 'p_rgh', eulerDim.p, p0, eulerWf('fixedFluxPressure', p0), eulerFv(p0), eulerWf('fixedFluxPressure', p0), title));
  const k = Math.max(1.5 * (0.05 * vm) ** 2, 1e-8), eps = (0.09 ** 0.75 * k ** 1.5) / (0.07 * c.D), mixIo = (v) => ['type            inletOutlet;', 'phi             phi;', `inletValue      uniform ${v};`, `value           uniform ${v};`];
  files.push(eulerFld('volScalarField', 'k', eulerDim.k, ff(k, 5), eulerFv(ff(k, 5)), mixIo(ff(k, 5)), eulerWf('kqRWallFunction', ff(k, 5)), title), eulerFld('volScalarField', 'epsilon', eulerDim.eps, ff(eps, 5), eulerFv(ff(eps, 5)), mixIo(ff(eps, 5)), eulerWf('epsilonWallFunction', ff(eps, 5)), title), eulerFld('volScalarField', 'nut', eulerDim.nu, '0', eulerCalc('0'), eulerCalc('0'), eulerWf('nutkWallFunction', '0'), title));
  const blended = ([a, b]) => [`    (${a} ${b})`, '    {', '        type blended;', ...[a, b].flatMap((n) => [`        ${n}`, '        {', '            type SchillerNaumann;', '            residualPhaseFraction 0;', '            residualSlip 0;', '        }']), '        residualPhaseFraction 1e-3;', '        residualSlip 1e-3;', '    }'];
  files.push({ path: 'constant/transportProperties', text: dictFile('transportProperties', 'constant', title, ['phases', '(', ...ph.flatMap(([n, , rho, mu, kap, cp, d]) => [`    ${n}`, '    {', `        rho             ${ff(rho, 6)};`, `        nu              ${ff(mu / rho, 6)};`, `        kappa           ${ff(kap, 5)};`, `        Cp              ${ff(cp, 5)};`, '        diameterModel   constant;', '        constantCoeffs', '        {', `            d               ${ff(d, 5)};`, '        }', '    }']), ');', '',
    'sigmas', '(', `    (gas oil)       ${ff(c.sigma, 5)}`, `    (gas water)     ${ff(sGW, 5)}`, `    (oil water)     ${ff(sOW, 5)}`, ');', '', 'interfaceCompression', '(', ...pairs.map(([a, b]) => `    (${a} ${b})     ${ff(cA, 3)}`), ');', '', 'virtualMass', '(', '    (gas oil)       0.5', '    (gas water)     0.5', ');', '', 'drag', '(', ...pairs.flatMap(blended), ');', '', '// single-phase transport entry read by the mixture turbulence model (laminar viscosity of its wall functions): the liquid of the case', 'transportModel  Newtonian;', `nu              ${ff(c.muL / c.rhoL, 6)};`]) });
  files.push({ path: 'constant/turbulenceProperties', text: dictFile('turbulenceProperties', 'constant', title, eulerRas('kEpsilon')) }, eulerGravity(c, title), { path: 'system/blockMeshDict', text: blockMeshDict(P.mesh, title, 'inlet') });
  files.push({ path: 'system/fvSchemes', text: dictFile('fvSchemes', 'system', title, ['ddtSchemes', '{', '    default         Euler;', '}', '', 'gradSchemes', '{', '    default         Gauss linear;', '}', '', 'divSchemes', '{', '    "div\\(phi,alpha.*\\)" Gauss vanLeer;', '    "div\\(phir,alpha.*,alpha.*\\)" Gauss vanLeer;', '    "div\\(alphaPhi.*,U.*\\)" Gauss limitedLinearV 1;', '    div(Rc)         Gauss linear;', '    "div\\(phi.*,U.*\\)" Gauss limitedLinearV 1;', '    "div\\(phi,(k|epsilon)\\)" Gauss upwind;', '    "div\\(\\(nuEff.*" Gauss linear;', '}', '', 'laplacianSchemes', '{', '    default         Gauss linear corrected;', '}', '', 'interpolationSchemes', '{', '    default         linear;', '}', '', 'snGradSchemes', '{', '    default         corrected;', '}']) });
  files.push({ path: 'system/fvSolution', text: dictFile('fvSolution', 'system', title, ['solvers', '{', '    "alpha.*"', '    {', '        nAlphaSubCycles 2;', '    }', '', '    p_rgh', '    {', '        solver          GAMG;', '        smoother        DIC;', '        tolerance       1e-8;', '        relTol          0.01;', '    }', '', '    p_rghFinal', '    {', '        $p_rgh;', '        tolerance       1e-9;', '        relTol          0;', '    }', '', '    "pcorr.*"', '    {', '        $p_rgh;', '        tolerance       1e-5;', '        relTol          0;', '    }', '', '    "(U|T|Theta|k|epsilon).*"', '    {', '        solver          smoothSolver;', '        smoother        symGaussSeidel;', '        tolerance       1e-7;', '        relTol          0;', '    }', '}', '', 'PIMPLE', '{', '    nOuterCorrectors 1;', '    nCorrectors      3;', '    nNonOrthogonalCorrectors 1;', '}', '', 'relaxationFactors', '{', '    equations', '    {', '        ".*"            1;', '    }', '}']) });
  const first = 'alpha.liquid', pre = ['    liquidFraction', '    {', '        // alpha.liquid = alpha.oil + alpha.water, registered for the function objects below', '        type            add;', '        libs            (fieldFunctionObjects);', '        fields          (alpha.oil alpha.water);', '        result          alpha.liquid;', '        executeControl  timeStep;', '        executeInterval 1;', '        writeControl    none;', '        log             false;', '    }'];
  files.push({ path: 'system/controlDict', text: eulerControl('multiphaseEulerFoam', P, title, eulerFunctions({ mesh: P.mesh, sections: P.sections, every: P.every, first, more: ['alpha.gas', 'alpha.oil', 'alpha.water'], pre })) }, ...eulerShell(P, title));
  const closures = [['Interfacial drag', 'Schiller–Naumann for each phase of a pair dispersed in the other, blended by phase fraction (`blended` / `SchillerNaumann`), for gas–oil, gas–water and oil–water'], ['Virtual mass', 'coefficient 0.5 for gas–oil and gas–water'], ['Surface tension', `gas–oil ${ff(c.sigma, 4)} N/m (case), gas–water ${ff(sGW, 4)} N/m, oil–water ${ff(sOW, 4)} N/m (options sigmaGasWater, sigmaOilWater)`], ['Interface compression', `coefficient ${ff(cA, 3)} for every pair (0 = dispersed multi-fluid; 1 sharpens the interfaces of separated layers; option interfaceCompression)`], ['Diameters', `gas bubbles ${ff(dB * 1000, 4)} mm, oil and water drops ${ff(dD * 1000, 4)} mm`], ['Turbulence', 'standard k–ε of the mixture (`kEpsilon`, one eddy viscosity shared by the three phases, wall functions); the molecular viscosity of each phase acts on that phase'], ['Phases', `gas ${ff(c.rhoG, 5)} kg/m³, ${ff(c.muG, 4)} Pa·s; oil ${ff(rhoO, 5)} kg/m³, ${ff(muO, 4)} Pa·s; water ${ff(rhoW, 5)} kg/m³, ${ff(muW, 4)} Pa·s (all incompressible)`]];
  return { name, title, files, P, U: vm, solver: 'multiphaseEulerFoam', label: 'multi-fluid gas–oil–water', first, invert: false, closures, extraRef: { inletFractions: { gas: aG, oil: aO, water: aW } },
    inlet: `all three phases at the mixture velocity ${ff(vm, 4)} m/s with the no-slip fractions gas ${ff(aG, 4)}, oil ${ff(aO, 4)}, water ${ff(aW, 4)} (water cut ${ff(wc, 3)})`, more: 'alpha.gas, alpha.oil, alpha.water',
    intro: ['The one-dimensional models of the application carry gas and one liquid; oil and water share a velocity and their split is the input water cut.', 'This case gives gas, oil and water each their own momentum equation: the water can lag behind or run ahead of the oil, settle into a layer at the bottom of the pipe, and the', 'local water fraction at the wall — what matters for corrosion and hydrate formation — comes out of the solution.'],
    notes: ['`alpha.liquid` in the result files is alpha.oil + alpha.water, built at run time by the `liquidFraction` function object.', 'The phases are incompressible in this solver; the outlet value of p_rgh is the pressure of the case so that p is absolute.'] };
}

function eulerAssemble(v, c) {
  const { P, name } = v, cells = P.mesh.nCells, files = [...v.files];
  const readme = text([`# ${name} — OpenFOAM hand-off`, '', `${v.title}: ${ff(P.L / c.D, 4)} diameters of pipe (${ff(P.L, 4)} m), solver \`${v.solver}\`.`, '', ...v.intro, '',
    '## Case data (taken from the HydraSlug case)', '', ...caseTable(c), ...kv([['Inlet', v.inlet], ['Outlet', `fixed pressure ${ff(c.P, 5)} bara (p_rgh = ${ff(c.P * 1e5)} Pa)`], ['Gravity', `resolved for the inclination of ${ff(c.angleDeg, 4)}°: x is the pipe axis, z the upward normal`]]),
    '## Closures', '', ...kv(v.closures), '## Mesh and cost', '', ...kv(eulerMeshRows(P, v.U)), (v.notes || []).filter(Boolean).map((n) => `- ${n}`), '',
    ...eulerRunBlock(name, v.solver, P), '## Results', '', `The function objects in system/controlDict write, every ${P.every} time steps:`, '',
    ...v.results, '', ...eulerImportNote, ...eulerFlavour]);
  return { name, files: [...files, { path: 'README.md', text: readme }], readme, summary: `${v.solver} ${v.label}: ${cells.toLocaleString('en-US')} cells, ${ff(P.endTime, 3)} s simulated, about ${hours(P.coreHours)}.`, mesh: P.mesh, sections: P.sections,
    plan: { cells, steps: P.steps, coreHours: P.coreHours, endTime: P.endTime, dt: P.dt, costFactor: P.factor, nProcs: P.nProcs }, commands: ['source /path/to/OpenFOAM-v2412/etc/bashrc', `unzip ${name}.zip && cd ${name}`, './Allrun'],
    reference: { holdup: c.holdup, dpdx: c.dpdx, tauW: c.tauW, slugFreq: null, forceN: null, sectionDistance: 0.5 * P.L, liquidField: v.first, invert: v.invert, ...(v.extraRef || {}) } };
}
const eulerResultLines = (first, more, firstText) => [`- \`postProcessing/holdupVolume/0/volFieldValue.dat\` — volume average of \`${first}\` (${firstText}) and of ${more};`, `- \`postProcessing/section1..3/0/surfaceFieldValue.dat\` — area average of \`${first}\`, the pressure \`p\` and ${more} on three cross-sections (25 %, 50 %, 75 % of the length);`, `- \`postProcessing/outletHoldup/0/surfaceFieldValue.dat\` — the same fractions averaged over the outlet;`, `- \`postProcessing/probes/0/p\`, \`p_rgh\`, \`${first}\` — point values near the centreline at the same sections.`];

/**
 * OpenFOAM Euler–Euler case of a straight pipe section: every phase is a continuum with its own velocity field.
 * opts: { case (from bridgeCase), phases: 2 | 3, dispersed: 'gas' | 'solids' (two phases only), lengthD, cellsPerDiameter, bubbleMm, dropMm,
 * particleUm, solidsFraction, solidsDensity, closures: { lift, virtualMass, wallLubrication, turbulentDispersion, bubbleInducedTurbulence },
 * sigmaGasWater, sigmaOilWater, interfaceCompression (three phases), endTime, flowThroughs, nProcs, name }.
 * Two phases: reactingTwoPhaseEulerFoam (gas–liquid two-fluid model, or liquid + particles with kinetic theory); three phases: multiphaseEulerFoam.
 * Returns { name, files, readme, summary, mesh, plan, commands, reference }.
 */
export function openfoamEulerCase(opts = {}) {
  const c = withCase(opts), v = +opts.phases === 3 ? eulerThreeFluid(c, opts) : opts.dispersed === 'solids' ? eulerSolids(c, opts) : eulerGasLiquid(c, opts);
  v.results = eulerResultLines(v.first, `\`${v.more.split(', ').join('`, `')}\``, 'the liquid fraction');
  return eulerAssemble(v, c);
}

/**
 * OpenFOAM population balance + CFD case (class method): reactingTwoPhaseEulerFoam with populationBalanceTwoPhaseSystem. The dispersed phase
 * (hydrate particles, or gas bubbles) is split into size groups that exchange volume by coalescence (agglomeration) and breakup while they are
 * carried through a straight pipe section by the liquid of the case.
 * opts: { case, dispersed: 'hydrate' | 'gas', nGroups, dMinUm, dMaxUm, dInUm, fraction, particleDensity, lengthD, cellsPerDiameter, endTime, flowThroughs, nProcs, name }.
 * Returns { name, files, readme, summary, mesh, plan, commands, reference (with groups: [{ name, field, d }]) }.
 */
export function openfoamPopulationCase(opts = {}) {
  const c = withCase(opts), gas = opts.dispersed === 'gas', d = gas ? 'gas' : 'particles', pop = gas ? 'bubbles' : 'agglomerates', U = c.vsl / c.holdup, TK = c.T + KEL;
  const n = Math.round(clamp(num(opts.nGroups, 12), 3, 40)), P = eulerPlan(c, opts, U, eulerCost.population(n)), dMin = clamp(num(opts.dMinUm, gas ? 500 : 20), 0.1, 1e5) * 1e-6, dMax = Math.max(clamp(num(opts.dMaxUm, gas ? 8000 : 2000), 0.2, 2e5) * 1e-6, 1.5 * dMin), dIn = clamp(num(opts.dInUm, gas ? 3000 : 100) * 1e-6, dMin, dMax);
  const eff = clamp(num(opts.collisionEfficiency, gas ? 1 : 0.01), 0, 100), dBreak = clamp(num(opts.dBreakUm, (dMax * 1e6) / 2) * 1e-6, dMin, dMax), kBreak = clamp(num(opts.breakupRate, 0.01), 0, 1e4), vol = (x) => (Math.PI * x ** 3) / 6, bExp = Math.min(1 / vol(dBreak), 20 / vol(dMax));
  const aD = clamp(num(opts.fraction, 0.05), 1e-4, 0.4), rhoP = clamp(num(opts.particleDensity, 920), 100, 20000), ratio = (dMax / dMin) ** (1 / (n - 1)), ds = Array.from({ length: n }, (_, i) => +ff(dMin * ratio ** i, 5));
  // inlet (and initial) distribution: log-normal weights around dIn, one size class wide, rounded so that they sum to exactly one
  const w = ds.map((di) => Math.exp(-0.5 * (Math.log(di / dIn) / Math.log(ratio)) ** 2)), ws = w.reduce((a, b) => a + b, 0), val = w.map((x) => Math.round((1e4 * x) / ws)), iMax = val.indexOf(Math.max(...val));
  val[iMax] += 1e4 - val.reduce((a, b) => a + b, 0);
  const groups = ds.map((di, i) => ({ name: 'f' + i, field: `f${i}.${d}.${pop}`, d: di, inlet: val[i] / 1e4 })), d32in = 1 / groups.reduce((s, g) => s + g.inlet / g.d, 0);
  const name = safeName(opts.name || `population_balance_${gas ? 'bubbles' : 'hydrate'}`), title = `Population balance (${n} size classes, ${gas ? 'gas bubbles' : 'hydrate particles'}) + Euler-Euler CFD, straight section at x = ${ff(c.x, 6)} m`, dl = `${d} in liquid`, cpL = pos(c.cpL, 2100), cpG = pos(c.cpG, 2400), mwEff = ((c.rhoG * R * TK) / (c.P * 1e5)) * 1000;
  const phaseProps = ['type            populationBalanceTwoPhaseSystem;', '', `phases          (${d} liquid);`, '', `populationBalances (${pop});`, '',
    ...eulerPhase(d, ['diameterModel   velocityGroup;', 'velocityGroupCoeffs', '{', `    populationBalance ${pop};`, '    formFactor      0.5235987756;', '    sizeGroups', '    (', ...groups.map((g) => `        ${g.name} { d ${ff(g.d, 5)}; value ${ff(g.inlet, 5)}; }`), '    );', '}'], ['residualAlpha   1e-6;']), ...eulerPhase('liquid', eulerConstD(1e-4), ['residualAlpha   1e-6;']),
    'populationBalanceCoeffs', '{', `    ${pop}`, '    {', '        continuousPhase liquid;', '', '        coalescenceModels', '        (', '            CoulaloglouTavlarides', '            {', `                C1              ${ff(2.8 * eff, 5)};`, '                C2              1.83e9;', '            }', '        );', '', '        binaryBreakupModels', '        (', '        );', '', '        breakupModels', '        (',
    ...(gas ? ['            LaakkonenAlopaeusAittamaa', '            {', '                daughterSizeDistributionModel LaakkonenAlopaeusAittamaa;', '            }'] : ['            exponential', '            {', `                C               ${ff(kBreak, 5)};`, `                exponent        ${ff(bExp, 6)};`, '                daughterSizeDistributionModel uniformBinary;', '            }']), '        );', '', '        driftModels', '        (', ...(gas ? ['            densityChange', '            {', '            }'] : []), '        );', '', '        nucleationModels', '        (', '        );', '    }', '}', '',
    'blending', '{', '    default', '    {', '        type            none;', '        continuousPhase liquid;', '    }', '}', '',
    ...eulerList('surfaceTension', [eulerPair(`${d} and liquid`, ['type            constant;', `sigma           ${ff(c.sigma, 5)};`])]), ...eulerList('aspectRatio', []),
    ...eulerList('drag', [eulerPair(dl, eulerDrag(gas ? 'IshiiZuber' : 'SchillerNaumann'))]), ...eulerList('virtualMass', [eulerPair(dl, ['type            constantCoefficient;', 'Cvm             0.5;'])]),
    ...eulerList('heatTransfer', [eulerPair(dl, ['type            RanzMarshall;', 'residualAlpha   1e-4;'])]), ...eulerList('phaseTransfer', []), ...eulerList('lift', []), ...eulerList('wallLubrication', []),
    ...eulerList('turbulentDispersion', [eulerPair(dl, ['type            Burns;', 'sigma           0.7;', 'Ctd             1;', 'residualAlpha   1e-3;'])]), '// Minimum allowable pressure', 'pMin            10000;'];
  const first = 'alpha.liquid', dF = 'd.' + d;
  const files = eulerTwoFluidFiles({ title, c, P, d, aD, U, TK, phaseProps, sizeGroups: groups, popName: pop,
    turb: { [d]: { model: gas ? 'continuousGasKEpsilon' : 'laminar', ke: gas }, liquid: { model: 'kEpsilon', ke: true } },
    thermo: { [d]: gas ? eulerThermo('perfectGas', mwEff, cpG, c.muG, (cpG * c.muG) / pos(c.kG, 0.04)) : ['// particle phase: constant density; its shear viscosity is set to that of the carrier liquid (dilute slurry)', ...eulerThermo('rhoConst', 100, 2100, c.muL, 1, rhoP)], liquid: eulerThermo('rhoConst', 100, cpL, c.muL, (cpL * c.muL) / pos(c.kL, 0.13), c.rhoL) },
    functions: eulerFunctions({ mesh: P.mesh, sections: P.sections, every: P.every, first, more: ['alpha.' + d, dF], outletMore: groups.map((g) => g.field) }) });
  const closures = [['Population balance', `class (size-group) method, \`populationBalanceTwoPhaseSystem\`, population \`${pop}\`, velocity group \`${d}\` with ${n} size groups from ${ff(dMin * 1e6, 4)} to ${ff(dMax * 1e6, 4)} µm (geometric, diameter ratio ${ff(ratio, 4)})`], ['Inlet size distribution', `log-normal around ${ff(dIn * 1e6, 4)} µm, one class wide; Sauter mean diameter ${ff(d32in * 1e6, 4)} µm; ${ff(100 * aD, 3)} % by volume`],
    [gas ? 'Coalescence' : 'Coalescence (agglomeration)', `Coulaloglou–Tavlarides (\`CoulaloglouTavlarides\`): turbulent collision frequency ∝ ε^(1/3) with a film-drainage (contact-time) efficiency; C1 = ${ff(2.8 * eff, 5)}, i.e. 2.8 × a collision efficiency of ${ff(eff, 4)} (option collisionEfficiency), C2 = 1.83e9 1/m²`], ['Breakup', gas ? 'Laakkonen–Alopaeus–Aittamaa (`LaakkonenAlopaeusAittamaa`) with its own daughter-size distribution (`daughterSizeDistributionModel LaakkonenAlopaeusAittamaa`)' : `exponential in the agglomerate volume (\`exponential\`): rate = ${ff(kBreak, 4)} 1/s × exp(v / v_c) with v_c the volume of a ${ff(dBreak * 1e6, 4)} µm agglomerate — negligible below that size, fast above it (options breakupRate, dBreakUm); binary daughters of uniform size probability (\`uniformBinary\`)`], gas ? ['Drift', 'bubble growth with falling pressure (`densityChange`)'] : null,
    ['Interfacial drag', gas ? 'Ishii–Zuber (`IshiiZuber`)' : 'Schiller–Naumann (`SchillerNaumann`) on the local Sauter mean diameter'], ['Virtual mass', '`constantCoefficient`, Cvm = 0.5'], ['Turbulent dispersion', 'Burns (`Burns`, σ = 0.7)'], ['Turbulence', gas ? 'liquid: standard k–ε (`kEpsilon`); gas: `continuousGasKEpsilon`' : 'liquid: standard k–ε (`kEpsilon`); particle phase: laminar stress with the viscosity of the carrier'], [gas ? 'Gas' : 'Particles', gas ? `perfect gas, ${ff(c.rhoG, 5)} kg/m³ at the case pressure` : `${ff(rhoP, 5)} kg/m³`]];
  const v = { name, title, files, P, U, solver: 'reactingTwoPhaseEulerFoam', label: `population balance (${n} size classes, ${gas ? 'bubbles' : 'hydrate particles'})`, first, invert: false, closures,
    extraRef: { dispersedPhase: d, fraction: aD, diameterField: dF, sauterInlet: d32in, groups: groups.map((g) => ({ name: g.name, field: g.field, d: g.d, inlet: g.inlet })) },
    inlet: `liquid and ${gas ? 'bubbles' : 'particles'} at the liquid velocity of the 1-D model, vsl/holdup = ${ff(U, 4)} m/s, dispersed fraction ${ff(aD, 4)}, size distribution as below`, more: `alpha.${d}, ${dF}`,
    intro: [gas ? 'The one-dimensional models of the application use one bubble size (or none) for the dispersed gas.' : 'The hydrate models of the application follow one mean particle or agglomerate size along the line.', `This case carries the whole size distribution: the ${gas ? 'bubbles' : 'particles'} are split into ${n} size classes that exchange volume by ${gas ? 'coalescence' : 'agglomeration'} and breakup at rates set by the local`, 'turbulence, while a three-dimensional two-fluid flow transports them — so the distribution at the outlet, and where the large agglomerates collect in the section, are results.'],
    notes: [gas ? null : 'The coalescence and breakup kernels of this OpenFOAM release were derived for fluid particles. For hydrate agglomerates they give the right dependences on dissipation rate and size, but their constants (cohesion, agglomerate strength) must be calibrated: edit `populationBalanceCoeffs` in constant/phaseProperties (`C1`, `C2` of CoulaloglouTavlarides, `C` and `exponent` of exponential) or use the options collisionEfficiency, breakupRate and dBreakUm.', 'Volume that coalesces beyond the largest size class leaves the distribution: keep the last class empty (the log prints `sizeGroup phase fraction first, last` and `sizeGroups-sum volume fraction` every step) by raising `dMaxUm` or the breakup rate.', 'The inlet and the initial field hold the same distribution (the `value` of each size group in constant/phaseProperties); the boundary types of all size groups come from `0.orig/f.' + d + '.' + pop + '`.', 'A residence time of a few seconds changes the distribution only if the kernels are fast; lengthen the section (`lengthD`) to follow the evolution along the line.'] };
  v.results = [...eulerResultLines(first, `\`alpha.${d}\`, \`${dF}\` (Sauter mean diameter, m)`, 'the liquid fraction').slice(0, 2), `- \`postProcessing/outletHoldup/0/surfaceFieldValue.dat\` — outlet averages of \`${first}\`, \`alpha.${d}\`, \`${dF}\` and of the size-group fractions \`${groups[0].field}\` … \`${groups[n - 1].field}\` (fraction of the dispersed volume in each class; class diameters in constant/phaseProperties and in the table below);`, `- \`postProcessing/probes/0/p\`, \`p_rgh\`, \`${first}\` — point values near the centreline at the same sections.`, '', '| Size group | Diameter (µm) | Inlet fraction |', '|---|---|---|', ...groups.map((g) => `| ${g.field} | ${ff(g.d * 1e6, 5)} | ${ff(g.inlet, 4)} |`)];
  return eulerAssemble(v, c);
}

// ---- reference equations of state: property scripts and property-table exchange ---------------------------------
export const TABLE_FIELDS = Object.freeze(['wG', 'rhoG', 'rhoO', 'muG', 'muO', 'cpG', 'cpO', 'kG', 'kO', 'hG', 'hO', 'jtG', 'jtO', 'sigma', 'zG', 'mwG', 'mwO']);
// CoolProp fluid names of the defined components; pseudo-components are mapped to the nearest n-alkane the mixture model covers
const HEOS_NAMES = Object.freeze({ N2: 'Nitrogen', CO2: 'CarbonDioxide', H2S: 'HydrogenSulfide', C1: 'Methane', C2: 'Ethane', C3: 'Propane', iC4: 'IsoButane', nC4: 'n-Butane', iC5: 'Isopentane', nC5: 'n-Pentane', C6: 'n-Hexane' });
const HEOS_ALKANES = Object.freeze([['n-Heptane', 100.202], ['n-Octane', 114.229], ['n-Nonane', 128.255], ['n-Decane', 142.282]]);
// PC-SAFT pure-component parameters m, σ (Å), ε/k (K): Gross & Sadowski, Ind. Eng. Chem. Res. 40 (2001) 1244; H2S from Tang & Gross, Fluid Phase Equilib. 293 (2010) 11
const PCSAFT = Object.freeze({ N2: [1.2053, 3.313, 90.96], CO2: [2.0729, 2.7852, 169.21], H2S: [1.6686, 3.0349, 229.0], C1: [1.0, 3.7039, 150.03], C2: [1.6069, 3.5206, 191.42], C3: [2.002, 3.6184, 208.11], iC4: [2.2616, 3.7574, 216.53], nC4: [2.3316, 3.7086, 222.88], iC5: [2.562, 3.8296, 230.75], nC5: [2.6896, 3.7729, 231.2], C6: [3.0576, 3.7983, 236.77] });
/** PC-SAFT parameters of a pseudo-component from its molar mass (n-alkane correlations of Tihic et al., Fluid Phase Equilib. 248 (2006) 29). */
export function pcsaftPseudo(MW) { const m = 0.02537 * MW + 0.9081; return [m, Math.cbrt((1.7284 * MW + 18.787) / m), (6.8311 * MW + 124.42) / m]; }

/**
 * Equilibrium compositions of the case fluid on a pressure–temperature grid (the grid of the application's property table),
 * for the property scripts. Returns { P[], T[], comps: [{ id, name, MW, cp[4], pseudo }], points: [{ i, j, beta, x[], y[] }], eosId }.
 */
export function propertyGrid(spec = DEFAULT_FLUID, { P, T } = {}) {
  const f = makeFluid({ ...DEFAULT_FLUID, ...spec, comp: { ...DEFAULT_FLUID.comp, ...(spec.comp || {}) } }), Ps = Array.isArray(P) && P.length ? P : linspace(Math.log10(1), Math.log10(600), 22).map((e) => 10 ** e), Ts = Array.isArray(T) && T.length ? T : linspace(-30, 170, 17), points = [];
  Ps.forEach((p, i) => Ts.forEach((tc, j) => { const r = flashPT(f, p, tc); points.push({ i, j, beta: r.phase === 'two' ? r.beta : r.phase === 'gas' ? 1 : 0, x: r.x.slice(), y: r.y.slice() }); }));
  return { P: Ps.slice(), T: Ts.slice(), comps: f.comps.map((c) => ({ id: c.id, name: c.name, MW: c.MW, cp: c.cp.slice(), pseudo: !!c.pseudo })), points, eosId: f.eosId };
}
const DEFAULT_GRID = () => propertyGrid(DEFAULT_FLUID, { P: [5, 20, 60, 150, 300], T: [0, 40, 80, 120] });
const pyList = (a, sig = 7) => '[' + a.map((v) => ff(v, sig)).join(', ') + ']';
const pyHeader = (title, lines) => ['#!/usr/bin/env python3', '"""', `HydraSlug hand-off: ${title}`, '', ...lines, '"""'];
const pyGrid = (g) => ['P_BAR = ' + pyList(g.P), 'T_C = ' + pyList(g.T), '# one row per grid point: (index of P, index of T, vapour mole fraction, liquid composition x, vapour composition y)', 'POINTS = [', ...g.points.map((p) => `    (${p.i}, ${p.j}, ${ff(p.beta, 6)}, ${pyList(p.x, 6)}, ${pyList(p.y, 6)}),`), ']'];

/**
 * Python script that evaluates the gas phase of the case with the GERG-2008 multi-fluid Helmholtz model of CoolProp on the grid of the
 * application's property table and writes `gerg2008_table.csv` (columns P, T, wG, rhoG, zG, cpG, jtG, wSound, source) for re-import.
 * opts: { grid (from propertyGrid), fluidName, liquid (bool: also evaluate the liquid phase — only meaningful for light condensates) }.
 */
export function coolpropScript(opts = {}) {
  const g = opts.grid && opts.grid.points ? opts.grid : DEFAULT_GRID(), heavy = [];
  const names = g.comps.map((c) => { if (HEOS_NAMES[c.id]) return HEOS_NAMES[c.id]; let best = HEOS_ALKANES[0]; for (const a of HEOS_ALKANES) if (Math.abs(a[1] - c.MW) < Math.abs(best[1] - c.MW)) best = a; heavy.push(`${c.name} (M = ${ff(c.MW, 4)} g/mol) → ${best[0]}`); return best[0]; });
  const out = 'gerg2008_table.csv', name = 'gerg2008_coolprop';
  const script = text([pyHeader('GERG-2008 (multi-fluid Helmholtz) properties with CoolProp', ['Evaluates the gas phase of the case fluid at the equilibrium vapour compositions computed by the application,', 'on the pressure-temperature grid of its property table, and writes ' + out + '.', 'Run:  pip install CoolProp ; python ' + name + '.py', 'Then import ' + out + ' on the External solvers page.']),
    'import csv', 'import sys', '', 'import CoolProp', 'import CoolProp.CoolProp as CP', '',
    'NAMES = ' + JSON.stringify(names), ...pyGrid(g), 'LIQUID = ' + (opts.liquid ? 'True' : 'False'), '', '',
    'def merged(z):', '    """Mole fractions per distinct CoolProp fluid (pseudo-components mapped to the same alkane are added)."""', '    tot = {}', '    for n, v in zip(NAMES, z):', '        tot[n] = tot.get(n, 0.0) + max(v, 0.0)', '    keys = [k for k in tot if tot[k] > 1e-9]', '    s = sum(tot[k] for k in keys)', '    return keys, [tot[k] / s for k in keys]', '', '',
    'def evaluate(z, p_pa, t_k, phase):', '    keys, frac = merged(z)', '    st = CP.AbstractState("HEOS", "&".join(keys))', '    st.set_mole_fractions(frac)', '    st.specify_phase(phase)', '    st.update(CP.PT_INPUTS, p_pa, t_k)', '    jt = st.first_partial_deriv(CP.iT, CP.iP, CP.iHmass)', '    return st.rhomass(), st.compressibility_factor(), st.cpmass(), jt, st.speed_sound()', '', '',
    'def main():', '    rows, failed = [], 0', '    for (i, j, beta, x, y) in POINTS:', '        p, t = P_BAR[i], T_C[j]', '        row = {"P": p, "T": t, "wG": "", "rhoG": "", "zG": "", "cpG": "", "jtG": "", "wSound": "", "rhoO": "", "cpO": "", "source": "GERG-2008 (CoolProp %s)" % CoolProp.__version__}', '        if beta > 1e-9:', '            try:', '                rho, zf, cp, jt, w = evaluate(y, p * 1e5, t + 273.15, CP.iphase_gas)', '                row.update(rhoG="%.8g" % rho, zG="%.8g" % zf, cpG="%.8g" % cp, jtG="%.8g" % jt, wSound="%.8g" % w)', '            except Exception as err:  # the reference model has no root here (far inside the two-phase region of the mapped mixture)', '                failed += 1', '                row["source"] = "failed: %s" % str(err)[:80].replace(",", ";")', '        if LIQUID and beta < 1 - 1e-9:', '            try:', '                rho, zf, cp, jt, w = evaluate(x, p * 1e5, t + 273.15, CP.iphase_liquid)', '                row.update(rhoO="%.8g" % rho, cpO="%.8g" % cp)', '            except Exception:', '                pass', '        rows.append(row)', `    with open("${out}", "w", newline="") as fh:`, '        wr = csv.DictWriter(fh, fieldnames=["P", "T", "wG", "rhoG", "zG", "cpG", "jtG", "wSound", "rhoO", "cpO", "source"])', '        wr.writeheader()', '        wr.writerows(rows)', `    print("wrote ${out}: %d grid points, %d without a gas-phase result" % (len(rows), failed))`, '    return 0', '', '', 'if __name__ == "__main__":', '    sys.exit(main())']);
  const readme = text([`# ${name} — CoolProp hand-off`, '', 'The application computes phase behaviour and properties with a cubic equation of state. GERG-2008 is the reference multi-fluid Helmholtz model for natural-gas mixtures', '(uncertainty of gas density about 0.1 %); it needs 21 pure-fluid reference equations and 210 binary departure functions and is evaluated here with CoolProp.', '',
    `Grid: ${g.P.length} pressures (${ff(g.P[0], 4)}–${ff(g.P[g.P.length - 1], 4)} bara) × ${g.T.length} temperatures (${ff(g.T[0], 4)}–${ff(g.T[g.T.length - 1], 4)} °C). At each point the script evaluates the gas phase at the vapour composition`, `found by the application's ${g.eosId} flash, so the comparison isolates the property model from the phase split.`, '',
    heavy.length ? ['GERG-2008 has no heptanes-plus fraction. The pseudo-components are mapped to the nearest n-alkane of the model:', '', ...heavy.map((h) => `- ${h}`), '', 'This is accurate for the gas phase, where the heavy end is a few tenths of a per cent; it is not a model of the oil phase, so liquid columns are left empty unless you set `LIQUID = True`.', ''] : null,
    '## Run', '', '```sh', 'python3 -m venv venv && . venv/bin/activate', 'pip install CoolProp', `python ${name}.py        # writes ${out}`, '```', '', '## Results', '', `\`${out}\` has one row per grid point: P (bara), T (°C), rhoG (kg/m³), zG, cpG (J/kg/K), jtG (K/Pa), wSound (m/s). Import it on the External solvers page:`, 'the application shows the deviation of its own table and can publish the GERG values as the gas properties of the case property table.', '']);
  return { name, files: [{ path: `${name}.py`, text: script }, { path: 'README.md', text: readme }], readme, summary: `CoolProp GERG-2008 script: gas-phase density, Z, cp, Joule–Thomson coefficient and speed of sound at ${g.points.length} grid points.`, commands: ['pip install CoolProp', `python ${name}.py`], output: out, plan: { cells: g.points.length, coreHours: 0.001 }, reference: {} };
}

/**
 * Python script that evaluates both phases of the case with PC-SAFT (teqp) at the equilibrium compositions of the application,
 * on the grid of its property table, and writes `pcsaft_table.csv` (P, T, rhoG, zG, cpG, jtG, rhoO, cpO, jtO) for re-import.
 * With model 'CPA' it also evaluates liquid water with the cubic-plus-association equation (column rhoAq).
 * opts: { grid (from propertyGrid), model: 'PCSAFT' | 'CPA', aqueous: { inhibitor: 'MeOH' | 'MEG' | 'none', wt } }.
 */
export function teqpScript(opts = {}) {
  const g = opts.grid && opts.grid.points ? opts.grid : DEFAULT_GRID(), cpa = opts.model === 'CPA', out = cpa ? 'cpa_table.csv' : 'pcsaft_table.csv', name = cpa ? 'cpa_teqp' : 'pcsaft_teqp';
  const coeffs = g.comps.map((c) => { const p = PCSAFT[c.id] || pcsaftPseudo(c.MW); return `    {"name": ${JSON.stringify(c.name)}, "m": ${ff(p[0], 6)}, "sigma_Angstrom": ${ff(p[1], 6)}, "epsilon_over_k": ${ff(p[2], 6)}, "BibTeXKey": "${PCSAFT[c.id] ? 'Gross-IECR-2001' : 'Tihic-FPE-2006'}"},`; });
  const aq = opts.aqueous || {}, inh = aq.inhibitor === 'MeOH' || aq.inhibitor === 'MEG' ? aq.inhibitor : 'none', wt = clamp(num(aq.wt, inh === 'none' ? 0 : 30), 0, 95);
  // CPA parameters (SRK + association): Kontogeorgis & Folas, Thermodynamic Models for Industrial Applications (2010), tables 9.1–9.3
  const CPA = { water: [0.12277, 1.4515e-5, 0.6736, 647.096, 16655, 0.0692, '4C', 18.015] };
  const cpaPure = (k) => `{"a0i / Pa m^6/mol^2": ${CPA[k][0]}, "bi / m^3/mol": ${CPA[k][1]}, "c1": ${CPA[k][2]}, "Tc / K": ${CPA[k][3]}, "epsABi / J/mol": ${CPA[k][4]}, "betaABi": ${CPA[k][5]}, "class": "${CPA[k][6]}"}`;
  const script = text([pyHeader(cpa ? 'cubic-plus-association (CPA) and PC-SAFT properties with teqp' : 'PC-SAFT properties with teqp', ['Evaluates gas and liquid of the case fluid with PC-SAFT at the equilibrium compositions computed by the application,', 'on the pressure-temperature grid of its property table, and writes ' + out + '.', cpa ? 'The aqueous phase (water and hydrate inhibitor) is evaluated with CPA: SRK plus Wertheim association.' : null, 'Run:  pip install teqp numpy ; python ' + name + '.py', 'Then import ' + out + ' on the External solvers page.'].filter(Boolean)),
    'import csv', 'import math', 'import sys', '', 'import numpy as np', 'import teqp', '', 'R = 8.314462618', '# PC-SAFT pure-component parameters (no binary interaction parameters: k_ij = 0)', 'COEFFS = [', ...coeffs, ']', 'MW = ' + pyList(g.comps.map((c) => c.MW)) + '  # g/mol', '# ideal-gas heat capacity cp0 = a + b T + c T^2 + d T^3 (J/mol/K), the same polynomials as the application', 'CP0 = [' + g.comps.map((c) => pyList(c.cp, 6)).join(', ') + ']', ...pyGrid(g),
    cpa ? ['# CPA parameters of water (4C association scheme)', `CPA_WATER = ${cpaPure('water')}`, `MW_AQ = [${CPA.water[7]}]`] : null, '', '',
    'def pressure(model, t, rho, z):', '    return rho * R * t * (1.0 + model.get_Ar01(t, rho, z))', '', '',
    'def density(model, t, p, z, liquid, rho_max):', '    """Molar density (mol/m3) of the gas-like or liquid-like root of p(rho) = p at fixed T and composition, or None."""', '    grid = np.geomspace(1e-2, rho_max, 240)', '    f = np.array([pressure(model, t, r, z) - p for r in grid])', '    roots = []', '    for k in range(len(grid) - 1):', '        if not (math.isfinite(f[k]) and math.isfinite(f[k + 1])) or f[k] * f[k + 1] > 0 or f[k + 1] < f[k]:', '            continue  # no sign change, or a mechanically unstable branch', '        a, b = grid[k], grid[k + 1]', '        for _ in range(80):', '            m = 0.5 * (a + b)', '            if (pressure(model, t, m, z) - p) * (pressure(model, t, a, z) - p) <= 0:', '                b = m', '            else:', '                a = m', '        roots.append(0.5 * (a + b))', '    if not roots:', '        return None', '    return roots[-1] if liquid else roots[0]', '', '',
    'def props(model, t, p, z, liquid, mw, cp0_poly, rho_max):', '    z = np.array(z, dtype=float)', '    z = np.maximum(z, 1e-12)', '    z = z / z.sum()', '    rho = density(model, t, p, z, liquid, rho_max)', '    if rho is None:', '        return None', '    m = float(np.dot(z, mw)) * 1e-3  # kg/mol', '    a01, a02 = model.get_Ar01(t, rho, z), model.get_Ar02(t, rho, z)', '    a11, a20 = model.get_Ar11(t, rho, z), model.get_Ar20(t, rho, z)', '    cp0 = sum(zi * (c[0] + c[1] * t + c[2] * t * t + c[3] * t ** 3) for zi, c in zip(z, cp0_poly))', '    cv = cp0 - R - R * a20', '    dpdrho = R * t * (1.0 + 2.0 * a01 + a02)', '    dpdt = rho * R * (1.0 + a01 - a11)', '    cp = cv + t * dpdt * dpdt / (rho * rho * dpdrho)', '    jt = (t * dpdt / (rho * dpdrho) - 1.0) / (rho * cp)  # K/Pa', '    return {"rho": rho * m, "Z": p / (rho * R * t), "cp": cp / m, "jt": jt}', '', '',
    'def main():', '    model = teqp.make_model({"kind": "PCSAFT", "model": {"coeffs": COEFFS}})', '    # densest packing of the hard-chain reference: eta = (pi/6) rho N_A sum(x m d^3) < 0.74', '    seg = [c["m"] * (c["sigma_Angstrom"] * 1e-10) ** 3 for c in COEFFS]', '    rows, failed = [], 0',
    cpa ? ['    aq_model, aq_z, aq_mw = None, np.array([1.0]), [MW_AQ[0]]', '    try:', '        aq_model = teqp.make_model({"kind": "CPA", "model": {"cubic": "SRK", "radial_dist": "KG", "pures": [CPA_WATER], "R_gas / J/mol/K": R}})', '    except Exception as err:', '        print("CPA model could not be built by this teqp version: %s" % err)'] : null,
    '    for (i, j, beta, x, y) in POINTS:', '        p, t = P_BAR[i] * 1e5, T_C[j] + 273.15', `        row = {"P": P_BAR[i], "T": T_C[j], "rhoG": "", "zG": "", "cpG": "", "jtG": "", "rhoO": "", "cpO": "", "jtO": "", ${cpa ? '"rhoAq": "", ' : ''}"source": "PC-SAFT (teqp %s)" % teqp.__version__}`, '        for liquid, z in ((False, y), (True, x)):', '            if (beta < 1e-9 and not liquid) or (beta > 1 - 1e-9 and liquid):', '                continue', '            zz = np.maximum(np.array(z, dtype=float), 1e-12)', '            zz = zz / zz.sum()', '            rho_max = 0.74 / (math.pi / 6.0 * 6.02214076e23 * float(np.dot(zz, seg)))', '            try:', '                r = props(model, t, p, z, liquid, MW, CP0, rho_max)', '            except Exception:', '                r = None', '            if r is None:', '                failed += 1', '                continue', '            if liquid:', '                row.update(rhoO="%.8g" % r["rho"], cpO="%.8g" % r["cp"], jtO="%.8g" % r["jt"])', '            else:', '                row.update(rhoG="%.8g" % r["rho"], zG="%.8g" % r["Z"], cpG="%.8g" % r["cp"], jtG="%.8g" % r["jt"])',
    cpa ? ['        if aq_model is not None:', '            try:', '                rho_aq = density(aq_model, t, p, aq_z, True, 0.98 / CPA_WATER["bi / m^3/mol"])', '                if rho_aq is not None:', '                    row["rhoAq"] = "%.8g" % (rho_aq * aq_mw[0] * 1e-3)', '            except Exception:', '                pass'] : null,
    '        rows.append(row)', `    with open("${out}", "w", newline="") as fh:`, `        wr = csv.DictWriter(fh, fieldnames=["P", "T", "rhoG", "zG", "cpG", "jtG", "rhoO", "cpO", "jtO", ${cpa ? '"rhoAq", ' : ''}"source"])`, '        wr.writeheader()', '        wr.writerows(rows)', `    print("wrote ${out}: %d grid points, %d phase evaluations without a root" % (len(rows), failed))`, '    return 0', '', '', 'if __name__ == "__main__":', '    sys.exit(main())']);
  const readme = text([`# ${name} — teqp hand-off`, '', `PC-SAFT (perturbed-chain statistical associating fluid theory) describes chain molecules from segment number, segment diameter and dispersion energy; it is more accurate than a cubic equation for liquid`, 'densities and derivative properties of heavy hydrocarbons. The script evaluates it with teqp at the equilibrium compositions of the application, so the comparison isolates the property model from the phase split.', '',
    'Parameters: Gross & Sadowski (2001) for the defined components, Tang & Gross (2010) for hydrogen sulphide, the n-alkane correlations of Tihic et al. (2006) with the molar mass of each pseudo-component;', 'binary interaction parameters are zero. The ideal-gas heat capacities are the polynomials of the application.', cpa ? ['', 'CPA (SRK + Wertheim association, 4C scheme, parameters of Kontogeorgis & Folas 2010) is evaluated for liquid water, the associating component of the produced fluids: column rhoAq.', inh !== 'none' ? `The case carries ${ff(wt, 4)} wt % ${inh}: cross-association of water with the inhibitor needs the association-site input of teqp and is not written by this script.` : null] : null, 'The n-alkane correlations underestimate the density of aromatic and naphthenic heavy ends; treat the oil columns as a cross-check, not as a replacement of a tuned model.', '',
    '## Run', '', '```sh', 'python3 -m venv venv && . venv/bin/activate', 'pip install teqp numpy', `python ${name}.py        # writes ${out}`, '```', '', '## Results', '', `\`${out}\`: P (bara), T (°C), rhoG, zG, cpG, jtG, rhoO, cpO, jtO (SI units, per kg)${cpa ? ', rhoAq' : ''}. Import it on the External solvers page to see the deviation of the cubic-EOS table and to publish the values.`, '']);
  return { name, files: [{ path: `${name}.py`, text: script }, { path: 'README.md', text: readme }], readme, summary: `teqp ${cpa ? 'CPA + ' : ''}PC-SAFT script: densities, heat capacities and Joule–Thomson coefficients of both phases at ${g.points.length} grid points.`, commands: ['pip install teqp numpy', `python ${name}.py`], output: out, plan: { cells: g.points.length, coreHours: 0.01 }, reference: {} };
}

/** CSV of a property table (the shape of buildTable): one row per grid point with every tabulated field. */
export function propertyTableCSV(table) {
  const rows = ['P,T,' + TABLE_FIELDS.join(',')];
  table.P.forEach((p, i) => table.T.forEach((t, j) => rows.push([p, t, ...TABLE_FIELDS.map((k) => table[k][i][j])].map((v) => (isNum(v) ? String(v) : '')).join(','))));
  return rows.join('\n') + '\n';
}
/** Minimal RFC-4180 reader for importers: returns { headers, rows: [[cell strings]] } with caps on size. */
function csvRows(src, maxRows = BRIDGE_LIMITS.rows) {
  if (typeof src !== 'string') throw new Error('Expected text.');
  if (src.length > BRIDGE_LIMITS.chars) throw new Error('File is too large to import.');
  const lines = src.replace(/^﻿/, '').split(/\r?\n/).filter((l) => l.trim() !== '' && !l.startsWith('#'));
  if (!lines.length) throw new Error('The file is empty.');
  if (lines.length > maxRows) throw new Error('Too many rows.');
  const delim = [',', ';', '\t'].map((d) => [d, lines[0].split(d).length]).sort((a, b) => b[1] - a[1])[0][0];
  const split = (l) => { const out = []; let cell = '', q = false; for (let i = 0; i < l.length; i++) { const ch = l[i]; if (q) { if (ch === '"') { if (l[i + 1] === '"') { cell += '"'; i++; } else q = false; } else cell += ch; } else if (ch === '"') q = true; else if (ch === delim) { out.push(cell); cell = ''; } else cell += ch; } out.push(cell); return out.map((s) => s.trim()); };
  return { headers: split(lines[0]), rows: lines.slice(1).map(split) };
}
const cellNum = (s) => { if (s === undefined || s === null || s === '') return null; const v = Number(s); return Number.isFinite(v) ? v : null; };
/**
 * Read a property table from CSV (columns P in bara, T in °C, then any of the tabulated fields) into the exact shape of
 * `buildTable` in thermo.js, so that it can be published as outputs.pvt.table. Fields or cells that the file does not give
 * are taken from `base` (the application's own table on the same grid); without a base every field must be present.
 * Returns { table, replaced: [fields], filled (cells taken from the file), missing (cells kept from the base), extra: { name: [iP][iT] } , source }.
 */
export function importPropertyTable(csv, base = null) {
  const { headers, rows } = csvRows(csv), col = (n) => headers.findIndex((h) => h.toLowerCase().replace(/[^a-z0-9]/g, '') === n.toLowerCase());
  const iP = [col('P'), col('Pbar'), col('Pbara')].find((k) => k >= 0), iT = [col('T'), col('TC'), col('TdegC')].find((k) => k >= 0);
  if (iP === undefined || iT === undefined) throw new Error('The table needs columns P (bara) and T (°C).');
  const uniq = (k) => [...new Set(rows.map((r) => cellNum(r[k])).filter((v) => v !== null))].sort((a, b) => a - b), P = uniq(iP), T = uniq(iT);
  if (P.length < 2 || T.length < 2) throw new Error('The table needs at least two pressures and two temperatures.');
  if (P.length * T.length > 40000 || P[0] <= 0) throw new Error('Unreasonable pressure–temperature grid.');
  const near = (arr, v) => { let b = 0; for (let k = 1; k < arr.length; k++) if (Math.abs(arr[k] - v) < Math.abs(arr[b] - v)) b = k; return b; };
  let useBase = null;
  if (base && Array.isArray(base.P) && Array.isArray(base.T)) { // the same grid within rounding of the CSV
    const same = base.P.length === P.length && base.T.length === T.length && base.P.every((p, k) => Math.abs(p - P[k]) <= 1e-4 * Math.max(1, Math.abs(p))) && base.T.every((t, k) => Math.abs(t - T[k]) <= 1e-4 * Math.max(1, Math.abs(t)));
    if (!same) throw new Error('The imported grid differs from the grid of the case property table; write the script again from the current case.');
    useBase = base;
  }
  const table = { P: useBase ? useBase.P.slice() : P, T: useBase ? useBase.T.slice() : T }, replaced = [], extra = {};
  table.lnP = table.P.map(Math.log); table.eosId = String(useBase?.eosId || 'imported').slice(0, 20);
  let filled = 0, missing = 0;
  const srcCol = col('source'), source = srcCol >= 0 ? String(rows.find((r) => r[srcCol] && !/^failed/.test(r[srcCol]))?.[srcCol] || '').slice(0, 80) : '';
  const grid = (k) => { const a = table.P.map(() => new Array(table.T.length).fill(null)); for (const r of rows) { const p = cellNum(r[iP]), t = cellNum(r[iT]), v = cellNum(r[k]); if (p === null || t === null || v === null) continue; a[near(P, p)][near(T, t)] = v; } return a; };
  for (const f of TABLE_FIELDS) {
    const k = col(f), a = k >= 0 ? grid(k) : null;
    if (!a && !useBase) throw new Error(`Column ${f} is missing and there is no case table to take it from.`);
    let n = 0;
    table[f] = table.P.map((_, i) => table.T.map((__, j) => { const v = a ? a[i][j] : null; if (v !== null) { n++; return v; } if (!useBase) throw new Error(`Column ${f} has an empty cell and there is no case table to take it from.`); missing++; return useBase[f][i][j]; }));
    if (n) { replaced.push(f); filled += n; }
  }
  // physical sanity of what will be interpolated by the flow solvers
  for (const f of ['rhoG', 'rhoO', 'muG', 'muO', 'cpG', 'cpO', 'zG']) if (table[f].some((r) => r.some((v) => !(v > 0)))) throw new Error(`Column ${f} must be positive everywhere.`);
  table.wG = table.wG.map((r) => r.map((v) => clamp(v, 0, 1)));
  headers.forEach((h, k) => { if (k !== iP && k !== iT && k !== srcCol && !TABLE_FIELDS.includes(h) && /^[A-Za-z][\w]{0,20}$/.test(h) && Object.keys(extra).length < 8) { const a = grid(k); if (a.some((r) => r.some((v) => v !== null))) extra[h] = a; } });
  if (useBase) { table.spec = JSON.parse(JSON.stringify(useBase.spec || {})); table.gasSG = useBase.gasSG; table.rates = JSON.parse(JSON.stringify(useBase.rates || {})); table.aq = JSON.parse(JSON.stringify(useBase.aq || {})); }
  else { table.spec = {}; table.gasSG = 0.7; table.rates = {}; table.aq = {}; }
  return { table, replaced, filled, missing, extra, source };
}
/** Deviation of an imported property table from the application's table, field by field: [{ field, n, meanPct, maxPct, bias }] over the cells the import supplied. */
export function comparePropertyTables(own, imported, fields = ['rhoG', 'zG', 'cpG', 'jtG', 'rhoO', 'cpO']) {
  const out = [];
  for (const f of fields) {
    if (!own[f] || !imported[f]) continue;
    let n = 0, s = 0, mx = 0, b = 0;
    own.P.forEach((_, i) => own.T.forEach((__, j) => { const a = own[f][i][j], c = imported[f][i]?.[j]; if (!isNum(a) || !isNum(c) || a === c || Math.abs(a) < 1e-300) return; if ((f.endsWith('G') && !(own.wG[i][j] > 0)) || (f.endsWith('O') && !(own.wG[i][j] < 1))) return; const d = (100 * (a - c)) / Math.abs(c); n++; s += Math.abs(d); b += d; mx = Math.max(mx, Math.abs(d)); }));
    if (n) out.push({ field: f, n, meanPct: s / n, maxPct: mx, bias: b / n });
  }
  return out;
}

// ---- importers: results of the external solvers --------------------------------------------------------------
const SEC_PER_YEAR = 31557600;
/** Split a text made by collect.sh ("### file: path" markers) into [{ name, text }]; a text without markers is one file. */
function splitCollected(src, name = '') {
  if (typeof src !== 'string') throw new Error('Expected text.');
  if (src.length > BRIDGE_LIMITS.chars) throw new Error('File is too large to import.');
  if (!/^### file: /m.test(src)) return [{ name, text: src }];
  const out = [], re = /^### file: (.*)$/gm; let m, last = null;
  while ((m = re.exec(src))) { if (last) out.push({ name: last.name, text: src.slice(last.end, m.index) }); last = { name: m[1].trim().slice(0, 300), end: re.lastIndex }; if (out.length > BRIDGE_LIMITS.files) throw new Error('Too many files in one import.'); }
  if (last) out.push({ name: last.name, text: src.slice(last.end) });
  return out;
}
/**
 * Parse the time-series tables written by OpenFOAM function objects (surfaceFieldValue, volFieldValue, forces, probes, cloud information)
 * or a plain CSV with a time column. src: text, or [{ name, text }]. Returns [{ file, object, columns: [names], t: [], data: [[...] per column] }].
 */
export function parseFoamTables(src) {
  const files = (Array.isArray(src) ? src.slice(0, BRIDGE_LIMITS.files).flatMap((f) => splitCollected(String(f.text ?? ''), String(f.name || ''))) : splitCollected(src)), out = [];
  let total = 0;
  for (const f of files) {
    total += f.text.length; if (total > BRIDGE_LIMITS.chars) throw new Error('Files are too large to import together.');
    const path = f.name.replace(/\\/g, '/'), parts = path.split('/').filter(Boolean), pp = parts.indexOf('postProcessing'), base = parts[parts.length - 1] || 'table', object = ((pp >= 0 && parts[pp + 1]) || base.replace(/\.(dat|csv|txt)$/i, '')) + (/^moment/i.test(base) ? '.moment' : '');
    const lines = f.text.split(/\r?\n/);
    let header = null, probes = 0; const rows = [];
    for (const raw of lines) {
      const l = raw.trim();
      if (!l) continue;
      if (l.startsWith('#')) { if (/^#\s*Probe\s+\d+/.test(l)) probes++; if (/^#\s*Time\b/i.test(l)) header = l.replace(/^#\s*/, ''); continue; }
      if (!header && !rows.length && /[A-Za-z]/.test(l.replace(/[eE][-+]?\d/g, '')) && /[,;\t]/.test(l)) { header = l; continue; } // CSV header line
      const cells = l.replace(/[()]/g, ' ').split(/[\s,;]+/).filter((c) => c !== '').map(Number);
      if (cells.length < 2 || !cells.every(Number.isFinite)) continue;
      rows.push(cells); if (rows.length > BRIDGE_LIMITS.rows) throw new Error('Too many rows in ' + (f.name || 'the file') + '.');
    }
    if (rows.length < 2) continue;
    const nCol = Math.min(BRIDGE_LIMITS.series + 1, Math.min(...rows.slice(0, 50).map((r) => r.length)));
    let names = header ? (/[,;]/.test(header) ? header.split(/[,;]/) : header.split(/\s+/)).map((s) => s.trim()).filter(Boolean) : [];
    if (probes) { const field = base, per = Math.max(1, Math.round((nCol - 1) / probes)); names = ['Time']; for (let p = 0; p < probes; p++) for (let k = 0; k < per; k++) names.push(`${field}.${p + 1}${per > 1 ? '.' + 'xyz'[k] : ''}`); }
    if (names.length !== nCol) names = ['Time', ...Array.from({ length: nCol - 1 }, (_, k) => names[k + 1] && names.length - 1 === nCol - 1 ? names[k + 1] : `col${k + 1}`)];
    const t = [], data = Array.from({ length: nCol - 1 }, () => []);
    let prev = -Infinity;
    for (const r of rows) { if (!(r[0] > prev)) continue; prev = r[0]; t.push(r[0]); for (let k = 1; k < nCol; k++) data[k - 1].push(r[k]); } // restarts repeat times: keep the first pass, strictly increasing
    out.push({ file: f.name, object: object.slice(0, 60), columns: names.slice(1).map((s) => String(s).slice(0, 80)), t, data });
  }
  if (!out.length) throw new Error('No time-series table was recognised. Choose hydraslug_results.txt or the .dat files under postProcessing.');
  return out;
}
/** Frequency of the dominant fluctuation of a signal: upward crossings of its mean with hysteresis, and the peak of its spectrum. */
export function fluctuation(t, y) {
  const n = Math.min(t.length, y.length);
  if (n < 8) return { mean: n ? mean(y.slice(0, n)) : 0, std: 0, min: n ? Math.min(...y.slice(0, n)) : 0, max: n ? Math.max(...y.slice(0, n)) : 0, crossings: 0, freqCrossing: 0, freqSpectrum: 0, amplitude: 0 };
  let mu = 0, lo = Infinity, hi = -Infinity; for (let i = 0; i < n; i++) { mu += y[i]; lo = Math.min(lo, y[i]); hi = Math.max(hi, y[i]); } mu /= n;
  let va = 0; for (let i = 0; i < n; i++) va += (y[i] - mu) ** 2; const sd = Math.sqrt(va / n), band = 0.5 * sd, dur = t[n - 1] - t[0];
  let state = y[0] > mu ? 1 : -1, up = 0, first = null, lastUp = null;
  for (let i = 1; i < n; i++) { if (state < 0 && y[i] > mu + band) { state = 1; up++; if (first === null) first = t[i]; lastUp = t[i]; } else if (state > 0 && y[i] < mu - band) state = -1; }
  const fCross = up >= 2 && lastUp > first ? (up - 1) / (lastUp - first) : 0;
  // spectrum of the signal resampled on a uniform grid (direct transform; at most 1024 points)
  const N = Math.min(1024, n), u = new Float64Array(N); for (let k = 0; k < N; k++) u[k] = interp1(t.slice(0, n), y.slice(0, n), t[0] + (dur * k) / (N - 1)) - mu;
  let best = 0, bestP = 0;
  for (let m = 1; m < N / 2; m++) { let re = 0, im = 0; const w = (2 * Math.PI * m) / N; for (let k = 0; k < N; k++) { const h = 0.5 - 0.5 * Math.cos((2 * Math.PI * k) / (N - 1)); re += h * u[k] * Math.cos(w * k); im -= h * u[k] * Math.sin(w * k); } const p = re * re + im * im; if (p > bestP) { bestP = p; best = m; } }
  const fSpec = sd > 1e-12 * Math.max(1, Math.abs(mu)) && dur > 0 && best >= 2 ? (best * (N - 1)) / (N * dur) : 0;
  return { mean: mu, std: sd, min: lo, max: hi, crossings: up, freqCrossing: fCross, freqSpectrum: fSpec, amplitude: 0.5 * (hi - lo) };
}
const tail = (t, y, frac) => { const t0 = t[0] + frac * (t[t.length - 1] - t[0]); let i = 0; while (i < t.length - 2 && t[i] < t0) i++; return { t: t.slice(i), y: y.slice(i) }; };
const slope = (t, y) => { const n = t.length; if (n < 2) return 0; const mt = mean(t), my = mean(y); let a = 0, b = 0; for (let i = 0; i < n; i++) { a += (t[i] - mt) * (y[i] - my); b += (t[i] - mt) ** 2; } return b > 0 ? a / b : 0; };

/**
 * Results of a hand-off OpenFOAM case. src: the text of hydraslug_results.txt (or of one function-object file), or [{ name, text }].
 * ref: the `reference` object of the generator that wrote the case (sectionDistance, wallArea, wallFaceArea, invert …); discard = fraction of the
 * simulated time left out of the averages as start-up (default 0.3).
 * Returns { t, series: [{ name, object, t, y }], metrics: { holdupMean, holdupSection, slugFrequency, slugFrequencySpectrum, holdupStd, pressureGradient,
 *   wallShear, forcePeak, forceMean, forceStd, pressureAmplitude, pressurePeriod, erosionRateMmY, erodedVolume, pressureDrop, sauterMean, duration, … (null when the files do not hold it) }, files }.
 */
export function importOpenfoamPostProcessing(src, ref = {}, { discard = 0.3 } = {}) {
  const tables = parseFoamTables(src), series = [], find = (obj, re) => series.find((s) => s.object === obj && (!re || re.test(s.column)));
  for (const tb of tables) tb.columns.forEach((cname, k) => { if (series.length < BRIDGE_LIMITS.series * 4) series.push({ name: `${tb.object}: ${cname}`, object: tb.object, column: cname, t: tb.t, y: tb.data[k] }); });
  const inv = (v) => (ref.invert ? 1 - v : v), frac = clamp(num(discard, 0.3), 0, 0.9), M = {}, stat = (s) => { if (!s) return null; const w = tail(s.t, s.y, frac); return fluctuation(w.t, w.y); };
  const alphaRe = /alpha\.|alpha$/i, hv = find('holdupVolume', alphaRe) || find('holdupVolume'), s1 = find('section1', alphaRe), s2 = find('section2', alphaRe), s3 = find('section3', alphaRe), ho = find('outletHoldup', alphaRe);
  const p1 = find('section1', /\(p\)|^p$/), p2 = find('section2', /\(p\)|^p$/), p3 = find('section3', /\(p\)|^p$/), primary = hv || s2 || series[0];
  M.duration = primary.t[primary.t.length - 1] - primary.t[0]; M.samples = primary.t.length;
  const a = stat(hv), b = stat(s2), cst = stat(s3 || ho || s2);
  M.holdupMean = a ? inv(a.mean) : b ? inv(b.mean) : null; M.holdupSection = b ? inv(b.mean) : null; M.holdupStd = cst ? cst.std : null; M.holdupMin = cst ? Math.min(inv(cst.min), inv(cst.max)) : null; M.holdupMax = cst ? Math.max(inv(cst.min), inv(cst.max)) : null;
  // slugs are counted only when the holdup really swings (a tenth of the full range)
  M.slugFrequency = cst && cst.max - cst.min > 0.1 ? cst.freqCrossing : cst ? 0 : null; M.slugFrequencySpectrum = cst && cst.max - cst.min > 0.1 ? cst.freqSpectrum : cst ? 0 : null;
  const q1 = stat(p1), q3 = stat(p3), q2 = stat(p2), dist = pos(ref.sectionDistance, 0);
  M.pressureGradient = q1 && q3 && dist > 0 ? (q1.mean - q3.mean) / dist : null;
  M.pressureAmplitude = q2 ? q2.amplitude : null; M.pressurePeriod = q2 && q2.freqCrossing > 0 ? 1 / q2.freqCrossing : null; M.pressureMean = q2 ? q2.mean : null;
  const fx = ['total_x', 'total_y', 'total_z'].map((n) => series.find((s) => s.object === 'forces' && s.column === n)), vx = ['viscous_x', 'viscous_y', 'viscous_z'].map((n) => series.find((s) => s.object === 'forces' && s.column === n));
  if (fx.every(Boolean)) {
    const w = fx.map((s) => tail(s.t, s.y, frac)), mag = w[0].y.map((_, i) => Math.hypot(w[0].y[i], w[1].y[i], w[2].y[i])), fs = fluctuation(w[0].t, mag);
    M.forcePeak = fs.max; M.forceMean = fs.mean; M.forceStd = fs.std; series.push({ name: 'forces: |total|', object: 'forces', column: 'magnitude', t: w[0].t, y: mag });
    // fluctuating part of the force about its mean vector: what excites the structure
    const mv = w.map((s) => mean(s.y)); M.forceFluctuationPeak = Math.max(...w[0].y.map((_, i) => Math.hypot(w[0].y[i] - mv[0], w[1].y[i] - mv[1], w[2].y[i] - mv[2])));
  } else { M.forcePeak = null; M.forceMean = null; M.forceStd = null; M.forceFluctuationPeak = null; }
  if (vx.every(Boolean) && pos(ref.wallArea, 0) > 0) { const w = vx.map((s) => tail(s.t, s.y, frac)); M.wallShear = mean(w[0].y.map((_, i) => Math.hypot(w[0].y[i], w[1].y[i], w[2].y[i]))) / ref.wallArea; } else M.wallShear = null;
  const eMax = find('erosionMax'), eTot = find('erosionTotal'), pin = find('inletPressure'), pout = find('outletPressure'), dS = series.find((s) => /^(sauter|diameter)/i.test(s.object) || /\(d\.[\w]+\)/.test(s.column));
  if (eMax && eMax.t.length > 3) { const w = tail(eMax.t, eMax.y, 0.5), rate = Math.max(slope(w.t, w.y), 0), area = pos(ref.wallFaceArea, 0); M.erodedVolumeMax = eMax.y[eMax.y.length - 1]; M.erosionRateMmY = area > 0 ? (rate / area) * 1000 * SEC_PER_YEAR : null; } else { M.erodedVolumeMax = null; M.erosionRateMmY = null; }
  M.erodedVolume = eTot ? eTot.y[eTot.y.length - 1] : null;
  M.pressureDrop = pin ? stat(pin).mean - (pout ? stat(pout).mean : 0) : q1 && q3 ? q1.mean - q3.mean : null;
  M.sauterMean = dS ? stat(dS).mean : null;
  for (const k of Object.keys(M)) if (typeof M[k] === 'number' && !Number.isFinite(M[k])) M[k] = null;
  return { t: primary.t.slice(), series: series.map((s) => ({ name: s.name, object: s.object, t: s.t, y: s.y })), metrics: M, files: tables.length };
}
/**
 * Set imported values beside the application's own: rows [{ quantity, unit, external, inApp, ratio (external / in-app), difference (external − in-app) }].
 * pairs: [[quantity, unit, external, inApp]]; rows without an external value are dropped.
 */
export function compareResults(pairs) {
  return pairs.filter((p) => isNum(p[2])).map(([quantity, unit, external, inApp]) => ({ quantity, unit, external, inApp: isNum(inApp) ? inApp : null, ratio: isNum(inApp) && Math.abs(inApp) > 1e-300 ? external / inApp : null, difference: isNum(inApp) ? external - inApp : null }));
}
/** Comparison rows of an imported flow case with the one-dimensional values of the case (ref = the generator's `reference`). */
export function compareFlow(metrics, ref = {}) {
  const m = metrics || {};
  return compareResults([['Liquid holdup (time and volume average)', '–', m.holdupMean, ref.holdup], ['Liquid holdup at the middle section', '–', m.holdupSection, ref.holdup], ['Pressure gradient', 'Pa/m', m.pressureGradient, ref.dpdx], ['Wall shear stress', 'Pa', m.wallShear, ref.tauW],
    ['Slug frequency (threshold crossings)', '1/s', m.slugFrequency, ref.slugFreq], ['Slug frequency (spectrum peak)', '1/s', m.slugFrequencySpectrum, ref.slugFreq], ['Peak force on the wall', 'N', m.forcePeak, ref.forceN], ['Peak fluctuating force', 'N', m.forceFluctuationPeak, ref.forceN],
    ['Pressure swing at the middle section (half range)', 'Pa', m.pressureAmplitude, ref.pressureAmplitude], ['Cycle period', 's', m.pressurePeriod, ref.slugFreq > 0 ? 1 / ref.slugFreq : null], ['Erosion rate at the worst wall face', 'mm/y', m.erosionRateMmY, ref.erosionRateMmY], ['Pressure drop over the section', 'Pa', m.pressureDrop, isNum(ref.dpdx) && isNum(ref.sectionLength) ? ref.dpdx * ref.sectionLength : null], ['Sauter mean diameter', 'm', m.sauterMean, ref.sauter]]);
}

const vonMises = (s) => Math.sqrt(0.5 * ((s[0] - s[1]) ** 2 + (s[1] - s[2]) ** 2 + (s[2] - s[0]) ** 2) + 3 * (s[3] * s[3] + s[4] * s[4] + s[5] * s[5]));
/**
 * CalculiX .dat file: eigenvalue table, printed nodal displacements, total forces and printed element stresses.
 * Returns { frequencies: [Hz], displacements: [{ set, time, n, max, mean }], forces: [{ set, time, f: [x, y, z] }], stresses: [{ set, time, n, maxMises }],
 *   metrics: { f1, maxDisplacement, firstDisplacement, maxMises, peakReaction } }.
 */
export function importCalculixDat(src) {
  if (typeof src !== 'string') throw new Error('Expected text.');
  if (src.length > BRIDGE_LIMITS.chars) throw new Error('File is too large to import.');
  const lines = src.split(/\r?\n/), frequencies = [], displacements = [], forces = [], stresses = [];
  let mode = null, cur = null;
  const close = () => { if (cur && cur.n) { if (mode === 'disp') displacements.push({ set: cur.set, time: cur.time, n: cur.n, max: cur.max, mean: cur.sum / cur.n }); else if (mode === 'stress') stresses.push({ set: cur.set, time: cur.time, n: cur.n, maxMises: cur.max }); } cur = null; };
  for (const raw of lines) {
    const l = raw.trim();
    if (!l) continue;
    let m;
    if (/E I G E N V A L U E\s+O U T P U T/.test(l)) { close(); mode = 'eig'; continue; }
    if (/P A R T I C I P A T I O N|E F F E C T I V E|S T E P/.test(l)) { close(); mode = null; continue; }
    if ((m = /^displacements .* for set (\S+) and time\s+([-+\d.eE]+)/i.exec(l))) { close(); mode = 'disp'; cur = { set: m[1].slice(0, 40), time: +m[2], n: 0, max: 0, sum: 0 }; continue; }
    if ((m = /^stresses .* for set (\S+) and time\s+([-+\d.eE]+)/i.exec(l))) { close(); mode = 'stress'; cur = { set: m[1].slice(0, 40), time: +m[2], n: 0, max: 0, sum: 0 }; continue; }
    if ((m = /^total force .* for set (\S+) and time\s+([-+\d.eE]+)/i.exec(l))) { close(); mode = 'force'; cur = { set: m[1].slice(0, 40), time: +m[2] }; continue; }
    if (/^[A-Za-z]/.test(l)) { if (mode !== 'eig') { close(); mode = null; } continue; }
    const v = l.split(/\s+/).map(Number);
    if (!v.every(Number.isFinite)) continue;
    if (mode === 'eig' && v.length >= 4 && Number.isInteger(v[0]) && frequencies.length < 2000) frequencies.push(v[3]);
    else if (mode === 'disp' && v.length >= 4 && cur) { const u = Math.hypot(v[1], v[2], v[3]); cur.n++; cur.sum += u; if (u > cur.max) cur.max = u; if (displacements.length > BRIDGE_LIMITS.rows) throw new Error('Too many result blocks.'); }
    else if (mode === 'stress' && v.length >= 8 && cur) { const s = vonMises(v.slice(2, 8)); cur.n++; if (s > cur.max) cur.max = s; }
    else if (mode === 'force' && v.length >= 3 && cur) { forces.push({ set: cur.set, time: cur.time, f: v.slice(0, 3) }); cur = null; mode = null; }
  }
  close();
  const maxU = displacements.reduce((a, d) => Math.max(a, d.max), 0);
  return { frequencies, displacements, forces, stresses, metrics: { f1: frequencies.length ? frequencies[0] : null, modes: frequencies.length, maxDisplacement: displacements.length ? maxU : null, firstDisplacement: displacements.length ? displacements[0].max : null, maxMises: stresses.length ? stresses.reduce((a, d) => Math.max(a, d.maxMises), 0) : null, peakReaction: forces.length ? forces.reduce((a, f) => Math.max(a, Math.hypot(...f.f)), 0) : null } };
}
/**
 * CalculiX .frd result file (ASCII): nodal displacement and stress blocks of every step.
 * Returns { nodes, steps: [{ time, kind: 'DISP' | 'STRESS', n, max (largest displacement magnitude or von Mises stress), node }], frequencies: [Hz] (modal files),
 *   metrics: { maxDisplacement, maxMises, f1 } }.
 */
export function importCalculixFrd(src) {
  if (typeof src !== 'string') throw new Error('Expected text.');
  if (src.length > 5 * BRIDGE_LIMITS.chars) throw new Error('File is too large to import (limit ' + (5 * BRIDGE_LIMITS.chars) / 1e6 + ' million characters); lower the output frequency in the deck.');
  const steps = [], frequencies = [];
  let pos0 = 0, nodes = 0, time = 0, kind = null, cur = null, inNodes = false, modal = false, lineNo = 0;
  const n = src.length;
  while (pos0 < n) {
    let e = src.indexOf('\n', pos0); if (e < 0) e = n;
    const c0 = src.charCodeAt(pos0 + 1), c1 = src.charCodeAt(pos0 + 2); lineNo++;
    if (c0 === 45 && c1 === 49 && (cur || inNodes)) { // " -1": a node line
      if (inNodes) nodes++;
      else {
        const l = src.slice(pos0, e > pos0 && src.charCodeAt(e - 1) === 13 ? e - 1 : e), v = []; for (let k = 13; k + 12 <= l.length && v.length < 6; k += 12) v.push(+l.slice(k, k + 12));
        if (v.every(Number.isFinite)) { const mag = kind === 'DISP' && v.length >= 3 ? Math.hypot(v[0], v[1], v[2]) : kind === 'STRESS' && v.length >= 6 ? vonMises(v) : null; if (mag !== null) { cur.n++; if (mag > cur.max) { cur.max = mag; cur.node = +l.slice(3, 13); } } }
      }
    } else {
      const l = src.slice(pos0, Math.min(e, pos0 + 120));
      if (/^\s+2C/.test(l)) inNodes = true;
      else if (/^\s-3/.test(l)) { inNodes = false; if (cur) { if (cur.n) steps.push(cur); cur = null; kind = null; } }
      else if (/^\s+100CL/.test(l)) { time = +l.slice(12, 24); modal = /MODAL/.test(l); if (modal && Number.isFinite(time) && !frequencies.includes(time) && frequencies.length < 2000) frequencies.push(time); }
      else if (/^\s-4\s+(DISP|STRESS)\b/.test(l)) { kind = /DISP/.test(l) ? 'DISP' : 'STRESS'; cur = { time: Number.isFinite(time) ? time : 0, kind, n: 0, max: 0, node: 0, modal }; if (steps.length > 20000) throw new Error('Too many result blocks.'); }
      else if (/^\s-4\s/.test(l)) { cur = null; kind = null; }
    }
    pos0 = e + 1;
  }
  if (!nodes && !steps.length) throw new Error('Not a CalculiX .frd result file.');
  if (nodes > BRIDGE_LIMITS.nodes) throw new Error('The model has more nodes than this importer accepts.');
  const real = steps.filter((s) => !s.modal), mx = (k) => { const a = real.filter((s) => s.kind === k); return a.length ? a.reduce((m, s) => Math.max(m, s.max), 0) : null; };
  return { nodes, steps: steps.map(({ modal: _m, ...s }) => s), frequencies, metrics: { maxDisplacement: mx('DISP'), maxMises: mx('STRESS'), f1: frequencies.length ? frequencies[0] : null } };
}
/** Comparison rows of imported finite-element results with the hand values of the case (ref = `reference` of calculixPipeCase). */
export function compareStructure(metrics, ref = {}) {
  const m = metrics || {};
  return compareResults([['First natural frequency', 'Hz', m.f1, ref.f1], ['Maximum von Mises stress', 'MPa', isNum(m.maxMises) ? m.maxMises / 1e6 : null, isNum(ref.vonMises) ? ref.vonMises / 1e6 : null], ['von Mises stress / SMYS', '–', isNum(m.maxMises) && ref.smys > 0 ? m.maxMises / ref.smys : null, isNum(ref.vmUtil) ? ref.vmUtil : null],
    ['Maximum displacement', 'mm', isNum(m.maxDisplacement) ? m.maxDisplacement * 1000 : null, isNum(ref.staticDisp) ? ref.staticDisp * 1000 : null], ['Peak support reaction', 'kN', isNum(m.peakReaction) ? m.peakReaction / 1000 : null, isNum(ref.forceN) ? ref.forceN / 1000 : null]]);
}
/**
 * Small, storable record of an import for the case store (`store.setOutputs('bridge', …)`): numbers and short text only.
 * kind: generator name; metrics: plain numbers; comparison: rows from compareFlow / compareStructure / comparePropertyTables.
 */
export function bridgeRecord(kind, solver, metrics, comparison, note = '') {
  const m = {}; for (const [k, v] of Object.entries(metrics || {})) if (isNum(v) || typeof v === 'string') m[k] = typeof v === 'string' ? v.slice(0, 80) : +v.toPrecision(6);
  return { kind: String(kind).slice(0, 40), solver: String(solver || '').slice(0, 80), metrics: m, comparison: (comparison || []).slice(0, 20).map((r) => Object.fromEntries(Object.entries(r).map(([k, v]) => [k, isNum(v) ? +v.toPrecision(6) : v === null || v === undefined ? null : String(v).slice(0, 80)]))), note: String(note).slice(0, 300) };
}

// ---- preCICE: two-way fluid–structure interaction of a jumper, bend or span (OpenFOAM + CalculiX) -----------------
/**
 * Coupled case in the layout of the preCICE tutorials: precice-config.xml, fluid-openfoam/ (pimpleFoam on a moving mesh with the
 * OpenFOAM adapter) and solid-calculix/ (C3D8I wall with the CalculiX adapter). The fluid is the homogeneous mixture at the slug-body density;
 * the slug is a velocity pulse at the inlet. opts: { case, geometry: 'jumper' | 'bend' | 'span', coupling: 'implicit' | 'explicit',
 * cellsPerDiameter, nC, nR, elemLen, endTime, timeStep, plus the geometry options of calculixPipeCase }.
 */
export function fsiCase(opts = {}) {
  const c = withCase(opts), D = c.D, geometry = FE_GEOMETRY[opts.geometry] ? opts.geometry : 'jumper', implicit = opts.coupling !== 'explicit', legs = FE_GEOMETRY[geometry](opts, D), st = c.steel, name = safeName(opts.name || `fsi_${geometry}_${implicit ? 'implicit' : 'explicit'}`);
  const title = `two-way FSI of a pipe ${geometry}, ${implicit ? 'implicit' : 'explicit'} coupling`, nPerD = Math.round(clamp(num(opts.cellsPerDiameter, 12), 6, 80)), nc = Math.max(3, Math.round(0.4 * nPerD)), nr = Math.max(2, Math.round((nPerD - nc) / 2)), dz = (2 * D) / nPerD;
  const fluid = ogridMesh({ D, legs, nc, nr, wallRatio: 0.4, dz }), solid = pipeFeMesh({ D, wt: c.wt, legs, element: 'solid8', nC: opts.nC || 4 * nc, nR: opts.nR || 2, elemLen: opts.elemLen || dz });
  const HS = clamp(c.slug.holdupBody, 0.3, 1), rho = HS * c.rhoL + (1 - HS) * c.rhoG, nu = c.muL / c.rhoL, vm = c.vsl + c.vsg, vt = Math.max(pos(c.slug.velocity, 1.2 * vm), vm * 1.05);
  const ri = solid.ri, ro = solid.ro, As = Math.PI * (ro * ro - ri * ri), I = (Math.PI / 4) * (ro ** 4 - ri ** 4), rhoWall = st.rho + (c.z < 0 ? (1025 * Math.PI * ro * ro) / As : 0), mLin = rhoWall * As + rho * Math.PI * ri * ri;
  const f1 = beamFrequency(geometry === 'jumper' ? pos(opts.jumperLength, 12) : fluid.length, st.E * I, mLin, 'fixed'), dt = pos(opts.timeStep, Math.min((0.5 * dz) / (1.5 * vt), 1 / (40 * f1))), ramp = Math.max(D / vt, 10 * dt), endTime = pos(opts.endTime, Math.min(Math.max(4 / f1, 2 * fluid.length / vt), 20000 * dt)), nWin = Math.ceil(endTime / dt);
  const theta = geometry === 'jumper' ? Math.PI / 2 : (c.angleDeg * Math.PI) / 180, mid = fluid.stations[Math.floor(fluid.stations.length / 2)].c, files = [], F = (p, t) => files.push({ path: p, text: t });
  const pulse = [[0, vm], [2 * ramp, vm], [3 * ramp, vt], [3 * ramp + c.slug.length / vt, vt], [4 * ramp + c.slug.length / vt, vm], [Math.max(endTime, 5 * ramp + c.slug.length / vt), vm]];
  F('precice-config.xml', text(['<?xml version="1.0" encoding="UTF-8" ?>', `<!-- HydraSlug hand-off: ${title} (preCICE v3) -->`, '<precice-configuration>', '  <log>', '    <sink', '      filter="%Severity% > debug and %Rank% = 0"', '      format="---[precice] %ColorizedSeverity% %Message%"', '      enabled="true" />', '  </log>', '',
    '  <data:vector name="Force" />', '  <data:vector name="DisplacementDelta" />', '', '  <mesh name="Fluid-Mesh-Nodes" dimensions="3">', '    <use-data name="DisplacementDelta" />', '  </mesh>', '', '  <mesh name="Fluid-Mesh-Faces" dimensions="3">', '    <use-data name="Force" />', '  </mesh>', '', '  <mesh name="Solid-Mesh" dimensions="3">', '    <use-data name="DisplacementDelta" />', '    <use-data name="Force" />', '  </mesh>', '',
    '  <participant name="Fluid">', '    <receive-mesh name="Solid-Mesh" from="Solid" />', '    <provide-mesh name="Fluid-Mesh-Nodes" />', '    <provide-mesh name="Fluid-Mesh-Faces" />', '    <write-data name="Force" mesh="Fluid-Mesh-Faces" />', '    <read-data name="DisplacementDelta" mesh="Fluid-Mesh-Nodes" />',
    '    <mapping:nearest-neighbor', '      direction="write"', '      from="Fluid-Mesh-Faces"', '      to="Solid-Mesh"', '      constraint="conservative" />', '    <mapping:rbf direction="read" from="Solid-Mesh" to="Fluid-Mesh-Nodes" constraint="consistent">', `      <basis-function:compact-polynomial-c6 support-radius="${ff(1.5 * D, 4)}" />`, '    </mapping:rbf>', '  </participant>', '',
    '  <participant name="Solid">', '    <provide-mesh name="Solid-Mesh" />', '    <write-data name="DisplacementDelta" mesh="Solid-Mesh" />', '    <read-data name="Force" mesh="Solid-Mesh" />', `    <watch-point mesh="Solid-Mesh" name="Midpoint" coordinate="${mid.map((v) => ff(v, 6)).join(';')}" />`, '  </participant>', '', '  <m2n:sockets acceptor="Fluid" connector="Solid" exchange-directory=".." />', '',
    implicit ? ['  <coupling-scheme:serial-implicit>', '    <participants first="Fluid" second="Solid" />', `    <max-time-windows value="${nWin}" />`, `    <time-window-size value="${ff(dt, 4)}" />`, '    <exchange data="Force" mesh="Solid-Mesh" from="Fluid" to="Solid" />', '    <exchange data="DisplacementDelta" mesh="Solid-Mesh" from="Solid" to="Fluid" />', '    <max-iterations value="30" />', '    <relative-convergence-measure limit="1e-3" data="DisplacementDelta" mesh="Solid-Mesh" />', '    <relative-convergence-measure limit="1e-3" data="Force" mesh="Solid-Mesh" />', '    <acceleration:IQN-ILS>', '      <data name="DisplacementDelta" mesh="Solid-Mesh" />', '      <preconditioner type="residual-sum" />', '      <filter type="QR2" limit="1e-3" />', '      <initial-relaxation value="0.1" />', '      <max-used-iterations value="50" />', '      <time-windows-reused value="10" />', '    </acceleration:IQN-ILS>', '  </coupling-scheme:serial-implicit>']
      : ['  <coupling-scheme:serial-explicit>', '    <participants first="Fluid" second="Solid" />', `    <max-time-windows value="${nWin}" />`, `    <time-window-size value="${ff(dt, 4)}" />`, '    <exchange data="Force" mesh="Solid-Mesh" from="Fluid" to="Solid" />', '    <exchange data="DisplacementDelta" mesh="Solid-Mesh" from="Solid" to="Fluid" />', '  </coupling-scheme:serial-explicit>'], '</precice-configuration>']));
  // fluid participant
  const ts = turbulenceSetup('kOmegaSST', { U: vm, D, nu, title }), fp = 'fluid-openfoam/';
  for (const f of ts.files) F(fp + f.path.replace('0.orig/', '0/'), f.text);
  F(fp + '0/U', fieldFile('volVectorField', 'U', '[0 1 -1 0 0 0 0]', "uniform (0 0 0)", { inlet: ['type            uniformFixedValue;', `uniformValue    table (${pulse.map(([t, u]) => `(${ff(t, 6)} (${ff(u, 6)} 0 0))`).join(' ')});`, `value           uniform (${ff(vm, 6)} 0 0);`], outlet: ['type            inletOutlet;', 'inletValue      uniform (0 0 0);', 'value           uniform (0 0 0);'], wall: ['type            movingWallVelocity;', 'value           uniform (0 0 0);'] }, title));
  F(fp + '0/p', fieldFile('volScalarField', 'p', '[0 2 -2 0 0 0 0]', 'uniform 0', { inlet: ['type            zeroGradient;'], outlet: ['type            fixedValue;', 'value           uniform 0;'], wall: ['type            zeroGradient;'] }, title));
  F(fp + '0/pointDisplacement', fieldFile('pointVectorField', 'pointDisplacement', '[0 1 0 0 0 0 0]', 'uniform (0 0 0)', { inlet: ['type            fixedValue;', 'value           uniform (0 0 0);'], outlet: ['type            fixedValue;', 'value           uniform (0 0 0);'], wall: ['type            fixedValue;', 'value           $internalField;'] }, title));
  F(fp + 'constant/transportProperties', dictFile('transportProperties', 'constant', title, ['transportModel  Newtonian;', '', `nu              ${ff(nu, 6)};`]));
  F(fp + 'constant/turbulenceProperties', dictFile('turbulenceProperties', 'constant', title, ts.props));
  F(fp + 'constant/dynamicMeshDict', dictFile('dynamicMeshDict', 'constant', title, ['dynamicFvMesh   dynamicMotionSolverFvMesh;', '', 'motionSolverLibs (fvMotionSolvers);', '', 'solver          displacementLaplacian;', '', 'displacementLaplacianCoeffs', '{', '    diffusivity     quadratic inverseDistance (wall);', '}']));
  F(fp + 'system/blockMeshDict', blockMeshDict(fluid, title, 'inlet'));
  F(fp + 'system/controlDict', dictFile('controlDict', 'system', title, ['application     pimpleFoam;', '', 'startFrom       startTime;', 'startTime       0;', 'stopAt          endTime;', `endTime         ${ff(nWin * dt, 8)};`, `deltaT          ${ff(dt, 4)};`, '', 'writeControl    timeStep;', `writeInterval   ${Math.max(1, Math.round(nWin / 50))};`, 'purgeWrite      0;', 'writeFormat     binary;', 'writePrecision  10;', 'writeCompression off;', 'timeFormat      general;', 'timePrecision   8;', 'runTimeModifiable false;', 'adjustTimeStep  no;', '', 'functions', '{',
    '    forces', '    {', '        type            forces;', '        libs            (forces);', '        patches         (wall);', '        rho             rhoInf;', `        rhoInf          ${ff(rho, 6)};`, '        CofR            (0 0 0);', '        log             false;', '        writeControl    timeStep;', '        writeInterval   1;', '    }', '    preCICE_Adapter', '    {', '        type            preciceAdapterFunctionObject;', '        libs            ("libpreciceAdapterFunctionObject.so");', '    }', '}']));
  F(fp + 'system/fvSchemes', dictFile('fvSchemes', 'system', title, ['ddtSchemes', '{', '    default         Euler;', '}', '', 'gradSchemes', '{', '    default         Gauss linear;', '}', '', 'divSchemes', '{', '    default         none;', '    div(phi,U)      Gauss linearUpwind grad(U);', '    div(phi,k)      Gauss upwind;', '    div(phi,omega)  Gauss upwind;', '    div((nuEff*dev2(T(grad(U))))) Gauss linear;', '}', '', 'laplacianSchemes', '{', '    default         Gauss linear corrected;', '}', '', 'interpolationSchemes', '{', '    default         linear;', '}', '', 'snGradSchemes', '{', '    default         corrected;', '}', '', 'wallDist', '{', '    method          meshWave;', '}']));
  F(fp + 'system/fvSolution', dictFile('fvSolution', 'system', title, ['solvers', '{', '    "(p|pcorr)"', '    {', '        solver          GAMG;', '        smoother        GaussSeidel;', '        tolerance       1e-7;', '        relTol          0.01;', '    }', '', '    "(p|pcorr)Final"', '    {', '        solver          GAMG;', '        smoother        GaussSeidel;', '        tolerance       1e-8;', '        relTol          0;', '    }', '', '    "(U|k|omega|cellDisplacement)"', '    {', '        solver          smoothSolver;', '        smoother        symGaussSeidel;', '        tolerance       1e-8;', '        relTol          0.01;', '    }', '', '    "(U|k|omega|cellDisplacement)Final"', '    {', '        solver          smoothSolver;', '        smoother        symGaussSeidel;', '        tolerance       1e-8;', '        relTol          0;', '    }', '}', '', 'PIMPLE', '{', '    nOuterCorrectors    1;', '    nCorrectors         3;', '    nNonOrthogonalCorrectors 1;', '    correctPhi          yes;', '}']));
  F(fp + 'system/decomposeParDict', decomposeDict(Math.max(1, Math.round(num(opts.nProcs, 4))), title));
  F(fp + 'system/preciceDict', dictFile('preciceDict', 'system', title, ['preciceConfig "../precice-config.xml";', '', 'participant Fluid;', '', 'modules (FSI);', '', 'interfaces', '{', '    Interface1', '    {', '        mesh              Fluid-Mesh-Nodes;', '        patches           (wall);', '        locations         faceNodes;', '        readData', '        (', '            DisplacementDelta', '        );', '        writeData', '        (', '        );', '    };', '    Interface2', '    {', '        mesh              Fluid-Mesh-Faces;', '        patches           (wall);', '        locations         faceCenters;', '        readData', '        (', '        );', '        writeData', '        (', '            Force', '        );', '    };', '};', '', 'FSI', '{', `    rho rho [1 -3 0 0 0 0 0] ${ff(rho, 6)};`, '}']));
  F(fp + 'run.sh', text(['#!/bin/sh', 'set -e', 'cd "${0%/*}"', 'blockMesh > log.blockMesh 2>&1', 'checkMesh > log.checkMesh 2>&1', 'pimpleFoam > log.pimpleFoam 2>&1']));
  F(fp + 'case.foam', '');
  // solid participant
  F('solid-calculix/mesh.inc', calculixMeshText(solid, `${geometry} wall, ${solid.type}`));
  F('solid-calculix/pipe.inp', text([`** HydraSlug hand-off — CalculiX participant of the coupled case (${title})`, '** units m, kg, s, N, Pa', '*INCLUDE, INPUT=mesh.inc', '*MATERIAL, NAME=STEEL', '*ELASTIC', `${ff(st.E, 7)}, ${ff(st.nu, 4)}`, '*DENSITY', ff(rhoWall, 7), '*SOLID SECTION, ELSET=Eall, MATERIAL=STEEL', `*STEP, NLGEOM, INC=${nWin + 100}`, '*DYNAMIC, DIRECT', `${ff(dt, 4)}, ${ff(nWin * dt, 8)}`, '*RESTART, WRITE, FREQUENCY=1', '*BOUNDARY', 'END0, 1, 3', 'END1, 1, 3', '*CLOAD', 'Ninterface, 1, 0.0', 'Ninterface, 2, 0.0', 'Ninterface, 3, 0.0', `*NODE FILE, FREQUENCY=${Math.max(1, Math.round(nWin / 50))}`, 'U', `*EL FILE, FREQUENCY=${Math.max(1, Math.round(nWin / 50))}`, 'S', '*NODE PRINT, NSET=MID, FREQUENCY=1', 'U', '*END STEP']));
  F('solid-calculix/config.yml', text(['participants:', '    Solid:', '        interfaces:', '        - nodes-mesh: Solid-Mesh', '          patch: interface', '          read-data: [Force]', '          write-data: [DisplacementDelta]', '', 'precice-config-file: ../precice-config.xml']));
  F('solid-calculix/run.sh', text(['#!/bin/sh', 'set -e', 'cd "${0%/*}"', 'export OMP_NUM_THREADS=1', 'export CCX_NPROC_EQUATION_SOLVER=1', 'ccx_preCICE -i pipe -precice-participant Solid > log.ccx 2>&1']));
  F('Allrun', text(['#!/bin/sh', '# Starts both participants; preCICE connects them through sockets in this directory.', 'cd "${0%/*}" || exit', './fluid-openfoam/run.sh &', './solid-calculix/run.sh &', 'wait', 'echo "Import fluid-openfoam/postProcessing/forces/0/force.dat and solid-calculix/pipe.dat on the External solvers page."']));
  F('Allclean', text(['#!/bin/sh', 'cd "${0%/*}" || exit', 'rm -rf precice-run precice-*.json precice-*.log fluid-openfoam/[1-9]* fluid-openfoam/0.[0-9]* fluid-openfoam/constant/polyMesh fluid-openfoam/postProcessing fluid-openfoam/log.* fluid-openfoam/processor*', 'cd solid-calculix && rm -f pipe.cvg pipe.dat pipe.frd pipe.sta pipe.12d pipe.rout spooles.out log.ccx']));
  const cost = (fluid.nCells * nWin * (implicit ? 6 : 1.5) * COST_PER_CELL_STEP) / 3600;
  const readme = text([`# ${name} — OpenFOAM + CalculiX through preCICE`, '', `Coupled simulation of the flow in a pipe ${geometry} and the motion of its wall: ${implicit ? 'implicit (strong) coupling with quasi-Newton acceleration, the fluid and the structure iterate to equilibrium in every time window' : 'explicit (loose) coupling, one exchange of force and displacement per time window — the wall load from the CFD drives the finite-element model and its motion is fed back with one window of lag'}.`, '',
    'The integrity suite applies a slug force formula to a beam model. Here the pressure and shear field on the wall come from the Navier–Stokes solution, the wall is a three-dimensional elastic solid, and', 'the wall motion changes the flow passage (added mass, damping and any lock-in follow from the solution).', '',
    '## Model', '', ...kv([['Fluid', `pimpleFoam (k–ω SST) on a moving mesh, ${fluid.nCells.toLocaleString('en-US')} cells; homogeneous mixture at the slug-body density ${ff(rho, 5)} kg/m³, kinematic viscosity ${ff(nu, 4)} m²/s`], ['Slug', `velocity pulse at the inlet from ${ff(vm, 4)} to ${ff(vt, 4)} m/s lasting ${ff(c.slug.length / vt, 3)} s (slug length ${ff(c.slug.length, 4)} m)`], ['Structure', `CalculiX ${solid.type}, ${solid.nodes.length.toLocaleString('en-US')} nodes, ${solid.elements.length.toLocaleString('en-US')} elements, both ends fixed; ${st.grade}, E = ${ff(st.E / 1e9, 4)} GPa; equivalent wall density ${ff(rhoWall, 5)} kg/m³${c.z < 0 ? ' (added mass of the sea water outside)' : ''}`], ['Interface', 'fluid patch `wall` ↔ solid node set `Ninterface` (inner wall); forces mapped conservatively (nearest neighbour), displacements by radial basis functions'], ['Time window', `${ff(dt, 3)} s × ${nWin.toLocaleString('en-US')} windows = ${ff(nWin * dt, 4)} s (beam estimate of the first natural frequency ${ff(f1, 3)} Hz)`], ['Estimated cost', hours(cost)]]),
    '- The fluid is single-phase: the density wave of a real slug is represented by the higher density and the velocity pulse, not by a resolved gas–liquid interface. For a resolved slug replace the fluid participant by', '  the interFoam case of the “pipe section” hand-off (same mesh generator) and keep this coupling set-up; the adapter then needs the density field (`FSI { }` without a constant rho).', '',
    '## Requirements', '', '- preCICE v3 (precice.org), the OpenFOAM adapter (`libpreciceAdapterFunctionObject.so`, v1.3 or later) and the CalculiX adapter (`ccx_preCICE`, built for CalculiX 2.20).', '',
    '## Run', '', '```sh', `unzip ${name}.zip && cd ${name}`, 'precice-tools check precice-config.xml      # optional: validates the coupling configuration', './Allrun                                    # or start fluid-openfoam/run.sh and solid-calculix/run.sh in two terminals', '```', '',
    '## Results', '', '- `fluid-openfoam/postProcessing/forces/0/force.dat` — force of the fluid on the wall in time;', '- `solid-calculix/pipe.dat` — displacement history of the mid ring (node set MID); `pipe.frd` — displacement and stress fields;', '- `precice-Solid-watchpoint-Midpoint.log` — force and displacement at the interface point nearest the middle of the path.', '', 'Import `force.dat`, `pipe.dat` and `pipe.frd` on the External solvers page.', '']);
  F('README.md', readme);
  return { name, files, readme, fluid, solid, summary: `preCICE ${implicit ? 'implicit' : 'explicit'} FSI: pimpleFoam ${fluid.nCells.toLocaleString('en-US')} cells + CalculiX ${solid.elements.length.toLocaleString('en-US')} ${solid.type} elements, ${nWin.toLocaleString('en-US')} time windows, about ${hours(cost)}.`,
    commands: [`unzip ${name}.zip && cd ${name}`, 'precice-tools check precice-config.xml', './Allrun'], plan: { cells: fluid.nCells, steps: nWin, coreHours: cost, endTime: nWin * dt, dt }, reference: { f1, forceN: c.slug.forceN, smys: st.smys, wallArea: Math.PI * D * fluid.length } };
}
