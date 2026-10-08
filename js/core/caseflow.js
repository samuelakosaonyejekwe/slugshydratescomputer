// Case-level steady flow used as the common fallback: when the flow suite has not been run yet, the suites
// downstream of it (solids, operations, integrity, economics) still need pressure, temperature, holdup and
// velocity along the line. This module builds that picture from whatever the case already holds — the network
// suite's geometry when published, otherwise the reference case — and caches it.
import { fluidModel } from './thermo.js';
import { marchSteady, seaTemperature } from './pipe.js';
import { interp1 } from './num.js';
import { BASE } from '../data/basecase.js';

const isArr = (a, n = 2) => Array.isArray(a) && a.length >= n && a.every((v) => typeof v === 'number' && Number.isFinite(v));
/**
 * Line description of the case: geometry published by the network suite, else the reference case.
 * Returns { profile: { x[], z[] }, id (m), wt (m), roughness (m), uValue (W/m²K), tSeabed, tSeaSurface (°C), riserBaseX (m), length (m), tIn (°C), pOut (bara), fromNet (bool) }.
 */
export function caseLine(ctx = {}, over = {}) {
  const n = ctx.outputs?.net || {}, ok = n.profile && isArr(n.profile.x) && isArr(n.profile.z) && n.profile.x.length === n.profile.z.length;
  const profile = over.profile || (ok ? { x: n.profile.x.slice(), z: n.profile.z.slice() } : { x: BASE.profile.map((p) => p.x), z: BASE.profile.map((p) => p.z) });
  const num = (a, b) => (typeof a === 'number' && Number.isFinite(a) && a > 0 ? a : b);
  let length = 0; for (let i = 1; i < profile.x.length; i++) length += Math.hypot(profile.x[i] - profile.x[i - 1], profile.z[i] - profile.z[i - 1]);
  return {
    profile, length, fromNet: !!ok,
    id: num(over.id, num(n.id, BASE.idMm / 1000)), wt: num(over.wt, num(n.wt, BASE.wtMm / 1000)), roughness: num(over.roughness, num(n.roughness, BASE.roughUm * 1e-6)), uValue: num(over.uValue, num(n.uValue, BASE.U)),
    tSeabed: over.tSeabed ?? n.tSeabed ?? ctx.site?.data?.seabedTemp ?? BASE.tSeabed, tSeaSurface: over.tSeaSurface ?? n.tSeaSurface ?? ctx.site?.data?.sst ?? BASE.tSeaSurface,
    riserBaseX: num(over.riserBaseX, num(n.riserBaseX, ok ? profile.x[profile.x.length - 1] : BASE.riserBaseX)),
    tIn: over.tIn ?? ctx.fluid?.Tin ?? BASE.tIn, pOut: over.pOut ?? ctx.fluid?.Pout ?? n.separatorP ?? BASE.pOut,
  };
}
/** Ambient temperature (°C) at an elevation z (m, negative below sea level): thermocline in the sea, air above it. */
export const ambientAt = (z, line, tAir = BASE.tAir) => (z >= 0 ? tAir : seaTemperature(-z, line.tSeaSurface, line.tSeabed));

const cache = new Map();
/**
 * Steady pressure/temperature/holdup solution of the case line (kernel marchSteady result plus `line` and `fm`).
 * opt: { mScale, model, n, idOf(s), roughOf(s), pIn, pOut, tIn, uValue, override (fluid override), ...line overrides }.
 * Throws when no steady solution exists (the message says why).
 */
export function steadyCase(ctx = {}, opt = {}) {
  const fm = fluidModel(ctx, opt.override), line = caseLine(ctx, opt), cacheable = typeof opt.idOf !== 'function' && typeof opt.roughOf !== 'function';
  const key = cacheable ? JSON.stringify([fm.spec, line, opt.mScale, opt.model, opt.n, opt.pIn, ctx.outputs?.pvt?.table ? 1 : 0]) : null;
  if (key && cache.has(key)) return cache.get(key);
  const r = marchSteady({ fm, profile: line.profile, id: line.id, rough: line.roughness, U: line.uValue, tAmbOf: (s, z) => ambientAt(z, line), tIn: line.tIn, ...(opt.pIn !== undefined ? { pIn: opt.pIn } : { pOut: line.pOut }), mScale: opt.mScale ?? 1, model: opt.model || 'beggsBrill', n: opt.n || 120, idOf: opt.idOf, roughOf: opt.roughOf });
  if (!r.ok) throw new Error(r.reason || 'No steady flow solution for this case.');
  r.line = line; r.fm = fm;
  if (key) { if (cache.size > 12) cache.clear(); cache.set(key, r); }
  return r;
}
/**
 * Flow picture for downstream suites: the flow suite's published profile when it exists, otherwise the kernel's
 * steady solution, both in the same shape:
 * { x[], z[], P[], T[], holdup[], vsl[], vsg[], vm[], rhoM[], dpdx[], tauW[], tAmb[], tHyd[], subcooling[], regime[], source: 'flow suite' | 'kernel estimate' }.
 */
export function flowPicture(ctx = {}, opt = {}) {
  const p = ctx.outputs?.flow?.profile, keys = ['x', 'z', 'P', 'T', 'holdup', 'vsl', 'vsg', 'vm', 'rhoM', 'dpdx', 'tauW', 'tAmb', 'tHyd', 'subcooling'];
  if (!opt.force && p && keys.every((k) => isArr(p[k]) && p[k].length === p.x.length) && Array.isArray(p.regime)) return { ...Object.fromEntries(keys.map((k) => [k, p[k].slice()])), regime: p.regime.slice(), source: 'flow suite' };
  const r = steadyCase(ctx, opt), out = { source: 'kernel estimate' };
  for (const k of keys) out[k] = r[k].slice();
  out.regime = r.regime.slice();
  return out;
}
/** Value of a profile array at distance x (linear interpolation). */
export const atX = (pic, key, x) => interp1(pic.x, pic[key], x);
