// Regression test runner: node tests/run.mjs [suiteId ...]
// For every suite: contract shape, default run, every preset, linked-data hooks, calibration model,
// verification checks, and a recursive scan for NaN / Infinity / undefined in everything that is shown.
import { SUITES, loadSuite } from '../js/suites/index.js';
import { DEFAULT_FLUID, makeFluid, flashPT, phaseProps, fluidModel, hydrateT0 } from '../js/core/thermo.js';
import { frictionFactor, marchSteady, slugVelocity } from '../js/core/pipe.js';
import { BASE } from '../js/data/basecase.js';
import * as P from '../js/core/props.js';
import { gci, brent, rk45, levenbergMarquardt, metrics } from '../js/core/num.js';
import { CATALOG } from '../js/data/catalog.js';

let fails = 0;
const ok = (cond, msg) => { if (!cond) { fails++; console.log('   ✗ ' + msg); } };
const near = (a, b, tol, msg) => ok(Math.abs(a - b) <= tol, `${msg}: got ${a}, expected ${b} ± ${tol}`);

function scan(x, path, out, depth = 0) {
  if (depth > 9) return;
  if (typeof x === 'number') { if (!Number.isFinite(x)) out.push(path + ' = ' + x); return; }
  if (x === undefined) { out.push(path + ' = undefined'); return; }
  if (x === null || typeof x !== 'object') return;
  if (ArrayBuffer.isView(x)) { for (let i = 0; i < x.length; i++) if (!Number.isFinite(x[i])) { out.push(path + `[${i}] = ${x[i]}`); break; } return; }
  if (Array.isArray(x)) { for (let i = 0; i < x.length; i++) { scan(x[i], `${path}[${i}]`, out, depth + 1); if (out.length > 6) return; } return; }
  for (const k of Object.keys(x)) { if (k.startsWith('_')) continue; scan(x[k], path + '.' + k, out, depth + 1); if (out.length > 6) return; }
}
function checkResult(res, label) {
  ok(res && Array.isArray(res.kpis) && res.kpis.length >= 6, `${label}: at least 6 KPIs`);
  ok(typeof res.summary === 'string' && res.summary.length > 20, `${label}: has a summary`);
  ok(Array.isArray(res.recommendations) && res.recommendations.length >= 1, `${label}: has recommendations`);
  try { structuredClone(res); } catch (e) { ok(false, `${label}: result cannot cross threads (${e.message})`); }
  ok(Array.isArray(res.plots) && res.plots.length >= 1, `${label}: has plots`);
  ok(res.outputs && typeof res.outputs === 'object', `${label}: has outputs`);
  const bad = [];
  scan(res.kpis, 'kpis', bad); scan(res.tables?.map((t) => t.rows), 'tables', bad); scan(res.outputs, 'outputs', bad); scan(res.balances || [], 'balances', bad);
  for (const [i, p] of (res.plots || []).entries()) {
    if (p.type === 'field') { for (let j = 0; j < p.z.length; j++) for (let k = 0; k < p.z[j].length; k++) if (!Number.isFinite(p.z[j][k]) && !(p.mask && p.mask[j][k])) { bad.push(`plots[${i}].z[${j}][${k}]`); j = 1e9; break; } ok(p.z.length === p.y.length && p.z[0].length === p.x.length, `${label}: field plot ${i} dimensions`); }
    else if (p.type === 'bar') scan(p.series, `plots[${i}].series`, bad);
    else for (const s of p.series || []) ok(s.x.length === s.y.length, `${label}: plot ${i} series "${s.name}" x/y length`);
  }
  ok(!bad.length, `${label}: non-finite values → ${bad.slice(0, 5).join('; ')}`);
  for (const b of res.balances || []) ok(Math.abs(b.in - b.out) <= 1e-3 * Math.max(1, Math.abs(b.in)), `${label}: balance "${b.name}" closes (in ${b.in}, out ${b.out})`);
}

console.log('Core');
near(P.density(25, 35), 1023.3, 0.6, 'seawater density 25 °C, 35 g/kg');
near(P.viscosity(25, 0) * 1000, 0.89, 0.01, 'water viscosity at 25 °C');
near(brent((x) => x * x - 2, 0, 2), Math.SQRT2, 1e-9, 'brent root');
near(rk45((t, y) => [-y[0]], [1], 0, 1).y.at(-1)[0], Math.exp(-1), 1e-5, 'rk45 exponential decay');
near(gci([1, 2, 4], [1 + 0.01, 1 + 0.04, 1 + 0.16]).p, 2, 1e-6, 'GCI observed order for a second-order sequence');
near(levenbergMarquardt((p) => [1, 2, 3, 4].map((x) => p[0] * x + p[1] - (3 * x + 1)), [1, 0]).p[0], 3, 1e-6, 'Levenberg–Marquardt line fit');
near(metrics([1, 2, 3], [1, 2, 3]).rmse, 0, 1e-12, 'metrics of perfect agreement');
{
  const c1 = makeFluid({ comp: { C1: 100 } }), c3 = makeFluid({ comp: { C3: 100 } });
  near(phaseProps(c1, [1], 100, 26.85, 'vapour').Z, 0.851, 0.01, 'methane compressibility at 100 bara, 300 K (Peng–Robinson)');
  near(phaseProps(c3, [1], 15, 26.85, 'liquid').rho, 489, 12, 'liquid propane density at 300 K');
  const f = makeFluid(DEFAULT_FLUID), fl = flashPT(f, 100, 60);
  ok(fl.phase === 'two' && fl.beta > 0.2 && fl.beta < 0.45, 'reference oil is two-phase at 100 bara, 60 °C');
  near(fl.x.reduce((s, v) => s + v, 0), 1, 1e-10, 'liquid mole fractions sum to one');
  near(fl.beta * fl.y[2] + (1 - fl.beta) * fl.x[2], f.z[2], 1e-9, 'methane material balance over the flash');
  near(frictionFactor(1e5, 1e-4), 0.01851, 2e-4, 'Colebrook friction factor at Re = 1e5, ε/D = 1e-4');
  near(frictionFactor(1000, 0), 0.064, 1e-12, 'laminar friction factor 64/Re');
  near(hydrateT0(70, 0.6), 15.5, 2, 'hydrate temperature of a 0.6-gravity gas at 70 bara (gas-gravity chart: about 15.5 °C)');
  near(slugVelocity(2, 0.254, 0).C0, 1.05, 1e-12, 'Bendiksen distribution coefficient, horizontal, low Froude number');
  const fm = fluidModel({}), r = marchSteady({ fm, profile: { x: BASE.profile.map((p) => p.x), z: BASE.profile.map((p) => p.z) }, id: BASE.idMm / 1000, U: BASE.U, tAmb: 4, tIn: BASE.tIn, pOut: BASE.pOut, n: 100 });
  ok(r.ok && r.pIn > 60 && r.pIn < 140, `reference tie-back inlet pressure is plausible (${r.ok ? r.pIn.toFixed(1) : r.reason} bara)`);
  near(r.pOut, BASE.pOut, 1e-3, 'shooting meets the outlet pressure');
  ok(r.tOut < BASE.tIn && r.tOut > 4, 'arrival temperature lies between ambient and inlet');
}
ok(Object.keys(CATALOG).length === 7, 'catalogue has 7 suites');

const only = process.argv.slice(2);
const ctx = () => ({ fluid: JSON.parse(JSON.stringify(DEFAULT_FLUID)), site: { data: {} }, outputs, progress() {}, tick: async () => {} });
const outputs = {};
for (const meta of SUITES) {
  if (only.length && !only.includes(meta.id)) continue;
  const t0 = Date.now();
  let s;
  try { s = await loadSuite(meta.id); } catch (e) { fails++; console.log(`${meta.num}. ${meta.id}: ✗ cannot load — ${e.message}`); continue; }
  console.log(`${meta.num}. ${s.title}`);
  ok(s.id === meta.id && s.num === meta.num, 'id/num match the registry');
  for (const k of ['title', 'tagline', 'inputs', 'run', 'verify']) ok(s[k], `declares ${k}`);
  const fields = s.inputs.flatMap((g) => g.fields), keys = new Set();
  for (const f of fields) { ok(f.key && f.label, `field has key and label (${f.key})`); ok(!keys.has(f.key), `duplicate field key ${f.key}`); keys.add(f.key); if (!f.type || f.type === 'number') ok(typeof f.value === 'number', `numeric default for ${f.key}`); }
  const d = Object.fromEntries(fields.map((f) => [f.key, f.type === 'table' ? JSON.parse(JSON.stringify(f.value)) : f.value ?? null]));
  try {
    const tr = Date.now(), res = await s.run({ ...d }, ctx());
    checkResult(res, 'default run'); ok(Date.now() - tr < 6000, `default run finishes within 6 s (${Date.now() - tr} ms)`);
    outputs[s.id] = res.outputs;
    for (const p of s.presets || []) { try { checkResult(await s.run({ ...d, ...JSON.parse(JSON.stringify(p.values)) }, ctx()), `preset "${p.name}"`); } catch (e) { ok(false, `preset "${p.name}" threw: ${e.message}`); } }
    // linked data must only offer known keys and usable values
    for (const hook of ['pull', 'site']) { if (!s[hook]) continue; const items = (hook === 'pull' ? s.pull(ctx()) : s.site({ data: { depth: 1350, elevation: 0, seabedTemp: 4, sst: 27, salinity: 35, currentSpeed: 0.3, currentDir: 90, waveHeight: 2, wavePeriod: 9, tideRange: 1.2, airTemp: 28, airTempMin: 22, airTempMax: 33, windSpeed: 6, groundTemp: 18, pressure: 1012, seabedSlope: 1.5, inflation: 3, lendingRate: 8, fxPerUSD: 1500, currency: 'NGN', electricityPrice: 0.1, gridCarbon: 0.4, oilPrice: 80, gasPrice: 3.5, carbonPrice: 40, taxRate: 30, country: 'Nigeria' } }) || []).filter(Boolean); for (const it of items) ok(keys.has(it.key), `${hook} offers unknown key "${it.key}"`); }
    // a second run with everything pulled from the other suites must still work
    const linked = { ...d }; for (const it of (s.pull?.(ctx()) || []).filter(Boolean)) if (keys.has(it.key) && it.value != null) linked[it.key] = it.value;
    try { checkResult(await s.run(linked, ctx()), 'linked run'); } catch (e) { ok(false, `linked run threw: ${e.message}`); }
    if (s.mesh) for (const st of Array.isArray(s.mesh) ? s.mesh : [s.mesh]) { for (const k of st.keys) ok(keys.has(k), `mesh key ${k} exists`); const r1 = await s.run({ ...d }, ctx()); for (const m of st.metrics) { let val; try { val = +m.get(r1); } catch (e) { if (/select the task/.test(e.message)) continue; throw e; } ok(Number.isFinite(val), `mesh metric "${m.label}" is finite`); } }
  } catch (e) { ok(false, `default run threw: ${e.stack?.split('\n').slice(0, 3).join(' | ')}`); }
  if (s.calibration) {
    const c = s.calibration;
    try { const m = c.model({ ...d }); for (const t of c.targets) ok(Number.isFinite(m[t.key]), `calibration model returns ${t.key}`); ok(c.sample.length >= 8, 'calibration sample has rows'); for (const p of c.params) ok(keys.has(p.key), `calibration parameter ${p.key} is an input`); }
    catch (e) { ok(false, `calibration threw: ${e.message}`); }
  }
  if (s.validationData) {
    ok(Array.isArray(s.validationData), 'validationData is an array');
    for (const ds of s.validationData || []) {
      const lab = `reference data "${ds.id}"`;
      ok(ds.id && ds.title && ds.target && Array.isArray(ds.rows) && ds.rows.length >= 8, `${lab}: id, title, target and ≥ 8 rows`);
      ok(ds.source && /^https?:\/\//.test(ds.source.url || '') && (ds.source.citation || '').length > 20 && ds.source.licence && /^\d{4}-\d\d-\d\d$/.test(ds.source.retrieved || ''), `${lab}: citation, address, licence and retrieval date`);
      try {
        const meas = [], pred = [];
        for (const r of ds.rows) { const m = +r[ds.target], q = +ds.model(r); if (Number.isFinite(m)) { meas.push(m); pred.push(q); } }
        ok(pred.every(Number.isFinite), `${lab}: model returns finite predictions`);
        const mt = metrics(meas, pred), t = ds.tolerance || {};
        const within = (t.mape === undefined || mt.mape <= t.mape) && (t.rmse === undefined || mt.rmse <= t.rmse) && (t.bias === undefined || Math.abs(mt.bias) <= t.bias) && (t.maxAbs === undefined || Math.max(...mt.residuals.map(Math.abs)) <= t.maxAbs);
        ok(within, `${lab}: outside its stated tolerance (MAPE ${mt.mape?.toFixed(2)} %, RMSE ${mt.rmse?.toPrecision(3)}, bias ${mt.bias?.toPrecision(3)})`);
        console.log(`   ${ds.id}: n = ${meas.length}, MAPE ${mt.mape?.toFixed(2)} %, bias ${mt.bias?.toPrecision(3)} ${ds.unit || ''}`);
      } catch (e) { ok(false, `${lab} threw: ${e.message}`); }
    }
  }
  try { const checks = await s.verify(); const failed = checks.filter((c) => !c.pass); ok(checks.length >= 12, `has ≥ 12 verification checks (${checks.length})`); ok(!failed.length, `verification failures: ${failed.map((c) => `${c.name} (got ${c.got}, expected ${c.expected})`).join('; ')}`); console.log(`   ${checks.length - failed.length}/${checks.length} verification checks pass · ${Date.now() - t0} ms`); }
  catch (e) { ok(false, `verify threw: ${e.message}`); }
}
console.log(fails ? `\n${fails} FAILURE(S)` : '\nAll tests passed');
process.exit(fails ? 1 : 0);
