// Coverage of the reference catalogue: for every suite and every list, how many items are solved in-app, how many
// are handed to an external solver, and which are covered by neither.   node tools/coverage.mjs [--strict]
import { CATALOG } from '../js/data/catalog.js';
import { SUITES, loadSuite } from '../js/suites/index.js';
const norm = (s) => String(s).toLowerCase().normalize('NFKD').replace(/[‐-―]/g, '-').replace(/[’'`]/g, '').replace(/s\b/g, '').replace(/[^a-z0-9]+/g, ' ').trim();
let bridge = null; try { bridge = await import('../js/core/bridge.js'); } catch { /* not built yet */ }
const LISTS = ['classical', 'hybrid', 'initial', 'boundary', 'inputs', 'outputs', 'calibration', 'verification', 'validation'];
let total = 0, inApp = 0, viaExt = 0; const open = [];
for (const m of SUITES) {
  const s = await loadSuite(m.id), k = (s.implemented || []).map(norm).filter(Boolean), r = (s.referenceOnly || []).map(norm).filter(Boolean), parts = [];
  for (const key of LISTS) {
    const items = CATALOG[m.num][key] || []; let a = 0, b = 0;
    for (const it of items) { const q = norm(it), on = !r.some((z) => q.includes(z)) && k.some((z) => q.includes(z) || z.includes(q)); if (on) a++; else if (bridge?.handoffFor?.(m.id, it)) b++; else open.push(`${m.id} · ${key} · ${it}`); }
    total += items.length; inApp += a; viaExt += b; parts.push(`${key} ${a}${b ? '+' + b : ''}/${items.length}`);
  }
  console.log(m.id.padEnd(7), parts.join(' · '), `| reference data sets ${s.validationData?.length || 0} | provenance entries ${s.provenance?.length || 0}`);
}
console.log(`\n${inApp} solved in-app + ${viaExt} via external solver = ${inApp + viaExt} of ${total} catalogue items; ${open.length} not covered`);
if (open.length) console.log(open.map((x) => '  - ' + x).join('\n'));
if (process.argv.includes('--strict') && open.length) process.exit(1);
