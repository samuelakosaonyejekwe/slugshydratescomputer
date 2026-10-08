// External solvers page: hands the formulations that cannot run in a browser at useful resolution to established open
// solvers. One card per hand-off — what it resolves, which solver, the estimated size and cost, a ready-to-run case
// as one .zip, the commands to run it, and an importer that sets the returned results beside the in-app values.
import { h, fill, btn, kpiGrid, dataTable, toast, badge, fieldRow } from '../core/ui.js';
import { store } from '../core/store.js';
import { plotCard } from '../core/plot.js';
import { download, checkFile } from '../core/io.js';
import { fmt, interp1 } from '../core/num.js';
import { fluidModel } from '../core/thermo.js';
import { flowPicture } from '../core/caseflow.js';
import * as B from '../core/bridge.js';

const card = (...kids) => h('section', { class: 'card' }, ...kids);
const pc = (spec) => plotCard(spec, { onDownload: download });
const opt = (pairs) => pairs.map(([value, label]) => ({ value, label }));
const N = (key, label, value, unit, min, max, help) => ({ key, label, value, unit, min, max, help });
const S = (key, label, value, options, help) => ({ key, label, type: 'select', value, options: opt(options), help });
const TURB = Object.entries(B.TURBULENCE).filter(([k]) => k !== 'laminar').map(([k, v]) => [k, v.label]);
const GEOM = [['span', 'Free span (straight)'], ['bend', 'Bend with two legs'], ['jumper', 'Rigid jumper (up – across – down)']];

// The hand-offs offered on this page. kind selects the importer: 'flow' (OpenFOAM time series), 'structure' (CalculiX), 'property' (CSV table).
const CARDS = [
  { id: 'pipe', icon: '🌀', title: 'Pipe section in three dimensions — turbulence and interface', generator: 'openfoamPipeCase', kind: 'flow', solver: 'OpenFOAM interFoam / interIsoFoam',
    resolves: 'The velocity field, the gas–liquid interface and the turbulence in a section of the line at the chosen location: RANS (k–ε family, k–ω SST, Spalart–Allmaras, Reynolds stress), LES, DES, IDDES or a DNS-resolution set-up, with an algebraic or a geometric volume-of-fluid interface. The in-app models give only section averages from closure laws.',
    fields: [S('turbulence', 'Turbulence treatment', 'kOmegaSST', TURB, 'RANS models average all turbulence; LES and DES resolve the large eddies and need a far finer mesh and time step.'), S('interface', 'Interface method', 'vof', [['vof', 'Algebraic VOF (interFoam, MULES)'], ['isoAdvector', 'Geometric VOF, reconstructed interface (interIsoFoam)'], ['plicRDF', 'Geometric VOF with reconstructed distance function (interIsoFoam, plicRDF)']]), S('geometry', 'Geometry', 'straight', [['straight', 'Straight section at the local inclination'], ['bend', '90° bend with two straight legs']]),
      N('lengthD', 'Section length', 20, 'diameters', 3, 2000, 'Developed length including the bend. Slug statistics need several slug units; short sections show initiation and the front structure.'), N('cellsPerDiameter', 'Cells across the diameter', 24, '', 8, 400), N('yPlus', 'Target y+ of the first cell', 50, '', 0.2, 300, 'About 30–100 with wall functions; about 1 for a wall-resolved simulation.'), N('flowThroughs', 'Simulated time', 5, 'flow-throughs', 1, 200)] },
  { id: 'coupled', icon: '🔁', title: 'One-dimensional transient driving a three-dimensional section', generator: 'openfoamCoupledCase', kind: 'flow', solver: 'OpenFOAM interFoam with flow-rate tables from the 1-D model',
    resolves: 'The three-dimensional section receives the liquid and gas flow-rate histories of the in-app transient model at its inlet (a slug train from the slug closures when no transient has been run) and returns holdup, pressure gradient, wall shear and forces to check the 1-D closures.',
    fields: [S('turbulence', 'Turbulence treatment', 'kOmegaSST', TURB), S('geometry', 'Geometry', 'straight', [['straight', 'Straight section at the local inclination'], ['bend', '90° bend with two straight legs']]), N('lengthD', 'Section length', 20, 'diameters', 3, 2000), N('cellsPerDiameter', 'Cells across the diameter', 24, '', 8, 400)] },
  { id: 'riser', icon: '⤴️', title: 'Pipeline–riser severe slugging', generator: 'openfoamRiserSlugCase', kind: 'flow', solver: 'OpenFOAM compressibleInterFoam (or interFoam)',
    resolves: 'The blocking, gas-compression, blow-out and fall-back cycle at the riser base, resolved in a flowline–bend–riser model with a compressible gas. The in-app model screens it with stability criteria and a lumped cycle.',
    fields: [S('dimension', 'Model', '2d', [['2d', 'Planar two-dimensional channel'], ['3d', 'Coarse three-dimensional pipe']]), S('solver', 'Gas phase', 'compressibleInterFoam', [['compressibleInterFoam', 'Compressible (needed for the slugging cycle)'], ['interFoam', 'Incompressible']]), N('flowlineLengthD', 'Flowline length modelled', 120, 'diameters', 10, 5000), N('riserHeightD', 'Riser height modelled', 60, 'diameters', 5, 5000), N('cellsPerDiameter', 'Cells across the diameter', 16, '', 6, 200), N('cycles', 'Simulated time', 4, 'riser fill times', 1, 100)] },
  { id: 'euler', icon: '🫧', title: 'Two-fluid and multi-fluid model in three dimensions', generator: 'openfoamEulerCase', kind: 'flow', solver: 'OpenFOAM reactingTwoPhaseEulerFoam / multiphaseEulerFoam',
    resolves: 'Euler–Euler transport with one momentum equation per phase and resolved interfacial forces (drag, lift, virtual mass, wall lubrication, turbulent dispersion, bubble-induced turbulence) — or solids as a continuum with kinetic-theory closures.',
    fields: [S('phases', 'Phases', '2', [['2', 'Two fluids'], ['3', 'Three fluids: gas, oil, water']]), S('dispersed', 'Dispersed phase (two fluids)', 'gas', [['gas', 'Gas in liquid'], ['solids', 'Solid particles in liquid']]), N('lengthD', 'Section length', 20, 'diameters', 3, 500), N('cellsPerDiameter', 'Cells across the diameter', 16, '', 8, 100)], coerce: (v) => ({ ...v, phases: +v.phases }) },
  { id: 'population', icon: '📊', title: 'Population balance resolved on the flow field', generator: 'openfoamPopulationCase', kind: 'flow', solver: 'OpenFOAM reactingTwoPhaseEulerFoam, size-group population balance',
    resolves: 'The particle (or bubble) size distribution transported in every cell of a three-dimensional turbulent flow, with agglomeration and breakage acting between size classes. The in-app population balance is zero- or one-dimensional.',
    fields: [S('dispersed', 'Dispersed phase', 'hydrate', [['hydrate', 'Hydrate particles'], ['gas', 'Gas bubbles']]), N('nGroups', 'Size classes', 12, '', 4, 40), N('fraction', 'Dispersed volume fraction', 0.05, '', 0.001, 0.4), N('lengthD', 'Section length', 20, 'diameters', 3, 500), N('cellsPerDiameter', 'Cells across the diameter', 16, '', 8, 100)] },
  { id: 'erosion', icon: '🪨', title: 'Particle tracking and erosion in a bend', generator: 'openfoamErosionCase', kind: 'flow', solver: 'OpenFOAM DPMFoam with a particle-erosion cloud function',
    resolves: 'Sand (or hydrate) parcels tracked individually through the resolved flow of a bend: where they hit, at which angle and speed, and the eroded-volume field on the wall. The in-app erosion models use one representative impact.',
    fields: [S('particles', 'Particles', 'sand', [['sand', 'Sand (rate and size from the case)'], ['hydrate', 'Hydrate particles']]), S('carrier', 'Carrier fluid', 'mixture', [['mixture', 'Gas–liquid mixture'], ['gas', 'Gas'], ['liquid', 'Liquid']]), N('bendRadiusD', 'Bend radius', 1.5, 'diameters', 1, 50), N('cellsPerDiameter', 'Cells across the diameter', 20, '', 8, 100)] },
  { id: 'fea', icon: '🧱', title: 'Shell or solid finite-element model of a span, bend or jumper', generator: 'calculixPipeCase', kind: 'structure', solver: 'CalculiX (S8R shells or C3D20R solids)',
    resolves: 'The stress field through the wall and around bends and supports, ovalisation, the mode shapes and the response in time to the slug force. The in-app integrity models are beam and thick-cylinder formulae.',
    fields: [S('geometry', 'Geometry', 'span', GEOM), S('element', 'Elements', 'shell', [['shell', 'Shell (S8R)'], ['solid', 'Solid (C3D20R)']]), S('ends', 'End conditions', 'fixed', [['fixed', 'Fixed'], ['pinned', 'Pinned']]), N('spanLength', 'Span length', 15, 'm', 1, 400), N('nC', 'Elements around the circumference', 16, '', 8, 96)] },
  { id: 'fsi', icon: '🔗', title: 'Flow and structure solved together (two-way coupling)', generator: 'fsiCase', kind: 'structure', solver: 'OpenFOAM + CalculiX coupled through preCICE',
    resolves: 'The Navier–Stokes flow and the dynamics of the pipe wall exchange force and displacement every time step: wall loads from the CFD on a finite-element model, and the wall motion back on the flow. The in-app model applies a slug-force formula to a beam.',
    fields: [S('geometry', 'Geometry', 'jumper', GEOM), S('coupling', 'Coupling', 'implicit', [['implicit', 'Implicit, iterated (two-way FSI)'], ['explicit', 'Explicit, one exchange per step (CFD loads on FEA)']]), N('cellsPerDiameter', 'Fluid cells across the diameter', 12, '', 6, 80)] },
  { id: 'gerg', icon: '🧪', title: 'Reference equation of state (GERG-2008)', generator: 'coolpropScript', kind: 'property', solver: 'CoolProp, multi-fluid Helmholtz model',
    resolves: 'Gas density, compressibility factor, heat capacity, Joule–Thomson coefficient and speed of sound from the GERG-2008 reference model on the grid of the case property table, to cross-check or replace the cubic-EOS values.', fields: [] },
  { id: 'saft', icon: '⚗️', title: 'Molecular-based equation of state (PC-SAFT, CPA)', generator: 'teqpScript', kind: 'property', solver: 'teqp',
    resolves: 'Densities, heat capacities and Joule–Thomson coefficients of both phases from PC-SAFT at the equilibrium compositions of the case; with CPA, liquid water with association.',
    fields: [S('model', 'Model', 'PCSAFT', [['PCSAFT', 'PC-SAFT'], ['CPA', 'PC-SAFT + CPA for water']])] },
];

const ctxOf = () => ({ fluid: store.case.fluid, site: store.case.site, outputs: store.case.outputs, inputs: store.case.inputs, name: store.case.name });
const itemsOf = (generator) => Object.entries(B.HANDOFF).flatMap(([suite, list]) => list.filter((e) => e.generator === generator).map((e) => ({ suite, item: e.item })));
async function readFiles(fileList) {
  const out = [];
  for (const f of [...fileList].slice(0, 40)) { checkFile(f); out.push({ name: f.name, text: await f.text() }); }
  return out;
}
/** Button that opens a file chooser accepting several files; handler receives [{ name, text }]. */
function importMany(label, handler, accept) {
  const inp = h('input', { type: 'file', accept, multiple: true, hidden: true });
  inp.addEventListener('change', async () => {
    const files = inp.files; if (!files || !files.length) return;
    try { await handler(await readFiles(files)); } catch (e) { toast(e.message || 'Could not read the file.', 'bad', 8000); }
    inp.value = '';
  });
  return h('span', null, inp, btn('⬆  Import results', () => inp.click(), '', 'Accepted: ' + accept));
}
const ratioStatus = (r) => (r === null ? '' : r > 0.8 && r < 1.25 ? 'ok' : r > 0.5 && r < 2 ? 'warn' : 'bad');
const comparisonView = (rows, title) => [
  kpiGrid(rows.slice(0, 8).map((r) => ({ label: r.quantity, value: r.external, unit: r.unit, status: ratioStatus(r.ratio), help: r.inApp === null ? 'No in-app value to compare with.' : `In-app: ${fmt(r.inApp)} ${r.unit} — ratio ${fmt(r.ratio, 3)}` }))),
  dataTable({ title, columns: ['Quantity', 'Unit', 'External solver', 'In-app / hand value', 'Ratio', 'Difference'], rows: rows.map((r) => [r.quantity, r.unit, r.external, r.inApp ?? '–', r.ratio ?? '–', r.difference ?? '–']), note: 'Ratio = external solver ÷ in-app value. Green tiles agree within 25 %, amber within a factor of two.' }),
];
const thin = (t, y, n = 600) => { if (t.length <= n) return { x: t, y }; const k = Math.ceil(t.length / n), x = [], yy = []; for (let i = 0; i < t.length; i += k) { x.push(t[i]); yy.push(y[i]); } return { x, y: yy }; };
const seriesPlot = (title, ylabel, list, hline) => (list.length ? pc({ type: 'line', title, xlabel: 'Time (s)', ylabel, series: list.slice(0, 8).map((s) => ({ name: s.name, ...thin(s.t, s.y) })), hlines: hline ? [hline] : [] }) : null);

function saveRecord(id, record, x) {
  const prev = store.case.outputs.bridge || {}, imports = { ...(prev.imports || {}) };
  imports[id] = record;
  store.setOutputs('bridge', { imports, count: Object.keys(imports).length, x });
}

/** Result view of an imported OpenFOAM case. */
function flowResults(def, res, ref) {
  const m = res.metrics, rows = B.compareFlow(m, ref), by = (re) => res.series.filter((s) => re.test(s.name));
  const kp = [['Simulated time', m.duration, 's'], ['Samples', m.samples, ''], ['Holdup swing (min – max)', m.holdupMin !== null ? `${fmt(m.holdupMin, 3)} – ${fmt(m.holdupMax, 3)}` : null, ''], ['Mean force on the wall', m.forceMean, 'N'], ['Eroded volume, whole wall', m.erodedVolume, 'm³']].filter((k) => k[1] !== null && k[1] !== undefined).map(([label, value, unit]) => ({ label, value, unit }));
  return h('div', null, h('h3', null, 'Imported results'), h('p', { class: 'summary' }, `${res.files} file${res.files === 1 ? '' : 's'}, ${res.series.length} time series. The first 30 % of the simulated time is left out of the averages as start-up.`),
    rows.length ? comparisonView(rows, `${def.title} — external solver beside the in-app values`) : h('p', { class: 'note' }, 'The files hold no quantity that has an in-app counterpart; the series are plotted below.'), kp.length ? kpiGrid(kp) : null,
    h('div', { class: 'plots' }, seriesPlot('Liquid holdup in time', 'Liquid holdup (–)', by(/^(holdupVolume|section\d|outletHoldup): .*alpha\.liquid/), ref.holdup ? { y: ref.holdup, label: 'in-app holdup', color: '#64748b' } : null), seriesPlot('Pressure on the sections', 'Pressure (Pa)', by(/^section\d: .*\(p\)|^(inlet|outlet)Pressure/)), seriesPlot('Force on the wall', 'Force (N)', by(/^forces: (\|total\||total_[xyz])/)), seriesPlot('Erosion', 'Eroded volume (m³)', by(/^erosion/)), seriesPlot('Mean particle diameter', 'Diameter (m)', by(/\(d\.[\w]+\)/))));
}
/** Result view of imported CalculiX files (.dat and/or .frd, optionally with OpenFOAM force files of a coupled run). */
function structureResults(def, files, ref) {
  const M = {}, plots = [], notes = [];
  let foam = [];
  for (const f of files) {
    if (/\.frd$/i.test(f.name)) { const r = B.importCalculixFrd(f.text); for (const k of ['maxDisplacement', 'maxMises', 'f1']) if (r.metrics[k] !== null) M[k] = Math.max(M[k] ?? -Infinity, r.metrics[k]); const d = r.steps.filter((s) => s.kind === 'DISP'), s = r.steps.filter((x) => x.kind === 'STRESS'); notes.push(`${f.name}: ${r.nodes.toLocaleString()} nodes, ${r.steps.length} result blocks`);
      if (!r.frequencies.length && d.length > 2) plots.push(pc({ type: 'line', title: `Largest displacement in the model — ${f.name}`, xlabel: 'Time (s)', ylabel: 'Displacement (mm)', series: [{ name: 'max |u|', x: d.map((x) => x.time), y: d.map((x) => x.max * 1000), mode: 'both' }] }));
      if (!r.frequencies.length && s.length > 2) plots.push(pc({ type: 'line', title: `Largest von Mises stress in the model — ${f.name}`, xlabel: 'Time (s)', ylabel: 'Stress (MPa)', series: [{ name: 'max von Mises', x: s.map((x) => x.time), y: s.map((x) => x.max / 1e6), mode: 'both' }], hlines: ref.smys ? [{ y: ref.smys / 1e6, label: 'SMYS', color: '#dc2626' }] : [] })); }
    else if (/E I G E N V A L U E|displacements \(vx|total force|stresses \(elem/.test(f.text.slice(0, 200000)) || /\.dat$/i.test(f.name) && !/^#/.test(f.text.trim())) { const r = B.importCalculixDat(f.text); if (r.metrics.f1 !== null) M.f1 = r.metrics.f1; if (r.metrics.maxDisplacement !== null) M.maxDisplacement = Math.max(M.maxDisplacement ?? 0, r.metrics.maxDisplacement); if (r.metrics.peakReaction !== null) M.peakReaction = r.metrics.peakReaction; notes.push(`${f.name}: ${r.frequencies.length} natural frequencies, ${r.displacements.length} displacement blocks`);
      if (r.frequencies.length) plots.push(pc({ type: 'bar', title: 'Natural frequencies', ylabel: 'Frequency (Hz)', categories: r.frequencies.slice(0, 20).map((_, i) => 'Mode ' + (i + 1)), series: [{ name: 'CalculiX', values: r.frequencies.slice(0, 20) }] }));
      if (r.displacements.length > 2) plots.push(pc({ type: 'line', title: `Displacement of node set ${r.displacements[0].set}`, xlabel: 'Time (s)', ylabel: 'Displacement (mm)', series: [{ name: 'largest in the set', x: r.displacements.map((d) => d.time), y: r.displacements.map((d) => d.max * 1000) }] })); }
    else foam.push(f);
  }
  let flow = null;
  if (foam.length) { flow = B.importOpenfoamPostProcessing(foam, ref); if (flow.metrics.forcePeak !== null) M.forcePeak = flow.metrics.forcePeak; plots.push(seriesPlot('Force of the fluid on the wall', 'Force (N)', flow.series.filter((s) => /^forces: (\|total\||total_[xyz])/.test(s.name)))); }
  if (!Object.keys(M).length) throw new Error('No CalculiX result was recognised. Choose the .dat and .frd files that ccx wrote (and force.dat of a coupled run).');
  const rows = [...B.compareStructure(M, ref), ...(flow ? B.compareResults([['Peak fluid force on the wall', 'N', flow.metrics.forcePeak, ref.forceN]]) : [])];
  return { view: h('div', null, h('h3', null, 'Imported results'), h('p', { class: 'summary' }, notes.join(' · ')), comparisonView(rows, `${def.title} — external solver beside the hand values`), h('div', { class: 'plots' }, plots)), metrics: M, rows };
}

/** One hand-off card. */
function handoffCard(def, state) {
  const values = Object.fromEntries(def.fields.map((f) => [f.key, f.value])), est = h('div'), cmdBox = h('div'), resBox = h('div', { role: 'region', 'aria-live': 'polite' }), items = itemsOf(def.generator);
  let built = null, table = null;
  const optsNow = () => {
    const v = def.coerce ? def.coerce({ ...values }) : { ...values };
    if (def.kind === 'property') { const fm = fluidModel(ctxOf()); table = fm.table; return { ...v, grid: B.propertyGrid(fm.spec, { P: fm.table.P, T: fm.table.T }), aqueous: { inhibitor: fm.spec.inhibitor, wt: fm.spec.inhWt } }; }
    return { ...v, case: state.snapshot() };
  };
  const build = () => { built = B[def.generator](optsNow()); return built; };
  const paint = () => {
    try {
      const r = build(), p = r.plan || {};
      fill(est, kpiGrid([p.cells ? { label: def.kind === 'structure' ? 'Elements / cells' : def.kind === 'property' ? 'Grid points' : 'Cells', value: p.cells, sig: 6 } : null, p.steps ? { label: 'Time steps', value: p.steps, sig: 6 } : null, p.coreHours !== undefined ? { label: 'Estimated cost', value: p.coreHours < 1 ? p.coreHours * 60 : p.coreHours, unit: p.coreHours < 1 ? 'core-min' : 'core-h', status: p.coreHours > 5000 ? 'bad' : p.coreHours > 100 ? 'warn' : 'ok', help: 'Order of magnitude: cells × time steps × cost per cell and step. Time the first steps of the run and rescale.' } : null, { label: 'Files in the case', value: r.files.length }].filter(Boolean)), h('p', { class: 'summary' }, r.summary));
      fill(cmdBox, h('details', { class: 'ref', open: true }, h('summary', null, 'Commands to run'), h('div', { class: 'tbl-scroll', tabindex: '0' }, h('pre', null, r.commands.join('\n')))), h('details', { class: 'ref' }, h('summary', null, 'Read-me of the case (also inside the .zip)'), h('div', { class: 'tbl-scroll', tabindex: '0' }, h('pre', null, r.readme))));
    } catch (e) { built = null; fill(est, h('p', { class: 'note' }, 'This case cannot be written with the present inputs: ' + (e.message || e))); fill(cmdBox); }
  };
  state.listeners.push(paint);
  const get = (k) => values[k.slice(def.id.length + 1)], set = (k, v) => { values[k.slice(def.id.length + 1)] = v; paint(); };
  const fieldsEl = def.fields.length ? h('div', { class: 'fields' }, def.fields.map((f) => fieldRow({ ...f, key: def.id + '_' + f.key }, get, set, {}))) : null;
  const dl = btn('⬇  Download case (.zip)', () => {
    try { const r = build(); download(new Blob([B.zipStore(r.files, { root: r.name })], { type: 'application/zip' }), r.name + '.zip', 'application/zip'); toast(`${r.name}.zip written: ${r.files.length} files.`, 'ok'); }
    catch (e) { toast(e.message || 'Could not write the case.', 'bad', 7000); }
  }, 'primary', 'A complete, ready-to-run case as one archive');
  const onImport = async (files) => {
    const r = built || build(), ref = { ...r.reference };
    if (def.kind === 'flow') {
      if (def.id === 'erosion' && typeof store.case.outputs.integ?.erosionRate === 'number') ref.erosionRateMmY = store.case.outputs.integ.erosionRate;
      const res = B.importOpenfoamPostProcessing(files, ref); fill(resBox, flowResults(def, res, ref));
      saveRecord(def.id, B.bridgeRecord(def.generator, def.solver, res.metrics, B.compareFlow(res.metrics, ref), `x = ${fmt(state.x)} m`), state.x);
    } else if (def.kind === 'structure') {
      const o = store.case.outputs.integ || {}; if (typeof o.vmUtil === 'number') ref.vmUtil = o.vmUtil;
      const s = structureResults(def, files, ref); fill(resBox, s.view); saveRecord(def.id, B.bridgeRecord(def.generator, def.solver, s.metrics, s.rows, `x = ${fmt(state.x)} m`), state.x);
    } else {
      const fm = fluidModel(ctxOf()), imp = B.importPropertyTable(files[0].text, fm.table), dev = B.comparePropertyTables(fm.table, imp.table), jMid = Math.floor(fm.table.T.length / 2), has = (f) => imp.replaced.includes(f);
      const canPublish = !!store.case.outputs.pvt, pub = btn('Use these values in the case property table', () => { store.setOutputs('pvt', { ...store.case.outputs.pvt, table: imp.table }); toast('Property table of the case replaced: ' + imp.replaced.join(', ') + '. Run the downstream suites again.', 'ok', 7000); }, 'ghost');
      pub.disabled = !canPublish;
      fill(resBox, h('h3', null, 'Imported property table'), h('p', { class: 'summary' }, `${imp.source || files[0].name}: ${imp.filled.toLocaleString()} values for ${imp.replaced.join(', ')}; ${imp.missing.toLocaleString()} cells kept from the case table.`),
        kpiGrid(dev.map((d) => ({ label: `${d.field}: mean deviation of the in-app table`, value: d.meanPct, unit: '%', status: d.meanPct < 2 ? 'ok' : d.meanPct < 10 ? 'warn' : 'bad', help: `Largest ${fmt(d.maxPct, 3)} %, bias ${fmt(d.bias, 3)} % over ${d.n} grid points` }))),
        dataTable({ title: 'In-app property table against the imported model', columns: ['Property', 'Grid points', 'Mean |deviation| %', 'Largest %', 'Bias %'], rows: dev.map((d) => [d.field, d.n, d.meanPct, d.maxPct, d.bias]), note: 'Deviation = (in-app − imported) ÷ imported, over the grid points where the phase exists.' }),
        h('div', { class: 'plots' }, ['rhoG', 'zG', 'rhoO'].filter(has).map((f) => pc({ type: 'line', title: `${f} along the ${fmt(fm.table.T[jMid], 3)} °C isotherm`, xlabel: 'Pressure (bara)', ylabel: f, logx: true, series: [{ name: 'in-app (cubic EOS)', x: fm.table.P, y: fm.table.P.map((_, i) => fm.table[f][i][jMid]) }, { name: 'imported', x: fm.table.P, y: fm.table.P.map((_, i) => imp.table[f][i][jMid]), mode: 'both', dash: true }] }))),
        h('div', { class: 'row-tools' }, pub, canPublish ? null : h('span', { class: 'note' }, 'Run the PVT suite once to publish a property table; the imported values can then replace its columns.')));
      saveRecord(def.id, B.bridgeRecord(def.generator, def.solver, Object.fromEntries(dev.map((d) => [d.field + 'MeanPct', d.meanPct])), dev, imp.source), state.x);
    }
    toast('Results imported and stored with the case.', 'ok');
  };
  const accept = def.kind === 'property' ? '.csv,.txt' : def.kind === 'structure' ? '.dat,.frd,.txt' : '.txt,.dat,.csv';
  const el = card(h('h2', { id: 'bridge-' + def.id, tabindex: '-1' }, def.icon + '  ' + def.title), h('p', null, h('b', null, 'What it resolves: '), def.resolves), h('p', { class: 'note' }, h('b', null, 'Solver: '), def.solver),
    items.length ? h('div', { class: 'chips', 'aria-label': 'Catalogue items covered by this hand-off' }, items.map((i) => h('a', { class: 'chip via', href: '#/suite/' + i.suite, title: 'Catalogue item handed off from suite ' + i.suite }, i.item))) : null,
    fieldsEl, est, h('div', { class: 'row-tools' }, dl, importMany('⬆  Import results', onImport, accept)), cmdBox, resBox);
  el.dataset.generator = def.generator;
  paint();
  return el;
}

/** Page "External solvers". */
export function bridgePage(root) {
  const ctx = ctxOf(), state = { x: 0, listeners: [], snapshot: () => B.bridgeCase(ctxOf(), { x: state.x }) };
  let pic;
  try { pic = flowPicture(ctx); } catch (e) { fill(root, h('header', { class: 'page-head' }, h('h1', null, 'External solvers')), card(h('p', null, 'The case has no steady flow solution yet, so there is nothing to hand off: ' + (e.message || e)))); return; }
  const x0 = pic.x[0], x1 = pic.x[pic.x.length - 1], prev = store.case.outputs.bridge?.x;
  state.x = typeof prev === 'number' && prev >= x0 && prev <= x1 ? prev : Math.round(0.5 * (x0 + x1));
  const kp = h('div'), plotBox = h('div', { class: 'plots' });
  const slider = h('input', { type: 'range', id: 'bridge-x', min: x0, max: x1, step: Math.max(1, Math.round((x1 - x0) / 400)), value: state.x, 'aria-label': 'Location along the line in metres' }), num = h('input', { type: 'number', id: 'bridge-xn', min: x0, max: x1, step: 'any', value: state.x, inputmode: 'decimal' });
  const paintLoc = () => {
    const c = state.snapshot();
    fill(kp, kpiGrid([{ label: 'Distance', value: c.x, unit: 'm', sig: 6 }, { label: 'Elevation', value: c.z, unit: 'm' }, { label: 'Inclination', value: c.angleDeg, unit: '°' }, { label: 'Pressure', value: c.P, unit: 'bara' }, { label: 'Temperature', value: c.T, unit: '°C' }, { label: 'Liquid holdup', value: c.holdup }, { label: 'Liquid / gas superficial velocity', value: `${fmt(c.vsl, 3)} / ${fmt(c.vsg, 3)}`, unit: 'm/s' }, { label: 'Flow regime', value: c.regime || '–' }]));
    fill(plotBox, pc({ type: 'line', title: 'Line profile and the hand-off location', xlabel: 'Distance (m)', ylabel: 'Elevation (m)', height: 220, series: [{ name: 'Elevation', x: pic.x, y: pic.z }], vlines: [{ x: state.x, label: 'hand-off' }] }), pc({ type: 'line', title: 'Liquid holdup along the line', xlabel: 'Distance (m)', ylabel: 'Holdup (–)', height: 220, series: [{ name: 'Holdup', x: pic.x, y: pic.holdup }], vlines: [{ x: state.x, label: fmt(interp1(pic.x, pic.holdup, state.x), 3) }] }));
  };
  let timer = 0;
  const moved = (v, from) => {
    if (!Number.isFinite(v)) return;
    state.x = Math.min(x1, Math.max(x0, v)); if (from !== slider) slider.value = state.x; if (from !== num) num.value = state.x;
    clearTimeout(timer); timer = setTimeout(() => { paintLoc(); state.listeners.forEach((fn) => fn()); }, 180);
  };
  slider.addEventListener('input', () => moved(+slider.value, slider)); num.addEventListener('change', () => moved(+num.value, num));
  const jump = (label, v) => btn(label, () => moved(v), 'mini');
  const maxAt = (key) => pic.x[pic[key].indexOf(Math.max(...pic[key]))], stored = store.case.outputs.bridge?.imports || {};
  const cards = CARDS.map((def) => handoffCard(def, state));
  fill(root,
    h('header', { class: 'page-head' }, h('h1', null, 'External solvers'), h('p', null, 'The engines of this application are one-dimensional and lumped, which is what lets them run on any device in seconds. Formulations that need a three-dimensional grid, a finite-element mesh or a reference property model are handed off here: the application writes a complete, ready-to-run case for an established open-source solver from the data the case already holds, and reads the results back so that they join the rest of the study.')),
    card(h('h2', null, 'Where along the line'), h('p', { class: 'note' }, `Every case below is written for this location, with pressure, temperature, velocities, holdup and fluid properties taken from the ${pic.source}.`),
      h('div', { class: 'fields' }, h('div', { class: 'field wide' }, h('label', { for: 'bridge-x' }, 'Location along the line'), slider), h('div', { class: 'field' }, h('label', { for: 'bridge-xn' }, 'Distance from the inlet'), h('div', { class: 'ctl' }, h('div', { class: 'inp' }, num, h('span', { class: 'unit' }, 'm'))))),
      h('div', { class: 'row-tools' }, jump('Inlet', x0), jump('Highest holdup', maxAt('holdup')), jump('Highest velocity', maxAt('vm')), jump('Steepest pressure gradient', maxAt('dpdx')), jump('Outlet', x1)), kp, plotBox),
    Object.keys(stored).length ? card(h('h2', null, 'Results already imported for this case'), h('div', { class: 'row-tools' }, Object.entries(stored).map(([id, r]) => badge(`${CARDS.find((c) => c.id === id)?.title || id}: ${r.solver}`, 'ok'))), h('p', { class: 'note' }, 'They are stored with the case and cited on the decision page. Import again to replace them.')) : null,
    card(h('h2', null, 'How a hand-off works'), h('ol', { class: 'steps' }, h('li', null, 'Choose the location and the options of a card; the size of the mesh and the cost of the run are estimated as you type.'), h('li', null, 'Download the case: one .zip with every input file, the run script and a read-me.'), h('li', null, 'Run it with the solver named on the card, on a workstation or a cluster. The solvers are free and open source; nothing is installed or sent by this application.'), h('li', null, 'Import the result file the run produced. The values appear beside the in-app results and are stored with the case.'))),
    cards);
  paintLoc();
  // deep link from the Equations tab: #/bridge?g=<generator>
  const want = /[?&]g=([\w]+)/.exec(globalThis.location?.hash || '')?.[1];
  if (want) { const el = cards.find((c) => c.dataset.generator === want); if (el) setTimeout(() => { el.scrollIntoView({ block: 'start' }); el.querySelector('h2')?.focus(); }, 60); }
}
export { CARDS as BRIDGE_CARDS };
