// Non-suite pages: home, case & fluid, global site data, data portal, integrated run, app & offline.
import { h, clear, btn, kpiGrid, dataTable, toast, badge, fieldRow, importBtn, help, fill } from '../core/ui.js';
import { store } from '../core/store.js';
import { SUITES, CHAIN, byId, loadSuite, downstream } from '../suites/index.js';
import { runSuite, applyLinks, linkItems, allFields, setInputValue } from '../core/suiteview.js';
import { readFiles, geometryCard, generatorPanel, formatCatalogue, attachGeometry, derive, ACCEPT } from '../core/geomview.js';
import { geometryLinks } from '../core/geomlinks.js';
import { SUITE_GEOMETRY, formatOf } from '../core/geom.js';
import { fetchSite, searchPlace, SOURCES, loadAtlas, mergeSiteData, ATLAS_LABELS } from '../core/live.js';
import { plotCard } from '../core/plot.js';
import { fmt } from '../core/num.js';
import { makeFluid, streams, saturationP, hydrateT, hydrateDepression, aqueous, INHIBITORS, EOS, DEFAULT_FLUID } from '../core/thermo.js';
import { download, readTable, extOf, checkFile } from '../core/io.js';
import { MIRRORS, APP } from '../data/app.js';
import { LAND } from '../data/atlas.js';
import { EXAMPLES } from '../data/examples.js';

const card = (...kids) => h('section', { class: 'card' }, ...kids);
const pc = (spec) => plotCard(spec, { onDownload: download });
const ago = (iso) => { if (!iso) return 'never'; const s = (Date.now() - new Date(iso)) / 1000; return s < 90 ? 'just now' : s < 5400 ? Math.round(s / 60) + ' min ago' : s < 172800 ? Math.round(s / 3600) + ' h ago' : Math.round(s / 86400) + ' d ago'; };
const N = SUITES.length;
const rateText = (f) => (f.rateBasis === 'gas' ? `${fmt(f.qGas)} MSm³/d gas` : f.rateBasis === 'mass' ? `${fmt(f.mdot)} kg/s` : `${fmt(f.qOil)} Sm³/d oil`);

// ---------------------------------------------------------------------------------------------- home
export function home(root) {
  const c = store.case, done = SUITES.filter((s) => c.outputs[s.id]).length;
  const cards = SUITES.map((s) => {
    const o = c.outputs[s.id], k = o?._kpis?.slice(0, 3) || [];
    return h('a', { class: 'suite-card' + (o ? ' done' : ''), href: '#/suite/' + s.id },
      h('div', { class: 'sc-top' }, h('span', { class: 'sc-ico', 'aria-hidden': 'true' }, s.icon), h('span', { class: 'sc-num' }, String(s.num).padStart(2, '0')), o ? badge('solved ' + ago(o._at), 'ok') : badge('not run', '')),
      h('h3', null, s.title), h('p', null, s.blurb),
      k.length ? h('ul', { class: 'sc-kpi' }, k.map((q) => h('li', null, q.label + ': ', h('b', null, typeof q.value === 'number' ? fmt(q.value, 3) : String(q.value)), ' ' + q.unit))) : null,
      h('div', { class: 'sc-links' }, s.uses.length ? '⛓ uses ' + s.uses.map((u) => byId(u).short).join(', ') : '⛓ starts the chain'));
  });
  const next = !c.site.fetchedAt && done === 0 ? ['Start with the fluid and rates of your case, or load a worked example.', '#/case', 'Define the case'] : done < N ? [`${done} of ${N} suites solved. Run the whole chain in one pass, with the couplings between suites iterated.`, '#/chain', 'Run everything'] : ['Every suite is solved. Read the ranked recommendations for this case.', '#/advisor', 'See the decisions'];
  fill(root,
    h('section', { class: 'hero' },
      h('div', null, h('h1', null, 'One connected workspace for multiphase flow assurance'),
        h('p', null, 'Seven engineering suites share one case: fluid and phase behaviour; geometry, wells, network and equipment; multiphase thermal-hydraulics and slugging; hydrates and other solids; operations and control; integrity, loads and risk; economics and decisions. Results pass from suite to suite — forwards and back — so one industrial case is analysed from the reservoir fluid to the net present value.'),
        h('div', { class: 'row-tools' }, h('a', { class: 'btn primary', href: '#/case' }, '1 · Define the case'), h('a', { class: 'btn', href: '#/site' }, '2 · Pull site data'), h('a', { class: 'btn', href: '#/portal' }, '3 · Bring geometry & data'), h('a', { class: 'btn', href: '#/chain' }, '4 · Run everything'), h('a', { class: 'btn', href: '#/advisor' }, '5 · Decide')),
        h('p', { class: 'note nextstep' }, h('b', null, 'Suggested next step: '), next[0] + ' ', h('a', { href: next[1] }, next[2] + ' →'))),
      h('div', { class: 'hero-stat' }, kpiGrid([{ label: 'Active case', value: c.name }, { label: 'Fluid', value: c.fluid.name || 'custom' }, { label: 'Rate', value: rateText(c.fluid) }, { label: 'Site', value: c.site.name || (c.site.lat !== null ? `${fmt(c.site.lat, 4)}, ${fmt(c.site.lon, 4)}` : 'not set') }, { label: 'Suites solved', value: `${done} / ${N}` }, { label: 'Site data', value: c.site.fetchedAt ? ago(c.site.fetchedAt) : 'not pulled' }]))),
    h('h2', { class: 'sect' }, `The ${N} suites`), h('div', { class: 'suite-grid' }, cards,
      [['🌍', 'Global site data', 'Live seabed depth and temperature, sea state, weather, prices and national figures for any location', '#/site', c.site.fetchedAt ? 'pulled ' + ago(c.site.fetchedAt) : 'not pulled'], ['📥', 'Data portal', 'Routes, well surveys, networks, CAD, meshes, terrain, inspection maps, logs and tables in one place', '#/portal', 'import'], ['🔗', 'Integrated run', 'Solve all suites in sequence and iterate the couplings between them', '#/chain', `${done} / ${N} solved`], ['🧭', 'Decision support', 'Ranked recommendations, sustainability scorecard and live evidence', '#/advisor', 'whole case']].map(([ico, t, b, href, st]) =>
        h('a', { class: 'suite-card whole', href }, h('div', { class: 'sc-top' }, h('span', { class: 'sc-ico', 'aria-hidden': 'true' }, ico), h('span', { class: 'sc-num' }, ''), badge(st, '')), h('h3', null, t), h('p', null, b), h('div', { class: 'sc-links' }, 'whole case')))),
    h('h2', { class: 'sect' }, 'How the suites are wired together'),
    card(h('div', { class: 'flowmap' }, CHAIN.map((id) => { const s = byId(id), d = downstream(id); return h('div', { class: 'fm-row' }, h('a', { class: 'chip on', href: '#/suite/' + id }, `${s.num}. ${s.short}`), h('span', { class: 'fm-arrow', 'aria-hidden': 'true' }, '→'), h('div', { class: 'fm-to' }, d.length ? d.map((x) => h('a', { class: 'chip', href: '#/suite/' + x.id }, `${x.num}. ${x.short}`)) : h('span', { class: 'note' }, 'final results — feeds the decision report'))); })),
      h('p', { class: 'note' }, 'Each row reads “this suite feeds →”. The fluid model supplies properties and the hydrate curve to everything; the network supplies the route, wall and equipment; the flow solution gives pressure, temperature, hold-up and slugs; solids narrow the bore and roughen the wall, which is sent back to the flow solution; operations change rates, chemicals and valves; integrity turns the loads into damage and risk; economics turns all of it into cost, value and a ranked decision.')));
}

// ------------------------------------------------------------------------------------- case & fluid
export function casePage(root) {
  const c = store.case, sumBox = h('div'), warnBox = h('div');
  let timer = 0;
  const paintSum = () => {
    clearTimeout(timer);
    timer = setTimeout(() => {
      try {
        const f = store.case.fluid, m = makeFluid(f), st = streams(f, m), sat = saturationP(m, f.Tres), aq = aqueous(f), sg = st.std.gasSG ?? m.MW / 28.9647, tot = Object.values(f.comp).reduce((a, b) => a + (+b || 0), 0);
        fill(sumBox, kpiGrid([
          { label: 'Gas–oil ratio', value: Number.isFinite(st.gor) ? st.gor : '∞ (no stock-tank liquid)', unit: Number.isFinite(st.gor) ? 'Sm³/Sm³' : '', help: 'Single-stage flash of the well stream to 1.01325 bara and 15 °C.' },
          st.std.api !== null ? { label: 'Stock-tank oil', value: st.std.api, unit: '°API' } : null, { label: 'Gas gravity', value: sg, unit: '(air = 1)' },
          sat.P ? { label: sat.type === 'dew' ? 'Dew point at reservoir T' : 'Bubble point at reservoir T', value: sat.P, unit: 'bara', status: sat.P > f.Pres ? 'warn' : 'ok', help: sat.P > f.Pres ? 'The fluid is already two-phase at reservoir conditions.' : 'The fluid is single-phase in the reservoir.' } : { label: 'Saturation point', value: 'single phase' },
          { label: 'Oil', value: st.qOilStd, unit: 'Sm³/d' }, { label: 'Gas', value: st.qGasStd / 1e6, unit: 'MSm³/d' }, { label: 'Water', value: st.qWaterStd, unit: 'Sm³/d' }, { label: 'Total mass rate', value: st.mHC + st.mW, unit: 'kg/s' },
          { label: 'Hydrate temperature at arrival pressure', value: hydrateT(f.Pout, sg, aq), unit: '°C', help: 'Screening value (gas-gravity correlation with the salt and inhibitor of the case). Suite 1 computes the rigorous curve.' },
          { label: 'Hydrate depression by salt + inhibitor', value: hydrateDepression(aq), unit: '°C' },
        ].filter(Boolean)));
        const w = [];
        if (Math.abs(tot - 100) > 0.05) w.push(`The composition adds up to ${fmt(tot, 5)} mol %; it is normalised to 100 % in every calculation.`);
        if (f.rateBasis === 'oil' && !(st.std.vOil > 0)) w.push('This fluid leaves no liquid at standard conditions, so an oil rate cannot define it: the gas rate is used instead. Switch the rate basis to gas.');
        if (f.Tin < hydrateT(Math.max(f.Pout, 50), sg, aq)) w.push('The inlet temperature is already below the hydrate temperature at line pressure: hydrates can form from the first metre.');
        fill(warnBox, w.length ? h('div', { class: 'warns' }, w.map((t) => h('div', { class: 'warn' }, h('b', null, 'Check'), ' ' + t))) : null);
      } catch (e) { fill(sumBox, h('p', { class: 'bad' }, 'This fluid cannot be characterised: ' + (e.message || e))); }
    }, 220);
  };
  const get = (k) => store.case.fluid[k];
  const set = (k, v, lib) => { if (k === 'comp') store.setFluid(lib ? { comp: v, c7MW: lib.c7MW, c7SG: lib.c7SG, name: lib.name, ...(lib.rateBasis ? { rateBasis: lib.rateBasis, qGas: lib.qGas, qWater: lib.qWater } : { rateBasis: 'oil' }) } : { comp: v }); else store.setFluid({ [k]: v }); paintSum(); if (lib || k === 'rateBasis' || k === 'inhibitor') casePage(root); };
  const F = (def) => fieldRow({ ...def, value: store.case.fluid[def.key] }, get, set, {});
  const f = c.fluid, basis = f.rateBasis;
  const lib = store.library(), libSel = h('select', { 'aria-label': 'Saved cases' }, h('option', { value: '' }, Object.keys(lib).length ? 'Open a saved case…' : 'No saved cases yet'), Object.keys(lib).map((n) => h('option', { value: n }, n)));
  libSel.addEventListener('change', () => { if (libSel.value) { store.loadFromLibrary(libSel.value); toast('Case loaded.', 'ok'); casePage(root); } });
  const exSel = h('select', { 'aria-label': 'Worked examples' }, h('option', { value: '' }, 'Load a worked example…'), EXAMPLES.map((e, i) => h('option', { value: i }, e.name)));
  exSel.addEventListener('change', () => { const e = EXAMPLES[+exSel.value]; if (e && (Object.keys(store.case.inputs).length === 0 || confirm('Replace the current case with this worked example? Unsaved inputs of the current case are discarded.'))) { store.importJSON(JSON.stringify(e.case)); toast(`Loaded: ${e.name}.`, 'ok'); casePage(root); } exSel.value = ''; });
  paintSum();
  fill(root,
    h('header', { class: 'page-head' }, h('h1', null, 'Case & fluid'), h('p', null, 'A case is one industrial study: its site, the produced fluid and rates, and the inputs and results of every suite. Everything is stored on this device only.')),
    card(h('h2', null, 'Case'),
      h('div', { class: 'fields' },
        h('div', { class: 'field' }, h('label', { for: 'c_name' }, 'Case name'), h('div', { class: 'ctl' }, h('input', { id: 'c_name', type: 'text', maxlength: 120, value: c.name, oninput: (e) => store.update({ name: e.target.value.slice(0, 120) || 'Untitled case' }) }))),
        h('div', { class: 'field' }, h('label', { for: 'c_auto' }, 'Auto-link suites', help('When on, each suite takes matching values from the case fluid, the site data and the other suites every time it runs.')), h('div', { class: 'ctl' }, h('label', { class: 'switch' }, h('input', { id: 'c_auto', type: 'checkbox', checked: c.autolink, onchange: (e) => store.update({ autolink: e.target.checked }) }), h('span', { class: 'slider' })))),
        h('div', { class: 'field wide' }, h('label', { for: 'c_notes' }, 'Notes'), h('div', { class: 'ctl wide' }, h('textarea', { id: 'c_notes', rows: 2, maxlength: 5000, oninput: (e) => store.update({ notes: e.target.value }) }, c.notes)))),
      h('div', { class: 'row-tools' },
        btn('Save to this device', () => toast(store.saveToLibrary() ? 'Case saved in the library on this device.' : 'Could not save — device storage is full.', 'ok'), 'primary'), libSel, exSel,
        btn('Export case file', () => download(store.exportJSON(), `${store.case.name}${APP.caseExt}`, 'application/json'), '', 'One portable JSON file with every input — share it or open it on another device'),
        importBtn('Import case file', async (file) => { checkFile(file); store.importJSON(await file.text()); toast('Case imported.', 'ok'); casePage(root); }, '.json'),
        btn('New blank case', () => { if (confirm('Start a new blank case? Unsaved inputs of the current case are discarded.')) { store.reset(); casePage(root); } }, 'ghost'))),
    card(h('h2', null, 'Produced fluid'), h('p', { class: 'note' }, 'The water-free well-stream composition. Pick a reference fluid, import a laboratory report (CSV / Excel with component names and mol %), or type the values. Every suite takes its fluid properties from this definition through the equation of state.'),
      h('div', { class: 'fields' },
        F({ key: 'name', label: 'Fluid name', type: 'text' }),
        F({ key: 'comp', label: 'Composition (mol %)', type: 'composition' }),
        F({ key: 'c7MW', label: 'C7+ molar mass', unit: 'g/mol', min: 96, max: 600, help: 'Average molar mass of the heptanes-plus fraction from the laboratory report.' }),
        F({ key: 'c7SG', label: 'C7+ specific gravity', unit: '(water = 1)', min: 0.7, max: 1.05 }),
        F({ key: 'eos', label: 'Equation of state', type: 'select', options: Object.entries(EOS).map(([value, e]) => ({ value, label: e.name })) }),
        F({ key: 'nPseudo', label: 'C7+ pseudo-components', type: 'select', options: [{ value: 1, label: '1 (single lump)' }, { value: 2, label: '2' }, { value: 3, label: '3 (recommended)' }] }))),
    card(h('h2', null, 'Rates, water and chemicals'),
      h('div', { class: 'fields' },
        F({ key: 'rateBasis', label: 'Rate specified as', type: 'select', options: [{ value: 'oil', label: 'Stock-tank oil rate (oil systems)' }, { value: 'gas', label: 'Sales-gas rate (gas and condensate systems)' }, { value: 'mass', label: 'Hydrocarbon mass rate' }] }),
        basis === 'oil' ? F({ key: 'qOil', label: 'Oil rate', unit: 'Sm³/d', min: 0, max: 1e6, help: 'Stock-tank oil. 1 Sm³/d = 6.29 bbl/d.' }) : basis === 'gas' ? F({ key: 'qGas', label: 'Gas rate', unit: 'MSm³/d', min: 0, max: 500, help: 'Million standard cubic metres per day. 1 MSm³/d = 35.3 MMscf/d.' }) : F({ key: 'mdot', label: 'Hydrocarbon mass rate', unit: 'kg/s', min: 0, max: 5000 }),
        basis === 'oil' ? F({ key: 'wc', label: 'Water cut', unit: '%', min: 0, max: 98, help: 'Water as a share of the total liquid at standard conditions.' }) : F({ key: 'qWater', label: 'Water rate', unit: 'Sm³/d', min: 0, max: 1e6 }),
        F({ key: 'salinity', label: 'Water salinity', unit: 'wt % NaCl eq.', min: 0, max: 26, help: 'Seawater is about 3.5 wt %. Salt lowers the hydrate temperature.' }),
        F({ key: 'inhibitor', label: 'Hydrate inhibitor in the water', type: 'select', options: Object.entries(INHIBITORS).map(([value, i]) => ({ value, label: i.name })) }),
        f.inhibitor !== 'none' ? F({ key: 'inhWt', label: 'Inhibitor concentration', unit: 'wt % of aqueous phase', min: 0, max: 95 }) : null)),
    card(h('h2', null, 'Operating conditions'),
      h('div', { class: 'fields' },
        F({ key: 'Tin', label: 'Flowline inlet temperature', unit: '°C', min: -20, max: 200, help: 'Temperature where the fluid enters the flowline (wellhead or manifold).' }),
        F({ key: 'Pout', label: 'Arrival pressure', unit: 'bara', min: 1, max: 600, help: 'Pressure at the receiving separator or terminal.' }),
        F({ key: 'Pres', label: 'Reservoir pressure', unit: 'bara', min: 1, max: 1500 }), F({ key: 'Tres', label: 'Reservoir temperature', unit: '°C', min: 0, max: 250 })),
      h('h3', null, 'What this fluid is'), sumBox, warnBox,
      h('div', { class: 'row-tools' }, h('a', { class: 'btn primary', href: '#/suite/pvt' }, 'Open suite 1 for the full phase behaviour →'), btn('Reset fluid to the reference case', () => { store.setFluid({ ...DEFAULT_FLUID, comp: { ...DEFAULT_FLUID.comp } }); casePage(root); }, 'ghost'))));
}

// ---------------------------------------------------------------------------------- global site data
// Map tiles are fetched (not hot-linked) so the HTTP status can be checked: a refused tile is never
// shown, the next provider is tried instead, and good tiles are kept in memory for the session.
const TILE_SOURCES = [(z, x, y) => `https://tile.openstreetmap.org/${z}/${x}/${y}.png`, (z, x, y) => `https://basemaps.cartocdn.com/rastertiles/voyager/${z}/${x}/${y}.png`];
const tileCache = new Map();
let tileSource = 0;
function loadTile(z, x, y) {
  const key = `${z}/${x}/${y}`;
  if (tileCache.has(key)) return tileCache.get(key);
  const p = (async () => {
    for (let k = tileSource; k < TILE_SOURCES.length; k++) {
      try {
        const r = await fetch(TILE_SOURCES[k](z, x, y), { mode: 'cors', credentials: 'omit', referrerPolicy: 'strict-origin-when-cross-origin' });
        if (r.ok && (r.headers.get('content-type') || '').startsWith('image/')) return URL.createObjectURL(await r.blob());
        if (r.status === 403 || r.status === 429) tileSource = Math.max(tileSource, k + 1); // provider refuses this app: stop asking it
      } catch { /* offline or blocked: try the next provider */ }
    }
    tileCache.delete(key);
    return null;
  })();
  if (tileCache.size > 600) { const first = tileCache.keys().next().value; tileCache.get(first).then((u) => u && URL.revokeObjectURL(u)); tileCache.delete(first); }
  tileCache.set(key, p);
  return p;
}

function slippyMap(lat, lon, onPick) {
  let zoom = lat === null || lat === undefined ? 2 : 8, cLat = lat ?? 22, cLon = lon ?? 30, mLat = lat, mLon = lon;
  const box = h('div', { class: 'map', tabindex: '0', role: 'application', 'aria-label': 'World map. Click or tap to choose the site.' }), layer = h('div', { class: 'map-layer' }), pin = h('div', { class: 'map-pin', hidden: true }, '📍');
  const baseCv = h('canvas', { class: 'map-base', 'aria-hidden': 'true' }); // built-in coastlines: always there, tiles are drawn over them when reachable
  let borders = null, bordersAsked = false; // country outlines of the built-in atlas: fetched only when map tiles cannot be loaded
  const needBorders = () => { if (bordersAsked) return; bordersAsked = true; loadAtlas().then((A) => A.atlasBorders()).then((b) => { borders = b; draw(); }).catch(() => { bordersAsked = false; }); };
  box.append(baseCv, layer, pin, h('div', { class: 'map-zoom' }, h('button', { type: 'button', 'aria-label': 'Zoom in', onclick: (e) => { e.stopPropagation(); zoom = Math.min(15, zoom + 1); draw(); } }, '+'), h('button', { type: 'button', 'aria-label': 'Zoom out', onclick: (e) => { e.stopPropagation(); zoom = Math.max(2, zoom - 1); draw(); } }, '−')),
    h('div', { class: 'map-attr' }, 'Coastlines and borders: Natural Earth · tiles © OpenStreetMap contributors · © CARTO'));
  const X = (lo, z) => ((lo + 180) / 360) * 2 ** z * 256, Y = (la, z) => { const r = (la * Math.PI) / 180; return ((1 - Math.log(Math.tan(r) + 1 / Math.cos(r)) / Math.PI) / 2) * 2 ** z * 256; };
  const invX = (x, z) => (x / (2 ** z * 256)) * 360 - 180, invY = (y, z) => (Math.atan(Math.sinh(Math.PI * (1 - (2 * y) / (2 ** z * 256)))) * 180) / Math.PI;
  function draw() {
    const w = box.clientWidth || 600, hgt = box.clientHeight || 360, cx = X(cLon, zoom), cy = Y(cLat, zoom), n = 2 ** zoom;
    { const dpr = Math.min(window.devicePixelRatio || 1, 2); baseCv.width = w * dpr; baseCv.height = hgt * dpr; const g = baseCv.getContext('2d'); g.setTransform(dpr, 0, 0, dpr, 0, 0);
      g.fillStyle = '#a9d3e8'; g.fillRect(0, 0, w, hgt); g.fillStyle = '#eef0e4'; g.strokeStyle = '#8fa3ad'; g.lineWidth = 0.8;
      const world = n * 256;
      for (const k of [-1, 0, 1]) { const ox = k * world - cx + w / 2; if (ox > w || ox + world < 0) continue;
        for (const ring of LAND) { g.beginPath(); for (let q = 0; q < ring.length; q += 2) { const px = X(ring[q], zoom) + ox, py = Y(Math.max(-85, Math.min(85, ring[q + 1])), zoom) - cy + hgt / 2; if (q) g.lineTo(px, py); else g.moveTo(px, py); } g.closePath(); g.fill(); if (!borders) g.stroke(); }
        if (borders) { // finer outlines with the country borders
          const lonL = invX(-ox, zoom), lonR = invX(w - ox, zoom), latT = invY(cy - hgt / 2, zoom), latB = invY(cy + hgt / 2, zoom), minSize = 3 / ((n * 256) / 360);
          for (const p of borders) {
            if (p.x1 < lonL || p.x0 > lonR || p.y1 < latB || p.y0 > latT || (p.x1 - p.x0 < minSize && p.y1 - p.y0 < minSize)) continue;
            g.beginPath(); for (let q = 0; q < p.xs.length; q++) { const px = X(p.xs[q], zoom) + ox, py = Y(Math.max(-85, Math.min(85, p.ys[q])), zoom) - cy + hgt / 2; if (q) g.lineTo(px, py); else g.moveTo(px, py); } g.closePath(); g.fill(); g.stroke();
          }
        } }
      g.strokeStyle = 'rgba(60,90,110,.18)'; g.lineWidth = 0.6; const stepDeg = zoom <= 3 ? 30 : zoom <= 5 ? 10 : zoom <= 8 ? 2 : 0.5;
      for (let lo = -180; lo <= 180; lo += stepDeg) { const px = X(lo, zoom) - cx + w / 2; if (px >= 0 && px <= w) { g.beginPath(); g.moveTo(px, 0); g.lineTo(px, hgt); g.stroke(); } }
      for (let la = -80; la <= 80; la += stepDeg) { const py = Y(la, zoom) - cy + hgt / 2; if (py >= 0 && py <= hgt) { g.beginPath(); g.moveTo(0, py); g.lineTo(w, py); g.stroke(); } } }
    clear(layer);
    for (let tx = Math.floor((cx - w / 2) / 256); tx <= Math.floor((cx + w / 2) / 256); tx++) for (let ty = Math.floor((cy - hgt / 2) / 256); ty <= Math.floor((cy + hgt / 2) / 256); ty++) {
      if (ty < 0 || ty >= n) continue;
      const img = h('img', { alt: '', draggable: 'false' });
      img.style.left = Math.round(tx * 256 - cx + w / 2) + 'px'; img.style.top = Math.round(ty * 256 - cy + hgt / 2) + 'px';
      layer.append(img);
      loadTile(zoom, ((tx % n) + n) % n, ty).then((url) => { if (url && img.isConnected) img.src = url; else { img.remove(); if (!url) needBorders(); } });
    }
    if (mLat !== null && mLat !== undefined) { pin.hidden = false; pin.style.left = X(mLon, zoom) - cx + w / 2 + 'px'; pin.style.top = Y(mLat, zoom) - cy + hgt / 2 + 'px'; } else pin.hidden = true;
  }
  let drag = null;
  box.addEventListener('pointerdown', (e) => { if (e.target.closest('.map-zoom')) return; drag = { x: e.clientX, y: e.clientY, cx: X(cLon, zoom), cy: Y(cLat, zoom), moved: false }; box.setPointerCapture(e.pointerId); });
  box.addEventListener('pointermove', (e) => { if (!drag) return; const dx = e.clientX - drag.x, dy = e.clientY - drag.y; if (Math.abs(dx) + Math.abs(dy) > 5) drag.moved = true; if (drag.moved) { cLon = invX(drag.cx - dx, zoom); cLat = Math.max(-84, Math.min(84, invY(drag.cy - dy, zoom))); draw(); } });
  box.addEventListener('pointerup', (e) => {
    if (drag && !drag.moved) { const r = box.getBoundingClientRect(), px = X(cLon, zoom) + (e.clientX - r.left - r.width / 2), py = Y(cLat, zoom) + (e.clientY - r.top - r.height / 2); mLat = invY(py, zoom); mLon = ((invX(px, zoom) + 540) % 360) - 180; draw(); onPick(mLat, mLon); }
    drag = null;
  });
  box.addEventListener('wheel', (e) => { e.preventDefault(); zoom = Math.max(2, Math.min(15, zoom + (e.deltaY < 0 ? 1 : -1))); draw(); }, { passive: false });
  box.addEventListener('keydown', (e) => { const s = 40 / 2 ** zoom; if (e.key === 'ArrowLeft') cLon -= s; else if (e.key === 'ArrowRight') cLon += s; else if (e.key === 'ArrowUp') cLat += s; else if (e.key === 'ArrowDown') cLat -= s; else if (e.key === '+') zoom++; else if (e.key === '-') zoom--; else return; e.preventDefault(); draw(); });
  new ResizeObserver(draw).observe(box);
  box.setView = (la, lo, z) => { cLat = la; cLon = lo; mLat = la; mLon = lo; if (z) zoom = z; draw(); };
  return box;
}

export function sitePage(root) {
  const s = store.case.site;
  const latI = h('input', { type: 'number', step: 'any', min: -90, max: 90, value: s.lat ?? '', id: 's_lat', placeholder: 'e.g. 3.50' }), lonI = h('input', { type: 'number', step: 'any', min: -180, max: 180, value: s.lon ?? '', id: 's_lon', placeholder: 'e.g. 5.75' });
  const srcBox = h('div', { class: 'sources' }), dataBox = h('div'), results = h('ul', { class: 'search-results' });
  const map = slippyMap(s.lat, s.lon, (la, lo) => { latI.value = la.toFixed(4); lonI.value = lo.toFixed(4); });
  const paintSources = (live = {}) => fill(srcBox, SOURCES.map((src) => { const st = live[src.id] || store.case.site.status?.[src.id]; return h('div', { class: 'source ' + (st === 'loading' ? 'loading' : st?.ok ? 'ok' : st ? 'fail' : '') }, h('b', null, src.name), h('span', null, src.gives), h('small', null, src.provider + ' · ' + src.host), h('em', { title: st?.atlas ? `Live source not reachable${st.reason ? ' (' + st.reason + ')' : ''}: its values are answered from the world atlas bundled with the app.` : null }, st === 'loading' ? 'fetching…' : st ? (st.ok ? (st.cached ? 'on this device · fetched ' : 'live · ') + ago(st.at) : st.atlas ? st.message + ' — live source not reachable' : 'unavailable — ' + st.message) : 'not fetched')); }));
  let plotCache = new Map();
  const paintData = () => {
    const site = store.case.site, d = site.data || {};
    clear(dataBox);
    if (!site.fetchedAt) return dataBox.append(h('p', { class: 'note' }, 'Pick the field location on the map — offshore or on land — (or search / type coordinates) and press “Fetch live site data”.'));
    const AF = new Set(d.atlasFields || []), notes = d.atlasNotes || {}, vint = d.atlasVintage || {};
    // a tile whose value came from the built-in atlas says so in its label, is tinted, and explains itself on hover
    const K = (label, v, unit, key, extra) => { if (v === null || v === undefined) return null; const a = key && AF.has(key); return { label: a ? label + ' · built-in atlas' : label, value: v, unit, status: a ? 'warn' : undefined, help: a ? `From the built-in world atlas, not live${notes[key] ? ' — ' + notes[key] : ''}.` : extra }; };
    dataBox.append(h('h2', { class: 'sect' }, `${site.name || 'Site'}${site.country ? ', ' + site.country : ''}${AF.has('country') ? ' (country from the built-in atlas)' : ''} — ${fmt(site.lat, 5)}°, ${fmt(site.lon, 5)}°`));
    if (AF.size) {
      const groups = [['sea temperature and salinity', ['sst', 'sstMonthly', 'salinity', 'salinityMonthly']], ['water depth and the seabed grid', ['depth', 'bathy', 'maxDepthNearby', 'elevation', 'seabedSlope']], ['seabed temperature and the temperature profile', ['seabedTemp', 'tempProfile']], ['tides, currents and waves', ['tideRange', 'tide', 'currentSpeed', 'currentMax', 'currentDir', 'waveHeight', 'wavePeriod', 'waveDir']],
        ['air-temperature and wind climate', ['airTempAnnual', 'airTemp', 'windSpeed', 'windAnnual', 'ghiAnnual', 'ghiDaily']], ['country and national economic and energy figures', ['country', 'currency', 'inflation', 'lendingRate', 'gdpPerCapita', 'gridCarbon', 'renewableShare', 'electricityPrice', 'electricityPriceWB', 'renewableElectricity']], ['exchange rate', ['fxPerUSD']], ['oil and gas prices', ['oilPrice', 'oilPriceWTI', 'gasPrice', 'gasPriceEurope', 'gasPriceAsia', 'gasPriceRegional', 'steelPrice']], ['cost indices', ['costIndex', 'costEscalation', 'steelIndex', 'usCpi', 'usInflation5y']], ['tax and interest rates', ['corporateTaxRate', 'bondYield', 'treasuryBillYield', 'policyRate']]];
      const single = [...new Set([...AF].map((k) => ATLAS_LABELS[k]).filter(Boolean))], names = single.length <= 6 ? single : groups.filter(([, ks]) => ks.some((k) => AF.has(k))).map(([nm]) => nm), dates = [['ocean', 'sea climatology'], ['relief', 'relief'], ['deep', 'temperature at depth'], ['coast', 'tides, currents and waves'], ['nations', 'national figures'], ['fx', 'exchange rates'], ['climate', 'climate'], ['prices', 'prices'], ['costs', 'cost indices'], ['markets', 'gas and steel markets'], ['fiscal', 'tax rates'], ['rates', 'interest rates'], ['power', 'electricity tariffs']].filter(([k]) => vint[k]).map(([k, nm]) => `${nm} ${vint[k]}`);
      dataBox.append(h('div', { class: 'warns', role: 'status' }, h('div', { class: 'warn' }, h('b', null, 'Built-in atlas in use.'), `Some values come from the built-in world atlas (either because a live service could not be reached, or because no open live service publishes them): ${names.join(single.length <= 6 ? ', ' : '; ')}. They are long-term or recent-period statistics bundled with the app${dates.length ? ' (' + dates.join('; ') + ')' : ''}, not current conditions; each one is marked “built-in atlas” below. Fetch again when the live services are reachable to replace them.`)));
    }
    const onLand = !(d.depth > 0);
    dataBox.append(
      h('h3', null, onLand ? 'Ground and climate at the site' : 'Seabed and sea at the site'),
      kpiGrid([K(/nearshore/.test(d.depthEstimated || '') ? 'Water depth (nearshore estimate)' : 'Water depth at point', d.depth, 'm', 'depth'), K('Seabed temperature', d.seabedTemp, '°C', 'seabedTemp', 'Annual mean at the seabed: the governing ambient temperature for hydrate and wax in a subsea line.'), K(onLand ? 'Ground slope' : 'Seabed slope', d.seabedSlope, '°', 'seabedSlope'), K('Deepest nearby', d.maxDepthNearby, 'm', 'maxDepthNearby'), K('Land elevation', onLand ? d.elevation : null, 'm', 'elevation'), K('Ground temperature (0.5 m)', onLand ? d.groundTemp : null, '°C', 'groundTemp', 'Soil temperature at about half a metre: the ambient temperature of a buried onshore line.'),
        K('Sea-surface temperature', d.sst, '°C', 'sst'), K(notes.salinity === 'regional estimate' ? 'Salinity (regional estimate)' : 'Salinity (climatology)', d.salinity, 'g/kg', 'salinity'), K('Mean current', d.currentSpeed, 'm/s', 'currentSpeed'), K('Peak current', d.currentMax, 'm/s', 'currentMax'),
        K(AF.has('tideRange') ? 'Tidal range (typical)' : 'Tidal range', d.tideRange, 'm', 'tideRange'), K(AF.has('waveHeight') ? 'Wave height (annual mean)' : 'Wave height', d.waveHeight, 'm', 'waveHeight'), AF.has('waveHeight') ? K('Wave height (95th percentile)', d.waveHeightP95, 'm', 'waveHeight') : null, K('Wave period', d.wavePeriod, 's', 'wavePeriod'),
        K(AF.has('airTemp') ? 'Air temperature (month mean)' : 'Air temperature', d.airTemp, '°C', 'airTemp'), K('Coldest air this fortnight', d.airTempMin, '°C'), K('Warmest air this fortnight', d.airTempMax, '°C'), K(AF.has('windSpeed') ? 'Wind speed (month mean)' : 'Wind speed', d.windSpeed, 'm/s', 'windSpeed'), K('Air temperature (long-term)', d.airTempAnnual, '°C', 'airTempAnnual')].filter(Boolean)),
      h('h3', null, 'Prices and markets'),
      kpiGrid([K('Brent crude' + (d.oilPriceDate ? ` (${d.oilPriceDate})` : ''), d.oilPrice, '$/bbl', 'oilPrice'), K('Brent, 30-day mean', d.oilPriceMean30, '$/bbl'), K('WTI crude' + (d.oilPriceWTIDate ? ` (${d.oilPriceWTIDate})` : ''), d.oilPriceWTI, '$/bbl', 'oilPriceWTI'), K('Henry Hub gas' + (d.gasPriceDate ? ` (${d.gasPriceDate})` : ''), d.gasPrice, '$/MMBtu', 'gasPrice'),
        K('European gas (TTF)' + (d.gasPriceEuropeDate ? ` (${d.gasPriceEuropeDate})` : ''), d.gasPriceEurope, '$/MMBtu', 'gasPriceEurope', d.gasPriceEuropeSource), K('Asian LNG (Japan)' + (d.gasPriceAsiaDate ? ` (${d.gasPriceAsiaDate})` : ''), d.gasPriceAsia, '$/MMBtu', 'gasPriceAsia', d.gasPriceAsiaSource), K(`Regional gas price${d.gasPriceRegionalName ? ' — ' + d.gasPriceRegionalName : ''}`, d.gasPriceRegional, '$/MMBtu', 'gasPriceRegional', 'The gas marker nearest this site; suite 7 uses it for gas sales.'),
        K('Brent volatility (1 year)', d.oilPriceVolatility != null ? 100 * d.oilPriceVolatility : null, '%/y', null, 'Annualised standard deviation of daily log returns over the last year; used by the real-options and Monte Carlo analyses of suite 7.'),
        K(`Oil-field equipment cost index${d.costIndexDate ? ' (' + d.costIndexDate + ')' : ''}`, d.costIndex, d.costIndexBase || '', 'costIndex', `${d.costIndexName || 'Producer-price index'} — ${d.costIndexSource || ''}. Suite 7 escalates its cost basis to today with it.`), K('Equipment cost escalation (5 years)', d.costEscalation, '%/y', 'costEscalation'),
        K(`Steel mill products index${d.steelIndexDate ? ' (' + d.steelIndexDate + ')' : ''}`, d.steelIndex, d.steelIndexBase || '', 'steelIndex', d.steelIndexName), K('Iron ore' + (d.steelPriceDate ? ` (${d.steelPriceDate})` : ''), d.steelPrice, d.steelPriceUnit || '$/t', 'steelPrice', d.steelPriceWhat), K('US inflation (5 years)', d.usInflation5y, '%/y', 'usInflation5y', 'Used to convert between money of the day and real terms for dollar costs.')].filter(Boolean)),
      h('h3', null, 'Fiscal terms and national figures'),
      kpiGrid([K('Tax on upstream profit' + (d.taxRateYear ? ` (${d.taxRateYear})` : ''), d.taxRate, '%', null, [d.taxRateBasis, d.taxRateSource].filter(Boolean).join(' — ') || 'Edit the fiscal terms in suite 7.'), K('Fiscal regime', d.fiscalRegime ? String(d.fiscalRegime).replace(/_/g, ' ') : null, '', null, d.fiscalNote), K('Royalty', d.royaltyRate, '%', null, d.fiscalSource), K('Petroleum / resource tax', d.petroleumTaxRate, '%', null, d.fiscalSource), K('Marginal government take', d.marginalTake, '%', null, d.fiscalSource),
        K('Statutory corporate tax' + (d.corporateTaxYear ? ` (${d.corporateTaxYear})` : ''), d.corporateTaxRate, '%', 'corporateTaxRate', d.corporateTaxSource), K('Carbon price' + (d.carbonPriceYear ? ` (${d.carbonPriceYear})` : ''), d.carbonPrice, '$/tCO₂', null, 'Emissions-weighted national carbon price.'),
        K('Government bond yield' + (d.bondYieldDate ? ` (${d.bondYieldDate})` : ''), d.bondYield, '%/y', 'bondYield', d.bondYieldSource), K('Treasury-bill yield' + (d.treasuryBillYieldDate ? ` (${d.treasuryBillYieldDate})` : ''), d.treasuryBillYield, '%/y', 'treasuryBillYield', d.treasuryBillYieldSource), K('Central-bank policy rate' + (d.policyRateDate ? ` (${d.policyRateDate})` : ''), d.policyRate, '%/y', 'policyRate', d.policyRateSource),
        K('Inflation' + (d.inflationYear ? ` (${d.inflationYear})` : ''), d.inflation, '%/y', 'inflation'), K('Lending rate' + (d.lendingRateYear ? ` (${d.lendingRateYear})` : ''), d.lendingRate, '%/y', 'lendingRate'), K(`${d.currency || ''} per USD` + (AF.has('fxPerUSD') && d.fxDate ? ` (${d.fxDate})` : ''), d.fxPerUSD, '', 'fxPerUSD'),
        K('Electricity for industry' + (d.electricityPriceDate ? ` (${d.electricityPriceDate})` : ''), d.electricityPrice, '$/kWh', 'electricityPrice', [d.electricityPriceBasis, d.electricityPriceSource].filter(Boolean).join(' — ')),
        K((d.gridCarbonLive || AF.has('gridCarbon')) && !notes.gridCarbon ? `Grid carbon (${d.gridCarbonYear || 'latest'})` : 'Grid carbon (indicative)', d.gridCarbon, 'kgCO₂/kWh', 'gridCarbon'), K('Renewable electricity' + (AF.has('renewableShare') && d.renewableShareYear ? ` (${d.renewableShareYear})` : ''), d.renewableShare, '%', 'renewableShare')].filter(Boolean)),
      d.fiscalUrl ? h('p', { class: 'note' }, 'Fiscal terms: ', d.fiscalSource || '', ' ', h('a', { href: d.fiscalUrl, target: '_blank', rel: 'noopener noreferrer' }, 'source'), d.fiscalNote ? ' — ' + d.fiscalNote : '') : null);
    const plots = [];
    if (d.bathy) {
      const flat = d.bathy.elev.flat(), zlo = Math.min(...flat), zhi = Math.max(...flat), allLand = zlo >= 0, allSea = zhi <= 0;
      const coarse = AF.has('bathy'), span = Math.round(Math.abs(d.bathy.lat[d.bathy.lat.length - 1] - d.bathy.lat[0]) * 55.6);
      plots.push({ type: 'field', title: (allLand ? 'Terrain around the site (elevation above sea level, m)' : allSea ? 'Seabed around the site (m, negative = below sea level)' : 'Seabed and terrain around the site (m; blue = sea, green to brown = land)') + (coarse ? ' — built-in atlas' : ''), xlabel: 'Longitude (°)', ylabel: 'Latitude (°)', zlabel: 'Elevation (m)', zunit: 'm', x: d.bathy.lon, y: d.bathy.lat, z: d.bathy.elev,
        cmap: allLand ? 'land' : allSea ? 'sea' : 'topo', zmid: allLand || allSea ? undefined : 0, contours: allLand ? 8 : [...[-1000, -500, -200, -100, -50, -20, -10].filter((q) => q > zlo && q < zhi).map((level) => ({ level, color: 'rgba(8,48,107,.45)', width: 0.7 })), ...(allSea ? [] : [{ level: 0, color: '#ffffff', width: 1.6 }])], equal: true, shade: 'geo', exaggeration: zhi - zlo < 60 ? 10 : zhi - zlo < 300 ? 4 : 1.6, markers: [{ x: site.lon, y: site.lat, label: 'site' }],
        onPick: (lo, la) => { latI.value = la.toFixed(4); lonI.value = lo.toFixed(4); map.setView(la, lo); toast(`Site moved to ${la.toFixed(4)}°, ${lo.toFixed(4)}° — fetching live data for the new point…`); doFetch(); },
        note: (coarse ? `Built-in atlas relief${vint.relief ? ' (' + vint.relief + ')' : ''}, not the live survey grid: about 5 km resolution at the coast and 25 km elsewhere, so canyons, scarps and channels are not resolved. ` : '') + 'Click anywhere on this panel to move the site to that point. ' + ( allLand ? `This point is inland: the ground is ${fmt(zlo, 3)}–${fmt(zhi, 3)} m above sea level and there is no sea within about ${coarse ? span : 13} km, so the marine figures above come from the nearest sea cell (or are unavailable). For a subsea study, click a point offshore.` : allSea ? 'Open water: the whole window is below sea level. Thin lines are depth contours at 10, 20, 50, 100, 200, 500 and 1000 m.' : 'The white line is the shoreline (0 m); thin blue lines are depth contours at 10, 20, 50, 100, 200, 500 and 1000 m.') });
    }
    if (d.tempProfile && d.depth > 0) plots.push({ type: 'line', title: 'Sea temperature with depth at the site — built-in atlas', xlabel: 'Temperature (°C)', ylabel: 'Elevation (m, sea level = 0)', series: [{ name: 'Annual-mean temperature', x: d.tempProfile.T, y: d.tempProfile.depth.map((z) => -z), mode: 'both' }], hlines: [{ y: -d.depth, label: `seabed ${fmt(d.depth, 4)} m`, color: '#b45309' }], note: `${notes.seabedTemp || 'World Ocean Atlas 2023 annual climatology'}. This profile is the ambient temperature seen by a riser from the seabed to the surface; suites 2, 3 and 5 take the seabed and surface values from it.` });
    for (const [key, name, unit] of [['oilPrice', 'Brent crude oil', '$/bbl'], ['gasPrice', 'Henry Hub natural gas', '$/MMBtu']]) { const sr = d[key + 'Series']; if (sr?.v?.length > 20) plots.push({ type: 'line', title: `${name} — daily spot price, last year`, xlabel: 'Trading days before the latest price', ylabel: unit, series: [{ name, x: sr.v.map((_, i) => i - sr.v.length + 1), y: sr.v }], hlines: [{ y: d[key + 'Mean30'], label: '30-day mean', color: '#64748b' }], note: `Latest ${fmt(d[key], 4)} ${unit} on ${d[key + 'Date']}; one-year range ${fmt(d[key + 'Min1y'], 4)}–${fmt(d[key + 'Max1y'], 4)}. US Energy Information Administration spot series.` }); }
    if (d.costIndexSeries?.v?.length > 12) plots.push({ type: 'line', title: `${d.costIndexName || 'Oil-field equipment cost index'} (${d.costIndexBase || 'index'})`, xlabel: 'Months before the latest value', ylabel: 'Index', series: [{ name: 'Cost index', x: d.costIndexSeries.v.map((_, i) => i - d.costIndexSeries.v.length + 1), y: d.costIndexSeries.v }], note: `${d.costIndexSource || ''}. Latest ${fmt(d.costIndex, 4)} in ${d.costIndexDate}; average escalation over the last five years ${fmt(d.costEscalation, 3)} %/y. Suite 7 uses this series to bring every cost from its basis year to today.` });
    if (d.tide) {
      const now = d.tide.nowHour ?? 72, m0 = d.seaLevelMean ?? d.tide.eta.reduce((a, b) => a + b, 0) / d.tide.eta.length;
      const synth = AF.has('tide') || d.tide.synthetic;
      plots.push({ type: 'line', title: synth ? 'Tide about mean sea level — built-in atlas (harmonic prediction)' : 'Tide about the local mean sea level — past 3 days and forecast', xlabel: 'Hours from now (negative = past)', ylabel: 'm about local mean level', series: [{ name: 'Sea level', x: d.tide.t.map((t) => t - now), y: d.tide.eta.map((e) => e - m0) }], vlines: [{ x: 0, label: 'now', color: '#f97316' }], hlines: [{ y: 0, label: 'local mean', color: '#64748b' }],
        note: synth ? `Synthetic series from the built-in atlas, not the live ocean model: a harmonic prediction from five tidal constituents (M2, S2, N2, K1, O1) fitted to modelled sea level${vint.seaPeriod ? ' of ' + vint.seaPeriod : ''} at the nearest coastal model point${notes.tideRange && /km away/.test(notes.tideRange) ? ' (' + notes.tideRange.replace(/^.*model point /, '') + ')' : ''}. Weather-driven surge and the smaller constituents are not included. Typical range ${fmt(d.tideRange, 3)} m${d.tideSpring ? `, spring range ${fmt(d.tideSpring, 3)} m` : ''}.` : `Tidal range over this window: ${fmt(d.tideRange, 3)} m. The local mean level is ${fmt(m0, 2)} m relative to the ocean model’s global mean-sea-level datum; that constant offset (regional sea-surface height, not tide) has been removed from the curve.` });
    }
    if (d.currents) { const now = d.currents.nowHour ?? 72; plots.push({ type: 'line', title: 'Ocean-current speed — past 3 days and forecast', xlabel: 'Hours from now (negative = past)', ylabel: 'm/s', zeroY: true, series: [{ name: 'Current speed', x: d.currents.t.map((t) => t - now), y: d.currents.speed, mode: 'step' }], vlines: [{ x: 0, label: 'now', color: '#f97316' }], hlines: [{ y: d.currentSpeed, label: 'mean', color: '#64748b' }],
      note: `The ocean model publishes currents in steps of 0.1 km/h (about 0.03 m/s), which is why the trace is stepped at low speeds.${d.currents.gaps ? ` ${d.currents.gaps} missing hours in the source were bridged by interpolation.` : ''} Mean ${fmt(d.currentSpeed, 2)} m/s, peak ${fmt(d.currentMax, 2)} m/s.` }); }
    if (d.salinityMonthly) plots.push({ type: 'line', title: 'Monthly climatology near the site' + (AF.has('salinityMonthly') || AF.has('sstMonthly') ? ' — built-in atlas' : ''), note: AF.has('salinityMonthly') || AF.has('sstMonthly') ? `Built-in atlas${vint.ocean ? ' (' + vint.ocean + ')' : ''}, not the live climatology service: salinity from the SeaDataCloud climatology on 1° cells (about 100 km), temperature from the NOAA CoralTemp satellite record on 0.5° cells near land; both as a smooth seasonal cycle (annual mean and two harmonics)${notes.salinity ? '; salinity: ' + notes.salinity : ''}.` : undefined, xlabel: 'Month', ylabel: 'Salinity (g/kg) · temperature (°C)', legendBelow: true, series: [{ name: 'Salinity', x: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], y: d.salinityMonthly, mode: 'both' }, { name: 'Temperature', x: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], y: d.sstMonthly, mode: 'both' }] });
    if (plots.length) { // reuse chart cards whose data are unchanged, so progressive updates do not redraw everything
      const sig = (p) => p.title + '|' + (p.type === 'field' ? p.z.length + ':' + p.z[0][0] + ':' + p.z[p.z.length - 1][p.z[0].length - 1] + ':' + (p.markers?.[0]?.x ?? '') + ':' + (p.markers?.[0]?.y ?? '') : p.series.map((q) => q.y.length + ':' + q.y[0] + ':' + q.y[q.y.length - 1] + ':' + q.x[0]).join(','));
      const next = new Map(); dataBox.append(h('div', { class: 'plots' }, plots.map((p) => { const k = sig(p), card = plotCache.get(k) || pc(p); next.set(k, card); return card; }))); plotCache = next;
    }
    dataBox.append(card(h('h2', null, 'Where these values go'), h('p', { class: 'note' }, 'Open any suite: matching inputs appear in its “linked data” bar and, with auto-link on, are applied when you run. Water depth, seabed and surface temperature and currents go to the network, flow, solids and operations suites; waves and currents to integrity; prices, tax, inflation, exchange rate, electricity and carbon figures to economics.'),
      h('div', { class: 'row-tools' }, h('a', { class: 'btn primary', href: '#/chain' }, 'Run every suite with these site data →'),
        btn('Download site data (JSON)', () => download(JSON.stringify(site, null, 1), 'site-data.json', 'application/json')))));
  };
  const doFetch = async (opt = {}) => {
    const la = +latI.value, lo = +lonI.value;
    if (!Number.isFinite(la) || !Number.isFinite(lo) || latI.value === '' || lonI.value === '' || Math.abs(la) > 90 || Math.abs(lo) > 180) return toast('Enter a valid latitude (−90…90) and longitude (−180…180), or click the map.', 'warn');
    if (!navigator.onLine) opt = { ...opt, atlasOnly: true }; // no connection: answer from the built-in atlas instead of refusing
    const t0 = performance.now(), live = {}, old = store.case.site, same = old.lat !== null && Math.abs(old.lat - la) < 0.02 && Math.abs(old.lon - lo) < 0.02;
    fetchBtn.disabled = true; fetchBtn.textContent = 'Fetching…';
    let raf = 0, latest = null, first = true;
    const commit = (site) => { // keep earlier answers of sources that are down right now, then repaint
      if (same) site.data = mergeSiteData(old.data, site.data); // a stored live value is not replaced by an atlas value
      store.setSite(site);
    };
    const paintSoon = (site) => { latest = site; if (raf) return; raf = requestAnimationFrame(() => { raf = 0; commit(latest); if (first) { map.setView(latest.lat, latest.lon); first = false; } paintSources(live); paintData(); }); };
    try {
      const site = await fetchSite(la, lo, (id, st, msg) => { live[id] = st === 'loading' ? 'loading' : { ok: st === 'ok', atlas: st === 'atlas', reason: st === 'atlas' ? live[id]?.message : undefined, message: msg, at: new Date().toISOString() }; paintSources(live); }, paintSoon, opt);
      if (raf) { cancelAnimationFrame(raf); raf = 0; }
      commit(site); map.setView(site.lat, site.lon);
      const ok = Object.values(site.status).filter((x) => x.ok).length, stored = Object.values(site.status).filter((x) => x.cached).length, atl = Object.values(site.status).filter((x) => x.atlas).length;
      toast(ok ? `Site data ready in ${((performance.now() - t0) / 1000).toFixed(1)} s: ${ok} of ${SOURCES.length} sources live${stored ? ` (${stored} from the copy kept on this device)` : ''}${atl ? `; ${atl} answered from the built-in atlas` : ''}.` : 'Live services could not be reached from here — showing the built-in world atlas (sea climatology, depth, tides, currents, waves, climate and national figures; every value is marked).', ok && !atl ? 'ok' : 'warn', ok && !atl ? 4200 : 9000);
    } catch (e) { toast('Could not fetch site data: ' + e.message, 'bad'); }
    fetchBtn.disabled = false; fetchBtn.textContent = 'Fetch live site data'; paintSources(); paintData();
  };
  const fetchBtn = btn('Fetch live site data', () => doFetch(), 'primary', 'Each source fills in as soon as it answers; answers already on this device appear instantly');
  const freshBtn = btn('Force refresh', () => doFetch({ fresh: true }), 'ghost', 'Ignore the copies kept on this device and ask every source again');
  const q = h('input', { type: 'search', placeholder: 'Search a city, port, terminal or basin…', 'aria-label': 'Search place', maxlength: 80 });
  const doSearch = async () => {
    if (!q.value.trim()) return;
    try { const r = await searchPlace(q.value.trim()); fill(results, r.length ? r.map((p) => h('li', null, h('button', { type: 'button', class: 'linklike', onclick: () => { latI.value = p.lat.toFixed(4); lonI.value = p.lon.toFixed(4); map.setView(p.lat, p.lon, 9); clear(results); } }, `${p.name}${p.admin ? ', ' + p.admin : ''}, ${p.country}`))) : h('li', { class: 'note' }, 'No match — try another spelling.')); }
    catch { toast('Search is unavailable (offline?). Type coordinates instead.', 'warn'); }
  };
  q.addEventListener('keydown', (e) => { if (e.key === 'Enter') doSearch(); });
  paintSources(); paintData();
  if (s.fetchedAt && s.data?.bathy && s.data.bathy.lat.length < 40 && navigator.onLine) setTimeout(() => doFetch(), 300); // stored by an earlier build with a coarse relief grid: refresh once
  fill(root, 
    h('header', { class: 'page-head' }, h('h1', null, 'Global site data'), h('p', null, 'Choose any location on Earth — land, shelf or deep water. Your browser pulls the seabed or terrain relief, water depth, seabed temperature and the temperature profile of the water column, sea state, tides and currents, weather and ground temperature, today’s oil and gas prices, exchange rate, and national economic, energy and carbon figures directly from open global data services, and offers them to every suite.')),
    card(
      h('div', { class: 'site-bar' },
        h('div', { class: 'site-search' }, q, btn('Search', doSearch)),
        h('label', { class: 'site-coord', for: 's_lat' }, h('span', null, 'Latitude °N'), latI), h('label', { class: 'site-coord', for: 's_lon' }, h('span', null, 'Longitude °E'), lonI),
        fetchBtn, freshBtn, btn('Use my location', () => navigator.geolocation?.getCurrentPosition((p) => { latI.value = p.coords.latitude.toFixed(4); lonI.value = p.coords.longitude.toFixed(4); map.setView(p.coords.latitude, p.coords.longitude, 9); }, () => toast('Location permission was not granted.', 'warn')), 'ghost')),
      results, map,
      h('p', { class: 'note' }, 'Click or tap the map to place the site. For a subsea development choose the field location (the deepest point of the route governs the seabed temperature); for an onshore line the marine sources report the nearest sea cell and the ground temperature is used instead.'),
      h('h3', null, 'Live sources'), srcBox),
    dataBox);
}

// --------------------------------------------------------------------------------------- data portal
export function portalPage(root) {
  const out = h('div');
  const TABULAR = ['csv', 'tsv', 'txt', 'dat', 'xlsx', 'xlsm'];
  // Show an imported or generated geometry and offer every suite that can use it.
  const showGeometry = async (g) => {
    clear(out);
    const d = derive(g), routes = [];
    for (const m of SUITES) { try { const su = await loadSuite(m.id), items = geometryLinks(su, g, d), accepts = (SUITE_GEOMETRY[m.id]?.accepts || []).includes(g.kind); if (items.length || accepts) routes.push({ su, items }); } catch { /* suite not available */ } }
    routes.sort((a, b) => b.items.length - a.items.length);
    const send = (su, items) => { attachGeometry(su.id, g); for (const it of items) setInputValue(su, it.key, it.value); store.pref('tab.' + su.id, 'geometry'); toast(`${g.name} sent to ${su.short}${items.length ? ` — ${items.length} input${items.length > 1 ? 's' : ''} set` : ''}.`, 'ok'); location.hash = '#/suite/' + su.id; };
    out.append(card(h('h2', null, `Geometry: ${g.name}`), geometryCard(g),
      h('h3', null, 'Send to a suite'),
      routes.length ? h('div', { class: 'tbl-scroll' }, h('table', { class: 'tbl' }, h('thead', null, h('tr', null, ['Suite', 'What it will take from this geometry', ''].map((x) => h('th', null, x)))),
        h('tbody', null, routes.map(({ su, items }) => h('tr', null, h('td', { class: 'lead' }, `${su.icon} ${su.num}. ${su.title}`), h('td', { class: 'lead wrap' }, items.length ? items.map((it) => `${allFields(su).find((f) => f.key === it.key)?.label || it.key} ← ${it.from}`).join('; ') : 'Attached for reference: preview and measured dimensions on its Geometry tab'), h('td', null, btn(items.length ? 'Send & apply' : 'Attach', () => send(su, items), items.length ? 'mini primary' : 'mini')))))))
        : h('p', { class: 'note' }, 'No suite uses this class of geometry directly. Its measured dimensions are shown above.')));
  };
  const handle = async (files) => {
    for (const f of files) checkFile(f);
    const file = files[0], ext = extOf(file.name);
    clear(out);
    if (files.length === 1 && ext === 'json' && /hydraslug/i.test(file.name)) { store.importJSON(await file.text()); toast('Case imported.', 'ok'); return out.append(card(h('h2', null, 'Case imported'), h('p', null, `“${store.case.name}” is now the active case.`), h('a', { class: 'btn primary', href: '#/case' }, 'Open the case'))); }
    if (files.length > 1 || !TABULAR.includes(ext)) return showGeometry(await readFiles(files));
    const t = await readTable(file);
    if (!t.records.length) throw new Error('The file contains no data rows.');
    const targets = [];
    for (const m of SUITES) { try { const su = await loadSuite(m.id); allFields(su).filter((f) => f.type === 'table').forEach((f) => targets.push({ su, f })); } catch { /* skip */ } }
    const sel = h('select', { 'aria-label': 'Destination' }, targets.map((x, i) => h('option', { value: i }, `${x.su.num}. ${x.su.short} → ${x.f.label}`)));
    const mapBox = h('div');
    const norm = (s) => String(s).toLowerCase().replace(/[^a-z0-9]/g, '');
    const paintMap = () => {
      const { f } = targets[+sel.value] || {};
      if (!f) return;
      fill(mapBox, h('table', { class: 'tbl' }, h('thead', null, h('tr', null, h('th', null, 'Destination column'), h('th', null, 'Take from file column'))),
        h('tbody', null, f.columns.map((c, i) => h('tr', null, h('td', { class: 'lead' }, c.label + (c.unit ? ` (${c.unit})` : '')), h('td', null, h('select', { dataset: { key: c.key } }, h('option', { value: '' }, '— leave empty —'), t.headers.map((hd, j) => h('option', { value: hd, selected: norm(hd) === norm(c.key) || norm(hd) === norm(c.label) || (j === i && !t.headers.some((x) => norm(x) === norm(c.key))) }, hd)))))))));
    };
    sel.addEventListener('change', paintMap); paintMap();
    out.append(card(h('h2', null, `Table: ${file.name}`), h('p', { class: 'note' }, `${t.records.length} rows × ${t.headers.length} columns detected.`),
      dataTable({ title: 'Preview (first 12 rows)', columns: t.headers, rows: t.records.slice(0, 12).map((r) => t.headers.map((hd) => r[hd])) }),
      h('div', { class: 'row-tools' }, btn('This file is geometry (coordinates, soundings, point cloud, network table) — read it as geometry', async () => { try { await showGeometry(await readFiles(files)); } catch (e) { toast(e.message, 'bad', 12000); } })),
      targets.length ? h('div', null, h('h3', null, 'Send to a suite as data'), h('div', { class: 'row-tools' }, sel), mapBox,
        h('div', { class: 'row-tools' }, btn('Load into suite', () => { const { su, f } = targets[+sel.value]; const m = Object.fromEntries([...mapBox.querySelectorAll('select')].map((s) => [s.dataset.key, s.value])); store.setInput(su.id, f.key, t.records.map((r) => Object.fromEntries(f.columns.map((c) => [c.key, m[c.key] ? r[m[c.key]] : null])))); toast(`${t.records.length} rows loaded into ${su.short}.`, 'ok'); location.hash = '#/suite/' + su.id; }, 'primary'))) : null,
      h('p', { class: 'note' }, 'A fluid composition (component names with mol %) is imported on the Case page with the “Import analysis” button of the composition editor.')));
  };
  const input = h('input', { type: 'file', hidden: true, multiple: true, accept: ACCEPT + ',.xlsx,.xlsm,.tsv' });
  const safe = async (files) => { try { await handle([...files]); } catch (e) { clear(out); const m = formatOf(files[0]?.name || ''); out.append(card(h('h2', null, `Could not read ${files[0]?.name || 'the file'} directly`), h('p', { class: m?.support === 'convert' ? 'summary' : 'bad' }, e.message || 'Unreadable file.'), m ? h('p', { class: 'note' }, `Recognised as ${m.name}.`) : h('p', { class: 'note' }, 'The file type was not recognised. See the list of supported formats below.'))); toast(m?.support === 'convert' ? 'This format needs one conversion step — see the instruction.' : 'Could not read the file.', m?.support === 'convert' ? 'warn' : 'bad', 8000); } };
  input.addEventListener('change', () => { if (input.files.length) safe([...input.files]); input.value = ''; });
  const drop = h('div', { class: 'drop', tabindex: '0', role: 'button', 'aria-label': 'Choose or drop files', onclick: () => input.click(), onkeydown: (e) => { if (e.key === 'Enter' || e.key === ' ') input.click(); } }, h('div', { class: 'drop-ico', 'aria-hidden': 'true' }, '⬆'), h('b', null, 'Drop files here or click to choose'), h('span', null, 'Files are read on this device only — nothing is uploaded. Select companion files together (e.g. .shp + .dbf, .mhd + .raw).'), input);
  drop.addEventListener('dragover', (e) => { e.preventDefault(); drop.classList.add('over'); });
  drop.addEventListener('dragleave', () => drop.classList.remove('over'));
  drop.addEventListener('drop', (e) => { e.preventDefault(); drop.classList.remove('over'); if (e.dataTransfer.files.length) safe([...e.dataTransfer.files]); });
  fill(root, 
    h('header', { class: 'page-head' }, h('h1', null, 'Data portal'), h('p', null, 'One place to bring case data in: pipeline routes and elevation profiles, well deviation surveys, network tables, CAD and meshes, seabed terrain, inspection and deposit maps, laboratory PVT tables, historian logs, cost tables and complete case files. Each import is recognised, previewed, measured and routed to the suites that can use it.')),
    card(drop, h('div', { class: 'formats' },
      h('div', null, h('b', null, 'Pipeline, riser and well geometry'), h('span', null, 'Chainage–elevation, XYZ centrelines, easting/northing/elevation, latitude/longitude/depth, MD–inclination–azimuth and MD–TVD surveys as CSV, TXT, XLSX, JSON')),
      h('div', null, h('b', null, 'Networks and topology'), h('span', null, 'Node–edge tables, JSON, YAML, XML, GraphML, PCF piping, IFC, AutomationML')),
      h('div', null, h('b', null, 'GIS, bathymetry, terrain'), h('span', null, 'Shapefile, GeoJSON, KML/KMZ, GPX, GML, LandXML, GeoPackage, GeoTIFF, DEM / ASCII grid, NetCDF (classic and 4), XYZ soundings, LAS/LAZ, E57, PTS/PTX, SEG-Y positions')),
      h('div', null, h('b', null, 'CAD and surfaces'), h('span', null, 'STEP (AP203/214/242), IGES, ACIS SAT, Parasolid X_T, DWG, DXF, DGN, STL, OBJ, PLY, OFF, 3MF, glTF/GLB, VRML/X3D, SVG')),
      h('div', null, h('b', null, 'CFD and FEA meshes'), h('span', null, 'CGNS, Exodus II, MED, VTK/VTU/VTP, Gmsh, Fluent, OpenFOAM, UNV, Nastran BDF, Abaqus INP, ANSYS CDB, LS-DYNA K, SU2, Tecplot')),
      h('div', null, h('b', null, 'Inspection and deposits'), h('span', null, 'Wall-thickness, corrosion and deposit maps (x, θ, value), defect lists, point clouds, voxel and image data')),
      h('div', null, h('b', null, 'Tables'), h('span', null, 'CSV, TSV, TXT, JSON, Excel .xlsx, MATLAB .mat, HDF5, NumPy — PVT reports, logs, curves, schedules, costs')),
      h('div', null, h('b', null, 'Cases'), h('span', null, `${APP.caseExt} — every input of all ${N} suites in one file`)))),
    out,
    card(h('h2', null, 'Generate geometry without a file'), h('p', { class: 'note' }, 'Parametric geometry for early studies before survey data exist: catenary and lazy-wave risers, undulating flowlines, build-and-hold well trajectories and jumpers, plus porous and packed structures for deposit studies.'), generatorPanel((g) => showGeometry(g).catch((e) => toast(e.message, 'bad')))),
    card(h('h2', null, 'Supported geometry formats'), h('p', { class: 'note' }, 'Open, documented formats are read directly on this device. That now includes formats usually regarded as closed: AutoCAD DWG (R13 to 2018, except 2007), MicroStation DGN (V7), ACIS SAT and Parasolid X_T text files — with pipe runs, bends and diameters recognised from their cylinders and tori — and the HDF5 family (CGNS, Exodus II, MED, NetCDF-4, MATLAB v7.3), GeoPackage, E57 and LAZ. The few that remain (JT, the binary SAB and X_B files, DWG 2007, DGN V8) are recognised and the portal says exactly which neutral format to export.'), formatCatalogue()));
}

// ------------------------------------------------------------------------------------ integrated run
// Quantities that tie the suites together in both directions. The chain is repeated until none of them moves by more
// than the tolerance, which is what turns seven calculators into one coupled model.
const COUPLED = [
  ['flow', 'pIn', 'Inlet pressure', 'bara'], ['flow', 'tOut', 'Arrival temperature', '°C'], ['flow', 'liquidInventory', 'Liquid inventory', 'm³'],
  ['solids', 'effectiveId', 'Smallest effective bore', 'm'], ['solids', 'inhibitorRequired', 'Inhibitor required', 'wt %'],
  ['ops', 'cooldownTime', 'Cooldown time', 'h'], ['ops', 'inhibitorRate', 'Inhibitor rate', 'm³/d'], ['integ', 'pof', 'Annual failure probability', '1/y'], ['integ', 'remainingLife', 'Remaining life', 'y'],
  ['econ', 'npv', 'Net present value', '$'], ['econ', 'eal', 'Expected annual loss', '$'],
];
export function chainPage(root) {
  const rows = h('tbody'), bar = h('div', { class: 'progress', hidden: true }, h('i')), summary = h('div'), stat = h('span', { class: 'run-status', role: 'status', 'aria-live': 'polite' });
  const sel = Object.fromEntries(CHAIN.map((id) => [id, h('input', { type: 'checkbox', checked: store.pref('chain.' + id) !== false, 'aria-label': 'Include ' + byId(id).title, onchange: (e) => store.pref('chain.' + id, e.target.checked) })]));
  const passes = h('select', { 'aria-label': 'Coupling passes' }, [[1, 'one forward pass (no feedback)'], [2, 'up to 2 passes'], [3, 'up to 3 passes'], [5, 'up to 5 passes']].map(([v, t]) => h('option', { value: v, selected: v === (store.pref('chain.passes') || 5) }, t)));
  passes.addEventListener('change', () => store.pref('chain.passes', +passes.value));
  const cell = {};
  CHAIN.forEach((id, i) => { const s = byId(id); cell[id] = { st: h('td', null, store.case.outputs[id] ? badge('solved ' + ago(store.case.outputs[id]._at), 'ok') : badge('waiting', '')), kp: h('td', { class: 'kp' }, kpText(store.case.outputs[id])), ln: h('td', null, '') }; rows.append(h('tr', null, h('td', null, sel[id]), h('td', { class: 'num' }, i + 1), h('td', { class: 'lead' }, h('a', { href: '#/suite/' + id }, `${s.icon} ${s.num}. ${s.title}`)), cell[id].ln, cell[id].st, cell[id].kp)); });
  function kpText(o) { return o?._kpis?.length ? o._kpis.slice(0, 4).map((k) => `${k.label} ${typeof k.value === 'number' ? fmt(k.value, 3) : k.value} ${k.unit}`).join(' · ') : '—'; }
  const snapshot = () => COUPLED.map(([id, key]) => { const v = store.case.outputs[id]?.[key]; return typeof v === 'number' && Number.isFinite(v) ? v : null; });
  let stop = false;
  const stopBtn = btn('Stop after this suite', () => { stop = true; stopBtn.disabled = true; }, 'ghost danger'); stopBtn.hidden = true;
  const go = btn('▶  Run the selected suites', async () => {
    go.disabled = true; bar.hidden = false; clear(summary); stop = false; stopBtn.hidden = false; stopBtn.disabled = false;
    const list = CHAIN.filter((id) => sel[id].checked), maxPass = Math.max(1, +passes.value || 1), tol = 0.01, history = [];
    let okN = 0, viol = 0, pass = 0, converged = false, failed = 0;
    for (pass = 1; pass <= maxPass && !stop; pass++) {
      // from the second pass on only the suites inside the feedback loops need to run again
      const todo = pass === 1 ? list : list.filter((id) => id !== 'pvt');
      viol = 0;
      for (let i = 0; i < todo.length && !stop; i++) {
        const id = todo[i]; bar.firstChild.style.width = Math.round((100 * ((pass - 1) * list.length + i)) / (maxPass * list.length)) + '%';
        stat.textContent = `Pass ${pass} of up to ${maxPass}: ${byId(id).title}…`;
        clear(cell[id].st).append(badge('running…', 'warn'));
        try {
          const su = await loadSuite(id), items = linkItems(su), n = applyLinks(su, items);
          cell[id].ln.textContent = n ? `${n} linked` : '–'; cell[id].ln.title = items.map((x) => `${x.key} ← ${x.from}`).join('\n');
          await new Promise((r) => setTimeout(r, 15));
          const res = await runSuite(su), bad = res.warnings.filter((w) => w.level === 'bad').length;
          viol += bad; if (pass === 1) okN++;
          clear(cell[id].st).append(badge(bad ? `${bad} limit issue${bad > 1 ? 's' : ''}` : 'solved', bad ? 'warn' : 'ok'), h('small', null, ` ${Math.round(res._ms)} ms${pass > 1 ? ' · pass ' + pass : ''}`));
          cell[id].kp.textContent = kpText(store.case.outputs[id]);
        } catch (e) { console.error(e); failed++; clear(cell[id].st).append(badge('failed', 'bad')); cell[id].kp.textContent = String(e.message || e).slice(0, 200); }
      }
      history.push(snapshot());
      if (history.length > 1) {
        const a = history[history.length - 2], b = history[history.length - 1];
        const change = Math.max(0, ...b.map((v, k) => (v === null || a[k] === null ? 0 : Math.abs(v - a[k]) / Math.max(Math.abs(v), Math.abs(a[k]), 1e-9))));
        if (change < tol) { converged = true; break; }
      }
      if (maxPass === 1 || failed) break;
    }
    const done = Math.min(pass, maxPass);
    bar.firstChild.style.width = '100%'; setTimeout(() => (bar.hidden = true), 600); go.disabled = false; stopBtn.hidden = true;
    stat.textContent = stop ? 'Stopped.' : `Finished after ${history.length} pass${history.length > 1 ? 'es' : ''}.`;
    const o = store.case.outputs, K = (cond, label, value, unit, status) => (cond !== undefined && cond !== null && value !== undefined && value !== null && !(typeof value === 'number' && !Number.isFinite(value)) ? { label, value, unit, status } : null);
    summary.append(h('h2', { class: 'sect' }, 'Case summary'), kpiGrid([
      K(o.pvt?.gor, 'Gas–oil ratio', o.pvt?.gor, 'Sm³/Sm³'), K(o.pvt?.psat, 'Saturation pressure', o.pvt?.psat, 'bara'),
      K(o.net?.length, 'Line length', (o.net?.length || 0) / 1000, 'km'), K(o.net?.uValue, 'U-value', o.net?.uValue, 'W/m²K'),
      K(o.flow?.pIn, 'Inlet pressure', o.flow?.pIn, 'bara'), K(o.flow?.tOut, 'Arrival temperature', o.flow?.tOut, '°C'), K(o.flow?.slug, 'Slugging', o.flow?.slug?.type, '', o.flow?.slug?.type === 'severe' ? 'bad' : o.flow?.slug?.type === 'none' ? 'ok' : 'warn'), K(o.flow?.slug?.surge, 'Liquid surge', o.flow?.slug?.surge, 'm³'),
      K(o.solids?.hydrateRisk, 'Hydrate risk index', o.solids?.hydrateRisk, '0–1', o.solids?.hydrateRisk > 0.5 ? 'bad' : o.solids?.hydrateRisk > 0.15 ? 'warn' : 'ok'), K(o.solids?.inhibitorRequired, 'Inhibitor required', o.solids?.inhibitorRequired, 'wt %'), K(o.solids?.piggingInterval, 'Pigging interval', o.solids?.piggingInterval, 'd'),
      K(o.ops?.cooldownTime, 'Cooldown time', o.ops?.cooldownTime, 'h'), K(o.ops?.maxShutdown, 'Longest safe shutdown', o.ops?.maxShutdown, 'h'), K(o.ops?.uptime, 'Uptime', 100 * (o.ops?.uptime || 0), '%'),
      K(o.integ?.remainingLife, 'Remaining life', o.integ?.remainingLife, 'y'), K(o.integ?.riskLevel, 'Integrity risk', o.integ?.riskLevel, '', /high/.test(o.integ?.riskLevel || '') ? 'bad' : 'ok'),
      K(o.econ?.npv, 'Net present value', (o.econ?.npv || 0) / 1e6, 'M$', o.econ?.npv < 0 ? 'bad' : 'ok'), K(o.econ?.irr, 'Rate of return', o.econ?.irr, '%/y'), K(o.econ?.utc, 'Unit technical cost', o.econ?.utc, '$/boe'), K(o.econ?.bestOption, 'Best strategy', o.econ?.bestOption, ''),
    ].filter(Boolean)), h('p', { class: 'summary' }, `${okN} of ${list.length} suites solved` + (failed ? `; ${failed} run${failed > 1 ? 's' : ''} failed (see the row for the reason)` : '') + (viol ? `, with ${viol} limit issue${viol > 1 ? 's' : ''} to review in the flagged suites.` : ', with no limit violations.') + (history.length > 1 ? (converged ? ` The couplings between the suites settled to within ${100 * tol} % after ${history.length} passes.` : ` The couplings were still moving after ${history.length} passes: allow more passes, or look at the table below to see which quantity keeps changing.`) : maxPass > 1 ? '' : ' One forward pass only: feedback from solids and operations to the flow solution was not iterated.')));
    if (history.length > 1) {
      const live = COUPLED.map((c, k) => ({ c, k })).filter(({ k }) => history.some((hh) => hh[k] !== null));
      summary.append(dataTable({ title: 'Coupling history — the quantities that pass between suites, pass by pass', columns: ['Quantity', 'From suite', ...history.map((_, i) => `Pass ${i + 1}`), 'Last change %'], rows: live.map(({ c, k }) => { const a = history[history.length - 2][k], b = history[history.length - 1][k]; return [`${c[2]} (${c[3]})`, byId(c[0]).short, ...history.map((hh) => hh[k]), a === null || b === null ? '–' : (100 * Math.abs(b - a)) / Math.max(Math.abs(a), Math.abs(b), 1e-9)]; }), note: 'Example of the loop: deposits from suite 4 narrow the bore → suite 3 recomputes pressure drop and slugging → suite 6 recomputes loads → suite 5 changes the inhibitor and operating plan → suite 7 prices the result.' }));
      summary.append(h('div', { class: 'plots' }, pc({ type: 'line', title: 'Convergence of the coupled quantities (each relative to its final value)', xlabel: 'Pass', ylabel: 'Value / final value', series: live.filter(({ k }) => history[history.length - 1][k]).map(({ c, k }) => ({ name: c[2], x: history.map((_, i) => i + 1), y: history.map((hh) => (hh[k] ?? history[history.length - 1][k]) / history[history.length - 1][k]), mode: 'both' })), hlines: [{ y: 1, label: 'final', color: '#64748b' }] })));
    }
    summary.append(h('div', { class: 'row-tools' }, btn('Export all results (JSON)', () => download(JSON.stringify({ case: store.case.name, site: { name: store.case.site.name, lat: store.case.site.lat, lon: store.case.site.lon }, fluid: store.case.fluid, outputs: store.case.outputs }, null, 1), `${store.case.name}_all_results.json`, 'application/json'))));
    summary.append(h('div', { class: 'linkbar' }, 'Next: ', h('a', { class: 'btn mini primary', href: '#/advisor' }, 'See ranked recommendations and the sustainability scorecard →')));
    toast(failed ? 'Integrated run finished with failures.' : 'Integrated run finished.', failed ? 'warn' : 'ok');
  }, 'primary');
  fill(root,
    h('header', { class: 'page-head' }, h('h1', null, 'Integrated run'), h('p', null, 'Solves the suites one after another in data-flow order — fluid → network → flow → solids → operations → integrity → economics — and then repeats the chain so that what the later suites found (deposits narrowing the bore, the inhibitor dose, the operating plan) is fed back to the earlier ones, until the quantities that pass between suites stop changing.')),
    card(h('div', { class: 'row-tools' }, go, stopBtn, h('label', { class: 'inline' }, 'Coupling ', passes), btn('Select all', () => CHAIN.forEach((id) => { sel[id].checked = true; store.pref('chain.' + id, true); }), 'ghost'), btn('Select none', () => CHAIN.forEach((id) => { sel[id].checked = false; store.pref('chain.' + id, false); }), 'ghost'), stat), bar,
      h('div', { class: 'tbl-scroll' }, h('table', { class: 'tbl chain' }, h('thead', null, h('tr', null, ['Run', 'Step', 'Suite', 'Linked inputs', 'Status', 'Key results'].map((x) => h('th', { scope: 'col' }, x)))), rows))),
    summary);
}

// ------------------------------------------------------------------------------- app, offline, mirrors
export function appPage(root, app) {
  const mirrorBox = h('div', { class: 'sources' }), stBox = h('div');
  const check = async () => {
    fill(mirrorBox, MIRRORS.map((m) => h('div', { class: 'source loading', dataset: { url: m.url, kind: m.kind || '' } }, h('b', null, m.name), h('span', null, m.note), h('small', null, m.url), h('em', null, 'checking…'))));
    for (const el of mirrorBox.children) {
      const url = el.dataset.url;
      const arch = el.dataset.kind === 'archive';
      try {
        const ctl = new AbortController(), tm = setTimeout(() => ctl.abort(), 12000);
        if (arch) { await fetch(url, { signal: ctl.signal, cache: 'no-store', credentials: 'omit', mode: 'no-cors' }); clearTimeout(tm); el.className = 'source ok'; el.lastChild.textContent = 'reachable'; }
        else { const r = await fetch(url + 'version.json?t=' + Date.now(), { signal: ctl.signal, cache: 'no-store', credentials: 'omit', mode: 'cors' }); clearTimeout(tm); const j = await r.json(); el.className = 'source ok'; el.lastChild.textContent = `online · build ${String(j.version).slice(0, 12)}`; }
      }
      catch { el.className = 'source fail'; el.lastChild.textContent = navigator.onLine ? 'not reachable from this network' : 'offline'; }
      el.append(h('a', { class: 'btn mini', href: arch ? url + (app.version && app.version !== 'single-file' ? '?v=' + app.version : '') : url + 'index.html', rel: 'noopener', title: arch ? 'Opens the capture of this build if the Archive holds one, otherwise its most recent capture' : null }, 'Open'));
    }
  };
  const paintStorage = async () => {
    let est = null; try { est = await navigator.storage?.estimate?.(); } catch { /* not supported */ }
    const persisted = await navigator.storage?.persisted?.().catch(() => false);
    fill(stBox, kpiGrid([{ label: 'App build', value: app.version || 'dev' }, { label: 'Connection', value: navigator.onLine ? 'online' : 'offline', status: navigator.onLine ? 'ok' : 'warn' }, { label: 'Offline copy', value: app.offlineReady() ? 'ready' : 'preparing…', status: app.offlineReady() ? 'ok' : 'warn' }, { label: 'Installed', value: app.installed() ? 'yes' : 'not yet' }, est ? { label: 'Storage used', value: (est.usage || 0) / 1e6, unit: 'MB' } : null, { label: 'Protected storage', value: persisted ? 'yes' : 'no', help: 'Protected storage is not cleared automatically by the browser when the device runs low on space.' }].filter(Boolean)));
  };
  paintStorage();
  const ios = /iphone|ipad|ipod/i.test(navigator.userAgent);
  fill(root, 
    h('header', { class: 'page-head' }, h('h1', null, 'Install, offline & availability'), h('p', null, `${APP.name} is a self-contained application: after the first visit it runs from this device, with or without a network, and keeps your cases locally.`)),
    card(h('h2', null, 'Install on this device'), stBox,
      h('div', { class: 'row-tools' }, btn('⬇  Install app', () => app.install(), 'primary'), btn('Check for updates', () => app.checkUpdate(true)), btn('Keep my data protected', async () => { const ok = await navigator.storage?.persist?.(); toast(ok ? 'Storage is now protected from automatic clean-up.' : 'The browser did not grant protected storage (install the app first).', ok ? 'ok' : 'warn'); paintStorage(); }, 'ghost'),
        h('a', { class: 'btn', href: 'standalone.html', download: 'HydraSlug-standalone.html', title: 'A single HTML file containing the whole application. Copy it to any computer or phone and open it — no server and no internet needed.' }, '⬇  Download single-file edition')),
      h('ul', { class: 'steps' },
        h('li', null, h('b', null, 'Windows · macOS · Linux · ChromeOS (Chrome, Edge): '), 'press “Install app” above, or the install icon at the right of the address bar.'),
        h('li', null, h('b', null, 'Android (Chrome, Edge, Samsung Internet): '), 'press “Install app”, or menu ⋮ → “Add to Home screen / Install app”.'),
        h('li', { class: ios ? 'hl' : '' }, h('b', null, 'iPhone · iPad (Safari): '), 'tap the Share button, then “Add to Home Screen”. Apple does not allow an install button inside web pages.'),
        h('li', null, h('b', null, 'macOS Safari: '), 'File → “Add to Dock”.'), h('li', null, h('b', null, 'Firefox desktop: '), 'no install prompt, but the app still works offline in a normal tab once loaded; or use the single-file edition.'))),
    card(h('h2', null, 'Works in aeroplane mode'), h('p', null, 'All seven calculation engines, the plotting, file import and your cases run entirely on the device. Only three things need a connection: pulling live site data, looking up published evidence on the decision page, and checking for a newer build. Site data already fetched remain stored with the case. Without a connection the Global site data page answers from a built-in world atlas instead: seabed depth and a local relief grid, seabed temperature and the temperature profile of the water column, sea-surface temperature and salinity, tidal range, currents, wave climate, country, national economic indicators, exchange rates, long-term wind and air temperature, and the oil and gas prices known when the build was made — each value labelled as atlas data with its date.'),
      h('p', { class: 'note' }, 'While you are online the app checks for a newer build in the background and refreshes stored site data that are more than six hours old; installed copies on supporting browsers also refresh periodically in the background.')),
    card(h('h2', null, 'Availability and mirrors'), h('p', null, MIRRORS.length > 1 ? 'The same build is published at more than one independent address. If one host is down, open another — or simply keep using the installed copy, which needs no host at all.' : 'Once installed (or saved as the single-file edition) the application needs no host at all.'), mirrorBox, h('div', { class: 'row-tools' }, btn('Check mirrors now', check))),
    card(h('h2', null, 'Security and privacy'), h('ul', { class: 'steps' },
      h('li', null, 'No account, no tracking, no analytics, no cookies. Cases, inputs and results never leave this device unless you export them.'),
      h('li', null, 'A strict content-security policy blocks inline and third-party scripts; the app loads no external code libraries.'),
      h('li', null, 'Imported files are parsed as data only, size-limited, and never executed; all text is displayed as plain text.'),
      h('li', null, 'Live data requests go over HTTPS only, to a fixed allow-list of public data services, without credentials or referrer.'),
      h('li', null, 'Calculations run in a background thread of your own browser; no simulation data are sent to any server, because there is no server.')),
      h('div', { class: 'row-tools' }, btn('Erase all data stored by this app on this device', async () => { if (!confirm('Erase every stored case, input and preference on this device? This cannot be undone.')) return; try { localStorage.clear(); } catch { /* ignore */ } toast('Local data erased. Reloading…'); setTimeout(() => location.reload(), 600); }, 'ghost danger'))),
    card(h('h2', null, 'Live data sources'), h('div', { class: 'sources' }, SOURCES.map((s) => h('div', { class: 'source' }, h('b', null, s.name), h('span', null, s.gives), h('small', null, s.provider + ' · ' + s.host))))));
  check();
}
