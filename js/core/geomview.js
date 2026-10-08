// Geometry hub UI: shared by the data portal and by the "Geometry" tab of every suite.
// Imported geometry of any pathway (CAD, surface, mesh, drawing, GIS / bathymetry, point cloud, well
// survey, network, tables, procedural) is previewed (pipeline profile, well plan and section, network
// layout, wall maps), measured, and routed to the inputs of the suite that can use it.
import { h, clear, btn, kpiGrid, dataTable, toast, importBtn, help } from './ui.js';
import { plotCard } from './plot.js';
import { fmt } from './num.js';
import { download, LIMITS } from './io.js';
import { FORMATS, PATHWAYS, SUITE_GEOMETRY, formatOf, importGeometry, sectionOf, gridOf, microstructure, networkSummary, dimensions, generate } from './geom.js';
import { geometryLinks, routeProfile } from './geomlinks.js';
import { mapFrom, networkFrom, checkProfile } from './route.js';

const attached = new Map(); // suite id -> Geometry (kept in memory for the session; files are never uploaded)
export const geometryOf = (id) => attached.get(id);
export function attachGeometry(id, g) { if (g) attached.set(id, g); else attached.delete(id); }

export const ACCEPT = [...new Set(FORMATS.flatMap((f) => f.ext.map((e) => '.' + e)))].join(',');
const pc = (spec) => plotCard(spec, { onDownload: download });
const safe = (fn, fb = null) => { try { return fn(); } catch { return fb; } };

/** Everything that can be measured from a geometry, computed once per import. */
export function derive(g) {
  if (g._derived) return g._derived;
  const d = { dims: safe(() => dimensions(g)), micro: null, net: null, grid: null };
  if (g.kind === 'voxels') d.micro = safe(() => microstructure(g));
  if (g.kind === 'network') d.net = safe(() => networkSummary(g));
  if (g.kind === 'grid') { let lo = Infinity, hi = -Infinity; for (const r of g.grid.z) for (const v of r) if (Number.isFinite(v)) { if (v < lo) lo = v; if (v > hi) hi = v; } d.grid = lo <= hi ? { min: lo, max: hi, nx: g.grid.x.length, ny: g.grid.y.length } : null; }
  d.profile = safe(() => routeProfile(g, {}));                         // pipeline / riser / flowline elevation profile, when the geometry is a route
  if (g.kind === 'table' && (g.stats?.role === 'thicknessMap' || g.stats?.role === 'depositMap')) d.map = safe(() => mapFrom(g));
  if (g.kind === 'network') d.network = safe(() => networkFrom(g));
  Object.defineProperty(g, '_derived', { value: d, enumerable: false });
  return d;
}

/** Read one or several files (companions such as .shp + .dbf or .mhd + .raw are matched by base name). */
export async function readFiles(fileList, opts = {}) {
  const files = [...fileList];
  if (!files.length) throw new Error('No file selected.');
  for (const f of files) if (f.size > LIMITS.fileBytes) throw new Error(`${f.name} is larger than ${LIMITS.fileBytes / 1e6} MB.`);
  const rank = (f) => { const m = formatOf(f.name); return m ? (m.support === 'convert' ? 1 : 0) : 2; };
  const companionExt = /\.(shx|dbf|prj|cpg|raw|hdr|img|bin|mtl)$/i;
  const main = files.filter((f) => !companionExt.test(f.name) || files.length === 1).sort((a, b) => rank(a) - rank(b))[0] || files[0];
  const companion = Object.fromEntries(files.filter((f) => f !== main).map((f) => [f.name, f]));
  return importGeometry(main, { ...opts, companion });
}

const km = (p) => (p.horizontalLength > 5000 ? 1e-3 : 1);
/** Elevation-versus-distance plot of a pipeline profile. */
export function profilePlot(p, title = 'Elevation profile') {
  const k = km(p);
  return pc({ type: 'line', title: `${title} — ${fmt(p.length, 5)} m along the pipe`, xlabel: k === 1 ? 'Horizontal distance (m)' : 'Horizontal distance (km)', ylabel: 'Elevation (m)', series: [{ name: 'Pipe centreline', x: p.x.map((v) => v * k), y: p.z }] });
}
/** Schematic node positions for a network without coordinates: columns by distance from the inlets. */
function schematic(net) {
  const rank = new Map(net.nodes.map((n) => [n.id, 0])), indeg = new Map(net.nodes.map((n) => [n.id, 0]));
  for (const e of net.edges) if (indeg.has(e.to) && e.from !== e.to) indeg.set(e.to, indeg.get(e.to) + 1);
  const q = net.nodes.filter((n) => !indeg.get(n.id)).map((n) => n.id);
  for (let i = 0, guard = 0; i < q.length && guard < 2e5; i++) for (const e of net.edges) { if (e.from !== q[i] || !rank.has(e.to) || e.to === e.from) continue; guard++; if (rank.get(e.to) < rank.get(q[i]) + 1 && rank.get(q[i]) < net.nodes.length) { rank.set(e.to, rank.get(q[i]) + 1); q.push(e.to); } }
  const col = new Map(), pos = new Map();
  for (const n of net.nodes) { const r = rank.get(n.id), k = col.get(r) || 0; col.set(r, k + 1); pos.set(n.id, { x: r, y: -k }); }
  return pos;
}
function preview(g) {
  const out = [], d = derive(g);
  if (g.survey) {
    const s = g.survey, hd = s.north.map((n, i) => Math.hypot(n - s.north[0], s.east[i] - s.east[0]));
    out.push(pc({ type: 'line', title: 'Well path — plan view', xlabel: 'East (m)', ylabel: 'North (m)', equal: true, series: [{ name: 'Well path', x: s.east, y: s.north }, { name: 'Wellhead', x: [s.east[0]], y: [s.north[0]], mode: 'points', size: 5 }] }));
    out.push(pc({ type: 'line', title: 'Well path — vertical section', xlabel: 'Horizontal departure (m)', ylabel: 'Elevation (m, TVD negative down)', series: [{ name: 'Well path', x: hd, y: s.tvd.map((v) => 0 - v) }] }));
    out.push(pc({ type: 'line', title: 'Inclination and dog-leg severity', xlabel: 'Measured depth (m)', ylabel: 'Inclination (°) / DLS (°/30 m)', series: [{ name: 'Inclination', x: s.md, y: s.inc }, { name: 'Dog-leg severity', x: s.md, y: s.dls, dash: true }] }));
    return out;
  }
  if (d.profile) out.push(profilePlot(d.profile, g.kind === 'network' ? 'Elevation profile of the longest route' : 'Elevation profile'));
  if (d.map) {
    const m = d.map, lo = Number.isFinite(m.min) ? m.min : 0;
    out.push(pc({ type: 'field', title: `${{ thickness: 'Wall thickness', corrosion: 'Metal-loss depth', deposit: 'Deposit thickness' }[m.kind]} map${m.times.length > 1 ? ` — latest of ${m.times.length} time steps` : ''}`, xlabel: 'Axial position', ylabel: m.thetaUnit === 'deg' ? 'Circumferential position (°)' : 'Circumferential position', zlabel: m.unit || m.kind, x: m.x, y: m.theta, z: m.t.map((r) => r.map((v) => (Number.isFinite(v) ? v : lo))), ...(m.missing ? { mask: m.t.map((r) => r.map((v) => !Number.isFinite(v))) } : {}), cmap: 'viridis', interpolate: false }));
  }
  if (g.kind === 'grid') out.push(pc({ type: 'field', title: `Elevation grid — ${g.name}`, xlabel: g.grid.geographic ? 'Longitude (°)' : 'x', ylabel: g.grid.geographic ? 'Latitude (°)' : 'y', zlabel: 'Elevation', x: g.grid.x, y: g.grid.y, z: g.grid.z.map((r) => Array.from(r)), cmap: d.grid && d.grid.max <= 0 ? 'sea' : d.grid && d.grid.min >= 0 ? 'land' : 'viridis', contours: 8 }));
  else if (g.kind === 'voxels') {
    const v = g.voxels, k = Math.floor(v.nz / 2), sx = Math.max(1, Math.ceil(v.nx / 220)), sy = Math.max(1, Math.ceil(v.ny / 160));
    const xs = [], ys = [], z = [];
    for (let i = 0; i < v.nx; i += sx) xs.push(i * v.spacing[0]);
    for (let j = 0; j < v.ny; j += sy) { ys.push(j * v.spacing[1]); const row = []; for (let i = 0; i < v.nx; i += sx) row.push(v.data[(k * v.ny + j) * v.nx + i] ? 1 : 0); z.push(row); }
    out.push(pc({ type: 'field', title: `Voxel volume${v.nz > 1 ? ` — slice ${k + 1} of ${v.nz}` : ''} (1 = solid)`, xlabel: 'x', ylabel: 'y', zlabel: 'solid', x: xs, y: ys, z, cmap: 'viridis', equal: true, zmin: 0, zmax: 1 }));
  } else if (g.kind === 'points') {
    const n = g.count, step = Math.max(1, Math.ceil(n / 4000)), x = [], y = [], route = g.stats?.role === 'route3d';
    for (let i = 0; i < n; i += step) { x.push(g.points[3 * i]); y.push(g.points[3 * i + 1]); }
    out.push(pc({ type: 'line', title: route ? 'Route — plan view' : `Point cloud — plan view (${Math.min(n, x.length).toLocaleString()} of ${n.toLocaleString()} points)`, xlabel: g.geographic ? 'Longitude (°)' : 'x', ylabel: g.geographic ? 'Latitude (°)' : 'y', series: [{ name: route ? 'Route' : 'Points', x, y, ...(route ? {} : { mode: 'points', size: 1.5 }) }] }));
    if (g.seabed) out.push(pc({ type: 'line', title: 'Seabed estimate along the seismic line (first strong arrival, 1500 m/s)', xlabel: 'Distance along line (m)', ylabel: 'Elevation (m)', series: [{ name: 'Estimated seabed', x: g.seabed.distance, y: g.seabed.depth.map((v) => (Number.isFinite(v) ? 0 - v : NaN)) }] }));
    const gr = route || g.seabed ? null : safe(() => gridOf(g, 60, 45));
    if (gr) out.push(pc({ type: 'field', title: 'Surface gridded from the points', xlabel: 'x', ylabel: 'y', zlabel: 'z', x: gr.x, y: gr.y, z: gr.z, cmap: 'viridis', contours: 8 }));
  } else if (g.kind === 'network') {
    const placed = g.network.nodes.some((n) => Number.isFinite(n.x) && Number.isFinite(n.y) && (n.x || n.y)), sch = placed ? null : schematic(g.network), at = (id) => (sch ? sch.get(id) : nodes.get(id));
    const nodes = new Map(g.network.nodes.map((n) => [n.id, n])), x = [], y = [];
    for (const e of g.network.edges) { const a = at(e.from), b = at(e.to); if (a && b && [a.x, a.y, b.x, b.y].every(Number.isFinite)) { x.push(a.x, b.x, NaN); y.push(a.y, b.y, NaN); } }
    const shown = g.network.nodes.map((n) => at(n.id)).filter((q) => q && Number.isFinite(q.x) && Number.isFinite(q.y));
    if (shown.length) out.push(pc({ type: 'line', title: placed ? 'Network layout (plan view)' : 'Network layout (schematic: the file gives no node coordinates)', xlabel: placed ? 'x' : 'Steps from the inlets', ylabel: placed ? 'y' : '', series: [{ name: 'Connections', x, y }, { name: 'Nodes / equipment', x: shown.map((q) => q.x), y: shown.map((q) => q.y), mode: 'points', size: 4 }] }));
  } else if (g.kind === 'mesh' || g.kind === 'polylines') {
    const flat = g.kind === 'polylines' && g.polylines.every((p) => p.z && p.y.every((v) => v === p.y[0]));    // x-z profile only: the plan view is a straight line
    const segs = flat ? [] : safe(() => sectionOf(g), []) || [], x = [], y = [];
    segs.slice(0, 8000).forEach((s) => { x.push(s[0], s[2], NaN); y.push(s[1], s[3], NaN); });
    if (segs.length) out.push(pc({ type: 'line', title: g.kind === 'mesh' ? 'Mid-plane section' : d.profile ? 'Route — plan view' : 'Outlines', xlabel: g.geographic ? 'Longitude (°)' : 'x', ylabel: g.geographic ? 'Latitude (°)' : 'y', series: [{ name: 'Geometry', x, y }] }));
  }
  return out;
}

function facts(g) {
  const d = derive(g), k = [{ label: 'Format', value: g.format || '—' }, { label: 'Pathway', value: PATHWAYS[g.pathway]?.title || g.pathway || '—' }, { label: 'Read as', value: g.survey ? 'Well trajectory' : { mesh: 'Triangulated surface', polylines: g.polylines?.some((p) => p.z) ? '3-D lines / route' : 'Lines / outlines', points: g.stats?.role === 'route3d' ? 'Ordered route points' : 'Point cloud', grid: 'Elevation grid', voxels: 'Voxel volume', network: 'Pipeline / equipment network', table: { profile: 'Pipeline profile table', route3d: 'Route table', thicknessMap: 'Wall-thickness / corrosion map', depositMap: 'Deposit map', timeseries: 'Time series', particles: 'Particle table' }[g.stats?.role] || 'Data table' }[g.kind] || g.kind }];
  if (g.count) k.push({ label: g.kind === 'points' ? 'Points' : 'Triangles', value: g.count });
  if (g.polylines) k.push({ label: 'Lines', value: g.polylines.length });
  for (const [key, val] of Object.entries(g.stats || {})) if ((typeof val === 'number' || typeof val === 'string') && !/^(nx|ny|nz|zmin|zmax|min|max|polylines|role)$/i.test(key)) k.push({ label: key.replace(/([A-Z])/g, ' $1').replace(/^./, (c) => c.toUpperCase()), value: val });
  if (d.dims?.size) k.push({ label: 'Extent', value: d.dims.size.filter((s) => Number.isFinite(s)).map((s) => fmt(s, 4)).join(' × ') });
  if (d.dims?.area && g.kind === 'mesh') k.push({ label: 'Surface area', value: d.dims.area }); if (d.dims?.volume) k.push({ label: 'Volume', value: d.dims.volume, help: d.dims.closed === false ? 'The surface is not closed, so the volume is approximate.' : '' });
  if (d.micro) k.push({ label: 'Porosity', value: d.micro.porosity }, { label: 'Tortuosity', value: (d.micro.tortuosity || []).filter(Number.isFinite).map((t) => fmt(t, 3)).join(' / ') || '—' }, { label: 'Mean pore size', value: d.micro.meanPoreSize }, { label: 'Specific surface', value: d.micro.specificSurface, unit: '1/length' }, { label: 'Percolates', value: (d.micro.percolates || []).map((p) => (p ? 'yes' : 'no')).join(' / ') });
  if (d.net) k.push({ label: 'Nodes', value: d.net.nodes }, { label: 'Connections', value: d.net.edges }, { label: 'Total pipe length', value: d.net.totalLength, unit: 'm' });
  if (d.network) k.push({ label: 'Disconnected nodes', value: d.network.issues.disconnected.length }, { label: 'Loops', value: d.network.issues.loops });
  if (d.profile) k.push({ label: 'Route length', value: d.profile.length, unit: 'm' }, { label: 'Horizontal length', value: d.profile.horizontalLength, unit: 'm' }, { label: 'Lowest point', value: d.profile.zMin, unit: 'm' }, { label: 'Highest point', value: d.profile.zMax, unit: 'm' }, { label: 'Steepest element', value: d.profile.incl.reduce((q, v) => Math.max(q, Math.abs(v)), 0), unit: '°' });
  if (d.map) k.push({ label: 'Map minimum', value: d.map.min, unit: d.map.unit || '' }, { label: 'Map maximum', value: d.map.max, unit: d.map.unit || '' }, { label: 'Map mean', value: d.map.mean, unit: d.map.unit || '' });
  if (d.grid) k.push({ label: 'Lowest point', value: d.grid.min, unit: 'm' }, { label: 'Highest point', value: d.grid.max, unit: 'm' }, { label: 'Grid', value: `${d.grid.nx} × ${d.grid.ny}` });
  return k.slice(0, 24);
}

/** Card describing one imported geometry: facts, warnings, preview and network tables. */
export function geometryCard(g) {
  const d = derive(g), box = h('div');
  box.append(kpiGrid(facts(g)));
  const notes = [...(g.warnings || []), ...(d.profile?.warnings || []).map((w) => 'Profile: ' + w)];
  if (notes.length) box.append(h('div', { class: 'warns' }, notes.slice(0, 10).map((w) => h('div', { class: 'warn warn' }, h('b', null, 'Note'), ' ', String(w)))));
  const pv = preview(g);
  if (pv.length) box.append(h('div', { class: 'plots' }, pv));
  if (d.profile) {
    const ck = safe(() => checkProfile(d.profile), []) || [];
    if (ck.length) box.append(dataTable({ title: `Profile checks — ${d.profile.source}`, columns: ['Check', 'Result', 'Got', 'Expected', 'What it means'], rows: ck.map((c) => [c.name, c.pass ? 'pass' : 'FAIL', c.got, c.expected, c.note]) }));
  }
  if (d.net) {
    box.append(dataTable({ title: 'Equipment inventory', columns: ['Type', 'Count'], rows: Object.entries(d.net.byType || {}).map(([t, n]) => [t, n]) }));
    if (d.net.pipes?.length) box.append(dataTable({ title: 'Connections', columns: ['From', 'To', 'Length (m)', 'Diameter (m)', 'Name'], rows: d.net.pipes.slice(0, 400).map((p) => [String(p.from), String(p.to), p.length ?? null, p.diameter ?? null, p.name || '']) }));
  }
  if (g.kind === 'table') box.append(dataTable({ title: 'Data preview', columns: g.headers, rows: g.records.slice(0, 15).map((r) => g.headers.map((hd) => r[hd])) }));
  return box;
}

/** Controls that build a procedural geometry (no file needed). */
export function generatorPanel(onGeometry) {
  const defs = {
    catenary: { label: 'Catenary riser (seabed to surface)', p: { waterDepth: [1000, 'Water depth (m)'], topAngle: [12, 'Hang-off angle from vertical (°)'], flowline: [250, 'Seabed lead-in length (m)'], n: [200, 'Points'] } },
    lazywave: { label: 'Lazy-wave riser', p: { waterDepth: [1000, 'Water depth (m)'], topAngle: [12, 'Hang-off angle from vertical (°)'], sagHeight: [120, 'Sag-bend height above seabed (m)'], hogHeight: [250, 'Hog-bend height above seabed (m)'], buoyancyRatio: [1, 'Net uplift / submerged weight (-)'], flowline: [250, 'Seabed lead-in length (m)'], n: [240, 'Points'] } },
    flowline: { label: 'Undulating flowline (random hills)', p: { length: [5000, 'Horizontal length (m)'], amplitude: [20, 'Hill amplitude (m)'], wavelength: [500, 'Hill wavelength (m)'], slope: [0, 'Overall slope (°, + uphill)'], waterDepth: [1000, 'Water depth at the start (m)'], seed: [1, 'Random seed'] } },
    well: { label: 'Well trajectory (build and hold)', p: { kickoff: [500, 'Kick-off depth (m MD)'], buildRate: [3, 'Build rate (°/30 m)'], inclination: [45, 'Hold inclination (°)'], md: [3000, 'Total measured depth (m)'], azimuth: [0, 'Azimuth (°)'], step: [30, 'Station spacing (m)'] } },
    jumper: { label: 'Jumper (M-shape)', p: { span: [30, 'Hub-to-hub span (m)'], height: [8, 'Leg height (m)'], dip: [3.2, 'Dip of the middle section (m)'], leg: [0.2, 'Shoulder length (fraction of span)'], waterDepth: [1000, 'Hub depth (m)'] } },
    spheres: { label: 'Random packed grains (sand or hydrate-particle bed)', p: { porosity: [0.45, 'Target porosity'], radius: [0.06, 'Grain radius (fraction of box)'], n: [96, 'Resolution (cells)'], seed: [7, 'Random seed'] } },
    spacer: { label: 'Channel with periodic cylindrical obstacles (2-D section)', p: { height: [0.05, 'Channel height (m)'], length: [0.8, 'Length (m)'], diameter: [0.02, 'Obstacle diameter (m)'], pitch: [0.2, 'Obstacle pitch (m)'] }, sel: { arrangement: ['zigzag', 'cavity', 'submerged'] } },
    lattice: { label: 'Crossed-wire screen (sand-screen unit cell)', p: { pitch: [0.003, 'Pitch (m)'], diameter: [0.0004, 'Wire diameter (m)'], angle: [90, 'Crossing angle (°)'], n: [96, 'Resolution (cells)'] } },
    voronoi: { label: 'Voronoi foam / cellular deposit structure', p: { cells: [30, 'Number of cells'], wall: [0.03, 'Wall thickness (fraction of box)'], n: [96, 'Resolution (cells)'], seed: [3, 'Random seed'] } },
    tpms: { label: 'Triply-periodic porous structure', p: { period: [1, 'Period'], porosity: [0.6, 'Target porosity'], n: [64, 'Resolution (cells)'] }, sel: { surface: ['gyroid', 'schwarzP', 'diamond'] } },
  };
  const sel = h('select', { 'aria-label': 'Generator type' }, Object.entries(defs).map(([k, v]) => h('option', { value: k }, v.label))), form = h('div', { class: 'fields' });
  const paint = () => { const d = defs[sel.value]; clear(form); for (const [k, [val, label]] of Object.entries(d.p)) form.append(h('div', { class: 'field' }, h('label', null, label), h('div', { class: 'ctl' }, h('input', { type: 'number', step: 'any', value: val, dataset: { k } })))); for (const [k, opts] of Object.entries(d.sel || {})) form.append(h('div', { class: 'field' }, h('label', null, k[0].toUpperCase() + k.slice(1)), h('div', { class: 'ctl' }, h('select', { dataset: { k } }, opts.map((o) => h('option', { value: o }, o)))))); };
  sel.addEventListener('change', paint); paint();
  return h('div', null, h('div', { class: 'row-tools' }, sel, btn('Generate geometry', () => {
    try { const spec = { type: sel.value }; for (const el of form.querySelectorAll('[data-k]')) spec[el.dataset.k] = el.tagName === 'SELECT' ? el.value : +el.value; const g = generate(spec); g.name ||= defs[sel.value].label; g.format ||= 'Procedural generator'; g.pathway ||= 'procedural'; onGeometry(g); }
    catch (e) { toast('Could not generate: ' + e.message, 'bad'); }
  }, 'primary')), form);
}

/** Supported-format catalogue grouped by pathway. */
export function formatCatalogue() {
  const badge = { full: ['reads directly', 'ok'], partial: ['reads (partial)', 'warn'], convert: ['convert first', ''] };
  return h('div', null, Object.entries(PATHWAYS).map(([id, p]) => {
    const list = FORMATS.filter((f) => f.pathway === id);
    if (!list.length) return null;
    return h('details', { class: 'ref' }, h('summary', null, `${p.title} — ${list.filter((f) => f.support !== 'convert').length} read directly, ${list.filter((f) => f.support === 'convert').length} via conversion`), h('p', null, p.blurb),
      h('div', { class: 'tbl-scroll' }, h('table', { class: 'tbl' }, h('thead', null, h('tr', null, ['Format', 'Extensions', 'Support', 'What happens'].map((x) => h('th', null, x)))),
        h('tbody', null, list.map((f) => h('tr', null, h('td', { class: 'lead' }, f.name), h('td', { class: 'lead' }, f.ext.map((e) => '.' + e).join(' ')), h('td', { class: 'lead' }, h('span', { class: 'badge ' + badge[f.support][1] }, badge[f.support][0])), h('td', { class: 'lead wrap' }, f.support === 'convert' ? f.convert : f.note)))))));
  }));
}

/**
 * The Geometry tab of a suite. api = { fields, values(), setValue(key, value), rerender() }.
 */
export function geometryTab(suite, api) {
  const info = SUITE_GEOMETRY[suite.id] || { classes: '', accepts: [] }, box = h('div', { class: 'groups' }), body = h('div');
  const paint = () => {
    const g = attached.get(suite.id);
    clear(body);
    if (!g) return body.append(h('p', { class: 'note' }, 'No geometry attached to this suite yet. Import a file here, generate one, or send one from the Data portal. Without geometry the suite uses the pipeline, well and equipment dimensions typed on its Inputs tab.'));
    const items = safe(() => geometryLinks(suite, g, derive(g)), []) || [];
    body.append(h('h3', null, g.name || 'Geometry'), geometryCard(g));
    if (items.length) body.append(h('fieldset', { class: 'group' }, h('legend', null, 'What this suite takes from the geometry'),
      h('ul', { class: 'linklist' }, items.map((it) => h('li', null, h('b', null, api.fields.find((f) => f.key === it.key)?.label || it.key), ' ← ', typeof it.value === 'number' ? fmt(it.value) : Array.isArray(it.value) ? `${it.value.length} rows` : 'geometry data', h('small', null, ' · ' + it.from)))),
      h('div', { class: 'row-tools' }, btn('Apply to inputs', () => { for (const it of items) api.setValue(it.key, it.value); toast(`${items.length} input${items.length > 1 ? 's' : ''} set from the geometry.`, 'ok'); }, 'primary'))));
    else body.append(h('p', { class: 'note' }, `This geometry was read correctly, but ${suite.short || suite.title} has no input that a “${g.kind}” geometry can set directly. Its measured dimensions above can be typed into the inputs, or send it from the Data portal to a suite that uses this class of geometry (routes and profiles: suites 2 and 3; deposit maps: suite 4; networks: suites 2 and 5; wall-thickness maps and structural meshes: suite 6).`));
    body.append(h('div', { class: 'row-tools' }, btn('Remove geometry', () => { attached.delete(suite.id); for (const f of api.fields.filter((x) => x.type === 'file')) api.setValue(f.key, null); paint(); }, 'ghost')));
  };
  const take = (g) => { attached.set(suite.id, g); const items = safe(() => geometryLinks(suite, g, derive(g)), []) || []; for (const it of items) api.setValue(it.key, it.value); toast(`${g.name} attached${items.length ? ` — ${items.length} input${items.length > 1 ? 's' : ''} set` : ''}.`, 'ok'); paint(); };
  const inp = h('input', { type: 'file', hidden: true, multiple: true, accept: ACCEPT });
  inp.addEventListener('change', async () => { const fl = [...inp.files]; inp.value = ''; if (!fl.length) return; try { take(await readFiles(fl)); } catch (e) { toast(e.message || 'Could not read the file.', 'bad', 12000); } });
  paint();
  box.append(
    h('p', { class: 'summary' }, 'Geometry for this suite: ', info.classes || 'parametric dimensions on the Inputs tab.'),
    h('fieldset', { class: 'group' }, h('legend', null, 'Attach geometry', help('Files are read on this device only. Select companion files together (for example .shp with .dbf and .prj, or an OpenFOAM points file with faces, owner and neighbour).')),
      h('div', { class: 'row-tools' }, inp, btn('Import geometry file…', () => inp.click(), 'primary'), h('a', { class: 'btn', href: '#/portal' }, 'Open the Data portal')),
      h('details', { class: 'ref' }, h('summary', null, 'Or generate a riser, flowline, well or jumper'), generatorPanel(take))),
    body,
    h('details', { class: 'ref' }, h('summary', null, 'All supported geometry formats'), formatCatalogue()));
  return box;
}
