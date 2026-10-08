// Suite registry: navigation metadata is static (so the shell renders instantly) and each engine
// is loaded on demand the first time its page is opened.
export const SUITES = [
  { id: 'pvt', num: 1, title: 'Fluid, PVT & Phase Behaviour', short: 'Fluid · PVT', icon: '🧪', blurb: 'Equation-of-state flash, phase envelope, properties, hydrate curve', load: () => import('./s01_pvt.js'), uses: [] },
  { id: 'net', num: 2, title: 'Geometry, Wells, Network & Equipment', short: 'Network', icon: '🛠️', blurb: 'Route and riser profile, wells, network hydraulics, equipment', load: () => import('./s02_net.js'), uses: ['pvt'] },
  { id: 'flow', num: 3, title: 'Multiphase Thermal-Hydraulics & Slugging', short: 'Flow · Slugs', icon: '🌊', blurb: 'Steady and transient multiphase flow, heat loss, slugging', load: () => import('./s03_flow.js'), uses: ['pvt', 'net', 'solids'] },
  { id: 'solids', num: 4, title: 'Hydrate & Multiphase Solids Flow Assurance', short: 'Hydrate · Solids', icon: '❄️', blurb: 'Hydrate kinetics and plugging, wax, scale, asphaltene, sand', load: () => import('./s04_solids.js'), uses: ['pvt', 'net', 'flow', 'ops'] },
  { id: 'ops', num: 5, title: 'Operations, Control & Flow-Assurance Management', short: 'Operations', icon: '🎛️', blurb: 'Shutdown, cooldown, restart, blowdown, pigging, control', load: () => import('./s05_ops.js'), uses: ['pvt', 'net', 'flow', 'solids'] },
  { id: 'integ', num: 6, title: 'Integrity, Loads, Risk & Engineering Assessment', short: 'Integrity · Risk', icon: '🛡️', blurb: 'Slug loads, stress, fatigue, corrosion, erosion, reliability', load: () => import('./s06_integ.js'), uses: ['pvt', 'net', 'flow', 'solids', 'ops'] },
  { id: 'econ', num: 7, title: 'Economics, Techno-Economics & Decision Analysis', short: 'Economics', icon: '💲', blurb: 'CAPEX, OPEX, cash flow, risk cost, uncertainty, decisions', load: () => import('./s07_econ.js'), uses: ['pvt', 'net', 'flow', 'solids', 'ops', 'integ'] },
];
export const byId = (id) => SUITES.find((s) => s.id === id);
/** Forward order of the main data flow: fluid → network → flow → solids → operations → integrity → economics. */
export const CHAIN = ['pvt', 'net', 'flow', 'solids', 'ops', 'integ', 'econ'];
export const downstream = (id) => SUITES.filter((s) => s.uses.includes(id));

const cache = new Map();
export async function loadSuite(id) {
  if (cache.has(id)) return cache.get(id);
  const meta = byId(id);
  if (!meta) throw new Error('Unknown suite: ' + id);
  const mod = (await meta.load()).default;
  cache.set(id, mod);
  return mod;
}
