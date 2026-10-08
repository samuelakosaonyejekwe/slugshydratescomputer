// Route / profile / network / wall-map tests: node tests/route.test.mjs
// Every export of js/core/route.js is checked against shapes with a closed-form answer.
import { profileFrom, drape, minimumCurvature, resample, simplify, segments, checkProfile, networkFrom, mapFrom, terrainFrom, haversine } from '../js/core/route.js';
import { importGeometry, generate } from '../js/core/geom.js';

let pass = 0, fails = 0;
const t0 = Date.now();
const ok = (cond, msg) => { if (cond) { pass++; console.log('ok   ' + msg); } else { fails++; console.log('FAIL ' + msg); } };
const near = (a, b, tol, msg) => ok(typeof a === 'number' && Math.abs(a - b) <= tol, `${msg} (got ${a}, expected ${b} ± ${tol})`);
const group = async (label, fn) => { try { await fn(); } catch (e) { fails++; console.log(`FAIL ${label} — threw: ${e && e.stack ? e.stack.split('\n').slice(0, 3).join(' | ') : e}`); } };
const throws = (label, fn, re) => { try { fn(); ok(false, `${label}: should have thrown`); } catch (e) { ok(e instanceof Error && e.message.length < 400 && (!re || re.test(e.message)), `${label}: clean error "${String(e.message).slice(0, 70)}"`); } };
const enc = new TextEncoder();
const F = (name, text) => { const u8 = enc.encode(text); return { name, size: u8.length, arrayBuffer: async () => u8.buffer.slice(u8.byteOffset, u8.byteOffset + u8.byteLength) }; };
const line = (x, y, z, extra = {}) => ({ kind: 'polylines', name: 'line', format: 'test', polylines: [{ x, y, ...(z ? { z } : {}), closed: false }], bbox: { min: [Math.min(...x), Math.min(...y)], max: [Math.max(...x), Math.max(...y)] }, stats: {}, warnings: [], ...extra });
const table = (headers, rows, stats = {}) => ({ kind: 'table', name: 'table', format: 'test', headers, records: rows.map((r) => Object.fromEntries(headers.map((h, i) => [h, r[i]]))), stats, warnings: [] });
const D = Math.PI / 180, DEG = D * 6371008.8;            // metres per degree of arc on the sphere used (mean Earth radius)

await group('minimumCurvature', async () => {
  const v = minimumCurvature([{ md: 0, inc: 0, azi: 0 }, { md: 1000, inc: 0, azi: 0 }, { md: 2500, inc: 0, azi: 0 }]);
  ok(v.tvd.join() === '0,1000,2500' && v.north.every((q) => q === 0) && v.east.every((q) => q === 0) && v.dls.every((q) => q === 0) && v.z.join() === '0,-1000,-2500', 'vertical well: TVD = MD, no departure, no dog-leg');
  const R = 400, n = 9, st = Array.from({ length: n + 1 }, (_, k) => [R * k * (90 / n) * D, k * (90 / n), 0]), b = minimumCurvature(st);
  near(b.tvd[n], R, 1e-9, '90° build of radius R heading north: TVD = R (closed form)');
  near(b.north[n], R, 1e-9, '90° build: northing = R');
  near(b.east[n], 0, 1e-9, '90° build: no easting');
  ok(b.dls.slice(1).every((q) => Math.abs(q - (30 / R) / D) < 1e-9) && b.dogleg.slice(1).every((q) => Math.abs(q - 10) < 1e-9) && b.dls[0] === 0, `90° build: dog-leg severity = 30 / R rad = ${((30 / R) / D).toFixed(4)}°/30 m at every station`);
  const one = minimumCurvature([[0, 0, 0], [(Math.PI / 2) * R, 90, 0]]);
  near(one.tvd[1], R, 1e-9, 'one 90° interval is still exact (circular arc)');
  const turn = minimumCurvature({ md: st.map((q) => q[0]), inc: st.map(() => 90), azi: st.map((q) => q[1]) });
  near(turn.north[n], R, 1e-9, 'horizontal 90° turn: northing = R');
  near(turn.east[n], R, 1e-9, 'horizontal 90° turn: easting = R');
  ok(Math.abs(turn.tvd[n]) < 1e-9 && turn.x[n] === turn.east[n] && turn.y[n] === turn.north[n], 'horizontal turn stays at constant TVD; x = east, y = north');
  const hold = minimumCurvature([[100, 30, 90], [300, 30, 90]], { tvd: 95, north: 5, east: 7 });
  near(hold.tvd[1], 95 + 200 * Math.cos(30 * D), 1e-9, 'straight hold from a tie-in: TVD advances by ΔMD · cos(inc)');
  near(hold.east[1], 7 + 200 * Math.sin(30 * D), 1e-9, 'straight hold: departure ΔMD · sin(inc) along the azimuth');
  ok(minimumCurvature([[50, 0, 0], [60, 0, 0]]).tvd[0] === 50, 'default tie-in: vertical hole above the first station');
  throws('decreasing MD', () => minimumCurvature([[10, 0, 0], [5, 0, 0]]), /decrease/);
  throws('non-numeric station', () => minimumCurvature([{ md: 0, inc: 'a' }]), /numeric/);
  throws('no stations', () => minimumCurvature([]), /stations/);
});

await group('profileFrom: polylines', async () => {
  const p = profileFrom(line([0, 3], [0, 0], [0, 4]));
  ok(p.x.join() === '0,3' && p.z.join() === '0,4' && p.s.join() === '0,5' && p.length === 5 && p.horizontalLength === 3 && p.zMin === 0 && p.zMax === 4 && p.incl.length === 1, '3-4-5 triangle: x = 3, rise = 4, length = 5');
  near(p.incl[0], Math.atan2(4, 3) / D, 1e-12, '3-4-5 triangle: inclination atan(4/3) = 53.13°, positive uphill');
  ok(typeof p.source === 'string' && p.source.length > 5 && Array.isArray(p.warnings) && p.warnings.length === 0, 'profile names its source and carries warnings');
  const q = profileFrom(line([0, 3, 3, 3], [0, 4, 4, 4], [0, 0, 0, 12]));
  ok(q.x.join() === '0,5,5' && q.s.join() === '0,5,17' && q.incl.join() === '0,90' && q.warnings.some((w) => /repeated/.test(w)), '3-D polyline: plan distance unrolled (3-4-5 in plan), vertical riser 90°, repeated point dropped');
  const down = profileFrom(line([0, 4], [0, 0], [0, -3]));
  near(down.incl[0], -Math.atan2(3, 4) / D, 1e-12, 'downhill element has a negative inclination');
  const rev = profileFrom(line([0, 4], [0, 0], [0, -3]), { reverse: true });
  ok(rev.z.join() === '-3,0' && rev.x.join() === '0,4' && rev.incl[0] > 0, 'reverse: the profile runs the other way and x restarts at 0');
  const xz = profileFrom(line([0, 3, 6], [0, 4, 0]));
  ok(xz.length === 10 && xz.z.join() === '0,4,0' && xz.zMax === 4 && /x–z/.test(xz.source) && xz.warnings.some((w) => /elevation profile/.test(w)), '2-D polyline with monotonic x read as an x–z elevation profile (length 5 + 5)');
  const plan = profileFrom(line([0, 3, 6], [0, 4, 0]), { plane: 'xy' });
  ok(plan.length === 10 && plan.zMax === 0 && plan.horizontalLength === 10 && plan.warnings.some((w) => /flat/.test(w)), 'plane: "xy" reads the same line as a flat plan-view route');
  const loop = profileFrom(line([0, 10, 10, 0], [0, 0, 10, 10]));
  ok(loop.zMax === 0 && loop.length === 30, '2-D polyline that turns back in x is a plan-view route');
  const geo = profileFrom(line([2, 3], [0, 0], [-100, -100], { geographic: true }));
  near(geo.length, DEG, 1e-6, 'geographic polyline: one degree of longitude on the equator in metres');
  near(profileFrom(line([2, 2], [59, 60], [-100, -100], { geographic: true })).length, DEG, 1e-6, 'geographic polyline: one degree of latitude');
  near(profileFrom(line([0, 1], [60, 60], [-100, -100], { geographic: true })).length, haversine(0, 60, 1, 60), 1e-9, 'geographic polyline: one degree of longitude at 60° N (great circle)');
  near(haversine(0, 60, 1, 60), DEG * Math.cos(60 * D), 0.0002 * DEG, 'haversine at 60° N ≈ cos(60°) of the equatorial degree');
  const gz = profileFrom(line([2, 2.01], [60, 60], [-300, -100], { geographic: true }));
  near(gz.length, Math.hypot(haversine(2, 60, 2.01, 60), 200), 1e-9, 'geographic polyline: elevations stay in metres');
  const joined = profileFrom({ kind: 'polylines', polylines: [{ x: [10, 20], y: [0, 0], z: [-5, -5], closed: false }, { x: [0, 10], y: [0, 0], z: [0, -5], closed: false }, { x: [30, 20], y: [0, 0], z: [0, -5], closed: false }, { x: [100, 101], y: [50, 50], z: [0, 0], closed: false }], stats: {} });
  ok(joined.horizontalLength === 30 && joined.x.length === 4 && joined.warnings.some((w) => /3 connected polylines/.test(w)) && joined.warnings.some((w) => /1 of 4 polylines/.test(w)), 'polylines touching end to end are chained (reversing where needed); a stray one is reported');
  const mm = profileFrom(line([0, 3000], [0, 0], [0, 4000], { stats: { units: 'mm' } }));
  ok(Math.abs(mm.length - 5) < 1e-12 && mm.warnings.some((w) => /mm to metres/.test(w)), 'coordinates in millimetres (unit given by the file) converted to metres');
  near(profileFrom(line([0, 3], [0, 0], [0, 4]), { unit: 'km' }).length, 5000, 1e-9, 'opts.unit = "km"');
  near(profileFrom(line([0, 3], [0, 0], [0, 4]), { unit: 'ft' }).length, 5 * 0.3048, 1e-12, 'opts.unit = "ft"');
  ok(profileFrom(line([0, 3], [0, 0], [10, 14]), { depthPositiveDown: true }).z.join() === '-10,-14', 'depthPositiveDown: true flips z');
  throws('unknown unit', () => profileFrom(line([0, 3], [0, 0], [0, 4]), { unit: 'furlong' }), /unit/);
  throws('no geometry', () => profileFrom(null), /geometry/);
  throws('single point', () => profileFrom(line([1, 1], [2, 2], [3, 3])), /two distinct points/);
  throws('elevation grid', () => profileFrom({ kind: 'grid', grid: { x: [0, 1], y: [0, 1], z: [[0, 0], [0, 0]] } }), /drape/);
  throws('voxels', () => profileFrom({ kind: 'voxels' }), /cannot give/);
});

await group('profileFrom: points, networks, tables, meshes', async () => {
  const pts = (rows, extra = {}) => ({ kind: 'points', points: Float64Array.from(rows.flat()), count: rows.length, bbox: { min: [0, 0, 0], max: [0, 0, 0] }, stats: {}, ...extra });
  const c = profileFrom(pts([[0, 0, -10], [30, 40, -10], [30, 40, 20], [60, 80, 20]]));
  ok(c.x.join() === '0,50,50,100' && c.s.join() === '0,50,80,130' && /file order/.test(c.source) && c.warnings.length === 0, 'ordered XYZ centreline followed in file order');
  const sc = profileFrom(pts([[40, 0, -4], [0, 0, 0], [30, 0, -3], [10, 0, -1], [20, 0, -2], [50, 0, -5], [5, 0, -0.5], [45, 0, -4.5], [15, 0, -1.5], [35, 0, -3.5], [25, 0, -2.5]]));
  ok(sc.horizontalLength === 50 && Math.abs(sc.length - Math.hypot(50, 5)) < 1e-9 && sc.warnings.some((w) => /sorted along/.test(w)), 'shuffled centreline points are sorted along their principal axis, with a warning');
  const gp = profileFrom(pts([[2, 60, -300], [2.01, 60, -300], [2.02, 60, -300]], { geographic: true, stats: { role: 'route3d' } }));
  near(gp.length, 2 * haversine(2, 60, 2.01, 60), 1e-6, 'geographic route points: great-circle distances');

  const net = { kind: 'network', network: { nodes: [{ id: 'WH', x: 0, y: 0, z: -1000 }, { id: 'PLEM', x: 3000, y: 0, z: -1000 }, { id: 'RB', x: 3000, y: 4000, z: -1000 }, { id: 'TOP', x: 3000, y: 4000, z: 0 }, { id: 'W2', x: 2900, y: 0, z: -1005 }], edges: [{ from: 'WH', to: 'PLEM', length: 3000 }, { from: 'PLEM', to: 'RB' }, { from: 'RB', to: 'TOP', length: 1100 }, { from: 'W2', to: 'PLEM', length: 120 }] } };
  const np = profileFrom(net);
  ok(np.nodes === undefined && np.s.join() === '0,3000,7000,8100' && np.z.join() === '-1000,-1000,-1000,0' && /WH → TOP over 3 connections/.test(np.source), 'network: longest source-to-sink route WH → PLEM → RB → TOP (declared lengths, coordinate length where none is given)');
  near(np.x[3] - np.x[2], Math.sqrt(1100 ** 2 - 1000 ** 2), 1e-9, 'network: a 1100 m riser over 1000 m of elevation keeps its length (horizontal offset √(L² − Δz²))');
  const ring = { kind: 'network', network: { nodes: ['a', 'b', 'c', 'd'].map((id, k) => ({ id, x: null, y: null, z: -k })), edges: [{ from: 'a', to: 'b', length: 10 }, { from: 'b', to: 'c', length: 10 }, { from: 'c', to: 'a', length: 10 }, { from: 'c', to: 'd', length: 50 }] } };
  const rp = profileFrom(ring);
  ok(rp.length === 60 && rp.warnings.some((w) => /loops/.test(w)), 'cyclic network: longest shortest path (10 + 50), with a warning');
  const noz = profileFrom({ kind: 'network', network: { nodes: [{ id: 'a' }, { id: 'b' }], edges: [{ from: 'a', to: 'b', length: 25 }] } });
  ok(noz.length === 25 && noz.zMin === 0 && noz.warnings.some((w) => /no elevation/.test(w)), 'network without elevations: flat, reported');
  throws('network without lengths', () => profileFrom({ kind: 'network', network: { nodes: [{ id: 'a' }, { id: 'b' }], edges: [{ from: 'a', to: 'b' }] } }), /length/);

  const ce = profileFrom(table(['chainage', 'elevation'], [[0, -100], [300, -100], [300.0001, -100], [600, 300]]));
  near(ce.length, 300 + 0.0001 + Math.hypot(299.9999, 400), 1e-9, 'table chainage + elevation');
  const kp = profileFrom(table(['KP', 'Depth (m)'], [[2, 500], [0, 100], [1, 100]]));
  ok(kp.x.join() === '0,1000,2000' && kp.z.join() === '-100,-100,-500' && kp.warnings.some((w) => /kilometre point/.test(w)) && kp.warnings.some((w) => /sorted by chainage/.test(w)) && kp.warnings.some((w) => /depth/.test(w)), 'table KP + depth: kilometres → metres, depth → negative elevation, rows sorted');
  ok(profileFrom(table(['KP', 'Depth (m)'], [[0, 100], [1, 500]]), { depthPositiveDown: false }).z.join() === '100,500', 'depthPositiveDown: false keeps a depth column as written');
  near(profileFrom(table(['Station (ft)', 'Elev (ft)'], [[0, 0], [300, 400]])).length, 500 * 0.3048, 1e-9, 'table with feet in the headers');
  const xyz = profileFrom(table(['x', 'y', 'z'], [[0, 0, 0], [3, 4, 0], [3, 4, 12]]));
  ok(xyz.s.join() === '0,5,17', 'table x / y / z');
  const en = profileFrom(table(['Easting', 'Northing', 'Elevation'], [[500000, 6500000, -300], [500300, 6500400, -300]]));
  ok(en.length === 500, 'table easting / northing / elevation');
  const lld = profileFrom(table(['Latitude', 'Longitude', 'Depth'], [[0, 2, 100], [0, 3, 100]]));
  near(lld.length, DEG, 1e-6, 'table latitude / longitude / depth: great-circle length');
  ok(lld.zMin === -100 && /geographic/.test(lld.source), 'table latitude / longitude / depth: depth below datum');
  const mt = profileFrom(table(['MD', 'TVD'], [[0, 0], [1000, 1000], [2000, 1600]]));
  ok(mt.length === 2000 && mt.z.join() === '-1600,-1000,0' && mt.horizontalLength === 800 && mt.warnings.some((w) => /bottom-hole to wellhead/.test(w)), 'table MD + TVD: length = MD, departure √(ΔMD² − ΔTVD²), reversed to production direction');
  ok(profileFrom(table(['MD', 'TVD'], [[0, 0], [1000, 1000]]), { reverse: false }).z.join() === '0,-1000', 'reverse: false keeps survey order');
  const R = 500, arc = (Math.PI / 2) * R, mia = profileFrom(table(['MD', 'Inc', 'Azi'], Array.from({ length: 19 }, (_, k) => [(arc * k) / 18, 5 * k, 30])));
  near(mia.length, arc, 1e-9, 'table MD + inclination + azimuth: profile length = measured depth');
  near(mia.zMin, -R, 1e-9, 'table MD / inc / azi: 90° build reaches TVD = R');
  near(mia.horizontalLength, R, 1, 'table MD / inc / azi: horizontal run of the quarter circle = R (to the 5° station spacing)');
  const pos = profileFrom(table(['a', 'b'], [[0, 0], [3, 4]]));
  ok(pos.length === 5 && pos.warnings.some((w) => /not recognised/.test(w)), 'unnamed two-column table: distance, elevation (with a warning)');
  ok(profileFrom(table(['p', 'q', 'r'], [[0, 0, 0], [3, 4, 0], [3, 4, 12]])).length === 17, 'unnamed three-column table: x, y, z');
  throws('table with no usable columns', () => profileFrom(table(['name', 'value'], [['a', 1], ['b', 2]])), /no columns/);
  throws('table with one row', () => profileFrom(table(['chainage', 'elevation'], [[0, 0]])), /two rows/);

  // straight pipe of radius 0.5 m from (0, 0, 0) to (30, 0, 40): 3-4-5 axis, length 50
  const tri = [], ax = [0.6, 0, 0.8], u = [0, 1, 0], v = [-0.8, 0, 0.6], nL = 100, nT = 16, P = (i, k) => { const a = (2 * Math.PI * k) / nT, s = (50 * i) / nL; return [0, 1, 2].map((c) => ax[c] * s + 0.5 * (Math.cos(a) * u[c] + Math.sin(a) * v[c])); };
  for (let i = 0; i < nL; i++) for (let k = 0; k < nT; k++) tri.push(...P(i, k), ...P(i + 1, k), ...P(i + 1, k + 1), ...P(i, k), ...P(i + 1, k + 1), ...P(i, k + 1));
  const mp = profileFrom({ kind: 'mesh', triangles: tri, count: tri.length / 9, bbox: { min: [0, 0, 0], max: [30, 1, 40] }, stats: {} });
  ok(mp.warnings.some((w) => /ESTIMATED/.test(w)) && /mesh/.test(mp.source), 'pipe mesh: the centreline is labelled as an estimate');
  near(mp.length, 50, 2.5, 'pipe mesh: centreline length ≈ 50 (ends shortened by up to half a slice)');
  ok(mp.incl.every((a) => Math.abs(a - Math.atan2(4, 3) / D) < 0.5), 'pipe mesh: inclination 53.13° ± 0.5° everywhere');
  near(mp.zMax - mp.zMin, 40, 2.5, 'pipe mesh: rise ≈ 40');
  const again = profileFrom(mp);
  ok(again.length === mp.length && again.x.length === mp.x.length, 'a profile object can be passed back in');
});

await group('profileFrom: size cap and generators', async () => {
  const n = 5000, x = Array.from({ length: n + 1 }, (_, i) => (20000 * i) / n), z = x.map((v) => -500 + 30 * Math.sin(v / 400) + 8 * Math.sin(v / 37));
  const full = profileFrom(line(x, x.map(() => 0), z), { maxPoints: 1e6 }), capped = profileFrom(line(x, x.map(() => 0), z)), small = profileFrom(line(x, x.map(() => 0), z), { maxPoints: 60 });
  ok(full.x.length === n + 1 && capped.x.length === 400 && small.x.length === 60 && capped.warnings.some((w) => /simplified from 5001 to 400/.test(w)), 'maxPoints: 5001 points capped to 400 (default) and to 60');
  ok(capped.zMin === full.zMin && capped.zMax === full.zMax && capped.horizontalLength === full.horizontalLength && capped.x[0] === 0, 'capping keeps both ends and the highest and lowest point');
  ok(capped.length <= full.length && capped.length > 0.985 * full.length && small.length > 0.95 * full.length, `capping loses little length (${full.length.toFixed(1)} → ${capped.length.toFixed(1)} → ${small.length.toFixed(1)})`);
  const Dw = 1500, th = 80 * D, t = Math.tan(th), a = Dw / (Math.sqrt(1 + t * t) - 1), cr = generate({ type: 'catenary riser', waterDepth: Dw, topAngle: 10, flowline: 0, n: 2000 }), cp = profileFrom(cr, { maxPoints: 5000 });
  near(cp.length, a * t, 1e-5 * a * t, 'catenary riser generator: profile length = a · sinh(X / a) = a · tan θ (closed form)');
  near(cp.horizontalLength, a * Math.asinh(t), 1e-9 * a, 'catenary riser generator: horizontal reach a · asinh(tan θ)');
  ok(cp.zMin === -Dw && cp.zMax === 0 && cp.incl.every((v, i) => i === 0 || v >= cp.incl[i - 1] - 1e-9) && cp.incl[0] < 0.2, 'catenary riser generator: leaves the seabed tangentially and steepens monotonically');
  near(cp.incl[cp.incl.length - 1], 80, 0.05, 'catenary riser generator: top inclination 80° (10° from vertical)');
  const mid = cp.x.findIndex((v) => v > 0.5 * cp.horizontalLength);
  near(cp.z[mid], -Dw + a * (Math.cosh(cp.x[mid] / a) - 1), 1e-6, 'catenary riser generator: z = a (cosh(x / a) − 1) above the seabed');
  const wl = generate({ type: 'well trajectory (build-hold)', kickoff: 400, buildRate: 2, inclination: 60, md: 2400 }), wp = profileFrom(wl);
  near(wp.length, 2400, 1e-9, 'well generator: profile length = total measured depth');
  ok(wp.z[0] === wp.zMin && wp.z[wp.z.length - 1] === 0 && Math.abs(wp.zMin + wl.stats.tvd) < 1e-9 && wp.incl[wp.incl.length - 1] === 90 && Math.abs(wp.incl[0] - 30) < 1e-6, 'well generator: runs bottom-hole → wellhead, vertical at the top, 90° − 60° at the bottom');
  const jp = profileFrom(generate({ type: 'jumper (M-shape)', span: 30, height: 8, dip: 3 }));
  near(jp.length, 52, 1e-9, 'jumper generator: length 2h + 2·dip + span');
  ok(jp.incl.join() === '90,0,-90,0,90,0,-90', 'jumper generator: vertical and horizontal elements only');
  const lz = profileFrom(generate({ type: 'lazy-wave riser' }));
  ok(lz.incl.some((v) => v < -1) && lz.incl[lz.incl.length - 1] > 70 && Math.abs(lz.length - generate({ type: 'lazy-wave riser' }).stats.length) < 0.5, 'lazy-wave generator: descends after the hog bend and hangs off steeply');
});

await group('drape', async () => {
  const xs = Array.from({ length: 21 }, (_, i) => i * 100), ys = Array.from({ length: 11 }, (_, j) => j * 100), plane = (x, y) => -100 + 0.01 * x + 0.02 * y;
  const grid = { kind: 'grid', name: 'plane', grid: { x: xs, y: ys, z: ys.map((y) => xs.map((x) => plane(x, y))), geographic: false }, bbox: { min: [0, 0, -100], max: [2000, 1000, -60] }, stats: {} };
  const route = line([0, 1000, 1000], [0, 0, 500]), d = drape(route, grid), p = d.polylines[0];
  ok(d.kind === 'polylines' && d.polylines.length === 1 && d.stats.draped && p.x.length === 16 && p.x[0] === 0 && p.x[10] === 1000 && p.y[15] === 500 && d.warnings.length === 0, 'drape: route densified to the grid spacing (10 + 5 steps)');
  ok(p.z.every((z, i) => Math.abs(z - plane(p.x[i], p.y[i])) < 1e-9), 'drape over a planar grid reproduces the plane exactly (bilinear)');
  const pr = profileFrom(d);
  near(pr.length, Math.hypot(1000, 10) + Math.hypot(500, 10), 1e-9, 'draped route profile: length of the two sloping legs');
  near(pr.incl[0], Math.atan(0.01) / D, 1e-9, 'draped route profile: first leg slope atan(0.01)');
  near(pr.incl[pr.incl.length - 1], Math.atan(0.02) / D, 1e-9, 'draped route profile: second leg slope atan(0.02)');
  const viaOpt = profileFrom(route, { grid });
  ok(Math.abs(viaOpt.length - pr.length) < 1e-12 && /draped/.test(viaOpt.source), 'profileFrom(route, { grid }) drapes a flat route first');
  const off = drape({ x: [50, 250], y: [-50, 150] }, { x: xs, y: ys, elev: grid.grid.z });
  ok(off.polylines[0].z[0] === plane(50, 0) && off.stats.outside === 1 && off.warnings.some((w) => /outside the grid/.test(w)) && Math.abs(off.polylines[0].z[off.polylines[0].z.length - 1] - plane(250, 150)) < 1e-9, 'drape: plain { x, y } route and { x, y, elev } terrain; points outside take the edge value');
  const holes = { x: [0, 1, 2], y: [0, 1], z: [[1, NaN, 3], [1, NaN, 3]] }, dh = drape({ x: [0, 2], y: [0.5, 0.5] }, holes);
  ok(dh.polylines[0].z.join() === '1,3' && dh.warnings.some((w) => /without data/.test(w)), 'drape: samples on cells without data are left out');
  const gg = { kind: 'grid', grid: { x: [2, 2.5, 3], y: [59, 60, 61], z: [[-100, -150, -200], [-100, -150, -200], [-100, -150, -200]], geographic: true } };
  const dg = drape(line([2, 3], [60, 60], null, { geographic: true }), gg), pg = profileFrom(dg);
  ok(dg.geographic && pg.zMin === -200 && pg.zMax === -100 && Math.abs(pg.horizontalLength - haversine(2, 60, 2.5, 60) - haversine(2.5, 60, 3, 60)) < 1e-6, 'drape: geographic route on a geographic grid, distances in metres');
  throws('projected route on a geographic grid', () => drape(line([0, 1], [0, 1]), gg), /coordinate system/);
  throws('drape without a grid', () => drape(route, { x: [0], y: [0], z: [[0]] }), /grid/);
  throws('drape without a route', () => drape({}, grid), /route/);
  throws('route off a hole-only grid', () => drape({ x: [0, 1], y: [0, 1] }, { x: [0, 1], y: [0, 1], z: [[NaN, NaN], [NaN, NaN]] }), /valid cell/);
});

await group('resample, simplify, segments', async () => {
  const p = profileFrom(line([0, 3], [0, 0], [0, 4])), r = resample(p, 6);
  ok(r.x.length === 6 && r.s.every((s, i) => Math.abs(s - i) < 1e-12) && Math.abs(r.length - 5) < 1e-12 && Math.abs(r.x[2] - 1.2) < 1e-12 && Math.abs(r.z[2] - 1.6) < 1e-12, 'resample: 6 nodes equally spaced along a 3-4-5 line (every 1 m of pipe)');
  const bend = profileFrom(line([0, 100, 100], [0, 0, 0], [0, 0, 100])), rb = resample(bend, 5);
  ok(rb.s.join() === '0,50,100,150,200' && rb.x.join() === '0,50,100,100,100' && rb.z.join() === '0,0,0,50,100', 'resample: spacing follows arc length around a corner');
  ok(resample(bend, 1).x.length === 2 && resample({ x: [0, 10], z: [0, 0] }, 3).x.join() === '0,5,10', 'resample: at least two nodes; accepts a bare { x, z }');
  const seg = segments(bend);
  ok(seg.length === 2 && JSON.stringify(seg[0]) === JSON.stringify({ x0: 0, x1: 100, z0: 0, z1: 0, length: 100, incl: 0 }) && seg[1].incl === 90 && seg[1].z1 === 100 && seg[1].length === 100, 'segments: rows { x0, x1, z0, z1, length, incl }');
  near(segments(p)[0].incl, 53.13010235415598, 1e-12, 'segments: 3-4-5 inclination');
  const zig = profileFrom(line([0, 10, 20, 30, 40, 50, 60], [0, 0, 0, 0, 0, 0, 0], [0, 0.05, 0, 5, 0, -0.02, 0]));
  ok(simplify(zig, 0).x.length === 7 && simplify(zig, 0.5).x.join() === '0,20,30,40,50,60' && simplify(zig, 0.5).z.includes(-0.02), 'simplify: bumps below the tolerance go, the peak stays, and the lowest point is kept even though it is tiny');
  const s6 = simplify(zig, 100);
  ok(s6.x.join() === '0,30,50,60' && s6.zMax === 5 && s6.zMin === -0.02, 'simplify with a huge tolerance keeps the ends, the highest and the lowest point');
  const col = profileFrom(line([0, 1, 2, 3, 4], [0, 0, 0, 0, 0], [0, 2, 4, 6, 8]));
  ok(simplify(col, 1e-9).x.join() === '0,4' && Math.abs(simplify(col, 1e-9).length - col.length) < 1e-12, 'simplify: collinear points removed without changing the length');
  throws('simplify without a profile', () => simplify({ x: [0] }, 1), /profile/);
  throws('segments without a profile', () => segments(null), /profile/);
});

await group('checkProfile', async () => {
  const p = profileFrom(generate({ type: 'undulating flowline', length: 3000, amplitude: 12, wavelength: 300 })), ck = checkProfile(p, { id: 0.254 }), byName = Object.fromEntries(ck.map((c) => [c.name, c]));
  ok(ck.length === 8 && ck.every((c) => c.pass === true && typeof c.name === 'string' && 'got' in c && 'expected' in c && typeof c.note === 'string'), `a generated flowline passes all ${ck.length} checks: ${ck.map((c) => c.name).join(', ')}`);
  near(byName['Internal volume'].got, (Math.PI / 4) * 0.254 ** 2 * p.length, 1e-9, 'internal volume = π/4 · id² · length');
  ok(byName['Length conservation'].got === byName['Length conservation'].expected || Math.abs(byName['Length conservation'].got - p.length) < 1e-9, 'length conservation: Σ element length = s at the last node');
  ok(Math.abs(byName['Round trip'].got - p.length) <= 1e-9 * p.length, 'round trip profile → table → profile reproduces the length to 1e-9');
  ok(checkProfile(p).length === 7 && !checkProfile(p).some((c) => c.name === 'Internal volume'), 'without a diameter the volume check is left out, not passed');
  const bad = { x: [0, 10, 10, 5, 20], z: [0, 0, 0, 1, NaN], s: [0, 10, 10, 15, 30], incl: [0, 0, 170, 0] }, cb = Object.fromEntries(checkProfile(bad, { id: 0.1 }).map((c) => [c.name, c.pass]));
  ok(!cb['No zero-length elements'] && !cb['No duplicate nodes'] && !cb['Monotonic chainage'] && !cb['Finite values'] && !cb['Length conservation'] && !cb['Inclination bounds'] && !cb['Round trip'], 'a broken profile (repeated node, back-tracking, NaN, wrong s) fails the matching checks');
  const steep = { x: [0, 0, 5], z: [0, -10, -10], s: [0, 10, 15], incl: [-90, 0] };
  ok(checkProfile(steep).every((c) => c.pass), 'vertical downward element is within bounds');
  throws('checkProfile without a profile', () => checkProfile({}), /profile/);
});

await group('networkFrom', async () => {
  const g = await importGeometry(F('field.csv', 'from,to,length,diameter,type,name\nWH1,MAN,2600,0.2,flowline,FL-1\nWH2,MAN,900,0.2,flowline,FL-2\nMAN,FPSO,7500,0.3,riser,R-1\nMAN,FPSO,7500,0.3,riser,R-2\nT1,T1,0,0.1,pipe,STUB\nX,Y,5,0.1,pipe,ISLAND\n'));
  const n = networkFrom(g);
  ok(n.nodes.length === 7 && n.edges.length === 6 && Object.keys(n.nodes[0]).join() === 'id,type,x,y,z' && Object.keys(n.edges[0]).join() === 'from,to,type,length,diameter,name' && n.edges[2].type === 'riser' && n.edges[0].diameter === 0.2, 'network geometry → nodes { id, type, x, y, z } and edges { from, to, type, length, diameter, name }');
  ok(n.issues.duplicates.join() === 'R-2' && n.issues.zeroLength.join() === 'STUB' && n.issues.disconnected.sort().join() === 'T1,X,Y' && n.issues.loops === 2, `issues: duplicate R-2, zero-length STUB, 3 disconnected nodes, 2 loops (${JSON.stringify(n.issues)})`);
  const pl = { kind: 'polylines', polylines: [{ x: [0, 50, 100], y: [0, 10, 0], z: [-10, -10, -10], closed: false }, { x: [100, 100], y: [0, 80], z: [-10, -20], closed: false }, { x: [100, 0], y: [80, 0], closed: false }, { x: [500, 510], y: [500, 500], closed: false }], stats: {} }, np = networkFrom(pl);
  ok(np.nodes.length === 5 && np.edges.length === 4 && np.edges[0].from === 'N1' && np.edges[0].to === 'N2' && np.edges[1].from === 'N2' && np.edges[2].to === 'N1' && np.issues.loops === 1 && np.issues.disconnected.length === 2, 'polylines → edges with coincident end points welded (triangle = 1 loop, stray line disconnected)');
  near(np.edges[0].length, 2 * Math.hypot(50, 10), 1e-12, 'polyline edge length follows the vertices');
  near(np.edges[1].length, Math.hypot(80, 10), 1e-12, 'polyline edge length is 3-D');
  ok(np.nodes[1].x === 100 && np.nodes[1].z === -10 && np.edges[3].name === 'L4' && np.edges[0].diameter === null, 'welded node keeps its coordinates; unknown diameter is null');
  const ring = networkFrom({ kind: 'polylines', polylines: [{ x: [0, 10, 10, 0], y: [0, 0, 10, 10], closed: true }], stats: {} });
  ok(ring.nodes.length === 1 && ring.edges[0].from === ring.edges[0].to && ring.edges[0].length === 40 && ring.issues.loops === 1, 'closed polyline: one node, one 40 m loop');
  const geo = networkFrom({ kind: 'polylines', geographic: true, polylines: [{ x: [2, 3], y: [0, 0], closed: false }], stats: {} });
  near(geo.edges[0].length, DEG, 1e-6, 'geographic polyline edge length in metres');
  const tb = networkFrom(table(['From', 'To', 'Length', 'Diameter', 'Type'], [['A', 'B', 100, 0.2, 'Pipe'], ['B', 'C', 50, 0.2, 'Valve'], ['', 'C', 1, 1, 'x']]));
  ok(tb.nodes.length === 3 && tb.edges.length === 2 && tb.edges[1].type === 'valve' && tb.edges[1].length === 50 && tb.issues.loops === 0 && tb.issues.disconnected.length === 0, 'table with from / to columns → network');
  throws('table without from / to', () => networkFrom(table(['a', 'b'], [[1, 2]])), /from and to/);
  throws('mesh', () => networkFrom({ kind: 'mesh' }), /network/);
});

await group('mapFrom', async () => {
  const wt = mapFrom(table(['x (m)', 'theta (deg)', 'wt (mm)'], [[0, 0, 12.7], [0, 90, 12.5], [0, 180, 9.5], [0, 270, 12.6], [2, 0, 12.7], [2, 90, 12.4], [2, 180, 10.1], [2, 270, 12.6]]));
  ok(wt.kind === 'thickness' && wt.x.join() === '0,2' && wt.theta.join() === '0,90,180,270' && wt.t.length === 4 && wt.t[0].length === 2 && wt.t[2][0] === 9.5 && wt.t[2][1] === 10.1 && wt.min === 9.5 && wt.max === 12.7 && wt.unit === 'mm' && wt.missing === 0, 'wall-thickness table → t[θ][x], min / max, unit from the header');
  near(wt.mean, (12.7 + 12.5 + 9.5 + 12.6 + 12.7 + 12.4 + 10.1 + 12.6) / 8, 1e-12, 'map mean');
  const ck = mapFrom(table(['KP', 'Clock Position', 'Pit Depth (mm)'], [[1, 12, 0.5], [1, 3, 1.5], [1.5, 6, 2.5]]));
  ok(ck.kind === 'corrosion' && ck.x.join() === '1000,1500' && ck.theta.join() === '0,90,180' && ck.t[2][1] === 2.5 && Number.isNaN(ck.t[2][0]) && ck.missing === 3 && ck.max === 2.5, 'corrosion map: clock positions → degrees (12 o’clock = 0°), KP → metres, gaps are NaN');
  const dep = mapFrom(table(['x', 'theta', 't', 'delta'], [[0, 0, 0, 0], [0, 180, 0, 0], [100, 0, 0, 0], [100, 180, 0, 0], [0, 0, 3600, 0.4], [0, 180, 3600, 1.2], [100, 0, 3600, 0.1], [100, 180, 3600, 2.5]]));
  ok(dep.kind === 'deposit' && dep.time === 3600 && dep.times.join() === '0,3600' && dep.max === 2.5 && dep.t[1][1] === 2.5 && dep.t[0][0] === 0.4, 'time-resolved deposit map (x, θ, t, δ): the latest time step by default');
  ok(mapFrom(table(['x', 'theta', 't', 'delta'], [[0, 0, 0, 0.1], [0, 0, 5, 0.9]]), { time: 0 }).max === 0.1, 'opts.time selects a time step');
  const sheet = mapFrom(table(['x', 'y', 'thickness'], [[0, 0, 10], [0, 0.5, 11], [1, 0, 12], [1, 0.5, 13]]));
  ok(sheet.kind === 'thickness' && sheet.thetaUnit === 'length' && sheet.theta.join() === '0,0.5' && sheet.t[1][1] === 13, 'unrolled sheet (x, y, thickness): y kept as a circumferential length');
  const fromGrid = mapFrom({ kind: 'grid', grid: { x: [0, 1, 2], y: [0, 180], z: [[5, 6, 7], [8, NaN, 4]] }, stats: { role: 'depositMap' } });
  ok(fromGrid.kind === 'deposit' && fromGrid.min === 4 && fromGrid.max === 8 && fromGrid.missing === 1 && fromGrid.t[1][0] === 8 && mapFrom({ kind: 'grid', grid: { x: [0, 1], y: [0, 1], z: [[1, 2], [3, 4]] } }, { kind: 'corrosion' }).kind === 'corrosion', 'grid → map (x axial, y circumferential); kind from the role or from opts');
  const fromPts = mapFrom({ kind: 'points', points: Float64Array.from([0, 0, 12, 0, 90, 11, 5, 0, 12.5, 5, 90, 10]), count: 4, stats: {} });
  ok(fromPts.x.join() === '0,5' && fromPts.theta.join() === '0,90' && fromPts.t[1][1] === 10 && fromPts.kind === 'thickness', 'points (x, θ, value) on a lattice → map');
  throws('table without a value column', () => mapFrom(table(['x', 'theta', 'colour'], [[0, 0, 1]])), /value column/);
  throws('mesh', () => mapFrom({ kind: 'mesh' }), /wall map/);
});

await group('terrainFrom', async () => {
  const xs = Array.from({ length: 5 }, (_, i) => i * 10), ys = [0, 10, 20], g = { kind: 'grid', name: 'bathy', grid: { x: xs, y: ys, z: ys.map((y) => xs.map((x) => -50 - x - 0.5 * y)), geographic: false }, bbox: { min: [0, 0, -100], max: [40, 20, -50] } };
  const t = terrainFrom(g);
  ok(t.x.join() === xs.join() && t.y.join() === ys.join() && t.elev[2][4] === -100 && t.elev[0][0] === -50 && t.geographic === false && t.filled === 0 && t.name === 'bathy' && t.elev !== g.grid.z, 'small grid passed through as { x, y, elev[y][x], geographic }');
  const big = { kind: 'grid', grid: { x: Array.from({ length: 201 }, (_, i) => i), y: Array.from({ length: 101 }, (_, j) => j), z: Array.from({ length: 101 }, (_, j) => Array.from({ length: 201 }, (_q, i) => -i - 2 * j)), geographic: true } }, tb = terrainFrom(big, 21);
  ok(tb.x.length === 21 && tb.y.length === 21 && tb.geographic === true && tb.x[20] === 200 && tb.y[20] === 100 && Math.abs(tb.elev[10][10] - (-100 - 100)) < 1e-9 && Math.abs(tb.elev[20][20] + 400) < 1e-9, 'large grid resampled bilinearly to n × n (a plane stays a plane)');
  const holed = terrainFrom({ kind: 'grid', grid: { x: [0, 1, 2], y: [0, 1, 2], z: [[-1, -1, -1], [-1, NaN, -1], [-3, NaN, -3]] } });
  ok(holed.filled === 2 && holed.elev[1][1] === -1 && holed.elev[2][1] === -3, 'cells without data are filled from their valid neighbours and counted');
  const pts = [], rnd = (() => { let s = 7; return () => (s = (s * 16807) % 2147483647) / 2147483647; })();
  for (let i = 0; i < 4000; i++) { const x = rnd() * 1000, y = rnd() * 500; pts.push(x, y, -200 + 0.05 * x - 0.1 * y); }
  const tp = terrainFrom({ kind: 'points', points: Float64Array.from(pts), count: 4000, bbox: { min: [0, 0, -250], max: [1000, 500, -150] } }, 20);
  let worst = 0; for (let j = 0; j < 20; j++) for (let i = 0; i < 20; i++) worst = Math.max(worst, Math.abs(tp.elev[j][i] - (-200 + 0.05 * tp.x[i] - 0.1 * tp.y[j])));
  ok(tp.x.length === 20 && tp.y.length === 20 && worst < 1.5, `scattered soundings on a plane gridded to 20 × 20 (worst error ${worst.toFixed(3)} m)`);
  throws('polylines', () => terrainFrom({ kind: 'polylines' }), /elevation grid/);
  throws('grid of holes', () => terrainFrom({ kind: 'grid', grid: { x: [0, 1], y: [0, 1], z: [[NaN, NaN], [NaN, NaN]] } }), /no valid/);
});

await group('imported files end to end', async () => {
  const prof = await importGeometry(F('profile.csv', 'KP (km),Elevation (m)\n0,-1000\n3,-1000\n3.0,-1000\n7,-1300\n7.4,-1000\n')), p = profileFrom(prof);
  near(p.length, 3000 + Math.hypot(4000, 300) + Math.hypot(400, 300), 1e-9, 'CSV KP / elevation → profile length in metres');
  ok(p.incl[p.incl.length - 1] > 36 && p.incl[1] < 0 && checkProfile(p, { id: 0.3 }).every((c) => c.pass), 'CSV profile: downhill then a steep climb; all checks pass');
  const rt = await importGeometry(F('route.csv', 'Easting,Northing,Elevation\n500000,6500000,-300\n500300,6500400,-300\n500300,6500400,-100\n')), pr = profileFrom(rt);
  ok(pr.s.join() === '0,500,700' && pr.incl[1] === 90, 'CSV easting / northing / elevation → 3-D route profile');
  const sv = await importGeometry(F('well.dev', 'MD INCL AZIM\n0 0 0\n1000 0 0\n1500 0 0\n')), ps = profileFrom(sv);
  ok(ps.length === 1500 && ps.z.join() === '-1500,-1000,0' && ps.incl.every((v) => v === 90), 'vertical well survey → 1500 m vertical profile, bottom-hole first');
  const inp = await importGeometry(F('riser.inp', '*NODE\n1, 0., 0., -100.\n2, 30., 40., -100.\n3, 30., 40., 0.\n*ELEMENT, TYPE=PIPE31\n1, 1, 2\n2, 2, 3\n')), pi = profileFrom(inp);
  ok(pi.s.join() === '0,50,150' && pi.zMin === -100, 'Abaqus PIPE31 string → profile');
  const nw = await importGeometry(F('net.json', JSON.stringify({ nodes: [{ id: 'A', x: 0, y: 0, z: -500 }, { id: 'B', x: 1200, y: 0, z: -500 }, { id: 'C', x: 1200, y: 0, z: 0 }], edges: [{ from: 'A', to: 'B' }, { from: 'B', to: 'C' }] }))), pn = profileFrom(nw);
  ok(pn.length === 1700 && pn.incl.join() === '0,90' && networkFrom(nw).issues.loops === 0, 'JSON network → longest-route profile');
  const dem = await importGeometry(F('bathy.asc', `ncols 5\nnrows 3\nxllcorner 0\nyllcorner 0\ncellsize 100\n${Array.from({ length: 3 }, () => [-300, -310, -320, -330, -340].join(' ')).join('\n')}\n`));
  const route = await importGeometry(F('route.xy', '50 150\n450 150\n')), pd = profileFrom(route, { grid: dem });
  near(pd.length, Math.hypot(400, 40), 1e-9, 'x-y route file draped over an ESRI ASCII grid: 400 m across a 10 % slope');
  ok(terrainFrom(dem).elev[1][2] === -320 && pd.zMin === -340 && pd.zMax === -300, 'the same grid as terrain; draped profile spans the sampled depths');
  const wt = await importGeometry(F('wt.csv', 'x (m),clock,wt (mm)\n0,12,12.7\n0,6,9.9\n1,12,12.6\n1,6,10.2\n')), m = mapFrom(wt);
  ok(m.theta.join() === '0,180' && m.min === 9.9 && m.kind === 'thickness', 'CSV wall-thickness map → mapFrom');
});

console.log(`\n${pass} checks passed, ${fails} failed, ${((Date.now() - t0) / 1000).toFixed(1)} s`);
process.exit(fails ? 1 : 0);
