// Neutral boundary-representation model shared by the ACIS and Parasolid readers, with the pipe recogniser.
//
// Model (plain objects, coordinates in model units):
//   { format, version, unitScale (metres per model unit or null),
//     faces: [{ surf: { type: 'plane' | 'cylinder' | 'cone' | 'sphere' | 'torus' | 'spline' | 'other',
//                       o, z, x  (origin / unit axis or normal / unit reference direction),
//                       r        (cylinder radius, cone radius at o, sphere radius, torus minor radius),
//                       R        (torus major radius), tanA (cone: d radius / d distance along +z), ratio (elliptical) },
//               sense (true: face normal = natural surface normal, i.e. away from the axis / centre / tube circle),
//               loops: [{ pts: [[x, y, z] …] (boundary polygon, curved edges sampled), single (a lone vertex) }] }],
//     edges: [{ pts, kind: 'line' | 'circle' | 'ellipse' | 'spline' | 'other', closed }],
//     vertices, counts, skipped, warnings }
//
// recognisePipes() reads pipe runs off the analytic surfaces: every cylinder gives an axis segment (extent of its face
// boundary along the axis), coaxial cylinders are merged into one run carrying all their radii (bore, wall, coating),
// every torus gives a bend (centre, bend radius, angle), and runs whose ends meet are chained into centrelines. Walls that
// do not close round their axis (fillets, grooves) are left out. A model counts as a pipe when its longest centreline is
// at least three diameters long, carries half of all run length and spans half of the model's diagonal.

export const TAU = 2 * Math.PI, ARC_N = 48;
export const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
export const add = (a, b, k = 1) => [a[0] + k * b[0], a[1] + k * b[1], a[2] + k * b[2]];
export const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
export const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
export const vlen = (a) => Math.hypot(a[0], a[1], a[2]);
export const unit = (a) => { const l = vlen(a); return l > 0 ? [a[0] / l, a[1] / l, a[2] / l] : [0, 0, 0]; };
/** A unit vector perpendicular to z, as close to `hint` as possible. */
export const perp = (z, hint) => { let x = hint ? sub(hint, [z[0] * dot(hint, z), z[1] * dot(hint, z), z[2] * dot(hint, z)]) : [0, 0, 0]; if (vlen(x) < 1e-9 * (hint ? vlen(hint) : 1) || !hint) x = cross(Math.abs(z[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0], z); return unit(x); };

/**
 * Points of an elliptical arc: centre c, major-axis vector M, minor-axis vector m, from A to B (points on the curve).
 * forward = true runs in the direction of increasing angle; A = B (or no B) gives the whole ellipse, closing point included.
 */
export function ellipseArc(c, M, m, A, B, forward = true) {
  const mm = dot(M, M), nn = dot(m, m), ang = (P) => { const d = sub(P, c); return Math.atan2(dot(d, m) / nn, dot(d, M) / mm); }, at = (t) => [c[0] + Math.cos(t) * M[0] + Math.sin(t) * m[0], c[1] + Math.cos(t) * M[1] + Math.sin(t) * m[1], c[2] + Math.cos(t) * M[2] + Math.sin(t) * m[2]];
  const ta = A ? ang(A) : 0, full = !A || !B || vlen(sub(A, B)) <= 1e-9 * Math.sqrt(mm);
  let sw = full ? TAU : ang(B) - ta;
  if (!full) { if (forward) { while (sw <= 1e-12) sw += TAU; } else while (sw >= -1e-12) sw -= TAU; } else if (!forward) sw = -TAU;
  const n = Math.max(2, Math.ceil((Math.abs(sw) / TAU) * ARC_N - 1e-9)), out = [];
  for (let i = 0; i <= n; i++) out.push(i === 0 && A ? A : i === n && B && !full ? B : i === n && A && full ? A : at(ta + (sw * i) / n));
  return out;
}

/** Angular interval covered by a set of angles: [start, sweep] with sweep = 2π when no gap wider than `gap` remains. */
export function angularExtent(angles, gap = TAU / 6) {
  if (!angles.length) return [0, TAU];
  const a = angles.map((t) => ((t % TAU) + TAU) % TAU).sort((p, q) => p - q);
  let best = a[0] + TAU - a[a.length - 1], bi = 0;
  for (let k = 1; k < a.length; k++) if (a[k] - a[k - 1] > best) { best = a[k] - a[k - 1]; bi = k; }
  return best < gap ? [0, TAU] : [a[bi], TAU - best];
}

/**
 * Pipe runs of a B-rep model. Returns
 *   runs: [{ kind: 'straight' | 'bend' | 'reducer' (a cone), a, b (centreline ends), radius (largest), radii (all coaxial
 *            radii, ascending; the two end radii of a reducer), length, axis, centre, bendRadius, angle (radians), normal }]
 *   centrelines: [{ pts, closed, radii, length }]  runs chained end to end, longest first (bends sampled as arcs)
 *   diameters: distinct diameters of the main centreline (ascending), isPipe, length (of the main centreline)
 */
export function recognisePipes(model, opts = {}) {
  const faces = model.faces || [], all = [];
  for (const f of faces) for (const l of f.loops) for (const p of l.pts) all.push(p);
  let size = 0;
  if (all.length) { const lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity]; for (const p of all) for (let k = 0; k < 3; k++) { if (p[k] < lo[k]) lo[k] = p[k]; if (p[k] > hi[k]) hi[k] = p[k]; } size = vlen(sub(hi, lo)); }
  const tol = (opts.tolerance ?? 1e-6) * (size || 1), key = (v, q) => v.map((c) => Math.round(c / q)).join(','), lines = new Map(), tori = new Map(), rnd = (v) => +v.toPrecision(10);
  for (const f of faces.slice(0, 2e5)) {
    const s = f.surf, pts = f.loops.flatMap((l) => l.pts);
    if (!s || !pts.length || !(s.r > 0)) continue;
    if ((s.type === 'cylinder' || s.type === 'cone') && (s.ratio === undefined || Math.abs(s.ratio - 1) < 1e-9)) {
      const flip = s.z[0] < -1e-9 || (Math.abs(s.z[0]) <= 1e-9 && (s.z[1] < -1e-9 || (Math.abs(s.z[1]) <= 1e-9 && s.z[2] < 0))), z = flip ? [-s.z[0], -s.z[1], -s.z[2]] : s.z, o = add(s.o, z, -dot(s.o, z));
      const taper = Math.abs(s.tanA || 0) < 1e-9 ? 0 : flip ? -s.tanA : s.tanA, r0 = s.r - dot(s.o, z) * taper, k = key(z, 1e-6) + '|' + key(o, 10 * tol) + (taper ? `|${Math.round(taper * 1e6)}|${Math.round(r0 / (10 * tol))}` : '');
      let v0 = Infinity, v1 = -Infinity;
      for (const p of pts) { const v = dot(p, z); if (v < v0) v0 = v; if (v > v1) v1 = v; }
      let g = lines.get(k);
      if (!g) lines.set(k, (g = { o, z, x: perp(z), parts: [], taper, r0 }));
      const y = cross(g.z, g.x);
      g.parts.push({ r: s.r, v0, v1, ang: pts.map((p) => { const d = sub(p, g.o); return Math.atan2(dot(d, y), dot(d, g.x)); }) });
    } else if (s.type === 'torus' && s.R > s.r) {
      const z = s.z, x = perp(z, s.x), y = cross(z, x), k = key(s.o, 10 * tol) + '|' + key(z.map(Math.abs), 1e-6) + '|' + Math.round(s.R / (10 * tol));
      let g = tori.get(k);
      if (!g) tori.set(k, (g = { o: s.o, z, x, y, R: s.R, tube: new Map(), ang: [] }));
      const tube = g.tube.get(rnd(s.r)) || [];
      for (const p of pts) { const d = sub(p, g.o), u = dot(d, g.x), w = dot(d, g.y); g.ang.push(Math.atan2(w, u)); tube.push(Math.atan2(dot(d, g.z), Math.hypot(u, w) - g.R)); }
      g.tube.set(rnd(s.r), tube);
    }
  }
  const runs = [];
  for (const g of lines.values()) {
    if (g.taper) {                                           // a cone closed round its axis: a reducer between two diameters
      const v0 = Math.min(...g.parts.map((q) => q.v0)), v1 = Math.max(...g.parts.map((q) => q.v1)), ra = Math.abs(g.r0 + v0 * g.taper), rb = Math.abs(g.r0 + v1 * g.taper);
      if (v1 - v0 > 10 * tol && angularExtent(g.parts.flatMap((q) => q.ang))[1] >= TAU) runs.push({ kind: 'reducer', a: add(g.o, g.z, v0), b: add(g.o, g.z, v1), axis: g.z, radius: Math.max(ra, rb), radii: [rnd(Math.min(ra, rb)), rnd(Math.max(ra, rb))], length: v1 - v0 });
      continue;
    }
    // union of the coaxial intervals; each connected stretch is one run with the radii that cover most of it
    const parts = g.parts.sort((a, b) => a.v0 - b.v0);
    for (let i = 0; i < parts.length;) {
      let j = i, v1 = parts[i].v1;
      while (j + 1 < parts.length && parts[j + 1].v0 <= v1 + 10 * tol) { j++; if (parts[j].v1 > v1) v1 = parts[j].v1; }
      const v0 = parts[i].v0, len = v1 - v0, cover = new Map();
      for (let k = i; k <= j; k++) { const r = rnd(parts[k].r), c = cover.get(r) || Object.assign([], { ang: [] }); c.push([parts[k].v0, parts[k].v1]); c.ang.push(...parts[k].ang); cover.set(r, c); }
      const radii = [], minor = [];
      for (const [r, iv] of cover) {
        if (angularExtent(iv.ang)[1] < TAU) continue;       // a fillet or groove: the wall does not close round the axis
        iv.sort((a, b) => a[0] - b[0]); let tot = 0, e = -Infinity; for (const [a, b] of iv) { if (b > e) { tot += b - Math.max(a, e); e = b; } }
        (tot >= 0.5 * len ? radii : minor).push(r);
      }
      radii.sort((a, b) => a - b);
      if (len > 10 * tol && radii.length) runs.push({ kind: 'straight', a: add(g.o, g.z, v0), b: add(g.o, g.z, v1), axis: g.z, radius: radii[radii.length - 1], radii, otherRadii: minor.sort((a, b) => a - b), length: len });
      i = j + 1;
    }
  }
  for (const g of tori.values()) {
    const [t0, sw] = angularExtent(g.ang), radii = [...g.tube].filter(([, a]) => angularExtent(a)[1] >= TAU).map(([r]) => r).sort((a, b) => a - b), at = (t) => [0, 1, 2].map((k) => g.o[k] + g.R * (Math.cos(t) * g.x[k] + Math.sin(t) * g.y[k]));
    if (radii.length) runs.push({ kind: 'bend', a: at(t0), b: at(t0 + sw), radius: radii[radii.length - 1], radii, length: g.R * sw, centre: g.o, bendRadius: g.R, angle: sw, normal: g.z, closed: sw >= TAU - 1e-9, t0, x: g.x, y: g.y });
  }
  // chain runs whose ends coincide (within a fraction of the pipe radius)
  const used = new Uint8Array(runs.length), near = (p, q, r) => vlen(sub(p, q)) <= Math.max(0.05 * r, 100 * tol);
  const sample = (run, rev) => {
    if (run.kind !== 'bend') return rev ? [run.b, run.a] : [run.a, run.b];
    const n = Math.max(2, Math.ceil((run.angle / TAU) * ARC_N - 1e-9)), pts = [];
    for (let i = 0; i <= n; i++) { const t = run.t0 + (run.angle * i) / n; pts.push([0, 1, 2].map((k) => run.centre[k] + run.bendRadius * (Math.cos(t) * run.x[k] + Math.sin(t) * run.y[k]))); }
    return rev ? pts.reverse() : pts;
  };
  const centrelines = [];
  const order = runs.map((_, k) => k).sort((p, q) => runs[q].length - runs[p].length);
  for (const k0 of order) {
    if (used[k0]) continue;
    used[k0] = 1;
    const chain = [[k0, false]];
    if (!runs[k0].closed) for (const front of [false, true]) for (let guard = 0; guard < runs.length; guard++) {
      const [ke, rev] = front ? chain[0] : chain[chain.length - 1], e = runs[ke], end = front ? (rev ? e.b : e.a) : rev ? e.a : e.b;
      let hit = -1, hrev = false, best = Infinity;
      for (let k = 0; k < runs.length; k++) {
        if (used[k] || runs[k].closed) continue;
        const r = Math.max(e.radius, runs[k].radius), da = vlen(sub(runs[k].a, end)), db = vlen(sub(runs[k].b, end));
        if (near(runs[k].a, end, r) && da < best) { hit = k; hrev = front; best = da; }
        if (near(runs[k].b, end, r) && db < best) { hit = k; hrev = !front; best = db; }
      }
      if (hit < 0) break;
      used[hit] = 1;
      if (front) chain.unshift([hit, hrev]); else chain.push([hit, hrev]);
    }
    const pts = [];
    for (const [k, rev] of chain) { const s = sample(runs[k], rev); for (let i = pts.length ? 1 : 0; i < s.length; i++) pts.push(s[i]); }
    const radii = [...new Set(chain.flatMap(([k]) => (runs[k].kind === 'reducer' && chain.length > 1 ? [] : runs[k].radii)))].sort((a, b) => a - b), length = chain.reduce((s, [k]) => s + runs[k].length, 0), closed = !!runs[k0].closed || (chain.length > 2 && near(pts[0], pts[pts.length - 1], runs[k0].radius));
    if (closed && pts.length > 2) pts.pop();
    centrelines.push({ pts, closed, radii, length, runs: chain.map(([k]) => k), maxRadius: Math.max(...chain.map(([k]) => runs[k].radius)) });
  }
  centrelines.sort((a, b) => b.length - a.length);
  const main = centrelines[0], total = runs.reduce((s, r) => s + r.length, 0);
  // a pipe: one chain clearly longer than its diameter that carries most of the cylinder / torus length of the model
  const isPipe = !!main && main.length >= 3 * 2 * main.maxRadius && main.length >= 0.5 * total && main.length >= 0.5 * size;
  for (const r of runs) { delete r.t0; delete r.x; delete r.y; }
  return { runs, centrelines, diameters: main ? main.radii.map((r) => rnd(2 * r)) : [], isPipe, length: main ? main.length : 0 };
}
