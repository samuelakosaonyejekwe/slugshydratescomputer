// Parasolid transmit-file reader (text .x_t) producing the neutral B-rep model of fmt_brep.js, written against the
// published "Parasolid XT Format Reference".
//
// Read:   the keyword header, the "T" text stream of modeller versions 14 and later (node records laid out by the base
//         schema SCH_13006 plus the field edits embedded in the file, so node types and fields added by newer kernels are
//         followed or skipped exactly), and schema-13006 files without embedded edits; the topology body → region →
//         shell → face → loop → fin → edge → vertex → point; edge curves LINE, CIRCLE, ELLIPSE, INTERSECTION (its chart
//         of points) and TRIMMED_CURVE; surfaces PLANE, CYLINDER, CONE, SPHERE and TORUS with face and surface senses.
//         Parasolid models are in metres.
// Not read: B_SURFACE, OFFSET_SURF, SWEPT_SURF, SPUN_SURF, BLENDED_EDGE and foreign surfaces (their faces are counted in
//         `skipped`, the boundary is still delivered), B_CURVE and SP_CURVE edges (chord between the vertices), instance
//         transforms of assemblies (bodies are given in their own coordinates), attributes, partition and delta files,
//         files of modellers before version 13 (their schemas are not embedded) and binary .x_b files.
// File content is untrusted: counts, lengths and indices are bounds-checked and the node total is capped.
import { sub, cross, vlen, unit, perp, ellipseArc } from './fmt_brep.js';

function fail(msg) { const e = new Error(msg); e.user = true; throw e; }

// Base schema SCH_13006: effective fields of each node type, as "name:type[count]" (count 1 = variable length).
// Types: d n u w p f t one number; c char; l logical; v h vector; i interval; b box.
const BASE = {
  10: 'highest_node_id:d attributes_groups:p attribute_chains:p list:p surface:p curve:p point:p key:p res_size:f res_linear:f ref_instance:p next:p previous:p state:u owner:p type:u sub_instance:p',
  11: 'node_id:d attributes_groups:p type:u part:p transform:p assembly:p next_in_part:p prev_in_part:p next_of_part:p prev_of_part:p',
  12: 'highest_node_id:d attributes_groups:p attribute_chains:p surface:p curve:p point:p key:p res_size:f res_linear:f ref_instance:p next:p previous:p state:u owner:p body_type:u nom_geom_state:u shell:p boundary_surface:p boundary_curve:p boundary_point:p region:p edge:p vertex:p',
  13: 'node_id:d attributes_groups:p body:p next:p face:p edge:p vertex:p region:p front_face:p',
  14: 'node_id:d attributes_groups:p tolerance:f next:p previous:p loop:p shell:p surface:p sense:c next_on_surface:p previous_on_surface:p next_front:p previous_front:p front_shell:p',
  15: 'node_id:d attributes_groups:p fin:p face:p next:p',
  16: 'node_id:d attributes_groups:p tolerance:f fin:p previous:p next:p curve:p next_on_curve:p previous_on_curve:p owner:p',
  17: 'attributes_groups:p loop:p forward:p backward:p vertex:p other:p edge:p curve:p next_at_vx:p sense:c',
  18: 'node_id:d attributes_groups:p fin:p previous:p next:p point:p tolerance:f owner:p',
  19: 'node_id:d attributes_groups:p body:p next:p previous:p shell:p type:c',
  29: 'node_id:d attributes_groups:p owner:p next:p previous:p pvec:v',
  30: 'node_id:d attributes_groups:p owner:p next:p previous:p geometric_owner:p sense:c pvec:v direction:v',
  31: 'node_id:d attributes_groups:p owner:p next:p previous:p geometric_owner:p sense:c centre:v normal:v x_axis:v radius:f',
  32: 'node_id:d attributes_groups:p owner:p next:p previous:p geometric_owner:p sense:c centre:v normal:v x_axis:v major_radius:f minor_radius:f',
  38: 'node_id:d attributes_groups:p owner:p next:p previous:p geometric_owner:p sense:c surface:p[2] chart:p start:p end:p',
  40: 'base_parameter:f base_scale:f chart_count:d chordal_error:f angular_error:f parameter_error:f[2] hvec:h[1]',
  41: 'type:c hvec:h[1]',
  45: 'vertices:f[1]',
  50: 'node_id:d attributes_groups:p owner:p next:p previous:p geometric_owner:p sense:c pvec:v normal:v x_axis:v',
  51: 'node_id:d attributes_groups:p owner:p next:p previous:p geometric_owner:p sense:c pvec:v axis:v radius:f x_axis:v',
  52: 'node_id:d attributes_groups:p owner:p next:p previous:p geometric_owner:p sense:c pvec:v axis:v radius:f sin_half_angle:f cos_half_angle:f x_axis:v',
  53: 'node_id:d attributes_groups:p owner:p next:p previous:p geometric_owner:p sense:c centre:v radius:f axis:v x_axis:v',
  54: 'node_id:d attributes_groups:p owner:p next:p previous:p geometric_owner:p sense:c centre:v axis:v major_radius:f minor_radius:f x_axis:v',
  56: 'node_id:d attributes_groups:p owner:p next:p previous:p geometric_owner:p sense:c blend_type:c surface:p[2] spine:p range:f[2] thumb_weight:f[2] boundary:p[2] start:p end:p',
  59: 'node_id:d attributes_groups:p owner:p next:p previous:p geometric_owner:p sense:c boundary:n blend:p',
  60: 'node_id:d attributes_groups:p owner:p next:p previous:p geometric_owner:p sense:c check:c true_offset:l surface:p offset:f scale:f',
  67: 'node_id:d attributes_groups:p owner:p next:p previous:p geometric_owner:p sense:c section:p sweep:v scale:f',
  68: 'node_id:d attributes_groups:p owner:p next:p previous:p geometric_owner:p sense:c profile:p base:v axis:v start:v end:v start_param:f end_param:f x_axis:v scale:f',
  70: 'node_id:d owner:p next:p previous:p list_type:d list_length:d block_length:d size_of_entry:d list_block:p finger_block:p finger_index:d notransmit:l',
  74: 'n_entries:d next_block:p entries:p[1]',
  79: 'string:c[1]',
  80: 'next:p identifier:p type_id:d actions:u[8] field_names:p legal_owners:l[14] fields:u[1]',
  81: 'node_id:d definition:p owner:p next:p previous:p next_of_type:p previous_of_type:p fields:p[1]',
  82: 'values:d[1]', 83: 'values:f[1]', 84: 'values:c[1]', 85: 'values:v[1]', 86: 'values:v[1]', 87: 'values:v[1]', 88: 'values:t[1]', 89: 'values:v[1]', 98: 'values:w[1]',
  90: 'node_id:d attributes_groups:p owner:p next:p previous:p type:u first_member:p',
  91: 'dummy_node_id:d owning_group:p owner:p next:p previous:p next_member:p previous_member:p',
  100: 'node_id:d owner:p next:p previous:p rotation_matrix:f[9] translation_vector:v scale:f flag:d perspective_vector:v',
  101: 'assembly:p attribute:p body:p transform:p surface:p curve:p point:p alive:l attrib_def:p highest_id:d current_id:d',
  102: 'string:c[1]',
  120: 'node_id:d attributes_groups:p owner:p next:p previous:p geometric_owner:p sense:c type:c data:p tf:p internal_geom:p[1]',
  121: 'geom_type:d real_array:p int_array:p', 122: 'key:p real_array:p int_array:p',
  124: 'node_id:d attributes_groups:p owner:p next:p previous:p geometric_owner:p sense:c nurbs:p data:p',
  125: 'original_uint:i original_vint:i extended_uint:i extended_vint:i self_int:u original_u_start:c original_u_end:c original_v_start:c original_v_end:c extended_u_start:c extended_u_end:c extended_v_start:c extended_v_end:c analytic_form_type:c swept_form_type:c spun_form_type:c blend_form_type:c analytic_form:p swept_form:p spun_form:p blend_form:p',
  126: 'u_periodic:l v_periodic:l u_degree:n v_degree:n n_u_vertices:d n_v_vertices:d u_knot_type:u v_knot_type:u n_u_knots:d n_v_knots:d rational:l u_closed:l v_closed:l surface_form:u vertex_dim:n bspline_vertices:p u_knot_mult:p v_knot_mult:p u_knots:p v_knots:p',
  127: 'mult:n[1]', 128: 'knots:f[1]',
  130: 'node_id:d attributes_groups:p owner:p next:p previous:p geometric_owner:p sense:c type:c data:p tf:p internal_geom:p[1]',
  133: 'node_id:d attributes_groups:p owner:p next:p previous:p geometric_owner:p sense:c basis_curve:p point_1:v point_2:v parm_1:f parm_2:f',
  134: 'node_id:d attributes_groups:p owner:p next:p previous:p geometric_owner:p sense:c nurbs:p data:p',
  135: 'self_int:u analytic_form:p',
  136: 'degree:n n_vertices:d vertex_dim:n n_knots:d knot_type:u periodic:l closed:l rational:l curve_form:u bspline_vertices:p knot_mult:p knots:p',
  137: 'node_id:d attributes_groups:p owner:p next:p previous:p geometric_owner:p sense:c surface:p b_curve:p original:p tolerance_to_original:f',
  141: 'owner:p next:p previous:p shared_geometry:p',
};
const fieldsOf = (spec) => spec.split(' ').map((s) => { const m = /^([^:]*):(\w)(?:\[(\d+)\])?$/.exec(s); return { name: m[1], type: m[2], n: m[3] ? +m[3] : 0 }; });
const SURF_NAME = { 50: 'plane', 51: 'cylinder', 52: 'cone', 53: 'sphere', 54: 'torus', 56: 'blended-edge', 59: 'blend-bound', 60: 'offset-surface', 67: 'swept-surface', 68: 'spun-surface', 120: 'foreign-surface', 124: 'b-surface' };

/** 'text' | 'binary' (bare) | 'neutral' (neutral binary) | null, judged by the flag that follows the keyword header. */
export function xtKind(u8) {
  const head = String.fromCharCode(...u8.subarray(0, Math.min(u8.length, 4096)));
  if (!/\*\*PARASOLID|TRANSMIT FILE/.test(head) && !/^\*\*ABCDEFGH/.test(head) && !(u8[0] === 0x50 && u8[1] === 0x53 && u8[2] === 0)) return null;
  const e = head.indexOf('**END_OF_HEADER');
  let p = e >= 0 ? head.indexOf('\n', e) + 1 : 0;
  while (p < u8.length && (u8[p] === 10 || u8[p] === 13)) p++;
  if (u8[p] === 0x50 && u8[p + 1] === 0x53 && u8[p + 2] === 0 && u8[p + 3] === 0) return 'neutral';
  if (u8[p] === 0x54) return 'text';
  return u8[p] === 0x42 ? 'binary' : /FORMAT=binary/.test(head) ? 'neutral' : 'text';
}
export function parseXB() { return fail('This is a binary Parasolid file (.x_b): the binary transmit encoding is not read (the text form is). Export the model as Parasolid text (.x_t) or as STEP AP242 from the CAD package.'); }

/** Node table of a text transmit file: { nodes: Map index → { t, f: { name: value } }, version, schema, header }. */
export function xtNodes(text) {
  const eoh = text.indexOf('**END_OF_HEADER');
  if (eoh < 0 && !/^\s*T\d+ : TRANSMIT/.test(text)) fail('Not a Parasolid transmit file (the header is missing).');
  const header = {}, hd = eoh >= 0 ? text.slice(0, eoh).replace(/[\r\n]/g, '') : '';
  for (const m of hd.matchAll(/([A-Z_0-9]+)=([^;]*);/g)) header[m[1]] = m[2];
  // the data stream: line breaks are layout only
  const s = (eoh >= 0 ? text.slice(text.indexOf('\n', eoh) + 1) : text).replace(/[\r\n]/g, ''), n = s.length;
  let p = 0;
  const bad = (what) => fail(`The Parasolid file is truncated or not laid out as its schema says (${what}).`);
  const num = () => {
    if (p >= n) bad('unexpected end');
    if (s.charCodeAt(p) === 63) { p++; return null; }
    const e = s.indexOf(' ', p), t = s.slice(p, e < 0 ? n : e);
    p = e < 0 ? n : e + 1;
    const v = +t;
    if (t === '' || v !== v) bad(`number expected, found "${t.slice(0, 20)}"`);
    return v;
  };
  const chr = () => { if (p >= n) bad('unexpected end'); let c = s[p++]; if (c === '\\' && p < n) { const d = s[p++]; c = d === '9' ? '         ' : d === 'n' ? '\n' : d === 'r' ? '\r' : d === '0' ? '\0' : d; } return c; };
  const str = () => { const len = num(); if (!(len >= 0 && len <= 4096)) bad('string length'); let o = ''; while (o.length < len) o += chr(); return o; };
  if (s[p] !== 'T') fail(s.charCodeAt(p) === 0x50 || s[p] === 'B' ? 'This Parasolid file is binary. Export the model as Parasolid text (.x_t) or STEP AP242.' : 'Not a Parasolid text transmit file.');
  p++;
  const modeller = str(), schema = str(), vm = /version (\d+)/.exec(modeller), ver = vm ? +vm[1] : 0, embedded = /^SCH_\d+_\d+_\d+$/.test(schema), schemaNo = +(/^SCH_\d+_(\d+)/.exec(schema) || [0, 0])[1];
  if (/partition|delta/i.test(modeller)) fail('This is a Parasolid partition or delta file, not a part file; it is not read. Export the part as Parasolid text (.x_t) or STEP.');
  if (!embedded && schemaNo !== 13006) fail(`This Parasolid file (modeller version ${ver ? Math.floor(ver / 100000) : '?'}, schema ${schemaNo}) carries no embedded schema description and its node layout is not known to this reader (files with an embedded schema, the default of current CAD packages, are read). Re-save it from a current CAD package, or export STEP AP242.`);
  if (embedded) num();                                      // highest node type in use
  const usfld = num();
  if (!(usfld >= 0 && usfld <= 16)) bad('user-field size');
  const layouts = new Map(), nodes = new Map();
  const fieldDef = () => { const name = str(), cls = num(), ne = num(), type = cls ? 'p' : str(); if (ne === 1) chr(); return { name, type, n: ne }; };
  for (;;) {
    if (p >= n) bad('terminator missing');
    const t = num();
    if (t === 1) break;
    if (!(t > 1 && t < 4096)) bad(`node type ${t}`);
    let lay = layouts.get(t);
    if (!lay) {
      const base = BASE[t] ? fieldsOf(BASE[t]) : null;
      if (!embedded) lay = base || bad(`node type ${t} is not in the base schema`);
      else {
        const nf = num();
        if (nf === 255) lay = base || bad(`node type ${t} is not in the base schema`);
        else if (/[CDIAZ]/.test(s[p])) {
          if (!base) bad(`field edits for unknown node type ${t}`);
          lay = [];
          for (let b = 0, guard = 0; ; guard++) {
            const c = s[p++];
            if (c === 'Z' || guard > 4096) break;
            if (c === 'C') { if (b >= base.length) bad(`node type ${t} has more base fields than known`); lay.push(base[b++]); }
            else if (c === 'D') b++;
            else if (c === 'I' || c === 'A') lay.push(fieldDef());
            else bad(`schema edit "${c}"`);
          }
        } else { str(); str(); lay = []; for (let k = 0; k < nf; k++) lay.push(fieldDef()); }
      }
      layouts.set(t, lay);
    }
    const variable = lay.some((f) => f.n === 1), vlen = variable ? num() : 0, idx = num();
    if (!(vlen >= 0 && vlen <= n) || !(idx > 0) || nodes.size > 4e6) bad('node index');
    const f = {};
    for (const fd of lay) {
      const cnt = fd.n === 1 ? vlen : fd.n || 1, ty = fd.type;
      let v;
      if (ty === 'c' && fd.n) { v = ''; for (let k = 0; k < cnt; k++) v += chr(); }
      else {
        const one = () => { if (ty === 'c') return chr(); if (ty === 'l') return chr() === 'T'; if (ty === 'v' || ty === 'h') { if (s.charCodeAt(p) === 63) { p++; return null; } return [num(), num(), num()]; } if (ty === 'i') return [num(), num()]; if (ty === 'b') return [num(), num(), num(), num(), num(), num()]; return num(); };
        if (!fd.n) v = one(); else { v = new Array(cnt); for (let k = 0; k < cnt; k++) v[k] = one(); }
      }
      f[fd.name] = v;
    }
    for (let k = 0; k < usfld; k++) num();                 // user field: written after every node
    nodes.set(idx, { t, f });
  }
  return { nodes, version: ver ? `${Math.floor(ver / 100000)}.${Math.floor(ver / 1000) % 100}` : '?', schema, header };
}

/** Parse a text .x_t file into the neutral B-rep model (see fmt_brep.js); coordinates in metres. */
export function parseXT(text) {
  const { nodes, version, schema, header } = xtNodes(text), get = (i) => (i ? nodes.get(i) : undefined), of = (i, t) => { const nd = get(i); return nd && nd.t === t ? nd.f : null; };
  const warnings = [], skipped = {}, counts = { body: 0, shell: 0, face: 0, loop: 0, fin: 0, edge: 0, vertex: 0, surfaces: {}, curves: {} }, vec = (v) => (Array.isArray(v) && v.length === 3 && v.every(Number.isFinite) ? v : null);
  for (const nd of nodes.values()) { const k = { 12: 'body', 13: 'shell', 14: 'face', 15: 'loop', 17: 'fin', 16: 'edge', 18: 'vertex' }[nd.t]; if (k) counts[k]++; }
  const pointOf = (vi) => { const v = of(vi, 18), pt = v && of(v.point, 29); return pt ? vec(pt.pvec) : null; };
  const edgeCache = new Map();
  /** Edge sampled from its start to its end vertex. */
  const edgePts = (ei) => {
    let e = edgeCache.get(ei);
    if (e) return e;
    e = { pts: [], kind: 'other', closed: false };
    edgeCache.set(ei, e);
    const ed = of(ei, 16);
    if (!ed) return e;
    // a positive fin ends at the edge's end vertex, a negative one at its start vertex
    let A = null, B = null, ring = true, fcurve = 0;
    for (let fi = ed.fin, k = 0; fi && k < 64; k++) { const fn = of(fi, 17); if (!fn) break; if (fn.vertex) { ring = false; const q = pointOf(fn.vertex); if (fn.sense === '+') B = B || q; else A = A || q; } fcurve = fcurve || fn.curve; fi = fn.other; if (fi === ed.fin) break; }
    let cv = get(ed.curve || fcurve), fwd = true;
    for (let d = 0; cv && cv.t === 133 && d < 4; d++) cv = get(cv.f.basis_curve);
    if (cv) fwd = cv.f.sense !== '-';
    const ct = cv ? cv.t : 0;
    e.kind = ct === 30 ? 'line' : ct === 31 ? 'circle' : ct === 32 ? 'ellipse' : ct === 38 || ct === 134 || ct === 137 ? 'spline' : 'other';
    counts.curves[e.kind] = (counts.curves[e.kind] || 0) + 1;
    if (!ring && (!A || !B)) return e;
    if (ct === 31 || ct === 32) {
      const c = vec(cv.f.centre), nrm = vec(cv.f.normal), x = vec(cv.f.x_axis), a = ct === 31 ? cv.f.radius : cv.f.major_radius, b = ct === 31 ? cv.f.radius : cv.f.minor_radius;
      if (c && nrm && x && a > 0 && b > 0) {
        const M = x.map((q) => q * a), m = cross(nrm, x).map((q) => q * b), closed = ring || vlen(sub(A, B)) <= 1e-12 * a;
        e.pts = ellipseArc(c, M, m, ring ? null : A, closed ? null : B, fwd); e.closed = closed;
        if (closed) e.pts.pop();
        return e;
      }
    }
    if (ring) return e;
    if (ct === 38) {                                         // intersection curve: the chart is its chordal approximation
      const ch = of(cv.f.chart, 40), hv = ch && Array.isArray(ch.hvec) ? ch.hvec.filter(vec) : [];
      if (hv.length > 1) {
        const near = (q) => { let bi = 0, bd = Infinity; hv.forEach((h, k) => { const d = vlen(sub(h, q)); if (d < bd) { bd = d; bi = k; } }); return bi; };
        const ia = near(A), ib = near(B), mid = ia <= ib ? hv.slice(ia + 1, ib) : hv.slice(ib + 1, ia).reverse();
        e.pts = [A, ...mid.filter((h) => vlen(sub(h, A)) > 0 && vlen(sub(h, B)) > 0), B];
        return e;
      }
    }
    if (ct && ct !== 30) skipped[ct === 134 ? 'b-curve edge' : ct === 137 ? 'sp-curve edge' : ct === 38 ? 'intersection edge' : `curve type ${ct}`] = (skipped[ct === 134 ? 'b-curve edge' : ct === 137 ? 'sp-curve edge' : ct === 38 ? 'intersection edge' : `curve type ${ct}`] || 0) + 1;
    e.pts = vlen(sub(A, B)) === 0 ? [A] : [A, B];
    return e;
  };
  const surfOf = (si) => {
    const nd = get(si);
    if (!nd) return { type: 'other', name: 'no surface' };
    const f = nd.f, plus = f.sense !== '-', z = vec(f.axis || f.normal), o = vec(f.pvec || f.centre), name = SURF_NAME[nd.t] || `surface type ${nd.t}`;
    if (nd.t < 50 || nd.t > 54 || !z || !o || !vlen(z)) return { type: nd.t === 124 ? 'spline' : 'other', name, plus };
    const base = { o, z: unit(z), x: perp(unit(z), vec(f.x_axis)), plus, name };
    if (nd.t === 50) return { type: 'plane', ...base };
    if (nd.t === 51) return f.radius > 0 ? { type: 'cylinder', ...base, r: f.radius, tanA: 0 } : { type: 'other', name, plus };
    // cone: the radius grows against the axis direction and the natural normal points towards the axis
    if (nd.t === 52) return f.radius >= 0 && f.cos_half_angle ? { type: 'cone', ...base, r: f.radius, tanA: f.sin_half_angle / f.cos_half_angle, inward: true } : { type: 'other', name, plus };
    if (nd.t === 53) return f.radius > 0 ? { type: 'sphere', ...base, r: f.radius } : { type: 'other', name, plus };
    return f.major_radius > 0 && f.minor_radius > 0 && f.major_radius > f.minor_radius ? { type: 'torus', ...base, R: f.major_radius, r: f.minor_radius } : { type: 'other', name: 'self-intersecting torus', plus };
  };
  const faces = [];
  let guard = 0;
  for (const [fi, nd] of nodes) {
    if (nd.t !== 14) continue;
    const { plus, name, inward, ...surf } = surfOf(nd.f.surface), loops = [];
    counts.surfaces[surf.type] = (counts.surfaces[surf.type] || 0) + 1;
    if (surf.type === 'spline' || surf.type === 'other') skipped[name] = (skipped[name] || 0) + 1;
    for (let li = nd.f.loop, nl = 0; li && nl < 1e5; nl++) {
      const lp = of(li, 15);
      if (!lp) break;
      const pts = [];
      let single = false;
      for (let ni = lp.fin, nf = 0; ni && nf < 1e5; nf++) {
        const fn = of(ni, 17);
        if (!fn || ++guard > 2e7) break;
        if (!fn.edge) { const q = pointOf(fn.vertex); if (q) { pts.push(q); single = true; } break; }
        const e = edgePts(fn.edge), back = fn.sense === '-', seq = back && !e.closed ? e.pts.slice().reverse() : back ? [e.pts[0], ...e.pts.slice(1).reverse()] : e.pts;
        for (let k = 0; k < seq.length - (e.closed ? 0 : 1) || (k === 0 && seq.length === 1); k++) pts.push(seq[k]);
        ni = fn.forward;
        if (ni === lp.fin) break;
      }
      if (pts.length) loops.push({ pts, single: single && pts.length === 1 });
      li = lp.next;
    }
    faces.push({ surf, sense: ((nd.f.sense !== '-') === (plus !== false)) !== !!inward, loops, id: fi });
  }
  const edges = [];
  for (const [ei, nd] of nodes) if (nd.t === 16) { const e = edgePts(ei); if (e.pts.length > 1) edges.push({ pts: e.pts, kind: e.kind, closed: e.closed }); }
  let inst = 0;
  for (const nd of nodes.values()) if (nd.t === 11) inst++;
  if (inst) warnings.push(`The file is an assembly with ${inst} instances; their placements are not applied and each body is shown once in its own coordinates.`);
  if (!faces.length && !edges.length) fail('The Parasolid file holds no faces or edges.');
  return { format: 'Parasolid X_T', version, schema, application: header.APPL || '', unitScale: 1, faces, edges, vertices: counts.vertex, counts, skipped, warnings };
}
