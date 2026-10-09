// Parasolid transmit-file reader (text .x_t and neutral binary .x_b) producing the neutral B-rep model of fmt_brep.js,
// written against the published "Parasolid XT Format Reference".
//
// Read:   the keyword header; the text stream ("T") and the neutral binary stream ("PS", big-endian, as written to .x_b
//         files and to the XT B-rep segments of JT files; checked on streams of Parasolid 17 and 25) of modeller versions
//         14 and later, whose node records are laid out by the base schema SCH_13006 plus the field edits embedded in the
//         file, so node types and fields added by newer kernels are followed or skipped exactly; files without an
//         embedded schema of the schemas 9008, 10004, 11004, 12006 / 13006 (layouts fitted to real files) and of fixed
//         later schemas (per node type the newest layout met in schema-embedding files: checked for 16100 and 32001);
//         the topology body → region → shell → face → loop → fin → edge → vertex → point; edge curves LINE, CIRCLE,
//         ELLIPSE, INTERSECTION (its chart of points), B_CURVE (sampled; rational forms included) and TRIMMED_CURVE;
//         surfaces PLANE, CYLINDER, CONE, SPHERE and TORUS, SWEPT_SURF of a line (a plane), of a circle along its normal
//         (a cylinder) and of a B-curve (a ruled B-spline strip), and B_SURFACE as its control net; face and surface
//         senses; assemblies: every instance places its body by its transform (rotation, translation, scale; nested).
//         Parasolid models are in metres.
// Not read: OFFSET_SURF, SPUN_SURF, BLENDED_EDGE and foreign surfaces (their faces are counted in `skipped`, the boundary
//         is still delivered), SP_CURVE edges (chord between the vertices), periodic B-spline knot sets, attributes,
//         partition and delta files, the machine-dependent "bare" binary form, schemas before 13006 other than those
//         named. No sample with a B_SURFACE face was at hand: that path follows the reference but is unverified.
// File content is untrusted: counts, lengths and indices are bounds-checked and the node total is capped.
import { sub, cross, dot, vlen, unit, perp, ellipseArc, nurbsCurve, nurbsSurfValid } from './fmt_brep.js';

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
  99: 'names:p[1]',
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
// Node layouts of schemas newer than the base, as met in files that embed their schema: type -> [[first schema seen, fields] …].
// A file written against a fixed schema without embedding it takes, per node type, the newest layout not above its schema.
const LATER = {
  10: [[28101, 'highest_node_id:d attributes_groups:p attribute_chains:p list:p surface:p curve:p point:p mesh:p polyline:p key:p res_size:f res_linear:f ref_instance:p next:p previous:p state:u owner:p type:u sub_instance:p'], [31001, 'highest_node_id:d attributes_groups:p attribute_chains:p list:p surface:p curve:p point:p mesh:p polyline:p key:p res_size:f res_linear:f ref_instance:p next:p previous:p state:u owner:p type:u sub_instance:p mesh_offset_data:p'], [34101, 'highest_node_id:d attributes_groups:p attribute_chains:p list:p lattice:p surface:p curve:p point:p mesh:p polyline:p key:p res_size:f res_linear:f ref_instance:p next:p previous:p state:u owner:p type:u sub_instance:p mesh_offset_data:p']],
  12: [[19008, 'highest_node_id:d attributes_groups:p attribute_chains:p surface:p curve:p point:p key:p res_size:f res_linear:f ref_instance:p next:p previous:p state:u owner:p body_type:u nom_geom_state:u shell:p boundary_surface:p boundary_curve:p boundary_point:p region:p edge:p vertex:p index_map_offset:d index_map:p node_id_index_map:p schema_embedding_map:p'], [25001, 'highest_node_id:d attributes_groups:p attribute_chains:p surface:p curve:p point:p key:p res_size:f res_linear:f ref_instance:p next:p previous:p state:u owner:p body_type:u nom_geom_state:u shell:p boundary_surface:p boundary_curve:p boundary_point:p region:p edge:p vertex:p index_map_offset:d index_map:p node_id_index_map:p schema_embedding_map:p child:p lowest_node_id:d'], [26105, 'highest_node_id:d attributes_groups:p attribute_chains:p surface:p curve:p point:p key:p res_size:f res_linear:f ref_instance:p next:p previous:p state:u owner:p body_type:u nom_geom_state:u shell:p boundary_surface:p boundary_curve:p boundary_point:p region:p edge:p vertex:p boundary_mesh:p index_map_offset:d index_map:p node_id_index_map:p schema_embedding_map:p child:p lowest_node_id:d'], [28101, 'highest_node_id:d attributes_groups:p attribute_chains:p surface:p curve:p point:p mesh:p polyline:p key:p res_size:f res_linear:f ref_instance:p next:p previous:p state:u owner:p body_type:u nom_geom_state:u shell:p boundary_surface:p boundary_curve:p boundary_point:p boundary_mesh:p boundary_polyline:p region:p edge:p vertex:p index_map_offset:d index_map:p node_id_index_map:p schema_embedding_map:p child:p lowest_node_id:d'], [31001, 'highest_node_id:d attributes_groups:p attribute_chains:p surface:p curve:p point:p mesh:p polyline:p key:p res_size:f res_linear:f ref_instance:p next:p previous:p state:u owner:p body_type:u nom_geom_state:u shell:p boundary_surface:p boundary_curve:p boundary_point:p boundary_mesh:p boundary_polyline:p region:p edge:p vertex:p index_map_offset:d index_map:p node_id_index_map:p schema_embedding_map:p child:p lowest_node_id:d mesh_offset_data:p'], [33103, 'highest_node_id:d attributes_groups:p attribute_chains:p lattice:p surface:p curve:p point:p mesh:p polyline:p key:p res_size:f res_linear:f ref_instance:p next:p previous:p state:u owner:p body_type:u nom_geom_state:u shell:p boundary_surface:p boundary_curve:p boundary_point:p boundary_mesh:p boundary_polyline:p region:p edge:p vertex:p index_map_offset:d index_map:p node_id_index_map:p schema_embedding_map:p child:p lowest_node_id:d mesh_offset_data:p'], [34101, 'highest_node_id:d attributes_groups:p attribute_chains:p lattice:p surface:p curve:p point:p mesh:p polyline:p key:p res_size:f res_linear:f ref_instance:p next:p previous:p state:u owner:p body_type:u nom_geom_state:u shell:p boundary_lattice:p boundary_surface:p boundary_curve:p boundary_point:p boundary_mesh:p boundary_polyline:p region:p edge:p vertex:p index_map_offset:d index_map:p node_id_index_map:p schema_embedding_map:p child:p lowest_node_id:d mesh_offset_data:p']],
  19: [[25001, 'node_id:d attributes_groups:p body:p next:p previous:p shell:p type:c owner:p'], [34101, 'node_id:d attributes_groups:p body:p next:p previous:p shell:p frame:p type:c owner:p']],
  38: [[30000, 'node_id:d attributes_groups:p owner:p next:p previous:p geometric_owner:p sense:c surface:p[2] chart:p start:p end:p intersection_data:p']],
  41: [[28101, 'type:c term_use:c hvec:h[1]']],
  70: [[16100, 'node_id:d list_type:u notransmit:l owner:p next:p previous:p list_length:d block_length:d finger_index:d finger_block:p list_block:p']],
  74: [[19008, 'n_entries:d index_map_offset:d next_block:p entries:p[1]']],
  80: [[34101, 'next:p identifier:p type_id:d actions:u[8] field_names:p legal_owners:l[16] fields:u[1]']],
  100: [[34101, 'node_id:d owner:p next:p previous:p rotation_matrix:f[9] translation_vector:v scale:f flag:d perspective_vector:v precision:p']],
  204: [[30000, 'uv_type:u values:f[1]']],
  // met only in files without an embedded schema: laid out from their contents (text form)
  176: [[30000, 'n_entries:d a:p b:p c:p d:p entries:p[1]']],
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
/** Value readers over the text form: every number is written out and followed by one space; "?" stands for a null value. */
function textStream(text) {
  const eoh = text.indexOf('**END_OF_HEADER');
  if (eoh < 0 && !/^\s*T\d+ : TRANSMIT/.test(text)) fail('Not a Parasolid transmit file (the header is missing).');
  const header = {}, hd = eoh >= 0 ? text.slice(0, eoh).replace(/[\r\n]/g, '') : '';
  for (const m of hd.matchAll(/([A-Z_0-9]+)=([^;]*);/g)) header[m[1]] = m[2];
  // the data stream: line breaks are layout only
  const s = (eoh >= 0 ? text.slice(text.indexOf('\n', eoh) + 1) : text).replace(/[\r\n]/g, ''), n = s.length;
  let p = 0;
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
  if (s[p] !== 'T') fail(s.charCodeAt(p) === 0x50 || s[p] === 'B' ? 'This Parasolid file is binary; read it as .x_b.' : 'Not a Parasolid text transmit file.');
  p++;
  return { header, size: n, head: [str, str], eof: () => p >= n, raw: () => s[p++], type: num, len: num, idx: num, byte: num, short: num, pos: num, sstr: str, chr, d: num, f: num, p: num, vec: () => { if (s.charCodeAt(p) === 63) { p++; return null; } return [num(), num(), num()]; } };
}
/** Value readers over the neutral binary form: big-endian, IEEE doubles, pointer indices in one or two short integers. */
function binaryStream(u8) {
  const N = u8.length, dv = new DataView(u8.buffer, u8.byteOffset, N), header = {};
  let p = 0;
  // an optional text header precedes the "PS\0\0" flag
  const eoh = String.fromCharCode(...u8.subarray(0, Math.min(N, 8192))).indexOf('**END_OF_HEADER');
  if (eoh >= 0) { const hd = String.fromCharCode(...u8.subarray(0, eoh)).replace(/[\r\n]/g, ''); for (const m of hd.matchAll(/([A-Z_0-9]+)=([^;]*);/g)) header[m[1]] = m[2]; p = eoh; while (p < N && u8[p] !== 10) p++; while (p < N && (u8[p] === 10 || u8[p] === 13)) p++; }
  if (u8[p] === 0x42) fail('This Parasolid file is in the machine-dependent "bare" binary form, which is not read (neutral binary and text are). Export the model as Parasolid text (.x_t) or neutral binary.');
  if (!(u8[p] === 0x50 && u8[p + 1] === 0x53 && u8[p + 2] === 0 && u8[p + 3] === 0)) fail('Not a Parasolid neutral binary file (the "PS" flag is missing).');
  p += 4;
  const need = (k) => { if (p + k > N) bad('unexpected end'); }, byte = () => { need(1); return u8[p++]; }, short = () => { need(2); p += 2; return dv.getInt16(p - 2); };
  const int = () => { need(4); p += 4; const v = dv.getInt32(p - 4); return v === -32764 ? null : v; }, dbl = () => { need(8); p += 8; const v = dv.getFloat64(p - 8); return v === -3.14158e13 ? null : v; };
  const ptr = () => { let r = short(), q = 0; if (r < 0) { q = short(); r = -r; } return q * 32767 + r - 1; };
  const chars = (k) => { if (!(k >= 0 && k <= 4096)) bad('string length'); need(k); let o = ''; for (let i = 0; i < k; i++) o += String.fromCharCode(u8[p + i]); p += k; return o; };
  return { header, size: N, binary: true, head: [() => chars(short()), () => chars(int())], eof: () => p >= N, raw: () => String.fromCharCode(byte()), type: short, len: int, idx: ptr, byte, short, pos: ptr, sstr: () => chars(byte()), chr: () => String.fromCharCode(byte()), d: int, f: dbl, p: ptr, n: () => { const v = short(); return v === -32764 ? null : v; }, u: byte, int,
    vec: () => { const v = [dbl(), dbl(), dbl()]; return v[0] === null ? null : v; } };
}
const bad = (what) => fail(`The Parasolid file is truncated or not laid out as its schema says (${what}).`);

/** Node table of a transmit stream: { nodes: Map index → { t, f: { name: value } }, version, schema, header }. */
function nodeTable(R) {
  const modeller = R.head[0](), schema = R.head[1](), vm = /version (\d+)/.exec(modeller), ver = vm ? +vm[1] : 0, embedded = /^SCH_\d+_\d+_\d+$/.test(schema), schemaNo = +(/^SCH_\d+_(\d+)/.exec(schema) || [0, 0])[1], n = R.size;
  if (/partition|delta/i.test(modeller)) fail('This is a Parasolid partition or delta file, not a part file; it is not read. Export the part as Parasolid text (.x_t) or STEP.');
  // schemas before the base are known for the versions listed in SCHEMA_EDITS only; later ones follow the layouts collected in LATER
  const edits = embedded ? null : schemaNo > 13006 ? Object.fromEntries(Object.entries(LATER).map(([t, l]) => [t, (l.filter((e) => e[0] <= schemaNo).pop() || [])[1]])) : SCHEMA_EDITS[schemaNo];
  if (!embedded && schemaNo !== 13006 && !edits) fail(`This Parasolid file (modeller version ${ver ? Math.floor(ver / 100000) : '?'}, schema ${schemaNo}) carries no embedded schema description and its node layout is not known to this reader (files with an embedded schema, the default since Parasolid 14, and the fixed schemas ${Object.keys(SCHEMA_EDITS).join(', ')}, 13006 and later are read). Export the model again from a current CAD package, or as STEP AP214 / AP242.`);
  if (embedded) R.short();                                  // highest node type in use
  const usfld = R.binary ? R.int() : R.d();
  if (!(usfld >= 0 && usfld <= 16)) bad('user-field size');
  const layouts = new Map(), nodes = new Map();
  const fieldDef = () => { const name = R.sstr(), cls = R.short(), ne = R.pos(), type = cls ? 'p' : R.sstr(); if (ne === 1) R.raw(); return { name, type, n: ne }; };
  for (;;) {
    if (R.eof()) bad('terminator missing');
    const t = R.type();
    if (t === 1) break;
    if (!(t > 1 && t < 4096)) bad(`node type ${t}`);
    let lay = layouts.get(t);
    if (!lay) {
      const base = BASE[t] ? fieldsOf(BASE[t]) : null;
      if (!embedded) { lay = edits && edits[t] ? fieldsOf(edits[t]) : base; if (!lay) bad(`node type ${t} is not in the schema ${schemaNo} layouts of this reader`); }
      else {
        const nf = R.byte();
        if (nf === 255) lay = base || bad(`node type ${t} is not in the base schema`);
        else if (base) {
          lay = [];
          for (let b = 0, guard = 0; ; guard++) {
            const c = R.raw();
            if (c === 'Z' || guard > 4096) break;
            if (c === 'C') { if (b >= base.length) bad(`node type ${t} has more base fields than known`); lay.push(base[b++]); }
            else if (c === 'D') b++;
            else if (c === 'I' || c === 'A') lay.push(fieldDef());
            else bad(`schema edit "${c}"`);
          }
        } else { R.sstr(); R.sstr(); lay = []; for (let k = 0; k < nf; k++) lay.push(fieldDef()); }
      }
      layouts.set(t, lay);
    }
    const variable = lay.some((f) => f.n === 1), vlen = variable ? R.len() : 0, idx = R.idx();
    if (!(vlen >= 0 && vlen <= n) || !(idx > 0) || nodes.size > 4e6) bad('node index');
    const f = {};
    for (const fd of lay) {
      const cnt = fd.n === 1 ? vlen : fd.n || 1, ty = fd.type;
      let v;
      if (ty === 'c' && fd.n) { v = ''; for (let k = 0; k < cnt; k++) v += R.chr(); }
      else {
        const one = () => {
          if (ty === 'c') return R.chr();
          if (ty === 'l') return R.binary ? R.byte() !== 0 : R.chr() === 'T';
          if (ty === 'v' || ty === 'h') return R.vec();
          if (ty === 'i') return [R.f(), R.f()];
          if (ty === 'b') return [R.f(), R.f(), R.f(), R.f(), R.f(), R.f()];
          if (!R.binary) return R.d();
          return ty === 'f' ? R.f() : ty === 'p' ? R.p() : ty === 'u' ? R.u() : ty === 'n' || ty === 'w' ? R.n() : ty === 'd' || ty === 't' ? R.d() : bad(`field type "${ty}"`);
        };
        if (!fd.n) v = one(); else { v = new Array(cnt); for (let k = 0; k < cnt; k++) v[k] = one(); }
      }
      f[fd.name] = v;
    }
    for (let k = 0; k < usfld; k++) R.d();                 // user field: written after every node
    nodes.set(idx, { t, f });
  }
  return { nodes, version: ver ? `${Math.floor(ver / 100000)}.${Math.floor(ver / 1000) % 100}` : '?', schema, header: R.header };
}
export const xtNodes = (text) => nodeTable(textStream(text));
// Layouts of schemas written without an embedded description, as the node types whose fields differ from SCH_13006.
const OLD_BODY = BASE[12].replace(' nom_geom_state:u', ''), OLD_FIN = BASE[17].replace('attributes_groups:p ', '');
const OLD_ATTDEF = 'next:p identifier:p type_id:d actions:u[8] legal_owners:l[13] fields:u[1]';
const SCHEMA_EDITS = { 9008: { 12: OLD_BODY, 17: OLD_FIN, 80: OLD_ATTDEF }, 10004: { 17: OLD_FIN, 80: OLD_ATTDEF }, 11004: {}, 12006: {} };

/** Parse a neutral binary transmit file (.x_b) into the neutral B-rep model; coordinates in metres. */
export function parseXB(u8) { return model(nodeTable(binaryStream(u8)), 'Parasolid X_B'); }
/** Parse a text .x_t file into the neutral B-rep model (see fmt_brep.js); coordinates in metres. */
export function parseXT(text) { return model(nodeTable(textStream(text)), 'Parasolid X_T'); }

function model({ nodes, version, schema, header }, format) {
  const get = (i) => (i ? nodes.get(i) : undefined), of = (i, t) => { const nd = get(i); return nd && nd.t === t ? nd.f : null; };
  const warnings = [], skipped = {}, counts = { body: 0, shell: 0, face: 0, loop: 0, fin: 0, edge: 0, vertex: 0, surfaces: {}, curves: {} }, vec = (v) => (Array.isArray(v) && v.length === 3 && v.every(Number.isFinite) ? v : null);
  for (const nd of nodes.values()) { const k = { 12: 'body', 13: 'shell', 14: 'face', 15: 'loop', 17: 'fin', 16: 'edge', 18: 'vertex' }[nd.t]; if (k) counts[k]++; }
  const pointOf = (vi) => { const v = of(vi, 18), pt = v && of(v.point, 29); return pt ? vec(pt.pvec) : null; };
  const edgeCache = new Map();
  const arr = (i, t, key) => { const nd = of(i, t), v = nd && nd[key]; return Array.isArray(v) ? v : null; };
  /** Full knot vector of a knot set with its multiplicities, or null. */
  const knotsOf = (ki, mi, nk) => { const k = arr(ki, 128, 'knots'), m = arr(mi, 127, 'mult'), out = []; if (!k || !m || k.length !== m.length || k.length !== nk) return null; for (let i = 0; i < k.length; i++) { if (!(m[i] >= 1 && m[i] <= 26) || out.length > 4e5) return null; for (let q = 0; q < m[i]; q++) out.push(k[i]); } return out; };
  /** Control vertices: vertex_dim values each; a rational form carries (w·x, w·y, w·z, w). */
  const polesOf = (vi, n, dim, rational) => {
    const v = arr(vi, 45, 'vertices'), cps = [], w = rational ? [] : null;
    if (!v || !(dim === 3 || (rational && dim === 4)) || v.length !== n * dim || n > 4e5) return null;
    for (let i = 0; i < n; i++) { const h = rational ? v[i * dim + 3] : 1; if (!(h > 0)) return null; cps.push([v[i * dim] / h, v[i * dim + 1] / h, v[i * dim + 2] / h]); if (w) w.push(h); }
    return { cps, w };
  };
  /** NURBS_CURVE node -> { deg, knots, cps, w } or null (periodic knot sets are not unwrapped). */
  const curveOf = (ni) => {
    const c = of(ni, 136);
    if (!c) return null;
    const knots = knotsOf(c.knots, c.knot_mult, c.n_knots), pl = polesOf(c.bspline_vertices, c.n_vertices, c.vertex_dim, !!c.rational);
    return knots && pl && knots.length === pl.cps.length + c.degree + 1 ? { deg: c.degree, knots, ...pl } : null;
  };
  /** NURBS_SURF node -> surface for nurbsSurfPoint or null; the vertices are stored with v running fastest. */
  const bsurfOf = (ni) => {
    const c = of(ni, 126);
    if (!c) return null;
    const ku = knotsOf(c.u_knots, c.u_knot_mult, c.n_u_knots), kv = knotsOf(c.v_knots, c.v_knot_mult, c.n_v_knots), nU = c.n_u_vertices, nV = c.n_v_vertices, pl = nU > 0 && nV > 0 ? polesOf(c.bspline_vertices, nU * nV, c.vertex_dim, !!c.rational) : null;
    if (!ku || !kv || !pl) return null;
    const g = { degU: c.u_degree, degV: c.v_degree, nU, nV, knotsU: ku, knotsV: kv, cps: [], w: pl.w ? [] : null };
    for (let b = 0; b < nV; b++) for (let a = 0; a < nU; a++) { g.cps.push(pl.cps[a * nV + b]); if (pl.w) g.w.push(pl.w[a * nV + b]); }
    return nurbsSurfValid(g) ? g : null;
  };
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
    if (ct === 134) {
      const c = curveOf(cv.f.nurbs), smp = c ? nurbsCurve(c.deg, c.knots, c.cps, c.w, NaN, NaN, Math.min(1000, 16 * c.cps.length)) : null;
      if (smp) {
        const near = (q) => { let bi = 0, bd = Infinity; smp.forEach((h, k) => { const d = vlen(sub(h, q)); if (d < bd) { bd = d; bi = k; } }); return bi; }, ia = near(A), ib = near(B);
        if (ia !== ib) { e.pts = [A, ...(ia < ib ? smp.slice(ia + 1, ib) : smp.slice(ib + 1, ia).reverse()), B]; return e; }
        if (vlen(sub(A, B)) === 0 && vlen(sub(smp[0], smp[smp.length - 1])) <= 1e-9 * (vlen(sub(smp[0], smp[smp.length >> 1])) || 1)) { e.pts = [A, ...smp.slice(ia + 1, -1), ...smp.slice(0, ia)]; e.closed = true; return e; }
      }
    }
    if (ct && ct !== 30) skipped[ct === 134 ? 'b-curve edge' : ct === 137 ? 'sp-curve edge' : ct === 38 ? 'intersection edge' : `curve type ${ct}`] = (skipped[ct === 134 ? 'b-curve edge' : ct === 137 ? 'sp-curve edge' : ct === 38 ? 'intersection edge' : `curve type ${ct}`] || 0) + 1;
    e.pts = vlen(sub(A, B)) === 0 ? [A] : [A, B];
    return e;
  };
  let size = 0;
  /** Largest coordinate magnitude of the model's points: bounds the sweep of an unbounded swept surface. */
  const extent = () => { if (!size) for (const nd of nodes.values()) if (nd.t === 29 && vec(nd.f.pvec)) size = Math.max(size, Math.abs(nd.f.pvec[0]), Math.abs(nd.f.pvec[1]), Math.abs(nd.f.pvec[2])); return size; };
  const surfOf = (si) => {
    const nd = get(si);
    if (!nd) return { type: 'other', name: 'no surface' };
    const f = nd.f, plus = f.sense !== '-', z = vec(f.axis || f.normal), o = vec(f.pvec || f.centre), name = SURF_NAME[nd.t] || `surface type ${nd.t}`;
    if (nd.t === 124) { const g = bsurfOf(f.nurbs); return g ? { type: 'spline', nurbs: g, name, plus } : { type: 'spline', name, plus }; }
    if (nd.t === 67) {
      // swept surface C(u) + v·D: a line gives a plane, a circle swept along its normal a cylinder, a B-curve a ruled B-spline strip
      const D = vec(f.sweep), sec = get(f.section), sf = sec ? sec.f : null;
      if (D && vlen(D) && sec && sec.t === 30 && vec(sf.pvec) && vec(sf.direction) && vlen(cross(sf.direction, D)) > 1e-9) { const n = unit(cross(sf.direction, D)); return { type: 'plane', o: sf.pvec, z: n, x: unit(sf.direction), plus, name }; }
      if (D && sec && sec.t === 31 && vec(sf.centre) && vec(sf.normal) && sf.radius > 0 && vlen(cross(sf.normal, D)) < 1e-9 * vlen(D)) { const n = unit(sf.normal); return { type: 'cylinder', o: sf.centre, z: n, x: perp(n, vec(sf.x_axis)), r: sf.radius, tanA: 0, plus: plus === (sf.sense !== '-'), name }; }
      const c = D && sec && sec.t === 134 ? curveOf(sf.nurbs) : null;
      if (c) {
        const d = unit(D), far = Math.max(1, 2 * extent()), lo = c.cps.map((q) => [q[0] - far * d[0], q[1] - far * d[1], q[2] - far * d[2]]), hi = c.cps.map((q) => [q[0] + far * d[0], q[1] + far * d[1], q[2] + far * d[2]]);
        const g = { degU: c.deg, degV: 1, nU: c.cps.length, nV: 2, knotsU: c.knots, knotsV: [-far, -far, far, far], cps: [...lo, ...hi], w: c.w ? [...c.w, ...c.w] : null };
        if (nurbsSurfValid(g)) return { type: 'spline', nurbs: g, plus: plus === (sf.sense !== '-'), name };
      }
    }
    if (nd.t < 50 || nd.t > 54 || !z || !o || !vlen(z)) return { type: 'other', name, plus };
    const base = { o, z: unit(z), x: perp(unit(z), vec(f.x_axis)), plus, name };
    if (nd.t === 50) return { type: 'plane', ...base };
    if (nd.t === 51) return f.radius > 0 ? { type: 'cylinder', ...base, r: f.radius, tanA: 0 } : { type: 'other', name, plus };
    // cone: the radius grows along the axis by tan(half angle) (checked on real files: boundary points lie on this surface)
    if (nd.t === 52) return f.radius >= 0 && f.cos_half_angle ? { type: 'cone', ...base, r: f.radius, tanA: f.sin_half_angle / f.cos_half_angle } : { type: 'other', name, plus };
    if (nd.t === 53) return f.radius > 0 ? { type: 'sphere', ...base, r: f.radius } : { type: 'other', name, plus };
    return f.major_radius > 0 && f.minor_radius > 0 && f.major_radius > f.minor_radius ? { type: 'torus', ...base, R: f.major_radius, r: f.minor_radius } : { type: 'other', name: 'self-intersecting torus', plus };
  };
  const faces = [];
  let guard = 0;
  for (const [fi, nd] of nodes) {
    if (nd.t !== 14) continue;
    const { plus, name, ...surf } = surfOf(nd.f.surface), loops = [];
    counts.surfaces[surf.type] = (counts.surfaces[surf.type] || 0) + 1;
    if ((surf.type === 'spline' && !surf.nurbs) || surf.type === 'other') skipped[name] = (skipped[name] || 0) + 1;
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
    const sh = of(nd.f.shell, 13);
    faces.push({ surf, sense: (nd.f.sense !== '-') === (plus !== false), loops, id: fi, body: sh ? sh.body : 0 });
  }
  const edges = [];
  const bodyOfEdge = (ed) => { const fn = of(ed.fin, 17), lp = fn && of(fn.loop, 15), fc = lp && of(lp.face, 14), sh = fc && of(fc.shell, 13); return sh ? sh.body : get(ed.owner) && get(ed.owner).t === 12 ? ed.owner : 0; };
  for (const [ei, nd] of nodes) if (nd.t === 16) { const e = edgePts(ei); if (e.pts.length > 1) edges.push({ pts: e.pts, kind: e.kind, closed: e.closed, body: bodyOfEdge(nd.f) }); }
  // assemblies: every instance places its part (a body or a sub-assembly) by a transform x' = scale · R·x + t
  const places = new Map();                                 // body index -> list of transforms (null = identity)
  let inst = 0, nPlaced = 0;
  const mul = (A, B) => (!A ? B : !B ? A : { R: [0, 1, 2].flatMap((i) => [0, 1, 2].map((j) => A.R[3 * i] * B.R[j] + A.R[3 * i + 1] * B.R[3 + j] + A.R[3 * i + 2] * B.R[6 + j])), t: [0, 1, 2].map((i) => A.k * (A.R[3 * i] * B.t[0] + A.R[3 * i + 1] * B.t[1] + A.R[3 * i + 2] * B.t[2]) + A.t[i]), k: A.k * B.k });
  const walk = (ai, T, depth, seen) => {
    const as = of(ai, 10);
    if (!as || depth > 16 || seen.has(ai)) return;
    const inner = new Set(seen).add(ai);
    for (let ii = as.sub_instance, g = 0; ii && g < 1e5; g++) {
      const ins = of(ii, 11);
      if (!ins || ++inst > 1e5) break;
      const tf = of(ins.transform, 100), R = tf && Array.isArray(tf.rotation_matrix) && tf.rotation_matrix.length === 9 && tf.rotation_matrix.every(Number.isFinite) ? tf.rotation_matrix : null;
      const own = R ? { R, t: vec(tf.translation_vector) || [0, 0, 0], k: tf.scale > 0 ? tf.scale : 1 } : null, M = mul(T, own), part = get(ins.part);
      if (part && part.t === 12 && nPlaced++ < 1e5) { const l = places.get(ins.part); if (l) l.push(M); else places.set(ins.part, [M]); }
      else if (part && part.t === 10) walk(ins.part, M, depth + 1, inner);
      ii = ins.next_in_part;
    }
  };
  const root = get(1);
  if (root && root.t === 10) walk(1, null, 0, new Set());
  if (places.size) {
    const apply = (M, q) => [0, 1, 2].map((i) => M.k * (M.R[3 * i] * q[0] + M.R[3 * i + 1] * q[1] + M.R[3 * i + 2] * q[2]) + M.t[i]), dir = (M, q) => [0, 1, 2].map((i) => M.R[3 * i] * q[0] + M.R[3 * i + 1] * q[1] + M.R[3 * i + 2] * q[2]);
    const det = (R) => R[0] * (R[4] * R[8] - R[5] * R[7]) - R[1] * (R[3] * R[8] - R[5] * R[6]) + R[2] * (R[3] * R[7] - R[4] * R[6]);
    const outF = [], outE = [];
    let total = 0;
    const put = (list, out, map) => { for (const it of list) { const Ms = places.get(it.body) || [null]; for (const M of Ms) { if (++total > 2e6) fail('The Parasolid assembly expands to too many faces and edges.'); out.push(M ? map(it, M) : it); } } };
    put(faces, outF, (f, M) => {
      const sf = f.surf, flip = det(M.R) < 0, t = { ...sf };
      if (sf.o) { t.o = apply(M, sf.o); t.z = unit(dir(M, sf.z)); t.x = unit(dir(M, sf.x)); }
      if (sf.r !== undefined) t.r = sf.r * M.k;
      if (sf.R !== undefined) t.R = sf.R * M.k;
      if (sf.nurbs) t.nurbs = { ...sf.nurbs, cps: sf.nurbs.cps.map((q) => apply(M, q)) };
      return { ...f, surf: t, sense: f.sense !== flip, loops: f.loops.map((l) => ({ ...l, pts: l.pts.map((q) => apply(M, q)) })) };
    });
    put(edges, outE, (e, M) => ({ ...e, pts: e.pts.map((q) => apply(M, q)) }));
    faces.length = 0; edges.length = 0;
    for (const f of outF) faces.push(f);
    for (const e of outE) edges.push(e);
    warnings.push(`The file is an assembly: ${nPlaced} placements of ${places.size} bodies were applied.`);
  }
  if (!faces.length && !edges.length) fail('The Parasolid file holds no faces or edges.');
  return { format, version, schema, application: header.APPL || '', unitScale: 1, faces, edges, vertices: counts.vertex, counts, skipped, warnings };
}
