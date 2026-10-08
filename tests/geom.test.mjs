// Geometry-import engine tests: node tests/geom.test.mjs
// Every reader gets a small file synthesised in memory; the analysis helpers and generators are checked on known shapes.
import zlib from 'node:zlib';
import { FORMATS, PATHWAYS, SUITE_GEOMETRY, formatOf, importGeometry, sectionOf, maskOf, gridOf, microstructure, networkSummary, dimensions, generate } from '../js/core/geom.js';

let pass = 0, fails = 0;
const t0 = Date.now();
const ok = (cond, msg) => { if (cond) { pass++; console.log('ok   ' + msg); } else { fails++; console.log('FAIL ' + msg); } };
const near = (a, b, tol, msg) => ok(typeof a === 'number' && Math.abs(a - b) <= tol, `${msg} (got ${a}, expected ${b} ± ${tol})`);
const enc = new TextEncoder();
const bytesOf = (d) => (typeof d === 'string' ? enc.encode(d) : d instanceof Uint8Array ? new Uint8Array(d) : new Uint8Array(d));
const F = (name, data) => { const u8 = bytesOf(data), ab = u8.buffer.slice(u8.byteOffset, u8.byteOffset + u8.byteLength); return { name, size: u8.length, arrayBuffer: async () => ab, text: async () => new TextDecoder().decode(u8) }; };
const SAMPLES = [];
const imp = (name, data, opts) => { SAMPLES.push([name, bytesOf(data), opts]); return importGeometry(F(name, data), opts); };
async function group(label, fn) { try { await fn(); } catch (e) { fails++; console.log(`FAIL ${label} — threw: ${e && e.stack ? e.stack.split('\n').slice(0, 3).join(' | ') : e}`); } }
async function throws(label, fn, re) {
  try { await fn(); ok(false, `${label}: should have thrown`); }
  catch (e) { ok(e instanceof Error && e.message.length > 0 && e.message.length < 500 && (!re || re.test(e.message)), `${label}: clean error "${String(e && e.message).slice(0, 70)}"`); }
}
const cat = (...parts) => { const bs = parts.map(bytesOf), out = new Uint8Array(bs.reduce((s, b) => s + b.length, 0)); let o = 0; for (const b of bs) { out.set(b, o); o += b.length; } return out; };
const le = (type, vals, little = true) => { const sz = { u8: 1, i8: 1, u16: 2, i16: 2, u32: 4, i32: 4, f32: 4, f64: 8, u64: 8 }[type], b = new Uint8Array(sz * vals.length), dv = new DataView(b.buffer); vals.forEach((v, i) => { const o = i * sz; if (type === 'u8') dv.setUint8(o, v); else if (type === 'i8') dv.setInt8(o, v); else if (type === 'u16') dv.setUint16(o, v, little); else if (type === 'i16') dv.setInt16(o, v, little); else if (type === 'u32') dv.setUint32(o, v, little); else if (type === 'i32') dv.setInt32(o, v, little); else if (type === 'f32') dv.setFloat32(o, v, little); else if (type === 'u64') dv.setBigUint64(o, BigInt(v), little); else dv.setFloat64(o, v, little); }); return b; };
const be = (type, vals) => le(type, vals, false);
const b64 = (u8) => Buffer.from(u8).toString('base64');

// stored / deflated zip writer
const CRC = (() => { const t = new Uint32Array(256); for (let n = 0; n < 256; n++) { let c = n; for (let k = 0; k < 8; k++) c = c & 1 ? 0xedb88320 ^ (c >>> 1) : c >>> 1; t[n] = c >>> 0; } return t; })();
const crc32 = (u8) => { let c = 0xffffffff; for (const b of u8) c = CRC[(c ^ b) & 255] ^ (c >>> 8); return (c ^ 0xffffffff) >>> 0; };
function zip(entries, deflate = false) {
  const locals = [], central = [];
  let off = 0;
  for (const [name, data] of entries) {
    const raw = bytesOf(data), nm = enc.encode(name), body = deflate ? new Uint8Array(zlib.deflateRawSync(raw)) : raw, m = deflate ? 8 : 0;
    const lh = cat(le('u32', [0x04034b50]), le('u16', [20, 0, m, 0, 0]), le('u32', [crc32(raw), body.length, raw.length]), le('u16', [nm.length, 0]), nm, body);
    central.push(cat(le('u32', [0x02014b50]), le('u16', [20, 20, 0, m, 0, 0]), le('u32', [crc32(raw), body.length, raw.length]), le('u16', [nm.length, 0, 0, 0, 0]), le('u32', [0, off]), nm));
    locals.push(lh); off += lh.length;
  }
  const cd = cat(...central);
  return cat(...locals, cd, le('u32', [0x06054b50]), le('u16', [0, 0, entries.length, entries.length]), le('u32', [cd.length, off]), le('u16', [0]));
}

// unit cube: 8 vertices, 6 outward quads, 12 outward triangles, and a 6-tetrahedron split
const V = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]];
const Q = [[0, 3, 2, 1], [4, 5, 6, 7], [0, 1, 5, 4], [1, 2, 6, 5], [2, 3, 7, 6], [3, 0, 4, 7]];
const T = Q.flatMap((q) => [[q[0], q[1], q[2]], [q[0], q[2], q[3]]]);
const TETS = [[0, 1, 2, 6], [0, 2, 3, 6], [0, 3, 7, 6], [0, 7, 4, 6], [0, 4, 5, 6], [0, 5, 1, 6]];
const vtxt = (sep = ' ') => V.map((v) => v.join(sep)).join('\n');
function isCube(g, label, { count = 12, lo = [0, 0, 0], hi = [1, 1, 1], area = 6, volume = 1 } = {}) {
  ok(g.kind === 'mesh' && g.count === count && g.triangles.length === count * 9 && Array.isArray(g.triangles), `${label}: mesh with ${count} triangles (got ${g.kind}, ${g.count})`);
  ok(g.bbox.min.every((v, k) => Math.abs(v - lo[k]) < 1e-6) && g.bbox.max.every((v, k) => Math.abs(v - hi[k]) < 1e-6), `${label}: bbox ${JSON.stringify(g.bbox.min)}–${JSON.stringify(g.bbox.max)}`);
  const d = dimensions(g);
  ok(Math.abs(d.area - area) < 1e-6 && d.closed && Math.abs(d.volume - volume) < 1e-6 && d.consistentNormals, `${label}: area ${d.area}, closed ${d.closed}, volume ${d.volume}, consistent normals ${d.consistentNormals}`);
  ok(typeof g.name === 'string' && typeof g.format === 'string' && Array.isArray(g.warnings) && g.stats && typeof g.pathway === 'string', `${label}: carries name / format "${g.format}" / pathway "${g.pathway}"`);
}

// ---------------------------------------------------------------------------------------------------------
await group('registry', async () => {
  ok(Array.isArray(FORMATS) && FORMATS.length > 60 && FORMATS.every((f) => Array.isArray(f.ext) && f.ext.length && f.name && f.pathway && ['full', 'partial', 'convert'].includes(f.support) && typeof f.note === 'string' && (f.support !== 'convert' || (typeof f.convert === 'string' && f.convert.length > 20))), 'FORMATS entries are well formed');
  ok(FORMATS.every((f) => PATHWAYS[f.pathway]) && ['cad', 'surface', 'mesh', 'drawing', 'gis', 'points', 'voxel', 'network', 'numeric'].every((p) => PATHWAYS[p] && PATHWAYS[p].title && PATHWAYS[p].blurb), 'every format pathway is described in PATHWAYS');
  const ids = ['pvt', 'net', 'flow', 'solids', 'ops', 'integ', 'econ'], kindsOk = ['mesh', 'polylines', 'points', 'grid', 'voxels', 'network', 'table', 'params'];
  ok(ids.every((s) => SUITE_GEOMETRY[s] && SUITE_GEOMETRY[s].classes.length > 10 && SUITE_GEOMETRY[s].accepts.length && SUITE_GEOMETRY[s].accepts.every((k) => kindsOk.includes(k))) && Object.keys(SUITE_GEOMETRY).join() === ids.join(), 'SUITE_GEOMETRY covers the 7 suites with known geometry kinds');
  ok(!/membrane|spacer|desalin|brine|reverse osmosis/i.test(JSON.stringify([PATHWAYS, SUITE_GEOMETRY])) && PATHWAYS.well && /riser/.test(PATHWAYS.procedural.blurb), 'pathway and suite texts speak about pipelines, wells and risers');
  const reads = { graphml: 'GraphML', landxml: 'LandXML', inp: 'Abaqus input deck', cdb: 'ANSYS CDB archive', k: 'LS-DYNA keyword deck', key: 'LS-DYNA keyword deck', dyn: 'LS-DYNA keyword deck', dev: 'Well deviation survey (ASCII)', wbt: 'Well deviation survey (ASCII)', survey: 'Well deviation survey (ASCII)', sgy: 'SEG-Y seismic', segy: 'SEG-Y seismic', mat: 'MATLAB MAT-file', dem: 'DEM / DTM / DSM elevation model', dtm: 'DEM / DTM / DSM elevation model', ifczip: 'IFCZIP' };
  ok(Object.entries(reads).every(([e, nm]) => { const f = formatOf('file.' + e); return f && f.name === nm && f.support === 'partial' && f.note.length > 60; }), 'new readable formats are registered as partial with a note of what is and is not read');
  // every format family of the client specification is either read or has a conversion instruction
  const spec = 'step stp iges igs x_t x_b sat sab jt dxf dwg dgn ifc ifczip svg stl obj ply 3mf off wrl gltf glb x3d cgns vtk vtu vtp vts vtr vti exo msh med unv cas foam bdf nas inp cdb key k xyz csv tsv txt dat json xml yaml yml graphml shp geojson gpkg kml kmz gml landxml tif tiff dem dtm dsm nc asc las laz e57 sgy segy h5 hdf5 mat dev'.split(' ');
  ok(spec.every((e) => { const f = formatOf('a.' + e); return f && (f.support !== 'convert' || f.convert.length > 40); }), `all ${spec.length} extensions of the format specification are registered`);
  ok(/AP203/.test(formatOf('a.step').note) && /AP242/.test(formatOf('a.x_t').convert) && /ODA File Converter/.test(formatOf('a.dwg').convert) && /laszip/.test(formatOf('a.laz').convert) && /cgnsconvert/.test(formatOf('a.cgns').convert) && /gdal_translate/.test(formatOf('a.gpkg').convert), 'conversion instructions name a concrete tool or command');
  ok(formatOf('a.STEP').name === 'STEP (ISO 10303)' && formatOf('scan.nii.gz').name === 'NIfTI-1' && formatOf('part.x_t').support === 'partial' && formatOf('part.jt').support === 'convert' && formatOf('noext') === null && formatOf('a.unknownext') === null && formatOf('polyMesh/points').name === 'OpenFOAM polyMesh', 'formatOf handles case, double extensions, underscores and unknown names');
  const convExt = FORMATS.filter((f) => f.support === 'convert').flatMap((f) => f.ext), readExt = FORMATS.filter((f) => f.support !== 'convert').flatMap((f) => f.ext);
  ok(convExt.every((e) => !readExt.includes(e)), 'no extension is both readable and convert-only');
});

await group('surface formats', async () => {
  const facet = (t) => `facet normal 0 0 0\nouter loop\n${t.map((i) => 'vertex ' + V[i].join(' ')).join('\n')}\nendloop\nendfacet`;
  isCube(await imp('cube.stl', `solid cube\n${T.map(facet).join('\n')}\nendsolid cube\n`), 'STL ascii');
  const bin = new Uint8Array(84 + 50 * 12), dv = new DataView(bin.buffer);
  dv.setUint32(80, 12, true);
  T.forEach((t, i) => t.forEach((vi, k) => V[vi].forEach((c, q) => dv.setFloat32(84 + i * 50 + 12 + k * 12 + q * 4, c, true))));
  isCube(await imp('cube.stl', bin), 'STL binary');
  isCube(await imp('padded.stl', cat(bin, new Uint8Array(2))), 'STL binary with trailing padding bytes');
  isCube(await imp('cube.obj', `# cube\n${V.map((v) => 'v ' + v.join(' ')).join('\n')}\n${Q.map((q) => 'f ' + q.map((i) => i + 1).join(' ')).join('\n')}\n`), 'OBJ');
  const plyHead = (fmt) => `ply\nformat ${fmt} 1.0\ncomment test\nelement vertex 8\nproperty float x\nproperty float y\nproperty float z\nproperty uchar red\nelement face 6\nproperty list uchar int vertex_indices\nend_header\n`;
  isCube(await imp('cube.ply', plyHead('ascii') + V.map((v) => v.join(' ') + ' 255').join('\n') + '\n' + Q.map((q) => '4 ' + q.join(' ')).join('\n') + '\n'), 'PLY ascii');
  for (const little of [true, false]) {
    const body = cat(...V.map((v) => cat(le('f32', v, little), le('u8', [200]))), ...Q.map((q) => cat(le('u8', [4]), le('i32', q, little))));
    isCube(await imp('cube.ply', cat(plyHead(little ? 'binary_little_endian' : 'binary_big_endian'), body)), `PLY binary ${little ? 'little' : 'big'} endian`);
  }
  const pc = await imp('pts.ply', 'ply\nformat ascii 1.0\nelement vertex 3\nproperty double x\nproperty double y\nproperty double z\nend_header\n0 0 0\n1 2 3\n4 5 6\n');
  ok(pc.kind === 'points' && pc.count === 3 && pc.points[5] === 3 && pc.bbox.max[2] === 6, 'PLY without faces → point cloud');
  isCube(await imp('cube.off', `OFF\n# comment\n8 6 12\n${vtxt()}\n${Q.map((q) => '4 ' + q.join(' ')).join('\n')}\n`), 'OFF');
  const model = `<?xml version="1.0" encoding="UTF-8"?>\n<model unit="millimeter" xmlns="http://schemas.microsoft.com/3dmanufacturing/core/2015/02"><resources><object id="1" type="model"><mesh><vertices>${V.map((v) => `<vertex x="${v[0]}" y="${v[1]}" z="${v[2]}"/>`).join('')}</vertices><triangles>${T.map((t) => `<triangle v1="${t[0]}" v2="${t[1]}" v3="${t[2]}"/>`).join('')}</triangles></mesh></object></resources><build><item objectid="1" transform="1 0 0 0 1 0 0 0 1 10 20 30"/></build></model>`;
  isCube(await imp('cube.3mf', zip([['[Content_Types].xml', '<Types/>'], ['3D/3dmodel.model', model]])), '3MF (stored zip, build transform)', { lo: [10, 20, 30], hi: [11, 21, 31] });
  isCube(await imp('cube.3mf', zip([['3D/3dmodel.model', model]], true)), '3MF (deflated zip)', { lo: [10, 20, 30], hi: [11, 21, 31] });
  const amf = `<?xml version="1.0"?><amf unit="millimeter"><object id="0"><mesh><vertices>${V.map((v) => `<vertex><coordinates><x>${v[0]}</x><y>${v[1]}</y><z>${v[2]}</z></coordinates></vertex>`).join('')}</vertices><volume>${T.map((t) => `<triangle><v1>${t[0]}</v1><v2>${t[1]}</v2><v3>${t[2]}</v3></triangle>`).join('')}</volume></mesh></object></amf>`;
  isCube(await imp('cube.amf', amf), 'AMF xml');
  isCube(await imp('cube.amf', zip([['cube.amf', amf]])), 'AMF zip-wrapped');
  const ci = Q.map((q) => q.join(' ') + ' -1').join(' ');
  isCube(await imp('cube.wrl', `#VRML V2.0 utf8\n# a cube\nShape { geometry IndexedFaceSet { coordIndex [ ${Q.map((q) => q.join(', ') + ', -1').join(', ')} ]\n coord Coordinate { point [ ${V.map((v) => v.join(' ')).join(', ')} ] } } }\n`), 'VRML 2 IndexedFaceSet');
  isCube(await imp('cube.iv', `#Inventor V2.1 ascii\nSeparator { Coordinate3 { point [ ${V.map((v) => v.join(' ')).join(', ')} ] } IndexedFaceSet { coordIndex [ ${Q.map((q) => q.join(', ') + ', -1').join(', ')} ] } }\n`), 'Open Inventor ascii');
  isCube(await imp('cube.x3d', `<?xml version="1.0"?><!DOCTYPE X3D PUBLIC "ISO//Web3D//DTD X3D 3.3//EN" "http://www.web3d.org/specifications/x3d-3.3.dtd"><X3D><Scene><Shape><IndexedFaceSet coordIndex="${ci}"><Coordinate point="${V.map((v) => v.join(' ')).join(', ')}"/></IndexedFaceSet></Shape></Scene></X3D>`), 'X3D');
  const daeSrc = `<source id="pos"><float_array id="pos-a" count="24">${V.flat().join(' ')}</float_array><technique_common><accessor source="#pos-a" count="8" stride="3"/></technique_common></source><source id="nrm"><float_array id="nrm-a" count="3">0 0 1</float_array></source><vertices id="vtx"><input semantic="POSITION" source="#pos"/></vertices>`;
  const dae = (prim) => `<?xml version="1.0"?><COLLADA version="1.4.1"><asset><unit meter="0.01" name="cm"/><up_axis>Z_UP</up_axis></asset><library_geometries><geometry id="g"><mesh>${daeSrc}${prim}</mesh></geometry></library_geometries></COLLADA>`;
  isCube(await imp('cube.dae', dae(`<triangles count="12"><input semantic="VERTEX" source="#vtx" offset="0"/><input semantic="NORMAL" source="#nrm" offset="1"/><p>${T.flat().map((i) => i + ' 0').join(' ')}</p></triangles>`)), 'COLLADA triangles (2 inputs)');
  isCube(await imp('cube.dae', dae(`<polylist count="6"><input semantic="VERTEX" source="#vtx" offset="0"/><vcount>4 4 4 4 4 4</vcount><p>${Q.flat().join(' ')}</p></polylist>`)), 'COLLADA polylist');
  // glTF: positions (float32) + indices (uint16) in one buffer
  const gbin = cat(le('f32', V.flat()), le('u16', T.flat()));
  const gj = (uri, extra = {}) => ({ asset: { version: '2.0' }, scene: 0, scenes: [{ nodes: [0] }], nodes: [{ mesh: 0, ...extra }], meshes: [{ primitives: [{ attributes: { POSITION: 0 }, indices: 1 }] }], buffers: [uri ? { uri, byteLength: gbin.length } : { byteLength: gbin.length }], bufferViews: [{ buffer: 0, byteOffset: 0, byteLength: 96 }, { buffer: 0, byteOffset: 96, byteLength: 72 }], accessors: [{ bufferView: 0, componentType: 5126, count: 8, type: 'VEC3' }, { bufferView: 1, componentType: 5123, count: 36, type: 'SCALAR' }] });
  isCube(await imp('cube.gltf', JSON.stringify(gj('data:application/octet-stream;base64,' + b64(gbin)))), 'glTF embedded base64');
  isCube(await imp('cube.gltf', JSON.stringify(gj('data:application/octet-stream;base64,' + b64(gbin), { translation: [5, 0, 0], scale: [2, 2, 2] }))), 'glTF node translation + scale', { lo: [5, 0, 0], hi: [7, 2, 2], area: 24, volume: 8 });
  isCube(await imp('cube.gltf', JSON.stringify(gj('cube.bin')), { companion: { 'cube.bin': F('cube.bin', gbin) } }), 'glTF with companion .bin');
  await throws('glTF without its external buffer', () => imp('cube.gltf', JSON.stringify(gj('cube.bin'))), /companion|\.glb/);
  let js = enc.encode(JSON.stringify(gj(null, { rotation: [0, 0, Math.SQRT1_2, Math.SQRT1_2] })));
  js = cat(js, ' '.repeat((4 - (js.length % 4)) % 4));
  const glb = cat(le('u32', [0x46546c67, 2, 12 + 8 + js.length + 8 + gbin.length]), le('u32', [js.length, 0x4e4f534a]), js, le('u32', [gbin.length, 0x004e4942]), gbin);
  isCube(await imp('cube.glb', glb), 'GLB binary, node rotated 90° about z', { lo: [-1, 0, 0], hi: [0, 1, 1] });
  // VTK legacy
  const vtkHead = (kind, mode = 'ASCII') => `# vtk DataFile Version 3.0\ncube\n${mode}\nDATASET ${kind}\n`;
  isCube(await imp('cube.vtk', vtkHead('POLYDATA') + `POINTS 8 float\n${vtxt()}\nPOLYGONS 6 30\n${Q.map((q) => '4 ' + q.join(' ')).join('\n')}\n`), 'VTK legacy ASCII POLYDATA');
  isCube(await imp('pv.vtk', vtkHead('POLYDATA') + `FIELD FieldData 1\nTIME 1 1 double\n0.5\nPOINTS 8 float\n${vtxt()}\n\nMETADATA\nINFORMATION 2\nNAME L2_NORM_RANGE LOCATION vtkDataArray\nDATA 2 0 1.73205\nNAME L2_NORM_FINITE_RANGE LOCATION vtkDataArray\nDATA 2 0 1.73205\n\nPOLYGONS 6 30\n${Q.map((q) => '4 ' + q.join(' ')).join('\n')}\n\nPOINT_DATA 8\nSCALARS p float\nLOOKUP_TABLE default\n0 1 2 3 4 5 6 7\n`), 'VTK legacy with FIELD and METADATA blocks (ParaView style)');
  isCube(await imp('cube.vtk', cat(vtkHead('POLYDATA', 'BINARY'), 'POINTS 8 float\n', be('f32', V.flat()), '\nPOLYGONS 6 30\n', be('i32', Q.flatMap((q) => [4, ...q])), '\n')), 'VTK legacy BINARY POLYDATA');
  const ug = await imp('hex.vtk', vtkHead('UNSTRUCTURED_GRID') + `POINTS 8 double\n${vtxt()}\nCELLS 1 9\n8 0 1 2 3 4 5 6 7\nCELL_TYPES 1\n12\n`);
  isCube(ug, 'VTK legacy UNSTRUCTURED_GRID hexahedron → boundary');
  ok(ug.stats.cells === 1 && ug.stats.nodes === 8 && ug.stats.cellTypes.hex === 1 && ug.stats.dimension === 3, 'VTK volume-mesh stats (nodes, cells, cellTypes, dimension)');
  isCube(await imp('cube51.vtk', `# vtk DataFile Version 5.1\ncube\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS 8 float\n${vtxt()}\nCELLS 7 24\nOFFSETS vtktypeint64\n0 4 8 12 16 20 24\nCONNECTIVITY vtktypeint64\n${TETS.flat().join(' ')}\nCELL_TYPES 6\n10 10 10 10 10 10\n`), 'VTK legacy 5.1 offsets/connectivity, six tetrahedra');
  const sp = await imp('img.vtk', vtkHead('STRUCTURED_POINTS') + 'DIMENSIONS 3 2 2\nORIGIN 0 0 0\nSPACING 0.5 0.5 0.5\nPOINT_DATA 12\nSCALARS phase unsigned_char 1\nLOOKUP_TABLE default\n1 0 0 1 0 0 1 0 0 1 0 0\n');
  ok(sp.kind === 'voxels' && sp.voxels.nx === 3 && sp.voxels.ny === 2 && sp.voxels.nz === 2 && sp.voxels.data[0] === 1 && sp.voxels.data[1] === 0 && Math.abs(sp.stats.porosity - 2 / 3) < 1e-12 && sp.voxels.spacing[0] === 0.5, 'VTK STRUCTURED_POINTS → voxels');
  isCube(await imp('grid.vtk', vtkHead('RECTILINEAR_GRID') + 'DIMENSIONS 3 2 2\nX_COORDINATES 3 float\n0 0.5 1\nY_COORDINATES 2 float\n0 1\nZ_COORDINATES 2 float\n0 1\n'), 'VTK RECTILINEAR_GRID → boundary', { count: 20 });
  // VTK XML
  const vtp = `<?xml version="1.0"?><VTKFile type="PolyData" version="0.1" byte_order="LittleEndian"><PolyData><Piece NumberOfPoints="8" NumberOfPolys="6"><Points><DataArray type="Float32" NumberOfComponents="3" format="ascii">${V.flat().join(' ')}</DataArray></Points><Polys><DataArray type="Int32" Name="connectivity" format="ascii">${Q.flat().join(' ')}</DataArray><DataArray type="Int32" Name="offsets" format="ascii">4 8 12 16 20 24</DataArray></Polys></Piece></PolyData></VTKFile>`;
  isCube(await imp('cube.vtp', vtp), 'VTP ascii');
  const hdr = (n, t = 'u32') => le(t, [n]);
  const vtu = (attrs, pts, con, off, typ, tail = '') => `<?xml version="1.0"?><VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian" ${attrs}><UnstructuredGrid><Piece NumberOfPoints="8" NumberOfCells="6"><Points>${pts}</Points><Cells>${con}${off}${typ}</Cells></Piece></UnstructuredGrid>${tail}</VTKFile>`;
  const da = (type, name, fmt, body, extra = '') => `<DataArray type="${type}" Name="${name}" ${name === 'Points' ? 'NumberOfComponents="3" ' : ''}format="${fmt}" ${extra}>${body}</DataArray>`;
  const P8 = le('f64', V.flat()), C8 = le('u64', TETS.flat()), O8 = le('u64', [4, 8, 12, 16, 20, 24]), Ty = le('u8', [10, 10, 10, 10, 10, 10]);
  const raw = (b, sep = true, ht = 'u32') => (sep ? b64(hdr(b.length, ht)) + b64(b) : b64(cat(hdr(b.length, ht), b)));
  isCube(await imp('cube.vtu', vtu('', da('Float64', 'Points', 'ascii', V.flat().join(' ')), da('Int64', 'connectivity', 'ascii', TETS.flat().join(' ')), da('Int64', 'offsets', 'ascii', '4 8 12 16 20 24'), da('UInt8', 'types', 'ascii', '10 10 10 10 10 10'))), 'VTU ascii (six tetrahedra)');
  isCube(await imp('cube.vtu', vtu('header_type="UInt32"', da('Float64', 'Points', 'binary', raw(P8)), da('Int64', 'connectivity', 'binary', raw(C8)), da('Int64', 'offsets', 'binary', raw(O8)), da('UInt8', 'types', 'binary', raw(Ty)))), 'VTU inline base64 (header encoded separately)');
  isCube(await imp('cube.vtu', vtu('header_type="UInt64"', da('Float64', 'Points', 'binary', raw(P8, false, 'u64')), da('Int64', 'connectivity', 'binary', raw(C8, false, 'u64')), da('Int64', 'offsets', 'binary', raw(O8, false, 'u64')), da('UInt8', 'types', 'binary', raw(Ty, false, 'u64')))), 'VTU inline base64 (UInt64 header, one stream)');
  const zl = (b) => { const c = new Uint8Array(zlib.deflateSync(b)); return b64(le('u32', [1, b.length, b.length, c.length])) + b64(c); };
  isCube(await imp('cube.vtu', vtu('compressor="vtkZLibDataCompressor"', da('Float64', 'Points', 'binary', zl(P8)), da('Int64', 'connectivity', 'binary', zl(C8)), da('Int64', 'offsets', 'binary', zl(O8)), da('UInt8', 'types', 'binary', zl(Ty)))), 'VTU inline base64, zlib-compressed');
  const blocks = [P8, C8, O8, Ty].map((b) => cat(hdr(b.length), b)), offs = blocks.map((_, k) => blocks.slice(0, k).reduce((s, b) => s + b.length, 0));
  const app = vtu('', da('Float64', 'Points', 'appended', '', `offset="${offs[0]}"`), da('Int64', 'connectivity', 'appended', '', `offset="${offs[1]}"`), da('Int64', 'offsets', 'appended', '', `offset="${offs[2]}"`), da('UInt8', 'types', 'appended', '', `offset="${offs[3]}"`), '<AppendedData encoding="raw">_@@@</AppendedData>');
  const [a1, a2] = app.split('@@@');
  isCube(await imp('cube.vtu', cat(a1, ...blocks, a2)), 'VTU appended raw data');
  const b64blocks = blocks.map((b) => b64(b)), boffs = b64blocks.map((_, k) => b64blocks.slice(0, k).reduce((s, b) => s + b.length, 0));
  isCube(await imp('cube.vtu', vtu('', da('Float64', 'Points', 'appended', '', `offset="${boffs[0]}"`), da('Int64', 'connectivity', 'appended', '', `offset="${boffs[1]}"`), da('Int64', 'offsets', 'appended', '', `offset="${boffs[2]}"`), da('UInt8', 'types', 'appended', '', `offset="${boffs[3]}"`), `<AppendedData encoding="base64">_${b64blocks.join('')}</AppendedData>`)), 'VTU appended base64 data');
  const vti = await imp('img.vti', `<VTKFile type="ImageData" byte_order="LittleEndian"><ImageData WholeExtent="0 2 0 1 0 1" Origin="0 0 0" Spacing="2 2 2"><Piece Extent="0 2 0 1 0 1"><PointData Scalars="s"><DataArray type="UInt8" Name="s" format="ascii">0 0 1 0 0 1 0 0 1 0 0 1</DataArray></PointData></Piece></ImageData></VTKFile>`);
  ok(vti.kind === 'voxels' && vti.voxels.nx === 3 && vti.voxels.data[2] === 1 && vti.voxels.spacing[1] === 2 && Math.abs(vti.stats.porosity - 2 / 3) < 1e-12, 'VTI ImageData → voxels');
  isCube(await imp('grid.vts', `<VTKFile type="StructuredGrid" byte_order="LittleEndian"><StructuredGrid WholeExtent="0 1 0 1 0 1"><Piece Extent="0 1 0 1 0 1"><Points><DataArray type="Float32" NumberOfComponents="3" format="ascii">0 0 0 1 0 0 0 1 0 1 1 0 0 0 1 1 0 1 0 1 1 1 1 1</DataArray></Points></Piece></StructuredGrid></VTKFile>`), 'VTS StructuredGrid → boundary');
  isCube(await imp('grid.vtr', `<VTKFile type="RectilinearGrid" byte_order="LittleEndian"><RectilinearGrid WholeExtent="0 1 0 1 0 1"><Piece Extent="0 1 0 1 0 1"><Coordinates><DataArray type="Float32" format="ascii">0 1</DataArray><DataArray type="Float32" format="ascii">0 1</DataArray><DataArray type="Float32" format="ascii">0 1</DataArray></Coordinates></Piece></RectilinearGrid></VTKFile>`), 'VTR RectilinearGrid → boundary');
  // GTS: vertices, edges, faces (all 1-based)
  const em = new Map(), edges = [], eid = (a, b) => { const k = a < b ? a + '_' + b : b + '_' + a; if (!em.has(k)) { edges.push([a, b]); em.set(k, edges.length); } return em.get(k); };
  const gf = T.map((t) => [eid(t[0], t[1]), eid(t[1], t[2]), eid(t[2], t[0])]);
  isCube(await imp('cube.gts', `8 ${edges.length} 12\n${vtxt()}\n${edges.map((e) => `${e[0] + 1} ${e[1] + 1}`).join('\n')}\n${gf.map((f) => f.join(' ')).join('\n')}\n`), 'GTS');
  isCube(await imp('cube.byu', `1 8 6 24\n1 6\n${vtxt()}\n${Q.map((q) => q.map((i, k) => (k === 3 ? -(i + 1) : i + 1)).join(' ')).join('\n')}\n`), 'BYU');
  isCube(await imp('cube.tri', `8 12\n${vtxt()}\n${T.map((t) => t.map((i) => i + 1).join(' ')).join('\n')}\n`), 'TRI (Cart3D-style counts line)');
  isCube(await imp('cube.fac', `8\n${V.map((v, k) => `${k + 1} ${v.join(' ')}`).join('\n')}\n12\n${T.map((t, k) => `${k + 1} ${t.map((i) => i + 1).join(' ')}`).join('\n')}\n`), 'FAC (separate counts, indexed rows)');
  await throws('unrecognised .surf dialect', () => imp('x.surf', 'surface data v2\nfoo bar\n'), /STL|OBJ|PLY/);
});

await group('computational meshes', async () => {
  const nodes22 = V.map((v, k) => `${k + 1} ${v.join(' ')}`).join('\n');
  const g2 = await imp('cube.msh', `$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n8\n${nodes22}\n$EndNodes\n$Elements\n8\n1 15 2 0 1 1\n2 3 2 0 1 1 4 3 2\n${TETS.map((t, k) => `${k + 3} 4 2 0 1 ${t.map((i) => i + 1).join(' ')}`).join('\n')}\n$EndElements\n`);
  isCube(g2, 'Gmsh 2.2 (six tetrahedra + boundary quad + point)');
  ok(g2.format === 'Gmsh MSH' && g2.stats.cellTypes.tet === 6 && g2.stats.cellTypes.quad === 1 && g2.stats.dimension === 3 && g2.stats.nodes === 8, 'Gmsh 2.2 stats');
  isCube(await imp('cube.msh', `$MeshFormat\n4.1 0 8\n$EndMeshFormat\n$Entities\n0 0 0 1\n1 0 0 0 1 1 1 0 0\n$EndEntities\n$Nodes\n1 8 1 8\n3 1 0 8\n1\n2\n3\n4\n5\n6\n7\n8\n${vtxt()}\n$EndNodes\n$Elements\n1 1 1 1\n3 1 5 1\n1 1 2 3 4 5 6 7 8\n$EndElements\n`), 'Gmsh 4.1 (one hexahedron)');
  isCube(await imp('cube.msh', `$MeshFormat\n4 0 8\n$EndMeshFormat\n$Nodes\n1 8\n1 3 0 8\n${nodes22}\n$EndNodes\n$Elements\n1 1\n1 3 5 1\n1 1 2 3 4 5 6 7 8\n$EndElements\n`), 'Gmsh 4.0 (one hexahedron)');
  const sq = await imp('square.msh', `$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n4\n1 0 0 0\n2 2 0 0\n3 2 2 0\n4 0 2 0\n$EndNodes\n$Elements\n2\n1 2 2 0 1 1 2 3\n2 2 2 0 1 1 3 4\n$EndElements\n`);
  ok(sq.kind === 'polylines' && sq.polylines.length === 1 && sq.polylines[0].closed && sq.polylines[0].x.length === 4 && sq.stats.dimension === 2 && Math.abs(dimensions(sq).area - 4) < 1e-12, '2-D Gmsh mesh → closed boundary outline (area 4)');
  await throws('binary Gmsh', () => imp('b.msh', '$MeshFormat\n2.2 1 8\n\x01\x00\x00\x00\n$EndMeshFormat\n'), /ASCII/);
  // Fluent ASCII: one hexahedral cell, six boundary quads (hex indices, c0 = cell 1, c1 = 0)
  const fl = `(0 "test mesh")\n(2 3)\n(10 (0 1 8 0 3))\n(12 (0 1 1 0))\n(13 (0 1 6 0))\n(10 (1 1 8 1 3)(\n${vtxt()}\n))\n(13 (2 1 6 3 4)(\n${Q.map((q) => q.map((i) => (i + 1).toString(16)).join(' ') + ' 1 0').join('\n')}\n))\n(12 (3 1 1 1 4))\n`;
  const gfl = await imp('cube.msh', fl);
  isCube(gfl, 'Fluent ASCII mesh (sniffed from .msh)');
  ok(gfl.format === 'Fluent mesh (ASCII)' && gfl.stats.boundaryFaces === 6 && gfl.stats.cells === 1, 'Fluent stats and format name');
  isCube(await imp('cube.su2', `% test\nNDIME= 3\nNELEM= 6\n${TETS.map((t, k) => `10 ${t.join(' ')} ${k}`).join('\n')}\nNPOIN= 8\n${V.map((v, k) => v.join(' ') + ' ' + k).join('\n')}\nNMARK= 1\nMARKER_TAG= wall\nMARKER_ELEMS= 2\n5 0 1 2\n5 0 2 3\n`), 'SU2');
  const s2 = await imp('sq.su2', 'NDIME= 2\nNELEM= 1\n9 0 1 2 3 0\nNPOIN= 4\n0 0 0\n1 0 1\n1 1 2\n0 1 3\n');
  ok(s2.kind === 'polylines' && s2.polylines[0].closed && Math.abs(dimensions(s2).area - 1) < 1e-12, 'SU2 2-D quad → outline');
  const pad = (v, w = 10) => String(v).padStart(w);
  const unv = `    -1\n  2411\n${V.map((v, k) => `${pad(k + 1)}${pad(1)}${pad(1)}${pad(11)}\n${v.map((c) => c.toExponential(16).replace('e', 'D').padStart(25)).join('')}`).join('\n')}\n    -1\n    -1\n  2412\n${pad(1)}${pad(115)}${pad(1)}${pad(1)}${pad(7)}${pad(8)}\n${[1, 2, 3, 4, 5, 6, 7, 8].map((i) => pad(i)).join('')}\n${pad(2)}${pad(94)}${pad(1)}${pad(1)}${pad(7)}${pad(4)}\n${[1, 4, 3, 2].map((i) => pad(i)).join('')}\n    -1\n`;
  const gu = await imp('cube.unv', unv);
  isCube(gu, 'UNV (2411 nodes, 2412 hexahedron + shell)');
  ok(gu.stats.cellTypes.hex === 1 && gu.stats.cellTypes.quad === 1, 'UNV element types');
  const f8 = (s) => String(s).padEnd(8).slice(0, 8), f16 = (s) => String(s).padEnd(16).slice(0, 16);
  const bdfFree = `$ free field\nBEGIN BULK\n${V.map((v, k) => `GRID,${k + 1},,${v.map((c) => c.toFixed(1)).join(',')}`).join('\n')}\nCHEXA,1,1,1,2,3,4,5,6,+C1\n+C1,7,8\nENDDATA\n`;
  isCube(await imp('cube.bdf', bdfFree), 'Nastran free field (CHEXA with continuation)');
  const bdfSmall = `SOL 101\nCEND\nBEGIN BULK\n${V.map((v, k) => f8('GRID') + f8(k + 1) + f8('') + v.map((c) => f8(c ? '1.0' : '0.')).join('')).join('\n')}\n${TETS.map((t, k) => f8('CTETRA') + f8(k + 1) + f8(1) + t.map((i) => f8(i + 1)).join('')).join('\n')}\nENDDATA\n`;
  isCube(await imp('cube.nas', bdfSmall), 'Nastran small field (CTETRA)');
  const bdfLarge = `BEGIN BULK\n${V.map((v, k) => f8('GRID*') + f16(k + 1) + f16('') + f16(v[0] ? '1.+0' : '0.0') + f16(v[1] ? '1.0D0' : '0.0') + f8('*G' + k) + '\n' + f8('*G' + k) + f16(v[2] ? '10.-1' : '0.0')).join('\n')}\n${f8('CHEXA')}${f8(1)}${f8(1)}${[1, 2, 3, 4, 5, 6].map(f8).join('')}${f8('+H')}\n${f8('+H')}${f8(7)}${f8(8)}\nENDDATA\n`;
  isCube(await imp('cube.bdf', bdfLarge), 'Nastran large-field GRID* with Nastran-style exponents');
  const tecFE = `TITLE = "cube"\nVARIABLES = "X", "Y", "Z", "P"\nZONE T="solid", N=8, E=1, DATAPACKING=POINT, ZONETYPE=FEBRICK\n${V.map((v) => v.join(' ') + ' 1.5').join('\n')}\n1 2 3 4 5 6 7 8\n`;
  isCube(await imp('cube.dat', tecFE), 'Tecplot FEBRICK (point packing, sniffed from .dat)');
  const tecBlock = `VARIABLES = "X" "Y" "Z"\nZONE N=8, E=6, F=FEBLOCK, ET=TETRAHEDRON\n${[0, 1, 2].map((c) => V.map((v) => v[c]).join(' ')).join('\n')}\n${TETS.map((t) => t.map((i) => i + 1).join(' ')).join('\n')}\n`;
  isCube(await imp('cube.tec', tecBlock), 'Tecplot FEBLOCK tetrahedra');
  const ijk = []; for (let k = 0; k < 2; k++) for (let j = 0; j < 2; j++) for (let i = 0; i < 3; i++) ijk.push([i / 2, j, k]);
  const tecO = await imp('cube.plt', `TITLE="ordered"\nVARIABLES="X","Y","Z"\nZONE I=3, J=2, K=2, F=POINT\n${ijk.map((p) => p.join(' ')).join('\n')}\n`);
  isCube(tecO, 'Tecplot ordered I×J×K zone (sniffed from .plt)', { count: 20 });
  ok(tecO.format === 'Tecplot ASCII' && tecO.stats.structured && tecO.stats.cells === 2, 'Tecplot ordered zone stats');
  await throws('Tecplot with an absurd VARLOCATION range', () => imp('cc.dat', 'VARIABLES = "X", "Y", "Z", "P"\nZONE N=8, E=1, DATAPACKING=BLOCK, ZONETYPE=FEBRICK, VARLOCATION=([1-999999999999]=CELLCENTERED)\n0 1 1 0 0 1 1 0\n0 0 1 1 0 0 1 1\n0 0 0 0 1 1 1 1\n1.5\n1 2 3 4 5 6 7 8\n'), /Cell-centred/);
  const tecCC = await imp('cc.dat', 'VARIABLES = "X", "Y", "Z", "P"\nZONE N=8, E=1, DATAPACKING=BLOCK, ZONETYPE=FEBRICK, VARLOCATION=([4]=CELLCENTERED)\n0 1 1 0 0 1 1 0\n0 0 1 1 0 0 1 1\n0 0 0 0 1 1 1 1\n1.5\n1 2 3 4 5 6 7 8\n');
  isCube(tecCC, 'Tecplot block packing with a cell-centred variable');
  const tec2 = await imp('sq.dat', 'VARIABLES = "X", "Y"\nZONE N=4, E=1, F=FEPOINT, ET=QUADRILATERAL\n0 0\n3 0\n3 1\n0 1\n1 2 3 4\n');
  ok(tec2.kind === 'polylines' && tec2.polylines[0].closed && Math.abs(dimensions(tec2).area - 3) < 1e-12, 'Tecplot 2-D quadrilateral → outline');
  await throws('binary Tecplot .plt', () => imp('b.plt', cat('#!TDV112', new Uint8Array(40))), /ASCII/);
  const p3 = [0, 1, 2].map((c) => ijk.map((p) => p[c]).join(' ')).join('\n');
  isCube(await imp('grid.p3d', `3 2 2\n${p3}\n`), 'Plot3D single block', { count: 20 });
  isCube(await imp('grid.x', `1\n3 2 2\n${p3}\n`), 'Plot3D multi-block header', { count: 20 });
  const px = await imp('grid.xyz', `3 2 2\n${p3}\n`);
  ok(px.kind === 'mesh' && px.format === 'Plot3D grid (ASCII)' && px.count === 20, '.xyz with a Plot3D header is read as a structured grid');
  isCube(await imp('mesh.json', JSON.stringify({ nodes: V, elements: [[0, 1, 2, 3, 4, 5, 6, 7]] })), 'JSON nodes + elements (hexahedron)');
  isCube(await imp('mesh.json', JSON.stringify({ vertices: V, faces: Q.map((q) => q.map((i) => i + 1)) })), 'JSON vertices + faces (1-based quads)');
  isCube(await imp('mesh.json', JSON.stringify({ nodes: V, cells: TETS })), 'JSON nodes + cells (tetrahedra recognised by shape)');
  // OpenFOAM polyMesh
  const foam = (cls, body, note = '') => `/*--------------------------------*- C++ -*----------------------------------*\\\n| OpenFOAM |\n\\*---------------------------------------------------------------------------*/\nFoamFile\n{\n    version 2.0;\n    format ascii;\n    class ${cls};\n    ${note}\n    object x;\n}\n// * * * //\n${body}\n`;
  const fp = F('points', foam('vectorField', `8\n(\n${V.map((v) => `(${v.join(' ')})`).join('\n')}\n)`)), ff = F('faces', foam('faceList', `6\n(\n${Q.map((q) => `4(${q.join(' ')})`).join('\n')}\n)`)), fo = F('owner', foam('labelList', '6\n(\n0\n0\n0\n0\n0\n0\n)', 'note "nPoints: 8 nCells: 1 nFaces: 6 nInternalFaces: 0";')), fn = F('neighbour', foam('labelList', '0\n(\n)'));
  const gfo = await importGeometry(fp, { companion: { faces: ff, owner: fo, neighbour: fn } });
  isCube(gfo, 'OpenFOAM polyMesh (points + faces + owner + neighbour)');
  ok(gfo.stats.cells === 1 && gfo.stats.internalFaces === 0, 'OpenFOAM stats');
  isCube(await importGeometry(ff, { companion: [fp, fo] }), 'OpenFOAM started from faces, companions as a list, nInternalFaces from the owner note');
  await throws('OpenFOAM points without companions', () => importGeometry(fp), /faces/);
});

// ---- CAD writers ----------------------------------------------------------------------------------------
const RN = (v) => { const s = String(v); return /[.e]/.test(s) ? s : s + '.'; };
function stepFile(build, schema = 'AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }') {
  let id = 0;
  const lines = [], E = (s) => { lines.push(`#${++id}=${s};`); return '#' + id; };
  build(E);
  return `ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('test; with ''quotes'''),'2;1');\nFILE_NAME('t.stp','2024-01-01T00:00:00',('a'),('o'),'p','s','');\nFILE_SCHEMA(('${schema}'));\nENDSEC;\nDATA;\n/* geometry */\n${lines.join('\n')}\nENDSEC;\nEND-ISO-10303-21;\n`;
}
const sP = (E, p) => E(`CARTESIAN_POINT('',(${p.map(RN).join(',')}))`);
const sD = (E, d) => E(`DIRECTION('',(${d.map(RN).join(',')}))`);
const sAX = (E, o, z, x) => E(`AXIS2_PLACEMENT_3D('',${sP(E, o)},${sD(E, z)},${sD(E, x)})`);
const sVX = (E, p) => E(`VERTEX_POINT('',${sP(E, p)})`);
const sLine = (E, va, vb, pa, pb) => E(`EDGE_CURVE('',${va},${vb},${E(`LINE('',${sP(E, pa)},${E(`VECTOR('',${sD(E, pb.map((v, k) => v - pa[k]))},1.)`)})`)},.T.)`);
const sLoop = (E, edges) => E(`EDGE_LOOP('',(${edges.map(([e, o]) => E(`ORIENTED_EDGE('',*,*,${e},${o ? '.T.' : '.F.'})`)).join(',')}))`);
const sUnits = (E) => E('(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.))');
const crossv = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
/** Planar quad faces (outward, counter-clockwise) -> ADVANCED_FACE refs. */
function sQuadFaces(E, pts, quads) {
  return quads.map((q) => {
    const p = q.map((i) => pts[i]), vs = p.map((x) => sVX(E, x)), es = p.map((x, k) => [sLine(E, vs[k], vs[(k + 1) % 4], x, p[(k + 1) % 4]), true]);
    const u = p[1].map((v, k) => v - p[0][k]), w = p[3].map((v, k) => v - p[0][k]), n = crossv(u, w), l = Math.hypot(...n), ul = Math.hypot(...u);
    return E(`ADVANCED_FACE('',(${E(`FACE_OUTER_BOUND('',${sLoop(E, es)},.T.)`)}),${E(`PLANE('',${sAX(E, p[0], n.map((v) => v / l), u.map((v) => v / ul))})`)},.T.)`);
  });
}
function igesFile(entities) {
  const row = (body, sec, n) => body.padEnd(72).slice(0, 72) + sec + String(n).padStart(7), f = (v) => String(v).padStart(8);
  const chunk = (s, w) => { const out = []; let cur = ''; for (const part of s.split(/(?<=[,;])/)) { if (cur.length + part.length > w) { out.push(cur); cur = ''; } cur += part; } if (cur) out.push(cur); return out; };
  const G = chunk('1H,,1H;,4HTEST,8HTEST.IGS,6HSYSTEM,3H1.0,32,38,6,308,15,4HTEST,1.,2,2HMM,1,0.01,15H20240101.000000,1E-06,100.,4HAUTH,3HORG,11,0,15H20240101.000000;', 72);
  const D = [], P = [];
  entities.forEach((e, i) => {
    const pl = chunk(e.params + ';', 64), start = P.length + 1;
    for (const s of pl) P.push(s.padEnd(64) + ' ' + String(2 * i + 1).padStart(7) + 'P' + String(P.length + 1).padStart(7));
    D.push(f(e.type) + f(start) + f(0) + f(0) + f(0) + f(0) + f(e.xf || 0) + f(0) + (e.status || '00000000') + 'D' + String(2 * i + 1).padStart(7));
    D.push(f(e.type) + f(0) + f(0) + f(pl.length) + f(e.form || 0) + f('') + f('') + f('') + f(0) + 'D' + String(2 * i + 2).padStart(7));
  });
  return [row('IGES test file', 'S', 1), ...G.map((g, k) => row(g, 'G', k + 1)), ...D, ...P, row(`S${String(1).padStart(7)}G${String(G.length).padStart(7)}D${String(D.length).padStart(7)}P${String(P.length).padStart(7)}`, 'T', 1)].join('\n') + '\n';
}
const A48 = 24 * Math.sin((2 * Math.PI) / 48), P48 = 96 * Math.sin(Math.PI / 48);   // area and perimeter of the 48-gon inscribed in a unit circle

await group('CAD: STEP and IGES', async () => {
  const B = V.map((v) => [v[0] * 2, v[1] * 3, v[2] * 4]);
  const box = await imp('box.step', stepFile((E) => { sUnits(E); const fs = sQuadFaces(E, B, Q); E(`MANIFOLD_SOLID_BREP('box',${E(`CLOSED_SHELL('',(${fs.join(',')}))`)})`); }));
  isCube(box, 'STEP box (ADVANCED_FACE on PLANE)', { hi: [2, 3, 4], area: 52, volume: 24 });
  ok(box.stats.units === 'mm' && box.stats.faces === 6 && box.stats.tessellatedFaces === 6 && box.format === 'STEP (ISO 10303)' && box.pathway === 'cad', 'STEP box: units mm, 6 of 6 faces, format and pathway');
  const part = await imp('part.stp', stepFile((E) => {
    sUnits(E); sQuadFaces(E, B, Q);
    const far = sP(E, [10, 11, 12]);
    E(`ADVANCED_FACE('',(),${E(`B_SPLINE_SURFACE_WITH_KNOTS('',1,1,((${far},${far}),(${far},${far})),.UNSPECIFIED.,.F.,.F.,.F.,(2,2),(2,2),(0.,1.),(0.,1.),.UNSPECIFIED.)`)},.T.)`);
  }));
  ok(part.kind === 'mesh' && part.count === 12 && part.warnings.some((w) => /6 of 7 faces tessellated/.test(w) && /b spline surface/.test(w)), 'STEP with a B-spline face: "6 of 7 faces tessellated" reported');
  ok(part.bbox.max[0] === 10 && part.bbox.max[2] === 12 && part.bbox.min[0] === 0, 'STEP bbox still covers every CARTESIAN_POINT when faces are skipped');
  // solid cylinder R = 1, h = 2
  const cyl = await imp('cyl.stp', stepFile((E) => {
    sUnits(E);
    const v0 = sVX(E, [1, 0, 0]), v1 = sVX(E, [1, 0, 2]), ax0 = sAX(E, [0, 0, 0], [0, 0, 1], [1, 0, 0]), ax1 = sAX(E, [0, 0, 2], [0, 0, 1], [1, 0, 0]);
    const e0 = E(`EDGE_CURVE('',${v0},${v0},${E(`CIRCLE('',${ax0},1.)`)},.T.)`), e1 = E(`EDGE_CURVE('',${v1},${v1},${E(`CIRCLE('',${ax1},1.)`)},.T.)`), e2 = sLine(E, v0, v1, [1, 0, 0], [1, 0, 2]);
    E(`ADVANCED_FACE('',(${E(`FACE_OUTER_BOUND('',${sLoop(E, [[e0, false]])},.T.)`)}),${E(`PLANE('',${ax0})`)},.F.)`);
    E(`ADVANCED_FACE('',(${E(`FACE_OUTER_BOUND('',${sLoop(E, [[e1, true]])},.T.)`)}),${E(`PLANE('',${ax1})`)},.T.)`);
    E(`ADVANCED_FACE('',(${E(`FACE_OUTER_BOUND('',${sLoop(E, [[e0, true], [e2, true], [e1, false], [e2, false]])},.T.)`)}),${E(`CYLINDRICAL_SURFACE('',${ax0},1.)`)},.T.)`);
  }));
  const dc = dimensions(cyl);
  ok(cyl.kind === 'mesh' && cyl.count === 46 * 2 + 96 && cyl.stats.tessellatedFaces === 3, `STEP cylinder: 3 faces → ${cyl.count} triangles`);
  ok(dc.closed && dc.consistentNormals && Math.abs(dc.volume - 2 * A48) < 1e-9 && Math.abs(dc.area - (2 * A48 + 2 * P48)) < 1e-9, `STEP cylinder: closed, volume ${dc.volume} (π·r²·h = ${2 * Math.PI}), area ${dc.area}`);
  ok(Math.abs(dc.volume - 2 * Math.PI) / (2 * Math.PI) < 0.005 && cyl.bbox.min.every((v, k) => Math.abs(v - [-1, -1, 0][k]) < 1e-9) && cyl.bbox.max.every((v, k) => Math.abs(v - [1, 1, 2][k]) < 1e-9), 'STEP cylinder: volume within 0.5 % of exact, bbox [-1,-1,0]–[1,1,2]');
  // plates of thickness 1 with N × N round holes of radius 1 on a pitch of 4: planar faces with holes (several bridged to the
  // same outline vertex) + inward-facing cylinders (same_sense .F.)
  for (const N of [1, 3]) {
    const Wd = 4 * N, plate = await imp(`plate${N}.stp`, stepFile((E) => {
      sUnits(E);
      const C = [[0, 0, 0], [Wd, 0, 0], [Wd, Wd, 0], [0, Wd, 0], [0, 0, 1], [Wd, 0, 1], [Wd, Wd, 1], [0, Wd, 1]], rings = [];
      sQuadFaces(E, C, Q.slice(2));
      for (let ia = 0; ia < N; ia++) for (let ib = 0; ib < N; ib++) {
        const cx = 2 + 4 * ia, cy = 2 + 4 * ib, w0 = sVX(E, [cx + 1, cy, 0]), w1 = sVX(E, [cx + 1, cy, 1]), ax0 = sAX(E, [cx, cy, 0], [0, 0, 1], [1, 0, 0]), ax1 = sAX(E, [cx, cy, 1], [0, 0, 1], [1, 0, 0]);
        const c0 = E(`EDGE_CURVE('',${w0},${w0},${E(`CIRCLE('',${ax0},1.)`)},.T.)`), c1 = E(`EDGE_CURVE('',${w1},${w1},${E(`CIRCLE('',${ax1},1.)`)},.T.)`), seam = sLine(E, w0, w1, [cx + 1, cy, 0], [cx + 1, cy, 1]);
        rings.push([c0, c1]);
        E(`ADVANCED_FACE('',(${E(`FACE_BOUND('',${sLoop(E, [[c0, true], [seam, true], [c1, false], [seam, false]])},.T.)`)}),${E(`CYLINDRICAL_SURFACE('',${ax0},1.)`)},.F.)`);
      }
      for (const [zq, q, sense] of [[0, Q[0], '.F.'], [1, Q[1], '.T.']]) {
        const p = q.map((i) => C[i]), vs = p.map((x) => sVX(E, x)), es = p.map((x, k) => [sLine(E, vs[k], vs[(k + 1) % 4], x, p[(k + 1) % 4]), true]);
        E(`ADVANCED_FACE('',(${[E(`FACE_OUTER_BOUND('',${sLoop(E, es)},.T.)`), ...rings.map((r) => E(`FACE_BOUND('',${sLoop(E, [[r[zq], zq === 0]])},.T.)`))].join(',')}),${E(`PLANE('',${sAX(E, [0, 0, zq], [0, 0, 1], [1, 0, 0])})`)},${sense})`);
      }
    }));
    const dp = dimensions(plate), vol = Wd * Wd - N * N * A48;
    ok(plate.kind === 'mesh' && plate.stats.tessellatedFaces === 6 + N * N && dp.closed && dp.consistentNormals, `STEP plate with ${N * N} hole(s): ${6 + N * N} faces, ${plate.count} triangles, closed with consistent normals (${dp.closed}, ${dp.consistentNormals})`);
    ok(Math.abs(dp.volume - vol) < 1e-9 && Math.abs(dp.area - (2 * vol + 4 * Wd + N * N * P48)) < 1e-9, `STEP plate with ${N * N} hole(s): volume ${dp.volume} (exact ${Wd * Wd - N * N * Math.PI}), area ${dp.area}`);
  }
  const broken = await imp('broken.stp', stepFile((E) => { sUnits(E); sQuadFaces(E, V, Q); E("PRODUCT('x',(1,2"); }));
  ok(broken.kind === 'mesh' && broken.count === 12, 'STEP with a malformed, unrelated entity still reads the solid');
  // AP242 tessellated geometry
  const cl = `COORDINATES_LIST('',8,(${V.map((v) => `(${v.map(RN).join(',')})`).join(',')}))`;
  isCube(await imp('tess.stp', stepFile((E) => { sUnits(E); E(`TRIANGULATED_FACE('',${E(cl)},8,((0.,0.,1.)),$,(),(${T.map((t) => `(${t.map((i) => i + 1).join(',')})`).join(',')}))`); }, 'AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }')), 'STEP AP242 TRIANGULATED_FACE');
  isCube(await imp('tess2.stp', stepFile((E) => {
    const strips = Q.slice(0, 3).map((q) => `(${[q[0], q[1], q[3], q[2]].map((i) => 8 - i).join(',')})`), fans = Q.slice(3).map((q) => `(${q.map((i) => 8 - i).join(',')})`);
    E(`COMPLEX_TRIANGULATED_FACE('',${E(cl)},8,((0.,0.,1.)),$,(8,7,6,5,4,3,2,1),(${strips.join(',')}),(${fans.join(',')}))`);
  })), 'STEP AP242 COMPLEX_TRIANGULATED_FACE (strips, fans, pnindex)');
  const fb = await imp('faceted.stp', stepFile((E) => { const ps = V.map((v) => sP(E, v)); const fs = Q.map((q) => E(`FACE('',(${E(`FACE_OUTER_BOUND('',${E(`POLY_LOOP('',(${q.map((i) => ps[i]).join(',')}))`)},.T.)`)}))`)); E(`FACETED_BREP('',${E(`CLOSED_SHELL('',(${fs.join(',')}))`)})`); }));
  isCube(fb, 'STEP FACETED_BREP / POLY_LOOP');
  ok(fb.warnings.some((w) => /length unit/i.test(w)), 'STEP without units says so');
  // IGES: two lines, a half-circle arc and a quadratic NURBS curve
  const ig = await imp('wire.igs', igesFile([{ type: 110, params: '110,0.,0.,0.,10.,0.,0.' }, { type: 110, params: '110,10.,0.,0.,10.,-1.,0.' }, { type: 100, params: '100,0.,5.,0.,10.,0.,0.,0.' }, { type: 126, params: '126,2,2,1,0,1,0,0.,0.,0.,1.,1.,1.,1.,1.,1.,0.,0.,0.,5.,-5.,0.,10.,0.,0.,0.,1.,0.,0.,1.' }, { type: 116, params: '116,5.,2.,0.,0' }, { type: 144, params: '144,1,0,0,0' }]));
  ok(ig.kind === 'polylines' && ig.polylines.length === 4 && ig.stats.units === 'mm' && ig.stats.droppedAxis === 'z' && ig.stats.surfaces === 1 && ig.stats.points === 1, `IGES: 4 curves, mm, z dropped (${ig.polylines.length}, ${ig.stats.units}, ${ig.stats.droppedAxis})`);
  ok(Math.abs(ig.bbox.min[0]) < 1e-9 && Math.abs(ig.bbox.max[0] - 10) < 1e-9 && Math.abs(ig.bbox.max[1] - 5) < 1e-9 && Math.abs(ig.bbox.min[1] + 2.5) < 1e-9, `IGES bbox: arc top y = 5, NURBS dip y = −2.5 (${ig.bbox.min}, ${ig.bbox.max})`);
  ok(ig.polylines[2].x.length === 25 && Math.abs(ig.polylines[2].x[12] - 5) < 1e-9 && Math.abs(ig.polylines[2].y[12] - 5) < 1e-9, 'IGES arc discretised counter-clockwise through (5, 5)');
  ok(ig.warnings.some((w) => /not tessellated/.test(w)) && ig.warnings.some((w) => /projected/.test(w)) && formatOf('a.igs').support === 'partial', 'IGES warns that surfaces are not tessellated');
  const ig2 = await imp('moved.iges', igesFile([{ type: 124, params: '124,1.,0.,0.,100.,0.,1.,0.,200.,0.,0.,1.,0.' }, { type: 110, params: '110,0.,0.,0.,1.,1.,0.', xf: 1 }, { type: 112, params: '112,1,1,3,1,0.,2.,0.,1.,0.,0.,0.,0.,1.,0.,0.,0.,0.,0.' }]));
  ok(ig2.polylines.length === 2 && ig2.bbox.max[0] === 101 && ig2.bbox.max[1] === 201, 'IGES transformation matrix (124) applied; parametric spline (112) evaluated');
});

await group('2-D drawings', async () => {
  const dxf = (ents) => `0\nSECTION\n2\nENTITIES\n${ents}0\nENDSEC\n0\nEOF\n`;
  const d1 = await imp('sq.dxf', dxf('0\nLWPOLYLINE\n8\n0\n90\n4\n70\n1\n10\n0\n20\n0\n10\n2\n20\n0\n10\n2\n20\n2\n10\n0\n20\n2\n0\nARC\n8\n0\n10\n5\n20\n0\n40\n1\n50\n0\n51\n180\n0\nELLIPSE\n10\n10\n20\n0\n11\n2\n21\n0\n40\n0.5\n41\n0\n42\n6.283185307179586\n'));
  ok(d1.kind === 'polylines' && d1.polylines.length === 3 && d1.polylines[0].closed && d1.polylines[0].x.length === 4, 'DXF: LWPOLYLINE square + ARC + ELLIPSE → 3 polylines');
  const arc = d1.polylines[1], ell = d1.polylines[2];
  ok(arc.x.length === 25 && Math.abs(arc.x[0] - 6) < 1e-9 && Math.abs(arc.y[12] - 1) < 1e-9 && Math.abs(arc.x[24] - 4) < 1e-9 && !arc.closed, 'DXF ARC from 0° to 180°');
  ok(ell.closed && Math.abs(Math.max(...ell.x) - 12) < 1e-9 && Math.abs(Math.max(...ell.y) - 1) < 1e-9, 'DXF ELLIPSE (major 2, ratio 0.5)');
  const d2 = await imp('face.dxf', dxf('0\n3DFACE\n10\n0\n20\n0\n30\n0\n11\n1\n21\n0\n31\n0\n12\n1\n22\n1\n32\n0\n13\n0\n23\n1\n33\n0\n'));
  ok(d2.kind === 'mesh' && d2.count === 2 && Math.abs(dimensions(d2).area - 1) < 1e-12 && !dimensions(d2).closed, 'DXF 3DFACE → mesh (open surface, area 1)');
  const d3 = await imp('spl.dxf', dxf('0\nSPLINE\n70\n8\n71\n2\n72\n6\n73\n3\n40\n0\n40\n0\n40\n0\n40\n1\n40\n1\n40\n1\n10\n0\n20\n0\n30\n0\n10\n5\n20\n10\n30\n0\n10\n10\n20\n0\n30\n0\n'));
  ok(d3.polylines.length === 1 && Math.abs(Math.max(...d3.polylines[0].y) - 5) < 1e-9, 'DXF SPLINE evaluated (quadratic peak at y = 5)');
  const dInf = await imp('odd.dxf', dxf('0\nARC\n10\n0\n20\n0\n40\n1\n50\n0\n51\n-1e308\n0\nARC\n10\n0\n20\n0\n40\n1\n50\n-Infinity\n51\n90\n0\nLINE\n10\n0\n20\n0\n11\n1\n21\n1\n'));
  ok(dInf.kind === 'polylines' && dInf.polylines.length >= 1 && dInf.polylines.every((p) => p.x.length <= 200), 'DXF arcs with absurd angles neither hang nor explode');
  await throws('DXF without supported entities', () => imp('e.dxf', dxf('0\nTEXT\n1\nhello\n')), /No LINE/);
  const svg = await imp('shapes.svg', `<?xml version="1.0"?><svg xmlns="http://www.w3.org/2000/svg" width="200" height="100"><defs><rect x="0" y="0" width="999" height="999"/></defs><rect x="10" y="10" width="10" height="10"/><g transform="translate(100,0) scale(2)"><path d="M0 0 h10 v10 h-10 z"/></g><circle cx="50" cy="50" r="5"/><polygon points="0,90 10,90 10,100"/><path d="M 60 0 A 5 5 0 0 1 70 0 M 80,20 C 80,10 90,10 90,20 q 5,5 10,0 L 100 30"/><line x1="0" y1="0" x2="1" y2="1" transform="matrix(1 0 0 1 150 50)"/></svg>`);
  ok(svg.kind === 'polylines' && svg.polylines.length === 7, `SVG: rect, group path, circle, polygon, 2 sub-paths, line → 7 polylines (${svg.polylines.length})`);
  const [r, gp, c] = svg.polylines, ar = (p) => { let s = 0; for (let i = 0; i < p.x.length; i++) { const j = (i + 1) % p.x.length; s += p.x[i] * p.y[j] - p.x[j] * p.y[i]; } return Math.abs(s) / 2; };
  ok(r.closed && r.x.length === 4 && ar(r) === 100 && gp.closed && Math.min(...gp.x) === 100 && Math.max(...gp.x) === 120 && ar(gp) === 400, 'SVG rect area 100; translate + scale group gives x 100–120, area 400');
  ok(c.closed && Math.abs(ar(c) - 25 * A48) < 1e-9 && svg.bbox.min[1] === -5 && svg.bbox.max[1] === 100, 'SVG circle flattened; elliptical arc bulges 5 units beyond y = 0');
  ok(Math.max(...r.y) === 85 && Math.min(...r.y) === 75, 'SVG y axis flipped within the bounding box');
  const hp = await imp('sq.hpgl', 'IN;SP1;PU0,0;PD4000,0,4000,4000,0,4000,0,0;PU;PA8000,0;CI400;PU;');
  ok(hp.kind === 'polylines' && hp.polylines.length === 2 && hp.polylines[0].closed && hp.polylines[0].x.length === 4 && Math.abs(dimensions(hp).area - (10000 + 100 * A48)) < 1e-6 && hp.stats.units === 'mm', 'HPGL: 100 mm square + circle of radius 10 mm');
  const hpA = await imp('arc.hpgl', 'IN;PU10,0;PD;AA0,0,1e15;PU;');
  ok(hpA.polylines[0].x.length <= 200, 'HPGL arc with an absurd sweep stays bounded');
  const hp2 = await imp('sq.plt', 'IN;PU10,10;PR;PD100,0,0,100,-100,0,0,-100;PU;');
  ok(hp2.format === 'HPGL plot file' && hp2.pathway === 'drawing' && hp2.polylines[0].closed && Math.abs(hp2.bbox.max[0] - 2.75) < 1e-12, 'HPGL sniffed from .plt, relative moves');
  const xy = await imp('profile.xy', '# x y\n0 0\n2 0\n2 1\n0 1\n0 0\n\n5 5\n6 6\n');
  ok(xy.kind === 'polylines' && xy.polylines.length === 2 && xy.polylines[0].closed && xy.polylines[0].x.length === 4 && !xy.polylines[1].closed, 'x-y file: blank line separates polylines, repeated first point closes');
});

// ---- raster writers ----------------------------------------------------------------------------------------
function lzwEncode(data) {
  const out = [];
  let buf = 0, cnt = 0, width = 9, next = 258, table = new Map(), w = '';
  const put = (code) => { buf = (buf << width) | code; cnt += width; while (cnt >= 8) { out.push((buf >> (cnt - 8)) & 255); cnt -= 8; } buf &= (1 << cnt) - 1; };
  const codeOf = (s) => (s.length === 1 ? s.charCodeAt(0) : table.get(s));
  put(256);
  for (const b of data) {
    const c = String.fromCharCode(b);
    if (w === '') { w = c; continue; }
    if (table.has(w + c)) { w += c; continue; }
    put(codeOf(w));
    table.set(w + c, next++);
    if (next === 4094) { put(256); table = new Map(); next = 258; width = 9; } else if (next > (1 << width) - 1) width++;
    w = c;
  }
  if (w !== '') put(codeOf(w));
  put(257);
  if (cnt) out.push((buf << (8 - cnt)) & 255);
  return new Uint8Array(out);
}
function tiffFile(pages, little = true) {
  const chunks = [new Uint8Array(8)], ifds = [];
  let pos = 8;
  const put = (b) => { const o = pos; chunks.push(b); pos += b.length; if (pos % 2) { chunks.push(new Uint8Array(1)); pos++; } return o; };
  for (const pg of pages) {
    const offs = pg.strips.map((s) => put(s)), tags = [...pg.tags, [pg.tiled ? 324 : 273, 4, offs], [pg.tiled ? 325 : 279, 4, pg.strips.map((s) => s.length)]].sort((a, b) => a[0] - b[0]);
    ifds.push(tags.map(([tag, type, vals]) => {
      const bytes = type === 2 ? enc.encode(vals + '\0') : le({ 3: 'u16', 4: 'u32', 12: 'f64' }[type], vals, little), count = type === 2 ? bytes.length : vals.length;
      let vf = new Uint8Array(4);
      if (bytes.length <= 4) vf.set(bytes); else vf = le('u32', [put(bytes)], little);
      return cat(le('u16', [tag, type], little), le('u32', [count], little), vf);
    }));
  }
  const at = [];
  let p = pos;
  for (const e of ifds) { at.push(p); p += 2 + 12 * e.length + 4; }
  ifds.forEach((e, k) => chunks.push(le('u16', [e.length], little), ...e, le('u32', [k + 1 < ifds.length ? at[k + 1] : 0], little)));
  const out = cat(...chunks);
  out.set(cat(little ? 'II' : 'MM', le('u16', [42], little), le('u32', [at[0]], little)), 0);
  return out;
}
/** One TIFF page: values row-major, type 'u8' | 'u16' | 'i16' | 'f32' | 'f64'. */
function tiffPage({ w, h, values, type = 'f32', comp = 1, pred = 1, tile = 0, little = true, geo = null, nodata, rowsPerStrip }) {
  const bp = { u8: 1, u16: 2, i16: 2, f32: 4, f64: 8 }[type], sf = type[0] === 'f' ? 3 : type[0] === 'i' ? 2 : 1;
  const encodeRows = (vals, cw, rows) => {
    let v = Array.from(vals);
    if (pred === 2) for (let r = 0; r < rows; r++) for (let i = cw - 1; i > 0; i--) v[r * cw + i] -= v[r * cw + i - 1];
    let b = le(type, pred === 2 && type === 'u16' ? v.map((x) => (x + 65536) % 65536) : pred === 2 && type === 'u8' ? v.map((x) => (x + 256) % 256) : v, little);
    if (pred === 3) {
      const src = le(type, v, false), o = new Uint8Array(src.length);
      for (let r = 0; r < rows; r++) { for (let s = 0; s < cw; s++) for (let k = 0; k < bp; k++) o[r * cw * bp + k * cw + s] = src[(r * cw + s) * bp + k]; for (let i = cw * bp - 1; i > 0; i--) o[r * cw * bp + i] = (o[r * cw * bp + i] - o[r * cw * bp + i - 1]) & 255; }
      b = o;
    }
    return comp === 5 ? lzwEncode(b) : comp === 8 ? new Uint8Array(zlib.deflateSync(b)) : comp === 32773 ? cat(...Array.from({ length: Math.ceil(b.length / 100) }, (_, k) => { const s = b.subarray(k * 100, k * 100 + 100); return cat(le('u8', [s.length - 1]), s); })) : b;
  };
  const strips = [], tags = [[256, 4, [w]], [257, 4, [h]], [258, 3, [bp * 8]], [259, 3, [comp]], [262, 3, [1]], [277, 3, [1]], [284, 3, [1]], [317, 3, [pred]], [339, 3, [sf]]];
  if (tile) {
    for (let ty = 0; ty < h; ty += tile) for (let tx = 0; tx < w; tx += tile) { const t = new Array(tile * tile).fill(0); for (let r = 0; r < tile && ty + r < h; r++) for (let c = 0; c < tile && tx + c < w; c++) t[r * tile + c] = values[(ty + r) * w + tx + c]; strips.push(encodeRows(t, tile, tile)); }
    tags.push([322, 3, [tile]], [323, 3, [tile]]);
  } else {
    const rps = rowsPerStrip || h;
    for (let r = 0; r < h; r += rps) strips.push(encodeRows(values.slice(r * w, Math.min(h, r + rps) * w), w, Math.min(rps, h - r)));
    tags.push([278, 4, [rps]]);
  }
  if (geo) tags.push([33550, 12, [geo.dx, geo.dy, 0]], [33922, 12, [0, 0, 0, geo.x0, geo.y1, 0]], [34735, 3, [1, 1, 0, 2, 1024, 0, 1, geo.geographic ? 2 : 1, 1025, 0, 1, 1]]);
  if (nodata !== undefined) tags.push([42113, 2, String(nodata)]);
  return { tags, strips, tiled: !!tile };
}
function netcdf(dims, vars, version = 1) {
  const pad4 = (b) => cat(b, new Uint8Array((4 - (b.length % 4)) % 4)), nm = (s) => cat(be('u32', [s.length]), pad4(enc.encode(s)));
  const NC = { i8: 1, char: 2, i16: 3, i32: 4, f32: 5, f64: 6 }, attList = (atts) => (atts && atts.length ? cat(be('u32', [12, atts.length]), ...atts.map(([n, t, v]) => cat(nm(n), be('u32', [NC[t], t === 'char' ? v.length : v.length]), pad4(t === 'char' ? enc.encode(v) : be(t, v))))) : be('u32', [0, 0]));
  const datas = vars.map((v) => pad4(be(v.type, v.data)));
  const header = (begins) => cat('CDF', le('u8', [version]), be('u32', [0]), be('u32', [10, dims.length]), ...dims.map(([n, l]) => cat(nm(n), be('u32', [l]))), attList([['title', 'char', 'test']]), be('u32', [11, vars.length]), ...vars.map((v, k) => cat(nm(v.name), be('u32', [v.dims.length, ...v.dims]), attList(v.atts), be('u32', [NC[v.type], datas[k].length]), version === 1 ? be('u32', [begins[k]]) : be('u64', [begins[k]]))));
  const hl = header(vars.map(() => 0)).length, begins = [];
  let p = hl;
  for (const d of datas) { begins.push(p); p += d.length; }
  return cat(header(begins), ...datas);
}

await group('GIS vectors', async () => {
  const gj = await imp('coast.geojson', JSON.stringify({ type: 'FeatureCollection', features: [{ type: 'Feature', properties: {}, geometry: { type: 'Polygon', coordinates: [[[10, 50], [11, 50], [11, 51], [10, 51], [10, 50]]] } }, { type: 'Feature', geometry: { type: 'LineString', coordinates: [[10, 50], [12, 52]] } }, { type: 'Feature', geometry: { type: 'Point', coordinates: [10.5, 50.5] } }] }));
  ok(gj.kind === 'polylines' && gj.polylines.length === 2 && gj.polylines[0].closed && gj.polylines[0].x.length === 4 && gj.geographic && gj.format === 'GeoJSON' && gj.markers.length === 1, 'GeoJSON polygon + line + point');
  const gp = await imp('pts.json', JSON.stringify({ type: 'MultiPoint', coordinates: [[1, 2, 3], [4, 5, 6]] }));
  ok(gp.kind === 'points' && gp.count === 2 && gp.points[5] === 6 && gp.format === 'GeoJSON', 'GeoJSON points only (sniffed from .json)');
  const kml = `<?xml version="1.0" encoding="UTF-8"?><kml xmlns="http://www.opengis.net/kml/2.2"><Document><Placemark><name>Area &amp; outfall</name><Polygon><outerBoundaryIs><LinearRing><coordinates>\n 54.0,24.0,0 54.5,24.0,0 54.5,24.4,0 54.0,24.4,0 54.0,24.0,0\n</coordinates></LinearRing></outerBoundaryIs></Polygon></Placemark><Placemark><LineString><coordinates>54.1,24.1 54.2,24.3</coordinates></LineString></Placemark><Placemark><Point><coordinates>54.25,24.2,-12</coordinates></Point></Placemark></Document></kml>`;
  const k1 = await imp('site.kml', kml);
  ok(k1.kind === 'polylines' && k1.polylines.length === 2 && k1.polylines[0].closed && k1.polylines[0].x.length === 4 && k1.bbox.max[0] === 54.5 && k1.bbox.min[1] === 24 && k1.geographic && k1.stats.points === 1, 'KML polygon ring, line string and point');
  const k2 = await imp('site.kmz', zip([['doc.kml', kml]], true));
  ok(k2.kind === 'polylines' && k2.polylines.length === 2 && k2.format === 'KML / KMZ', 'KMZ (deflated zip)');
  const gpx = await imp('track.gpx', '<?xml version="1.0"?><gpx version="1.1"><wpt lat="24.5" lon="54.5"><ele>3</ele></wpt><trk><trkseg><trkpt lat="24.0" lon="54.0"><ele>1</ele></trkpt><trkpt lat="24.1" lon="54.2"/><trkpt lat="24.3" lon="54.3"/></trkseg></trk><rte><rtept lat="1" lon="2"/><rtept lat="3" lon="4"/></rte></gpx>');
  ok(gpx.kind === 'polylines' && gpx.polylines.length === 2 && gpx.polylines[0].x[1] === 54.2 && gpx.polylines[0].y[2] === 24.3 && gpx.markers[0][0] === 54.5, 'GPX track, route and waypoint (x = lon, y = lat)');
  const gml = await imp('zone.gml', '<gml:FeatureCollection xmlns:gml="http://www.opengis.net/gml"><gml:featureMember><gml:Polygon><gml:exterior><gml:LinearRing><gml:posList srsDimension="2">0 0 4 0 4 3 0 3 0 0</gml:posList></gml:LinearRing></gml:exterior></gml:Polygon></gml:featureMember><gml:LineString><gml:coordinates>1,1 2,2 3,1</gml:coordinates></gml:LineString><gml:Point><gml:pos>9 9</gml:pos></gml:Point></gml:FeatureCollection>');
  ok(gml.kind === 'polylines' && gml.polylines.length === 2 && gml.polylines[0].closed && Math.abs(dimensions(gml).area - 12) < 1e-12 && gml.polylines[1].x.length === 3 && gml.markers.length === 1, 'GML posList ring, coordinates line string and pos point');
  // shapefile: one polygon record (ring of 5 points) and one polyline with two parts
  const ring = [[0, 0], [0, 2], [3, 2], [3, 0], [0, 0]], rec = (n, content) => cat(be('i32', [n, content.length / 2]), content);
  const polyRec = (type, parts, pts) => cat(le('i32', [type]), le('f64', [0, 0, 3, 2]), le('i32', [parts.length, pts.length, ...parts]), le('f64', pts.flat()));
  const shp = (type, recs) => { const body = cat(...recs), h = new Uint8Array(100), dv = new DataView(h.buffer); dv.setInt32(0, 9994, false); dv.setInt32(24, (100 + body.length) / 2, false); dv.setInt32(28, 1000, true); dv.setInt32(32, type, true); return cat(h, body); };
  const dbf = cat(le('u8', [3, 124, 1, 1]), le('u32', [1]), le('u16', [65, 11]), new Uint8Array(20), cat('NAME', new Uint8Array(7), 'C', new Uint8Array(4), le('u8', [10]), new Uint8Array(15)), le('u8', [0x0d]), ' intake    ');
  const s1 = await imp('zone.shp', shp(5, [rec(1, polyRec(5, [0], ring)), rec(2, polyRec(3, [0, 2], [[5, 5], [6, 6], [7, 7], [8, 9]]))]), { companion: { 'zone.dbf': F('zone.dbf', dbf), 'zone.prj': F('zone.prj', 'GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984"]]') } });
  ok(s1.kind === 'polylines' && s1.polylines.length === 3 && s1.polylines[0].closed && s1.polylines[0].x.length === 4 && !s1.polylines[1].closed && s1.polylines[2].y[1] === 9, 'Shapefile polygon ring + two-part polyline');
  ok(s1.stats.records === 2 && s1.stats.attributes[0] === 'NAME' && s1.stats.attributeRecords === 1 && s1.geographic === true && /GEOGCS/.test(s1.stats.projection), 'Shapefile .dbf field names and .prj read from companions');
  const s2 = await imp('wells.shp', shp(11, [rec(1, cat(le('i32', [11]), le('f64', [500000, 2700000, -15, 0]))), rec(2, cat(le('i32', [11]), le('f64', [500010, 2700020, -18, 0])))]));
  ok(s2.kind === 'points' && s2.count === 2 && s2.points[2] === -15 && s2.points[4] === 2700020 && s2.geographic === false, 'Shapefile PointZ, lone .shp');
  const mif = await imp('zone.mif', 'Version 300\nCharset "WindowsLatin1"\nDelimiter ","\nColumns 1\n  Name Char(10)\nData\n\nRegion 1\n  5\n0 0\n4 0\n4 4\n0 4\n0 0\n    Pen (1,2,0)\n    Brush (2,16777215,16777215)\nPline 3\n1 1\n2 2\n3 1\nLine 0 0 9 9\nPoint 5 5\n');
  ok(mif.kind === 'polylines' && mif.polylines.length === 3 && mif.polylines[0].closed && mif.polylines[0].x.length === 4 && mif.polylines[1].x.length === 3 && mif.bbox.max[0] === 9, 'MapInfo MIF region, pline and line');
});

await group('rasters and DEMs', async () => {
  // 5 × 4 DEM, value = 100·row(from the top) + col; lower-left corner (1000, 2000), cell 10
  const W = 5, H = 4, dem = Array.from({ length: W * H }, (_, i) => 100 * Math.floor(i / W) + (i % W) + 0.5);
  const isDem = (g, label, { nod = 0, x0 = 1005, y0 = 2005, d = 10 } = {}) => {
    ok(g.kind === 'grid' && g.grid.x.length === W && g.grid.y.length === H && g.grid.z.length === H && g.grid.z[0].length === W, `${label}: 5 × 4 grid (kind ${g.kind})`);
    ok(Math.abs(g.grid.x[0] - x0) < 1e-9 && Math.abs(g.grid.x[4] - (x0 + 4 * d)) < 1e-9 && Math.abs(g.grid.y[0] - y0) < 1e-9 && Math.abs(g.grid.y[3] - (y0 + 3 * d)) < 1e-9 && g.grid.y[1] > g.grid.y[0], `${label}: cell-centre coordinates ascending (${g.grid.x[0]}, ${g.grid.y[0]})`);
    ok(g.grid.z[0][1] === 301.5 && g.grid.z[3][4] === 4.5 && (nod ? Number.isNaN(g.grid.z[3][0]) : g.grid.z[3][0] === 0.5) && g.grid.nodata === nod, `${label}: values (top row of the file is the highest y), ${g.grid.nodata} nodata`);
  };
  const ascBody = (vals) => Array.from({ length: H }, (_, r) => vals.slice(r * W, r * W + W).join(' ')).join('\n');
  isDem(await imp('dem.asc', `ncols 5\nnrows 4\nxllcorner 1000\nyllcorner 2000\ncellsize 10\nNODATA_value -9999\n${ascBody(dem)}\n`), 'ESRI ASCII grid');
  const withHole = dem.slice(); withHole[0] = -9999;
  isDem(await imp('dem.asc', `NCOLS 5\nNROWS 4\nXLLCENTER 1005\nYLLCENTER 2005\nCELLSIZE 10\nNODATA_VALUE -9999\n${ascBody(withHole)}\n`), 'ESRI ASCII grid (xllcenter, nodata)', { nod: 1 });
  const geo = { dx: 10, dy: 10, x0: 1000, y1: 2040 };
  const t1 = await imp('dem.tif', tiffFile([tiffPage({ w: W, h: H, values: dem, geo })]));
  isDem(t1, 'GeoTIFF float32 uncompressed');
  ok(t1.format === 'GeoTIFF / TIFF' && t1.pathway === 'gis' && t1.grid.geographic === false, 'GeoTIFF format / pathway / projected');
  isDem(await imp('dem.tif', tiffFile([tiffPage({ w: W, h: H, values: dem, geo, comp: 8 })])), 'GeoTIFF Deflate');
  isDem(await imp('dem.tif', tiffFile([tiffPage({ w: W, h: H, values: dem, geo, comp: 5 })])), 'GeoTIFF LZW');
  isDem(await imp('dem.tif', tiffFile([tiffPage({ w: W, h: H, values: dem, geo, comp: 32773, rowsPerStrip: 1 })])), 'GeoTIFF PackBits, one row per strip');
  isDem(await imp('dem.tiff', tiffFile([tiffPage({ w: W, h: H, values: dem, geo, comp: 8, pred: 3, little: false })], false)), 'GeoTIFF big-endian, Deflate, floating-point predictor');
  isDem(await imp('dem.tif', tiffFile([tiffPage({ w: W, h: H, values: withHole, geo, comp: 5, pred: 3, type: 'f64', tile: 16, nodata: -9999 })])), 'GeoTIFF float64, tiled, LZW + predictor 3, GDAL_NODATA', { nod: 1 });
  const i16 = dem.map((v) => Math.round(v * 10)), ti = await imp('dem16.tif', tiffFile([tiffPage({ w: W, h: H, values: i16, type: 'i16', geo: { ...geo, geographic: true }, comp: 8, pred: 2 })]));
  ok(ti.kind === 'grid' && ti.grid.z[0][1] === 3015 && ti.grid.z[3][4] === 45 && ti.grid.geographic === true, 'GeoTIFF int16 with horizontal predictor, geographic key');
  // LZW across the 9 → 10 → 11-bit code-width changes, and a multi-page stack
  const big = Array.from({ length: 64 * 64 }, (_, i) => (i * 7919 + ((i * i) % 251)) % 256 > 127 ? 255 : (i * 31) % 97);
  const pg = (vals, comp) => tiffPage({ w: 64, h: 64, values: vals, type: 'u8', comp });
  const st = await imp('stack.tif', tiffFile([pg(big, 5), pg(big.map((v) => 255 - v), 8), pg(big, 32773)]), { threshold: 127 });
  const solid0 = big.filter((v) => v > 127).length;
  ok(st.kind === 'voxels' && st.voxels.nz === 3 && st.voxels.nx === 64 && st.format === 'TIFF stack' && st.pathway === 'voxel', 'multi-page 8-bit TIFF → 3-D voxel stack');
  let s0 = 0, s1 = 0, same = 0;
  for (let i = 0; i < 4096; i++) { s0 += st.voxels.data[i]; s1 += st.voxels.data[4096 + i]; same += st.voxels.data[i] === st.voxels.data[8192 + i] ? 1 : 0; }
  ok(s0 === solid0 && s1 === 4096 - big.filter((v) => 255 - v <= 127).length && same === 4096 && st.voxels.data[63 * 64] === (big[0] > 127 ? 1 : 0), `TIFF stack: LZW page decodes exactly across code-width changes (${s0} = ${solid0} solid), rows flipped to y-up`);
  const im = await imp('mask.tif', tiffFile([pg(big, 1)]));
  ok(im.kind === 'voxels' && im.voxels.nz === 1 && im.format === 'TIFF image' && im.stats.method === 'otsu', 'single 8-bit TIFF without georeferencing → 2-D voxels (Otsu)');
  await throws('BigTIFF', () => imp('big.tif', cat('II', le('u16', [43, 8, 0]), new Uint8Array(16))), /BigTIFF/);
  await throws('JPEG-compressed TIFF', () => imp('j.tif', tiffFile([{ ...tiffPage({ w: 2, h: 2, values: [1, 2, 3, 4], type: 'u8' }), tags: tiffPage({ w: 2, h: 2, values: [1, 2, 3, 4], type: 'u8' }).tags.map((t) => (t[0] === 259 ? [259, 3, [7]] : t)) }])), /JPEG/);
  // NetCDF classic: y descending in the file, packed shorts with scale / offset / fill
  const ys = [2035, 2025, 2015, 2005], xs = [1005, 1015, 1025, 1035, 1045];
  const nc1 = await imp('dem.nc', netcdf([['y', 4], ['x', 5]], [{ name: 'x', dims: [1], type: 'f64', data: xs, atts: [['units', 'char', 'm']] }, { name: 'y', dims: [0], type: 'f64', data: ys }, { name: 'elevation', dims: [0, 1], type: 'f32', data: dem, atts: [['units', 'char', 'm']] }]));
  isDem(nc1, 'NetCDF classic (CDF-1) float variable with coordinate variables');
  ok(nc1.stats.variable === 'elevation' && nc1.stats.units === 'm' && nc1.format === 'NetCDF classic', 'NetCDF stats (variable, units)');
  const packed = withHole.map((v) => (v === -9999 ? -32768 : Math.round((v - 100) / 0.5)));
  isDem(await imp('dem.grd', netcdf([['lat', 4], ['lon', 5], ['t', 1]], [{ name: 'lon', dims: [1], type: 'f32', data: xs }, { name: 'lat', dims: [0], type: 'f32', data: ys }, { name: 'z', dims: [2, 0, 1], type: 'i16', data: packed, atts: [['scale_factor', 'f64', [0.5]], ['add_offset', 'f64', [100]], ['_FillValue', 'i16', [-32768]]] }], 2)), 'NetCDF CDF-2 (.grd), 3-D packed shorts with scale_factor / add_offset / _FillValue', { nod: 1 });
  await throws('truncated NetCDF-4 / HDF5 in a .nc file', () => imp('v4.nc', cat(le('u8', [0x89, 0x48, 0x44, 0x46, 0x0d, 0x0a, 0x1a, 0x0a]), new Uint8Array(64))), /HDF5 file is truncated or corrupt/);
  // Surfer grids (rows from the lowest y upwards)
  const up = Array.from({ length: H }, (_, r) => dem.slice((H - 1 - r) * W, (H - r) * W)).flat();
  isDem(await imp('dem.grd', `DSAA\n5 4\n1005 1045\n2005 2035\n0.5 304.5\n${ascBody(up)}\n`), 'Surfer ASCII grid (DSAA)');
  isDem(await imp('dem.grd', cat('DSBB', le('i16', [5, 4]), le('f64', [1005, 1045, 2005, 2035, 0.5, 304.5]), le('f32', up))), 'Surfer 6 binary grid (DSBB)');
  isDem(await imp('dem.grd', cat('DSRB', le('i32', [4, 1]), 'GRID', le('i32', [72, 4, 5]), le('f64', [1005, 2005, 10, 10, 0.5, 304.5, 0, 1.70141e38]), 'DATA', le('i32', [160]), le('f64', up))), 'Surfer 7 binary grid (DSRB)');
  // band-interleaved rasters with header companions
  const bil = cat(...Array.from({ length: H }, (_, r) => cat(le('f32', dem.slice(r * W, r * W + W)), le('f32', new Array(W).fill(7)))));
  isDem(await imp('dem.bil', bil, { companion: { 'dem.hdr': F('dem.hdr', 'BYTEORDER I\nLAYOUT BIL\nNROWS 4\nNCOLS 5\nNBANDS 2\nNBITS 32\nPIXELTYPE FLOAT\nULXMAP 1005\nULYMAP 2035\nXDIM 10\nYDIM 10\n') } }), 'ESRI BIL (2 bands) + .hdr');
  isDem(await imp('dem.bsq', cat(new Uint8Array(8), be('f64', dem)), { companion: { 'dem.hdr': F('dem.hdr', 'ENVI\nsamples = 5\nlines = 4\nbands = 1\nheader offset = 8\ndata type = 5\ninterleave = bsq\nbyte order = 1\nmap info = {UTM, 1, 1, 1000, 2040, 10, 10, 40, North}\n') } }), 'ENVI BSQ float64 big-endian + .hdr');
  await throws('BIL without its header', () => imp('dem.bil', bil), /\.hdr/);
  // soundings
  const rows = []; for (let j = 0; j < 4; j++) for (let i = 0; i < 5; i++) rows.push(`${54 + i * 0.01},${24 + j * 0.01},${10 + i + j}`);
  const snd = await imp('soundings.csv', 'Longitude,Latitude,Depth (m)\n' + rows.join('\n') + '\n');
  ok(snd.kind === 'grid' && snd.grid.x.length === 5 && snd.grid.y.length === 4 && snd.grid.z[0][0] === -10 && snd.grid.z[3][4] === -17 && snd.grid.geographic && snd.warnings.some((w) => /positive down/.test(w)), 'lon, lat, depth CSV on a lattice → grid with depths as negative elevations');
  let seed = 7; const rnd = () => ((seed = (seed * 16807) % 2147483647) / 2147483647);
  const scat = Array.from({ length: 600 }, () => { const x = rnd() * 100, y = rnd() * 50; return `${x.toFixed(3)} ${y.toFixed(3)} ${(5 + 0.1 * x).toFixed(3)}`; });
  const sc = await imp('survey.xyz', 'easting northing depth\n' + scat.join('\n') + '\n');
  ok(sc.kind === 'grid' && sc.count === 600 && sc.points.length === 1800 && sc.warnings.some((w) => /binned/.test(w)) && sc.grid.z.flat().every((v) => v <= -4.9 && v >= -15.1), `scattered easting / northing / depth → binned grid ${sc.stats.nx} × ${sc.stats.ny}, points kept`);
});

await group('point clouds', async () => {
  const hd = new Uint8Array(227), dv = new DataView(hd.buffer);
  hd.set(enc.encode('LASF')); hd[24] = 1; hd[25] = 2; dv.setUint16(94, 227, true); dv.setUint32(96, 227, true); hd[104] = 1; dv.setUint16(105, 28, true); dv.setUint32(107, 10, true);
  [0.01, 0.01, 0.001].forEach((v, k) => dv.setFloat64(131 + 8 * k, v, true)); [100, 200, 0].forEach((v, k) => dv.setFloat64(155 + 8 * k, v, true));
  const recs = cat(...Array.from({ length: 10 }, (_, i) => cat(le('i32', [i * 100, i * 50, -i * 1000]), new Uint8Array(16))));
  const las = await imp('scan.las', cat(hd, recs));
  ok(las.kind === 'points' && las.count === 10 && las.points[0] === 100 && las.points[1] === 200 && Math.abs(las.points[27] - 109) < 1e-9 && Math.abs(las.points[28] - 204.5) < 1e-9 && Math.abs(las.points[29] + 9) < 1e-9 && las.stats.pointFormat === 1, 'LAS 1.2 format 1: 10 points scaled and offset');
  const cutLas = await imp('cut.las', cat(hd, recs.subarray(0, 28 * 6 + 5)));
  ok(cutLas.count === 6 && cutLas.warnings.some((w) => /truncated/.test(w)), 'truncated LAS keeps the complete records and warns');
  hd[104] = 0x80 | 1;
  await throws('LAZ-compressed point format inside .las', () => imp('c.las', cat(hd, recs)), /laszip|Decompress/);
  const pts = await imp('scan.pts', '3\n1 2 3 100 255 0 0\n4 5 6 90 0 255 0\n7 8 9 80 0 0 255\n');
  ok(pts.kind === 'points' && pts.count === 3 && pts.points[3] === 4 && pts.points[8] === 9 && pts.warnings.some((w) => /ignored/.test(w)), 'PTS with count line, intensity and colour');
  const ptx = await imp('scan.ptx', '2\n2\n0 0 0\n1 0 0\n0 1 0\n0 0 1\n1 0 0 0\n0 1 0 0\n0 0 1 0\n10 20 30 1\n1 1 1 0.5\n0 0 0 0.5\n2 2 2 0.5\n3 3 3 0.5\n');
  ok(ptx.kind === 'points' && ptx.count === 3 && ptx.points[0] === 11 && ptx.points[1] === 21 && ptx.points[8] === 33 && ptx.warnings.some((w) => /empty/.test(w)), 'PTX: registration transform applied, empty return dropped');
  const xyz = await imp('cloud.xyz', '0.5 1.5 2.5\n1 2 3\n-4 5 6\n');
  ok(xyz.kind === 'points' && xyz.count === 3 && xyz.bbox.min[0] === -4 && xyz.bbox.max[2] === 6 && xyz.format === 'XYZ points / soundings', 'XYZ point list');
  const asc = await imp('cloud.asc', '// X Y Z I\n1 2 3 9\n4 5 6 9\n');
  ok(asc.kind === 'points' && asc.count === 2 && asc.pathway === 'points', '.asc without a grid header → point list');
  const csv = await imp('pts.csv', 'x,y,z\n0,0,1\n1,0,2\n0,1,3\n1,1,5\n2,2,2\n');
  ok(csv.kind === 'points' && csv.count === 5 && csv.points[14] === 2 && csv.pathway === 'points', 'CSV with x, y, z header → points');
  const txt = await imp('pts.txt', '1\t2\t3\n4\t5\t6\n');
  ok(txt.kind === 'points' && txt.count === 2, 'tab-separated TXT → points');
});

// ---- voxel volumes ---------------------------------------------------------------------------------------
function npy(descr, shape, bytes, fortran = false, major = 1) {
  const pre = major === 1 ? 10 : 12;
  let h = `{'descr': '${descr}', 'fortran_order': ${fortran ? 'True' : 'False'}, 'shape': (${shape.join(', ')}${shape.length === 1 ? ',' : ''}), }`;
  h = h.padEnd(Math.ceil((pre + h.length + 1) / 64) * 64 - pre - 1) + '\n';
  return cat(le('u8', [0x93]), 'NUMPY', le('u8', [major, 0]), major === 1 ? le('u16', [h.length]) : le('u32', [h.length]), h, bytes);
}
await group('image / voxel volumes', async () => {
  // 8 × 6 × 4 block (x, y, z), solid where x < 2 → porosity 0.75
  const NX = 8, NY = 6, NZ = 4, blk = new Uint8Array(NX * NY * NZ);
  for (let k = 0; k < NZ; k++) for (let j = 0; j < NY; j++) for (let i = 0; i < NX; i++) blk[(k * NY + j) * NX + i] = i < 2 ? 1 : 0;
  const isBlock = (g, label, { sp } = {}) => {
    ok(g.kind === 'voxels' && g.voxels.nx === NX && g.voxels.ny === NY && g.voxels.nz === NZ && g.voxels.data instanceof Uint8Array && g.voxels.data.length === blk.length, `${label}: ${NX} × ${NY} × ${NZ} voxels (kind ${g.kind}, ${g.voxels && [g.voxels.nx, g.voxels.ny, g.voxels.nz]})`);
    ok(g.voxels.data.every((v, i) => v === blk[i]) && Math.abs(g.stats.porosity - 0.75) < 1e-12 && Math.abs(microstructure(g).porosity - 0.75) < 1e-12, `${label}: solid block reproduced, porosity ${g.stats.porosity}`);
    if (sp) ok(g.voxels.spacing.every((v, k) => Math.abs(v - sp[k]) < 1e-6), `${label}: spacing ${g.voxels.spacing}`);
  };
  isBlock(await imp('block.npy', npy('|b1', [NZ, NY, NX], blk)), 'NPY bool, C order (z, y, x)');
  const grey = Float32Array.from(blk, (v, i) => (v ? 200 : 20) + (i % 7) * 0.25);
  isBlock(await imp('block.npy', npy('<f4', [NZ, NY, NX], le('f32', Array.from(grey)), false, 2)), 'NPY float32 v2 header, Otsu threshold');
  const fort = new Uint8Array(blk.length);            // Fortran order with shape (nx, ny, nz) has the same memory layout
  fort.set(blk);
  isBlock(await imp('block.npy', npy('>i2', [NX, NY, NZ], be('i16', Array.from(fort)), true, 3)), 'NPY big-endian int16, Fortran order, v3 header');
  const inv = await imp('block.npy', npy('|u1', [NZ, NY, NX], blk), { invert: true });
  ok(Math.abs(inv.stats.porosity - 0.25) < 1e-12 && inv.voxels.data[0] === 0, 'opts.invert swaps solid and pore');
  const np = await imp('pts.npy', npy('<f8', [5, 3], le('f64', [0, 0, 0, 1, 0, 0.5, 0, 1, 0.25, 1, 1, 2, 3, 3, 3])));
  ok(np.kind === 'points' && np.count === 5 && np.points[5] === 0.5 && np.bbox.max[2] === 3, 'NPY N × 3 float array → points');
  const img = await imp('img.npy', npy('|u1', [3, 4], le('u8', [1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0])));
  ok(img.kind === 'voxels' && img.voxels.nx === 4 && img.voxels.ny === 3 && img.voxels.nz === 1 && img.stats.porosity === 0.5, 'NPY 2-D integer array → 2-D voxels');
  const dm = await imp('dem.npy', npy('<f8', [2, 4], le('f64', [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5])));
  ok(dm.kind === 'grid' && dm.grid.z[1][2] === 6.5 && dm.grid.x.length === 4 && dm.grid.y.length === 2, 'NPY 2-D real array → grid');
  const lab = await imp('labels.npy', npy('|u1', [1, 2, 4], le('u8', [0, 1, 2, 3, 0, 0, 5, 5]))), two = await imp('two.npy', npy('|i1', [1, 2, 4], le('i8', [-1, 1, 1, -1, -1, -1, 1, -1])));
  ok(lab.stats.method === 'labels' && lab.stats.porosity === 3 / 8 && two.stats.method === 'binary' && two.stats.porosity === 5 / 8, 'labelled volume: non-zero = solid; two-valued volume: the larger value = solid');
  const z1 = await imp('arrays.npz', zip([['volume.npy', npy('|b1', [NZ, NY, NX], blk)], ['other.npy', npy('<f8', [2], le('f64', [1, 2]))]], true));
  isBlock(z1, 'NPZ (deflated), first array');
  ok(z1.stats.array === 'volume' && z1.stats.arrays === 2 && z1.warnings.some((w) => /2 arrays/.test(w)), 'NPZ reports which array was used');
  const nh = (enc2, extra = '') => `NRRD0004\n# comment\ntype: unsigned char\ndimension: 3\nsizes: ${NX} ${NY} ${NZ}\nspacings: 0.5 0.5 2\nencoding: ${enc2}\nendian: little\n${extra}\n`;
  isBlock(await imp('block.nrrd', cat(nh('raw'), blk)), 'NRRD raw', { sp: [0.5, 0.5, 2] });
  isBlock(await imp('block.nrrd', cat(nh('gzip'), new Uint8Array(zlib.gzipSync(blk)))), 'NRRD gzip');
  isBlock(await imp('block.nrrd', nh('ascii') + Array.from(blk).join(' ') + '\n'), 'NRRD ascii');
  isBlock(await imp('block.nhdr', `NRRD0004\ntype: uint8\ndimension: 3\nsizes: ${NX} ${NY} ${NZ}\nspace directions: (0.1,0,0) (0,0.2,0) (0,0,0.3)\nencoding: raw\ndata file: block.raw\n`, { companion: { 'block.raw': F('block.raw', blk) } }), 'NRRD detached header + companion', { sp: [0.1, 0.2, 0.3] });
  const mh = (file, extra = '') => `ObjectType = Image\nNDims = 3\nBinaryData = True\nBinaryDataByteOrderMSB = False\n${extra}ElementSpacing = 0.25 0.25 0.5\nDimSize = ${NX} ${NY} ${NZ}\nElementType = MET_USHORT\nElementDataFile = ${file}\n`;
  const u16 = le('u16', Array.from(blk, (v) => v * 1000));
  isBlock(await imp('block.mha', cat(mh('LOCAL'), u16)), 'MetaImage .mha inline', { sp: [0.25, 0.25, 0.5] });
  isBlock(await imp('block.mhd', mh('block.zraw', 'CompressedData = True\n'), { companion: { 'block.zraw': F('block.zraw', new Uint8Array(zlib.deflateSync(u16))) } }), 'MetaImage .mhd + zlib-compressed companion');
  await throws('MetaImage .mhd without its data file', () => imp('block.mhd', mh('block.raw')), /block\.raw/);
  const nii = (little) => { const h = new Uint8Array(352), dv = new DataView(h.buffer); dv.setInt32(0, 348, little); [3, NX, NY, NZ, 1, 1, 1, 1].forEach((v, k) => dv.setInt16(40 + 2 * k, v, little)); dv.setInt16(70, 16, little); dv.setInt16(72, 32, little); [1, 0.5, 0.5, 1.5].forEach((v, k) => dv.setFloat32(76 + 4 * k, v, little)); dv.setFloat32(108, 352, little); dv.setFloat32(112, 2, little); dv.setFloat32(116, -1, little); h.set(enc.encode('n+1'), 344); return cat(h, le('f32', Array.from(blk), little)); };
  isBlock(await imp('block.nii', nii(true)), 'NIfTI-1 float32 with scl_slope / scl_inter', { sp: [0.5, 0.5, 1.5] });
  isBlock(await imp('block.nii.gz', new Uint8Array(zlib.gzipSync(nii(false)))), 'NIfTI-1 big-endian, gzip (.nii.gz)');
  const cube = new Uint8Array(27).map((_, i) => (i % 3 === 0 ? 255 : 0)), rc = await imp('cube.raw', cube);
  ok(rc.kind === 'voxels' && rc.voxels.nx === 3 && rc.voxels.nz === 3 && Math.abs(rc.stats.porosity - 2 / 3) < 1e-12 && rc.warnings.some((w) => /inferred/.test(w)), 'RAW 27 bytes → 3 × 3 × 3 cube inferred');
  isBlock(await imp('block.vol', u16, { dims: [NX, NY, NZ], dtype: 'uint16', spacing: [2, 2, 2] }), 'RAW with dims + dtype', { sp: [2, 2, 2] });
  await throws('RAW without dimensions', () => imp('odd.raw', new Uint8Array(30)), /dims/);
  await throws('RAW smaller than its dimensions', () => imp('small.raw', u16, { dims: [100, 100, 100], dtype: 'uint16' }), /too small/);
  // DICOM: explicit VR little endian, one 4 × 3 slice of unsigned shorts, with an undefined-length sequence to skip
  const tag = (g, e) => le('u16', [g, e]), even = (b) => (b.length % 2 ? cat(b, ' ') : bytesOf(b));
  const el = (g, e, vr, val) => { const v = even(val); return ['OB', 'OW', 'SQ', 'UN'].includes(vr) ? cat(tag(g, e), vr, le('u16', [0]), le('u32', [v.length]), v) : cat(tag(g, e), vr, le('u16', [v.length]), v); };
  const us = (v) => le('u16', [v]), px = [0, 0, 900, 900, 0, 0, 900, 900, 10, 10, 800, 800];
  const seq = cat(tag(0x0008, 0x1140), 'SQ', le('u16', [0]), le('u32', [0xffffffff]), tag(0xfffe, 0xe000), le('u32', [0xffffffff]), el(0x0008, 0x1150, 'UI', '1.2.3'), tag(0xfffe, 0xe00d), le('u32', [0]), tag(0xfffe, 0xe0dd), le('u32', [0]));
  const dcm = (ts, pixel) => cat(new Uint8Array(128), 'DICM', el(0x0002, 0x0010, 'UI', ts), seq, el(0x0018, 0x0050, 'DS', '2.5'), el(0x0028, 0x0002, 'US', us(1)), el(0x0028, 0x0004, 'CS', 'MONOCHROME2'), el(0x0028, 0x0010, 'US', us(3)), el(0x0028, 0x0011, 'US', us(4)), el(0x0028, 0x0030, 'DS', '0.5\\0.25'), el(0x0028, 0x0100, 'US', us(16)), el(0x0028, 0x0103, 'US', us(0)), pixel);
  const dc = await imp('slice.dcm', dcm('1.2.840.10008.1.2.1', el(0x7fe0, 0x0010, 'OW', le('u16', px))), { threshold: 400 });
  ok(dc.kind === 'voxels' && dc.voxels.nx === 4 && dc.voxels.ny === 3 && dc.voxels.nz === 1 && dc.stats.porosity === 0.5 && dc.voxels.spacing[0] === 0.25 && dc.voxels.spacing[1] === 0.5 && dc.voxels.spacing[2] === 2.5, 'DICOM explicit-VR slice: 4 × 3, pixel spacing and thickness, sequence skipped');
  ok(dc.voxels.data[2] === 1 && dc.voxels.data[0] === 0 && dc.stats.method === 'threshold', 'DICOM thresholded at the given value with rows flipped to y-up');
  await throws('compressed DICOM transfer syntax', () => imp('j.dcm', dcm('1.2.840.10008.1.2.4.70', el(0x7fe0, 0x0010, 'OW', le('u16', px)))), /gdcmconv|Decompress|TIFF/);
  await throws('PNG outside a browser', () => imp('photo.png', cat(le('u8', [0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a]), new Uint8Array(32))), /browser/);
  const mtx = await imp('mask.txt', Array.from({ length: 4 }, (_, j) => Array.from({ length: 10 }, (_q, i) => (i < 5 ? 1 : 0)).join(' ')).join('\n') + '\n');
  ok(mtx.kind === 'voxels' && mtx.voxels.nx === 10 && mtx.voxels.ny === 4 && mtx.stats.porosity === 0.5, 'numeric 0/1 matrix in a text file → 2-D voxels');
  const jv = await imp('vox.json', JSON.stringify({ voxels: [[[1, 0], [1, 0]], [[1, 0], [0, 0]]], spacing: [2, 2, 2] }));
  ok(jv.kind === 'voxels' && jv.voxels.nx === 2 && jv.voxels.nz === 2 && jv.stats.porosity === 5 / 8 && jv.voxels.spacing[0] === 2, 'JSON nested 3-D array → voxels');
});

// ---- networks -----------------------------------------------------------------------------------------------
await group('plant / piping networks', async () => {
  const net = { nodes: [{ id: 'feed', type: 'source', x: 0, y: 0, z: 0 }, { id: 'P1', type: 'Pump', name: 'HP pump', x: 10, y: 0, z: 0 }, { id: 'RO', kind: 'membrane', position: [10, 0, 5] }, { tag: 'brine', class: 'sink', x: 20, y: 0, elevation: 5 }], pipes: [{ from: 'feed', to: 'P1', L: 10, D: 0.2 }, { source: 'P1', target: 'RO', dn: 150 }, { upstream: 'RO', downstream: 'brine', length: 12, diameter: 0.1, type: 'brine line', name: 'B-1' }] };
  const isNet = (g, label, { coords = true } = {}) => {
    const s = networkSummary(g);
    ok(g.kind === 'network' && s.nodes === 4 && s.edges === 3 && g.network.nodes.every((n) => typeof n.id === 'string' && typeof n.type === 'string' && typeof n.name === 'string'), `${label}: network of 4 nodes, 3 edges`);
    ok(Math.abs(s.totalLength - 27) < 1e-9 && s.byType.pump === 1 && s.byType.pipe === 2 && s.byType['brine line'] === 1 && s.order.join('>') === 'feed>P1>RO>brine' && s.acyclic, `${label}: total length ${s.totalLength}, byType, topological order ${s.order.join('>')}`);
    ok(s.pipes.length === 3 && s.pipes[0].diameter === 0.2 && Math.abs(s.pipes[1].diameter - 0.15) < 1e-12 && s.pipes[2].name === 'B-1' && (!coords || s.maxElevationChange === 5), `${label}: pipe list (DN 150 → 0.15 m), elevation change ${s.maxElevationChange}`);
  };
  const nj = await imp('plant.json', JSON.stringify(net));
  isNet(nj, 'JSON network');
  ok(nj.warnings.some((w) => /computed from node coordinates/.test(w)) && nj.pathway === 'numeric' && nj.network.nodes[1].name === 'HP pump' && nj.network.nodes[2].z === 5, 'JSON network: liberal keys (kind, position, tag, elevation), computed length noted');
  const yaml = `# desalination train\nplant:\n  name: "Train A"   # comment\n  units:\n    - id: feed\n      type: source\n      x: 0\n      y: 0\n      z: 0\n    - {id: P1, type: Pump, name: 'HP pump', x: 10, y: 0, z: 0}\n    - id: RO\n      kind: membrane\n      position: [10, 0, 5]\n    - tag: brine\n      class: sink\n      x: 20\n      y: 0\n      elevation: 5\n  streams:\n  - from: feed\n    to: P1\n    L: 10\n    D: 0.2\n  - {source: P1, target: RO, dn: 150}\n  - upstream: RO\n    downstream: brine\n    length: 12.0\n    diameter: 1.0e-1\n    type: brine line\n    name: B-1\n`;
  isNet(await imp('plant.yaml', yaml), 'YAML network (block + inline maps, nested under "plant")');
  isNet(await imp('plant.xml', `<?xml version="1.0"?><plant><equipment><unit id="feed" type="source" x="0" y="0" z="0"/><unit id="P1" type="Pump" x="10" y="0" z="0"><name>HP pump</name></unit><unit id="RO" kind="membrane" x="10" y="0" z="5"/><unit tag="brine" class="sink" x="20" y="0" elevation="5"/></equipment><piping><pipe from="feed" to="P1" L="10" D="0.2"/><pipe source="P1" target="RO" dn="150"/><pipe upstream="RO" downstream="brine" length="12" diameter="0.1" type="brine line" name="B-1"/></piping></plant>`), 'XML network');
  const csvNet = await imp('links.csv', 'from,to,length,diameter,type,name\nfeed,P1,10,0.2,pipe,L1\nP1,RO,5,0.15,pipe,L2\nRO,brine,12,0.1,brine line,B-1\n');
  ok(csvNet.kind === 'network' && csvNet.pathway === 'network' && networkSummary(csvNet).totalLength === 27 && networkSummary(csvNet).order.join('>') === 'feed>P1>RO>brine' && networkSummary(csvNet).byType['brine line'] === 1, 'CSV connectivity table → network');
  const loop = networkSummary(await imp('loop.json', JSON.stringify({ edges: [['a', 'b', 1], ['b', 'c', 1], ['c', 'a', 1]] })));
  ok(loop.nodes === 3 && loop.order.length === 0 && loop.acyclic === false && loop.totalLength === 3, 'cyclic network: no topological order');
  await imp('evil.json', '{"__proto__": {"polluted": 1}, "constructor": {"prototype": {"polluted": 1}}, "nodes": [{"id": "a", "__proto__": {"polluted": 1}}], "edges": [{"from": "a", "to": "b"}]}');
  await imp('evil.yaml', '__proto__:\n  polluted: 1\nnodes:\n  - id: a\n    __proto__: {polluted: 1}\nedges:\n  - {from: a, to: b}\n');
  ok(({}).polluted === undefined && Object.prototype.polluted === undefined, '__proto__ / constructor keys in JSON and YAML do not pollute prototypes');
  await throws('YAML with an anchor', () => imp('a.yaml', 'base: &b\n  x: 1\nuse: *b\n'), /anchors/);
  // PCF: three pipes, an elbow and a valve (millimetres)
  const pcf = `ISOGEN-FILES ISOGEN.FLS\nUNITS-BORE MM\nUNITS-CO-ORDS MM\nUNITS-WEIGHT KGS\nPIPELINE-REFERENCE LINE-100\n    PIPING-SPEC CS150\nPIPE\n    END-POINT 0 0 0 100\n    END-POINT 1000 0 0 100\n    ITEM-CODE P1\nELBOW\n    END-POINT 1000 0 0 100\n    END-POINT 1150 150 0 100\n    CENTRE-POINT 1150 0 0\n    SKEY ELBW\nPIPE\n    END-POINT 1150.0004 150 0 100\n    END-POINT 1150 2150 0 100\nVALVE\n    END-POINT 1150 2150 0 100\n    END-POINT 1150 2350 0 100\n    TAG XV-101\nPIPE\n    END-POINT 1150 2350 0 100\n    END-POINT 1150 2350 3000 100\nMATERIALS\nITEM-CODE P1\n    DESCRIPTION PIPE\n`;
  const gp = await imp('line.pcf', pcf), sp = networkSummary(gp);
  ok(gp.kind === 'network' && sp.nodes === 6 && sp.edges === 5 && sp.byType.pipe === 3 && sp.byType.elbow === 1 && sp.byType.valve === 1, `PCF: 6 merged nodes, 5 components (${sp.nodes}, ${sp.edges})`);
  ok(Math.abs(sp.totalLength - (1 + 0.15 * Math.PI / 2 + 2 + 0.2 + 3)) < 1e-6 && sp.pipes.every((p) => p.diameter === 0.1) && sp.maxElevationChange === 3 && gp.stats.units === 'm' && sp.pipes[3].name === 'XV-101', `PCF: lengths in metres incl. elbow arc (${sp.totalLength}), bore 0.1 m, rise 3 m`);
  ok(sp.order.length === 6 && sp.order[0] === 'N1' && sp.sources.length === 1 && sp.sinks.length === 1, 'PCF: line is a single chain from source to sink');
  // IFC: pump, two pipe segments (extrusion depth / length quantity), a tank, ports and connections; millimetre model
  const ifc = `ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('ViewDefinition [CoordinationView]'),'2;1');\nFILE_NAME('p.ifc','2024-01-01T00:00:00',(''),(''),'','','');\nFILE_SCHEMA(('IFC4'));\nENDSEC;\nDATA;\n#1=IFCSIUNIT(*,.LENGTHUNIT.,.MILLI.,.METRE.);\n#2=IFCUNITASSIGNMENT((#1));\n#10=IFCCARTESIANPOINT((0.,0.,0.));\n#11=IFCAXIS2PLACEMENT3D(#10,$,$);\n#12=IFCLOCALPLACEMENT($,#11);\n#13=IFCCARTESIANPOINT((1000.,0.,500.));\n#14=IFCAXIS2PLACEMENT3D(#13,$,$);\n#15=IFCLOCALPLACEMENT(#12,#14);\n#16=IFCDIRECTION((0.,1.,0.));\n#17=IFCDIRECTION((1.,0.,0.));\n#18=IFCAXIS2PLACEMENT3D(#13,#17,#16);\n#19=IFCLOCALPLACEMENT(#15,#18);\n#21=IFCLOCALPLACEMENT(#19,#14);\n#20=IFCPUMP('guidPump',$,'Feed pump P-101',$,$,#15,$,'P-101',.CIRCULATOR.);\n#30=IFCCIRCLEPROFILEDEF(.AREA.,$,$,50.);\n#31=IFCDIRECTION((0.,0.,1.));\n#32=IFCEXTRUDEDAREASOLID(#30,#11,#31,2500.);\n#33=IFCSHAPEREPRESENTATION($,'Body','SweptSolid',(#32));\n#34=IFCPRODUCTDEFINITIONSHAPE($,$,(#33));\n#40=IFCPIPESEGMENT('guidPipe1',$,'Pipe 1',$,$,#12,#34,$,.RIGIDSEGMENT.);\n#41=IFCPIPESEGMENT('guidPipe2',$,'Pipe 2',$,$,#21,$,$,.RIGIDSEGMENT.);\n#50=IFCQUANTITYLENGTH('Length',$,$,4000.,$);\n#51=IFCELEMENTQUANTITY('q',$,'Qto_PipeSegmentBaseQuantities',$,$,(#50));\n#52=IFCRELDEFINESBYPROPERTIES('r',$,$,$,(#41),#51);\n#60=IFCDISTRIBUTIONPORT('p1',$,'out',$,$,$,$,.SOURCE.,.PIPE.,.PIPING.);\n#61=IFCDISTRIBUTIONPORT('p2',$,'in',$,$,$,$,.SINK.,.PIPE.,.PIPING.);\n#62=IFCDISTRIBUTIONPORT('p3',$,'out',$,$,$,$,.SOURCE.,.PIPE.,.PIPING.);\n#63=IFCDISTRIBUTIONPORT('p4',$,'in',$,$,$,$,.SINK.,.PIPE.,.PIPING.);\n#70=IFCRELNESTS('n1',$,$,$,#20,(#60));\n#71=IFCRELNESTS('n2',$,$,$,#40,(#61,#62));\n#72=IFCRELNESTS('n3',$,$,$,#41,(#63));\n#80=IFCRELCONNECTSPORTS('c1',$,$,$,#60,#61,$);\n#81=IFCRELCONNECTSPORTS('c2',$,$,$,#63,#62,$);\n#90=IFCTANK('guidTank',$,'Permeate tank',$,$,$,$,$,$);\n#91=IFCPUMPTYPE('guidType',$,'type',$,$,$,$,$,$,.CIRCULATOR.);\nENDSEC;\nEND-ISO-10303-21;\n`;
  const gi = await imp('plant.ifc', ifc), si = networkSummary(gi), byId = Object.fromEntries(gi.network.nodes.map((n) => [n.id, n]));
  const gz = await imp('plant.ifczip', zip([['readme.txt', 'x'], ['model/plant.ifc', ifc]], true));
  ok(gz.kind === 'network' && gz.format === 'IFCZIP' && gz.network.nodes.length === gi.network.nodes.length && gz.network.edges.length === gi.network.edges.length, 'IFCZIP: the .ifc inside the archive is read like a plain IFC file');
  await throws('IFCZIP holding only ifcXML', () => importGeometry(F('a.ifczip', zip([['a.ifcxml', '<ifcXML/>']]))), /ifcXML/);
  ok(gi.kind === 'network' && si.nodes === 4 && si.edges === 2 && si.byType.pump === 1 && si.byType.pipe === 2 && si.byType.tank === 1 && gi.format === 'IFC (STEP text)', `IFC: pump + 2 pipe segments + tank, 2 port connections (${si.nodes}, ${si.edges})`);
  ok(byId.guidPump.name === 'Feed pump P-101' && byId.guidPump.tag === 'P-101' && Math.abs(byId.guidPump.x - 1) < 1e-12 && Math.abs(byId.guidPump.z - 0.5) < 1e-12 && byId.guidTank.x === null, 'IFC: names, tags and placements in metres (nested IfcLocalPlacement)');
  ok(Math.abs(byId.guidPipe2.x - 2.5) < 1e-12 && Math.abs(byId.guidPipe2.y - 1) < 1e-12 && Math.abs(byId.guidPipe2.z - 1) < 1e-12,`IFC: rotated child placement composed with its parent (${byId.guidPipe2.x}, ${byId.guidPipe2.y}, ${byId.guidPipe2.z})`);
  ok(Math.abs(byId.guidPipe1.length - 2.5) < 1e-12 && Math.abs(byId.guidPipe1.diameter - 0.1) < 1e-12 && Math.abs(byId.guidPipe2.length - 4) < 1e-12 && Math.abs(si.totalLength - 6.5) < 1e-12, 'IFC: pipe lengths from extrusion depth and Qto length quantity');
  ok(gi.network.edges[0].from === 'guidPump' && gi.network.edges[0].to === 'guidPipe1' && gi.network.edges[1].from === 'guidPipe1' && gi.network.edges[1].to === 'guidPipe2' && gi.warnings.some((w) => /2 of 2 segments/.test(w)), 'IFC: connectivity through nested ports, directed from SOURCE to SINK; recovery reported');
  const aml = `<?xml version="1.0" encoding="utf-8"?><CAEXFile FileName="plant.aml" SchemaVersion="2.15"><InstanceHierarchy Name="Plant"><InternalElement Name="Train1" ID="t1"><InternalElement Name="HP pump" ID="e1"><ExternalInterface Name="out" ID="i1"/><RoleRequirements RefBaseRoleClassPath="ProcessRoles/Pump"/></InternalElement><InternalElement Name="RO vessel" ID="e2"><ExternalInterface Name="in" ID="i2"/><ExternalInterface Name="brine" ID="i3"/></InternalElement><InternalElement Name="ERD" ID="e3" RefBaseSystemUnitPath="Lib/EnergyRecovery"><ExternalInterface Name="in" ID="i4"/></InternalElement><InternalLink Name="feed" RefPartnerSideA="i1" RefPartnerSideB="i2"/><InternalLink Name="brine" RefPartnerSideA="e2:brine" RefPartnerSideB="{e3}:in"/></InternalElement></InstanceHierarchy></CAEXFile>`;
  const ga = await imp('plant.aml', aml), sa = networkSummary(ga);
  ok(ga.kind === 'network' && sa.nodes === 4 && sa.edges === 2 && sa.byType.pump === 1 && sa.byType.energyrecovery === 1 && ga.network.nodes[1].parent === 't1', 'AutomationML: InternalElement hierarchy → 4 nodes with roles');
  ok(ga.network.edges[0].from === 'e1' && ga.network.edges[0].to === 'e2' && ga.network.edges[1].from === 'e2' && ga.network.edges[1].to === 'e3' && ga.network.edges[1].name === 'brine', 'AutomationML: InternalLinks resolved by interface ID and by element:interface');
  const viaXml = await imp('plant.xml', aml);
  ok(viaXml.kind === 'network' && viaXml.format === 'AutomationML / CAEX', 'generic .xml routed by its root element');
  const prm = await imp('element.json', JSON.stringify({ element: { diameter_in: 8, length_m: 1.016, area_m2: 37 }, spacer: { thickness_mil: 34, angle: 90 } }));
  ok(prm.kind === 'table' && prm.records.length === 1 && prm.records[0]['spacer.angle'] === 90 && prm.params.element.area_m2 === 37, 'JSON parameter set (parametric membrane definition) → one-row table + params');
  const tb = await imp('notes.csv', 'name,value,unit\nrecovery,0.45,-\npressure,62,bar\n');
  ok(tb.kind === 'table' && tb.headers.length === 3 && tb.records.length === 2 && tb.records[1].value === 62, 'unrecognised CSV columns → plain table');
});

// ---- analysis helpers -------------------------------------------------------------------------------------
await group('sectionOf / maskOf / gridOf', async () => {
  const cube = await imp('cube.obj', `${V.map((v) => 'v ' + v.join(' ')).join('\n')}\n${Q.map((q) => 'f ' + q.map((i) => i + 1).join(' ')).join('\n')}\n`);
  for (const axis of [0, 1, 2]) {
    const s = sectionOf(cube, { axis });
    const len = s.reduce((a, q) => a + Math.hypot(q[2] - q[0], q[3] - q[1]), 0), flat = s.flat();
    ok(s.length === 8 && Math.abs(len - 4) < 1e-12 && Math.min(...flat) === 0 && Math.max(...flat) === 1, `sectionOf(cube, axis ${axis}): unit-square outline, perimeter ${len}`);
  }
  ok(sectionOf(cube, { axis: 2, value: 5 }).length === 0, 'sectionOf outside the mesh gives no segments');
  const frac = (m) => m.flat().filter(Boolean).length / (m.length * m[0].length);
  const m1 = maskOf(cube, 0, 2, 0, 2, 40, 40, { margin: 0.25 });
  ok(m1.length === 40 && m1[0].length === 40 && frac(m1) === 0.25 && m1[20][20] === true && m1[2][2] === false && typeof m1[0][0] === 'boolean', `maskOf(cube) fitted with a 25 % margin: area fraction ${frac(m1)}`);
  ok(frac(maskOf(cube, 0, 2, 0, 2, 40, 40)) === 1 && Math.abs(frac(maskOf(cube, 0, 4, 0, 2, 80, 40)) - 0.5) < 1e-12 && frac(maskOf(cube, 0, 4, 0, 2, 80, 40, { stretch: true })) === 1, 'maskOf fit: uniform scale by default, stretch fills the window');
  ok(Math.abs(frac(maskOf(cube, -1, 2, -1, 2, 30, 30, { fit: false })) - 1 / 9) < 1e-12 && Math.abs(frac(maskOf(cube, -1, 2, -1, 2, 30, 30, { fit: false, invert: true })) - 8 / 9) < 1e-12, 'maskOf without fit keeps model coordinates; invert flips');
  const sq = await imp('sq.xy', '0 0\n2 0\n2 1\n0 1\n0 0\n');
  ok(sectionOf(sq).length === 4 && frac(maskOf(sq, 0, 4, 0, 4, 40, 40)) === 0.5, 'sectionOf / maskOf on polylines (2 × 1 rectangle fills half of a square window)');
  const sp = generate({ type: 'spacer', height: 1, length: 8, diameter: 0.5, pitch: 2, n: 320 });
  const mv = maskOf(sp, 0, 8, 0, 1, sp.voxels.nx, sp.voxels.ny);
  let same = 0;
  for (let j = 0; j < sp.voxels.ny; j++) for (let i = 0; i < sp.voxels.nx; i++) same += mv[j][i] === (sp.voxels.data[j * sp.voxels.nx + i] === 1) ? 1 : 0;
  ok(same === sp.voxels.nx * sp.voxels.ny, 'maskOf(voxels) on the native grid reproduces the voxel data');
  const sv = sectionOf(sp), per = sv.reduce((a, q) => a + Math.hypot(q[2] - q[0], q[3] - q[1]), 0);
  ok(sv.length > 100 && Math.abs(per - 4 * Math.PI * 0.5) / (2 * Math.PI) < 0.08, `sectionOf(voxels): marching-squares outline of 4 filaments, perimeter ${per.toFixed(3)} vs ${(2 * Math.PI).toFixed(3)}`);
  const blk = generate({ type: 'csg', tree: { type: 'box', min: [0.25, 0.25, 0.25], max: [0.75, 0.75, 0.75] }, bounds: { min: [0, 0, 0], max: [1, 1, 1] }, n: 20 });
  const s3 = sectionOf(blk, { axis: 0, value: 0.5 }), f3 = s3.flat();
  ok(s3.length > 0 && Math.abs(Math.min(...f3) - 0.25) < 0.03 && Math.abs(Math.max(...f3) - 0.75) < 0.03 && sectionOf(blk, { axis: 0, value: 0.05 }).length === 0, 'sectionOf(3-D voxels, axis, value) cuts the chosen slice');
  // gridOf: scattered samples of the plane z = 2x + 3y + 1, with an empty patch that must be filled
  let seed = 11; const rnd = () => ((seed = (seed * 48271) % 2147483647) / 2147483647), pts = [];
  while (pts.length < 3 * 6000) { const x = rnd(), y = rnd(); if (x > 0.42 && x < 0.58 && y > 0.42 && y < 0.58) continue; pts.push(x, y, 2 * x + 3 * y + 1); }
  const pc = await imp('plane.npy', npy('<f8', [pts.length / 3, 3], le('f64', pts))), gr = gridOf(pc, 10, 10);
  let worst = 0;
  for (let j = 0; j < 10; j++) for (let i = 0; i < 10; i++) worst = Math.max(worst, Math.abs(gr.z[j][i] - (2 * gr.x[i] + 3 * gr.y[j] + 1)));
  ok(gr.x.length === 10 && gr.y.length === 10 && gr.z.length === 10 && gr.z[0].length === 10 && gr.z.flat().every(Number.isFinite) && worst < 1e-9, `gridOf(scattered points) recovers the plane, worst cell error ${worst.toExponential(2)}`);
  const hole = Array.from(pts); for (let k = 0; k < hole.length; k += 3) if (hole[k] > 0.3 && hole[k] < 0.7 && hole[k + 1] > 0.3 && hole[k + 1] < 0.7) hole[k + 2] = NaN;
  const gh = gridOf({ kind: 'points', points: Float64Array.from(hole) }, 10, 10);
  ok(Math.abs(gh.z[4][4] - (2 * gh.x[4] + 3 * gh.y[4] + 1)) < 1e-9 && Math.abs(gh.z[5][5] - (2 * gh.x[5] + 3 * gh.y[5] + 1)) < 1e-9, 'gridOf fills empty cells from the surrounding bins');
  const dem = await imp('dem.asc', 'ncols 3\nnrows 2\nxllcorner 0\nyllcorner 0\ncellsize 1\n3 4 5\n0 1 2\n'), rg = gridOf(dem, 5, 3);
  ok(rg.z[0][0] === 0 && rg.z[2][4] === 5 && rg.z[1][2] === 2.5 && rg.z[0][1] === 0.5 && rg.x[4] === 2.5, 'gridOf(grid) resamples bilinearly');
  const dg = dimensions(dem), cs = sectionOf(dem, { value: 2.5 });
  ok(dg.area === 2 && cs.length >= 2 && cs.every((s) => s.every(Number.isFinite)), 'grid: extent area and iso-line section');
  const gm = gridOf(cube, 2, 2);
  ok(gm.z.flat().every((v) => v === 0.5), 'gridOf(mesh) bins the mesh vertices');
});

await group('microstructure', async () => {
  // slabs normal to y: channels run along x and z
  const nx = 12, ny = 8, nz = 6, data = new Uint8Array(nx * ny * nz);
  for (let k = 0; k < nz; k++) for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) data[(k * ny + j) * nx + i] = j % 4 < 2 ? 1 : 0;
  const m = microstructure({ kind: 'voxels', voxels: { nx, ny, nz, data, spacing: [1, 1, 1] } });
  ok(m.porosity === 0.5 && m.solidFraction === 0.5 && m.connectedPorosity === 0.5, 'parallel channels: porosity exactly 0.5, all of it connected');
  ok(m.tortuosity.length === 3 && m.tortuosity[0] === 1 && m.tortuosity[2] === 1 && m.tortuosity[1] === null && m.percolates.join() === 'true,false,true', `parallel channels: tortuosity 1 along the channels, no percolation across (${m.tortuosity})`);
  ok(Math.abs(m.specificSurface - 0.375) < 1e-12 && Math.abs(m.meanPoreSize - 4) < 1e-12 && Math.abs(m.meanSolidSize - 4) < 1e-12, `parallel channels: specific surface ${m.specificSurface}, mean chord ${m.meanPoreSize}`);
  const m2 = microstructure({ kind: 'voxels', voxels: { nx, ny, nz, data, spacing: [2, 2, 2] } });
  ok(Math.abs(m2.specificSurface - 0.1875) < 1e-12 && m2.meanPoreSize === 8 && m2.tortuosity[0] === 1, 'spacing scales lengths and specific surface');
  // serpentine channel in a 2-D mask (true = solid)
  const W = 21, H = 9, mask = Array.from({ length: H }, () => new Array(W).fill(false));
  for (const [x, gapTop] of [[5, true], [10, false], [15, true]]) for (let j = 0; j < H; j++) if (gapTop ? j < H - 1 : j > 0) mask[j][x] = true;
  const ms = microstructure(mask);
  ok(ms.tortuosity.length === 2 && ms.percolates[0] === true && ms.tortuosity[0] > 1.5 && ms.tortuosity[0] < 3 && Math.abs(ms.porosity - (1 - 24 / 189)) < 1e-12, `serpentine 2-D mask: tortuosity ${ms.tortuosity[0].toFixed(3)} along x`);
  const closed = mask.map((r) => r.slice()); closed[H - 1][5] = true;
  const mc = microstructure(closed);
  ok(mc.percolates[0] === false && mc.tortuosity[0] === null && mc.percolates[1] === true && mc.tortuosity[1] >= 1 && mc.tortuosity[1] < 1.05 && mc.poreClusters === 2, 'blocked serpentine: no percolation along x, still open along y');
  const ring = Array.from({ length: 7 }, (_, j) => Array.from({ length: 7 }, (_q, i) => Math.max(Math.abs(i - 3), Math.abs(j - 3)) === 2)), mr = microstructure(ring);
  ok(Math.abs(mr.porosity - 33 / 49) < 1e-12 && Math.abs(mr.connectedPorosity - 24 / 49) < 1e-12 && mr.poreClusters === 2, 'sealed cavity: connected porosity counts only clusters that span the sample');
  const diag = Array.from({ length: 30 }, (_, j) => Array.from({ length: 30 }, (_q, i) => Math.abs(i - j) > 2));
  const md = microstructure(diag);
  ok(md.tortuosity[0] > 1.3 && md.tortuosity[0] < 1.5, `45° channel: chamfer tortuosity ${md.tortuosity[0].toFixed(3)} ≈ √2`);
  const dv = dimensions({ kind: 'voxels', voxels: { nx, ny, nz, data, spacing: [0.5, 0.5, 0.5] }, bbox: { min: [0, 0, 0], max: [6, 4, 3] } });
  ok(dv.volume === 36 && dv.size.join() === '6,4,3' && Math.abs(dv.area - 0.75 * 72) < 1e-9, 'dimensions(voxels): solid volume by voxel count');
});

await group('generate()', async () => {
  const por = (g) => microstructure(g, { tortuosity: false }).porosity, sameData = (a, b) => a.voxels.data.length === b.voxels.data.length && a.voxels.data.every((v, i) => v === b.voxels.data[i]);
  const shape = (g, label) => ok(g.kind === 'voxels' && g.voxels.data instanceof Uint8Array && g.voxels.data.length === g.voxels.nx * g.voxels.ny * g.voxels.nz && g.voxels.spacing.length === 3 && g.pathway === 'procedural' && Array.isArray(g.warnings) && Math.abs(g.stats.porosity - por(g)) < 1e-12, `${label}: voxel geometry ${g.voxels.nx} × ${g.voxels.ny} × ${g.voxels.nz}, porosity ${g.stats.porosity.toFixed(4)}`);
  const gy = generate({ type: 'tpms', surface: 'gyroid', porosity: 0.6, n: 32 });
  shape(gy, 'tpms gyroid');
  ok(Math.abs(por(gy) - 0.6) < 0.002 && sameData(gy, generate({ type: 'tpms', surface: 'gyroid', porosity: 0.6, n: 32 })) && microstructure(gy).percolates.every(Boolean), 'gyroid hits its porosity target, is deterministic and percolates in x, y, z');
  for (const s of ['schwarzP', 'diamond']) { const g = generate({ type: 'tpms', surface: s, n: 32, period: 0.5 }); ok(Math.abs(por(g) - 0.5) < 0.02 && g.voxels.nx === 32, `tpms ${s} at level 0: porosity ${por(g).toFixed(3)} ≈ 0.5`); }
  const sheet = generate({ type: 'tpms', surface: 'gyroid', sheet: true, porosity: 0.8, n: 32 });
  ok(Math.abs(por(sheet) - 0.8) < 0.002, `sheet gyroid hits its porosity target (${por(sheet)})`);
  const vo = generate({ type: 'voronoi', cells: 12, porosity: 0.8, seed: 3, n: 32 });
  shape(vo, 'voronoi foam');
  ok(Math.abs(por(vo) - 0.8) < 0.002 && sameData(vo, generate({ type: 'voronoi', cells: 12, porosity: 0.8, seed: 3, n: 32 })) && !sameData(vo, generate({ type: 'voronoi', cells: 12, porosity: 0.8, seed: 4, n: 32 })), 'voronoi: porosity target met, same seed → same structure, other seed → different');
  const vw = generate({ type: 'voronoi', cells: 8, wall: 0.06, seed: 2, n: [96, 96, 1] });
  ok(vw.voxels.nz === 1 && por(vw) > 0.5 && por(vw) < 0.95 && microstructure(vw).tortuosity.length === 2, `2-D voronoi with wall thickness: porosity ${por(vw).toFixed(3)}`);
  const la = generate({ type: 'lattice', pitch: 1, diameter: 0.3, angle: 90, n: 96 });
  shape(la, 'filament lattice (spacer unit cell)');
  ok(Math.abs(por(la) - (1 - (Math.PI * 0.3) / 4)) < 0.02 && Math.abs(la.bbox.max[2] - 0.6) < 1e-12 && Math.abs(la.bbox.max[0] - 2 * Math.SQRT2) < 1e-9, `lattice: porosity ${por(la).toFixed(4)} vs analytic ${(1 - (Math.PI * 0.3) / 4).toFixed(4)} for touching crossed filaments`);
  const sph = generate({ type: 'spheres', porosity: 0.45, radius: 0.08, seed: 5, n: 48 });
  shape(sph, 'sphere pack');
  ok(Math.abs(por(sph) - 0.45) < 0.01 && sameData(sph, generate({ type: 'spheres', porosity: 0.45, radius: 0.08, seed: 5, n: 48 })) && !sameData(sph, generate({ type: 'spheres', porosity: 0.45, radius: 0.08, seed: 6, n: 48 })), `spheres: porosity ${por(sph).toFixed(4)} within 0.01 of 0.45, seeded`);
  const disc = generate({ type: 'spheres', porosity: 0.7, radius: 0.05, seed: 1, n: [128, 128, 1], overlap: false });
  ok(disc.voxels.nz === 1 && Math.abs(por(disc) - 0.7) < 0.01 && disc.stats.overlap === false, `non-overlapping discs: porosity ${por(disc).toFixed(4)}`);
  const csg = generate({ type: 'csg', tree: { op: 'difference', children: [{ type: 'box', min: [0, 0, 0], max: [1, 1, 1] }, { type: 'sphere', center: [0.5, 0.5, 0.5], radius: 0.4 }] }, bounds: { min: [0, 0, 0], max: [1, 1, 1] }, n: 40 });
  shape(csg, 'CSG box minus sphere');
  ok(Math.abs(1 - por(csg) - (1 - (4 / 3) * Math.PI * 0.064)) < 0.01, `CSG: solid fraction ${(1 - por(csg)).toFixed(4)} vs ${(1 - (4 / 3) * Math.PI * 0.064).toFixed(4)}`);
  const csg2 = generate({ type: 'csg', tree: { op: 'union', children: [{ type: 'cylinder', center: [0, 0, 0], radius: 0.5, height: 2, axis: 'z' }, { op: 'intersection', children: [{ type: 'box', center: [2, 0, 0], size: [1, 1, 1] }, { type: 'sphere', center: [2, 0, 0], radius: 0.6 }] }] }, n: 80 });
  const dv = dimensions(csg2);
  ok(csg2.bbox.min[0] < -0.5 && csg2.bbox.max[0] > 2.5 && csg2.voxels.origin[0] < -0.5 && Math.abs(dv.volume - (Math.PI * 0.5 + 0.9)) < 0.12, `CSG union with automatic bounds: volume ${dv.volume.toFixed(3)} (cylinder π/2 + clipped box)`);
  const im = generate({ type: 'implicit', fn: (x, y, z) => Math.hypot(x, y, z) - 0.4, bounds: { min: [-0.5, -0.5, -0.5], max: [0.5, 0.5, 0.5] }, n: 40 });
  shape(im, 'implicit sphere');
  ok(Math.abs(1 - por(im) - (4 / 3) * Math.PI * 0.064) < 0.01 && im.bbox.min[0] === -0.5, 'implicit function sampled inside its bounds (solid where fn < 0)');
  let threw = null;
  try { generate({ type: 'implicit', fn: 'Math.hypot(x,y,z)-0.4' }); } catch (e) { threw = e; }
  ok(threw instanceof Error && /function/.test(threw.message), 'implicit geometry refuses a string instead of a function');
  const sp = generate({ type: 'spacer', height: 1, length: 8, diameter: 0.5, pitch: 2, arrangement: 'zigzag', n: 400 });
  shape(sp, 'spacer-filled channel');
  const msp = microstructure(sp);
  ok(sp.voxels.nz === 1 && sp.voxels.nx === 400 && sp.voxels.ny === 50 && Math.abs(1 - msp.porosity - (4 * Math.PI * 0.0625) / 8) < 0.003 && msp.percolates[0] && msp.tortuosity[0] > 1 && msp.tortuosity[0] < 1.2 && sp.stats.filaments === 4, `spacer: 4 filaments, solid fraction ${(1 - msp.porosity).toFixed(4)}, tortuosity ${msp.tortuosity[0].toFixed(3)}`);
  const cav = generate({ type: 'spacer', height: 1, length: 8, diameter: 0.5, pitch: 2, arrangement: 'cavity', n: 400 });
  ok(cav.voxels.data.subarray(49 * 400).every((v) => v === 0) && cav.voxels.data.subarray(0, 400).some((v) => v === 1), 'cavity arrangement keeps every filament on the bottom wall');
  for (const bad of [{ type: 'nonsense' }, { type: 'csg', tree: { type: 'torus' } }, { type: 'tpms', surface: 'nope' }, { type: 'tpms', n: 5000 }]) { let e = null; try { generate(bad); } catch (x) { e = x; } ok(e instanceof Error && e.message.length < 300, `generate(${JSON.stringify(bad)}) throws: ${e && e.message.slice(0, 50)}`); }
});

await group('GraphML and LandXML', async () => {
  const gml = `<?xml version="1.0" encoding="UTF-8"?><graphml xmlns="http://graphml.graphdrawing.org/xmlns"><key id="d0" for="node" attr.name="x" attr.type="double"/><key id="d1" for="node" attr.name="y" attr.type="double"/><key id="d2" for="node" attr.name="elevation" attr.type="double"/><key id="d3" for="node" attr.name="type" attr.type="string"/><key id="d4" for="edge" attr.name="length" attr.type="double"/><key id="d5" for="edge" attr.name="diameter" attr.type="double"/><key id="d6" for="edge" attr.name="name" attr.type="string"/><key id="d7" for="edge" attr.name="type" attr.type="string"/>
<graph id="G" edgedefault="directed"><node id="WH1"><data key="d0">0</data><data key="d1">0</data><data key="d2">-1200</data><data key="d3">wellhead</data></node><node id="WH2"><data key="d0">0</data><data key="d1">800</data><data key="d2">-1210</data><data key="d3">wellhead</data></node><node id="MAN"><data key="d0">2500</data><data key="d1">400</data><data key="d2">-1180</data><data key="d3">manifold</data></node><node id="FPSO"><data key="d0">9000</data><data key="d1">400</data><data key="d2">0</data><data key="d3">topsides</data></node>
<edge id="e1" source="WH1" target="MAN"><data key="d4">2600</data><data key="d5">0.2</data><data key="d6">FL-1</data><data key="d7">flowline</data></edge><edge id="e2" source="WH2" target="MAN"><data key="d5">0.2</data></edge><edge id="e3" source="MAN" target="FPSO"><data key="d4">7500</data><data key="d5">0.3</data><data key="d6">RISER</data></edge><hyperedge/></graph></graphml>`;
  const g = await imp('field.graphml', gml), byId = Object.fromEntries(g.network.nodes.map((n) => [n.id, n])), s = networkSummary(g);
  ok(g.kind === 'network' && g.format === 'GraphML' && g.pathway === 'network' && g.network.nodes.length === 4 && g.network.edges.length === 3, 'GraphML: 4 nodes, 3 edges');
  ok(byId.WH1.z === -1200 && byId.MAN.type === 'manifold' && byId.FPSO.x === 9000 && g.network.edges[0].length === 2600 && g.network.edges[0].diameter === 0.2 && g.network.edges[0].name === 'FL-1' && g.network.edges[0].type === 'flowline', 'GraphML data keys: x, y, elevation, type, length, diameter, name');
  near(g.network.edges[1].length, Math.hypot(2500, 400, 30), 1e-9, 'GraphML edge without a length gets it from the node coordinates');
  ok(s.acyclic && s.sources.length === 2 && s.sinks[0] === 'FPSO' && g.warnings.some((w) => /hyperedge/.test(w)) && g.network.edges[1].name === 'e2', 'GraphML: directed graph summary, hyperedge warning, edge id as name');
  const viaXml = await imp('field.xml', gml);
  ok(viaXml.kind === 'network' && viaXml.format === 'GraphML' && viaXml.network.edges.length === 3, 'GraphML recognised inside a generic .xml');
  const yed = await imp('yed.graphml', '<graphml xmlns:y="http://www.yworks.com/xml/graphml"><key id="d6" for="node" yfiles.type="nodegraphics"/><graph><node id="n0"><data key="d6"><y:ShapeNode><y:Geometry height="30" width="30" x="10.5" y="20"/><y:NodeLabel>Choke</y:NodeLabel></y:ShapeNode></data></node><node id="n1"/><edge source="n0" target="n1"/></graph></graphml>');
  ok(yed.network.nodes[0].x === 10.5 && yed.network.nodes[0].y === 20 && yed.network.nodes[0].name === 'Choke', 'GraphML: yEd node geometry and label');
  await throws('GraphML without content', () => importGeometry(F('a.graphml', '<graphml><graph/></graphml>')), /no nodes|No nodes/);

  // LandXML: straight 300 m east, then a 100 m radius quarter circle turning north; profile rises 10 m over the first 300 m
  const arc = (Math.PI / 2) * 100, lx = `<?xml version="1.0"?><LandXML xmlns="http://www.landxml.org/schema/LandXML-1.2" version="1.2"><Units><Metric linearUnit="meter" areaUnit="squareMeter" diameterUnit="millimeter"/></Units>
<CgPoints><CgPoint name="P1">10 20 -3</CgPoint><CgPoint name="P2">30 40 -4</CgPoint></CgPoints>
<Alignments><Alignment name="FL-A" length="${300 + arc}" staStart="1000"><CoordGeom><Line><Start>0 0</Start><End>0 300</End></Line><Curve rot="ccw" radius="100"><Start>0 300</Start><Center>100 300</Center><End>100 400</End></Curve></CoordGeom>
<Profile><ProfAlign name="design"><PVI>1000 -100</PVI><PVI>1150 -95</PVI><PVI>1300 -90</PVI><PVI>${1300 + arc} -90</PVI></ProfAlign></Profile></Alignment></Alignments>
<PlanFeatures><PlanFeature name="umbilical"><CoordGeom><IrregularLine><Start>0 0 -5</Start><PntList3D>0 10 -6 0 20 -7</PntList3D><End>0 30 -8</End></IrregularLine></CoordGeom></PlanFeature></PlanFeatures>
<Surfaces><Surface name="seabed"><Definition surfType="TIN"><Pnts><P id="1">0 0 -100</P><P id="2">0 500 -90</P><P id="3">500 500 -80</P><P id="4">500 0 -95</P></Pnts><Faces><F>1 2 3</F><F>1 3 4</F><F i="1">1 2 4</F></Faces></Definition></Surface></Surfaces>
<PipeNetworks><PipeNetwork name="N"><Structs><Struct name="S1" elevRim="-50" elevSump="-52"><Center>0 0</Center></Struct><Struct name="S2" elevRim="-40"><Center>0 100</Center></Struct></Structs><Pipes><Pipe name="PIPE-1" refStart="S1" refEnd="S2" length="100.5"><CircPipe diameter="300"/></Pipe></Pipes></PipeNetwork></PipeNetworks></LandXML>`;
  const al = await imp('route.landxml', lx), p0 = al.polylines[0], n0 = p0.x.length;
  ok(al.kind === 'polylines' && al.format === 'LandXML' && al.pathway === 'gis' && al.polylines.length === 2 && al.stats.units === 'm', 'LandXML: alignment + plan feature as polylines');
  ok(p0.x[0] === 0 && p0.y[0] === 0 && Math.abs(p0.x[n0 - 1] - 400) < 1e-9 && Math.abs(p0.y[n0 - 1] - 100) < 1e-9 && p0.z[0] === -100 && p0.z[n0 - 1] === -90, 'LandXML: northing / easting swapped to x = east, y = north; profile gives z at both ends');
  const i150 = p0.x.findIndex((v) => Math.abs(v - 150) < 1e-9), i300 = p0.x.findIndex((v, k) => Math.abs(v - 300) < 1e-9 && Math.abs(p0.y[k]) < 1e-9);
  ok(i150 > 0 && Math.abs(p0.z[i150] + 95) < 1e-9 && Math.abs(p0.z[i300] + 90) < 1e-9, 'LandXML: a vertex is inserted at every PVI station (staStart honoured)');
  let l2 = 0; for (let i = 1; i < n0; i++) l2 += Math.hypot(p0.x[i] - p0.x[i - 1], p0.y[i] - p0.y[i - 1]);
  near(l2, 300 + arc, 0.2, 'LandXML: line + tessellated ccw curve has the alignment length');
  ok(al.polylines[1].z.join() === '-5,-6,-7,-8' && al.polylines[1].x.join() === '0,10,20,30' && al.bbox.min[2] === -100 && al.warnings.some((w) => /also holds/.test(w)), 'LandXML: 3-D IrregularLine plan feature; other content is announced');
  const sf = await imp('surface.landxml', lx, { prefer: 'surface' });
  ok(sf.kind === 'mesh' && sf.count === 2 && sf.bbox.max[0] === 500 && sf.bbox.min[2] === -100 && dimensions(sf).area > 250000 && dimensions(sf).area < 250200 && sf.stats.surfaces === 1, 'LandXML: TIN surface (Pnts + Faces, invisible face skipped) as mesh via prefer');
  const cp = await imp('points.landxml', lx, { prefer: 'points' });
  ok(cp.kind === 'points' && cp.count === 2 && cp.points[0] === 20 && cp.points[1] === 10 && cp.points[2] === -3, 'LandXML: CgPoints as points (east, north, elevation)');
  const pn = await imp('pipes.landxml', lx, { prefer: 'network' });
  ok(pn.kind === 'network' && pn.network.nodes.length === 2 && pn.network.nodes[0].z === -52 && pn.network.edges[0].length === 100.5 && Math.abs(pn.network.edges[0].diameter - 0.3) < 1e-12 && pn.network.edges[0].name === 'PIPE-1', 'LandXML: pipe network with diameters converted from millimetres');
  const lxml = await imp('route.xml', lx);
  ok(lxml.kind === 'polylines' && lxml.format === 'LandXML' && lxml.polylines.length === 2, 'LandXML recognised inside a generic .xml');
  const cw = await imp('cw.landxml', '<LandXML><Units><Imperial linearUnit="USSurveyFoot"/></Units><Alignments><Alignment name="a"><CoordGeom><Curve rot="cw"><Start>0 0</Start><Center>0 100</Center><End>100 100</End></Curve><Spiral><Start>100 100</Start><End>150 200</End></Spiral></CoordGeom></Alignment></Alignments></LandXML>');
  const pc = cw.polylines[0], mid = Math.floor((pc.x.length - 2) / 2);
  ok(cw.stats.units === 'ft' && !pc.z && Math.abs(Math.hypot(pc.x[mid] - 100, pc.y[mid]) - 100) < 1e-9 && pc.y[mid] > 0 && pc.x[mid] < 100 && cw.warnings.some((w) => /spiral/i.test(w)) && cw.warnings.some((w) => /feet/.test(w)), 'LandXML: clockwise curve on its circle, spiral replaced by its chord, feet reported');
  await throws('LandXML without geometry', () => importGeometry(F('a.landxml', '<LandXML><Units/></LandXML>')), /no alignments/);
  await throws('.landxml that is not LandXML', () => importGeometry(F('a.landxml', '<kml/>')), /LandXML/);
});

await group('structural decks: Abaqus, ANSYS CDB, LS-DYNA', async () => {
  const nodeLines = V.map((v, i) => `${i + 1}, ${v.join(', ')}`).join('\n');
  const hexInp = `*HEADING\n cube\n** comment line\n*NODE, NSET=ALL\n${nodeLines}\n*ELEMENT, TYPE=C3D8R, ELSET=BLOCK\n1, 1, 2, 3, 4, 5, 6, 7, 8\n*SOLID SECTION, ELSET=BLOCK, MATERIAL=STEEL\n1.,\n*END STEP\n`;
  const ai = await imp('cube.inp', hexInp);
  isCube(ai, 'Abaqus C3D8R');
  ok(ai.format === 'Abaqus input deck' && ai.pathway === 'mesh' && ai.stats.cellTypes.hex === 1 && ai.stats.nodes === 8, 'Abaqus: one hexahedron, 8 nodes');
  isCube(await imp('tets.inp', `*Node\n${nodeLines}\n*Element, type=C3D4\n${TETS.map((t, i) => `${i + 1}, ${t.map((k) => k + 1).join(', ')}`).join('\n')}\n`), 'Abaqus C3D4');
  const sh = await imp('shell.inp', `*NODE\n${nodeLines}\n*ELEMENT, TYPE=S4R\n${Q.map((q, i) => `${i + 1}, ${q.map((k) => k + 1).join(', ')}`).join('\n')}\n*ELEMENT, TYPE=SPRINGA\n90, 1, 2\n`);
  isCube(sh, 'Abaqus S4R shells');
  ok(sh.stats.surfaceOnly && sh.warnings.some((w) => /1 × SPRINGA/.test(w)), 'Abaqus: shell surface, unsupported SPRINGA reported');
  const t10 = await imp('tet10.inp', '*NODE\n1,0,0,0\n2,1,0,0\n3,0,1,0\n4,0,0,1\n5,.5,0,0\n6,.5,.5,0\n7,0,.5,0\n8,0,0,.5\n9,.5,0,.5\n10,0,.5,.5\n*ELEMENT, TYPE=C3D10\n1, 1,2,3,4,5,6,7,\n8,9,10\n*ELEMENT, TYPE=S8R\n2, 1,2,3,4,5,6,7,8\n');
  ok(t10.kind === 'mesh' && t10.stats.cellTypes.tet === 1 && Math.abs(dimensions(t10).volume - 1 / 6) < 1e-9, 'Abaqus C3D10 over two lines reduced to its corner tetrahedron');
  const pipe = await imp('pipe.inp', '*NODE\n1, 0., 0., 0.\n2, 3., 0., 0.\n3, 3., 0., 4.\n4, 3., 5., 4.\n*ELEMENT, TYPE=PIPE31, ELSET=RISER\n1, 1, 2\n2, 2, 3\n3, 3, 4\n*BEAM SECTION, SECTION=PIPE, ELSET=RISER, MATERIAL=X65\n0.1, 0.01\n');
  ok(pipe.kind === 'polylines' && pipe.polylines.length === 1 && pipe.polylines[0].z.join() === '0,0,4,4' && pipe.bbox.max.join() === '3,5,4' && pipe.stats.dimension === 1, 'Abaqus PIPE31 only → one 3-D polyline');
  near(dimensions(pipe).length, 12, 1e-12, 'Abaqus pipe string length 3 + 4 + 5 in 3-D');
  const b32 = await imp('b32.inp', '*NODE\n1,0,0,0\n2,1,0,0\n3,2,0,1\n*ELEMENT,TYPE=B32\n1,1,2,3\n');
  ok(b32.kind === 'polylines' && b32.polylines[0].x.join() === '0,1,2' && b32.polylines[0].z.join() === '0,0,1', 'Abaqus B32: mid node kept in order');
  const parts = await imp('parts.inp', `*PART, NAME=A\n*NODE\n${nodeLines}\n*ELEMENT, TYPE=C3D8\n1,1,2,3,4,5,6,7,8\n*END PART\n*PART, NAME=B\n*NODE\n${V.map((v, i) => `${i + 1}, ${v[0] + 5}, ${v[1]}, ${v[2]}`).join('\n')}\n*ELEMENT, TYPE=C3D8\n1,1,2,3,4,5,6,7,8\n*END PART\n*ASSEMBLY, NAME=ASM\n*INSTANCE, NAME=A-1, PART=A\n*END INSTANCE\n*INSTANCE, NAME=B-1, PART=B\n10., 0., 0.\n*END INSTANCE\n*END ASSEMBLY\n*NGEN\n1,8\n*INCLUDE, INPUT=more.inp\n`);
  ok(parts.count === 24 && parts.bbox.max[0] === 6 && parts.stats.parts === 2 && parts.warnings.some((w) => /instance/.test(w)) && parts.warnings.some((w) => /\*NGEN/.test(w)) && parts.warnings.some((w) => /INCLUDE/.test(w)), 'Abaqus parts keep separate node numbering; instance placement, *NGEN and *INCLUDE are reported as not applied');
  await throws('Abaqus deck without nodes', () => importGeometry(F('a.inp', '*HEADING\n*MATERIAL, NAME=A\n')), /NODE/);
  await throws('.inp that is not Abaqus', () => importGeometry(F('a.inp', 'hello\nworld\n')), /Abaqus/);

  // ANSYS CDB
  const i9 = (...a) => a.map((v) => String(v).padStart(9)).join(''), e21 = (v) => v.toExponential(13).replace(/e([+-])(\d+)$/, (m, s, d) => 'E' + s + d.padStart(3, '0')).padStart(21);
  const nblock = `NBLOCK,6,SOLID,       8,       8\n(3i9,6e21.13e3)\n${V.map((v, i) => i9(i + 1, 0, 0) + v.filter((_, k) => k < 2 || v[2] !== 0).map(e21).join('')).join('\n')}\nN,R5.3,LOC,       -1,\n`;
  const eline = (type, num, nodes) => i9(1, type, 1, 1, 0, 0, 0, 0, nodes.length, 0, num, ...nodes.slice(0, 8)) + (nodes.length > 8 ? '\n' + i9(...nodes.slice(8)) : '');
  const cdbHex = `/COM,ANSYS RELEASE\n/PREP7\nET,       1,185\nET,2,SHELL181\n${nblock}EBLOCK,19,SOLID,       1,       1\n(19i9)\n${eline(1, 1, [1, 2, 3, 4, 5, 6, 7, 8])}\n       -1\nFINISH\n`;
  const ch = await imp('cube.cdb', cdbHex);
  isCube(ch, 'ANSYS CDB SOLID185');
  ok(ch.format === 'ANSYS CDB archive' && ch.stats.cellTypes.hex === 1, 'CDB: NBLOCK with omitted trailing zeros + EBLOCK');
  const cdbTet = `ET,1,185\n${nblock}EBLOCK,19,SOLID,6,6\n(19i9)\n${TETS.map((t, i) => eline(1, i + 1, [t[0], t[1], t[2], t[2], t[3], t[3], t[3], t[3]].map((k) => k + 1))).join('\n')}\n-1\n`;
  const ct = await imp('tets.cdb', cdbTet);
  isCube(ct, 'ANSYS CDB degenerate tetrahedra');
  ok(ct.stats.cellTypes.tet === 6, 'CDB: collapsed hexahedra recognised as tetrahedra');
  const cs = await imp('shell.cdb', `ET,7,181\n${nblock}EBLOCK,19,SOLID,6,6\n(19i9)\n${Q.map((q, i) => eline(7, i + 1, q.map((k) => k + 1))).join('\n')}\n-1\nEBLOCK,10,,2,2\n(15i9)\n${i9(1, 2, 1, 1, 7, 1, 1, 2, 3, 4)}\n-1\n`);
  isCube(cs, 'ANSYS CDB SHELL181');
  ok(cs.warnings.some((w) => /without the SOLID key/.test(w)), 'CDB: non-SOLID EBLOCK reported as not read');
  const c187 = await imp('t10.cdb', `ET,1,187\nET,3,188\nNBLOCK,6,SOLID,10,10\n(3i9,6e21.13e3)\n${[[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [0.5, 0, 0], [0.5, 0.5, 0], [0, 0.5, 0], [0, 0, 0.5], [0.5, 0, 0.5], [0, 0.5, 0.5]].map((v, i) => i9(i + 1, 0, 0) + v.map(e21).join('')).join('\n')}\nN,R5.3,LOC,-1,\nEBLOCK,19,SOLID,2,2\n(19i9)\n${eline(1, 1, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10])}\n${eline(3, 2, [1, 2, 3])}\n-1\n`);
  ok(c187.kind === 'mesh' && c187.stats.cellTypes.tet === 1 && c187.stats.cellTypes.line === 1 && Math.abs(dimensions(c187).volume - 1 / 6) < 1e-9, 'CDB: SOLID187 on a continuation line, BEAM188 as a line element');
  await throws('CDB with unblocked nodes', () => importGeometry(F('a.cdb', 'N,1,0,0,0\nN,2,1,0,0\nEN,1,1,2\n')), /unblocked|CDWRITE/);

  // LS-DYNA
  const i8 = (...a) => a.map((v) => String(v).padStart(8)).join(''), e16 = (v) => v.toFixed(6).padStart(16);
  const dnodes = `*NODE\n$#   nid               x               y               z      tc      rc\n${V.map((v, i) => i8(i + 1) + v.map(e16).join('') + i8(0, 0)).join('\n')}\n`;
  const dk = await imp('cube.k', `*KEYWORD\n*TITLE\ncube\n${dnodes}*ELEMENT_SOLID\n$#   eid     pid      n1      n2      n3      n4      n5      n6      n7      n8\n${i8(1, 1, 1, 2, 3, 4, 5, 6, 7, 8)}\n*END\n`);
  isCube(dk, 'LS-DYNA *ELEMENT_SOLID (one card)');
  ok(dk.format === 'LS-DYNA keyword deck' && dk.stats.cellTypes.hex === 1, 'LS-DYNA: fixed-format nodes and solid');
  const d2 = await imp('tets.key', `*KEYWORD\n${dnodes}*ELEMENT_SOLID\n${TETS.map((t, i) => `${i8(i + 1, 1)}\n${i8(...[t[0], t[1], t[2], t[3], t[3], t[3], t[3], t[3]].map((k) => k + 1))}`).join('\n')}\n*END\n`);
  isCube(d2, 'LS-DYNA *ELEMENT_SOLID (two cards, degenerate tetrahedra)');
  const d3 = await imp('shell.dyn', `*KEYWORD\n*NODE\n${V.map((v, i) => `${i + 1},${v.join(',')}`).join('\n')}\n*ELEMENT_SHELL_THICKNESS\n${Q.map((q, i) => `${i8(i + 1, 1, ...q.map((k) => k + 1))}\n${[0.01, 0.01, 0.01, 0.01].map(e16).join('')}`).join('\n')}\n*ELEMENT_MASS\n${i8(99, 1)}             1.0\n*INCLUDE\nother.k\n*END\n`);
  isCube(d3, 'LS-DYNA *ELEMENT_SHELL_THICKNESS (comma nodes)');
  ok(d3.warnings.some((w) => /option cards/.test(w)) && d3.warnings.some((w) => /ELEMENT_MASS/.test(w)) && d3.warnings.some((w) => /INCLUDE/.test(w)), 'LS-DYNA: thickness cards skipped, *ELEMENT_MASS and *INCLUDE reported');
  const db = await imp('pipe.k', `*KEYWORD\n*NODE\n${[[0, 0, 0], [3, 0, 0], [3, 0, 4]].map((v, i) => i8(i + 1) + v.map(e16).join('')).join('\n')}\n*ELEMENT_BEAM\n${i8(1, 1, 1, 2, 3)}\n${i8(2, 1, 2, 3, 1)}\n*END\n`);
  ok(db.kind === 'polylines' && db.polylines[0].z.join() === '0,0,4' && Math.abs(dimensions(db).length - 7) < 1e-12, 'LS-DYNA *ELEMENT_BEAM only → 3-D polyline of length 7');
  const dl = await imp('long.k', `*KEYWORD LONG=Y\n*NODE\n${V.map((v, i) => String(i + 1).padStart(20) + v.map((c) => c.toFixed(8).padStart(20)).join('')).join('\n')}\n*ELEMENT_SOLID\n${[1, 1, 1, 2, 3, 4, 5, 6, 7, 8].map((v) => String(v).padStart(20)).join('')}\n*END\n`);
  isCube(dl, 'LS-DYNA LONG=Y format');
  await throws('LS-DYNA deck without nodes', () => importGeometry(F('a.k', '*KEYWORD\n*PART\nx\n*END\n')), /NODE/);
});

await group('well deviation surveys', async () => {
  const vert = await imp('vertical.csv', 'MD,Inclination,Azimuth\n0,0,0\n100,0,0\n250,0,0\n');
  ok(vert.kind === 'polylines' && vert.pathway === 'well' && vert.format === 'Well deviation survey (ASCII)' && vert.stats.wellSurvey === true && vert.polylines.length === 1, 'CSV with MD / Inclination / Azimuth detected as a well survey');
  ok(vert.polylines[0].z.join() === '0,-100,-250' && vert.polylines[0].x.every((v) => v === 0) && vert.stats.md === 250 && vert.stats.tvd === 250 && vert.stats.maxInclination === 0 && vert.stats.maxDLS === 0 && vert.survey.tvd[2] === 250, 'vertical well: TVD = MD, z negative down, no dog-leg');
  // quarter-circle build to horizontal heading east, radius R: closed form TVD = R, east = R, DLS = 30 / R rad
  const R = 300, st = Array.from({ length: 10 }, (_, k) => [(R * k * 10 * Math.PI) / 180, k * 10, 90]);
  const build = await imp('build.dev', `# WELL: A-1\n# UNITS: METRES\n#    MD        INCL      AZIM\n${st.map((r) => r.map((v) => v.toFixed(9).padStart(16)).join('')).join('\n')}\n`);
  const pb = build.polylines[0], e = pb.x.length - 1;
  near(pb.x[e], R, 1e-6, 'minimum curvature: 90° build, east departure = R');
  near(pb.z[e], -R, 1e-6, 'minimum curvature: 90° build, TVD = R');
  ok(Math.abs(pb.y[e]) < 1e-9 && build.format === 'Well deviation survey (ASCII)' && build.stats.stations === 10 && build.stats.maxInclination === 90 && Math.abs(build.stats.maxDLS - (30 / R) * (180 / Math.PI)) < 1e-6 && build.survey.dls[0] === 0 && build.survey.md.length === 10, `.dev with commented header: stats ${JSON.stringify(build.stats)}`);
  const ft = await imp('feet.csv', 'Measured Depth (ft),Incl (deg),Azim (deg),TVD (ft)\n0,0,0,0\n1000,0,0,1000\n2000,0,0,2000\n');
  ok(ft.stats.wellSurvey && Math.abs(ft.stats.md - 609.6) < 1e-9 && Math.abs(ft.polylines[0].z[2] + 609.6) < 1e-9 && ft.warnings.some((w) => /metres/.test(w)), 'multi-word headers with feet: depths converted to metres');
  const bare = await imp('plain.survey', '0 0 0\n500 0 0\n1000 30 45\n');
  ok(bare.stats.wellSurvey && bare.polylines[0].x.length === 3 && bare.polylines[0].x[2] > 0 && Math.abs(bare.polylines[0].x[2] - bare.polylines[0].y[2]) < 1e-9 && bare.warnings.some((w) => /No column header/.test(w)), 'header-less .survey read as MD, inclination, azimuth (azimuth 45°: east = north)');
  const mt = await imp('mdtvd.txt', 'MD TVD\n0 0\n100 100\n200 180\n300 240\n');
  ok(mt.stats.wellSurvey && mt.polylines[0].z.join() === '0,-100,-180,-240' && Math.abs(mt.polylines[0].x[2] - 60) < 1e-9 && Math.abs(mt.polylines[0].x[3] - 140) < 1e-9 && mt.polylines[0].y.every((v) => v === 0) && mt.warnings.some((w) => /due east/.test(w)), 'MD + TVD: departure from √(ΔMD² − ΔTVD²), laid out east, with a warning');
  near(mt.survey.inc[3], Math.acos(0.6) * 180 / Math.PI, 1e-9, 'MD + TVD: interval inclination acos(ΔTVD / ΔMD)');
  const ne = await imp('offsets.wbt', 'MD,TVD,NS,EW\n0,0,0,0\n100,100,0,0\n200,180,60,0\n300,240,60,80\n');
  ok(ne.polylines[0].y.join() === '0,0,60,60' && ne.polylines[0].x.join() === '0,0,0,80' && Math.abs(ne.survey.azi[3] - 90) < 1e-9 && Math.abs(ne.survey.azi[2]) < 1e-9, 'MD + TVD + NS / EW offsets used directly; azimuth derived');
  const js = await imp('survey.json', JSON.stringify([{ md: 0, inc: 0, azi: 0 }, { md: 500, inc: 0, azi: 0 }, { md: 800, inc: 30, azi: 180 }]));
  ok(js.kind === 'polylines' && js.stats.wellSurvey && js.polylines[0].y[2] < 0 && js.format === 'JSON data' && js.pathway === 'well', 'JSON survey records (md, inc, azi) → well path heading south');
  const dropped = await imp('gaps.dev', 'MD INC AZI\n0 0 0\n100 0 0\n100 1 0\n90 1 0\n200 2 0\nx y z\n');
  ok(dropped.stats.stations === 3 && dropped.warnings.some((w) => /dropped/.test(w)), 'non-increasing and non-numeric survey rows dropped with a warning');
  await throws('survey with one station', () => importGeometry(F('a.dev', 'MD INC AZI\n0 0 0\n')), /two stations/);
  await throws('survey with impossible inclination', () => importGeometry(F('a.dev', '0 0 0\n10 500 0\n')), /between 0 and 180/);
});

await group('pipeline and inspection tables', async () => {
  const prof = await imp('profile.csv', 'KP (km),Elevation (m),Comment\n0,-1000,PLEM\n1.5,-990,\n3,-1005,low point\n4,-200,riser base\n');
  ok(prof.kind === 'table' && prof.stats.role === 'profile' && prof.headers.join() === 'KP (km),Elevation (m),Comment' && prof.records[2]['Elevation (m)'] === -1005 && prof.records[2].Comment === 'low point' && prof.records.length === 4, 'chainage + elevation table keeps headers / records, role "profile"');
  const ws = await imp('profile.txt', 'chainage   depth\n0     100\n500   120\n900   95\n');
  ok(ws.kind === 'table' && ws.stats.role === 'profile' && ws.headers.join() === 'chainage,depth' && ws.records[1].depth === 120 && ws.records[2].chainage === 900, 'whitespace-separated station / depth table');
  const rt = await imp('route.csv', 'Easting,Northing,Elevation\n500000,6500000,-300\n500400,6500300,-310\n500800,6500600,-305\n501200,6500900,-250\n');
  ok(rt.kind === 'points' && rt.count === 4 && rt.stats.role === 'route3d' && rt.points[2] === -300 && rt.points[9] === 501200, 'easting / northing / elevation in route order → ordered points, role "route3d"');
  const ll = await imp('route_ll.csv', 'Latitude,Longitude,Depth (m)\n60.0,2.0,300\n60.01,2.02,320\n60.02,2.05,310\n60.03,2.06,150\n60.05,2.07,20\n');
  ok(ll.kind === 'points' && ll.geographic && ll.stats.role === 'route3d' && ll.points[0] === 2 && ll.points[1] === 60 && ll.points[2] === -300 && ll.warnings.some((w) => /single line/.test(w)) && ll.warnings.some((w) => /depth/.test(w)), 'latitude / longitude / depth along a line → geographic route with negative elevations');
  const kp = await imp('route_kp.csv', 'KP,x,y,z\n0,0,0,-100\n1,800,600,-90\n2,1600,1200,-95\n');
  ok(kp.kind === 'table' && kp.stats.role === 'route3d' && kp.records[1].KP === 1 && kp.records[2].z === -95, 'x / y / z with a KP column stays a table, role "route3d"');
  const wt = await imp('wt.csv', 'x (m),Clock Position,Wall Thickness (mm)\n0,12,12.7\n0,3,12.5\n0,6,9.8\n0,9,12.6\n0.5,12,12.7\n0.5,3,12.4\n0.5,6,10.1\n0.5,9,12.6\n');
  ok(wt.kind === 'table' && wt.stats.role === 'thicknessMap' && wt.headers[2] === 'Wall Thickness (mm)' && wt.records[2]['Wall Thickness (mm)'] === 9.8 && wt.records.length === 8, 'wall-thickness map (x, clock position, thickness), multi-word headers');
  const co = await imp('corrosion.csv', 'chainage,theta,pit depth\n10,0,0.5\n10,90,1.5\n20,0,0.2\n20,90,0.1\n');
  ok(co.stats.role === 'thicknessMap' && co.records[1]['pit depth'] === 1.5, 'corrosion map (chainage, theta, pit depth)');
  const dp = await imp('deposit.csv', 'x,theta,t,delta\n0,0,0,0\n0,180,0,0\n100,0,0,0\n100,180,0,0\n0,0,3600,0.4\n0,180,3600,1.2\n100,0,3600,0.1\n100,180,3600,2.5\n');
  ok(dp.kind === 'table' && dp.stats.role === 'depositMap' && dp.records.length === 8 && dp.records[7].delta === 2.5 && dp.records[7].t === 3600, 'time-resolved deposit map (x, θ, t, δ)');
  const pa = await imp('particles.csv', 'x,y,z,diameter,velocity\n0,0,0,0.001,1.5\n1,0.1,0,0.002,1.4\n2,0.3,0.1,0.004,1.2\n0.5,0.2,0.05,0.003,1.1\n');
  ok(pa.kind === 'points' && pa.stats.role === 'particles' && pa.count === 4 && pa.attributes.diameter[2] === 0.004 && pa.attributes.velocity[3] === 1.1 && Object.keys(pa.attributes).length === 2, 'particle cloud: positions as points, size and velocity kept as attributes');
  const ts = await imp('ops.csv', 'time,choke opening,pressure\n0,0.5,120\n60,0.55,118\n120,0.6,115\n');
  ok(ts.kind === 'table' && ts.stats.role === 'timeseries' && ts.records[2].pressure === 115, 'operating time series, role "timeseries"');
  const tsd = await imp('ops_dates.csv', 'timestamp,valve,state\n2024-01-01T00:00,V-1,open\n2024-01-01T01:00,V-1,closed\n');
  ok(tsd.kind === 'table' && tsd.stats.role === 'timeseries' && tsd.records[1].state === 'closed', 'time series with text cells');
  const jr = await imp('wt.json', JSON.stringify([{ x: 0, theta: 0, wt: 12 }, { x: 0, theta: 180, wt: 11 }]));
  ok(jr.kind === 'table' && jr.stats.role === 'thicknessMap', 'JSON record table gets the same role');
  const plain = await imp('plain.csv', 'name,value\nrecovery,0.45\n');
  ok(plain.kind === 'table' && plain.stats.role === undefined, 'a parameter table has no role');
});

await group('SEG-Y', async () => {
  const text = (ebc) => { const lines = Array.from({ length: 40 }, (_, k) => `C${String(k + 1).padStart(2)} ${k === 0 ? 'CLIENT HYDRASLUG TEST LINE 2D-01' : ''}`.padEnd(80)).join(''), u = new Uint8Array(3200); for (let i = 0; i < 3200; i++) { const c = lines.charCodeAt(i); u[i] = !ebc ? c : c === 32 ? 0x40 : c === 45 ? 0x60 : c >= 48 && c <= 57 ? 0xf0 + c - 48 : c >= 65 && c <= 73 ? 0xc1 + c - 65 : c >= 74 && c <= 82 ? 0xd1 + c - 74 : c >= 83 && c <= 90 ? 0xe2 + c - 83 : 0x40; } return u; };
  const segy = ({ fmt = 5, ebc = true, wd = 0, cdp = false, n = 5, ns = 50, coords = true } = {}) => {
    const bh = new Uint8Array(400), bv = new DataView(bh.buffer), bps = { 1: 4, 5: 4, 3: 2 }[fmt], parts = [text(ebc), bh];
    bv.setUint16(16, 4000); bv.setUint16(20, ns); bv.setUint16(24, fmt);
    for (let t = 0; t < n; t++) {
      const th = new Uint8Array(240), tv = new DataView(th.buffer), sm = new Uint8Array(ns * bps), sv = new DataView(sm.buffer), k = 10 + t;
      tv.setInt32(0, t + 1); tv.setInt16(68, 1); tv.setInt16(70, -10); tv.setInt32(60, wd); tv.setInt16(114, ns); tv.setInt16(116, 4000);
      if (coords) { tv.setInt32(cdp ? 180 : 72, 5000000 + 1000 * t); tv.setInt32(cdp ? 184 : 76, 65000000); }
      for (let s = 0; s < ns; s++) { const a = s === k ? 1 : s === k + 1 ? -0.4 : s === 3 ? 0.05 : 0; if (fmt === 5) sv.setFloat32(4 * s, a); else if (fmt === 3) sv.setInt16(2 * s, Math.round(a * 1000)); else sv.setUint32(4 * s, a === 1 ? 0x41100000 : a === -0.4 ? 0xc0666666 : a === 0.05 ? 0x3fcccccd : 0); }
      parts.push(th, sm);
    }
    return cat(...parts);
  };
  const g = await imp('line.sgy', segy());
  ok(g.kind === 'points' && g.format === 'SEG-Y seismic' && g.count === 5 && g.stats.traces === 5 && g.stats.samplesPerTrace === 50 && g.stats.sampleInterval === 0.004 && g.stats.sampleFormat === 'IEEE float32' && g.stats.textHeader === 'EBCDIC', `SEG-Y headers: ${JSON.stringify(g.stats)}`);
  ok(g.points[0] === 500000 && g.points[1] === 6500000 && g.points[3] === 500100 && /CLIENT HYDRASLUG TEST LINE 2D-01/.test(g.textHeader) && g.textHeader.startsWith('C 1'), 'SEG-Y: source coordinates with scalar −10 applied, EBCDIC text header decoded');
  ok(g.seabed && g.seabed.twt.every((v, t) => Math.abs(v - (10 + t) * 0.004) < 1e-12) && g.seabed.depth.every((v, t) => Math.abs(v - 3 * (10 + t)) < 1e-9) && Math.abs(g.seabed.distance[4] - 400) < 1e-9, 'SEG-Y: first strong arrival picked per trace, 1500 m/s → depth 3 m per sample');
  ok(Math.abs(g.points[2] + 30) < 1e-9 && Math.abs(g.points[14] + 42) < 1e-9 && g.warnings.some((w) => /ESTIMATE/.test(w)) && g.warnings.some((w) => /no water depth/.test(w)), 'SEG-Y: without header water depths z is the estimate, labelled as such');
  const ibm = await imp('ibm.segy', segy({ fmt: 1, ebc: false, wd: 1500, cdp: true }));
  ok(ibm.stats.sampleFormat === 'IBM float' && ibm.stats.textHeader === 'ASCII' && ibm.points[2] === -1500 && ibm.points[0] === 500000 && Math.abs(ibm.seabed.depth[2] - 36) < 1e-9 && ibm.warnings.some((w) => /CDP/.test(w)) && /CLIENT/.test(ibm.textHeader), 'SEG-Y: IBM floats decoded, ASCII header, CDP coordinates, z from water depth at source');
  const i16s = await imp('int16.sgy', segy({ fmt: 3 }));
  ok(i16s.count === 5 && !i16s.seabed && i16s.points[2] === 0 && i16s.warnings.some((w) => /not decoded/.test(w)), 'SEG-Y: int16 samples give positions only');
  const nopos = await imp('nopos.sgy', segy({ coords: false }));
  ok(nopos.points[3] === 1 && nopos.points[6] === 2 && nopos.warnings.some((w) => /no source or CDP/.test(w)), 'SEG-Y without coordinates: traces placed by index, with a warning');
  const cutf = segy(), tr = await imp('cut.sgy', cutf.subarray(0, cutf.length - 100));
  ok(tr.count === 4 && tr.warnings.some((w) => /ends inside trace 5/.test(w)), 'SEG-Y: truncated last trace reported');
  await throws('SEG-Y that is too short', () => importGeometry(F('a.sgy', new Uint8Array(1000))), /SEG-Y/);
  await throws('SEG-Y with an invalid format code', () => importGeometry(F('a.sgy', new Uint8Array(5000))), /format code/);
});

await group('MATLAB MAT-files', async () => {
  const pad8 = (u) => cat(u, new Uint8Array((8 - (u.length % 8)) % 8));
  const el = (type, data) => pad8(cat(le('u32', [type, data.length]), data));
  const small = (type, data) => cat(le('u32', [(data.length << 16) | type]), data, new Uint8Array(4 - data.length));
  const mx = (name, dims, vals, { cls = 6, dt = 9, enc = 'f64', flags = 0 } = {}) => el(14, cat(el(6, le('u32', [cls | flags, 0])), el(5, le('i32', dims)), name.length <= 4 ? small(1, enc8(name)) : el(1, enc8(name)), el(dt, le(enc, vals))));
  const enc8 = (s) => new TextEncoder().encode(s);
  const head = (txt = 'MATLAB 5.0 MAT-file, Platform: GLNXA64, Created on: Mon Jan  1 00:00:00 2024') => cat(enc8(txt.padEnd(116)), new Uint8Array(8), le('u16', [0x0100]), enc8('IM'));
  const z = (u) => { const d = new Uint8Array(zlib.deflateSync(u)); return cat(le('u32', [15, d.length]), d); };
  const md = [0, 500, 1000, 1500], inc = [0, 0, 20, 40], azi = [10, 10, 10, 10];
  const sv = await imp('survey.mat', cat(head(), z(mx('md', [4, 1], md)), z(mx('inc', [4, 1], inc)), z(mx('azi', [1, 4], azi)), mx('note', [1, 3], [72, 105, 33], { cls: 4, dt: 4, enc: 'u16' })));
  ok(sv.kind === 'polylines' && sv.stats.wellSurvey && sv.format === 'MATLAB MAT-file' && sv.stats.md === 1500 && sv.survey.inc[3] === 40 && sv.stats.variables.join() === 'md,inc,azi', 'MAT v7 (compressed): equal-length vectors md, inc, azi combined into a survey');
  ok(sv.warnings.some((w) => /3 vectors of length 4/.test(w)) && sv.warnings.some((w) => /note \(char\)/.test(w)), 'MAT: vector grouping and skipped character array reported');
  const P = [[0, -100, 7], [10, -98, 8], [20, -97, 9], [30, -99, 10], [40, -90, 11]], colMajor = [0, 1, 2].flatMap((c) => P.map((r) => r[c]));
  const tb = await imp('matrix.mat', cat(head(), mx('P', [5, 3], colMajor)));
  ok(tb.kind === 'table' && tb.headers.join() === 'P_1,P_2,P_3' && tb.records.length === 5 && tb.records[3].P_2 === -99 && tb.records[4].P_3 === 11 && tb.stats.variable === 'P' && tb.stats.shape.join() === '5,3', 'MAT v6 (uncompressed): 5 × 3 matrix as a table (column-major order honoured)');
  const nr = 12, nc = 20, dem = []; for (let c = 0; c < nc; c++) for (let r = 0; r < nr; r++) dem.push(-100 + 0.5 * c + 0.25 * r);
  const gr = await imp('dem.mat', cat(head(), z(mx('bathy', [nr, nc], dem, { cls: 7, dt: 7, enc: 'f32' })), mx('k', [1, 1], [3])));
  ok(gr.kind === 'grid' && gr.grid.x.length === nc && gr.grid.y.length === nr && Math.abs(gr.grid.z[2][4] - (-100 + 2 + 0.5)) < 1e-6 && gr.stats.variable === 'bathy' && gr.warnings.some((w) => /largest/.test(w)), 'MAT: single-precision 12 × 20 matrix as a grid; largest variable chosen');
  const pick = await imp('pick.mat', cat(head(), z(mx('bathy', [nr, nc], dem)), mx('xz', [2, 9], Array.from({ length: 18 }, (_, i) => (i % 2 ? -50 - i : i * 5)))), { variable: 'xz' });
  ok(pick.kind === 'table' && pick.records.length === 9 && pick.headers.length === 2 && pick.records[1].xz_1 === 10 && pick.records[1].xz_2 === -53 && pick.warnings.some((w) => /rows as table columns/.test(w)), 'MAT: opts.variable selects an array; a 2 × N matrix is transposed into columns');
  const vx = await imp('vox.mat', cat(head(), mx('lab', [2, 2, 2], [1, 0, 0, 0, 1, 1, 0, 1], { cls: 9, dt: 2, enc: 'u8' })));
  ok(vx.kind === 'voxels' && vx.voxels.nx === 2 && vx.voxels.nz === 2 && vx.stats.porosity === 0.5, 'MAT: uint8 3-D array as voxels');
  await throws('MAT v7.3 (HDF5)', () => importGeometry(F('a.mat', cat(head('MATLAB 7.3 MAT-file, Platform: GLNXA64, Created on: Mon Jan  1 00:00:00 2024 HDF5 schema 1.00 .'), new Uint8Array(600)))), /7\.3.*HDF5|HDF5/);
  await throws('MAT v4 or other data', () => importGeometry(F('a.mat', new Uint8Array(400))), /Level 5/);
  await throws('MAT with only text', () => importGeometry(F('a.mat', cat(head(), mx('note', [1, 3], [72, 105, 33], { cls: 4, dt: 4, enc: 'u16' })))), /no numeric arrays/);
  await throws('MAT with a missing variable', () => importGeometry(F('a.mat', cat(head(), mx('P', [5, 3], colMajor))), { variable: 'nope' }), /no numeric variable/);
});

await group('DEM / DTM text models and IFCZIP', async () => {
  const asc = await imp('seabed.dem', 'ncols 3\nnrows 2\nxllcorner 1000\nyllcorner 2000\ncellsize 10\nNODATA_value -9999\n-50 -51 -52\n-60 -61 -9999\n');
  ok(asc.kind === 'grid' && asc.format === 'ESRI ASCII grid' && asc.grid.x.length === 3 && asc.grid.z[0][0] === -60 && asc.grid.z[1][2] === -52 && asc.grid.nodata === 1, '.dem holding an ESRI ASCII grid');
  const fld = (s, w) => String(s).padStart(w), A = new Array(1024).fill(' ');
  const putA = (o, s) => { for (let k = 0; k < s.length; k++) A[o + k] = s[k]; };
  putA(0, 'SYNTHETIC QUAD'); putA(156, fld(1, 6)); putA(162, fld(31, 6)); putA(528, fld(2, 6)); putA(534, fld(2, 6)); putA(540, fld(4, 6)); putA(816, '0.300000E+020.300000E+020.100000E+01'); putA(852, fld(1, 6) + fld(3, 6));
  const prof = (c, y0, z) => `${fld(1, 6)}${fld(c, 6)}${fld(z.length, 6)}${fld(1, 6)}   0.${String(5000 + 300 * (c - 1)).padEnd(15, '0')}D+03   0.${String(y0).padEnd(15, '0')}D+04   0.000000000000000D+00   0.100000000000000D+02   0.900000000000000D+02${z.map((v) => fld(v, 6)).join('')}`.padEnd(1024);
  const us = await imp('quad.dem', A.join('') + prof(1, 1000, [10, 20, 30, 40]) + prof(2, 1030, [21, 31, 41]) + prof(3, 1000, [12, 22, -32767, 42]));
  ok(us.kind === 'grid' && us.format === 'DEM / DTM / DSM elevation model' && us.stats.demFormat === 'USGS ASCII DEM' && us.grid.x.join() === '500,530,560' && us.grid.y.join() === '1000,1030,1060,1090', `USGS ASCII DEM: profiles placed on a regular grid (${us.grid.x.join()} × ${us.grid.y.join()})`);
  ok(us.grid.z[0][0] === 10 && us.grid.z[3][0] === 40 && Number.isNaN(us.grid.z[0][1]) && us.grid.z[1][1] === 21 && us.grid.z[3][1] === 41 && Number.isNaN(us.grid.z[2][2]) && us.grid.z[3][2] === 42 && us.bbox.max[2] === 42, 'USGS ASCII DEM: staggered profile start and void value handled');
  const xyz = await imp('scan.dtm', '0 0 -10\n10 0 -11\n0 10 -12\n10 10 -13\n5 5 -11.5\n');
  ok(xyz.kind === 'points' && xyz.count === 5 && xyz.pathway === 'gis', '.dtm holding x y z rows');
  await throws('binary DEM', () => importGeometry(F('a.dem', new Uint8Array(4000).map((_, i) => (i * 37) & 255))), /binary.*gdal_translate/);
  await throws('DEM that is a plain table', () => importGeometry(F('a.dsm', 'name,value\na,1\n')), /neither/);
});

await group('line generators', async () => {
  const D = 1000, th = (78 * Math.PI) / 180, t = Math.tan(th), a = D / (Math.sqrt(1 + t * t) - 1), len3 = (g) => { const p = g.polylines[0]; let l = 0; for (let i = 1; i < p.x.length; i++) l += Math.hypot(p.x[i] - p.x[i - 1], p.y[i] - p.y[i - 1], p.z[i] - p.z[i - 1]); return l; };
  const cr = generate({ type: 'catenary riser', waterDepth: D, topAngle: 12, flowline: 0, n: 400 });
  ok(cr.kind === 'polylines' && cr.pathway === 'procedural' && cr.format === 'Procedural: catenary' && cr.polylines[0].z[0] === -D && cr.polylines[0].z[cr.polylines[0].z.length - 1] === 0 && cr.polylines[0].y.every((v) => v === 0), 'catenary riser: seabed to surface, x–z plane');
  near(cr.stats.suspendedLength, a * t, 1e-9, 'catenary riser: suspended length = a · tan(θ_top) (closed form)');
  near(len3(cr), a * t, 1e-3 * a * t, 'catenary riser: polyline length matches the closed form');
  near(cr.stats.horizontalLength, a * Math.asinh(t), 1e-9, 'catenary riser: horizontal reach = a · asinh(tan θ_top)');
  const pz = cr.polylines[0], e = pz.x.length - 1;
  near(Math.atan2(pz.z[e] - pz.z[e - 1], pz.x[e] - pz.x[e - 1]) * 180 / Math.PI, 78, 0.2, 'catenary riser: top angle 12° from vertical');
  const lw = generate({ type: 'lazy-wave riser', waterDepth: D, topAngle: 10, sagHeight: 100, hogHeight: 260, buoyancyRatio: 1.5, flowline: 200 });
  const zl = lw.polylines[0].z, ext = []; for (let i = 1; i + 1 < zl.length; i++) if ((zl[i] - zl[i - 1]) * (zl[i + 1] - zl[i]) < 0) ext.push(zl[i] + D);
  ok(lw.format === 'Procedural: lazywave' && zl[0] === -D && zl[1] === -D && lw.polylines[0].x[1] === 200 && zl[zl.length - 1] === 0 && ext.length === 2 && Math.abs(ext[0] - 260) < 1 && Math.abs(ext[1] - 100) < 1, `lazy-wave riser: lead-in on the seabed, hog at 260 m and sag at 100 m above it (${ext.map((v) => v.toFixed(1))})`);
  near(len3(lw), lw.stats.length, 2e-3 * lw.stats.length, 'lazy-wave riser: polyline length matches the sum of the catenary arcs');
  const fl = generate({ type: 'undulating flowline', length: 4000, amplitude: 15, wavelength: 400, slope: 1, waterDepth: 500, seed: 4 }), fl2 = generate({ type: 'flowline', length: 4000, amplitude: 15, wavelength: 400, slope: 1, waterDepth: 500, seed: 4 }), fl3 = generate({ type: 'flowline', length: 4000, amplitude: 15, wavelength: 400, slope: 1, waterDepth: 500, seed: 5 });
  const fp = fl.polylines[0], dev = fp.x.map((x, i) => fp.z[i] + 500 - x * Math.tan(Math.PI / 180));
  ok(fp.x[fp.x.length - 1] === 4000 && dev.every((v) => Math.abs(v) <= 15 + 1e-9) && Math.max(...dev) - Math.min(...dev) > 10 && fp.z.join() === fl2.polylines[0].z.join() && fp.z.join() !== fl3.polylines[0].z.join(), 'undulating flowline: hills bounded by the amplitude about the mean slope, reproducible per seed');
  const R = 30 / ((3 * Math.PI) / 180), wl = generate({ type: 'well trajectory (build-hold)', kickoff: 500, buildRate: 3, inclination: 45, md: 3000, azimuth: 90 });
  ok(wl.stats.wellSurvey && wl.survey.md[wl.survey.md.length - 1] === 3000 && wl.survey.md.includes(950) && wl.stats.maxInclination === 45 && Math.abs(wl.stats.maxDLS - 3) < 1e-9, 'build-and-hold well: stations to TD, end of build at 950 m MD, DLS = build rate');
  near(wl.stats.tvd, 500 + R * Math.sin(Math.PI / 4) + 2050 * Math.cos(Math.PI / 4), 1e-6, 'build-and-hold well: TVD = KOP + R sin α + hold · cos α (closed form)');
  near(wl.polylines[0].x[wl.polylines[0].x.length - 1], R * (1 - Math.cos(Math.PI / 4)) + 2050 * Math.sin(Math.PI / 4), 1e-6, 'build-and-hold well: departure = R (1 − cos α) + hold · sin α, due east');
  const jm = generate({ type: 'jumper (M-shape)', span: 30, height: 8, dip: 3, waterDepth: 1200 });
  ok(jm.polylines[0].x.length === 8 && jm.bbox.min[2] === -1200 && jm.bbox.max[2] === -1192 && Math.abs(len3(jm) - (2 * 8 + 2 * 3 + 30)) < 1e-9 && jm.stats.length === 52 && jm.polylines[0].z[3] === -1195, 'M-shaped jumper: length 2h + 2·dip + span');
  ok(generate({ type: 'jumper', dip: 0 }).polylines[0].x.length === 6 && Math.abs(dimensions(jm).length - 52) < 1e-9, 'jumper without a dip is an inverted U; dimensions() measures 3-D polyline length');
});

await group('suite links (geomlinks)', async () => {
  const { geometryLinks, profileTable, routeProfile } = await import('../js/core/geomlinks.js');
  const suite = (id, keys, extra = {}) => ({ id, inputs: [{ fields: keys.map((k) => (typeof k === 'string' ? { key: k, label: k } : k)) }], ...extra });
  const all = ['profile', 'length', 'waterDepth', 'riserHeight', 'id', 'wt', 'survey', 'wellDepth', 'wellMD', 'network', 'terrain', 'wtMap', 'minWt', 'depositMap', 'depositMax', 'structMesh', 'spanLength'].map((k) => (['profile', 'survey', 'network'].includes(k) ? { key: k, label: k, type: 'table' } : ['terrain', 'wtMap', 'depositMap', 'structMesh'].includes(k) ? { key: k, label: k, type: 'file' } : k));
  const get = (items, k) => (items.find((it) => it.key === k) || {}).value, keysOf = (items) => items.map((it) => it.key).sort().join();
  const riser = generate({ type: 'catenary riser', waterDepth: 1000, topAngle: 12, flowline: 500, n: 300 });
  const net = geometryLinks(suite('net', all), riser, {}), prof = get(net, 'profile');
  ok(Array.isArray(prof) && prof.length <= 60 && prof.length > 10 && prof[0].x === 0 && prof[0].z === -1000 && prof[prof.length - 1].z === 0 && Object.keys(prof[0]).join() === 'x,z', `net: riser profile offered as ≤ 60 rows of { x, z } (${prof.length})`);
  near(get(net, 'length'), riser.stats.length, 1e-3 * riser.stats.length, 'net: length along the pipe');
  ok(get(net, 'waterDepth') === 1000 && get(net, 'riserHeight') > 950 && get(net, 'riserHeight') < 1000 && keysOf(net) === 'length,profile,riserHeight,waterDepth' && net.every((it) => typeof it.from === 'string' && it.from.length > 5), `net: water depth and riser height (${get(net, 'riserHeight')}); every offer names its source`);
  ok(keysOf(geometryLinks(suite('flow', all), riser, {})) === 'length,profile,riserHeight,waterDepth' && keysOf(geometryLinks(suite('flow', ['length']), riser, {})) === 'length' && keysOf(geometryLinks(suite('solids', all), riser, {})) === 'length,riserHeight,waterDepth', 'flow takes the profile; only declared inputs are offered; solids gets no profile');
  const well = generate({ type: 'well', md: 2500 }), wo = geometryLinks(suite('net', all), well, {});
  ok(keysOf(wo) === 'survey,wellDepth,wellMD' && get(wo, 'survey').length === well.survey.md.length && Object.keys(get(wo, 'survey')[0]).join() === 'md,inc,azi' && get(wo, 'wellMD') === 2500 && Math.abs(get(wo, 'wellDepth') - well.stats.tvd) < 1e-9, 'net: well survey table, well depth and MD (no pipeline profile from a well)');
  ok(get(geometryLinks(suite('econ', all), well, {}), 'wellDepth') > 1000, 'econ: well depth');
  const nw = await importGeometry(F('field.csv', 'from,to,length,diameter,type\nWH1,MAN,2600,0.2,flowline\nWH2,MAN,900,0.2,flowline\nMAN,FPSO,7500,0.3,riser\n')), no = geometryLinks(suite('ops', all), nw, {});
  ok(get(no, 'network').length === 3 && Object.keys(get(no, 'network')[0]).join() === 'from,to,type,length,diameter' && get(no, 'network')[2].type === 'riser' && get(no, 'id') === 200 && Math.abs(get(no, 'length') - 10100) < 1e-9, 'ops: network table, median bore in mm, longest route length');
  ok(get(geometryLinks(suite('net', all), nw, {}), 'network').length === 3, 'net: network table');
  const dem = await importGeometry(F('bathy.asc', `ncols 4\nnrows 3\nxllcorner 0\nyllcorner 0\ncellsize 100\n${[-300, -310, -320, -330, -305, -315, -325, -335, -310, -320, -330, -340].join(' ')}\n`)), to = geometryLinks(suite('net', all), dem, {}), ter = get(to, 'terrain');
  ok(ter && ter.x.length === 4 && ter.y.length === 3 && ter.elev[0][3] === -340 && ter.geographic === false && get(to, 'waterDepth') === 340 && keysOf(to) === 'terrain,waterDepth', 'net: bathymetry grid offered as terrain with the deepest point as water depth');
  const wtg = await importGeometry(F('wt.csv', 'x (m),theta (deg),wt (mm)\n0,0,12.7\n0,180,9.8\n1,0,12.6\n1,180,10.4\n')), io = geometryLinks(suite('integ', all), wtg, {});
  ok(get(io, 'wtMap').kind === 'thickness' && get(io, 'wtMap').t[1][0] === 9.8 && get(io, 'minWt') === 9.8 && keysOf(io) === 'minWt,wtMap', 'integ: wall-thickness map and minimum wall thickness');
  const dpg = await importGeometry(F('dep.csv', 'x,theta,delta (m)\n0,0,0.001\n0,180,0.004\n50,0,0.002\n50,180,0.003\n')), so = geometryLinks(suite('solids', all), dpg, {});
  ok(get(so, 'depositMap').kind === 'deposit' && Math.abs(get(so, 'depositMax') - 4) < 1e-12 && keysOf(geometryLinks(suite('integ', all), dpg, {})) === '', 'solids: deposit map and maximum thickness converted to mm; integ ignores a deposit map');
  const cube = await importGeometry(F('cube.inp', `*NODE\n${V.map((v, i) => `${i + 1}, ${v.join(', ')}`).join('\n')}\n*ELEMENT, TYPE=C3D8\n1,1,2,3,4,5,6,7,8\n`)), mo = geometryLinks(suite('integ', all), cube, {});
  ok(get(mo, 'structMesh') === cube && keysOf(mo) === 'structMesh' && keysOf(geometryLinks(suite('flow', all), cube, {})) === '', 'integ: structural mesh passed through; flow gets nothing from a solid mesh');
  const span = { kind: 'polylines', name: 'span', format: 'test', polylines: [{ x: [0, 20, 40], y: [0, 0, 0], z: [-100, -100.8, -100], closed: false }], bbox: { min: [0, 0, -100.8], max: [40, 0, -100] }, stats: {} };
  ok(get(geometryLinks(suite('integ', all), span, {}), 'spanLength') === 40, 'integ: free-span length from a short, level polyline');
  const par = await importGeometry(F('params.csv', 'name,value\nReservoir pressure,250\nwt,12.7\n'));
  const po = geometryLinks(suite('pvt', [{ key: 'pRes', label: 'Reservoir pressure' }, 'wt']), par, {});
  ok(get(po, 'pRes') === 250 && get(po, 'wt') === 12.7, 'parameter table matched to inputs by key or label (fromTable)');
  const custom = geometryLinks(suite('flow', all, { geometry: () => [{ key: 'length', value: 42, from: 'custom hook' }, { key: 'notDeclared', value: 1, from: 'x' }] }), riser, {});
  ok(get(custom, 'length') === 42 && custom.filter((it) => it.key === 'length').length === 1 && !custom.some((it) => it.key === 'notDeclared'), 'suite.geometry() hook wins over the built-in rule and is filtered to declared inputs');
  ok(geometryLinks(suite('unknown', all), riser, {}).length === 0 && geometryLinks(suite('net', all), riser).length === 4, 'unknown suite id offers nothing; derived data are optional');
  const rp = routeProfile(riser, {}), pt = profileTable(rp, 20);
  ok(pt.length <= 20 && pt[0].z === -1000 && pt[pt.length - 1].z === 0 && routeProfile(well, {}) === null && routeProfile(dem, {}) === null, 'profileTable caps the row count and keeps both ends; wells and grids are not pipeline routes');
});

// ---- Fixtures of the binary container formats: gzip + base64 of files written by independent libraries, or of small real files
// (origin and licence beside each). The expected numbers in the checks come from those libraries, not from this reader.
const gunzip64 = (s) => new Uint8Array(zlib.gunzipSync(Buffer.from(s, 'base64')));
const FX = {
  // HDF5, superblock 0 / old-style groups / v1 B-tree chunks; written for this test with h5py 3.16 (libver earliest)
  h5old: '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',
  // HDF5, superblock 3 / link messages / extensible-array chunk index; written for this test with h5py 3.16 (libver latest)
  h5new: 'H4sIAAAAAAACA+1aC1QUVRi+82BddldYdAElH4uCiQ9akJZ8tDsLy7oIiYh1xMfippC4kDxWXK14SCZS4QOPouQDpDKJ8lWKdhJRyAcqlXXI6iBaZqUkPijKtDs7d2DY0YIepzzNt+ee79z//nPnznf//849d3a5UW/oKfGSEGIxYHEHocWPqauQfXjoL/Ioo34yDmx4L64jCwkOMADkcgAsA5lLpqJLRa4ANpHJ5pT0KryTyWp+Kv1tEWNy7wVIfq+gY2RgW0yzcxeGIFLQfTunzV9gjY+zLVqsxtD1FPjLsJyaVWUfAgYwzAcAHMewAaiNQOxYJxQAc1Nq0ZNQIEgKgu7ARic4WkIq9QAYjtMDx+wemP0HlQDiXnYfAhMTYkISgpQDL4D/CEprr1406Ix6gCkkJLKZ0Cj7+9WcNuj0IQBrHzeayNWcHlYj22aOLcCZ4Q0c2wlkK+DYHrVGP4GCgY457wcYcwKabBQKTJgtNCclhaLpEMnsGtvWkZ3qi7Y6ccNwWT+mdqAHup3bXQRoabk0xDEYZKjN2SGdnLnB4D0QMHdXAk8x8MwCbDCIPQAg6Sl3IklSSZJySRT5e9khC6bDTAQUdOOCpxOt6aAvhgGoh71gSHh6CugSK2E6s9sj/tZgqI9qHRJGBwMQK0k5KXFsZ3NwrYM+CsSH0ZrxqkRT0K4pbISSYs4O1/xuQknotCGwWC/G18gOYAe476AJGJvoqAUbX7JuaFGCtJh1H2uRX5T5hqMW7MqCc7S4W/jbtcAZLTC0TijvYy0O9zuoC5kCc82e2fTKNGeopT8bF2xkhGRoJhpiYsJAX9SE+ePgI7sm/cqw+TTbTDOPj7QsMYOlTdKXpfmDPGID5436XjqxfpDYQ7s9TH/oyecyj5ce25ekylNNrd+RpM6Lrxuk+fDwRn+XyLbzy/CsN4sWNVev97W0Xt5neGiX+dDzxMKLH07f49dyYuzmOaWbL9VM3zx/WFnpwCXDTp47mpXSW52ja1Nve3XZzb1N4zPOHHt5f1vqQ9qE4qB3w6+dJc4f8JtmM408GRCBm0F201MT5oa/QqxrlR/6JJf6QOljrLxulWTlXNgrPjstJosiX59b4TLDarXMdav+LFafH2qr9TSlb1r+2fXa4JgEReVXuXGbGi7qCp0GP70tpebBFWv2bT2+wbf3PJ/3rh9fqSa8bv1wKjZ15JrcjSoXk36vtyJeG12uaIw74rq0LD1ydXLwl9Znlu3JWmuIzph8+rmvd112LTcElNlMyXUVE5dE08MrHDJ7dtwFPElRuCfGZzdWTAQ/Iku9VTSi9u0q2XfuecUnZ3mfaEj8glw/XOnVcDRp/uLq8sPuPT/18RjXsrDaYmrb8oH0xxn1c449ODNz/fVTBerTZyrLaj0bp/zUVECsci/++NkTyQNa/bIKb+fOG7Dd7+Cqaeap32nX/LJ1hCYh9ZWa8NETX8yNNcsnj/nc7RmPAIvN9HBdiQTX9TnSohW5FdZkus4aV5wS93HKzFURB8L1M90jX6g4+6y0sbrnhrdMgW24yL3ufPMV1bfBqj5Fob3fSnk+f4zvS965Y79SLTdVPdCqyl0bF7ppZbXY5WpOzA0na5Ns+p3GFzXvBKZ89E1OQ8ljGf5Xr5U/ERxvOht0Y39FcP/+M/ybG+uWn+l5mxhxylczPjQq0h5p8s6rP5upyckOmYs4YlJA561Tu53ZhrLbDPatnJBosy5Ii2/PjJWuQIAAAQIECBAgQIAAAQIECBDw74BiSI64D2IlYh/EQxGPQKxC3M/Bb5hD+yjEasSjEY9DrEWsYzg7BNX1FH2+wB6KE6iwdfpsIirUGHH3T01Ajwnz+jfAUpyZbZhsnGQ/wMbZ8yKu4v7okIg9xmw/BCI79yTr4h1JzpXMHA6F3OLL2qQitZchhv6Swd4ZA525o6dJwAavZlua0LEUOmJu59v5tx4PmRJpgM/oQdd96YeDaRDxWgjIwzoMLrXW8cwhPTK8X3JzPKjgeIxa7DQFjOV4uK+IfRzc4BgO7K9IA0s5hpE7H8sBGk4fd/QXloDZHEPowmsFoC/HMKa5tgps5/Txc27lRXCaYxhUnHazblf5BCF+uwp6JRG50ZEiSl+QlhG/aDcbttxPeNZEa1L8PT/h9eCsTnQIX/F883KYLjwE9o4i2nsww6mIDyF282H4irJzbv1fsH31RmMY/Tm6XSkaZqTKTsQYWgOiEZf4/rNquf5KlQqZ8V+HdgvcwMC8oVKYjQw1CZaVsOyDJRuWbbDUM7sMqgoADbRpIVMqyHQ7TGMKrupUESx0X0dhuQrbzkE+1/X+qZ+Y/qlLHf1TlMNGa4tD/dxdN2Id9WyHOv0MebDsgKWhY7xawIy3u3rwnv8e/XdZj6rO49UB/vi7pccfzFe3N7p/MF9/cuN8z/nqLgxGuAay/22x4/XdlYEEfJvI6DeKWaUCcjFo/zMBznHkOAXwnDC+UyDPCec7jeI5EXynIJ4TyXd6mOfkxHdS85xEfKdgnlMPvtMjPCcx32k0z8mZ5xTAV1zCd+IrLhUWZgHdxm9U2vp68SkAAA==',
  // CGNS/HDF5 structured zone: VTK testing data (Testing/Data/bc_struct.cgns, BSD-3-Clause, Kitware)
  cgnsStruct: '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',
  // CGNS/HDF5 unstructured HEXA_8 zone: VTK testing data Example_grids.cgns (BSD-3-Clause, Kitware), converted from ADF to HDF5 with cgnsconvert 4.5
  cgnsHex: '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',
  // Salome MED 4 mesh of 6 tetrahedra; written for this test with meshio 5.3 (h5py)
  medTet: 'H4sIAAAAAAACA+1bMU/bQBQ+B0KiAFVQaOtSVfJSESFQA4qqCFE1bnIhqRLbsoMELKFqkcJQFbVIXTvSf9CxY0b+QRkzdmRkZOzYjTrxnZN3jmMnBKG27zGczveez7777vvuHc7XcrE0n1hKkI7F42SaJEm/XTNTXsA6bz9gpcTKU1a2Ivy61G2T2fUFdn/Rr25S2vG+Foz3cxl1yjhB+x+tTFWjU+6yeoaV7Qj0q2gl3WpsU42aapVahGpWo6aqpKRa3fZCWQV4TQf0K5EowG+S4Tcdc+rnpB+/xBe/aIjffvxyHjuNinhzzJhQv5amFyWb1TnODzYFvmb11pZTcpxfxaBfkj1wxlEKorDr+WfQr8nqp+vwPZUp6Fdm9Wb075g/rlvic88xJuqM8Zz9tjX1tTtezrvnwXuzcG9cRRsaR3ziTFoNFeenw+K8Bulwcx51GHmsx2N8t2gkoN/7w09NsN5v2K9EZgfq8Nk9puMh95HfUohfxG8PvzmOp3not5bx2Pqga6JuK0N1WHJx10rAuJ+s/msO5yiMDovjJepisVIbS4epZYwVZ1IjpH6n3bi4/bej1bv1lOT4BvtXAvzhc1lmfehzSW6cAfrR9JrbT6bPb1eBtifU94X6HFth/L5Farn3lfvuW7P1Qnn78fDNyeE75fPRSVPpKMjRB/F96nvhxlkiMwP14nIRvniQXnyXUS9QL3p6wfGUT0E/TacOjtUKwGP+RnlbTy8uknC92kTo6E0K5yiMXojjJfJKYTssT8I4rTg8jvOJJ043Q8b1+HPR9jRYfwvKSzcuO0uy10J+zdTxh9/51dUjdn5GwvHg+WPkQeRBb95XlqFfQdfJJPM+hwcjbn/xpxDnGVY3liE+L+9Dv/gDlifK/zcPiuMwKR4sbFsh4+B+0ChV3f3gw/79IC02NL1hmHqpUqWNilanpqZW7faY4yPx+zovIkkS0GfZ3a/LXVc/vozaVzvtkUhEcvDjRE5J8SeM+Ybur8fWjVfD+b/Hs/JI+4SzJeSqCVg+YPsW1M5nMQ/LUdtFS+b9ShFvQWVEyGenBdzxY/mZEf34fWMjlmhoaGhoaHdhkq1wg/LkNssrjJB58sUK5smYJ3vz5ONl6Fen2YnmyTy/E/PhcfMVv3MjZXW0c6Pfa7gecD1418OXFehnp7G3cG7Uy4vlNMR3ntUPVnCOwvCKOF7//rlRsusqPpf/eVH7Ts+LZtwelbHWx9UqYn28fePg7+rE7yAD/7+SRZ1EnfTqZGsd+k3+u7rB+E1vsPtnw+G3nUP8In6930k0s9CvpNYqVVuY96mpg/3YhPZ5ZAPqqIhjtNHOP5Kbg84/8HcuaOF4gK/L45yIN4irm5rIA/HnsJ8Mqxs5nKNQ50jCeHnykZ3h39cG7S/E3x0hr6CNwyvtrdvllT93fmZnID0AAA==',
  // Exodus II in NetCDF-4 / HDF5; written for this test with meshio 5.3 (netCDF4 1.7)
  exoTet: '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',
  // Exodus II in classic NetCDF (64-bit offset): VTK testing data box-noglom.ex2 written by vtkExodusIIWriter (BSD-3-Clause, Kitware)
  exoClassic: 'H4sIAAAAAAACA+1XO1PbQBCWsTEGhAlGgAMplBkKihQIh4RhmAkM4BkmTR6eyaNRZHzGN8iSR5bBpkqRIkV+QIoUyaRJk0mTJh1l6lT5ASn4ESlyK+/JF/mBgBThsTPnW+1+393u3urh9Y1snyRJETaG2BiB2SSWXnUdau1IIDfZSIDNpBZh+kM2YkW75sDMxqBLy4QRSEXyxcNbRpkgf8CqlfUCLYMvChy4tuwCqTbxM8ABGzEJgOYgjua1njd3PQzEmEQbtcCseXbAjuJ6eoU4zO85IDaZr6nvGRBvnI1hjm2avHhk7jMqVN8jTpXaFvj6Yd/Vj5kM5NDFni6atuGyaukVm1quvm87Bb1KDwjGAHEnIOciNQl3CD4Y/S51TbTDeSytO8RwSUHNN9Q9d3ezbhdq1a2tJw51iXNLzZVq6lrFUTVN1e4uLy4uz8+rC/NaBo5QGi8bdVqGHFn9dXYQO25J3E+FPNmYgBp7Z7dfsnHzCNaDSxxjh9Z4ATmQPDtpw63huUVafj+nGFJfgY/hK45d0QSsjHoM+wPyjW5ttPFfw8+2zapZ9/eKB2L7ieobH9s4DjtQ8LEHx2ETK5gDhFrFWBOtGP3zmm2qgwrUFNb2GU1/fwdOATlLsOa2bVlk29XQN4j3I6/VoNfDbqMi9Eg8t5l7tObXve8lrvcODss7e6HxkTPcIY7vTVW+wX7G9gyz6vM0ktdauCjGJdQpsoq9UWrjLoTlJj+1cTNhuaNHbdzbYbljc23cxbDcVHu+d8JyxyFfmdN018gLfTUcuJfSeDa/oMjefvjgCuzVpX+VoyBvIRRv4jDIy4TiTb71nrnQe60HbB/2cZf+l38L7x9PVhWcZ6W/5Lzbg+LjlM684Hxe1jmpdN33pPPsv5nDyNPjAM+OAzwPs00E75shvGdGUAdbEt8rfajH8PmRRBvncD7Hc38c8ZwvB76FwD6KuhzgcfwAG9dwjwTGF0ObuB9/Z46hngzwRDz3xxHP+Sl8psZRT2B8KSGmuMBPYTyTQv7jqINNwfwV9IF9CnUlwBPx3C8jnvPTGI+MehK/sdJoU4Q6Kmgfwdwn0X4ddbBNC/lPo20K9aTA4XyOnxbynxL4PL40+pL43Z1Cm8jz8FnqVF11nZimuuY4RqOtOx8T9tlS6IHIlajTyy9l2d8It9QDkaXFnn7pMa33Bvj3Udto5vcAPto7czG/7ohmfj1WuPBy7+vn+53ms/JPu965q9+HbzOHHeaz8k+73lX9Llf9ruRKLqOs/Ih+eX8B7Vfyf8sfVrdYe3AXAAA=',
  // NetCDF-4 bathymetry grid (packed int16, deflate + shuffle, _FillValue); written for this test with netCDF4 1.7
  nc4: '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',
  // MATLAB v7.3 MAT-file with md / inc / azi vectors; written for this test with hdf5storage 0.2
  mat73: 'H4sIAAAAAAACA+1W32vTUBS+SddZ11Y2RBziw3kRRGbIYrKVgphurUywLrgO9Mldm3QtpElJU1h8UQTBgSCK4usexb/CR/8EH/1PapJ77kxSUHxpH9rz8uV859cN3z0hzVrrYW0HtqU70Ky1bnd6trUBhk39juv1q9A1O9rQdz16bIEsKZKyAbueRX3LBNepQqs7gv22D3IFNrerqlrVNFBkZQv26vc1GLa7Vp/CpiTLIMEfI9zEB02ysJnaaahUeeXaSuwUCmSJrKI2DMZoZ5fSPq8/QhQQ3yJ+FTkvxLF15Newfzav9bjRiLLHGeNz5CWGhYVkc2l7jZoR4RP0+Ufkh5jO65sMe047Rvqil7qflf+cK5AiqxW4f5OhKKT6rmb2gPu58OZH3Brcw7urE7VI1HGOxfMhG8VFkTUs4KbkBCiyjFep85SIgXPK4fM6MQK/6zrSQZcOLD6Xz4nzBCG5zuc4eV7GnOQzgczGHTzarwvxW+H+/2Mh+ZyT5b/n8TaQm4/7XOL6kOvhuxe4js6oPwikQ8e0PDvoOcetgMlKLod6XEjI0rFd6m+pYZ8K9rk62WfXdXzacyxvst4xqefRIKyXz+9TVN+M/weetW06HPK5y4k60x09ty2SmFtKzk2dt5ycFx+IT434i3Hs59T36myqe7XQeVY6Q2mh8/R0ng97quPfBcPvHvoBQ/0l8q8Z/nqD8VOG8A7z3jOUP2D+R4bGJ6z7zHDwReeTD2PL6wyvIA/o39C/xXYLeRl9BeMa8hXkq8jfRV5Hv4bxHeTreJ6jmeFvquzI8aAOAAA=',
  // GeoPackage with a 3-D line layer and a point layer; written for this test with GDAL (geopandas / pyogrio)
  gpkgRoute: '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',
  // GeoPackage polygon layer: GDAL autotest data ogr/data/gpkg/poly_non_conformant.gpkg (MIT, GDAL project)
  gpkgPoly: 'H4sIAAAAAAACA+1ZTW/jRBgeO92mS1UpB0IQqGhkCW0jnCbOJz4g6rYmtE3sbJqIraLiDskktUjsYDtLs0LahK/bXpAQFw7LH0ArsUIrcVmO7IUrF6Q9UWmFhDjBkbHH2SYk0h7Rsn409szzzjMz7zuemVjO4dWS7mDYNq0ecmAGLAGGAVsQAgBWyRUGF7hErqUpzkwKxYqQAouxChJ/3ltjvwPLzCsgdI39gfmd7TNHhDwVEJeXY9Eo846D3uvivtkduld4pypLNRnWpO2SDF0L3GjrLR52sNmDFbV0VFSV+DhyKRwrFJjxm17jTv/9juYKsGMNtabZHfQMe6Fxeab7hRK44fWpGaiHeWqj5cdCZ9j3q23L1ohzN/hefMyEwrFkkhkLFy6ZHYv0ajjYcOw5w9K8K9PVU17AmnytBhWVXPVSCVaqe2WpegQP5CO+jZEzsDBpNTCcuLS0HKsnGaAbLXxmf9Ali09DA8f0uDY3hCbMmdxluOI+nPFLbDh25Qozti/CmQllQkLzYUyqZmayhRzkzRyvt0il3taxxbew3bT0vqObBt9FtqM1T5HRwXxPN7Qz7z7ke+jMLZP70J/v+PhVJhzL5Zhx78I5u48cHXU1C7c1e2gvsrHzrv5LATfcEaYfrWl1kKHfQJ6L04REaVot0sZVtXBbN3RPMxVSfM3dyGwSkBQgwJPwicKsqGyU/bwI1MPiNhTETB4m4bZFVpZ9ChVv3aEuLFp6K3KuqIocOa9U1f2dwwb3xAYcX5TV4oyU43elWr1MLRq1HFbelqvq3m6Dk3Ry1AmvZ1Icn88UCrl8ZjMj5vm0KG5m0lkxn80fH/PuWSSTHooWxsaHevOU41PHfF3ZqzW4XdwhVmLYTAmFbC6TFtM5QRSzpJDzmqr78k5tT1UaXM1Chn0dWzbWythqIse0OKKQqlJZrsnVBtcloTiDFtbMtmZaekc3OD4rzkia5FSxyFbuYRKtjoggkZ4R2E1EDqM2arq9E69EUcynhHRBKMzI2qhL3MDkNNKNDhkl5WKBwjAt59STJARfQ+MuYweTAYTjY2//h34EJAUIEOBpwxvManQ9urYWJmDdl1H/fct2y9K6cMcVSdXXtrx8Hd+lPPFoRPLIOd3/fwGSAgQI8GwgzISil90T4jlv//8NSAoQIMCzcgD4+5/+/v8GSAoQIMD/Ei8y4Wj0ZXZlxd3y7kdi/yu1//7PfgpIChDgP8HHQ2bp8mfXixUQOo9Q00dvrUsku/9C1cvhZqLq5rcfkZwJAfrn0/NUe5Lco6KbVHR/o+jxh19QPuns5FvKRz/Revg95fDrEtX/QvnJLcpv/+rX31TWHw/u8rjq8dEDv/+7V2n7uj++7/RIpvxhTJkJAkplqi/59QeUn7xL+VbTr3f8oP/w/bvl6xWff+PzOwe0/T0/vq/26fiTeH+m8zP60o9ver7+AX61UvgAHAAA',
  // ACIS 1.6 SAT part: NIST Design, Process Planning and Assembly Repository, Allied-Signal part09.sat (U.S. Government work, public domain)
  satPart: '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',
  // Parasolid 26 text file of a spar-buoy hull (two cylinders joined by a cone): github.com/yang7857854/OpenF2A CADFiles/HywindSpar.x_t (Apache-2.0)
  xtSpar: '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',
  // MicroStation V7 2-D design file: GDAL autotest data ogr/data/dgn/smalltest.dgn (MIT, GDAL project)
  dgnSmall: 'H4sIAAAAAAACA+1ZeUhUQRj/zXu77ZHVhtllxRqY0EkJXSiuq9l2Z0pFUZZmEaUFFYQd+7Iw/wmsKKI/IiqCKKjEwCJi+6MoCiqiA4kO6L6QoEAtp5k3b3Pd3F2X1kya3+Pbub75rveYb2bWamtWiBVtwomyKf76W9UOmL9X019Y/JXS5gY+RnD/FkHOCsAEzdW2rFD9oTAu/aKaU+EYzustWpvpu2ZKx/bivUmoGT4ARdWx1Hp48jDXqY+Hknl9a4CviS8oPf4BGoECu5ngYBXgMMVK65G0wbmLjt5o5JEM1Kq+ofTOE6H13RMFdRNiqbUgbWnBl+2+p8Faq59TevkW5zCj8a4JJ3fEMsKL0xyvq67cSAmOcLCvTTH2dXTF7sL62kgRbuwUrR3na8N5YghhD3RioF5mhFfYEklrJO1tjXdD7R7Arq8aJZs2bAolM1TJopEiqDKeWlTmhyVAemUlHjOvJul8Lb4w+uWjFjlKhCwoMWptjTLqXubKnbjQ1eNIotLu2Gva9TLEEGasNJ6OgIn/WM+Fi1Lr+LQjvoGxWvZefyudhH/HEthn4l8FYfk/3NuPvA5Fuz51hAQJCQkJiaizpHl9qzYNAFCo93nm5hdkZxdMm5qdNzlnzuz8mdPco51R6PB2KrZ7Le6RJLSPQ4i+ETJFyJNGNkxHBtgGXs+Z1ETVcDkzA0vwUn2lBPYtRCZ+qDcV+eV1TQSfFmzu+Urob8uj8O+qJitYimqQwCzX53mjdsGk+TraeqlX6pV6/xj94Yzm3BLIZ2NP67Npe2FlD9Hn+uWp+r0I/Q0wxhOSoKQgGQm2uPkXqLhDIe1Yr8wtFgeV/v5lmX6b2Bz4UMUoi5GHhZ63Pf663oYgFjOPg5NPkMEXFq7oykhvZAazYbpDZjEJia6Bvu7rYc4ntca+3b/DCtxpWdhqWf9Y3Po5Uv623U6Dfs9XXRsD3ZfC5I8zev5AvPxuJSQkJGKFfu5vIdddgreKSC0JaGLtLc6s5YVri51jnPnFrNgmoychISHRdUB6JwHatUHQZu/V29rOddBKn4l9drx+haRaBmBNKScxh3Hj4QkbsotL1rGlf/NGKL2GMikP4qg347uQUlpGtEerhZQ+upTxHsfte3VHz/pLzpflwTF+u5PqQXfeVrqNYnImnoa2X/wlpCVdhZZaLuT0A3oiz8y6R5wH9u0AGCcOrAIYF/gM32Ugtdw/bnFkQrES5Bm+Kmp/JvFwIbTyZ0J6XRG05E9Cup1byVhZH/g44wMfi7W0/xGU/gQzVdHgACoAAA==',
  // MicroStation V7 3-D design file with a global origin; written for this test with the GDAL DGN driver
  dgn3d: 'H4sIAAAAAAACA8WVD0yUdRjHn4d/onje4fHvOOV37UCmLnfkllhJRyOzrfRytVQY3rD4k2NnKgmyuvcoElMTN1eo5UEyaMwaCmm5srM2jHLIFiP5k1wYcRC+HAaEy+7X83LAAeGi7abP7fb+3t/7/X2+v/d5nvd9v5vr8glEmDE0IOjHx4rYeQD+7bV8LAq1uduflWH/Nf5U+RawHV8L9SvWw87znhUzx2yvd64unjuiCnw8Fj5P9rju1j6dkSRbaVjPFbkfwp8P50FmUxamn/CWa8HqoiTNKdvyM+DclzLhWqUduJIu23eY8/kugMpogIxFAA5fb7keSLBX7zqXmQIwiB5XFzfmcZ7scLvuMQAkJ3jT1bUmvTB1xYh2qivnX/cCnEhyu/J8gCCjN11jEqzLLupKnpzu+mkzwAaV27XteYBLid50HVhzZ3NDzUjuVFeAY0szuJDldm0l12+87NqUEnc98NpU1xIQ3gPYFnsvXG+fOaXIfAkgQABBklBFyZWbaRtm6f695zo5EmHVfokNgJAN2dk1SaPuE9r/OtLT7+f+O3w5RsCrU15QRQVQmxDnJ035O+E1beTGTbLLJs4DXPJoU81+WV9+AJgPXWzQfCmpzzamNkrH3d9DtHT0LbQ3ABxVIr6Y7WFWFvz7LoL26p+L36SXnVT70JlcDrMIQajbO/0M76KazTxC+tjPnQ/vxtuUDgysvut1X9DZZ9rn6LfhCfhf8/c+fCZGxl6p77n5vm1l3jMzTqNw/7MkzI+CFovNX2bVN4lmjWnY4qxj1rSrohl6AKjlAgGWoieXcxRpqZZJz4LgNzRpvTFnyGLQR1ljzrrXd49qlkALkoKRgkNPCyrqFzCeJnJ7VysaDy5g/aLIbZ1tWL1Bzg7liLy0ox1BKWfioMiF9nY0NCvY4HaRG1t+xtJjCpbVKnJ983V0pgQzWaLINT92oF4TzPLKiN/YgcU3gtnw38S/Ykd7xUL22Gbi1/+CulcWsiXniF/XicKDSrY1nPjfdmKjU8nUZuLbbqCmVsmaaD/6r37FHXtCWEwG8S90oe2REPZZN/HPd6HCFcqSTMSv+Q2Nl0LZD38Rv7obq98KYxfeIf4nDoR14exRFfGretAQFM7eLCd+ZQ+WNkSwmGXEL+9F55EItvFj4pf9jvoXVCwoivgn+7A4QsUSDxD/eB/a21TsJyfx37+JOmskKzAQ/6iIQlok++g08Y/0Y2OMmt2eQ/x3+1HTo2ZbXyb+ASfuOK1muy4Tv2gAbaZFzKAjfuEtVKykMeXKbrmFxuHFrDeS+G/8gdVfLGZlHxA/fxDh9SgWEkP83CGkijKqKDfmCH5a8FR4vFu2XZUqDxA69taYvGIYS30Ze4iqpN85jKRmpOYaE/gtB4cF1oVZ4Y5odphuWnRWrRXiRLNQNUoKlEhSHsMYKSiPUi7UzH5Y5H0lUi60jNTcYZrekRCwClotxoNyq32L02yoaLeAMtjaWO40xxP5gVGyP0idJ2ek4PFV7m7zjIMZqbmhQtJMHo9rpjty/g9u8VIERgsAAA==',
  // MicroStation V7 3-D design file with complex shapes of lines and arcs: NIST Design Repository, Bentley/MODELER/HOUSING.DGN (U.S. Government repository, public domain)
  dgnHousing: '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',
  // AutoCAD 2000 drawing (AC1015); drawn for this test with ezdxf 1.4 and written as DWG with LibreDWG 0.13.3 dxf2dwg
  dwg2000: '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',
  // ASTM E57 scan of 60 points with a pose; written for this test with pye57 0.4 (libE57Format)
  e57: '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',
  // LAZ, LAS 1.2 point format 1 (150 points); written for this test with laspy 2.7 + lazrs
  laz1: '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',
  // LAZ, LAS 1.2 point format 3 with extra bytes (200 points); written for this test with laspy 2.7 + lazrs
  laz3: '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',
  // LAZ, LAS 1.4 point format 6 (layered compression, rejected); written for this test with laspy 2.7 + lazrs
  laz6: 'H4sIAAAAAAACA/NxDHZjwAIYWfxDPFyDGHCDnMTigkoFIz1zPQOs8pKMr9jLGa8yAs0C8trkcBpULbLO/WFViz06zYAPNHg4kMRHA0YMAw2AwVeVWaCQmpecn5KawsCwJ0yDIaOkpMBKXx8ipZdflI6mhxmImZgg7IDDDAz/0QAjAxeDHFDVXiZiXAAKg0ggdiSksKJR3Sws7e9UwSq9S3uo5H1QmrDxN02S1jmUenPOt2mfN5oFGe76Yn/84+1YrwUPy679FJzNqck8xUNMcYusxoM7DluO2DhHpx96+Oi5K781B6uJJpJBLheABADbMXLdygIAAA==',
};
// Hand-made ACIS 7.0 text file of a pipe spool (two straights and a 90° bend, bore 0.16 m, outside 0.2 m), laid out like
// the entity records of real SAT exports; no CAD-written pipe spool with a redistribution licence was found.
const SPOOL_SAT = '700 0 1 0\n@9 hand-made @8 ACIS 7.0 @24 Thu Jan 01 00:00:00 2026\n1000 9.9999999999999995e-007 1e-010\nbody $-1 -1 $-1 $1 $-1 $-1 #\nlump $-1 -1 $-1 $-1 $2 $0 #\nshell $-1 -1 $-1 $-1 $-1 $35 $-1 $1 #\npoint $-1 -1 $-1 0 0 0.08 #\nvertex $-1 -1 $-1 $6 $3 #\nellipse-curve $-1 -1 $-1 0 0 0 1 0 0 0 0 0.08 1 I I #\nedge $-1 -1 $-1 $4 0 $4 6.2831853071795862 $40 $5 forward @7 unknown #\npoint $-1 -1 $-1 0 0 0.1 #\nvertex $-1 -1 $-1 $10 $7 #\nellipse-curve $-1 -1 $-1 0 0 0 1 0 0 0 0 0.1 1 I I #\nedge $-1 -1 $-1 $8 0 $8 6.2831853071795862 $38 $9 forward @7 unknown #\npoint $-1 -1 $-1 1 0 0.08 #\nvertex $-1 -1 $-1 $14 $11 #\nellipse-curve $-1 -1 $-1 1 0 0 1 0 0 0 0 0.08 1 I I #\nedge $-1 -1 $-1 $12 0 $12 6.2831853071795862 $52 $13 forward @7 unknown #\npoint $-1 -1 $-1 1 0 0.1 #\nvertex $-1 -1 $-1 $18 $15 #\nellipse-curve $-1 -1 $-1 1 0 0 1 0 0 0 0 0.1 1 I I #\nedge $-1 -1 $-1 $16 0 $16 6.2831853071795862 $46 $17 forward @7 unknown #\npoint $-1 -1 $-1 1.3 0.3 0.08 #\nvertex $-1 -1 $-1 $22 $19 #\nellipse-curve $-1 -1 $-1 1.3 0.3 0 0 1 0 0 0 0.08 1 I I #\nedge $-1 -1 $-1 $20 0 $20 6.2831853071795862 $62 $21 forward @7 unknown #\npoint $-1 -1 $-1 1.3 0.3 0.1 #\nvertex $-1 -1 $-1 $26 $23 #\nellipse-curve $-1 -1 $-1 1.3 0.3 0 0 1 0 0 0 0.1 1 I I #\nedge $-1 -1 $-1 $24 0 $24 6.2831853071795862 $56 $25 forward @7 unknown #\npoint $-1 -1 $-1 1.3 1.3 0.08 #\nvertex $-1 -1 $-1 $30 $27 #\nellipse-curve $-1 -1 $-1 1.3 1.3 0 0 1 0 0 0 0.08 1 I I #\nedge $-1 -1 $-1 $28 0 $28 6.2831853071795862 $64 $29 forward @7 unknown #\npoint $-1 -1 $-1 1.3 1.3 0.1 #\nvertex $-1 -1 $-1 $34 $31 #\nellipse-curve $-1 -1 $-1 1.3 1.3 0 0 1 0 0 0 0.1 1 I I #\nedge $-1 -1 $-1 $32 0 $32 6.2831853071795862 $58 $33 forward @7 unknown #\nface $-1 -1 $-1 $41 $37 $2 $-1 $36 forward single #\nplane-surface $-1 -1 $-1 0 0 0 -1 0 0 0 0 1 forward_v I I I I #\nloop $-1 -1 $-1 $39 $38 $35 #\ncoedge $-1 -1 $-1 $38 $38 $-1 $10 forward $37 $-1 #\nloop $-1 -1 $-1 $-1 $40 $35 #\ncoedge $-1 -1 $-1 $40 $40 $-1 $6 forward $39 $-1 #\nface $-1 -1 $-1 $47 $43 $2 $-1 $42 forward single #\ncone-surface $-1 -1 $-1 0 0 0 1 0 0 0 0 0.1 1 I I 0 1 0.1 forward I I I I #\nloop $-1 -1 $-1 $45 $44 $41 #\ncoedge $-1 -1 $-1 $44 $44 $-1 $10 forward $43 $-1 #\nloop $-1 -1 $-1 $-1 $46 $41 #\ncoedge $-1 -1 $-1 $46 $46 $-1 $18 forward $45 $-1 #\nface $-1 -1 $-1 $53 $49 $2 $-1 $48 reversed single #\ncone-surface $-1 -1 $-1 0 0 0 1 0 0 0 0 0.08 1 I I 0 1 0.08 forward I I I I #\nloop $-1 -1 $-1 $51 $50 $47 #\ncoedge $-1 -1 $-1 $50 $50 $-1 $6 forward $49 $-1 #\nloop $-1 -1 $-1 $-1 $52 $47 #\ncoedge $-1 -1 $-1 $52 $52 $-1 $14 forward $51 $-1 #\nface $-1 -1 $-1 $59 $55 $2 $-1 $54 forward single #\ncone-surface $-1 -1 $-1 1.3 0.3 0 0 1 0 0 0 0.1 1 I I 0 1 0.1 forward I I I I #\nloop $-1 -1 $-1 $57 $56 $53 #\ncoedge $-1 -1 $-1 $56 $56 $-1 $26 forward $55 $-1 #\nloop $-1 -1 $-1 $-1 $58 $53 #\ncoedge $-1 -1 $-1 $58 $58 $-1 $34 forward $57 $-1 #\nface $-1 -1 $-1 $65 $61 $2 $-1 $60 reversed single #\ncone-surface $-1 -1 $-1 1.3 0.3 0 0 1 0 0 0 0.08 1 I I 0 1 0.08 forward I I I I #\nloop $-1 -1 $-1 $63 $62 $59 #\ncoedge $-1 -1 $-1 $62 $62 $-1 $22 forward $61 $-1 #\nloop $-1 -1 $-1 $-1 $64 $59 #\ncoedge $-1 -1 $-1 $64 $64 $-1 $30 forward $63 $-1 #\nface $-1 -1 $-1 $71 $67 $2 $-1 $66 forward single #\ntorus-surface $-1 -1 $-1 1 0.3 0 0 0 1 0.3 0.1 1 0 0 forward_v I I I I #\nloop $-1 -1 $-1 $69 $68 $65 #\ncoedge $-1 -1 $-1 $68 $68 $-1 $18 forward $67 $-1 #\nloop $-1 -1 $-1 $-1 $70 $65 #\ncoedge $-1 -1 $-1 $70 $70 $-1 $26 forward $69 $-1 #\nface $-1 -1 $-1 $77 $73 $2 $-1 $72 reversed single #\ntorus-surface $-1 -1 $-1 1 0.3 0 0 0 1 0.3 0.08 1 0 0 forward_v I I I I #\nloop $-1 -1 $-1 $75 $74 $71 #\ncoedge $-1 -1 $-1 $74 $74 $-1 $14 forward $73 $-1 #\nloop $-1 -1 $-1 $-1 $76 $71 #\ncoedge $-1 -1 $-1 $76 $76 $-1 $22 forward $75 $-1 #\nface $-1 -1 $-1 $-1 $79 $2 $-1 $78 forward single #\nplane-surface $-1 -1 $-1 1.3 1.3 0 0 1 0 0 0 1 forward_v I I I I #\nloop $-1 -1 $-1 $81 $80 $77 #\ncoedge $-1 -1 $-1 $80 $80 $-1 $34 forward $79 $-1 #\nloop $-1 -1 $-1 $-1 $82 $77 #\ncoedge $-1 -1 $-1 $82 $82 $-1 $30 forward $81 $-1 #\nEnd-of-ACIS-data\n';
const len3 = (g) => g.polylines.reduce((s, p) => { const n = p.x.length; let l = 0; for (let i = 1; i < n + (p.closed ? 1 : 0); i++) l += Math.hypot(p.x[i % n] - p.x[i - 1], p.y[i % n] - p.y[i - 1], p.z ? p.z[i % n] - p.z[i - 1] : 0); return s + l; }, 0);
const fx = (k) => gunzip64(FX[k]);

await group('HDF5 containers: plain HDF5, CGNS, MED, Exodus II, NetCDF-4, MATLAB v7.3', async () => {
  for (const [key, label] of [['h5old', 'HDF5 (superblock 0, v1 B-trees)'], ['h5new', 'HDF5 (superblock 3, link messages, extensible array)']]) {
    const g = await imp(key + '.h5', fx(key));
    ok(g.kind === 'polylines' && g.format === 'HDF5' && g.pathway === 'numeric' && g.polylines.length === 1 && g.polylines[0].x.length === 30 && g.polylines[0].z[29] === -158 && g.stats.dataset === '/survey/route_xyz', `${label}: N × 3 data set (gzip + shuffle chunks) read as a 3-D route`);
    near(len3(g), 296.91684766613804, 1e-9, `${label}: route length as computed by numpy`);
    ok(Array.isArray(g.contents) && g.contents.some((c) => c.path === '/maps/wall' && c.kind === 'dataset' && c.shape.join() === '12,9' && c.type === 'float32') && g.contents.some((c) => c.path === '/tags' && c.type === 'string') && g.contents.some((c) => c.path === '/maps' && c.kind === 'group'), `${label}: g.contents lists groups and data sets with shape and type`);
    const t = await importGeometry(F(key + '.hdf5', fx(key)), { dataset: 'maps/wall' });
    let sum = 0;
    for (const r of t.records) for (const h of t.headers) sum += r[h];
    ok(t.kind === 'table' && t.records.length === 12 && t.headers.length === 9 && sum === 426 && t.records[3].wall_5 === 4 && t.records[11].wall_9 === 4, `${label}: opts.dataset picks a big-endian float32 chunked matrix (sum ${sum})`);
    await throws(`${label}: unknown data set`, () => importGeometry(F('a.h5', fx(key)), { dataset: 'nope' }), /no numeric data set "nope".*route_xyz/);
  }
  await throws('HDF4 file', () => importGeometry(F('a.hdf', cat(le('u8', [0x0e, 0x03, 0x13, 0x01]), new Uint8Array(64)))), /HDF4.*h4toh5/);
  // CGNS: counts, surface area and volume as read by VTK's CGNS reader
  const cs = await imp('bc_struct.cgns', fx('cgnsStruct')), ds = dimensions(cs);
  ok(cs.kind === 'mesh' && cs.format === 'CGNS' && cs.pathway === 'mesh' && cs.stats.structured && cs.stats.nodes === 125 && cs.stats.cells === 64 && ds.closed && Math.abs(ds.area - 600) < 1e-9 && Math.abs(ds.volume - 1000) < 1e-9, `CGNS structured zone: 125 nodes, 64 cells, area ${ds.area}, volume ${ds.volume}`);
  const ch = await imp('grids.cgns', fx('cgnsHex')), dh = dimensions(ch);
  ok(ch.kind === 'mesh' && ch.stats.nodes === 8 && ch.stats.cellTypes.hex === 1 && dh.closed && Math.abs(dh.area - 24) < 1e-9 && Math.abs(dh.volume - 8) < 1e-9 && ch.bbox.min.join() === '-1,-1,-1', 'CGNS unstructured HEXA_8 section: one cell of area 24 and volume 8');
  await throws('CGNS in the ADF container', () => importGeometry(F('old.cgns', cat('@(#)ADF Database Version B02012>', new Uint8Array(200)))), /ADF.*cgnsconvert/);
  await throws('HDF5 that is not CGNS', () => importGeometry(F('x.cgns', fx('h5old'))), /not a CGNS tree/);
  // MED and Exodus II: the same 6 tetrahedra of a 2 × 3 × 4 box written by meshio
  for (const [nm, key, fmt] of [['tet.med', 'medTet', 'Salome MED'], ['tet.exo', 'exoTet', 'Exodus II']]) {
    const g = await imp(nm, fx(key)), d = dimensions(g);
    ok(g.kind === 'mesh' && g.format === fmt && g.stats.nodes === 8 && g.stats.cellTypes.tet === 6 && d.closed && Math.abs(d.volume - 24) < 1e-9 && Math.abs(d.area - 52) < 1e-9 && g.bbox.max.join() === '2,3,4', `${fmt} (HDF5): 8 nodes, 6 tetrahedra, volume ${d.volume}, area ${d.area}`);
  }
  const ec = await imp('box.ex2', fx('exoClassic')), dc = dimensions(ec);
  ok(ec.kind === 'mesh' && ec.format === 'Exodus II' && ec.stats.container === 'CDF-2' && ec.stats.nodes === 27 && ec.stats.elementBlocks.TETRA === 40 && dc.closed && Math.abs(dc.volume - 1000) < 1e-9 && Math.abs(dc.area - 600) < 1e-9, 'Exodus II (classic NetCDF, 64-bit offset): 27 nodes, 40 tetrahedra, volume 1000');
  await throws('Exodus II that is neither NetCDF flavour', () => importGeometry(F('a.exo', 'not a mesh at all, just text')), /NetCDF/);
  // NetCDF-4: sum and count of valid cells as read back by the netCDF4 library
  const nc = await imp('bathy.nc', fx('nc4'));
  let s4 = 0, n4 = 0;
  for (const r of nc.grid.z) for (const v of r) if (v === v) { s4 += v; n4++; }
  ok(nc.kind === 'grid' && nc.format === 'NetCDF-4 (HDF5)' && nc.pathway === 'gis' && nc.grid.x.length === 41 && nc.grid.y.length === 31 && nc.grid.geographic && nc.grid.x[0] === 2 && nc.grid.x[40] === 4 && nc.grid.y[0] === 60 && nc.grid.y[30] === 61.5 && n4 === 1266 && nc.grid.nodata === 5 && nc.stats.variable === 'elevation' && nc.stats.units === 'm', `NetCDF-4 grid: 41 × 31 cells on lon / lat, ${n4} valid, y re-ordered to ascend`);
  near(s4, -379862.9, 1e-6, 'NetCDF-4 grid: scale_factor / add_offset / _FillValue applied to deflated int16 chunks');
  ok((await importGeometry(F('bathy.nc4', fx('nc4')))).kind === 'grid' && formatOf('a.nc4').name === 'NetCDF-4 (HDF5)', '.nc4 extension is registered');
  // MATLAB v7.3
  const m7 = await imp('survey.mat', fx('mat73'));
  ok(m7.kind === 'polylines' && m7.format === 'MATLAB MAT-file' && m7.stats.wellSurvey && m7.stats.md === 1500 && m7.stats.stations === 16 && m7.stats.maxInclination === 40 && m7.stats.container === 'HDF5 (v7.3)', 'MATLAB v7.3: md / inc / azi vectors become a well survey like a Level 5 file');
});

await group('GeoPackage', async () => {
  const g = await imp('route.gpkg', fx('gpkgRoute'));
  ok(g.kind === 'polylines' && g.format === 'GeoPackage' && g.pathway === 'gis' && g.polylines.length === 2 && g.polylines[0].x.length === 40 && g.polylines[1].x.length === 3 && g.markers.length === 3 && g.geographic === false && g.stats.layers.length === 2 && g.stats.layers[0].srs === 'WGS 84 / UTM zone 31N', 'GeoPackage: two feature layers merged, 43 vertices and 3 point markers in UTM');
  near(len3(g), 3211.060270900347, 1e-6, 'GeoPackage: 3-D length of the LineString Z features as read by GDAL');
  ok(g.polylines[0].name === 'Export line A' && g.polylines[0].layer === 'pipelines' && g.polylines[0].attrs.od_mm === 323.9 && g.polylines[1].attrs.kp0 === 12 && g.polylines[1].z[2] === -151, 'GeoPackage: polylines carry layer, name and attribute values');
  const one = await importGeometry(F('route.gpkg', fx('gpkgRoute')), { layer: 'structures' });
  ok(one.kind === 'points' && one.count === 3 && one.points[2] === -120 && one.attributes.length === 3 && one.attributes[1].tag === 'MAN-2' && one.attributes[2].depth === 151, 'GeoPackage: opts.layer reads the point layer with its attributes');
  await throws('GeoPackage: unknown layer', () => importGeometry(F('route.gpkg', fx('gpkgRoute')), { layer: 'zz' }), /no feature layer "zz".*pipelines/);
  const p = await imp('poly.gpkg', fx('gpkgPoly'));
  ok(p.kind === 'polylines' && p.polylines.length === 1 && p.polylines[0].closed && p.polylines[0].x.length === 19 && p.stats.layers[0].type === 'POLYGON', 'GeoPackage (GDAL autotest): one polygon ring of 19 vertices');
  near(len3(p), 2700.247416, 1e-5, 'GeoPackage: ring perimeter as read by GDAL');
  await throws('SQLite that is not a GeoPackage header', () => importGeometry(F('a.gpkg', cat('SQLite format 3\0', new Uint8Array(200)))), /corrupt|GeoPackage/);
});

await group('Kernel B-rep files: ACIS SAT and Parasolid X_T with the pipe recogniser', async () => {
  // hand-made spool: centreline 1 + 0.3·π/2 + 1 m, bore 0.16 m, outside 0.2 m
  const sp = await imp('spool.sat', SPOOL_SAT), L = 2 + 0.15 * Math.PI;
  ok(sp.kind === 'polylines' && sp.format === 'ACIS SAT/SAB' && sp.pathway === 'cad' && sp.pipe.isPipe && sp.polylines.length === 1 && sp.polylines[0].z && sp.stats.units === 'm' && sp.stats.diameters.join() === '0.16,0.2', 'ACIS spool: recognised as a pipe, centreline polyline with z, diameters 0.16 and 0.2');
  near(sp.pipe.length, L, 1e-9, 'ACIS spool: centreline length 2 + 0.3·π/2');
  const bend = sp.stats.pipeRuns.find((r) => r.kind === 'bend'), str = sp.stats.pipeRuns.filter((r) => r.kind === 'straight');
  ok(str.length === 2 && str.every((r) => Math.abs(r.length - 1) < 1e-12) && bend && Math.abs(bend.bendRadius - 0.3) < 1e-12 && Math.abs(bend.angle - 90) < 1e-9 && sp.polylines[0].x[0] === 0 && Math.abs(sp.polylines[0].y[sp.polylines[0].y.length - 1] - 1.3) < 1e-12, 'ACIS spool: two 1 m straights and a 90° bend of radius 0.3, chained end to end');
  const sm = await importGeometry(F('spool.sat', SPOOL_SAT), { prefer: 'mesh' }), dm = dimensions(sm), vol = Math.PI * (0.1 ** 2 - 0.08 ** 2) * L;
  ok(sm.kind === 'mesh' && dm.closed && dm.consistentNormals && Math.abs(dm.volume / vol - 1) < 5e-3 && sm.stats.tessellatedFaces === 8 && sm.pipe.diameters.length === 2, `ACIS spool: opts.prefer = "mesh" gives a closed wall of volume ${dm.volume.toFixed(5)} (exact ${vol.toFixed(5)})`);
  // real part: area and volume of the same body as computed by OpenCASCADE from its STEP twin
  const pt = await imp('part09.sat', fx('satPart')), dp = dimensions(pt);
  ok(pt.kind === 'mesh' && pt.stats.version === '1.6' && pt.stats.faces === 31 && pt.stats.tessellatedFaces === 31 && pt.stats.surfaces.plane === 16 && pt.stats.surfaces.cylinder === 15 && dp.closed && dp.consistentNormals, 'ACIS 1.6 part (NIST): 16 planar and 15 cylindrical faces tessellated into a closed surface');
  ok(Math.abs(dp.volume / 0.1016208407 - 1) < 1e-3 && Math.abs(dp.area / 2.89343471 - 1) < 1e-3 && Math.abs(pt.bbox.min[0] + 0.334) < 1e-9 && Math.abs(pt.bbox.max[2] - 0.091) < 1e-9, `ACIS part: volume ${dp.volume.toFixed(6)} and area ${dp.area.toFixed(5)} within 0.1 % of the exact B-rep values`);
  await throws('binary ACIS (.sab)', () => importGeometry(F('a.sab', cat('ACIS BinaryFile', new Uint8Array(64)))), /binary ACIS.*\.sat/);
  await throws('SAT without entities', () => importGeometry(F('a.sat', '700 0 1 0\n@4 test @8 ACIS 7.0 @4 date\n1 1e-6 1e-10\nEnd-of-ACIS-data\n')), /no entity records/);
  // Parasolid: spar hull = cylinder Ø 9.4 × 108, cone 12 long, cylinder Ø 6.5 × 10 (dimensions of the published spar concept)
  const xt = await imp('spar.x_t', fx('xtSpar'));
  ok(xt.kind === 'polylines' && xt.format === 'Parasolid' && xt.pipe.isPipe && xt.stats.units === 'm' && xt.stats.version === '26.0' && xt.stats.diameters.join() === '6.5,9.4' && xt.stats.pipeRuns.length === 3 && xt.stats.pipeRuns.some((r) => r.kind === 'reducer' && Math.abs(r.length - 12) < 1e-9), 'Parasolid: two cylinders and a cone recognised as one run with diameters 6.5 and 9.4 and a reducer');
  ok(Math.abs(xt.pipe.length - 130) < 1e-9 && xt.polylines.length === 1 && Math.min(...xt.polylines[0].z) === -120 && Math.max(...xt.polylines[0].z) === 10, 'Parasolid: centreline from z = -120 to z = 10 (130 m)');
  const xm = await importGeometry(F('spar.x_t', fx('xtSpar')), { prefer: 'mesh' }), dx = dimensions(xm), vx = Math.PI * (4.7 ** 2 * 108 + 3.25 ** 2 * 10 + 4 * (4.7 ** 2 + 4.7 * 3.25 + 3.25 ** 2));
  ok(xm.kind === 'mesh' && dx.closed && dx.consistentNormals && Math.abs(dx.volume / vx - 1) < 5e-3 && xm.stats.surfaces.cone === 1 && xm.stats.surfaces.cylinder === 2, `Parasolid: tessellated hull is closed with volume ${dx.volume.toFixed(1)} (exact ${vx.toFixed(1)})`);
  const xhead = '**ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz**************************\n**PARASOLID !"#$%&\'()*+,-./:;<=>?@[\\]^_`{|}~0123456789**************************\n**PART1;MC=x;FORMAT=text;\n**PART2;SCH=SCH_900000_9008;USFLD_SIZE=0;\n**PART3;\n**END_OF_HEADER*****************************************************************\n';
  await throws('Parasolid file of modeller 9 (schema not embedded)', () => importGeometry(F('old.x_t', xhead + 'T50 : TRANSMIT FILE created by modeller version 90000016 SCH_900000_90080 12 1 0 0\n')), /schema 9008.*embedded/);
  await throws('binary Parasolid (.x_b)', () => importGeometry(F('a.x_b', cat(xhead.replace('FORMAT=text', 'FORMAT=binary'), 'PS\0\0\0', new Uint8Array(40)))), /binary Parasolid.*\.x_t/);
});

await group('MicroStation DGN V7', async () => {
  // values asserted by GDAL's own test of this file (ogr_dgn.py) and by its DGN driver
  const g = await imp('smalltest.dgn', fx('dgnSmall')), el = g.polylines.find((p) => p.type === 15), sh = g.polylines.find((p) => p.type === 6), ln = g.polylines.find((p) => p.type === 3);
  ok(g.kind === 'polylines' && g.format === 'MicroStation DGN' && g.pathway === 'drawing' && g.stats.dimension === 2 && g.polylines.length === 3 && g.stats.elements.ellipse === 1 && g.stats.levels.join() === '2', 'DGN V7 2-D (GDAL autotest): line, shape and ellipse on level 2');
  ok(Math.abs(el.x[0] - 9.68780658389143) < 1e-9 && Math.abs(el.y[0] - 4.5835) < 1e-9 && el.closed && el.x.length === 48 && sh.closed && [4.5355, 4.3832, 4.9441, 4.832].every((v, k) => Math.abs(sh.x[k] - v) < 1e-12) && Math.abs(sh.y[2] - 2.5235) < 1e-12 && [2.5562, 5.7218, 2.5242, 6.0709].every((v, k) => Math.abs([ln.x[0], ln.y[0], ln.x[1], ln.y[1]][k] - v) < 1e-12), 'DGN: VAX doubles, middle-endian integers and working units decoded (coordinates as GDAL reports them)');
  const g3 = await imp('route3d.dgn', fx('dgn3d')), ch = g3.polylines.find((p) => p.type === 12);
  ok(g3.stats.dimension === 3 && g3.stats.units === 'm' && g3.stats.subUnit === 'mm' && g3.polylines.length === 3 && ch.x.length === 40 && ch.x[0] === 431000 && ch.y[0] === 6521000 && Math.abs(ch.z[39] + 151.2) < 1e-9 && g3.polylines.some((p) => p.closed && p.x.length === 4 && p.z[2] === -6), 'DGN V7 3-D (written by GDAL): global origin applied, complex chain of two line strings joined into 40 vertices');
  near(len3(g3), 3711.0667203, 1e-5, 'DGN 3-D: total length of the route, spur and shape as written');
  const h = await imp('housing.dgn', fx('dgnHousing')), shp = h.polylines.find((p) => p.type === 6);
  ok(h.polylines.length === 4 && h.polylines.every((p) => p.closed && p.z) && h.stats.elements['complex shape'] === 3 && h.stats.elements.arc === 9 && h.stats.elements.line === 9 && h.stats.levels.join() === '10,11', 'DGN (NIST repository): three complex shapes of lines, line strings and arcs joined into closed 3-D outlines');
  near(len3({ polylines: [shp] }), 450, 1e-9, 'DGN: perimeter of the 3-D shape element as reported by GDAL');
  ok(h.polylines.filter((p) => p.type === 14).every((p) => p.z.every((z) => Math.abs(z - p.z[0]) < 1e-6)) && h.polylines.filter((p) => p.type === 14).map((p) => p.z[0]).join() === '-62.5,12.5,12.5', 'DGN: 3-D arcs are turned by their quaternion into the plane of the complex shape they belong to');
  await throws('DGN V8 (OLE2 compound document)', () => importGeometry(F('v8.dgn', cat(le('u8', [0xd0, 0xcf, 0x11, 0xe0, 0xa1, 0xb1, 0x1a, 0xe1]), new Uint8Array(600)))), /V8.*V7/);
  await throws('DGN cell library', () => importGeometry(F('lib.dgn', cat(le('u8', [0x08, 0x05, 0x17, 0x00]), new Uint8Array(600)))), /cell library/);
});

await group('AutoCAD DWG', async () => {
  // counts, lengths and extents of the same drawing as read by ezdxf from its DXF source
  const g = await imp('route.dwg', fx('dwg2000')), by = (t) => g.polylines.filter((p) => p.type === t);
  ok(g.kind === 'polylines' && g.format === 'AutoCAD DWG' && g.pathway === 'drawing' && g.stats.version === 'AC1015' && g.stats.units === 'm' && g.stats.insunits === 6 && ['ROUTE', 'STRUCT'].every((l) => g.stats.layers.includes(l)), 'DWG 2000: version, $INSUNITS = 6 (m) and layer table');
  ok(g.stats.entities.LINE === 9 && g.stats.entities.LWPOLYLINE === 1 && g.stats.entities.POLYLINE === 1 && g.stats.entities.CIRCLE === 3 && g.stats.entities.ARC === 1 && g.stats.entities.ELLIPSE === 1 && g.stats.entities.INSERT === 2 && g.stats.entities.POINT === 1 && g.markers.length === 1, 'DWG: entity counts with two block references expanded');
  near(len3({ polylines: by('LINE') }), 270.062, 1e-3, 'DWG: total length of LINE entities (one in model space, eight from the scaled and rotated block)');
  near(len3({ polylines: by('POLYLINE') }), 170.543, 1e-3, 'DWG: 3-D polyline length through its vertex chain');
  ok(Math.abs(len3({ polylines: by('LWPOLYLINE') }) / 256.988 - 1) < 2e-3 && Math.abs(len3({ polylines: by('ARC') }) / 47.129 - 1) < 2e-3 && by('LWPOLYLINE')[0].z.every((z) => z === -104) && by('LWPOLYLINE')[0].layer === 'ROUTE', 'DWG: LWPOLYLINE bulges become arcs at the polyline elevation; ARC length within the facet tolerance');
  ok(Math.abs(g.bbox.min[0] + 10.2606) < 1e-3 && Math.abs(g.bbox.min[1] + 20) < 1e-9 && g.bbox.min[2] === -130 && Math.abs(g.bbox.max[0] - 416.1603) < 1e-3 && Math.abs(g.bbox.max[1] - 217.9904) < 1e-3, 'DWG: extents equal those of the DXF source');
  await throws('DWG 2007', () => importGeometry(F('a.dwg', cat('AC1021', new Uint8Array(300)))), /2007.*AC1021.*DXF/);
  await throws('DWG R12', () => importGeometry(F('a.dwg', cat('AC1009', new Uint8Array(300)))), /R11 \/ R12.*AC1009.*DXF/);
  await throws('DWG 2018 without sections', () => importGeometry(F('a.dwg', cat('AC1032', new Uint8Array(600)))), /truncated or corrupt/);
});

await group('Point clouds: E57 and LAZ', async () => {
  // sums and extents as read back by pye57 (pose applied) and laspy
  const e = await imp('scan.e57', fx('e57'));
  let sx = 0, sz = 0;
  for (let i = 0; i < e.points.length; i += 3) { sx += e.points[i]; sz += e.points[i + 2]; }
  ok(e.kind === 'points' && e.format === 'ASTM E57' && e.pathway === 'points' && e.count === 60 && e.stats.scans === 1 && e.stats.scanNames[0] === 'scan0', 'E57: one scan of 60 points');
  ok(Math.abs(sx - 59948.33731213721) < 1e-6 && Math.abs(sz + 8996.28059387207) < 1e-6 && Math.abs(e.bbox.min[0] - 980.4752958819739) < 1e-9 && Math.abs(e.bbox.max[1] - 2014.0070811086382) < 1e-9, 'E57: scan pose (rotation + translation) applied; coordinate sums equal those of libE57Format');
  await throws('E57 without scans', () => importGeometry(F('a.e57', cat('ASTM-E57', new Uint8Array(100)))), /E57/);
  for (const [key, n, sum, z0] of [['laz1', 150, [64652080.36, 978152248.44, -17761.348], -121.046], ['laz3', 200, [86204935.53, 1304203614.6, -23948.005], -122.696]]) {
    const g = await imp(key + '.laz', fx(key)), s = [0, 0, 0];
    for (let i = 0; i < g.points.length; i++) s[i % 3] += g.points[i];
    ok(g.kind === 'points' && g.format === 'LAZ (compressed LAS)' && g.count === n && g.stats.pointFormat === (key === 'laz1' ? 1 : 3) && /LASzip/.test(g.stats.compression) && s.every((v, k) => Math.abs(v - sum[k]) < 1e-4) && Math.abs(g.points[g.points.length - 1] - z0) < 1e-9, `LAZ point format ${g.stats.pointFormat}: ${n} arithmetic-decoded points, coordinate sums equal those of laspy`);
  }
  const viaLas = await importGeometry(F('cloud.las', fx('laz1')));
  ok(viaLas.count === 150 && viaLas.format === 'LAZ (compressed LAS)', 'compressed points under a .las name are passed to the LAZ reader');
  await throws('LAZ of LAS 1.4 point format 6', () => imp('f6.laz', fx('laz6')), /format 6.*layered.*laszip/);
});

await group('JT container summary', async () => {
  const seg = (id, type, body) => cat(le('u8', new Array(16).fill(id)), le('u32', [type, 24 + body]), new Uint8Array(body));
  const s1 = seg(1, 1, 40), s2 = seg(2, 7, 60), toc = cat(le('u32', [2]), le('u8', new Array(16).fill(1)), le('u32', [0, s1.length, 1 << 24]), le('u8', new Array(16).fill(2)), le('u32', [0, s2.length, 7 << 24]));
  const head = cat('Version 8.1 JT  DM 8.3.0.0'.padEnd(79, ' ') + '\n', le('u8', [0]), le('u32', [0, 105]), le('u8', new Array(16).fill(1)));
  await throws('JT file: message lists its segments', () => importGeometry(F('a.jt', cat(head, toc, s1, s2))), /JT 8\.1 file holds 2 segments.*logical scene graph.*shape LOD 0.*STEP AP242/);
  await throws('.jt that is not JT', () => importGeometry(F('a.jt', 'x')), /JT Open Toolkit/);
});

await group('convert-only formats', async () => {
  const conv = FORMATS.filter((f) => f.support === 'convert');
  let n = 0, good = 0;
  for (const f of conv) for (const ext of f.ext) {
    n++;
    try { await importGeometry(F('model.' + ext, cat('\x89HDF\r\n\x1a\n binary payload ', new Uint8Array(64)))); }
    catch (e) { if (e instanceof Error && e.message === f.convert) good++; else console.log(`     .${ext}: ${e && e.message}`); }
  }
  ok(n >= 40 && good === n,`all ${n} convert-only extensions (${conv.length} formats) throw their conversion instruction`);
  for (const ext of ['sldprt', 'fbx', 'CATPart', '3dm', 'pdf', 'parquet', '000', 'rvt', 'jt', 'bag', 'grib2', 'prt', 'ipt']) await throws(`.${ext}`, () => importGeometry(F('file.' + ext, 'x')), /Export|export|Convert|convert|Decompress|Save|save|use /);
  for (const ext of ['x_t', 'x_b', 'sat', 'sab', 'dwg', 'dgn', 'cgns', 'exo', 'med', 'h5', 'hdf5', 'nc4', 'gpkg', 'e57', 'laz']) await throws(`.${ext} holding something else`, () => importGeometry(F('file.' + ext, 'x'.repeat(300))), /^Not a|not a|Not an|is not|does not|neither|missing/);
  ok(['x_t', 'sat', 'dwg', 'dgn', 'cgns', 'exo', 'med', 'h5', 'gpkg', 'e57', 'laz'].every((e) => { const f = formatOf('a.' + e); return f.support === 'partial' && f.note.length > 150 && /not read|rejected|skipped|ignored/.test(f.note); }), 'formats read through the new container readers are registered as partial, with what is not read spelled out');
  await throws('unknown extension', () => importGeometry(F('file.qqq', 'x')), /Unsupported/);
  await throws('shapefile side-car selected alone', () => importGeometry(F('zone.dbf', 'x')), /\.shp/);
  await throws('empty file', () => importGeometry(F('a.stl', '')), /empty/);
  await throws('no file', () => importGeometry(null), /No file/);
  await throws('oversized file', () => importGeometry({ name: 'a.stl', size: 1e9, arrayBuffer: async () => new ArrayBuffer(0) }), /larger/);
});

await group('malformed and truncated input', async () => {
  const seen = new Set(), settle = (p) => Promise.race([p.then((g) => ['ok', g], (e) => ['err', e]), new Promise((r) => setTimeout(() => r(['hang']), 4000))]);
  const kinds = ['mesh', 'polylines', 'points', 'grid', 'voxels', 'network', 'table'];
  let cases = 0;
  for (const [name, bytes, opts] of SAMPLES) {
    const key = name.replace(/^[^.]*\./, '') + ':' + bytes.length;
    if (seen.has(key) || bytes.length < 8) continue;
    seen.add(key);
    const variants = [bytes.subarray(0, bytes.length >> 1), bytes.subarray(0, Math.floor(bytes.length / 3)), bytes.subarray(0, Math.min(bytes.length, 12)), bytes.subarray(0, bytes.length - 3)];
    for (const stride of [7, 31]) { const c = bytes.slice(); for (let i = stride; i < c.length; i += stride * 3) c[i] = (c[i] * 31 + 101) & 255; variants.push(c); }
    const flip = bytes.slice(); for (let i = 0; i < flip.length; i++) if (flip[i] >= 0x31 && flip[i] <= 0x38 && i % 5 === 0) flip[i] = 0x39; variants.push(flip);
    let bad = '';
    for (const v of variants) {
      cases++;
      const [st, r] = await settle(importGeometry(F(name, v), opts));
      if (st === 'hang') bad = 'hang';
      else if (st === 'err' && !(r instanceof Error && r.message.length > 0 && r.message.length < 500 && !/undefined|is not a function|Cannot read|Invalid typed array|out of bounds|RangeError|TypeError/.test(r.message))) bad = `raw error: ${r && r.message}`;
      else if (st === 'ok' && !(r && kinds.includes(r.kind) && Array.isArray(r.warnings))) bad = 'invalid geometry returned';
      if (bad) break;
    }
    ok(!bad, `damaged ${name} (${variants.length} variants) → clean result or clean error${bad ? ' — ' + bad : ''}`);
  }
  ok(cases > 500, `${cases} damaged inputs processed without a crash or hang`);
  await throws('STEP that is not STEP', () => importGeometry(F('a.step', 'hello world')), /ISO-10303-21/);
  await throws('XML bomb-style nesting', () => importGeometry(F('a.kml', '<a>'.repeat(5000))), /nested|No line|truncated/);
  await throws('garbage glTF', () => importGeometry(F('a.glb', cat('glTF', le('u32', [2, 9999, 500, 0x4e4f534a]), '{}'))), /truncated|glTF/);
  await throws('PLY with an absurd vertex count', () => importGeometry(F('a.ply', 'ply\nformat binary_little_endian 1.0\nelement vertex 4000000000\nproperty float x\nproperty float y\nproperty float z\nend_header\n')), /impossible|truncated|malformed/);
  await throws('NPY with an absurd shape', () => importGeometry(F('a.npy', npy('<f8', [100000, 100000, 100], new Uint8Array(8)))), /truncated|large|too many/);
  await throws('corrupt gzip', () => importGeometry(F('a.nii.gz', cat(le('u8', [0x1f, 0x8b, 8, 0, 0, 0, 0, 0, 0, 3]), new Uint8Array(40).fill(0xaa)))), /corrupt|truncated|NIfTI/);
  const bomb = `<model><resources>${Array.from({ length: 9 }, (_, k) => `<object id="${k}"><components>${Array.from({ length: 60 }, () => `<component objectid="${k + 1}"/>`).join('')}</components></object>`).join('')}</resources><build><item objectid="0"/></build></model>`;
  await throws('3MF component bomb', () => importGeometry(F('bomb.3mf', zip([['3D/3dmodel.model', bomb]]))), /too many|No triangles/);
  await throws('zip without directory', () => importGeometry(F('a.3mf', cat('PK\x03\x04', new Uint8Array(40)))), /zip/);
});

console.log(`\n${pass} checks passed, ${fails} failed, ${((Date.now() - t0) / 1000).toFixed(1)} s`);
process.exit(fails ? 1 : 0);
