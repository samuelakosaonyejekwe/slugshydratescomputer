// Geometry-import engine. One entry point (importGeometry) turns CAD, surface, mesh, drawing, GIS, point-cloud,
// voxel/image, pipeline-network, well-survey, seismic and plain numeric files into a single Geometry object, and a few
// analysis helpers (sectionOf, maskOf, gridOf, microstructure, networkSummary, dimensions, generate) make that object
// usable by every suite; js/core/route.js turns it into pipeline profiles, networks, terrain and wall maps.
// Everything is parsed locally; file content is treated as data only and is never executed.
import { LIMITS as IO_LIMITS, checkFile, parseSTL, parseOBJ, parseDXF, parseGeoJSON, parseDelimited, tableFromRows, sliceMesh, polylinesToSegments, rasterize } from './io.js';
import { rng } from './num.js';
import { openHDF5, hdf5Offset, isNumeric, typeName } from './hdf5.js';
import { openSQLite } from './sqlite.js';
import { parseDGN } from './fmt_dgn.js';
import { parseE57 } from './fmt_e57.js';
import { parseLAZ, isLAZ } from './fmt_laz.js';
import { recognisePipes, angularExtent, nurbsSurfPoint } from './fmt_brep.js';
import { parseSAT, parseSAB, isSAB } from './fmt_acis.js';
import { parseXT, parseXB, xtKind } from './fmt_xt.js';
import { parseDWG } from './fmt_dwg.js';
import { parseJT, isJT } from './fmt_jt.js';

export const GEOM_LIMITS = { ...IO_LIMITS, points: 2e6, voxels: 64e6, grid: 4e6, cells: 2e6, xmlNodes: 2e6, inflated: 512e6 };
const L = GEOM_LIMITS;

// ---- Format registry -----------------------------------------------------------------------------
const CAD_OUT = 'Export as STEP (AP214 or AP242) or STL from your CAD package.';
const F = (ext, name, pathway, support, note, convert) => Object.freeze({ ext: Object.freeze(ext.split(' ')), name, pathway, support, note, ...(convert ? { convert } : {}) });
const C = (ext, name, pathway, convert) => F(ext, name, pathway, 'convert', 'Closed or kernel-dependent format: not read directly.', convert);

export const FORMATS = Object.freeze([
  // CAD / solid
  F('step stp p21', 'STEP (ISO 10303)', 'cad', 'partial', 'AP203, AP214 and AP242 files: reads AP242 tessellated faces, faceted B-reps and ADVANCED_FACEs on planes, cylinders and cones (line, circle, ellipse and B-spline edges); other surface types are skipped and counted, and assembly placements are not applied.', 'For free-form surfaces export AP242 with tessellation, or STL.'),
  F('iges igs', 'IGES', 'cad', 'partial', 'Reads the curve wireframe (entities 100, 106, 110, 112, 116, 126 with their 124 transformations) projected onto its dominant plane; trimmed surfaces and solids are not tessellated.', 'For solids export STEP (AP214 or AP242) or STL.'),
  F('x_t x_b xmt_txt xmt_bin', 'Parasolid', 'cad', 'partial', 'Text transmit files (.x_t) that carry an embedded schema (the default since Parasolid 14; checked on files of versions 19 to 35), read through the base schema and the field edits embedded in each file: bodies with their region / shell / face / loop / fin / edge / vertex topology in metres; line, circle, ellipse and intersection-curve edges; faces on planes, cylinders and cones are tessellated, toroidal and spherical faces are filled over the parametric rectangle of their boundary. Cylinders and tori are also read as pipe runs and bends: a model recognised as a pipe, jumper or spool is returned as its 3-D centreline with stats.pipeRuns and stats.diameters (opts.prefer = "mesh" | "centreline" | "wireframe" overrides). B-spline, offset, swept, spun and blend surfaces are skipped and counted, B-curve edges become chords, and instance placements of assemblies are not applied. Binary files (.x_b) and text files without an embedded schema (modellers before version 14, or written against a fixed newer schema) are rejected.', 'Export a binary or schema-less file as Parasolid text (.x_t) from a current CAD package, and free-form shapes as STEP AP242 with tessellation or as STL.'),
  F('sat sab', 'ACIS SAT/SAB', 'cad', 'partial', 'Text ACIS files (.sat) of releases 1.x to 3x, including the ShapeManager flavour of Autodesk products: the entity table with body / lump / shell / face / loop / coedge / edge / vertex topology and the body transform; straight and elliptical edges; faces on planes, cylinders and cones are tessellated, toroidal and spherical faces are filled over the parametric rectangle of their boundary. Cylinders and tori are also read as pipe runs and bends (axis, length, diameters, bend radius and angle): a model recognised as a pipe, jumper or spool is returned as its 3-D centreline with stats.pipeRuns and stats.diameters (opts.prefer = "mesh" | "centreline" | "wireframe" overrides). Spline surfaces are skipped and counted, edges on procedural curves become chords, and a model without analytic faces is returned as its wireframe. Binary .sab files are not read.', 'Save a binary .sab model as text ACIS (.sat), and export free-form shapes as STEP AP242 with tessellation or as STL.'),
  F('jt', 'JT', 'cad', 'convert', 'Recognised, not read: the header, table of contents and segment headers of JT 8, 9 and 10 files are parsed, so the error message lists what the file holds (scene graph, shape LODs, B-rep, PMI), but the tessellation itself sits in the JT codecs (bit-length, Huffman and arithmetic coded index streams, quantised or topologically compressed vertices) and the exact shape in JT or Parasolid B-rep segments, none of which is decoded.', 'JT keeps its tessellation in proprietary-style compression codecs (ISO 14306) that are not decoded here. Export the JT model as STEP AP242 or STL from NX or Teamcenter Visualization, or convert it with the JT Open Toolkit.'),
  C('3dxml', '3DXML', 'cad', 'Export from CATIA / 3DEXPERIENCE as STEP AP214/AP242 or STL.'),
  C('prc', 'PRC', 'cad', 'Export the PRC model as STEP or STL from the authoring CAD package.'),
  C('vda', 'VDA-FS', 'cad', 'Convert VDA-FS to STEP or IGES in your CAD package.'),
  C('cgr', 'CATIA CGR', 'cad', 'Export from CATIA as STEP AP214/AP242 or STL.'),
  C('catpart catproduct model', 'CATIA V4/V5', 'cad', 'Export from CATIA as STEP AP214/AP242 or STL.'),
  C('sldprt sldasm', 'SOLIDWORKS', 'cad', 'In SOLIDWORKS use Save As: STEP AP214/AP242 or STL.'),
  C('ipt iam', 'Autodesk Inventor', 'cad', 'In Inventor use Export: CAD Format, STEP or STL.'),
  C('prt asm', 'Creo / NX / Solid Edge part or assembly', 'cad', CAD_OUT),
  C('par psm', 'Solid Edge', 'cad', 'In Solid Edge use Save As: STEP or STL.'),
  C('3dm', 'Rhino 3DM', 'cad', 'In Rhino use Export: STEP, STL, OBJ or PLY.'),
  C('rvt rfa', 'Revit', 'cad', 'Export from Revit as IFC (equipment and connectivity) or as STL / DXF (shape).'),
  C('f3d f3z', 'Fusion design', 'cad', 'In Fusion use Export: STEP or STL.'),
  C('_pd xcgm', 'CADDS / XCGM', 'cad', CAD_OUT),
  // Surface / tessellated
  F('stl', 'STL', 'surface', 'full', 'ASCII and binary triangle facets.'),
  F('obj', 'Wavefront OBJ', 'surface', 'full', 'Vertices and polygonal faces (fan-triangulated); materials and curves are ignored.'),
  F('ply', 'PLY', 'surface', 'full', 'ASCII and binary (little/big endian); faces become triangles, vertex-only files become a point cloud.'),
  F('off', 'OFF', 'surface', 'full', 'OFF / COFF / NOFF polygon lists, fan-triangulated.'),
  F('3mf', '3MF', 'surface', 'full', 'Mesh objects of the core specification with build-item and component transforms.'),
  F('amf', 'AMF', 'surface', 'full', 'Plain or zip-wrapped XML; straight-edged triangles (curved-edge data is ignored).'),
  F('wrl vrml iv', 'VRML / Open Inventor', 'surface', 'partial', 'IndexedFaceSet coordinates and indices of text files; Transform nodes are not applied.', 'Export as STL, OBJ or PLY if the model relies on transforms.'),
  F('x3d', 'X3D', 'surface', 'partial', 'IndexedFaceSet / IndexedTriangleSet of the XML encoding; Transform nodes are not applied.', 'Export as STL, OBJ or PLY if the model relies on transforms.'),
  F('dae', 'COLLADA', 'surface', 'partial', 'triangles / polylist / polygons with POSITION accessors; scene-node transforms are not applied.', 'Export as STL, OBJ or glTF if the model relies on node transforms.'),
  F('gltf glb', 'glTF 2.0', 'surface', 'partial', 'POSITION accessors and indices of triangle primitives with node transforms; embedded, GLB or companion buffers. Draco / meshopt compression, sparse accessors, skins and morph targets are not read.', 'Re-export the glTF without Draco or meshopt compression.'),
  C('fbx', 'Autodesk FBX', 'surface', 'Convert FBX to glTF/GLB, OBJ or STL (for example with Blender or the Autodesk FBX converter).'),
  F('vtk', 'VTK legacy', 'mesh', 'full', 'ASCII and binary POLYDATA, UNSTRUCTURED_GRID, STRUCTURED_GRID, RECTILINEAR_GRID and STRUCTURED_POINTS.'),
  F('vtp vtu vti vts vtr', 'VTK XML', 'mesh', 'full', 'ASCII, inline base64 and appended (raw or base64) arrays, uncompressed or zlib-compressed; parallel (.pvtu …) index files are not followed.'),
  F('byu', 'BYU surface', 'surface', 'full', 'Movie.BYU polygon parts in free-format ASCII.'),
  F('gts', 'GTS', 'surface', 'full', 'GNU Triangulated Surface vertex / edge / face lists.'),
  F('tri fac surf', 'Generic ASCII node-face file', 'surface', 'partial', 'Reads the common "counts, vertex rows, triangle rows" layouts; other dialects of these extensions are rejected.', 'Export the surface as STL, OBJ or PLY.'),
  // Computational meshes
  F('msh', 'Gmsh MSH', 'mesh', 'partial', 'ASCII versions 2.2, 4.0 and 4.1 (binary files are rejected); a .msh that is a Fluent ASCII mesh is recognised and read instead.', 'Save the mesh from Gmsh as ASCII (Mesh.Binary = 0), version 2.2 or 4.1.'),
  F('cas msh', 'Fluent mesh (ASCII)', 'mesh', 'partial', 'ASCII node (10) and face (13) sections; boundary faces are those with one neighbour cell. Binary sections are rejected.', 'Write the case/mesh file from Fluent in ASCII, or export CGNS-free formats such as Gmsh .msh, VTK .vtu, SU2 or UNV.'),
  F('su2', 'SU2 mesh', 'mesh', 'full', 'NDIME / NELEM / NPOIN sections of the native ASCII format.'),
  F('unv', 'I-DEAS Universal', 'mesh', 'full', 'Datasets 2411 (nodes) and 2412 (elements); parabolic elements are reduced to their corner nodes.'),
  F('bdf nas nastran', 'Nastran bulk data', 'mesh', 'partial', 'GRID, CTRIA3/6, CQUAD4/8, CTETRA, CHEXA, CPENTA and CPYRAM in small, large and free field; local coordinate systems are not applied.'),
  F('inp', 'Abaqus input deck', 'mesh', 'partial', '*NODE and *ELEMENT blocks: C3D4/5/6/8/10/15/20 solids, S3/S4/S8R, STRI, M3D, CPS/CPE/CAX and R3D shells and planes, PIPE/B/T/ELBOW line elements (a deck of line elements only becomes 3-D polylines). Second-order elements are reduced to their corner nodes; *INCLUDE files, *NGEN/*NFILL/*ELGEN generators, *SYSTEM transforms and instance placements of *ASSEMBLY are not applied.', 'Write a flat (non-parts) input deck with all nodes and elements listed explicitly (Abaqus/CAE: model attribute "Do not use parts and assemblies in input files").'),
  F('cdb', 'ANSYS CDB archive', 'mesh', 'partial', 'Blocked NBLOCK / EBLOCK (SOLID format) written by CDWRITE: solid, shell, plane and beam/pipe elements classified by their ET type, degenerate tetrahedra, wedges and pyramids recognised, mid-side nodes dropped. Unblocked N / EN commands, the non-SOLID EBLOCK layout and local coordinate systems are not read.', 'Write the archive with CDWRITE,DB,file,cdb with blocked output (CDOPT default, ,,,,BLOCKED).'),
  F('k key dyn', 'LS-DYNA keyword deck', 'mesh', 'partial', '*NODE, *ELEMENT_SHELL, *ELEMENT_SOLID (one- and two-card forms, degenerate tetrahedra / wedges / pyramids), *ELEMENT_TSHELL and *ELEMENT_BEAM in fixed, comma-separated and LONG=Y format; option cards (thickness, beta, ortho) are skipped. *INCLUDE files, *NODE_TRANSFORM and *DEFINE_TRANSFORMATION are not applied.', 'Merge the includes into one deck (LS-PrePost: File, Save Keyword As with "include files merged").'),
  F('dat tec plt', 'Tecplot ASCII', 'mesh', 'partial', 'Finite-element zones (triangle, quadrilateral, tetrahedron, brick) and ordered I/J/K zones with nodal coordinates; binary .plt is not read.', 'Save the data set from Tecplot as ASCII (.dat).'),
  F('p3d x', 'Plot3D grid (ASCII)', 'mesh', 'partial', 'Formatted single- or multi-block structured grids without IBLANK; a .xyz file is read as Plot3D only when its header looks like one.', 'Write the Plot3D grid formatted (ASCII), whole, without IBLANK.'),
  F('foam', 'OpenFOAM polyMesh', 'mesh', 'partial', 'ASCII points + faces + owner (+ neighbour) files passed together (opts.companion); the boundary surface is built from faces without a neighbour.', 'Select the polyMesh "points" file together with "faces", "owner" and "neighbour" (ASCII; run foamFormatConvert first if binary).'),
  F('cgns', 'CGNS', 'mesh', 'partial', 'HDF5-flavoured CGNS files: every zone of every base with its GridCoordinates (Cartesian, or cylindrical R / Theta converted); structured zones give their block boundaries, unstructured Elements_t sections of BAR, TRI, QUAD, TETRA, PYRA, PENTA and HEXA elements (higher-order ones reduced to their corner nodes), MIXED sections and NGON_n / NFACE_n polyhedra (reduced to their boundary faces) give the boundary surface, outline or line work. Internal links are followed. ADF-flavoured files, links to other files, solutions, boundary-condition patches and zone connectivity are not read.', 'Convert an ADF-flavoured file to HDF5 with the CGNS tools: cgnsconvert -h in.cgns out.cgns (or adf2hdf in.cgns out.cgns).'),
  F('e exo ex2 exii g', 'Exodus II', 'mesh', 'partial', 'Classic NetCDF (CDF-1, 64-bit-offset CDF-2, CDF-5) and NetCDF-4 / HDF5 files: nodal coordinates (coordx / coordy / coordz or coord) and every connect<N> element block typed by its elem_type attribute (HEX, TETRA, WEDGE, PYRAMID, QUAD / SHELL, TRI, BAR / BEAM / TRUSS; higher-order elements reduced to their corner nodes). Polyhedral NSIDED / NFACED blocks, sphere elements, node and side sets and result variables are not read.', 'For a mesh that relies on polyhedral blocks convert it to VTK .vtu: meshio convert in.exo out.vtu, or ParaView File > Save Data.'),
  F('med rmed', 'Salome MED', 'mesh', 'partial', 'MED 2.x to 4.x (HDF5): the node coordinates and element connectivity of one unstructured mesh (opts.mesh = name, default the first; the first computation step): SE2, TR3, QU4, TE4, PY5, PE6, HE8, their quadratic variants reduced to corner nodes, and polygons (POG). Structured grids, polyhedra (POE), families / groups and result fields are not read.', 'For a structured or polyhedral MED mesh export UNV or Gmsh .msh from Salome, or run meshio convert in.med out.vtu.'),
  C('case', 'EnSight Gold', 'mesh', 'EnSight cases span several binary files. Export as VTK .vtu, Gmsh .msh or Tecplot ASCII.'),
  // 2-D drawings
  F('dxf', 'DXF', 'drawing', 'partial', 'ASCII DXF: LINE, LWPOLYLINE, POLYLINE, CIRCLE, ARC, ELLIPSE, SPLINE and 3DFACE; bulges, blocks/inserts and hatches are ignored.', 'Save the drawing as ASCII DXF with blocks exploded.'),
  F('svg', 'SVG', 'drawing', 'partial', 'path, polygon, polyline, rect, circle, ellipse and line with nested transform attributes; curves are flattened and the y axis is flipped to point up. <use> references, viewBox scaling, CSS transforms and text are ignored.', 'Convert text and <use> clones to plain paths (Inkscape: Object to Path) before exporting.'),
  F('hpgl hpg plt', 'HPGL plot file', 'drawing', 'partial', 'Pen moves PU / PD / PA / PR plus CI circles and AA arcs, converted from plotter units to millimetres.'),
  F('xy', 'x-y polyline', 'drawing', 'full', 'Two numeric columns per row; blank lines separate polylines.'),
  F('dwg', 'AutoCAD DWG', 'drawing', 'partial', 'AutoCAD R13, R14, 2000, 2004, 2010, 2013 and 2018 drawings (AC1012, AC1014, AC1015, AC1018, AC1024, AC1027, AC1032): model-space LINE, LWPOLYLINE and POLYLINE (bulges as arcs, 2-D and 3-D), ARC, CIRCLE, ELLIPSE, SPLINE, POINT and 3DFACE, block references (INSERT / MINSERT with scale, rotation and extrusion, nested), layer names, 3-D coordinates, and the drawing unit $INSUNITS of 2000 drawings. Text, dimensions, hatches, leaders, solids and ACIS bodies, meshes, attributes, paper space and invisible entities are skipped; the drawing unit of 2004 and later drawings is not read. AutoCAD 2007 drawings (AC1021, a Reed-Solomon coded layout of their own) and releases before R13 are rejected.', 'Save a 2007 or pre-R13 drawing as ASCII DXF (AutoCAD: SAVEAS, DXF; or the free ODA File Converter, output "ASCII DXF"), or re-save it in the 2010 or a later DWG format.'),
  F('dgn', 'MicroStation DGN', 'drawing', 'partial', 'MicroStation V7 design files (ISFF), 2-D and 3-D: lines, line strings, shapes, curves, arcs and ellipses (flattened), complex chains and shapes joined into one polyline, B-spline curves as their control polygon, point strings; coordinates in master units from the design-file header with the global origin applied, element levels kept, cell components read in place. Text, tags, dimensions, surfaces and solids are skipped; the placement of shared-cell instances is not decoded, so a shared-cell definition is shown once in its own coordinates. V8 design files (OLE2 compound documents with an unpublished element stream) and cell libraries are rejected with an explanation.', 'Save a V8 design file as V7 (MicroStation: File > Save As, "MicroStation V7 DGN") or export it as DXF.'),
  C('eps ai pdf', 'EPS / AI / PDF vector drawing', 'drawing', 'Convert the vector drawing to SVG or DXF (for example with Inkscape).'),
  C('gbr ger', 'Gerber', 'drawing', 'Export the layer as DXF or SVG from your PCB/CAM tool.'),
  C('idf emn', 'IDF board outline', 'drawing', 'Export the outline as DXF or STEP.'),
  // GIS / terrain
  F('geojson', 'GeoJSON', 'gis', 'full', 'Polygon, LineString and Point geometries of any nesting.'),
  F('kml kmz', 'KML / KMZ', 'gis', 'full', 'Polygon rings, LineStrings and Points (longitude, latitude).'),
  F('gpx', 'GPX', 'gis', 'full', 'Tracks, routes and waypoints (longitude, latitude, elevation).'),
  F('gml', 'GML', 'gis', 'partial', 'posList / pos / coordinates of rings, line strings and points; axis order is taken as written.'),
  F('shp', 'ESRI Shapefile', 'gis', 'full', 'Point, MultiPoint, PolyLine and Polygon including Z and M variants (MultiPatch records are skipped); .dbf field names and .prj text are summarised when supplied as companions.'),
  F('asc', 'ESRI ASCII grid', 'gis', 'full', 'Header + row-major elevations; an .asc without that header is read as an ASCII point list.'),
  F('grd', 'Surfer grid', 'gis', 'full', 'Surfer ASCII (DSAA), Surfer 6 binary (DSBB) and Surfer 7 binary (DSRB); NetCDF .grd files are passed to the NetCDF reader.'),
  F('tif tiff', 'GeoTIFF / TIFF', 'gis', 'partial', 'Classic TIFF, strips or tiles, 1/8/16/32-bit integer and 32/64-bit float, uncompressed, Deflate, LZW or PackBits with predictors 1-3; a georeferenced or float page becomes a DEM, 8/16-bit pages become a 2-D or (multi-page) 3-D voxel image. BigTIFF and JPEG-in-TIFF are not read.', 'Rewrite the raster as a classic (non-BigTIFF) TIFF with Deflate or LZW compression, e.g. gdal_translate -co BIGTIFF=NO -co COMPRESS=DEFLATE.'),
  F('nc cdf', 'NetCDF classic', 'gis', 'partial', 'CDF-1, CDF-2 and CDF-5 files: the first 2-D (or first slice of a higher-dimensional) numeric variable with its coordinate variables becomes a grid. A NetCDF-4 file under these extensions is passed to the NetCDF-4 reader.'),
  F('nc4', 'NetCDF-4 (HDF5)', 'gis', 'partial', 'NetCDF-4 files (HDF5 container, also under .nc): the first 2-D (or first slice of a higher-dimensional) numeric variable with the coordinate variables of its dimension scales becomes a grid, with scale_factor / add_offset / _FillValue applied; contiguous, chunked, deflated and shuffled storage. Without a gridded variable the file is read as plain HDF5. Szip and third-party compression filters, user-defined compound variables and groups other than through their full path are not read.', 'Rewrite a file compressed with another filter: nccopy -d 4 in.nc out.nc (deflate), or nccopy -k classic in.nc out.nc.'),
  F('bil bip bsq', 'ENVI / ESRI band raster', 'gis', 'full', 'First band of a BIL/BIP/BSQ raster described by its .hdr companion (opts.companion).', 'Supply the .hdr header file together with the raster.'),
  F('mif', 'MapInfo MIF', 'gis', 'partial', 'REGION, PLINE, LINE, RECT and POINT objects.'),
  F('landxml', 'LandXML', 'gis', 'partial', 'LandXML 1.x (also recognised in .xml): Alignments and PlanFeatures CoordGeom (Line, Curve, IrregularLine; Spiral as its chord) as polylines with elevations from Profile / ProfAlign PVIs or 3-D points, Surfaces (Pnts + Faces TIN) as a mesh, CgPoints as points and PipeNetworks (Structs + Pipes) as a network. One content class is returned per import (opts.prefer = "alignment" | "surface" | "points" | "network"); vertical curves are joined by straight grades and cross-sections, parcels and grade models are ignored.'),
  F('dem dtm dsm', 'DEM / DTM / DSM elevation model', 'gis', 'partial', 'Text elevation models only: an ESRI ASCII grid, a USGS ASCII DEM (record A + profiles, resampled onto a regular grid) or rows of x y z. Binary DEM flavours (SRTM .hgt, BIL, ERDAS, DTED, SDTS) are rejected.', 'Convert a binary elevation model to GeoTIFF or ESRI ASCII grid: gdal_translate -of AAIGrid in.dem out.asc (or -of GTiff -co COMPRESS=DEFLATE).'),
  F('sgy segy', 'SEG-Y seismic', 'points', 'partial', 'Textual header (EBCDIC or ASCII), binary header and per-trace CDP / source coordinates with the coordinate scalar become trace positions; water depth at source gives z. For IBM-float (1) and IEEE-float (5) samples of small files the first strong arrival of each trace is picked as a seabed two-way time and converted with 1500 m/s: an estimate, not an interpreted horizon. Other sample formats give positions only; SEG-Y rev 2 extended trace headers are skipped.', 'For an interpreted seabed or horizon export it from the interpretation package as XYZ, ESRI ASCII grid or GeoTIFF.'),
  F('las', 'LAS point cloud', 'points', 'full', 'LAS 1.0-1.4 uncompressed, point formats 0-10, scaled to real coordinates and evenly sub-sampled above the point limit.'),
  F('laz', 'LAZ (compressed LAS)', 'points', 'partial', 'LASzip-compressed LAS 1.0-1.4 with point formats 0-3: the arithmetic-coded "pointwise" and "pointwise chunked" compressors with the version-2 POINT10, GPSTIME11, RGB12 and extra-byte codecs; x, y, z scaled to real coordinates and evenly sub-sampled above the point limit (of very large clouds only every n-th chunk is decompressed). LAS 1.4 point formats 6-10 (layered compression), wave-packet formats 4-5, variable-size chunks and files of LASzip releases before 2.0 are rejected.', 'Decompress such a file to LAS: laszip -i in.laz -o out.las, or pdal translate in.laz out.las.'),
  F('e57', 'ASTM E57', 'points', 'partial', 'Point coordinates of every scan: the XML section and the CompressedVector binary sections with the bit-pack codec for float, double, integer and scaled-integer fields; Cartesian coordinates, or spherical ones converted; points flagged invalid are dropped, each scan pose is applied and clouds above the point limit are evenly sub-sampled. Images, intensity, colour and other per-point fields are ignored and page checksums are not verified.'),
  F('gpkg', 'GeoPackage', 'gis', 'partial', 'Vector feature tables of the SQLite container read in place: gpkg_contents / gpkg_geometry_columns / gpkg_spatial_ref_sys, GeoPackage binary headers and WKB Point, LineString, Polygon, their Multi and collection forms with Z / M (ISO or EWKB codes, either byte order); every feature layer is merged (opts.layer = table name reads one) and polylines carry their layer, name and attribute values. Circular-arc curve types are reduced to their control points. Tile pyramids and gridded-coverage elevation rasters, the RTree index and extensions are not read, coordinates are not re-projected, and changes still in a write-ahead log (-wal file) are not replayed.', 'Export a raster layer as GeoTIFF: gdal_translate in.gpkg out.tif. Checkpoint a database in WAL mode first (sqlite3 in.gpkg "PRAGMA wal_checkpoint(TRUNCATE)").'),
  C('grib grb grib2 grb2', 'GRIB / GRIB2', 'gis', 'Convert the GRIB field to NetCDF classic or GeoTIFF (cdo -f nc copy, or gdal_translate).'),
  C('000', 'S-57 / S-101 chart', 'gis', 'Export the chart objects as Shapefile or GeoJSON and soundings as XYZ (ogr2ogr reads S-57).'),
  C('bag', 'BAG bathymetry', 'gis', 'BAG is an HDF5 container. Convert the surface to GeoTIFF or XYZ (gdal_translate reads BAG).'),
  C('tab mid', 'MapInfo TAB / MID', 'gis', 'Export the table as MapInfo MIF, Shapefile or GeoJSON.'),
  // Point clouds
  F('pts', 'PTS point cloud', 'points', 'full', 'Optional count line then x y z [intensity r g b] rows.'),
  F('ptx', 'PTX scan', 'points', 'full', 'One or more scans with their 4×4 registration transform applied; empty returns are dropped.'),
  F('xyz xyzi xyzrgb', 'XYZ points / soundings', 'points', 'full', 'x y z rows with optional header and extra columns; a Plot3D-style header switches to the structured-grid reader.'),
  // Image / voxel / microstructure
  F('png jpg jpeg bmp webp gif', 'Raster image', 'voxel', 'partial', 'Decoded by the browser, converted to greyscale and thresholded (given value or Otsu) into a 2-D solid/pore mask. Not available outside a browser.'),
  F('raw vol bin', 'RAW voxel volume', 'voxel', 'partial', 'Headerless volume with opts.dims and opts.dtype; a uint8 file whose size is a perfect cube or square is inferred.'),
  F('nrrd nhdr', 'NRRD', 'voxel', 'full', 'raw, ascii and gzip encodings; detached data through opts.companion.'),
  F('mha mhd', 'MetaImage', 'voxel', 'full', 'Inline (.mha) or companion (.mhd + raw) data, optionally zlib-compressed, single channel.'),
  F('nii nii.gz', 'NIfTI-1', 'voxel', 'full', 'Single-file .nii and gzip-compressed .nii.gz volumes (first 3-D volume of 4-D data).'),
  F('npy npz', 'NumPy array', 'numeric', 'full', 'Format versions 1-3, bool / integer / float types in C or Fortran order; 2-D and 3-D arrays become voxels, N×2 and N×3 arrays become points; .npz uses its first array.'),
  F('dcm dicom', 'DICOM slice', 'voxel', 'partial', 'Uncompressed little-endian (explicit or implicit VR) single- or multi-frame monochrome images.', 'Decompress the DICOM file (gdcmconv --raw, or dcmdjpeg) or export the stack as TIFF / NRRD / NIfTI.'),
  C('jp2 j2k jpx', 'JPEG 2000', 'voxel', 'Convert JPEG 2000 images to PNG or TIFF.'),
  F('h5 hdf5 hdf he5 h4', 'HDF5', 'numeric', 'partial', 'HDF5 files read in place (superblock versions 0-3, old- and new-style groups, contiguous / compact / chunked data with every chunk index, deflate and shuffle filters, integer, floating-point, string, enum, array and compound types). The content decides the reading: a CGNS tree, a MED mesh, an Exodus II or NetCDF-4 data set and a MATLAB v7.3 file go to those readers; otherwise the most geometry-like numeric data set (opts.dataset = path overrides) becomes points or a polyline (N × 2, N × 3), a grid or table (2-D), voxels (3-D) or a table of equal-length vectors, and g.contents lists what the file holds. HDF4 files, external and virtual data sets, szip / n-bit / scale-offset and third-party filters are not read.', 'Rewrite data compressed with another filter: h5repack -f GZIP=4 in.h5 out.h5; convert HDF4 with h4toh5.'),
  // Plant / piping / network topology
  F('dev wbt survey', 'Well deviation survey (ASCII)', 'well', 'partial', 'Column text with a header naming MD + inclination + azimuth (positions by the minimum-curvature method) or MD + TVD (+ northing / easting offsets); common aliases and feet are recognised, and a header-less file is taken as MD, inclination, azimuth. The same tables are detected in .csv / .txt. Casing, tubing and completion records and binary survey databases are not read.', 'Export the survey from the well-planning package as ASCII or CSV with MD, inclination and azimuth columns.'),
  F('graphml', 'GraphML', 'network', 'partial', 'Nodes and edges with their <data> keys (x, y, z / elevation, type, name, length, diameter) and yEd node geometry; nested graphs are flattened, hyperedges and ports are skipped.'),
  F('ifczip', 'IFCZIP', 'network', 'partial', 'Zip archive holding one IFC STEP file, read exactly like .ifc (equipment inventory, connectivity, pipe lengths; no solid geometry). An archive that holds only ifcXML is rejected.', 'Save the model as IFC (STEP physical file, .ifc) rather than ifcXML.'),
  F('pcf', 'Piping Component File', 'network', 'partial', 'Components with END-POINT / CENTRE-POINT / BRANCH1-POINT become a centre-line network in metres; material lists and supports are not interpreted.'),
  F('ifc', 'IFC (STEP text)', 'network', 'partial', 'Builds an equipment inventory (pumps, tanks, valves, pipe segments, fittings, exchangers, filters, meters …) with placements, port connectivity and pipe lengths from quantities or extrusions; no solid geometry is tessellated.', 'For the building shape export STL or OBJ from the BIM tool.'),
  F('aml', 'AutomationML / CAEX', 'network', 'partial', 'InternalElement hierarchy as nodes and InternalLink as edges; attributes other than names and roles are ignored.'),
  F('yaml yml', 'YAML definition', 'network', 'partial', 'A safe subset (maps, lists, scalars, inline lists/maps) holding a network, mesh, points or parameter table.'),
  // Numeric
  F('json', 'JSON data', 'numeric', 'full', 'Sniffed as GeoJSON, glTF, network (nodes + edges), mesh (nodes + elements / vertices + faces), array (grid, voxels, points) or record table.'),
  F('xml', 'XML data', 'numeric', 'partial', 'Routed by root element to LandXML, GraphML, KML, GPX, GML, VTK, COLLADA, X3D, AMF, SVG or CAEX; otherwise read as a network (nodes + edges) definition.'),
  F('csv tsv txt', 'Delimited / column text', 'numeric', 'full', 'Sniffed as node-edge table, well survey (MD / inclination / azimuth or MD / TVD), chainage-elevation profile, 3-D route, wall-thickness / corrosion / deposit map (x, θ, value[, time]), particle cloud, lon/lat/depth soundings, x y z points, x y polyline, numeric matrix, time series or a plain table.'),
  F('mat', 'MATLAB MAT-file', 'numeric', 'partial', 'Level 5 MAT-files (v5, v6, v7; plain or zlib-compressed) and v7.3 MAT-files (HDF5): real numeric and logical arrays. Equal-length vectors become one table with the variable names as headers, a matrix becomes a table (up to 16 columns), a grid or a voxel volume. Structs, cells, sparse and character arrays and imaginary parts are skipped; v4 files are rejected.', 'Write struct or cell contents as plain numeric variables, or as CSV (writematrix / writetable).'),
  C('parquet', 'Apache Parquet', 'numeric', 'Export the table as CSV or JSON (pandas: read_parquet(...).to_csv(...)).'),
]);

export const PATHWAYS = Object.freeze({
  cad: { title: 'CAD / solid geometry', blurb: 'Exact B-rep exchange files (STEP AP203 / AP214 / AP242, IGES) and the kernel text formats Parasolid X_T and ACIS SAT for trees, manifolds, jumpers, spools and vessels; a pipe, jumper or spool in X_T or SAT is recognised from its cylinders and bends and returned as a centreline with diameters. JT and native CAD files must be exported to STEP or STL first; a long pipeline is better described by its centreline than by a solid.' },
  surface: { title: 'Surface / tessellated geometry', blurb: 'Faceted approximations of equipment, pipe walls and deposit surfaces: STL, OBJ, PLY, OFF, 3MF, AMF, VRML/X3D, COLLADA, glTF, GTS, BYU.' },
  mesh: { title: 'CFD and structural meshes', blurb: 'Volume, shell and beam meshes from CFD and FEA (Gmsh, VTK, CGNS, Exodus II, Salome MED, SU2, UNV, Nastran, Abaqus, ANSYS, LS-DYNA, Fluent, Tecplot, Plot3D, OpenFOAM); the boundary surface, outline or pipe centreline is extracted.' },
  drawing: { title: '2-D drawings and profiles', blurb: 'DXF, DWG (R13 to 2018, except 2007), MicroStation DGN V7, SVG, HPGL and x-y polylines for alignment sheets, elevation profiles, field layouts and cross-sections; 3-D line work keeps its elevations.' },
  gis: { title: 'GIS, terrain and bathymetry', blurb: 'Pipeline routes and seabed: GeoJSON, KML, GPX, GML, Shapefile, LandXML, ASCII / Surfer grids, DEM / DTM, GeoTIFF, NetCDF (classic and NetCDF-4), GeoPackage vector layers, band rasters, MIF. A route can be draped over a bathymetry grid to give its elevation profile.' },
  points: { title: 'Point clouds, soundings and seismic lines', blurb: 'LAS, LAZ, E57, PTS, PTX, XYZ and CSV points from bathymetric surveys, laser scans and inspection; SEG-Y trace positions with an estimated seabed. Gridded on demand.' },
  voxel: { title: 'Image and voxel data', blurb: 'Images, TIFF stacks, RAW, NRRD, MetaImage, NIfTI, DICOM and NumPy volumes thresholded into solid / void voxels: deposit and plug scans, sand packs, core samples.' },
  network: { title: 'Pipeline, equipment and network topology', blurb: 'Node-edge definitions of gathering networks, wells, manifolds, valves, chokes and sensors (JSON, YAML, XML, GraphML, CSV), PCF piping, IFC / IFCZIP equipment and AutomationML hierarchies.' },
  well: { title: 'Wells and trajectories', blurb: 'Deviation surveys (MD + inclination + azimuth by minimum curvature, or MD + TVD) from ASCII, CSV, JSON or MATLAB tables become 3-D well paths with dog-leg severity.' },
  numeric: { title: 'Tabular and scientific data', blurb: 'CSV, text, JSON, XML, YAML, NumPy, MATLAB (incl. v7.3) and HDF5 data sniffed into profiles (chainage / elevation), 3-D routes, wall-thickness, corrosion and deposit maps, particle clouds, time series, grids, networks or parameter tables.' },
  procedural: { title: 'Parametric and procedural geometry', blurb: 'Generated without a file: catenary and lazy-wave risers, undulating flowlines, build-and-hold wells, M-shaped jumpers, and porous structures (sphere packs, Voronoi foams, TPMS, lattices, CSG trees, implicit functions).' },
});

export const SUITE_GEOMETRY = Object.freeze({
  pvt: { classes: 'Mainly tabular data: CSV, TXT/DAT, JSON, XML, YAML, XLSX-exported CSV, NetCDF and MAT tables of composition, PVT, phase-equilibrium, hydrate-curve and inhibitor data. STEP, IGES, STL and OBJ are accepted from other modules but drive no thermodynamic calculation', accepts: ['table', 'params', 'grid'] },
  net: { classes: 'STEP, IGES, STL, OBJ, 3MF, glTF, DXF, IFC; meshes (VTK, Gmsh, Nastran, Abaqus, Fluent, OpenFOAM, UNV); pipeline centrelines (XYZ, chainage / elevation, easting / northing / elevation, lat / lon / depth), node-edge networks (CSV, JSON, XML, GraphML), well surveys (MD / TVD, MD / inclination / azimuth), GIS routes (Shapefile, GeoJSON, KML/KMZ, LandXML) and bathymetry (XYZ, DEM/DTM, GeoTIFF, NetCDF, ESRI ASCII, LAS, SEG-Y)', accepts: ['polylines', 'network', 'points', 'grid', 'mesh', 'table', 'params'] },
  flow: { classes: '1-D pipeline and network data first: CSV / XYZ centrelines, chainage or distance versus elevation, MD / TVD, node-edge networks (CSV, JSON, XML, GraphML), GIS polylines (Shapefile, GeoJSON, KML/KMZ, LandXML) draped on DEM / GeoTIFF / NetCDF / XYZ bathymetry; STEP, STL and CFD meshes (VTK, Gmsh, Fluent, OpenFOAM, UNV) for local detail', accepts: ['polylines', 'network', 'points', 'grid', 'table', 'mesh', 'params'] },
  solids: { classes: 'The pipe and network geometry of the flow suite plus deposit data: CSV / XYZ deposition maps (x, θ, t, δ), structured grids, voxel data, surface meshes (STL, OBJ, PLY, VTP) and CFD meshes (VTK, Gmsh, Fluent, OpenFOAM, UNV), particle and point clouds with size, velocity and density', accepts: ['table', 'grid', 'points', 'voxels', 'mesh', 'polylines', 'network', 'params'] },
  ops: { classes: 'Network topology and equipment connectivity: CSV, JSON, XML, YAML and GraphML node-edge tables with valve, choke, sensor and actuator locations, GIS pipeline routes, time-series operating data; STEP, IGES and DXF only where physical equipment geometry is needed', accepts: ['network', 'table', 'polylines', 'params'] },
  integ: { classes: 'STEP, IGES, STL, OBJ, DXF, IFC; structural / FEA meshes (Nastran BDF, Abaqus INP, ANSYS CDB, LS-DYNA KEY, UNV, Gmsh, VTK); inspection data: wall-thickness maps, corrosion and pit-depth grids, dent and free-span geometry, laser-scan point clouds (CSV, XYZ, VTK, PLY, LAS); GIS and bathymetric formats', accepts: ['mesh', 'table', 'grid', 'points', 'polylines', 'network', 'params'] },
  econ: { classes: 'No mandatory CAD geometry. Network and equipment identifiers and engineering results as CSV, TSV, TXT/DAT, JSON, XML, YAML, NetCDF, HDF5 and other tabular or scientific tables; optional GIS or network references where costs or risks depend on location (route length, water depth, well depth, equipment counts)', accepts: ['network', 'polylines', 'table', 'params'] },
});

const extName = (name) => { const s = String(name).toLowerCase(), m = s.match(/\.([a-z0-9_]+)$/); return /\.nii\.gz$/.test(s) ? 'nii.gz' : m ? m[1] : ''; };
const baseName = (name) => String(name).replace(/^.*[\\/]/, '').toLowerCase();
const FOAM_FILES = ['points', 'faces', 'owner', 'neighbour'];

/** File name -> FORMATS entry (the first one registered for the extension) or null. */
export function formatOf(filename) {
  const ext = extName(filename);
  if (!ext) return FOAM_FILES.includes(baseName(filename).replace(/\.gz$/, '')) ? FORMATS.find((f) => f.ext[0] === 'foam') : null;
  return FORMATS.find((f) => f.ext.includes(ext)) || null;
}
const fmtNamed = (name) => FORMATS.find((f) => f.name === name);

// ---- Small shared helpers ------------------------------------------------------------------------
function fail(msg) { const e = new Error(msg); e.user = true; throw e; }
const latin1 = new TextDecoder('latin1'), utf8 = new TextDecoder();
const splitLines = (t) => t.split(/\r\n|\n|\r/);
const numsOf = (s) => { const m = s.match(/[^\s,;]+/g) || [], out = new Float64Array(m.length); for (let i = 0; i < m.length; i++) out[i] = +m[i]; return out; };
const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const vlen = (a) => Math.hypot(a[0], a[1], a[2]);
const unit = (a) => { const l = vlen(a); return l > 0 ? [a[0] / l, a[1] / l, a[2] / l] : [0, 0, 0]; };
const TAU = 2 * Math.PI, ARC_N = 48;

function bboxOf(a, stride = 3, dim = stride) {
  const lo = new Array(dim).fill(Infinity), hi = new Array(dim).fill(-Infinity);
  for (let i = 0; i + dim <= a.length; i += stride) for (let c = 0; c < dim; c++) { const v = a[i + c]; if (v < lo[c]) lo[c] = v; if (v > hi[c]) hi[c] = v; }
  return { min: lo, max: hi };
}
function meshGeom(tri, extra = {}) {
  if (tri.length > L.triangles * 9) fail(`Mesh has too many triangles (limit ${L.triangles}).`);
  const t = Array.isArray(tri) ? tri : Array.from(tri);
  t.length -= t.length % 9;
  if (!t.length) fail('No triangles could be read from the file.');
  for (let i = 0; i < t.length; i++) if (!Number.isFinite(t[i])) fail('The mesh contains non-numeric coordinates.');
  return { kind: 'mesh', triangles: t, count: t.length / 9, bbox: bboxOf(t), ...extra };
}
function polyGeom(polys, extra = {}) {
  const good = polys.filter((p) => p.x.length > 0 && p.x.length === p.y.length && p.x.every(Number.isFinite) && p.y.every(Number.isFinite));
  if (!good.length) fail('No line or outline geometry could be read from the file.');
  const lo = [Infinity, Infinity], hi = [-Infinity, -Infinity];
  let zl = Infinity, zh = -Infinity;                     // optional p.z makes a polyline 3-D (routes, risers, well paths)
  for (const p of good) {
    if (p.z !== undefined && !(p.z && p.z.length === p.x.length && p.z.every(Number.isFinite))) delete p.z;
    for (let i = 0; i < p.x.length; i++) { lo[0] = Math.min(lo[0], p.x[i]); hi[0] = Math.max(hi[0], p.x[i]); lo[1] = Math.min(lo[1], p.y[i]); hi[1] = Math.max(hi[1], p.y[i]); if (p.z) { zl = Math.min(zl, p.z[i]); zh = Math.max(zh, p.z[i]); } }
  }
  if (zl <= zh) { lo.push(zl); hi.push(zh); }
  return { kind: 'polylines', polylines: good, bbox: { min: lo, max: hi }, ...extra };
}
/** xyz: flat array of triples. Rows with a non-finite coordinate are dropped; large clouds are sub-sampled evenly. */
function pointsGeom(xyz, extra = {}) {
  const warnings = extra.warnings || (extra.warnings = []);
  let n = Math.floor(xyz.length / 3), bad = 0;
  const step = Math.max(1, Math.ceil(n / L.points)), out = new Float64Array(Math.ceil(n / step) * 3);
  let m = 0;
  for (let i = 0; i < n; i += step) { const x = xyz[3 * i], y = xyz[3 * i + 1], z = xyz[3 * i + 2]; if (Number.isFinite(x) && Number.isFinite(y) && Number.isFinite(z)) { out[m++] = x; out[m++] = y; out[m++] = z; } else bad++; }
  if (!m) fail('No points could be read from the file.');
  if (step > 1) warnings.push(`Point cloud of ${n} points evenly sub-sampled to ${m / 3} (every ${step}th point).`);
  if (bad) warnings.push(`${bad} rows with non-numeric coordinates were dropped.`);
  const points = out.subarray(0, m);
  return { kind: 'points', points, count: m / 3, bbox: bboxOf(points), ...extra, stats: { points: m / 3, sourcePoints: n, ...(extra.stats || {}) } };
}
/** Raster -> grid. zf is row-major with row j at y = yOf(j); rows/columns are re-ordered so x and y ascend. */
function gridGeom(zf, nx, ny, xOf, yOf, { nodata, geographic = false, warnings = [], stats = {} } = {}) {
  if (!(nx > 0 && ny > 0) || zf.length < nx * ny) fail('Raster dimensions do not match the amount of data in the file.');
  const s = Math.max(1, Math.ceil(Math.sqrt((nx * ny) / L.grid))), ix = [], iy = [];
  for (let i = 0; i < nx; i += s) ix.push(i);
  for (let j = 0; j < ny; j += s) iy.push(j);
  if (s > 1) warnings.push(`Raster of ${nx} × ${ny} cells down-sampled by ${s} in each direction.`);
  if (ix.length > 1 && xOf(ix[0]) > xOf(ix[1])) ix.reverse();
  if (iy.length > 1 && yOf(iy[0]) > yOf(iy[1])) iy.reverse();
  const isNo = typeof nodata === 'function' ? nodata : (v) => v === nodata;
  let nod = 0, lo = Infinity, hi = -Infinity;
  const z = iy.map((j) => ix.map((i) => { const v = zf[j * nx + i]; if (!Number.isFinite(v) || isNo(v)) { nod++; return NaN; } if (v < lo) lo = v; if (v > hi) hi = v; return v; }));
  if (nod === ix.length * iy.length) fail('The raster holds no valid values.');
  const x = ix.map(xOf), y = iy.map(yOf);
  return { kind: 'grid', grid: { x, y, z, nodata: nod, geographic }, bbox: { min: [x[0], y[0], lo], max: [x[x.length - 1], y[y.length - 1], hi] }, warnings, stats: { nx: x.length, ny: y.length, zmin: lo, zmax: hi, ...stats } };
}
function voxGeom(nx, ny, nz, data, spacing = [1, 1, 1], extra = {}, origin = [0, 0, 0]) {
  const sp = [0, 1, 2].map((k) => (Number.isFinite(+spacing[k]) && +spacing[k] > 0 ? +spacing[k] : 1));
  let solid = 0;
  for (let i = 0; i < data.length; i++) solid += data[i];
  return { kind: 'voxels', voxels: { nx, ny, nz, data, spacing: sp, origin }, bbox: { min: origin.slice(), max: [origin[0] + nx * sp[0], origin[1] + ny * sp[1], origin[2] + nz * sp[2]] }, ...extra, stats: { nx, ny, nz, porosity: 1 - solid / data.length, ...(extra.stats || {}) } };
}

function otsu(v, n, lo, hi) {
  const nb = 256, h = new Float64Array(nb), sc = (nb - 1) / (hi - lo || 1);
  for (let i = 0; i < n; i++) { const x = v[i]; if (x === x) h[Math.round((x - lo) * sc)]++; }
  let tot = 0, sumAll = 0;
  for (let b = 0; b < nb; b++) { tot += h[b]; sumAll += b * h[b]; }
  let w0 = 0, s0 = 0, best = -1, bi = 0, bj = 0;
  for (let b = 0; b < nb - 1; b++) {
    w0 += h[b]; s0 += b * h[b];
    const w1 = tot - w0;
    if (!w0 || !w1) continue;
    const d = s0 / w0 - (sumAll - s0) / w1, v2 = w0 * w1 * d * d;
    if (v2 > best * (1 + 1e-12)) { best = v2; bi = bj = b; } else if (Math.abs(v2 - best) <= best * 1e-12) bj = b;
  }
  return lo + (0.5 * (bi + bj) + 0.5) / sc;
}
/** Scalar volume (x fastest, then y, then z) -> solid/pore voxels. Labels: nonzero = solid; grey values: Otsu or opts.threshold. */
function voxelsFromValues(v, nx, ny, nz, spacing, opts = {}, extra = {}) {
  const n = nx * ny * nz;
  if (![nx, ny, nz].every((k) => Number.isInteger(k) && k > 0)) fail('Volume dimensions are missing or invalid.');
  if (n > L.voxels) fail(`Volume has too many voxels (limit ${L.voxels / 1e6} million).`);
  if (v.length < n) fail('Volume data is shorter than its dimensions require (file truncated?).');
  let lo = Infinity, hi = -Infinity, isInt = true;
  const seen = new Set();
  for (let i = 0; i < n; i++) { const x = +v[i]; if (x < lo) lo = x; if (x > hi) hi = x; if (isInt) { if (x !== Math.floor(x)) isInt = false; else if (seen.size <= 16) seen.add(x); } }
  if (!(hi >= lo)) fail('The volume holds no numeric values.');
  const warnings = extra.warnings || [], data = new Uint8Array(n), inv = !!opts.invert;
  let thr, method;
  if (opts.threshold !== undefined && opts.threshold !== null && Number.isFinite(+opts.threshold)) { thr = +opts.threshold; method = 'threshold'; }
  else if (isInt && seen.size === 2) { thr = 0.5 * (lo + hi); method = 'binary'; }
  else if (isInt && seen.size <= 16) { thr = 0; method = 'labels'; }
  else { thr = otsu(v, n, lo, hi); method = 'otsu'; }
  if (method === 'labels') for (let i = 0; i < n; i++) data[i] = (v[i] != 0) !== inv ? 1 : 0;
  else for (let i = 0; i < n; i++) data[i] = (v[i] > thr) !== inv ? 1 : 0;
  warnings.push(method === 'labels' ? `Labelled volume: ${inv ? 'zero' : 'non-zero'} voxels taken as solid.` : method === 'binary' ? `Two-valued volume: voxels equal to ${inv ? lo : hi} taken as solid.` : `Grey values ${method === 'otsu' ? 'split by Otsu threshold' : 'split at the given threshold'} ${+thr.toPrecision(6)}: ${inv ? 'darker' : 'brighter'} voxels taken as solid.`);
  return voxGeom(nx, ny, nz, data, spacing, { ...extra, warnings, stats: { threshold: thr, method, min: lo, max: hi, ...(extra.stats || {}) } });
}

// ---- Binary helpers ------------------------------------------------------------------------------
const DT = { u1: [Uint8Array, 1], i1: [Int8Array, 1], u2: [Uint16Array, 2], i2: [Int16Array, 2], u4: [Uint32Array, 4], i4: [Int32Array, 4], f4: [Float32Array, 4], f8: [Float64Array, 8], i8: [BigInt64Array, 8], u8: [BigUint64Array, 8] };
const DT_ALIAS = { bool: 'u1', b1: 'u1', uint8: 'u1', uchar: 'u1', 'unsigned char': 'u1', uint8_t: 'u1', byte: 'u1', ubyte: 'u1', met_uchar: 'u1', int8: 'i1', char: 'i1', 'signed char': 'i1', int8_t: 'i1', met_char: 'i1', uint16: 'u2', ushort: 'u2', 'unsigned short': 'u2', 'unsigned short int': 'u2', uint16_t: 'u2', met_ushort: 'u2', int16: 'i2', short: 'i2', 'short int': 'i2', 'signed short': 'i2', int16_t: 'i2', met_short: 'i2', uint32: 'u4', uint: 'u4', 'unsigned int': 'u4', uint32_t: 'u4', met_uint: 'u4', int32: 'i4', int: 'i4', 'signed int': 'i4', int32_t: 'i4', met_int: 'i4', float32: 'f4', float: 'f4', single: 'f4', met_float: 'f4', float64: 'f8', double: 'f8', met_double: 'f8', int64: 'i8', longlong: 'i8', 'long long': 'i8', int64_t: 'i8', met_long_long: 'i8', vtktypeint64: 'i8', uint64: 'u8', ulonglong: 'u8', 'unsigned long long': 'u8', uint64_t: 'u8', met_ulong_long: 'u8', vtktypeuint64: 'u8', vtktypeint32: 'i4', vtkidtype: 'i4', unsigned_char: 'u1', unsigned_short: 'u2', unsigned_int: 'u4', long: 'i8', unsigned_long: 'u8', bit: 'u1' };
const dtypeOf = (s) => { const k = String(s ?? '').trim().toLowerCase(); return DT[k] ? k : DT_ALIAS[k] || null; };
const HOST_LE = new Uint8Array(new Uint16Array([1]).buffer)[0] === 1;
/** Bounds-checked copy of `count` values of type dt starting at byte `off`; 64-bit integers come back as doubles. */
function typed(u8, off, count, dt, little = true) {
  const d = DT[dt];
  if (!d) fail(`Unsupported data type "${dt}".`);
  const nb = count * d[1];
  if (!Number.isInteger(count) || count < 0 || !Number.isSafeInteger(nb) || !(off >= 0) || off + nb > u8.length) fail('The file is truncated: its data section is shorter than the header declares.');
  const b = u8.slice(off, off + nb);
  if (d[1] > 1 && little !== HOST_LE) for (let i = 0; i < nb; i += d[1]) for (let a = 0, z = d[1] - 1; a < z; a++, z--) { const t = b[i + a]; b[i + a] = b[i + z]; b[i + z] = t; }
  const arr = new d[0](b.buffer, 0, count);
  if (dt === 'i8' || dt === 'u8') { const f = new Float64Array(count); for (let i = 0; i < count; i++) f[i] = Number(arr[i]); return f; }
  return arr;
}
function b64(s) {
  let bin;
  try { bin = atob(s.replace(/\s+/g, '')); } catch { fail('Invalid base64 data in the file.'); }
  const u = new Uint8Array(bin.length);
  for (let i = 0; i < bin.length; i++) u[i] = bin.charCodeAt(i);
  return u;
}
/** Size-capped decompression: 'deflate' (zlib), 'deflate-raw' or 'gzip'. */
async function inflate(bytes, format = 'deflate', max = L.inflated) {
  if (typeof DecompressionStream !== 'function') fail('This browser cannot decompress data (DecompressionStream is missing).');
  const chunks = [];
  let n = 0, over = false;
  try {
    const reader = new Blob([bytes]).stream().pipeThrough(new DecompressionStream(format)).getReader();
    for (;;) {
      const { done, value } = await reader.read();
      if (done) break;
      n += value.length;
      if (n > max) { over = true; await reader.cancel().catch(() => {}); break; }
      chunks.push(value);
    }
  } catch { fail('Compressed data in the file is corrupt or truncated.'); }
  if (over) fail('Decompressed data is too large.');
  const out = new Uint8Array(n);
  let p = 0;
  for (const c of chunks) { out.set(c, p); p += c.length; }
  return out;
}
/** Zip central-directory reader -> Map(name -> bytes) for the members accepted by want(name). */
async function unzipFiles(buf, want = () => true) {
  const u8 = new Uint8Array(buf), dv = new DataView(buf), out = new Map();
  let eocd = -1;
  for (let i = u8.length - 22; i >= Math.max(0, u8.length - 66000); i--) if (dv.getUint32(i, true) === 0x06054b50) { eocd = i; break; }
  if (eocd < 0) fail('Not a valid zip archive.');
  const n = dv.getUint16(eocd + 10, true);
  let p = dv.getUint32(eocd + 16, true);
  for (let k = 0; k < n && p + 46 <= u8.length; k++) {
    if (dv.getUint32(p, true) !== 0x02014b50) break;
    const method = dv.getUint16(p + 10, true), csize = dv.getUint32(p + 20, true), usize = dv.getUint32(p + 24, true);
    const nl = dv.getUint16(p + 28, true), el = dv.getUint16(p + 30, true), cl = dv.getUint16(p + 32, true), lho = dv.getUint32(p + 42, true);
    const name = utf8.decode(u8.subarray(p + 46, p + 46 + nl));
    p += 46 + nl + el + cl;
    if (name.endsWith('/') || !want(name)) continue;
    if (csize === 0xffffffff || usize === 0xffffffff || lho === 0xffffffff) fail('Zip64 archives are not supported.');
    if (usize > L.inflated) fail('A file inside the archive is too large.');
    if (lho + 30 > u8.length) fail('The zip archive is truncated.');
    const start = lho + 30 + dv.getUint16(lho + 26, true) + dv.getUint16(lho + 28, true);
    if (start + csize > u8.length) fail('The zip archive is truncated.');
    const data = u8.subarray(start, start + csize);
    if (method !== 0 && method !== 8) fail('The zip archive uses an unsupported compression method.');
    out.set(name, method === 0 ? data : await inflate(data, 'deflate-raw'));
  }
  return out;
}
const isZip = (u8) => u8.length > 4 && u8[0] === 0x50 && u8[1] === 0x4b && (u8[2] === 3 || u8[2] === 5);
const isGzip = (u8) => u8.length > 2 && u8[0] === 0x1f && u8[1] === 0x8b;
const SKIP_KEYS = new Set(['__proto__', 'constructor', 'prototype']);
function safeJSON(text) {
  try { return JSON.parse(text.replace(/^﻿/, ''), (k, v) => (SKIP_KEYS.has(k) ? undefined : v)); } catch { fail('The file is not valid JSON.'); }
}

// ---- Minimal XML tree builder (no DTD processing, no custom entities) ----------------------------------
const xLocal = (s) => { const k = s.lastIndexOf(':'); return (k >= 0 ? s.slice(k + 1) : s).toLowerCase(); };
const XENT = { lt: '<', gt: '>', amp: '&', quot: '"', apos: "'" };
const unent = (s) => (s.indexOf('&') < 0 ? s : s.replace(/&(#x[0-9a-f]{1,6}|#\d{1,7}|\w{2,5});/gi, (m, c) => {
  if (c[0] !== '#') return XENT[c] ?? m;
  const cp = c[1] === 'x' || c[1] === 'X' ? parseInt(c.slice(2), 16) : parseInt(c.slice(1), 10);
  return cp > 0 && cp <= 0x10ffff && !(cp >= 0xd800 && cp <= 0xdfff) ? String.fromCodePoint(cp) : '';
}));
/** Returns the document node; element and attribute names are lower-cased local names. */
function parseXML(text) {
  const mk = (name) => ({ name, attrs: Object.create(null), children: [], text: '' });
  const root = mk('#document'), stack = [root], n = text.length;
  let i = 0, count = 0;
  while (i < n) {
    const lt = text.indexOf('<', i), top = stack[stack.length - 1];
    if (lt < 0) break;
    if (lt > i && stack.length > 1) { const s = text.slice(i, lt); if (/\S/.test(s)) top.text += unent(s); }
    const c1 = text[lt + 1];
    if (c1 === '!' || c1 === '?') {
      let e;
      if (text.startsWith('<!--', lt)) e = text.indexOf('-->', lt + 4) + 2;
      else if (text.startsWith('<![CDATA[', lt)) { e = text.indexOf(']]>', lt + 9); if (e >= 0) { top.text += text.slice(lt + 9, e); e += 2; } }
      else if (c1 === '?') e = text.indexOf('?>', lt + 2) + 1;
      else { e = text.indexOf('>', lt + 2); const b = text.indexOf('[', lt + 2); if (b >= 0 && b < e) { e = text.indexOf(']>', b); if (e >= 0) e += 1; } }
      if (e < lt + 2) break;
      i = e + 1; continue;
    }
    if (c1 === '/') {
      const e = text.indexOf('>', lt + 2);
      if (e < 0) break;
      const nm = xLocal(text.slice(lt + 2, e).trim());
      for (let k = stack.length - 1; k > 0; k--) if (stack[k].name === nm) { stack.length = k; break; }
      i = e + 1; continue;
    }
    let e = lt + 1, q = 0;
    for (; e < n; e++) { const c = text.charCodeAt(e); if (q) { if (c === q) q = 0; } else if (c === 34 || c === 39) q = c; else if (c === 62) break; }
    if (e >= n) fail('The XML document is truncated.');
    const selfClose = text.charCodeAt(e - 1) === 47, tag = text.slice(lt + 1, selfClose ? e - 1 : e), m = /^[^\s]+/.exec(tag);
    i = e + 1;
    if (!m) continue;
    const node = mk(xLocal(m[0])), re = /([^\s=]+)\s*=\s*(?:"([^"]*)"|'([^']*)')/g;
    let a;
    while ((a = re.exec(tag))) { const key = xLocal(a[1]); if (!SKIP_KEYS.has(key)) node.attrs[key] = unent(a[2] ?? a[3]); }
    top.children.push(node);
    if (++count > L.xmlNodes) fail('The XML document has too many elements.');
    if (!selfClose) { if (stack.length > 300) fail('The XML document is nested too deeply.'); stack.push(node); }
  }
  if (!root.children.length) fail('No XML elements were found in the file.');
  return root;
}
/** All descendants with one of the given names, in document order. */
function xAll(node, ...names) {
  const out = [], st = [node];
  while (st.length) { const nd = st.pop(); if (nd !== node && names.includes(nd.name)) out.push(nd); for (let k = nd.children.length - 1; k >= 0; k--) st.push(nd.children[k]); }
  return out;
}
const xKids = (node, name) => node.children.filter((c) => c.name === name);
const xKid = (node, name) => node.children.find((c) => c.name === name);

// ---- Cell meshes -> boundary surface / outline --------------------------------------------------------
const CELL_FACES = {
  tet: [[0, 2, 1], [0, 1, 3], [1, 2, 3], [0, 3, 2]],
  hex: [[0, 3, 2, 1], [4, 5, 6, 7], [0, 1, 5, 4], [1, 2, 6, 5], [2, 3, 7, 6], [3, 0, 4, 7]],
  wedge: [[0, 2, 1], [3, 4, 5], [0, 1, 4, 3], [1, 2, 5, 4], [2, 0, 3, 5]],
  pyramid: [[0, 3, 2, 1], [0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]],
};
function cellStore() {
  const t = [], v = [];
  return { t, v, add(type, verts) { if (t.length >= L.cells) fail(`Mesh has too many cells (limit ${L.cells / 1e6} million).`); t.push(type); v.push(verts); } };
}
/** Chain undirected index pairs into open or closed vertex chains. */
function chainEdges(edges) {
  const adj = new Map(), used = new Uint8Array(edges.length), out = [];
  const put = (a, k) => { const l = adj.get(a); if (l) l.push(k); else adj.set(a, [k]); };
  edges.forEach(([a, b], k) => { put(a, k); put(b, k); });
  const walk = (start) => {
    const ids = [start];
    for (let cur = start; ;) {
      const k = adj.get(cur).find((q) => !used[q]);
      if (k === undefined) return { ids, closed: false };
      used[k] = 1;
      cur = edges[k][0] === cur ? edges[k][1] : edges[k][0];
      if (cur === start) return { ids, closed: true };
      ids.push(cur);
    }
  };
  for (const [node, ks] of adj) if (ks.length % 2 === 1 && ks.some((k) => !used[k])) out.push(walk(node));
  for (let k = 0; k < edges.length; k++) if (!used[k]) out.push(walk(edges[k][0]));
  return out.filter((c) => c.ids.length > 1);
}
/**
 * Nodes (flat xyz) + cells -> Geometry. Volume cells give their boundary surface (faces owned by one cell, oriented
 * outwards); a flat tri/quad mesh gives its boundary outline as polylines; a curved tri/quad mesh is a surface.
 */
function cellMesh(xyz, store, extra = {}) {
  const nn = Math.floor(xyz.length / 3), cellTypes = {}, P = (i) => [xyz[3 * i], xyz[3 * i + 1], xyz[3 * i + 2]];
  const warnings = extra.warnings || [];
  let has3 = false, has2 = false;
  for (let c = 0; c < store.t.length; c++) {
    const t = store.t[c];
    cellTypes[t] = (cellTypes[t] || 0) + 1;
    for (const i of store.v[c]) if (!Number.isInteger(i) || i < 0 || i >= nn) fail('The mesh connectivity refers to nodes that do not exist.');
    if (CELL_FACES[t]) has3 = true; else if (t !== 'line') has2 = true;
  }
  const stats = { nodes: nn, cells: store.t.length, cellTypes, ...(extra.stats || {}) };
  if (!nn) fail('No mesh nodes were found in the file.');
  if (!has3 && !has2) {
    if (cellTypes.line) {
      const ch = chainEdges(store.v.filter((_, c) => store.t[c] === 'line'));
      let z3 = false;
      for (let i = 0; i < nn && !z3; i++) if (xyz[3 * i + 2] !== xyz[2]) z3 = true;       // out-of-plane line meshes (pipes, beams) keep z
      return polyGeom(ch.map((c) => ({ x: c.ids.map((i) => xyz[3 * i]), y: c.ids.map((i) => xyz[3 * i + 1]), ...(z3 ? { z: c.ids.map((i) => xyz[3 * i + 2]) } : {}), closed: c.closed })), { ...extra, warnings, stats: { ...stats, dimension: 1 } });
    }
    warnings.push('The file holds nodes but no surface or volume cells; the nodes are returned as points.');
    return pointsGeom(xyz, { ...extra, warnings, stats });
  }
  const tri = [], fan = (ids) => { for (let k = 1; k + 1 < ids.length; k++) for (const q of [ids[0], ids[k], ids[k + 1]]) tri.push(xyz[3 * q], xyz[3 * q + 1], xyz[3 * q + 2]); };
  if (has3) {
    const faces = new Map(), wide = nn > 2e5;    // key = the three smallest node ids of a face; a face met twice is interior
    for (let c = 0; c < store.t.length; c++) {
      const cf = CELL_FACES[store.t[c]], v = store.v[c];
      if (!cf) continue;
      for (let f = 0; f < cf.length; f++) {
        let a = Infinity, b = Infinity, d = Infinity;
        for (const k of cf[f]) { const x = v[k]; if (x < a) { d = b; b = a; a = x; } else if (x < b) { d = b; b = x; } else if (x < d) d = x; }
        const key = wide ? a + ',' + b + ',' + d : (a * nn + b) * nn + d;
        if (!faces.delete(key)) faces.set(key, c * 8 + f);
      }
    }
    for (const code of faces.values()) {
      const c = Math.floor(code / 8), v = store.v[c], ids = CELL_FACES[store.t[c]][code % 8].map((k) => v[k]), cc = [0, 0, 0], fc = [0, 0, 0], nrm = [0, 0, 0], pts = ids.map(P);
      for (const i of v) { cc[0] += xyz[3 * i] / v.length; cc[1] += xyz[3 * i + 1] / v.length; cc[2] += xyz[3 * i + 2] / v.length; }
      for (let k = 0; k < pts.length; k++) { const a = pts[k], b = pts[(k + 1) % pts.length]; for (let d = 0; d < 3; d++) fc[d] += a[d] / pts.length; nrm[0] += (a[1] - b[1]) * (a[2] + b[2]); nrm[1] += (a[2] - b[2]) * (a[0] + b[0]); nrm[2] += (a[0] - b[0]) * (a[1] + b[1]); }
      fan(dot(nrm, sub(fc, cc)) < 0 ? ids.reverse() : ids);
      if (tri.length > L.triangles * 9) fail(`The mesh boundary has too many triangles (limit ${L.triangles}).`);
    }
    if (!tri.length) fail('The volume mesh has no boundary faces.');
    return meshGeom(tri, { ...extra, warnings, stats: { ...stats, dimension: 3, boundaryTriangles: tri.length / 9 } });
  }
  const bb = bboxOf(xyz), ext = [0, 1, 2].map((k) => bb.max[k] - bb.min[k]), big = Math.max(...ext), flat = ext.findIndex((e) => e <= 1e-9 * big);
  if (extra.flat === undefined ? flat >= 0 : extra.flat) {
    const drop = flat >= 0 ? flat : 2, ax = [0, 1, 2].filter((k) => k !== drop), em = new Map();
    for (let c = 0; c < store.t.length; c++) {
      const v = store.v[c];
      if (store.t[c] === 'line') continue;
      for (let k = 0; k < v.length; k++) { const a = v[k], b = v[(k + 1) % v.length], key = a < b ? a + ',' + b : b + ',' + a, e = em.get(key); if (e) e.n++; else em.set(key, { a, b, n: 1 }); }
    }
    const be = [];
    for (const e of em.values()) if (e.n === 1) be.push([e.a, e.b]);
    const polys = chainEdges(be).map((c) => ({ x: c.ids.map((i) => xyz[3 * i + ax[0]]), y: c.ids.map((i) => xyz[3 * i + ax[1]]), closed: c.closed }));
    if (drop !== 2) warnings.push(`Planar mesh lies in a plane of constant ${'xyz'[drop]}; outline given in the ${'xyz'[ax[0]]}-${'xyz'[ax[1]]} plane.`);
    const { flat: _f, ...rest } = extra;
    return polyGeom(polys, { ...rest, warnings, stats: { ...stats, dimension: 2, boundaryEdges: be.length } });
  }
  for (let c = 0; c < store.t.length; c++) if (store.t[c] !== 'line') { fan(store.v[c]); if (tri.length > L.triangles * 9) fail(`Mesh has too many triangles (limit ${L.triangles}).`); }
  const { flat: _f, ...rest } = extra;
  return meshGeom(tri, { ...rest, warnings, stats: { ...stats, dimension: 3, surfaceOnly: true } });
}
/** Structured blocks [{ xyz, ni, nj, nk }] -> boundary surface (3-D blocks) or outline / surface (2-D blocks). */
function structuredGeom(blocks, extra = {}) {
  const warnings = extra.warnings || [], tri = [];
  let nodes = 0, cells = 0, any3 = false;
  for (const b of blocks) {
    const { xyz, ni, nj, nk } = b;
    if (![ni, nj, nk].every((k) => Number.isInteger(k) && k > 0) || xyz.length < 3 * ni * nj * nk) fail('Structured grid dimensions do not match its coordinates.');
    nodes += ni * nj * nk; cells += Math.max(1, ni - 1) * Math.max(1, nj - 1) * Math.max(1, nk - 1);
    if (ni > 1 && nj > 1 && nk > 1) any3 = true;
  }
  if (!any3) {
    const store = cellStore(), all = new Float64Array(nodes * 3);
    let base = 0;
    for (const { xyz, ni, nj, nk } of blocks) {
      all.set(xyz.subarray(0, 3 * ni * nj * nk), base * 3);
      const d = [ni, nj, nk], ax = [0, 1, 2].filter((k) => d[k] > 1), st = [1, ni, ni * nj];
      if (ax.length === 2) for (let b = 0; b < d[ax[1]] - 1; b++) for (let a = 0; a < d[ax[0]] - 1; a++) { const o = base + a * st[ax[0]] + b * st[ax[1]]; store.add('quad', [o, o + st[ax[0]], o + st[ax[0]] + st[ax[1]], o + st[ax[1]]]); }
      else if (ax.length === 1) for (let a = 0; a < d[ax[0]] - 1; a++) store.add('line', [base + a * st[ax[0]], base + (a + 1) * st[ax[0]]]);
      base += ni * nj * nk;
    }
    return cellMesh(all, store, { ...extra, warnings, stats: { structured: true, blocks: blocks.length, ...(extra.stats || {}) } });
  }
  for (const { xyz, ni, nj, nk } of blocks) {
    if (!(ni > 1 && nj > 1 && nk > 1)) continue;
    if (2 * ((ni - 1) * (nj - 1) + (nj - 1) * (nk - 1) + (ni - 1) * (nk - 1)) * 2 + tri.length / 9 > L.triangles) fail(`The structured grid boundary has too many triangles (limit ${L.triangles}).`);
    const id = (i, j, k) => 3 * (i + ni * (j + nj * k)), p = (o) => [xyz[o], xyz[o + 1], xyz[o + 2]], o0 = p(0);
    const rh = dot(cross(sub(p(id(1, 0, 0)), o0), sub(p(id(0, 1, 0)), o0)), sub(p(id(0, 0, 1)), o0)) >= 0;
    const quad = (a, b, c, d) => { const q = rh ? [a, b, c, a, c, d] : [a, c, b, a, d, c]; for (const o of q) tri.push(xyz[o], xyz[o + 1], xyz[o + 2]); };
    for (let k = 0; k < nk - 1; k++) for (let j = 0; j < nj - 1; j++) { quad(id(0, j, k), id(0, j, k + 1), id(0, j + 1, k + 1), id(0, j + 1, k)); quad(id(ni - 1, j, k), id(ni - 1, j + 1, k), id(ni - 1, j + 1, k + 1), id(ni - 1, j, k + 1)); }
    for (let k = 0; k < nk - 1; k++) for (let i = 0; i < ni - 1; i++) { quad(id(i, 0, k), id(i + 1, 0, k), id(i + 1, 0, k + 1), id(i, 0, k + 1)); quad(id(i, nj - 1, k), id(i, nj - 1, k + 1), id(i + 1, nj - 1, k + 1), id(i + 1, nj - 1, k)); }
    for (let j = 0; j < nj - 1; j++) for (let i = 0; i < ni - 1; i++) { quad(id(i, j, 0), id(i, j + 1, 0), id(i + 1, j + 1, 0), id(i + 1, j, 0)); quad(id(i, j, nk - 1), id(i + 1, j, nk - 1), id(i + 1, j + 1, nk - 1), id(i, j + 1, nk - 1)); }
  }
  if (blocks.length > 1) warnings.push(`${blocks.length} structured blocks: every block boundary is kept, including block-to-block interfaces.`);
  return meshGeom(tri, { ...extra, warnings, stats: { nodes, cells, cellTypes: { hex: cells }, dimension: 3, structured: true, blocks: blocks.length, ...(extra.stats || {}) } });
}

// ---- Curves and planar triangulation ---------------------------------------------------------------------
/** Rational B-spline point by de Boor's algorithm (w may be null). */
function nurbsPoint(deg, knots, cps, w, u) {
  const n = cps.length - 1;
  let k = deg;
  if (u >= knots[n + 1]) k = n; else while (k < n && u >= knots[k + 1]) k++;
  while (k > deg && knots[k] === knots[k + 1]) k--;
  const d = [];
  for (let j = 0; j <= deg; j++) { const i = j + k - deg, wi = w ? w[i] : 1, c = cps[i]; d.push([c[0] * wi, c[1] * wi, c[2] * wi, wi]); }
  for (let r = 1; r <= deg; r++) for (let j = deg; j >= r; j--) {
    const i = j + k - deg, den = knots[i + deg - r + 1] - knots[i], a = den > 0 ? (u - knots[i]) / den : 0;
    for (let c = 0; c < 4; c++) d[j][c] = (1 - a) * d[j - 1][c] + a * d[j][c];
  }
  const h = d[deg][3] || 1;
  return [d[deg][0] / h, d[deg][1] / h, d[deg][2] / h];
}
/** Sampled B-spline curve, or null when the definition is inconsistent. */
function nurbsCurve(deg, knots, cps, w, u0, u1, nSeg) {
  if (!Number.isInteger(deg) || deg < 1 || deg > 25 || cps.length <= deg || knots.length !== cps.length + deg + 1) return null;
  if (cps.some((c) => !c || !c.every(Number.isFinite)) || knots.some((k) => !Number.isFinite(k)) || (w && (w.length !== cps.length || w.some((x) => !(x > 0))))) return null;
  for (let i = 1; i < knots.length; i++) if (knots[i] < knots[i - 1]) return null;
  const a = Math.max(knots[deg], Number.isFinite(u0) ? u0 : -Infinity), b = Math.min(knots[cps.length], Number.isFinite(u1) ? u1 : Infinity);
  if (!(b > a)) return null;
  const n = Math.min(2000, Math.max(8, nSeg | 0)), out = [];
  for (let i = 0; i <= n; i++) out.push(nurbsPoint(deg, knots, cps, w, i === n ? b : a + ((b - a) * i) / n));
  return out;
}
/** Ear clipping of a CCW outer ring with CW holes; points are [u, v, id]; returns id triples. */
function triangulate(outer, holes) {
  const X = outer.map((p) => p[0]), Y = outer.map((p) => p[1]), ID = outer.map((p) => p[2]);
  const insert = (arr, at, items) => { for (let k = 0; k < items.length; k += 8192) arr.splice(at + k, 0, ...items.slice(k, k + 8192)); };
  // bridge every hole (rightmost first) to the vertex of the current outline that its rightmost point can see
  const hs = holes.filter((h) => h.length >= 3).map((h) => { let m = 0; for (let k = 1; k < h.length; k++) if (h[k][0] > h[m][0]) m = k; return { h, m }; }).sort((a, b) => b.h[b.m][0] - a.h[a.m][0]);
  for (const { h, m } of hs) {
    const mx = h[m][0], my = h[m][1], np = X.length;
    let best = Infinity, bi = -1, pick = 0;
    for (let i = 0, j = np - 1; i < np; j = i++) {              // edge j -> i against the ray from M towards +x
      const ay = Y[j], by = Y[i];
      if ((ay > my) === (by > my)) continue;
      const x = X[j] + ((my - ay) * (X[i] - X[j])) / (by - ay);
      if (x >= mx && x < best) { best = x; bi = j; }
    }
    if (bi < 0) { let d = Infinity; for (let k = 0; k < np; k++) { const q = (X[k] - mx) ** 2 + (Y[k] - my) ** 2; if (q < d) { d = q; pick = k; } } }
    else {
      const bj = (bi + 1) % np;
      pick = X[bi] > X[bj] ? bi : bj;
      const cx = X[pick], cy = Y[pick], xmax = Math.max(best, cx);
      let tanMin = Infinity;
      for (let k = 0; k < np; k++) {                              // a vertex inside triangle (M, hit point, C) hides C: take the one closest to the ray
        const qx = X[k], qy = Y[k];
        if (qx < mx || qx > xmax || (qx === cx && qy === cy)) continue;
        const d1 = (best - mx) * (qy - my), d2 = (cx - best) * (qy - my) - (cy - my) * (qx - best), d3 = (mx - cx) * (qy - cy) - (my - cy) * (qx - cx);
        if (!((d1 >= 0 && d2 >= 0 && d3 >= 0) || (d1 <= 0 && d2 <= 0 && d3 <= 0))) continue;
        const t = Math.abs(qy - my) / (qx - mx || 1e-300);
        if (t < tanMin) { tanMin = t; pick = k; }
      }
    }
    // a vertex already used by earlier bridges exists in several copies: take the copy whose interior wedge faces the hole
    const px = X[pick], py = Y[pick];
    for (let k = 0; k < np; k++) {
      if (X[k] !== px || Y[k] !== py) continue;
      const a = (k + np - 1) % np, c = (k + 1) % np, l1 = (px - X[a]) * (my - Y[a]) - (py - Y[a]) * (mx - X[a]), l2 = (X[c] - px) * (my - py) - (Y[c] - py) * (mx - px);
      if ((px - X[a]) * (Y[c] - Y[a]) - (py - Y[a]) * (X[c] - X[a]) > 0 ? l1 > 0 && l2 > 0 : l1 > 0 || l2 > 0) { pick = k; break; }
    }
    const ord = [];
    for (let k = 0; k <= h.length; k++) ord.push(h[(m + k) % h.length]);
    insert(X, pick + 1, [...ord.map((p) => p[0]), px]); insert(Y, pick + 1, [...ord.map((p) => p[1]), py]); insert(ID, pick + 1, [...ord.map((p) => p[2]), ID[pick]]);
  }
  const n = X.length, out = [];
  if (n < 3) return out;
  let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity;
  for (let k = 0; k < n; k++) { if (X[k] < x0) x0 = X[k]; if (X[k] > x1) x1 = X[k]; if (Y[k] < y0) y0 = Y[k]; if (Y[k] > y1) y1 = Y[k]; }
  const eps = 1e-12 * ((x1 - x0) ** 2 + (y1 - y0) ** 2 + 1e-300), c2 = (a, b, c) => (X[b] - X[a]) * (Y[c] - Y[a]) - (Y[b] - Y[a]) * (X[c] - X[a]);
  // uniform bucket grid so that the "no other vertex inside the ear" test only looks at nearby vertices
  const G = Math.max(1, Math.min(1024, Math.ceil(Math.sqrt(n / 2)))), gw = (x1 - x0) / G || 1, gh = (y1 - y0) / G || 1, cells = Array.from({ length: G * G }, () => []);
  const gx = (x) => Math.max(0, Math.min(G - 1, Math.floor((x - x0) / gw))), gy = (y) => Math.max(0, Math.min(G - 1, Math.floor((y - y0) / gh)));
  for (let k = 0; k < n; k++) cells[gy(Y[k]) * G + gx(X[k])].push(k);
  const alive = new Uint8Array(n).fill(1), prev = new Int32Array(n), next = new Int32Array(n);
  for (let k = 0; k < n; k++) { prev[k] = (k + n - 1) % n; next[k] = (k + 1) % n; }
  const blocked = (a, i, c) => {
    const ax = X[a], ay = Y[a], bx = X[i], by = Y[i], cx = X[c], cy = Y[c];
    const ca = gx(Math.min(ax, bx, cx)), cb = gx(Math.max(ax, bx, cx)), ra = gy(Math.min(ay, by, cy)), rb = gy(Math.max(ay, by, cy));
    for (let yy = ra; yy <= rb; yy++) for (let xx = ca; xx <= cb; xx++) {
      const list = cells[yy * G + xx];
      for (let q = 0; q < list.length; q++) {
        const k = list[q];
        if (!alive[k] || k === a || k === i || k === c) continue;
        const qx = X[k], qy = Y[k];
        if ((qx === ax && qy === ay) || (qx === bx && qy === by) || (qx === cx && qy === cy)) continue;
        if ((bx - ax) * (qy - ay) - (by - ay) * (qx - ax) >= -eps && (cx - bx) * (qy - by) - (cy - by) * (qx - bx) >= -eps && (ax - cx) * (qy - cy) - (ay - cy) * (qx - cx) >= -eps) return true;
      }
    }
    return false;
  };
  let left = n, i = 0;
  while (left > 3) {
    let found = false, flat = i, flatV = Infinity;
    for (let t = 0; t < left; t++, i = next[i]) {
      const cr = c2(prev[i], i, next[i]);
      if (Math.abs(cr) < flatV) { flatV = Math.abs(cr); flat = i; }
      if (cr > eps && !blocked(prev[i], i, next[i])) { found = true; break; }
    }
    if (!found) i = flat;                       // degenerate ring: clip the flattest corner so the loop always ends
    const a = prev[i], c = next[i];
    out.push(ID[a], ID[i], ID[c]);
    alive[i] = 0; next[a] = c; prev[c] = a; left--; i = a;
  }
  out.push(ID[prev[i]], ID[i], ID[next[i]]);
  return out;
}
const ringArea = (r) => { let s = 0; for (let k = 0; k < r.length; k++) { const a = r[k], b = r[(k + 1) % r.length]; s += a[0] * b[1] - b[0] * a[1]; } return s / 2; };
function newell(pts) {
  const n = [0, 0, 0];
  for (let k = 0; k < pts.length; k++) { const a = pts[k], b = pts[(k + 1) % pts.length]; n[0] += (a[1] - b[1]) * (a[2] + b[2]); n[1] += (a[2] - b[2]) * (a[0] + b[0]); n[2] += (a[0] - b[0]) * (a[1] + b[1]); }
  return n;
}
const frameFromNormal = (z) => { const x = unit(cross(Math.abs(z[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0], z)); return { x, y: cross(z, x), z }; };
/** Planar face (3-D loops, optional plane frame) -> triangles appended to tri. */
function planarFace(loops, pl, sense, tri) {
  loops = loops.filter((l) => l.pts.length >= 3);
  if (!loops.length) return;
  let fr = pl, flip = !sense;
  if (!fr) {
    const l0 = loops.find((l) => l.outer) || loops[0], z = unit(newell(l0.pts));
    if (!vlen(z)) return;
    fr = { o: l0.pts[0], ...frameFromNormal(z) }; flip = false;
  }
  const P3 = [], rings = loops.map((l) => { const r = l.pts.map((p) => { const d = sub(p, fr.o); P3.push(p); return [dot(d, fr.x), dot(d, fr.y), P3.length - 1]; }); return { r, a: ringArea(r), outer: !!l.outer }; });
  let oi = rings.findIndex((r) => r.outer);
  if (oi < 0) { oi = 0; rings.forEach((r, k) => { if (Math.abs(r.a) > Math.abs(rings[oi].a)) oi = k; }); }
  const outer = rings[oi].a < 0 ? rings[oi].r.reverse() : rings[oi].r;
  const holes = rings.filter((_, k) => k !== oi).map((h) => (h.a > 0 ? h.r.reverse() : h.r));
  const idx = triangulate(outer, holes);
  for (let k = 0; k < idx.length; k += 3) for (const q of flip ? [idx[k], idx[k + 2], idx[k + 1]] : [idx[k], idx[k + 1], idx[k + 2]]) tri.push(P3[q][0], P3[q][1], P3[q][2]);
}
/**
 * Face on a surface of revolution about the frame's z axis with radius r(v) = R + v·tanA (cylinder or cone).
 * Boundary loops are unwrapped into (angle, height), wrapped into one 2π window and filled strip by strip.
 */
function revolvedFace(loops, pl, R, tanA, sense, tri) {
  const { o, x, y, z } = pl, eps = 1e-9, edges = [];
  const at = (t, v) => { const r = Math.abs(R + v * tanA), c = Math.cos(t) * r, s = Math.sin(t) * r; return [o[0] + c * x[0] + s * y[0] + v * z[0], o[1] + c * x[1] + s * y[1] + v * z[1], o[2] + c * x[2] + s * y[2] + v * z[2]]; };
  let T0 = null, periodic = 0;
  const add = (t1, v1, P1, t2, v2, P2) => {
    if (Math.abs(t1 - t2) < eps) return;
    if (t1 > t2) [t1, v1, P1, t2, v2, P2] = [t2, v2, P2, t1, v1, P1];
    const k = Math.floor((t1 - T0 + eps) / TAU), top = T0 + TAU;
    t1 -= k * TAU; t2 -= k * TAU;
    if (t2 > top + eps) { const vm = v1 + ((v2 - v1) * (top - t1)) / (t2 - t1); edges.push([t1, v1, P1, top, vm, null], [T0, vm, null, t2 - TAU, v2, P2]); }
    else edges.push([t1, v1, P1, Math.min(t2, top), v2, P2]);
  };
  for (const l of loops) {
    if (l.single || l.pts.length < 2) continue;
    const prm = [];
    let prev = null;
    for (const p of l.pts) {
      const d = sub(p, o), a = dot(d, x), b = dot(d, y), v = dot(d, z);
      let t = Math.hypot(a, b) <= 1e-9 * (Math.abs(R) + Math.abs(v * tanA)) ? NaN : Math.atan2(b, a);
      if (t !== t) t = prev ?? NaN; else if (prev !== null) t += TAU * Math.round((prev - t) / TAU);
      if (t === t) prev = t;
      prm.push([t, v, p]);
    }
    const f = prm.find((q) => q[0] === q[0]);
    if (!f) continue;
    for (const q of prm) if (q[0] !== q[0]) q[0] = f[0];
    if (T0 === null) T0 = prm[0][0];
    const last = prm[prm.length - 1], close = prm[0][0] + TAU * Math.round((last[0] - prm[0][0]) / TAU);
    if (Math.abs(close - prm[0][0]) > 1) periodic++;
    for (let k = 0; k + 1 < prm.length; k++) add(prm[k][0], prm[k][1], prm[k][2], prm[k + 1][0], prm[k + 1][1], prm[k + 1][2]);
    add(last[0], last[1], last[2], close, prm[0][1], prm[0][2]);
  }
  if (T0 === null) return;
  if (periodic % 2 === 1) {                     // a single closed parallel: the band must end at an apex (vertex loop)
    const s = loops.find((l) => l.single);
    if (!s) return;
    const va = dot(sub(s.pts[0], o), z);
    edges.push([T0, va, null, T0 + TAU, va, null]);
  }
  const all = [];
  for (const e of edges) all.push(e[0], e[3]);
  all.sort((a, b) => a - b);
  const cols = [];
  for (const c of all) {
    const p = cols.length ? cols[cols.length - 1] : null;
    if (p === null) cols.push(c);
    else if (c - p > eps) { const m = Math.ceil((c - p) / (TAU / ARC_N) - 1e-6); for (let q = 1; q < m; q++) cols.push(p + ((c - p) * q) / m); cols.push(c); }
  }
  const eq = (a, b) => a[0] === b[0] && a[1] === b[1] && a[2] === b[2];
  const emit = (a, b, c) => { if (eq(a, b) || eq(b, c) || eq(a, c)) return; for (const p of sense ? [a, b, c] : [a, c, b]) tri.push(p[0], p[1], p[2]); };
  for (let k = 0; k + 1 < cols.length; k++) {
    const ca = cols[k], cb = cols[k + 1], sp = [];
    for (const e of edges) {
      if (!(e[0] <= ca + eps && e[3] >= cb - eps)) continue;
      const hgt = (c) => e[1] + ((e[4] - e[1]) * (c - e[0])) / (e[3] - e[0]);
      const pnt = (c) => (Math.abs(c - e[0]) < eps && e[2] ? e[2] : Math.abs(c - e[3]) < eps && e[5] ? e[5] : at(c, hgt(c)));
      sp.push({ m: hgt(0.5 * (ca + cb)), pa: pnt(ca), pb: pnt(cb) });
    }
    sp.sort((p, q) => p.m - q.m);
    for (let q = 0; q + 1 < sp.length; q += 2) { const lo = sp[q], hi = sp[q + 1]; emit(lo.pa, lo.pb, hi.pb); emit(lo.pa, hi.pb, hi.pa); }
  }
}

// ---- STEP physical file (ISO 10303-21) parsing; shared by STEP and IFC -------------------------------
function p21Entity(s, id) {
  const n = s.length, numRe = /[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?/y, idRe = /[A-Za-z_][A-Za-z0-9_]*/y;
  let p = 0;
  const ws = () => { while (p < n && s.charCodeAt(p) <= 32) p++; };
  const bad = () => fail(`Malformed STEP entity #${id}.`);
  const list = (depth) => {
    if (depth > 60) fail('A STEP entity is nested too deeply.');
    p++;
    const out = [];
    ws();
    if (s[p] === ')') { p++; return out; }
    for (;;) { out.push(value(depth)); ws(); const c = s[p++]; if (c === ')') return out; if (c !== ',') bad(); }
  };
  const value = (depth) => {
    ws();
    const c = s[p];
    if (c === '(') return list(depth + 1);
    if (c === '#') { numRe.lastIndex = p + 1; const m = numRe.exec(s); if (!m) bad(); p = numRe.lastIndex; return { r: +m[0] }; }
    if (c === "'") { let q = p + 1, out = ''; for (;;) { const e = s.indexOf("'", q); if (e < 0) bad(); out += s.slice(q, e); if (s[e + 1] === "'") { out += "'"; q = e + 2; } else { p = e + 1; return out; } } }
    if (c === '$' || c === '*') { p++; return null; }
    if (c === '"') { const e = s.indexOf('"', p + 1); if (e < 0) bad(); p = e + 1; return null; }
    if (c === '.' && /[A-Za-z_]/.test(s[p + 1] || '')) { const e = s.indexOf('.', p + 1); if (e < 0) bad(); const v = { e: s.slice(p + 1, e).toUpperCase() }; p = e + 1; return v; }
    numRe.lastIndex = p;
    let m = numRe.exec(s);
    if (m) { p = numRe.lastIndex; return +m[0]; }
    idRe.lastIndex = p;
    m = idRe.exec(s);
    if (!m) bad();
    p = idRe.lastIndex; ws();
    return s[p] === '(' ? { t: m[0].toUpperCase(), a: list(depth + 1) } : m[0];
  };
  const ident = () => { ws(); idRe.lastIndex = p; const m = idRe.exec(s); if (!m) bad(); p = idRe.lastIndex; ws(); if (s[p] !== '(') bad(); return m[0].toUpperCase(); };
  ws();
  if (s[p] === '(') {
    p++;
    const parts = Object.create(null);
    for (ws(); p < n && s[p] !== ')'; ws()) { const nm = ident(); parts[nm] = list(0); }
    return { id, type: 'COMPLEX', args: null, parts };
  }
  const type = ident();
  return { id, type, args: list(0), parts: null };
}
function parseP21(text) {
  if (!/ISO-10303-21/i.test(text.slice(0, 600))) fail('Not a STEP physical file (the ISO-10303-21 header is missing).');
  if (text.indexOf('/*') >= 0) text = text.replace(/\/\*[\s\S]*?\*\//g, ' ');
  const span = new Map(), byType = new Map(), cache = new Map(), stat = { malformed: 0 }, n = text.length, typeRe = /\s*([A-Za-z0-9_]+)\s*\(/y, partRe = /([A-Za-z_][A-Za-z0-9_]*)\s*\(/g;
  let i = text.indexOf('ENDSEC');
  if (i < 0) fail('The STEP file has no DATA section.');
  const idx = (t, id) => { const l = byType.get(t); if (l) l.push(id); else byType.set(t, [id]); };
  for (;;) {
    const h = text.indexOf('#', i);
    if (h < 0) break;
    let j = h + 1, c;
    while ((c = text.charCodeAt(j)) >= 48 && c <= 57) j++;
    const id = +text.slice(h + 1, j);
    while (text.charCodeAt(j) <= 32) j++;
    if (j === h + 1 || text[j] !== '=') { i = j + (j === h + 1 ? 0 : 1); if (i <= h) i = h + 1; continue; }
    let k = j + 1, q = false;
    for (; k < n; k++) { const ch = text.charCodeAt(k); if (ch === 39) q = !q; else if (ch === 59 && !q) break; }
    if (k >= n) break;
    span.set(id, [j + 1, k]);
    typeRe.lastIndex = j + 1;
    const m = typeRe.exec(text);
    if (m) idx(m[1].toUpperCase(), id);
    else { const body = text.slice(j + 1, Math.min(k, j + 4000)); let pm; partRe.lastIndex = 0; while ((pm = partRe.exec(body))) idx(pm[1].toUpperCase(), id); }
    i = k + 1;
  }
  if (!span.size) fail('The STEP file holds no entities.');
  const get = (id) => {
    let e = cache.get(id);
    if (e !== undefined) return e;
    const sp = span.get(id);
    try { e = sp ? p21Entity(text.slice(sp[0], sp[1]), id) : null; } catch (err) { if (!err || !err.user) throw err; e = null; stat.malformed++; }
    cache.set(id, e);
    return e;
  };
  const ids = (...types) => types.flatMap((t) => byType.get(t) || []);
  const ref = (v) => (v && typeof v === 'object' && 'r' in v ? get(v.r) : null);
  const part = (e, name) => (!e ? null : e.type === name ? e.args : e.parts ? e.parts[name] || null : null);
  return { get, ids, ref, part, byType, stat, count: span.size, schema: (/FILE_SCHEMA\s*\(\s*\(\s*'([^']*)'/i.exec(text.slice(0, i)) || [, ''])[1] };
}
const P21_PREFIX = { EXA: 1e18, PETA: 1e15, TERA: 1e12, GIGA: 1e9, MEGA: 1e6, KILO: 1e3, HECTO: 100, DECA: 10, DECI: 0.1, CENTI: 0.01, MILLI: 1e-3, MICRO: 1e-6, NANO: 1e-9 };
const UNIT_NAME = { 1: 'm', 0.001: 'mm', 0.01: 'cm', 0.1: 'dm', 1000: 'km', 0.000001: 'µm' };

async function readSTEP(ctx) {
  const db = parseP21(await ctx.text()), { get, ids, ref, part } = db, warnings = [], stats = { entities: db.count, schema: db.schema };
  const flag = (v, dflt = true) => (v && v.e ? v.e === 'T' : dflt);
  const lastList = (a) => (a ? a.filter(Array.isArray).pop() : null);
  const coords = (v) => { const c = lastList(part(ref(v), 'CARTESIAN_POINT')); return c && c.length >= 2 && c.every(Number.isFinite) ? [c[0], c[1], c[2] || 0] : null; };
  const ptCache = new Map();
  const pt = (v) => { if (!v || v.r === undefined) return null; let p = ptCache.get(v.r); if (p === undefined) { p = coords(v); ptCache.set(v.r, p); } return p; };
  const dir = (v) => { const c = lastList(part(ref(v), 'DIRECTION')); return c && c.length >= 2 && c.every(Number.isFinite) ? [c[0], c[1], c[2] || 0] : null; };
  const vertex = (v) => { const a = part(ref(v), 'VERTEX_POINT'); return a ? pt(a[a.length - 1]) : null; };
  const placement = (v) => {
    const a = part(ref(v), 'AXIS2_PLACEMENT_3D');
    if (!a) return null;
    const o = pt(a[a.length - 3]), z = unit(dir(a[a.length - 2]) || [0, 0, 1]);
    if (!o || !vlen(z)) return null;
    const xr = dir(a[a.length - 1]) || (Math.abs(z[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0]), k = dot(xr, z);
    let x = unit([xr[0] - k * z[0], xr[1] - k * z[1], xr[2] - k * z[2]]);
    if (!vlen(x)) x = frameFromNormal(z).x;
    return { o, x, y: cross(z, x), z };
  };
  // length unit of the model
  for (const id of ids('LENGTH_UNIT')) {
    const e = get(id);
    if (!e || !e.parts) continue;
    if (e.parts.CONVERSION_BASED_UNIT) { stats.units = String(e.parts.CONVERSION_BASED_UNIT[0] || '').toLowerCase() || stats.units; break; }
    const si = e.parts.SI_UNIT;
    if (si && si[1] && si[1].e === 'METRE' && !stats.units) { const f = si[0] && si[0].e ? P21_PREFIX[si[0].e] : 1; stats.units = UNIT_NAME[f] || `${f} m`; }
  }
  const tri = [];
  // (a) AP242 tessellated geometry
  let tess = 0, badIdx = 0;
  for (const id of ids('TRIANGULATED_FACE', 'COMPLEX_TRIANGULATED_FACE', 'TRIANGULATED_SURFACE_SET', 'COMPLEX_TRIANGULATED_SURFACE_SET')) {
    const e = get(id), a = e && e.args;
    if (!a || a.length < 6) continue;
    const cl = lastList(part(ref(a[1]), 'COORDINATES_LIST'));
    if (!cl) continue;
    const cx = e.type.startsWith('COMPLEX'), pn = cx ? a[a.length - 3] : a[a.length - 2], usePn = Array.isArray(pn) && pn.length > 0;
    const vtx = (i) => { const c = cl[(usePn ? pn[i - 1] : i) - 1]; return Array.isArray(c) && c.length >= 3 ? c : null; };
    const put = (i, j, k) => { const A = vtx(i), B = vtx(j), Cc = vtx(k); if (A && B && Cc) tri.push(A[0], A[1], A[2], B[0], B[1], B[2], Cc[0], Cc[1], Cc[2]); else badIdx++; };
    const lists = (v) => (Array.isArray(v) ? v.filter(Array.isArray) : []);
    if (!cx) for (const t of lists(a[a.length - 1])) put(t[0], t[1], t[2]);
    else {
      for (const s of lists(a[a.length - 2])) for (let k = 0; k + 2 < s.length; k++) if (k % 2) put(s[k + 1], s[k], s[k + 2]); else put(s[k], s[k + 1], s[k + 2]);
      for (const f of lists(a[a.length - 1])) for (let k = 1; k + 1 < f.length; k++) put(f[0], f[k], f[k + 1]);
    }
    tess++;
    if (tri.length > L.triangles * 9) fail(`The STEP tessellation has too many triangles (limit ${L.triangles}).`);
  }
  // (b, c) boundary-representation faces
  const faceIds = ids('ADVANCED_FACE', 'FACE_SURFACE', 'FACE'), skipped = {};
  let done = 0, chordEdges = 0;
  if (tess) {
    stats.tessellatedFaces = tess;
    if (faceIds.length) warnings.push(`Used the ${tess} AP242 tessellated faces stored in the file; its ${faceIds.length} exact B-rep faces were not re-tessellated.`);
    if (badIdx) warnings.push(`${badIdx} tessellated triangles referred to missing coordinates and were dropped.`);
  } else {
    const edgeCache = new Map();
    const curveOf = (e) => { for (let d = 0; e && e.args && d < 4 && ['SURFACE_CURVE', 'SEAM_CURVE', 'INTERSECTION_CURVE', 'BOUNDED_SURFACE_CURVE', 'TRIMMED_CURVE'].includes(e.type); d++) e = ref(e.args[1]); return e; };
    const edgePts = (id) => {
      let pts = edgeCache.get(id);
      if (pts) return pts;
      edgeCache.set(id, (pts = []));
      const a = part(get(id), 'EDGE_CURVE');
      if (!a) return pts;
      const A = vertex(a[1]), B = vertex(a[2]), c = curveOf(ref(a[3])), same = flag(a[4]);
      if (!A || !B) return pts;
      pts.push(A, B);
      const con = part(c, 'CIRCLE') || part(c, 'ELLIPSE'), bw = part(c, 'B_SPLINE_CURVE_WITH_KNOTS');
      if (con) {
        const pl = placement(con[1]), rx = +con[2], ry = c.type === 'ELLIPSE' ? +con[3] : rx;
        if (!pl || !(rx > 0) || !(ry > 0)) { chordEdges++; return pts; }
        const ang = (Q) => { const d = sub(Q, pl.o); return Math.atan2(dot(d, pl.y) / ry, dot(d, pl.x) / rx); }, ta = ang(A);
        let sw;
        if ((a[1] && a[2] && a[1].r === a[2].r) || vlen(sub(A, B)) <= 1e-9 * (rx + ry)) sw = same ? TAU : -TAU;
        else { sw = ang(B) - ta; if (same) { while (sw <= 0) sw += TAU; } else while (sw >= 0) sw -= TAU; }
        const n = Math.max(2, Math.ceil((Math.abs(sw) / TAU) * ARC_N - 1e-9));
        pts.length = 0;
        for (let i = 0; i <= n; i++) { const t = ta + (sw * i) / n, cs = rx * Math.cos(t), sn = ry * Math.sin(t); pts.push(i === 0 ? A : i === n ? B : [pl.o[0] + cs * pl.x[0] + sn * pl.y[0], pl.o[1] + cs * pl.x[1] + sn * pl.y[1], pl.o[2] + cs * pl.x[2] + sn * pl.y[2]]); }
      } else if (bw) {
        const simple = c.type === 'B_SPLINE_CURVE_WITH_KNOTS', base = simple ? bw.slice(1) : part(c, 'B_SPLINE_CURVE'), kn = simple ? bw.slice(6) : bw, wts = part(c, 'RATIONAL_B_SPLINE_CURVE');
        let smp = null;
        if (base && Array.isArray(base[1]) && Array.isArray(kn[0]) && Array.isArray(kn[1]) && kn[0].length === kn[1].length) {
          const cps = base[1].map(pt), knots = [];
          kn[0].forEach((m, q) => { for (let r = 0; r < m && knots.length < 1e5; r++) knots.push(kn[1][q]); });
          smp = nurbsCurve(base[0], knots, cps, wts ? wts[wts.length - 1] : null, NaN, NaN, 4 * cps.length);
        }
        if (smp && !same) smp.reverse();
        const tol = smp ? 1e-3 * (vlen(sub(A, B)) + vlen(sub(smp[0], smp[smp.length >> 1]))) + 1e-12 : 0;
        if (smp && vlen(sub(smp[0], A)) <= tol && vlen(sub(smp[smp.length - 1], B)) <= tol) { smp[0] = A; smp[smp.length - 1] = B; pts.length = 0; pts.push(...smp); } else chordEdges++;
      } else if (c && c.type !== 'LINE') chordEdges++;
      return pts;
    };
    const loopPts = (lp) => {
      const pl = part(lp, 'POLY_LOOP'), vl = part(lp, 'VERTEX_LOOP'), el = part(lp, 'EDGE_LOOP');
      if (pl) return { pts: (lastList(pl) || []).map(pt).filter(Boolean) };
      if (vl) { const p = vertex(vl[vl.length - 1]); return p ? { pts: [p], single: true } : null; }
      if (!el) return null;
      const pts = [];
      for (const oe of lastList(el) || []) {
        const o = ref(oe);
        if (!o || !o.args) continue;
        const isO = o.type === 'ORIENTED_EDGE', ep = edgePts(isO ? (o.args[3] || {}).r : o.id), seq = !isO || flag(o.args[4]) ? ep : ep.slice().reverse();
        for (let k = 0; k + 1 < seq.length; k++) pts.push(seq[k]);
      }
      return { pts };
    };
    for (const id of faceIds) {
      const e = get(id), a = e && e.args;
      if (!a) continue;
      const loops = [];
      for (const b of Array.isArray(a[1]) ? a[1] : []) {
        const be = ref(b);
        if (!be || !be.args) continue;
        const lp = loopPts(ref(be.args[1]));
        if (!lp || !lp.pts.length) continue;
        if (!flag(be.args[2])) lp.pts.reverse();
        lp.outer = be.type === 'FACE_OUTER_BOUND';
        loops.push(lp);
      }
      const surf = ref(a[2]), sense = flag(a[3]), before = tri.length, st = surf ? surf.type : 'PLANE';
      if (st === 'PLANE') planarFace(loops, surf ? placement(surf.args[1]) : null, sense, tri);
      else if (st === 'CYLINDRICAL_SURFACE' || st === 'CONICAL_SURFACE') {
        const pl = placement(surf.args[1]), R = +surf.args[2], ta = st === 'CONICAL_SURFACE' ? Math.tan(+surf.args[3]) : 0;
        if (pl && Number.isFinite(R) && Number.isFinite(ta)) revolvedFace(loops, pl, R, ta, sense, tri);
      } else { skipped[st] = (skipped[st] || 0) + 1; continue; }
      if (tri.length > before) done++; else skipped['degenerate ' + st] = (skipped['degenerate ' + st] || 0) + 1;
      if (tri.length > L.triangles * 9) fail(`The STEP model tessellates to too many triangles (limit ${L.triangles}).`);
    }
    stats.faces = faceIds.length; stats.tessellatedFaces = done;
    const nSkip = faceIds.length - done;
    if (nSkip) warnings.push(`${done} of ${faceIds.length} faces tessellated; ${nSkip} skipped (${Object.entries(skipped).map(([k, v]) => `${v} × ${k.toLowerCase().replace(/_/g, ' ')}`).join(', ')}). Export AP242 with tessellation or STL for a complete surface.`);
    if (chordEdges) warnings.push(`${chordEdges} curved edges of unsupported type were replaced by straight chords.`);
  }
  const cp = [];
  for (const id of ids('CARTESIAN_POINT')) { const p = pt({ r: id }); if (p) cp.push(p[0], p[1], p[2]); }
  const pb = cp.length ? bboxOf(cp) : null;
  if (ids('NEXT_ASSEMBLY_USAGE_OCCURRENCE').length) warnings.push('The file is an assembly: component placements are not applied, parts are shown in their own coordinates.');
  if (!stats.units) warnings.push('No length unit was found in the file; coordinates are returned as written.');
  if (db.stat.malformed) warnings.push(`${db.stat.malformed} malformed entities were ignored.`);
  if (!tri.length) {
    if (!cp.length) fail('The STEP file holds no geometry that can be read.');
    warnings.push('No faces could be tessellated; the Cartesian points of the model are returned as a point set (extent only).');
    return pointsGeom(cp, { warnings, stats });
  }
  const g = meshGeom(tri, { warnings, stats });
  if (pb && (tess ? false : done < faceIds.length)) g.bbox = { min: g.bbox.min.map((v, k) => Math.min(v, pb.min[k])), max: g.bbox.max.map((v, k) => Math.max(v, pb.max[k])) };
  if (pb) stats.pointBbox = pb;
  return g;
}

// ---- IGES (fixed 80-column ASCII) --------------------------------------------------------------------
function igesTokens(s, pd, rd) {
  const out = [], n = s.length;
  let p = 0;
  while (p < n) {
    while (p < n && s[p] === ' ') p++;
    const m = /^(\d+)H/.exec(s.slice(p, p + 12));
    if (m) { const st = p + m[0].length; out.push(s.slice(st, st + +m[1])); p = st + +m[1]; while (p < n && s[p] !== pd && s[p] !== rd) p++; }
    else { let q = p; while (q < n && s[q] !== pd && s[q] !== rd) q++; out.push(s.slice(p, q).trim()); p = q; }
    if (p >= n || s[p] === rd) break;
    p++;
  }
  return out;
}
const IGES_UNITS = { 1: 'inch', 2: 'mm', 4: 'ft', 5: 'mile', 6: 'm', 7: 'km', 8: 'mil', 9: 'µm', 10: 'cm', 11: 'µin' };
const IGES_SURF = new Set([108, 114, 118, 120, 122, 128, 140, 143, 144, 150, 152, 154, 156, 158, 160, 162, 164, 168, 180, 182, 184, 186, 190, 192, 194, 196, 198, 510]);
async function readIGES(ctx) {
  const sec = { S: [], G: [], D: [], P: [] }, text = await ctx.text();
  for (const ln of splitLines(text)) {
    if (ln.length < 73) continue;
    const c = ln[72];
    if (sec[c]) sec[c].push(ln); else if (c === 'C' || c === 'B') fail('Compressed or binary IGES is not supported; save as fixed-format ASCII IGES.');
  }
  if (!sec.D.length || !sec.P.length) fail('Not a fixed-format ASCII IGES file (no directory or parameter section found).');
  const g = sec.G.map((l) => l.slice(0, 72)).join('');
  let pd = ',', rd = ';', gp = 0;
  if (g.startsWith('1H')) { pd = g[2]; gp = 4; } else if (g[0] === ',') gp = 1;
  if (g.startsWith('1H', gp)) rd = g[gp + 2];
  const gt = igesTokens(g, pd, rd), warnings = [], stats = { entities: sec.D.length >> 1, units: IGES_UNITS[parseInt(gt[13], 10)] || (parseInt(gt[13], 10) === 3 ? String(gt[14] || '').toLowerCase() : undefined) };
  const pBySeq = new Map();
  for (const ln of sec.P) pBySeq.set(parseInt(ln.slice(73, 80), 10), ln.slice(0, 64));
  const fld = (l, k) => l.slice(8 * k, 8 * k + 8).trim();
  const params = (de) => {
    const l1 = sec.D[de], l2 = sec.D[de + 1], p0 = parseInt(fld(l1, 1), 10), cnt = Math.min(Math.max(1, parseInt(fld(l2, 3), 10) || 1), 200000);
    let s = '';
    for (let k = 0; k < cnt; k++) s += pBySeq.get(p0 + k) || '';
    return igesTokens(s, pd, rd).map((t) => parseFloat(String(t).replace(/[dD]/, 'e')));
  };
  const xform = (p, seq, depth = 0) => {
    const de = seq - 1;
    if (!(seq > 0) || depth > 8 || de + 1 >= sec.D.length || parseInt(fld(sec.D[de], 0), 10) !== 124) return p;
    const v = params(de);
    if (v.length < 13 || !v.slice(1, 13).every(Number.isFinite)) return p;
    const q = [v[1] * p[0] + v[2] * p[1] + v[3] * p[2] + v[4], v[5] * p[0] + v[6] * p[1] + v[7] * p[2] + v[8], v[9] * p[0] + v[10] * p[1] + v[11] * p[2] + v[12]];
    return xform(q, parseInt(fld(sec.D[de], 6), 10) || 0, depth + 1);
  };
  const curves = [], marks = [], counts = {};
  let surfaces = 0, bad = 0, total = 0;
  for (let de = 0; de + 1 < sec.D.length; de += 2) {
    const l1 = sec.D[de], type = parseInt(fld(l1, 0), 10), status = l1.slice(64, 72), xf = parseInt(fld(l1, 6), 10) || 0, form = parseInt(fld(sec.D[de + 1], 4), 10) || 0;
    if (IGES_SURF.has(type)) { surfaces++; continue; }
    if (![100, 106, 110, 112, 116, 126].includes(type) || status.slice(0, 2) === '01' || status.slice(4, 6) === '05') continue;
    const v = params(de);
    let pts = null;
    if (type === 110 && v.length >= 7) pts = [[v[1], v[2], v[3]], [v[4], v[5], v[6]]];
    else if (type === 100 && v.length >= 8) {
      const [, zt, cx, cy, sx, sy, ex, ey] = v, r = Math.hypot(sx - cx, sy - cy), a0 = Math.atan2(sy - cy, sx - cx);
      let sw = Math.atan2(ey - cy, ex - cx) - a0;
      if (sw <= 1e-12) sw += TAU;
      const n = Math.max(2, Math.ceil((sw / TAU) * ARC_N - 1e-9));
      pts = [];
      for (let i = 0; i <= n; i++) pts.push(i === n && sw < TAU - 1e-9 ? [ex, ey, zt] : [cx + r * Math.cos(a0 + (sw * i) / n), cy + r * Math.sin(a0 + (sw * i) / n), zt]);
    } else if (type === 106 && v.length >= 3) {
      const ip = v[1], n = Math.min(v[2], v.length);
      pts = [];
      if (ip === 1) for (let i = 0; i < n; i++) pts.push([v[4 + 2 * i], v[5 + 2 * i], v[3]]);
      else { const w = ip === 3 ? 6 : 3; for (let i = 0; i < n; i++) pts.push([v[3 + w * i], v[4 + w * i], v[5 + w * i]]); }
      if (form >= 1 && form <= 3) { for (const p of pts) marks.push(xform(p, xf)); pts = null; counts[type] = (counts[type] || 0) + 1; }
    } else if (type === 116 && v.length >= 4) { marks.push(xform([v[1], v[2], v[3]], xf)); counts[116] = (counts[116] || 0) + 1; }
    else if (type === 112 && v.length >= 6) {
      const ns = v[4], t0 = 5, c0 = 6 + ns;
      pts = [];
      if (Number.isInteger(ns) && ns > 0 && v.length >= c0 + 12 * ns) for (let s = 0; s < ns; s++) {
        const h = v[t0 + s + 1] - v[t0 + s], m = 8;
        for (let i = s ? 1 : 0; i <= m; i++) { const u = (h * i) / m, c = c0 + 12 * s, ev = (o) => v[c + o] + u * (v[c + o + 1] + u * (v[c + o + 2] + u * v[c + o + 3])); pts.push([ev(0), ev(4), ev(8)]); }
      }
    } else if (type === 126 && v.length >= 8) {
      const K = v[1], M = v[2], nk = K + M + 2, k0 = 7, w0 = k0 + nk, c0 = w0 + K + 1, e0 = c0 + 3 * (K + 1);
      if (Number.isInteger(K) && Number.isInteger(M) && K >= M && M >= 1 && v.length >= e0 + 2) {
        const cps = [];
        for (let i = 0; i <= K; i++) cps.push([v[c0 + 3 * i], v[c0 + 3 * i + 1], v[c0 + 3 * i + 2]]);
        const w = v.slice(w0, w0 + K + 1);
        pts = nurbsCurve(M, v.slice(k0, k0 + nk), cps, w.every((x) => x === w[0]) ? null : w, v[e0], v[e0 + 1], 8 * (K + 1));
      }
    }
    if (pts === null && (type === 116 || (type === 106 && form >= 1 && form <= 3))) continue;
    if (!pts || pts.length < 2 || pts.some((p) => !p.every(Number.isFinite))) { bad++; continue; }
    counts[type] = (counts[type] || 0) + 1;
    total += pts.length;
    if (total > 4e6) fail('The IGES wireframe is too large.');
    curves.push(pts.map((p) => xform(p, xf)));
  }
  if (!curves.length) fail('The IGES file holds no readable curve entities' + (surfaces ? ` (only ${surfaces} surface or solid entities, which are not tessellated)` : '') + '. For solids export STEP (AP214 or AP242) or STL.');
  const flat = curves.flat().concat(marks).flat(), bb = bboxOf(flat), ext = bb.max.map((v, k) => v - bb.min[k]);
  const drop = ext[2] <= ext[0] && ext[2] <= ext[1] ? 2 : ext[1] <= ext[0] ? 1 : 0, ax = [0, 1, 2].filter((k) => k !== drop);
  const polys = curves.map((c) => { const cl = c.length > 2 && vlen(sub(c[0], c[c.length - 1])) <= 1e-9 * (vlen(ext) + 1e-300), q = cl ? c.slice(0, -1) : c; return { x: q.map((p) => p[ax[0]]), y: q.map((p) => p[ax[1]]), closed: cl }; });
  warnings.push(`Wire-frame of ${curves.length} curves projected onto the ${'xyz'[ax[0]]}-${'xyz'[ax[1]]} plane (the ${'xyz'[drop]} extent of ${+ext[drop].toPrecision(4)} was dropped).`);
  warnings.push(surfaces ? `${surfaces} surface or solid entities (trimmed, NURBS or analytic surfaces) are not tessellated; export STEP or STL for the solid shape.` : 'Trimmed surfaces are not tessellated by this reader; only curves are returned.');
  if (bad) warnings.push(`${bad} curve entities were malformed and skipped.`);
  return polyGeom(polys, { warnings, markers: marks.map((p) => [p[ax[0]], p[ax[1]]]), stats: { ...stats, curves: curves.length, points: marks.length, surfaces, byEntity: counts, droppedAxis: 'xyz'[drop], bbox3: bb } });
}

// ---- Surface formats -----------------------------------------------------------------------------------
const pushFan = (tri, xyz, ids) => {
  const n = xyz.length / 3;
  for (const i of ids) if (!Number.isInteger(i) || i < 0 || i >= n) return false;
  for (let k = 1; k + 1 < ids.length; k++) for (const q of [ids[0], ids[k], ids[k + 1]]) tri.push(xyz[3 * q], xyz[3 * q + 1], xyz[3 * q + 2]);
  if (tri.length > L.triangles * 9) fail(`Mesh has too many triangles (limit ${L.triangles}).`);
  return true;
};
/** Index list with -1 separators -> polygons. */
const splitPolys = (idx) => { const out = []; let cur = []; for (const v of idx) { if (v < 0) { if (cur.length > 2) out.push(cur); cur = []; } else cur.push(v); } if (cur.length > 2) out.push(cur); return out; };

async function readPLY(ctx) {
  const u8 = await ctx.bytes(), head = latin1.decode(u8.subarray(0, Math.min(u8.length, 1 << 17))), m = /end_header[ \t]*(\r\n|\n|\r)/.exec(head);
  if (!/^ply\s/.test(head) || !m) fail('Not a PLY file (header not found).');
  let off = m.index + m[0].length, fmt = '';
  const elems = [];
  for (const ln of splitLines(head.slice(0, m.index))) {
    const t = ln.trim().split(/\s+/);
    if (t[0] === 'format') fmt = t[1];
    else if (t[0] === 'element') elems.push({ name: t[1], count: +t[2], props: [] });
    else if (t[0] === 'property' && elems.length) {
      const p = t[1] === 'list' ? { list: true, ct: dtypeOf(t[2]), it: dtypeOf(t[3]), name: t[4] } : { t: dtypeOf(t[1]), name: t[2] };
      if (!(p.list ? p.ct && p.it : p.t)) fail('The PLY header uses an unknown property type.');
      elems[elems.length - 1].props.push(p);
    }
  }
  const ascii = fmt === 'ascii', le = fmt === 'binary_little_endian';
  if (!ascii && !le && fmt !== 'binary_big_endian') fail('The PLY header has no valid format line.');
  const dv = new DataView(u8.buffer, u8.byteOffset, u8.byteLength);
  const RD = { u1: (o) => dv.getUint8(o), i1: (o) => dv.getInt8(o), u2: (o) => dv.getUint16(o, le), i2: (o) => dv.getInt16(o, le), u4: (o) => dv.getUint32(o, le), i4: (o) => dv.getInt32(o, le), f4: (o) => dv.getFloat32(o, le), f8: (o) => dv.getFloat64(o, le), i8: (o) => Number(dv.getBigInt64(o, le)), u8: (o) => Number(dv.getBigUint64(o, le)) };
  const toks = ascii ? latin1.decode(u8.subarray(off)).split(/\s+/).filter(Boolean) : null;
  let tp = 0;
  const read = ascii ? () => +toks[tp++] : (t) => { if (off + DT[t][1] > u8.length) fail('The PLY data is truncated.'); const v = RD[t](off); off += DT[t][1]; return v; };
  let verts = new Float64Array(0);
  const tri = [];
  let faces = 0, badFaces = 0;
  for (const e of elems) {
    if (!Number.isInteger(e.count) || e.count < 0 || e.count > u8.length) fail('The PLY header declares an impossible element count.');
    const isV = e.name === 'vertex', isF = e.name === 'face', ci = ['x', 'y', 'z'].map((c) => e.props.findIndex((p) => p.name === c));
    if (isV) verts = new Float64Array(e.count * 3);
    for (let r = 0; r < e.count; r++) for (let k = 0; k < e.props.length; k++) {
      const p = e.props[k];
      if (!p.list) { const v = read(p.t); if (isV) { if (k === ci[0]) verts[3 * r] = v; else if (k === ci[1]) verts[3 * r + 1] = v; else if (k === ci[2]) verts[3 * r + 2] = v; } continue; }
      const n = read(p.ct);
      if (!(n >= 0 && n <= 65536)) fail('The PLY data is corrupt or truncated.');
      const ids = [];
      for (let i = 0; i < n; i++) ids.push(read(p.it));
      if (isF && /^vertex_ind/.test(p.name)) { faces++; if (!pushFan(tri, verts, ids)) badFaces++; }
    }
  }
  if (ascii && tp > toks.length) fail('The PLY data is truncated.');
  const warnings = badFaces ? [`${badFaces} faces referred to missing vertices and were dropped.`] : [];
  if (!tri.length) { warnings.push('The PLY file has no faces; its vertices are returned as a point cloud.'); return pointsGeom(verts, { warnings }); }
  return meshGeom(tri, { warnings, stats: { vertices: verts.length / 3, faces } });
}

async function readOFF(ctx) {
  const ls = splitLines(await ctx.text()).map((l) => l.replace(/#.*/, '').trim()).filter(Boolean), hm = /^(?:ST)?(?:C)?(?:N)?(4)?(n)?OFF\b\s*(.*)$/.exec(ls[0] || '');
  if (!hm) fail('Not an OFF file (the OFF keyword is missing).');
  if (hm[1] || hm[2]) fail('Only 3-D OFF files are supported.');
  let i = hm[3] ? 1 : 2;
  const [nv, nf] = (hm[3] || ls[1] || '').split(/\s+/).map(Number);
  if (!Number.isInteger(nv) || !Number.isInteger(nf) || nv < 0 || nf < 0 || nv + nf > ls.length) fail('The OFF file is truncated or its counts are invalid.');
  const xyz = new Float64Array(nv * 3), tri = [];
  for (let k = 0; k < nv; k++) { const t = ls[i++].split(/\s+/); xyz[3 * k] = +t[0]; xyz[3 * k + 1] = +t[1]; xyz[3 * k + 2] = +t[2]; }
  let bad = 0;
  for (let k = 0; k < nf; k++) { const t = ls[i++].split(/\s+/).map(Number), n = t[0]; if (!(n >= 3) || !pushFan(tri, xyz, t.slice(1, 1 + n))) bad++; }
  return meshGeom(tri, { warnings: bad ? [`${bad} faces were malformed and dropped.`] : [], stats: { vertices: nv, faces: nf } });
}

async function read3MF(ctx) {
  const files = await unzipFiles(await ctx.buf(), (n) => /\.model$/i.test(n)), key = [...files.keys()].find((n) => /3dmodel\.model$/i.test(n)) || [...files.keys()][0];
  if (!key) fail('No 3-D model part was found in the 3MF archive.');
  const model = xAll(parseXML(utf8.decode(files.get(key))), 'model')[0];
  if (!model) fail('The 3MF model part has no <model> element.');
  const objs = new Map(), tri = [], mat = (s) => { const m = numsOf(s || ''); return m.length === 12 && m.every(Number.isFinite) ? m : null; };
  for (const o of xAll(model, 'object')) {
    const mesh = xKid(o, 'mesh'), ob = { v: null, t: null, comps: xAll(o, 'component').map((c) => ({ id: c.attrs.objectid, m: mat(c.attrs.transform) })) };
    if (mesh) {
      const vs = xAll(mesh, 'vertex'), ts = xAll(mesh, 'triangle');
      ob.v = new Float64Array(vs.length * 3); ob.t = new Int32Array(ts.length * 3);
      vs.forEach((v, k) => { ob.v[3 * k] = +v.attrs.x; ob.v[3 * k + 1] = +v.attrs.y; ob.v[3 * k + 2] = +v.attrs.z; });
      ts.forEach((t, k) => { ob.t[3 * k] = +t.attrs.v1; ob.t[3 * k + 1] = +t.attrs.v2; ob.t[3 * k + 2] = +t.attrs.v3; });
    }
    objs.set(o.attrs.id, ob);
  }
  let emits = 0;
  const emit = (id, chain, depth) => {
    const ob = objs.get(id);
    if (!ob || depth > 16) return;
    if (++emits > 5e4) fail('The 3MF model instantiates too many components.');
    if (ob.v) {
      if (tri.length + ob.t.length * 3 > L.triangles * 9) fail(`Mesh has too many triangles (limit ${L.triangles}).`);
      for (let k = 0; k < ob.t.length; k++) {
        const q = ob.t[k];
        if (q < 0 || q >= ob.v.length / 3) fail('A 3MF triangle refers to a missing vertex.');
        let x = ob.v[3 * q], y = ob.v[3 * q + 1], z = ob.v[3 * q + 2];
        for (const m of chain) [x, y, z] = [x * m[0] + y * m[3] + z * m[6] + m[9], x * m[1] + y * m[4] + z * m[7] + m[10], x * m[2] + y * m[5] + z * m[8] + m[11]];
        tri.push(x, y, z);
      }
    }
    for (const c of ob.comps) emit(c.id, c.m ? [c.m, ...chain] : chain, depth + 1);
  };
  const items = xAll(model, 'item');
  if (items.length) for (const it of items) { const m = mat(it.attrs.transform); emit(it.attrs.objectid, m ? [m] : [], 0); }
  else for (const id of objs.keys()) emit(id, [], 0);
  return meshGeom(tri, { stats: { objects: objs.size, units: model.attrs.unit || 'millimeter' } });
}

async function readAMF(ctx) {
  let u8 = await ctx.bytes();
  if (isZip(u8)) { const f = await unzipFiles(await ctx.buf(), (n) => /\.amf$/i.test(n)); u8 = f.values().next().value; if (!u8) fail('No .amf document was found inside the archive.'); }
  const doc = parseXML(utf8.decode(u8)), tri = [], root = xAll(doc, 'amf')[0];
  if (!root) fail('Not an AMF file (no <amf> element).');
  const num = (n, k) => { const c = n && xKid(n, k); return c ? +c.text : NaN; };
  for (const mesh of xAll(root, 'mesh')) {
    const vs = xAll(mesh, 'vertex'), xyz = new Float64Array(vs.length * 3);
    vs.forEach((v, k) => { const c = xKid(v, 'coordinates'); xyz[3 * k] = num(c, 'x'); xyz[3 * k + 1] = num(c, 'y'); xyz[3 * k + 2] = num(c, 'z'); });
    for (const t of xAll(mesh, 'triangle')) if (!pushFan(tri, xyz, [num(t, 'v1'), num(t, 'v2'), num(t, 'v3')])) fail('An AMF triangle refers to a missing vertex.');
  }
  return meshGeom(tri, { stats: { units: root.attrs.unit || 'millimeter' } });
}

async function readVRML(ctx) {
  const raw = await ctx.text();
  if (/^#Inventor V[\d.]+ binary/i.test(raw)) fail('Binary Open Inventor files are not supported; save as ASCII .iv, VRML, STL or OBJ.');
  if (!/^#(VRML|Inventor)/i.test(raw)) fail('Not a VRML or Open Inventor text file.');
  const src = raw.replace(/#[^\n\r]*/g, ' '), coords = [], defs = new Map(), tri = [];
  let m;
  const reC = /(?:DEF\s+(\S+)\s+)?Coordinate3?\s*\{\s*point\s*\[([^\]]*)\]/g, reF = /IndexedFaceSet\s*\{/g;
  while ((m = reC.exec(src))) { const c = { at: m.index, pts: numsOf(m[2]) }; coords.push(c); if (m[1]) defs.set(m[1], c); }
  let sets = 0;
  while ((m = reF.exec(src))) {
    let e = reF.lastIndex, depth = 1;
    for (; e < src.length && depth; e++) { const ch = src.charCodeAt(e); if (ch === 123) depth++; else if (ch === 125) depth--; }
    const block = src.slice(m.index, e), ci = /coordIndex\s*\[([^\]]*)\]/.exec(block), use = /coord\s+USE\s+(\S+)/.exec(block);
    const c = coords.find((q) => q.at > m.index && q.at < e) || (use && defs.get(use[1])) || coords.filter((q) => q.at < m.index).pop();
    if (!ci || !c) continue;
    sets++;
    for (const poly of splitPolys(numsOf(ci[1]))) pushFan(tri, c.pts, poly);
  }
  if (!sets) fail('No IndexedFaceSet with coordinates was found in the file.');
  const warnings = /\b(translation|rotation|scale|scaleFactor|MatrixTransform)\b/.test(src) ? ['Transform nodes are present but not applied; shapes are shown in their local coordinates.'] : [];
  return meshGeom(tri, { warnings, stats: { faceSets: sets } });
}

async function readX3D(ctx, doc) {
  doc = doc || parseXML(await ctx.text());
  const defs = new Map(), tri = [];
  for (const c of xAll(doc, 'coordinate')) if (c.attrs.def && c.attrs.point) defs.set(c.attrs.def, c);
  let sets = 0;
  for (const fs of xAll(doc, 'indexedfaceset', 'indexedtriangleset')) {
    let c = xAll(fs, 'coordinate')[0];
    if (c && c.attrs.use) c = defs.get(c.attrs.use);
    const idx = fs.attrs.coordindex ?? fs.attrs.index;
    if (!c || !c.attrs.point || idx === undefined) continue;
    const pts = numsOf(c.attrs.point), id = numsOf(idx);
    sets++;
    if (fs.name === 'indexedtriangleset') for (let k = 0; k + 2 < id.length; k += 3) pushFan(tri, pts, [id[k], id[k + 1], id[k + 2]]);
    else for (const poly of splitPolys(id)) pushFan(tri, pts, poly);
  }
  if (!sets) fail('No IndexedFaceSet or IndexedTriangleSet was found in the X3D file.');
  const warnings = xAll(doc, 'transform').some((t) => t.attrs.translation || t.attrs.rotation || t.attrs.scale) ? ['Transform nodes are present but not applied; shapes are shown in their local coordinates.'] : [];
  return meshGeom(tri, { warnings, stats: { faceSets: sets } });
}

async function readDAE(ctx, doc) {
  doc = doc || parseXML(await ctx.text());
  const tri = [];
  for (const mesh of xAll(doc, 'mesh')) {
    const src = new Map(), vmap = new Map();
    for (const s of xKids(mesh, 'source')) {
      const fa = xKid(s, 'float_array'), acc = xAll(s, 'accessor')[0];
      if (fa) src.set('#' + s.attrs.id, { d: numsOf(fa.text), stride: Math.max(1, +(acc && acc.attrs.stride) || 3), off: +(acc && acc.attrs.offset) || 0 });
    }
    for (const v of xKids(mesh, 'vertices')) { const inp = xKids(v, 'input').find((q) => String(q.attrs.semantic).toUpperCase() === 'POSITION'); if (inp) vmap.set('#' + v.attrs.id, inp.attrs.source); }
    for (const prim of mesh.children) {
      if (!['triangles', 'polylist', 'polygons', 'trifans', 'tristrips'].includes(prim.name)) continue;
      const inputs = xKids(prim, 'input'), vin = inputs.find((q) => String(q.attrs.semantic).toUpperCase() === 'VERTEX');
      const pos = vin && src.get(vmap.get(vin.attrs.source) || vin.attrs.source);
      if (!pos) continue;
      const step = Math.max(...inputs.map((q) => +q.attrs.offset || 0)) + 1, vo = +vin.attrs.offset || 0, nv = Math.floor((pos.d.length - pos.off) / pos.stride);
      const poly = (ids) => {
        if (ids.some((i) => !Number.isInteger(i) || i < 0 || i >= nv)) return;
        for (let k = 1; k + 1 < ids.length; k++) for (const q of [ids[0], ids[k], ids[k + 1]]) { const o = pos.off + q * pos.stride; tri.push(pos.d[o], pos.d[o + 1], pos.stride > 2 ? pos.d[o + 2] : 0); }
        if (tri.length > L.triangles * 9) fail(`Mesh has too many triangles (limit ${L.triangles}).`);
      };
      const ps = xKids(prim, 'p').map((p) => numsOf(p.text)), vidx = (p) => { const out = []; for (let k = vo; k < p.length; k += step) out.push(p[k]); return out; };
      if (prim.name === 'triangles') for (const p of ps) { const v = vidx(p); for (let k = 0; k + 2 < v.length; k += 3) poly([v[k], v[k + 1], v[k + 2]]); }
      else if (prim.name === 'polylist') { const vc = numsOf((xKid(prim, 'vcount') || { text: '' }).text), v = vidx(ps[0] || []); let q = 0; for (const n of vc) { if (!(n > 0) || q + n > v.length) break; poly(v.slice(q, q + n)); q += n; } }
      else if (prim.name === 'tristrips') for (const p of ps) { const v = vidx(p); for (let k = 0; k + 2 < v.length; k++) poly(k % 2 ? [v[k + 1], v[k], v[k + 2]] : [v[k], v[k + 1], v[k + 2]]); }
      else for (const p of ps) poly(vidx(p));
    }
  }
  const unitNode = xAll(doc, 'unit')[0], up = xAll(doc, 'up_axis')[0], warnings = [];
  if (xAll(doc, 'node').some((n) => n.children.some((c) => ['matrix', 'translate', 'rotate', 'scale'].includes(c.name)))) warnings.push('Scene-node transforms are present but not applied; geometry is shown in its local coordinates.');
  return meshGeom(tri, { warnings, stats: { unitMetres: unitNode ? +unitNode.attrs.meter || 1 : 1, upAxis: up ? up.text.trim() : 'Y_UP' } });
}

// glTF 2.0 (JSON with embedded or companion buffers, or binary GLB)
const m4mul = (a, b) => { const o = new Array(16); for (let c = 0; c < 4; c++) for (let r = 0; r < 4; r++) o[c * 4 + r] = a[r] * b[c * 4] + a[4 + r] * b[c * 4 + 1] + a[8 + r] * b[c * 4 + 2] + a[12 + r] * b[c * 4 + 3]; return o; };
const M4I = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1];
async function readGLTF(ctx, pre) {
  const u8 = await ctx.bytes();
  let json = pre, bin = null;
  if (!json && u8.length >= 12 && u8[0] === 0x67 && u8[1] === 0x6c && u8[2] === 0x54 && u8[3] === 0x46) {
    const dv = new DataView(u8.buffer, u8.byteOffset, u8.byteLength);
    if (dv.getUint32(4, true) !== 2) fail('Only glTF version 2.0 binary files are supported.');
    for (let p = 12; p + 8 <= u8.length;) {
      const len = dv.getUint32(p, true), type = dv.getUint32(p + 4, true);
      if (p + 8 + len > u8.length) fail('The GLB file is truncated.');
      if (type === 0x4e4f534a) json = safeJSON(utf8.decode(u8.subarray(p + 8, p + 8 + len))); else if (type === 0x004e4942 && !bin) bin = u8.subarray(p + 8, p + 8 + len);
      p += 8 + len;
    }
  } else if (!json) json = safeJSON(utf8.decode(u8));
  if (!json || typeof json !== 'object' || !json.asset || !Array.isArray(json.meshes)) fail('Not a glTF 2.0 file (no asset / meshes).');
  if ((json.extensionsRequired || []).some((e) => /draco|meshopt/i.test(e))) fail(fmtNamed('glTF 2.0').convert);
  const buffers = [];
  for (const [i, b] of (json.buffers || []).entries()) {
    if (!b || typeof b.uri !== 'string') buffers.push(i === 0 ? bin : null);
    else if (/^data:/i.test(b.uri)) buffers.push(b64(b.uri.slice(b.uri.indexOf(',') + 1)));
    else { const c = ctx.companion((n) => baseName(n) === baseName(decodeURIComponentSafe(b.uri))); buffers.push(c ? new Uint8Array(await c.arrayBuffer()) : null); }
  }
  const CT = { 5120: 'i1', 5121: 'u1', 5122: 'i2', 5123: 'u2', 5125: 'u4', 5126: 'f4' }, NC = { SCALAR: 1, VEC2: 2, VEC3: 3, VEC4: 4 }, cache = new Map();
  const accessor = (i) => {
    if (cache.has(i)) return cache.get(i);
    const a = (json.accessors || [])[i], bv = a && (json.bufferViews || [])[a.bufferView], buf = bv && buffers[bv.buffer], ct = a && CT[a.componentType], nc = a && NC[a.type];
    if (!a || !bv || !ct || !nc) fail('The glTF file has an accessor that cannot be read.');
    if (!buf) fail('The glTF file refers to an external buffer; supply the .bin file as a companion or export a self-contained .glb.');
    const sz = DT[ct][1], stride = bv.byteStride || nc * sz, base = (bv.byteOffset || 0) + (a.byteOffset || 0), count = a.count;
    if (!Number.isInteger(count) || count < 0 || !Number.isInteger(base) || base < 0 || (count && base + (count - 1) * stride + nc * sz > buf.length)) fail('A glTF accessor points outside its buffer (file truncated?).');
    const dv = new DataView(buf.buffer, buf.byteOffset, buf.byteLength), out = new Float64Array(count * nc);
    const get = { i1: (o) => dv.getInt8(o), u1: (o) => dv.getUint8(o), i2: (o) => dv.getInt16(o, true), u2: (o) => dv.getUint16(o, true), u4: (o) => dv.getUint32(o, true), f4: (o) => dv.getFloat32(o, true) }[ct];
    for (let k = 0; k < count; k++) for (let c = 0; c < nc; c++) out[k * nc + c] = get(base + k * stride + c * sz);
    cache.set(i, out);
    return out;
  };
  const tri = [], warnings = [];
  let skipped = 0;
  const emitMesh = (mi, M) => {
    const mesh = json.meshes[mi];
    for (const prim of (mesh && mesh.primitives) || []) {
      const mode = prim.mode ?? 4;
      if (![4, 5, 6].includes(mode) || !prim.attributes || prim.attributes.POSITION === undefined) { skipped++; continue; }
      const pos = accessor(prim.attributes.POSITION), nv = pos.length / 3, idx = prim.indices !== undefined ? accessor(prim.indices) : null, n = idx ? idx.length : nv, at = (k) => (idx ? idx[k] : k);
      const put = (a, b, c) => {
        if (a >= nv || b >= nv || c >= nv) fail('A glTF index refers to a missing vertex.');
        for (const q of [a, b, c]) { const x = pos[3 * q], y = pos[3 * q + 1], z = pos[3 * q + 2]; tri.push(M[0] * x + M[4] * y + M[8] * z + M[12], M[1] * x + M[5] * y + M[9] * z + M[13], M[2] * x + M[6] * y + M[10] * z + M[14]); }
      };
      if (tri.length / 9 + n / 3 > L.triangles) fail(`Mesh has too many triangles (limit ${L.triangles}).`);
      if (mode === 4) for (let k = 0; k + 2 < n; k += 3) put(at(k), at(k + 1), at(k + 2));
      else if (mode === 5) for (let k = 0; k + 2 < n; k++) if (k % 2) put(at(k + 1), at(k), at(k + 2)); else put(at(k), at(k + 1), at(k + 2));
      else for (let k = 1; k + 1 < n; k++) put(at(0), at(k), at(k + 1));
    }
  };
  const local = (nd) => {
    if (Array.isArray(nd.matrix) && nd.matrix.length === 16) return nd.matrix.map(Number);
    const t = nd.translation || [0, 0, 0], [x, y, z, w] = nd.rotation || [0, 0, 0, 1], s = nd.scale || [1, 1, 1];
    return [(1 - 2 * (y * y + z * z)) * s[0], 2 * (x * y + z * w) * s[0], 2 * (x * z - y * w) * s[0], 0, 2 * (x * y - z * w) * s[1], (1 - 2 * (x * x + z * z)) * s[1], 2 * (y * z + x * w) * s[1], 0, 2 * (x * z + y * w) * s[2], 2 * (y * z - x * w) * s[2], (1 - 2 * (x * x + y * y)) * s[2], 0, t[0], t[1], t[2], 1];
  };
  const nodes = Array.isArray(json.nodes) ? json.nodes : [];
  let visits = 0;
  const walk = (ni, M, depth) => {
    const nd = nodes[ni];
    if (!nd || typeof nd !== 'object' || depth > 64 || ++visits > 200000) return;
    const W = m4mul(M, local(nd));
    if (Number.isInteger(nd.mesh)) emitMesh(nd.mesh, W);
    for (const c of Array.isArray(nd.children) ? nd.children : []) walk(c, W, depth + 1);
  };
  if (nodes.some((nd) => nd && Number.isInteger(nd.mesh))) {
    const sc = (json.scenes || [])[json.scene ?? 0], kids = new Set(nodes.flatMap((nd) => (nd && Array.isArray(nd.children) ? nd.children : [])));
    for (const r of sc && Array.isArray(sc.nodes) ? sc.nodes : nodes.map((_, k) => k).filter((k) => !kids.has(k))) walk(r, M4I, 0);
  }
  if (!tri.length) json.meshes.forEach((_, k) => emitMesh(k, M4I));
  if (skipped) warnings.push(`${skipped} point or line primitives were ignored.`);
  return meshGeom(tri, { warnings, stats: { meshes: json.meshes.length, nodes: nodes.length } });
}
const decodeURIComponentSafe = (s) => { try { return decodeURIComponent(s); } catch { return s; } };

// VTK cell type -> [cell kind, corner count]; higher-order cells keep their corner nodes only
const VTK_CELL = { 3: ['line', 2], 21: ['line', 2], 5: ['tri', 3], 22: ['tri', 3], 9: ['quad', 4], 23: ['quad', 4], 28: ['quad', 4], 10: ['tet', 4], 24: ['tet', 4], 12: ['hex', 8], 25: ['hex', 8], 29: ['hex', 8], 13: ['wedge', 6], 26: ['wedge', 6], 14: ['pyramid', 5], 27: ['pyramid', 5] };
function vtkCell(store, type, ids) {
  if (type === 7) store.add(ids.length === 3 ? 'tri' : ids.length === 4 ? 'quad' : 'poly', ids);
  else if (type === 8 && ids.length >= 4) store.add('quad', [ids[0], ids[1], ids[3], ids[2]]);
  else if (type === 11 && ids.length >= 8) store.add('hex', [ids[0], ids[1], ids[3], ids[2], ids[4], ids[5], ids[7], ids[6]]);
  else if (type === 6) for (let k = 0; k + 2 < ids.length; k++) store.add('tri', k % 2 ? [ids[k + 1], ids[k], ids[k + 2]] : [ids[k], ids[k + 1], ids[k + 2]]);
  else { const c = VTK_CELL[type]; if (!c || ids.length < c[1]) return false; store.add(c[0], ids.length === c[1] ? ids : ids.slice(0, c[1])); }
  return true;
}
/** Regular scalar image -> DEM grid (2-D, real-valued) or solid/pore voxels. */
function imageGeom(vals, nx, ny, nz, spacing, origin, opts, extra = {}) {
  const n = nx * ny * nz;
  if (nz === 1 && opts.threshold === undefined && n > 0 && n <= vals.length) {
    let real = false;
    for (let i = 0; i < n && !real; i++) if (vals[i] !== Math.floor(vals[i])) real = true;
    if (real) return gridGeom(vals, nx, ny, (i) => origin[0] + i * spacing[0], (j) => origin[1] + j * spacing[1], extra);
  }
  const g = voxelsFromValues(vals, nx, ny, nz, spacing, opts, extra);
  if (origin.some((v) => v)) { g.voxels.origin = origin.slice(0, 3); g.bbox = { min: origin.slice(0, 3), max: [0, 1, 2].map((k) => origin[k] + [nx, ny, nz][k] * g.voxels.spacing[k]) }; }
  return g;
}
const rectilinear = (cx, cy, cz) => {
  const ni = cx.length, nj = cy.length, nk = cz.length;
  if (ni * nj * nk > 3 * L.cells) fail('The rectilinear grid is too large.');
  const xyz = new Float64Array(3 * ni * nj * nk);
  let o = 0;
  for (let k = 0; k < nk; k++) for (let j = 0; j < nj; j++) for (let i = 0; i < ni; i++) { xyz[o++] = cx[i]; xyz[o++] = cy[j]; xyz[o++] = cz[k]; }
  return { xyz, ni, nj, nk };
};

async function readVTK(ctx) {
  const u8 = await ctx.bytes(), n = u8.length;
  let p = 0;
  const line = () => { let e = p; while (e < n && u8[e] !== 10) e++; const s = latin1.decode(u8.subarray(p, e)); p = e + 1; return s; };
  if (!/^# vtk DataFile/i.test(line())) fail('Not a legacy VTK file.');
  line();
  const binary = /BINARY/i.test(line());
  const word = () => { while (p < n && u8[p] <= 32) p++; let e = p; while (e < n && u8[e] > 32) e++; const s = latin1.decode(u8.subarray(p, e)); p = e; return s; };
  const int = () => { const v = +word(); if (!Number.isInteger(v) || v < 0) fail('The VTK file has an invalid count (truncated or corrupt).'); return v; };
  const arr = (count, type) => {
    const dt = dtypeOf(type);
    if (!dt) fail(`The VTK file uses an unsupported data type "${type}".`);
    if (!Number.isInteger(count) || count < 0 || count > n) fail('The VTK file declares more data than it contains.');
    if (binary) {
      if (String(type).toLowerCase() === 'bit') fail('Bit-packed binary VTK arrays are not supported.');
      while (p < n && (u8[p] === 32 || u8[p] === 9)) p++;
      if (u8[p] === 13) p++;
      if (u8[p] === 10) p++;
      const a = typed(u8, p, count, dt, false);
      p += count * DT[dt][1];
      return a;
    }
    const out = new Float64Array(count);
    for (let i = 0; i < count; i++) { const w = word(); if (!w) fail('The VTK file is truncated.'); out[i] = +w; }
    return out;
  };
  const cells = () => {
    const a = int(), b = int(), save = p, out = [];
    if (word().toUpperCase() === 'OFFSETS') {
      const off = arr(a, word());
      if (word().toUpperCase() !== 'CONNECTIVITY') fail('The VTK cell section is malformed.');
      const con = arr(b, word());
      for (let i = 0; i + 1 < off.length; i++) out.push(Array.from(con.subarray(off[i], off[i + 1])));
      return out;
    }
    p = save;
    const d = arr(b, 'int');
    for (let q = 0, c = 0; c < a && q < d.length; c++) { const k = d[q]; if (!(k >= 0) || q + 1 + k > d.length) fail('The VTK cell section is malformed.'); out.push(Array.from(d.subarray(q + 1, q + 1 + k))); q += 1 + k; }
    return out;
  };
  const store = cellStore(), cxyz = [null, null, null], warnings = [];
  let dataset = '', pts = null, dims = null, origin = [0, 0, 0], spacing = [1, 1, 1], cellList = null, types = null, dataN = 0, dataKind = '', scal = null, scalKind = '', skipped = 0;
  parse: for (;;) {
    const w = word().toUpperCase();
    switch (w) {
      case '': break parse;
      case 'DATASET': dataset = word().toUpperCase(); break;
      case 'POINTS': { const np = int(); pts = Float64Array.from(arr(np * 3, word())); break; }
      case 'VERTICES': case 'LINES': cells(); break;
      case 'POLYGONS': for (const c of cells()) store.add(c.length === 3 ? 'tri' : c.length === 4 ? 'quad' : 'poly', c); break;
      case 'TRIANGLE_STRIPS': for (const c of cells()) vtkCell(store, 6, c); break;
      case 'CELLS': cellList = cells(); break;
      case 'CELL_TYPES': types = arr(int(), 'int'); break;
      case 'DIMENSIONS': dims = [int(), int(), int()]; break;
      case 'ORIGIN': origin = [+word(), +word(), +word()]; break;
      case 'SPACING': case 'ASPECT_RATIO': spacing = [+word(), +word(), +word()]; break;
      case 'X_COORDINATES': case 'Y_COORDINATES': case 'Z_COORDINATES': { const k = int(); cxyz['XYZ'.indexOf(w[0])] = Float64Array.from(arr(k, word())); break; }
      case 'POINT_DATA': case 'CELL_DATA': dataN = int(); dataKind = w; break;
      case 'SCALARS': {
        word();
        const type = word();
        let save = p, t = word(), nc = 1;
        if (/^\d+$/.test(t)) { nc = +t; save = p; t = word(); }
        if (t.toUpperCase() === 'LOOKUP_TABLE') word(); else p = save;
        const a = arr(dataN * nc, type);
        if (!scal && nc === 1) { scal = a; scalKind = dataKind; }
        break;
      }
      case 'VECTORS': case 'NORMALS': { word(); arr(dataN * 3, word()); break; }
      case 'METADATA': { line(); while (p < n && line().trim() !== ''); break; }
      case 'FIELD': {
        word();
        const na = int();
        let known = true;
        for (let q = 0; q < na && known; q++) { word(); const nc = int(), nt = int(), ty = word(); if (dtypeOf(ty)) arr(nc * nt, ty); else known = false; }
        if (!known) break parse;
        break;
      }
      default: break parse;
    }
  }
  const title = 'VTK ' + dataset.toLowerCase().replace(/_/g, ' ');
  if (dataset === 'STRUCTURED_POINTS') {
    if (!dims || !scal) fail('The VTK image has no scalar data to threshold.');
    const d = scalKind === 'CELL_DATA' ? dims.map((k) => Math.max(1, k - 1)) : dims;
    return imageGeom(scal, d[0], d[1], d[2], spacing, origin, ctx.opts, { warnings, stats: { dataset: title } });
  }
  if (dataset === 'STRUCTURED_GRID') { if (!dims || !pts) fail('The VTK structured grid has no points.'); return structuredGeom([{ xyz: pts, ni: dims[0], nj: dims[1], nk: dims[2] }], { warnings }); }
  if (dataset === 'RECTILINEAR_GRID') { if (!cxyz.every(Boolean)) fail('The VTK rectilinear grid has no coordinates.'); return structuredGeom([rectilinear(cxyz[0], cxyz[1], cxyz[2])], { warnings }); }
  if (!pts) fail('The VTK file has no POINTS section.');
  if (cellList) cellList.forEach((c, k) => { if (!vtkCell(store, types ? types[k] : 7, c)) skipped++; });
  if (skipped) warnings.push(`${skipped} cells of unsupported type (vertices, polyhedra …) were ignored.`);
  return cellMesh(pts, store, { warnings });
}

async function readVTKXML(ctx, pre) {
  const u8 = await ctx.bytes();
  let text = latin1.decode(u8), appRaw = null, appText = null;
  const ai = text.indexOf('<AppendedData');
  if (ai >= 0) {
    const gt = text.indexOf('>', ai), us = gt < 0 ? -1 : text.indexOf('_', gt);
    if (us >= 0) { if (/encoding\s*=\s*["']raw/i.test(text.slice(ai, gt))) appRaw = u8.subarray(us + 1); else { const e = text.indexOf('<', us); appText = text.slice(us + 1, e < 0 ? undefined : e).replace(/\s+/g, ''); } }
    text = text.slice(0, ai) + '</VTKFile>';
  }
  const root = xAll(ai >= 0 || !pre ? parseXML(text) : pre, 'vtkfile')[0];
  if (!root) fail('Not a VTK XML file (no <VTKFile> element).');
  const kind = String(root.attrs.type || '').toLowerCase(), little = String(root.attrs.byte_order || 'LittleEndian').toLowerCase() !== 'bigendian', ht = String(root.attrs.header_type || 'UInt32').toLowerCase() === 'uint64' ? 'u8' : 'u4', hs = DT[ht][1], comp = root.attrs.compressor || '';
  if (comp && !/zlib/i.test(comp)) fail('The VTK XML data uses an unsupported compressor (LZ4 / LZMA); re-save it with zlib or without compression.');
  const b64Reader = (s, pos) => {
    let q = new Uint8Array(0);
    return (nb) => {
      if (!(nb >= 0) || nb > L.inflated) fail('A VTK data array is too large.');
      const out = new Uint8Array(nb);
      let o = 0;
      while (o < nb) {
        if (!q.length) {
          let chunk = s.substr(pos, Math.ceil((nb - o) / 3) * 4);
          const eqi = chunk.indexOf('=');
          if (eqi >= 0) chunk = chunk.slice(0, (Math.floor(eqi / 4) + 1) * 4);
          if (chunk.length < 4) fail('A VTK data array is truncated.');
          pos += chunk.length;
          q = b64(chunk);
          if (!q.length) fail('A VTK data array is truncated.');
        }
        const k = Math.min(q.length, nb - o);
        out.set(q.subarray(0, k), o); o += k; q = q.subarray(k);
      }
      return out;
    };
  };
  const rawReader = (b, pos) => (nb) => { if (!(nb >= 0) || pos + nb > b.length) fail('A VTK data array is truncated.'); const o = b.subarray(pos, pos + nb); pos += nb; return o; };
  const da = async (node) => {
    if (!node) return null;
    const dt = dtypeOf(node.attrs.type), fmt = String(node.attrs.format || 'ascii').toLowerCase();
    if (!dt) fail(`The VTK file uses an unsupported data type "${node.attrs.type}".`);
    if (fmt === 'ascii') return numsOf(node.text);
    const off = +node.attrs.offset || 0;
    let src;
    if (fmt === 'binary') src = b64Reader(node.text.replace(/\s+/g, ''), 0);
    else if (appRaw) src = rawReader(appRaw, off);
    else if (appText !== null) src = b64Reader(appText, off);
    else fail('The VTK file refers to appended data that is missing.');
    const hd = (k) => typed(src(k * hs), 0, k, ht, little);
    let bytes;
    if (!comp) bytes = src(hd(1)[0]);
    else {
      const [nb] = hd(3);
      if (!(nb >= 0 && nb <= 1e6)) fail('A compressed VTK array has an invalid header.');
      const cs = nb ? hd(nb) : [], parts = [];
      let tot = 0;
      for (const c of cs) { const part = await inflate(src(c), 'deflate', L.inflated - tot); tot += part.length; parts.push(part); }
      bytes = new Uint8Array(tot);
      let o = 0;
      for (const part of parts) { bytes.set(part, o); o += part.length; }
    }
    return typed(bytes, 0, Math.floor(bytes.length / DT[dt][1]), dt, little);
  };
  const named = (parent, name) => parent && parent.children.find((c) => c.name === 'dataarray' && String(c.attrs.name || '').toLowerCase() === name);
  const first = (parent) => parent && parent.children.find((c) => c.name === 'dataarray');
  const body = root.children.find((c) => c.name === kind), pieces = body ? xKids(body, 'piece') : [];
  if (!body || !pieces.length) fail('The VTK XML file has no data pieces (parallel index files such as .pvtu are not followed).');
  const ext = (s) => { const e = numsOf(s || ''); return e.length >= 6 ? [e[1] - e[0] + 1, e[3] - e[2] + 1, e[5] - e[4] + 1] : null; }, warnings = [];
  if (kind === 'imagedata') {
    const pc = pieces[0], d = ext(pc.attrs.extent || body.attrs.wholeextent), pd = first(xKid(pc, 'pointdata')), cd = first(xKid(pc, 'celldata'));
    if (!d || (!pd && !cd)) fail('The VTK image has no extent or no data array.');
    const dd = pd ? d : d.map((k) => Math.max(1, k - 1)), sp = numsOf(body.attrs.spacing || '1 1 1'), or = numsOf(body.attrs.origin || '0 0 0');
    return imageGeom(await da(pd || cd), dd[0], dd[1], dd[2], Array.from(sp), Array.from(or), ctx.opts, { warnings });
  }
  if (kind === 'structuredgrid' || kind === 'rectilineargrid') {
    const blocks = [];
    for (const pc of pieces) {
      const d = ext(pc.attrs.extent || body.attrs.wholeextent);
      if (!d) fail('The VTK grid piece has no extent.');
      if (kind === 'structuredgrid') blocks.push({ xyz: Float64Array.from(await da(first(xKid(pc, 'points'))) || []), ni: d[0], nj: d[1], nk: d[2] });
      else { const cs = xKid(pc, 'coordinates'), arrs = []; for (const c of cs ? xKids(cs, 'dataarray') : []) arrs.push(await da(c)); if (arrs.length < 3) fail('The VTK rectilinear grid has no coordinates.'); blocks.push(rectilinear(arrs[0], arrs[1], arrs[2])); }
    }
    return structuredGeom(blocks, { warnings });
  }
  if (kind !== 'polydata' && kind !== 'unstructuredgrid') fail(`VTK XML data sets of type "${root.attrs.type}" are not supported.`);
  const store = cellStore(), chunks = [];
  let base = 0, skipped = 0;
  for (const pc of pieces) {
    const pts = await da(first(xKid(pc, 'points')));
    if (!pts) continue;
    chunks.push(pts);
    for (const sec of kind === 'polydata' ? ['polys', 'strips'] : ['cells']) {
      const s = xKid(pc, sec), con = await da(named(s, 'connectivity')), off = await da(named(s, 'offsets')), types = sec === 'cells' ? await da(named(s, 'types')) : null;
      if (!con || !off) continue;
      for (let c = 0, a = 0; c < off.length; c++) {
        const b = off[c];
        if (!(b >= a && b <= con.length)) fail('The VTK connectivity offsets are invalid.');
        const ids = new Array(b - a);
        for (let k = a; k < b; k++) ids[k - a] = con[k] + base;
        if (!vtkCell(store, sec === 'strips' ? 6 : types ? types[c] : 7, ids)) skipped++;
        a = b;
      }
    }
    base += Math.floor(pts.length / 3);
  }
  const xyz = new Float64Array(base * 3);
  let o = 0;
  for (const c of chunks) { xyz.set(c.subarray(0, Math.floor(c.length / 3) * 3), o); o += Math.floor(c.length / 3) * 3; }
  if (skipped) warnings.push(`${skipped} cells of unsupported type (vertices, polyhedra …) were ignored.`);
  return cellMesh(xyz, store, { warnings });
}

async function readGTS(ctx) {
  const ls = splitLines(await ctx.text()).map((l) => l.replace(/[#!].*/, '').trim()).filter(Boolean), h = (ls[0] || '').split(/\s+/).map(Number), [nv, ne, nf] = h;
  if (![nv, ne, nf].every((k) => Number.isInteger(k) && k >= 0) || 1 + nv + ne + nf > ls.length) fail('Not a GTS file, or the file is truncated.');
  const xyz = new Float64Array(nv * 3), ed = new Int32Array(ne * 2), tri = [];
  for (let k = 0; k < nv; k++) { const t = ls[1 + k].split(/\s+/); xyz[3 * k] = +t[0]; xyz[3 * k + 1] = +t[1]; xyz[3 * k + 2] = +t[2]; }
  for (let k = 0; k < ne; k++) { const t = ls[1 + nv + k].split(/\s+/); ed[2 * k] = +t[0] - 1; ed[2 * k + 1] = +t[1] - 1; }
  for (let k = 0; k < nf; k++) {
    const t = ls[1 + nv + ne + k].split(/\s+/), e1 = +t[0] - 1, e2 = +t[1] - 1;
    if (!(e1 >= 0 && e1 < ne && e2 >= 0 && e2 < ne)) fail('A GTS face refers to a missing edge.');
    const a = ed[2 * e1], b = ed[2 * e1 + 1], c = ed[2 * e2], d = ed[2 * e2 + 1];
    if (!pushFan(tri, xyz, b === c ? [a, b, d] : b === d ? [a, b, c] : a === c ? [b, a, d] : [b, a, c])) fail('A GTS edge refers to a missing vertex.');
  }
  return meshGeom(tri, { stats: { vertices: nv, edges: ne, faces: nf } });
}

async function readBYU(ctx) {
  const v = numsOf(await ctx.text()), [np, nj, npoly, nconn] = v;
  if (![np, nj, npoly, nconn].every((k) => Number.isInteger(k) && k >= 0) || v.length < 4 + 2 * np + 3 * nj + nconn) fail('Not a BYU surface file, or the file is truncated.');
  const xyz = v.subarray(4 + 2 * np, 4 + 2 * np + 3 * nj), con = v.subarray(4 + 2 * np + 3 * nj, 4 + 2 * np + 3 * nj + nconn), tri = [];
  let cur = [];
  for (const c of con) { cur.push(Math.abs(c) - 1); if (c < 0) { if (!pushFan(tri, xyz, cur)) fail('A BYU polygon refers to a missing vertex.'); cur = []; } }
  return meshGeom(tri, { stats: { parts: np, vertices: nj, polygons: npoly } });
}

/** Generic "counts, vertex rows, triangle rows" ASCII surfaces (.tri / .fac / .surf). */
async function readNodeFace(ctx) {
  const rows = splitLines(await ctx.text()).map((l) => l.replace(/[#!%].*/, '').trim()).filter(Boolean).map((l) => l.split(/[\s,]+/).map(Number));
  const bad = () => fail('This .tri/.fac/.surf dialect is not recognised. ' + fmtNamed('Generic ASCII node-face file').convert);
  if (!rows.length || rows.some((r) => r.some((x) => !Number.isFinite(x)))) bad();
  const ints = (r) => r.every(Number.isInteger);
  let nv, nt, i = 1;
  if (rows[0].length === 2 && ints(rows[0])) [nv, nt] = rows[0];
  else if (rows[0].length === 1 && ints(rows[0])) { nv = rows[0][0]; const r = rows[1 + nv]; if (!r || r.length !== 1 || !ints(r)) bad(); nt = r[0]; }
  else bad();
  if (!(nv >= 3 && nt >= 1)) bad();
  const two = rows[0].length === 2, t0 = two ? 1 + nv : 2 + nv;
  if (rows.length < t0 + nt) bad();
  const xyz = new Float64Array(nv * 3), w = rows[1].length, vo = w === 4 ? 1 : 0;
  if (w !== 3 && w !== 4) bad();
  for (let k = 0; k < nv; k++) { const r = rows[i + k]; if (r.length !== w) bad(); xyz[3 * k] = r[vo]; xyz[3 * k + 1] = r[vo + 1]; xyz[3 * k + 2] = r[vo + 2]; }
  const fw = rows[t0].length, fo = fw === 4 && rows[t0][0] === 1 && rows[t0 + nt - 1][0] === nt ? 1 : 0, tris = [];
  if (fw < 3 || fw > 5) bad();
  let lo = Infinity;
  for (let k = 0; k < nt; k++) { const r = rows[t0 + k]; if (r.length < fo + 3 || !ints(r)) bad(); const t = [r[fo], r[fo + 1], r[fo + 2]]; lo = Math.min(lo, ...t); tris.push(t); }
  const tri = [];
  for (const t of tris) if (!pushFan(tri, xyz, t.map((q) => q - (lo >= 1 ? 1 : 0)))) bad();
  return meshGeom(tri, { stats: { vertices: nv, faces: nt } });
}

// ---- Computational meshes ------------------------------------------------------------------------------
const GMSH_CELL = { 1: ['line', 2], 8: ['line', 2], 2: ['tri', 3], 9: ['tri', 3], 20: ['tri', 3], 21: ['tri', 3], 3: ['quad', 4], 10: ['quad', 4], 16: ['quad', 4], 4: ['tet', 4], 11: ['tet', 4], 29: ['tet', 4], 5: ['hex', 8], 12: ['hex', 8], 17: ['hex', 8], 6: ['wedge', 6], 13: ['wedge', 6], 18: ['wedge', 6], 7: ['pyramid', 5], 14: ['pyramid', 5], 19: ['pyramid', 5] };
const intsOf = (l) => l.trim().split(/\s+/).map(Number);
async function readGmsh(ctx) {
  const text = await ctx.text();
  const sect = (name) => { const a = text.indexOf('$' + name); if (a < 0) return null; const s = text.indexOf('\n', a) + 1, e = text.indexOf('$End' + name, s); return splitLines(text.slice(s, e < 0 ? undefined : e)).filter((l) => l.trim()); };
  const mf = sect('MeshFormat');
  if (!mf) fail(text.includes('$NOD') ? 'Gmsh version 1 files are not supported; save the mesh as version 2.2 or 4.1 ASCII.' : 'Not a Gmsh .msh file ($MeshFormat is missing).');
  const [ver, ftype] = intsOf(mf[0] || '');
  if (ftype === 1) fail('Binary Gmsh files are not supported. ' + fmtNamed('Gmsh MSH').convert);
  const nl = sect('Nodes'), el = sect('Elements');
  if (!nl || !el || !nl.length || !el.length) fail('The Gmsh file has no $Nodes or $Elements section.');
  const idmap = new Map(), store = cellStore();
  let xyz, skipped = 0;
  const addEl = (type, tags) => {
    const c = GMSH_CELL[type];
    if (!c) { skipped++; return; }
    const ids = [];
    for (let k = 0; k < c[1]; k++) { const q = idmap.get(tags[k]); if (q === undefined) fail('A Gmsh element refers to a missing node.'); ids.push(q); }
    store.add(c[0], ids);
  };
  const chk = (n, have) => { if (!Number.isInteger(n) || n < 0 || n > have) fail('The Gmsh file is truncated or its counts are invalid.'); return n; };
  if (ver < 3) {
    const nn = chk(+nl[0], nl.length - 1);
    xyz = new Float64Array(nn * 3);
    for (let k = 0; k < nn; k++) { const t = nl[1 + k].trim().split(/\s+/); idmap.set(+t[0], k); xyz[3 * k] = +t[1]; xyz[3 * k + 1] = +t[2]; xyz[3 * k + 2] = +t[3]; }
    const ne = chk(+el[0], el.length - 1);
    for (let k = 0; k < ne; k++) { const t = intsOf(el[1 + k]); addEl(t[1], t.slice(3 + (t[2] || 0))); }
  } else {
    const v41 = ver >= 4.1, [nb, nn] = intsOf(nl[0]);
    chk(nn, nl.length); chk(nb, nl.length);
    xyz = new Float64Array(nn * 3);
    let i = 1, q = 0;
    for (let b = 0; b < nb; b++) {
      const h = intsOf(nl[i++] || ''), cnt = chk(h[3], nl.length - i);
      if (q + cnt > nn) fail('The Gmsh node blocks hold more nodes than declared.');
      if (v41) { for (let k = 0; k < cnt; k++) idmap.set(+nl[i + k], q + k); i += cnt; for (let k = 0; k < cnt; k++, q++) { const t = nl[i + k].trim().split(/\s+/); xyz[3 * q] = +t[0]; xyz[3 * q + 1] = +t[1]; xyz[3 * q + 2] = +t[2]; } i += cnt; }
      else for (let k = 0; k < cnt; k++, q++) { const t = nl[i++].trim().split(/\s+/); idmap.set(+t[0], q); xyz[3 * q] = +t[1]; xyz[3 * q + 1] = +t[2]; xyz[3 * q + 2] = +t[3]; }
    }
    const [eb] = intsOf(el[0]);
    i = 1;
    for (let b = 0; b < chk(eb, el.length); b++) { const h = intsOf(el[i++] || ''), cnt = chk(h[3], el.length - i); for (let k = 0; k < cnt; k++) { const t = intsOf(el[i++]); addEl(h[2], t.slice(1)); } }
  }
  return cellMesh(xyz, store, { warnings: skipped ? [`${skipped} elements of unsupported type (points, high-order …) were ignored.`] : [], stats: { version: ver } });
}

function signedVolume(t) {
  let v = 0;
  for (let i = 0; i < t.length; i += 9) v += t[i] * (t[i + 4] * t[i + 8] - t[i + 5] * t[i + 7]) - t[i + 1] * (t[i + 3] * t[i + 8] - t[i + 5] * t[i + 6]) + t[i + 2] * (t[i + 3] * t[i + 7] - t[i + 4] * t[i + 6]);
  return v / 6;
}
const flipTriangles = (t) => { for (let i = 0; i < t.length; i += 9) for (let c = 0; c < 3; c++) { const s = t[i + 3 + c]; t[i + 3 + c] = t[i + 6 + c]; t[i + 6 + c] = s; } };

async function readFluent(ctx) {
  const text = await ctx.text(), HX = '([0-9a-fA-F]+)', hex = (s) => parseInt(s, 16);
  if (/\(\s*[23]0(10|12|13)\s*\(/.test(text)) fail('This Fluent file has binary sections. ' + fmtNamed('Fluent mesh (ASCII)').convert);
  const dm = /\(\s*2\s+(\d)\s*\)/.exec(text), dim = dm ? +dm[1] : 3, cm = new RegExp(`\\(\\s*12\\s*\\(\\s*0\\s+${HX}\\s+${HX}`).exec(text);
  const reN = new RegExp(`\\(\\s*10\\s*\\(\\s*${HX}\\s+${HX}\\s+${HX}\\s+${HX}\\s*([0-9a-fA-F]*)\\s*\\)\\s*(\\()?`, 'g'), reF = new RegExp(`\\(\\s*13\\s*\\(\\s*${HX}\\s+${HX}\\s+${HX}\\s+${HX}\\s*([0-9a-fA-F]*)\\s*\\)\\s*(\\()?`, 'g');
  let m, xyz = null, nn = 0;
  const bodyOf = (re) => { const s = re.lastIndex, e = text.indexOf(')', s); if (e < 0) fail('The Fluent mesh file is truncated.'); re.lastIndex = e; return text.slice(s, e); };
  while ((m = reN.exec(text))) {
    const first = hex(m[2]), last = hex(m[3]), nd = m[5] ? hex(m[5]) : dim;
    if (hex(m[1]) === 0 || !m[6]) { nn = Math.max(nn, last); continue; }
    nn = Math.max(nn, last);
    if (nn > text.length) fail('The Fluent mesh declares more nodes than the file can hold.');
    if (!xyz || xyz.length < nn * 3) { const g = new Float64Array(nn * 3); if (xyz) g.set(xyz); xyz = g; }
    const v = numsOf(bodyOf(reN));
    if (!(nd === 2 || nd === 3) || v.length < (last - first + 1) * nd) fail('A Fluent node section is shorter than declared.');
    for (let k = 0; k <= last - first; k++) { xyz[3 * (first - 1 + k)] = v[k * nd]; xyz[3 * (first - 1 + k) + 1] = v[k * nd + 1]; xyz[3 * (first - 1 + k) + 2] = nd === 3 ? v[k * nd + 2] : 0; }
  }
  if (!xyz) fail('No ASCII node section "(10 (…" was found: not a Fluent ASCII mesh.');
  const tri = [], edges = [];
  let nFaces = 0, nBnd = 0;
  while ((m = reF.exec(text))) {
    const first = hex(m[2]), last = hex(m[3]), et = m[5] ? hex(m[5]) : 0;
    if (hex(m[1]) === 0 || !m[6]) continue;
    const tk = bodyOf(reF).match(/\S+/g) || [];
    let q = 0;
    for (let f = first; f <= last; f++) {
      const n = et === 0 || et === 5 ? hex(tk[q++]) : et;
      if (!(n >= 2 && n <= 256) || q + n + 2 > tk.length) fail('A Fluent face section is shorter than declared.');
      const ids = [];
      for (let k = 0; k < n; k++) ids.push(hex(tk[q++]) - 1);
      const c0 = hex(tk[q++]), c1 = hex(tk[q++]);
      nFaces++;
      if (c0 !== 0 && c1 !== 0) continue;
      nBnd++;
      if (c0 === 0) ids.reverse();
      if (n === 2) edges.push(ids); else if (!pushFan(tri, xyz, ids)) fail('A Fluent face refers to a missing node.');
    }
  }
  const stats = { nodes: nn, cells: cm ? hex(cm[2]) : undefined, faces: nFaces, boundaryFaces: nBnd, cellTypes: {}, dimension: dim };
  if (edges.length && !tri.length) {
    for (const e of edges) for (const i of e) if (!(i >= 0 && i < nn)) fail('A Fluent face refers to a missing node.');
    return polyGeom(chainEdges(edges).map((c) => ({ x: c.ids.map((i) => xyz[3 * i]), y: c.ids.map((i) => xyz[3 * i + 1]), closed: c.closed })), { stats: { ...stats, dimension: 2 } });
  }
  if (!tri.length) fail('The Fluent mesh has no boundary faces.');
  if (signedVolume(tri) < 0) flipTriangles(tri);
  return meshGeom(tri, { stats });
}

async function readSU2(ctx) {
  const ls = splitLines(await ctx.text()), store = cellStore();
  let ndim = 3, xyz = null, skipped = 0;
  const need = (i, n) => { if (!Number.isInteger(n) || n < 0 || i + n > ls.length - 1) fail('The SU2 file is truncated or its counts are invalid.'); };
  for (let i = 0; i < ls.length; i++) {
    const m = /^\s*([A-Za-z_]+)\s*=\s*(-?\d+)/.exec(ls[i].replace(/%.*/, ''));
    if (!m) continue;
    const key = m[1].toUpperCase(), val = +m[2];
    if (key === 'NDIME') ndim = val;
    else if (key === 'NELEM') { need(i, val); for (let k = 0; k < val; k++) { const t = intsOf(ls[++i]), c = VTK_CELL[t[0]]; if (!c || t.length < 1 + c[1]) skipped++; else store.add(c[0], t.slice(1, 1 + c[1])); } }
    else if (key === 'NPOIN') { need(i, val); xyz = new Float64Array(val * 3); for (let k = 0; k < val; k++) { const t = ls[++i].trim().split(/\s+/); xyz[3 * k] = +t[0]; xyz[3 * k + 1] = +t[1]; xyz[3 * k + 2] = ndim === 3 ? +t[2] : 0; } }
    else if (key === 'MARKER_ELEMS') { need(i, val); i += val; }
  }
  if (!xyz) fail('Not an SU2 mesh (NPOIN section is missing).');
  return cellMesh(xyz, store, { flat: ndim === 2 ? true : undefined, warnings: skipped ? [`${skipped} elements of unsupported type were ignored.`] : [] });
}

async function readUNV(ctx) {
  const ls = splitLines(await ctx.text()), sep = (l) => /^\s*-1\s*$/.test(l), idmap = new Map(), coords = [], raw = [];
  let i = 0;
  while (i < ls.length) {
    if (!sep(ls[i])) { i++; continue; }
    const id = parseInt(ls[i + 1], 10);
    i += 2;
    if (id === 2411 || id === 781) while (i + 1 < ls.length && !sep(ls[i])) {
      const c = ls[i + 1].replace(/[dD]/g, 'e').trim().split(/\s+/).map(Number);
      idmap.set(parseInt(ls[i], 10), coords.length / 3); coords.push(c[0], c[1], c[2] || 0); i += 2;
    } else if (id === 2412) while (i < ls.length && !sep(ls[i])) {
      const r = intsOf(ls[i++]), fe = r[0 + 1], nn = r[5], nodes = [];
      if (!(nn > 0 && nn <= 64)) fail('The UNV element data set is malformed.');
      if ([11, 21, 22, 23, 24].includes(fe)) i++;
      while (nodes.length < nn && i < ls.length && !sep(ls[i])) nodes.push(...intsOf(ls[i++]));
      raw.push([fe, nodes]);
    } else while (i < ls.length && !sep(ls[i])) i++;
    i++;
  }
  if (!coords.length) fail('No node data set (2411) was found in the UNV file.');
  const store = cellStore(), pick = (n, ix) => ix.map((k) => n[k]);
  let skipped = 0;
  for (const [fe, tags] of raw) {
    const n = tags.map((t) => idmap.get(t)), k = n.length;
    if (n.some((q) => q === undefined)) fail('A UNV element refers to a missing node.');
    if ([11, 21, 22, 23, 24].includes(fe)) store.add('line', [n[0], n[k - 1]]);
    else if (fe < 100) { if (k === 3) store.add('tri', n); else if (k === 6) store.add('tri', pick(n, [0, 2, 4])); else if (k === 4) store.add('quad', n); else if (k === 8) store.add('quad', pick(n, [0, 2, 4, 6])); else skipped++; }
    else if (k === 4) store.add('tet', n); else if (k === 10) store.add('tet', pick(n, [0, 2, 4, 9]));
    else if (k === 6) store.add('wedge', n); else if (k === 15) store.add('wedge', pick(n, [0, 2, 4, 9, 11, 13]));
    else if (k === 8) store.add('hex', n); else if (k === 20) store.add('hex', pick(n, [0, 2, 4, 6, 12, 14, 16, 18]));
    else if (k === 5) store.add('pyramid', n); else skipped++;
  }
  return cellMesh(Float64Array.from(coords), store, { warnings: skipped ? [`${skipped} elements of unsupported type were ignored.`] : [] });
}

const NAS_CELL = { CTRIA3: ['tri', 3], CTRIA6: ['tri', 3], CTRIAR: ['tri', 3], CQUAD4: ['quad', 4], CQUAD8: ['quad', 4], CQUADR: ['quad', 4], CTETRA: ['tet', 4], CHEXA: ['hex', 8], CPENTA: ['wedge', 6], CPYRAM: ['pyramid', 5], CBAR: ['line', 2], CBEAM: ['line', 2], CROD: ['line', 2] };
async function readNastran(ctx) {
  const ls = splitLines(await ctx.text());
  let start = ls.findIndex((l) => /^\s*BEGIN\s+BULK/i.test(l));
  start = start < 0 ? 0 : start + 1;
  const num = (s) => { let t = String(s ?? '').trim().replace(/[dD]/, 'e'); if (!t) return NaN; if (!/e/i.test(t)) t = t.replace(/([0-9.])([+-]\d+)$/, '$1e$2'); return Number(t); };
  const cards = [];
  let cur = null;
  for (let i = start; i < ls.length; i++) {
    let l = ls[i];
    const ci = l.indexOf('$');
    if (ci >= 0) l = l.slice(0, ci);
    if (!l.trim()) continue;
    if (/^\s*ENDDATA/i.test(l)) break;
    if (l.includes('\t')) { let o = ''; for (const ch of l) o += ch === '\t' ? ' '.repeat(8 - (o.length % 8)) : ch; l = o; }
    let name, data;
    if (l.includes(',')) {
      const f = l.split(',').map((s) => s.trim());
      name = f[0];
      const cnt = /\*/.test(name) ? 4 : 8;
      data = f.slice(1, 1 + cnt);
      if (data.length && /^[+*]/.test(data[data.length - 1])) data.pop();
      while (data.length < cnt) data.push('');
    } else {
      name = l.slice(0, 8).trim();
      const big = name.endsWith('*') || name.startsWith('*'), w = big ? 16 : 8, cnt = big ? 4 : 8;
      data = [];
      for (let k = 0; k < cnt; k++) data.push(l.slice(8 + k * w, 8 + (k + 1) * w).trim());
    }
    if (name === '' || name[0] === '+' || name[0] === '*') { if (cur) cur.f.push(...data); }
    else { cur = { name: name.replace(/\*$/, '').trim().toUpperCase(), f: data }; if (cur.name === 'GRID' || NAS_CELL[cur.name]) cards.push(cur); if (cards.length > 3 * L.cells) fail('The Nastran deck is too large.'); }
  }
  const idmap = new Map(), coords = [], store = cellStore();
  let local = 0;
  for (const c of cards) if (c.name === 'GRID') { idmap.set(parseInt(c.f[0], 10), coords.length / 3); if (parseInt(c.f[1], 10) > 0) local++; coords.push(num(c.f[2]) || 0, num(c.f[3]) || 0, num(c.f[4]) || 0); }
  if (!coords.length) fail('No GRID cards were found: not a Nastran bulk data file.');
  for (const c of cards) {
    const d = NAS_CELL[c.name];
    if (!d) continue;
    const ids = c.f.slice(2, 2 + d[1]).map((s) => idmap.get(parseInt(s, 10)));
    if (ids.length < d[1] || ids.some((q) => q === undefined)) fail(`A ${c.name} card refers to a missing GRID point.`);
    store.add(d[0], ids);
  }
  return cellMesh(Float64Array.from(coords), store, { warnings: local ? [`${local} GRID points use a local coordinate system, which is not applied.`] : [] });
}

const TEC_ET = { TRIANGLE: 3, QUADRILATERAL: 4, TETRAHEDRON: 4, BRICK: 8, LINESEG: 2 };
async function readTecplot(ctx) {
  const zones = [];
  let varsText = '', mode = '', cur = null, hdr = '';
  for (const rawLine of splitLines(await ctx.text())) {
    const l = rawLine.trim();
    if (!l || l[0] === '#') continue;
    if (/^[-+.\d]/.test(l) && !l.includes('=')) { if (mode === 'zone') { cur = { hdr, rows: [] }; zones.push(cur); mode = ''; } if (cur) cur.rows.push(l); continue; }
    const up = l.toUpperCase();
    if (/^ZONE\b/.test(up)) { mode = 'zone'; hdr = l.slice(4); }
    else if (/^VARIABLES\b/.test(up)) { mode = 'vars'; varsText = l.slice(9); }
    else if (/^(TITLE|TEXT|GEOMETRY|DATASETAUXDATA|CUSTOMLABELS|FILETYPE)\b/.test(up)) mode = 'other';
    else if (mode === 'zone') hdr += ' ' + l; else if (mode === 'vars') varsText += ' ' + l;
  }
  if (!zones.length) fail('No ZONE with data was found: not a Tecplot ASCII file.');
  let names = (varsText.match(/"[^"]*"/g) || []).map((s) => s.slice(1, -1));
  if (!names.length) names = varsText.replace(/=/, ' ').split(/[\s,]+/).filter(Boolean);
  const find = (re) => names.findIndex((s) => re.test(s.trim()));
  let xi = find(/^(x|coordinatex)(\s*[[(].*)?$/i), yi = find(/^(y|coordinatey)(\s*[[(].*)?$/i), zi = find(/^(z|coordinatez)(\s*[[(].*)?$/i);
  const store = cellStore(), blocks = [], parts = [], warnings = [];
  let base = 0, any2 = false;
  for (const z of zones) {
    const kv = {}, re = /([A-Za-z]+)\s*=\s*("[^"]*"|\([^)]*\)|[^\s,]+)/g;
    let m;
    while ((m = re.exec(z.hdr))) kv[m[1].toUpperCase()] = m[2].replace(/"/g, '');
    const F = (kv.F || '').toUpperCase(), zt = (kv.ZONETYPE || '').toUpperCase(), v = numsOf(z.rows.join(' '));
    if (/POLY/.test(zt)) fail('Tecplot polygonal / polyhedral zones are not supported; export tetrahedral or brick zones, or VTK.');
    const fe = zt.startsWith('FE') || F.startsWith('FE'), block = (kv.DATAPACKING || '').toUpperCase() === 'BLOCK' || F === 'FEBLOCK' || F === 'BLOCK';
    const etName = fe ? (zt.startsWith('FE') ? zt.slice(2) : (kv.ET || 'QUADRILATERAL').toUpperCase()) : '', npe = TEC_ET[etName];
    if (fe && !npe) fail(`Tecplot element type "${etName}" is not supported.`);
    const I = parseInt(kv.I || '1', 10), J = parseInt(kv.J || '1', 10), K = parseInt(kv.K || '1', 10);
    const nNode = fe ? parseInt(kv.N || kv.NODES, 10) : I * J * K, nEl = fe ? parseInt(kv.E || kv.ELEMENTS, 10) : Math.max(1, I - 1) * Math.max(1, J - 1) * Math.max(1, K - 1);
    if (!(nNode > 0) || !(nEl >= 0) || nNode > v.length) fail('A Tecplot zone header is missing its node / element counts, or the data is truncated.');
    let nv = names.length;
    if (!nv) { nv = block ? Math.round((v.length - (fe ? nEl * npe : 0)) / nNode) : z.rows[0].split(/[\s,]+/).length; }
    if (xi < 0) { xi = 0; yi = 1; zi = nv >= 3 && (K > 1 || npe === 8 || etName === 'TETRAHEDRON' || names.length === 3) ? 2 : -1; }
    if (yi < 0 || xi >= nv || yi >= nv || zi >= nv) fail('The Tecplot file has no recognisable X / Y coordinate variables.');
    const cc = new Set();
    if (kv.VARLOCATION && /CELLCENTERED/i.test(kv.VARLOCATION)) for (const g of kv.VARLOCATION.matchAll(/\[([^\]]*)\]\s*=\s*CELLCENTERED/gi)) for (const r of g[1].split(',')) { const [a, b] = r.split('-').map(Number); for (let q = a; q <= Math.min(b || a, 4096); q++) cc.add(q - 1); }
    if (cc.has(xi) || cc.has(yi) || cc.has(zi)) fail('Cell-centred coordinate variables are not supported.');
    const off = [];
    let tot = 0;
    for (let q = 0; q < nv; q++) { off.push(tot); tot += block ? (cc.has(q) ? nEl : nNode) : 1; }
    const nodeVals = block ? tot : nNode * nv;
    if (v.length < nodeVals + (fe ? nEl * npe : 0)) fail('A Tecplot zone holds less data than its header declares.');
    const val = block ? (n, q) => v[off[q] + n] : (n, q) => v[n * nv + q], xyz = new Float64Array(nNode * 3);
    for (let n = 0; n < nNode; n++) { xyz[3 * n] = val(n, xi); xyz[3 * n + 1] = val(n, yi); xyz[3 * n + 2] = zi >= 0 ? val(n, zi) : 0; }
    if (zi < 0) any2 = true;
    parts.push(xyz);
    if (!fe) blocks.push({ xyz, ni: I, nj: J, nk: K, base });
    else for (let e = 0; e < nEl; e++) {
      const ids = [];
      for (let q = 0; q < npe; q++) ids.push(v[nodeVals + e * npe + q] - 1 + base);
      const u = ids.filter((x, q) => ids.indexOf(x) === q);
      if (etName === 'BRICK') { if (u.length === 8) store.add('hex', ids); else if (u.length === 4) store.add('tet', u); else if (u.length === 5) store.add('pyramid', u); else if (u.length === 6) store.add('wedge', u); }
      else if (etName === 'TETRAHEDRON') store.add('tet', ids);
      else if (npe === 2) store.add('line', ids);
      else if (u.length >= 3) store.add(u.length === 3 ? 'tri' : 'quad', u);
    }
    base += nNode;
  }
  const stats = { zones: zones.length, variables: names };
  if (!store.t.length) return structuredGeom(blocks, { warnings, stats });
  const all = new Float64Array(base * 3);
  let o = 0;
  for (const p of parts) { all.set(p, o); o += p.length; }
  for (const b of blocks) {                    // mixed files: ordered zones join the finite-element cells
    const id = (i, j, k) => b.base + i + b.ni * (j + b.nj * k);
    if (b.ni > 1 && b.nj > 1 && b.nk > 1) { for (let k = 0; k < b.nk - 1; k++) for (let j = 0; j < b.nj - 1; j++) for (let i = 0; i < b.ni - 1; i++) store.add('hex', [id(i, j, k), id(i + 1, j, k), id(i + 1, j + 1, k), id(i, j + 1, k), id(i, j, k + 1), id(i + 1, j, k + 1), id(i + 1, j + 1, k + 1), id(i, j + 1, k + 1)]); }
    else if (b.nk === 1 && b.ni > 1 && b.nj > 1) for (let j = 0; j < b.nj - 1; j++) for (let i = 0; i < b.ni - 1; i++) store.add('quad', [id(i, j, 0), id(i + 1, j, 0), id(i + 1, j + 1, 0), id(i, j + 1, 0)]);
  }
  return cellMesh(all, store, { flat: any2 ? true : undefined, warnings, stats });
}

/** Plot3D formatted grid tokens -> blocks, or null when the token count fits no single/multi-block 2-D/3-D layout. */
function plot3dBlocks(v) {
  const ok = (a) => a.every((k) => Number.isInteger(k) && k > 0 && k < 1e7);
  const build = (dimsAt, nb, nd) => {
    const dims = [];
    let tot = dimsAt + nb * nd;
    for (let b = 0; b < nb; b++) { const d = Array.from(v.subarray(dimsAt + b * nd, dimsAt + (b + 1) * nd)); if (d.length < nd || !ok(d)) return null; if (nd === 2) d.push(1); dims.push(d); tot += nd * d[0] * d[1] * d[2]; }
    if (tot !== v.length) return null;
    const blocks = [];
    let p = dimsAt + nb * nd;
    for (const [ni, nj, nk] of dims) {
      const n = ni * nj * nk, xyz = new Float64Array(3 * n);
      for (let q = 0; q < n; q++) { xyz[3 * q] = v[p + q]; xyz[3 * q + 1] = v[p + n + q]; xyz[3 * q + 2] = nd === 3 ? v[p + 2 * n + q] : 0; }
      p += nd * n;
      blocks.push({ xyz, ni, nj, nk });
    }
    return blocks;
  };
  if (v.length < 6) return null;
  const nb = v[0];
  return build(0, 1, 3) || build(0, 1, 2) || (Number.isInteger(nb) && nb > 0 && nb <= 4096 ? build(1, nb, 3) || build(1, nb, 2) : null);
}
async function readPlot3D(ctx) {
  const u8 = await ctx.bytes();
  if (u8.subarray(0, 4096).some((b) => b === 0)) fail('Unformatted (binary) Plot3D grids are not supported. ' + fmtNamed('Plot3D grid (ASCII)').convert);
  const blocks = plot3dBlocks(numsOf(await ctx.text()));
  if (!blocks) fail('The file does not match a formatted Plot3D grid (whole, no IBLANK). ' + fmtNamed('Plot3D grid (ASCII)').convert);
  return structuredGeom(blocks, {});
}

async function readFoam(ctx) {
  const part = {}, strip = (n) => baseName(n).replace(/\.gz$/, '');
  for (const k of FOAM_FILES) part[k] = strip(ctx.name) === k ? ctx.file : ctx.companion((n) => strip(n) === k);
  if (!part.points || !part.faces) fail(fmtNamed('OpenFOAM polyMesh').convert);
  const load = async (f) => {
    let u = new Uint8Array(await f.arrayBuffer());
    if (isGzip(u)) u = await inflate(u, 'gzip');
    let t = utf8.decode(u);
    const head = t.slice(0, 4000), note = /nInternalFaces:\s*(\d+)/.exec(head);
    if (/format\s+binary/.test(head)) fail('Binary OpenFOAM meshes are not supported; run foamFormatConvert to write ASCII.');
    t = t.replace(/\/\*[\s\S]*?\*\//g, ' ').replace(/\/\/[^\n]*/g, ' ').replace(/FoamFile\s*\{[^}]*\}/, ' ');
    return { v: numsOf(t.replace(/[()]/g, ' ')), nInt: note ? +note[1] : null };
  };
  const P = await load(part.points), Fc = await load(part.faces), O = part.owner ? await load(part.owner) : null, N = part.neighbour ? await load(part.neighbour) : null;
  const np = P.v[0], nf = Fc.v[0];
  if (!Number.isInteger(np) || P.v.length < 1 + 3 * np || !Number.isInteger(nf)) fail('The OpenFOAM points / faces files are truncated or malformed.');
  const xyz = P.v.subarray(1, 1 + 3 * np), warnings = [], tri = [];
  let nInt = N ? N.v[0] : (O && O.nInt) ?? Fc.nInt;
  if (!Number.isInteger(nInt)) { nInt = 0; warnings.push('No neighbour file or nInternalFaces note was supplied: every face (internal ones included) is kept.'); }
  let q = 1;
  for (let f = 0; f < nf; f++) {
    const k = Fc.v[q];
    if (!(k >= 3) || q + 1 + k > Fc.v.length) fail('The OpenFOAM faces file is truncated or malformed.');
    if (f >= nInt && !pushFan(tri, xyz, Array.from(Fc.v.subarray(q + 1, q + 1 + k)))) fail('An OpenFOAM face refers to a missing point.');
    q += 1 + k;
  }
  let cells;
  if (O) { cells = 0; for (let k = 1; k < O.v.length; k++) if (O.v[k] >= cells) cells = O.v[k] + 1; }
  return meshGeom(tri, { warnings, stats: { nodes: np, faces: nf, internalFaces: nInt, cells, cellTypes: { polyhedron: cells }, dimension: 3 } });
}

// ---- 2-D drawings --------------------------------------------------------------------------------------
function linesOrPoints(polys, pts, extra = {}) {
  if (polys.length) { const g = polyGeom(polys, extra); if (pts.length) { g.markers = []; for (let i = 0; i < pts.length; i += 3) g.markers.push([pts[i], pts[i + 1]]); g.stats = { ...(g.stats || {}), points: pts.length / 3 }; } return g; }
  if (pts.length) return pointsGeom(pts, extra);
  return fail('No line, polygon or point geometry was found in the file.');
}
const closeRing = (x, y) => { const n = x.length; if (n > 2 && x[0] === x[n - 1] && y[0] === y[n - 1]) { x.pop(); y.pop(); return true; } return false; };
const arcPts = (cx, cy, rx, ry, t0, sweep, rot = 0) => {
  const n = Math.min(4 * ARC_N, Math.max(2, Math.ceil((Math.abs(sweep) / TAU) * ARC_N - 1e-9) || 2)), c = Math.cos(rot), s = Math.sin(rot), out = [];
  for (let i = 0; i <= n; i++) { const t = t0 + (sweep * i) / n, u = rx * Math.cos(t), v = ry * Math.sin(t); out.push([cx + u * c - v * s, cy + u * s + v * c]); }
  return out;
};

async function readDXF(ctx) {
  const text = await ctx.text();
  if (text.startsWith('AutoCAD Binary DXF')) fail('Binary DXF is not supported. ' + fmtNamed('DXF').convert);
  let polys = [];
  try { polys = parseDXF(text).polylines.slice(); } catch { polys = []; }
  const ls = splitLines(text), tri = [], count = {}, warnings = [];
  let d = null, type = '';
  const flush = () => {
    if (!d) return;
    const g = (c, k = 0, dflt = NaN) => (d[c] && d[c][k] !== undefined ? d[c][k] : dflt);
    let pts = null, closed = false;
    if (type === 'ARC') { const a0 = (g(50, 0, 0) * Math.PI) / 180; let sw = ((g(51, 0, 360) * Math.PI) / 180 - a0) % TAU; if (sw <= 0) sw += TAU; pts = arcPts(g(10), g(20), g(40), g(40), a0, sw); }
    else if (type === 'ELLIPSE') {
      const mx = g(11), my = g(21), a = Math.hypot(mx, my), t0 = g(41, 0, 0);
      let sw = g(42, 0, TAU) - t0;
      if (sw <= 0 || sw > TAU) sw = ((sw % TAU) + TAU) % TAU || TAU;
      pts = arcPts(g(10), g(20), a, a * g(40, 0, 1), t0, sw, Math.atan2(my, mx));
      if (sw >= TAU - 1e-9) { pts.pop(); closed = true; }
    } else if (type === 'SPLINE') {
      const cps = (d[10] || []).map((x, k) => [x, g(20, k, 0), 0]), deg = g(71, 0, 3), fit = (d[11] || []).map((x, k) => [x, g(21, k, 0)]);
      const smp = cps.length ? nurbsCurve(deg, d[40] || [], cps, d[41] && d[41].length === cps.length ? d[41] : null, NaN, NaN, 8 * cps.length) : null;
      pts = smp || (fit.length > 1 ? fit : cps);
      closed = (g(70, 0, 0) & 1) === 1 && !smp;
    } else if (type === '3DFACE') {
      const p = [0, 1, 2, 3].map((k) => [g(10 + k), g(20 + k), g(30 + k, 0, 0)]);
      if (p.slice(0, 3).every((q) => q.every(Number.isFinite))) { tri.push(...p[0], ...p[1], ...p[2]); if (p[3].every(Number.isFinite) && p[3].some((v, k) => v !== p[2][k])) tri.push(...p[0], ...p[2], ...p[3]); }
    }
    count[type] = (count[type] || 0) + 1;
    if (pts && pts.length > 1 && pts.every((q) => Number.isFinite(q[0]) && Number.isFinite(q[1]))) polys.push({ x: pts.map((q) => q[0]), y: pts.map((q) => q[1]), closed });
    d = null;
  };
  for (let i = 0; i + 1 < ls.length; i += 2) {
    const code = parseInt(ls[i], 10), val = ls[i + 1].trim();
    if (code === 0) { flush(); type = val; d = ['ARC', 'ELLIPSE', 'SPLINE', '3DFACE'].includes(val) ? {} : null; }
    else if (d && code > 0 && code < 300) (d[code] || (d[code] = [])).push(+val);
  }
  flush();
  if (!polys.length && !tri.length) fail('No LINE, POLYLINE, CIRCLE, ARC, ELLIPSE, SPLINE or 3DFACE entities were found in the DXF file.');
  if (/\n\s*INSERT\s*\r?\n/.test(text)) warnings.push('Block references (INSERT) are not expanded; explode blocks before exporting.');
  if (!polys.length) return meshGeom(tri, { warnings, stats: { faces3d: count['3DFACE'] } });
  if (tri.length) warnings.push(`${count['3DFACE']} 3DFACE entities were ignored in favour of the 2-D outlines.`);
  return polyGeom(polys, { warnings, stats: { entities: count } });
}

/** SVG path data -> [{ pts, closed }] with curves flattened. */
function svgPath(dstr) {
  const d = String(dstr || ''), n = d.length, out = [], numRe = /[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?/y;
  let i = 0, cur = null, cx = 0, cy = 0, sx = 0, sy = 0, px = 0, py = 0, prev = '', guard = 0;
  const skip = () => { while (i < n && (d.charCodeAt(i) <= 32 || d[i] === ',')) i++; };
  const num = () => { skip(); numRe.lastIndex = i; const m = numRe.exec(d); if (!m) return NaN; i = numRe.lastIndex; return +m[0]; };
  const flag = () => { skip(); const c = d[i]; if (c === '0' || c === '1') { i++; return +c; } return NaN; };
  const more = () => { skip(); return i < n && /[-+.\d]/.test(d[i]); };
  const lineTo = (x, y) => { if (!cur) { cur = { pts: [[cx, cy]], closed: false }; out.push(cur); } cur.pts.push([x, y]); cx = x; cy = y; };
  const bez = (p0, p1, p2, p3) => { const N = p3 ? 16 : 12; for (let k = 1; k <= N; k++) { const t = k / N, u = 1 - t; if (p3) lineTo(u * u * u * p0[0] + 3 * u * u * t * p1[0] + 3 * u * t * t * p2[0] + t * t * t * p3[0], u * u * u * p0[1] + 3 * u * u * t * p1[1] + 3 * u * t * t * p2[1] + t * t * t * p3[1]); else lineTo(u * u * p0[0] + 2 * u * t * p1[0] + t * t * p2[0], u * u * p0[1] + 2 * u * t * p1[1] + t * t * p2[1]); } };
  while (i < n) {
    skip();
    if (i >= n || !/[a-zA-Z]/.test(d[i])) break;
    const cmd = d[i++], rel = cmd === cmd.toLowerCase(), ox = () => (rel ? cx : 0), oy = () => (rel ? cy : 0);
    let C = cmd.toUpperCase(), first = true;
    if (C === 'Z') { if (cur) cur.closed = true; cur = null; cx = sx; cy = sy; prev = 'Z'; continue; }
    while (first || more()) {
      first = false;
      if (++guard > 2e6) fail('An SVG path is too long.');
      const bx = ox(), by = oy();
      let v;
      if (C === 'M') { v = [num() + bx, num() + by]; if (!v.every(Number.isFinite)) return out; cur = null; cx = sx = v[0]; cy = sy = v[1]; cur = { pts: [[cx, cy]], closed: false }; out.push(cur); prev = 'M'; C = 'L'; continue; }
      if (C === 'L') v = [num() + bx, num() + by];
      else if (C === 'H') v = [num() + bx, cy];
      else if (C === 'V') v = [cx, num() + by];
      else if (C === 'C') v = [num() + bx, num() + by, num() + bx, num() + by, num() + bx, num() + by];
      else if (C === 'S') v = [num() + bx, num() + by, num() + bx, num() + by];
      else if (C === 'Q') v = [num() + bx, num() + by, num() + bx, num() + by];
      else if (C === 'T') v = [num() + bx, num() + by];
      else if (C === 'A') v = [num(), num(), num(), flag(), flag(), num() + bx, num() + by];
      else return out;
      if (!v.every(Number.isFinite)) return out;
      if (C === 'L' || C === 'H' || C === 'V') lineTo(v[0], v[1]);
      else if (C === 'C') { bez([cx, cy], [v[0], v[1]], [v[2], v[3]], [v[4], v[5]]); px = v[2]; py = v[3]; }
      else if (C === 'S') { const r = prev === 'C' || prev === 'S' ? [2 * cx - px, 2 * cy - py] : [cx, cy]; bez([cx, cy], r, [v[0], v[1]], [v[2], v[3]]); px = v[0]; py = v[1]; }
      else if (C === 'Q') { bez([cx, cy], [v[0], v[1]], [v[2], v[3]]); px = v[0]; py = v[1]; }
      else if (C === 'T') { const r = prev === 'Q' || prev === 'T' ? [2 * cx - px, 2 * cy - py] : [cx, cy]; bez([cx, cy], r, [v[0], v[1]]); px = r[0]; py = r[1]; }
      else {
        let rx = Math.abs(v[0]), ry = Math.abs(v[1]);
        const phi = (v[2] * Math.PI) / 180, x2 = v[5], y2 = v[6];
        if (!rx || !ry || (cx === x2 && cy === y2)) lineTo(x2, y2);
        else {
          const cp = Math.cos(phi), sp = Math.sin(phi), dx = (cx - x2) / 2, dy = (cy - y2) / 2, x1p = cp * dx + sp * dy, y1p = -sp * dx + cp * dy, lam = (x1p * x1p) / (rx * rx) + (y1p * y1p) / (ry * ry);
          if (lam > 1) { rx *= Math.sqrt(lam); ry *= Math.sqrt(lam); }
          const den = rx * rx * y1p * y1p + ry * ry * x1p * x1p;
          let co = Math.sqrt(Math.max(0, (rx * rx * ry * ry - den) / den));
          if (v[3] === v[4]) co = -co;
          const cxp = (co * rx * y1p) / ry, cyp = (-co * ry * x1p) / rx, ccx = cp * cxp - sp * cyp + (cx + x2) / 2, ccy = sp * cxp + cp * cyp + (cy + y2) / 2;
          const ang = (ux, uy, vx, vy) => Math.atan2(ux * vy - uy * vx, ux * vx + uy * vy), t1 = ang(1, 0, (x1p - cxp) / rx, (y1p - cyp) / ry);
          let dt = ang((x1p - cxp) / rx, (y1p - cyp) / ry, (-x1p - cxp) / rx, (-y1p - cyp) / ry);
          if (!v[4] && dt > 0) dt -= TAU; else if (v[4] && dt < 0) dt += TAU;
          const pts = arcPts(ccx, ccy, rx, ry, t1, dt, phi);
          for (let k = 1; k < pts.length; k++) if (k === pts.length - 1) lineTo(x2, y2); else lineTo(pts[k][0], pts[k][1]);
        }
      }
      prev = C;
    }
  }
  return out;
}
const SVG_SKIP = new Set(['defs', 'clippath', 'mask', 'pattern', 'symbol', 'marker', 'metadata', 'style', 'script', 'title', 'desc', 'text', 'lineargradient', 'radialgradient', 'filter', 'image']);
async function readSVG(ctx, doc) {
  doc = doc || parseXML(await ctx.text());
  const svg = xAll(doc, 'svg')[0];
  if (!svg) fail('Not an SVG file (no <svg> element).');
  const mul = (m, k) => [m[0] * k[0] + m[2] * k[1], m[1] * k[0] + m[3] * k[1], m[0] * k[2] + m[2] * k[3], m[1] * k[2] + m[3] * k[3], m[0] * k[4] + m[2] * k[5] + m[4], m[1] * k[4] + m[3] * k[5] + m[5]];
  const xform = (s) => {
    let M = [1, 0, 0, 1, 0, 0], m;
    const re = /(matrix|translate|scale|rotate|skewx|skewy)\s*\(([^)]*)\)/gi;
    while ((m = re.exec(s || ''))) {
      const a = numsOf(m[2]), k = m[1].toLowerCase(), r = (a[0] * Math.PI) / 180;
      let T = null;
      if (k === 'matrix' && a.length >= 6) T = Array.from(a.subarray(0, 6));
      else if (k === 'translate') T = [1, 0, 0, 1, a[0] || 0, a[1] || 0];
      else if (k === 'scale' && a.length) T = [a[0], 0, 0, a.length > 1 ? a[1] : a[0], 0, 0];
      else if (k === 'rotate' && a.length) { T = [Math.cos(r), Math.sin(r), -Math.sin(r), Math.cos(r), 0, 0]; if (a.length >= 3) T = mul(mul([1, 0, 0, 1, a[1], a[2]], T), [1, 0, 0, 1, -a[1], -a[2]]); }
      else if (k === 'skewx' && a.length) T = [1, 0, Math.tan(r), 1, 0, 0];
      else if (k === 'skewy' && a.length) T = [1, Math.tan(r), 0, 1, 0, 0];
      if (T && T.every(Number.isFinite)) M = mul(M, T);
    }
    return M;
  };
  const polys = [], A = (nd, k, dflt = 0) => { const v = parseFloat(nd.attrs[k]); return Number.isFinite(v) ? v : dflt; };
  const add = (pts, closed, M) => {
    if (pts.length < 2) return;
    const x = pts.map((p) => M[0] * p[0] + M[2] * p[1] + M[4]), y = pts.map((p) => M[1] * p[0] + M[3] * p[1] + M[5]);
    if (closed) closeRing(x, y);
    polys.push({ x, y, closed });
    if (polys.length > 2e5) fail('The SVG holds too many shapes.');
  };
  const walk = (node, M) => {
    for (const c of node.children) {
      if (SVG_SKIP.has(c.name)) continue;
      const T = c.attrs.transform ? mul(M, xform(c.attrs.transform)) : M;
      if (c.name === 'path') for (const s of svgPath(c.attrs.d)) add(s.pts, s.closed, T);
      else if (c.name === 'polygon' || c.name === 'polyline') { const v = numsOf(c.attrs.points || ''), pts = []; for (let k = 0; k + 1 < v.length; k += 2) pts.push([v[k], v[k + 1]]); add(pts, c.name === 'polygon', T); }
      else if (c.name === 'rect') { const x = A(c, 'x'), y = A(c, 'y'), w = A(c, 'width'), h = A(c, 'height'); if (w > 0 && h > 0) add([[x, y], [x + w, y], [x + w, y + h], [x, y + h]], true, T); }
      else if (c.name === 'circle' || c.name === 'ellipse') { const rx = A(c, c.name === 'circle' ? 'r' : 'rx'), ry = A(c, c.name === 'circle' ? 'r' : 'ry'); if (rx > 0 && ry > 0) add(arcPts(A(c, 'cx'), A(c, 'cy'), rx, ry, 0, TAU).slice(0, -1), true, T); }
      else if (c.name === 'line') add([[A(c, 'x1'), A(c, 'y1')], [A(c, 'x2'), A(c, 'y2')]], false, T);
      if (c.children.length) walk(c, T);
    }
  };
  walk(svg, svg.attrs.transform ? xform(svg.attrs.transform) : [1, 0, 0, 1, 0, 0]);
  const g = polyGeom(polys, { warnings: ['SVG user units are returned as written, with the y axis flipped to point upwards.'] });
  const s = g.bbox.min[1] + g.bbox.max[1];
  for (const p of g.polylines) p.y = p.y.map((v) => s - v);
  return g;
}

async function readHPGL(ctx) {
  const text = (await ctx.text()).replace(/\x1b\.[^A-Za-z]*[A-Za-z][^:;A-Za-z]*:?/g, ''), re = /([A-Za-z]{2})([^A-Za-z;]*)/g, out = [], U = 0.025;
  let m, down = false, rel = false, x = 0, y = 0, cur = null, cmds = 0;
  const end = () => { if (cur && cur.length > 1) out.push(cur); cur = null; };
  const move = (nx, ny) => { if (down) { if (!cur) cur = [[x, y]]; cur.push([nx, ny]); } else end(); x = nx; y = ny; };
  const pairs = (a) => { for (let k = 0; k + 1 < a.length; k += 2) move(rel ? x + a[k] : a[k], rel ? y + a[k + 1] : a[k + 1]); };
  while ((m = re.exec(text))) {
    const c = m[1].toUpperCase(), a = numsOf(m[2]);
    if (a.some((v) => !Number.isFinite(v))) continue;
    cmds++;
    if (c === 'PU') { end(); down = false; pairs(a); }
    else if (c === 'PD') { down = true; pairs(a); }
    else if (c === 'PA') { rel = false; pairs(a); }
    else if (c === 'PR') { rel = true; pairs(a); }
    else if (c === 'CI' && a.length) { end(); const p = arcPts(x, y, Math.abs(a[0]), Math.abs(a[0]), 0, TAU); p.pop(); p.closed = true; out.push(p); }
    else if ((c === 'AA' || c === 'AR') && a.length >= 3) { const ccx = c === 'AR' ? x + a[0] : a[0], ccy = c === 'AR' ? y + a[1] : a[1], r = Math.hypot(x - ccx, y - ccy), p = arcPts(ccx, ccy, r, r, Math.atan2(y - ccy, x - ccx), (a[2] * Math.PI) / 180), wasDown = down; down = true; for (let k = 1; k < p.length; k++) move(p[k][0], p[k][1]); down = wasDown; if (!down) end(); }
    else if ((c === 'EA' || c === 'ER') && a.length >= 2) { end(); const x2 = c === 'ER' ? x + a[0] : a[0], y2 = c === 'ER' ? y + a[1] : a[1], p = [[x, y], [x2, y], [x2, y2], [x, y2]]; p.closed = true; out.push(p); }
    else if (c === 'IN' || c === 'DF') { end(); down = false; rel = false; }
    if (out.length > 5e5) fail('The plot file holds too many strokes.');
  }
  end();
  if (!out.length) fail('No pen-down strokes were found: not an HPGL plot file.');
  const polys = out.map((p) => { const xs = p.map((q) => q[0] * U), ys = p.map((q) => q[1] * U); return { x: xs, y: ys, closed: !!p.closed || closeRing(xs, ys) }; });
  return polyGeom(polys, { warnings: ['Plotter units (0.025 mm) converted to millimetres.'], stats: { units: 'mm', commands: cmds } });
}

// ---- GIS vector formats ------------------------------------------------------------------------------
function kmlGeom(doc) {
  const polys = [], pts = [];
  const coords = (nd) => { const v = []; for (const t of nd.text.match(/\S+/g) || []) { const c = t.split(',').map(Number); if (c.length >= 2 && Number.isFinite(c[0]) && Number.isFinite(c[1])) v.push(c); } return v; };
  for (const nd of xAll(doc, 'linearring', 'linestring', 'point', 'track')) {
    const v = nd.name === 'track' ? xKids(nd, 'coord').map((c) => Array.from(numsOf(c.text))).filter((c) => c.length >= 2) : coords(xAll(nd, 'coordinates')[0] || { text: '' });
    if (nd.name === 'point') { if (v[0]) pts.push(v[0][0], v[0][1], v[0][2] || 0); continue; }
    if (v.length < 2) continue;
    const x = v.map((c) => c[0]), y = v.map((c) => c[1]), closed = nd.name === 'linearring';
    if (closed) closeRing(x, y);
    polys.push(!closed && v.some((c) => c[2]) && v.every((c) => Number.isFinite(c[2] ?? 0)) ? { x, y, z: v.map((c) => c[2] ?? 0), closed } : { x, y, closed });
  }
  return linesOrPoints(polys, pts, { geographic: true, warnings: ['Coordinates are geographic (longitude, latitude in degrees).'] });
}
async function readKML(ctx, doc) {
  if (!doc) {
    let u8 = await ctx.bytes();
    if (isZip(u8)) { const f = await unzipFiles(await ctx.buf(), (n) => /\.kml$/i.test(n)); u8 = f.get('doc.kml') || f.values().next().value; if (!u8) fail('No .kml document was found inside the KMZ archive.'); }
    doc = parseXML(utf8.decode(u8));
  }
  return kmlGeom(doc);
}
async function readGPX(ctx, doc) {
  doc = doc || parseXML(await ctx.text());
  const polys = [], pts = [], ll = (p) => [+p.attrs.lon, +p.attrs.lat, +((xKid(p, 'ele') || {}).text || 0) || 0];
  for (const seg of [...xAll(doc, 'trkseg'), ...xAll(doc, 'rte')]) { const v = seg.children.filter((c) => c.name === 'trkpt' || c.name === 'rtept').map(ll); if (v.length > 1) polys.push({ x: v.map((c) => c[0]), y: v.map((c) => c[1]), z: v.map((c) => c[2]), closed: false }); }
  for (const w of xAll(doc, 'wpt')) pts.push(...ll(w));
  return linesOrPoints(polys, pts, { geographic: true, warnings: ['Coordinates are geographic (longitude, latitude in degrees).'] });
}
async function readGML(ctx, doc) {
  doc = doc || parseXML(await ctx.text());
  const polys = [], pts = [];
  for (const nd of xAll(doc, 'linearring', 'linestring', 'linestringsegment', 'point')) {
    const pl = xKid(nd, 'poslist'), ps = xKids(nd, 'pos'), co = xKid(nd, 'coordinates');
    let v = [];
    if (pl) { const a = numsOf(pl.text), dim = +(pl.attrs.srsdimension || nd.attrs.srsdimension) || (a.length % 2 === 0 ? 2 : 3); for (let k = 0; k + dim <= a.length; k += dim) v.push([a[k], a[k + 1], dim > 2 ? a[k + 2] : 0]); }
    else if (ps.length) v = ps.map((p) => Array.from(numsOf(p.text)));
    else if (co) { const cs = co.attrs.cs || ',', ts = co.attrs.ts || ' '; v = co.text.trim().split(ts === ' ' ? /\s+/ : ts).map((t) => t.split(cs).map(Number)); }
    v = v.filter((c) => c.length >= 2 && Number.isFinite(c[0]) && Number.isFinite(c[1]));
    if (nd.name === 'point') { if (v[0]) pts.push(v[0][0], v[0][1], v[0][2] || 0); continue; }
    if (v.length < 2) continue;
    const x = v.map((c) => c[0]), y = v.map((c) => c[1]), closed = nd.name === 'linearring';
    if (closed) closeRing(x, y);
    polys.push({ x, y, closed });
  }
  return linesOrPoints(polys, pts, { warnings: ['GML axis order is taken as written (first value = x); geographic CRSs such as EPSG:4326 may list latitude first.'] });
}

async function readSHP(ctx) {
  const u8 = await ctx.bytes(), dv = new DataView(u8.buffer, u8.byteOffset, u8.byteLength);
  if (u8.length < 100 || dv.getInt32(0, false) !== 9994) fail('Not an ESRI shapefile (.shp main file).');
  const polys = [], pts = [], warnings = [], f64 = (o) => dv.getFloat64(o, true), i32 = (o) => dv.getInt32(o, true);
  let p = 100, nrec = 0, skipped = 0, total = 0, cut = false;
  while (p + 12 <= u8.length) {
    const len = dv.getInt32(p + 4, false) * 2, s = p + 8;
    if (len < 4 || s + len > u8.length) { cut = true; break; }
    const t = i32(s), base = t % 10, hasZ = t >= 10 && t < 20;
    nrec++;
    if (t === 0) { /* null shape */ }
    else if (base === 1 && len >= 20) pts.push(f64(s + 4), f64(s + 12), hasZ && len >= 28 ? f64(s + 20) : 0);
    else if (base === 8 && len >= 40) {
      const n = i32(s + 36), zo = s + 40 + 16 * n + 16;
      if (n < 0 || 40 + 16 * n > len) { cut = true; break; }
      for (let k = 0; k < n; k++) pts.push(f64(s + 40 + 16 * k), f64(s + 48 + 16 * k), hasZ && zo + 8 * n <= s + len ? f64(zo + 8 * k) : 0);
    } else if ((base === 3 || base === 5) && len >= 44) {
      const np = i32(s + 36), n = i32(s + 40), po = s + 44, co = po + 4 * np;
      if (np < 0 || n < 0 || 44 + 4 * np + 16 * n > len) { cut = true; break; }
      total += n;
      if (total > 5e6) fail('The shapefile holds too many vertices.');
      for (let k = 0; k < np; k++) {
        const a = i32(po + 4 * k), b = k + 1 < np ? i32(po + 4 * k + 4) : n, x = [], y = [], zo = co + 16 * n + 16, z = hasZ && base === 3 && zo + 8 * n <= s + len ? [] : null;
        if (a < 0 || b > n || b < a) continue;
        for (let q = a; q < b; q++) { x.push(f64(co + 16 * q)); y.push(f64(co + 16 * q + 8)); if (z) z.push(f64(zo + 8 * q)); }
        if (base === 5) closeRing(x, y);
        if (x.length > 1) polys.push(z ? { x, y, z, closed: false } : { x, y, closed: base === 5 });
      }
    } else skipped++;
    p = s + len;
  }
  if (cut) warnings.push(`The shapefile is truncated; ${nrec} records were recovered.`);
  if (skipped) warnings.push(`${skipped} records of unsupported shape type (e.g. MultiPatch) were skipped.`);
  const stats = { records: nrec, shapeType: i32(32) }, dbf = ctx.companion((n) => /\.dbf$/i.test(n)), prj = ctx.companion((n) => /\.prj$/i.test(n));
  if (dbf) {
    const b = new Uint8Array(await dbf.arrayBuffer());
    if (b.length >= 32) { const hl = b[8] | (b[9] << 8), names = []; for (let o = 32; o + 32 <= Math.min(hl, b.length) && b[o] !== 0x0d; o += 32) names.push(latin1.decode(b.subarray(o, o + 11)).replace(/\0.*$/, '').trim()); stats.attributes = names; stats.attributeRecords = (b[4] | (b[5] << 8) | (b[6] << 16) | (b[7] << 24)) >>> 0; }
  }
  let geographic;
  if (prj) { const t = utf8.decode(new Uint8Array(await prj.arrayBuffer())).trim(); stats.projection = t.slice(0, 160); geographic = /^GEOGCS|^GEOGCRS/i.test(t); }
  const g = linesOrPoints(polys, pts, { warnings, stats });
  g.geographic = geographic ?? (g.bbox.min[0] >= -180 && g.bbox.max[0] <= 180 && g.bbox.min[1] >= -90 && g.bbox.max[1] <= 90);
  if (!prj) warnings.push(`No .prj companion supplied: coordinates are assumed ${g.geographic ? 'geographic (degrees), judging by their range' : 'projected'}.`);
  return g;
}

async function readMIF(ctx) {
  const text = await ctx.text(), di = text.search(/^\s*DATA\b/im), tk = (di >= 0 ? text.slice(di).replace(/^\s*DATA/i, '') : text).match(/"[^"]*"|\([^)]*\)|[^\s,]+/g) || [], polys = [], pts = [];
  let i = 0;
  const num = () => +tk[i++];
  const ring = (closed) => {
    const n = num(), x = [], y = [];
    if (!(n >= 0) || i + 2 * n > tk.length) fail('The MIF file is truncated or malformed.');
    for (let k = 0; k < n; k++) { x.push(num()); y.push(num()); }
    if (closed) closeRing(x, y);
    if (x.length > 1) polys.push({ x, y, closed });
  };
  while (i < tk.length) {
    const w = tk[i++].toUpperCase();
    if (w === 'REGION') { const n = num(); if (!(n >= 0 && n <= tk.length)) fail('The MIF file is malformed.'); for (let k = 0; k < n; k++) ring(true); }
    else if (w === 'PLINE') { if (String(tk[i]).toUpperCase() === 'MULTIPLE') { i++; const n = num(); if (!(n >= 0 && n <= tk.length)) fail('The MIF file is malformed.'); for (let k = 0; k < n; k++) ring(false); } else ring(false); }
    else if (w === 'LINE' && i + 4 <= tk.length) { const v = [num(), num(), num(), num()]; polys.push({ x: [v[0], v[2]], y: [v[1], v[3]], closed: false }); }
    else if ((w === 'RECT' || w === 'ROUNDRECT') && i + 4 <= tk.length) { const v = [num(), num(), num(), num()]; polys.push({ x: [v[0], v[2], v[2], v[0]], y: [v[1], v[1], v[3], v[3]], closed: true }); }
    else if (w === 'POINT' && i + 2 <= tk.length) pts.push(num(), num(), 0);
  }
  return linesOrPoints(polys, pts, {});
}

// ---- Kernel B-rep files: ACIS SAT and Parasolid X_T (neutral model of fmt_brep.js) ---------------------------------------
const BREP_UNIT = { 1: 'm', 0.001: 'mm', 0.01: 'cm', 0.0254: 'in', 0.3048: 'ft' };
/** Parametric rectangle [u0, u0 + du] × [v0, v0 + dv] of a surface P(u, v) whose ∂u × ∂v points outwards -> triangles. */
function patchFace(P, u0, du, v0, dv, sense, tri) {
  const nu = Math.max(2, Math.ceil((Math.abs(du) / TAU) * ARC_N - 1e-9)), nv = Math.max(2, Math.ceil((Math.abs(dv) / TAU) * ARC_N - 1e-9)), g = [];
  for (let j = 0; j <= nv; j++) for (let i = 0; i <= nu; i++) g.push(P(u0 + (du * i) / nu, v0 + (dv * j) / nv));
  const eq = (a, b) => a[0] === b[0] && a[1] === b[1] && a[2] === b[2], put = (a, b, c) => { if (eq(a, b) || eq(b, c) || eq(a, c)) return; for (const q of sense ? [a, b, c] : [a, c, b]) tri.push(q[0], q[1], q[2]); };
  for (let j = 0; j < nv; j++) for (let i = 0; i < nu; i++) { const a = g[j * (nu + 1) + i], b = g[j * (nu + 1) + i + 1], c = g[(j + 1) * (nu + 1) + i + 1], d = g[(j + 1) * (nu + 1) + i]; put(a, b, c); put(a, c, d); }
}
/**
 * B-spline face: the surface is filled over the parametric rectangle spanned by its boundary (each boundary point is
 * projected by a grid search refined by halving steps). Returns false when the boundary does not lie on the surface.
 */
function splineFace(f, tri) {
  const s = f.surf.nurbs, U0 = s.knotsU[s.degU], U1 = s.knotsU[s.nU], V0 = s.knotsV[s.degV], V1 = s.knotsV[s.nV], pts = f.loops.flatMap((l) => l.pts);
  let a = U0, b = U1, c = V0, d = V1;
  if (pts.length) {
    const G = 12, grid = [], step = Math.max(1, Math.ceil(pts.length / 24)), dist = (p, q) => Math.hypot(p[0] - q[0], p[1] - q[1], p[2] - q[2]);
    for (let j = 0; j <= G; j++) for (let i = 0; i <= G; i++) { const u = U0 + ((U1 - U0) * i) / G, v = V0 + ((V1 - V0) * j) / G; grid.push([u, v, nurbsSurfPoint(s, u, v)]); }
    let lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity], worst = 0;
    for (const q of pts) for (let k = 0; k < 3; k++) { lo[k] = Math.min(lo[k], q[k]); hi[k] = Math.max(hi[k], q[k]); }
    a = c = Infinity; b = d = -Infinity;
    for (let n = 0; n < pts.length; n += step) {
      const q = pts[n];
      let best = Infinity, bu = U0, bv = V0, su = (U1 - U0) / G, sv = (V1 - V0) / G;
      for (const g of grid) { const e = dist(g[2], q); if (e < best) { best = e; bu = g[0]; bv = g[1]; } }
      for (let it = 0; it < 60 && (su > 1e-9 * (U1 - U0) || sv > 1e-9 * (V1 - V0)); it++) {
        let moved = false;
        for (const [du, dv] of [[su, 0], [-su, 0], [0, sv], [0, -sv]]) { const u = Math.min(U1, Math.max(U0, bu + du)), v = Math.min(V1, Math.max(V0, bv + dv)), e = dist(nurbsSurfPoint(s, u, v), q); if (e < best) { best = e; bu = u; bv = v; moved = true; } }
        if (!moved) { su /= 2; sv /= 2; }
      }
      worst = Math.max(worst, best); a = Math.min(a, bu); b = Math.max(b, bu); c = Math.min(c, bv); d = Math.max(d, bv);
    }
    if (worst > 1e-3 * (dist(lo, hi) || 1)) return false;
  }
  if (!(b > a) || !(d > c)) return false;
  // one strip per knot span inside the rectangle, refined for curved (degree > 1) directions
  const cuts = (deg, knots, lo, hi) => { const ks = [lo, ...knots.filter((k, i) => k > lo && k < hi && k !== knots[i - 1]), hi], out = [], n = deg > 1 ? Math.max(2, Math.min(8, Math.ceil(96 / ks.length))) : 1; for (let i = 0; i + 1 < ks.length; i++) for (let j = 0; j < n; j++) out.push(ks[i] + ((ks[i + 1] - ks[i]) * j) / n); out.push(hi); return out.length > 400 ? out.filter((_, i) => i % Math.ceil(out.length / 400) === 0 || i === out.length - 1) : out; };
  const us = cuts(s.degU, s.knotsU, a, b), vs = cuts(s.degV, s.knotsV, c, d), g = [];
  for (const v of vs) for (const u of us) g.push(nurbsSurfPoint(s, u, v));
  const nu = us.length, eq = (p, q) => p[0] === q[0] && p[1] === q[1] && p[2] === q[2], put = (p, q, r) => { if (eq(p, q) || eq(q, r) || eq(p, r)) return; for (const w of f.sense ? [p, q, r] : [p, r, q]) tri.push(w[0], w[1], w[2]); };
  for (let j = 0; j + 1 < vs.length; j++) for (let i = 0; i + 1 < nu; i++) { const p = g[j * nu + i], q = g[j * nu + i + 1], r = g[(j + 1) * nu + i + 1], t = g[(j + 1) * nu + i]; put(p, q, r); put(p, r, t); }
  return true;
}
/**
 * Neutral B-rep model -> Geometry. A model recognised as a pipe gives its centreline (3-D polyline) with g.pipe = { runs,
 * diameters, length }; any other model gives the tessellation of its analytic faces, or its wireframe when no face can be
 * tessellated. opts.prefer = 'centreline' | 'mesh' | 'wireframe' overrides the choice.
 */
function brepGeom(m, opts) {
  const warnings = m.warnings.slice(), tri = [], pr = recognisePipes(m), left = {}, units = BREP_UNIT[m.unitScale], prefer = String(opts.prefer || '');
  let done = 0, nSpline = 0, spl = 0;
  for (const f of m.faces) {
    const s = f.surf, before = tri.length, pl = s.o ? { o: s.o, x: s.x, y: cross(s.z, s.x), z: s.z } : null, pts = f.loops.flatMap((l) => l.pts);
    let why = s.type === 'spline' ? 'free-form (spline)' : s.type;
    if (s.type === 'spline' && s.nurbs) { if (nSpline++ < 20000 && splineFace(f, tri)) spl++; else why = 'B-spline surface whose boundary could not be located'; }
    else if (s.type === 'plane') planarFace(f.loops.filter((l) => !l.single).map((l) => ({ pts: l.pts })), pl, f.sense, tri);
    else if ((s.type === 'cylinder' || s.type === 'cone') && s.ratio === undefined) revolvedFace(f.loops, pl, s.r, s.tanA, f.sense, tri);
    else if (s.type === 'cylinder' || s.type === 'cone') why = 'elliptical ' + s.type;
    else if (s.type === 'torus') {
      // the face is taken as the parametric rectangle spanned by its boundary: exact for elbows and untrimmed patches
      const th = [], ph = [];
      for (const q of pts) { const d = sub(q, s.o), u = dot(d, pl.x), w = dot(d, pl.y); th.push(Math.atan2(w, u)); ph.push(Math.atan2(dot(d, s.z), Math.hypot(u, w) - s.R)); }
      const [t0, dt] = angularExtent(th), [p0, dp] = angularExtent(ph);
      patchFace((t, q) => { const c = s.R + s.r * Math.cos(q), a = c * Math.cos(t), b = c * Math.sin(t), h = s.r * Math.sin(q); return [s.o[0] + a * pl.x[0] + b * pl.y[0] + h * s.z[0], s.o[1] + a * pl.x[1] + b * pl.y[1] + h * s.z[1], s.o[2] + a * pl.x[2] + b * pl.y[2] + h * s.z[2]]; }, t0, dt, p0, dp, f.sense, tri);
    } else if (s.type === 'sphere') {
      const lat = pts.map((q) => Math.asin(Math.max(-1, Math.min(1, dot(sub(q, s.o), s.z) / s.r)))), band = f.loops.length >= 2 || !f.loops.length;
      if (band) { const a = f.loops.length ? Math.min(...lat) : -Math.PI / 2, b = f.loops.length ? Math.max(...lat) : Math.PI / 2; patchFace((t, q) => { const c = s.r * Math.cos(q), a2 = c * Math.cos(t), b2 = c * Math.sin(t), h = s.r * Math.sin(q); return [s.o[0] + a2 * pl.x[0] + b2 * pl.y[0] + h * s.z[0], s.o[1] + a2 * pl.x[1] + b2 * pl.y[1] + h * s.z[1], s.o[2] + a2 * pl.x[2] + b2 * pl.y[2] + h * s.z[2]]; }, 0, TAU, a, b - a, f.sense, tri); }
      else why = 'trimmed sphere';
    }
    if (tri.length > before) done++; else left[why] = (left[why] || 0) + 1;
    if (tri.length > L.triangles * 9) fail(`The model tessellates to too many triangles (limit ${L.triangles}).`);
  }
  const runs = pr.runs.slice(0, 500).map((r) => ({ kind: r.kind, length: r.length, diameters: r.radii.map((q) => 2 * q), ...(r.kind === 'bend' ? { bendRadius: r.bendRadius, angle: (r.angle * 180) / Math.PI } : {}), from: r.a, to: r.b }));
  const stats = { version: m.version, ...(m.schema ? { schema: m.schema } : {}), ...(units ? { units } : {}), faces: m.faces.length, tessellatedFaces: done, edges: m.edges.length, surfaces: m.counts.surfaces, pipeRuns: runs, diameters: pr.diameters, ...(pr.centrelines.length ? { centrelineLength: pr.length } : {}) };
  const nLeft = m.faces.length - done;
  if (m.unitScale === null || m.unitScale === undefined) warnings.push('The file gives no length unit; coordinates are returned as written.');
  else if (!units) { stats.unitScale = m.unitScale; warnings.push(`One model unit is ${m.unitScale} m; coordinates are returned in model units.`); }
  const poly = (list) => list.map((c) => ({ x: c.pts.map((q) => q[0]), y: c.pts.map((q) => q[1]), z: c.pts.map((q) => q[2]), closed: !!c.closed }));
  const want = prefer === 'centreline' && pr.centrelines.length ? 'centreline' : prefer === 'mesh' && tri.length ? 'mesh' : prefer === 'wireframe' && m.edges.length ? 'wireframe' : pr.isPipe ? 'centreline' : tri.length ? 'mesh' : m.edges.length ? 'wireframe' : pr.centrelines.length ? 'centreline' : '';
  const pipe = { isPipe: pr.isPipe, runs, diameters: pr.diameters, length: pr.length };
  if (want === 'centreline') {
    const main = pr.centrelines.filter((c, k) => k === 0 || c.length >= 0.02 * pr.length);
    warnings.push(`${pr.isPipe ? 'Recognised as a pipe' : 'Pipe-like surfaces found'}: ${pr.runs.filter((r) => r.kind === 'straight').length} straight runs, ${pr.runs.filter((r) => r.kind === 'bend').length} bends and ${pr.runs.filter((r) => r.kind === 'reducer').length} reducers give a centreline of ${+pr.length.toPrecision(6)} model units with diameters ${pr.diameters.join(', ')}. Pass opts.prefer = "mesh" for the tessellated faces.`);
    return Object.assign(polyGeom(poly(main), { warnings, stats: { ...stats, representation: 'pipe centreline' } }), { pipe });
  }
  if (want === 'mesh') {
    if (nLeft) warnings.push(`${done} of ${m.faces.length} faces tessellated; ${nLeft} left out (${Object.entries(left).map(([k, v]) => `${v} × ${k}`).join(', ')}). Export STEP AP242 with tessellation or STL for a complete surface.`);
    if (m.counts.surfaces.torus || m.counts.surfaces.sphere || spl) warnings.push(`Toroidal, spherical${spl ? ' and B-spline' : ''} faces are filled over the parametric rectangle of their boundary (exact for elbows and untrimmed patches${spl ? `; ${spl} B-spline faces` : ''}).`);
    return Object.assign(meshGeom(tri, { warnings, stats: { ...stats, representation: 'tessellated faces' } }), pr.runs.length ? { pipe } : {});
  }
  if (want !== 'wireframe') fail('The model holds no face or edge geometry that can be shown.');
  warnings.push(`No face could be tessellated (${Object.entries(left).map(([k, v]) => `${v} × ${k}`).join(', ') || 'no faces'}); the edge wireframe is returned.`);
  return Object.assign(polyGeom(poly(m.edges), { warnings, stats: { ...stats, representation: 'wireframe' } }), pr.runs.length ? { pipe } : {});
}
async function readSAT(ctx) {
  const u8 = await ctx.bytes();
  if (isSAB(u8)) return brepGeom(parseSAB(u8), ctx.opts);
  if (ctx.ext === 'sab') fail('The .sab file does not start with an ACIS binary header ("ACIS BinaryFile" or "ASM BinaryFile").');
  return brepGeom(parseSAT(latin1.decode(u8)), ctx.opts);
}

async function readXT(ctx) {
  const u8 = await ctx.bytes(), kind = xtKind(u8);
  if (kind === 'binary' || kind === 'neutral') return brepGeom(parseXB(u8), ctx.opts);
  if (!kind) fail('Not a Parasolid transmit file (the "**PARASOLID" keyword header is missing).');
  return brepGeom(parseXT(latin1.decode(u8)), ctx.opts);
}

// ---- AutoCAD DWG (fmt_dwg.js) -------------------------------------------------------------------------------------------
const DWG_UNIT = { 1: 'in', 2: 'ft', 4: 'mm', 5: 'cm', 6: 'm', 7: 'km', 10: 'yd', 14: 'dm' };
async function readDWG(ctx) {
  const d = parseDWG(await ctx.bytes(), { space: ctx.opts.space }), warnings = d.warnings, polys = [];
  for (const p of d.polylines) {
    if (!p.spline) { polys.push(p); continue; }
    const s = p.spline, smp = s.points.length ? nurbsCurve(s.degree, s.knots, s.points, s.weights, NaN, NaN, 8 * s.points.length) : null, pts = smp || (s.fit.length > 1 ? s.fit : s.points);
    if (pts.length > 1) polys.push({ x: pts.map((q) => q[0]), y: pts.map((q) => q[1]), ...(pts.some((q) => q[2] !== 0) ? { z: pts.map((q) => q[2]) } : {}), closed: !smp && s.closed, layer: p.layer, type: 'SPLINE' });
  }
  const sk = Object.entries(d.skipped), unit = DWG_UNIT[d.insunits], labels = d.labels.slice(0, 5000);
  const stats = { version: d.version, release: `AutoCAD ${d.release}`, entities: d.counts, layers: d.layers.slice(0, 200), ...(unit ? { units: unit } : {}), ...(d.insunits !== null ? { insunits: d.insunits } : {}), ...(d.labels.length ? { labels: d.labels.length } : {}), ...(d.paperSpace ? { paperSpaceEntities: d.paperSpace } : {}), ...(d.hidden ? { invisibleEntities: d.hidden } : {}) };
  const offLayers = d.layerInfo.filter((l) => l.off || l.frozen).map((l) => l.name);
  if (offLayers.length) { stats.layersOff = offLayers.slice(0, 200); warnings.push(`Layers switched off or frozen are still read: ${offLayers.slice(0, 8).join(', ')}${offLayers.length > 8 ? ' …' : ''}.`); }
  if (sk.length) warnings.push(`Entities that carry no line work were skipped: ${sk.slice(0, 12).map(([k, v]) => `${v} × ${k}`).join(', ')}.`);
  if (d.insunits === null) warnings.push(/^AC101[24]$/.test(d.version) ? 'R13 and R14 drawings carry no drawing unit ($INSUNITS came with AutoCAD 2000); coordinates are returned as drawn.' : 'The drawing unit ($INSUNITS) could not be read; coordinates are returned as drawn.');
  else if (!unit) warnings.push('The drawing is unitless ($INSUNITS = 0 or an unusual unit); coordinates are returned as drawn.');
  if (d.paperSpace && (ctx.opts.space || 'model') === 'model') warnings.push(`${d.paperSpace} paper-space entities (sheet layouts) were left out; pass opts.space = "paper" or "all" to read them.`);
  if (d.hidden) warnings.push(`${d.hidden} invisible entities were left out.`);
  // embedded ACIS bodies (3DSOLID, REGION, BODY): read through the SAT / SAB reader, tessellated in drawing coordinates
  const solidTri = [], solidPolys = [];
  let nSolid = 0, badSolid = 0, pipe = null;
  for (const so of d.solids) {
    try {
      const m = so.format === 'sab' ? parseSAB(so.data) : parseSAT(so.data), g = brepGeom(m, { prefer: ctx.opts.prefer || 'mesh' });
      if (g.kind === 'mesh') for (const v of g.triangles) solidTri.push(v); else for (const q of g.polylines) solidPolys.push({ ...q, layer: so.layer, type: so.type });
      if (g.pipe && g.pipe.isPipe && !pipe) pipe = g.pipe;
      nSolid++;
    } catch (e) { if (!(e && e.user) && !(e instanceof RangeError) && !(e instanceof TypeError)) throw e; badSolid++; }
    if (solidTri.length > L.triangles * 9) fail(`The solids of the drawing tessellate to too many triangles (limit ${L.triangles}).`);
  }
  if (nSolid) stats.solids = nSolid;
  if (badSolid) warnings.push(`${badSolid} embedded ACIS bodies could not be read.`);
  const wantSolids = String(ctx.opts.prefer || '') === 'solids' || polys.length + d.points.length / 3 < nSolid;   // a drawing of solids with a stray line or two is a 3-D model
  let g;
  if (solidTri.length && wantSolids) g = Object.assign(meshGeom(solidTri, { warnings, stats: { ...stats, representation: 'tessellated ACIS solids' } }), pipe ? { pipe } : {});
  else if (!polys.length && !d.points.length && !solidPolys.length && d.faces.length) g = meshGeom(d.faces, { warnings, stats: { ...stats, faces3d: d.counts['3DFACE'] } });
  else {
    if (solidTri.length) warnings.push(`${nSolid} ACIS solids were read but left out in favour of the line work; pass opts.prefer = "solids" for their tessellation.`);
    if (d.faces.length) warnings.push(`${d.counts['3DFACE']} 3DFACE entities were ignored in favour of the line work.`);
    if (!polys.length && !d.points.length && !solidPolys.length) fail(`No LINE, POLYLINE, CIRCLE, ARC, ELLIPSE, SPLINE, POINT, 3DFACE or readable 3DSOLID entities were found in the ${ctx.opts.space === 'paper' ? 'paper' : 'model'} space of this AutoCAD ${d.release} drawing${sk.length ? ` (it holds ${sk.slice(0, 6).map(([k, v]) => `${v} × ${k}`).join(', ')})` : ''}${labels.length ? `; it carries ${d.labels.length} text labels` : ''}.`);
    g = linesOrPoints(polys.length || d.points.length ? polys : solidPolys, d.points, { warnings, stats });
  }
  if (labels.length) g.labels = labels;
  return g;
}

// ---- MicroStation DGN V7 (fmt_dgn.js) ------------------------------------------------------------------------------
async function readDGN(ctx) {
  const d = parseDGN(await ctx.bytes()), names = Object.fromEntries(Object.entries(d.counts).map(([t, n]) => [({ 2: 'cell', 3: 'line', 4: 'line string', 6: 'shape', 11: 'curve', 12: 'complex chain', 14: 'complex shape', 15: 'ellipse', 16: 'arc', 21: 'B-spline poles', 22: 'point string' })[t] || `type ${t}`, n]));
  const sk = Object.entries(d.skipped).filter(([k]) => k !== 'deleted'), unit = d.units.master.toLowerCase();
  if (sk.length) d.warnings.push(`Elements that carry no line work were skipped: ${sk.map(([k, v]) => `${v} × ${k}`).join(', ')}.`);
  return linesOrPoints(d.polylines, d.points, { warnings: d.warnings, stats: { version: 'V7', dimension: d.is3d ? 3 : 2, ...(unit ? { units: unit } : {}), subUnit: d.units.sub || undefined, elements: names, levels: d.levels } });
}

// ---- GeoPackage (SQLite container, GeoPackage binary + WKB geometries) ------------------------------------------------
/**
 * One WKB geometry (OGC, ISO Z / M / ZM codes or EWKB flags, either byte order) starting at b[p0] -> out.polys / out.pts.
 * Curve types are flattened to their control points (counted in out.curves). Returns false when the blob is malformed.
 */
function wkbGeom(b, p0, out, extra) {
  const dv = new DataView(b.buffer, b.byteOffset, b.byteLength), n = b.length;
  let p = p0;
  const put = (s, mode) => {
    if (mode === 2) return s;
    const m = s.x.length;
    if (m === 1) out.pts.push(s.x[0], s.y[0], s.z ? s.z[0] : 0);
    if (m < 2) return null;
    if (mode === 1 && m > 2 && s.x[0] === s.x[m - 1] && s.y[0] === s.y[m - 1]) { s.x.pop(); s.y.pop(); if (s.z) s.z.pop(); }
    out.vertices += m;
    out.polys.push({ x: s.x, y: s.y, ...(s.z ? { z: s.z } : {}), closed: mode === 1, ...extra });
    return null;
  };
  const geom = (depth, mode) => {
    if (depth > 16 || p + 5 > n || b[p] > 1) throw 0;
    const le = b[p] === 1;
    let t = dv.getUint32(p + 1, le), hasZ = false, hasM = false;
    p += 5;
    if (t & 0x80000000) hasZ = true;
    if (t & 0x40000000) hasM = true;
    if (t & 0x20000000) p += 4;
    t &= 0x0fffffff;
    if (t >= 3000) { hasZ = hasM = true; t -= 3000; } else if (t >= 2000) { hasM = true; t -= 2000; } else if (t >= 1000) { hasZ = true; t -= 1000; }
    const dim = 2 + (hasZ ? 1 : 0) + (hasM ? 1 : 0);
    const count = (unit) => { if (p + 4 > n) throw 0; const c = dv.getUint32(p, le); p += 4; if (c * unit > n - p) throw 0; return c; };
    const seq = () => {
      const c = count(8 * dim), x = new Array(c), y = new Array(c), z = hasZ ? new Array(c) : null;
      for (let k = 0; k < c; k++, p += 8 * dim) { x[k] = dv.getFloat64(p, le); y[k] = dv.getFloat64(p + 8, le); if (z) z[k] = dv.getFloat64(p + 16, le); }
      return { x, y, z };
    };
    if (t === 1) {
      if (p + 8 * dim > n) throw 0;
      const x = dv.getFloat64(p, le), y = dv.getFloat64(p + 8, le), z = hasZ ? dv.getFloat64(p + 16, le) : 0;
      p += 8 * dim;
      if (x === x && y === y) { if (mode === 2) return { x: [x], y: [y], z: hasZ ? [z] : null }; out.pts.push(x, y, z === z ? z : 0); }
      return null;
    }
    if (t === 2 || t === 8) { if (t === 8) out.curves++; return put(seq(), mode); }
    if (t === 3 || t === 17) { for (let k = 0, c = count(4); k < c; k++) put(seq(), 1); return null; }
    if (t === 9) {
      const acc = { x: [], y: [], z: hasZ ? [] : null };
      for (let k = 0, c = count(5); k < c; k++) { const s = geom(depth + 1, 2); if (!s) continue; const skip = acc.x.length && s.x[0] === acc.x[acc.x.length - 1] && s.y[0] === acc.y[acc.y.length - 1] ? 1 : 0; for (let i = skip; i < s.x.length; i++) { acc.x.push(s.x[i]); acc.y.push(s.y[i]); if (acc.z) acc.z.push(s.z ? s.z[i] : 0); } }
      return put(acc, mode);
    }
    if (t === 10) { for (let k = 0, c = count(5); k < c; k++) geom(depth + 1, 1); return null; }
    if ([4, 5, 6, 7, 11, 12, 15, 16].includes(t)) { for (let k = 0, c = count(5); k < c; k++) geom(depth + 1, mode === 2 ? 0 : mode); return null; }
    throw 0;
  };
  try { geom(0, 0); return true; } catch (e) { if (e === 0) return false; throw e; }
}
const GPKG_ENV = [0, 32, 48, 48, 64];
async function readGPKG(ctx) {
  const db = openSQLite(await ctx.bytes()), warnings = [], opts = ctx.opts;
  if (!db.tables.has('gpkg_contents')) fail('The file is an SQLite database but not a GeoPackage (table gpkg_contents is missing).');
  if (db.wal) warnings.push('The GeoPackage is in write-ahead-log mode; changes still held in a "-wal" file beside it are not read.');
  const contents = db.rows('gpkg_contents', { max: 10000 }), gcols = db.rows('gpkg_geometry_columns', { max: 10000 }).filter((r) => db.tables.has(String(r.table_name)));
  const srsOf = new Map(db.rows('gpkg_spatial_ref_sys', { max: 100000 }).map((r) => [r.srs_id, r])), rasters = contents.filter((r) => /tiles|gridded/i.test(String(r.data_type))).map((r) => String(r.table_name));
  const rasterNote = rasters.length ? `Tile pyramids and gridded-coverage (elevation) rasters are not read (${rasters.slice(0, 4).join(', ')}); export them as GeoTIFF: gdal_translate in.gpkg out.tif.` : '';
  if (!gcols.length) fail(`The GeoPackage holds no vector feature table. ${rasterNote || 'Its attribute tables carry no geometry.'}`);
  const want = opts.layer !== undefined ? gcols.filter((r) => String(r.table_name) === String(opts.layer)) : gcols;
  if (!want.length) fail(`The GeoPackage has no feature layer "${String(opts.layer).slice(0, 60)}" (it holds ${gcols.slice(0, 10).map((r) => r.table_name).join(', ')}).`);
  if (rasterNote) warnings.push(rasterNote);
  const out = { polys: [], pts: [], vertices: 0, curves: 0 }, layers = [], attributes = [];
  let nullGeom = 0, badGeom = 0, total = 0, geographic, cut = false;
  for (const lay of want.slice(0, 64)) {
    const name = String(lay.table_name), gc = String(lay.column_name), t = db.tables.get(name), cols = t.cols.filter((c) => c !== gc).slice(0, 32), nameCol = cols.find((c) => /^(name|label|title|line_?name|pipeline|route|tag|ref|descr(iption)?)$/i.test(c)), srs = srsOf.get(lay.srs_id);
    let nf = 0;
    if (srs && geographic === undefined && !/^\s*undefined\s*$/i.test(String(srs.definition || 'undefined'))) geographic = /^\s*GEOG(CS|CRS|RAPHICCRS)\b/i.test(String(srs.definition)) || (/^epsg$/i.test(String(srs.organization || '')) && +srs.organization_coordsys_id === 4326);
    for (const r of db.rows(name, { max: 2e6 })) {
      const g = r[gc];
      nf++;
      if (!(g instanceof Uint8Array) || g.length < 8) { nullGeom++; continue; }
      if (g[0] !== 0x47 || g[1] !== 0x50) { badGeom++; continue; }
      const flags = g[3], env = GPKG_ENV[(flags >> 1) & 7];
      if (flags & 0x10) { nullGeom++; continue; }           // empty geometry
      if (env === undefined || 8 + env + 5 > g.length) { badGeom++; continue; }
      const attrs = {};
      for (const c of cols) { const v = r[c]; if (v !== null && v !== undefined && !(v instanceof Uint8Array)) attrs[c] = typeof v === 'string' ? v.slice(0, 200) : v; }
      const np = out.pts.length, tag = { layer: name, ...(nameCol && r[nameCol] !== null ? { name: String(r[nameCol]).slice(0, 120) } : {}), ...(Object.keys(attrs).length ? { attrs } : {}) };
      if (!wkbGeom(g, 8 + env, out, tag)) badGeom++;
      for (let k = np; k < out.pts.length && attributes.length < L.rows; k += 3) attributes.push(attrs);
      if (out.vertices > 5e6 || out.pts.length > 3 * 8e6) { cut = true; break; }
    }
    total += nf;
    layers.push({ name, type: String(lay.geometry_type_name || ''), features: nf, srs: srs ? String(srs.srs_name || lay.srs_id).slice(0, 80) : lay.srs_id, hasZ: lay.z === 1 || lay.z === 2 });
    if (cut) { warnings.push('The GeoPackage holds more vertices than can be shown; the remaining features were not read.'); break; }
  }
  if (nullGeom) warnings.push(`${nullGeom} features without geometry were skipped.`);
  if (badGeom) warnings.push(`${badGeom} geometries are not valid GeoPackage binary / WKB and were skipped.`);
  if (out.curves) warnings.push(`${out.curves} circular-arc strings were reduced to their control points.`);
  if (gcols.length > want.length) warnings.push(`Only layer "${want[0].table_name}" of ${gcols.length} feature layers was read.`);
  else if (layers.length > 1) warnings.push(`${layers.length} feature layers were merged (${layers.map((l) => l.name).slice(0, 8).join(', ')}); pass opts.layer to read one.`);
  if (!out.polys.length && !out.pts.length) fail(`The GeoPackage feature ${layers.length > 1 ? 'layers hold' : `layer "${layers[0].name}" holds`} no readable geometry (${total} features${nullGeom ? `, ${nullGeom} without geometry` : ''}${badGeom ? `, ${badGeom} not valid GeoPackage binary` : ''}).`);
  const g = linesOrPoints(out.polys, out.pts, { warnings, stats: { layers, features: total } });
  if (g.kind === 'points' && attributes.length === g.count && attributes.some((a) => Object.keys(a).length)) g.attributes = attributes;
  g.geographic = geographic ?? (g.bbox.min[0] >= -180 && g.bbox.max[0] <= 180 && g.bbox.min[1] >= -90 && g.bbox.max[1] <= 90);
  return g;
}

// ---- Rasters: ESRI ASCII, Surfer, TIFF, NetCDF classic, band-interleaved ----------------------------------
function readASCGrid(text) {
  const h = Object.create(null);
  let pos = 0, m;
  const re = /\s*([A-Za-z_]+)[ \t]+([-+.\deE]+)[ \t]*(?:\r\n|\n|\r)/y;
  while ((m = re.exec(text))) { h[m[1].toLowerCase()] = +m[2]; pos = re.lastIndex; }
  const nx = h.ncols, ny = h.nrows, dx = h.cellsize ?? h.dx, dy = h.cellsize ?? h.dy;
  if (!Number.isInteger(nx) || !Number.isInteger(ny) || nx < 1 || ny < 1 || !(dx > 0) || !(dy > 0)) fail('The ESRI ASCII grid header is incomplete (ncols, nrows, cellsize).');
  if (nx * ny > L.voxels) fail('The raster is too large.');
  const x0 = h.xllcorner ?? (h.xllcenter ?? 0) - dx / 2, y0 = h.yllcorner ?? (h.yllcenter ?? 0) - dy / 2, v = numsOf(text.slice(pos));
  if (v.length < nx * ny) fail(`The ESRI ASCII grid holds ${v.length} values but its header declares ${nx * ny}.`);
  const geo = x0 >= -360 && x0 + nx * dx <= 360.0001 && y0 >= -90.0001 && y0 + ny * dy <= 90.0001 && dx < 1;
  return gridGeom(v, nx, ny, (i) => x0 + (i + 0.5) * dx, (j) => y0 + (ny - 1 - j + 0.5) * dy, { nodata: h.nodata_value ?? -9999, geographic: geo, stats: { cellsize: [dx, dy] } });
}
async function readSurfer(ctx) {
  const u8 = await ctx.bytes(), sig = latin1.decode(u8.subarray(0, 4)), blank = (v) => v >= 1.70141e38, lin = (lo, hi, n) => (i) => (n > 1 ? lo + ((hi - lo) * i) / (n - 1) : lo);
  if (sig === 'CDF\x01' || sig === 'CDF\x02' || sig === 'CDF\x05' || sig === '\x89HDF') return readNetCDF(ctx);
  if (/^ncols/i.test(latin1.decode(u8.subarray(0, 64)).trim())) return Object.assign(readASCGrid(await ctx.text()), { format: 'ESRI ASCII grid' });
  if (sig === 'DSAA') {
    const v = numsOf((await ctx.text()).replace(/^\s*DSAA/, '')), nx = v[0], ny = v[1];
    if (!Number.isInteger(nx) || !Number.isInteger(ny) || nx < 1 || ny < 1 || v.length < 8 + nx * ny) fail('The Surfer ASCII grid is truncated or malformed.');
    return gridGeom(v.subarray(8), nx, ny, lin(v[2], v[3], nx), lin(v[4], v[5], ny), { nodata: blank });
  }
  const dv = new DataView(u8.buffer, u8.byteOffset, u8.byteLength);
  if (sig === 'DSBB') {
    const nx = dv.getInt16(4, true), ny = dv.getInt16(6, true), f = (o) => dv.getFloat64(o, true);
    if (nx < 1 || ny < 1) fail('The Surfer grid has invalid dimensions.');
    return gridGeom(typed(u8, 56, nx * ny, 'f4'), nx, ny, lin(f(8), f(16), nx), lin(f(24), f(32), ny), { nodata: blank });
  }
  if (sig === 'DSRB') {
    let p = 0, hd = null;
    while (p + 8 <= u8.length) {
      const id = latin1.decode(u8.subarray(p, p + 4)), size = dv.getInt32(p + 4, true), s = p + 8, f = (o) => dv.getFloat64(s + o, true);
      if (size < 0 || s + size > u8.length) break;
      if (id === 'GRID' && size >= 72) hd = { ny: dv.getInt32(s, true), nx: dv.getInt32(s + 4, true), x0: f(8), y0: f(16), dx: f(24), dy: f(32), blank: f(64) };
      else if (id === 'DATA' && hd) { if (hd.nx < 1 || hd.ny < 1 || hd.nx * hd.ny > L.voxels) fail('The Surfer grid has invalid dimensions.'); return gridGeom(typed(u8, s, hd.nx * hd.ny, 'f8'), hd.nx, hd.ny, (i) => hd.x0 + i * hd.dx, (j) => hd.y0 + j * hd.dy, { nodata: (v) => v >= hd.blank || v >= 1.70141e38 }); }
      p = s + size;
    }
    fail('The Surfer 7 grid has no GRID / DATA sections.');
  }
  return fail('Not a Surfer grid (DSAA, DSBB or DSRB signature missing).');
}

function tiffLZW(src, expected) {
  const out = new Uint8Array(expected), prefix = new Int32Array(4097), suffix = new Uint8Array(4097), len = new Int32Array(4097), firstc = new Uint8Array(4097), nbits = src.length * 8;
  for (let k = 0; k < 256; k++) { len[k] = 1; firstc[k] = k; suffix[k] = k; }
  let o = 0, bitPos = 0, width = 9, next = 258, prev = -1;
  while (bitPos + width <= nbits && o < expected) {
    const bp = bitPos >> 3, v = (src[bp] << 16) | ((src[bp + 1] || 0) << 8) | (src[bp + 2] || 0), code = (v >> (24 - (bitPos & 7) - width)) & ((1 << width) - 1);
    bitPos += width;
    if (code === 257) break;
    if (code === 256) { width = 9; next = 258; prev = -1; continue; }
    if (prev < 0) { if (code > 255) break; out[o++] = code; prev = code; continue; }
    if (code > next || (code >= 256 && code < 258)) break;
    if (next < 4096) { prefix[next] = prev; firstc[next] = firstc[prev]; suffix[next] = code === next ? firstc[prev] : firstc[code]; len[next] = len[prev] + 1; }
    else if (code === next) break;
    const n = len[code];
    for (let c = code, q = o + n - 1; q >= o; q--) { if (q < expected) out[q] = suffix[c]; c = prefix[c]; }
    o += n;
    if (next < 4096) next++;
    if (next >= (1 << width) - 1 && width < 12) width++;
    prev = code;
  }
  return out;
}
function packBits(src, expected) {
  const out = new Uint8Array(expected);
  let o = 0, p = 0;
  while (p < src.length && o < expected) {
    const n = src[p++] << 24 >> 24;
    if (n >= 0) { for (let k = 0; k <= n && p < src.length && o < expected; k++) out[o++] = src[p++]; }
    else if (n !== -128) { const b = src[p++]; for (let k = 0; k < 1 - n && o < expected; k++) out[o++] = b; }
  }
  return out;
}
const TIFF_SZ = [0, 1, 1, 2, 4, 8, 1, 1, 2, 4, 8, 4, 8];
const TIFF_COMP = { 2: 'CCITT', 3: 'CCITT Group 3', 4: 'CCITT Group 4', 6: 'old JPEG', 7: 'JPEG', 34712: 'JPEG 2000', 34887: 'LERC', 34925: 'LZMA', 50000: 'Zstandard', 50001: 'WebP' };
async function readTIFF(ctx) {
  const u8 = await ctx.bytes(), dv = new DataView(u8.buffer, u8.byteOffset, u8.byteLength), le = u8[0] === 0x49 && u8[1] === 0x49, conv = fmtNamed('GeoTIFF / TIFF').convert;
  if (u8.length < 8 || (!le && !(u8[0] === 0x4d && u8[1] === 0x4d))) fail('Not a TIFF file.');
  const magic = dv.getUint16(2, le);
  if (magic === 43) fail('BigTIFF files are not supported. ' + conv);
  if (magic !== 42) fail('Not a TIFF file.');
  const pages = [], seen = new Set();
  for (let off = dv.getUint32(4, le); off && pages.length < 20000; ) {
    if (seen.has(off) || off + 2 > u8.length) break;
    seen.add(off);
    const n = dv.getUint16(off, le), tags = new Map();
    if (off + 2 + 12 * n + 4 > u8.length) fail('The TIFF directory is truncated.');
    for (let k = 0; k < n; k++) {
      const e = off + 2 + 12 * k, type = dv.getUint16(e + 2, le), count = dv.getUint32(e + 4, le), sz = TIFF_SZ[type] || 0, vo = sz * count <= 4 ? e + 8 : dv.getUint32(e + 8, le);
      if (sz && vo + sz * count <= u8.length) tags.set(dv.getUint16(e, le), { type, count, vo });
    }
    pages.push(tags);
    off = dv.getUint32(off + 2 + 12 * n, le);
  }
  const val = (tags, tag, dflt = []) => {
    const t = tags.get(tag);
    if (!t) return dflt;
    const out = [], n = Math.min(t.count, 4e6);
    for (let k = 0; k < n; k++) {
      const o = t.vo + k * TIFF_SZ[t.type];
      out.push(t.type === 3 ? dv.getUint16(o, le) : t.type === 4 ? dv.getUint32(o, le) : t.type === 12 ? dv.getFloat64(o, le) : t.type === 11 ? dv.getFloat32(o, le) : t.type === 8 ? dv.getInt16(o, le) : t.type === 9 ? dv.getInt32(o, le) : t.type === 6 ? dv.getInt8(o) : t.type === 5 ? dv.getUint32(o, le) / (dv.getUint32(o + 4, le) || 1) : t.type === 10 ? dv.getInt32(o, le) / (dv.getInt32(o + 4, le) || 1) : u8[o]);
    }
    return out;
  };
  const decode = async (tags) => {
    const w = val(tags, 256)[0], h = val(tags, 257)[0], bps = val(tags, 258, [1]), spp = val(tags, 277, [1])[0], comp = val(tags, 259, [1])[0], pred = val(tags, 317, [1])[0], sf = val(tags, 339, [1])[0], bits = bps[0];
    if (!(w > 0 && h > 0)) fail('The TIFF page has no image dimensions.');
    if (w * h > L.voxels) fail('The TIFF image is too large.');
    if (![1, 5, 8, 32946, 32773].includes(comp)) fail(`TIFF compression ${TIFF_COMP[comp] || comp} is not supported. ` + conv);
    if (![1, 8, 16, 32, 64].includes(bits) || bps.some((b) => b !== bits) || (bits === 64 && sf !== 3) || (sf === 3 && bits < 32) || (bits === 1 && spp !== 1)) fail(`TIFF images with ${bits}-bit ${sf === 3 ? 'float' : 'integer'} samples are not supported. ` + conv);
    if (spp > 1 && val(tags, 284, [1])[0] === 2) fail('Planar-separated multi-band TIFF images are not supported. ' + conv);
    if (![1, 2, 3].includes(pred) || (pred === 3 && sf !== 3)) fail(`TIFF predictor ${pred} is not supported. ` + conv);
    const tiled = tags.has(322), tw = tiled ? val(tags, 322)[0] : w, th = tiled ? val(tags, 323)[0] : Math.min(val(tags, 278, [h])[0] || h, h), offs = val(tags, tiled ? 324 : 273), cnts = val(tags, tiled ? 325 : 279);
    if (!(tw > 0 && th > 0)) fail('The TIFF strip / tile layout is invalid.');
    const across = Math.ceil(w / tw), nChunks = across * Math.ceil(h / th), dt = sf === 3 ? (bits === 32 ? 'f4' : 'f8') : sf === 2 ? { 8: 'i1', 16: 'i2', 32: 'i4' }[bits] : { 1: 'u1', 8: 'u1', 16: 'u2', 32: 'u4' }[bits];
    if (!dt || offs.length < nChunks) fail('The TIFF strip / tile table is incomplete.');
    const bp = bits >> 3, ns = tw * spp, rowBytes = bits === 1 ? Math.ceil(tw / 8) : ns * bp, out = spp === 1 ? new DT[dt][0](w * h) : new Float32Array(w * h), white0 = val(tags, 262, [1])[0] === 0;
    for (let c = 0; c < nChunks; c++) {
      const y0 = Math.floor(c / across) * th, x0 = (c % across) * tw, rows = tiled ? th : Math.min(th, h - y0), need = rowBytes * rows, cnt = cnts[c] ?? need;
      if (offs[c] + cnt > u8.length) fail('The TIFF file is truncated: a strip or tile lies beyond the end of the file.');
      const raw = u8.subarray(offs[c], offs[c] + cnt);
      let dec = comp === 1 ? raw : comp === 5 ? tiffLZW(raw, need) : comp === 32773 ? packBits(raw, need) : await inflate(raw, 'deflate', need + 65536);
      if (dec.length < need) { if (comp === 1) fail('The TIFF file is truncated: a strip is shorter than the image needs.'); const t = new Uint8Array(need); t.set(dec); dec = t; }
      let vals;
      if (bits === 1) { vals = new Uint8Array(tw * rows); for (let r = 0; r < rows; r++) for (let x = 0; x < tw; x++) { const b = (dec[r * rowBytes + (x >> 3)] >> (7 - (x & 7))) & 1; vals[r * tw + x] = white0 ? 1 - b : b; } }
      else if (pred === 3) {
        const t = new Uint8Array(need);
        for (let r = 0; r < rows; r++) { const o = r * rowBytes; for (let i = spp; i < rowBytes; i++) dec[o + i] = (dec[o + i] + dec[o + i - spp]) & 255; for (let s = 0; s < ns; s++) for (let k = 0; k < bp; k++) t[o + s * bp + (bp - 1 - k)] = dec[o + k * ns + s]; }
        vals = typed(t, 0, ns * rows, dt, true);
      } else {
        vals = typed(dec, 0, ns * rows, dt, le);
        if (pred === 2) for (let r = 0; r < rows; r++) for (let i = spp, o = r * ns; i < ns; i++) vals[o + i] += vals[o + i - spp];
      }
      const m = Math.min(3, spp);
      for (let r = 0; r < rows && y0 + r < h; r++) for (let x = 0; x < tw && x0 + x < w; x++) {
        if (spp === 1) out[(y0 + r) * w + x0 + x] = vals[r * tw + x];
        else { let s = 0; for (let q = 0; q < m; q++) s += vals[(r * tw + x) * spp + q]; out[(y0 + r) * w + x0 + x] = s / m; }
      }
    }
    return { w, h, data: out, real: sf === 3, bits };
  };
  const full = pages.filter((t) => t.has(256) && t.has(257) && !(val(t, 254, [0])[0] & 1));
  if (!full.length) fail('The TIFF file holds no full-resolution image.');
  const t0 = full[0], first = await decode(t0), { w, h } = first, warnings = [];
  const scale = val(t0, 33550), tie = val(t0, 33922), xfm = val(t0, 34264), keys = val(t0, 34735), nd = tags2str(t0, 42113, u8), geoKey = (id) => { for (let k = 4; k + 3 < keys.length; k += 4) if (keys[k] === id && keys[k + 1] === 0) return keys[k + 3]; return undefined; };
  const same = full.filter((t) => val(t, 256)[0] === w && val(t, 257)[0] === h), geo = scale.length >= 2 || tie.length >= 6 || xfm.length >= 16;
  const flip = (src, dst, at) => { for (let r = 0; r < h; r++) for (let x = 0; x < w; x++) dst[at + (h - 1 - r) * w + x] = src[r * w + x]; };
  if (same.length > 1 && !geo && ctx.opts.as !== 'grid') {
    if (w * h * same.length > L.voxels) fail(`The TIFF stack has too many voxels (limit ${L.voxels / 1e6} million).`);
    const vol = new first.data.constructor(w * h * same.length);
    flip(first.data, vol, 0);
    for (let k = 1; k < same.length; k++) flip((await decode(same[k])).data, vol, k * w * h);
    if (same.length < full.length) warnings.push(`${full.length - same.length} pages of a different size were ignored.`);
    return voxelsFromValues(vol, w, h, same.length, ctx.opts.spacing || [1, 1, 1], ctx.opts, { warnings, format: 'TIFF stack', pathway: 'voxel', stats: { pages: same.length, bits: first.bits } });
  }
  if ((geo || first.real || ctx.opts.as === 'grid') && ctx.opts.as !== 'voxels') {
    const half = geoKey(1025) === 2 ? 0 : 0.5, nod = nd === null ? undefined : parseFloat(nd);
    let xOf = (i) => i + 0.5, yOf = (j) => h - j - 0.5;
    if (scale.length >= 2 && tie.length >= 6) { xOf = (i) => tie[3] + (i - tie[0] + half) * scale[0]; yOf = (j) => tie[4] - (j - tie[1] + half) * scale[1]; }
    else if (xfm.length >= 16) { xOf = (i) => xfm[3] + (i + half) * xfm[0]; yOf = (j) => xfm[7] + (j + half) * xfm[5]; if (xfm[1] || xfm[4]) warnings.push('The GeoTIFF is rotated; the rotation terms of its transformation are ignored.'); }
    else warnings.push('The TIFF carries no georeferencing; pixel indices are used as coordinates.');
    return gridGeom(first.data, w, h, xOf, yOf, { nodata: Number.isFinite(nod) ? (v) => v === nod || Math.abs(v - nod) <= 1e-6 * Math.abs(nod) : (v) => v !== v, geographic: geoKey(1024) === 2, warnings, stats: { bits: first.bits, pages: full.length } });
  }
  const img = new first.data.constructor(w * h);
  flip(first.data, img, 0);
  return voxelsFromValues(img, w, h, 1, ctx.opts.spacing || [1, 1, 1], ctx.opts, { warnings, format: 'TIFF image', pathway: 'voxel', stats: { bits: first.bits } });
}
const tags2str = (tags, tag, u8) => { const t = tags.get(tag); return t && t.type === 2 ? latin1.decode(u8.subarray(t.vo, t.vo + t.count)).replace(/\0.*$/, '').trim() : null; };

const NC_TYPE = { 1: 'i1', 3: 'i2', 4: 'i4', 5: 'f4', 6: 'f8', 7: 'u1', 8: 'u2', 9: 'u4', 10: 'i8', 11: 'u8' };
/**
 * NetCDF classic header -> variable accessor { vars: [{ name, dims, shape, at, numeric, all(), slice2() }], at } shared by the
 * grid reader and Exodus II. all() gives every value in C order, slice2() the first 2-D slice over the last two dimensions.
 */
function cdfClassic(u8) {
  const dv = new DataView(u8.buffer, u8.byteOffset, u8.byteLength), ver = u8[3];
  if (u8.length < 32 || u8[0] !== 0x43 || u8[1] !== 0x44 || u8[2] !== 0x46 || ![1, 2, 5].includes(ver)) fail('Not a NetCDF classic file (CDF signature missing).');
  let p = 4;
  const bad = () => fail('The NetCDF header is truncated or corrupt.');
  const i4 = () => { if (p + 4 > u8.length) bad(); const v = dv.getInt32(p); p += 4; return v; };
  const nn = () => { if (ver !== 5) return i4() >>> 0; if (p + 8 > u8.length) bad(); const v = Number(dv.getBigUint64(p)); p += 8; return v; };
  const offs = () => { if (ver === 1) return i4() >>> 0; if (p + 8 > u8.length) bad(); const v = Number(dv.getBigUint64(p)); p += 8; return v; };
  const pad = (n) => (4 - (n % 4)) % 4;
  const name = () => { const n = nn(); if (n > 4096 || p + n > u8.length) bad(); const s = utf8.decode(u8.subarray(p, p + n)); p += n + pad(n); return s; };
  const list = (want, item) => { const tag = i4(), n = nn(); if (tag === 0 && n === 0) return []; if (tag !== want || n > 1e5) bad(); const out = []; for (let k = 0; k < n; k++) out.push(item()); return out; };
  const atts = () => {
    const a = Object.create(null);
    list(12, () => {
      const nm = name(), t = i4(), n = nn(), sz = t === 2 ? 1 : NC_TYPE[t] ? DT[NC_TYPE[t]][1] : 0;
      if (!sz || p + n * sz > u8.length) bad();
      a[nm] = t === 2 ? utf8.decode(u8.subarray(p, p + n)).replace(/\0+$/, '') : Array.from(typed(u8, p, Math.min(n, 1e5), NC_TYPE[t], false));
      p += n * sz + pad(n * sz);
    });
    return a;
  };
  const numrecs = nn(), dims = list(10, () => ({ name: name(), len: nn() })), gat = atts();
  const raw = list(11, () => { const nm = name(), nd = nn(); if (nd > 32) bad(); const dimids = []; for (let k = 0; k < nd; k++) { const d = nn(); if (d >= dims.length) bad(); dimids.push(d); } const at = atts(), type = i4(), vsize = nn(), begin = offs(); return { name: nm, dimids, at, type, vsize, begin }; });
  const isRec = (v) => v.dimids.length > 0 && dims[v.dimids[0]].len === 0, shape = (v) => v.dimids.map((d, k) => (k === 0 && dims[d].len === 0 ? numrecs : dims[d].len));
  const recVars = raw.filter(isRec);
  const recSize = recVars.length === 1 && NC_TYPE[recVars[0].type] ? shape(recVars[0]).slice(1).reduce((a, b) => a * b, 1) * DT[NC_TYPE[recVars[0].type]][1] : recVars.reduce((a, v) => a + v.vsize, 0);
  /** First n values of a variable in C order (record variables are gathered record by record). */
  const first = (v, n) => {
    const dt = NC_TYPE[v.type];
    if (!dt) fail(`NetCDF variable "${v.name}" is not numeric.`);
    if (!(n <= u8.length)) bad();
    if (!isRec(v)) return typed(u8, v.begin, n, dt, false);
    const per = shape(v).slice(1).reduce((a, b) => a * b, 1), out = new Float64Array(n);
    for (let r = 0, o = 0; o < n; r++, o += per) out.set(typed(u8, v.begin + r * recSize, Math.min(per, n - o), dt, false), o);
    return out;
  };
  const vars = raw.map((v) => { const sh = shape(v); return { name: v.name, dims: v.dimids.map((d) => dims[d].name), shape: sh, at: v.at, numeric: !!NC_TYPE[v.type], all: async () => first(v, sh.reduce((a, b) => a * b, 1)), slice2: async () => first(v, sh.slice(-2).reduce((a, b) => a * b, 1)) }; });
  return { vars, at: gat, flavour: `CDF-${ver}` };
}
/** The same accessor over a NetCDF-4 file (HDF5 datasets; dimensions resolved through DIMENSION_LIST references). */
function cdfHDF5(h5) {
  const byAddr = new Map(), vars = [], norm = (at) => { const o = Object.create(null); for (const k in at) { const v = at[k]; o[k] = typeof v === 'number' ? [v] : ArrayBuffer.isView(v) ? Array.from(v) : v; } return o; };
  h5.walk((path, addr, inf) => { if (inf.kind === 'dataset') byAddr.set(addr, path.slice(1)); }, { maxDepth: 8 });
  for (const [addr, name] of byAddr) {
    const inf = h5.info(addr), at = h5.attrs(addr), dl = at.DIMENSION_LIST, sh = inf.shape, rank = sh.length;
    if (typeof at.NAME === 'string' && at.NAME.startsWith('This is a netCDF dimension but not a netCDF variable')) continue;
    const dimNames = Array.isArray(dl) && dl.length === rank ? dl.map((r, k) => (r && r.length ? byAddr.get(r[0]) : undefined) ?? `dim${k}`) : rank === 1 ? [name] : sh.map((_, k) => `dim${k}`);
    vars.push({ name, dims: dimNames, shape: sh, at: norm(at), numeric: isNumeric(inf.type) && inf.type.cls !== 10, all: async () => (await h5.read(addr)).data, slice2: async () => (await h5.read(addr, { count: sh.map((n, k) => (k < rank - 2 ? 1 : n)) })).data });
  }
  return { vars, at: norm(h5.attrs(h5.root)), flavour: 'NetCDF-4' };
}
/** First gridded variable of a NetCDF data set (classic or NetCDF-4) with its coordinate variables -> grid. */
async function cdfGrid(nc) {
  const { vars } = nc, warnings = [], byName = (nm) => vars.find((v) => v.name === nm);
  // legacy GMT layout: z(xysize) with x_range / y_range / dimension
  const gz = byName('z'), gd = byName('dimension'), gx = byName('x_range'), gy = byName('y_range');
  if (gz && gd && gx && gy && gz.shape.length === 1 && gz.numeric) {
    const [nx, ny] = await gd.all(), xr = await gx.all(), yr = await gy.all();
    if (!(nx > 0 && ny > 0) || nx * ny > L.voxels || nx * ny > gz.shape[0]) fail('The NetCDF header is truncated or corrupt.');
    const z = (await gz.all()).subarray(0, nx * ny), sf = (gz.at.scale_factor || [1])[0], ao = (gz.at.add_offset || [0])[0], zz = sf === 1 && ao === 0 ? z : Float64Array.from(z, (q) => q * sf + ao);
    return gridGeom(zz, nx, ny, (i) => xr[0] + ((xr[1] - xr[0]) * (i + 0.5)) / nx, (j) => yr[1] - ((yr[1] - yr[0]) * (j + 0.5)) / ny, { nodata: (q) => q !== q, warnings, stats: { variable: 'z', layout: 'GMT v3' } });
  }
  const cand = vars.filter((v) => v.numeric && v.shape.length >= 2 && v.shape.slice(-2).every((n) => n > 1) && !/bnds|bounds|^crs$/i.test(v.name));
  const v = cand.find((q) => /^(z|elev|elevation|topo|bathy|bathymetry|depth|height|band1|dem|altitude)/i.test(q.name)) || cand[0];
  if (!v) fail('No 2-D numeric variable was found in the NetCDF file.');
  const sh = v.shape, ny = sh[sh.length - 2], nx = sh[sh.length - 1];
  if (nx * ny > L.voxels) fail('The NetCDF variable is too large.');
  const raw = await v.slice2();
  if (sh.length > 2) warnings.push(`Variable "${v.name}" has ${sh.length} dimensions (${sh.join(' × ')}); its first 2-D slice is used.`);
  const sf = (v.at.scale_factor || [1])[0], ao = (v.at.add_offset || [0])[0], fills = [v.at._FillValue, v.at.missing_value].filter(Array.isArray).map((a) => a[0]);
  const z = new Float64Array(nx * ny);
  for (let k = 0; k < z.length; k++) { const q = raw[k]; z[k] = fills.includes(q) || q !== q ? NaN : q * sf + ao; }
  const coord = async (dim, n) => { const cv = vars.find((q) => q.name === dim && q.shape.length === 1 && q.numeric); if (!cv) return { f: (i) => i, v: null }; const a = await cv.all(); return a.length >= n ? { f: (i) => a[i], v: cv } : { f: (i) => i, v: null }; };
  const cx = await coord(v.dims[v.dims.length - 1], nx), cy = await coord(v.dims[v.dims.length - 2], ny);
  if (!cx.v || !cy.v) warnings.push('No coordinate variables were found for the grid; cell indices are used as coordinates.');
  const geographic = !!cx.v && (/degree/i.test(String(cx.v.at.units || '')) || /^lon/i.test(cx.v.name));
  return gridGeom(z, nx, ny, cx.f, cy.f, { nodata: (q) => q !== q, geographic, warnings, stats: { variable: v.name, units: typeof v.at.units === 'string' ? v.at.units : undefined, dimensions: v.dims } });
}
async function readNetCDF(ctx) {
  const u8 = await ctx.bytes();
  if (hdf5Offset(u8) === 0) return readHDF(ctx, 'netcdf');
  return Object.assign(await cdfGrid(cdfClassic(u8)), { format: 'NetCDF classic' });
}
async function readExodus(ctx) {
  const u8 = await ctx.bytes();
  if (hdf5Offset(u8) === 0) return readHDF(ctx, 'exodus');
  if (!(u8[0] === 0x43 && u8[1] === 0x44 && u8[2] === 0x46)) fail('Not an Exodus II file: it is neither a NetCDF classic nor a NetCDF-4 / HDF5 container.');
  return Object.assign(await exodusGeom(cdfClassic(u8)), { format: 'Exodus II' });
}

async function readBandRaster(ctx) {
  const hf = ctx.companion((n) => /\.hdr$/i.test(n));
  if (!hf) fail(fmtNamed('ENVI / ESRI band raster').convert);
  const ht = utf8.decode(new Uint8Array(await hf.arrayBuffer())), kv = Object.create(null), u8 = await ctx.bytes();
  let nx, ny, nb, off, dt, il, le, x0, y0, dx, dy, nod;
  if (/^\s*ENVI/i.test(ht)) {
    for (const m of ht.matchAll(/^[ \t]*([^=\n\r]+?)[ \t]*=[ \t]*(\{[^}]*\}|[^\n\r]*)/gm)) kv[m[1].toLowerCase()] = m[2].trim();
    nx = +kv.samples; ny = +kv.lines; nb = +(kv.bands || 1); off = +(kv['header offset'] || 0); il = (kv.interleave || 'bsq').toLowerCase(); le = +(kv['byte order'] || 0) === 0;
    dt = { 1: 'u1', 2: 'i2', 3: 'i4', 4: 'f4', 5: 'f8', 12: 'u2', 13: 'u4', 14: 'i8', 15: 'u8' }[+kv['data type']];
    const mi = (kv['map info'] || '').replace(/[{}]/g, '').split(',').map((s) => s.trim());
    if (mi.length >= 7) { dx = +mi[5]; dy = +mi[6]; x0 = +mi[3] - (+mi[1] - 1) * dx; y0 = +mi[4] + (+mi[2] - 1) * dy; }
    nod = kv['data ignore value'] !== undefined ? +kv['data ignore value'] : undefined;
  } else {
    for (const l of splitLines(ht)) { const t = l.trim().split(/\s+/); if (t.length >= 2) kv[t[0].toLowerCase()] = t[1]; }
    nx = +kv.ncols; ny = +kv.nrows; nb = +(kv.nbands || 1); off = +(kv.skipbytes || 0); il = (kv.layout || 'bil').toLowerCase(); le = (kv.byteorder || 'I').toUpperCase() !== 'M';
    const bits = +(kv.nbits || 8), pt = (kv.pixeltype || 'UNSIGNEDINT').toUpperCase();
    dt = pt.startsWith('FLOAT') ? { 32: 'f4', 64: 'f8' }[bits] : pt.startsWith('SIGNED') ? { 8: 'i1', 16: 'i2', 32: 'i4' }[bits] : { 8: 'u1', 16: 'u2', 32: 'u4' }[bits];
    dx = +(kv.xdim || 1); dy = +(kv.ydim || 1);
    if (kv.ulxmap !== undefined) { x0 = +kv.ulxmap - dx / 2; y0 = +kv.ulymap + dy / 2; }
    nod = kv.nodata !== undefined ? +kv.nodata : undefined;
  }
  if (![nx, ny, nb].every((k) => Number.isInteger(k) && k > 0) || !dt || !['bsq', 'bil', 'bip'].includes(il)) fail('The raster header (.hdr) is incomplete or uses an unsupported data type.');
  if (nx * ny > L.voxels) fail('The raster is too large.');
  let z;
  if (il === 'bsq' || nb === 1) z = typed(u8, off, nx * ny, dt, le);
  else {
    if (nx * ny * nb > 4 * L.voxels) fail('The multi-band raster is too large.');
    const all = typed(u8, off, nx * ny * nb, dt, le);
    z = new Float64Array(nx * ny);
    for (let r = 0; r < ny; r++) for (let c = 0; c < nx; c++) z[r * nx + c] = il === 'bil' ? all[r * nx * nb + c] : all[(r * nx + c) * nb];
  }
  if (!Number.isFinite(x0) || !Number.isFinite(y0)) { x0 = 0; dx = dx > 0 ? dx : 1; dy = dy > 0 ? dy : 1; y0 = ny * dy; }
  return gridGeom(z, nx, ny, (i) => x0 + (i + 0.5) * dx, (j) => y0 - (j + 0.5) * dy, { nodata: nod, warnings: nb > 1 ? [`Band 1 of ${nb} was read.`] : [], stats: { bands: nb, interleave: il } });
}

// ---- Column text, point clouds ---------------------------------------------------------------------------
/** Numeric column text -> { headers, ncol, rows, data (row-major, NaN for non-numeric cells), breaks, textCells }. */
function readColumns(text, maxRows = 8e6) {
  if (text.indexOf('\n') < 0 && text.indexOf('\r') >= 0) text = text.replace(/\r/g, '\n');
  const n = text.length, breaks = [], isSep = (c) => c === 32 || c === 9 || c === 44 || c === 59 || c === 13;
  let headers = null, ncol = 0, rows = 0, data = new Float64Array(4096), textCells = 0, started = false, comment = null, pos = text.charCodeAt(0) === 0xfeff ? 1 : 0;
  while (pos < n) {
    let e = text.indexOf('\n', pos), a = pos;
    if (e < 0) e = n;
    pos = e + 1;
    while (a < e && isSep(text.charCodeAt(a))) a++;
    if (a >= e) { if (rows && breaks[breaks.length - 1] !== rows) breaks.push(rows); continue; }
    const c0 = text.charCodeAt(a);
    if (c0 === 35 || c0 === 37 || c0 === 33 || (c0 === 47 && text.charCodeAt(a + 1) === 47)) { if (!started && !comment) comment = text.slice(a, e).replace(/^(#|\/\/|%|!)+\s*/, '').trim(); continue; }
    if (!started || !ncol) {
      const parts = text.slice(a, e).trim().split(/[,;\t ]+/).filter((s) => s !== '');
      if (!started) { started = true; if (parts.some((s) => !Number.isFinite(+s.replace(/^"|"$/g, '')))) { headers = parts.map((s) => s.replace(/^"|"$/g, '')); continue; } }
      ncol = parts.length;
    }
    if (rows >= maxRows) fail('The text file holds too many rows.');
    if ((rows + 1) * ncol > data.length) { const g = new Float64Array(Math.max(data.length * 2, (rows + 1) * ncol)); g.set(data); data = g; }
    let k = 0;
    for (let p = a; p < e && k < ncol; k++) {
      let q = p;
      while (q < e && !isSep(text.charCodeAt(q))) q++;
      const v = +text.slice(p, q);
      if (v !== v) textCells++;
      data[rows * ncol + k] = v;
      p = q;
      while (p < e && isSep(text.charCodeAt(p))) p++;
    }
    for (; k < ncol; k++) { data[rows * ncol + k] = NaN; textCells++; }
    rows++;
  }
  if (!headers && comment) { const h = comment.split(/[,;\t ]+/).filter(Boolean); if (h.length === ncol && h.every((s) => !Number.isFinite(+s))) headers = h; }
  if (headers && ncol && headers.length !== ncol) headers = headers.length > ncol ? headers.slice(0, ncol) : null;
  return { headers, ncol, rows, data: data.subarray(0, rows * ncol), breaks: breaks.filter((b) => b < rows), textCells };
}
const normHeader = (h) => String(h).toLowerCase().replace(/[[(].*$/, '').replace(/[_\s]+(m|mm|km|ft|usft|deg|degrees|dd|s|h|hr)$/, '').replace(/[^a-z0-9]/g, '');
const ROLE = { x: /^(x|lon|long|longitude|lng|easting|east|e|xcoord|xm)$/, y: /^(y|lat|latitude|northing|north|n|ycoord|ym)$/, z: /^(z|elev|elevation|height|alt|altitude|h|level|zcoord|zm)$/, depth: /^(depth|d|sounding|soundings|bathy|bathymetry|waterdepth)$/ };
function columnRoles(headers) {
  if (!headers) return null;
  const hn = headers.map(normHeader), find = (re) => hn.findIndex((h) => re.test(h)), r = { x: find(ROLE.x), y: find(ROLE.y), z: find(ROLE.z), depth: find(ROLE.depth) };
  if (r.x < 0 || r.y < 0) return null;
  r.geo = /^(lon|long|longitude|lng)$/.test(hn[r.x]) && /^(lat|latitude)$/.test(hn[r.y]);
  return r;
}
/** Column table -> xyz triples, honouring header names; depth columns (positive down) become negative elevations. */
function tableToXYZ(tab, warnings) {
  const r = columnRoles(tab.headers) || { x: 0, y: 1, z: tab.ncol > 2 ? 2 : -1, depth: -1, geo: false }, zc = r.z >= 0 ? r.z : r.depth, sign = r.z < 0 && r.depth >= 0 ? -1 : 1;
  if (tab.ncol < 2) fail('The file needs at least two numeric columns (x, y).');
  if (sign < 0) warnings.push(`Column "${tab.headers[r.depth]}" is treated as depth (positive down) and converted to negative elevation.`);
  const xyz = new Float64Array(tab.rows * 3), d = tab.data, nc = tab.ncol;
  for (let i = 0; i < tab.rows; i++) { xyz[3 * i] = d[i * nc + r.x]; xyz[3 * i + 1] = d[i * nc + r.y]; xyz[3 * i + 2] = zc >= 0 ? sign * d[i * nc + zc] : 0; }
  return { xyz, roles: r, hasZ: zc >= 0 };
}
/** Scattered xyz that actually form a complete regular lattice -> grid Geometry, else null. */
function latticeGrid(xyz, extra) {
  const n = xyz.length / 3;
  if (n < 4 || n > L.grid) return null;
  const uniq = (o) => { const s = new Set(); for (let i = 0; i < n; i++) { s.add(xyz[3 * i + o]); if (s.size > 4096) return null; } return [...s].sort((a, b) => a - b); };
  const xs = uniq(0), ys = xs && uniq(1);
  if (!xs || !ys || xs.length < 2 || ys.length < 2 || xs.length * ys.length !== n) return null;
  const xi = new Map(xs.map((v, k) => [v, k])), yi = new Map(ys.map((v, k) => [v, k])), z = new Float64Array(n).fill(NaN);
  for (let i = 0; i < n; i++) z[yi.get(xyz[3 * i + 1]) * xs.length + xi.get(xyz[3 * i])] = xyz[3 * i + 2];
  return gridGeom(z, xs.length, ys.length, (i) => xs[i], (j) => ys[j], { ...extra, nodata: (v) => v !== v });
}
async function readPoints(ctx, text) {
  let t = text ?? (await ctx.text());
  if (ctx.ext === 'pts') t = t.replace(/^\s*\d+\s*(\r\n|\n|\r)/, '');
  const tab = readColumns(t), warnings = [];
  if (tab.rows < 1 || tab.ncol < 3) fail('The point file needs rows of at least three numbers (x y z).');
  if (tab.ncol > 3) warnings.push(`Columns beyond x, y, z (intensity / colour, ${tab.ncol - 3} per row) are ignored.`);
  const { xyz, roles } = tableToXYZ(tab, warnings);
  const g = pointsGeom(xyz, { warnings, stats: { columns: tab.headers || undefined } });
  if (roles.geo) g.geographic = true;
  return g;
}
async function readPTX(ctx) {
  const ls = splitLines(await ctx.text()), out = [];
  let i = 0, scans = 0, empty = 0;
  while (i + 10 <= ls.length) {
    if (!ls[i].trim()) { i++; continue; }
    const nc = parseInt(ls[i], 10), nr = parseInt(ls[i + 1], 10), M = [6, 7, 8, 9].map((k) => Array.from(numsOf(ls[i + k])));
    if (!(nc > 0 && nr > 0) || M.some((r) => r.length < 3 || r.some((v) => !Number.isFinite(v)))) { if (!scans) fail('Not a PTX scan file (header of column / row counts and matrices missing).'); break; }
    i += 10;
    const n = Math.min(nc * nr, ls.length - i);
    if (out.length / 3 + n > 12e6) fail('The PTX file holds too many points.');
    for (let k = 0; k < n; k++, i++) {
      const p = ls[i].trim().split(/\s+/), x = +p[0], y = +p[1], z = +p[2];
      if (x === 0 && y === 0 && z === 0) { empty++; continue; }
      out.push(x * M[0][0] + y * M[1][0] + z * M[2][0] + M[3][0], x * M[0][1] + y * M[1][1] + z * M[2][1] + M[3][1], x * M[0][2] + y * M[1][2] + z * M[2][2] + M[3][2]);
    }
    scans++;
  }
  if (!scans) fail('Not a PTX scan file.');
  return pointsGeom(out, { warnings: empty ? [`${empty} empty returns (0 0 0) were dropped.`] : [], stats: { scans } });
}
async function readLAS(ctx) {
  const u8 = await ctx.bytes(), dv = new DataView(u8.buffer, u8.byteOffset, u8.byteLength);
  if (u8.length < 227 || latin1.decode(u8.subarray(0, 4)) !== 'LASF') fail('Not a LAS file (LASF signature missing).');
  const major = u8[24], minor = u8[25], hsize = dv.getUint16(94, true), off = dv.getUint32(96, true), fmtRaw = u8[104], fmt = fmtRaw & 0x3f, rec = dv.getUint16(105, true);
  if (fmtRaw & 0xc0) return Object.assign(await readLAZ(ctx), { format: 'LAZ (compressed LAS)' });
  if (major !== 1 || minor > 4 || fmt > 10 || rec < 20 || off < hsize || off > u8.length) fail(`Unsupported or corrupt LAS file (version ${major}.${minor}, point format ${fmt}).`);
  let n = dv.getUint32(107, true);
  if (minor >= 4 && hsize >= 375 && u8.length >= 255) { const n64 = Number(dv.getBigUint64(247, true)); if (n64 > 0) n = n64; }
  const avail = Math.floor((u8.length - off) / rec), warnings = [], f = (o) => dv.getFloat64(o, true), sc = [f(131), f(139), f(147)], of = [f(155), f(163), f(171)];
  if (n > avail) { warnings.push(`The header declares ${n} points but the file holds ${avail}; the file appears truncated.`); n = avail; }
  if (!(n > 0)) fail('The LAS file holds no points.');
  const step = Math.max(1, Math.ceil(n / L.points)), m = Math.ceil(n / step), xyz = new Float64Array(m * 3);
  for (let k = 0, i = 0; i < n; i += step, k++) { const o = off + i * rec; xyz[3 * k] = dv.getInt32(o, true) * sc[0] + of[0]; xyz[3 * k + 1] = dv.getInt32(o + 4, true) * sc[1] + of[1]; xyz[3 * k + 2] = dv.getInt32(o + 8, true) * sc[2] + of[2]; }
  if (step > 1) warnings.push(`Point cloud of ${n} points evenly sub-sampled to ${m} (every ${step}th point).`);
  return pointsGeom(xyz, { warnings, stats: { version: `${major}.${minor}`, pointFormat: fmt, sourcePoints: n } });
}

async function readLAZ(ctx) {
  const u8 = await ctx.bytes();
  if (u8.length > 110 && !isLAZ(u8)) return readLAS(ctx);
  const d = parseLAZ(u8, { maxPoints: L.points });
  return pointsGeom(d.xyz, { warnings: d.warnings, stats: { version: d.version, pointFormat: d.pointFormat, sourcePoints: d.total, compression: `LASzip ${d.lazVersion}, ${d.compressor}`, chunks: d.chunks } });
}
async function readE57(ctx) {
  const d = parseE57(await ctx.bytes(), { maxPoints: L.points });
  return pointsGeom(d.xyz, { warnings: d.warnings, stats: { version: d.version, scans: d.scans.length, sourcePoints: d.total, invalidPoints: d.invalid, scanNames: d.scans.map((s) => s.name).filter(Boolean).slice(0, 20) } });
}

// ---- Images and voxel volumes -----------------------------------------------------------------------------
const flipRows = (src, w, h, nz = 1) => { const out = new src.constructor(src.length); for (let k = 0; k < nz; k++) for (let r = 0; r < h; r++) out.set(src.subarray(k * w * h + r * w, k * w * h + (r + 1) * w), k * w * h + (h - 1 - r) * w); return out; };
async function readImage(ctx) {
  if (typeof createImageBitmap !== 'function' || (typeof OffscreenCanvas !== 'function' && typeof document === 'undefined')) fail('Image decoding (PNG / JPEG / BMP / WebP) needs a browser. Outside a browser convert the image to TIFF, NRRD or NumPy .npy first.');
  let bmp;
  try { bmp = await createImageBitmap(new Blob([await ctx.buf()])); } catch { fail('The image could not be decoded by this browser.'); }
  const w = bmp.width, h = bmp.height;
  if (!(w > 0 && h > 0) || w * h > L.voxels) fail('The image is empty or too large.');
  const cv = typeof OffscreenCanvas === 'function' ? new OffscreenCanvas(w, h) : Object.assign(document.createElement('canvas'), { width: w, height: h }), c2 = cv.getContext('2d', { willReadFrequently: true });
  c2.drawImage(bmp, 0, 0);
  const px = c2.getImageData(0, 0, w, h).data, grey = new Uint8Array(w * h);
  if (bmp.close) bmp.close();
  for (let i = 0; i < w * h; i++) { const a = px[4 * i + 3] / 255; grey[i] = Math.round((0.299 * px[4 * i] + 0.587 * px[4 * i + 1] + 0.114 * px[4 * i + 2]) * a + 255 * (1 - a)); }
  return voxelsFromValues(flipRows(grey, w, h), w, h, 1, ctx.opts.spacing || [1, 1, 1], ctx.opts, { warnings: ['Image rows are flipped so that y points up; spacing is one unit per pixel unless opts.spacing is given.'] });
}

async function readRAW(ctx) {
  const u8 = await ctx.bytes(), o = ctx.opts, dt = dtypeOf(o.dtype || 'uint8'), off = Math.max(0, o.headerBytes | 0), warnings = [];
  if (!dt) fail(`Unknown data type "${o.dtype}" (use uint8, int16, uint16, int32, float32 or float64).`);
  const n = (u8.length - off) / DT[dt][1];
  let dims = Array.isArray(o.dims) ? o.dims.map(Number) : null;
  if (!dims) {
    const c = Math.round(Math.cbrt(n)), s = Math.round(Math.sqrt(n));
    if (dt === 'u1' && c > 1 && c * c * c === n) dims = [c, c, c];
    else if (dt === 'u1' && s > 1 && s * s === n) dims = [s, s, 1];
    else fail(`RAW volumes have no header: give the dimensions (dims = [nx, ny, nz]) and data type (dtype). This file holds ${u8.length - off} bytes, which is not a perfect cube or square of 8-bit voxels.`);
    warnings.push(`No dimensions given: inferred ${dims.join(' × ')} 8-bit voxels from the file size.`);
  }
  if (dims.length === 2) dims.push(1);
  const [nx, ny, nz] = dims;
  if (dims.length !== 3 || !dims.every((k) => Number.isInteger(k) && k > 0)) fail('RAW dimensions must be three positive integers [nx, ny, nz].');
  if (nx * ny * nz > L.voxels) fail(`Volume has too many voxels (limit ${L.voxels / 1e6} million).`);
  if (nx * ny * nz > n) fail(`The file is too small for ${nx} × ${ny} × ${nz} voxels of type ${o.dtype || 'uint8'} (needs ${nx * ny * nz * DT[dt][1] + off} bytes, has ${u8.length}).`);
  if (nx * ny * nz < Math.floor(n)) warnings.push('The file is larger than the given dimensions need; trailing bytes are ignored.');
  return voxelsFromValues(typed(u8, off, nx * ny * nz, dt, o.littleEndian !== false), nx, ny, nz, o.spacing || [1, 1, 1], o, { warnings });
}

async function readNRRD(ctx) {
  const u8 = await ctx.bytes(), head = latin1.decode(u8.subarray(0, Math.min(u8.length, 1 << 16)));
  if (!/^NRRD\d{4}/.test(head)) fail('Not an NRRD file (NRRD magic missing).');
  const m = /\r?\n\r?\n/.exec(head), kv = Object.create(null);
  for (const l of splitLines(m ? head.slice(0, m.index) : head)) { const q = /^([^:#=]+):\s+(.*)$/.exec(l); if (q) kv[q[1].trim().toLowerCase()] = q[2].trim(); }
  const dt = dtypeOf(kv.type), sizes = Array.from(numsOf(kv.sizes || '')), enc = (kv.encoding || 'raw').toLowerCase(), file = kv['data file'] || kv.datafile;
  if (!dt || sizes.length < 2 || sizes.length > 3 || +kv.dimension !== sizes.length) fail(sizes.length > 3 ? 'NRRD files with more than three axes (vector or time data) are not supported.' : 'The NRRD header is incomplete (type, dimension, sizes).');
  const [nx, ny, nz = 1] = sizes, n = nx * ny * nz;
  if (![nx, ny, nz].every((k) => Number.isInteger(k) && k > 0) || n > L.voxels) fail('The NRRD sizes are invalid or the volume is too large.');
  let data;
  if (file) { const c = ctx.companion((nm) => baseName(nm) === baseName(file)); if (!c) fail(`This NRRD header refers to a detached data file; supply "${file}" as a companion.`); data = new Uint8Array(await c.arrayBuffer()); }
  else { if (!m) fail('The NRRD file has no data section.'); data = u8.subarray(m.index + m[0].length); }
  let vals;
  if (enc === 'gzip' || enc === 'gz') data = await inflate(data, 'gzip', n * DT[dt][1] + (1 << 20));
  else if (!['raw', 'ascii', 'text', 'txt'].includes(enc)) fail(`NRRD encoding "${enc}" is not supported; re-save the volume with raw or gzip encoding.`);
  if (enc === 'ascii' || enc === 'text' || enc === 'txt') vals = numsOf(latin1.decode(data));
  else { const skip = +kv['byte skip'] || 0; vals = typed(data, skip === -1 ? data.length - n * DT[dt][1] : Math.max(0, skip), n, dt, (kv.endian || 'little').toLowerCase() !== 'big'); }
  let sp = Array.from(numsOf(kv.spacings || ''));
  if (!sp.length && kv['space directions']) sp = (kv['space directions'].match(/\([^)]*\)/g) || []).map((v) => Math.hypot(...numsOf(v.replace(/[()]/g, ''))));
  return voxelsFromValues(vals, nx, ny, nz, sp, ctx.opts, { stats: { encoding: enc, type: kv.type } });
}

async function readMHA(ctx) {
  const u8 = await ctx.bytes(), head = latin1.decode(u8.subarray(0, Math.min(u8.length, 1 << 16))), kv = Object.create(null);
  let p = 0;
  while (p < head.length) {
    const e = head.indexOf('\n', p), line = head.slice(p, e < 0 ? undefined : e), q = /^\s*(\w+)\s*=\s*(.*?)\s*$/.exec(line);
    p = e < 0 ? head.length : e + 1;
    if (!q) break;
    kv[q[1].toLowerCase()] = q[2];
    if (q[1].toLowerCase() === 'elementdatafile') break;
  }
  if (!kv.elementdatafile || !kv.dimsize) fail('Not a MetaImage file (DimSize / ElementDataFile missing).');
  const dims = Array.from(numsOf(kv.dimsize)), dt = dtypeOf(kv.elementtype), yes = (s) => /^true$/i.test(s || ''), [nx, ny, nz = 1] = dims, n = nx * ny * nz;
  if (!dt || dims.length < 2 || dims.length > 3 || ![nx, ny, nz].every((k) => Number.isInteger(k) && k > 0) || n > L.voxels) fail('The MetaImage header is invalid, has more than three dimensions, or the volume is too large.');
  if (+(kv.elementnumberofchannels || 1) !== 1) fail('Multi-channel MetaImage volumes are not supported.');
  let data;
  if (/^local$/i.test(kv.elementdatafile)) data = u8.subarray(p);
  else if (/^list/i.test(kv.elementdatafile) || /%/.test(kv.elementdatafile)) fail('MetaImage file lists (one file per slice) are not supported; save the volume as a single .mha.');
  else { const c = ctx.companion((nm) => baseName(nm) === baseName(kv.elementdatafile)); if (!c) fail(`This MetaImage header refers to a separate data file; supply "${kv.elementdatafile}" as a companion.`); data = new Uint8Array(await c.arrayBuffer()); }
  if (yes(kv.compresseddata)) data = await inflate(data, 'deflate', n * DT[dt][1] + (1 << 20));
  const hs = +kv.headersize || 0, vals = /^false$/i.test(kv.binarydata || '') ? numsOf(latin1.decode(data)) : typed(data, hs === -1 ? data.length - n * DT[dt][1] : Math.max(0, hs), n, dt, !(yes(kv.binarydatabyteordermsb) || yes(kv.elementbyteordermsb)));
  return voxelsFromValues(vals, nx, ny, nz, Array.from(numsOf(kv.elementspacing || kv.elementsize || '')), ctx.opts, { stats: { type: kv.elementtype } });
}

const NII_TYPE = { 2: 'u1', 4: 'i2', 8: 'i4', 16: 'f4', 64: 'f8', 256: 'i1', 512: 'u2', 768: 'u4', 1024: 'i8', 1280: 'u8' };
async function readNIfTI(ctx) {
  let u8 = await ctx.bytes();
  if (isGzip(u8)) u8 = await inflate(u8, 'gzip');
  if (u8.length < 348) fail('Not a NIfTI-1 file (header too short).');
  const dv = new DataView(u8.buffer, u8.byteOffset, u8.byteLength), a = dv.getInt32(0, true), b = dv.getInt32(0, false);
  if (a === 540 || b === 540) fail('NIfTI-2 files are not supported; convert to NIfTI-1, NRRD or MetaImage.');
  if (a !== 348 && b !== 348) fail('Not a NIfTI-1 file (header size field is not 348).');
  const le = a === 348, magic = latin1.decode(u8.subarray(344, 347)), dim = [];
  for (let k = 0; k < 8; k++) dim.push(dv.getInt16(40 + 2 * k, le));
  const nx = dim[1], ny = Math.max(1, dim[2]), nz = dim[0] >= 3 ? Math.max(1, dim[3]) : 1, dt = NII_TYPE[dv.getInt16(70, le)], warnings = [];
  if (!dt) fail('The NIfTI data type (complex, RGB or 128-bit) is not supported.');
  if (!(nx > 0) || nx * ny * nz > L.voxels) fail('The NIfTI dimensions are invalid or the volume is too large.');
  let src = u8, off = Math.round(dv.getFloat32(108, le));
  if (magic === 'ni1') { const c = ctx.companion((nm) => /\.img(\.gz)?$/i.test(nm)); if (!c) fail('This NIfTI header (.hdr) needs its .img data file as a companion.'); src = new Uint8Array(await c.arrayBuffer()); if (isGzip(src)) src = await inflate(src, 'gzip'); off = 0; }
  else if (magic !== 'n+1') fail('Not a NIfTI-1 file (magic string missing).');
  if (dim[0] > 3 && dim[4] > 1) warnings.push(`4-D NIfTI with ${dim[4]} volumes: the first volume is used.`);
  let vals = typed(src, Math.max(off, magic === 'n+1' ? 348 : 0), nx * ny * nz, dt, le);
  const slope = dv.getFloat32(112, le), inter = dv.getFloat32(116, le);
  if (Number.isFinite(slope) && slope !== 0 && (slope !== 1 || inter !== 0)) vals = Float32Array.from(vals, (v) => v * slope + inter);
  return voxelsFromValues(vals, nx, ny, nz, [1, 2, 3].map((k) => Math.abs(dv.getFloat32(76 + 4 * k, le))), ctx.opts, { warnings });
}

function parseNPY(u8) {
  if (u8.length < 12 || u8[0] !== 0x93 || latin1.decode(u8.subarray(1, 6)) !== 'NUMPY') fail('Not a NumPy .npy file.');
  const major = u8[6], hlen = major === 1 ? u8[8] | (u8[9] << 8) : (u8[8] | (u8[9] << 8) | (u8[10] << 16) | (u8[11] << 24)) >>> 0, hoff = major === 1 ? 10 : 12;
  if (major < 1 || major > 3 || hoff + hlen > u8.length) fail('The .npy header is invalid or truncated.');
  const h = (major === 3 ? utf8 : latin1).decode(u8.subarray(hoff, hoff + hlen)), dm = /'descr'\s*:\s*'([^']*)'/.exec(h), fm = /'fortran_order'\s*:\s*(True|False)/.exec(h), sm = /'shape'\s*:\s*\(([^)]*)\)/.exec(h);
  if (!dm || !sm) fail(/'descr'\s*:\s*\[/.test(h) ? 'Structured NumPy arrays are not supported; save a plain numeric array.' : 'The .npy header is invalid.');
  const d = /^([<>|=])([bBuif])(\d)$/.exec(dm[1]), dt = d ? dtypeOf((d[2] === 'b' ? 'u' : d[2] === 'B' ? 'u' : d[2]) + d[3]) : null;
  if (!dt || (d[2] === 'f' && d[3] === '2')) fail(`NumPy dtype "${dm[1]}" is not supported (use bool, integer, float32 or float64).`);
  const shape = sm[1].split(',').map((s) => s.trim()).filter(Boolean).map(Number);
  if (shape.some((k) => !Number.isInteger(k) || k < 0)) fail('The .npy shape is invalid.');
  const count = shape.reduce((a, k) => a * k, 1);
  return { shape, fortran: !!fm && fm[1] === 'True', dtype: dm[1], isBool: d[2] === 'b', data: typed(u8, hoff + hlen, count, dt, d[1] !== '>') };
}
function npyGeom(a, opts) {
  const { shape, fortran, data } = a, nd = shape.length, stats = { shape, dtype: a.dtype, order: fortran ? 'Fortran' : 'C' };
  if (nd === 3) { const [nx, ny, nz] = fortran ? shape : [shape[2], shape[1], shape[0]]; return voxelsFromValues(data, nx, ny, nz, opts.spacing || [1, 1, 1], opts, { stats, warnings: [`Array axes read as ${fortran ? '(x, y, z)' : '(z, y, x)'}.`] }); }
  if (nd === 2 && (shape[1] === 2 || shape[1] === 3) && !a.isBool && (shape[0] > 4 || data.some((v) => v !== Math.floor(v)))) {
    const n = shape[0], w = shape[1], xyz = new Float64Array(n * 3);
    for (let i = 0; i < n; i++) for (let c = 0; c < w; c++) xyz[3 * i + c] = fortran ? data[c * n + i] : data[i * w + c];
    return pointsGeom(xyz, { stats });
  }
  if (nd === 2) { const [nx, ny] = fortran ? shape : [shape[1], shape[0]]; return imageGeom(data, nx, ny, 1, opts.spacing || [1, 1, 1], [0, 0, 0], opts, { stats }); }
  if (nd === 1) return { kind: 'table', headers: ['value'], records: Array.from(data.subarray(0, L.rows), (v) => ({ value: v })), stats };
  return fail(`NumPy arrays with ${nd} dimensions are not supported (use 2-D or 3-D arrays).`);
}
async function readNPY(ctx) {
  let u8 = await ctx.bytes();
  if (isZip(u8)) {
    const f = await unzipFiles(await ctx.buf(), (n) => /\.npy$/i.test(n)), first = f.entries().next().value;
    if (!first) fail('The .npz archive holds no arrays.');
    const g = npyGeom(parseNPY(first[1]), ctx.opts);
    g.stats = { ...(g.stats || {}), array: first[0].replace(/\.npy$/i, ''), arrays: f.size };
    if (f.size > 1) (g.warnings || (g.warnings = [])).push(`The archive holds ${f.size} arrays; the first ("${g.stats.array}") is used.`);
    return g;
  }
  return npyGeom(parseNPY(u8), ctx.opts);
}

const DCM_LONG = new Set(['OB', 'OW', 'OF', 'SQ', 'UT', 'UN', 'OD', 'OL', 'UC', 'UR', 'OV', 'SV', 'UV']);
async function readDICOM(ctx) {
  const u8 = await ctx.bytes(), dv = new DataView(u8.buffer, u8.byteOffset, u8.byteLength), conv = fmtNamed('DICOM slice').convert;
  if (u8.length < 140 || latin1.decode(u8.subarray(128, 132)) !== 'DICM') fail('Not a DICOM Part-10 file (DICM marker missing).');
  let explicit = true;
  const head = (p) => {                          // -> [group, element, value offset, length, vr]
    const g = dv.getUint16(p, true), e = dv.getUint16(p + 2, true);
    if (g === 0xfffe) return [g, e, p + 8, dv.getUint32(p + 4, true), ''];
    if (g === 2 || explicit) { const vr = latin1.decode(u8.subarray(p + 4, p + 6)); return DCM_LONG.has(vr) ? [g, e, p + 12, dv.getUint32(p + 8, true), vr] : [g, e, p + 8, dv.getUint16(p + 6, true), vr]; }
    return [g, e, p + 8, dv.getUint32(p + 4, true), ''];
  };
  const skipUndef = (p, depth) => {              // skip a sequence / item of undefined length; returns the offset after its delimiter
    if (depth > 24) fail('The DICOM file is nested too deeply.');
    while (p + 8 <= u8.length) {
      const [g, e, s, len] = head(p);
      if (g === 0xfffe && (e === 0xe0dd || e === 0xe00d)) return s;
      p = len === 0xffffffff ? skipUndef(s, depth + 1) : s + len;
    }
    return u8.length;
  };
  const str = (s, len) => latin1.decode(u8.subarray(s, s + len)).replace(/[\0 ]+$/, '').trim(), t = {};
  let p = 132, pix = null;
  while (p + 8 <= u8.length) {
    const [g, e, s, len] = head(p), key = g * 65536 + e;
    if (len === 0xffffffff) { if (key === 0x7fe00010) fail('The DICOM pixel data is compressed (encapsulated). ' + conv); p = skipUndef(s, 0); continue; }
    if (s + len > u8.length) { if (key === 0x7fe00010) fail('The DICOM file is truncated.'); break; }
    if (key === 0x00020010) { const ts = str(s, len); if (ts === '1.2.840.10008.1.2') explicit = false; else if (ts !== '1.2.840.10008.1.2.1') fail(`DICOM transfer syntax ${ts} (compressed or big-endian) is not supported. ` + conv); }
    else if (key === 0x00280010) t.rows = dv.getUint16(s, true); else if (key === 0x00280011) t.cols = dv.getUint16(s, true);
    else if (key === 0x00280100) t.bits = dv.getUint16(s, true); else if (key === 0x00280103) t.signed = dv.getUint16(s, true) === 1;
    else if (key === 0x00280002) t.spp = dv.getUint16(s, true); else if (key === 0x00280008) t.frames = parseInt(str(s, len), 10);
    else if (key === 0x00280030) t.ps = str(s, len).split('\\').map(Number); else if (key === 0x00180050) t.thick = parseFloat(str(s, len));
    else if (key === 0x00281052) t.inter = parseFloat(str(s, len)); else if (key === 0x00281053) t.slope = parseFloat(str(s, len));
    else if (key === 0x00280004) t.photo = str(s, len);
    else if (key === 0x7fe00010) { pix = [s, len]; break; }
    p = s + len;
  }
  const nz = t.frames > 0 ? t.frames : 1, dt = { 8: t.signed ? 'i1' : 'u1', 16: t.signed ? 'i2' : 'u2', 32: t.signed ? 'i4' : 'u4' }[t.bits];
  if (!pix || !(t.rows > 0 && t.cols > 0) || !dt) fail('The DICOM file has no readable uncompressed pixel data. ' + conv);
  if ((t.spp || 1) !== 1) fail('Colour DICOM images are not supported. ' + conv);
  if (t.rows * t.cols * nz > L.voxels) fail('The DICOM image is too large.');
  let vals = flipRows(typed(u8, pix[0], t.rows * t.cols * nz, dt, true), t.cols, t.rows, nz);
  const slope = Number.isFinite(t.slope) ? t.slope : 1, inter = Number.isFinite(t.inter) ? t.inter : 0, warnings = [];
  if (slope !== 1 || inter !== 0) vals = Float32Array.from(vals, (v) => v * slope + inter);
  if (t.photo === 'MONOCHROME1') warnings.push('MONOCHROME1 image: low values are bright on display; thresholding uses the stored values.');
  const ps = t.ps && t.ps.length === 2 && t.ps.every((v) => v > 0) ? t.ps : [1, 1];
  return voxelsFromValues(vals, t.cols, t.rows, nz, [ps[1], ps[0], t.thick > 0 ? t.thick : 1], ctx.opts, { warnings, stats: { units: t.ps ? 'mm' : undefined, frames: nz, bits: t.bits } });
}

// ---- Plant / piping / network topology ---------------------------------------------------------------------
const lcKeys = (o) => { const m = Object.create(null); for (const k of Object.keys(o)) m[k.toLowerCase().replace(/[_\-\s]/g, '')] = o[k]; return m; };
const firstOf = (m, keys) => { for (const k of keys) if (m[k] !== undefined && m[k] !== null && m[k] !== '') return m[k]; return undefined; };
const toNum = (v) => { if (typeof v === 'number') return Number.isFinite(v) ? v : null; if (typeof v === 'string' && v.trim() !== '') { const f = parseFloat(v); return Number.isFinite(f) ? f : null; } return null; };
const NODE_KEYS = ['nodes', 'units', 'equipment', 'equipments', 'components', 'unitops', 'blocks', 'vertices', 'items'], EDGE_KEYS = ['edges', 'links', 'pipes', 'streams', 'connections', 'lines', 'segments', 'arcs', 'piping'];
const FROM_KEYS = ['from', 'source', 'upstream', 'src', 'start', 'inlet', 'fromnode', 'node1', 'u'], TO_KEYS = ['to', 'target', 'downstream', 'dst', 'end', 'outlet', 'tonode', 'node2', 'v'];
function netGeom(nodes, edges, extra = {}) {
  if (!nodes.length) fail('The network definition holds no nodes or connections.');
  const lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
  for (const n of nodes) [n.x, n.y, n.z].forEach((v, k) => { if (typeof v === 'number') { lo[k] = Math.min(lo[k], v); hi[k] = Math.max(hi[k], v); } });
  for (let k = 0; k < 3; k++) if (lo[k] > hi[k]) lo[k] = hi[k] = 0;
  return { kind: 'network', network: { nodes, edges }, bbox: { min: lo, max: hi }, ...extra, stats: { nodes: nodes.length, edges: edges.length, ...(extra.stats || {}) } };
}
/** Liberal node / edge records -> normalised network Geometry. */
function buildNetwork(rawNodes, rawEdges, warnings = [], extra = {}) {
  const nodes = [], byId = new Map();
  const node = (id) => { id = String(id); let n = byId.get(id); if (!n) { n = { id, type: 'node', name: id, x: null, y: null, z: null }; byId.set(id, n); nodes.push(n); if (nodes.length > 2e5) fail('The network has too many nodes.'); } return n; };
  const list = Array.isArray(rawNodes) ? rawNodes.map((r, k) => [k + 1, r]) : rawNodes && typeof rawNodes === 'object' ? Object.entries(rawNodes) : [];
  for (const [key, r] of list) {
    if (r === null || typeof r !== 'object') { node(Array.isArray(rawNodes) ? r : key); continue; }
    const m = lcKeys(r), n = node(firstOf(m, ['id', 'tag', 'name', 'key', 'label']) ?? key), pos = firstOf(m, ['position', 'pos', 'coords', 'coordinates', 'xyz', 'location']);
    const p = Array.isArray(pos) ? pos : [];
    n.type = String(firstOf(m, ['type', 'kind', 'class', 'category']) ?? 'node').toLowerCase();
    n.name = String(firstOf(m, ['name', 'label', 'tag', 'description']) ?? n.id);
    n.x = toNum(m.x ?? p[0]); n.y = toNum(m.y ?? p[1]); n.z = toNum(firstOf(m, ['z', 'elevation', 'elev', 'height']) ?? p[2]);
  }
  const declared = nodes.length, edges = [];
  let bad = 0, dn = 0, computed = 0;
  for (const r of Array.isArray(rawEdges) ? rawEdges : []) {
    const m = Array.isArray(r) ? { from: r[0], to: r[1], length: r[2], diameter: r[3] } : r && typeof r === 'object' ? lcKeys(r) : null;
    let a = m && firstOf(m, FROM_KEYS), b = m && firstOf(m, TO_KEYS);
    if (a && typeof a === 'object') a = firstOf(lcKeys(a), ['id', 'tag', 'name']);
    if (b && typeof b === 'object') b = firstOf(lcKeys(b), ['id', 'tag', 'name']);
    if (a === undefined || b === undefined) { bad++; continue; }
    const na = node(a), nb = node(b);
    let length = toNum(firstOf(m, ['length', 'l', 'len'])), diameter = toNum(firstOf(m, ['diameter', 'd', 'dia', 'bore', 'innerdiameter']));
    if (diameter === null && toNum(m.dn) !== null) { diameter = toNum(m.dn) / 1000; dn++; }
    if (length === null && [na.x, na.y, nb.x, nb.y].every((v) => v !== null)) { length = Math.hypot(nb.x - na.x, nb.y - na.y, (nb.z || 0) - (na.z || 0)); computed++; }
    edges.push({ from: na.id, to: nb.id, type: String(firstOf(m, ['type', 'kind', 'class']) ?? 'pipe').toLowerCase(), name: String(firstOf(m, ['name', 'id', 'tag', 'label']) ?? `${na.id}-${nb.id}`), length, diameter });
    if (edges.length > 5e5) fail('The network has too many connections.');
  }
  if (!edges.length && !declared) fail('No nodes or connections were recognised in the network definition.');
  if (nodes.length > declared && declared) warnings.push(`${nodes.length - declared} nodes are referenced by connections but not declared; they were added.`);
  if (bad) warnings.push(`${bad} connections without a from / to pair were skipped.`);
  if (dn) warnings.push(`${dn} nominal diameters (DN, millimetres) were converted to metres.`);
  if (computed) warnings.push(`${computed} connection lengths were computed from node coordinates.`);
  return netGeom(nodes, edges, { ...extra, warnings });
}
function findNet(obj) {
  if (!obj || typeof obj !== 'object' || Array.isArray(obj)) return null;
  const m = lcKeys(obj), nk = NODE_KEYS.find((k) => m[k] && typeof m[k] === 'object'), ek = EDGE_KEYS.find((k) => Array.isArray(m[k]));
  if (ek && m[ek].some((e) => (Array.isArray(e) && e.length >= 2) || (e && typeof e === 'object' && firstOf(lcKeys(e), FROM_KEYS) !== undefined))) return { nodes: nk ? m[nk] : [], edges: m[ek] };
  for (const k of ['network', 'plant', 'model', 'graph', 'topology', 'flowsheet', 'process']) { const r = findNet(m[k]); if (r) return r; }
  return null;
}

/** Small safe YAML subset: block maps and lists, scalars, inline [..] and {..}; no anchors, tags or block scalars. */
function parseYAML(text) {
  const lines = [];
  for (const raw of splitLines(text.replace(/^﻿/, ''))) {
    if (/^\s*(#.*)?$/.test(raw) || /^(---|\.\.\.)(\s|$)/.test(raw) || /^%/.test(raw)) continue;
    const ind = raw.length - raw.trimStart().length;
    if (raw.slice(0, ind).includes('\t')) fail('YAML indentation must use spaces, not tabs.');
    let s = raw.trim(), q = '';
    for (let k = 0; k < s.length; k++) { const c = s[k]; if (q) { if (c === q) q = ''; } else if (c === '"' || c === "'") q = c; else if (c === '#' && (k === 0 || s[k - 1] === ' ')) { s = s.slice(0, k).trim(); break; } }
    if (s) lines.push({ ind, s });
    if (lines.length > 5e5) fail('The YAML file is too long.');
  }
  if (!lines.length) fail('The YAML file is empty.');
  const plain = (t) => {
    if (t === '' || t === '~' || /^null$/i.test(t)) return null;
    if (/^(true|yes)$/i.test(t)) return true;
    if (/^(false|no)$/i.test(t)) return false;
    if (/^[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$/.test(t)) return Number(t);
    if (/^[&*!|>]/.test(t)) fail('YAML anchors, aliases, tags and block scalars are not supported.');
    return t;
  };
  const flow = (s) => {
    let p = 0;
    const ws = () => { while (s[p] === ' ') p++; };
    const quoted = () => { const q = s[p], e = s.indexOf(q, p + 1); if (e < 0) fail('A YAML string is not closed.'); const v = s.slice(p + 1, e); p = e + 1; return v; };
    const val = (depth) => {
      ws();
      if (depth > 40) fail('YAML inline collections are nested too deeply.');
      if (s[p] === '[') { p++; const a = []; ws(); if (s[p] === ']') { p++; return a; } for (;;) { a.push(val(depth + 1)); ws(); if (s[p] === ',') { p++; continue; } if (s[p] === ']') { p++; return a; } fail('A YAML inline list is malformed.'); } }
      if (s[p] === '{') {
        p++;
        const o = {};
        ws();
        if (s[p] === '}') { p++; return o; }
        for (;;) {
          ws();
          let k;
          if (s[p] === '"' || s[p] === "'") k = quoted(); else { const e = s.indexOf(':', p); if (e < 0) fail('A YAML inline map is malformed.'); k = s.slice(p, e).trim(); p = e; }
          ws();
          if (s[p] !== ':') fail('A YAML inline map is malformed.');
          p++;
          const v = val(depth + 1);
          if (!SKIP_KEYS.has(k)) o[k] = v;
          ws();
          if (s[p] === ',') { p++; continue; }
          if (s[p] === '}') { p++; return o; }
          fail('A YAML inline map is malformed.');
        }
      }
      if (s[p] === '"' || s[p] === "'") return quoted();
      let e = p;
      while (e < s.length && !',]}'.includes(s[e])) e++;
      const t = s.slice(p, e).trim();
      p = e;
      return plain(t);
    };
    const out = val(0);
    ws();
    if (p < s.length) fail('A YAML inline value has trailing text.');
    return out;
  };
  const scalar = (t) => (t[0] === '[' || t[0] === '{' ? flow(t) : (t[0] === '"' || t[0] === "'") && t[t.length - 1] === t[0] && t.length > 1 ? t.slice(1, -1) : plain(t));
  const keyRe = /^("[^"]*"|'[^']*'|[^\s:#[\]{},"'&*!|>-][^:]*?|-[^\s:][^:]*?)\s*:(?:\s+(.*))?$/;
  let i = 0;
  const block = (ind, depth) => {
    if (depth > 60) fail('The YAML document is nested too deeply.');
    if (/^-(\s|$)/.test(lines[i].s)) {
      const arr = [];
      while (i < lines.length && lines[i].ind === ind && /^-(\s|$)/.test(lines[i].s)) {
        const s = lines[i].s, rest = s.slice(1).trim();
        if (!rest) { i++; arr.push(i < lines.length && lines[i].ind > ind ? block(lines[i].ind, depth + 1) : null); }
        else if (keyRe.test(rest) || /^-(\s|$)/.test(rest)) { const col = ind + s.length - rest.length; lines[i] = { ind: col, s: rest }; arr.push(block(col, depth + 1)); }
        else { arr.push(scalar(rest)); i++; }
      }
      return arr;
    }
    const obj = {};
    while (i < lines.length && lines[i].ind === ind && !/^-(\s|$)/.test(lines[i].s)) {
      const m = keyRe.exec(lines[i].s);
      if (!m) { if (depth === 0 && i === 0 && lines.length === 1) return scalar(lines[i++].s); fail(`YAML line not understood: "${lines[i].s.slice(0, 60)}".`); }
      const k = m[1].replace(/^(["'])(.*)\1$/, '$2'), v = (m[2] || '').trim();
      i++;
      let val;
      if (v !== '') val = scalar(v);
      else if (i < lines.length && (lines[i].ind > ind || (lines[i].ind === ind && /^-(\s|$)/.test(lines[i].s)))) val = block(lines[i].ind, depth + 1);
      else val = null;
      if (!SKIP_KEYS.has(k)) obj[k] = val;
    }
    return obj;
  };
  const root = block(lines[0].ind, 0);
  if (i < lines.length) fail(`YAML indentation is inconsistent near "${lines[i].s.slice(0, 60)}".`);
  return root;
}

async function readPCF(ctx) {
  const UNIT = { MM: 1e-3, CM: 0.01, M: 1, METRE: 1, METER: 1, INCH: 0.0254, INCHES: 0.0254, FEET: 0.3048, FT: 0.3048, 'MM-HUNDREDTHS': 1e-5 }, comps = [], warnings = [];
  let cu = 1e-3, bu = 1e-3, cuName = 'MM', buName = 'MM', cur = null, pipeline = null, seen = false;
  for (const raw of splitLines(await ctx.text())) {
    if (!raw.trim()) continue;
    const t = raw.trim().split(/\s+/), key = t[0].toUpperCase();
    const pt = () => { const v = t.slice(1, 5).map(Number); return v.slice(0, 3).every(Number.isFinite) ? { p: [v[0] * cu, v[1] * cu, v[2] * cu], bore: Number.isFinite(v[3]) ? v[3] * bu : null } : null; };
    if (!/^\s/.test(raw)) {
      if (key === 'UNITS-CO-ORDS') { cuName = (t[1] || 'MM').toUpperCase(); cu = UNIT[cuName] || 1e-3; seen = true; }
      else if (key === 'UNITS-BORE') { buName = (t[1] || 'MM').toUpperCase(); bu = UNIT[buName] || 1e-3; seen = true; }
      else if (key.startsWith('UNITS-') || key === 'ISOGEN-FILES') { seen = true; cur = null; }
      else if (key === 'PIPELINE-REFERENCE') { pipeline = t.slice(1).join(' '); cur = null; seen = true; }
      else if (key === 'MATERIALS') break;
      else { cur = { type: key, ends: [], branch: [], centre: null, co: null, at: {} }; comps.push(cur); if (comps.length > 5e5) fail('The PCF file holds too many components.'); }
    } else if (cur) {
      const p = /POINT$|^CO-ORDS$/.test(key) ? pt() : null;
      if (key === 'END-POINT' && p) cur.ends.push(p);
      else if ((key === 'CENTRE-POINT' || key === 'CENTER-POINT') && p) cur.centre = p;
      else if (/^BRANCH\d*-POINT$/.test(key) && p) cur.branch.push(p);
      else if (key === 'CO-ORDS' && p) cur.co = p;
      else if (['ITEM-CODE', 'SKEY', 'NAME', 'TAG', 'COMPONENT-IDENTIFIER'].includes(key)) cur.at[key] = t.slice(1).join(' ');
    }
  }
  const real = comps.filter((c) => c.ends.length || c.co);
  if (!real.length) fail('No piping components with END-POINT data were found: not a PCF file.');
  const nodes = [], edges = [], byKey = new Map(), tol = 1e-3, count = {};
  const node = (p, type) => {
    const k = p.p.map((v) => Math.round(v / tol)).join(',');
    let n = byKey.get(k);
    if (!n) { n = { id: 'N' + (nodes.length + 1), type: type || 'junction', name: 'N' + (nodes.length + 1), x: p.p[0], y: p.p[1], z: p.p[2] }; byKey.set(k, n); nodes.push(n); }
    else if (type && n.type === 'junction') n.type = type;
    return n;
  };
  const dist = (a, b) => vlen(sub(a.p, b.p));
  for (const c of real) {
    const type = c.type.toLowerCase().split('-')[0];
    count[type] = (count[type] || 0) + 1;
    const name = c.at.TAG || c.at.NAME || c.at['COMPONENT-IDENTIFIER'] || `${c.type}-${count[type]}`, [a, b] = c.ends;
    const edge = (p, q, length) => edges.push({ from: node(p).id, to: node(q).id, type, name, length, diameter: p.bore ?? q.bore ?? null, ...(c.at['ITEM-CODE'] ? { itemCode: c.at['ITEM-CODE'] } : {}) });
    if (a && b) {
      if (c.centre && c.branch.length) { c.centre.bore = a.bore; node(c.centre, 'tee'); edge(a, c.centre, dist(a, c.centre)); edge(c.centre, b, dist(c.centre, b)); for (const br of c.branch) edge(c.centre, br, dist(c.centre, br)); }
      else if (c.centre) {
        const la = dist(a, c.centre), lb = dist(c.centre, b), cosv = la && lb ? dot(sub(a.p, c.centre.p), sub(b.p, c.centre.p)) / (la * lb) : -1, phi = Math.acos(Math.max(-1, Math.min(1, cosv)));
        edge(a, b, phi > 1e-6 && phi < Math.PI - 1e-6 ? la * Math.tan(phi / 2) * (Math.PI - phi) : la + lb);
      } else edge(a, b, dist(a, b));
      for (const extra of c.ends.slice(2)) edge(a, extra, dist(a, extra));
    } else node(a || c.co, type).name = name;
  }
  if (!seen) warnings.push('No UNITS lines were found; millimetres are assumed.');
  warnings.push(`Coordinates (${cuName}) and bores (${buName}) converted to metres; nodes within 1 mm are merged.`);
  return netGeom(nodes, edges, { warnings, stats: { units: 'm', components: count, pipeline } });
}

const IFC_KIND = [[/^IFCPUMP$/, 'pump'], [/^IFCTANK$/, 'tank'], [/^IFCVALVE$/, 'valve'], [/^IFC(PIPE|FLOW)SEGMENT$/, 'pipe'], [/^IFC(PIPE|FLOW)FITTING$/, 'fitting'], [/^IFCHEATEXCHANGER$/, 'heat exchanger'], [/^IFCEVAPORATOR$/, 'evaporator'], [/^IFCCONDENSER$/, 'condenser'], [/^IFCFILTER$/, 'filter'], [/^IFCFLOWMETER$/, 'meter'], [/^IFCBOILER$/, 'boiler'], [/^IFCCOOLINGTOWER$/, 'cooling tower'], [/^IFCCOMPRESSOR$/, 'compressor'], [/^IFCFAN$/, 'fan'], [/^IFCINTERCEPTOR$/, 'interceptor'], [/^IFCELECTRICMOTOR$/, 'motor'], [/^IFCFLOWSTORAGEDEVICE$/, 'tank'], [/^IFCFLOWMOVINGDEVICE$/, 'pump'], [/^IFCFLOWCONTROLLER$/, 'valve'], [/^IFCFLOWTREATMENTDEVICE$/, 'filter'], [/^IFC(FLOW|SANITARY|WASTE)TERMINAL$/, 'terminal'], [/^IFCENERGYCONVERSIONDEVICE$/, 'energy conversion'], [/^IFCDISTRIBUTIONCHAMBERELEMENT$/, 'chamber'], [/^IFCUNITARYEQUIPMENT$/, 'equipment'], [/^IFCTUBEBUNDLE$/, 'tube bundle'], [/^IFCCHILLER$/, 'chiller'], [/^IFCBURNER$/, 'burner'], [/^IFCENGINE$/, 'engine'], [/^IFCSTACKTERMINAL$/, 'terminal']];
async function readIFC(ctx) {
  const db = parseP21(await ctx.text()), { get, ids, ref, byType } = db, warnings = [];
  if (![...byType.keys()].some((k) => k.startsWith('IFC'))) fail('Not an IFC file (no IFC entities found).');
  const en = (v) => (v && v.e) || '', numOf = (v) => (typeof v === 'number' ? v : v && Array.isArray(v.a) && typeof v.a[0] === 'number' ? v.a[0] : null), listOf = (v) => (Array.isArray(v) ? v : []);
  let lu = 1, unitName = 'm';
  const units = ids('IFCUNITASSIGNMENT').flatMap((id) => listOf((get(id) || {}).args && get(id).args[0]).map(ref)).filter(Boolean);
  for (const u of units.length ? units : ids('IFCSIUNIT').map(get)) {
    if (!u || !u.args || en(u.args[1]) !== 'LENGTHUNIT') continue;
    if (u.type === 'IFCSIUNIT') { lu = en(u.args[2]) ? P21_PREFIX[en(u.args[2])] || 1 : 1; unitName = UNIT_NAME[lu] || `${lu} m`; }
    else if (u.type === 'IFCCONVERSIONBASEDUNIT') { const nm = String(u.args[2] || '').toUpperCase(); lu = /INCH/.test(nm) ? 0.0254 : /F(OO|EE)T/.test(nm) ? 0.3048 : 1; unitName = nm.toLowerCase(); }
    break;
  }
  const vec = (v, T) => { const e = ref(v); const c = e && e.type === T && Array.isArray(e.args[0]) ? e.args[0] : null; return c && c.every((q) => typeof q === 'number') ? [c[0], c[1] || 0, c[2] || 0] : null; };
  const I3 = { o: [0, 0, 0], x: [1, 0, 0], y: [0, 1, 0], z: [0, 0, 1] }, pcache = new Map();
  const axis = (v) => {
    const e = ref(v);
    if (!e || e.type !== 'IFCAXIS2PLACEMENT3D') return I3;
    const o = vec(e.args[0], 'IFCCARTESIANPOINT') || [0, 0, 0], z = unit(vec(e.args[1], 'IFCDIRECTION') || [0, 0, 1]), xr = vec(e.args[2], 'IFCDIRECTION') || [1, 0, 0], k = dot(xr, z);
    let x = unit([xr[0] - k * z[0], xr[1] - k * z[1], xr[2] - k * z[2]]);
    if (!vlen(x)) x = frameFromNormal(z).x;
    return { o, x, y: cross(z, x), z };
  };
  const place = (v, depth = 0) => {
    const e = ref(v);
    if (!e || e.type !== 'IFCLOCALPLACEMENT' || depth > 64) return I3;
    if (pcache.has(e.id)) return pcache.get(e.id);
    const P = place(e.args[0], depth + 1), A = axis(e.args[1]), rot = (q) => [P.x[0] * q[0] + P.y[0] * q[1] + P.z[0] * q[2], P.x[1] * q[0] + P.y[1] * q[1] + P.z[1] * q[2], P.x[2] * q[0] + P.y[2] * q[1] + P.z[2] * q[2]];
    const ro = rot(A.o), out = { o: [P.o[0] + ro[0], P.o[1] + ro[1], P.o[2] + ro[2]], x: rot(A.x), y: rot(A.y), z: rot(A.z) };
    pcache.set(e.id, out);
    return out;
  };
  // lengths from quantity sets
  const qLen = new Map();
  for (const id of ids('IFCRELDEFINESBYPROPERTIES')) {
    const a = (get(id) || {}).args, q = a && ref(a[5]);
    if (!q || q.type !== 'IFCELEMENTQUANTITY') continue;
    for (const qr of listOf(q.args[5])) { const ql = ref(qr); if (ql && ql.type === 'IFCQUANTITYLENGTH' && /length/i.test(String(ql.args[0])) && numOf(ql.args[3]) !== null) { for (const o of listOf(a[4])) if (o && o.r !== undefined && !qLen.has(o.r)) qLen.set(o.r, numOf(ql.args[3]) * lu); break; } }
  }
  const shapeOf = (v) => {                       // -> { length, diameter } from extrusions / swept discs
    const out = { length: null, diameter: null }, items = [];
    const reps = (pd, depth) => { const e = ref(pd); if (!e || depth > 4) return; if (e.type === 'IFCPRODUCTDEFINITIONSHAPE') for (const r of listOf(e.args[2])) reps(r, depth + 1); else if (e.type === 'IFCSHAPEREPRESENTATION') for (const it of listOf(e.args[3])) { const ie = ref(it); if (ie && ie.type === 'IFCMAPPEDITEM') { const src = ref(ie.args[0]); if (src) reps(src.args[1], depth + 1); } else if (ie) items.push(ie); } };
    reps(v, 0);
    for (const it of items) {
      if (it.type === 'IFCEXTRUDEDAREASOLID') {
        const d = numOf(it.args[3]), pr = ref(it.args[0]);
        if (d !== null && (out.length === null || d * lu > out.length)) out.length = d * lu;
        if (pr && /^IFCCIRCLE(HOLLOW)?PROFILEDEF$/.test(pr.type) && numOf(pr.args[3]) !== null) out.diameter = 2 * numOf(pr.args[3]) * lu;
      } else if (it.type === 'IFCSWEPTDISKSOLID') {
        const dx = ref(it.args[0]), r = numOf(it.args[1]);
        if (r !== null) out.diameter = 2 * r * lu;
        if (dx && dx.type === 'IFCPOLYLINE') { const ps = listOf(dx.args[0]).map((q) => vec(q, 'IFCCARTESIANPOINT')).filter(Boolean); let s = 0; for (let k = 1; k < ps.length; k++) s += vlen(sub(ps[k], ps[k - 1])); if (s > 0) out.length = s * lu; }
      }
    }
    return out;
  };
  const nodes = [], byEnt = new Map();
  let pipes = 0, pipesWithLength = 0, placed = 0;
  const addNode = (id, kind) => {
    if (byEnt.has(id)) return byEnt.get(id);
    const e = get(id), a = e && e.args;
    if (!a) return null;
    const pl = a[5] ? place(a[5]) : null, n = { id: typeof a[0] === 'string' && a[0] ? a[0] : '#' + id, type: kind || e.type.replace(/^IFC/, '').toLowerCase(), name: String(a[2] || a[7] || `${e.type} #${id}`), x: null, y: null, z: null, ifcClass: e.type, step: id };
    if (pl && ref(a[5])) { n.x = pl.o[0] * lu; n.y = pl.o[1] * lu; n.z = pl.o[2] * lu; placed++; }
    if (typeof a[7] === 'string' && a[7]) n.tag = a[7];
    if (n.type === 'pipe') {
      pipes++;
      const sh = shapeOf(a[6]);
      n.length = qLen.get(id) ?? sh.length; n.diameter = sh.diameter;
      if (n.length !== null && n.length !== undefined) pipesWithLength++; else n.length = null;
    }
    byEnt.set(id, n); nodes.push(n);
    if (nodes.length > 2e5) fail('The IFC model holds too many distribution elements.');
    return n;
  };
  for (const [type, list] of byType) { const k = IFC_KIND.find((q) => q[0].test(type)); if (k) for (const id of list) addNode(id, k[1]); }
  // port -> owning element
  const owner = new Map();
  for (const id of ids('IFCRELCONNECTSPORTTOELEMENT')) { const a = (get(id) || {}).args; if (a && a[4] && a[5] && a[4].r !== undefined && a[5].r !== undefined) owner.set(a[4].r, a[5].r); }
  for (const id of ids('IFCRELNESTS')) { const a = (get(id) || {}).args; if (a && a[4] && a[4].r !== undefined) for (const o of listOf(a[5])) { const pe = ref(o); if (pe && pe.type === 'IFCDISTRIBUTIONPORT') owner.set(o.r, a[4].r); } }
  const edges = [];
  let dangling = 0;
  const flowDir = (pid) => { const e = get(pid); return e && e.args ? e.args.map(en).find((s) => s === 'SOURCE' || s === 'SINK' || s === 'SOURCEANDSINK') || '' : ''; };
  for (const id of ids('IFCRELCONNECTSPORTS')) {
    const a = (get(id) || {}).args;
    if (!a || !a[4] || !a[5] || a[4].r === undefined || a[5].r === undefined) continue;
    const ea = owner.get(a[4].r), eb = owner.get(a[5].r), na = ea !== undefined ? addNode(ea) : null, nb = eb !== undefined ? addNode(eb) : null;
    if (!na || !nb) { dangling++; continue; }
    const rev = flowDir(a[4].r) === 'SINK' && flowDir(a[5].r) !== 'SINK', [f, t] = rev ? [nb, na] : [na, nb];
    edges.push({ from: f.id, to: t.id, type: 'connection', name: `${f.name} → ${t.name}`, length: null, diameter: f.diameter ?? t.diameter ?? null });
  }
  if (!nodes.length) fail('The IFC file holds no pumps, tanks, valves, pipe segments or other distribution elements.');
  const kinds = {};
  for (const n of nodes) kinds[n.type] = (kinds[n.type] || 0) + 1;
  warnings.push(`Recovered ${nodes.length} elements (${Object.entries(kinds).map(([k, v]) => `${v} ${k}`).join(', ')}); ${placed} with a placement.`);
  warnings.push(edges.length ? `${edges.length} port connections recovered${dangling ? `; ${dangling} connections to elements outside the inventory were dropped` : ''}.` : 'No port connectivity (IfcRelConnectsPorts) was found: the result is an equipment inventory without connections.');
  if (pipes) warnings.push(`Pipe-segment lengths recovered for ${pipesWithLength} of ${pipes} segments (from length quantities or extrusion depth).`);
  warnings.push('Solid geometry is not tessellated; for shape export STL or OBJ from the BIM tool.');
  if (db.stat.malformed) warnings.push(`${db.stat.malformed} malformed entities were ignored.`);
  return netGeom(nodes, edges, { warnings, stats: { schema: db.schema, units: 'm', sourceLengthUnit: unitName, byType: kinds, pipes, pipesWithLength } });
}

function amlGeom(doc) {
  const root = xAll(doc, 'caexfile')[0];
  if (!root) fail('Not an AutomationML / CAEX file (no <CAEXFile> element).');
  const nodes = [], edges = [], owner = new Map(), strip = (s) => String(s || '').replace(/[{}]/g, ''), last = (s) => String(s || '').split('/').pop();
  const walk = (el, parent, depth) => {
    if (depth > 200) return;
    for (const ie of xKids(el, 'internalelement')) {
      const id = strip(ie.attrs.id) || ie.attrs.name || 'E' + (nodes.length + 1), rr = xKid(ie, 'rolerequirements'), sr = xKid(ie, 'supportedroleclass');
      const type = last((rr && rr.attrs.refbaseroleclasspath) || (sr && sr.attrs.refroleclasspath) || ie.attrs.refbasesystemunitpath) || 'element';
      nodes.push({ id, type: type.toLowerCase(), name: ie.attrs.name || id, x: null, y: null, z: null, parent });
      if (nodes.length > 2e5) fail('The AutomationML file holds too many elements.');
      for (const xi of xKids(ie, 'externalinterface')) { if (xi.attrs.id) owner.set(strip(xi.attrs.id), id); if (xi.attrs.name) owner.set(id + ':' + xi.attrs.name, id); }
      walk(ie, id, depth + 1);
    }
  };
  for (const ih of xAll(root, 'instancehierarchy')) walk(ih, null, 0);
  let bad = 0;
  for (const l of xAll(root, 'internallink')) {
    const res = (s) => { const k = strip(s); return owner.get(k) ?? owner.get(k.split(':')[0]) ?? (nodes.some((n) => n.id === k.split(':')[0]) ? k.split(':')[0] : undefined); };
    const a = res(l.attrs.refpartnersidea), b = res(l.attrs.refpartnersideb);
    if (a === undefined || b === undefined) { bad++; continue; }
    edges.push({ from: a, to: b, type: 'link', name: l.attrs.name || `${a}-${b}`, length: null, diameter: null });
  }
  return netGeom(nodes, edges, { warnings: bad ? [`${bad} InternalLinks could not be resolved to elements and were skipped.`] : [], stats: { hierarchies: xAll(root, 'instancehierarchy').length } });
}
function xmlNetwork(doc) {
  const obj = (el) => { const o = {}; for (const k in el.attrs) o[k] = el.attrs[k]; for (const c of el.children) if (!c.children.length && c.text.trim() && !(c.name in o)) o[c.name] = c.text.trim(); return o; };
  const NT = ['node', 'unit', 'equipment', 'component', 'unitop', 'block', 'vertex'], nodes = xAll(doc, ...NT).filter((el) => !xAll(el, ...NT).length && (Object.keys(el.attrs).length || el.children.length)).map(obj), edges = xAll(doc, 'edge', 'link', 'pipe', 'stream', 'connection', 'segment', 'arc').map(obj);
  if (!edges.some((e) => firstOf(lcKeys(e), FROM_KEYS) !== undefined)) fail('The XML file is not a recognised geometry or network definition (expected nodes / units / equipment and edges / links / pipes / streams with from and to).');
  return buildNetwork(nodes, edges, []);
}

// ---- Numeric sniffing: JSON / YAML objects, XML, delimited text -----------------------------------------------
const isNum = (v) => typeof v === 'number' && Number.isFinite(v);
/** Nested numeric arrays -> { shape, data } (row-major) or null when ragged / non-numeric. */
function nestedArray(a) {
  const shape = [];
  for (let q = a; Array.isArray(q); q = q[0]) { shape.push(q.length); if (shape.length > 3) return null; }
  if (!shape.length || shape.some((n) => !n)) return null;
  const n = shape.reduce((p, q) => p * q, 1);
  if (n > L.voxels) fail('The array is too large.');
  const data = new Float64Array(n);
  let o = 0, ok = true;
  const fill = (q, d) => { if (!ok) return; if (!Array.isArray(q) || q.length !== shape[d]) { ok = false; return; } if (d === shape.length - 1) for (const v of q) { if (typeof v === 'boolean') data[o++] = +v; else if (isNum(v)) data[o++] = v; else if (v === null) data[o++] = NaN; else { ok = false; return; } } else for (const s of q) fill(s, d + 1); };
  fill(a, 0);
  return ok ? { shape, data } : null;
}
function jsonMesh(rawNodes, rawElems, key, warnings) {
  const nn = rawNodes.length, xyz = new Float64Array(nn * 3), store = cellStore();
  let dim = 2;
  for (let i = 0; i < nn; i++) { const p = Array.isArray(rawNodes[i]) ? rawNodes[i] : rawNodes[i] && typeof rawNodes[i] === 'object' ? [rawNodes[i].x, rawNodes[i].y, rawNodes[i].z] : []; if (!isNum(p[0]) || !isNum(p[1])) fail('Mesh node coordinates must be numeric [x, y, z] arrays.'); xyz[3 * i] = p[0]; xyz[3 * i + 1] = p[1]; if (isNum(p[2])) { xyz[3 * i + 2] = p[2]; if (p[2] !== 0) dim = 3; } }
  const conn = rawElems.map((e) => (Array.isArray(e) ? e : e && typeof e === 'object' ? firstOf(lcKeys(e), ['nodes', 'n', 'conn', 'connectivity', 'vertices', 'v']) : null));
  let lo = Infinity, hi = -Infinity, skipped = 0;
  for (const c of conn) if (Array.isArray(c)) for (const v of c) { if (v < lo) lo = v; if (v > hi) hi = v; }
  const base = lo >= 1 && hi === nn ? 1 : 0, P = (i) => [xyz[3 * i], xyz[3 * i + 1], xyz[3 * i + 2]];
  if (base) warnings.push('Element connectivity is 1-based and was shifted to 0-based.');
  for (const c of conn) {
    if (!Array.isArray(c) || c.some((v) => !Number.isInteger(v))) { skipped++; continue; }
    const ids = c.map((v) => v - base), n = ids.length;
    if (ids.some((v) => v < 0 || v >= nn)) fail('The mesh connectivity refers to nodes that do not exist.');
    let type = n === 2 ? 'line' : n === 3 ? 'tri' : n === 8 ? 'hex' : n === 6 ? 'wedge' : n === 5 ? 'pyramid' : n === 4 ? 'quad' : n > 8 ? 'poly' : null;
    if (n === 4 && dim === 3 && !/face|tri|quad|poly/i.test(key)) { const [a, b, cc, d] = ids.map(P), vol = Math.abs(dot(cross(sub(b, a), sub(cc, a)), sub(d, a))), s = vlen(sub(b, a)) * vlen(sub(cc, a)) * vlen(sub(d, a)); if (vol > 1e-6 * s) type = 'tet'; }
    if ((n === 5 || n === 6 || n === 8) && (dim === 2 || /face|poly/i.test(key))) type = 'poly';
    if (!type) { skipped++; continue; }
    store.add(type, ids);
  }
  if (skipped) warnings.push(`${skipped} elements were malformed or of unsupported size and were ignored.`);
  return cellMesh(xyz, store, { warnings });
}
/** Generic object (from JSON or YAML) -> network, mesh, voxels / grid, points, polylines or table. */
function objectGeom(j, opts) {
  const warnings = [], net = findNet(j);
  if (net) return buildNetwork(net.nodes, net.edges, warnings);
  const m = j && typeof j === 'object' && !Array.isArray(j) ? lcKeys(j) : null;
  if (m) {
    const nk = ['nodes', 'vertices', 'points', 'coordinates', 'coords'].find((k) => Array.isArray(m[k]) && m[k].length), ek = ['elements', 'cells', 'faces', 'triangles', 'tris', 'connectivity', 'polygons'].find((k) => Array.isArray(m[k]) && m[k].length);
    if (nk && ek) return jsonMesh(m[nk], m[ek], ek, warnings);
    if (Array.isArray(m.polylines) && m.polylines.length) return polyGeom(m.polylines.map((p) => { const q = p && typeof p === 'object' ? lcKeys(p) : {}; return Array.isArray(p) ? { x: p.map((v) => +v[0]), y: p.map((v) => +v[1]), closed: false } : { x: (q.x || []).map(Number), y: (q.y || []).map(Number), closed: !!q.closed }; }), { warnings });
    if (Array.isArray(m.x) && Array.isArray(m.y) && m.x.length === m.y.length && m.x.length > 1 && m.x.every(isNum)) {
      if (Array.isArray(m.z) && Array.isArray(m.z[0])) { const z = nestedArray(m.z); if (z && z.shape.length === 2 && z.shape[0] === m.y.length && z.shape[1] === m.x.length) return gridGeom(z.data, m.x.length, m.y.length, (i) => m.x[i], (k) => m.y[k], { nodata: (v) => v !== v, warnings }); }
      if (Array.isArray(m.z) && m.z.length === m.x.length) { const xyz = new Float64Array(m.x.length * 3); m.x.forEach((v, i) => { xyz[3 * i] = v; xyz[3 * i + 1] = +m.y[i]; xyz[3 * i + 2] = +m.z[i]; }); return pointsGeom(xyz, { warnings }); }
      const x = m.x.slice(), y = m.y.map(Number);
      return polyGeom([{ x, y, closed: closeRing(x, y) || !!m.closed }], { warnings });
    }
  }
  const ak = m ? ['voxels', 'mask', 'labels', 'data', 'grid', 'values', 'array', 'points', 'z', 'elevation', 'sdf', 'levelset'].find((k) => Array.isArray(m[k]) && Array.isArray(m[k][0])) : null, arr = Array.isArray(j) ? j : ak ? m[ak] : null;
  if (arr && Array.isArray(arr[0])) {
    const a = nestedArray(arr);
    if (!a) fail('The array is ragged or holds non-numeric values.');
    const sp = (m && (m.spacing || m.voxelsize)) || opts.spacing || [1, 1, 1], spacing = Array.isArray(sp) ? sp : [sp, sp, sp], sdf = ak === 'sdf' || ak === 'levelset';
    if (a.shape.length === 3) { const g = voxelsFromValues(sdf ? a.data.map((v) => -v) : a.data, a.shape[2], a.shape[1], a.shape[0], spacing, sdf ? { ...opts, threshold: opts.threshold ?? 0 } : opts, { warnings }); if (sdf) g.warnings.push('Signed-distance / level-set field: negative values taken as solid.'); return g; }
    const [n0, n1] = a.shape;
    if ((n1 === 2 || n1 === 3) && ak !== 'voxels' && ak !== 'mask' && ak !== 'grid' && (n0 > 4 || ak === 'points' || a.data.some((v) => v !== Math.floor(v)))) {
      if (n1 === 2 && ak !== 'points') { const x = [], y = []; for (let i = 0; i < n0; i++) { x.push(a.data[2 * i]); y.push(a.data[2 * i + 1]); } return polyGeom([{ x, y, closed: closeRing(x, y) }], { warnings }); }
      const xyz = new Float64Array(n0 * 3);
      for (let i = 0; i < n0; i++) for (let c = 0; c < n1; c++) xyz[3 * i + c] = a.data[i * n1 + c];
      return pointsGeom(xyz, { warnings });
    }
    return imageGeom(sdf ? a.data.map((v) => -v) : a.data, n1, n0, 1, spacing, [0, 0, 0], sdf ? { ...opts, threshold: opts.threshold ?? 0 } : opts, { warnings });
  }
  const recs = Array.isArray(j) ? j : m && Array.isArray(m.records) ? m.records : m && Array.isArray(m.data) ? m.data : null;
  if (recs && recs.length && recs.every((r) => r && typeof r === 'object' && !Array.isArray(r))) {
    const headers = [...new Set(recs.slice(0, 200).flatMap((r) => Object.keys(r)))].filter((k) => !SKIP_KEYS.has(k));
    if (headers.some((h) => FROM_KEYS.includes(normHeader(h))) && headers.some((h) => TO_KEYS.includes(normHeader(h)))) return buildNetwork([], recs, warnings);
    const tc = tableColumns(headers), records = recs.slice(0, L.rows).map((r) => Object.fromEntries(headers.map((h) => [h, r[h] ?? ''])));
    if (tc.role === 'survey') return surveyFromRecords(headers, records, tc, warnings);
    return { kind: 'table', headers, records, warnings, stats: { rows: recs.length, ...(tc.role ? { role: tc.role } : {}) } };
  }
  if (m && Object.keys(j).length) {             // flat parameter set, e.g. a parametric membrane or spacer definition
    const flat = {}, put = (o, pre, depth) => { for (const [k, v] of Object.entries(o)) { if (v && typeof v === 'object' && !Array.isArray(v) && depth < 4) put(v, pre + k + '.', depth + 1); else flat[pre + k] = Array.isArray(v) ? JSON.stringify(v).slice(0, 200) : v; } };
    put(j, '', 0);
    warnings.push('The file holds parameters rather than shape data; it is returned as a one-row parameter table.');
    return { kind: 'table', headers: Object.keys(flat), records: [flat], params: j, warnings, stats: { parameters: Object.keys(flat).length } };
  }
  return fail('The file holds no nodes/edges, mesh, array, point list or record table that can be used as geometry.');
}
async function readJSON(ctx) {
  const text = await ctx.text(), j = safeJSON(text);
  if (j && typeof j === 'object' && !Array.isArray(j)) {
    if (typeof j.type === 'string' && /^(FeatureCollection|Feature|GeometryCollection|(Multi)?(Polygon|LineString|Point))$/.test(j.type)) return geojsonGeom(text, j);
    if (j.asset && Array.isArray(j.meshes)) return Object.assign(await readGLTF(ctx, j), { format: 'glTF 2.0' });
  }
  return objectGeom(j, ctx.opts);
}
function geojsonGeom(text, j) {
  const pts = [], warnings = [];
  const walk = (o, depth) => {
    if (!o || typeof o !== 'object' || depth > 32) return;
    if (o.type === 'FeatureCollection') (o.features || []).forEach((f) => walk(f, depth + 1));
    else if (o.type === 'Feature') walk(o.geometry, depth + 1);
    else if (o.type === 'GeometryCollection') (o.geometries || []).forEach((q) => walk(q, depth + 1));
    else if (o.type === 'Point' && Array.isArray(o.coordinates)) pts.push(+o.coordinates[0], +o.coordinates[1], +o.coordinates[2] || 0);
    else if (o.type === 'MultiPoint' && Array.isArray(o.coordinates)) for (const c of o.coordinates) if (Array.isArray(c)) pts.push(+c[0], +c[1], +c[2] || 0);
  };
  walk(j, 0);
  let g = null;
  try { g = parseGeoJSON(text); } catch { g = null; }
  if (!g && !pts.length) fail('No polygon, line or point geometry was found in the GeoJSON file.');
  const polys = g ? g.polylines.map((p) => { const x = p.x.slice(), y = p.y.slice(); if (p.closed) closeRing(x, y); return { x, y, closed: p.closed }; }) : [];
  const out = linesOrPoints(polys, pts, { warnings, format: 'GeoJSON', pathway: 'gis' });
  out.geographic = out.bbox.min[0] >= -180 && out.bbox.max[0] <= 180 && out.bbox.min[1] >= -90 && out.bbox.max[1] <= 90;
  return out;
}
async function readYAML(ctx) { return objectGeom(parseYAML(await ctx.text()), ctx.opts); }
async function readXML(ctx) {
  const text = await ctx.text(), head = text.slice(0, 4000);
  if (/<VTKFile/i.test(head)) return Object.assign(await readVTKXML(ctx), { format: 'VTK XML' });
  const doc = parseXML(text), root = doc.children[0].name, as = (g, name) => Object.assign(g, { format: name });
  if (root === 'landxml') return as(landxmlGeom(doc, ctx.opts), 'LandXML');
  if (root === 'graphml') return as(graphmlGeom(doc), 'GraphML');
  if (root === 'kml') return as(kmlGeom(doc), 'KML / KMZ');
  if (root === 'gpx') return as(await readGPX(ctx, doc), 'GPX');
  if (root === 'svg') return as(await readSVG(ctx, doc), 'SVG');
  if (root === 'collada') return as(await readDAE(ctx, doc), 'COLLADA');
  if (root === 'x3d') return as(await readX3D(ctx, doc), 'X3D');
  if (root === 'caexfile') return as(amlGeom(doc), 'AutomationML / CAEX');
  if (root === 'amf') return as(await readAMF(ctx), 'AMF');
  if (xAll(doc, 'poslist', 'linearring', 'linestring').length) return as(await readGML(ctx, doc), 'GML');
  return xmlNetwork(doc);
}
/** Ordered points that trace one line (a route or centreline) rather than a scattered cloud: short hop-to-hop path. */
function pathLike(xyz) {
  const n = Math.floor(xyz.length / 3);
  if (n < 3 || n > 2e5) return false;
  const bb = bboxOf(xyz), diag = Math.hypot(bb.max[0] - bb.min[0], bb.max[1] - bb.min[1]);
  let sum = 0;
  for (let i = 1; i < n && sum <= 3 * diag; i++) sum += Math.hypot(xyz[3 * i] - xyz[3 * i - 3], xyz[3 * i + 1] - xyz[3 * i - 2]);
  return diag > 0 && sum <= 3 * diag;
}
/**
 * CSV / TSV / TXT / DAT / XYZ / ASC text: node-edge table, well survey, profile / route / wall-map table (stats.role),
 * particle cloud, soundings, points, polyline, matrix or plain table.
 */
async function readDelimited(ctx) {
  const text = await ctx.text(), tab = readColumns(text), warnings = [], hc = headerCells(text);
  if (hc && tab.ncol && hc.length === tab.ncol) tab.headers = hc;
  const asTable = () => {
    const t = tableFromRows(parseDelimited(text, ctx.ext === 'tsv' ? '\t' : undefined));
    if (!t.records.length) fail('The file holds no rows of data.');
    return t;
  };
  const heads = hc || tab.headers, tc = heads ? tableColumns(heads) : null, role = tc ? tc.role : null;
  if (heads) {
    const hn = heads.map(normHeader);
    if (hn.some((h) => FROM_KEYS.includes(h)) && hn.some((h) => TO_KEYS.includes(h))) return Object.assign(buildNetwork([], asTable().records, warnings), { pathway: 'network' });
    if (role === 'survey') { const sv = surveyFromText(text, true); if (sv) return sv; }
  }
  if (!tab.rows || !tab.ncol) fail('The file holds no numeric rows.');
  const tbl = (t) => ({ kind: 'table', ...t, warnings, stats: { rows: t.records.length, ...(role ? { role } : {}) } });
  const numTable = () => {            // numeric rows keyed by header; also covers whitespace-separated tables
    if (!tab.headers || tab.textCells) return asTable();
    const records = [], nc = tab.ncol;
    for (let i = 0; i < Math.min(tab.rows, L.rows); i++) { const r = {}; for (let c = 0; c < nc; c++) r[tab.headers[c]] = tab.data[i * nc + c]; records.push(r); }
    if (tab.rows > L.rows) warnings.push(`Only the first ${L.rows} of ${tab.rows} rows are kept.`);
    return { headers: tab.headers.slice(), records };
  };
  if (role === 'thicknessMap' || role === 'depositMap' || role === 'profile' || (role === 'route3d' && tc.c.chainage !== undefined)) return tbl(numTable());
  const roles = columnRoles(tab.headers);
  if (role === 'particles' && roles && !tab.textCells) {
    const { xyz } = tableToXYZ(tab, warnings), g = pointsGeom(xyz, { warnings, stats: { columns: tab.headers, role } });
    if (g.count === tab.rows) { g.attributes = {}; tab.headers.forEach((h, c) => { if (![roles.x, roles.y, roles.z, roles.depth].includes(c)) g.attributes[h] = Float64Array.from({ length: tab.rows }, (_, i) => tab.data[i * tab.ncol + c]); }); }
    else warnings.push('Per-particle attributes were dropped because rows were removed or sub-sampled.');
    return Object.assign(g, { pathway: 'points' });
  }
  if (tab.textCells > 0.05 * tab.rows * tab.ncol && !roles) return tbl(asTable());
  if (roles && (roles.z >= 0 || roles.depth >= 0)) {
    const { xyz } = tableToXYZ(tab, warnings), ex = { warnings, geographic: roles.geo, stats: { columns: tab.headers } };
    if (roles.geo || roles.depth >= 0) {
      const lg = latticeGrid(xyz, ex);
      if (lg) return Object.assign(lg, { pathway: 'gis' });
      if (pathLike(xyz)) { warnings.push('The rows trace a single line, so they are kept as an ordered 3-D route rather than gridded as soundings.'); return Object.assign(pointsGeom(xyz, { ...ex, stats: { ...ex.stats, role: 'route3d' } }), { pathway: 'gis' }); }
      const pg = pointsGeom(xyz, { warnings }), n = pg.count, side = Math.max(8, Math.min(512, Math.round(Math.sqrt(n)))), w = pg.bbox.max[0] - pg.bbox.min[0], h = pg.bbox.max[1] - pg.bbox.min[1];
      const nx = Math.max(4, Math.round(side * Math.sqrt((w || 1) / (h || 1)))), ny = Math.max(4, Math.round((side * side) / nx)), gr = scatterGrid(pg.points, nx, ny);
      warnings.push(`${n} scattered soundings binned onto a ${nx} × ${ny} grid (empty cells filled by inverse-distance weighting); the original points are kept in "points".`);
      return { kind: 'grid', grid: { ...gr, nodata: 0, geographic: roles.geo }, points: pg.points, count: n, bbox: pg.bbox, pathway: 'gis', warnings, stats: { nx, ny, points: n, columns: tab.headers } };
    }
    if (pathLike(xyz)) ex.stats.role = 'route3d';
    const g = pointsGeom(xyz, ex);
    return Object.assign(g, { pathway: 'points' });
  }
  if (tab.headers && !roles) return tbl(numTable());
  if (tab.ncol === 2) {
    const cuts = [0, ...tab.breaks, tab.rows], polys = [];
    for (let k = 0; k + 1 < cuts.length; k++) { const x = [], y = []; for (let i = cuts[k]; i < cuts[k + 1]; i++) { x.push(tab.data[2 * i]); y.push(tab.data[2 * i + 1]); } if (x.length > 1) polys.push({ x, y, closed: closeRing(x, y) || !!ctx.opts.closed }); }
    return Object.assign(polyGeom(polys, { warnings }), { pathway: 'drawing' });
  }
  if (tab.ncol >= 3 && tab.ncol <= 7) {
    if (tab.ncol > 3) warnings.push(`Columns beyond x, y, z (${tab.ncol - 3} per row: intensity, colour or radius) are ignored.`);
    const { xyz } = tableToXYZ(tab, warnings);
    return Object.assign(pointsGeom(xyz, { warnings, stats: { columns: tab.headers || undefined, ...(pathLike(xyz) ? { role: 'route3d' } : {}) } }), { pathway: 'points' });
  }
  if (tab.ncol >= 8 && tab.rows >= 2 && !tab.textCells) { warnings.push(`Numeric matrix of ${tab.rows} rows × ${tab.ncol} columns read as a regular array (row 1 at y = 0).`); return imageGeom(tab.data, tab.ncol, tab.rows, 1, ctx.opts.spacing || [1, 1, 1], [0, 0, 0], ctx.opts, { warnings }); }
  return tbl(asTable());
}

// ---- Table roles: profiles, routes, well surveys, wall maps, particles --------------------------------------
const LEN_UNIT = { m: 1, km: 1000, mm: 1e-3, cm: 0.01, ft: 0.3048, feet: 0.3048, usft: 1200 / 3937, in: 0.0254, inch: 0.0254 };
/** Unit named in a column header such as "KP (km)", "MD [ft]" or "wt_mm"; null when none is given. */
const headerUnit = (h) => { const m = /[[(]\s*([a-zµ°%/]+)\s*[\])]/i.exec(String(h)) || /[_\s](mm|km|cm|m|ft|usft|in|deg|rad|s|h|hr|min|d)$/i.exec(String(h)); return m ? m[1].toLowerCase() : null; };
const COL = {
  md: /^(md|mdepth|measureddepth|measdepth|alongholedepth|ahd|dept|mdrkb|mdkb|mdbrt|mdft)$/,
  inc: /^(inc|incl|inclin|inclination|incdeg|incldeg|deviation|devi|dev|drift|holeangle)$/,
  azi: /^(azi|azim|azimuth|az|azideg|azimdeg|hazi|azimuthtrue|azimuthgrid|azitrue|azigrid|bearing|direction)$/,
  tvd: /^(tvd|tvdrkb|tvdkb|tvdbrt|tvdss|tvdmsl|trueverticaldepth|verticaldepth|tvdft)$/,
  north: /^(ns|north|northing|n|dn|deltan|nsoffset|northoffset|localn|ynorth|dy|y|yoffset)$/,
  east: /^(ew|east|easting|e|de|deltae|ewoffset|eastoffset|locale|xeast|dx|x|xoffset)$/,
  chainage: /^(chainage|ch|kp|kilometrepoint|kilometerpoint|station|sta|stationing|distance|dist|horizontaldistance|pipelinedistance|routedistance|alongroute|arclength)$/,
  lon: /^(lon|long|longitude|lng)$/, lat: /^(lat|latitude)$/,
  x: ROLE.x, y: ROLE.y, z: ROLE.z,
  depth: /^(depth|d|sounding|soundings|bathy|bathymetry|waterdepth|wd|seabeddepth)$/,
  theta: /^(theta|angle|ang|clock|clockposition|clockpos|oclock|circumferential|circ|circumferentialposition|orientation|phi)$/,
  thickness: /^(wt|wallthickness|thickness|thk|wallthk|remainingwall|remainingthickness|remainingwt|twall|tmeas|measuredthickness)$/,
  loss: /^(metalloss|wallloss|corrosion|corrosiondepth|pitdepth|defectdepth|lossdepth|dentdepth|loss)$/,
  deposit: /^(deposit|depositthickness|delta|hydrate|hydratethickness|wax|waxthickness|scale|scalethickness|depositheight|layerthickness|film|filmthickness)$/,
  time: /^(t|time|timestamp|date|datetime|elapsed|elapsedtime|times|hours|days|seconds|minutes)$/,
  diameter: /^(dp|diameter|diam|size|particlediameter|particlesize|radius|agglomeratesize|d50)$/,
  velocity: /^(u|v|w|vx|vy|vz|velocity|speed|density|rho|hydratefraction|massfraction)$/,
};
/**
 * Column roles of a header row. Returns { c: { md, inc, azi, tvd, north, east, chainage, lon, lat, x, y, z, depth, theta,
 * thickness, loss, deposit, time, diameter, velocity → column index }, unit: { same keys → unit text or null }, role } where
 * role is 'survey' | 'thicknessMap' | 'depositMap' | 'particles' | 'profile' | 'route3d' | 'timeseries' | null.
 */
export function tableColumns(headers) {
  const hn = (headers || []).map(normHeader), c = {}, unit = {};
  for (const k of Object.keys(COL)) { const i = hn.findIndex((h) => COL[k].test(h)); if (i >= 0) { c[k] = i; unit[k] = headerUnit(headers[i]); } }
  if (c.md === undefined && c.inc !== undefined && c.depth !== undefined) { c.md = c.depth; unit.md = unit.depth; }
  const has = (k) => c[k] !== undefined, value = has('deposit') || has('thickness') || has('loss'), axial = has('chainage') || has('x');
  let role = null;
  if (has('md') && (has('inc') || has('tvd'))) role = 'survey';
  else if (value && axial && (has('theta') || (has('y') && !has('z')))) role = has('deposit') ? 'depositMap' : 'thicknessMap';
  else if (has('x') && has('y') && (has('diameter') || has('velocity'))) role = 'particles';
  else if (has('chainage') && (has('z') || has('depth')) && !(has('x') && has('y'))) role = 'profile';
  else if (has('x') && has('y') && (has('z') || has('depth'))) role = 'route3d';
  else if (has('time') && hn.length > 1) role = 'timeseries';
  return { c, unit, role };
}
/** Header cells of the first non-comment line when that line is a header (holds a non-numeric cell), else null. */
function headerCells(text) {
  for (let pos = text.charCodeAt(0) === 0xfeff ? 1 : 0, k = 0; pos < text.length && k < 200; k++) {
    let e = text.indexOf('\n', pos);
    if (e < 0) e = text.length;
    const l = text.slice(pos, e).trim();
    pos = e + 1;
    if (!l || /^(#|%|!|\/\/)/.test(l)) continue;
    const cells = (/[,;\t]/.test(l) ? l.split(/[,;\t]/) : l.split(/\s+/)).map((s) => s.trim().replace(/^"|"$/g, '').trim());
    return cells.some((s) => s !== '' && !Number.isFinite(+s)) ? cells : null;
  }
  return null;
}

// ---- Well deviation surveys ----------------------------------------------------------------------------------
/**
 * Minimum-curvature positions of a directional survey. stations: [{ md, inc, azi }] (or [md, inc, azi] triples, or
 * { md[], inc[], azi[] }) in metres and degrees; start: { tvd, north, east } of the first station (default tvd = md of the
 * first station, i.e. a vertical hole above it). Returns { md, inc, azi, tvd, north, east, dogleg (deg per interval),
 * dls (deg / 30 m) } as arrays; dls[0] = dogleg[0] = 0.
 */
export function minimumCurvature(stations, start = {}) {
  const src = stations && !Array.isArray(stations) && Array.isArray(stations.md) ? stations.md.map((m, i) => [m, stations.inc ? stations.inc[i] : 0, stations.azi ? stations.azi[i] : 0]) : stations;
  if (!Array.isArray(src) || !src.length) throw new Error('minimumCurvature needs survey stations with md, inc and azi.');
  const md = [], inc = [], azi = [], D = Math.PI / 180;
  for (const q of src) {
    const a = Array.isArray(q) ? q : q && typeof q === 'object' ? [q.md, q.inc ?? q.incl ?? q.inclination, q.azi ?? q.azim ?? q.azimuth ?? 0] : [];
    if (![a[0], a[1], a[2]].every((v) => typeof v === 'number' && Number.isFinite(v))) throw new Error('Survey stations must hold numeric md, inc and azi.');
    if (md.length && a[0] < md[md.length - 1]) throw new Error('Measured depth must not decrease along the survey.');
    md.push(a[0]); inc.push(a[1]); azi.push(a[2]);
  }
  const fin0 = (v, d) => (typeof v === 'number' && Number.isFinite(v) ? v : d);
  const tvd = [fin0(start.tvd, md[0])], north = [fin0(start.north, 0)], east = [fin0(start.east, 0)], dogleg = [0], dls = [0];
  for (let i = 1; i < md.length; i++) {
    const dm = md[i] - md[i - 1], i1 = inc[i - 1] * D, i2 = inc[i] * D, a1 = azi[i - 1] * D, a2 = azi[i] * D;
    const dl = Math.acos(Math.max(-1, Math.min(1, Math.cos(i2 - i1) - Math.sin(i1) * Math.sin(i2) * (1 - Math.cos(a2 - a1))))), rf = dl < 1e-6 ? 1 + (dl * dl) / 12 : (2 / dl) * Math.tan(dl / 2);
    north.push(north[i - 1] + 0.5 * dm * (Math.sin(i1) * Math.cos(a1) + Math.sin(i2) * Math.cos(a2)) * rf);
    east.push(east[i - 1] + 0.5 * dm * (Math.sin(i1) * Math.sin(a1) + Math.sin(i2) * Math.sin(a2)) * rf);
    tvd.push(tvd[i - 1] + 0.5 * dm * (Math.cos(i1) + Math.cos(i2)) * rf);
    dogleg.push(dl / D); dls.push(dm > 0 ? ((dl / D) * 30) / dm : 0);
  }
  return { md, inc, azi, tvd, north, east, dogleg, dls };
}
/** Survey columns (metres, degrees; tvd / north / east optional) -> 3-D well path: x east, y north, z elevation (negative down). */
function surveyGeom(col, extra = {}) {
  const warnings = extra.warnings || [], rows = [], n0 = col.md.length, opt = ['inc', 'azi', 'tvd', 'north', 'east'].filter((k) => col[k]);
  let dropped = 0;
  for (let i = 0; i < n0; i++) {
    const r = { md: col.md[i] };
    for (const k of opt) r[k] = col[k][i];
    if (!Number.isFinite(r.md) || opt.some((k) => !Number.isFinite(r[k])) || (rows.length && r.md <= rows[rows.length - 1].md)) { dropped++; continue; }
    rows.push(r);
    if (rows.length > L.rows) fail('The survey holds too many stations.');
  }
  if (rows.length < 2) fail('A well survey needs at least two stations with increasing measured depth.');
  if (dropped) warnings.push(`${dropped} survey rows with missing values or non-increasing measured depth were dropped.`);
  const md = rows.map((r) => r.md), D = Math.PI / 180;
  let s;
  if (col.inc) {
    if (!col.azi) warnings.push('The survey has no azimuth column: the well is laid out in one vertical plane heading north.');
    if (rows.some((r) => r.inc < 0 || r.inc > 180)) fail('Survey inclinations must lie between 0 and 180 degrees.');
    s = minimumCurvature(rows.map((r) => [r.md, r.inc, col.azi ? r.azi : 0]), { tvd: col.tvd ? rows[0].tvd : rows[0].md, north: col.north ? rows[0].north : 0, east: col.east ? rows[0].east : 0 });
    if (!col.tvd && rows[0].md > 0) warnings.push(`The first station is at MD ${+rows[0].md.toFixed(3)}; the hole above it is taken as vertical.`);
    if (col.tvd) { const e = Math.abs(s.tvd[s.tvd.length - 1] - rows[rows.length - 1].tvd); if (e > 0.01 * Math.max(1, md[md.length - 1] - md[0])) warnings.push(`TVD recomputed by minimum curvature differs from the file's TVD column by ${+e.toFixed(2)} at the last station; the computed path is used.`); }
  } else {
    const tvd = rows.map((r) => r.tvd), ne = col.north && col.east, north = [ne ? rows[0].north : 0], east = [ne ? rows[0].east : 0], inc = [], azi = [];
    let bad = 0;
    for (let i = 1; i < rows.length; i++) {
      const dm = md[i] - md[i - 1], dz = tvd[i] - tvd[i - 1];
      if (Math.abs(dz) > dm * (1 + 1e-6)) bad++;
      if (ne) { north.push(rows[i].north); east.push(rows[i].east); } else { north.push(0); east.push(east[i - 1] + Math.sqrt(Math.max(0, dm * dm - dz * dz))); }
      inc.push(Math.acos(Math.max(-1, Math.min(1, dz / dm))) / D);
      azi.push(ne ? ((Math.atan2(east[i] - east[i - 1], north[i] - north[i - 1]) / D) + 360) % 360 : 90);
    }
    inc.unshift(inc[0]); azi.unshift(azi[0]);
    const dogleg = [0], dls = [0];
    for (let i = 1; i < rows.length; i++) { const i1 = inc[i - 1] * D, i2 = inc[i] * D, dl = Math.acos(Math.max(-1, Math.min(1, Math.cos(i2 - i1) - Math.sin(i1) * Math.sin(i2) * (1 - Math.cos((azi[i] - azi[i - 1]) * D))))) / D; dogleg.push(dl); dls.push((dl * 30) / (md[i] - md[i - 1])); }
    if (bad) warnings.push(`${bad} intervals have a TVD change larger than their MD change; check the depth columns.`);
    warnings.push(ne ? 'No inclination column: inclination, azimuth and dog-leg severity are derived from the MD, TVD and offset columns (interval averages).' : 'Only MD and TVD are given: the horizontal departure is laid out due east (azimuth unknown) and inclination is the interval average.');
    s = { md, inc, azi, tvd, north, east, dogleg, dls };
  }
  const g = polyGeom([{ x: s.east.slice(), y: s.north.slice(), z: s.tvd.map((v) => 0 - v), closed: false }], { ...extra, warnings });
  g.survey = { md: s.md, inc: s.inc, azi: s.azi, tvd: s.tvd, north: s.north, east: s.east, dls: s.dls };
  g.stats = { ...(extra.stats || {}), md: md[md.length - 1], tvd: Math.max(...s.tvd), maxInclination: Math.max(...s.inc), maxDLS: Math.max(...s.dls), stations: md.length, wellSurvey: true };
  g.pathway = 'well';
  return g;
}
/** Survey table from records keyed by header. */
function surveyFromRecords(headers, records, tc, warnings) {
  const col = {}, toM = {};
  for (const k of ['md', 'inc', 'azi', 'tvd', 'north', 'east']) if (tc.c[k] !== undefined) { const h = headers[tc.c[k]]; col[k] = records.map((r) => (typeof r[h] === 'number' ? r[h] : parseFloat(r[h]))); toM[k] = LEN_UNIT[tc.unit[k]] ?? null; }
  return surveyGeom(surveyUnits(col, toM, warnings), { warnings });
}
function surveyUnits(col, toM, warnings) {
  const f = toM.md ?? toM.tvd ?? 1;
  if (f !== 1) { for (const k of ['md', 'tvd', 'north', 'east']) if (col[k]) col[k] = col[k].map((v) => v * (toM[k] ?? f)); warnings.push(`Survey depths converted to metres (× ${+f.toPrecision(6)}).`); }
  if (!(col.north && col.east)) { delete col.north; delete col.east; }
  return col;
}
/**
 * Deviation-survey text: the last header line before the data that names MD + inclination (+ azimuth) or MD + TVD decides
 * the columns. Returns null when no such header exists and `strict` is set; otherwise a header-less file is MD, inc, azi.
 */
function surveyFromText(text, strict) {
  const ls = splitLines(text), warnings = [];
  let tc = null, feet = false, first = -1;
  for (let i = 0; i < ls.length && i < 500; i++) {
    const raw = ls[i].replace(/^[\s#%!;*/]+/, '').trim();
    if (!raw) continue;
    const v = numsOf(raw);
    if (v.length >= 2 && v.every(Number.isFinite)) { first = i; break; }
    if (/\b(ft|feet|foot)\b/i.test(raw) && /\b(unit|units|depth|md|tvd)\b/i.test(raw)) feet = true;
    for (const cells of [/[,;\t]/.test(raw) ? raw.split(/[,;\t]/) : raw.split(/\s{2,}/), raw.split(/[\s,;]+/)]) { const q = tableColumns(cells.map((s) => s.trim().replace(/^"|"$/g, ''))); if (q.role === 'survey') { tc = q; break; } }
  }
  if (first < 0) return strict ? null : fail('The survey file holds no numeric rows.');
  if (!tc) {
    if (strict) return null;
    const nc = numsOf(ls[first].trim()).length;
    tc = nc >= 3 ? { c: { md: 0, inc: 1, azi: 2 }, unit: {} } : { c: { md: 0, tvd: 1 }, unit: {} };
    warnings.push(nc >= 3 ? 'No column header found: columns taken as MD, inclination, azimuth.' : 'No column header found: the two columns are taken as MD and TVD.');
  }
  const keys = ['md', 'inc', 'azi', 'tvd', 'north', 'east'].filter((k) => tc.c[k] !== undefined), col = Object.fromEntries(keys.map((k) => [k, []])), need = Math.max(...keys.map((k) => tc.c[k])) + 1, toM = {};
  let skipped = 0;
  for (let i = first; i < ls.length; i++) {
    const l = ls[i].trim();
    if (!l || /^(#|%|!|\/\/)/.test(l)) continue;
    const v = numsOf(l);
    if (v.length < need) { skipped++; continue; }
    for (const k of keys) col[k].push(v[tc.c[k]]);
    if (col.md.length > L.rows) fail('The survey holds too many stations.');
  }
  for (const k of keys) toM[k] = LEN_UNIT[tc.unit[k]] ?? (feet ? 0.3048 : null);
  if (skipped) warnings.push(`${skipped} lines with too few columns were skipped.`);
  return surveyGeom(surveyUnits(col, toM, warnings), { warnings, format: 'Well deviation survey (ASCII)' });
}
async function readSurvey(ctx) { return surveyFromText(await ctx.text(), false); }

// ---- GraphML ---------------------------------------------------------------------------------------------------
function graphmlGeom(doc) {
  const keys = new Map(), warnings = [];
  for (const k of xAll(doc, 'key')) if (k.attrs.id) keys.set(k.attrs.id, String(k.attrs['attr.name'] ?? k.attrs.id).toLowerCase());
  const rec = (el) => { const o = {}; for (const d of xKids(el, 'data')) { const name = keys.get(d.attrs.key) ?? String(d.attrs.key ?? '').toLowerCase(), t = d.text.trim(); if (name && !SKIP_KEYS.has(name) && t !== '' && !d.children.length) o[name] = t; } return o; };
  const nodes = xAll(doc, 'node').map((n) => {
    const o = rec(n), yg = xAll(n, 'geometry')[0], lab = xAll(n, 'nodelabel')[0];
    if (yg && o.x === undefined && o.y === undefined) { o.x = yg.attrs.x; o.y = yg.attrs.y; }
    if (lab && o.name === undefined && lab.text.trim()) o.name = lab.text.trim();
    o.id = n.attrs.id;
    return o;
  }).filter((o) => o.id !== undefined);
  const edges = xAll(doc, 'edge').map((e) => { const o = rec(e); for (const k of [...FROM_KEYS, ...TO_KEYS]) delete o[k]; if (o.name === undefined && e.attrs.id !== undefined) o.name = e.attrs.id; return { ...o, from: e.attrs.source, to: e.attrs.target }; });
  const hyper = xAll(doc, 'hyperedge').length;
  if (hyper) warnings.push(`${hyper} hyperedges were skipped (only two-ended edges are read).`);
  if (!nodes.length && !edges.length) fail('The GraphML file holds no nodes or edges.');
  return buildNetwork(nodes, edges, warnings);
}

// ---- LandXML ---------------------------------------------------------------------------------------------------
function landxmlGeom(doc, opts = {}) {
  const root = doc.children[0], warnings = [], un = xAll(root, 'metric', 'imperial')[0], lu = un ? String(un.attrs.linearunit || '').toLowerCase() : '', du = un ? String(un.attrs.diameterunit || '').toLowerCase() : '';
  const units = /foot|feet/.test(lu) ? 'ft' : lu.startsWith('millim') ? 'mm' : lu.startsWith('kilom') ? 'km' : 'm', stats = { units };
  const named = new Map(), cg = [];
  const pt = (el) => {
    if (!el) return null;
    if (el.attrs.pntref !== undefined && named.has(el.attrs.pntref)) return named.get(el.attrs.pntref);
    const v = numsOf(el.text);
    return v.length >= 2 && Number.isFinite(v[0]) && Number.isFinite(v[1]) ? [v[1], v[0], v.length > 2 && Number.isFinite(v[2]) ? v[2] : NaN] : null;   // LandXML order is northing, easting, elevation
  };
  for (const p of xAll(root, 'cgpoint')) { const q = pt(p); if (q) { cg.push(q); if (p.attrs.name !== undefined) named.set(p.attrs.name, q); } }
  let spirals = 0, badGeom = 0;
  const chain = (coordGeom) => {
    const P = [], push = (q) => { if (!q) { badGeom++; return; } const l = P[P.length - 1]; if (!l || l[0] !== q[0] || l[1] !== q[1]) P.push(q); };
    for (const el of coordGeom.children) {
      const a = pt(xKid(el, 'start')), b = pt(xKid(el, 'end'));
      if (el.name === 'line') { push(a); push(b); }
      else if (el.name === 'irregularline') { push(a); const pl = xKid(el, 'pntlist3d') || xKid(el, 'pntlist2d'), dim = pl && pl.name === 'pntlist3d' ? 3 : 2, v = pl ? numsOf(pl.text) : []; for (let k = 0; k + dim <= v.length; k += dim) push([v[k + 1], v[k], dim === 3 ? v[k + 2] : NaN]); push(b); }
      else if (el.name === 'curve') {
        const c = pt(xKid(el, 'center'));
        if (!a || !b || !c) { badGeom++; continue; }
        const r = Math.hypot(a[0] - c[0], a[1] - c[1]), t0 = Math.atan2(a[1] - c[1], a[0] - c[0]), cw = String(el.attrs.rot || '').toLowerCase() === 'cw';
        let sw = Math.atan2(b[1] - c[1], b[0] - c[0]) - t0;
        if (cw && sw >= 0) sw -= TAU; else if (!cw && sw <= 0) sw += TAU;
        const n = Math.max(2, Math.ceil((Math.abs(sw) / TAU) * ARC_N)), za = a[2], zb = b[2];
        push(a);
        for (let k = 1; k < n; k++) { const t = t0 + (sw * k) / n; push([c[0] + r * Math.cos(t), c[1] + r * Math.sin(t), za + ((zb - za) * k) / n]); }
        push(b);
      } else if (el.name === 'spiral') { spirals++; push(a); push(b); }
    }
    return P;
  };
  const polys = [];
  let vcurves = 0, flat = 0;
  const addPoly = (P, prof, sta0) => {
    if (P.length < 2 || polys.length > 2e4) return;
    let x = P.map((q) => q[0]), y = P.map((q) => q[1]), z = P.every((q) => Number.isFinite(q[2])) ? P.map((q) => q[2]) : null;
    if (prof && prof.length >= 2) {                    // vertical alignment: station -> elevation, vertices added at every PVI
      const st = [sta0];
      for (let i = 1; i < x.length; i++) st.push(st[i - 1] + Math.hypot(x[i] - x[i - 1], y[i] - y[i - 1]));
      const all = [...new Set([...st, ...prof.map((q) => q[0]).filter((s) => s > st[0] && s < st[st.length - 1])])].sort((p, q) => p - q), lerp = (xs, ys, s) => { let k = 1; while (k < xs.length - 1 && xs[k] < s) k++; const t = xs[k] === xs[k - 1] ? 0 : Math.max(0, Math.min(1, (s - xs[k - 1]) / (xs[k] - xs[k - 1]))); return ys[k - 1] + t * (ys[k] - ys[k - 1]); };
      const ps = prof.map((q) => q[0]), pz = prof.map((q) => q[1]), nx = all.map((s) => lerp(st, x, s)), ny = all.map((s) => lerp(st, y, s));
      z = all.map((s) => lerp(ps, pz, s)); x = nx; y = ny;
    }
    if (!z) flat++;
    polys.push(z ? { x, y, z, closed: false } : { x, y, closed: false });
  };
  for (const al of xAll(root, 'alignment')) {
    const g = xKid(al, 'coordgeom');
    if (!g) continue;
    const pa = xAll(al, 'profalign')[0], ps = xAll(al, 'profsurf')[0];
    let prof = null;
    if (pa) { prof = []; for (const el of pa.children) { const v = numsOf(el.text); if (['pvi', 'paracurve', 'circcurve', 'unsymparacurve'].includes(el.name) && v.length >= 2) { prof.push([v[0], v[1]]); if (el.name !== 'pvi') vcurves++; } } }
    else if (ps) { const pl = xKid(ps, 'pntlist2d'), v = pl ? numsOf(pl.text) : []; prof = []; for (let k = 0; k + 2 <= v.length; k += 2) prof.push([v[k], v[k + 1]]); }
    if (prof) prof = prof.filter((q) => Number.isFinite(q[0]) && Number.isFinite(q[1])).sort((p, q) => p[0] - q[0]);
    addPoly(chain(g), prof, Number.isFinite(+al.attrs.stastart) ? +al.attrs.stastart : 0);
  }
  for (const pf of xAll(root, 'planfeature')) { const g = xKid(pf, 'coordgeom'); if (g) addPoly(chain(g), null, 0); }
  // pipe networks
  const structs = xAll(root, 'struct'), pipes = xAll(root, 'pipe'), dscale = du.startsWith('millim') ? 0.001 : du.startsWith('inch') ? 0.0254 : du.startsWith('centim') ? 0.01 : 1;
  // surfaces
  const tri = [];
  let surfaces = 0;
  const readSurfaces = () => {
    for (const sf of xAll(root, 'surface')) {
      const def = xKid(sf, 'definition');
      if (!def) continue;
      const P = new Map();
      for (const p of xAll(def, 'p')) { const q = pt(p); if (q && p.attrs.id !== undefined && Number.isFinite(q[2])) P.set(p.attrs.id, q); }
      let n = 0;
      for (const f of xAll(def, 'f')) {
        if (f.attrs.i === '1') continue;
        const ids = f.text.trim().split(/\s+/).map((k) => P.get(k));
        if (ids.length < 3 || ids.some((q) => !q)) { badGeom++; continue; }
        for (let k = 1; k + 1 < ids.length; k++) for (const q of [ids[0], ids[k], ids[k + 1]]) tri.push(q[0], q[1], q[2]);
        n++;
        if (tri.length > L.triangles * 9) fail(`The LandXML surface has too many triangles (limit ${L.triangles}).`);
      }
      if (n) surfaces++;
    }
  };
  const have = { alignment: polys.length > 0, network: pipes.length > 0 && structs.length > 0, surface: xAll(root, 'surface').length > 0, points: cg.length > 0 };
  const pref = ['alignment', 'network', 'surface', 'points'], want = String(opts.prefer || '').toLowerCase(), pick = have[want] ? want : pref.find((k) => have[k]);
  if (!pick) fail('The LandXML file holds no alignments, plan features, surfaces, pipe networks or survey points.');
  const other = pref.filter((k) => k !== pick && have[k]).map((k) => ({ alignment: `${polys.length} alignments / plan features`, network: `a pipe network of ${pipes.length} pipes`, surface: 'a TIN surface', points: `${cg.length} survey points` })[k]);
  if (other.length) warnings.push(`The file also holds ${other.join(', ')}; only the ${{ alignment: 'alignments and plan features', network: 'pipe network', surface: 'surface', points: 'survey points' }[pick]} are returned (importGeometry option prefer selects another class).`);
  if (badGeom) warnings.push(`${badGeom} geometry elements with missing or unresolved points were skipped.`);
  if (units !== 'm') warnings.push(`Coordinates are in ${units === 'ft' ? 'feet' : units}; they are kept as written (stats.units = "${units}").`);
  if (pick === 'alignment') {
    if (spirals) warnings.push(`${spirals} spirals were replaced by their chords.`);
    if (vcurves) warnings.push(`${vcurves} vertical curves were replaced by straight grades between their intersection points.`);
    if (flat) warnings.push(`${flat} of ${polys.length} lines carry no elevations (no Profile and no 3-D points).`);
    return polyGeom(polys, { warnings, stats: { ...stats, alignments: polys.length } });
  }
  if (pick === 'network') {
    const nodes = structs.map((s) => { const c = pt(xKid(s, 'center')), z = [s.attrs.elevsump, s.attrs.elevrim].map(Number).find(Number.isFinite); return { id: s.attrs.name ?? s.attrs.id, type: xKid(s, 'inlet') ? 'inlet' : xKid(s, 'outlet') ? 'outlet' : xKid(s, 'connection') ? 'connection' : 'structure', name: s.attrs.desc || s.attrs.name, x: c ? c[0] : null, y: c ? c[1] : null, z: z ?? (c && Number.isFinite(c[2]) ? c[2] : null) }; }).filter((n) => n.id !== undefined);
    const edges = pipes.map((p) => { const cp = xKid(p, 'circpipe') || xKid(p, 'ellippipe') || xKid(p, 'rectpipe'), d = cp ? +(cp.attrs.diameter ?? cp.attrs.span ?? cp.attrs.width) : NaN, o = { from: p.attrs.refstart, to: p.attrs.refend, name: p.attrs.name, type: 'pipe' }; if (Number.isFinite(+p.attrs.length)) o.length = +p.attrs.length; if (Number.isFinite(d)) o.diameter = d * dscale; return o; });
    return buildNetwork(nodes, edges, warnings, { stats });
  }
  if (pick === 'surface') {
    readSurfaces();
    if (!tri.length) fail('The LandXML surfaces hold no usable faces (Pnts + Faces are required).');
    if (surfaces > 1) warnings.push(`${surfaces} surfaces were merged into one mesh.`);
    return meshGeom(tri, { warnings, stats: { ...stats, surfaces, surfaceOnly: true } });
  }
  const xyz = new Float64Array(cg.length * 3);
  cg.forEach((q, i) => { xyz[3 * i] = q[0]; xyz[3 * i + 1] = q[1]; xyz[3 * i + 2] = Number.isFinite(q[2]) ? q[2] : 0; });
  return pointsGeom(xyz, { warnings, stats });
}
async function readLandXML(ctx) {
  const doc = parseXML(await ctx.text());
  if (doc.children[0].name !== 'landxml') fail('Not a LandXML file (root element LandXML is missing).');
  return landxmlGeom(doc, ctx.opts);
}

// ---- Structural decks: Abaqus, ANSYS CDB, LS-DYNA -----------------------------------------------------------------
/** Eight node ids of a (possibly collapsed) hexahedron -> [cell type, ids] or null. */
function degenerateSolid(v) {
  const u = [];
  for (const q of v) if (!u.includes(q)) u.push(q);
  if (u.length === 8) return ['hex', v.slice(0, 8)];
  if (u.length === 4) return ['tet', u];
  if (u.length === 5 && v[4] === v[5] && v[5] === v[6] && v[6] === v[7]) return ['pyramid', [v[0], v[1], v[2], v[3], v[4]]];
  if (u.length === 6 && v[2] === v[3] && v[6] === v[7]) return ['wedge', [v[0], v[1], v[2], v[4], v[5], v[6]]];
  if (u.length === 6 && v[4] === v[5] && v[6] === v[7]) return ['wedge', [v[0], v[1], v[4], v[3], v[2], v[6]]];
  return null;
}
const SOLID_N = { 4: ['tet', 4], 10: ['tet', 4], 8: ['hex', 8], 20: ['hex', 8], 27: ['hex', 8], 6: ['wedge', 6], 15: ['wedge', 6], 18: ['wedge', 6], 5: ['pyramid', 5], 13: ['pyramid', 5] };
/** Abaqus element type name -> { type, n (corner nodes used), total (nodes listed), mid } or null. */
function abaqusCell(name) {
  const t = String(name).toUpperCase();
  let m = /^(?:DC|AC|Q|EMC)?C3D(\d+)/.exec(t);
  if (m) { const c = SOLID_N[+m[1]]; return c ? { type: c[0], n: c[1], total: +m[1] } : null; }
  if ((m = /^(?:SC|CSS)([68])/.exec(t))) return { type: m[1] === '8' ? 'hex' : 'wedge', n: +m[1], total: +m[1] };
  if ((m = /^(?:STRI|SFM3D|M3D|CPEG|CPS|CPE|CAX|DC2D|DCAX|R3D|DS|S)(\d)/.exec(t))) { const k = +m[1]; return k === 3 || k === 6 ? { type: 'tri', n: 3, total: k } : k === 4 || k === 8 || k === 9 ? { type: 'quad', n: 4, total: k } : null; }
  if ((m = /^T(\d)D(\d)/.exec(t))) return +m[2] === 3 ? { type: 'line', n: 3, total: 3, mid: true } : { type: 'line', n: 2, total: 2 };
  if ((m = /^(?:PIPE|ELBOW|B)(\d)(\d)/.exec(t))) return m[2] === '2' ? { type: 'line', n: 3, total: 3, mid: true } : { type: 'line', n: 2, total: 2 };
  if (/^FRAME[23]D/.test(t)) return { type: 'line', n: 2, total: 2 };
  return null;
}
async function readAbaqus(ctx) {
  const ls = splitLines(await ctx.text()), coords = [], store = cellStore(), warnings = [], skippedTypes = Object.create(null), notApplied = new Set();
  let idmap = new Map(), mode = null, cell = null, etype = '', pending = null, parts = 0, instances = 0, includes = 0, bad = 0, sawKeyword = false;
  const emit = (f) => {
    if (!cell) { skippedTypes[etype] = (skippedTypes[etype] || 0) + 1; return; }
    const ids = f.slice(1, 1 + cell.n).map((q) => idmap.get(q));
    if (ids.length < cell.n || ids.some((q) => q === undefined)) { bad++; return; }
    if (cell.mid) { store.add('line', [ids[0], ids[1]]); store.add('line', [ids[1], ids[2]]); } else store.add(cell.type, ids);
  };
  for (const line of ls) {
    if (line.startsWith('**')) continue;
    const t = line.trim();
    if (!t) continue;
    if (t[0] === '*') {
      const kw = t.slice(1).split(',')[0].trim().toUpperCase().replace(/\s+/g, ' ');
      sawKeyword = true; mode = null; pending = null;
      if (kw === 'NODE') { mode = 'node'; if (/SYSTEM\s*=\s*[CS]/i.test(t)) notApplied.add('*NODE, SYSTEM'); }
      else if (kw === 'ELEMENT') { const m = /TYPE\s*=\s*([^,\s]+)/i.exec(t); etype = m ? m[1].toUpperCase() : '?'; cell = abaqusCell(etype); mode = 'elem'; }
      else if (kw === 'PART') { parts++; idmap = new Map(); }
      else if (kw === 'INSTANCE') instances++;
      else if (kw === 'INCLUDE') includes++;
      else if (['NGEN', 'NFILL', 'NCOPY', 'NMAP', 'ELGEN', 'ELCOPY', 'SYSTEM', 'TRANSFORM'].includes(kw)) notApplied.add('*' + kw);
      continue;
    }
    if (mode === 'node') {
      const f = t.split(','), id = parseInt(f[0], 10), x = parseFloat(f[1]);
      if (!Number.isInteger(id) || !Number.isFinite(x)) { bad++; continue; }
      if (coords.length > 9 * L.cells) fail('The Abaqus deck holds too many nodes.');
      idmap.set(id, coords.length / 3); coords.push(x, parseFloat(f[2]) || 0, parseFloat(f[3]) || 0);
    } else if (mode === 'elem') {
      const f = t.split(',').map((s) => s.trim()).filter((s) => s !== '').map((s) => parseInt(s, 10));
      pending = pending ? pending.concat(f) : f;
      if (cell ? pending.length - 1 < cell.total && /,$/.test(t) : /,$/.test(t) && pending.length < 64) continue;
      emit(pending); pending = null;
    }
  }
  if (pending) emit(pending);
  if (!coords.length) fail(sawKeyword ? 'No *NODE data were found in the Abaqus input deck.' : 'Not an Abaqus input deck (no keyword lines).');
  const sk = Object.entries(skippedTypes);
  if (sk.length) warnings.push(`Elements of unsupported type were ignored: ${sk.slice(0, 8).map(([k, n]) => `${n} × ${k}`).join(', ')}.`);
  if (bad) warnings.push(`${bad} malformed node or element lines (or elements referring to missing nodes) were skipped.`);
  if (includes) warnings.push(`${includes} *INCLUDE files are not followed; only the data in this file are read.`);
  if (notApplied.size) warnings.push(`Not applied: ${[...notApplied].join(', ')} (generated nodes / elements and coordinate transforms are missing).`);
  if (instances > 1 || (parts > 1 && instances)) warnings.push(`${parts} parts / ${instances} instances: instance translations and rotations of the assembly are not applied, parts are shown in their own coordinates.`);
  return cellMesh(Float64Array.from(coords), store, { warnings, stats: parts ? { parts } : {} });
}

const CDB_FACE = new Set([13, 25, 28, 41, 42, 43, 55, 57, 63, 75, 77, 78, 82, 83, 93, 131, 132, 157, 163, 181, 182, 183, 208, 209, 223, 230, 233, 281]);
const CDB_SOLID = new Set([5, 45, 62, 64, 65, 69, 70, 87, 90, 92, 95, 96, 97, 98, 117, 122, 123, 164, 168, 185, 186, 187, 190, 226, 227, 231, 232, 236, 237, 278, 279, 285, 291]);
const CDB_LINE = new Set([1, 3, 4, 8, 10, 11, 12, 16, 17, 18, 20, 23, 24, 33, 44, 59, 180, 188, 189, 288, 289, 290]);
async function readCDB(ctx) {
  const ls = splitLines(await ctx.text()), et = new Map(), idmap = new Map(), coords = [], raw = [], warnings = [];
  let i = 0, unblocked = 0, oldBlocks = 0;
  const ints = (r, w) => { const out = []; for (let o = 0; o + 1 <= r.length; o += w) { const s = r.slice(o, o + w).trim(); if (s === '') break; out.push(parseInt(s, 10)); } return out; };
  while (i < ls.length) {
    const u = ls[i].trim().toUpperCase();
    if (/^ET\s*,/.test(u)) { const f = u.split(','), m = /(\d+)\s*$/.exec(f[2] || ''); if (m) et.set(parseInt(f[1], 10), +m[1]); i++; }
    else if (u.startsWith('NBLOCK')) {
      const fm = /\(\s*(\d+)i(\d+)\s*,\s*(\d+)e(\d+)/i.exec(ls[i + 1] || '');
      if (!fm) fail('The NBLOCK format line of the CDB file is missing or unusual.');
      const ni = +fm[1], wi = +fm[2], wr = +fm[4], o = ni * wi;
      for (i += 2; i < ls.length; i++) {
        const r = ls[i];
        if (/^\s*-1\s*$/.test(r) || /^[A-Za-z*/!]/.test(r.trim())) break;
        const id = parseInt(r.slice(0, wi), 10);
        if (!(id > 0)) continue;
        if (coords.length > 9 * L.cells) fail('The CDB file holds too many nodes.');
        idmap.set(id, coords.length / 3); coords.push(parseFloat(r.slice(o, o + wr)) || 0, parseFloat(r.slice(o + wr, o + 2 * wr)) || 0, parseFloat(r.slice(o + 2 * wr, o + 3 * wr)) || 0);
      }
    } else if (u.startsWith('EBLOCK')) {
      const solid = /,\s*SOLID/.test(u), fm = /\(\s*(\d+)i(\d+)/i.exec(ls[i + 1] || '');
      if (!fm) fail('The EBLOCK format line of the CDB file is missing or unusual.');
      const w = +fm[2];
      if (!solid) oldBlocks++;
      for (i += 2; i < ls.length; i++) {
        const r = ls[i];
        if (/^\s*-1\s*$/.test(r) || /^[A-Za-z*/!]/.test(r.trim())) break;
        if (!solid) continue;
        const f = ints(r, w);
        if (f.length < 12) continue;
        const nn = f[8], nodes = f.slice(11);
        if (!(nn > 0 && nn <= 64)) continue;
        while (nodes.length < nn && i + 1 < ls.length && !/^\s*-1\s*$/.test(ls[i + 1])) nodes.push(...ints(ls[++i], w));
        nodes.length = Math.min(nodes.length, nn);
        raw.push([f[1], nodes]);
        if (raw.length > L.cells) fail(`Mesh has too many cells (limit ${L.cells / 1e6} million).`);
      }
    } else { if (/^(N|EN|E)\s*,\s*\d/.test(u)) unblocked++; i++; }
  }
  if (!coords.length) fail(unblocked ? 'The CDB file uses unblocked N / EN commands. ' + fmtNamed('ANSYS CDB archive').convert : 'No NBLOCK node data were found: not an ANSYS CDB archive.');
  const store = cellStore(), xyz = Float64Array.from(coords);
  let skipped = 0, guessed = 0, missing = 0;
  for (const [type, nodes] of raw) {
    const en = et.get(type), nn = nodes.length, ids = nodes.map((q) => idmap.get(q)), ok = (k) => ids.slice(0, k).every((q) => q !== undefined);
    let kind = CDB_SOLID.has(en) ? 'solid' : CDB_FACE.has(en) ? 'face' : CDB_LINE.has(en) ? 'line' : null;
    if (!kind) {
      guessed++;
      if (nn <= 3) kind = 'line';
      else if (nn === 4 && ok(4)) { const P = ids.map((q) => [xyz[3 * q], xyz[3 * q + 1], xyz[3 * q + 2]]), vol = Math.abs(dot(cross(sub(P[1], P[0]), sub(P[2], P[0])), sub(P[3], P[0]))), sc = vlen(sub(P[1], P[0])) * vlen(sub(P[2], P[0])) * vlen(sub(P[3], P[0])); kind = vol > 1e-6 * sc ? 'solid' : 'face'; }
      else kind = 'solid';
    }
    if (kind === 'line') { if (nn >= 2 && ok(2)) store.add('line', [ids[0], ids[1]]); else missing++; }
    else if (kind === 'face') {
      const k = nn === 6 || nn === 3 ? 3 : 4;
      if (nn < 3 || !ok(k)) { missing++; continue; }
      if (k === 3 || ids[2] === ids[3]) store.add('tri', ids.slice(0, 3)); else store.add('quad', ids.slice(0, 4));
    } else if (nn === 4 || nn === 10) { if (ok(4)) store.add('tet', ids.slice(0, 4)); else missing++; }
    else if (nn >= 8) { const c = ok(8) ? degenerateSolid(ids.slice(0, 8)) : null; if (c) store.add(c[0], c[1]); else if (ok(8)) skipped++; else missing++; }
    else if (nn === 6 && ok(6)) store.add('wedge', ids); else if (nn === 5 && ok(5)) store.add('pyramid', ids); else skipped++;
  }
  if (oldBlocks) warnings.push(`${oldBlocks} EBLOCK sections without the SOLID key use a different layout and were not read.`);
  if (unblocked) warnings.push(`${unblocked} unblocked N / EN / E commands were ignored.`);
  if (guessed) warnings.push(`${guessed} elements have no recognised ET type; they were classified by node count and shape.`);
  if (skipped) warnings.push(`${skipped} elements of unsupported shape were ignored.`);
  if (missing) warnings.push(`${missing} elements refer to nodes that are not in NBLOCK and were skipped.`);
  return cellMesh(xyz, store, { warnings });
}

async function readDyna(ctx) {
  const ls = splitLines(await ctx.text()), idmap = new Map(), coords = [], raw = [], warnings = [], other = Object.create(null);
  let mode = null, longAll = false, long = false, twoCard = null, includes = 0, transforms = 0, sawKeyword = false, optCards = 0;
  const cut = (l, widths) => { if (l.includes(',')) return l.split(',').map((s) => s.trim()); const out = []; let o = 0; for (const w of widths) { if (o >= l.length) break; out.push(l.slice(o, o + w).trim()); o += w; } return out; };
  const intsOf8 = (l) => { const w = long ? 20 : 8, f = l.includes(',') ? l.split(',').map((s) => s.trim()) : cut(l, new Array(Math.ceil(l.length / w)).fill(w)); while (f.length && f[f.length - 1] === '') f.pop(); return f.map((s) => (s === '' ? 0 : parseInt(s, 10))); };
  for (const l of ls) {
    if (l[0] === '$' || !l.trim()) continue;
    if (l[0] === '*') {
      const kw = l.trim().toUpperCase(), name = kw.split(/\s+/)[0];
      sawKeyword = true; twoCard = null;
      if (name === '*KEYWORD') { longAll = /LONG\s*=\s*Y/.test(kw); mode = null; continue; }
      long = /\s\+\s*$/.test(kw) ? true : /\s-\s*$/.test(kw) ? false : longAll;
      if (name === '*END') break;
      if (name === '*NODE' || name === '*NODE_MERGE') mode = 'node';
      else if (name.startsWith('*ELEMENT_SHELL')) mode = 'shell';
      else if (name.startsWith('*ELEMENT_SOLID')) mode = 'solid';
      else if (name.startsWith('*ELEMENT_TSHELL')) mode = 'tshell';
      else if (name.startsWith('*ELEMENT_BEAM')) mode = 'beam';
      else { mode = null; if (name.startsWith('*ELEMENT_')) other[name] = 1; if (name.startsWith('*INCLUDE')) includes++; if (name === '*NODE_TRANSFORM' || name.startsWith('*DEFINE_TRANSFORMATION')) transforms++; }
      continue;
    }
    if (!mode) continue;
    if (mode === 'node') {
      const f = cut(l, long ? [20, 20, 20, 20] : [8, 16, 16, 16]), id = parseInt(f[0], 10), x = parseFloat(f[1]);
      if (!Number.isInteger(id) || !Number.isFinite(x)) continue;
      if (coords.length > 9 * L.cells) fail('The LS-DYNA deck holds too many nodes.');
      idmap.set(id, coords.length / 3); coords.push(x, parseFloat(f[2]) || 0, parseFloat(f[3]) || 0);
      continue;
    }
    if (/[.eE]/.test(l) && !twoCard) { optCards++; continue; }      // thickness / beta / ortho option cards hold reals
    const f = intsOf8(l);
    if (f.some((v) => !Number.isInteger(v))) continue;
    if (mode === 'solid') {
      if (twoCard) { raw.push(['solid', f.slice(0, 10)]); twoCard = null; }
      else if (f.length <= 2) twoCard = f;
      else raw.push(['solid', f.slice(2, 12)]);
    } else if (f.length >= 4) raw.push([mode, f.slice(2, 10)]);
    if (raw.length > L.cells) fail(`Mesh has too many cells (limit ${L.cells / 1e6} million).`);
  }
  if (!coords.length) fail(sawKeyword ? 'No *NODE data were found in the LS-DYNA keyword deck.' : 'Not an LS-DYNA keyword deck (no keyword lines).');
  const store = cellStore();
  let skipped = 0, missing = 0;
  for (const [kind, nodes] of raw) {
    const ids = nodes.map((q) => (q > 0 ? idmap.get(q) : -1));
    if (ids.some((q) => q === undefined)) { missing++; continue; }
    if (kind === 'beam') { if (ids[0] >= 0 && ids[1] >= 0) store.add('line', [ids[0], ids[1]]); else skipped++; }
    else if (kind === 'shell') { if (ids.length < 3 || ids[0] < 0 || ids[1] < 0 || ids[2] < 0) skipped++; else if (ids[3] === undefined || ids[3] < 0 || ids[3] === ids[2]) store.add('tri', ids.slice(0, 3)); else store.add('quad', ids.slice(0, 4)); }
    else {
      if (ids.length >= 10 && ids[8] >= 0 && ids[9] >= 0) { store.add('tet', ids.slice(0, 4)); continue; }
      const v = ids.slice(0, 8).filter((q) => q >= 0);
      if (v.length === 4) { store.add('tet', v); continue; }
      while (v.length && v.length < 8) v.push(v[v.length - 1]);
      const c = v.length === 8 ? degenerateSolid(v) : null;
      if (c) store.add(c[0], c[1]); else skipped++;
    }
  }
  const ot = Object.keys(other);
  if (ot.length) warnings.push(`Element keywords that are not read: ${ot.slice(0, 8).join(', ')}.`);
  if (optCards) warnings.push(`${optCards} element option cards (thickness, beta, orientation) were skipped.`);
  if (skipped) warnings.push(`${skipped} elements of unsupported shape were ignored.`);
  if (missing) warnings.push(`${missing} elements refer to nodes that are not defined in this file and were skipped.`);
  if (includes) warnings.push(`${includes} *INCLUDE keywords are not followed; only the data in this file are read.`);
  if (transforms) warnings.push('*NODE_TRANSFORM / *DEFINE_TRANSFORMATION are not applied.');
  return cellMesh(Float64Array.from(coords), store, { warnings });
}

// ---- SEG-Y --------------------------------------------------------------------------------------------------------
const EBCDIC = (() => { const t = new Array(256).fill(' '), put = (o, s) => { for (let k = 0; k < s.length; k++) t[o + k] = s[k]; }; put(0x81, 'abcdefghi'); put(0x91, 'jklmnopqr'); put(0xa2, 'stuvwxyz'); put(0xc1, 'ABCDEFGHI'); put(0xd1, 'JKLMNOPQR'); put(0xe2, 'STUVWXYZ'); put(0xf0, '0123456789'); for (const [o, c] of [[0x4b, '.'], [0x4c, '<'], [0x4d, '('], [0x4e, '+'], [0x50, '&'], [0x5b, '$'], [0x5c, '*'], [0x5d, ')'], [0x5e, ';'], [0x60, '-'], [0x61, '/'], [0x6b, ','], [0x6c, '%'], [0x6d, '_'], [0x6e, '>'], [0x6f, '?'], [0x7a, ':'], [0x7b, '#'], [0x7c, '@'], [0x7d, "'"], [0x7e, '='], [0x7f, '"']]) t[o] = c; return t; })();
const SEGY_BPS = { 1: 4, 2: 4, 3: 2, 5: 4, 6: 8, 8: 1, 9: 8, 10: 4, 11: 2, 12: 8, 16: 1 };
const SEGY_FMT = { 1: 'IBM float', 2: 'int32', 3: 'int16', 5: 'IEEE float32', 6: 'IEEE float64', 8: 'int8', 9: 'int64', 10: 'uint32', 11: 'uint16', 12: 'uint64', 16: 'uint8' };
async function readSEGY(ctx) {
  const u8 = await ctx.bytes(), dv = new DataView(u8.buffer, u8.byteOffset, u8.byteLength), warnings = [];
  if (u8.length < 3600 + 240) fail('Not a SEG-Y file (shorter than its 3600-byte file header plus one trace header).');
  let hi = 0;
  for (let k = 0; k < 3200; k++) if (u8[k] >= 0x80 || u8[k] === 0x40) hi++;
  const ebcdic = hi > 1600, head = ebcdic ? Array.from(u8.subarray(0, 3200), (b) => EBCDIC[b]).join('') : latin1.decode(u8.subarray(0, 3200)).replace(/[^\x20-\x7e]/g, ' ');
  let be = true;
  if (!SEGY_BPS[dv.getUint16(3224, false)]) { if (SEGY_BPS[dv.getUint16(3224, true)]) be = false; else fail('Not a SEG-Y file (the binary header holds no valid sample format code).'); }
  const u16 = (o) => dv.getUint16(o, !be), i16 = (o) => dv.getInt16(o, !be), i32 = (o) => dv.getInt32(o, !be);
  const fmt = u16(3224), bps = SEGY_BPS[fmt], dtUs = u16(3216), nsBin = u16(3220), ext = i16(3504), extraTH = u8.length >= 3508 ? i16(3506) : 0;
  if (extraTH > 0 && extraTH < 100 && (u8[3500] === 2 || u8[3501] === 2)) fail('SEG-Y revision 2 files with additional trace-header blocks are not read. ' + fmtNamed('SEG-Y seismic').convert);
  if (ext < 0 || ext > 100) warnings.push('The number of extended textual headers is not fixed; none are assumed.');
  let off = 3600 + (ext > 0 && ext <= 100 ? ext * 3200 : 0), nTr = 0, samples = 0, cut = false, capped = false;
  const pos = [], wd = [], tr = [];
  let nzC = 0, nzS = 0, arcsec = false, degUnits = false;
  while (off + 240 <= u8.length) {
    const ns = u16(off + 114) || nsBin, size = 240 + ns * bps;
    if (!ns) fail('The SEG-Y file declares no samples per trace.');
    if (off + size > u8.length) { cut = true; break; }
    if (nTr >= L.points) { capped = true; break; }
    const sc = i16(off + 70), k = sc > 0 ? sc : sc < 0 ? -1 / sc : 1, es = i16(off + 68), ke = es > 0 ? es : es < 0 ? -1 / es : 1, cu = i16(off + 88);
    const sx = i32(off + 72) * k, sy = i32(off + 76) * k, cx = i32(off + 180) * k, cy = i32(off + 184) * k;
    if (cx || cy) nzC++; if (sx || sy) nzS++;
    if (cu === 2) arcsec = true; else if (cu === 3) degUnits = true;
    pos.push(sx, sy, cx, cy); wd.push(i32(off + 60) * ke); tr.push(off, ns, i16(off + 108), u16(off + 116) || dtUs);
    samples += ns; nTr++; off += size;
  }
  if (!nTr) fail('The SEG-Y file holds no complete trace.');
  if (cut) warnings.push(`The file ends inside trace ${nTr + 1}; ${nTr} complete traces were read.`);
  if (capped) warnings.push(`Only the first ${nTr} traces were read.`);
  const useC = nzC >= nzS && nzC > 0, o = useC ? 2 : 0, sA = arcsec ? 1 / 3600 : 1, xyz = new Float64Array(nTr * 3), hasWD = wd.some((v) => v !== 0), located = nzC > 0 || nzS > 0;
  for (let t = 0; t < nTr; t++) { xyz[3 * t] = located ? pos[4 * t + o] * sA : t; xyz[3 * t + 1] = located ? pos[4 * t + o + 1] * sA : 0; xyz[3 * t + 2] = hasWD ? -Math.abs(wd[t]) : 0; }
  if (!located) warnings.push('The trace headers carry no source or CDP coordinates; traces are placed by their index along x.');
  else warnings.push(`Trace positions taken from the ${useC ? 'CDP (bytes 181-188)' : 'source (bytes 73-80)'} coordinates with the coordinate scalar applied${arcsec ? ', converted from arc-seconds to degrees' : ''}.`);
  const stats = { traces: nTr, samplesPerTrace: nsBin || tr[1], sampleInterval: dtUs / 1e6, sampleFormat: SEGY_FMT[fmt], textHeader: ebcdic ? 'EBCDIC' : 'ASCII', byteOrder: be ? 'big-endian' : 'little-endian' };
  let seabed = null;
  if ((fmt === 1 || fmt === 5) && samples <= 2e7 && dtUs > 0) {
    const ibm = (q) => { const b = dv.getUint32(q, !be), f = b & 0xffffff; return f ? (b >>> 31 ? -1 : 1) * (f / 16777216) * 16 ** (((b >>> 24) & 0x7f) - 64) : 0; }, val = fmt === 1 ? ibm : (q) => dv.getFloat32(q, !be);
    const twt = new Array(nTr), dist = [0];
    let picked = 0;
    for (let t = 0; t < nTr; t++) {
      const p = tr[4 * t] + 240, ns = tr[4 * t + 1];
      let peak = 0, first = -1;
      for (let s = 0; s < ns; s++) { const a = Math.abs(val(p + 4 * s)); if (a > peak && Number.isFinite(a)) peak = a; }
      if (peak > 0) for (let s = 0; s < ns; s++) if (Math.abs(val(p + 4 * s)) >= 0.5 * peak) { first = s; break; }
      twt[t] = first >= 0 ? Math.max(0, tr[4 * t + 2]) / 1000 + (first * tr[4 * t + 3]) / 1e6 : NaN;
      if (first >= 0) picked++;
      if (t) dist.push(dist[t - 1] + (arcsec || degUnits ? Math.hypot((xyz[3 * t] - xyz[3 * t - 3]) * 111320 * Math.cos((xyz[3 * t + 1] * Math.PI) / 180), (xyz[3 * t + 1] - xyz[3 * t - 2]) * 110540) : Math.hypot(xyz[3 * t] - xyz[3 * t - 3], xyz[3 * t + 1] - xyz[3 * t - 2])));
    }
    if (picked) {
      seabed = { distance: dist, twt, depth: twt.map((v) => (v * 1500) / 2) };
      stats.seabedPicked = picked;
      warnings.push(`Seabed ESTIMATE: the first sample reaching half of each trace's peak amplitude was taken as the seabed two-way time and converted with 1500 m/s (${picked} of ${nTr} traces). It is not an interpreted horizon: direct arrivals, noise and processing mutes shift it.`);
      if (!hasWD) { for (let t = 0; t < nTr; t++) xyz[3 * t + 2] = Number.isFinite(twt[t]) ? -seabed.depth[t] : 0; warnings.push('The trace headers hold no water depth, so the z of each trace position is that estimate.'); }
    }
  } else warnings.push(fmt === 1 || fmt === 5 ? 'The file is too large for the seabed estimate; only trace positions are returned.' : `Sample format ${SEGY_FMT[fmt]} is not decoded; only trace positions are returned.`);
  const g = pointsGeom(xyz, { warnings, stats });
  g.textHeader = head.match(/.{1,80}/g).map((s) => s.trimEnd()).filter(Boolean).slice(0, 40).join('\n');
  if (seabed && g.count === nTr) g.seabed = seabed;
  if (arcsec || degUnits) g.geographic = true;
  return g;
}

// ---- MATLAB Level 5 MAT-file ----------------------------------------------------------------------------------------
const MAT_DT = { 1: 'i1', 2: 'u1', 3: 'i2', 4: 'u2', 5: 'i4', 6: 'u4', 7: 'f4', 9: 'f8', 12: 'i8', 13: 'u8' };
const MAT_CLASS = { 1: 'cell', 2: 'struct', 3: 'object', 4: 'char', 5: 'sparse', 16: 'function handle', 17: 'opaque' };
async function readMAT(ctx) {
  const u8 = await ctx.bytes(), head = latin1.decode(u8.subarray(0, Math.min(116, u8.length))), how = fmtNamed('MATLAB MAT-file').convert, opts = ctx.opts, warnings = [];
  if (/^MATLAB 7\.3 MAT-file/.test(head) || hdf5Offset(u8) >= 0) return readHDF(ctx, 'matlab');
  if (u8.length < 136 || !/^MATLAB 5\.0 MAT-file/.test(head)) fail('Not a Level 5 MAT-file (version 4 files and other data are not read). ' + how);
  const little = u8[126] === 0x49 && u8[127] === 0x4d;
  if (!little && !(u8[126] === 0x4d && u8[127] === 0x49)) fail('The MAT-file header has no valid byte-order mark.');
  const vars = [], skipped = [];
  let budget = L.inflated, count = 0;
  const element = (b, p) => {
    if (p + 8 > b.length) return null;
    const dv = new DataView(b.buffer, b.byteOffset, b.byteLength), w = dv.getUint32(p, little);
    if (w >>> 16) { const n = w >>> 16; if (n > 4) fail('The MAT-file is corrupt (bad small data element).'); return { type: w & 0xffff, data: b.subarray(p + 4, p + 4 + n), next: p + 8 }; }
    const n = dv.getUint32(p + 4, little);
    if (p + 8 + n > b.length) fail('The MAT-file is truncated.');
    return { type: w, data: b.subarray(p + 8, p + 8 + n), next: p + 8 + n + (w === 15 ? 0 : (8 - (n % 8)) % 8) };
  };
  const matrix = (b) => {
    const fl = element(b, 0), dm = fl && element(b, fl.next), nm = dm && element(b, dm.next), pr = nm && element(b, nm.next);
    if (!fl || !dm || !nm || fl.data.length < 8) return;
    const flags = new DataView(fl.data.buffer, fl.data.byteOffset, 8).getUint32(0, little), cls = flags & 0xff, name = latin1.decode(nm.data).replace(/\0.*$/, '') || `var${vars.length + skipped.length + 1}`;
    const dims = Array.from(typed(dm.data, 0, Math.floor(dm.data.length / 4), 'i4', little));
    if (cls < 6 || cls > 15) { skipped.push(`${name} (${MAT_CLASS[cls] || 'class ' + cls})`); return; }
    const n = dims.reduce((a, k) => a * k, 1), dt = pr && MAT_DT[pr.type];
    if (!dt || !(n > 0) || dims.some((k) => !(k >= 0))) { skipped.push(`${name} (empty or unsupported storage)`); return; }
    if (n > L.voxels) fail(`Array "${name}" is too large.`);
    vars.push({ name, dims, count: n, data: typed(pr.data, 0, n, dt, little), complex: !!(flags & 0x0800), logical: !!(flags & 0x0200) });
  };
  for (let p = 128; p + 8 <= u8.length && count < 4096 && vars.length < 256; count++) {
    const el = element(u8, p);
    if (!el) break;
    p = el.next;
    if (el.type === 15) { const inf = await inflate(el.data, 'deflate', budget); budget -= inf.length; if (budget <= 0) fail('Decompressed data is too large.'); const inner = element(inf, 0); if (inner && inner.type === 14) matrix(inner.data); }
    else if (el.type === 14) matrix(el.data);
  }
  return matGeom(vars, skipped, opts, warnings);
}
/** MATLAB variables [{ name, dims, count, data (column-major), complex }] -> table, survey, grid, voxels (Level 5 and v7.3 alike). */
function matGeom(vars, skipped, opts, warnings) {
  if (!vars.length) fail(`The MAT-file holds no numeric arrays${skipped.length ? ` (skipped: ${skipped.slice(0, 6).join(', ')})` : ''}.`);
  if (skipped.length) warnings.push(`Variables that are not plain numeric arrays were skipped: ${skipped.slice(0, 8).join(', ')}.`);
  if (vars.some((v) => v.complex)) warnings.push('Imaginary parts of complex arrays are ignored.');
  const isVec = (v) => v.dims.length === 2 && (v.dims[0] === 1 || v.dims[1] === 1) && v.count > 1, byLen = new Map();
  for (const v of vars) if (isVec(v)) { const l = byLen.get(v.count); if (l) l.push(v); else byLen.set(v.count, [v]); }
  const groups = [...byLen.values()].filter((l) => l.length > 1).sort((a, b) => b.length * b[0].count - a.length * a[0].count), want = opts.variable !== undefined ? vars.find((v) => v.name === String(opts.variable)) : null, big = vars.slice().sort((a, b) => b.count - a.count)[0];
  if (opts.variable !== undefined && !want) fail(`The MAT-file has no numeric variable "${String(opts.variable).slice(0, 60)}" (it holds ${vars.slice(0, 12).map((v) => v.name).join(', ')}).`);
  const stats = { variables: vars.map((v) => v.name).slice(0, 40) };
  const table = (headers, col, n) => {
    const records = [];
    for (let r = 0; r < Math.min(n, L.rows); r++) records.push(Object.fromEntries(headers.map((h, c) => [h, col(r, c)])));
    if (n > L.rows) warnings.push(`Only the first ${L.rows} of ${n} rows are kept.`);
    const tc = tableColumns(headers);
    if (tc.role === 'survey') return surveyFromRecords(headers, records, tc, warnings);
    return { kind: 'table', headers, records, warnings, stats: { ...stats, rows: n, ...(tc.role ? { role: tc.role } : {}) } };
  };
  let g;
  if (!want && groups.length && groups[0].length * groups[0][0].count >= big.count) {
    const grp = groups[0].slice(0, 64);
    warnings.push(`${grp.length} vectors of length ${grp[0].count} were combined into one table (columns ${grp.map((v) => v.name).join(', ')}).`);
    g = table(grp.map((v) => v.name), (r, c) => grp[c].data[r], grp[0].count);
  } else {
    const v = want || big, d = v.dims.slice();
    while (d.length > 2 && d[d.length - 1] === 1) d.pop();
    if (vars.length > 1) warnings.push(`The file holds ${vars.length} numeric variables; "${v.name}" (${v.dims.join(' × ')}) is used${want ? '' : ' because it is the largest'}.`);
    stats.variable = v.name; stats.shape = v.dims;
    if (d.length === 3) g = npyGeom({ shape: d, fortran: true, data: v.data, dtype: 'MATLAB' }, opts);
    else if (d.length !== 2) fail(`MATLAB arrays with ${d.length} dimensions are not supported (use vectors, matrices or 3-D arrays).`);
    else {
      const [nr, nc] = d;
      if (nr === 1 || nc === 1) g = table([v.name], (r) => v.data[r], v.count);
      else if (nc <= 16 && nr > nc) g = table(Array.from({ length: nc }, (_, c) => `${v.name}_${c + 1}`), (r, c) => v.data[c * nr + r], nr);
      else if (nr <= 4 && nc > 4 * nr) { warnings.push(`Matrix "${v.name}" of ${nr} × ${nc} is read with its rows as table columns.`); g = table(Array.from({ length: nr }, (_, c) => `${v.name}_${c + 1}`), (r, c) => v.data[r * nr + c], nc); }
      else { const rm = new Float64Array(nr * nc); for (let r = 0; r < nr; r++) for (let c = 0; c < nc; c++) rm[r * nc + c] = v.data[c * nr + r]; warnings.push(`Matrix "${v.name}" of ${nr} × ${nc} read as a regular array (row 1 at y = 0, column 1 at x = 0).`); g = imageGeom(rm, nc, nr, 1, opts.spacing || [1, 1, 1], [0, 0, 0], opts, {}); }
    }
  }
  g.stats = { ...stats, ...(g.stats || {}) };
  g.warnings = [...new Set([...warnings, ...(g.warnings || [])])];
  return g;
}

// ---- HDF5 containers: CGNS, Salome MED, Exodus II, NetCDF-4, MATLAB v7.3 and plain HDF5 data sets --------------------
const CGNS_ELEM = { 2: ['node', 1], 3: ['line', 2], 4: ['line', 3], 5: ['tri', 3], 6: ['tri', 6], 7: ['quad', 4], 8: ['quad', 8], 9: ['quad', 9], 10: ['tet', 4], 11: ['tet', 10], 12: ['pyramid', 5], 13: ['pyramid', 14], 14: ['wedge', 6], 15: ['wedge', 15], 16: ['wedge', 18], 17: ['hex', 8], 18: ['hex', 20], 19: ['hex', 27], 21: ['pyramid', 13], 24: ['line', 4], 25: ['tri', 9], 26: ['tri', 10], 27: ['quad', 12], 28: ['quad', 16], 29: ['tet', 16], 30: ['tet', 20], 31: ['pyramid', 21], 32: ['pyramid', 29], 33: ['pyramid', 30], 34: ['wedge', 24], 35: ['wedge', 38], 36: ['wedge', 40], 37: ['hex', 32], 38: ['hex', 56], 39: ['hex', 64], 40: ['line', 5], 41: ['tri', 12], 42: ['tri', 15], 43: ['quad', 16], 44: ['quad', 25], 45: ['tet', 22], 46: ['tet', 34], 47: ['tet', 35], 48: ['pyramid', 29], 49: ['pyramid', 50], 50: ['pyramid', 55], 51: ['wedge', 33], 52: ['wedge', 66], 53: ['wedge', 75], 54: ['hex', 44], 55: ['hex', 98], 56: ['hex', 125] };
const CORNERS = { node: 1, line: 2, tri: 3, quad: 4, tet: 4, pyramid: 5, wedge: 6, hex: 8 };
const h5Text = (v) => (typeof v === 'string' ? v : ArrayBuffer.isView(v) ? latin1.decode(Uint8Array.from(v, (c) => c & 255)).replace(/\0[\s\S]*$/, '') : Array.isArray(v) && typeof v[0] === 'string' ? v[0] : '').trim();
/** CGNS (HDF5 flavour): every zone of every base; structured zones as blocks, unstructured Elements_t sections as cells. */
async function cgnsGeom(h5) {
  const warnings = [], stats = { bases: 0, zones: 0 };
  const target = (addr) => (h5Text(h5.attrs(addr).type) === 'LK' ? h5.resolve(' link', addr) : addr);
  const kidsOf = (addr) => { const out = []; for (const [name, l] of h5.kids(addr)) { if (name[0] === ' ' || l.addr === undefined) continue; const a = target(l.addr); if (a < 0) { stats.externalLinks = (stats.externalLinks || 0) + 1; continue; } out.push({ name, addr: a, label: h5Text(h5.attrs(a).label) }); } return out; };
  const data = async (addr) => { const d = h5.resolve(' data', addr); return d < 0 || h5.info(d).kind !== 'dataset' ? null : (await h5.read(d)).data; };
  const bases = kidsOf(h5.root).filter((n) => n.label === 'CGNSBase_t');
  if (!bases.length) fail('The CGNS file holds no CGNSBase_t node.');
  const blocks = [], store = cellStore(), xyzParts = [], skipped = {};
  let nodeBase = 0, sNodes = 0;
  for (const base of bases) {
    stats.bases++;
    const bk = kidsOf(base.addr), du = bk.find((n) => n.label === 'DimensionalUnits_t');
    if (du && !stats.units) { const t = await data(du.addr); if (t && t.length >= 64) stats.units = h5Text(t.subarray(32, 64)) || undefined; }
    for (const zone of bk.filter((n) => n.label === 'Zone_t')) {
      stats.zones++;
      const zk = kidsOf(zone.addr), zd = await data(zone.addr), zt = zk.find((n) => n.label === 'ZoneType_t'), ztype = zt ? h5Text(await data(zt.addr)) : 'Structured';
      const gc = zk.find((n) => n.label === 'GridCoordinates_t' && n.name === 'GridCoordinates') || zk.find((n) => n.label === 'GridCoordinates_t');
      if (!zd || !gc) { warnings.push(`Zone "${zone.name}" has no size or no grid coordinates and was skipped.`); continue; }
      const idim = Math.floor(zd.length / 3), size = Array.from(zd.subarray(0, idim)), nv = size.reduce((a, b) => a * b, 1), ck = kidsOf(gc.addr), comp = [];
      for (const nm of ['CoordinateX', 'CoordinateY', 'CoordinateZ']) { const c = ck.find((n) => n.name === nm); comp.push(c ? await data(c.addr) : null); }
      const cr = ck.find((n) => n.name === 'CoordinateR'), ct = ck.find((n) => n.name === 'CoordinateTheta');
      if (!comp[0] && cr && ct && !ck.some((n) => n.name === 'CoordinatePhi')) {
        const R = await data(cr.addr), T = await data(ct.addr);
        if (R && T && T.length >= R.length) { comp[0] = Float64Array.from(R, (q, i) => q * Math.cos(T[i])); comp[1] = Float64Array.from(R, (q, i) => q * Math.sin(T[i])); warnings.push(`Zone "${zone.name}": cylindrical coordinates converted to Cartesian.`); }
      }
      if (!comp[0] || !(nv > 0) || comp.some((c) => c && c.length < nv)) { warnings.push(`Zone "${zone.name}" has no usable CoordinateX / Y / Z arrays and was skipped.`); continue; }
      if (nv > 3 * L.cells) fail('The CGNS zone has too many nodes.');
      const xyz = new Float64Array(3 * nv);
      for (let d = 0; d < 3; d++) if (comp[d]) for (let i = 0; i < nv; i++) xyz[3 * i + d] = comp[d][i];
      if (/^Structured/i.test(ztype)) { blocks.push({ xyz, ni: size[0] || 1, nj: size[1] || 1, nk: size[2] || 1 }); sNodes += nv; continue; }
      xyzParts.push(xyz);
      // element sections in ElementRange order; NGON_n / NFACE_n are resolved once all sections are known
      const secs = [];
      for (const s of zk.filter((n) => n.label === 'Elements_t')) {
        const sd = await data(s.addr), sk = kidsOf(s.addr), by = (nm) => sk.find((n) => n.name === nm), rg = by('ElementRange') && (await data(by('ElementRange').addr)), cn = by('ElementConnectivity') && (await data(by('ElementConnectivity').addr));
        if (!sd || !rg || !cn) continue;
        secs.push({ type: sd[0], start: rg[0], end: rg[1], conn: cn, off: by('ElementStartOffset') ? await data(by('ElementStartOffset').addr) : null, pe: sd[0] === 22 && by('ParentElements') ? await data(by('ParentElements').addr) : null, name: s.name });
      }
      secs.sort((a, b) => a.start - b.start);
      const ngon = new Map(), nface = [], add = (kind, c, o) => { const n = CORNERS[kind], ids = new Array(n); for (let k = 0; k < n; k++) ids[k] = nodeBase + c[o + k] - 1; store.add(kind, ids); };
      for (const s of secs) {
        const ne = s.end - s.start + 1, c = s.conn;
        if (!(ne > 0)) continue;
        if (s.type === 20) {                                  // MIXED: each element is preceded by its type
          for (let e = 0, o = 0; e < ne && o < c.length; e++) { if (s.off) o = s.off[e]; const el = CGNS_ELEM[c[o]]; if (!el) { skipped['mixed type ' + c[o]] = (skipped['mixed type ' + c[o]] || 0) + 1; if (!s.off) break; continue; } if (o + 1 + el[1] > c.length) break; if (el[0] !== 'node') add(el[0], c, o + 1); o += 1 + el[1]; }
        } else if (s.type === 22 || s.type === 23) {         // polygon faces / polyhedral cells: offsets (CGNS 4) or inline counts (CGNS 3)
          for (let e = 0, o = 0; e < ne && o < c.length; e++) {
            let a, b;
            if (s.off) { a = s.off[e]; b = s.off[e + 1]; } else { a = o + 1; b = a + c[o]; o = b; }
            if (!(a >= 0 && b > a && b <= c.length)) break;
            if (s.type === 22) ngon.set(s.start + e, [c, a, b, s.pe && s.pe.length >= 2 * ne ? !s.pe[e] || !s.pe[ne + e] : null]); else nface.push([c, a, b]);
          }
        } else {
          const el = CGNS_ELEM[s.type];
          if (!el) { skipped[`element type ${s.type}`] = (skipped[`element type ${s.type}`] || 0) + ne; continue; }
          if (el[0] === 'node') continue;
          if (c.length < ne * el[1]) fail(`CGNS section "${s.name}" has fewer connectivity entries than its element range requires.`);
          for (let e = 0; e < ne; e++) add(el[0], c, e * el[1]);
        }
      }
      if (ngon.size) {
        const use = new Map();
        for (const [c, a, b] of nface) for (let k = a; k < b; k++) { const f = Math.abs(c[k]); use.set(f, (use.get(f) || 0) + 1); }
        let nb = 0;                                           // boundary = faces of one NFACE_n cell, or with one parent element
        for (const [id, [c, a, b, edge]] of ngon) { if (nface.length ? use.get(id) !== 1 : edge === false) continue; const ids = []; for (let k = a; k < b; k++) ids.push(nodeBase + c[k] - 1); store.add('poly', ids); nb++; }
        stats.polyhedra = (stats.polyhedra || 0) + nface.length; stats.polygonFaces = (stats.polygonFaces || 0) + ngon.size;
        if (nb < ngon.size) warnings.push(`Zone "${zone.name}": polyhedral cells (${nface.length ? 'NFACE_n' : 'ParentElements'}) reduced to their ${nb} boundary faces of ${ngon.size}.`);
      }
      nodeBase += nv;
    }
  }
  if (stats.externalLinks) warnings.push(`${stats.externalLinks} links to other CGNS files were not followed.`);
  if (stats.bases > 1) warnings.push(`The file holds ${stats.bases} bases (${bases.slice(0, 6).map((b) => b.name).join(', ')}); their zones are shown together.`);
  const sk = Object.entries(skipped);
  if (sk.length) warnings.push(`Unsupported element sections were skipped: ${sk.map(([k, v]) => `${v} × ${k}`).join(', ')}.`);
  const uNodes = nodeBase;
  if (blocks.length && (!store.t.length || sNodes >= uNodes)) {
    if (store.t.length) warnings.push('The file mixes structured and unstructured zones; the structured zones are shown.');
    return structuredGeom(blocks, { warnings, stats });
  }
  if (!xyzParts.length) fail('No zone with grid coordinates could be read from the CGNS file.');
  if (blocks.length) warnings.push('The file mixes structured and unstructured zones; the unstructured zones are shown.');
  const all = new Float64Array(3 * uNodes);
  let o = 0;
  for (const part of xyzParts) { all.set(part, o); o += part.length; }
  return cellMesh(all, store, { warnings, stats, ...(store.t.includes('poly') && !store.t.some((t) => CELL_FACES[t]) ? { flat: false } : {}) });
}

const MED_ELEM = { SE2: 'line', SE3: 'line', SE4: 'line', TR3: 'tri', TR6: 'tri', TR7: 'tri', QU4: 'quad', QU8: 'quad', QU9: 'quad', TE4: 'tet', T10: 'tet', PY5: 'pyramid', P13: 'pyramid', PE6: 'wedge', P15: 'wedge', P18: 'wedge', HE8: 'hex', H20: 'hex', H27: 'hex' };
/** Salome MED (2.x - 4.x): nodes and element connectivity of one unstructured mesh, quadratic cells reduced to corners. */
async function medGeom(h5, opts) {
  const maa = h5.resolve('/ENS_MAA');
  if (maa < 0) fail('The MED file holds no mesh (group ENS_MAA is missing); MED files with fields only are not read.');
  const names = [...h5.kids(maa)].filter(([, l]) => l.addr !== undefined), pick = names.find(([n]) => n === String(opts.mesh ?? '')) || names[0], warnings = [];
  if (!pick) fail('The MED file holds no mesh.');
  if (names.length > 1) warnings.push(`The file holds ${names.length} meshes (${names.slice(0, 8).map(([n]) => n).join(', ')}); "${pick[0]}" is used.`);
  let m = pick[1].addr;
  if (h5.resolve('NOE/COO', m) < 0) { const step = [...h5.kids(m)].find(([, l]) => l.addr !== undefined && h5.resolve('NOE/COO', l.addr) >= 0); if (!step) fail('The MED mesh has no node coordinates (NOE/COO); structured MED grids are not read.'); m = step[1].addr; }
  const coo = h5.resolve('NOE/COO', m), cd = (await h5.read(coo)).data, nn = +h5.attrs(coo).NBR || 0, dim = nn ? Math.round(cd.length / nn) : 0;
  if (!(nn > 0) || dim < 1 || dim > 3 || cd.length < nn * dim) fail('The MED node coordinates are inconsistent.');
  const xyz = new Float64Array(3 * nn), store = cellStore(), skipped = {};
  for (let d = 0; d < dim; d++) for (let i = 0; i < nn; i++) xyz[3 * i + d] = cd[d * nn + i];
  const mai = h5.resolve('MAI', m);
  for (const [type, l] of mai < 0 ? [] : h5.kids(mai)) {
    const nod = l.addr === undefined ? -1 : h5.resolve('NOD', l.addr), kind = MED_ELEM[type];
    if (nod < 0) continue;
    const c = (await h5.read(nod)).data;
    if (type === 'POG') {                                    // polygons: INN holds the 1-based start of each polygon in NOD
      const inn = h5.resolve('INN', l.addr), ix = inn < 0 ? null : (await h5.read(inn)).data;
      for (let e = 0; ix && e + 1 < ix.length; e++) { const ids = []; for (let k = ix[e] - 1; k < ix[e + 1] - 1 && k < c.length; k++) ids.push(c[k] - 1); if (ids.length > 2) store.add('poly', ids); }
      continue;
    }
    const ne = +h5.attrs(nod).NBR || 0, npe = ne ? Math.round(c.length / ne) : 0;
    if (!kind || !(ne > 0) || npe < CORNERS[kind]) { if (type !== 'PO1') skipped[type] = (skipped[type] || 0) + ne; continue; }
    for (let e = 0; e < ne; e++) { const ids = new Array(CORNERS[kind]); for (let k = 0; k < ids.length; k++) ids[k] = c[k * ne + e] - 1; store.add(kind, ids); }
  }
  const sk = Object.entries(skipped);
  if (sk.length) warnings.push(`Unsupported MED element types were skipped: ${sk.map(([k, v]) => `${v} × ${k}`).join(', ')}.`);
  return cellMesh(xyz, store, { warnings, stats: { mesh: pick[0], spaceDimension: dim } });
}

/** Exodus II over either NetCDF flavour: coordx / coordy / coordz (or coord) and the connect<N> element blocks. */
async function exodusGeom(nc) {
  const by = (n) => nc.vars.find((v) => v.name === n), cx = by('coordx'), co = by('coord'), warnings = [], skipped = {}, blocks = {};
  let xyz, nn;
  if (cx) {
    const c = [await cx.all(), by('coordy') ? await by('coordy').all() : null, by('coordz') ? await by('coordz').all() : null];
    nn = c[0].length; xyz = new Float64Array(3 * nn);
    for (let d = 0; d < 3; d++) if (c[d] && c[d].length >= nn) for (let i = 0; i < nn; i++) xyz[3 * i + d] = c[d][i];
  } else if (co && co.shape.length === 2 && co.shape[0] <= 3) {
    const a = await co.all(), nd = co.shape[0];
    nn = co.shape[1]; xyz = new Float64Array(3 * nn);
    for (let d = 0; d < nd; d++) for (let i = 0; i < nn; i++) xyz[3 * i + d] = a[d * nn + i];
  } else fail('Not an Exodus II mesh: the coordinate variables (coordx, coordy, coordz or coord) are missing.');
  const store = cellStore(), txt = (v) => (typeof v === 'string' ? v : Array.isArray(v) && typeof v[0] === 'string' ? v[0] : '').trim().toUpperCase();
  for (const v of nc.vars.filter((q) => /^connect\d+$/.test(q.name)).sort((a, b) => +a.name.slice(7) - +b.name.slice(7))) {
    if (v.shape.length !== 2 || !v.numeric) continue;
    const [ne, npe] = v.shape, et = txt(v.at.elem_type);
    const kind = /^HEX/.test(et) ? 'hex' : /^TET/.test(et) ? 'tet' : /^(WEDGE|PENTA|PRISM)/.test(et) ? 'wedge' : /^PYR/.test(et) ? 'pyramid' : /^TRI/.test(et) ? 'tri' : /^(QUAD|SHELL|RECT)/.test(et) ? (npe === 3 || npe === 6 ? 'tri' : 'quad') : /^(BAR|BEAM|TRUSS|EDGE|ROD|LINE)/.test(et) ? 'line' : null;
    if (!kind || npe < CORNERS[kind]) { skipped[et || 'untyped'] = (skipped[et || 'untyped'] || 0) + ne; continue; }
    const c = await v.all(), n = CORNERS[kind];
    blocks[et] = (blocks[et] || 0) + ne;
    for (let e = 0; e < ne; e++) { const ids = new Array(n); for (let k = 0; k < n; k++) ids[k] = c[e * npe + k] - 1; store.add(kind, ids); }
  }
  const sk = Object.entries(skipped);
  if (sk.length) warnings.push(`Element blocks of unsupported type were skipped: ${sk.map(([k, v]) => `${v} × ${k}`).join(', ')}.`);
  return cellMesh(xyz, store, { warnings, stats: { elementBlocks: blocks, title: typeof nc.at.title === 'string' ? nc.at.title.slice(0, 80) : undefined, container: nc.flavour } });
}

const MAT73_NUM = new Set(['double', 'single', 'int8', 'uint8', 'int16', 'uint16', 'int32', 'uint32', 'int64', 'uint64', 'logical']);
/** MATLAB v7.3: root data sets with a numeric MATLAB_class -> the variable records of the Level 5 reader. */
async function mat73Vars(h5) {
  const vars = [], skipped = [];
  for (const [name, l] of h5.kids(h5.root)) {
    if (l.addr === undefined || name[0] === '#') continue;
    const inf = h5.info(l.addr), at = h5.attrs(l.addr), cls = h5Text(at.MATLAB_class) || (inf.kind === 'group' ? 'struct' : 'unknown');
    const cplx = inf.kind === 'dataset' && inf.type && inf.type.cls === 6 && inf.type.members.length === 2 && inf.type.members.every((m) => isNumeric(m.type));
    if (inf.kind !== 'dataset' || !MAT73_NUM.has(cls) || !(isNumeric(inf.type) || cplx) || at.MATLAB_empty || at.MATLAB_sparse !== undefined) { skipped.push(`${name} (${at.MATLAB_sparse !== undefined ? 'sparse' : at.MATLAB_empty ? 'empty' : cls})`); continue; }
    if (inf.count > L.voxels) fail(`Array "${name}" is too large.`);
    if (vars.length >= 256) break;
    const r = await h5.read(l.addr), dims = inf.shape.slice().reverse();
    while (dims.length < 2) dims.push(1);
    vars.push({ name, dims, count: inf.count, data: cplx ? r.data[inf.type.members[0].name] : r.data, complex: cplx, logical: cls === 'logical' });
  }
  return { vars, skipped };
}

/** Plain HDF5: the most geometry-like numeric data set (opts.dataset = path overrides) -> points, polyline, grid, voxels or table. */
async function h5Generic(h5, opts, contents) {
  const warnings = [], ds = contents.filter((c) => c.kind === 'dataset' && c.numeric && c.count > 1), short = (c) => c.path.replace(/^.*\//, '');
  if (!ds.length) fail(`The HDF5 file holds no numeric data set${contents.length ? ` (it holds ${contents.slice(0, 6).map((c) => c.path).join(', ')})` : ''}.`);
  const isXYZ = (c) => c.shape.length === 2 && (c.shape[1] === 2 || c.shape[1] === 3) && c.shape[0] > c.shape[1], load = async (c) => (await h5.read(c.addr)).data;
  const big = (list) => list.slice().sort((a, b) => b.count - a.count)[0];
  let c = opts.dataset !== undefined ? ds.find((q) => q.path === String(opts.dataset) || q.path === '/' + String(opts.dataset)) : null;
  if (opts.dataset !== undefined && !c) fail(`The HDF5 file has no numeric data set "${String(opts.dataset).slice(0, 80)}" (it holds ${ds.slice(0, 10).map((q) => q.path).join(', ')}).`);
  const vecs = ds.filter((q) => q.shape.length === 1), byLen = new Map();
  for (const v of vecs) { const k = v.path.replace(/[^/]*$/, '') + '|' + v.count; if (byLen.has(k)) byLen.get(k).push(v); else byLen.set(k, [v]); }
  const groups = [...byLen.values()].filter((l) => l.length > 1).sort((a, b) => b.length * b[0].count - a.length * a[0].count);
  if (!c) c = big(ds.filter(isXYZ).filter((q) => /coord|point|xyz|vert|node|pos|route|path|centre|center|line|track/i.test(q.path))) || big(ds.filter(isXYZ));
  const stats = { datasets: ds.length };
  const table = (headers, col, n) => {
    const records = [];
    for (let r = 0; r < Math.min(n, L.rows); r++) records.push(Object.fromEntries(headers.map((h, k) => [h, col(r, k)])));
    if (n > L.rows) warnings.push(`Only the first ${L.rows} of ${n} rows are kept.`);
    const tc = tableColumns(headers);
    if (tc.role === 'survey') return surveyFromRecords(headers, records, tc, warnings);
    return { kind: 'table', headers, records, warnings, stats: { ...stats, rows: n, ...(tc.role ? { role: tc.role } : {}) } };
  };
  if (!c && groups.length && groups[0].length * groups[0][0].count >= big(ds).count) {
    const grp = groups[0].slice(0, 64), cols = [];
    for (const v of grp) cols.push(await load(v));
    warnings.push(`${grp.length} vectors of length ${grp[0].count} were combined into one table (columns ${grp.map(short).join(', ')}).`);
    return table(grp.map(short), (r, k) => cols[k][r], grp[0].count);
  }
  if (!c) c = big(ds.filter((q) => q.shape.length === 2)) || big(ds.filter((q) => q.shape.length === 3)) || big(ds);
  if (c.count > L.voxels) fail(`Data set "${c.path}" is too large.`);
  if (ds.length > 1) warnings.push(`The file holds ${ds.length} numeric data sets; "${c.path}" (${c.shape.join(' × ')}) is used. Pass opts.dataset to choose another.`);
  stats.dataset = c.path; stats.shape = c.shape; stats.dtype = c.type;
  const d = await load(c), sh = c.shape.filter((n, k) => n > 1 || k >= c.shape.length - 2);
  let g;
  if (isXYZ(c)) {
    const n = c.shape[0], w = c.shape[1], xyz = new Float64Array(3 * n);
    for (let i = 0; i < n; i++) for (let k = 0; k < w; k++) xyz[3 * i + k] = d[i * w + k];
    if (pathLike(xyz) && !/point|cloud|scan|vert|node/i.test(c.path)) { const x = [], y = [], z = []; for (let i = 0; i < n; i++) { x.push(xyz[3 * i]); y.push(xyz[3 * i + 1]); z.push(xyz[3 * i + 2]); } g = polyGeom([{ x, y, ...(w === 3 ? { z } : {}), closed: false }], { warnings, stats }); }
    else g = pointsGeom(xyz, { warnings, stats });
  } else if (sh.length === 1) g = table([short(c)], (r) => d[r], c.count);
  else if (sh.length === 2) {
    const [nr, ncol] = sh;
    if (ncol <= 16 && nr > ncol) g = table(Array.from({ length: ncol }, (_, k) => `${short(c)}_${k + 1}`), (r, k) => d[r * ncol + k], nr);
    else {
      // 1-D data sets of matching length beside the array serve as its coordinates
      const dir = c.path.replace(/[^/]*$/, ''), sib = vecs.filter((v) => v.path.startsWith(dir) && !v.path.slice(dir.length).includes('/')), cx = sib.find((v) => v.count === ncol && /^(x|lon|longitude|easting|east)$/i.test(short(v))), cy = sib.find((v) => v.count === nr && /^(y|lat|latitude|northing|north)$/i.test(short(v)));
      if (cx && cy) { const X = await load(cx), Y = await load(cy); g = gridGeom(Float64Array.from(d), ncol, nr, (i) => X[i], (j) => Y[j], { nodata: (q) => q !== q, geographic: /^lon/i.test(short(cx)), warnings, stats }); }
      else { warnings.push(`Array "${c.path}" of ${nr} × ${ncol} read as a regular array (row 1 at y = 0, column 1 at x = 0).`); g = imageGeom(Float64Array.from(d), ncol, nr, 1, opts.spacing || [1, 1, 1], [0, 0, 0], opts, { warnings, stats }); }
    }
  } else if (sh.length === 3) g = npyGeom({ shape: sh, fortran: false, data: d, dtype: c.type }, opts);
  else fail(`HDF5 data sets with ${sh.length} dimensions are not supported (use vectors, matrices or 3-D arrays).`);
  g.stats = { ...stats, ...(g.stats || {}) };
  g.warnings = [...new Set([...warnings, ...(g.warnings || [])])];
  return g;
}

/** Any HDF5 container. hint: 'cgns' | 'med' | 'exodus' | 'netcdf' | 'matlab' forces that reading; otherwise the content decides. */
async function readHDF(ctx, hint) {
  const u8 = await ctx.bytes(), opts = ctx.opts;
  if (u8.length > 4 && u8[0] === 0x0e && u8[1] === 0x03 && u8[2] === 0x13 && u8[3] === 0x01) fail('This is an HDF4 file, which is not read (HDF5 is). Convert it with h4toh5 in.hdf out.h5, or export the data set as GeoTIFF / CSV.');
  if (hint === 'cgns' && /ADF Database Version/.test(latin1.decode(u8.subarray(0, 64)))) fail('This CGNS file uses the older ADF container, which is not read (the HDF5 container is). Convert it with the CGNS tools: cgnsconvert -h in.cgns out.cgns (or adf2hdf in.cgns out.cgns).');
  const h5 = openHDF5(u8, { maxBytes: L.inflated }), rootKids = h5.kids(h5.root), rat = h5.attrs(h5.root), contents = [];
  const has = (n) => rootKids.has(n);
  let dimScales = false;
  h5.walk((path, addr, inf) => {
    if (path === '/') return;
    if (contents.length < 400) contents.push({ path, kind: inf.kind, addr, ...(inf.kind === 'dataset' ? { shape: inf.shape, type: typeName(inf.type), count: inf.count, numeric: isNumeric(inf.type) && inf.type.cls !== 10 } : {}) });
    if (inf.kind === 'dataset' && !dimScales && contents.length < 60) { const a = h5.attrs(addr); if (a.DIMENSION_LIST !== undefined || a.CLASS === 'DIMENSION_SCALE') dimScales = true; }
  }, { maxDepth: 12, maxNodes: 50000 });
  const isCgns = has('CGNSLibraryVersion') || [...rootKids].some(([n, l]) => n[0] !== ' ' && l.addr !== undefined && h5Text(h5.attrs(l.addr).label) === 'CGNSBase_t');
  const isMat = h5.userBlock === 512 && /^MATLAB 7\.3 MAT-file/.test(latin1.decode(u8.subarray(0, 24)));
  const kind = hint && hint !== 'netcdf' ? hint : isCgns ? 'cgns' : has('ENS_MAA') ? 'med' : isMat ? 'matlab' : has('coordx') || (has('coord') && has('connect1')) ? 'exodus' : rat._NCProperties !== undefined || dimScales || hint === 'netcdf' ? 'netcdf' : 'hdf5';
  let g;
  if (kind === 'cgns') { if (!isCgns) fail('The file is an HDF5 container but not a CGNS tree (no CGNSBase_t node).'); g = Object.assign(await cgnsGeom(h5), { format: 'CGNS' }); }
  else if (kind === 'med') g = Object.assign(await medGeom(h5, opts), { format: 'Salome MED' });
  else if (kind === 'exodus') g = Object.assign(await exodusGeom(cdfHDF5(h5)), { format: 'Exodus II' });
  else if (kind === 'matlab') { const { vars, skipped } = await mat73Vars(h5); g = Object.assign(matGeom(vars, skipped, opts, []), { format: 'MATLAB MAT-file' }); g.stats.container = 'HDF5 (v7.3)'; }
  else {
    let err = null;
    if (kind === 'netcdf') try { g = Object.assign(await cdfGrid(cdfHDF5(h5)), { format: 'NetCDF-4 (HDF5)' }); } catch (e) { if (!e || !e.user) throw e; err = e; }
    if (!g) try { g = Object.assign(await h5Generic(h5, opts, contents), { format: kind === 'netcdf' ? 'NetCDF-4 (HDF5)' : 'HDF5' }); } catch (e) { if (err && e && e.user) throw err; throw e; }
  }
  g.contents = contents.slice(0, 200).map(({ path, kind: k, shape, type }) => ({ path, kind: k, ...(shape ? { shape, type } : {}) }));
  if (contents.length >= 400) g.warnings.push('The file holds more objects than are listed in the contents (first 200 shown).');
  return g;
}

// ---- DEM / DTM / DSM text elevation models ---------------------------------------------------------------------------
/** USGS ASCII DEM (1024-byte record A, then one profile record per column) -> grid, or null when the header does not fit. */
function usgsDEM(text, warnings) {
  if (text.length < 1100) return null;
  const A = text.slice(0, 1024), num = (s) => +s.trim().replace(/[dD]/, 'e'), int = (a, b) => parseInt(A.slice(a, b), 10);
  const guni = int(528, 534), euni = int(534, 540), ncols = int(858, 864), res = [num(A.slice(816, 828)), num(A.slice(828, 840)), num(A.slice(840, 852))];
  if (!(ncols > 0 && ncols <= 20000) || !(res[0] > 0 && res[1] > 0) || ![0, 1, 2, 3].includes(guni) || ![1, 2].includes(euni)) return null;
  const re = /[-+]?(?:\d+\.?\d*|\.\d+)(?:[dDeE][-+]?\d+)?/g, zres = res[2] > 0 ? res[2] : 1, hs = guni === 3 ? 1 / 3600 : guni === 1 ? 0.3048 : guni === 0 ? 180 / Math.PI : 1, zs = euni === 1 ? 0.3048 : 1;
  re.lastIndex = 1024;
  const next = () => { const m = re.exec(text); return m ? +m[0].replace(/[dD]/, 'e') : NaN; }, profs = [];
  let total = 0;
  for (let p = 0; p < ncols; p++) {
    next(); const col = next(), m = next(), n = next(), x0 = next(), y0 = next(), datum = next();
    next(); next();
    if (!(m > 0 && m <= 1e5 && n === 1) || !Number.isFinite(col) || !Number.isFinite(datum)) break;
    total += m;
    if (total > L.voxels) fail('The DEM is too large.');
    const z = new Float64Array(m);
    for (let k = 0; k < m; k++) { const v = next(); z[k] = v === -32767 || v === -32768 || v !== v ? NaN : (datum + v * zres) * zs; }
    profs.push({ x: x0 * hs, y: y0 * hs, z });
  }
  if (!profs.length) return null;
  if (profs.length < ncols) warnings.push(`The DEM header declares ${ncols} profiles but ${profs.length} could be read; the file appears truncated.`);
  profs.sort((a, b) => a.x - b.x);
  const dy = res[1] * hs, y0 = Math.min(...profs.map((q) => q.y)), ny = Math.round(Math.max(...profs.map((q) => q.y + (q.z.length - 1) * dy - y0)) / dy) + 1, nx = profs.length;
  if (!(ny > 0) || nx * ny > L.voxels) fail('The DEM is too large.');
  const zf = new Float64Array(nx * ny).fill(NaN);
  profs.forEach((q, i) => { const j0 = Math.round((q.y - y0) / dy); for (let k = 0; k < q.z.length; k++) if (j0 + k >= 0 && j0 + k < ny) zf[(j0 + k) * nx + i] = q.z[k]; });
  if (euni === 1) warnings.push('Elevations converted from feet to metres.');
  if (guni === 1) warnings.push('Ground coordinates converted from feet to metres.');
  return gridGeom(zf, nx, ny, (i) => profs[i].x, (j) => y0 + j * dy, { nodata: (v) => v !== v, geographic: guni === 3 || guni === 0, warnings, stats: { cellsize: [res[0] * hs, dy], demFormat: 'USGS ASCII DEM' } });
}
async function readDEM(ctx) {
  const u8 = await ctx.bytes(), n = Math.min(u8.length, 2048), how = fmtNamed('DEM / DTM / DSM elevation model').convert;
  let odd = 0;
  for (let k = 0; k < n; k++) if ((u8[k] < 9 || (u8[k] > 13 && u8[k] < 32) || u8[k] > 126)) odd++;
  if (odd > 0.02 * n) fail('This elevation model is a binary file, which is not read. ' + how);
  const text = await ctx.text(), warnings = [];
  if (/^\s*ncols\s/i.test(text.slice(0, 200))) return Object.assign(readASCGrid(text), { format: 'ESRI ASCII grid' });
  const us = usgsDEM(text, warnings);
  if (us) return us;
  let g = null;
  try { g = await readDelimited(ctx); } catch (e) { if (!(e && e.user)) throw e; }
  if (!g || (g.kind !== 'grid' && g.kind !== 'points')) fail('The elevation model is neither an ESRI ASCII grid, a USGS ASCII DEM nor rows of x y z. ' + how);
  return Object.assign(g, { pathway: 'gis' });
}

async function readIFCZIP(ctx) {
  const u8 = await ctx.bytes();
  if (!isZip(u8)) return readIFC(ctx);
  const f = await unzipFiles(await ctx.buf(), (n) => /\.(ifc|ifcxml)$/i.test(n)), ifc = [...f.keys()].find((n) => /\.ifc$/i.test(n));
  if (!ifc) fail(f.size ? 'The IFCZIP archive holds ifcXML only, which is not read. ' + fmtNamed('IFCZIP').convert : 'No .ifc file was found inside the IFCZIP archive.');
  const text = utf8.decode(f.get(ifc));
  return readIFC({ ...ctx, text: async () => text });
}

// ---- Dispatch ------------------------------------------------------------------------------------------
const sniffHead = async (ctx, n = 4096) => latin1.decode((await ctx.bytes()).subarray(0, n));
const named = (name, fn) => async (ctx) => Object.assign(await fn(ctx), { format: name });
const READERS = {
  step: readSTEP, stp: readSTEP, p21: readSTEP, iges: readIGES, igs: readIGES,
  stl: async (ctx) => {
    let b = await ctx.buf();
    const u = await ctx.bytes(), need = u.length > 84 ? 84 + new DataView(b).getUint32(80, true) * 50 : 0;   // binary STL padded with a few trailing bytes
    if (need > 84 && need < u.length && u.length - need <= 1024 && !/^\s*solid[^\0]*facet/.test(latin1.decode(u.subarray(0, 1024)))) b = b.slice(0, need);
    return parseSTL(b);
  }, obj: async (ctx) => parseOBJ(await ctx.text()),
  ply: readPLY, off: readOFF, '3mf': read3MF, amf: readAMF, wrl: readVRML, vrml: readVRML, iv: readVRML, x3d: (ctx) => readX3D(ctx), dae: (ctx) => readDAE(ctx), gltf: (ctx) => readGLTF(ctx), glb: (ctx) => readGLTF(ctx),
  vtk: readVTK, vtp: (ctx) => readVTKXML(ctx), vtu: (ctx) => readVTKXML(ctx), vti: (ctx) => readVTKXML(ctx), vts: (ctx) => readVTKXML(ctx), vtr: (ctx) => readVTKXML(ctx), byu: readBYU, gts: readGTS, tri: readNodeFace, fac: readNodeFace, surf: readNodeFace,
  msh: async (ctx) => { const h = await sniffHead(ctx); if (/\$MeshFormat|\$NOD/.test(h)) return readGmsh(ctx); if (/^\s*\(/.test(h)) return named('Fluent mesh (ASCII)', readFluent)(ctx); return fail('The .msh file is neither a Gmsh mesh ($MeshFormat) nor a Fluent ASCII mesh. ' + fmtNamed('Gmsh MSH').convert); },
  cas: readFluent, su2: readSU2, unv: readUNV, bdf: readNastran, nas: readNastran, nastran: readNastran, tec: readTecplot, p3d: readPlot3D, x: readPlot3D, foam: readFoam,
  dat: async (ctx) => (/^\s*(TITLE|VARIABLES|ZONE)\b/im.test(await sniffHead(ctx)) ? readTecplot(ctx) : named('Delimited / column text', readDelimited)(ctx)),
  plt: async (ctx) => { const h = await sniffHead(ctx); if (h.startsWith('#!TDV')) fail('Binary Tecplot .plt files are not supported. ' + fmtNamed('Tecplot ASCII').convert); return /\b(ZONE|VARIABLES|TITLE)\b/i.test(h) ? readTecplot(ctx) : named('HPGL plot file', readHPGL)(ctx); },
  dxf: readDXF, svg: (ctx) => readSVG(ctx), hpgl: readHPGL, hpg: readHPGL, xy: named('x-y polyline', readDelimited),
  geojson: readJSON, json: readJSON, kml: (ctx) => readKML(ctx), kmz: (ctx) => readKML(ctx), gpx: (ctx) => readGPX(ctx), gml: (ctx) => readGML(ctx), shp: readSHP, mif: readMIF,
  asc: async (ctx) => (/^\s*ncols\s/i.test(await sniffHead(ctx, 200)) ? readASCGrid(await ctx.text()) : Object.assign(await readDelimited(ctx), { format: 'XYZ points / soundings', pathway: 'points' })),
  grd: readSurfer, tif: readTIFF, tiff: readTIFF, nc: readNetCDF, cdf: readNetCDF, nc4: readNetCDF, bil: readBandRaster, bip: readBandRaster, bsq: readBandRaster,
  las: readLAS, pts: (ctx) => readPoints(ctx), ptx: readPTX, xyzi: readDelimited, xyzrgb: readDelimited,
  xyz: async (ctx) => { const h = await sniffHead(ctx, 400), l1 = (splitLines(h).find((l) => l.trim()) || '').trim().split(/\s+/); if (l1.length <= 3 && l1.every((t) => /^\d+$/.test(t))) { const b = plot3dBlocks(numsOf(await ctx.text())); if (b) return Object.assign(structuredGeom(b, {}), { format: 'Plot3D grid (ASCII)' }); return readPoints(ctx, (await ctx.text()).replace(/^\s*\d+\s*(\r\n|\n|\r)/, '')); } return readDelimited(ctx); },
  png: readImage, jpg: readImage, jpeg: readImage, bmp: readImage, webp: readImage, gif: readImage,
  raw: readRAW, vol: readRAW, bin: readRAW, nrrd: readNRRD, nhdr: readNRRD, mha: readMHA, mhd: readMHA, nii: readNIfTI, 'nii.gz': readNIfTI, npy: readNPY, npz: readNPY, dcm: readDICOM, dicom: readDICOM,
  graphml: async (ctx) => graphmlGeom(parseXML(await ctx.text())), landxml: readLandXML, inp: readAbaqus, cdb: readCDB, k: readDyna, key: readDyna, dyn: readDyna,
  dev: readSurvey, wbt: readSurvey, survey: readSurvey, sgy: readSEGY, segy: readSEGY, mat: readMAT, dem: readDEM, dtm: readDEM, dsm: readDEM, ifczip: readIFCZIP,
  x_t: readXT, x_b: readXT, xmt_txt: readXT, xmt_bin: readXT, sat: readSAT, sab: readSAT, dwg: readDWG, dgn: readDGN, gpkg: readGPKG, e57: readE57, laz: readLAZ, cgns: (ctx) => readHDF(ctx, 'cgns'), med: (ctx) => readHDF(ctx, 'med'), rmed: (ctx) => readHDF(ctx, 'med'), e: readExodus, exo: readExodus, ex2: readExodus, exii: readExodus, g: readExodus,
  h5: (ctx) => readHDF(ctx), hdf5: (ctx) => readHDF(ctx), hdf: (ctx) => readHDF(ctx), he5: (ctx) => readHDF(ctx), h4: (ctx) => readHDF(ctx),
  pcf: readPCF, ifc: readIFC, aml: async (ctx) => amlGeom(parseXML(await ctx.text())), yaml: readYAML, yml: readYAML, xml: readXML, csv: readDelimited, tsv: readDelimited, txt: readDelimited,
};

/**
 * Read any registered geometry file into a Geometry object. `file` needs { name, size, arrayBuffer() } (text() optional).
 * opts: { dims, dtype, threshold, invert, spacing, companion: { name: file } } — companions serve multi-file formats.
 */
export async function importGeometry(file, opts = {}) {
  checkFile(file);
  opts = opts && typeof opts === 'object' ? opts : {};
  const name = String(file.name || ''), fmt = formatOf(name);
  let ext = extName(name);
  if (!fmt) {
    if (['shx', 'dbf', 'prj'].includes(ext)) fail('Select the .shp file of the shapefile; its .dbf and .prj can be supplied as companions.');
    if (ext === 'hdr') fail('Select the raster data file (.bil / .bip / .bsq) and supply this .hdr as its companion.');
    fail(`Unsupported geometry format "${ext ? '.' + ext : name}". See the format list for what can be read or how to convert it.`);
  }
  if (fmt.support === 'convert') {
    if (ext === 'jt' && typeof file.arrayBuffer === 'function') { const b = new Uint8Array(await file.arrayBuffer()); if (isJT(b)) parseJT(b); }   // says what the JT file holds
    fail(fmt.convert);
  }
  if (!ext) ext = 'foam';
  let buf, u8, txt;
  const comp = opts.companion && typeof opts.companion === 'object' ? (Array.isArray(opts.companion) ? opts.companion.map((f) => [f && f.name, f]) : Object.entries(opts.companion)) : [];
  const ctx = {
    file, name, ext, opts,
    buf: async () => { if (!buf) { const b = await file.arrayBuffer(); buf = b instanceof ArrayBuffer ? b : b.buffer.slice(b.byteOffset, b.byteOffset + b.byteLength); if (buf.byteLength > L.fileBytes) fail(`File is larger than ${L.fileBytes / 1e6} MB.`); } return buf; },
    bytes: async () => u8 || (u8 = new Uint8Array(await ctx.buf())),
    text: async () => { if (txt === undefined) txt = utf8.decode(await ctx.bytes()).replace(/^﻿/, ''); return txt; },
    companion: (pred) => { for (const [k, f] of comp) if (f && typeof f.arrayBuffer === 'function' && (pred(String(k || '')) || pred(String(f.name || '')))) { if (f.size > L.fileBytes) fail(`Companion file is larger than ${L.fileBytes / 1e6} MB.`); return f; } return null; },
  };
  const reader = READERS[ext];
  if (!reader) fail(`Unsupported geometry format ".${ext}".`);
  let g;
  try {
    if (!(await ctx.bytes()).length) fail('The file is empty.');
    g = await reader(ctx);
  } catch (e) {
    if (e && e.user) throw e;
    if (e instanceof Error && /Too many rows|too many triangles|Could not read any triangles|No faces found/i.test(e.message)) throw e;
    throw new Error(`Could not read "${name}" as ${fmt.name}: the file is truncated or malformed.`, { cause: e });
  }
  const used = fmtNamed(g.format) || fmt, { kind, format, pathway, warnings, stats, bbox, ...rest } = g;
  const out = { name, format: format || used.name, pathway: pathway || used.pathway, kind, warnings: warnings || [], stats: { ...(stats || {}) }, bbox: bbox || { min: [], max: [] }, ...rest };
  if (kind === 'mesh') out.stats.triangles = out.count;
  else if (kind === 'polylines') { out.stats.polylines = out.polylines.length; out.stats.vertices = out.polylines.reduce((s, p) => s + p.x.length, 0); }
  return out;
}

// ---- Analysis helpers ------------------------------------------------------------------------------------
/** Marching squares over a point-sampled field val(i, j) -> segments [x1, y1, x2, y2] of the iso-line at `level`. */
function contour(val, nx, ny, level, xOf, yOf) {
  const segs = [];
  for (let j = 0; j + 1 < ny; j++) for (let i = 0; i + 1 < nx; i++) {
    const v = [val(i, j), val(i + 1, j), val(i + 1, j + 1), val(i, j + 1)];
    if (v.some((q) => q !== q)) continue;
    const X = [xOf(i), xOf(i + 1), xOf(i + 1), xOf(i)], Y = [yOf(j), yOf(j), yOf(j + 1), yOf(j + 1)], cr = [];
    for (let e = 0; e < 4; e++) { const a = e, b = (e + 1) % 4; if ((v[a] > level) !== (v[b] > level)) { const t = (level - v[a]) / (v[b] - v[a]); cr.push([X[a] + t * (X[b] - X[a]), Y[a] + t * (Y[b] - Y[a]), e]); } }
    if (cr.length === 2) segs.push([cr[0][0], cr[0][1], cr[1][0], cr[1][1]]);
    else if (cr.length === 4) {
      const mid = (v[0] + v[1] + v[2] + v[3]) / 4 > level, pair = mid === v[0] > level ? [[0, 1], [2, 3]] : [[3, 0], [1, 2]];
      for (const [p, q] of pair) segs.push([cr[p][0], cr[p][1], cr[q][0], cr[q][1]]);
    }
  }
  return segs;
}
const asVoxels = (g) => {
  if (g && g.kind === 'voxels') return g.voxels;
  if (g && g.data && g.nx > 0 && g.ny > 0) return { nz: 1, spacing: [1, 1, 1], origin: [0, 0, 0], ...g };
  const m = Array.isArray(g) ? g : g && Array.isArray(g.mask) ? g.mask : null;
  if (m && m.length && m[0] && m[0].length !== undefined) { const ny = m.length, nx = m[0].length, data = new Uint8Array(nx * ny); for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) data[j * nx + i] = m[j][i] ? 1 : 0; return { nx, ny, nz: 1, data, spacing: [1, 1, 1], origin: [0, 0, 0] }; }
  throw new Error('A voxel geometry or a 2-D mask is required.');
};
/** One axis-aligned slice of a voxel volume as a 2-D solid sampler in physical coordinates. */
function voxelSlice(g, axis, value) {
  const v = asVoxels(g), d = [v.nx, v.ny, v.nz], sp = v.spacing || [1, 1, 1], or = v.origin || [0, 0, 0], ax = v.nz === 1 ? 2 : axis === 0 || axis === 1 ? axis : 2, [a, b] = [0, 1, 2].filter((k) => k !== ax);
  const k = value === undefined || value === null ? d[ax] >> 1 : Math.max(0, Math.min(d[ax] - 1, Math.floor((value - or[ax]) / sp[ax]))), st = [1, v.nx, v.nx * v.ny];
  return { w: d[a], h: d[b], dx: sp[a], dy: sp[b], x0: or[a], y0: or[b], index: k, get: (i, j) => (i >= 0 && j >= 0 && i < d[a] && j < d[b] ? v.data[k * st[ax] + i * st[a] + j * st[b]] : 0) };
}

/**
 * 2-D outline of a geometry as segments [x1, y1, x2, y2]: mesh cut by the plane {axis = value}, polylines as they are,
 * solid/pore boundary of a voxel slice, iso-line of a grid at z = value, or the plan view of a network.
 */
export function sectionOf(g, { axis = 2, value } = {}) {
  if (!g) return [];
  const fin = (s) => s.every(Number.isFinite);
  if (g.kind === 'mesh') return sliceMesh(g, axis, value).filter(fin);
  if (g.kind === 'polylines') return polylinesToSegments(g.polylines || []).filter(fin);
  if (g.kind === 'voxels') { const s = voxelSlice(g, axis, value); return contour((i, j) => s.get(i - 1, j - 1), s.w + 2, s.h + 2, 0.5, (i) => s.x0 + (i - 0.5) * s.dx, (j) => s.y0 + (j - 0.5) * s.dy); }
  if (g.kind === 'grid') { const { x, y, z } = g.grid, lv = value ?? 0.5 * (g.bbox.min[2] + g.bbox.max[2]); return contour((i, j) => z[j][i], x.length, y.length, lv, (i) => x[i], (j) => y[j]); }
  if (g.kind === 'network') { const at = new Map(g.network.nodes.map((n) => [n.id, n])); return g.network.edges.map((e) => { const a = at.get(e.from), b = at.get(e.to); return a && b ? [a.x, a.y, b.x, b.y] : [NaN]; }).filter((s) => s.length === 4 && s.every((q) => typeof q === 'number' && Number.isFinite(q))); }
  return [];
}

/**
 * Solid mask[ny][nx] (true = solid) of a geometry over the window [x0, x1] × [y0, y1].
 * opt: axis / value (section plane), fit (default true: scale the outline's bounding box uniformly into the window,
 * centred), margin (fraction of the window left free on each side), stretch (fill the window anisotropically), invert.
 */
export function maskOf(g, x0, x1, y0, y1, nx, ny, opt = {}) {
  const { axis = 2, value, fit = true, invert = false, margin = 0, stretch = false } = opt;
  if (!(nx > 0 && ny > 0) || !(x1 > x0) || !(y1 > y0)) throw new Error('maskOf needs a non-empty window and grid.');
  const vox = g && g.kind === 'voxels' ? voxelSlice(g, axis, value) : null, grd = g && g.kind === 'grid' ? g.grid : null;
  let segs = null, bb;
  if (vox) bb = [vox.x0, vox.y0, vox.x0 + vox.w * vox.dx, vox.y0 + vox.h * vox.dy];
  else if (grd) bb = [grd.x[0], grd.y[0], grd.x[grd.x.length - 1], grd.y[grd.y.length - 1]];
  else {
    segs = sectionOf(g, { axis, value });
    if (!segs.length) throw new Error('The geometry gives no outline in the chosen section plane.');
    bb = [Infinity, Infinity, -Infinity, -Infinity];
    for (const s of segs) for (const q of [0, 2]) { bb[0] = Math.min(bb[0], s[q]); bb[2] = Math.max(bb[2], s[q]); bb[1] = Math.min(bb[1], s[q + 1]); bb[3] = Math.max(bb[3], s[q + 1]); }
  }
  let sx = 1, sy = 1, ox = 0, oy = 0;                                  // window = geometry * s + o
  if (fit) {
    const w = Math.max(bb[2] - bb[0], 1e-300), h = Math.max(bb[3] - bb[1], 1e-300), W = (x1 - x0) * (1 - 2 * margin), H = (y1 - y0) * (1 - 2 * margin);
    sx = W / w; sy = H / h;
    if (!stretch) sx = sy = Math.min(sx, sy);
    ox = 0.5 * (x0 + x1) - 0.5 * (bb[0] + bb[2]) * sx; oy = 0.5 * (y0 + y1) - 0.5 * (bb[1] + bb[3]) * sy;
  }
  let mask;
  if (segs) mask = rasterize(segs.map((s) => [s[0] * sx + ox, s[1] * sy + oy, s[2] * sx + ox, s[3] * sy + oy]), x0, x1, y0, y1, nx, ny);
  else {
    const dx = (x1 - x0) / nx, dy = (y1 - y0) / ny, lv = value ?? 0, near = (a, q) => { let lo = 0, hi = a.length - 1; while (hi - lo > 1) { const m = (lo + hi) >> 1; if (a[m] <= q) lo = m; else hi = m; } return q - a[lo] <= a[hi] - q ? lo : hi; };
    mask = Array.from({ length: ny }, (_, j) => Array.from({ length: nx }, (_q, i) => {
      const gx = (x0 + (i + 0.5) * dx - ox) / sx, gy = (y0 + (j + 0.5) * dy - oy) / sy;
      if (vox) return vox.get(Math.floor((gx - vox.x0) / vox.dx), Math.floor((gy - vox.y0) / vox.dy)) === 1;
      if (gx < bb[0] || gx > bb[2] || gy < bb[1] || gy > bb[3]) return false;
      return grd.z[near(grd.y, gy)][near(grd.x, gx)] >= lv;
    }));
  }
  if (invert) for (const row of mask) for (let i = 0; i < nx; i++) row[i] = !row[i];
  return mask;
}

/**
 * Scattered xyz -> cell-centred elevation grid. Samples are binned (mean position and elevation per cell); each cell value
 * comes from an inverse-distance-weighted plane through the bins around it, so sloping data is not biased by where the
 * samples happen to fall inside a cell, and empty cells are filled from the nearest bins.
 */
function scatterGrid(xyz, nx, ny) {
  const n = Math.floor(xyz.length / 3), bb = bboxOf(xyz), w = bb.max[0] - bb.min[0] || 1, h = bb.max[1] - bb.min[1] || 1, dx = w / nx, dy = h / ny;
  const sx = new Float64Array(nx * ny), sy = new Float64Array(nx * ny), sz = new Float64Array(nx * ny), cnt = new Uint32Array(nx * ny);
  for (let p = 0; p < n; p++) {
    const xv = xyz[3 * p], yv = xyz[3 * p + 1], zv = xyz[3 * p + 2];
    if (zv !== zv) continue;
    const q = Math.min(ny - 1, Math.max(0, Math.floor((yv - bb.min[1]) / dy))) * nx + Math.min(nx - 1, Math.max(0, Math.floor((xv - bb.min[0]) / dx)));
    sx[q] += xv; sy[q] += yv; sz[q] += zv; cnt[q]++;
  }
  const z = Array.from({ length: ny }, () => new Array(nx)), rmax = Math.max(nx, ny), eps = 0.0025 * (dx * dx + dy * dy);
  let budget = 3e8;
  for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {
    const xc = bb.min[0] + (i + 0.5) * dx, yc = bb.min[1] + (j + 0.5) * dy, own = cnt[j * nx + i] > 0;
    let S = 0, Sx = 0, Sy = 0, Sxx = 0, Sxy = 0, Syy = 0, Sz = 0, Sxz = 0, Syz = 0, found = 0, lo = Infinity, hi = -Infinity;
    for (let r = 0; r <= rmax && budget > 0; r++) {
      for (let jj = j - r; jj <= j + r; jj++) {
        if (jj < 0 || jj >= ny) continue;
        const edge = jj === j - r || jj === j + r;
        for (let ii = i - r; ii <= i + r; ii += edge || r === 0 ? 1 : 2 * r) {
          budget--;
          const q = jj * nx + ii, c = cnt[q];
          if (ii < 0 || ii >= nx || !c) continue;
          const u = sx[q] / c - xc, v = sy[q] / c - yc, zm = sz[q] / c, wt = c / (u * u + v * v + eps);
          S += wt; Sx += wt * u; Sy += wt * v; Sxx += wt * u * u; Sxy += wt * u * v; Syy += wt * v * v; Sz += wt * zm; Sxz += wt * u * zm; Syz += wt * v * zm;
          found++; if (zm < lo) lo = zm; if (zm > hi) hi = zm;
        }
      }
      if (r >= 1 && (own || found >= 3)) break;
    }
    if (!found) { z[j][i] = NaN; continue; }
    const det = S * (Sxx * Syy - Sxy * Sxy) - Sx * (Sx * Syy - Sxy * Sy) + Sy * (Sx * Sxy - Sxx * Sy);
    let a = Sz / S;
    if (found >= 3 && det > 1e-9 * S * Sxx * Syy) a = (Sz * (Sxx * Syy - Sxy * Sxy) - Sx * (Sxz * Syy - Sxy * Syz) + Sy * (Sxz * Sxy - Sxx * Syz)) / det;
    z[j][i] = Math.max(lo - (hi - lo), Math.min(hi + (hi - lo), a));
  }
  return { x: Array.from({ length: nx }, (_, i) => bb.min[0] + (i + 0.5) * dx), y: Array.from({ length: ny }, (_, j) => bb.min[1] + (j + 0.5) * dy), z };
}
/**
 * Elevation grid { x[nx], y[ny], z[ny][nx] }: a DEM grid is resampled bilinearly, scattered points / soundings and
 * mesh vertices are binned and fitted locally (empty cells filled from the nearest bins), a voxel volume gives the height of its top solid.
 */
export function gridOf(g, nx, ny) {
  nx = Math.max(1, nx | 0); ny = Math.max(1, ny | 0);
  if (!g) throw new Error('No geometry given.');
  if (g.kind === 'grid') {
    const { x, y, z } = g.grid, lin = (a, n) => Array.from({ length: n }, (_, i) => (n === 1 ? 0.5 * (a[0] + a[a.length - 1]) : a[0] + ((a[a.length - 1] - a[0]) * i) / (n - 1)));
    const X = lin(x, nx), Y = lin(y, ny), loc = (a, q) => { let lo = 0, hi = a.length - 1; if (hi === 0) return [0, 0, 0]; while (hi - lo > 1) { const m = (lo + hi) >> 1; if (a[m] <= q) lo = m; else hi = m; } return [lo, hi, Math.max(0, Math.min(1, (q - a[lo]) / (a[hi] - a[lo] || 1)))]; };
    const ix = X.map((q) => loc(x, q)), iy = Y.map((q) => loc(y, q));
    return { x: X, y: Y, z: iy.map(([j0, j1, tj]) => ix.map(([i0, i1, ti]) => { let s = 0, w = 0; for (const [j, i, wt] of [[j0, i0, (1 - tj) * (1 - ti)], [j0, i1, (1 - tj) * ti], [j1, i0, tj * (1 - ti)], [j1, i1, tj * ti]]) { const v = z[j][i]; if (v === v && wt > 0) { s += wt * v; w += wt; } } return w > 0 ? s / w : NaN; })) };
  }
  if (g.kind === 'points') return scatterGrid(g.points, nx, ny);
  if (g.kind === 'mesh') return scatterGrid(g.triangles, nx, ny);
  if (g.kind === 'voxels') {
    const v = g.voxels, or = v.origin || [0, 0, 0], pts = new Float64Array(v.nx * v.ny * 3);
    for (let j = 0; j < v.ny; j++) for (let i = 0; i < v.nx; i++) { let top = 0; for (let k = v.nz - 1; k >= 0; k--) if (v.data[(k * v.ny + j) * v.nx + i]) { top = k + 1; break; } const o = 3 * (j * v.nx + i); pts[o] = or[0] + (i + 0.5) * v.spacing[0]; pts[o + 1] = or[1] + (j + 0.5) * v.spacing[1]; pts[o + 2] = or[2] + top * v.spacing[2]; }
    return scatterGrid(pts, nx, ny);
  }
  throw new Error(`A ${g.kind} geometry has no elevations to grid.`);
}

/**
 * Pore-scale statistics of a solid/pore voxel geometry (or a 2-D boolean mask, true = solid).
 * specificSurface is the voxel-face interface area per unit bulk volume; chord sizes are mean intercept lengths along the
 * axes; tortuosity[a] is the mean shortest pore path (5-7-9 chamfer metric, no corner cutting) from the inlet to the outlet
 * face along axis a divided by the straight distance, or null when the pore space does not percolate along that axis;
 * connectedPorosity counts only pore clusters that span the sample along at least one axis.
 * opt.tortuosity = false skips the path search; above 8 million voxels it only runs when opt.tortuosity = true.
 */
export function microstructure(g, opt = {}) {
  const v = asVoxels(g), { nx, ny, nz, data } = v, sp = v.spacing || [1, 1, 1], n = nx * ny * nz, d = [nx, ny, nz], st = [1, nx, nx * ny], axes = nz > 1 ? [0, 1, 2] : [0, 1];
  if (!n || data.length < n) throw new Error('The voxel geometry is empty.');
  let solid = 0, area = 0, poreLen = 0, poreRuns = 0, solidLen = 0, solidRuns = 0;
  for (let i = 0; i < n; i++) solid += data[i] ? 1 : 0;
  const coord = (i, a) => (a === 0 ? i % nx : a === 1 ? Math.floor(i / nx) % ny : Math.floor(i / st[2]));
  for (const a of axes) {                       // walk every voxel line along axis a: interface faces and chord runs
    const fa = (sp[0] * sp[1] * sp[2]) / sp[a], [b, c] = [0, 1, 2].filter((k) => k !== a), step = st[a], da = d[a];
    let faces = 0, solidVox = 0, runsS = 0, runsP = 0;
    for (let q = 0; q < d[c]; q++) for (let p = 0; p < d[b]; p++) {
      let prev = -1;
      for (let t = 0, i = p * st[b] + q * st[c]; t < da; t++, i += step) {
        const s = data[i] ? 1 : 0;
        if (s !== prev) { if (prev >= 0) faces++; if (s) runsS++; else runsP++; prev = s; }
        solidVox += s;
      }
    }
    area += faces * fa; solidRuns += runsS; poreRuns += runsP; solidLen += solidVox * sp[a]; poreLen += (n - solidVox) * sp[a];
  }
  // 6-connected pore clusters and the axes they span
  const label = new Uint8Array(n), queue = new Int32Array(n - solid), percolates = axes.map(() => false), nxy = st[2];
  let nlab = 0, connected = 0;
  for (let s0 = 0; s0 < n; s0++) {
    if (data[s0] || label[s0]) continue;
    nlab++;
    let head = 0, tail = 0, i0 = nx, i1 = -1, j0 = ny, j1 = -1, k0 = nz, k1 = -1;
    queue[tail++] = s0; label[s0] = 1;
    while (head < tail) {
      const i = queue[head++], ck = (i / nxy) | 0, r = i - ck * nxy, cj = (r / nx) | 0, ci = r - cj * nx;
      if (ci < i0) i0 = ci; if (ci > i1) i1 = ci; if (cj < j0) j0 = cj; if (cj > j1) j1 = cj; if (ck < k0) k0 = ck; if (ck > k1) k1 = ck;
      if (ci > 0 && !data[i - 1] && !label[i - 1]) { label[i - 1] = 1; queue[tail++] = i - 1; }
      if (ci + 1 < nx && !data[i + 1] && !label[i + 1]) { label[i + 1] = 1; queue[tail++] = i + 1; }
      if (cj > 0 && !data[i - nx] && !label[i - nx]) { label[i - nx] = 1; queue[tail++] = i - nx; }
      if (cj + 1 < ny && !data[i + nx] && !label[i + nx]) { label[i + nx] = 1; queue[tail++] = i + nx; }
      if (ck > 0 && !data[i - nxy] && !label[i - nxy]) { label[i - nxy] = 1; queue[tail++] = i - nxy; }
      if (ck + 1 < nz && !data[i + nxy] && !label[i + nxy]) { label[i + nxy] = 1; queue[tail++] = i + nxy; }
    }
    const span = [i0 === 0 && i1 === nx - 1, j0 === 0 && j1 === ny - 1, k0 === 0 && k1 === nz - 1];
    let spans = false;
    axes.forEach((a, q) => { if (span[a]) { percolates[q] = true; spans = true; } });
    if (spans) connected += tail;
  }
  // chamfer geodesic distance from the inlet face (Dial's bucket queue)
  const smin = Math.min(...axes.map((a) => sp[a])), mv = [];
  for (let dz = nz > 1 ? -1 : 0; dz <= (nz > 1 ? 1 : 0); dz++) for (let dy = -1; dy <= 1; dy++) for (let dx = -1; dx <= 1; dx++) {
    if (!dx && !dy && !dz) continue;
    const side = [dx, dy * st[1], dz * st[2]].filter(Boolean);
    mv.push([dx, dy, dz, dx + dy * st[1] + dz * st[2], Math.max(1, Math.round((5 * Math.hypot(dx * sp[0], dy * sp[1], dz * sp[2])) / smin)), ...(side.length > 1 ? side : []), 0, 0, 0]);
  }
  const nm = mv.length, col = (k) => Int32Array.from(mv, (m) => m[k]), mdx = col(0), mdy = col(1), mdz = col(2), moff = col(3), mw = col(4), ms0 = col(5), ms1 = col(6), ms2 = col(7);
  const nb = Math.max(...mw) + 1, INF = 0x7fffffff;
  const skipPaths = opt.tortuosity === false || (n > 8e6 && opt.tortuosity !== true);   // large volumes: only on request
  const tortuosity = axes.map((a, q) => {
    if (!percolates[q] || skipPaths) return null;
    if (d[a] < 2) return 1;
    const dist = new Int32Array(n).fill(INF), buckets = Array.from({ length: nb }, () => []);
    let pending = 0;
    for (let i = 0; i < n; i++) if (!data[i] && coord(i, a) === 0) { dist[i] = 0; buckets[0].push(i); pending++; }
    for (let lev = 0; pending > 0; lev++) {
      const cur = buckets[lev % nb];
      if (!cur.length) continue;
      buckets[lev % nb] = [];
      pending -= cur.length;
      for (let c = 0; c < cur.length; c++) {
        const i = cur[c];
        if (dist[i] !== lev) continue;
        const ck = (i / nxy) | 0, r = i - ck * nxy, cj = (r / nx) | 0, ci = r - cj * nx;
        for (let m = 0; m < nm; m++) {
          const ii = ci + mdx[m], jj = cj + mdy[m], kk = ck + mdz[m];
          if (ii < 0 || jj < 0 || kk < 0 || ii >= nx || jj >= ny || kk >= nz) continue;
          const u = i + moff[m], nd = lev + mw[m];
          if (data[u] || nd >= dist[u]) continue;
          const s0 = ms0[m];
          if (s0 !== 0 && (data[i + s0] || data[i + ms1[m]] || (ms2[m] !== 0 && data[i + ms2[m]]))) continue;   // no squeezing through solid corners
          dist[u] = nd; buckets[nd % nb].push(u); pending++;
        }
      }
    }
    let sum = 0, cnt = 0, wa = 5;
    for (let i = 0; i < n; i++) if (!data[i] && coord(i, a) === d[a] - 1 && dist[i] < INF) { sum += dist[i]; cnt++; }
    for (let m = 0; m < nm; m++) if ([mdx[m], mdy[m], mdz[m]][a] === 1 && Math.abs(mdx[m]) + Math.abs(mdy[m]) + Math.abs(mdz[m]) === 1) wa = mw[m];
    return cnt ? sum / cnt / (wa * (d[a] - 1)) : null;
  });
  const bulk = n * sp[0] * sp[1] * sp[2];
  return { porosity: 1 - solid / n, solidFraction: solid / n, specificSurface: area / bulk, meanPoreSize: poreRuns ? poreLen / poreRuns : 0, meanSolidSize: solidRuns ? solidLen / solidRuns : 0, tortuosity, percolates, connectedPorosity: connected / n, poreClusters: nlab, tortuositySkipped: skipPaths };
}

/** Totals, pipe list and (for acyclic graphs) a topological order of a network geometry. */
export function networkSummary(g) {
  if (!g || g.kind !== 'network') throw new Error('A network geometry is required.');
  const { nodes, edges } = g.network, byType = {}, num = (v) => (typeof v === 'number' && Number.isFinite(v) ? v : 0);
  let totalLength = 0, zlo = Infinity, zhi = -Infinity;
  for (const n of nodes) { byType[n.type] = (byType[n.type] || 0) + 1; totalLength += num(n.length); if (typeof n.z === 'number' && Number.isFinite(n.z)) { zlo = Math.min(zlo, n.z); zhi = Math.max(zhi, n.z); } }
  for (const e of edges) { if (!nodes.length || e.type !== 'connection') byType[e.type] = (byType[e.type] || 0) + 1; totalLength += num(e.length); }
  const pipes = edges.map((e) => ({ from: e.from, to: e.to, length: e.length ?? null, diameter: e.diameter ?? null, name: e.name }));
  for (const n of nodes) if (n.type === 'pipe' && n.length !== undefined) pipes.push({ from: n.id, to: n.id, length: n.length ?? null, diameter: n.diameter ?? null, name: n.name });
  const indeg = new Map(nodes.map((n) => [n.id, 0])), out = new Map(nodes.map((n) => [n.id, []]));
  for (const e of edges) if (indeg.has(e.to) && out.has(e.from)) { indeg.set(e.to, indeg.get(e.to) + 1); out.get(e.from).push(e.to); }
  const ready = nodes.filter((n) => !indeg.get(n.id)).map((n) => n.id), order = [];
  for (let q = 0; q < ready.length; q++) { const id = ready[q]; order.push(id); for (const t of out.get(id)) { indeg.set(t, indeg.get(t) - 1); if (!indeg.get(t)) ready.push(t); } }
  const acyclic = order.length === nodes.length;
  return { nodes: nodes.length, edges: edges.length, totalLength, byType, pipes, maxElevationChange: zhi >= zlo ? zhi - zlo : 0, order: acyclic ? order : [], acyclic, sources: nodes.filter((n) => !edges.some((e) => e.to === n.id)).map((n) => n.id), sinks: nodes.filter((n) => !edges.some((e) => e.from === n.id)).map((n) => n.id) };
}

/**
 * Overall size, surface area, enclosed volume and closedness of a geometry, with a unit hint when the file gave one.
 * Meshes: area of all triangles, volume by the divergence theorem when every edge is shared by exactly two triangles.
 * Voxels: solid volume by count and the solid/pore interface area inside the sample (faces on the sample boundary excluded).
 * Polylines: enclosed area of the closed outlines and total length. Networks: total pipe length.
 */
export function dimensions(g) {
  if (!g || !g.bbox) throw new Error('No geometry given.');
  const size = [0, 1, 2].map((k) => (Number.isFinite(g.bbox.max[k]) && Number.isFinite(g.bbox.min[k]) ? g.bbox.max[k] - g.bbox.min[k] : 0));
  const out = { size, area: null, volume: null, closed: false, units: (g.stats && g.stats.units) || (g.geographic || (g.grid && g.grid.geographic) ? 'degrees' : null) };
  if (g.kind === 'mesh') {
    const t = g.triangles, nt = t.length / 9, nv3 = nt * 3, diag = Math.hypot(...size) || 1, inv = 1e7 / diag, o = g.bbox.min, vid = new Int32Array(nv3), edges = new Map();
    // weld vertices on a 1e-7·diagonal lattice with an open-addressing hash table
    let cap = 16, nvert = 0, area = 0;
    while (cap < nv3 * 2) cap *= 2;
    const slot = new Int32Array(cap).fill(-1), qa = new Int32Array(nv3), qb = new Int32Array(nv3), qc = new Int32Array(nv3);
    for (let i = 0; i < nv3; i++) {
      const a = Math.round((t[3 * i] - o[0]) * inv + 0.3183) | 0, b = Math.round((t[3 * i + 1] - o[1]) * inv + 0.3183) | 0, c = Math.round((t[3 * i + 2] - o[2]) * inv + 0.3183) | 0;
      for (let h = (Math.imul(a, 73856093) ^ Math.imul(b, 19349663) ^ Math.imul(c, 83492791)) & (cap - 1); ; h = (h + 1) & (cap - 1)) {
        const s = slot[h];
        if (s < 0) { slot[h] = nvert; qa[nvert] = a; qb[nvert] = b; qc[nvert] = c; vid[i] = nvert++; break; }
        if (qa[s] === a && qb[s] === b && qc[s] === c) { vid[i] = s; break; }
      }
    }
    for (let f = 0; f < nt; f++) {
      const o = 9 * f, a = [t[o + 3] - t[o], t[o + 4] - t[o + 1], t[o + 5] - t[o + 2]], b = [t[o + 6] - t[o], t[o + 7] - t[o + 1], t[o + 8] - t[o + 2]];
      area += 0.5 * vlen(cross(a, b));
      for (let e = 0; e < 3; e++) { const p = vid[3 * f + e], q = vid[3 * f + ((e + 1) % 3)]; if (p === q) continue; const k = p < q ? p * nvert + q : q * nvert + p; edges.set(k, (edges.get(k) || 0) + (p < q ? 1 : 65536)); }
    }
    let closed = edges.size > 0, oriented = true;
    for (const c of edges.values()) { const fwd = c & 65535, bwd = c >> 16; if (fwd + bwd !== 2) { closed = false; break; } if (fwd !== 1) oriented = false; }
    out.area = area; out.closed = closed; out.vertices = nvert;
    if (closed) { out.volume = Math.abs(signedVolume(t)); out.consistentNormals = oriented; }
  } else if (g.kind === 'voxels') {
    const v = g.voxels, m = microstructure(g, { tortuosity: false }), bulk = v.nx * v.ny * v.nz * v.spacing[0] * v.spacing[1] * v.spacing[2];
    out.volume = m.solidFraction * bulk; out.area = m.specificSurface * bulk; out.closed = true; out.porosity = m.porosity;
  } else if (g.kind === 'polylines') {
    let area = 0, len = 0, allClosed = g.polylines.length > 0;
    for (const p of g.polylines) {
      const n = p.x.length;
      for (let i = 0; i + 1 < n; i++) len += Math.hypot(p.x[i + 1] - p.x[i], p.y[i + 1] - p.y[i], p.z ? p.z[i + 1] - p.z[i] : 0);
      if (p.closed && n > 2) { len += Math.hypot(p.x[0] - p.x[n - 1], p.y[0] - p.y[n - 1]); let s = 0; for (let i = 0; i < n; i++) { const j = (i + 1) % n; s += p.x[i] * p.y[j] - p.x[j] * p.y[i]; } area += Math.abs(s) / 2; } else allClosed = false;
    }
    out.area = area; out.length = len; out.closed = allClosed;
  } else if (g.kind === 'grid') out.area = size[0] * size[1];
  else if (g.kind === 'network') out.length = networkSummary(g).totalLength;
  return out;
}

// ---- Procedural geometry -----------------------------------------------------------------------------------
const TPMS = {
  gyroid: (x, y, z) => Math.sin(x) * Math.cos(y) + Math.sin(y) * Math.cos(z) + Math.sin(z) * Math.cos(x),
  schwarzp: (x, y, z) => Math.cos(x) + Math.cos(y) + Math.cos(z),
  diamond: (x, y, z) => Math.sin(x) * Math.sin(y) * Math.sin(z) + Math.sin(x) * Math.cos(y) * Math.cos(z) + Math.cos(x) * Math.sin(y) * Math.cos(z) + Math.cos(x) * Math.cos(y) * Math.sin(z),
};
const vec3 = (v, dflt) => (Array.isArray(v) ? [0, 1, 2].map((k) => (Number.isFinite(+v[k]) ? +v[k] : dflt[k])) : Number.isFinite(+v) && v !== null && v !== undefined && v !== '' ? [+v, +v, +v] : dflt.slice());
/** Signed-distance function of a JSON CSG tree (negative inside). */
function csgCompile(node, boxes, depth = 0, count = { n: 0 }) {
  if (!node || typeof node !== 'object' || depth > 64 || ++count.n > 20000) throw new Error('The CSG tree is missing, too deep or too large.');
  const kind = String(node.op || node.type || node.shape || '').toLowerCase(), kids = node.children || node.args || node.items;
  if (['union', 'intersection', 'difference', 'subtract', 'intersect', 'add'].includes(kind)) {
    if (!Array.isArray(kids) || !kids.length) throw new Error(`CSG operation "${kind}" needs children.`);
    const neg = kind === 'difference' || kind === 'subtract', fs = kids.map((k, q) => csgCompile(k, neg && q > 0 ? [] : boxes, depth + 1, count));
    if (kind === 'union' || kind === 'add') return (x, y, z) => { let m = Infinity; for (const f of fs) m = Math.min(m, f(x, y, z)); return m; };
    if (neg) return (x, y, z) => { let m = fs[0](x, y, z); for (let q = 1; q < fs.length; q++) m = Math.max(m, -fs[q](x, y, z)); return m; };
    return (x, y, z) => { let m = -Infinity; for (const f of fs) m = Math.max(m, f(x, y, z)); return m; };
  }
  if (kind === 'box' || kind === 'cuboid' || kind === 'cube') {
    const lo = node.min ? vec3(node.min, [0, 0, 0]) : null, sz = vec3(node.size, [1, 1, 1]), c = lo ? vec3(node.max, [1, 1, 1]).map((v, k) => 0.5 * (v + lo[k])) : vec3(node.center || node.centre, [0, 0, 0]), h = lo ? vec3(node.max, [1, 1, 1]).map((v, k) => 0.5 * (v - lo[k])) : sz.map((v) => v / 2);
    boxes.push([c.map((v, k) => v - h[k]), c.map((v, k) => v + h[k])]);
    return (x, y, z) => { const qx = Math.abs(x - c[0]) - h[0], qy = Math.abs(y - c[1]) - h[1], qz = Math.abs(z - c[2]) - h[2]; return Math.hypot(Math.max(qx, 0), Math.max(qy, 0), Math.max(qz, 0)) + Math.min(Math.max(qx, qy, qz), 0); };
  }
  if (kind === 'sphere' || kind === 'ball') {
    const c = vec3(node.center || node.centre, [0, 0, 0]), r = +(node.radius ?? node.r ?? (node.diameter ?? 1) / 2);
    boxes.push([c.map((v) => v - r), c.map((v) => v + r)]);
    return (x, y, z) => Math.hypot(x - c[0], y - c[1], z - c[2]) - r;
  }
  if (kind === 'cylinder') {
    const r = +(node.radius ?? node.r ?? (node.diameter ?? 1) / 2), hgt = +(node.height ?? node.h ?? node.length ?? 1), a = Math.max(0, 'xyz'.indexOf(String(node.axis ?? 'z').toLowerCase())), c = node.base ? vec3(node.base, [0, 0, 0]).map((v, k) => (k === a ? v + hgt / 2 : v)) : vec3(node.center || node.centre, [0, 0, 0]), [u, w] = [0, 1, 2].filter((k) => k !== a);
    boxes.push([c.map((v, k) => v - (k === a ? hgt / 2 : r)), c.map((v, k) => v + (k === a ? hgt / 2 : r))]);
    return (x, y, z) => { const p = [x, y, z], qr = Math.hypot(p[u] - c[u], p[w] - c[w]) - r, qa = Math.abs(p[a] - c[a]) - hgt / 2; return Math.hypot(Math.max(qr, 0), Math.max(qa, 0)) + Math.min(Math.max(qr, qa), 0); };
  }
  throw new Error(`Unknown CSG node "${kind}" (use union, difference, intersection, box, sphere or cylinder).`);
}
/**
 * Line generators -> 3-D polylines with x = horizontal distance (m), y = 0 and z = elevation (m, negative below sea level),
 * ordered in the production direction except the well, which runs from the wellhead down like a survey.
 *   catenary  waterDepth, topAngle (hang-off from vertical, deg), flowline (seabed lead-in length), n
 *   lazywave  as catenary plus sagHeight, hogHeight (above seabed) and buoyancyRatio (net uplift / submerged weight)
 *   flowline  length, amplitude, wavelength, slope (deg, + uphill), waterDepth (at the start), seed, n
 *   well      kickoff (m MD), buildRate (deg / 30 m), inclination (hold angle), md (total), azimuth, step
 *   jumper    span, height, dip, leg (fraction of the span before the first drop), waterDepth
 */
function generateLine(kind, spec, seed) {
  const pos = (v, d) => (Number.isFinite(+v) && +v > 0 ? +v : d), num = (v, d) => (v !== null && v !== undefined && v !== '' && Number.isFinite(+v) ? +v : d), warnings = [], stats = { type: kind }, D = Math.PI / 180;
  const depth = pos(spec.waterDepth ?? spec.depth, 1000), x = [], z = [];
  let g = null;
  if (kind === 'catenary' || kind === 'lazywave') {
    const hang = Math.max(1, Math.min(75, num(spec.topAngle ?? spec.hangoff, 12))), t = Math.tan((90 - hang) * D), ct = Math.sqrt(1 + t * t), lead = Math.max(0, num(spec.flowline ?? spec.lead, 0.25 * depth)), n = Math.max(8, Math.min(4000, Math.round(pos(spec.n, 200))));
    let pieces;               // catenary arcs { a, sa, sb }: slope runs from sa to sb, concave up when sb > sa
    if (kind === 'catenary') { const a = depth / (ct - 1); pieces = [{ a, sa: 0, sb: t }]; stats.catenaryParameter = a; }
    else {
      const hSag = Math.min(0.8 * depth, pos(spec.sagHeight, 0.12 * depth)), hHog = Math.min(0.95 * depth, Math.max(pos(spec.hogHeight, 0.25 * depth), hSag * 1.05)), r = pos(spec.buoyancyRatio, 1);
      const a = (depth - hSag) / (ct - 1), a2 = a / r, c1 = 1 + hHog / (a + a2), c3 = 1 + (hHog - hSag) / (a + a2), s1 = Math.sqrt(c1 * c1 - 1), s3 = Math.sqrt(c3 * c3 - 1);
      pieces = [{ a, sa: 0, sb: s1 }, { a: a2, sa: s1, sb: 0 }, { a: a2, sa: 0, sb: -s3 }, { a, sa: -s3, sb: 0 }, { a, sa: 0, sb: t }];
      Object.assign(stats, { catenaryParameter: a, buoyantParameter: a2, sagHeight: hSag, hogHeight: hHog, buoyantLength: a2 * (s1 + s3) });
    }
    const susp = pieces.reduce((q, p) => q + p.a * Math.abs(p.sb - p.sa), 0);
    if (lead > 0) { x.push(0); z.push(-depth); }
    let x0 = lead, z0 = -depth;
    x.push(x0); z.push(z0);
    for (const p of pieces) {
      const sg = p.sb > p.sa ? 1 : -1, w0 = sg * Math.asinh(p.sa), w1 = sg * Math.asinh(p.sb), m = Math.max(2, Math.round((n * p.a * Math.abs(p.sb - p.sa)) / susp));
      for (let k = 1; k <= m; k++) { const w = w0 + ((w1 - w0) * k) / m; x.push(x0 + p.a * (w - w0)); z.push(z0 + sg * p.a * (Math.cosh(w) - Math.cosh(w0))); }
      x0 = x[x.length - 1]; z0 = z[z.length - 1];
    }
    z[z.length - 1] = 0;
    Object.assign(stats, { waterDepth: depth, topAngle: hang, touchdownX: lead, suspendedLength: susp, length: lead + susp, horizontalLength: x0 });
  } else if (kind === 'flowline') {
    const len = pos(spec.length, 5000), A = Math.abs(num(spec.amplitude, 20)), lam = pos(spec.wavelength, 500), slope = Math.max(-45, Math.min(45, num(spec.slope, 0))), n = Math.max(50, Math.min(4000, Math.round(pos(spec.n, (len / lam) * 24)))), R = rng(seed);
    const hs = [[1, 0.5], [0.61, 0.25], [1.7, 0.15], [0.37, 0.1]].map(([f, w]) => ({ k: TAU / (lam * f * R.uniform(0.9, 1.1)), w, ph: R.uniform(0, TAU) }));
    for (let i = 0; i <= n; i++) { const xi = (len * i) / n; x.push(xi); z.push(-depth + xi * Math.tan(slope * D) + A * hs.reduce((q, h) => q + h.w * Math.sin(h.k * xi + h.ph), 0)); }
    Object.assign(stats, { seed, horizontalLength: len, amplitude: A, wavelength: lam, slope, waterDepth: depth });
  } else if (kind === 'well') {
    const kop = Math.max(0, num(spec.kickoff ?? spec.kop, 500)), bur = pos(spec.buildRate ?? spec.bur, 3), hold = Math.max(0, Math.min(95, num(spec.inclination ?? spec.hold, 45))), mdT = pos(spec.md ?? spec.length, 3000), az = num(spec.azimuth, 0), step = Math.max(mdT / 2000, pos(spec.step, 30)), eob = kop + (hold / bur) * 30;
    const mds = new Set([0, mdT, ...(kop < mdT ? [kop] : []), ...(eob < mdT ? [eob] : [])]);
    for (let m = step; m < mdT; m += step) mds.add(m);
    const md = [...mds].sort((a, b) => a - b).filter((m, i, a) => !i || m - a[i - 1] > 1e-9), inc = md.map((m) => (m <= kop ? 0 : Math.min(hold, ((m - kop) * bur) / 30)));
    g = surveyGeom({ md, inc, azi: md.map(() => az) }, { warnings });
    Object.assign(g.stats, { type: kind, kickoff: kop, buildRate: bur, endOfBuild: eob, azimuth: az });
    if (eob > mdT) warnings.push(`The hold angle is not reached: the build section would end at MD ${+eob.toFixed(1)} m, beyond the total depth.`);
  } else {
    const span = pos(spec.span, 30), h = pos(spec.height, 8), dip = Math.max(0, Math.min(h, num(spec.dip, 0.4 * h))), w1 = span * Math.max(0.02, Math.min(0.49, num(spec.leg, 0.2)));
    x.push(0, 0, w1, w1, span - w1, span - w1, span, span); z.push(-depth, -depth + h, -depth + h, -depth + h - dip, -depth + h - dip, -depth + h, -depth + h, -depth);
    if (!dip) { x.splice(3, 2); z.splice(3, 2); }
    Object.assign(stats, { span, height: h, dip, length: 2 * h + 2 * dip + span, waterDepth: depth });
    warnings.push('Bends are drawn as sharp corners; the bend radius is not modelled.');
  }
  if (!g) { g = polyGeom([{ x, y: x.map(() => 0), z, closed: false }], { warnings }); g.stats = stats; }
  g.stats.polylines = 1; g.stats.vertices = g.polylines[0].x.length;
  return { name: spec.name || `Procedural ${kind}`, format: 'Procedural: ' + kind, pathway: 'procedural', ...g, warnings };
}
/** Names accepted by generate() for the line generators, e.g. "catenary riser", "lazy-wave riser", "undulating flowline". */
const lineKind = (type) => { const k = type.replace(/[^a-z]/g, ''); return /^catenary/.test(k) ? 'catenary' : /^lazy/.test(k) ? 'lazywave' : /^(undulating|flowline)/.test(k) ? 'flowline' : /^well/.test(k) ? 'well' : /^jumper/.test(k) ? 'jumper' : null; };

/**
 * Procedural geometry. Voxel types (2-D when nz = 1): 'tpms' | 'voronoi' | 'lattice' | 'spheres' | 'csg' | 'implicit' |
 * 'spacer' with common fields n (voxels along x, or [nx, ny, nz]), size ([lx, ly, lz]) and seed. Line types (3-D polylines,
 * see generateLine): 'catenary riser' | 'lazy-wave riser' | 'undulating flowline' | 'well trajectory (build-hold)' |
 * 'jumper (M-shape)'. spec.fn (type 'implicit') must be a function supplied by application code; strings are never evaluated.
 */
export function generate(spec = {}) {
  if (!spec || typeof spec !== 'object') throw new Error('generate() needs a specification object.');
  const type = String(spec.type || '').toLowerCase(), seed = Number.isFinite(+spec.seed) ? +spec.seed : 1, R = rng(seed), warnings = [], stats = { type, seed };
  if (lineKind(type)) return generateLine(lineKind(type), spec, seed);
  const pos = (v, dflt) => (Number.isFinite(+v) && +v > 0 ? +v : dflt), clampPhi = (p) => Math.max(0.005, Math.min(0.995, +p));
  let size = vec3(spec.size, [1, 1, 1]), origin = [0, 0, 0], dimsGiven = Array.isArray(spec.n) ? spec.n.map((k) => Math.max(1, k | 0)) : null, n0 = pos(Array.isArray(spec.n) ? NaN : spec.n ?? spec.resolution, 0), two = !!spec.twoD || (dimsGiven && dimsGiven[2] === 1);
  let solidAt = null, field = null, target = null;      // solidAt(x, y, z) -> bool, or field(x, y, z) with a quantile target
  if (type === 'tpms') {
    const f = TPMS[String(spec.surface || spec.kind || 'gyroid').toLowerCase().replace(/[^a-z]/g, '')];
    if (!f) throw new Error('Unknown TPMS surface (use gyroid, schwarzP or diamond).');
    const per = vec3(spec.period, [size[0], size[0], size[0]]), k = per.map((p) => TAU / p), sheet = !!spec.sheet, lvl = Number.isFinite(+spec.level) ? +spec.level : sheet ? 0.3 : 0;
    stats.surface = String(spec.surface || spec.kind || 'gyroid'); stats.period = per;
    field = sheet ? (x, y, z) => Math.abs(f(k[0] * x, k[1] * y, k[2] * z)) : (x, y, z) => f(k[0] * x, k[1] * y, k[2] * z);
    if (spec.porosity !== undefined) target = 1 - clampPhi(spec.porosity); else { const fl = field; field = null; solidAt = (x, y, z) => fl(x, y, z) < lvl; stats.level = lvl; }
  } else if (type === 'voronoi') {
    const nc = Math.max(2, Math.min(5000, (spec.cells ?? spec.count ?? 20) | 0)), wall = pos(spec.wall ?? spec.thickness, 0.03 * size[0]), s = [];
    for (let q = 0; q < nc; q++) s.push([R.uniform(0, size[0]), R.uniform(0, size[1]), R.uniform(0, size[2])]);
    stats.cells = nc;
    field = (x, y, z) => {
      let d1 = Infinity, d2 = Infinity, i1 = 0, i2 = 0;
      for (let q = 0; q < nc; q++) { const c = s[q], dq = (x - c[0]) ** 2 + (y - c[1]) ** 2 + (two ? 0 : (z - c[2]) ** 2); if (dq < d1) { d2 = d1; i2 = i1; d1 = dq; i1 = q; } else if (dq < d2) { d2 = dq; i2 = q; } }
      const a = s[i1], b = s[i2], sep = Math.hypot(a[0] - b[0], a[1] - b[1], two ? 0 : a[2] - b[2]);
      return sep > 0 ? (d2 - d1) / (2 * sep) : 0;
    };
    if (spec.porosity !== undefined) target = 1 - clampPhi(spec.porosity); else { const fl = field; field = null; solidAt = (x, y, z) => fl(x, y, z) < wall / 2; stats.wall = wall; }
    stats.work = nc;
  } else if (type === 'lattice') {
    const p = pos(spec.pitch, 1), dia = pos(spec.diameter, 0.25 * p), ang = ((Number.isFinite(+spec.angle) ? +spec.angle : 90) * Math.PI) / 180, H = pos(spec.height, 2 * dia), hs = Math.sin(ang / 2), hc = Math.cos(ang / 2);
    if (spec.size === undefined) size = [hs > 1e-6 ? (2 * p) / hs : 2 * p, hc > 1e-6 ? (2 * p) / hc : 2 * p, H]; else size[2] = spec.size[2] ?? H;
    const fam = (x, y, z, sgn, zc) => { const sC = -x * hs * sgn + y * hc, ds = sC - p * Math.round(sC / p); return Math.hypot(ds, z - zc); };
    solidAt = (x, y, z) => fam(x, y, z, 1, dia / 2) < dia / 2 || fam(x, y, z, -1, H - dia / 2) < dia / 2;
    Object.assign(stats, { pitch: p, diameter: dia, angleDeg: (ang * 180) / Math.PI, height: H });
  } else if (type === 'spacer') {
    const H = pos(spec.height, 1), Lc = pos(spec.length, 8 * H), dia = pos(spec.diameter, 0.5 * H), p = pos(spec.pitch, 4 * dia), arr = String(spec.arrangement || 'zigzag').toLowerCase();
    size = [Lc, H, Lc / pos(n0, 400)]; two = true;
    if (!dimsGiven) { const nx = Math.round(pos(n0, 400)); dimsGiven = [nx, Math.max(16, Math.round((nx * H) / Lc)), 1]; }
    solidAt = (x, y) => { const k = Math.floor(x / p), xc = (k + 0.5) * p, yc = arr === 'submerged' ? H / 2 : arr === 'cavity' || arr === 'bottom' ? dia / 2 : arr === 'top' ? H - dia / 2 : k % 2 ? H - dia / 2 : dia / 2; return Math.hypot(x - xc, y - yc) < dia / 2; };
    Object.assign(stats, { height: H, length: Lc, diameter: dia, pitch: p, arrangement: arr, filaments: Math.floor(Lc / p) });
  } else if (type === 'csg') {
    const boxes = [], sdf = csgCompile(spec.tree || spec.root || spec.csg, boxes), b = spec.bounds;
    let lo, hi;
    if (b && b.min && b.max) { lo = vec3(b.min, [0, 0, 0]); hi = vec3(b.max, [1, 1, 1]); }
    else { if (!boxes.length) throw new Error('The CSG tree has no primitives.'); lo = [0, 1, 2].map((k) => Math.min(...boxes.map((q) => q[0][k]))); hi = [0, 1, 2].map((k) => Math.max(...boxes.map((q) => q[1][k]))); const pad = 0.04 * Math.max(...hi.map((v, k) => v - lo[k])); lo = lo.map((v) => v - pad); hi = hi.map((v) => v + pad); }
    origin = lo; size = hi.map((v, k) => v - lo[k]);
    if (!size.every((v) => v > 0 && Number.isFinite(v))) throw new Error('The CSG bounds are empty.');
    solidAt = (x, y, z) => sdf(x, y, z) < 0;
  } else if (type === 'implicit') {
    if (typeof spec.fn !== 'function') throw new Error('An implicit geometry needs spec.fn as a function (x, y, z) => number; text is never evaluated.');
    const b = spec.bounds;
    if (b && b.min && b.max) { origin = vec3(b.min, [0, 0, 0]); size = vec3(b.max, [1, 1, 1]).map((v, k) => v - origin[k]); }
    if (!size.every((v) => v > 0 && Number.isFinite(v))) throw new Error('The implicit-geometry bounds are empty.');
    const fn = spec.fn;
    solidAt = (x, y, z) => { const v = fn(x, y, z); return typeof v === 'number' && v < 0; };
  } else if (type !== 'spheres') throw new Error(`Unknown procedural geometry type "${spec.type}" (use tpms, voronoi, lattice, spheres, csg, implicit, spacer, catenary, lazywave, flowline, well or jumper).`);
  // grid
  let nx, ny, nz;
  if (dimsGiven) { [nx, ny] = dimsGiven; nz = two ? 1 : dimsGiven[2] || 1; ny = ny || nx; }
  else { nx = Math.round(pos(n0, two ? 256 : 64)); ny = Math.max(1, Math.round((nx * size[1]) / size[0])); nz = two ? 1 : Math.max(1, Math.round((nx * size[2]) / size[0])); }
  const n = nx * ny * nz;
  if (!(n > 0) || n > L.voxels) throw new Error(`The requested grid of ${nx} × ${ny} × ${nz} voxels is empty or too large (limit ${L.voxels / 1e6} million).`);
  if (stats.work && n * stats.work > 6e8) throw new Error('The Voronoi request is too heavy; lower the resolution or the cell count.');
  delete stats.work;
  const sp = [size[0] / nx, size[1] / ny, nz === 1 ? (two && type !== 'spacer' ? Math.min(size[0] / nx, size[1] / ny) : size[2]) : size[2] / nz], data = new Uint8Array(n);
  if (nz === 1 && type !== 'spacer') sp[2] = Math.min(sp[0], sp[1]);
  const zOf = (k) => (nz === 1 ? origin[2] + (two ? 0 : size[2] / 2) : origin[2] + (k + 0.5) * sp[2]);
  if (type === 'spheres') {
    const phi = clampPhi(spec.porosity ?? 0.4), r = pos(spec.radius ?? (spec.diameter ? spec.diameter / 2 : NaN), Math.min(size[0], size[1]) / 12), want = Math.round((1 - phi) * n), overlap = spec.overlap !== false, centres = [];
    let solid = 0, tries = 0, placed = 0;
    const maxTries = overlap ? 2e6 : 2e5;
    while (solid < want && tries++ < maxTries) {
      const c = [R.uniform(-r, size[0] + r), R.uniform(-r, size[1] + r), nz === 1 ? 0 : R.uniform(-r, size[2] + r)];
      if (!overlap) { if (centres.some((q) => Math.hypot(q[0] - c[0], q[1] - c[1], q[2] - c[2]) < 2 * r)) continue; centres.push(c); }
      placed++;
      const i0 = Math.max(0, Math.floor((c[0] - r) / sp[0])), i1 = Math.min(nx - 1, Math.ceil((c[0] + r) / sp[0])), j0 = Math.max(0, Math.floor((c[1] - r) / sp[1])), j1 = Math.min(ny - 1, Math.ceil((c[1] + r) / sp[1])), k0 = nz === 1 ? 0 : Math.max(0, Math.floor((c[2] - r) / sp[2])), k1 = nz === 1 ? 0 : Math.min(nz - 1, Math.ceil((c[2] + r) / sp[2]));
      for (let k = k0; k <= k1; k++) for (let j = j0; j <= j1; j++) for (let i = i0; i <= i1; i++) {
        const q = (k * ny + j) * nx + i;
        if (!data[q] && ((i + 0.5) * sp[0] - c[0]) ** 2 + ((j + 0.5) * sp[1] - c[1]) ** 2 + (nz === 1 ? 0 : ((k + 0.5) * sp[2] - c[2]) ** 2) < r * r) { data[q] = 1; solid++; }
      }
    }
    if (solid < want) warnings.push(`Target porosity ${phi} was not reached${overlap ? '' : ' without overlap (random sequential packing saturates)'}; porosity is ${+(1 - solid / n).toFixed(4)}.`);
    Object.assign(stats, { radius: r, spheres: placed, targetPorosity: phi, overlap });
  } else if (field) {
    const f = new Float32Array(n);
    let q = 0;
    for (let k = 0; k < nz; k++) for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) f[q++] = field(origin[0] + (i + 0.5) * sp[0], origin[1] + (j + 0.5) * sp[1], zOf(k));
    const sorted = Float32Array.from(f).sort(), kk = Math.max(1, Math.min(n - 1, Math.round(target * n))), lvl = sorted[kk - 1];
    let quota = kk;                              // values tied at the level are admitted in index order until the target count is met
    for (let i = 0; i < n; i++) if (f[i] < lvl) { data[i] = 1; quota--; }
    for (let i = 0; i < n && quota > 0; i++) if (f[i] === lvl) { data[i] = 1; quota--; }
    Object.assign(stats, { level: lvl, targetPorosity: 1 - target });
  } else {
    let q = 0;
    try { for (let k = 0; k < nz; k++) for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) data[q++] = solidAt(origin[0] + (i + 0.5) * sp[0], origin[1] + (j + 0.5) * sp[1], zOf(k)) ? 1 : 0; }
    catch (e) { throw new Error('The geometry function failed while being sampled: ' + (e && e.message ? e.message : 'error')); }
  }
  const g = voxGeom(nx, ny, nz, data, sp, { warnings, stats }, origin);
  return { name: spec.name || `Procedural ${type}`, format: 'Procedural: ' + type, pathway: 'procedural', kind: g.kind, warnings, stats: g.stats, bbox: g.bbox, voxels: g.voxels };
}
