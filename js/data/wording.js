// Suite-specific wording for the workspace tabs and the run button, so each discipline reads in its
// own vocabulary. Order of the labels: overview, data, model settings, numerical resolution, geometry,
// results, calibration/validation, verification, theory.
const W = (run, guide, inputs, setup, mesh, geometry, results, cal, verify, theory) => ({ run, guide, inputs, setup, mesh, geometry, results, cal, verify, theory });
export const WORDING = {
  default: W('Run simulation', 'Guide', 'Inputs', 'Model setup', 'Mesh', 'Geometry', 'Results', 'Calibrate & validate', 'Verify', 'Theory'),
  pvt: W('Run PVT', 'Overview', 'Fluid & lab data', 'EOS & characterisation', 'P–T grid', 'Data files', 'Phase behaviour', 'Tune & validate', 'Verify', 'Equations'),
  net: W('Build & solve network', 'Overview', 'Route, wells & equipment', 'Hydraulic models', 'Discretisation', 'Geometry import', 'Network & equipment', 'Calibrate & validate', 'Verify', 'Equations'),
  flow: W('Run flow simulation', 'Overview', 'Rates & boundaries', 'Models & conditions', 'Mesh & time step', 'Route & terrain', 'Profiles & slugging', 'Calibrate & validate', 'Verify', 'Equations'),
  solids: W('Run solids model', 'Overview', 'Solids data', 'Kinetics & deposition', 'Grid & classes', 'Deposit geometry', 'Hydrates & deposits', 'Calibrate & validate', 'Verify', 'Equations'),
  ops: W('Run operation', 'Overview', 'Schedules & set points', 'Procedures & control', 'Time step', 'Topology & tags', 'Transient response', 'Replay & validate', 'Verify', 'Equations'),
  integ: W('Run assessment', 'Overview', 'Loads & materials', 'Codes & limit states', 'Mesh & samples', 'Structure & inspection', 'Margins & risk', 'Calibrate & validate', 'Verify', 'Equations'),
  econ: W('Run economics', 'Briefing', 'Costs & prices', 'Fiscal & scenarios', 'Sample size', 'Asset lists & layouts', 'Cash flow & decisions', 'Benchmark & back-cast', 'Audit checks', 'Methods'),
};
export const wording = (id) => WORDING[id] || WORDING.default;
