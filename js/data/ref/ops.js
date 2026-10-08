// Reference data sets for suite 5 (Operations, Control & Flow-Assurance Management).
// Every number below was read from the cited address on the retrieval date; nothing is estimated or "typical".
// The suite (js/suites/s05_ops.js) attaches its own blind prediction to each set; this module holds data only.
const RETRIEVED = '2026-10-08';

// ---- 1. Skogestad (2003): sensitivity peaks of PI / series-PID loops with the SIMC settings (Tables 3 and 4) ----
// plant: k Π(T0 s + 1) e^(-delay s) / (s^ints Π(τ s + 1)); negative T0 = right-half-plane zero (the sign follows from the
// half-rule approximation printed in the same table row). Controller: Kc, τI, τD in the series form Kc (τI s + 1)(τD s + 1)/(τI s);
// ki = pure integral controller. Ms is the published sensitivity peak.
const SIMC_ROWS = [
  { id: 'T3-1', process: 'e^-s', k: 1, lags: [], zeros: [], ints: 0, delay: 1, ki: 0.5, kc: 0, ti: 0, td: 0, Ms: 1.59 },
  { id: 'T3-2', process: 'e^-s / s', k: 1, lags: [], zeros: [], ints: 1, delay: 1, kc: 0.5, ti: 8, td: 0, Ms: 1.70 },
  { id: 'T3-3', process: 'e^-s / (s (4s+1))', k: 1, lags: [4], zeros: [], ints: 1, delay: 1, kc: 0.5, ti: 8, td: 4, Ms: 1.70 },
  { id: 'T3-4', process: 'e^-s / s²', k: 1, lags: [], zeros: [], ints: 2, delay: 1, kc: 0.0625, ti: 8, td: 8, Ms: 1.96 },
  { id: 'T3-5', process: 'e^-s / (4s+1)', k: 1, lags: [4], zeros: [], ints: 0, delay: 1, kc: 2, ti: 4, td: 0, Ms: 1.59 },
  { id: 'E1 PI', process: '1 / ((s+1)(0.2s+1))', k: 1, lags: [1, 0.2], zeros: [], ints: 0, delay: 0, kc: 5.5, ti: 0.8, td: 0, Ms: 1.56 },
  { id: 'E2 PI', process: '(-0.3s+1)(0.08s+1) / ((2s+1)(s+1)(0.4s+1)(0.2s+1)(0.05s+1)³)', k: 1, lags: [2, 1, 0.4, 0.2, 0.05, 0.05, 0.05], zeros: [-0.3, 0.08], ints: 0, delay: 0, kc: 0.85, ti: 2.5, td: 0, Ms: 1.66 },
  { id: 'E2 PID', process: '(-0.3s+1)(0.08s+1) / ((2s+1)(s+1)(0.4s+1)(0.2s+1)(0.05s+1)³)', k: 1, lags: [2, 1, 0.4, 0.2, 0.05, 0.05, 0.05], zeros: [-0.3, 0.08], ints: 0, delay: 0, kc: 1.30, ti: 2, td: 1.2, Ms: 1.73 },
  { id: 'E3 PI', process: '2(15s+1) / ((20s+1)(s+1)(0.1s+1)²)', k: 2, lags: [20, 1, 0.1, 0.1], zeros: [15], ints: 0, delay: 0, kc: 2.33, ti: 1.05, td: 0, Ms: 1.55 },
  { id: 'E3 PID', process: '2(15s+1) / ((20s+1)(s+1)(0.1s+1)²)', k: 2, lags: [20, 1, 0.1, 0.1], zeros: [15], ints: 0, delay: 0, kc: 6.67, ti: 0.4, td: 0.15, Ms: 1.47 },
  { id: 'E4 PI', process: '1 / (s+1)⁴', k: 1, lags: [1, 1, 1, 1], zeros: [], ints: 0, delay: 0, kc: 0.3, ti: 1.5, td: 0, Ms: 1.46 },
  { id: 'E4 PID', process: '1 / (s+1)⁴', k: 1, lags: [1, 1, 1, 1], zeros: [], ints: 0, delay: 0, kc: 0.5, ti: 1.5, td: 1, Ms: 1.43 },
  { id: 'E5 PI', process: '1 / ((s+1)(0.2s+1)(0.04s+1)(0.0008s+1))', k: 1, lags: [1, 0.2, 0.04, 0.0008], zeros: [], ints: 0, delay: 0, kc: 3.72, ti: 1.1, td: 0, Ms: 1.59 },
  { id: 'E5 PID', process: '1 / ((s+1)(0.2s+1)(0.04s+1)(0.0008s+1))', k: 1, lags: [1, 0.2, 0.04, 0.0008], zeros: [], ints: 0, delay: 0, kc: 17.9, ti: 0.224, td: 0.22, Ms: 1.58 },
  { id: 'E6 PI', process: '(0.17s+1)² / (s (s+1)² (0.028s+1))', k: 1, lags: [1, 1, 0.028], zeros: [0.17, 0.17], ints: 1, delay: 0, kc: 0.296, ti: 13.5, td: 0, Ms: 1.48 },
  { id: 'E6 PID', process: '(0.17s+1)² / (s (s+1)² (0.028s+1))', k: 1, lags: [1, 1, 0.028], zeros: [0.17, 0.17], ints: 1, delay: 0, kc: 1.40, ti: 2.86, td: 1.33, Ms: 1.23 },
  { id: 'E7 PI', process: '(-2s+1) / (s+1)³', k: 1, lags: [1, 1, 1], zeros: [-2], ints: 0, delay: 0, kc: 0.214, ti: 1.5, td: 0, Ms: 1.66 },
  { id: 'E7 PID', process: '(-2s+1) / (s+1)³', k: 1, lags: [1, 1, 1], zeros: [-2], ints: 0, delay: 0, kc: 0.3, ti: 1.5, td: 1, Ms: 1.85 },
  { id: 'E8 PI', process: '1 / (s (s+1)²)', k: 1, lags: [1, 1], zeros: [], ints: 1, delay: 0, kc: 0.33, ti: 12, td: 0, Ms: 1.76 },
  { id: 'E8 PID', process: '1 / (s (s+1)²)', k: 1, lags: [1, 1], zeros: [], ints: 1, delay: 0, kc: 1.5, ti: 4, td: 1.5, Ms: 1.79 },
  { id: 'E9 PI', process: 'e^-s / (s+1)²', k: 1, lags: [1, 1], zeros: [], ints: 0, delay: 1, kc: 0.5, ti: 1.5, td: 0, Ms: 1.61 },
  { id: 'E9 PID', process: 'e^-s / (s+1)²', k: 1, lags: [1, 1], zeros: [], ints: 0, delay: 1, kc: 0.5, ti: 1, td: 1, Ms: 1.59 },
  { id: 'E10 PI', process: 'e^-s / ((20s+1)(2s+1))', k: 1, lags: [20, 2], zeros: [], ints: 0, delay: 1, kc: 5.25, ti: 16, td: 0, Ms: 1.72 },
  { id: 'E10 PID', process: 'e^-s / ((20s+1)(2s+1))', k: 1, lags: [20, 2], zeros: [], ints: 0, delay: 1, kc: 10, ti: 8, td: 2, Ms: 1.65 },
  { id: 'E11 PI', process: '(-s+1) e^-s / ((6s+1)(2s+1)²)', k: 1, lags: [6, 2, 2], zeros: [-1], ints: 0, delay: 1, kc: 0.7, ti: 7, td: 0, Ms: 1.63 },
  { id: 'E11 PID', process: '(-s+1) e^-s / ((6s+1)(2s+1)²)', k: 1, lags: [6, 2, 2], zeros: [-1], ints: 0, delay: 1, kc: 1, ti: 6, td: 3, Ms: 1.66 },
  { id: 'E12 PI', process: '(6s+1)(3s+1) e^-0.3s / ((10s+1)(8s+1)(s+1))', k: 1, lags: [10, 8, 1], zeros: [6, 3], ints: 0, delay: 0.3, kc: 7.41, ti: 1, td: 0, Ms: 1.66 },
  { id: 'E13 PI', process: '(2s+1) e^-s / ((10s+1)(0.5s+1))', k: 1, lags: [10, 0.5], zeros: [2], ints: 0, delay: 1, kc: 2.88, ti: 4.5, td: 0, Ms: 1.74 },
  { id: 'E14 PI', process: '(-s+1) / s', k: 1, lags: [], zeros: [-1], ints: 1, delay: 0, kc: 0.5, ti: 8, td: 0, Ms: 2 },
  { id: 'E15 PI', process: '(-s+1) / (s+1)', k: 1, lags: [1], zeros: [-1], ints: 0, delay: 0, kc: 0.5, ti: 1, td: 0, Ms: 2 },
];

// ---- 2. Nitrogen vessel blowdown, experiment I1 of Haque et al. (1992), as digitised in the HydDown validation suite ----
// Vessel: flat-ended vertical cylinder, length 1.524 m, inner diameter 0.273 m, wall 25 mm (steel 7800 kg/m³, 500 J/kg/K);
// nitrogen at 150 bara and 288 K through a 6.35 mm orifice (discharge coefficient 0.8 in the input file) to 1.013 bara.
const N2 = {
  vessel: { length: 1.524, diameter: 0.273, thickness: 0.025, cp: 500, rho: 7800, T0: 288, P0: 150e5, orifice: 0.00635, cd: 0.8, pBack: 101300, tAmb: 288 },
  pt: [0.28869, 5.2776, 10.214, 15.131, 19.77, 24.674, 29.847, 34.747, 39.644, 44.541, 49.436, 54.331, 59.225, 64.119, 69.012, 73.905, 78.798, 83.69, 88.583, 93.475, 98.367],
  p: [150.02, 92.559, 65.72, 50.581, 39.226, 31.656, 25.806, 20.989, 17.548, 14.108, 12.043, 9.9785, 8.2581, 6.5376, 5.5054, 4.4731, 3.7849, 3.0968, 2.4086, 2.0645, 1.7204],
  thT: [0.050285, 5.0799, 10.09, 15.068, 20.011, 24.94, 30.124, 35.029, 39.933, 45.109, 50.004, 54.896, 60.059, 64.947, 70.106, 74.997, 79.885, 85.048, 89.94, 95.102, 99.994],
  th: [288.93, 261.4, 238.25, 222.05, 213.58, 208.2, 206.68, 206.7, 206.73, 207.01, 209.1, 211.96, 215.07, 218.7, 222.59, 225.71, 229.34, 232.46, 235.31, 238.43, 241.29],
  tlT: [0.32393, 5.3957, 10.408, 15.37, 20.328, 25.264, 30.191, 35.381, 40.292, 45.194, 50.096, 55.264, 60.16, 65.054, 70.22, 75.115, 80.279, 85.171, 90.06, 95.223, 100.11],
  tl: [288.67, 251.87, 228.21, 215.62, 203.8, 196.87, 192, 189.2, 187.68, 188.22, 188.76, 190.84, 192.67, 195.01, 197.61, 199.7, 202.56, 205.42, 208.79, 211.91, 215.28],
};

// ---- 3. Pipeline–riser test case: OLGA reference values and the small-scale rig (Jahanshahi & Skogestad) ----
// Parameters of the case (thesis Table 2.1) and of the laboratory rig (Table 2.3); values of the comparison table (paper Table 1 / thesis Table 2.2).
const RISER = {
  olga: { D: 0.12, Dr: 0.1, Lp: 4300, Vp: 48.63, Lr: 300, Lh: 100, Vr: 3.14, theta: Math.PI / 180, rhoL: 832.2, mwG: 0.02, Z: 1, Tp: 337, Tr: 298.3, muL: 1.426e-4, wG: 0.36, wL: 8.64, Ps: 50.1e5, kH: 0.7, kG: 3.49e-2, kL: 2.81e-1, Kpc: 1.16e-2, fric: 'dkm' },
  rig: { D: 0.02, Dr: 0.02, Lp: 69.71, Vp: 0.0219, Lr: 3, Lh: 0.2, Vr: 0.001, theta: (15 * Math.PI) / 180, rhoL: 832.2, mwG: 0.018, Z: 1, Tp: 288, Tr: 288, muL: 1.426e-4, qWaterLmin: 4, qAirLmin: 4.5, PsAtm: 1.01325e5, kH: 1, kG: 2.07e-2, kL: 1.57e-1, Kpc: 2.21e-4, fric: 'dkm' },
  rows: [
    { q: 'zCrit', quantity: 'Critical valve opening (OLGA case)', unit: '%', value: 5 },
    { q: 'period', quantity: 'Oscillation period at the critical opening', unit: 'min', value: 15.6 },
    { q: 'ssPin', quantity: 'Steady-state inlet pressure, valve fully open', unit: 'bar', value: 68.22 },
    { q: 'ssPrb', quantity: 'Steady-state riser-base pressure, valve fully open', unit: 'bar', value: 66.76 },
    { q: 'ssPrt', quantity: 'Steady-state riser-top pressure, valve fully open', unit: 'bar', value: 50.10 },
    { q: 'minPin', quantity: 'Minimum inlet pressure of the slug cycle, valve fully open', unit: 'bar', value: 63.50 },
    { q: 'minPrb', quantity: 'Minimum riser-base pressure of the slug cycle', unit: 'bar', value: 62.08 },
    { q: 'minPrt', quantity: 'Minimum riser-top pressure of the slug cycle', unit: 'bar', value: 50.09 },
    { q: 'minW', quantity: 'Minimum outlet mass rate of the slug cycle', unit: 'kg/s', value: 0.791 },
    { q: 'maxPin', quantity: 'Maximum inlet pressure of the slug cycle, valve fully open', unit: 'bar', value: 75.83 },
    { q: 'maxPrb', quantity: 'Maximum riser-base pressure of the slug cycle', unit: 'bar', value: 74.55 },
    { q: 'maxPrt', quantity: 'Maximum riser-top pressure of the slug cycle', unit: 'bar', value: 50.14 },
    { q: 'maxW', quantity: 'Maximum outlet mass rate of the slug cycle', unit: 'kg/s', value: 31.18 },
    { q: 'zCritRig', quantity: 'Critical valve opening, small-scale rig (experiment)', unit: '%', value: 15 },
  ],
};

// ---- 4. Hock–Schittkowski test problems: best known objective values from the SIF files of the CUTEst collection ----
const HS_ROWS = [
  { problem: 'HS12', n: 2, m: 1, fStar: -30.0 }, { problem: 'HS21', n: 2, m: 1, fStar: -99.96 }, { problem: 'HS22', n: 2, m: 2, fStar: 1.0 }, { problem: 'HS24', n: 2, m: 3, fStar: -1.0 },
  { problem: 'HS29', n: 3, m: 1, fStar: -22.6274169 }, { problem: 'HS30', n: 3, m: 1, fStar: 1.0 }, { problem: 'HS34', n: 3, m: 2, fStar: -0.83403245 }, { problem: 'HS35', n: 3, m: 1, fStar: 0.1111111111 },
  { problem: 'HS39', n: 4, m: 2, fStar: -1.0 }, { problem: 'HS43', n: 4, m: 3, fStar: -44.0 }, { problem: 'HS65', n: 3, m: 1, fStar: 0.9535288567 }, { problem: 'HS71', n: 4, m: 2, fStar: 17.0140173 }, { problem: 'HS100', n: 7, m: 4, fStar: 680.6300573 },
];

// ---- 5. Complementary error function (the travelling-front profile of advection–dispersion is C/C0 = ½ erfc ξ) ----
const ERFC_ROWS = [[0, 1], [0.1, 0.887537084], [0.2, 0.777297411], [0.3, 0.671373241], [0.4, 0.571607645], [0.5, 0.479500122], [0.6, 0.396143909], [0.7, 0.322198806], [0.8, 0.257899035], [0.9, 0.203091788], [1, 0.157299207], [1.1, 0.119794930], [1.2, 0.089686022], [1.3, 0.065992055], [1.4, 0.047714880], [1.5, 0.033894854]].map(([xi, erfc]) => ({ xi, erfc }));

export const REF = Object.freeze({
  retrieved: RETRIEVED,
  simc: {
    id: 'simc-sensitivity-peaks', title: 'Sensitivity peak Ms of 30 PI / PID loops tuned by the SIMC rule', quantity: 'Ms', unit: '–', kind: 'benchmark',
    source: { citation: 'S. Skogestad (2003), Simple analytic rules for model reduction and PID controller tuning, Journal of Process Control 13, 291–309, Tables 3 and 4 (author\'s posted copy)', url: 'https://folk.ntnu.no/skoge/publications/2003/tuningPID/finalpaper.pdf', licence: 'Published numerical results (facts) quoted with citation; paper © Elsevier, copy posted by the author', retrieved: RETRIEVED },
    rows: SIMC_ROWS,
  },
  n2: {
    idP: 'haque-n2-blowdown-pressure', idT: 'haque-n2-blowdown-gas-temperature', vessel: N2.vessel,
    source: { citation: 'M.A. Haque, S.M. Richardson, G. Saville, G. Chamberlain, L. Shirvill (1992), Blowdown of pressure vessels II: experimental validation of computer model and case studies, Trans IChemE 70B, 10–17, experiment I1; points as digitised in the validation suite of HydDown (A. Andreasen, MIT licence), file validation/N2_Exp_I1.yml', url: 'https://raw.githubusercontent.com/andr1976/HydDown/main/validation/N2_Exp_I1.yml', licence: 'MIT (HydDown repository); experimental points are facts from the cited paper', retrieved: RETRIEVED },
    pressure: N2.pt.map((t, i) => ({ t, P: N2.p[i] })),
    temperature: [...N2.thT.map((t, i) => ({ t, sensor: 'upper', T: N2.th[i] })), ...N2.tlT.map((t, i) => ({ t, sensor: 'lower', T: N2.tl[i] }))],
  },
  riser: {
    id: 'jahanshahi-riser-olga', olga: RISER.olga, rig: RISER.rig, rows: RISER.rows,
    source: { citation: 'E. Jahanshahi, S. Skogestad (2011), Simplified dynamical models for control of severe slugging in multiphase risers, 18th IFAC World Congress, 1634–1639, Table 1; parameters from E. Jahanshahi (2013), Control solutions for multiphase flow, PhD thesis, NTNU, Tables 2.1–2.3', url: 'https://folk.ntnu.no/skoge/publications/thesis/2013_jahanshahi/PhDThesis_Jahanshahi_2013.pdf', licence: 'Published numerical results (facts) quoted with citation; thesis and preprint posted openly by NTNU', retrieved: RETRIEVED },
  },
  hs: {
    id: 'hock-schittkowski-sqp', rows: HS_ROWS,
    source: { citation: 'W. Hock, K. Schittkowski (1981), Test examples for nonlinear programming codes, Lecture Notes in Economics and Mathematical Systems 187, Springer; objective values (LO SOLTN) and problem data from the SIF files of the CUTEst test set', url: 'https://bitbucket.org/optrove/sif/raw/HEAD/HS71.SIF', licence: 'CUTEst SIF collection, freely distributed for research and testing; values are mathematical facts', retrieved: RETRIEVED },
  },
  erfc: {
    id: 'advection-dispersion-front', rows: ERFC_ROWS,
    source: { citation: 'Table of values of the error function and its complement, Wikipedia, "Error function" (values agree with Abramowitz & Stegun, Handbook of Mathematical Functions, Table 7.1)', url: 'https://en.wikipedia.org/wiki/Error_function', licence: 'CC BY-SA 4.0', retrieved: RETRIEVED },
  },
  // single published worked example used as a check in verify(): methanol K-value at 1000 psia, 10 °F, 25 wt % methanol
  methanolExample: { psia: 1000, degF: 10, wt: 25, K: 0.00093, x: 0.1579, y: 0.000147, lbPerMMscf: 12.4, url: 'https://jmcampbell.com/tip-of-the-month/2011/08/a-simple-model-for-estimation-of-methanol-loss-to-vapor-phase' },
  // published margins of the SIMC PI settings (τc = θ), Skogestad (2003) Table 2
  simcMargins: { fopdt: { gm: 3.14, pm: 61.4, ms: 1.59, mt: 1.00 }, integrating: { gm: 2.96, pm: 46.9, ms: 1.70, mt: 1.30 } },
});
