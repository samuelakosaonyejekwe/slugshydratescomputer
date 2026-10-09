// Reference data sets of suite 2 (Geometry, Wells, Network & Equipment).
// Every number below was read from the cited address on the retrieval date given with it (measurements and published tables) or is the
// output of the cited open-source library run on that date (kind: 'benchmark'); nothing is typed from memory and nothing is fitted.
// Licences: public domain, CC0, CC BY 4.0 and MIT-licensed software and its output; one manufacturer table is quoted as twelve cited numbers.
// The engine predictions are attached in js/suites/s02_net.js (validationData); this module holds data only.

/** Nikuradse (1933) sand-roughened pipes: friction factor against Reynolds number for six relative roughnesses r/k. */
export const NIKURADSE_ROUGH = {
  id: 'nikuradse-rough-pipe', title: 'Sand-roughened pipes, Darcy friction factor (Nikuradse 1933), r/k = 15 … 507', quantity: 'Darcy friction factor λ', unit: '–', kind: 'experiment',
  source: { citation: 'Nikuradse, J. (1933) Strömungsgesetze in rauhen Rohren, VDI-Forschungsheft 361 (English: Laws of flow in rough pipes, NACA TM 1292, 1950); tabulated data file NikRough_f_vs_Re.xls, sheet rawNikdata, distributed by the Princeton Gas Dynamics Laboratory (A. J. Smits)',
    url: 'https://smits.princeton.edu/wp-content/uploads/sites/356/2016/03/NikRough_f_vs_Re.xls', licence: 'Measured data of a 1933 report republished by NACA (US government, public domain); data file openly distributed by Princeton University', retrieved: '2026-10-08' },
  columns: [{ key: 'rk', label: 'Pipe radius / sand-grain size r/k' }, { key: 'logRe', label: 'log₁₀ Re' }, { key: 'log100f', label: 'log₁₀(100 λ)' }, { key: 'f', label: 'Measured λ' }],
  // seven points per roughness, evenly spaced through each series of the source sheet (after sorting by Reynolds number); logRe and log100f exactly as tabulated
  rows: [
    {"rk":507.0,"logRe":4.114,"log100f":0.456},
    {"rk":507.0,"logRe":4.568,"log100f":0.365},
    {"rk":507.0,"logRe":4.851,"log100f":0.292},
    {"rk":507.0,"logRe":5.143,"log100f":0.26},
    {"rk":507.0,"logRe":5.377,"log100f":0.255},
    {"rk":507.0,"logRe":5.709,"log100f":0.272},
    {"rk":507.0,"logRe":5.987,"log100f":0.286},
    {"rk":252.0,"logRe":4.21,"log100f":0.4506},
    {"rk":252.0,"logRe":4.644,"log100f":0.3434},
    {"rk":252.0,"logRe":4.916,"log100f":0.3222},
    {"rk":252.0,"logRe":5.173,"log100f":0.3197},
    {"rk":252.0,"logRe":5.494,"log100f":0.3504},
    {"rk":252.0,"logRe":5.748,"log100f":0.3608},
    {"rk":252.0,"logRe":5.991,"log100f":0.3716},
    {"rk":126.0,"logRe":3.63,"log100f":0.594},
    {"rk":126.0,"logRe":4.196,"log100f":0.451},
    {"rk":126.0,"logRe":4.76,"log100f":0.4},
    {"rk":126.0,"logRe":5.196,"log100f":0.43},
    {"rk":126.0,"logRe":5.432,"log100f":0.447},
    {"rk":126.0,"logRe":5.736,"log100f":0.452},
    {"rk":126.0,"logRe":5.987,"log100f":0.447},
    {"rk":60.0,"logRe":3.653,"log100f":0.593},
    {"rk":60.0,"logRe":4.236,"log100f":0.487},
    {"rk":60.0,"logRe":4.74,"log100f":0.517},
    {"rk":60.0,"logRe":5.083,"log100f":0.545},
    {"rk":60.0,"logRe":5.408,"log100f":0.55},
    {"rk":60.0,"logRe":5.659,"log100f":0.551},
    {"rk":60.0,"logRe":5.962,"log100f":0.555},
    {"rk":30.6,"logRe":3.672,"log100f":0.592},
    {"rk":30.6,"logRe":3.978,"log100f":0.578},
    {"rk":30.6,"logRe":4.425,"log100f":0.637},
    {"rk":30.6,"logRe":4.9,"log100f":0.656},
    {"rk":30.6,"logRe":5.223,"log100f":0.656},
    {"rk":30.6,"logRe":5.473,"log100f":0.657},
    {"rk":30.6,"logRe":6.0,"log100f":0.659},
    {"rk":15.0,"logRe":3.77,"log100f":0.696},
    {"rk":15.0,"logRe":4.196,"log100f":0.754},
    {"rk":15.0,"logRe":4.5,"log100f":0.777},
    {"rk":15.0,"logRe":4.865,"log100f":0.777},
    {"rk":15.0,"logRe":5.139,"log100f":0.78},
    {"rk":15.0,"logRe":5.515,"log100f":0.781},
    {"rk":15.0,"logRe":6.008,"log100f":0.78},
  ].map((r) => ({ ...r, Re: 10 ** r.logRe, f: 10 ** r.log100f / 100 })),
  target: 'f',
};

/** Smooth-pipe friction factors from the Oregon (liquid helium, gases) and Princeton Superpipe experiments. */
export const SMOOTH_PIPE = {
  id: 'mckeon-smooth-pipe', title: 'Smooth pipe, Darcy friction factor, Re = 11 … 1.05 × 10⁶ (Oregon and Princeton data)', quantity: 'Darcy friction factor', unit: '–', kind: 'experiment',
  source: { citation: 'McKeon, B. J., Swanson, C. J., Zagarola, M. V., Donnelly, R. J. & Smits, A. J. (2004) Friction factors for smooth pipe flow, J. Fluid Mech. 511, 41–44, doi:10.1017/S0022112004009796; values as tabulated in oregon_smooth_data of the open-source fluids library (C. Bell)',
    url: 'https://raw.githubusercontent.com/CalebBell/fluids/master/fluids/friction.py', licence: 'MIT (fluids library); measured values are facts from the cited paper', retrieved: '2026-10-08' },
  columns: [{ key: 'Re', label: 'Reynolds number' }, { key: 'f', label: 'Measured friction factor' }],
  rows: [11.21,20.22,29.28,43.19,57.73,64.58,86.05,113.3,135.3,157.5,179.4,206.4,228.0,270.9,315.2,358.9,402.9,450.2,522.5,583.1,671.8,789.8,891.0,1013.0,1197.0,1300.0,1390.0,1669.0,1994.0,2227.0,2554.0,2868.0,2903.0,2926.0,2955.0,2991.0,2997.0,3047.0,3080.0,3264.0,3980.0,4835.0,5959.0,8162.0,10900.0,13650.0,18990.0,29430.0,40850.0,59220.0,84760.0,120000.0,176000.0,237700.0,298200.0,467800.0,587500.0,824200.0,1050000.0].map((Re, i) => ({ Re, f: [5.537,3.492,2.329,1.523,1.173,0.9863,0.7826,0.5709,0.4815,0.4182,0.3655,0.3237,0.2884,0.2433,0.2077,0.1834,0.1656,0.1475,0.1245,0.1126,0.09917,0.08501,0.07722,0.06707,0.0588,0.05328,0.04815,0.04304,0.03739,0.03405,0.03091,0.02804,0.03182,0.03846,0.03363,0.04124,0.035,0.03875,0.04285,0.0426,0.03995,0.03797,0.0361,0.03364,0.03088,0.02903,0.0267,0.02386,0.02086,0.02,0.01805,0.01686,0.01594,0.01511,0.01462,0.01365,0.01313,0.01244,0.01198][i] })),
  target: 'f',
};

/** Crane fully turbulent friction factor f_T by pipe size as computed by the fluids library (ft_Crane: Colebrook at Re = 7.5 × 10⁶ D with a size-dependent commercial-steel roughness, which reproduces the rounded Crane TP-410 table). */
export const CRANE_FT = {
  id: 'crane-ft', title: 'Fully turbulent friction factor f_T of clean commercial steel pipe by size (fluids ft_Crane)', quantity: 'f_T', unit: '–', kind: 'benchmark',
  source: { citation: 'Bell, C. and contributors, fluids: fluid dynamics component of Chemical Engineering Design Library (ChEDL), version 1.3.1 (Python package), function fluids.friction.ft_Crane(D); output computed on 2026-10-09, seven significant figures kept', url: 'https://pypi.org/project/fluids/1.3.1/', licence: 'MIT (software and its output)', retrieved: '2026-10-09' },
  columns: [{ key: 'd', label: 'Inner diameter', unit: 'm' }, { key: 'fT', label: 'f_T (fluids)' }],
  rows: [
    { d: 0.01576, fT: 0.02556221 },
    { d: 0.02664, fT: 0.02225221 },
    { d: 0.04094, fT: 0.01999874 },
    { d: 0.05248, fT: 0.0188491 },
    { d: 0.07792, fT: 0.01721361 },
    { d: 0.10226, fT: 0.01620925 },
    { d: 0.154, fT: 0.01485458 },
    { d: 0.20274, fT: 0.01403791 },
    { d: 0.25446, fT: 0.01341281 },
    { d: 0.33334, fT: 0.01272276 },
    { d: 0.381, fT: 0.01240077 },
    { d: 0.53994, fT: 0.01161601 },
    { d: 0.8759, fT: 0.01064477 },
  ],
  target: 'fT',
};

/** Field production tests through wellhead chokes (liquid rate against wellhead pressure, gas–liquid ratio and choke size). */
export const CHOKE_FIELD = {
  id: 'choke-field-tests', title: 'Wellhead choke production tests: liquid rate, 36 of 565 field points', quantity: 'Liquid rate', unit: 'STB/d', kind: 'field',
  source: { citation: 'Dabiri, M.-S., Hadavimoghaddam, F., Ashoorian, S., Schaffie, M. & Hemmati-Sarapardeh, A. (2024) Modeling liquid rate through wellhead chokes using machine learning techniques, Scientific Reports 14, 6945, doi:10.1038/s41598-024-54010-2, Supplementary Table S1',
    url: 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41598-024-54010-2/MediaObjects/41598_2024_54010_MOESM1_ESM.docx', licence: 'CC BY 4.0', retrieved: '2026-10-08' },
  columns: [{ key: 'n', label: 'Row number in Table S1' }, { key: 'pwh', label: 'Wellhead pressure', unit: 'psi' }, { key: 'glr', label: 'Gas–liquid ratio', unit: 'scf/STB' }, { key: 'bean', label: 'Choke size', unit: '1/64 in' }, { key: 'q', label: 'Measured liquid rate', unit: 'STB/d' }],
  // every sixteenth row of the 565-row table (first and last included), copied unchanged
  rows: [
    {"n":1,"pwh":784.0,"glr":69.0,"bean":58.88,"q":25878.0},
    {"n":17,"pwh":1484.0,"glr":868.0,"bean":72.0,"q":17829.0},
    {"n":33,"pwh":1430.0,"glr":767.0,"bean":72.0,"q":17290.0},
    {"n":49,"pwh":603.0,"glr":142.0,"bean":64.0,"q":14760.0},
    {"n":65,"pwh":430.0,"glr":119.0,"bean":64.0,"q":11900.0},
    {"n":82,"pwh":2285.0,"glr":1695.2,"bean":50.0,"q":9584.7},
    {"n":98,"pwh":549.0,"glr":188.0,"bean":64.0,"q":11000.0},
    {"n":114,"pwh":2920.0,"glr":2010.9,"bean":48.0,"q":9499.6},
    {"n":130,"pwh":512.0,"glr":201.0,"bean":64.0,"q":9800.0},
    {"n":146,"pwh":502.0,"glr":210.0,"bean":64.0,"q":9320.0},
    {"n":162,"pwh":2073.0,"glr":1665.0,"bean":64.0,"q":11457.3},
    {"n":178,"pwh":543.0,"glr":280.0,"bean":64.0,"q":8300.0},
    {"n":194,"pwh":2385.0,"glr":1800.2,"bean":46.0,"q":7292.3},
    {"n":210,"pwh":1830.0,"glr":1916.4,"bean":68.0,"q":10180.3},
    {"n":227,"pwh":2500.0,"glr":1040.0,"bean":47.0,"q":9610.0},
    {"n":243,"pwh":2279.0,"glr":1646.7,"bean":50.0,"q":7888.6},
    {"n":259,"pwh":910.0,"glr":738.0,"bean":48.0,"q":5012.0},
    {"n":275,"pwh":1450.0,"glr":759.0,"bean":32.0,"q":3812.0},
    {"n":291,"pwh":2970.0,"glr":1964.3,"bean":48.0,"q":8120.1},
    {"n":307,"pwh":2400.0,"glr":1272.5,"bean":31.0,"q":4148.3},
    {"n":323,"pwh":1453.0,"glr":737.0,"bean":32.0,"q":3647.0},
    {"n":339,"pwh":455.0,"glr":364.0,"bean":61.44,"q":5118.0},
    {"n":356,"pwh":1230.0,"glr":816.0,"bean":48.0,"q":5430.0},
    {"n":372,"pwh":420.0,"glr":410.0,"bean":61.44,"q":4330.0},
    {"n":388,"pwh":1330.0,"glr":753.0,"bean":32.0,"q":2907.0},
    {"n":404,"pwh":905.0,"glr":699.0,"bean":24.0,"q":1324.0},
    {"n":420,"pwh":433.0,"glr":661.0,"bean":61.44,"q":3160.0},
    {"n":436,"pwh":2185.0,"glr":3421.0,"bean":40.0,"q":2804.0},
    {"n":452,"pwh":342.0,"glr":67.0,"bean":30.72,"q":2039.0},
    {"n":468,"pwh":311.0,"glr":176.0,"bean":30.72,"q":927.0},
    {"n":484,"pwh":507.0,"glr":390.0,"bean":38.4,"q":1545.0},
    {"n":501,"pwh":1703.0,"glr":828.1,"bean":44.0,"q":4582.3},
    {"n":517,"pwh":1756.0,"glr":833.0,"bean":44.0,"q":4427.9},
    {"n":533,"pwh":2150.0,"glr":1824.3,"bean":64.0,"q":7589.0},
    {"n":549,"pwh":2250.0,"glr":1824.3,"bean":64.0,"q":7047.0},
    {"n":565,"pwh":2825.0,"glr":1820.4,"bean":68.0,"q":10258.0},
  ],
  target: 'q',
};

/** Darcy friction factor from the fluids library (Colebrook equation solved exactly by Clamond’s algorithm; 64/Re when laminar). */
export const BENCH_FRICTION = {
  id: 'fluids-friction', title: 'Code-to-code: Darcy friction factor, Re = 500 … 10⁸, ε/D = 0 … 0.05 (fluids)', quantity: 'Darcy friction factor', unit: '–', kind: 'benchmark',
  source: { citation: 'Bell, C. and contributors, fluids: fluid dynamics component of Chemical Engineering Design Library (ChEDL), version 1.3.1 (Python package), function fluids.friction.friction_factor(Re, eD); output computed on 2026-10-09, seven significant figures kept', url: 'https://pypi.org/project/fluids/1.3.1/', licence: 'MIT (software and its output)', retrieved: '2026-10-09' },
  columns: [{ key: 'Re', label: 'Reynolds number' }, { key: 'eD', label: 'Relative roughness' }, { key: 'f', label: 'Friction factor (fluids)' }],
  rows: [
    { Re: 500, eD: 0, f: 0.128 },
    { Re: 1500, eD: 0.001, f: 0.04266667 },
    { Re: 4000, eD: 0, f: 0.03990701 },
    { Re: 10000, eD: 0, f: 0.03088295 },
    { Re: 10000, eD: 0.001, f: 0.03238181 },
    { Re: 100000, eD: 0, f: 0.01798977 },
    { Re: 100000, eD: 0.0001, f: 0.01851387 },
    { Re: 100000, eD: 0.01, f: 0.03850354 },
    { Re: 1000000, eD: 0, f: 0.01164504 },
    { Re: 1000000, eD: 0.00001, f: 0.01186954 },
    { Re: 1000000, eD: 0.001, f: 0.01994347 },
    { Re: 1000000, eD: 0.05, f: 0.07157375 },
    { Re: 10000000, eD: 0.00001, f: 0.008995712 },
    { Re: 10000000, eD: 0.0001, f: 0.01216608 },
    { Re: 10000000, eD: 0.01, f: 0.03790983 },
    { Re: 100000000, eD: 0, f: 0.005940466 },
    { Re: 100000000, eD: 0.00001, f: 0.008187559 },
    { Re: 100000000, eD: 0.001, f: 0.01963863 },
  ],
  target: 'f',
};

/** Crane TP-410 resistance coefficient of pipe bends of 90° and 180° as coded in the fluids library. */
export const BENCH_BENDS = {
  id: 'fluids-crane-bends', title: 'Code-to-code: Crane bend resistance coefficient, 90° and 180°, r/D = 1.5 … 10 (fluids)', quantity: 'K', unit: '–', kind: 'benchmark',
  source: { citation: 'Bell, C. and contributors, fluids: fluid dynamics component of Chemical Engineering Design Library (ChEDL), version 1.3.1 (Python package), function fluids.fittings.bend_rounded(Di, angle, bend_diameters, method=\'Crane\'); output computed on 2026-10-09, seven significant figures kept', url: 'https://pypi.org/project/fluids/1.3.1/', licence: 'MIT (software and its output)', retrieved: '2026-10-09' },
  columns: [{ key: 'Di', label: 'Inner diameter', unit: 'm' }, { key: 'rD', label: 'Bend radius / diameter' }, { key: 'angle', label: 'Angle', unit: '°' }, { key: 'K', label: 'K (fluids)' }],
  rows: [
    { Di: 0.1023, rD: 1.5, angle: 90, K: 0.2278894 },
    { Di: 0.1023, rD: 3, angle: 90, K: 0.1964592 },
    { Di: 0.1023, rD: 5, angle: 90, K: 0.2485084 },
    { Di: 0.1023, rD: 10, angle: 90, K: 0.4867525 },
    { Di: 0.1023, rD: 1.5, angle: 180, K: 0.3609285 },
    { Di: 0.1023, rD: 5, angle: 180, K: 0.4364107 },
    { Di: 0.2545, rD: 1.5, angle: 90, K: 0.1885838 },
    { Di: 0.2545, rD: 3, angle: 90, K: 0.1625746 },
    { Di: 0.2545, rD: 5, angle: 90, K: 0.2056465 },
    { Di: 0.2545, rD: 10, angle: 90, K: 0.4027991 },
    { Di: 0.2545, rD: 1.5, angle: 180, K: 0.2986768 },
    { Di: 0.2545, rD: 5, angle: 180, K: 0.3611401 },
    { Di: 0.3874, rD: 1.5, angle: 90, K: 0.1738077 },
    { Di: 0.3874, rD: 3, angle: 90, K: 0.1498364 },
    { Di: 0.3874, rD: 5, angle: 90, K: 0.1895335 },
    { Di: 0.3874, rD: 10, angle: 90, K: 0.3712386 },
    { Di: 0.3874, rD: 1.5, angle: 180, K: 0.2752746 },
    { Di: 0.3874, rD: 5, angle: 180, K: 0.3328437 },
  ],
  target: 'K',
};

/** The same function below 90°, where the library applies Crane’s multi-turn formula with n = angle/90 < 1. */
export const BENCH_BENDS_LOW = {
  id: 'fluids-crane-bends-below-90', title: 'Code-to-code: bend resistance coefficient below 90° (fluids, Crane formula continued to n < 1)', quantity: 'K', unit: '–', kind: 'benchmark',
  source: { citation: 'Bell, C. and contributors, fluids: fluid dynamics component of Chemical Engineering Design Library (ChEDL), version 1.3.1 (Python package), function fluids.fittings.bend_rounded(Di, angle, bend_diameters, method=\'Crane\'); output computed on 2026-10-09, seven significant figures kept', url: 'https://pypi.org/project/fluids/1.3.1/', licence: 'MIT (software and its output)', retrieved: '2026-10-09' },
  columns: [{ key: 'Di', label: 'Inner diameter', unit: 'm' }, { key: 'rD', label: 'Bend radius / diameter' }, { key: 'angle', label: 'Angle', unit: '°' }, { key: 'K', label: 'K (fluids)' }],
  rows: [
    { Di: 0.1023, rD: 1.5, angle: 45, K: 0.1613698 },
    { Di: 0.1023, rD: 3, angle: 45, K: 0.12825 },
    { Di: 0.1023, rD: 5, angle: 60, K: 0.1858743 },
    { Di: 0.1023, rD: 10, angle: 30, K: 0.2396375 },
    { Di: 0.2545, rD: 1.5, angle: 45, K: 0.1335373 },
    { Di: 0.2545, rD: 3, angle: 45, K: 0.1061298 },
    { Di: 0.2545, rD: 5, angle: 60, K: 0.1538153 },
    { Di: 0.2545, rD: 10, angle: 30, K: 0.1983056 },
    { Di: 0.3874, rD: 1.5, angle: 45, K: 0.1230743 },
    { Di: 0.3874, rD: 3, angle: 45, K: 0.09781427 },
    { Di: 0.3874, rD: 5, angle: 60, K: 0.1417634 },
    { Di: 0.3874, rD: 10, angle: 30, K: 0.1827678 },
  ],
  target: 'K',
};

/** Crane TP-410 resistance coefficients of full-bore valves as coded in the fluids library. */
export const BENCH_VALVES = {
  id: 'fluids-crane-valves', title: 'Code-to-code: Crane valve resistance coefficients, gate, globe, ball, swing check, butterfly, 2 … 16 in (fluids)', quantity: 'K', unit: '–', kind: 'benchmark',
  source: { citation: 'Bell, C. and contributors, fluids: fluid dynamics component of Chemical Engineering Design Library (ChEDL), version 1.3.1 (Python package), function fluids.fittings.K_gate_valve_Crane, K_globe_valve_Crane, K_ball_valve_Crane, K_swing_check_valve_Crane, K_butterfly_valve_Crane; output computed on 2026-10-09, seven significant figures kept', url: 'https://pypi.org/project/fluids/1.3.1/', licence: 'MIT (software and its output)', retrieved: '2026-10-09' },
  columns: [{ key: 'type', label: 'Valve' }, { key: 'fn', label: 'Library call' }, { key: 'D', label: 'Inner diameter', unit: 'm' }, { key: 'K', label: 'K (fluids)' }],
  rows: [
    { type: "gate", fn: "K_gate_valve_Crane(D, D, 0)", D: 0.0525, K: 0.1507793 },
    { type: "globe", fn: "K_globe_valve_Crane(D, D)", D: 0.0525, K: 6.40812 },
    { type: "ball", fn: "K_ball_valve_Crane(D, D, 0)", D: 0.0525, K: 0.05654224 },
    { type: "check", fn: "K_swing_check_valve_Crane(D, angled=False)", D: 0.0525, K: 0.9423707 },
    { type: "butterfly", fn: "K_butterfly_valve_Crane(D, style=0)", D: 0.0525, K: 0.8481336 },
    { type: "gate", fn: "K_gate_valve_Crane(D, D, 0)", D: 0.1023, K: 0.1296629 },
    { type: "globe", fn: "K_globe_valve_Crane(D, D)", D: 0.1023, K: 5.510674 },
    { type: "ball", fn: "K_ball_valve_Crane(D, D, 0)", D: 0.1023, K: 0.0486236 },
    { type: "check", fn: "K_swing_check_valve_Crane(D, angled=False)", D: 0.1023, K: 0.8103933 },
    { type: "butterfly", fn: "K_butterfly_valve_Crane(D, style=0)", D: 0.1023, K: 0.729354 },
    { type: "gate", fn: "K_gate_valve_Crane(D, D, 0)", D: 0.2027, K: 0.1123077 },
    { type: "globe", fn: "K_globe_valve_Crane(D, D)", D: 0.2027, K: 4.773079 },
    { type: "ball", fn: "K_ball_valve_Crane(D, D, 0)", D: 0.2027, K: 0.0421154 },
    { type: "check", fn: "K_swing_check_valve_Crane(D, angled=False)", D: 0.2027, K: 0.7019234 },
    { type: "butterfly", fn: "K_butterfly_valve_Crane(D, style=0)", D: 0.2027, K: 0.631731 },
    { type: "gate", fn: "K_gate_valve_Crane(D, D, 0)", D: 0.2545, K: 0.1072991 },
    { type: "globe", fn: "K_globe_valve_Crane(D, D)", D: 0.2545, K: 4.560212 },
    { type: "ball", fn: "K_ball_valve_Crane(D, D, 0)", D: 0.2545, K: 0.04023717 },
    { type: "check", fn: "K_swing_check_valve_Crane(D, angled=False)", D: 0.2545, K: 0.6706194 },
    { type: "butterfly", fn: "K_butterfly_valve_Crane(D, style=0)", D: 0.2545, K: 0.4694336 },
    { type: "gate", fn: "K_gate_valve_Crane(D, D, 0)", D: 0.4286, K: 0.09701708 },
    { type: "globe", fn: "K_globe_valve_Crane(D, D)", D: 0.4286, K: 4.123226 },
    { type: "ball", fn: "K_ball_valve_Crane(D, D, 0)", D: 0.4286, K: 0.03638141 },
    { type: "check", fn: "K_swing_check_valve_Crane(D, angled=False)", D: 0.4286, K: 0.6063568 },
    { type: "butterfly", fn: "K_butterfly_valve_Crane(D, style=0)", D: 0.4286, K: 0.3031784 },
  ],
  target: 'K',
};

/** IEC 60534-2-1 liquid sizing by the fluids library (no attached fittings); the first row is the worked example of the standard (Kv = 165). */
export const BENCH_VALVE_LIQUID = {
  id: 'fluids-iec-valve-liquid', title: 'Code-to-code: IEC 60534 control-valve sizing, liquids incl. choked flow (fluids)', quantity: 'Flow coefficient Cv', unit: 'US gpm/psi^½', kind: 'benchmark',
  source: { citation: 'Bell, C. and contributors, fluids: fluid dynamics component of Chemical Engineering Design Library (ChEDL), version 1.3.1 (Python package), function fluids.control_valve.size_control_valve_l(rho, Psat, Pc, mu, P1, P2, Q, FL=FL, Fd=1), Kv_to_Cv; output computed on 2026-10-09, seven significant figures kept', url: 'https://pypi.org/project/fluids/1.3.1/', licence: 'MIT (software and its output)', retrieved: '2026-10-09' },
  columns: [{ key: 'rho', label: 'Density', unit: 'kg/m³' }, { key: 'Psat', label: 'Vapour pressure', unit: 'Pa' }, { key: 'Pc', label: 'Critical pressure', unit: 'Pa' }, { key: 'P1', label: 'Inlet pressure', unit: 'Pa' }, { key: 'P2', label: 'Outlet pressure', unit: 'Pa' }, { key: 'Q', label: 'Volume flow', unit: 'm³/s' }, { key: 'FL', label: 'FL' }, { key: 'choked', label: 'Choked' }, { key: 'Cv', label: 'Cv (fluids)' }],
  rows: [
    { rho: 965.4, Psat: 70100, Pc: 22120000, mu: 0.00031472, P1: 680000, P2: 220000, Q: 0.1, FL: 0.9, choked: false, Kv: 164.9955, Cv: 190.7511 },
    { rho: 965.4, Psat: 70100, Pc: 22120000, mu: 0.00031472, P1: 680000, P2: 600000, Q: 0.1, FL: 0.9, choked: false, Kv: 395.6453, Cv: 457.4052 },
    { rho: 965.4, Psat: 70100, Pc: 22120000, mu: 0.00031472, P1: 680000, P2: 80000, Q: 0.1, FL: 0.9, choked: true, Kv: 158.7054, Cv: 183.4792 },
    { rho: 965.4, Psat: 70100, Pc: 22120000, mu: 0.00031472, P1: 680000, P2: 220000, Q: 0.1, FL: 0.6, choked: true, Kv: 238.0582, Cv: 275.2189 },
    { rho: 998.2, Psat: 2340, Pc: 22064000, mu: 0.001, P1: 1000000, P2: 700000, Q: 0.02, FL: 0.9, choked: false, Kv: 41.55042, Cv: 48.03641 },
    { rho: 998.2, Psat: 2340, Pc: 22064000, mu: 0.001, P1: 5000000, P2: 1000000, Q: 0.05, FL: 0.85, choked: true, Kv: 29.94122, Cv: 34.61502 },
    { rho: 820, Psat: 300000, Pc: 2500000, mu: 0.002, P1: 9000000, P2: 8000000, Q: 0.03, FL: 0.9, choked: false, Kv: 30.94037, Cv: 35.77014 },
    { rho: 820, Psat: 300000, Pc: 2500000, mu: 0.002, P1: 9000000, P2: 4000000, Q: 0.03, FL: 0.9, choked: false, Kv: 13.83695, Cv: 15.99689 },
    { rho: 700, Psat: 1500000, Pc: 3000000, mu: 0.0004, P1: 4000000, P2: 2500000, Q: 0.08, FL: 0.9, choked: false, Kv: 62.24302, Cv: 71.95911 },
    { rho: 1100, Psat: 1000, Pc: 22000000, mu: 0.0015, P1: 300000, P2: 250000, Q: 0.004, FL: 0.9, choked: false, Kv: 21.36823, Cv: 24.7038 },
  ],
  target: 'Cv',
};

/** IEC 60534-2-1 gas sizing by the fluids library (no attached fittings); the first row uses the conditions of the standard’s carbon-dioxide example. */
export const BENCH_VALVE_GAS = {
  id: 'fluids-iec-valve-gas', title: 'Code-to-code: IEC 60534 control-valve sizing, gases incl. choked flow (fluids)', quantity: 'Flow coefficient Cv', unit: 'US gpm/psi^½', kind: 'benchmark',
  source: { citation: 'Bell, C. and contributors, fluids: fluid dynamics component of Chemical Engineering Design Library (ChEDL), version 1.3.1 (Python package), function fluids.control_valve.size_control_valve_g(T, MW, mu, gamma, Z, P1, P2, Q, xT=xT), Kv_to_Cv; output computed on 2026-10-09, seven significant figures kept', url: 'https://pypi.org/project/fluids/1.3.1/', licence: 'MIT (software and its output)', retrieved: '2026-10-09' },
  columns: [{ key: 'T', label: 'Temperature', unit: 'K' }, { key: 'MW', label: 'Molar mass', unit: 'g/mol' }, { key: 'gamma', label: 'Isentropic exponent' }, { key: 'Z', label: 'Compressibility' }, { key: 'P1', label: 'Inlet pressure', unit: 'Pa' }, { key: 'P2', label: 'Outlet pressure', unit: 'Pa' }, { key: 'Q', label: 'Flow at 0 °C, 101.325 kPa', unit: 'm³/s' }, { key: 'xT', label: 'xT' }, { key: 'choked', label: 'Choked' }, { key: 'Cv', label: 'Cv (fluids)' }],
  rows: [
    { T: 433, MW: 44.01, mu: 0.00014665, gamma: 1.3, Z: 0.988, P1: 680000, P2: 310000, Q: 1.055556, xT: 0.6, w: 2.072591, rho1: 8.413588, choked: false, Kv: 62.65206, Cv: 72.432 },
    { T: 433, MW: 44.01, mu: 0.00014665, gamma: 1.3, Z: 0.988, P1: 680000, P2: 600000, Q: 1.055556, xT: 0.6, w: 2.072591, rho1: 8.413588, choked: false, Kv: 97.75636, Cv: 113.0161 },
    { T: 433, MW: 44.01, mu: 0.00014665, gamma: 1.3, Z: 0.988, P1: 680000, P2: 100000, Q: 1.055556, xT: 0.6, w: 2.072591, rho1: 8.413588, choked: true, Kv: 62.63912, Cv: 72.41704 },
    { T: 320, MW: 16.04, mu: 0.000011, gamma: 1.31, Z: 0.998, P1: 130000, P2: 120000, Q: 0.1, xT: 0.7, w: 0.07156251, rho1: 0.7852953, choked: false, Kv: 30.23299, Cv: 34.95233 },
    { T: 300, MW: 19, mu: 0.000012, gamma: 1.28, Z: 0.85, P1: 10000000, P2: 9000000, Q: 30, xT: 0.7, w: 25.43057, rho1: 89.6147, choked: false, Kv: 101.9449, Cv: 117.8584 },
    { T: 300, MW: 19, mu: 0.000012, gamma: 1.28, Z: 0.85, P1: 10000000, P2: 6000000, Q: 30, xT: 0.7, w: 25.43057, rho1: 89.6147, choked: false, Kv: 61.03278, Cv: 70.55996 },
    { T: 300, MW: 19, mu: 0.000012, gamma: 1.28, Z: 0.85, P1: 10000000, P2: 2500000, Q: 30, xT: 0.7, w: 25.43057, rho1: 89.6147, choked: true, Kv: 57.29765, Cv: 66.24177 },
    { T: 350, MW: 28.97, mu: 0.00002, gamma: 1.4, Z: 1, P1: 1000000, P2: 700000, Q: 5, xT: 0.5, w: 6.462488, rho1: 9.955115, choked: false, Kv: 168.1493, Cv: 194.3973 },
    { T: 350, MW: 28.97, mu: 0.00002, gamma: 1.4, Z: 1, P1: 1000000, P2: 200000, Q: 5, xT: 0.5, w: 6.462488, rho1: 9.955115, choked: true, Kv: 156.2975, Cv: 180.6954 },
    { T: 290, MW: 17.4, mu: 0.000011, gamma: 1.27, Z: 0.9, P1: 6000000, P2: 4500000, Q: 12, xT: 0.75, w: 9.315619, rho1: 48.10894, choked: false, Kv: 44.95392, Cv: 51.97119 },
  ],
  target: 'Cv',
};

/** Orifice mass flow with the ISO 5167 expansibility by the fluids library (k = null: liquid, expansibility 1). */
export const BENCH_ORIFICE = {
  id: 'fluids-orifice', title: 'Code-to-code: orifice mass flow, liquids and gases, β = 0.2 … 0.75 (fluids)', quantity: 'Mass flow', unit: 'kg/s', kind: 'benchmark',
  source: { citation: 'Bell, C. and contributors, fluids: fluid dynamics component of Chemical Engineering Design Library (ChEDL), version 1.3.1 (Python package), function fluids.flow_meter.flow_meter_discharge(D, Do, P1, P2, rho, C, expansibility=orifice_expansibility(D, Do, P1, P2, k)); output computed on 2026-10-09, seven significant figures kept', url: 'https://pypi.org/project/fluids/1.3.1/', licence: 'MIT (software and its output)', retrieved: '2026-10-09' },
  columns: [{ key: 'D', label: 'Pipe diameter', unit: 'm' }, { key: 'Do', label: 'Orifice diameter', unit: 'm' }, { key: 'P1', label: 'Upstream pressure', unit: 'Pa' }, { key: 'P2', label: 'Downstream pressure', unit: 'Pa' }, { key: 'rho', label: 'Upstream density', unit: 'kg/m³' }, { key: 'C', label: 'Discharge coefficient' }, { key: 'k', label: 'Isentropic exponent' }, { key: 'm', label: 'Mass flow (fluids)', unit: 'kg/s' }],
  rows: [
    { D: 0.1, Do: 0.05, P1: 200000, P2: 180000, rho: 998, C: 0.61, k: null, eps: 1, m: 7.815726 },
    { D: 0.1, Do: 0.02, P1: 1000000, P2: 500000, rho: 998, C: 0.61, k: null, eps: 1, m: 6.058885 },
    { D: 0.2545, Do: 0.1, P1: 5000000, P2: 4900000, rho: 850, C: 0.6, k: null, eps: 1, m: 62.18768 },
    { D: 0.2545, Do: 0.19, P1: 5000000, P2: 4950000, rho: 850, C: 0.62, k: null, eps: 1, m: 195.1985 },
    { D: 0.05, Do: 0.0127, P1: 10000000, P2: 4000000, rho: 700, C: 0.85, k: null, eps: 1, m: 9.88921 },
    { D: 0.1, Do: 0.05, P1: 200000, P2: 180000, rho: 2.4, C: 0.61, k: 1.4, eps: 0.9731308, m: 0.372976 },
    { D: 0.1, Do: 0.07, P1: 1000000, P2: 800000, rho: 12, C: 0.6, k: 1.3, eps: 0.9264877, m: 5.376743 },
    { D: 0.2545, Do: 0.1, P1: 5000000, P2: 4500000, rho: 45, C: 0.6, k: 1.28, eps: 0.9717414, m: 31.09116 },
    { D: 0.2545, Do: 0.19, P1: 5000000, P2: 4000000, rho: 45, C: 0.62, k: 1.28, eps: 0.9167667, m: 184.1396 },
    { D: 0.0779, Do: 0.03, P1: 2000000, P2: 1600000, rho: 16, C: 0.61, k: 1.31, eps: 0.9440737, m: 1.472661 },
  ],
  target: 'm',
};

/** Compression of an ideal gas with a compressibility factor by the fluids library: polytropic exponent from the polytropic efficiency, polytropic and isentropic work per mole divided by the molar mass, isentropic efficiency and discharge temperature. */
export const BENCH_COMPRESSOR = {
  id: 'fluids-compressor-head', title: 'Code-to-code: compressor polytropic head, ratio 1.5 … 4.4 (fluids)', quantity: 'Polytropic head', unit: 'J/kg', kind: 'benchmark',
  source: { citation: 'Bell, C. and contributors, fluids: fluid dynamics component of Chemical Engineering Design Library (ChEDL), version 1.3.1 (Python package), function fluids.compressible.polytropic_exponent(k, eta_p), isentropic_work_compression(T1, n, Z, P1, P2), isentropic_efficiency(P1, P2, k, eta_p), isentropic_T_rise_compression(T1, P1, P2, k, eta_s); output computed on 2026-10-09, seven significant figures kept', url: 'https://pypi.org/project/fluids/1.3.1/', licence: 'MIT (software and its output)', retrieved: '2026-10-09' },
  columns: [{ key: 'T1', label: 'Suction temperature', unit: 'K' }, { key: 'k', label: 'Isentropic exponent' }, { key: 'Z', label: 'Compressibility' }, { key: 'P1', label: 'Suction pressure', unit: 'Pa' }, { key: 'P2', label: 'Discharge pressure', unit: 'Pa' }, { key: 'etaP', label: 'Polytropic efficiency' }, { key: 'MW', label: 'Molar mass', unit: 'g/mol' }, { key: 'headPoly', label: 'Polytropic head (fluids)', unit: 'J/kg' }, { key: 'headIsen', label: 'Isentropic head (fluids)', unit: 'J/kg' }, { key: 'etaIsen', label: 'Isentropic efficiency (fluids)' }, { key: 'T2', label: 'Discharge temperature (fluids)', unit: 'K' }],
  rows: [
    { T1: 303.15, k: 1.3, Z: 0.95, P1: 1000000, P2: 2500000, etaP: 0.78, MW: 20, n: 1.420168, headPoly: 126012.8, headIsen: 122164.1, etaIsen: 0.7561773, T2: 397.5496 },
    { T1: 303.15, k: 1.3, Z: 0.95, P1: 1000000, P2: 1500000, etaP: 0.78, MW: 20, n: 1.420168, headPoly: 51576.07, headIsen: 50888.01, etaIsen: 0.7695943, T2: 341.787 },
    { T1: 313.15, k: 1.28, Z: 0.9, P1: 2500000, P2: 7000000, etaP: 0.75, MW: 19.5, n: 1.411765, headPoly: 144315.3, headIsen: 138770.4, etaIsen: 0.721183, T2: 422.8375 },
    { T1: 313.15, k: 1.28, Z: 0.9, P1: 2500000, P2: 11000000, etaP: 0.75, MW: 19.5, n: 1.411765, headPoly: 222708.9, headIsen: 210282.7, etaIsen: 0.7081532, T2: 482.4209 },
    { T1: 288.15, k: 1.4, Z: 1, P1: 100000, P2: 300000, etaP: 0.82, MW: 28.97, n: 1.534759, headPoly: 110692.8, headIsen: 106731, etaIsen: 0.7906513, T2: 422.5353 },
    { T1: 288.15, k: 1.4, Z: 1, P1: 100000, P2: 150000, etaP: 0.7, MW: 28.97, n: 1.689655, headPoly: 36466.18, headIsen: 35551.38, etaIsen: 0.6824397, T2: 340.0107 },
    { T1: 320, k: 1.25, Z: 0.85, P1: 6000000, P2: 15000000, etaP: 0.8, MW: 22, n: 1.333333, headPoly: 105853.5, headIsen: 103374.9, etaIsen: 0.7812677, T2: 402.3787 },
    { T1: 300, k: 1.31, Z: 0.98, P1: 500000, P2: 1750000, etaP: 0.76, MW: 16.04, n: 1.452159, headPoly: 233506.6, headIsen: 222233.6, etaIsen: 0.7233094, T2: 443.1265 },
    { T1: 330, k: 1.15, Z: 0.8, P1: 800000, P2: 2000000, etaP: 0.72, MW: 44.1, n: 1.221239, headPoly: 49610.84, headIsen: 48444.43, etaIsen: 0.7030719, T2: 389.5871 },
    { T1: 298.15, k: 1.67, Z: 1, P1: 200000, P2: 500000, etaP: 0.85, MW: 4.003, n: 1.893929, headPoly: 709920, headIsen: 685781.4, etaIsen: 0.8210984, T2: 459.4745 },
  ],
  target: 'headPoly',
};

/** Souders–Brown coefficient of a vessel with a mesh mist eliminator against pressure by the fluids library (York curve). */
export const BENCH_DEMISTER = {
  id: 'fluids-demister-k', title: 'Code-to-code: mist-eliminator Souders–Brown K, 1.5 … 150 bar (fluids)', quantity: 'K', unit: 'm/s', kind: 'benchmark',
  source: { citation: 'Bell, C. and contributors, fluids: fluid dynamics component of Chemical Engineering Design Library (ChEDL), version 1.3.1 (Python package), function fluids.separator.K_separator_demister_York(P); output computed on 2026-10-09, seven significant figures kept', url: 'https://pypi.org/project/fluids/1.3.1/', licence: 'MIT (software and its output)', retrieved: '2026-10-09' },
  columns: [{ key: 'P', label: 'Pressure', unit: 'Pa' }, { key: 'K', label: 'K (fluids)', unit: 'm/s' }],
  rows: [
    { P: 150000, K: 0.10668 },
    { P: 500000, K: 0.1010325 },
    { P: 1000000, K: 0.09617328 },
    { P: 2500000, K: 0.08974972 },
    { P: 4000000, K: 0.0864548 },
    { P: 6720000, K: 0.08281785 },
    { P: 10000000, K: 0.08003124 },
    { P: 15000000, K: 0.07718877 },
  ],
  target: 'K',
};

/** Beij (1938): measured loss coefficients of nine 90° bends in 4-inch steel pipe (water, Re about 4 × 10⁴ … 4 × 10⁵). eta is the loss of the bend in excess of the same axial length of straight pipe, theta the excess loss in the downstream tangent. */
export const BEIJ_BENDS = {
  id: 'beij-90-degree-bends', title: 'Measured bend-loss coefficients of nine 90° bends, 4-inch steel pipe, R/d = 0.97 … 19.96 (Beij 1938)', quantity: 'Bend coefficient η', unit: '–', kind: 'experiment',
  source: { citation: 'Beij, K. H. (1938) Pressure losses for fluid flow in 90° pipe bends, Journal of Research of the National Bureau of Standards 21(1), 1–18, Research Paper RP1110, Table 1 (p. 11)', url: 'https://nvlpubs.nist.gov/nistpubs/jres/21/jresv21n1p1_A1b.pdf', licence: 'public domain (work of the US National Bureau of Standards)', retrieved: '2026-10-09' },
  columns: [{ key: 'bend', label: 'Bend number' }, { key: 'rD', label: 'R/d' }, { key: 'eta', label: 'Measured bend coefficient η' }, { key: 'theta', label: 'Measured tangent coefficient θ' }],
  rows: [
    { bend: 1, rD: 0.97, eta: 0.36, theta: 0.18 },
    { bend: 2, rD: 1.47, eta: 0.214, theta: 0.18 },
    { bend: 3, rD: 3.35, eta: 0.145, theta: 0.18 },
    { bend: 4, rD: 4.97, eta: 0.174, theta: 0.18 },
    { bend: 5, rD: 7.97, eta: 0.27, theta: 0.18 },
    { bend: 6, rD: 11.93, eta: 0.347, theta: 0.09 },
    { bend: 7, rD: 4.04, eta: 0.178, theta: 0.18 },
    { bend: 8, rD: 9.93, eta: 0.165, theta: 0.09 },
    { bend: 9, rD: 19.96, eta: 0.411, theta: 0.18 },
  ],
  target: 'eta',
  pipeId: 0.1023, // m, actual inside diameter
};

/** Manufacturer flow coefficients of a 3-inch cage-guided globe control valve with an equal-percentage trim against travel (twelve numbers quoted: ten Cv values and FL = 0.85, xT = 0.70 at full travel). */
export const VALVE_CV_TRAVEL = {
  id: 'valve-cv-travel', title: 'Control-valve Cv against travel, DN 80 cage-guided globe valve, equal-percentage trim (manufacturer table)', quantity: 'Cv', unit: 'US gpm/psi^½', kind: 'experiment',
  source: { citation: 'Valmet Flow Control (Neles), “Cv-values: GBAE” — GB series globe control valve, balanced standard trim, equal percentage; sizing-coefficient table dated 13.1.2025, row DN 80 / 3 in / reduction number 0', url: 'https://www.valmet.com/globalassets/flow-control/services/GB_GBAE.html', licence: 'manufacturer data, all rights reserved: twelve individual numbers quoted as cited facts', retrieved: '2026-10-09' },
  columns: [{ key: 'travel', label: 'Travel', unit: '%' }, { key: 'cv', label: 'Published Cv' }],
  rows: [
    { travel: 10, cv: 4 },
    { travel: 20, cv: 8.5 },
    { travel: 30, cv: 14.5 },
    { travel: 40, cv: 22 },
    { travel: 50, cv: 33.5 },
    { travel: 60, cv: 51 },
    { travel: 70, cv: 75 },
    { travel: 80, cv: 108 },
    { travel: 90, cv: 128 },
    { travel: 100, cv: 138 },
  ],
  target: 'cv',
  ratedCv: 138, FL: 0.85, xT: 0.7, // at 100 % travel; the page states no rangeability
};

/** Centrifugal compressor map of 55 points on six speed lines (speed, mass flow, pressure ratio, efficiency; the file states no units and no geometry — a map from a design manual, not a documented test). The 740 line is the comparison target; the other five lines are the map handed to the interpolation. */
export const COMPRESSOR_MAP = {
  id: 'compressor-map-speed-line', title: 'Compressor map: speed of the held-out 740 line from the five other speed lines (open map data)', quantity: 'Speed', unit: 'as tabulated', kind: 'experiment',
  source: { citation: 'Li, X., Yang, C., Wang, Y., Wang, H., Zu, X. & Sun, Y. (2018) Data from: Compressor map regression modelling based on partial least squares, Dryad, doi:10.5061/dryad.3g0p68m (file “original data.xlsx”, sheet 1); article: Royal Society Open Science 5, 172454', url: 'https://zenodo.org/api/records/4993647/files/original%20data.xlsx/content', licence: 'CC0 1.0', retrieved: '2026-10-09' },
  columns: [{ key: 'n', label: 'Speed' }, { key: 'm', label: 'Mass flow' }, { key: 'pr', label: 'Pressure ratio' }, { key: 'eta', label: 'Efficiency', unit: '%' }],
  rows: [
    { n: 740, m: 2.4443, pr: 2.8774, eta: 63.74 },
    { n: 740, m: 2.4441, pr: 3.0653, eta: 67.883 },
    { n: 740, m: 2.4347, pr: 3.2529, eta: 72.009 },
    { n: 740, m: 2.3884, pr: 3.4375, eta: 75.535 },
    { n: 740, m: 2.3455, pr: 3.5284, eta: 77.061 },
    { n: 740, m: 2.2791, pr: 3.6176, eta: 78.372 },
    { n: 740, m: 2.2346, pr: 3.6615, eta: 78.91 },
    { n: 740, m: 2.2007, pr: 3.6879, eta: 79.194 },
    { n: 740, m: 2.1722, pr: 3.7056, eta: 79.354 },
    { n: 740, m: 2.0557, pr: 3.7494, eta: 79.531 },
    { n: 740, m: 2.0277, pr: 3.7534, eta: 79.477 },
  ],
  target: 'n',
  map: [ { n: 500, m: 1.4292, pr: 1.5347, eta: 52.253 }, { n: 500, m: 1.407, pr: 1.7201, eta: 66.968 }, { n: 500, m: 1.2501, pr: 1.9033, eta: 78.574 }, { n: 500, m: 1.137, pr: 1.9461, eta: 80.209 }, { n: 500, m: 1.0664, pr: 1.9629, eta: 80.347 }, { n: 500, m: 1.0196, pr: 1.9706, eta: 80.171 }, { n: 600, m: 1.8012, pr: 1.9528, eta: 58.014 }, { n: 600, m: 1.8009, pr: 2.1379, eta: 66.767 }, { n: 600, m: 1.7907, pr: 2.2329, eta: 70.847 }, { n: 600, m: 1.7581, pr: 2.3257, eta: 74.572 }, { n: 600, m: 1.6932, pr: 2.4159, eta: 77.734 }, { n: 600, m: 1.5685, pr: 2.5039, eta: 80.127 }, { n: 600, m: 1.4008, pr: 2.5444, eta: 80.301 }, { n: 700, m: 2.2575, pr: 2.7597, eta: 67.65 }, { n: 700, m: 2.2568, pr: 2.8536, eta: 70.21 }, { n: 700, m: 2.2477, pr: 2.9437, eta: 72.629 }, { n: 700, m: 2.225, pr: 3.0397, eta: 74.837 }, { n: 700, m: 2.1858, pr: 3.1309, eta: 76.826 }, { n: 700, m: 2.1228, pr: 3.2205, eta: 78.548 }, { n: 700, m: 2.0076, pr: 3.3072, eta: 79.784 }, { n: 700, m: 1.8626, pr: 3.3486, eta: 79.908 }, { n: 700, m: 1.7395, pr: 3.3533, eta: 79.261 }, { n: 800, m: 2.6836, pr: 3.4748, eta: 65.923 }, { n: 800, m: 2.6786, pr: 3.6626, eta: 69.211 }, { n: 800, m: 2.6599, pr: 3.8488, eta: 72.154 }, { n: 800, m: 2.6447, pr: 3.9414, eta: 73.453 }, { n: 800, m: 2.624, pr: 4.0335, eta: 74.638 }, { n: 800, m: 2.596, pr: 4.1252, eta: 75.721 }, { n: 800, m: 2.5621, pr: 4.2172, eta: 76.716 }, { n: 800, m: 2.487, pr: 4.3097, eta: 77.559 }, { n: 800, m: 2.3307, pr: 4.401, eta: 78.158 }, { n: 800, m: 2.2452, pr: 4.4274, eta: 78.157 }, { n: 840, m: 2.788, pr: 3.6882, eta: 61.79 }, { n: 840, m: 2.7873, pr: 3.8741, eta: 64.674 }, { n: 840, m: 2.7827, pr: 4.0593, eta: 67.949 }, { n: 840, m: 2.7769, pr: 4.1508, eta: 69.589 }, { n: 840, m: 2.7685, pr: 4.2429, eta: 71.023 }, { n: 840, m: 2.7582, pr: 4.3352, eta: 72.234 }, { n: 840, m: 2.7458, pr: 4.4281, eta: 73.278 }, { n: 840, m: 2.7288, pr: 4.5216, eta: 74.2 }, { n: 840, m: 2.7034, pr: 4.6166, eta: 75.053 }, { n: 840, m: 2.6654, pr: 4.7119, eta: 75.833 }, { n: 840, m: 2.6029, pr: 4.8067, eta: 76.523 }, { n: 840, m: 2.401, pr: 4.897, eta: 76.731 }, ],
};

/** Seven gas wells of the Yan’an field: binomial deliverability coefficients A and B (p_R² − p_wf² = A q + B q², MPa² and 10⁴ m³/d) from isochronal tests (Table 3), the stabilised test point of each well (Table 4) and the absolute open-flow potential listed with the coefficients (Table 3; reservoir pressure there from the tabulated Δp² at 0.101 MPa). */
export const GAS_WELL_TESTS = {
  id: 'gas-well-deliverability', title: 'Gas-well tests: stabilised rates and open-flow potentials of seven wells from their binomial coefficients (Yan’an field)', quantity: 'Gas rate', unit: '10⁴ m³/d', kind: 'field',
  source: { citation: 'Liu, E., Liu, Y., Gao, L., Zhou, D., Liu, X. & Xu, J. (2021) On the one-point model for the productivity evaluation in Jingbian sector of Yan’an gas field, Frontiers in Earth Science 9, 793293, doi:10.3389/feart.2021.793293, Tables 3 and 4', url: 'https://www.frontiersin.org/journals/earth-science/articles/10.3389/feart.2021.793293/xml', licence: 'CC BY 4.0', retrieved: '2026-10-09' },
  columns: [{ key: 'well', label: 'Well' }, { key: 'point', label: 'Point' }, { key: 'pRes', label: 'Reservoir pressure', unit: 'MPa' }, { key: 'pwf', label: 'Flowing bottom-hole pressure', unit: 'MPa' }, { key: 'A', label: 'Laminar coefficient A' }, { key: 'B', label: 'Turbulence coefficient B' }, { key: 'q', label: 'Reported gas rate', unit: '10⁴ m³/d' }],
  rows: [
    { well: "J32-1", point: "stabilised test", pRes: 33.38, pwf: 29.407, A: 34.4174, B: 0.1619, q: 6.2109 },
    { well: "J44", point: "stabilised test", pRes: 30.858, pwf: 8.831, A: 184.0653, B: 1.7577, q: 3.857 },
    { well: "J53-1", point: "stabilised test", pRes: 31.972, pwf: 13.296, A: 273.1753, B: 1.588, q: 2.3937 },
    { well: "Y924", point: "stabilised test", pRes: 33.159, pwf: 30.103, A: 52.038, B: 0.7448, q: 3.2295 },
    { well: "Y942-3", point: "stabilised test", pRes: 32.564, pwf: 26.438, A: 87.4788, B: 0.4263, q: 3.7263 },
    { well: "Y865", point: "stabilised test", pRes: 33.184, pwf: 24.315, A: 46.9808, B: 0.4088, q: 9.0216 },
    { well: "Y313-1", point: "stabilised test", pRes: 31.039, pwf: 24.55, A: 64.7204, B: 0.4913, q: 4.5254 },
    { well: "J32-1", point: "absolute open flow", pRes: 33.38, pwf: 0.101, A: 34.4174, B: 0.1619, q: 28.5416 },
    { well: "J44", point: "absolute open flow", pRes: 30.858, pwf: 0.101, A: 184.0653, B: 1.7577, q: 4.9401 },
    { well: "J53-1", point: "absolute open flow", pRes: 31.872, pwf: 0.101, A: 273.1753, B: 1.588, q: 3.6415 },
    { well: "Y924", point: "absolute open flow", pRes: 33.159, pwf: 0.101, A: 52.038, B: 0.7448, q: 16.995 },
    { well: "Y942-3", point: "absolute open flow", pRes: 32.564, pwf: 0.101, A: 87.4788, B: 0.4263, q: 11.4796 },
    { well: "Y865", point: "absolute open flow", pRes: 33.184, pwf: 0.101, A: 46.9808, B: 0.4088, q: 19.9688 },
    { well: "Y313-1", point: "absolute open flow", pRes: 31.039, pwf: 0.101, A: 64.7204, B: 0.4913, q: 13.5019 },
  ],
  target: 'q',
};

/** Deviation survey of geothermal well Utah FORGE 16A(78)-32 (build-and-hold to 65–69°): all 422 accepted stations with the positions reported by the survey contractor (minimum curvature, grid north, feet below the kelly bushing). */
export const SURVEY_FORGE = {
  id: 'deviation-survey-forge-16a', title: 'Deviation survey of Utah FORGE well 16A(78)-32: true vertical depth at 24 of 422 stations', quantity: 'True vertical depth', unit: 'ft', kind: 'field',
  source: { citation: 'Utah FORGE: Well 16A(78)-32 drilling data, Geothermal Data Repository submission 1283 (University of Utah / US Department of Energy), file “16A(78)-32 Survey.xlsx” (survey report: 422 accepted surveys with measured depth, inclination, grid azimuth, TVD, northing and easting offsets)',
    url: 'https://gdr.openei.org/files/1283/16A(78)-32%20Survey.xlsx', licence: 'CC BY 4.0', retrieved: '2026-10-09' },
  columns: [{ key: 'md', label: 'Measured depth', unit: 'ft' }, { key: 'inc', label: 'Inclination', unit: '°' }, { key: 'azi', label: 'Grid azimuth', unit: '°' }, { key: 'tvd', label: 'Reported TVD', unit: 'ft' }, { key: 'north', label: 'Reported offset north', unit: 'ft' }, { key: 'east', label: 'Reported offset east', unit: 'ft' }],
  // all 422 stations (needed to integrate the path) and 24 evenly spaced stations used for the comparison, copied unchanged (the source prints two decimals)
  survey: { md: [0,80,175,267,359,453,548,643,739,833,928,1022,1117,1211,1306,1402,1497,1592,1687,1783,1878,1974,2069,2164,2258,2352,2446,2540,2634,2729,2823,2918,3013,3107,3201,3295,3390,3484,3578,3673,3766,3862,3956,4051,4146,4241,4335,4430,4524,4618,4713,4807,4902,4997,5007,5017,5027,5037,5047,5057,5067,5077,5087,5097,5107,5117,5127,5137,5147,5157,5167,5177,5187,5197,5207,5217,5227,5237,5247,5257,5267,5277,5287,5297,5307,5317,5327,5337,5347,5357,5367,5377,5387,5397,5407,5417,5427,5437,5447,5457,5467,5477,5487,5497,5507,5517,5527,5537,5547,5557,5567,5577,5587,5597,5607,5617,5627,5637,5647,5657,5667,5677,5687,5697,5707,5717,5727,5737,5747,5757,5767,5777,5787,5797,5807,5817,5827,5837,5847,5857,5867,5877,5887,5897,5907,5917,5927,5937,5947,5957,5967,5977,5987,5997,6007,6017,6027,6037,6047,6057,6067,6077,6087,6097,6107,6117,6127,6137,6147,6157,6167,6177,6187,6197,6207,6217,6227,6237,6247,6257,6267,6277,6287,6297,6307,6317,6327,6337,6347,6357,6367,6377,6387,6397,6407,6417,6427,6437,6447,6457,6467,6477,6487,6497,6507,6517,6527,6537,6547,6557,6567,6577,6587,6597,6607,6617,6627,6637,6647,6657,6667,6677,6687,6697,6707,6717,6727,6737,6747,6757,6767,6777,6787,6797,6807,6817,6827,6837,6847,6857,6867,6877,6887,6897,6907,6917,6927,6937,6947,6957,6967,6977,6987,6997,7007,7017,7027,7037,7047,7057,7067,7077,7087,7097,7107,7117,7127,7137,7147,7157,7167,7177,7187,7197,7207,7217,7227,7237,7247,7257,7267,7277,7287,7297,7307,7317,7327,7337,7347,7357,7367,7377,7387,7397,7407,7417,7427,7437,7447,7457,7467,7477,7487,7497,7507,7517,7527,7537,7547,7557,7567,7577,7587,7597,7607,7617,7627,7637,7647,7657,7667,7677,7687,7697,7707,7717,7727,7737,7747,7757,7767,7777,7787,7797,7807,7817,7827,7837,7847,7857,7867,7877,7887,7897,7907,7917,7927,7937,7947,7957,7967,7977,7987,7997,8007,8017,8027,8037,8047,8057,8067,8077,8087,8097,8107,8117,8127,8137,8147,8157,8167,8177,8187,8197,8207,8217,8227,8237,8247,8257,8267,8277,8287,8297,8307,8317,8327,8337,8347,8357,8367,8377,8387,8397,8406.1,8466,8573,8668,8763,8859,8953,9048,9144,9239,9333,9429,9524,9618,9741,9837,9933,10028,10124,10220,10316,10411,10505,10604,10694,10790,10886,10955], inc: [0,0.34,0.33,0.47,0.36,0.03,0.56,0.59,0.51,0.84,1.24,1.62,1.93,2.07,1.68,1.19,0.7,0.72,0.82,0.87,1,1.14,1.16,1.04,0.87,0.71,0.74,0.76,1.13,1.21,1.01,0.63,1.04,1.38,1.54,1.3,0.99,0.81,0.48,0.75,1.46,1.71,1.66,1.35,1.12,0.84,1.07,1.11,1.45,1.64,1.45,1.13,1.21,1.6,1.57,1.54,1.51,1.49,1.48,1.46,1.46,1.45,1.46,1.46,1.47,1.49,1.5,1.53,1.55,1.58,1.58,1.56,1.55,1.54,1.53,1.52,1.51,1.51,1.5,1.5,1.49,1.48,1.48,1.49,1.5,1.52,1.54,1.56,1.59,1.53,1.46,1.4,1.33,1.28,1.23,1.18,1.14,1.11,1.03,0.9,0.76,0.62,0.49,0.35,0.21,0.08,0.06,0.19,0.27,0.38,0.5,0.62,0.75,0.87,1,1.12,1.25,1.36,1.44,1.52,1.61,1.69,1.78,1.87,1.95,2.04,2.11,2.07,2.04,2,1.96,1.93,1.9,1.86,1.83,1.79,1.86,1.99,2.12,2.26,2.4,2.54,2.68,2.83,2.98,3.22,3.83,4.44,5.05,5.67,6.29,6.9,7.52,8.14,8.75,9.37,9.99,10.6,11.22,11.84,12.46,13.07,13.69,14.31,14.85,15.13,15.43,15.73,16.04,16.36,16.69,17.03,17.37,17.72,17.86,17.99,18.13,18.26,18.38,18.39,18.39,18.29,18.21,18.14,17.95,17.7,17.59,17.45,17.39,17.38,17.32,17.55,17.87,18.84,19.62,20.26,20.69,20.94,21.1,21.45,21.94,22.43,22.82,23.07,23.29,23.56,24.08,24.83,25.53,26.28,27.19,28.16,28.6,28.74,28.98,29.48,30.08,30.88,31.63,32.41,33.19,34.4,35.07,35.42,35.81,36.14,36.58,37.15,37.72,38.24,39.14,39.56,39.85,40.2,40.35,40.87,41.59,42.46,43.24,44.22,45.11,45.59,45.94,46.42,46.97,47.58,48.11,48.66,49.11,49.78,50.57,51.25,52,52.69,53.14,54.01,54.65,55.19,55.93,56.46,56.92,57.3,58.09,58.37,58.62,58.79,58.85,59.18,59.83,60.49,61.16,61.71,62.07,62.08,62.64,63.42,63.66,64.27,64.67,64.94,65.18,65.19,65.36,65.84,66.12,66.28,66.6,66.77,66.89,67.06,66.46,67.49,66.81,66.16,65.23,64.75,64.15,64.43,64.66,65.22,65.75,66.18,66.15,65.6,65.17,65.72,65.71,64.48,65.76,65.63,66.17,66.87,66.29,65.05,65.18,65.42,64.36,64.46,64.1,63.14,62.36,61.19,60.25,60.16,59.67,58.31,59.08,57.74,57.54,56.93,57.3,56.69,56.55,56.88,56.27,56.84,55.56,56.05,55.7,55.44,56.5,56.72,56.33,57.31,57.39,57.7,57.09,57.82,57.78,58.69,58.35,58.05,59.55,59.63,60.18,59.61,60.79,61.73,61.68,62.77,63.65,63.47,63.23,64.29,64.56,64.96,64.51,65.4,65.45,65.39,64.46,64.5,65.25,64.59,65.19,64.26,65.61,65.27,65.35,65.99,66.53,66.37,66.65,66.08,66.72,66.25,66.07,66.3,66.31,66.11,66.69,65.39,66.37,66.78,66.17,65.6,65.9,67.92,68.32,69.41,65.25,66.34,66.13,66.64,64.64,63.62,62.73,62.36,62.89,65.42,67.55,65.93,65.17,64.32,63.32,64.68,64.86,68.46,69.27,68,68.6,68.6], azi: [0,324.91,322.79,312.36,310.55,309.87,157.85,172.48,198.76,191.72,186.85,188.31,182.52,175.74,171.17,173.45,142.32,125.01,118.93,95.89,65.55,14.63,321.31,278.62,236.76,163.66,107.66,48.82,359.64,312.16,262.95,188.09,123.23,79.63,28.28,349.04,291.49,250.73,184.81,95.07,59.45,33.75,3.11,327.79,281.82,230.53,183.14,116.98,32.77,341.93,284.79,233.6,185.35,134.6,131.93,129.15,126.28,123.33,120.3,117.22,114.09,110.94,107.78,104.64,101.54,98.48,95.49,92.58,89.76,87.04,85.49,84.02,82.52,81.01,79.47,77.91,76.33,74.74,73.14,71.38,68.16,64.92,61.67,58.44,55.24,52.1,49.03,46.04,43.16,40.89,38.32,35.51,32.44,29.08,25.43,21.46,17.2,12.65,9.93,10.09,10.32,10.64,11.14,12.02,14.04,22.99,170.93,181.76,157.9,146.55,140.39,136.61,134.07,132.26,130.9,129.84,129,127.86,125.27,122.96,120.88,119.01,117.32,115.78,114.39,113.12,112.36,112.92,113.5,114.1,114.72,115.37,116.03,116.73,117.45,118.19,116.32,113.39,110.82,108.55,106.55,104.77,103.18,101.75,100.47,99.27,97.79,96.72,95.91,95.27,94.76,94.34,93.99,93.69,93.43,93.21,93.01,92.84,92.68,92.54,92.42,92.3,92.2,92.11,92.48,94.05,95.56,97.02,98.42,99.77,101.07,102.32,103.52,104.68,104.94,105.1,105.27,105.44,105.62,105.91,106.2,106.44,106.73,106.96,107.64,108.18,108.98,109.58,110.07,110.41,110.72,111,111.26,111.6,111.9,112.17,112.57,112.89,113.11,113.28,113.43,113.57,113.85,114.51,115.06,115.51,115.61,115.58,115.62,115.65,115.58,115.37,115.2,115.11,114.98,114.72,114.42,114.12,113.78,113.19,112.51,111.33,110.58,109.97,109.5,109.05,108.38,107.37,106.36,105.52,104.36,103.1,102.29,101.75,101.44,101.13,100.87,100.59,100.39,100.23,100.02,100.03,100.12,100.11,100,99.93,99.88,99.89,99.97,100.06,100.25,100.41,100.57,100.72,100.86,100.93,101.02,101.15,101.38,101.78,102.04,102.22,102.39,102.47,102.55,102.58,102.58,102.6,102.58,102.54,102.49,102.45,102.41,102.25,101.98,101.76,101.49,101.43,101.36,101.41,101.42,101.43,101.42,101.36,101.32,101.3,101.3,101.29,101.29,101.32,101.55,101.64,101.86,102.17,102.39,102.57,102.63,102.69,102.72,102.75,102.8,102.85,103.07,103.43,103.88,104.35,104.71,105.01,105.24,105.45,105.64,105.77,105.95,106.17,106.32,106.55,106.78,106.88,107.09,107.76,108.2,108.52,108.65,108.88,109.11,109.2,109.31,109.39,109.51,109.63,109.74,109.9,110.04,110.15,110.21,110.22,110.23,110.26,110.17,110.13,110.1,110,109.9,109.81,109.66,109.57,109.46,109.45,109.45,109.48,109.52,109.53,109.52,109.46,109.46,109.42,109.42,109.39,109.34,109.37,109.37,109.39,109.44,109.57,109.61,109.64,109.72,109.88,109.94,110.02,110.12,110.2,110.26,110.36,110.44,110.51,110.53,110.45,110.31,109.92,109.62,109.34,109.26,109.26,109.23,109.27,109.2,109.09,108.82,108.5,108.25,108.05,107.96,107.99,108.02,108.41,106.23,103.08,100.55,99.11,98.06,98,99.91,103.7,103.67,105.78,106.61,105.28,106.07,105.64,105.62,105.44,106.84,108.36,108.3,105.97,105.02,103.44,101.13,101.18,101.44,101.44] },
  rows: [
    { md: 1687, inc: 0.82, azi: 118.93, tvd: 1686.72, north: -20.15, east: 0.68 },
    { md: 3295, inc: 1.3, azi: 349.04, tvd: 3294.48, north: -13.08, east: 4.42 },
    { md: 4997, inc: 1.6, azi: 134.6, tvd: 4996.13, north: -6.5, east: 1.25 },
    { md: 5167, inc: 1.58, azi: 85.49, tvd: 5166.07, north: -7.93, east: 5.31 },
    { md: 5347, inc: 1.59, azi: 43.16, tvd: 5346, north: -6.18, east: 9.63 },
    { md: 5517, inc: 0.08, azi: 22.99, tvd: 5515.97, north: -3.57, east: 10.84 },
    { md: 5697, inc: 1.87, azi: 115.78, tvd: 5695.94, north: -5.35, east: 13.23 },
    { md: 5867, inc: 2.4, azi: 106.55, tvd: 5865.84, north: -7.74, east: 18.62 },
    { md: 6047, inc: 11.22, azi: 92.68, tvd: 6044.64, north: -9.44, east: 37.24 },
    { md: 6217, inc: 17.99, azi: 105.1, tvd: 6208.52, north: -15.55, east: 81.5 },
    { md: 6397, inc: 18.84, azi: 111.6, tvd: 6379.77, north: -32.67, east: 134.17 },
    { md: 6577, inc: 28.16, azi: 115.37, tvd: 6545.43, north: -61.38, east: 198.05 },
    { md: 6747, inc: 37.72, azi: 106.36, tvd: 6688.12, north: -95.54, east: 283.29 },
    { md: 6927, inc: 48.11, azi: 99.88, tvd: 6820.02, north: -119.45, east: 402.97 },
    { md: 7097, inc: 58.37, azi: 102.47, tvd: 6921.06, north: -145.52, east: 536.84 },
    { md: 7277, inc: 65.19, azi: 101.43, tvd: 7005.91, north: -178.74, east: 691.9 },
    { md: 7447, inc: 64.66, azi: 102.72, tvd: 7075.15, north: -210.41, east: 843.86 },
    { md: 7627, inc: 64.36, azi: 106.78, tvd: 7149.5, north: -252.01, east: 1002.37 },
    { md: 7797, inc: 56.88, azi: 110.15, tvd: 7235, north: -299.32, east: 1141.23 },
    { md: 7977, inc: 58.05, azi: 109.53, tvd: 7333.16, north: -350.55, east: 1283.11 },
    { md: 8147, inc: 65.45, azi: 109.94, tvd: 7411.72, north: -400.87, east: 1425.09 },
    { md: 8327, inc: 66.07, azi: 109.2, tvd: 7486.17, north: -456.71, east: 1579.14 },
    { md: 9239, inc: 66.64, azi: 103.7, tvd: 7845.66, north: -637.67, east: 2395.44 },
    { md: 10955, inc: 68.6, azi: 101.44, tvd: 8558.83, north: -1039.83, east: 3901.17 },
  ],
  target: 'tvd',
};
