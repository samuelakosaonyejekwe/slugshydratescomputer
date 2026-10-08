# Suite module contract

Every suite is one ES module in `js/suites/` that `export default`s a plain object. The generic
workspace (`js/core/suiteview.js`) turns that object into the full workflow (guide, inputs, model
setup, mesh, geometry, results, mesh-sensitivity, calibration, validation, verification, theory). A suite
module therefore contains **physics and declarations only** — no DOM code, no network access, no storage.

Suite modules may import only from `../core/num.js`, `../core/thermo.js`, `../core/pipe.js`, `../core/caseflow.js`, `../core/props.js`
(water/brine properties), `../core/route.js` (pure geometry helpers), `../data/basecase.js` and other suite
modules' **named** engine exports. They must load and run under plain Node (`node tests/run.mjs <id>`) as well
as in the browser and in a module Web Worker (so: no `window`, no `document`, results must be structured-cloneable —
plain objects, arrays, numbers, strings, booleans; no functions, no class instances, no typed arrays in `outputs`).

```js
export default {
  id: 'flow',                // fixed id (table below)
  num: 3,                    // 1..7
  title: 'Multiphase Thermal-Hydraulics & Slugging',
  short: 'Flow · Slugs',     // <= 18 characters, used in navigation
  icon: '🌊',
  tagline: 'One sentence shown under the title.',
  description: 'Two to four sentences: what is solved and how.',
  guide: ['Step 1 …', 'Step 2 …'],          // optional
  implemented: ['darcy–weisbach', …],       // lower-case fragments matched (substring, accents/dashes ignored) against the
                                            // reference catalogue (js/data/catalog.js) to tick what is really solved
  referenceOnly: ['dns'],                    // optional: fragments that must stay unticked even if matched
  equationsNote: 'optional note on model scope and validity limits',

  inputs: [                                 // groups of fields
    { group: 'Rates', tab: 'inputs' | 'setup' | 'mesh', help: '…', showIf: (v) => true,
      fields: [
        { key: 'tIn', label: 'Inlet temperature', unit: '°C', value: 70, min: -20, max: 200, typical: [40, 120], help: '…' },
        { key: 'model', label: 'Holdup model', type: 'select', value: 'beggsBrill', options: [{ value: 'beggsBrill', label: 'Beggs & Brill' }] },
        { key: 'transient', label: 'Run transient', type: 'bool', value: true },
        { key: 'note', label: 'Tag', type: 'text', value: '' },
        { key: 'profile', label: 'Elevation profile', type: 'table', columns: [{ key: 'x', label: 'Distance', unit: 'm' }, { key: 'z', label: 'Elevation', unit: 'm' }], value: [{ x: 0, z: -1250 }] },
        { key: 'terrain', label: 'Seabed terrain', type: 'file', value: null },   // parsed geometry object attached on the Geometry tab
        // table columns may have type: 'text'; any field may have showIf: (v) => boolean
      ] },
  ],
  presets: [{ name: 'Gas-condensate trunk line', values: { tIn: 55 } }],   // >= 3 realistic industrial presets

  // Values offered from the case fluid and from other suites' outputs (all optional chaining!)
  pull: ({ fluid, site, outputs }) => [{ key: 'tIn', value: fluid?.Tin, from: 'Case fluid: inlet temperature' }],
  // Values offered from the Global Site Data page (site.data fields are listed below)
  site: (site) => [{ key: 'tSeabed', value: site.data.seabedTemp, from: 'Seabed temperature at site' }],

  // The engine. May be async. v = current input values keyed by field key.
  // ctx = { fluid, site, outputs, progress(fraction, message), tick() -> Promise (yield to the UI) }
  run(v, ctx) {
    return {
      summary: 'One plain-language sentence describing the outcome.',
      kpis: [{ label: 'Inlet pressure', value: 92.4, unit: 'bara', status: 'ok' | 'warn' | 'bad', help: '…' }],   // >= 6
      warnings: [{ level: 'bad' | 'warn' | 'info', msg: 'Erosional velocity ratio 1.3 exceeds 1.0 at the riser top.' }],
      recommendations: ['Plain-language, specific next action with the number that justifies it …'],
      plots: [ /* plot specs, see below */ ],
      tables: [{ title: 'Segment summary', columns: ['Segment', 'Length (m)'], rows: [[1, 550]], note: '…' }],
      balances: [{ name: 'Total mass', in: 44.95, out: 44.9499 }],   // shown on the Verify tab; must close within 0.1 %
      outputs: { /* machine-readable values consumed by other suites — see table below */ },
    };
  },

  // Numerical-uncertainty study. keys are integer resolution inputs multiplied by the refinement ratio
  // (or, with refine: 'divide', step sizes divided by it). May be an array of several studies.
  mesh: { name: 'Axial grid', keys: ['nCells'], min: 20, note: '…', metrics: [{ label: 'Inlet pressure', unit: 'bara', get: (res) => res.outputs.pIn }] },

  // Parameter estimation + validation. model must be synchronous and fast (< 50 ms per call).
  calibration: {
    note: '…',
    params: [{ key: 'roughUm', label: 'Wall roughness', lo: 5, hi: 500 }],      // keys of numeric input fields
    columns: [{ key: 'qOil', label: 'Oil rate', unit: 'Sm³/d' }, { key: 'dp', label: 'Pressure drop', unit: 'bar' }],
    targets: [{ key: 'dp', label: 'Pressure drop', unit: 'bar' }],              // subset of columns that are measurements
    model: (v) => ({ dp: 12.3 }),                                               // predictions for one operating point
    sample: [ { qOil: 2000, dp: 30 }, … ],        // realistic synthetic "measured" data (>= 8 rows, generated from the
                                                  // model with different parameter values plus a few % noise, hard-coded)
    validationSample: [ … ],                     // different operating points (>= 5 rows)
  },

  // Code/equation verification: conservation, limiting cases, analytical or hand calculations, benchmark problems.
  // >= 12 checks, every one must pass, each with an independent expected value (not the code compared with itself).
  verify() { return [{ name: 'Hagen–Poiseuille limit', expected: 0.064, got: 0.064, tol: 1e-9, pass: true, note: 'f = 64/Re at Re = 1000' }]; },

  // Optional live feed: the shell lets the user follow a local export file and reloads this table input as it grows.
  live: { key: 'log', label: 'Operating log', help: '…' },

  // Optional custom geometry mapping: g is an imported geometry object, d its derived measurements.
  geometry: (g, d) => [{ key: 'length', value: 19500, from: 'Imported route length' }],

  // Optional extra tabs with custom content. `el` is an empty container; build DOM only with api.h(...)
  views: [{ id: 'envelope', label: 'Operating envelope', tip: '…', render(el, api) { /* api = { h, values(), set(k,v), result(), run(), plotCard(spec), dataTable(spec), kpiGrid(items), toast, download, store } */ } }],
};
```

`kpis`, `tables[].rows`, `outputs` and `balances` must never contain `NaN`, `Infinity` or `undefined` (use `null` in
outputs when a value does not exist, and text such as `'—'` in tables). Every numeric field needs a numeric default.
The default inputs must describe the reference case in `js/data/basecase.js` so that all suites agree out of the box,
and a default `run` should finish in about 2 s or less under Node (hard limit 6 s; a mesh study runs it three times).
Use `ctx.progress(f, 'message')` and `await ctx.tick()` inside long loops.

## Plot specs

```js
{ type: 'line', title, xlabel, ylabel, logx, logy, ymin, ymax, xmin, xmax, zeroY, height,
  series: [{ name, x: [], y: [], mode: 'line' | 'points' | 'both' | 'step', dash: true, color }],
  hlines: [{ y, label, color }], vlines: [{ x, label }], note }
{ type: 'bar', title, ylabel, categories: [], series: [{ name, values: [] }], stacked: true }
{ type: 'field', title, xlabel, ylabel, zlabel, zunit, x: [nx], y: [ny], z: [ny][nx] /* rows */, zmin, zmax,
  cmap: 'viridis' | 'turbo' | 'coolwarm' | 'salinity' | 'thermal', contours: 8, equal: true,
  u: [ny][nx], v: [ny][nx], stream: true, vectors: true, mask: [ny][nx] /* true = solid */,
  shapes: [{ x: [], y: [], closed, color, dash, fill }], markers: [{ x, y, label }] }
```
Series x and y must have equal length. A tornado chart is a `bar` plot; a risk matrix or a P–T map is a `field`
plot; a distance–time map (holdup waves, temperature during cooldown) is a `field` plot with x = distance, y = time.

## Units

Engines compute in SI. At the interfaces (inputs, outputs, tables): pressure **bara**, temperature **°C**,
length/elevation **m**, pipe diameter and wall thickness inputs in **mm** (outputs `id`, `wt` in **m**), roughness
input in **µm**, rates at standard conditions (1.01325 bara, 15 °C) in **Sm³/d**, mass rates kg/s, time in s inside
engines and h or d in user-facing fields (say so in the unit), money in US$ (millions as `M$` where labelled).
Elevation `z` is relative to mean sea level, negative below; a profile is `{ x: [], z: [] }` with x the horizontal
distance from the inlet.

## The case fluid (`ctx.fluid`)

```js
{ name, comp: { N2, CO2, H2S, C1, C2, C3, iC4, nC4, iC5, nC5, C6, C7p } /* mol %, water-free well stream */, c7MW, c7SG,
  rateBasis: 'oil' | 'gas' | 'mass', qOil (Sm³/d stock-tank oil), qGas (million Sm³/d), mdot (kg/s hydrocarbon),
  wc (% water cut of the standard liquid, oil basis), qWater (Sm³/d, gas and mass basis), salinity (wt % NaCl equivalent),
  inhibitor: 'none' | 'MeOH' | 'MEG' | 'DEG' | 'TEG' | 'EtOH', inhWt (wt % of the aqueous phase),
  eos: 'PR' | 'SRK' | 'RK' | 'vdW', nPseudo (1–3), Tin (°C at the flowline inlet), Pout (bara arrival pressure),
  Tres (°C), Pres (bara) }
```
`ctx.fluid` may be missing or partial: always go through `fluidModel(ctx)` from `core/thermo.js`, which fills
defaults (`DEFAULT_FLUID`), reuses the property table published by the PVT suite when it belongs to the same fluid
and otherwise builds and caches one.

## Shared kernel (read the JSDoc in the files for exact signatures)

- `core/thermo.js`: `COMPONENTS, COMP_IDS, COMP_LABELS, INHIBITORS, EOS, DEFAULT_FLUID, R, P_STD, T_STD, VM_STD, MW_AIR,
  makeFluid(spec, opts), eosPhase(f, x, Pbar, TK, kind), rachfordRice(z, K), stability(f, Pbar, TK), flashPT(f, Pbar, Tc),
  phaseProps(f, x, Pbar, Tc, kind), enthalpyMolar, interfacialTension, props(f, Pbar, Tc), saturationP(f, Tc), lowerDewP(f, Tc),
  phaseEnvelope(f, opt), stdFlash(f), aqueous(spec), streams(spec, f), waterProps(Pbar, Tc, aq), waterContent(Pbar, Tc),
  hydrateDepression(aq), inhibitorFor(dT, inhId, S), hydrateT0(Pbar, sg), hydrateT(Pbar, sg, aq), hydrateP(Tc, sg, aq),
  pseudoProps(M, SG), tbSoreide(M, SG), buildTable(spec, opt), lookup(table, Pbar, Tc), fluidModel(ctx, override)`.
  `fluidModel(ctx).at(P, T, mScale)` returns every phase property plus actual phase flow rates at the case rate × mScale.
- `core/pipe.js`: `G, frictionFactor(Re, rel, model), hInside, hOutside, uValue({...}), seaTemperature(depth, tSurf, tBed, scale),
  stratifiedLevel, flowPattern, slugVelocity, slugBodyHoldup, slugFrequency, slugLength, slugUnit, severeSlugging,
  gradient(p, model), discretise(profile, n), marchSteady({...})`.
- `core/caseflow.js`: `caseLine(ctx, over)` (line geometry of the case: the network suite's when published, else the reference case),
  `ambientAt(z, line)`, `steadyCase(ctx, opt)` (cached kernel steady solution of the case line), `flowPicture(ctx)` (the flow suite's
  published profile when present, otherwise the kernel estimate, in one shape) and `atX(pic, key, x)`. Suites 4–7 must use
  `flowPicture(ctx)` / `caseLine(ctx)` so that they work both before and after the upstream suites have been run.
- `core/num.js`: `brent, solve1, newton1, newtonN, solveLinear, tridiag, rk4, rk45, nelderMead, diffEvolution,
  levenbergMarquardt, lstsq, linfit, trapz, interp1, linspace, logspace(a, b, n) (a and b are the end VALUES, not exponents),
  metrics, gci, lhs, rng, histogram, mean, std, variance, quantile, sum, clamp, fmt, isNum`.
- `core/props.js` (water and brine; T in °C, S in g/kg): `density, viscosity, cp, conductivityThermal, psat, tsat, latentHeat, …`.
- `data/basecase.js`: `BASE` — the reference deep-water tie-back (profile, diameters, insulation, temperatures, steel, prices).

Kernel files are frozen: do not edit them. If you find a defect or need something added, work around it inside your
own suite file and report it.

## Fixed ids and the outputs each suite publishes

Consumers must read every value with optional chaining and fall back to their own inputs/defaults when it is absent.
Publish every key listed for your suite (use `null` when not applicable); you may add more.

| num | id | outputs |
|---|---|---|
| 1 | `pvt` | `table` (the property table from `buildTable`, for the tuned model), `eos`, `gor` (Sm³/Sm³), `api`, `gasSG`, `mwGas`, `rhoOilStd` (kg/m³), `psat` (bara at reservoir temperature), `psatType` ('bubble' \| 'dew'), `bo` (m³/Sm³ at reservoir conditions), `envelope` `{ T: [], P: [], type: [] }`, `critical` `{ T, P }` \| null, `cricondenbar`, `cricondentherm`, `hydrateCurve` `{ P: [] (bara, ascending), T0: [] (°C, fresh water), T: [] (°C, with the case salt and inhibitor) }`, `hydrateStructure` ('sI' \| 'sII' \| 'sH'), `hydrateDepression` (°C), `waterContent` (kg/Sm³ at arrival conditions), `wat` (wax appearance temperature °C), `waxContent` (wt %), `inhibitorWt` (wt % needed for the coldest condition), `components` `[{ id, z, MW, Tc, Pc, w }]`, `rates` `{ mHC, mW, qOilStd, qGasStd, qWaterStd }` |
| 2 | `net` | `profile` `{ x: [], z: [] }` (≤ 80 points), `length` (m), `id` (m), `wt` (m), `od` (m incl. coatings), `roughness` (m), `uValue` (W/m²K on ID), `layers` `[{ name, t, k }]`, `volume` (m³), `waterDepth` (m), `riserHeight` (m), `riserBaseX` (m), `tSeabed`, `tSeaSurface` (°C), `burial` (m or 0), `kLoss` (total minor-loss coefficient), `bends` `[{ x, angle (deg), radius (m) }]`, `spans` `[{ x, length (m) }]`, `ipr` `{ type, pRes (bara), tRes, pi (Sm³/d/bar), qMax }`, `wellTVD`, `wellMD` (m), `tubingId` (m), `whp` (wellhead pressure bara at the case rate), `operatingRate` (Sm³/d liquid at the nodal solution), `chokeCv`, `chokeOpening` (%), `separatorP` (bara), `slugCatcherVol` (m³), `network` `{ nodes: [{ id, type, p (bara), z }], edges: [{ from, to, type, length, id, q (kg/s), dp (bar) }] }`, `pumpPower`, `compressorPower` (kW), `material` `{ grade, smys, smts, E (MPa), poisson, alphaT, rho }`, `designPressure` (bara), `designTemp` (°C) |
| 3 | `flow` | `profile` `{ x, z, P (bara), T (°C), holdup, vsl, vsg, vm (m/s), rhoM (kg/m³), dpdx (Pa/m), tauW (Pa), tAmb, tHyd, subcooling (°C), regime: [] (text) }` (≤ 120 points), `pIn`, `pOut` (bara), `tIn`, `tOut` (°C), `dpTotal`, `dpFric`, `dpGrav` (bar), `mdot` (kg/s), `qLiq`, `qGas` (actual m³/s at the outlet), `liquidInventory`, `volume` (m³), `residence` (h), `heatLoss` (kW), `uValue`, `slug` `{ type: 'none' \| 'hydrodynamic' \| 'terrain' \| 'severe', freq (1/s), period (s), length (m), lengthMax (m), velocity (m/s), holdupBody, volume (m³ liquid per slug), surge (m³ to accommodate at the receiving vessel), x (m where slugging is strongest) }`, `severeSlugging` (bool), `boe`, `pots`, `erosionalRatio` (max v/v_erosional), `maxVelocity` (m/s), `maxSubcooling` (°C), `hydrateLength` (m of pipe inside the hydrate region), `series` `{ t: [] (s), pIn: [], qLiqOut: [], qGasOut: [] (actual m³/s), holdupOut: [] }` (≤ 400 points, transient), `pInAmplitude` (bar), `turndownRate` (minimum stable rate fraction) |
| 4 | `solids` | `hydrateRisk` (0–1), `maxSubcooling` (°C), `onsetX` (m) \| null, `onsetTime` (h) \| null, `hydrateFraction` (max volume fraction in the liquid), `hydrateRate` (kg/s formed), `depositProfile` `{ x: [], hydrate: [], wax: [], scale: [], total: [] }` (m thickness, same x), `effectiveId` (m, minimum), `roughnessEff` (m), `blockage` (max area fraction lost), `plugTime` (h) \| null, `plugX` (m) \| null, `plugProbability` (0–1), `inhibitorRequired` (wt %), `inhibitorRate` (m³/d), `wat` (°C), `waxRate` (mm/d max), `waxMass` (kg after the simulated time), `piggingInterval` (d), `scaleSI` (max saturation index), `scaleMineral`, `asphalteneRisk` ('low' \| 'medium' \| 'high'), `sandCriticalVelocity` (m/s), `sandBed` (bool), `slurryViscosityFactor` |
| 5 | `ops` | `cooldownTime` (h until the first point reaches the hydrate temperature), `noTouchTime` (h), `maxShutdown` (h), `coldSpotX` (m), `restartPressure` (bara), `restartTime` (h to steady state), `blowdownTime` (h), `blowdownMinT` (°C), `blowdownEndP` (bara), `pigTransit` (h), `pigSurge` (m³), `pigDp` (bar), `inhibitorDose` (wt %), `inhibitorRate` (m³/d), `inhibitorCostPerDay` ($), `heatingPower` (kW), `uptime` (fraction of the year), `deferredVolume` (Sm³ oil per year from planned and unplanned stops), `envelope` `{ qMin, qMax (fraction of the case rate), limits: [text] }`, `controller` `{ kc, ti, td }`, `chokeOpening` (%), `slugSuppressed` (bool), `alarms` `[{ tag, level, msg }]`, `eventsPerYear` `{ shutdowns, pigRuns, blowdowns }` |
| 6 | `integ` | `slugForce` (kN, peak on a 90° bend), `hoopUtil`, `vmUtil`, `collapseUtil` (utilisation 0–1+), `mawp` (bara), `minWallRequired` (mm), `corrosionRate`, `erosionRate` (mm/y), `wallLossRate` (mm/y combined), `remainingLife` (y), `fatigueDamagePerYear`, `fatigueLife` (y), `fivRatio` (excitation frequency / first natural frequency), `naturalFrequency` (Hz), `pof` (annual probability of failure), `reliabilityIndex`, `riskLevel` ('low' \| 'medium' \| 'high' \| 'very high'), `riskCostPerYear` ($), `consequence` ($ per failure), `inspectionInterval` (y), `criticalLocations` `[{ x, mechanism, utilisation }]`, `violations` (count), `upheavalUtil`, `erosionalRatio` |
| 7 | `econ` | `capex`, `opexPerYear`, `npv`, `riskedNpv` (US$), `irr` (%/y), `mirr`, `payback`, `discountedPayback` (y), `pi` (profitability index), `utc` ($/boe), `breakevenPrice` ($/bbl), `eal` (expected annual loss $), `deferredCost` ($/y), `chemicalCost`, `energyCost` ($/y), `carbon` (tCO₂e/y), `carbonIntensity` (kgCO₂e/boe), `npvP10`, `npvP50`, `npvP90`, `probLoss` (probability NPV < 0), `bestOption` (name of the top-ranked flow-assurance strategy), `ranking` `[{ option, npv, capex, risk, score }]` |

## Coupling (forward and backward)

The required input and output data of each suite are the `inputs` and `outputs` lists of its entry in `js/data/catalog.js`
(generated from the corrected module specification); a suite must accept every listed input and compute every listed output,
or leave the item unticked.

Forward: 1 → 2 → 3 → 4 → 5 → 6 → 7. Backward: suite 4 publishes `effectiveId`, `roughnessEff` and `depositProfile`; suites 2
and 3 offer them through `pull` (as optional inputs `depositProfile`-aware effective diameter / roughness) so that a second pass
recomputes pressure drop and slugging with the restricted bore; suite 5 publishes `inhibitorDose` and `chokeOpening`, which
suite 4 and suite 3 may pull; suites 6 and 7 always read the latest of everything. The Integrated-run page repeats the chain
until the coupled quantities stop changing, so engines must be deterministic for identical inputs (seed any random numbers).

## Site data (`site.data`, any field may be missing)

`depth` (m water depth at the site, positive down; 0 on land), `elevation` (m land elevation), `seabedTemp` (°C at the seabed),
`sst` (°C sea-surface temperature), `salinity` (g/kg), `currentSpeed` (m/s), `currentDir` (°), `waveHeight` (m significant),
`wavePeriod` (s), `tideRange` (m), `airTemp` (°C), `airTempMin`, `airTempMax` (°C climatological), `windSpeed` (m/s),
`groundTemp` (°C, soil at about 1 m), `pressure` (hPa), `bathy` `{ lat[], lon[], elev[][] }` (local grid, m, negative below sea),
`seabedSlope` (deg), `inflation` (%/y), `lendingRate` (%/y), `fxPerUSD`, `currency`, `electricityPrice` ($/kWh),
`gridCarbon` (kgCO₂/kWh), `oilPrice` ($/bbl Brent), `gasPrice` ($/MMBtu), `carbonPrice` ($/tCO₂), `taxRate` (%), `country`.

## Sourced parameters and reference validation data

Numbers that come from the literature (model constants, correlation coefficients, property data, default costs) must be
traceable. A suite may export, next to its default object:

```js
// Every constant set the engine relies on, with where it was checked.
export const PROVENANCE = [
  { item: 'Calcite solubility product, log K(T)', used: 'scaleIndices()', source: 'USGS PHREEQC database phreeqc.dat', url: 'https://…',
    retrieved: '2026-10-08', status: 'verified' | 'corrected' | 'unverified', note: 'what was compared and the largest difference found' },
];
```

and the default object may carry

```js
validationData: [
  { id: 'nist-methane-z', title: 'Methane compressibility factor, 250–400 K, 1–600 bar', quantity: 'Z', unit: '–',
    kind: 'reference-fluid' | 'experiment' | 'field' | 'benchmark' | 'market',
    source: { citation: 'Full citation of the publication or database', url: 'https://… (the address the numbers were actually read from)',
              licence: 'public domain / CC BY 4.0 / …', retrieved: '2026-10-08' },
    columns: [{ key: 'T', label: 'Temperature', unit: 'K' }, { key: 'P', label: 'Pressure', unit: 'bar' }, { key: 'Z', label: 'Measured Z' }],
    rows: [{ T: 250, P: 10, Z: 0.9712 }, …],          // >= 8 rows, copied from the source, never invented or "typical"
    target: 'Z',                                       // the measured column
    model: (row) => 0.9705,                            // the engine's blind prediction for that row (no fitting to this data set)
    tolerance: { mape: 2 },                            // acceptance limit the engine is expected to meet (optional: bias, rmse, maxAbs)
    note: 'range of applicability, what the comparison shows' },
],
```

Data sets live in `js/data/ref/<suite id>.js` (one module per suite, imported by that suite) so that the engine file stays
readable. Rules: only numbers that were actually read from the cited address during development; keep the citation and the
address; respect the licence (facts and public-domain or openly licensed data only — no copied copyrighted tables beyond a
few cited data points); never tune the model to a validation set; if the engine misses the tolerance, say so in `note` and
leave the tolerance honest rather than widening it silently. The workspace shows every data set on the calibration tab with
a parity plot and error metrics, and `tests/run.mjs` evaluates them all.

## Hand-off to external open-source solvers

Formulations that cannot run inside a browser at useful resolution (three-dimensional LES/DNS, interface-capturing CFD,
shell and solid finite elements, two-way fluid–structure interaction, multi-parameter reference equations of state) are
covered by `js/core/bridge.js`: it writes ready-to-run case files for open solvers from the case and the suite results and
reads their results back. `bridge.js` exports `HANDOFF = { <suite id>: [{ match: 'lower-case fragment of the catalogue item',
solver: 'OpenFOAM interFoam', generator: 'openfoamVof' }] }`; the Equations tab shows those items as “via external solver”
— a third state between “solved in-app” and “reference”. Suites must not tick such items as implemented.

### Further site fields (sourced, each with a `…Date`/`…Year` and `…Source` companion)

`corporateTaxRate` (%), `taxRate` (% headline tax on upstream profit, with `taxRateBasis`), `fiscalRegime`, `royaltyRate`,
`petroleumTaxRate`, `marginalTake` (%), `costIndex` (with `costIndexBase`, `costIndexSeries { t, v }`), `costEscalation` (%/y),
`steelIndex`, `steelPrice` (US$/t iron ore), `usCpi`, `usInflation5y` (%/y), `gasPriceEurope`, `gasPriceAsia`,
`gasPriceRegional` (US$/MMBtu), `bondYield`, `treasuryBillYield`, `policyRate` (%/y), `electricityPriceSource`.
`site.countryCode` (ISO 3166-1 alpha-2) and `site.country` sit beside `site.data`.
