# HydraSlug

**Integrated multiphase flow-assurance, integrity and techno-economic simulation suite** for oil and gas
production systems — an installable web application that runs entirely on the user's device, online or offline.

**Open it:** https://samuelakosaonyejekwe.github.io/slugshydratescomputer/

Seven engineering suites share one case and pass results to each other, forwards and backwards:

| # | Suite | What it solves |
|---|---|---|
| 1 | Fluid, PVT & Phase Behaviour | Equation-of-state flash, phase envelope, properties, laboratory experiments, hydrate curve, wax and asphaltene screening |
| 2 | Geometry, Wells, Network & Equipment | Route and riser profile, wall and insulation, wells and nodal analysis, chokes, pumps, compressors, separators, network hydraulics |
| 3 | Multiphase Thermal-Hydraulics & Slugging | Steady and transient multiphase flow, heat loss, flow regimes, hydrodynamic, terrain and severe slugging |
| 4 | Hydrate & Multiphase Solids Flow Assurance | Hydrate kinetics, particles, deposition and plugging; wax, scale, asphaltene and sand |
| 5 | Operations, Control & Flow-Assurance Management | Shutdown and cooldown, restart, blowdown, pigging, chemical injection, slug control, operating envelope |
| 6 | Integrity, Loads, Risk & Engineering Assessment | Stress and collapse, slug loads and vibration, fatigue, corrosion, erosion, reliability and risk |
| 7 | Economics, Techno-Economics & Decision Analysis | CAPEX, OPEX, cash flow, risk cost, uncertainty, decision analysis and strategy optimisation |

Around the suites:

* **Case & fluid** — one fluid composition, rates and operating conditions used by every suite; worked examples.
* **Global site data** — for any point on Earth (land, shelf or deep water) the browser pulls relief and water depth,
  seabed temperature and the temperature profile of the water column, sea state, weather, today's oil and gas prices,
  exchange rate and national economic, energy and carbon figures from open data services. A bundled world atlas
  answers when a service cannot be reached or the device is offline; every value says where it came from.
* **Data portal** — pipeline routes, well surveys, network tables, CAD, CFD/FEA meshes, terrain, inspection and
  deposit maps and plain tables are recognised, previewed, measured and routed to the suites that use them.
  This is not a CAD program: it reads geometry made elsewhere.
* **Integrated run** — solves the chain 1 → 7 and repeats it until the quantities that couple the suites
  (deposits narrowing the bore, inhibitor dose, operating plan, risk and cost) stop changing.
* **Decision support** — ranked, quantified recommendations and a sustainability scorecard, benchmarked
  against published practice and the live context of the site, with live literature look-up.

Each suite has the same workflow: overview → inputs → model set-up → mesh / step and a three-level
grid-convergence study (Richardson extrapolation, GCI) → geometry → results (KPIs, plots, tables, reports)
→ calibration against measured data and validation on independent data → verification checks → equations.
The *Equations* tab ticks exactly what the built-in engine solves and lists the rest as reference.

## Install and offline use

* Any computer (Chrome, Edge): press **Install app** in the top bar.
* Android: **Install app**, or browser menu → *Add to Home screen*.
* iPhone / iPad: Safari → Share → *Add to Home Screen*.
* A service worker stores the whole application, so every calculation works in aeroplane mode.
* `standalone.html` is the entire application in one file — copy it anywhere and open it.

## Availability

The static build is host-independent (relative paths only). The addresses it is published at are listed in
`js/data/app.js` and checked live on the app's *Install & offline* page. Installed copies and the single-file
edition keep working with no host at all, so an outage of any web address does not stop existing users.
Whenever a browser first sees a new build on the primary address it asks the Internet Archive to capture the
single-file edition, so an independent copy stays current without anyone's computer being on.
`tools/deploy-mirrors.sh` publishes the same build to Cloudflare Pages, Netlify, GitLab or Codeberg.

## Live data

All requests go straight from the user's browser to public HTTPS services on a fixed allow-list — nothing passes
through a server of this application, so the data are as fresh as the services wherever the app is opened:
Open-Meteo (weather, marine, geocoding), NOAA NCEI (global relief), EMODnet / SeaDataCloud (sea climatology),
World Bank (national indicators), Our World in Data (grid carbon, carbon price), NASA POWER (climate),
US EIA daily spot series (Brent, WTI, Henry Hub), an open exchange-rate service, and OpenAlex /
Crossref (research literature for the decision page). Stored copies are refreshed in the background.

## Security

* Strict Content-Security-Policy: no inline scripts, no third-party code, no `eval`.
* No accounts, cookies, analytics or tracking; cases stay in the browser's local storage.
* Imported files are size-limited and parsed as data only; all text is rendered as plain text.
* Network requests are HTTPS-only to a fixed allow-list, without credentials or referrer.

## Development

No dependencies and no build step are needed to run the app:

```bash
python3 -m http.server 8080      # then open http://localhost:8080
npm test                         # suites, geometry importers, route helpers, external-solver bridge
node tests/live.test.mjs         # live data connectors against the real services
node tools/coverage.mjs --strict # every catalogue item solved in-app or handed off
node tools/build.mjs             # stamp version.json + sw.js (and standalone.html if esbuild is installed)
```

`docs/SUITE_CONTRACT.md` describes how a suite module is written and what each suite publishes to the others.
`tools/catalog.py` regenerates the reference catalogue, `tools/atlas_deep.py` the temperature-at-depth atlas and
`tools/prices.mjs` the stand-by commodity prices.

## Scope, evidence and limits

* **What is solved where.** Every item of the specification catalogue (1,531 across the seven suites) is solved
  by the built-in engines on the user's device. The engines are one-dimensional, lumped or small two- and
  three-dimensional models, so the three-dimensional items (LES, DES, DNS, interface capturing, fluid–structure
  coupling) run at coarse, stated resolution; a fully resolved turbulent channel DNS is beyond a browser. For
  production resolution the *External solvers* page writes a complete, ready-to-run case (OpenFOAM, CalculiX,
  preCICE, CoolProp, teqp) from the study and reads the results back. `node tools/coverage.mjs --strict` lists
  the status of every item.
* **Where the numbers come from.** Each suite's *Equations* tab carries a provenance table: every set of literature
  constants, the source it was compared with, and whether it was verified, corrected or could not be checked
  against an openly readable source. Unverified sets are labelled as such.
* **Reference data.** Each suite's calibration tab holds published measurements, reference-fluid tables and
  benchmarks (NIST, PHREEQC test data, public test databases, regulator statistics and open-access papers), each
  with citation, address, licence and retrieval date. The engine predicts them blind and the errors — including
  the cases it misses — are shown and are part of the test suite.
* **Costs and fiscal terms.** The cost basis is dated and escalated with a producer-price index; each entry is
  marked as sourced, derived or an engineering estimate. Tax and fiscal terms, cost indices, interest rates and
  energy prices are read live where an open service exists, with a dated bundled snapshot offline.
* **File formats.** Besides open formats the portal reads DWG (R13–2018, with labels and embedded solids),
  DGN V7 and the basic elements of V8, ACIS SAT/SAB and Parasolid X_T/X_B (recognising pipe runs, bends and
  diameters), CGNS (HDF5 and ADF), Exodus II, MED, NetCDF-4, HDF5, MATLAB v7.3, BAG, GeoPackage, S-57 charts,
  E57 and LAZ/COPC. JT is read when it carries a Parasolid body; pre-R13 DWG and native CAD part files are
  recognised but not read.
* **Shared data cache.** Rate-limited public services (US cost indices, OECD tax rates and others) are fetched
  once per period for all users by a small function on the mirror (`functions/api/feed.js`), with the direct
  source and a dated bundled snapshot as fallbacks. `docs/DATA_SOURCES.md` lists every bundled data set and
  constant set with its source and licence.
* **Use.** These are engineering models for design studies, screening, teaching and decision support. Results
  for a real asset should be calibrated and validated against that asset's data with the built-in tools, and
  checked against the governing design code, before they are relied on.
