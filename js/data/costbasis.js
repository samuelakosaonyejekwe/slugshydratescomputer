// Sourced cost basis, cost indices, fiscal terms and emission factors of the economics suite.
// Every entry carries the year its number refers to and where it was read; entries that could not be traced to an
// open publication say so (`source: null`, status 'engineering estimate'). Costs are stored in the money of their own
// basis year and are moved to the common basis year (and from there to the evaluation year) with a cost index.
const RETRIEVED = '2026-10-08';

/** Year in which the suite's default cost inputs are expressed. */
export const BASIS_YEAR = 2024;

/** Annual means of the monthly producer-price indices used for escalation (2026 is the mean of January–August, preliminary). */
export const COST_INDEX = Object.freeze({
  id: 'PCU333132333132',
  label: 'US producer-price index, oil and gas field machinery and equipment manufacturing',
  source: { citation: 'U.S. Bureau of Labor Statistics, Producer Price Index by Industry: Oil and Gas Field Machinery and Equipment Manufacturing, series PCU333132333132, monthly, not seasonally adjusted (annual means of the monthly values)', url: 'https://api.bls.gov/publicAPI/v1/timeseries/data/PCU333132333132', licence: 'public domain (U.S. government work)', retrieved: RETRIEVED },
  years: Object.freeze([1987, 1988, 1989, 1990, 1991, 1992, 1993, 1994, 1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026]),
  values: Object.freeze([119.6, 124.7, 126.8, 130.1, 136.5, 137.0, 138.4, 141.1, 146.3, 152.7, 156.8, 160.2, 161.4, 162.0, 168.4, 169.6, 171.6, 177.3, 192.8, 210.3, 225.5, 243.8, 248.4, 247.0, 253.4, 260.8, 263.8, 267.8, 267.9, 266.1, 265.4, 268.0, 271.4, 269.5, 273.8, 292.0, 310.0, 323.2, 328.7, 334.6]),
});
export const STEEL_INDEX = Object.freeze({
  id: 'WPU101706',
  label: 'US producer-price index, steel pipe and tube',
  source: { citation: 'U.S. Bureau of Labor Statistics, Producer Price Index by Commodity: Metals and Metal Products: Steel Pipe and Tube, series WPU101706, monthly, not seasonally adjusted (annual means of the monthly values)', url: 'https://api.bls.gov/publicAPI/v1/timeseries/data/WPU101706', licence: 'public domain (U.S. government work)', retrieved: RETRIEVED },
  years: COST_INDEX.years,
  values: Object.freeze([90.9, 99.2, 102.6, 102.6, 100.8, 94.1, 92.8, 96.9, 104.4, 103.2, 106.9, 109.4, 102.5, 106.6, 104.0, 106.6, 113.4, 166.3, 193.3, 200.8, 202.4, 251.7, 215.5, 241.8, 277.1, 279.4, 252.8, 253.0, 228.7, 218.8, 245.1, 280.6, 281.7, 267.1, 387.1, 491.3, 407.0, 372.2, 376.1, 409.1]),
});
/** US consumer price index (2010 = 100), annual, for the inflation fit and the real escalation of the cost index. */
export const CPI_INDEX = Object.freeze({
  id: 'FP.CPI.TOTL', label: 'Consumer price index, United States (2010 = 100)',
  source: { citation: 'World Bank, World Development Indicators, indicator FP.CPI.TOTL "Consumer price index (2010 = 100)", United States (source: International Monetary Fund, International Financial Statistics); last updated 13 July 2026', url: 'https://api.worldbank.org/v2/country/US/indicator/FP.CPI.TOTL?format=json&per_page=100&date=1987:2025', licence: 'CC BY 4.0', retrieved: RETRIEVED },
  years: Object.freeze([1987, 1988, 1989, 1990, 1991, 1992, 1993, 1994, 1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024]),
  values: Object.freeze([52.11, 54.23, 56.85, 59.92, 62.46, 64.35, 66.25, 67.98, 69.88, 71.93, 73.61, 74.76, 76.39, 78.97, 81.2, 82.49, 84.36, 86.62, 89.56, 92.45, 95.09, 98.74, 98.39, 100, 103.16, 105.29, 106.83, 108.57, 108.7, 110.07, 112.41, 115.16, 117.24, 118.69, 124.27, 134.21, 139.74, 143.86]),
});
const at = (ix, year) => { const i = ix.years.indexOf(Math.round(year)); return i >= 0 ? ix.values[i] : year < ix.years[0] ? ix.values[0] : ix.values[ix.values.length - 1]; };
/** Factor that moves money of year `from` to year `to` on one of the bundled indices ('machinery' | 'steel' | 'none'). */
export const bundledFactor = (index, from, to) => (index === 'none' || from === to ? 1 : at(index === 'steel' ? STEEL_INDEX : COST_INDEX, to) / at(index === 'steel' ? STEEL_INDEX : COST_INDEX, from));

// ---- sources (each address was opened and read on the retrieval date) -----------------------------------------------
const SRC = {
  eia: { citation: 'U.S. Energy Information Administration (prepared by IHS Global Inc.), Trends in U.S. Oil and Natural Gas Upstream Costs, March 2016, section IX "Deepwater Gulf of Mexico"', url: 'https://www.eia.gov/analysis/studies/drilling/pdf/upstream.pdf', licence: 'U.S. government publication (analysis by IHS)', retrieved: RETRIEVED },
  kaiser: { citation: 'Kaiser, M.J., A review of deepwater pipeline construction in the U.S. Gulf of Mexico — contracts, cost, and installation methods, Journal of Marine Science and Application 15(3), 288–306, 2016, doi:10.1007/s11804-016-1373-7', url: 'https://html.rhhz.net/jmsa/html/20160308.html', licence: 'copyright of the publisher; free to read; a few cited data points', retrieved: RETRIEVED },
  boem: { citation: 'Kaiser, M.J., The Offshore Pipeline Construction Industry and Activity Modeling in the US Gulf of Mexico, OCS Study BOEM 2019-070, U.S. Bureau of Ocean Energy Management, October 2019, section 11.3 and appendix K', url: 'https://espis.boem.gov/final%20reports/BOEM_2019-070.pdf', licence: 'U.S. government-funded study, publicly released', retrieved: RETRIEVED },
  bsee: { citation: 'ICF Incorporated for the U.S. Bureau of Safety and Environmental Enforcement, Decommissioning Methodology and Cost Evaluation, TAP 738, 2015, chapter 8', url: 'https://www.bsee.gov/sites/bsee.gov/files/2023-03/738aa1%5B1%5D.pdf', licence: 'U.S. government contractor report published by BSEE', retrieved: RETRIEVED },
  kaiserPA: { citation: 'Kaiser, M.J., Review of analytical models can help operators properly determine well P&A costs, Offshore magazine, 11 December 2023 (BSEE data of May 2022)', url: 'https://offshore-mag.com/decommissioning/article/14301713/center-for-energy-studies-louisiana-state-university-review-of-analytical-models-can-help-operators-properly-determine-well-pa-costs', licence: 'copyright of the publisher; single cited fact', retrieved: RETRIEVED },
  parker: { citation: 'Parker, N., Using Natural Gas Transmission Pipeline Costs to Estimate Hydrogen Pipeline Costs, Institute of Transportation Studies, University of California, Davis, report UCD-ITS-RR-04-35, 2004', url: 'https://itspubs.ucdavis.edu/download_pdf.php?id=197', licence: 'open research report; cited facts', retrieved: RETRIEVED },
  ingaa: { citation: 'ICF International for the INGAA Foundation, North America Midstream Infrastructure through 2035: Capitalizing on Our Energy Abundance, March 2014', url: 'https://ingaa.org/wp-content/uploads/2014/03/21527.pdf', licence: 'copyright INGAA Foundation; cited fact', retrieved: RETRIEVED },
  steel: { citation: 'SteelBenchmarker (World Steel Dynamics), Price History: Tables and Charts, report #491, 30 September 2026', url: 'https://steelbenchmarker.com/history.pdf', licence: 'copyright SteelBenchmarker; public report; cited facts', retrieved: RETRIEVED },
  methanex: { citation: 'Methanex Corporation, Methanex Methanol Price Sheet, 30 September 2026 (prices valid for October 2026)', url: 'https://www.methanex.com/wp-content/uploads/Mx-Price-Sheet-Sep-2026-1.pdf', licence: 'publicly posted prices; cited facts', retrieved: RETRIEVED },
  meglobal: { citation: 'MEGlobal, press releases "MEGlobal announces ACP" for arrival June to October 2026 (Asian Contract Price of monoethylene glycol)', url: 'https://www.meglobal.biz/news-and-media/', licence: 'company press release; cited facts', retrieved: RETRIEVED },
  riviera: { citation: 'Wingrove, M., Market on fire for contractors, vessel owners, Riviera Maritime Media, 23 July 2024', url: 'https://www.rivieramm.com/news-content-hub/market-on-fire-for-contractors-vessel-owners-81587', licence: 'copyright of the publisher; single cited fact', retrieved: RETRIEVED },
  redlinger: { citation: 'Redlinger, M., Drilling Down the Bakken Learning Curve, working paper, United States Association for Energy Economics, 11 May 2015', url: 'https://usaee.org/aws/USAEE/asset_manager/get_file/527969', licence: 'open working paper; cited fact', retrieved: RETRIEVED },
  nstaUoc: { citation: 'North Sea Transition Authority, "North Sea cost efficiency improved in 2025" (UKCS Unit Operating Cost report), 29 September 2026', url: 'https://www.nstauthority.co.uk/news-publications/north-sea-cost-efficiency-improved-in-2025/', licence: 'UK public body (Open Government Licence v3.0 as published by the NSTA)', retrieved: RETRIEVED },
  ogaUoc: { citation: 'Oil and Gas Authority, UKCS Operating Costs in 2018, 2019', url: 'https://www.nstauthority.co.uk/media/4iagopkn/cost-report-pdf-versioncommsv2.pdf', licence: 'UK public body (Open Government Licence v3.0 as published by the OGA)', retrieved: RETRIEVED },
  ogaDecom: { citation: 'Oil and Gas Authority, UKCS Decommissioning Benchmarking Report, November 2021, table 1', url: 'https://www.nstauthority.co.uk/media/8172/decom_benchmarking-report-2021_finalv2.pdf', licence: 'UK public body (Open Government Licence v3.0 as published by the OGA)', retrieved: RETRIEVED },
  nstaWells: { citation: 'North Sea Transition Authority, "Reinstated wells produced 16 million boe in 2025" (UKCS Wells Insights Report 2026), 11 August 2026', url: 'https://www.nstauthority.co.uk/news-publications/reinstated-wells-produced-16-million-boe-in-2025/', licence: 'UK public body (Open Government Licence v3.0 as published by the NSTA)', retrieved: RETRIEVED },
  jpt: { citation: 'Journal of Petroleum Technology, 2026 Offshore Challenge: Softening Demand Puts the Brakes on Day Rates, 21 November 2025 (Westwood data)', url: 'https://jpt.spe.org/2026-offshore-challenge-softening-demand-puts-the-brakes-on-day-rates', licence: 'copyright SPE; single cited fact', retrieved: RETRIEVED },
};
const E = (key, label, value, unit, basisYear, low, high, source, note, extra = {}) => Object.freeze({ key, label, value, unit, basisYear, low, high, source: source || null, status: source ? extra.derived ? 'derived from sourced numbers' : 'sourced' : 'engineering estimate', note, index: 'machinery', ...extra });
const EST = (key, label, value, unit, input, note = 'No open publication with this unit cost was found; order-of-magnitude value that should be replaced by a quotation.') => E(key, label, value, unit, BASIS_YEAR, null, null, null, note, { input });

/**
 * The cost basis. `value` is the published number in the money of `basisYear`; `index` names the index that moves it to the
 * common basis year ('machinery', 'steel' or 'none' for ratios and current market prices); `input` is the suite input it feeds;
 * `opex: true` marks evaluation-year operating-cost inputs, which are not escalated a second time.
 */
export const COST_ENTRIES = Object.freeze([
  // ---- wells and subsea
  E('wellCost', 'Deep-water well, drilling and completion (Miocene play, Gulf of Mexico)', 70, 'M$ per well', 2015, 70, 165, SRC.eia, 'Published range for Miocene wells: 70–165 M$ (all deep-water plays: 60–240 M$), rising with well depth; rig and related cost are about 89 % of it. The reference well is 3.4 km long below the mudline, at the short end of that play, so the low end of the range is the default and the top of the range is carried as the high value.', { input: 'wellCost' }),
  E('treeRef', 'Subsea tree, wellhead and controls', 9, 'M$ per well', 2015, 6, 12, SRC.eia, 'Published: production and wellhead equipment including the electric submersible pump 11–15 M$ per well, pump alone 3–5 M$. The difference of the midpoints is used; no tree-only price was found in an open source.', { input: 'costBasis: tree', derived: true }),
  E('subseaSystem5', 'Subsea system for two satellite wells, 5-mile tie-back', 200, 'M$', 2015, null, null, SRC.eia, 'Lower end of figure 9-41 ("near MM$ 200 to over MM$ 500 for a 5 mile to 65 mile tie-in distance"), normal conditions, no gas lift, water injection or chemical treatment. Used to calibrate the subsea scope (see PIPE_COST_FIT).', { input: 'calibration target' }),
  E('subseaSystem65', 'Subsea system for two satellite wells, 65-mile tie-back', 500, 'M$', 2015, null, null, SRC.eia, 'Upper end of the same curve ("over MM$ 500").', { input: 'calibration target' }),
  E('subseaChemUplift', 'Uplift of the subsea system cost for chemical injection and acid-gas service', 37.5, '%', 2015, 30, 45, SRC.eia, 'Midpoint of the published 30–45 %; shown for reference, not applied.', { index: 'none', input: 'reference' }),
  E('hostMod', 'Host modification for a subsea tie-back (processing, umbilical and control system, riser tube)', 60, 'M$', 2015, null, null, SRC.eia, 'One project (Kodiak to Devils Tower); compare with the topsides scope of the estimate.', { input: 'cross-check' }),
  // ---- flowlines and pipelines
  E('flowlineInfield', 'Deep-water infield flowline systems, procurement and installation (mean of 41 projects, 1979–2015)', 3.61, 'M$ per mile', 2014, 0.42, 6.8, SRC.kaiser, 'Table 2: mean 3.61, standard deviation 3.19 M$ per mile in 2014 dollars; low and high are the mean ∓ one standard deviation. 2.24 M$ per km.', { input: 'cross-check and validation' }),
  E('flowlineJulia', 'Insulated 10.75-inch flowlines with steel catenary risers and terminations, 30 miles, 2,200 m water depth (Julia)', 6.0, 'M$ per mile', 2014, 6.0, 7.76, SRC.kaiser, 'The closest published analogue of the reference line. The paper gives 6.0 in its table and 7.76 in its text, from a contract announced as 133–333 M$.', { input: 'cross-check and validation' }),
  E('flowlinePip', 'Pipe-in-pipe flowlines (80 miles) with 56 miles of umbilicals, EPCI (Big Bend, Dantzler, Gunflint)', 2.21, 'M$ per mile', 2014, null, null, SRC.kaiser, 'Three contracts worth 300 M$ in all.', { input: 'cross-check' }),
  E('flowlineHeated', 'Electrically heated 6 × 10-inch pipe-in-pipe flowlines, material and installation (Serrano, Oregano, 2001)', 2.05, 'M$ per mile', 2001, 1.7, 2.4, SRC.kaiser, 'The only published unit cost of a heated flowline that was found; midpoint of the two projects.', { input: 'cross-check', derived: true }),
  E('offshoreMaterials', 'Offshore gas pipeline material cost (pipe, coatings, anodes, buckle arrestors), FERC filings 1995–2014', 814000, '$ per mile', 2014, 368000, 1600000, SRC.boem, 'Mostly 14–24-inch export lines in shallow water.', { input: 'cross-check' }),
  E('offshoreCoating', 'Offshore pipe coating (corrosion and concrete weight coat), FERC filings', 288000, '$ per mile', 2014, 232000, 382000, SRC.boem, 'Not thermal insulation.', { input: 'cross-check' }),
  E('offshoreInchMile', 'Offshore gas pipelines, engineer–procure–install cost per inch of diameter and mile, FERC filings 1995–2014', 136000, '$ per inch-mile', 2014, 21000, 251000, SRC.boem, 'Mean 136, standard deviation 115 thousand $ per inch-mile.', { input: 'cross-check' }),
  E('onshoreInchMile', 'Onshore transmission pipelines, US and Canada, planning average', 155000, '$ per inch-mile', 2012, null, null, SRC.ingaa, 'For comparison with the offshore figure.', { input: 'cross-check' }),
  E('depthUplift', 'Pipeline cost increase beyond 7,000 ft of water', 50, '%', 2015, null, null, SRC.eia, '"Minor cost increase" from 1,000 to 6,500 ft, "over 50 %" beyond 7,000 ft (2,134 m): 0.5 ÷ 2.134 = 0.23 per 1,000 m is used for the depth factor on the lay spread.', { index: 'none', input: 'depthCoef', scale: 0.0046, derived: true }),
  // ---- steel and vessels
  E('steelPrice', 'Steel plate, United States, FOB mill', 1630, '$ per tonne', 2026, 510, 1630, SRC.steel, 'Benchmark of 30 September 2026. Line pipe is made from plate or coil; the mill conversion to X65 pipe is not in this price. Low: world export hot-rolled band.', { index: 'steel', input: 'steelPrice' }),
  E('vesselRate', 'Pipelay vessel charter (three- to four-year fixtures, Brazil)', 263, 'k$ per day', 2024, 242, 284, SRC.riviera, 'Midpoint of the published 242,000–284,000 $ per day; vessel only.', { index: 'none', input: 'vesselRate', derived: true }),
  E('spreadFactor', 'Total offshore spread ÷ vessel charter (non-drilling vessel)', 1.79, '–', 2015, null, null, SRC.bsee, 'Published: multi-service vessel 260,000 $ per day, total spread with services 466,000 $ per day. The ratio is applied to the pipelay charter as the day-rate factor.', { index: 'none', input: 'layFactor', derived: true }),
  E('spreadRate', 'Intervention spread: construction support vessel × spread factor', 268, 'k$ per day', 2024, null, null, SRC.riviera, 'Published charter "up to around US$150,000" per day for subsea construction support vessels (2024), multiplied by the spread factor 1.79.', { index: 'none', input: 'spreadRate', opex: true, derived: true }),
  E('drillship', 'Drillship day rate, near-term fixtures', 400, 'k$ per day', 2025, null, 400, SRC.jpt, 'Published as "below 400,000 $ per day"; the 2015 study quotes 378,708 (new fixtures) and 436,482 $ per day (earned).', { index: 'none', input: 'reference' }),
  E('learnRate', 'Learning: drilling time for each doubling of rig experience', 95, '% of the previous', 2014, 92.8, 98.1, SRC.redlinger, 'Onshore horizontal wells (Bakken): −5.0 % per doubling of rig experience, −1.9 % per doubling of operator experience; an analogue, since no open subsea learning rate was found.', { index: 'none', input: 'learnRate' }),
  // ---- chemicals
  E('meohPrice', 'Methanol, posted reference price, US Gulf Coast', 1148, '$ per m³', 2026, 554, 1148, SRC.methanex, 'Published 1,450 $ per tonne (non-discounted reference, October 2026) × 0.792 t/m³; Asia-Pacific posted 700 $ per tonne gives the low value. Delivery offshore is not included.', { index: 'none', input: 'meohPrice', opex: true, derived: true }),
  E('megPrice', 'Monoethylene glycol, Asian contract price, CFR main ports', 979, '$ per m³', 2026, 857, 979, SRC.meglobal, 'Published 880 $ per tonne for arrival October 2026 (770–880 since June 2026) × 1.113 t/m³. A seller nomination; delivery offshore is not included.', { index: 'none', input: 'megPrice', opex: true, derived: true }),
  // ---- decommissioning
  E('abandonWell', 'Abandonment of a deep-water subsea well, Gulf of Mexico (838 wells with probabilistic estimates)', 17.2, 'M$ per well', 2022, 8, 17.2, SRC.kaiserPA, 'Regulator estimates, not out-turn cost. The 2015 BSEE study gives 8–16 M$ per wet-tree well.', { input: 'abandon' }),
  E('abandonLine', 'Pipeline decommissioning: flush, cut, plug and abandon in place', 27.5, '$ per foot', 2015, 15, 40, SRC.bsee, 'Midpoint of the published 15–40 $ per foot ("highly variable").', { input: 'abandon', derived: true }),
  E('abandonUmbilical', 'Umbilical removal', 6, '$ per foot', 2015, 2, 10, SRC.bsee, 'Midpoint of the published 2–10 $ per foot.', { input: 'abandon', derived: true }),
  E('abandonAddOn', 'Decommissioning add-ons: engineering 8 %, work provision 10 % (wet tree), weather 20 %', 38, '%', 2015, null, null, SRC.bsee, 'Applied to the pipeline and umbilical items; the well figure is an all-in estimate.', { index: 'none', input: 'abandon' }),
  E('abandonUkWell', 'Decommissioning of a subsea development well, UK central and northern North Sea (P50)', 7.8, '£M per well', 2020, 6.3, 10.0, SRC.ogaDecom, 'P25–P75 range; for comparison (not converted).', { index: 'none', input: 'reference' }),
  // ---- operating cost benchmarks
  E('uocUk2025', 'Unit operating cost, UK continental shelf, 2025', 17.81, '£ per boe', 2025, null, null, SRC.nstaUoc, '19.60 £ per boe in 2024 (2025 prices).', { index: 'none', input: 'benchmark' }),
  E('uocUk2018', 'Unit operating cost, UK continental shelf, 2018', 15.5, '$ per boe', 2018, null, null, SRC.ogaUoc, 'Published as 11.6 £ per boe and 15.5 $ per boe; used as the benchmark for the lifting cost of the case.', { index: 'none', input: 'benchmark' }),
  E('ukDrilling', 'Development drilling, UK continental shelf, 2025', 10000, '£ per metre drilled', 2025, null, null, SRC.nstaWells, 'For comparison.', { index: 'none', input: 'reference' }),
  // ---- entries for which no open source was found
  EST('insPrice', 'Wet insulation applied (syntactic or solid polyurethane / polypropylene)', 5000, '$ per m³', 'insPrice'),
  EST('pipPremium', 'Pipe-in-pipe premium: annulus insulation, centralisers, bulkheads, assembly', 650, '$ per m', 'pipPremium', 'No open unit price; the published pipe-in-pipe project costs above are the cross-check.'),
  EST('dehCable', 'Direct electrical heating cable and anodes', 450, '$ per m', 'dehCable', 'No open unit price; the 2001 heated pipe-in-pipe projects above are the only published reference.'),
  EST('coatPrice', 'Anti-corrosion coating', 60, '$ per m²', 'coatPrice', 'No open price per square metre; the published offshore coating cost per mile above includes concrete weight coat.'),
  EST('fabPerM', 'Welding, non-destructive testing and field joints', 120, '$ per m', 'fabPerM'),
  EST('craFactor', 'CRA-clad line pipe ÷ carbon-steel line pipe', 4.5, '–', 'craFactor'),
  EST('riserFactor', 'Steel catenary riser ÷ flowline cost per metre', 2.5, '–', 'riserFactor', 'No stand-alone riser cost was found: published contracts bundle the risers with the flowlines.'),
  EST('layRate', 'Lay rate of a 10-inch line', 2.5, 'km per day', 'layRate'),
  EST('mobCost', 'Mobilisation and demobilisation', 6, 'M$', 'mobCost'),
  EST('manifoldRef', 'Production manifold with foundation, four slots', 16, 'M$', 'costBasis: manifold'),
  EST('jumperRef', 'Rigid jumper with connectors', 1.4, 'M$ each', 'costBasis: jumper'),
  EST('pletRef', 'Pipeline end termination', 3, 'M$ each', 'costBasis: plet'),
  EST('umbilicalRef', 'Electro-hydraulic control umbilical, supply', 1.1, 'M$ per km', 'costBasis: umbilical', 'No open supply price; the published installation-only cost is 476,000 $ per mile (Kaiser 2016).'),
  EST('chemlineRef', 'Chemical line in the umbilical', 0.14, 'M$ per km', 'costBasis: chemline'),
  EST('slugcatcherRef', 'Slug catcher / inlet separator vessel, 60 m³ (purchased)', 1.6, 'M$', 'costBasis: slugcatcher'),
  EST('megregenRef', 'MEG regeneration and reclamation package, 300 m³/d (purchased)', 14, 'M$', 'costBasis: megregen'),
  EST('cheminjRef', 'Chemical-injection skid, 10 m³/d (purchased)', 1.2, 'M$', 'costBasis: cheminj'),
  EST('pumpRef', 'Booster or export pump with driver, 1,000 kW (purchased)', 2.4, 'M$', 'costBasis: pump'),
  EST('compressorRef', 'Gas compressor train with driver, 5,000 kW (purchased)', 8.5, 'M$', 'costBasis: compressor'),
  EST('dehpowerRef', 'Heating power unit, riser cable and feeder, 2,000 kW (purchased)', 7, 'M$', 'costBasis: dehpower'),
  EST('pigtrapRef', 'Pig launcher and receiver pair, 10-inch (purchased)', 1.1, 'M$', 'costBasis: pigtrap'),
  EST('arrestorCost', 'Buckle arrestor: ring forging, two girth welds, coating', 45, 'k$ each', 'arrestorCost'),
  EST('sleeperCost', 'Sleeper, buckle initiator or span support, installed', 350, 'k$ each', 'sleeperCost'),
  EST('labourRate', 'Fully loaded cost of one position', 150, 'k$ per year', 'labourRate'),
  EST('ldhiPrice', 'Low-dosage hydrate inhibitor', 9000, '$ per m³', 'ldhiPrice'),
]);

/** Abandonment cost of the reference tie-back (M$, evaluation-year money) built from the sourced unit costs above. */
export function abandonmentEstimate({ wells = 2, lineM = 20300, umbilicalM = 21300, year = 2026 } = {}) {
  const e = Object.fromEntries(COST_ENTRIES.map((x) => [x.key, x])), f = (k) => bundledFactor('machinery', e[k].basisYear, year), ft = 1 / 0.3048;
  const well = wells * e.abandonWell.value * f('abandonWell'), line = ((lineM * ft * e.abandonLine.value) / 1e6) * f('abandonLine') * (1 + e.abandonAddOn.value / 100), umb = ((umbilicalM * ft * e.abandonUmbilical.value) / 1e6) * f('abandonUmbilical') * (1 + e.abandonAddOn.value / 100);
  return { wells: well, line, umbilical: umb, total: well + line + umb };
}

// ---- petroleum fiscal terms by country ------------------------------------------------------------------------------
const EY_URL = 'https://ualberta.scholaris.ca/bitstreams/87b03a19-5372-4a2c-8b6d-2f8436500b17/download';
const ey = (chapter) => ({ citation: `EY, Global oil and gas tax guide 2019 (EYGM Limited, May 2019; law as at 1 January 2019), chapter "${chapter}"`, url: EY_URL, year: 2019 });
const R = (r, t, note, more = {}) => ({ regime: 'tax', royalty: r, taxRate: t, note, ...more });
const P = (r, t, cap, split, note, more = {}) => ({ regime: 'psc', royalty: r, taxRate: t, costOilCap: cap, profitSplit: split, note, ...more });
const OLD = 'Terms are those of the 2019 edition of the guide and may have changed since.';
/**
 * Upstream fiscal terms by ISO 3166 alpha-2 code. `terms` quotes what the sources state; `model` is how the single-field
 * cash-flow engine represents them (royalty % of revenue, one marginal tax rate on profit, and for production sharing the
 * cost-oil cap and contractor profit share — `undefined` where no opened source publishes the number, in which case the
 * value entered by the user is kept). `status` says what is sourced and what is a modelling simplification.
 */
export const FISCAL = Object.freeze({
  NO: { country: 'Norway', label: 'continental shelf, company tax and special tax', year: 2026, model: R(0, 78, 'One marginal rate of 78 %. The real system writes investments off at once against the special tax and over six years against company tax; the engine uses the depreciation method entered, so early tax relief is understated.'),
    terms: { royalty: 'none', tax: 'company tax 22 %; special tax 71.8 % on a base net of the calculated company tax; combined marginal rate 78 %', ringFence: 'no field ring fence: shelf income is consolidated; financing costs are not deductible against the special tax' },
    sources: [{ citation: 'Norwegian Petroleum (Ministry of Energy and Norwegian Offshore Directorate), "The petroleum tax system"', url: 'https://www.norskpetroleum.no/en/economy/petroleum-tax/', year: 2026 }], status: 'sourced (official); single-rate representation' },
  GB: { country: 'United Kingdom', label: 'ring-fence corporation tax, supplementary charge and Energy Profits Levy', year: 2026, model: R(0, 78, 'One marginal rate of 78 % for the whole life. The Energy Profits Levy (38 %) is legislated to end on 31 March 2030 or earlier under the price floor, after which the rate is 40 %; first-year allowances and the investment allowance are not modelled.'),
    terms: { royalty: 'none (petroleum revenue tax 0 %)', tax: 'ring-fence corporation tax 30 % + supplementary charge 10 % + Energy Profits Levy 38 % = 78 % marginal', ringFence: 'ring fence: losses from other activities cannot shelter extraction profit; 100 % first-year allowances' },
    sources: [{ citation: 'North Sea Transition Authority, "Taxation" (overview of the UK oil and gas fiscal regime)', url: 'https://www.nstauthority.co.uk/exploration-production/taxation/', year: 2026 }], status: 'sourced (official); single-rate representation' },
  US: { country: 'United States', label: 'federal offshore leases (Gulf), royalty and corporate income tax', year: 2025, model: R(12.5, 21, 'Royalty of leases sold from December 2025; leases of 2008 and later carry 18.75 %, some earlier ones 16.667 %. State taxes do not apply offshore.'),
    terms: { royalty: '12.5 % for new offshore leases (the lowest rate permitted by statute); 18.75 % and 16.667 % for earlier lease sales', tax: 'federal corporate income tax 21 % (2019 guide; not re-confirmed)', ringFence: 'none' },
    sources: [{ citation: 'Bureau of Ocean Energy Management, press release "BOEM Advances First Two OBBBA Offshore Lease Sales", 7 November 2025', url: 'https://www.boem.gov/newsroom/press-releases/boem-advances-first-two-obbba-offshore-lease-sales', year: 2025 }, ey('United States of America')], status: 'royalty sourced (official); tax rate from the 2019 guide' },
  BR: { country: 'Brazil', label: 'concession contracts', year: 2019, model: R(10, 34, 'Special participation (10–40 % of net revenue of large fields) is not modelled. Pre-salt production-sharing blocks pay a 15 % royalty and a bid profit-oil share instead.'),
    terms: { royalty: '10 % under concession (reducible to 5 %); 15 % under production sharing', tax: 'corporate income tax 34 %; special participation 10–40 % for large production volumes (concessions)', ringFence: 'not stated' },
    sources: [{ citation: 'ANP (Agência Nacional do Petróleo, Gás Natural e Biocombustíveis), press kit "Participações governamentais"', url: 'https://www.gov.br/anp/pt-br/canais_atendimento/imprensa/kits-de-imprensa-1/kit-imprensa-part-govern.pdf', year: 2019 }, ey('Brazil')], status: 'royalty sourced (official); tax from the 2019 guide; special participation omitted' },
  NG: { country: 'Nigeria', label: 'deep offshore, Petroleum Industry Act 2021', year: 2021, model: R(7.5, 30, 'Deep-offshore royalty (5 % for fields at or below 50,000 bbl/d) and companies income tax. The price-based royalty (0–10 %) and the 65 % cost-price-ratio limit are not modelled; no opened source gives a hydrocarbon tax for deep offshore, so none is applied.'),
    terms: { royalty: 'onshore 15 %, shallow water 12.5 %, deep offshore 7.5 % (5 % up to 50,000 bbl/d), plus a price-based royalty of 0–10 %', tax: 'companies income tax 30 %; hydrocarbon tax 30 % / 15 % onshore and shallow water', ringFence: 'cost-price-ratio limit of 65 % of gross revenue' },
    sources: [{ citation: 'EY Global Tax Alert, "Nigerian Government signs Petroleum Industry Bill 2020 into law", 2 September 2021', url: 'https://globaltaxnews.ey.com/news/2021-5913-nigerian-government-signs-petroleum-industry-bill-2020-into-law', year: 2021 }], status: 'sourced (secondary); price-based royalty omitted' },
  AO: { country: 'Angola', label: 'production-sharing agreements', year: 2019, model: P(0, 50, 50, undefined, `Cost-oil cap 50 % (up to 65 % if development cost is not recovered in four to five years). The profit split slides with the contractor's rate of return and is specific to each agreement: the entered split is kept. ${OLD}`),
    terms: { royalty: 'none under production-sharing agreements', tax: 'petroleum income tax 50 %', costOil: '50 % (may rise to 65 %)', profit: '', ringFence: 'each block assessed separately' }, sources: [ey('Angola')], status: 'sourced (2019 guide); profit split not published' },
  GH: { country: 'Ghana', label: 'petroleum agreements', year: 2019, model: R(7.75, 35, `Royalty is set in each agreement between 3 % and 12.5 %: the midpoint is used. State participation (carried interest of at least 15 %) is not modelled. ${OLD}`),
    terms: { royalty: '3–12.5 % by agreement', tax: 'income tax 35 %', ringFence: 'by petroleum agreement' }, sources: [ey('Ghana')], status: 'sourced (2019 guide); royalty is the midpoint of the published range' },
  GY: { country: 'Guyana', label: 'Stabroek block, 2016 petroleum agreement', year: 2023, model: P(2, 0, 75, 50, 'No income tax is charged to the contractor in the model: the 2019 guide states that agreements may have the government settle it from its share of profit oil. Blocks under the 2023 model agreement: royalty 10 %, cost recovery 65 %, 50 % split, 10 % corporate tax.'),
    terms: { royalty: '2 % (10 % in the 2023 model agreement)', tax: 'settled by the government from its profit-oil share (10 % corporate tax in the 2023 model)', costOil: '75 % (65 % in the 2023 model)', profit: '50 % to the contractor', ringFence: 'by petroleum agreement' },
    sources: [{ citation: 'Argus Media, "Guyana extends deadline for first oil block bids", 2023', url: 'https://www.argusmedia.com/ja/news-and-insights/latest-market-news/2438577-guyana-extends-deadline-for-first-oil-block-bids', year: 2023 }, ey('Guyana')], status: 'sourced (trade press and 2019 guide); the agreement text was not read' },
  AU: { country: 'Australia', label: 'offshore Commonwealth waters, income tax and petroleum resource rent tax', year: 2019, model: R(0, 58, `Combined rate 40 % + 30 % × (1 − 40 %) = 58 %, taking the resource rent tax as deductible for income tax; the carry-forward uplift of resource-rent-tax deductions is not modelled. ${OLD}`),
    terms: { royalty: 'none offshore (10–12.5 % for onshore and state projects)', tax: 'corporate income tax 30 %; petroleum resource rent tax 40 %', ringFence: 'no project ring fence for income tax; the resource rent tax is assessed by project' }, sources: [ey('Australia')], status: 'rates sourced (2019 guide); combined rate derived' },
  CA: { country: 'Canada', label: 'Newfoundland and Labrador offshore, generic royalty regime', year: 2017, model: R(7.5, 30, 'Top basic royalty of 7.5 % (it steps from 1 % with the R-factor); the net royalty of 10–50 % on R between 1 and 3 is not modelled. Income tax: federal 15 % plus a provincial rate inside the published 11.5–16 % range (15 % assumed).'),
    terms: { royalty: 'basic royalty 1–7.5 % stepping with the R-factor; net royalty 10–50 % sliding for R between 1 and 3', tax: 'federal 15 % plus provincial 11.5–16 %', ringFence: 'by project for the royalty' },
    sources: [{ citation: 'Government of Newfoundland and Labrador, Industry, Energy and Technology, "Generic Offshore Oil Royalty Regime" table, 2017', url: 'https://www.gov.nl.ca/em/files/royalties-2017-gorr-table.pdf', year: 2017 }, ey('Canada')], status: 'royalty sourced (official); net royalty omitted; provincial tax rate assumed inside the published range' },
  MY: { country: 'Malaysia', label: 'production-sharing contracts with PETRONAS', year: 2019, model: P(10, 38, undefined, undefined, `Cost-oil ceiling and profit split follow the revenue-over-cost ratio of each contract and are not published: the entered values are kept. ${OLD}`),
    terms: { royalty: '10 %', tax: 'petroleum income tax 38 %', costOil: '', profit: '', ringFence: 'each contract is a separate chargeable person' }, sources: [ey('Malaysia')], status: 'royalty and tax sourced (2019 guide); sharing terms not published' },
  ID: { country: 'Indonesia', label: 'gross-split production-sharing contracts (2017 rules)', year: 2019, model: P(0, 40, 0, 43, `Gross split: the contractor receives 43 % of oil production before any cost recovery (48 % of gas; the oil split is used for both). Tax 25 % plus 20 % branch profits tax on the remainder = 40 %. Variable and progressive adjustments of the split are not modelled. ${OLD}`),
    terms: { royalty: 'none', tax: 'corporate income tax 25 % and branch profits tax 20 %', costOil: 'no cost recovery under the gross split', profit: 'base split: contractor 43 % of oil, 48 % of gas', ringFence: 'by working area' }, sources: [ey('Indonesia')], status: 'sourced (2019 guide); the gross-split rules were reissued in 2024 and not re-read' },
  SA: { country: 'Saudi Arabia', label: 'crude-oil concession', year: 2025, model: R(15, 50, 'Royalty on the first 70 $/bbl; the marginal royalty of 45 % between 70 and 100 $/bbl and 80 % above is not modelled, so the take is understated at high prices.'),
    terms: { royalty: '15 % up to 70 $/bbl; 45 % of the increment from 70 to 100 $/bbl; 80 % above 100 $/bbl', tax: 'income tax 50 % on the upstream business', ringFence: 'not stated' },
    sources: [{ citation: 'Arab Gulf States Institute, "Aramco and the Saudi Government Budget"', url: 'https://agsi.org/analysis/aramco-and-the-saudi-government-budget/', year: 2025 }, ey('Saudi Arabia')], status: 'sourced (secondary); sliding royalty simplified to its first band' },
  TT: { country: 'Trinidad and Tobago', label: 'exploration and production licences', year: 2019, model: R(12.5, 55, `Petroleum profits tax 50 % plus unemployment levy 5 % (deep-water blocks: 35 %). The supplemental petroleum tax on gross crude income, which depends on the oil price, is not modelled. ${OLD}`),
    terms: { royalty: '12.5 %', tax: 'petroleum profits tax 50 % (35 % deep water) and unemployment levy 5 %; supplemental petroleum tax by price band', ringFence: 'production-sharing operations are ring-fenced' }, sources: [ey('Trinidad and Tobago')], status: 'sourced (2019 guide); supplemental petroleum tax omitted' },
  SN: { country: 'Senegal', label: 'deep-water production-sharing contracts, 2019 Petroleum Code', year: 2019, model: P(8, 30, 65, undefined, `Royalty 6–10 % (midpoint used); cost-oil cap 65 % in deep water (55 % onshore, 60 % shallow, 70 % ultra-deep). The split follows a ratio R set in each contract: the entered split is kept. ${OLD}`),
    terms: { royalty: '6–10 %', tax: 'corporate tax 30 %', costOil: '55 / 60 / 65 / 70 % by water depth', profit: '', ringFence: 'by contract' }, sources: [ey('Senegal')], status: 'sourced (2019 guide); profit split not published' },
  MZ: { country: 'Mozambique', label: 'exploration and production concession contracts', year: 2019, model: P(8, 32, 60, undefined, `Petroleum production tax 6–10 % (midpoint used); the profit split follows an R-factor set in each contract: the entered split is kept. ${OLD}`),
    terms: { royalty: '6–10 % (petroleum production tax)', tax: 'corporate income tax 32 %', costOil: '60 %', profit: '', ringFence: 'not stated' }, sources: [ey('Mozambique')], status: 'sourced (2019 guide); profit split not published' },
  NA: { country: 'Namibia', label: 'petroleum licences', year: 2019, model: R(5, 35, `The additional profits tax (three tranches triggered at after-tax rates of return of 15 % and 20–25 %) is not modelled. ${OLD}`),
    terms: { royalty: '5 %', tax: 'petroleum income tax 35 %; additional profits tax on after-tax net cash flow', ringFence: 'by licence area' }, sources: [ey('Namibia')], status: 'sourced (2019 guide); additional profits tax omitted' },
  SR: { country: 'Suriname', label: 'offshore production-sharing contracts with Staatsolie', year: 2026, model: P(6.25, 36, undefined, undefined, 'The cost-oil ceiling and the R-factor split are contractual and not published: the entered values are kept. The published government take is 60–70 % after costs.'),
    terms: { royalty: '6.25 %', tax: 'income tax 36 %', costOil: '', profit: '', ringFence: 'not stated' },
    sources: [{ citation: 'Staatsolie Maatschappij Suriname N.V., Staatsolie Hydrocarbon Institute, "FAQ"', url: 'https://www.staatsolie.com/en/shi/faq/', year: 2026 }], status: 'royalty and tax sourced (official); sharing terms not published' },
  NL: { country: 'Netherlands', label: 'Mining Act licences', year: 2019, model: R(0, 50, `State profit share of 50 % with corporate income tax credited against it; royalty 0–7 % (0 used). ${OLD}`),
    terms: { royalty: '0–7 %', tax: 'corporate income tax 25 %; state profit share 50 %', ringFence: 'state profit share ring-fenced to mineral production' }, sources: [ey('Netherlands')], status: 'sourced (2019 guide)' },
  DK: { country: 'Denmark', label: 'hydrocarbon tax, all licences', year: 2019, model: R(0, 64, `Combined rate 64 % (corporate tax 25 % deductible from the 52 % hydrocarbon tax). ${OLD}`),
    terms: { royalty: 'none', tax: 'corporate income tax 25 % and hydrocarbon tax 52 %; combined 64 %', ringFence: 'no field ring fence since 2014' }, sources: [ey('Denmark')], status: 'sourced (2019 guide)' },
  KW: { country: 'Kuwait', label: 'foreign companies under the income tax law', year: 2019, model: R(15, 15, OLD), terms: { royalty: '15 %', tax: 'corporate income tax 15 %', ringFence: 'not stated' }, sources: [ey('Kuwait')], status: 'sourced (2019 guide)' },
  GA: { country: 'Gabon', label: 'exploitation and production-sharing contracts, 2014 code', year: 2019, model: P(9, 35, 70, 50, `Royalty 6–12 % and cost-oil cap 65–75 % (midpoints used); the contractor receives at most 50 % of profit oil. ${OLD}`),
    terms: { royalty: '6–12 %', tax: 'corporate income tax 35 %', costOil: '65–75 %', profit: 'at most 50 % to the contractor', ringFence: 'not stated' }, sources: [ey('Gabon')], status: 'sourced (2019 guide); midpoints of published ranges' },
  GQ: { country: 'Equatorial Guinea', label: 'production-sharing contracts, Hydrocarbon Law 8/2006', year: 2019, model: P(13, 35, undefined, undefined, `Royalty of not less than 13 %; cost-oil cap and split are set in each contract: the entered values are kept. ${OLD}`),
    terms: { royalty: 'not less than 13 %', tax: 'corporate income tax 35 %', costOil: '', profit: '', ringFence: 'by contract' }, sources: [ey('Equatorial Guinea')], status: 'royalty and tax sourced (2019 guide); sharing terms not published' },
  CN: { country: 'China', label: 'offshore production-sharing contracts with foreign contractors', year: 2019, model: P(6, 25, 56, undefined, `Resource tax of 6 % of sales in place of royalty; offshore cost-recovery limit 50–62.5 % (midpoint used). The special oil gain levy above 65 $/bbl is not modelled and the allocation factor is contract-specific: the entered split is kept. ${OLD}`),
    terms: { royalty: 'resource tax 6 % of sales (royalty 0–12.5 % for contracts before November 2011)', tax: 'corporate income tax 25 %; special oil gain levy above 65 $/bbl', costOil: '50–62.5 % offshore', profit: '', ringFence: 'no clear ring-fencing rule' }, sources: [ey('China')], status: 'sourced (2019 guide); profit split not published' },
  OM: { country: 'Oman', label: 'production-sharing agreements', year: 2019, model: P(0, 55, undefined, undefined, `Cost-oil cap and split are set in each agreement: the entered values are kept. ${OLD}`), terms: { royalty: 'none', tax: 'income tax 55 % on petroleum income', costOil: '', profit: '', ringFence: 'not stated' }, sources: [ey('Oman')], status: 'tax sourced (2019 guide); sharing terms not published' },
  AR: { country: 'Argentina', label: 'federal and provincial concessions', year: 2019, model: R(12, 30, `Corporate income tax 30 % for 2018–2019 (25 % was legislated from 2020). The 12 % export duty of 2019 is not modelled. ${OLD}`), terms: { royalty: '12 %', tax: 'corporate income tax 30 %; export duty 12 % (2019)', ringFence: 'not stated' }, sources: [ey('Argentina')], status: 'sourced (2019 guide)' },
});
export const FISCAL_NOTE = 'Also researched, but not offered because the governing terms are contract-specific, price-dependent or not published in an opened source: Mexico, Egypt, Qatar, United Arab Emirates, Kazakhstan, Azerbaijan, Russia, India, Colombia, Iraq, Algeria, Libya, Republic of the Congo, Thailand, Vietnam, Brunei. Without a country choice the terms entered by hand apply; they are generic and belong to no country.';
/** Headline marginal rate on UK extraction profit by year: ring-fence tax 30 % + supplementary charge 10 %, plus the Energy Profits Levy of 25 % (2022), 35 % (2023 to 31 October 2024) and 38 % since. */
export const FISCAL_HISTORY = Object.freeze([{ year: 2019, royalty: 0, tax: 40 }, { year: 2022, royalty: 0, tax: 65 }, { year: 2023, royalty: 0, tax: 75 }, { year: 2025, royalty: 0, tax: 78 }, { year: 2026, royalty: 0, tax: 78 }].map(Object.freeze));
export const FISCAL_HISTORY_SOURCE = Object.freeze({ citation: 'North Sea Transition Authority, "Taxation": ring-fence corporation tax 30 %, supplementary charge 10 %, Energy Profits Levy 25 % until 31 December 2022, 35 % until 31 October 2024, 38 % since; 2019 rates from the EY guide', url: 'https://www.nstauthority.co.uk/exploration-production/taxation/', retrieved: RETRIEVED });
// ---- emission factors and unit definitions (official tables) --------------------------------------------------------
const IPCC_E = { citation: '2006 IPCC Guidelines for National Greenhouse Gas Inventories, volume 2 (Energy), chapter 1, table 1.4 "Default CO2 emission factors for combustion"', url: 'https://www.ipcc-nggip.iges.or.jp/public/2006gl/pdf/2_Volume2/V2_1_Ch1_Introduction.pdf', retrieved: RETRIEVED };
const IPCC_C = { citation: '2006 IPCC Guidelines for National Greenhouse Gas Inventories, volume 3 (Industrial Processes and Product Use), chapter 3, tables 3.12, 3.14 and 3.20', url: 'https://www.ipcc-nggip.iges.or.jp/public/2006gl/pdf/3_Volume3/V3_3_Ch3_Chemical_Industry.pdf', retrieved: RETRIEVED };
const DESNZ = { citation: 'UK Department for Energy Security and Net Zero, Greenhouse gas reporting: conversion factors 2026, full set, worksheet "Fuels"', url: 'https://assets.publishing.service.gov.uk/media/6a29392bade52dc0882218a8/ghg-conversion-factors-2026-full-set.xlsx', licence: 'Open Government Licence v3.0', retrieved: RETRIEVED };
const EPA = { citation: 'US Environmental Protection Agency, Emission Factors for Greenhouse Gas Inventories (GHG Emission Factors Hub), 15 January 2025', url: 'https://www.epa.gov/system/files/other-files/2025-01/ghg-emission-factors-hub-2025.xlsx', licence: 'public domain (U.S. government work)', retrieved: RETRIEVED };
const WSTEEL = { citation: 'World Steel Association, Sustainability Indicators Report 2025: CO2 emissions intensity of crude steel production', url: 'https://worldsteel.org/wider-sustainability/sustainability-indicators/', retrieved: RETRIEVED };
const NIST = { citation: 'Thompson, A. and Taylor, B.N., Guide for the Use of the International System of Units (SI), NIST Special Publication 811, 2008, appendix B.8', url: 'https://www.nist.gov/pml/special-publication-811/nist-guide-si-appendix-b-conversion-factors/nist-guide-si-appendix-b8', licence: 'public domain (U.S. government work)', retrieved: RETRIEVED };
/** value, unit, source and how the engine uses it. `estimate` marks the part of a factor that is not from the cited table. */
export const EMISSION_FACTORS = Object.freeze({
  naturalGas: { value: 56100, unit: 'kg CO2 per TJ (net calorific value)', source: IPCC_E, use: 'fuel gas for glycol regeneration, gas-turbine power and flaring' },
  gasOil: { value: 74100, unit: 'kg CO2 per TJ (net calorific value)', source: IPCC_E, use: 'reference' },
  netToGross: { value: 0.18231 / 0.20199, unit: 'net ÷ gross calorific value of natural gas', source: DESNZ, use: 'converts the gross heating value entered for the gas to the net basis of the IPCC factor (0.18231 and 0.20199 kgCO2e per kWh on the gross and net basis)' },
  naturalGasEpa: { value: 53.06, unit: 'kg CO2 per MMBtu (higher heating value)', source: EPA, use: 'cross-check of the gas factor' },
  gridUS: { value: 771.523 * 0.45359237 / 1000, unit: 'kg CO2 per kWh (771.523 lb per MWh, US average, eGRID2023)', source: EPA, use: 'default grid carbon intensity' },
  gridUK: { value: 0.13096, unit: 'kg CO2e per kWh (2026)', source: DESNZ, use: 'reference' },
  marineGasOil: { value: 3245.30441, unit: 'kg CO2e per tonne of marine gas oil', source: DESNZ, use: 'vessel days' },
  vesselFuel: { value: 30, unit: 't of fuel per vessel day', source: null, estimate: true, use: 'vessel days (no open figure for construction-vessel fuel use was found)' },
  turbineEfficiency: { value: 0.337, unit: 'electrical efficiency of an offshore gas turbine', source: null, estimate: true, use: 'gas-turbine power: 56.1 kg/GJ × 3.6 MJ/kWh ÷ 0.337 = 0.60 kg CO2 per kWh' },
  methanol: { value: 0.67, unit: 't CO2 per t of methanol (conventional steam reforming without primary reformer, default)', source: IPCC_C, use: 'embodied emissions of methanol' },
  ethyleneOxide: { value: 0.863, unit: 't CO2 per t of ethylene oxide (air process, default)', source: IPCC_C, use: 'embodied emissions of MEG' },
  ethyleneEthane: { value: 0.95, unit: 't CO2 per t of ethylene (steam cracking of ethane)', source: IPCC_C, use: 'embodied emissions of MEG' },
  ethyleneNaphtha: { value: 1.73, unit: 't CO2 per t of ethylene (steam cracking of naphtha)', source: IPCC_C, use: 'embodied emissions of MEG' },
  ldhi: { value: 3, unit: 't CO2e per t of low-dosage inhibitor', source: null, estimate: true, use: 'embodied emissions of the low-dosage inhibitor (no open factor found)' },
  steel: { value: 1.92, unit: 't CO2 per t of crude steel (2021–2024)', source: WSTEEL, use: 'embodied emissions of the line pipe' },
});
/** Embodied CO2 of monoethylene glycol from the IPCC factors: 0.7097 t ethylene oxide per t MEG, and 0.9097 t ethylene per t oxide at the default 70 % selectivity; mean of ethane and naphtha feed. */
export const megEmbodied = () => { const e = EMISSION_FACTORS, eo = 44.05 / 62.07, eth = (28.05 / 44.05 / 0.7) * eo; return { ethane: e.ethyleneOxide.value * eo + e.ethyleneEthane.value * eth, naphtha: e.ethyleneOxide.value * eo + e.ethyleneNaphtha.value * eth, value: e.ethyleneOxide.value * eo + 0.5 * (e.ethyleneEthane.value + e.ethyleneNaphtha.value) * eth }; };
export const UNIT_DEFS = Object.freeze({
  barrel: { value: 0.1589873, unit: 'm³ per barrel of petroleum (42 US gallons)', source: NIST },
  btu: { value: 1055.056, unit: 'J per British thermal unit (International Table)', source: NIST },
  boe: { value: 5.8, unit: 'million Btu per barrel of oil equivalent', source: { citation: 'United States Code, title 26, section 45K(d)(5), "Barrel-of-oil equivalent"', url: 'https://www.law.cornell.edu/uscode/text/26/45K', licence: 'public domain (U.S. federal statute)', retrieved: RETRIEVED } },
  toe: { value: 41.868, unit: 'GJ per tonne of oil equivalent', source: { citation: 'INSEE, Definitions: Ton of oil equivalent (toe)', url: 'https://www.insee.fr/en/metadonnees/definition/c1355', retrieved: RETRIEVED } },
});


const BY_KEY = () => Object.fromEntries(COST_ENTRIES.map((e) => [e.key, e]));
/** Default of a suite cost input: the sourced value moved to the basis year on its index (rounded to three significant figures), or the fallback. */
export function costDefault(key, fallback) {
  const e = BY_KEY()[key];
  if (!e || !(e.value >= 0)) return fallback;
  const v = e.value * (e.scale ?? 1) * (e.opex ? 1 : bundledFactor(e.index || 'machinery', e.basisYear, BASIS_YEAR));
  return +v.toPrecision(3);
}
export const FAILURE_RECORDS = [];
