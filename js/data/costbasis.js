// Sourced cost basis, cost indices, fiscal terms and emission factors of the economics suite.
// Every entry carries the year its number refers to and where it was read; entries that could not be traced to an
// open publication say so (`source: null`, status 'engineering estimate'). Costs are stored in the money of their own
// basis year and are moved to the common basis year (and from there to the evaluation year) with a cost index.
import { FISCAL as FISCAL_TABLE, fiscalOf } from './fiscal.js';

const RETRIEVED = '2026-10-08', RETRIEVED2 = '2026-10-09';

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
  borehole: { citation: 'U.S. Bureau of Safety and Environmental Enforcement, BSEE Data Center, Borehole raw data (file mv_boreholes_all.txt of 8 October 2026: spud date, total-depth date, status date, measured depth, rotary-table elevation and water depth of every Gulf of Mexico borehole); regressions computed for this cost basis', url: 'https://www.data.bsee.gov/Well/Files/BoreholeRawData.zip', licence: 'public domain (U.S. government data)', retrieved: RETRIEVED2 },
  rigs24: { citation: 'Edralin, C. (Westwood RigLogix), Offshore Drilling 2025: 3 Things to Watch During a Year of Market Corrections, Offshore Engineer, 24 December 2024', url: 'https://www.oedigital.com/news/520595-offshore-drilling-2025-3-things-to-watch-during-a-year-of-market-corrections', licence: 'copyright of the publisher; single cited fact', retrieved: RETRIEVED2 },
  nexansUmb: { citation: 'Nexans, press release "Nexans wins a 33.5 million Euro subsea umbilical contract from Petrobras for Brazil\'s Tambau and Urugua deepwater gas and oil fields", Paris, 27 July 2009 (regulated-information filing, Autorité des marchés financiers open data)', url: 'https://echanges.dila.gouv.fr/OPENDATA/AMF/BWR/2009/07/FCBWR051507_20090727.pdf', licence: 'company regulatory release in a French government open-data repository; cited facts', retrieved: RETRIEVED2 },
  prysmianUmb: { citation: 'Offshore magazine, Petrobras awards large-scale umbilicals order to Prysmian, 20 July 2021', url: 'https://www.offshore-mag.com/subsea/article/14207154/petrobras-awards-large-scale-umbilicals-order-to-prysmian', licence: 'copyright of the publisher; single cited fact', retrieved: RETRIEVED2 },
  dehTyrihans: { citation: 'Offshore magazine, Nexans to heat Tyrihans flowlines, 23 April 2006', url: 'https://offshore-mag.com/pipelines/article/16792721/nexans-to-heat-tyrihans-flowlines', licence: 'copyright of the publisher; single cited fact', retrieved: RETRIEVED2 },
  dehFossekall: { citation: 'Nexans, press release "Nexans wins a 20 million euro DEH system contract for Statoil\'s new Fossekall Dompap oil and gas field on the Norwegian Continental Shelf", Paris, 3 May 2011 (regulated-information filing, Autorité des marchés financiers open data)', url: 'https://echanges.dila.gouv.fr/OPENDATA/AMF/BWR/2011/05/FCBWR066614_20110503.pdf', licence: 'company regulatory release in a French government open-data repository; cited facts', retrieved: RETRIEVED2 },
  insJsm: { citation: 'ShawCor Ltd., press release "Shawcor Announces Contract To Provide Gulf Of Mexico Pipe Coating And Insulation Services For The Jack/St. Malo Project", Toronto, 22 December 2010', url: 'https://prnewswire.com/news-releases/shawcor-announces-contract-to-provide-gulf-of-mexico-pipe-coating-and-insulation-services-for-the-jackst-malo-project-112340514.html', licence: 'company press release; cited facts', retrieved: RETRIEVED2 },
  insVega: { citation: 'Offshore magazine, Bredero Shaw wins contract to coat pipe for Vega project, 20 June 2008', url: 'https://www.offshore-mag.com/pipelines/article/16779044/bredero-shaw-wins-contract-to-coat-pipe-for-vega-project', licence: 'copyright of the publisher; single cited fact', retrieved: RETRIEVED2 },
  towler: { citation: 'Tamarona, P.B., Vlugt, T.J.H., Ramdin, M., OpenPyTEA: an open-source Python toolkit for techno-economic assessment of chemical process plants and energy systems, SoftwareX 35, 102816, 2026, doi:10.1016/j.softx.2026.102816 — data files cost_correlations.csv (purchased-cost correlations of Towler & Sinnott, Chemical Engineering Design, table 7.2, 2010 basis, as reproduced there) and cepci_values.csv, package version 3.0.0', url: 'https://pypi.org/project/openpytea/3.0.0/', licence: 'MIT licence (open-source reproduction of the published correlations)', retrieved: RETRIEVED2 },
  treeAker: { citation: 'Offshore magazine, Aker Kvaerner awarded 18 subsea christmas trees by Petrobras, 11 April 2007', url: 'https://offshore-mag.com/subsea/article/16798807/aker-kvaerner-awarded-18-subsea-christmas-trees-by-petrobras', licence: 'copyright of the publisher; single cited fact', retrieved: RETRIEVED },
  boostVigdis: { citation: 'Equinor ASA, news release "Boosting Vigdis", 5 December 2018', url: 'https://www.equinor.com/en/news/2018-12-05-vigdis.html', licence: 'company news release; cited facts', retrieved: RETRIEVED2 },
  boostDraugen: { citation: 'Offshore magazine, Shell contracts Framo for subsea booster pump at Draugen, 27 August 2012', url: 'https://www.offshore-mag.com/subsea/article/16784248/shell-contracts-framo-for-subsea-booster-pump-at-draugen', licence: 'copyright of the publisher; single cited fact', retrieved: RETRIEVED2 },
  fx: { citation: 'European Central Bank, euro foreign exchange reference rates, US dollar per euro, annual averages (series EXR.A.USD.EUR.SP00.A); Norges Bank, exchange rates, Norwegian krone per US dollar, annual averages (series EXR/A.USD.NOK.SP)', url: 'https://data-api.ecb.europa.eu/service/data/EXR/A.USD.EUR.SP00.A?format=csvdata', urlNok: 'https://data.norges-bank.no/api/data/EXR/A.USD.NOK.SP?format=csv', licence: 'ECB and Norges Bank statistics, free re-use with attribution', retrieved: RETRIEVED2 },
  jpt: { citation: 'Journal of Petroleum Technology, 2026 Offshore Challenge: Softening Demand Puts the Brakes on Day Rates, 21 November 2025 (Westwood data)', url: 'https://jpt.spe.org/2026-offshore-challenge-softening-demand-puts-the-brakes-on-day-rates', licence: 'copyright SPE; single cited fact', retrieved: RETRIEVED },
};
const E = (key, label, value, unit, basisYear, low, high, source, note, extra = {}) => Object.freeze({ key, label, value, unit, basisYear, low, high, source: source || null, status: source ? extra.derived ? 'derived from sourced numbers' : 'sourced' : 'engineering estimate', note, index: 'machinery', ...extra });
// an engineering estimate: the basis it rests on and a low / high range that is deliberately wide (default × 0.5 … × 2)
const EST = (key, label, value, unit, input, basis, low = null, high = null) => E(key, label, value, unit, BASIS_YEAR, low ?? +(value * 0.5).toPrecision(2), high ?? +(value * 2).toPrecision(2), null, `Engineering estimate. Basis: ${basis} Searched without result: contract-award releases, regulator benchmarks, open cost reports and theses. Replace by a quotation.`, { input });

// ---- deep-water well-cost model ------------------------------------------------------------------------------------------
/**
 * Coefficients of the well-cost model, each from an open source (see the matching COST_ENTRIES rows).
 * Drilling time: ln(days from spud to total depth) = const + perKmBml × (measured depth below the mudline, km) + perKmWater × (water depth, km)
 * + development (for development wells), least squares on the original holes of 908 Gulf of Mexico wells in at least 600 m of water spudded 2005–2025.
 * Cost: rig charter = day rate × (drilling + completion days); rig and related cost = charter ÷ rigShare; drilling and completion = that ÷ rigRelatedShare.
 */
export const WELL_MODEL = Object.freeze({
  time: Object.freeze({ const: 2.4561, perKmBml: 0.21372, perKmWater: 0.31034, development: 0.22319, sdLog: 0.6855, n: 908, nDevelopment: 199, r2: 0.31, se: Object.freeze({ const: 0.0959, perKmBml: 0.01106, perKmWater: 0.04097, development: 0.05518 }) }),
  completion: Object.freeze({ p10: 26, p50: 54, p90: 93, n: 221 }), // days from total depth to the completed status, development wells in at least 600 m of water
  rigShare: 0.43, rigRelatedShare: 0.89,
  rigRate: Object.freeze({ value: 400, low: 380, high: 500, year: 2025 }), // k$ per day
  equipment: Object.freeze({ value: 13, low: 11, high: 15, year: 2015 }), // production and wellhead equipment including the tree and a down-hole pump, M$
  esp: Object.freeze({ value: 4, low: 3, high: 5, year: 2015 }),
  wellCorrelation: 0.5, // correlation of the cost of wells drilled in one campaign (same rig, same geology): an assumption, no open figure
  published: Object.freeze({ year: 2015, miocene: Object.freeze({ low: 70, high: 165, mean: 120, depthLo: 6100, depthHi: 7300 }), all: Object.freeze({ low: 60, high: 240 }), lucius: 103 }),
});
const ZS = Object.freeze([-1.8339, -1.2816, -0.9674, -0.7279, -0.5244, -0.3407, -0.1679, 0, 0.1679, 0.3407, 0.5244, 0.7279, 0.9674, 1.2816, 1.8339]), Z90C = 1.2815515655446004;
/**
 * Drilling and completion cost of one deep-water well (M$, money of `year`).
 * { waterDepth (m), mdBml (measured depth below the mudline, m), type: 'subsea' | 'dry', development (bool), rigRate (k$/d), year, nWells }.
 * Returns the build-up at the median durations ({ drillDays, complDays, rigDays, rigCharter, rigRelated, gross, equipment, tree, atMedians }), where `gross`
 * includes production and wellhead equipment with tree and down-hole pump; the distribution for one well of what the suite charges as drilling and
 * completion ({ p10, p50, p90, mean, sdLog } from the scatter of drilling time, completion time and rig rate; the subsea tree and wellhead equipment are
 * costed with the subsea scope, a dry-tree well keeps its surface tree and wellhead in the well cost); `cost` = the P50, which is the estimate; and the
 * low / high multipliers of a campaign of nWells wells ({ lo, hi } = P10 and P90 of campaign cost ÷ median).
 */
export function wellCostModel({ waterDepth = 1350, mdBml = 3400, type = 'subsea', development = true, rigRate = WELL_MODEL.rigRate.value, year = BASIS_YEAR, nWells = 1 } = {}) {
  const W = WELL_MODEL, t = W.time, k = W.rigShare * W.rigRelatedShare, fe = bundledFactor('machinery', W.equipment.year, year), rate = Math.max(rigRate, 1);
  const drillDays = Math.exp(t.const + (t.perKmBml * Math.max(mdBml, 0)) / 1000 + (t.perKmWater * Math.max(waterDepth, 0)) / 1000 + (development ? t.development : 0)), complDays = W.completion.p50;
  const equipment = W.equipment.value * fe, esp = W.esp.value * fe, tree = equipment - esp, off = type === 'dry' ? esp : equipment, grossOf = (days, r) => (r * days) / 1000 / k;
  const gross = grossOf(drillDays + complDays, rate), sC = Math.log(W.completion.p90 / W.completion.p10) / (2 * Z90C), rs = [[W.rigRate.low / W.rigRate.value, 0.3], [1, 0.4], [W.rigRate.high / W.rigRate.value, 0.3]], pts = [];
  for (const zi of ZS) for (const zj of ZS) for (const [f, w] of rs) pts.push([Math.max(grossOf(drillDays * Math.exp(t.sdLog * zi) + complDays * Math.exp(sC * zj), rate * f) - off, 0.05 * gross), w]);
  pts.sort((a, b) => a[0] - b[0]);
  const wSum = pts.reduce((a, p) => a + p[1], 0), qOf = (p) => { let acc = 0; for (const [v, w] of pts) { acc += w / wSum; if (acc >= p) return v; } return pts[pts.length - 1][0]; };
  const p10 = qOf(0.1), p50 = qOf(0.5), p90 = qOf(0.9), mean = pts.reduce((a, p) => a + p[0] * p[1], 0) / wSum, sdLog = Math.log(p90 / p10) / (2 * Z90C), n = Math.max(1, Math.round(nWells)), sdN = sdLog * Math.sqrt((1 + (n - 1) * W.wellCorrelation) / n);
  return { drillDays, complDays, rigDays: drillDays + complDays, rigCharter: grossOf(drillDays + complDays, rate) * k, rigRelated: gross * W.rigRelatedShare, gross, equipment, tree, cost: p50, atMedians: Math.max(gross - off, 0), p10, p50, p90, mean, sdLog, campaign: { lo: Math.exp(-Z90C * sdN), hi: Math.exp(Z90C * sdN), sdLog: sdN }, perMetre: (p50 * 1e6) / Math.max(mdBml, 1), publishedLow: WELL_MODEL.published.miocene.low * fe, publishedHigh: WELL_MODEL.published.miocene.high * fe, withEquipment: p50 + off };
}
const WELL_REF = wellCostModel();


/**
 * The cost basis. `value` is the published number in the money of `basisYear`; `index` names the index that moves it to the
 * common basis year ('machinery', 'steel' or 'none' for ratios and current market prices); `input` is the suite input it feeds;
 * `opex: true` marks evaluation-year operating-cost inputs, which are not escalated a second time.
 */
export const COST_ENTRIES = Object.freeze([
  // ---- wells and subsea
  E('wellCost', 'Reference well, drilling and completion without tree and wellhead equipment: 3,400 m below the mudline in 1,350 m of water (well-cost model)', +WELL_REF.cost.toPrecision(3), 'M$ per well', BASIS_YEAR, +WELL_REF.p10.toPrecision(3), +WELL_REF.p90.toPrecision(3), SRC.borehole, `Model: rig charter = day rate × (drilling days + completion days); ÷ ${WELL_MODEL.rigShare} (share of the floating rig in rig and related cost) ÷ ${WELL_MODEL.rigRelatedShare} (share of rig and related cost in drilling and completion); less the production and wellhead equipment, which the suite costs as the subsea tree. Reference well at the median durations: ${WELL_REF.drillDays.toFixed(0)} + ${WELL_REF.complDays} rig days at ${WELL_MODEL.rigRate.value} k$/d = ${WELL_REF.atMedians.toFixed(0)} M$. The estimate is the P50 of the distribution that follows from the scatter of drilling time, completion time and day rate; low and high are its P10 and P90. At the depth of the published Miocene wells the same model gives 128–135 M$ at the median durations against a published average of about 120 M$ (2015).`, { input: 'wellCost', derived: true, index: 'none' }),
  E('wellPublished', 'Published range: deep-water well, drilling and completion with wet tree and down-hole pump (Miocene play, Gulf of Mexico, 20,000–24,000 ft below sea level)', 120, 'M$ per well', 2015, 70, 165, SRC.eia, 'Average about 120 M$; all deep-water plays 60–240 M$; Lucius (six subsea wells of about 19,000 ft) 103 M$ each. Not used as an input: the well-cost model is compared with it. The reference well ends 3,950 m (13,000 ft) below sea level, far shallower than this play.', { input: 'comparison with the well-cost model' }),
  E('wellTime', 'Drilling time of deep-water wells against depth: ln(days) = 2.456 + 0.2137 × depth below mudline (km) + 0.3103 × water depth (km) + 0.2232 (development wells)', 45.8, 'days (reference well, median)', 2025, 19, 110, SRC.borehole, 'Least squares on 908 original holes (199 development wells) in at least 600 m of water, spudded 2005–2025, 5–400 days from spud to total depth; R² 0.31, standard deviation of ln(days) 0.686, all four coefficients more than four standard errors from zero. Low and high are the P10 and P90 of the fitted scatter.', { index: 'none', input: 'well-cost model: drilling days', derived: true }),
  E('wellCompletion', 'Completion time of deep-water development wells: days from total depth to the completed status', 54, 'days', 2025, 26, 93, SRC.borehole, 'Median, P10 and P90 of 221 development wells in at least 600 m of water spudded 2005–2025 whose completed status follows total depth by 3–120 days.', { index: 'none', input: 'well-cost model: completion days', derived: true }),
  E('rigShare', 'Floating rig charter as a share of rig and related cost of a deep-water well', 43, '%', 2015, null, null, SRC.eia, '"Almost 43 %" for the floating rig and over 33 % for support and supply vessels; the rest is helicopters, logging, cementing and testing services.', { index: 'none', input: 'well-cost model' }),
  E('rigRelatedShare', 'Rig and related cost as a share of deep-water drilling and completion cost', 89, '%', 2015, null, null, SRC.eia, 'Figure 9-21 and its text.', { index: 'none', input: 'well-cost model' }),
  E('rigRate', 'Drillship day rate, near-term fixtures', 400, 'k$ per day', 2025, 380, 500, SRC.jpt, 'Published as "below 400,000 $ per day" (2025). Low: 378,708 $ per day, new fixtures of the first quarter of 2015 (earned rates then 436,482). High: clean floating-rig rates "above $500,000 per day" in 2024, with one fixture over 600,000 (Westwood RigLogix).', { index: 'none', input: 'rigRate' }),
  E('wellEquipment', 'Production and wellhead equipment of a deep-water well, including tree and down-hole pump', 13, 'M$ per well', 2015, 11, 15, SRC.eia, 'Midpoint of the published 11–15 M$; deducted from the modelled drilling and completion cost because the subsea tree is costed with the subsea scope. Down-hole pump alone 3–5 M$.', { input: 'well-cost model', derived: true }),
  E('treeRef', 'Subsea tree, wellhead and controls', 9, 'M$ per well', 2015, 6, 12, SRC.eia, 'Published: production and wellhead equipment including the electric submersible pump 11–15 M$ per well, pump alone 3–5 M$. The difference of the midpoints is used.', { input: 'costBasis: tree', derived: true }),
  E('treeContract', 'Subsea trees only: 18 dual-bore trees for 2,000 m of water, 50 M$ contract (Petrobras, 2007)', 2.78, 'M$ per tree', 2007, null, null, SRC.treeAker, '50 M$ ÷ 18. Trees without wellhead, controls and installation: about 4 M$ in 2024 money, inside the 9 M$ used for tree, wellhead and controls.', { input: 'cross-check', derived: true }),
  E('boostRef', 'Subsea multiphase boosting station (pump module, template, protection structure): supply contract', 86, 'M$', 2018, 86, 103, SRC.boostVigdis, 'Vigdis: contract "estimated at NOK 700 million" ÷ 8.134 NOK per US$ (2018 average); the partners invest "some NOK 1.4 billion" in the whole project, with the modifications on two platforms and the power umbilical. High: Draugen, 100 M$ in 2012 for a complete system with topside power and control and the umbilical (103 M$ in 2018 money).', { input: 'costBasis: boost', derived: true }),
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
  // ---- entries sourced in the second search (contract awards with published scope and value, an open reproduction of published equipment-cost correlations)
  E('insPrice', 'Wet insulation applied (syntactic polypropylene on a 10-inch deep-water flowline)', 5900, '$ per m³', BASIS_YEAR, 5500, 7500, SRC.insJsm, 'Jack/St. Malo: contract "in excess of US$40 million" for about 92 km of 10-inch pipe with three-layer polypropylene and syntactic polypropylene insulation = at least 435 $ per metre in 2010, 569 $ in 2024 money, or 663 $ per m² of pipe surface. Vega (Offshore magazine, 20 June 2008): 30 M$ for 51 km of 12- and 14-inch pipe = 780 $ per metre and 680 $ per m² in 2024 money. At 670 $ per m², less the anti-corrosion coating, over the 0.093 m³ per metre of the reference 80 mm layer: 5,900 $ per m³. The thickness of the two contracts is not published, which is the main uncertainty of the conversion.', { input: 'insPrice', derived: true, index: 'none', source2: SRC.insVega }),
  E('dehCable', 'Direct electrical heating system: piggyback cable and associated equipment', 1170, '$ per m', BASIS_YEAR, 965, 1420, SRC.dehTyrihans, 'Tyrihans: 23 M€ for 46 km of cable and associated equipment (2006) = 628 $ per metre at 1.256 $ per €, 965 $ in 2024 money. Fossekall Dompap: 20 M€ for a complete system with riser cable, feeder cables and 25 km of piggyback cable (2011) = 1,114 $ per metre at 1.392 $ per €, 1,420 $ in 2024 money. Geometric mean used.', { input: 'dehCable', derived: true, index: 'none', source2: SRC.dehFossekall }),
  E('umbilicalRef', 'Electro-hydraulic control umbilical, supply', 0.58, 'M$ per km', BASIS_YEAR, 0.37, 0.94, SRC.nexansUmb, 'Tambau and Urugua (1,500 m of water): 33.5 M€ for over 65 km in 16 lengths (2009) = 0.72 M$ per km at 1.395 $ per €, 0.94 in 2024 money. Petrobras frame order (Offshore magazine, 20 July 2021): about 92 M€ (108.6 M$) for 350 km of steel-tube and thermoplastic umbilicals with services = 0.31 M$ per km, 0.37 in 2024 money. Geometric mean used; the published installation-only cost is 476,000 $ per mile (Kaiser 2016).', { input: 'costBasis: umbilical', derived: true, index: 'none', source2: SRC.prysmianUmb }),
  E('layRate', 'Lay rate of a pipelay vessel', 2.5, 'km per day', 2019, 2, 4, SRC.boem, 'Published: "Typical lay rates are 2–4 km/day" for S-lay and "2–3 km/day" for J-lay; the midpoint of the J-lay range, which applies in deep water, is used for the 10-inch reference line.', { input: 'layRate', index: 'none', derived: true }),
  E('slugcatcherRef', 'Slug catcher / inlet separator: horizontal carbon-steel pressure vessel, 60 m³, 30 bar design', 0.27, 'M$', BASIS_YEAR, 0.19, 0.4, SRC.towler, 'Published correlation (2010 US$): cost = 10,200 + 31 × (shell mass, kg)^0.85, valid for 160–50,000 kg. Shell mass 25,800 kg for a 3 m × 8.5 m vessel with two heads and a 33 mm wall (thin-wall formula at 30 bar and 138 MPa allowable stress): 184,500 $, × 800 ÷ 550.8 (plant cost index 2024 ÷ 2010) = 268,000 $. Nozzles and internals are not in the correlation. Low and high: −30 % / +50 % (accuracy of a correlation of this kind).', { input: 'costBasis: slugcatcher', derived: true, index: 'none' }),
  E('compressorRef', 'Centrifugal gas compressor with explosion-proof electric motor, 5,000 kW', 6.2, 'M$', BASIS_YEAR, 4.3, 9.3, SRC.towler, 'Published correlations (2010 US$): compressor 580,000 + 20,000 × (driver power, kW)^0.6, valid for 75–30,000 kW; motor −1,100 + 2,100 × kW^0.6, valid to 2,500 kW and extrapolated. At 5,000 kW: 3.89 + 0.35 = 4.24 M$, × 800 ÷ 550.8 = 6.2 M$. Local exponent of the cost with power: 0.51. Low and high: −30 % / +50 %.', { input: 'costBasis: compressor', derived: true, index: 'none' }),
  // ---- entries for which no open source was found: the basis of the estimate and a wide range are stated
  EST('pipPremium', 'Pipe-in-pipe premium: annulus insulation, centralisers, bulkheads, assembly', 650, '$ per m', 'pipPremium', 'above the insulated single pipe (about 600 $ per metre for insulation and coating from the contracts above); the published pipe-in-pipe project costs (2.2 M$ per mile for flowlines with umbilicals, 2014) are the cross-check.', 300, 1500),
  EST('coatPrice', 'Anti-corrosion coating', 60, '$ per m²', 'coatPrice', 'the published offshore coating cost of 288,000 $ per mile for 14–24-inch lines (2014) is about 120 $ per m² including the concrete weight coat; half of it is taken for the corrosion coat alone.', 25, 130),
  EST('fabPerM', 'Welding, non-destructive testing and field joints', 120, '$ per m', 'fabPerM', 'one girth weld with inspection and field-joint coating for each 12.2 m joint, about 1,500 $ a joint.', 50, 300),
  EST('craFactor', 'CRA-clad line pipe ÷ carbon-steel line pipe', 4.5, '–', 'craFactor', 'a 3 mm alloy 625 or 316L layer on a carbon-steel pipe; alloy prices are several times that of carbon steel.', 2.5, 8),
  EST('riserFactor', 'Steel catenary riser ÷ flowline cost per metre', 2.5, '–', 'riserFactor', 'thicker wall, fatigue-class welds, strakes and flex joint; published contracts bundle the risers with the flowlines, so no stand-alone figure exists.', 1.5, 5),
  EST('mobCost', 'Mobilisation and demobilisation', 6, 'M$', 'mobCost', 'about two weeks of the installation spread (471 k$ per day) for transit and set-up.', 2, 15),
  EST('manifoldRef', 'Production manifold with foundation, four slots', 16, 'M$', 'costBasis: manifold', 'no award with a stated value for a manifold alone was found; the whole subsea scope is calibrated to the published cost of a two-well tie-back, which bounds this item together with the others.', 6, 35),
  EST('jumperRef', 'Rigid jumper with connectors', 1.4, 'M$ each', 'costBasis: jumper', 'two collet connectors and 30–50 m of insulated pipe; bounded by the calibration of the subsea scope.', 0.5, 4),
  EST('pletRef', 'Pipeline end termination', 3, 'M$ each', 'costBasis: plet', 'an 80-tonne structure with a hub and a mudmat (two 86-tonne units were installed at Thunder Hawk); bounded by the calibration of the subsea scope.', 1, 8),
  EST('chemlineRef', 'Chemical line in the umbilical', 0.14, 'M$ per km', 'costBasis: chemline', 'about a quarter of the sourced umbilical supply cost for each added tube.', 0.05, 0.4),
  EST('megregenRef', 'MEG regeneration and reclamation package, 300 m³/d (purchased)', 14, 'M$', 'costBasis: megregen', 'a distillation and vacuum-reclaimer package; no award with a stated capacity and value was found.', 5, 40),
  EST('cheminjRef', 'Chemical-injection skid, 10 m³/d (purchased)', 1.2, 'M$', 'costBasis: cheminj', 'storage tanks and high-pressure metering pumps for subsea injection.', 0.4, 3.5),
  EST('pumpRef', 'Booster or export pump with driver, 1,000 kW (purchased)', 2.4, 'M$', 'costBasis: pump', 'the open correlations for general-purpose centrifugal pumps stop at 300 kW and give far lower values; a multistage high-pressure offshore pump set is outside their range.', 0.6, 5),
  EST('dehpowerRef', 'Heating power unit, riser cable and feeder, 2,000 kW (purchased)', 7, 'M$', 'costBasis: dehpower', 'topside transformer, variable-speed drive and compensation; the riser and feeder cables are partly inside the sourced heating-system cost per metre.', 2.5, 15),
  EST('pigtrapRef', 'Pig launcher and receiver pair, 10-inch (purchased)', 1.1, 'M$', 'costBasis: pigtrap', 'two short pressure vessels with quick-opening closures, valves and instruments.', 0.4, 3),
  EST('arrestorCost', 'Buckle arrestor: ring forging, two girth welds, coating', 45, 'k$ each', 'arrestorCost', 'a forged ring of about 0.5 tonne and two extra girth welds.', 15, 120),
  EST('sleeperCost', 'Sleeper, buckle initiator or span support, installed', 350, 'k$ each', 'sleeperCost', 'a 20–40 tonne fabricated structure and about half a day of the installation spread.', 120, 900),
  EST('labourRate', 'Fully loaded cost of one position', 150, 'k$ per year', 'labourRate', 'salary, rotation, social charges and training of an offshore position; the official wage tables could not be opened.', 90, 280),
  EST('ldhiPrice', 'Low-dosage hydrate inhibitor', 9000, '$ per m³', 'ldhiPrice', 'a speciality surfactant formulation; prices are not published.', 4000, 20000),
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
 * Fallback of the shared fiscal table (js/data/fiscal.js, the single source of fiscal terms). The country table this file used to
 * hold has been merged into the shared one and deleted here; the one row left is what the shared table does not carry in a form the regime model can read
 * (the royalty of the production-sharing variant of a mixed regime). A value present in the shared table always wins; a row is deleted
 * once the shared table has it. `model` is how the single-field cash-flow engine represents a regime (royalty %
 * of revenue, one marginal tax rate on profit, cost-oil cap and contractor profit share — `undefined` where no opened source
 * publishes the number, in which case the value entered by the user is kept); `range` holds published low / high values.
 */
export const FISCAL_SUPPLEMENT = Object.freeze({
  AO: { country: 'Angola', label: 'production-sharing agreements', year: 2019, model: P(0, 50, 50, undefined, `Cost-oil cap 50 % (up to 65 % if development cost is not recovered in four to five years). The profit split slides with the contractor's rate of return and is specific to each agreement: the entered split is kept. ${OLD}`),
    terms: { royalty: 'none under production-sharing agreements', tax: 'petroleum income tax 50 %', costOil: '50 % (may rise to 65 %)', profit: '', ringFence: 'each block assessed separately' }, sources: [ey('Angola')], status: 'sourced (2019 guide); profit split not published' },
});
export const FISCAL_NOTE = 'The shared fiscal table of the application lists the countries for which an open publication states the terms. Where the governing terms are contract-specific or price-dependent the table records only what its sources state, and the values entered by hand are kept for the rest. Without a country choice the terms entered by hand apply; they are generic and belong to no country.';

// ---- mapping of the shared fiscal table to the regime model of the cash-flow engine -----------------------------------------
const num = (x) => (typeof x === 'number' && Number.isFinite(x) ? x : null);
/** low / headline / high of one field of a shared-table row, whichever way the row carries a range: { low, headline | mid | value, high } objects, `<key>Low` / `<key>High`, `<key>Range: [low, high]` or `ranges[key]`. */
function triple(row, key) {
  const v = row?.[key], o = v && typeof v === 'object' && !Array.isArray(v) ? v : null, arr = Array.isArray(v) ? v : row?.[`${key}Range`] ?? row?.ranges?.[key] ?? row?.range?.[key];
  let mid = num(o ? o.headline ?? o.mid ?? o.value ?? o.rate : v), lo = num(o?.low ?? o?.lo ?? o?.min ?? row?.[`${key}Low`] ?? row?.[`${key}Min`]), hi = num(o?.high ?? o?.hi ?? o?.max ?? row?.[`${key}High`] ?? row?.[`${key}Max`]);
  if (Array.isArray(arr) && arr.length >= 2) { lo ??= num(arr[0]); hi ??= num(arr[arr.length - 1]); }
  else if (arr && typeof arr === 'object') { lo ??= num(arr.low ?? arr.lo ?? arr.min); hi ??= num(arr.high ?? arr.hi ?? arr.max); }
  if (lo !== null && hi !== null && hi < lo) [lo, hi] = [hi, lo];
  if (mid === null && lo !== null && hi !== null) mid = 0.5 * (lo + hi); // decision: where a source gives a range, the headline is its midpoint
  return { mid, lo: lo !== null && hi !== null && hi > lo ? lo : null, hi: lo !== null && hi !== null && hi > lo ? hi : null };
}
/** ISO codes offered to the user: every country of the shared table plus those only the supplement has. */
export const fiscalCodes = () => [...new Set([...Object.keys(FISCAL_TABLE), ...Object.keys(FISCAL_SUPPLEMENT)])];
/**
 * Fiscal terms of a country in the form the cash-flow engine uses, from the shared table (js/data/fiscal.js) with the fallback rows
 * filling what it lacks. opt: { waterDepth (m, for terms set by water depth) }. Returns null for an unknown code, else { code, country, label, year, regime: 'tax' | 'psc', royalty, taxRate,
 * costOilCap, profitSplit (% or undefined = keep the entered value), ranges: { royalty, taxRate, costOilCap, profitSplit: [low, high] | undefined },
 * note, basis, sources: [{ citation, url }], from: { field: 'shared table' | 'supplement' }, status }.
 */
export function fiscalTerms(code, opt = {}) {
  const c = String(code || '').toUpperCase(), f = opt.row || fiscalOf(c, opt), sp = FISCAL_SUPPLEMENT[c] || null; // opt.waterDepth selects the terrain class of the shared table; opt.row maps a given row (used by the tests)
  if (!f && !sp) return null;
  const m = sp?.model || {}, from = {}, ranges = {}, pick = (key, shared, fallback) => { if (shared !== null && shared !== undefined) { from[key] = 'shared table'; return shared; } if (fallback !== undefined && fallback !== null) { from[key] = 'supplement'; return fallback; } return undefined; };
  const roy = f ? triple(f, 'royalty') : { mid: null }, tax = f ? triple(f, 'headline') : { mid: null }, alt = f ? num(f.marginalTake) ?? num(f.corporateTax) : null, cap = f ? triple(f, 'costOilCap') : { mid: null }, split = f ? triple(f, 'profitSplit') : { mid: null };
  // a production-sharing regime needs sharing terms: the shared table names the regime, the supplement supplies the terms where it has them
  const psc = f?.regime === 'psc' || (f?.regime === 'mixed' && (cap.mid !== null || split.mid !== null)) || (f?.regime !== 'tax-royalty' && m.regime === 'psc') || (!f && m.regime === 'psc'), regime = psc && (m.regime === 'psc' || cap.mid !== null || split.mid !== null || f?.regime === 'psc') ? 'psc' : 'tax';
  from.regime = f?.regime ? 'shared table' : 'supplement';
  // a mixed regime: the shared royalty belongs to the licence variant, the supplement's to the production-sharing one
  const out = { code: c, country: f?.country || sp.country, label: sp?.label || f?.basis || '', year: f?.year ?? sp?.year, regime,
    royalty: f?.regime === 'mixed' && regime === 'psc' && m.royalty !== undefined ? ((from.royalty = 'supplement'), m.royalty) : pick('royalty', roy.mid, m.royalty), taxRate: pick('taxRate', tax.mid ?? alt, m.taxRate), costOilCap: pick('costOilCap', cap.mid, m.costOilCap), profitSplit: pick('profitSplit', split.mid, m.profitSplit), scale: m.scale };
  for (const [key, t] of [['royalty', roy], ['taxRate', tax], ['costOilCap', cap], ['profitSplit', split]]) { const r = t.lo !== null && t.lo !== undefined && from[key] === 'shared table' ? [t.lo, t.hi] : from[key] === 'supplement' || (from[key] === undefined) ? sp?.range?.[key] : sp?.range?.[key] && t.mid !== null && t.mid >= sp.range[key][0] && t.mid <= sp.range[key][1] ? sp.range[key] : undefined; if (r && r[1] > r[0]) ranges[key] = [r[0], r[1]]; }
  out.ranges = ranges;
  out.note = [f?.note, f?.basis && f.basis !== out.label ? `Headline rate: ${f.basis}.` : '', sp && Object.values(from).includes('supplement') ? m.note : ''].filter(Boolean).join(' ');
  out.terms = sp?.terms || { royalty: f?.royalty === null || f?.royalty === undefined ? 'not stated by the source' : `${f.royalty} %`, tax: f?.basis || '', ringFence: '' };
  out.sources = [f?.source, f?.source2, f?.source3, ...(sp && Object.values(from).includes('supplement') ? sp.sources : [])].filter(Boolean).map((x) => ({ citation: x.citation, url: x.url }));
  out.from = from; out.status = f ? (Object.values(from).includes('supplement') ? 'shared fiscal table, completed from the supplement' : 'shared fiscal table') : `supplement only (not yet in the shared fiscal table): ${sp.status}`;
  return out;
}
/** What the supplement holds that the shared table lacks today: { countries: [{ code, country }], fields: [{ code, country, fields: [] }] }. */
export function fiscalMergeList() {
  const countries = [], fields = [];
  for (const [code, sp] of Object.entries(FISCAL_SUPPLEMENT)) {
    const f = fiscalOf(code);
    if (!f) { countries.push({ code, country: sp.country }); continue; }
    const miss = [];
    if ((f.royalty === null || f.royalty === undefined) && sp.model.royalty !== undefined) miss.push(`royalty ${sp.model.royalty} %`);
    if (sp.model.regime === 'psc' && sp.model.costOilCap !== undefined && f.costOilCap === undefined) miss.push(`cost-oil cap ${sp.model.costOilCap} %`);
    if (sp.model.regime === 'psc' && sp.model.profitSplit !== undefined && f.profitSplit === undefined) miss.push(`contractor profit share ${sp.model.profitSplit} %`);
    if (!f.regime && sp.model.regime) miss.push(`regime ${sp.model.regime === 'psc' ? 'production sharing' : 'tax and royalty'}`);
    for (const [k, r] of Object.entries(sp.range || {})) if (triple(f, k === 'taxRate' ? 'headline' : k).lo === null) miss.push(`${k} range ${r[0]}–${r[1]} %`);
    if (miss.length) fields.push({ code, country: sp.country, fields: miss });
  }
  return { countries, fields };
}
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
/**
 * The cost and price basis the suite carried before its numbers were traced to sources (hand-entered and undated), kept so that the
 * results can walk from it to the present basis driver by driver. `recorded` is what the reference case showed on that basis.
 */
export const EARLIER_BASIS = Object.freeze({
  recorded: Object.freeze({ capex: 361, npv: 1271, breakeven: 21 }), // M$, M$, $/bbl
  wells: 2, wellCost: 70, treeInWell: true, subseaCal: 1, vesselRate: 350, layFactor: 1, depthCoef: 0.2, steelPrice: 1800, learnRate: 90, abandon: 60,
  meohPrice: 550, megPrice: 1100, escFactor: 830 / 800, gridCarbon: 0.45, priceModel: 'ou', priceVol: 25, priceKappa: 0.3,
  capexDist: Object.freeze({ dist: 'triangular', lo: 0.9, mode: 1, hi: 1.5 }), costBasis: Object.freeze({ umbilical: { ref: 1.1 }, slugcatcher: { ref: 1.6, exp: 0.6 }, compressor: { ref: 8.5, exp: 0.75 } }), insPrice: 5000, dehCable: 450, boosting: false,
});
