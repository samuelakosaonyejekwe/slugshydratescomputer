// Published reference data for the economics suite. Every number below was read from the cited address on the
// retrieval date; nothing is interpolated, rounded differently or "typical". The suite attaches its blind
// predictions to these rows (see `validationData` in js/suites/s07_econ.js).
const RETRIEVED = '2026-10-08';

/** Annual average spot prices: Brent ($/bbl) and Henry Hub ($/MMBtu). The Henry Hub page showed no 2013 value on the retrieval date. */
export const PRICE_HISTORY = Object.freeze({
  id: 'eia-annual-prices',
  title: 'Annual average spot prices: Europe Brent FOB and Henry Hub natural gas, 1987–2025',
  source: {
    citation: 'U.S. Energy Information Administration, Petroleum & Other Liquids / Natural Gas data: "Europe Brent Spot Price FOB (Dollars per Barrel)", annual, series RBRTE, and "Henry Hub Natural Gas Spot Price (Dollars per Million Btu)", annual, series RNGWHHD; release of 7 October 2026',
    url: 'https://www.eia.gov/dnav/pet/hist/LeafHandler.ashx?n=PET&s=RBRTE&f=A',
    urlGas: 'https://www.eia.gov/dnav/ng/hist/rngwhhdA.htm',
    licence: 'public domain (U.S. government work)', retrieved: RETRIEVED,
  },
  rows: Object.freeze([
    [1987, 18.53, null], [1988, 14.91, null], [1989, 18.23, null], [1990, 23.76, null], [1991, 20.04, null], [1992, 19.32, null], [1993, 17.01, null], [1994, 15.86, null], [1995, 17.02, null], [1996, 20.64, null],
    [1997, 19.11, 2.49], [1998, 12.76, 2.09], [1999, 17.9, 2.27], [2000, 28.66, 4.31], [2001, 24.46, 3.96], [2002, 24.99, 3.38], [2003, 28.85, 5.47], [2004, 38.26, 5.89], [2005, 54.57, 8.69], [2006, 65.16, 6.73],
    [2007, 72.44, 6.97], [2008, 96.94, 8.86], [2009, 61.74, 3.94], [2010, 79.61, 4.37], [2011, 111.26, 4.0], [2012, 111.63, 2.75], [2013, 108.56, null], [2014, 98.97, 4.37], [2015, 52.32, 2.62], [2016, 43.64, 2.52],
    [2017, 54.13, 2.99], [2018, 71.34, 3.15], [2019, 64.3, 2.56], [2020, 41.96, 2.03], [2021, 70.86, 3.89], [2022, 100.93, 6.45], [2023, 82.49, 2.53], [2024, 80.52, 2.19], [2025, 69.14, 3.52],
  ].map((r) => Object.freeze({ year: r[0], oil: r[1], gas: r[2] }))),
});

/** Annual net oil production of two Norwegian fields with a long decline (million Sm³ per year). The partial year 2026 is left out. */
export const FIELD_PRODUCTION = Object.freeze({
  id: 'sodir-field-production',
  title: 'Annual net oil production of the Draugen and Norne fields, Norwegian continental shelf',
  source: {
    citation: 'Norwegian Offshore Directorate (Sokkeldirektoratet), FactPages, table "Field – Production – Yearly – by field" (column prfPrdOilNetMillSm3), fields DRAUGEN and NORNE',
    url: 'https://factpages.sodir.no/public?/Factpages/external/tableview/field_production_yearly&rs:Command=Render&rc:Toolbar=false&rc:Parameters=f&IpAddress=not_used&CultureCode=en&rs:Format=CSV&Top100=false',
    licence: 'Norwegian Licence for Open Government Data (NLOD)', retrieved: RETRIEVED,
  },
  fields: Object.freeze({
    Draugen: Object.freeze({ first: 1993, oil: Object.freeze([0.123, 3.887, 6.971, 8.444, 10.452, 11.194, 12.139, 11.744, 11.859, 11.067, 7.434, 7.807, 6.018, 4.549, 4.124, 4.04, 3.398, 2.409, 2.173, 2.024, 1.469, 1.805, 1.721, 1.277, 1.339, 1.16, 1.13, 0.986, 0.922, 0.881, 0.782, 0.825, 0.694]) }),
    Norne: Object.freeze({ first: 1997, oil: Object.freeze([0.415, 6.316, 8.305, 10.42, 11.297, 10.27, 8.552, 7.008, 5.267, 4.849, 4.274, 3.482, 2.214, 1.931, 1.615, 0.743, 0.802, 0.706, 0.631, 0.43, 0.57, 0.615, 0.573, 0.546, 0.586, 0.422, 0.261, 0.179, 0.124]) }),
  }),
});

/** Onshore US transmission pipeline construction cost by nominal diameter, mean of the filed projects, year-2000 dollars per mile. */
export const PARKER_TABLE = Object.freeze({
  id: 'parker-pipeline-cost',
  title: 'US onshore transmission pipeline construction cost by diameter (FERC filings compiled by Oil & Gas Journal, 1991–2003)',
  source: {
    citation: 'Parker, N., Using Natural Gas Transmission Pipeline Costs to Estimate Hydrogen Pipeline Costs, Institute of Transportation Studies, University of California, Davis, research report UCD-ITS-RR-04-35, 2004, table 1 (average cost per mile by diameter, year-2000 dollars)',
    url: 'https://itspubs.ucdavis.edu/download_pdf.php?id=197', licence: 'open research report; factual data points with citation', retrieved: RETRIEVED,
  },
  // [diameter in, materials, labour, miscellaneous, right of way, total] — mean $ per mile
  rows: Object.freeze([[4, 60017, 268585, 101668, 56222, 486492], [6, 57863, 239916, 115264, 54364, 467407], [8, 93436, 208658, 139034, 36947, 478076], [10, 102258, 246771, 110033, 43427, 503489], [12, 113981, 404051, 174573, 63389, 755993], [16, 150324, 407615, 214930, 82542, 855411],
    [20, 210178, 491082, 273170, 81100, 1055529], [24, 245372, 574579, 297635, 99112, 1210092], [30, 395461, 637608, 349755, 86631, 1469456], [36, 519622, 764100, 398088, 86900, 1768710], [42, 713651, 998242, 492774, 96377, 2301044]].map((r) => Object.freeze({ d: r[0], materials: r[1], labour: r[2], misc: r[3], row: r[4], total: r[5] }))),
  // the report's own regression of the same data: cost = [a·d² + b·d + c]·miles + fixed
  equations: Object.freeze({ materials: [330.5, 687, 26960, 35000], labour: [343, 2074, 170013, 185000], total: [674, 11754, 234085, 405000] }),
});

/**
 * Deep-water Gulf of Mexico pipeline and flowline contracts with a stated diameter and length (2014 dollars per mile).
 * `d` (inches; the mean when two sizes are given) and `miles` (length of one line) are read from the published description;
 * `ins` and `riser` mark the contracts whose description names insulation or risers.
 */
export const KAISER_PROJECTS = Object.freeze({
  id: 'kaiser-deepwater-pipelines',
  title: 'Deep-water Gulf of Mexico pipeline and flowline contract cost per mile, 1998–2015',
  source: {
    citation: 'Kaiser, M.J., A review of deepwater pipeline construction in the U.S. Gulf of Mexico — contracts, cost, and installation methods, Journal of Marine Science and Application 15(3), 288–306, 2016, doi:10.1007/s11804-016-1373-7, table 3 (press-release cost data, 2014 dollars)',
    url: 'https://html.rhhz.net/jmsa/html/20160308.html', licence: 'copyright of the publisher, free to read; eleven cited data points', retrieved: RETRIEVED,
  },
  rows: Object.freeze([
    ['Ursa', 1998, '18 inch, 47 mi oil; 20 inch, 47 mi gas', 19, 47, 0, 0, 1.45], ['Brutus', 2002, '20 inch, 26 mi oil; 20 inch 24 mi gas', 20, 25, 0, 0, 2.46], ['Falcon', 2002, '32 mi, 10 inch flowline and umbilical', 10, 32, 0, 0, 0.8], ['Glider', 2003, '6 inch, 6 mi flowlines', 6, 6, 0, 0, 10.0],
    ['Droshky', 2008, '8 inch, 36 mi flowline', 8, 36, 0, 0, 1.26], ['Big Foot', 2009, '20 inch, 40 mi oil', 20, 40, 0, 0, 7.05], ['Keathley Canyon', 2013, '20 inch, 215 mi', 20, 215, 0, 0, 2.81], ['Julia', 2014, '30 mi insulated flowlines, risers, PLETs (10.75 inch)', 10.75, 30, 1, 1, 6.0],
    ['Shell', 2014, '27 miles of 8 inch flowlines, SCR and PLETs', 8, 27, 0, 1, 2.77], ['Walker Ridge', 2015, '8 and 10 inch, 170 miles gas', 9, 170, 0, 0, 2.94], ['Stampede', 2015, '18 inch, 16 mi oil', 18, 16, 0, 0, 8.1],
  ].map((r) => Object.freeze({ project: r[0], year: r[1], description: r[2], d: r[3], miles: r[4], ins: r[5], riser: r[6], cost: r[7] }))),
  // quoted in the text of the same paper for the Independence Hub area (ultra-deep-water gas flowlines, dollars of 2007): total M$ and M$ per mile
  text2007: Object.freeze([['two 9 inch flowlines running 56 miles', 9, 56, 0, 136, 2.3], ['8 inch, 25 mile flowline with riser', 8, 25, 1, 103, 4.0], ['10 inch, 22 mile flowline with riser', 10, 22, 1, 64, 2.8]].map((r) => Object.freeze({ description: r[0], d: r[1], miles: r[2], riser: r[3], total: r[4], cost: r[5] }))),
});

const NCS_ROWS = [["2020-2021", "Sverdrup Byggetrinn I", 2015, 132.519, 100.212], ["2020-2021", "Utgard", 2017, 3.466, 2.782], ["2020-2021", "Valhall Flanke Vest", 2018, 5.841, 5.98], ["2020-2021", "Skogul", 2018, 1.591, 2.321], ["2020-2021", "Ærfugl", 2018, 8.849, 8.594], ["2020-2021", "Gullfaks Shetland/Lista fase 2", 2019, 2.275, 2.275], ["2021-2022", "Martin Linge", 2012, 31.434, 63.091], ["2021-2022", "Dvalin", 2017, 11.327, 10.193], ["2021-2022", "Snorre Expansion Project (SEP)", 2018, 20.894, 18.269], ["2021-2022", "Troll Fase 3 Trinn 1", 2018, 8.076, 8.557], ["2021-2022", "Tor II", 2019, 6.3, 6.176], ["2021-2022", "Duva", 2019, 5.675, 5.197], ["2022-2023", "Yme New Development", 2018, 9.357, 12.7], ["2022-2023", "Nova", 2018, 10.665, 11.424], ["2022-2023", "Solveig", 2019, 7.032, 6.235], ["2022-2023", "Hod Nyutvikling", 2020, 6.138, 6.935], ["2023-2024", "Njord Future", 2017, 18.4, 33.2], ["2023-2024", "Bauge", 2017, 4.8, 4.6], ["2023-2024", "Fenja", 2018, 12.5, 12.6], ["2023-2024", "Sverdrup Byggetrinn II", 2019, 50.3, 52.2], ["2023-2024", "Frosk", 2022, 2.3, 2.0], ["2023-2024", "Hywind Tampen", 2020, 5.5, 8.0], ["2024-2025", "Breidablikk", 2021, 22.4, 22.9], ["2024-2025", "Sleipner Kraft frå land", 2021, 1.0, 1.2], ["2024-2025", "Tommeliten A", 2022, 14.5, 13.0], ["2024-2025", "Kobra East and Gekko", 2022, 9.2, 8.5], ["2024-2025", "Eldfisk Nord", 2022, 11.8, 13.7], ["2024-2025", "Kristin Sør", 2022, 7.8, 8.5], ["2025-2026", "Johan Castberg", 2018, 61.4, 87.8], ["2025-2026", "Balder Future", 2020, 24.0, 55.0], ["2025-2026", "Gina Krog – alternativ oljeeksportløysing", 2022, 1.4, 1.5], ["2025-2026", "Ormen Lange fase 3", 2022, 13.8, 12.3], ["2025-2026", "Halten Øst", 2023, 9.9, 10.5], ["2025-2026", "Tyrving", 2023, 6.9, 6.5], ["2025-2026", "Maria fase 2", 2023, 4.5, 4.6]];
const NCS_URLS = {
 "2020-2021": "https://web.archive.org/web/20210613081925id_/https://www.regjeringen.no/contentassets/83dae08df4cb4794879c3b93843f62b5/no/pdfs/prp202020210001oeddddpdfs.pdf",
 "2021-2022": "https://web.archive.org/web/20211104154617id_/https://www.regjeringen.no/contentassets/f2da0c393fb24b1cb2ed0254d76521da/no/pdfs/prp202120220001oeddddpdfs.pdf",
 "2022-2023": "https://web.archive.org/web/20221013165821id_/https://www.regjeringen.no/contentassets/eed690487f7846c397104ee2b984de1b/no/pdfs/prp202220230001oeddddpdfs.pdf",
 "2023-2024": "https://web.archive.org/web/20240509193552id_/https://www.regjeringen.no/contentassets/9fd440e350c8423ab46e9aa67ee674ca/nn-no/pdfs/prp202320240001oeddddpdfs.pdf",
 "2024-2025": "https://web.archive.org/web/20241121024335id_/https://www.regjeringen.no/contentassets/5174ef1ed30a4785bc433090c2db85e1/nn-no/pdfs/prp202420250001_eddddpdfs.pdf",
 "2025-2026": "https://web.archive.org/web/20260909095215id_/https://www.regjeringen.no/contentassets/93021cd505cc48b3996ab338e9e9cb14/nn-no/pdfs/prp202520260001_eddddpdfs.pdf"
};
const DALLAS_NEW = [[2026, "Other U.S. (Shale)", 62.38], [2026, "Eagle Ford", 63], [2026, "Permian Basin – Delaware", 63.21], [2026, "Other U.S. (Nonshale)", 67.82], [2026, "Permian Basin – Midland", 68.57], [2026, "Permian Basin – Other", 70]];

/**
 * Norwegian continental shelf: investment estimate in the plan for development and operation (PDO) against the last estimate
 * reported by the ministry for projects completed in the reporting year. Both numbers of a row are in the fixed prices of
 * the edition year (billion NOK; the editions up to 2022–2023 print million NOK, divided by 1,000 here).
 */
export const NCS_PROJECTS = Object.freeze({
  id: 'ncs-capex-outcomes',
  title: 'Norwegian shelf projects completed 2020–2025: investment estimate at approval against the final estimate',
  source: {
    citation: 'Norwegian Ministry of Petroleum and Energy / Ministry of Energy, Prop. 1 S (budget proposition to the Storting), editions 2020–2021 to 2025–2026, chapter on projects under development on the Norwegian continental shelf, table "Investeringsanslag, prosjekt som er ferdigstilte"; read from the Internet Archive copies of the regjeringen.no documents (the address of each edition is in `urls`)',
    url: NCS_URLS['2025-2026'], licence: 'Norwegian government publication (Norwegian Licence for Open Government Data as generally applied by regjeringen.no; not stated on the pages read)', retrieved: RETRIEVED,
  },
  urls: Object.freeze(NCS_URLS),
  rows: Object.freeze(NCS_ROWS.map((r) => Object.freeze({ edition: r[0], project: r[1], approved: r[2], pdo: r[3], final: r[4] }))),
});

/** UK continental shelf decommissioning: total cost estimate for 2023 onwards as reported by successive annual surveys. */
export const UKCS_DECOM = Object.freeze({
  id: 'ukcs-decommissioning-estimates',
  title: 'UK continental shelf decommissioning cost for 2023 onwards as estimated by the surveys of 2021 to 2025',
  source: {
    citation: 'North Sea Transition Authority, UKCS Decommissioning Cost and Performance Update 2026, figures 2 and 16 (actual and forecast decommissioning costs by survey year, in-year prices and 2025 prices)',
    url: 'https://www.nstauthority.co.uk/media/w1fa3xvv/ukcs-decommissioning-cost-and-performance-update-2026.pdf', licence: 'NSTA copyright; re-use free of charge under the NSTA user agreement', retrieved: RETRIEVED,
  },
  // total = spend since the start of 2023 + remaining forecast, £ billion
  rows: Object.freeze([[2021, 'real', 44.0], [2022, 'real', 45.7], [2023, 'real', 48.8], [2024, 'real', 49.7], [2025, 'real', 50.5], [2021, 'nominal', 36.3], [2022, 'nominal', 39.9], [2023, 'nominal', 45.3], [2024, 'nominal', 48.0], [2025, 'nominal', 50.5]].map((r) => Object.freeze({ survey: r[0], basis: r[1], total: r[2] }))),
  // same report, table 1 and figure 1: actual spend in-year (£ billion) and the abandonment expenditure forecast and out-turn of 2021 (NSTA, UKCS Decommissioning Cost Estimate 2022)
  spend: Object.freeze({ 2023: 2.0, 2024: 2.4, 2025: 2.6 }), abex2021: Object.freeze({ forecast: 1.4, actual: 1.2, url: 'https://www.nstauthority.co.uk/media/8907/decom_cost-estimate-2022_020822_final_v4.pdf' }),
});

/** WTI price needed to drill a new well profitably, by play: mean answer of the latest survey ($/bbl). Earlier years and the ranges of answers are not bundled: the Bank allows reproduction with credit when it is not for private gain, which is narrower than an open licence, so only these cited figures are kept. */
export const DALLAS_BREAKEVEN = Object.freeze({
  id: 'dallasfed-breakeven',
  title: 'Break-even oil price of new wells by US play, Dallas Fed Energy Survey, first quarter of 2026',
  source: {
    citation: 'Federal Reserve Bank of Dallas, Dallas Fed Energy Survey, first quarter 2026, special question "What WTI oil price does your firm need to profitably drill a new well?", chart-data workbook',
    url: 'https://www.dallasfed.org/-/media/Documents/research/surveys/DES/2026/2601/des26q1_charts.xlsx', page: 'https://www.dallasfed.org/research/surveys/des', terms: 'https://www.dallasfed.org/fed/disclaimer',
    licence: 'Federal Reserve Bank of Dallas: reproduction permitted with credit when not distributed for private gain; a few cited figures', retrieved: RETRIEVED,
  },
  newWell: Object.freeze(DALLAS_NEW.map((r) => Object.freeze({ year: r[0], play: r[1], mean: r[2] }))),
});

/** Upstream carbon intensity of crude oil production, reference year 2015 (g CO2-equivalent per MJ of crude): the figures printed in the text of the accepted manuscript. The country table of the paper's supplement was available only in a secondary reproduction without a licence statement and is not bundled. */
export const UPSTREAM_CI = Object.freeze({
  id: 'masnadi-upstream-ci',
  title: 'Upstream carbon intensity of crude oil production, global statistics (2015)',
  source: {
    citation: 'Masnadi, M.S., El-Houjeiri, H.M., Schunack, D., et al., Global carbon intensity of crude oil production, Science 361(6405), 851–853, 2018, doi:10.1126/science.aar6859; accepted manuscript deposited at the U.S. Department of Energy Office of Scientific and Technical Information',
    url: 'https://www.osti.gov/servlets/purl/1485127', licence: 'accepted manuscript in the public repository of the U.S. Department of Energy; seven cited figures', retrieved: RETRIEVED,
  },
  globalMean: 10.3, errorPlus: 16.8, errorMinus: 8.6, // volume-weighted average and its error bar
  countryMin: 3.3, countryMax: 20.3, // range of the country-level volume-weighted averages
  percentiles: Object.freeze({ p5: 4.7, p25: 7.3, p50: 9.1, p75: 11.2, p95: 19.5 }),
});

/**
 * Development cost of subsea tie-backs on the Norwegian shelf that started production in 2010–2024, computed from the open field tables:
 * investments (million NOK of each year) up to and including the year after first production, converted at the annual average exchange
 * rate of each year and moved to 2024 US$ with the US consumer price index, over the original recoverable volume (2025 reserves version).
 * Columns: field, first production, recoverable volume (million Sm³ oil equivalent), oil share, investment (million NOK, nominal),
 * investment (million US$ of 2024), water depth (m; null where the description does not state it).
 */
const TIEBACK_ROWS = [["Morvin", 2010, 13.76, 0.67, 7994, 1943, 360], ["Vega", 2010, 12.07, 0.0, 8503, 2070, 370], ["Trym", 2011, 5.76, 0.0, 2956, 710, 65], ["Atla", 2012, 1.72, 0.19, 1403, 334, 120], ["Gaupe", 2012, 5.14, 0.24, 2329, 572, 90], ["Marulk", 2012, 12.24, 0.06, 3520, 852, 370], ["Oselvar", 2012, 8.56, 0.47, 4433, 1069, 70], ["Visund Sør", 2012, 15.65, 0.24, 3821, 904, 290], ["Hyme", 2013, 4.03, 0.8, 4436, 1049, 250], ["Jette", 2013, 1.98, 0.82, 3562, 837, 127], ["Skuld", 2013, 14.56, 0.92, 9525, 2230, 340], ["Brynhild", 2014, 3.2, 1.0, 8088, 1757, 80], ["Svalin", 2014, 12.08, 1.0, 4664, 1034, 120], ["Bøyla", 2015, 3.71, 0.92, 5396, 1156, 120], ["Flyndre", 2017, 0.48, 0.88, 601, 119, 70], ["Oda", 2019, 7.53, 0.95, 4652, 710, 65], ["Trestakk", 2019, 11.73, 0.92, 6577, 947, 300], ["Utgard", 2019, 4.4, 0.32, 1619, 244, null], ["Dvalin", 2020, 18.79, 0.0, 9589, 1392, null], ["Skogul", 2020, 1.49, 0.9, 2253, 317, 110], ["Duva", 2021, 14.04, 0.26, 5558, 742, 350], ["Nova", 2022, 12.51, 0.73, 11438, 1448, 370], ["Bauge", 2023, 11.64, 0.68, 3465, 475, 280], ["Breidablikk", 2023, 30.46, 1.0, 19573, 2104, 130], ["Fenja", 2023, 15.36, 0.71, 10914, 1469, 325], ["Tommeliten A", 2023, 21.34, 0.27, 11703, 1230, 75], ["Hanz", 2024, 3.07, 0.84, 5398, 574, 115]];
export const NCS_TIEBACKS = Object.freeze({
  id: 'sodir-tieback-development-cost',
  title: 'Development cost per barrel of oil equivalent of Norwegian subsea tie-backs, first production 2010–2024',
  source: {
    citation: 'Norwegian Offshore Directorate, FactPages, tables "Field – Reserves", "Field – Investments (yearly)" and "Field – Description" (fields whose development text names a subsea tie-in to a host), synchronised 8 October 2026; exchange rates: Norges Bank, annual average NOK per US$; cost per barrel computed for this suite',
    url: 'https://factpages.sodir.no/public?/Factpages/external/tableview/field_investment_yearly&rs:Command=Render&rc:Toolbar=false&rc:Parameters=f&IpAddress=not_used&CultureCode=en&rs:Format=CSV&Top100=false', urlReserves: 'https://factpages.sodir.no/public?/Factpages/external/tableview/field_reserves&rs:Command=Render&rc:Toolbar=false&rc:Parameters=f&IpAddress=not_used&CultureCode=en&rs:Format=CSV&Top100=false', urlFx: 'https://data.norges-bank.no/api/data/EXR/A.USD.NOK.SP?format=csv',
    licence: 'Norwegian Licence for Open Government Data (NLOD); derived values computed here', retrieved: '2026-10-09',
  },
  boePerSm3: 6.2898,
  rows: Object.freeze(TIEBACK_ROWS.map((r) => Object.freeze({ field: r[0], first: r[1], oe: r[2], oilShare: r[3], nok: r[4], usd2024: r[5], waterDepth: r[6], perBoe: +(r[5] / (r[2] * 6.2898)).toFixed(1) }))),
});

/** Break-even prices of new supply by segment as published by a market analyst, and the published cost range of the project class of the reference case. */
export const BREAKEVEN_PUBLISHED = Object.freeze({
  source: { citation: 'Rystad Energy, news release "Shale project economics still reign supreme as cost of new oil production rises further", 2024 (average break-even Brent price of not-yet-producing fields by supply segment)', url: 'https://www.rystadenergy.com/news/upstream-breakeven-shale-oil-inflation', licence: 'copyright Rystad Energy; six cited figures', retrieved: '2026-10-09' },
  year: 2024, nonOpecAverage: 47,
  segments: Object.freeze([{ name: 'Onshore Middle East', value: 27 }, { name: 'Offshore shelf', value: 37 }, { name: 'Offshore deep water', value: 43 }, { name: 'North American shale', value: 45 }, { name: 'Oil sands', value: 57 }].map(Object.freeze)),
  tiebackClass: Object.freeze({ low: 100, high: 1500, year: 2015, what: 'total project cost of the subsea tie-backs selected for most Miocene fields of the deep-water Gulf of Mexico', citation: 'U.S. Energy Information Administration (prepared by IHS Global Inc.), Trends in U.S. Oil and Natural Gas Upstream Costs, March 2016, section IX', url: 'https://www.eia.gov/analysis/studies/drilling/pdf/upstream.pdf' }),
});

/** UK continental shelf emission intensity as printed in the regulator's monitoring reports. */
export const UKCS_INTENSITY = Object.freeze({
  source: { citation: 'North Sea Transition Authority, Emissions Monitoring Report 2023', url: 'https://www.nstauthority.co.uk/media/bicn5tva/nsta-emissions-monitoring-report-2023-final-accessible.pdf', licence: 'NSTA copyright; cited data points', retrieved: RETRIEVED },
  year: 2022, total: 28.8, offshoreCO2: 19.8, // kgCO2e/boe (offshore facilities and terminals); kgCO2/boe (offshore facilities only)
});

/** Worked examples of an open textbook and published results of a numerical library, used as reference cases of the financial formulas. */
export const FINANCE_CASES = Object.freeze({
  id: 'finance-reference-cases',
  title: 'Published worked examples: net present value, annuity, terminal value and Black–Scholes–Merton call prices',
  source: {
    citation: 'Dahlquist, J., Knight, R., et al., Principles of Finance, OpenStax, Rice University, 2022, sections 8.2, 16.2 and 16.4 (worked examples); Numerical Algorithms Group, NAG Library Manual Mark 27, routine S30AAF, example program results',
    url: 'https://openstax.org/books/principles-finance/pages/16-2-net-present-value-npv-method', urlOptions: 'https://support.nag.com/numeric/nl/nagdoc_27/examples/baseresults/s30aafe.r.html', licence: 'OpenStax: Creative Commons Attribution-NonCommercial-ShareAlike 4.0 — two worked results cited as facts, no text reproduced; NAG documentation: six values cited', retrieved: RETRIEVED,
  },
  rows: Object.freeze([
    { case: 'NPV of −16,000; 2,000; 4,000; 5,000 × 4 at 9 %', kind: 'npv', value: 2835.63 }, { case: 'Terminal value of the same inflows at 9 %, end of year 6', kind: 'tv', value: 31595.22 }, { case: 'Future value of 1,000 a year for 5 years at 7 %', kind: 'fva', value: 5750.74 },
    { case: 'Call: S 55, K 58, T 0.7, r 10 %, σ 30 %', kind: 'bs', K: 58, T: 0.7, value: 5.9198 }, { case: 'Call: S 55, K 58, T 0.8', kind: 'bs', K: 58, T: 0.8, value: 6.5506 }, { case: 'Call: S 55, K 60, T 0.7', kind: 'bs', K: 60, T: 0.7, value: 5.0809 },
    { case: 'Call: S 55, K 60, T 0.8', kind: 'bs', K: 60, T: 0.8, value: 5.6992 }, { case: 'Call: S 55, K 62, T 0.7', kind: 'bs', K: 62, T: 0.7, value: 4.3389 }, { case: 'Call: S 55, K 62, T 0.8', kind: 'bs', K: 62, T: 0.8, value: 4.9379 },
  ].map(Object.freeze)),
  irrRounded: 0.14, mirrRounded: 0.12, // as printed for the same project (sections 16.3 and 16.4)
});

/** Default rows of the predicted-against-actual table: the Norwegian project outcomes (billion NOK of the edition year) and the UK decommissioning estimate of 2021 against 2025 (£ billion, 2025 prices). */
export const REF_SETS = Object.freeze({
  records: Object.freeze([
    ...NCS_ROWS.map((r) => Object.freeze({ type: 'capex', year: +r[0].slice(5), region: '', actual: r[4], predicted: r[3], note: `${r[1]} (approved ${r[2]}), bn NOK` })),
    Object.freeze({ type: 'abandonment', year: 2025, region: '', actual: 50.5, predicted: 44.0, note: 'UK shelf, 2021 survey against 2025, £bn' }),
  ]),
});
