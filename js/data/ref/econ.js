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
const DALLAS_NEW = [[2017, "Permian (Midland)", 46, 25, 65, 13], [2017, "SCOOP/STACK", 46.88, 35, 75, 8], [2017, "Eagle Ford", 47.5, 40, 55, 4], [2017, "Permian (Delaware)", 48, 30, 60, 10], [2017, "Permian (Central Platform)", 49.62, 35, 65, 13], [2017, "Other U.S. (non-shale)", 52.53, 20, 100, 40], [2017, "Other U.S. (Shale)", 55.25, 45, 65, 8], [2018, "Permian Basin – Midland", 47.33, 20, 70, 15], [2018, "Permian Basin – Delaware", 48.92, 38, 70, 13], [2018, "Bakken", 50, 40, 60, 4], [2018, "Permian Basin – Other", 52.33, 40, 75, 18], [2018, "Oklahoma – SCOOP/STACK", 52.86, 38, 65, 7], [2018, "Other U.S. (Shale)", 53.57, 30, 70, 7], [2018, "Other U.S. (Non-shale)", 54.76, 20, 75, 34], [2019, "Permian Basin – Midland", 47.5, 23, 65, 17], [2019, "Other U.S. (Shale)", 48.5, 35, 60, 12], [2019, "Permian Basin – Delaware", 48.69, 40, 65, 13], [2019, "Other U.S. (Non-shale)", 49.41, 20, 75, 45], [2019, "Eagle Ford", 51.27, 40, 75, 11], [2019, "Oklahoma – SCOOP/STACK", 52.6, 48, 60, 5], [2019, "Permian Basin – Other", 53.63, 40, 70, 20], [2020, "Permian Basin – Midland", 45.57, 30, 60, 22], [2020, "Eagle Ford", 46, 40, 55, 7], [2020, "Permian Basin – Other", 50, 30, 70, 23], [2020, "Other U.S. (Non-shale)", 50.34, 15, 70, 44], [2020, "Bakken", 51, 40, 60, 5], [2020, "Other U.S. (Shale)", 51, 45, 65, 8], [2020, "Permian Basin – Delaware", 51.56, 35, 70, 18], [2021, "Permian Basin – Midland", 45.8, 30, 75, 20], [2021, "Eagle Ford", 46.25, 32.5, 65, 6], [2021, "Permian Basin – Delaware", 49.25, 30, 75, 18], [2021, "Other U.S. (Nonshale)", 53.17, 20, 100, 50], [2021, "Permian Basin – Other", 53.17, 35, 100, 23], [2021, "Other U.S. (Shale)", 58.06, 35, 100, 18], [2022, "Eagle Ford", 48, 40, 70, 5], [2022, "Permian Basin – Delaware", 49.74, 35, 70, 19], [2022, "Permian Basin – Midland", 51.29, 30, 85, 17], [2022, "Permian Basin – Other", 54.13, 30, 85, 23], [2022, "Other U.S. (Nonshale)", 60.26, 27.5, 150, 47], [2022, "Other U.S. (Shale)", 68.7, 15, 150, 10], [2023, "Eagle Ford", 56.43, 45, 70, 7], [2023, "Permian Basin – Midland", 58.33, 40, 80, 18], [2023, "Permian Basin – Delaware", 60.71, 40, 75, 17], [2023, "Other U.S. (Shale)", 60.91, 40, 75, 11], [2023, "Other U.S. (Non-shale)", 63.28, 35, 100, 47], [2023, "Permian Basin – Other", 65.81, 48, 90, 16], [2024, "Other U.S. (Shale)", 59.41, 35, 75, 17], [2024, "Permian Basin – Midland", 61.61, 40, 85, 23], [2024, "Permian Basin – Delaware", 63.81, 30, 95, 21], [2024, "Oklahoma – SCOOP/STACK", 65.2, 36, 75, 5], [2024, "Other U.S. (Nonshale)", 65.72, 30, 95, 45], [2024, "Permian Basin – Other", 70, 45, 90, 18], [2025, "Permian Basin – Midland", 60.56, 45, 75, 9], [2025, "Eagle Ford", 61.67, 45, 80, 6], [2025, "Permian Basin – Delaware", 62.22, 45, 90, 18], [2025, "Other U.S. (Shale)", 63, 35, 80, 16], [2025, "Other U.S. (Nonshale)", 65.77, 30, 90, 48], [2025, "Permian Basin – Other", 69.75, 50, 90, 20], [2026, "Other U.S. (Shale)", 62.38, 35, 85, 17], [2026, "Eagle Ford", 63, 50, 70, 5], [2026, "Permian Basin – Delaware", 63.21, 50, 85, 19], [2026, "Other U.S. (Nonshale)", 67.82, 25, 100, 53], [2026, "Permian Basin – Midland", 68.57, 50, 90, 14], [2026, "Permian Basin – Other", 70, 50, 85, 12]];
const DALLAS_OPEX = [[2025, "Eagle Ford", 25.83, 10, 50, 6], [2025, "Permian Basin – Delaware", 33.03, 9, 70, 18], [2025, "Permian Basin – Midland", 35.44, 14, 70, 9], [2025, "Other U.S. (Shale)", 40.81, 6, 75, 16], [2025, "Other U.S. (Nonshale)", 44.62, 10, 75, 50], [2025, "Permian Basin – Other", 45.15, 15, 75, 20], [2026, "Permian Basin – Delaware", 33.82, 3.5, 70, 19], [2026, "Eagle Ford", 40, 15, 55, 5], [2026, "Permian Basin – Midland", 41.54, 15, 65, 13], [2026, "Permian Basin – Other", 43.63, 3.5, 70, 12], [2026, "Other U.S. (Shale)", 44.12, 15, 75, 17], [2026, "Other U.S. (Nonshale)", 47, 10, 75, 52]];
const CI_ROWS = [["Syria", 29.8], ["Democratic Republic of Congo", 29.2], ["Uzbekistan", 27.4], ["Yemen", 26.9], ["Albania", 23.7], ["Algeria", 20.3], ["Venezuela", 20.3], ["Myanmar", 20.2], ["Cameroon", 18.4], ["Canada", 17.6], ["Iran", 17.1], ["Turkmenistan", 15.9], ["Tunisia", 15.4], ["Indonesia", 15.3], ["Georgia", 15.2], ["Sudan", 14.9], ["Mauritania", 14.8], ["Trinidad and Tobago", 14.3], ["Iraq", 14.1], ["Gabon", 13.2], ["Malaysia", 12.9], ["Nigeria", 12.6], ["Pakistan", 12.2], ["Ukraine", 11.8], ["Oman", 11.7], ["Philippines", 11.6], ["Niger", 11.3], ["United States", 11.3], ["Chile", 11.2], ["Libya", 11.0], ["Peru", 10.9], ["Republic of Congo", 10.6], ["Egypt", 10.6], ["Brazil", 10.3], ["Chad", 10.2], ["Mexico", 9.9], ["Guatemala", 9.8], ["Lithuania", 9.7], ["Russian Federation", 9.7], ["Kazakhstan", 9.7], ["Kyrgyzstan", 9.4], ["Tajikistan", 9.4], ["Morocco", 9.3], ["Ecuador", 9.3], ["Barbados", 9.3], ["Argentina", 9.1], ["Australia", 9.1], ["Cuba", 9.0], ["Bolivia", 9.0], ["Latvia", 8.9], ["Vietnam", 8.8], ["Belize", 8.8], ["Bulgaria", 8.6], ["India", 8.6], ["Papua New Guinea", 8.5], ["Turkey", 8.4], ["Colombia", 8.3], ["Afghanistan", 8.3], ["Suriname", 8.2], ["Poland", 8.2], ["New Zealand", 8.2], ["United Kingdom", 7.9], ["Hungary", 7.9], ["Croatia", 7.8], ["Germany", 7.7], ["Japan", 7.7], ["Serbia", 7.7], ["Austria", 7.6], ["France", 7.5], ["Angola", 7.5], ["Romania", 7.4], ["United Arab Emirates", 7.1], ["China", 7.0], ["Kuwait", 6.9], ["Qatar", 6.5], ["Equatorial Guinea", 6.4], ["Jordan", 6.3], ["Azerbaijan", 6.3], ["Cote d'Ivoire", 6.1], ["Italy", 6.1], ["Greece", 5.9], ["Brunei", 5.7], ["Norway", 5.6], ["Ghana", 5.2], ["Thailand", 5.1], ["Bahrain", 5.0], ["Saudi Arabia", 4.6], ["Spain", 4.1], ["Netherlands", 3.9], ["Denmark", 3.3]];

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

/** WTI price needed to drill a new well profitably, by play: survey mean, lowest and highest answer, number of answers ($/bbl). */
export const DALLAS_BREAKEVEN = Object.freeze({
  id: 'dallasfed-breakeven',
  title: 'Break-even oil price of new wells by US play, Dallas Fed Energy Survey, first quarters of 2017–2026',
  source: {
    citation: 'Federal Reserve Bank of Dallas, Dallas Fed Energy Survey, first-quarter special questions 2017–2026 ("What WTI oil price does your firm need to profitably drill a new well?" and "… to cover operating expenses for existing wells?"), chart-data workbooks',
    url: 'https://www.dallasfed.org/-/media/Documents/research/surveys/DES/2026/2601/des26q1_charts.xlsx', page: 'https://www.dallasfed.org/research/surveys/des', licence: 'publicly released survey data of the Federal Reserve Bank of Dallas; cited', retrieved: RETRIEVED,
  },
  newWell: Object.freeze(DALLAS_NEW.map((r) => Object.freeze({ year: r[0], play: r[1], mean: r[2], min: r[3], max: r[4], n: r[5] }))),
  operating: Object.freeze(DALLAS_OPEX.map((r) => Object.freeze({ year: r[0], play: r[1], mean: r[2], min: r[3], max: r[4], n: r[5] }))),
});

/** Volume-weighted upstream carbon intensity of crude oil production by country, reference year 2015 (g CO2-equivalent per MJ of crude). */
export const UPSTREAM_CI = Object.freeze({
  id: 'masnadi-upstream-ci',
  title: 'Upstream carbon intensity of crude oil production by country (2015)',
  source: {
    citation: 'Masnadi, M.S., El-Houjeiri, H.M., Schunack, D., et al., Global carbon intensity of crude oil production, Science 361(6405), 851–853, 2018, doi:10.1126/science.aar6859; country values as reproduced in Science Based Targets initiative, Annex D — Fuel specific calculations for the WTW indicator, 10 August 2020, table D.Oil.3; global average from the accepted manuscript',
    url: 'https://files.sciencebasedtargets.org/production/legacy/2020/08/OG-Annex-D.pdf', urlGlobal: 'https://www.osti.gov/servlets/purl/1485127', licence: 'data points cited from a secondary reproduction of the supplementary table; no licence statement read', retrieved: RETRIEVED,
  },
  globalMean: 10.3, // g CO2-eq/MJ, with an error bar of +16.8 / −8.6
  rows: Object.freeze(CI_ROWS.map((r) => Object.freeze({ country: r[0], ci: r[1] }))),
  // ISO codes of the producing countries offered by the fiscal table and the site page
  codes: Object.freeze({ NO: 'Norway', GB: 'United Kingdom', US: 'United States', BR: 'Brazil', NG: 'Nigeria', AO: 'Angola', GH: 'Ghana', AU: 'Australia', CA: 'Canada', MY: 'Malaysia', ID: 'Indonesia', SA: 'Saudi Arabia', TT: 'Trinidad and Tobago', SR: 'Suriname', NL: 'Netherlands', DK: 'Denmark', KW: 'Kuwait', GA: 'Gabon', GQ: 'Equatorial Guinea', CN: 'China', OM: 'Oman', AR: 'Argentina', DZ: 'Algeria', VE: 'Venezuela', IR: 'Iran', IQ: 'Iraq', RU: 'Russian Federation', KZ: 'Kazakhstan', MX: 'Mexico', EG: 'Egypt', LY: 'Libya', QA: 'Qatar', AE: 'United Arab Emirates', IN: 'India', CO: 'Colombia', TH: 'Thailand', VN: 'Vietnam', AZ: 'Azerbaijan', EC: 'Ecuador', CG: 'Republic of Congo', BN: 'Brunei', IT: 'Italy', TN: 'Tunisia', CM: 'Cameroon', TM: 'Turkmenistan', PK: 'Pakistan', PE: 'Peru', TR: 'Turkey', RO: 'Romania', DE: 'Germany', FR: 'France', NZ: 'New Zealand', PG: 'Papua New Guinea', CI: "Cote d'Ivoire", TD: 'Chad', SD: 'Sudan', YE: 'Yemen', SY: 'Syria', BH: 'Bahrain', JP: 'Japan' }),
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
    url: 'https://openstax.org/books/principles-finance/pages/16-2-net-present-value-npv-method', urlOptions: 'https://support.nag.com/numeric/nl/nagdoc_27/examples/baseresults/s30aafe.r.html', licence: 'CC BY-NC-SA 4.0 (OpenStax); NAG documentation, values cited', retrieved: RETRIEVED,
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
