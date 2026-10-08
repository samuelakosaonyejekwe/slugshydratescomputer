// Live global data connectors. Every request goes straight from the user's browser to a public,
// key-less HTTPS service on the allow-list below; nothing passes through any server of this app, so the
// data stay fresh wherever the app is opened. Responses are treated as untrusted numbers/text only.
import { clamp, mean, quantile } from './num.js';
import { ATLAS } from '../data/atlas.js';
import { PRICES, MARKETS, COSTS, CORPORATE_TAX, ELECTRICITY, RATES, ISO3 } from '../data/prices.js';
import { fiscalOf } from '../data/fiscal.js';

export const SOURCES = [
  { id: 'place', name: 'Place and country', host: 'api.bigdatacloud.net', provider: 'BigDataCloud reverse geocoding', gives: 'Locality, country' },
  { id: 'weather', name: 'Weather and solar', host: 'api.open-meteo.com', provider: 'Open-Meteo forecast (national weather models)', gives: 'Air temperature, wind, humidity, pressure, solar irradiance, land elevation' },
  { id: 'marine', name: 'Sea state, currents and tides', host: 'marine-api.open-meteo.com', provider: 'Open-Meteo marine (wave and ocean models)', gives: 'Sea-surface temperature, waves, ocean currents, sea-level / tide series' },
  { id: 'bathy', name: 'Seabed and terrain', host: 'gis.ngdc.noaa.gov', provider: 'NOAA NCEI global DEM mosaic (best available resolution); fallback SRTM30+ via PacIOOS', gives: 'Bathymetry / topography grid, water depth, seabed slope; seabed temperature from the bundled World Ocean Atlas 2023 profile' },
  { id: 'salinity', name: 'Sea temperature and salinity', host: 'erddap.emodnet-physics.eu', provider: 'SeaDataCloud global T–S climatology (EMODnet Physics ERDDAP)', gives: 'Monthly near-surface salinity and temperature climatology' },
  { id: 'economy', name: 'Inflation and interest', host: 'api.worldbank.org', provider: 'World Bank Open Data', gives: 'Consumer-price inflation, lending interest rate' },
  { id: 'fx', name: 'Currency', host: 'open.er-api.com', provider: 'Open exchange-rate API', gives: 'Local currency per US dollar' },
  { id: 'energy', name: 'Grid carbon, renewables and carbon price', host: 'ourworldindata.org', provider: 'Our World in Data (Ember / Energy Institute / World Bank carbon-pricing series)', gives: 'Carbon intensity of electricity, renewable share of generation, national carbon price' },
  { id: 'prices', name: 'Oil and gas prices', host: 'raw.githubusercontent.com', provider: 'US EIA daily spot series (Brent, WTI, Henry Hub) from the open “datasets” collection; fallback: the EIA open-data service', gives: 'Latest Brent and WTI crude and Henry Hub gas prices, 30-day mean and one-year range' },
  { id: 'climate', name: 'Solar and wind climatology', host: 'power.larc.nasa.gov', provider: 'NASA POWER long-term climatology', gives: 'Monthly and annual solar irradiation, wind speed, air temperature' },
  { id: 'fiscal', name: 'Corporate tax', host: 'sdmx.oecd.org', provider: 'OECD Tax Database, statutory corporate income tax rates (146 jurisdictions); fallback: Tax Foundation worldwide corporate tax rates', gives: 'Statutory corporate income tax rate and its year; the petroleum fiscal terms (royalty, petroleum tax, headline rate) come from the sourced table bundled with the app' },
  { id: 'costs', name: 'Cost escalation', host: 'api.bls.gov', provider: 'US Bureau of Labor Statistics public data API (producer and consumer price indexes)', gives: 'Oil and gas field machinery and equipment price index (20 years, monthly) and its 5-year escalation, steel mill products index, US consumer price index' },
  { id: 'markets', name: 'Regional gas and iron ore prices', host: 'api.imf.org', provider: 'IMF Primary Commodity Price System (monthly); fallback for European gas: Our World in Data (Energy Institute, annual)', gives: 'European (TTF) and Asian (Japan LNG) natural-gas prices, iron ore price as the steel-cost proxy' },
  { id: 'rates', name: 'Bond yield and policy rate', host: 'api.imf.org', provider: 'IMF Monetary and Financial Statistics, interest rates; fallbacks: OECD long-term interest rates, BIS central-bank policy rates', gives: 'Government bond yield, Treasury bill yield, central-bank policy rate' },
  { id: 'power', name: 'Industrial electricity price', host: 'ec.europa.eu', provider: 'Eurostat electricity prices for non-household consumers (Europe); US EIA average industrial retail price (United States)', gives: 'Electricity price for industry, where an open series exists' },
];
export const EVIDENCE_SOURCES = [
  { name: 'OpenAlex', host: 'api.openalex.org', gives: 'Open index of the global research literature' },
  { name: 'Crossref', host: 'api.crossref.org', gives: 'DOI registry of scholarly and technical publications' },
];
const GEOCODE_HOST = 'geocoding-api.open-meteo.com';
const ALLOWED = new Set([...SOURCES.map((s) => s.host), ...EVIDENCE_SOURCES.map((s) => s.host), GEOCODE_HOST, 'pae-paha.pacioos.hawaii.edu', 'api.eia.gov', 'stats.bis.org', 'api.db.nomics.world']);

/** GET text (CSV) from an allow-listed HTTPS host with a timeout and a size cap. */
export async function getText(url, ms = 20000, maxChars = 400000) {
  const u = new URL(url);
  if (u.protocol !== 'https:' || !ALLOWED.has(u.hostname)) throw new Error('Blocked request to a host that is not on the allow-list.');
  const ctl = new AbortController(), timer = setTimeout(() => ctl.abort(), ms);
  try { const r = await fetch(u.href, { signal: ctl.signal, credentials: 'omit', referrerPolicy: 'no-referrer', cache: 'no-store' }); if (!r.ok) throw new Error('HTTP ' + r.status); return (await r.text()).slice(0, maxChars); }
  finally { clearTimeout(timer); }
}

/** GET JSON from an allow-listed HTTPS host with a timeout. */
export async function getJSON(url, ms = 20000, form = null) {
  const u = new URL(url);
  if (u.protocol !== 'https:' || !ALLOWED.has(u.hostname)) throw new Error('Blocked request to a host that is not on the allow-list.');
  const ctl = new AbortController(), timer = setTimeout(() => ctl.abort(), ms);
  try {
    const opt = { signal: ctl.signal, credentials: 'omit', referrerPolicy: 'no-referrer', cache: 'no-store' };
    const r = await fetch(u.href, form ? { ...opt, method: 'POST', body: new URLSearchParams(form) } : opt);
    if (!r.ok) throw new Error('HTTP ' + r.status);
    return await r.json();
  } finally { clearTimeout(timer); }
}
const num = (x) => (typeof x === 'number' && Number.isFinite(x) ? x : null);
const txt = (x, n = 80) => (typeof x === 'string' ? x.slice(0, n) : '');

/** Place search by name -> [{ name, country, lat, lon }]. */
export async function searchPlace(q) {
  const j = await getJSON(`https://${GEOCODE_HOST}/v1/search?name=${encodeURIComponent(String(q).slice(0, 80))}&count=8&language=en&format=json`);
  return (j.results || []).map((r) => ({ name: txt(r.name), admin: txt(r.admin1), country: txt(r.country), code: txt(r.country_code, 3), lat: num(r.latitude), lon: num(r.longitude) })).filter((r) => r.lat !== null && r.lon !== null);
}

// Stand-by industrial electricity tariff ($/kWh) and grid carbon intensity (kgCO₂/kWh) by country: hand-entered
// ESTIMATES, used only where no sourced figure exists and always labelled as estimates in the site data. The tariff
// of every country with an open industrial series (Eurostat for Europe, US EIA for the United States — live, or the
// dated snapshot in js/data/prices.js) comes from that series, so those countries carry null here. For the others
// no open, cross-origin-readable industrial tariff series exists; the World Bank business tariff (2019, small
// commercial connection) is used for countries missing from this table. The grid-carbon figure is replaced by the
// Our World in Data series (live, or the copy in the built-in atlas) whenever that is available.
// CURRENCY is a quick table for the common cases; the built-in atlas carries the full ISO 4217 list.
const ENERGY = {
  AE: [0.08, 0.40], SA: [0.05, 0.57], QA: [0.04, 0.49], KW: [0.03, 0.57], BH: [0.07, 0.49], OM: [0.06, 0.44], IL: [0.11, 0.45], EG: [0.05, 0.43], MA: [0.10, 0.62], DZ: [0.04, 0.48], TN: [0.09, 0.47], LY: [0.03, 0.55],
  ES: [null, 0.17], IT: [null, 0.30], GR: [null, 0.34], CY: [null, 0.60], MT: [null, 0.39], TR: [null, 0.42], PT: [null, 0.15], FR: [null, 0.06], GB: [0.24, 0.21], DE: [null, 0.36], NL: [null, 0.27],
  US: [null, 0.37], MX: [0.11, 0.42], CL: [0.13, 0.30], PE: [0.08, 0.20], BR: [0.13, 0.10], AR: [0.07, 0.31], CA: [0.09, 0.12],
  AU: [0.14, 0.55], NZ: [0.11, 0.11], CN: [0.09, 0.56], IN: [0.10, 0.71], PK: [0.14, 0.40], BD: [0.09, 0.57], SG: [0.17, 0.41], JP: [0.19, 0.46], KR: [0.11, 0.43], ID: [0.07, 0.68], MY: [0.09, 0.59], TH: [0.12, 0.47], VN: [0.08, 0.47], PH: [0.15, 0.61],
  ZA: [0.09, 0.71], NG: [0.07, 0.37], GH: [0.13, 0.30], KE: [0.16, 0.09], NA: [0.11, 0.06], TZ: [0.10, 0.34], SN: [0.17, 0.53], DJ: [0.25, 0.50], CV: [0.28, 0.55], IR: [0.02, 0.49], IQ: [0.06, 0.60], JO: [0.11, 0.39], YE: [0.15, 0.60],
};
const CURRENCY = {
  AE: 'AED', SA: 'SAR', QA: 'QAR', KW: 'KWD', BH: 'BHD', OM: 'OMR', IL: 'ILS', EG: 'EGP', MA: 'MAD', DZ: 'DZD', TN: 'TND', LY: 'LYD', ES: 'EUR', IT: 'EUR', GR: 'EUR', CY: 'EUR', MT: 'EUR', PT: 'EUR', FR: 'EUR', DE: 'EUR', NL: 'EUR', TR: 'TRY', GB: 'GBP',
  US: 'USD', MX: 'MXN', CL: 'CLP', PE: 'PEN', BR: 'BRL', AR: 'ARS', CA: 'CAD', AU: 'AUD', NZ: 'NZD', CN: 'CNY', IN: 'INR', PK: 'PKR', BD: 'BDT', SG: 'SGD', JP: 'JPY', KR: 'KRW', ID: 'IDR', MY: 'MYR', TH: 'THB', VN: 'VND', PH: 'PHP',
  ZA: 'ZAR', NG: 'NGN', GH: 'GHS', KE: 'KES', NA: 'NAD', TZ: 'TZS', SN: 'XOF', DJ: 'DJF', CV: 'CVE', IR: 'IRR', IQ: 'IQD', JO: 'JOD', YE: 'YER',
};

// ---- economic feeds ------------------------------------------------------------------------------------
// One function per open service; each returns plain numbers keyed by country or series. The connectors below
// ask for one country, tools/prices.mjs asks for all of them to rewrite the bundled snapshot (js/data/prices.js),
// so the live values and the offline fallback always come from the same series.
const ym = (p) => txt(String(p ?? ''), 12).replace(/^(\d{4})-M(\d\d)$/, '$1-$2'); // "2026-M09" -> "2026-09"
const csvRows = (csv) => { const [head, ...lines] = csv.trim().split(/\r?\n/).map((l) => l.split(',')); return lines.map((c) => Object.fromEntries(head.map((h, i) => [h.trim(), (c[i] ?? '').trim()]))); };
const BLS_COST = 'PCU333132333132', BLS_STEEL = 'WPU1017', BLS_CPI = 'CUUR0000SA0';
const IMF_RATES = { S13BOND_RT_PT_A_PT: 'bondYield', GSTBILY_RT_PT_A_PT: 'treasuryBillYield', MFS166_RT_PT_A_PT: 'policyRate' };
const IMF_COMMODITIES = { PNGASEU: 'gasEurope', PNGASJP: 'gasAsia', PIORECR: 'ironOre' };
/** Countries reported in the Eurostat non-household electricity price table (ISO-2 -> Eurostat code). */
const EUROSTAT_GEO = Object.fromEntries('BE BG CZ DK DE EE IE ES FR HR IT CY LV LT LU HU MT NL AT PL PT RO SI SK FI SE IS LI NO BA ME MD MK GE AL RS TR UA XK'.split(' ').map((c) => [c, c]).concat([['GR', 'EL']]));
export const SOURCE_NAMES = {
  oecdCit: 'OECD Tax Database, combined statutory corporate income tax rate', taxFoundation: 'Tax Foundation, Worldwide Corporate Tax Rates', bls: 'US Bureau of Labor Statistics', imfPcps: 'IMF Primary Commodity Price System',
  imfMfs: 'IMF Monetary and Financial Statistics', oecdRates: 'OECD long-term interest rates (10-year government bonds)', bis: 'BIS central-bank policy rates', eurostat: 'Eurostat, electricity prices for non-household consumers (nrg_pc_205)',
  eia: 'US EIA, average retail price of electricity, industrial sector', owidGas: 'Our World in Data / Energy Institute Statistical Review, annual average',
};
export const feeds = {
  /** US BLS public API (one request for several series): { seriesId: [[yyyy-mm, value], …] }, oldest first. */
  async bls(ids, y0, y1) {
    const j = await getJSON('https://api.bls.gov/publicAPI/v2/timeseries/data/', 20000, { seriesid: ids.join(','), startyear: String(y0), endyear: String(y1) });
    if (j?.status !== 'REQUEST_SUCCEEDED') throw new Error(txt(String(j?.message?.[0] || 'Request not processed'), 90));
    const out = {};
    for (const s of j.Results?.series || []) out[txt(s.seriesID, 30)] = (s.data || []).filter((r) => /^M(0[1-9]|1[0-2])$/.test(r.period) && /^\d{4}$/.test(String(r.year)) && Number.isFinite(parseFloat(r.value))).map((r) => [`${r.year}-${r.period.slice(1)}`, parseFloat(r.value)]).sort((a, b) => (a[0] < b[0] ? -1 : 1));
    return out;
  },
  /** The same BLS series from the DBnomics mirror (complete history, refreshed less often than the BLS service itself). */
  async blsMirror() {
    const one = async (path) => { const j = await getJSON(`https://api.db.nomics.world/v22/series/BLS/${path}?observations=1&metadata=false`, 25000), d = j?.series?.docs?.[0] || {}; return (d.period || []).map((t, i) => [txt(String(t), 7), parseFloat(d.value?.[i])]).filter((r) => /^\d{4}-\d\d$/.test(r[0]) && Number.isFinite(r[1])); };
    const [cost, steel, cpi] = await Promise.all([one('pc/' + BLS_COST), one('wp/' + BLS_STEEL), one('cu/' + BLS_CPI)]);
    if (cost.length < 24) throw new Error('No index series returned');
    return { [BLS_COST]: cost, [BLS_STEEL]: steel, [BLS_CPI]: cpi };
  },
  /**
   * Twenty years of the three price indexes: { cost, steel, cpi } as [[yyyy-mm, value], …]. The key-less BLS service
   * answers ten years per request and a limited number of requests a day, so the older decade comes from a second
   * request or, failing that, from the mirror.
   */
  async indexes() {
    const y = new Date().getUTCFullYear(), ids = [BLS_COST, BLS_STEEL, BLS_CPI];
    const [a, b] = await Promise.allSettled([feeds.bls(ids, y - 9, y), feeds.bls(ids, y - 19, y - 10)]);
    if (a.status !== 'fulfilled') throw a.reason;
    const old = b.status === 'fulfilled' ? b.value : await feeds.blsMirror().catch(() => ({}));
    const join = (id) => { const rec = a.value[id] || [], from = rec[0]?.[0] || '9999', lo = `${y - 19}-01`; return [...(old[id] || []).filter((r) => r[0] < from && r[0] >= lo), ...rec]; };
    return { cost: join(BLS_COST), steel: join(BLS_STEEL), cpi: join(BLS_CPI) };
  },
  /** IMF data service (SDMX 2.1, XML): [{ country, indicator, obs: [[period, value], …] }]. */
  async imf(flow, key, query) {
    const xml = await getText(`https://api.imf.org/external/sdmx/2.1/data/${flow}/${key}?${query}`, 25000, 4e6), out = [];
    const attr = (s, k) => (new RegExp(`\\s${k}="([^"]*)"`).exec(s) || [])[1] || '';
    for (const m of xml.matchAll(/<Series(\s[^>]*[^/>])>([\s\S]*?)<\/Series>/g)) {
      const obs = [...m[2].matchAll(/<Obs(\s[^>]*)>/g)].map((o) => [ym(attr(o[1], 'TIME_PERIOD')), parseFloat(attr(o[1], 'OBS_VALUE'))]).filter((o) => /^\d{4}(-\d\d)?$/.test(o[0]) && Number.isFinite(o[1])).sort((a, b) => (a[0] < b[0] ? -1 : 1));
      if (obs.length) out.push({ country: txt(attr(m[1], 'COUNTRY'), 8), indicator: txt(attr(m[1], 'INDICATOR'), 40), obs });
    }
    return out;
  },
  /** Monthly benchmark prices of the last ten years: { gasEurope, gasAsia, ironOre } as [[yyyy-mm, value], …] (US$/MMBtu, US$/MMBtu, US$/t). */
  async commodities() {
    const out = {};
    for (const s of await feeds.imf('PCPS', `G001.${Object.keys(IMF_COMMODITIES).join('+')}.USD.M`, `startPeriod=${new Date().getUTCFullYear() - 10}-M01`)) if (IMF_COMMODITIES[s.indicator]) out[IMF_COMMODITIES[s.indicator]] = s.obs;
    if (!Object.keys(out).length) throw new Error('No commodity series returned');
    return out;
  },
  /** Latest government bond yield, Treasury bill yield and policy rate: { ISO3: { bondYield: [value, yyyy-mm], … } } (iso3 '' = every country). */
  async imfRates(iso3 = '') {
    const out = {};
    for (const s of await feeds.imf('MFS_IR', `${iso3}.${Object.keys(IMF_RATES).join('+')}.M`, 'lastNObservations=1')) { const k = IMF_RATES[s.indicator], o = s.obs[s.obs.length - 1]; if (k && s.country && Math.abs(o[1]) < 200) (out[s.country] ??= {})[k] = [+o[1].toFixed(3), o[0]]; }
    return out;
  },
  /** OECD data service, CSV rows as objects. */
  async oecd(flow, key, query = 'lastNObservations=1') { return csvRows(await getText(`https://sdmx.oecd.org/public/rest/data/${flow}/${key}?${query}&format=csvfile`, 25000, 2e6)); },
  /** Combined (central + sub-central) statutory corporate income tax rate: { ISO3: [rate %, year] } (iso3 '' = every jurisdiction). */
  async corporateTax(iso3 = '') {
    const out = {};
    for (const r of await feeds.oecd('OECD.CTP.TPS,DSD_TAX_CIT@DF_CIT,', `${iso3}.A.CIT_C.ST.....`)) { const v = parseFloat(r.OBS_VALUE), y = parseInt(r.TIME_PERIOD, 10); if (/^[A-Z]{3}$/.test(r.REF_AREA) && v >= 0 && v <= 90 && y > 1990 && !(out[r.REF_AREA]?.[1] > y)) out[r.REF_AREA] = [v, y]; }
    return out;
  },
  /** Tax Foundation table (2023 edition): [rate %, year] of one country (ISO-2) or null. */
  async corporateTaxTF(code) {
    const csv = await getText('https://raw.githubusercontent.com/TaxFoundation/worldwide-corporate-tax-rates/master/final_data/final_data_long.csv', 25000, 3e6); let hit = null;
    for (const line of csv.split('\n')) { const c = line.split(',').map((x) => x.replace(/^"|"$/g, '')); const v = parseFloat(c[c.length - 7]), y = parseInt(c[c.length - 8], 10); if (c[0] === code && v >= 0 && v <= 90 && y > 1990 && !(hit?.[1] > y)) hit = [v, y]; }
    return hit;
  },
  /** OECD long-term (10-year government bond) interest rate: { ISO3: [value %, yyyy-mm] }. */
  async oecdBond(iso3 = '') {
    const out = {};
    for (const r of await feeds.oecd('OECD.SDD.STES,DSD_STES@DF_FINMARK,4.0', `${iso3}.M.IRLT.PA.....`)) { const v = parseFloat(r.OBS_VALUE); if (/^[A-Z]{3}$/.test(r.REF_AREA) && /^\d{4}-\d\d$/.test(r.TIME_PERIOD) && Math.abs(v) < 200 && !(out[r.REF_AREA]?.[1] > r.TIME_PERIOD)) out[r.REF_AREA] = [+v.toFixed(3), r.TIME_PERIOD]; }
    return out;
  },
  /** BIS central-bank policy rate of one country (ISO-2): [value %, yyyy-mm] or null. */
  async bisPolicy(code) {
    const csv = await getText(`https://stats.bis.org/api/v2/data/dataflow/BIS/WS_CBPOL/1.0/M.${encodeURIComponent(code)}?lastNObservations=1&format=csv`, 20000, 20000);
    const m = /,(\d{4}-\d\d),(-?\d+(?:\.\d+)?),/.exec(csv.split(/\r?\n/)[1] || '');
    return m && Math.abs(+m[2]) < 200 ? [+m[2], m[1]] : null;
  },
  /** Eurostat half-yearly electricity price for industry (band ID: 2 000–19 999 MWh a year, excluding VAT and other recoverable taxes), EUR/kWh: { ISO-2: [value, period] } (code '' = every reporting country). */
  async eurostatPower(code = '') {
    if (code && !EUROSTAT_GEO[code]) return {};
    const j = await getJSON(`https://ec.europa.eu/eurostat/api/dissemination/statistics/1.0/data/nrg_pc_205?format=JSON&lang=EN&nrg_cons=MWH2000-19999&tax=X_VAT&currency=EUR&lastTimePeriod=4${code ? '&geo=' + EUROSTAT_GEO[code] : ''}`, 25000);
    const gi = j?.dimension?.geo?.category?.index || {}, ti = j?.dimension?.time?.category?.index || {}, nt = Object.keys(ti).length, back = Object.fromEntries(Object.entries(EUROSTAT_GEO).map(([a, b]) => [b, a])), out = {};
    if (!nt || (j.size || []).reduce((a, b) => a * b, 1) !== Object.keys(gi).length * nt) throw new Error('Unexpected table layout');
    for (const [g, i] of Object.entries(gi)) for (const [t, k] of Object.entries(ti)) { const v = num(j.value?.[String(i * nt + k)]); if (back[g] && v !== null && v > 0 && v < 2 && !(out[back[g]]?.[1] > t)) out[back[g]] = [v, txt(t, 8)]; }
    return out;
  },
  /** US average industrial retail electricity price, mean of the last twelve published months: [US$/kWh, 'yyyy-mm' of the latest month]. */
  async eiaPower() {
    const j = await getJSON('https://api.eia.gov/v2/electricity/retail-sales/data/?api_key=DEMO_KEY&frequency=monthly&data[0]=price&facets[sectorid][]=IND&facets[stateid][]=US&sort[0][column]=period&sort[0][direction]=desc&length=12');
    const rows = (j?.response?.data || []).map((r) => [txt(r.period, 7), parseFloat(r.price)]).filter((r) => /^\d{4}-\d\d$/.test(r[0]) && r[1] > 0 && r[1] < 100);
    if (rows.length < 6) throw new Error('No industrial price returned');
    return [+(mean(rows.map((r) => r[1])) / 100).toFixed(4), rows[0][0]];
  },
  /** Annual average Netherlands TTF gas price (Energy Institute via Our World in Data), converted from US$/MWh: [US$/MMBtu, year]. */
  async owidGasEurope() {
    const csv = await getText('https://ourworldindata.org/grapher/natural-gas-prices.csv?csvType=full'); let hit = null;
    for (const line of csv.split('\n')) { const c = line.split(','), v = parseFloat(c[c.length - 1]), y = parseInt(c[c.length - 2], 10); if (c[0] === 'Netherlands TTF' && v > 0 && y > 2000 && !(hit?.[1] > y)) hit = [+(v / 3.412142).toFixed(3), y]; }
    return hit;
  },
};

/** Site-data fields of the cost-escalation source from the three index series ([[yyyy-mm, value], …]); the indexes are re-based to a stated year. */
export function costStats(cost = [], steel = [], cpi = []) {
  const data = {}; if (cost.length < 24) return data;
  const baseYear = +cost[cost.length - 1][0].slice(0, 4) - 9;
  const rebase = (rows) => { const b = rows.filter((r) => +r[0].slice(0, 4) === baseYear).map((r) => r[1]); return b.length >= 6 ? rows.map((r) => [r[0], +((100 * r[1]) / mean(b)).toFixed(2)]) : null; };
  const growth = (rows) => { const n = rows.length; return n > 60 && rows[n - 61][1] > 0 ? +(100 * ((rows[n - 1][1] / rows[n - 61][1]) ** 0.2 - 1)).toFixed(2) : null; }; // %/y over the last five years
  const c = rebase(cost), s = rebase(steel), p = rebase(cpi), last = (rows) => rows[rows.length - 1], series = (rows) => ({ t: rows.map((r) => r[0]), v: rows.map((r) => r[1]) });
  if (!c) return data;
  Object.assign(data, { costIndex: last(c)[1], costIndexDate: last(c)[0], costIndexBase: `${baseYear} = 100`, costIndexName: `US producer price index, oil and gas field machinery and equipment (BLS series ${BLS_COST})`, costIndexSeries: series(c) });
  const ce = growth(c); if (ce !== null) data.costEscalation = ce;
  if (s) { Object.assign(data, { steelIndex: last(s)[1], steelIndexDate: last(s)[0], steelIndexBase: `${baseYear} = 100`, steelIndexName: `US producer price index, steel mill products (BLS series ${BLS_STEEL})` }); const g = growth(s); if (g !== null) data.steelEscalation = g; }
  if (p) { Object.assign(data, { usCpi: last(p)[1], usCpiDate: last(p)[0], usCpiBase: `${baseYear} = 100`, usCpiName: `US consumer price index, all urban consumers, all items (BLS series ${BLS_CPI})`, usCpiSeries: series(p) }); const g = growth(p); if (g !== null) data.usInflation5y = g; }
  return data;
}
const MARKET_INFO = {
  gasEurope: ['gasPriceEurope', 'US$/MMBtu', 'Netherlands TTF day-ahead natural gas, monthly average'], gasAsia: ['gasPriceAsia', 'US$/MMBtu', 'LNG in Japan (Indonesian LNG import price), monthly average'],
  ironOre: ['steelPrice', 'US$/t', 'iron ore fines 62% Fe, spot, CFR China (a proxy for steel cost, not a steel price)'],
};
/** Site-data fields of the regional-prices source from { gasEurope, gasAsia, ironOre } series. */
export function marketStats(m = {}, source = SOURCE_NAMES.imfPcps) {
  const data = {};
  for (const [k, [key, unit, what]] of Object.entries(MARKET_INFO)) {
    const rows = m[k]; if (!rows?.length) continue;
    const last = rows[rows.length - 1]; if (!(last[1] > 0)) continue;
    Object.assign(data, { [key]: +last[1].toFixed(3), [key + 'Date']: String(last[0]), [key + 'Unit']: unit, [key + 'What']: what, [key + 'Source']: source });
    if (rows.length > 12) data[key + 'Series'] = { t: rows.map((r) => r[0]), v: rows.map((r) => +r[1].toFixed(3)) };
  }
  return data;
}
/** Benchmark gas region of a point: the Americas price off Henry Hub, Europe–Africa–Middle East off the European hub, Asia–Pacific off LNG delivered to Japan. */
export const gasRegion = (lon) => (lon < -30 ? 'americas' : lon < 60 ? 'europe' : 'asia');
const freshMonth = (p, months = 36) => { const t = Date.parse(String(p).length === 4 ? p + '-12-01' : p + '-01'); return Number.isFinite(t) && Date.now() - t < months * 30.5 * 24 * 3600e3; };
/** Interest-rate fields from { bondYield: [v, period], treasuryBillYield, policyRate }; observations older than three years are ignored. */
function rateFields(r = {}, source) {
  const data = {};
  for (const k of ['bondYield', 'treasuryBillYield', 'policyRate']) if (r[k] && Number.isFinite(r[k][0]) && freshMonth(r[k][1])) { data[k] = r[k][0]; data[k + 'Date'] = String(r[k][1]); data[k + 'Source'] = source; }
  return data;
}
const POWER_BASIS = { eurostat: 'industrial consumers using 2 000–19 999 MWh a year, excluding VAT and other recoverable taxes', eia: 'average revenue per kWh sold to industrial customers, mean of twelve months' };
const connectors = {
  async place(lat, lon) {
    const j = await getJSON(`https://api.bigdatacloud.net/data/reverse-geocode-client?latitude=${lat}&longitude=${lon}&localityLanguage=en`);
    let code = txt(j.countryCode, 2).toUpperCase(), country = txt(j.countryName), offshoreKm = null;
    if (!code) { // open sea: the nearest country (from the bundled outlines) supplies the national context — prices, tax, currency
      try { const p = await (await loadAtlas()).atlasPlace(lat, lon); if (p?.countryCode) { code = p.countryCode; country = p.country; offshoreKm = p.offshoreKm ?? null; } } catch { /* stays without a country */ }
    }
    const e = ENERGY[code]; // the tariff and the tax rate are no longer taken from hand-entered tables here: see the power and fiscal sources
    return { meta: { name: txt(j.city || j.locality || j.principalSubdivision || ''), country, countryCode: code, ...(offshoreKm ? { offshoreKm } : {}) },
      data: { ...(e ? { gridCarbon: e[1] } : {}), currency: CURRENCY[code] || 'USD' } };
  },
  async weather(lat, lon) {
    const j = await getJSON(`https://api.open-meteo.com/v1/forecast?latitude=${lat}&longitude=${lon}&current=temperature_2m,relative_humidity_2m,surface_pressure,wind_speed_10m,wind_direction_10m,shortwave_radiation&daily=shortwave_radiation_sum,temperature_2m_max,temperature_2m_min&hourly=soil_temperature_54cm&wind_speed_unit=ms&past_days=7&forecast_days=7&timezone=GMT`);
    const c = j.current || {}, d = j.daily || {};
    const ghi = (d.shortwave_radiation_sum || []).filter((x) => num(x) !== null).map((x) => x / 3.6); // MJ/m² -> kWh/m²
    const soil = (j.hourly?.soil_temperature_54cm || []).filter((x) => num(x) !== null);
    return { data: { airTemp: num(c.temperature_2m), humidity: num(c.relative_humidity_2m), pressure: num(c.surface_pressure), windSpeed: num(c.wind_speed_10m), windDir: num(c.wind_direction_10m), solar: num(c.shortwave_radiation),
      ghiDaily: ghi.length ? mean(ghi) : null, airTempMax: d.temperature_2m_max?.length ? Math.max(...d.temperature_2m_max.filter((x) => num(x) !== null)) : null, airTempMin: d.temperature_2m_min?.length ? Math.min(...d.temperature_2m_min.filter((x) => num(x) !== null)) : null, elevation: num(j.elevation), groundTemp: soil.length ? mean(soil) : null } };
  },
  async marine(lat, lon) {
    const j = await getJSON(`https://marine-api.open-meteo.com/v1/marine?latitude=${lat}&longitude=${lon}&current=wave_height,wave_direction,wave_period,sea_surface_temperature,ocean_current_velocity,ocean_current_direction&hourly=sea_level_height_msl,ocean_current_velocity,ocean_current_direction,sea_surface_temperature&past_days=3&forecast_days=5&timezone=GMT&cell_selection=sea`);
    const c = j.current || {}, hh = j.hourly || {}, n = (hh.time || []).length;
    const keep = (a) => (a || []).map((x) => num(x));
    const eta = keep(hh.sea_level_height_msl), sp = keep(hh.ocean_current_velocity).map((x) => (x === null ? null : x / 3.6)), dir = keep(hh.ocean_current_direction), sst = keep(hh.sea_surface_temperature).filter((x) => x !== null);
    const etaOk = eta.filter((x) => x !== null), spOk = sp.filter((x) => x !== null);
    // gaps in the source series are bridged by linear interpolation (never by zeros)
    const bridge = (a) => { const o = a.slice(); let last = -1; for (let i = 0; i < o.length; i++) { if (o[i] === null) continue; if (last < 0) for (let k = 0; k < i; k++) o[k] = o[i]; else for (let k = last + 1; k < i; k++) o[k] = o[last] + ((o[i] - o[last]) * (k - last)) / (i - last); last = i; } if (last >= 0) for (let k = last + 1; k < o.length; k++) o[k] = o[last]; return o.map((x) => x ?? 0); };
    const t0 = Date.parse(String(hh.time?.[0] || '') + 'Z'), nowHour = Number.isFinite(t0) ? clamp((Date.now() - t0) / 3600e3, 0, Math.max(0, n - 1)) : 72;
    const t = Array.from({ length: n }, (_, i) => i);
    const data = { waveHeight: num(c.wave_height), waveDir: num(c.wave_direction), wavePeriod: num(c.wave_period), sst: num(c.sea_surface_temperature) ?? (sst.length ? mean(sst) : null),
      currentSpeed: spOk.length ? mean(spOk) : c.ocean_current_velocity != null ? num(c.ocean_current_velocity) / 3.6 : null, currentMax: spOk.length ? Math.max(...spOk) : null, currentDir: num(c.ocean_current_direction),
      sstMin: sst.length ? Math.min(...sst) : null, sstMax: sst.length ? Math.max(...sst) : null };
    if (etaOk.length > 24) { data.tideRange = quantile(etaOk, 0.98) - quantile(etaOk, 0.02); data.seaLevelMean = mean(etaOk); data.tide = { t, eta: bridge(eta), nowHour }; }
    if (spOk.length > 24) data.currents = { t, speed: bridge(sp), dir: bridge(dir), nowHour, gaps: sp.length - spOk.length };
    if (data.sst === null && !etaOk.length && !spOk.length) throw new Error('No marine data at this point (inland site?)');
    return { data };
  },
  async bathy(lat, lon) {
    const half = 0.12, n = 57; // 57 × 57 samples ≈ the 15-arc-second native resolution of the mosaic over this window
    const lats = Array.from({ length: n }, (_, i) => clamp(lat - half + (2 * half * i) / (n - 1), -89.9, 89.9)), lons = Array.from({ length: n }, (_, i) => clamp(lon - half + (2 * half * i) / (n - 1), -179.9, 179.9));
    let elev;
    try { // primary: one multipoint sample request against the NOAA global DEM mosaic
      const points = lats.flatMap((la) => lons.map((lo) => [+lo.toFixed(5), +la.toFixed(5)]));
      const flat = new Array(n * n).fill(null), CH = 1000, jobs = []; // the service returns at most 1000 samples per request
      for (let o = 0; o < points.length; o += CH) jobs.push(getJSON('https://gis.ngdc.noaa.gov/arcgis/rest/services/DEM_mosaics/DEM_global_mosaic/ImageServer/getSamples', 12000,
        { geometry: JSON.stringify({ points: points.slice(o, o + CH), spatialReference: { wkid: 4326 } }), geometryType: 'esriGeometryMultipoint', returnFirstValueOnly: 'true', f: 'json' })
        .then((j) => { for (const sm of j.samples || []) { const v = parseFloat(sm.value), id = o + sm.locationId; if (Number.isInteger(sm.locationId) && id >= 0 && id < n * n && Number.isFinite(v) && Math.abs(v) < 12000) flat[id] = v; } }));
      await Promise.all(jobs);
      if (flat.filter((v) => v !== null).length < 0.6 * n * n) throw new Error('Incomplete relief grid');
      elev = lats.map((_, a) => lons.map((_, b) => flat[a * n + b] ?? 0));
    } catch { // fallback: SRTM30+ 1 km relief resampled onto the same grid
      const j = await getJSON(`https://pae-paha.pacioos.hawaii.edu/erddap/griddap/srtm30plus_v11_bathy.json?elev%5B(${lats[0].toFixed(4)}):1:(${lats[n - 1].toFixed(4)})%5D%5B(${lons[0].toFixed(4)}):1:(${lons[n - 1].toFixed(4)})%5D`, 20000);
      const rows = (j.table?.rows || []).filter((r) => num(r[0]) !== null && num(r[1]) !== null && num(r[2]) !== null);
      if (rows.length < 9) throw new Error('No relief data returned');
      elev = lats.map((la) => lons.map((lo) => { let best = rows[0], bd = Infinity; for (const r of rows) { const d = (r[0] - la) ** 2 + (r[1] - lo) ** 2; if (d < bd) { bd = d; best = r; } } return best[2]; }));
    }
    const mid = (n - 1) / 2, here = elev[mid][mid], flat = elev.flat();
    return { data: { bathy: { lat: lats, lon: lons, elev }, depth: here < 0 ? -here : 0, elevationRelief: here, maxDepthNearby: Math.max(0, -Math.min(...flat)), seaFraction: flat.filter((z) => z < 0).length / flat.length } };
  },
  async salinity(lat, lon) {
    const g = (v) => Math.round(v * 4) / 4, la = clamp(g(lat), -79.5, 79.5), lo = clamp(g(lon), -179.5, 179.25);
    const q = `%5B0:1:11%5D%5B0%5D%5B(${la - 0.25}):1:(${la + 0.25})%5D%5B(${lo - 0.25}):1:(${lo + 0.25})%5D`;
    const url = `https://erddap.emodnet-physics.eu/erddap/griddap/SDC_GLO_CLIM_TS_V2_1.json?Salinity${q},Temperature${q}`;
    const j = await getJSON(url, 20000); // slow server: its answer is cached on the device for 90 days, so this wait happens once per place
    const rows = (j.table?.rows || []).filter((r) => num(r[4]) !== null && r[4] > 0 && r[4] < 60);
    if (!rows.length) throw new Error('No ocean cell near this point');
    const byMonth = Array.from({ length: 12 }, () => ({ s: [], t: [] }));
    for (const r of rows) { const m = new Date(r[0]).getUTCMonth(); if (m >= 0) { byMonth[m].s.push(r[4]); if (num(r[5]) !== null) byMonth[m].t.push(r[5]); } }
    const sM = byMonth.map((b) => (b.s.length ? mean(b.s) : null)), tM = byMonth.map((b) => (b.t.length ? mean(b.t) : null)), sOk = sM.filter((x) => x !== null);
    return { data: { salinity: quantile(sOk, 0.5), salinityMin: Math.min(...sOk), salinityMax: Math.max(...sOk), salinityMonthly: sM.map((x) => x ?? 0), sstMonthly: tM.map((x) => x ?? 0) } };
  },
  async economy(lat, lon, site) {
    const code = site.countryCode;
    if (!code) throw new Error('Country unknown');
    const one = async (ind) => { const j = await getJSON(`https://api.worldbank.org/v2/country/${encodeURIComponent(code)}/indicator/${ind}?format=json&mrnev=1`); const r = j?.[1]?.[0]; return r && num(r.value) !== null ? { value: r.value, year: txt(String(r.date), 6), iso3: txt(r.countryiso3code, 3) } : null; };
    // key -> World Bank indicator
    const IND = { inflation: 'FP.CPI.TOTL.ZG', lendingRate: 'FR.INR.LEND', waterStress: 'ER.H2O.FWST.ZS', freshwaterPerCapita: 'ER.H2O.INTR.PC', renewableElectricity: 'EG.ELC.RNEW.ZS', tariffWB: 'IC.ELC.PRI.KH.DB1619', gdpPerCapita: 'NY.GDP.PCAP.CD', safeWaterAccess: 'SH.H2O.SMDW.ZS' };
    const keys = Object.keys(IND), res = await Promise.allSettled(keys.map((k) => one(IND[k]))), data = {};
    res.forEach((r, i) => { const v = r.status === 'fulfilled' ? r.value : null; if (!v) return; data[keys[i]] = v.value; data[keys[i] + 'Year'] = v.year; if (v.iso3) data.iso3 = v.iso3; });
    if (data.tariffWB != null) { data.electricityPriceWB = data.tariffWB / 100; delete data.tariffWB; data.electricityPriceWBYear = data.tariffWBYear; delete data.tariffWBYear; }
    if (!Object.keys(data).length) throw new Error('No indicators published for this country');
    return { data };
  },
  async energy(lat, lon, site) {
    const iso3 = site.data?.iso3, name = site.country;
    if (!iso3 && !name) throw new Error('Country unknown');
    // The chart addresses are requested directly: the short aliases answer with a redirect that browsers refuse cross-origin.
    const owid = (direct, alias) => getText('https://ourworldindata.org/grapher/' + direct).catch(() => getText('https://ourworldindata.org/grapher/' + alias));
    const pick = (csv) => { for (const line of csv.split('\n')) { const c = line.split(','); if ((iso3 && c[1] === iso3) || (!iso3 && c[0] === name)) { const v = parseFloat(c[3]); if (Number.isFinite(v)) return { value: v, year: txt(c[4] || c[2], 6) }; } } return null; };
    const [ci, rn] = await Promise.allSettled([owid('electricity-mix.csv?frequency=annual&metric=carbon_intensity&source=total&csvType=filtered&time=latest', 'carbon-intensity-electricity.csv?csvType=filtered&time=latest'), owid('electricity-mix.csv?frequency=annual&metric=share_of_generation&source=renewables&csvType=filtered&time=latest', 'share-electricity-renewables.csv?csvType=filtered&time=latest')]);
    const data = {}, a = ci.status === 'fulfilled' ? pick(ci.value) : null, b = rn.status === 'fulfilled' ? pick(rn.value) : null;
    if (a && a.value >= 0 && a.value < 1500) { data.gridCarbon = a.value / 1000; data.gridCarbonYear = a.year; data.gridCarbonLive = true; }
    if (b && b.value >= 0 && b.value <= 100) { data.renewableShare = b.value; data.renewableShareYear = b.year; }
    try { // national emissions-weighted carbon price (US$ per tonne CO₂); the series lists every country by year, so the last row of the country is the latest
      const csv = await getText('https://ourworldindata.org/grapher/emissions-weighted-carbon-price.csv?v=1&csvType=filtered', 20000, 900000); let hit = null;
      for (const line of csv.split('\n')) { const c = line.split(','); if ((iso3 && c[1] === iso3) || (!iso3 && c[0] === name)) { const v = parseFloat(c[3]); if (Number.isFinite(v)) hit = { value: v, year: txt(c[2], 6) }; } }
      if (hit && hit.value >= 0 && hit.value < 1000) { data.carbonPrice = hit.value; data.carbonPriceYear = hit.year; }
    } catch { /* optional */ }
    if (!Object.keys(data).length) throw new Error('Country not listed');
    return { data };
  },
  async climate(lat, lon) {
    const j = await getJSON(`https://power.larc.nasa.gov/api/temporal/climatology/point?parameters=ALLSKY_SFC_SW_DWN,WS10M,T2M&community=RE&longitude=${lon.toFixed(3)}&latitude=${lat.toFixed(3)}&format=JSON`, 25000);
    const p = j?.properties?.parameter || {}, M = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC'];
    const series = (o) => (o ? M.map((m) => (num(o[m]) !== null && o[m] > -900 ? o[m] : 0)) : null), ann = (o) => (o && num(o.ANN) !== null && o.ANN > -900 ? o.ANN : null);
    const data = { ghiAnnual: ann(p.ALLSKY_SFC_SW_DWN), ghiMonthly: series(p.ALLSKY_SFC_SW_DWN), windAnnual: ann(p.WS10M), windMonthly: series(p.WS10M), airTempAnnual: ann(p.T2M), airTempMonthly: series(p.T2M) };
    if (data.ghiAnnual === null) throw new Error('No climatology for this point');
    return { data };
  },
  async prices() {
    // Daily spot series as two-column CSV (date, price). Several independent copies are tried in turn.
    const parse = (csv) => csv.trim().split('\n').map((l) => l.split(',')).filter((c) => /^\d{4}-\d\d-\d\d$/.test(c[0]) && num(parseFloat(c[1])) !== null).map((c) => [c[0], parseFloat(c[1])]);
    const series = async (paths) => { let err; for (const u of paths) { try { const rows = parse(await getText(u, 20000, 4e6)); if (rows.length > 30) return rows; } catch (e) { err = e; } } throw err || new Error('No price series'); };
    const eia = async (path, id) => { const j = await getJSON(`https://api.eia.gov/v2/${path}/data/?api_key=DEMO_KEY&frequency=daily&data[0]=value&facets[series][]=${id}&sort[0][column]=period&sort[0][direction]=desc&length=260`); return (j?.response?.data || []).map((r) => [txt(r.period, 10), parseFloat(r.value)]).filter((r) => Number.isFinite(r[1])).reverse(); };
    const get = (dh, gh, path, id) => series([`https://raw.githubusercontent.com/datasets/${gh}`]).catch(() => eia(path, id));
    const [b, w, g] = await Promise.allSettled([get('oil-prices/r/brent-daily.csv', 'oil-prices/main/data/brent-daily.csv', 'petroleum/pri/spt', 'RBRTE'), get('oil-prices/r/wti-daily.csv', 'oil-prices/main/data/wti-daily.csv', 'petroleum/pri/spt', 'RWTC'), get('natural-gas/r/daily.csv', 'natural-gas/main/data/daily.csv', 'natural-gas/pri/fut', 'RNGWHHD')]);
    const data = {}, stat = (rows, key) => { if (!rows?.length) return; const last = rows[rows.length - 1], y = rows.slice(-252).map((r) => r[1]), m = rows.slice(-22).map((r) => r[1]); data[key] = last[1]; data[key + 'Date'] = last[0]; data[key + 'Mean30'] = +mean(m).toFixed(2); data[key + 'Min1y'] = Math.min(...y); data[key + 'Max1y'] = Math.max(...y); data[key + 'Series'] = { t: rows.slice(-252).map((r) => r[0]), v: y }; data[key + 'Volatility'] = +(Math.sqrt(252) * Math.sqrt(mean(y.slice(1).map((v, i) => Math.log(v / y[i]) ** 2)))).toFixed(3); };
    stat(b.status === 'fulfilled' ? b.value : null, 'oilPrice'); stat(w.status === 'fulfilled' ? w.value : null, 'oilPriceWTI'); stat(g.status === 'fulfilled' ? g.value : null, 'gasPrice');
    if (data.oilPrice == null && data.oilPriceWTI != null) { data.oilPrice = data.oilPriceWTI; data.oilPriceDate = data.oilPriceWTIDate; }
    if (data.oilPrice == null && data.gasPrice == null) throw new Error('No price series reachable');
    return { data };
  },
  async fx(lat, lon, site) {
    let cur = CURRENCY[site.countryCode] || site.data?.currency || 'USD';
    if (site.countryCode && !CURRENCY[site.countryCode]) { // not in the short table above: the bundled ISO 4217 list knows every country
      try { cur = (await (await loadAtlas()).atlasNational(site.countryCode))?.currency || cur; } catch { /* keep the default */ }
    }
    const j = await getJSON('https://open.er-api.com/v6/latest/USD');
    const rate = num(j?.rates?.[cur]);
    if (rate === null) throw new Error('Currency not listed');
    return { data: { currency: cur, fxPerUSD: rate, fxDate: txt(j.time_last_update_utc, 40) } };
  },
  async fiscal(lat, lon, site) {
    const code = site.countryCode, iso3 = site.data?.iso3 || ISO3[code];
    if (!code) throw new Error('Country unknown');
    let hit = null, source = SOURCE_NAMES.oecdCit;
    try { if (iso3) hit = (await feeds.corporateTax(iso3))[iso3] || null; } catch { /* the second table is tried below */ }
    if (!hit) { hit = await feeds.corporateTaxTF(code); source = SOURCE_NAMES.taxFoundation; }
    if (!hit) throw new Error('Country not listed');
    if (source === SOURCE_NAMES.taxFoundation && CORPORATE_TAX.rows[code]?.[1] > hit[1]) throw new Error('OECD table not reachable; the bundled copy is newer than the second source'); // the dated snapshot then stands in
    return { data: { corporateTaxRate: hit[0], corporateTaxYear: String(hit[1]), corporateTaxSource: source } };
  },
  async costs() {
    let x, source = SOURCE_NAMES.bls;
    try { x = await feeds.indexes(); }
    catch (e) { // daily request limit of the key-less service reached, or unreachable: the mirror, unless the bundled snapshot is newer
      const m = await feeds.blsMirror().catch(() => null), lo = `${new Date().getUTCFullYear() - 19}-01`, cut = (id) => (m[id] || []).filter((r) => r[0] >= lo), end = (rows) => rows?.[rows.length - 1]?.[0] || '';
      if (!m || end(m[BLS_COST]) < end(COSTS.cost)) throw e;
      x = { cost: cut(BLS_COST), steel: cut(BLS_STEEL), cpi: cut(BLS_CPI) }; source = SOURCE_NAMES.bls + ' (DBnomics mirror)';
    }
    const data = costStats(x.cost, x.steel, x.cpi);
    if (data.costIndex == null) throw new Error('No index series returned');
    return { data: { ...data, costIndexSource: source } };
  },
  async markets() {
    let data = {};
    try { data = marketStats(await feeds.commodities()); } catch { /* the annual series below still gives the European price */ }
    if (data.gasPriceEurope == null) { const g = await feeds.owidGasEurope(); if (g) Object.assign(data, { gasPriceEurope: g[0], gasPriceEuropeDate: String(g[1]), gasPriceEuropeUnit: 'US$/MMBtu', gasPriceEuropeWhat: 'Netherlands TTF natural gas, annual average', gasPriceEuropeSource: SOURCE_NAMES.owidGas }); }
    if (!Object.keys(data).length) throw new Error('No price series reachable');
    return { data };
  },
  async rates(lat, lon, site) {
    const code = site.countryCode, iso3 = site.data?.iso3 || ISO3[code];
    if (!code) throw new Error('Country unknown');
    let data = {};
    try { if (iso3) data = rateFields((await feeds.imfRates(iso3))[iso3], SOURCE_NAMES.imfMfs); } catch { /* fallbacks below */ }
    const [b, p] = await Promise.allSettled([data.bondYield == null && iso3 ? feeds.oecdBond(iso3) : null, data.policyRate == null ? feeds.bisPolicy(code) : null]);
    if (b.status === 'fulfilled' && b.value?.[iso3]) data = { ...rateFields({ bondYield: b.value[iso3] }, SOURCE_NAMES.oecdRates), ...data };
    if (p.status === 'fulfilled' && p.value) data = { ...rateFields({ policyRate: p.value }, SOURCE_NAMES.bis), ...data };
    if (!Object.keys(data).length) throw new Error('No interest rates published for this country');
    return { data };
  },
  async power(lat, lon, site) {
    const code = site.countryCode;
    if (!code) throw new Error('Country unknown');
    if (code === 'US') { const r = await feeds.eiaPower(); return { data: { electricityPrice: r[0], electricityPriceDate: r[1], electricityPriceSource: SOURCE_NAMES.eia, electricityPriceBasis: POWER_BASIS.eia } }; }
    if (!EUROSTAT_GEO[code]) throw new Error('No open industrial tariff series for this country');
    const r = (await feeds.eurostatPower(code))[code];
    if (!r) throw new Error('No price published for this country');
    let eur = null; // euro per US dollar: today's rate, or the one bundled with the snapshot
    try { eur = num((await getJSON('https://open.er-api.com/v6/latest/USD'))?.rates?.EUR); } catch { /* bundled rate */ }
    if (!(eur > 0.3 && eur < 3)) eur = ELECTRICITY.eurPerUSD;
    return { data: { electricityPrice: +(r[0] / eur).toFixed(4), electricityPriceDate: r[1], electricityPriceEUR: r[0], electricityPriceSource: SOURCE_NAMES.eurostat, electricityPriceBasis: POWER_BASIS.eurostat } };
  },
};
/** One live source, straight from the network (no device cache, no atlas): { data, meta? }. Used by the tests and tools. */
export const fetchSource = (id, lat, lon, site = {}) => connectors[id](lat, lon, { data: {}, ...site });

// ---- response cache ---------------------------------------------------------------------------------
// Slow-changing answers (relief, climatologies, national indicators) are kept on the device so that a
// repeat fetch of the same place is instant; fast-changing ones (weather, sea state) are kept briefly.
const CACHE_KEY = 'hydraslug.live.v1', H = 3600e3;
const TTL = { place: 30 * 24 * H, weather: 0.5 * H, marine: 0.5 * H, bathy: 90 * 24 * H, salinity: 90 * 24 * H, economy: 7 * 24 * H, fx: 12 * H, energy: 7 * 24 * H, climate: 90 * 24 * H, prices: 6 * H, fiscal: 30 * 24 * H, costs: 7 * 24 * H, markets: 24 * H, rates: 7 * 24 * H, power: 30 * 24 * H };
let memCache = null;
function cacheAll() {
  if (memCache) return memCache;
  try { memCache = JSON.parse(localStorage.getItem(CACHE_KEY) || '{}') || {}; } catch { memCache = {}; }
  return memCache;
}
function cacheGet(key, ttl) { const e = cacheAll()[key]; return e && Date.now() - e.t < ttl ? e : null; }
function cachePut(key, v) {
  const all = cacheAll(); all[key] = { t: Date.now(), v };
  const keys = Object.keys(all);
  if (keys.length > 160) keys.sort((a, b) => all[a].t - all[b].t).slice(0, keys.length - 120).forEach((k) => delete all[k]);
  try { localStorage.setItem(CACHE_KEY, JSON.stringify(all)); }
  catch { try { for (const k of Object.keys(all)) if (k.startsWith('bathy')) delete all[k]; localStorage.setItem(CACHE_KEY, JSON.stringify(all)); } catch { /* storage unavailable: memory cache only */ } }
}
/** Cache key of a connector: grid-cell for site data, country for national data. */
function keyOf(id, lat, lon, site) {
  if (id === 'economy' || id === 'energy' || id === 'fiscal' || id === 'rates' || id === 'power') return site.countryCode ? `${id}:${site.countryCode}` : null;
  if (id === 'costs' || id === 'markets') return `${id}:world`;
  if (id === 'place') return `place2:${Math.round(lat * 100)}:${Math.round(lon * 100)}`; // answers stored by earlier builds carried hand-entered tariff and tax figures
  if (id === 'fx') return `fx:${site.data?.currency || site.countryCode || 'USD'}`;
  if (id === 'prices') return 'prices:world';
  const r = id === 'salinity' ? 4 : id === 'climate' ? 2 : id === 'bathy' ? 250 : id === 'place' ? 100 : 20; // cells per degree
  return `${id}:${Math.round(lat * r)}:${Math.round(lon * r)}`;
}

// ---- built-in world atlas -----------------------------------------------------------------------------
// When a live source cannot be reached (no connection, a blocked host, a service that is down) the fields
// it would have supplied are answered from data bundled with the app: js/data/atlas_*.js, built from open
// global data sets by tools/atlas_build.py. Those modules are large, so they are imported only when needed.
let atlasLib = null;
/** The look-up module of the built-in atlas (js/data/atlas_lookup.js), loaded on first use. */
export const loadAtlas = () => (atlasLib ??= import('../data/atlas_lookup.js').catch((e) => { atlasLib = null; throw e; }));

/** Fields by which a source counts as "answered by the atlas", and the atlas data set that stands in for it. */
const ATLAS_SOURCE = {
  place: ['nations', ['country']], economy: ['nations', ['inflation', 'lendingRate', 'gdpPerCapita', 'waterStress', 'renewableElectricity', 'electricityPriceWB', 'freshwaterPerCapita', 'safeWaterAccess']],
  fx: ['fx', ['fxPerUSD']], energy: ['nations', ['gridCarbon', 'renewableShare']], weather: ['climate', ['airTemp', 'windSpeed', 'ghiDaily', 'elevation']],
  marine: ['coast', ['tideRange', 'waveHeight', 'currentSpeed', 'tide', 'sst']], bathy: ['relief', ['bathy', 'depth']], salinity: ['ocean', ['salinity', 'salinityMonthly']], climate: ['climate', ['ghiAnnual', 'windAnnual', 'airTempAnnual']], prices: ['prices', ['oilPrice', 'gasPrice']],
  fiscal: ['fiscal', ['corporateTaxRate']], costs: ['costs', ['costIndex']], markets: ['markets', ['gasPriceEurope', 'gasPriceAsia', 'steelPrice']], rates: ['rates', ['bondYield', 'policyRate', 'treasuryBillYield']], power: ['power', ['electricityPrice']],
};
/** Plain-language name of every site-data field the atlas can supply (used for labels and the notice on the site page). */
export const ATLAS_LABELS = {
  country: 'country', currency: 'currency', electricityPrice: 'electricity price', inflation: 'inflation', lendingRate: 'lending rate', gdpPerCapita: 'GDP per capita', waterStress: 'water stress', renewableElectricity: 'renewable electricity (World Bank)',
  electricityPriceWB: 'business electricity tariff', freshwaterPerCapita: 'freshwater per capita', safeWaterAccess: 'safe water access', fxPerUSD: 'exchange rate', gridCarbon: 'grid carbon', renewableShare: 'renewable share',
  salinity: 'salinity', salinityMonthly: 'monthly salinity', sstMonthly: 'monthly sea temperature', sst: 'sea temperature', bathy: 'seabed and terrain grid', depth: 'water depth', maxDepthNearby: 'deepest point nearby', elevation: 'land elevation',
  tideRange: 'tidal range', tide: 'tide series', currentSpeed: 'mean current', currentMax: 'peak current', currentDir: 'current direction', waveHeight: 'wave height', wavePeriod: 'wave period', waveDir: 'wave direction',
  oilPrice: 'oil price', oilPriceWTI: 'WTI oil price', gasPrice: 'gas price', seabedTemp: 'seabed temperature', tempProfile: 'sea temperature profile',
  ghiAnnual: 'solar resource', ghiDaily: 'solar resource', windAnnual: 'long-term wind', airTempAnnual: 'long-term air temperature', airTemp: 'air temperature', windSpeed: 'wind speed',
  corporateTaxRate: 'corporate tax rate', taxRate: 'tax rate', costIndex: 'equipment cost index', costEscalation: 'cost escalation', steelIndex: 'steel price index', usCpi: 'US consumer price index', usInflation5y: 'US inflation (5-year)',
  gasPriceEurope: 'European gas price', gasPriceAsia: 'Asian LNG price', gasPriceRegional: 'regional gas price', steelPrice: 'iron ore price', bondYield: 'government bond yield', treasuryBillYield: 'Treasury bill yield', policyRate: 'central-bank policy rate',
};
const month = () => new Date().getUTCMonth();

/**
 * Fill every site-data field that is still missing from the built-in atlas. Live values are never overwritten.
 *   site    { lat, lon, data, country, countryCode, name }  — changed in place
 *   status  per-source status of the live fetch; sources answered by the atlas get { atlas: true, message: 'built-in atlas (yyyy-mm)' }
 *   ctx     { mark: Map(field -> atlas data set), notes: {}, done: Set(source ids that have settled), vintage: {} }
 * A field is filled only once every live source that could still deliver it has settled (ctx.done; default: all).
 * Afterwards site.data.atlasFields lists the fields taken from the atlas, site.data.atlasVintage the dates of the
 * bundled data sets in use and site.data.atlasNotes any caveat per field.
 */
export async function atlasFill(site, status = {}, ctx = {}) {
  const d = (site.data ??= {}), { lat, lon } = site;
  const mark = (ctx.mark ??= new Map()), notes = (ctx.notes ??= {}), done = (ctx.done ??= new Set(SOURCES.map((x) => x.id))), vintage = (ctx.vintage ??= {});
  const ready = (...ids) => ids.every((id) => done.has(id)), missing = (...keys) => keys.some((k) => d[k] == null);
  const put = (set, k, v, note) => { if (v == null || (d[k] != null && !mark.has(k))) return; d[k] = v; mark.set(k, set); if (note) notes[k] = note; else delete notes[k]; };
  let A = null;
  try { A = await loadAtlas(); } catch { /* the atlas modules could not be loaded: the small table shipped with the core is used below */ }
  const part = async (fn) => { try { await fn(); } catch { /* this data set is unavailable; the others still apply */ } };
  if (A) {
    await part(async () => { // country and national figures
      if (!ready('place')) return;
      if (!site.countryCode) { const p = await A.atlasPlace(lat, lon); if (p) { site.country = p.country; site.countryCode = p.countryCode; mark.set('country', 'nations'); if (p.offshoreKm > 0) notes.country = `nearest country, ${p.offshoreKm} km from its outline`; } }
      if (!site.countryCode) return;
      const eco = ['inflation', 'lendingRate', 'gdpPerCapita', 'waterStress', 'renewableElectricity', 'freshwaterPerCapita', 'safeWaterAccess', 'electricityPriceWB'];
      const wantEco = ready('economy') && missing(...eco), wantEn = ready('energy') && (!d.gridCarbonLive || missing('renewableShare')), wantFx = ready('fx') && missing('fxPerUSD');
      if (!wantEco && !wantEn && !wantFx && !missing('currency')) return;
      const n = await A.atlasNational(site.countryCode);
      if (!n) return;
      if (!site.country) site.country = n.countryName;
      if (n.currency && (d.currency == null || (d.currency === 'USD' && !CURRENCY[site.countryCode] && !status.fx?.ok))) { d.currency = n.currency; mark.set('currency', 'nations'); }
      if (wantEco) { for (const k of eco) { put('nations', k, n[k]); if (mark.has(k)) d[k + 'Year'] = n[k + 'Year']; } d.iso3 ??= n.iso3; }
      if (wantEn) {
        if (!d.gridCarbonLive && n.gridCarbon != null) { d.gridCarbon = n.gridCarbon; d.gridCarbonYear = n.gridCarbonYear; mark.set('gridCarbon', 'nations'); }
        put('nations', 'renewableShare', n.renewableShare); if (mark.has('renewableShare')) d.renewableShareYear = n.renewableShareYear;
      }
      if (wantFx) { const r = d.currency === n.currency ? n : await A.atlasRate(d.currency); if (r?.fxPerUSD > 0) { put('fx', 'fxPerUSD', r.fxPerUSD); if (mark.has('fxPerUSD')) d.fxDate = r.fxDate; } }
      const e = ENERGY[site.countryCode];
      if (e && ready('energy') && d.gridCarbon == null) put('nations', 'gridCarbon', e[1], 'estimate: hand-entered planning default');
    });
    await part(async () => { // sea-surface salinity and temperature
      const wantS = ready('salinity') && missing('salinity', 'salinityMonthly', 'sstMonthly'), wantT = ready('salinity', 'marine') && missing('sst');
      if (!wantS && !wantT) return;
      const o = await A.atlasOcean(lat, lon);
      if (!o) return;
      const note = o.oceanCellKm > 0 ? `nearest sea cell, ${o.oceanCellKm} km away` : '';
      if (wantS) { put('ocean', 'salinity', o.salinity, note); if (mark.has('salinity')) { d.salinityMin = o.salinityMin; d.salinityMax = o.salinityMax; } put('ocean', 'salinityMonthly', o.salinityMonthly); put('ocean', 'sstMonthly', o.sstMonthly); }
      if (wantT && (d.sstMonthly || o.sstMonthly)) { const own = mark.has('sstMonthly') || d.sstMonthly == null, m = own ? o.sstMonthly : d.sstMonthly; put(own ? 'ocean' : 'derived', 'sst', m[month()], 'climatological mean for this month'); if (mark.has('sst')) { d.sstMin = Math.min(...m); d.sstMax = Math.max(...m); } }
    });
    await part(async () => { // seabed and terrain
      const wantB = ready('bathy') && missing('bathy', 'depth'), wantE = ready('bathy', 'weather') && missing('elevation');
      if (!wantB && !wantE) return;
      const r = await A.atlasRelief(lat, lon);
      if (wantB) {
        put('relief', 'bathy', r.bathy, 'coarse grid (about 5 km at the coast)'); put('relief', 'depth', r.depth, r.depthEstimated || '');
        if (mark.has('depth')) { d.elevationRelief = r.elevationRelief; d.seaFraction = r.seaFraction; d.depthEstimated = r.depthEstimated || 'built-in atlas relief'; }
        put('relief', 'maxDepthNearby', r.maxDepthNearby, 'within about 50 km');
      }
      if (wantE) put('relief', 'elevation', r.elevation);
    });
    await part(async () => { // tides, currents, waves
      if (!ready('marine') || !missing('tideRange', 'tide', 'currentSpeed', 'waveHeight', 'wavePeriod')) return;
      const c = await A.atlasCoast(lat, lon);
      if (!c) return;
      const near = c.coastPointKm != null ? `, model point ${c.coastPointKm} km away` : '';
      put('coast', 'tideRange', c.tideRange, `typical range ${c.seaPeriod}${near}`);
      if (mark.has('tideRange')) { d.tideSpring = c.tideSpring; d.tideNeap = c.tideNeap; }
      if (d.tide == null || mark.has('tide')) { put('coast', 'tide', c.tide, 'synthetic: harmonic prediction from five tidal constituents'); if (mark.has('tide')) d.seaLevelMean = 0; }
      put('coast', 'currentSpeed', c.currentSpeed, `mean ${c.seaPeriod}`); put('coast', 'currentMax', c.currentMax, `99th percentile ${c.seaPeriod}`); put('coast', 'currentDir', c.currentDir, 'predominant direction');
      const wn = c.wavePointKm != null ? `, model point ${c.wavePointKm} km away` : '';
      put('coast', 'waveHeight', c.waveHeight, `annual mean ${c.wavePeriodOfRecord.slice(0, 4)}${wn}`);
      if (mark.has('waveHeight')) { d.waveHeightP95 = c.waveHeightP95; d.waveHeightMax = c.waveHeightMax; d.waveHeightMonthly = c.waveHeightMonthly; d.wavePeriodStorm = c.wavePeriodStorm; }
      put('coast', 'wavePeriod', c.wavePeriod, 'annual mean'); put('coast', 'waveDir', c.waveDir, 'predominant direction');
    });
    await part(async () => { // long-term solar, wind, air temperature
      const wantC = ready('climate') && missing('ghiAnnual', 'windAnnual', 'airTempAnnual'), wantW = ready('climate', 'weather') && missing('airTemp', 'windSpeed', 'ghiDaily');
      if (!wantC && !wantW) return;
      const c = await A.atlasClimate(lat, lon);
      if (!c) return;
      if (wantC) for (const k of ['ghiAnnual', 'ghiMonthly', 'windAnnual', 'windMonthly', 'airTempAnnual', 'airTempMonthly']) put('climate', k, c[k], k.endsWith('Annual') ? `long-term mean ${c.climatePeriod}` : '');
      if (wantW) {
        const own = (k) => (mark.has(k) || d[k] == null ? 'climate' : 'derived'), at = d.airTempMonthly || c.airTempMonthly, wm = d.windMonthly || c.windMonthly;
        if (at) put(own('airTempMonthly'), 'airTemp', at[month()], 'long-term mean for this month'); if (wm) put(own('windMonthly'), 'windSpeed', wm[month()], 'long-term mean for this month');
        put(own('ghiAnnual'), 'ghiDaily', d.ghiAnnual ?? c.ghiAnnual, 'long-term mean');
      }
    });
  }
  if (ready('prices') && missing('oilPrice', 'gasPrice')) { // commodity prices known when this build was made
    const note = `price on ${PRICES.date}, bundled with this build`;
    put('prices', 'oilPrice', PRICES.brent, note); if (mark.has('oilPrice')) d.oilPriceDate = PRICES.date;
    put('prices', 'oilPriceWTI', PRICES.wti, note); put('prices', 'gasPrice', PRICES.henryHub, note); if (mark.has('gasPrice')) d.gasPriceDate = PRICES.date;
    vintage.prices = PRICES.date;
  }
  // ---- dated snapshot bundled with the build (js/data/prices.js, rewritten by tools/prices.mjs from the live series) ----
  const code = site.countryCode, tag = (k) => mark.has(k);
  if (ready('costs') && missing('costIndex')) {
    const c = costStats(COSTS.cost, COSTS.steel, COSTS.cpi), note = `index of ${c.costIndexDate}, snapshot of ${COSTS.retrieved} bundled with this build`;
    put('costs', 'costIndex', c.costIndex, note);
    if (tag('costIndex')) { for (const k of ['costIndexDate', 'costIndexBase', 'costIndexName', 'costIndexSeries']) d[k] = c[k]; d.costIndexSource = SOURCE_NAMES.bls; put('costs', 'costEscalation', c.costEscalation, `five years to ${c.costIndexDate}`); vintage.costs = COSTS.retrieved; }
    put('costs', 'steelIndex', c.steelIndex, note); if (tag('steelIndex')) for (const k of ['steelIndexDate', 'steelIndexBase', 'steelIndexName', 'steelEscalation']) d[k] = c[k];
    put('costs', 'usCpi', c.usCpi, note); if (tag('usCpi')) { for (const k of ['usCpiDate', 'usCpiBase', 'usCpiName', 'usCpiSeries']) d[k] = c[k]; put('costs', 'usInflation5y', c.usInflation5y, `five years to ${c.usCpiDate}`); }
  }
  if (ready('markets') && missing('gasPriceEurope', 'gasPriceAsia', 'steelPrice')) {
    const m = marketStats(Object.fromEntries(Object.entries(MARKETS.rows).map(([k, v]) => [k, [[v[1], v[0]]]])));
    for (const key of ['gasPriceEurope', 'gasPriceAsia', 'steelPrice']) { put('markets', key, m[key], `monthly average of ${m[key + 'Date']}, bundled with this build`); if (tag(key)) { for (const x of ['Date', 'Unit', 'What', 'Source']) d[key + x] = m[key + x]; delete d[key + 'Series']; vintage.markets = MARKETS.retrieved; } }
  }
  if (code && ready('place', 'fiscal') && missing('corporateTaxRate')) {
    const r = CORPORATE_TAX.rows[code];
    if (r) { put('fiscal', 'corporateTaxRate', r[0], `statutory rate of ${r[1]}, snapshot of ${CORPORATE_TAX.retrieved} bundled with this build`); if (tag('corporateTaxRate')) { d.corporateTaxYear = String(r[1]); d.corporateTaxSource = SOURCE_NAMES.oecdCit; vintage.fiscal = CORPORATE_TAX.retrieved; } }
  }
  if (code && ready('place', 'rates') && missing('bondYield', 'policyRate', 'treasuryBillYield')) {
    const r = rateFields(RATES.rows[code], SOURCE_NAMES.imfMfs);
    for (const k of ['bondYield', 'treasuryBillYield', 'policyRate']) { put('rates', k, r[k], `value of ${r[k + 'Date']}, bundled with this build`); if (tag(k)) { d[k + 'Date'] = r[k + 'Date']; d[k + 'Source'] = r[k + 'Source']; vintage.rates = RATES.retrieved; } }
  }
  if (code && ready('place', 'power') && (d.electricityPrice == null || mark.has('electricityPrice'))) { // sourced snapshot, then the hand-entered estimate, then the World Bank business tariff
    const r = ELECTRICITY.rows[code], e = ENERGY[code], src = (date, source, basis) => { if (date) d.electricityPriceDate = date; else delete d.electricityPriceDate; d.electricityPriceSource = source; d.electricityPriceBasis = basis; delete d.electricityPriceEUR; };
    if (r) { put('power', 'electricityPrice', r[0], `industrial price of ${r[1]}, snapshot of ${ELECTRICITY.retrieved} bundled with this build`); if (mark.get('electricityPrice') === 'power') { src(r[1], code === 'US' ? SOURCE_NAMES.eia : SOURCE_NAMES.eurostat, code === 'US' ? POWER_BASIS.eia : POWER_BASIS.eurostat); vintage.power = ELECTRICITY.retrieved; } }
    else if (e && e[0] != null) { put('estimate', 'electricityPrice', e[0], 'ESTIMATE: hand-entered planning default — no open industrial tariff series exists for this country'); if (mark.get('electricityPrice') === 'estimate') src('', 'estimate (hand-entered planning default)', 'indicative industrial tariff'); }
    else if (ready('economy') && d.electricityPriceWB != null) { put(mark.has('electricityPriceWB') ? 'nations' : 'derived', 'electricityPrice', d.electricityPriceWB, `World Bank business tariff of ${d.electricityPriceWBYear} (small commercial connection, not an industrial rate)`); src(String(d.electricityPriceWBYear || ''), 'World Bank, Doing Business: price of electricity', 'small commercial connection, not an industrial rate'); }
  }
  // last resort (also when the atlas modules are unavailable): the small table shipped with the core
  if (ready('salinity') && d.salinity == null) {
    const at = A ? null : atlasSite(lat, lon), reg = regionalSalinity(lat, lon), enclosed = reg !== 35.5 && reg !== 34.3;
    put('core', 'salinity', at && !enclosed ? at.salinity : reg, at && !enclosed ? 'coarse 4° table' : 'regional estimate');
    if (at) { put('core', 'salinityMonthly', at.salinityMonthly); put('core', 'sstMonthly', at.sstMonthly); if (ready('marine')) put('core', 'sst', at.sst, 'coarse 4° table'); }
    else if (d.sstMonthly && d.salinityMonthly == null) put('core', 'salinityMonthly', new Array(12).fill(d.salinity), 'regional estimate, no seasonal cycle');
  }
  if (ready('climate', 'weather') && d.ghiAnnual == null && d.ghiDaily == null) put('core', 'ghiDaily', solarEstimate(lat), 'estimate from latitude');
  // dates of the bundled data, and the per-source status
  const sets = new Set(mark.values());
  if (A) for (const s of sets) if (s !== 'core' && s !== 'derived' && s !== 'estimate' && vintage[s] === undefined) Object.assign(vintage, await A.atlasVintage([s === 'fx' ? 'nations' : s]).catch(() => ({})));
  for (const k of [...mark.keys()]) if (mark.get(k) === 'derived') mark.delete(k); // computed from live climatology: not atlas data
  for (const [id, [set, keys]] of Object.entries(ATLAS_SOURCE)) {
    const st = status[id];
    if (!st || st.ok || st.atlas || !keys.some((k) => mark.has(k) && mark.get(k) !== 'estimate')) continue;
    status[id] = { ok: false, atlas: true, message: `built-in atlas (${vintage[set] || vintage.nations || 'bundled'})`, reason: st.message, at: st.at };
    ctx.onStatus?.(id, 'atlas', status[id].message);
  }
  atlasStamp(site, ctx);
  return site;
}
/** Write the atlas book-keeping (atlasFields, atlasNotes, atlasVintage and the older flags) into site.data. */
function atlasStamp(site, ctx) {
  const d = site.data, mark = ctx.mark, sets = new Set(mark.values());
  d.atlasFields = [...mark.keys()].sort();
  d.atlasNotes = Object.fromEntries(Object.entries(ctx.notes).filter(([k]) => mark.has(k)));
  const v = Object.fromEntries(Object.entries(ctx.vintage || {}).filter(([k]) => sets.has(k) || (k === 'seaPeriod' || k === 'wavePeriod' ? sets.has('coast') : k === 'climatePeriod' ? sets.has('climate') : false)));
  if (Object.keys(v).length) d.atlasVintage = v; else delete d.atlasVintage;
  d.atlas = mark.size > 0; d.salinityEstimated = mark.has('salinity'); d.sstEstimated = mark.has('sst'); d.solarEstimated = mark.get('ghiDaily') === 'core';
  if (!mark.has('depth')) delete d.depthEstimated;
}
/**
 * Merge a fresh fetch into the stored data of the same place: fresh live values replace anything, but a fresh
 * atlas value never replaces a stored live one. Returns the merged data object (with a consistent atlasFields).
 */
export function mergeSiteData(old = {}, fresh = {}) {
  const oldAtlas = new Set(old.atlasFields || []), newAtlas = new Set(fresh.atlasFields || []), out = { ...old }, fields = new Set();
  for (const k of oldAtlas) if (old[k] != null) fields.add(k);
  for (const [k, v] of Object.entries(fresh)) {
    if (v === null || v === undefined || k.startsWith('atlas')) continue;
    if (newAtlas.has(k) && old[k] != null && !oldAtlas.has(k)) continue; // stored live value stays
    out[k] = v; if (newAtlas.has(k)) fields.add(k); else fields.delete(k);
  }
  out.atlasFields = [...fields].sort();
  out.atlasNotes = Object.fromEntries(Object.entries({ ...(old.atlasNotes || {}), ...(fresh.atlasNotes || {}) }).filter(([k]) => fields.has(k)));
  out.atlasVintage = { ...(old.atlasVintage || {}), ...(fresh.atlasVintage || {}) }; if (!fields.size) delete out.atlasVintage;
  out.atlas = fields.size > 0; out.salinityEstimated = fields.has('salinity'); out.sstEstimated = fields.has('sst'); out.solarEstimated = fields.has('ghiDaily') && fresh.solarEstimated === true;
  if (!fields.has('depth')) delete out.depthEstimated;
  return out;
}

/**
 * Pull everything for one location, as fast as the sources allow:
 *  - every connector starts at once (national data start the moment the country is known);
 *  - answers already on the device are used immediately;
 *  - onData(site) is called after each answer so the page fills in progressively instead of waiting for the slowest source;
 *  - whatever a source fails to deliver is filled from the built-in atlas (atlasFill), field by field, and listed in
 *    data.atlasFields; a live answer always replaces an atlas value.
 * opt: { fresh: ignore the device cache, atlasOnly: do not use the network at all }.
 * Returns { meta, data, status: { id: { ok, message, at, cached, atlas } } }.
 */
export async function fetchSite(lat, lon, onStatus = () => {}, onData = () => {}, opt = {}) {
  const { fresh = false } = opt;
  lat = clamp(+lat, -90, 90); lon = ((((+lon + 180) % 360) + 360) % 360) - 180;
  const site = { lat, lon, data: {}, name: '', country: '', countryCode: '' }, status = {};
  const ctx = { mark: new Map(), notes: {}, done: new Set(), vintage: {}, onStatus };
  const snapshot = () => { atlasStamp(site, ctx); return { ...site, data: { ...site.data }, status: { ...status }, fetchedAt: new Date().toISOString() }; };
  const absorb = (r) => { // live (or stored live) answer: it replaces any atlas value of the same field
    for (const [k, v] of Object.entries(r.data || {})) if (v !== null && v !== undefined) { site.data[k] = v; ctx.mark.delete(k); }
    if (r.meta) { Object.assign(site, r.meta); if (r.meta.countryCode) ctx.mark.delete('country'); }
  };
  const run = async (id) => {
    const key = keyOf(id, lat, lon, site), hit = !fresh && key ? cacheGet(key, TTL[id]) : null;
    if (hit) { absorb(hit.v); ctx.done.add(id); status[id] = { ok: true, message: 'Live', at: new Date(hit.t).toISOString(), cached: true }; onStatus(id, 'ok', 'Live'); finish(); onData(snapshot()); return; }
    onStatus(id, 'loading');
    try {
      const r = await connectors[id](lat, lon, site);
      absorb(r);
      if (key) cachePut(key, { data: r.data, meta: r.meta });
      status[id] = { ok: true, message: 'Live', at: new Date().toISOString() };
    } catch (e) {
      const stale = key ? cacheAll()[key] : null; // an older stored answer beats no answer
      if (stale) { absorb(stale.v); status[id] = { ok: true, message: 'Stored copy', at: new Date(stale.t).toISOString(), cached: true }; }
      else status[id] = { ok: false, message: e.name === 'AbortError' ? 'Timed out' : txt(e.message || 'Unavailable', 90), at: new Date().toISOString() };
    }
    ctx.done.add(id);
    onStatus(id, status[id].ok ? 'ok' : 'fail', status[id].message);
    if (!status[id].ok) await atlasFill(site, status, ctx); // stand in for this source straight away (the country, in particular, unlocks the national sources)
    finish(); onData(snapshot());
  };
  const finish = () => { // values derived from other fields
    const d = site.data, mk = ctx.mark, derive = (k, from, v) => { d[k] = v; if (mk.has(from)) mk.set(k, mk.get(from)); else mk.delete(k); };
    if (d.sst == null && d.sstMonthly && d.sstMonthly[month()]) derive('sst', 'sstMonthly', d.sstMonthly[month()]);
    if (d.ghiAnnual != null && (d.ghiDaily == null || mk.has('ghiDaily') || !mk.has('ghiAnnual'))) derive('ghiDaily', 'ghiAnnual', d.ghiAnnual); // the long-term mean is the better design basis than this week's weather
    // Headline tax rate: the sourced petroleum fiscal table bundled with the app, else the statutory corporate income tax rate.
    const f = fiscalOf(site.countryCode);
    if (f) {
      for (const [k, v] of [['fiscalRegime', f.regime], ['royaltyRate', f.royalty], ['petroleumTaxRate', f.petroleumTax], ['upstreamCorporateTax', f.corporateTax], ['marginalTake', f.marginalTake], ['fiscalNote', f.note], ['fiscalYear', String(f.year)], ['fiscalSource', f.source.citation + (f.source2 ? '; ' + f.source2.citation : '')], ['fiscalUrl', f.source.url]]) if (v != null && v !== '') d[k] = v;
    }
    if (f && f.headline != null) { d.taxRate = f.headline; d.taxRateYear = String(f.year); d.taxRateSource = f.source.citation; d.taxRateBasis = f.basis; mk.delete('taxRate'); }
    else if (d.corporateTaxRate != null) { d.taxRate = d.corporateTaxRate; mk.delete('taxRate'); d.taxRateYear = d.corporateTaxYear; d.taxRateSource = d.corporateTaxSource; d.taxRateBasis = 'statutory corporate income tax (no petroleum-specific rate on file for this country)' + (mk.has('corporateTaxRate') ? ', from the snapshot bundled with this build' : ''); }
    else if (d.taxRate != null && !d.taxRateSource) delete d.taxRate; // a figure stored by an earlier build, without a source
    // Regional gas benchmark nearest the site.
    const reg = gasRegion(site.lon), gk = reg === 'americas' ? 'gasPrice' : reg === 'europe' ? 'gasPriceEurope' : 'gasPriceAsia';
    if (d[gk] != null) { derive('gasPriceRegional', gk, d[gk]); d.gasPriceRegionalDate = d[gk + 'Date']; d.gasPriceRegionalName = reg === 'americas' ? 'US Henry Hub' : reg === 'europe' ? 'Europe (Netherlands TTF)' : 'Asia (LNG delivered to Japan)'; }
    if (d.bathy && d.bathy.elev?.length > 4) { // seabed / ground slope at the site from the relief grid (central differences over about 1 km)
      const b = d.bathy, n = b.lat.length, m = b.lon.length, i = Math.floor(n / 2), j = Math.floor(m / 2), k = Math.max(1, Math.round(n / 28));
      const dy = (b.lat[i + k] - b.lat[i - k]) * 110540, dx = (b.lon[j + k] - b.lon[j - k]) * 111320 * Math.cos((site.lat * Math.PI) / 180);
      if (dy > 0 && dx > 0) derive('seabedSlope', 'bathy', +((Math.atan(Math.hypot((b.elev[i + k][j] - b.elev[i - k][j]) / dy, (b.elev[i][j + k] - b.elev[i][j - k]) / dx)) * 180) / Math.PI).toFixed(2));
    }
  };
  // Seabed temperature: bundled World Ocean Atlas profile at the site depth; in shallow water the coldest month at the surface governs.
  const deep = async () => {
    const d = site.data;
    if (!(d.depth > 0)) { if (ctx.mark.get('seabedTemp') === 'deep') { delete d.seabedTemp; delete d.tempProfile; ctx.mark.delete('seabedTemp'); ctx.mark.delete('tempProfile'); } return; }
    try {
      const { deepTemperature } = await import('../data/deep_lookup.js'), r = deepTemperature(site.lat, site.lon, d.depth);
      if (!r || r.cellKm > 450) return;
      const shallow = d.depth < 40 && (d.sstMin != null || d.sstMonthly), t = shallow ? Math.min(r.T, d.sstMin ?? Math.min(...d.sstMonthly.filter((x) => x > -3))) : r.T;
      d.seabedTemp = +t.toFixed(2); d.tempProfile = r.profile; ctx.mark.set('seabedTemp', 'deep'); ctx.mark.set('tempProfile', 'deep');
      ctx.notes.seabedTemp = `${r.vintage}, annual mean at ${Math.round(Math.min(d.depth, r.levelDepth))} m in the 2° cell ${r.cellKm} km away${shallow ? '; shallow water, so the coldest surface month is used when it is lower' : ''}${r.extrapolated ? '; the cell is shallower than the site, deepest level used' : ''}`;
      ctx.vintage.deep = r.vintage;
    } catch { /* atlas module unavailable */ }
  };
  if (opt.atlasOnly) { // no network at all: answer from the built-in atlas
    for (const src of SOURCES) { status[src.id] = { ok: false, message: 'No connection', at: new Date().toISOString() }; ctx.done.add(src.id); onStatus(src.id, 'fail', 'No connection'); }
    await atlasFill(site, status, ctx); finish(); await deep();
    return snapshot();
  }
  const national = run('place').then(() => Promise.all([run('fx'), run('economy').then(() => run('energy')), run('fiscal'), run('rates'), run('power')])); // these need the country
  await Promise.all([national, run('prices'), run('costs'), run('markets'), ...['weather', 'marine', 'salinity', 'climate'].map(run), run('bathy').then(deep).then(() => onData(snapshot()))]);
  await atlasFill(site, status, ctx); finish(); await deep();
  return snapshot();
}

/**
 * Coarse stand-by value for a point: nearest ocean cell of the 4° climatology shipped with the core (searched outwards
 * up to three cells). Returns null far inland. Used only if the modules of the full built-in atlas cannot be loaded.
 */
export function atlasSite(lat, lon) {
  const A = ATLAS, i0 = Math.round((lat - A.lat0) / A.step), j0 = Math.round((((lon - A.lon0) % 360) + 360) % 360 / A.step);
  let best = null;
  for (let r = 0; r <= 3 && !best; r++) for (let di = -r; di <= r; di++) for (let dj = -r; dj <= r; dj++) {
    if (Math.max(Math.abs(di), Math.abs(dj)) !== r) continue;
    const i = i0 + di, j = (((j0 + dj) % A.nlon) + A.nlon) % A.nlon;
    if (i < 0 || i >= A.nlat) continue;
    const k = i * A.nlon + j;
    if (A.sal[k] < 0) continue;
    const d = Math.hypot(di, dj * Math.cos((lat * Math.PI) / 180));
    if (!best || d < best.d) best = { d, k };
  }
  if (!best) return null;
  const k = best.k, t0 = A.t0[k] / 10, ta = A.ta[k] / 10, tp = A.tp[k] / 10, month = new Date().getUTCMonth();
  const sstMonthly = Array.from({ length: 12 }, (_, m) => +(t0 + ta * Math.cos((2 * Math.PI * (m - tp)) / 12)).toFixed(2));
  return { salinity: A.sal[k] / 10, sst: sstMonthly[month], sstMin: +(t0 - ta).toFixed(2), sstMax: +(t0 + ta).toFixed(2), sstMonthly, salinityMonthly: new Array(12).fill(A.sal[k] / 10) };
}
/** Clear-sky-based estimate of the long-term mean solar irradiation (kWh/m²·d) from latitude alone (used outside the atlas grid, beyond 60° S / 80° N). */
export const solarEstimate = (lat) => +clamp(6.4 * Math.cos((Math.abs(lat) * Math.PI) / 180) ** 1.15 + 0.6, 1.5, 6.8).toFixed(2);

/** Regional salinity estimate, used only where neither the live service nor the built-in atlas has a value (e.g. the Caspian Sea, far inland). */
export function regionalSalinity(lat, lon) {
  const box = (a, b, c, d) => lat >= a && lat <= b && lon >= c && lon <= d;
  if (box(23.5, 30.5, 47.5, 56.5)) return 42;
  if (box(12, 30, 32, 43.5)) return 40;
  if (box(30, 46, -6, 36.5)) return 38.3;
  if (box(53, 66, 9, 30)) return 7.5;
  if (box(40.5, 47, 27, 42)) return 18;
  if (box(36, 47, 46.5, 55)) return 12.5;
  return Math.abs(lat) < 35 ? 35.5 : 34.3;
}

/**
 * Live literature evidence for a decision topic: recent, well-cited peer-reviewed work from OpenAlex,
 * with Crossref as a fallback. Returns [{ title, year, venue, cited, url }] (plain text only).
 */
export async function evidence(topic, n = 5) {
  const q = encodeURIComponent(String(topic).slice(0, 160));
  const year = new Date().getUTCFullYear() - 6;
  const safeUrl = (doi) => { const d = txt(String(doi || '').replace(/^https?:\/\/(dx\.)?doi\.org\//i, ''), 200); return /^10\.\d{4,9}\/\S+$/.test(d) ? 'https://doi.org/' + encodeURI(d) : ''; };
  try {
    const j = await getJSON(`https://api.openalex.org/works?filter=title_and_abstract.search:${q},from_publication_date:${year}-01-01,type:article|review,cited_by_count:%3E4&sort=relevance_score:desc&per-page=${n}&select=title,publication_year,doi,cited_by_count,primary_location`, 15000);
    const out = (j.results || []).map((w) => ({ title: txt(w.title, 260), year: num(w.publication_year), venue: txt(w.primary_location?.source?.display_name, 120), cited: num(w.cited_by_count), url: safeUrl(w.doi) })).filter((w) => w.title);
    if (out.length) return { source: 'OpenAlex', items: out };
  } catch { /* fall through to Crossref */ }
  const j = await getJSON(`https://api.crossref.org/works?query=${q}&rows=${n}&select=title,DOI,issued,container-title,is-referenced-by-count&filter=from-pub-date:${year}`, 15000);
  return { source: 'Crossref', items: (j.message?.items || []).map((w) => ({ title: txt(w.title?.[0], 260), year: num(w.issued?.['date-parts']?.[0]?.[0]), venue: txt(w['container-title']?.[0], 120), cited: num(w['is-referenced-by-count']), url: safeUrl(w.DOI) })).filter((w) => w.title) };
}
