// Live-data test: node tests/live.test.mjs   (Node 20+, needs the network for the live part)
//  1. every connector of js/core/live.js is called against the real service for two sites and its fields are range-checked;
//  2. the offline path (fetchSite … { atlasOnly: true }) must fill the same economic fields from the bundled snapshot, with dates;
//  3. every address the connectors requested is checked for the cross-origin header a browser needs
//     (Access-Control-Allow-Origin: * or the page origin, on the final response of any redirect chain, HTTPS, allow-listed host).
// A source that cannot be reached (no network, time-out, HTTP 429/5xx, daily request limit) is reported as
// "skipped (unreachable)"; a source that answers with unusable data fails the test.
import { readFileSync } from 'node:fs';
import { SOURCES, ATLAS_LABELS, fetchSource, fetchSite, feeds, costStats, gasRegion } from '../js/core/live.js';
import { FISCAL, fiscalOf } from '../js/data/fiscal.js';
import { PRICES, MARKETS, COSTS, CORPORATE_TAX, RATES, ELECTRICITY, ISO3 } from '../js/data/prices.js';

const ORIGIN = 'https://samuelakosaonyejekwe.github.io';
let fails = 0, skips = 0, passes = 0;
const ok = (cond, msg) => { if (cond) passes++; else { fails++; console.log('   ✗ ' + msg); } return cond; };
const between = (v, a, b) => typeof v === 'number' && Number.isFinite(v) && v >= a && v <= b;
const monthsOld = (p) => { const s = String(p || ''), t = Date.parse(/^\d{4}$/.test(s) ? s + '-12-31' : /^\d{4}-\d\d$/.test(s) ? s + '-01' : /^\d{4}-S[12]$/.test(s) ? `${s.slice(0, 4)}-${s.endsWith('1') ? '06' : '12'}-01` : s); return Number.isFinite(t) ? (Date.now() - t) / (30.44 * 24 * 3600e3) : NaN; };

// ---- network recorder: sends the page origin, notes the answer of every address the connectors use ------------------
const seen = new Map(); let netFailures = 0, calls = 0;
const realFetch = globalThis.fetch;
globalThis.fetch = async (url, opt = {}) => {
  calls++;
  const u = new URL(String(url)), key = `${opt.method || 'GET'} ${u.origin}${u.pathname.replace(/\/country\/[^/]+\/indicator\/.*$/, '/country/…/indicator/…')}`;
  try {
    const r = await realFetch(url, { ...opt, headers: { ...(opt.headers || {}), Origin: ORIGIN } });
    const acao = r.headers.get('access-control-allow-origin'), fin = new URL(r.url || String(url));
    if (!seen.has(key) || r.ok) seen.set(key, { status: r.status, acao, finalHost: fin.hostname, https: fin.protocol === 'https:', redirected: r.redirected, sample: u.href.slice(0, 150) });
    if (r.status === 429 || r.status >= 500 || r.status === 403) netFailures++;
    return r;
  } catch (e) { netFailures++; if (!seen.has(key)) seen.set(key, { status: 0, acao: null, error: e.cause?.code || e.name, sample: u.href.slice(0, 150) }); throw e; }
};
const unreachable = (e, failed) => failed > 0 || e?.name === 'AbortError' || e?.name === 'TypeError' || /^HTTP (403|429|5\d\d)|daily threshold|not processed|fetch failed/i.test(e?.message || '');

/** Run one connector; returns its data, or null when skipped/failed (already counted). `expectError` = a refusal that is the correct answer. */
async function source(id, label, lat, lon, site, check, expectError) {
  const before = netFailures;
  try {
    const r = await fetchSource(id, lat, lon, structuredClone(site));
    if (expectError) { ok(false, `${id} @ ${label}: expected "${expectError}" but got data`); return null; }
    const n = fails; check(r.data || {}, r.meta || {});
    console.log(` ${fails === n ? '✓' : '✗'} ${id.padEnd(8)} ${label}`);
    return r;
  } catch (e) {
    if (expectError && String(e.message).includes(expectError)) { passes++; console.log(` ✓ ${id.padEnd(8)} ${label}: correctly reports “${e.message}”`); return null; }
    if (unreachable(e, netFailures - before)) { skips++; console.log(` – ${id.padEnd(8)} ${label}: skipped (unreachable) — ${e.message}`); return null; }
    ok(false, `${id} @ ${label}: reachable but unusable — ${e.message}`); return null;
  }
}

const SITES = [
  { label: 'Nigeria offshore 3.5 N 5.75 E', lat: 3.5, lon: 5.75, site: { countryCode: 'NG', country: 'Nigeria', data: { iso3: 'NGA', currency: 'NGN' } }, tax: 30, industrialPower: false },
  { label: 'Norway 61.2 N 1.8 E', lat: 61.2, lon: 1.8, site: { countryCode: 'NO', country: 'Norway', data: { iso3: 'NOR', currency: 'NOK' } }, tax: 22, industrialPower: true },
];
const year = new Date().getUTCFullYear();
const CHECKS = {
  place: (d, m) => ok(typeof m.countryCode === 'string' && typeof d.currency === 'string', 'place: country code and currency are text'),
  weather: (d) => { ok(between(d.airTemp, -70, 60), `weather: air temperature ${d.airTemp}`); ok(between(d.windSpeed, 0, 80), `weather: wind speed ${d.windSpeed}`); },
  marine: (d) => { ok(d.sst == null || between(d.sst, -3, 38), `marine: sea-surface temperature ${d.sst}`); ok(d.waveHeight == null || between(d.waveHeight, 0, 25), `marine: wave height ${d.waveHeight}`); },
  bathy: (d) => { ok(between(d.depth, 0, 11500), `bathy: depth ${d.depth}`); ok(d.bathy?.elev?.length === 57, 'bathy: 57 × 57 grid'); },
  salinity: (d) => ok(between(d.salinity, 1, 45) && d.salinityMonthly?.length === 12, `salinity: ${d.salinity}`),
  economy: (d) => { ok(between(d.inflation, -30, 2000) && /^\d{4}$/.test(d.inflationYear), `economy: inflation ${d.inflation} (${d.inflationYear})`); ok(/^[A-Z]{3}$/.test(d.iso3 || ''), 'economy: ISO-3 code'); },
  fx: (d) => ok(d.fxPerUSD > 0 && typeof d.fxDate === 'string', `fx: ${d.fxPerUSD}`),
  energy: (d) => ok(between(d.gridCarbon, 0, 1.5) && /^\d{4}$/.test(String(d.gridCarbonYear)), `energy: grid carbon ${d.gridCarbon} (${d.gridCarbonYear})`),
  prices: (d) => { ok(between(d.oilPrice, 5, 400) && /^\d{4}-\d\d-\d\d$/.test(d.oilPriceDate), `prices: Brent ${d.oilPrice} (${d.oilPriceDate})`); ok(between(d.gasPrice, 0.2, 100), `prices: Henry Hub ${d.gasPrice}`); },
  climate: (d) => ok(between(d.ghiAnnual, 0.2, 10), `climate: solar ${d.ghiAnnual}`),
  fiscal: (d, m, s) => { ok(between(d.corporateTaxRate, 0, 60), `fiscal: corporate tax ${d.corporateTaxRate}`); ok(Math.abs(d.corporateTaxRate - s.tax) <= 3, `fiscal: ${s.label} statutory rate ${d.corporateTaxRate}, expected about ${s.tax}`); ok(+d.corporateTaxYear >= 2022 && +d.corporateTaxYear <= year, `fiscal: year ${d.corporateTaxYear}`); ok(/OECD|Tax Foundation/.test(d.corporateTaxSource || ''), 'fiscal: source named'); },
  costs: (d) => {
    ok(between(d.costIndex, 60, 400), `costs: index ${d.costIndex}`); ok(/^\d{4} = 100$/.test(d.costIndexBase || ''), `costs: base stated (${d.costIndexBase})`); ok(monthsOld(d.costIndexDate) < 9, `costs: latest month ${d.costIndexDate}`);
    ok(d.costIndexSeries?.t?.length >= 120 && d.costIndexSeries.t.length === d.costIndexSeries.v.length, `costs: series of ${d.costIndexSeries?.t?.length} months (at least 120)`); ok(between(d.costEscalation, -15, 25), `costs: escalation ${d.costEscalation} %/y`);
    ok(between(d.steelIndex, 40, 500), `costs: steel index ${d.steelIndex}`); ok(between(d.usCpi, 80, 300) && between(d.usInflation5y, -5, 20) && d.usCpiSeries?.t?.length >= 100, `costs: US CPI ${d.usCpi}, ${d.usInflation5y} %/y`);
  },
  markets: (d) => {
    ok(between(d.gasPriceEurope, 0.5, 150) && monthsOld(d.gasPriceEuropeDate) < 24, `markets: European gas ${d.gasPriceEurope} (${d.gasPriceEuropeDate})`);
    if (/IMF/.test(d.gasPriceEuropeSource || '')) { ok(between(d.gasPriceAsia, 0.5, 150) && monthsOld(d.gasPriceAsiaDate) < 6, `markets: Asian LNG ${d.gasPriceAsia} (${d.gasPriceAsiaDate})`); ok(between(d.steelPrice, 20, 400) && d.steelPriceUnit === 'US$/t' && /iron ore/.test(d.steelPriceWhat), `markets: iron ore ${d.steelPrice}`); ok(d.gasPriceEuropeSeries?.t?.length >= 100, 'markets: ten-year monthly series'); }
    else console.log('     (IMF series not reachable: annual European price from the fallback only)');
  },
  rates: (d) => { const ks = ['bondYield', 'treasuryBillYield', 'policyRate'].filter((k) => d[k] != null); ok(ks.length > 0, 'rates: at least one rate'); for (const k of ks) ok(between(d[k], -2, 90) && monthsOld(d[k + 'Date']) < 37 && d[k + 'Source'], `rates: ${k} ${d[k]} (${d[k + 'Date']})`); },
  power: (d) => ok(between(d.electricityPrice, 0.01, 0.8) && monthsOld(d.electricityPriceDate) < 30 && /Eurostat|EIA/.test(d.electricityPriceSource), `power: ${d.electricityPrice} $/kWh (${d.electricityPriceDate})`),
};

console.log('— bundled data —');
{ // the fiscal table and the snapshot are internally sound
  for (const [c, f] of Object.entries(FISCAL)) {
    ok(/^[A-Z]{2}$/.test(c) && f.source?.citation && /^https:\/\//.test(f.source.url) && /^\d{4}-\d\d-\d\d$/.test(f.source.retrieved), `fiscal table ${c}: citation, address and retrieval date`);
    for (const k of ['royalty', 'petroleumTax', 'corporateTax', 'marginalTake', 'headline']) ok(f[k] === null || between(f[k], 0, 95), `fiscal table ${c}: ${k} = ${f[k]}`);
    ok([null, 'tax-royalty', 'psc', 'mixed'].includes(f.regime) && between(f.year, 2020, year), `fiscal table ${c}: regime and year`);
    ok(f.headline === null || (typeof f.basis === 'string' && f.basis.length > 10), `fiscal table ${c}: the headline rate says what it is`);
  }
  ok(fiscalOf('no').marginalTake === 78 && fiscalOf('GB').headline === 78 && fiscalOf('ZZ') === null, 'fiscal table: look-up');
  ok(Object.keys(CORPORATE_TAX.rows).length > 100 && Object.keys(RATES.rows).length > 40 && Object.keys(ELECTRICITY.rows).length > 25 && Object.keys(ISO3).length > 180, 'snapshot: coverage');
  for (const s of [PRICES.date, MARKETS.retrieved, COSTS.retrieved, CORPORATE_TAX.retrieved, RATES.retrieved, ELECTRICITY.retrieved]) ok(/^\d{4}-\d\d-\d\d$/.test(s) && monthsOld(s) < 18, `snapshot: dated (${s})`);
  const c = costStats(COSTS.cost, COSTS.steel, COSTS.cpi); CHECKS.costs({ ...c, costIndexDate: c.costIndexDate && monthsOld(c.costIndexDate) < 24 ? new Date().toISOString().slice(0, 7) : c.costIndexDate });
  ok(gasRegion(-90) === 'americas' && gasRegion(5) === 'europe' && gasRegion(120) === 'asia', 'regional gas benchmark by longitude');
  for (const k of ['corporateTaxRate', 'costIndex', 'costEscalation', 'steelIndex', 'usCpi', 'gasPriceEurope', 'gasPriceAsia', 'gasPriceRegional', 'steelPrice', 'bondYield', 'policyRate', 'treasuryBillYield', 'electricityPrice']) ok(typeof ATLAS_LABELS[k] === 'string', `plain-language label for ${k}`);
  console.log(` ${fails ? '✗' : '✓'} fiscal table (${Object.keys(FISCAL).length} countries) and snapshot`);
}

console.log('— live connectors (two sites) —');
const live = {};
for (const s of SITES) {
  const res = await Promise.all(SOURCES.map((src) => source(src.id, s.label, s.lat, s.lon, s.site, (d, m) => CHECKS[src.id](d, m, s), src.id === 'power' && !s.industrialPower ? 'No open industrial tariff series' : null)));
  live[s.label] = Object.assign({}, ...res.filter(Boolean).map((r) => r.data));
}
console.log('— fallback services —');
for (const [name, fn, check] of [
  ['Tax Foundation corporate tax (NO)', () => feeds.corporateTaxTF('NO'), (r) => ok(r && between(r[0], 15, 30) && r[1] >= 2022, `Tax Foundation: ${r}`)],
  ['OECD long-term interest rate (NOR)', () => feeds.oecdBond('NOR'), (r) => ok(r.NOR && between(r.NOR[0], -1, 20) && monthsOld(r.NOR[1]) < 12, `OECD bond yield: ${r.NOR}`)],
  ['BIS policy rate (NO)', () => feeds.bisPolicy('NO'), (r) => ok(r && between(r[0], -1, 30) && monthsOld(r[1]) < 12, `BIS policy rate: ${r}`)],
  ['Our World in Data European gas', () => feeds.owidGasEurope(), (r) => ok(r && between(r[0], 0.5, 150) && r[1] >= year - 2, `annual TTF: ${r}`)],
  ['BLS mirror (DBnomics)', () => feeds.blsMirror(), (r) => ok(Object.values(r).every((x) => x.length > 200), 'mirror: three long series')],
  ['US industrial electricity (EIA)', () => feeds.eiaPower(), (r) => ok(between(r[0], 0.03, 0.3) && monthsOld(r[1]) < 12, `EIA: ${r}`)],
]) {
  const before = netFailures, n = fails;
  try { check(await fn()); console.log(` ${fails === n ? '✓' : '✗'} ${name}`); }
  catch (e) { if (unreachable(e, netFailures - before)) { skips++; console.log(` – ${name}: skipped (unreachable) — ${e.message}`); } else ok(false, `${name}: reachable but unusable — ${e.message}`); }
}

console.log('— offline path (built-in atlas and snapshot, no network) —');
const NEW = ['taxRate', 'taxRateYear', 'taxRateSource', 'taxRateBasis', 'corporateTaxRate', 'corporateTaxYear', 'corporateTaxSource', 'costIndex', 'costIndexDate', 'costIndexBase', 'costIndexSeries', 'costEscalation', 'steelIndex', 'steelIndexDate', 'usCpi', 'usCpiDate', 'usCpiSeries', 'usInflation5y',
  'gasPriceEurope', 'gasPriceEuropeDate', 'gasPriceAsia', 'gasPriceAsiaDate', 'gasPriceRegional', 'gasPriceRegionalDate', 'gasPriceRegionalName', 'steelPrice', 'steelPriceDate', 'steelPriceUnit', 'steelPriceWhat', 'policyRate', 'policyRateDate', 'electricityPrice', 'electricityPriceSource'];
for (const [label, lat, lon, code, liveKey] of [['Nigeria offshore 3.5 N 5.75 E', 3.5, 5.75, 'NG', SITES[0].label], ['Norway coast 60.4 N 5.0 E', 60.4, 5.0, 'NO', SITES[1].label], ['Norway 61.2 N 1.8 E (open sea, no country within 100 km)', 61.2, 1.8, '', null]]) {
  const before = calls, n = fails, s = await fetchSite(lat, lon, () => {}, () => {}, { atlasOnly: true }), d = s.data, AF = new Set(d.atlasFields), notes = d.atlasNotes || {};
  ok(calls === before, `${label}: the offline path made ${calls - before} network requests`);
  ok(s.countryCode === code, `${label}: country ${s.countryCode || '(none)'}, expected ${code || '(none)'}`);
  const global = ['costIndex', 'costEscalation', 'steelIndex', 'usCpi', 'gasPriceEurope', 'gasPriceAsia', 'gasPriceRegional', 'steelPrice', 'oilPrice', 'gasPrice'];
  for (const k of global) ok(d[k] != null && AF.has(k), `${label}: ${k} filled and marked built-in atlas`);
  for (const k of ['costIndex', 'gasPriceEurope', 'steelPrice']) ok(/\d{4}-\d\d/.test(notes[k] || ''), `${label}: ${k} carries its date in the note (${notes[k]})`);
  CHECKS.costs({ ...d, costIndexDate: new Date().toISOString().slice(0, 7) }); ok(/^\d{4}-\d\d$/.test(d.costIndexDate), `${label}: cost index date ${d.costIndexDate}`);
  ok(between(d.gasPriceEurope, 0.5, 150) && between(d.gasPriceAsia, 0.5, 150) && between(d.steelPrice, 20, 400), `${label}: snapshot prices in range`);
  if (code) {
    const f = fiscalOf(code);
    ok(d.taxRate === f.headline && d.taxRateSource === f.source.citation && +d.taxRateYear === f.year && d.taxRateBasis, `${label}: tax rate ${d.taxRate} from the sourced fiscal table, with source and year`);
    ok(between(d.corporateTaxRate, 0, 60) && AF.has('corporateTaxRate') && +d.corporateTaxYear >= 2024, `${label}: corporate tax ${d.corporateTaxRate} (${d.corporateTaxYear}) from the snapshot`);
    ok(['bondYield', 'policyRate', 'treasuryBillYield'].some((k) => d[k] != null && AF.has(k) && /^\d{4}-\d\d$/.test(d[k + 'Date'])), `${label}: an interest rate from the snapshot, dated`);
    ok(between(d.electricityPrice, 0.01, 0.8) && AF.has('electricityPrice'), `${label}: electricity price ${d.electricityPrice}`);
    if (code === 'NO') ok(/Eurostat/.test(d.electricityPriceSource) && /^\d{4}-S[12]$/.test(d.electricityPriceDate) && !/ESTIMATE/.test(notes.electricityPrice), `${label}: electricity price is the sourced Eurostat value (${d.electricityPriceDate})`);
    if (code === 'NG') ok(/^ESTIMATE/.test(notes.electricityPrice || '') && /estimate/.test(d.electricityPriceSource), `${label}: electricity price is labelled an estimate`);
    const lv = liveKey ? live[liveKey] : null;
    if (lv && Object.keys(lv).length) { // every new field the live sources delivered must also exist offline
      const liveAll = { ...lv, taxRate: 1, taxRateYear: 1, taxRateSource: 1, taxRateBasis: 1, ...(lv.gasPriceEurope != null ? { gasPriceRegional: 1, gasPriceRegionalDate: 1, gasPriceRegionalName: 1 } : {}) };
      const missing = NEW.filter((k) => liveAll[k] != null && d[k] == null);
      ok(!missing.length, `${label}: fields delivered live but missing offline → ${missing.join(', ')}`);
    }
  } else ok(d.taxRate == null && d.corporateTaxRate == null, `${label}: no national figures without a country`);
  console.log(` ${fails === n ? '✓' : '✗'} ${label}: ${d.atlasFields.length} fields from the atlas and snapshot`);
}

console.log('— cross-origin readability of every address used (Origin: ' + ORIGIN + ') —');
const csp = readFileSync(new URL('../index.html', import.meta.url), 'utf8').match(/connect-src([^;]*)/)?.[1] || '', sw = readFileSync(new URL('../tools/sw.template.js', import.meta.url), 'utf8').match(/DATA_HOSTS = \[([^\]]*)\]/)?.[1] || '';
for (const [key, v] of [...seen].sort()) {
  const reach = v.status >= 200 && v.status < 300, cors = v.acao === '*' || v.acao === ORIGIN;
  const verdict = !v.status ? `unreachable (${v.error})` : !reach ? `HTTP ${v.status} — not judged` : cors && v.https ? 'readable' : 'BLOCKED for browsers';
  console.log(` ${!reach ? '–' : cors ? '✓' : '✗'} ${key}  →  ${v.status || '—'}  allow-origin: ${v.acao ?? 'none'}${v.redirected ? '  (redirected to ' + v.finalHost + ')' : ''}  ${verdict}`);
  if (reach) { ok(cors && v.https, `${key}: no usable Access-Control-Allow-Origin on the final response (${v.sample})`); ok(csp.includes('https://' + v.finalHost) && sw.includes(`'${v.finalHost}'`), `${v.finalHost}: missing from the connect-src list of index.html or from DATA_HOSTS of the service worker`); }
}

console.log(`\n${passes} checks passed, ${skips} sources skipped (unreachable), ${fails} failed.`);
process.exit(fails ? 1 : 0);
