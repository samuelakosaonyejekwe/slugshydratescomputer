// Live-data test: node tests/live.test.mjs   (Node 20+, needs the network for the live part)
//  1. the bundled tables are sound: fiscal terms (ranges, terrain by water depth, merge notes), snapshot, national
//     electricity table, upstream carbon intensity;
//  2. the shared cache (functions/api/feed.js) is driven with stubs: allow-list, GET only, one upstream request per
//     period, last good answer when the upstream fails, size cap, shared store;
//  3. the order shared cache → service itself → bundled snapshot is checked with stubbed transports;
//  4. every connector of js/core/live.js is called against the real service for two sites (shared cache switched off)
//     and its fields are range-checked; then the same economic connectors run through the handler of the shared cache
//     in-process against the real upstream services;
//  5. the offshore jurisdiction look-up (exclusive economic zone, else nearest country) is checked on six basins;
//  6. the offline path (fetchSite … { atlasOnly: true }) must fill the same economic fields from the bundled snapshot, with dates;
//  7. every address the connectors requested is checked for the cross-origin header a browser needs
//     (Access-Control-Allow-Origin: * or the page origin, on the final response of any redirect chain, HTTPS, allow-listed host).
// A source that cannot be reached (no network, time-out, HTTP 429/5xx, daily request limit) is reported as
// "skipped (unreachable)"; a source that answers with unusable data fails the test.
import { readFileSync } from 'node:fs';
import { SOURCES, ATLAS_LABELS, EDGE, EDGE_INFO, fetchSource, fetchSite, feeds, feedInfo, costStats, gasRegion, setEdgeTransport, offshoreState, nearestCountry } from '../js/core/live.js';
import { FISCAL, FISCAL_TERMS, fiscalOf } from '../js/data/fiscal.js';
import { FISCAL_SUPPLEMENT as ECON_FISCAL } from '../js/data/costbasis.js'; // the economics suite's former table, now a supplement of the shared one
import { PRICES, MARKETS, COSTS, CORPORATE_TAX, RATES, ELECTRICITY, ELECTRICITY_NATIONAL, UPSTREAM_CI, ISO3 } from '../js/data/prices.js';
import { handle, FEEDS } from '../functions/api/feed.js';

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
const edgeOff = () => setEdgeTransport(() => { throw new TypeError('shared cache switched off for this part of the test'); });
const unreachable = (e, failed) => failed > 0 || e?.name === 'AbortError' || e?.name === 'TypeError' || /^HTTP (403|429|5\d\d)|daily threshold|not processed|fetch failed|bundled copy is newer|not reachable/i.test(e?.message || '');

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
  { label: 'Norway 61.2 N 2.2 E', lat: 61.2, lon: 2.2, site: { countryCode: 'NO', country: 'Norway', data: { iso3: 'NOR', currency: 'NOK' } }, tax: 22, industrialPower: true },
];
const year = new Date().getUTCFullYear();
const CHECKS = {
  place: (d, m, s) => { ok(typeof m.countryCode === 'string' && typeof d.currency === 'string', 'place: country code and currency are text'); ok(m.countryCode === s.site.countryCode, `place: ${s.label} → ${m.countryCode || '(none)'}, expected ${s.site.countryCode}`); if (d.jurisdictionBasis) ok(/exclusive economic zone|nearest country/.test(d.jurisdictionBasis) && d.jurisdictionSource, `place: jurisdiction ${d.jurisdictionBasis}`); },
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
  power: (d) => ok(between(d.electricityPrice, 0.01, 0.8) && monthsOld(d.electricityPriceDate) < 30 && /Eurostat|EIA|Energy Security|ElCom|Singapore/.test(d.electricityPriceSource) && d.electricityPriceKind === 'industrial', `power: ${d.electricityPrice} $/kWh (${d.electricityPriceDate})`),
};

console.log('— bundled data —');
{ // the fiscal table and the snapshot are internally sound
  const src = (x) => x?.citation && /^https:\/\//.test(x.url) && /^\d{4}-\d\d-\d\d$/.test(x.retrieved);
  for (const [c, f] of Object.entries(FISCAL)) {
    ok(/^[A-Z]{2}$/.test(c) && src(f.source) && ['source2', 'source3', 'source4'].every((k) => f[k] === undefined || src(f[k])), `fiscal table ${c}: citation, address and retrieval date`);
    for (const k of FISCAL_TERMS) {
      ok(f[k] === null || between(f[k], 0, 95), `fiscal table ${c}: ${k} = ${f[k]}`);
      const lo = f[k + 'Low'], hi = f[k + 'High'];
      if (lo !== undefined || hi !== undefined) ok(between(lo, 0, 95) && between(hi, 0, 95) && lo < hi && Math.abs(f[k] - (lo + hi) / 2) < 1e-9, `fiscal table ${c}: ${k} ${f[k]} is the midpoint of ${lo}–${hi}`);
    }
    ok([null, 'tax-royalty', 'psc', 'mixed'].includes(f.regime) && between(f.year, 2019, year), `fiscal table ${c}: regime and year`);
    ok(f.headline === null || (typeof f.basis === 'string' && f.basis.length > 10), `fiscal table ${c}: the headline rate says what it is`);
    if (f.headlineLow !== undefined) ok(/midpoint/.test(f.basis) && f.basis.includes(`${f.headlineLow}–${f.headlineHigh}%`), `fiscal table ${c}: the basis names the midpoint and the range (${f.basis})`);
    for (const t of f.terrain || []) ok(typeof t.name === 'string' && (t.minDepth === undefined || t.minDepth >= 0) && (t.maxDepth === undefined || t.maxDepth >= (t.minDepth ?? 0)) && typeof f.terrainBasis === 'string' && f.terrainBasis.length > 20, `fiscal table ${c}: terrain class ${t.name}`);
  }
  ok(fiscalOf('no').marginalTake === 78 && fiscalOf('GB').headline === 78 && fiscalOf('ZZ') === null, 'fiscal table: look-up');
  // ranges: midpoint as the headline, both bounds kept (owner's decision)
  for (const [c, k, lo, hi] of [['SA', 'headline', 50, 85], ['SA', 'petroleumTax', 50, 85], ['VN', 'headline', 25, 50], ['CO', 'headline', 35, 50], ['KZ', 'royalty', 5, 18], ['GH', 'royalty', 3, 12.5], ['GA', 'costOilCap', 65, 75], ['NL', 'royalty', 0, 7]]) { const f = fiscalOf(c); ok(f[k] === (lo + hi) / 2 && f[k + 'Low'] === lo && f[k + 'High'] === hi, `fiscal range ${c} ${k}: ${f[k]} (${f[k + 'Low']}–${f[k + 'High']}), expected the midpoint of ${lo}–${hi}`); }
  ok(fiscalOf('KZ').headline === 20 && fiscalOf('KZ').headlineLow === undefined && fiscalOf('NO').royaltyLow === undefined, 'fiscal ranges: single figures carry no bounds');
  // terms by water depth
  const ng = (d) => fiscalOf('NG', { waterDepth: d });
  ok(ng(0).royalty === 15 && ng(0).headline === 60 && ng(0).terrainName === 'onshore', 'Nigeria onshore: royalty 15%, 60%');
  ok(ng(150).royalty === 12.5 && ng(200).headline === 60 && /shallow/.test(ng(200).terrainName), 'Nigeria shallow water (to 200 m): royalty 12.5%, 60%');
  ok(ng(201).royalty === 7.5 && ng(1500).headline === 30 && ng(1500).petroleumTax === 0 && /deep offshore/.test(ng(1500).terrainName) && /hydrocarbon tax/.test(ng(1500).basis) && ng(1500).terrainBasis.includes('200 m'), 'Nigeria deep offshore (beyond 200 m): royalty 7.5%, 30%, basis explains');
  ok(fiscalOf('NG').headline === 60 && fiscalOf('NG').terrainName === undefined && fiscalOf('NG', {}).headline === 60 && fiscalOf('NG', { waterDepth: NaN }).headline === 60 && fiscalOf('NG', { waterDepth: -5 }).headline === 60, 'Nigeria without a usable depth: the row as tabled');
  ok([0, 300, 900, 2000].map((d) => fiscalOf('IN', { waterDepth: d }).royalty).join() === '12.5,7.5,5,2', 'India: royalty by depth class (400 m and 1 500 m isobaths)');
  ok(fiscalOf('SN', { waterDepth: 0 }).costOilCap === 55 && fiscalOf('SN', { waterDepth: 0 }).costOilCapLow === undefined && fiscalOf('SN', { waterDepth: 900 }).costOilCapLow === 60 && fiscalOf('SN', { waterDepth: 900 }).costOilCapHigh === 70 && fiscalOf('SN').costOilCapLow === 55, 'Senegal: cost-oil cap onshore 55%, offshore 60–70%');
  ok(fiscalOf('TT', { waterDepth: 900 }).headline === 35 && fiscalOf('TT', { waterDepth: 100 }).headline === 55 && fiscalOf('NO', { waterDepth: 300 }) === fiscalOf('NO'), 'Trinidad deep water 35%; a country without terrain classes is unchanged');
  // single source: every country of the economics table is here, with the fields that suite needs
  const need = ['regime', 'royalty', 'petroleumTax', 'corporateTax', 'marginalTake', 'costOilCap', 'profitSplit'];
  for (const c of Object.keys(ECON_FISCAL)) ok(FISCAL[c] && need.every((k) => k in FISCAL[c]), `economics table country ${c} is in the fiscal table with all fields`);
  for (const [c, e] of Object.entries(ECON_FISCAL)) { const f = FISCAL[c], m = e.model; if (!f) continue; const same = (a, b) => a == null || b == null || Math.abs(a - b) < 0.01; if (!(same(m.royalty, f.royalty) && same(m.costOilCap, f.costOilCap) && same(m.profitSplit, f.profitSplit)) || (f.headline != null && !same(m.taxRate, f.headline))) ok(typeof f.merged === 'string' && f.merged.length > 30, `fiscal table ${c}: differs from the economics table (royalty ${m.royalty}/${f.royalty}, tax ${m.taxRate}/${f.headline}, cost oil ${m.costOilCap}/${f.costOilCap}) without a note on the choice`); }
  // national electricity table and upstream carbon intensity
  for (const [c, r] of Object.entries(ELECTRICITY_NATIONAL.rows)) ok(/^[A-Z]{2}$/.test(c) && between(r.usd, 0.01, 0.8) && /^\d{4}$/.test(r.period) && ['industrial', 'proxy'].includes(r.kind) && r.source.length > 20 && /^https:\/\//.test(r.url) && r.licence && r.basis.length > 20 && (r.kind === 'proxy' ? /not an industrial tariff/.test(r.basis) : between(r.price, 0.001, 1000) && /^[A-Z]{3}$/.test(r.currency) && Math.abs(r.price / r.fxRate - r.usd) < 0.001) && !ELECTRICITY.rows[c], `national electricity row ${c}`);
  ok(/doi/.test(UPSTREAM_CI.citation) && UPSTREAM_CI.licence === 'CC BY 4.0' && Object.entries(UPSTREAM_CI.rows).every(([c, v]) => /^[A-Z]{2}$/.test(c) && v[1] <= v[0] && v[0] <= v[2] && between(v[0], 1, 200)), 'upstream carbon intensity table: citation, licence, low ≤ value ≤ high');
  ok(Object.keys(CORPORATE_TAX.rows).length > 100 && Object.keys(RATES.rows).length > 40 && Object.keys(ELECTRICITY.rows).length > 25 && Object.keys(ISO3).length > 180, 'snapshot: coverage');
  for (const s of [PRICES.date, MARKETS.retrieved, COSTS.retrieved, CORPORATE_TAX.retrieved, RATES.retrieved, ELECTRICITY.retrieved, ELECTRICITY_NATIONAL.retrieved, UPSTREAM_CI.retrieved]) ok(/^\d{4}-\d\d-\d\d$/.test(s) && monthsOld(s) < 18, `snapshot: dated (${s})`);
  const c = costStats(COSTS.cost, COSTS.steel, COSTS.cpi); CHECKS.costs({ ...c, costIndexDate: c.costIndexDate && monthsOld(c.costIndexDate) < 24 ? new Date().toISOString().slice(0, 7) : c.costIndexDate });
  ok(gasRegion(-90) === 'americas' && gasRegion(5) === 'europe' && gasRegion(120) === 'asia', 'regional gas benchmark by longitude');
  for (const k of ['upstreamCarbonIntensity', 'royaltyRate', 'costOilCap', 'profitSplit', 'jurisdictionBasis', 'corporateTaxRate', 'costIndex', 'costEscalation', 'steelIndex', 'usCpi', 'gasPriceEurope', 'gasPriceAsia', 'gasPriceRegional', 'steelPrice', 'bondYield', 'policyRate', 'treasuryBillYield', 'electricityPrice']) ok(typeof ATLAS_LABELS[k] === 'string', `plain-language label for ${k}`);
  console.log(` ${fails ? '✗' : '✓'} fiscal table (${Object.keys(FISCAL).length} countries) and snapshot`);
}

console.log('— shared cache: handler with stubbed upstream and cache —');
{
  const n0 = fails;
  let clock = Date.UTC(2026, 9, 9, 12);
  const makeCache = () => { const m = new Map(); return { m, async match(req) { const e = m.get(req.url); if (!e || e.until <= clock) return undefined; return e.res.clone(); }, async put(req, res) { const s = /s-maxage=(\d+)/.exec(res.headers.get('cache-control') || ''); if (res.headers.get('vary') === '*') throw new TypeError('Vary: *'); m.set(req.url, { until: clock + (s ? +s[1] : 0) * 1000, res }); } }; };
  const ctx = { now: () => clock, waitUntil() {} };
  const CIT = 'DATAFLOW,REF_AREA,FREQ,MEASURE,TIME_PERIOD,OBS_VALUE\n' + Array.from({ length: 70 }, (_, i) => `X,${String.fromCharCode(65 + (i % 26))}${String.fromCharCode(65 + Math.floor(i / 26))}Q,A,CIT_C,2026,${10 + (i % 25)}`).join('\n') + '\nX,NOR,A,CIT_C,2026,22\nX,NOR,A,CIT_C,2025,21\n';
  let up = [], mode = 'ok';
  const upstream = async (url, opt = {}) => { up.push({ host: new URL(url).hostname, url: String(url), method: opt.method || 'GET' }); if (mode === 'down') throw new TypeError('fetch failed'); if (mode === '500') return new Response('oops', { status: 500 }); if (mode === 'huge') return new Response('x', { headers: { 'content-length': '99999999' } }); return new Response(CIT, { headers: { 'content-type': 'text/csv' } }); };
  const call = (q, cache, env = {}, init) => handle(new Request('https://hydraslug.pages.dev/api/feed' + q, init), env, ctx, upstream, cache);
  const cache = makeCache();
  let r = await call('?id=cit', cache), j = await r.json();
  ok(r.status === 200 && r.headers.get('x-feed-cache') === 'miss' && up.length === 1 && j.id === 'cit' && j.data.rows.NOR[0] === 22 && j.data.rows.NOR[1] === 2026, 'first request reads the upstream once and returns the stripped rows');
  ok(r.headers.get('access-control-allow-origin') === '*' && r.headers.get('vary') === '*' && /s-maxage=\d+/.test(r.headers.get('cache-control')) && r.headers.get('x-content-type-options') === 'nosniff', 'answer carries the cross-origin and cache headers');
  ok(Object.keys(j).sort().join() === 'data,fetched,id,refreshedEverySeconds,source,terms' && !JSON.stringify(j).includes('DATAFLOW'), 'answer holds only the stripped payload, its source and terms');
  r = await call('?id=cit&junk=1&url=https://evil.example/x&callback=a', cache);
  ok(r.headers.get('x-feed-cache') === 'hit' && up.length === 1, 'second request (with stray parameters) is served from the cache without an upstream request');
  ok(up.every((u) => u.host === 'sdmx.oecd.org') && !up.some((u) => /evil/.test(u.url)), 'no part of the request reaches an upstream address');
  clock += 3 * 24 * 3600e3; r = await call('?id=cit', cache);
  ok(r.headers.get('x-feed-cache') === 'hit' && +r.headers.get('x-feed-age') === 3 * 24 * 3600 && up.length === 1, 'within the period the stored answer is served with its age');
  clock += 5 * 24 * 3600e3; mode = 'down'; r = await call('?id=cit', cache); j = await r.json();
  ok(r.status === 200 && r.headers.get('x-feed-stale') === '1' && r.headers.get('x-feed-cache') === 'stale' && +r.headers.get('x-feed-age') === 8 * 24 * 3600 && j.data.rows.NOR[0] === 22 && up.length === 2, 'after the period, with the upstream down, the last good answer is served, marked stale with its age');
  r = await call('?id=cit', cache); ok(r.headers.get('x-feed-stale') === '1' && up.length === 2, 'a failing upstream is not asked again straight away');
  clock += 20 * 60e3; mode = 'ok'; r = await call('?id=cit', cache);
  ok(r.headers.get('x-feed-cache') === 'miss' && !r.headers.get('x-feed-stale') && up.length === 3 && +r.headers.get('x-feed-age') === 0, 'once the upstream answers again the copy is refreshed');
  // no stored answer at all
  const empty = makeCache(); mode = '500'; up = [];
  r = await call('?id=cit', empty); ok(r.status === 502 && (await r.json()).error && r.headers.get('access-control-allow-origin') === '*' && up.length === 1, 'upstream failure without a stored answer: 502 with a reason');
  r = await call('?id=cit', empty); ok(r.status === 502 && up.length === 1, 'the failure is remembered for a few minutes');
  mode = 'huge'; clock += 10 * 60e3; r = await call('?id=cit', empty); ok(r.status === 502 && /too large/.test((await r.json()).reason), 'an oversized upstream answer is refused');
  // request validation
  mode = 'ok'; up = [];
  r = await call('?id=nonsense', cache); ok(r.status === 404, 'unknown feed: 404');
  r = await call('?id=__proto__', cache); ok(r.status === 404, 'object-prototype names are not feeds');
  r = await call('?id=cit', cache, {}, { method: 'POST' }); ok(r.status === 405 && r.headers.get('allow') === 'GET', 'POST: 405');
  r = await call('?id=cit', cache, {}, { method: 'DELETE' }); ok(r.status === 405, 'DELETE: 405');
  r = await call('?id=eez&lat=91&lon=0', cache); ok(r.status === 400, 'eez: latitude out of range refused');
  r = await call('?id=eez&lat=1;drop&lon=0', cache); ok(r.status === 400, 'eez: non-numeric coordinate refused');
  r = await call('?id=eez', cache); ok(r.status === 400 && up.length === 0, 'eez: missing coordinates refused, nothing requested upstream');
  r = await call('', cache); j = await r.json(); ok(r.status === 200 && j.feeds.length === Object.keys(FEEDS).length, 'no id: the list of feeds');
  // shared store between locations (optional KV binding)
  const kvm = new Map(), kv = { async get(k, t) { const v = kvm.get(k); return v == null ? null : t === 'json' ? JSON.parse(v) : v; }, async put(k, v) { kvm.set(k, v); } };
  up = []; const a = makeCache(), b = makeCache();
  r = await call('?id=cit', a, { FEED_STORE: kv }); await new Promise((f) => setTimeout(f, 5));
  r = await call('?id=cit', b, { FEED_STORE: kv });
  ok(up.length === 1 && r.status === 200 && r.headers.get('x-feed-cache') === 'store' && (await r.json()).data.rows.NOR[0] === 22, 'with the shared store a second location answers without asking the upstream');
  // twenty simultaneous first requests
  up = []; const c = makeCache(); const rs = await Promise.all(Array.from({ length: 20 }, () => call('?id=cit', c)));
  ok(rs.every((x) => x.status === 200) && up.length === 1, `twenty simultaneous first requests cause one upstream request (${up.length})`);
  // the eez feed builds its address from two rounded numbers only
  up = []; const z = async (url) => { up.push(String(url)); return new Response(JSON.stringify({ features: [{ properties: { mrgid: 5686, geoname: 'Norwegian Exclusive Economic Zone', pol_type: '200NM', territory1: 'Norway', iso_ter1: 'NOR', sovereign1: 'Norway', iso_sov1: 'NOR' } }] })); };
  r = await handle(new Request('https://hydraslug.pages.dev/api/feed?id=eez&lat=61.2049999&lon=2.2000001'), {}, ctx, z, makeCache()); j = await r.json();
  ok(j.data.territory === 'NOR' && j.data.none === false && up.length === 1 && up[0].startsWith('https://geo.vliz.be/') && up[0].includes('POINT(2.20%2061.20)'), 'eez: coordinates are rounded and placed in a fixed address');
  const g = async (url) => { up.push(String(url)); if (String(url).startsWith('https://geo.vliz.be/')) return new Response('down', { status: 503 }); return new Response(JSON.stringify([{ MRGID: 5696, placeType: 'EEZ', preferredGazetteerName: 'British Exclusive Economic Zone' }, { MRGID: 2353, placeType: 'IHO Sea Area', preferredGazetteerName: 'Norwegian Sea' }])); };
  r = await handle(new Request('https://hydraslug.pages.dev/api/feed?id=eez&lat=61.2&lon=1.8'), {}, ctx, g, makeCache()); j = await r.json();
  ok(j.data.territory === 'GBR' && j.data.sovereign === 'GBR' && j.data.mrgid === 5696, 'eez: the gazetteer stands in for the boundary service (zone number → country code)');
  console.log(` ${fails > n0 ? '✗' : '✓'} allow-list, GET only, one upstream request per period, last good answer, size cap, shared store`);
}

console.log('— order of sources: shared cache, the service itself, bundled snapshot (stubbed transports) —');
{
  const recorder = globalThis.fetch, n0 = fails, NO = { countryCode: 'NO', country: 'Norway', data: { iso3: 'NOR', currency: 'NOK' } };
  const citRows = Object.fromEntries(Object.values(ISO3).slice(0, 90).map((c) => [c, [20, 2026]])); citRows.NOR = [22, 2026];
  const answer = (id, data, headers = {}) => new Response(JSON.stringify({ id, fetched: new Date().toISOString(), source: 'stub', data }), { headers: { 'content-type': 'application/json', ...headers } });
  let direct = [], directDown = false, asked = [];
  globalThis.fetch = async (url) => { const u = new URL(String(url)); direct.push(u.hostname); if (directDown) throw new TypeError('fetch failed'); if (u.hostname === 'sdmx.oecd.org') return new Response('DATAFLOW,REF_AREA,TIME_PERIOD,OBS_VALUE\nX,NOR,2026,23\n'); return new Response('not stubbed', { status: 503 }); };
  const edge = (fn) => { asked = []; direct = []; setEdgeTransport(async (url) => { asked.push(String(url)); return fn(new URL(url).searchParams.get('id'), new URL(url).searchParams); }); };
  const down = () => { throw new TypeError('fetch failed'); };
  try {
    edge((id) => answer(id, { rows: citRows })); directDown = false;
    let r = await fetchSource('fiscal', 61, 2, NO);
    ok(r.data.corporateTaxRate === 22 && r.data.corporateTaxVia === 'shared cache' && direct.length === 0 && asked.length === 1 && asked[0] === EDGE + '?id=cit', `step 1: the shared cache answers and the service is not asked (rate ${r.data.corporateTaxRate}, ${r.data.corporateTaxVia}, ${direct.length} direct requests)`);
    edge(down); r = await fetchSource('fiscal', 61, 2, NO);
    ok(r.data.corporateTaxRate === 23 && r.data.corporateTaxVia === 'direct' && asked.length === 1 && direct.join() === 'sdmx.oecd.org', `step 2: shared cache unreachable → the service itself (rate ${r.data.corporateTaxRate}, ${r.data.corporateTaxVia})`);
    edge(() => new Response('<!doctype html><title>HydraSlug</title>', { headers: { 'content-type': 'text/html' } })); r = await fetchSource('fiscal', 61, 2, NO);
    ok(r.data.corporateTaxRate === 23 && r.data.corporateTaxVia === 'direct', 'step 2: a host without the feed (answers with a page) → the service itself');
    edge(() => answer('rates', { rows: citRows })); r = await fetchSource('fiscal', 61, 2, NO);
    ok(r.data.corporateTaxRate === 23, 'step 2: an answer for another feed is not accepted');
    edge((id) => answer(id, { rows: { NOR: [22, 2026] } })); r = await fetchSource('fiscal', 61, 2, NO);
    ok(r.data.corporateTaxRate === 23, 'step 2: an implausibly short table is not accepted');
    edge((id) => answer(id, { rows: citRows }, { 'x-feed-stale': '1', 'x-feed-age': '7200' })); r = await fetchSource('fiscal', 61, 2, NO);
    ok(r.data.corporateTaxRate === 23 && r.data.corporateTaxVia === 'direct', 'a last-good answer of the shared cache yields to a fresh answer of the service');
    edge((id) => answer(id, { rows: citRows }, { 'x-feed-stale': '1', 'x-feed-age': '7200' })); directDown = true; r = await fetchSource('fiscal', 61, 2, NO);
    ok(r.data.corporateTaxRate === 22 && r.data.corporateTaxVia === 'shared cache, last good answer (2 h old)' && direct.includes('sdmx.oecd.org'), `service down too → the last good answer of the shared cache, labelled (${r.data.corporateTaxVia})`);
    edge(down); directDown = true; let threw = null; try { await fetchSource('fiscal', 61, 2, NO); } catch (e) { threw = e; }
    ok(threw && asked.length === 1 && direct.length >= 1, 'both unreachable: the connector fails after asking the shared cache, then the service');
    const s = await fetchSite(60.4, 5.0, () => {}, () => {}, { fresh: true }), d = s.data;
    ok(s.countryCode === 'NO' && d.corporateTaxRate === CORPORATE_TAX.rows.NO[0] && d.atlasFields.includes('corporateTaxRate') && s.status.fiscal.atlas === true && d.corporateTaxVia === undefined, `step 3: both unreachable → the bundled snapshot (${d.corporateTaxRate}, ${s.status.fiscal?.message})`);
    ok(d.costIndex != null && d.atlasFields.includes('costIndex') && d.oilPrice === PRICES.brent && d.electricityPrice === ELECTRICITY.rows.NO[0] && d.taxRate === 78, 'step 3: cost index, oil price, electricity and tax rate from the bundled tables');
    // other feeds through the same order
    const ymNow = new Date().toISOString().slice(0, 7), months = (n, end) => Array.from({ length: n }, (_, i) => { const t = new Date(Date.UTC(+end.slice(0, 4), +end.slice(5, 7) - n + i, 1)); return [t.toISOString().slice(0, 7), 100 + i / 10]; });
    const lastSnap = COSTS.cost[COSTS.cost.length - 1][0], old = months(200, '2024-06'), cur = months(230, lastSnap);
    edge((id) => answer(id, { cost: old, steel: old, cpi: old })); directDown = true; threw = null; try { await fetchSource('costs', 0, 0, {}); } catch (e) { threw = e; }
    ok(threw && /bundled copy is newer/.test(threw.message), `costs: a series older than the bundled snapshot is refused (${threw?.message})`);
    edge((id) => answer(id, { cost: cur, steel: cur, cpi: cur })); r = await fetchSource('costs', 0, 0, {});
    ok(r.data.costIndexVia === 'shared cache' && r.data.costIndexDate === lastSnap && direct.length === 0, 'costs: twenty years of indexes from the shared cache in one request');
    edge((id) => answer(id, { imf: { USA: { policyRate: [3.9, ymNow] } }, oecd: { USA: [4.9, ymNow] }, bis: { NO: [4.25, ymNow] }, us10y: [5.22, ymNow + '-08'] })); directDown = false;
    r = await fetchSource('rates', 0, 0, { countryCode: 'US', data: { iso3: 'USA' } });
    ok(r.data.bondYield === 5.22 && r.data.bondYieldDate === ymNow + '-08' && /Treasury/.test(r.data.bondYieldSource) && r.data.policyRate === 3.9 && r.data.ratesVia === 'shared cache' && direct.length === 0, `rates: United States from the shared cache, daily Treasury yield preferred (${r.data.bondYield}, ${r.data.bondYieldSource})`);
    r = await fetchSource('rates', 0, 0, NO);
    ok(r.data.policyRate === 4.25 && /BIS/.test(r.data.policyRateSource) && r.data.bondYield === undefined, 'rates: a country only the BIS table lists');
    const gb = ELECTRICITY_NATIONAL.rows.GB;
    edge((id) => (id === 'tariffs' ? answer(id, { national: { GB: [gb.price, 'GBP', gb.period] } }) : down()));
    r = await fetchSource('power', 0, 0, { countryCode: 'GB', data: {} });
    ok(Math.abs(r.data.electricityPrice - gb.usd) < 2e-4 && r.data.electricityPriceCurrency === 'GBP' && r.data.electricityPriceKind === 'industrial' && r.data.electricityPriceVia === 'shared cache' && /Energy Security/.test(r.data.electricityPriceSource), `power: United Kingdom national table through the shared cache (${r.data.electricityPrice} $/kWh)`);
    edge((id, q) => (id === 'eez' && q.get('lat') === '61.20' && q.get('lon') === '2.20' ? answer(id, { zones: [{ name: 'Norwegian Exclusive Economic Zone' }], name: 'Norwegian Exclusive Economic Zone', territory: 'NOR', sovereign: 'NOR', territoryName: 'Norway', sovereignName: 'Norway', disputed: false, none: false }) : down()));
    r = await fetchSource('place', 61.2, 2.2, {});
    ok(r.meta.countryCode === 'NO' && r.data.jurisdictionBasis === 'exclusive economic zone' && r.data.jurisdictionZone === 'Norwegian Exclusive Economic Zone' && r.data.jurisdictionVia === 'shared cache' && r.data.currency === 'USD' && direct.join() === 'api.bigdatacloud.net', 'place: offshore point takes the country of its exclusive economic zone from the shared cache (coordinates rounded to 0.01°)');
    edge((id) => answer(id, { zones: [{ name: 'Overlapping claim Falkland / Malvinas Islands: United Kingdom / Argentina' }], name: 'Overlapping claim Falkland / Malvinas Islands: United Kingdom / Argentina', territory: 'FLK', sovereign: 'GBR', territoryName: 'Falkland / Malvinas Islands', sovereignName: 'United Kingdom', disputed: true, none: false }));
    r = await offshoreState(-51, -58.5);
    ok(r.countryCode === 'GB' && /overlapping claim/.test(r.data.jurisdictionNote), 'place: a territory without national data falls to the sovereign state; a disputed zone is flagged');
  } finally { globalThis.fetch = recorder; setEdgeTransport(); }
  console.log(` ${fails === n0 ? '✓' : '✗'} three-step order for tax, costs, rates, electricity and place`);
}

console.log('— live connectors (two sites, shared cache switched off: every request goes to the service itself) —');
edgeOff();
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

console.log('— shared cache: handler in-process, real upstream services —');
const mem = new Map(), memCache = { async match(q) { return mem.get(q.url)?.clone(); }, async put(q, x) { mem.set(q.url, x); } };
let upstreamFailures = 0;
const serverFetch = async (url, opt) => { try { const r = await realFetch(url, opt); if (!r.ok) upstreamFailures++; return r; } catch (e) { upstreamFailures++; throw e; } }; // server side: no browser rules apply
const rowsOk = (rows, n, re) => Array.isArray(rows) && rows.length >= n && rows.every((r) => re.test(r[0]) && Number.isFinite(r[1]));
const FEED_CHECKS = {
  indexes: (d) => ok(rowsOk(d.cost, 200, /^\d{4}-\d\d$/) && rowsOk(d.steel, 200, /^\d{4}-\d\d$/) && rowsOk(d.cpi, 200, /^\d{4}-\d\d$/) && monthsOld(d.cost[d.cost.length - 1][0]) < 9, 'indexes: three series of twenty years, recent'),
  cit: (d) => ok(Object.keys(d.rows).length > 100 && d.rows.NOR?.[0] === 22 && d.rows.NOR[1] >= 2024, `cit: ${Object.keys(d.rows).length} jurisdictions, Norway ${d.rows.NOR}`),
  rates: (d) => { ok(d.imf || d.oecd || d.bis, 'rates: at least one table'); if (d.imf) ok(Object.keys(d.imf).length > 60, 'rates: IMF table'); if (d.oecd) ok(between(d.oecd.NOR?.[0], -1, 20), `rates: OECD Norway ${d.oecd.NOR}`); if (d.bis) ok(between(d.bis.NO?.[0], -1, 30) && between(d.bis.US?.[0], -1, 30), `rates: BIS Norway ${d.bis.NO}`); if (d.us10y) ok(between(d.us10y[0], 0, 20) && monthsOld(d.us10y[1]) < 2, `rates: US 10-year ${d.us10y}`); },
  commodities: (d) => ok(rowsOk(d.series.gasEurope, 100, /^\d{4}-\d\d$/) && rowsOk(d.series.gasAsia, 100, /^\d{4}-\d\d$/) && rowsOk(d.series.ironOre, 100, /^\d{4}-\d\d$/), 'commodities: three ten-year monthly series'),
  spot: (d) => ok(rowsOk(d.brent, 200, /^\d{4}-\d\d-\d\d$/) && rowsOk(d.henryHub, 200, /^\d{4}-\d\d-\d\d$/) && monthsOld(d.brent[d.brent.length - 1][0]) < 2, 'spot: daily Brent and Henry Hub of about a year'),
  power: (d) => { ok(d.eurostat || d.us, 'power: Eurostat or EIA'); if (d.eurostat) ok(Object.keys(d.eurostat).length > 25 && between(d.eurostat.NO?.[0], 0.01, 1) && !d.eurostat.EA, `power: Eurostat ${Object.keys(d.eurostat).length} countries`); },
  tariffs: (d) => { for (const [c, cur] of [['GB', 'GBP'], ['CH', 'CHF'], ['SG', 'SGD']]) { if (d.national?.[c]) ok(between(d.national[c][0], 0.03, 1) && d.national[c][1] === cur && +d.national[c][2] >= 2024 && ELECTRICITY_NATIONAL.rows[c]?.feed === true, `tariffs: national table ${c} ${d.national[c]}`); else console.log(`     (national table of ${c} not reachable: the bundled row stands in)`); } },
  eez: (d) => ok(d.territory === 'NOR' && d.none === false && /Norwegian/.test(d.name), `eez: 61.2 N 2.2 E → ${d.name}`),
};
ok(Object.keys(FEED_CHECKS).sort().join() === Object.keys(FEEDS).sort().join() && Object.keys(EDGE_INFO.feeds).sort().join() === Object.keys(FEEDS).sort().join(), 'every feed of the shared cache is tested and known to the app');
const feedUp = {};
for (const id of Object.keys(FEEDS)) {
  const q = EDGE + '?id=' + id + (id === 'eez' ? '&lat=61.2&lon=2.2' : ''), n = fails;
  const r = await handle(new Request(q), {}, {}, serverFetch, memCache), j = await r.json().catch(() => ({}));
  if (r.status !== 200) { skips++; console.log(` – ${id.padEnd(11)} skipped (upstream unreachable) — ${j.reason || r.status}`); continue; }
  FEED_CHECKS[id](j.data); feedUp[id] = true;
  ok(JSON.stringify(j).length < 300e3 && j.source && j.terms && r.headers.get('x-feed-cache') === 'miss', `${id}: size, source, terms`);
  const again = await handle(new Request(q + '&again=1'), {}, {}, () => { throw new Error('no upstream request expected'); }, memCache);
  ok(again.status === 200 && again.headers.get('x-feed-cache') === 'hit', `${id}: the second request is answered from the cache`);
  console.log(` ${fails === n ? '✓' : '✗'} ${id.padEnd(11)} ${String(JSON.stringify(j).length).padStart(6)} bytes  ${j.source.slice(0, 96)}${j.note ? '  [' + j.note.slice(0, 60) + ']' : ''}`);
}
setEdgeTransport((url) => handle(new Request(url), {}, {}, serverFetch, memCache));
for (const [id, feed, viaKey, site] of [['fiscal', 'cit', 'corporateTaxVia', SITES[1]], ['costs', 'indexes', 'costIndexVia', SITES[1]], ['markets', 'commodities', 'marketsVia', SITES[1]], ['rates', 'rates', 'ratesVia', SITES[1]], ['power', 'power', 'electricityPriceVia', SITES[1]], ['prices', 'spot', 'pricesVia', SITES[1]],
  ['power', 'tariffs', 'electricityPriceVia', { label: 'United Kingdom', lat: 57, lon: 1, site: { countryCode: 'GB', country: 'United Kingdom', data: { iso3: 'GBR' } }, industrialPower: true }]]) {
  const r = await source(id, `${site.label} through the shared cache`, site.lat, site.lon, site.site, (d, m) => CHECKS[id](d, m, site));
  if (r && feedUp[feed]) ok(r.data[viaKey] === 'shared cache', `${id}: answered by the shared cache (${r.data[viaKey]})`);
}

console.log('— offshore jurisdiction: exclusive economic zone, else the nearest country —');
const BASINS = [['North Sea, Norwegian sector 61.2 N 2.2 E', 61.2, 2.2, 'NO'], ['North Sea, UK sector 61.2 N 1.8 E', 61.2, 1.8, 'GB'], ['Gulf of Mexico 28.0 N 90.0 W', 28.0, -90.0, 'US'], ['Santos Basin 25.5 S 43.0 W', -25.5, -43.0, 'BR'], ['Gulf of Guinea 3.5 N 5.75 E', 3.5, 5.75, 'NG'], ['North-West Shelf 19.6 S 116.1 E', -19.6, 116.1, 'AU'], ['Barents Sea 72.5 N 21.0 E', 72.5, 21.0, 'NO']];
for (const [mode, prepare] of [['shared cache', () => setEdgeTransport((url) => handle(new Request(url), {}, {}, serverFetch, memCache))], ['boundary service', edgeOff]]) {
  prepare();
  for (const [name, lat, lon, code] of BASINS) {
    const before = netFailures + upstreamFailures, r = await offshoreState(lat, lon), zone = r?.data?.jurisdictionBasis === 'exclusive economic zone';
    if (!zone && netFailures + upstreamFailures > before) { skips++; console.log(` – ${name}: skipped (${mode} unreachable); nearest country gives ${r?.countryCode}`); continue; }
    ok(zone && r.countryCode === code && r.data.jurisdictionZone && /Marine Regions/.test(r.data.jurisdictionSource), `${name} via ${mode}: ${r?.countryCode} (${r?.data?.jurisdictionBasis}), expected ${code}`);
    if (zone) console.log(` ${r.countryCode === code ? '✓' : '✗'} ${name.padEnd(42)} ${r.countryCode}  ${r.data.jurisdictionZone}  [${r.data.jurisdictionVia}]`);
  }
}
for (const [name, lat, lon, code] of BASINS) { const before = calls, r = await offshoreState(lat, lon, false); ok(calls === before && r?.countryCode === code && /^nearest country, \d+ km$/.test(r.data.jurisdictionBasis) && r.offshoreKm > 50 && r.offshoreKm <= 400 && !r.data.jurisdictionZone, `${name} without a connection: ${r?.countryCode} (${r?.data?.jurisdictionBasis}), expected ${code}, labelled as the nearest country`); }
{ const far = await offshoreState(40, -150, false), near = await nearestCountry(61.2, 1.8, 100); ok(far === null && near === null && (await nearestCountry(61.2, 1.8, 400))?.countryCode === 'GB', 'mid-ocean: no country; the 400 km search reaches where the 100 km one does not'); console.log(' ✓ without a connection: nearest country within 400 km for the seven points, none in mid-ocean'); }

console.log('— the deployed shared cache —');
setEdgeTransport();
try {
  const r = await realFetch(EDGE + '?id=list', { headers: { Origin: ORIGIN } }), j = await r.json().catch(() => null);
  if (!j?.feeds) { skips++; console.log(` – ${EDGE}: not deployed yet (the address answers with the app page) — skipped; the app then uses the services directly`); }
  else {
    const c = await realFetch(EDGE + '?id=cit', { headers: { Origin: ORIGIN } }), cj = await c.json().catch(() => ({}));
    if (c.status !== 200) { skips++; console.log(` – deployed feed cit: HTTP ${c.status} — ${cj.reason || ''}`); }
    else { ok(c.headers.get('access-control-allow-origin') === '*' && cj.data?.rows?.NOR && c.headers.get('x-feed-cache'), 'deployed shared cache: readable cross-origin, rows present'); console.log(` ✓ deployed: ${j.feeds.length} feeds, cit ${c.headers.get('x-feed-cache')}, age ${c.headers.get('x-feed-age')} s`); }
  }
} catch (e) { skips++; console.log(` – ${EDGE}: unreachable — ${e.message}`); }
edgeOff();

console.log('— offline path (built-in atlas and snapshot, no network) —');
const NEW = ['taxRate', 'taxRateYear', 'taxRateSource', 'taxRateBasis', 'corporateTaxRate', 'corporateTaxYear', 'corporateTaxSource', 'costIndex', 'costIndexDate', 'costIndexBase', 'costIndexSeries', 'costEscalation', 'steelIndex', 'steelIndexDate', 'usCpi', 'usCpiDate', 'usCpiSeries', 'usInflation5y',
  'gasPriceEurope', 'gasPriceEuropeDate', 'gasPriceAsia', 'gasPriceAsiaDate', 'gasPriceRegional', 'gasPriceRegionalDate', 'gasPriceRegionalName', 'steelPrice', 'steelPriceDate', 'steelPriceUnit', 'steelPriceWhat', 'policyRate', 'policyRateDate', 'electricityPrice', 'electricityPriceSource'];
for (const [label, lat, lon, code, liveKey] of [['Nigeria offshore 3.5 N 5.75 E', 3.5, 5.75, 'NG', SITES[0].label], ['Norway coast 60.4 N 5.0 E', 60.4, 5.0, 'NO', SITES[1].label], ['North Sea 61.2 N 1.8 E (open sea, UK sector, more than 100 km from any coast)', 61.2, 1.8, 'GB', null], ['mid-Pacific 40 N 150 W (no country within 400 km)', 40, -150, '', null]]) {
  const before = calls, n = fails, s = await fetchSite(lat, lon, () => {}, () => {}, { atlasOnly: true }), d = s.data, AF = new Set(d.atlasFields), notes = d.atlasNotes || {};
  ok(calls === before, `${label}: the offline path made ${calls - before} network requests`);
  ok(s.countryCode === code, `${label}: country ${s.countryCode || '(none)'}, expected ${code || '(none)'}`);
  const global = ['costIndex', 'costEscalation', 'steelIndex', 'usCpi', 'gasPriceEurope', 'gasPriceAsia', 'gasPriceRegional', 'steelPrice', 'oilPrice', 'gasPrice'];
  for (const k of global) ok(d[k] != null && AF.has(k), `${label}: ${k} filled and marked built-in atlas`);
  for (const k of ['costIndex', 'gasPriceEurope', 'steelPrice']) ok(/\d{4}-\d\d/.test(notes[k] || ''), `${label}: ${k} carries its date in the note (${notes[k]})`);
  CHECKS.costs({ ...d, costIndexDate: new Date().toISOString().slice(0, 7) }); ok(/^\d{4}-\d\d$/.test(d.costIndexDate), `${label}: cost index date ${d.costIndexDate}`);
  ok(between(d.gasPriceEurope, 0.5, 150) && between(d.gasPriceAsia, 0.5, 150) && between(d.steelPrice, 20, 400), `${label}: snapshot prices in range`);
  if (code) {
    const f = fiscalOf(code, { waterDepth: d.depth });
    if (code === 'NG') ok(d.depth > 200 && d.taxRate === 30 && d.royaltyRate === 7.5 && /deep offshore/.test(d.fiscalTerrain) && d.fiscalTerrainBasis, `${label}: deep-offshore terms chosen by the water depth (${Math.round(d.depth)} m: tax ${d.taxRate}%, royalty ${d.royaltyRate}%, ${d.fiscalTerrain})`);
    if (lat !== 60.4) ok(/^nearest country, \d+ km$/.test(d.jurisdictionBasis) && /nearest country/.test(notes.country || ''), `${label}: the country is labelled as the nearest one (${d.jurisdictionBasis})`);
    if (UPSTREAM_CI.rows[code]) ok(d.upstreamCarbonIntensity === UPSTREAM_CI.rows[code][0] && d.upstreamCarbonIntensityUnit === 'kgCO₂e/bbl' && /Dixit/.test(d.upstreamCarbonIntensitySource) && d.upstreamCarbonIntensityLow < d.upstreamCarbonIntensityHigh, `${label}: upstream carbon intensity ${d.upstreamCarbonIntensity} with source and range`); else ok(d.upstreamCarbonIntensity === undefined, `${label}: no upstream carbon intensity on file`);
    ok(d.taxRate === f.headline && d.taxRateSource === f.source.citation && +d.taxRateYear === f.year && d.taxRateBasis, `${label}: tax rate ${d.taxRate} from the sourced fiscal table, with source and year`);
    ok(between(d.corporateTaxRate, 0, 60) && AF.has('corporateTaxRate') && +d.corporateTaxYear >= 2024, `${label}: corporate tax ${d.corporateTaxRate} (${d.corporateTaxYear}) from the snapshot`);
    if (RATES.rows[code]) ok(['bondYield', 'policyRate', 'treasuryBillYield'].some((k) => d[k] != null && AF.has(k) && /^\d{4}-\d\d$/.test(d[k + 'Date'])), `${label}: an interest rate from the snapshot, dated`);
    ok(between(d.electricityPrice, 0.01, 0.8) && AF.has('electricityPrice'), `${label}: electricity price ${d.electricityPrice}`);
    if (code === 'NO') ok(/Eurostat/.test(d.electricityPriceSource) && /^\d{4}-S[12]$/.test(d.electricityPriceDate) && !/ESTIMATE/.test(notes.electricityPrice), `${label}: electricity price is the sourced Eurostat value (${d.electricityPriceDate})`);
    if (code === 'NG') ok(/^ESTIMATE/.test(notes.electricityPrice || '') && /estimate/.test(d.electricityPriceSource) && d.electricityPriceKind === 'estimate', `${label}: electricity price is labelled an estimate`);
    if (code === 'GB') ok(/Energy Security/.test(d.electricityPriceSource) && d.electricityPriceKind === 'industrial' && d.electricityPriceCurrency === 'GBP' && d.electricityPriceDate === ELECTRICITY_NATIONAL.rows.GB.period && !/ESTIMATE|PROXY/.test(notes.electricityPrice), `${label}: electricity price is the sourced national value (${d.electricityPrice}, ${d.electricityPriceDate})`);
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

// The shared cache is reached at the Cloudflare mirror: allowed by the page policy, and deliberately NOT in DATA_HOSTS —
// its answers carry `Vary: *`, which the service worker cannot store (the app keeps its own dated copy of every answer),
// and listing the host would also store the mirror's version checks.
ok(csp.includes('https://' + EDGE_INFO.host) && EDGE.startsWith(`https://${EDGE_INFO.host}/`), 'shared cache host is in the connect-src list of index.html');
ok(csp.includes('https://geo.vliz.be') && sw.includes("'geo.vliz.be'"), 'boundary service host is in the connect-src list and in DATA_HOSTS');
console.log(`\n${passes} checks passed, ${skips} sources skipped (unreachable), ${fails} failed.`);
process.exit(fails ? 1 : 0);
