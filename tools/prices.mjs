// Refreshes js/data/prices.js — the dated stand-by snapshot bundled with a build — from the same open series the
// app reads live (the fetch and parse code is shared: js/core/live.js `feeds`):  node tools/prices.mjs
// A source that cannot be reached keeps its previous snapshot, so a partial outage never empties the file.
import { writeFileSync } from 'node:fs';
import { join, dirname } from 'node:path';
import { fileURLToPath } from 'node:url';
import { feeds, getJSON, setEdgeTransport } from '../js/core/live.js';
import { handle } from '../functions/api/feed.js';
import * as OLD from '../js/data/prices.js';

const today = new Date().toISOString().slice(0, 10), log = [];
// The shared cache of the app is run in-process (same code as the deployed function, real upstream services), so the
// snapshot is built along the path the app itself uses: shared feed first, then the service directly.
const mem = new Map(), memCache = { async match(r) { return mem.get(r.url)?.clone(); }, async put(r, x) { mem.set(r.url, x); } };
setEdgeTransport((url) => handle(new Request(url), process.env, {}, fetch, memCache)); // BLS_KEY / EIA_KEY in the environment are honoured
const keep = async (name, old, fn) => { try { const v = await fn(); if (v == null) throw new Error('empty'); log.push(`${name}: refreshed`); return v; } catch (e) { log.push(`${name}: KEPT previous snapshot (${e.message})`); return old; } };
const PRICES = await keep('oil and gas spot prices', OLD.PRICES, async () => {
  const x = await feeds.spot(), end = (rows) => rows[rows.length - 1];
  return x && x.brent.length && x.wti.length && x.henryHub.length ? { date: end(x.brent)[0], brent: end(x.brent)[1], wti: end(x.wti)[1], henryHub: end(x.henryHub)[1] } : null;
});
// ISO-2 -> ISO-3 codes of every country (World Bank country list), needed to address the OECD and IMF services
const ISO3 = await keep('country codes', OLD.ISO3, async () => {
  const j = await getJSON('https://api.worldbank.org/v2/country?format=json&per_page=400');
  const out = Object.fromEntries(j[1].filter((c) => c.region?.id !== 'NA' && /^[A-Z]{2}$/.test(c.iso2Code) && /^[A-Z]{3}$/.test(c.id)).map((c) => [c.iso2Code, c.id]).sort());
  return Object.keys(out).length > 150 ? out : null;
});
const to2 = Object.fromEntries(Object.entries(ISO3).map(([a, b]) => [b, a]));
const byIso2 = (rows) => Object.fromEntries(Object.entries(rows).filter(([k]) => to2[k]).map(([k, v]) => [to2[k], v]).sort());

const MARKETS = await keep('regional gas and iron ore (IMF)', OLD.MARKETS, async () => {
  const m = await feeds.commodities(), rows = {};
  for (const [k, s] of Object.entries(m)) { const r = s[s.length - 1]; rows[k] = [+r[1].toFixed(3), r[0]]; }
  return { retrieved: today, rows };
});
const COSTS = await keep('price indexes (BLS)', OLD.COSTS, async () => { const x = await feeds.indexes(); return x.cost.length > 100 ? { retrieved: today, ...x } : null; });
const CORPORATE_TAX = await keep('corporate tax (OECD)', OLD.CORPORATE_TAX, async () => { const rows = byIso2(await feeds.corporateTax('')); return Object.keys(rows).length > 80 ? { retrieved: today, rows } : null; });
const RATES = await keep('interest rates (IMF)', OLD.RATES, async () => {
  const cut = String(new Date().getUTCFullYear() - 2), rows = byIso2(await feeds.imfRates(''));
  for (const r of Object.values(rows)) for (const k of Object.keys(r)) if (r[k][1] < cut) delete r[k]; // stale series are not bundled
  for (const k of Object.keys(rows)) if (!Object.keys(rows[k]).length) delete rows[k];
  return Object.keys(rows).length > 40 ? { retrieved: today, rows } : null;
});
const power = await feeds.powerTable('').catch(() => null), tariffs = await feeds.tariffs().catch(() => null);
const ELECTRICITY = await keep('industrial electricity (Eurostat, EIA)', OLD.ELECTRICITY, async () => {
  const eur = (await getJSON('https://open.er-api.com/v6/latest/USD')).rates.EUR, eu = power?.eurostat || {}, rows = {};
  if (!(eur > 0.3 && eur < 3)) return null;
  for (const [k, v] of Object.entries(eu).sort()) rows[k] = [+(v[0] / eur).toFixed(4), v[1]];
  if (power?.us) rows.US = power.us; else if (OLD.ELECTRICITY?.rows?.US) { rows.US = OLD.ELECTRICITY.rows.US; log.push('  US industrial price kept (not reachable)'); }
  return Object.keys(rows).length > 20 ? { retrieved: today, eurPerUSD: +eur.toFixed(4), rows } : null;
});
// National rows: those marked feed: true are re-read (price and period); the conversion uses the World Bank annual
// average exchange rate of the period, or today's rate while that year is not yet published. The other rows are
// maintained by hand from the publication cited on each row and are written back unchanged.
const ELECTRICITY_NATIONAL = await keep('industrial electricity (national tables)', OLD.ELECTRICITY_NATIONAL, async () => {
  const rows = structuredClone(OLD.ELECTRICITY_NATIONAL.rows); let n = 0;
  for (const [c, v] of Object.entries(tariffs || {})) {
    const r = rows[c]; if (!r?.feed || r.currency !== v[1]) continue;
    if (r.price === v[0] && r.period === v[2]) { n++; continue; }
    let rate = null, what = '';
    try { const j = await getJSON(`https://api.worldbank.org/v2/country/${c}/indicator/PA.NUS.FCRF?format=json&date=${v[2]}`); rate = j?.[1]?.[0]?.value; what = `${v[1]} per US$, ${v[2]} average (World Bank, official exchange rate PA.NUS.FCRF)`; } catch { /* today's rate */ }
    if (!(rate > 0)) { rate = (await getJSON('https://open.er-api.com/v6/latest/USD')).rates[v[1]]; what = `${v[1]} per US$ on ${today} (open.er-api.com)`; }
    if (!(rate > 0)) continue;
    Object.assign(r, { price: v[0], period: v[2], fxRate: +rate.toFixed(4), fx: what, usd: +(v[0] / rate).toFixed(4) }); n++;
  }
  return n ? { retrieved: today, rows } : null;
});
const UPSTREAM_CI = OLD.UPSTREAM_CI; // from the publication cited in the table; not a live series

const J = (x) => JSON.stringify(x);
const lines = (o) => { const { rows, ...head } = o; return J(head).slice(0, -1) + ',"rows":{\n' + Object.entries(rows).map(([k, v]) => `  ${J(k)}:${J(v)}`).join(',\n') + '\n}}'; };
const pairs = (rows) => '[' + rows.map((r) => `['${r[0]}',${r[1]}]`).join(',') + ']';
const out = `// Dated stand-by snapshot bundled with a build; shown (with its date, and marked "built-in atlas") only when the
// live services cannot be reached. Rewritten by tools/prices.mjs from the same open series the app reads live.

// Last daily spot prices (US EIA series): Brent and WTI in US$/bbl, Henry Hub in US$/MMBtu.
export const PRICES = { date: '${PRICES.date}', brent: ${PRICES.brent}, wti: ${PRICES.wti}, henryHub: ${PRICES.henryHub} };

// IMF Primary Commodity Price System, latest monthly average [value, 'yyyy-mm']: European gas (Netherlands TTF) and
// LNG in Japan in US$/MMBtu, iron ore 62% Fe CFR China in US$/t.
export const MARKETS = ${J(MARKETS)};

// US Bureau of Labor Statistics index levels as published, [['yyyy-mm', value], …]: cost = PCU333132333132 (oil and gas
// field machinery and equipment), steel = WPU1017 (steel mill products), cpi = CUUR0000SA0 (consumer prices, all items).
export const COSTS = {
  retrieved: '${COSTS.retrieved}',
  cost: ${pairs(COSTS.cost)},
  steel: ${pairs(COSTS.steel)},
  cpi: ${pairs(COSTS.cpi)},
};

// OECD Tax Database: combined (central and sub-central) statutory corporate income tax rate, ISO-2 -> [rate %, year].
export const CORPORATE_TAX = ${J(CORPORATE_TAX)};

// IMF Monetary and Financial Statistics, latest month, ISO-2 -> { bondYield, treasuryBillYield, policyRate } as [% a year, 'yyyy-mm'].
export const RATES = ${J(RATES)};

// Electricity price for industry in US$/kWh, ISO-2 -> [price, period]: Eurostat nrg_pc_205 (consumption band 2 000–19 999 MWh
// a year, excluding VAT and other recoverable taxes; converted at eurPerUSD) and, for US, the EIA twelve-month mean
// industrial retail price ending in the month shown.
export const ELECTRICITY = ${J(ELECTRICITY)};

// Electricity price for industry of countries outside the two tables above, ISO-2 -> one dated row with its source.
// kind 'industrial': a published industrial / non-domestic price (price in the national currency per kWh, converted at
// fxRate); rows with feed: true are re-read through the shared cache of the app. kind 'proxy': the utility's average
// revenue per kWh over all customer classes (published in US$) — sourced, but not an industrial tariff, and labelled so.
export const ELECTRICITY_NATIONAL = ${lines(ELECTRICITY_NATIONAL)};

// Upstream carbon intensity of crude oil production by country, ISO-2 -> [value, low, high] in the unit stated.
export const UPSTREAM_CI = ${lines(UPSTREAM_CI)};

// ISO-2 -> ISO-3 country codes (World Bank country list).
export const ISO3 = ${J(ISO3)};
`;
writeFileSync(join(dirname(fileURLToPath(import.meta.url)), '../js/data/prices.js'), out);
console.log(log.join('\n'));
console.log(`prices.js written: ${out.length} characters; ${Object.keys(CORPORATE_TAX.rows).length} tax rows, ${Object.keys(RATES.rows).length} rate rows, ${Object.keys(ELECTRICITY.rows).length} electricity rows, ${Object.keys(ELECTRICITY_NATIONAL.rows).length} national electricity rows, ${COSTS.cost.length} index months`);
