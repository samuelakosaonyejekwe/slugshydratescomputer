// Refreshes js/data/prices.js — the dated stand-by snapshot bundled with a build — from the same open series the
// app reads live (the fetch and parse code is shared: js/core/live.js `feeds`):  node tools/prices.mjs
// A source that cannot be reached keeps its previous snapshot, so a partial outage never empties the file.
import { writeFileSync } from 'node:fs';
import { join, dirname } from 'node:path';
import { fileURLToPath } from 'node:url';
import { feeds, getJSON } from '../js/core/live.js';
import * as OLD from '../js/data/prices.js';

const today = new Date().toISOString().slice(0, 10), log = [];
const keep = async (name, old, fn) => { try { const v = await fn(); if (v == null) throw new Error('empty'); log.push(`${name}: refreshed`); return v; } catch (e) { log.push(`${name}: KEPT previous snapshot (${e.message})`); return old; } };
const last = async (url) => { const rows = (await (await fetch(url)).text()).trim().split('\n').map((l) => l.split(',')).filter((c) => /^\d{4}-\d\d-\d\d$/.test(c[0]) && Number.isFinite(+c[1])); return rows[rows.length - 1]; };

const PRICES = await keep('oil and gas spot prices', OLD.PRICES, async () => {
  const [b, w, g] = await Promise.all(['oil-prices/main/data/brent-daily.csv', 'oil-prices/main/data/wti-daily.csv', 'natural-gas/main/data/daily.csv'].map((p) => last('https://raw.githubusercontent.com/datasets/' + p)));
  return { date: b[0], brent: +b[1], wti: +w[1], henryHub: +g[1] };
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
const ELECTRICITY = await keep('industrial electricity (Eurostat, EIA)', OLD.ELECTRICITY, async () => {
  const eur = (await getJSON('https://open.er-api.com/v6/latest/USD')).rates.EUR, eu = await feeds.eurostatPower(''), rows = {};
  if (!(eur > 0.3 && eur < 3)) return null;
  for (const [k, v] of Object.entries(eu).sort()) rows[k] = [+(v[0] / eur).toFixed(4), v[1]];
  try { const us = await feeds.eiaPower(); rows.US = us; } catch (e) { if (OLD.ELECTRICITY?.rows?.US) rows.US = OLD.ELECTRICITY.rows.US; log.push('  US industrial price kept (' + e.message + ')'); }
  return Object.keys(rows).length > 20 ? { retrieved: today, eurPerUSD: +eur.toFixed(4), rows } : null;
});

const J = (x) => JSON.stringify(x);
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

// ISO-2 -> ISO-3 country codes (World Bank country list).
export const ISO3 = ${J(ISO3)};
`;
writeFileSync(join(dirname(fileURLToPath(import.meta.url)), '../js/data/prices.js'), out);
console.log(log.join('\n'));
console.log(`prices.js written: ${out.length} characters; ${Object.keys(CORPORATE_TAX.rows).length} tax rows, ${Object.keys(RATES.rows).length} rate rows, ${Object.keys(ELECTRICITY.rows).length} electricity rows, ${COSTS.cost.length} index months`);
