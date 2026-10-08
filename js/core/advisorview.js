// Decision-support UI: recommendation cards with live literature evidence, and the
// decision-support & sustainability page.
import { h, clear, btn, kpiGrid, dataTable, toast, badge, fill } from './ui.js';
import { plotCard } from './plot.js';
import { fmt } from './num.js';
import { download } from './io.js';
import { store } from './store.js';
import { advise, sustainability, BENCHMARKS } from './advisor.js';
import { evidence, SOURCES, EVIDENCE_SOURCES } from './live.js';
import { byId, SUITES } from '../suites/index.js';

const PRI = { 1: ['Act now', 'bad'], 2: ['Improve', 'warn'], 3: ['Consider', 'ok'] };

/** One recommendation card. The evidence button queries the open research indexes live. */
export function recCard(r) {
  const ev = h('div', { class: 'evidence' });
  const find = btn('Find current evidence', async () => {
    if (!navigator.onLine) return toast('Evidence look-up needs a connection; the recommendation itself does not.', 'warn');
    find.disabled = true; fill(ev, h('p', { class: 'note' }, 'Searching the open research indexes…'));
    try {
      const res = await evidence(r.topic, 5);
      fill(ev, res.items.length ? h('ul', { class: 'linklist' }, res.items.map((w) => h('li', null, w.url ? h('a', { href: w.url, target: '_blank', rel: 'noopener noreferrer' }, w.title) : w.title, h('small', null, ` · ${[w.venue, w.year, w.cited !== null ? w.cited + ' citations' : ''].filter(Boolean).join(' · ')}`)))) : h('p', { class: 'note' }, 'No recent indexed work matched this topic.'),
        h('p', { class: 'note' }, `Source: ${res.source}, queried ${new Date().toLocaleString()}. Titles are shown as published; read the papers before relying on them.`));
    } catch { fill(ev, h('p', { class: 'note' }, 'The research indexes could not be reached right now.')); }
    find.disabled = false;
  }, 'mini');
  const s = r.goTo && byId(r.goTo);
  return h('article', { class: 'rec p' + r.priority },
    h('header', null, badge(PRI[r.priority][0], PRI[r.priority][1]), h('span', { class: 'rec-area' }, r.area), h('h3', null, r.title)),
    h('dl', null, h('dt', null, 'Finding'), h('dd', null, r.why), h('dt', null, 'Recommended action'), h('dd', null, r.action), h('dt', null, 'Expected benefit'), h('dd', null, r.benefit), h('dt', null, 'Judged against'), h('dd', { class: 'note' }, r.basis)),
    h('div', { class: 'row-tools' }, s ? h('a', { class: 'btn mini primary', href: '#/suite/' + s.id }, `Act in ${s.num}. ${s.short} →`) : null, r.topic ? find : null), ev);
}

/** Panel shown under the results of a suite. */
export function advicePanel(suite, ownRecs = []) {
  const c = store.case, recs = advise(c.outputs, c.site, c.fluid, suite.id);
  if (!recs.length && !ownRecs.length) return h('div', { class: 'recs' }, h('h3', null, 'Decision support'), h('p', { class: 'note' }, 'No benchmark concerns were found for this suite with the results available. Run the linked suites to widen the assessment, or open the Decision support page for the whole case.'), h('a', { class: 'btn mini', href: '#/advisor' }, 'Open decision support →'));
  return h('section', { class: 'advice' }, h('div', { class: 'tbl-head' }, h('h3', null, 'Decision support — what to do next'), h('a', { class: 'btn mini', href: '#/advisor' }, 'Whole-case view →')),
    ownRecs.length ? h('div', { class: 'recs' }, h('b', null, 'From this run'), h('ul', null, ownRecs.map((x) => h('li', null, x)))) : null,
    recs.length ? h('div', { class: 'rec-grid' }, recs.map(recCard)) : null,
    h('p', { class: 'note' }, c.site.fetchedAt ? `Benchmarked with live site context (${c.site.country || 'site'}; data pulled ${new Date(c.site.fetchedAt).toLocaleDateString()}).` : 'Tip: pull Global site data to benchmark against the real seabed temperature, water depth, oil and gas prices, carbon price and financing conditions of your location.'));
}

function reportHTML(recs, S) {
  const doc = document.implementation.createHTMLDocument('Decision report'), add = (p, tag, text) => { const e = doc.createElement(tag); if (text != null) e.textContent = text; p.append(e); return e; };
  add(doc.head, 'style', 'body{font:14px/1.5 system-ui,sans-serif;max-width:900px;margin:24px auto;padding:0 16px;color:#0f172a}h2{border-bottom:2px solid #0d9488;padding-bottom:4px;margin-top:26px}table{border-collapse:collapse;width:100%}td,th{border:1px solid #cbd5e1;padding:5px 8px;text-align:left;font-size:13px}th{background:#f1f5f9}.r{border-left:4px solid #0d9488;padding:6px 12px;margin:10px 0;background:#f8fafc}');
  const b = doc.body, c = store.case;
  add(b, 'h1', `Decision report — ${c.name}`); add(b, 'p', `Site: ${c.site.name || 'not set'}${c.site.country ? ', ' + c.site.country : ''} · generated ${new Date().toLocaleString()}`);
  add(b, 'h2', 'Sustainability scorecard'); add(b, 'p', S.overall === null ? 'Not yet assessed.' : `Overall ${Math.round(S.overall)} / 100 (${S.rating}), ${S.assessed} of ${S.pillars.length} pillars assessed.`);
  const t = add(b, 'table'), hr = add(t, 'tr'); ['Pillar', 'Goal', 'Indicator', 'Score', 'Rating', 'Evidence'].forEach((x) => add(hr, 'th', x));
  for (const p of S.pillars) { const tr = add(t, 'tr'); [p.title, p.sdg, p.value === null ? '—' : `${fmt(p.value, 3)} ${p.unit}`, p.score === null ? '—' : String(Math.round(p.score)), p.rating, p.detail || p.need].forEach((x) => add(tr, 'td', x)); }
  add(b, 'h2', `Recommendations (${recs.length})`);
  for (const r of recs) { const d = add(b, 'div'); d.className = 'r'; add(d, 'h3', `[${PRI[r.priority][0]}] ${r.area}: ${r.title}`); add(d, 'p', 'Finding: ' + r.why); add(d, 'p', 'Action: ' + r.action); add(d, 'p', 'Benefit: ' + r.benefit); add(d, 'p', 'Judged against: ' + r.basis); }
  return '<!doctype html>' + doc.documentElement.outerHTML;
}

/** Whole-case decision support and sustainability page. */
export function advisorPage(root) {
  const c = store.case, recs = advise(c.outputs, c.site, c.fluid), S = sustainability(c.outputs, c.site, c.fluid), d = c.site.data || {};
  const solved = SUITES.filter((s) => c.outputs[s.id]).length;
  const bar = (p) => h('div', { class: 'pillar ' + (p.score === null ? 'na' : p.score >= 80 ? 'ok' : p.score >= 55 ? 'mid' : 'low') },
    h('div', { class: 'pillar-top' }, h('b', null, p.title), h('span', { class: 'chip' }, p.sdg), h('span', { class: 'pillar-score' }, p.score === null ? '—' : Math.round(p.score))),
    h('div', { class: 'meter', role: 'meter', 'aria-valuemin': 0, 'aria-valuemax': 100, 'aria-valuenow': p.score === null ? 0 : Math.round(p.score), 'aria-label': p.title }, (() => { const i = h('i'); i.style.width = (p.score === null ? 0 : Math.max(3, p.score)) + '%'; return i; })()),
    h('p', { class: 'note' }, p.score === null ? p.need + '.' : `${p.value !== null ? fmt(p.value, 3) + ' ' + p.unit + ' — ' : ''}${p.rating}. ${p.detail}`));
  const live = [
    ['Water depth', d.depth > 0 ? d.depth : null, 'm', 'Global relief model (NOAA NCEI)'], ['Seabed temperature', d.seabedTemp, '°C', 'World Ocean Atlas 2023 (bundled)'], ['Sea-surface temperature', d.sst, '°C', 'Open-Meteo marine, now'], ['Mean current', d.currentSpeed, 'm/s', 'Open-Meteo marine'], ['Significant wave height', d.waveHeight, 'm', 'Open-Meteo marine'],
    ['Air temperature', d.airTemp, '°C', 'Open-Meteo weather, now'], ['Ground temperature (0.5 m)', d.depth > 0 ? null : d.groundTemp, '°C', 'Open-Meteo weather'],
    ['Brent crude', d.oilPrice, '$/bbl', d.oilPriceDate ? `US EIA spot series, ${d.oilPriceDate}` : ''], ['Henry Hub gas', d.gasPrice, '$/MMBtu', d.gasPriceDate ? `US EIA spot series, ${d.gasPriceDate}` : ''], ['Carbon price', d.carbonPrice, '$/tCO₂', d.carbonPriceYear ? `Our World in Data, ${d.carbonPriceYear}` : ''],
    ['Grid carbon intensity', d.gridCarbon, 'kgCO₂/kWh', d.gridCarbonLive ? `Our World in Data, ${d.gridCarbonYear || 'latest'}` : d.gridCarbon != null ? 'Bundled planning default' : ''], ['Renewable share of electricity', d.renewableShare, '%', d.renewableShareYear ? `Our World in Data, ${d.renewableShareYear}` : ''],
    ['Electricity price', d.electricityPrice, '$/kWh', d.electricityPriceWB != null && d.electricityPrice === d.electricityPriceWB ? `World Bank, ${d.electricityPriceWBYear}` : 'Bundled planning default (edit in Economics)'],
    ['Inflation', d.inflation, '%/y', d.inflationYear ? `World Bank, ${d.inflationYear}` : ''], ['Lending interest rate', d.lendingRate, '%/y', d.lendingRateYear ? `World Bank, ${d.lendingRateYear}` : ''], [`${d.currency || 'Local currency'} per US dollar`, d.fxPerUSD, '', 'Open exchange-rate service'], ['Petroleum profit tax (indicative)', d.taxRate, '%', 'Bundled planning default'],
  ].filter((r) => r[1] !== null && r[1] !== undefined);
  const plots = [];
  const sc = S.pillars.filter((p) => p.score !== null);
  if (sc.length) plots.push(plotCard({ type: 'bar', title: 'Sustainability scores by pillar (0–100)', ylabel: 'Score', categories: sc.map((p) => p.title), series: [{ name: 'Score', values: sc.map((p) => Math.round(p.score)) }], colors: sc.map((p) => (p.score >= 80 ? '#10b981' : p.score >= 55 ? '#eab308' : '#ef4444')) }, { onDownload: download }));
  if (d.oilPriceSeries?.v?.length > 20) plots.push(plotCard({ type: 'line', title: 'Brent crude — daily spot price over the last year', xlabel: 'Trading days before the latest price', ylabel: '$/bbl', series: [{ name: 'Brent', x: d.oilPriceSeries.v.map((_, i) => i - d.oilPriceSeries.v.length + 1), y: d.oilPriceSeries.v }], hlines: c.outputs.econ?.breakevenPrice ? [{ y: c.outputs.econ.breakevenPrice, label: 'break-even of this case', color: '#ef4444' }] : [] }, { onDownload: download }));
  fill(root, 
    h('header', { class: 'page-head' }, h('h1', null, 'Decision support & sustainability'), h('p', null, 'Every result of this case is benchmarked against published industry practice and the live context of the site — seabed, sea state, prices, carbon and financing — then turned into ranked actions and a sustainability scorecard. Recommendations update each time a suite is run.')),
    kpiGrid([{ label: 'Suites solved', value: `${solved} / ${SUITES.length}` }, { label: 'Recommendations', value: recs.length }, { label: 'Act now', value: recs.filter((r) => r.priority === 1).length, status: recs.some((r) => r.priority === 1) ? 'bad' : 'ok' },
      { label: 'Sustainability score', value: S.overall === null ? '—' : Math.round(S.overall), unit: S.overall === null ? '' : '/ 100', status: S.overall === null ? '' : S.overall >= 70 ? 'ok' : S.overall >= 45 ? 'warn' : 'bad' }, { label: 'Pillars assessed', value: `${S.assessed} / ${S.pillars.length}` }, { label: 'Site context', value: c.site.fetchedAt ? (c.site.country || 'pulled') : 'not pulled', status: c.site.fetchedAt ? 'ok' : 'warn' }]),
    !solved ? h('div', { class: 'linkbar' }, 'Nothing has been solved yet. ', h('a', { class: 'btn mini primary', href: '#/chain' }, 'Run the whole case →')) : null,
    !c.site.fetchedAt ? h('div', { class: 'linkbar quiet' }, 'Site data have not been pulled, so the context of your location (seabed temperature, prices, carbon and financing) is not part of the assessment. ', h('a', { class: 'btn mini', href: '#/site' }, 'Pull global site data →')) : null,
    h('div', { class: 'row-tools' }, btn('Export decision report (HTML)', () => download(reportHTML(recs, S), `${c.name}_decision_report.html`, 'text/html'), 'primary'), btn('Export as JSON', () => download(JSON.stringify({ case: c.name, recommendations: recs, sustainability: { overall: S.overall, pillars: S.pillars } }, null, 1), `${c.name}_decisions.json`, 'application/json'))),
    h('h2', { class: 'sect' }, 'Ranked recommendations'),
    recs.length ? h('div', { class: 'rec-grid' }, recs.map(recCard)) : h('p', { class: 'note' }, solved ? 'No benchmark concerns were found with the suites solved so far.' : 'Run suites to receive recommendations.'),
    h('h2', { class: 'sect' }, 'Sustainability scorecard'), h('div', { class: 'pillars' }, S.pillars.map(bar)),
    plots.length ? h('div', { class: 'plots' }, plots) : null,
    live.length ? dataTable({ title: 'Live context used for these decisions', columns: ['Indicator', 'Value', 'Unit', 'Source'], rows: live }) : null,
    dataTable({ title: 'Benchmarks the results are judged against', columns: ['Indicator', 'Best practice', 'Typical', 'Weak', 'Unit', 'Basis'], rows: Object.values(BENCHMARKS).map((b) => [b.label, b.good, b.typical, b.poor, b.unit, b.basis]), note: 'Benchmarks are planning guides drawn from published industry practice; the applicable design code, company standards, local regulations and contracts always take precedence.' }),
    h('section', { class: 'card' }, h('h2', null, 'Where the knowledge comes from'), h('div', { class: 'sources' }, [...SOURCES, ...EVIDENCE_SOURCES.map((s) => ({ ...s, provider: s.name, name: 'Research evidence' }))].map((s) => h('div', { class: 'source' }, h('b', null, s.name), h('span', null, s.gives), h('small', null, `${s.provider} · ${s.host}`))))));
}
