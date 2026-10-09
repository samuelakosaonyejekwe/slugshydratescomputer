// HydraSlug shared data relay (Cloudflare Pages Function):  GET /api/feed?id=<feed>[&lat=..&lon=..]
//
// Several public statistical services limit the number of requests per network address (the key-less US BLS service
// stops after a handful of requests a day), so a browser on a busy network often gets no answer. This function asks
// each service at most once per period for ALL users of the app and hands out the stored answer:
//   - it is NOT a proxy: `id` selects one of the feeds below, whose upstream addresses are written in this file;
//     no part of a request is ever used as an address (the only inputs are a feed id and, for `eez`, two numbers);
//   - GET only, no cookies, no credentials, `Access-Control-Allow-Origin: *`;
//   - the answer is cut down to the numbers the app needs (compact JSON, size-capped);
//   - fresh copy: Cache API entry kept for the feed's period (`s-maxage`); last good copy: a second entry kept for
//     30 days and served, marked `X-Feed-Stale: 1` with its age in `X-Feed-Age`, whenever the upstream fails;
//   - a failure is remembered for a few minutes so a broken upstream is not asked again by every visitor;
//   - optional bindings (Pages project → Settings → Functions): KV namespace `FEED_STORE` shares the stored copy
//     between Cloudflare locations (the Cache API is per location); secrets `BLS_KEY`, `EIA_KEY` (free registration
//     keys) raise the quotas of those two services, which are counted per network address and would otherwise be
//     shared with every other user of the same Cloudflare addresses. None is required; `BLS_KEY` is recommended.
// The handler is a pure function of its arguments so that tests can drive it with stubs:
//   handle(request, env, ctx, fetchImpl, cacheImpl) -> Response
const H = 3600, DAY = 24 * H;
const MAX_UPSTREAM = 6e6, MAX_OUT = 300e3, STALE_KEEP = 30 * DAY, RETRY_AFTER_FAILURE = 300, RETRY_WHILE_STALE = 900;
const UA = 'HydraSlug-feed/1 (shared cache of public statistics; https://hydraslug.pages.dev)';
const BLS_COST = 'PCU333132333132', BLS_STEEL = 'WPU1017', BLS_CPI = 'CUUR0000SA0';

// ---- small parsers ---------------------------------------------------------------------------------------------------
const year = () => new Date().getUTCFullYear();
const fin = (x) => typeof x === 'number' && Number.isFinite(x);
const csvRows = (csv) => { const [head, ...lines] = csv.trim().split(/\r?\n/).map((l) => l.split(',')); return lines.map((c) => Object.fromEntries(head.map((h, i) => [h.trim(), (c[i] ?? '').trim()]))); };
const xmlAttr = (s, k) => (new RegExp(`\\s${k}="([^"]*)"`).exec(s) || [])[1] || '';
/** IMF SDMX 2.1 structure-specific XML -> [{ country, indicator, obs: [[period, value], …] }]. */
function imfSeries(xml) {
  const out = [];
  for (const m of xml.matchAll(/<Series(\s[^>]*[^/>])>([\s\S]*?)<\/Series>/g)) {
    const obs = [...m[2].matchAll(/<Obs(\s[^>]*)>/g)].map((o) => [xmlAttr(o[1], 'TIME_PERIOD').replace(/^(\d{4})-M(\d\d)$/, '$1-$2'), parseFloat(xmlAttr(o[1], 'OBS_VALUE'))]).filter((o) => /^\d{4}(-\d\d)?$/.test(o[0]) && Number.isFinite(o[1])).sort((a, b) => (a[0] < b[0] ? -1 : 1));
    if (obs.length) out.push({ country: xmlAttr(m[1], 'COUNTRY').slice(0, 8), indicator: xmlAttr(m[1], 'INDICATOR').slice(0, 40), obs });
  }
  return out;
}
const inflate = async (raw) => new Uint8Array(await new Response(new Blob([raw]).stream().pipeThrough(new DecompressionStream('deflate-raw'))).arrayBuffer());
/** Members of a ZIP archive (Uint8Array) whose name passes `want`: { name: text }. Stored and deflated members only. */
export async function unzipText(buf, want) {
  const dv = new DataView(buf.buffer, buf.byteOffset, buf.byteLength), dec = new TextDecoder(), out = {};
  let e = buf.length - 22; while (e >= 0 && dv.getUint32(e, true) !== 0x06054b50) e--;
  if (e < 0) throw new Error('not a ZIP archive');
  let p = dv.getUint32(e + 16, true);
  for (let i = dv.getUint16(e + 10, true); i > 0 && p + 46 <= buf.length && dv.getUint32(p, true) === 0x02014b50; i--) {
    const method = dv.getUint16(p + 10, true), csize = dv.getUint32(p + 20, true), usize = dv.getUint32(p + 24, true), nl = dv.getUint16(p + 28, true), el = dv.getUint16(p + 30, true), cl = dv.getUint16(p + 32, true), off = dv.getUint32(p + 42, true);
    const name = dec.decode(buf.subarray(p + 46, p + 46 + nl));
    if (want(name)) {
      if (usize > 12e6 || (method !== 0 && method !== 8)) throw new Error('archive member not readable: ' + name);
      const start = off + 30 + dv.getUint16(off + 26, true) + dv.getUint16(off + 28, true), raw = buf.subarray(start, start + csize);
      out[name] = dec.decode(method === 0 ? raw : await inflate(raw));
    }
    p += 46 + nl + el + cl;
  }
  return out;
}
const unXml = (s) => s.replace(/<[^>]+>/g, '').replace(/&lt;/g, '<').replace(/&gt;/g, '>').replace(/&quot;/g, '"').replace(/&apos;/g, "'").replace(/&amp;/g, '&');
/**
 * Rows of one worksheet of an XLSX workbook (Uint8Array): the first `head` and the last `tail` rows as { column letter: text }.
 * `sheet` is the visible name of the worksheet.
 */
export async function xlsxRows(buf, sheet, head = 12, tail = 140) {
  const base = await unzipText(buf, (n) => n === 'xl/workbook.xml' || n === 'xl/_rels/workbook.xml.rels' || n === 'xl/sharedStrings.xml');
  const tag = [...(base['xl/workbook.xml'] || '').matchAll(/<sheet\s[^>]*>/g)].map((m) => m[0]).find((t) => unXml(xmlAttr(t, 'name')) === sheet);
  if (!tag) throw new Error(`worksheet “${sheet}” not found`);
  const rel = [...(base['xl/_rels/workbook.xml.rels'] || '').matchAll(/<Relationship\s[^>]*>/g)].map((m) => m[0]).find((t) => xmlAttr(t, 'Id') === xmlAttr(tag, 'r:id'));
  const path = 'xl/' + xmlAttr(rel || '', 'Target').replace(/^\/?xl\//, '').replace(/^\//, '');
  const xml = (await unzipText(buf, (n) => n === path))[path];
  if (!xml) throw new Error('worksheet file missing');
  const shared = (base['xl/sharedStrings.xml'] || '').split('<si>').slice(1).map((s) => unXml(s.slice(0, s.indexOf('</si>')).replace(/<rPh[\s\S]*?<\/rPh>/g, '')));
  const cut = []; let a = 0;
  for (let i = 0; i < head; i++) { const s = xml.indexOf('<row ', a); if (s < 0) break; const t = xml.indexOf('</row>', s); if (t < 0) break; cut.push(xml.slice(s, t)); a = t + 6; }
  const last = []; let b = xml.length;
  for (let i = 0; i < tail; i++) { const s = xml.lastIndexOf('<row ', b - 1); if (s < a) break; const t = xml.indexOf('</row>', s); if (t > 0) last.unshift(xml.slice(s, t)); b = s; }
  return [...cut, ...last].map((r) => {
    const o = { _row: +xmlAttr(r.slice(0, r.indexOf('>') + 1), 'r') };
    for (const m of r.matchAll(/<c\s([^>]*?)(?:\/>|>([\s\S]*?)<\/c>)/g)) {
      const col = (/r="([A-Z]+)\d+"/.exec(m[1]) || [])[1], body = m[2] || '', v = (/<v>([\s\S]*?)<\/v>/.exec(body) || [])[1];
      if (!col) continue;
      if (/t="s"/.test(m[1])) { if (v !== undefined) o[col] = shared[+v] ?? ''; } else if (/t="inlineStr"/.test(m[1])) o[col] = unXml(body); else if (v !== undefined) o[col] = unXml(v);
    }
    return o;
  });
}

// ---- upstream readers (each address is written out here; `get` only adds the time-out and the size cap) -------------
async function blsIndexes(get, env) {
  const y = year(), ids = [BLS_COST, BLS_STEEL, BLS_CPI];
  const ask = async (a, b) => {
    const j = JSON.parse(await get('https://api.bls.gov/publicAPI/v2/timeseries/data/', { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify({ seriesid: ids, startyear: String(a), endyear: String(b), ...(env.BLS_KEY ? { registrationkey: String(env.BLS_KEY) } : {}) }) }));
    if (j?.status !== 'REQUEST_SUCCEEDED') throw new Error('BLS: ' + String(j?.message?.[0] || 'request not processed').slice(0, 80));
    const out = {};
    for (const s of j.Results?.series || []) out[s.seriesID] = (s.data || []).filter((r) => /^M(0[1-9]|1[0-2])$/.test(r.period) && /^\d{4}$/.test(String(r.year)) && Number.isFinite(parseFloat(r.value))).map((r) => [`${r.year}-${r.period.slice(1)}`, parseFloat(r.value)]);
    return out;
  };
  const parts = env.BLS_KEY ? [await ask(y - 19, y)] : await Promise.all([ask(y - 19, y - 10), ask(y - 9, y)]); // the key-less service answers ten years per request
  const join = (id) => parts.flatMap((p) => p[id] || []).sort((a, b) => (a[0] < b[0] ? -1 : 1));
  return { cost: join(BLS_COST), steel: join(BLS_STEEL), cpi: join(BLS_CPI) };
}
/** The same three series from the DBnomics mirror of the BLS files (complete history; refreshed a little later than BLS itself). */
async function mirrorIndexes(get) {
  const lo = `${year() - 19}-01`;
  const one = async (path) => { const d = JSON.parse(await get(`https://api.db.nomics.world/v22/series/BLS/${path}?observations=1&metadata=false`, { ms: 9000, max: 2e6 }))?.series?.docs?.[0] || {}; return (d.period || []).map((t, i) => [String(t).slice(0, 7), parseFloat(d.value?.[i])]).filter((r) => /^\d{4}-\d\d$/.test(r[0]) && r[0] >= lo && Number.isFinite(r[1])); };
  const [cost, steel, cpi] = await Promise.all([one('pc/' + BLS_COST), one('wp/' + BLS_STEEL), one('cu/' + BLS_CPI)]);
  const end = cost[cost.length - 1]?.[0] || '', limit = new Date(Date.now() - 200 * DAY * 1000).toISOString().slice(0, 7);
  if (end < limit) throw new Error(`mirror ends ${end || 'nowhere'}: too old to stand in`); // an old mirror must not displace a newer stored or bundled copy
  return { cost, steel, cpi };
}
async function oecdCsv(get, flow, key) { return csvRows(await get(`https://sdmx.oecd.org/public/rest/data/${flow}/${key}?lastNObservations=1&format=csvfile`)); }
// National industrial electricity prices outside the Eurostat and EIA tables: [price in the national currency per kWh, ISO 4217 code, period].
const NATIONAL_POWER = {
  /** United Kingdom: DESNZ Quarterly Energy Prices table 5.4.2, medium consumers (2 000–19 999 MWh a year), including taxes and levies, excluding VAT. */
  async GB(get) {
    const j = JSON.parse(await get('https://www.gov.uk/api/content/government/statistical-data-sets/international-non-domestic-energy-prices', { max: 5e5 }));
    const url = (j?.details?.attachments || []).map((a) => String(a?.url || '')).find((u) => /^https:\/\/assets\.publishing\.service\.gov\.uk\/media\/[0-9a-f]{16,40}\/table_541\.xlsx$/.test(u)); // the file moves with every release; only this exact shape of address is followed
    if (!url) throw new Error('DESNZ: table 5.4.1 not listed');
    const rows = await xlsxRows(await get(url, { binary: true, max: 3e6, ms: 9000 }), 'Annual incl tax', 24, 110), head = rows.find((r) => r.A === 'Customer size'), col = head && Object.keys(head).find((c) => head[c] === 'UK');
    const best = rows.filter((r) => r.A === 'Medium' && /^\d{4}/.test(r.B || '') && parseFloat(r[col]) > 0).sort((a, b) => parseInt(a.B, 10) - parseInt(b.B, 10)).pop();
    if (!col || !best || !(parseFloat(best[col]) < 200)) throw new Error('DESNZ: UK column not found');
    return [+(parseFloat(best[col]) / 100).toFixed(4), 'GBP', String(parseInt(best.B, 10))];
  },
  /** Switzerland: ElCom median regulated tariff of the current tariff year, consumption profile C7 (7.5 GWh a year), standard product, excluding VAT. */
  async CH(get) {
    const q = 'PREFIX cube: <https://cube.link/> PREFIX d: <https://energy.ld.admin.ch/elcom/electricityprice/dimension/> SELECT ?period ?total WHERE { <https://energy.ld.admin.ch/elcom/electricityprice-swiss> cube:observationSet/cube:observation ?o . ?o d:period ?period ; d:category ?category ; d:product ?product ; d:total ?total . FILTER(STR(?period) >= "2024" && REGEX(STR(?category),"/C7$") && REGEX(STR(?product),"/standard$")) } ORDER BY ?period';
    const csv = await get('https://lindas.admin.ch/query', { method: 'POST', headers: { Accept: 'text/csv', 'Content-Type': 'application/x-www-form-urlencoded' }, body: 'query=' + encodeURIComponent(q), max: 2e5 }), y = year();
    const rows = csv.trim().split(/\r?\n/).slice(1).map((l) => l.split(',')).map((c) => [parseInt(c[0], 10), parseFloat(c[1])]).filter((r) => r[0] <= y && r[1] > 0 && r[1] < 200).sort((a, b) => a[0] - b[0]), r = rows.pop(); // next year's tariffs are published in advance: not used
    if (!r) throw new Error('ElCom: no tariff');
    return [+(r[1] / 100).toFixed(5), 'CHF', String(r[0])];
  },
  /** Singapore: SP Group regulated tariff, high-tension large supplies, energy charge weighted over the published peak (07:00–23:00) and off-peak (23:00–07:00) periods; capacity charges and GST excluded. */
  async SG(get) {
    const rec = JSON.parse(await get('https://data.gov.sg/api/action/datastore_search?resource_id=d_d610f8ed1864daa6c7e790318bbc3323&limit=50', { max: 2e5 }))?.result?.records || [];
    const row = (re) => rec.find((r) => re.test(String(r?.DataSeries || ''))), pk = row(/^High Tension Large Supplies - Peak Period \(7\.00 AM To 11\.00 PM\)/), op = row(/^High Tension Large Supplies - Off-Peak Period \(11\.00 PM To 7\.00 AM\)/);
    const yr = pk && Object.keys(pk).filter((k) => /^\d{4}$/.test(k) && parseFloat(pk[k]) > 0 && parseFloat(op?.[k]) > 0).sort().pop();
    if (!yr) throw new Error('data.gov.sg: tariff rows not found');
    return [+((16 * parseFloat(pk[yr]) + 8 * parseFloat(op[yr])) / 24 / 100).toFixed(5), 'SGD', yr];
  },
};

// Exclusive economic zones of the Maritime Boundaries Geodatabase (version 12): gazetteer record number -> ISO-3 code of
// the territory (and of the sovereign state where that differs). Needed only when the boundary service answers by name.
const EEZ_ISO = '3293:BEL 5668:NLD 5669:DEU 5670:ALB 5672:BGR 5673:HRV 5674:DNK 5675:EST 5676:FIN 5677:FRA 5678:GEO 5679:GRC 5680:ISL 5681:IRL 5682:ITA 5683:LVA 5684:LTU 5685:MLT 5686:NOR 5687:POL 5688:PRT 5689:ROU 5690:RUS 5691:MNE 5692:SVN 5693:ESP 5694:SWE 5695:UKR 5696:GBR 5697:TUR 8308:CCK/AUS 8309:CXR/AUS 8310:NFK/AUS 8311:/AUS 8312:NCL/FRA 8313:VUT 8314:SLB 8315:PLW 8316:FSM 8317:NRU 8318:MHL 8319:UMI/USA 8321:TWN 8322:PHL 8323:AUS 8324:PNG 8325:FJI 8326:TUV 8327:KOR 8328:PRK 8331:KHM 8332:THA 8333:/IND 8334:COM 8337:SYC 8338:REU/FRA 8339:ATF/FRA 8340:ATF/FRA 8341:ATF/FRA 8343:MUS 8345:MDV 8346:LKA 8347:MOZ 8348:MDG 8349:KEN 8350:SOM 8351:ERI 8352:DJI 8353:YEM 8354:OMN 8355:SDN 8356:SAU 8357:KWT 8358:BHR 8359:PAK 8360:ARE 8361:/PRT 8362:CPV 8363:/PRT 8364:/ESP 8365:GIB 8366:TUN 8367:MAR 8368:ESH 8369:MRT 8370:GMB 8371:SEN 8372:LBY 8373:SYR 8374:LBN 8375:ISR 8376:CYP 8378:DZA 8379:SHN/GBR 8380:SHN/GBR 8381:/BRA 8382:SHN/GBR 8383:SGS/GBR 8384:/ZAF 8385:ATF/FRA 8386:ATF/FRA 8387:ATF/FRA 8388:HMD/AUS 8389:FLK/GBR 8390:SLE 8391:LBR 8392:TGO 8393:BEN 8394:COG 8395:NAM 8396:ZAF 8397:STP 8398:GNQ 8400:GHA 8401:/FRA 8402:BMU/GBR 8403:/ECU 8404:BHS 8405:TCA/GBR 8406:CUB 8407:CYM/GBR 8408:HTI 8409:DOM 8411:VGB/GBR 8412:AIA/GBR 8413:KNA 8414:ATG 8415:MSR/GBR 8416:LCA 8417:DMA 8418:BRB 8419:GRD 8420:TTO 8421:VCT 8423:PAN 8424:CRI 8425:NIC 8426:COL 8427:HND 8428:SLV 8429:MEX 8430:GTM 8431:ECU 8432:PER 8433:VEN 8435:FRO/DNK 8437:SJM/NOR 8438:GRL/DNK 8439:PCN/GBR 8440:PYF/FRA 8441:KIR 8442:UMI/USA 8443:UMI/USA 8444:ASM/USA 8445:WSM 8446:COK/NZL 8447:NIU/NZL 8448:TON 8449:TKL/NZL 8450:KIR 8451:UMI/USA 8452:UMI/USA 8453:/USA 8454:WLF/FRA 8455:NZL 8456:USA 8457:BLZ 8459:JAM 8460:GUY 8461:SUR 8462:GUF/FRA 8463:/USA 8464:BRA 8465:CHL 8466:ARG 8467:URY 8468:QAT 8469:IRN 8470:IRQ 8471:GNB 8472:GIN 8473:CIV 8474:NGA 8475:CMR 8476:GAB 8477:COD 8478:AGO 8479:TZA 8480:IND 8481:BGD 8482:MMR 8483:MYS 8484:VNM 8485:SGP 8486:CHN 8487:JPN 8488:KIR 8490:EGY 8491:JOR 8492:IDN 8493:CAN 8494:SPM/FRA 8495:MAF/FRA 8758:TLS 21787:/CHL 21788:GGY/GBR 21789:JEY/GBR 21790:MCO 21791:/TLS 21792:JAM 21796:KOR 21797:STP 21798:PNG 21803:SXM/NLD 22491:BIH 22756:CHL 26517:CUW/NLD 26518:BES/NLD 26519:ABW/NLD 26520:BES/NLD 26521:BRN 26522:KAZ 26523:TKM 26524:AZE 26526:BES/NLD 26582:/SDN 33177:GLP/FRA 33178:MTQ/FRA 33179:PRI/USA 33180:VIR/USA 33181:SJM/NOR 33185:GUY 48943:CAN 48944:MYT/FRA 48945:/MDG 48946:/FRA 48947:/ARE 48948:/FRA 48950:/JPN 48951:/HTI 48952:BLM/FRA 48953:PSE 48954:/TWN 48955:/JPN 48956:/ERI 48957:GUM/USA 48961:URY 48962:PER 48964:SEN 48965:VEN 48966:FRA 48967:GBR 48968:CRI 48969:BRB 48970:DOM 48971:COL 48972:HND 48973:FRO/DNK 48974:ECU 48975:ISL 48976:FRA 48977:NOR 48978:USA 48980:MNP/USA 48982:PRI/USA 48984:/COL 48985:/COL 48997:/ESP 48998:/ESP 48999:/ESP 49000:/ESP 49001:/ESP 49002:/ESP 49003:CHN 50167:HRV 50170:QAT 62589:/MUS 62596:/COL 62598:/COL 64430:EGY 64431:BLZ 64440:MYS 64446:IMN/GBR 64459:NOR 64460:IRQ';
const eezTable = () => Object.fromEntries(EEZ_ISO.split(' ').map((s) => { const [id, c] = s.split(':'), [t, v] = (c || '').split('/'); return [id, { territory: t || v || '', sovereign: v || t || '' }]; }));

const TERMS = {
  bls: 'US Bureau of Labor Statistics — public domain; BLS.gov cannot vouch for the data or analyses derived from these data after the data have been retrieved from BLS.gov',
  treasury: 'US Department of the Treasury, daily par yield curve rates — public domain',
  oecd: 'OECD — data may be shared and adapted for any purpose with credit to the OECD', imf: 'International Monetary Fund — data may be distributed with attribution to the IMF as the source (commercial reuse needs permission)', bis: 'Bank for International Settlements — use unrestricted with the BIS cited as the source',
  eurostat: 'Eurostat — reuse authorised with the source acknowledged', eia: 'US Energy Information Administration — public domain, acknowledgement requested',
  vliz: 'Flanders Marine Institute (2023), Maritime Boundaries Geodatabase: Maritime Boundaries and Exclusive Economic Zones (200NM), version 12, https://doi.org/10.14284/632 — Creative Commons Attribution 4.0',
  pddl: 'Open “datasets” collection (US EIA series) — Open Data Commons Public Domain Dedication and License',
};

/** The feeds. `load(get, env, params)` returns { data, source, terms } or throws; `params(url)` validates the query. */
export const FEEDS = {
  indexes: {
    ttl: 12 * H, what: 'US price indexes, 20 years monthly: oil and gas field machinery (PCU333132333132), steel mill products (WPU1017), consumer prices (CUUR0000SA0)',
    async load(get, env) {
      let data, source = 'US Bureau of Labor Statistics', terms = TERMS.bls, note;
      try { data = await blsIndexes(get, env); if (data.cost.length < 200) throw new Error('BLS: short answer'); }
      catch (e) { try { data = await mirrorIndexes(get); } catch (m) { throw new Error(`${String(e.message || e).slice(0, 70)}; ${String(m.message || m).slice(0, 70)}`); } source = 'US Bureau of Labor Statistics series, read from the DBnomics mirror'; terms = TERMS.bls + '; DBnomics mirror: data under the terms of the original provider'; note = String(e.message || e).slice(0, 120); }
      if (data.cost.length < 200 || data.steel.length < 200 || data.cpi.length < 200) throw new Error('index series too short');
      return { data, source, terms, ...(note ? { note } : {}) };
    },
  },
  cit: {
    ttl: 7 * DAY, what: 'OECD Tax Database: combined statutory corporate income tax rate of every jurisdiction, latest year',
    async load(get) {
      const rows = {};
      for (const r of await oecdCsv(get, 'OECD.CTP.TPS,DSD_TAX_CIT@DF_CIT,', '.A.CIT_C.ST.....')) { const v = parseFloat(r.OBS_VALUE), y = parseInt(r.TIME_PERIOD, 10); if (/^[A-Z]{3}$/.test(r.REF_AREA) && v >= 0 && v <= 90 && y > 1990 && !(rows[r.REF_AREA]?.[1] > y)) rows[r.REF_AREA] = [v, y]; }
      if (Object.keys(rows).length < 60) throw new Error('OECD: short answer');
      return { data: { rows }, source: 'OECD Tax Database, combined statutory corporate income tax rate', terms: TERMS.oecd };
    },
  },
  rates: {
    ttl: 24 * H, what: 'Interest rates of every country, latest month: IMF government bond yield, Treasury bill yield and policy rate; OECD 10-year yield; BIS policy rate; US 10-year Treasury par yield (daily)',
    async load(get) {
      const IMF = { S13BOND_RT_PT_A_PT: 'bondYield', GSTBILY_RT_PT_A_PT: 'treasuryBillYield', MFS166_RT_PT_A_PT: 'policyRate' };
      const [a, b, c, d] = await Promise.allSettled([
        (async () => { const out = {}; for (const s of imfSeries(await get(`https://api.imf.org/external/sdmx/2.1/data/MFS_IR/.${Object.keys(IMF).join('+')}.M?lastNObservations=1`, { ms: 9000 }))) { const k = IMF[s.indicator], o = s.obs[s.obs.length - 1]; if (k && /^[A-Z]{3}$/.test(s.country) && Math.abs(o[1]) < 200) (out[s.country] ??= {})[k] = [+o[1].toFixed(3), o[0]]; } if (Object.keys(out).length < 30) throw new Error('IMF: short answer'); return out; })(),
        (async () => { const out = {}; for (const r of await oecdCsv(get, 'OECD.SDD.STES,DSD_STES@DF_FINMARK,4.0', '.M.IRLT.PA.....')) { const v = parseFloat(r.OBS_VALUE); if (/^[A-Z]{3}$/.test(r.REF_AREA) && /^\d{4}-\d\d$/.test(r.TIME_PERIOD) && Math.abs(v) < 200 && !(out[r.REF_AREA]?.[1] > r.TIME_PERIOD)) out[r.REF_AREA] = [+v.toFixed(3), r.TIME_PERIOD]; } if (Object.keys(out).length < 10) throw new Error('OECD: short answer'); return out; })(),
        (async () => { const out = {}; for (const r of csvRows((await get('https://stats.bis.org/api/v2/data/dataflow/BIS/WS_CBPOL/1.0/M.?lastNObservations=1&format=csv')).replace(/"[^"]*"/g, ''))) { const v = parseFloat(r.OBS_VALUE); if (/^[A-Z0-9]{2}$/.test(r.REF_AREA) && /^\d{4}-\d\d$/.test(r.TIME_PERIOD) && Math.abs(v) < 200) out[r.REF_AREA] = [v, r.TIME_PERIOD]; } if (Object.keys(out).length < 10) throw new Error('BIS: short answer'); return out; })(),
        (async () => { // the Treasury publishes one small file per month; early in a month the previous one holds the latest day
          for (const back of [0, 1]) { const d = new Date(); d.setUTCDate(1); d.setUTCMonth(d.getUTCMonth() - back); const xml = await get(`https://home.treasury.gov/resource-center/data-chart-center/interest-rates/pages/xml?data=daily_treasury_yield_curve&field_tdr_date_value_month=${d.toISOString().slice(0, 7).replace('-', '')}`, { max: 4e5 }); const rows = [...xml.matchAll(/<d:NEW_DATE[^>]*>(\d{4}-\d\d-\d\d)[^<]*<[\s\S]*?<d:BC_10YEAR[^>]*>(-?[\d.]+)</g)].map((m) => [+m[2], m[1]]).filter((r) => Math.abs(r[0]) < 50).sort((a, b) => (a[1] < b[1] ? -1 : 1)); if (rows.length) return rows[rows.length - 1]; }
          throw new Error('Treasury: no yield');
        })(),
      ]);
      const ok = (x) => (x.status === 'fulfilled' ? x.value : null), data = { imf: ok(a), oecd: ok(b), bis: ok(c), us10y: ok(d) };
      if (!data.imf && !data.oecd && !data.bis) throw new Error('no interest-rate service answered');
      return { data, source: 'IMF Monetary and Financial Statistics; OECD long-term interest rates; BIS central-bank policy rates; US Department of the Treasury daily par yield curve', terms: [TERMS.imf, TERMS.oecd, TERMS.bis, TERMS.treasury].join('; ') };
    },
  },
  commodities: {
    ttl: 24 * H, what: 'Monthly benchmark prices, ten years: European gas (TTF), LNG in Japan, iron ore',
    async load(get) {
      const IDS = { PNGASEU: 'gasEurope', PNGASJP: 'gasAsia', PIORECR: 'ironOre' }, series = {};
      for (const s of imfSeries(await get(`https://api.imf.org/external/sdmx/2.1/data/PCPS/G001.${Object.keys(IDS).join('+')}.USD.M?startPeriod=${year() - 10}-M01`, { ms: 9000 }))) if (IDS[s.indicator]) series[IDS[s.indicator]] = s.obs.map((o) => [o[0], +o[1].toFixed(3)]);
      if (!series.gasEurope?.length) throw new Error('IMF: no commodity series');
      return { data: { provider: 'imf', series }, source: 'IMF Primary Commodity Price System', terms: TERMS.imf };
    },
  },
  spot: {
    ttl: 6 * H, what: 'Daily spot prices, about one year: Brent and WTI crude, Henry Hub natural gas (US EIA series)',
    async load(get, env) {
      const keep = (rows) => rows.filter((r) => r[1] > 0).slice(-270), full = (d) => d.brent?.length > 30 && d.henryHub?.length > 30;
      const tries = [
        ['US EIA daily spot series from the open “datasets” collection', TERMS.pddl, async () => { const one = async (p) => keep((await get('https://raw.githubusercontent.com/datasets/' + p)).trim().split('\n').slice(-300).map((l) => l.split(',')).filter((c) => /^\d{4}-\d\d-\d\d$/.test(c[0]) && Number.isFinite(parseFloat(c[1]))).map((c) => [c[0], parseFloat(c[1])])); const [b, w, g] = await Promise.all(['oil-prices/main/data/brent-daily.csv', 'oil-prices/main/data/wti-daily.csv', 'natural-gas/main/data/daily.csv'].map(one)); return { brent: b, wti: w, henryHub: g }; }],
        ['US EIA open-data service', TERMS.eia, async () => { const one = async (path, id) => keep((JSON.parse(await get(`https://api.eia.gov/v2/${path}/data/?api_key=${encodeURIComponent(env.EIA_KEY || 'DEMO_KEY')}&frequency=daily&data[0]=value&facets[series][]=${id}&sort[0][column]=period&sort[0][direction]=desc&length=270`))?.response?.data || []).map((r) => [String(r.period).slice(0, 10), parseFloat(r.value)]).filter((r) => Number.isFinite(r[1])).reverse()); const [b, w, g] = await Promise.all([one('petroleum/pri/spt', 'RBRTE'), one('petroleum/pri/spt', 'RWTC'), one('natural-gas/pri/fut', 'RNGWHHD')]); return { brent: b, wti: w, henryHub: g }; }],
      ];
      let err;
      for (const [source, terms, fn] of tries) { try { const data = await fn(); if (full(data)) return { data, source, terms }; err = new Error('short answer'); } catch (e) { err = e; } }
      throw err;
    },
  },
  power: {
    ttl: 7 * DAY, what: 'Electricity price for industry: Eurostat (Europe, EUR/kWh, band 2 000–19 999 MWh a year, excluding VAT and recoverable taxes), US EIA (US$/kWh, twelve-month mean)',
    async load(get, env) {
      const [eu, us] = await Promise.allSettled([
        (async () => {
          const j = JSON.parse(await get('https://ec.europa.eu/eurostat/api/dissemination/statistics/1.0/data/nrg_pc_205?format=JSON&lang=EN&nrg_cons=MWH2000-19999&tax=X_VAT&currency=EUR&lastTimePeriod=4', { ms: 9000 }));
          const gi = j?.dimension?.geo?.category?.index || {}, ti = j?.dimension?.time?.category?.index || {}, nt = Object.keys(ti).length, out = {};
          if (!nt || (j.size || []).reduce((a, b) => a * b, 1) !== Object.keys(gi).length * nt) throw new Error('Eurostat: unexpected table layout');
          for (const [g, i] of Object.entries(gi)) for (const [t, k] of Object.entries(ti)) { const v = j.value?.[String(i * nt + k)], c = g === 'EL' ? 'GR' : g; if (/^[A-Z]{2}$/.test(c) && c !== 'EA' && fin(v) && v > 0 && v < 2 && !(out[c]?.[1] > t)) out[c] = [v, String(t).slice(0, 8)]; }
          if (Object.keys(out).length < 15) throw new Error('Eurostat: short answer');
          return out;
        })(),
        (async () => {
          const j = JSON.parse(await get(`https://api.eia.gov/v2/electricity/retail-sales/data/?api_key=${encodeURIComponent(env.EIA_KEY || 'DEMO_KEY')}&frequency=monthly&data[0]=price&facets[sectorid][]=IND&facets[stateid][]=US&sort[0][column]=period&sort[0][direction]=desc&length=12`));
          const rows = (j?.response?.data || []).map((r) => [String(r.period).slice(0, 7), parseFloat(r.price)]).filter((r) => /^\d{4}-\d\d$/.test(r[0]) && r[1] > 0 && r[1] < 100);
          if (rows.length < 6) throw new Error('EIA: no industrial price');
          return [+(rows.reduce((a, r) => a + r[1], 0) / rows.length / 100).toFixed(4), rows[0][0]];
        })(),
      ]);
      const data = { eurostat: eu.status === 'fulfilled' ? eu.value : null, us: us.status === 'fulfilled' ? us.value : null };
      if (!data.eurostat && !data.us) throw new Error('no electricity-price service answered');
      return { data, source: 'Eurostat, electricity prices for non-household consumers (nrg_pc_205); US EIA, average retail price of electricity, industrial sector', terms: [TERMS.eurostat, TERMS.eia].join('; ') };
    },
  },
  tariffs: { // kept apart from `power`: reading the UK workbook is the heaviest job here and must not hold up the two big tables
    ttl: 7 * DAY, what: 'Electricity price for industry from national tables, national currency per kWh: United Kingdom (DESNZ, medium non-domestic consumers), Switzerland (ElCom median tariff, 7.5 GWh a year), Singapore (SP Group high-tension large supplies)',
    async load(get) {
      const national = {};
      await Promise.all(Object.entries(NATIONAL_POWER).map(async ([c, fn]) => { try { national[c] = await fn(get); } catch { /* this country then comes from the table bundled with the app */ } }));
      if (!Object.keys(national).length) throw new Error('no national tariff table answered');
      return { data: { national }, source: 'UK DESNZ Quarterly Energy Prices table 5.4.2; Swiss Federal Electricity Commission (ElCom) median tariffs; SP Group tariffs published by the Singapore Department of Statistics', terms: 'DESNZ: Open Government Licence v3.0; ElCom on opendata.swiss: open use; Singapore Open Data Licence v1.0' };
    },
  },
  eez: {
    ttl: 90 * DAY, what: 'Exclusive economic zone at a point (lat, lon in degrees): zone name and the ISO-3 code of its territory and sovereign state',
    params(q) {
      const lat = Number(q.get('lat')), lon = Number(q.get('lon'));
      if (q.get('lat') === null || q.get('lon') === null || !(Math.abs(lat) <= 90) || !(Math.abs(lon) <= 180)) return null;
      return { lat: (Math.round(lat * 100) / 100).toFixed(2), lon: (Math.round(lon * 100) / 100).toFixed(2) }; // 0.01° ≈ 1 km: one stored answer per cell
    },
    async load(get, env, p) {
      const pick = (zones) => ({ zones, ...(zones.find((z) => z.type === '200NM') || zones[0] || {}), disputed: zones.length > 1 || (zones[0] && zones[0].type !== '200NM') || false, none: !zones.length });
      try {
        const j = JSON.parse(await get(`https://geo.vliz.be/geoserver/MarineRegions/wfs?service=WFS&version=1.0.0&request=GetFeature&typeName=MarineRegions:eez&cql_filter=INTERSECTS(the_geom,POINT(${p.lon}%20${p.lat}))&outputFormat=application/json&propertyName=mrgid,geoname,pol_type,territory1,iso_ter1,sovereign1,iso_sov1,iso_sov2&maxFeatures=6`, { ms: 8000 }));
        if (!Array.isArray(j?.features)) throw new Error('boundary service: unexpected answer');
        const zones = j.features.map((f) => f.properties || {}).filter((z) => Number.isInteger(z.mrgid)).map((z) => ({ mrgid: z.mrgid, name: String(z.geoname || '').slice(0, 120), type: String(z.pol_type || '').slice(0, 30), territory: String(z.iso_ter1 || z.iso_sov1 || '').slice(0, 3), sovereign: String(z.iso_sov1 || z.iso_ter1 || '').slice(0, 3), territoryName: String(z.territory1 || '').slice(0, 60), sovereignName: String(z.sovereign1 || '').slice(0, 60), ...(z.iso_sov2 ? { otherClaimant: String(z.iso_sov2).slice(0, 3) } : {}) }));
        return { data: pick(zones), source: 'Marine Regions, Maritime Boundaries Geodatabase version 12 (Flanders Marine Institute)', terms: TERMS.vliz };
      } catch (e) { // second service of the same institute: the gazetteer answers by name, the table above supplies the codes
        const j = JSON.parse(await get(`https://www.marineregions.org/rest/getGazetteerRecordsByLatLong.json/${p.lat}/${p.lon}/`, { ms: 8000 })), T = eezTable();
        if (!Array.isArray(j)) throw e;
        const zones = j.filter((r) => r.placeType === 'EEZ' && T[r.MRGID]).map((r) => ({ mrgid: r.MRGID, name: String(r.preferredGazetteerName || '').slice(0, 120), type: /^Overlapping/.test(r.preferredGazetteerName) ? 'Overlapping claim' : /^Joint/.test(r.preferredGazetteerName) ? 'Joint regime' : '200NM', ...T[r.MRGID], territoryName: '', sovereignName: '' }));
        return { data: pick(zones), source: 'Marine Regions gazetteer, Maritime Boundaries Geodatabase version 12 (Flanders Marine Institute)', terms: TERMS.vliz };
      }
    },
  },
};

// ---- the handler -----------------------------------------------------------------------------------------------------
const CORS = { 'Access-Control-Allow-Origin': '*', 'Access-Control-Expose-Headers': 'Age, X-Feed-Cache, X-Feed-Age, X-Feed-Stale, X-Feed-Fetched', 'X-Content-Type-Options': 'nosniff', 'Cross-Origin-Resource-Policy': 'cross-origin', 'Referrer-Policy': 'no-referrer' };
const json = (status, obj, headers = {}) => new Response(typeof obj === 'string' ? obj : JSON.stringify(obj), { status, headers: { 'Content-Type': 'application/json; charset=utf-8', ...CORS, ...headers } });
const inflight = new Map(); // one upstream request per feed at a time within this instance

/** Bounded upstream GET/POST: time-out (short enough for the whole feed to answer before the app gives up on it), size cap, text or bytes. Only called with addresses written in this file. */
function reader(fetchImpl) {
  return async (url, { method = 'GET', headers = {}, body, ms = 9000, max = MAX_UPSTREAM, binary = false } = {}) => {
    const ctl = new AbortController(), timer = setTimeout(() => ctl.abort(), ms);
    try {
      const r = await fetchImpl(url, { method, body, headers: { 'User-Agent': UA, Accept: '*/*', ...headers }, signal: ctl.signal, redirect: 'follow' });
      if (!r.ok) throw new Error(`${new URL(url).hostname}: HTTP ${r.status}`);
      if (+(r.headers.get('content-length') || 0) > max) throw new Error(`${new URL(url).hostname}: answer too large`);
      const buf = new Uint8Array(await r.arrayBuffer());
      if (buf.length > max) throw new Error(`${new URL(url).hostname}: answer too large`);
      return binary ? buf : new TextDecoder().decode(buf);
    } catch (e) { throw e?.name === 'AbortError' ? new Error(`${new URL(url).hostname}: timed out`) : e; }
    finally { clearTimeout(timer); }
  };
}

export async function handle(request, env = {}, ctx = {}, fetchImpl = globalThis.fetch, cacheImpl = globalThis.caches?.default) {
  env = env || {}; ctx = ctx || {};
  if (request.method !== 'GET') return json(405, { error: 'GET only' }, { Allow: 'GET', 'Cache-Control': 'no-store' });
  const url = new URL(request.url), id = String(url.searchParams.get('id') || '').slice(0, 24), feed = Object.hasOwn(FEEDS, id) ? FEEDS[id] : null;
  if (!id || id === 'list') return json(200, { feeds: Object.entries(FEEDS).map(([k, f]) => ({ id: k, refreshedEverySeconds: f.ttl, gives: f.what, ...(f.params ? { needs: 'lat, lon' } : {}) })) }, { 'Cache-Control': 'public, max-age=3600' });
  if (!feed) return json(404, { error: 'Unknown feed', feeds: Object.keys(FEEDS) }, { 'Cache-Control': 'public, max-age=3600' });
  const params = feed.params ? feed.params(url.searchParams) : {};
  if (!params) return json(400, { error: 'This feed needs lat (−90…90) and lon (−180…180)' }, { 'Cache-Control': 'public, max-age=3600' });
  const canon = Object.entries(params).map(([k, v]) => `${k}=${v}`).join('&'), now = () => (typeof ctx.now === 'function' ? ctx.now() : Date.now());
  const freshKey = new Request(`${url.origin}/api/feed/__fresh/${id}?${canon}`), staleKey = new Request(`${url.origin}/api/feed/__last/${id}?${canon}`), kvKey = `feed:${id}:${canon}`;
  const later = (p) => { const q = Promise.resolve(p).catch(() => {}); if (typeof ctx.waitUntil === 'function') ctx.waitUntil(q); return q; };
  const store = (key, body, fetched, seconds, extra = {}) => cacheImpl ? later(cacheImpl.put(key, new Response(body, { headers: { 'Content-Type': 'application/json; charset=utf-8', 'Cache-Control': `public, max-age=${seconds}, s-maxage=${seconds}`, 'X-Feed-Fetched': String(fetched), ...extra } }))) : Promise.resolve();
  // `Vary: *` keeps browsers and the app's service worker from storing an answer (the Cache API refuses such responses), so a
  // page served from this same host never sees a copy older than the one held here.
  const reply = (body, fetched, state, status = 200) => {
    const age = Math.max(0, Math.round((now() - fetched) / 1000)), stale = state === 'stale';
    return json(status, body, { 'Cache-Control': status === 200 && !stale ? `public, max-age=300, s-maxage=${Math.max(60, feed.ttl - age)}` : 'no-store', Vary: '*', Age: String(age), 'X-Feed-Cache': state, 'X-Feed-Age': String(age), 'X-Feed-Fetched': new Date(fetched).toISOString(), ...(stale ? { 'X-Feed-Stale': '1', Warning: '110 - "upstream unavailable; last good answer"' } : {}) });
  };

  // 1. fresh copy at this location
  const hit = cacheImpl ? await cacheImpl.match(freshKey).catch(() => null) : null;
  if (hit) {
    const fetched = +hit.headers.get('X-Feed-Fetched') || now(), failed = hit.headers.get('X-Feed-Failed');
    if (failed) return reply(await hit.text(), fetched, 'error', 502);
    return reply(await hit.text(), fetched, hit.headers.get('X-Feed-Stale') ? 'stale' : 'hit');
  }
  // 2. copy shared between locations (optional KV binding)
  const kv = env.FEED_STORE && typeof env.FEED_STORE.get === 'function' ? env.FEED_STORE : null;
  let shared = null;
  if (kv) { try { shared = await kv.get(kvKey, 'json'); } catch { shared = null; } if (shared && !(typeof shared.body === 'string' && fin(shared.t))) shared = null; }
  if (shared && now() - shared.t < feed.ttl * 1000) {
    await store(freshKey, shared.body, shared.t, Math.max(60, Math.round(feed.ttl - (now() - shared.t) / 1000)));
    return reply(shared.body, shared.t, 'store');
  }
  // 3. the upstream, once
  const flight = `${id}?${canon}`;
  let job = inflight.get(flight);
  if (!job) {
    job = (async () => {
      const r = await feed.load(reader(fetchImpl), env, params), fetched = now();
      const body = JSON.stringify({ id, fetched: new Date(fetched).toISOString(), refreshedEverySeconds: feed.ttl, source: r.source, terms: r.terms, ...(r.note ? { note: r.note } : {}), data: r.data });
      if (body.length > MAX_OUT) throw new Error('answer larger than the size cap');
      return { body, fetched };
    })();
    inflight.set(flight, job); job.finally(() => inflight.delete(flight)).catch(() => {});
  }
  try {
    const { body, fetched } = await job;
    await Promise.all([store(freshKey, body, fetched, feed.ttl), store(staleKey, body, fetched, STALE_KEEP)]);
    if (kv) later(kv.put(kvKey, JSON.stringify({ t: fetched, body }), { expirationTtl: STALE_KEEP }));
    return reply(body, fetched, 'miss');
  } catch (e) {
    // 4. the last good answer, with its age; the failure is remembered briefly so the upstream is left alone
    const old = cacheImpl ? await cacheImpl.match(staleKey).catch(() => null) : null;
    const last = old ? { body: await old.text(), t: +old.headers.get('X-Feed-Fetched') || 0 } : shared;
    if (last) { await store(freshKey, last.body, last.t, RETRY_WHILE_STALE, { 'X-Feed-Stale': '1' }); return reply(last.body, last.t, 'stale'); }
    const body = JSON.stringify({ id, error: 'Upstream unavailable and no stored answer yet', reason: String(e?.message || e).slice(0, 160) });
    await store(freshKey, body, now(), RETRY_AFTER_FAILURE, { 'X-Feed-Failed': '1' });
    return reply(body, now(), 'error', 502);
  }
}

export const onRequest = (c) => handle(c.request, c.env, c, fetch, caches.default);
