// Upstream petroleum fiscal terms by country: the single table of the app (the site data, the advisor and the economics
// suite read it through `fiscalOf`). Every figure is read from the publication cited on its row (addresses opened on the
// `retrieved` date). Only what a cited page states is recorded: a term no page gives is null, never a guess, and
// countries for which no open publication was found are absent (the site data then fall back to the statutory corporate
// income tax rate of the OECD table, live or from the bundled snapshot).
//
//   regime        'tax-royalty' (licence / concession), 'psc' (production sharing), 'mixed', or null where the source does not say
//   royalty       % of gross production (headline rate; alternatives in `note`)
//   petroleumTax  % special petroleum / resource-rent / hydrocarbon tax (on top of, or in place of, corporate income tax — see `note`)
//   corporateTax  % income tax applying to upstream companies as stated by the source
//   marginalTake  % combined marginal rate on upstream profit, only where the source publishes it or it follows from the
//                 stated rates by the formula shown in `note`
//   costOilCap    % of production from which cost may be recovered (production sharing)
//   profitSplit   % of profit oil going to the contractor (production sharing)
//   headline      % rate offered to the economics suite as `taxRate` (what `basis` says it is); null = use the statutory
//                 corporate income tax rate
//   year          year of the source's last review (of the headline's source where a row has several)
// RANGES. Where a source gives a range, the row is written [low, high]; the table then holds the MIDPOINT under the
// plain name and the bounds as `<name>Low` and `<name>High` (royaltyLow, headlineHigh, …), and `basis` says so.
// TERRAIN. Where the law sets terms by water depth, `terrain` lists the classes (minDepth exclusive, maxDepth
// inclusive, metres; maxDepth 0 = on land) with the terms that differ; `fiscalOf(code, { waterDepth })` applies them.
// MERGE. The economics suite kept a second table (js/data/costbasis.js, 26 countries). Its countries and fields are
// merged here; `merged` on a row says what came from it and how a disagreement was settled (the more recent or the
// more official source wins).
// Real terms depend on the licence or contract, the terrain and the vintage: every term stays editable in the economics suite.
const PWC = (slug, page, country, reviewed) => ({ citation: `PwC Worldwide Tax Summaries, ${country}, Corporate — ${page === 'other-taxes' ? 'Other taxes' : 'Taxes on corporate income'} (last reviewed ${reviewed})`, url: `https://taxsummaries.pwc.com/${slug}/corporate/${page}`, retrieved: '2026-10-08' });
const CIT = 'taxes-on-corporate-income', OTH = 'other-taxes';
const P2 = (slug, country, reviewed) => ({ ...PWC(slug, CIT, country, reviewed), retrieved: '2026-10-09' });
const EY = (chapter) => ({ citation: `EY, Global oil and gas tax guide 2019 (law as at 1 January 2019), chapter “${chapter}”`, url: 'https://www.ey.com/content/dam/ey-unified-site/ey-com/en-gl/technical/tax/documents/ey-global-oil-and-gas-tax-guide-2019.pdf', retrieved: '2026-10-09' });

const ROWS = {
  NO: { country: 'Norway', regime: 'tax-royalty', royalty: 0, petroleumTax: 71.8, corporateTax: 22, marginalTake: 78, headline: 78, basis: 'combined marginal rate: 22% company tax plus special petroleum tax', year: 2026,
    note: 'Special tax 71.8% is levied on the base after deduction of the 22% company tax (56% effective), giving 78% combined; royalties are no longer part of the system.',
    source: { citation: 'Norwegian Offshore Directorate and Ministry of Energy, Norwegian Petroleum: “The petroleum tax system” (page updated 8 October 2026)', url: 'https://www.norskpetroleum.no/en/economy/petroleum-tax/', retrieved: '2026-10-08' } },
  GB: { country: 'United Kingdom', regime: 'tax-royalty', royalty: 0, petroleumTax: 48, corporateTax: 30, marginalTake: 78, headline: 78, basis: 'marginal rate on UK and UKCS extraction income: ring fence corporation tax 30% + supplementary charge 10% + energy profits levy 38%', year: 2026,
    note: 'Petroleum tax shown = supplementary charge 10% + energy profits levy 38%; petroleum revenue tax is permanently 0%. No royalty is charged (recorded in the economics table from the same authority’s taxation page). The energy profits levy is legislated to end on 31 March 2030 or earlier, after which the rate is 40%.',
    merged: 'royalty (none) and the end date of the levy from the economics table; rates identical in both tables',
    source: { citation: 'North Sea Transition Authority, “Taxation — overview” (2026)', url: 'https://www.nstauthority.co.uk/exploration-production/taxation/overview/', retrieved: '2026-10-08' } },
  DK: { country: 'Denmark', regime: 'tax-royalty', royalty: 0, petroleumTax: 52, corporateTax: 25, marginalTake: 64, headline: 64, basis: 'effective rate: 25% ring-fenced company tax, deductible against the 52% hydrocarbon tax', year: 2026,
    note: 'The ordinary 22% rate does not apply upstream; 25% + 52% × (1 − 0.25) = 64%. No royalty (EY guide 2019).', source: PWC('denmark', CIT, 'Denmark', '31 July 2026'), source2: EY('Denmark'),
    merged: 'royalty (none) from the economics table (EY guide 2019); rates identical in both tables' },
  NG: { country: 'Nigeria', regime: 'tax-royalty', royalty: 12.5, petroleumTax: 30, corporateTax: 30, marginalTake: 60, headline: 60, basis: 'highest headline rate under the Petroleum Industry Act: 30% companies income tax + 30% hydrocarbon tax (onshore and shallow-water mining leases)', year: 2026,
    terrain: [
      { name: 'onshore', maxDepth: 0, royalty: 15, petroleumTax: 30, marginalTake: 60, headline: 60 },
      { name: 'shallow water (to 200 m)', minDepth: 0, maxDepth: 200, royalty: 12.5, petroleumTax: 30, marginalTake: 60, headline: 60 },
      { name: 'deep offshore (beyond 200 m)', minDepth: 200, royalty: 7.5, petroleumTax: 0, marginalTake: 30, headline: 30, basis: 'companies income tax 30%; the hydrocarbon tax of the Petroleum Industry Act does not apply to deep offshore (the 2026 source adds that the Nigeria Tax Act may extend it there at a rate not yet clear)' },
    ],
    terrainBasis: 'Petroleum Industry Act 2021: deep offshore is water deeper than 200 m (s. 318, as quoted by Ole & Herbert 2022); the 200 m upper limit of shallow water follows from that definition',
    note: 'Royalty by terrain: onshore 15%, shallow water 12.5%, deep offshore and frontier 7.5%, plus a price-based royalty of 0–10% (0 below US$50/bbl, 5% at US$100, 10% above US$150). Hydrocarbon tax 30% for mining leases and 15% for prospecting licences and marginal fields onshore and in shallow water; it does not apply to deep offshore, where the headline is the 30% companies income tax. Licences not converted to the Act stay under petroleum profit tax (50% PSC, 65.75%/85% other). Royalty also slides with production: deep offshore 5% up to 50,000 bbl/d, onshore and shallow water 5% up to 5,000 bbl/d and 7.5% for the next 5,000 bbl/d (Petroleum Royalty Regulations 2022); frontier basins 7.5%. A cost-price-ratio limit caps deductible cost at 65% of gross revenue.',
    merged: 'the economics table holds the deep-offshore case (royalty 7.5%, tax 30%), this table the onshore and shallow-water case (12.5%, 60%): both are kept and selected by water depth; the cost-price-ratio limit comes from the economics table',
    source: PWC('nigeria', CIT, 'Nigeria', '29 May 2026'),
    source2: { citation: 'EY Global Tax News, “Nigerian Government signs Petroleum Industry Bill 2020 into law” (2 September 2021) — royalty rates', url: 'https://globaltaxnews.ey.com/news/2021-5913-nigerian-government-signs-petroleum-industry-bill-2020-into-law', retrieved: '2026-10-08' },
    source3: { citation: 'Federal Republic of Nigeria Official Gazette No. 205 (22 November 2022), S.I. No. 73, Petroleum Royalty Regulations 2022 (regulator’s copy, Internet Archive)', url: 'https://web.archive.org/web/20240909172300id_/https://www.nuprc.gov.ng/wp-content/uploads/2022/11/Petroleum-Royalty-Regulations-2022-pdf-1.pdf', retrieved: '2026-10-09' },
    source4: { citation: 'N. C. Ole and E. B. Herbert, Studia Iuridica Lublinensia 31(3), 2022, doi:10.17951/sil.2022.31.3.143-163 — quotes the deep-offshore definition of s. 318 of the Act; Chartered Institute of Taxation of Nigeria, “An Overview of the Petroleum Industry Act 2021” — hydrocarbon tax not applied to deep offshore', url: 'https://www.journals.umcs.pl/sil/article/download/13580/pdf', retrieved: '2026-10-09' } },
  AO: { country: 'Angola', regime: 'mixed', royalty: 20, petroleumTax: 50, corporateTax: null, marginalTake: null, costOilCap: 50, headline: 50, basis: 'petroleum income tax, regular rate under production sharing agreements', year: 2024,
    note: 'Regular rates: petroleum income tax 50% under production sharing agreements and 65.75% under association agreements; petroleum production tax (royalty) 20% under association agreements. Reduced rates apply to incremental production (Presidential Decree 8/24). The State profit-oil share of a production sharing agreement comes on top. Production sharing agreements: no royalty, cost-oil cap 50% (it may rise to 65% if development cost is not recovered in four to five years); the profit split slides with the contractor’s rate of return and is set in each agreement (EY guide 2019).',
    merged: 'cost-oil cap from the economics table (EY guide 2019); the petroleum income tax rate is the same in both tables, here from the more recent source', source2: EY('Angola'),
    source: { citation: 'Mayer Brown, “Angola: Incremental Production Decree and Other Ongoing Developments” (25 November 2024)', url: 'https://www.mayerbrown.com/en/insights/publications/2024/11/angola-incremental-production-decree-and-other-ongoing-developments', retrieved: '2026-10-08' } },
  AU: { country: 'Australia', regime: 'tax-royalty', royalty: 0, petroleumTax: 40, corporateTax: 30, marginalTake: 58, headline: 58, basis: 'petroleum resource rent tax 40%, deductible against the 30% company tax', year: 2026,
    note: 'Offshore projects: 40% + (1 − 0.40) × 30% = 58% (computed from the two stated rates; the source states that PRRT payments are deductible for income tax). No royalty offshore in Commonwealth waters; onshore and state projects pay 10–12.5% (EY guide 2019).',
    merged: 'offshore royalty (none) from the economics table (EY guide 2019); the 58% combined rate is the same in both tables',
    source: PWC('australia', OTH, 'Australia', '1 September 2026'), source2: PWC('australia', CIT, 'Australia', '1 September 2026') },
  US: { country: 'United States', regime: 'tax-royalty', royalty: 12.5, petroleumTax: null, corporateTax: 21, marginalTake: null, headline: null, basis: '', year: 2025,
    note: 'Federal offshore leases (Gulf lease sales from December 2025): 12.5% royalty for shallow and deep water, the lowest rate permitted by statute. Federal corporate tax 21%; state taxes come on top onshore (the combined statutory rate is taken from the OECD table); leases of 2008 and later carry 18.75%, some earlier ones 16.667% (economics table, same BOEM release).',
    merged: 'royalty of earlier lease sales from the economics table; same source and rates in both tables',
    source: { citation: 'US Bureau of Ocean Energy Management, press release “BOEM Advances First Two OBBBA Offshore Lease Sales” (November 2025)', url: 'https://www.boem.gov/newsroom/press-releases/boem-advances-first-two-obbba-offshore-lease-sales', retrieved: '2026-10-08' },
    source2: PWC('united-states', CIT, 'United States', '4 September 2026') },
  BR: { country: 'Brazil', regime: 'tax-royalty', royalty: 10, petroleumTax: null, corporateTax: 34, marginalTake: null, headline: 34, basis: 'corporate income tax 15% + 10% surcharge + 9% social contribution on net income', year: 2026,
    note: 'Concession regime of Law 9.478/1997 (pre-salt production sharing is governed by a separate law and is not covered by this row). Royalty: 10% of production (art. 47; the regulator may reduce it to 5% in the tender). The special participation on large fields (10–40% of net revenue, by production volume) is not included here; pre-salt production sharing pays a 15% royalty and a bid profit-oil share (ANP press kit 2019).',
    merged: 'special participation range and production-sharing royalty from the economics table (ANP, “Participações governamentais”, 2019); royalty and tax identical in both tables',
    source: { citation: 'Brazil, Lei nº 9.478 de 6 de agosto de 1997 (Lei do Petróleo), art. 47 — Presidência da República, consolidated text', url: 'https://www.planalto.gov.br/ccivil_03/leis/l9478.htm', retrieved: '2026-10-08' },
    source2: PWC('brazil', CIT, 'Brazil', '23 September 2026') },
  NL: { country: 'Netherlands', regime: 'tax-royalty', royalty: [0, 7], petroleumTax: 50, corporateTax: null, marginalTake: null, headline: 50, basis: 'State profit share under the Mining Act', year: 2024,
    note: 'Holders of a production licence for oil or natural gas owe a State profit share of 50% of the profit, besides royalty, surface rights fee and provincial contributions (rates of those are not given by the source). Royalty 0–7% (EY guide 2019; midpoint shown); corporate income tax is credited against the State profit share.',
    merged: 'royalty range from the economics table (EY guide 2019). Its 25% corporate income tax of 2019 is not carried over: the OECD table has the current rate',
    source: { citation: 'Loyens & Loeff, “Out now: Levies under the Mining Act” (26 June 2024)', url: 'https://www.loyensloeff.com/insights/news--events/news/out-now-levies-under-the-mining-act/', retrieved: '2026-10-08' }, source2: EY('Netherlands') },
  MY: { country: 'Malaysia', regime: 'psc', royalty: 10, petroleumTax: 38, corporateTax: null, marginalTake: null, headline: 38, basis: 'petroleum income tax (in place of corporate income tax)', year: 2026,
    note: 'An effective rate of 25% applies to marginal fields; no other taxes are imposed on income from petroleum operations. Production sharing contracts with PETRONAS: royalty 10%; the cost-oil ceiling and the profit split follow the revenue-over-cost ratio of each contract and are not published (EY guide 2019).', source: PWC('malaysia', CIT, 'Malaysia', '16 June 2026'), source2: EY('Malaysia'),
    merged: 'regime and royalty from the economics table (EY guide 2019); petroleum income tax 38% in both tables, here from the more recent source' },
  TT: { country: 'Trinidad and Tobago', regime: 'tax-royalty',
    terrain: [{ name: 'deep water (beyond 400 m)', minDepth: 400, petroleumTax: 30, marginalTake: 35, headline: 35, basis: 'petroleum profits tax for deep-sea production 30% + unemployment levy 5% of taxable profits' }],
    terrainBasis: 'deep-water blocks are those in more than 400 m of water (EY guide 2019); the 30% rate is that of the 2026 source', royalty: 12.5, petroleumTax: 50, corporateTax: null, marginalTake: 55, headline: 55, basis: 'petroleum profits tax 50% + unemployment levy 5% of taxable profits', year: 2026,
    note: 'Deep-sea production is taxed at 30%. A supplementary petroleum tax on gross crude income (less royalties) also applies and is deductible for petroleum profits tax. Royalty 12.5% (EY guide 2019).', source: PWC('trinidad-and-tobago', CIT, 'Trinidad and Tobago', '2 June 2026'), source2: EY('Trinidad and Tobago'),
    merged: 'regime and royalty from the economics table (EY guide 2019). Conflict: that table has 35% petroleum profits tax for deep water (2019); the 2026 source states 30%, which is kept as the more recent' },
  TH: { country: 'Thailand', regime: 'mixed', royalty: null, petroleumTax: 50, corporateTax: null, marginalTake: null, headline: 50, basis: 'petroleum income tax on concessionaires (in place of corporate income tax)', year: 2026,
    note: 'Concession, production sharing contract or service contract; concessionaires pay 50% of annual net profit from petroleum operations, with royalties deductible.', source: PWC('thailand', CIT, 'Thailand', '24 August 2026') },
  BN: { country: 'Brunei Darussalam', regime: null, royalty: null, petroleumTax: 55, corporateTax: null, marginalTake: null, headline: 55, basis: 'petroleum tax on exploration and production profit', year: 2026, note: '', source: PWC('brunei-darussalam', CIT, 'Brunei Darussalam', '3 August 2026') },
  OM: { country: 'Oman', regime: 'psc', royalty: 0, petroleumTax: 55, corporateTax: null, marginalTake: null, headline: 55, basis: 'petroleum income tax rate specified for companies selling petroleum', year: 2026,
    note: 'Applied as set out in each exploration and production sharing agreement; the government pays the company’s tax out of its own share of production, so the tax is not actually borne by the company. No royalty; cost-oil cap and split are set in each agreement (EY guide 2019).', source: PWC('oman', CIT, 'Oman', '7 July 2026'), source2: EY('Oman'),
    merged: 'royalty (none) from the economics table (EY guide 2019); 55% in both tables' },
  SA: { country: 'Saudi Arabia', regime: 'tax-royalty', royalty: 15, petroleumTax: [50, 85], corporateTax: null, marginalTake: null, headline: [50, 85], basis: 'income tax on oil and hydrocarbon production: midpoint of the published 50–85% range', year: 2026,
    note: 'Income from oil and hydrocarbon production is taxed at 50% to 85%; natural-gas investment falls under the general 20% rate. The source does not list the steps of the range (the rate depends on the company’s capital investment). Crude-oil royalty 15% up to US$70/bbl, 45% of the increment from US$70 to US$100 and 80% above US$100 (first band shown).', source: PWC('saudi-arabia', CIT, 'Saudi Arabia', '29 July 2026'),
    source2: { citation: 'Arab Gulf States Institute, “Aramco and the Saudi Government Budget” (2025) — royalty bands', url: 'https://agsi.org/analysis/aramco-and-the-saudi-government-budget/', retrieved: '2026-10-08' },
    merged: 'regime and royalty bands from the economics table (Arab Gulf States Institute, 2025). Conflict: that table uses 50% as the single tax rate; the published range is 50–85% and, by the owner’s decision for ranges, the midpoint is the headline here with both bounds kept' },
  QA: { country: 'Qatar', regime: null, royalty: null, petroleumTax: 35, corporateTax: null, marginalTake: null, headline: 35, basis: 'minimum income tax rate for oil operations', year: 2026, note: 'The rate for oil operations may not be less than 35%; the agreement with the State can set a higher one.', source: PWC('qatar', CIT, 'Qatar', '17 September 2026') },
  IQ: { country: 'Iraq', regime: null, royalty: null, petroleumTax: 35, corporateTax: null, marginalTake: null, headline: 35, basis: 'income tax on foreign oil companies and their subcontractors', year: 2026, note: 'The general corporate rate is 15%.', source: PWC('iraq', CIT, 'Iraq', '24 June 2026') },
  EG: { country: 'Egypt', regime: null, royalty: null, petroleumTax: 40.55, corporateTax: null, marginalTake: null, headline: 40.55, basis: 'income tax on oil exploration companies', year: 2026, note: '', source: PWC('egypt', CIT, 'Egypt', '17 August 2026') },
  GH: { country: 'Ghana', regime: 'tax-royalty', royalty: [3, 12.5], petroleumTax: null, corporateTax: 35, marginalTake: null, headline: 35, basis: 'corporate income tax rate for upstream petroleum companies', year: 2026, note: 'The general corporate rate is 25%. Royalty is set in each petroleum agreement between 3% and 12.5% (EY guide 2019; midpoint shown); the State’s carried interest of at least 15% is not a tax and is not included.', source: PWC('ghana', CIT, 'Ghana', '11 March 2026'), source2: EY('Ghana'),
    merged: 'regime and royalty range from the economics table (EY guide 2019); 35% in both tables' },
  GA: { country: 'Gabon', regime: 'psc', royalty: [6, 12], petroleumTax: null, corporateTax: 35, marginalTake: null, costOilCap: [65, 75], profitSplit: 50, headline: 35, basis: 'corporate income tax rate for the oil and mining sectors', year: 2026, note: 'The general corporate rate is 30%. Exploitation and production sharing contracts under the 2014 code: royalty 6–12%, cost-oil cap 65–75% (midpoints shown), at most 50% of profit oil to the contractor (EY guide 2019).', source: PWC('gabon', CIT, 'Gabon', '6 August 2026'), source2: EY('Gabon'),
    merged: 'regime, royalty, cost-oil cap and profit split from the economics table (EY guide 2019); 35% in both tables' },
  MZ: { country: 'Mozambique', regime: 'psc', royalty: 10, petroleumTax: null, corporateTax: 32, marginalTake: null, costOilCap: 60, headline: 32, basis: 'corporate income tax', year: 2026,
    note: 'Petroleum production tax under a concession agreement: 10% for crude oil and condensate, 6% for natural gas and LNG. Concession contracts recover cost from at most 60% of production and share the remaining profit petroleum by an R-factor set in each contract (EY guide 2019).', source: PWC('mozambique', OTH, 'Mozambique', '11 August 2026'), source2: PWC('mozambique', CIT, 'Mozambique', '11 August 2026'), source3: EY('Mozambique'),
    merged: 'cost-recovery cap from the economics table (EY guide 2019). Conflict: that table uses the midpoint 8% of a 6–10% production tax; the 2026 source gives 10% for oil and 6% for gas, so the oil rate is kept. Regime changed from licence to production sharing because the contracts share profit petroleum' },
  PG: { country: 'Papua New Guinea', regime: 'tax-royalty', royalty: 2, petroleumTax: null, corporateTax: 30, marginalTake: null, headline: 30, basis: 'corporate income tax', year: 2026,
    note: 'Royalty 2% of the wellhead value; new petroleum projects also pay a tax-deductible development levy at the same 2%.', source: PWC('papua-new-guinea', OTH, 'Papua New Guinea', '23 September 2026'), source2: PWC('papua-new-guinea', CIT, 'Papua New Guinea', '23 September 2026') },
  VN: { country: 'Vietnam', regime: null, royalty: null, petroleumTax: null, corporateTax: [25, 50], marginalTake: null, headline: [25, 50], basis: 'corporate income tax for the oil and gas industry: midpoint of the published 25–50% range (the rate is set contract by contract)', year: 2026, note: 'The rate is set contract by contract between 25% and 50%; the standard rate is 20%.', source: PWC('vietnam', CIT, 'Vietnam', '23 September 2026') },
  KZ: { country: 'Kazakhstan', regime: null, royalty: [5, 18], petroleumTax: null, corporateTax: 20, marginalTake: null, headline: 20, basis: 'general corporate income tax', year: 2026,
    note: 'Mineral extraction tax on crude oil and condensate 5–18% depending on annual production (shown as royalty: midpoint, with both bounds); excess profit tax is progressive, 10–60% of the net income above 25% of deductions (not shown as a rate: it applies to a slice of income only); an alternative tax of 0–42% (3–14% for complex marine projects) can replace it and the extraction tax.', source: PWC('kazakhstan', OTH, 'Kazakhstan', '23 July 2026'), source2: PWC('kazakhstan', CIT, 'Kazakhstan', '23 July 2026') },
  CO: { country: 'Colombia', regime: null, royalty: null, petroleumTax: [0, 15], corporateTax: 35, marginalTake: null, headline: [35, 50], basis: 'corporate income tax 35% plus the price-dependent surtax on oil extraction (0 to 15 points): midpoint of the resulting 35–50% range', year: 2026,
    note: 'A surtax of 5% to 15% applies to oil extraction when the average price of the year is at or above 65% of the average of the preceding 120 months, so the rate can reach 50%.', source: PWC('colombia', CIT, 'Colombia', '21 July 2026') },
  // ---- rows added from the economics table and from sources opened on 9 October 2026 ---------------------------------
  CA: { country: 'Canada', regime: 'tax-royalty', royalty: [1, 7.5], petroleumTax: null, corporateTax: [26.5, 31], marginalTake: null, headline: [26.5, 31], basis: 'federal corporate income tax 15% plus a provincial rate of 11.5–16%: midpoint of the resulting 26.5–31% range', year: 2026,
    note: 'Newfoundland and Labrador offshore, generic royalty regime: basic royalty stepping from 1% to 7.5% with the R-factor (midpoint shown), and a net royalty of 10–50% for R between 1 and 3, which is not included. Other provinces have their own Crown royalty schedules. Provincial income tax range from the EY guide 2019 (2018 rates).',
    source: P2('canada', 'Canada', '12 June 2026'), source2: { citation: 'Government of Newfoundland and Labrador, Industry, Energy and Technology, “Generic Offshore Oil Royalty Regime” table (2017)', url: 'https://www.gov.nl.ca/em/files/royalties-2017-gorr-table.pdf', retrieved: '2026-10-08' }, source3: EY('Canada'),
    merged: 'royalty from the economics table (official provincial table). Conflict: that table assumes a single 30% income tax (15% + an assumed 15% provincial rate); the published provincial range is kept here with its midpoint' },
  MX: { country: 'Mexico', regime: 'mixed', royalty: null, petroleumTax: null, corporateTax: 30, marginalTake: null, costOilCap: 60, headline: 30, basis: 'federal corporate income tax', year: 2026,
    note: 'Licence contracts and production or profit sharing contracts coexist. Contracts add a price-linked royalty on revenue (floor 7.5% for oil) and a biddable additional royalty or State profit share; the 60% cost-recovery limit is that of bid rounds R1.01, R1.02, R2.01 and R3.01. Deep-water activities (more than 500 m) are ring-fenced for income tax without a different rate. Contract terms from the EY guide 2019.',
    source: P2('mexico', 'Mexico', '6 August 2026'), source2: EY('Mexico') },
  ID: { country: 'Indonesia', regime: 'psc', royalty: 0, petroleumTax: null, corporateTax: 22, marginalTake: null, costOilCap: 0, profitSplit: 43, headline: null, basis: '', year: 2026,
    note: 'Gross-split production sharing contracts (2017 rules): no royalty and no cost recovery; the contractor’s base share is 43% of oil and 48% of gas before variable and progressive adjustments (EY guide 2019; the rules were reissued in 2024 and not re-read). Income tax is computed under each contract; the general rate is 22%.',
    source: P2('indonesia', 'Indonesia', '11 June 2026'), source2: EY('Indonesia'),
    merged: 'regime, cost-recovery (none) and base split from the economics table (EY guide 2019). Conflict: that table uses 40% tax (25% plus 20% branch profits tax, 2019); the general rate is now 22% (2026 source) and no current source gives the upstream rate, so no petroleum-specific headline is offered and the statutory rate stands in' },
  KW: { country: 'Kuwait', regime: null, royalty: null, petroleumTax: null, corporateTax: 15, marginalTake: null, headline: 15, basis: 'flat corporate income tax on foreign companies', year: 2026,
    note: 'No separate oil and gas tax law is described by the source.', source: P2('kuwait', 'Kuwait', '22 July 2026'),
    merged: 'Conflict: the economics table carries a 15% royalty from the EY guide 2019; in that guide the 15% line is income tax on royalty income, not a production royalty, so no royalty is recorded here' },
  GY: { country: 'Guyana', regime: 'psc', royalty: 10, petroleumTax: null, corporateTax: 10, marginalTake: null, costOilCap: 65, profitSplit: 50, headline: 10, basis: 'corporate tax under the model production sharing agreement used since the 2022 licensing round', year: 2025,
    note: 'Model agreement: 10% royalty, 10% corporate tax, cost recovery up to 65%, profit oil shared 50/50. The 2016 Stabroek agreement differs: 2% royalty, cost recovery up to 75%, 50/50 split, and the contractor’s income tax settled by the government out of its share (economics table: Argus 2023 and EY guide 2019).',
    source: { citation: 'News Room Guyana, “Guyana to sign new offshore oil deal with Cybele for S7 shallow-water block” (5 December 2025) — terms of the model agreement as reported', url: 'https://newsroom.gy/2025/12/05/guyana-to-sign-new-offshore-oil-deal-with-cybele-for-s7-shallow-water-block/', retrieved: '2026-10-09' },
    source2: { citation: 'Department of Public Information, Guyana, “Gov’t reaffirms new model PSA will accrue more economic benefits for Guyana” — 2% royalty of the Stabroek agreement', url: 'https://dpi.gov.gy/govt-reaffirms-new-model-psa-will-accrue-more-economic-benefits-for-guyana-maintains-2016-stabroek-psa-will-remain-uncontended/', retrieved: '2026-10-09' },
    source3: { citation: 'Argus Media, “Guyana extends deadline for first oil block bids” (2023)', url: 'https://www.argusmedia.com/ja/news-and-insights/latest-market-news/2438577-guyana-extends-deadline-for-first-oil-block-bids', retrieved: '2026-10-08' },
    merged: 'Conflict: the economics table models the 2016 Stabroek agreement (royalty 2%, no tax on the contractor, cost oil 75%); this row carries the current model agreement, which applies to new blocks, and keeps the Stabroek terms in the note' },
  AR: { country: 'Argentina', regime: 'tax-royalty', royalty: 12, petroleumTax: null, corporateTax: [25, 35], marginalTake: null, headline: [25, 35], basis: 'corporate income tax on a progressive 25–35% scale: midpoint of the range (large companies pay 35% on the top slice)', year: 2026,
    note: 'Royalty 12% on federal and provincial concessions and a 12% export duty in 2019 (EY guide 2019).', source: P2('argentina', 'Argentina', '16 June 2026'), source2: EY('Argentina'),
    merged: 'royalty from the economics table (EY guide 2019). Conflict: that table has 30% income tax (the 2018–2019 rate); the 2026 scale of 25–35% is kept' },
  IN: { country: 'India', regime: 'mixed', royalty: 12.5, petroleumTax: null, corporateTax: null, marginalTake: null, headline: null, basis: '', year: 2025,
    terrain: [
      { name: 'onland', maxDepth: 0, royalty: 12.5 },
      { name: 'shallow water (to the 400 m isobath)', minDepth: 0, maxDepth: 400, royalty: 7.5 },
      { name: 'deep water (400 m to 1 500 m)', minDepth: 400, maxDepth: 1500, royalty: 5 },
      { name: 'ultra-deep water (beyond 1 500 m)', minDepth: 1500, royalty: 2 },
    ],
    terrainBasis: 'depth classes and crude-oil royalty of the Model Revenue Sharing Contract, Open Acreage Licensing round X (December 2025); deep and ultra-deep water pay no royalty for the first seven years',
    note: 'Revenue sharing contracts of the Hydrocarbon Exploration and Licensing Policy (no cost recovery; the government share of revenue is bid); older blocks are production sharing contracts. Onland natural gas pays 10%. A news report of 12 May 2026 describes a revised schedule (deep water 5% then 10%, ultra-deep 0% then 5%); the notification was not read, so the contract values stand.',
    source: { citation: 'Government of India, Directorate General of Hydrocarbons, Model Revenue Sharing Contract, Open Acreage Licensing Policy bid round X (December 2025; Internet Archive copy)', url: 'https://web.archive.org/web/20251231161600id_/https://online.dghindia.org/oalp/Files/pdf/MRSC_OALP-X.pdf', retrieved: '2026-10-09' },
    source2: { citation: 'Indian Economic Service, Arthapedia, “Hydrocarbon Exploration and Licensing Policy (HELP)” — onland royalty', url: 'https://ies.gov.in/arthapedia/concept/hydrocarbon-exploration-and-licensing-policy-help', retrieved: '2026-10-09' } },
  CN: { country: 'China', regime: 'psc', royalty: 6, petroleumTax: null, corporateTax: 25, marginalTake: null, costOilCap: [50, 62.5], headline: 25, basis: 'corporate income tax', year: 2026,
    note: 'Offshore production sharing contracts with foreign contractors: a resource tax of 6% of sales stands in place of royalty (the royalty of 0–12.5% applies only to contracts concluded before November 2011), cost recovery is limited to 50–62.5% offshore (midpoint shown), and a special oil gain levy applies above US$65/bbl; the allocation of profit oil is contract-specific (EY guide 2019).',
    source: P2('peoples-republic-of-china', 'China', '22 July 2026'), source2: EY('China'), merged: 'regime, resource tax and cost-recovery range from the economics table (EY guide 2019); 25% in both tables' },
  DZ: { country: 'Algeria', regime: 'mixed', royalty: 10, petroleumTax: [10, 50], corporateTax: 30, marginalTake: null, headline: 30, basis: 'income tax under hydrocarbon law 19-13 (the hydrocarbon revenue tax of 10–50% comes on top)', year: 2025,
    note: 'Law 19-13: hydrocarbon royalty 10%, hydrocarbon revenue tax 10–50% (midpoint shown), income tax 30%; reductions are possible with ministerial approval (royalty not below 5%). Under production sharing the foreign contractor’s share is capped at 49% of total production.',
    source: { citation: 'International Bar Association, “Algeria Bid Round 2024: a strategic shift under the new Algerian hydrocarbon law” (8 June 2025)', url: 'https://www.ibanet.org/algeria-bid-round-2024-hydrocarbon-law', retrieved: '2026-10-09' } },
  LY: { country: 'Libya', regime: 'psc', royalty: 16.67, petroleumTax: null, corporateTax: 20, marginalTake: 65, headline: 20, basis: 'statutory corporate income tax (under exploration and production sharing agreements the taxes are notional and deemed paid by the national oil company)', year: 2026,
    note: 'Notional royalty of 16.67% of production under production sharing; under the concession rules a surtax brings tax, royalty adjustment and other deductions to a composite 65% of profits (EY guide 2019; shown as the marginal take).', source: P2('libya', 'Libya', '31 May 2026'), source2: EY('Libya') },
  KE: { country: 'Kenya', regime: 'psc', royalty: null, petroleumTax: null, corporateTax: 30, marginalTake: null, headline: 30, basis: 'corporate income tax', year: 2026,
    note: 'The State’s take comes mainly through the sharing of profit oil under each contract; the opened sources give no percentages.', source: P2('kenya', 'Kenya', '17 July 2026'), source2: EY('Kenya') },
  TZ: { country: 'Tanzania', regime: 'psc', royalty: [7.5, 12.5], petroleumTax: [25, 35], corporateTax: 30, marginalTake: null, headline: 30, basis: 'corporate income tax (the additional profits tax of 25% or 35% comes on top on its own bases)', year: 2026,
    note: 'Royalty 12.5% for onshore and shelf areas and 7.5% for offshore areas (the source gives no depth limit, so the range is shown with its midpoint); additional profits tax at 25% or 35% (Petroleum Act 2015, EY guide 2019).', source: P2('tanzania', 'Tanzania', '9 September 2026'), source2: EY('Tanzania') },
  SN: { country: 'Senegal', regime: 'psc', royalty: [7, 10], petroleumTax: null, corporateTax: 30, marginalTake: null, costOilCap: [55, 70], headline: 30, basis: 'corporate income tax', year: 2026,
    terrain: [{ name: 'onshore', maxDepth: 0, costOilCap: 55 }, { name: 'offshore (shallow 60%, deep 65%, ultra-deep 70%)', minDepth: 0, costOilCap: [60, 70] }],
    terrainBasis: 'Petroleum Code 2019: recoverable cost is capped at 55% onshore, 60% in shallow, 65% in deep and 70% in ultra-deep offshore; the opened text gives no depth limits for the three offshore classes, so offshore sites show their range',
    note: 'Petroleum Code of 2019: royalty 7–10% on liquids (6% on gas), by zone. The split of profit oil follows a ratio set in each contract.',
    source: { citation: 'Sénégal, Loi n° 2019-03 du 1er février 2019 portant Code pétrolier — summary sheet on vie-publique.sn', url: 'https://www.vie-publique.sn/documents/10994/code-petrolier-du-senegal-loi-2019-03-2019', retrieved: '2026-10-09' }, source2: P2('senegal', 'Senegal', '7 August 2026'), source3: EY('Senegal'),
    merged: 'Conflict: the economics table has royalty 6–10% (EY guide 2019, midpoint 8) and one cost-oil cap of 65%; the summary of the law gives 7–10% for liquids and the caps by terrain, kept here as the more official source' },
  MR: { country: 'Mauritania', regime: 'psc', royalty: 0, petroleumTax: null, corporateTax: 25, marginalTake: null, headline: 25, basis: 'corporate income tax (the contract rate, at least the common rate)', year: 2026,
    note: 'Hydrocarbons Code of 2010: no royalty; the State may take a participation of not less than 10% (EY guide 2019). Permits under the 1988 ordinance pay a royalty of not less than 10%.', source: P2('mauritania', 'Mauritania', '7 August 2026'), source2: EY('Mauritania') },
  NA: { country: 'Namibia', regime: 'tax-royalty', royalty: 5, petroleumTax: null, corporateTax: 35, marginalTake: null, headline: 35, basis: 'petroleum income tax', year: 2019,
    note: 'An additional profits tax in three tranches is triggered at after-tax rates of return of 15% and 20–25%; its rates are not given by the source. Terms are those of the 2019 guide and may have changed.', source: EY('Namibia'),
    merged: 'row taken from the economics table (EY guide 2019), re-read in the guide' },
  SR: { country: 'Suriname', regime: 'psc', royalty: 6.25, petroleumTax: null, corporateTax: 36, marginalTake: null, headline: 36, basis: 'income tax under offshore production sharing contracts', year: 2026,
    note: 'Royalty 6.25% of gross production; the cost-oil ceiling is set by contract and profit oil is shared by an R-factor (no figures published). The state company states a total government take of 60–70% after costs.',
    source: { citation: 'Staatsolie Maatschappij Suriname N.V., Staatsolie Hydrocarbon Institute, “FAQ”', url: 'https://www.staatsolie.com/en/shi/faq/', retrieved: '2026-10-09' }, merged: 'same source and figures as the economics table' },
  GQ: { country: 'Equatorial Guinea', regime: 'psc', royalty: 13, petroleumTax: null, corporateTax: 25, marginalTake: null, headline: 25, basis: 'general corporate income tax (2025 source; an oil and gas guide of 2019 gave 35%)', year: 2025,
    note: 'Production sharing contracts under Hydrocarbon Law 8/2006: royalty of not less than 13% (the minimum is shown) and a State participation of not less than 20%; cost-oil cap and split are set in each contract (EY guide 2019).', source: P2('equatorial-guinea', 'Equatorial Guinea', '21 November 2025'), source2: EY('Equatorial Guinea'),
    merged: 'regime and royalty from the economics table (EY guide 2019). Conflict: that table has 35% income tax (2019 guide, oil and gas chapter); the 2025 source gives a general rate of 25% without a petroleum exception, kept as the more recent — the 35% may still apply under existing contracts' },
  CG: { country: 'Republic of the Congo', regime: 'psc', royalty: 15, petroleumTax: null, corporateTax: 30, marginalTake: null, costOilCap: 50, headline: 30, basis: 'standard corporate income tax', year: 2026,
    note: 'Hydrocarbon Code of 2016: mineral fee 15%, reduced to 12% in specific zones such as deep water; cost recovery limited to 50% of net production, up to 70% for especially difficult work; the State’s share of profit oil is not less than 35% (EY guide 2019; no depth limits stated).', source: P2('republic-of-congo', 'Republic of the Congo', '7 August 2026'), source2: EY('Republic of the Congo') },
  CM: { country: 'Cameroon', regime: 'mixed', royalty: null, petroleumTax: null, corporateTax: [33, 50], marginalTake: null, headline: [33, 50], basis: 'corporate income tax set in each petroleum contract: midpoint of the 33–50% range', year: 2019,
    note: 'Production sharing contract or concession; a royalty applies to concession holders at a contract rate and profit oil is shared by an R-factor (no figures in the source). Terms are those of the 2019 guide and may have changed.', source: EY('Cameroon') },
  CI: { country: 'Côte d’Ivoire', regime: 'psc', royalty: null, petroleumTax: null, corporateTax: 25, marginalTake: null, costOilCap: [70, 80], headline: 25, basis: 'corporate income tax (in recent contracts it is settled out of the State’s share of profit oil)', year: 2026,
    note: 'Cost oil is usually capped at 70–80% of production (midpoint shown); the profit split is contract-specific — the guide’s example gives the holder 55% up to 100,000 bbl/d falling to 40% above 300,000 bbl/d (EY guide 2019).', source: P2('ivory-coast', 'Côte d’Ivoire', '9 September 2026'), source2: EY('Côte d’Ivoire') },
  AZ: { country: 'Azerbaijan', regime: 'psc', royalty: 0, petroleumTax: null, corporateTax: [20, 32], marginalTake: null, headline: [20, 32], basis: 'profit tax negotiated in each production sharing agreement: midpoint of the 20–32% range (statutory rate 20%)', year: 2026,
    note: 'Contractors under production sharing agreements pay profit tax only and are exempt from royalties (EY guide 2019); each agreement has its own tax regime and the statutory profit tax is 20%.', source: P2('azerbaijan', 'Azerbaijan', '28 September 2026'), source2: EY('Azerbaijan') },
  TM: { country: 'Turkmenistan', regime: 'psc', royalty: null, petroleumTax: null, corporateTax: 20, marginalTake: null, headline: 20, basis: 'corporate income tax for entities operating under the Petroleum Law', year: 2021,
    note: 'Royalty and the split of production are set in each production sharing agreement (the published model agreement leaves the percentages blank). The source dates from January 2021.',
    source: { citation: 'Deloitte, “International Tax — Turkmenistan Highlights 2021” (January 2021)', url: 'https://www.deloitte.com/content/dam/assets-shared/docs/services/consulting/2024/dttl-tax-turkmenistanhighlights-2021.pdf', retrieved: '2026-10-09' } },
};

export const FISCAL_TERMS = ['royalty', 'petroleumTax', 'corporateTax', 'marginalTake', 'costOilCap', 'profitSplit', 'headline'];
/** [low, high] -> midpoint under the plain name, bounds as <name>Low / <name>High; a missing term becomes null (rows only, not terrain classes). */
function expand(row, fill) {
  const o = { ...row };
  for (const k of FISCAL_TERMS) {
    const v = row[k];
    if (Array.isArray(v)) { o[k] = +((v[0] + v[1]) / 2).toFixed(3); o[k + 'Low'] = v[0]; o[k + 'High'] = v[1]; }
    else if (v === undefined) { if (fill) o[k] = null; }
    else if (!fill) { o[k + 'Low'] = undefined; o[k + 'High'] = undefined; } // a single figure for a terrain class replaces the range of the row
  }
  if (row.terrain) o.terrain = row.terrain.map((t) => Object.freeze(expand(t, false)));
  return Object.freeze(o);
}
export const FISCAL = Object.freeze(Object.fromEntries(Object.entries(ROWS).map(([c, r]) => [c, expand(r, true)])));

/**
 * Fiscal terms of a country (ISO-2 code) or null.
 *   fiscalOf('NG')                       the row as tabled (its `terrain` lists the depth classes, if the law has any)
 *   fiscalOf('NG', { waterDepth: 1200 }) the row with the terms of the terrain class holding that water depth (metres;
 *                                        0 = on land) applied, plus `terrainName` and `terrainBasis`
 * Fields: regime, royalty, petroleumTax, corporateTax, marginalTake, costOilCap, profitSplit, headline (each a number
 * or null, with <name>Low / <name>High where the source gives a range), basis, year, note, source(s), merged.
 */
export function fiscalOf(code, opt) {
  const f = FISCAL[String(code || '').toUpperCase()] || null, depth = opt?.waterDepth;
  if (!f || !f.terrain || typeof depth !== 'number' || !(depth >= 0)) return f;
  const t = f.terrain.find((x) => depth > (x.minDepth ?? -1) && depth <= (x.maxDepth ?? Infinity));
  if (!t) return f;
  const { name, minDepth, maxDepth, ...terms } = t, o = { ...f, ...terms, terrainName: name, terrainBasis: f.terrainBasis || '' };
  for (const k of Object.keys(o)) if (o[k] === undefined) delete o[k]; // a class with a single figure clears the bounds of the row's range
  return o;
}
/** Plain-language name of a regime code. */
export const REGIME_LABELS = { 'tax-royalty': 'tax and royalty (licence)', psc: 'production sharing contract', mixed: 'licences and production sharing contracts' };
