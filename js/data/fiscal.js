// Upstream petroleum fiscal terms by country, read from the publications cited on each row (addresses opened on the
// `retrieved` date). Only what the cited page states is recorded: a term the page does not give is null, never a guess,
// and countries for which no open publication was found are simply absent (the site data then fall back to the statutory
// corporate income tax rate of the OECD table, live or from the bundled snapshot).
//
//   regime        'tax-royalty' (licence / concession), 'psc' (production sharing), 'mixed', or null where the source does not say
//   royalty       % of gross production (headline rate; alternatives in `note`)
//   petroleumTax  % special petroleum / resource-rent / hydrocarbon tax (on top of, or in place of, corporate income tax — see `note`)
//   corporateTax  % income tax applying to upstream companies as stated by the source
//   marginalTake  % combined marginal rate on upstream profit, only where the source publishes it or it follows from the
//                 stated rates by the formula shown in `note`
//   headline      % rate offered to the economics suite as `taxRate` (what `basis` says it is); null = use the statutory
//                 corporate income tax rate
//   year          year of the source's last review
// Real terms depend on the licence or contract, the terrain and the vintage: every term stays editable in the economics suite.
const PWC = (slug, page, country, reviewed) => ({ citation: `PwC Worldwide Tax Summaries, ${country}, Corporate — ${page === 'other-taxes' ? 'Other taxes' : 'Taxes on corporate income'} (last reviewed ${reviewed})`, url: `https://taxsummaries.pwc.com/${slug}/corporate/${page}`, retrieved: '2026-10-08' });
const CIT = 'taxes-on-corporate-income', OTH = 'other-taxes';

export const FISCAL = {
  NO: { country: 'Norway', regime: 'tax-royalty', royalty: 0, petroleumTax: 71.8, corporateTax: 22, marginalTake: 78, headline: 78, basis: 'combined marginal rate: 22% company tax plus special petroleum tax', year: 2026,
    note: 'Special tax 71.8% is levied on the base after deduction of the 22% company tax (56% effective), giving 78% combined; royalties are no longer part of the system.',
    source: { citation: 'Norwegian Offshore Directorate and Ministry of Energy, Norwegian Petroleum: “The petroleum tax system” (page updated 8 October 2026)', url: 'https://www.norskpetroleum.no/en/economy/petroleum-tax/', retrieved: '2026-10-08' } },
  GB: { country: 'United Kingdom', regime: 'tax-royalty', royalty: null, petroleumTax: 48, corporateTax: 30, marginalTake: 78, headline: 78, basis: 'marginal rate on UK and UKCS extraction income: ring fence corporation tax 30% + supplementary charge 10% + energy profits levy 38%', year: 2026,
    note: 'Petroleum tax shown = supplementary charge 10% + energy profits levy 38%; petroleum revenue tax is permanently 0%. The source does not mention a royalty.',
    source: { citation: 'North Sea Transition Authority, “Taxation — overview” (2026)', url: 'https://www.nstauthority.co.uk/exploration-production/taxation/overview/', retrieved: '2026-10-08' } },
  DK: { country: 'Denmark', regime: 'tax-royalty', royalty: null, petroleumTax: 52, corporateTax: 25, marginalTake: 64, headline: 64, basis: 'effective rate: 25% ring-fenced company tax, deductible against the 52% hydrocarbon tax', year: 2026,
    note: 'The ordinary 22% rate does not apply upstream; 25% + 52% × (1 − 0.25) = 64%.', source: PWC('denmark', CIT, 'Denmark', '31 July 2026') },
  NG: { country: 'Nigeria', regime: 'tax-royalty', royalty: 12.5, petroleumTax: 30, corporateTax: 30, marginalTake: 60, headline: 60, basis: 'highest headline rate under the Petroleum Industry Act: 30% companies income tax + 30% hydrocarbon tax (onshore and shallow-water mining leases)', year: 2026,
    note: 'Royalty by terrain: onshore 15%, shallow water 12.5%, deep offshore and frontier 7.5%, plus a price-based royalty of 0–10% (0 below US$50/bbl, 5% at US$100, 10% above US$150). Hydrocarbon tax 30% for mining leases and 15% for prospecting licences and marginal fields onshore and in shallow water; it does not apply to deep offshore, where the headline is the 30% companies income tax. Licences not converted to the Act stay under petroleum profit tax (50% PSC, 65.75%/85% other).',
    source: PWC('nigeria', CIT, 'Nigeria', '29 May 2026'),
    source2: { citation: 'EY Global Tax News, “Nigerian Government signs Petroleum Industry Bill 2020 into law” (2 September 2021) — royalty rates', url: 'https://globaltaxnews.ey.com/news/2021-5913-nigerian-government-signs-petroleum-industry-bill-2020-into-law', retrieved: '2026-10-08' } },
  AO: { country: 'Angola', regime: 'mixed', royalty: 20, petroleumTax: 50, corporateTax: null, marginalTake: null, headline: 50, basis: 'petroleum income tax, regular rate under production sharing agreements', year: 2024,
    note: 'Regular rates: petroleum income tax 50% under production sharing agreements and 65.75% under association agreements; petroleum production tax (royalty) 20% under association agreements. Reduced rates apply to incremental production (Presidential Decree 8/24). The State profit-oil share of a production sharing agreement comes on top.',
    source: { citation: 'Mayer Brown, “Angola: Incremental Production Decree and Other Ongoing Developments” (25 November 2024)', url: 'https://www.mayerbrown.com/en/insights/publications/2024/11/angola-incremental-production-decree-and-other-ongoing-developments', retrieved: '2026-10-08' } },
  AU: { country: 'Australia', regime: 'tax-royalty', royalty: null, petroleumTax: 40, corporateTax: 30, marginalTake: 58, headline: 58, basis: 'petroleum resource rent tax 40%, deductible against the 30% company tax', year: 2026,
    note: 'Offshore projects: 40% + (1 − 0.40) × 30% = 58% (computed from the two stated rates; the source states that PRRT payments are deductible for income tax).',
    source: PWC('australia', OTH, 'Australia', '1 September 2026'), source2: PWC('australia', CIT, 'Australia', '1 September 2026') },
  US: { country: 'United States', regime: 'tax-royalty', royalty: 12.5, petroleumTax: null, corporateTax: 21, marginalTake: null, headline: null, basis: '', year: 2025,
    note: 'Federal offshore leases (Gulf lease sales from December 2025): 12.5% royalty for shallow and deep water, the lowest rate permitted by statute. Federal corporate tax 21%; state taxes come on top (the combined statutory rate is taken from the OECD table).',
    source: { citation: 'US Bureau of Ocean Energy Management, press release “BOEM Advances First Two OBBBA Offshore Lease Sales” (November 2025)', url: 'https://www.boem.gov/newsroom/press-releases/boem-advances-first-two-obbba-offshore-lease-sales', retrieved: '2026-10-08' },
    source2: PWC('united-states', CIT, 'United States', '4 September 2026') },
  BR: { country: 'Brazil', regime: 'tax-royalty', royalty: 10, petroleumTax: null, corporateTax: 34, marginalTake: null, headline: 34, basis: 'corporate income tax 15% + 10% surcharge + 9% social contribution on net income', year: 2026,
    note: 'Concession regime of Law 9.478/1997 (pre-salt production sharing is governed by a separate law and is not covered by this row). Royalty: 10% of production (art. 47; the regulator may reduce it to 5% in the tender). The special participation on large fields is not included here.',
    source: { citation: 'Brazil, Lei nº 9.478 de 6 de agosto de 1997 (Lei do Petróleo), art. 47 — Presidência da República, consolidated text', url: 'https://www.planalto.gov.br/ccivil_03/leis/l9478.htm', retrieved: '2026-10-08' },
    source2: PWC('brazil', CIT, 'Brazil', '23 September 2026') },
  NL: { country: 'Netherlands', regime: 'tax-royalty', royalty: null, petroleumTax: 50, corporateTax: null, marginalTake: null, headline: 50, basis: 'State profit share under the Mining Act', year: 2024,
    note: 'Holders of a production licence for oil or natural gas owe a State profit share of 50% of the profit, besides royalty, surface rights fee and provincial contributions (rates of those are not given by the source).',
    source: { citation: 'Loyens & Loeff, “Out now: Levies under the Mining Act” (26 June 2024)', url: 'https://www.loyensloeff.com/insights/news--events/news/out-now-levies-under-the-mining-act/', retrieved: '2026-10-08' } },
  MY: { country: 'Malaysia', regime: null, royalty: null, petroleumTax: 38, corporateTax: null, marginalTake: null, headline: 38, basis: 'petroleum income tax (in place of corporate income tax)', year: 2026,
    note: 'An effective rate of 25% applies to marginal fields; no other taxes are imposed on income from petroleum operations.', source: PWC('malaysia', CIT, 'Malaysia', '16 June 2026') },
  TT: { country: 'Trinidad and Tobago', regime: null, royalty: null, petroleumTax: 50, corporateTax: null, marginalTake: 55, headline: 55, basis: 'petroleum profits tax 50% + unemployment levy 5% of taxable profits', year: 2026,
    note: 'Deep-sea production is taxed at 30%. A supplementary petroleum tax on gross crude income (less royalties) also applies and is deductible for petroleum profits tax.', source: PWC('trinidad-and-tobago', CIT, 'Trinidad and Tobago', '2 June 2026') },
  TH: { country: 'Thailand', regime: 'mixed', royalty: null, petroleumTax: 50, corporateTax: null, marginalTake: null, headline: 50, basis: 'petroleum income tax on concessionaires (in place of corporate income tax)', year: 2026,
    note: 'Concession, production sharing contract or service contract; concessionaires pay 50% of annual net profit from petroleum operations, with royalties deductible.', source: PWC('thailand', CIT, 'Thailand', '24 August 2026') },
  BN: { country: 'Brunei Darussalam', regime: null, royalty: null, petroleumTax: 55, corporateTax: null, marginalTake: null, headline: 55, basis: 'petroleum tax on exploration and production profit', year: 2026, note: '', source: PWC('brunei-darussalam', CIT, 'Brunei Darussalam', '3 August 2026') },
  OM: { country: 'Oman', regime: 'psc', royalty: null, petroleumTax: 55, corporateTax: null, marginalTake: null, headline: 55, basis: 'petroleum income tax rate specified for companies selling petroleum', year: 2026,
    note: 'Applied as set out in each exploration and production sharing agreement; the government pays the company’s tax out of its own share of production, so the tax is not actually borne by the company.', source: PWC('oman', CIT, 'Oman', '7 July 2026') },
  SA: { country: 'Saudi Arabia', regime: null, royalty: null, petroleumTax: 50, corporateTax: null, marginalTake: null, headline: 50, basis: 'income tax on oil and hydrocarbon production, lower bound of the published 50–85% range', year: 2026,
    note: 'Income from oil and hydrocarbon production is taxed at 50% to 85%; natural-gas investment falls under the general 20% rate.', source: PWC('saudi-arabia', CIT, 'Saudi Arabia', '29 July 2026') },
  QA: { country: 'Qatar', regime: null, royalty: null, petroleumTax: 35, corporateTax: null, marginalTake: null, headline: 35, basis: 'minimum income tax rate for oil operations', year: 2026, note: 'The rate for oil operations may not be less than 35%; the agreement with the State can set a higher one.', source: PWC('qatar', CIT, 'Qatar', '17 September 2026') },
  IQ: { country: 'Iraq', regime: null, royalty: null, petroleumTax: 35, corporateTax: null, marginalTake: null, headline: 35, basis: 'income tax on foreign oil companies and their subcontractors', year: 2026, note: 'The general corporate rate is 15%.', source: PWC('iraq', CIT, 'Iraq', '24 June 2026') },
  EG: { country: 'Egypt', regime: null, royalty: null, petroleumTax: 40.55, corporateTax: null, marginalTake: null, headline: 40.55, basis: 'income tax on oil exploration companies', year: 2026, note: '', source: PWC('egypt', CIT, 'Egypt', '17 August 2026') },
  GH: { country: 'Ghana', regime: null, royalty: null, petroleumTax: null, corporateTax: 35, marginalTake: null, headline: 35, basis: 'corporate income tax rate for upstream petroleum companies', year: 2026, note: 'The general corporate rate is 25%.', source: PWC('ghana', CIT, 'Ghana', '11 March 2026') },
  GA: { country: 'Gabon', regime: null, royalty: null, petroleumTax: null, corporateTax: 35, marginalTake: null, headline: 35, basis: 'corporate income tax rate for the oil and mining sectors', year: 2026, note: 'The general corporate rate is 30%.', source: PWC('gabon', CIT, 'Gabon', '6 August 2026') },
  MZ: { country: 'Mozambique', regime: 'tax-royalty', royalty: 10, petroleumTax: null, corporateTax: 32, marginalTake: null, headline: 32, basis: 'corporate income tax', year: 2026,
    note: 'Petroleum production tax under a concession agreement: 10% for crude oil and condensate, 6% for natural gas and LNG.', source: PWC('mozambique', OTH, 'Mozambique', '11 August 2026'), source2: PWC('mozambique', CIT, 'Mozambique', '11 August 2026') },
  PG: { country: 'Papua New Guinea', regime: 'tax-royalty', royalty: 2, petroleumTax: null, corporateTax: 30, marginalTake: null, headline: 30, basis: 'corporate income tax', year: 2026,
    note: 'Royalty 2% of the wellhead value; new petroleum projects also pay a tax-deductible development levy at the same 2%.', source: PWC('papua-new-guinea', OTH, 'Papua New Guinea', '23 September 2026'), source2: PWC('papua-new-guinea', CIT, 'Papua New Guinea', '23 September 2026') },
  VN: { country: 'Vietnam', regime: null, royalty: null, petroleumTax: null, corporateTax: 25, marginalTake: null, headline: 25, basis: 'corporate income tax for the oil and gas industry, lower bound of the published 25–50% range', year: 2026, note: 'The rate is set contract by contract between 25% and 50%; the standard rate is 20%.', source: PWC('vietnam', CIT, 'Vietnam', '23 September 2026') },
  KZ: { country: 'Kazakhstan', regime: null, royalty: 5, petroleumTax: null, corporateTax: 20, marginalTake: null, headline: 20, basis: 'general corporate income tax', year: 2026,
    note: 'Mineral extraction tax on crude oil and condensate 5–18% depending on annual production (lower bound shown as royalty); excess profit tax is progressive, 10–60% of the net income above 25% of deductions.', source: PWC('kazakhstan', OTH, 'Kazakhstan', '23 July 2026'), source2: PWC('kazakhstan', CIT, 'Kazakhstan', '23 July 2026') },
  CO: { country: 'Colombia', regime: null, royalty: null, petroleumTax: null, corporateTax: 35, marginalTake: null, headline: 35, basis: 'general corporate income tax, before the price-dependent surtax on oil extraction', year: 2026,
    note: 'A surtax of 5% to 15% applies to oil extraction when the average price of the year is at or above 65% of the average of the preceding 120 months, so the rate can reach 50%.', source: PWC('colombia', CIT, 'Colombia', '21 July 2026') },
};

/** Fiscal terms of a country (ISO-2 code) or null. */
export const fiscalOf = (code) => FISCAL[String(code || '').toUpperCase()] || null;
/** Plain-language name of a regime code. */
export const REGIME_LABELS = { 'tax-royalty': 'tax and royalty (licence)', psc: 'production sharing contract', mixed: 'licences and production sharing contracts' };
