// Refreshes js/data/prices.js (the stand-by commodity prices bundled with a build) from the same open
// series the app reads live:  node tools/prices.mjs
import { writeFileSync } from 'node:fs';
import { join, dirname } from 'node:path';
import { fileURLToPath } from 'node:url';
const last = async (url) => { const rows = (await (await fetch(url)).text()).trim().split('\n').map((l) => l.split(',')).filter((c) => /^\d{4}-\d\d-\d\d$/.test(c[0]) && Number.isFinite(+c[1])); return rows[rows.length - 1]; };
const [b, w, g] = await Promise.all(['oil-prices/main/data/brent-daily.csv', 'oil-prices/main/data/wti-daily.csv', 'natural-gas/main/data/daily.csv'].map((p) => last('https://raw.githubusercontent.com/datasets/' + p)));
writeFileSync(join(dirname(fileURLToPath(import.meta.url)), '../js/data/prices.js'), `// Last commodity prices known when this build was made; shown (with their date) only when the live price\n// services cannot be reached. Refreshed by tools/prices.mjs.\nexport const PRICES = { date: '${b[0]}', brent: ${+b[1]}, wti: ${+w[1]}, henryHub: ${+g[1]} };\n`);
console.log('prices', b, w, g);
