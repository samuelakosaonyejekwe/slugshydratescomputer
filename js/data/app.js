// Application identity and the addresses the build is published at. Add each further host here (and to the
// connect-src list in index.html) after publishing to it with tools/deploy-mirrors.sh.
export const APP = { name: 'HydraSlug', tagline: 'Integrated multiphase flow-assurance, integrity and techno-economic simulation suite', caseExt: '.hydraslug.json' };
const PRIMARY = 'https://samuelakosaonyejekwe.github.io/slugshydratescomputer/';
/** Address the Internet Archive is asked to re-capture whenever a browser first sees a new build. */
export const ARCHIVE_SOURCE = PRIMARY + 'standalone.html';
export const MIRRORS = [
  { name: 'Primary — GitHub Pages', url: PRIMARY, note: 'Main address' },
  { name: 'Mirror — Cloudflare Pages', url: 'https://hydraslug.pages.dev/', note: 'The same build on an independent network and company: use it if the main address is unreachable. Everything works here, including live data and installation.' },
  { name: 'Independent copy — Internet Archive', kind: 'archive', url: 'https://web.archive.org/web/29991231235959id_/' + PRIMARY + 'standalone.html', note: 'The single-file edition preserved by the Internet Archive, on infrastructure unrelated to GitHub (the address resolves to the most recent capture). All seven engines, file import and decision support work there. That host blocks calls to outside data services, so site data come from the built-in atlas and the map uses the built-in coastlines.' },
];
