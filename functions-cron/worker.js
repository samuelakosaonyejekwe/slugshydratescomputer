// Daily warm-up of the shared data cache (Cloudflare Worker with a cron trigger; Pages Functions have none).
// It only asks the public feed address of the app for each feed, so that the first visitor of the day does not wait
// for the upstream services and a failing upstream is noticed while the last good answer is still young.
// Deploy (once, and after changes):   npx wrangler deploy --config functions-cron/wrangler.toml
// Try it without waiting for the clock: open the Worker's own address — it runs the same warm-up and lists the result.
// Note: the Cache API of the feed is per Cloudflare location, so this warms the location the Worker runs in; with the
// optional KV binding FEED_STORE on the Pages project (see functions/api/feed.js) the warmed copy serves every location.
const FEEDS = ['indexes', 'cit', 'rates', 'commodities', 'spot', 'power', 'tariffs'];

async function warm(env) {
  const base = String(env.FEED_URL || 'https://hydraslug.pages.dev/api/feed');
  return Promise.all(FEEDS.map(async (id) => {
    try {
      const r = await fetch(`${base}?id=${id}`, { headers: { 'User-Agent': 'HydraSlug-feed-warmer/1' } });
      await r.arrayBuffer();
      return { id, status: r.status, cache: r.headers.get('x-feed-cache'), ageSeconds: Number(r.headers.get('x-feed-age')), stale: r.headers.get('x-feed-stale') === '1' };
    } catch (e) { return { id, status: 0, error: String(e?.message || e).slice(0, 120) }; }
  }));
}

export default {
  async scheduled(event, env, ctx) { ctx.waitUntil(warm(env).then((r) => console.log(JSON.stringify(r)))); },
  async fetch(request, env) {
    if (request.method !== 'GET') return new Response('GET only', { status: 405, headers: { Allow: 'GET' } });
    return new Response(JSON.stringify(await warm(env), null, 1), { headers: { 'Content-Type': 'application/json; charset=utf-8', 'Cache-Control': 'no-store' } });
  },
};
