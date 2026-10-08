// Look-up of the bundled temperature-at-depth atlas (js/data/atlas_deep.js): seabed temperature and the
// temperature profile of the water column for any ocean location, without a connection.
import { DEEP } from './atlas_deep.js';

let bytes = null;
function grid() {
  if (bytes) return bytes;
  const bin = atob(DEEP.data); bytes = new Uint8Array(bin.length);
  for (let i = 0; i < bin.length; i++) bytes[i] = bin.charCodeAt(i);
  return bytes;
}
const cell = (k, i, j) => grid()[(k * DEEP.nlat + i) * DEEP.nlon + j];

/**
 * Temperature profile of the nearest ocean cell and the temperature at a given depth.
 * Returns { T (°C at `depth`, or at the deepest level with water in that cell), profile: { depth[], T[] }, cellKm (distance to
 * the cell centre used), levelDepth (deepest atlas level with water there), extrapolated (true when `depth` is below it) } or null far inland.
 */
export function deepTemperature(lat, lon, depth = 0) {
  const D = DEEP, i0 = Math.max(0, Math.min(D.nlat - 1, Math.round((lat - D.lat0) / D.step))), j0 = (((Math.round((lon - D.lon0) / D.step) % D.nlon) + D.nlon) % D.nlon);
  // which level do we need? the deepest atlas level not below the site depth (at least the surface)
  let need = 0; for (let k = 0; k < D.depths.length; k++) if (D.depths[k] <= depth) need = k;
  let best = null;
  for (let r = 0; r <= 4; r++) {
    for (let di = -r; di <= r; di++) for (let dj = -r; dj <= r; dj++) {
      if (Math.max(Math.abs(di), Math.abs(dj)) !== r) continue;
      const i = i0 + di, j = (((j0 + dj) % D.nlon) + D.nlon) % D.nlon;
      if (i < 0 || i >= D.nlat || !cell(0, i, j)) continue;
      let deepest = 0; for (let k = 0; k < D.depths.length; k++) if (cell(k, i, j)) deepest = k; else break;
      const la = D.lat0 + i * D.step, lo = D.lon0 + j * D.step, km = 111.2 * Math.hypot(la - lat, (((lo - lon + 540) % 360) - 180) * Math.cos((lat * Math.PI) / 180));
      // prefer cells that reach the required level, then the nearest
      const score = (deepest >= need ? 0 : 1e5 * (need - deepest)) + km;
      if (!best || score < best.score) best = { score, i, j, deepest, km };
    }
    if (best && best.deepest >= need && r >= 1) break;
  }
  if (!best) return null;
  const depths = [], T = [];
  for (let k = 0; k <= best.deepest; k++) { depths.push(D.depths[k]); T.push(+(D.offset + D.scale * cell(k, best.i, best.j)).toFixed(2)); }
  let t;
  if (depth >= depths[depths.length - 1]) t = T[T.length - 1];
  else { let k = 0; while (k < depths.length - 2 && depths[k + 1] < depth) k++; t = T[k] + ((T[k + 1] - T[k]) * (Math.max(depth, 0) - depths[k])) / (depths[k + 1] - depths[k]); }
  return { T: +t.toFixed(2), profile: { depth: depths, T }, cellKm: Math.round(best.km), levelDepth: depths[depths.length - 1], extrapolated: depth > depths[depths.length - 1] + 250, vintage: D.vintage };
}
