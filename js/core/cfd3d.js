// Three-dimensional incompressible Navier–Stokes solver in plain JavaScript (typed arrays; browser, Web Worker and Node).
//
// Grid       staggered Cartesian box (Harlow–Welch): u, v, w on the cell faces, pressure, liquid fraction and eddy viscosity at the cell
//            centres; uniform in x and z, uniform or tanh-stretched in y. Each direction is periodic, a no-slip wall or a free-slip wall.
//            A circular pipe along x is represented inside the box by a volume-fraction immersed boundary (the velocity of a face is
//            multiplied by the fluid fraction of its control volume after the predictor; Kajishima-type direct forcing).
// Advection  second-order central differences in divergence form with the volume-consistent averages of Verstappen & Veldman, which
//            conserve kinetic energy on the staggered grid (single phase); a van Leer limited upwind flux is used with the interface.
// Time       explicit: low-storage third-order Runge–Kutta of Williamson (single phase, one projection per stage) or second-order
//            Adams–Bashforth (two phases).
// Pressure   fractional-step projection; the Poisson equation is solved directly: Fourier (periodic) or cosine (walls) transforms in
//            x and z with the modified wavenumbers of the second-order stencil and a tridiagonal solve in y (or a third transform when
//            y is periodic). With two phases the variable-coefficient equation is reduced to this constant-coefficient one by the
//            pressure-splitting of Dodd & Ferrante (2014), so the direct solver is kept.
// Turbulence none (DNS when the grid resolves the flow, otherwise an under-resolved simulation), Smagorinsky with van Driest damping,
//            WALE, and the Spalart–Allmaras family: DES97, DDES and IDDES (length-scale switches of Spalart et al. 1997, 2006 and
//            Shur et al. 2008). Constants are listed in CFD3D_CONSTANTS with their sources.
// Interface  volume of fluid: direction-split THINC/WLIC with the dilatation correction of Weymouth & Yue (volume conserved to
//            round-off), variable density and viscosity, gravity, no surface tension.
//
// Everything here is pure (no DOM, no I/O) and deterministic. Units are whatever consistent set the caller uses.

/** Model constants with the open source each one was read from (see PROVENANCE of js/suites/s03_flow.js). */
export const CFD3D_CONSTANTS = Object.freeze({
  rk3: Object.freeze({ a: [0, -5 / 9, -153 / 128], b: [1 / 3, 15 / 16, 8 / 15], c: [1 / 3, 5 / 12, 1 / 4] }), // Williamson (1980) low-storage scheme; c = stage time increments
  smagorinsky: Object.freeze({ cs: 0.1, aPlus: 26 }), // channel value of the Smagorinsky constant; van Driest damping 1 − exp(−y⁺/A⁺)
  wale: Object.freeze({ cw: 0.325 }), // Nicoud & Ducros (1999); OpenFOAM WALE.H
  sa: Object.freeze({ cb1: 0.1355, cb2: 0.622, sigma: 2 / 3, kappa: 0.41, cw2: 0.3, cw3: 2, cv1: 7.1 }), // Spalart & Allmaras (1992); NASA Turbulence Modeling Resource
  des: Object.freeze({ cDes: 0.65 }), // Shur et al. (1999)
  ddes: Object.freeze({ cd1: 8, cd2: 3 }), // f_d = 1 − tanh[(8 r_d)³], Spalart et al. (2006)
  iddes: Object.freeze({ cw: 0.15, ct: 1.63, cl: 3.55, fwStar: 0.424 }), // Shur et al. (2008), Spalart–Allmaras branch
  thincBeta: 2.3,
});

// ---- transforms ---------------------------------------------------------------------------------------------------------------
const isPow2 = (n) => n > 1 && (n & (n - 1)) === 0;
// In-place radix-2 complex FFT with precomputed bit-reversal and twiddle tables (tb = { rev, cos, sin } from fftTables); inverse is unnormalised.
function fftTables(n) { const rev = new Int32Array(n), cos = new Float64Array(n >> 1), sin = new Float64Array(n >> 1); for (let i = 1, j = 0; i < n; i++) { let bit = n >> 1; for (; j & bit; bit >>= 1) j ^= bit; j ^= bit; rev[i] = j; } for (let k = 0; k < n >> 1; k++) { cos[k] = Math.cos((2 * Math.PI * k) / n); sin[k] = Math.sin((2 * Math.PI * k) / n); } return { rev, cos, sin }; }
function fftRadix2(re, im, n, inverse, tb) {
  const rev = tb.rev, cs = tb.cos, sn = tb.sin, sg = inverse ? 1 : -1;
  for (let i = 1; i < n; i++) { const j = rev[i]; if (i < j) { const tr = re[i]; re[i] = re[j]; re[j] = tr; const ti = im[i]; im[i] = im[j]; im[j] = ti; } }
  for (let len = 2, stp = n >> 1; len <= n; len <<= 1, stp >>= 1) {
    const half = len >> 1;
    for (let i = 0; i < n; i += len) for (let k = 0, q = 0; k < half; k++, q += stp) { const a = i + k, b = a + half, cr = cs[q], ci = sg * sn[q], xr = re[b] * cr - im[b] * ci, xi = re[b] * ci + im[b] * cr; re[b] = re[a] - xr; im[b] = im[a] - xi; re[a] += xr; im[a] += xi; }
  }
}
/**
 * Eigen-decomposition of the three-point second-difference operator on n equal cells of size d.
 * type 'periodic' (Fourier modes; fast transform when n is a power of two) or 'neumann' (cosine modes of the staggered grid).
 * Returns { n, lam[] (eigenvalues, ≤ 0), fwd(a, off, stride), inv(a, off, stride) } transforming a strided line in place.
 */
export function makeTransform(n, d, type = 'periodic') {
  const lam = new Float64Array(n);
  if (type === 'periodic' && isPow2(n)) {
    const re = new Float64Array(n), im = new Float64Array(n), h = n >> 1, tb = fftTables(n);
    for (let m = 0; m < n; m++) { const k = m <= h ? m : n - m; lam[m] = -(2 - 2 * Math.cos((2 * Math.PI * k) / n)) / (d * d); }
    return { n, lam, fast: true,
      fwd(a, off, s) { for (let i = 0; i < n; i++) { re[i] = a[off + i * s]; im[i] = 0; } fftRadix2(re, im, n, false, tb); const q = 1 / n; a[off] = re[0] * q; a[off + h * s] = re[h] * q; for (let k = 1; k < h; k++) { a[off + k * s] = re[k] * q; a[off + (n - k) * s] = im[k] * q; } },
      inv(a, off, s) { re[0] = a[off]; im[0] = 0; re[h] = a[off + h * s]; im[h] = 0; for (let k = 1; k < h; k++) { const r = a[off + k * s], q = a[off + (n - k) * s]; re[k] = r; im[k] = q; re[n - k] = r; im[n - k] = -q; } fftRadix2(re, im, n, true, tb); for (let i = 0; i < n; i++) a[off + i * s] = re[i]; } };
  }
  const Q = new Float64Array(n * n), t = new Float64Array(n);
  if (type === 'periodic') {
    let m = 0; const put = (f, ev) => { for (let i = 0; i < n; i++) Q[m * n + i] = f(i); lam[m++] = ev; };
    put(() => 1 / Math.sqrt(n), 0);
    for (let k = 1; 2 * k < n; k++) { const ev = -(2 - 2 * Math.cos((2 * Math.PI * k) / n)) / (d * d), c = Math.sqrt(2 / n); put((i) => c * Math.cos((2 * Math.PI * k * i) / n), ev); put((i) => c * Math.sin((2 * Math.PI * k * i) / n), ev); }
    if (n % 2 === 0) put((i) => (i % 2 ? -1 : 1) / Math.sqrt(n), -4 / (d * d));
  } else for (let k = 0; k < n; k++) { const c = Math.sqrt((k ? 2 : 1) / n); lam[k] = -(2 - 2 * Math.cos((Math.PI * k) / n)) / (d * d); for (let i = 0; i < n; i++) Q[k * n + i] = c * Math.cos((Math.PI * k * (i + 0.5)) / n); }
  return { n, lam, fast: false,
    fwd(a, off, s) { for (let k = 0; k < n; k++) { let x = 0; const r = k * n; for (let i = 0; i < n; i++) x += Q[r + i] * a[off + i * s]; t[k] = x; } for (let k = 0; k < n; k++) a[off + k * s] = t[k]; },
    inv(a, off, s) { for (let i = 0; i < n; i++) { let x = 0; for (let k = 0; k < n; k++) x += Q[k * n + i] * a[off + k * s]; t[i] = x; } for (let i = 0; i < n; i++) a[off + i * s] = t[i]; } };
}

/**
 * Direct solver of the discrete Poisson equation ∇²p = r on the staggered box with homogeneous Neumann conditions on walls.
 * g: { nx, ny, nz, dx, dz, dy[1..ny], dyc[0..ny], perX, perY, perZ }. Returns { solve(W) } where W (Float64Array nx·ny·nz, index
 * i + nx (j + ny l), zero-based) holds r on entry and p on return (mean of the constant mode removed by pinning one value).
 */
export function makePoisson(g) {
  const { nx, ny, nz } = g, Tx = makeTransform(nx, g.dx, g.perX ? 'periodic' : 'neumann'), Tz = makeTransform(nz, g.dz, g.perZ ? 'periodic' : 'neumann');
  const Ty = g.perY ? makeTransform(ny, g.dy[1], 'periodic') : null, sxy = nx * ny;
  const a = new Float64Array(ny), c = new Float64Array(ny), cp = new Float64Array(ny), dp = new Float64Array(ny);
  if (!Ty) for (let j = 0; j < ny; j++) { a[j] = j > 0 ? 1 / (g.dy[j + 1] * g.dyc[j]) : 0; c[j] = j < ny - 1 ? 1 / (g.dy[j + 1] * g.dyc[j + 1]) : 0; }
  return { fast: Tx.fast && Tz.fast,
    solve(W) {
      for (let l = 0; l < nz; l++) for (let j = 0; j < ny; j++) Tx.fwd(W, nx * (j + ny * l), 1);
      for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) Tz.fwd(W, i + nx * j, sxy);
      if (Ty) {
        for (let l = 0; l < nz; l++) for (let i = 0; i < nx; i++) { const off = i + sxy * l; Ty.fwd(W, off, nx); for (let j = 0; j < ny; j++) { const e = Tx.lam[i] + Ty.lam[j] + Tz.lam[l]; W[off + j * nx] = e < -1e-14 ? W[off + j * nx] / e : 0; } Ty.inv(W, off, nx); }
      } else {
        for (let l = 0; l < nz; l++) for (let i = 0; i < nx; i++) {
          const e = Tx.lam[i] + Tz.lam[l], off = i + sxy * l, sing = e > -1e-14;
          let b = -(a[0] + c[0]) + e; if (sing) b *= 2; // pins the first value of the constant mode
          cp[0] = c[0] / b; dp[0] = W[off] / b;
          for (let j = 1; j < ny; j++) { const m = -(a[j] + c[j]) + e - a[j] * cp[j - 1]; cp[j] = c[j] / m; dp[j] = (W[off + j * nx] - a[j] * dp[j - 1]) / m; }
          W[off + (ny - 1) * nx] = dp[ny - 1];
          for (let j = ny - 2; j >= 0; j--) W[off + j * nx] = dp[j] - cp[j] * W[off + (j + 1) * nx];
        }
      }
      for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) Tz.inv(W, i + nx * j, sxy);
      for (let l = 0; l < nz; l++) for (let j = 0; j < ny; j++) Tx.inv(W, nx * (j + ny * l), 1);
    } };
}

// ---- volume of fluid ------------------------------------------------------------------------------------------------------------
// THINC flux through the downstream face of a donor cell (fraction of the donor volume) for a Courant number c ≥ 0; pm, p, pp: upstream, donor, downstream.
function thincFlux(pm, p, pp, c, beta) {
  let fl = c * p;
  if (p > 1e-8 && p < 1 - 1e-8 && (pp - p) * (p - pm) > 0) { const gm = pp > pm ? 1 : -1, qq = (beta * (2 * p - 1)) / gm, w = (Math.exp(beta) - Math.exp(qq)) / (Math.exp(qq) - Math.exp(-beta)); if (w > 0 && Number.isFinite(w)) { const xt = Math.log(w) / (2 * beta); fl = 0.5 * (c + (gm / beta) * (Math.log(Math.cosh(beta * (1 - xt))) - Math.log(Math.cosh(beta * (1 - c - xt))))); } }
  const lo = Math.max(0, c - (1 - p)), hi = Math.min(c, p);
  return fl < lo ? lo : fl > hi ? hi : fl;
}
const vanLeer = (qU, qC, qD) => { const a = qC - qU, b = qD - qC; return a * b > 0 ? qC + (a * b) / (a + b) : qC; }; // limited face value from the upwind side (harmonic limiter)

/**
 * Create a flow solver. Options (all optional):
 *  nx, ny, nz, lx, ly, lz        cells and box size; stretchY > 0 clusters the y grid at both walls (tanh, factor ≈ 1.5–2.5)
 *  bc: { x, y, z }               'periodic' | 'wall' (no slip) | 'slip'; default x and z periodic, y walls
 *  nu                            kinematic viscosity (single phase)
 *  force: [fx, fy, fz]           body acceleration; bulk: target mean x-velocity (the uniform forcing is then adjusted every step)
 *  sgs                           'none' | 'smagorinsky' | 'wale' | 'des' | 'ddes' | 'iddes'
 *  scheme                        'rk3' (default, single phase) | 'ab2'; cfl (default 0.8 / 0.2); dtMax; advection 'central' | 'vanleer'
 *  wallOrder                     3 (default): wall shear from the parabola through the wall and two points (plane Poiseuille flow exact); 2: ghost cell
 *  pipe: { radius, cy, cz }      circular pipe along x by immersed boundary (bc.y and bc.z should be walls or slip)
 *  vof: { rhoL, rhoG, muL, muG, gravity: [gx, gy, gz], c0(x, y, z) → liquid fraction, pIter }   two-phase flow
 *  init(x, y, z) → [u, v, w]     initial velocity (evaluated at each face); nuTilde0(d) initial Spalart–Allmaras variable
 * Returns the solver object: { grid, t, steps, u, v, w, p, c, nut, step(), stats(), profile(), slice(), … } (see the members below).
 */
export function createFlow3D(o = {}) {
  const nx = Math.max(2, Math.round(o.nx ?? 16)), ny = Math.max(2, Math.round(o.ny ?? 16)), nz = Math.max(2, Math.round(o.nz ?? 16)), lx = o.lx ?? 1, ly = o.ly ?? 1, lz = o.lz ?? 1;
  const bc = { x: 'periodic', y: 'wall', z: 'periodic', ...(o.bc || {}) }, perX = bc.x === 'periodic', perY = bc.y === 'periodic', perZ = bc.z === 'periodic';
  const sgnX = bc.x === 'slip' ? 1 : -1, sgnY = bc.y === 'slip' ? 1 : -1, sgnZ = bc.z === 'slip' ? 1 : -1;
  const sx = nx + 2, sy = ny + 2, sz = nz + 2, N = sx * sy * sz, Y = sx, Z = sx * sy, dx = lx / nx, dz = lz / nz, rdx = 1 / dx, rdz = 1 / dz;
  const K = CFD3D_CONSTANTS, vof = o.vof || null, sgs = vof && !o.sgs ? 'none' : o.sgs || 'none', saFamily = sgs === 'des' || sgs === 'ddes' || sgs === 'iddes';
  const scheme = vof ? 'ab2' : o.scheme === 'ab2' ? 'ab2' : 'rk3', cfl = o.cfl ?? (scheme === 'rk3' ? 0.8 : 0.2), limited = (o.advection || (vof ? 'vanleer' : 'central')) === 'vanleer', wall3 = (o.wallOrder ?? 3) === 3 && bc.y === 'wall' && ny >= 3;
  // ---- y grid
  const yf = new Float64Array(ny + 1), yc = new Float64Array(ny + 2), dy = new Float64Array(ny + 2), dyc = new Float64Array(ny + 1), st = perY ? 0 : o.stretchY || 0;
  for (let j = 0; j <= ny; j++) yf[j] = st > 0 ? 0.5 * ly * (1 + Math.tanh(st * ((2 * j) / ny - 1)) / Math.tanh(st)) : (ly * j) / ny;
  for (let j = 1; j <= ny; j++) { yc[j] = 0.5 * (yf[j] + yf[j - 1]); dy[j] = yf[j] - yf[j - 1]; }
  dy[0] = perY ? dy[ny] : dy[1]; dy[ny + 1] = perY ? dy[1] : dy[ny]; yc[0] = yf[0] - 0.5 * dy[0]; yc[ny + 1] = yf[ny] + 0.5 * dy[ny + 1];
  for (let j = 0; j <= ny; j++) dyc[j] = yc[j + 1] - yc[j];
  const xc = new Float64Array(nx + 2), zc = new Float64Array(nz + 2); for (let i = 0; i <= nx + 1; i++) xc[i] = (i - 0.5) * dx; for (let l = 0; l <= nz + 1; l++) zc[l] = (l - 0.5) * dz;
  const dyMin = Math.min(...Array.from(dy.subarray(1, ny + 1)));
  // one-sided wall-gradient weights (exact for a parabola through the wall)
  const y1s = yc[1] - yf[0], y2s = yc[2] - yf[0], gs1 = y2s / (y1s * (y2s - y1s)), gs2 = -y1s / (y2s * (y2s - y1s)), y1n = yf[ny] - yc[ny], y2n = yf[ny] - yc[ny - 1], gn1 = y2n / (y1n * (y2n - y1n)), gn2 = -y1n / (y2n * (y2n - y1n));
  // ---- fields
  const u = new Float64Array(N), v = new Float64Array(N), w = new Float64Array(N), p = new Float64Array(N), qu = new Float64Array(N), qv = new Float64Array(N), qw = new Float64Array(N), ru = new Float64Array(N), rv = new Float64Array(N), rw = new Float64Array(N);
  const nut = sgs !== 'none' ? new Float64Array(N) : null, mu = new Float64Array(N), W = new Float64Array(nx * ny * nz), poisson = makePoisson({ nx, ny, nz, dx, dz, dy, dyc, perX, perY, perZ });
  const C = vof ? new Float64Array(N) : null, rho = vof ? new Float64Array(N) : null, pOld = vof ? new Float64Array(N) : null, pHat = vof ? new Float64Array(N) : null, flux = vof ? new Float64Array(N) : null, cc = vof ? new Uint8Array(N) : null;
  const nuT = saFamily ? new Float64Array(N) : null, nuT2 = saFamily ? new Float64Array(N) : null, dWall = new Float64Array(N), lDes = saFamily ? new Float64Array(N) : null;
  const nu0 = vof ? 0 : o.nu ?? 1e-3, rho0 = vof ? Math.min(vof.rhoL, vof.rhoG) : 1, grav = vof ? vof.gravity || [0, -9.80665, 0] : [0, 0, 0], force = (o.force || [0, 0, 0]).slice();
  const idx = (i, j, l) => i + sx * (j + sy * l);
  // ---- immersed pipe
  const pipe = o.pipe || null; let fu = null, fv = null, fw = null, fcell = null, fluidVolume = lx * ly * lz;
  if (pipe) {
    const R = pipe.radius, cy = pipe.cy ?? ly / 2, cz = pipe.cz ?? lz / 2, ns = pipe.samples || 6, frac = (ya, yb, za, zb) => { let s = 0; for (let a = 0; a < ns; a++) for (let b = 0; b < ns; b++) { const yy = ya + ((a + 0.5) * (yb - ya)) / ns - cy, zz = za + ((b + 0.5) * (zb - za)) / ns - cz; if (yy * yy + zz * zz < R * R) s++; } return s / (ns * ns); };
    fu = new Float64Array(N); fv = new Float64Array(N); fw = new Float64Array(N); fcell = new Float64Array(N); fluidVolume = 0;
    for (let l = 0; l <= nz + 1; l++) for (let j = 0; j <= ny + 1; j++) {
      const ya = yc[j] - 0.5 * dy[j], yb = yc[j] + 0.5 * dy[j], za = zc[l] - 0.5 * dz, zb = zc[l] + 0.5 * dz, fc = frac(ya, yb, za, zb), fvv = j <= ny ? frac(yc[j], yc[Math.min(j + 1, ny + 1)], za, zb) : 0, fww = frac(ya, yb, zc[l], zc[l] + dz);
      for (let i = 0; i <= nx + 1; i++) { const k = idx(i, j, l); fu[k] = fc; fcell[k] = fc; fv[k] = fvv; fw[k] = fww; if (i >= 1 && i <= nx && j >= 1 && j <= ny && l >= 1 && l <= nz) fluidVolume += fc * dx * dy[j] * dz; dWall[k] = Math.max(R - Math.hypot(yc[j] - cy, zc[l] - cz), 1e-9 * R); }
    }
  } else {
    const useY = bc.y === 'wall', useZ = bc.z === 'wall', useX = bc.x === 'wall';
    for (let l = 0; l <= nz + 1; l++) for (let j = 0; j <= ny + 1; j++) for (let i = 0; i <= nx + 1; i++) { let d = Infinity; if (useY) d = Math.min(d, yc[j] - yf[0], yf[ny] - yc[j]); if (useZ) d = Math.min(d, zc[l], lz - zc[l]); if (useX) d = Math.min(d, xc[i], lx - xc[i]); dWall[idx(i, j, l)] = Number.isFinite(d) ? Math.max(d, 1e-12) : 1e30; }
    if (typeof o.wallDistance === 'function') for (let l = 1; l <= nz; l++) for (let j = 1; j <= ny; j++) for (let i = 1; i <= nx; i++) dWall[idx(i, j, l)] = Math.max(o.wallDistance(xc[i], yc[j], zc[l]), 1e-12);
  }
  // ---- ghost cells
  const nxU = perX ? nx : nx - 1, nyV = perY ? ny : ny - 1, nzW = perZ ? nz : nz - 1;
  // kind: 0 scalar (zero gradient at walls), 1/2/3 velocity component normal to x/y/z, 4 scalar that vanishes at no-slip walls
  function fill(f, kind) {
    for (let l = 0; l < sz; l++) for (let j = 0; j < sy; j++) { const r = sx * (j + sy * l);
      if (perX) { f[r] = f[r + nx]; f[r + nx + 1] = f[r + 1]; } else if (kind === 1) { f[r] = 0; f[r + nx] = 0; f[r + nx + 1] = 0; } else { const s = kind === 0 ? 1 : kind === 4 ? (bc.x === 'wall' ? -1 : 1) : sgnX; f[r] = s * f[r + 1]; f[r + nx + 1] = s * f[r + nx]; } }
    for (let l = 0; l < sz; l++) for (let i = 0; i < sx; i++) { const r = i + Z * l;
      if (perY) { f[r] = f[r + ny * Y]; f[r + (ny + 1) * Y] = f[r + Y]; } else if (kind === 2) { f[r] = 0; f[r + ny * Y] = 0; f[r + (ny + 1) * Y] = 0; } else { const s = kind === 0 ? 1 : kind === 4 ? (bc.y === 'wall' ? -1 : 1) : sgnY; f[r] = s * f[r + Y]; f[r + (ny + 1) * Y] = s * f[r + ny * Y]; } }
    for (let j = 0; j < sy; j++) for (let i = 0; i < sx; i++) { const r = i + Y * j;
      if (perZ) { f[r] = f[r + nz * Z]; f[r + (nz + 1) * Z] = f[r + Z]; } else if (kind === 3) { f[r] = 0; f[r + nz * Z] = 0; f[r + (nz + 1) * Z] = 0; } else { const s = kind === 0 ? 1 : kind === 4 ? (bc.z === 'wall' ? -1 : 1) : sgnZ; f[r] = s * f[r + Z]; f[r + (nz + 1) * Z] = s * f[r + nz * Z]; } }
  }
  const fillVel = () => { fill(u, 1); fill(v, 2); fill(w, 3); };
  // ---- initial state
  const S = { t: 0, steps: 0, dt: 0, forcing: 0, forcingMean: 0, forcingTime: 0, uTau: 0, poissonSolves: 0, maxDiv: 0 };
  if (typeof o.init === 'function') for (let l = 1; l <= nz; l++) for (let j = 1; j <= ny; j++) for (let i = 1; i <= nx; i++) { const k = idx(i, j, l); if (i <= nxU) u[k] = o.init(xc[i] + 0.5 * dx, yc[j], zc[l])[0]; if (j <= nyV) v[k] = o.init(xc[i], yf[j], zc[l])[1]; if (l <= nzW) w[k] = o.init(xc[i], yc[j], zc[l] + 0.5 * dz)[2]; }
  if (vof) { for (let l = 1; l <= nz; l++) for (let j = 1; j <= ny; j++) for (let i = 1; i <= nx; i++) C[idx(i, j, l)] = Math.min(1, Math.max(0, typeof vof.c0 === 'function' ? vof.c0(xc[i], yc[j], zc[l], dx, dy[j], dz) : 0)); fill(C, 0); }
  if (saFamily) { const f0 = typeof o.nuTilde0 === 'function' ? o.nuTilde0 : () => 3 * nu0; for (let k = 0; k < N; k++) nuT[k] = Math.max(f0(dWall[k]), 0); fill(nuT, 4); }
  const cellVolume = (j) => dx * dy[j] * dz;
  function properties() { // density and dynamic viscosity (two phases) or kinematic viscosity (single phase) at the cell centres
    if (vof) { for (let k = 0; k < N; k++) { const c = C[k]; rho[k] = vof.rhoG + (vof.rhoL - vof.rhoG) * c; mu[k] = vof.muG + (vof.muL - vof.muG) * c + (nut ? rho[k] * nut[k] : 0); } }
    else if (nut) for (let k = 0; k < N; k++) mu[k] = nu0 + nut[k];
  }
  // ---- sub-grid viscosity
  function eddyViscosity(dt) {
    if (!nut) return;
    const cs = o.cs ?? K.smagorinsky.cs, aPlus = K.smagorinsky.aPlus, cw = o.cw ?? K.wale.cw, sa = K.sa, cw1 = sa.cb1 / (sa.kappa * sa.kappa) + (1 + sa.cb2) / sa.sigma, cDes = o.cDes ?? K.des.cDes, k2 = sa.kappa * sa.kappa, uTau = S.uTau, hx = dx, hz = dz;
    for (let l = 1; l <= nz; l++) for (let j = 1; j <= ny; j++) {
      const rdyj = 1 / dy[j], rdy2 = 1 / (yc[j + 1] - yc[j - 1]), delta = Math.cbrt(dx * dy[j] * dz), hMax = Math.max(dx, dy[j], dz);
      for (let i = 1; i <= nx; i++) {
        const k = i + sx * (j + sy * l);
        const g11 = (u[k] - u[k - 1]) * rdx, g22 = (v[k] - v[k - Y]) * rdyj, g33 = (w[k] - w[k - Z]) * rdz;
        const g12 = 0.5 * (u[k + Y] + u[k + Y - 1] - u[k - Y] - u[k - Y - 1]) * rdy2, g13 = 0.25 * (u[k + Z] + u[k + Z - 1] - u[k - Z] - u[k - Z - 1]) * rdz;
        const g21 = 0.25 * (v[k + 1] + v[k + 1 - Y] - v[k - 1] - v[k - 1 - Y]) * rdx, g23 = 0.25 * (v[k + Z] + v[k + Z - Y] - v[k - Z] - v[k - Z - Y]) * rdz;
        const g31 = 0.25 * (w[k + 1] + w[k + 1 - Z] - w[k - 1] - w[k - 1 - Z]) * rdx, g32 = 0.5 * (w[k + Y] + w[k + Y - Z] - w[k - Y] - w[k - Y - Z]) * rdy2;
        const s12 = 0.5 * (g12 + g21), s13 = 0.5 * (g13 + g31), s23 = 0.5 * (g23 + g32), ss = g11 * g11 + g22 * g22 + g33 * g33 + 2 * (s12 * s12 + s13 * s13 + s23 * s23), d = dWall[k];
        if (sgs === 'smagorinsky') { const damp = uTau > 0 && d < 1e29 ? 1 - Math.exp((-d * uTau) / (nu0 * aPlus)) : 1, ls = cs * delta * damp; nut[k] = ls * ls * Math.sqrt(2 * ss); }
        else if (sgs === 'wale') {
          const a11 = g11 * g11 + g12 * g21 + g13 * g31, a12 = g11 * g12 + g12 * g22 + g13 * g32, a13 = g11 * g13 + g12 * g23 + g13 * g33, a21 = g21 * g11 + g22 * g21 + g23 * g31, a22 = g21 * g12 + g22 * g22 + g23 * g32, a23 = g21 * g13 + g22 * g23 + g23 * g33, a31 = g31 * g11 + g32 * g21 + g33 * g31, a32 = g31 * g12 + g32 * g22 + g33 * g32, a33 = g31 * g13 + g32 * g23 + g33 * g33;
          const tr = (a11 + a22 + a33) / 3, d11 = a11 - tr, d22 = a22 - tr, d33 = a33 - tr, d12 = 0.5 * (a12 + a21), d13 = 0.5 * (a13 + a31), d23 = 0.5 * (a23 + a32), sd = d11 * d11 + d22 * d22 + d33 * d33 + 2 * (d12 * d12 + d13 * d13 + d23 * d23);
          const den = ss ** 2.5 + sd ** 1.25; nut[k] = den > 1e-300 ? (cw * delta) ** 2 * (sd ** 1.5 / den) : 0;
        } else { // Spalart–Allmaras family: length scale first
          const gg = g11 * g11 + g22 * g22 + g33 * g33 + g12 * g12 + g13 * g13 + g21 * g21 + g23 * g23 + g31 * g31 + g32 * g32, sMag = Math.max(Math.sqrt(gg), 1e-10);
          const om = Math.sqrt((g12 - g21) ** 2 + (g13 - g31) ** 2 + (g23 - g32) ** 2), nt = nuT[k], chi = nt / nu0, c3 = chi * chi * chi, fv1 = c3 / (c3 + sa.cv1 ** 3), fv2 = 1 - chi / (1 + chi * fv1);
          let len;
          if (sgs === 'des') len = Math.min(d, cDes * hMax);
          else if (sgs === 'ddes') { const rd = (nt * fv1 + nu0) / (sMag * k2 * d * d), fd = 1 - Math.tanh((K.ddes.cd1 * rd) ** K.ddes.cd2); len = d - fd * Math.max(0, d - cDes * hMax); }
          else {
            const I = K.iddes, del = Math.min(Math.max(I.cw * d, I.cw * hMax, dy[j]), hMax), alpha = 0.25 - d / hMax, fB = Math.min(2 * Math.exp(-9 * alpha * alpha), 1), rdt = (nt * fv1) / (sMag * k2 * d * d), rdl = nu0 / (sMag * k2 * d * d);
            const fdt = 1 - Math.tanh((8 * rdt) ** 3), fdTilde = Math.max(1 - fdt, fB), fe1 = alpha >= 0 ? 2 * Math.exp(-11.09 * alpha * alpha) : 2 * Math.exp(-9 * alpha * alpha), ft = Math.tanh((I.ct * I.ct * rdt) ** 3), fl = Math.tanh((I.cl * I.cl * rdl) ** 10);
            const psi = Math.sqrt(Math.min(100, (1 - (sa.cb1 / (cw1 * k2 * I.fwStar)) * fv2) / Math.max(fv1, 1e-10))), fe = Math.max(fe1 - 1, 0) * psi * (1 - Math.max(ft, fl));
            len = fdTilde * (1 + fe) * d + (1 - fdTilde) * cDes * psi * del;
          }
          lDes[k] = len;
          // transport of the working variable (explicit Euler, first-order upwind advection, destruction treated implicitly)
          const uc = 0.5 * (u[k] + u[k - 1]), vc = 0.5 * (v[k] + v[k - Y]), wc = 0.5 * (w[k] + w[k - Z]), dn = dyc[j], ds = dyc[j - 1];
          const adv = (uc > 0 ? uc * (nt - nuT[k - 1]) : uc * (nuT[k + 1] - nt)) * rdx + (vc > 0 ? (vc * (nt - nuT[k - Y])) / ds : (vc * (nuT[k + Y] - nt)) / dn) + (wc > 0 ? wc * (nt - nuT[k - Z]) : wc * (nuT[k + Z] - nt)) * rdz;
          const de = nu0 + 0.5 * (nt + nuT[k + 1]), dw_ = nu0 + 0.5 * (nt + nuT[k - 1]), dnn = nu0 + 0.5 * (nt + nuT[k + Y]), dss = nu0 + 0.5 * (nt + nuT[k - Y]), dt_ = nu0 + 0.5 * (nt + nuT[k + Z]), db = nu0 + 0.5 * (nt + nuT[k - Z]);
          const diff = (de * (nuT[k + 1] - nt) - dw_ * (nt - nuT[k - 1])) * rdx * rdx + ((dnn * (nuT[k + Y] - nt)) / dn - (dss * (nt - nuT[k - Y])) / ds) * rdyj + (dt_ * (nuT[k + Z] - nt) - db * (nt - nuT[k - Z])) * rdz * rdz;
          const gx = 0.5 * (nuT[k + 1] - nuT[k - 1]) * rdx, gy = (nuT[k + Y] - nuT[k - Y]) * rdy2, gz = 0.5 * (nuT[k + Z] - nuT[k - Z]) * rdz;
          const sT = Math.max(om + (nt * fv2) / (k2 * len * len), 0.3 * om), r = Math.min(nt / (Math.max(sT, 1e-12) * k2 * len * len), 10), g = r + sa.cw2 * (r ** 6 - r), fwv = g * ((1 + sa.cw3 ** 6) / (g ** 6 + sa.cw3 ** 6)) ** (1 / 6);
          const src = -adv + sa.cb1 * sT * nt + (diff + sa.cb2 * (gx * gx + gy * gy + gz * gz)) / sa.sigma;
          nuT2[k] = Math.max((nt + dt * src) / (1 + (dt * cw1 * fwv * nt) / (len * len)), 0) * (fcell ? fcell[k] : 1);
        }
      }
    }
    if (saFamily) { if (dt > 0) nuT.set(nuT2); fill(nuT, 4); for (let k = 0; k < N; k++) { const c = nuT[k] / nu0, c3 = c * c * c; nut[k] = nuT[k] > 0 ? (nuT[k] * c3) / (c3 + sa.cv1 ** 3) : 0; } }
    fill(nut, 4);
  }
  // ---- right-hand sides (advection, viscous stress, body force); pressure is added by the projection
  const variableMu = !!nut || !!vof;
  function rhs() {
    const fx = force[0] + S.forcing + grav[0], fy = force[1] + grav[1], fz = force[2] + grav[2], nu = nu0;
    for (let l = 1; l <= nz; l++) for (let j = 1; j <= ny; j++) {
      const rdyj = 1 / dy[j], rn = 1 / dyc[j], rs = 1 / dyc[j - 1], wj = dy[j] / (dy[j] + dy[j + 1]), wj1 = 1 - wj, southWall = wall3 && j === 1, northWall = wall3 && j === ny, row = sx * (j + sy * l);
      const rdyn = j < ny || perY ? 1 / dy[j + 1] : 0;
      for (let i = 1; i <= nx; i++) {
        const k = row + i, uk = u[k], vk = v[k], wk = w[k];
        // ---------------- u at (i+½, j, l)
        if (i <= nxU) {
          const ue = 0.5 * (uk + u[k + 1]), uw = 0.5 * (uk + u[k - 1]), vn = 0.5 * (vk + v[k + 1]), vs = 0.5 * (v[k - Y] + v[k - Y + 1]), wt = 0.5 * (wk + w[k + 1]), wb = 0.5 * (w[k - Z] + w[k - Z + 1]);
          let qe, qw_, qn, qs, qt, qb;
          if (limited) {
            qe = ue >= 0 ? vanLeer(u[k - 1], uk, u[k + 1]) : vanLeer(i + 2 <= nx + 1 ? u[k + 2] : u[k + 1], u[k + 1], uk); qw_ = uw >= 0 ? vanLeer(i - 2 >= 0 ? u[k - 2] : u[k - 1], u[k - 1], uk) : vanLeer(u[k + 1], uk, u[k - 1]);
            qn = vn >= 0 ? vanLeer(u[k - Y], uk, u[k + Y]) : vanLeer(j + 2 <= ny + 1 ? u[k + 2 * Y] : u[k + Y], u[k + Y], uk); qs = vs >= 0 ? vanLeer(j - 2 >= 0 ? u[k - 2 * Y] : u[k - Y], u[k - Y], uk) : vanLeer(u[k + Y], uk, u[k - Y]);
            qt = wt >= 0 ? vanLeer(u[k - Z], uk, u[k + Z]) : vanLeer(l + 2 <= nz + 1 ? u[k + 2 * Z] : u[k + Z], u[k + Z], uk); qb = wb >= 0 ? vanLeer(l - 2 >= 0 ? u[k - 2 * Z] : u[k - Z], u[k - Z], uk) : vanLeer(u[k + Z], uk, u[k - Z]);
          } else { qe = ue; qw_ = uw; qn = 0.5 * (uk + u[k + Y]); qs = 0.5 * (uk + u[k - Y]); qt = 0.5 * (uk + u[k + Z]); qb = 0.5 * (uk + u[k - Z]); }
          let r = -((ue * qe - uw * qw_) * rdx + (vn * qn - vs * qs) * rdyj + (wt * qt - wb * qb) * rdz), gN = (u[k + Y] - uk) * rn, gS = (uk - u[k - Y]) * rs;
          if (southWall) gS = gs1 * uk + gs2 * u[k + Y]; if (northWall) gN = -(gn1 * uk + gn2 * u[k - Y]);
          if (variableMu) {
            const mN = 0.25 * (mu[k] + mu[k + 1] + mu[k + Y] + mu[k + Y + 1]), mS = 0.25 * (mu[k] + mu[k + 1] + mu[k - Y] + mu[k - Y + 1]), mT = 0.25 * (mu[k] + mu[k + 1] + mu[k + Z] + mu[k + Z + 1]), mB = 0.25 * (mu[k] + mu[k + 1] + mu[k - Z] + mu[k - Z + 1]);
            const vis = (2 * mu[k + 1] * (u[k + 1] - uk) - 2 * mu[k] * (uk - u[k - 1])) * rdx * rdx + (mN * (gN + (v[k + 1] - vk) * rdx) - mS * (gS + (v[k + 1 - Y] - v[k - Y]) * rdx)) * rdyj + (mT * ((u[k + Z] - uk) * rdz + (w[k + 1] - wk) * rdx) - mB * ((uk - u[k - Z]) * rdz + (w[k + 1 - Z] - w[k - Z]) * rdx)) * rdz;
            r += vof ? (2 * vis) / (rho[k] + rho[k + 1]) : vis;
          } else r += nu * ((u[k + 1] - 2 * uk + u[k - 1]) * rdx * rdx + (gN - gS) * rdyj + (u[k + Z] - 2 * uk + u[k - Z]) * rdz * rdz);
          ru[k] = r + fx;
        }
        // ---------------- v at (i, j+½, l)
        if (j <= nyV) {
          const ue = wj * uk + wj1 * u[k + Y], uw = wj * u[k - 1] + wj1 * u[k - 1 + Y], vn = 0.5 * (vk + v[k + Y]), vs = 0.5 * (vk + v[k - Y]), wt = wj * wk + wj1 * w[k + Y], wb = wj * w[k - Z] + wj1 * w[k - Z + Y];
          let qe, qw_, qn, qs, qt, qb;
          if (limited) {
            qe = ue >= 0 ? vanLeer(v[k - 1], vk, v[k + 1]) : vanLeer(i + 2 <= nx + 1 ? v[k + 2] : v[k + 1], v[k + 1], vk); qw_ = uw >= 0 ? vanLeer(i - 2 >= 0 ? v[k - 2] : v[k - 1], v[k - 1], vk) : vanLeer(v[k + 1], vk, v[k - 1]);
            qn = vn >= 0 ? vanLeer(v[k - Y], vk, v[k + Y]) : vanLeer(j + 2 <= ny + 1 ? v[k + 2 * Y] : v[k + Y], v[k + Y], vk); qs = vs >= 0 ? vanLeer(j - 2 >= 0 ? v[k - 2 * Y] : v[k - Y], v[k - Y], vk) : vanLeer(v[k + Y], vk, v[k - Y]);
            qt = wt >= 0 ? vanLeer(v[k - Z], vk, v[k + Z]) : vanLeer(l + 2 <= nz + 1 ? v[k + 2 * Z] : v[k + Z], v[k + Z], vk); qb = wb >= 0 ? vanLeer(l - 2 >= 0 ? v[k - 2 * Z] : v[k - Z], v[k - Z], vk) : vanLeer(v[k + Z], vk, v[k - Z]);
          } else { qe = 0.5 * (vk + v[k + 1]); qw_ = 0.5 * (vk + v[k - 1]); qn = vn; qs = vs; qt = 0.5 * (vk + v[k + Z]); qb = 0.5 * (vk + v[k - Z]); }
          let r = -((ue * qe - uw * qw_) * rdx + (vn * qn - vs * qs) * rn + (wt * qt - wb * qb) * rdz);
          if (variableMu) {
            const mE = 0.25 * (mu[k] + mu[k + 1] + mu[k + Y] + mu[k + Y + 1]), mW = 0.25 * (mu[k] + mu[k - 1] + mu[k + Y] + mu[k + Y - 1]), mT = 0.25 * (mu[k] + mu[k + Z] + mu[k + Y] + mu[k + Y + Z]), mB = 0.25 * (mu[k] + mu[k - Z] + mu[k + Y] + mu[k + Y - Z]);
            const vis = (mE * ((u[k + Y] - uk) * rn + (v[k + 1] - vk) * rdx) - mW * ((u[k - 1 + Y] - u[k - 1]) * rn + (vk - v[k - 1]) * rdx)) * rdx + (2 * mu[k + Y] * (v[k + Y] - vk) * rdyn - 2 * mu[k] * (vk - v[k - Y]) * rdyj) * rn + (mT * ((v[k + Z] - vk) * rdz + (w[k + Y] - wk) * rn) - mB * ((vk - v[k - Z]) * rdz + (w[k + Y - Z] - w[k - Z]) * rn)) * rdz;
            r += vof ? (2 * vis) / (rho[k] + rho[k + Y]) : vis;
          } else r += nu * ((v[k + 1] - 2 * vk + v[k - 1]) * rdx * rdx + ((v[k + Y] - vk) * rdyn - (vk - v[k - Y]) * rdyj) * rn + (v[k + Z] - 2 * vk + v[k - Z]) * rdz * rdz);
          rv[k] = r + fy;
        }
        // ---------------- w at (i, j, l+½)
        if (l <= nzW) {
          const ue = 0.5 * (uk + u[k + Z]), uw = 0.5 * (u[k - 1] + u[k - 1 + Z]), vn = 0.5 * (vk + v[k + Z]), vs = 0.5 * (v[k - Y] + v[k - Y + Z]), wt = 0.5 * (wk + w[k + Z]), wb = 0.5 * (wk + w[k - Z]);
          let qe, qw_, qn, qs, qt, qb;
          if (limited) {
            qe = ue >= 0 ? vanLeer(w[k - 1], wk, w[k + 1]) : vanLeer(i + 2 <= nx + 1 ? w[k + 2] : w[k + 1], w[k + 1], wk); qw_ = uw >= 0 ? vanLeer(i - 2 >= 0 ? w[k - 2] : w[k - 1], w[k - 1], wk) : vanLeer(w[k + 1], wk, w[k - 1]);
            qn = vn >= 0 ? vanLeer(w[k - Y], wk, w[k + Y]) : vanLeer(j + 2 <= ny + 1 ? w[k + 2 * Y] : w[k + Y], w[k + Y], wk); qs = vs >= 0 ? vanLeer(j - 2 >= 0 ? w[k - 2 * Y] : w[k - Y], w[k - Y], wk) : vanLeer(w[k + Y], wk, w[k - Y]);
            qt = wt >= 0 ? vanLeer(w[k - Z], wk, w[k + Z]) : vanLeer(l + 2 <= nz + 1 ? w[k + 2 * Z] : w[k + Z], w[k + Z], wk); qb = wb >= 0 ? vanLeer(l - 2 >= 0 ? w[k - 2 * Z] : w[k - Z], w[k - Z], wk) : vanLeer(w[k + Z], wk, w[k - Z]);
          } else { qe = 0.5 * (wk + w[k + 1]); qw_ = 0.5 * (wk + w[k - 1]); qn = 0.5 * (wk + w[k + Y]); qs = 0.5 * (wk + w[k - Y]); qt = wt; qb = wb; }
          let r = -((ue * qe - uw * qw_) * rdx + (vn * qn - vs * qs) * rdyj + (wt * qt - wb * qb) * rdz), gN = (w[k + Y] - wk) * rn, gS = (wk - w[k - Y]) * rs;
          if (southWall) gS = gs1 * wk + gs2 * w[k + Y]; if (northWall) gN = -(gn1 * wk + gn2 * w[k - Y]);
          if (variableMu) {
            const mE = 0.25 * (mu[k] + mu[k + 1] + mu[k + Z] + mu[k + Z + 1]), mW = 0.25 * (mu[k] + mu[k - 1] + mu[k + Z] + mu[k + Z - 1]), mN = 0.25 * (mu[k] + mu[k + Y] + mu[k + Z] + mu[k + Z + Y]), mS = 0.25 * (mu[k] + mu[k - Y] + mu[k + Z] + mu[k + Z - Y]);
            const vis = (mE * ((w[k + 1] - wk) * rdx + (u[k + Z] - uk) * rdz) - mW * ((wk - w[k - 1]) * rdx + (u[k - 1 + Z] - u[k - 1]) * rdz)) * rdx + (mN * (gN + (v[k + Z] - vk) * rdz) - mS * (gS + (v[k - Y + Z] - v[k - Y]) * rdz)) * rdyj + (2 * mu[k + Z] * (w[k + Z] - wk) - 2 * mu[k] * (wk - w[k - Z])) * rdz * rdz;
            r += vof ? (2 * vis) / (rho[k] + rho[k + Z]) : vis;
          } else r += nu * ((w[k + 1] - 2 * wk + w[k - 1]) * rdx * rdx + (gN - gS) * rdyj + (w[k + Z] - 2 * wk + w[k - Z]) * rdz * rdz);
          rw[k] = r + fz;
        }
      }
    }
  }
  // ---- projection: makes (u, v, w) solenoidal; cdt = time increment of the (sub)step
  function project(cdt) {
    if (fu) for (let k = 0; k < N; k++) { u[k] *= fu[k]; v[k] *= fv[k]; w[k] *= fw[k]; }
    fillVel();
    const two = !!vof, s0 = two ? rho0 / cdt : 1 / cdt;
    if (two) { for (let k = 0; k < N; k++) pHat[k] = 2 * p[k] - pOld[k]; pOld.set(p); }
    const its = two ? Math.max(1, vof.pIter | 0 || 1) : 1;
    for (let it = 0; it < its; it++) {
      if (two && it > 0) pHat.set(p);
      let m = 0;
      for (let l = 1; l <= nz; l++) for (let j = 1; j <= ny; j++) { const rdyj = 1 / dy[j], row = sx * (j + sy * l); for (let i = 1; i <= nx; i++) { const k = row + i; let r = s0 * ((u[k] - u[k - 1]) * rdx + (v[k] - v[k - Y]) * rdyj + (w[k] - w[k - Z]) * rdz);
        if (two) { const rc = rho[k], fe = i < nx || perX ? (1 - (2 * rho0) / (rc + rho[k + 1])) * (pHat[k + 1] - pHat[k]) * rdx : 0, fw_ = i > 1 || perX ? (1 - (2 * rho0) / (rc + rho[k - 1])) * (pHat[k] - pHat[k - 1]) * rdx : 0, fn = j < ny || perY ? ((1 - (2 * rho0) / (rc + rho[k + Y])) * (pHat[k + Y] - pHat[k])) / dyc[j] : 0, fs = j > 1 || perY ? ((1 - (2 * rho0) / (rc + rho[k - Y])) * (pHat[k] - pHat[k - Y])) / dyc[j - 1] : 0, ft = l < nz || perZ ? (1 - (2 * rho0) / (rc + rho[k + Z])) * (pHat[k + Z] - pHat[k]) * rdz : 0, fb = l > 1 || perZ ? (1 - (2 * rho0) / (rc + rho[k - Z])) * (pHat[k] - pHat[k - Z]) * rdz : 0; r += (fe - fw_) * rdx + (fn - fs) * rdyj + (ft - fb) * rdz; }
        W[m++] = r; } }
      // W is ordered i + nx (j + ny l): reorder index m accordingly (loops above run l, j, i → m = i + nx (j + ny l))
      poisson.solve(W); S.poissonSolves++;
      m = 0; for (let l = 1; l <= nz; l++) for (let j = 1; j <= ny; j++) { const row = sx * (j + sy * l); for (let i = 1; i <= nx; i++) p[row + i] = W[m++]; }
      fill(p, 0);
    }
    for (let l = 1; l <= nz; l++) for (let j = 1; j <= ny; j++) { const rn = 1 / dyc[j], row = sx * (j + sy * l); for (let i = 1; i <= nx; i++) { const k = row + i;
      if (two) { const r0 = 1 / rho0; if (i <= nxU) { const ir = 2 / (rho[k] + rho[k + 1]); u[k] -= cdt * (r0 * (p[k + 1] - p[k]) + (ir - r0) * (pHat[k + 1] - pHat[k])) * rdx; } if (j <= nyV) { const ir = 2 / (rho[k] + rho[k + Y]); v[k] -= cdt * (r0 * (p[k + Y] - p[k]) + (ir - r0) * (pHat[k + Y] - pHat[k])) * rn; } if (l <= nzW) { const ir = 2 / (rho[k] + rho[k + Z]); w[k] -= cdt * (r0 * (p[k + Z] - p[k]) + (ir - r0) * (pHat[k + Z] - pHat[k])) * rdz; } }
      else { if (i <= nxU) u[k] -= cdt * (p[k + 1] - p[k]) * rdx; if (j <= nyV) v[k] -= cdt * (p[k + Y] - p[k]) * rn; if (l <= nzW) w[k] -= cdt * (p[k + Z] - p[k]) * rdz; } } }
    fillVel();
  }
  // ---- interface advection (one split step with the face velocities)
  let sweepOrder = 0;
  function advectInterface(dt) {
    const beta = K.thincBeta;
    const sweep = (dir) => {
      const vel = dir === 0 ? u : dir === 1 ? v : w, s = dir === 0 ? 1 : dir === 1 ? Y : Z; flux.fill(0);
      for (let l = dir === 2 ? 0 : 1; l <= nz; l++) for (let j = dir === 1 ? 0 : 1; j <= ny; j++) for (let i = dir === 0 ? 0 : 1; i <= nx; i++) {
        const k = i + sx * (j + sy * l), q = vel[k]; if (q === 0) continue;
        const d = q > 0 ? k : k + s, jd = dir === 1 ? (q > 0 ? j : j + 1) : j, h = dir === 0 ? dx : dir === 1 ? dy[jd] : dz, c = (Math.abs(q) * dt) / h, sg = q > 0 ? s : -s, cd = C[d];
        let f;
        if (cd <= 1e-8 || cd >= 1 - 1e-8) f = c * cd;
        else { const zp = d + Z < N ? C[d + Z] : cd, zm = d - Z >= 0 ? C[d - Z] : cd, gx = Math.abs(C[d + 1] - C[d - 1]) * rdx, gy = Math.abs((d + Y < N ? C[d + Y] : cd) - (d - Y >= 0 ? C[d - Y] : cd)) / dy[Math.min(Math.max(jd, 1), ny)], gz = Math.abs(zp - zm) * rdz, gsum = gx + gy + gz, wgt = gsum < 1e-12 ? 1 / 3 : (dir === 0 ? gx : dir === 1 ? gy : gz) / gsum; const iu = d - sg, id = d + sg; f = wgt * thincFlux(iu >= 0 && iu < N ? C[iu] : cd, cd, id >= 0 && id < N ? C[id] : cd, c, beta) + (1 - wgt) * c * cd; }
        flux[k] = (q > 0 ? f : -f) * h;
      }
      if (dir === 0 && perX) for (let l = 1; l <= nz; l++) for (let j = 1; j <= ny; j++) { const r = sx * (j + sy * l); flux[r] = flux[r + nx]; } // the face shared across a periodic boundary carries one flux
      if (dir === 1 && perY) for (let l = 1; l <= nz; l++) for (let i = 1; i <= nx; i++) { const r = i + Z * l; flux[r] = flux[r + ny * Y]; }
      if (dir === 2 && perZ) for (let j = 1; j <= ny; j++) for (let i = 1; i <= nx; i++) { const r = i + Y * j; flux[r] = flux[r + nz * Z]; }
      for (let l = 1; l <= nz; l++) for (let j = 1; j <= ny; j++) { const h = dir === 0 ? dx : dir === 1 ? dy[j] : dz; for (let i = 1; i <= nx; i++) { const k = i + sx * (j + sy * l); C[k] += (-(flux[k] - flux[k - s]) + cc[k] * (vel[k] - vel[k - s]) * dt) / h; } }
      for (let l = 1; l <= nz; l++) for (let j = 1; j <= ny; j++) for (let i = 1; i <= nx; i++) { const k = i + sx * (j + sy * l); if (C[k] < 0 && C[k] > -1e-12) C[k] = 0; else if (C[k] > 1 && C[k] < 1 + 1e-12) C[k] = 1; }
      fill(C, 0);
    };
    for (let k = 0; k < N; k++) cc[k] = C[k] > 0.5 ? 1 : 0; // one indicator for the three sweeps: their dilatation terms then sum to cc × div u = 0
    const ord = [[0, 1, 2], [1, 2, 0], [2, 0, 1], [0, 2, 1], [2, 1, 0], [1, 0, 2]][sweepOrder++ % 6]; for (const d of ord) sweep(d);
  }
  // ---- diagnostics
  function maxRates() { // advective and viscous rates that limit the explicit step
    let a = 0, d = 0; const nuBase = vof ? Math.max(vof.muL / vof.rhoL, vof.muG / vof.rhoG) : nu0;
    for (let l = 1; l <= nz; l++) for (let j = 1; j <= ny; j++) { const rdyj = 1 / dy[j], h2 = rdx * rdx + rdyj * rdyj + rdz * rdz; for (let i = 1; i <= nx; i++) { const k = i + sx * (j + sy * l), r = Math.abs(u[k]) * rdx + Math.abs(v[k]) * rdyj + Math.abs(w[k]) * rdz; if (r > a) a = r; const q = (nuBase + (nut ? nut[k] : 0)) * h2; if (q > d) d = q; } }
    return { adv: a, diff: d };
  }
  function divergence() { let m = 0; for (let l = 1; l <= nz; l++) for (let j = 1; j <= ny; j++) for (let i = 1; i <= nx; i++) { const k = i + sx * (j + sy * l), d = Math.abs((u[k] - u[k - 1]) * rdx + (v[k] - v[k - Y]) / dy[j] + (w[k] - w[k - Z]) * rdz) * (fcell ? (fcell[k] >= 1 ? 1 : 0) : 1); if (d > m) m = d; } return m; }
  function bulkVelocity() { let s = 0; for (let l = 1; l <= nz; l++) for (let j = 1; j <= ny; j++) { const a = dy[j] * dz * dx; for (let i = 1; i <= nx; i++) s += u[i + sx * (j + sy * l)] * a; } return s / fluidVolume; }
  function wallShear() { // mean du/dy at the two y walls (positive = drag on a flow in +x), averaged over the walls
    if (bc.y !== 'wall') return 0; let s = 0;
    for (let l = 1; l <= nz; l++) for (let i = 1; i <= nx; i++) { const a = idx(i, 1, l), b = idx(i, ny, l); s += wall3 ? gs1 * u[a] + gs2 * u[a + Y] + gn1 * u[b] + gn2 * u[b - Y] : u[a] / y1s + u[b] / y1n; }
    return s / (2 * nx * nz);
  }
  function kinetic() { // volume-averaged kinetic energy per unit mass from the face velocities (the quantity the advection scheme conserves); per unit volume from cell-centre values with two phases
    let e = 0, vol = 0;
    if (!vof) { for (let l = 1; l <= nz; l++) for (let j = 1; j <= ny; j++) { const a = dy[j], b = j <= nyV ? dyc[j] : 0; for (let i = 1; i <= nx; i++) { const k = i + sx * (j + sy * l); e += 0.5 * ((i <= nxU ? u[k] * u[k] : 0) * a + v[k] * v[k] * b + (l <= nzW ? w[k] * w[k] : 0) * a); vol += a; } } return e / vol; }
    for (let l = 1; l <= nz; l++) for (let j = 1; j <= ny; j++) { const a = cellVolume(j); for (let i = 1; i <= nx; i++) { const k = i + sx * (j + sy * l), uc = 0.5 * (u[k] + u[k - 1]), vc = 0.5 * (v[k] + v[k - Y]), wc = 0.5 * (w[k] + w[k - Z]); e += 0.5 * (uc * uc + vc * vc + wc * wc) * a * rho[k]; vol += a; } } return e / vol;
  }
  function enstrophy() { // volume average of ½ |ω|² (vorticity from the face velocities, evaluated on the cell edges and averaged)
    let e = 0, vol = 0; for (let l = 1; l <= nz; l++) for (let j = 1; j <= ny; j++) { const a = cellVolume(j), rn = 1 / dyc[j]; for (let i = 1; i <= nx; i++) { const k = i + sx * (j + sy * l), ox = (w[k + Y] - w[k]) * rn - (v[k + Z] - v[k]) * rdz, oy = (u[k + Z] - u[k]) * rdz - (w[k + 1] - w[k]) * rdx, oz = (v[k + 1] - v[k]) * rdx - (u[k + Y] - u[k]) * rn; e += 0.5 * (ox * ox + oy * oy + oz * oz) * a; vol += a; } } return e / vol;
  }
  function liquidVolume() { if (!vof) return 0; let s = 0; for (let l = 1; l <= nz; l++) for (let j = 1; j <= ny; j++) { const a = cellVolume(j); for (let i = 1; i <= nx; i++) s += C[i + sx * (j + sy * l)] * a; } return s; }
  function updateUTau() { if (bc.y === 'wall' && !pipe && !vof) S.uTau = Math.sqrt(Math.abs(nu0 * wallShear())); else if (pipe) S.uTau = Math.sqrt(Math.abs((force[0] + S.forcingMeanInst) * pipe.radius * 0.5)); }
  S.forcingMeanInst = 0;
  // ---- one time step
  let first = true, volume0 = 0;
  function step(dtFixed) {
    if (first) { fillVel(); if (vof) { properties(); volume0 = liquidVolume(); if (o.vof.hydrostaticStart !== false) { if (grav[1] < 0 && grav[0] === 0 && grav[2] === 0 && !perY) { for (let l = 0; l < sz; l++) for (let i = 0; i < sx; i++) { let ph = 0; for (let j = ny; j >= 1; j--) { const k = idx(i, j, l); ph += -grav[1] * (j === ny ? 0.5 * rho[k] * dy[j] : 0.5 * (rho[k] * dy[j] + rho[k + Y] * dy[j + 1])); p[k] = ph; } } fill(p, 0); pOld.set(p); } const keep = [u.slice(), v.slice(), w.slice()], n0 = Math.max(0, vof.startIterations ?? 20); for (let q = 0; q < n0; q++) { u.set(keep[0]); v.set(keep[1]); w.set(keep[2]); const g = 1e-3; for (let k = 0; k < N; k++) { u[k] += g * grav[0]; v[k] += g * grav[1]; w[k] += g * grav[2]; } pOld.set(p); project(g); } u.set(keep[0]); v.set(keep[1]); w.set(keep[2]); pOld.set(p); fillVel(); } } updateUTau(); eddyViscosity(0); properties(); first = false; }
    const rate = maxRates(), lim = scheme === 'rk3' ? 0.55 : 0.2; let dt = Math.min(cfl / Math.max(rate.adv, 1e-30), lim / Math.max(rate.diff, 1e-30), o.dtMax ?? Infinity);
    if (vof) { const gm = Math.hypot(grav[0], grav[1], grav[2]); if (gm > 0) dt = Math.min(dt, 0.3 * Math.sqrt(Math.min(dx, dyMin, dz) / gm)); }
    if (dtFixed > 0) dt = dtFixed;
    if (scheme === 'rk3') {
      const A = K.rk3.a, B = K.rk3.b, Cs = K.rk3.c;
      for (let s = 0; s < 3; s++) {
        rhs();
        for (let k = 0; k < N; k++) { qu[k] = A[s] * qu[k] + dt * ru[k]; qv[k] = A[s] * qv[k] + dt * rv[k]; qw[k] = A[s] * qw[k] + dt * rw[k]; u[k] += B[s] * qu[k]; v[k] += B[s] * qv[k]; w[k] += B[s] * qw[k]; }
        project(Cs[s] * dt);
      }
    } else {
      if (vof) { advectInterface(dt); properties(); }
      rhs();
      const c1 = S.steps === 0 ? 1 : 1.5, c0 = S.steps === 0 ? 0 : -0.5;
      for (let k = 0; k < N; k++) { u[k] += dt * (c1 * ru[k] + c0 * qu[k]); v[k] += dt * (c1 * rv[k] + c0 * qv[k]); w[k] += dt * (c1 * rw[k] + c0 * qw[k]); qu[k] = ru[k]; qv[k] = rv[k]; qw[k] = rw[k]; }
      project(dt);
    }
    if (o.bulk !== undefined && o.bulk !== null && perX) { // hold the bulk velocity: uniform shift of the fluid, equivalent to a body force
      const ub = bulkVelocity(), du = o.bulk - ub, n = fu || null; let q = 1;
      if (n) { let s = 0; for (let l = 1; l <= nz; l++) for (let j = 1; j <= ny; j++) { const a = dy[j] * dz * dx; for (let i = 1; i <= nx; i++) s += fu[i + sx * (j + sy * l)] * a; } q = fluidVolume / s; }
      for (let l = 1; l <= nz; l++) for (let j = 1; j <= ny; j++) for (let i = 1; i <= nx; i++) { const k = i + sx * (j + sy * l); u[k] += du * q * (n ? n[k] : 1); }
      fill(u, 1); S.forcingMeanInst = S.forcing + du / dt; S.forcing = S.forcingMeanInst; // carried as a body force in the next step so that the correction stays small
    } else S.forcingMeanInst = 0;
    S.t += dt; S.steps++; S.dt = dt; S.forcingMean += (force[0] + S.forcing) * dt; S.forcingTime += dt;
    updateUTau(); eddyViscosity(dt); if (!vof) properties();
    return dt;
  }
  /** Plane averages over x and z: { y[], u[], v[], w[], uu[], vv[], ww[], uv[], nut[] } (second moments about zero, at the cell centres). */
  function profile() {
    const out = { y: [], u: [], v: [], w: [], uu: [], vv: [], ww: [], uv: [], nut: [] }, n = nx * nz;
    for (let j = 1; j <= ny; j++) { let a = 0, b = 0, c = 0, aa = 0, bb = 0, cq = 0, ab = 0, e = 0; for (let l = 1; l <= nz; l++) for (let i = 1; i <= nx; i++) { const k = i + sx * (j + sy * l), uc = 0.5 * (u[k] + u[k - 1]), vc = 0.5 * (v[k] + v[k - Y]), wc = 0.5 * (w[k] + w[k - Z]); a += uc; b += vc; c += wc; aa += uc * uc; bb += vc * vc; cq += wc * wc; ab += uc * vc; if (nut) e += nut[k]; }
      out.y.push(yc[j]); out.u.push(a / n); out.v.push(b / n); out.w.push(c / n); out.uu.push(aa / n); out.vv.push(bb / n); out.ww.push(cq / n); out.uv.push(ab / n); out.nut.push(e / n); }
    return out;
  }
  /** A plane of a field as nested arrays: field 'u' | 'v' | 'w' | 'p' | 'c' | 'nut' | 'speed'; normal 'x' | 'y' | 'z'; index = cell index (1-based) of the plane. Returns { a[], b[], z[][] } (rows along b). */
  function slice(field, normal = 'z', index) {
    const get = (i, j, l) => { const k = idx(i, j, l); if (field === 'u') return 0.5 * (u[k] + u[k - 1]); if (field === 'v') return 0.5 * (v[k] + v[k - Y]); if (field === 'w') return 0.5 * (w[k] + w[k - Z]); if (field === 'p') return p[k]; if (field === 'c') return C ? C[k] : 0; if (field === 'nut') return nut ? nut[k] : 0; return Math.hypot(0.5 * (u[k] + u[k - 1]), 0.5 * (v[k] + v[k - Y]), 0.5 * (w[k] + w[k - Z])); };
    const X = Array.from(xc.subarray(1, nx + 1)), Yv = Array.from(yc.subarray(1, ny + 1)), Zv = Array.from(zc.subarray(1, nz + 1)), z = [];
    if (normal === 'z') { const l = index ?? Math.max(1, nz >> 1); for (let j = 1; j <= ny; j++) { const r = []; for (let i = 1; i <= nx; i++) r.push(get(i, j, l)); z.push(r); } return { a: X, b: Yv, z }; }
    if (normal === 'x') { const i = index ?? Math.max(1, nx >> 1); for (let j = 1; j <= ny; j++) { const r = []; for (let l = 1; l <= nz; l++) r.push(get(i, j, l)); z.push(r); } return { a: Zv, b: Yv, z }; }
    const j = index ?? Math.max(1, ny >> 1); for (let l = 1; l <= nz; l++) { const r = []; for (let i = 1; i <= nx; i++) r.push(get(i, j, l)); z.push(r); } return { a: X, b: Zv, z };
  }
  return { grid: { nx, ny, nz, lx, ly, lz, dx, dz, yc, yf, dy, xc, zc, sx, sy, sz, idx, cells: nx * ny * nz, fastPoisson: poisson.fast, fluidVolume }, state: S, u, v, w, p, c: C, nut, nuTilde: nuT, lengthScale: lDes, fluidFraction: fcell, wallDistance: dWall, force,
    step, project, fillVel, profile, slice, kinetic, enstrophy, divergence, bulkVelocity, wallShear, liquidVolume, initialVolume: () => volume0, scheme, sgs, nu: nu0 };
}

// ---- cases -----------------------------------------------------------------------------------------------------------------------
const now = () => (typeof performance !== 'undefined' && performance.now ? performance.now() : Date.now());
/** Wraps a solver in a steppable run: { done, advance(maxSteps or ms budget), progress, result() }. */
function runner(sim, tEnd, maxSteps, each, finish) {
  const R = { sim, done: false, wall: 0, progress: 0,
    advance(n = 50, budgetMs = Infinity) { // up to n steps or budgetMs of wall time, whichever comes first
      const t0 = now(); let k = 0;
      while (!R.done && k < n && now() - t0 < budgetMs) {
        if (tEnd - sim.state.t <= 1e-12 * Math.max(tEnd, 1)) { R.done = true; break; }
        const dt = sim.step(); k++; if (each) each(sim, dt);
        if (!Number.isFinite(dt) || !Number.isFinite(sim.u[sim.grid.idx(1, 1, 1)])) { R.failed = true; R.done = true; }
        if (sim.state.steps >= maxSteps || sim.state.t >= tEnd * (1 - 1e-12)) R.done = true;
      }
      R.wall += now() - t0; R.progress = Math.min(1, Math.max(sim.state.t / tEnd, sim.state.steps / maxSteps)); return R;
    },
    result() { return finish(sim, R); } };
  return R;
}
const toArr = (a) => Array.from(a);

/**
 * Taylor–Green vortex in a periodic box of side 2π L.
 * dim 2: u = U sin x cos y, v = −U cos x sin y (the third direction is uniform); exact decay exp(−2 ν t) of the velocity.
 * dim 3: u = U sin x cos y cos z, v = −U cos x sin y cos z, w = 0 (Brachet et al. 1983); Re = U L / ν.
 * o: { n (cells per side), re, dim, tEnd (in L/U), cfl, sgs }. Returns a run (see runner): result() =
 *    { t[], ek[] (kinetic energy / U²), dissipation[] (−dEk/dt from the energy history), enstrophyDissipation[] (2 ν × enstrophy),
 *      exactEk[] (2-D), errorMax (2-D: largest relative velocity-amplitude error), peak: { t, value }, steps, cells, msPerStep, divergenceMax }.
 */
export function taylorGreen(o = {}) {
  const n = o.n || 32, dim = o.dim === 2 ? 2 : 3, re = o.re || 100, tEnd = o.tEnd ?? (dim === 2 ? 1 : 10), L = 2 * Math.PI, nu = 1 / re, nz = dim === 2 ? Math.max(2, o.nz || 2) : n;
  const sim = createFlow3D({ nx: n, ny: n, nz, lx: L, ly: L, lz: dim === 2 ? (L * nz) / n : L, bc: { x: 'periodic', y: 'periodic', z: 'periodic' }, nu, sgs: o.sgs || 'none', cfl: o.cfl ?? 0.6, dtMax: o.dtMax,
    init: (x, y, z) => (dim === 2 ? [Math.sin(x) * Math.cos(y), -Math.cos(x) * Math.sin(y), 0] : [Math.sin(x) * Math.cos(y) * Math.cos(z), -Math.cos(x) * Math.sin(y) * Math.cos(z), 0]) });
  sim.fillVel(); const hist = { t: [0], ek: [sim.kinetic()], ens: [sim.enstrophy()] }; let div = 0;
  return runner(sim, tEnd, o.maxSteps || 200000, (s) => { hist.t.push(s.state.t); hist.ek.push(s.kinetic()); hist.ens.push(s.enstrophy()); if (s.state.steps % 10 === 0) div = Math.max(div, s.divergence()); }, (s, R) => {
    const m = hist.t.length, diss = new Array(m).fill(0); for (let i = 1; i < m - 1; i++) diss[i] = -(hist.ek[i + 1] - hist.ek[i - 1]) / (hist.t[i + 1] - hist.t[i - 1]); diss[0] = 2 * nu * hist.ens[0]; if (m > 1) diss[m - 1] = diss[m - 2];
    let pk = 0; for (let i = 1; i < m; i++) if (diss[i] > diss[pk]) pk = i;
    const out = { dim, n, re, t: hist.t, ek: hist.ek, dissipation: diss, enstrophyDissipation: hist.ens.map((e) => 2 * nu * e), peak: { t: hist.t[pk], value: diss[pk] }, steps: s.state.steps, cells: s.grid.cells, msPerStep: R.wall / Math.max(s.state.steps, 1), divergenceMax: Math.max(div, s.divergence()), fastPoisson: s.grid.fastPoisson, failed: !!R.failed };
    if (dim === 2) { out.exactEk = hist.t.map((t) => hist.ek[0] * Math.exp(-4 * nu * t)); out.errorMax = Math.max(...hist.t.map((t, i) => Math.abs(Math.sqrt(hist.ek[i] / hist.ek[0]) / Math.exp(-2 * nu * t) - 1))); }
    return out;
  });
}

/**
 * Plane channel between two walls at y = 0 and 2 h (h = 1), periodic in x and z; wall units from Re_τ = u_τ h / ν.
 * o: { reTau (180), nx, ny, nz, lx (2π), lz (π), stretchY (1.9), sgs, mode: 'gradient' (constant forcing u_τ² / h, default) | 'bulk' (constant
 *      flow rate; bulkPlus = U_b / u_τ target), laminar: true (parabolic start, no perturbation), tEnd and tStats (in h / u_τ), cfl, seed }
 * Returns a run; result() = { yPlus[], uPlus[], uRms[], vRms[], wRms[], uvPlus[], nutPlus[] (lower half, averaged over both halves and over
 *   the sampling time), reTauActual, uBulkPlus, uCentrePlus, cf (2 τ_w / ρ U_b²), reBulk (U_b 2h / ν), resolution: { dxPlus, dzPlus,
 *   dyPlusWall, dyPlusMax }, history: { t[], uBulk[], reTau[], ek[] }, samples, steps, cells, msPerStep, usPerCellStep, balance }.
 */
export function channel3D(o = {}) {
  const reTau = o.reTau || 180, nx = o.nx || 32, ny = o.ny || 32, nz = o.nz || 32, lx = o.lx ?? 2 * Math.PI, lz = o.lz ?? Math.PI, nu = 1 / reTau, laminar = !!o.laminar, mode = o.mode === 'bulk' ? 'bulk' : 'gradient';
  const tEnd = o.tEnd ?? 12, tStats = Math.min(o.tStats ?? 0.5 * tEnd, tEnd), ubLog = Math.log(reTau) / 0.41 + 5.2 - 1 / 0.41 - 0.6, ub0 = laminar ? reTau / 3 : o.bulkPlus || ubLog;
  let seed = (o.seed || 12345) >>> 0; const rnd = () => { seed = (Math.imul(seed, 1664525) + 1013904223) >>> 0; return seed / 4294967296 - 0.5; };
  const mean = (y) => { const d = Math.min(y, 2 - y); if (laminar) return 0.5 * reTau * y * (2 - y); const yp = d * reTau; return (1 / 0.41) * Math.log(1 + 0.41 * yp) + 7.8 * (1 - Math.exp(-yp / 11) - (yp / 11) * Math.exp(-yp / 3)); }; // Reichardt's profile
  const amp = laminar ? 0 : o.perturbation ?? 0.2 * ub0, ax = (2 * Math.PI) / lx, az = (2 * Math.PI) / lz;
  const sim = createFlow3D({ nx, ny, nz, lx, ly: 2, lz, stretchY: o.stretchY ?? (laminar ? 0 : 1.9), nu, sgs: o.sgs || 'none', cfl: o.cfl ?? 0.8, dtMax: o.dtMax, wallOrder: o.wallOrder, cs: o.cs, cw: o.cw, cDes: o.cDes,
    force: mode === 'gradient' ? [1, 0, 0] : [0, 0, 0], bulk: mode === 'bulk' ? ub0 : undefined, nuTilde0: (d) => Math.min(0.41 * d, 0.09) * (1 - Math.exp(-(d * reTau) / 19)) ** 2,
    init: (x, y, z) => { const sh = Math.sin(Math.PI * y * 0.5), s2 = sh * sh; return [mean(y) + amp * (Math.cos(az * z) * s2 * 0.5 + 0.25 * Math.sin(ax * x) * Math.sin(2 * az * z) * s2) + amp * 0.1 * rnd() * s2, amp * 0.4 * Math.sin(ax * x) * Math.cos(az * z) * s2 * 0 + amp * 0.1 * rnd() * s2, amp * (0.5 * Math.sin(ax * x) * s2 + 0.3 * Math.sin(2 * ax * x) * Math.cos(az * z) * s2) + amp * 0.1 * rnd() * s2]; } });
  sim.project(1); // the perturbed start is made solenoidal
  const acc = { n: 0, time: 0, u: new Float64Array(ny), uu: new Float64Array(ny), vv: new Float64Array(ny), ww: new Float64Array(ny), uv: new Float64Array(ny), nut: new Float64Array(ny), tau: 0 }, hist = { t: [], uBulk: [], reTau: [], ek: [] };
  const every = Math.max(1, o.sampleEvery || 5);
  return runner(sim, tEnd, o.maxSteps || 400000, (s, dt) => {
    const st = s.state;
    if (st.steps % every === 0) { const tw = nu * s.wallShear(); hist.t.push(st.t); hist.uBulk.push(s.bulkVelocity()); hist.reTau.push(Math.sqrt(Math.abs(tw)) * reTau); hist.ek.push(s.kinetic());
      if (st.t >= tEnd - tStats) { const pr = s.profile(); for (let j = 0; j < ny; j++) { acc.u[j] += pr.u[j]; acc.uu[j] += pr.uu[j]; acc.vv[j] += pr.vv[j]; acc.ww[j] += pr.ww[j]; acc.uv[j] += pr.uv[j]; acc.nut[j] += pr.nut[j]; } acc.tau += tw; acc.n++; } }
  }, (s, R) => {
    const n = Math.max(acc.n, 1), g = s.grid, half = ny >> 1, tauW = acc.n ? acc.tau / n : nu * s.wallShear(), ut = Math.sqrt(Math.abs(tauW)), U = Array.from(acc.u, (x) => x / n);
    if (!acc.n) { const pr = s.profile(); for (let j = 0; j < ny; j++) { U[j] = pr.u[j]; acc.uu[j] = pr.uu[j]; acc.vv[j] = pr.vv[j]; acc.ww[j] = pr.ww[j]; acc.uv[j] = pr.uv[j]; acc.nut[j] = pr.nut[j]; } }
    let ub = 0; for (let j = 0; j < ny; j++) ub += U[j] * g.dy[j + 1]; ub /= 2;
    const sym = (a, sgn = 1) => { const r = []; for (let j = 0; j < half; j++) r.push(0.5 * (a[j] + sgn * a[ny - 1 - j])); return r; }, Um = sym(U), rms = (a) => sym(Array.from(a, (x, j) => x / n)).map((x, j) => x);
    const uu = rms(acc.uu).map((x, j) => Math.sqrt(Math.max(x - Um[j] * Um[j], 0)) / ut), vv = rms(acc.vv).map((x) => Math.sqrt(Math.max(x, 0)) / ut), ww = rms(acc.ww).map((x) => Math.sqrt(Math.max(x, 0)) / ut), uv = sym(Array.from(acc.uv, (x) => x / n), -1).map((x) => -x / (ut * ut));
    const reA = ut * reTau, yPlus = []; for (let j = 0; j < half; j++) yPlus.push(g.yc[j + 1] * reA);
    const cells = g.cells, ms = R.wall / Math.max(s.state.steps, 1);
    return { reTau, sgs: s.sgs, mode, nx, ny, nz, lx, lz, y: toArr(g.yc.subarray(1, half + 1)), yPlus, uPlus: Um.map((x) => x / ut), uRms: uu, vRms: vv, wRms: ww, uvPlus: uv, nutPlus: sym(Array.from(acc.nut, (x) => x / n)).map((x) => x * reTau),
      reTauActual: reA, uTau: ut, uBulk: ub, uBulkPlus: ub / ut, uCentrePlus: Um[half - 1] / ut, cf: (2 * tauW) / (ub * ub), reBulk: 2 * ub * reTau, forcingMean: s.state.forcingMean / Math.max(s.state.forcingTime, 1e-30),
      slice: (() => { const q = s.slice('u', 'x'); return { z: q.a, y: q.b, u: q.z.map((row) => row.map((x) => x / ut)) }; })(), resolution: { dxPlus: g.dx * reA, dzPlus: g.dz * reA, dyPlusWall: g.dy[1] * reA, dyPlusMax: g.dy[half] * reA, firstCentrePlus: g.yc[1] * reA }, history: hist, samples: acc.n, time: s.state.t, steps: s.state.steps, cells, msPerStep: ms, usPerCellStep: (1000 * ms) / cells, divergenceMax: s.divergence(), failed: !!R.failed,
      balance: { wallShear: tauW, forcing: s.state.forcingMean / Math.max(s.state.forcingTime, 1e-30) } };
  });
}

/**
 * Laminar or turbulent flow in a circular pipe along x (periodic), represented by the immersed boundary in a square box.
 * o: { radius (0.5), n (cells across the box side), nx, lx, nu, bulk (mean velocity, held) or gradient (driving acceleration −(1/ρ) dp/dx),
 *      tEnd, sgs, init: 'parabola' | 'rest', steady (stop when the forcing / bulk velocity changes by less than this per unit time) }
 * Returns a run; result() = { bulk, gradient (−(1/ρ) dp/dx), friction (Darcy), reynolds, frictionLaminar (64/Re), error, profile: { r[], u[] }, … }.
 */
export function pipe3D(o = {}) {
  const R0 = o.radius ?? 0.5, n = o.n || 24, nx = o.nx || 4, margin = o.margin ?? 1.5, side = 2 * R0 + (2 * margin * 2 * R0) / n, lx = o.lx ?? (nx * side) / n, nu = o.nu ?? 0.01, useBulk = o.gradient === undefined, ub = o.bulk ?? 1, tEnd = o.tEnd ?? 2;
  const cy = side / 2, cz = side / 2, uInit = o.init === 'rest' ? () => [0, 0, 0] : (x, y, z) => { const r2 = ((y - cy) ** 2 + (z - cz) ** 2) / (R0 * R0); return [r2 < 1 ? 2 * (useBulk ? ub : (o.gradient * R0 * R0) / (8 * nu)) * (1 - r2) : 0, 0, 0]; };
  const sim = createFlow3D({ nx, ny: n, nz: n, lx, ly: side, lz: side, bc: { x: 'periodic', y: 'slip', z: 'slip' }, nu, sgs: o.sgs || 'none', cfl: o.cfl ?? 0.6, pipe: { radius: R0, cy, cz, samples: o.samples || 8 }, force: useBulk ? [0, 0, 0] : [o.gradient, 0, 0], bulk: useBulk ? ub : undefined, init: uInit, wallOrder: 2 });
  const hist = { t: [], forcing: [], bulk: [] }; let last = null, steadyAt = null;
  const run = runner(sim, tEnd, o.maxSteps || 400000, (s) => { const st = s.state; if (st.steps % 5 === 0) { const q = useBulk ? st.forcingMeanInst : s.bulkVelocity(); hist.t.push(st.t); hist.forcing.push(st.forcingMeanInst + s.force[0]); hist.bulk.push(s.bulkVelocity()); if (o.steady && last !== null && st.steps > 20 && Math.abs(q - last.q) / Math.max(Math.abs(q), 1e-30) / Math.max(st.t - last.t, 1e-30) < o.steady) { steadyAt = st.t; run.done = true; } last = { q, t: st.t }; } }, (s, Rr) => {
    const g = s.grid, ubA = s.bulkVelocity(), grad = s.force[0] + (useBulk ? s.state.forcingMeanInst : 0), D = 2 * R0, re = (ubA * D) / nu, f = (2 * D * grad) / (ubA * ubA), prof = { r: [], u: [] }, l = Math.max(1, Math.round(n / 2));
    for (let j = 1; j <= n; j++) { const k = g.idx(1, j, l); prof.r.push((g.yc[j] - cy) / R0); prof.u.push(0.5 * (s.u[k] + s.u[k + g.sx * g.sy]) ); }
    return { n, nx, cells: g.cells, cellsAcross: (2 * R0) / g.dz, radius: R0, nu, bulk: ubA, gradient: grad, friction: f, reynolds: re, frictionLaminar: 64 / re, error: f / (64 / re) - 1, profile: prof, history: hist, steadyAt, time: s.state.t, steps: s.state.steps, msPerStep: Rr.wall / Math.max(s.state.steps, 1), areaError: g.fluidVolume / (Math.PI * R0 * R0 * lx) - 1, failed: !!Rr.failed };
  });
  return run;
}

/**
 * Two-phase cases with the volume-of-fluid interface (water-like liquid under a gas, gravity in −y).
 * kind 'dambreak': a liquid column a × 2a × lz collapses in a closed box 8a × 4a × lz (3-D: the column may be given a spanwise shape).
 * kind 'stratified': a liquid layer under a gas in a periodic channel driven by a body force, with an initial interface wave.
 * o: { kind, n (cells over the height), a, rhoL, rhoG, muL, muG, g, tEnd, sgs, spanCells, shape3d, level, wave, force, lx, lz }.
 * Returns a run; result() = { t[], front[] (dam break: position of the leading edge on the floor), height[] (column height at the wall),
 *   volume: { initial, final, error }, holdup, kinetic[], slice (liquid fraction on the mid-plane), steps, cells, msPerStep, … }.
 */
export function twoPhase3D(o = {}) {
  const kind = o.kind === 'stratified' ? 'stratified' : 'dambreak', g = o.g ?? 9.80665, n = o.n || 24, rhoL = o.rhoL ?? 998.2, rhoG = o.rhoG ?? 1.204, muL = o.muL ?? 1.0e-3, muG = o.muG ?? 1.82e-5;
  let sim, a, ly, lx, lz, nz, nx, tEnd;
  if (kind === 'dambreak') {
    a = o.a ?? 0.1; lx = 8 * a; ly = 4 * a; nx = 2 * n; nz = Math.max(2, o.spanCells || 4); lz = o.lz ?? (ly * nz) / n; tEnd = o.tEnd ?? 3 * Math.sqrt(a / g);
    const shape = o.shape3d ? (z) => a * (1 + 0.25 * Math.cos((2 * Math.PI * z) / lz)) : () => a;
    sim = createFlow3D({ nx, ny: n, nz, lx, ly, lz, bc: { x: 'wall', y: 'wall', z: o.shape3d || nz > 4 ? 'slip' : 'periodic' }, sgs: o.sgs || 'none', cfl: o.cfl ?? 0.2, cs: o.cs,
      vof: { rhoL, rhoG, muL, muG, gravity: [0, -g, 0], pIter: o.pIter || 2, c0: (x, y, z, hx, hy) => { const w = shape(z), fx = Math.min(1, Math.max(0, (w - (x - hx / 2)) / hx)), fy = Math.min(1, Math.max(0, (2 * a - (y - hy / 2)) / hy)); return fx * fy; } } });
  } else {
    const h = o.height ?? 0.05; ly = h; lx = o.lx ?? 4 * h; nx = Math.round((n * lx) / ly); nz = Math.max(2, o.spanCells || Math.max(4, n >> 1)); lz = o.lz ?? (ly * nz) / n; a = h; tEnd = o.tEnd ?? 1;
    const lev = (o.level ?? 0.4) * h, wv = (o.wave ?? 0.08) * h;
    sim = createFlow3D({ nx, ny: n, nz, lx, ly, lz, bc: { x: 'periodic', y: 'wall', z: 'periodic' }, sgs: o.sgs || 'none', cfl: o.cfl ?? 0.2, cs: o.cs, force: [o.force ?? 0, 0, 0],
      init: (x, y) => [(o.uGas ?? 0) * (y > lev ? 1 : 0) + (o.uLiquid ?? 0) * (y <= lev ? 1 : 0), 0, 0],
      vof: { rhoL, rhoG, muL, muG, gravity: [0, -g, 0], pIter: o.pIter || 2, c0: (x, y, z, hx, hy) => { const s = lev + wv * Math.sin((2 * Math.PI * x) / lx) * (1 + 0.3 * Math.cos((2 * Math.PI * z) / lz)); return Math.min(1, Math.max(0, (s - (y - hy / 2)) / hy)); } } });
  }
  const hist = { t: [], front: [], height: [], kinetic: [], volume: [] }, G = sim.grid;
  const sample = (s) => { let front = 0, hgt = 0; for (let i = 1; i <= nx; i++) { let c = 0; for (let l = 1; l <= nz; l++) c += s.c[G.idx(i, 1, l)]; c /= nz; if (c > 0.5) front = Math.max(front, i * G.dx); } for (let j = 1; j <= n; j++) { let c = 0; for (let l = 1; l <= nz; l++) c += s.c[G.idx(1, j, l)]; hgt += (c / nz) * G.dy[j]; } hist.t.push(s.state.t); hist.front.push(front); hist.height.push(hgt); hist.kinetic.push(s.kinetic()); hist.volume.push(s.liquidVolume()); };
  let v0 = null;
  return runner(sim, tEnd, o.maxSteps || 200000, (s) => { if (v0 === null) v0 = s.initialVolume(); if (s.state.steps % 4 === 0) sample(s); }, (s, R) => {
    sample(s); const vf = s.liquidVolume(), vi = v0 ?? vf, sl = s.slice('c', 'z'), su = s.slice('u', 'z'), sv = s.slice('v', 'z'); let cmin = 0, cmax = 1; for (let k = 0; k < s.c.length; k++) { if (s.c[k] < cmin) cmin = s.c[k]; if (s.c[k] > cmax) cmax = s.c[k]; }
    let qL = 0, qG = 0; for (let l = 1; l <= nz; l++) for (let j = 1; j <= n; j++) { const k = G.idx(1, j, l), cf = 0.5 * (s.c[k] + s.c[k + 1]), ar = G.dy[j] * G.dz; qL += s.u[k] * cf * ar; qG += s.u[k] * (1 - cf) * ar; }
    return { kind, a, nx, ny: n, nz, lx, ly, lz, cells: G.cells, t: hist.t, front: hist.front, height: hist.height, kinetic: hist.kinetic, volume: { initial: vi, final: vf, error: (vf - vi) / vi, history: hist.volume }, holdup: vf / (lx * ly * lz), bounds: { min: cmin, max: cmax - 1 }, vsl: qL / (ly * lz), vsg: qG / (ly * lz),
      slice: { x: sl.a, y: sl.b, c: sl.z, u: su.z, v: sv.z }, time: s.state.t, steps: s.state.steps, msPerStep: R.wall / Math.max(s.state.steps, 1), divergenceMax: s.divergence(), failed: !!R.failed, densityRatio: rhoL / rhoG };
  });
}

/**
 * Two-way coupling of a one-dimensional line model with a three-dimensional pipe section.
 * The 1-D model is a line of nCells cells between a fixed inlet and outlet pressure (incompressible liquid, the flow rate follows from the
 * momentum balance ρ L dU/dt = Δp − Σ (dp/dx)_i Δx_i with the laminar / Colebrook closure in every cell). In the coupled cell the closure is
 * replaced by the pressure gradient returned by the 3-D section, which receives the flow rate (mean velocity) and the fluid properties
 * from the 1-D model. The exchange is repeated until the flow rate and the 3-D gradient are consistent.
 * o: { D, length, nCells, cell (index of the coupled cell), rho, mu, dp (inlet − outlet pressure), n (3-D cells across), nx, closure(Re) → Darcy f,
 *      tol, maxIter, relax, tEnd3d (in D²/ν units, default 0.6) }.
 * Returns { iterations: [{ u, f1d, f3d, gradient3d }], velocity, velocityAnalytic (Hagen–Poiseuille), velocity1d (closure only), error3d, converged, … }.
 */
export function coupled1D3D(o = {}) {
  const D = o.D ?? 0.1, L = o.length ?? 100, nc = o.nCells || 20, ic = Math.min(nc - 1, Math.max(0, o.cell ?? nc >> 1)), rho = o.rho ?? 850, mu = o.mu ?? 0.5, dp = o.dp ?? 2000, nu = mu / rho, dxc = L / nc;
  const closure = o.closure || ((re) => 64 / Math.max(re, 1e-9)), grad1 = (U) => (closure((U * D) / nu) * rho * U * U) / (2 * D);
  const solve1d = (g3) => { // steady state of the 1-D momentum balance with the gradient of the coupled cell given as g3(U)
    const res = (U) => dp - (nc - 1) * dxc * grad1(U) - dxc * g3(U); let lo = 0, hi = 1e-6; while (res(hi) > 0 && hi < 1e4) hi *= 2; for (let k = 0; k < 200; k++) { const m = 0.5 * (lo + hi); if (res(m) > 0) lo = m; else hi = m; } return 0.5 * (lo + hi); };
  const uAnalytic = (dp * D * D) / (32 * mu * L), u1 = solve1d(grad1), its = [], tol = o.tol ?? 1e-4, relax = o.relax ?? 1;
  let U = u1, ratio = 1, conv = false, last3 = null;
  for (let it = 0; it < (o.maxIter || 8); it++) {
    // 3-D section in units of D and U: Re = U D / ν; it returns the Darcy friction factor at that Reynolds number
    const re = (U * D) / nu, run = pipe3D({ radius: 0.5, n: o.n || 20, nx: o.nx || 2, nu: 1 / re, bulk: 1, tEnd: (o.tEnd3d ?? 0.6) * re, steady: o.steady ?? 1e-5 * 1, cfl: 0.6 }); while (!run.done) run.advance(200);
    const r3 = run.result(), f3 = r3.friction, g3 = (f3 * rho * U * U) / (2 * D); last3 = r3; ratio = f3 / closure(re);
    its.push({ u: U, reynolds: re, f1d: closure(re), f3d: f3, gradient3d: g3, gradient1d: grad1(U), steps3d: r3.steps });
    const Un = solve1d((V) => ratio * grad1(V)), change = Math.abs(Un - U) / Math.max(Un, 1e-30); U = U + relax * (Un - U);
    if (change < tol) { conv = true; break; }
  }
  return { D, length: L, nCells: nc, cell: ic, rho, mu, dp, iterations: its, velocity: U, velocityAnalytic: uAnalytic, velocity1d: u1, error1d: u1 / uAnalytic - 1, errorCoupled: U / uAnalytic - 1, expectedCoupled: 1 / (1 + (ratio - 1) / nc) - 1, ratio3d: ratio, error3d: ratio - 1, converged: conv, cells3d: last3 ? last3.cells : 0, cellsAcross: last3 ? last3.cellsAcross : 0, msPerStep3d: last3 ? last3.msPerStep : 0, pressureDrop3dCell: its.length ? its[its.length - 1].gradient3d * dxc : 0, pressureDropBalance: { given: dp, sum: (nc - 1) * dxc * grad1(U) + dxc * ratio * grad1(U) } };
}
