// Reference industrial case used for every suite's default inputs: a deep-water oil tie-back, 18.4 km of
// insulated 10-inch flowline on an undulating seabed at about 1,300 m water depth, a 1,370 m steel catenary
// riser and a topsides separator. All suites start from the same numbers so that they agree out of the box.
export const BASE = Object.freeze({
  name: 'Deep-water oil tie-back (reference case)',
  // elevation profile: horizontal distance from the wellhead/manifold (m), elevation relative to mean sea level (m)
  profile: Object.freeze([
    { x: 0, z: -1250 }, { x: 3000, z: -1280 }, { x: 6000, z: -1265 }, { x: 9000, z: -1310 }, { x: 12000, z: -1295 }, { x: 15000, z: -1340 }, { x: 18000, z: -1350 }, { x: 18400, z: -1348 },
    { x: 18700, z: -1300 }, { x: 19000, z: -1100 }, { x: 19250, z: -800 }, { x: 19450, z: -400 }, { x: 19560, z: 0 }, { x: 19570, z: 25 },
  ].map(Object.freeze)),
  riserBaseX: 18400, // m: where the flowline meets the riser
  idMm: 254, wtMm: 15.9, roughUm: 45, // 10-inch flowline, carbon steel
  insulation: Object.freeze({ t: 0.08, k: 0.17, name: 'Wet insulation (GSPU)' }), // 80 mm
  U: 3.0, // W/m²/K referred to the inner diameter
  waterDepth: 1350, riserHeight: 1375,
  tSeabed: 4, tSeaSurface: 24, tAir: 26, currentSpeed: 0.3,
  tIn: 70, pOut: 25, // °C at the flowline inlet, bara at the topsides separator
  // reservoir and well
  pRes: 300, tRes: 90, pi: 25, // bara, °C, productivity index Sm³/d/bar (stock-tank liquid)
  wellTVD: 2600, wellMD: 3400, tubingIdMm: 114.3, // below mudline
  wells: 3, // producing wells sharing the case rate equally: the smallest count that holds it within a 70 bar drawdown with the chokes about 60 % open and subsea boosting (well-count study of the network suite)
  // steel: API 5L X65
  smys: 450, smts: 535, E: 207000, poisson: 0.3, alphaT: 1.17e-5, rhoSteel: 7850, corrosionAllowanceMm: 3, designPressure: 345, designTemp: 110, // MPa, MPa, MPa, -, 1/K, kg/m³, mm, bara, °C
  // topsides
  separatorP: 25, slugCatcherVol: 60, // bara, m³
  // economics
  oilPrice: 75, gasPrice: 4, discountRate: 10, projectLife: 20, // $/bbl, $/MMBtu, %/y, y
});
