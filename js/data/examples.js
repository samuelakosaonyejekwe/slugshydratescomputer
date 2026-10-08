// Worked example cases offered on the Case page. Each is a complete case file: a real producing province for
// the site (so that the Global site data page can pull its seabed, sea state and national figures) and a fluid
// typical of that kind of development. Suite inputs start from each suite's defaults and linked data.
import { FLUID_LIBRARY } from './fluids.js';
import { DEFAULT_FLUID } from '../core/thermo.js';

const fluid = (id, extra) => ({ ...DEFAULT_FLUID, ...FLUID_LIBRARY[id], comp: { ...FLUID_LIBRARY[id].comp }, kind: undefined, ...extra });
const mk = (name, notes, lat, lon, place, f) => ({ name, case: { schema: 1, name, notes, site: { name: place, lat, lon, country: '', countryCode: '', data: {}, fetchedAt: null }, fluid: f, inputs: {}, outputs: {}, autolink: true } });
export const EXAMPLES = [
  mk('Deep-water oil tie-back, Gulf of Guinea', 'Reference case: 18 km insulated 10-inch flowline at about 1,300 m water depth, steel catenary riser, 3,000 Sm³/d of light oil at 20 % water cut.', 3.5, 5.75, 'Offshore Niger Delta', fluid('lightOil', { qOil: 3000, wc: 20, Tin: 70, Pout: 25, Pres: 300, Tres: 90 })),
  mk('Ultra-deep long oil tie-back, Gulf of Mexico', 'Black oil from a 2,000 m deep field: long cold flowline, high hydrostatic head in the riser, cooldown and restart govern.', 27.3, -90.1, 'Green Canyon, Gulf of Mexico', fluid('blackOil', { qOil: 4500, wc: 10, Tin: 85, Pout: 30, Pres: 550, Tres: 110 })),
  mk('Gas-condensate subsea-to-shore, North-West Shelf', 'Wet gas with condensate and continuous MEG injection through a long trunk line to an onshore plant: liquid hold-up, turndown and slug-catcher sizing govern.', -19.6, 115.9, 'Carnarvon Basin, Australia', fluid('gasCondensate', { rateBasis: 'gas', qGas: 14, qWater: 120, inhibitor: 'MEG', inhWt: 45, Tin: 60, Pout: 70, Pres: 320, Tres: 105 })),
  mk('Wet-gas tie-back with methanol, northern North Sea', 'Cold shallow water, short tie-back to a host platform: hydrate margin at low rates and during shutdown governs.', 61.2, 1.8, 'Northern North Sea', fluid('wetGas', { rateBasis: 'gas', qGas: 5, qWater: 30, inhibitor: 'MeOH', inhWt: 20, Tin: 45, Pout: 60, Pres: 250, Tres: 95 })),
  mk('Pre-salt volatile oil with CO₂, Santos Basin', 'High-GOR oil with a high CO₂ content in 2,100 m of water: corrosion, wax and hydrate management together.', -25.3, -42.8, 'Santos Basin, Brazil', fluid('volatileOil', { comp: { ...FLUID_LIBRARY.volatileOil.comp, CO2: 12, C1: 50.4 }, qOil: 4000, wc: 15, Tin: 60, Pout: 28, Pres: 550, Tres: 65 })),
  mk('Onshore hilly-terrain multiphase gathering line', 'Buried onshore line over hills carrying oil, gas and water to a central facility: terrain slugging and liquid accumulation at low points govern.', 31.9, -102.1, 'Permian Basin, Texas', fluid('lightOil', { qOil: 1200, wc: 45, Tin: 55, Pout: 12, Pres: 180, Tres: 75 })),
  mk('Heavy waxy oil in late life, shallow water', 'Viscous oil at a high water cut and a low rate: wax deposition, sand transport and corrosion govern.', 4.3, 7.1, 'Shallow-water Niger Delta', fluid('heavyOil', { qOil: 900, wc: 65, Tin: 60, Pout: 10, Pres: 140, Tres: 80 })),
];
