// Reference reservoir fluids offered on the Case page: water-free well-stream compositions in mol % with the
// molar mass and specific gravity of the heptanes-plus fraction. Typical published compositions by fluid class.
const F = (name, kind, comp, c7MW, c7SG, extra = {}) => Object.freeze({ name, kind, comp: Object.freeze(comp), c7MW, c7SG, ...extra });
export const FLUID_LIBRARY = Object.freeze({
  lightOil: F('Light oil with associated gas', 'Black / light oil', { N2: 0.5, CO2: 2, H2S: 0, C1: 45, C2: 7, C3: 5, iC4: 1, nC4: 2.5, iC5: 1, nC5: 1.5, C6: 2.5, C7p: 32 }, 210, 0.84),
  blackOil: F('Black oil (low gas–oil ratio)', 'Black oil', { N2: 0.2, CO2: 0.9, H2S: 0, C1: 33, C2: 5.5, C3: 4.2, iC4: 0.9, nC4: 2.2, iC5: 1.1, nC5: 1.5, C6: 3.5, C7p: 47 }, 245, 0.875),
  heavyOil: F('Heavy oil', 'Heavy oil', { N2: 0.2, CO2: 0.5, H2S: 0, C1: 20, C2: 2.5, C3: 1.8, iC4: 0.5, nC4: 1.2, iC5: 0.8, nC5: 1, C6: 2.5, C7p: 69 }, 330, 0.93),
  volatileOil: F('Volatile oil', 'Volatile oil', { N2: 0.9, CO2: 2.4, H2S: 0, C1: 60, C2: 8, C3: 4.5, iC4: 1, nC4: 2, iC5: 0.9, nC5: 1.1, C6: 1.7, C7p: 17.5 }, 185, 0.815),
  gasCondensate: F('Gas condensate', 'Retrograde gas', { N2: 0.6, CO2: 2.5, H2S: 0, C1: 76, C2: 7.5, C3: 3.5, iC4: 0.7, nC4: 1.3, iC5: 0.6, nC5: 0.6, C6: 1, C7p: 5.7 }, 155, 0.79, { rateBasis: 'gas', qGas: 8, qWater: 60 }),
  wetGas: F('Wet gas', 'Wet gas', { N2: 0.8, CO2: 1.8, H2S: 0, C1: 86, C2: 5.5, C3: 2.4, iC4: 0.5, nC4: 0.8, iC5: 0.3, nC5: 0.3, C6: 0.4, C7p: 1.2 }, 125, 0.765, { rateBasis: 'gas', qGas: 10, qWater: 40 }),
  dryGas: F('Lean dry gas', 'Dry gas', { N2: 1.2, CO2: 1, H2S: 0, C1: 94, C2: 2.6, C3: 0.7, iC4: 0.15, nC4: 0.2, iC5: 0.05, nC5: 0.05, C6: 0.05, C7p: 0 }, 110, 0.75, { rateBasis: 'gas', qGas: 12, qWater: 20 }),
  sourGas: F('Sour CO₂-rich gas', 'Sour gas', { N2: 1.5, CO2: 9, H2S: 4, C1: 74, C2: 5.5, C3: 2.5, iC4: 0.5, nC4: 0.9, iC5: 0.3, nC5: 0.3, C6: 0.5, C7p: 1 }, 125, 0.765, { rateBasis: 'gas', qGas: 6, qWater: 50 }),
});
