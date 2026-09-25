import { dilutionAfterSpike } from "./calc.mjs";

const EXTRACT_ML = 0.5;

function partLoad(name, massG, waterFactor, inactivatedFraction, dnaScale, waterConcentration) {
  return {
    name,
    massG,
    infectiveTotal: waterConcentration * waterFactor * EXTRACT_ML,
    inactivatedFraction,
    extractMl: EXTRACT_ML,
    dnaNgPerUl: 1.8 * dnaScale,
    includeHdp: false,
  };
}

/** Illustrative Ocimum basilicum sheet. Part totals follow examples/ms2_t4_worked_example.py. */
export function ms2WorkedExample() {
  const water = dilutionAfterSpike(2e8, 1, 50) * (1 - 0.646);
  return {
    note: "Illustrative MS2 sheet for Ocimum basilicum at 4 h. Part totals follow the worked example. HDp is shown and is not in the closure.",
    input: {
      herb: "Ocimum basilicum",
      virus: "ms2",
      unitSystem: "PFU",
      hours: 4,
      notes:
        "Illustrative MS2 sheet. Part totals follow examples/ms2_t4_worked_example.py. HDp is shown and is not in the closure.",
      stockConcentration: 2e8,
      stockVolumeMl: 1,
      diluentVolumeMl: 50,
      blankMode: "percent",
      blankDecayPercent: 2.3,
      waterConcentrationPerMl: water,
      genomeLengthBp: 3569,
      molecularWeight: 330,
      targetParticles: 4.5e5,
      rootParticles: 3e5,
      waterParticles: water,
      parts: [
        partLoad("root", 1.2, 0.3, 0.35, 1, water),
        partLoad("leaf", 2.1, 0.45, 0.56, 1.4, water),
        partLoad("shoot", 0.9, 0.12, 0.3, 0.5, water),
      ],
    },
  };
}

/** T4 spike from the worked example. Water is back-calculated from the reported 65.8% removal. */
export function t4ReportedWater() {
  const water = dilutionAfterSpike(2.37e7, 1, 50) * (1 - 0.658);
  return {
    note: "T4 spike from the worked example. The water titer is back-calculated from the reported 65.8% removal. Parts are empty.",
    input: {
      virus: "t4",
      unitSystem: "PFU",
      hours: 4,
      notes:
        "T4 spike. Water titer is back-calculated from the reported 65.8% removal. Parts are empty.",
      stockConcentration: 2.37e7,
      stockVolumeMl: 1,
      diluentVolumeMl: 50,
      blankMode: "percent",
      blankDecayPercent: 2.3,
      waterConcentrationPerMl: water,
      genomeLengthBp: 169000,
      molecularWeight: 660,
      parts: [{ name: "root" }, { name: "leaf" }, { name: "shoot" }],
    },
  };
}
