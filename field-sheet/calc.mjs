/**
 * Eq. 1 to Eq. 9 and the pre-harvest closure.
 * Keep this file in step with antiviral_mass_balance.py.
 */

export const MS2_GENOME_BP = 3569;
export const MS2_MOLECULAR_WEIGHT_EQ9 = 330;
export const T4_GENOME_BP = 169_000;
export const T4_MOLECULAR_WEIGHT_EQ9 = 660;
export const AVOGADRO = 6.022e23;

export function viralInactivationPercent(initialInfectivity, finalInfectivity) {
  if (initialInfectivity <= 0) {
    throw new Error("initial_infectivity must be positive");
  }
  return ((initialInfectivity - finalInfectivity) / initialInfectivity) * 100.0;
}

export function initialInfectivityPfuMl(plaqueCount, dilutionFactor, virusDilutedVolumeMl) {
  const den = dilutionFactor * virusDilutedVolumeMl;
  if (den <= 0) {
    throw new Error("dilution_factor and virus_diluted_volume_ml must be positive");
  }
  return plaqueCount / den;
}

export function finalInfectivityPfuPerG(aPartPfuMl, crudeExtractVolumeMl, freshWeightG) {
  if (freshWeightG <= 0) {
    throw new Error("fresh_weight_g must be positive");
  }
  return (aPartPfuMl * crudeExtractVolumeMl) / freshWeightG;
}

export function finalCombinedInfectivityWeighted(piByPart, massByPart) {
  const keys = Object.keys(piByPart).filter((key) =>
    Object.prototype.hasOwnProperty.call(massByPart, key),
  );
  if (keys.length === 0) {
    throw new Error("pi_by_part and mass_by_part must share at least one key");
  }
  const massTotal = keys.reduce((sum, key) => sum + massByPart[key], 0);
  if (massTotal <= 0) {
    throw new Error("sum of masses must be positive");
  }
  return keys.reduce((sum, key) => sum + (massByPart[key] / massTotal) * piByPart[key], 0);
}

export function finalMixtureInfectivityPfuPerG(
  aMixturePfuMl,
  crudeExtractVolumeMl,
  sumFwtPartsG,
  fwtWholeHerbG,
) {
  if (sumFwtPartsG <= 0) {
    throw new Error("sum_fwt_parts_g must be positive");
  }
  return ((aMixturePfuMl * crudeExtractVolumeMl) / sumFwtPartsG) * fwtWholeHerbG;
}

export function phytoremediationEfficiencyPercent(initialParticles, finalParticles) {
  if (initialParticles <= 0) {
    throw new Error("initial_particles must be positive");
  }
  return ((initialParticles - finalParticles) / initialParticles) * 100.0;
}

export function translocationFactor(viralParticlesTarget, viralParticlesRoot) {
  if (viralParticlesRoot <= 0) {
    throw new Error("viral_particles_root must be positive");
  }
  return viralParticlesTarget / viralParticlesRoot;
}

export function bioconcentrationFactor(viralParticlesTarget, viralParticlesWater) {
  if (viralParticlesWater <= 0) {
    throw new Error("viral_particles_water must be positive");
  }
  return viralParticlesTarget / viralParticlesWater;
}

export function highlyDegradedVpPerMl(dnaNgPerUl, targetLengthBp, molecularWeightEq9) {
  if (targetLengthBp <= 0 || molecularWeightEq9 <= 0) {
    throw new Error("target_length_bp and molecular_weight_eq9 must be positive");
  }
  return (dnaNgPerUl * AVOGADRO) / (targetLengthBp * 1e9 * molecularWeightEq9);
}

export function naturalDecaySubtractedInitial(
  initialInfectivityBulk,
  naturalDecayFractionOrPercent,
  asPercent = true,
) {
  const loss = asPercent
    ? initialInfectivityBulk * (naturalDecayFractionOrPercent / 100.0)
    : initialInfectivityBulk * naturalDecayFractionOrPercent;
  return initialInfectivityBulk - loss;
}

export function dilutionAfterSpike(stockConcentrationPerMl, stockVolumeMl, diluentVolumeMl) {
  const total = stockVolumeMl + diluentVolumeMl;
  if (total <= 0) {
    throw new Error("total volume must be positive");
  }
  return (stockConcentrationPerMl * stockVolumeMl) / total;
}

export function preharvestMassBalance({
  cultureVolumeMl,
  initialConcentrationPerMl,
  naturalDecayPercentFromBlank,
  waterConcentrationPerMlAtT,
  parts,
}) {
  const n0Eff = initialConcentrationPerMl * cultureVolumeMl;
  const n0EffAdj = naturalDecaySubtractedInitial(n0Eff, naturalDecayPercentFromBlank, true);
  const nWater = waterConcentrationPerMlAtT * cultureVolumeMl;
  const nPlantTotal = parts.reduce(
    (sum, part) => sum + part.infectiveAp + part.inactivatedIap + part.highlyDegradedHdp,
    0,
  );
  const nResidual = n0EffAdj - nWater - nPlantTotal;
  return {
    n0EffAdj,
    nWater,
    nPlantTotal,
    nResidual,
    get closureFraction() {
      if (n0EffAdj <= 0) {
        throw new Error("n0_eff_adj must be positive");
      }
      return (nWater + nPlantTotal) / n0EffAdj;
    },
    get residualFraction() {
      if (n0EffAdj <= 0) {
        throw new Error("n0_eff_adj must be positive");
      }
      return nResidual / n0EffAdj;
    },
  };
}

export function readNumber(value) {
  if (value === null || value === undefined) return null;
  if (typeof value === "number") return Number.isFinite(value) ? value : Number.NaN;
  const text = String(value).trim();
  if (text === "") return null;
  const parsed = Number(text);
  return Number.isFinite(parsed) ? parsed : Number.NaN;
}

function invalid(errors, id, message) {
  errors.push({ id, message });
}

function requireFinite(errors, id, value, label) {
  const parsed = readNumber(value);
  if (parsed === null) return null;
  if (Number.isNaN(parsed)) {
    invalid(errors, id, `${label} is not a number.`);
    return null;
  }
  return parsed;
}

function partIsTouched(part) {
  const fields = [
    part.massG,
    part.plaqueCount,
    part.dilutionFactor,
    part.inoculumMl,
    part.extractMl,
    part.infectiveTotal,
    part.inactivatedTotal,
    part.inactivatedFraction,
    part.dnaNgPerUl,
  ];
  return fields.some((value) => readNumber(value) !== null);
}

function resolvePart(part, index, genome, errors) {
  const id = `part-${index}`;
  const name = String(part.name || "").trim() || `Part ${index + 1}`;
  const row = {
    name,
    massG: null,
    a: null,
    pi: null,
    infective: null,
    inactivated: null,
    hdpPerMl: null,
    hdpTotal: null,
    includedHdp: 0,
    used: false,
  };
  if (!partIsTouched(part)) return row;

  row.used = true;
  row.massG = requireFinite(errors, `${id}-mass`, part.massG, `${name} fresh weight`);
  const plaqueCount = requireFinite(errors, `${id}-plaque`, part.plaqueCount, `${name} plaque count`);
  const dilution = requireFinite(errors, `${id}-dilution`, part.dilutionFactor, `${name} dilution`);
  const inoculum = requireFinite(errors, `${id}-inoculum`, part.inoculumMl, `${name} inoculum`);
  const extract = requireFinite(errors, `${id}-extract`, part.extractMl, `${name} extract volume`);
  const manualInfective = requireFinite(
    errors,
    `${id}-infective`,
    part.infectiveTotal,
    `${name} infective total`,
  );
  const manualInactivated = requireFinite(
    errors,
    `${id}-inactivated`,
    part.inactivatedTotal,
    `${name} inactivated total`,
  );
  const fraction = requireFinite(
    errors,
    `${id}-fraction`,
    part.inactivatedFraction,
    `${name} inactivated fraction`,
  );
  const dna = requireFinite(errors, `${id}-dna`, part.dnaNgPerUl, `${name} DNA`);

  const plaqueFields = [plaqueCount, dilution, inoculum];
  const plaqueStarted = plaqueFields.some((value) => value !== null);
  const plaqueComplete = plaqueFields.every((value) => value !== null);

  if (plaqueStarted && !plaqueComplete) {
    invalid(errors, `${id}-plaque`, `${name} needs a plaque count, a dilution and an inoculum volume.`);
  }

  if (plaqueComplete) {
    if (row.massG === null || extract === null) {
      invalid(errors, `${id}-mass`, `${name} needs a fresh weight and an extract volume to turn plaques into a total.`);
    } else {
      try {
        row.a = initialInfectivityPfuMl(plaqueCount, dilution, inoculum);
        row.pi = finalInfectivityPfuPerG(row.a, extract, row.massG);
        row.infective = row.pi * row.massG;
      } catch (caught) {
        invalid(errors, `${id}-plaque`, caught instanceof Error ? caught.message : String(caught));
      }
    }
  } else if (manualInfective !== null) {
    row.infective = manualInfective;
    if (row.massG !== null && row.massG > 0) {
      row.pi = manualInfective / row.massG;
    }
  }

  if (row.used && row.infective === null && !plaqueStarted) {
    invalid(errors, `${id}-infective`, `${name} has no infective total.`);
  }

  if (fraction !== null && fraction < 0) {
    invalid(errors, `${id}-fraction`, `${name} inactivated fraction cannot be negative.`);
  }
  if (manualInactivated !== null && manualInactivated < 0) {
    invalid(errors, `${id}-inactivated`, `${name} inactivated total cannot be negative.`);
  }

  if (manualInactivated !== null) {
    row.inactivated = manualInactivated;
  } else if (fraction !== null && row.infective !== null && fraction >= 0) {
    row.inactivated = row.infective * fraction;
  } else if (row.infective !== null) {
    row.inactivated = 0;
  }

  if (dna !== null) {
    if (genome.lengthBp === null || genome.molecularWeight === null) {
      invalid(errors, "genomeLength", "Enter a genome length and a mass term before using a NanoDrop reading.");
    } else {
      try {
        row.hdpPerMl = highlyDegradedVpPerMl(dna, genome.lengthBp, genome.molecularWeight);
        const volume = extract === null ? null : extract;
        if (volume === null) {
          invalid(errors, `${id}-extract`, `${name} needs an extract volume to turn HDp per mL into a total.`);
        } else {
          row.hdpTotal = row.hdpPerMl * volume;
        }
      } catch (caught) {
        invalid(errors, `${id}-dna`, caught instanceof Error ? caught.message : String(caught));
      }
    }
  }

  if (part.includeHdp) {
    if (row.hdpTotal === null) {
      invalid(errors, `${id}-dna`, `${name} is marked to enter the closure, but HDp has no total yet.`);
    } else {
      row.includedHdp = row.hdpTotal;
    }
  }

  return row;
}

export function computeSheet(input) {
  const errors = [];
  const stockConcentration = requireFinite(
    errors,
    "stockConcentration",
    input.stockConcentration,
    "Stock concentration",
  );
  const stockVolume = requireFinite(errors, "stockVolumeMl", input.stockVolumeMl, "Stock volume");
  const diluentVolume = requireFinite(errors, "diluentVolumeMl", input.diluentVolumeMl, "Diluent volume");
  const cultureOverride = requireFinite(
    errors,
    "cultureVolumeMl",
    input.cultureVolumeMl,
    "Culture volume",
  );

  let cultureVolumeMl = null;
  if (cultureOverride !== null) {
    cultureVolumeMl = cultureOverride;
  } else if (stockVolume !== null && diluentVolume !== null) {
    cultureVolumeMl = stockVolume + diluentVolume;
  }

  let c0 = null;
  if (stockConcentration !== null && stockVolume !== null && diluentVolume !== null) {
    try {
      c0 = dilutionAfterSpike(stockConcentration, stockVolume, diluentVolume);
    } catch (caught) {
      invalid(errors, "stockVolumeMl", caught instanceof Error ? caught.message : String(caught));
    }
  }

  let blankDecayPercent = null;
  if (input.blankMode === "pair") {
    const initial = requireFinite(errors, "blankInitial", input.blankInitial, "Blank initial titer");
    const finalTiter = requireFinite(errors, "blankFinal", input.blankFinal, "Blank final titer");
    if (initial !== null && finalTiter !== null) {
      try {
        blankDecayPercent = viralInactivationPercent(initial, finalTiter);
      } catch (caught) {
        invalid(errors, "blankInitial", caught instanceof Error ? caught.message : String(caught));
      }
    } else if (initial !== null || finalTiter !== null) {
      invalid(errors, "blankInitial", "Blank titers need both the start and the end.");
    }
  } else {
    blankDecayPercent = requireFinite(
      errors,
      "blankDecayPercent",
      input.blankDecayPercent,
      "Blank decay",
    );
  }

  const water = requireFinite(
    errors,
    "waterConcentration",
    input.waterConcentrationPerMl,
    "Water concentration",
  );

  let pePercent = null;
  if (c0 !== null && water !== null) {
    try {
      pePercent = phytoremediationEfficiencyPercent(c0, water);
    } catch (caught) {
      invalid(errors, "waterConcentration", caught instanceof Error ? caught.message : String(caught));
    }
  }

  const genomeLength = requireFinite(errors, "genomeLength", input.genomeLengthBp, "Genome length");
  const molecularWeight = requireFinite(
    errors,
    "molecularWeight",
    input.molecularWeight,
    "Molecular weight term",
  );
  const genome = { lengthBp: genomeLength, molecularWeight };

  const parts = (input.parts || []).map((part, index) => resolvePart(part, index, genome, errors));

  const weightedRows = parts.filter((part) => part.pi !== null && part.massG !== null && part.massG > 0);
  let combinedPi = null;
  if (weightedRows.length > 0) {
    const piByPart = {};
    const massByPart = {};
    weightedRows.forEach((part, index) => {
      piByPart[`row-${index}`] = part.pi;
      massByPart[`row-${index}`] = part.massG;
    });
    combinedPi = finalCombinedInfectivityWeighted(piByPart, massByPart);
  }

  const mixtureA = requireFinite(errors, "mixtureA", input.mixtureA, "Mixture titer");
  const mixtureExtract = requireFinite(
    errors,
    "mixtureExtract",
    input.mixtureExtractMl,
    "Mixture extract volume",
  );
  const mixtureSum = requireFinite(errors, "mixtureSum", input.mixtureSumFwtG, "Sum of fresh weights");
  const mixtureHerb = requireFinite(
    errors,
    "mixtureHerb",
    input.mixtureWholeHerbG,
    "Whole herb fresh weight",
  );
  const mixtureStarted = [mixtureA, mixtureExtract, mixtureSum, mixtureHerb].some((value) => value !== null);
  let mixturePi = null;
  if (mixtureStarted && [mixtureA, mixtureExtract, mixtureSum, mixtureHerb].every((value) => value !== null)) {
    try {
      mixturePi = finalMixtureInfectivityPfuPerG(mixtureA, mixtureExtract, mixtureSum, mixtureHerb);
    } catch (caught) {
      invalid(errors, "mixtureSum", caught instanceof Error ? caught.message : String(caught));
    }
  } else if (mixtureStarted) {
    invalid(errors, "mixtureA", "The mixture block needs all four values.");
  }

  const target = requireFinite(errors, "targetParticles", input.targetParticles, "Target particles");
  const rootParticles = requireFinite(errors, "rootParticles", input.rootParticles, "Root particles");
  const waterParticles = requireFinite(errors, "waterParticles", input.waterParticles, "Water particles");
  let tf = null;
  if (target !== null && rootParticles !== null) {
    try {
      tf = translocationFactor(target, rootParticles);
    } catch (caught) {
      invalid(errors, "rootParticles", caught instanceof Error ? caught.message : String(caught));
    }
  }
  let bcf = null;
  if (target !== null && waterParticles !== null) {
    try {
      bcf = bioconcentrationFactor(target, waterParticles);
    } catch (caught) {
      invalid(errors, "waterParticles", caught instanceof Error ? caught.message : String(caught));
    }
  }

  const usedParts = parts.filter((part) => part.used);
  const partErrors = errors.some(
    (error) => error.id.startsWith("part-") || error.id === "genomeLength",
  );
  let balance = null;
  if (c0 !== null && cultureVolumeMl !== null && blankDecayPercent !== null && water !== null && !partErrors) {
    const ready = usedParts.every((part) => part.infective !== null && part.inactivated !== null);
    if (ready) {
      try {
        balance = preharvestMassBalance({
          cultureVolumeMl,
          initialConcentrationPerMl: c0,
          naturalDecayPercentFromBlank: blankDecayPercent,
          waterConcentrationPerMlAtT: water,
          parts: usedParts.map((part) => ({
            infectiveAp: part.infective,
            inactivatedIap: part.inactivated,
            highlyDegradedHdp: part.includedHdp,
          })),
        });
        if (balance.n0EffAdj <= 0) {
          balance = null;
          invalid(errors, "blankDecayPercent", "Adjusted spike is not positive. Check the blank decay.");
        }
      } catch (caught) {
        invalid(errors, "blankDecayPercent", caught instanceof Error ? caught.message : String(caught));
      }
    }
  }

  const missing = [];
  if (c0 === null) missing.push("spike");
  if (blankDecayPercent === null) missing.push("blank");
  if (water === null) missing.push("water");

  return {
    errors,
    missing,
    cultureVolumeMl,
    c0,
    blankDecayPercent,
    pePercent,
    tf,
    bcf,
    parts,
    combinedPi,
    mixturePi,
    balance,
    hdpIncluded: usedParts.some((part) => part.includedHdp !== 0),
  };
}
