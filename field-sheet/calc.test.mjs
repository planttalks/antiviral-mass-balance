import assert from "node:assert/strict";
import test from "node:test";

import * as amb from "./calc.mjs";
import { ms2WorkedExample, t4ReportedWater } from "./examples.mjs";

test("viral inactivation percent", () => {
  assert.ok(Math.abs(amb.viralInactivationPercent(100, 30) - 70) < 1e-9);
  assert.throws(() => amb.viralInactivationPercent(0, 0), /initial_infectivity must be positive/);
});

test("initial infectivity", () => {
  assert.ok(Math.abs(amb.initialInfectivityPfuMl(60, 10, 0.1) - 60) < 1e-9);
  assert.throws(
    () => amb.initialInfectivityPfuMl(1, 0, 1),
    /dilution_factor and virus_diluted_volume_ml must be positive/,
  );
});

test("final infectivity per gram", () => {
  assert.ok(Math.abs(amb.finalInfectivityPfuPerG(1e6, 0.5, 2) - 250_000) < 1e-6);
  assert.throws(() => amb.finalInfectivityPfuPerG(1, 1, 0), /fresh_weight_g must be positive/);
});

test("combined infectivity is mass weighted", () => {
  const value = amb.finalCombinedInfectivityWeighted({ root: 100, leaf: 300 }, { root: 3, leaf: 1 });
  assert.ok(Math.abs(value - 150) < 1e-9);
  assert.throws(
    () => amb.finalCombinedInfectivityWeighted({ a: 1 }, { b: 1 }),
    /share at least one key/,
  );
});

test("mixture infectivity", () => {
  assert.ok(Math.abs(amb.finalMixtureInfectivityPfuPerG(1e5, 1, 4, 2) - 50_000) < 1e-6);
});

test("phytoremediation efficiency and ratios", () => {
  assert.ok(Math.abs(amb.phytoremediationEfficiencyPercent(100, 40) - 60) < 1e-9);
  assert.ok(Math.abs(amb.translocationFactor(20, 10) - 2) < 1e-9);
  assert.ok(Math.abs(amb.bioconcentrationFactor(5, 10) - 0.5) < 1e-9);
});

test("highly degraded MS2 shape", () => {
  const value = amb.highlyDegradedVpPerMl(1, amb.MS2_GENOME_BP, amb.MS2_MOLECULAR_WEIGHT_EQ9);
  assert.ok(value > 0);
  assert.throws(
    () => amb.highlyDegradedVpPerMl(1, 0, 330),
    /target_length_bp and molecular_weight_eq9 must be positive/,
  );
});

test("natural decay and dilution", () => {
  assert.ok(Math.abs(amb.naturalDecaySubtractedInitial(100, 2.3, true) - 97.7) < 1e-9);
  assert.ok(Math.abs(amb.naturalDecaySubtractedInitial(100, 0.05, false) - 95) < 1e-9);
  assert.ok(Math.abs(amb.dilutionAfterSpike(2e8, 1, 50) - 2e8 / 51) < 1e-6);
});

test("perfect closure", () => {
  const c0 = amb.dilutionAfterSpike(2e8, 1, 50);
  const volume = 51;
  const n0 = c0 * volume;
  const half = n0 / 2;
  const balance = amb.preharvestMassBalance({
    cultureVolumeMl: volume,
    initialConcentrationPerMl: c0,
    naturalDecayPercentFromBlank: 0,
    waterConcentrationPerMlAtT: half / volume,
    parts: [{ infectiveAp: half, inactivatedIap: 0, highlyDegradedHdp: 0 }],
  });
  assert.ok(Math.abs(balance.closureFraction - 1) < 1e-9);
  assert.ok(Math.abs(balance.nResidual) < 1e-6);
});

test("closure fraction rejects a zero adjusted spike", () => {
  const balance = amb.preharvestMassBalance({
    cultureVolumeMl: 1,
    initialConcentrationPerMl: 0,
    naturalDecayPercentFromBlank: 0,
    waterConcentrationPerMlAtT: 0,
    parts: [],
  });
  assert.throws(() => balance.closureFraction, /n0_eff_adj must be positive/);
});

test("sheet closure keeps HDp out until it is marked", () => {
  const sheet = amb.computeSheet({
    stockConcentration: 2e8,
    stockVolumeMl: 1,
    diluentVolumeMl: 50,
    blankMode: "percent",
    blankDecayPercent: 0,
    waterConcentrationPerMl: 1e6,
    genomeLengthBp: amb.MS2_GENOME_BP,
    molecularWeight: amb.MS2_MOLECULAR_WEIGHT_EQ9,
    parts: [
      {
        name: "root",
        massG: 1.2,
        plaqueCount: 36,
        dilutionFactor: 1e-3,
        inoculumMl: 0.1,
        extractMl: 0.5,
        dnaNgPerUl: 1.8,
        includeHdp: false,
      },
    ],
  });
  assert.equal(sheet.errors.length, 0);
  assert.equal(sheet.hdpIncluded, false);
  assert.ok(sheet.parts[0].hdpPerMl > 0);
  assert.equal(sheet.parts[0].includedHdp, 0);
  const infective = sheet.parts[0].infective;
  const expectedPlant = infective;
  assert.ok(Math.abs(sheet.balance.nPlantTotal - expectedPlant) < 1e-6);
});

test("marked HDp enters the plant total", () => {
  const sheet = amb.computeSheet({
    stockConcentration: 100,
    stockVolumeMl: 1,
    diluentVolumeMl: 1,
    blankMode: "percent",
    blankDecayPercent: 0,
    waterConcentrationPerMl: 10,
    genomeLengthBp: 1000,
    molecularWeight: 100,
    parts: [
      {
        name: "leaf",
        massG: 2,
        infectiveTotal: 50,
        inactivatedFraction: 0.5,
        extractMl: 2,
        dnaNgPerUl: 1,
        includeHdp: true,
      },
    ],
  });
  assert.equal(sheet.errors.length, 0);
  const hdp = amb.highlyDegradedVpPerMl(1, 1000, 100) * 2;
  assert.ok(Math.abs(sheet.balance.nPlantTotal - (50 + 25 + hdp)) < 1e-3);
  assert.equal(sheet.hdpIncluded, true);
});

test("a half-filled part blocks closure", () => {
  const sheet = amb.computeSheet({
    stockConcentration: 100,
    stockVolumeMl: 1,
    diluentVolumeMl: 9,
    blankMode: "percent",
    blankDecayPercent: 0,
    waterConcentrationPerMl: 1,
    parts: [{ name: "root", massG: 1 }],
  });
  assert.equal(sheet.balance, null);
  assert.ok(sheet.errors.some((error) => error.id === "part-0-infective"));
});

test("blank pair matches Eq. 1", () => {
  const sheet = amb.computeSheet({
    blankMode: "pair",
    blankInitial: 100,
    blankFinal: 30,
    parts: [],
  });
  assert.ok(Math.abs(sheet.blankDecayPercent - 70) < 1e-9);
});

test("MS2 worked example closes without HDp", () => {
  const sheet = amb.computeSheet(ms2WorkedExample().input);
  assert.equal(sheet.errors.length, 0);
  assert.equal(sheet.hdpIncluded, false);
  assert.ok(Math.abs(sheet.pePercent - 64.6) < 1e-9);
  assert.ok(Math.abs(sheet.tf - 1.5) < 1e-9);
  assert.ok(Math.abs(sheet.balance.closureFraction - 0.36682021795412123) < 1e-9);
  assert.ok(sheet.parts[0].hdpPerMl > 1e8);
});

test("T4 reported removal closes on water alone", () => {
  const sheet = amb.computeSheet(t4ReportedWater().input);
  assert.equal(sheet.errors.length, 0);
  assert.ok(Math.abs(sheet.pePercent - 65.8) < 1e-9);
  assert.ok(Math.abs(sheet.balance.closureFraction - 0.342 / 0.977) < 1e-9);
  assert.equal(sheet.balance.nPlantTotal, 0);
});

test("culture volume override replaces stock plus diluent", () => {
  const sheet = amb.computeSheet({
    stockConcentration: 10,
    stockVolumeMl: 1,
    diluentVolumeMl: 1,
    cultureVolumeMl: 8,
    blankMode: "percent",
    blankDecayPercent: 0,
    waterConcentrationPerMl: 1,
    parts: [],
  });
  assert.equal(sheet.cultureVolumeMl, 8);
  assert.equal(sheet.c0, 5);
  assert.equal(sheet.balance.nWater, 8);
});
