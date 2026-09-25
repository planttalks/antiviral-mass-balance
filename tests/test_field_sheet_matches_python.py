"""The browser sheet must stay on the same equations as the Python module."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
from pathlib import Path

import pytest

import antiviral_mass_balance as amb

ROOT = Path(__file__).resolve().parents[1]


def test_field_sheet_js_matches_python() -> None:
    node = shutil.which("node")
    if node is None:
        if os.environ.get("CI") == "true" and os.environ.get("FIELD_SHEET_JS") == "1":
            pytest.fail("node is missing in CI")
        pytest.skip("node is required to compare the field sheet")

    script = """
import * as sheet from "./field-sheet/calc.mjs";
const c0 = sheet.dilutionAfterSpike(2e8, 1, 50);
const volume = 51;
const half = (c0 * volume) / 2;
const balance = sheet.preharvestMassBalance({
  cultureVolumeMl: volume,
  initialConcentrationPerMl: c0,
  naturalDecayPercentFromBlank: 0,
  waterConcentrationPerMlAtT: half / volume,
  parts: [{ infectiveAp: half, inactivatedIap: 0, highlyDegradedHdp: 0 }],
});
console.log(JSON.stringify({
  vi: sheet.viralInactivationPercent(100, 30),
  a: sheet.initialInfectivityPfuMl(60, 10, 0.1),
  pi: sheet.finalInfectivityPfuPerG(1e6, 0.5, 2),
  ci: sheet.finalCombinedInfectivityWeighted({ root: 100, leaf: 300 }, { root: 3, leaf: 1 }),
  mi: sheet.finalMixtureInfectivityPfuPerG(1e5, 1, 4, 2),
  pe: sheet.phytoremediationEfficiencyPercent(100, 40),
  tf: sheet.translocationFactor(20, 10),
  bcf: sheet.bioconcentrationFactor(5, 10),
  hdp: sheet.highlyDegradedVpPerMl(1, sheet.MS2_GENOME_BP, sheet.MS2_MOLECULAR_WEIGHT_EQ9),
  decay: sheet.naturalDecaySubtractedInitial(100, 2.3, true),
  c0,
  ms2: [sheet.MS2_GENOME_BP, sheet.MS2_MOLECULAR_WEIGHT_EQ9],
  t4: [sheet.T4_GENOME_BP, sheet.T4_MOLECULAR_WEIGHT_EQ9],
  closure: balance.closureFraction,
  residual: balance.nResidual,
}));
"""
    completed = subprocess.run(
        [node, "--input-type=module", "-e", script],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    got = json.loads(completed.stdout)
    assert got["vi"] == pytest.approx(amb.viral_inactivation_percent(100.0, 30.0))
    assert got["a"] == pytest.approx(amb.initial_infectivity_pfu_ml(60.0, 10.0, 0.1))
    assert got["pi"] == pytest.approx(amb.final_infectivity_pfu_per_g(1e6, 0.5, 2.0))
    combined = amb.final_combined_infectivity_weighted(
        {"root": 100.0, "leaf": 300.0},
        {"root": 3.0, "leaf": 1.0},
    )
    assert got["ci"] == pytest.approx(combined)
    assert got["mi"] == pytest.approx(amb.final_mixture_infectivity_pfu_per_g(1e5, 1.0, 4.0, 2.0))
    assert got["pe"] == pytest.approx(amb.phytoremediation_efficiency_percent(100.0, 40.0))
    assert got["tf"] == pytest.approx(amb.translocation_factor(20.0, 10.0))
    assert got["bcf"] == pytest.approx(amb.bioconcentration_factor(5.0, 10.0))
    assert got["hdp"] == pytest.approx(
        amb.highly_degraded_vp_per_ml(1.0, amb.MS2_GENOME_BP, amb.MS2_MOLECULAR_WEIGHT_EQ9)
    )
    decay = amb.natural_decay_subtracted_initial(100.0, 2.3, as_percent=True)
    assert got["decay"] == pytest.approx(decay)
    assert got["c0"] == pytest.approx(amb.dilution_after_spike(2e8, 1.0, 50.0))
    assert got["ms2"] == [amb.MS2_GENOME_BP, amb.MS2_MOLECULAR_WEIGHT_EQ9]
    assert got["t4"] == [amb.T4_GENOME_BP, amb.T4_MOLECULAR_WEIGHT_EQ9]
    assert got["closure"] == pytest.approx(1.0)
    assert got["residual"] == pytest.approx(0.0, abs=1e-6)
