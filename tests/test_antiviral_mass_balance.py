"""Tests for antiviral_mass_balance (manuscript Eq.1–Eq.9 and closure helpers)."""

from __future__ import annotations

import pytest

import antiviral_mass_balance as amb


def test_viral_inactivation_percent() -> None:
    assert amb.viral_inactivation_percent(100.0, 30.0) == pytest.approx(70.0)
    with pytest.raises(ValueError):
        amb.viral_inactivation_percent(0.0, 0.0)


def test_initial_infectivity_pfu_ml() -> None:
    # plaque_count / (dilution * vol)
    assert amb.initial_infectivity_pfu_ml(60.0, 10.0, 0.1) == pytest.approx(60.0)
    with pytest.raises(ValueError):
        amb.initial_infectivity_pfu_ml(1.0, 0.0, 1.0)


def test_final_infectivity_pfu_per_g() -> None:
    assert amb.final_infectivity_pfu_per_g(1e6, 0.5, 2.0) == pytest.approx(250_000.0)
    with pytest.raises(ValueError):
        amb.final_infectivity_pfu_per_g(1.0, 1.0, 0.0)


def test_final_combined_infectivity_weighted() -> None:
    pi = {"root": 100.0, "leaf": 300.0}
    m = {"root": 3.0, "leaf": 1.0}
    # (3/4)*100 + (1/4)*300 = 75 + 75 = 150
    assert amb.final_combined_infectivity_weighted(pi, m) == pytest.approx(150.0)
    with pytest.raises(ValueError):
        amb.final_combined_infectivity_weighted({"a": 1.0}, {"b": 1.0})


def test_final_mixture_infectivity_pfu_per_g() -> None:
    out = amb.final_mixture_infectivity_pfu_per_g(1e5, 1.0, 4.0, 2.0)
    assert out == pytest.approx(50_000.0)


def test_phytoremediation_efficiency_percent() -> None:
    assert amb.phytoremediation_efficiency_percent(100.0, 40.0) == pytest.approx(60.0)


def test_tf_bcf() -> None:
    assert amb.translocation_factor(20.0, 10.0) == pytest.approx(2.0)
    assert amb.bioconcentration_factor(5.0, 10.0) == pytest.approx(0.5)


def test_highly_degraded_vp_per_ml_ms2_shape() -> None:
    v = amb.highly_degraded_vp_per_ml(1.0, amb.MS2_GENOME_BP, amb.MS2_MOLECULAR_WEIGHT_EQ9)
    assert v > 0
    with pytest.raises(ValueError):
        amb.highly_degraded_vp_per_ml(1.0, 0.0, 330.0)


def test_natural_decay_subtracted_initial() -> None:
    assert amb.natural_decay_subtracted_initial(100.0, 2.3, as_percent=True) == pytest.approx(97.7)
    adj = amb.natural_decay_subtracted_initial(100.0, 0.05, as_percent=False)
    assert adj == pytest.approx(95.0)


def test_dilution_after_spike() -> None:
    c = amb.dilution_after_spike(2e8, 1.0, 50.0)
    assert c == pytest.approx(2e8 / 51.0)


def test_preharvest_mass_balance_closure() -> None:
    c0 = amb.dilution_after_spike(2e8, 1.0, 50.0)
    v_ml = 51.0
    n0 = c0 * v_ml
    n0_adj = amb.natural_decay_subtracted_initial(n0, 0.0, as_percent=True)
    # Perfect bookkeeping: half in water, half in one plant part
    half = n0_adj / 2.0
    parts = [
        amb.PlantPartLoads(
            mass_g=1.0,
            infective_ap=half,
            inactivated_iap=0.0,
            highly_degraded_hdp=0.0,
        )
    ]
    bal = amb.preharvest_mass_balance(
        culture_volume_ml=v_ml,
        initial_concentration_per_ml=c0,
        natural_decay_percent_from_blank=0.0,
        water_concentration_per_ml_at_t=half / v_ml,
        parts=parts,
    )
    assert bal.n0_eff_adj == pytest.approx(n0_adj)
    assert bal.closure_fraction == pytest.approx(1.0)
    assert bal.n_residual == pytest.approx(0.0, abs=1e-6)


def test_preharvest_mass_balance_closure_fraction_raises_on_zero_n0_adj() -> None:
    bal = amb.PreharvestMassBalance(n0_eff_adj=0.0, n_water=0.0, n_plant_total=0.0, n_residual=0.0)
    with pytest.raises(ValueError):
        _ = bal.closure_fraction
