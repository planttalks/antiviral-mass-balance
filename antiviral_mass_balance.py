"""
Mass-balance and antiviral metrics for pre-harvest plants spiked with surrogate viruses.

Equations follow the Materials and methods in Zure, Kuo, & Drizo, *Chemosphere*
(2024; DOI https://doi.org/10.1016/j.chemosphere.2023.141101;
PII https://www.sciencedirect.com/science/article/abs/pii/S0045653523033714): Eq.1–Eq.9.

Use one consistent infectivity or genome unit per virus (PFU, calibrated gene copies, or
Eq.9-derived Vp/ml) before closing a mass balance; mix units only after explicit conversion.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

__version__ = "0.1.0"
__all__ = [
    "viral_inactivation_percent",
    "initial_infectivity_pfu_ml",
    "final_infectivity_pfu_per_g",
    "final_combined_infectivity_weighted",
    "final_mixture_infectivity_pfu_per_g",
    "phytoremediation_efficiency_percent",
    "translocation_factor",
    "bioconcentration_factor",
    "highly_degraded_vp_per_ml",
    "natural_decay_subtracted_initial",
    "preharvest_mass_balance",
    "dilution_after_spike",
    "PlantPartLoads",
    "PreharvestMassBalance",
    "MS2_GENOME_BP",
    "MS2_MOLECULAR_WEIGHT_EQ9",
    "T4_GENOME_BP",
    "T4_MOLECULAR_WEIGHT_EQ9",
    "AVOGADRO",
]

# MS2 / T4 genome and mass constants from the manuscript (Eq.9 narrative)
MS2_GENOME_BP = 3569
# Manuscript average mass term for MS2 (inside Eq.9 as published)
MS2_MOLECULAR_WEIGHT_EQ9 = 330
T4_GENOME_BP = 169_000
# Manuscript average mass term for T4 (inside Eq.9 as published)
T4_MOLECULAR_WEIGHT_EQ9 = 660

AVOGADRO = 6.022e23


def viral_inactivation_percent(initial_infectivity: float, final_infectivity: float) -> float:
    """Eq.1: VI (%) = (initial - final) / initial * 100."""
    if initial_infectivity <= 0:
        raise ValueError("initial_infectivity must be positive")
    return (initial_infectivity - final_infectivity) / initial_infectivity * 100.0


def initial_infectivity_pfu_ml(
    plaque_count: float,
    dilution_factor: float,
    virus_diluted_volume_ml: float,
) -> float:
    """Eq.2: A (PFU/mL) = plaque_count / (dilution_factor * vol_virus_diluted)."""
    den = dilution_factor * virus_diluted_volume_ml
    if den <= 0:
        raise ValueError("dilution_factor and virus_diluted_volume_ml must be positive")
    return plaque_count / den


def final_infectivity_pfu_per_g(
    a_part_pfu_ml: float, crude_extract_volume_ml: float, fresh_weight_g: float
) -> float:
    """Eq.3: PI (PFU/g) = (A_part * vol_extract) / fwt_part."""
    if fresh_weight_g <= 0:
        raise ValueError("fresh_weight_g must be positive")
    return (a_part_pfu_ml * crude_extract_volume_ml) / fresh_weight_g


def final_combined_infectivity_weighted(
    pi_by_part_pfu_per_g: Mapping[str, float], mass_by_part_g: Mapping[str, float]
) -> float:
    """
    Mass-weighted mean PFU/g for an herbal combination (recommended form for manuscript Eq.4).

    Published Eq.4 repeats PI_root and PI_part/fwt_herb does not yield PFU/g; use
    CI = sum_i (m_i / M_herb) * PI_i with PI_i from Eq.3.
    """
    keys = set(pi_by_part_pfu_per_g) & set(mass_by_part_g)
    if not keys:
        raise ValueError("pi_by_part and mass_by_part must share at least one key")
    m_tot = sum(mass_by_part_g[k] for k in keys)
    if m_tot <= 0:
        raise ValueError("sum of masses must be positive")
    return sum((mass_by_part_g[k] / m_tot) * pi_by_part_pfu_per_g[k] for k in keys)


def final_mixture_infectivity_pfu_per_g(
    a_mixture_pfu_ml: float,
    crude_extract_volume_ml: float,
    sum_fwt_parts_g: float,
    fwt_whole_herb_g: float,
) -> float:
    """Eq.5: MI (PFU/g) = (A_mixture * vol_extract) / (sum fwt parts) * fwt_herb."""
    if sum_fwt_parts_g <= 0:
        raise ValueError("sum_fwt_parts_g must be positive")
    return (a_mixture_pfu_ml * crude_extract_volume_ml) / sum_fwt_parts_g * fwt_whole_herb_g


def phytoremediation_efficiency_percent(initial_particles: float, final_particles: float) -> float:
    """Eq.6: PE (%) = (initial - final) / initial * 100 (bulk water / culture)."""
    if initial_particles <= 0:
        raise ValueError("initial_particles must be positive")
    return (initial_particles - final_particles) / initial_particles * 100.0


def translocation_factor(viral_particles_target: float, viral_particles_root: float) -> float:
    """Eq.7: TF = particles_target / particles_root."""
    if viral_particles_root <= 0:
        raise ValueError("viral_particles_root must be positive")
    return viral_particles_target / viral_particles_root


def bioconcentration_factor(viral_particles_target: float, viral_particles_water: float) -> float:
    """Eq.8: BCF = particles_target / particles_water."""
    if viral_particles_water <= 0:
        raise ValueError("viral_particles_water must be positive")
    return viral_particles_target / viral_particles_water


def highly_degraded_vp_per_ml(
    dna_ng_per_ul: float, target_length_bp: float, molecular_weight_eq9: float
) -> float:
    """
    Eq.9 as published: HDp (Vp/ml) = (DNA (ng/ul) * 6.022e23) / (target_length (bp) * 10^9 * MW).

    Use manuscript constants for MW: 330 (MS2) and 660 (T4). Same symbols as in the paper.
    """
    if target_length_bp <= 0 or molecular_weight_eq9 <= 0:
        raise ValueError("target_length_bp and molecular_weight_eq9 must be positive")
    return (dna_ng_per_ul * AVOGADRO) / (target_length_bp * 1e9 * molecular_weight_eq9)


@dataclass(frozen=True)
class PlantPartLoads:
    """Per-part inventories in total units (e.g., PFU per part = PFU/g * mass_g)."""

    mass_g: float
    infective_ap: float
    inactivated_iap: float
    highly_degraded_hdp: float

    @property
    def total_accounted(self) -> float:
        """Sum of Ap + IAp + HDp in the same unit system."""
        return self.infective_ap + self.inactivated_iap + self.highly_degraded_hdp


@dataclass(frozen=True)
class PreharvestMassBalance:
    """
    Closed-culture mass balance after subtracting natural decay (manuscript: IAe from blank).

    N0_eff = initial total in culture after spike (single virus type, one unit system).
    At time t:
        N0_eff_adj = N0_eff - N_natural_decay   (manuscript: subtract IAe-style loss from blank)
        N_water = infective (+ optional non-infective) remaining in bulk water
        N_plant = sum_p (Ap_p + IAp_p + HDp_p) in same units
        N_residual = N0_eff_adj - N_water - N_plant   (undetected degradation / bias)

    closure_fraction and residual_fraction use N0_eff_adj as the denominator.
    """

    n0_eff_adj: float
    n_water: float
    n_plant_total: float
    n_residual: float

    @property
    def closure_fraction(self) -> float:
        """Fraction of adjusted spike recovered in water + plant (1 = perfect closure)."""
        if self.n0_eff_adj <= 0:
            raise ValueError("n0_eff_adj must be positive")
        return (self.n_water + self.n_plant_total) / self.n0_eff_adj

    @property
    def residual_fraction(self) -> float:
        """Fraction of adjusted spike unaccounted for (negative = over-recovery from assay bias)."""
        if self.n0_eff_adj <= 0:
            raise ValueError("n0_eff_adj must be positive")
        return self.n_residual / self.n0_eff_adj


def natural_decay_subtracted_initial(
    initial_infectivity_bulk: float,
    natural_decay_fraction_or_percent: float,
    as_percent: bool = True,
) -> float:
    """
    Manuscript: IAe from blank; subtract from initial before interpreting IAp and HDp in plant.

    If as_percent, pass VI_blank from Eq.1 style (percent). Else pass fractional loss in (0,1).
    """
    if as_percent:
        loss = initial_infectivity_bulk * (natural_decay_fraction_or_percent / 100.0)
    else:
        loss = initial_infectivity_bulk * natural_decay_fraction_or_percent
    return initial_infectivity_bulk - loss


def preharvest_mass_balance(
    culture_volume_ml: float,
    initial_concentration_per_ml: float,
    natural_decay_percent_from_blank: float,
    water_concentration_per_ml_at_t: float,
    parts: Sequence[PlantPartLoads],
) -> PreharvestMassBalance:
    """
    Build N0 from C0 * V, subtract blank decay, aggregate plant compartments.

    Concentrations must be in the same unit system as inventory (e.g., PFU/ml vs PFU total).
    Plant loads should be totals per part: concentration_per_g * mass_g for each category.
    """
    n0_eff = initial_concentration_per_ml * culture_volume_ml
    iae_pct = natural_decay_percent_from_blank
    n0_adj = natural_decay_subtracted_initial(n0_eff, iae_pct, as_percent=True)
    n_water = water_concentration_per_ml_at_t * culture_volume_ml
    n_plant = sum(p.total_accounted for p in parts)
    n_res = n0_adj - n_water - n_plant
    return PreharvestMassBalance(
        n0_eff_adj=n0_adj,
        n_water=n_water,
        n_plant_total=n_plant,
        n_residual=n_res,
    )


def dilution_after_spike(
    stock_concentration_per_ml: float, stock_volume_ml: float, diluent_volume_ml: float
) -> float:
    """Helper: C_initial in mixed culture after adding stock to water (e.g., 1 ml into 50 ml)."""
    v_tot = stock_volume_ml + diluent_volume_ml
    if v_tot <= 0:
        raise ValueError("total volume must be positive")
    return stock_concentration_per_ml * stock_volume_ml / v_tot
