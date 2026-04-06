# Antiviral mass-balance

[![CI](https://github.com/planttalks/antiviral-mass-balance/actions/workflows/ci.yml/badge.svg)](https://github.com/planttalks/antiviral-mass-balance/actions/workflows/ci.yml)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

Companion Python module for Eq.1-Eq.9 and pre-harvest mass-balance calculations described in:

Zure D, Kuo H-W D, Drizo A. *Insights of phytoremediation mechanisms for viruses based on in-vitro, in-vivo and in-silico assessments of selected herbal plants*. Chemosphere. 2024;351:141101. DOI: [10.1016/j.chemosphere.2023.141101](https://doi.org/10.1016/j.chemosphere.2023.141101).

## Scope

- Implements Eq.1-Eq.9 from the manuscript methods
- Supports infectivity metrics from plaque assay inputs
- Supports RT-PCR-derived particle metrics
- Supports HDp conversion and pre-harvest closure accounting

## Equation map

- Eq.1 -> `viral_inactivation_percent`
- Eq.2 -> `initial_infectivity_pfu_ml`
- Eq.3 -> `final_infectivity_pfu_per_g`
- Eq.4 -> `final_combined_infectivity_weighted` (mass-weighted implementation)
- Eq.5 -> `final_mixture_infectivity_pfu_per_g`
- Eq.6 -> `phytoremediation_efficiency_percent`
- Eq.7 -> `translocation_factor`
- Eq.8 -> `bioconcentration_factor`
- Eq.9 -> `highly_degraded_vp_per_ml`

## Required inputs

Three data sources feed the mass-balance. All come directly from the bench.

| Source | What you record | Code variable |
|---|---|---|
| **Culture setup** | Total vessel volume (mL) | `culture_volume_ml` |
| | Stock virus titer (PFU/mL) | `stock_concentration_per_ml` |
| | Volume of stock spiked in (mL) | `stock_volume_ml` |
| **Plaque assay** | Plaque count, dilution factor, inoculated volume per plate (blank, bulk water, each part) | `plaque_count`, `dilution_factor`, `virus_diluted_volume_ml` |
| | Fresh weight of each harvested part (g) | `fresh_weight_g` |
| | Volume of crude tissue extract (mL) | `crude_extract_volume_ml` |
| **RT-PCR + NanoDrop** | Gene copies/mL from bulk water at harvest | `viral_particles_water` |
| | Gene copies/mL from root and target-part extracts | `viral_particles_root`, `viral_particles_target` |
| | NanoDrop DNA reading from tissue extract (ng/uL) | `dna_ng_per_ul` |

Genome length and molecular weight constants for MS2 and T4 are built in (`MS2_GENOME_BP`, `MS2_MOLECULAR_WEIGHT_EQ9`, etc.) and require no measurement.

## Quick start

```python
import antiviral_mass_balance as amb

# 1. Initial concentration after spiking 1 mL stock into 50 mL water
c0 = amb.dilution_after_spike(stock_concentration_per_ml=2e8,
                               stock_volume_ml=1.0, diluent_volume_ml=50.0)

# 2. Blank decay from herb-free control (Eq.1)
blank_decay_pct = amb.viral_inactivation_percent(c0, c0_blank_at_t)

# 3. Infectivity per plant part from plaque counts (Eq.2 + Eq.3)
a_root  = amb.initial_infectivity_pfu_ml(plaque_count=36, dilution_factor=1e-3,
                                          virus_diluted_volume_ml=0.1)
pi_root = amb.final_infectivity_pfu_per_g(a_root, crude_extract_volume_ml=0.5,
                                           fresh_weight_g=1.2)

# 4. Phytoremediation efficiency, translocation and bioconcentration (Eq.6-Eq.8)
pe  = amb.phytoremediation_efficiency_percent(c0, c_water_at_t)
tf  = amb.translocation_factor(gene_copies_leaf, gene_copies_root)
bcf = amb.bioconcentration_factor(gene_copies_leaf, gene_copies_water)

# 5. Highly degraded viral particles from NanoDrop reading (Eq.9)
hdp = amb.highly_degraded_vp_per_ml(dna_ng_per_ul=1.8,
                                     target_length_bp=amb.MS2_GENOME_BP,
                                     molecular_weight_eq9=amb.MS2_MOLECULAR_WEIGHT_EQ9)

# 6. Close the mass balance
bal = amb.preharvest_mass_balance(
    culture_volume_ml=51.0,
    initial_concentration_per_ml=c0,
    natural_decay_percent_from_blank=blank_decay_pct,
    water_concentration_per_ml_at_t=c_water_at_t,
    parts=[
        amb.PlantPartLoads(mass_g=1.2, infective_ap=pi_root*1.2,
                           inactivated_iap=pi_root*1.2*0.35,
                           highly_degraded_hdp=hdp*0.5),
        # add one PlantPartLoads per harvested part
    ],
)
print(f"Closure: {bal.closure_fraction:.1%}  |  Residual: {bal.residual_fraction:.1%}")
```

## Minimal workflow

1. Spike aqueous medium with surrogate viruses and run blank control.
2. Harvest plant parts and record fresh mass.
3. Run plaque assay and compute infectivity metrics (Eq.1-Eq.5).
4. Run RT-PCR and compute PE/TF/BCF (Eq.6-Eq.8).
5. Estimate HDp from nucleic-acid mass (Eq.9).
6. Subtract blank decay and close the mass balance with `preharvest_mass_balance`.

See [`examples/ms2_t4_worked_example.py`](examples/ms2_t4_worked_example.py) for a fully annotated run using manuscript numbers, and `tests/test_antiviral_mass_balance.py` for unit tests.

## Install

```bash
pip install -e ".[dev]"
```

## Verify

```bash
ruff check .
ruff format --check .
pytest --cov=antiviral_mass_balance --cov-fail-under=85
```

## Citation

Please cite:

Zure D, Kuo H-W D, Drizo A. Insights of phytoremediation mechanisms for viruses based on in-vitro, in-vivo and in-silico assessments of selected herbal plants. *Chemosphere*. 2024;351:141101. doi:[10.1016/j.chemosphere.2023.141101](https://doi.org/10.1016/j.chemosphere.2023.141101).

Also cite the software release recorded in `CITATION.cff`. See [`CONTRIBUTING.md`](CONTRIBUTING.md) before opening a pull request.

## License

MIT (see `LICENSE`).
