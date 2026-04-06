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

## Minimal workflow

1. Spike aqueous medium with surrogate viruses and run blank control.
2. Harvest plant parts and record fresh mass.
3. Run plaque assay and compute infectivity metrics (Eq.1-Eq.5).
4. Run RT-PCR and compute PE/TF/BCF (Eq.6-Eq.8).
5. Estimate HDp from nucleic-acid mass (Eq.9).
6. Subtract blank decay and close the mass balance with `preharvest_mass_balance`.

See `manuscript.md` for experimental details and `tests/test_antiviral_mass_balance.py` for executable examples.

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

Also cite the software release recorded in `CITATION.cff`.

## License

MIT (see `LICENSE`).
