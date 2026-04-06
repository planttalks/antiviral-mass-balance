# Changelog

All notable changes to this project are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [0.1.0] – 2026-04-06

### Added
- `antiviral_mass_balance.py`: implementations of Eq.1–Eq.9 from Zure et al. (2024)
- `PlantPartLoads` and `PreharvestMassBalance` dataclasses for structured mass-balance closure
- `preharvest_mass_balance()` helper combining culture volume, blank decay and per-part loads
- `dilution_after_spike()` helper for initial concentration after spiking
- Constants `MS2_GENOME_BP`, `MS2_MOLECULAR_WEIGHT_EQ9`, `T4_GENOME_BP`, `T4_MOLECULAR_WEIGHT_EQ9`
- Unit test suite (`tests/test_antiviral_mass_balance.py`, 12 tests, 86 % branch coverage)
- CI workflow (GitHub Actions, Python 3.10–3.13, Ruff lint + format, pytest-cov)
- `CITATION.cff` with preferred citation for the *Chemosphere* article
- `examples/ms2_t4_worked_example.py` reproducing key manuscript results
