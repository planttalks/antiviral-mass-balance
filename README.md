# Antiviral mass-balance (Chemosphere companion code)

**Repository:** [github.com/planttalks/antiviral-mass-balance](https://github.com/planttalks/antiviral-mass-balance) · **Maintainer:** [@planttalks](https://github.com/planttalks)

Pure-Python helpers for **Eq.1–Eq.9** and a **pre-harvest culture closure** described in the Materials and methods of:

**Diaiti Zure**, **Hsion-Wen David Kuo**, **Aleksandra Drizo** (2024). Insights of phytoremediation mechanisms for viruses based on in-vitro, in-vivo and in-silico assessments of selected herbal plants. *Chemosphere*, 351, 141101.

- DOI: [https://doi.org/10.1016/j.chemosphere.2023.141101](https://doi.org/10.1016/j.chemosphere.2023.141101)
- ScienceDirect (PII `S0045653523033714`): [article link](https://www.sciencedirect.com/science/article/abs/pii/S0045653523033714)

Author names and order match the [Crossref record](https://api.crossref.org/works/10.1016/j.chemosphere.2023.141101) for that DOI (the PII above resolves to the same article).

Use **one consistent unit system** (PFU, calibrated gene copies, or Eq.9-derived particle estimates) before interpreting `closure_fraction` and `residual_fraction`.

## Pre-harvest workflow: aqueous spike → assays → mass balance

The steps below follow the *Materials and methods* in the paper. Protocol detail, instrument settings, Text S1 and figures sit in `manuscript.md` in this folder (use that file as the full lab narrative).

### 1. Culture and spike (in-vivo)

Grow plants under your chosen pre-harvest system (the paper used greenhouse Kratky-style non-circulating hydroponics). Move intact plants into an aqueous volume and spike with surrogate bacteriophages (e.g. MS2 and T4 stock titers as in the manuscript). Run controls in parallel: virus-free water and herb-free spiked water. The herb-free spiked control is the **blank** used later for **natural decay (IAe)**.

### 2. Sampling time points and fresh mass

Sample bulk water, defined plant parts and controls at fixed times after exposure (e.g. 0, 2 and 4 h). Record **fresh weight (fwt)** for each harvested part. You need mass to convert extract-based plaque readouts to a per-gram basis and to weight multi-part combinations.

### 3. Tissue preparation

Prepare plant material for molecular and infectivity work using your validated prep (the paper cites Li & Uyttendaele (2018) with study-specific handling). The goal is reproducible lysates or extracts whose volume you will pair with plaque counts and PCR targets.

### 4. Double-layer agar plaque assay (infectivity)

Use the plaque method referenced in the paper (Cock & Kalt (2010) style) to obtain **plaque counts** at known **dilution** and **virus dilution volume**. Feed those into **Eq.2** (`initial_infectivity_pfu_ml`) to get **A** in PFU/mL for water, extracts, or part-specific supernatants.

For a given plant part supernatant, combine **A**, **extract volume** and **fwt** through **Eq.3** (`final_infectivity_pfu_per_g`) to get **PI** (PFU/g). **Eq.1** (`viral_inactivation_percent`) compares initial versus final infectivity for **VI (%)** when you substitute the correct initial and final pair.

**In-vitro** screens on crude extracts use the same plaque pipeline. For multi-part composites, the paper’s **Eq.4** has a typo in print; this repo implements the mass-weighted form via `final_combined_infectivity_weighted`. **Eq.5** maps mixture plaque results to **MI** (`final_mixture_infectivity_pfu_per_g`) when you track **sum of part masses** and **whole-herb mass** as defined in the methods.

### 5. Real-time PCR (particle counts in water and tissues)

Quantify viral genomes or targets in water and in each part using your RT-PCR workflow (the paper used StepOne; see Text S1 in `manuscript.md`). Use those concentrations with **Eq.6–Eq.8** for **PE (%)** (`phytoremediation_efficiency_percent`), **TF** (`translocation_factor`) and **BCF** (`bioconcentration_factor`) once you have compatible “particle” units in water, root and target part.

Treat **infective (Ap)** and **inactivated (IAp)** tiers in the plant with the plaque-centered logic described in the paper: **Ap** ties to **Eq.3** style readouts; **IAp** follows the same measurement family with **Eq.1** applied to distinguish loss of infectivity versus the starting plaque-based level. Calibrate RT-PCR to PFU if your study applies a PFU anchor (the manuscript discusses MS2 PFU versus RT-PCR bias).

### 6. Highly degraded fraction (HDp) from nucleic acid mass

For material that is no longer captured well by plaque or PCR, measure released nucleic acid (e.g. NanoDrop ng/µL) and convert with **Eq.9** (`highly_degraded_vp_per_ml`) using the **genome length** and **manuscript molecular-weight terms** bundled as constants for MS2 and T4 in `antiviral_mass_balance.py`.

### 7. Blank correction (IAe) before interpreting plant pools

From the herb-free spiked control, estimate **IAe** the same way you estimate inactivation-related loss in bulk liquid. **Subtract that decay from the effective initial load** before you interpret **IAp** and **HDp** inside plant tissues (`natural_decay_subtracted_initial`).

### 8. Closing the pre-harvest mass balance in code

Convert every compartment to **total** counts in one unit system (multiply per-mL or per-g values by **culture volume** or **part mass** as appropriate). Build one `PlantPartLoads` per harvested part with totals for **Ap**, **IAp** and **HDp**. Call `preharvest_mass_balance` with (1) culture volume, (2) initial concentration in the liquid, (3) **IAe** percent from the blank, (4) water concentration at the sampled time and (5) `parts` as a sequence of `PlantPartLoads`. Inspect `closure_fraction` and `n_residual` against known assay limits (the paper discusses NanoDrop bias, host DNA and PCR recovery).

For more information on rationale, statistics and extended methods, read **`manuscript.md`**.

## Install (local / editable)

```bash
pip install -e ".[dev]"
```

## Run tests and linters

```bash
ruff check .
ruff format --check .
pytest --cov=antiviral_mass_balance --cov-fail-under=85
```

## CI

[![CI](https://github.com/planttalks/antiviral-mass-balance/actions/workflows/ci.yml/badge.svg)](https://github.com/planttalks/antiviral-mass-balance/actions/workflows/ci.yml)

Workflow `.github/workflows/ci.yml` runs Ruff and pytest on Python 3.10–3.12 on each push and pull request to `main`.

## Citation

See `CITATION.cff` for software metadata. Cite the *Chemosphere* article for the experimental and modeling context.

## License

MIT — see `LICENSE`.
