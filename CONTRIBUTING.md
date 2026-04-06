# Contributing

Thanks for your interest. This is a small companion package for a published study, so contributions are focused on correctness, reproducibility and usability for the environmental-engineering research community.

## Scope

Accepted contributions:
- Bug fixes in equation implementations (Eq.1–Eq.9)
- Additional unit tests or worked examples
- Documentation improvements
- CI or packaging improvements

Out of scope:
- New scientific models unrelated to Zure et al. (2024)
- Breaking changes to public API without discussion

## Setup

```bash
git clone https://github.com/planttalks/antiviral-mass-balance.git
cd antiviral-mass-balance
pip install -e ".[dev]"
```

## Before submitting a pull request

```bash
ruff check .
ruff format .
mypy antiviral_mass_balance.py
pytest --cov=antiviral_mass_balance --cov-fail-under=85
```

All four commands must pass with no errors.

## Reporting equation errors

If you believe an equation implementation differs from the published paper, open an issue and cite the specific equation number and manuscript line. Include your expected and actual values with units.

## Citation

If this code contributes to published work, cite both the software (see `CITATION.cff`) and the original article:

Zure D, Kuo H-W D, Drizo A. Insights of phytoremediation mechanisms for viruses based on in-vitro, in-vivo and in-silico assessments of selected herbal plants. *Chemosphere*. 2024;351:141101. doi:[10.1016/j.chemosphere.2023.141101](https://doi.org/10.1016/j.chemosphere.2023.141101).
