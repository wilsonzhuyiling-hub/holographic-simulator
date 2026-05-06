# Contributing to Holographic Simulator

Thanks for your interest. This project is in early alpha — the dynamics engine and statistics layer are not yet implemented. Contributions are welcome but please read this guide first.

---

## Branching model

- `main` — stable, only updated when a stage is complete and reviewed
- `dev` — integration branch for in-progress work
- `feature/<name>` — individual feature branches, branched off `dev`

Pull requests should target `dev`, not `main`.

---

## Three types of contributions

### Physics

Questions, critiques, or extensions of the theoretical framework. Use the **Physics** issue template. Reference the relevant paper (P0–P4) and section/equation number.

If you're proposing a new physical relation or correcting an existing one, please attach:

- The mathematical derivation (LaTeX in the issue body, or a linked PDF)
- The expected impact on `physics_to_code_map.md`
- Whether it changes any model-endogenous quantity vs. a bridge/posterior quantity

### Engineering

Code, tests, CI, infrastructure, performance improvements. Use the **Engineering** issue template.

Conventions:

- Python: type hints required, `ruff` for lint, `pytest` for tests
- New numerical code must include a unit test that locks the numerical value (so future refactors don't silently change physics)
- Public API surface (anything in `backend/api/`) requires a Pydantic schema

### Data

Additional observational datasets, improved loaders, covariance matrices. Use the **Data** issue template.

Data files go in `backend/observation/data/<source>/` with a small `README.md` describing:

- Origin (URL, DOI, paper citation)
- Date retrieved
- Format and units
- Whether covariance is included

---

## Development setup

```bash
git clone https://github.com/wilsonzhuyiling-hub/holographic-simulator.git
cd holographic-simulator
git checkout dev
git checkout -b feature/your-feature-name

# Backend
cd backend
pip install -e ".[dev]"
pytest
```

Frontend setup will be documented once Stage 4 begins.

---

## Pull request checklist

- [ ] Branched off `dev`, targeting `dev`
- [ ] Tests added or updated
- [ ] `pytest` passes locally
- [ ] If touching physics: `physics_to_code_map.md` updated
- [ ] If touching API: schema updated, OpenAPI docs regenerate cleanly
- [ ] Commit messages are descriptive (not "fix" or "update")

---

## Code of conduct

Be respectful. Disagreement on physics is welcome and expected; ad hominem is not. We're trying to build something honest about the limits of current cosmological data analysis — that requires charitable engagement with critique.
