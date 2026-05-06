# Stage 0 — Scaffolding

**Status**: Complete

This stage sets up the repository structure, CI, documentation skeletons, and migrates the v1.0 prototype code into the new layered architecture. No new physics or algorithms are introduced.

## What's in place

- Repository skeleton: `backend/{core,dynamics,observation,statistics,api,tests}`, `frontend/`, `docs/`, `.github/`
- License (MIT), README, CONTRIBUTING, .gitignore
- `backend/pyproject.toml` with dev dependencies (pytest, ruff, mypy)
- `backend/core/universe_simulator.py` — legacy v1.0 derived-constants generator (unchanged)
- `backend/core/render_universe.py` — legacy renderer (unchanged, will be deprecated by frontend)
- `backend/tests/test_core_constants.py` — locks numerical values of all framework constants
- `backend/tests/test_dynamics_placeholders.py` — skipped tests for Stage 1 milestones
- GitHub Actions CI: `pytest` + `ruff` on push/PR to `main` and `dev`
- Issue templates: physics, engineering, data
- Documentation: `architecture.md`, `physics_to_code_map.md`

## Verifying Stage 0

Local check:

```bash
cd backend
pip install -e ".[dev]"
pytest                  # all core tests pass; dynamics tests skip
ruff check .            # no lint errors
```

CI check: push triggers `tests` workflow on GitHub Actions. Green badge expected.

## What's NOT in place yet

- ODE solver (Stage 1)
- DESI loader (Stage 2)
- Statistics layer (Stage 2)
- API service (Stage 3)
- Frontend (Stage 4)
- Docker / deployment (Stage 5)

## Next: Stage 1

Open the `feature/dynamics-engine` branch off `dev`. Implement:

- `dynamics/ode_system.py`
- `dynamics/solver.py`
- `dynamics/critical_points.py`

Validation: the four currently-skipped tests in `test_dynamics_placeholders.py` should pass without modification of their expected values.
