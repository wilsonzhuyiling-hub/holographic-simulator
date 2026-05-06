# Stage 0 — Scaffolding

**Status**: Complete

This stage sets up the repository structure, CI, documentation skeletons, and migrates the v1.0 prototype code into the new layered architecture. No new physics or algorithms are introduced.

## What's in place

- Repository skeleton: `backend/{core,dynamics,observation,statistics,api,tests}`, `frontend/`, `docs/`, `.github/`
- License (MIT), README, CONTRIBUTING, .gitignore
- `backend/pyproject.toml` with runtime/dev dependencies and editable-install support
- `backend/core/universe_simulator.py` — legacy v1.0 derived-constants generator (unchanged)
- `backend/core/render_universe.py` — legacy renderer (unchanged, will be deprecated by frontend)
- `backend/tests/test_core_constants.py` — locks numerical values of all framework constants
- `backend/tests/test_dynamics_placeholders.py` — Stage 1 validation tests for the dynamics branch
- GitHub Actions CI: `pytest` + `ruff` on push/PR to `main` and `dev`
- Issue templates: physics, engineering, data
- Documentation: `architecture.md`, `physics_to_code_map.md`

## Verifying Stage 0

Local check:

```bash
cd backend
pip install -e ".[dev]"
pytest                  # core constants and current dynamics tests pass
ruff check .            # no lint errors
```

CI check: push triggers `tests` workflow on GitHub Actions. Green badge expected.

Current local verification on the Stage 1 branch:

```text
28 passed
```

## Stage boundary

- Stage 0 is complete: repository structure, packaging, CI, docs skeleton, and legacy core migration are in place.
- The current `feature/dynamics-engine` branch already contains the Stage 1 ODE solver and critical-point detection work.
- Stage 0 introduced no new physics or algorithms; Stage 1 is where the dynamics equations become executable.

## What's NOT in place yet

- DESI loader (Stage 2)
- Statistics layer (Stage 2)
- API service (Stage 3)
- Frontend (Stage 4)
- Docker / deployment (Stage 5)

## Next: Stage 1

Stage 1 lives on the `feature/dynamics-engine` branch. Its implementation targets are:

- `dynamics/ode_system.py`
- `dynamics/solver.py`
- `dynamics/critical_points.py`

Validation: the dynamics tests in `test_dynamics_placeholders.py` should recover the paper target values without hard-coded critical points.
