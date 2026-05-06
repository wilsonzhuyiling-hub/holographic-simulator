# Architecture

This document describes the module structure, data contracts, and design rationale of the holographic-simulator codebase.

---

## Layered architecture

```
┌──────────────────────────────────────────────────────────┐
│  Layer 4 — Frontend (React + D3 + Tailwind)              │
│  Visual exploration, residual panels, paper provenance   │
└──────────────────────────────────────────────────────────┘
                             ↑ HTTP/JSON
┌──────────────────────────────────────────────────────────┐
│  Layer 3 — API (FastAPI)                                  │
│  Stateless endpoints, Pydantic schemas, OpenAPI docs     │
└──────────────────────────────────────────────────────────┘
                             ↑ Python imports
┌──────────────────────────────────────────────────────────┐
│  Layer 2 — Statistics (chi2, GP, Bayes factor)           │
│  Continuous-domain comparison, no binning loss           │
└──────────────────────────────────────────────────────────┘
                             ↑
┌──────────────────────────────────────────────────────────┐
│  Layer 1 — Dynamics (ODE solver, lightcone projection)   │
│  Self-consistent matter-deformed evolution               │
└──────────────────────────────────────────────────────────┘
                             ↑
┌──────────────────────────────────────────────────────────┐
│  Layer 0 — Core (axioms, derived constants, events)      │
│  Pure-function constant generators (legacy from v1.0)    │
└──────────────────────────────────────────────────────────┘

                             ↑
┌──────────────────────────────────────────────────────────┐
│  Observation — DESI / Planck loaders                      │
│  External data, used only for validation                  │
└──────────────────────────────────────────────────────────┘
```

Each layer depends only on layers below it. There are no circular imports.

---

## Module responsibilities

### `backend/core/` — Layer 0

**Status**: Inherited from v1.0 prototype, working.

Pure-function generators for the framework's derived constants. Computes `C₁`, `Ω_Λ_fp`, `Ω_DM_base`, `f_sph`, `κ_geom`, etc. directly from axioms (π, e, integer dimensions). Also generates particle event tables, force validation tables, and CMB/W-Z bridges.

**Key principle**: Functions in this layer are *deterministic* and *parameter-free*. The only inputs are physical constants used as unit bridges (CODATA `m_e`, `c`, `ℏ`, etc.) and validation targets (PDG, Planck observed values). No fitting.

### `backend/dynamics/` — Layer 1

**Status**: Stage 1 deliverable, not yet implemented.

Solves the matter-deformed self-consistent ODE system from P2 Sec. 3:

```
dτ/dz = -(1 - ε(z)) / (1 + z)
ε(z)  = (3/2) Ω_m(z)
Ω_m(z) = Ω_m0 (1+z)³ / [Ω_m0 (1+z)³ + Ω_Λ(z)]
Ω_Λ(z) = Ω_Λ0 exp[-2(τ(z) - 1)]
```

with boundary `τ(z=0) = 1`, `Ω_Λ0 = π³/(6e²)`, `Ω_m0 = 0.3` (validation input, hardcoded).

Submodules:

- `ode_system.py` — defines the ODE right-hand side
- `solver.py` — `scipy.integrate.solve_ivp` wrapper, RK45 with adaptive stepping
- `critical_points.py` — automatic detection of `τ_min`, `C_peak`, overflow window edges (no human-placed anchors)
- `lightcone.py` — projects physical history (τ-domain) to observer coordinates (z-domain)

**Key principle**: The output is two continuous curves, not discrete event markers. Critical points are *discovered* numerically, not declared.

### `backend/observation/` — supporting layer

**Status**: Stage 2 deliverable.

Loaders for public datasets. Each dataset lives in `data/<source>/` with a small README documenting origin, retrieval date, format, and license. The MVP target is DESI DR2 BAO measurements with full covariance matrix.

### `backend/statistics/` — Layer 2

**Status**: Stage 2 deliverable.

Three improvements over standard ΛCDM analysis:

- `chi2.py` — full-covariance χ²: `χ² = ΔH^T · C⁻¹ · ΔH`. Includes a diagonal-approximation version for direct comparison.
- `gp_reconstruction.py` — Gaussian Process regression on discrete observations to produce a continuous reconstruction comparable point-for-point with the framework's continuous prediction.
- `bayes_factor.py` — model comparison via Bayes factor, rewarding the framework's zero-parameter prior simplicity.
- `binning_loss.py` — quantifies, in bits or KL divergence, how much information is lost by ΛCDM-style binning of a signal that has fine structure within bins.

### `backend/api/` — Layer 3

**Status**: Stage 3 deliverable.

FastAPI endpoints. Route plan:

```
GET  /api/theory       — physical history curve in τ
GET  /api/lightcone    — observer's view in z
GET  /api/observation  — DESI data + covariance
GET  /api/residual     — three-way comparison (model / ΛCDM / data)
GET  /api/critical     — auto-detected critical points
GET  /api/statistics   — χ², Bayes factor, binning loss
GET  /api/explore/{z}  — full diagnostic at a chosen redshift
```

All responses are validated against Pydantic schemas in `api/schemas.py` and reflected in `/docs` via FastAPI's automatic OpenAPI generation.

### `frontend/` — Layer 4

**Status**: Stage 4 deliverable.

React + Vite + D3 + Tailwind. Designed for migration from a static single-page prototype to a Flipbook-style nested-navigation app without backend changes — the API contract is the stable interface.

---

## Data contract principles

Three rules govern data flow between layers:

1. **No layer mutates data from a lower layer.** All transformations produce new objects.
2. **All quantities carry provenance.** Every numerical value should be traceable to its formula and paper reference (see `physics_to_code_map.md`).
3. **Parameter classification is explicit.** Each quantity is tagged as one of:
   - `axiom` — π, e, integer dimensions, fixed-point τ
   - `endogenous` — derived from axioms only
   - `bridge` — CODATA-style unit conversion or experimental anchor
   - `posterior` — diagnostic comparing model to observation

---

## Performance considerations

The MVP does not optimize. After Stage 4, profiling targets:

- ODE solve at default precision: < 100ms (well within FastAPI request budget)
- GP reconstruction on DESI DR2 (~50 points): < 500ms
- Cache theory and lightcone curves at the API layer; recompute only on parameter change (only Ω_m0 is a knob, and it's hardcoded for now)
