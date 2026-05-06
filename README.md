# Holographic Simulator

> ⚠️ **Status: Alpha (work in progress)** — Stage 0 scaffolding is complete. The dynamics engine is under active development on the Stage 1 branch, and the statistics/API/frontend layers are not yet ready for use.

Open-source simulator for the **Holographic Geometric Framework**, a zero-free-parameter theoretical framework for cosmological observables. The simulator solves the framework's self-consistent ODE system and compares its predictions to public observational data (DESI, Planck) using continuous-domain statistical methods that avoid the binning, diagonal-covariance, and free-parameter limitations of standard ΛCDM analyses.

---

## Why this exists

Standard cosmological data analysis has three structural limitations:

1. **Discrete binning loss** — observations are aggregated into wide redshift bins (Δz ≈ 0.1–0.2), erasing fine structure that a continuous theory could otherwise predict.
2. **Diagonal-covariance χ²** — most analyses assume bin-to-bin independence, systematically underestimating significance when the true signal has correlated structure.
3. **Free-parameter absorption** — multi-parameter models (w₀wₐCDM and similar) absorb genuine theoretical signatures into nuisance parameters.

The Holographic Geometric Framework predicts a specific, **rigid** functional form for `H(z)`, `Ω_Λ(z)`, and `fσ₈(z)` from the geometric coefficient `C(τ) = π²/16 · e^(−2τ)` and the self-consistent matter-deformed ODE system. Because the framework has zero free parameters, its predictions can be tested directly against data without fitting.

This simulator does two things:

- **Generates the theoretical history curve** — what the universe actually does, in continuous τ.
- **Generates the lightcone observation curve** — what an observer at z = 0 sees today, in continuous z.

The two curves can then be compared to each other, to ΛCDM, and to real observations using full-covariance χ², Gaussian Process reconstruction, and Bayes factor analysis.

---

## Architecture

```
┌─────────────────────────────────────────────────┐
│              Frontend (React + D3)              │
│   Dual-curve view · Residual panel · Statistics │
└─────────────────────────────────────────────────┘
                       ↑ JSON over HTTP
┌─────────────────────────────────────────────────┐
│              API Layer (FastAPI)                │
└─────────────────────────────────────────────────┘
                       ↑ Python imports
┌─────────────────────────────────────────────────┐
│   Computational Core (Python · scipy · sklearn) │
│   • core/         — derived constants & events  │
│   • dynamics/     — ODE solver, critical points │
│   • observation/  — DESI/Planck data loaders    │
│   • statistics/   — χ², GP, Bayes factor        │
└─────────────────────────────────────────────────┘
```

See [`docs/architecture.md`](docs/architecture.md) for module-level detail.

---

## Roadmap

| Stage | Status | Deliverable |
|---|---|---|
| 0 — Scaffolding | ✅ complete | Repo structure, CI, docs skeleton, editable backend install |
| 1 — Dynamics engine | 🟡 in progress | RK45 ODE solver, automatic critical-point detection (recovers τ_min ≈ 0.808, C_peak/C₁ ≈ 1.467) |
| 2 — Statistics layer | 🔲 planned | DESI DR2 loader, full-covariance χ², GP reconstruction, Bayes factor |
| 3 — API service | 🔲 planned | FastAPI endpoints, Pydantic schemas, OpenAPI docs |
| 4 — Frontend prototype | 🔲 planned | Dual-curve view, residual panel, paper provenance popover |
| 5 — Polish & deploy | 🔲 planned | Docker, Hugging Face Spaces deployment |

---

## Theoretical foundation

This simulator implements the framework described in the **Holographic Geometric Framework paper series** (Zhu 2026):

- **P0** — Reader's introduction
- **P1** — Geometric fixed-point reading of the CKN coefficient: `C₁ = π²/(16e²)`
- **P2** — From the de Sitter fixed point to the dark universe: relaxation dynamics and geometric overflow
- **P3** — Holographic constraint at laboratory and nuclear scales
- **P4** — Real to complex: holographic double-helix and Standard Model structure

Code-to-paper provenance is tracked in [`docs/physics_to_code_map.md`](docs/physics_to_code_map.md).

---

## Quick start (developers)

> The API and frontend do not exist yet. Stage 0 setup is complete; the legacy event/constant generator remains available:

```bash
cd backend
python core/universe_simulator.py
```

This produces `outputs/*.csv` and `outputs/*.json` containing derived constants, particle events, force validations, and the dark-matter relaxation table.

Developer setup:

```bash
cd backend
pip install -e ".[dev]"
pytest
```

---

## Contributing

Issues and discussion welcome once Stage 1 is complete. See [`CONTRIBUTING.md`](CONTRIBUTING.md) for guidelines.

Three types of contributions are explicitly invited:

- **Physics review** — questions or critiques of the theoretical framework (cite the relevant paper section)
- **Engineering** — code, tests, infrastructure
- **Data** — additional observational datasets or improved loaders

Use the corresponding Issue template.

---

## License

[MIT](LICENSE) — free for academic and commercial use, with attribution.

---

## Author

Yiling (Wilson) Zhu — Independent researcher, London
[wilsonzhuyiling@gmail.com](mailto:wilsonzhuyiling@gmail.com)
