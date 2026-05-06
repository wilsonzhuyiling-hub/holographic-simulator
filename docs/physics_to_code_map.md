# Physics-to-code map

This document maps each formula in the Holographic Geometric Framework paper series (P0–P4) to its implementation in the codebase. It exists for two reasons:

1. **Traceability** — when a number changes, you can find which paper section justifies it.
2. **Patching** — when a paper is updated, you can find which code needs revision.

The table is the source of truth. Entries marked **(stage X)** belong to a future deliverable.

---

## Layer 0 — Core constants

| Quantity | Formula | Paper ref | Code location |
|---|---|---|---|
| `C₀` | `π² / 16` | P1 Eq. (8) | `core/universe_simulator.py` → `derive_C0_axis()` |
| `C₁` | `π² / (16 e²)` | P1 main result, Eq. (12) | `derive_C1()` |
| `Ω_Λ_fp` | `π³ / (6 e²)` | P2 Eq. (9) | `derive_omega_lambda_fp()` |
| `Ω_DM_base` | `π⁴(√(4π) − e) / (6 e⁴)` | P2 Sec. 5 | `derive_omega_dm()` |
| `f_sph` | `√(4π) / e` | P2 Eq. (15), P3 Sec. 4 | `derive_f_sph()` |
| `C_sphere` | `C₁ · f_sph = π² √(4π) / (16 e³)` | P2 Eq. (16) | `DerivedConstants.C_sphere` |
| `λ_capacity` | `1 / (f_sph · e²)` | P2 Sec. 5.2 | `derive_lambda_capacity()` |
| `κ_geom` | `1 / √(1 − λ)` | P2 Sec. 5.2, P5.1 | `derive_kappa_geom()` |
| `Ω_DM_fp` | `κ · Ω_DM_base` | P2 Eq. (17) | `derive_omega_dm_fp()` |

---

## τ ↔ z mapping

| Quantity | Formula | Paper ref | Code location |
|---|---|---|---|
| `τ(z)` (de Sitter limit) | `1 − ln(1 + z)` | P2 Eq. (5) | `tau_from_redshift()` |
| `z(τ)` | `e^(1−τ) − 1` | P2 inverse of (5) | `redshift_from_tau()` |
| Scale factor | `a(τ) = e^(τ−1)` | P2 Eq. (4) | `scale_factor_from_tau()` |
| τ-temperature | `T(τ) = T₀ · e^(−(τ−τ₀))` | P2 Sec. 6 | `temperature_from_tau()` |
| **CMB bridge** | `T_CMB = T_rec · e^(−(1−τ_rec))` | P2 Sec. 6 | `T_CMB_from_recombination()` |

---

## Geometric coefficient

| Quantity | Formula | Paper ref | Code location |
|---|---|---|---|
| `C(R, τ_R, τ_I)` | `(π²/16)(l_UV/R)² e^(−2τ_R) e^(−2iτ_I)` | P4 Sec. 3 | `C_complex()` |
| `C̄(R, τ_R, τ_I)` | conjugate-mirror partner | P4 Sec. 3.2 | `C_conjugate_partner()` |
| `Z` (coherence) | `\|C + C̄\|²` | P4 Sec. 4 | `coherence_Z()` |
| Conservation product | `C · C̄ = const` | P4 Sec. 3.3 | `conservation_product()` |

---

## Dynamics — Stage 1 deliverables

These are not yet implemented. The expected mapping:

| Quantity | Formula | Paper ref | Planned location |
|---|---|---|---|
| ODE: `dτ/dz` | `−(1 − ε(z)) / (1 + z)` | P2 Eq. (10) | `dynamics/ode_system.py::dtau_dz` |
| ε(z) | `(3/2) Ω_m(z)` | P2 Eq. (11) | `dynamics/ode_system.py::epsilon` |
| Ω_m(z) | matter fraction closure | P2 Eq. (12) | `dynamics/ode_system.py::omega_m_of_z` |
| Ω_Λ(z) | `Ω_Λ0 · e^(−2(τ(z)−1))` | P2 Eq. (13) | `dynamics/ode_system.py::omega_lambda_of_z` |
| τ_min detection | numerical, ≈ 0.808 at z ≈ 0.898 | P2 Sec. 3.3 | `dynamics/critical_points.py::find_tau_minimum` |
| C_peak detection | numerical, C/C₁ ≈ 1.467 | P2 Sec. 3.3 | `dynamics/critical_points.py::find_C_peak` |
| Overflow window | `[z₁, z₂]` where `C > C_sphere` | P2 Sec. 4.1 | `dynamics/critical_points.py::overflow_window` |

---

## Bridges and validation targets

These quantities are *not* model predictions. They are CODATA values, experimental anchors, or PDG/Planck targets used only for unit conversion or validation.

| Quantity | Source | Code location |
|---|---|---|
| `c, h, ℏ, k_B, α` | CODATA | `PhysicalInputs` |
| `m_e_kg` | CODATA electron mass | `PhysicalInputs.m_e_kg` |
| `m_p_MeV, m_n_MeV` | PDG | `PhysicalInputs` |
| `T_rec_K, z_rec` | recombination physics | `PhysicalInputs` (until internally closed) |
| `T_CMB_observed_K` | Fixsen 2009 | validation only |
| `m_W, m_Z observed` | PDG | validation only |
| `Ω_Λ_PDG_2024` | PDG 2024 | validation only |
| `Ω_DM_Planck_derived` | Planck 2018 | validation only |
| **`Ω_m0 = 0.3`** | observational input | hardcoded in `dynamics/` (Stage 1) |

---

## Quadrant assignments (P4)

| Quadrant | Sign(τ_R, τ_I) | Physics | Code |
|---|---|---|---|
| QI | (+, +) | stable matter / atomic binding | `quadrant_of()` returns "QI stable matter" |
| QII | (−, +) | Z⁰ neutral weak saddle | "QII Z0 saddle" |
| QIII | (−, −) | W⁻ charged weak saddle | "QIII W- saddle" |
| QIV | (+, −) | cosmological radiation / CMB | "QIV radiation background" |

---

## Open derivations (not yet closed)

These are explicitly marked in the framework as boundaries, not failures. They become bridges/posteriors until closed:

| Quantity | Status | Plan |
|---|---|---|
| `α` itself | P4 records relation `α = (4/3)(m_e/m_s)`; not derived from geometry alone | Future paper |
| `m_u` (up-quark) | not closed; no free parameter introduced | Awaiting QI island integral closure |
| Recombination temperature `T_rec` | external until internal closure found | — |
| Z⁰ absolute mass anchor | external | — |

---

## How to update this map

When you add or change a quantity:

1. Add a row in the appropriate table
2. Cite the paper section/equation
3. Link to the code symbol (function or attribute path)
4. If the quantity is bridge or posterior, say so explicitly

If a paper revision changes a formula, the corresponding code change must update this map in the same PR.
