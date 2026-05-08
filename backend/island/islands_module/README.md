# Black-Hole Island Emergence Module

This module extends the Holographic Geometric Framework (HGF) to model
black-hole evaporation and the emergence of entanglement islands.

## Physical basis

Two inputs are derived directly from the HGF framework (no new free parameters):

**1. Log-horizon time map**
```
τ(t) = ln(1 − t / t_evap)
```
Derived from the native area law A ∝ e^{2τ}. Replaces the earlier linear
parametrisation.

**2. Double-helix macro/micro mirror map**
```
τ_micro = −1 + (1 − τ_macro) / e
```
From Paper III / IV. Acts as the coherence source for the matter term in S_gen.

**3. Area-law Hawking branch**
```
S_Hawking ∝ 1 − (1 − t/t_evap)²
```
Follows from the same area law; produces the correct upward-curving shape.

## Key results (all phase-response variants)

| Quantity | Value |
|---|---|
| t_Page | 0.096 (stable across all 5 variants) |
| τ_macro at Page | 0.504 |
| τ_micro at Page | −0.818 |
| island depth at Page | 0.210 |
| island depth (late) | ≈ 0.97 (interior) |
| coherence evolution | 0.000 → 1.000 |

Page time is an output of the fixed-point structure, not a fitted parameter.

## Files

- `entropy_model.py` — core simulation: S_gen extremisation, island scan, CSV/JSON output
- `render_island.py` — PIL rendering: Page-curve competition + S_gen landscape heatmap
- `README.md` — this file

## Usage

```bash
cd backend
python -m island.islands_module.entropy_model      # runs scan, writes CSV + JSON to outputs/islands/
python -m island.islands_module.render_island      # renders PNG to outputs/islands/
```

## Connection to Wang et al. (2021)

Wang, Li, Wang — *Page curves for a family of exactly solvable evaporating
black holes*, Phys. Rev. D 103, 126026 (2021) — derive analytically that the
Page transition occurs at t_Page = t_evap / 3 for the RST-BPP family.

The structural correspondences with this module are:

1. Page time is geometry-determined, not fitted (both frameworks)
2. Island boundary stabilises inside the horizon post-transition (both frameworks)
3. Branch switch is discrete and irreversible (both frameworks)

The quantitative difference (t_Page/t_evap ≈ 0.10 here vs 0.33 in Wang et al.)
is under investigation. The HGF identity G·ρ_Λ^fp / (c²H₀²) = C₁ fixes the
gravitational coupling side of a potential dictionary; the correspondence for
the CFT central charge N is an open question.
