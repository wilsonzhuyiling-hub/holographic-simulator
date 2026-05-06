"""BAO distance observables for the holographic framework.

Converts the dimensionless ODE output (E(z), τ(z)) into the physical distance
quantities that DESI and other BAO surveys measure.

The conversion requires Tier 3 foundational data (Planck 2018):
    H₀ = 67.36 km/s/Mpc   — sets the Hubble distance D_H0 = c/H₀
    r_d = 147.09 Mpc       — baryon drag sound horizon (BAO standard ruler)

These values do not enter the ODE system and have no effect on the intrinsic
model curves (τ, E, Ω_m, Ω_Λ, C/C₁).  They serve only as unit-conversion
anchors that map the dimensionless prediction onto the observable axis.

Quantities computed
-------------------
D_C(z)  comoving distance     = D_H0 · ∫₀ᶻ dz'/E(z')        [Mpc]
D_M(z)  transverse comoving   = D_C(z)   (flat universe)     [Mpc]
D_H(z)  Hubble distance       = c/H(z) = D_H0/E(z)           [Mpc]
D_V(z)  spherically-averaged  = [z · D_M² · D_H]^(1/3)       [Mpc]

BAO ratio observables (dimensionless, direct DESI comparison):
    D_M(z) / r_d
    D_H(z) / r_d
    D_V(z) / r_d
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.interpolate import CubicSpline

from dynamics.solver import SolverResult


# ---------------------------------------------------------------------------
# Tier 3 — Foundational model data (Planck 2018)
# Unit-conversion anchors only; do not enter the ODE system.
# ---------------------------------------------------------------------------

C_KMS: float = 299_792.458
"""Speed of light [km/s] (exact, CODATA 2018)."""

H0_PLANCK: float = 67.36
"""Hubble constant today [km/s/Mpc] — Planck 2018 TT,TE,EE+lowE+lensing."""

R_D_PLANCK: float = 147.09
"""Sound horizon at baryon drag epoch [Mpc] — Planck 2018 TT,TE,EE+lowE+lensing."""

R_D_DESI: float = 147.05
"""
Sound horizon at baryon drag epoch [Mpc] — DESI DR2 fiducial value (arXiv:2503.14738).
Use this when computing BAO ratios for direct DESI DR2 comparison.
P2 Sec. 9 explicitly uses this value: c/(H₀ r_d) = 30.27 with r_d = 147.05 Mpc.
Difference from R_D_PLANCK: 0.04 Mpc (0.027%).
"""

DH0: float = C_KMS / H0_PLANCK
"""Hubble distance D_H0 = c/H₀ ≈ 4451 Mpc."""


# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------

@dataclass
class ObservablesResult:
    """
    BAO distance observables computed from a SolverResult.

    All distance arrays are in Mpc.  BAO ratio observables are dimensionless.
    """
    z_arr: np.ndarray

    D_C_arr: np.ndarray        # comoving distance [Mpc]
    D_M_arr: np.ndarray        # transverse comoving distance [Mpc] (= D_C, flat)
    D_H_arr: np.ndarray        # Hubble distance c/H(z) [Mpc]
    D_V_arr: np.ndarray        # spherically-averaged BAO distance [Mpc]

    # Continuous interpolants
    D_C_of_z: CubicSpline = field(repr=False)
    D_M_of_z: CubicSpline = field(repr=False)
    D_H_of_z: CubicSpline = field(repr=False)
    D_V_of_z: CubicSpline = field(repr=False)

    # BAO ratio observables (dimensionless, direct DESI comparison)
    DM_over_rd_of_z: CubicSpline = field(repr=False)
    DH_over_rd_of_z: CubicSpline = field(repr=False)
    DV_over_rd_of_z: CubicSpline = field(repr=False)

    # Tier 3 anchors recorded for provenance
    H0: float = H0_PLANCK
    r_d: float = R_D_PLANCK


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def compute_observables(
    result: SolverResult,
    H0: float = H0_PLANCK,
    r_d: float = R_D_PLANCK,
) -> ObservablesResult:
    """
    Compute BAO distance observables from a SolverResult.

    Parameters
    ----------
    result : SolverResult from dynamics.solver.solve_history()
    H0     : Hubble constant [km/s/Mpc] — Tier 3 anchor, default Planck 2018
    r_d    : Sound horizon at drag epoch [Mpc] — Tier 3 anchor, default Planck 2018

    Returns
    -------
    ObservablesResult with distance arrays, interpolants, and BAO ratios.
    """
    dh0 = C_KMS / H0
    z_arr = result.z_arr

    # E(z) = H(z)/H₀ on the solver grid
    E_arr = np.array([result.H_normalized(float(z)) for z in z_arr])

    # D_H(z) = c/H(z) = D_H0 / E(z)  [Mpc]
    D_H_arr = dh0 / E_arr

    # D_C(z) = D_H0 · ∫₀ᶻ dz'/E(z')  [Mpc]
    # Build a CubicSpline of 1/E(z) and use its antiderivative for accuracy.
    inv_E_spline = CubicSpline(z_arr, 1.0 / E_arr)
    _antideriv = inv_E_spline.antiderivative()
    D_C_arr = dh0 * (_antideriv(z_arr) - float(_antideriv(0.0)))

    # D_M(z) = D_C(z)  (flat universe)
    D_M_arr = D_C_arr.copy()

    # D_V(z) = [z · D_M(z)² · D_H(z)]^(1/3)  [Mpc]
    # Undefined at z=0; set to 0.0 there.
    D_V_arr = np.where(
        z_arr > 0.0,
        (z_arr * D_M_arr ** 2 * D_H_arr) ** (1.0 / 3.0),
        0.0,
    )

    # Continuous interpolants
    D_C_of_z = CubicSpline(z_arr, D_C_arr)
    D_M_of_z = CubicSpline(z_arr, D_M_arr)
    D_H_of_z = CubicSpline(z_arr, D_H_arr)
    D_V_of_z = CubicSpline(z_arr, D_V_arr)

    # BAO ratio observables
    DM_over_rd_of_z = CubicSpline(z_arr, D_M_arr / r_d)
    DH_over_rd_of_z = CubicSpline(z_arr, D_H_arr / r_d)
    DV_over_rd_of_z = CubicSpline(z_arr, D_V_arr / r_d)

    return ObservablesResult(
        z_arr=z_arr,
        D_C_arr=D_C_arr,
        D_M_arr=D_M_arr,
        D_H_arr=D_H_arr,
        D_V_arr=D_V_arr,
        D_C_of_z=D_C_of_z,
        D_M_of_z=D_M_of_z,
        D_H_of_z=D_H_of_z,
        D_V_of_z=D_V_of_z,
        DM_over_rd_of_z=DM_over_rd_of_z,
        DH_over_rd_of_z=DH_over_rd_of_z,
        DV_over_rd_of_z=DV_over_rd_of_z,
        H0=H0,
        r_d=r_d,
    )
