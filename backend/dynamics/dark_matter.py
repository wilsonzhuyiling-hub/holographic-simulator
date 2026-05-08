"""Stage 2 dark-matter diagnostic curve.

This module keeps the dark-matter curve separate from the Stage 1 ODE system.
The Stage 1 solver still computes the background trajectory.  Stage 2 reads
that trajectory, applies the model-endogenous kappa factor, and marks the
relaxation interval discovered from the ODE output.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.interpolate import CubicSpline

from dynamics.critical_points import detect_all
from dynamics.solver import SolverResult

_PI = math.pi
_E = math.e

OMEGA_DM_BASE: float = _PI**4 * (math.sqrt(4.0 * _PI) - _E) / (6.0 * _E**4)
"""Base dark-matter density from P2 Sec. 5."""

F_SPH: float = math.sqrt(4.0 * _PI) / _E
"""First spherical release factor."""

LAMBDA_CAPACITY: float = 1.0 / (F_SPH * _E**2)
"""Geometric capacity factor."""

KAPPA_GEOM: float = 1.0 / math.sqrt(1.0 - LAMBDA_CAPACITY)
"""Model-endogenous geometric compression factor."""

OMEGA_DM_FP: float = KAPPA_GEOM * OMEGA_DM_BASE
"""Kappa-corrected dark-matter fixed-point density."""


@dataclass(frozen=True)
class DarkMatterCurveResult:
    """Continuous Omega_DM diagnostic curve."""

    z_arr: np.ndarray
    omega_dm_arr: np.ndarray
    relaxation_loss_fraction_arr: np.ndarray
    omega_dm_of_z: CubicSpline
    relaxation_loss_fraction_of_z: CubicSpline
    omega_dm_base: float
    kappa_geom: float
    omega_dm_fp: float
    z_peak: float
    z_relaxation_end: float


@dataclass(frozen=True)
class OmegaDMEffDemoResult:
    """DESI-style redshift-shell sampling of Omega_DM_eff(z)."""

    z_eff_arr: np.ndarray
    z_lo_arr: np.ndarray
    z_hi_arr: np.ndarray
    tau_eff_arr: np.ndarray
    C_ratio_eff_arr: np.ndarray
    omega_dm_eff_arr: np.ndarray
    relaxation_loss_fraction_arr: np.ndarray
    continuous_curve: DarkMatterCurveResult


def relaxation_loss_fraction(
    result: SolverResult,
    z: float,
    z_peak: float,
    z_relaxation_end: float,
) -> float:
    """
    Point diagnostic for the C_peak -> C/C1=1 relaxation interval.

    This is kept as a local source kernel for debugging.  The production
    Omega_DM_eff curve uses the causal integral in _relaxation_loss_array(),
    because relaxation is history-dependent rather than point-symmetric.

    The returned value is zero outside the relaxation interval and proportional
    to the current supersaturation Delta C/C1 inside it.
    """
    if z <= z_peak or z >= z_relaxation_end:
        return 0.0

    c_peak_ratio = float(result.C_ratio_of_z(z_peak))
    c_ratio = float(result.C_ratio_of_z(z))
    delta_actual = max(0.0, c_ratio - 1.0)
    return min(delta_actual / (2.0 * c_peak_ratio), 1.0)


def _relaxation_loss_array(
    result: SolverResult,
    z_peak: float,
    z_relaxation_end: float,
) -> np.ndarray:
    """
    Causal relaxation loss over the solver grid.

    The loss state L starts at zero at C_peak and evolves only forward along
    the observer's redshift shells:

        dL/dtau = max(C/C1 - 1, 0) - 2L

    This implements a one-sided relaxation response: the source is set by the
    current supersaturation, while the -2L term erases memory only gradually.
    The final loss fraction is normalized by 2*C_peak to match the dimensionless
    scale used by the legacy relaxation diagnostic.
    """
    z_arr = result.z_arr
    tau_arr = result.tau_arr
    c_ratio_arr = result.C_ratio_arr
    c_peak_ratio = float(result.C_ratio_of_z(z_peak))

    state = 0.0
    raw_loss = np.zeros_like(z_arr, dtype=float)

    start_idx = int(np.searchsorted(z_arr, z_peak, side="left"))
    end_idx = int(np.searchsorted(z_arr, z_relaxation_end, side="right"))

    for idx in range(start_idx + 1, end_idx):
        dtau = max(0.0, float(tau_arr[idx] - tau_arr[idx - 1]))
        source_prev = max(0.0, float(c_ratio_arr[idx - 1] - 1.0))
        source_now = max(0.0, float(c_ratio_arr[idx] - 1.0))
        source_mid = 0.5 * (source_prev + source_now)
        state += dtau * (source_mid - 2.0 * state)
        state = max(0.0, state)
        raw_loss[idx] = state

    return np.clip(raw_loss / (2.0 * c_peak_ratio), 0.0, 1.0)


def compute_omega_dm_eff_curve(result: SolverResult) -> DarkMatterCurveResult:
    """
    Build the lightcone Omega_DM_eff(z) curve with the relaxation interval.

    Boundaries are detected from the solved trajectory:
    - start: C_peak, same point as tau_min
    - end: weak overflow upper boundary, where C/C1 returns to 1
    """
    critical = detect_all(result)
    if critical.overflow_weak is None:
        raise ValueError("Cannot build relaxation curve: weak overflow window was not detected.")

    z_peak = critical.z_peak
    z_relaxation_end = critical.overflow_weak.z_hi
    z_arr = result.z_arr

    loss_arr = _relaxation_loss_array(result, z_peak, z_relaxation_end)
    omega_dm_arr = OMEGA_DM_FP * (1.0 - loss_arr)

    return DarkMatterCurveResult(
        z_arr=z_arr,
        omega_dm_arr=omega_dm_arr,
        relaxation_loss_fraction_arr=loss_arr,
        omega_dm_of_z=CubicSpline(z_arr, omega_dm_arr),
        relaxation_loss_fraction_of_z=CubicSpline(z_arr, loss_arr),
        omega_dm_base=OMEGA_DM_BASE,
        kappa_geom=KAPPA_GEOM,
        omega_dm_fp=OMEGA_DM_FP,
        z_peak=z_peak,
        z_relaxation_end=z_relaxation_end,
    )


def compute_dark_matter_curve(result: SolverResult) -> DarkMatterCurveResult:
    """Backward-compatible alias for compute_omega_dm_eff_curve()."""
    return compute_omega_dm_eff_curve(result)


def compute_omega_dm_eff_demo(
    result: SolverResult,
    z_eff: np.ndarray | None = None,
    shell_width: float = 0.10,
) -> OmegaDMEffDemoResult:
    """
    Sample Omega_DM_eff(z) as DESI-style redshift shells.

    DESI reports observables by effective redshift for finite redshift regions.
    This demo follows that observation principle: the continuous model is not
    refit or re-anchored at those points; it is only evaluated on the observer's
    past lightcone at each shell's z_eff.

    Parameters
    ----------
    result      : solved Stage 1 ODE trajectory
    z_eff       : effective redshifts for the shells.  If omitted, a small
                  DESI-like demo grid is used without importing data values.
    shell_width : displayed shell width in redshift units
    """
    curve = compute_omega_dm_eff_curve(result)
    if z_eff is None:
        z_eff_arr = np.array([0.30, 0.51, 0.71, 0.93, 1.32, 1.49, 2.33])
    else:
        z_eff_arr = np.asarray(z_eff, dtype=float)

    if np.any(z_eff_arr < result.z_arr[0]) or np.any(z_eff_arr > result.z_arr[-1]):
        raise ValueError("All effective redshifts must lie inside the solver range.")
    if shell_width <= 0.0:
        raise ValueError("shell_width must be positive.")

    half_width = 0.5 * shell_width
    z_lo_arr = np.maximum(result.z_arr[0], z_eff_arr - half_width)
    z_hi_arr = np.minimum(result.z_arr[-1], z_eff_arr + half_width)

    return OmegaDMEffDemoResult(
        z_eff_arr=z_eff_arr,
        z_lo_arr=z_lo_arr,
        z_hi_arr=z_hi_arr,
        tau_eff_arr=np.array([float(result.tau_of_z(z)) for z in z_eff_arr]),
        C_ratio_eff_arr=np.array([float(result.C_ratio_of_z(z)) for z in z_eff_arr]),
        omega_dm_eff_arr=np.array([float(curve.omega_dm_of_z(z)) for z in z_eff_arr]),
        relaxation_loss_fraction_arr=np.array([
            float(curve.relaxation_loss_fraction_of_z(z))
            for z in z_eff_arr
        ]),
        continuous_curve=curve,
    )
