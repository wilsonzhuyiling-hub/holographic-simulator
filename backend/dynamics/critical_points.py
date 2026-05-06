"""Automatic critical-point detection for the holographic dynamics.

All critical points are *discovered* numerically from the ODE solution.
None are hard-coded from the paper — the paper values serve only as
validation targets in the test suite.

Critical points detected
------------------------
tau_min     : minimum of τ(z), coincides with acceleration crossing (q=0)
C_peak      : maximum of C(τ)/C₁, same z as tau_min
z_accel     : acceleration transition z where q(z) = 0  (ε = 1, Ω_m = 2/3)
overflow_weak  : [z_lo, z_hi] where C/C₁ > 1  (weak supersaturation window)
overflow_strong: [z_lo, z_hi] where C/C₁ > C_sphere/C₁  (strong window)

Paper reference values (P2 Sec. 3.3, used only in tests)
---------------------------------------------------------
τ_min ≈ 0.808 at z ≈ 0.898
C_peak / C₁ ≈ 1.467
Strong overflow window: z ∈ [0.345, 1.721]
Weak overflow window: z ∈ [0.001, 2.846]
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

from dynamics.ode_system import C1, C_SPHERE, OMEGA_M0_BRIDGE
from dynamics.solver import SolverResult


# ---------------------------------------------------------------------------
# Result containers
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class TauMinimum:
    """Location and value of the τ(z) minimum."""
    z: float
    tau: float
    C_ratio: float       # C/C₁ at the minimum (= C_peak)
    omega_m: float       # Ω_m at the minimum  (≈ 2/3 at q=0)
    omega_lambda: float


@dataclass(frozen=True)
class OverflowWindow:
    """
    A redshift interval where C(τ)/C₁ exceeds a threshold.

    z_lo : lower boundary (C/C₁ crosses threshold going up, toward higher z)
    z_hi : upper boundary (C/C₁ crosses threshold going back down)
    threshold : C/C₁ value defining the window
    label : 'weak'  (threshold = 1.0)   or
            'strong'  (threshold = C_sphere/C₁)
    """
    z_lo: float
    z_hi: float
    threshold: float
    label: str

    @property
    def width(self) -> float:
        return self.z_hi - self.z_lo


@dataclass
class CriticalPoints:
    """All automatically-detected critical points for one solver run."""
    tau_minimum: TauMinimum
    overflow_weak: OverflowWindow | None
    overflow_strong: OverflowWindow | None

    # convenience
    @property
    def z_peak(self) -> float:
        return self.tau_minimum.z

    @property
    def C_peak_ratio(self) -> float:
        return self.tau_minimum.C_ratio


# ---------------------------------------------------------------------------
# Detection functions
# ---------------------------------------------------------------------------

def find_tau_minimum(result: SolverResult) -> TauMinimum:
    """
    Locate the global minimum of τ(z) in the solver result.

    τ(z) is non-monotone: it decreases from 1 at z=0, reaches τ_min at
    the acceleration crossing, then increases.  The minimum is detected
    by scanning the discrete array and then refining with a root-finder
    on dτ/dz = 0 (i.e. the spline derivative equals zero).

    Returns
    -------
    TauMinimum with the refined z, τ, C/C₁, Ω_m, Ω_Λ.
    """
    tau_arr = result.tau_arr
    z_arr = result.z_arr

    # Coarse location from array minimum
    idx_min = int(np.argmin(tau_arr))

    # Guard: if minimum is at the boundary, return boundary value
    if idx_min == 0 or idx_min == len(z_arr) - 1:
        z_min = float(z_arr[idx_min])
    else:
        # Refine: find z where d(tau)/dz = 0 using the cubic spline derivative
        dtau_spline = result.tau_of_z.derivative()

        # Search bracket around coarse index
        z_lo = float(z_arr[max(0, idx_min - 50)])
        z_hi = float(z_arr[min(len(z_arr) - 1, idx_min + 50)])

        try:
            z_min = brentq(dtau_spline, z_lo, z_hi, xtol=1e-8, rtol=1e-10)
        except ValueError:
            # Brentq failed (no sign change) — fall back to coarse value
            z_min = float(z_arr[idx_min])

    tau_min = float(result.tau_of_z(z_min))
    state = result.state_at(z_min)

    return TauMinimum(
        z=z_min,
        tau=tau_min,
        C_ratio=state.C_ratio,
        omega_m=state.omega_m,
        omega_lambda=state.omega_lambda,
    )


def find_overflow_window(
    result: SolverResult,
    threshold_ratio: float,
    label: str,
) -> OverflowWindow | None:
    """
    Find the redshift interval where C(τ)/C₁ > threshold_ratio.

    Uses the cubic spline interpolant for sub-grid accuracy.

    Parameters
    ----------
    result          : SolverResult from solve_history()
    threshold_ratio : C/C₁ threshold (1.0 for weak, C_sphere/C₁ for strong)
    label           : 'weak' or 'strong'

    Returns
    -------
    OverflowWindow or None if the threshold is never exceeded.
    """
    C_ratio_arr = result.C_ratio_arr
    z_arr = result.z_arr

    # Check whether the threshold is ever exceeded
    if np.max(C_ratio_arr) <= threshold_ratio:
        return None

    C_ratio_spline = result.C_ratio_of_z

    def f(z: float) -> float:
        return float(C_ratio_spline(z)) - threshold_ratio

    # Find z_lo: first crossing (C/C₁ rising above threshold)
    # C/C₁ starts at 1.0 at z=0 for the weak window, or below C_sphere/C1
    # for the strong window — scan upward from z=0
    z_lo = _find_first_crossing(f, z_arr, direction="up")
    z_hi = _find_first_crossing(f, z_arr, direction="down")

    if z_lo is None or z_hi is None:
        return None

    return OverflowWindow(
        z_lo=z_lo,
        z_hi=z_hi,
        threshold=threshold_ratio,
        label=label,
    )


def _find_first_crossing(
    f: callable,
    z_arr: np.ndarray,
    direction: str,
) -> float | None:
    """
    Find the z where f(z) changes sign.

    direction='up'  : find z where f goes from negative to positive
                      (C/C₁ rising above threshold)
    direction='down': find z where f goes from positive to negative
                      (C/C₁ falling below threshold after the peak)
    """
    f_vals = np.array([f(z) for z in z_arr])

    if direction == "up":
        # Look for first sign change from - to +
        for i in range(len(z_arr) - 1):
            if f_vals[i] <= 0 and f_vals[i + 1] > 0:
                try:
                    return brentq(f, z_arr[i], z_arr[i + 1], xtol=1e-8)
                except ValueError:
                    continue
    else:
        # Look for sign change from + to - after the peak
        # Find the peak first
        peak_idx = int(np.argmax(f_vals))
        for i in range(peak_idx, len(z_arr) - 1):
            if f_vals[i] >= 0 and f_vals[i + 1] < 0:
                try:
                    return brentq(f, z_arr[i], z_arr[i + 1], xtol=1e-8)
                except ValueError:
                    continue

    return None


def detect_all(result: SolverResult) -> CriticalPoints:
    """
    Run all critical-point detectors on a SolverResult.

    This is the main entry point for the detection pipeline.

    Returns
    -------
    CriticalPoints with tau_minimum, overflow_weak, overflow_strong.
    """
    tau_min = find_tau_minimum(result)

    C_sphere_ratio = C_SPHERE / C1   # ≈ 1.304 (f_sph = √(4π)/e)

    overflow_weak = find_overflow_window(
        result,
        threshold_ratio=1.0,
        label="weak",
    )
    overflow_strong = find_overflow_window(
        result,
        threshold_ratio=C_sphere_ratio,
        label="strong",
    )

    return CriticalPoints(
        tau_minimum=tau_min,
        overflow_weak=overflow_weak,
        overflow_strong=overflow_strong,
    )
