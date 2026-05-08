"""ODE solver for the matter-deformed holographic framework.

Integrates the self-consistent system from P2 Sec. 3 using scipy's
RK45 adaptive-step solver. Produces two continuous curves:

    Physical history  : τ(z), Ω_Λ(z), Ω_m(z), C/C₁  vs. z
    (z runs forward from 0 to z_max; τ starts at 1 and evolves)

The solver returns a dense output (continuous interpolant) so the curves
can be evaluated at arbitrary z without re-solving.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

from dynamics.ode_system import (
    OMEGA_LAMBDA_FP,
    OMEGA_M1,
    C1,
    C_SPHERE,
    ODEState,
    build_state,
    dtau_dz,
    C_ratio,
    _omega_lambda,
    _omega_m,
)


# ---------------------------------------------------------------------------
# Default integration range
# ---------------------------------------------------------------------------

Z_MIN: float = 1e-4     # avoid z=0 singularity; τ(z=0) = 1 is the BC
Z_MAX: float = 5.0      # well past the overflow window (z ≈ 2.8 upper edge)
N_EVAL: int = 5_000     # number of evaluation points for the dense output


@dataclass
class SolverResult:
    """
    Output of solve_history().

    Contains the raw scipy result plus convenience arrays and an
    interpolant for continuous evaluation.
    """
    z_arr: np.ndarray           # redshift grid
    tau_arr: np.ndarray         # τ(z)
    omega_lambda_arr: np.ndarray
    omega_m_arr: np.ndarray
    C_ratio_arr: np.ndarray     # C(τ)/C₁

    # Continuous interpolants (CubicSpline, callable)
    tau_of_z: CubicSpline = field(repr=False)
    C_ratio_of_z: CubicSpline = field(repr=False)
    omega_lambda_of_z: CubicSpline = field(repr=False)
    omega_m_of_z: CubicSpline = field(repr=False)

    # Solver metadata
    success: bool = True
    message: str = ""
    n_steps: int = 0
    omega_m0: float = OMEGA_M1   # matter fraction used during integration

    def state_at(self, z: float) -> ODEState:
        """Evaluate all diagnostics at an arbitrary redshift z."""
        tau = float(self.tau_of_z(z))
        return build_state(z, tau, self.omega_m0)

    def H_normalized(self, z: float, omega_m0: float | None = None) -> float:
        """
        E(z) = H(z)/H₀ — the normalized Hubble parameter.

            E²(z) = Ω_m1(1+z)³ + Ω_Λ(z)

        Uses the omega_m0 stored at solve time by default.
        Pass omega_m0 explicitly only for sensitivity comparisons.
        """
        if omega_m0 is None:
            omega_m0 = self.omega_m0
        omega_lam = float(self.omega_lambda_of_z(z))
        return math.sqrt(omega_m0 * (1.0 + z) ** 3 + omega_lam)


def solve_history(
    z_max: float = Z_MAX,
    n_eval: int = N_EVAL,
    omega_m0: float = OMEGA_M1,
    rtol: float = 1e-9,
    atol: float = 1e-11,
) -> SolverResult:
    """
    Solve dτ/dz from z=0 to z=z_max with boundary condition τ(0) = 1.

    Parameters
    ----------
    z_max   : maximum redshift to integrate to
    n_eval  : number of points in the output grid
    omega_m0: matter fraction bridge input (P2 convention)
    rtol    : relative tolerance for RK45
    atol    : absolute tolerance for RK45

    Returns
    -------
    SolverResult with continuous interpolants and diagnostic arrays
    """
    tau0 = 1.0   # τ(z=0) = 1, the de Sitter fixed-point anchor  [P2 BC]

    def rhs(z: float, y: list[float]) -> list[float]:
        return [dtau_dz(z, y[0], omega_m0)]

    z_span = (Z_MIN, z_max)
    z_eval = np.linspace(Z_MIN, z_max, n_eval)

    # Propagate the IC from z=0 to z=Z_MIN analytically (pure dS limit)
    # τ(z) = 1 - ln(1+z) for small z — negligible error
    tau_at_zmin = 1.0 - math.log1p(Z_MIN)

    sol = solve_ivp(
        rhs,
        z_span,
        [tau_at_zmin],
        method="RK45",
        t_eval=z_eval,
        rtol=rtol,
        atol=atol,
        dense_output=False,
    )

    if not sol.success:
        raise RuntimeError(f"ODE solver failed: {sol.message}")

    tau_arr = sol.y[0]
    z_arr = sol.t

    # Prepend the z=0 anchor point
    z_arr = np.concatenate([[0.0], z_arr])
    tau_arr = np.concatenate([[1.0], tau_arr])

    # Derived arrays
    omega_lambda_arr = np.array([
        float(_omega_lambda(tau)) for tau in tau_arr
    ])
    omega_m_arr = np.array([
        float(_omega_m(z, ol, omega_m0))
        for z, ol in zip(z_arr, omega_lambda_arr)
    ])
    C_ratio_arr = np.exp(-2.0 * (tau_arr - 1.0))

    # Build continuous interpolants (cubic spline, monotone z grid)
    tau_of_z = CubicSpline(z_arr, tau_arr)
    C_ratio_of_z = CubicSpline(z_arr, C_ratio_arr)
    omega_lambda_of_z = CubicSpline(z_arr, omega_lambda_arr)
    omega_m_of_z = CubicSpline(z_arr, omega_m_arr)

    return SolverResult(
        z_arr=z_arr,
        tau_arr=tau_arr,
        omega_lambda_arr=omega_lambda_arr,
        omega_m_arr=omega_m_arr,
        C_ratio_arr=C_ratio_arr,
        tau_of_z=tau_of_z,
        C_ratio_of_z=C_ratio_of_z,
        omega_lambda_of_z=omega_lambda_of_z,
        omega_m_of_z=omega_m_of_z,
        success=sol.success,
        message=sol.message,
        n_steps=sol.t.shape[0],
        omega_m0=omega_m0,
    )
