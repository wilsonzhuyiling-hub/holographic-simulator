"""Dynamics engine — Stage 1/2.

ODE-based integration of the matter-deformed framework equations.
Implements P2 Sec. 3 self-consistent system:

    dτ/dz = -(1 - ε(z)) / (1 + z)       [P2 Eq. 10]
    ε(z)  = (3/2) Ω_m(z)                 [P2 Eq. 11]
    Ω_m(z) = Ω_m0(1+z)³ / D(z)          [P2 Eq. 12]
    Ω_Λ(z) = Ω_Λ0 exp[-2(τ(z) - 1)]    [P2 Eq. 13]

with boundary condition τ(z=0) = 1 and Ω_m0 = 0.3 (validation bridge input).

Public API
----------
    solve_history()        → SolverResult      (continuous interpolants)
    detect_all()           → CriticalPoints    (τ_min, C_peak, overflow windows)
    compute_observables()  → ObservablesResult (BAO distances, ratio observables)
"""

from dynamics.solver import solve_history, SolverResult
from dynamics.critical_points import detect_all, CriticalPoints, TauMinimum, OverflowWindow
from dynamics.observables import (
    compute_observables,
    ObservablesResult,
    DH0,
    H0_PLANCK,
    R_D_PLANCK,
    R_D_DESI,
    C_KMS,
)
from dynamics.ode_system import (
    C1,
    C_SPHERE,
    OMEGA_LAMBDA_FP,
    OMEGA_M1,
    OMEGA_M0_BRIDGE,
    OMEGA_M0_PLANCK,
    ODEState,
    build_state,
    C_geom,
    C_ratio,
)

__all__ = [
    "solve_history",
    "SolverResult",
    "detect_all",
    "CriticalPoints",
    "TauMinimum",
    "OverflowWindow",
    "compute_observables",
    "ObservablesResult",
    "DH0",
    "H0_PLANCK",
    "R_D_PLANCK",
    "R_D_DESI",
    "C_KMS",
    "C1",
    "C_SPHERE",
    "OMEGA_LAMBDA_FP",
    "OMEGA_M1",
    "OMEGA_M0_BRIDGE",
    "OMEGA_M0_PLANCK",
    "ODEState",
    "build_state",
    "C_geom",
    "C_ratio",
]
