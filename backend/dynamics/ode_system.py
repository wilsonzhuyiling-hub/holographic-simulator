"""ODE system for the matter-deformed holographic framework.

Implements the self-consistent system from P2 Sec. 3 (Eqs. 10–13):

    dτ/dz = -(1 - ε(z)) / (1 + z)          [P2 Eq. 10]
    ε(z)  = (3/2) Ω_m(z)                    [P2 Eq. 11]
    Ω_m(z) = Ω_m0(1+z)³ / D(z)             [P2 Eq. 12]
    Ω_Λ(z) = Ω_Λ0 · exp[-2(τ(z) - 1)]     [P2 Eq. 13]
    D(z)   = Ω_m0(1+z)³ + Ω_Λ(z)

Boundary condition: τ(z=0) = 1  (de Sitter fixed point anchor)

Parameter classification
------------------------
Model-endogenous (zero free parameters):
    Ω_Λ0 = π³/(6e²)  — derived from geometry (P1/P2)

Validation input (hardcoded, not a fit knob):
    Ω_m0 = 0.3  — observational background input, labelled as bridge
"""

from __future__ import annotations

import math
from typing import NamedTuple


# ---------------------------------------------------------------------------
# Model-endogenous constants (from core, duplicated here for zero-import
# dependency so dynamics/ can be tested standalone)
# ---------------------------------------------------------------------------

_PI = math.pi
_E = math.e

C1: float = _PI**2 / (16 * _E**2)
"""de Sitter fixed-point coefficient C₁ = π²/(16e²)  [P1 main result]"""

OMEGA_LAMBDA_FP: float = _PI**3 / (6 * _E**2)
"""Dark-energy fixed-point fraction Ω_Λ0 = π³/(6e²)  [P2 Eq. 9]"""

OMEGA_M1: float = 1.0 - OMEGA_LAMBDA_FP
"""
Matter fraction at the C₁ fixed point: Ω_m1 = 1 − π³/(6e²) ≈ 0.3006.

Tier 1 geometric constant.  Derived from flat-universe closure at τ = 1:
    Ω_m1 + Ω_Λ0 = 1   →   Ω_m1 = 1 − Ω_Λ0

Guarantees E(z=0) = 1 exactly.  This is the correct default matter input
for the ODE system; it replaces the earlier OMEGA_M0_BRIDGE = 0.3 round
value, which caused a 0.063% Friedmann closure error.

The name carries the subscript '1' because it is the matter complement of
C₁, analogous to Ω_Λ0 = π³/(6e²) being the dark-energy reading at C₁.
"""

C_SPHERE: float = C1 * math.sqrt(4 * _PI) / _E
"""First spherical release threshold C_sphere = C₁·f_sph  [P2 Eq. 16]"""

# ---------------------------------------------------------------------------
# Validation input — observational bridge, NOT a fit knob
# ---------------------------------------------------------------------------

OMEGA_M0_BRIDGE: float = 0.3
"""
Deprecated — use OMEGA_M1 instead.
Legacy round value used in early drafts of P2 (when the fixed point was
labelled C₀ rather than C₁).  Retained for reference only; causes a 0.063%
Friedmann closure error.  Not used as a default anywhere in the codebase.
"""

OMEGA_M0_PLANCK: float = 0.315
"""
Planck 2018 observed total matter fraction (Tier 3 reference).
Planck 2018 TT,TE,EE+lowE+lensing: Ω_m = 0.315 ± 0.007.
The framework's geometric prediction OMEGA_M1 ≈ 0.3006 differs from this
by ~0.5%, providing a falsifiable comparison.
"""


# ---------------------------------------------------------------------------
# State vector layout
# ---------------------------------------------------------------------------

class ODEState(NamedTuple):
    """State vector at a single redshift z."""
    z: float
    tau: float
    omega_lambda: float
    omega_m: float
    C_ratio: float       # C(τ)/C₁  — convenience, derived from tau
    epsilon: float       # ε(z) = (3/2)Ω_m(z)  — deceleration departure


# ---------------------------------------------------------------------------
# ODE right-hand side
# ---------------------------------------------------------------------------

def dtau_dz(z: float, tau: float, omega_m0: float = OMEGA_M1) -> float:
    """
    dτ/dz = -(1 - ε(z)) / (1 + z)    [P2 Eq. 10]

    where ε(z) = (3/2) Ω_m(z) and Ω_m(z) is evaluated self-consistently
    using the current τ to compute Ω_Λ(z).

    Parameters
    ----------
    z       : redshift
    tau     : current τ value (the ODE state variable)
    omega_m0: matter fraction today (validation bridge input)

    Returns
    -------
    dτ/dz at this point
    """
    omega_lam = _omega_lambda(tau)
    om = _omega_m(z, omega_lam, omega_m0)
    eps = _epsilon(om)
    return -(1.0 - eps) / (1.0 + z)


def _omega_lambda(tau: float) -> float:
    """Ω_Λ(z) = Ω_Λ0 · exp[-2(τ - 1)]  [P2 Eq. 13]"""
    return OMEGA_LAMBDA_FP * math.exp(-2.0 * (tau - 1.0))


def _omega_m(z: float, omega_lambda: float, omega_m0: float) -> float:
    """
    Ω_m(z) = Ω_m0(1+z)³ / [Ω_m0(1+z)³ + Ω_Λ(z)]  [P2 Eq. 12]

    Clamped to [0, 1] to guard against numerical excursions at extreme z.
    """
    matter_term = omega_m0 * (1.0 + z) ** 3
    denom = matter_term + omega_lambda
    if denom <= 0.0:
        return 0.0
    return min(matter_term / denom, 1.0)


def _epsilon(omega_m: float) -> float:
    """ε(z) = (3/2) Ω_m(z)  [P2 Eq. 11]"""
    return 1.5 * omega_m


def C_geom(tau: float) -> float:
    """
    C(τ) = (π²/16) · e^(-2τ)

    The geometric coefficient at the real axis (τ_I = 0).
    """
    return (_PI**2 / 16.0) * math.exp(-2.0 * tau)


def C_ratio(tau: float) -> float:
    """C(τ) / C₁ — supersaturation level relative to de Sitter fixed point."""
    return math.exp(-2.0 * (tau - 1.0))


def build_state(z: float, tau: float, omega_m0: float = OMEGA_M1) -> ODEState:
    """Compute the full diagnostic state at a given (z, τ) point."""
    omega_lam = _omega_lambda(tau)
    om = _omega_m(z, omega_lam, omega_m0)
    eps = _epsilon(om)
    return ODEState(
        z=z,
        tau=tau,
        omega_lambda=omega_lam,
        omega_m=om,
        C_ratio=C_ratio(tau),
        epsilon=eps,
    )
