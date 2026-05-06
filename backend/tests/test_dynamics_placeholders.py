"""Placeholder tests for the dynamics engine — Stage 1.

These tests are skipped until the dynamics module is implemented. They lock
the expected behaviour of the ODE solver against P2 Sec. 3 numerical results.

Expected milestones (from P2 paper):
- τ_min ≈ 0.808 at z ≈ 0.898
- C_peak / C₁ ≈ 1.467
- Acceleration crossing at the same z (q = 0)
- Strong overflow window: z ∈ [0.345, 1.721]
"""

import pytest

pytestmark = pytest.mark.skip(reason="Dynamics engine not yet implemented (Stage 1)")


def test_tau_minimum_recovered():
    """τ(z) should reach minimum ≈ 0.808 at z ≈ 0.898."""
    # from dynamics.solver import solve_history
    # history = solve_history()
    # tau_min, z_at_min = find_tau_minimum(history)
    # assert tau_min == pytest.approx(0.808, abs=0.005)
    # assert z_at_min == pytest.approx(0.898, abs=0.01)


def test_C_peak_ratio_recovered():
    """C/C₁ should peak at ≈ 1.467 at the same redshift as τ_min."""
    # ...


def test_acceleration_crossing_coincides_with_tau_min():
    """q(z) = 0 should occur at the same z as τ_min."""
    # ...


def test_strong_overflow_window_boundaries():
    """C > C_sphere region should span z ∈ [0.345, 1.721]."""
    # ...
