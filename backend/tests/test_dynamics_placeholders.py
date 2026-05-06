"""Tests for the dynamics engine — Stage 1.

These tests validate that the ODE solver recovers the critical-point values
reported in P2 Sec. 3.3.  The paper values are the *targets*; the code must
discover them numerically without any hard-coding.

Expected values (P2 paper)
--------------------------
τ_min ≈ 0.808          at z ≈ 0.898
C_peak / C₁ ≈ 1.467   at the same z
Acceleration crossing: Ω_m(z_min) ≈ 2/3  (q = 0 condition)
Strong overflow window: z ∈ [0.345, 1.721]
Weak overflow window:   z ∈ [0.001, 2.846]
"""

import math
import pytest
import numpy as np

from dynamics import solve_history, detect_all
from dynamics.ode_system import C1, C_SPHERE, OMEGA_LAMBDA_FP


@pytest.fixture(scope="module")
def result():
    return solve_history(z_max=5.0, n_eval=10_000)


@pytest.fixture(scope="module")
def critical(result):
    return detect_all(result)


class TestSolverBasics:

    def test_solver_succeeds(self, result):
        assert result.success

    def test_boundary_condition(self, result):
        tau_at_zero = float(result.tau_of_z(0.0))
        assert tau_at_zero == pytest.approx(1.0, abs=1e-4)

    def test_omega_lambda_at_z0(self, result):
        ol_at_zero = float(result.omega_lambda_of_z(0.0))
        assert ol_at_zero == pytest.approx(OMEGA_LAMBDA_FP, rel=1e-4)

    def test_tau_decreases_then_increases(self, result):
        tau_arr = result.tau_arr
        idx_min = int(np.argmin(tau_arr))
        assert 0 < idx_min < len(tau_arr) - 1

    def test_C_ratio_peaks_above_one(self, result):
        assert float(np.max(result.C_ratio_arr)) > 1.0

    def test_interpolants_are_callable(self, result):
        for z in [0.1, 0.5, 0.898, 1.5, 3.0]:
            assert math.isfinite(float(result.tau_of_z(z)))
            assert math.isfinite(float(result.C_ratio_of_z(z)))
            assert math.isfinite(float(result.omega_lambda_of_z(z)))
            assert math.isfinite(float(result.omega_m_of_z(z)))


class TestTauMinimum:

    def test_tau_minimum_recovered(self, critical):
        tau_min = critical.tau_minimum.tau
        z_at_min = critical.tau_minimum.z
        assert tau_min == pytest.approx(0.808, abs=0.005), (
            f"tau_min = {tau_min:.4f}, expected ~0.808"
        )
        assert z_at_min == pytest.approx(0.898, abs=0.015), (
            f"z(tau_min) = {z_at_min:.4f}, expected ~0.898"
        )

    def test_C_peak_ratio_recovered(self, critical):
        C_peak = critical.C_peak_ratio
        assert C_peak == pytest.approx(1.467, abs=0.01), (
            f"C_peak/C1 = {C_peak:.4f}, expected ~1.467"
        )

    def test_acceleration_crossing_coincides_with_tau_min(self, critical):
        omega_m_at_min = critical.tau_minimum.omega_m
        assert omega_m_at_min == pytest.approx(2.0 / 3.0, abs=0.02), (
            f"Omega_m at tau_min = {omega_m_at_min:.4f}, expected ~2/3"
        )


class TestOverflowWindows:

    def test_strong_overflow_window_boundaries(self, critical):
        strong = critical.overflow_strong
        assert strong is not None, "Strong overflow window not detected"
        assert strong.z_lo == pytest.approx(0.345, abs=0.02)
        assert strong.z_hi == pytest.approx(1.721, abs=0.03)

    def test_weak_overflow_window_boundaries(self, critical):
        weak = critical.overflow_weak
        assert weak is not None, "Weak overflow window not detected"
        assert weak.z_lo < 0.05
        assert weak.z_hi == pytest.approx(2.846, abs=0.05)

    def test_strong_window_inside_weak_window(self, critical):
        strong = critical.overflow_strong
        weak = critical.overflow_weak
        assert strong is not None and weak is not None
        assert strong.z_lo > weak.z_lo
        assert strong.z_hi < weak.z_hi

    def test_C_sphere_threshold_correct(self, critical):
        C_sphere_ratio = C_SPHERE / C1
        assert critical.overflow_strong.threshold == pytest.approx(C_sphere_ratio, rel=1e-6)


class TestHubbleNormalized:

    def test_H_at_z0_equals_one(self, result):
        E0 = result.H_normalized(0.0)
        assert E0 == pytest.approx(1.0, abs=0.01)

    def test_H_increases_with_z(self, result):
        z_grid = [0.0, 0.3, 0.6, 1.0, 2.0]
        E_vals = [result.H_normalized(z) for z in z_grid]
        for i in range(len(E_vals) - 1):
            assert E_vals[i] < E_vals[i + 1]
