"""Tests for the Stage 2 dark-matter diagnostic curve."""

import numpy as np
import pytest

from dynamics import (
    KAPPA_GEOM,
    OMEGA_DM_BASE,
    OMEGA_DM_FP,
    compute_dark_matter_curve,
    compute_omega_dm_eff_demo,
    solve_history,
)


@pytest.fixture(scope="module")
def dm_curve():
    result = solve_history(z_max=3.5, n_eval=5_000)
    return compute_dark_matter_curve(result)


class TestDarkMatterConstants:
    def test_omega_dm_base(self):
        assert pytest.approx(0.24580, rel=1e-4) == OMEGA_DM_BASE

    def test_kappa_geom(self):
        assert pytest.approx(1.0563113, rel=1e-5) == KAPPA_GEOM

    def test_omega_dm_fp(self):
        assert pytest.approx(KAPPA_GEOM * OMEGA_DM_BASE, rel=1e-12) == OMEGA_DM_FP


class TestRelaxationWindow:
    def test_boundaries_are_detected_from_ode(self, dm_curve):
        assert dm_curve.z_peak == pytest.approx(0.898, rel=2e-3)
        assert 2.83 < dm_curve.z_relaxation_end < 2.85

    def test_loss_zero_outside_relaxation_interval(self, dm_curve):
        assert float(dm_curve.relaxation_loss_fraction_of_z(0.1)) == pytest.approx(0.0, abs=1e-8)
        assert float(dm_curve.relaxation_loss_fraction_of_z(3.0)) == pytest.approx(0.0, abs=1e-8)

    def test_loss_is_positive_inside_relaxation_interval(self, dm_curve):
        z_mid = 0.5 * (dm_curve.z_peak + dm_curve.z_relaxation_end)
        assert float(dm_curve.relaxation_loss_fraction_of_z(z_mid)) > 0.0

    def test_omega_dm_returns_to_kappa_corrected_value(self, dm_curve):
        assert float(dm_curve.omega_dm_of_z(0.1)) == pytest.approx(OMEGA_DM_FP, rel=1e-8)
        assert float(dm_curve.omega_dm_of_z(3.0)) == pytest.approx(OMEGA_DM_FP, rel=1e-8)

    def test_relaxation_is_a_downward_correction(self, dm_curve):
        assert np.min(dm_curve.omega_dm_arr) < OMEGA_DM_FP


class TestOmegaDMEffDemo:
    def test_demo_samples_effective_redshifts(self):
        result = solve_history(z_max=3.5, n_eval=5_000)
        z_eff = np.array([0.3, 0.93, 1.49, 2.33])
        demo = compute_omega_dm_eff_demo(result, z_eff=z_eff, shell_width=0.2)

        np.testing.assert_array_equal(demo.z_eff_arr, z_eff)
        assert len(demo.omega_dm_eff_arr) == len(z_eff)
        assert np.all(demo.z_lo_arr < demo.z_eff_arr)
        assert np.all(demo.z_hi_arr > demo.z_eff_arr)

    def test_demo_values_are_lightcone_samples_of_continuous_curve(self):
        result = solve_history(z_max=3.5, n_eval=5_000)
        demo = compute_omega_dm_eff_demo(result, z_eff=np.array([0.93, 1.49]))

        for z, omega_eff in zip(demo.z_eff_arr, demo.omega_dm_eff_arr, strict=True):
            assert omega_eff == pytest.approx(
                float(demo.continuous_curve.omega_dm_of_z(z)),
                rel=1e-12,
            )
