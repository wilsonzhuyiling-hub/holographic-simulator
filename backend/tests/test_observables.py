"""Tests for BAO distance observables — Stage 2.

Validates that compute_observables() produces physically consistent distance
curves from the dynamics engine output.  Paper target values are NOT tested
here — those belong in the statistics layer once DESI data is loaded.

What is tested
--------------
- Boundary conditions (D_C=0 and D_H=DH0 at z=0)
- Flat-universe identity (D_M = D_C)
- Monotonicity of D_C, D_M
- Positivity of D_V for z > 0
- Consistency of BAO ratio interpolants with the distance arrays
- Tier 3 anchors are recorded on the result object
"""

import pytest
import numpy as np

from dynamics import solve_history, compute_observables
from dynamics.observables import DH0, H0_PLANCK, R_D_PLANCK, C_KMS


@pytest.fixture(scope="module")
def obs():
    result = solve_history(z_max=3.5, n_eval=5_000)
    return compute_observables(result)


class TestConstants:

    def test_DH0_definition(self):
        assert DH0 == pytest.approx(C_KMS / H0_PLANCK, rel=1e-9)

    def test_DH0_magnitude(self):
        assert 4400 < DH0 < 4510


class TestBoundaryConditions:

    def test_D_C_zero_at_z0(self, obs):
        assert float(obs.D_C_of_z(0.0)) == pytest.approx(0.0, abs=1.0)

    def test_D_M_zero_at_z0(self, obs):
        assert float(obs.D_M_of_z(0.0)) == pytest.approx(0.0, abs=1.0)

    def test_D_H_at_z0_equals_DH0(self, obs):
        assert float(obs.D_H_of_z(0.0)) == pytest.approx(DH0, rel=1e-3)

    def test_D_V_zero_at_z0(self, obs):
        assert obs.D_V_arr[0] == pytest.approx(0.0, abs=1e-6)


class TestMonotonicity:

    def test_D_C_monotone_increasing(self, obs):
        assert np.all(np.diff(obs.D_C_arr) > 0)

    def test_D_M_monotone_increasing(self, obs):
        assert np.all(np.diff(obs.D_M_arr) > 0)

    def test_D_H_decreasing_at_low_z(self, obs):
        # H(z) increases → D_H = c/H(z) decreases, at least up to z ~ 1
        mask = obs.z_arr <= 1.0
        D_H_low = obs.D_H_arr[mask]
        assert np.all(np.diff(D_H_low) < 0)

    def test_D_V_positive_for_z_gt_0(self, obs):
        assert np.all(obs.D_V_arr[obs.z_arr > 0.0] > 0)


class TestFlatUniverse:

    def test_D_M_equals_D_C(self, obs):
        np.testing.assert_array_equal(obs.D_M_arr, obs.D_C_arr)


class TestBAORatioConsistency:

    def test_DM_rd_interpolant_matches_array(self, obs):
        z = 1.0
        assert float(obs.DM_over_rd_of_z(z)) == pytest.approx(
            float(obs.D_M_of_z(z)) / R_D_PLANCK, rel=1e-4
        )

    def test_DH_rd_interpolant_matches_array(self, obs):
        z = 1.0
        assert float(obs.DH_over_rd_of_z(z)) == pytest.approx(
            float(obs.D_H_of_z(z)) / R_D_PLANCK, rel=1e-4
        )

    def test_DV_rd_interpolant_matches_array(self, obs):
        z = 0.5
        assert float(obs.DV_over_rd_of_z(z)) == pytest.approx(
            float(obs.D_V_of_z(z)) / R_D_PLANCK, rel=1e-4
        )


class TestBAORatioSanity:

    def test_DM_rd_at_z1_cosmological_range(self, obs):
        # Rough sanity: D_M/r_d at z=1 should be in [15, 25] for any
        # reasonable flat cosmology with H₀ ~ 67, Ω_m ~ 0.3
        ratio = float(obs.DM_over_rd_of_z(1.0))
        assert 15 < ratio < 25

    def test_DH_rd_at_z1_cosmological_range(self, obs):
        ratio = float(obs.DH_over_rd_of_z(1.0))
        assert 10 < ratio < 20

    def test_DV_rd_at_z05_cosmological_range(self, obs):
        ratio = float(obs.DV_over_rd_of_z(0.5))
        assert 6 < ratio < 15


class TestTierThreeProvenance:

    def test_H0_recorded(self, obs):
        assert obs.H0 == pytest.approx(H0_PLANCK)

    def test_rd_recorded(self, obs):
        assert obs.r_d == pytest.approx(R_D_PLANCK)
