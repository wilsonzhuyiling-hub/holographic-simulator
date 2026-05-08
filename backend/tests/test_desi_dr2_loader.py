"""Tests for the DESI DR2 BAO loader."""

import numpy as np
import pytest

from observation import load_desi_dr2_bao


@pytest.fixture(scope="module")
def desi():
    return load_desi_dr2_bao()


class TestDESIDR2Loader:
    def test_loads_all_consensus_measurements(self, desi):
        assert len(desi.measurements) == 13
        assert desi.covariance.shape == (13, 13)

    def test_covariance_diagonal_matches_sigmas(self, desi):
        sigmas = np.sqrt(np.diag(desi.covariance))
        row_sigmas = np.array([row.sigma for row in desi.measurements])
        np.testing.assert_allclose(row_sigmas, sigmas)

    def test_quantities_are_bao_observables(self, desi):
        assert set(desi.quantities) == {"DM_over_rs", "DH_over_rs", "DV_over_rs"}

    def test_bgs_dv_is_single_low_redshift_point(self, desi):
        dv_rows = [row for row in desi.measurements if row.quantity == "DV_over_rs"]
        assert len(dv_rows) == 1
        assert dv_rows[0].z == pytest.approx(0.295)
        assert dv_rows[0].value == pytest.approx(7.94167639)

    def test_unique_redshifts_match_desi_shells(self, desi):
        np.testing.assert_allclose(
            desi.unique_redshifts,
            np.array([0.295, 0.510, 0.706, 0.934, 1.321, 1.484, 2.33]),
        )
