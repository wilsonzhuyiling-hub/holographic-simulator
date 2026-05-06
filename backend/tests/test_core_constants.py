"""Sanity tests for the legacy core module.

These tests lock the numerical values of the framework's primary derived
constants. If a refactor changes any of these values, the test fails — this
is intentional. The framework is zero-free-parameter; numerical values are
predictions, not configuration.
"""

import math

import pytest

from core.universe_simulator import (
    Axioms,
    PhysicalInputs,
    derive_C0_axis,
    derive_C1,
    derive_constants,
    derive_f_sph,
    derive_kappa_geom,
    derive_lambda_capacity,
    derive_omega_dm,
    derive_omega_lambda_fp,
    tau_from_redshift,
    redshift_from_tau,
    T_CMB_from_recombination,
    wz_ratio_from_stability_inversion,
)


class TestFundamentalConstants:
    """Lock the geometric constants from P1."""

    def test_C0_axis(self):
        # C₀ = π²/16, the unsuppressed mirror-axis value
        assert derive_C0_axis() == pytest.approx(math.pi**2 / 16, rel=1e-12)

    def test_C1_fixed_point(self):
        # C₁ = π²/(16e²), de Sitter fixed point (P1 main result)
        expected = math.pi**2 / (16 * math.e**2)
        assert derive_C1() == pytest.approx(expected, rel=1e-12)
        assert derive_C1() == pytest.approx(0.083485, rel=1e-4)

    def test_omega_lambda_fixed_point(self):
        # Ω_Λ_fp = π³/(6e²) ≈ 0.6994
        assert derive_omega_lambda_fp() == pytest.approx(0.69937, rel=1e-4)

    def test_f_sph(self):
        # f_sph = √(4π)/e ≈ 1.304
        assert derive_f_sph() == pytest.approx(math.sqrt(4 * math.pi) / math.e, rel=1e-12)

    def test_omega_dm_base(self):
        # Ω_DM_base = π⁴(√(4π) - e)/(6e⁴) ≈ 0.2458
        assert derive_omega_dm() == pytest.approx(0.24580, rel=1e-4)

    def test_kappa_geom(self):
        # κ_geom ≈ 1.0563 from P5.1 (recorded in README v1.0)
        assert derive_kappa_geom() == pytest.approx(1.0563113, rel=1e-5)

    def test_lambda_capacity(self):
        # λ = 1/(f_sph·e²)
        f_sph = math.sqrt(4 * math.pi) / math.e
        expected = 1 / (f_sph * math.e**2)
        assert derive_lambda_capacity() == pytest.approx(expected, rel=1e-12)


class TestRedshiftMapping:
    """Pure de Sitter τ ↔ z relation from P2 Eq. (5)."""

    def test_z_zero_corresponds_to_tau_one(self):
        # τ(z=0) = 1 by anchor convention
        assert tau_from_redshift(0.0) == pytest.approx(1.0, rel=1e-12)

    def test_inverse_consistency(self):
        # τ → z → τ round-trip
        for tau in [0.1, 0.5, 1.0, 1.5, 2.0]:
            z = redshift_from_tau(tau)
            assert tau_from_redshift(z) == pytest.approx(tau, rel=1e-10)

    def test_pure_desitter_formula(self):
        # τ(z) = 1 - ln(1+z)
        for z in [0.1, 1.0, 10.0, 100.0]:
            expected = 1 - math.log(1 + z)
            assert tau_from_redshift(z) == pytest.approx(expected, rel=1e-12)


class TestObservationalBridges:
    """Bridge quantities that compare against external measurements."""

    def test_T_CMB_close_to_observed(self):
        inputs = PhysicalInputs()
        T_predicted = T_CMB_from_recombination(inputs.T_rec_K, inputs.z_rec)
        # Within 0.1% of Fixsen 2009 value
        assert T_predicted == pytest.approx(2.7255, rel=1e-3)

    def test_wz_ratio_within_2_percent(self):
        # m_W/m_Z prediction should be within ~2% of PDG
        ratio = wz_ratio_from_stability_inversion()
        assert ratio == pytest.approx(0.88136, rel=2e-2)


class TestDerivedConstantsBundle:
    """The full derived-constants object should self-assemble correctly."""

    def test_bundle_assembles(self):
        constants = derive_constants(Axioms(), PhysicalInputs())
        assert constants.C1 > 0
        assert constants.C0_axis > constants.C1  # C₀ > C₁ since e² > 1
        assert constants.C_sphere > constants.C1  # spherical opens above fixed point
        assert 0 < constants.lambda_capacity < 1
        assert constants.kappa_geom > 1
