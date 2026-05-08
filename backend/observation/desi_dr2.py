"""DESI DR2 BAO data loader.

The bundled files are the DESI DR2 all-sample Gaussian BAO consensus from
CobayaSampler/bao_data, sourced from the DESI DR2 BAO publications.
They contain BAO observables as functions of effective redshift:

    D_M / r_s, D_H / r_s, and D_V / r_s

They are not direct Omega_DM measurements.  Downstream comparisons should keep
the BAO observables on their own axis and use the redshifts only to sample
Omega_DM_eff(z).
"""

from __future__ import annotations

from dataclasses import dataclass
from importlib.resources import files

import numpy as np

_DATA_DIR = files("observation") / "data" / "desi_dr2"
_MEAN_FILE = "desi_gaussian_bao_ALL_GCcomb_mean.txt"
_COV_FILE = "desi_gaussian_bao_ALL_GCcomb_cov.txt"


@dataclass(frozen=True)
class DESIBAOMeasurement:
    """One DESI DR2 BAO measurement at an effective redshift."""

    z: float
    value: float
    quantity: str
    sigma: float


@dataclass(frozen=True)
class DESIDR2BAOData:
    """DESI DR2 BAO mean vector and covariance matrix."""

    measurements: tuple[DESIBAOMeasurement, ...]
    covariance: np.ndarray
    source: str

    @property
    def z_eff_arr(self) -> np.ndarray:
        """Effective redshifts in the same order as the mean vector."""
        return np.array([row.z for row in self.measurements])

    @property
    def values(self) -> np.ndarray:
        """BAO mean vector in the same order as the covariance matrix."""
        return np.array([row.value for row in self.measurements])

    @property
    def quantities(self) -> tuple[str, ...]:
        """Observable labels in the same order as the mean vector."""
        return tuple(row.quantity for row in self.measurements)

    @property
    def unique_redshifts(self) -> np.ndarray:
        """Sorted unique effective redshifts."""
        return np.array(sorted({row.z for row in self.measurements}))


def load_desi_dr2_bao() -> DESIDR2BAOData:
    """Load the bundled DESI DR2 all-sample BAO mean vector and covariance."""
    mean_path = _DATA_DIR / _MEAN_FILE
    cov_path = _DATA_DIR / _COV_FILE

    z_vals: list[float] = []
    values: list[float] = []
    quantities: list[str] = []

    with mean_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            z_str, value_str, quantity = stripped.split()
            z_vals.append(float(z_str))
            values.append(float(value_str))
            quantities.append(quantity)

    covariance = np.loadtxt(cov_path)
    if covariance.shape != (len(values), len(values)):
        raise ValueError(
            f"DESI covariance has shape {covariance.shape}, expected {(len(values), len(values))}."
        )

    sigmas = np.sqrt(np.diag(covariance))
    measurements = tuple(
        DESIBAOMeasurement(z=z, value=value, quantity=quantity, sigma=float(sigma))
        for z, value, quantity, sigma in zip(z_vals, values, quantities, sigmas, strict=True)
    )

    return DESIDR2BAOData(
        measurements=measurements,
        covariance=covariance,
        source="CobayaSampler/bao_data desi_bao_dr2; DESI DR2 BAO arXiv:2503.14738, 2503.14739",
    )
