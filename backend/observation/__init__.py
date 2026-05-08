"""Observational data layer — Stage 2.

Loaders for public cosmological datasets:

- DESI DR2 BAO measurements (H(z), fσ₈(z)) with full covariance matrix
- Planck 2018 derived parameters (validation only, not for fitting)

Not yet implemented — Stage 2 deliverable.
"""

from observation.desi_dr2 import DESIDR2BAOData, DESIBAOMeasurement, load_desi_dr2_bao

__all__ = [
    "DESIDR2BAOData",
    "DESIBAOMeasurement",
    "load_desi_dr2_bao",
]
