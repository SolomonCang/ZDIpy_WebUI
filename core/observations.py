"""core/observations.py — Observation data container dataclasses."""

from dataclasses import dataclass
import numpy as np


@dataclass
class SpectralObservation:
    """Data container for a single LSD line profile observation file."""
    filename: str
    jdate: float
    vel_center_kms: float
    wl: np.ndarray  # (N_vel,)  velocity grid (km/s)
    stokes_I: np.ndarray  # (N_vel,)
    stokes_I_sigma: np.ndarray
    stokes_V: np.ndarray  # (N_vel,)
    stokes_V_sigma: np.ndarray


@dataclass
class LightCurveObservation:
    """Data container for a single photometric light curve observation point."""
    jdate: float
    flux: float
    flux_sigma: float
    band: str = "V"
