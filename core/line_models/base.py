"""Abstract base class for spectral line physical models."""
from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np


class LineModel(ABC):
    """Spectral line physical model interface.

    Implementing classes provide local line profiles and their derivatives
    with respect to model parameters, evaluated at each surface element.
    """
    @abstractmethod
    def compute_profile(
        self,
        vel_grid: np.ndarray,  # (N_vel,)  velocity grid in km/s
        B_los: np.ndarray,  # (N_cells,) line-of-sight B field (kG)
        brightness: np.ndarray,  # (N_cells,) brightness scale
        view_angle: np.ndarray,  # (N_cells,) angle to line-of-sight (rad)
    ) -> tuple[np.ndarray, np.ndarray]:
        """Return (StokesI, StokesV) local profiles, each shape (N_cells, N_vel)."""

    @abstractmethod
    def compute_derivatives(
            self,
            vel_grid: np.ndarray,  # (N_vel,)
            B_los: np.ndarray,  # (N_cells,)
            dB_los_d_coeff: np.ndarray,  # (N_coeff, N_cells)
            brightness: np.ndarray,  # (N_cells,)
            view_angle: np.ndarray,  # (N_cells,)
    ) -> dict[str, np.ndarray]:
        """Return dict with keys: dI_dBright, dV_dBlos, dV_dCoeff."""
