"""VoigtLineModel — wraps the Humlicek Voigt profile implementation.

This module preserves the original numerical algorithm from
``core.line_models.profile`` while exposing the ``LineModel`` ABC interface.
Disk integration lives in ``core.line_models.disk_integration``.
"""
from __future__ import annotations

import numpy as np

from core.line_models.base import LineModel

_C_KMS = 2.99792458e5  # speed of light in km/s


def _humlicek_voigt(wl_gauss_norm: np.ndarray,
                    width_lorentz: float) -> np.ndarray:
    """Humlicek (1982) complex error function w4, accurate to ~1e-4 everywhere.

    Parameters
    ----------
    wl_gauss_norm:
        Distance from line centre in Gaussian-width units, shape (...,).
    width_lorentz:
        Lorentzian-to-Gaussian width ratio (scalar).

    Returns
    -------
    w4 : complex ndarray, same shape as *wl_gauss_norm*.
        Real part = Voigt profile H; imag part = 2× Faraday-Voigt function.
    """
    z = width_lorentz - 1j * wl_gauss_norm
    zz = z * z
    s = np.abs(wl_gauss_norm) + width_lorentz
    w4 = np.zeros(wl_gauss_norm.shape, dtype=complex)

    con1 = s >= 15.0
    zt = z[con1]
    w4[con1] = 0.56418958355 * zt / (0.5 + zt * zt)

    con2 = (s >= 5.5) & (s < 15.0)
    zt = z[con2]
    zzt = zz[con2]
    w4[con2] = zt * (1.4104739589 + 0.56418958355 * zzt) / (
        (3.0 + zzt) * zzt + 0.7499999999)

    con3 = (width_lorentz >= 0.195 * np.abs(wl_gauss_norm) - 0.176) & (s < 5.5)
    zt = z[con3]
    w4[con3] = (16.4954955 + zt *
                (20.2093334 + zt *
                 (11.9648172 + zt * (3.77898687 + zt * 0.564223565)))) / (
                     16.4954955 + zt * (38.8236274 + zt *
                                        (39.2712051 + zt *
                                         (21.6927370 + zt *
                                          (6.69939801 + zt)))))

    con4 = w4 == 0. + 0j
    zt = z[con4]
    zzt = zz[con4]
    w4[con4] = np.exp(zzt) - zt * (36183.30536 - zzt *
                                   (3321.990492 - zzt *
                                    (1540.786893 - zzt *
                                     (219.0312964 - zzt *
                                      (35.76682780 - zzt *
                                       (1.320521697 - zzt * 0.5641900381)))))
                                   ) / (32066.59372 - zzt *
                                        (24322.84021 - zzt *
                                         (9022.227659 - zzt *
                                          (2186.181081 - zzt *
                                           (364.2190727 - zzt *
                                            (61.57036588 - zzt *
                                             (1.841438936 - zzt)))))))

    return w4


class VoigtLineModel(LineModel):
    """Voigt line profile model in the weak-field approximation.

    Encapsulates ``localProfileAndDeriv`` and ``diskIntProfAndDeriv`` logic.
    Accepts per-cell wavelength grids of shape ``(N_vel, N_cells)`` (the same
    convention as the original code) and returns profiles without disk-summing.

    Parameters
    ----------
    line_data:
        A ``lineprofileVoigt.lineData`` instance (or any object with the same
        attributes: ``wl0``, ``str``, ``widthGauss``, ``widthLorentz``, ``g``,
        ``limbDark``, ``gravDark``, ``instRes``).
    """
    def __init__(self, line_data) -> None:
        self._ld = line_data
        # Pre-compute combined Gaussian+instrument width and adjusted line strength
        ld = line_data
        if ld.instRes > 0.0:
            vel_inst = _C_KMS / ld.instRes * 0.6005612043932249
            self._width_gauss = float(
                np.sqrt(ld.widthGauss[0]**2 + vel_inst**2))
            self._width_lorentz = float(ld.widthLorentz[0] *
                                        (ld.widthGauss[0] / self._width_gauss))
            self._line_str = float(ld.str[0] *
                                   (ld.widthGauss[0] / self._width_gauss))
        else:
            self._width_gauss = float(ld.widthGauss[0])
            self._width_lorentz = float(ld.widthLorentz[0])
            self._line_str = float(ld.str[0])
        self._g_coeff = -4.6686e-12 * float(ld.wl0[0])**2 * float(ld.g[0])

    # ------------------------------------------------------------------
    # Public helpers
    # ------------------------------------------------------------------

    def local_unscaled(self,
                       wl_cells: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Compute unscaled Stokes-I and dI/dλ profiles for all cells.

        Parameters
        ----------
        wl_cells : (N_vel, N_cells)
            Per-cell wavelength grids in nm.

        Returns
        -------
        I_unscaled : (N_vel, N_cells)
        dI_dWlG   : (N_vel, N_cells)
        """
        wl0 = float(self._ld.wl0[0])
        width_gauss_wl = self._width_gauss / _C_KMS * wl0
        wl_gauss_norm = (wl0 - wl_cells) / width_gauss_wl  # (N_vel, N_cells)
        w4 = _humlicek_voigt(wl_gauss_norm, self._width_lorentz)
        dVoigt_dWl = (self._line_str / width_gauss_wl * 2.0) * (
            -wl_gauss_norm * w4.real + self._width_lorentz * w4.imag)
        I_unscaled = 1.0 - self._line_str * w4.real
        dI_dWlG = self._g_coeff * dVoigt_dWl
        return I_unscaled, dI_dWlG

    # ------------------------------------------------------------------
    # LineModel ABC implementation
    # ------------------------------------------------------------------

    def compute_profile(
        self,
        vel_grid: np.ndarray,
        B_los: np.ndarray,
        brightness: np.ndarray,
        view_angle: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Return local StokesI, StokesV per cell, shape (N_cells, N_vel).

        Note: This does NOT include limb-darkening or projected-area weighting;
        those are applied by the caller (``disk_integrate``).
        """
        wl0 = float(self._ld.wl0[0])
        # Convert velocity grid to wavelength for each cell (no Doppler shift here)
        wl = wl0 * (1.0 + vel_grid / _C_KMS)  # (N_vel,)
        wl_cells = wl[:, np.newaxis] * np.ones(
            (1, B_los.shape[0]))  # (N_vel, N_cells)
        I_unscaled, dI_dWlG = self.local_unscaled(wl_cells)
        # Scale by brightness (N_cells,) — transposed to (N_cells, N_vel)
        stokes_I = (I_unscaled *
                    brightness[np.newaxis, :]).T  # (N_cells, N_vel)
        stokes_V = (B_los[np.newaxis, :] * brightness[np.newaxis, :] *
                    dI_dWlG).T
        return stokes_I, stokes_V

    def compute_derivatives(
        self,
        vel_grid: np.ndarray,
        B_los: np.ndarray,
        dB_los_d_coeff: np.ndarray,
        brightness: np.ndarray,
        view_angle: np.ndarray,
    ) -> dict[str, np.ndarray]:
        """Return dict: dI_dBright (N_cells, N_vel), dV_dBlos (N_cells, N_vel),
        dV_dCoeff (N_coeff, N_cells, N_vel)."""
        stokes_I, stokes_V = self.compute_profile(vel_grid, B_los, brightness,
                                                  view_angle)
        wl0 = float(self._ld.wl0[0])
        wl = wl0 * (1.0 + vel_grid / _C_KMS)
        wl_cells = wl[:, np.newaxis] * np.ones((1, B_los.shape[0]))
        _, dI_dWlG = self.local_unscaled(wl_cells)  # (N_vel, N_cells)

        # dI/dBright_k: (N_cells, N_vel)
        dI_dBright = dI_dWlG.T  # just the unscaled I per cell per vel (scaled below)
        I_unscaled = (1.0 - self._line_str * _humlicek_voigt(
            (wl0 - wl_cells) /
            (self._width_gauss / _C_KMS * wl0), self._width_lorentz).real)
        dI_dBright = I_unscaled.T  # (N_cells, N_vel) — derivative of I·b wrt b_k at cell k

        # dV/dBlos_k: (N_cells, N_vel)
        dV_dBlos = (brightness[np.newaxis, :] * dI_dWlG).T  # (N_cells, N_vel)

        # dV/dCoeff: (N_coeff, N_cells, N_vel)
        # dV_dcoeff[i,k,j] = dBlos_d_coeff[i,k] * brightness[k] * dI_dWlG[j,k]
        # Use einsum: 'ik,kj->ikj' where k=cells, j=vel
        dV_dCoeff = np.einsum('ik,kj->ikj',
                              dB_los_d_coeff * brightness[np.newaxis, :],
                              dI_dWlG.T)

        return {
            "dI_dBright": dI_dBright,
            "dV_dBlos": dV_dBlos,
            "dV_dCoeff": dV_dCoeff
        }
