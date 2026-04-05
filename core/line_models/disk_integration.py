"""Standalone disk integration for spectral line profiles.

Decoupled from any specific line model so it can be reused with any
``LineModel`` implementation.
"""
from __future__ import annotations

import numpy as np
from scipy.ndimage import gaussian_filter1d


def disk_integrate(
        profiles: np.ndarray,  # (N_phases, N_cells, N_vel) or (N_cells, N_vel)
        proj_areas: np.ndarray,  # (N_phases, N_cells) or (N_cells,)
        limb_weights: np.ndarray,  # same shape as proj_areas
        inst_fwhm_pix:
    float = 0.0,  # instrumental FWHM in pixels (0 = no convolution)
) -> np.ndarray:
    """Sum over surface elements, optionally convolve with instrumental profile.

    Parameters
    ----------
    profiles : (N_phases, N_cells, N_vel) or (N_cells, N_vel)
        Local Stokes I or V profiles, already brightness-scaled.
    proj_areas : (N_phases, N_cells) or (N_cells,)
        Projected area weights (includes visibility mask).
    limb_weights : same shape as *proj_areas*
        Limb-darkening (and gravity-darkening) multipliers.
    inst_fwhm_pix : float
        Instrumental FWHM expressed in pixels along the velocity axis.
        Convolution is performed with a Gaussian (sigma = fwhm / 2.355).
        Skipped when <= 0.

    Returns
    -------
    ndarray, shape (N_phases, N_vel) or (N_vel,)
    """
    squeeze = profiles.ndim == 2
    if squeeze:
        profiles = profiles[np.newaxis]  # (1, N_cells, N_vel)
        proj_areas = proj_areas[np.newaxis]
        limb_weights = limb_weights[np.newaxis]

    weights = proj_areas * limb_weights  # (N_phases, N_cells)
    # Weighted sum over cells: (N_phases, N_vel)
    integrated = np.einsum('pc,pcv->pv', weights, profiles)

    if inst_fwhm_pix > 0.0:
        sigma = inst_fwhm_pix / (2.0 * np.sqrt(2.0 * np.log(2.0)))
        integrated = gaussian_filter1d(integrated, sigma=sigma, axis=-1)

    return integrated[0] if squeeze else integrated


def normalize_by_continuum(
        stokes_I: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (I/Ic, Ic) where Ic is the per-phase continuum (sum of weights).

    This matches the normalisation convention used in the original
    ``localProfileAndDeriv.updateProfDeriv``.

    Parameters
    ----------
    stokes_I : (N_phases, N_vel)

    Returns
    -------
    I_norm : (N_phases, N_vel) — I/Ic evaluated at continuum level
    Ic     : (N_phases,)
    """
    # Continuum is approximated as the mean of the first and last 5 velocity points
    n = stokes_I.shape[-1]
    edge = max(1, n // 10)
    Ic = 0.5 * (stokes_I[..., :edge].mean(axis=-1) +
                stokes_I[..., -edge:].mean(axis=-1))
    return stokes_I / Ic[..., np.newaxis], Ic
