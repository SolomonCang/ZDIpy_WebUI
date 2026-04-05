"""Synthetic photometry (light curve) calculation for ZDI inversion.

Provides ``compute_synthetic_photometry`` which computes the disk-integrated
flux for each rotation phase given a brightness map and observing geometry.
The function signature and normalisation follow Section 7 of the architecture
specification.
"""
from __future__ import annotations

import numpy as np
from core.line_models.line_utils import limbDarkening  # noqa: F401


def compute_synthetic_photometry(
    star_grid,
    visible_grid,
    bright_map,
    limb_dark_coeff: float,
    grav_dark_coeff: float,
) -> float | np.ndarray:
    """Compute normalised photometric flux for one or more rotation phases.

    The flux is normalised so that a perfectly uniform, fully bright star
    (``bright_map.bright = 1`` everywhere) returns 1.0 per phase.

    Parameters
    ----------
    star_grid :
        ``core.geometryStellar.starGrid`` instance.
    visible_grid :
        Either a ``_SinglePhaseView`` (single phase) or ``BatchVisibleGrid``
        (N_phases).  The function detects the type from attribute shapes.
    bright_map :
        ``core.brightnessGeom.brightMap`` instance.  ``bright`` is shape
        ``(N_cells,)``.
    limb_dark_coeff : float
        Linear limb-darkening coefficient η.
    grav_dark_coeff : float
        Gravity-darkening coefficient (0 for spherical stars).

    Returns
    -------
    float or (N_phases,) ndarray
        Normalised photometric flux.

    Notes
    -----
    The light-curve Jacobian ``dFlux/dBright_i`` is returned by
    ``compute_lc_jacobian`` below.
    """
    bright = bright_map.bright  # (N_cells,)
    is_batch = hasattr(visible_grid,
                       'proj_area') and visible_grid.proj_area.ndim == 2

    if is_batch:
        # visible_grid.proj_area: (N_phases, N_cells)
        # visible_grid.visible: (N_phases, N_cells)
        # visible_grid.view_angle: (N_phases, N_cells)
        proj_area = visible_grid.proj_area  # (N_phases, N_cells)
        visible = visible_grid.visible.astype(float)  # (N_phases, N_cells)
        view_angle = visible_grid.view_angle  # (N_phases, N_cells)

        limb = limbDarkening(limb_dark_coeff,
                             view_angle)  # (N_phases, N_cells)
        grav = star_grid.gravityDarkening(grav_dark_coeff)  # (N_cells,)
        weights = proj_area * limb * grav[
            np.newaxis, :] * visible  # (N_phases, N_cells)

        flux = np.einsum('pc,c->p', weights, bright)  # (N_phases,)

        # Normalise by uniform-bright reference
        ref_flux = np.sum(weights, axis=1)  # (N_phases,)
        # Guard against zero (non-visible phases)
        ref_flux = np.where(ref_flux > 0.0, ref_flux, 1.0)
        return flux / ref_flux

    else:
        # Single-phase _SinglePhaseView
        proj_area = visible_grid.projArea  # (N_cells,)
        visible = visible_grid.visible.astype(float)
        view_angle = visible_grid.viewAngle  # (N_cells,)

        limb = limbDarkening(limb_dark_coeff, view_angle)  # (N_cells,)
        grav = visible_grid.gravityDarkening(grav_dark_coeff)  # (N_cells,)
        weights = proj_area * limb * grav * visible  # (N_cells,)

        flux = np.dot(weights, bright)
        ref_flux = np.sum(weights)
        return (flux / ref_flux) if ref_flux > 0.0 else 0.0


def compute_lc_jacobian(
    star_grid,
    visible_grid,
    bright_map,
    limb_dark_coeff: float,
    grav_dark_coeff: float,
) -> np.ndarray:
    """Compute dFlux/dBright_i for each brightness surface element.

    Used to assemble the Jacobian rows for the light-curve data constraint
    in the MEM response matrix.

    Parameters
    ----------
    (same as ``compute_synthetic_photometry``)

    Returns
    -------
    dFlux_dBright : (N_phases, N_cells) or (N_cells,) ndarray
        Partial derivative of normalised flux with respect to each cell
        brightness value.
    """
    bright = bright_map.bright  # (N_cells,)
    is_batch = hasattr(visible_grid,
                       'proj_area') and visible_grid.proj_area.ndim == 2

    if is_batch:
        proj_area = visible_grid.proj_area
        visible = visible_grid.visible.astype(float)
        view_angle = visible_grid.view_angle

        limb = limbDarkening(limb_dark_coeff, view_angle)
        grav = star_grid.gravityDarkening(grav_dark_coeff)
        weights = proj_area * limb * grav[
            np.newaxis, :] * visible  # (N_phases, N_cells)

        flux = np.einsum('pc,c->p', weights, bright)  # (N_phases,)
        ref_flux = np.einsum('pc->p', weights)  # (N_phases,)
        ref_flux = np.where(ref_flux > 0.0, ref_flux, 1.0)

        # d(flux/ref)/dBright_i = weights_i / ref
        # (assuming ref does not depend on bright when bright=const reference)
        jac = weights / ref_flux[:, np.newaxis]  # (N_phases, N_cells)
        return jac

    else:
        proj_area = visible_grid.projArea
        visible = visible_grid.visible.astype(float)
        view_angle = visible_grid.viewAngle

        limb = limbDarkening(limb_dark_coeff, view_angle)
        grav = visible_grid.gravityDarkening(grav_dark_coeff)
        weights = proj_area * limb * grav * visible

        ref_flux = np.sum(weights)
        return weights / ref_flux if ref_flux > 0.0 else weights
