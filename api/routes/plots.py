"""api/routes/plots.py — GET /api/plots/{profiles|surface_map|light_curve}.

Endpoints that return Plotly-ready JSON dicts for the frontend to render
directly with ``Plotly.newPlot(div, resp.data, resp.layout)``.
"""

import sys
from pathlib import Path
from typing import Any, Dict

import numpy as np
from fastapi import APIRouter, HTTPException, Query

_ROOT = str(Path(__file__).resolve().parent.parent.parent)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

router = APIRouter(tags=["plots"])


def _get_result():
    """Read ZDIResult from shared run state; raise 404 when absent."""
    from api.state import get_state  # noqa: PLC0415
    result = get_state().get("result")
    if result is None:
        raise HTTPException(
            status_code=404,
            detail="No result available. Run ZDI inversion first.",
        )
    return result


# ---------------------------------------------------------------------------
# GET /api/plots/profiles
# ---------------------------------------------------------------------------
@router.get("/plots/profiles")
def get_profiles_plot() -> Dict[str, Any]:
    """Return Stokes I/V profile comparison as a Plotly JSON dict."""
    from core.plotting.plotly_backend import PlotlyBackend  # noqa: PLC0415
    from core.plotting.data import ProfilePlotData  # noqa: PLC0415

    result = _get_result()

    obs_I, mod_I, obs_V, mod_V = [], [], [], []
    obs_I_sigma, obs_V_sigma, phases = [], [], []
    vel_grid = np.array([])

    for obs, syn in zip(result["observed_profiles"],
                        result["synthetic_profiles"]):
        phases.append(obs.get("phase", syn.get("phase", 0.0)))
        vel_grid = np.array(obs["vel"])
        obs_I.append(np.array(obs["I_obs"]))
        mod_I.append(np.array(syn["I_mod"]))
        obs_V.append(np.array(obs["V_obs"]))
        mod_V.append(np.array(syn["V_mod"]))
        obs_I_sigma.append(np.array(obs.get("I_sig", [])))
        obs_V_sigma.append(np.array(obs.get("V_sig", [])))

    plot_data = ProfilePlotData(
        phases=phases,
        vel_grid=vel_grid,
        obs_I=obs_I,
        mod_I=mod_I,
        obs_V=obs_V,
        mod_V=mod_V,
        obs_I_sigma=obs_I_sigma,
        obs_V_sigma=obs_V_sigma,
    )
    return PlotlyBackend().plot_profiles(plot_data)


# ---------------------------------------------------------------------------
# GET /api/plots/surface_map
# ---------------------------------------------------------------------------
_VALID_MAP_TYPES = {"brightness", "radial_B", "meridional_B", "azimuthal_B"}


@router.get("/plots/surface_map")
def get_surface_map_plot(
    map_type: str = Query(
        "brightness",
        description="brightness | radial_B | meridional_B | azimuthal_B"),
) -> Dict[str, Any]:
    """Return a stellar surface map as a Plotly JSON dict."""
    if map_type not in _VALID_MAP_TYPES:
        raise HTTPException(
            status_code=422,
            detail=f"map_type must be one of {sorted(_VALID_MAP_TYPES)}",
        )

    from core.plotting.plotly_backend import PlotlyBackend  # noqa: PLC0415
    from core.plotting.data import SurfaceMapData  # noqa: PLC0415

    result = _get_result()

    _meta = result.get("metadata", {})
    clat = np.array(_meta.get("clat", []))
    lon = np.array(_meta.get("lon", []))
    if clat.size == 0 or lon.size == 0:
        raise HTTPException(
            status_code=503,
            detail=
            "Grid coordinates (clat/lon) not available in result metadata. "
            "Re-run the inversion with the updated pipeline.",
        )

    _mag_coeffs = result.get("mag_coeffs", {})

    def _to_complex(lst):
        """Reconstruct complex array from [[real, imag], ...] pairs."""
        return np.array([r + 1j * i for r, i in lst], dtype=complex)

    if map_type == "brightness":
        values = np.array(result["bright_map"])
        vmin, vmax = float(np.min(values)), float(np.max(values))
    else:
        # Magnetic field components require recomputing from coefficients + grid
        try:
            import core.magneticGeom as mG  # noqa: PLC0415
            mag_geom = mG.magSphHarmonics(len(_mag_coeffs.get("alpha", [])))
            mag_geom.alpha = _to_complex(_mag_coeffs["alpha"])
            mag_geom.beta = _to_complex(_mag_coeffs["beta"])
            mag_geom.gamma = _to_complex(_mag_coeffs["gamma"])
            mag_geom.initMagGeom(clat, lon)
            vec_b = mag_geom.getAllMagVectors()  # (3, N_cells)
            component_map = {
                "radial_B": 0,
                "meridional_B": 1,
                "azimuthal_B": 2
            }
            values = vec_b[component_map[map_type], :]
            abs_max = float(np.max(np.abs(values)))
            vmin, vmax = -abs_max, abs_max
        except Exception as exc:
            raise HTTPException(
                status_code=500,
                detail=f"Failed to compute magnetic field map: {exc}",
            ) from exc

    plot_data = SurfaceMapData(
        clat=clat,
        lon=lon,
        values=np.asarray(values),
        map_type=map_type,
        vmin=vmin,
        vmax=vmax,
    )
    return PlotlyBackend().plot_surface_map(plot_data)


# ---------------------------------------------------------------------------
# GET /api/plots/light_curve
# ---------------------------------------------------------------------------
@router.get("/plots/light_curve")
def get_light_curve_plot() -> Dict[str, Any]:
    """Return light curve comparison as a Plotly JSON dict.

    Returns HTTP 404 when no light curve data is present in the result.
    """
    from core.plotting.plotly_backend import PlotlyBackend  # noqa: PLC0415
    from core.plotting.data import LightCurvePlotData  # noqa: PLC0415

    result = _get_result()
    if result.get("light_curve_synthetic") is None:
        raise HTTPException(
            status_code=404,
            detail=
            "No light curve data in result (fit_light_curve=0 or not yet fitted).",
        )

    _meta = result.get("metadata", {})
    lc_obs = _meta.get("lc_obs", {})
    plot_data = LightCurvePlotData(
        jdates=np.array(lc_obs.get("jdates", [])),
        obs_flux=np.array(lc_obs.get("flux", [])),
        mod_flux=np.array(result["light_curve_synthetic"]),
        sigma=np.array(lc_obs.get("sigma", [])),
    )
    return PlotlyBackend().plot_light_curve(plot_data)
