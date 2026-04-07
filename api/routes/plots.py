"""api/routes/plots.py — GET /api/plots/{profiles|surface_map|light_curve}.

Endpoints that return Plotly-ready JSON dicts for the frontend to render
directly with ``Plotly.newPlot(div, resp.data, resp.layout)``.
"""

import io
import math
import sys
from pathlib import Path
from typing import Any, Dict

import numpy as np
from fastapi import APIRouter, HTTPException, Query
from fastapi.responses import Response

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
            # nTot = nl*(nl+1)//2 + nl; invert to get nl from nTot
            _nTot = len(_mag_coeffs.get("alpha", []))
            _nl = int(round((-3.0 + math.sqrt(9.0 + 8.0 * _nTot)) / 2.0))
            mag_geom = mG.magSphHarmonics(_nl)
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


# ---------------------------------------------------------------------------
# GET /api/plots/magnetic_polar
# ---------------------------------------------------------------------------
@router.get("/plots/magnetic_polar")
def get_magnetic_polar_plot():
    """Render a 3-panel polar magnetic field map via Matplotlib and return PNG.

    The image is also saved to ``results/magnetic_polar.png`` in the project
    root so the user can retrieve it later.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from core.plotting.matplotlib_backend import MatplotlibBackend  # noqa: PLC0415
    from core.plotting.data import MagneticPolarData  # noqa: PLC0415
    import core.magneticGeom as mG  # noqa: PLC0415

    result = _get_result()

    _mag_coeffs = result.get("mag_coeffs", {})
    if not _mag_coeffs.get("alpha"):
        raise HTTPException(
            status_code=404,
            detail=
            "No magnetic coefficients in result. Run ZDI inversion first.",
        )

    def _to_complex(lst):
        return np.array([r + 1j * i for r, i in lst], dtype=complex)

    # ---- build dense grid (reference: 90 × 180) -------------------------
    _EPS = 1e-6
    npClat, npLon = 90, 180
    lon_grid = np.linspace(0.0, 2.0 * np.pi, npLon, endpoint=False)
    clat_grid = np.linspace(_EPS, np.pi - _EPS, npClat)
    full_clat = np.repeat(clat_grid, npLon)
    full_lon = np.tile(lon_grid, npClat)

    # ---- reconstruct magnetic field on dense grid -----------------------
    _nTot = len(_mag_coeffs["alpha"])
    _nl = int(round((-3.0 + math.sqrt(9.0 + 8.0 * _nTot)) / 2.0))
    try:
        mag_geom = mG.magSphHarmonics(_nl)
        mag_geom.alpha = _to_complex(_mag_coeffs["alpha"])
        mag_geom.beta = _to_complex(_mag_coeffs["beta"])
        mag_geom.gamma = _to_complex(_mag_coeffs["gamma"])
        mag_geom.initMagGeom(full_clat, full_lon)
        vec_b = mag_geom.getAllMagVectors()  # (3, npClat * npLon)
    except Exception as exc:
        raise HTTPException(
            status_code=500,
            detail=f"Failed to compute magnetic field: {exc}",
        ) from exc

    Br = vec_b[0].reshape(npClat, npLon)
    Blat = vec_b[1].reshape(npClat, npLon)
    Blon = vec_b[2].reshape(npClat, npLon)

    # ---- observed rotation phases (cycle % 1) ---------------------------
    syn_profs = result.get("synthetic_profiles", [])
    obs_phases = (np.array([p["phase"] % 1.0 for p in syn_profs])
                  if syn_profs else np.array([]))

    # ---- build data container and render --------------------------------
    data = MagneticPolarData(
        Br=Br,
        Blon=Blon,
        Blat=Blat,
        lon_grid=lon_grid,
        clat_grid=clat_grid,
        obs_phases=obs_phases,
    )
    fig = MatplotlibBackend().plot_magnetic_polar(data)

    # ---- save to results/ -----------------------------------------------
    results_dir = Path(_ROOT) / "results"
    results_dir.mkdir(exist_ok=True)
    save_path = results_dir / "magnetic_polar.png"
    fig.savefig(str(save_path), dpi=150, bbox_inches='tight')

    # ---- return as PNG response -----------------------------------------
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    buf.seek(0)
    return Response(content=buf.read(), media_type="image/png")


# ---------------------------------------------------------------------------
# GET /api/plots/brightness_polar
# ---------------------------------------------------------------------------
@router.get("/plots/brightness_polar")
def get_brightness_polar_plot():
    """Render a polar projection brightness map via Matplotlib and return PNG.

    The image is also saved to ``results/brightness_polar.png``.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from core.plotting.matplotlib_backend import MatplotlibBackend  # noqa: PLC0415
    from core.plotting.data import BrightnessPolarData  # noqa: PLC0415

    result = _get_result()

    bright_map = result.get("bright_map")
    if bright_map is None:
        raise HTTPException(
            status_code=404,
            detail="No brightness map in result. Run ZDI inversion first.",
        )

    # ---- build dense grid (90 × 180) ------------------------------------
    _EPS = 1e-6
    npClat, npLon = 90, 180
    lon_grid = np.linspace(0.0, 2.0 * np.pi, npLon, endpoint=False)
    clat_grid = np.linspace(_EPS, np.pi - _EPS, npClat)

    # ---- interpolate / reshape bright_map onto dense grid ----------------
    # The bright_map is stored on the stellar grid (nRings*2 cells).
    # Retrieve the original clat/lon from result metadata to interpolate.
    _meta = result.get("metadata", {})
    orig_clat = np.array(_meta.get("clat", []))
    orig_lon = np.array(_meta.get("lon", []))
    orig_bri = np.array(bright_map)

    if orig_clat.size == 0 or orig_lon.size == 0:
        raise HTTPException(
            status_code=503,
            detail="Grid coordinates (clat/lon) not in result metadata. "
            "Re-run the inversion with the updated pipeline.",
        )

    # Nearest-neighbour interpolation onto dense grid
    from scipy.interpolate import griddata  # noqa: PLC0415
    orig_lat = np.pi / 2.0 - orig_clat
    orig_lon_plot = orig_lon

    full_clat = np.repeat(clat_grid, npLon)
    full_lon = np.tile(lon_grid, npClat)
    dense_lat = np.pi / 2.0 - full_clat
    dense_bri = griddata(
        np.column_stack([orig_lon_plot, orig_lat]),
        orig_bri,
        np.column_stack([full_lon, dense_lat]),
        method='nearest',
    )
    brightness_2d = dense_bri.reshape(npClat, npLon)

    # ---- observed rotation phases (cycle % 1) ---------------------------
    syn_profs = result.get("synthetic_profiles", [])
    obs_phases = (np.array([p["phase"] % 1.0 for p in syn_profs])
                  if syn_profs else np.array([]))

    # ---- render ---------------------------------------------------------
    data = BrightnessPolarData(
        brightness=brightness_2d,
        lon_grid=lon_grid,
        clat_grid=clat_grid,
        obs_phases=obs_phases,
    )
    fig = MatplotlibBackend().plot_brightness_polar(data)

    # ---- save to results/ -----------------------------------------------
    results_dir = Path(_ROOT) / "results"
    results_dir.mkdir(exist_ok=True)
    save_path = results_dir / "brightness_polar.png"
    fig.savefig(str(save_path), dpi=150, bbox_inches='tight')

    # ---- return as PNG --------------------------------------------------
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    buf.seek(0)
    return Response(content=buf.read(), media_type="image/png")
