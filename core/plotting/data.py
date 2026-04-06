"""core/plotting/data.py — Plot data container dataclasses."""

from dataclasses import dataclass, field
import numpy as np


@dataclass
class ProfilePlotData:
    """Stokes I/V line profile observation vs. model comparison data."""
    phases: list[float]  # rotation phase for each observation
    vel_grid: np.ndarray  # (N_vel,) km/s
    obs_I: list[np.ndarray]  # len=N_phases, each (N_vel,)
    mod_I: list[np.ndarray]
    obs_V: list[np.ndarray]
    mod_V: list[np.ndarray]
    obs_I_sigma: list[np.ndarray] = field(default_factory=list)
    obs_V_sigma: list[np.ndarray] = field(default_factory=list)


@dataclass
class SurfaceMapData:
    """Stellar surface map (equal-area projection) data. map_type determines the physical quantity."""
    clat: np.ndarray  # (N_cells,) colatitude (rad)
    lon: np.ndarray  # (N_cells,) longitude (rad)
    values: np.ndarray  # (N_cells,) brightness or magnetic field component
    map_type: str  # "brightness" | "radial_B" | "meridional_B" | "azimuthal_B"
    vmin: float | None = None
    vmax: float | None = None


@dataclass
class LightCurvePlotData:
    """Light curve observation vs. model comparison data."""
    jdates: np.ndarray  # (N_obs,)
    obs_flux: np.ndarray  # (N_obs,)
    mod_flux: np.ndarray  # (N_obs,)
    sigma: np.ndarray  # (N_obs,)


@dataclass
class MagneticPolarData:
    """Three-component polar-projection magnetic field map data.

    The 2-D arrays (Br, Blon, Blat) are on a regular (npClat × npLon) grid
    that spans latitudes from ~-30° to ~+90° using a polar-coordinate view.
    """
    Br: np.ndarray  # (npClat, npLon) radial component [G]
    Blon: np.ndarray  # (npClat, npLon) azimuthal component [G]
    Blat: np.ndarray  # (npClat, npLon) meridional component [G]
    lon_grid: np.ndarray  # (npLon,) longitude [rad], 0..2π
    clat_grid: np.ndarray  # (npClat,) colatitude [rad], eps..π-eps
    obs_phases: np.ndarray  # (N_obs,) rotation phases in [0, 1)
    discrete_levels: int = 10
