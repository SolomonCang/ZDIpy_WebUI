"""api/routes/results.py — GET /api/results/{profiles|magnetic|brightness}."""

import json
import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
from fastapi import APIRouter

_ROOT = str(Path(__file__).resolve().parent.parent.parent)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

router = APIRouter(tags=["results"])


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _load_dat(path: str) -> Optional[np.ndarray]:
    """Load whitespace-delimited data, skip comment/blank lines."""
    if not os.path.isfile(path):
        return None
    rows: List[List[float]] = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                try:
                    rows.append([float(x) for x in line.split()])
                except ValueError:
                    pass
    return np.array(rows) if rows else None


def _parse_phases(path: str) -> List[np.ndarray]:
    """Parse a multi-phase profile file (phases separated by blank/comment lines)."""
    phases: List[np.ndarray] = []
    current: List[List[float]] = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                if current:
                    phases.append(np.array(current))
                    current = []
            else:
                try:
                    current.append([float(x) for x in line.split()])
                except ValueError:
                    pass
    if current:
        phases.append(np.array(current))
    return phases


def _col(arr: np.ndarray, idx: int) -> Optional[List[float]]:
    """Return column idx as list, or None if the column doesn't exist."""
    if arr.shape[1] > idx:
        return arr[:, idx].tolist()
    return None


# ---------------------------------------------------------------------------
# GET /api/results/profiles
# ---------------------------------------------------------------------------
@router.get("/results/profiles")
def get_profiles() -> Dict[str, Any]:
    """Parse outObserved.dat + outLineModels.dat and return phase data for Plotly."""
    obs_file = os.path.join(_ROOT, "outObserved.dat")
    mod_file = os.path.join(_ROOT, "outLineModels.dat")

    if not os.path.isfile(obs_file) or not os.path.isfile(mod_file):
        return {"available": False, "phases": []}

    obs_phases = _parse_phases(obs_file)
    mod_phases = _parse_phases(mod_file)
    n = min(len(obs_phases), len(mod_phases))
    if n == 0:
        return {"available": False, "phases": []}

    phases = []
    for i in range(n):
        obs = obs_phases[i]
        mod = mod_phases[i]
        phase: Dict[str, Any] = {"index": i}
        phase["vel_obs"] = _col(obs, 0) or []
        phase["vel_mod"] = _col(mod, 0) or []
        phase["stokes_I_obs"] = _col(obs, 1)
        phase["stokes_I_err"] = _col(obs, 2)
        phase["stokes_V_obs"] = _col(obs, 3)
        phase["stokes_V_err"] = _col(obs, 4)
        phase["stokes_I_mod"] = _col(mod, 1)
        phase["stokes_V_mod"] = _col(mod, 2)
        phases.append(phase)

    return {"available": True, "phases": phases}


# ---------------------------------------------------------------------------
# GET /api/results/magnetic
# ---------------------------------------------------------------------------
@router.get("/results/magnetic")
def get_magnetic() -> Dict[str, Any]:
    """Reconstruct Br/Bclat/Blon surface maps from outMagCoeff.dat."""
    mag_file = os.path.join(_ROOT, "outMagCoeff.dat")
    if not os.path.isfile(mag_file):
        return {"available": False}

    try:
        # Load config for lMax and nRings
        cfg_path = os.path.join(_ROOT, "config", "_run_config.json")
        if not os.path.isfile(cfg_path):
            cfg_path = os.path.join(_ROOT, "config", "config.json")
        with open(cfg_path) as f:
            cfg = json.load(f)
        l_max = int(cfg.get("magnetic", {}).get("l_max", 15))
        n_rings = int(cfg.get("grid", {}).get("nRings", 30))

        import core.geometryStellar as gS  # noqa: PLC0415
        import core.magneticGeom as mG  # noqa: PLC0415

        s_grid = gS.starGrid(n_rings, verbose=0)
        mag_geom = mG.SetupMagSphHarmoics(s_grid,
                                          1,
                                          mag_file,
                                          l_max,
                                          "full",
                                          verbose=0)
        vec_b = mag_geom.getAllMagVectors()  # shape (3, nPoints)

        clat_deg = np.degrees(s_grid.clat).tolist()
        lon_deg = np.degrees(s_grid.long).tolist()
        lat_deg = (90.0 - np.degrees(s_grid.clat)).tolist()

        return {
            "available": True,
            "lat": lat_deg,
            "lon": lon_deg,
            "clat": clat_deg,
            "Br": vec_b[0, :].tolist(),
            "Bclat": vec_b[1, :].tolist(),
            "Blon": vec_b[2, :].tolist(),
        }
    except Exception as exc:
        return {"available": False, "error": str(exc)}


# ---------------------------------------------------------------------------
# GET /api/results/brightness
# ---------------------------------------------------------------------------
@router.get("/results/brightness")
def get_brightness() -> Dict[str, Any]:
    """Parse outBrightMap.dat (columns: clat_rad, lon_rad, brightness)."""
    bri_file = os.path.join(_ROOT, "outBrightMap.dat")
    data = _load_dat(bri_file)
    if data is None or data.shape[1] < 3:
        return {"available": False}

    clat_deg = np.degrees(data[:, 0]).tolist()
    lon_deg = np.degrees(data[:, 1]).tolist()
    lat_deg = (90.0 - np.degrees(data[:, 0])).tolist()
    brightness = data[:, 2].tolist()

    return {
        "available": True,
        "lat": lat_deg,
        "lon": lon_deg,
        "clat": clat_deg,
        "brightness": brightness,
    }
