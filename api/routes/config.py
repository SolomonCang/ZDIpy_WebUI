"""api/routes/config.py — GET/PUT /api/config endpoints."""

import copy
import json
import sys
from pathlib import Path
from typing import Any, Dict, Optional

from fastapi import APIRouter, HTTPException

_ROOT = str(Path(__file__).resolve().parent.parent.parent)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

_APP_ROOT = Path(_ROOT).resolve()
_REPO_ROOT = _APP_ROOT.parent.parent
_DEFAULT_CONFIG = _APP_ROOT / "config.json"
_CONFIG_LIBRARY = _REPO_ROOT / "zdi_configs"

router = APIRouter(tags=["config"])


def _library_relpath(path: Path) -> str:
    return path.relative_to(_CONFIG_LIBRARY).as_posix()


def _resolve_config_path(config_path: Optional[str] = None) -> Path:
    """Resolve a frontend config path to either config.json or zdi_configs/."""
    if not config_path:
        return _DEFAULT_CONFIG

    raw = config_path.strip()
    if not raw or raw in {"default", "config.json"}:
        return _DEFAULT_CONFIG

    rel = Path(raw)
    if rel.is_absolute():
        raise HTTPException(status_code=400,
                            detail="Config path must be relative")
    parts = rel.parts
    if parts and parts[0] == "zdi_configs":
        rel = Path(*parts[1:])
    if rel.suffix.lower() != ".json":
        raise HTTPException(status_code=400,
                            detail="Config path must point to a JSON file")

    target = (_CONFIG_LIBRARY / rel).resolve()
    if not str(target).startswith(str(_CONFIG_LIBRARY.resolve())):
        raise HTTPException(status_code=400,
                            detail="Config path escapes zdi_configs")
    return target


def _load_raw(config_path: Optional[str] = None) -> Dict[str, Any]:
    """Load config.json and return as dict (raises HTTPException on error)."""
    target = _resolve_config_path(config_path)
    try:
        with open(target, "r", encoding="utf-8") as f:
            return json.load(f)
    except FileNotFoundError:
        raise HTTPException(status_code=404,
                            detail=f"Config file not found: {target}")
    except json.JSONDecodeError as exc:
        raise HTTPException(status_code=500,
                            detail=f"Config JSON malformed: {exc}")


def _to_frontend_shape(raw: Dict[str, Any]) -> Dict[str, Any]:
    cfg = copy.deepcopy(raw)

    bri = cfg.get("brightness", {})
    ev = bri.get("epoch_variation", {})
    if ev:
        if "enabled" in ev:
            bri["epoch_var_enabled"] = ev["enabled"]
        if "var_mode" in ev:
            bri["epoch_var_mode"] = ev["var_mode"]
        if "lifetime_cycles" in ev:
            bri["epoch_var_lifetime"] = ev["lifetime_cycles"]
        if "kk" in ev:
            bri["epoch_var_kk"] = ev["kk"]
        if "dphi" in ev:
            bri["epoch_var_dphi"] = ev["dphi"]
        if "dT" in ev:
            bri["epoch_var_dt"] = ev["dT"]

    return cfg


@router.get("/config/library")
def list_config_library() -> Dict[str, Any]:
    """Return all selectable ZDI config files under the root zdi_configs/."""
    _CONFIG_LIBRARY.mkdir(parents=True, exist_ok=True)
    files = []
    for path in sorted(_CONFIG_LIBRARY.rglob("*.json")):
        if not path.is_file():
            continue
        rel = _library_relpath(path)
        parts = Path(rel).parts
        if len(parts) < 5:
            continue
        target, epoch, data_type, stokes = parts[:4]
        try:
            raw = json.loads(path.read_text(encoding="utf-8"))
        except Exception:
            raw = {}
        meta = raw.get("zdi_config_meta", {}) if isinstance(raw, dict) else {}
        obs = raw.get("observations", {}) if isinstance(raw, dict) else {}
        files.append({
            "path": rel,
            "target": meta.get("target", target),
            "epoch": meta.get("epoch", epoch),
            "data_type": meta.get("data_type", data_type),
            "stokes": meta.get("stokes", stokes),
            "name": path.name,
            "n_observations": len(obs.get("files", [])),
            "updated": path.stat().st_mtime,
        })

    return {
        "root": "zdi_configs",
        "files": files,
        "targets": sorted({f["target"]
                           for f in files}),
    }


@router.get("/config")
def get_config(config_path: Optional[str] = None) -> Dict[str, Any]:
    """Return a config JSON object.

    Converts nested epoch_variation structure to flat frontend field names for
    compatibility. When config_path is provided, it must point inside the root
    zdi_configs/ library.
    """
    return _to_frontend_shape(_load_raw(config_path))


@router.put("/config")
def put_config(body: Dict[str, Any],
               config_path: Optional[str] = None) -> Dict[str, Any]:
    """Validate and persist an updated config to config.json.

    Returns the saved config on success, or raises 422 on validation error.
    """
    target = _resolve_config_path(config_path)

    # Minimal structural validation — ensure required top-level keys
    required = {
        "star", "grid", "inversion", "magnetic", "brightness", "line_model",
        "instrument", "velocity_grid", "observations", "output"
    }
    missing = required - body.keys()
    if missing:
        raise HTTPException(
            status_code=422,
            detail=f"Missing config sections: {sorted(missing)}")

    # Convert flat frontend epoch_variation fields to nested structure
    bri = body.get("brightness", {})
    if "epoch_var_enabled" in bri or "epoch_var_mode" in bri or \
         "epoch_var_lifetime" in bri or "epoch_var_kk" in bri or "epoch_var_dphi" in bri or "epoch_var_dt" in bri:
        # Ensure nested epoch_variation dict exists
        if "epoch_variation" not in bri:
            bri["epoch_variation"] = {}
        # Map flat keys to nested structure
        if "epoch_var_enabled" in bri:
            bri["epoch_variation"]["enabled"] = bri.pop("epoch_var_enabled")
        if "epoch_var_mode" in bri:
            bri["epoch_variation"]["var_mode"] = bri.pop("epoch_var_mode")
        if "epoch_var_lifetime" in bri:
            bri["epoch_variation"]["lifetime_cycles"] = bri.pop(
                "epoch_var_lifetime")
        if "epoch_var_kk" in bri:
            bri["epoch_variation"]["kk"] = bri.pop("epoch_var_kk")
        if "epoch_var_dphi" in bri:
            bri["epoch_variation"]["dphi"] = bri.pop("epoch_var_dphi")
        if "epoch_var_dt" in bri:
            bri["epoch_variation"]["dT"] = bri.pop("epoch_var_dt")
        body["brightness"] = bri

    # Attempt a lightweight ZDIConfig parse to catch value errors early
    try:
        import config_loader as cl
        if body.get("observations", {}).get("files"):
            cl.ZDIConfig.from_dict(body, str(target.parent))
    except Exception as exc:
        raise HTTPException(status_code=422,
                            detail=f"Config validation error: {exc}")

    # Persist to config.json or a selected zdi_configs file.
    try:
        target.parent.mkdir(parents=True, exist_ok=True)
        with open(target, "w", encoding="utf-8") as f:
            json.dump(body, f, indent=2, ensure_ascii=False)
    except OSError as exc:
        raise HTTPException(status_code=500,
                            detail=f"Failed to write config: {exc}")

    return _to_frontend_shape(body)
