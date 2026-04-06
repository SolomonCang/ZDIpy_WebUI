"""api/routes/config.py — GET/PUT /api/config endpoints."""

import json
import sys
from pathlib import Path
from typing import Any, Dict

from fastapi import APIRouter, HTTPException

_ROOT = str(Path(__file__).resolve().parent.parent.parent)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

_DEFAULT_CONFIG = Path(_ROOT) / "config.json"

router = APIRouter(tags=["config"])


def _load_raw() -> Dict[str, Any]:
    """Load config.json and return as dict (raises HTTPException on error)."""
    try:
        with open(_DEFAULT_CONFIG, "r", encoding="utf-8") as f:
            return json.load(f)
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail="config.json not found")
    except json.JSONDecodeError as exc:
        raise HTTPException(status_code=500,
                            detail=f"Config JSON malformed: {exc}")


@router.get("/config")
def get_config() -> Dict[str, Any]:
    """Return the current config.json as a JSON object."""
    return _load_raw()


@router.put("/config")
def put_config(body: Dict[str, Any]) -> Dict[str, Any]:
    """Validate and persist an updated config to config.json.

    Returns the saved config on success, or raises 422 on validation error.
    """
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

    # Attempt a lightweight ZDIConfig parse to catch value errors early
    try:
        import config_loader as cl
        if body.get("observations", {}).get("files"):
            cl.ZDIConfig.from_dict(body, config_path=str(_DEFAULT_CONFIG))
    except Exception as exc:
        raise HTTPException(status_code=422,
                            detail=f"Config validation error: {exc}")

    # Persist to config.json
    try:
        with open(_DEFAULT_CONFIG, "w", encoding="utf-8") as f:
            json.dump(body, f, indent=2, ensure_ascii=False)
    except OSError as exc:
        raise HTTPException(status_code=500,
                            detail=f"Failed to write config: {exc}")

    return body
