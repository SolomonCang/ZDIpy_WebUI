"""api/routes/observations.py — Manage LSD profile files in LSDprof/."""

import os
import re
import sys
from pathlib import Path
from typing import Dict, List

from fastapi import APIRouter, File, HTTPException, UploadFile
from pydantic import BaseModel

_ROOT = str(Path(__file__).resolve().parent.parent.parent)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

_LSD_DIR = Path(_ROOT) / "LSDprof"
_ALLOWED_EXTS = {".prof", ".norm", ".lsd", ".txt", ".dat"}
_SAFE_FNAME = re.compile(r"^[\w\-. ]+$")  # alphanums + _ - . space only

router = APIRouter(tags=["observations"])


def _validate_fname(raw: str) -> str:
    """Return sanitised basename or raise HTTPException."""
    basename = os.path.basename(raw)
    if not basename or "\x00" in basename:
        raise HTTPException(status_code=400, detail="Invalid filename")
    if not _SAFE_FNAME.match(basename):
        raise HTTPException(status_code=400,
                            detail=f"Unsafe filename: {basename!r}")
    ext = Path(basename).suffix.lower()
    if ext not in _ALLOWED_EXTS:
        raise HTTPException(
            status_code=400,
            detail=f"Extension {ext!r} not allowed. Allowed: {sorted(_ALLOWED_EXTS)}",
        )
    # Path traversal guard
    dest = (_LSD_DIR / basename).resolve()
    if not str(dest).startswith(str(_LSD_DIR.resolve())):
        raise HTTPException(status_code=400, detail="Path traversal detected")
    return basename


@router.get("/observations")
def list_observations() -> List[dict]:
    """List all LSD profile files in LSDprof/."""
    _LSD_DIR.mkdir(parents=True, exist_ok=True)
    files = []
    for p in sorted(_LSD_DIR.iterdir()):
        if p.is_file() and p.suffix.lower() in _ALLOWED_EXTS:
            files.append({"name": p.name, "size_bytes": p.stat().st_size})
    return files


@router.post("/observations/upload")
async def upload_observations(files: List[UploadFile] = File(...)) -> dict:
    """Upload one or more LSD profile files to LSDprof/."""
    _LSD_DIR.mkdir(parents=True, exist_ok=True)
    saved: List[str] = []
    for uf in files:
        basename = _validate_fname(uf.filename or "")
        dest = _LSD_DIR / basename
        content = await uf.read()
        dest.write_bytes(content)
        saved.append(basename)
    return {"saved": saved, "count": len(saved)}


@router.delete("/observations/{fname}")
def delete_observation(fname: str) -> dict:
    """Delete a specific LSD profile file from LSDprof/."""
    basename = _validate_fname(fname)
    target = _LSD_DIR / basename
    if not target.exists():
        raise HTTPException(status_code=404, detail=f"{basename!r} not found")
    target.unlink()
    return {"deleted": basename}


class ValidatePathsRequest(BaseModel):
    paths: List[str]


class PathStatus(BaseModel):
    exists: bool
    stokes: str  # "I" | "I+V+N" | "unknown"


@router.post("/observations/validate")
def validate_observation_paths(
        body: ValidatePathsRequest) -> Dict[str, PathStatus]:
    """Check which observation file paths exist and detect Stokes columns."""
    root = Path(_ROOT).resolve()
    results: Dict[str, PathStatus] = {}
    for p in body.paths:
        if not p:
            results[p] = PathStatus(exists=False, stokes="unknown")
            continue
        try:
            target = (root / p).resolve()
        except Exception:
            results[p] = PathStatus(exists=False, stokes="unknown")
            continue
        # Path traversal guard
        if not str(target).startswith(str(root)):
            results[p] = PathStatus(exists=False, stokes="unknown")
            continue
        if not target.is_file():
            results[p] = PathStatus(exists=False, stokes="unknown")
            continue
        # Detect number of columns from second header line (Donati LSD format)
        stokes = "unknown"
        try:
            with open(target, "r", encoding="utf-8", errors="replace") as fh:
                fh.readline()
                second = fh.readline()
            num_cols = int(second.split()[1])
            stokes = "I" if num_cols == 2 else "I+V+N"
        except Exception:
            pass
        results[p] = PathStatus(exists=True, stokes=stokes)
    return results
