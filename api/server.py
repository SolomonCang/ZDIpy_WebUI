"""api/server.py — FastAPI application factory for ZDIpy WebUI."""

import os
import sys
from pathlib import Path

from fastapi import FastAPI
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles

# Ensure project root is importable from every context
_ROOT = str(Path(__file__).resolve().parent.parent)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from api.routes import config as cfg_routes  # noqa: E402
from api.routes import observations as obs_routes  # noqa: E402
from api.routes import results as res_routes  # noqa: E402
from api.routes import run as run_routes  # noqa: E402

_FRONTEND_DIR = str(Path(_ROOT) / "frontend")

# ---------------------------------------------------------------------------
# Application
# ---------------------------------------------------------------------------
app = FastAPI(
    title="ZDIpy WebUI API",
    description=
    ("REST + SSE backend for the ZDIpy Zeeman Doppler Imaging web interface. "
     "Interactive API docs available at /docs."),
    version="2.0.0",
)

# --- API routes (must be registered BEFORE the static-file fallback) --------
app.include_router(cfg_routes.router, prefix="/api")
app.include_router(obs_routes.router, prefix="/api")
app.include_router(run_routes.router, prefix="/api")
app.include_router(res_routes.router, prefix="/api")


# --- Root: serve index.html explicitly so /docs still works via FastAPI -----
@app.get("/", include_in_schema=False)
async def root() -> FileResponse:
    return FileResponse(os.path.join(_FRONTEND_DIR, "index.html"))


# --- Static assets: CSS / JS (registered last to avoid shadowing /api) ------
os.makedirs(os.path.join(_FRONTEND_DIR, "css"), exist_ok=True)
os.makedirs(os.path.join(_FRONTEND_DIR, "js"), exist_ok=True)

app.mount("/css",
          StaticFiles(directory=os.path.join(_FRONTEND_DIR, "css")),
          name="css")
app.mount("/js",
          StaticFiles(directory=os.path.join(_FRONTEND_DIR, "js")),
          name="js")
