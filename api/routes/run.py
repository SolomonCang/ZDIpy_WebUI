"""api/routes/run.py — POST /api/run  |  GET /api/run/status  |  GET /api/run/stream (SSE)."""

import asyncio
import json
import logging
import sys
import threading
import traceback
from pathlib import Path
from typing import AsyncGenerator


class _StateLogHandler(logging.Handler):
    """Append zdipy log records to the shared run state (SSE stream)."""

    def emit(self, record: logging.LogRecord) -> None:
        from api.state import append_log  # deferred to avoid circular import
        append_log(self.format(record))


_log_handler = _StateLogHandler()
_log_handler.setFormatter(
    logging.Formatter("%(levelname)s %(name)s: %(message)s"))
logging.getLogger("zdipy").addHandler(_log_handler)

from fastapi import APIRouter
from fastapi.responses import StreamingResponse

from api.models import RunRequest, RunStatus
from api.state import _state, _state_lock, append_log, extend_log, update_state  # noqa: F401

_ROOT = str(Path(__file__).resolve().parent.parent.parent)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

router = APIRouter(tags=["run"])

# Prevents concurrent ZDI runs.  Acquired in start_run(), released in _run_thread().
_run_lock = threading.Lock()
# Signals the running pipeline to stop at the next iteration boundary.
_stop_event = threading.Event()


# ---------------------------------------------------------------------------
# Background run thread
# ---------------------------------------------------------------------------
def _run_thread(config_path: str, forward_only: bool, verbose: int) -> None:

    def _callback(msg: str) -> None:
        append_log(msg)

    _stop_event.clear()
    try:
        from zdi_runner import run_zdi  # noqa: PLC0415
        result = run_zdi(
            config_path=config_path,
            forward_only=forward_only,
            verbose=verbose,
            progress_callback=_callback,
            stop_event=_stop_event,
        )
        # Extract H-alpha init plot before serialising the full result
        if result.halpha_init_plot is not None:
            update_state(halpha_init_plot=result.halpha_init_plot)
        _summary = {
            "iterations": result.iterations,
            "entropy": f"{result.entropy:.5f}",
            "chi2": f"{result.chi2:.6f}",
            "test": f"{result.test:.6f}",
            "converged": result.converged,
            "mean_bright": result.metadata.get("mean_bright"),
            "mean_bright_diff": result.metadata.get("mean_bright_diff"),
            "mean_mag": f"{result.metadata.get('mean_mag', 0.0):.4f} G",
        }
        update_state(status="done", result=result.to_serializable())
        extend_log([
            "",
            "=" * 60,
            "Run complete",
            "=" * 60,
            *[f"  {k:<25s}: {v}" for k, v in _summary.items()],
        ])
    except Exception:
        tb = traceback.format_exc()
        update_state(status="error", error=tb)
        append_log(f"\n\u274c Error:\n{tb}")
    finally:
        _run_lock.release()


# ---------------------------------------------------------------------------
# Endpoints
# ---------------------------------------------------------------------------
@router.post("/run")
def start_run(req: RunRequest) -> dict:
    """Start ZDI pipeline in a background thread.  Returns 'started' or 'busy'."""
    if not _run_lock.acquire(blocking=False):
        return {
            "status": "busy",
            "message": "Another run is already in progress."
        }

    config_path = req.config_path or str(Path(_ROOT) / "config.json")
    update_state(status="running",
                 log_lines=[],
                 result=None,
                 error=None,
                 halpha_init_plot=None)

    t = threading.Thread(
        target=_run_thread,
        args=(config_path, req.forward_only, req.verbose),
        daemon=True,
    )
    t.start()
    return {"status": "started", "message": "ZDI pipeline started."}


@router.get("/run/status", response_model=RunStatus)
def get_status() -> RunStatus:
    """Return current run state and last 50 log lines."""
    from api.state import get_state  # noqa: PLC0415
    st = get_state()
    return RunStatus(
        status=st["status"],
        log_tail=st["log_lines"][-50:],
        result=st["result"],
        error=st["error"],
    )


@router.get("/run/stream")
async def stream_log() -> StreamingResponse:
    """Server-Sent Events endpoint: streams log lines in real time."""

    async def _generator() -> AsyncGenerator[str, None]:
        from api.state import get_state  # noqa: PLC0415
        sent = 0
        while True:
            st = get_state()
            lines = st["log_lines"]
            status = st["status"]
            result = st["result"]
            error = st["error"]

            # Emit any new lines
            for line in lines[sent:]:
                yield f"data: {json.dumps({'type': 'log', 'text': line})}\n\n"
            sent = len(lines)

            # Emit "done" / "error" sentinel when the run finishes
            if status in ("done", "error") and sent >= len(lines):
                if status == "done":
                    yield f"data: {json.dumps({'type': 'complete', 'result': result})}\n\n"
                else:
                    yield f"data: {json.dumps({'type': 'error', 'message': error})}\n\n"
                break

            await asyncio.sleep(0.25)

    return StreamingResponse(
        _generator(),
        media_type="text/event-stream",
        headers={
            "Cache-Control": "no-cache",
            "X-Accel-Buffering": "no"
        },
    )


# ---------------------------------------------------------------------------
# GET /api/run/halpha_init_plot
# ---------------------------------------------------------------------------
@router.get("/run/halpha_init_plot")
def get_halpha_init_plot():
    """Return the H-alpha pre-processing Plotly JSON dict (if available)."""
    from api.state import get_state  # noqa: PLC0415
    from fastapi import HTTPException  # noqa: PLC0415
    st = get_state()
    plot = st.get("halpha_init_plot")
    if plot is None:
        raise HTTPException(
            status_code=404,
            detail=
            "No H-alpha init plot available. Run with halpha_compound model first.",
        )
    return plot


@router.delete("/run")
def stop_run() -> dict:
    """Signal the running ZDI pipeline to stop gracefully at the next iteration."""
    from api.state import get_state  # noqa: PLC0415
    status = get_state()["status"]
    if status != "running":
        return {
            "ok": False,
            "message": f"No run in progress (status={status!r})"
        }
    _stop_event.set()
    update_state(status="cancelled")
    return {
        "ok": True,
        "message": "Stop signal sent. Pipeline will halt at next iteration."
    }
