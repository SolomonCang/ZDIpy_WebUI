"""api/routes/run.py — POST /api/run  |  GET /api/run/status  |  GET /api/run/stream (SSE)."""

import asyncio
import json
import sys
import threading
import traceback
from pathlib import Path
from typing import AsyncGenerator

from fastapi import APIRouter
from fastapi.responses import StreamingResponse

from api.models import RunRequest, RunStatus

_ROOT = str(Path(__file__).resolve().parent.parent.parent)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

router = APIRouter(tags=["run"])

# ---------------------------------------------------------------------------
# Shared mutable run state (protected by _state_lock)
# ---------------------------------------------------------------------------
_state_lock = threading.Lock()
_state: dict = {
    "status": "idle",  # idle | running | done | error
    "log_lines": [],
    "result": None,
    "error": None,
}

# Prevents concurrent ZDI runs.  Acquired in start_run(), released in _run_thread().
_run_lock = threading.Lock()


# ---------------------------------------------------------------------------
# Stdout/stderr capture
# ---------------------------------------------------------------------------
class _StreamCapture:
    """Redirect writes to both _state["log_lines"] and the original stream."""
    def __init__(self, original):
        self._original = original

    def write(self, s: str) -> None:
        if s.strip():
            with _state_lock:
                _state["log_lines"].append(s.rstrip())
        self._original.write(s)

    def flush(self) -> None:
        self._original.flush()


# ---------------------------------------------------------------------------
# Background run thread
# ---------------------------------------------------------------------------
def _run_thread(config_path: str, forward_only: bool, verbose: int) -> None:
    old_stdout = sys.stdout
    old_stderr = sys.stderr
    sys.stdout = _StreamCapture(old_stdout)
    sys.stderr = _StreamCapture(old_stderr)
    try:
        from zdi_runner import run_zdi  # noqa: PLC0415
        result = run_zdi(
            config_path=config_path,
            forward_only=forward_only,
            verbose=verbose,
        )
        with _state_lock:
            _state["status"] = "done"
            _state["result"] = result
            _state["log_lines"] += [
                "",
                "=" * 60,
                "Run complete",
                "=" * 60,
                *[f"  {k:<25s}: {v}" for k, v in result.items()],
            ]
    except Exception:
        tb = traceback.format_exc()
        with _state_lock:
            _state["status"] = "error"
            _state["error"] = tb
            _state["log_lines"].append(f"\n❌ Error:\n{tb}")
    finally:
        sys.stdout = old_stdout
        sys.stderr = old_stderr
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

    config_path = req.config_path or str(
        Path(_ROOT) / "config" / "config.json")
    with _state_lock:
        _state["status"] = "running"
        _state["log_lines"] = []
        _state["result"] = None
        _state["error"] = None

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
    with _state_lock:
        return RunStatus(
            status=_state["status"],
            log_tail=_state["log_lines"][-50:],
            result=_state["result"],
            error=_state["error"],
        )


@router.get("/run/stream")
async def stream_log() -> StreamingResponse:
    """Server-Sent Events endpoint: streams log lines in real time."""
    async def _generator() -> AsyncGenerator[str, None]:
        sent = 0
        while True:
            with _state_lock:
                lines = list(_state["log_lines"])
                status = _state["status"]
                result = _state["result"]
                error = _state["error"]

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
