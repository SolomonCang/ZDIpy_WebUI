"""api/state.py — Shared mutable run state for ZDIpy WebUI API.

All route modules must read/write run state exclusively through the
functions in this module.  Never import _state or _state_lock
directly from api.routes.run.
"""
import threading
from typing import Any, Dict, List

_state_lock: threading.Lock = threading.Lock()
_state: Dict[str, Any] = {
    "status": "idle",
    "log_lines": [],
    "result": None,
    "error": None,
    "halpha_init_plot": None,
}


def get_state() -> Dict[str, Any]:
    with _state_lock:
        return dict(_state)


def update_state(**kwargs: Any) -> None:
    with _state_lock:
        _state.update(kwargs)


def append_log(msg: str) -> None:
    with _state_lock:
        _state["log_lines"].append(msg)


def extend_log(lines: List[str]) -> None:
    with _state_lock:
        _state["log_lines"].extend(lines)


def reset_state() -> None:
    with _state_lock:
        _state["status"] = "idle"
        _state["log_lines"] = []
        _state["result"] = None
        _state["error"] = None
        _state["halpha_init_plot"] = None
