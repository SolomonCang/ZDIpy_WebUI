"""api/models.py — Pydantic request/response models for ZDIpy WebUI API."""

from typing import Any, Dict, List, Optional
from pydantic import BaseModel, field_validator


class RunRequest(BaseModel):
    config_path: Optional[str] = None  # None → use default config/config.json
    forward_only: bool = False
    verbose: int = 1

    @field_validator("verbose")
    @classmethod
    def _check_verbose(cls, v: int) -> int:
        if v not in (0, 1, 2):
            raise ValueError("verbose must be 0, 1, or 2")
        return v


class RunStatus(BaseModel):
    status: str  # idle | running | done | error
    log_tail: List[str]  # last N log lines
    result: Optional[Dict[str, Any]] = None
    error: Optional[str] = None


class ObservationFile(BaseModel):
    filename: str
    jdate: float
    vel_center_kms: float


class ObservationsList(BaseModel):
    files: List[ObservationFile]


class ConfigSaveResponse(BaseModel):
    ok: bool
    message: str


class ObsFileInfo(BaseModel):
    name: str
    size_bytes: int
