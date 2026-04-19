"""
LSD Interactive Web — FastAPI 后端。

提供两个 API 端点：
  GET  /api/init     — 加载 config.json 和 LSD 观测数据，返回初始参数与观测谱
  POST /api/compute  — 根据参数实时计算 Voigt 或 UR 模型轮廓

启动：
    cd ZDIpy_WebUI
    python scripts/lsd_interactive_web/server.py
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np
import uvicorn
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

# ---------------------------------------------------------------------------
# Path setup
# ---------------------------------------------------------------------------
_SCRIPT_DIR = Path(__file__).resolve().parent
_PROJECT_ROOT = _SCRIPT_DIR.parent.parent  # ZDIpy_WebUI/
sys.path.insert(0, str(_PROJECT_ROOT))

from core import readObs  # noqa: E402

# ---------------------------------------------------------------------------
# Line profile computation (identical to lsd_interactive.py)
# ---------------------------------------------------------------------------
_C_KMS: float = 2.99792458e5


def _humlicek_w4(wl_gauss_norm: np.ndarray,
                 width_lorentz: float) -> np.ndarray:
    z = width_lorentz - 1j * wl_gauss_norm
    zz = z * z
    s = np.abs(wl_gauss_norm) + width_lorentz
    w4 = np.zeros(wl_gauss_norm.shape, dtype=complex)

    c1 = s >= 15.0
    zt = z[c1]
    w4[c1] = 0.56418958355 * zt / (0.5 + zt * zt)

    c2 = (s >= 5.5) & (s < 15.0)
    zt = z[c2]
    zzt = zz[c2]
    w4[c2] = zt * (1.4104739589 + 0.56418958355 * zzt) / (
        (3.0 + zzt) * zzt + 0.7499999999)

    c3 = (width_lorentz >= 0.195 * np.abs(wl_gauss_norm) - 0.176) & (s < 5.5)
    zt = z[c3]
    w4[c3] = (16.4954955 + zt *
              (20.2093334 + zt *
               (11.9648172 + zt *
                (3.77898687 + zt * 0.564223565)))) / (16.4954955 + zt *
                                                      (38.8236274 + zt *
                                                       (39.2712051 + zt *
                                                        (21.6927370 + zt *
                                                         (6.69939801 + zt)))))

    c4 = w4 == 0.0 + 0j
    zt = z[c4]
    zzt = zz[c4]
    w4[c4] = (np.exp(zzt) - zt * (36183.30536 - zzt *
                                  (3321.990492 - zzt *
                                   (1540.786893 - zzt *
                                    (219.0312964 - zzt *
                                     (35.76682780 - zzt *
                                      (1.320521697 - zzt * 0.5641900381)))))) /
              (32066.59372 - zzt * (24322.84021 - zzt *
                                    (9022.227659 - zzt *
                                     (2186.181081 - zzt *
                                      (364.2190727 - zzt *
                                       (61.57036588 - zzt *
                                        (1.841438936 - zzt))))))))
    return w4


def voigt_profile(vel: np.ndarray, v_center: float, line_strength: float,
                  gauss_width: float, lorentz_frac: float, wl0_nm: float,
                  inst_res: float) -> np.ndarray:
    if inst_res > 0:
        vi = _C_KMS / inst_res * 0.6005612043932249
        eg = np.sqrt(gauss_width**2 + vi**2)
        el = lorentz_frac * (gauss_width / eg)
        es = line_strength * (gauss_width / eg)
    else:
        eg, el, es = gauss_width, lorentz_frac, line_strength

    wgw = eg / _C_KMS * wl0_nm
    wl = (vel - v_center) / _C_KMS * wl0_nm + wl0_nm
    wgn = (wl0_nm - wl) / wgw

    w4 = _humlicek_w4(wgn, el)
    return (1.0 - es * w4.real).astype(float)


def unno_profile(vel: np.ndarray, v_center: float, line_strength: float,
                 gauss_width: float, lorentz_frac: float, wl0_nm: float,
                 inst_res: float, beta: float, fI: float) -> np.ndarray:
    if inst_res > 0:
        vi = _C_KMS / inst_res * 0.6005612043932249
        eg = np.sqrt(gauss_width**2 + vi**2)
        el = lorentz_frac * (gauss_width / eg)
    else:
        eg, el = gauss_width, lorentz_frac

    wgw = eg / _C_KMS * wl0_nm
    wl = (vel - v_center) / _C_KMS * wl0_nm + wl0_nm
    wgn = (wl - wl0_nm) / wgw
    norm = 1.0 / (np.sqrt(np.pi) * eg)

    W = _humlicek_w4(-wgn, el)
    eta_comp = W.real * norm
    etaI = 1.0 + line_strength * eta_comp
    beta_mu = beta
    sI = (1.0 + beta_mu / etaI) / (1.0 + beta_mu)
    return sI.astype(float)


# ---------------------------------------------------------------------------
# Load observation data once at startup
# ---------------------------------------------------------------------------
_CONFIG_PATH = _PROJECT_ROOT / "config.json"

with open(_CONFIG_PATH, "r") as f:
    _cfg = json.load(f)

_obs_files = _cfg["observations"]["files"]
_file_list = [str(_PROJECT_ROOT / entry["filename"]) for entry in _obs_files]
_vel_rs = np.array([entry["vel_center_kms"] for entry in _obs_files])
_line_cfg = _cfg["line_model"]
_inst_res = float(_cfg["instrument"]["spectral_resolution"])

print("Loading LSD observations…")
_obsSet = readObs.obsProfSet(_file_list)

# Build median Stokes I
_vel_ref = _obsSet[0].wl.copy()
_aligned_I = []
for _i, _obs in enumerate(_obsSet):
    _shifted = _obs.wl - _vel_rs[_i]
    _aligned_I.append(np.interp(_vel_ref, _shifted, _obs.specI))
_aligned_I = np.array(_aligned_I)
_median_I = np.median(_aligned_I, axis=0)
_v_dense = np.linspace(float(_vel_ref[0]), float(_vel_ref[-1]), 800)

# Initial parameters from config
_wl0_nm = float(_line_cfg["wavelength_nm"])
_beta_cfg = float(_line_cfg.get("unno_beta", -1))
if _beta_cfg <= 0:
    _ld = float(_line_cfg["limb_darkening"])
    if abs(_ld - 1.0) < 1e-10:
        _ld = 1.0 - 1e-6
    _beta_cfg = _ld / (1.0 - _ld)
_fI_cfg = float(_line_cfg.get("unno_filling_factor_I", 1.0))
_model_type = _line_cfg.get("model_type", "voigt").lower()
if _model_type not in ("voigt", "unno"):
    _model_type = "voigt"

print(f"Ready — {len(_obsSet)} observations, model_type={_model_type}")

# ---------------------------------------------------------------------------
# FastAPI application
# ---------------------------------------------------------------------------
app = FastAPI(title="LSD Interactive Fit")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)


class ComputeRequest(BaseModel):
    model_type: str = "voigt"
    line_strength: float = 0.4
    gauss_width_kms: float = 20.0
    lorentz_width_fraction: float = 0.1
    v_center: float = 0.0
    inst_resolution: float = 65000.0
    limb_darkening: float = 0.3
    beta: float = 0.4286
    fI: float = 1.0


@app.get("/api/init")
def api_init():
    """返回初始参数和观测数据。"""
    return {
        "vel_obs": _vel_ref.tolist(),
        "median_I": _median_I.tolist(),
        "individual_I": _aligned_I.tolist(),
        "vel_model": _v_dense.tolist(),
        "params": {
            "model_type": _model_type,
            "line_strength": float(_line_cfg["line_strength"]),
            "gauss_width_kms": float(_line_cfg["gauss_width_kms"]),
            "lorentz_width_fraction":
            float(_line_cfg["lorentz_width_fraction"]),
            "v_center": float(_vel_rs[0]),
            "inst_resolution": _inst_res,
            "limb_darkening": float(_line_cfg["limb_darkening"]),
            "beta": _beta_cfg,
            "fI": _fI_cfg,
        },
        "slider_ranges": {
            "line_strength": {
                "min": 0.01,
                "max": 2.0,
                "step": 0.01,
                "max_unno": 100.0,
                "step_unno": 0.5,
            },
            "gauss_width_kms": {
                "min": 1.0,
                "max": 80.0,
                "step": 0.5
            },
            "lorentz_width_fraction": {
                "min": 0.0,
                "max": 1.5,
                "step": 0.01
            },
            "v_center": {
                "min": -80.0,
                "max": 80.0,
                "step": 0.5
            },
            "inst_resolution": {
                "min": 5000,
                "max": 200000,
                "step": 1000
            },
            "limb_darkening": {
                "min": 0.0,
                "max": 1.0,
                "step": 0.01
            },
            "beta": {
                "min": 0.01,
                "max": 10.0,
                "step": 0.01
            },
            "fI": {
                "min": 0.01,
                "max": 1.0,
                "step": 0.01
            },
        },
    }


@app.post("/api/compute")
def api_compute(req: ComputeRequest):
    """根据参数计算模型轮廓，返回模型 I 和 RMS。"""
    if req.model_type == "unno":
        I_model = unno_profile(_v_dense, req.v_center, req.line_strength,
                               req.gauss_width_kms, req.lorentz_width_fraction,
                               _wl0_nm, req.inst_resolution, req.beta, req.fI)
    else:
        I_model = voigt_profile(_v_dense, req.v_center, req.line_strength,
                                req.gauss_width_kms,
                                req.lorentz_width_fraction, _wl0_nm,
                                req.inst_resolution)

    # RMS in ±3σ
    I_at_obs = np.interp(_vel_ref, _v_dense, I_model)
    mask = np.abs(_vel_ref - req.v_center) <= req.gauss_width_kms * 3.0
    rms = None
    if np.any(mask):
        rms = float(np.sqrt(np.mean((_median_I[mask] - I_at_obs[mask])**2)))

    return {
        "I_model": I_model.tolist(),
        "rms": rms,
    }


class AutoEstimateRequest(BaseModel):
    model_type: str = "voigt"
    gauss_width_kms: float = 20.0
    lorentz_width_fraction: float = 0.1
    v_center: float = 0.0
    inst_resolution: float = 65000.0
    limb_darkening: float = 0.3
    beta: float = 0.4286
    fI: float = 1.0


@app.post("/api/auto_estimate")
def api_auto_estimate(req: AutoEstimateRequest):
    """在 ±3σ 窗口内，通过最小化 RMS 自动估计 line_strength。"""
    from scipy.optimize import minimize_scalar

    mask = np.abs(_vel_ref - req.v_center) <= req.gauss_width_kms * 3.0
    if not np.any(mask):
        return {"line_strength": 0.4, "I_model": None, "rms": None}

    obs_masked = _median_I[mask]

    def rms_func(ls):
        if req.model_type == "unno":
            I_m = unno_profile(_v_dense, req.v_center, ls, req.gauss_width_kms,
                               req.lorentz_width_fraction, _wl0_nm,
                               req.inst_resolution, req.beta, req.fI)
        else:
            I_m = voigt_profile(_v_dense, req.v_center, ls,
                                req.gauss_width_kms,
                                req.lorentz_width_fraction, _wl0_nm,
                                req.inst_resolution)
        I_at_obs = np.interp(_vel_ref, _v_dense, I_m)
        return float(np.sqrt(np.mean((obs_masked - I_at_obs[mask])**2)))

    # Use the slider range as search bounds
    ls_range = api_init()["slider_ranges"]["line_strength"]
    ls_max = ls_range.get("max_unno", ls_range["max"]) \
        if req.model_type == "unno" else ls_range["max"]
    result = minimize_scalar(rms_func,
                             bounds=(ls_range["min"], ls_max),
                             method="bounded")
    best_ls = float(result.x)

    # Compute final model at best line_strength
    if req.model_type == "unno":
        I_model = unno_profile(_v_dense, req.v_center, best_ls,
                               req.gauss_width_kms, req.lorentz_width_fraction,
                               _wl0_nm, req.inst_resolution, req.beta, req.fI)
    else:
        I_model = voigt_profile(_v_dense, req.v_center, best_ls,
                                req.gauss_width_kms,
                                req.lorentz_width_fraction, _wl0_nm,
                                req.inst_resolution)

    I_at_obs = np.interp(_vel_ref, _v_dense, I_model)
    rms = float(np.sqrt(np.mean((obs_masked - I_at_obs[mask])**2)))

    return {
        "line_strength": round(best_ls, 5),
        "I_model": I_model.tolist(),
        "rms": rms,
    }


class SaveConfigRequest(BaseModel):
    model_type: str = "voigt"
    line_strength: float = 0.4
    gauss_width_kms: float = 20.0
    lorentz_width_fraction: float = 0.1
    limb_darkening: float = 0.3
    unno_beta: float = 0.4286
    unno_filling_factor_I: float = 1.0


@app.post("/api/save_config")
def api_save_config(req: SaveConfigRequest):
    """将当前参数写回 config.json 的 line_model 段。"""
    with open(_CONFIG_PATH, "r") as f:
        cfg = json.load(f)

    lm = cfg["line_model"]
    lm["model_type"] = req.model_type
    lm["line_strength"] = round(req.line_strength, 5)
    lm["gauss_width_kms"] = round(req.gauss_width_kms, 3)
    lm["lorentz_width_fraction"] = round(req.lorentz_width_fraction, 5)
    lm["limb_darkening"] = round(req.limb_darkening, 4)
    lm["unno_beta"] = round(req.unno_beta, 5)
    lm["unno_filling_factor_I"] = round(req.unno_filling_factor_I, 4)

    with open(_CONFIG_PATH, "w") as f:
        json.dump(cfg, f, indent=2)
        f.write("\n")

    return {"ok": True}


# Serve static frontend files
_STATIC_DIR = _SCRIPT_DIR / "static"
app.mount("/",
          StaticFiles(directory=str(_STATIC_DIR), html=True),
          name="static")

# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="LSD Interactive Fit — Web UI")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8051)
    args = parser.parse_args()
    print(f"Starting server at http://{args.host}:{args.port}")
    uvicorn.run(app, host=args.host, port=args.port, log_level="warning")
