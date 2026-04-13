"""
halpha_interactive.py — H-alpha 双 Voigt 模型交互式参数调整工具。

用法（在项目根目录下）：
    python tests/halpha_interactive.py

启动后会弹出 matplotlib 窗口：
  - 上方主图：中值 Stokes I（蓝）、发射成分（橙虚）、自吸收成分（绿虚）、
    复合模型（红实）
  - 下方 7 个滑轨：实时控制所有模型参数
  - [Reset] 按钮：恢复自动估算的初始值
  - [Save PNG] 按钮：将当前状态保存为 PNG 文件（与本脚本同目录）
"""
from __future__ import annotations

import json
import os
import sys

import matplotlib.pyplot as plt
import matplotlib.widgets as mwidgets
import numpy as np

# 把项目根目录加入 sys.path，确保 core/ 可以 import
sys.path.insert(0,
                os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from core import readObs
from core.line_models import halpha_preproc
from core.line_models.halpha_preproc import _build_median_spectrum, _voigt_real

# ---------------------------------------------------------------------------
# 1. 加载观测数据 & 自动估算初始参数
# ---------------------------------------------------------------------------
_BASE = os.path.dirname(__file__)
_CONFIG = os.path.join(_BASE, "..", "config.json")

with open(_CONFIG, "r") as _f:
    _cfg = json.load(_f)

_obs_files = _cfg["observations"]["files"]
_file_list = [os.path.join(_BASE, "..", f["filename"]) for f in _obs_files]
_vel_rs = np.array([f["vel_center_kms"] for f in _obs_files])
_jdates = np.array([f["jdate"] for f in _obs_files])
_jdate_ref = float(_cfg["observations"]["jdate_ref"])

_star = _cfg["star"]
_line_cfg = _cfg["line_model"]
_inc_rad = float(_star["inclination_deg"]) / 180.0 * np.pi
_vsini = float(_star["vsini_kms"])
_vel_eq = _vsini / np.sin(_inc_rad)

_stellar_params = {
    "vel_eq_kms": _vel_eq,
    "inc_rad": _inc_rad,
    "period_days": float(_star["period_days"]),
    "d_omega": float(_star["differential_rotation_rad_per_day"]),
    "jdates": _jdates,
    "jdate_ref": _jdate_ref,
    "wl0_nm": float(_line_cfg["wavelength_nm"]),
    "lande_g": float(_line_cfg["lande_g"]),
    "limb_dark": float(_line_cfg.get("limb_darkening", 0.0)),
    "grav_dark": float(_line_cfg.get("gravity_darkening", 0.0)),
    "fV": float(_line_cfg.get("filling_factor_V", 1.0)),
    "inst_res": float(_cfg["instrument"]["spectral_resolution"]),
}

print("Loading observations…")
_obsSet = readObs.obsProfSet(_file_list)

print("Running auto parameter estimation (disk-integrated forward model)…")
_result = halpha_preproc.auto_estimate_halpha_params(
    _obsSet, _vel_rs, None, stellar_params=_stellar_params)

_p0 = _result["params"]
_v_ctr0 = _result["v_center"]

# 中值谱（用于绘图背景）
_vel, _median_I, _indiv_I, _ = _build_median_spectrum(_obsSet, _vel_rs)

# 提取盘积分前向模型均值曲线（用于静态参考线）
_vel_fwd: np.ndarray | None = None
_I_fwd: np.ndarray | None = None
for _trace in _result.get("plot_data", {}).get("data", []):
    if "前向模型均值" in _trace.get("name", ""):
        _vel_fwd = np.array(_trace["x"])
        _I_fwd = np.array(_trace["y"])
        break

print("Auto-estimated parameters:")
for k, v in _p0.items():
    print(f"  {k}: {v:.5g}")
print(f"  v_center: {_v_ctr0:.3f} km/s")
print(f"  fit_ok: {_result['fit_ok']}")
if _vel_fwd is not None:
    print("  disk-integrated forward model mean: loaded")
else:
    print(
        "  disk-integrated forward model mean: not available (forward model may have failed)"
    )

# ---------------------------------------------------------------------------
# 2. 复合模型计算（供实时更新调用）
# ---------------------------------------------------------------------------
_v_dense = np.linspace(float(_vel[0]), float(_vel[-1]), 800)


def _compute_model(A_em, sig_em, a_em, A_abs, sig_abs, a_abs, v_ctr):
    u_em = (_v_dense - v_ctr) / max(sig_em, 0.1)
    u_abs = (_v_dense - v_ctr) / max(sig_abs, 0.1)
    H_em = _voigt_real(u_em, max(0.0, a_em))
    H_abs = _voigt_real(u_abs, max(0.0, a_abs))
    I_em = 1.0 + A_em * H_em
    I_abs_comp = A_abs * H_abs
    I_fit = I_em - I_abs_comp
    I_abs_trace = 1.0 - I_abs_comp
    return I_em, I_abs_trace, I_fit


# ---------------------------------------------------------------------------
# 3. 构建图形布局
# ---------------------------------------------------------------------------
_SLIDER_DEFS = [
    # (label, key,          init_val, vmin,  vmax,  vstep,  fmt)
    ("A_em  (emission strength)", "emission_strength",
     _p0["emission_strength"], 0.01, 6.0, 0.01, "%.3f"),
    ("σ_em  (emission width km/s)", "emission_gauss_kms",
     _p0["emission_gauss_kms"], 1.0, 120.0, 0.5, "%.1f"),
    ("a_em  (Lorentz ratio)", "emission_lorentz_ratio",
     _p0["emission_lorentz_ratio"], 0.0, 1.5, 0.01, "%.3f"),
    ("A_abs (absorption strength)", "absorption_strength",
     _p0["absorption_strength"], 0.0, 4.0, 0.01, "%.3f"),
    ("σ_abs (absorption width km/s)", "absorption_gauss_kms",
     _p0["absorption_gauss_kms"], 1.0, 80.0, 0.5, "%.1f"),
    ("a_abs (abs Lorentz ratio)", "absorption_lorentz_ratio",
     _p0["absorption_lorentz_ratio"], 0.0, 1.5, 0.01, "%.3f"),
    ("v_center  (km/s)", "_v_center", _v_ctr0, -80.0, 80.0, 0.5, "%.1f"),
]

_N_SLIDERS = len(_SLIDER_DEFS)
_SLIDER_H = 0.032  # height of each slider ax
_SLIDER_GAP = 0.012  # gap between sliders
_BTN_H = 0.04
_BTN_W = 0.10
_BOTTOM_MARGIN = 0.04
_TOTAL_SLIDER_H = _N_SLIDERS * (_SLIDER_H +
                                _SLIDER_GAP) + _BTN_H + _BOTTOM_MARGIN + 0.02

fig = plt.figure(figsize=(10, 8))
fig.patch.set_facecolor("#1e1e2e")

# Main plot axes (leaves room at bottom for sliders)
_ax_main_top = 0.97
_ax_main_bot = _TOTAL_SLIDER_H + 0.04
ax = fig.add_axes([0.10, _ax_main_bot, 0.88, _ax_main_top - _ax_main_bot])
ax.set_facecolor("#1e1e2e")
for spine in ax.spines.values():
    spine.set_edgecolor("#585b70")
ax.tick_params(colors="#cdd6f4")
ax.xaxis.label.set_color("#cdd6f4")
ax.yaxis.label.set_color("#cdd6f4")
ax.title.set_color("#cdd6f4")

# Plot static elements
_line_obs, = ax.plot(_vel,
                     _median_I,
                     color="#89b4fa",
                     lw=1.8,
                     label="Median Stokes I",
                     zorder=3)

# Dynamic model lines (initial values)
_I_em0, _I_abs0, _I_fit0 = _compute_model(
    _p0["emission_strength"],
    _p0["emission_gauss_kms"],
    _p0["emission_lorentz_ratio"],
    _p0["absorption_strength"],
    _p0["absorption_gauss_kms"],
    _p0["absorption_lorentz_ratio"],
    _v_ctr0,
)
_line_em, = ax.plot(_v_dense,
                    _I_em0,
                    "--",
                    color="#fab387",
                    lw=1.5,
                    label="Emission Component",
                    zorder=4)
_line_abs, = ax.plot(_v_dense,
                     _I_abs0,
                     "--",
                     color="#a6e3a1",
                     lw=1.5,
                     label="Absorption Component",
                     zorder=4)
_line_fit, = ax.plot(_v_dense,
                     _I_fit0,
                     "-",
                     color="#f38ba8",
                     lw=2.2,
                     label="Local Voigt Model",
                     zorder=5)

# 静态盘积分前向模型均值线（橙黄色实线）
_line_fwd = None
if _vel_fwd is not None and _I_fwd is not None:
    _label_fwd = ("Disk-Integrated Mean (fit)"
                  if _result["fit_ok"] else "Disk-Integrated Mean (est.)")
    (_line_fwd, ) = ax.plot(_vel_fwd,
                            _I_fwd,
                            "-",
                            color="#f9e2af",
                            lw=2.0,
                            alpha=0.85,
                            label=_label_fwd,
                            zorder=6)

ax.axhline(1.0, color="#585b70", lw=0.8, ls="--", zorder=1)
ax.set_xlabel("Velocity (km/s)")
ax.set_ylabel("Stokes I")
_fit_mode = "disk-integrated" if _result["fit_ok"] else "morphological"
ax.set_title(f"H-alpha Dual-Voigt Interactive Fit  [{_fit_mode}]")
leg = ax.legend(facecolor="#313244",
                edgecolor="#585b70",
                labelcolor="#cdd6f4",
                fontsize=9)

# Residual info text
_resid_label = f"fit_ok={_result['fit_ok']}"
_txt_status = ax.text(0.02,
                      0.97,
                      _resid_label,
                      transform=ax.transAxes,
                      va="top",
                      ha="left",
                      fontsize=8,
                      color="#a6adc8")

# ---------------------------------------------------------------------------
# 4. 创建滑轨
# ---------------------------------------------------------------------------
_slider_ax_list = []
_slider_list = []

_SLIDER_LEFT = 0.15
_SLIDER_WIDTH = 0.72

for i, (label, key, init, vmin, vmax, vstep, fmt) in enumerate(_SLIDER_DEFS):
    # bottom of this slider (count from top: slider 0 is highest)
    _bot = _BOTTOM_MARGIN + _BTN_H + 0.02 + (_N_SLIDERS - 1 -
                                             i) * (_SLIDER_H + _SLIDER_GAP)
    _sax = fig.add_axes([_SLIDER_LEFT, _bot, _SLIDER_WIDTH, _SLIDER_H],
                        facecolor="#313244")
    _sax.tick_params(colors="#cdd6f4")
    sl = mwidgets.Slider(
        _sax,
        label,
        vmin,
        vmax,
        valinit=float(np.clip(init, vmin, vmax)),
        valstep=vstep,
        color="#89dceb",
    )
    sl.label.set_color("#cdd6f4")
    sl.label.set_fontsize(8)
    sl.valtext.set_color("#f5c2e7")
    sl.valtext.set_fontsize(8)
    _slider_ax_list.append(_sax)
    _slider_list.append(sl)

# ---------------------------------------------------------------------------
# 5. 按钮
# ---------------------------------------------------------------------------
_btn_bot = _BOTTOM_MARGIN
_ax_btn_reset = fig.add_axes([0.15, _btn_bot, _BTN_W, _BTN_H],
                             facecolor="#45475a")
_ax_btn_save = fig.add_axes([0.27, _btn_bot, _BTN_W, _BTN_H],
                            facecolor="#45475a")

_btn_reset = mwidgets.Button(_ax_btn_reset,
                             "Reset",
                             color="#45475a",
                             hovercolor="#585b70")
_btn_save = mwidgets.Button(_ax_btn_save,
                            "Save PNG",
                            color="#45475a",
                            hovercolor="#585b70")
_btn_reset.label.set_color("#cdd6f4")
_btn_save.label.set_color("#cdd6f4")

# CheckButtons: 控制盘积分均值线的显示/隐藏
_ax_chk = fig.add_axes([0.39, _btn_bot, 0.22, _BTN_H], facecolor="#313244")
_chk_fwd = mwidgets.CheckButtons(
    _ax_chk,
    ["Disk-Integrated Mean"],
    actives=[_line_fwd is not None],
)
_chk_fwd.labels[0].set_color("#f9e2af")
_chk_fwd.labels[0].set_fontsize(8)


# ---------------------------------------------------------------------------
# 6. 回调函数
# ---------------------------------------------------------------------------
def _get_slider_vals():
    vals = {}
    for sl, (_, key, *_rest) in zip(_slider_list, _SLIDER_DEFS):
        vals[key] = sl.val
    return vals


def _update(_=None):
    v = _get_slider_vals()
    A_em = v["emission_strength"]
    sig_em = v["emission_gauss_kms"]
    a_em = v["emission_lorentz_ratio"]
    A_abs = v["absorption_strength"]
    sig_abs = v["absorption_gauss_kms"]
    a_abs = v["absorption_lorentz_ratio"]
    v_ctr = v["_v_center"]

    I_em, I_abs_tr, I_fit = _compute_model(A_em, sig_em, a_em, A_abs, sig_abs,
                                           a_abs, v_ctr)
    _line_em.set_ydata(I_em)
    _line_abs.set_ydata(I_abs_tr)
    _line_fit.set_ydata(I_fit)

    # Update legend labels with current values
    _line_em.set_label(
        f"Emission  A={A_em:.2f}, σ={sig_em:.0f} km/s, a={a_em:.2f}")
    _line_abs.set_label(f"Absorption  A={A_abs:.2f}, σ={sig_abs:.0f} km/s")
    _line_fit.set_label(f"Model  v_c={v_ctr:+.1f} km/s")
    ax.legend(facecolor="#313244",
              edgecolor="#585b70",
              labelcolor="#cdd6f4",
              fontsize=8)

    # Residual RMS in the emission region
    I_fit_at_obs = np.interp(_vel, _v_dense, I_fit)
    mask = np.abs(_vel - v_ctr) <= sig_em * 3.0
    if np.any(mask):
        rms = float(np.sqrt(np.mean(
            (_median_I[mask] - I_fit_at_obs[mask])**2)))
        _txt_status.set_text(f"RMS in ±3σ = {rms:.4f}")
    else:
        _txt_status.set_text("")

    fig.canvas.draw_idle()


def _reset(_=None):
    for sl, (_, key, init, vmin, vmax, *_) in zip(_slider_list, _SLIDER_DEFS):
        sl.set_val(float(np.clip(init, vmin, vmax)))
    # _update is triggered automatically via slider observers


def _toggle_fwd(_=None):
    if _line_fwd is not None:
        _line_fwd.set_visible(_chk_fwd.get_status()[0])
        ax.legend(facecolor="#313244",
                  edgecolor="#585b70",
                  labelcolor="#cdd6f4",
                  fontsize=8)
        fig.canvas.draw_idle()


def _save_png(_=None):
    save_path = os.path.join(_BASE, "halpha_interactive_result.png")
    fig.savefig(save_path, dpi=200, facecolor=fig.get_facecolor())
    print(f"Saved: {save_path}")


for sl in _slider_list:
    sl.on_changed(_update)

_btn_reset.on_clicked(_reset)
_btn_save.on_clicked(_save_png)
_chk_fwd.on_clicked(_toggle_fwd)

# Initial draw
_update()

plt.show()
