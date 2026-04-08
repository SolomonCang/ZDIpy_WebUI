"""core.line_models.halpha_preproc — H-alpha 发射线预处理工具。

两个辅助函数，在加载观测数据之后、主反演循环之前调用：

- ``normalize_halpha_emission``  : 将所有历元的 Stokes I 发射峰高度
  等比例标定到中值水平，并对 V/N 及所有噪声数组施加相同缩放。

- ``auto_estimate_halpha_params``: 对归一化后的 Stokes I 求中值得到
  平均发射谱，然后对双 Voigt 模型做最小二乘拟合，估算初始参数，
  并返回 Plotly 格式的绘图 dict。
"""
from __future__ import annotations

from typing import Callable, Optional

import numpy as np

_C_KMS: float = 2.99792458e5

# ---------------------------------------------------------------------------
# 公共 Voigt 辅助（直接借调 unno 中的实现）
# ---------------------------------------------------------------------------


def _voigt_real(u: np.ndarray, a: float) -> np.ndarray:
    """Humlicek Voigt 轮廓（实部），接受 1-D 数组。"""
    from core.line_models.unno import _voigt_faraday_humlicek  # noqa: PLC0415
    u_arr = np.asarray(u, dtype=float).ravel()
    W = _voigt_faraday_humlicek(u_arr, float(a))
    return W.real.reshape(np.asarray(u).shape)


# ===========================================================================
# 模块 1：发射强度归一化
# ===========================================================================


def normalize_halpha_emission(
    obsSet: list,
    vel_rs: np.ndarray,
    log_fn: Optional[Callable[[str], None]] = None,
) -> np.ndarray:
    """对所有历元的 H-alpha 发射峰高度进行归一化。

    针对 H-alpha 发射线性质（specI 在线心大于连续谱），计算每个历元的
    超出连续谱幅度 ``h_i = max(specI) − 1``，以所有历元 h_i 的**中值**
    作为参考，各历元乘以缩放因子 ``s_i = h_ref / h_i``。

    Stokes V、Null（N）及对应的噪声数组（specIsig / specVsig / specNsig）
    均乘以相同的 s_i，保持信噪比不变。

    Parameters
    ----------
    obsSet : list of obsProf
        已读取的观测数组（in-place 修改）。
    vel_rs : (N_obs,) ndarray
        各历元恒星视向速度（km/s），仅用于日志输出。
    log_fn : callable, optional
        字符串消息回调，与 ZDIPipeline._log 接口相同。

    Returns
    -------
    scale_factors : (N_obs,) ndarray
        各历元的乘法缩放因子。
    """

    def _log(msg: str) -> None:
        if log_fn is not None:
            log_fn(msg)

    # --- 计算各历元发射峰高度 ---
    heights = np.array(
        [max(0.0,
             float(np.max(obs.specI)) - 1.0) for obs in obsSet])
    valid = heights > 1e-6
    if not np.any(valid):
        _log("  [normalize_halpha] 未检测到有效发射峰，跳过归一化。")
        return np.ones(len(obsSet))

    ref_h = float(np.median(heights[valid]))
    scale_factors = np.where(valid, ref_h / heights, 1.0)

    _log("H-alpha 发射强度归一化：")
    _log(f"  参考峰高（中值）: {ref_h:.4f}")
    for i, (obs, sf) in enumerate(zip(obsSet, scale_factors)):
        _log(f"  历元 {i:02d}  峰高={heights[i]:.4f}  "
             f"v_r={float(vel_rs[i]):+.2f} km/s  缩放因子={sf:.4f}")

    # --- 在原数组上施加缩放 ---
    for obs, sf in zip(obsSet, scale_factors):
        if abs(sf - 1.0) < 1e-9:
            continue
        obs.specI = 1.0 + (obs.specI - 1.0) * sf
        obs.specIsig = obs.specIsig * sf
        obs.specV = obs.specV * sf
        obs.specVsig = obs.specVsig * sf
        obs.specN = obs.specN * sf
        obs.specNsig = obs.specNsig * sf

    return scale_factors


# ===========================================================================
# 模块 2：自动参数估算
# ===========================================================================

# ---------------------------------------------------------------------------
# 内部辅助
# ---------------------------------------------------------------------------


def _build_median_spectrum(
    obsSet: list,
    vel_rs: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, list]:
    """将所有观测插值到公共速度网格，返回中值谱。

    Returns
    -------
    vel_common : (M,) ndarray
        公共速度格点（km/s，恒星静止系）。
    median_I : (M,) ndarray
        各格点 Stokes I 中值。
    indiv_I : list of (M,) ndarray
        各历元插值后的 Stokes I（归一）。
    """
    # 各历元移到恒星静止系
    rest_vels = [obs.wl - float(vel_rs[i]) for i, obs in enumerate(obsSet)]

    # 公共速度范围：取所有历元的交集（保守）
    v_min = max(float(np.min(v)) for v in rest_vels)
    v_max = min(float(np.max(v)) for v in rest_vels)

    # 速度步长：取最细采样的步长
    dv_arr = [
        float(np.mean(np.abs(np.diff(v)))) for v in rest_vels if len(v) > 1
    ]
    dv = min(dv_arr) if dv_arr else 1.0

    if v_max <= v_min:
        # 退化情况：使用第一个观测的速度网格
        vel_common = rest_vels[0]
        indiv_I = [obs.specI.copy() for obs in obsSet]
        median_I = np.median(np.vstack(indiv_I), axis=0)
        return vel_common, median_I, indiv_I

    vel_common = np.arange(v_min, v_max + 0.5 * dv, dv)
    indiv_I = []
    for i, obs in enumerate(obsSet):
        v_i = rest_vels[i]
        # 边界值采用最近端点
        I_interp = np.interp(vel_common,
                             v_i,
                             obs.specI,
                             left=obs.specI[0],
                             right=obs.specI[-1])
        indiv_I.append(I_interp)

    median_I = np.median(np.vstack(indiv_I), axis=0)
    return vel_common, median_I, indiv_I


# ---------------------------------------------------------------------------
# 主函数
# ---------------------------------------------------------------------------


def auto_estimate_halpha_params(
    obsSet: list,
    vel_rs: np.ndarray,
    current_lineData,
    log_fn: Optional[Callable[[str], None]] = None,
) -> dict:
    """从中值发射谱估算 H-alpha 双 Voigt 初始参数。

    算法步骤：

    1. 对所有历元 Stokes I 在恒星静止系插值至公共速度格点。
    2. 取中值得到代表性发射谱。
    3. 用简单形态法估算初值（峰高、半宽、中心速度、吸收深度）。
    4. 用 ``scipy.optimize.curve_fit`` 对双 Voigt 模型做非线性最小二乘拟合。
    5. 构建 Plotly 格式绘图数据（个别轮廓 + 中值 + 分量 + 拟合模型）并返回。

    Parameters
    ----------
    obsSet : list of obsProf
        （已归一化的）观测数组。
    vel_rs : (N_obs,) ndarray
        各历元恒星视向速度（km/s）。
    current_lineData : lineDataHalpha
        现有线参数，用于拟合失败时的回退值。
    log_fn : callable, optional

    Returns
    -------
    dict
        ``params``    — 最优参数 dict（键名与 config.json ``line_model`` 块一致）
        ``plot_data`` — Plotly JSON dict（``data`` + ``layout``）
        ``fit_ok``    — bool，拟合是否收敛
    """
    from scipy.optimize import curve_fit  # noqa: PLC0415

    def _log(msg: str) -> None:
        if log_fn is not None:
            log_fn(msg)

    _log("自动估算 H-alpha 复合模型初始参数…")

    # --- 构建中值谱 ---
    vel, median_I, indiv_I = _build_median_spectrum(obsSet, vel_rs)
    _log(f"  速度范围: [{float(vel[0]):.1f}, {float(vel[-1]):.1f}] km/s，"
         f"格点数: {len(vel)}")

    # --- 初值估算 ---
    i_peak = int(np.argmax(median_I))
    A_em_0 = max(0.01, float(median_I[i_peak]) - 1.0)
    v_ctr_0 = float(vel[i_peak])

    # 半高宽 → Gaussian 宽度
    half_max = 1.0 + A_em_0 * 0.5
    left_hm = np.where(median_I[:i_peak] >= half_max)[0]
    right_hm = np.where(median_I[i_peak:] >= half_max)[0]
    fwhm_pts = (len(left_hm) + len(right_hm))
    dv_mean = float(np.mean(np.abs(np.diff(vel)))) if len(vel) > 1 else 1.0
    fwhm_est = fwhm_pts * dv_mean
    sigma_em_0 = max(5.0, fwhm_est / 2.3548)

    # 自吸收：在线心附近 ±sigma/4 范围内找强度最低点
    half_win = max(1, int(sigma_em_0 * 0.25 / dv_mean))
    iL = max(0, i_peak - half_win)
    iR = min(len(vel) - 1, i_peak + half_win)
    i_dip = iL + int(np.argmin(median_I[iL:iR + 1]))
    # 用纯发射模型在 dip 处估算"理论值"与观测差
    u_dip = (vel[i_dip] - v_ctr_0) / max(sigma_em_0, 0.1)
    em_at_dip = 1.0 + A_em_0 * _voigt_real(np.array([u_dip]), 0.15)[0]
    A_abs_0 = max(0.0, float(em_at_dip) - float(median_I[i_dip]))
    sigma_abs_0 = max(2.0, sigma_em_0 / 5.0)

    _log(f"  初值估算: A_em={A_em_0:.3f}, σ_em={sigma_em_0:.1f} km/s, "
         f"A_abs={A_abs_0:.3f}, σ_abs={sigma_abs_0:.1f} km/s, "
         f"v_center={v_ctr_0:.2f} km/s")

    # --- 非线性拟合（两阶段）---
    # Stage 1: 仅拟合发射成分（固定 A_abs=0），得到稳定的发射参数
    def _emission_only(vel, A_em, sigma_em, a_em, v_center):
        u_em = (vel - v_center) / max(float(sigma_em), 0.1)
        H_em = _voigt_real(u_em, max(0.0, float(a_em)))
        return 1.0 + float(A_em) * H_em

    p0_em = [A_em_0, sigma_em_0, 0.15, v_ctr_0]
    lo_em = [0.0, 1.0, 0.0, v_ctr_0 - 30.0]
    hi_em = [
        max(10.0, A_em_0 * 6.0),
        max(300.0, sigma_em_0 * 6.0), 2.0, v_ctr_0 + 30.0
    ]
    p0_em_safe = [
        float(np.clip(p0_em[k], lo_em[k], hi_em[k])) for k in range(4)
    ]

    A_em_s1, sig_em_s1, a_em_s1, v_ctr_s1 = p0_em_safe  # fallback
    fit_ok = False
    try:
        popt_em, _ = curve_fit(
            _emission_only,
            vel,
            median_I,
            p0=p0_em_safe,
            bounds=(lo_em, hi_em),
            maxfev=20000,
            ftol=1e-10,
            xtol=1e-10,
        )
        A_em_s1, sig_em_s1, a_em_s1, v_ctr_s1 = popt_em
        fit_ok = True
    except Exception as exc:
        _log(f"  [警告] Stage-1 curve_fit 未收敛: {exc}，使用形态估算。")

    # Stage 2: 拟合自吸收成分（扣除 Stage-1 发射后求残差）
    residual = median_I - (
        _emission_only(vel, A_em_s1, sig_em_s1, a_em_s1, v_ctr_s1) - 1.0)

    def _absorption_only(vel, A_abs, sigma_abs, a_abs):
        u_abs = (vel - v_ctr_s1) / max(float(sigma_abs), 0.1)
        H_abs = _voigt_real(u_abs, max(0.0, float(a_abs)))
        return 1.0 - float(A_abs) * H_abs

    # 只在峰值附近寻找吸收（中心窗口）
    win_half = max(1, int(sig_em_s1 * 0.6 / dv_mean))
    i_ctr = int(np.argmin(np.abs(vel - v_ctr_s1)))
    iL2 = max(0, i_ctr - win_half)
    iR2 = min(len(vel) - 1, i_ctr + win_half)
    A_abs_dip = max(0.0, float(1.0 - np.min(residual[iL2:iR2 + 1])))
    sigma_abs_0_s2 = max(2.0, sig_em_s1 / 4.0)

    A_abs_f, sig_abs_f, a_abs_f = A_abs_0, sigma_abs_0_s2, 0.10  # fallback

    if A_abs_dip > 0.01:
        p0_abs = [A_abs_dip, sigma_abs_0_s2, 0.10]
        lo_abs = [0.0, 0.5, 0.0]
        hi_abs = [
            min(A_em_s1, A_abs_dip * 5.0 + 0.5),
            max(80.0, sigma_abs_0_s2 * 4.0), 2.0
        ]
        p0_abs_safe = [
            float(np.clip(p0_abs[k], lo_abs[k], hi_abs[k])) for k in range(3)
        ]
        try:
            popt_abs, _ = curve_fit(
                _absorption_only,
                vel[iL2:iR2 + 1],
                residual[iL2:iR2 + 1],
                p0=p0_abs_safe,
                bounds=(lo_abs, hi_abs),
                maxfev=10000,
            )
            A_abs_f, sig_abs_f, a_abs_f = popt_abs
        except Exception as exc2:
            _log(f"  [警告] Stage-2 curve_fit 未收敛: {exc2}，自吸收强度设为 0。")
            A_abs_f, sig_abs_f, a_abs_f = 0.0, sigma_abs_0_s2, 0.10

    A_em_f = float(A_em_s1)
    sig_em_f = float(sig_em_s1)
    a_em_f = float(a_em_s1)
    A_abs_f = float(A_abs_f)
    sig_abs_f = float(sig_abs_f)
    a_abs_f = float(a_abs_f)
    v_ctr_f = float(v_ctr_s1)

    _log("")
    _log("H-alpha 复合模型自动初始参数：")
    _log(
        f"  发射成分:  A_em={A_em_f:.4f},  σ_G={sig_em_f:.2f} km/s,  a={a_em_f:.4f}"
    )
    _log(
        f"  吸收成分:  A_abs={A_abs_f:.4f}, σ_G={sig_abs_f:.2f} km/s,  a={a_abs_f:.4f}"
    )
    _log(f"  v_center = {v_ctr_f:.3f} km/s  {'（拟合）' if fit_ok else '（估算）'}")
    _log("")

    params = {
        "emission_strength": float(A_em_f),
        "emission_gauss_kms": float(sig_em_f),
        "emission_lorentz_ratio": float(a_em_f),
        "absorption_strength": float(A_abs_f),
        "absorption_gauss_kms": float(sig_abs_f),
        "absorption_lorentz_ratio": float(a_abs_f),
    }

    # --- 构建拟合曲线（高密度）---
    v_dense = np.linspace(float(vel[0]), float(vel[-1]), 600)
    u_dense_em = (v_dense - v_ctr_f) / max(sig_em_f, 0.1)
    u_dense_abs = (v_dense - v_ctr_f) / max(sig_abs_f, 0.1)
    I_em_only = 1.0 + A_em_f * _voigt_real(u_dense_em, max(0.0, a_em_f))
    I_abs_comp = A_abs_f * _voigt_real(u_dense_abs, max(0.0, a_abs_f))
    I_fit = I_em_only - I_abs_comp
    # 1 + (I_em-1) 即为仅发射分量（含连续谱）；吸收作为单独参考线
    I_abs_trace = 1.0 - I_abs_comp  # 仅自吸收分量的"轮廓"（最低点<1）

    # --- Plotly 可视化数据 ---
    traces = []

    # 各历元（灰色细线）
    for i, I_obs in enumerate(indiv_I):
        traces.append({
            "type": "scatter",
            "x": vel.tolist(),
            "y": I_obs.tolist(),
            "mode": "lines",
            "name": f"Obs {i + 1:02d}",
            "line": {
                "color": "rgba(150,150,200,0.40)",
                "width": 1
            },
            "showlegend": False,
        })

    # 中值谱
    traces.append({
        "type": "scatter",
        "x": vel.tolist(),
        "y": median_I.tolist(),
        "mode": "lines",
        "name": "中值谱",
        "line": {
            "color": "#cdd6f4",
            "width": 2.5
        },
    })

    # 发射成分（点线，橙色）
    traces.append({
        "type": "scatter",
        "x": v_dense.tolist(),
        "y": I_em_only.tolist(),
        "mode": "lines",
        "name": f"发射成分 (A={A_em_f:.2f}, σ={sig_em_f:.0f} km/s)",
        "line": {
            "color": "#fab387",
            "width": 1.5,
            "dash": "dot"
        },
    })

    # 自吸收成分（虚线，绿色），仅当 A_abs >阈值时显示
    if A_abs_f > 0.01:
        traces.append({
            "type": "scatter",
            "x": v_dense.tolist(),
            "y": I_abs_trace.tolist(),
            "mode": "lines",
            "name": f"自吸收成分 (A={A_abs_f:.2f}, σ={sig_abs_f:.0f} km/s)",
            "line": {
                "color": "#a6e3a1",
                "width": 1.5,
                "dash": "dash"
            },
        })

    # 拟合复合模型（蓝色实线）
    fit_label = "拟合模型" if fit_ok else "初始估算模型"
    traces.append({
        "type": "scatter",
        "x": v_dense.tolist(),
        "y": I_fit.tolist(),
        "mode": "lines",
        "name": fit_label,
        "line": {
            "color": "#89dceb",
            "width": 2.5
        },
    })

    layout = {
        "title": {
            "text": "Hα 复合模型自动初始参数估算",
            "font": {
                "size": 14,
                "color": "#cdd6f4"
            },
        },
        "xaxis": {
            "title": "速度 (km/s)",
            "zeroline": True,
            "zerolinecolor": "rgba(200,200,200,0.25)",
            "gridcolor": "rgba(200,200,200,0.10)",
            "color": "#cdd6f4",
        },
        "yaxis": {
            "title": "I / I_c",
            "gridcolor": "rgba(200,200,200,0.10)",
            "color": "#cdd6f4",
        },
        "plot_bgcolor": "#1e1e2e",
        "paper_bgcolor": "#1e1e2e",
        "legend": {
            "bgcolor": "rgba(0,0,0,0)",
            "font": {
                "size": 11,
                "color": "#cdd6f4"
            },
        },
        "margin": {
            "l": 60,
            "r": 20,
            "t": 50,
            "b": 50
        },
        "hovermode": "x unified",
    }

    return {
        "params": params,
        "plot_data": {
            "data": traces,
            "layout": layout
        },
        "fit_ok": fit_ok,
    }
