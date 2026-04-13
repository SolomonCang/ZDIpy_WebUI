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
    # 兼容 numpy 1.24+，用 np.real 替代 .real 属性
    return np.real(W).reshape(np.asarray(u).shape)


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
) -> tuple[np.ndarray, np.ndarray, list, np.ndarray]:
    """将所有观测插值到公共速度网格，返回中值谱。

    Returns
    -------
    vel_common : (M,) ndarray
        公共速度格点（km/s，恒星静止系）。
    median_I : (M,) ndarray
        各格点 Stokes I 中值。
    indiv_I : list of (M,) ndarray
        各历元插值后的 Stokes I（归一）。
    median_Isig : (M,) ndarray
        各格点 Stokes I 噪声中值。
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
        indiv_Isig = [obs.specIsig.copy() for obs in obsSet]
        median_I = np.median(np.vstack(indiv_I), axis=0)
        median_Isig = np.median(np.vstack(indiv_Isig), axis=0)
        return vel_common, median_I, indiv_I, median_Isig

    vel_common = np.arange(v_min, v_max + 0.5 * dv, dv)
    indiv_I = []
    indiv_Isig = []
    for i, obs in enumerate(obsSet):
        v_i = rest_vels[i]
        # 边界值采用最近端点
        I_interp = np.interp(vel_common,
                             v_i,
                             obs.specI,
                             left=obs.specI[0],
                             right=obs.specI[-1])
        Isig_interp = np.interp(vel_common,
                                v_i,
                                obs.specIsig,
                                left=obs.specIsig[0],
                                right=obs.specIsig[-1])
        indiv_I.append(I_interp)
        indiv_Isig.append(Isig_interp)

    median_I = np.median(np.vstack(indiv_I), axis=0)
    median_Isig = np.median(np.vstack(indiv_Isig), axis=0)
    return vel_common, median_I, indiv_I, median_Isig


# ===========================================================================
# 模块 2：严格前向模型拟合（盘积分 χ² 最小化）
# ===========================================================================


def _fit_via_forward_model(
    obsSet: list,
    vel_rs: np.ndarray,
    stellar_params: dict,
    p0: list,
    lo: list,
    hi: list,
    log_fn: Optional[Callable[[str], None]] = None,
) -> tuple:
    """通过盘积分前向模型对局域 H-alpha 线参数做最小二乘拟合。

    对每次函数评估，为每个观测历元调用 ``diskIntProfAndDerivHalpha``，
    将旋转展宽内禀地纳入正向模型，然后对各历元观测谱 Stokes I 计算总
    chi-squared 并最小化。

    与直接拟合中值盘积分谱相比，本函数得到的 ``widthGauss_em`` /
    ``widthGauss_abs`` 是真正的局域热/湍动宽度，不包含旋转展宽。

    Parameters
    ----------
    obsSet : list of obsProf
        观测数组（各历元独立）。
    vel_rs : (N_obs,) ndarray
        各历元恒星视向速度（km/s）。
    stellar_params : dict
        必须包含以下键：
        ``vel_eq_kms``, ``inc_rad``, ``period_days``, ``jdates``,
        ``jdate_ref``, ``wl0_nm``, ``lande_g``, ``limb_dark``,
        ``grav_dark``, ``fV``, ``inst_res``。
        可选键：``d_omega`` (默认 0.0)、``n_rings`` (默认 15)。
    p0 : list of 5 floats
        初值 [A_em, sigma_em, a_em, A_abs, sigma_abs]（来自形态估算）。
    lo, hi : list of 5 floats
        参数下界和上界。
    log_fn : callable, optional

    Returns
    -------
    popt : list of 5 floats
        最优局域参数 [A_em, sigma_em, a_em, A_abs, sigma_abs]。
    fit_ok : bool
    vel_model : (M,) ndarray or None
        平均模型轮廓的速度轴（km/s，恒星静止系），用于可视化。
    I_model_mean : (M,) ndarray or None
        各相位前向模型 IIc 的算术平均，用于可视化。
    """
    from scipy.optimize import minimize  # noqa: PLC0415
    from core.line_models.halpha import (  # noqa: PLC0415
        lineDataHalpha as _HalphaData, diskIntProfAndDerivHalpha as _DiskInt,
    )
    from core.geometry.stellar_grid import starGrid  # noqa: PLC0415
    from core.geometry.visibility import BatchVisibleGrid  # noqa: PLC0415
    from core.geometry.differential_rotation import getCyclesClat  # noqa: PLC0415
    from core.brightnessGeom import brightMap as _BrightMap  # noqa: PLC0415

    def _log(msg: str) -> None:
        if log_fn is not None:
            log_fn(msg)

    # --- 解包恒星参数 ---
    vel_eq = float(stellar_params["vel_eq_kms"])
    inc_rad = float(stellar_params["inc_rad"])
    period = float(stellar_params["period_days"])
    d_omega = float(stellar_params.get("d_omega", 0.0))
    jdates = np.asarray(stellar_params["jdates"])
    jdate_ref = float(stellar_params["jdate_ref"])
    n_rings = int(stellar_params.get("n_rings", 15))
    wl0 = float(stellar_params.get("wl0_nm", 656.28))
    lande_g = float(stellar_params.get("lande_g", 1.048))
    limb_dark = float(stellar_params.get("limb_dark", 0.0))
    grav_dark = float(stellar_params.get("grav_dark", 0.0))
    fV = float(stellar_params.get("fV", 1.0))
    inst_res = float(stellar_params.get("inst_res", -1.0))

    _log(f"  [前向模型] 建立格点: nRings={n_rings}, "
         f"v_eq={vel_eq:.1f} km/s, incl={float(np.degrees(inc_rad)):.1f}°, "
         f"P={period:.3f} d")

    # --- 建立恒星表面格点（粗网格，仅用于初始化，不影响主反演）---
    s_grid = starGrid(n_rings, verbose=0)
    n_cells = s_grid.numPoints

    # --- 建立批量可见性（所有观测相位）---
    cycles_at_clat = getCyclesClat(period, d_omega, jdates, jdate_ref,
                                   s_grid.clat)
    batch_vis = BatchVisibleGrid(
        star_grid=s_grid,
        inclination=inc_rad,
        cycles_at_clat=cycles_at_clat,
        period=period,
        dOmega=d_omega,
    )

    # --- 均匀亮度图 + 零磁场（纯 Stokes I 拟合，无需磁场）---
    bri_map = _BrightMap(s_grid.clat, s_grid.long)  # bright = ones
    v_mag_cart = np.zeros((3, n_cells))

    # --- 为各历元构建波长网格（nm，恒星静止系）---
    wl_syn_list = [(obs.wl - float(vel_rs[i])) / _C_KMS * wl0 + wl0
                   for i, obs in enumerate(obsSet)]

    # 预缓存噪声权重（避免在目标函数内重复 clip）
    sig_list = [np.clip(obs.specIsig, 1e-6, None) for obs in obsSet]

    # 原地可修改的 lineDataHalpha 工作对象（避免每次评估构造新对象）
    ldata_work = _HalphaData.from_parameters(
        wavelength_nm=wl0,
        lande_g=lande_g,
        limb_darkening=limb_dark,
        gravity_darkening=grav_dark,
        emission_strength=float(p0[0]),
        emission_gauss_kms=float(p0[1]),
        emission_lorentz_ratio=float(p0[2]),
        absorption_strength=float(p0[3]),
        absorption_gauss_kms=float(p0[4]),
        absorption_lorentz_ratio=0.15,
        filling_factor_V=fV,
        instRes=inst_res,
    )

    def _update_ldata(params: list) -> None:
        ldata_work.A_em[0] = params[0]
        ldata_work.widthGauss_em[0] = params[1]
        ldata_work.widthLorentz_em[0] = params[2]
        ldata_work.A_abs[0] = params[3]
        ldata_work.widthGauss_abs[0] = params[4]
        ldata_work.str[0] = params[0]
        ldata_work.widthGauss[0] = params[1]
        ldata_work.widthLorentz[0] = params[2]

    def _forward_chi2(params_raw: np.ndarray) -> float:
        params = [
            float(np.clip(params_raw[k], lo[k], hi[k])) for k in range(5)
        ]
        _update_ldata(params)
        chi2 = 0.0
        for i_obs, obs in enumerate(obsSet):
            vis = batch_vis[i_obs]
            spec = _DiskInt(vis, v_mag_cart, 0, bri_map, ldata_work, vel_eq,
                            wl_syn_list[i_obs], 0, 0)
            spec.convolveIGnumpy(inst_res)
            diff = obs.specI - spec.IIc
            chi2 += float(np.sum((diff / sig_list[i_obs])**2))
        return chi2

    # --- 运行 L-BFGS-B 优化 ---
    bounds = [(lo[k], hi[k]) for k in range(5)]
    fit_ok = False
    popt = list(p0)
    chi2_p0 = _forward_chi2(p0)
    _log(f"  [前向模型] 初值 χ²={chi2_p0:.2f}，开始 L-BFGS-B 优化 "
         f"({len(obsSet)} 历元 × {n_cells} 格点)…")

    try:
        result = minimize(
            _forward_chi2,
            p0,
            method='L-BFGS-B',
            bounds=bounds,
            options={
                'maxiter': 300,
                'ftol': 1e-10,
                'gtol': 1e-6
            },
        )
        chi2_final = float(result.fun)
        if chi2_final < chi2_p0 * 0.999:
            popt = [
                float(np.clip(result.x[k], lo[k], hi[k])) for k in range(5)
            ]
            fit_ok = True
            _log(f"  [前向模型] 收敛: χ²={chi2_final:.2f}, 迭代={result.nit}")
        else:
            _log(f"  [前向模型] 警告: 未能改善初值 "
                 f"(χ²={chi2_final:.2f} vs {chi2_p0:.2f})，回退形态估算。")
    except Exception as exc:
        _log(f"  [前向模型] 优化异常: {exc}，回退形态估算。")

    # --- 构建平均模型轮廓（可视化用）---
    vel_model: Optional[np.ndarray] = None
    I_model_mean: Optional[np.ndarray] = None
    try:
        _update_ldata(popt)
        rest_vels_m = [
            obs.wl - float(vel_rs[i]) for i, obs in enumerate(obsSet)
        ]
        v_min_m = max(float(np.min(v)) for v in rest_vels_m)
        v_max_m = min(float(np.max(v)) for v in rest_vels_m)
        dv_m = min(
            float(np.mean(np.abs(np.diff(v)))) for v in rest_vels_m
            if len(v) > 1)
        vel_model = np.arange(v_min_m, v_max_m + 0.5 * dv_m, dv_m)
        wl_common_m = vel_model / _C_KMS * wl0 + wl0
        I_models = []
        for i_obs in range(len(obsSet)):
            vis = batch_vis[i_obs]
            spec = _DiskInt(vis, v_mag_cart, 0, bri_map, ldata_work, vel_eq,
                            wl_common_m, 0, 0)
            spec.convolveIGnumpy(inst_res)
            I_models.append(spec.IIc.copy())
        I_model_mean = np.mean(np.vstack(I_models), axis=0)
    except Exception as exc2:
        _log(f"  [前向模型] 可视化轮廓构建失败: {exc2}")

    return popt, fit_ok, vel_model, I_model_mean


# ===========================================================================
# 模块 3：自动参数估算（入口）
# ===========================================================================

# ---------------------------------------------------------------------------
# 主函数
# ---------------------------------------------------------------------------


def auto_estimate_halpha_params(
    obsSet: list,
    vel_rs: np.ndarray,
    current_lineData,
    stellar_params: Optional[dict] = None,
    log_fn: Optional[Callable[[str], None]] = None,
) -> dict:
    """估算 H-alpha 双 Voigt 初始参数。

    算法步骤：

    1. 对所有历元 Stokes I 在恒星静止系插值至公共速度格点，取中值。
    2. 用简单形态法估算初值（峰高、半宽、中心速度、吸收深度）。
    3a. 若提供 ``stellar_params``（严格方案）：调用
        ``_fit_via_forward_model``，以恒星几何正向建模盘积分轮廓，
        最小化各历元 chi-squared，得到真正的局域热/湍动线宽。
    3b. 否则（旧方案）：用 ``curve_fit`` 直接拟合中值盘积分谱，
        所得线宽包含旋转展宽，对快速自转星存在系统偏差。
    4. 构建 Plotly 格式绘图数据并返回。

    Parameters
    ----------
    obsSet : list of obsProf
        （已归一化的）观测数组。
    vel_rs : (N_obs,) ndarray
        各历元恒星视向速度（km/s）。
    current_lineData : lineDataHalpha or None
        现有线参数（当前未使用，保留接口兼容性）。
    stellar_params : dict, optional
        恒星几何参数字典，启用严格前向模型拟合。
        必须包含键：
        ``vel_eq_kms``, ``inc_rad``, ``period_days``, ``jdates``,
        ``jdate_ref``, ``wl0_nm``, ``lande_g``, ``limb_dark``,
        ``grav_dark``, ``fV``, ``inst_res``。
        可选键：``d_omega`` (默认 0.0)、``n_rings`` (默认 15)。
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
    vel, median_I, indiv_I, median_Isig = _build_median_spectrum(
        obsSet, vel_rs)
    _log(f"  速度范围: [{float(vel[0]):.1f}, {float(vel[-1]):.1f}] km/s，"
         f"格点数: {len(vel)}")

    # --- 初值估算 ---
    i_peak = int(np.argmax(median_I))
    A_em_0 = max(0.01, float(median_I[i_peak]) - 1.0)
    dv_mean = float(np.mean(np.abs(np.diff(vel)))) if len(vel) > 1 else 1.0

    # 发射线中心：用超出连续谱部分(I-1)的加权质心，对自反转/双峰轮廓更稳健
    em_excess = np.clip(median_I - 1.0, 0.0, None)
    total_weight = float(np.sum(em_excess))
    if total_weight > 1e-9:
        v_ctr_0 = float(np.sum(vel * em_excess) / total_weight)
    else:
        v_ctr_0 = float(vel[i_peak])
    i_ctr_0 = int(np.argmin(np.abs(vel - v_ctr_0)))

    # v_center 预先由质心确定后固定（不作为拟合自由参数）
    _log(f"  固定 v_center = {v_ctr_0:.2f} km/s（发射加权质心）")

    # 半高宽 → Gaussian 宽度
    half_max = 1.0 + A_em_0 * 0.5
    above_hm = np.where(median_I >= half_max)[0]
    if len(above_hm) >= 2:
        fwhm_est = (int(above_hm[-1]) - int(above_hm[0])) * dv_mean
    else:
        fwhm_est = 10.0 * dv_mean
    sigma_em_0 = max(5.0, fwhm_est / 2.3548)

    # --- 初值修正：用观测峰到质心的距离对 A_em 做高斯补偿 ---
    # 在质心为中心的Voigt中，观测峰(可能在v_ctr左右)并非峰顶，需要向上修正
    v_peak_offset = abs(float(vel[i_peak]) - v_ctr_0)
    H_at_peak_offset = float(
        np.exp(-0.5 * (v_peak_offset / max(sigma_em_0, 1.0))**2))
    if H_at_peak_offset > 0.1:
        A_em_0 = max(A_em_0, float(median_I[i_peak] - 1.0) / H_at_peak_offset)

    # 吸收初值：用质心处的理论发射值与观测值之差估算
    I_em_at_ctr = 1.0 + A_em_0  # 如无自吸收时质心处的理论强度
    I_data_at_ctr = float(median_I[i_ctr_0])
    A_abs_0 = max(0.0, I_em_at_ctr - I_data_at_ctr)
    sigma_abs_0 = max(2.0, sigma_em_0 / 2.0)

    _log(f"  初值估算: A_em={A_em_0:.3f}, σ_em={sigma_em_0:.1f} km/s, "
         f"A_abs={A_abs_0:.3f}, σ_abs={sigma_abs_0:.1f} km/s, "
         f"v_center={v_ctr_0:.2f} km/s")

    # --- 非线性拟合（复合模型，v_center 固定）---
    # 直接拟合 emission - absorption 复合模型；a_em 作为自由参数以适配翼部形状
    def _composite_model(vel, A_em, sigma_em, a_em, A_abs, sigma_abs):
        u_em = (vel - v_ctr_0) / max(float(sigma_em), 0.1)
        u_abs = (vel - v_ctr_0) / max(float(sigma_abs), 0.1)
        H_em = _voigt_real(u_em, max(0.0, float(a_em)))
        H_abs = _voigt_real(u_abs, 0.15)
        return 1.0 + float(A_em) * H_em - float(A_abs) * H_abs

    p0_c = [A_em_0, sigma_em_0, 0.15, A_abs_0, sigma_abs_0]
    lo_c = [0.0, sigma_em_0 * 0.3, 0.0, 0.0, 1.0]
    hi_c = [
        max(5.0, A_em_0 * 2.5),  # A_em: 不超过初值2.5倍
        sigma_em_0 * 2.5,  # sigma_em: 不超过初值2.5倍
        1.5,  # a_em: Lorentz比例（0=纯Gaussian，1.5接近Lorentzian翼）
        max(A_abs_0 * 3.0, 1.0),  # A_abs: 不超过A_abs_0的3倍
        sigma_em_0 * 0.9,  # sigma_abs < sigma_em（自吸收更窄）
    ]
    p0_c_safe = [float(np.clip(p0_c[k], lo_c[k], hi_c[k])) for k in range(5)]

    # 使用形态估算作为 fallback
    A_em_f, sig_em_f = p0_c_safe[0], p0_c_safe[1]
    a_em_f = 0.15
    A_abs_f, sig_abs_f = p0_c_safe[3], p0_c_safe[4]
    a_abs_f = 0.15
    fit_ok = False
    _vel_model: Optional[np.ndarray] = None
    _I_model_mean: Optional[np.ndarray] = None

    if stellar_params is not None:
        # --- 严格方案：盘积分前向模型拟合 ---
        _log("  使用严格前向模型（盘积分）估算局域线参数…")
        try:
            _popt, fit_ok, _vel_model, _I_model_mean = _fit_via_forward_model(
                obsSet, vel_rs, stellar_params, p0_c_safe, lo_c, hi_c, log_fn)
            A_em_f, sig_em_f, a_em_f, A_abs_f, sig_abs_f = _popt
        except Exception as exc:
            _log(f"  [前向模型] 意外错误: {exc}，回退形态估算。")
    else:
        # --- 旧方案：直接拟合中值盘积分谱（忽略旋转展宽）---
        try:
            # 权重：发射峰越高权重越大；宽松的中心衰减（2.5σ）保留翼部数据权重
            sigma_noise = np.clip(median_Isig, 1e-6, None)
            w_peak = np.clip(median_I - 1.0, 0.05, None)
            w_center = np.exp(-0.5 * ((vel - v_ctr_0) / (sigma_em_0 * 2.5))**2)
            sigma_eff = sigma_noise / np.clip(w_peak * w_center, 1e-4, None)

            # 拟合窗口：质心±3.5σ（覆盖发射翼部，使优化器拟合陡翼）
            fit_mask = np.abs(vel - v_ctr_0) <= sigma_em_0 * 3.5
            popt, _ = curve_fit(
                _composite_model,
                vel[fit_mask],
                median_I[fit_mask],
                p0=p0_c_safe,
                bounds=(lo_c, hi_c),
                sigma=sigma_eff[fit_mask],
                absolute_sigma=True,
                maxfev=30000,
                ftol=1e-9,
                xtol=1e-9,
            )
            A_em_f, sig_em_f, a_em_f, A_abs_f, sig_abs_f = popt

            # 合理性检验：发射宽度不超过初值2.5倍、吸收强度不超过发射的75%
            if sig_em_f > sigma_em_0 * 2.5 or A_abs_f > A_em_f * 0.75:
                _log(f"  [调整] 拟合结果异常(σ={sig_em_f:.1f}/{sigma_em_0:.1f})，"
                     f"回退到形态估算。")
                A_em_f, sig_em_f = p0_c_safe[0], p0_c_safe[1]
                a_em_f = 0.15
                A_abs_f, sig_abs_f = p0_c_safe[3], p0_c_safe[4]
            else:
                fit_ok = True
        except Exception as exc:
            _log(f"  [警告] curve_fit 未收敛: {exc}，使用形态估算。")

    v_ctr_f = v_ctr_0

    _method = "前向模型（盘积分）" if stellar_params is not None else "直接拟合中值谱"
    _log("")
    _log(f"H-alpha 复合模型自动初始参数 [{_method}]：")
    _log(
        f"  发射成分:  A_em={A_em_f:.4f},  σ_G={sig_em_f:.2f} km/s,  a_em={a_em_f:.3f}"
    )
    _log(f"  吸收成分:  A_abs={A_abs_f:.4f}, σ_G={sig_abs_f:.2f} km/s")
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
    # --- 构建拟合曲线（高密度速度网格）---
    v_dense = np.linspace(float(vel[0]), float(vel[-1]), 600)
    u_dense_em = (v_dense - v_ctr_f) / max(sig_em_f, 0.1)
    u_dense_abs = (v_dense - v_ctr_f) / max(sig_abs_f, 0.1)
    # 局域 Voigt 分量轮廓（不含旋转展宽，用于组分示意）
    I_em_only = 1.0 + A_em_f * _voigt_real(u_dense_em, max(0.0, a_em_f))
    I_abs_comp = A_abs_f * _voigt_real(u_dense_abs, max(0.0, a_abs_f))
    I_fit_local = I_em_only - I_abs_comp  # 局域复合轮廓
    I_abs_trace = 1.0 - I_abs_comp  # 仅自吸收分量的"轮廓"（最低点<1）

    # 决定"拟合模型"曲线来源：
    # - 严格方案：使用前向模型平均 IIc（包含旋转展宽）
    # - 旧方案：使用局域复合 Voigt（不含旋转展宽）
    if _vel_model is not None and _I_model_mean is not None:
        fit_x = _vel_model.tolist()
        fit_y = _I_model_mean.tolist()
        fit_label = "前向模型均值（盘积分）" if fit_ok else "前向模型均值（估算）"
    else:
        fit_x = v_dense.tolist()
        fit_y = I_fit_local.tolist()
        fit_label = "拟合模型" if fit_ok else "初始估算模型"

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

    # 局域发射成分（点线，橙色）
    em_label = f"局域发射分量 (A={A_em_f:.2f}, σ_local={sig_em_f:.0f} km/s)"
    traces.append({
        "type": "scatter",
        "x": v_dense.tolist(),
        "y": I_em_only.tolist(),
        "mode": "lines",
        "name": em_label,
        "line": {
            "color": "#fab387",
            "width": 1.5,
            "dash": "dot"
        },
    })

    # 局域自吸收成分（虚线，绿色），仅当 A_abs > 阈值时显示
    if A_abs_f > 0.01:
        abs_label = f"局域自吸收分量 (A={A_abs_f:.2f}, σ_local={sig_abs_f:.0f} km/s)"
        traces.append({
            "type": "scatter",
            "x": v_dense.tolist(),
            "y": I_abs_trace.tolist(),
            "mode": "lines",
            "name": abs_label,
            "line": {
                "color": "#a6e3a1",
                "width": 1.5,
                "dash": "dash"
            },
        })

    # 拟合/前向模型曲线（蓝色实线）
    traces.append({
        "type": "scatter",
        "x": fit_x,
        "y": fit_y,
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
        "v_center": float(v_ctr_f),
        "plot_data": {
            "data": traces,
            "layout": layout
        },
        "fit_ok": fit_ok,
    }
