"""core.line_models.line_utils — 谱线轮廓辅助函数。

提供轮廓拟合和等效宽度计算中的独立无状态工具函数：
- limbDarkening      : 线性临边昏暗
- calcSynEW          : 合成 Stokes I 轮廓等效宽度（numpy 向量化）
- equivWidComp2      : 等效宽度残差目标函数（scipy.optimize 兼容）
- fitLineStrength    : 通过等效宽度匹配拟合谱线强度

这些函数不依赖 ``profile.py`` 中的任何大型类，可被 ``profile.py`` 安全导入。
"""
from __future__ import annotations

import numpy as np


def limbDarkening(coeffEta: float, angle: np.ndarray) -> np.ndarray:
    """线性临边昏暗定律（Gray 2005, eq. 17.11）。

    .. math::
        I(\\mu) = 1 - \\varepsilon + \\varepsilon \\cos(\\theta)

    Parameters
    ----------
    coeffEta : float or array-like
        临边昏暗系数 ε（可以是 lineData.limbDark[0]）。
    angle : ndarray
        视线与格点法向量夹角（弧度），支持任意形状。

    Returns
    -------
    limbD : ndarray, same shape as *angle*
        临边昏暗权重（无量纲）。
    """
    return 1.0 - coeffEta + coeffEta * np.cos(angle)


def calcSynEW(spec) -> float:
    """计算合成 Stokes I 轮廓的等效宽度（向量化梯形积分）。

    Parameters
    ----------
    spec
        具有 ``wl`` (N_vel,) 和 ``IIc`` (N_vel,) 属性的谱线对象
        （如 ``diskIntProfAndDeriv``）。

    Returns
    -------
    float
        等效宽度（与 ``wl`` 同单位，通常为 nm）。
    """
    dw = spec.wl[1:] - spec.wl[:-1]
    integrand = (1.0 - spec.IIc[:-1]) * dw
    equivWidth = float(np.sum(integrand))
    # 末尾点用最后一个区间宽度近似
    equivWidth += (1.0 - spec.IIc[-1]) * (spec.wl[-1] - spec.wl[-2])
    return equivWidth


def equivWidComp2(lineStr: float, meanEquivWidObs: float, setSynSpec: list,
                  lineData) -> float:
    """返回给定线强下模型等效宽度与观测值的绝对差，供标量最优化使用。

    Parameters
    ----------
    lineStr : float
        待试验的谱线强度（深度）缩放系数。
    meanEquivWidObs : float
        观测平均等效宽度（比较基准；传 0 得到模型 EW 本身）。
    setSynSpec : list
        ``diskIntProfAndDeriv`` 对象列表（每个观测相位一个）。
    lineData
        ``lineData`` 对象，提供 ``str[0]`` 基准线强。

    Returns
    -------
    float
        |EW_model - EW_obs|。
    """
    lineStr0 = lineData.str[0]
    meanEW = 0.0
    nObs = 0
    for spec in setSynSpec:
        scaleI = 1.0 - (1.0 - spec.IIc) * (lineStr / lineStr0)
        dw = spec.wl[1:] - spec.wl[:-1]
        ew = float(np.sum((1.0 - scaleI[:-1]) * dw))
        ew += (1.0 - scaleI[-1]) * (spec.wl[-1] - spec.wl[-2])
        meanEW += ew
        nObs += 1
    meanEW /= float(nObs)
    return float(np.abs(meanEquivWidObs - meanEW))


def fitLineStrength(meanEquivWidObs: float,
                    par,
                    listGridView,
                    vecMagCart: np.ndarray,
                    dMagCart0,
                    briMap,
                    lineData,
                    wlSynSet: list,
                    verbose: int = 1) -> None:
    """通过使模型等效宽度匹配观测值，原地更新 ``lineData.str[0]``。

    .. warning::
        本函数仅适用于 Voigt 吸收线模型（``lineData`` 为 ``lineData`` 或
        ``lineDataUnno`` 类型）。对于 H-alpha 复合发射线模型
        (``lineDataHalpha``/``halpha_compound``)，EW 为负（发射），
        内部 Voigt 吸收模型无法匹配负 EW，会错误地将线强推至近零，
        导致 Stokes V 振幅被彻底压制。
        请在调用方（``ZDIPipeline``）跳过对该模型类型的调用。
    """
    from scipy.optimize import minimize_scalar  # noqa: PLC0415

    # 根据 lineData 实际类型选择正确的盘积分类，避免 UR 模式下用 Voigt EW 校准
    from core.line_models.unno import lineDataUnno as _lineDataUnno  # noqa: PLC0415
    if isinstance(lineData, _lineDataUnno):
        from core.line_models.unno import diskIntProfAndDerivUnno as _DiskIntCls  # noqa: PLC0415
    else:
        from core.line_models.profile import diskIntProfAndDeriv as _DiskIntCls  # noqa: PLC0415

    if verbose == 1:
        print('fitting line strength (by equivalent width)')

    kL0 = lineData.str[0]

    if isinstance(lineData, _lineDataUnno):
        # UR 模型中 EW 与 kL 的关系是非线性的（存在饱和效应），
        # equivWidComp2 中的线性缩放假设仅对 Voigt 弱场模型有效。
        # 此处对每个试验 kL 值重新计算完整盘积分以得到准确 EW。
        n_obs = len(par.cycleList)

        def _ur_ew_residual(kL_trial: float) -> float:
            lineData.str[0] = float(kL_trial)
            total_ew = 0.0
            for nObs in range(n_obs):
                spec = _DiskIntCls(listGridView[nObs], vecMagCart, dMagCart0,
                                   briMap, lineData, par.velEq, wlSynSet[nObs],
                                   0, 0)
                spec.convolveIGnumpy(par.instrumentRes)
                total_ew += calcSynEW(spec)
            return float(np.abs(meanEquivWidObs - total_ew / n_obs))

        fitResult = minimize_scalar(
            _ur_ew_residual,
            bounds=(kL0 * 1e-4, kL0 * 1e4),
            method='bounded',
            options={'xatol': 1e-5},
        )
        lineData.str[0] = fitResult.x

        if verbose == 1:
            # 重新计算以报告最终 EW
            lineData.str[0] = fitResult.x
            final_ew = meanEquivWidObs - _ur_ew_residual(
                fitResult.x) + meanEquivWidObs
            # _ur_ew_residual 返回 |EW_model - EW_obs|，EW_model = meanEquivWidObs ± residual
            # 重新算一次以便打印真实 EW
            lineData.str[0] = fitResult.x
            total_ew = 0.0
            for nObs in range(n_obs):
                spec = _DiskIntCls(listGridView[nObs], vecMagCart, dMagCart0,
                                   briMap, lineData, par.velEq, wlSynSet[nObs],
                                   0, 0)
                spec.convolveIGnumpy(par.instrumentRes)
                total_ew += calcSynEW(spec)
            print('best match line strength is {:f} (ew {:f})'.format(
                fitResult.x, total_ew / n_obs))
    else:
        # Voigt 模型：深度与 kL 线性相关，equivWidComp2 的线性缩放有效
        setSynSpec = []
        for nObs, _phase in enumerate(par.cycleList):
            spec = _DiskIntCls(listGridView[nObs], vecMagCart, dMagCart0,
                               briMap, lineData, par.velEq, wlSynSet[nObs], 0,
                               0)
            spec.convolveIGnumpy(par.instrumentRes)
            setSynSpec.append(spec)

        # 使用有界搜索，搜索区间 [kL*1e-4, kL*1e4]
        fitResult = minimize_scalar(
            equivWidComp2,
            bounds=(kL0 * 1e-4, kL0 * 1e4),
            method='bounded',
            options={'xatol': 1e-5},
            args=(meanEquivWidObs, setSynSpec, lineData),
        )

        if verbose == 1:
            ew = equivWidComp2(fitResult.x, 0, setSynSpec, lineData)
            print('best match line strength is {:f} (ew {:f})'.format(
                fitResult.x, ew))

        lineData.str[0] = fitResult.x
