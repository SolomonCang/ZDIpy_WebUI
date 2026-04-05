"""core.geometry.differential_rotation — 差分自转速率辅助函数。

提供太阳型差分自转律（赤道-极差）的各类计算：
- calcOmegaClat : 各余纬对应的角速度
- getCyclesClat : 各余纬对应的旋转圈数（相位），形状 (N_phases, N_cells)
- calcVelDiffrotFactor : 归一化速度比例因子 Omega(clat)/Omega_eq

这些函数不依赖其他 geometry 子模块，可被 stellar_grid.py 和
visibility.py 安全导入。
"""
from __future__ import annotations

import numpy as np


def calcOmegaClat(period: float, dOmega: float,
                  clat: np.ndarray) -> np.ndarray:
    """按余纬 clat 计算角速度 Omega(clat)。

    采用太阳型差分自转律::

        Omega(lat)  = Omega_eq - dOmega * sin^2(lat)
        Omega(clat) = Omega_eq - dOmega * cos^2(clat)

    Parameters
    ----------
    period : float
        恒星赤道自转周期（天）。
    dOmega : float
        赤道-极角速度差 (rad/day)。
    clat : (N_cells,) ndarray
        表面元余纬坐标（弧度）。

    Returns
    -------
    omega_clat : (N_cells,) ndarray
        各余纬处的角速度 (rad/day)。
    """
    omega_eq = 2.0 * np.pi / period
    omega_clat = omega_eq - dOmega * np.cos(clat)**2
    return omega_clat


def getCyclesClat(period: float, dOmega: float, jDates: np.ndarray,
                  jDateRef: float, clat: np.ndarray) -> np.ndarray:
    """计算各余纬处对应各观测时刻的旋转圈数（相位）。

    Parameters
    ----------
    period : float
        赤道自转周期（天）。
    dOmega : float
        赤道-极角速度差 (rad/day)。
    jDates : (N_phases,) ndarray
        观测儒略日列表。
    jDateRef : float
        参考儒略日（相位零点）。
    clat : (N_cells,) ndarray
        表面元余纬坐标（弧度）。

    Returns
    -------
    cycleAtClat : (N_phases, N_cells) ndarray
        各相位、各余纬处的旋转圈数。
    """
    omega_clat = calcOmegaClat(period, dOmega, clat)
    period_clat = 2.0 * np.pi / omega_clat
    cycleAtClat = (jDates[:, np.newaxis] -
                   jDateRef) / period_clat[np.newaxis, :]
    return cycleAtClat


def calcVelDiffrotFactor(period: float, dOmega: float,
                         clat: np.ndarray) -> np.ndarray:
    """返回各余纬相对赤道旋转速度的归一化比例因子。

    说明：factor(clat) = Omega(clat) / Omega_eq，用于将赤道线速度
    （由 vsini 估计）推广到任意余纬。

    Parameters
    ----------
    period : float
        赤道自转周期（天）。
    dOmega : float
        赤道-极角速度差 (rad/day)。
    clat : (N_cells,) ndarray
        表面元余纬坐标（弧度）。

    Returns
    -------
    norm_clat : (N_cells,) ndarray
        归一化速度比例因子。
    """
    omega_eq = 2.0 * np.pi / period
    omega_clat = calcOmegaClat(period, dOmega, clat)
    return omega_clat / omega_eq
