"""core.line_models.halpha — H-alpha 双峰发射线复合 Voigt 模型。

基于弱场近似（Weak-Field Approximation），将 H-alpha 轮廓建模为
**宽发射成分 + 窄自吸收成分**的双 Voigt 叠加：

    I(λ) = 1 + A_em · H_em(λ) − A_abs · H_abs(λ)
    V(λ) = −Δλ_Z · B_LOS · [A_em · dH_em/dλ − A_abs · dH_abs/dλ]

其中 Voigt 函数导数由 Humlicek w4 的虚部（Faraday-Voigt 函数）解析给出：
    dH/dλ = (2A·norm/Δλ_D) · (−u·H_raw + a·F_raw)
    u = (λ₀ − λ) / Δλ_D,  a = widthLorentz_ratio

对外接口与 ``core.line_models.unno`` 对等，可直接在 pipeline 中替换：

- ``lineDataHalpha``              : H-alpha 线参数容器
- ``localProfileAndDerivHalpha``  : 单格点预计算轮廓及导数核
- ``diskIntProfAndDerivHalpha``   : 盘积分轮廓及全部一阶导数
- ``getAllProfDirivHalpha``        : 批量初始化辅助函数

设计参考：``docs/future/halpha_compound_model.md``
"""
from __future__ import annotations

import numpy as np

from core.line_models.line_utils import limbDarkening
from core.line_models.unno import _voigt_faraday_humlicek

_C_KMS: float = 2.99792458e5  # km/s
_ZEEMAN_CONST: float = 4.66864e-12  # nm/G (e/4πm_e c²)

# ---------------------------------------------------------------------------
# 单分量 Voigt 贡献及导数核的计算辅助
# ---------------------------------------------------------------------------


def _voigt_component(
    wl0: float,
    width_gauss_kms: float,
    width_lorentz: float,
    amplitude: float,
    wls: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """计算单个 Voigt 分量的轮廓值和波长导数。

    与 ``core.line_models.profile`` 完全相同的符号约定：
    ``u = (wl0 − wl) / Δλ_D``（正值对应蓝翼），``w4.real`` 在线心处≈1
    （峰值归一化，不除以 sqrt(π)·Δλ_D），振幅参数 ``amplitude`` 直接
    控制峰高（连续谱归一单位）。

    Voigt 导数公式（与 profile.py 相同）：
        d(amplitude · w4.real)/dλ = (2·amplitude/Δλ_D) · (−u·W.real + a·W.imag)

    Parameters
    ----------
    wl0 : float
        谱线线心波长（nm）。
    width_gauss_kms : float
        Gaussian 多普勒宽度（km/s）。
    width_lorentz : float
        Lorentzian/Gaussian 宽度比（无量纲）。
    amplitude : float
        分量强度（A_em 或 A_abs）。A_em 即峰值高于连续谱的分数；
        A_abs 即中心自吸收深度，均为连续谱归一单位。
    wls : (N_vel, N_cells) ndarray
        各格点已 Doppler 移位的波长网格（nm）。

    Returns
    -------
    comp_I : (N_vel, N_cells)
        amplitude × w4.real(u, a)  — 轮廓贡献
    dcomp_dlambda : (N_vel, N_cells)
        d(comp_I)/dλ  — Stokes V 弱场核的贡献项
    """
    delta_lambda_D = width_gauss_kms / _C_KMS * wl0  # Gaussian 宽度（nm）

    # profile.py 符号约定：u = (wl0 - wl) / Δλ_D
    u = (wl0 - wls) / delta_lambda_D  # (N_vel, N_cells)

    W = _voigt_faraday_humlicek(u, width_lorentz)  # (N_vel, N_cells) complex

    # 峰值归一化（w4.real 在 u=0, a≈0 处 ≈1），与 profile.py 相同约定
    comp_I = amplitude * W.real  # (N_vel, N_cells)

    # d(amplitude · w4.real)/dλ = (2·amplitude/Δλ_D) · (u·W.real − a·W.imag)
    # 推导（u = (wl0-wl)/Δλ_D，du/dλ = -1/Δλ_D，w4.real 导数经 profile.py 验证）：
    #   d(w4.real)/dλ = d(w4.real)/du · (−1/Δλ_D) = 2(−u·H + a·F) · (−1/Δλ_D)
    #                 = (2/Δλ_D) · (u·H − a·F)
    dcomp_dlambda = (amplitude * 2.0 / delta_lambda_D) * (
        u * W.real - width_lorentz * W.imag)  # (N_vel, N_cells)

    return comp_I, dcomp_dlambda


# ---------------------------------------------------------------------------
# lineDataHalpha
# ---------------------------------------------------------------------------


class lineDataHalpha:
    """H-alpha 双峰复合 Voigt 线参数容器。

    字段设计兼容 ``lineData`` / ``lineDataUnno`` 接口：

    - ``wl0``          : 线心波长（nm），通常 656.28
    - ``g``            : 有效 Landé 因子，通常 1.048
    - ``limbDark``     : 临边昏暗系数（色球层建议取 0.0）
    - ``gravDark``     : 重力昏暗系数（通常 0.0）
    - ``str``          : 等同 A_em（保持 interface 兼容）
    - ``widthGauss``   : 等同 widthGauss_em（保持 interface 兼容）
    - ``widthLorentz`` : 等同 widthLorentz_em（保持 interface 兼容）

    H-alpha 专有字段：

    - ``A_em``              : 发射峰强度（> 0，连续谱归一单位）
    - ``widthGauss_em``     : 发射成分 Gaussian 半宽（km/s）
    - ``widthLorentz_em``   : 发射成分 Lorentz/Gauss 比
    - ``A_abs``             : 自吸收强度（≥ 0；= 0 则退化为单发射峰）
    - ``widthGauss_abs``    : 自吸收成分 Gaussian 半宽（km/s）
    - ``widthLorentz_abs``  : 自吸收成分 Lorentz/Gauss 比
    - ``fV``                : 整体 Stokes V 填充因子（默认 1.0）
    """

    def __init__(self) -> None:
        self.wl0 = np.array([])
        self.g = np.array([])
        self.limbDark = np.array([])
        self.gravDark = np.array([])
        # 发射成分
        self.A_em = np.array([])
        self.widthGauss_em = np.array([])
        self.widthLorentz_em = np.array([])
        # 自吸收成分
        self.A_abs = np.array([])
        self.widthGauss_abs = np.array([])
        self.widthLorentz_abs = np.array([])
        # 填充因子
        self.fV = np.array([])
        # 接口兼容字段（alias to emission component）
        self.str = np.array([])
        self.widthGauss = np.array([])
        self.widthLorentz = np.array([])
        self.numLines = 0
        self.instRes: float = -1.0

    @classmethod
    def from_parameters(
        cls,
        wavelength_nm: float,
        lande_g: float,
        limb_darkening: float,
        gravity_darkening: float,
        emission_strength: float,
        emission_gauss_kms: float,
        emission_lorentz_ratio: float,
        absorption_strength: float,
        absorption_gauss_kms: float,
        absorption_lorentz_ratio: float,
        filling_factor_V: float = 1.0,
        instRes: float = -1.0,
    ) -> "lineDataHalpha":
        """由 config.json 数值字段构造 ``lineDataHalpha``。

        Parameters
        ----------
        wavelength_nm : float
            线心波长，H-alpha = 656.28 nm。
        emission_strength : float
            发射峰强度 A_em（连续谱归一，应 > 0）。
        emission_gauss_kms : float
            发射成分 Gaussian 半宽（km/s）。
        emission_lorentz_ratio : float
            发射成分 Lorentz/Gauss 宽度比。
        absorption_strength : float
            自吸收强度 A_abs（≥ 0；置 0 则为无自吸收的单峰发射）。
        absorption_gauss_kms : float
            自吸收成分 Gaussian 半宽（km/s，通常 < emission_gauss_kms）。
        absorption_lorentz_ratio : float
            自吸收成分 Lorentz/Gauss 宽度比。
        filling_factor_V : float
            Stokes V 整体填充因子，默认 1.0。
        instRes : float
            仪器分辨率 R = λ/Δλ；< 0 表示不施加额外仪器展宽。
        """
        obj = cls()
        obj.instRes = float(instRes)
        obj.numLines = 1

        obj.wl0 = np.array([float(wavelength_nm)])
        obj.g = np.array([float(lande_g)])
        obj.limbDark = np.array([float(limb_darkening)])
        obj.gravDark = np.array([float(gravity_darkening)])

        A_em_v = float(emission_strength)
        wG_em_v = float(emission_gauss_kms)
        aL_em_v = float(emission_lorentz_ratio)
        A_abs_v = max(0.0, float(absorption_strength))
        wG_abs_v = float(absorption_gauss_kms)
        aL_abs_v = float(absorption_lorentz_ratio)

        obj.A_em = np.array([A_em_v])
        obj.widthGauss_em = np.array([wG_em_v])
        obj.widthLorentz_em = np.array([aL_em_v])
        obj.A_abs = np.array([A_abs_v])
        obj.widthGauss_abs = np.array([wG_abs_v])
        obj.widthLorentz_abs = np.array([aL_abs_v])
        obj.fV = np.array([max(1e-6, float(filling_factor_V))])

        # 接口兼容（指向发射成分参数）
        obj.str = np.array([A_em_v])
        obj.widthGauss = np.array([wG_em_v])
        obj.widthLorentz = np.array([aL_em_v])

        return obj


# ---------------------------------------------------------------------------
# localProfileAndDerivHalpha
# ---------------------------------------------------------------------------


class localProfileAndDerivHalpha:
    """预计算 H-alpha 复合 Voigt 轮廓及 Stokes V 核。

    与 ``localProfileAndDeriv`` 的关键区别：

    1. 轮廓形状为 ``1 + A_em·H_em − A_abs·H_abs``（发射 + 自吸收）。
    2. Stokes V 核（``VsigKernel``）为两分量导数的线性叠加：
       ``gCoeff · (dI_em/dλ − dI_abs/dλ)``。
    3. 模型关于 Blos 线性（弱场近似），无需数值导数，响应矩阵可预计算。

    Parameters
    ----------
    ldata : lineDataHalpha
    numPts : int
        波长格点数（保留兼容性）。
    wl_grid : (N_vel, N_cells) ndarray
        各格点已 Doppler 移位的波长网格（nm）。
    """

    def __init__(
        self,
        ldata: lineDataHalpha,
        numPts: int,
        wl_grid: np.ndarray,
    ) -> None:
        self.wl = wl_grid  # (N_vel, N_cells)

        wl0 = float(ldata.wl0[0])
        g = float(ldata.g[0])
        fV = float(ldata.fV[0])

        # --- 预计算两个 Voigt 分量的轮廓及导数 ---
        I_em, dI_em = _voigt_component(
            wl0,
            float(ldata.widthGauss_em[0]),
            float(ldata.widthLorentz_em[0]),
            float(ldata.A_em[0]),
            wl_grid,
        )
        I_abs, dI_abs = _voigt_component(
            wl0,
            float(ldata.widthGauss_abs[0]),
            float(ldata.widthLorentz_abs[0]),
            float(ldata.A_abs[0]),
            wl_grid,
        )

        # I_local = 1 + A_em*H_em − A_abs*H_abs
        self.Iunscaled = 1.0 + I_em - I_abs  # (N_vel, N_cells)

        # Stokes V 核：gCoeff × dI/dλ  (gCoeff = −Δλ_Z < 0)
        # 此核在 B 和亮度不变时为常数，可在 __init__ 处完全预计算。
        gCoeff = -_ZEEMAN_CONST * wl0**2 * g  # nm/G，负数
        self.VsigKernel = fV * gCoeff * (dI_em - dI_abs)  # (N_vel, N_cells)

    def updateProfDeriv(
        self,
        ldata: lineDataHalpha,
        Blos: np.ndarray,
        dBlos_d,
        brightMap,
        numPts: int,
        wl_grid: np.ndarray,
        surfaceScale: np.ndarray,
        calcDI: int,
        calcDV: int,
    ) -> None:
        """计算盘积分轮廓及全部一阶导数。

        接口与 ``localProfileAndDeriv.updateProfDeriv`` 完全对等。

        Parameters
        ----------
        Blos : (N_cells,)
            视线磁场分量（G）。
        dBlos_d : (n_types, nTot, N_cells) 或 0
            Blos 对球谐系数的导数（``calcDV=0`` 时传 0）。
        brightMap : object
            亮度图对象（``bright`` 属性为 (N_cells,)）。
        surfaceScale : (N_cells,)
            limb darkening × gravity darkening × projArea × visible。
        calcDI, calcDV : int
            亮度/磁场导数开关（0=不计算，1=计算）。
        """
        contin = surfaceScale * brightMap.bright  # (N_cells,)

        self.Ic = np.sum(contin)
        invSumIc0 = 1.0 / self.Ic

        # V = Blos × contin × VsigKernel  （线性于 Blos）
        V_cell = (Blos *
                  contin)[np.newaxis, :] * self.VsigKernel  # (N_vel, N_cells)

        # 盘积分并归一化
        self.I = np.sum(contin[np.newaxis, :] * self.Iunscaled,
                        axis=1) * invSumIc0
        self.V = np.sum(V_cell, axis=1) * invSumIc0

        # --- dV/dCoeff（=dBlos/dCoeff 链式法则）---
        if calcDV == 1:
            VsigKernelCont = self.VsigKernel * contin[
                np.newaxis, :]  # (N_vel, N_cells)
            # 与 profile.py 完全相同的 einsum 模式
            self.dVsum = np.einsum('ijl,kl->ikj', dBlos_d, VsigKernelCont)
            self.dVsum = np.swapaxes(self.dVsum, 1,
                                     2)  # (n_types, nTot, N_vel)
            self.dVsum *= invSumIc0
        else:
            self.dVsum = 0

        # --- dI/dBright ---
        if calcDI == 1:
            self.dIcsum = (surfaceScale * invSumIc0)[:, np.newaxis] * (
                self.Iunscaled.T - self.I[np.newaxis, :])  # (N_cells, N_vel)
        else:
            self.dIcsum = 0

        # --- dV/dBright ---
        if calcDI == 1 and calcDV == 1:
            V_unscaled = Blos[
                np.newaxis, :] * self.VsigKernel  # (N_vel, N_cells)
            self.dVdBright = (surfaceScale * invSumIc0)[:, np.newaxis] * (
                V_unscaled.T - self.V[np.newaxis, :])  # (N_cells, N_vel)
        else:
            self.dVdBright = 0

    def dopplerShift(self, vel_shift: float) -> None:
        """对波长网格施加 Doppler 偏移（+ 为远离观测者方向）。"""
        self.wl = self.wl + self.wl * vel_shift / _C_KMS


# ---------------------------------------------------------------------------
# diskIntProfAndDerivHalpha
# ---------------------------------------------------------------------------


class diskIntProfAndDerivHalpha:
    """H-alpha 盘积分轮廓及全部一阶导数。

    接口与 ``diskIntProfAndDeriv`` / ``diskIntProfAndDerivUnno`` 完全对等：

    - ``IIc``    : 归一化 Stokes I，形状 (N_vel,)
    - ``VIc``    : 归一化 Stokes V，形状 (N_vel,)
    - ``dVIc``   : dV/d(mag coeff)，形状 (n_types, nTot, N_vel)
    - ``dIIc``   : dI/d(brightness)，形状 (N_cells, N_vel)
    - ``dVdBri`` : dV/d(brightness)，形状 (N_cells, N_vel)
    - ``dImag``  : dI/d(mag coeff)（弱场近似下设为 0）

    与 Voigt 版本的核心区别：仪器展宽**始终使用显式卷积**（不折入 Gaussian 宽度）。
    """

    def __init__(
        self,
        visible_grid,
        v_mag_cart: np.ndarray,
        d_mag_cart,
        bright_map,
        ldata: lineDataHalpha,
        vel_eq: float,
        wl_grid: np.ndarray,
        calcDI: int,
        calcDV: int,
    ) -> None:
        self.wl = wl_grid
        self.numPts = wl_grid.shape[0]
        self.wlStart = wl_grid[0]
        self.wlEnd = wl_grid[-1]

        self.IIc = np.zeros(self.numPts)
        self.QIc = np.zeros(self.numPts)
        self.UIc = np.zeros(self.numPts)
        self.VIc = np.zeros(self.numPts)

        # 各格点的旋转 Doppler 移位
        vr_cell = vel_eq * visible_grid.velRotProj  # (N_cells,)
        self.wlCells = np.outer(self.wl, 1.0 / (1.0 + vr_cell / _C_KMS))

        self.prof = localProfileAndDerivHalpha(ldata, self.numPts,
                                               self.wlCells)
        self.updateIntProfDeriv(visible_grid, v_mag_cart, d_mag_cart,
                                bright_map, ldata, calcDI, calcDV)

    def updateIntProfDeriv(
        self,
        visible_grid,
        v_mag_cart: np.ndarray,
        d_mag_cart,
        bright_map,
        ldata: lineDataHalpha,
        calcDI: int,
        calcDV: int,
    ) -> None:
        """按新几何/亮度/磁场重新计算盘积分轮廓及导数。"""
        self.calcDI = calcDI
        self.calcDV = calcDV

        # --- 视线磁场分量 Blos 及其导数 ---
        Blos, dBlos_d = self._BlosProjected(visible_grid.vViewCart, v_mag_cart,
                                            d_mag_cart)

        # --- 表面缩放 ---
        limb_d = limbDarkening(ldata.limbDark, visible_grid.viewAngle)
        g_dark = visible_grid.gravityDarkening(ldata.gravDark)
        surface_scale = limb_d * g_dark * visible_grid.projArea * visible_grid.visible

        self.prof.updateProfDeriv(
            ldata,
            Blos,
            dBlos_d,
            bright_map,
            self.numPts,
            self.wlCells,
            surface_scale,
            calcDI,
            calcDV,
        )

        self.IIc = self.prof.I
        self.VIc = self.prof.V
        self.dVIc = self.prof.dVsum
        self.dIIc = self.prof.dIcsum
        self.dImag = 0  # 弱场近似下不使用
        self.dVdBri = self.prof.dVdBright

    def _BlosProjected(
            self,
            v_view: np.ndarray,  # (3, N_cells)
            v_mag_cart: np.ndarray,  # (3, N_cells)
            d_mag_cart,  # (3, n_types, nTot, N_cells) 或 0
    ) -> tuple:
        """计算视线磁场分量 Blos 及其对球谐系数的导数。

        Returns
        -------
        Blos    : (N_cells,)
        dBlos_d : (n_types, nTot, N_cells) 或 0
        """
        # Blos = dot(v_view, v_mag_cart)，逐格点
        Blos = np.einsum('ij,ij->j', v_view, v_mag_cart)  # (N_cells,)

        if self.calcDV == 1:
            dBlos_d = np.einsum('il,ijkl->jkl',
                                v_view,
                                d_mag_cart,
                                optimize=True)
        else:
            dBlos_d = 0

        return Blos, dBlos_d

    def convolveIGnumpy(self, fwhm: float) -> None:
        """用高斯仪器轮廓对盘积分谱显式卷积（分辨率 R = fwhm）。

        H-alpha 模型始终使用显式卷积（仪器展宽不折入 Gaussian 宽度中），
        fwhm ≤ 0 时不卷积。
        """
        if fwhm <= 0.0:
            return

        wl_center = (self.wlStart + self.wlEnd) / 2.0
        wl_step = (self.wlEnd - self.wlStart) / (self.numPts - 1)
        fwhm_wl = wl_center / fwhm
        wl_end_g = 3.0 * fwhm_wl
        num_pts_g = 2 * int(wl_end_g / wl_step) + 1

        if fwhm_wl < wl_step:
            return  # 分辨率极高，卷积核小于一个像素，跳过

        wl_g = np.linspace(-wl_end_g, wl_end_g, num_pts_g)
        prof_g = (0.939437 / fwhm_wl) * np.exp(-2.772589 * (wl_g / fwhm_wl)**2)
        prof_g /= np.sum(prof_g)
        pad = num_pts_g // 2

        def _pad_conv(arr: np.ndarray) -> np.ndarray:
            padded = np.concatenate(
                [np.repeat(arr[:1], pad), arr,
                 np.repeat(arr[-1:], pad)])
            return np.convolve(padded, prof_g, mode='valid')

        self.IIc = _pad_conv(self.IIc)
        self.VIc = _pad_conv(self.VIc)

        if self.calcDI == 1:
            n_pad = pad
            tmpdI = np.concatenate([
                np.tile(self.dIIc[:, :1], (1, n_pad)),
                self.dIIc,
                np.tile(self.dIIc[:, -1:], (1, n_pad)),
            ],
                                   axis=1)
            flat = np.convolve(tmpdI.ravel(), prof_g, mode='same')
            self.dIIc = flat.reshape(tmpdI.shape)[:, n_pad:-n_pad]

        if self.calcDV == 1:
            n_pad = pad
            tmpdV = np.concatenate([
                np.tile(self.dVIc[:, :, :1], (1, 1, n_pad)),
                self.dVIc,
                np.tile(self.dVIc[:, :, -1:], (1, 1, n_pad)),
            ],
                                   axis=2)
            flat = np.convolve(tmpdV.ravel(), prof_g, mode='same')
            self.dVIc = flat.reshape(tmpdV.shape)[:, :, n_pad:-n_pad]

        if self.calcDI == 1 and self.calcDV == 1:
            n_pad = pad
            tmpdV = np.concatenate([
                np.tile(self.dVdBri[:, :1], (1, n_pad)),
                self.dVdBri,
                np.tile(self.dVdBri[:, -1:], (1, n_pad)),
            ],
                                   axis=1)
            flat = np.convolve(tmpdV.ravel(), prof_g, mode='same')
            self.dVdBri = flat.reshape(tmpdV.shape)[:, n_pad:-n_pad]


# ---------------------------------------------------------------------------
# 批量初始化辅助
# ---------------------------------------------------------------------------


def getAllProfDirivHalpha(
    par,
    list_grid_view,
    vec_mag_cart: np.ndarray,
    d_mag_cart0,
    bri_map,
    ldata: lineDataHalpha,
    wl_syn_set: list,
) -> list:
    """为各观测相位构造 ``diskIntProfAndDerivHalpha`` 列表。

    接口与 ``getAllProfDiriv`` / ``getAllProfDirivUnno`` 完全对等，
    可在 pipeline 中直接替换。

    Parameters
    ----------
    par : ZDIConfig
        参数对象（需有 ``cycleList``, ``velEq``, ``instrumentRes``）。
    list_grid_view : list
        各相位可见性几何视图列表。
    vec_mag_cart : (3, N_cells) ndarray
        初始磁场笛卡尔分量。
    d_mag_cart0 : ndarray 或 0
        磁场对球谐系数的导数；``fitMag=0`` 时传 0。
    bri_map : object
        亮度图对象。
    ldata : lineDataHalpha
        H-alpha 线参数对象。
    wl_syn_set : list[ndarray]
        各相位合成波长网格列表。

    Returns
    -------
    list[diskIntProfAndDerivHalpha]
    """
    set_syn_spec = []
    for nobs, _phase in enumerate(par.cycleList):
        spec = diskIntProfAndDerivHalpha(
            list_grid_view[nobs],
            vec_mag_cart,
            d_mag_cart0,
            bri_map,
            ldata,
            par.velEq,
            wl_syn_set[nobs],
            0,  # calcDI=0 in initial build（由 pipeline 按需开启）
            par.calcDV,
        )
        spec.convolveIGnumpy(par.instrumentRes)
        set_syn_spec.append(spec)
    return set_syn_spec
