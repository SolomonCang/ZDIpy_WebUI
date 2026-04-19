"""core.line_models.unno — Unno-Rachkovsky (1967) 谱线轮廓模型。

基于 Milne-Eddington 大气的完整偏振辐射转移（Stokes I/V），
适用于强磁场或需要全 Unno-Rachkovsky 解的场景。

参考：Landi Degl'Innocenti & Landolfi (2004)
      《Polarization in Spectral Lines》Eqs 9.32 / 9.110–9.112
参考实现：J. Morin / J.F. Donati (lineprofileUnnoR.py)

对外接口与 ``core.line_models.profile`` 对等，可作为可选模型替换进 pipeline：

- ``lineDataUnno``            : Unno-Rachkovsky 谱线参数容器
  （含 β 源函数坡度、fI / fV 填充因子）
- ``localProfileAndDerivUnno``: 单格点 Unno 轮廓及其导数（预计算安静区 sI0）
- ``diskIntProfAndDerivUnno`` : 盘积分轮廓及全部一阶导数
- ``getAllProfDirivUnno``      : 批量初始化辅助函数

与 Voigt 弱场模型的主要区别
----------------------------
- Stokes I 同样依赖于磁场（通过 etaI/Delta），因此存在 dI/dB 导数
  （``dImag`` 属性），但现有 MEM 优化器不使用该项（已设为 0 以保持兼容）。
- 仪器展宽默认使用**显式卷积**（积分后处理），而非折入 Gaussian 宽度。
- 需要完整笛卡尔磁场向量以计算 Bmod（场强）和 Btheta（场方向角）。
"""
from __future__ import annotations

import numpy as np

from core.line_models.line_utils import limbDarkening

_C_KMS: float = 2.99792458e5  # speed of light (km/s)
_ZEEMAN_CONST: float = 4.66864e-12  # e/(4π m_e c²), B in G, λ in nm

# ---------------------------------------------------------------------------
# Humlicek (1982) complex error function — Voigt + Faraday-Voigt
# ---------------------------------------------------------------------------


def _voigt_faraday_humlicek(
    wl_gauss_norm: np.ndarray,
    width_lorentz: float,
) -> np.ndarray:
    """Humlicek (1982) 近似 Faddeeva 函数，精度约 10^{-4}。

    Parameters
    ----------
    wl_gauss_norm:
        距线心距离（以 Gaussian 宽度为单位），任意形状。
    width_lorentz:
        Lorentzian 宽度 / Gaussian 宽度的比值（标量）。

    Returns
    -------
    w4 : complex ndarray, 与 ``wl_gauss_norm`` 同形状。
        实部 = Voigt 轮廓 H；虚部 = 2 × Faraday-Voigt 函数。
    """
    z = width_lorentz - 1j * wl_gauss_norm
    zz = z * z
    s = np.abs(wl_gauss_norm) + width_lorentz
    w4 = np.zeros(wl_gauss_norm.shape, dtype=complex)

    con1 = s >= 15.0
    zt = z[con1]
    w4[con1] = 0.56418958355 * zt / (0.5 + zt * zt)

    con2 = (s >= 5.5) & (s < 15.0)
    zt = z[con2]
    zzt = zz[con2]
    w4[con2] = zt * (1.4104739589 + 0.56418958355 * zzt) / (
        (3.0 + zzt) * zzt + 0.7499999999)

    con3 = (width_lorentz >= 0.195 * np.abs(wl_gauss_norm) - 0.176) & (s < 5.5)
    zt = z[con3]
    w4[con3] = (16.4954955 + zt *
                (20.2093334 + zt *
                 (11.9648172 + zt * (3.77898687 + zt * 0.564223565)))) / (
                     16.4954955 + zt * (38.8236274 + zt *
                                        (39.2712051 + zt *
                                         (21.6927370 + zt *
                                          (6.69939801 + zt)))))

    con4 = w4 == 0.0 + 0j
    zt = z[con4]
    zzt = zz[con4]
    w4[con4] = np.exp(zzt) - zt * (36183.30536 - zzt *
                                   (3321.990492 - zzt *
                                    (1540.786893 - zzt *
                                     (219.0312964 - zzt *
                                      (35.76682780 - zzt *
                                       (1.320521697 - zzt * 0.5641900381)))))
                                   ) / (32066.59372 - zzt *
                                        (24322.84021 - zzt *
                                         (9022.227659 - zzt *
                                          (2186.181081 - zzt *
                                           (364.2190727 - zzt *
                                            (61.57036588 - zzt *
                                             (1.841438936 - zzt)))))))

    return w4


# ---------------------------------------------------------------------------
# Unno-Rachkovsky local profile
# ---------------------------------------------------------------------------


def _unno_profile(
    width_gauss: float,
    width_lorentz: float,
    ldata: "lineDataUnno",
    wls: np.ndarray,
    view_angle: np.ndarray,
    Bmod: np.ndarray,
    Btheta: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """单格点 Unno-Rachkovsky Stokes I 和 V 轮廓。

    Parameters
    ----------
    width_gauss:
        Gaussian 多普勒宽度（km/s）。
    width_lorentz:
        Lorentzian 宽度 / Gaussian 宽度（无量纲比值）。
    ldata:
        ``lineDataUnno`` 对象（取索引 [0] 的参数）。
    wls: (N_vel, N_cells)
        各格点的波长网格（nm），已按 Doppler 移位修正。
    view_angle: (N_cells,)
        视线与格点法线的夹角（rad）。
    Bmod: (N_cells,)
        磁场强度（G）。
    Btheta: (N_cells,)
        磁场方向与视线方向的夹角（rad）。

    Returns
    -------
    sI : (N_vel, N_cells)
        连续谱归一化 Stokes I。
    sV : (N_vel, N_cells)
        连续谱归一化 Stokes V。
    """
    kL = float(ldata.str[0])
    beta = float(ldata.beta[0])
    fV_val = float(ldata.fV[0])

    # mu = cos(view_angle)，用于 Milne-Eddington 源函数
    mu = np.cos(view_angle)  # (N_cells,)
    theta = Btheta  # (N_cells,) — 磁场与视线夹角
    chi = 0.0  # 磁场方位角方位角固定为 0（U=0）

    wl0 = float(ldata.wl0[0])
    g_eff = float(ldata.g[0])

    # Zeeman 分裂量（nm），sigma 组分偏移
    wl_b_sig = _ZEEMAN_CONST * g_eff * wl0**2 * Bmod / fV_val  # (N_cells,)

    # Gaussian 宽度换算为波长单位
    width_gauss_wl = width_gauss / _C_KMS * wl0

    # 各组分的线心偏移（Gaussian 宽度为单位）
    wl_gauss_norm = (wls - wl0) / width_gauss_wl  # (N_vel, N_cells)
    wl_gauss_norm_L = (wls - wl0 + wl_b_sig) / width_gauss_wl
    wl_gauss_norm_R = (wls - wl0 - wl_b_sig) / width_gauss_wl

    # --- Voigt + Faraday-Voigt 函数，各 Zeeman 组分 ---
    norm = 1.0 / (np.sqrt(np.pi) * width_gauss)

    W = _voigt_faraday_humlicek(wl_gauss_norm, width_lorentz)
    eta_P = W.real * norm
    rho_P = W.imag * norm

    W = _voigt_faraday_humlicek(wl_gauss_norm_L, width_lorentz)
    eta_L = W.real * norm
    rho_L = W.imag * norm

    W = _voigt_faraday_humlicek(wl_gauss_norm_R, width_lorentz)
    eta_R = W.real * norm
    rho_R = W.imag * norm

    # --- UR 吸收—色散系数（Landi Degl'Innocenti & Landolfi 2004 式 9.32） ---
    # theta 为 (N_cells,)，eta_P 等为 (N_vel, N_cells)，广播正确
    sin2 = np.sin(theta)**2  # (N_cells,)
    cos2 = np.cos(theta)**2  # (N_cells,)
    eta_half = kL * 0.5

    # "1 + etaI" in Landi notation
    etaI = (1.0 + eta_half * (eta_P * sin2 + 0.5 * (eta_L + eta_R) *
                              (1.0 + cos2)))
    etaQ = eta_half * (eta_P - 0.5 *
                       (eta_L + eta_R)) * sin2 * np.cos(2.0 * chi)
    etaV = eta_half * (eta_R - eta_L) * np.cos(theta)

    rhoQ = eta_half * (rho_P - 0.5 *
                       (rho_L + rho_R)) * sin2 * np.cos(2.0 * chi)
    rhoV = eta_half * (rho_R - rho_L) * np.cos(theta)

    # --- 行列式 Delta（chi=0，故 etaU=rhoU=0）---
    Delta = (etaI**2) * (etaI**2 + rhoQ**2 + rhoV**2 - etaQ**2 - etaV**2) \
            - (etaQ * rhoQ + etaV * rhoV)**2

    # --- Stokes I 和 V（连续谱归一，式 9.110 & 9.112）---
    beta_mu = beta * mu  # (N_cells,)
    one_plus_bmu = 1.0 + beta_mu  # (N_cells,)

    # 防止除零（mu≈0 的临边区域）
    safe_denom = np.where(np.abs(one_plus_bmu) > 1e-14, one_plus_bmu, 1e-14)

    sI = (1.0 + beta_mu / Delta * etaI *
          (etaI**2 + rhoQ**2 + rhoV**2)) / safe_denom
    sV = -(beta_mu / safe_denom) / Delta * (etaV * etaI**2 + rhoV *
                                            (etaQ * rhoQ + etaV * rhoV))

    return sI, sV


# ---------------------------------------------------------------------------
# lineDataUnno
# ---------------------------------------------------------------------------


class lineDataUnno:
    """Unno-Rachkovsky 谱线参数容器。

    在 ``core.line_models.profile.lineData`` 的基础上增加三个 UR 特定字段：

    - ``beta`` : Planck 函数随光深的线性坡度；若 ≤0 则从 ``limbDark`` 推导。
    - ``fI``   : Stokes I 填充因子（磁活动区面积分数，默认 1.0）。
    - ``fV``   : Stokes V 填充因子 / Zeeman 分裂修正因子（默认 1.0）。

    通常通过 ``from_parameters()`` 类方法由 config.json 参数构造，
    也可通过 ``from_file()`` 从扩展格式文件读取。
    """
    def __init__(self) -> None:
        self.wl0 = np.array([])
        self.str = np.array([])
        self.beta = np.array([])
        self.widthGauss = np.array([])
        self.widthLorentz = np.array([])
        self.g = np.array([])
        self.limbDark = np.array([])
        self.fI = np.array([])
        self.fV = np.array([])
        self.gravDark = np.array([])
        self.numLines = 0
        self.instRes: float = -1.0
        self.macro_turb_kms: float = 0.0  # 宏观湍流 Gaussian FWHM（km/s），仅展宽 Stokes I

    @classmethod
    def from_parameters(
        cls,
        wavelength_nm: float,
        line_strength: float,
        gauss_width_kms: float,
        lorentz_width_fraction: float,
        lande_g: float,
        limb_darkening: float,
        gravity_darkening: float,
        instRes: float = -1.0,
        beta: float = -1.0,
        filling_factor_I: float = 1.0,
        filling_factor_V: float = 1.0,
        macro_turb_kms: float = 0.0,
    ) -> "lineDataUnno":
        """由 config.json 数值字段构造 ``lineDataUnno``。

        Parameters
        ----------
        beta:
            Planck 函数坡度。≤0 表示自动由 ``limb_darkening`` 推导：
            ``beta = limbd / (1 - limbd)``。
        filling_factor_I:
            Stokes I 填充因子（磁活动区分数），默认 1.0。
        filling_factor_V:
            Stokes V 填充因子，默认 1.0。其他参数定义同 ``lineData.from_parameters``。
        """
        obj = cls()
        obj.instRes = float(instRes)
        obj.numLines = 1
        obj.wl0 = np.array([float(wavelength_nm)])
        obj.str = np.array([float(line_strength)])
        obj.widthGauss = np.array([float(gauss_width_kms)])
        obj.widthLorentz = np.array([float(lorentz_width_fraction)])
        obj.g = np.array([float(lande_g)])
        obj.limbDark = np.array([float(limb_darkening)])
        obj.gravDark = np.array([float(gravity_darkening)])

        # 填充因子
        fI_val = float(filling_factor_I)
        fV_val = float(filling_factor_V)
        obj.fI = np.array([max(0.0, min(1.0, fI_val))])
        obj.fV = np.array([max(1e-6, fV_val)])  # 防止除零

        # Beta（Planck 函数坡度）
        beta_val = float(beta)
        if beta_val <= 0.0:
            limbd = float(limb_darkening)
            if abs(limbd - 1.0) < 1e-10:
                limbd = 1.0 - 1e-6  # 防止除零
            beta_val = limbd / (1.0 - limbd)
            print(
                f'lineDataUnno: using beta derived from limb darkening: {beta_val:.4f}'
            )
        obj.beta = np.array([beta_val])

        # 宏观湍流（仅展宽 Stokes I，不影响 V）
        obj.macro_turb_kms = max(0.0, float(macro_turb_kms))

        return obj

    @classmethod
    def from_file(cls, filename: str, instRes: float = -1.0) -> "lineDataUnno":
        """从文件读取 Unno-Rachkovsky 谱线参数。

        文件每行格式（'#' 开头为注释）::

            wl0  lineStr  beta  widthGauss  widthLorentz  g  limbDark  fI  fV  [gravDark]
        """
        obj = cls()
        obj.instRes = float(instRes)

        with open(filename, 'r') as fh:
            for line in fh:
                stripped = line.strip()
                if not stripped or stripped[0] == '#':
                    continue
                parts = stripped.split()
                obj.numLines += 1
                obj.wl0 = np.append(obj.wl0, float(parts[0]))
                obj.str = np.append(obj.str, float(parts[1]))

                beta_raw = float(parts[2])
                if beta_raw <= 0.0:
                    limbd = float(parts[6])
                    if abs(limbd - 1.0) < 1e-10:
                        limbd = 1.0 - 1e-6
                    beta_raw = limbd / (1.0 - limbd)
                    print(
                        f'lineDataUnno: auto-computing beta from limb darkening: {beta_raw:.4f}'
                    )
                obj.beta = np.append(obj.beta, beta_raw)

                obj.widthGauss = np.append(obj.widthGauss, float(parts[3]))
                obj.widthLorentz = np.append(obj.widthLorentz,
                                             float(parts[4]) / float(parts[3]))
                obj.g = np.append(obj.g, float(parts[5]))
                obj.limbDark = np.append(obj.limbDark, float(parts[6]))

                fI_raw = float(parts[7]) if len(parts) > 7 else 1.0
                fV_raw = float(parts[8]) if len(parts) > 8 else 1.0
                obj.fI = np.append(obj.fI, max(0.0, min(1.0, fI_raw)))
                obj.fV = np.append(obj.fV, max(1e-6, fV_raw))

                try:
                    obj.gravDark = np.append(obj.gravDark, float(parts[9]))
                except (ValueError, IndexError):
                    obj.gravDark = np.append(obj.gravDark, 0.0)

        if obj.numLines > 1:
            print(
                'ERROR: multi-line input not currently supported for Unno model.'
            )
        return obj


# ---------------------------------------------------------------------------
# localProfileAndDerivUnno
# ---------------------------------------------------------------------------


class localProfileAndDerivUnno:
    """预计算安静区（B=0）Stokes I，供 ``updateProfDeriv`` 重复使用。

    与 ``localProfileAndDeriv`` 的关键区别：

    1. 构造时需要 ``view_angle``（视线与格点法线夹角），用于 Unno 源函数。
    2. 磁场参数为 ``Bmod``（场强）和 ``Btheta``（场方向角），而非弱场近似的 Blos。
    3. 计算安静区 ``sI0``（fI<1 的混合归一化）。

    Parameters
    ----------
    ldata:
        ``lineDataUnno`` 对象。
    numPts:
        波长点数（保留兼容性，实际从 wl_grid 推断）。
    wl_grid: (N_vel, N_cells)
        各格点已 Doppler 移位后的波长网格（nm）。
    view_angle: (N_cells,)
        视线与格点法线夹角（rad）。
    """
    def __init__(
        self,
        ldata: lineDataUnno,
        numPts: int,
        wl_grid: np.ndarray,
        view_angle: np.ndarray,
    ) -> None:
        self.wl = wl_grid  # (N_vel, N_cells)
        self.view_angle = view_angle  # (N_cells,)

        # 对 Unno 模型始终使用原始 Gaussian 宽度（不折入仪器展宽）
        self.width_gauss = float(ldata.widthGauss[0])
        self.width_lorentz = float(ldata.widthLorentz[0])

        # 预计算安静区（B=0）轮廓，用于磁/非磁混合（填充因子 fI < 1）
        b_zero = np.zeros_like(view_angle)  # (N_cells,)
        self.sI0, _ = _unno_profile(
            self.width_gauss,
            self.width_lorentz,
            ldata,
            self.wl,
            self.view_angle,
            b_zero,
            b_zero,
        )  # (N_vel, N_cells)

    def updateProfDeriv(
        self,
        ldata: lineDataUnno,
        Bmod: np.ndarray,
        Btheta: np.ndarray,
        dBmod_d,
        dBtheta_d,
        brightMap,
        numPts: int,
        wl_grid: np.ndarray,
        surface_scale: np.ndarray,
        calcDI: int,
        calcDV: int,
    ) -> None:
        """计算盘积分轮廓及导数（求和与连续谱归一化均在此完成）。

        Parameters
        ----------
        Bmod: (N_cells,)
            各格点磁场强度（G）。
        Btheta: (N_cells,)
            各格点磁场与视线夹角。
        dBmod_d: (n_types, nTot, N_cells) 或 0
            Bmod 对球谐系数的导数。
        dBtheta_d: (n_types, nTot, N_cells) 或 0
            Btheta 对球谐系数的导数。
        brightMap:
            亮度图对象（``bright`` 属性为 (N_cells,)）。
        surface_scale: (N_cells,)
            表面缩放（limb darkening × gravity darkening × projArea × visible）。
        calcDI, calcDV:
            是否计算亮度/磁场导数的开关。
        """
        sI, sV = _unno_profile(
            self.width_gauss,
            self.width_lorentz,
            ldata,
            self.wl,
            self.view_angle,
            Bmod,
            Btheta,
        )  # (N_vel, N_cells) each

        # --- 局部连续谱缩放：limb darkening + 亮度图 ---
        contin = surface_scale * brightMap.bright  # (N_cells,)

        # --- 含填充因子的混合 Stokes I 和 V ---
        fI_val = float(ldata.fI[0])
        fV_val = float(ldata.fV[0])
        I_unscaled = sI * fI_val + self.sI0 * (1.0 - fI_val
                                               )  # (N_vel, N_cells)
        V_unscaled = sV * fV_val  # (N_vel, N_cells)

        # --- 连续谱总量（相同对所有波长点）---
        self.Ic = np.sum(contin)
        inv_ic = 1.0 / self.Ic

        # --- 盘积分并归一化 ---
        self.I = np.sum(contin[np.newaxis, :] * I_unscaled, axis=1) * inv_ic
        self.V = np.sum(contin[np.newaxis, :] * V_unscaled, axis=1) * inv_ic

        # --- 磁场导数（dV/dCoeff）---
        if calcDV == 1:
            delta_bmod = 1e-8
            delta_btheta = 1e-4
            bmean = float(np.average(Bmod))
            if bmean > 10.0:
                delta_bmod = bmean * 1e-8

            # 数值导数：对 Bmod
            sIp, sVp = _unno_profile(
                self.width_gauss,
                self.width_lorentz,
                ldata,
                self.wl,
                self.view_angle,
                Bmod + delta_bmod,
                Btheta,
            )
            dV_dBmod = contin[np.newaxis, :] * fV_val * (sVp - sV) / delta_bmod
            if calcDI == 1:
                dI_dBmod = contin[np.newaxis, :] * fI_val * (sIp -
                                                             sI) / delta_bmod

            # 数值导数：对 Btheta
            sIp, sVp = _unno_profile(
                self.width_gauss,
                self.width_lorentz,
                ldata,
                self.wl,
                self.view_angle,
                Bmod,
                Btheta + delta_btheta,
            )
            dV_dBtheta = contin[np.newaxis, :] * fV_val * (sVp -
                                                           sV) / delta_btheta
            if calcDI == 1:
                dI_dBtheta = contin[np.newaxis, :] * fI_val * (
                    sIp - sI) / delta_btheta

            # 链式法则：dV/dCoeff = dV/dBmod * dBmod/dCoeff + dV/dBtheta * dBtheta/dCoeff
            # dBmod_d, dBtheta_d: (n_types, nTot, N_cells)  → indices i,j,l
            # dV_dBmod:           (N_vel, N_cells)           → indices k,l
            # einsum 'ijl,kl->ijk': output (n_types, nTot, N_vel)
            self.dVsum = (
                np.einsum('ijl,kl->ijk', dBmod_d, dV_dBmod, optimize=True) +
                np.einsum('ijl,kl->ijk', dBtheta_d, dV_dBtheta,
                          optimize=True)) * inv_ic

            if calcDI == 1:
                self.dIdBsum = (
                    np.einsum('ijl,kl->ijk', dBmod_d, dI_dBmod, optimize=True)
                    + np.einsum(
                        'ijl,kl->ijk', dBtheta_d, dI_dBtheta,
                        optimize=True)) * inv_ic
            else:
                self.dIdBsum = 0
        else:
            self.dVsum = 0
            self.dIdBsum = 0

        # --- 亮度导数（dI/dBright）---
        if calcDI == 1:
            self.dIcsum = (surface_scale * inv_ic)[:, np.newaxis] * (
                I_unscaled.T - self.I[np.newaxis, :])
        else:
            self.dIcsum = 0

        # --- 交叉导数（dV/dBright，同时 calcDI 和 calcDV 为真时才有效）---
        if calcDI == 1 and calcDV == 1:
            self.dVdBright = (surface_scale * inv_ic)[:, np.newaxis] * (
                V_unscaled.T - self.V[np.newaxis, :])
        else:
            self.dVdBright = 0

    def dopplerShift(self, vel_shift: float) -> None:
        """对波长网格施加 Doppler 偏移（+ 为远离观测者方向）。"""
        self.wl = self.wl + self.wl * vel_shift / _C_KMS


# ---------------------------------------------------------------------------
# diskIntProfAndDerivUnno
# ---------------------------------------------------------------------------


class diskIntProfAndDerivUnno:
    """Unno-Rachkovsky 盘积分轮廓及全部一阶导数。

    接口与 ``diskIntProfAndDeriv`` 对等，可直接用于 pipeline 和 MEM 反演：

    - ``IIc``    : 归一化 Stokes I
    - ``VIc``    : 归一化 Stokes V
    - ``dVIc``   : dV/d(mag coeff)，形状 (n_types, nTot, N_vel)
    - ``dIIc``   : dI/d(brightness)，形状 (N_cells, N_vel)（供 ``packResponseMatrix`` 用）
    - ``dVdBri`` : dV/d(brightness)，形状 (N_cells, N_vel)
    - ``dImag``  : dI/d(mag coeff)（现有 MEM 优化器不使用，设为 0）

    与 Voigt 版本的核心区别：``updateIntProfDeriv`` 使用 ``BdBprojected()``
    计算 Bmod 和 Btheta，而非仅取 B 的视线分量。
    """
    def __init__(
        self,
        visible_grid,
        v_mag_cart: np.ndarray,
        d_mag_cart,
        bright_map,
        ldata: lineDataUnno,
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
        # 每格点的波长网格（逆 Doppler 移位）
        self.wlCells = np.outer(self.wl, 1.0 / (1.0 + vr_cell / _C_KMS))

        self.prof = localProfileAndDerivUnno(ldata, self.numPts, self.wlCells,
                                             visible_grid.viewAngle)
        self.updateIntProfDeriv(visible_grid, v_mag_cart, d_mag_cart,
                                bright_map, ldata, calcDI, calcDV)

    def updateIntProfDeriv(
        self,
        visible_grid,
        v_mag_cart: np.ndarray,
        d_mag_cart,
        bright_map,
        ldata: lineDataUnno,
        calcDI: int,
        calcDV: int,
    ) -> None:
        """按新几何/亮度/磁场重新计算盘积分轮廓及导数。"""
        self.calcDI = calcDI
        self.calcDV = calcDV

        # --- B 向量向视线方向的投影及其导数 ---
        Bmod, Btheta, dBmod_d, dBtheta_d = self._BdBprojected(
            visible_grid.vViewCart, v_mag_cart, d_mag_cart)

        # --- 表面缩放（limb darkening × gravity darkening × projArea × visible）---
        limb_d = limbDarkening(ldata.limbDark, visible_grid.viewAngle)
        g_dark = visible_grid.gravityDarkening(ldata.gravDark)
        surface_scale = limb_d * g_dark * visible_grid.projArea * visible_grid.visible

        self.prof.updateProfDeriv(
            ldata,
            Bmod,
            Btheta,
            dBmod_d,
            dBtheta_d,
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
        self.dImag = self.prof.dIdBsum  # 存储但 MEM 优化器不使用
        self.dVdBri = self.prof.dVdBright

    def _BdBprojected(
            self,
            v_view: np.ndarray,  # (3, N_cells)
            v_mag_cart: np.ndarray,  # (3, N_cells)
            d_mag_cart,  # (3, n_types, nTot, N_cells) 或 0
    ) -> tuple:
        """计算 Bmod、Btheta 及其对球谐系数的导数。

        Returns
        -------
        Bmod    : (N_cells,)
        Btheta  : (N_cells,)
        dBmod_d : (n_types, nTot, N_cells) 或 0
        dBtheta_d: (n_types, nTot, N_cells) 或 0
        """
        # B 沿视线方向的投影（点积，逐格点）
        proj_blos = np.einsum('ij,ij->j', v_view, v_mag_cart)  # (N_cells,)

        # 磁场强度
        Bmod = np.linalg.norm(v_mag_cart, axis=0)  # (N_cells,)
        # 防止除零（B=0 的格点）
        Bmod_safe = np.where(Bmod > 0.0, Bmod, 1e-10)

        cos_theta = proj_blos / Bmod_safe  # (N_cells,)
        # 夹紧至 [-1, 1] 防止数值误差导致 arccos 失败
        cos_theta = np.clip(cos_theta, -1.0, 1.0)
        Btheta = np.arccos(cos_theta)  # (N_cells,)

        if self.calcDV == 1:
            # d(B_los)/d(coeff)，形状 (n_types, nTot, N_cells)
            dBlos_d = np.einsum('il,ijkl->jkl',
                                v_view,
                                d_mag_cart,
                                optimize=True)
            # d(|B|)/d(coeff)
            dBmod_d = np.einsum(
                'il,ijkl->jkl', v_mag_cart, d_mag_cart,
                optimize=True) / Bmod_safe
            # d(theta)/d(coeff)：通过链式法则
            sin_theta_sq = 1.0 - cos_theta**2
            sin_theta_safe = np.where(sin_theta_sq > 1e-20,
                                      np.sqrt(np.maximum(sin_theta_sq, 0.0)),
                                      1e-10)
            safe_denom = Bmod_safe**2 * sin_theta_safe  # (N_cells,)
            dBtheta_d = (proj_blos * dBmod_d - Bmod * dBlos_d) / safe_denom
        else:
            dBmod_d = 0
            dBtheta_d = 0

        return Bmod, Btheta, dBmod_d, dBtheta_d

    def convolveIGnumpy(self, fwhm: float) -> None:
        """用高斯仪器轮廓对盘积分谱进行显式卷积（分辨率 R = fwhm）。

        Unno 模型始终使用显式卷积（不折入 Gaussian 宽度），
        因此该方法在 fwhm > 0 时总是执行，除非 FWHM 小于一个像素。
        """
        if fwhm <= 0.0:
            return

        wl_center = (self.wlStart + self.wlEnd) / 2.0
        wl_step = (self.wlEnd - self.wlStart) / (self.numPts - 1)
        fwhm_wl = wl_center / fwhm
        wl_end_g = 3.0 * fwhm_wl
        num_pts_g = 2 * int(wl_end_g / wl_step) + 1

        if fwhm_wl < wl_step:
            return  # 分辨率极高，跳过卷积

        wl_g = np.linspace(-wl_end_g, wl_end_g, num_pts_g)
        prof_g = (0.939437 / fwhm_wl) * np.exp(-2.772589 * (wl_g / fwhm_wl)**2)
        prof_g /= np.sum(prof_g)
        pad = num_pts_g // 2

        def _pad_conv(arr: np.ndarray) -> np.ndarray:
            """对 1D 数组边缘填充后卷积。"""
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

    def convolveMacroTurbI(self, macro_turb_kms: float) -> None:
        """仅对 Stokes I 施加宏观湍流 Gaussian 展宽，Stokes V 完全不受影响。

        宏观湍流（macroturbulence）作用于盘积分后的涌现 I 轮廓（large-scale
        velocity fields average over the disk），而 Stokes V 的特征由局域磁场
        Zeeman 分裂决定，不被宏观湍流展宽。这是解决"宽 I / 锐 V"矛盾的
        物理正确途径。

        Parameters
        ----------
        macro_turb_kms:
            宏观湍流 Gaussian FWHM（km/s）。≤0 时跳过。
        """
        if macro_turb_kms <= 0.0:
            return

        wl_center = (self.wlStart + self.wlEnd) / 2.0
        wl_step = (self.wlEnd - self.wlStart) / (self.numPts - 1)
        fwhm_wl = wl_center * macro_turb_kms / _C_KMS
        wl_end_g = 3.0 * fwhm_wl
        num_pts_g = 2 * int(wl_end_g / wl_step) + 1

        if fwhm_wl < wl_step or num_pts_g < 3:
            return  # 宏观湍流宽度远小于一个像素，跳过

        wl_g = np.linspace(-wl_end_g, wl_end_g, num_pts_g)
        prof_g = (0.939437 / fwhm_wl) * np.exp(-2.772589 * (wl_g / fwhm_wl)**2)
        prof_g /= np.sum(prof_g)
        pad = num_pts_g // 2

        def _pad_conv(arr: np.ndarray) -> np.ndarray:
            padded = np.concatenate(
                [np.repeat(arr[:1], pad), arr,
                 np.repeat(arr[-1:], pad)])
            return np.convolve(padded, prof_g, mode='valid')

        # 仅卷积 IIc — VIc 保持不变
        self.IIc = _pad_conv(self.IIc)

        # 同步卷积 I 的亮度导数（dI/d(brightness)），保持 Jacobian 与展宽后 I 一致
        if self.calcDI == 1 and isinstance(self.dIIc, np.ndarray):
            n_pad = pad
            tmpdI = np.concatenate([
                np.tile(self.dIIc[:, :1], (1, n_pad)),
                self.dIIc,
                np.tile(self.dIIc[:, -1:], (1, n_pad)),
            ],
                                   axis=1)
            flat = np.convolve(tmpdI.ravel(), prof_g, mode='same')
            self.dIIc = flat.reshape(tmpdI.shape)[:, n_pad:-n_pad]

        # VIc, dVIc, dVdBri 不做卷积 — 宏观湍流不影响 Stokes V


# ---------------------------------------------------------------------------
# 批量初始化辅助
# ---------------------------------------------------------------------------


def getAllProfDirivUnno(
    par,
    list_grid_view,
    vec_mag_cart: np.ndarray,
    d_mag_cart0,
    bri_map,
    ldata: lineDataUnno,
    wl_syn_set: list,
) -> list:
    """为各观测相位构造初始化 ``diskIntProfAndDerivUnno`` 列表。

    接口与 ``getAllProfDiriv`` 对等，可在 pipeline 中直接替换。

    Parameters
    ----------
    par:
        参数对象（需有 ``cycleList``, ``velEq``, ``instrumentRes``）。
    list_grid_view:
        各相位可见性几何视图（支持索引访问）。
    vec_mag_cart: (3, N_cells)
        初始磁场笛卡尔向量。
    d_mag_cart0:
        磁场对球谐系数的导数；``fitMag=0`` 时传 0。
    bri_map:
        亮度图对象。
    ldata:
        ``lineDataUnno`` 对象。
    wl_syn_set:
        各相位合成波长网格列表。

    Returns
    -------
    list[diskIntProfAndDerivUnno]
    """
    set_syn_spec = []
    for nobs, _phase in enumerate(par.cycleList):
        spec = diskIntProfAndDerivUnno(
            list_grid_view[nobs],
            vec_mag_cart,
            d_mag_cart0,
            bri_map,
            ldata,
            par.velEq,
            wl_syn_set[nobs],
            0,
            par.calcDV,
        )
        spec.convolveIGnumpy(par.instrumentRes)
        macro_turb = getattr(ldata, 'macro_turb_kms', 0.0)
        if macro_turb > 0.0:
            spec.convolveMacroTurbI(macro_turb)
        set_syn_spec.append(spec)
    return set_syn_spec
