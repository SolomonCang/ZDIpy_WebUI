"""core.line_models.profile — ZDI 谱线轮廓主类。

提供基于 Humlicek (1982) 加速 Voigt 轮廓算法的弱场近似 Stokes I/V 合成：

- ``lineData``              : 谱线参数文件读取器
- ``localProfileAndDeriv``  : 单格点轮廓及其亮度导数（预计算 Voigt + Stokes V 核）
- ``diskIntProfAndDeriv``   : 盘积分轮廓及全部一阶导数（亮度图和磁场系数）
- ``getAllProfDiriv``        : 批量构造各相位 diskIntProfAndDeriv 列表

.. note::
    新代码应优先使用 ``core.line_models.base.LineModel`` 抽象接口。
    本模块中的类保留以保证与既有流水线（``pipeline/pipeline.py``）的向后兼容性。
"""
from __future__ import annotations

import numpy as np

from core.line_models.line_utils import limbDarkening

# When False, the instrumental profile is folded into the local Voigt width
# (more efficient). When True, an explicit post-integration convolution is done.
explicitConvolution: bool = False

_C_KMS = 2.99792458e5  # speed of light in km/s

# ---------------------------------------------------------------------------
# lineData
# ---------------------------------------------------------------------------


class lineData:
    """从文件读取谱线参数（Donati LSD 输入格式）。

    文件格式（每参数行，'#' 开头为注释）::

        wl0  lineStr  widthGauss  widthLorentz  g  limbDark  [gravDark]

    Parameters
    ----------
    inFileName : str
        谱线参数文件路径。
    instRes : float, optional
        仪器分辨率 R = λ/Δλ；< 0 表示不施加仪器展宽（默认 -1）。
    """
    def __init__(self, inFileName: str, instRes: float = -1.0) -> None:
        self.wl0 = np.array([])
        self.str = np.array([])
        self.widthGauss = np.array([])
        self.widthLorentz = np.array([])
        self.g = np.array([])
        self.limbDark = np.array([])
        self.gravDark = np.array([])
        self.numLines = 0
        self.instRes = instRes

        with open(inFileName, 'r') as inFile:
            for line in inFile:
                if not line.strip() or line.strip()[0] == '#':
                    continue
                self.numLines += 1
                parts = line.split()
                self.wl0 = np.append(self.wl0, float(parts[0]))
                self.str = np.append(self.str, float(parts[1]))
                self.widthGauss = np.append(self.widthGauss, float(parts[2]))
                self.widthLorentz = np.append(self.widthLorentz,
                                              float(parts[3]))
                self.g = np.append(self.g, float(parts[4]))
                self.limbDark = np.append(self.limbDark, float(parts[5]))
                try:
                    self.gravDark = np.append(self.gravDark, float(parts[6]))
                except (ValueError, IndexError):
                    self.gravDark = np.append(self.gravDark, 0.0)


# ---------------------------------------------------------------------------
# localProfileAndDeriv
# ---------------------------------------------------------------------------


class localProfileAndDeriv:
    """预计算单格点（逐波长格点）Voigt 核及 Stokes V 核。

    将与亮度图和磁场无关的量预计算一次，供 ``updateProfDeriv`` 高效重复使用。

    Parameters
    ----------
    lineData
        ``lineData`` 对象。
    numPts : int
        波长格点数（未使用，保留兼容性）。
    wlGrid : (N_vel, N_cells) ndarray
        各格点 Doppler 移位后的波长网格。
    """
    def __init__(self, lineData, numPts: int, wlGrid: np.ndarray) -> None:
        self.wl = wlGrid
        self.Q = np.zeros(self.wl.shape)
        self.U = np.zeros(self.wl.shape)

        # 确定有效 Gaussian 宽度（折叠仪器展宽或保持原始）
        if not explicitConvolution and lineData.instRes > 0.0:
            velInstRes = _C_KMS / lineData.instRes * 0.6005612043932249
            widthGauss = float(
                np.sqrt(lineData.widthGauss[0]**2 + velInstRes**2))
            widthLorentz = float(lineData.widthLorentz[0] *
                                 (lineData.widthGauss[0] / widthGauss))
            lineStr = float(lineData.str[0] *
                            (lineData.widthGauss[0] / widthGauss))
        elif not explicitConvolution and lineData.instRes < 0.0:
            print('Warning: convolution with an instrumental profile '
                  'may not be performed')
            widthGauss = float(lineData.widthGauss[0])
            widthLorentz = float(lineData.widthLorentz[0])
            lineStr = float(lineData.str[0])
        else:
            widthGauss = float(lineData.widthGauss[0])
            widthLorentz = float(lineData.widthLorentz[0])
            lineStr = float(lineData.str[0])

        widthGaussWl = widthGauss / _C_KMS * lineData.wl0[0]
        wlGaussNorm = (lineData.wl0[0] - self.wl) / widthGaussWl

        # Humlicek (1982) 分段 w4 计算
        z = widthLorentz - 1j * wlGaussNorm
        zz = z * z
        s = np.abs(wlGaussNorm) + widthLorentz
        w4 = np.zeros(self.wl.shape, dtype=complex)

        con1 = s >= 15.0
        zt = z[con1]
        w4[con1] = 0.56418958355 * zt / (0.5 + zt * zt)

        con2 = (s >= 5.5) & (s < 15.0)
        zt = z[con2]
        zzt = zz[con2]
        w4[con2] = zt * (1.4104739589 + 0.56418958355 * zzt) / (
            (3.0 + zzt) * zzt + 0.7499999999)

        con3 = (widthLorentz >= 0.195 * np.abs(wlGaussNorm) - 0.176) & (s <
                                                                        5.5)
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
        w4[con4] = (np.exp(zzt) - zt *
                    (36183.30536 - zzt *
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

        # Voigt 轮廓对波长的导数（用于 Stokes V 核）
        dVoigtdWl = (lineStr / widthGaussWl * 2.0) * (-wlGaussNorm * w4.real +
                                                      widthLorentz * w4.imag)

        self.Iunscaled = 1.0 - lineStr * w4.real
        # Stokes V 核：gCoeff = e/(4π m_e c²) × λ₀² × g_eff  (cgs, λ in nm)
        gCoeff = -4.6686e-12 * lineData.wl0[0]**2 * lineData.g[0]
        self.dIdWlG = gCoeff * dVoigtdWl

    def updateProfDeriv(self, lineData, Blos, dBlos_d, brightMap, numPts,
                        wlGrid, surfaceScale, calcDI: int,
                        calcDV: int) -> None:
        """按当前亮度图和磁场更新盘积分轮廓及导数。"""
        contin = surfaceScale * brightMap.bright
        self.I = contin[np.newaxis, :] * self.Iunscaled  # noqa: E741
        self.V = (Blos * contin)[np.newaxis, :] * self.dIdWlG

        self.Ic = np.sum(contin)
        invSumIc0 = 1.0 / np.sum(self.Ic)

        # 盘积分并按连续谱归一化
        self.I = np.sum(self.I, axis=1) * invSumIc0  # noqa: E741

        if calcDV == 1:
            dIdWlGcont = self.dIdWlG * contin[np.newaxis, :]
            self.dVsum = np.einsum('ijl,kl->ikj', dBlos_d, dIdWlGcont)
            self.dVsum = np.swapaxes(self.dVsum, 1, 2)
            self.dVsum *= invSumIc0
        else:
            self.dVsum = 0

        if calcDI == 1:
            self.dIcsum = (surfaceScale * invSumIc0)[:, np.newaxis] * (
                self.Iunscaled.T - self.I[np.newaxis, :])
        else:
            self.dIcsum = 0

        if calcDI == 1 and calcDV == 1:
            self.dVdBright = (surfaceScale * invSumIc0)[:, np.newaxis] * (
                (Blos[np.newaxis, :] * self.dIdWlG).T - self.V[np.newaxis, :])
        else:
            self.dVdBright = 0

    def dopplerShift(self, velShift: float) -> None:
        """对轮廓施加 Doppler 偏移（+为远离观测者）。"""
        self.wl = self.wl + self.wl * velShift / _C_KMS


# ---------------------------------------------------------------------------
# diskIntProfAndDeriv
# ---------------------------------------------------------------------------


class diskIntProfAndDeriv:
    """盘积分轮廓及其关于亮度图和磁场球谐系数的全部一阶导数。

    Parameters
    ----------
    visibleGrid : _SinglePhaseView
        单相位可见性几何视图。
    vMagCart : (3, N_cells) ndarray
        笛卡尔磁场分量。
    dMagCart : ndarray 或 0
        磁场对球谐系数的导数（``fitMag=0`` 时传 0）。
    brightMap
        ``brightMap`` 亮度图对象。
    lData : lineData
        谱线参数对象。
    velEq : float
        赤道线速度（km/s）。
    wlGrid : (N_vel,) ndarray
        合成波长网格（来自观测 LSD 格点）。
    calcDI, calcDV : int
        分别控制是否计算亮度导数和磁场导数（0=否，1=是）。
    """
    def __init__(self, visibleGrid, vMagCart: np.ndarray, dMagCart, brightMap,
                 lData, velEq: float, wlGrid: np.ndarray, calcDI: int,
                 calcDV: int) -> None:
        self.wl = wlGrid
        self.numPts = wlGrid.shape[0]
        self.wlStart = wlGrid[0]
        self.wlEnd = wlGrid[-1]

        self.IIc = np.zeros(self.numPts)
        self.QIc = np.zeros(self.numPts)
        self.UIc = np.zeros(self.numPts)
        self.VIc = np.zeros(self.numPts)

        # 各格点的赤道旋转 Doppler 速度
        vrCell = velEq * visibleGrid.velRotProj
        # 每个格点的波长网格（Dopppler 偏移的逆移位）
        self.wlCells = np.outer(self.wl, 1.0 / (1.0 + vrCell / _C_KMS))

        self.prof = localProfileAndDeriv(lData, self.numPts, self.wlCells)
        self.updateIntProfDeriv(visibleGrid, vMagCart, dMagCart, brightMap,
                                lData, calcDI, calcDV)

    def updateIntProfDeriv(self, visibleGrid, vMagCart: np.ndarray, dMagCart,
                           brightMap, lData, calcDI: int, calcDV: int) -> None:
        """按新几何/亮度/磁场重新计算盘积分轮廓及导数。"""
        self.calcDI = calcDI
        self.calcDV = calcDV

        vView = visibleGrid.vViewCart
        projBlos = np.einsum('ij,ij->j', vView, vMagCart)

        if calcDV == 1:
            dBlos_d = np.einsum('il,ijkl->jkl', vView, dMagCart)
        else:
            dBlos_d = 0

        limbD = limbDarkening(lData.limbDark, visibleGrid.viewAngle)
        gDark = visibleGrid.gravityDarkening(lData.gravDark)
        surfaceScale = (limbD * gDark * visibleGrid.projArea *
                        visibleGrid.visible)

        self.prof.updateProfDeriv(lData, projBlos, dBlos_d, brightMap,
                                  self.numPts, self.wlCells, surfaceScale,
                                  calcDI, calcDV)

        self.IIc = self.prof.I
        self.VIc = self.prof.V
        self.dVIc = self.prof.dVsum
        self.dIIc = self.prof.dIcsum
        self.dVdBri = self.prof.dVdBright

    def convolveIGnumpy(self, fwhm: float) -> None:
        """对盘积分轮廓施加 Gaussian 仪器展宽（分辨率 R = fwhm）。

        若 ``explicitConvolution=False``（默认），则直接返回无操作。
        假设输入波长网格均匀间隔。
        """
        if not explicitConvolution:
            return

        wlCenter = (self.wlStart + self.wlEnd) / 2.0
        wlGstep = (self.wlEnd - self.wlStart) / (self.numPts - 1)
        fwhmWl = wlCenter / fwhm
        wlEndG = 3.0 * fwhmWl
        numPtsG = 2 * int(wlEndG / wlGstep) + 1

        if fwhmWl < wlGstep:
            return

        # 检查是否均匀间隔
        steps = self.wl[1:-1] - self.wl[:-2]
        if not np.allclose(steps, wlGstep, rtol=1e-3):
            print(
                'Error: uneven spacing in wavelength pixels. '
                'Convolution with the instrumental profile may produce errors.'
            )

        wlG = np.linspace(-wlEndG, wlEndG, numPtsG)
        profG = (0.939437 / fwhmWl) * np.exp(-2.772589 * (wlG / fwhmWl)**2)
        profG /= np.sum(profG)
        pad = numPtsG // 2

        def _pad_and_convolve(arr: np.ndarray) -> np.ndarray:
            padded = np.concatenate(
                [np.repeat(arr[:1], pad), arr,
                 np.repeat(arr[-1:], pad)])
            return np.convolve(padded, profG, mode='valid')

        self.IIc = _pad_and_convolve(self.IIc)
        self.QIc = _pad_and_convolve(self.QIc)
        self.UIc = _pad_and_convolve(self.UIc)
        self.VIc = _pad_and_convolve(self.VIc)

        if self.calcDI == 1:
            tmpdI = np.concatenate([
                np.tile(self.dIIc[:, :1], (1, pad)),
                self.dIIc,
                np.tile(self.dIIc[:, -1:], (1, pad)),
            ],
                                   axis=1)
            flat = np.convolve(tmpdI.ravel(), profG, mode='same')
            self.dIIc = flat.reshape(tmpdI.shape)[:, pad:-pad]

        if self.calcDV == 1:
            tmpdV = np.concatenate([
                np.tile(self.dVIc[:, :, :1], (1, 1, pad)),
                self.dVIc,
                np.tile(self.dVIc[:, :, -1:], (1, 1, pad)),
            ],
                                   axis=2)
            flat = np.convolve(tmpdV.ravel(), profG, mode='same')
            self.dVIc = flat.reshape(tmpdV.shape)[:, :, pad:-pad]

        if self.calcDI == 1 and self.calcDV == 1:
            tmpdV = np.concatenate([
                np.tile(self.dVdBri[:, :1], (1, pad)),
                self.dVdBri,
                np.tile(self.dVdBri[:, -1:], (1, pad)),
            ],
                                   axis=1)
            flat = np.convolve(tmpdV.ravel(), profG, mode='same')
            self.dVdBri = flat.reshape(tmpdV.shape)[:, pad:-pad]


# ---------------------------------------------------------------------------
# 批量初始化辅助
# ---------------------------------------------------------------------------


def getAllProfDiriv(par, listGridView, vecMagCart: np.ndarray, dMagCart0,
                    briMap, lineData, wlSynSet: list) -> list:
    """为各观测相位构造初始化 ``diskIntProfAndDeriv`` 列表。

    Parameters
    ----------
    par
        参数对象（需有 ``cycleList``, ``velEq``, ``instrumentRes``）。
    listGridView
        BatchVisibleGrid 或列表，按索引返回单相位视图。
    vecMagCart : (3, N_cells) ndarray
        初始磁场笛卡尔向量。
    dMagCart0
        磁场导数（``fitMag=0`` 时传 0）。
    briMap
        亮度图对象。
    lineData
        谱线参数对象。
    wlSynSet : list
        各相位合成波长网格列表。

    Returns
    -------
    list[diskIntProfAndDeriv]
        各相位的盘积分轮廓对象。
    """
    setSynSpec = []
    for nObs, _phase in enumerate(par.cycleList):
        spec = diskIntProfAndDeriv(listGridView[nObs], vecMagCart, dMagCart0,
                                   briMap, lineData, par.velEq, wlSynSet[nObs],
                                   0, par.calcDV)
        spec.convolveIGnumpy(par.instrumentRes)
        setSynSpec.append(spec)
    return setSynSpec
