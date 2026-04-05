"""core.geometry.stellar_grid — 恒星表面格点 starGrid 类。

starGrid 将恒星表面离散为近等面积的球坐标 (余纬, 经度) 网格元，
支持球形和扁球形（Roche 模型）两种几何。

主要功能
--------
- 生成近等面积网格，格点数由 numPoints（余纬圈数）决定。
- 支持扁球形快速旋转星（Roche 面）：半径 r(clat)、面积积分、表面法向量。
- 预计算 Cartesian 法向量、旋转速度向量和微分自转缓存，
  供 BatchVisibleGrid 高效调用。
- gravityDarkening：逐格点重力昏暗因子（仅扁球形有效）。
"""
from __future__ import annotations

from math import pi, acos

import numpy as np

from core.geometry.differential_rotation import calcVelDiffrotFactor


class starGrid:
    """恒星表面格点对象。

    Parameters
    ----------
    numPoints : int
        余纬方向的格点（环）数。总格点数由此决定，约为 2*numPoints^2。
    period : float, optional
        自转周期（天）。用于扁球形的 Roche 半径计算。
    mass : float, optional
        恒星质量（太阳质量）。> 0 时启用扁球形几何。
    radiusEq : float, optional
        赤道半径（太阳半径）。> 0 时启用扁球形几何。
    verbose : int, optional
        0 = 静默；1 = 打印网格信息和扁度信息（默认 1）。
    """
    def __init__(self,
                 numPoints: int,
                 period: float = 0.0,
                 mass: float = 0.0,
                 radiusEq: float = 0.0,
                 verbose: int = 1) -> None:
        # 余纬方向格点数
        self.numPtsClat = numPoints
        # 赤道处经度格点数（约 2 倍余纬格点数）
        self.numPtsLongEq = 2 * self.numPtsClat
        # 各余纬环的经度格点数
        self.numPtsLong_Clat = np.zeros(self.numPtsClat, dtype=int)

        # 计算总格点数（近等面积格点）
        self.numPoints = 0
        for i in range(self.numPtsClat):
            _clat = pi * ((float(i) + 0.5) / float(self.numPtsClat))
            lnumPtsLong = int(round(np.sin(_clat) * float(self.numPtsLongEq)))
            self.numPtsLong_Clat[i] = lnumPtsLong
            self.numPoints += lnumPtsLong

        self.clat = np.zeros(self.numPoints)
        self.long = np.zeros(self.numPoints)
        self.radius = np.zeros(self.numPoints)
        self.dClat = np.zeros(self.numPoints)
        self.dLong = np.zeros(self.numPoints)
        self.dRadius = np.zeros(self.numPoints)
        self.gdark_local = None

        # 建立余纬-经度网格
        n = 0
        for i in range(self.numPtsClat):
            _clat = pi * ((float(i) + 0.5) / float(self.numPtsClat))
            _dClat = pi * (1.0 / float(self.numPtsClat))
            lnumPtsLong = int(self.numPtsLong_Clat[i])
            for j in range(lnumPtsLong):
                self.clat[n] = _clat
                self.dClat[n] = _dClat
                _long = 2.0 * pi * ((float(j) + 0.5) / float(lnumPtsLong))
                self.long[n] = _long
                _dLong = 2.0 * pi * (1.0 / float(lnumPtsLong))
                self.dLong[n] = _dLong
                n += 1

        # 默认球形几何：r = 1 everywhere
        self.radius[:] = 1.0
        self.dRadius[:] = 0.0
        self.fracOmegaCrit = 0.0

        # 如给出质量和半径，计算扁球形 Roche 模型
        if mass > 0.0 and radiusEq > 0.0:
            self.radius, self.fracOmegaCrit = self.calcRclatOblate(
                period, mass, radiusEq, verbose)
            self.dRadius = self.getOblateRadiusStep(period, mass, radiusEq)

        # 格点面积
        self.area = self.GetSurfaceArea()

        # 预计算相位无关量，供 BatchVisibleGrid 复用
        self._normals: np.ndarray = self.getCartesianNormals()  # (3, N_cells)
        self._rot_vel: np.ndarray = self.GetCartesianRotVel()  # (3, N_cells)
        self._diffrot_cache_key: tuple | None = None
        self._diffrot_cache: np.ndarray | None = None

        if verbose == 1:
            print("initiated a {:} point stellar grid  "
                  "( {:} clat, {:} long_equator)".format(
                      self.numPoints, self.numPtsClat, self.numPtsLongEq))

    # ------------------------------------------------------------------
    # 格点几何辅助方法
    # ------------------------------------------------------------------

    def GetCellCorners(self, i: int) -> np.ndarray:
        """返回格点 i 的 4 个角点坐标 (long, clat, radius)，形状 (4, 3)。"""
        long1 = self.long[i] - 0.5 * self.dLong[i]
        long2 = self.long[i] + 0.5 * self.dLong[i]
        clat1 = self.clat[i] - 0.5 * self.dClat[i]
        clat2 = self.clat[i] + 0.5 * self.dClat[i]
        r1 = self.radius[i] - 0.5 * self.dRadius[i]
        r2 = self.radius[i] + 0.5 * self.dRadius[i]
        return np.array([[long1, clat1, r1], [long2, clat1, r1],
                         [long1, clat2, r2], [long2, clat2, r2]])

    def GetCartesianCells(self) -> np.ndarray:
        """返回所有格点中心的笛卡尔坐标，形状 (3, N_cells)。"""
        x = self.radius * np.sin(self.clat) * np.cos(self.long)
        y = self.radius * np.sin(self.clat) * np.sin(self.long)
        z = self.radius * np.cos(self.clat)
        return np.array([x, y, z])

    def GetCartesianCellCorners(self, i: int) -> np.ndarray:
        """返回格点 i 的 4 个角点笛卡尔坐标，形状 (4, 3)。"""
        cartCor = np.zeros((4, 3))
        cor = self.GetCellCorners(i)
        for j in range(4):
            x = cor[j, 2] * np.sin(cor[j, 1]) * np.cos(cor[j, 0])
            y = cor[j, 2] * np.sin(cor[j, 1]) * np.sin(cor[j, 0])
            z = cor[j, 2] * np.cos(cor[j, 1])
            cartCor[j, 0] = x
            cartCor[j, 1] = y
            cartCor[j, 2] = z
        return cartCor

    def GetSurfaceArea(self) -> np.ndarray:
        """计算各格点的表面积，返回 (N_cells,) 数组。

        球形星采用解析精确公式；扁球形星用数值积分（scipy.integrate.quad）。
        """
        area = np.zeros(self.numPoints)

        if self.fracOmegaCrit > 0.0:
            # 扁球形：通过对势函数梯度构造面法向，再做数值积分
            def fThetaP(theta: float, wo: float) -> float:
                dS = (9.0 / wo**2 / np.sin(theta) * np.cos(
                    (np.pi + np.arccos(wo * np.sin(theta))) / 3.0)**2)
                rs = (3.0 / (wo * np.sin(theta)) * np.cos(
                    (np.pi + np.arccos(wo * np.sin(theta))) / 3.0))
                gradPhi_r = 1.0 / rs**2 - 8.0 / 27.0 * rs * (wo *
                                                             np.sin(theta))**2
                gradPhi_theta = -8.0 / 27.0 * rs * wo**2 * np.sin(
                    theta) * np.cos(theta)
                cosGamma = gradPhi_r / np.sqrt(gradPhi_r**2 + gradPhi_theta**2)
                return dS / cosGamma

            from scipy.integrate import quad
            for i in range(self.numPoints):
                if self.clat[i] != self.clat[i - 1]:
                    intTheta, _ = quad(fThetaP,
                                       self.clat[i] - 0.5 * self.dClat[i],
                                       self.clat[i] + 0.5 * self.dClat[i],
                                       args=(self.fracOmegaCrit, ))
                    area[i] = self.dLong[i] * intTheta
        else:
            # 球形：精确解析面积
            area = self.dLong * (np.cos(self.clat - 0.5 * self.dClat) -
                                 np.cos(self.clat + 0.5 * self.dClat))
        return area

    def GetCartesianRotVel(self) -> np.ndarray:
        """返回所有格点以赤道线速度为单位的旋转速度向量，形状 (3, N_cells)。

        假设固体旋转（solid body rotation）。
        """
        velX = -self.radius * np.sin(self.clat) * np.sin(self.long)
        velY = self.radius * np.sin(self.clat) * np.cos(self.long)
        velZ = np.zeros(len(self.clat))
        velVec = np.array([velX, velY, velZ])
        # 对扁球形星按赤道半径归一化
        iEquator = np.argmin(np.abs(self.clat - 0.5 * np.pi))
        velVec /= self.radius[iEquator]
        return velVec

    def getCartesianNormals(self) -> np.ndarray:
        """计算所有格点的单位外法向量（笛卡尔），形状 (3, N_cells)。

        球形星：法向量即归一化半径向量。
        扁球形星：使用等势面梯度解析表达式。
        """
        if self.fracOmegaCrit <= 0.0:
            vNormals = self.GetCartesianCells()
            lenThisPoint = np.sqrt(vNormals[0, :]**2 + vNormals[1, :]**2 +
                                   vNormals[2, :]**2)
            vNormals /= lenThisPoint
        else:
            wo = self.fracOmegaCrit
            sr = self.radius
            sclat = self.clat
            slong = self.long
            # 利用势函数梯度得到表面法向量（解析形式）
            gradPot = np.array([
                -8.0 / 27.0 * wo**2 * sr * np.sin(sclat)**2 + 1.0 / sr**2,
                -8.0 / 27.0 * wo**2 * sr * np.sin(sclat) * np.cos(sclat),
                np.zeros_like(sclat),
            ])
            gpotnorm = gradPot / np.sqrt(gradPot[0]**2 + gradPot[1]**2)
            # 球坐标 → 笛卡尔旋转矩阵
            R = np.array([
                [
                    np.sin(sclat) * np.cos(slong),
                    np.cos(sclat) * np.cos(slong), -np.sin(slong)
                ],
                [
                    np.sin(sclat) * np.sin(slong),
                    np.cos(sclat) * np.sin(slong),
                    np.cos(slong)
                ],
                [np.cos(sclat), -np.sin(sclat),
                 np.zeros_like(sclat)],
            ])
            vNormals = np.einsum('ijk,jk->ik', R, gpotnorm)
        return vNormals

    # ------------------------------------------------------------------
    # 扁球形 Roche 几何
    # ------------------------------------------------------------------

    def calcRclatOblate(self,
                        period: float,
                        mass: float,
                        radius: float,
                        verbose: int = 1) -> tuple[np.ndarray, float]:
        """按 Roche 模型计算各余纬处的归一化半径 r(clat)。

        参考 Tassoul 1978；Collins 1963, 1965, 1966。
        返回 (x_clat, wo)：归一化半径数组和临界角速度比率 Omega/Omega_crit。
        """
        G = 6.67408e-11  # 引力常数 (m^3 kg^-1 s^-2)
        Msun = 1.98892e30  # 太阳质量 (kg)
        Rsun = 6.955e8  # 太阳半径 (m)
        Omega = 2.0 * np.pi / (period * 86400.0)
        Omega_break = np.sqrt(
            (8.0 / 27.0) * G * (mass * Msun) / (radius * Rsun)**3)
        wo = Omega / Omega_break

        if wo > 1.0:
            print("Error: star above breakup velocity "
                  "(Omega/Omega_break {:}),\n"
                  "   for given mass {:} Msun, radius {:} Rsun, "
                  "and period {:} d".format(wo, mass, radius, period))
            print("Assuming breakup velocity for geometry calculations")
            wo = 1.0
        elif wo < 1e-8:
            wo = 1e-8

        obl = (3.0 / wo) * np.cos((np.pi + acos(wo)) / 3.0)
        if verbose:
            print("The Oblateness of the star: "
                  "Req/Rp = {:} ({:} of breakup)".format(obl, wo))

        x_clat = (3.0 / (wo * np.sin(self.clat)) * np.cos(
            (np.pi + np.arccos(wo * np.sin(self.clat))) / 3.0))
        return x_clat, wo

    def getOblateRadiusStep(self, period: float, mass: float,
                            radiusEq: float) -> np.ndarray:
        """计算扁球形格点角点间的半径跨度 dR(clat)，用于面积积分。"""
        if mass > 0.0 and radiusEq > 0.0:
            _refClat = self.clat
            self.clat = _refClat + 0.499999999 * self.dClat
            radiusP, _ = self.calcRclatOblate(period,
                                              mass,
                                              radiusEq,
                                              verbose=False)
            self.clat = _refClat - 0.499999999 * self.dClat
            radiusM, _ = self.calcRclatOblate(period,
                                              mass,
                                              radiusEq,
                                              verbose=False)
            self.clat = _refClat  # 恢复正确中心值
            dRadius = radiusP - radiusM
        else:
            dRadius = np.zeros_like(self.clat)
        return dRadius

    # ------------------------------------------------------------------
    # 差分自转和重力昏暗
    # ------------------------------------------------------------------

    def get_diffrot_factor(self, period: float, dOmega: float) -> np.ndarray:
        """返回归一化差分自转速度比例因子 (N_cells,)。

        结果按 (period, dOmega) 键缓存；键变化时自动重算。
        """
        key = (period, dOmega)
        if self._diffrot_cache_key != key:
            self._diffrot_cache_key = key
            self._diffrot_cache = calcVelDiffrotFactor(period, dOmega,
                                                       self.clat)
        return self._diffrot_cache  # type: ignore[return-value]

    def gravityDarkening(self, gravDarkCoeff: float) -> np.ndarray:
        """计算逐格点重力昏暗因子（仅扁球形有效，球形星返回全 1）。

        结果按 gravDarkCoeff 缓存。
        """
        bCalc = True
        if self.gdark_local is not None and self.gravDarkCoeff == gravDarkCoeff:
            bCalc = False

        if bCalc:
            wo = self.fracOmegaCrit
            if wo > 0.0 and gravDarkCoeff > 0.0:
                gr = (-1.0 / self.radius**2 + 8.0 / 27.0 * self.radius *
                      (wo * np.sin(self.clat))**2)
                gtheta = (8.0 / 27.0 * self.radius * np.sin(self.clat) *
                          np.cos(self.clat) * wo**2)
                gdark_local = np.sqrt(gr**2 + gtheta**2)**gravDarkCoeff
            else:
                gdark_local = np.ones_like(self.clat)
            self.gravDarkCoeff = gravDarkCoeff
            self.gdark_local = gdark_local

        return self.gdark_local

    # ------------------------------------------------------------------
    # 辅助（其他用途）
    # ------------------------------------------------------------------

    def GetDistRotAxis(self) -> np.ndarray:
        """返回各格点到自转轴的距离（当前未使用）。"""
        return self.radius * np.sin(self.clat)
