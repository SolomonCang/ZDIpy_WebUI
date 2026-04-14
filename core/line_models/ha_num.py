"""core.line_models.ha_num — H-alpha 数值模板缩放模型

直接使用数值Hα发射线I(λ)轮廓模板，通过缩放和弱场近似线性响应反演Stokes V。

接口与halpha.py兼容：
- lineDataHaNum
- localProfileAndDerivHaNum
- diskIntProfAndDerivHaNum
- getAllProfDirivHaNum

"""

import numpy as np


class lineDataHaNum:
    """H-alpha 数值模板参数容器。"""

    def __init__(self,
                 wl_template,
                 I_template,
                 g=1.048,
                 fV=1.0,
                 limb_darkening=0.0,
                 gravity_darkening=0.0):
        self.wl_template = np.array(wl_template)  # (N_wl,)
        self.I_template = np.array(I_template)  # (N_wl,)
        self.g = float(g)
        self.fV = float(fV)
        self.limb_darkening = float(limb_darkening)
        self.gravity_darkening = float(gravity_darkening)
        self.numLines = 1
        self.instRes = -1.0
        # --- pipeline 兼容属性（对齐 lineData / lineDataHalpha 接口） ---
        wl_center = float(self.wl_template[len(self.wl_template) // 2])
        self.wl0 = np.array([wl_center])
        self.limbDark = np.array([limb_darkening])
        self.gravDark = np.array([gravity_darkening])
        # 等效线强（供 pipeline metadata 用）：模板轮廓最大吸收深度
        self.str = np.array([float(1.0 - np.min(self.I_template))])

    @classmethod
    def from_arrays(cls,
                    wl_template,
                    I_template,
                    g=1.048,
                    fV=1.0,
                    limb_darkening=0.0,
                    gravity_darkening=0.0):
        return cls(wl_template, I_template, g, fV, limb_darkening,
                   gravity_darkening)


class localProfileAndDerivHaNum:
    """预计算数值Hα轮廓及Stokes V核。"""
    _C_KMS = 2.99792458e5
    _ZEEMAN_CONST = 4.66864e-12

    def __init__(self, ldata: lineDataHaNum, numPts: int, wl_grid: np.ndarray):
        self.wl = wl_grid  # (N_vel, N_cells)
        # 插值I(λ)
        I_interp = np.interp(wl_grid, ldata.wl_template, ldata.I_template)
        self.Iunscaled = I_interp  # (N_vel, N_cells)
        # 数值导数：wl_grid 为 2D，用链式法则 dI/dλ = (dI/dk)/(dλ/dk)
        dI_dlambda = np.gradient(I_interp, axis=0) / np.gradient(wl_grid,
                                                                 axis=0)
        gCoeff = -self._ZEEMAN_CONST * ldata.wl_template[len(ldata.wl_template)
                                                         // 2]**2 * ldata.g
        self.VsigKernel = ldata.fV * gCoeff * dI_dlambda

    def updateProfDeriv(self, ldata, Blos, dBlos_d, brightMap, numPts, wl_grid,
                        surfaceScale, calcDI, calcDV):
        contin = surfaceScale * brightMap.bright
        self.Ic = np.sum(contin)
        invSumIc0 = 1.0 / self.Ic
        V_cell = (Blos * contin)[np.newaxis, :] * self.VsigKernel
        # 避免歧义，改为 IIc
        self.IIc = np.sum(contin[np.newaxis, :] * self.Iunscaled,
                          axis=1) * invSumIc0
        self.V = np.sum(V_cell, axis=1) * invSumIc0
        if calcDV == 1:
            VsigKernelCont = self.VsigKernel * contin[np.newaxis, :]
            self.dVsum = np.einsum('ijl,kl->ikj', dBlos_d, VsigKernelCont)
            self.dVsum = np.swapaxes(self.dVsum, 1, 2) * invSumIc0
        else:
            self.dVsum = 0
        if calcDI == 1:
            self.dIcsum = (surfaceScale * invSumIc0)[:, np.newaxis] * (
                self.Iunscaled.T - self.IIc[np.newaxis, :])
        else:
            self.dIcsum = 0
        if calcDI == 1 and calcDV == 1:
            V_unscaled = Blos[np.newaxis, :] * self.VsigKernel
            self.dVdBright = (surfaceScale * invSumIc0)[:, np.newaxis] * (
                V_unscaled.T - self.V[np.newaxis, :])
        else:
            self.dVdBright = 0

    def dopplerShift(self, vel_shift: float):
        self.wl = self.wl + self.wl * vel_shift / self._C_KMS


class diskIntProfAndDerivHaNum:
    """H-alpha 数值模板盘积分轮廓及全部一阶导数。"""

    def __init__(self, visible_grid, v_mag_cart, d_mag_cart, bright_map, ldata,
                 vel_eq, wl_grid, calcDI, calcDV):
        self.wl = wl_grid
        self.numPts = wl_grid.shape[0]
        self.wlStart = wl_grid[0]
        self.wlEnd = wl_grid[-1]
        self.IIc = np.zeros(self.numPts)
        self.QIc = np.zeros(self.numPts)
        self.UIc = np.zeros(self.numPts)
        self.VIc = np.zeros(self.numPts)
        vr_cell = vel_eq * visible_grid.velRotProj
        self.wlCells = np.outer(
            self.wl, 1.0 / (1.0 + vr_cell / localProfileAndDerivHaNum._C_KMS))
        self.prof = localProfileAndDerivHaNum(ldata, self.numPts, self.wlCells)
        self.updateIntProfDeriv(visible_grid, v_mag_cart, d_mag_cart,
                                bright_map, ldata, calcDI, calcDV)

    def updateIntProfDeriv(self, visible_grid, v_mag_cart, d_mag_cart,
                           bright_map, ldata, calcDI, calcDV):
        self.calcDI = calcDI
        self.calcDV = calcDV
        Blos, dBlos_d = self._BlosProjected(visible_grid.vViewCart, v_mag_cart,
                                            d_mag_cart)
        from core.line_models.line_utils import limbDarkening
        limb_d = limbDarkening(ldata.limb_darkening, visible_grid.viewAngle)
        g_dark = visible_grid.gravityDarkening(ldata.gravity_darkening)
        surface_scale = limb_d * g_dark * visible_grid.projArea * visible_grid.visible
        self.prof.updateProfDeriv(ldata, Blos, dBlos_d, bright_map,
                                  self.numPts, self.wlCells, surface_scale,
                                  calcDI, calcDV)
        self.IIc = self.prof.IIc
        self.VIc = self.prof.V
        self.dVIc = self.prof.dVsum
        self.dIIc = self.prof.dIcsum
        self.dImag = 0
        self.dVdBri = self.prof.dVdBright

    def _BlosProjected(self, v_view, v_mag_cart, d_mag_cart):
        Blos = np.einsum('ij,ij->j', v_view, v_mag_cart)
        if self.calcDV == 1:
            dBlos_d = np.einsum('il,ijkl->jkl',
                                v_view,
                                d_mag_cart,
                                optimize=True)
        else:
            dBlos_d = 0
        return Blos, dBlos_d

    def convolveIGnumpy(self, fwhm: float):
        if fwhm <= 0.0:
            return
        wl_center = (self.wlStart + self.wlEnd) / 2.0
        wl_step = (self.wlEnd - self.wlStart) / (self.numPts - 1)
        fwhm_wl = wl_center / fwhm
        wl_end_g = 3.0 * fwhm_wl
        num_pts_g = 2 * int(wl_end_g / wl_step) + 1
        if fwhm_wl < wl_step:
            return
        wl_g = np.linspace(-wl_end_g, wl_end_g, num_pts_g)
        prof_g = (0.939437 / fwhm_wl) * np.exp(-2.772589 * (wl_g / fwhm_wl)**2)
        prof_g /= np.sum(prof_g)
        pad = num_pts_g // 2

        def _pad_conv(arr):
            padded = np.concatenate(
                [np.repeat(arr[:1], pad), arr,
                 np.repeat(arr[-1:], pad)])
            return np.convolve(padded, prof_g, mode="valid")

        self.IIc = _pad_conv(self.IIc)
        self.VIc = _pad_conv(self.VIc)
        # Convolve each row (per-cell velocity profile) independently to avoid
        # cross-cell contamination from a naive ravel+convolve on a 2-D array.
        if self.calcDI == 1:
            self.dIIc = np.asarray([_pad_conv(row) for row in self.dIIc])
        if self.calcDV == 1:
            # dVIc shape: (n_comp, n_cells, n_vel)
            self.dVIc = np.asarray([[_pad_conv(row) for row in comp]
                                    for comp in self.dVIc])
        if self.calcDI == 1 and self.calcDV == 1:
            self.dVdBri = np.asarray([_pad_conv(row) for row in self.dVdBri])


def getAllProfDirivHaNum(par, list_grid_view, vec_mag_cart, d_mag_cart0,
                         bri_map, ldata, wl_syn_set):
    set_syn_spec = []
    for nobs, _phase in enumerate(par.cycleList):
        spec = diskIntProfAndDerivHaNum(
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
        set_syn_spec.append(spec)
    return set_syn_spec
