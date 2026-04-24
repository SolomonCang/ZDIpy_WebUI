"""core.brightnessGeom — 恒星表面亮度图处理。

提供目标：
- brightMap     : 像素亮度图（limb/重力暗化加权投影）
- saveMap       : 将亮度图写入文件
- readMap       : 从文件读取亮度图（吼坐标对比校验）
- SetupBrightMap: 初始化亮度图（从文件或默认常数）
"""

import numpy as np


# Saves the brightness map as a set of pixels,
# using the grid of colatitude and longitude supplied at initialization.
class brightMap:

    def __init__(self, clat, lon):
        self.clat = clat
        self.lon = lon
        # assume we initialize to a uniform (normalized) brightness of one
        self.npt = clat.shape[0]
        self.bright = np.ones(self.npt)

    def makeRoundSpot(
        self,
        clat: float,
        lon: float,
        radius: float,
        brigthness: float,
    ) -> None:
        """Place a round spot on the brightness map (vectorised, no Python loop)."""
        dist = np.arccos(
            np.sin(self.clat) * np.sin(clat) * np.cos(self.lon - lon) +
            np.cos(self.clat) * np.cos(clat))
        self.bright = np.where(dist < radius, brigthness, self.bright)

    def projected(
            self,
            visible_batch,
            limb_dark_coeff: float,
            grav_dark_factor: np.ndarray,  # (N_cells,)
    ) -> np.ndarray:
        """Return visibility/limb/gravity-darkening weighted brightness.

        Parameters
        ----------
        visible_batch : BatchVisibleGrid
            Batch geometry for all phases.
        limb_dark_coeff : float
            Linear limb-darkening coefficient (epsilon in I(mu)=1-eps*(1-mu)).
        grav_dark_factor : (N_cells,)
            Gravity darkening factor from starGrid.gravityDarkening().

        Returns
        -------
        weighted : np.ndarray, shape (N_phases, N_cells)
            bright * limb_weight * grav_dark * proj_area, zero for hidden cells.
        """
        mu = np.cos(visible_batch.view_angle)  # (N_phases, N_cells)
        limb_weight = 1.0 - limb_dark_coeff * (1.0 - mu)  # (N_phases, N_cells)
        weighted = (
            self.bright[np.newaxis, :]  # (1, N_cells)
            * limb_weight  # (N_phases, N_cells)
            * grav_dark_factor[np.newaxis, :]  # (1, N_cells)
            * visible_batch.proj_area  # (N_phases, N_cells)
        )
        return np.where(visible_batch.visible, weighted, 0.0)


# Save a brightness map (as in the class above), to a specified file name
def saveMap(brightMap, fileName):
    with open(fileName, 'w') as outFile:
        for i in range(brightMap.npt):
            outFile.write(
                '%8.6f %8.6f %11.8f\n' %
                (brightMap.clat[i], brightMap.lon[i], brightMap.bright[i]))


# Read a brightness map in from a file (same format as saveMap())
# Check the input map's stellar grid coordinates against the grid coordinates the user supplies.
# If the coordinate systems are not consistent return an error.
def readMap(fileName, userClat, userLon, verbose=1):

    clat, lon, bright = np.loadtxt(fileName, unpack=True)

    if ((clat.shape[0] != userClat.shape[0]) |
        (lon.shape[0] != userLon.shape[0])):
        print(
            'ERROR reading map {:} requested coordinate grid shape does not match input coordinate grid shape'  # noqa: E501
            .format(fileName))
        return

    # small should be set to account for the precision of the coordinates used in the saved map.
    small = 1e-5
    if ((np.abs(clat - userClat) < small).all() &
        (np.abs(lon - userLon) < small).all()):
        bMap = brightMap(clat, lon)
        bMap.bright = bright
        if (verbose == 1):
            print('Initialized brightness map from {:}'.format(fileName))
        return bMap
    else:
        print(
            'ERROR reading map {:} requested coordinate grid does not match input coordinate grid'
            .format(fileName))
        return


def SetupBrightMap(sGrid,
                   initBrightFromFile,
                   initBrightFile,
                   defaultBright,
                   verbose=1):
    # initialize the brightness map from a file, or set to the 'default brightness'
    if (initBrightFromFile == 1):
        briMap = readMap(initBrightFile, sGrid.clat, sGrid.long, verbose)
    else:
        briMap = brightMap(sGrid.clat, sGrid.long)
        briMap.bright[:] = defaultBright

    return briMap


def epoch_brightness_scale(cq: np.ndarray,
                           phase: float,
                           var: int,
                           t: float = 100.0,
                           kk: float = 10.0,
                           dphi: float = 0.0,
                           dT: float = 1.0) -> tuple:
    """
    Apply phase-dependent brightness modulation to grid-cell brightness (Cq).

    This implements the epoch-dependent temperature variation algorithm from CTTSzdi2,
    where brightness at each grid cell is modulated by a phase-dependent envelope function.
    Cooler regions (spots) have reduced brightness; hotter regions have increased brightness.
    dT ∝ -dCq (temperature decrease = brightness decrease).

    Parameters
    ----------
    cq : np.ndarray
        Base brightness for each grid cell, shape (N_cells,).
    phase : float
        Current rotational phase (cycles).
    var : int
        Variation model selector:
        - 1: Gaussian profile (smooth peak with exponential tails)
        - 2: Cosine² profile (even smoother, bounded)
        - 3: Tanh profile (prevents unphysical brightness bounds)
        - 4: Exponential profile (rapid evolution)
    t : float, optional
        Typical spot lifetime in rotation cycles (default 100.0).
    kk : float, optional
        Scaling parameter for phase-to-time conversion (default 10.0).
    dphi : float, optional
        Phase offset: center of feature at dphi (default 0.0).
    dT : float, optional
        Amplitude scale of epoch variation. 1.0 keeps the original modulation,
        0.0 disables modulation, values >1 amplify modulation (default 1.0).

    Returns
    -------
    ccq : np.ndarray
        Modified brightness after epoch variation, shape (N_cells,).
    dcq : np.ndarray
        Partial derivative ∂ccq/∂cq for Jacobian computation, shape (N_cells,).
    """
    # Normalized phase offset
    s = t / (np.sqrt(np.log(2.0)) * 2.0)
    dd = (phase - dphi) / s if s > 0 else 0.0

    def _apply_dt_scale(ccq_base: np.ndarray, dcq_base: np.ndarray) -> tuple:
        # Scale the modulation amplitude while preserving backward compatibility at dT=1.
        ccq = cq + dT * (ccq_base - cq)
        dcq = 1.0 + dT * (dcq_base - 1.0)
        return ccq, dcq

    if var == 1:
        # Gaussian profile: f(dd) = exp(-dd²)
        f = np.exp(-(dd**2))
        ccq_base = 1.0 - (1.0 - cq) * f
        dcq_base = np.full_like(cq, f, dtype=np.float64)
        return _apply_dt_scale(ccq_base, dcq_base)

    elif var == 2:
        # Cosine² profile: f(dd) = cos²(π·dd/2) for |dd| ≤ 1.0, else 0
        if abs(dd) <= 1.0:
            f = np.cos(np.pi * dd / 2.0)**2
        else:
            f = 0.0
        ccq_base = 1.0 - (1.0 - cq) * f
        dcq_base = np.full_like(cq, f, dtype=np.float64)
        return _apply_dt_scale(ccq_base, dcq_base)

    elif var == 3:
        # Tanh profile: smooth phase transition, prevents negative brightness
        ph0 = kk * np.log(kk)
        # Handle edge case when cq ≈ 1.0 (xx ≈ 0, division by zero)
        safe_xx = np.where(
            np.abs(1.0 / cq - 1.0) < 1e-10, 1e-10, 1.0 / cq - 1.0)
        zz = np.tanh(dd * ph0 * safe_xx) / safe_xx
        fac = 1.0 + zz
        ccq_base = 1.0 - (1.0 - cq) * fac
        dcq_base = np.full_like(cq, 1.0,
                                dtype=np.float64)  # Conservative Jacobian
        return _apply_dt_scale(ccq_base, dcq_base)

    elif var == 4:
        # Exponential profile: multiplicative scaling
        ph0 = kk * np.log(kk)
        fac = np.exp(dd * ph0)
        ccq_base = cq * fac
        dcq_base = np.full_like(cq, fac, dtype=np.float64)
        return _apply_dt_scale(ccq_base, dcq_base)

    else:
        # Default: no variation
        return cq.copy(), np.ones_like(cq, dtype=np.float64)
