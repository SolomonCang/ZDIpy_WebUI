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
