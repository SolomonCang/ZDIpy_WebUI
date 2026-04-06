"""core.magneticGeom — 球谐展开磁场计算。

提供：
- magSphHarmonics      : 球谐系数对象（含预计算 Legendre 多项式批量向量化）
- SetupMagSphHarmoics  : 初始化磁场（从文件或零常数）
- magSphHarmoicsFromFile: 从 JF Donati 格式实系数文件读入
- magSphHarmoics       : magSphHarmonics 的向后兼容别名
- dLegendre_dTheta     : 单点标量 Legendre 导数（备用）
"""
# Notes:
# Spherical harmonic calculations of (X,Y,Z) moved to an initialization step.
# Batch-vectorised associated Legendre polynomials use scipy.special.lpmn.

import numpy as np
import scipy.special
try:
    from scipy.special import factorial
except ImportError:  # For some older scipy versions
    from scipy.misc import factorial


def _compute_legendre_batch(
        m_max: int,
        n_max: int,
        x: np.ndarray,  # (N_cells,)  cos(colatitude)
) -> tuple[np.ndarray, np.ndarray]:
    """Batch-vectorised associated Legendre polynomials and derivatives.

     Returns
     -------
     p  : (m_max+1, n_max+1, N_cells)   P_n^m(x_i)
     pd : (m_max+1, n_max+1, N_cells)   dP_n^m/dx evaluated at x_i

     The outer loops are over (m, n) harmonic indices — O(n_harmonics^2)
     iterations total, each vectorised over N_cells — rather than the old
     O(N_cells) Python loop with O(n_harmonics^2) work per cell.
     """
    N = len(x)
    p = np.zeros((m_max + 1, n_max + 1, N))
    pd = np.zeros((m_max + 1, n_max + 1, N))

    # Polynomial values: O(n_harmonics^2) vectorised lpmv calls.
    # lpmv includes the Condon-Shortley phase (-1)^m, but the original
    # scipy.special.lpmn (used in the reference code) did NOT.
    # Multiply by (-1)^m to remove the CS phase and match lpmn's convention.
    for m_val in range(m_max + 1):
        cs_phase = (-1.0)**m_val
        for n_val in range(m_val, n_max + 1):
            p[m_val, n_val, :] = scipy.special.lpmv(m_val, n_val, x) * cs_phase

    # Derivatives via recurrence (vectorised over N_cells):
    #   dP_n^m/dx = (n*x*P_n^m - (n+m)*P_{n-1}^m) / (x^2 - 1)
    denom = x * x - 1.0  # (N_cells,)
    safe_denom = np.where(
        np.abs(denom) > 1e-14, denom, np.full_like(denom, np.inf))

    for m_val in range(m_max + 1):
        for n_val in range(m_val, n_max + 1):
            pm1 = p[m_val, n_val - 1, :] if n_val > m_val else np.zeros(N)
            pd[m_val, n_val, :] = (n_val * x * p[m_val, n_val, :] -
                                   (n_val + m_val) * pm1) / safe_denom

    # At exact poles (|x|=1) apply limiting form for m=0 (higher m → 0)
    pole_mask = np.abs(denom) <= 1e-14
    if np.any(pole_mask):
        s = np.where(x[pole_mask] > 0, 1.0, -1.0)
        for n_val in range(n_max + 1):
            pd[0, n_val,
               pole_mask] = s**(n_val + 1) * n_val * (n_val + 1) / 2.0

    return p, pd


# Saves a set of spherical harmonic coefficients for magnetic fields,
# and has functions to return magnetic vectors based on those coefficients
class magSphHarmonics:
    def __init__(self, nHarmics):

        self.nl = nHarmics
        # for a given l, m goes from 0 to l
        # so the total number of values is: sum of i+1 from i=1 to nl
        self.nTot = self.nl * (self.nl + 1) // 2 + self.nl

        # setup arrays of l and m values, to go with the subsequent arrays of coefficients
        self.l = np.zeros(self.nTot)  # noqa: E741
        self.m = np.zeros(self.nTot)
        index = 0
        for i in range(self.nl):
            for j in range(i + 2):
                self.l[index] = i + 1
                self.m[index] = j
                index += 1

        # initialize coefficients alpha, beta, and gamma.
        # actual values must be set elsewhere
        self.magGeomType = 'full'
        self.alpha = np.zeros(self.nTot, dtype=complex)
        self.beta = np.zeros(self.nTot, dtype=complex)
        self.gamma = np.zeros(self.nTot, dtype=complex)

        # XTerm, YTerm, and ZTerm (and clat, lon) are properly initialized by initMagGeom
        # This is done so that alpha, beta and gamma can be stored
        # even if we don't care about the stellar geometry.
        self.clat = []
        self.lon = []
        self.XTerm = []
        self.YTerm = []
        self.ZTerm = []
        # Hash of the last grid passed to initMagGeom; None means uninitialised
        self._init_geom_hash: int | None = None

    def setMagGeomType(self, magGeomType):
        # Set the magnetic geometry type flag to one of the allowed values.
        # Error and exit if this does not receive an allowed value.
        test = magGeomType.lower()
        if not (test == 'full' or test == 'poloidal' or test == 'pottor'
                or test == 'potential'):
            print((
                'ERROR: read an unrecognized magnetic geometry type ({:}).  ' +
                'Use one of: Full, Poloidal, PotTor, Potential').format(
                    self.magGeomType))
            import sys
            sys.exit()
        else:
            self.magGeomType = test
            # Generally the type should be set before coefficient values
            # are set, but just in case modify the coefficient accordingly.
            if self.magGeomType == 'poloidal':
                self.gamma[:] = 0.0
            elif self.magGeomType == 'pottor':
                self.beta = self.alpha
            elif self.magGeomType == 'potential':
                self.beta = self.alpha
                self.gamma[:] = 0.0

    def initMagGeom(self, clat: np.ndarray, lon: np.ndarray) -> None:
        """Batch-precompute SH basis matrices XTerm, YTerm, ZTerm.

          Replaces the former O(N_cells * n_harmonics^2) Python loop with
          O(n_harmonics^2) vectorised scipy calls over all N_cells at once.
          Results are cached by grid hash; re-calling with the same grid is a no-op.
          """
        geom_hash = hash((clat.tobytes(), lon.tobytes()))
        if geom_hash == self._init_geom_hash:
            return  # grid unchanged – reuse cached XTerm/YTerm/ZTerm
        self._init_geom_hash = geom_hash

        cos_clat = np.cos(clat)
        self.clat = clat
        self.lon = lon

        # 1. Batch computation of P_n^m and dP_n^m/dx for all cells at once
        #    p_all:  (nl+1, nl+1, N_cells)
        #    dp_all: (nl+1, nl+1, N_cells)
        p_all, dp_all = _compute_legendre_batch(self.nl, self.nl, cos_clat)

        # 2. Extract the upper-triangle (l,m) pairs in the order used by self.l / self.m
        #    np.tril_indices(nl+1) returns (row, col) for lower triangle
        #    After swapping and slicing off l=0,m=0 we get (l_idx, m_idx) for nTot entries
        ind_tri = np.tril_indices(self.nl + 1)
        ind_l = ind_tri[0][1:]  # l indices, shape (nTot,)
        ind_m = ind_tri[1][1:]  # m indices, shape (nTot,)

        p_legendre = p_all[ind_m, ind_l, :]  # (nTot, N_cells)
        dp_dcos = dp_all[ind_m, ind_l, :]  # (nTot, N_cells)

        # 3. dP/dθ = dP/d(cosθ) * d(cosθ)/dθ = dP/dx * (-sinθ)
        sin_clat = np.sin(clat)  # (N_cells,)
        dp_dth = dp_dcos * (-sin_clat[np.newaxis, :])  # (nTot, N_cells)

        # 4. Azimuthal exponential: exp(i*m*phi), shape (nTot, N_cells)
        exp_term = np.exp(1j * np.outer(self.m, lon))  # (nTot, N_cells)

        # 5. Normalisation scalars (shape (nTot,))
        c_term = np.sqrt(
            (2.0 * self.l + 1.0) / (4.0 * np.pi) * factorial(self.l - self.m) /
            factorial(self.l + self.m))
        c2_term = c_term / (self.l + 1.0)
        c3_term = 1j * c2_term * self.m
        inv_sin = 1.0 / sin_clat  # (N_cells,)

        # 6. Build and store SH basis matrices — pure matrix multiplications
        self.YTerm = c_term[:, np.
                            newaxis] * p_legendre * exp_term  # (nTot, N_cells)
        self.XTerm = (c3_term[:, np.newaxis] * inv_sin[np.newaxis, :] *
                      p_legendre * exp_term)  # (nTot, N_cells)
        self.ZTerm = c2_term[:,
                             np.newaxis] * dp_dth * exp_term  # (nTot, N_cells)

    def getAllMagVectors(self):
        # Returns an array of magnetic vectors in spherical coordinates,
        # for positions in colatitude and longitude given to initMagGeom.
        # Requires initMagGeom to have been called.

        # Matching Pascal Petit's definitions, for axes of colatitude and rotation phase
        # Br = np.real(np.sum(self.alpha*YTerm))
        # Bclat = np.real(np.sum(self.beta*ZTerm - self.gamma*XTerm))  #Btheta
        # Blon = np.real(np.sum(self.beta*XTerm + self.gamma*ZTerm))  #Bphi

        # Matching JF Donati's line profiles, making some unusual assumptions:
        # for left handed longitude in Donati's equations and -ve latitude vs colatitude term
        # Br = np.real(np.sum(self.alpha*YTerm))
        # Bclat -ve of Pascal and what is codded in Donati's ZDI, because those both use
        # latitude internally as opposed to colatitude, and the direction of increasing latitude
        # is opposite to my direction of increasing colatitude. Therefore the +ve Bclat
        # direction is in the opposite direction, hence the -ve sign here.
        # Bclat = -np.real(np.sum(self.beta*ZTerm - self.gamma*XTerm))
        # Blon for left handed longitude in Donati's equations
        # Blon = -np.real(np.sum(self.beta*XTerm + self.gamma*ZTerm))

        # Matching JF Donati's line profiles, using the complex conjugates of
        # alpha, beta and gamma, and thereby getting closest to the published equations.
        # This is closest to what appears in JF Donati ZDI code, but he has a +ve clat term
        # due to Donati using latitude increasing in the opposite direction to the
        # colatitude used here

        Br = np.real(np.dot(self.alpha, self.YTerm))
        Bclat = -np.real(
            np.dot(self.beta, self.ZTerm) + np.dot(self.gamma, self.XTerm))
        Blon = -np.real(
            np.dot(self.beta, self.XTerm) - np.dot(self.gamma, self.ZTerm))

        vMagAll = np.array([Br, Bclat, Blon])
        return vMagAll

    def getAllMagDerivs(self):
        # Returns an array of derivatives of B with respect to the spherical harmonic
        # coefficients at positions in colatitude and longitude given to initMagGeom.
        # Only calculate for components that are used in the fit
        # based on the magGeomType flag (other components are just set to 0).
        # Requires initMagGeom to have been called.

        dBr_dAlpha = np.conj(self.YTerm)
        if self.magGeomType == 'full':  # (free alpha, beta gamma)
            dBclat_dBeta = -np.conj(self.ZTerm)
            dBclat_dGamma = -np.conj(self.XTerm)
            dBlon_dBeta = -np.conj(self.XTerm)
            dBlon_dGamma = np.conj(self.ZTerm)
            derivMagAll = np.concatenate(
                (dBr_dAlpha, dBclat_dBeta, dBclat_dGamma, dBlon_dBeta,
                 dBlon_dGamma),
                axis=0)
            derivMagAll = derivMagAll.reshape(
                (5, self.nTot, self.YTerm.shape[1]))
        elif self.magGeomType == 'poloidal':  # (gamma = 0)
            dBclat_dBeta = -np.conj(self.ZTerm)
            dBlon_dBeta = -np.conj(self.XTerm)
            derivMagAll = np.concatenate(
                (dBr_dAlpha, dBclat_dBeta, dBlon_dBeta), axis=0)
            derivMagAll = derivMagAll.reshape(
                (3, self.nTot, self.YTerm.shape[1]))
        elif self.magGeomType == 'pottor':  # (beta = alpha)
            dBclat_dAlpha = -np.conj(self.ZTerm)
            dBclat_dGamma = -np.conj(self.XTerm)
            dBlon_dAlpha = -np.conj(self.XTerm)
            dBlon_dGamma = np.conj(self.ZTerm)
            derivMagAll = np.concatenate(
                (dBr_dAlpha, dBclat_dAlpha, dBclat_dGamma, dBlon_dAlpha,
                 dBlon_dGamma),
                axis=0)
            derivMagAll = derivMagAll.reshape(
                (5, self.nTot, self.YTerm.shape[1]))
        elif self.magGeomType == 'potential':  # (beta = alpha & gamma = 0)
            dBclat_dAlpha = -np.conj(self.ZTerm)
            dBlon_dAlpha = -np.conj(self.XTerm)
            derivMagAll = np.concatenate(
                (dBr_dAlpha, dBclat_dAlpha, dBlon_dAlpha))
            derivMagAll = derivMagAll.reshape(
                (3, self.nTot, self.YTerm.shape[1]))

        return derivMagAll

    def getAllMagVectorsCart(self):
        # Returns an array of magnetic vectors in Cartesian coordinates,
        # for the positions in colatitude and longitude given to initMagGeom
        # Requires initMagGeom to have been called.

        sinClat = np.sin(self.clat)
        cosClat = np.cos(self.clat)
        sinLon = np.sin(self.lon)
        cosLon = np.cos(self.lon)

        # First get the magnetic vector in spherical coordinates
        vecB = self.getAllMagVectors()
        Br = vecB[0, :]
        Bclat = vecB[1, :]
        Blon = vecB[2, :]

        # Then convert the magnetic vector from spherical to Cartesian coordinates
        # vecB = Br*r^ + Btheta*theta^ + Bphi*phi^   (where x^ denotes a unit vector)
        # and in cartesian coordinates [x^, y^, z^]  (note: theta=clat, phi=long)
        # r^ = sin(theta)cos(phi)*x^ + sin(theta)sin(phi)*y^ + cos(theta)*z^
        # theta^ = cos(theta)cos(phi)*x^ + cos(theta)sin(phi)*y^ - sin(theta)*z^
        # phi^ = -sin(phi)*x^ + cos(phi)*y^
        # then the x^, y^ and z^ components of B_vec are:
        Bx = Br * sinClat * cosLon + Bclat * cosClat * cosLon - Blon * sinLon
        By = Br * sinClat * sinLon + Bclat * cosClat * sinLon + Blon * cosLon
        Bz = Br * cosClat - Bclat * sinClat
        Bcart = np.concatenate(
            (Bx[np.newaxis], By[np.newaxis], Bz[np.newaxis]), axis=0)

        return Bcart

    def getAllMagDerivsCart(self):
        # Returns an array of derivatives of B with respect to spherical harmonic coefficients,
        # in Cartesian coordinates, for the positions given to initMagGeom
        # Requires initMagGeom to have been called.

        sinClat = np.sin(self.clat)
        cosClat = np.cos(self.clat)
        sinLon = np.sin(self.lon)
        cosLon = np.cos(self.lon)

        dB = self.getAllMagDerivs()

        # Unpack the array of derivatives in spherical coordinates
        dBr_dAlpha = dB[0, ...]
        if self.magGeomType == 'full':
            dBclat_dBeta = dB[1, ...]
            dBclat_dGamma = dB[2, ...]
            dBlon_dBeta = dB[3, ...]
            dBlon_dGamma = dB[4, ...]
        elif self.magGeomType == 'poloidal':
            dBclat_dBeta = dB[1, ...]
            dBlon_dBeta = dB[2, ...]
        elif self.magGeomType == 'pottor':
            dBclat_dAlpha = dB[1, ...]
            dBclat_dGamma = dB[2, ...]
            dBlon_dAlpha = dB[3, ...]
            dBlon_dGamma = dB[4, ...]
        elif self.magGeomType == 'potential':
            dBclat_dAlpha = dB[1, ...]
            dBlon_dAlpha = dB[2, ...]

        # Convert from spherical coordinates to Cartesian coordinates
        if self.magGeomType == 'full':
            dBx_dAlpha = dBr_dAlpha * (sinClat * cosLon)[np.newaxis, :]
            dBy_dAlpha = dBr_dAlpha * (sinClat * sinLon)[np.newaxis, :]
            dBz_dAlpha = dBr_dAlpha * cosClat[np.newaxis, :]
            dBx_dBeta = dBclat_dBeta * (cosClat * cosLon)[np.newaxis, :] \
                         - dBlon_dBeta * sinLon[np.newaxis, :]
            dBy_dBeta = dBclat_dBeta * (cosClat * sinLon)[np.newaxis, :] \
                         + dBlon_dBeta * cosLon[np.newaxis, :]
            dBz_dBeta = -dBclat_dBeta * sinClat[np.newaxis, :]
            dBx_dGamma = dBclat_dGamma * (cosClat * cosLon)[np.newaxis, :] \
                         - dBlon_dGamma * sinLon[np.newaxis, :]
            dBy_dGamma = dBclat_dGamma * (cosClat * sinLon)[np.newaxis, :] \
                         + dBlon_dGamma * cosLon[np.newaxis, :]
            dBz_dGamma = -dBclat_dGamma * sinClat[np.newaxis, :]
        elif self.magGeomType == 'poloidal':
            dBx_dAlpha = dBr_dAlpha * (sinClat * cosLon)[np.newaxis, :]
            dBy_dAlpha = dBr_dAlpha * (sinClat * sinLon)[np.newaxis, :]
            dBz_dAlpha = dBr_dAlpha * cosClat[np.newaxis, :]
            dBx_dBeta = dBclat_dBeta * (cosClat * cosLon)[np.newaxis, :] \
                         - dBlon_dBeta * sinLon[np.newaxis, :]
            dBy_dBeta = dBclat_dBeta * (cosClat * sinLon)[np.newaxis, :] \
                         + dBlon_dBeta * cosLon[np.newaxis, :]
            dBz_dBeta = -dBclat_dBeta * sinClat[np.newaxis, :]
        elif self.magGeomType == 'pottor':
            dBx_dAlpha = dBr_dAlpha * (sinClat * cosLon)[np.newaxis, :] \
                         + dBclat_dAlpha * (cosClat * cosLon)[np.newaxis, :] \
                         - dBlon_dAlpha * sinLon[np.newaxis, :]
            dBy_dAlpha = dBr_dAlpha * (sinClat * sinLon)[np.newaxis, :] \
                         + dBclat_dAlpha * (cosClat * sinLon)[np.newaxis, :] \
                         + dBlon_dAlpha * cosLon[np.newaxis, :]
            dBz_dAlpha = dBr_dAlpha * cosClat[np.newaxis, :] \
                         - dBclat_dAlpha * sinClat[np.newaxis, :]
            dBx_dGamma = dBclat_dGamma * (cosClat * cosLon)[np.newaxis, :] \
                         - dBlon_dGamma * sinLon[np.newaxis, :]
            dBy_dGamma = dBclat_dGamma * (cosClat * sinLon)[np.newaxis, :] \
                         + dBlon_dGamma * cosLon[np.newaxis, :]
            dBz_dGamma = -dBclat_dGamma * sinClat[np.newaxis, :]
        elif self.magGeomType == 'potential':
            dBx_dAlpha = dBr_dAlpha * (sinClat * cosLon)[np.newaxis, :] \
                         + dBclat_dAlpha * (cosClat * cosLon)[np.newaxis, :] \
                         - dBlon_dAlpha * sinLon[np.newaxis, :]
            dBy_dAlpha = dBr_dAlpha * (sinClat * sinLon)[np.newaxis, :] \
                         + dBclat_dAlpha * (cosClat * sinLon)[np.newaxis, :] \
                         + dBlon_dAlpha * cosLon[np.newaxis, :]
            dBz_dAlpha = dBr_dAlpha * cosClat[np.newaxis, :] \
                         - dBclat_dAlpha * sinClat[np.newaxis, :]

        # dBcart has dimensions of (x-y-z, alpha-beta-gamma, harmonic coefficient, surface element)
        # Joining the arrays here is rather slow, a single concatenation
        # follows by reshaping is an attempt to be as efficient as possible.
        if self.magGeomType == 'full':
            dBcart = np.concatenate(
                (dBx_dAlpha, dBx_dBeta, dBx_dGamma, dBy_dAlpha, dBy_dBeta,
                 dBy_dGamma, dBz_dAlpha, dBz_dBeta, dBz_dGamma))
            dBcart = np.reshape(dBcart, (3, 3, dB.shape[1], dB.shape[2]))
        elif self.magGeomType == 'poloidal':
            dBcart = np.concatenate((dBx_dAlpha, dBx_dBeta, dBy_dAlpha,
                                     dBy_dBeta, dBz_dAlpha, dBz_dBeta))
            dBcart = np.reshape(dBcart, (3, 2, dB.shape[1], dB.shape[2]))
        elif self.magGeomType == 'pottor':
            dBcart = np.concatenate((dBx_dAlpha, dBx_dGamma, dBy_dAlpha,
                                     dBy_dGamma, dBz_dAlpha, dBz_dGamma))
            dBcart = np.reshape(dBcart, (3, 2, dB.shape[1], dB.shape[2]))
        elif self.magGeomType == 'potential':
            dBcart = np.concatenate((dBx_dAlpha, dBy_dAlpha, dBz_dAlpha))
            dBcart = np.reshape(dBcart, (3, 1, dB.shape[1], dB.shape[2]))

        return dBcart

    def saveToFile(self, fName, compatibility=False):
        """将球谐系数写入文件（JF Donati 格式）。

        Parameters
        ----------
        fName : str
            输出文件路径。
        compatibility : bool, optional
            True = 完全兼容 Donati 格式（取共轭并写 -3 标志），默认 False。
        """
        with open(fName, 'w') as fOut:
            if compatibility is False:
                fOut.write('ZDIpy: general poloidal plus toroidal field\n')
                fOut.write('%i %i %i\n' % (self.nTot, 3, -30))
                for i in range(self.nTot):
                    fOut.write('%2i %2i %14e %14e\n' %
                               (self.l[i], self.m[i], np.real(
                                   self.alpha[i]), np.imag(self.alpha[i])))
                fOut.write('\n')
                for i in range(self.nTot):
                    fOut.write('%2i %2i %14e %14e\n' %
                               (self.l[i], self.m[i], np.real(
                                   self.beta[i]), np.imag(self.beta[i])))
                fOut.write('\n')
                for i in range(self.nTot):
                    fOut.write('%2i %2i %14e %14e\n' %
                               (self.l[i], self.m[i], np.real(
                                   self.gamma[i]), np.imag(self.gamma[i])))
                fOut.write('\n')
            elif compatibility is True:
                fOut.write('General poloidal plus toroidal field\n')
                fOut.write('%i %i %i\n' % (self.nTot, 3, -3))
                for i in range(self.nTot):
                    fOut.write('%2i %2i %13e %13e\n' %
                               (self.l[i], self.m[i], np.real(
                                   self.alpha[i]), -np.imag(self.alpha[i])))
                fOut.write('\n')
                for i in range(self.nTot):
                    fOut.write('%2i %2i %13e %13e\n' %
                               (self.l[i], self.m[i], np.real(
                                   self.beta[i]), -np.imag(self.beta[i])))
                fOut.write('\n')
                for i in range(self.nTot):
                    fOut.write('%2i %2i %13e %13e\n' %
                               (self.l[i], self.m[i], np.real(
                                   self.gamma[i]), -np.imag(self.gamma[i])))
                fOut.write('\n')
            else:
                print('File writing error, unknown type')


def dLegendre_dTheta(vl, vm, cosTheta):
    # assumes you are passing cos(theta),
    # returns dP(l,m)(cos(theta))/dtheta, where P(l,m) is the associated Legendre polynomial.
    # identity taken from mathworld.wolfram.com  Find a better reference?

    Plm = scipy.special.lpmv(vm, vl, cosTheta)
    Pl2m = scipy.special.lpmv(vm, vl - 1., cosTheta)
    dPdTh = (vl * cosTheta * Plm -
             (vl + vm) * Pl2m) / np.sqrt(1. - cosTheta**2)
    # if (cosTheta >= 1.) | (cosTheta <= -1.):
    #      print('ERROR in dLegendre_dTheta: value of cosTheta >= 1 or <= 1. Setting to 0.')
    #      dPdTh = 0.

    # alternately use dPdTheta formulation from mappot_rect_inc1.2.c
    # Plmplus1 = scipy.special.lpmv(vm+1, vl, cosTheta)
    # dPdTh2 = cosTheta/np.sqrt(1. - cosTheta**2)*vm*Plm + Plmplus1

    return dPdTh


def magSphHarmoicsFromFile(fname, lmax=0, verbose=1):
    """从 JF Donati 格式系数文件读取 magSphHarmonics 对象。"""
    with open(fname, 'r') as inFile:
        inFile.readline()  # skip header comment
        line = inFile.readline()
        nValues = int(line.split()[0])
        nPotential = int(line.split()[2])

        nl = int(0.5 * (np.sqrt(9 + 8 * nValues) - 3))
        if (lmax == 0):
            lmax = nl
        mSphHar = magSphHarmonics(lmax)

        for i in range(nValues):
            line = inFile.readline()
            if (i < mSphHar.nTot):
                mSphHar.l[i] = int(line.split()[0])
                mSphHar.m[i] = int(line.split()[1])
                mSphHar.alpha[i] = complex(float(line.split()[2]),
                                           float(line.split()[3]))

        inFile.readline()
        for i in range(nValues):
            line = inFile.readline()
            if (i < mSphHar.nTot):
                mSphHar.beta[i] = complex(float(line.split()[2]),
                                          float(line.split()[3]))

        inFile.readline()
        for i in range(nValues):
            line = inFile.readline()
            if (i < mSphHar.nTot):
                mSphHar.gamma[i] = complex(float(line.split()[2]),
                                           float(line.split()[3]))

    # In JF Donati's ZDI code, alpha beta and gamma are read in as real and imaginary parts,
    # but the imaginary part is treated as -ve of the imaginary part in calculating Br Btheta and Bphi.  # noqa: E501
    # Specifically, in the terms like Br = alpha*Y = alpha*P_lm(theta)*exp(i*m*phi) that code uses:
    # real(alpha)*real(Y) + imag(alpha)*imag(Y)
    # = real(alpha)*P_lm(theta)*cos(m*phi) + imag(alpha)**P_lm(theta)*sin(m*phi)
    # where as this should actually be:
    # real(alpha)*real(Y) - imag(alpha)*imag(Y)
    #  (The real part of the product of two complex numbers)
    # = real(alpha)*P_lm(theta)*cos(m*phi) - imag(alpha)**P_lm(theta)*sin(m*phi)
    # Consequently we take the complex conjugate of alpha beta and gamma
    # if we are using coefficients output by Donati's ZDI code (in potential = -3 mode).
    if (nPotential == -3):
        mSphHar.alpha = np.conj(mSphHar.alpha)
        mSphHar.beta = np.conj(mSphHar.beta)
        mSphHar.gamma = np.conj(mSphHar.gamma)
        if (verbose == 1):
            print('Treating spherical harmonics in Donati\'s "-3" ZDI fashion')

    return mSphHar


def magSphHarmoicsFromMean(lMax, Binit):
    # Misleading name/not fully fully functional (the mean field is not actually Binit)
    # need some description of the distribution of power in harmonics
    # (e.g. proportional to 1/l, or lmax-l)

    mSphHar = magSphHarmonics(lMax)

    coeffBinit = Binit / float(mSphHar.nTot)
    for i in range(mSphHar.nTot):
        mSphHar.alpha[i] = coeffBinit * (1. + 1.j)
        mSphHar.beta[i] = coeffBinit * (1. + 1.j)
        mSphHar.gamma[i] = coeffBinit * (1. + 1.j)

    return mSphHar


def SetupMagSphHarmoics(sGrid,
                        initMagFromFile,
                        initMagGeomFile,
                        lMax,
                        magGeomType='full',
                        verbose=1):
    # initialize the magnetic geometry spherical harmonics from a file of coefficients, or to a constant value (of 0).  # noqa: E501
    # Here the magGeomType flag is only used to limit what derivatives are calculated.
    # The alpha beta and gamma values all still exist and are used to calculate
    # magnetic vectors (even if some are not free parameters in the fit).
    if (initMagFromFile == 1):
        magGeom = magSphHarmoicsFromFile(initMagGeomFile, lMax, verbose)
    else:
        magGeom = magSphHarmoicsFromMean(lMax, 0.0)

    magGeom.setMagGeomType(magGeomType)
    # Save the stellar grid into magnetic geometry, and calculate spherical harmonics
    magGeom.initMagGeom(sGrid.clat, sGrid.long)

    return magGeom


# Backward-compatible alias for the old (misspelled) class name
magSphHarmoics = magSphHarmonics
