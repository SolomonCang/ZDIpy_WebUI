# core/mem/zdi_adapter.py - ZDI adapter layer for the MEM inversion engine
#
# Contains ZDI-specific image-vector packing/unpacking, response-matrix
# assembly, entropy weight setup, and the mem_iter / get_s_grads wrappers.
# This is the primary interface between the ZDI physical model
# (mainFuncs.py, zdi_runner.py) and the generic MEM core.
#
# Backward-compatible re-export: core.memSimple3
#
# Public API
# ----------
# mem_iter        -- one MEM iteration step (Skilling & Bryan 1984)
# get_s_grads     -- entropy value, gradients, and image norm
# updateImg       -- apply step and enforce boundary constraints
# packDataVector  -- flatten obs data into 1-D arrays (Data, sig2)
# packModelVector -- flatten synthetic spectra into 1-D model vector
# packImageVector -- flatten briMap + magGeom into 1-D image vector
# unpackImageVector -- write 1-D image vector back into briMap + magGeom
# packResponseMatrix -- assemble dModel/dImage response matrix
# constantsMEM    -- helper class: dimension bookkeeping for mem_iter
# setEntropyWeights  -- compute per-element entropy weights
#

import numpy as np

# ---------------------------------------------------------------------------
# Import algorithm primitives from the generic MEM engine.
# Backward-compatible aliases keep old names working if any caller uses them.
# ---------------------------------------------------------------------------
from core.mem.generic import (
    _get_c_gradc,
    _search_dir,
    _diag_dir,
    _get_cmu_smu,
    _get_l0_squared,
    _get_alpha_min,
    _get_caim_quad,
    _get_saim_quad,
    _chop_down,
    _chop_up,
    _control_chi2,
    _control_entropy,
    _get_test,
)

# Legacy names (kept for any code that imports them directly)
get_c_gradc = _get_c_gradc
getCmuSmu = _get_cmu_smu
getL0squared = _get_l0_squared
getAlphaMin = _get_alpha_min
getCaimQuad = _get_caim_quad
chopDown = _chop_down
chopUp = _chop_up
get_test = _get_test


def searchDir(ntot, maxDir, Resp, sig2, fsupi, gradC, gradS, gradgradS):
    """Backward-compatible wrapper for _search_dir."""
    return _search_dir(ntot, maxDir, Resp, sig2, fsupi, gradC, gradS,
                       gradgradS)


def diagDir(edir, nedir, gradgradS, sig2, Resp):
    """Backward-compatible wrapper for _diag_dir."""
    return _diag_dir(edir, nedir, gradgradS, sig2, Resp)


def control(C0, gamma, Cmu, Smu, Caimq, L02, alphaMin):
    """Backward-compatible wrapper for _control_chi2 (convTol=1e-5)."""
    return _control_chi2(C0, gamma, Cmu, Smu, Caimq, L02, alphaMin, 1e-5)


def mem_iter(n1, n2, ntot, Img, Data, Fmodel, sig2, Resp, weights, defImg,
             defIpm, ffIMax, targetAim, fixEntropy):
    """
    Maximum entropy method (MEM) image reconstruction, using the algorithm
    of Skilling & Bryan (1984, MNRAS 211, 111).
    This allows for different forms of entropy to be applied to different parts of
    the image vector, for simultainous magnetic and brightness mapping.

    :param n1: use standard image entropy (+ve, no upper limit) for the first n1 elements of Img.
    :param n2: use filling factor entropy (+ve, limiting upper value) for the n1:n2 elements of Img.  # noqa: E501
      Filling factor entropy follows A. Collier Cameron (1992, in Surface Inhomogeneities
      on Late-Type Stars, p. 33) and Y.C. Unruh & A. Collier Cameron (1995, MNRAS, 273, 1) Eq 4.
    :param ntot: total number of elements in the image vector Img.
      Positive/negative (magnetic) entropy (allows -ve & +ve, trends towards 0, no upper/lower limits)  # noqa: E501
      is used for the n2:ntot elements.  This form of entropy is from
      Hobson & Lasenby (1998, MNRAS, 298, 905) (e.g. Eq 8)
    :param Img: array of current image values (ntot long), corresponds to f in Skilling and Bryan (Eq 1).  # noqa: E501
    :param Data: array of data values the model is fit to (Eq 2).
    :param Fmodel: array of model, i.e. simulated data, calculated from the current input image Img.  # noqa: E501
      (Same dimensions as Data, and for the same points.)  For ZDI this is the model spectrum.
    :param sig2: array standard deviations (1 sigma errors) squared, for data values in Data.
      (Same dimensions as Data.)
    :param Resp: response matrix, derivatives of the model (Fmodel) with respect to the image (Img).  # noqa: E501
      That is, each element of the Resp is defined as R(k,j) = dF(k)/dI(j).
      (Dimensions of [Data.shape, Img.shape])
    :param weights: array of additional weights, which are multiplied with the entropy terms
      for each image pixel.  (Same Dimension as Img)
    :param defImg: default image intensity (Eq 6), for standard entropy.
    :param defIpm: "default" value for the positive/negative (magnetic) entropy.
    :param ffIMax: maximum image value allowed for filling factor entropy.
    :param targetAim: the desired chi-squared (Eq 4) value to converge to,
      or the desired entropy to converge to.
    :param fixEntropy: A flag for whether to fit to target chi^2 or target entropy.

    :return: entropy, chi2, test, Img, entStand, entFF, entMag
    """

    # Set the maximum number of search directions used
    maxDir = 10

    # L_fac limits the step size of the iteration in parameter space.
    L_fac = 0.3

    # for S (entropy) and C (chi^2).
    C0, gradC = _get_c_gradc(Data, Fmodel, sig2, Resp)
    S0, gradS, gradgradS, fsupi, Itot, entStand, entFF, entMag = \
        get_s_grads(n1, n2, ntot, Img, weights, defImg, defIpm, ffIMax)
    test = _get_test(gradC, gradS)

    # Calculate the normalized search directions.
    edir, nedir, gamma = _search_dir(ntot, maxDir, Resp, sig2, fsupi, gradC,
                                     gradS, gradgradS)

    # Calculate some values needed by the control subroutine.
    Cmu, Smu = _get_cmu_smu(gradC, gradS, edir)
    L02 = _get_l0_squared(L_fac, Itot)
    alphaMin = _get_alpha_min(gamma)

    if (fixEntropy == 1):
        Saim = targetAim
        Saimq = _get_saim_quad(Saim, S0, gamma, Smu)
        xq = _control_entropy(S0, gamma, Cmu, Smu, Saimq, L02, alphaMin, 1e-5)
    else:
        chiAim = targetAim
        Caimq = _get_caim_quad(chiAim, C0, gamma, Cmu)
        xq = _control_chi2(C0, gamma, Cmu, Smu, Caimq, L02, alphaMin, 1e-5)

    # Update the image vector along the search directions.
    Img = updateImg(xq, edir, Img, n1, n2, ntot, ffIMax)

    entropy = S0
    chi2 = C0

    return entropy, chi2, test, Img, entStand, entFF, entMag


def get_s_grads(n1, n2, ntot, Img, weights, defImg, defIpm, maxIff):
    # Finds the current, entropy, gradient in entropy,
    # and some additional entropy/image based quantities.

    gradS = np.zeros(ntot)
    gradgradS = np.zeros(ntot)
    fsupi = np.zeros(ntot)
    # fsupi = f^i = -1/g_ii = -1/gradgradS from eq 18 and following paragraph

    # Calculate the current image entropy S (S0) (Eq 6), gradient of S (gradS) (Eq 7),
    # and diagonal elements of the second derivative matrix of S (gradgradS)
    # (Eq 7) (off diagonal elements are zero).

    # Entropy for brightness image
    gradS[0:n1] = weights[0:n1] * (np.log(defImg) - np.log(Img[0:n1])
                                   )  # (Eq 7a)
    gradgradS[0:n1] = -weights[0:n1] / Img[0:n1]  # (Eq 7b)
    fsupi[0:n1] = -1.0 / gradgradS[0:n1]  # fsupi # f^i from eq 18
    # entStand = np.sum(-weights[0:n1]*(Img[0:n1]*(np.log(Img[0:n1]/defImg) - 1.0))) #(eq 6)
    entStand = np.sum(
        -weights[0:n1] *
        (Img[0:n1] *
         (np.log(Img[0:n1] / defImg) - 1.0) + defImg))  # normalized entropy
    Itot = np.sum(weights[0:n1] * Img[0:n1])

    # Entropy for filling factors (image with limited brightness)
    # from Collier Cameron (1992), Unruh & Collier Cameron (1995)  (Eq 4.)
    gradS[n1:n2] = weights[n1:n2] * (np.log(defImg / Img[n1:n2]) - np.log(
        (maxIff - defImg) / (maxIff - Img[n1:n2])))
    gradgradS[n1:n2] = -weights[n1:n2] * maxIff / (Img[n1:n2] *
                                                   (maxIff - Img[n1:n2]))
    fsupi[n1:n2] = -1.0 / gradgradS[n1:n2]
    entFF = np.sum(-weights[n1:n2] *
                   (Img[n1:n2] * np.log(Img[n1:n2] / defImg) +
                    (maxIff - Img[n1:n2]) * np.log(
                        (maxIff - Img[n1:n2]) / (maxIff - defImg))))
    Itot += np.sum(weights[n1:n2] * Img[n1:n2])

    # Entropy for magnetic spherical harmonic coefficients (positive and negative values)
    # from Hobson & Lasenby (1998) (Eq. 8)
    tmpPsi = np.sqrt(Img[n2:ntot]**2 + 4.0 * defIpm**2)
    gradS[n2:ntot] = -weights[n2:ntot] * np.log(
        (tmpPsi + Img[n2:ntot]) / (2.0 * defIpm))
    gradgradS[n2:ntot] = -weights[n2:ntot] / tmpPsi
    fsupi[n2:ntot] = -1.0 / gradgradS[n2:ntot]
    entMag = np.sum(weights[n2:ntot] *
                    (tmpPsi - 2.0 * defIpm - Img[n2:ntot] * np.log(
                        (tmpPsi + Img[n2:ntot]) / (2.0 * defIpm))))
    Itot += np.sum(weights[n2:ntot] * np.maximum(np.abs(Img[n2:ntot]), defIpm))

    S0 = entStand + entFF + entMag

    return S0, gradS, gradgradS, fsupi, Itot, entStand, entFF, entMag


def updateImg(xq, edir, Img, n1, n2, ntot, maxIff):
    # Update the image array (eq 25).
    # Image(new) = I + deltaI = I + x^nu*e_nu
    Img += np.inner(xq, edir)

    # Protect against stray negative image values.
    # (only for the first n1 image elements. For brightness with regular entropy)
    intMod = np.where(Img[:n1] <= 0.0)
    Img[intMod] = 1.0E-6

    # (for the next n1-n2 elements, protect against negative or too large values. For filling factor entropy)  # noqa: E501
    intMod = np.where(Img[n1:n2] <= 0.0)
    Img[intMod] = 1.0E-6 * maxIff
    intMod = np.where(Img[n1:n2] >= maxIff)
    Img[intMod] = maxIff * (1.0 - 1.0E-6)
    # elements between n2 and ntot can have any value (usually the SHD coefficients for ZDI)
    return Img


# Convert the observed data, and uncertainties, into 1D format used by mem_iter
def packDataVector(obsSet, fitBri, fitMag):

    Data = np.empty(0)
    sig2 = np.empty(0)
    if (fitBri == 1):
        for tmpObs in obsSet:
            Data = np.append(Data, tmpObs.specI)
            sig2 = np.append(sig2, tmpObs.specIsig**2)
    if (fitMag == 1):
        for tmpObs in obsSet:
            Data = np.append(Data, tmpObs.specV)
            sig2 = np.append(sig2, tmpObs.specVsig**2)

    return Data, sig2


def packResponseMatrix(setSynSpec, nDataTot, npBriMap, magGeom, magGeomType,
                       fitBri, fitMag):

    # Setup the array for I in brightness
    if (fitBri == 1):
        nDataUsed = 0
        allModeldI = np.zeros([nDataTot, npBriMap])
        for spec in setSynSpec:
            # save the derivatives of I wrt brightness
            allModeldI[nDataUsed:nDataUsed + spec.numPts, :] = spec.dIIc.T
            nDataUsed += spec.numPts

    # Setup the array for V in magnetic coefficients
    if (fitMag == 1):
        # A little sanity check:
        if magGeom.magGeomType != magGeomType:
            print('ERROR: miss-match in magGeomType flags!')
            import sys
            sys.exit()

        nDataUsed = 0
        nMagCoeff = 0
        if magGeom.magGeomType == 'full':
            nMagCoeff = 3 * 2 * magGeom.nTot
        elif magGeom.magGeomType == 'poloidal' or magGeom.magGeomType == 'pottor':
            nMagCoeff = 2 * 2 * magGeom.nTot
        elif magGeom.magGeomType == 'potential':
            nMagCoeff = 1 * 2 * magGeom.nTot
        allModeldV_lin = np.zeros((nDataTot, nMagCoeff))

        for spec in setSynSpec:
            # save the derivatives of V wrt magnetic coefficients, for the used coefficients
            nDataUsedNext = nDataUsed + spec.numPts
            if magGeom.magGeomType == 'full':
                for a in range(3):  # loop over alpha, beta, gamma
                    allModeldV_lin[nDataUsed:nDataUsedNext,
                                   a * 2 * magGeom.nTot: (a * 2 + 1) * magGeom.nTot] \
                                = np.real(spec.dVIc[a, :, :]).T
                    allModeldV_lin[nDataUsed:nDataUsedNext,
                                   (a * 2 + 1) * magGeom.nTot: (a * 2 + 2) * magGeom.nTot] \
                                = np.imag(spec.dVIc[a, :, :]).T
            elif magGeom.magGeomType == 'poloidal' or magGeom.magGeomType == 'pottor':
                for a in range(2):
                    allModeldV_lin[nDataUsed:nDataUsedNext,
                                   a * 2 * magGeom.nTot: (a * 2 + 1) * magGeom.nTot] \
                                = np.real(spec.dVIc[a, :, :]).T
                    allModeldV_lin[nDataUsed:nDataUsedNext,
                                   (a * 2 + 1) * magGeom.nTot: (a * 2 + 2) * magGeom.nTot] \
                                = np.imag(spec.dVIc[a, :, :]).T
            elif magGeom.magGeomType == 'potential':
                allModeldV_lin[nDataUsed:nDataUsedNext, 0: magGeom.nTot] \
                    = np.real(spec.dVIc[0, :, :]).T
                allModeldV_lin[nDataUsed:nDataUsedNext, magGeom.nTot: 2 * magGeom.nTot] \
                    = np.imag(spec.dVIc[0, :, :]).T

            nDataUsed += spec.numPts

    # Setup the array for V in brightness changes
    if ((fitBri == 1) & (fitMag == 1)):
        nDataUsed = 0
        allModeldVdBri = np.zeros([nDataTot, npBriMap])
        for spec in setSynSpec:
            # save the derivatives of V wrt brightness
            allModeldVdBri[nDataUsed:nDataUsed +
                           spec.numPts, :] = spec.dVdBri.T
            nDataUsed += spec.numPts

    # concatenate the two sets of derivatives (or just return one set)
    if ((fitBri == 1) & (fitMag == 1)):

        allModeldIdV = np.zeros((nDataTot * 2, npBriMap + nMagCoeff))
        # dI by dBrightness
        allModeldIdV[0:nDataTot, 0:npBriMap] = allModeldI
        # dI by dMagnetic
        # allModeldIdV[0:nDataTot, npBriMap:npBriMap+nMagCoeff] = 0.
        # dV by dBrightness
        allModeldIdV[nDataTot:nDataTot * 2, 0:npBriMap] = allModeldVdBri
        # dV by dMagnetic
        allModeldIdV[nDataTot:nDataTot * 2,
                     npBriMap:npBriMap + nMagCoeff] = allModeldV_lin
    elif (fitBri == 1):
        allModeldIdV = allModeldI
    elif (fitMag == 1):
        allModeldIdV = allModeldV_lin
    else:
        print('ERROR: got no fittable parameters in packResponseMatrix()')

    return allModeldIdV


# Merge the set of observed spectra (I and/or V) into one 1D array
def packModelVector(setSynSpec, fitBri, fitMag):
    allModelI = np.empty(0)
    allModelV = np.empty(0)
    if (fitBri == 1):
        for spec in setSynSpec:
            allModelI = np.append(allModelI, spec.IIc)
    if (fitMag == 1):
        for spec in setSynSpec:
            allModelV = np.append(allModelV, spec.VIc)

    allModelIV = np.concatenate((allModelI, allModelV))

    return allModelIV


# Convert the set of model (fitting) parameters into a 1D "image" vector for passing to mem_iter
def packImageVector(briMap, magGeom, magGeomType, fitBri, fitMag):
    # A little sanity check:
    if magGeom.magGeomType != magGeomType:
        print('ERROR: miss-match in magGeomType flags!')
        import sys
        sys.exit()

    Image = np.empty(0)
    if (fitBri == 1):
        Image = briMap.bright
    # All three alpha beta and gamma exist in the magnetic geometry object,
    # but here we only use the ones that are fit, then in upackImageVector
    # the coefficients that are not fit are set based on the type of geometry assumed.
    if (fitMag == 1):
        if magGeomType == 'full':
            Image = np.concatenate([
                Image,
                np.real(magGeom.alpha),
                np.imag(magGeom.alpha),
                np.real(magGeom.beta),
                np.imag(magGeom.beta),
                np.real(magGeom.gamma),
                np.imag(magGeom.gamma)
            ])
        elif magGeomType == 'poloidal':
            Image = np.concatenate([
                Image,
                np.real(magGeom.alpha),
                np.imag(magGeom.alpha),
                np.real(magGeom.beta),
                np.imag(magGeom.beta)
            ])
        elif magGeomType == 'pottor':
            Image = np.concatenate([
                Image,
                np.real(magGeom.alpha),
                np.imag(magGeom.alpha),
                np.real(magGeom.gamma),
                np.imag(magGeom.gamma)
            ])
        elif magGeomType == 'potential':
            Image = np.concatenate(
                [Image, np.real(magGeom.alpha),
                 np.imag(magGeom.alpha)])

    return Image


# Convert a 1D "image" vector from mem_iter back into the normal format
# for storing parameters used elsewhere.
def unpackImageVector(Image, briMap, magGeom, magGeomType, fitBri, fitMag):
    # A little sanity check:
    if magGeom.magGeomType != magGeomType:
        print('ERROR: miss-match in second magGeomType flags!')
        import sys
        sys.exit()

    npBriMap = briMap.bright.shape[0]

    endBri = 0
    if (fitBri == 1):
        briMap.bright[:] = Image[0:npBriMap]
        endBri = npBriMap

    # If only some magnetic geometry coefficients are free
    # (due to an assumed type of geometry)
    # then calculate the fixed values from the free values.
    if (fitMag == 1):
        if magGeomType == 'full':
            magGeom.alpha = Image[endBri:endBri + magGeom.nTot] + \
                            Image[endBri + magGeom.nTot:endBri + magGeom.nTot * 2] * 1.j
            magGeom.beta = Image[endBri + magGeom.nTot * 2:endBri + magGeom.nTot * 3] + \
                            Image[endBri + magGeom.nTot * 3:endBri + magGeom.nTot * 4] * 1.j
            magGeom.gamma = Image[endBri + magGeom.nTot * 4:endBri + magGeom.nTot * 5] + \
                            Image[endBri + magGeom.nTot * 5:endBri + magGeom.nTot * 6] * 1.j
        elif magGeomType == 'poloidal':
            magGeom.alpha = Image[endBri:endBri + magGeom.nTot] + \
                            Image[endBri + magGeom.nTot:endBri + magGeom.nTot * 2] * 1.j
            magGeom.beta = Image[endBri + magGeom.nTot * 2:endBri + magGeom.nTot * 3] + \
                            Image[endBri + magGeom.nTot * 3:endBri + magGeom.nTot * 4] * 1.j
            magGeom.gamma[:] = 0. + 0.j
        elif magGeomType == 'pottor':
            magGeom.alpha = Image[endBri:endBri + magGeom.nTot] + \
                            Image[endBri + magGeom.nTot:endBri + magGeom.nTot * 2] * 1.j
            magGeom.beta = Image[endBri:endBri + magGeom.nTot] + \
                            Image[endBri + magGeom.nTot:endBri + magGeom.nTot * 2] * 1.j
            magGeom.gamma = Image[endBri + magGeom.nTot * 2:endBri + magGeom.nTot * 3] + \
                            Image[endBri + magGeom.nTot * 3:endBri + magGeom.nTot * 4] * 1.j
        elif magGeomType == 'potential':
            magGeom.alpha = Image[endBri:endBri + magGeom.nTot] + \
                            Image[endBri + magGeom.nTot:endBri + magGeom.nTot * 2] * 1.j
            magGeom.beta = Image[endBri:endBri + magGeom.nTot] + \
                            Image[endBri + magGeom.nTot:endBri + magGeom.nTot * 2] * 1.j
            magGeom.gamma[:] = 0. + 0.j

    return


class constantsMEM:
    # Hold some control constants for the mem_iter routine
    def __init__(self, par, briMap, magGeom, nDataTot):
        # setup control constants for input to mem_iter
        self.npBriMap = 0
        self.npMagGeom = 0
        self.nDataTotIV = 0
        if (par.fitBri == 1):
            self.npBriMap = briMap.bright.shape[0]
            self.nDataTotIV += nDataTot
        if (par.fitMag == 1):
            if magGeom.magGeomType == 'full':
                self.npMagGeom = magGeom.nTot * 6
            elif magGeom.magGeomType == 'poloidal' or magGeom.magGeomType == 'pottor':
                self.npMagGeom = magGeom.nTot * 4
            elif magGeom.magGeomType == 'potential':
                self.npMagGeom = magGeom.nTot * 2
            else:
                self.npMagGeom = magGeom.nTot * 6
            # chi2_scale_V upweights V residuals; nDataTotIV is scaled accordingly
            # so that chi_aim = chiTarget * nDataTotIV corresponds to per-dof
            # chi2 == chiTarget for *both* I and V at convergence.
            chi_scale_v = float(getattr(par, 'chiScaleV', 1.0))
            self.nDataTotIV += nDataTot * chi_scale_v
        self.n1Model = self.npBriMap
        self.n2Model = self.n1Model
        if (par.fEntropyBright == 2):
            # If using "filling factor" entropy (limits brightness to >0, <1)
            self.n1Model = 0
            self.n2Model = self.npBriMap
        self.nModelTot = self.n2Model + self.npMagGeom


# Weight applied to the entropy terms (and their derivatives) in mem_iter
# typically it is set to the l order of the harmonic coefficient
def setEntropyWeights(par, magGeom, sGrid):
    weightEntropy = 0

    # weightEntropyV = np.ones(magGeom.nTot*6)  #alternately use no weighting (everything = 1)
    if magGeom.magGeomType == 'poloidal':
        weightEntropyV = np.tile(magGeom.l, 4)
    elif magGeom.magGeomType == 'pottor':
        weightEntropyV = np.tile(magGeom.l, 4)
    elif magGeom.magGeomType == 'potential':
        weightEntropyV = np.tile(magGeom.l, 2)
    else:
        weightEntropyV = np.tile(magGeom.l, 6)

    # Allow for a difference in relative scaling of entropy between
    # magnetic field and brightness.
    weightEntropyI = sGrid.area * par.brightEntScale
    if (par.fitBri == 1):
        weightEntropy = weightEntropyI
    if (par.fitMag == 1):
        weightEntropy = weightEntropyV
    if ((par.fitMag == 1) & (par.fitBri == 1)):
        weightEntropy = np.concatenate((weightEntropyI, weightEntropyV))
    return weightEntropy
