"""
core/mem/generic.py — Generic Maximum Entropy Method (MEM) Optimization Engine
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Based on the algorithm of Skilling & Bryan (1984, MNRAS 211, 111).

This module contains **pure mathematical optimization algorithms**, independent
of specific physical models and parameterizations.  Users customize via
callback functions:
  - Entropy and its derivatives          (compute_entropy_callback)
  - Constraint statistic and derivative  (compute_constraint_callback)
  - Boundary enforcement                 (boundary_constraint_callback, optional)

The chi-squared-target control follows Skilling & Bryan directly.
The entropy-target control follows the SAIM formulation used in
memSaim3.py (fitting to a target entropy S >= Saim rather than chi^2).

Usage example: see core/mem/zdi_adapter.py  (ZDI adapter)
"""

import numpy as np
from scipy import linalg
from typing import Callable, Optional, Tuple

# ---------------------------------------------------------------------------
# Public class
# ---------------------------------------------------------------------------


class MEMOptimizer:
    """
    Generic MEM optimizer.

    Provide callbacks to define entropy forms and data fitting;
    call iterate() once per inversion step.

    Parameters
    ----------
    compute_entropy_callback : callable
        Signature (Image, weights, n1, n2, ntot) -> (S0, gradS, gradgradS, Itot)
        where Itot is the total 'image norm' used to set the step-length limit.
    compute_constraint_callback : callable
        Signature (Data, Fmodel, sig2, Resp) -> (C0, gradC)
    boundary_constraint_callback : callable, optional
        Signature (Image, n1, n2, ntot) -> Image_corrected
    max_search_dirs : int
        Maximum search directions (Skilling & Bryan recommend 3-10).
    step_length_factor : float
        L_fac in l0^2 = L_fac * Itot  (typical range 0.1–0.5).
    convergence_tol : float
        Bisection convergence tolerance for alpha and P.
    """

    def __init__(
        self,
        compute_entropy_callback: Callable,
        compute_constraint_callback: Callable,
        boundary_constraint_callback: Optional[Callable] = None,
        max_search_dirs: int = 10,
        step_length_factor: float = 0.3,
        convergence_tol: float = 1e-5,
    ):
        self.compute_entropy = compute_entropy_callback
        self.compute_constraint = compute_constraint_callback
        self.apply_boundary = boundary_constraint_callback
        self.max_search_dirs = max_search_dirs
        self.step_length_factor = step_length_factor
        self.convergence_tol = convergence_tol

    def iterate(
        self,
        Image: np.ndarray,
        Fmodel: np.ndarray,
        Data: np.ndarray,
        sig2: np.ndarray,
        Resp: np.ndarray,
        weights: np.ndarray,
        n1: int,
        n2: int,
        ntot: int,
        fixEntropy: int = 0,
        targetAim: Optional[float] = None,
    ) -> Tuple[float, float, float, np.ndarray]:
        """
        Execute one MEM iteration step.

        Parameters
        ----------
        Image   : current parameter vector (length ntot)
        Fmodel  : current model prediction  (length ndata)
        Data    : observational data         (length ndata)
        sig2    : noise variance             (length ndata)
        Resp    : response matrix dModel/dImage  (ndata × ntot)
        weights : entropy term weights       (length ntot)
        n1      : end index of standard-entropy segment
        n2      : end index of filling-factor-entropy segment  (n2 >= n1)
        ntot    : total length of Image
        fixEntropy : 0 = fit to chi-squared target; 1 = fit to entropy target
        targetAim  : target chi-squared (fixEntropy=0) or target entropy (fixEntropy=1)

        Returns
        -------
        (entropy, chi2, test, Image_new)
        """
        # ------------------------------------------------------------------
        # Constraint statistic (chi-squared) and its gradient
        # ------------------------------------------------------------------
        C0, gradC = self.compute_constraint(Data, Fmodel, sig2, Resp)

        # ------------------------------------------------------------------
        # Entropy, gradients, and image norm (Itot)
        # ------------------------------------------------------------------
        entropy_result = self.compute_entropy(Image, weights, n1, n2, ntot)
        S0, gradS, gradgradS, Itot = entropy_result

        # ------------------------------------------------------------------
        # Skilling-Bryan test statistic (Eq. 37)
        # ------------------------------------------------------------------
        test = _get_test(gradC, gradS)

        # ------------------------------------------------------------------
        # Search directions (Skilling & Bryan Eqs. 20-24)
        # ------------------------------------------------------------------
        fsupi = -1.0 / gradgradS
        edir, nedir, gamma = _search_dir(ntot, self.max_search_dirs, Resp,
                                         sig2, fsupi, gradC, gradS, gradgradS)

        # ------------------------------------------------------------------
        # Projections of gradients onto search basis (Eqs. 24a, 24c)
        # ------------------------------------------------------------------
        Cmu, Smu = _get_cmu_smu(gradC, gradS, edir)

        # Step-length limit and alpha floor
        L02 = _get_l0_squared(self.step_length_factor, Itot)
        alphaMin = _get_alpha_min(gamma)

        # ------------------------------------------------------------------
        # Control procedure: find optimal step coefficients xq
        # ------------------------------------------------------------------
        aim = targetAim if targetAim is not None else (
            S0 if fixEntropy else C0)
        if fixEntropy == 1:
            Saimq = _get_saim_quad(aim, S0, gamma, Smu)
            xq = _control_entropy(S0, gamma, Cmu, Smu, Saimq, L02, alphaMin,
                                  self.convergence_tol)
        else:
            Caimq = _get_caim_quad(aim, C0, gamma, Cmu)
            xq = _control_chi2(C0, gamma, Cmu, Smu, Caimq, L02, alphaMin,
                               self.convergence_tol)

        # ------------------------------------------------------------------
        # Update image (Eq. 25)
        # ------------------------------------------------------------------
        Image_new = Image + np.dot(edir, xq)

        # Apply boundary constraints (e.g. positivity)
        if self.apply_boundary is not None:
            Image_new = self.apply_boundary(Image_new, n1, n2, ntot)

        return S0, C0, test, Image_new


# ===========================================================================
# Private helper functions
# (use leading underscore; exported so memSimple3.py can re-use them)
# ===========================================================================


def _get_c_gradc(
    Data: np.ndarray,
    Fmodel: np.ndarray,
    sig2: np.ndarray,
    Resp: np.ndarray,
) -> Tuple[float, np.ndarray]:
    """Chi-squared and its gradient (Skilling & Bryan Eqs. 4, 8)."""
    C0 = float(np.sum((Fmodel - Data)**2 / sig2))
    gradC = 2.0 * np.sum(Resp.T * (Fmodel - Data) / sig2, axis=1)
    return C0, gradC


def _search_dir(
    ntot: int,
    maxDir: int,
    Resp: np.ndarray,
    sig2: np.ndarray,
    fsupi: np.ndarray,
    gradC: np.ndarray,
    gradS: np.ndarray,
    gradgradS: np.ndarray,
) -> Tuple[np.ndarray, int, np.ndarray]:
    """
    Generate up to maxDir linearly independent search directions
    (Skilling & Bryan Eqs. 18-20), then diagonalise the search subspace
    (Sect. 3.7.1).
    """
    edir = np.zeros((ntot, maxDir))

    # 1st direction: e_1 = f^i * grad(C)
    edir[:, 0] = gradC * fsupi
    norm_sq = np.dot(edir[:, 0], edir[:, 0])
    err_gradCis0 = 0
    if norm_sq > 0.0:
        edir[:, 0] /= np.sqrt(norm_sq)
    else:
        err_gradCis0 = 1

    if maxDir == 1:
        return _diag_dir(edir[:, :1], 1, gradgradS, sig2, Resp)

    # 2nd direction: e_2 = f^i * grad(S)
    edir[:, 1] = gradS * fsupi
    norm_sq = np.dot(edir[:, 1], edir[:, 1])
    err_gradSis0 = 0
    if norm_sq > 0.0:
        edir[:, 1] /= np.sqrt(norm_sq)
    else:
        err_gradSis0 = 1

    if (err_gradCis0 == 1) and (err_gradSis0 == 1):
        print(
            "Error: f(gradS) and f(gradC) are both zero; problem not constrained."
        )
        import sys
        sys.exit()

    if maxDir > 2:
        # Higher directions: e_n = f(grad(grad(C))) . e_{n-2}
        # grad(grad(C)) = 2 * R^T . diag(1/sigma^2) . R
        sig2_inv = 1.0 / sig2
        for i in range(2, maxDir):
            tempDot = np.dot(Resp, edir[:, i - 2])
            edir[:, i] = 2.0 * fsupi * np.dot(Resp.T, tempDot * sig2_inv)
            norm_sq = np.dot(edir[:, i], edir[:, i])
            if norm_sq > 0.0:
                edir[:, i] /= np.sqrt(norm_sq)

    return _diag_dir(edir, maxDir, gradgradS, sig2, Resp)


def _diag_dir(
    edir: np.ndarray,
    nedir: int,
    gradgradS: np.ndarray,
    sig2: np.ndarray,
    Resp: np.ndarray,
) -> Tuple[np.ndarray, int, np.ndarray]:
    """
    Diagonalise the search subspace (Skilling & Bryan Sect. 3.7.1).

    1. Build metric g_mu,nu = e^T_mu . (-gradgradS) . e_nu  (Eq. 24b).
    2. Diagonalise g; discard near-zero eigenvalue directions.
    3. Build M_mu,nu = e^T_mu . grad(grad(C)) . e_nu  (Eq. 24d).
    4. Diagonalise M; project directions onto M eigenbasis.
    """
    # Metric g (eq 24b) – only diagonal of grad(grad(S)) is used
    g = np.dot((-gradgradS * edir.T), edir)
    gamma_g, eiVec_g = linalg.eigh(g)

    # Drop near-zero eigenvalues (linear dependence)
    maxGamma = np.max(np.abs(gamma_g))
    iok = gamma_g > 1e-8 * maxGamma

    # Project and normalise (Cartesian metric)
    edir = np.dot(edir, eiVec_g[:, iok]) / np.sqrt(gamma_g[iok])
    nedir = edir.shape[1]

    # M matrix (eq 24d): M = 2 * (R.e)^T . diag(1/sig2) . (R.e)
    tmpRe = np.dot(Resp, edir)
    MM = 2.0 * np.dot((tmpRe.T / sig2), tmpRe)

    gammaM, eiVec_M = linalg.eigh(MM)
    edir = np.dot(edir, eiVec_M)

    return edir, nedir, gammaM


def _get_cmu_smu(
    gradC: np.ndarray,
    gradS: np.ndarray,
    edir: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """Projections of gradients onto search basis (Eqs. 24a, 24c)."""
    return np.dot(gradC, edir), np.dot(gradS, edir)


def _get_l0_squared(L_fac: float, Itot: float) -> float:
    """Step-length limit l0^2 = L_fac * Itot (above Eq. 28)."""
    return L_fac * Itot


def _get_alpha_min(gamma: np.ndarray) -> float:
    """
    Minimum alpha: alpha_min = max(0, -min(gamma))
    (paragraph below Eq. 32).
    """
    minGamma = float(np.min(gamma))
    return -minGamma if minGamma < 0.0 else 0.0


def _get_caim_quad(chiAim: float, C0: float, gamma: np.ndarray,
                   Cmu: np.ndarray) -> float:
    """
    Quadratic approximation of chi-squared target (Skilling & Bryan Eq. 29).

    C_aimq = (2/3) * C_minq + (1/3) * C0,  but no less than chiAim.
    """
    Cminq = C0 - 0.5 * np.sum(Cmu * Cmu / gamma)
    Caimq = 0.66666667 * Cminq + 0.33333333 * C0
    if chiAim > Caimq:
        Caimq = chiAim
    return Caimq


def _get_saim_quad(Saim: float, S0: float, gamma: np.ndarray,
                   Smu: np.ndarray) -> float:
    """
    Quadratic approximation of entropy target (SAIM formulation).

    S_maxq = S0 + 0.5 * |Smu|^2  (maximum of the entropy quadratic in the
    normalised metric basis, analogous to Eq. 28 but for entropy).
    S_aimq = (2/3) * S_maxq + (1/3) * S0,  but no more than Saim.
    """
    Smaxq = S0 + 0.5 * np.sum(Smu * Smu)
    Saimq = 0.66666667 * Smaxq + 0.33333333 * S0
    if Saim < Saimq:
        Saimq = Saim
    return Saimq


def _chop_down(alpha: float, alphaLow: float) -> Tuple[float, float]:
    """Bisect alpha down toward alphaLow."""
    alphaHigh = alpha
    alpha = 0.5 * (alphaLow + alpha)
    return alpha, alphaHigh


def _chop_up(alpha: float, alphaHigh: float) -> Tuple[float, float]:
    """Bisect alpha up toward alphaHigh (or double if unconstrained)."""
    alphaLow = alpha
    if alphaHigh > 0.0:
        alpha = 0.5 * (alphaHigh + alpha)
    else:
        alpha = 2.0 * alpha + 0.1
    return alpha, alphaLow


def _control_chi2(
    C0: float,
    gamma: np.ndarray,
    Cmu: np.ndarray,
    Smu: np.ndarray,
    Caimq: float,
    L02: float,
    alphaMin: float,
    convTol: float,
) -> np.ndarray:
    """
    Skilling-Bryan control procedure targeting chi-squared
    (Fig. 3, Sects. 3.7.2-3.7.3).

    Finds alpha and P by bisection so that Caimq <= ~C(x) < ~Cp(x) <= C0
    while |x|^2 <= l0^2, then returns the step coefficients xqp.
    """
    # Allow chi^2 to relax slightly around the target
    Caimr = Caimq * 1.001 if C0 < Caimq * 1.001 else Caimq

    P = Plow = 0.0
    Phigh = 0.0
    Pfinished = 0

    while Pfinished == 0:
        alphaLow = alphaMin
        alphaHigh = -1.0
        alpha = alphaMin + 1.0
        afinished = 0

        while afinished == 0:
            asuccess = 0

            xqp = (alpha * Smu - Cmu) / (P + gamma + alpha)
            Cq = C0 + np.dot(Cmu, xqp) + 0.5 * np.sum(gamma * xqp * xqp)
            Cqp = C0 + np.dot(Cmu, xqp) + 0.5 * np.sum((P + gamma) * xqp * xqp)
            L2 = np.dot(xqp, xqp)

            if (Cqp > C0) and (Cqp > Caimr):
                alpha, alphaHigh = _chop_down(alpha, alphaLow)
            elif (Cq < Caimq) and (Cq < C0):
                alpha, alphaLow = _chop_up(alpha, alphaHigh)
            elif L2 > L02:
                if Cqp < C0:
                    alpha, alphaLow = _chop_up(alpha, alphaHigh)
                else:
                    alpha, alphaHigh = _chop_down(alpha, alphaLow)
            else:
                asuccess = 1
                if Cq < Caimq:
                    alpha, alphaLow = _chop_up(alpha, alphaHigh)
                else:
                    alpha, alphaHigh = _chop_down(alpha, alphaLow)

            if (alphaHigh > 0.0 and
                    abs(alphaHigh - alphaLow) < convTol * alphaHigh + 1.0e-10):
                afinished = 1
            if alpha > 1.0e20:
                afinished = 1

        if asuccess != 1:
            P, Plow = _chop_up(P, Phigh)
        else:
            if P == 0.0 or abs(Phigh - Plow) < convTol * Phigh + 1.0e-10:
                Pfinished = 1
            else:
                P, Phigh = _chop_down(P, Plow)

        if (asuccess == 0) and (P > 1.0e20):
            Pfinished = 1
            print("P chop blow up encountered.")

    return xqp


def _control_entropy(
    S0: float,
    gamma: np.ndarray,
    Cmu: np.ndarray,
    Smu: np.ndarray,
    Saimq: float,
    L02: float,
    alphaMin: float,
    convTol: float,
) -> np.ndarray:
    """
    Modified control procedure targeting entropy (SAIM formulation).

    Finds alpha and P by bisection so that S0 <= ~Sp(x) < ~S(x) <= Saimq
    while |x|^2 <= l0^2, then returns the step coefficients xqp.

    The quadratic entropy approximation in the normalised (Cartesian) metric
    basis is:  ~S(x) = S0 + Smu.x - 0.5 * |x|^2
    and  ~Sp(x) = S0 + Smu.x - 0.5*(1 + P/alpha)*|x|^2

    Special case: when S0 >= 0 (initial zero-entropy start), fall back to
    alpha = alphaMin so chi-squared fitting dominates.
    """
    # Allow entropy to relax slightly around the target
    Saimr = Saimq * 1.001 if S0 > Saimq * 1.001 else Saimq

    P = Plow = 0.0
    Phigh = 0.0
    Pfinished = 0

    while Pfinished == 0:
        alphaLow = alphaMin
        alphaHigh = -1.0
        alpha = alphaMin + 1.0
        afinished = 0

        while afinished == 0:
            asuccess = 0

            xqp = (alpha * Smu - Cmu) / (P + gamma + alpha)
            L2 = np.dot(xqp, xqp)
            Sq = S0 + np.dot(Smu, xqp) - 0.5 * np.dot(xqp, xqp)
            Sqp = S0 + np.dot(
                Smu, xqp) - 0.5 * (1.0 + P / alpha) * np.dot(xqp, xqp)

            if (Sqp < S0) and (Sqp < Saimr):
                alpha, alphaLow = _chop_up(alpha, alphaHigh)
            elif (Sq > Saimq) and (Sq > S0):
                alpha, alphaHigh = _chop_down(alpha, alphaLow)
            elif L2 > L02:
                if Sqp > S0:
                    alpha, alphaHigh = _chop_down(alpha, alphaLow)
                else:
                    alpha, alphaLow = _chop_up(alpha, alphaHigh)
            else:
                asuccess = 1
                if Sq < Saimq:
                    alpha, alphaLow = _chop_up(alpha, alphaHigh)
                else:
                    alpha, alphaHigh = _chop_down(alpha, alphaLow)

            if (alphaHigh > 0.0 and
                    abs(alphaHigh - alphaLow) < convTol * alphaHigh + 1.0e-10):
                afinished = 1

            # Special case: initialised to zero entropy → defer to chi^2 fitting
            if S0 >= 0.0:
                alpha = alphaMin
                asuccess = 1
                afinished = 1

            if alpha > 1.0e20:
                afinished = 1

        if asuccess != 1:
            P, Plow = _chop_up(P, Phigh)
        else:
            if P == 0.0 or abs(Phigh - Plow) < convTol * Phigh + 1.0e-10:
                Pfinished = 1
            else:
                P, Phigh = _chop_down(P, Plow)

        if (asuccess == 0) and (P > 1.0e20):
            Pfinished = 1
            print("P chop blow up encountered.")

    return xqp


def _get_test(gradC: np.ndarray, gradS: np.ndarray) -> float:
    """
    Skilling-Bryan convergence test statistic (Eq. 37).

    Measures anti-parallelism of the entropy and chi-squared gradients.
    Range [0, 0.5]; close to 0 means near optimal.
    """
    mag_gradS = np.sqrt(np.sum(gradS**2))
    mag_gradC = np.sqrt(np.sum(gradC**2))
    inv_mag_gradS = 1.0 / mag_gradS if mag_gradS != 0.0 else 0.0
    inv_mag_gradC = 1.0 / mag_gradC if mag_gradC != 0.0 else 0.0
    return 0.5 * float(
        np.sum((gradS * inv_mag_gradS - gradC * inv_mag_gradC)**2))
