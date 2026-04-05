#The modified functions needed to run the maximum entropy method fitting (mem_*),
#but fitting with the constraint that entropy is >= Saim, rather than chi^2 <= chiAim.
#
# This module now delegates to core/mem_generic.py for the algorithmic parts.
# The public names getSaimQuad() and control() are kept as backward-compatible
# wrappers so that any code that imports them directly continues to work.
#

import numpy as np
from core.mem_generic import _get_saim_quad, _control_entropy, _chop_down, _chop_up


def getSaimQuad(Saim, S0, gamma, Smu):
    """
    Find the quadratic-approximation target entropy for this iteration.
    Backward-compatible wrapper for core.mem_generic._get_saim_quad.
    """
    return _get_saim_quad(Saim, S0, gamma, Smu)


# Expose chop helpers under old names for any direct callers.
chopDown = _chop_down
chopUp = _chop_up


def control(S0, gamma, Cmu, Smu, Saimq, L02, alphaMin):
    """
    SAIM control procedure (fitting to target entropy).
    Backward-compatible wrapper for core.mem_generic._control_entropy.
    """
    return _control_entropy(S0,
                            gamma,
                            Cmu,
                            Smu,
                            Saimq,
                            L02,
                            alphaMin,
                            convTol=1e-5)
