# core/mem/saim_adapter.py - SAIM (entropy-target) control wrappers
#
# Backward-compatible wrappers for the entropy-target control procedure
# originally in memSaim3.py.  Delegates to core.mem.generic.
#
# Backward-compatible re-export: core.memSaim3
#

from core.mem.generic import _get_saim_quad, _control_entropy, _chop_down, _chop_up


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
