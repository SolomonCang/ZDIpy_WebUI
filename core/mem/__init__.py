"""
core/mem — Maximum Entropy Method (MEM) package for ZDIpy
══════════════════════════════════════════════════════════

Sub-modules
-----------
generic            Pure MEM algorithm (Skilling & Bryan 1984):
                   MEMOptimizer class and private numerical helpers.
monitoring         Per-iteration recording and real-time progress:
                   IterationHistory, ProgressMonitor.
optimization       Caching and numerical-stability tools:
                   ResponseMatrixCache, StabilityMonitor, CacheStats.
iteration_manager  High-level inversion-loop controller:
                   ConvergenceChecker, IterationManager,
                   create_iteration_manager_from_config().
zdi_adapter        ZDI-specific pack/unpack functions and mem_iter wrapper
                   (public API of the former core/memSimple3.py).
saim_adapter       SAIM entropy-target control wrappers
                   (public API of the former core/memSaim3.py).

Quick-start
-----------
    # High-level optimizer (generic, model-agnostic)
    from core.mem import MEMOptimizer, IterationManager

    # ZDI pipeline helpers
    from core.mem.zdi_adapter import (
        mem_iter, packDataVector, packModelVector,
        packImageVector, unpackImageVector,
        packResponseMatrix, constantsMEM, setEntropyWeights,
    )

    # Monitoring / caching utilities
    from core.mem import IterationHistory, ResponseMatrixCache

Backward compatibility
----------------------
The original flat-module names are still importable:

    import core.memSimple3 as memSimple   # → core.mem.zdi_adapter
    import core.memSaim3   as memSaim     # → core.mem.saim_adapter
    from core.mem_generic import ...      # → core.mem.generic
"""

# ---------------------------------------------------------------------------
# Core algorithm
# ---------------------------------------------------------------------------
from core.mem.generic import MEMOptimizer

# ---------------------------------------------------------------------------
# Monitoring
# ---------------------------------------------------------------------------
from core.mem.monitoring import IterationHistory, ProgressMonitor

# ---------------------------------------------------------------------------
# Optimization utilities
# ---------------------------------------------------------------------------
from core.mem.optimization import ResponseMatrixCache, StabilityMonitor, CacheStats

# ---------------------------------------------------------------------------
# Iteration management
# ---------------------------------------------------------------------------
from core.mem.iteration_manager import (
    ConvergenceChecker,
    IterationManager,
    create_iteration_manager_from_config,
)

# ---------------------------------------------------------------------------
# ZDI adapter (high-level pipeline functions)
# ---------------------------------------------------------------------------
from core.mem.zdi_adapter import (
    mem_iter,
    get_s_grads,
    updateImg,
    packDataVector,
    packModelVector,
    packImageVector,
    unpackImageVector,
    packResponseMatrix,
    constantsMEM,
    setEntropyWeights,
)

__all__ = [
    # generic
    "MEMOptimizer",
    # monitoring
    "IterationHistory",
    "ProgressMonitor",
    # optimization
    "ResponseMatrixCache",
    "StabilityMonitor",
    "CacheStats",
    # iteration_manager
    "ConvergenceChecker",
    "IterationManager",
    "create_iteration_manager_from_config",
    # zdi_adapter
    "mem_iter",
    "get_s_grads",
    "updateImg",
    "packDataVector",
    "packModelVector",
    "packImageVector",
    "unpackImageVector",
    "packResponseMatrix",
    "constantsMEM",
    "setEntropyWeights",
]
