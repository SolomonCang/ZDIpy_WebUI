"""
core/mem_iteration_manager.py - MEM Inversion Iteration Process Manager
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Classes
-------
ConvergenceChecker
    Checks relative chi-squared change and detects stagnation.

IterationManager
    Central controller for the inversion loop: counts iterations,
    integrates ProgressMonitor and IterationHistory, and decides
    when to stop.

Factory
-------
create_iteration_manager_from_config(config, ...) -> IterationManager
"""

import time
import numpy as np
from typing import Any, Dict, Optional, Tuple

from core.mem_monitoring import IterationHistory, ProgressMonitor

# ===========================================================================
# ConvergenceChecker
# ===========================================================================


class ConvergenceChecker:
    """
    Incremental convergence checker based on relative chi-squared change.

    Calls check(chi2) each iteration; returns True once
    ``consecutive_small_changes >= stall_threshold``.

    Parameters
    ----------
    convergence_threshold : float
        Relative |Δχ²| / |χ²| threshold (default 1e-6).
    stall_threshold : int
        Number of consecutive sub-threshold changes before declaring
        convergence (default 5).
    verbose : int
        Verbosity level.
    """
    def __init__(
        self,
        convergence_threshold: float = 1e-6,
        stall_threshold: int = 5,
        verbose: int = 0,
    ):
        self.threshold = convergence_threshold
        self.stall_threshold = stall_threshold
        self.verbose = verbose

        self._prev_chi2: Optional[float] = None
        self._stall_count: int = 0

    def check(self, chi2: float) -> Tuple[bool, str]:
        """
        Return (should_stop, reason_str).

        On the first call (no previous value) always returns False.
        """
        if self._prev_chi2 is None:
            self._prev_chi2 = chi2
            return False, ""

        rel_change = abs(chi2 - self._prev_chi2) / (abs(self._prev_chi2) +
                                                    1e-20)

        if rel_change < self.threshold:
            self._stall_count += 1
        else:
            self._stall_count = 0

        self._prev_chi2 = chi2

        if self._stall_count >= self.stall_threshold:
            return True, f"convergence (Δχ²/χ²={rel_change:.3e})"
        return False, ""

    def reset(self) -> None:
        """Reset internal state."""
        self._prev_chi2 = None
        self._stall_count = 0


# ===========================================================================
# IterationManager
# ===========================================================================


class IterationManager:
    """
    Central manager for the MEM inversion loop.

    Integrates
    ----------
    * ConvergenceChecker – chi-squared stagnation
    * ProgressMonitor    – timing and ETA (optional)
    * IterationHistory   – per-iteration record   (optional)

    Typical usage::

        manager = IterationManager(max_iterations=200, verbose=1)
        for iIter in range(par.numIterations):
            manager.start_iteration()
            # … compute entropy, chi2 …
            manager.record_iteration(chi2=Chi2, entropy=S, ...)
            should_stop, reason = manager.should_stop(Chi2)
            if should_stop:
                break
        summary = manager.get_summary()

    Parameters
    ----------
    max_iterations : int
    config : dict, optional
        May contain 'convergence_threshold' and 'stall_threshold'.
    use_progress_monitor : bool
    use_iteration_history : bool
    convergence_threshold : float
    verbose : int
    """
    def __init__(
        self,
        max_iterations: int,
        config: Optional[Dict[str, Any]] = None,
        use_progress_monitor: bool = True,
        use_iteration_history: bool = True,
        convergence_threshold: float = 1e-6,
        verbose: int = 0,
    ):
        self.max_iterations = max_iterations
        self.verbose = verbose
        self.iteration: int = 0
        self.convergence_reason: str = ""
        self._start_time: float = 0.0
        self.total_elapsed: float = 0.0

        cfg = config or {}
        thr = cfg.get('convergence_threshold', convergence_threshold)
        stall = cfg.get('stall_threshold', 5)

        self.convergence_checker = ConvergenceChecker(
            convergence_threshold=thr,
            stall_threshold=stall,
            verbose=verbose,
        )

        self.progress_monitor: Optional[ProgressMonitor] = (ProgressMonitor(
            total_iterations=max_iterations,
            verbose=verbose,
            check_convergence=True,
        ) if use_progress_monitor else None)

        self.iteration_history: Optional[IterationHistory] = (IterationHistory(
            verbose=verbose) if use_iteration_history else None)

    # ------------------------------------------------------------------

    def start_iteration(self) -> None:
        """Mark the start of a new iteration (call before computation)."""
        if self.iteration == 0:
            self._start_time = time.time()
        if self.progress_monitor is not None:
            self.progress_monitor.on_iteration_start()

    def record_iteration(
        self,
        chi2: float,
        entropy: float,
        grad_S_norm: float = 0.0,
        grad_C_norm: float = 0.0,
        alpha: float = 0.0,
        param_delta: float = 0.0,
        diagnostics: Optional[Dict[str, Any]] = None,
    ) -> None:
        """
        Record results of the current iteration and advance the counter.

        Raises ValueError for any non-finite scalar.
        """
        for name, val in [('chi2', chi2), ('entropy', entropy),
                          ('grad_S_norm', grad_S_norm),
                          ('grad_C_norm', grad_C_norm)]:
            if not np.isfinite(val):
                raise ValueError(
                    f"record_iteration: {name}={val} is not finite")

        if self.progress_monitor is not None:
            self.progress_monitor.on_iteration_complete(entropy)

        if self.iteration_history is not None:
            self.iteration_history.record_iteration(
                iteration=self.iteration,
                entropy=entropy,
                chisq=chi2,
                Q=chi2 - entropy,
                grad_S_norm=grad_S_norm,
                grad_C_norm=grad_C_norm,
                alpha=alpha,
                param_delta=param_delta,
                diagnostics=diagnostics or {},
            )

        self.iteration += 1

    def should_stop(self, chi2: Optional[float] = None) -> Tuple[bool, str]:
        """
        Return (should_stop, reason).

        Checks maximum iteration count first, then convergence.
        """
        if self.iteration >= self.max_iterations:
            reason = (f"max_iterations reached "
                      f"({self.iteration}/{self.max_iterations})")
            self.convergence_reason = reason
            return True, reason

        if chi2 is not None:
            converged, reason = self.convergence_checker.check(chi2)
            if converged:
                self.convergence_reason = f"converged: {reason}"
                return True, reason

        return False, ""

    def get_summary(self) -> Dict[str, Any]:
        """Return a dictionary summarising the completed run."""
        self.total_elapsed = (time.time() - self._start_time
                              if self._start_time > 0 else 0.0)
        summary: Dict[str, Any] = {
            'iterations_completed':
            self.iteration,
            'max_iterations':
            self.max_iterations,
            'converged': (len(self.convergence_reason) > 0
                          and 'max_iterations' not in self.convergence_reason),
            'convergence_reason':
            self.convergence_reason,
            'total_elapsed_seconds':
            self.total_elapsed,
            'avg_time_per_iteration':
            (self.total_elapsed / max(1, self.iteration)),
        }
        if self.progress_monitor is not None:
            summary['progress_summary'] = self.progress_monitor.get_summary()
        if self.iteration_history is not None:
            summary['history_summary'] = self.iteration_history.get_summary()
        return summary

    def get_iteration_history(self) -> Optional[IterationHistory]:
        """Return the IterationHistory object (None if disabled)."""
        return self.iteration_history

    def reset(self) -> None:
        """Reset to initial state (keeps configuration)."""
        self.iteration = 0
        self.convergence_reason = ""
        self._start_time = 0.0
        self.total_elapsed = 0.0
        self.convergence_checker.reset()
        if self.progress_monitor is not None:
            self.progress_monitor = ProgressMonitor(
                total_iterations=self.max_iterations,
                verbose=self.verbose,
                check_convergence=True,
            )
        if self.iteration_history is not None:
            self.iteration_history = IterationHistory(verbose=self.verbose)


# ===========================================================================
# Factory helper
# ===========================================================================


def create_iteration_manager_from_config(
    config: Any,
    use_progress_monitor: bool = True,
    use_iteration_history: bool = True,
    verbose: int = 0,
) -> IterationManager:
    """
    Build an IterationManager from a config object or dict.

    Supports both dict-style (``config['max_iterations']``) and
    attribute-style (``config.numIterations``) access.
    """
    if isinstance(config, dict):
        max_iters = config.get('max_iterations') or config.get(
            'numIterations', 100)
        thr = config.get('convergence_threshold', 1e-6)
        stall = config.get('stall_threshold', 5)
    else:
        max_iters = getattr(config, 'max_iterations', None) or \
                    getattr(config, 'numIterations', 100)
        thr = getattr(config, 'convergence_threshold', 1e-6)
        stall = getattr(config, 'stall_threshold', 5)

    return IterationManager(
        max_iterations=int(max_iters),
        config={
            'convergence_threshold': thr,
            'stall_threshold': stall
        },
        use_progress_monitor=use_progress_monitor,
        use_iteration_history=use_iteration_history,
        verbose=verbose,
    )
