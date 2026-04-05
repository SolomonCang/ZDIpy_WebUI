"""
core/mem/monitoring.py - MEM Inversion Iteration Monitoring
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Classes
-------
IterationHistory
    Records the full state (entropy, chi2, gradients, timing) for every
    iteration step.  Provides a structured summary and NPZ checkpoint.

ProgressMonitor
    Tracks iteration timing, estimates ETA, and flags stagnation.
"""

import time
import numpy as np
from datetime import datetime, timedelta
from typing import Any, Dict, List, Optional

# ===========================================================================
# IterationHistory
# ===========================================================================


class IterationHistory:
    """
    Records detailed state for each MEM inversion iteration.

    Stored per-iteration
    --------------------
    iteration, entropy, chisq, Q (= chisq - entropy),
    grad_S_norm, grad_C_norm, alpha, param_delta,
    elapsed_iter, cum_elapsed, timestamp, diagnostics (dict)

    Parameters
    ----------
    verbose : int
        0 = silent, 1 = print one line per iteration.
    """
    def __init__(self, verbose: int = 0):
        self.iterations: List[Dict[str, Any]] = []
        self.start_time: float = time.time()
        self.verbose = verbose
        self.niter: int = 0

    # ------------------------------------------------------------------
    def record_iteration(
        self,
        iteration: int,
        entropy: float,
        chisq: float,
        Q: float,
        grad_S_norm: float,
        grad_C_norm: float,
        alpha: float = 0.0,
        param_delta: float = 0.0,
        diagnostics: Optional[Dict[str, Any]] = None,
        elapsed_iter: float = 0.0,
    ) -> None:
        """
        Record results of a single iteration.

        All scalar values must be finite; raises ValueError otherwise.
        """
        for name, val in [('entropy', entropy), ('chisq', chisq), ('Q', Q),
                          ('grad_S_norm', grad_S_norm),
                          ('grad_C_norm', grad_C_norm)]:
            if not np.isfinite(val):
                raise ValueError(
                    f"Invalid iteration record: {name}={val} (must be finite)")

        record: Dict[str, Any] = {
            'iteration': iteration,
            'entropy': entropy,
            'chisq': chisq,
            'Q': Q,
            'grad_S_norm': grad_S_norm,
            'grad_C_norm': grad_C_norm,
            'alpha': alpha,
            'param_delta': param_delta,
            'elapsed_iter': elapsed_iter,
            'cum_elapsed': time.time() - self.start_time,
            'timestamp': datetime.now(),
            'diagnostics': diagnostics or {},
        }
        self.iterations.append(record)
        self.niter = len(self.iterations)

        if self.verbose >= 1:
            print(f"iter {iteration:4d}: "
                  f"S={entropy:12.6e}  χ²={chisq:12.6e}  "
                  f"α={alpha:8.2e}  Δp={param_delta:8.2e}")

    # ------------------------------------------------------------------
    def get_history(self) -> Dict[str, List[float]]:
        """Return history as dict of lists (keys: entropy, chisq, Q, …)."""
        if not self.iterations:
            return {}
        keys = [
            'entropy', 'chisq', 'Q', 'grad_S_norm', 'grad_C_norm', 'alpha',
            'param_delta'
        ]
        history = {k: [it[k] for it in self.iterations] for k in keys}
        # Convenience aliases
        history['chi2'] = history['chisq']
        history['regularization'] = history['Q']
        return history

    # ------------------------------------------------------------------
    def get_last_iteration(self) -> Optional[Dict[str, Any]]:
        """Return the most recent iteration record."""
        return self.iterations[-1] if self.niter > 0 else None

    # ------------------------------------------------------------------
    def get_summary(self) -> str:
        """Return a human-readable summary string."""
        if self.niter == 0:
            return "IterationHistory: no iterations recorded."

        S_arr = np.array([it['entropy'] for it in self.iterations])
        chi_arr = np.array([it['chisq'] for it in self.iterations])
        alpha_arr = np.array([it['alpha'] for it in self.iterations])
        delta_arr = np.array([it['param_delta'] for it in self.iterations])

        lines = [
            f"IterationHistory (n_iter={self.niter})",
            f"  Entropy : {S_arr[0]:.6e} → {S_arr[-1]:.6e}  (Δ={S_arr[-1]-S_arr[0]:.6e})",
            f"  Chi2    : {chi_arr[0]:.6e} → {chi_arr[-1]:.6e}  (Δ={chi_arr[-1]-chi_arr[0]:.6e})",
            f"  alpha   : [{np.min(alpha_arr):.2e}, {np.max(alpha_arr):.2e}]",
            f"  Δparam  : [{np.min(delta_arr):.2e}, {np.max(delta_arr):.2e}]",
        ]
        if self.niter > 0:
            lines.append(
                f"  Total time: {self.iterations[-1]['cum_elapsed']:.2f}s")
        return "\n".join(lines)

    # ------------------------------------------------------------------
    def save_checkpoint(self, filepath: str) -> None:
        """Save history to an NPZ file."""
        np.savez(
            filepath,
            iterations=np.array(self.iterations, dtype=object),
            niter=np.array(self.niter),
            start_time=np.array(self.start_time),
            allow_pickle=True,
        )

    def load_checkpoint(self, filepath: str) -> None:
        """Load history from an NPZ file."""
        data = np.load(filepath, allow_pickle=True)
        self.iterations = data['iterations'].tolist()
        self.niter = int(data['niter'])
        self.start_time = float(data['start_time'])


# ===========================================================================
# ProgressMonitor
# ===========================================================================


class ProgressMonitor:
    """
    Real-time iteration progress tracker with ETA estimation.

    Parameters
    ----------
    total_iterations : int
        Expected total number of iterations.
    verbose : int
        0 = silent, 1 = per-iteration messages.
    check_convergence : bool
        If True, track stagnation (entropy change < 1e-6 relative).
    """
    def __init__(
        self,
        total_iterations: int,
        verbose: int = 0,
        check_convergence: bool = True,
    ):
        self.total_iterations = total_iterations
        self.verbose = verbose
        self.check_convergence = check_convergence

        self.iter_times: List[float] = []
        self.iter_count: int = 0
        self._t_iter_start: float = 0.0
        self.stall_count: int = 0
        self._last_entropy: Optional[float] = None

    # ------------------------------------------------------------------
    def on_iteration_start(self) -> None:
        """Call at the beginning of each iteration."""
        self._t_iter_start = time.time()

    def on_iteration_complete(self, entropy: float) -> None:
        """
        Call at the end of each iteration.

        Parameters
        ----------
        entropy : float
            Current entropy value (used for stagnation detection).
        """
        elapsed = time.time() - self._t_iter_start
        self.iter_times.append(elapsed)
        self.iter_count += 1

        if self.check_convergence and self._last_entropy is not None:
            rel_change = abs(entropy - self._last_entropy) / (
                abs(self._last_entropy) + 1e-20)
            self.stall_count = (self.stall_count +
                                1) if rel_change < 1e-6 else 0
        self._last_entropy = entropy

        if self.verbose >= 1:
            eta = self._estimate_eta()
            print(f"  iter {self.iter_count}/{self.total_iterations}  "
                  f"elapsed={elapsed:.2f}s  ETA={eta}")

    # ------------------------------------------------------------------
    def get_eta_seconds(self) -> float:
        """Estimated remaining time in seconds; -1 if unknown."""
        if len(self.iter_times) < 2:
            return -1.0
        avg = float(np.mean(self.iter_times[-5:]))
        return avg * (self.total_iterations - self.iter_count)

    def _estimate_eta(self) -> str:
        eta_sec = self.get_eta_seconds()
        if eta_sec < 0:
            return "estimating…"
        td = timedelta(seconds=eta_sec)
        return str(td).split('.')[0]

    # ------------------------------------------------------------------
    def is_stalled(self, stall_threshold: int = 3) -> bool:
        """Return True if stall_count >= stall_threshold."""
        return self.stall_count >= stall_threshold

    # ------------------------------------------------------------------
    def get_summary(self) -> str:
        """Return a human-readable progress summary."""
        if self.iter_count == 0:
            return "ProgressMonitor: no iterations completed."
        total = float(np.sum(self.iter_times))
        avg = total / self.iter_count
        pct = 100.0 * self.iter_count / self.total_iterations
        lines = [
            "ProgressMonitor",
            f"  Progress : {self.iter_count}/{self.total_iterations} ({pct:.1f}%)",
            f"  Total t  : {total:.2f}s",
            f"  Avg/iter : {avg:.3f}s",
            f"  ETA      : {self._estimate_eta()}",
            f"  Stall    : {self.stall_count} consecutive stalled iters",
        ]
        return "\n".join(lines)
