"""
core/mem/optimization.py - MEM Optimization Toolkit
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Provides caching and numerical-stability tools for the MEM inversion loop.

Classes
-------
ResponseMatrixCache
    LRU cache for response matrices keyed on magnetic-geometry hash.
StabilityMonitor
    Lightweight numerical-stability checker for gradients and step sizes.

Usage example::

    from core.mem_optimization import ResponseMatrixCache, StabilityMonitor
    cache = ResponseMatrixCache(max_size=5)
    Resp  = cache.get_or_compute(magGeom, obsSet, compute_fn)
    stats = cache.get_stats()

    monitor = StabilityMonitor(verbose=1)
    monitor.check_gradient(gradC, gradS)
"""

import hashlib
import warnings
import numpy as np
from typing import Any, Callable, Dict, List
from dataclasses import dataclass

# ===========================================================================
# ResponseMatrixCache
# ===========================================================================


@dataclass
class CacheStats:
    """Cache statistics snapshot."""
    hits: int = 0
    misses: int = 0
    hit_rate: float = 0.0
    size: int = 0
    max_size: int = 0
    memory_usage: int = 0  # bytes


class ResponseMatrixCache:
    """
    LRU cache for ZDI response matrices.

    Cache keys are based on a content-hash of the magnetic geometry
    coefficients (alpha, beta, gamma arrays) and the Python object id of
    the observation set (cheap, avoids re-hashing large data arrays).

    Parameters
    ----------
    max_size : int
        Maximum number of cached matrices before LRU eviction.
    verbose : int
        0 = silent, 1 = miss/evict messages, 2 = full debug.
    """
    def __init__(self, max_size: int = 5, verbose: int = 0):
        self._cache: Dict[str, np.ndarray] = {}
        self.max_size = max_size
        self._order: List[str] = []  # access order, oldest first
        self.hits = 0
        self.misses = 0
        self.verbose = verbose

    # ------------------------------------------------------------------
    def get_or_compute(
        self,
        magGeom: Any,
        obsSet: Any,
        compute_fn: Callable[[Any], np.ndarray],
    ) -> np.ndarray:
        """
        Return a cached response matrix or compute and cache it.

        Parameters
        ----------
        magGeom    : magnetic geometry object (must have alpha, beta, gamma attrs)
        obsSet     : observation set object
        compute_fn : callable(magGeom) -> np.ndarray

        Returns
        -------
        np.ndarray  response matrix (ndata × nparam)
        """
        key = self._make_key(magGeom, obsSet)

        if key in self._cache:
            self.hits += 1
            self._order.remove(key)
            self._order.append(key)
            if self.verbose >= 2:
                print(f"[Cache HIT]  key={key[:16]}…  hits={self.hits}")
            return self._cache[key]

        self.misses += 1
        Resp = compute_fn(magGeom)

        if not isinstance(Resp, np.ndarray):
            raise ValueError(
                f"compute_fn must return np.ndarray, got {type(Resp)}")
        Resp = np.asarray(Resp, dtype=float)

        if len(self._cache) >= self.max_size:
            self._evict_lru()

        self._cache[key] = Resp
        self._order.append(key)

        if self.verbose >= 1:
            print(f"[Cache MISS] key={key[:16]}…  shape={Resp.shape}  "
                  f"size={len(self._cache)}/{self.max_size}")
        return Resp

    # ------------------------------------------------------------------
    def _make_key(self, magGeom: Any, obsSet: Any) -> str:
        """Hash magnetic geometry coefficients + observation set content.

        Uses content hashing for both magGeom and obsSet so that cache
        correctness is not sensitive to Python object identity (id()).
        """
        try:

            def _h(arr: Any) -> str:
                return hashlib.md5(np.asarray(
                    arr, dtype=complex).tobytes()).hexdigest()[:8]

            obs_hash = self._hash_obs_set(obsSet)
            key = f"{_h(magGeom.alpha)}_{_h(magGeom.beta)}_{_h(magGeom.gamma)}_{obs_hash}"
        except AttributeError as exc:
            raise ValueError(
                f"magGeom must have alpha, beta, gamma attributes: {exc}")
        return key

    @staticmethod
    def _hash_obs_set(obsSet: Any) -> str:
        """Content hash for an observation set (array of obsProf objects).

        Hashes the velocity grid of each observation so that new objects
        with the same data produce the same key.  Using .tobytes() is
        efficient even for large arrays.
        """
        h = hashlib.md5()
        try:
            for obs in obsSet:
                wl = getattr(obs, "wl", None)
                if wl is not None and hasattr(wl, "tobytes"):
                    h.update(wl.tobytes())
                spec_i = getattr(obs, "specI", None)
                if spec_i is not None and hasattr(spec_i, "tobytes"):
                    h.update(spec_i.tobytes())
        except TypeError:
            # obsSet not iterable — fall back to repr hash
            h.update(repr(obsSet).encode())
        return h.hexdigest()[:12]

    def _evict_lru(self) -> None:
        if self._order:
            old = self._order.pop(0)
            del self._cache[old]
            if self.verbose >= 1:
                print(f"[Cache LRU]  evicted {old[:16]}…  "
                      f"remaining={len(self._cache)}")

    def clear(self) -> None:
        """Clear all cache entries and reset counters."""
        self._cache.clear()
        self._order.clear()
        self.hits = 0
        self.misses = 0

    def get_stats(self) -> CacheStats:
        """Return a snapshot of cache statistics."""
        total = self.hits + self.misses
        return CacheStats(
            hits=self.hits,
            misses=self.misses,
            hit_rate=self.hits / total if total > 0 else 0.0,
            size=len(self._cache),
            max_size=self.max_size,
            memory_usage=sum(a.nbytes for a in self._cache.values()),
        )

    def __repr__(self) -> str:
        s = self.get_stats()
        return (f"ResponseMatrixCache(size={s.size}/{s.max_size}, "
                f"hit_rate={s.hit_rate:.1%}, "
                f"memory={s.memory_usage / 1e6:.1f}MB)")


# ===========================================================================
# StabilityMonitor
# ===========================================================================


class StabilityMonitor:
    """
    Lightweight numerical-stability checker for MEM iterations.

    Checks
    ------
    - Gradient near-zero (convergence plateau)
    - NaN / Inf in gradients
    - Gradient norm excessively large
    - Step length out of reasonable range

    Parameters
    ----------
    verbose : int
        0 = silent, 1 = warnings on detection, 2 = all checks.
    """
    def __init__(self, verbose: int = 0):
        self.verbose = verbose
        self.warnings: List[str] = []

    # ------------------------------------------------------------------
    def check_gradient(
        self,
        gradC: np.ndarray,
        gradS: np.ndarray,
        tol: float = 1e-10,
    ) -> bool:
        """
        Check health of constraint and entropy gradients.

        Returns True if all checks pass.
        """
        gradC = np.asarray(gradC, dtype=float)
        gradS = np.asarray(gradS, dtype=float)

        if np.allclose(gradC, 0, atol=tol) and np.allclose(gradS, 0, atol=tol):
            return self._warn(
                "Gradient near-zero (possible convergence plateau)")

        if np.any(np.isnan(gradC)) or np.any(np.isnan(gradS)):
            return self._warn("NaN detected in gradients")

        if np.any(np.isinf(gradC)) or np.any(np.isinf(gradS)):
            return self._warn("Inf detected in gradients")

        max_grad = max(float(np.max(np.abs(gradC))),
                       float(np.max(np.abs(gradS))))
        if max_grad > 1e10:
            return self._warn(f"Gradient norm very large: max={max_grad:.3e}")

        return True

    # ------------------------------------------------------------------
    def check_step_length(
        self,
        step_length: float,
        min_step: float = 1e-15,
        max_step: float = 1.0,
    ) -> bool:
        """Return True if step_length is within [min_step, max_step]."""
        if step_length < min_step:
            return self._warn(f"Step size too small: {step_length:.3e}")
        if step_length > max_step:
            return self._warn(f"Step size too large: {step_length:.3e}")
        return True

    # ------------------------------------------------------------------
    def _warn(self, msg: str) -> bool:
        self.warnings.append(msg)
        if self.verbose >= 1:
            warnings.warn(msg, stacklevel=3)
        return False

    def get_summary(self) -> str:
        if not self.warnings:
            return "[StabilityMonitor] All checks passed."
        lines = [f"[StabilityMonitor] {len(self.warnings)} warning(s):"]
        for i, w in enumerate(self.warnings, 1):
            lines.append(f"  {i}. {w}")
        return "\n".join(lines)

    def clear(self) -> None:
        """Reset warning list."""
        self.warnings.clear()
