"""pipeline/result.py — ZDI inversion result dataclass."""

from __future__ import annotations

from dataclasses import dataclass, field
import numpy as np


@dataclass
class ZDIResult:
    """
    Standard result container produced by a completed ZDI inversion run.

    All fields are plain Python types or numpy arrays so the object can be
    directly JSON-serialised (via ``to_serializable()``) or saved to a
    compressed NumPy archive (via ``save()``).
    """
    iterations: int
    entropy: float
    chi2: float
    test: float
    converged: bool
    bright_map: np.ndarray  # (N_cells,)
    mag_coeffs: dict  # {"alpha": [...], "beta": [...], "gamma": [...]}
    synthetic_profiles: list[
        dict]  # per-phase {"phase": float, "vel": [...], "I_mod": [...], "V_mod": [...]}
    observed_profiles: list[
        dict]  # per-phase {"phase": float, "vel": [...], "I_obs": [...], "V_obs": [...], "I_sig": [...], "V_sig": [...]}
    light_curve_synthetic: np.ndarray | None = None  # (N_lc_obs,), None when not fitted
    metadata: dict = field(
        default_factory=dict
    )  # free key-value pairs (run time, config path, grid coords, …)

    def to_serializable(self) -> dict:
        """Return a dict that can be passed directly to ``json.dumps``.

        All numpy arrays are converted to plain Python lists.
        """
        return {
            "iterations":
            self.iterations,
            "entropy":
            self.entropy,
            "chi2":
            self.chi2,
            "test":
            self.test,
            "converged":
            self.converged,
            "bright_map":
            self.bright_map.tolist(),
            "mag_coeffs": {
                k: [[c.real, c.imag] if isinstance(c, complex) else
                    (c.tolist() if isinstance(c, np.ndarray) else c)
                    for c in (v.tolist() if isinstance(v, np.ndarray) else v)]
                for k, v in self.mag_coeffs.items()
            },
            "synthetic_profiles":
            self.synthetic_profiles,
            "observed_profiles":
            self.observed_profiles,
            "light_curve_synthetic":
            (self.light_curve_synthetic.tolist()
             if self.light_curve_synthetic is not None else None),
            "metadata": {
                k: (v.tolist() if isinstance(v, np.ndarray) else v)
                for k, v in self.metadata.items()
            },
        }

    def save(self, path: str | "Path") -> None:
        """Save numerical results to a compressed ``.npz`` file.

        Reload with ``numpy.load(path, allow_pickle=False)``.
        """
        save_kwargs: dict[str, np.ndarray] = {
            "bright_map": self.bright_map,
            "mag_alpha": np.array(self.mag_coeffs.get("alpha", []),
                                  dtype=complex),
            "mag_beta": np.array(self.mag_coeffs.get("beta", []),
                                 dtype=complex),
            "mag_gamma": np.array(self.mag_coeffs.get("gamma", []),
                                  dtype=complex),
        }
        if self.light_curve_synthetic is not None:
            save_kwargs["lc_synthetic"] = self.light_curve_synthetic
        np.savez_compressed(path, **save_kwargs)
