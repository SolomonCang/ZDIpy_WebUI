"""
config_loader.py - Load and convert config.json to ZDIpy parameter objects.

Bridges the standardized config.json format and the legacy ZDIpy parameter
object (readParamsZDI) so that both the WebUI and the CLI can share the same
configuration.
"""

import json
import os
from pathlib import Path
import numpy as np

c_kms = 2.99792458e5  # speed of light in km/s


class ZDIConfig:
    """
    Unified ZDI configuration object loaded from config.json.

    Mirrors the interface of core.mainFuncs.readParamsZDI so that the rest of
    the ZDIpy pipeline works without modification.
    """

    def __init__(self, config_path: str):
        self.config_path = str(config_path)
        with open(config_path, "r", encoding="utf-8") as f:
            raw = json.load(f)
        self._raw = raw
        self._parse(raw)

    # ------------------------------------------------------------------
    # Parsing helpers
    # ------------------------------------------------------------------

    def _parse(self, cfg: dict):
        star = cfg["star"]
        grid = cfg["grid"]
        inv = cfg["inversion"]
        mag = cfg["magnetic"]
        bri = cfg["brightness"]
        line = cfg["line_model"]
        inst = cfg["instrument"]
        vel = cfg["velocity_grid"]
        obs_block = cfg["observations"]
        out = cfg.get("output", {})

        # --- Stellar parameters ------------------------------------------
        self.inclination = float(star["inclination_deg"])
        self.vsini = float(star["vsini_kms"])
        self.period = float(star["period_days"])
        self.dOmega = float(star["differential_rotation_rad_per_day"])
        self.mass = float(star["mass_msun"])
        self.radius = float(star["radius_rsun"])

        self.incRad = self.inclination / 180.0 * np.pi
        self.velEq = self.vsini / np.sin(self.incRad)

        # --- Grid ----------------------------------------------------------
        self.nRingsStellarGrid = int(grid["nRings"])

        # --- Inversion control ---------------------------------------------
        self.targetForm = str(inv["target_form"]).upper()
        self.targetValue = float(inv["target_value"])
        self.numIterations = int(inv["num_iterations"])
        self.test_aim = float(inv["test_aim"])

        # --- Magnetic ------------------------------------------------------
        self.fitMag = int(mag["fit_magnetic"])
        self.lMax = int(mag["l_max"])
        self.defaultBent = float(mag["default_bent"])
        self.magGeomType = str(mag["geometry_type"]).lower()
        self.initMagFromFile = int(mag["init_from_file"])
        self.initMagGeomFile = str(mag["init_file"])

        # Validate magnetic geometry type
        valid_mag_types = ("full", "poloidal", "pottor", "potential")
        if self.magGeomType not in valid_mag_types:
            raise ValueError(
                f"Invalid magnetic geometry type: {self.magGeomType!r}. "
                f"Must be one of: {valid_mag_types}"
            )

        # --- Brightness ----------------------------------------------------
        self.fitBri = int(bri["fit_brightness"])
        self.chiScaleI = float(bri["chi2_scale_I"])
        self.brightEntScale = float(bri["entropy_scale"])
        self.fEntropyBright = int(bri["entropy_form"])
        self.defaultBright = float(bri["default_bright"])
        self.maximumBright = float(bri["max_bright"])
        self.initBrightFromFile = int(bri["init_from_file"])
        self.initBrightFile = str(bri["init_file"])

        # --- Line model ----------------------------------------------------
        self.estimateStrenght = int(line["estimate_strength"])
        self.model_file = str(line.get("model_file", "model-voigt-line.dat"))

        # --- Instrument ----------------------------------------------------
        self.instrumentRes = float(inst["spectral_resolution"])

        # --- Velocity grid -------------------------------------------------
        self.velStart = float(vel["vel_start_kms"])
        self.velEnd = float(vel["vel_end_kms"])

        # --- Observations --------------------------------------------------
        self.jDateRef = float(obs_block["jdate_ref"])
        obs_files = obs_block["files"]

        self.fnames = np.array([o["filename"] for o in obs_files])
        self.jDates = np.array([float(o["jdate"]) for o in obs_files])
        self.velRs = np.array([float(o["vel_center_kms"]) for o in obs_files])
        self.numObs = len(obs_files)

        # Sanity warnings for observations
        for i in range(self.numObs):
            if abs(self.jDateRef - self.jDates[i]) > 500.0:
                print(
                    f"Warning: possible date mismatch between jDateRef "
                    f"({self.jDateRef}) and obs {i} ({self.jDates[i]})"
                )
            if abs(self.velRs[i]) > 500.0:
                print(f"Warning: extreme Vr for obs {i}: {self.velRs[i]}")

        # --- Output filenames ----------------------------------------------
        self.outMagCoeffFile = str(out.get("mag_coeff_file", "outMagCoeff.dat"))
        self.outBrightMapFile = str(out.get("bright_map_file", "outBrightMap.dat"))
        self.outBrightMapGDarkFile = str(
            out.get("bright_map_gdark_file", "outBrightMapGDark.dat")
        )
        self.outLineModelsFile = str(out.get("line_models_file", "outLineModels.dat"))
        self.outObservedFile = str(out.get("observed_used_file", "outObserved.dat"))
        self.outFitSummaryFile = str(out.get("fit_summary_file", "outFitSummary.txt"))

        # Derived / computed after parsing
        self.cycleList = None
        self.fixedEntropy = None
        self.chiTarget = None
        self.ent_aim = None
        self.calcDI = None
        self.calcDV = None

    # ------------------------------------------------------------------
    # Methods matching readParamsZDI interface
    # ------------------------------------------------------------------

    def calcCycles(self, verbose: int = 1):
        """Calculate rotation cycle/phase from Julian dates."""
        self.cycleList = (self.jDates - self.jDateRef) / self.period
        if self.dOmega != 0.0 and verbose == 1:
            lap = 2.0 * np.pi / self.dOmega
            print(
                f"Equator-pole lap time: {lap:.4f} days, "
                f"or {lap/self.period:.4f} rotation cycles"
            )
            span = np.max(self.jDates) - np.min(self.jDates)
            print(
                f"Observations span: {span:.4f} days, "
                f"or {np.max(self.cycleList)-np.min(self.cycleList):.4f} cycles"
            )

    def setTarget(self):
        """Set chi-squared or entropy target from config."""
        if self.targetForm == "C":
            self.fixedEntropy = 0
            self.chiTarget = self.targetValue
            self.ent_aim = -1e6
        elif self.targetForm == "E":
            self.fixedEntropy = 1
            self.ent_aim = self.targetValue
            self.chiTarget = 1.0
        else:
            raise ValueError(
                f"Unknown target form: {self.targetForm!r}. Must be 'C' or 'E'."
            )

    def setCalcdIdV(self, verbose: int = 1):
        """Determine which line profile components to fit."""
        self.calcDI = 0
        self.calcDV = 0
        if self.fitBri == 1:
            self.calcDI = 1
        elif self.fitBri not in (0, 1):
            raise ValueError(f"Invalid fitBri: {self.fitBri}")
        if self.fitMag == 1:
            self.calcDV = 1
        elif self.fitMag not in (0, 1):
            raise ValueError(f"Invalid fitMag: {self.fitMag}")

        if self.calcDI == 0 and self.calcDV == 0:
            if verbose == 1:
                print("Warning: no parameters to fit!")
            self.numIterations = 0

        if self.fEntropyBright not in (1, 2):
            raise ValueError(
                f"Unrecognized brightness entropy flag: {self.fEntropyBright}"
            )

    # ------------------------------------------------------------------
    # Serialization
    # ------------------------------------------------------------------

    def to_dict(self) -> dict:
        """Return the raw parsed dict (suitable for JSON serialization)."""
        return self._raw

    def save(self, path: str = None):
        """Save current configuration back to a JSON file."""
        save_path = path or self.config_path
        with open(save_path, "w", encoding="utf-8") as f:
            json.dump(self._raw, f, indent=2, ensure_ascii=False)

    @classmethod
    def from_dict(cls, cfg: dict, config_path: str = "config.json") -> "ZDIConfig":
        """Create a ZDIConfig from a plain dictionary (e.g. from the WebUI)."""
        # Write to a temp file, then load
        tmp_path = config_path
        with open(tmp_path, "w", encoding="utf-8") as f:
            json.dump(cfg, f, indent=2, ensure_ascii=False)
        return cls(tmp_path)

    @classmethod
    def default_config(cls) -> dict:
        """Return a deep copy of a minimal default configuration dictionary."""
        return {
            "star": {
                "inclination_deg": 45.0,
                "vsini_kms": 67.7,
                "period_days": 0.4232,
                "differential_rotation_rad_per_day": 0.02,
                "mass_msun": 0.66,
                "radius_rsun": 0.72,
            },
            "grid": {"nRings": 30},
            "inversion": {
                "target_form": "C",
                "target_value": 1.0,
                "num_iterations": 20,
                "test_aim": 1e-4,
            },
            "magnetic": {
                "fit_magnetic": 1,
                "l_max": 15,
                "default_bent": 100.0,
                "geometry_type": "Full",
                "init_from_file": 0,
                "init_file": "outMagCoeff.dat",
            },
            "brightness": {
                "fit_brightness": 0,
                "chi2_scale_I": 1.0,
                "entropy_scale": 1.0,
                "entropy_form": 1,
                "default_bright": 1.0,
                "max_bright": 1.01,
                "init_from_file": 0,
                "init_file": "outBrightMap.dat",
            },
            "line_model": {
                "estimate_strength": 1,
                "model_file": "model-voigt-line.dat",
                "wavelength_nm": 650.0,
                "line_strength": 0.6306,
                "gauss_width_kms": 2.41,
                "lorentz_width_fraction": 0.89,
                "lande_g": 1.195,
                "limb_darkening": 0.66,
                "gravity_darkening": 0.5,
            },
            "instrument": {"spectral_resolution": 65000.0},
            "velocity_grid": {"vel_start_kms": -80.0, "vel_end_kms": 80.0},
            "observations": {
                "jdate_ref": 2456892.015,
                "files": [],
            },
            "output": {
                "mag_coeff_file": "outMagCoeff.dat",
                "bright_map_file": "outBrightMap.dat",
                "bright_map_gdark_file": "outBrightMapGDark.dat",
                "line_models_file": "outLineModels.dat",
                "observed_used_file": "outObserved.dat",
                "fit_summary_file": "outFitSummary.txt",
            },
        }
