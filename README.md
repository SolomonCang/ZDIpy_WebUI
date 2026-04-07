# ZDIpy_WebUI

A modern, web-based Zeeman Doppler Imaging (ZDI) application combining:

- **Physical models** from [ZDIpy](https://github.com/SolomonCang/ZDIpy) (Folsom et al. 2018) — stellar geometry, magnetic spherical harmonics, Voigt line profiles
- **MEM optimization engine** from [pyZeeTom](https://github.com/SolomonCang/pyZeeTom) — Skilling & Bryan (1984) Maximum Entropy Method
- **Standardized `config.json`** input format replacing legacy `.dat` files
- **FastAPI + native JS WebUI** for interactive parameter editing, model execution, and result visualization

---

## Quick Start

### 1. Install dependencies

```bash
pip install -r requirements.txt
```

### 2. Launch the WebUI

```bash
python app.py
```

Open your browser at `http://localhost:7860`.

### 2.1 One-click launcher (macOS)

```bash
chmod +x start_webui.command
```

Then double-click `start_webui.command` in Finder.

The launcher will:
- create `.venv` automatically (first run only)
- install/update dependencies from `requirements.txt`
- start the app and open your browser at `http://127.0.0.1:7860`

### 3. CLI mode (no WebUI)

```bash
python app.py --cli --config config.json
# or
python zdi_runner.py --config config.json
```

---

## Project Structure

```
ZDIpy_WebUI/
├── app.py                   # Main entry point (WebUI + CLI)
├── zdi_runner.py            # CLI runner accepting config.json
├── config_loader.py         # JSON config loader / ZDIConfig class
├── requirements.txt
│
├── config.json              # Standardized input configuration
│
├── core/                    # ZDI physics engine (from ZDIpy)
│   ├── __init__.py
│   ├── mainFuncs.py         # Main fitting loop, parameter I/O
│   ├── geometryStellar.py   # Stellar surface grid
│   ├── magneticGeom.py      # Magnetic spherical harmonics
│   ├── brightnessGeom.py    # Brightness map
│   ├── lineprofileVoigt.py  # Voigt line profile model
│   ├── memSimple3.py        # MEM algorithm (Skilling & Bryan 1984)
│   ├── memSaim3.py          # MEM entropy-fitting variant
│   └── readObs.py           # LSD profile reader
│
└── LSDprof/                 # Directory for observed LSD profiles
```

> **Legacy WebUI removed**: The Gradio-based `webui_legacy/` directory was removed in favour of
> the native JS + FastAPI frontend (`frontend/` + `api/`). The new WebUI supports offline
> operation (Plotly served locally), task cancellation, and a structured logging pipeline.

---

## Configuration (`config.json`)

The `config.json` file replaces the legacy `inzdi.dat` format with a
structured, human-readable JSON schema:

```json
{
  "star": {
    "inclination_deg": 45.0,
    "vsini_kms": 67.7,
    "period_days": 0.4232,
    "differential_rotation_rad_per_day": 0.02,
    "mass_msun": 0.66,
    "radius_rsun": 0.72
  },
  "grid": { "nRings": 30 },
  "inversion": {
    "target_form": "C",
    "target_value": 1.0,
    "num_iterations": 20,
    "test_aim": 1e-4
  },
  "magnetic": {
    "fit_magnetic": 1,
    "l_max": 15,
    "default_bent": 100.0,
    "geometry_type": "Full"
  },
  "brightness": { "fit_brightness": 0, "default_bright": 1.0 },
  "line_model": { "wavelength_nm": 650.0, "lande_g": 1.195 },
  "instrument": { "spectral_resolution": 65000 },
  "velocity_grid": { "vel_start_kms": -80.0, "vel_end_kms": 80.0 },
  "observations": {
    "jdate_ref": 2456892.015,
    "files": [
      { "filename": "LSDprof/obs1.prof", "jdate": 2456886.4, "vel_center_kms": -19.8 }
    ]
  }
}
```

See [config.json](config.json) for the full schema with inline documentation.

---

## WebUI Features

| Tab | Description |
|-----|-------------|
| ⚙️ Configuration | Edit parameters with sliders/fields or raw JSON editor |
| 📂 Observations | Upload LSD profiles and manage the observation list |
| ▶️ Run Model | Run forward model or MEM inversion |
| 📊 Results | View line profile fits, magnetic maps, brightness maps |
| ℹ️ About | Documentation and references |

---

## Input File Format (LSD Profiles)

```
# comment
nPoints  nColumns
vel_kms  specI  sigI  specV  sigV  specN  sigN
...
```

---

## Output Files

All output files are written to the `results/` subdirectory (created automatically).

| File | Description |
|------|-------------|
| `results/outMagCoeff.dat` | Magnetic spherical harmonic coefficients |
| `results/outBrightMap.dat` | Brightness map (colatitude, longitude, brightness) |
| `results/outBrightMapGDark.dat` | Gravity-darkening-weighted brightness map |
| `results/outLineModels.dat` | Synthetic line profiles for all phases |
| `results/outObserved.dat` | Observed profiles used in fit |
| `results/outFitSummary.txt` | Per-iteration fit summary |

Output paths can be customised via the `"output"` block in `config.json`.

---

## Gravity Darkening

The `gravity_darkening` coefficient in `config.json → line_model` is used in the
form $g^{\beta}$. Set it to `0.0` to disable gravity darkening entirely.
For a fully radiative star use `1.0`; for convective stars consult
Claret & Bloemen (2011, A&A 529, A75) for wavelength-dependent values derived
from a wide range of model atmospheres.

---

## References

- Folsom et al. (2018) MNRAS 474, 4956 — ZDIpy
- Donati et al. (2006) MNRAS 370, 629 — ZDI magnetic mapping
- Skilling & Bryan (1984) MNRAS 211, 111 — Maximum Entropy Method
- Cang et al. (2020) A&A 643 A39 — Oblate rapid rotators

---

## License

MIT License — see [LICENSE](LICENSE).
