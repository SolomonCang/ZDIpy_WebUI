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

### 2.2 Python launcher (same bootstrap rules as `.command`)

```bash
python start_webui.py
```

This launcher applies the same startup rules as `start_webui.command`:

- create `.venv` automatically (first run only)
- install/update dependencies from `requirements.txt`
- open your browser automatically
- start `app.py` on port `7860` by default

You can also override the port:

```bash
python start_webui.py --port 8080
```

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
  "line_model": {
    "model_type": "voigt",
    "wavelength_nm": 650.0,
    "lande_g": 1.195,
    "...": "(voigt/unno common fields)",
    "emission_strength": 2.5,
    "emission_gauss_kms": 80.0,
    "emission_lorentz_ratio": 0.15,
    "absorption_strength": 1.2,
    "absorption_gauss_kms": 25.0,
    "absorption_lorentz_ratio": 0.10,
    "filling_factor_V": 1.0
  },
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

| Tab                | Description                                            |
| ------------------ | ------------------------------------------------------ |
| ⚙️ Configuration | Edit parameters with sliders/fields or raw JSON editor |
| 📂 Observations    | Upload LSD profiles and manage the observation list    |
| ▶️ Run Model     | Run forward model or MEM inversion                     |
| 📊 Results         | View line profile fits, magnetic maps, brightness maps |
| ℹ️ About         | Documentation and references                           |

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

| File                              | Description                                        |
| --------------------------------- | -------------------------------------------------- |
| `results/outMagCoeff.dat`       | Magnetic spherical harmonic coefficients           |
| `results/outBrightMap.dat`      | Brightness map (colatitude, longitude, brightness) |
| `results/outBrightMapGDark.dat` | Gravity-darkening-weighted brightness map          |
| `results/outLineModels.dat`     | Synthetic line profiles for all phases             |
| `results/outObserved.dat`       | Observed profiles used in fit                      |
| `results/outFitSummary.txt`     | Per-iteration fit summary                          |

Output paths can be customised via the `"output"` block in `config.json`.

---

## Spectral Line Models

The line profile model is selected via `line_model.model_type` in `config.json`:

| `model_type`      | Description                                                                                                                                                                                                                                               |
| ------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `voigt`           | Voigt weak-field approximation (default). Fast; suitable for most LSD ZDI.                                                                                                                                                                                |
| `unno`            | Unno-Rachkovsky full polarised radiative transfer (Milne-Eddington). More accurate for very strong fields or large filling factors.                                                                                                                       |
| `halpha_compound` | H-alpha double-Voigt compound model (emission + self-absorption, weak-field). For active stars showing chromospheric H-alpha emission. Response matrix is linear and can be pre-computed (same efficiency as `voigt`). See Cang et al. (2026, in prep). |

For `halpha_compound`, the extra `line_model` fields are:

| Field                        | Default | Description                                                  |
| ---------------------------- | ------- | ------------------------------------------------------------ |
| `emission_strength`        | 2.5     | Peak amplitude of the broad emission component               |
| `emission_gauss_kms`       | 80.0    | Gaussian half-width of the emission component (km/s)         |
| `emission_lorentz_ratio`   | 0.15    | Lorentzian fraction of the emission Voigt profile            |
| `absorption_strength`      | 1.2     | Depth of the narrow self-absorption component (0 = disabled) |
| `absorption_gauss_kms`     | 25.0    | Gaussian half-width of the absorption component (km/s)       |
| `absorption_lorentz_ratio` | 0.10    | Lorentzian fraction of the absorption Voigt profile          |
| `filling_factor_V`         | 1.0     | Stokes V filling factor (< 1 dilutes the Zeeman signal)      |

---

## Gravity Darkening

The `gravity_darkening` coefficient in `config.json → line_model` is used in the
form $g^{\beta}$. Set it to `0.0` to disable gravity darkening entirely.
For a fully radiative star use `1.0`; for convective stars consult
Claret & Bloemen (2011, A&A 529, A75) for wavelength-dependent values derived
from a wide range of model atmospheres.

---

## Module Documentation

Detailed physics, API reference, and algorithm notes for each submodule:

- [`docs/stellar_geometry.md`](docs/stellar_geometry.md) — `core.geometry`: 恒星表面球面等面积网格、扁球体（Roche 外形）几何、可见性与投影速度批量计算、差分自转相位修正。
- [`docs/magnetic_map.md`](docs/magnetic_map.md) — `core.magneticGeom`: 矢量磁场球谐展开（$\alpha_{\ell m}$、$\beta_{\ell m}$、$\gamma_{\ell m}$）、Legendre 多项式批量预计算、四种磁场拓扑约束（Full / Poloidal / PotTor / Potential）。
- [`docs/brightness_map.md`](docs/brightness_map.md) — `core.brightnessGeom`: 逐格点表面亮度图、临边昏暗与重力昏暗加权、MEM 熵正则化参数接口。
- [`docs/line_models.md`](docs/line_models.md) — `core.line_models`: Voigt 弱场模型与 Unno-Rachkovsky 完整偏振辐射转移两种谱线轮廓模型、盘积分与仪器卷积、响应矩阵预计算。
- [`docs/halpha_compound_model.md`](docs/halpha_compound_model.md) — `core.line_models.halpha`: H-alpha 双峰复合 Voigt 模型（宽发射 + 窄自吸收，弱场近似），专为活跃恒星色球 H-alpha 轮廓设计，线性响应矩阵可预计算。
- [`docs/MEM.md`](docs/MEM.md) — `core.mem`: Skilling & Bryan (1984) MEM 优化引擎，向量打包/解包、ZDI 专用响应矩阵组装、收敛监控。

---

## License

MIT License — see [LICENSE](LICENSE).
