"""
webui/app.py - Gradio WebUI for ZDIpy_WebUI

Provides a web-based interface for:
1. Editing ZDI parameters (config.json)
2. Managing observation files
3. Running forward model and MEM inversion
4. Visualizing results (line profile fits, magnetic maps, brightness maps)
"""

import os
import sys
import json
import io
import traceback
import threading
from pathlib import Path

# Ensure the project root is in the Python path
_ROOT = str(Path(__file__).resolve().parent.parent)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import gradio as gr
import numpy as np
import matplotlib
matplotlib.use("Agg")  # non-interactive backend for server use
import matplotlib.pyplot as plt

from config_loader import ZDIConfig

# ---------------------------------------------------------------------------
# Global state
# ---------------------------------------------------------------------------
_run_lock = threading.Lock()
_DEFAULT_CONFIG_PATH = os.path.join(_ROOT, "config", "config.json")


# ---------------------------------------------------------------------------
# Helper: load config JSON as string
# ---------------------------------------------------------------------------
def _load_config_str(path: str) -> str:
    try:
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
        # Remove comment-only keys for display
        clean = {k: v for k, v in data.items() if not k.startswith("_")}
        return json.dumps(clean, indent=2, ensure_ascii=False)
    except Exception as e:
        return f"// Error loading config: {e}"


def _save_config_str(json_str: str, path: str) -> str:
    try:
        data = json.loads(json_str)
        with open(path, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2, ensure_ascii=False)
        return f"✅ Config saved to {path}"
    except json.JSONDecodeError as e:
        return f"❌ JSON parse error: {e}"
    except Exception as e:
        return f"❌ Save error: {e}"


# ---------------------------------------------------------------------------
# Parameter widget helpers
# ---------------------------------------------------------------------------
def _cfg_get(cfg: dict, *keys, default=None):
    """Navigate nested dict with dot-separated keys."""
    d = cfg
    for k in keys:
        if not isinstance(d, dict) or k not in d:
            return default
        d = d[k]
    return d


def _cfg_set(cfg: dict, value, *keys):
    """Set nested dict value, creating dicts as needed."""
    d = cfg
    for k in keys[:-1]:
        d = d.setdefault(k, {})
    d[keys[-1]] = value


# ---------------------------------------------------------------------------
# Run ZDI pipeline
# ---------------------------------------------------------------------------
def run_zdi_pipeline(
    config_json_str: str,
    forward_only: bool,
    verbose_level: int,
    progress: gr.Progress = gr.Progress()
) -> tuple:
    """
    Run the ZDI pipeline and return (log_text, figure_profiles, figure_mag, figure_bright).
    """
    if not _run_lock.acquire(blocking=False):
        return ("⚠️ Another run is already in progress. Please wait.", None, None, None)

    log_lines = []

    class _StreamCapture:
        """Capture stdout/stderr into log_lines."""
        def __init__(self, original):
            self.original = original

        def write(self, s):
            if s.strip():
                log_lines.append(s.rstrip())
            self.original.write(s)

        def flush(self):
            self.original.flush()

    old_stdout = sys.stdout
    old_stderr = sys.stderr
    sys.stdout = _StreamCapture(old_stdout)
    sys.stderr = _StreamCapture(old_stderr)

    try:
        progress(0.05, desc="Parsing config…")

        # Parse & save config
        try:
            cfg_dict = json.loads(config_json_str)
        except json.JSONDecodeError as e:
            return (f"❌ JSON parse error: {e}", None, None, None)

        tmp_config_path = os.path.join(_ROOT, "config", "_run_config.json")
        with open(tmp_config_path, "w", encoding="utf-8") as f:
            json.dump(cfg_dict, f, indent=2, ensure_ascii=False)

        progress(0.1, desc="Loading observations…")

        # Run pipeline
        from zdi_runner import run_zdi
        result = run_zdi(
            config_path=tmp_config_path,
            forward_only=forward_only,
            verbose=verbose_level,
        )

        progress(0.85, desc="Generating plots…")

        log_lines.append("")
        log_lines.append("=" * 60)
        log_lines.append("Run Complete")
        log_lines.append("=" * 60)
        for k, v in result.items():
            log_lines.append(f"  {k:25s}: {v}")

        fig_profiles = _plot_line_profiles(result)
        fig_mag = _plot_magnetic_map(result)
        fig_bright = _plot_brightness_map(result)

        progress(1.0, desc="Done")
        return ("\n".join(log_lines), fig_profiles, fig_mag, fig_bright)

    except Exception:
        tb = traceback.format_exc()
        log_lines.append("\n❌ Error during run:\n" + tb)
        return ("\n".join(log_lines), None, None, None)

    finally:
        sys.stdout = old_stdout
        sys.stderr = old_stderr
        _run_lock.release()


# ---------------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------------
def _load_dat(fname: str):
    """Load a whitespace-delimited data file, skipping comment lines."""
    if not os.path.isfile(fname):
        return None
    rows = []
    with open(fname, "r") as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                try:
                    rows.append([float(x) for x in line.split()])
                except ValueError:
                    pass
    if not rows:
        return None
    return np.array(rows)


def _plot_line_profiles(result: dict):
    """Plot observed vs model line profiles (Stokes I and V)."""
    obs_file = os.path.join(_ROOT, "outObserved.dat")
    mod_file = os.path.join(_ROOT, "outLineModels.dat")

    if not os.path.isfile(obs_file) or not os.path.isfile(mod_file):
        return None

    # Parse multi-phase files (phases separated by blank lines)
    def _parse_phases(path):
        phases = []
        current = []
        with open(path, "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    if current:
                        phases.append(np.array(current))
                        current = []
                else:
                    try:
                        current.append([float(x) for x in line.split()])
                    except ValueError:
                        pass
        if current:
            phases.append(np.array(current))
        return phases

    obs_phases = _parse_phases(obs_file)
    mod_phases = _parse_phases(mod_file)

    n_phases = min(len(obs_phases), len(mod_phases), 16)
    if n_phases == 0:
        return None

    ncols = min(4, n_phases)
    nrows = (n_phases + ncols - 1) // ncols
    fig, axes = plt.subplots(
        nrows * 2, ncols,
        figsize=(4 * ncols, 3 * nrows * 2),
        squeeze=False,
    )
    fig.suptitle("Line Profile Fits (top: Stokes I, bottom: Stokes V)", fontsize=13)

    for idx in range(n_phases):
        row_I = (idx // ncols) * 2
        row_V = row_I + 1
        col = idx % ncols

        obs = obs_phases[idx]
        mod = mod_phases[idx]

        vel_obs = obs[:, 0] if obs.shape[1] > 0 else []
        vel_mod = mod[:, 0] if mod.shape[1] > 0 else []

        ax_I = axes[row_I][col]
        ax_V = axes[row_V][col]

        if obs.shape[1] >= 3:
            ax_I.errorbar(vel_obs, obs[:, 1], yerr=obs[:, 2],
                          fmt="k.", ms=2, lw=0.5, label="Obs")
        if obs.shape[1] >= 2:
            ax_I.plot(vel_obs, obs[:, 1], "k.", ms=2)

        if mod.shape[1] >= 2:
            ax_I.plot(vel_mod, mod[:, 1], "r-", lw=1.5, label="Model")

        ax_I.set_xlabel("Velocity (km/s)", fontsize=8)
        ax_I.set_ylabel("Stokes I/Ic", fontsize=8)
        ax_I.set_title(f"Phase {idx+1}", fontsize=9)
        ax_I.tick_params(labelsize=7)
        if idx == 0:
            ax_I.legend(fontsize=7)

        if obs.shape[1] >= 5:
            ax_V.errorbar(vel_obs, obs[:, 3], yerr=obs[:, 4],
                          fmt="k.", ms=2, lw=0.5)
        if obs.shape[1] >= 4:
            ax_V.plot(vel_obs, obs[:, 3], "k.", ms=2)

        if mod.shape[1] >= 3:
            ax_V.plot(vel_mod, mod[:, 2], "b-", lw=1.5)

        ax_V.axhline(0, color="gray", lw=0.5, ls="--")
        ax_V.set_xlabel("Velocity (km/s)", fontsize=8)
        ax_V.set_ylabel("Stokes V/Ic", fontsize=8)
        ax_V.tick_params(labelsize=7)

    # Hide unused axes
    for idx in range(n_phases, nrows * ncols):
        row_I = (idx // ncols) * 2
        row_V = row_I + 1
        col = idx % ncols
        axes[row_I][col].set_visible(False)
        axes[row_V][col].set_visible(False)

    fig.tight_layout()
    return fig


def _plot_magnetic_map(result: dict):
    """Plot magnetic map from outMagCoeff.dat using a simple surface plot."""
    mag_file = os.path.join(_ROOT, "outMagCoeff.dat")
    if not os.path.isfile(mag_file):
        return None

    # We need to reconstruct the surface magnetic field from coefficients.
    # For a simple visualization, compute Br on a colatitude/longitude grid.
    try:
        import core.magneticGeom as magneticGeom
        import core.geometryStellar as geometryStellar

        # Load config to get lMax
        cfg_path = os.path.join(_ROOT, "config", "_run_config.json")
        if not os.path.isfile(cfg_path):
            cfg_path = _DEFAULT_CONFIG_PATH
        with open(cfg_path) as f:
            cfg = json.load(f)
        lMax = int(cfg.get("magnetic", {}).get("l_max", 15))
        nRings = int(cfg.get("grid", {}).get("nRings", 30))

        sGrid = geometryStellar.starGrid(nRings, verbose=0)
        magGeom = magneticGeom.SetupMagSphHarmoics(
            sGrid, 1, mag_file, lMax, "full", verbose=0
        )
        vecB = magGeom.getAllMagVectors()  # shape (3, nPoints)
        Br = vecB[0, :]

        clat_deg = np.degrees(sGrid.clat)
        lon_deg = np.degrees(sGrid.long)

        fig, axes = plt.subplots(1, 3, figsize=(15, 5))
        labels = ["Br (G)", "B_clat (G)", "B_lon (G)"]
        for ax, comp, label in zip(axes, vecB, labels):
            sc = ax.scatter(
                lon_deg, 90 - clat_deg, c=comp,
                cmap="RdBu_r", s=4, vmin=-np.max(np.abs(comp)),
                vmax=np.max(np.abs(comp))
            )
            plt.colorbar(sc, ax=ax, label=label)
            ax.set_xlabel("Longitude (°)")
            ax.set_ylabel("Latitude (°)")
            ax.set_title(label)
        fig.suptitle("Magnetic Map", fontsize=13)
        fig.tight_layout()
        return fig

    except Exception as e:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.text(0.5, 0.5, f"Magnetic map unavailable:\n{e}",
                ha="center", va="center", transform=ax.transAxes, fontsize=10)
        ax.set_axis_off()
        return fig


def _plot_brightness_map(result: dict):
    """Plot brightness map from outBrightMap.dat."""
    bri_file = os.path.join(_ROOT, "outBrightMap.dat")
    if not os.path.isfile(bri_file):
        return None

    data = _load_dat(bri_file)
    if data is None or data.shape[1] < 3:
        return None

    clat_deg = np.degrees(data[:, 0])
    lon_deg = np.degrees(data[:, 1])
    bright = data[:, 2]

    fig, ax = plt.subplots(figsize=(8, 5))
    sc = ax.scatter(lon_deg, 90 - clat_deg, c=bright, cmap="hot_r", s=5,
                    vmin=0, vmax=max(bright.max(), 1.5))
    plt.colorbar(sc, ax=ax, label="Relative Brightness")
    ax.set_xlabel("Longitude (°)")
    ax.set_ylabel("Latitude (°)")
    ax.set_title("Brightness Map")
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Build the Gradio interface
# ---------------------------------------------------------------------------
def build_ui() -> gr.Blocks:

    default_config_str = _load_config_str(_DEFAULT_CONFIG_PATH)

    with gr.Blocks(title="ZDIpy WebUI") as demo:

        gr.Markdown(
            """
# ZDIpy WebUI — Zeeman Doppler Imaging

**Based on** [ZDIpy](https://github.com/SolomonCang/ZDIpy) (Folsom et al. 2018) physical models &
[pyZeeTom](https://github.com/SolomonCang/pyZeeTom) MEM optimization engine.

Edit parameters, upload observations, run forward/inversion models and visualize results.
"""
        )

        # ----------------------------------------------------------------
        # Tab 1: Configuration Editor
        # ----------------------------------------------------------------
        with gr.Tab("⚙️ Configuration"):
            gr.Markdown("### Edit `config.json` directly or use the parameter panels below.")

            with gr.Row():
                with gr.Column(scale=1):
                    gr.Markdown("#### Quick Parameters")

                    with gr.Accordion("🌟 Stellar Parameters", open=True):
                        inc = gr.Slider(0, 90, value=45.0, step=0.5,
                                        label="Inclination (°)")
                        vsini = gr.Slider(0, 300, value=67.7, step=0.1,
                                          label="v sin i (km/s)")
                        period = gr.Number(value=0.4232, label="Period (days)")
                        domega = gr.Number(value=0.02,
                                           label="Differential rotation Ω (rad/day)")
                        mass = gr.Number(value=0.66, label="Mass (M☉)")
                        radius = gr.Number(value=0.72, label="Radius (R☉)")

                    with gr.Accordion("🔭 Grid & Inversion", open=True):
                        nrings = gr.Slider(10, 120, value=30, step=5,
                                           label="Grid rings (nRings)")
                        target_form = gr.Radio(
                            choices=["C", "E"],
                            value="C",
                            label="Fit to: C = chi² target, E = entropy target",
                        )
                        target_val = gr.Number(value=1.0, label="Target value")
                        num_iter = gr.Slider(0, 200, value=20, step=5,
                                             label="Max iterations")
                        test_aim = gr.Number(value=1e-4,
                                             label="Convergence test_aim")

                    with gr.Accordion("🧲 Magnetic Field", open=True):
                        fit_mag = gr.Checkbox(value=True, label="Fit magnetic map")
                        l_max = gr.Slider(1, 30, value=15, step=1,
                                          label="ℓ_max (spherical harmonic order)")
                        default_bent = gr.Number(value=100.0,
                                                 label="Default B entropy slope (G)")
                        mag_geom = gr.Dropdown(
                            choices=["Full", "Poloidal", "PotTor", "Potential"],
                            value="Full",
                            label="Magnetic geometry type",
                        )

                    with gr.Accordion("☀️ Brightness Map", open=False):
                        fit_bri = gr.Checkbox(value=False, label="Fit brightness map")
                        default_bright = gr.Number(value=1.0,
                                                   label="Default brightness")

                    with gr.Accordion("📡 Instrument & Line", open=False):
                        spectral_res = gr.Number(value=65000,
                                                 label="Spectral resolution R")
                        vel_start = gr.Number(value=-80.0,
                                              label="Velocity grid start (km/s)")
                        vel_end = gr.Number(value=80.0,
                                            label="Velocity grid end (km/s)")
                        estimate_str = gr.Checkbox(
                            value=True, label="Auto-estimate line strength from EW"
                        )

                    apply_btn = gr.Button("↩ Apply to JSON Editor", variant="secondary")

                with gr.Column(scale=2):
                    gr.Markdown("#### JSON Editor")
                    config_editor = gr.Code(
                        value=default_config_str,
                        language="json",
                        label="config.json",
                        lines=40,
                        elem_classes=["code-editor"],
                    )
                    with gr.Row():
                        load_btn = gr.Button("📂 Load from disk", variant="secondary")
                        save_btn = gr.Button("💾 Save to disk", variant="primary")
                    save_status = gr.Textbox(label="Save status", interactive=False)

            # Wire quick-param panel → JSON editor
            def _apply_quick_params(
                _inc, _vsini, _period, _domega, _mass, _radius,
                _nrings, _tf, _tv, _ni, _ta,
                _fit_mag, _lmax, _bent, _geom,
                _fit_bri, _dbright,
                _sres, _vs, _ve, _est,
                current_json,
            ):
                try:
                    cfg = json.loads(current_json)
                except json.JSONDecodeError:
                    cfg = ZDIConfig.default_config()

                cfg.setdefault("star", {})
                cfg["star"]["inclination_deg"] = _inc
                cfg["star"]["vsini_kms"] = _vsini
                cfg["star"]["period_days"] = _period
                cfg["star"]["differential_rotation_rad_per_day"] = _domega
                cfg["star"]["mass_msun"] = _mass
                cfg["star"]["radius_rsun"] = _radius

                cfg.setdefault("grid", {})
                cfg["grid"]["nRings"] = int(_nrings)

                cfg.setdefault("inversion", {})
                cfg["inversion"]["target_form"] = _tf
                cfg["inversion"]["target_value"] = _tv
                cfg["inversion"]["num_iterations"] = int(_ni)
                cfg["inversion"]["test_aim"] = _ta

                cfg.setdefault("magnetic", {})
                cfg["magnetic"]["fit_magnetic"] = 1 if _fit_mag else 0
                cfg["magnetic"]["l_max"] = int(_lmax)
                cfg["magnetic"]["default_bent"] = _bent
                cfg["magnetic"]["geometry_type"] = _geom

                cfg.setdefault("brightness", {})
                cfg["brightness"]["fit_brightness"] = 1 if _fit_bri else 0
                cfg["brightness"]["default_bright"] = _dbright

                cfg.setdefault("instrument", {})
                cfg["instrument"]["spectral_resolution"] = _sres

                cfg.setdefault("velocity_grid", {})
                cfg["velocity_grid"]["vel_start_kms"] = _vs
                cfg["velocity_grid"]["vel_end_kms"] = _ve

                cfg.setdefault("line_model", {})
                cfg["line_model"]["estimate_strength"] = 1 if _est else 0

                return json.dumps(cfg, indent=2, ensure_ascii=False)

            apply_btn.click(
                fn=_apply_quick_params,
                inputs=[
                    inc, vsini, period, domega, mass, radius,
                    nrings, target_form, target_val, num_iter, test_aim,
                    fit_mag, l_max, default_bent, mag_geom,
                    fit_bri, default_bright,
                    spectral_res, vel_start, vel_end, estimate_str,
                    config_editor,
                ],
                outputs=[config_editor],
            )

            load_btn.click(
                fn=lambda: _load_config_str(_DEFAULT_CONFIG_PATH),
                outputs=[config_editor],
            )
            save_btn.click(
                fn=lambda txt: _save_config_str(txt, _DEFAULT_CONFIG_PATH),
                inputs=[config_editor],
                outputs=[save_status],
            )

        # ----------------------------------------------------------------
        # Tab 2: Observations
        # ----------------------------------------------------------------
        with gr.Tab("📂 Observations"):
            gr.Markdown(
                """
### Upload LSD Profiles

Upload your LSD profile files. They will be saved to the `LSDprof/` directory.
Then add their entries to the **Observation list** in the JSON editor on the
Configuration tab.

**Expected file format** (Donati's LSD format):
```
# header line
nPoints  nColumns
vel_kms  specI  sigI  specV  sigV  specN  sigN
...
```
"""
            )
            obs_upload = gr.File(
                label="Upload LSD profile files",
                file_count="multiple",
                file_types=[".prof", ".norm", ".lsd", ".txt", ".dat"],
            )
            obs_upload_status = gr.Textbox(label="Upload status", interactive=False)

            def _upload_obs(files):
                if not files:
                    return "No files uploaded."
                lsd_dir = os.path.join(_ROOT, "LSDprof")
                os.makedirs(lsd_dir, exist_ok=True)
                saved = []
                for f in files:
                    dest = os.path.join(lsd_dir, Path(f.name).name)
                    import shutil
                    shutil.copy(f.name, dest)
                    saved.append(dest)
                return "Saved:\n" + "\n".join(saved)

            obs_upload.change(fn=_upload_obs, inputs=[obs_upload],
                              outputs=[obs_upload_status])

            gr.Markdown("### Current Observation List")
            obs_table_btn = gr.Button("🔄 Refresh observation list")
            obs_table = gr.Textbox(
                label="Observations (JSON array, one entry per line):\n"
                      '[ {"filename": "LSDprof/obs.prof", "jdate": 2456886.4, "vel_center_kms": -19.8}, ... ]',
                lines=10,
                interactive=True,
                placeholder='[{"filename": "LSDprof/obs.prof", "jdate": 2456886.4, "vel_center_kms": -19.8}]',
            )

            save_obs_btn = gr.Button("💾 Save observation list to config", variant="primary")
            obs_save_status = gr.Textbox(label="", interactive=False)

            def _refresh_obs():
                try:
                    with open(_DEFAULT_CONFIG_PATH) as f:
                        cfg = json.load(f)
                    files = cfg.get("observations", {}).get("files", [])
                    return json.dumps(files, indent=2)
                except Exception as e:
                    return f"// Error: {e}"

            def _save_obs_table(obs_json_str):
                try:
                    obs_list = json.loads(obs_json_str)
                    if not isinstance(obs_list, list):
                        return "❌ Expected a JSON array."
                    with open(_DEFAULT_CONFIG_PATH) as f:
                        cfg = json.load(f)
                    cfg.setdefault("observations", {})["files"] = obs_list
                    with open(_DEFAULT_CONFIG_PATH, "w") as f:
                        json.dump(cfg, f, indent=2)
                    return f"✅ Saved {len(obs_list)} observations to config."
                except json.JSONDecodeError as e:
                    return f"❌ JSON error: {e}"
                except Exception as e:
                    return f"❌ Error: {e}"

            obs_table_btn.click(fn=_refresh_obs, outputs=[obs_table])
            save_obs_btn.click(fn=_save_obs_table, inputs=[obs_table],
                               outputs=[obs_save_status])

        # ----------------------------------------------------------------
        # Tab 3: Run
        # ----------------------------------------------------------------
        with gr.Tab("▶️ Run Model"):
            gr.Markdown("### Run ZDI Forward Model and/or MEM Inversion")

            with gr.Row():
                with gr.Column():
                    run_config_display = gr.Code(
                        value=default_config_str,
                        language="json",
                        label="Active config.json (read-only preview)",
                        lines=25,
                        interactive=False,
                    )
                    refresh_run_cfg = gr.Button("🔄 Reload config preview")
                    refresh_run_cfg.click(
                        fn=lambda: _load_config_str(_DEFAULT_CONFIG_PATH),
                        outputs=[run_config_display],
                    )

                with gr.Column():
                    forward_only_cb = gr.Checkbox(
                        value=False, label="Forward model only (no MEM inversion)"
                    )
                    verbose_sel = gr.Radio(
                        choices=[0, 1, 2], value=1, label="Verbosity"
                    )
            run_btn = gr.Button("🚀 Run ZDI", variant="primary", size="lg")

            run_log = gr.Textbox(
                label="Run Log",
                lines=20,
                max_lines=50,
                interactive=False,
                placeholder="Run log will appear here…",
            )

            run_fig_profiles = gr.Plot(label="Line Profile Fits", visible=False)
            run_fig_mag = gr.Plot(label="Magnetic Map", visible=False)
            run_fig_bright = gr.Plot(label="Brightness Map", visible=False)

            run_btn.click(
                fn=run_zdi_pipeline,
                inputs=[run_config_display, forward_only_cb, verbose_sel],
                outputs=[run_log, run_fig_profiles, run_fig_mag, run_fig_bright],
            )

        # ----------------------------------------------------------------
        # Tab 4: Results
        # ----------------------------------------------------------------
        with gr.Tab("📊 Results"):
            gr.Markdown("### Visualization of ZDI Results")
            gr.Markdown(
                "_After a successful run, click **Refresh** to load the latest plots._"
            )

            refresh_results_btn = gr.Button("🔄 Refresh Results", variant="primary")

            with gr.Row():
                fig_profiles = gr.Plot(label="Line Profile Fits (Stokes I & V)")
            with gr.Row():
                fig_mag = gr.Plot(label="Magnetic Field Map")
            with gr.Row():
                fig_bright = gr.Plot(label="Brightness Map")

            results_status = gr.Textbox(label="Status", interactive=False)

            def _refresh_results():
                dummy_result = {}
                fp = _plot_line_profiles(dummy_result)
                fm = _plot_magnetic_map(dummy_result)
                fb = _plot_brightness_map(dummy_result)
                msgs = []
                for name, fig in [("profiles", fp), ("magnetic", fm), ("brightness", fb)]:
                    msgs.append(f"{'✅' if fig else '⚠️'} {name}")
                return fp, fm, fb, "  |  ".join(msgs)

            refresh_results_btn.click(
                fn=_refresh_results,
                outputs=[fig_profiles, fig_mag, fig_bright, results_status],
            )

        # ----------------------------------------------------------------
        # Tab 5: About
        # ----------------------------------------------------------------
        with gr.Tab("ℹ️ About"):
            gr.Markdown(
                """
## ZDIpy WebUI

**ZDIpy WebUI** is a modern web-based interface for Zeeman Doppler Imaging (ZDI)
of stellar magnetic fields and brightness distributions.

### Physics Engine
- **Physical model**: ZDIpy by [C. P. Folsom](https://github.com/SolomonCang/ZDIpy),
  based on Donati et al. (2001, 2006) and Folsom et al. (2018, MNRAS 474, 4956).
- **MEM solver**: Enhanced Skilling & Bryan (1984) algorithm from
  [pyZeeTom](https://github.com/SolomonCang/pyZeeTom).

### Key Features
- JSON-based configuration (`config.json`) replacing legacy `.dat` files
- Web interface for parameter editing, file management, and result visualization
- Forward modelling (synthetic Stokes I + V profiles)
- MEM inversion recovering magnetic field spherical harmonic coefficients
  and stellar brightness maps
- Support for oblate rapid rotators (via gravity darkening)
- Supports magnetic geometry: Full, Poloidal, PotTor, Potential

### Quick Start
1. Edit stellar parameters in the **⚙️ Configuration** tab
2. Upload your LSD profiles in the **📂 Observations** tab
3. Click **🚀 Run ZDI** in the **▶️ Run Model** tab
4. View results in the **📊 Results** tab

### References
- Folsom et al. (2018) MNRAS 474, 4956
- Donati et al. (2006) MNRAS 370, 629
- Skilling & Bryan (1984) MNRAS 211, 111
- Cang et al. (2020) A&A 643 A39
"""
            )

    return demo


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
def launch(share: bool = False, server_port: int = 7860, **kwargs):
    """Launch the Gradio WebUI server."""
    demo = build_ui()
    demo.launch(share=share, server_port=server_port, **kwargs)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="ZDIpy WebUI server")
    parser.add_argument("--share", action="store_true",
                        help="Create a public Gradio link")
    parser.add_argument("--port", type=int, default=7860,
                        help="Server port (default: 7860)")
    args = parser.parse_args()
    launch(share=args.share, server_port=args.port)
