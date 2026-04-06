#!/usr/bin/env python3
"""
scripts/compare_model_profiles.py
Compare synthetic Stokes V profiles from this project (outLineModels.dat)
against the reference .model files from ZDIpy-main, and overlay the
observed LSD profiles from the LSDprof directory.

Usage:
    python scripts/compare_model_profiles.py
"""
import re
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

# ── Paths ────────────────────────────────────────────────────────────────────
NEW_MODELS = Path(__file__).parent.parent / "outLineModels.dat"
OLD_MODELS_DIR = Path(
    "/Users/tianqi/Documents/Codes_collection/ZDI_and/ZDIpy-main/LSDprof")
OBS_DIR = Path(__file__).parent.parent / "LSDprof"
PLOT_OUT = Path(__file__).parent.parent / "compare_profiles.png"

# 16 observations in the same order as config.json / inzdi.dat
OBS_FILES = [
    "lopeg_16aug14_v_02.prof.norm",
    "lopeg_19aug14_v_01.prof.norm",
    "lopeg_19aug14_v_05.prof.norm",
    "lopeg_19aug14_v_10.prof.norm",
    "lopeg_20aug14_v_01.prof.norm",
    "lopeg_20aug14_v_05.prof.norm",
    "lopeg_23aug14_v_01.prof.norm",
    "lopeg_23aug14_v_07.prof.norm",
    "lopeg_25aug14_v_01.prof.norm",
    "lopeg_25aug14_v_05.prof.norm",
    "lopeg_25aug14_v_10.prof.norm",
    "lopeg_27aug14_v_01.prof.norm",
    "lopeg_27aug14_v_07.prof.norm",
    "lopeg_31aug14_v_01.prof.norm",
    "lopeg_31aug14_v_05.prof.norm",
    "lopeg_31aug14_v_10.prof.norm",
]


def parse_new_models(path: Path):
    """Parse outLineModels.dat -> list of (phase, vel[], I[], V[])"""
    blocks = []
    phase = None
    vel, Imod, Vmod = [], [], []

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#cycle"):
                if vel:
                    blocks.append(
                        (phase, np.array(vel), np.array(Imod), np.array(Vmod)))
                    vel, Imod, Vmod = [], [], []
                m = re.search(r"[-\d.]+$", line)
                phase = float(m.group()) if m else None
            else:
                parts = line.split()
                if len(parts) >= 3:
                    vel.append(float(parts[0]))
                    Imod.append(float(parts[1]))
                    Vmod.append(float(parts[2]))
    if vel:
        blocks.append((phase, np.array(vel), np.array(Imod), np.array(Vmod)))
    return blocks


def parse_old_model(path: Path):
    """Parse a single .model file -> (phase, vel[], I[], V[])"""
    vel, Imod, Vmod = [], [], []
    phase = None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                m = re.search(r"[-\d.]+$", line)
                if m:
                    phase = float(m.group())
                continue
            parts = line.split()
            if len(parts) == 2:  # header line e.g. "89 6"
                continue
            if len(parts) >= 4:
                vel.append(float(parts[0]))
                Imod.append(float(parts[1]))
                Vmod.append(float(parts[3]))
    return phase, np.array(vel), np.array(Imod), np.array(Vmod)


def parse_obs(path: Path):
    """Parse a Donati-format LSD .prof.norm file -> (vel[], I[], Isig[], V[], Vsig[])"""
    vel, I, Isig, V, Vsig = [], [], [], [], []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) == 2:  # header "224 6"
                continue
            if len(parts) >= 5:
                vel.append(float(parts[0]))
                I.append(float(parts[1]))
                Isig.append(float(parts[2]))
                V.append(float(parts[3]))
                Vsig.append(float(parts[4]))
    return (np.array(vel), np.array(I), np.array(Isig), np.array(V),
            np.array(Vsig))


# ── Load data ────────────────────────────────────────────────────────────────
print(f"Loading new models: {NEW_MODELS}")
new_blocks = parse_new_models(NEW_MODELS)
print(f"  -> {len(new_blocks)} observation phases")

old_blocks = []
for obs in OBS_FILES:
    model_path = OLD_MODELS_DIR / (obs + ".model")
    if model_path.exists():
        old_blocks.append(parse_old_model(model_path))
    else:
        print(f"  [WARNING] not found: {model_path}")
        old_blocks.append((None, None, None, None))
print(
    f"Loaded {len([b for b in old_blocks if b[1] is not None])} reference .model files"
)

obs_data = []
for obs in OBS_FILES:
    obs_path = OBS_DIR / obs
    if obs_path.exists():
        obs_data.append(parse_obs(obs_path))
    else:
        print(f"  [WARNING] observed file not found: {obs_path}")
        obs_data.append((None, ) * 5)
print(
    f"Loaded {len([o for o in obs_data if o[0] is not None])} observed LSD profiles"
)

# ── Numerical comparison ──────────────────────────────────────────────────────
print(
    "\n=== Stokes V profile comparison (phase, rms_new, rms_ref, rms_diff) ==="
)
v_rms_new, v_rms_old, v_rms_diff = [], [], []

for i, (obs_name, (nphase, nvel, nI, nV),
        (ophase, ovel, oI,
         oV)) in enumerate(zip(OBS_FILES, new_blocks, old_blocks)):
    if ovel is None:
        print(f"  [{i+1:2d}] {obs_name:40s}  -- reference file missing")
        continue
    # Align velocity grids (ranges may differ slightly)
    vel_min = max(nvel[0], ovel[0])
    vel_max = min(nvel[-1], ovel[-1])
    mask_n = (nvel >= vel_min - 0.01) & (nvel <= vel_max + 0.01)
    mask_o = (ovel >= vel_min - 0.01) & (ovel <= vel_max + 0.01)
    nV_trim = nV[mask_n]
    oV_trim = oV[mask_o]
    if len(nV_trim) != len(oV_trim):
        oV_trim = np.interp(nvel[mask_n], ovel[mask_o], oV_trim)
    rms_n = float(np.sqrt(np.mean(nV_trim**2)))
    rms_o = float(np.sqrt(np.mean(oV_trim**2)))
    rms_d = float(np.sqrt(np.mean((nV_trim - oV_trim)**2)))
    v_rms_new.append(rms_n)
    v_rms_old.append(rms_o)
    v_rms_diff.append(rms_d)
    print(f"  [{i+1:2d}] phase={nphase:9.4f}  rms_new={rms_n:.3e}  "
          f"rms_ref={rms_o:.3e}  rms_diff={rms_d:.3e}")

print(f"\nMean Stokes V rms diff : {np.mean(v_rms_diff):.3e}")
print(f"Mean Stokes V rms (new): {np.mean(v_rms_new):.3e}")
print(f"Mean Stokes V rms (ref): {np.mean(v_rms_old):.3e}")

# ── Plot ──────────────────────────────────────────────────────────────────────
ncols = 4
nrows = int(np.ceil(len(OBS_FILES) / ncols))
fig, axes = plt.subplots(nrows, ncols, figsize=(16, nrows * 3), sharex=False)
axes = axes.flatten()

for i, (obs_name, (nphase, nvel, nI, nV), (ophase, ovel, oI, oV),
        (avel, aI, aIsig, aV,
         aVsig)) in enumerate(zip(OBS_FILES, new_blocks, old_blocks,
                                  obs_data)):
    ax = axes[i]
    ax.set_title(f"phase = {nphase:.3f}", fontsize=8)
    ax.axhline(0, color="gray", lw=0.5, ls="--")
    # Observed data with error bars (plot every 2nd point to avoid clutter)
    if avel is not None:
        step = max(1, len(avel) // 50)
        ax.errorbar(avel[::step],
                    aV[::step] * 1e4,
                    yerr=aVsig[::step] * 1e4,
                    fmt="none",
                    ecolor="#999999",
                    elinewidth=0.8,
                    capsize=1.5,
                    zorder=1)
        ax.plot(avel,
                aV * 1e4,
                color="#888888",
                lw=0.8,
                alpha=0.7,
                label="Observed",
                zorder=2)
    ax.plot(nvel, nV * 1e4, "b-", lw=1.3, label="ZDIpy_WebUI", zorder=3)
    if ovel is not None:
        ax.plot(ovel,
                oV * 1e4,
                "r--",
                lw=1.3,
                label="ZDIpy-main ref",
                zorder=4)
    ax.set_xlabel("v (km/s)", fontsize=7)
    ax.set_ylabel(r"$V/I_c$ ($\times 10^{-4}$)", fontsize=7)
    ax.tick_params(labelsize=6)
    if i == 0:
        ax.legend(fontsize=6, loc="upper right")

for j in range(len(OBS_FILES), len(axes)):
    axes[j].set_visible(False)

fig.suptitle(
    "Stokes V profiles: observed (grey) / ZDIpy_WebUI (blue) / ZDIpy-main ref (red)",
    fontsize=11)
plt.tight_layout()
plt.savefig(PLOT_OUT, dpi=120)
print(f"\nComparison plot saved to: {PLOT_OUT}")
