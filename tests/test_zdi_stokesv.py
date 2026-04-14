"""Stokes V 振幅对照调试脚本。

基于同一份 ``config.json`` 运行多组参数案例，输出：

1) 每组的 ``chi2/dof``；
2) 每个相位的振幅比（max / p95 / rms）；
3) 自动诊断：振幅偏小主要来自 ``filling_factor_V`` 还是 ``num_iterations``。

Usage:
    python tests/test_zdi_stokesv.py
"""

import os
import sys

sys.path.insert(0,
                os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from config_loader import ZDIConfig
from pipeline.pipeline import ZDIPipeline


def _run_case(config_path,
              name,
              filling_factor_v,
              chi2_scale_v=None,
              num_iterations=None,
              halpha_auto_init=None):
    """Run one pipeline case with temporary parameter overrides."""
    par = ZDIConfig(config_path)
    par.halpha_filling_factor_V = float(filling_factor_v)
    if chi2_scale_v is not None:
        par.chiScaleV = float(chi2_scale_v)
    if num_iterations is not None:
        par.numIterations = int(num_iterations)
    if halpha_auto_init is not None:
        par.halpha_auto_init = bool(halpha_auto_init)

    print(f"\n=== Case: {name} ===")
    print(f"  halpha_filling_factor_V = {par.halpha_filling_factor_V}")
    print(f"  chi2_scale_V            = {getattr(par, 'chiScaleV', 1.0)}")
    print(f"  num_iterations          = {par.numIterations}")
    print(
        f"  halpha_auto_init        = {getattr(par, 'halpha_auto_init', True)}"
    )

    result = ZDIPipeline(par, verbose=0).run()

    obs_list = sorted(result.observed_profiles, key=lambda d: d["phase"])
    syn_list = sorted(result.synthetic_profiles, key=lambda d: d["phase"])

    amp_ratios_max = []
    amp_ratios_p95 = []
    amp_ratios_rms = []

    for obs, syn in zip(obs_list, syn_list):
        v_obs = np.asarray(obs["V_obs"])
        v_mod = np.asarray(syn["V_mod"])
        v_obs_abs = np.abs(v_obs)
        v_mod_abs = np.abs(v_mod)

        obs_max = float(np.max(v_obs_abs))
        mod_max = float(np.max(v_mod_abs))
        ratio_max = np.nan
        if obs_max > 0.0:
            ratio_max = mod_max / obs_max
        amp_ratios_max.append(ratio_max)

        obs_p95 = float(np.percentile(v_obs_abs, 95))
        mod_p95 = float(np.percentile(v_mod_abs, 95))
        ratio_p95 = np.nan
        if obs_p95 > 0.0:
            ratio_p95 = mod_p95 / obs_p95
        amp_ratios_p95.append(ratio_p95)

        obs_rms = float(np.sqrt(np.mean(v_obs_abs**2)))
        mod_rms = float(np.sqrt(np.mean(v_mod_abs**2)))
        ratio_rms = np.nan
        if obs_rms > 0.0:
            ratio_rms = mod_rms / obs_rms
        amp_ratios_rms.append(ratio_rms)

    amp_ratios_max = np.asarray(amp_ratios_max, dtype=float)
    amp_ratios_p95 = np.asarray(amp_ratios_p95, dtype=float)
    amp_ratios_rms = np.asarray(amp_ratios_rms, dtype=float)

    med_ratio_max = float(np.nanmedian(amp_ratios_max))
    med_ratio_p95 = float(np.nanmedian(amp_ratios_p95))
    med_ratio_rms = float(np.nanmedian(amp_ratios_rms))

    print(f"  Iterations             : {result.iterations}")
    print(f"  chi2/dof               : {result.chi2:.5f}")
    print(f"  Converged              : {result.converged}")
    print(f"  Median max-ratio       : {med_ratio_max:.4f}")
    print(f"  Median p95-ratio       : {med_ratio_p95:.4f}")
    print(f"  Median rms-ratio       : {med_ratio_rms:.4f}")

    return {
        "name": name,
        "filling_factor_v": float(par.halpha_filling_factor_V),
        "chi2_scale_v": float(getattr(par, "chiScaleV", 1.0)),
        "result": result,
        "obs_list": obs_list,
        "syn_list": syn_list,
        "median_ratio_max": med_ratio_max,
        "median_ratio_p95": med_ratio_p95,
        "median_ratio_rms": med_ratio_rms,
    }


def _plot_case(case_data, save_path):
    """Create stacked Stokes V plot for one case."""
    obs_list = case_data["obs_list"]
    syn_list = case_data["syn_list"]
    result = case_data["result"]

    n_phases = len(obs_list)
    all_v_obs = np.concatenate([np.asarray(d["V_obs"]) for d in obs_list])
    v_span = float(np.max(np.abs(all_v_obs))) * 2.0
    offset_step = max(v_span * 1.2, 0.002)

    fig, ax = plt.subplots(figsize=(8, 2.5 + n_phases * 0.9))

    for idx, (obs, syn) in enumerate(zip(obs_list, syn_list)):
        phase = obs["phase"]
        offset = idx * offset_step

        vel_obs = np.asarray(obs["vel"])
        v_obs = np.asarray(obs["V_obs"])
        v_sig = np.asarray(obs["V_sig"])

        vel_syn = np.asarray(syn["vel"])
        v_mod = np.asarray(syn["V_mod"])

        ax.fill_between(vel_obs,
                        offset + v_obs - v_sig,
                        offset + v_obs + v_sig,
                        color="steelblue",
                        alpha=0.25,
                        linewidth=0)
        ax.plot(vel_obs,
                offset + v_obs,
                color="steelblue",
                linewidth=1.0,
                label="Observed" if idx == 0 else "")

        ax.plot(vel_syn,
                offset + v_mod,
                color="tomato",
                linewidth=1.2,
                linestyle="--",
                label="Model" if idx == 0 else "")

        ax.text(vel_obs[-1] + 3,
                offset,
                f"phi = {phase:.4f}",
                va="center",
                ha="left",
                fontsize=7.5,
                color="0.3")

    ax.axhline(0, color="0.85", linewidth=0.5, linestyle=":")
    ax.set_xlabel("Velocity (km/s)")
    ax.set_ylabel("Stokes V / Ic (+ phase offset)")
    ax.set_title(f"{case_data['name']} ({n_phases} phases, "
                 f"chi2/dof = {result.chi2:.4f})")
    ax.legend(loc="upper left", fontsize=8)
    ax.set_xlim(vel_obs[0] - 5, vel_obs[-1] + 60)
    plt.tight_layout()

    fig.savefig(save_path, dpi=150)
    print(f"Plot saved to: {save_path}")
    plt.close(fig)


def run_and_plot():
    config_path = os.path.join(os.path.dirname(__file__), "../config.json")
    print("Running comparative Stokes V debug cases...")

    cases = [
        _run_case(config_path,
                  "baseline_fV0.5_it6",
                  filling_factor_v=0.5,
                  num_iterations=6),
        _run_case(config_path,
                  "test_fV1.0_it6",
                  filling_factor_v=1.0,
                  num_iterations=6),
        _run_case(config_path,
                  "test_fV1.0_chiV5_it6",
                  filling_factor_v=1.0,
                  chi2_scale_v=5.0,
                  num_iterations=6),
        _run_case(config_path,
                  "test_fV1.0_it20",
                  filling_factor_v=1.0,
                  num_iterations=20),
        _run_case(config_path,
                  "test_fV1.0_it80",
                  filling_factor_v=1.0,
                  num_iterations=80),
    ]

    print("\n=== Summary ===")
    for case in cases:
        print(f"  {case['name']}: max={case['median_ratio_max']:.4f}, "
              f"p95={case['median_ratio_p95']:.4f}, "
              f"rms={case['median_ratio_rms']:.4f} "
              f"(fV={case['filling_factor_v']}, "
              f"it={case['result'].iterations}, "
              f"chi2_scale_V={case['chi2_scale_v']})")

    case_map = {c["name"]: c for c in cases}
    base = case_map["baseline_fV0.5_it6"]
    f1 = case_map["test_fV1.0_it6"]
    i80 = case_map["test_fV1.0_it80"]

    d_fv = f1["median_ratio_max"] - base["median_ratio_max"]
    d_it = i80["median_ratio_max"] - f1["median_ratio_max"]
    d_it_p95 = i80["median_ratio_p95"] - f1["median_ratio_p95"]

    print("\n=== Auto Diagnosis ===")
    print(f"  fV提升(0.5->1.0, it=6) max-ratio增量: {d_fv:+.4f}")
    print(f"  迭代提升(it=6->80, fV=1.0) max-ratio增量: {d_it:+.4f}")
    print(f"  迭代提升(it=6->80, fV=1.0) p95-ratio增量: {d_it_p95:+.4f}")

    if d_it > 0.15:
        print("  [Diagnosis] 主要瓶颈是 num_iterations 过低，而非 filling_factor_V 本身。")
    else:
        print("  [Diagnosis] 主瓶颈不在迭代步数，需继续检查线型参数或速度中心。")

    out_dir = os.path.dirname(__file__)
    for case in cases:
        fig_name = f"test_zdi_stokesv_{case['name']}.png"
        _plot_case(case, os.path.join(out_dir, fig_name))


if __name__ == "__main__":
    run_and_plot()
