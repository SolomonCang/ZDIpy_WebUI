#!/usr/bin/env python3
"""
scripts/compare_voigt_vs_ur.py
对比 Voigt 弱场模型与 Unno-Rachkovsky 模型（fI=1, fV=0.16）的 ZDI 反演结果。

用法:
    cd /path/to/ZDIpy_WebUI
    python scripts/compare_voigt_vs_ur.py [--iterations N]

输出:
    - results/voigt/         : Voigt 运行输出文件
    - results/ur_fV016/      : UR 运行输出文件（fI=1, fV=0.16）
    - compare_voigt_vs_ur.png: 合成轮廓对比图（Stokes I 和 V）
    - compare_voigt_vs_ur_summary.txt: 数值指标对比摘要
"""

import argparse
import copy
import json
import logging
import os
import sys
import tempfile
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ── 路径设置 ─────────────────────────────────────────────────────────────────
ROOT = Path(__file__).parent.parent.resolve()
sys.path.insert(0, str(ROOT))

CONFIG_PATH = ROOT / "config.json"
OUT_VOIGT = ROOT / "results" / "voigt"
OUT_UR = ROOT / "results" / "ur_fV016"
PLOT_OUT = ROOT / "compare_voigt_vs_ur.png"
SUMMARY_OUT = ROOT / "compare_voigt_vs_ur_summary.txt"

# 配置 logging：只显示 zdipy 相关日志
logging.basicConfig(
    level=logging.WARNING,
    format="%(message)s",
)
logging.getLogger("zdipy").setLevel(logging.INFO)
logging.getLogger("zdipy.pipeline").setLevel(logging.INFO)

_fh = logging.StreamHandler(sys.stdout)
_fh.setLevel(logging.INFO)
logging.getLogger("zdipy.pipeline").addHandler(_fh)

# ── 辅助函数 ─────────────────────────────────────────────────────────────────


def _load_base_config() -> dict:
    """从 config.json 读取基础配置字典（深拷贝，不会污染原始对象）。"""
    with open(CONFIG_PATH, encoding="utf-8") as f:
        return json.load(f)


def _make_voigt_config(base: dict, iterations: int) -> dict:
    """克隆基础配置，切换为 Voigt 弱场模型，重写输出路径。"""
    cfg = copy.deepcopy(base)
    cfg["line_model"]["model_type"] = "voigt"
    cfg["inversion"]["num_iterations"] = iterations
    _set_output_paths(cfg, str(OUT_VOIGT))
    return cfg


def _make_ur_config(base: dict,
                    iterations: int,
                    f_I: float = 1.0,
                    f_V: float = 0.16) -> dict:
    """克隆基础配置，切换为 UR 模型并设定 fI, fV，重写输出路径。"""
    cfg = copy.deepcopy(base)
    cfg["line_model"]["model_type"] = "unno"
    cfg["line_model"]["unno_filling_factor_I"] = f_I
    cfg["line_model"]["unno_filling_factor_V"] = f_V
    cfg["inversion"]["num_iterations"] = iterations
    _set_output_paths(cfg, str(OUT_UR))
    return cfg


def _set_output_paths(cfg: dict, out_dir: str) -> None:
    """将 output 块中的所有路径重写到指定目录（保留文件名）。"""
    os.makedirs(out_dir, exist_ok=True)
    for key, val in cfg.get("output", {}).items():
        fname = Path(val).name
        cfg["output"][key] = str(Path(out_dir) / fname)


def _run_pipeline(cfg_dict: dict, label: str):
    """
    根据配置字典构造 ZDIConfig 并运行完整反演。

    Returns
    -------
    ZDIResult
    """
    import config_loader as cl
    from pipeline.pipeline import ZDIPipeline

    print(f"\n{'='*60}")
    print(f"  运行: {label}")
    print(f"{'='*60}")
    t0 = time.time()

    # 写临时配置文件到项目根目录，确保相对路径（LSDprof/, results/ 等）解析正确
    with tempfile.NamedTemporaryFile(mode="w",
                                     suffix=".json",
                                     delete=False,
                                     dir=str(ROOT),
                                     encoding="utf-8") as tf:
        json.dump(cfg_dict, tf, indent=2, ensure_ascii=False)
        tmp_path = tf.name
    try:
        par = cl.ZDIConfig(tmp_path)
        pipeline = ZDIPipeline(par, verbose=1)
        result = pipeline.run()
    finally:
        os.unlink(tmp_path)

    elapsed = time.time() - t0
    print(f"  >> 完成，耗时 {elapsed:.1f} 秒")
    return result


def _mag_energy(result) -> dict:
    """从 result.metadata 中提取磁场能量分析项（若已计算）。"""
    return result.metadata.get("mag_energy", {})


def _mean_mag(result) -> float:
    """从 metadata 中提取平均磁场强度（单位 G）。"""
    return float(result.metadata.get("mean_mag", np.nan))


def _chi2_I_V(result) -> tuple[float, float]:
    """从 metadata 中获取 chi2_I 和 chi2_V（若单独存储）。"""
    meta = result.metadata
    chi2_I = float(meta.get("chi2_I", np.nan))
    chi2_V = float(meta.get("chi2_V", np.nan))
    return chi2_I, chi2_V


# ── 绘图 ─────────────────────────────────────────────────────────────────────


def _plot_comparison(res_voigt, res_ur, f_V: float = 0.16):
    """
    绘制 Stokes I 和 Stokes V 合成轮廓对比图（Voigt vs UR）。
    每行一个观测相位，左列 I，右列 V。
    """
    n_obs = len(res_voigt.synthetic_profiles)
    if n_obs == 0:
        print("警告: 无合成轮廓数据，跳过绘图")
        return

    # 限制最多 16 个相位以保持图幅合理
    n_show = min(n_obs, 16)
    ncols = 2
    nrows = n_show

    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(10, max(nrows * 1.2, 4)),
        sharex="col",
        gridspec_kw={
            "hspace": 0.08,
            "wspace": 0.3
        },
    )
    if nrows == 1:
        axes = axes[np.newaxis, :]  # 保证二维索引

    ur_label = f"UR (fI=1, fV={f_V})"
    col_labels = ["Stokes I / Ic", "Stokes V / Ic (%)"]

    for col, col_label in enumerate(col_labels):
        axes[0, col].set_title(col_label, fontsize=9, pad=4)

    for row in range(n_show):
        sp_v = res_voigt.synthetic_profiles[row]
        sp_u = res_ur.synthetic_profiles[row]
        op = res_voigt.observed_profiles[row]

        vel_v = np.asarray(sp_v["vel"])
        vel_u = np.asarray(sp_u["vel"])
        vel_o = np.asarray(op["vel"])

        phase_str = f"φ={sp_v['phase']:.4f}"

        # --- Stokes I ---
        ax_I = axes[row, 0]
        ax_I.plot(vel_o,
                  op["I_obs"],
                  "k.",
                  ms=1.5,
                  alpha=0.6,
                  label="obs" if row == 0 else None)
        ax_I.plot(vel_v,
                  sp_v["I_mod"],
                  "b-",
                  lw=1.0,
                  label="Voigt" if row == 0 else None)
        ax_I.plot(vel_u,
                  sp_u["I_mod"],
                  "r--",
                  lw=1.0,
                  label=ur_label if row == 0 else None)
        ax_I.set_xlim(-200, 200)
        ax_I.text(0.97,
                  0.80,
                  phase_str,
                  transform=ax_I.transAxes,
                  ha="right",
                  va="top",
                  fontsize=6,
                  color="0.4")
        ax_I.tick_params(labelsize=6)
        ax_I.yaxis.set_tick_params(labelleft=True)

        # --- Stokes V ---
        ax_V = axes[row, 1]
        V_obs = np.asarray(op["V_obs"]) * 100
        V_voi = np.asarray(sp_v["V_mod"]) * 100
        V_ur = np.asarray(sp_u["V_mod"]) * 100
        ax_V.plot(vel_o, V_obs, "k.", ms=1.5, alpha=0.6)
        ax_V.plot(vel_v, V_voi, "b-", lw=1.0)
        ax_V.plot(vel_u, V_ur, "r--", lw=1.0)
        ax_V.axhline(0, lw=0.4, color="0.6")
        ax_V.set_xlim(-200, 200)
        ax_V.tick_params(labelsize=6)

    # 共享 x 轴标签
    for col in range(ncols):
        axes[-1, col].set_xlabel("v (km/s)", fontsize=8)

    # 图例放在右上角子图外
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles,
               labels,
               loc="upper right",
               fontsize=7,
               framealpha=0.9,
               bbox_to_anchor=(1.0, 1.0))

    fig.suptitle(
        f"ZDI Profile Comparison: Voigt vs UR (fI=1, fV={f_V})\n"
        f"Voigt: chi2/dof={res_voigt.chi2:.4f}   UR: chi2/dof={res_ur.chi2:.4f}",
        fontsize=9,
    )
    fig.savefig(PLOT_OUT, dpi=150, bbox_inches="tight")
    print(f"\n对比图已保存: {PLOT_OUT}")


# ── 摘要输出 ──────────────────────────────────────────────────────────────────


def _print_and_save_summary(res_v, res_u, f_I: float, f_V: float):
    """打印并保存数值对比摘要。"""
    en_v = _mag_energy(res_v)
    en_u = _mag_energy(res_u)

    lines = []
    lines.append("=" * 60)
    lines.append("ZDI 反演对比摘要: Voigt  vs  UR (fI={}, fV={})".format(f_I, f_V))
    lines.append("=" * 60)
    lines.append(f"{'指标':<28}{'Voigt':>12}{'UR':>12}")
    lines.append("-" * 52)

    def _row(label, v1, v2, fmt="{:.5f}"):
        s1 = fmt.format(v1) if not np.isnan(v1) else "  N/A"
        s2 = fmt.format(v2) if not np.isnan(v2) else "  N/A"
        lines.append(f"  {label:<26}{s1:>12}{s2:>12}")

    _row("迭代次数",
         float(res_v.iterations),
         float(res_u.iterations),
         fmt="{:.0f}")
    _row("熵 S", res_v.entropy, res_u.entropy)
    _row("χ²/dof", res_v.chi2, res_u.chi2)
    _row("test 收敛量", res_v.test, res_u.test)
    _row("收敛?", float(res_v.converged), float(res_u.converged), fmt="{:.0f}")
    _row("平均亮度", float(res_v.metadata.get("mean_bright", np.nan)),
         float(res_u.metadata.get("mean_bright", np.nan)))
    _row("平均 |B| (G)", _mean_mag(res_v), _mean_mag(res_u), fmt="{:.2f}")

    if en_v and en_u:
        lines.append("")
        lines.append("  -- 磁场能量分析 --")
        _row("极向场占比 (%)",
             en_v.get("pct_pol", np.nan) * 100,
             en_u.get("pct_pol", np.nan) * 100,
             fmt="{:.1f}")
        _row("环向场占比 (%)",
             en_v.get("pct_tor", np.nan) * 100,
             en_u.get("pct_tor", np.nan) * 100,
             fmt="{:.1f}")
        _row("轴对称占比 (%)",
             en_v.get("pct_axi", np.nan) * 100,
             en_u.get("pct_axi", np.nan) * 100,
             fmt="{:.1f}")

    lines.append("=" * 60)

    summary_str = "\n".join(lines)
    print("\n" + summary_str)

    with open(SUMMARY_OUT, "w", encoding="utf-8") as f:
        f.write(summary_str + "\n")
    print(f"\n摘要已保存: {SUMMARY_OUT}")


# ── CLI 入口 ──────────────────────────────────────────────────────────────────


def main():
    parser = argparse.ArgumentParser(
        description="对比 Voigt 与 UR (fI=1, fV=0.16) 模式 ZDI 反演结果")
    parser.add_argument(
        "--iterations",
        "-n",
        type=int,
        default=None,
        help="MEM 迭代次数（默认使用 config.json 中的 num_iterations）",
    )
    parser.add_argument(
        "--fV",
        type=float,
        default=0.16,
        help="UR 模式 Stokes V 填充因子（默认 0.16）",
    )
    parser.add_argument(
        "--fI",
        type=float,
        default=1.0,
        help="UR 模式 Stokes I 填充因子（默认 1.0）",
    )
    args = parser.parse_args()

    base = _load_base_config()
    default_iters = int(base["inversion"]["num_iterations"])
    iterations = args.iterations if args.iterations is not None else default_iters

    print(f"迭代次数: {iterations}")
    print(f"UR 参数: fI={args.fI}, fV={args.fV}")

    # --- 运行两个模型 ---
    cfg_voigt = _make_voigt_config(base, iterations)
    cfg_ur = _make_ur_config(base, iterations, f_I=args.fI, f_V=args.fV)

    res_voigt = _run_pipeline(cfg_voigt, label=f"Voigt  (弱场近似)")
    res_ur = _run_pipeline(cfg_ur,
                           label=f"UR     (fI={args.fI}, fV={args.fV})")

    # --- 对比输出 ---
    _print_and_save_summary(res_voigt, res_ur, f_I=args.fI, f_V=args.fV)
    _plot_comparison(res_voigt, res_ur, f_V=args.fV)

    print("\n全部完成。")


if __name__ == "__main__":
    main()
