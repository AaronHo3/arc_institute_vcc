"""Render the report figures from Arc's ST-HVG-Replogle paired evaluations.

Figure 1: fewshot vs zeroshot per metric, per held-out cell line (dumbbells).
Figure 2: the same contrast stratified by true effect size.

Usage:
    python make_figures.py [--evals-dir DIR] [--out-dir figures]
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from degradation_analysis import EFFECT_BINS, EFFECT_LABELS, LINES, load_paired

BLUE = "#2a78d6"    # fewshot (categorical slot 1; palette validated)
ORANGE = "#eb6834"  # zeroshot (slot 2)
INK = "#0b0b0b"
INK2 = "#52514e"
GRID = "#e4e3e0"

FIG1_METRICS = [
    ("pearson_delta", "Delta correlation", (0, 1), None),
    ("de_direction_match", "DE direction match", (0.4, 1), 0.5),
    ("overlap_at_100", "DE overlap@100", (0, 0.6), None),
    ("de_sig_genes_recall", "DE sig-gene recall", (0, 0.6), None),
    ("discrimination_score_l1", "Discrimination (L1)", (0.4, 1), 0.5),
    ("mae", "Mean absolute error", (0, 0.12), None),
]


def style_axis(ax):
    ax.grid(axis="x", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=9)


def fig_degradation(paired, out_dir: Path):
    fig, axes = plt.subplots(2, 3, figsize=(10.5, 5.6), constrained_layout=True)
    display = {"k562": "K562", "rpe1": "RPE1", "hepg2": "HepG2", "jurkat": "Jurkat"}
    names = [display[ln] for ln in LINES]
    ys = range(len(LINES))

    for ax, (col, title, xlim, chance) in zip(axes.flat, FIG1_METRICS):
        for y, line in zip(ys, LINES):
            fs, zs = paired[line]
            f, z = fs[col].mean(), zs[col].mean()
            ax.plot([z, f], [y, y], color=GRID, linewidth=2, zorder=1)
            ax.plot(z, y, "o", color=ORANGE, markersize=8, zorder=2)
            ax.plot(f, y, "o", color=BLUE, markersize=8, zorder=2)
        if chance is not None:
            ax.axvline(chance, color=INK2, linewidth=1, linestyle=":", zorder=0)
            ax.text(chance, len(LINES) - 0.35, " chance", color=INK2,
                    fontsize=8, va="bottom")
        ax.set_yticks(list(ys), names)
        ax.set_xlim(*xlim)
        ax.set_ylim(-0.6, len(LINES) - 0.4)
        ax.invert_yaxis()
        ax.set_title(title, fontsize=10, color=INK, loc="left")
        style_axis(ax)

    handles = [plt.Line2D([], [], marker="o", linestyle="", color=BLUE,
                          markersize=8, label="few-shot (cell type seen perturbed)"),
               plt.Line2D([], [], marker="o", linestyle="", color=ORANGE,
                          markersize=8, label="zero-shot (cell type never seen perturbed)")]
    fig.legend(handles=handles, loc="upper center", ncol=2, frameon=False,
               bbox_to_anchor=(0.5, 1.09), fontsize=9)
    fig.suptitle("STATE performance on held-out perturbations, by held-out cell line\n",
                 fontsize=11, color=INK, y=1.16)
    for ext in ("png", "svg"):
        fig.savefig(out_dir / f"fig1_degradation.{ext}", dpi=200,
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)


def fig_stratified(paired, out_dir: Path):
    frames = []
    for line in LINES:
        fs, zs = paired[line]
        frames.append(pd.DataFrame({
            "effect": pd.cut(zs["de_nsig_counts_real"], EFFECT_BINS,
                             labels=EFFECT_LABELS, right=False),
            "fs_pearson_delta": fs["pearson_delta"],
            "zs_pearson_delta": zs["pearson_delta"],
            "fs_dir": fs["de_direction_match"],
            "zs_dir": zs["de_direction_match"],
            "fs_ov": fs["overlap_at_100"],
            "zs_ov": zs["overlap_at_100"],
        }))
    d = pd.concat(frames)
    g = d.groupby("effect", observed=True).mean(numeric_only=True)
    n = d.groupby("effect", observed=True).size()
    x = range(len(EFFECT_LABELS))

    panels = [("Delta correlation", "fs_pearson_delta", "zs_pearson_delta", None),
              ("DE direction match", "fs_dir", "zs_dir", 0.5),
              ("DE overlap@100", "fs_ov", "zs_ov", None)]
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.4), constrained_layout=True)
    for ax, (title, fc, zc, chance) in zip(axes, panels):
        ax.plot(x, g[fc], "-o", color=BLUE, linewidth=2, markersize=8)
        ax.plot(x, g[zc], "-o", color=ORANGE, linewidth=2, markersize=8)
        if chance is not None:
            ax.axhline(chance, color=INK2, linewidth=1, linestyle=":")
        labels = [f"{lbl}\nn={n[lbl]}" for lbl in EFFECT_LABELS]
        ax.set_xticks(list(x), labels)
        ax.set_ylim(0, 1)
        ax.set_title(title, fontsize=10, color=INK, loc="left")
        ax.grid(axis="y", color=GRID, linewidth=0.8)
        ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color(GRID)
        ax.tick_params(colors=INK2, labelsize=9)
    axes[0].set_ylabel("metric value", color=INK2, fontsize=9)
    axes[1].set_xlabel("true effect size (n significant DE genes in real data)",
                       color=INK2, fontsize=9)
    handles = [plt.Line2D([], [], marker="o", color=BLUE, markersize=8,
                          linewidth=2, label="few-shot"),
               plt.Line2D([], [], marker="o", color=ORANGE, markersize=8,
                          linewidth=2, label="zero-shot")]
    fig.legend(handles=handles, loc="upper center", ncol=2, frameon=False,
               bbox_to_anchor=(0.5, 1.13), fontsize=9)
    for ext in ("png", "svg"):
        fig.savefig(out_dir / f"fig2_stratified.{ext}", dpi=200,
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--evals-dir", type=Path,
                    default=Path("C:/Users/aaron/vcc_data/st_hvg_replogle_evals"))
    ap.add_argument("--out-dir", type=Path, default=Path("figures"))
    args = ap.parse_args()
    args.out_dir.mkdir(exist_ok=True)
    paired = load_paired(args.evals_dir)
    fig_degradation(paired, args.out_dir)
    fig_stratified(paired, args.out_dir)
    print(f"wrote figures to {args.out_dir}/")


if __name__ == "__main__":
    main()
