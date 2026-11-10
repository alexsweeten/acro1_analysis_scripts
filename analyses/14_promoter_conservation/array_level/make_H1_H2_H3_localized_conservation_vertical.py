#!/usr/bin/env python3
"""Render the strict-majority H1-H3 array experiment as a vertical figure."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SOURCE_SCRIPT = Path(__file__).with_name("make_H1_H2_H3_localized_conservation.py")
TABLES = ROOT / "tables"
OUT = ROOT / "figures"


def load_source_module():
    spec = importlib.util.spec_from_file_location("h1_h2_h3_array_source", SOURCE_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {SOURCE_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def draw_profile(axis, source, base, profile: pd.DataFrame) -> None:
    local = profile[
        profile.alignment_column.between(base.PROFILE_START, base.PROFILE_END)
    ].copy()
    axis.axvspan(
        base.PROMOTER_START - 0.5,
        base.PROMOTER_END + 0.5,
        color="#F4D03F",
        alpha=0.25,
        zorder=0,
    )
    axis.axvline(
        base.PROMOTER_START - 0.5,
        color="#B7950B",
        linewidth=1.0,
        linestyle=":",
    )
    axis.axvline(
        base.PROMOTER_END + 0.5,
        color="#B7950B",
        linewidth=1.0,
        linestyle=":",
    )
    for group in source.GROUP_ORDER:
        subset = local[local.dominant_cluster.eq(group)]
        x = subset.alignment_column.to_numpy()
        mean = 100 * subset.individual_weighted_mean_rolling11_divergence.to_numpy()
        low = 100 * subset.ci95_low.to_numpy()
        high = 100 * subset.ci95_high.to_numpy()
        axis.fill_between(x, low, high, color=source.COLORS[group], alpha=0.12, linewidth=0)
        axis.plot(
            x,
            mean,
            color=source.COLORS[group],
            linewidth=2.5,
            label=source.LABELS[group],
        )

    upper = 1.12 * 100 * local.ci95_high.max()
    axis.set_ylim(0, upper)
    axis.set_xlim(base.PROFILE_START, base.PROFILE_END)
    axis.text(
        (base.PROMOTER_START + base.PROMOTER_END) / 2,
        upper * 0.965,
        "50-bp candidate interval",
        ha="center",
        va="top",
        fontsize=12.5,
        fontweight="bold",
        color="#6E5A00",
    )
    axis.set_xlabel("Human-only L1PA10 alignment column")
    axis.set_ylabel("Mean rolling pairwise divergence (%)")
    axis.set_title("Localized divergence profile", loc="left", fontweight="bold")
    axis.legend(frameon=False, loc="upper right", bbox_to_anchor=(0.99, 0.90))
    axis.grid(axis="y", color="#DDE2E7", linewidth=0.9)
    axis.set_axisbelow(True)
    axis.spines[["top", "right"]].set_visible(False)


def draw_effects(
    axis,
    source,
    individual: pd.DataFrame,
    summary: pd.DataFrame,
    pairwise: pd.DataFrame,
) -> None:
    rng = np.random.default_rng(20261108)
    positions = np.arange(len(source.GROUP_ORDER))
    for position, group in zip(positions, source.GROUP_ORDER):
        values = 100 * individual.loc[
            individual.dominant_cluster.eq(group), "mean_paired_effect"
        ].to_numpy()
        violin = axis.violinplot(
            values,
            positions=[position],
            widths=0.70,
            showmeans=False,
            showmedians=False,
            showextrema=False,
            bw_method=0.28,
        )
        for body in violin["bodies"]:
            body.set_facecolor(source.COLORS[group])
            body.set_edgecolor(source.COLORS[group])
            body.set_alpha(0.25)
            body.set_linewidth(1.1)
        jitter = rng.normal(0, 0.065, size=len(values))
        axis.scatter(
            position + jitter,
            values,
            s=12,
            alpha=0.26,
            color=source.COLORS[group],
            edgecolors="none",
            rasterized=True,
            zorder=2,
        )
        axis.boxplot(
            [values],
            positions=[position],
            widths=0.15,
            showfliers=False,
            patch_artist=True,
            medianprops={"color": "#202020", "linewidth": 1.8},
            boxprops={"facecolor": "white", "edgecolor": "#202020", "linewidth": 1.2},
            whiskerprops={"color": "#202020", "linewidth": 1.2},
            capprops={"color": "#202020", "linewidth": 1.2},
            zorder=3,
        )
        row = summary[summary.dominant_cluster.eq(group)].iloc[0]
        mean = 100 * row.mean_effect
        low = 100 * row.bootstrap_ci95_low
        high = 100 * row.bootstrap_ci95_high
        axis.errorbar(
            position,
            mean,
            yerr=[[mean - low], [high - mean]],
            fmt="D",
            markersize=7.6,
            markerfacecolor="#F2C14E",
            markeredgecolor="#171717",
            markeredgewidth=1.1,
            ecolor="#171717",
            elinewidth=1.4,
            capsize=4,
            zorder=4,
        )
        axis.text(
            position,
            mean + 0.24,
            f"{mean:.2f} pp",
            ha="center",
            va="bottom",
            fontsize=10.8,
            fontweight="bold",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78, "pad": 0.8},
        )

    axis.axhline(0, color="#222222", linewidth=1.1, linestyle="--", zorder=1)
    y_min = 100 * individual.mean_paired_effect.min()
    y_max = 100 * individual.mean_paired_effect.max()
    base_y = y_max + 0.28
    step = 0.53
    height = 0.12
    bracket_specs = [
        (0, 1, "H2_minus_H1", base_y),
        (1, 2, "H2_minus_H3", base_y + step),
        (0, 2, "H1_minus_H3", base_y + 2 * step),
    ]
    for x0, x1, comparison, y in bracket_specs:
        row = pairwise[pairwise.comparison.eq(comparison)].iloc[0]
        source.add_bracket(axis, x0, x1, y, height, f"q = {row.BH_q_value:.1e}")

    axis.set_ylim(y_min - 0.34, base_y + 2 * step + height + 0.48)
    axis.set_xlim(-0.55, 2.55)
    tick_labels = []
    for group in source.GROUP_ORDER:
        row = summary[summary.dominant_cluster.eq(group)].iloc[0]
        tick_labels.append(
            f"{source.LABELS[group]}\n"
            f"(n={int(row.individuals)} individuals)\n"
            f"{100 * row.fraction_positive:.1f}% positive"
        )
    axis.set_xticks(positions, tick_labels)
    axis.set_ylabel("Paired conservation effect\n(control - promoter, percentage points)")
    axis.set_title("Paired promoter-control effect", loc="left", fontweight="bold")
    axis.grid(axis="y", color="#DDE2E7", linewidth=0.9)
    axis.set_axisbelow(True)
    axis.spines[["top", "right"]].set_visible(False)


def main() -> None:
    source = load_source_module()
    base = source.load_base_module()
    base.configure_style()
    OUT.mkdir(parents=True, exist_ok=True)

    profile = pd.read_csv(TABLES / "H1_H2_H3_local_divergence_profile.tsv", sep="\t")
    individual = pd.read_csv(TABLES / "H1_H2_H3_individual_paired_effects.tsv.gz", sep="\t")
    summary = pd.read_csv(TABLES / "H1_H2_H3_effect_summary.tsv", sep="\t")
    summary = summary[summary.analysis.eq("strict_majority")].copy()
    pairwise = pd.read_csv(TABLES / "H1_H2_H3_pairwise_permutation_tests.tsv", sep="\t")

    figure = plt.figure(figsize=(10.8, 13.8))
    grid = figure.add_gridspec(2, 1, height_ratios=[1.0, 1.12], hspace=0.34)
    axis_profile = figure.add_subplot(grid[0])
    axis_effect = figure.add_subplot(grid[1])
    draw_profile(axis_profile, source, base, profile)
    draw_effects(axis_effect, source, individual, summary, pairwise)

    for label, axis in zip(("A", "B"), (axis_profile, axis_effect)):
        axis.text(
            -0.105,
            1.08,
            label,
            transform=axis.transAxes,
            fontsize=27,
            fontweight="bold",
            va="top",
            ha="left",
        )

    figure.suptitle(
        "Localized L1PA10 promoter conservation across H1-H3 arrays",
        fontsize=24,
        fontweight="bold",
        y=0.988,
    )
    figure.subplots_adjust(left=0.12, right=0.975, top=0.925, bottom=0.095)

    pdf = OUT / "Figure6_L1PA10_promoter_conservation_array_dominant_vertical.pdf"
    png = OUT / "Figure6_L1PA10_promoter_conservation_array_dominant_vertical.png"
    figure.savefig(pdf, bbox_inches="tight")
    figure.savefig(png, dpi=400, bbox_inches="tight", facecolor="white")
    plt.close(figure)
    print(pdf)
    print(png)


if __name__ == "__main__":
    main()
