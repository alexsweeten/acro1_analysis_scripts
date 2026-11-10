#!/usr/bin/env python3
"""Build Figure 6 from the direct H1-H3 ACRE-level promoter analysis."""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "figures"
ANALYSIS = ROOT
GROUP_VALUES = ANALYSIS / "tables" / "array_H_cluster_promoter_metrics.tsv.gz"
SUMMARY = ANALYSIS / "tables" / "H_cluster_summary_statistics.tsv"
ONE_SAMPLE = ANALYSIS / "tables" / "H_cluster_effect_vs_zero_tests.tsv"
PAIRWISE = ANALYSIS / "tables" / "within_array_H_cluster_pairwise_tests.tsv"
PROFILE = ANALYSIS / "tables" / "H1_H2_H3_ACRE_level_localized_divergence_profile.tsv"

GROUPS = ["H1", "H2", "H3"]
CLUSTER_NUMBERS = {"H1": 1, "H2": 2, "H3": 3}
COLORS = {"H1": "#4E79A7", "H2": "#F28E2B", "H3": "#59A14F"}
PROMOTER_START = 74
PROMOTER_END = 123


def configure_style() -> None:
    helvetica = Path("/System/Library/Fonts/Helvetica.ttc")
    if helvetica.exists():
        font_manager.fontManager.addfont(str(helvetica))
    mpl.rcParams.update(
        {
            "font.family": "Helvetica, Arial, sans-serif",
            "font.size": 14,
            "axes.titlesize": 20,
            "axes.labelsize": 17,
            "xtick.labelsize": 13,
            "ytick.labelsize": 13,
            "legend.fontsize": 12.5,
            "axes.linewidth": 1.2,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def individual_means(values: pd.DataFrame, metric: str, group: str) -> np.ndarray:
    cluster = CLUSTER_NUMBERS[group]
    return (
        values.loc[values.cluster.eq(cluster), ["individual_id", metric]]
        .dropna()
        .groupby("individual_id", sort=False)[metric]
        .mean()
        .to_numpy(float)
    )


def localized_profile_panel(axis: plt.Axes, profile: pd.DataFrame) -> None:
    """Draw individual-balanced localized divergence for direct H1-H3 ACRE groups."""
    x_column = "display_alignment_column"
    mean_column = "individual_balanced_mean_rolling11_divergence"
    axis.axvspan(
        PROMOTER_START - 0.5,
        PROMOTER_END + 0.5,
        color="#F4D03F",
        alpha=0.26,
        zorder=0,
    )
    for group in GROUPS:
        subset = profile.loc[profile.cluster.eq(group)].sort_values(x_column)
        x = subset[x_column].to_numpy(float)
        mean = 100 * subset[mean_column].to_numpy(float)
        low = 100 * subset.ci95_low.to_numpy(float)
        high = 100 * subset.ci95_high.to_numpy(float)
        axis.fill_between(x, low, high, color=COLORS[group], alpha=0.16, linewidth=0)
        axis.plot(x, mean, color=COLORS[group], linewidth=2.7, label=group)

    axis.text(
        (PROMOTER_START + PROMOTER_END) / 2,
        0.95,
        "50-bp candidate interval",
        transform=axis.get_xaxis_transform(),
        ha="center",
        va="top",
        fontsize=12.5,
        fontweight="bold",
        color="#725A00",
    )
    axis.set_xlim(profile[x_column].min(), profile[x_column].max())
    maximum = 100 * profile.ci95_high.max()
    axis.set_ylim(0, np.ceil(maximum + 0.5))
    axis.set_xlabel("Human-only L1PA10 alignment column")
    axis.set_ylabel("Mean rolling pairwise divergence (%)")
    axis.set_title("Localized divergence profile", loc="left", fontweight="bold")
    axis.legend(
        title="NTRiage cluster",
        frameon=False,
        loc="upper right",
        fontsize=12.5,
        title_fontsize=12.5,
    )
    axis.grid(axis="y", color="#DDE2E7", linewidth=0.9)
    axis.set_axisbelow(True)
    axis.spines[["top", "right"]].set_visible(False)


def violin_panel(
    axis: plt.Axes,
    values: pd.DataFrame,
    summary: pd.DataFrame,
    tests: pd.DataFrame,
) -> None:
    rng = np.random.default_rng(20261108)
    positions = np.arange(1, 4)
    estimable_groups = {"H1": 548, "H2": 1177, "H3": 750}
    estimable_people = {"H1": 206, "H2": 230, "H3": 216}
    estimable_fraction = {"H1": 30.9, "H2": 100.0, "H3": 99.5}

    for position, group in zip(positions, GROUPS):
        data = 100 * individual_means(values, "local_effect", group)
        violin = axis.violinplot(
            data,
            positions=[position],
            widths=0.82,
            showmeans=False,
            showmedians=False,
            showextrema=False,
            bw_method=0.24,
        )
        for body in violin["bodies"]:
            body.set_facecolor(COLORS[group])
            body.set_edgecolor(COLORS[group])
            body.set_alpha(0.42)
            body.set_linewidth(1.1)

        jitter = rng.uniform(-0.14, 0.14, len(data))
        axis.scatter(
            np.full(len(data), position) + jitter,
            data,
            s=14,
            alpha=0.24,
            color=COLORS[group],
            edgecolors="none",
            rasterized=True,
            zorder=2,
        )
        axis.boxplot(
            [data],
            positions=[position],
            widths=0.16,
            patch_artist=True,
            showfliers=False,
            boxprops={"facecolor": "white", "edgecolor": "#202020", "linewidth": 1.0},
            medianprops={"color": "#202020", "linewidth": 1.4},
            whiskerprops={"color": "#202020", "linewidth": 1.0},
            capprops={"color": "#202020", "linewidth": 1.0},
            zorder=3,
        )

        row = summary[
            summary.cohort.eq("all_arrays")
            & summary.cluster.eq(group)
            & summary.metric.eq("local_effect")
        ].iloc[0]
        mean = 100 * float(row.individual_weighted_mean)
        low = 100 * float(row.individual_boot_ci95_lo)
        high = 100 * float(row.individual_boot_ci95_hi)
        axis.errorbar(
            position,
            mean,
            yerr=np.array([[mean - low], [high - mean]]),
            fmt="D",
            markersize=8.2,
            markerfacecolor="#F2C14E",
            markeredgecolor="#161616",
            markeredgewidth=1.0,
            ecolor="#161616",
            elinewidth=1.4,
            capsize=4,
            zorder=5,
        )
        offset = 0.22 if group != "H3" else -0.34
        va = "bottom" if group != "H3" else "top"
        axis.text(
            position,
            mean + offset,
            f"{mean:+.2f}",
            ha="center",
            va=va,
            fontsize=11.5,
            fontweight="bold",
            color="#171717",
            zorder=6,
        )

    q_lookup = tests[
        tests.cohort.eq("all_arrays")
        & tests.metric.eq("local_effect")
    ].set_index("cluster").q_bh_within_cohort_metric
    axis.text(
        2,
        3.05,
        r"H1 and H2: $q = 7.5 \times 10^{-6}$; H3: opposite direction",
        ha="center",
        va="bottom",
        fontsize=10.5,
        fontweight="bold",
    )

    axis.axhline(0, color="#202020", linewidth=1.1, linestyle="--", zorder=1)
    axis.set_ylim(-3.55, 3.42)
    axis.set_xlim(0.45, 3.55)
    axis.set_xticks(
        positions,
        [
            f"{group}\n{estimable_groups[group]:,} groups\n{estimable_people[group]} people"
            for group in GROUPS
        ],
    )
    axis.set_ylabel("Matched-divergence effect\n(control - candidate, percentage points)")
    axis.set_title("Within-cluster matched divergence", loc="left", fontweight="bold")
    axis.text(
        0.01,
        0.99,
        "Matched-control availability: "
        + ", ".join(f"{g} {estimable_fraction[g]:g}%" for g in GROUPS),
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=10.5,
        color="#4B5563",
    )
    axis.grid(axis="y", color="#DDE2E7", linewidth=0.9)
    axis.set_axisbelow(True)
    axis.spines[["top", "right"]].set_visible(False)


def add_bracket(axis: plt.Axes, x0: float, x1: float, y: float, height: float, label: str) -> None:
    axis.plot(
        [x0, x0, x1, x1],
        [y, y + height, y + height, y],
        color="#202020",
        linewidth=1.15,
        clip_on=False,
    )
    axis.text((x0 + x1) / 2, y + height + 0.8, label, ha="center", va="bottom", fontsize=11.5)


def normalized_panel(axis: plt.Axes, summary: pd.DataFrame, pairwise: pd.DataFrame) -> None:
    metric = "promoter_depletion_scaled_by_l1pa10_nonpromoter_divergence"
    positions = [1, 2]
    groups = ["H1", "H2"]
    for position, group in zip(positions, groups):
        row = summary[
            summary.cohort.eq("all_arrays")
            & summary.cluster.eq(group)
            & summary.metric.eq(metric)
        ].iloc[0]
        mean = 100 * float(row.individual_weighted_mean)
        low = 100 * float(row.individual_boot_ci95_lo)
        high = 100 * float(row.individual_boot_ci95_hi)
        axis.errorbar(
            position,
            mean,
            yerr=np.array([[mean - low], [high - mean]]),
            fmt="D",
            markersize=13,
            markerfacecolor=COLORS[group],
            markeredgecolor="#161616",
            markeredgewidth=1.1,
            ecolor=COLORS[group],
            elinewidth=2.3,
            capsize=7,
            capthick=1.7,
            zorder=4,
        )
        axis.text(
            position,
            high + 2.2,
            f"{mean:.1f}%",
            ha="center",
            va="bottom",
            fontsize=14,
            fontweight="bold",
        )

    test = pairwise[
        pairwise.cohort.eq("all_arrays")
        & pairwise.metric.eq(metric)
        & pairwise.cluster_a.eq("H1")
        & pairwise.cluster_b.eq("H2")
    ].iloc[0]
    add_bracket(
        axis,
        1,
        2,
        71.5,
        1.6,
        rf"paired $q$ = {float(test.q_bh_within_cohort_metric):.3f}",
    )
    axis.axhline(0, color="#202020", linewidth=1.0, linestyle="--")
    axis.set_xlim(0.45, 2.55)
    axis.set_ylim(0, 79)
    axis.set_xticks(positions, ["H1\n231 people", "H2\n229 people"])
    axis.set_ylabel("Fractional promoter depletion (%)")
    axis.set_title("Background-scaled H1-H2 comparison", loc="left", fontweight="bold")
    axis.text(
        0.5,
        0.04,
        "Paired test: 586 arrays from 219 individuals",
        transform=axis.transAxes,
        ha="center",
        va="bottom",
        fontsize=10.5,
        color="#4B5563",
    )
    axis.grid(axis="y", color="#DDE2E7", linewidth=0.9)
    axis.set_axisbelow(True)
    axis.spines[["top", "right"]].set_visible(False)


def main() -> None:
    configure_style()
    OUT.mkdir(parents=True, exist_ok=True)

    profile = pd.read_csv(PROFILE, sep="\t")
    values = pd.read_csv(GROUP_VALUES, sep="\t")
    summary = pd.read_csv(SUMMARY, sep="\t")
    one_sample = pd.read_csv(ONE_SAMPLE, sep="\t")
    pairwise = pd.read_csv(PAIRWISE, sep="\t")

    figure = plt.figure(figsize=(19.2, 8.9))
    outer = figure.add_gridspec(1, 3, width_ratios=[1.58, 0.92, 0.72], wspace=0.26)
    axis_profile = figure.add_subplot(outer[0])
    axis_effect = figure.add_subplot(outer[1])
    axis_scaled = figure.add_subplot(outer[2])

    localized_profile_panel(axis_profile, profile)
    violin_panel(axis_effect, values, summary, one_sample)
    normalized_panel(axis_scaled, summary, pairwise)

    axis_profile.text(
        -0.13, 1.08, "A", transform=axis_profile.transAxes,
        fontsize=28, fontweight="bold", va="top",
    )
    axis_effect.text(
        -0.17, 1.08, "B", transform=axis_effect.transAxes,
        fontsize=28, fontweight="bold", va="top",
    )
    axis_scaled.text(
        -0.22, 1.08, "C", transform=axis_scaled.transAxes,
        fontsize=28, fontweight="bold", va="top",
    )

    figure.suptitle(
        "Candidate L1PA10 promoter conservation across H1-H3 ACREs",
        fontsize=26,
        fontweight="bold",
        y=0.992,
    )
    figure.subplots_adjust(left=0.064, right=0.988, top=0.875, bottom=0.13)

    pdf = OUT / "Figure6_L1PA10_promoter_conservation.pdf"
    png = OUT / "Figure6_L1PA10_promoter_conservation.png"
    figure.savefig(pdf, bbox_inches="tight")
    figure.savefig(png, dpi=400, bbox_inches="tight", facecolor="white")
    plt.close(figure)
    print(pdf)
    print(png)


if __name__ == "__main__":
    main()
