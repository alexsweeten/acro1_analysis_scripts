#!/usr/bin/env python3
"""Compare localized L1PA10 promoter conservation across H1-H3 dominant arrays."""

from __future__ import annotations

import importlib.util
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "figures"
TABLES = ROOT / "tables"
BASE_SCRIPT = Path(
    os.environ.get("ACRO1_PROMOTER_ARRAY_BASE", ROOT / "base_candidate_analysis.py")
)
if not BASE_SCRIPT.exists():
    raise FileNotFoundError(f"Could not find the candidate analysis base: {BASE_SCRIPT}")

GROUP_ORDER = ["H1", "H2", "H3"]
COLORS = {"H1": "#4E79A7", "H2": "#F28E2B", "H3": "#59A14F"}
LABELS = {group: f"{group}-dominant" for group in GROUP_ORDER}
N_PERMUTATIONS = 500_000


def load_base_module():
    spec = importlib.util.spec_from_file_location("h2_h3_base", BASE_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {BASE_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def summarize_effects(base, frame: pd.DataFrame, analysis: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    individual = base.individual_effects(frame)
    rows = []
    for group in GROUP_ORDER:
        arrays = frame[frame.dominant_cluster.eq(group)]
        values = individual.loc[
            individual.dominant_cluster.eq(group), "mean_paired_effect"
        ].to_numpy(float)
        low, high = base.S.bootstrap_mean_ci(values, f"H1-H3|{analysis}|{group}")
        rows.append(
            {
                "analysis": analysis,
                "dominant_cluster": group,
                "display_group": LABELS[group],
                "arrays": len(arrays),
                "individuals": len(values),
                "median_dominant_fraction": arrays.dominant_fraction.median(),
                "mean_effect": values.mean(),
                "median_effect": np.median(values),
                "fraction_positive": (values > 0).mean(),
                "bootstrap_ci95_low": low,
                "bootstrap_ci95_high": high,
            }
        )
    return individual, pd.DataFrame(rows)


def pairwise_tests(base, frame: pd.DataFrame, individual: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    comparisons = [
        ("H2", "H1"),
        ("H1", "H3"),
        ("H2", "H3"),
    ]
    permutation_rows = []
    paired_rows = []
    for group1, group0 in comparisons:
        subset = frame[frame.dominant_cluster.isin([group1, group0])].copy()
        group = subset.dominant_cluster.eq(group1).to_numpy()
        strata = pd.factorize(subset.individual_id)[0]
        result = base.S.perm_diff_means(
            subset.effect.to_numpy(float),
            group,
            N_PERMUTATIONS,
            f"H1-H3-dominant|{group1}-minus-{group0}",
            strata=strata,
        )
        permutation_rows.append(
            {
                "comparison": f"{group1}_minus_{group0}",
                "group1": group1,
                "group0": group0,
                "arrays_group1": int(group.sum()),
                "arrays_group0": int((~group).sum()),
                "difference_in_array_means": result["diff"],
                "p_two_sided": result["p_two_sided"],
                "p_group1_greater": result["p_group1_greater"],
                "permutations": result["draws"],
            }
        )

        wide = individual[
            individual.dominant_cluster.isin([group1, group0])
        ].pivot(index="individual_id", columns="dominant_cluster", values="mean_paired_effect")
        wide = wide.dropna(subset=[group1, group0])
        differences = wide[group1] - wide[group0]
        paired_rows.append(
            {
                "comparison": f"{group1}_minus_{group0}",
                "group1": group1,
                "group0": group0,
                "paired_individuals": len(wide),
                "mean_paired_difference": differences.mean(),
                "median_paired_difference": differences.median(),
                "paired_wilcoxon_p": wilcoxon(differences).pvalue,
            }
        )

    permutations = pd.DataFrame(permutation_rows)
    permutations["BH_q_value"] = base.S.bh(permutations.p_two_sided.to_numpy())
    return permutations, pd.DataFrame(paired_rows)


def add_bracket(axis, x0: float, x1: float, y: float, height: float, label: str) -> None:
    axis.plot(
        [x0, x0, x1, x1],
        [y, y + height, y + height, y],
        color="#222222",
        linewidth=1.05,
        clip_on=False,
    )
    axis.text(
        (x0 + x1) / 2,
        y + height + 0.045,
        label,
        ha="center",
        va="bottom",
        fontsize=10.5,
        fontweight="bold",
    )


def make_figure(
    base,
    profile: pd.DataFrame,
    individual: pd.DataFrame,
    summary: pd.DataFrame,
    pairwise: pd.DataFrame,
) -> tuple[Path, Path]:
    figure, (axis_profile, axis_effect) = plt.subplots(
        1,
        2,
        figsize=(15.4, 6.9),
        gridspec_kw={"width_ratios": [1.36, 1.08]},
    )

    local = profile[
        profile.alignment_column.between(base.PROFILE_START, base.PROFILE_END)
    ].copy()
    axis_profile.axvspan(
        base.PROMOTER_START - 0.5,
        base.PROMOTER_END + 0.5,
        color="#F4D03F",
        alpha=0.25,
        zorder=0,
    )
    axis_profile.axvline(base.PROMOTER_START - 0.5, color="#B7950B", linewidth=1.0, linestyle=":")
    axis_profile.axvline(base.PROMOTER_END + 0.5, color="#B7950B", linewidth=1.0, linestyle=":")

    for group in GROUP_ORDER:
        subset = local[local.dominant_cluster.eq(group)]
        x = subset.alignment_column.to_numpy()
        mean = 100 * subset.individual_weighted_mean_rolling11_divergence.to_numpy()
        low = 100 * subset.ci95_low.to_numpy()
        high = 100 * subset.ci95_high.to_numpy()
        axis_profile.fill_between(x, low, high, color=COLORS[group], alpha=0.12, linewidth=0)
        axis_profile.plot(x, mean, color=COLORS[group], linewidth=2.35, label=LABELS[group])

    upper = 1.12 * 100 * local.ci95_high.max()
    axis_profile.set_ylim(0, upper)
    axis_profile.set_xlim(base.PROFILE_START, base.PROFILE_END)
    axis_profile.text(
        (base.PROMOTER_START + base.PROMOTER_END) / 2,
        upper * 0.965,
        "50-bp candidate interval",
        ha="center",
        va="top",
        fontsize=12,
        fontweight="bold",
        color="#6E5A00",
    )
    axis_profile.set_xlabel("Human-only L1PA10 alignment column")
    axis_profile.set_ylabel("Mean rolling pairwise divergence (%)")
    axis_profile.set_title("Localized divergence profile", loc="left", fontweight="bold")
    axis_profile.legend(frameon=False, loc="upper right", bbox_to_anchor=(0.99, 0.90))
    axis_profile.grid(axis="y", color="#DDE2E7", linewidth=0.9)
    axis_profile.set_axisbelow(True)
    axis_profile.spines[["top", "right"]].set_visible(False)

    rng = np.random.default_rng(20261108)
    positions = np.arange(len(GROUP_ORDER))
    for position, group in zip(positions, GROUP_ORDER):
        values = 100 * individual.loc[
            individual.dominant_cluster.eq(group), "mean_paired_effect"
        ].to_numpy()
        violin = axis_effect.violinplot(
            values,
            positions=[position],
            widths=0.70,
            showmeans=False,
            showmedians=False,
            showextrema=False,
            bw_method=0.28,
        )
        for body in violin["bodies"]:
            body.set_facecolor(COLORS[group])
            body.set_edgecolor(COLORS[group])
            body.set_alpha(0.25)
            body.set_linewidth(1.1)
        jitter = rng.normal(0, 0.065, size=len(values))
        axis_effect.scatter(
            position + jitter,
            values,
            s=12,
            alpha=0.26,
            color=COLORS[group],
            edgecolors="none",
            rasterized=True,
            zorder=2,
        )
        axis_effect.boxplot(
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
        axis_effect.errorbar(
            position,
            mean,
            yerr=[[mean - low], [high - mean]],
            fmt="D",
            markersize=7.2,
            markerfacecolor="#F2C14E",
            markeredgecolor="#171717",
            markeredgewidth=1.1,
            ecolor="#171717",
            elinewidth=1.4,
            capsize=4,
            zorder=4,
        )
        axis_effect.text(
            position,
            mean + 0.24,
            f"{mean:.2f} pp",
            ha="center",
            va="bottom",
            fontsize=10.8,
            fontweight="bold",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78, "pad": 0.8},
        )

    axis_effect.axhline(0, color="#222222", linewidth=1.1, linestyle="--", zorder=1)
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
        add_bracket(axis_effect, x0, x1, y, height, f"q = {row.BH_q_value:.1e}")

    axis_effect.set_ylim(y_min - 0.34, base_y + 2 * step + height + 0.48)
    axis_effect.set_xlim(-0.55, 2.55)
    tick_labels = []
    for group in GROUP_ORDER:
        row = summary[summary.dominant_cluster.eq(group)].iloc[0]
        tick_labels.append(
            f"{LABELS[group]}\n(n={int(row.individuals)} individuals)\n{100 * row.fraction_positive:.1f}% positive"
        )
    axis_effect.set_xticks(positions, tick_labels)
    axis_effect.set_ylabel(
        "Paired conservation effect\n(control - promoter, percentage points)"
    )
    axis_effect.set_title("Paired promoter-control effect", loc="left", fontweight="bold")
    axis_effect.grid(axis="y", color="#DDE2E7", linewidth=0.9)
    axis_effect.set_axisbelow(True)
    axis_effect.spines[["top", "right"]].set_visible(False)

    axis_profile.text(
        -0.13, 1.08, "A", transform=axis_profile.transAxes,
        fontsize=25, fontweight="bold", va="top", ha="left"
    )
    axis_effect.text(
        -0.13, 1.08, "B", transform=axis_effect.transAxes,
        fontsize=25, fontweight="bold", va="top", ha="left"
    )
    figure.suptitle(
        "Localized L1PA10 promoter conservation across H1-H3 arrays",
        fontsize=25,
        fontweight="bold",
        x=0.07,
        y=0.985,
        ha="left",
    )
    figure.subplots_adjust(left=0.072, right=0.988, top=0.82, bottom=0.205, wspace=0.27)

    pdf = OUT / "Figure6_candidate_H1_H2_H3_localized_promoter_conservation.pdf"
    png = OUT / "Figure6_candidate_H1_H2_H3_localized_promoter_conservation.png"
    figure.savefig(pdf, bbox_inches="tight")
    figure.savefig(png, dpi=350, bbox_inches="tight", facecolor="white")
    plt.close(figure)
    return pdf, png


def main() -> None:
    base = load_base_module()
    base.configure_style()
    OUT.mkdir(parents=True, exist_ok=True)
    TABLES.mkdir(parents=True, exist_ok=True)

    classified, _ = base.classify_arrays()
    excluded_ties = classified[classified.dominant_fraction <= 0.5].copy()
    selected = classified[classified.dominant_fraction > 0.5].copy()

    base.GROUP_ORDER = GROUP_ORDER
    base.COLORS = COLORS
    base.LABELS = LABELS
    profile_metadata, profiles = base.load_array_profiles(selected)
    profile_summary = base.summarize_profiles(profile_metadata, profiles)
    interval_summary = base.summarize_profile_intervals(profile_summary)

    individual, summary = summarize_effects(base, selected, "strict_majority")
    pairwise, paired = pairwise_tests(base, selected, individual)

    high_purity = selected[selected.dominant_fraction >= 0.8].copy()
    _, sensitivity_summary = summarize_effects(
        base, high_purity, "dominant_fraction_ge_0.8"
    )

    selected.to_csv(TABLES / "H1_H2_H3_strict_majority_array_values.tsv.gz", sep="\t", index=False, compression="gzip")
    excluded_ties.to_csv(TABLES / "excluded_tied_arrays.tsv", sep="\t", index=False)
    profile_summary.to_csv(TABLES / "H1_H2_H3_local_divergence_profile.tsv", sep="\t", index=False)
    interval_summary.to_csv(TABLES / "H1_H2_H3_profile_interval_summary.tsv", sep="\t", index=False)
    individual.to_csv(TABLES / "H1_H2_H3_individual_paired_effects.tsv.gz", sep="\t", index=False, compression="gzip")
    summary.to_csv(TABLES / "H1_H2_H3_effect_summary.tsv", sep="\t", index=False)
    pairwise.to_csv(TABLES / "H1_H2_H3_pairwise_permutation_tests.tsv", sep="\t", index=False)
    paired.to_csv(TABLES / "H1_H2_H3_paired_individual_tests.tsv", sep="\t", index=False)
    sensitivity_summary.to_csv(TABLES / "H1_H2_H3_high_purity_effect_summary.tsv", sep="\t", index=False)

    pdf, png = make_figure(base, profile_summary, individual, summary, pairwise)
    print(pdf)
    print(png)
    print(summary.to_string(index=False))
    print(pairwise.to_string(index=False))
    print(paired.to_string(index=False))
    print(sensitivity_summary.to_string(index=False))


if __name__ == "__main__":
    main()
