#!/usr/bin/env python3
"""Compare localized L1PA10 promoter conservation in H2- and H3-dominant arrays."""

from __future__ import annotations

import json
import os
from pathlib import Path
import sys
import warnings

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon


ROOT = Path(__file__).resolve().parent
DATA_ROOT = Path(os.environ.get("ACRO1_WORKSPACE", ROOT / "inputs"))
ANALYSIS = Path(os.environ.get("ACRO1_PROMOTER_HUMAN_ONLY", DATA_ROOT / "promoter_human_only"))
FIG5_DATA = Path(os.environ.get("ACRO1_FIGURE5_DATA", DATA_ROOT / "figure5_data"))
OUT = ROOT / "figures"
TABLES = ROOT / "tables"

PP_SCRIPTS = ROOT.parent / "pipeline"
sys.path.insert(0, str(PP_SCRIPTS))
import pp_common as C  # noqa: E402
import pp_stats as S  # noqa: E402

CACHE = ANALYSIS / "cache"
ARRAY_VALUES = Path(os.environ.get("ACRO1_ARRAY_EFFECTS", FIG5_DATA / "Figure5B_array_values.tsv.gz"))
NTRIAGE_UNITS = Path(os.environ.get("ACRO1_HPRC_UNITS", DATA_ROOT / "hprc.clusters.tsv"))
NTRIAGE_ASSIGNMENTS = Path(os.environ.get("ACRO1_HPRC_ASSIGNMENTS", DATA_ROOT / "assignments.tsv.gz"))

PROMOTER_START = 74
PROMOTER_END = 123
FLANK_WIDTH = 50
PROFILE_START = PROMOTER_START - FLANK_WIDTH
PROFILE_END = PROMOTER_END + FLANK_WIDTH
GROUP_ORDER = ["H2", "H3"]
COLORS = {"H2": "#F28E2B", "H3": "#59A14F"}
LABELS = {"H2": "H2-dominant", "H3": "H3-dominant"}
N_PERMUTATIONS = 500_000


def configure_style() -> None:
    helvetica = Path("/System/Library/Fonts/Helvetica.ttc")
    if helvetica.exists():
        font_manager.fontManager.addfont(str(helvetica))
    mpl.rcParams.update(
        {
            "font.family": "Helvetica",
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


def classify_arrays() -> tuple[pd.DataFrame, pd.DataFrame]:
    effects = pd.read_csv(ARRAY_VALUES, sep="\t")
    units = pd.read_csv(
        NTRIAGE_UNITS,
        sep="\t",
        usecols=["acre_id", "array_id", "source_name"],
    )
    assignments = pd.read_csv(
        NTRIAGE_ASSIGNMENTS,
        sep="\t",
        usecols=["acre_id", "selected_cluster"],
    )
    merged = units.merge(assignments, on="acre_id", validate="one_to_one")
    if set(merged.selected_cluster.unique()) != {1, 2, 3}:
        raise ValueError("Expected selected clusters 1-3")

    counts = (
        merged.groupby(["array_id", "source_name", "selected_cluster"])
        .size()
        .unstack(fill_value=0)
    )
    for cluster in (1, 2, 3):
        if cluster not in counts.columns:
            counts[cluster] = 0
    counts = counts[[1, 2, 3]].rename(
        columns={1: "H1_count", 2: "H2_count", 3: "H3_count"}
    )
    counts["total_ACREs"] = counts[["H1_count", "H2_count", "H3_count"]].sum(axis=1)
    for cluster in ("H1", "H2", "H3"):
        counts[f"{cluster}_fraction"] = counts[f"{cluster}_count"] / counts.total_ACREs
    counts["dominant_cluster"] = counts[["H1_count", "H2_count", "H3_count"]].idxmax(axis=1).str[:2]
    counts["dominant_fraction"] = counts.apply(
        lambda row: row[f"{row.dominant_cluster}_fraction"], axis=1
    )
    counts = counts.reset_index().rename(columns={"array_id": "ntriage_array_id", "source_name": "array_id"})

    classified = effects.merge(counts, on="array_id", how="left", validate="one_to_one")
    if classified.dominant_cluster.isna().any():
        raise ValueError("Some promoter arrays did not map to an HPRC NTRiage array")
    selected = classified[classified.dominant_cluster.isin(GROUP_ORDER)].copy()
    if (selected.dominant_fraction <= 0.5).any():
        raise ValueError("H2/H3 dominant assignment did not have a strict majority")
    return classified, selected


def load_array_profiles(selected: pd.DataFrame) -> tuple[pd.DataFrame, np.ndarray]:
    selected_by_id = selected.set_index("array_id")
    records: list[dict[str, object]] = []
    profiles: list[np.ndarray] = []

    for cache_path in sorted(CACHE.glob("*.npz")):
        with np.load(cache_path, allow_pickle=False) as data:
            metadata = json.loads(str(data["result"]))
            array_id = str(metadata["array_id"])
            if array_id not in selected_by_id.index:
                continue
            rolling = data["rolling"].astype(float)
            callable_fraction = data["ncall"].astype(float) / int(metadata["n_anchors"])
            rolling[callable_fraction < 0.8] = np.nan

        row = selected_by_id.loc[array_id]
        records.append(
            {
                "array_id": array_id,
                "individual_id": str(metadata["individual_id"]),
                "dominant_cluster": str(row.dominant_cluster),
                "dominant_fraction": float(row.dominant_fraction),
            }
        )
        profiles.append(rolling)

    metadata = pd.DataFrame(records)
    profile_matrix = np.vstack(profiles)
    if len(metadata) != len(selected):
        missing = set(selected.array_id) - set(metadata.array_id)
        raise ValueError(f"Missing cached profiles for {len(missing)} selected arrays")
    return metadata, profile_matrix


def summarize_profiles(metadata: pd.DataFrame, profiles: np.ndarray) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    positions = np.arange(1, profiles.shape[1] + 1)

    for group in GROUP_ORDER:
        group_mask = metadata.dominant_cluster.eq(group).to_numpy()
        group_profiles = profiles[group_mask]
        group_ids = metadata.loc[group_mask, "individual_id"].to_numpy()
        individual_ids = sorted(set(group_ids))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            individual_profiles = np.vstack(
                [np.nanmean(group_profiles[group_ids == individual], axis=0) for individual in individual_ids]
            )
            mean_profile = np.nanmean(individual_profiles, axis=0)
            n_callable = np.sum(np.isfinite(individual_profiles), axis=0)
            standard_error = np.nanstd(individual_profiles, axis=0, ddof=1) / np.sqrt(n_callable)

        ci_low = np.maximum(0.0, mean_profile - 1.96 * standard_error)
        ci_high = mean_profile + 1.96 * standard_error
        for index, position in enumerate(positions):
            rows.append(
                {
                    "alignment_column": int(position),
                    "dominant_cluster": group,
                    "display_group": LABELS[group],
                    "arrays": int(group_mask.sum()),
                    "individuals": len(individual_ids),
                    "individuals_callable": int(n_callable[index]),
                    "individual_weighted_mean_rolling11_divergence": mean_profile[index],
                    "ci95_low": ci_low[index],
                    "ci95_high": ci_high[index],
                }
            )
    return pd.DataFrame(rows)


def summarize_profile_intervals(profile: pd.DataFrame) -> pd.DataFrame:
    intervals = {
        "left_50bp_flank": (PROFILE_START, PROMOTER_START - 1),
        "candidate_50bp_interval": (PROMOTER_START, PROMOTER_END),
        "right_50bp_flank": (PROMOTER_END + 1, PROFILE_END),
    }
    rows = []
    for group in GROUP_ORDER:
        subset = profile[profile.dominant_cluster.eq(group)]
        values = {}
        for label, (start, end) in intervals.items():
            interval = subset[subset.alignment_column.between(start, end)]
            values[label] = interval.individual_weighted_mean_rolling11_divergence.mean()
        flank_mean = np.mean([values["left_50bp_flank"], values["right_50bp_flank"]])
        rows.append(
            {
                "dominant_cluster": group,
                "display_group": LABELS[group],
                **values,
                "mean_of_two_flanks": flank_mean,
                "flank_minus_candidate": flank_mean - values["candidate_50bp_interval"],
            }
        )
    return pd.DataFrame(rows)


def individual_effects(frame: pd.DataFrame) -> pd.DataFrame:
    return (
        frame.groupby(["individual_id", "dominant_cluster"], as_index=False)
        .agg(
            mean_promoter_divergence=("promoter_divergence", "mean"),
            mean_control_median_divergence=("control_median_divergence", "mean"),
            mean_paired_effect=("effect", "mean"),
            arrays=("array_id", "size"),
            mean_cluster_fraction=("dominant_fraction", "mean"),
        )
    )


def effect_summary(frame: pd.DataFrame, label: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.Series]:
    individual = individual_effects(frame)
    rows = []
    for group in GROUP_ORDER:
        subset_arrays = frame[frame.dominant_cluster.eq(group)]
        values = individual.loc[
            individual.dominant_cluster.eq(group), "mean_paired_effect"
        ].to_numpy(float)
        low, high = S.bootstrap_mean_ci(values, f"H2-H3|{label}|{group}")
        rows.append(
            {
                "analysis": label,
                "dominant_cluster": group,
                "display_group": LABELS[group],
                "arrays": len(subset_arrays),
                "individuals": len(values),
                "median_dominant_fraction": subset_arrays.dominant_fraction.median(),
                "mean_effect": values.mean(),
                "median_effect": np.median(values),
                "fraction_positive": (values > 0).mean(),
                "bootstrap_ci95_low": low,
                "bootstrap_ci95_high": high,
            }
        )
    summary = pd.DataFrame(rows)

    wide = individual.pivot(
        index="individual_id", columns="dominant_cluster", values="mean_paired_effect"
    ).dropna(subset=GROUP_ORDER)
    differences = wide.H2 - wide.H3
    paired = pd.Series(
        {
            "analysis": label,
            "paired_individuals": len(wide),
            "mean_H2_minus_H3": differences.mean(),
            "median_H2_minus_H3": differences.median(),
            "paired_wilcoxon_p": wilcoxon(differences).pvalue if len(wide) else np.nan,
        }
    )
    return individual, summary, paired


def within_individual_permutation(frame: pd.DataFrame) -> pd.Series:
    group = frame.dominant_cluster.eq("H2").to_numpy()
    strata = pd.factorize(frame.individual_id)[0]
    result = S.perm_diff_means(
        frame.effect.to_numpy(float),
        group,
        N_PERMUTATIONS,
        "H2-H3-dominant-array-promoter-effect",
        strata=strata,
    )
    return pd.Series(
        {
            "comparison": "H2_minus_H3",
            "arrays_H2": int(group.sum()),
            "arrays_H3": int((~group).sum()),
            "individuals": frame.individual_id.nunique(),
            "difference_in_array_means": result["diff"],
            "p_two_sided": result["p_two_sided"],
            "p_H2_greater": result["p_group1_greater"],
            "permutations": result["draws"],
        }
    )


def make_figure(
    profile: pd.DataFrame,
    individual: pd.DataFrame,
    summary: pd.DataFrame,
    permutation: pd.Series,
) -> tuple[Path, Path]:
    figure, (axis_profile, axis_effect) = plt.subplots(
        1,
        2,
        figsize=(14.4, 6.6),
        gridspec_kw={"width_ratios": [1.42, 1.0]},
    )

    local = profile[profile.alignment_column.between(PROFILE_START, PROFILE_END)].copy()
    axis_profile.axvspan(
        PROMOTER_START - 0.5,
        PROMOTER_END + 0.5,
        color="#F4D03F",
        alpha=0.25,
        zorder=0,
    )
    axis_profile.axvline(PROMOTER_START - 0.5, color="#B7950B", linewidth=1.0, linestyle=":")
    axis_profile.axvline(PROMOTER_END + 0.5, color="#B7950B", linewidth=1.0, linestyle=":")

    for group in GROUP_ORDER:
        subset = local[local.dominant_cluster.eq(group)]
        x = subset.alignment_column.to_numpy()
        mean = 100 * subset.individual_weighted_mean_rolling11_divergence.to_numpy()
        low = 100 * subset.ci95_low.to_numpy()
        high = 100 * subset.ci95_high.to_numpy()
        axis_profile.fill_between(x, low, high, color=COLORS[group], alpha=0.13, linewidth=0)
        axis_profile.plot(
            x,
            mean,
            color=COLORS[group],
            linewidth=2.4,
            label=LABELS[group],
        )

    upper = 1.12 * 100 * local.ci95_high.max()
    axis_profile.set_ylim(0, upper)
    axis_profile.set_xlim(PROFILE_START, PROFILE_END)
    axis_profile.text(
        (PROMOTER_START + PROMOTER_END) / 2,
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
            widths=0.78,
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

        jitter = rng.normal(0, 0.075, size=len(values))
        axis_effect.scatter(
            position + jitter,
            values,
            s=13,
            alpha=0.28,
            color=COLORS[group],
            edgecolors="none",
            rasterized=True,
            zorder=2,
        )
        axis_effect.boxplot(
            [values],
            positions=[position],
            widths=0.17,
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
            markersize=7.5,
            markerfacecolor="#F2C14E",
            markeredgecolor="#171717",
            markeredgewidth=1.1,
            ecolor="#171717",
            elinewidth=1.4,
            capsize=4,
            zorder=4,
        )
        axis_effect.text(
            position + 0.16,
            mean,
            f"mean {mean:.2f} pp",
            ha="left",
            va="center",
            fontsize=11.5,
            fontweight="bold",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78, "pad": 1.0},
        )

    axis_effect.axhline(0, color="#222222", linewidth=1.1, linestyle="--", zorder=1)
    y_min = 100 * individual.mean_paired_effect.min()
    y_max = 100 * individual.mean_paired_effect.max()
    bracket_y = y_max + 0.33
    bracket_height = 0.14
    axis_effect.plot(
        [0, 0, 1, 1],
        [bracket_y, bracket_y + bracket_height, bracket_y + bracket_height, bracket_y],
        color="#222222",
        linewidth=1.1,
    )
    delta = 100 * permutation.difference_in_array_means
    axis_effect.text(
        0.5,
        bracket_y + bracket_height + 0.07,
        f"Difference in means = {delta:.2f} pp; P = {permutation.p_two_sided:.1e}",
        ha="center",
        va="bottom",
        fontsize=11.5,
        fontweight="bold",
    )
    axis_effect.set_ylim(y_min - 0.34, bracket_y + bracket_height + 0.53)
    axis_effect.set_xlim(-0.55, 1.55)
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
        "Localized L1PA10 promoter conservation in H2 and H3 arrays",
        fontsize=25,
        fontweight="bold",
        x=0.07,
        y=0.985,
        ha="left",
    )
    figure.subplots_adjust(left=0.075, right=0.985, top=0.83, bottom=0.20, wspace=0.28)

    pdf = OUT / "Figure6_candidate_H2_H3_localized_promoter_conservation.pdf"
    png = OUT / "Figure6_candidate_H2_H3_localized_promoter_conservation.png"
    figure.savefig(pdf, bbox_inches="tight")
    figure.savefig(png, dpi=350, bbox_inches="tight", facecolor="white")
    plt.close(figure)
    return pdf, png


def main() -> None:
    configure_style()
    OUT.mkdir(parents=True, exist_ok=True)
    TABLES.mkdir(parents=True, exist_ok=True)

    classified, selected = classify_arrays()
    profile_metadata, profiles = load_array_profiles(selected)
    profile_summary = summarize_profiles(profile_metadata, profiles)
    interval_summary = summarize_profile_intervals(profile_summary)

    individual, summary, paired = effect_summary(selected, "dominant_cluster")
    permutation = within_individual_permutation(selected)
    high_purity = selected[selected.dominant_fraction >= 0.8].copy()
    sensitivity_individual, sensitivity_summary, sensitivity_paired = effect_summary(
        high_purity, "dominant_fraction_ge_0.8"
    )

    classified.to_csv(TABLES / "all_estimable_array_cluster_classification.tsv.gz", sep="\t", index=False, compression="gzip")
    selected.to_csv(TABLES / "H2_H3_dominant_array_values.tsv.gz", sep="\t", index=False, compression="gzip")
    profile_summary.to_csv(TABLES / "H2_H3_local_divergence_profile.tsv", sep="\t", index=False)
    interval_summary.to_csv(TABLES / "H2_H3_profile_interval_summary.tsv", sep="\t", index=False)
    individual.to_csv(TABLES / "H2_H3_individual_paired_effects.tsv.gz", sep="\t", index=False, compression="gzip")
    summary.to_csv(TABLES / "H2_H3_effect_summary.tsv", sep="\t", index=False)
    pd.DataFrame([paired]).to_csv(TABLES / "H2_H3_paired_individual_test.tsv", sep="\t", index=False)
    pd.DataFrame([permutation]).to_csv(TABLES / "H2_H3_within_individual_permutation.tsv", sep="\t", index=False)
    sensitivity_individual.to_csv(TABLES / "H2_H3_high_purity_individual_effects.tsv.gz", sep="\t", index=False, compression="gzip")
    sensitivity_summary.to_csv(TABLES / "H2_H3_high_purity_effect_summary.tsv", sep="\t", index=False)
    pd.DataFrame([sensitivity_paired]).to_csv(TABLES / "H2_H3_high_purity_paired_test.tsv", sep="\t", index=False)

    pdf, png = make_figure(profile_summary, individual, summary, permutation)
    print(pdf)
    print(png)
    print(summary.to_string(index=False))
    print(pd.DataFrame([permutation]).to_string(index=False))
    print(pd.DataFrame([paired]).to_string(index=False))
    print(sensitivity_summary.to_string(index=False))


if __name__ == "__main__":
    main()
