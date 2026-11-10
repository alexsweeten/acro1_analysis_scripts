#!/usr/bin/env python3
"""Split the weakest HG002 acrocentric K=3 cluster and compare with global K=4."""

from __future__ import annotations

import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent
(ROOT / ".mplconfig").mkdir(parents=True, exist_ok=True)
(ROOT / ".cache" / "fontconfig").mkdir(parents=True, exist_ok=True)
os.environ["LOKY_MAX_CPU_COUNT"] = str(os.cpu_count() or 1)
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".mplconfig"))
os.environ.setdefault("XDG_CACHE_HOME", str(ROOT / ".cache"))

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "font.size": 12,
        "axes.titlesize": 14,
        "axes.labelsize": 12,
        "axes.unicode_minus": False,
        "legend.fontsize": 11,
        "legend.title_fontsize": 12,
        "pdf.use14corefonts": True,
        "ps.useafm": True,
    }
)
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, silhouette_score


SOURCE = ROOT.parent / "hg002_acrocentric_k3"
FEATURES = SOURCE / "native" / "hg002_acrocentric_k3.features.tsv"
K3_ASSIGNMENTS = SOURCE / "results" / "unit_assignments.tsv"
K3_SUPPORT = SOURCE / "results" / "fixed_k3_cluster_support.tsv"
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
PROVENANCE = ROOT / "provenance"

CHROMOSOMES = ("chr13", "chr14", "chr15", "chr21", "chr22")
FINAL_LABELS = ("H1", "H2", "H3", "H4")
SOURCE_TO_HUMAN = {"C1": "H1", "C2": "H2", "C3": "H3"}
COLORS = {
    "H1": "#4C78A8",
    "H2": "#F58518",
    "H3": "#54A24B",
    "H4": "#E45756",
}
TARGET_INTERVALS = (
    {
        "target_interval": "chr15_PATERNAL:2254000-2550000",
        "original_chrom": "chr15_PATERNAL",
        "start": 2_254_000,
        "end": 2_550_000,
        "color": "#6A3D9A",
    },
    {
        "target_interval": "chr14_PATERNAL:110000-376000",
        "original_chrom": "chr14_PATERNAL",
        "start": 110_000,
        "end": 376_000,
        "color": "#111111",
    },
)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def scale_features(features: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    raw = features.iloc[:, 3:].to_numpy(dtype=np.float64, copy=True)
    means = raw.mean(axis=0, dtype=np.float64)
    scales = raw.std(axis=0, ddof=1, dtype=np.float64)
    usable = np.isfinite(scales) & (scales > 1e-12)
    x = (raw[:, usable] - means[usable]) / scales[usable]
    norms = np.linalg.norm(x, axis=1)
    nonzero = norms > 0
    x[nonzero] /= norms[nonzero, None]
    x[~nonzero] = 0.0
    return x, usable


def stable_integer_labels(labels: np.ndarray) -> np.ndarray:
    first = [
        (int(np.flatnonzero(labels == label)[0]), int(label))
        for label in np.unique(labels)
    ]
    remap = {old: new for new, (_, old) in enumerate(sorted(first))}
    return np.fromiter((remap[int(label)] for label in labels), dtype=np.int32)


def align_integer_labels(reference: np.ndarray, labels: np.ndarray, k: int) -> np.ndarray:
    contingency = np.zeros((k, k), dtype=np.int64)
    for ref, label in zip(reference, labels, strict=True):
        contingency[int(ref), int(label)] += 1
    rows, columns = linear_sum_assignment(-contingency)
    remap = {int(column): int(row) for row, column in zip(rows, columns, strict=True)}
    return np.fromiter((remap[int(label)] for label in labels), dtype=np.int32)


def fit_multistart(
    x: np.ndarray, k: int, seeds: range
) -> tuple[np.ndarray, np.ndarray, pd.DataFrame, dict]:
    fitted = []
    for seed in seeds:
        model = KMeans(
            n_clusters=k,
            init="k-means++",
            n_init=1,
            max_iter=300,
            tol=1e-4,
            algorithm="lloyd",
            random_state=seed,
        ).fit(x)
        fitted.append((seed, model))
    best_seed, best_model = min(fitted, key=lambda value: value[1].inertia_)
    selected = stable_integer_labels(best_model.labels_)
    aligned = np.column_stack(
        [align_integer_labels(selected, model.labels_, k) for _, model in fitted]
    )
    agreement = (aligned == selected[:, None]).mean(axis=1)
    runs = []
    for (seed, model), labels in zip(fitted, aligned.T, strict=True):
        runs.append(
            {
                "seed": seed,
                "k": k,
                "inertia": model.inertia_,
                "iterations": model.n_iter_,
                "selected_best_inertia": seed == best_seed,
                "ari_to_selected": adjusted_rand_score(selected, labels),
            }
        )
    pairwise_ari = [
        adjusted_rand_score(aligned[:, left], aligned[:, right])
        for left in range(aligned.shape[1])
        for right in range(left + 1, aligned.shape[1])
    ]
    metadata = {
        "k": k,
        "starts": len(fitted),
        "best_seed": best_seed,
        "best_inertia": float(best_model.inertia_),
        "pairwise_ari_mean": float(np.mean(pairwise_ari)),
        "pairwise_ari_min": float(np.min(pairwise_ari)),
    }
    return selected, agreement, pd.DataFrame(runs), metadata


def align_named_labels(
    reference: np.ndarray, labels: np.ndarray, names: tuple[str, ...]
) -> np.ndarray:
    contingency = np.zeros((len(names), len(names)), dtype=np.int64)
    for ref, label in zip(reference, labels, strict=True):
        contingency[names.index(str(ref)), int(label)] += 1
    rows, columns = linear_sum_assignment(-contingency)
    remap = {
        int(column): names[int(row)]
        for row, column in zip(rows, columns, strict=True)
    }
    return np.asarray([remap[int(label)] for label in labels], dtype=object)


def summarize_clusters(assignments: pd.DataFrame, label_column: str) -> pd.DataFrame:
    rows = []
    for cluster in FINAL_LABELS:
        group = assignments.loc[assignments[label_column].eq(cluster)]
        s1 = int(group["subfamily_1"].sum())
        s2 = int(group["subfamily_2"].sum())
        s3 = int(group["subfamily_3"].sum())
        classified = s1 + s2 + s3
        rows.append(
            {
                "cluster": cluster,
                "units": len(group),
                "arrays": group["array_id"].nunique(),
                "chromosome_haplotypes": group["original_chrom"].nunique(),
                "mean_unit_length_bp": group["length"].mean(),
                "median_unit_length_bp": group["length"].median(),
                "min_unit_length_bp": group["length"].min(),
                "max_unit_length_bp": group["length"].max(),
                "mean_classified_acro1_copies": group["classified_total"].mean(),
                "subfamily_1": s1,
                "subfamily_2": s2,
                "subfamily_3": s3,
                "classified_total": classified,
                "fraction_subfamily_1": s1 / classified if classified else np.nan,
                "fraction_subfamily_2": s2 / classified if classified else np.nan,
                "fraction_subfamily_3": s3 / classified if classified else np.nan,
            }
        )
    return pd.DataFrame(rows)


def summarize_by_chromosome_haplotype(
    assignments: pd.DataFrame, label_column: str
) -> pd.DataFrame:
    output = (
        assignments.groupby(
            ["chromosome", "haplotype", label_column], observed=True
        )
        .agg(
            units=("acre_id", "size"),
            arrays=("array_id", "nunique"),
            mean_unit_length_bp=("length", "mean"),
        )
        .reset_index()
        .rename(columns={label_column: "cluster"})
    )
    totals = output.groupby(["chromosome", "haplotype"])["units"].transform("sum")
    output["fraction_within_chromosome_haplotype"] = output["units"] / totals
    return output


def support_table(
    labels: np.ndarray,
    agreement: np.ndarray,
    names: tuple[str, ...],
    support_source: str,
) -> pd.DataFrame:
    rows = []
    for index, name in enumerate(names):
        values = agreement[labels == index]
        rows.append(
            {
                "cluster": name,
                "units": len(values),
                "mean_agreement": values.mean(),
                "agreement_q05": np.quantile(values, 0.05),
                "fraction_agreement_below_0.75": (values < 0.75).mean(),
                "support_source": support_source,
            }
        )
    return pd.DataFrame(rows)


def annotate_target_intervals(assignments: pd.DataFrame) -> pd.DataFrame:
    assignments = assignments.copy()
    assignments["target_interval"] = ""
    for target in TARGET_INTERVALS:
        overlap = (
            assignments["original_chrom"].eq(target["original_chrom"])
            & assignments["genome_start"].lt(target["end"])
            & assignments["genome_end"].gt(target["start"])
        )
        assignments.loc[overlap, "target_interval"] = target["target_interval"]
    assignments["target_interval_overlap"] = assignments["target_interval"].ne("")
    return assignments


def draw_fig3_style_pca(
    ax: plt.Axes,
    assignments: pd.DataFrame,
    explained: tuple[float, float],
    *,
    standalone: bool,
) -> None:
    point_size = 34 if standalone else 25
    for cluster in FINAL_LABELS:
        part = assignments.loc[assignments["constrained_k4_cluster"].eq(cluster)]
        ax.scatter(
            part["pca1"],
            part["pca2"],
            s=point_size,
            marker="o",
            color=COLORS[cluster],
            alpha=0.67,
            edgecolors="none",
            rasterized=True,
            zorder=2,
        )

    centroid_size = 230 if standalone else 115
    label_size = 18 if standalone else 11
    y_offset = 0.035 if standalone else 0.026
    for cluster in FINAL_LABELS:
        part = assignments.loc[assignments["constrained_k4_cluster"].eq(cluster)]
        x_center = part["pca1"].mean()
        y_center = part["pca2"].mean()
        ax.scatter(
            [x_center],
            [y_center],
            s=centroid_size,
            marker="D",
            facecolor=COLORS[cluster],
            edgecolor="#111111",
            linewidth=1.4,
            zorder=5,
        )
        ax.text(
            x_center,
            y_center + y_offset,
            cluster,
            ha="center",
            va="bottom",
            fontsize=label_size,
            fontweight="bold",
            color="#111111",
            zorder=6,
        )

    target_handles = []
    for target in TARGET_INTERVALS:
        part = assignments.loc[assignments["target_interval"].eq(target["target_interval"])]
        ax.scatter(
            part["pca1"],
            part["pca2"],
            s=point_size * 2.2,
            marker="o",
            facecolors="none",
            edgecolors=target["color"],
            linewidths=1.35 if standalone else 1.0,
            zorder=4,
        )
        target_handles.append(
            plt.Line2D(
                [0], [0], marker="o", linestyle="none", markerfacecolor="none",
                markeredgecolor=target["color"], markeredgewidth=1.5,
                markersize=8, label=f"{target['target_interval']} (n={len(part)})"
            )
        )

    ax.set_axisbelow(True)
    ax.grid(color="#DFE3E8", linewidth=0.85)
    for spine in ax.spines.values():
        spine.set_color("#333333")
        spine.set_linewidth(0.9)
    ax.tick_params(axis="both", labelsize=18 if standalone else 10, width=1.0)
    ax.set_xlabel(
        f"PC1 ({explained[0] * 100:.1f}% variance)",
        fontsize=24 if standalone else 12,
        fontweight="bold",
    )
    ax.set_ylabel(
        f"PC2 ({explained[1] * 100:.1f}% variance)",
        fontsize=24 if standalone else 12,
        fontweight="bold",
    )

    if standalone:
        ax.set_title(
            "HG002 acrocentric ACRO-unit PCA",
            loc="left",
            fontsize=32,
            fontweight="bold",
            pad=16,
        )
        cluster_handles = [
            plt.Line2D(
                [0], [0], marker="o", linestyle="none",
                markerfacecolor=COLORS[cluster], markeredgecolor="none",
                markersize=13,
                label=(
                    f"{cluster}  (n="
                    f"{assignments['constrained_k4_cluster'].eq(cluster).sum()})"
                ),
            )
            for cluster in FINAL_LABELS
        ]
        cluster_legend = ax.legend(
            handles=cluster_handles,
            title="NTRiage cluster",
            frameon=False,
            loc="upper left",
            bbox_to_anchor=(1.02, 0.78),
            fontsize=18,
            title_fontsize=22,
            labelspacing=0.65,
            borderaxespad=0,
        )
        cluster_legend.get_title().set_fontweight("bold")
        ax.add_artist(cluster_legend)
        target_legend = ax.legend(
            handles=target_handles,
            title="Circled intervals",
            frameon=False,
            loc="upper left",
            bbox_to_anchor=(1.02, 0.34),
            fontsize=11.5,
            title_fontsize=14,
            labelspacing=0.8,
            borderaxespad=0,
        )
        target_legend.get_title().set_fontweight("bold")
    else:
        ax.set_title(
            "A  HG002 acrocentric ACRO-unit PCA",
            loc="left",
            fontsize=13,
            fontweight="bold",
        )
        ax.legend(
            handles=target_handles,
            title="Circled intervals",
            frameon=False,
            loc="lower right",
            fontsize=7.5,
            title_fontsize=8.5,
        )


def make_standalone_pca(
    assignments: pd.DataFrame, explained: tuple[float, float]
) -> None:
    fig, ax = plt.subplots(figsize=(12.2, 9.2))
    draw_fig3_style_pca(ax, assignments, explained, standalone=True)
    fig.subplots_adjust(left=0.12, right=0.70, bottom=0.14, top=0.88)
    fig.savefig(
        FIGURES / "hg002_acrocentric_k4_pca_fig3_style.png",
        dpi=300,
        bbox_inches="tight",
        facecolor="white",
    )
    fig.savefig(
        FIGURES / "hg002_acrocentric_k4_pca_fig3_style.pdf",
        bbox_inches="tight",
        facecolor="white",
    )
    plt.close(fig)


def make_figure(
    assignments: pd.DataFrame,
    by_chrom_hap: pd.DataFrame,
    cluster_summary: pd.DataFrame,
    explained: tuple[float, float],
) -> None:
    fig, axes = plt.subplots(
        1, 3,
        figsize=(17.0, 6.1),
        gridspec_kw={"width_ratios": [1.25, 1.05, 0.9]},
    )

    draw_fig3_style_pca(axes[0], assignments, explained, standalone=False)

    ax = axes[1]
    row_order = [
        (chromosome, haplotype)
        for chromosome in CHROMOSOMES
        for haplotype in ("MATERNAL", "PATERNAL")
    ]
    index = pd.MultiIndex.from_tuples(row_order, names=["chromosome", "haplotype"])
    pivot = (
        by_chrom_hap.pivot(
            index=["chromosome", "haplotype"],
            columns="cluster",
            values="fraction_within_chromosome_haplotype",
        )
        .reindex(index)
        .fillna(0)
    )
    bottom = np.zeros(len(pivot))
    for cluster in FINAL_LABELS:
        values = pivot.get(cluster, pd.Series(0, index=pivot.index)).to_numpy()
        bars = ax.barh(
            np.arange(len(pivot)), values, left=bottom, color=COLORS[cluster],
            edgecolor="white", linewidth=0.5, label=cluster
        )
        for bar, value, left in zip(bars, values, bottom, strict=True):
            if value >= 0.12:
                ax.text(
                    left + value / 2,
                    bar.get_y() + bar.get_height() / 2,
                    f"{value * 100:.0f}%",
                    ha="center",
                    va="center",
                    color="white",
                    fontsize=7.5,
                    fontweight="bold",
                )
        bottom += values
    counts = (
        assignments.groupby(["chromosome", "haplotype"])
        .size()
        .reindex(index)
        .astype(int)
    )
    labels = [
        f"{chromosome} {'M' if haplotype == 'MATERNAL' else 'P'}  "
        f"(n={counts.loc[(chromosome, haplotype)]})"
        for chromosome, haplotype in row_order
    ]
    ax.set_yticks(np.arange(len(index)))
    ax.set_yticklabels(labels, fontsize=8.3)
    ax.set_xlim(0, 1)
    ax.set_xlabel("Fraction of complete ACRE units", labelpad=5)
    ax.set_title("B  Chromosome/haplotype composition", loc="left", fontweight="bold")
    ax.invert_yaxis()
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.legend(frameon=False, ncol=4, loc="lower center", bbox_to_anchor=(0.5, -0.22))

    ax = axes[2]
    composition = cluster_summary.set_index("cluster")[[
        "fraction_subfamily_1", "fraction_subfamily_2", "fraction_subfamily_3"
    ]]
    composition.columns = ["S1", "S2", "S3"]
    composition = composition.reindex(FINAL_LABELS)
    bottom = np.zeros(len(composition))
    allele_colors = {"S1": "#4C78A8", "S2": "#F58518", "S3": "#B279A2"}
    for allele in ("S1", "S2", "S3"):
        values = composition[allele].to_numpy()
        bars = ax.bar(
            np.arange(len(composition)), values, bottom=bottom,
            color=allele_colors[allele], edgecolor="white", linewidth=0.6,
            label=allele,
        )
        for bar, value, lower in zip(bars, values, bottom, strict=True):
            if value >= 0.10:
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    lower + value / 2,
                    f"{value * 100:.0f}%",
                    ha="center",
                    va="center",
                    color="white",
                    fontsize=8,
                    fontweight="bold",
                )
        bottom += values
    ax.set_xticks(np.arange(len(composition)), FINAL_LABELS)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Fraction of classified ACRO1 copies")
    ax.set_title("C  ACRO1 composition", loc="left", fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, ncol=3, loc="lower center", bbox_to_anchor=(0.5, -0.16))

    fig.suptitle(
        "HG002 acrocentric K=4 sensitivity analysis: split the weakest K=3 cluster",
        x=0.01,
        ha="left",
        fontsize=15,
        fontweight="bold",
    )
    fig.tight_layout(w_pad=2.3, rect=(0, 0, 1, 0.95))
    fig.savefig(
        FIGURES / "hg002_acrocentric_k4_weakest_split.png",
        dpi=300,
        bbox_inches="tight",
        facecolor="white",
    )
    fig.savefig(
        FIGURES / "hg002_acrocentric_k4_weakest_split.pdf",
        bbox_inches="tight",
        facecolor="white",
    )
    plt.close(fig)


def write_results_markdown(
    cluster_summary: pd.DataFrame,
    by_chrom_hap: pd.DataFrame,
    constrained_support: pd.DataFrame,
    global_support: pd.DataFrame,
    crosswalk: pd.DataFrame,
    metadata: dict,
) -> None:
    cluster_lines = []
    for row in cluster_summary.itertuples(index=False):
        cluster_lines.append(
            f"| {row.cluster} | {row.units} | {row.mean_unit_length_bp:,.1f} | "
            f"{row.median_unit_length_bp:,.0f} | {row.fraction_subfamily_1 * 100:.1f}% | "
            f"{row.fraction_subfamily_2 * 100:.1f}% | "
            f"{row.fraction_subfamily_3 * 100:.1f}% |"
        )

    composition = by_chrom_hap.pivot_table(
        index=["chromosome", "haplotype"],
        columns="cluster",
        values="units",
        aggfunc="sum",
        fill_value=0,
    )
    composition_lines = []
    for chromosome in CHROMOSOMES:
        for haplotype in ("MATERNAL", "PATERNAL"):
            row = composition.loc[(chromosome, haplotype)]
            composition_lines.append(
                f"| {chromosome} | {haplotype.title()} | "
                f"{int(row.get('H1', 0))} | {int(row.get('H2', 0))} | "
                f"{int(row.get('H3', 0))} | {int(row.get('H4', 0))} |"
            )

    support_lines = []
    for row in constrained_support.itertuples(index=False):
        support_lines.append(
            f"| {row.cluster} | {row.units} | {row.mean_agreement:.3f} | "
            f"{row.agreement_q05:.3f} | {row.support_source} |"
        )

    global_support_lines = []
    for row in global_support.itertuples(index=False):
        global_support_lines.append(
            f"| {row.cluster} | {row.units} | {row.mean_agreement:.3f} | "
            f"{row.agreement_q05:.3f} |"
        )

    matrix = crosswalk.pivot(
        index="constrained_k4_cluster", columns="global_k4_cluster", values="units"
    ).reindex(index=FINAL_LABELS, columns=FINAL_LABELS, fill_value=0).fillna(0)
    crosswalk_lines = []
    for cluster in FINAL_LABELS:
        row = matrix.loc[cluster]
        crosswalk_lines.append(
            f"| {cluster} | " + " | ".join(str(int(row[name])) for name in FINAL_LABELS) + " |"
        )

    target_lines = []
    for target in metadata["target_intervals"]:
        target_lines.append(
            f"| {target['target_interval']} | {target['acres']} | {target['arrays']} | "
            f"{target['H1']} | {target['H2']} | {target['H3']} | {target['H4']} |"
        )

    text = f"""# HG002 acrocentric K=4 sensitivity analysis

This analysis starts from the 487 complete ACRE units in the acrocentric-only
HG002 K=3 run. The weakest K=3 cluster was {metadata['weakest_k3_cluster']}
(mean assignment agreement {metadata['weakest_k3_mean_agreement']:.3f}). H1 and
H2 were held fixed, and the {metadata['weakest_k3_units']} H3 units were split
with 50 K-means++ starts at K=2. The stable S3-rich daughter retains the H3
label, while the smaller S1-rich daughter is labeled H4.

The conditional H3 split was {metadata['daughter_sizes'][0]}/{metadata['daughter_sizes'][1]} units.
It reduced within-H3 inertia by {metadata['split_inertia_reduction_fraction'] * 100:.1f}%
and had silhouette {metadata['split_silhouette']:.3f}. H3 is the stable,
S3-rich daughter; H4 is smaller, more S1-rich, and retains most of the
multi-start ambiguity.

## Constrained split summary

| Cluster | Units | Mean ACRE length (bp) | Median (bp) | S1 | S2 | S3 |
|---|---:|---:|---:|---:|---:|---:|
{chr(10).join(cluster_lines)}

| Chromosome | Haplotype | H1 | H2 | H3 | H4 |
|---|---|---:|---:|---:|---:|
{chr(10).join(composition_lines)}

## Circled PCA intervals

| Interval | ACREs | Arrays | H1 | H2 | H3 | H4 |
|---|---:|---:|---:|---:|---:|---:|
{chr(10).join(target_lines)}

## Multi-start support

The H1/H2 values below are inherited from the original K=3 fit because those
assignments were held fixed. H3/H4 support is conditional on membership in
the original H3 parent.

| Cluster | Units | Mean agreement | 5th percentile | Support source |
|---|---:|---:|---:|---|
{chr(10).join(support_lines)}

## Unconstrained global K=4 check

An independent 50-start global K=4 fit agreed closely with the constrained
split (ARI {metadata['global_vs_constrained_ari']:.3f};
{metadata['global_vs_constrained_exact_fraction'] * 100:.1f}% exact assignment).
It preserved every original H1 and H2 assignment. Among the 125 original H3
units, {metadata['global_k4_changed_from_constrained']} differed from the
constrained daughter assignment: four moved to global H1 and seven moved
between H3 and H4.

| Constrained \\ Global | H1 | H2 | H3 | H4 |
|---|---:|---:|---:|---:|
{chr(10).join(crosswalk_lines)}

| Global cluster | Units | Mean agreement | 5th percentile |
|---|---:|---:|---:|
{chr(10).join(global_support_lines)}

The split therefore isolates a reproducible S3-rich H3, but it does not
eliminate the weak boundary: H4 remains the least stable component.
"""
    (ROOT / "RESULTS.md").write_text(text)


def main() -> None:
    for path in (FEATURES, K3_ASSIGNMENTS, K3_SUPPORT):
        if not path.is_file():
            raise FileNotFoundError(path)
    for directory in (RESULTS, FIGURES, PROVENANCE, ROOT / ".mplconfig"):
        directory.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".mplconfig"))

    features = pd.read_csv(FEATURES, sep="\t")
    assignments = pd.read_csv(K3_ASSIGNMENTS, sep="\t")
    k3_support = pd.read_csv(K3_SUPPORT, sep="\t")
    if not features["acre_id"].is_unique or not assignments["acre_id"].is_unique:
        raise RuntimeError("ACRE identifiers are not unique")
    assignments = (
        features[["acre_id"]]
        .merge(assignments, on="acre_id", validate="one_to_one")
    )
    if len(assignments) != len(features):
        raise RuntimeError("feature and K=3 assignment rows do not match")
    x, usable = scale_features(features)

    weakest = k3_support.sort_values("mean_agreement").iloc[0]
    weakest_source_cluster = str(weakest["cluster"])
    weakest_cluster = SOURCE_TO_HUMAN[weakest_source_cluster]
    parent_mask = assignments["cluster"].eq(weakest_source_cluster).to_numpy()
    parent_x = x[parent_mask]
    daughter_labels, daughter_agreement, daughter_runs, daughter_fit = fit_multistart(
        parent_x, 2, range(2001, 2051)
    )
    # The first-occurring daughter is the smaller, S1-rich residual (H4), while
    # the second is the stable, S3-rich daughter that retains the H3 label.
    daughter_names = ("H4", "H3")
    constrained = assignments["cluster"].map(SOURCE_TO_HUMAN).to_numpy(dtype=object)
    if pd.isna(constrained).any():
        raise RuntimeError("one or more source K=3 labels lack a human H-label mapping")
    constrained[parent_mask] = np.asarray(
        [daughter_names[value] for value in daughter_labels], dtype=object
    )
    assignments["k3_cluster"] = assignments["cluster"].map(SOURCE_TO_HUMAN)
    assignments["constrained_k4_cluster"] = constrained
    assignments["conditional_split_agreement"] = np.nan
    assignments.loc[parent_mask, "conditional_split_agreement"] = daughter_agreement
    assignments = annotate_target_intervals(assignments)

    one_cluster_inertia = float(np.sum((parent_x - parent_x.mean(axis=0)) ** 2))
    split_inertia_reduction = 1.0 - daughter_fit["best_inertia"] / one_cluster_inertia
    split_silhouette = float(silhouette_score(parent_x, daughter_labels))

    daughter_support = support_table(
        daughter_labels,
        daughter_agreement,
        daughter_names,
        "conditional K=2 split of original H3",
    )
    inherited_support = k3_support.loc[k3_support["cluster"].isin(("C1", "C2"))].copy()
    inherited_support["cluster"] = inherited_support["cluster"].map(SOURCE_TO_HUMAN)
    inherited_support["support_source"] = "inherited from original K=3 fit"
    constrained_support = pd.concat(
        [
            inherited_support[
                [
                    "cluster", "units", "mean_agreement", "agreement_q05",
                    "fraction_agreement_below_0.75", "support_source",
                ]
            ],
            daughter_support,
        ],
        ignore_index=True,
    )
    constrained_support["cluster"] = pd.Categorical(
        constrained_support["cluster"], categories=FINAL_LABELS, ordered=True
    )
    constrained_support = constrained_support.sort_values("cluster").reset_index(drop=True)
    constrained_support["cluster"] = constrained_support["cluster"].astype(str)

    constrained_summary = summarize_clusters(assignments, "constrained_k4_cluster")
    constrained_by_chrom_hap = summarize_by_chromosome_haplotype(
        assignments, "constrained_k4_cluster"
    )

    global_labels_int, _, global_runs_initial, global_fit = fit_multistart(
        x, 4, range(3001, 3051)
    )
    global_labels = align_named_labels(constrained, global_labels_int, FINAL_LABELS)

    # Refit the agreement calculation against the named best global solution.
    global_models = []
    for seed in range(3001, 3051):
        model = KMeans(
            n_clusters=4,
            init="k-means++",
            n_init=1,
            max_iter=300,
            tol=1e-4,
            algorithm="lloyd",
            random_state=seed,
        ).fit(x)
        global_models.append((seed, model))
    global_aligned = np.column_stack(
        [align_named_labels(global_labels, model.labels_, FINAL_LABELS) for _, model in global_models]
    )
    global_agreement = (global_aligned == global_labels[:, None]).mean(axis=1)
    assignments["global_k4_cluster"] = global_labels
    assignments["global_k4_agreement"] = global_agreement

    global_support_rows = []
    for cluster in FINAL_LABELS:
        values = global_agreement[global_labels == cluster]
        global_support_rows.append(
            {
                "cluster": cluster,
                "units": len(values),
                "mean_agreement": values.mean(),
                "agreement_q05": np.quantile(values, 0.05),
                "fraction_agreement_below_0.75": (values < 0.75).mean(),
            }
        )
    global_support = pd.DataFrame(global_support_rows)

    crosswalk = (
        assignments.groupby(
            ["constrained_k4_cluster", "global_k4_cluster"], observed=True
        )
        .size()
        .rename("units")
        .reset_index()
    )
    k3_global_crosswalk = (
        assignments.groupby(["k3_cluster", "global_k4_cluster"], observed=True)
        .size()
        .rename("units")
        .reset_index()
    )

    target_rows = []
    for target in TARGET_INTERVALS:
        part = assignments.loc[
            assignments["target_interval"].eq(target["target_interval"])
        ]
        counts = part["constrained_k4_cluster"].value_counts()
        target_rows.append(
            {
                "target_interval": target["target_interval"],
                "original_chrom": target["original_chrom"],
                "interval_start": target["start"],
                "interval_end": target["end"],
                "acres": len(part),
                "arrays": part["array_id"].nunique(),
                **{cluster: int(counts.get(cluster, 0)) for cluster in FINAL_LABELS},
            }
        )
    target_summary = pd.DataFrame(target_rows)

    assignments.to_csv(
        RESULTS / "constrained_and_global_k4_unit_assignments.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    constrained_summary.to_csv(
        RESULTS / "constrained_k4_cluster_summary.tsv", sep="\t", index=False
    )
    constrained_by_chrom_hap.to_csv(
        RESULTS / "constrained_k4_cluster_by_chromosome_haplotype.tsv",
        sep="\t",
        index=False,
    )
    constrained_support.to_csv(
        RESULTS / "constrained_k4_support.tsv", sep="\t", index=False
    )
    daughter_runs.to_csv(RESULTS / "conditional_h3_split_runs.tsv", sep="\t", index=False)
    global_runs_initial.to_csv(RESULTS / "global_k4_runs.tsv", sep="\t", index=False)
    global_support.to_csv(RESULTS / "global_k4_support.tsv", sep="\t", index=False)
    crosswalk.to_csv(
        RESULTS / "global_vs_constrained_k4_crosswalk.tsv", sep="\t", index=False
    )
    k3_global_crosswalk.to_csv(
        RESULTS / "original_k3_vs_global_k4_crosswalk.tsv", sep="\t", index=False
    )
    target_summary.to_csv(
        RESULTS / "target_interval_overlaps.tsv", sep="\t", index=False
    )

    # The PCA coordinates and exact explained fractions are inherited from K=3 provenance.
    with (SOURCE / "provenance" / "analysis_metadata.json").open() as handle:
        k3_metadata = json.load(handle)
    explained = tuple(float(value) for value in k3_metadata["pca_explained_variance"][:2])
    make_standalone_pca(assignments, explained)
    make_figure(assignments, constrained_by_chrom_hap, constrained_summary, explained)

    exact_fraction = float(np.mean(constrained == global_labels))
    metadata = {
        "created_utc": utc_now(),
        "source_analysis": str(SOURCE),
        "source_units": int(len(assignments)),
        "raw_features": int(features.shape[1] - 3),
        "nonconstant_features": int(usable.sum()),
        "preprocessing": "sample-SD feature standardization followed by row L2 normalization",
        "weakest_k3_cluster": weakest_cluster,
        "source_weakest_k3_cluster": weakest_source_cluster,
        "weakest_k3_units": int(parent_mask.sum()),
        "weakest_k3_mean_agreement": float(weakest["mean_agreement"]),
        "daughter_labels": list(daughter_names),
        "daughter_sizes": [int((daughter_labels == value).sum()) for value in range(2)],
        "daughter_fit": daughter_fit,
        "split_one_cluster_inertia": one_cluster_inertia,
        "split_inertia_reduction_fraction": float(split_inertia_reduction),
        "split_silhouette": split_silhouette,
        "global_fit": global_fit,
        "global_vs_constrained_ari": float(adjusted_rand_score(constrained, global_labels)),
        "global_vs_constrained_exact_fraction": exact_fraction,
        "global_k4_changed_from_constrained": int(np.sum(constrained != global_labels)),
        "target_intervals": target_summary.to_dict(orient="records"),
        "source_sha256": {
            "features": sha256_file(FEATURES),
            "k3_assignments": sha256_file(K3_ASSIGNMENTS),
            "k3_support": sha256_file(K3_SUPPORT),
        },
        "label_rule": (
            "Human labels use H1-H4. The first-occurring conditional daughter "
            "(49 units; S1-rich) is H4 and the second daughter (76 units; S3-rich) "
            "retains H3; global K=4 labels "
            "maximize overlap with the constrained labels"
        ),
    }
    (PROVENANCE / "analysis_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    write_results_markdown(
        constrained_summary,
        constrained_by_chrom_hap,
        constrained_support,
        global_support,
        crosswalk,
        metadata,
    )

    print(constrained_summary.to_string(index=False))
    print("\nConstrained support")
    print(constrained_support.to_string(index=False))
    print("\nGlobal versus constrained K=4")
    print(crosswalk.to_string(index=False))
    print(f"\nCompleted: {ROOT}")


if __name__ == "__main__":
    main()
