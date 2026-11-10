#!/usr/bin/env python3
"""Compare gorilla-only clusters with active-centromere and rDNA/NOR distance."""

from __future__ import annotations

import itertools
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import seaborn as sns


ROOT = Path(__file__).resolve().parents[1]
OUT = Path(__file__).resolve().parent
ACRES = ROOT / "acro_primate_individual/gorilla_all_acres/all_gorilla_annotations.tsv"
METADATA = ROOT / "acro_gorilla_only/gorilla_tip_metadata.base.tsv"
CLUSTERS = OUT / "gorilla_cluster_assignments.tsv.gz"
CENTROMERES = ROOT / "spatial_centromere/centromere_intervals.tsv"
RDNA = ROOT / "spatial_centromere/rdna_nor_intervals.tsv"

ORDER = ["G4", "G5", "G6"]
COLORS = {"G4": "#E15759", "G5": "#B07AA1", "G6": "#76B7B2"}


def interval_distance(start: int, end: int, feature_start: int, feature_end: int) -> int:
    """Zero-based half-open edge-to-edge distance between two intervals."""
    if end <= feature_start:
        return feature_start - end
    if feature_end <= start:
        return start - feature_end
    return 0


def annotate_distances() -> pd.DataFrame:
    acres = pd.read_csv(ACRES, sep="\t")
    metadata = pd.read_csv(METADATA, sep="\t")[["acre_id", "synteny_group", "nor_status"]]
    clusters = pd.read_csv(CLUSTERS, sep="\t")[["acre_id", "gorilla_cluster"]]
    frame = acres.merge(metadata, on="acre_id", validate="one_to_one").merge(
        clusters, on="acre_id", validate="one_to_one"
    )

    centromeres = pd.read_csv(CENTROMERES, sep="\t")
    centromeres = centromeres[centromeres["species_id"] == "mGorGor1"]
    centromere_lookup = {
        row.original_chrom: (int(row.feature_start), int(row.feature_end))
        for row in centromeres.itertuples(index=False)
    }
    rdna = pd.read_csv(RDNA, sep="\t")
    rdna = rdna[rdna["species_id"] == "mGorGor1"]
    rdna_lookup = {
        chrom: list(zip(group["feature_start"].astype(int), group["feature_end"].astype(int)))
        for chrom, group in rdna.groupby("original_chrom")
    }

    centromere_distance = []
    rdna_distance = []
    rdna_site_count = []
    for row in frame.itertuples(index=False):
        centromere = centromere_lookup.get(row.original_chrom)
        if centromere is None:
            centromere_distance.append(np.nan)
        else:
            centromere_distance.append(
                interval_distance(
                    int(row.original_acre_start), int(row.original_acre_end), *centromere
                )
            )
        sites = rdna_lookup.get(row.original_chrom, [])
        rdna_site_count.append(len(sites))
        if not sites:
            rdna_distance.append(np.nan)
        else:
            rdna_distance.append(
                min(
                    interval_distance(
                        int(row.original_acre_start), int(row.original_acre_end), *site
                    )
                    for site in sites
                )
            )

    frame["centromere_distance_bp"] = centromere_distance
    frame["rdna_distance_bp"] = rdna_distance
    frame["rdna_sites_on_sequence"] = rdna_site_count
    frame["centromere_distance_mb"] = frame["centromere_distance_bp"] / 1e6
    frame["rdna_distance_mb"] = frame["rdna_distance_bp"] / 1e6
    return frame


def cluster_summary(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for cluster, group in frame.groupby("gorilla_cluster", sort=True):
        row: dict[str, object] = {
            "gorilla_cluster": cluster,
            "acres": len(group),
            "parent_arrays": group["array_id"].nunique(),
            "synteny_groups": group["synteny_group"].nunique(),
        }
        for metric in ["centromere_distance_mb", "rdna_distance_mb"]:
            prefix = metric.removesuffix("_distance_mb")
            values = group[metric].dropna()
            row.update(
                {
                    f"{prefix}_mapped_acres": len(values),
                    f"{prefix}_mean_distance_mb": values.mean(),
                    f"{prefix}_median_distance_mb": values.median(),
                    f"{prefix}_q25_distance_mb": values.quantile(0.25),
                    f"{prefix}_q75_distance_mb": values.quantile(0.75),
                    f"{prefix}_min_distance_mb": values.min(),
                    f"{prefix}_max_distance_mb": values.max(),
                }
            )
        rows.append(row)
    return pd.DataFrame(rows)


def array_summary(frame: pd.DataFrame) -> pd.DataFrame:
    return (
        frame.groupby(
            ["gorilla_cluster", "array_id", "original_chrom", "synteny_group", "nor_status"],
            dropna=False,
        )
        .agg(
            acres=("acre_id", "size"),
            mean_centromere_distance_mb=("centromere_distance_mb", "mean"),
            median_centromere_distance_mb=("centromere_distance_mb", "median"),
            mean_rdna_distance_mb=("rdna_distance_mb", "mean"),
            median_rdna_distance_mb=("rdna_distance_mb", "median"),
        )
        .reset_index()
    )


def exact_array_permutation(frame: pd.DataFrame, metric: str) -> dict[str, float | int | str]:
    subset = frame[frame["gorilla_cluster"].isin(ORDER)]
    arrays = subset.groupby(["gorilla_cluster", "array_id"])[metric].mean().dropna().reset_index()
    values = arrays[metric].to_numpy(float)
    labels = arrays["gorilla_cluster"].to_numpy(str)
    group_sizes = {cluster: int((labels == cluster).sum()) for cluster in ORDER}
    if len(values) != sum(group_sizes.values()):
        raise ValueError("Unexpected missing array mean in exact permutation test")

    def between_group_ss(test_labels: np.ndarray) -> float:
        grand_mean = values.mean()
        return float(
            sum(
                np.sum(test_labels == cluster)
                * (values[test_labels == cluster].mean() - grand_mean) ** 2
                for cluster in ORDER
            )
        )

    observed = between_group_ss(labels)
    indices = set(range(len(values)))
    null = []
    for g4_indices in itertools.combinations(range(len(values)), group_sizes["G4"]):
        remaining = sorted(indices - set(g4_indices))
        for g5_indices in itertools.combinations(remaining, group_sizes["G5"]):
            permuted = np.full(len(values), "G6", dtype=object)
            permuted[list(g4_indices)] = "G4"
            permuted[list(g5_indices)] = "G5"
            null.append(between_group_ss(permuted))
    null = np.asarray(null)
    return {
        "metric": metric,
        "arrays": len(values),
        "observed_between_group_ss": observed,
        "exact_permutations": len(null),
        "permutations_at_least_observed": int(np.sum(null >= observed - 1e-12)),
        "exact_p_value": float(np.mean(null >= observed - 1e-12)),
    }


def plot_violin(frame: pd.DataFrame, metric: str, title: str, stem: str) -> None:
    subset = frame[frame["gorilla_cluster"].isin(ORDER)].dropna(subset=[metric]).copy()
    array_means = (
        subset.groupby(["gorilla_cluster", "array_id"], as_index=False)[metric].mean()
    )
    sns.set_theme(style="whitegrid", context="talk")
    fig, axis = plt.subplots(figsize=(9.5, 8.2), constrained_layout=True)
    sns.violinplot(
        data=subset,
        x="gorilla_cluster",
        y=metric,
        order=ORDER,
        palette=COLORS,
        hue="gorilla_cluster",
        hue_order=ORDER,
        legend=False,
        inner=None,
        cut=0,
        density_norm="width",
        linewidth=1.5,
        saturation=0.9,
        ax=axis,
    )
    sns.stripplot(
        data=subset,
        x="gorilla_cluster",
        y=metric,
        order=ORDER,
        color="#17202A",
        alpha=0.18,
        jitter=0.17,
        size=2.8,
        ax=axis,
    )

    for cluster_index, cluster in enumerate(ORDER):
        cluster_arrays = array_means[array_means["gorilla_cluster"] == cluster]
        offsets = np.linspace(-0.10, 0.10, len(cluster_arrays)) if len(cluster_arrays) > 1 else np.array([0.0])
        axis.scatter(
            cluster_index + offsets,
            cluster_arrays[metric],
            marker="D",
            s=70,
            facecolor="white",
            edgecolor="black",
            linewidth=1.25,
            zorder=6,
        )
        cluster_mean = subset.loc[subset["gorilla_cluster"] == cluster, metric].mean()
        axis.hlines(cluster_mean, cluster_index - 0.22, cluster_index + 0.22, color="black", lw=3.2, zorder=7)
        axis.annotate(
            f"mean {cluster_mean:.2f} Mb",
            (cluster_index, cluster_mean),
            xytext=(0, 10),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=11,
            fontweight="bold",
        )

    axis.set_xlabel("Gorilla-only N-TRIAGE cluster")
    axis.set_ylabel("Nearest-edge distance (Mb)")
    axis.set_title(title, fontsize=17, fontweight="bold", pad=14)
    axis.legend(
        handles=[
            Line2D([0], [0], marker="D", linestyle="", markerfacecolor="white", markeredgecolor="black", markersize=8, label="Parent-array mean"),
            Line2D([0], [0], color="black", lw=3.2, label="Cluster ACRE-weighted mean"),
        ],
        loc="upper right",
        frameon=True,
        fontsize=10,
    )
    axis.text(
        0.01,
        0.01,
        "Violin and small points: individual ACREs; distance is interval edge-to-edge.",
        transform=axis.transAxes,
        ha="left",
        va="bottom",
        fontsize=9.5,
        color="#4B5563",
    )
    sns.despine(ax=axis)
    fig.savefig(OUT / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    frame = annotate_distances()
    summary = cluster_summary(frame)
    arrays = array_summary(frame)
    tests = [
        exact_array_permutation(frame, "centromere_distance_mb"),
        exact_array_permutation(frame, "rdna_distance_mb"),
    ]
    frame.to_csv(OUT / "gorilla_cluster_spatial_distances.tsv.gz", sep="\t", index=False, compression="gzip")
    summary.to_csv(OUT / "gorilla_cluster_spatial_summary.tsv", sep="\t", index=False)
    arrays.to_csv(OUT / "gorilla_cluster_array_spatial_summary.tsv", sep="\t", index=False)
    pd.DataFrame(tests).to_csv(OUT / "gorilla_G4_G6_spatial_permutation_tests.tsv", sep="\t", index=False)
    (OUT / "gorilla_cluster_spatial_summary.json").write_text(
        json.dumps({"g4_g6_exact_array_permutation_tests": tests}, indent=2) + "\n"
    )
    plot_violin(
        frame,
        "centromere_distance_mb",
        "Gorilla G4–G6 distance from the active centromere",
        "gorilla_G4_G6_centromere_distance_violin",
    )
    plot_violin(
        frame,
        "rdna_distance_mb",
        "Gorilla G4–G6 distance from the nearest rDNA/NOR",
        "gorilla_G4_G6_rdna_distance_violin",
    )
    print(summary.to_string(index=False))
    print(pd.DataFrame(tests).to_string(index=False))


if __name__ == "__main__":
    main()
