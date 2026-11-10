#!/usr/bin/env python3
"""Characterize the explicit HG002 k=4 NTRiage solution at six target loci."""

from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "hg002" / "k4_target_analysis"
OUT.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".mplconfig"))
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
import matplotlib.pyplot as plt  # noqa: E402


FEATURES = ROOT / "runs" / "hg002" / "native" / "hg002.features.tsv"
UNITS = ROOT / "runs" / "hg002" / "native" / "hg002.clusters.tsv"
COMPOSITION = ROOT / "runs" / "hg002" / "native" / "hg002.acro1_composition.tsv"
SCALING = ROOT / "results" / "hg002" / "k_selection" / "feature_scaling.tsv"
ASSIGNMENTS = (
    ROOT / "results" / "hg002" / "k_selection" / "candidate_refits"
    / "array_balanced_K04_rejected" / "assignments.tsv.gz"
)
SUPPORT = (
    ROOT / "results" / "hg002" / "k_selection" / "candidate_refits"
    / "array_balanced_K04_rejected" / "final_cluster_support.tsv"
)
TARGETS = ROOT / "paired_hg002_gorilla" / "inputs" / "human_target_loci.tsv"

COLORS = {1: "#4E79A7", 2: "#F28E2B", 3: "#59A14F", 4: "#E15759"}


def main() -> None:
    features = pd.read_csv(FEATURES, sep="\t")
    units = pd.read_csv(UNITS, sep="\t")
    assignments = pd.read_csv(ASSIGNMENTS, sep="\t")
    support = pd.read_csv(SUPPORT, sep="\t")
    targets = pd.read_csv(TARGETS, sep="\t")
    composition = pd.read_csv(COMPOSITION, sep="\t")
    scaling = pd.read_csv(SCALING, sep="\t")

    data = units.merge(
        assignments[
            [
                "acre_id", "selected_cluster", "selected_cluster_label", "agreement",
                "vote_margin", "distance_to_centroid",
            ]
        ],
        on="acre_id",
        validate="one_to_one",
        suffixes=("_native", ""),
    ).merge(
        composition[
            ["acre_id", "subfamily_1", "subfamily_2", "subfamily_3", "classified_total"]
        ],
        on="acre_id",
        how="left",
        validate="one_to_one",
    )
    target_chromosomes = set(targets["analysis_chrom"])
    data["specified_target_locus"] = data["chrom"].isin(target_chromosomes)
    data["target_category"] = np.where(
        data["chrom"].eq("PF000029"),
        "chr15 paternal ~2.5 Mb",
        np.where(data["specified_target_locus"], "five telomeric loci", "other HG002"),
    )
    data["ge_7000bp"] = data["length"].ge(7000)
    data["ge_8000bp"] = data["length"].ge(8000)

    array_cluster = (
        data.groupby(["array_id", "selected_cluster"], as_index=False)
        .agg(array_cluster_mean_length_bp=("length", "mean"))
    )
    array_weighted = (
        array_cluster.groupby("selected_cluster")["array_cluster_mean_length_bp"]
        .agg(array_weighted_mean_unit_length_bp="mean",
             array_weighted_median_unit_length_bp="median")
    )
    summary = (
        data.groupby("selected_cluster")
        .agg(
            units=("acre_id", "size"),
            arrays=("array_id", "nunique"),
            mean_unit_length_bp=("length", "mean"),
            median_unit_length_bp=("length", "median"),
            sd_unit_length_bp=("length", "std"),
            min_unit_length_bp=("length", "min"),
            max_unit_length_bp=("length", "max"),
            units_ge_7000bp=("ge_7000bp", "sum"),
            fraction_units_ge_7000bp=("ge_7000bp", "mean"),
            units_ge_8000bp=("ge_8000bp", "sum"),
            fraction_units_ge_8000bp=("ge_8000bp", "mean"),
            target_units=("specified_target_locus", "sum"),
            target_fraction_within_cluster=("specified_target_locus", "mean"),
        )
        .join(array_weighted)
        .reset_index()
    )
    classified = data.groupby("selected_cluster")[["subfamily_1", "subfamily_2", "subfamily_3"]].sum()
    classified["classified_total"] = classified.sum(axis=1)
    for allele in ("subfamily_1", "subfamily_2", "subfamily_3"):
        classified[f"fraction_{allele}"] = classified[allele] / classified["classified_total"]
    summary = summary.merge(classified.reset_index(), on="selected_cluster", validate="one_to_one")
    summary = summary.merge(
        support[
            [
                "selected_cluster", "mean_agreement", "agreement_q05",
                "passes_mean_agreement",
            ]
        ],
        on="selected_cluster",
        validate="one_to_one",
    )
    summary.to_csv(OUT / "hg002_k4_cluster_summary.tsv", sep="\t", index=False)

    target_summary = (
        data.loc[data["specified_target_locus"]]
        .groupby(["chrom", "selected_cluster"], as_index=False)
        .agg(
            units=("acre_id", "size"),
            total_unit_bp=("length", "sum"),
            mean_unit_length_bp=("length", "mean"),
            median_unit_length_bp=("length", "median"),
            min_unit_length_bp=("length", "min"),
            max_unit_length_bp=("length", "max"),
            units_ge_7000bp=("ge_7000bp", "sum"),
            units_ge_8000bp=("ge_8000bp", "sum"),
        )
        .merge(targets, left_on="chrom", right_on="analysis_chrom", validate="many_to_one")
    )
    target_summary.to_csv(OUT / "target_locus_by_k4_cluster.tsv", sep="\t", index=False)

    contrast = (
        data.loc[data["specified_target_locus"]]
        .groupby("target_category", as_index=False)
        .agg(
            loci=("chrom", "nunique"),
            units=("acre_id", "size"),
            total_bp=("length", "sum"),
            mean_unit_length_bp=("length", "mean"),
            median_unit_length_bp=("length", "median"),
            sd_unit_length_bp=("length", "std"),
            min_unit_length_bp=("length", "min"),
            max_unit_length_bp=("length", "max"),
            units_ge_7000bp=("ge_7000bp", "sum"),
            units_ge_8000bp=("ge_8000bp", "sum"),
        )
    )
    contrast.to_csv(OUT / "internal_chr15_vs_telomeric_target_lengths.tsv", sep="\t", index=False)

    used = scaling["used_nonconstant"].astype(bool)
    names = scaling.loc[used, "feature"].tolist()
    x = features[names].to_numpy(dtype=np.float64)
    x = (x - scaling.loc[used, "mean"].to_numpy()) / scaling.loc[used, "sample_sd"].to_numpy()
    norms = np.linalg.norm(x, axis=1)
    x[norms > 0] /= norms[norms > 0, None]
    pca = PCA(n_components=2, svd_solver="full")
    coordinates = pca.fit_transform(x)
    coordinate_table = features[["acre_id", "array_id"]].copy()
    coordinate_table["PC1"] = coordinates[:, 0]
    coordinate_table["PC2"] = coordinates[:, 1]
    coordinate_table = coordinate_table.merge(
        data[
            [
                "acre_id", "chrom", "length", "selected_cluster", "agreement",
                "specified_target_locus", "target_category",
            ]
        ],
        on="acre_id",
        validate="one_to_one",
    )
    coordinate_table.to_csv(OUT / "hg002_k4_pca_coordinates.tsv.gz", sep="\t", index=False,
                            compression="gzip")
    data.to_csv(OUT / "hg002_k4_unit_assignments.tsv.gz", sep="\t", index=False,
                compression="gzip")

    fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.6), constrained_layout=True)
    ax = axes[0]
    for cluster in range(1, 5):
        subset = coordinate_table.loc[coordinate_table["selected_cluster"].eq(cluster)]
        ax.scatter(subset["PC1"], subset["PC2"], s=25, c=COLORS[cluster], alpha=0.68,
                   edgecolors="white", linewidths=0.3, label=f"C{cluster} (n={len(subset)})",
                   rasterized=True)
    target = coordinate_table.loc[coordinate_table["specified_target_locus"]]
    ax.scatter(target["PC1"], target["PC2"], s=48, facecolors="none", edgecolors="#B21875",
               linewidths=1.0, label="Specified loci", rasterized=True)
    ax.set_title("A  HG002 NTRiage PCA at $k$ = 4", loc="left", fontweight="bold")
    ax.set_xlabel(f"PC1 ({100*pca.explained_variance_ratio_[0]:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({100*pca.explained_variance_ratio_[1]:.1f}% variance)")
    ax.legend(frameon=True, loc="best")

    ax = axes[1]
    rng = np.random.default_rng(20261103)
    values = [data.loc[data["selected_cluster"].eq(c), "length"].to_numpy() for c in range(1, 5)]
    boxes = ax.boxplot(values, positions=range(1, 5), widths=0.58, patch_artist=True,
                       showfliers=False, medianprops={"color": "#222222", "linewidth": 1.3})
    for cluster, box in enumerate(boxes["boxes"], start=1):
        box.set_facecolor(COLORS[cluster])
        box.set_alpha(0.38)
        box.set_edgecolor(COLORS[cluster])
    for cluster in range(1, 5):
        subset = data.loc[data["selected_cluster"].eq(cluster)]
        x_all = cluster + rng.uniform(-0.18, 0.18, size=len(subset))
        ax.scatter(x_all, subset["length"], s=10, c=COLORS[cluster], alpha=0.22,
                   edgecolors="none", rasterized=True)
        marked = subset.loc[subset["specified_target_locus"]]
        x_target = cluster + rng.uniform(-0.13, 0.13, size=len(marked))
        ax.scatter(x_target, marked["length"], s=31, facecolors="#B21875",
                   edgecolors="white", linewidths=0.35, alpha=0.9, rasterized=True)
    ax.axhline(7000, color="#555555", linestyle="--", linewidth=0.9, label="7 kb")
    ax.axhline(8000, color="#888888", linestyle=":", linewidth=0.9, label="8 kb")
    ax.set_xticks(range(1, 5), [f"C{c}" for c in range(1, 5)])
    ax.set_ylabel("ACRE unit length (bp)")
    ax.set_title("B  Unit length by cluster\nMagenta: specified loci", loc="left",
                 fontweight="bold")
    ax.legend(frameon=True, loc="upper right")
    for axis in axes:
        axis.grid(color="#EBEBEB", linewidth=0.65)
        axis.set_axisbelow(True)
        for spine in ("top", "right"):
            axis.spines[spine].set_visible(False)
    fig.suptitle("HG002 four-cluster NTRiage sensitivity analysis", x=0.01, ha="left",
                 fontsize=15, fontweight="bold")
    fig.savefig(OUT / "hg002_k4_pca_and_unit_length.png", dpi=600, bbox_inches="tight",
                facecolor="white")
    fig.savefig(OUT / "hg002_k4_pca_and_unit_length.pdf", bbox_inches="tight",
                facecolor="white")
    plt.close(fig)

    telomeric = contrast.loc[contrast["target_category"].eq("five telomeric loci")].iloc[0]
    internal = contrast.loc[contrast["target_category"].eq("chr15 paternal ~2.5 Mb")].iloc[0]
    metadata_out = {
        "units": int(len(data)),
        "k": 4,
        "target_loci": sorted(target_chromosomes),
        "features_used": len(names),
        "pca_explained_variance_percent": {
            "PC1": float(100 * pca.explained_variance_ratio_[0]),
            "PC2": float(100 * pca.explained_variance_ratio_[1]),
        },
        "internal_minus_telomeric_mean_length_bp": float(
            internal["mean_unit_length_bp"] - telomeric["mean_unit_length_bp"]
        ),
        "stability_caveat": (
            "The prespecified k=4 refit was rejected because C3 mean aligned-start "
            "agreement was 0.643; C1, C2, and C4 agreements were 0.955, 0.919, "
            "and 0.982. No specified-locus unit was assigned to C3."
        ),
    }
    (OUT / "analysis_metadata.json").write_text(json.dumps(metadata_out, indent=2) + "\n")
    print(summary.to_string(index=False))
    print("\nTarget contrast")
    print(contrast.to_string(index=False))


if __name__ == "__main__":
    main()
