#!/usr/bin/env python3
"""Create a PCA figure for the authoritative HG002 NTRiage clustering."""

from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd


ANALYSIS = Path(__file__).resolve().parents[1]
OUTPUT = ANALYSIS / "results" / "hg002" / "pca"
OUTPUT.mkdir(parents=True, exist_ok=True)

os.environ.setdefault("MPLCONFIGDIR", str(ANALYSIS / ".mplconfig"))
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib.pyplot as plt  # noqa: E402
from sklearn.decomposition import PCA  # noqa: E402


FEATURES = ANALYSIS / "runs" / "hg002" / "native" / "hg002.features.tsv"
SCALING = ANALYSIS / "results" / "hg002" / "k_selection" / "feature_scaling.tsv"
ASSIGNMENTS = (
    ANALYSIS / "results" / "hg002" / "k_selection" / "assignments.tsv.gz"
)
UNITS = ANALYSIS / "runs" / "hg002" / "native" / "hg002.clusters.tsv"

COLORS = {1: "#4E79A7", 2: "#F28E2B"}


def main() -> None:
    features = pd.read_csv(FEATURES, sep="\t")
    scaling = pd.read_csv(SCALING, sep="\t")
    assignments = pd.read_csv(ASSIGNMENTS, sep="\t")
    units = pd.read_csv(
        UNITS,
        sep="\t",
        usecols=["acre_id", "unit_in_array", "chrom", "start", "end", "length"],
    )

    if len(features) != len(assignments):
        raise ValueError("feature and assignment row counts differ")
    if features["acre_id"].duplicated().any() or assignments["acre_id"].duplicated().any():
        raise ValueError("duplicate ACRE identifier")

    used = scaling.loc[scaling["used_nonconstant"].astype(bool)].copy()
    names = used["feature"].tolist()
    missing = sorted(set(names) - set(features.columns))
    if missing:
        raise ValueError(f"features missing scaling columns: {missing[:5]}")

    raw = features[names].to_numpy(dtype=np.float64)
    means = used["mean"].to_numpy(dtype=np.float64)
    sample_sd = used["sample_sd"].to_numpy(dtype=np.float64)
    transformed = (raw - means) / sample_sd
    row_norm = np.linalg.norm(transformed, axis=1)
    nonzero = row_norm > 0
    transformed[nonzero] /= row_norm[nonzero, None]
    transformed[~nonzero] = 0.0

    pca = PCA(n_components=2, svd_solver="full")
    coordinates = pca.fit_transform(transformed)

    plot_data = features[["acre_id", "array_id"]].copy()
    plot_data["PC1"] = coordinates[:, 0]
    plot_data["PC2"] = coordinates[:, 1]
    plot_data = plot_data.merge(
        assignments[
            [
                "acre_id",
                "selected_cluster",
                "selected_cluster_label",
                "agreement",
                "vote_margin",
            ]
        ],
        on="acre_id",
        how="left",
        validate="one_to_one",
    ).merge(units, on="acre_id", how="left", validate="one_to_one")
    if plot_data["selected_cluster"].isna().any():
        raise ValueError("unmatched authoritative cluster assignments")

    plot_data.to_csv(OUTPUT / "hg002_ntriage_k2_pca_coordinates.tsv", sep="\t", index=False)

    explained = 100.0 * pca.explained_variance_ratio_
    counts = plot_data["selected_cluster"].value_counts().sort_index()

    fig, ax = plt.subplots(figsize=(7.2, 5.7), constrained_layout=True)
    for cluster in sorted(counts.index):
        subset = plot_data.loc[plot_data["selected_cluster"].eq(cluster)]
        ax.scatter(
            subset["PC1"],
            subset["PC2"],
            s=25,
            alpha=0.76,
            c=COLORS[int(cluster)],
            edgecolors="white",
            linewidths=0.35,
            label=f"C{int(cluster)} (n = {len(subset):,})",
            rasterized=True,
            zorder=2,
        )
        centroid = subset[["PC1", "PC2"]].mean().to_numpy()
        ax.scatter(
            centroid[0],
            centroid[1],
            marker="X",
            s=105,
            c=COLORS[int(cluster)],
            edgecolors="black",
            linewidths=0.7,
            zorder=4,
        )
        ax.annotate(
            f"C{int(cluster)}",
            xy=centroid,
            xytext=(6, 6),
            textcoords="offset points",
            fontsize=9,
            fontweight="bold",
            color="#222222",
            zorder=5,
        )

    ax.axhline(0, color="#D8D8D8", linewidth=0.7, zorder=0)
    ax.axvline(0, color="#D8D8D8", linewidth=0.7, zorder=0)
    ax.grid(color="#EAEAEA", linewidth=0.65, zorder=0)
    ax.set_axisbelow(True)
    ax.set_xlabel(f"PC1 ({explained[0]:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({explained[1]:.1f}% variance)")
    ax.set_title(
        "HG002 NTRiage feature-space PCA\n"
        "Authoritative $k$ = 2; 533 complete ACRE units",
        loc="left",
        fontsize=13,
        fontweight="bold",
        pad=10,
    )
    legend = ax.legend(
        title="NTRiage cluster",
        frameon=True,
        facecolor="white",
        edgecolor="#CFCFCF",
        loc="best",
    )
    legend._legend_box.align = "left"
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.spines["left"].set_color("#888888")
    ax.spines["bottom"].set_color("#888888")
    ax.text(
        0.995,
        0.01,
        "PCA is for visualization; clustering used the full normalized feature space.",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=7.5,
        color="#666666",
    )

    png = OUTPUT / "hg002_ntriage_k2_pca.png"
    pdf = OUTPUT / "hg002_ntriage_k2_pca.pdf"
    fig.savefig(png, dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    metadata = {
        "dataset": "HG002",
        "units": int(len(plot_data)),
        "arrays": int(plot_data["array_id"].nunique()),
        "authoritative_k": 2,
        "features_before_constant_filter": int(len(scaling)),
        "features_used": int(len(names)),
        "preprocessing": (
            "saved full-data sample-SD column standardization followed by "
            "row L2 normalization"
        ),
        "pca_fit": "all 533 HG002 complete units; sklearn PCA svd_solver=full",
        "explained_variance_percent": {
            "PC1": float(explained[0]),
            "PC2": float(explained[1]),
        },
        "cluster_counts": {f"C{int(k)}": int(v) for k, v in counts.items()},
        "interpretation": (
            "PCA was fitted for visualization and was not used to assign clusters. "
            "Colors are the authoritative stability-filtered k=2 assignments; the "
            "rejected k=4 candidate is not shown."
        ),
    }
    (OUTPUT / "hg002_ntriage_k2_pca_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n"
    )
    (OUTPUT / "CAPTION.md").write_text(
        "**Principal-component projection of the authoritative HG002 NTRiage "
        "partition.** Complete ACRE units (n = 533) were projected from the same "
        f"{len(names)} nonconstant, standardized and row-L2-normalized features used for "
        "clustering. Points are colored by the stability-filtered two-cluster "
        "assignment; X symbols mark cluster centroids. PC1 and PC2 explain "
        f"{explained[0]:.1f}% and {explained[1]:.1f}% of feature-space variance, "
        "respectively. PCA was used only for visualization and did not contribute "
        "to cluster assignment.\n"
    )

    print(f"Wrote {png}")
    print(f"Wrote {pdf}")
    print(f"Explained variance: PC1={explained[0]:.3f}%, PC2={explained[1]:.3f}%")


if __name__ == "__main__":
    main()
