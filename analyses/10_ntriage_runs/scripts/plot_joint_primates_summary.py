#!/usr/bin/env python3
"""Draw the focused summary figure for the six-primate joint NTRIAGE run."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
SPECIES_ORDER = [
    "hg002",
    "gorilla",
    "chimpanzee",
    "bonobo",
    "sumatran_orangutan",
    "bornean_orangutan",
]
SPECIES_LABELS = {
    "hg002": "HG002",
    "gorilla": "Gorilla",
    "chimpanzee": "Chimpanzee",
    "bonobo": "Bonobo",
    "sumatran_orangutan": "Sumatran orangutan",
    "bornean_orangutan": "Bornean orangutan",
}
CLUSTER_COLORS = {1: "#4C78A8", 2: "#F58518"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=PROJECT)
    parser.add_argument("--dataset", default="joint_primates")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/six_primates_no_hprc/figures"),
    )
    return parser.parse_args()


def resolve(root: Path, path: Path) -> Path:
    return path.resolve() if path.is_absolute() else (root / path).resolve()


def main() -> None:
    args = parse_args()
    root = args.root.resolve()
    analysis = root / "results" / args.dataset / "cluster_analysis"
    species = pd.read_csv(analysis / "species_cluster_summary.tsv", sep="\t")
    summary = pd.read_csv(analysis / "cluster_summary.tsv", sep="\t").sort_values(
        "cluster"
    )
    intervals = pd.read_csv(
        analysis / "cluster_array_bootstrap_intervals.tsv", sep="\t"
    )
    intervals = intervals.loc[
        intervals.metric.eq("mean_acro1_copy_number")
    ].sort_values("cluster")
    if summary.cluster.astype(int).tolist() != [1, 2]:
        raise ValueError("This focused figure expects the authoritative K=2 result")

    species_counts = (
        species.pivot(index="species_id", columns="cluster", values="units")
        .reindex(SPECIES_ORDER)
        .fillna(0)
    )
    for cluster in (1, 2):
        if cluster not in species_counts:
            species_counts[cluster] = 0
    species_fraction = species_counts.div(species_counts.sum(axis=1), axis=0)

    plt.rcParams.update({"font.size": 10, "axes.titlesize": 12})
    fig, axes = plt.subplots(1, 3, figsize=(15.2, 5.2), constrained_layout=True)

    y = np.arange(len(species_fraction))
    left = np.zeros(len(species_fraction))
    for cluster in (1, 2):
        values = species_fraction[cluster].to_numpy()
        axes[0].barh(
            y,
            values,
            left=left,
            color=CLUSTER_COLORS[cluster],
            label=f"Cluster {cluster}",
        )
        left += values
    axes[0].set_yticks(y, [SPECIES_LABELS[s] for s in SPECIES_ORDER])
    axes[0].invert_yaxis()
    axes[0].set_xlim(0, 1)
    axes[0].set_xlabel("Fraction of complete units")
    axes[0].set_title("A  Cluster membership by cohort", loc="left", fontweight="bold")
    axes[0].legend(frameon=False, loc="lower right")

    x = np.arange(2)
    means = intervals.array_weighted_mean.to_numpy()
    lower = means - intervals.bootstrap_95ci_low.to_numpy()
    upper = intervals.bootstrap_95ci_high.to_numpy() - means
    axes[1].bar(x, means, color=[CLUSTER_COLORS[1], CLUSTER_COLORS[2]], width=0.62)
    axes[1].errorbar(
        x,
        means,
        yerr=np.vstack([lower, upper]),
        fmt="none",
        color="black",
        capsize=4,
        linewidth=1.2,
    )
    axes[1].set_xticks(x, ["Cluster 1", "Cluster 2"])
    axes[1].set_ylabel("Raw ACRO1 anchors per unit")
    axes[1].set_title(
        "B  Array-weighted ACRO1 copy proxy", loc="left", fontweight="bold"
    )
    axes[1].text(
        0.02,
        0.98,
        "Mean and 95% array-bootstrap CI",
        transform=axes[1].transAxes,
        va="top",
        color="#555555",
        fontsize=9,
    )

    components = [
        ("mean_broad_acro1_annotation_bp", "Broad ACRO1", "#B279A2"),
        ("mean_explicit_backbone_bp", "Explicit backbone", "#59A14F"),
        ("mean_unclassified_bp", "Unclassified", "#BAB0AC"),
    ]
    bottom = np.zeros(2)
    for column, label, color in components:
        values = summary[column].to_numpy()
        axes[2].bar(x, values, bottom=bottom, label=label, color=color, width=0.62)
        bottom += values
    axes[2].set_xticks(x, ["Cluster 1", "Cluster 2"])
    axes[2].set_ylabel("Mean bp per complete unit")
    axes[2].set_title("C  Annotated sequence composition", loc="left", fontweight="bold")
    axes[2].legend(frameon=False, loc="upper right")

    for axis in axes:
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(axis="y", color="#dddddd", linewidth=0.7, alpha=0.7)
        axis.set_axisbelow(True)
    fig.suptitle(
        "Six-primate joint NTRIAGE result (HG002 + five nonhuman primate assemblies; no HPRC)",
        fontsize=15,
        fontweight="bold",
    )

    output = resolve(root, args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    png = output / "joint_primates_focused_summary.png"
    pdf = output / "joint_primates_focused_summary.pdf"
    fig.savefig(png, dpi=220, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {png}")
    print(f"Wrote {pdf}")


if __name__ == "__main__":
    main()
