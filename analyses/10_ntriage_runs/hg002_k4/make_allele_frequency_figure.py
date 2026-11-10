#!/usr/bin/env python3
"""Plot HG002 H1-H4 ACRO1 allele fractions in the Figure 1A style."""

from __future__ import annotations

import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent
MPLCONFIG = ROOT / ".mplconfig"
CACHE = ROOT / ".cache"
MPLCONFIG.mkdir(parents=True, exist_ok=True)
(CACHE / "fontconfig").mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIG))
os.environ.setdefault("XDG_CACHE_HOME", str(CACHE))

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "font.size": 13,
        "axes.unicode_minus": False,
        "pdf.use14corefonts": True,
        "ps.useafm": True,
    }
)
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


ASSIGNMENTS = ROOT / "results" / "constrained_and_global_k4_unit_assignments.tsv.gz"
SUMMARY = ROOT / "results" / "hg002_cluster_acro1_allele_fractions.tsv"
PDF_OUTPUT = ROOT / "output" / "pdf" / "hg002_cluster_acro1_allele_fractions.pdf"
PNG_OUTPUT = ROOT / "figures" / "hg002_cluster_acro1_allele_fractions.png"

CLUSTERS = ("H1", "H2", "H3", "H4")
ALLELES = (
    ("canonical", "subfamily_1", "Canonical (147 bp)", "#4E79A7"),
    ("internal_CCC_deletion", "subfamily_2", "Internal CCC deletion (144 bp)", "#E41A1C"),
    ("terminal_3bp_deletion", "subfamily_3", "Terminal 3-bp deletion (144 bp)", "#FF9D00"),
)


def build_summary(assignments: pd.DataFrame) -> pd.DataFrame:
    if set(assignments["constrained_k4_cluster"]) != set(CLUSTERS):
        raise RuntimeError("expected constrained H1-H4 assignments")
    rows = []
    for cluster in CLUSTERS:
        part = assignments.loc[assignments["constrained_k4_cluster"].eq(cluster)]
        counts = {
            key: int(part[source_column].sum())
            for key, source_column, _, _ in ALLELES
        }
        total = sum(counts.values())
        if total <= 0:
            raise RuntimeError(f"cluster {cluster} has no classified ACRO1 monomers")
        row = {
            "cluster": cluster,
            "acre_units": len(part),
            "classified_acro1_monomers": total,
        }
        for key, _, _, _ in ALLELES:
            row[f"{key}_count"] = counts[key]
            row[f"{key}_fraction"] = counts[key] / total
        rows.append(row)
    return pd.DataFrame(rows)


def label_segments(ax: plt.Axes, summary: pd.DataFrame, bar_width: float) -> None:
    x = np.arange(len(summary), dtype=float)
    bottoms = np.zeros(len(summary), dtype=float)
    for key, _, _, color in ALLELES:
        fractions = summary[f"{key}_fraction"].to_numpy(dtype=float)
        for index, fraction in enumerate(fractions):
            if fraction <= 0:
                continue
            label = f"{fraction * 100:.1f}%"
            midpoint = bottoms[index] + fraction / 2
            if fraction >= 0.045:
                ax.text(
                    x[index], midpoint, label,
                    ha="center", va="center", color="white", fontsize=14,
                )
            elif key == "terminal_3bp_deletion":
                ax.plot(
                    [x[index], x[index]], [1.0, 1.025],
                    color=color, linewidth=1.1, clip_on=False,
                )
                ax.text(
                    x[index], 1.035, label,
                    ha="center", va="bottom", color=color, fontsize=13,
                    clip_on=False,
                )
            else:
                x_edge = x[index] + bar_width / 2
                ax.plot(
                    [x_edge, x_edge + 0.055], [midpoint, midpoint],
                    color=color, linewidth=1.1, clip_on=False,
                )
                ax.text(
                    x_edge + 0.065, midpoint, label,
                    ha="left", va="center", color=color, fontsize=12,
                    clip_on=False,
                )
        bottoms += fractions


def make_figure(summary: pd.DataFrame) -> None:
    PDF_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    PNG_OUTPUT.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(12.0, 9.0))
    x = np.arange(len(summary), dtype=float)
    bar_width = 0.76
    bottoms = np.zeros(len(summary), dtype=float)
    for key, _, _, color in ALLELES:
        fractions = summary[f"{key}_fraction"].to_numpy(dtype=float)
        ax.bar(
            x,
            fractions,
            width=bar_width,
            bottom=bottoms,
            color=color,
            edgecolor="white",
            linewidth=0.9,
            zorder=3,
        )
        bottoms += fractions

    label_segments(ax, summary, bar_width)

    ax.set_ylim(0, 1.0)
    ax.set_xlim(-0.55, len(summary) - 0.45)
    ax.set_yticks(np.linspace(0, 1, 6))
    ax.set_yticklabels([f"{value:.0%}" for value in np.linspace(0, 1, 6)], fontsize=13)
    ax.set_xticks(x)
    ax.set_xticklabels(
        [
            f"{row.cluster}\nn = {int(row.classified_acro1_monomers):,}"
            for row in summary.itertuples(index=False)
        ],
        fontsize=13,
        linespacing=1.25,
    )
    ax.tick_params(axis="x", pad=14, length=0)
    ax.tick_params(axis="y", pad=7, width=0.9, length=6)
    ax.set_ylabel(
        "Fraction of classified ACRO1 monomers",
        fontsize=17,
        labelpad=16,
    )
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, color="#D8E0E5", linewidth=0.9)
    ax.xaxis.grid(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#222222")
    ax.spines["bottom"].set_color("#222222")
    ax.spines["left"].set_linewidth(0.9)
    ax.spines["bottom"].set_linewidth(0.9)

    fig.suptitle(
        "ACRO1 allele fractions across HG002 clusters",
        x=0.5,
        y=0.965,
        ha="center",
        fontsize=23,
        fontweight="normal",
    )
    fig.text(
        0.5,
        0.917,
        "Strictly classified 144/147-bp monomers; other spacings and ambiguous "
        "144-bp calls excluded",
        ha="center",
        va="center",
        fontsize=13,
        color="#65717B",
    )
    legend_handles = [
        Line2D(
            [0], [0], marker="s", linestyle="none", markersize=8,
            markerfacecolor=color, markeredgecolor=color, label=label,
        )
        for _, _, label, color in ALLELES
    ]
    fig.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.86),
        ncol=3,
        frameon=False,
        fontsize=13,
        handletextpad=0.6,
        columnspacing=2.2,
    )

    fig.subplots_adjust(left=0.13, right=0.97, bottom=0.14, top=0.76)
    fig.savefig(PDF_OUTPUT, facecolor="white")
    fig.savefig(PNG_OUTPUT, dpi=300, facecolor="white")
    plt.close(fig)


def main() -> None:
    if not ASSIGNMENTS.is_file():
        raise FileNotFoundError(ASSIGNMENTS)
    assignments = pd.read_csv(ASSIGNMENTS, sep="\t")
    summary = build_summary(assignments)
    summary.to_csv(SUMMARY, sep="\t", index=False)
    make_figure(summary)
    print(summary.to_string(index=False))
    print(f"PDF: {PDF_OUTPUT}")
    print(f"PNG: {PNG_OUTPUT}")


if __name__ == "__main__":
    main()
