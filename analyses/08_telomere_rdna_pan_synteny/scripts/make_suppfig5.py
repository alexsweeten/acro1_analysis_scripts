#!/usr/bin/env python3
"""Build Supplementary Figure 5 from the audited Pan synteny summaries."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
ANALYSIS = ROOT / "acro_telomere_rdna_position_analysis"
OUT_ROOT = ROOT / "pan_rdna_inversion_suppfig5"
FIG_DIR = OUT_ROOT / "figures"
RESULT_DIR = OUT_ROOT / "results"

ORIENTATION_PATH = ANALYSIS / "pan_synteny_orientation_by_homolog.tsv"
SIDE_PATH = ANALYSIS / "pan_rdna_side_summary.tsv"

FIG_DIR.mkdir(parents=True, exist_ok=True)
RESULT_DIR.mkdir(parents=True, exist_ok=True)


def configure_plotting() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
            "font.size": 12,
            "axes.titlesize": 17,
            "axes.titleweight": "bold",
            "axes.labelsize": 13,
            "xtick.labelsize": 11.5,
            "ytick.labelsize": 11.5,
            "legend.fontsize": 11.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def load_orientation() -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = pd.read_csv(ORIENTATION_PATH, sep="\t")
    raw["total_matching_bp"] = raw["forward_matching_bp"] + raw["inverted_matching_bp"]
    raw["forward_fraction"] = raw["forward_matching_bp"] / raw["total_matching_bp"]
    raw["inverted_fraction"] = raw["inverted_matching_bp"] / raw["total_matching_bp"]

    pooled = (
        raw.groupby("species", as_index=False)[
            [
                "forward_blocks",
                "inverted_blocks",
                "forward_matching_bp",
                "inverted_matching_bp",
                "total_matching_bp",
            ]
        ]
        .sum()
    )
    pooled["forward_fraction"] = pooled["forward_matching_bp"] / pooled["total_matching_bp"]
    pooled["inverted_fraction"] = pooled["inverted_matching_bp"] / pooled["total_matching_bp"]
    pooled["homolog"] = "Pooled"

    raw.to_csv(RESULT_DIR / "orientation_fractions_by_homolog.tsv", sep="\t", index=False)
    pooled.to_csv(RESULT_DIR / "orientation_fractions_pooled.tsv", sep="\t", index=False)
    return raw, pooled


def load_repeat_sides() -> pd.DataFrame:
    side = pd.read_csv(SIDE_PATH, sep="\t")
    side = side.copy()
    side["distal_fraction"] = side["distal_bp"] / side["total_classified_bp"]
    side["feature_label"] = side["feature"].map(
        {
            "ACRO": "ACRO1-annotated sequence",
            "HSat1A": "HSat1A",
            "ACRO_plus_HSat1A": "ACRO1 + HSat1A",
        }
    )
    side.to_csv(RESULT_DIR / "repeat_burden_relative_to_rdna.tsv", sep="\t", index=False)
    return side


def draw_arrow_array(ax, x0, x1, y, direction, color, edge="#3A2632") -> None:
    """Draw a stylized tandem array with repeated directional arrow units."""
    n = 8
    unit = (x1 - x0) / n
    for j in range(n):
        left = x0 + j * unit
        right = left + unit * 0.93
        if direction == "right":
            start, end = (left + 0.08, y), (right - 0.08, y)
        else:
            start, end = (right - 0.08, y), (left + 0.08, y)
        arrow = FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=17,
            linewidth=1.0,
            facecolor=color,
            edgecolor=edge,
        )
        ax.add_patch(arrow)


def panel_a(ax) -> None:
    ax.set_xlim(0, 100)
    ax.set_ylim(-0.15, 2.15)
    ax.axis("off")

    green = "#2E8B57"
    purple = "#6F3FA0"
    rdna = "#9B3A53"
    other = "#D7DDE1"
    outline = "#20282E"

    ax.text(-0.015, 1.08, "A", transform=ax.transAxes, fontsize=23, fontweight="bold", va="top")
    ax.text(
        0,
        2.02,
        "rDNA orientation and rDNA-relative repeat organization",
        fontsize=17,
        fontweight="bold",
        ha="left",
        va="bottom",
    )
    ax.text(5, 1.83, "p-arm telomere", fontsize=11.5, color="#4F5F68", ha="center")
    ax.text(94, 1.83, "centromere", fontsize=11.5, color="#4F5F68", ha="center")
    ax.text(18, 1.65, "distal", fontsize=11.5, color="#4F5F68", ha="center")
    ax.text(78, 1.65, "proximal", fontsize=11.5, color="#4F5F68", ha="center")

    rows = [("Human", 1.17), (r"$\it{Pan}$", 0.38)]
    for label, y in rows:
        ax.plot([5, 95], [y, y], color=outline, linewidth=2.2, zorder=1)
        ax.plot([5, 5], [y - 0.16, y + 0.16], color=outline, linewidth=2.2)
        ax.add_patch(Rectangle((92.5, y - 0.22), 3.5, 0.44, facecolor="#9CA7AD", edgecolor=outline, linewidth=1.7))
        ax.text(1.2, y, label, ha="left", va="center", fontsize=14.5, fontweight="bold")

    # Human: a distal ACRO domain and centromere-directed rDNA transcription.
    ax.add_patch(Rectangle((10, 1.02), 18, 0.30, facecolor=green, edgecolor=outline, linewidth=1.1))
    ax.text(19, 1.17, "ACRO", color="white", fontsize=11.5, fontweight="bold", ha="center", va="center")
    ax.add_patch(Rectangle((30, 1.02), 8, 0.30, facecolor=other, edgecolor=outline, linewidth=1.0))
    draw_arrow_array(ax, 39, 61, 1.17, "right", rdna)
    ax.text(50, 1.45, "rDNA transcription toward centromere", fontsize=11.2, ha="center", color=rdna, fontweight="bold")
    ax.add_patch(Rectangle((62, 1.02), 28, 0.30, facecolor=other, edgecolor=outline, linewidth=1.0))

    # Pan: reversed rDNA and predominantly proximal interspersed ACRO/HSat1A burden.
    ax.add_patch(Rectangle((8, 0.23), 5, 0.30, facecolor=other, edgecolor=outline, linewidth=1.0))
    draw_arrow_array(ax, 14, 36, 0.38, "left", rdna)
    ax.text(25, 0.66, "rDNA transcription toward telomere", fontsize=11.2, ha="center", color=rdna, fontweight="bold")
    ax.add_patch(Rectangle((38, 0.23), 5, 0.30, facecolor=other, edgecolor=outline, linewidth=1.0))
    x = 44
    widths = [6.5, 5.2, 7.2, 4.7, 6.8, 5.0, 6.2, 5.4]
    colors = [green, purple, green, purple, green, purple, green, purple]
    for w, c in zip(widths, colors):
        ax.add_patch(Rectangle((x, 0.23), w, 0.30, facecolor=c, edgecolor="white", linewidth=0.8))
        x += w
    ax.add_patch(Rectangle((44, 0.23), x - 44, 0.30, facecolor="none", edgecolor=outline, linewidth=1.1))
    ax.text((44 + x) / 2, -0.02, "predominantly proximal ACRO-HSat1A burden", fontsize=11.5, ha="center", fontweight="bold")

    model = FancyArrowPatch(
        (25, 0.91),
        (59, 0.63),
        connectionstyle="arc3,rad=-0.23",
        arrowstyle="-|>",
        mutation_scale=18,
        linewidth=1.8,
        linestyle=(0, (4, 3)),
        color="#59656B",
    )
    ax.add_patch(model)
    ax.text(43, 0.84, "inversion/repositioning model", fontsize=11, color="#59656B", ha="center", fontstyle="italic")

def panel_b(ax, raw: pd.DataFrame, pooled: pd.DataFrame) -> None:
    forward = "#76A9D2"
    inverted = "#E76F51"
    no_data = "#D3D7DA"

    species_order = ["Chimpanzee", "Bonobo"]
    homolog_order = ["hsa13", "hsa14", "hsa21", "hsa22", "Pooled"]
    rows = []
    for species in species_order:
        for homolog in homolog_order:
            if homolog == "Pooled":
                z = pooled[pooled["species"] == species]
            else:
                z = raw[(raw["species"] == species) & (raw["homolog"] == homolog)]
            if z.empty:
                rows.append((species, homolog, np.nan, np.nan, 0, 0, 0))
            else:
                r = z.iloc[0]
                rows.append(
                    (
                        species,
                        homolog,
                        float(r["forward_fraction"]),
                        float(r["inverted_fraction"]),
                        int(r["total_matching_bp"]),
                        int(r["forward_blocks"]),
                        int(r["inverted_blocks"]),
                    )
                )

    y = np.arange(len(rows))[::-1]
    for yi, (species, homolog, fwd, inv, total_bp, n_fwd, n_inv) in zip(y, rows):
        if np.isnan(fwd):
            ax.barh(yi, 100, color=no_data, edgecolor="white", height=0.68)
            ax.text(50, yi, "no retained reciprocal block", ha="center", va="center", fontsize=9.7, color="#5F696E")
        else:
            lw = 1.8 if homolog == "Pooled" else 0.7
            ax.barh(yi, fwd * 100, color=forward, edgecolor="white", height=0.68, linewidth=lw)
            ax.barh(yi, inv * 100, left=fwd * 100, color=inverted, edgecolor="white", height=0.68, linewidth=lw)
            ax.text(fwd * 100 + inv * 50, yi, f"{inv * 100:.1f}%", ha="center", va="center", fontsize=10.2, fontweight="bold", color="white")
            ax.text(101.5, yi, f"{total_bp/1000:.0f} kb", ha="left", va="center", fontsize=9.8, color="#56636A")

    labels = [f"{species}  {homolog}" for species, homolog, *_ in rows]
    ax.set_yticks(y, labels)
    for tick, (_, homolog, *_rest) in zip(ax.get_yticklabels(), rows):
        if homolog == "Pooled":
            tick.set_fontweight("bold")
    ax.axvline(50, color="#49565C", linestyle=(0, (3, 3)), linewidth=1.1, zorder=0)
    ax.set_xlim(0, 112)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xlabel("Fraction of aligned matching bases (%)")
    ax.set_title("Reciprocal unique-sequence block orientation", loc="left", pad=12)
    ax.text(-0.09, 1.08, "B", transform=ax.transAxes, fontsize=23, fontweight="bold", va="top")
    ax.grid(axis="x", color="#E7EBED", linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    ax.legend(
        handles=[
            Rectangle((0, 0), 1, 1, facecolor=forward, label="Forward"),
            Rectangle((0, 0), 1, 1, facecolor=inverted, label="Inverted"),
        ],
        loc="upper right",
        bbox_to_anchor=(1.00, 1.12),
        frameon=False,
        ncol=2,
    )


def panel_c(ax, side: pd.DataFrame) -> None:
    distal = "#DCE5E6"
    proximal = "#236B6B"
    order = [
        ("Chimpanzee", "ACRO"),
        ("Chimpanzee", "HSat1A"),
        ("Chimpanzee", "ACRO_plus_HSat1A"),
        ("Bonobo", "ACRO"),
        ("Bonobo", "HSat1A"),
        ("Bonobo", "ACRO_plus_HSat1A"),
    ]
    label_map = {
        "ACRO": "ACRO1",
        "HSat1A": "HSat1A",
        "ACRO_plus_HSat1A": "ACRO1 + HSat1A",
    }
    rows = []
    for species, feature in order:
        r = side[(side["species"] == species) & (side["feature"] == feature)].iloc[0]
        rows.append((species, feature, float(r["distal_fraction"]), float(r["proximal_fraction"]), int(r["n_rdna_chromosomes"])))

    y = np.array([6, 5, 4, 2.5, 1.5, 0.5])
    for yi, (species, feature, dist, prox, n_chr) in zip(y, rows):
        ax.barh(yi, dist * 100, color=distal, edgecolor="white", height=0.70)
        ax.barh(yi, prox * 100, left=dist * 100, color=proximal, edgecolor="white", height=0.70)
        text_color = "white" if prox >= 0.28 else "#1D3034"
        ax.text(dist * 100 + prox * 50, yi, f"{prox * 100:.1f}% proximal", ha="center", va="center", fontsize=10.0, fontweight="bold", color=text_color)

    labels = [f"{species}  {label_map[feature]}" for species, feature, *_ in rows]
    ax.set_yticks(y, labels)
    ax.set_xlim(0, 100)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xlabel("Fraction of classified bases outside rDNA (%)")
    ax.set_title("Repeat burden relative to the rDNA block", loc="left", pad=12)
    ax.text(-0.09, 1.08, "C", transform=ax.transAxes, fontsize=23, fontweight="bold", va="top")
    ax.grid(axis="x", color="#E7EBED", linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    ax.legend(
        handles=[
            Rectangle((0, 0), 1, 1, facecolor=distal, label="Distal / telomeric side"),
            Rectangle((0, 0), 1, 1, facecolor=proximal, label="Proximal / centromeric side"),
        ],
        loc="upper right",
        bbox_to_anchor=(1.00, 1.12),
        frameon=False,
        ncol=2,
    )
def main() -> None:
    configure_plotting()
    orientation, pooled = load_orientation()
    side = load_repeat_sides()

    fig = plt.figure(figsize=(16, 11.5), constrained_layout=False)
    gs = fig.add_gridspec(
        2,
        2,
        height_ratios=[0.82, 1.58],
        width_ratios=[1.18, 1.0],
        hspace=0.34,
        wspace=0.34,
        left=0.15,
        right=0.965,
        top=0.88,
        bottom=0.12,
    )
    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[1, 1])

    panel_a(ax_a)
    panel_b(ax_b, orientation, pooled)
    panel_c(ax_c, side)

    fig.suptitle(
        r"$\it{Pan}$ rDNA orientation, short-arm synteny, and repeat position",
        fontsize=25,
        fontweight="bold",
        x=0.15,
        ha="left",
        y=0.965,
    )
    fig.text(
        0.15,
        0.925,
        "Observed orientation and repeat distributions are consistent with, but do not resolve, inversion-mediated repositioning.",
        fontsize=14,
        color="#526169",
        ha="left",
    )
    fig.text(
        0.15,
        0.022,
        "Synteny fractions are weighted by matching bases in CenSat-masked reciprocal unique-sequence blocks from hsa13, hsa14, hsa21, and hsa22; bonobo hsa21 had no retained block.\n"
        "Repeat positions were pooled across 8 chimpanzee and 7 bonobo rDNA-bearing haplotypes. The schematic is an evolutionary model, not a breakpoint-resolved reconstruction.",
        fontsize=10.5,
        color="#56636A",
        ha="left",
        va="bottom",
        wrap=True,
    )

    stem = FIG_DIR / "Supplementary_Figure_5_Pan_rDNA_inversion_ACRO_HSat1A"
    fig.savefig(stem.with_suffix(".pdf"), dpi=300, facecolor="white")
    fig.savefig(stem.with_suffix(".png"), dpi=300, facecolor="white")
    fig.savefig(stem.with_suffix(".svg"), dpi=300, facecolor="white")
    plt.close(fig)

    summary = {
        row.species: {
            "inverted_fraction": row.inverted_fraction,
            "forward_bp": int(row.forward_matching_bp),
            "inverted_bp": int(row.inverted_matching_bp),
            "forward_blocks": int(row.forward_blocks),
            "inverted_blocks": int(row.inverted_blocks),
        }
        for row in pooled.itertuples()
    }
    with (RESULT_DIR / "figure_values.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["species", "inverted_fraction", "forward_matching_bp", "inverted_matching_bp", "forward_blocks", "inverted_blocks"])
        for species in ["Chimpanzee", "Bonobo"]:
            d = summary[species]
            writer.writerow([species, f"{d['inverted_fraction']:.12f}", d["forward_bp"], d["inverted_bp"], d["forward_blocks"], d["inverted_blocks"]])

    print(f"Wrote {stem.with_suffix('.pdf')}")
    print(f"Wrote {stem.with_suffix('.png')}")
    print(f"Wrote {stem.with_suffix('.svg')}")
    for species in ["Chimpanzee", "Bonobo"]:
        print(f"{species}: {summary[species]['inverted_fraction']:.6%} inverted matching bases")


if __name__ == "__main__":
    main()
