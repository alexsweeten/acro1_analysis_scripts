#!/usr/bin/env python3
"""Build Supplementary Figure 1: full-length ACRO1 logos across great apes.

The figure uses every sequence in the three strict ACRO1 classes from the
whole-primate callset. Each class is plotted in its complete native,
RepeatMasker/Dfam-hit-phased coordinate system (147 or 144 bp), retaining
every observed nucleotide and avoiding lossy projection of complex indels.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parent
CACHE = ROOT / ".matplotlib_cache"
CACHE.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(CACHE))
os.environ.setdefault("XDG_CACHE_HOME", str(CACHE))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Patch, PathPatch
from matplotlib.textpath import TextPath
from matplotlib.transforms import Affine2D, blended_transform_factory
import numpy as np
import pandas as pd


matplotlib.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    }
)


SOURCE = PROJECT / "orf_preserving/inputs/all_primate_acro1_subfamily_locations.tsv.gz"
SAVED_CONSENSUS = (
    PROJECT
    / "acro1_primate_allele_protein_comparison/results/"
    "species_allele_consensus_and_ctd.tsv"
)
PANELS = ROOT / "panels"
DATA = ROOT / "data"

SPECIES_ORDER = [
    "hg002",
    "mPanTro3",
    "mPanPan1",
    "mGorGor1",
    "mPonAbe1",
    "mPonPyg2",
]
SPECIES_LABELS = {
    "hg002": "Human (HG002)",
    "mPanTro3": "Chimpanzee",
    "mPanPan1": "Bonobo",
    "mGorGor1": "Gorilla",
    "mPonAbe1": "Sumatran orangutan",
    "mPonPyg2": "Bornean orangutan",
}
CLASS_ORDER = ["canonical_ACRO1", "internal_CCC", "boundary_terminal"]
CLASS_LABELS = {
    "canonical_ACRO1": "Canonical ACRO1 (147 bp)",
    "internal_CCC": "Internal-deletion ACRO1 (144 bp)",
    "boundary_terminal": "Terminal-deletion ACRO1 (144 bp)",
}
CLASS_SHORT = {
    "canonical_ACRO1": "canonical_147bp",
    "internal_CCC": "internal_deletion_144bp",
    "boundary_terminal": "terminal_deletion_144bp",
}
CLASS_PANEL = {
    "canonical_ACRO1": "A",
    "internal_CCC": "B",
    "boundary_terminal": "C",
}
CLASS_ACCENT = {
    "canonical_ACRO1": "#4E79A7",
    "internal_CCC": "#E41A1C",
    "boundary_terminal": "#FF9D00",
}
NATIVE_LENGTHS = {
    "canonical_ACRO1": 147,
    "internal_CCC": 144,
    "boundary_terminal": 144,
}
BASES = "ACGT"
BASE_COLORS = {
    "A": "#2CA02C",
    "C": "#377EB8",
    "G": "#F2A900",
    "T": "#D62728",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def plurality_base(counts: Counter[str]) -> str:
    """Return the deterministic A/C/G/T plurality base (alphabetic tie break)."""
    best = max(counts[base] for base in BASES)
    return next(base for base in BASES if counts[base] == best)


def native_summary(sequences: list[str]) -> tuple[str, list[float], list[int]]:
    length = len(sequences[0])
    if any(len(sequence) != length for sequence in sequences):
        raise RuntimeError("A species/class group contains unequal sequence lengths")
    consensus: list[str] = []
    modal_fractions: list[float] = []
    tied_positions: list[int] = []
    for position in range(length):
        counts = Counter(sequence[position] for sequence in sequences)
        valid = sum(counts[base] for base in BASES)
        if valid == 0:
            consensus.append("N")
            modal_fractions.append(0.0)
            continue
        modal = max(counts[base] for base in BASES)
        winners = [base for base in BASES if counts[base] == modal]
        consensus.append(winners[0])
        modal_fractions.append(modal / valid)
        if len(winners) > 1:
            tied_positions.append(position + 1)
    return "".join(consensus), modal_fractions, tied_positions


def native_pwm(sequences: list[str]) -> pd.DataFrame:
    """Build a full information-content PWM in native hit-phased coordinates."""
    length = len(sequences[0])
    if any(len(sequence) != length for sequence in sequences):
        raise RuntimeError("A species/class group contains unequal sequence lengths")
    rows: list[dict[str, float | int | str]] = []
    for position in range(length):
        counts = Counter(sequence[position] for sequence in sequences)
        valid = sum(counts[base] for base in BASES)
        probabilities = {
            base: counts[base] / valid if valid else 0.0 for base in BASES
        }
        entropy = -sum(
            probability * math.log2(probability)
            for probability in probabilities.values()
            if probability > 0
        )
        rows.append(
            {
                "native_position_1based": position + 1,
                **{f"p_{base}": probabilities[base] for base in BASES},
                "ambiguous_fraction": sum(
                    count
                    for base, count in counts.items()
                    if base not in BASES
                )
                / len(sequences),
                "information_bits_non_gap": max(0.0, 2.0 - entropy)
                if valid
                else 0.0,
                "n_sequences": len(sequences),
                "n_non_gap": valid,
            }
        )
    return pd.DataFrame(rows)


def frame_consensus(frame: pd.DataFrame) -> str:
    consensus: list[str] = []
    for _, row in frame.iterrows():
        probabilities = Counter(
            {base: float(row[f"p_{base}"]) for base in BASES}
        )
        consensus.append(plurality_base(probabilities))
    return "".join(consensus)


def draw_letter(axis: plt.Axes, base: str, x: float, y: float, height: float) -> None:
    if height <= 0.002:
        return
    path = TextPath(
        (0, 0),
        base,
        size=1,
        prop=FontProperties(family="Helvetica", weight="bold"),
    )
    box = path.get_extents()
    transform = (
        Affine2D()
        .translate(-box.x0, -box.y0)
        .scale(0.86 / box.width, height / box.height)
        .translate(x + 0.07, y)
    )
    axis.add_patch(
        PathPatch(
            path,
            transform=transform + axis.transData,
            color=BASE_COLORS[base],
            linewidth=0,
        )
    )


def draw_logo_axis(
    axis: plt.Axes,
    frame: pd.DataFrame,
    species_id: str,
    motif_type: str,
    show_x_axis: bool,
    compact: bool,
) -> None:
    for position, (_, row) in enumerate(frame.iterrows()):
        information = float(row["information_bits_non_gap"])
        heights = {
            base: float(row[f"p_{base}"]) * information for base in BASES
        }
        bottom = 0.0
        for base, height in sorted(heights.items(), key=lambda item: item[1]):
            draw_letter(axis, base, position, bottom, height)
            bottom += height

    native_length = len(frame)
    axis.set_xlim(0, native_length)
    axis.set_ylim(0, 2.04)
    axis.set_yticks([0, 1, 2])
    axis.tick_params(
        axis="y",
        labelsize=5.8 if compact else 7.2,
        length=2.2,
        width=0.6,
        pad=1.5,
    )
    axis.set_ylabel("bits", fontsize=6.2 if compact else 7.5, labelpad=1.5)

    interior_ticks = list(range(10, native_length, 10))
    if interior_ticks and native_length - interior_ticks[-1] < 6:
        interior_ticks.pop()
    ticks = np.asarray([1] + interior_ticks + [native_length])
    axis.set_xticks(ticks - 0.5)
    if show_x_axis:
        axis.set_xticklabels(ticks, fontsize=23.2 if compact else 7.0)
        axis.tick_params(
            axis="x",
            length=5.0 if compact else 2.5,
            width=1.0 if compact else 0.6,
            pad=7 if compact else 2,
        )
        axis.set_xlabel(
            f"Position in the {native_length}-bp native, hit-phased class coordinate",
            fontsize=27.2 if compact else 8.2,
            labelpad=11 if compact else 3,
        )
    else:
        axis.set_xticklabels([])
        axis.tick_params(axis="x", length=2.0, width=0.5)

    n_sequences = int(frame["n_sequences"].iloc[0])
    n_arrays = int(frame["n_arrays"].iloc[0])
    axis.text(
        -0.018,
        0.52,
        f"{SPECIES_LABELS[species_id]}\n"
        f"{native_length} bp  ·  {n_sequences:,} copies  ·  {n_arrays:,} arrays",
        transform=axis.transAxes,
        ha="right",
        va="center",
        fontsize=28.8 if compact else 9.1,
        fontweight="bold",
        color="#17202A",
        clip_on=False,
        linespacing=1.22,
    )

    consensus = frame_consensus(frame)
    consensus_transform = blended_transform_factory(axis.transData, axis.transAxes)
    consensus_font = 4.0 if compact else 5.2
    for position, base in enumerate(consensus):
        axis.text(
            position + 0.5,
            1.055,
            base,
            transform=consensus_transform,
            ha="center",
            va="bottom",
            fontsize=consensus_font,
            fontfamily="Helvetica",
            fontweight="bold",
            color=BASE_COLORS.get(base, "#7B8188"),
            clip_on=False,
        )

    axis.spines[["top", "right"]].set_visible(False)
    axis.spines[["left", "bottom"]].set_color("#53616B")
    axis.spines[["left", "bottom"]].set_linewidth(0.65)


def add_base_legend(figure: plt.Figure, x: float, y: float, fontsize: float) -> None:
    handles = [
        Patch(facecolor=BASE_COLORS[base], edgecolor="none", label=base)
        for base in BASES
    ]
    figure.legend(
        handles=handles,
        loc="upper right",
        bbox_to_anchor=(x, y),
        ncol=4,
        frameon=False,
        fontsize=fontsize,
        handlelength=1.0,
        handleheight=0.7,
        columnspacing=0.9,
        handletextpad=0.35,
    )


def make_class_figure(
    motif_type: str, frames: dict[tuple[str, str], pd.DataFrame]
) -> plt.Figure:
    figure = plt.figure(figsize=(18, 12.2), facecolor="white")
    grid = figure.add_gridspec(
        6,
        1,
        left=0.18,
        right=0.987,
        top=0.855,
        bottom=0.085,
        hspace=0.63,
    )
    for index, species_id in enumerate(SPECIES_ORDER):
        axis = figure.add_subplot(grid[index, 0])
        draw_logo_axis(
            axis,
            frames[(species_id, motif_type)],
            species_id,
            motif_type,
            show_x_axis=index == len(SPECIES_ORDER) - 1,
            compact=False,
        )

    panel = CLASS_PANEL[motif_type]
    figure.text(
        0.025,
        0.965,
        panel,
        ha="left",
        va="top",
        fontsize=48,
        fontweight="bold",
        color=CLASS_ACCENT[motif_type],
    )
    figure.text(
        0.065,
        0.965,
        CLASS_LABELS[motif_type],
        ha="left",
        va="top",
        fontsize=36,
        fontweight="bold",
        color="#17202A",
    )
    figure.text(
        0.065,
        0.925,
        "Full information-content logo; the native-coordinate plurality consensus is shown above each species row.",
        ha="left",
        va="top",
        fontsize=10,
        color="#48545D",
    )
    figure.text(
        0.065,
        0.895,
        "All strict-class copies are included; copy and contributing-array counts are shown at left.",
        ha="left",
        va="top",
        fontsize=9,
        color="#66727B",
    )
    add_base_legend(figure, 0.985, 0.915, 8.5)
    return figure


def make_overview_figure(
    frames: dict[tuple[str, str], pd.DataFrame]
) -> plt.Figure:
    figure = plt.figure(figsize=(30, 27), facecolor="white")
    outer = figure.add_gridspec(
        3,
        1,
        left=0.40,
        right=0.99,
        top=0.895,
        bottom=0.035,
        hspace=0.18,
    )
    for class_index, motif_type in enumerate(CLASS_ORDER):
        inner = outer[class_index].subgridspec(6, 1, hspace=0.60)
        for species_index, species_id in enumerate(SPECIES_ORDER):
            axis = figure.add_subplot(inner[species_index, 0])
            draw_logo_axis(
                axis,
                frames[(species_id, motif_type)],
                species_id,
                motif_type,
                show_x_axis=species_index == len(SPECIES_ORDER) - 1,
                compact=True,
            )
        bounds = outer[class_index].get_position(figure)
        figure.text(
            0.025,
            bounds.y1 + 0.021,
            CLASS_PANEL[motif_type],
            ha="left",
            va="bottom",
            fontsize=38,
            fontweight="bold",
            color=CLASS_ACCENT[motif_type],
        )
        figure.text(
            0.053,
            bounds.y1 + 0.021,
            CLASS_LABELS[motif_type],
            ha="left",
            va="bottom",
            fontsize=26.4,
            fontweight="bold",
            color="#17202A",
        )

    figure.text(
        0.025,
        0.985,
        "Supplementary Figure 1 | ACRO1 consensus sequences and full-length sequence logos across the great apes",
        ha="left",
        va="top",
        fontsize=18,
        fontweight="bold",
        color="#17202A",
    )
    figure.text(
        0.025,
        0.958,
        "The colored plurality consensus is shown above each information-content logo; each class uses its complete native hit-phased coordinate.",
        ha="left",
        va="top",
        fontsize=9.5,
        color="#48545D",
    )
    add_base_legend(figure, 0.988, 0.958, 33.2)
    return figure


def write_fasta(summary: pd.DataFrame, path: Path) -> None:
    with path.open("w") as handle:
        for row in summary.itertuples(index=False):
            safe_species = row.species_label.replace(" ", "_").replace("(", "").replace(")", "")
            handle.write(
                f">{safe_species}|{row.species_id}|{row.motif_type}|"
                f"n={row.n_sequences}|length={row.native_length_bp}\n"
            )
            handle.write(f"{row.native_plurality_consensus}\n")


def main() -> None:
    PANELS.mkdir(exist_ok=True)
    DATA.mkdir(exist_ok=True)
    labels = pd.read_csv(SOURCE, sep="\t", low_memory=False)
    strict = labels.loc[
        labels["species_id"].isin(SPECIES_ORDER)
        & labels["subfamily_label"].isin(CLASS_ORDER)
    ].copy()
    strict["sequence"] = strict["sequence"].astype(str).str.upper()

    unexpected_species = sorted(set(strict["species_id"]) - set(SPECIES_ORDER))
    if unexpected_species:
        raise RuntimeError(f"Unexpected species in strict subset: {unexpected_species}")

    sequences: dict[tuple[str, str], list[str]] = {}
    native_rows: list[dict[str, object]] = []
    for species_id in SPECIES_ORDER:
        for motif_type in CLASS_ORDER:
            group_table = strict.loc[
                (strict["species_id"] == species_id)
                & (strict["subfamily_label"] == motif_type)
            ]
            group = group_table["sequence"].tolist()
            if not group:
                raise RuntimeError(f"Empty species/class group: {species_id}/{motif_type}")
            observed_lengths = sorted(set(map(len, group)))
            if observed_lengths != [NATIVE_LENGTHS[motif_type]]:
                raise RuntimeError(
                    f"Unexpected lengths for {species_id}/{motif_type}: {observed_lengths}"
                )
            invalid = sum(any(base not in BASES for base in sequence) for sequence in group)
            if invalid:
                raise RuntimeError(
                    f"Found {invalid} non-ACGT sequences in {species_id}/{motif_type}"
                )
            consensus, modal_fractions, ties = native_summary(group)
            sequences[(species_id, motif_type)] = group
            native_rows.append(
                {
                    "species_id": species_id,
                    "species_label": SPECIES_LABELS[species_id],
                    "motif_type": motif_type,
                    "motif_label": CLASS_LABELS[motif_type],
                    "n_sequences": len(group),
                    "n_arrays": int(group_table["array_id"].dropna().nunique()),
                    "native_length_bp": NATIVE_LENGTHS[motif_type],
                    "plus_strand_copies": int((group_table["strand"] == "+").sum()),
                    "minus_strand_copies": int((group_table["strand"] == "-").sum()),
                    "native_plurality_consensus": consensus,
                    "mean_modal_base_fraction": float(np.mean(modal_fractions)),
                    "minimum_modal_base_fraction": float(np.min(modal_fractions)),
                    "tied_consensus_positions_1based": ",".join(map(str, ties)),
                }
            )

    summary = pd.DataFrame(native_rows)

    frames: dict[tuple[str, str], pd.DataFrame] = {}
    pwm_tables: list[pd.DataFrame] = []
    for species_id in SPECIES_ORDER:
        for motif_type in CLASS_ORDER:
            frame = native_pwm(sequences[(species_id, motif_type)])
            n_arrays = int(
                summary.loc[
                    (summary["species_id"] == species_id)
                    & (summary["motif_type"] == motif_type),
                    "n_arrays",
                ].iloc[0]
            )
            frame.insert(0, "n_arrays", n_arrays)
            frame.insert(0, "motif_type", motif_type)
            frame.insert(0, "species_label", SPECIES_LABELS[species_id])
            frame.insert(0, "species_id", species_id)
            frames[(species_id, motif_type)] = frame
            pwm_tables.append(frame)
            rebuilt_consensus = frame_consensus(frame)
            expected_consensus = summary.loc[
                (summary["species_id"] == species_id)
                & (summary["motif_type"] == motif_type),
                "native_plurality_consensus",
            ].iloc[0]
            if rebuilt_consensus != expected_consensus:
                raise RuntimeError(
                    f"PWM consensus mismatch for {species_id}/{motif_type}"
                )

    # Confirm that this rebuild reproduces the previously reported 18 native
    # species/class consensuses and source counts exactly.
    saved = pd.read_csv(SAVED_CONSENSUS, sep="\t")
    saved = saved.rename(
        columns={
            "allele": "motif_type",
            "unit_count": "saved_n_sequences",
            "consensus_nt": "saved_native_consensus",
        }
    )[["species_id", "motif_type", "saved_n_sequences", "saved_native_consensus"]]
    validation = summary.merge(saved, on=["species_id", "motif_type"], validate="one_to_one")
    validation["count_matches_saved"] = (
        validation["n_sequences"] == validation["saved_n_sequences"]
    )
    validation["native_consensus_matches_saved"] = (
        validation["native_plurality_consensus"]
        == validation["saved_native_consensus"]
    )
    if not validation["count_matches_saved"].all():
        raise RuntimeError("At least one source count does not match the saved analysis")
    if not validation["native_consensus_matches_saved"].all():
        failed = validation.loc[
            ~validation["native_consensus_matches_saved"],
            ["species_id", "motif_type"],
        ]
        raise RuntimeError(f"Native consensus validation failed:\n{failed}")

    summary.to_csv(
        DATA / "great_ape_ACRO1_consensus_sequences.tsv", sep="\t", index=False
    )
    write_fasta(summary, DATA / "great_ape_ACRO1_consensus_sequences.fasta")
    pd.concat(pwm_tables, ignore_index=True).to_csv(
        DATA / "great_ape_ACRO1_native_coordinate_pwm.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    validation.to_csv(DATA / "consensus_validation.tsv", sep="\t", index=False)

    manifest = pd.DataFrame(
        [
            {
                "figure": "Supplementary Figure 1",
                "source_file": str(SOURCE),
                "source_sha256": sha256(SOURCE),
                "source_total_rows": len(labels),
                "strict_rows_plotted": len(strict),
                "species": len(SPECIES_ORDER),
                "classes": len(CLASS_ORDER),
                "coordinate_system": (
                    "native RepeatMasker/Dfam-hit phase within each 147/144-bp class"
                ),
                "bases_omitted_from_logos": 0,
            }
        ]
    )
    manifest.to_csv(DATA / "source_manifest.tsv", sep="\t", index=False)

    qc = {
        "status": "PASS",
        "groups_validated": int(len(validation)),
        "source_rows_plotted": int(len(strict)),
        "saved_counts_reproduced": bool(validation["count_matches_saved"].all()),
        "saved_native_consensuses_reproduced": bool(
            validation["native_consensus_matches_saved"].all()
        ),
        "bases_omitted_from_logos": 0,
        "native_lengths_by_class": NATIVE_LENGTHS,
    }
    with (DATA / "validation_summary.json").open("w") as handle:
        json.dump(qc, handle, indent=2)
        handle.write("\n")

    multi_page = ROOT / "Supplementary_Figure_1_ACRO1_consensus_logos_multipage.pdf"
    with PdfPages(multi_page) as pdf:
        for motif_type in CLASS_ORDER:
            figure = make_class_figure(motif_type, frames)
            pdf.savefig(figure, facecolor="white")
            stem = (
                f"Supplementary_Figure_1{CLASS_PANEL[motif_type]}_"
                f"{CLASS_SHORT[motif_type]}"
            )
            figure.savefig(PANELS / f"{stem}.pdf", facecolor="white")
            figure.savefig(PANELS / f"{stem}.svg", facecolor="white")
            figure.savefig(PANELS / f"{stem}.png", dpi=400, facecolor="white")
            plt.close(figure)

    overview = make_overview_figure(frames)
    overview.savefig(
        ROOT / "Supplementary_Figure_1_ACRO1_consensus_logos.pdf",
        facecolor="white",
    )
    overview.savefig(
        ROOT / "Supplementary_Figure_1_ACRO1_consensus_logos.png",
        dpi=240,
        facecolor="white",
    )
    plt.close(overview)

    counts = (
        summary.pivot(index="species_label", columns="motif_type", values="n_sequences")
        .reindex([SPECIES_LABELS[species] for species in SPECIES_ORDER])
        .reindex(columns=CLASS_ORDER)
    )
    print("Built Supplementary Figure 1 from 18 validated species/class groups")
    print(counts.to_string())
    print("Every observed nucleotide was retained in native class coordinates")


if __name__ == "__main__":
    main()
