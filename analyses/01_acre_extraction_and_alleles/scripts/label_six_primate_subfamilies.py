#!/usr/bin/env python3
"""Label canonical, internal-CCC, and boundary/terminal ACRO1 across six primates."""

from __future__ import annotations

import bisect
import gzip
import json
import math
import re
from collections import Counter
from pathlib import Path

import edlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "primate_whole_subfamily_analysis"
ARRAY_MAPPING = ROOT / "primate_full/six_primate_array_mapping.tsv"
ACRE_METADATA = ROOT / "nor_analysis/acre_nor_annotations_full.tsv"
MOTIF_SUMMARY = ROOT / "gorilla_acro1_three_motifs/acro1_three_type_motif_summary.tsv"
GORILLA_CLASSIFICATIONS = ROOT / "gorilla_acro1_144/gorilla_acro1_144_sequence_classifications.tsv"

SPECIES_ORDER = ["hg002", "mPanPan1", "mPanTro3", "mGorGor1", "mPonPyg2", "mPonAbe1"]
SPECIES_CONFIG = {
    "hg002": {
        "name": "Human (HG002)",
        "assembly": "HG002v1.1",
        "genome": ROOT / "primate_genomes/hg002v1.1.fasta.gz",
        "bed": ROOT / "primate_acro_bedfiles/HG002/acro_monomers_full_HG002.P.bed",
    },
    "mPanPan1": {
        "name": "Bonobo",
        "assembly": "mPanPan1.dip.cur.20231122",
        "genome": ROOT / "primate_genomes/mPanPan1.dip.cur.20231122.fasta.gz",
        "bed": ROOT / "primate_acro_bedfiles/mPanPan1/acro_monomers_full_mPanPan1.P.bed",
    },
    "mPanTro3": {
        "name": "Chimpanzee",
        "assembly": "mPanTro3.dip.cur.20231122",
        "genome": ROOT / "primate_genomes/mPanTro3.dip.cur.20231122.fasta.gz",
        "bed": ROOT / "primate_acro_bedfiles/mPanTro3/acro_monomers_full_mPanTro3.P.bed",
    },
    "mGorGor1": {
        "name": "Gorilla",
        "assembly": "mGorGor1.dip.cur.20231122",
        "genome": ROOT / "primate_genomes/mGorGor1.dip.cur.20231122.fasta.gz",
        "bed": ROOT / "primate_acro_bedfiles/mGorGor1/acro_monomers_full_mGorGor1.P.bed",
    },
    "mPonPyg2": {
        "name": "Bornean orangutan",
        "assembly": "mPonPyg2.dip.cur.20231122",
        "genome": ROOT / "primate_genomes/mPonPyg2.dip.cur.20231122.fasta.gz",
        "bed": ROOT / "primate_acro_bedfiles/mPonPyg2/acro_monomers_full_mPonPyg2.P.bed",
    },
    "mPonAbe1": {
        "name": "Sumatran orangutan",
        "assembly": "mPonAbe1.dip.cur.20231205",
        "genome": ROOT / "primate_genomes/mPonAbe1.dip.cur.20231205.fasta.gz",
        "bed": ROOT / "primate_acro_bedfiles/mPonAbe1/acro_monomers_full_mPonAbe1.P.bed",
    },
}

TYPE_ORDER = ["canonical_ACRO1", "internal_CCC", "boundary_terminal"]
ALL_LABELS = TYPE_ORDER + ["unclassified_144", "other_spacing"]
TYPE_COLORS = {
    "canonical_ACRO1": "#4E79A7",
    "internal_CCC": "#E41A1C",
    "boundary_terminal": "#FF9D00",
    "unclassified_144": "#8C8C8C",
    "other_spacing": "#CFCFCF",
}
BED_RGB = {
    "canonical_ACRO1": "78,121,167",
    "internal_CCC": "228,26,28",
    "boundary_terminal": "255,157,0",
    "unclassified_144": "140,140,140",
    "other_spacing": "207,207,207",
}
NEAR_MIN = 139
NEAR_MAX = 158

DISPLAY_NAMES = {
    "canonical_ACRO1": "canonical (147 bp)",
    "internal_CCC": "CCC-missing (144 bp; deletion at ACRO1 positions 6-8)",
    "boundary_terminal": "terminal-missing (144 bp; deletion at ACRO1 positions 145-147)",
    "unclassified_144": "unclassified 144 bp",
    "other_spacing": "other 139-158 bp spacing",
}


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def consensus(sequences: list[str]) -> str:
    if not sequences:
        raise RuntimeError("Cannot form a consensus from zero sequences")
    length = len(sequences[0])
    if any(len(sequence) != length for sequence in sequences):
        raise RuntimeError("Consensus sequences have unequal lengths")
    result = []
    for position in range(length):
        counts = Counter(sequence[position] for sequence in sequences if sequence[position] in "ACGT")
        result.append(counts.most_common(1)[0][0] if counts else "N")
    return "".join(result)


def best_circular_distance(sequence: str, target: str) -> tuple[int, int]:
    best_distance = math.inf
    best_rotation = 0
    for rotation in range(len(sequence)):
        rotated = sequence[rotation:] + sequence[:rotation]
        distance = edlib.align(rotated, target, mode="NW", task="distance")["editDistance"]
        if distance < best_distance:
            best_distance = distance
            best_rotation = rotation
    return int(best_distance), best_rotation


def gap_signature(sequence: str, reference: str) -> dict[str, object]:
    alignment = edlib.align(sequence, reference, mode="NW", task="path")
    nice = edlib.getNiceAlignment(alignment, sequence, reference)
    reference_position = 0
    query_gaps = []
    target_gaps = 0
    for query_base, target_base in zip(nice["query_aligned"], nice["target_aligned"]):
        if target_base == "-":
            target_gaps += query_base != "-"
        else:
            if query_base == "-":
                query_gaps.append(reference_position)
            reference_position += 1
    return {
        "edit_distance_to_species_147_consensus": int(alignment["editDistance"]),
        "query_gap_positions_zero_based_on_species_147": ",".join(map(str, query_gaps)),
        "query_gap_count": len(query_gaps),
        "target_gap_count": target_gaps,
        "has_internal_CCC_signature": {5, 6, 7}.issubset(query_gaps),
        "has_boundary_terminal_signature": {144, 145, 146}.issubset(query_gaps),
        "simple_internal_CCC": target_gaps == 0 and tuple(query_gaps) == (5, 6, 7),
        "simple_boundary_terminal": target_gaps == 0 and tuple(query_gaps) == (144, 145, 146),
    }


def load_array_mapping() -> dict[tuple[str, str], dict[str, object]]:
    frame = pd.read_csv(ARRAY_MAPPING, sep="\t")
    frame["array_id"] = frame["input_index"].map(lambda value: f"ARRAY_{int(value):06d}")
    frame["array_label"] = frame["original_bed_name"].str.split("|").str[0]
    return {
        (str(row.species_id), str(row.array_label)): row._asdict()
        for row in frame.itertuples(index=False)
    }


def load_acre_index() -> tuple[pd.DataFrame, dict[tuple[str, str], tuple[list[int], list[dict[str, object]]]]]:
    metadata = pd.read_csv(ACRE_METADATA, sep="\t")
    result: dict[tuple[str, str], tuple[list[int], list[dict[str, object]]]] = {}
    for (species, chrom), group in metadata.groupby(["species_id", "original_chrom"], sort=False):
        group = group.sort_values("original_acre_start")
        records = group.to_dict(orient="records")
        result[(str(species), str(chrom))] = (
            [int(record["original_acre_start"]) for record in records], records
        )
    return metadata, result


def locate_acre(
    species: str,
    chrom: str,
    midpoint: int,
    index: dict[tuple[str, str], tuple[list[int], list[dict[str, object]]]],
) -> dict[str, object] | None:
    key = (species, chrom)
    if key not in index:
        return None
    starts, records = index[key]
    position = bisect.bisect_right(starts, midpoint) - 1
    if position < 0:
        return None
    record = records[position]
    if midpoint < int(record["original_acre_end"]):
        return record
    return None


def extract_species_links(
    species: str,
    config: dict[str, object],
    array_mapping: dict[tuple[str, str], dict[str, object]],
    acre_index: dict[tuple[str, str], tuple[list[int], list[dict[str, object]]]],
) -> tuple[pd.DataFrame, dict[str, int]]:
    columns = ["chrom", "start", "end", "name", "score", "strand"]
    hits = pd.read_csv(config["bed"], sep="\t", names=columns, usecols=range(6), header=None)
    hits["array_core"] = hits["name"].str.replace(r"_monomer_[0-9]+.*$", "", regex=True)
    hits["array_label"] = hits["array_core"].str.split("|").str[0]
    hits = hits.sort_values(["array_core", "strand", "start", "end"]).reset_index(drop=True)
    genome = pysam.FastaFile(str(config["genome"]))
    rows: list[dict[str, object]] = []
    total_adjacent = 0
    skipped_non_near = 0
    groups = 0
    for (array_core, strand), group in hits.groupby(["array_core", "strand"], sort=False):
        groups += 1
        group = group.sort_values(["start", "end"]).reset_index(drop=True)
        array_label = str(group.iloc[0]["array_label"])
        mapping = array_mapping.get((species, array_label))
        if mapping is None:
            raise RuntimeError(f"No array mapping for {species}/{array_label}")
        for position in range(len(group) - 1):
            total_adjacent += 1
            left = group.iloc[position]
            right = group.iloc[position + 1]
            spacing = int(right.start - left.start)
            if not (NEAR_MIN <= spacing <= NEAR_MAX):
                skipped_non_near += 1
                continue
            start = int(left.start)
            end = int(right.start)
            sequence = genome.fetch(str(left.chrom), start, end).upper()
            if str(strand) == "-":
                sequence = reverse_complement(sequence)
            if len(sequence) != spacing:
                raise RuntimeError(f"Fetched length mismatch: {species} {left.chrom}:{start}-{end}")
            midpoint = (start + end) // 2
            acre = locate_acre(species, str(left.chrom), midpoint, acre_index)
            rows.append(
                {
                    "species_id": species,
                    "species_name": config["name"],
                    "assembly": config["assembly"],
                    "chrom": str(left.chrom),
                    "start": start,
                    "end": end,
                    "spacing": spacing,
                    "strand": str(strand),
                    "source_hit_name": str(left["name"]),
                    "next_hit_name": str(right["name"]),
                    "array_label": array_label,
                    "array_id": str(mapping["array_id"]),
                    "position_in_hit_chain": position + 1,
                    "acre_id": str(acre["acre_id"]) if acre else "",
                    "unit_in_array": int(acre["unit_in_array"]) if acre else np.nan,
                    "cluster8": int(acre["cluster8"]) if acre else np.nan,
                    "synteny_group": str(acre["synteny_group"]) if acre else "unmapped",
                    "nor_status": str(acre["nor_status"]) if acre else "unmapped",
                    "sequence": sequence,
                }
            )
    genome.close()
    return pd.DataFrame(rows), {
        "input_hits": len(hits),
        "strand_array_groups": groups,
        "adjacent_hit_pairs": total_adjacent,
        "near_spacing_links": len(rows),
        "non_near_adjacent_pairs": skipped_non_near,
    }


def classify_species(
    frame: pd.DataFrame,
    species_reference: str,
    gorilla_motifs: dict[str, str],
) -> pd.DataFrame:
    rows = []
    for row in frame.itertuples(index=False):
        sequence = str(row.sequence)
        result = {
            "subfamily_label": "other_spacing",
            "classification_method": "spacing_outside_144_147",
            "pattern_class": "other_spacing",
            "edit_distance_to_species_147_consensus": np.nan,
            "query_gap_positions_zero_based_on_species_147": "",
            "query_gap_count": np.nan,
            "target_gap_count": np.nan,
            "has_internal_CCC_signature": False,
            "has_boundary_terminal_signature": False,
            "edit_distance_to_gorilla_internal_motif": np.nan,
            "edit_distance_to_gorilla_boundary_motif": np.nan,
            "nearest_144_gorilla_motif": "",
            "nearest_144_motif_distance_margin": np.nan,
        }
        if row.spacing == 147:
            result.update(
                {
                    "subfamily_label": "canonical_ACRO1",
                    "classification_method": "exact_147_spacing",
                    "pattern_class": "canonical_147",
                    "edit_distance_to_species_147_consensus": int(
                        edlib.align(sequence, species_reference, mode="NW", task="distance")["editDistance"]
                    ),
                }
            )
        elif row.spacing == 144:
            signature = gap_signature(sequence, species_reference)
            internal_distance = int(
                edlib.align(sequence, gorilla_motifs["internal_CCC"], mode="NW", task="distance")["editDistance"]
            )
            boundary_distance = int(
                edlib.align(sequence, gorilla_motifs["boundary_terminal"], mode="NW", task="distance")["editDistance"]
            )
            if internal_distance < boundary_distance:
                nearest, margin = "internal_CCC", boundary_distance - internal_distance
            elif boundary_distance < internal_distance:
                nearest, margin = "boundary_terminal", internal_distance - boundary_distance
            else:
                nearest, margin = "tie", 0
            if signature["has_internal_CCC_signature"] and signature["has_boundary_terminal_signature"]:
                label = "unclassified_144"
                method = "both_gap_signatures"
                pattern = "both_signatures"
            elif signature["has_internal_CCC_signature"]:
                label = "internal_CCC"
                method = "gap_signature_6_8"
                pattern = "simple_internal_CCC" if signature["simple_internal_CCC"] else "complex_internal_CCC"
            elif signature["has_boundary_terminal_signature"]:
                label = "boundary_terminal"
                method = "gap_signature_145_147"
                pattern = "simple_boundary_terminal" if signature["simple_boundary_terminal"] else "complex_boundary_terminal"
            else:
                label = "unclassified_144"
                method = "no_defining_gap_signature"
                pattern = "other_144_deletion_pattern"
            result.update(signature)
            result.update(
                {
                    "subfamily_label": label,
                    "classification_method": method,
                    "pattern_class": pattern,
                    "edit_distance_to_gorilla_internal_motif": internal_distance,
                    "edit_distance_to_gorilla_boundary_motif": boundary_distance,
                    "nearest_144_gorilla_motif": nearest,
                    "nearest_144_motif_distance_margin": margin,
                }
            )
        rows.append(result)
    classifications = pd.DataFrame(rows)
    result = pd.concat([frame.reset_index(drop=True), classifications], axis=1)
    result["valid_base_fraction"] = result["sequence"].map(
        lambda sequence: sum(base in "ACGT" for base in sequence) / len(sequence)
    )
    return result


def write_bed(frame: pd.DataFrame, path: Path) -> None:
    bed = pd.DataFrame(
        {
            "chrom": frame["chrom"],
            "start": frame["start"].astype(int),
            "end": frame["end"].astype(int),
            "name": [
                f"{species}|{label}|{acre or 'no_ACRE'}|{array_id}|unit{position}"
                for species, label, acre, array_id, position in zip(
                    frame["species_id"], frame["subfamily_label"], frame["acre_id"],
                    frame["array_id"], frame["position_in_hit_chain"],
                )
            ],
            "score": 0,
            "strand": frame["strand"],
            "thickStart": frame["start"].astype(int),
            "thickEnd": frame["end"].astype(int),
            "itemRgb": frame["subfamily_label"].map(BED_RGB),
        }
    )
    bed.to_csv(path, sep="\t", header=False, index=False)


def label_count_summary(frame: pd.DataFrame, groups: list[str]) -> pd.DataFrame:
    counts = (
        frame.groupby(groups + ["subfamily_label"], observed=True)
        .size().rename("monomers").reset_index()
        .pivot(index=groups, columns="subfamily_label", values="monomers")
        .fillna(0).astype(int).reset_index()
    )
    counts.columns.name = None
    for label in ALL_LABELS:
        if label not in counts:
            counts[label] = 0
    counts["three_subfamily_total"] = counts[TYPE_ORDER].sum(axis=1)
    counts["exact_144_total"] = counts[["internal_CCC", "boundary_terminal", "unclassified_144"]].sum(axis=1)
    counts["strict_144_classification_rate"] = np.where(
        counts["exact_144_total"] > 0,
        (counts["internal_CCC"] + counts["boundary_terminal"]) / counts["exact_144_total"],
        np.nan,
    )
    for label in TYPE_ORDER:
        counts[f"fraction_{label}_among_three"] = np.where(
            counts["three_subfamily_total"] > 0,
            counts[label] / counts["three_subfamily_total"], np.nan,
        )
    counts["fraction_internal_among_classified_144"] = np.where(
        (counts["internal_CCC"] + counts["boundary_terminal"]) > 0,
        counts["internal_CCC"] / (counts["internal_CCC"] + counts["boundary_terminal"]),
        np.nan,
    )
    return counts


def build_acre_summary(frame: pd.DataFrame, acre_metadata: pd.DataFrame) -> pd.DataFrame:
    mapped = frame[frame["acre_id"] != ""]
    counts = label_count_summary(mapped, ["acre_id"])
    result = acre_metadata.merge(counts, on="acre_id", how="left", validate="one_to_one")
    numeric = ALL_LABELS + ["three_subfamily_total", "exact_144_total"]
    result[numeric] = result[numeric].fillna(0).astype(int)
    result["strict_144_classification_rate"] = np.where(
        result["exact_144_total"] > 0,
        (result["internal_CCC"] + result["boundary_terminal"]) / result["exact_144_total"],
        np.nan,
    )
    for label in TYPE_ORDER:
        result[f"fraction_{label}_among_three"] = np.where(
            result["three_subfamily_total"] > 0,
            result[label] / result["three_subfamily_total"], np.nan,
        )
    result["fraction_internal_among_classified_144"] = np.where(
        (result["internal_CCC"] + result["boundary_terminal"]) > 0,
        result["internal_CCC"] / (result["internal_CCC"] + result["boundary_terminal"]),
        np.nan,
    )
    return result


def build_array_summary(frame: pd.DataFrame, array_mapping: pd.DataFrame) -> pd.DataFrame:
    counts = label_count_summary(frame, ["array_id"])
    metadata = array_mapping.copy()
    metadata["array_id"] = metadata["input_index"].map(lambda value: f"ARRAY_{int(value):06d}")
    result = metadata.merge(counts, on="array_id", how="left", validate="one_to_one")
    numeric = ALL_LABELS + ["three_subfamily_total", "exact_144_total"]
    result[numeric] = result[numeric].fillna(0).astype(int)
    return result


def gorilla_validation(frame: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, object]]:
    observed = frame[(frame["species_id"] == "mGorGor1") & (frame["spacing"] == 144)].copy()
    expected = pd.read_csv(GORILLA_CLASSIFICATIONS, sep="\t")
    expected["expected_label"] = expected["broad_signature"].map(
        {
            "internal_CCC_signature": "internal_CCC",
            "boundary_terminal_signature": "boundary_terminal",
            "other": "unclassified_144",
            "both_signatures": "unclassified_144",
        }
    )
    merged = observed.merge(
        expected[["array_id", "start", "end", "expected_label"]],
        on=["array_id", "start", "end"], how="inner", validate="one_to_one",
    )
    confusion = pd.crosstab(merged["expected_label"], merged["subfamily_label"]).reset_index()
    return confusion, {
        "expected_rows": len(expected),
        "matched_rows": len(merged),
        "exact_label_matches": int((merged["expected_label"] == merged["subfamily_label"]).sum()),
        "accuracy": float((merged["expected_label"] == merged["subfamily_label"]).mean()),
    }


def plot_results(species_summary: pd.DataFrame, cluster_summary: pd.DataFrame) -> None:
    names = [SPECIES_CONFIG[species]["name"] for species in SPECIES_ORDER]
    summary = species_summary.set_index("species_id").reindex(SPECIES_ORDER)
    fig, axes = plt.subplots(1, 2, figsize=(17, 7.5), constrained_layout=True)

    ax = axes[0]
    bottom = np.zeros(len(summary))
    for label in TYPE_ORDER:
        values = summary[f"fraction_{label}_among_three"].fillna(0).to_numpy() * 100
        ax.bar(
            names, values, bottom=bottom, color=TYPE_COLORS[label], edgecolor="white",
            label=DISPLAY_NAMES[label].split(" (")[0],
        )
        bottom += values
    ax.set_ylabel("Percent of classified three-subfamily monomers")
    ax.set_ylim(0, 100)
    ax.set_title("A. ACRO1 subfamily composition")
    ax.tick_params(axis="x", rotation=30)
    ax.legend(frameon=False, fontsize=8, bbox_to_anchor=(0.5, -0.30), loc="upper center", ncol=3)
    ax.spines[["top", "right"]].set_visible(False)

    ax = axes[1]
    classified = summary["internal_CCC"] + summary["boundary_terminal"]
    internal = np.where(classified > 0, summary["internal_CCC"] / classified * 100, 0)
    boundary = np.where(classified > 0, summary["boundary_terminal"] / classified * 100, 0)
    unclassified = np.where(summary["exact_144_total"] > 0, summary["unclassified_144"] / summary["exact_144_total"] * 100, 0)
    ax.bar(names, internal, color=TYPE_COLORS["internal_CCC"], label="CCC-missing among classified 144")
    ax.bar(names, boundary, bottom=internal, color=TYPE_COLORS["boundary_terminal"], label="Terminal-missing among classified 144")
    for index, value in enumerate(unclassified):
        ax.text(index, 102, f"{value:.1f}% unclassified", ha="center", va="bottom", fontsize=7, rotation=25)
    ax.set_ylabel("Percent of classified 144-bp monomers")
    ax.set_ylim(0, 116)
    ax.set_title("B. Split of the two 144-bp forms")
    ax.tick_params(axis="x", rotation=30)
    ax.legend(frameon=False, fontsize=8, bbox_to_anchor=(1.01, 1), loc="upper left")
    ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Three ACRO1 subfamilies across six primates", fontsize=16, fontweight="bold")
    fig.savefig(OUT / "six_primate_subfamily_composition.png", dpi=220, facecolor="white")
    fig.savefig(OUT / "six_primate_subfamily_composition.pdf", facecolor="white")
    plt.close(fig)

    matrix = cluster_summary.pivot(index="species_id", columns="cluster8", values="fraction_internal_among_classified_144")
    matrix = matrix.reindex(index=SPECIES_ORDER, columns=range(1, 9))
    fig, ax = plt.subplots(figsize=(10, 5.8), constrained_layout=True)
    image = ax.imshow(matrix.to_numpy() * 100, aspect="auto", cmap="RdYlBu_r", vmin=0, vmax=100)
    ax.set_xticks(range(8), [f"C{i}" for i in range(1, 9)])
    ax.set_yticks(range(6), names)
    ax.set_xlabel("Joint six-primate ACRE cluster")
    ax.set_title("Internal-CCC share among classified 144-bp ACRO1 monomers")
    for row in range(matrix.shape[0]):
        for column in range(matrix.shape[1]):
            value = matrix.iloc[row, column]
            if pd.notna(value):
                ax.text(column, row, f"{value * 100:.0f}%", ha="center", va="center", fontsize=8,
                        color="white" if value < 0.18 or value > 0.82 else "black")
    colorbar = fig.colorbar(image, ax=ax, shrink=0.86)
    colorbar.set_label("Internal-CCC (%)")
    fig.savefig(OUT / "species_cluster_internal_fraction.png", dpi=220, facecolor="white")
    fig.savefig(OUT / "species_cluster_internal_fraction.pdf", facecolor="white")
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    motif_summary = pd.read_csv(MOTIF_SUMMARY, sep="\t").set_index("motif_type")
    gorilla_motifs = motif_summary["native_consensus"].to_dict()
    array_mapping_dict = load_array_mapping()
    array_mapping_frame = pd.read_csv(ARRAY_MAPPING, sep="\t")
    acre_metadata, acre_index = load_acre_index()

    extracted_frames = []
    coverage_rows = []
    for species in SPECIES_ORDER:
        frame, coverage = extract_species_links(
            species, SPECIES_CONFIG[species], array_mapping_dict, acre_index
        )
        extracted_frames.append(frame)
        coverage_rows.append({"species_id": species, "species_name": SPECIES_CONFIG[species]["name"], **coverage})
    extracted = pd.concat(extracted_frames, ignore_index=True)

    species_consensuses = {}
    consensus_rows = []
    classified_frames = []
    for species in SPECIES_ORDER:
        frame = extracted[extracted["species_id"] == species].copy()
        sequences_147 = frame.loc[frame["spacing"] == 147, "sequence"].tolist()
        species_reference = consensus(sequences_147)
        species_consensuses[species] = species_reference
        direct_distance = int(
            edlib.align(species_reference, gorilla_motifs["canonical_ACRO1"], mode="NW", task="distance")["editDistance"]
        )
        circular_distance, rotation = best_circular_distance(
            species_reference, gorilla_motifs["canonical_ACRO1"]
        )
        consensus_rows.append(
            {
                "species_id": species,
                "species_name": SPECIES_CONFIG[species]["name"],
                "canonical_147_sequences": len(sequences_147),
                "direct_edit_distance_to_gorilla_canonical": direct_distance,
                "best_circular_edit_distance_to_gorilla_canonical": circular_distance,
                "best_circular_rotation": rotation,
                "canonical_147_consensus": species_reference,
            }
        )
        classified_frames.append(classify_species(frame, species_reference, gorilla_motifs))
    classified = pd.concat(classified_frames, ignore_index=True)

    with (OUT / "six_primate_species_canonical_consensus.fasta").open("w") as handle:
        for species in SPECIES_ORDER:
            handle.write(f">{species}_canonical_ACRO1_147\n{species_consensuses[species]}\n")
    pd.DataFrame(consensus_rows).to_csv(
        OUT / "species_canonical_consensus_summary.tsv", sep="\t", index=False
    )
    pd.DataFrame(coverage_rows).to_csv(OUT / "input_coverage_summary.tsv", sep="\t", index=False)

    classified.to_csv(
        OUT / "all_primate_acro1_subfamily_labels.tsv.gz", sep="\t", index=False,
        compression={"method": "gzip", "compresslevel": 6},
    )
    for species in SPECIES_ORDER:
        subset = classified[classified["species_id"] == species].copy()
        subset.to_csv(
            OUT / f"{species}_acro1_subfamily_labels.tsv.gz", sep="\t", index=False,
            compression={"method": "gzip", "compresslevel": 6},
        )
        write_bed(subset, OUT / f"{species}_acro1_subfamilies.bed")

    species_summary = label_count_summary(classified, ["species_id", "species_name"])
    species_summary.to_csv(OUT / "species_subfamily_summary.tsv", sep="\t", index=False)
    label_count_summary(
        classified[classified["acre_id"] != ""],
        ["species_id", "species_name", "nor_status"],
    ).to_csv(OUT / "species_nor_subfamily_summary.tsv", sep="\t", index=False)
    label_count_summary(
        classified[classified["acre_id"] != ""],
        ["species_id", "species_name", "synteny_group"],
    ).to_csv(OUT / "species_synteny_subfamily_summary.tsv", sep="\t", index=False)
    cluster_summary = label_count_summary(
        classified[classified["acre_id"] != ""],
        ["species_id", "species_name", "cluster8"],
    )
    cluster_summary.to_csv(OUT / "species_cluster_subfamily_summary.tsv", sep="\t", index=False)

    acre_summary = build_acre_summary(classified, acre_metadata)
    acre_summary.to_csv(OUT / "acre_subfamily_summary.tsv.gz", sep="\t", index=False, compression="gzip")
    array_summary = build_array_summary(classified, array_mapping_frame)
    array_summary.to_csv(OUT / "array_subfamily_summary.tsv", sep="\t", index=False)

    confusion, validation = gorilla_validation(classified)
    confusion.to_csv(OUT / "gorilla_classifier_recovery_confusion.tsv", sep="\t", index=False)
    plot_results(species_summary, cluster_summary)

    summary = {
        "species": SPECIES_ORDER,
        "near_spacing_links": int(len(classified)),
        "three_subfamily_classified_links": int(classified["subfamily_label"].isin(TYPE_ORDER).sum()),
        "mapped_to_acre": int((classified["acre_id"] != "").sum()),
        "labels": {str(key): int(value) for key, value in classified["subfamily_label"].value_counts().items()},
        "gorilla_recovery_validation": validation,
        "classification_definition": {
            "canonical_ACRO1": "exact 147-bp start-to-start spacing",
            "internal_CCC": "exact 144-bp spacing with query gaps at one-based ACRO1 positions 6-8",
            "boundary_terminal": "exact 144-bp spacing with query gaps at one-based ACRO1 positions 145-147",
            "unclassified_144": "exact 144-bp spacing lacking exactly one defining broad gap signature",
            "other_spacing": f"near spacing {NEAR_MIN}-{NEAR_MAX} bp other than 144/147",
        },
    }
    with (OUT / "analysis_summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
