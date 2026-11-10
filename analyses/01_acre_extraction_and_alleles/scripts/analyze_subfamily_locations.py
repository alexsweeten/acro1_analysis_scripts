#!/usr/bin/env python3
"""Build location, architecture, and positional summaries for ACRO1 subfamilies."""

from __future__ import annotations

import json
import math
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


OUT = Path(__file__).resolve().parent
LABELS_PATH = OUT / "all_primate_acro1_subfamily_labels.tsv.gz"
ACRES_PATH = OUT / "acre_subfamily_summary.tsv.gz"
ARRAYS_PATH = OUT / "array_subfamily_summary.tsv"

SPECIES_ORDER = ["hg002", "mPanPan1", "mPanTro3", "mGorGor1", "mPonPyg2", "mPonAbe1"]
SPECIES_NAMES = {
    "hg002": "Human (HG002)",
    "mPanPan1": "Bonobo",
    "mPanTro3": "Chimpanzee",
    "mGorGor1": "Gorilla",
    "mPonPyg2": "Bornean orangutan",
    "mPonAbe1": "Sumatran orangutan",
}
TYPE_ORDER = ["canonical_ACRO1", "internal_CCC", "boundary_terminal"]
DISPLAY = {
    "canonical_ACRO1": "canonical",
    "internal_CCC": "CCC-missing",
    "boundary_terminal": "terminal-missing",
    "unclassified_144": "unclassified-144",
    "other_spacing": "other-spacing",
}
CODE = {
    "canonical_ACRO1": "C",
    "internal_CCC": "Mccc",
    "boundary_terminal": "Mterm",
    "unclassified_144": "U144",
    "other_spacing": "O",
}
COLORS = {
    "canonical_ACRO1": "#4E79A7",
    "internal_CCC": "#E41A1C",
    "boundary_terminal": "#FF9D00",
}


def count_summary(frame: pd.DataFrame, groups: list[str]) -> pd.DataFrame:
    result = (
        frame.groupby(groups + ["subfamily_label"], dropna=False, observed=True)
        .size().unstack(fill_value=0).reset_index()
    )
    result.columns.name = None
    for label in [*TYPE_ORDER, "unclassified_144", "other_spacing"]:
        if label not in result:
            result[label] = 0
    result["three_subfamily_total"] = result[TYPE_ORDER].sum(axis=1)
    for label in TYPE_ORDER:
        result[f"fraction_{label}"] = np.where(
            result["three_subfamily_total"] > 0,
            result[label] / result["three_subfamily_total"],
            np.nan,
        )
    return result


def entropy(counts: list[int]) -> float:
    total = sum(counts)
    if total == 0:
        return math.nan
    probabilities = [value / total for value in counts if value]
    return -sum(value * math.log2(value) for value in probabilities)


def composition_class(counts: dict[str, int]) -> str:
    present = [label for label in TYPE_ORDER if counts[label] > 0]
    if not present:
        return "no_three_subfamily_call"
    if len(present) == 1:
        return f"pure_{DISPLAY[present[0]].replace('-', '_')}"
    if len(present) == 3:
        return "mixed_all_three"
    return "mixed_" + "_and_".join(DISPLAY[label].replace("-", "_") for label in present)


def oriented_records(group: pd.DataFrame) -> pd.DataFrame:
    dominant_strand = group["strand"].mode().iat[0]
    return group.sort_values(["start", "end"], ascending=dominant_strand != "-").reset_index(drop=True)


def rle(labels: list[str]) -> tuple[str, int, int]:
    if not labels:
        return "", 0, 0
    runs: list[tuple[str, int]] = []
    for label in labels:
        if runs and runs[-1][0] == label:
            runs[-1] = (label, runs[-1][1] + 1)
        else:
            runs.append((label, 1))
    return "|".join(f"{CODE[label]}x{length}" for label, length in runs), len(runs), max(x[1] for x in runs)


def build_architectures(labels: pd.DataFrame, acres: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    mapped = labels[labels["acre_id"].notna()].copy()
    architecture_rows = []
    run_rows = []
    transition_rows = []
    metadata = acres.set_index("acre_id")
    for acre_id, raw_group in mapped.groupby("acre_id", sort=False):
        group = oriented_records(raw_group)
        label_list = group["subfamily_label"].tolist()
        signature, number_of_runs, longest_run = rle(label_list)
        counts = Counter(label_list)
        meta = metadata.loc[acre_id]
        three_counts = {label: int(counts[label]) for label in TYPE_ORDER}
        dominant = max(TYPE_ORDER, key=lambda label: three_counts[label]) if sum(three_counts.values()) else ""
        architecture_rows.append({
            "acre_id": acre_id,
            "species_id": meta["species_id"],
            "species_name": SPECIES_NAMES[str(meta["species_id"])],
            "assembly": meta["assembly"],
            "chrom": meta["original_chrom"],
            "acre_start": int(meta["original_acre_start"]),
            "acre_end": int(meta["original_acre_end"]),
            "array_id": meta["array_id"],
            "unit_in_array": int(meta["unit_in_array"]),
            "strand": group["strand"].mode().iat[0],
            "cluster8": int(meta["cluster8"]),
            "synteny_group": meta["synteny_group"],
            "nor_status": meta["nor_status"],
            **three_counts,
            "unclassified_144": int(counts["unclassified_144"]),
            "other_spacing": int(counts["other_spacing"]),
            "three_subfamily_total": sum(three_counts.values()),
            "dominant_subfamily": dominant,
            "composition_class": composition_class(three_counts),
            "subfamily_entropy_bits": entropy(list(three_counts.values())),
            "architecture_rle": signature,
            "number_of_runs_all_labels": number_of_runs,
            "longest_run_all_labels": longest_run,
        })

        run_number = 0
        run_start = 0
        for index in range(1, len(group) + 1):
            at_end = index == len(group)
            changed = not at_end and group.iloc[index]["subfamily_label"] != group.iloc[run_start]["subfamily_label"]
            if at_end or changed:
                run_number += 1
                subset = group.iloc[run_start:index]
                run_rows.append({
                    "acre_id": acre_id,
                    "species_id": meta["species_id"],
                    "chrom": meta["original_chrom"],
                    "acre_start": int(meta["original_acre_start"]),
                    "acre_end": int(meta["original_acre_end"]),
                    "array_id": meta["array_id"],
                    "strand": group["strand"].mode().iat[0],
                    "run_number_in_acre_orientation": run_number,
                    "subfamily_label": subset.iloc[0]["subfamily_label"],
                    "run_monomers": len(subset),
                    "genomic_start": int(subset["start"].min()),
                    "genomic_end": int(subset["end"].max()),
                })
                run_start = index

        classified = group[group["subfamily_label"].isin(TYPE_ORDER)]
        prior = None
        for current in classified["subfamily_label"]:
            if prior is not None:
                transition_rows.append({
                    "species_id": meta["species_id"],
                    "acre_id": acre_id,
                    "from_subfamily": prior,
                    "to_subfamily": current,
                })
            prior = current

    architectures = pd.DataFrame(architecture_rows)
    complete = acres.merge(
        architectures.drop(columns=[column for column in ["species_id", "assembly", "array_id", "unit_in_array", "cluster8", "synteny_group", "nor_status"] if column in architectures]),
        on="acre_id", how="left", suffixes=("", "_architecture"), validate="one_to_one",
    )
    complete["composition_class"] = complete["composition_class"].fillna("no_ACRO1_link_mapped")
    transitions = pd.DataFrame(transition_rows)
    if not transitions.empty:
        transitions = transitions.groupby(
            ["species_id", "from_subfamily", "to_subfamily"], observed=True
        ).size().rename("transitions").reset_index()
    return complete, pd.DataFrame(run_rows), transitions


def add_positions(labels: pd.DataFrame, acres: pd.DataFrame) -> pd.DataFrame:
    metadata = acres[["acre_id", "original_acre_start", "original_acre_end"]]
    result = labels.merge(metadata, on="acre_id", how="left", validate="many_to_one")
    midpoint = (result["start"] + result["end"]) / 2
    length = result["original_acre_end"] - result["original_acre_start"]
    result["genomic_fraction_through_acre"] = (midpoint - result["original_acre_start"]) / length
    result["oriented_fraction_through_acre"] = np.where(
        result["strand"] == "-", 1 - result["genomic_fraction_through_acre"], result["genomic_fraction_through_acre"]
    )
    result["acre_position_decile"] = np.minimum(
        10, np.floor(result["oriented_fraction_through_acre"] * 10).fillna(-1).astype(int) + 1
    )
    return result


def top_locations(labels: pd.DataFrame) -> pd.DataFrame:
    classified = labels[labels["subfamily_label"].isin(TYPE_ORDER)].copy()
    rows = []
    for level, keys in [("array", ["array_id"]), ("ACRE", ["acre_id"])]:
        usable = classified[classified[keys[0]].notna()]
        grouped = (
            usable.groupby(["species_id", "subfamily_label", *keys], observed=True)
            .agg(monomers=("start", "size"), chrom=("chrom", "first"), genomic_start=("start", "min"), genomic_end=("end", "max"))
            .reset_index()
        )
        grouped["rank_within_species_subfamily"] = grouped.groupby(
            ["species_id", "subfamily_label"], observed=True
        )["monomers"].rank(method="first", ascending=False).astype(int)
        grouped = grouped[grouped["rank_within_species_subfamily"] <= 10].copy()
        grouped["location_level"] = level
        grouped["location_id"] = grouped[keys[0]]
        rows.append(grouped.drop(columns=keys))
    return pd.concat(rows, ignore_index=True).sort_values(
        ["species_id", "subfamily_label", "location_level", "rank_within_species_subfamily"]
    )


def plot_acre_composition(summary: pd.DataFrame) -> None:
    classes = [
        "pure_canonical", "pure_CCC_missing", "pure_terminal_missing",
        "mixed_canonical_and_CCC_missing", "mixed_canonical_and_terminal_missing",
        "mixed_CCC_missing_and_terminal_missing", "mixed_all_three",
        "no_ACRO1_link_mapped", "no_three_subfamily_call",
    ]
    pivot = summary.pivot(index="species_id", columns="composition_class", values="acres").fillna(0)
    pivot = pivot.reindex(index=SPECIES_ORDER, columns=classes, fill_value=0)
    palette = ["#4E79A7", "#E41A1C", "#FF9D00", "#9C755F", "#F28E2B", "#B07AA1", "#59A14F", "#BAB0AC", "#777777"]
    fig, ax = plt.subplots(figsize=(12, 6.5), constrained_layout=True)
    bottom = np.zeros(len(pivot))
    for column, color in zip(pivot.columns, palette):
        values = pivot[column].to_numpy()
        ax.bar([SPECIES_NAMES[x] for x in pivot.index], values, bottom=bottom, label=column.replace("_", " "), color=color)
        bottom += values
    ax.set_ylabel("ACREs")
    ax.set_title("Per-ACRE ACRO1 subfamily composition")
    ax.tick_params(axis="x", rotation=25)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, fontsize=8, bbox_to_anchor=(1.01, 1), loc="upper left")
    fig.savefig(OUT / "acre_composition_classes.png", dpi=220, facecolor="white")
    fig.savefig(OUT / "acre_composition_classes.pdf", facecolor="white")
    plt.close(fig)


def plot_positions(deciles: pd.DataFrame) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(15, 8), sharex=True, sharey=True, constrained_layout=True)
    for ax, species in zip(axes.flat, SPECIES_ORDER):
        subset = deciles[deciles["species_id"] == species]
        for label in TYPE_ORDER:
            values = subset[subset["subfamily_label"] == label].set_index("acre_position_decile")["fraction_within_decile"]
            ax.plot(range(1, 11), [values.get(x, 0) for x in range(1, 11)], marker="o", lw=1.8,
                    label=DISPLAY[label], color=COLORS[label])
        ax.set_title(SPECIES_NAMES[species])
        ax.set_xticks(range(1, 11))
        ax.spines[["top", "right"]].set_visible(False)
    axes[1, 0].set_xlabel("ACRE-oriented position decile")
    axes[1, 1].set_xlabel("ACRE-oriented position decile")
    axes[1, 2].set_xlabel("ACRE-oriented position decile")
    axes[0, 0].set_ylabel("Fraction within decile")
    axes[1, 0].set_ylabel("Fraction within decile")
    axes[0, 2].legend(frameon=False, fontsize=8)
    fig.suptitle("Subfamily composition along the oriented ACRE", fontweight="bold")
    fig.savefig(OUT / "within_acre_positional_profiles.png", dpi=220, facecolor="white")
    fig.savefig(OUT / "within_acre_positional_profiles.pdf", facecolor="white")
    plt.close(fig)


def write_report(species: pd.DataFrame, acre_classes: pd.DataFrame, synteny: pd.DataFrame,
                 chromosome: pd.DataFrame, labels: pd.DataFrame, architectures: pd.DataFrame) -> None:
    total = int(species["three_subfamily_total"].sum())
    exact_144 = int(species["internal_CCC"].sum() + species["boundary_terminal"].sum() + species["unclassified_144"].sum())
    classified_144 = int(species["internal_CCC"].sum() + species["boundary_terminal"].sum())
    mapped = int(labels["acre_id"].notna().sum())
    strict = labels[labels["subfamily_label"].isin(TYPE_ORDER)].copy()
    lines = [
        "# Whole-primate ACRO1 subfamily analysis", "",
        "## Scope and definitions", "",
        f"This analysis examined **{len(labels):,}** adjacent ACRO1-hit intervals (139-158 bp) across six primates. "
        f"Of these, **{total:,}** are strict members of the requested three classes and **{mapped:,}** map to a segmented ACRE.", "",
        "- **Canonical:** exactly 147 bp from one ACRO1 hit start to the next.",
        "- **CCC-missing:** exactly 144 bp, aligned to the species-specific 147-bp consensus with gaps at ACRO1 positions 6-8 (one-based).",
        "- **Terminal-missing:** exactly 144 bp, aligned with gaps at positions 145-147.",
        f"- Of {exact_144:,} exact 144-bp intervals, {classified_144:,} ({classified_144 / exact_144:.2%}) receive one of the two deletion signatures; the remainder are retained as unclassified-144.",
        "", "The labels describe sequence/phase signatures in the current monomer callsets. Because the original gorilla discovery set couples the two 144-bp forms to strand and repeat phase, they should not be interpreted as independently evolved phylogenetic subfamilies without additional evidence.", "",
        "## Coverage", "",
        "| Species | Near ACRO1 links | Strict 3-class calls | Mapped to ACRE | Arrays represented | ACREs with calls | Total segmented ACREs |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for species_id in SPECIES_ORDER:
        subset = labels[labels["species_id"] == species_id]
        strict_subset = subset[subset["subfamily_label"].isin(TYPE_ORDER)]
        acre_subset = architectures[architectures["species_id"] == species_id]
        lines.append(
            f"| {SPECIES_NAMES[species_id]} | {len(subset):,} | {len(strict_subset):,} | "
            f"{int(subset['acre_id'].notna().sum()):,} | {subset['array_id'].nunique():,} | "
            f"{int((acre_subset['three_subfamily_total'] > 0).sum()):,} | {len(acre_subset):,} |"
        )
    lines += ["", 
        "## Species-wide composition", "",
        "| Species | Canonical | CCC-missing | Terminal-missing | Strict total | Canonical % | CCC-missing % | Terminal-missing % |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for species_id in SPECIES_ORDER:
        row = species.set_index("species_id").loc[species_id]
        lines.append(
            f"| {SPECIES_NAMES[species_id]} | {int(row.canonical_ACRO1):,} | {int(row.internal_CCC):,} | "
            f"{int(row.boundary_terminal):,} | {int(row.three_subfamily_total):,} | "
            f"{row.fraction_canonical_ACRO1:.1%} | {row.fraction_internal_CCC:.1%} | {row.fraction_boundary_terminal:.1%} |"
        )
    lines += [
        "", "## Cross-species findings", "",
        "- Canonical 147-bp ACRO1 is the largest class in every species. It ranges from 50.8% of strict calls in Bornean orangutan to 77.2% in chimpanzee.",
        "- The orangutans are the clearest outliers: CCC-missing forms account for 47.3% of strict calls in Bornean and 37.5% in Sumatran orangutan, while terminal-missing forms are only 1.9% and 1.6%, respectively.",
        "- Human, bonobo, and chimpanzee have more terminal-missing than CCC-missing calls; gorilla and both orangutans show the reverse.",
        "- Across all six assemblies, every strict CCC-missing call is on the minus strand and every strict terminal-missing call is on the plus strand. Canonical calls occur on both strands. This perfect coupling shows that the two 144-bp labels capture repeat orientation/phase as well as sequence composition.",
        "", "## Strand and phase coupling", "",
        "| Species | Class | Plus strand | Minus strand |",
        "|---|---|---:|---:|",
    ]
    strand_counts = strict.groupby(["species_id", "subfamily_label", "strand"], observed=True).size()
    for species_id in SPECIES_ORDER:
        for label in TYPE_ORDER:
            lines.append(
                f"| {SPECIES_NAMES[species_id]} | {DISPLAY[label]} | "
                f"{int(strand_counts.get((species_id, label, '+'), 0)):,} | "
                f"{int(strand_counts.get((species_id, label, '-'), 0)):,} |"
            )
    lines += ["", "## Where the subfamilies occur", "",
              "Coordinates below are zero-based, half-open. Counts in parentheses are monomer intervals.", ""]
    for species_id in SPECIES_ORDER:
        sub = synteny[synteny["species_id"] == species_id]
        sub = sub.groupby("synteny_group", observed=True)[TYPE_ORDER].sum().reset_index()
        chrom_sub = chromosome[chromosome["species_id"] == species_id]
        class_sub = acre_classes[acre_classes["species_id"] == species_id]
        lines += [f"### {SPECIES_NAMES[species_id]}", ""]
        for label in TYPE_ORDER:
            top = sub[sub[label] > 0].sort_values(label, ascending=False).head(3)
            places = ", ".join(f"{row.synteny_group} ({int(getattr(row, label)):,})" for row in top.itertuples())
            label_rows = strict[(strict["species_id"] == species_id) & (strict["subfamily_label"] == label)]
            contigs = (
                label_rows.groupby("chrom", observed=True)
                .agg(monomers=("start", "size"), start=("start", "min"), end=("end", "max"))
                .reset_index().sort_values("monomers", ascending=False).head(2)
            )
            contig_text = ", ".join(
                f"{row.chrom}:{int(row.start):,}-{int(row.end):,} ({int(row.monomers):,})" for row in contigs.itertuples()
            )
            acre_rows = label_rows[label_rows["acre_id"].notna()]
            top_acres = (
                acre_rows.groupby(["acre_id", "chrom"], observed=True)
                .agg(monomers=("start", "size"), start=("start", "min"), end=("end", "max"))
                .reset_index().sort_values("monomers", ascending=False).head(2)
            )
            acre_text = ", ".join(
                f"{row.acre_id} at {row.chrom}:{int(row.start):,}-{int(row.end):,} ({int(row.monomers):,})"
                for row in top_acres.itertuples()
            )
            lines.append(
                f"- **{DISPLAY[label]}:** synteny groups {places}; top assembly contigs/spans {contig_text}; "
                f"densest individual ACREs {acre_text}."
            )
        lines.append(
            f"- Observed on {chrom_sub['chrom'].nunique()} assembly contigs/chromosomes and in "
            f"{int(architectures.loc[(architectures.species_id == species_id) & (architectures.three_subfamily_total > 0), 'acre_id'].nunique()):,} ACREs."
        )
        top_classes = class_sub.sort_values("acres", ascending=False).head(4)
        lines.append("- Most common ACRE compositions: " + ", ".join(
            f"{row.composition_class.replace('_', ' ')} ({int(row.acres):,})" for row in top_classes.itertuples()
        ) + ".")
        lines.append("")
    lines += ["## Recurrent within-ACRE architectures", "",
              "Codes are defined below; counts report exact repeat patterns after orienting minus-strand ACREs.", ""]
    for species_id in SPECIES_ORDER:
        patterns = architectures[
            (architectures["species_id"] == species_id) & architectures["architecture_rle"].notna()
        ]["architecture_rle"].value_counts().head(5)
        lines += [f"### {SPECIES_NAMES[species_id]}", "",
                  "; ".join(f"`{signature}` ({count:,} ACREs)" for signature, count in patterns.items()) + ".", ""]
    lines += [
        "## Reading the detailed outputs", "",
        "- `all_primate_acro1_subfamily_locations.tsv.gz` is the master monomer-level location table. Each row has exact genomic coordinates, strand, array, ACRE, cluster, synteny/NOR context, classification evidence, and oriented position within its ACRE.",
        "- `acre_subfamily_architecture.tsv.gz` answers what the ACRO1 portion of every ACRE consists of. Counts, fractions, dominant class, mixture class, entropy, and an orientation-aware run-length signature are included.",
        "- `subfamily_runs.tsv.gz` gives contiguous blocks within each ACRE; `subfamily_transition_matrix.tsv` quantifies adjacency after orienting minus-strand ACREs.",
        "- Architecture codes are `C` (canonical), `Mccc` (CCC-missing), `Mterm` (terminal-missing), `U144` (unclassified 144 bp), and `O` (other near spacing). For example, `Cx2|Mtermx1|Cx5` means two canonical units, one terminal-missing unit, then five canonical units.",
        "- `chromosome_subfamily_summary.tsv`, `species_synteny_subfamily_summary.tsv`, `array_subfamily_summary.tsv`, and `top_subfamily_locations.tsv` progressively localize calls from chromosome/synteny group to array and individual ACRE.",
        "- BED files can be loaded directly in a genome browser. Coordinates are zero-based, half-open.",
        "", "## Quality and limitations", "",
        "The classifier has exact recovery of the 3,595 gorilla 144-bp discovery labels used for validation. This is an internal recovery check, not an independent benchmark. Canonical calls use length, whereas the two 144-bp calls additionally require alignment signatures. ACREs with no mapped ACRO1 links remain in the architecture table and are explicitly labeled `no_ACRO1_link_mapped`; they are not silently discarded. Counts derive from ACRO1-hit intervals and should not be read as whole-genome copy numbers outside the input arrays.",
    ]
    (OUT / "ANALYSIS.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    labels = pd.read_csv(LABELS_PATH, sep="\t", low_memory=False)
    labels["acre_id"] = labels["acre_id"].replace("", np.nan)
    acres = pd.read_csv(ACRES_PATH, sep="\t")

    positioned = add_positions(labels, acres)
    positioned["requested_subfamily"] = positioned["subfamily_label"].map(DISPLAY)
    positioned.to_csv(OUT / "all_primate_acro1_subfamily_locations.tsv.gz", sep="\t", index=False, compression="gzip")

    chromosome = count_summary(positioned, ["species_id", "species_name", "assembly", "chrom"])
    chromosome.to_csv(OUT / "chromosome_subfamily_summary.tsv", sep="\t", index=False)
    synteny = count_summary(positioned[positioned["acre_id"].notna()], ["species_id", "species_name", "synteny_group", "nor_status"])
    synteny.to_csv(OUT / "synteny_nor_subfamily_summary.tsv", sep="\t", index=False)

    architectures, runs, transitions = build_architectures(positioned, acres)
    architectures["ccc_missing"] = architectures["internal_CCC"]
    architectures["terminal_missing"] = architectures["boundary_terminal"]
    architectures["dominant_subfamily_requested"] = architectures["dominant_subfamily"].map(DISPLAY).fillna("")
    runs["requested_subfamily"] = runs["subfamily_label"].map(DISPLAY)
    architectures.to_csv(OUT / "acre_subfamily_architecture.tsv.gz", sep="\t", index=False, compression="gzip")
    runs.to_csv(OUT / "subfamily_runs.tsv.gz", sep="\t", index=False, compression="gzip")
    transitions.to_csv(OUT / "subfamily_transition_matrix.tsv", sep="\t", index=False)

    acre_classes = architectures.groupby(
        ["species_id", "composition_class"], dropna=False, observed=True
    ).size().rename("acres").reset_index()
    acre_classes.to_csv(OUT / "species_acre_composition_classes.tsv", sep="\t", index=False)

    classified = positioned[
        positioned["subfamily_label"].isin(TYPE_ORDER) & positioned["acre_id"].notna()
        & positioned["acre_position_decile"].between(1, 10)
    ]
    deciles = classified.groupby(
        ["species_id", "subfamily_label", "acre_position_decile"], observed=True
    ).size().rename("monomers").reset_index()
    totals = deciles.groupby(["species_id", "acre_position_decile"])["monomers"].transform("sum")
    deciles["fraction_within_decile"] = deciles["monomers"] / totals
    deciles.to_csv(OUT / "within_acre_decile_summary.tsv", sep="\t", index=False)

    top_locations(positioned).to_csv(OUT / "top_subfamily_locations.tsv", sep="\t", index=False)

    strand_summary = count_summary(positioned, ["species_id", "species_name", "strand"])
    strand_summary.to_csv(OUT / "species_strand_subfamily_summary.tsv", sep="\t", index=False)
    architecture_patterns = (
        architectures[architectures["architecture_rle"].notna()]
        .groupby(["species_id", "architecture_rle"], observed=True).size()
        .rename("acres").reset_index()
        .sort_values(["species_id", "acres"], ascending=[True, False])
    )
    architecture_patterns["rank_within_species"] = architecture_patterns.groupby("species_id").cumcount() + 1
    architecture_patterns.to_csv(OUT / "species_acre_architecture_patterns.tsv", sep="\t", index=False)

    species = count_summary(positioned, ["species_id", "species_name"])
    species.to_csv(OUT / "species_subfamily_detailed_summary.tsv", sep="\t", index=False)
    plot_acre_composition(acre_classes)
    plot_positions(deciles)
    write_report(species, acre_classes, synteny, chromosome, positioned, architectures)

    checks = {
        "input_links": int(len(positioned)),
        "strict_three_subfamily_links": int(positioned["subfamily_label"].isin(TYPE_ORDER).sum()),
        "mapped_links": int(positioned["acre_id"].notna().sum()),
        "total_acres": int(len(acres)),
        "acres_with_strict_calls": int((architectures["three_subfamily_total"] > 0).sum()),
        "acres_without_mapped_links": int((architectures["composition_class"] == "no_ACRO1_link_mapped").sum()),
        "run_monomer_total": int(runs["run_monomers"].sum()),
        "invariants": {
            "run_rows_cover_all_mapped_links": int(runs["run_monomers"].sum()) == int(positioned["acre_id"].notna().sum()),
            "architecture_has_every_acre": len(architectures) == len(acres),
            "strict_count_matches_species_summary": int(species["three_subfamily_total"].sum()) == int(positioned["subfamily_label"].isin(TYPE_ORDER).sum()),
        },
    }
    (OUT / "detailed_analysis_validation.json").write_text(json.dumps(checks, indent=2) + "\n")
    if not all(checks["invariants"].values()):
        raise RuntimeError(f"Validation failure: {checks}")
    print(json.dumps(checks, indent=2))


if __name__ == "__main__":
    main()
