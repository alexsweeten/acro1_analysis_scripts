#!/usr/bin/env python3
"""Extract every gorilla ACRE and its feature/coordinate metadata."""

from __future__ import annotations

import csv
import re
from pathlib import Path

import pandas as pd


HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / "primate_full"
NOR_POSITIVE = {"hsa21", "hsa22"}
NOR_NEGATIVE = {"hsa2a", "hsa2b", "hsa9", "hsa13", "hsa15", "hsa18", "hsaY"}


def fasta_records(path: Path):
    name: str | None = None
    parts: list[str] = []
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts).upper()
                name = line[1:].strip().split()[0]
                parts = []
            else:
                parts.append(line.strip())
    if name is not None:
        yield name, "".join(parts).upper()


def synteny(chromosome: str) -> str:
    match = re.search(r"hsa(2[ab]|[0-9]+|[XY])", chromosome)
    if not match:
        raise ValueError(f"Cannot resolve HSA synteny from {chromosome}")
    return f"hsa{match.group(1)}"


def nor_status(group: str) -> str:
    if group in NOR_POSITIVE:
        return "NOR+"
    if group in NOR_NEGATIVE:
        return "NOR-"
    if group == "hsa14":
        return "ancestral"
    return "non-acro"


def main() -> None:
    metadata = pd.read_csv(SOURCE / "six_primate_acre_metadata.tsv", sep="\t")
    metadata = metadata[metadata["species_id"] == "mGorGor1"].copy()
    if metadata["acre_id"].duplicated().any():
        raise RuntimeError("Duplicate gorilla ACRE IDs in source metadata")

    metadata["tip_id"] = metadata["acre_id"]
    metadata["synteny_group"] = metadata["original_chrom"].map(synteny)
    metadata["haplotype"] = metadata["original_chrom"].str.extract(
        r"_(mat|pat)_", expand=False
    )
    metadata["nor_status"] = metadata["synteny_group"].map(nor_status)
    # Deliberately exclude all prior joint-clustering fields from the gorilla-only
    # metadata. New cluster assignments are added only after the de novo sweep.
    prior_cluster_fields = {
        "cluster8",
        "best_k",
        "best_cluster",
        "cluster8_centroid_distance",
    }
    ordered = [
        "tip_id",
        "acre_id",
        "array_id",
        "unit_in_array",
        "species_id",
        "species_name",
        "assembly",
        "original_chrom",
        "synteny_group",
        "haplotype",
        "nor_status",
        "original_acre_start",
        "original_acre_end",
        "length",
        "estimated_period",
        "periodicity_score",
        "top_component",
        "component_coverage",
        "valid_kmers",
        "distinct_kmers",
        "exclusive_kmers",
        "exclusive_fraction",
        "synthetic_chrom",
        "start",
        "end",
        "source_name",
    ]
    ordered += [
        column
        for column in metadata.columns
        if column not in ordered and column not in prior_cluster_fields
    ]
    metadata = metadata[ordered]
    metadata.to_csv(HERE / "gorilla_tip_metadata.base.tsv", sep="\t", index=False)

    features = pd.read_csv(SOURCE / "six_primate_full.features.tsv", sep="\t")
    features = features[features["acre_id"].isin(metadata["acre_id"])].copy()
    features = metadata[["acre_id", "array_id"]].merge(
        features.drop(columns=["array_id", "cluster"]),
        on="acre_id",
        how="left",
        validate="one_to_one",
    )
    if features.isna().any().any():
        raise RuntimeError("Missing gorilla features after extraction")
    features.to_csv(HERE / "gorilla.features.tsv", sep="\t", index=False)

    wanted = set(metadata["acre_id"])
    sequences: dict[str, str] = {}
    for raw_name, sequence in fasta_records(SOURCE / "six_primate_full.acres.fasta"):
        name = raw_name.split("|", 1)[0]
        if name in wanted:
            sequences[name] = sequence
    if set(sequences) != wanted:
        missing = sorted(wanted - set(sequences))[:10]
        raise RuntimeError(f"Missing {len(wanted - set(sequences))} sequences: {missing}")

    lengths = metadata.set_index("acre_id")["length"].astype(int)
    invalid = []
    with (HERE / "gorilla_acres.fasta").open("w") as output:
        for name in metadata["acre_id"]:
            sequence = sequences[name]
            if len(sequence) != lengths[name]:
                raise RuntimeError(
                    f"Length mismatch for {name}: FASTA={len(sequence)}, metadata={lengths[name]}"
                )
            if not set(sequence) & set("ACGT"):
                invalid.append(name)
            output.write(f">{name}\n")
            for offset in range(0, len(sequence), 80):
                output.write(sequence[offset : offset + 80] + "\n")
    if invalid:
        raise RuntimeError(f"All-ambiguous gorilla sequences: {invalid[:10]}")

    chromosome_counts = (
        metadata.groupby(["synteny_group", "original_chrom", "haplotype", "nor_status"], observed=True)
        .agg(acres=("acre_id", "size"), parent_arrays=("array_id", "nunique"))
        .reset_index()
    )
    chromosome_counts.to_csv(HERE / "chromosome_counts.tsv", sep="\t", index=False)

    with (HERE / "extraction_summary.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerow(["complete_gorilla_acres", len(metadata)])
        writer.writerow(["parent_arrays", metadata["array_id"].nunique()])
        writer.writerow(["original_chromosome_contigs", metadata["original_chrom"].nunique()])
        writer.writerow(["hsa_synteny_chromosomes", metadata["synteny_group"].nunique()])
        writer.writerow(["total_sequence_bp", int(metadata["length"].sum())])
        writer.writerow(["numeric_features", len(features.columns) - 2])

    print(
        f"Extracted {len(metadata)} gorilla ACREs from "
        f"{metadata['array_id'].nunique()} parent arrays and "
        f"{metadata['synteny_group'].nunique()} chromosome groups"
    )


if __name__ == "__main__":
    main()
