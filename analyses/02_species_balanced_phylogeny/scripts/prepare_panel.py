#!/usr/bin/env python3
"""Prepare a 50-unit-per-species cross-primate ACRO phylogeny panel."""

from __future__ import annotations

import csv
import os
from collections import defaultdict
from pathlib import Path

import edlib


PROJECT = Path(__file__).resolve().parents[1]
ACREAGE = Path(os.environ.get("ACREAGE_ROOT", PROJECT / "inputs"))
ORIGIN = Path(os.environ.get("HSA19_ORTHOLOGY_ROOT", PROJECT / "hsa19_inputs"))

TARGET_FASTA = ACREAGE / "acro_washu_cn/recombination_haplotype_screen/PAN010_chr13_distal_7kb_reference.fasta"
ARRAY_FASTA = ACREAGE / "primate_full/six_primate_arrays.fasta"
MOTIF_TABLE = ACREAGE / "primate_recombination_haplotype/motif_to_motif_units.tsv"
WINDOW_MANIFEST = ORIGIN / "results/ortholog_window_manifest.tsv"
WINDOW_DIR = ORIGIN / "intermediate/ortholog_windows"
PAF_DIR = ORIGIN / "intermediate/per_window_component_mappings"

GREAT_APE_PER_SPECIES = 50
SIAMANG_PER_HAPLOTYPE = 25
SPECIES_ORDER = [
    "human",
    "gorilla",
    "bonobo",
    "chimpanzee",
    "sumatran_orangutan",
    "bornean_orangutan",
]


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    name: str | None = None
    chunks: list[str] = []
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(chunks).upper()
                name = line[1:].strip().split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        records[name] = "".join(chunks).upper()
    return records


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def wrap(sequence: str, width: int = 80) -> str:
    return "\n".join(sequence[i : i + width] for i in range(0, len(sequence), width))


def global_identity(reference: str, query: str) -> float:
    result = edlib.align(reference, query, mode="NW", task="distance")
    return 1 - int(result["editDistance"]) / max(len(reference), len(query))


def parse_paf(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            if len(fields) < 12:
                continue
            tags = {
                field[:2]: field[5:]
                for field in fields[12:]
                if len(field) > 5 and field[2:5] == ":i:"
            }
            rows.append(
                {
                    "query_length": int(fields[1]),
                    "query_start": int(fields[2]),
                    "query_end": int(fields[3]),
                    "strand": fields[4],
                    "target_start": int(fields[7]),
                    "target_end": int(fields[8]),
                    "matching_bases": int(fields[9]),
                    "alignment_score": int(tags.get("AS", 0)),
                }
            )
    return rows


def nonoverlapping_full_unit_hits(path: Path) -> list[dict[str, object]]:
    candidates = []
    for row in parse_paf(path):
        qcov = (int(row["query_end"]) - int(row["query_start"])) / int(row["query_length"])
        span = int(row["target_end"]) - int(row["target_start"])
        if qcov >= 0.97 and 5_500 <= span <= 7_500:
            candidates.append(row)
    candidates.sort(
        key=lambda row: (
            -int(row["alignment_score"]),
            -int(row["matching_bases"]),
            int(row["target_start"]),
        )
    )
    kept: list[dict[str, object]] = []
    for row in candidates:
        start, end = int(row["target_start"]), int(row["target_end"])
        if any(
            max(0, min(end, int(old["target_end"])) - max(start, int(old["target_start"]))) > 500
            for old in kept
        ):
            continue
        kept.append(row)
    return sorted(kept, key=lambda row: int(row["target_start"]))


def main() -> None:
    target = next(iter(read_fasta(TARGET_FASTA).values()))
    arrays = read_fasta(ARRAY_FASTA)
    motif_rows = read_tsv(MOTIF_TABLE)
    window_manifest = {row["assembly_key"]: row for row in read_tsv(WINDOW_MANIFEST)}

    rows: list[dict[str, object]] = []
    sequences: dict[str, str] = {}

    # Great apes: the 50 closest complete family-38 motif-bounded units in
    # each species. PAN010 is used only as the fixed ranking reference and is
    # not included as a tip.
    for species in SPECIES_ORDER:
        candidates = [row for row in motif_rows if row["species_name"] == species]
        candidates.sort(
            key=lambda row: float(row["global_identity_to_pan010_h2_chr13_7kb"]),
            reverse=True,
        )
        if len(candidates) < GREAT_APE_PER_SPECIES:
            raise RuntimeError(f"Only {len(candidates)} units available for {species}")
        for rank, source in enumerate(candidates[:GREAT_APE_PER_SPECIES], start=1):
            rel_start = int(source["start"]) - int(source["array_start"])
            rel_end = int(source["end"]) - int(source["array_start"])
            sequence = arrays[source["synthetic_chrom"]][rel_start:rel_end]
            if source["orientation"] == "-":
                sequence = reverse_complement(sequence)
            species_key = species.replace("_", "-")
            label = "__".join(
                [
                    species_key,
                    source["chromosome"].replace("_", "-"),
                    source["start"],
                    f"rank{rank:02d}",
                ]
            )
            sequences[label] = sequence
            rows.append(
                {
                    "tip_label": label,
                    "species_key": species_key,
                    "species_display": species.replace("_", " ").title(),
                    "assembly": source["assembly"],
                    "chromosome": source["chromosome"],
                    "start": source["start"],
                    "end": source["end"],
                    "length_bp": len(sequence),
                    "selection_source": "top50_family38_motif_to_motif",
                    "completeness": "complete",
                    "identity_to_pan010_ranking_reference": source["global_identity_to_pan010_h2_chr13_7kb"],
                    "within_species_rank": rank,
                }
            )

    # Siamang: select 25 highest-identity near-complete non-overlapping units
    # from each phased HSA19-orthologous array.
    for assembly_key, species_key, display in [
        ("Siamang_hap1", "siamang-hap1", "Siamang haplotype 1"),
        ("Siamang_hap2", "siamang-hap2", "Siamang haplotype 2"),
    ]:
        window = next(iter(read_fasta(WINDOW_DIR / f"{assembly_key}.fa").values()))
        hits = nonoverlapping_full_unit_hits(PAF_DIR / f"{assembly_key}.unit.paf")
        scored: list[tuple[float, dict[str, object], str]] = []
        for hit in hits:
            start, end = int(hit["target_start"]), int(hit["target_end"])
            sequence = window[start:end]
            if hit["strand"] == "-":
                sequence = reverse_complement(sequence)
            scored.append((global_identity(target, sequence), hit, sequence))
        scored.sort(key=lambda item: item[0], reverse=True)
        if len(scored) < SIAMANG_PER_HAPLOTYPE:
            raise RuntimeError(f"Only {len(scored)} full mappings found for {assembly_key}")
        manifest = window_manifest[assembly_key]
        for rank, (identity, hit, sequence) in enumerate(scored[:SIAMANG_PER_HAPLOTYPE], start=1):
            genomic_start = int(manifest["window_start"]) + int(hit["target_start"])
            genomic_end = int(manifest["window_start"]) + int(hit["target_end"])
            label = "__".join([species_key, manifest["chrom"], str(genomic_start), f"rank{rank:02d}"])
            sequences[label] = sequence
            rows.append(
                {
                    "tip_label": label,
                    "species_key": species_key,
                    "species_display": display,
                    "assembly": assembly_key,
                    "chromosome": manifest["chrom"],
                    "start": genomic_start,
                    "end": genomic_end,
                    "length_bp": len(sequence),
                    "selection_source": "near_complete_reference_unit_mapping",
                    "completeness": "complete",
                    "identity_to_pan010_ranking_reference": identity,
                    "within_species_rank": rank,
                }
            )

    # Marmoset lacks a complete unit. Retain its best composite fragment as a
    # visibly partial rooting outgroup.
    window = next(iter(read_fasta(WINDOW_DIR / "Marmoset.fa").values()))
    best = max(
        parse_paf(PAF_DIR / "Marmoset.unit.paf"),
        key=lambda row: (int(row["matching_bases"]), int(row["alignment_score"])),
    )
    sequence = window[int(best["target_start"]):int(best["target_end"])]
    if best["strand"] == "-":
        sequence = reverse_complement(sequence)
    manifest = window_manifest["Marmoset"]
    genomic_start = int(manifest["window_start"]) + int(best["target_start"])
    genomic_end = int(manifest["window_start"]) + int(best["target_end"])
    label = "__".join(["marmoset-partial", manifest["chrom"], str(genomic_start), "best-fragment"])
    sequences[label] = sequence
    rows.append(
        {
            "tip_label": label,
            "species_key": "marmoset-partial",
            "species_display": "Marmoset (partial tract)",
            "assembly": "GCA_049354715.1",
            "chromosome": manifest["chrom"],
            "start": genomic_start,
            "end": genomic_end,
            "length_bp": len(sequence),
            "selection_source": "best_partial_reference_unit_mapping",
            "completeness": "partial",
            "identity_to_pan010_ranking_reference": global_identity(target, sequence),
            "within_species_rank": 1,
        }
    )

    write_tsv(PROJECT / "inputs/panel_manifest.tsv", rows)
    with (PROJECT / "inputs/panel_sequences.fasta").open("w") as handle:
        for row in rows:
            label = str(row["tip_label"])
            handle.write(f">{label}\n{wrap(sequences[label])}\n")

    by_species: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        by_species[str(row["species_key"])].append(row)
    summary = []
    for species, group in by_species.items():
        summary.append(
            {
                "species_key": species,
                "n_sequences": len(group),
                "complete_sequences": sum(row["completeness"] == "complete" for row in group),
                "partial_sequences": sum(row["completeness"] == "partial" for row in group),
                "minimum_length_bp": min(int(row["length_bp"]) for row in group),
                "maximum_length_bp": max(int(row["length_bp"]) for row in group),
            }
        )
    write_tsv(PROJECT / "results/panel_summary.tsv", summary)
    print(f"Wrote {len(rows)} sequences")
    for row in summary:
        print(row)


if __name__ == "__main__":
    main()
