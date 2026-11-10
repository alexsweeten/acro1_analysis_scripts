#!/usr/bin/env python3
"""Build a cross-assembly alignment of de novo ACRO1 monomer consensuses."""

from __future__ import annotations

import csv
import math
import subprocess
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "primate_acro1_orf_test"
RC_TABLE = str.maketrans("ACGTN", "TGCAN")
CODONS = [left + middle + right for left in "TCAG" for middle in "TCAG" for right in "TCAG"]
AMINO_ACIDS = "FFLLSSSSYY**CC*W" "LLLLPPPPHHQQRRRR" "IIIMTTTTNNKKSSRR" "VVVVAAAADDEEGGGG"
GENETIC_CODE = dict(zip(CODONS, AMINO_ACIDS))


def reverse_complement(sequence: str) -> str:
    return sequence.translate(RC_TABLE)[::-1]


def phase_to_motif(sequence: str, motif: str) -> tuple[str, str, int, int]:
    """Orient and circularly phase a monomer to a de novo reference motif."""
    best: tuple[int, str, int, str] | None = None
    for strand, oriented in [("+", sequence), ("-", reverse_complement(sequence))]:
        circular = oriented + oriented[: len(motif) - 1]
        for position in range(len(oriented)):
            window = circular[position : position + len(motif)]
            mismatches = sum(left != right for left, right in zip(window, motif))
            candidate = (mismatches, strand, position, oriented[position:] + oriented[:position])
            if best is None or candidate[:3] < best[:3]:
                best = candidate
    assert best is not None
    return best[3], best[1], best[2], best[0]


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    name = None
    parts: list[str] = []
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(parts)
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line.upper())
    if name is not None:
        records[name] = "".join(parts)
    return records


rows = [
    row for row in csv.DictReader((OUT / "array_orf_results.tsv").open(), delimiter="\t")
    if row["status"] == "ok"
]
reference_row = max(
    (row for row in rows if row["key"] == "CHM13" and len(row["monomer_consensus"]) == 147),
    key=lambda row: int(row["representative_score"]),
)
reference = reference_row["monomer_consensus"]
motif = reference[:21]
metadata: dict[str, tuple[str, str, str]] = {}
phase_rows = []
unaligned = OUT / "monomer_consensuses.fa"
with unaligned.open("w") as handle:
    for index, row in enumerate(rows, start=1):
        name = f"m{index:04d}"
        metadata[name] = (row["key"], row["species"], row["array_id"])
        phased, strand, rotation, mismatches = phase_to_motif(row["monomer_consensus"], motif)
        phase_rows.append([name, row["key"], row["array_id"], strand, rotation, mismatches, len(phased)])
        handle.write(f">{name}\n{phased}\n")

with (OUT / "monomer_phase_map.tsv").open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["sequence_id", "key", "array_id", "orientation", "rotation_bp", "motif_mismatches", "monomer_length"])
    writer.writerows(phase_rows)

aligned = OUT / "monomer_consensuses.aligned.fa"
with aligned.open("w") as handle:
    subprocess.run(
        ["mafft", "--quiet", "--auto", str(unaligned)],
        check=True,
        stdout=handle,
    )

alignment = read_fasta(aligned)
normalized: dict[str, str] = {}
for name, sequence in alignment.items():
    clean = name.removeprefix("_R_")
    normalized[clean] = sequence

widths = {len(sequence) for sequence in normalized.values()}
if len(widths) != 1:
    raise SystemExit("MAFFT output has inconsistent widths")
width = widths.pop()

by_key: dict[str, list[str]] = defaultdict(list)
for name, sequence in normalized.items():
    by_key[metadata[name][0]].append(sequence)

species_consensus: dict[str, str] = {}
consensus_rows = []
for key, sequences in by_key.items():
    consensus_alignment: list[str] = []
    informative = 0
    for column in range(width):
        counts = Counter(sequence[column] for sequence in sequences if sequence[column] in "ACGT")
        if sum(counts.values()) < math.ceil(len(sequences) / 2):
            consensus_alignment.append("-")
            continue
        informative += 1
        consensus_alignment.append(counts.most_common(1)[0][0])
    aligned_sequence = "".join(consensus_alignment)
    ungapped = aligned_sequence.replace("-", "")
    species_consensus[key] = ungapped
    example = next(row for row in rows if row["key"] == key)
    consensus_rows.append([
        key,
        example["species"],
        len(sequences),
        width,
        informative,
        len(ungapped),
        ungapped,
    ])

with (OUT / "assembly_monomer_consensuses.tsv").open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["key", "species", "array_consensuses", "alignment_width", "informative_columns", "consensus_length", "consensus_sequence"])
    writer.writerows(consensus_rows)

with (OUT / "assembly_monomer_consensuses.fa").open("w") as handle:
    for key, sequence in species_consensus.items():
        handle.write(f">{key}\n{sequence}\n")

peptide_rows = []
with (OUT / "assembly_candidate_orf_peptides.fa").open("w") as handle:
    for key, sequence in species_consensus.items():
        peptide = "".join(GENETIC_CODE.get(sequence[index:index + 3], "X") for index in range(0, len(sequence) - 2, 3))
        example = next(row for row in rows if row["key"] == key)
        peptide_rows.append([
            key,
            example["species"],
            len(sequence),
            len(peptide),
            peptide.count("*"),
            peptide.find("M") + 1 if "M" in peptide else "NA",
            sum(residue in "PST" for residue in peptide) / len(peptide),
            peptide,
        ])
        handle.write(f">{key}_candidate_plus0\n{peptide}\n")

with (OUT / "assembly_candidate_orf_peptides.tsv").open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["key", "species", "monomer_bp", "peptide_aa", "internal_stops", "first_methionine_aa", "PST_fraction", "peptide"])
    writer.writerows(peptide_rows)

with (OUT / "assembly_candidate_orf_peptides.aligned.fa").open("w") as handle:
    subprocess.run(
        ["mafft", "--quiet", "--auto", str(OUT / "assembly_candidate_orf_peptides.fa")],
        check=True,
        stdout=handle,
    )

print(aligned)
