#!/usr/bin/env python3
"""Validate and combine the 2,794 telomeric ACRO-array FASTAs for N-TRIAGE."""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import re
from pathlib import Path


VALID_BASES = set("ACGTRYSWKMBDHVN")
COORD_RE = re.compile(r"^(.*)_([0-9]+)-([0-9]+)$")


def read_one_fasta(path: Path) -> tuple[str, str]:
    header = ""
    sequence_parts: list[str] = []
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    raise ValueError(f"Multiple FASTA records in {path}")
                header = line[1:].split()[0]
            else:
                if not header:
                    raise ValueError(f"Sequence before header in {path}")
                sequence_parts.append(line.upper())
    sequence = "".join(sequence_parts)
    if not header or not sequence:
        raise ValueError(f"Missing FASTA record or sequence in {path}")
    invalid = sorted(set(sequence) - VALID_BASES)
    if invalid:
        raise ValueError(f"Invalid sequence symbols in {path}: {invalid}")
    return header, sequence


def wrap(sequence: str, width: int = 80) -> str:
    return "\n".join(sequence[i : i + width] for i in range(0, len(sequence), width))


def filename_for_array(array_id: str) -> str:
    return array_id.replace("#", "_") + ".fa"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--fasta-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    combined_fasta = args.out_dir / "telomeric_2794.fasta"
    combined_bed = args.out_dir / "telomeric_2794.bed"
    core_bed = args.out_dir / "telomeric_2794.core.bed"
    audit_tsv = args.out_dir / "input_audit.tsv"
    acro1_fasta = args.out_dir / "ACRO1_147.fasta"

    with args.manifest.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 2794:
        raise ValueError(f"Expected 2,794 manifest rows; observed {len(rows):,}")
    array_ids = [row["array_id"] for row in rows]
    if len(set(array_ids)) != len(array_ids):
        raise ValueError("Manifest array_id values are not unique")

    observed_fastas = sorted(args.fasta_dir.glob("*.fa"))
    if len(observed_fastas) != 2794:
        raise ValueError(f"Expected 2,794 .fa files; observed {len(observed_fastas):,}")
    expected_paths = {args.fasta_dir / filename_for_array(aid) for aid in array_ids}
    observed_paths = set(observed_fastas)
    if expected_paths != observed_paths:
        missing = sorted(str(p) for p in expected_paths - observed_paths)
        extra = sorted(str(p) for p in observed_paths - expected_paths)
        raise ValueError(f"Manifest/FASTA mismatch; missing={missing[:10]}, extra={extra[:10]}")

    total_bp = 0
    audit_columns = [
        "task_index", "array_id", "individual_id", "haplotype", "population",
        "superpopulation", "has_recombination", "recombination_class", "source_fasta",
        "original_header", "sequence_length_bp", "manifest_fasta_length",
        "length_matches_manifest", "coordinate_closed_span_bp", "length_matches_header_span",
        "core_local_start", "core_local_end", "core_length_bp",
        "left_flank_bp", "right_flank_bp", "sequence_sha256", "ambiguous_bases"
    ]

    with (
        combined_fasta.open("w") as fasta_out,
        combined_bed.open("w") as bed_out,
        core_bed.open("w") as core_bed_out,
        audit_tsv.open("w", newline="") as audit_handle,
    ):
        writer = csv.DictWriter(audit_handle, fieldnames=audit_columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for index, row in enumerate(rows, start=1):
            array_id = row["array_id"]
            source_path = args.fasta_dir / filename_for_array(array_id)
            original_header, sequence = read_one_fasta(source_path)
            sequence_length = len(sequence)
            manifest_length = int(row["fasta_length"])
            match = COORD_RE.match(array_id)
            closed_span = -1
            if match:
                closed_span = int(match.group(3)) - int(match.group(2)) + 1
            if sequence_length != manifest_length:
                raise ValueError(
                    f"Length mismatch for {array_id}: FASTA={sequence_length}, manifest={manifest_length}"
                )
            if closed_span != sequence_length:
                raise ValueError(
                    f"Coordinate-span mismatch for {array_id}: span={closed_span}, FASTA={sequence_length}"
                )

            original_start = int(row["original_start"])
            original_end = int(row["original_end"])
            flanked_start = int(row["flanked_start"])
            flanked_end = int(row["flanked_end"])
            # original_start/original_end come from a BED interval (0-based,
            # half-open), whereas flanked_start/flanked_end name the region
            # passed to samtools faidx (1-based, inclusive).  FASTA offset 0
            # therefore corresponds to genomic position flanked_start, and
            # the BED start corresponds to 1-based position original_start+1.
            core_local_start = original_start - flanked_start + 1
            core_local_end = original_end - flanked_start + 1
            core_length = core_local_end - core_local_start
            left_flank = core_local_start
            right_flank = sequence_length - core_local_end
            if not (0 <= core_local_start < core_local_end <= sequence_length):
                raise ValueError(f"Invalid local core interval for {array_id}")
            if core_length != original_end - original_start:
                raise ValueError(f"Core length mismatch for {array_id}")
            if flanked_end - flanked_start + 1 != sequence_length:
                raise ValueError(f"Flanked manifest span mismatch for {array_id}")

            fasta_out.write(f">{array_id}\n{wrap(sequence)}\n")
            bed_out.write(f"{array_id}\t0\t{sequence_length}\t{array_id}\t0\t.\n")
            core_bed_out.write(
                f"{array_id}\t{core_local_start}\t{core_local_end}\t{array_id}\t0\t.\n"
            )
            total_bp += sequence_length
            writer.writerow({
                "task_index": index,
                "array_id": array_id,
                "individual_id": row["individual_id"],
                "haplotype": row["haplotype"],
                "population": row["population"],
                "superpopulation": row["superpopulation"],
                "has_recombination": row["has_recombination"],
                "recombination_class": row["recombination_class"],
                "source_fasta": str(source_path),
                "original_header": original_header,
                "sequence_length_bp": sequence_length,
                "manifest_fasta_length": manifest_length,
                "length_matches_manifest": 1,
                "coordinate_closed_span_bp": closed_span,
                "length_matches_header_span": 1,
                "core_local_start": core_local_start,
                "core_local_end": core_local_end,
                "core_length_bp": core_length,
                "left_flank_bp": left_flank,
                "right_flank_bp": right_flank,
                "sequence_sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest(),
                "ambiguous_bases": sum(base not in "ACGT" for base in sequence),
            })

    # Canonical ACRO1 sequence from the validated N-TRIAGE reference collection.
    refs = Path(os.environ.get("ACRO1_COMPONENT_REFERENCES", "component_references.fasta"))
    current = ""
    pieces: list[str] = []
    with refs.open() as handle:
        for raw in handle:
            line = raw.strip()
            if line.startswith(">"):
                if current == "ACRO1":
                    break
                current = line[1:].split()[0]
            elif current == "ACRO1":
                pieces.append(line.upper())
    acro1 = "".join(pieces)
    if len(acro1) != 147:
        raise ValueError(f"Expected a 147-bp ACRO1 reference; observed {len(acro1)}")
    acro1_fasta.write_text(f">ACRO1\n{wrap(acro1)}\n")

    empty_placeholders = sum(1 for path in args.fasta_dir.iterdir() if path.is_file() and path.suffix != ".fa" and path.stat().st_size == 0)
    print(f"Validated and combined {len(rows):,} arrays totaling {total_bp:,} bp")
    print(f"Ignored {empty_placeholders:,} zero-byte non-.fa placeholders")
    print(f"Combined FASTA: {combined_fasta}")
    print(f"Combined BED: {combined_bed}")
    print(f"Core-only BED: {core_bed}")


if __name__ == "__main__":
    main()
