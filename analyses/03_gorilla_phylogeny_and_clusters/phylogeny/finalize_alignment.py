#!/usr/bin/env python3
"""Normalize MAFFT orientation labels and create occupancy-filtered alignments."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
FILTER_THRESHOLDS = [0.02, 0.05, 0.10]


def records(path: Path):
    name: str | None = None
    parts: list[str] = []
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts)
                name = line[1:].strip().split()[0]
                parts = []
            else:
                parts.append(line.strip())
    if name is not None:
        yield name, "".join(parts)


def write_fasta(path: Path, names: list[str], matrix: np.ndarray) -> None:
    with path.open("w") as handle:
        for row, name in enumerate(names):
            sequence = matrix[row].tobytes().decode("ascii")
            handle.write(f">{name}\n")
            for offset in range(0, len(sequence), 80):
                handle.write(sequence[offset : offset + 80] + "\n")


def main() -> None:
    source = HERE / "gorilla_acres.aligned.raw.fasta"
    output = HERE / "gorilla_acres.aligned.fasta"
    metadata = HERE / "gorilla_tip_metadata.tsv"
    expected = set()
    with metadata.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            expected.add(row["tip_id"])

    names: list[str] = []
    sequences: list[str] = []
    orientations: list[str] = []
    seen: set[str] = set()
    alignment_length: int | None = None
    for raw_name, raw_sequence in records(source):
        reverse = raw_name.startswith("_R_")
        name = raw_name[3:] if reverse else raw_name
        if name in seen:
            raise RuntimeError(f"Duplicate aligned ID: {name}")
        sequence = raw_sequence.upper()
        if alignment_length is None:
            alignment_length = len(sequence)
        elif len(sequence) != alignment_length:
            raise RuntimeError("Alignment rows have inconsistent lengths")
        names.append(name)
        sequences.append(sequence)
        orientations.append("reverse_complement" if reverse else "forward")
        seen.add(name)
    if seen != expected:
        raise RuntimeError(
            f"Alignment/metadata mismatch: {len(seen)} aligned, {len(expected)} expected"
        )
    if alignment_length is None:
        raise RuntimeError("Empty alignment")

    matrix = np.frombuffer("".join(sequences).encode("ascii"), dtype="S1").reshape(
        len(sequences), alignment_length
    )
    write_fasta(output, names, matrix)
    with (HERE / "gorilla_orientation.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["tip_id", "mafft_orientation"])
        writer.writerows(zip(names, orientations))

    valid = np.isin(matrix, np.array([b"A", b"C", b"G", b"T"], dtype="S1"))
    occupancy = valid.mean(axis=0)
    filter_rows: list[list[object]] = []
    masks: dict[float, np.ndarray] = {}
    for threshold in FILTER_THRESHOLDS:
        keep = occupancy >= threshold
        masks[threshold] = keep
        label = f"{int(round(threshold * 100)):02d}"
        filtered = matrix[:, keep]
        write_fasta(HERE / f"gorilla_acres.aligned.occupancy{label}.fasta", names, filtered)
        valid_per_sequence = valid[:, keep].sum(axis=1)
        gap_or_n = np.isin(filtered, np.array([b"-", b"N"], dtype="S1")).mean()
        filter_rows.append(
            [
                threshold,
                int(keep.sum()),
                keep.mean(),
                gap_or_n,
                int(valid_per_sequence.min()),
                float(np.median(valid_per_sequence)),
                int((valid_per_sequence == 0).sum()),
            ]
        )

    with (HERE / "alignment_column_occupancy.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            ["alignment_column_1based", "valid_base_count", "valid_base_fraction", "retained_occupancy10"]
        )
        for index, fraction in enumerate(occupancy, 1):
            writer.writerow(
                [index, int(valid[:, index - 1].sum()), fraction, int(masks[0.10][index - 1])]
            )

    with (HERE / "alignment_filter_summary.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "occupancy_threshold",
                "retained_columns",
                "retained_fraction",
                "gap_or_N_fraction",
                "minimum_valid_bases_per_sequence",
                "median_valid_bases_per_sequence",
                "all_gap_sequences",
            ]
        )
        writer.writerows(filter_rows)

    gap_or_n_full = np.isin(matrix, np.array([b"-", b"N"], dtype="S1")).mean()
    with (HERE / "alignment_summary.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sequences", "alignment_columns", "reverse_complemented", "gap_or_N_fraction"])
        writer.writerow([len(names), alignment_length, orientations.count("reverse_complement"), gap_or_n_full])
    print(
        f"Validated {len(names)} sequences x {alignment_length} columns; "
        f"reverse-complemented={orientations.count('reverse_complement')}; "
        f"gap/N={gap_or_n_full:.3%}"
    )


if __name__ == "__main__":
    main()
