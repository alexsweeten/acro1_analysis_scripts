#!/usr/bin/env python3
"""Compute ACRE-scale higher-order-repeat (HOR) scores.

The implementation generalizes the workflow used for the gorilla and HPRC
analyses. Complete ACREs are ordered within each parent array and aligned with
MAFFT. At every non-zero inter-unit offset, a qualifying block is a run of at
least ``--min-block-length`` consecutive aligned ACRE pairs meeting the chosen
identity threshold. Each unit's distinct partners across all qualifying blocks
are counted once, bidirectionally, and normalized by the number of complete
ACREs in its parent array.

Identity is implemented as in the original analysis: an aligned pair passes a
threshold when its number of unequal alignment characters is no greater than
``floor((1 - identity) * mean_unaligned_array_unit_length)``. Consequently,
the stated percentages are nominal identity thresholds rather than estimates
from gap-filtered pairwise alignments.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import random
import statistics
import subprocess
import tempfile
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Iterable


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True, type=Path, help="Unit metadata TSV")
    parser.add_argument("--fasta", required=True, type=Path, help="Complete-unit FASTA")
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--id-column", default="acre_id")
    parser.add_argument("--array-column", default="array_id")
    parser.add_argument("--order-column", default="unit_in_array")
    parser.add_argument("--cluster-column", default="cluster")
    parser.add_argument(
        "--identity-thresholds",
        default="90,92,96",
        help="Comma-separated nominal identity percentages (default: 90,92,96)",
    )
    parser.add_argument("--min-block-length", type=int, default=4)
    parser.add_argument("--workers", type=int, default=max(1, os.cpu_count() or 1))
    parser.add_argument("--mafft", default="mafft", help="MAFFT executable")
    parser.add_argument("--bootstrap-reps", type=int, default=20_000)
    parser.add_argument("--seed", type=int, default=20261105)
    return parser.parse_args()


def parse_identity_thresholds(value: str) -> tuple[float, ...]:
    try:
        thresholds = tuple(float(item.strip()) for item in value.split(",") if item.strip())
    except ValueError as exc:
        raise argparse.ArgumentTypeError("identity thresholds must be numeric") from exc
    if not thresholds or any(item <= 0 or item > 100 for item in thresholds):
        raise argparse.ArgumentTypeError("identity thresholds must fall in (0, 100]")
    if len(set(thresholds)) != len(thresholds):
        raise argparse.ArgumentTypeError("identity thresholds must be unique")
    return tuple(sorted(thresholds))


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"No rows to write: {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_fasta(path: Path, wanted: set[str]) -> dict[str, str]:
    """Read wanted FASTA records using the legacy first-token/pipe-prefix IDs."""
    sequences: dict[str, str] = {}
    name: str | None = None
    parts: list[str] = []
    keep = False
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None and keep:
                    sequences[name] = "".join(parts).upper()
                name = line[1:].split()[0].split("|", 1)[0]
                if name in sequences:
                    raise ValueError(f"Duplicate FASTA identifier: {name}")
                keep = name in wanted
                parts = []
            elif keep:
                parts.append(line)
        if name is not None and keep:
            sequences[name] = "".join(parts).upper()
    return sequences


def parse_alignment(text: str) -> dict[str, str]:
    sequences: dict[str, str] = {}
    name: str | None = None
    parts: list[str] = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                sequences[name] = "".join(parts).upper()
            name = line[1:].split()[0]
            parts = []
        else:
            parts.append(line)
    if name is not None:
        sequences[name] = "".join(parts).upper()
    lengths = {len(sequence) for sequence in sequences.values()}
    if len(lengths) != 1:
        raise ValueError(f"MAFFT returned unequal alignment lengths: {sorted(lengths)}")
    return sequences


def mafft_align(records: list[tuple[str, str]], executable: str) -> dict[str, str]:
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as handle:
        temporary_path = Path(handle.name)
        for name, sequence in records:
            handle.write(f">{name}\n")
            for offset in range(0, len(sequence), 80):
                handle.write(sequence[offset : offset + 80] + "\n")
    try:
        completed = subprocess.run(
            [
                executable,
                "--quiet",
                "--inputorder",
                "--kimura",
                "1",
                "--retree",
                "1",
                str(temporary_path),
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        return parse_alignment(completed.stdout)
    finally:
        temporary_path.unlink(missing_ok=True)


def score_diagonal_runs(
    variant_offsets: list[list[int]],
    mismatch_threshold: int,
    clusters: list[str],
    min_block_length: int = 4,
) -> tuple[list[int], list[int], list[int]]:
    """Return all-partner counts, same-cluster counts, and block lengths."""
    if min_block_length < 1:
        raise ValueError("min_block_length must be positive")
    partners = [set() for _ in clusters]
    same_cluster_partners = [set() for _ in clusters]
    block_lengths: list[int] = []

    for unit_offset, mismatch_counts in enumerate(variant_offsets, start=1):
        flags = [value <= mismatch_threshold for value in mismatch_counts]
        run_start: int | None = None
        for position, passes in enumerate(flags + [False]):
            if passes and run_start is None:
                run_start = position
            elif not passes and run_start is not None:
                run_length = position - run_start
                if run_length >= min_block_length:
                    block_lengths.append(run_length)
                    for left_index in range(run_start, position):
                        right_index = left_index + unit_offset
                        partners[left_index].add(right_index)
                        partners[right_index].add(left_index)
                        if clusters[left_index] == clusters[right_index]:
                            same_cluster_partners[left_index].add(right_index)
                            same_cluster_partners[right_index].add(left_index)
                run_start = None

    return (
        [len(values) for values in partners],
        [len(values) for values in same_cluster_partners],
        block_lengths,
    )


def numeric_order(value: str) -> tuple[int, float | str]:
    try:
        return 0, float(value)
    except ValueError:
        return 1, value


def analyze_array(
    array_id: str,
    rows: list[dict[str, str]],
    sequences: dict[str, str],
    id_column: str,
    cluster_column: str,
    identity_thresholds: tuple[float, ...],
    min_block_length: int,
    mafft_executable: str,
) -> tuple[str, float, dict[float, dict]]:
    records = [(row[id_column], sequences[row[id_column]]) for row in rows]
    mean_length = statistics.fmean(len(sequence) for _, sequence in records)
    clusters = [row[cluster_column] for row in rows]

    if len(records) >= 2:
        alignment = mafft_align(records, mafft_executable)
        aligned = [alignment[name] for name, _ in records]
        variant_offsets = [
            [
                sum(left != right for left, right in zip(aligned[index], aligned[index + offset]))
                for index in range(len(records) - offset)
            ]
            for offset in range(1, len(records))
        ]
    else:
        variant_offsets = []

    results: dict[float, dict] = {}
    for identity in identity_thresholds:
        mismatch_fraction = 1.0 - identity / 100.0
        mismatch_threshold = math.floor(mismatch_fraction * mean_length + 1e-12)
        partner_counts, same_cluster_counts, block_lengths = score_diagonal_runs(
            variant_offsets,
            mismatch_threshold,
            clusters,
            min_block_length=min_block_length,
        )
        results[identity] = {
            "mismatch_threshold": mismatch_threshold,
            "partner_counts": partner_counts,
            "same_cluster_counts": same_cluster_counts,
            "block_count": len(block_lengths),
            "maximum_block_length": max(block_lengths) if block_lengths else 0,
        }
    return array_id, mean_length, results


def quantile(values: Iterable[float], probability: float) -> float:
    ordered = sorted(values)
    if not ordered:
        return math.nan
    position = (len(ordered) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    return ordered[lower] * (upper - position) + ordered[upper] * (position - lower)


def bootstrap_ci(
    values_by_array: dict[str, float], reps: int, seed: int
) -> tuple[float, float]:
    array_ids = sorted(values_by_array)
    if not array_ids or reps < 1:
        return math.nan, math.nan
    rng = random.Random(seed)
    estimates = [
        statistics.fmean(
            values_by_array[array_ids[rng.randrange(len(array_ids))]] for _ in array_ids
        )
        for _ in range(reps)
    ]
    return quantile(estimates, 0.025), quantile(estimates, 0.975)


def summarize_clusters(
    unit_rows: list[dict], bootstrap_reps: int, seed: int
) -> list[dict]:
    output: list[dict] = []
    identities = sorted({float(row["identity_threshold_pct"]) for row in unit_rows})
    clusters = sorted({str(row["cluster"]) for row in unit_rows})
    for identity_index, identity in enumerate(identities):
        for cluster_index, cluster in enumerate(clusters):
            selected = [
                row
                for row in unit_rows
                if float(row["identity_threshold_pct"]) == identity
                and row["cluster"] == cluster
            ]
            scores = [float(row["hor_score_pct"]) for row in selected]
            same_scores = [float(row["same_cluster_hor_score_pct"]) for row in selected]
            by_array: dict[str, list[float]] = defaultdict(list)
            for row in selected:
                by_array[row["array_id"]].append(float(row["hor_score_pct"]))
            array_means = {
                array_id: statistics.fmean(values) for array_id, values in by_array.items()
            }
            ci_low, ci_high = bootstrap_ci(
                array_means,
                bootstrap_reps,
                seed + identity_index * 10_000 + cluster_index,
            )
            total_partners = sum(int(row["distinct_partner_count"]) for row in selected)
            total_same = sum(int(row["same_cluster_partner_count"]) for row in selected)
            output.append(
                {
                    "cluster": cluster,
                    "identity_threshold_pct": identity,
                    "units": len(selected),
                    "arrays": len(by_array),
                    "pooled_mean_score_pct": statistics.fmean(scores),
                    "median_score_pct": statistics.median(scores),
                    "q1_score_pct": quantile(scores, 0.25),
                    "q3_score_pct": quantile(scores, 0.75),
                    "pct_units_nonzero": 100.0 * sum(score > 0 for score in scores) / len(scores),
                    "array_balanced_mean_score_pct": statistics.fmean(array_means.values()),
                    "array_balanced_bootstrap95_low": ci_low,
                    "array_balanced_bootstrap95_high": ci_high,
                    "mean_same_cluster_score_pct": statistics.fmean(same_scores),
                    "cross_cluster_partner_pct": (
                        100.0 * (total_partners - total_same) / total_partners
                        if total_partners
                        else 0.0
                    ),
                }
            )
    return output


def main() -> None:
    args = parse_args()
    identity_thresholds = parse_identity_thresholds(args.identity_thresholds)
    if args.min_block_length < 1:
        raise ValueError("--min-block-length must be positive")
    if args.workers < 1:
        raise ValueError("--workers must be positive")

    metadata = read_tsv(args.metadata)
    if not metadata:
        raise ValueError("Metadata table is empty")
    required_columns = {
        args.id_column,
        args.array_column,
        args.order_column,
        args.cluster_column,
    }
    missing_columns = required_columns - set(metadata[0])
    if missing_columns:
        raise ValueError(f"Missing metadata columns: {sorted(missing_columns)}")

    seen_ids: set[str] = set()
    by_array: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in metadata:
        unit_id = row[args.id_column]
        if unit_id in seen_ids:
            raise ValueError(f"Duplicate unit identifier in metadata: {unit_id}")
        seen_ids.add(unit_id)
        by_array[row[args.array_column]].append(row)
    for rows in by_array.values():
        rows.sort(key=lambda row: numeric_order(row[args.order_column]))

    sequences = read_fasta(args.fasta, seen_ids)
    missing_sequences = sorted(seen_ids - set(sequences))
    if missing_sequences:
        raise ValueError(
            f"FASTA is missing {len(missing_sequences)} metadata IDs; "
            f"first examples: {missing_sequences[:5]}"
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    start_time = time.time()
    print(
        f"INPUT\tarrays={len(by_array)}\tunits={len(seen_ids)}\tworkers={args.workers}",
        flush=True,
    )

    unit_output: list[dict] = []
    array_output: list[dict] = []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(
                analyze_array,
                array_id,
                rows,
                sequences,
                args.id_column,
                args.cluster_column,
                identity_thresholds,
                args.min_block_length,
                args.mafft,
            ): array_id
            for array_id, rows in by_array.items()
        }
        for completed_index, future in enumerate(as_completed(futures), start=1):
            array_id, mean_length, results = future.result()
            rows = by_array[array_id]
            array_size = len(rows)
            for identity in identity_thresholds:
                result = results[identity]
                array_output.append(
                    {
                        "array_id": array_id,
                        "identity_threshold_pct": identity,
                        "units": array_size,
                        "mean_acre_length_bp": mean_length,
                        "mismatch_threshold_columns": result["mismatch_threshold"],
                        "qualifying_blocks": result["block_count"],
                        "maximum_block_length": result["maximum_block_length"],
                    }
                )
                for index, row in enumerate(rows):
                    partner_count = result["partner_counts"][index]
                    same_count = result["same_cluster_counts"][index]
                    unit_output.append(
                        {
                            "unit_id": row[args.id_column],
                            "array_id": array_id,
                            "unit_order": row[args.order_column],
                            "cluster": row[args.cluster_column],
                            "array_units": array_size,
                            "identity_threshold_pct": identity,
                            "distinct_partner_count": partner_count,
                            "same_cluster_partner_count": same_count,
                            "hor_score_pct": 100.0 * partner_count / array_size,
                            "same_cluster_hor_score_pct": 100.0 * same_count / array_size,
                        }
                    )
            if completed_index % 100 == 0 or completed_index == len(futures):
                print(f"PROGRESS\t{completed_index}/{len(futures)}", flush=True)

    unit_output.sort(
        key=lambda row: (
            row["array_id"],
            numeric_order(str(row["unit_order"])),
            float(row["identity_threshold_pct"]),
        )
    )
    array_output.sort(key=lambda row: (row["array_id"], row["identity_threshold_pct"]))
    cluster_output = summarize_clusters(unit_output, args.bootstrap_reps, args.seed)

    write_tsv(args.output_dir / "hor_score_per_unit.tsv", unit_output)
    write_tsv(args.output_dir / "hor_score_per_array.tsv", array_output)
    write_tsv(args.output_dir / "hor_score_cluster_summary.tsv", cluster_output)
    manifest = {
        "metadata": str(args.metadata.resolve()),
        "fasta": str(args.fasta.resolve()),
        "id_column": args.id_column,
        "array_column": args.array_column,
        "order_column": args.order_column,
        "cluster_column": args.cluster_column,
        "identity_thresholds_pct": list(identity_thresholds),
        "min_block_length": args.min_block_length,
        "workers": args.workers,
        "mafft": args.mafft,
        "bootstrap_reps": args.bootstrap_reps,
        "seed": args.seed,
        "arrays": len(by_array),
        "units": len(seen_ids),
        "wall_time_seconds": time.time() - start_time,
    }
    with (args.output_dir / "run_manifest.json").open("w") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(
        f"DONE\tarrays={len(by_array)}\tunits={len(seen_ids)}\t"
        f"seconds={manifest['wall_time_seconds']:.1f}",
        flush=True,
    )


if __name__ == "__main__":
    main()
