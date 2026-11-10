#!/usr/bin/env python3
"""Compute publication-style HOR scores for the six gorilla NTRiage clusters.

This mirrors the relaxed human analysis:
  * complete ACREs are treated as monomers;
  * ACREs are ordered within their original parent array;
  * sequences are aligned per array with MAFFT;
  * qualifying matching blocks contain at least four consecutive ACRE pairs;
  * the primary pairwise identity threshold is 90%, with 92% and 96%
    sensitivity analyses;
  * partners are counted bidirectionally, self-matches are excluded, and each
    ACRE score is normalized by the number of ACREs in its parent array.
"""

from __future__ import annotations

import csv
import math
import os
import random
import statistics
import subprocess
import tempfile
import threading
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


METADATA = Path(os.environ.get("GORILLA_ACRE_METADATA", "gorilla_tip_metadata.tsv"))
FASTA = Path(os.environ.get("GORILLA_ACRE_FASTA", "gorilla_acres.fasta"))
OUTDIR = Path(__file__).resolve().parent
WORKERS = int(os.environ.get("WORKERS", "8"))

# Divergence fractions corresponding to 96%, 92%, and 90% identity.
CUTOFFS = (0.04, 0.08, 0.10)
PRIMARY = 0.10
CLUSTERS = tuple(f"G{i}" for i in range(1, 7))


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_fasta(path: Path, wanted: set[str]) -> dict[str, str]:
    sequences: dict[str, str] = {}
    name = None
    parts: list[str] = []
    keep = False
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if line.startswith(">"):
                if name is not None and keep:
                    sequences[name] = "".join(parts).upper()
                name = line[1:].split("|", 1)[0].split()[0]
                keep = name in wanted
                parts = []
            elif keep and line:
                parts.append(line)
        if name is not None and keep:
            sequences[name] = "".join(parts).upper()
    return sequences


def parse_alignment(text: str) -> dict[str, str]:
    sequences: dict[str, str] = {}
    name = None
    parts: list[str] = []
    for raw in text.splitlines():
        line = raw.strip()
        if line.startswith(">"):
            if name is not None:
                sequences[name] = "".join(parts).upper()
            name = line[1:].split()[0]
            parts = []
        elif line:
            parts.append(line)
    if name is not None:
        sequences[name] = "".join(parts).upper()
    lengths = {len(seq) for seq in sequences.values()}
    if len(lengths) != 1:
        raise ValueError(f"MAFFT returned unequal alignment lengths: {lengths}")
    return sequences


def mafft_align(records: list[tuple[str, str]]) -> dict[str, str]:
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as handle:
        path = handle.name
        for name, sequence in records:
            handle.write(f">{name}\n")
            for i in range(0, len(sequence), 80):
                handle.write(sequence[i : i + 80] + "\n")
    try:
        completed = subprocess.run(
            [
                "mafft",
                "--quiet",
                "--inputorder",
                "--kimura",
                "1",
                "--retree",
                "1",
                path,
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        return parse_alignment(completed.stdout)
    finally:
        os.unlink(path)


def score_diagonal_runs(
    variant_offsets: list[list[int]], threshold: int, clusters: list[str]
) -> tuple[list[int], list[int], list[int]]:
    n = len(clusters)
    partners = [set() for _ in range(n)]
    same_cluster_partners = [set() for _ in range(n)]
    block_lengths: list[int] = []

    for offset, variants in enumerate(variant_offsets, start=1):
        flags = [value <= threshold for value in variants]
        run_start = None
        for position, flag in enumerate(flags + [False]):
            if flag and run_start is None:
                run_start = position
            elif not flag and run_start is not None:
                run_length = position - run_start
                if run_length >= 4:
                    block_lengths.append(run_length)
                    for i in range(run_start, position):
                        j = i + offset
                        partners[i].add(j)
                        partners[j].add(i)
                        if clusters[i] == clusters[j]:
                            same_cluster_partners[i].add(j)
                            same_cluster_partners[j].add(i)
                run_start = None

    return (
        [len(values) for values in partners],
        [len(values) for values in same_cluster_partners],
        block_lengths,
    )


def analyze_array(item: tuple[str, list[dict[str, str]], dict[str, str]]):
    array_id, rows, sequences = item
    n = len(rows)
    records = [(row["acre_id"], sequences[row["acre_id"]]) for row in rows]
    mean_length = statistics.fmean(len(sequence) for _, sequence in records)
    clusters = [row["cluster"] for row in rows]

    if n >= 2:
        alignment = mafft_align(records)
        aligned = [alignment[name] for name, _ in records]
        variant_offsets = []
        for offset in range(1, n):
            variant_offsets.append(
                [
                    sum(a != b for a, b in zip(aligned[i], aligned[i + offset]))
                    for i in range(n - offset)
                ]
            )
    else:
        variant_offsets = []

    results = {}
    for cutoff in CUTOFFS:
        threshold = math.floor(cutoff * mean_length + 1e-12)
        partners, same_cluster, block_lengths = score_diagonal_runs(
            variant_offsets, threshold, clusters
        )
        results[cutoff] = {
            "threshold": threshold,
            "partners": partners,
            "same_cluster": same_cluster,
            "block_count": len(block_lengths),
            "max_block": max(block_lengths) if block_lengths else 0,
        }
    return array_id, n, mean_length, results


def quantile(values: list[float], q: float) -> float:
    if not values:
        return math.nan
    ordered = sorted(values)
    position = (len(ordered) - 1) * q
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    return ordered[lower] * (upper - position) + ordered[upper] * (position - lower)


def array_balanced_bootstrap_ci(
    array_means: dict[str, float], reps: int = 20_000, seed: int = 20261105
) -> tuple[float, float]:
    array_ids = sorted(array_means)
    if not array_ids:
        return math.nan, math.nan
    rng = random.Random(seed)
    estimates = []
    for _ in range(reps):
        estimates.append(
            statistics.fmean(
                array_means[array_ids[rng.randrange(len(array_ids))]]
                for _ in array_ids
            )
        )
    return quantile(estimates, 0.025), quantile(estimates, 0.975)


def cluster_summary(
    units: list[dict], cutoff: float, eligible_only: bool = False
) -> list[dict]:
    rows_out = []
    score_field = f"score_{cutoff}"
    partner_field = f"partners_{cutoff}"
    same_field = f"same_partners_{cutoff}"
    same_score_field = f"same_score_{cutoff}"

    for cluster in CLUSTERS:
        selected = [
            row
            for row in units
            if row["cluster"] == cluster
            and (not eligible_only or row["array_n"] >= 4)
        ]
        scores = [row[score_field] for row in selected]
        same_scores = [row[same_score_field] for row in selected]
        by_array: dict[str, list[float]] = defaultdict(list)
        for row in selected:
            by_array[row["array_id"]].append(row[score_field])
        array_means = {
            array_id: statistics.fmean(values) for array_id, values in by_array.items()
        }
        lo, hi = array_balanced_bootstrap_ci(
            array_means,
            seed=20261105 + int(cluster[1:]) + int(cutoff * 1000) + eligible_only,
        )
        total_partners = sum(row[partner_field] for row in selected)
        total_same = sum(row[same_field] for row in selected)
        cross_partner_pct = (
            100 * (total_partners - total_same) / total_partners
            if total_partners
            else 0.0
        )
        rows_out.append(
            {
                "cluster": cluster,
                "identity_threshold_pct": int(round(100 * (1 - cutoff))),
                "scope": "arrays_n_ge_4" if eligible_only else "all_arrays",
                "units": len(selected),
                "arrays": len(by_array),
                "arrays_n_ge_4": len(
                    {row["array_id"] for row in selected if row["array_n"] >= 4}
                ),
                "pooled_mean_score_pct": statistics.fmean(scores),
                "median_score_pct": statistics.median(scores),
                "q1_score_pct": quantile(scores, 0.25),
                "q3_score_pct": quantile(scores, 0.75),
                "pct_units_nonzero": 100 * sum(value > 0 for value in scores) / len(scores),
                "array_balanced_mean_score_pct": statistics.fmean(array_means.values()),
                "array_balanced_bootstrap95_lo": lo,
                "array_balanced_bootstrap95_hi": hi,
                "mean_same_cluster_score_pct": statistics.fmean(same_scores),
                "cross_cluster_partner_pct": cross_partner_pct,
            }
        )
    return rows_out


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"No rows to write to {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    start_time = time.time()
    metadata = read_tsv(METADATA)
    if len(metadata) != 1232:
        raise ValueError(f"Expected 1,232 gorilla ACREs, found {len(metadata)}")

    by_array: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in metadata:
        row["cluster"] = f'G{row["de_novo_cluster"]}'
        by_array[row["array_id"]].append(row)
    for rows in by_array.values():
        rows.sort(key=lambda row: int(row["unit_in_array"]))

    wanted = {row["acre_id"] for rows in by_array.values() for row in rows}
    sequences = read_fasta(FASTA, wanted)
    if set(sequences) != wanted:
        missing = sorted(wanted - set(sequences))
        raise ValueError(f"Missing {len(missing)} FASTA records: {missing[:5]}")

    array_ids = sorted(by_array, key=lambda value: int(value.split("_")[-1]))
    print(
        f"INPUT\tarrays={len(array_ids)}\tunits={len(wanted)}\tworkers={WORKERS}",
        flush=True,
    )

    unit_results: list[dict] = []
    array_results: list[dict] = []
    tasks = [(array_id, by_array[array_id], sequences) for array_id in array_ids]
    completed_arrays = 0

    with ThreadPoolExecutor(max_workers=WORKERS) as pool:
        futures = {pool.submit(analyze_array, task): task[0] for task in tasks}
        for future in as_completed(futures):
            array_id, n, mean_length, results = future.result()
            rows = by_array[array_id]
            for index, row in enumerate(rows):
                unit = {
                    "acre_id": row["acre_id"],
                    "array_id": array_id,
                    "unit_in_array": int(row["unit_in_array"]),
                    "cluster": row["cluster"],
                    "original_chrom": row["original_chrom"],
                    "synteny_group": row["synteny_group"],
                    "haplotype": row["haplotype"],
                    "nor_status": row["nor_status"],
                    "array_n": n,
                }
                for cutoff in CUTOFFS:
                    partners = results[cutoff]["partners"][index]
                    same = results[cutoff]["same_cluster"][index]
                    unit[f"partners_{cutoff}"] = partners
                    unit[f"same_partners_{cutoff}"] = same
                    unit[f"score_{cutoff}"] = 100 * partners / n
                    unit[f"same_score_{cutoff}"] = 100 * same / n
                unit_results.append(unit)

            array_row = {
                "array_id": array_id,
                "units": n,
                "mean_acre_length_bp": mean_length,
                "original_chrom": rows[0]["original_chrom"],
                "synteny_group": rows[0]["synteny_group"],
                "haplotype": rows[0]["haplotype"],
                "nor_status": rows[0]["nor_status"],
                "clusters_present": ",".join(sorted({row["cluster"] for row in rows})),
            }
            for cutoff in CUTOFFS:
                array_row[f"variant_threshold_{cutoff}"] = results[cutoff]["threshold"]
                array_row[f"qualifying_blocks_{cutoff}"] = results[cutoff]["block_count"]
                array_row[f"max_block_length_{cutoff}"] = results[cutoff]["max_block"]
            array_results.append(array_row)

            completed_arrays += 1
            if completed_arrays % 20 == 0 or completed_arrays == len(tasks):
                elapsed = time.time() - start_time
                print(
                    f"PROGRESS\t{completed_arrays}/{len(tasks)}\telapsed_s={elapsed:.1f}",
                    flush=True,
                )

    unit_results.sort(
        key=lambda row: (int(row["array_id"].split("_")[-1]), row["unit_in_array"])
    )
    array_results.sort(key=lambda row: int(row["array_id"].split("_")[-1]))

    summaries: list[dict] = []
    for cutoff in CUTOFFS:
        summaries.extend(cluster_summary(unit_results, cutoff, eligible_only=False))
        summaries.extend(cluster_summary(unit_results, cutoff, eligible_only=True))

    write_tsv(OUTDIR / "gorilla_hor_score_per_unit.tsv", unit_results)
    write_tsv(OUTDIR / "gorilla_hor_score_per_array.tsv", array_results)
    write_tsv(OUTDIR / "gorilla_hor_score_cluster_summary.tsv", summaries)

    primary = [
        row
        for row in summaries
        if row["identity_threshold_pct"] == 90 and row["scope"] == "all_arrays"
    ]
    print("PRIMARY_90_PERCENT_IDENTITY")
    print(
        "cluster\tunits\tarrays\tpooled_mean\tmedian\tpct_nonzero\t"
        "array_balanced_mean\tarray_ci_lo\tarray_ci_hi\tcross_cluster_partner_pct"
    )
    for row in primary:
        print(
            f'{row["cluster"]}\t{row["units"]}\t{row["arrays"]}\t'
            f'{row["pooled_mean_score_pct"]:.6f}\t{row["median_score_pct"]:.6f}\t'
            f'{row["pct_units_nonzero"]:.6f}\t'
            f'{row["array_balanced_mean_score_pct"]:.6f}\t'
            f'{row["array_balanced_bootstrap95_lo"]:.6f}\t'
            f'{row["array_balanced_bootstrap95_hi"]:.6f}\t'
            f'{row["cross_cluster_partner_pct"]:.6f}'
        )

    print(
        f"QC\tarrays={len(array_results)}\tunits={len(unit_results)}\t"
        f"arrays_n_ge_4={sum(row['units'] >= 4 for row in array_results)}\t"
        f"seconds={time.time() - start_time:.1f}"
    )


if __name__ == "__main__":
    main()
