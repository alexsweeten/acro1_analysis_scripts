#!/usr/bin/env python3
"""Estimate within-array ACRE divergence with Mash distance.

Mash distance is an alignment-free estimate of the per-site mutation rate from
the Jaccard similarity of canonical 21-mers. Each complete ACRE is sketched
individually (5,000 hashes), every pair within a parent array is compared, and
same-cluster pairs are averaged to one value per array and H cluster.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import os
import statistics
import subprocess
import tempfile
from collections import OrderedDict, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


H_LABEL = {1: "H1", 2: "H2", 3: "H3"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--assignments", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--limit-arrays", type=int, default=0)
    parser.add_argument("--kmer", type=int, default=21)
    parser.add_argument("--sketch-size", type=int, default=5000)
    return parser.parse_args()


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def read_assignments(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    with open_text(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            cluster = int(row["selected_cluster"])
            if cluster not in H_LABEL:
                raise ValueError(f"Unexpected selected_cluster={cluster}")
            result[row["acre_id"]] = H_LABEL[cluster]
    return result


def read_metadata(
    path: Path, assignments: dict[str, str], limit_arrays: int
) -> tuple[dict[str, dict[str, str]], dict[str, list[str]]]:
    with open_text(path) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    selected_arrays: set[str] | None = None
    if limit_arrays:
        ordered = list(dict.fromkeys(row["array_id"] for row in rows))
        selected_arrays = set(ordered[:limit_arrays])
    unit_meta: dict[str, dict[str, str]] = {}
    array_units: dict[str, list[str]] = defaultdict(list)
    for row in rows:
        if selected_arrays is not None and row["array_id"] not in selected_arrays:
            continue
        acre_id = row["acre_id"]
        meta = {
            "acre_id": acre_id,
            "array_id": row["array_id"],
            "source_name": row["source_name"],
            "unit_in_array": row["unit_in_array"],
            "h_cluster": assignments[acre_id],
        }
        unit_meta[acre_id] = meta
        array_units[row["array_id"]].append(acre_id)
    return unit_meta, dict(array_units)


class HandleCache:
    def __init__(self, max_open: int = 128):
        self.max_open = max_open
        self.handles: OrderedDict[Path, object] = OrderedDict()

    def get(self, path: Path):
        handle = self.handles.pop(path, None)
        if handle is None:
            handle = path.open("a")
        self.handles[path] = handle
        if len(self.handles) > self.max_open:
            _, old_handle = self.handles.popitem(last=False)
            old_handle.close()
        return handle

    def close(self) -> None:
        for handle in self.handles.values():
            handle.close()
        self.handles.clear()


def split_fasta_by_array(
    fasta: Path,
    unit_meta: dict[str, dict[str, str]],
    temp_root: Path,
) -> dict[str, Path]:
    array_paths: dict[str, Path] = {}
    cache = HandleCache()

    def write_record(name: str | None, chunks: list[str]) -> None:
        if name is None or name not in unit_meta:
            return
        array_id = unit_meta[name]["array_id"]
        path = temp_root / f"{array_id}.fa"
        array_paths[array_id] = path
        cache.get(path).write(f">{name}\n{''.join(chunks).upper()}\n")

    name: str | None = None
    sequence_chunks: list[str] = []
    with open_text(fasta) as handle:
        for line in handle:
            if line.startswith(">"):
                write_record(name, sequence_chunks)
                name = line[1:].split("|", 1)[0].split()[0]
                sequence_chunks = []
            else:
                sequence_chunks.append(line.strip())
        write_record(name, sequence_chunks)
    cache.close()
    expected = {meta["array_id"] for meta in unit_meta.values()}
    missing = expected - set(array_paths)
    if missing:
        raise ValueError(f"Missing sequence files for {len(missing)} arrays")
    return array_paths


def quantile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    if not ordered:
        return math.nan
    if len(ordered) == 1:
        return ordered[0]
    position = (len(ordered) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1 - weight) + ordered[upper] * weight


def analyse_array(
    array_id: str,
    fasta_path: Path,
    acre_ids: list[str],
    unit_meta: dict[str, dict[str, str]],
    kmer: int,
    sketch_size: int,
) -> tuple[list[dict[str, object]], list[dict[str, object]], dict[str, object]]:
    command = [
        "mash", "triangle", "-i", "-E", "-p", "1",
        "-k", str(kmer), "-s", str(sketch_size), str(fasta_path),
    ]
    completed = subprocess.run(
        command,
        check=True,
        text=True,
        capture_output=True,
        env={**os.environ, "LC_ALL": "C"},
    )
    pairs: dict[tuple[str, str], dict[str, object]] = {}
    for line in completed.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) < 5:
            continue
        name_1 = fields[0].split()[0]
        name_2 = fields[1].split()[0]
        if name_1 == name_2:
            continue
        pair = tuple(sorted((name_1, name_2)))
        shared, denominator = fields[4].split("/", 1)
        distance = float(fields[2])
        record = {
            "array_id": array_id,
            "source_name": unit_meta[name_1]["source_name"],
            "acre_id_1": pair[0],
            "acre_id_2": pair[1],
            "h_cluster_1": unit_meta[pair[0]]["h_cluster"],
            "h_cluster_2": unit_meta[pair[1]]["h_cluster"],
            "mash_distance": distance,
            "mash_divergence_pct": 100.0 * distance,
            "p_value": float(fields[3]),
            "shared_hashes": int(shared),
            "sketch_hashes": int(denominator),
            "shared_hash_fraction": int(shared) / int(denominator),
        }
        pairs[pair] = record

    expected_all_pairs = len(acre_ids) * (len(acre_ids) - 1) // 2
    if len(pairs) != expected_all_pairs:
        raise ValueError(
            f"{array_id}: Mash returned {len(pairs)} of {expected_all_pairs} pairs"
        )

    same_cluster_pairs = [
        row for row in pairs.values()
        if row["h_cluster_1"] == row["h_cluster_2"]
    ]
    cluster_units: dict[str, list[str]] = defaultdict(list)
    for acre_id in acre_ids:
        cluster_units[unit_meta[acre_id]["h_cluster"]].append(acre_id)

    array_rows: list[dict[str, object]] = []
    for cluster in ("H1", "H2", "H3"):
        units = cluster_units.get(cluster, [])
        if len(units) < 2:
            continue
        expected_pairs = len(units) * (len(units) - 1) // 2
        observed = [
            row for row in same_cluster_pairs
            if row["h_cluster_1"] == cluster
        ]
        if len(observed) != expected_pairs:
            raise ValueError(
                f"{array_id}/{cluster}: observed {len(observed)} of "
                f"{expected_pairs} expected same-cluster pairs"
            )
        divergence = [float(row["mash_divergence_pct"]) for row in observed]
        shared_fraction = [float(row["shared_hash_fraction"]) for row in observed]
        array_rows.append({
            "array_id": array_id,
            "source_name": unit_meta[units[0]]["source_name"],
            "h_cluster": cluster,
            "cluster_units": len(units),
            "pairwise_comparisons": expected_pairs,
            "mean_mash_divergence_pct": statistics.fmean(divergence),
            "median_mash_divergence_pct": statistics.median(divergence),
            "q25_pairwise_divergence_pct": quantile(divergence, 0.25),
            "q75_pairwise_divergence_pct": quantile(divergence, 0.75),
            "mean_shared_hash_fraction": statistics.fmean(shared_fraction),
        })
    qc = {
        "array_id": array_id,
        "units": len(acre_ids),
        "all_pairwise_comparisons": expected_all_pairs,
        "same_cluster_pairwise_comparisons": len(same_cluster_pairs),
    }
    return same_cluster_pairs, array_rows, qc


def make_writer(handle, fieldnames: list[str]):
    writer = csv.DictWriter(
        handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore"
    )
    writer.writeheader()
    return writer


def write_table(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = make_writer(handle, fields)
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    assignments = read_assignments(args.assignments)
    unit_meta, array_units = read_metadata(
        args.metadata, assignments, args.limit_arrays
    )
    pair_fields = [
        "array_id", "source_name", "acre_id_1", "acre_id_2",
        "h_cluster_1", "h_cluster_2", "mash_distance",
        "mash_divergence_pct", "p_value", "shared_hashes", "sketch_hashes",
        "shared_hash_fraction",
    ]
    array_fields = [
        "array_id", "source_name", "h_cluster", "cluster_units",
        "pairwise_comparisons", "mean_mash_divergence_pct",
        "median_mash_divergence_pct", "q25_pairwise_divergence_pct",
        "q75_pairwise_divergence_pct", "mean_shared_hash_fraction",
    ]
    qc_fields = [
        "array_id", "units", "all_pairwise_comparisons",
        "same_cluster_pairwise_comparisons",
    ]
    array_rows: list[dict[str, object]] = []
    qc_rows: list[dict[str, object]] = []
    pair_count = 0
    pair_path = args.out_dir / "pairwise_within_array_mash_divergence.tsv.gz"
    with tempfile.TemporaryDirectory(prefix="hprc_mash_") as temp_name:
        temp_root = Path(temp_name)
        array_paths = split_fasta_by_array(args.fasta, unit_meta, temp_root)
        with gzip.open(pair_path, "wt", newline="") as pair_handle:
            pair_writer = make_writer(pair_handle, pair_fields)
            with ThreadPoolExecutor(max_workers=args.workers) as executor:
                future_map = {
                    executor.submit(
                        analyse_array,
                        array_id,
                        array_paths[array_id],
                        acre_ids,
                        unit_meta,
                        args.kmer,
                        args.sketch_size,
                    ): array_id
                    for array_id, acre_ids in array_units.items()
                }
                completed_count = 0
                for future in as_completed(future_map):
                    pairs, arrays, qc = future.result()
                    pair_writer.writerows(pairs)
                    pair_count += len(pairs)
                    array_rows.extend(arrays)
                    qc_rows.append(qc)
                    completed_count += 1
                    if completed_count % 250 == 0:
                        print(
                            f"Completed {completed_count}/{len(future_map)} arrays",
                            flush=True,
                        )

    array_rows.sort(key=lambda row: (row["h_cluster"], row["array_id"]))
    qc_rows.sort(key=lambda row: row["array_id"])
    write_table(args.out_dir / "array_cluster_mash_divergence.tsv", array_rows, array_fields)
    write_table(args.out_dir / "mash_qc.tsv", qc_rows, qc_fields)
    parameter_rows = [
        {"parameter": "metric", "value": "Mash distance x 100"},
        {"parameter": "kmer", "value": args.kmer},
        {"parameter": "sketch_size", "value": args.sketch_size},
        {"parameter": "array_count", "value": len(array_units)},
        {"parameter": "unit_count", "value": len(unit_meta)},
        {"parameter": "same_cluster_pair_count", "value": pair_count},
        {"parameter": "cluster_label_mapping", "value": "1=H1;2=H2;3=H3"},
    ]
    write_table(
        args.out_dir / "analysis_parameters.tsv",
        parameter_rows,
        ["parameter", "value"],
    )
    print(
        f"Wrote {pair_count:,} same-cluster pairs and "
        f"{len(array_rows):,} array-cluster summaries"
    )


if __name__ == "__main__":
    main()
