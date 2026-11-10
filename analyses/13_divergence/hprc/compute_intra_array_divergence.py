#!/usr/bin/env python3
"""Compute within-array ACRE sequence divergence by HPRC H cluster.

For each parent array, complete ACRE sequences are compared all-vs-all with
minimap2. Only pairs whose two units have the same final H1-H3 label are
summarized. Arrays, rather than individual pairs, are the independent units in
the downstream cluster comparison.
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
    parser.add_argument("--min-coverage", type=float, default=0.80)
    parser.add_argument("--min-units", type=int, default=2)
    return parser.parse_args()


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def read_assignments(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            cluster = int(row["selected_cluster"])
            if cluster not in H_LABEL:
                raise ValueError(f"Unexpected selected_cluster={cluster}")
            result[row["acre_id"]] = H_LABEL[cluster]
    return result


def read_metadata(
    path: Path, assignments: dict[str, str], limit_arrays: int
) -> tuple[dict[str, dict[str, str]], dict[str, list[str]]]:
    unit_meta: dict[str, dict[str, str]] = {}
    array_units: dict[str, list[str]] = defaultdict(list)
    selected_arrays: set[str] | None = None
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    if limit_arrays:
        ordered_arrays = list(dict.fromkeys(row["array_id"] for row in rows))
        selected_arrays = set(ordered_arrays[:limit_arrays])
    for row in rows:
        if selected_arrays is not None and row["array_id"] not in selected_arrays:
            continue
        acre_id = row["acre_id"]
        if acre_id not in assignments:
            raise ValueError(f"No final cluster assignment for {acre_id}")
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
            path.parent.mkdir(parents=True, exist_ok=True)
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
        handle = cache.get(path)
        handle.write(f">{name}\n{''.join(chunks).upper()}\n")

    current_name: str | None = None
    chunks: list[str] = []
    with open_text(fasta) as handle:
        for line in handle:
            if line.startswith(">"):
                write_record(current_name, chunks)
                current_name = line[1:].split("|", 1)[0].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        write_record(current_name, chunks)
    cache.close()
    missing = set(meta["array_id"] for meta in unit_meta.values()) - set(array_paths)
    if missing:
        raise ValueError(f"Missing FASTA sequences for {len(missing)} arrays")
    return array_paths


def quantile(values: list[float], probability: float) -> float:
    if not values:
        return math.nan
    ordered = sorted(values)
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
    min_coverage: float,
) -> tuple[list[dict[str, object]], list[dict[str, object]], dict[str, object]]:
    command = [
        "minimap2", "-x", "ava-ont", "-X", "-t", "1",
        str(fasta_path), str(fasta_path),
    ]
    completed = subprocess.run(
        command, check=True, text=True, capture_output=True,
        env={**os.environ, "LC_ALL": "C"},
    )
    best: dict[tuple[str, str], dict[str, object]] = {}
    total_paf = 0
    for line in completed.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) < 12:
            continue
        total_paf += 1
        qname, tname = fields[0], fields[5]
        if qname == tname or qname not in unit_meta or tname not in unit_meta:
            continue
        pair = tuple(sorted((qname, tname)))
        qlen = int(fields[1])
        qstart = int(fields[2])
        qend = int(fields[3])
        tlen = int(fields[6])
        tstart = int(fields[7])
        tend = int(fields[8])
        nmatch = int(fields[9])
        block_length = int(fields[10])
        if block_length <= 0:
            continue
        query_coverage = (qend - qstart) / qlen
        target_coverage = (tend - tstart) / tlen
        symmetric_coverage = min(query_coverage, target_coverage)
        ani_pct = 100.0 * nmatch / block_length
        full_length_identity_pct = 100.0 * nmatch / max(qlen, tlen)
        record = {
            "array_id": array_id,
            "source_name": unit_meta[qname]["source_name"],
            "acre_id_1": pair[0],
            "acre_id_2": pair[1],
            "h_cluster_1": unit_meta[pair[0]]["h_cluster"],
            "h_cluster_2": unit_meta[pair[1]]["h_cluster"],
            "length_1": qlen if pair[0] == qname else tlen,
            "length_2": tlen if pair[1] == tname else qlen,
            "aligned_block_bp": block_length,
            "matching_bp": nmatch,
            "alignment_coverage": symmetric_coverage,
            "ani_pct": ani_pct,
            "aligned_divergence_pct": 100.0 - ani_pct,
            "full_length_identity_pct": full_length_identity_pct,
            "full_length_divergence_pct": 100.0 - full_length_identity_pct,
            "qualifies": int(symmetric_coverage >= min_coverage),
        }
        score = (symmetric_coverage, nmatch, ani_pct, block_length)
        previous = best.get(pair)
        if previous is None or score > previous["_score"]:
            record["_score"] = score
            best[pair] = record

    pair_rows: list[dict[str, object]] = []
    for record in best.values():
        record.pop("_score", None)
        if record["h_cluster_1"] == record["h_cluster_2"]:
            pair_rows.append(record)

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
            row for row in pair_rows
            if row["h_cluster_1"] == cluster
        ]
        qualifying = [row for row in observed if row["qualifies"]]
        divergences = [float(row["aligned_divergence_pct"]) for row in qualifying]
        full_divergences = [float(row["full_length_divergence_pct"]) for row in qualifying]
        coverages = [float(row["alignment_coverage"]) for row in qualifying]
        array_rows.append({
            "array_id": array_id,
            "source_name": unit_meta[units[0]]["source_name"],
            "h_cluster": cluster,
            "cluster_units": len(units),
            "expected_pairs": expected_pairs,
            "observed_pairs": len(observed),
            "qualifying_pairs": len(qualifying),
            "observed_pair_fraction": len(observed) / expected_pairs,
            "qualifying_pair_fraction": len(qualifying) / expected_pairs,
            "mean_ani_pct": statistics.fmean(
                float(row["ani_pct"]) for row in qualifying
            ) if qualifying else math.nan,
            "mean_aligned_divergence_pct": statistics.fmean(divergences) if divergences else math.nan,
            "median_aligned_divergence_pct": statistics.median(divergences) if divergences else math.nan,
            "q25_aligned_divergence_pct": quantile(divergences, 0.25),
            "q75_aligned_divergence_pct": quantile(divergences, 0.75),
            "mean_full_length_divergence_pct": statistics.fmean(full_divergences) if full_divergences else math.nan,
            "mean_alignment_coverage": statistics.fmean(coverages) if coverages else math.nan,
        })
    qc = {
        "array_id": array_id,
        "units": len(acre_ids),
        "paf_records": total_paf,
        "unique_pairs_any_cluster": len(best),
        "same_cluster_pairs": len(pair_rows),
    }
    return pair_rows, array_rows, qc


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if path.suffix == ".gz" else open
    kwargs = {"mode": "wt", "newline": ""} if path.suffix == ".gz" else {"mode": "w", "newline": ""}
    with opener(path, **kwargs) as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    assignments = read_assignments(args.assignments)
    unit_meta, array_units = read_metadata(
        args.metadata, assignments, args.limit_arrays
    )
    with tempfile.TemporaryDirectory(prefix="hprc_divergence_") as temp_name:
        temp_root = Path(temp_name)
        array_paths = split_fasta_by_array(args.fasta, unit_meta, temp_root)
        pair_rows: list[dict[str, object]] = []
        array_rows: list[dict[str, object]] = []
        qc_rows: list[dict[str, object]] = []
        with ThreadPoolExecutor(max_workers=args.workers) as executor:
            future_map = {
                executor.submit(
                    analyse_array,
                    array_id,
                    array_paths[array_id],
                    acre_ids,
                    unit_meta,
                    args.min_coverage,
                ): array_id
                for array_id, acre_ids in array_units.items()
            }
            completed_count = 0
            for future in as_completed(future_map):
                pairs, arrays, qc = future.result()
                pair_rows.extend(pairs)
                array_rows.extend(arrays)
                qc_rows.append(qc)
                completed_count += 1
                if completed_count % 250 == 0:
                    print(f"Completed {completed_count}/{len(future_map)} arrays", flush=True)

    pair_rows.sort(key=lambda row: (row["array_id"], row["h_cluster_1"], row["acre_id_1"], row["acre_id_2"]))
    array_rows.sort(key=lambda row: (row["h_cluster"], row["array_id"]))
    qc_rows.sort(key=lambda row: row["array_id"])
    pair_fields = [
        "array_id", "source_name", "acre_id_1", "acre_id_2",
        "h_cluster_1", "h_cluster_2", "length_1", "length_2",
        "aligned_block_bp", "matching_bp", "alignment_coverage", "ani_pct",
        "aligned_divergence_pct", "full_length_identity_pct",
        "full_length_divergence_pct", "qualifies",
    ]
    array_fields = [
        "array_id", "source_name", "h_cluster", "cluster_units",
        "expected_pairs", "observed_pairs", "qualifying_pairs",
        "observed_pair_fraction", "qualifying_pair_fraction", "mean_ani_pct",
        "mean_aligned_divergence_pct", "median_aligned_divergence_pct",
        "q25_aligned_divergence_pct", "q75_aligned_divergence_pct",
        "mean_full_length_divergence_pct", "mean_alignment_coverage",
    ]
    qc_fields = [
        "array_id", "units", "paf_records", "unique_pairs_any_cluster",
        "same_cluster_pairs",
    ]
    write_tsv(args.out_dir / "pairwise_within_array_divergence.tsv.gz", pair_rows, pair_fields)
    write_tsv(args.out_dir / "array_cluster_divergence.tsv", array_rows, array_fields)
    write_tsv(args.out_dir / "alignment_qc.tsv", qc_rows, qc_fields)
    print(
        f"Wrote {len(pair_rows):,} same-cluster pairs and "
        f"{len(array_rows):,} array-cluster summaries"
    )


if __name__ == "__main__":
    main()
