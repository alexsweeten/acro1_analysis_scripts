#!/usr/bin/env python3
"""Compute direct H1-H3 localized L1PA10 divergence profiles.

Each profile is calculated from ACREs assigned to the same NTRiage cluster
within a parent array.  Array-cluster groups require at least three ACREs.
Profiles are then averaged within individuals before the cohort mean and 95%
confidence interval are calculated, matching the analysis unit in Figure 6.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import warnings

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
PIPE = Path(os.environ.get("ACRO1_PROMOTER_PIPELINE_DATA", ROOT / "inputs" / "promoter_purifying_2794"))
TABLES = ROOT / "tables"
MAPPINGS = TABLES / "hmm_to_current_H_cluster_mapping.tsv.gz"
INPUT_FASTA = PIPE / "input_fasta"
PER_ARRAY = PIPE / "per_array"

sys.path.insert(0, str(ROOT.parent / "pipeline"))
import pp_common as C  # noqa: E402
import pp_core as K  # noqa: E402


MIN_GROUP_SIZE = 3
FRAGMENT_BASES = 320
PROFILE_START_300 = C.PROMOTER_300[0] - 50
PROFILE_END_300 = C.PROMOTER_300[1] + 50
DISPLAY_SHIFT = C.PROMOTER_300[0] - C.PROMOTER_FULL[0]
GROUP_ORDER = [1, 2, 3]


def extract_queries(
    mapping: pd.DataFrame,
) -> tuple[pd.DataFrame, list[tuple[str, str]], dict[str, str]]:
    """Load saved anchor alignments or extract prefixes that still need mapping."""
    metadata: list[dict[str, object]] = []
    queries: list[tuple[str, str]] = []
    prealigned: dict[str, str] = {}
    query_index = 0

    for safe_id, frame in mapping.groupby("safe_id", sort=False):
        array_dir = PER_ARRAY / safe_id
        aligned_path = array_dir / "l1pa10_anchors.aligned.fa"
        qc_path = array_dir / "units_qc.tsv"
        saved_alignment: dict[str, str] | None = None
        unit_to_anchor: dict[str, str] = {}
        array_sequence = ""
        if aligned_path.exists() and qc_path.exists():
            saved_alignment = {
                name.split()[0]: sequence.upper()
                for name, sequence in C.read_fasta(aligned_path)
            }
            qc = pd.read_csv(qc_path, sep="\t", usecols=["unit_id", "start_anchor_id"])
            unit_to_anchor = dict(zip(qc.unit_id.astype(str), qc.start_anchor_id.astype(str)))
        else:
            fasta = INPUT_FASTA / f"{safe_id}.fa"
            records = C.read_fasta(fasta)
            if len(records) != 1:
                raise ValueError(f"Expected one record in {fasta}; found {len(records)}")
            array_sequence = records[0][1].upper()

        for row in frame.itertuples(index=False):
            query_index += 1
            query_id = f"Q{query_index:08d}"
            unit_id = str(row.unit_id)
            if saved_alignment is not None:
                anchor_id = unit_to_anchor.get(unit_id)
                if anchor_id is None or anchor_id not in saved_alignment:
                    raise ValueError(f"Saved anchor missing for {safe_id}/{unit_id}")
                prealigned[query_id] = saved_alignment[anchor_id]
            else:
                unit = array_sequence[int(row.hmm_start) - 1 : int(row.hmm_end)]
                if unit_id.endswith("m"):
                    unit = C.revcomp(unit)
                fragment = unit[:FRAGMENT_BASES]
                if len(fragment) < PROFILE_END_300:
                    raise ValueError(
                        f"Unit {safe_id}/{unit_id} is too short for the localized profile"
                    )
                queries.append((query_id, fragment))
            metadata.append(
                {
                    "query_id": query_id,
                    "safe_id": safe_id,
                    "source_name": row.source_name,
                    "individual_id": row.individual_id,
                    "cluster": int(row.cluster),
                    "unit_id": unit_id,
                }
            )
    return pd.DataFrame(metadata), queries, prealigned


def align_batch(batch: list[tuple[str, str]]) -> dict[str, str]:
    """Map one batch of unit prefixes to the fixed 300-column L1PA10 MSA."""
    wanted = {name for name, _ in batch}
    with tempfile.TemporaryDirectory(prefix="figure6_l1pa10_") as temp:
        work = Path(temp)
        query_fasta = work / "queries.fa"
        output_fasta = work / "queries_plus_reference.fa"
        error_log = work / "mafft.stderr"
        C.write_fasta(query_fasta, batch)
        with open(output_fasta, "w") as out, open(error_log, "w") as err:
            subprocess.run(
                [
                    "mafft",
                    "--quiet",
                    "--thread",
                    "1",
                    "--addfragments",
                    str(query_fasta),
                    "--keeplength",
                    str(C.REF_ALN300),
                ],
                stdout=out,
                stderr=err,
                check=True,
            )
        aligned = {
            name.split()[0]: sequence.upper()
            for name, sequence in C.read_fasta(output_fasta)
            if name.split()[0] in wanted
        }
    if set(aligned) != wanted:
        missing = sorted(wanted - set(aligned))
        raise ValueError(f"MAFFT omitted {len(missing)} queries; example: {missing[:1]}")
    bad = [name for name, sequence in aligned.items() if len(sequence) != C.ALN_COLS]
    if bad:
        raise ValueError(f"MAFFT returned non-{C.ALN_COLS}-column sequence: {bad[0]}")
    return aligned


def group_profiles(metadata: pd.DataFrame, aligned: dict[str, str]) -> tuple[pd.DataFrame, np.ndarray]:
    """Calculate one callable-masked rolling profile per array-cluster group."""
    lo = PROFILE_START_300 - 1
    hi = PROFILE_END_300
    rows: list[dict[str, object]] = []
    profiles: list[np.ndarray] = []

    for (safe_id, cluster), frame in metadata.groupby(["safe_id", "cluster"], sort=False):
        sequences = [aligned[query_id] for query_id in frame.query_id]
        encoded = K.encode(sequences)
        stats = K.col_stats(encoded)
        rolling = K.rolling_pooled(stats["diffs"], stats["pairs"], C.ROLL_DIV // 2)
        profiles.append(rolling[lo:hi])
        first = frame.iloc[0]
        rows.append(
            {
                "safe_id": safe_id,
                "source_name": first.source_name,
                "individual_id": first.individual_id,
                "cluster": int(cluster),
                "n_ACREs": len(frame),
            }
        )
    return pd.DataFrame(rows), np.vstack(profiles)


def summarize_profiles(groups: pd.DataFrame, profiles: np.ndarray) -> pd.DataFrame:
    """Average arrays within individuals, then summarize equally weighted individuals."""
    anchor_columns = np.arange(PROFILE_START_300, PROFILE_END_300 + 1)
    display_columns = anchor_columns - DISPLAY_SHIFT
    output: list[dict[str, object]] = []

    for cluster in GROUP_ORDER:
        mask = groups.cluster.eq(cluster).to_numpy()
        cluster_profiles = profiles[mask]
        cluster_people = groups.loc[mask, "individual_id"].to_numpy(str)
        people = sorted(set(cluster_people))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            individual_profiles = np.vstack(
                [
                    np.nanmean(cluster_profiles[cluster_people == person], axis=0)
                    for person in people
                ]
            )
            mean = np.nanmean(individual_profiles, axis=0)
            n_callable = np.sum(np.isfinite(individual_profiles), axis=0)
            se = np.nanstd(individual_profiles, axis=0, ddof=1) / np.sqrt(n_callable)
        low = np.maximum(0.0, mean - 1.96 * se)
        high = mean + 1.96 * se

        for index, display_column in enumerate(display_columns):
            output.append(
                {
                    "display_alignment_column": int(display_column),
                    "anchor_alignment_column": int(anchor_columns[index]),
                    "cluster": f"H{cluster}",
                    "array_cluster_groups": int(mask.sum()),
                    "individuals": len(people),
                    "individuals_callable": int(n_callable[index]),
                    "individual_balanced_mean_rolling11_divergence": mean[index],
                    "ci95_low": low[index],
                    "ci95_high": high[index],
                }
            )
    return pd.DataFrame(output)


def interval_summary(profile: pd.DataFrame) -> pd.DataFrame:
    intervals = {
        "left_50bp_flank": (C.PROMOTER_FULL[0] - 50, C.PROMOTER_FULL[0] - 1),
        "candidate_50bp_interval": C.PROMOTER_FULL,
        "right_50bp_flank": (C.PROMOTER_FULL[1] + 1, C.PROMOTER_FULL[1] + 50),
    }
    rows = []
    metric = "individual_balanced_mean_rolling11_divergence"
    for cluster in ("H1", "H2", "H3"):
        subset = profile[profile.cluster.eq(cluster)]
        values = {}
        for label, (start, end) in intervals.items():
            values[label] = subset.loc[
                subset.display_alignment_column.between(start, end), metric
            ].mean()
        rows.append({"cluster": cluster, **values})
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--batch-size", type=int, default=2_000)
    args = parser.parse_args()
    TABLES.mkdir(parents=True, exist_ok=True)

    mapping = pd.read_csv(MAPPINGS, sep="\t")
    mapping = mapping[mapping.mapping_pass].copy()
    sizes = mapping.groupby(["safe_id", "cluster"]).size()
    retained = sizes[sizes >= MIN_GROUP_SIZE].index
    keep = pd.MultiIndex.from_frame(mapping[["safe_id", "cluster"]]).isin(retained)
    mapping = mapping.loc[keep].copy()
    metadata, queries, aligned = extract_queries(mapping)
    print(
        f"Prepared {len(metadata):,} ACREs from "
        f"{metadata[['safe_id', 'cluster']].drop_duplicates().shape[0]:,} array-cluster groups; "
        f"reused {len(aligned):,} saved alignments and will map {len(queries):,} prefixes",
        flush=True,
    )

    batches = [queries[i : i + args.batch_size] for i in range(0, len(queries), args.batch_size)]
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(align_batch, batch) for batch in batches]
        for index, future in enumerate(concurrent.futures.as_completed(futures), 1):
            aligned.update(future.result())
            print(f"Aligned batch {index}/{len(batches)}", flush=True)

    groups, profiles = group_profiles(metadata, aligned)
    summary = summarize_profiles(groups, profiles)
    summary.to_csv(
        TABLES / "H1_H2_H3_ACRE_level_localized_divergence_profile.tsv",
        sep="\t",
        index=False,
    )
    groups.to_csv(
        TABLES / "H1_H2_H3_ACRE_level_profile_groups.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    interval_summary(summary).to_csv(
        TABLES / "H1_H2_H3_ACRE_level_profile_interval_summary.tsv",
        sep="\t",
        index=False,
    )
    print(groups.groupby("cluster").agg(groups=("safe_id", "size"), individuals=("individual_id", "nunique")))


if __name__ == "__main__":
    main()
