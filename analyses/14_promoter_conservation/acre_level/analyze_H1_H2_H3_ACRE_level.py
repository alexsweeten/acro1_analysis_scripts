#!/usr/bin/env python3
"""Recompute fixed-promoter divergence within current H1-H3 ACRE classes.

Individual HMM-defined complete ACREs are mapped to the current HPRC k=3
NTRiage units.  Divergence is then calculated separately among same-H-cluster
ACREs within each parent array.  The inferential observation is an
array-by-cluster group, not an individual ACRE.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import importlib.util
import json
import os
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent
DATA_ROOT = Path(os.environ.get("ACRO1_WORKSPACE", ROOT / "inputs"))
BASE_PATH = Path(os.environ.get("ACRO1_PROMOTER_CLUSTER_BASE", ROOT / "base_promoter_by_cluster.py"))
CURRENT_UNITS = Path(os.environ.get("ACRO1_HPRC_UNITS", DATA_ROOT / "hprc.clusters.tsv"))
CURRENT_ASSIGNMENTS = Path(os.environ.get("ACRO1_HPRC_ASSIGNMENTS", DATA_ROOT / "assignments.tsv.gz"))


def load_base():
    spec = importlib.util.spec_from_file_location("promoter_cluster_base_h3", BASE_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not import {BASE_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    module.OUT = ROOT
    module.TABLES = ROOT / "tables"
    module.FIGURES = ROOT / "figures"
    module.LOGS = ROOT / "logs"
    module.CACHE = ROOT / "cache"
    return module


def current_native_table() -> pd.DataFrame:
    units = pd.read_csv(
        CURRENT_UNITS,
        sep="\t",
        usecols=["acre_id", "source_name", "unit_in_array", "start", "end"],
    )
    assignments = pd.read_csv(
        CURRENT_ASSIGNMENTS,
        sep="\t",
        usecols=["acre_id", "selected_cluster"],
    )
    native = units.merge(assignments, on="acre_id", how="inner", validate="one_to_one")
    if len(native) != len(units):
        raise ValueError(f"Mapped {len(native):,} of {len(units):,} current NTRiage units")
    if set(native.selected_cluster.unique()) != {1, 2, 3}:
        raise ValueError("Current assignments are not exactly k=3")
    return native.rename(columns={"selected_cluster": "cluster"})


def relabel(frame: pd.DataFrame) -> pd.DataFrame:
    output = frame.copy()
    for column in ("cluster", "cluster_a", "cluster_b"):
        if column in output:
            output[column] = output[column].astype(str).str.replace(r"^C", "H", regex=True)
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=12)
    args = parser.parse_args()
    base = load_base()
    for directory in (base.TABLES, base.FIGURES, base.LOGS, base.CACHE):
        directory.mkdir(parents=True, exist_ok=True)

    manifest = pd.read_csv(base.PIPE / "manifest.validated.tsv", sep="\t")
    native = current_native_table()
    native_groups = {key: value for key, value in native.groupby("source_name", sort=False)}
    tasks = []
    for _, row in manifest.iterrows():
        source = row.array_id
        if source not in native_groups:
            continue
        tasks.append(
            {
                "source_name": source,
                "safe_id": row.safe_id,
                "individual_id": row.individual_id,
                "haplotype": row.haplotype,
                "has_recombination": row.has_recombination,
                "native_units": native_groups[source][
                    ["unit_in_array", "start", "end", "cluster"]
                ].to_dict("records"),
            }
        )

    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(base.analyze_array, task) for task in tasks]
        for index, future in enumerate(concurrent.futures.as_completed(futures), 1):
            results.append(future.result())
            if index % 200 == 0 or index == len(futures):
                errors = sum(row["status"] == "error" for row in results)
                print(f"completed {index}/{len(futures)}; errors={errors}", flush=True)

    errors = pd.DataFrame(
        [
            {key: value for key, value in row.items() if key not in ("groups", "mappings")}
            for row in results
            if row["status"] == "error"
        ],
        columns=[
            "source_name", "safe_id", "individual_id", "haplotype",
            "has_recombination", "status", "error",
        ],
    )
    errors.to_csv(base.LOGS / "errors.tsv", sep="\t", index=False)
    if len(errors):
        raise SystemExit(f"{len(errors)} arrays failed; see {base.LOGS / 'errors.tsv'}")

    mappings = pd.DataFrame([record for row in results for record in row["mappings"]])
    groups = pd.DataFrame([record for row in results for record in row["groups"]])
    groups["distal_backbone_absolute_depletion"] = (
        groups["distal_backbone_divergence"] - groups["raw_promoter_divergence"]
    )
    groups["l1pa10_nonpromoter_absolute_depletion"] = (
        groups["l1pa10_nonpromoter_divergence"] - groups["raw_promoter_divergence"]
    )
    mappings.to_csv(
        base.TABLES / "hmm_to_current_H_cluster_mapping.tsv.gz",
        sep="\t", index=False, compression="gzip",
    )
    groups.to_csv(
        base.TABLES / "array_H_cluster_promoter_metrics.tsv.gz",
        sep="\t", index=False, compression="gzip",
    )

    group_counts = groups.groupby("source_name").cluster.nunique()
    mixed_ids = set(group_counts[group_counts > 1].index)
    single_ids = set(group_counts[group_counts == 1].index)
    cohorts = {
        "all_arrays": groups,
        "mixed_H_cluster_arrays": groups[groups.source_name.isin(mixed_ids)].copy(),
        "single_H_cluster_arrays": groups[groups.source_name.isin(single_ids)].copy(),
    }
    metrics = [
        "raw_promoter_divergence",
        "distal_backbone_divergence",
        "l1pa10_nonpromoter_divergence",
        "distal_backbone_absolute_depletion",
        "l1pa10_nonpromoter_absolute_depletion",
        "promoter_depletion_scaled_by_distal_backbone_divergence",
        "promoter_depletion_scaled_by_l1pa10_nonpromoter_divergence",
        "local_effect",
        "backbone_effect",
    ]
    summaries = []
    pairwise = []
    one_sample = []
    estimability = []
    correlations = []
    for name, frame in cohorts.items():
        summaries.append(base.cluster_summary(frame, name, metrics))
        pairwise.append(base.pairwise_tests(frame, name, metrics))
        one_sample.append(
            base.one_sample_effect_tests(frame, name, ["local_effect", "backbone_effect"])
        )
        estimability.append(base.estimability_summary(frame, name))
        correlations.append(base.correlation_summary(frame, name))

    relabel(pd.concat(summaries, ignore_index=True)).to_csv(
        base.TABLES / "H_cluster_summary_statistics.tsv", sep="\t", index=False
    )
    relabel(pd.concat(pairwise, ignore_index=True)).to_csv(
        base.TABLES / "within_array_H_cluster_pairwise_tests.tsv", sep="\t", index=False
    )
    relabel(pd.concat(one_sample, ignore_index=True)).to_csv(
        base.TABLES / "H_cluster_effect_vs_zero_tests.tsv", sep="\t", index=False
    )
    relabel(pd.concat(estimability, ignore_index=True)).to_csv(
        base.TABLES / "H_cluster_control_estimability.tsv", sep="\t", index=False
    )
    relabel(pd.concat(correlations, ignore_index=True)).to_csv(
        base.TABLES / "H_cluster_promoter_background_correlations.tsv", sep="\t", index=False
    )

    metadata = {
        "arrays_processed": len(results),
        "current_ntriage_units": len(native),
        "hmm_units": len(mappings),
        "mapped_hmm_units": int(mappings.mapping_pass.sum()),
        "mapping_rate": float(mappings.mapping_pass.mean()),
        "array_H_cluster_groups_min_3_units": len(groups),
        "arrays_with_multiple_analyzable_H_clusters": len(mixed_ids),
        "arrays_with_one_analyzable_H_cluster": len(single_ids),
        "minimum_units_per_array_cluster_group": base.MIN_CLUSTER_UNITS,
        "cluster_mapping": "selected_cluster 1=H1; 2=H2; 3=H3",
        "fixed_promoter_full_alignment_columns_1based_inclusive": list(base.C.PROMOTER_FULL),
        "signflip_draws": base.N_SIGNFLIP,
        "bootstrap_draws": base.N_BOOT,
        "seed": base.SEED,
    }
    (ROOT / "analysis_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
