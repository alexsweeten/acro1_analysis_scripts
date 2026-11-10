#!/usr/bin/env python3
"""Build compact run- and cluster-level summaries from completed results."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=PROJECT)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument(
        "--datasets",
        nargs="*",
        help="Dataset IDs to summarize (default: all rows in the manifest)",
    )
    parser.add_argument(
        "--combined-output-dir",
        type=Path,
        default=Path("results/combined"),
        help=(
            "Directory containing aggregate analysis tables and receiving compact "
            "summaries, relative to --root unless absolute (default: results/combined)"
        ),
    )
    return parser.parse_args()


def resolved(value: object, base: Path) -> Path:
    path = Path(str(value)).expanduser()
    return path if path.is_absolute() else (base / path).resolve()


def require_columns(frame: pd.DataFrame, columns: list[str], source: Path) -> None:
    missing = sorted(set(columns) - set(frame.columns))
    if missing:
        raise ValueError(f"{source}: missing columns {missing}")


def main() -> None:
    args = parse_args()
    root = args.root.resolve()
    manifest_path = (args.manifest or root / "inputs/dataset_manifest.tsv").resolve()
    manifest = pd.read_csv(manifest_path, sep="\t")
    require_columns(
        manifest,
        [
            "dataset_id",
            "display_name",
            "ntriage_prefix",
            "n_arrays",
            "analyzed_bp",
            "min_period_bp",
            "max_period_bp",
            "period_target_bp",
        ],
        manifest_path,
    )
    if manifest.dataset_id.duplicated().any():
        raise ValueError(f"{manifest_path}: duplicate dataset_id")
    manifest["dataset_id"] = manifest.dataset_id.astype(str)
    if args.datasets:
        unknown = sorted(set(args.datasets) - set(manifest.dataset_id))
        if unknown:
            raise ValueError(f"Unknown --datasets: {', '.join(unknown)}")
        order = {dataset: index for index, dataset in enumerate(args.datasets)}
        manifest = manifest.loc[manifest.dataset_id.isin(args.datasets)].copy()
        manifest = manifest.sort_values(
            "dataset_id", key=lambda values: values.map(order)
        )

    combined_output = resolved(args.combined_output_dir, root)
    combined_path = combined_output / "all_cluster_summaries.tsv"
    clusters = pd.read_csv(combined_path, sep="\t")
    required_cluster_columns = [
        "dataset_id",
        "cluster",
        "cluster_label",
        "selected_k",
        "units",
        "arrays",
        "fraction_of_units",
        "mean_unit_length_bp",
        "median_unit_length_bp",
        "min_unit_length_bp",
        "max_unit_length_bp",
        "array_weighted_mean_unit_length_bp",
        "mean_acro1_copy_number",
        "median_acro1_copy_number",
        "array_weighted_mean_acro1_copy_number",
        "mean_classified_acro1_copy_number",
        "raw_acro1_anchor_total",
        "classified_acro1_total",
        "fraction_raw_acro1_anchors_classified",
        "subfamily_1_count",
        "subfamily_2_count",
        "subfamily_3_count",
        "fraction_subfamily_1_pooled",
        "fraction_subfamily_2_pooled",
        "fraction_subfamily_3_pooled",
        "mean_raw_acro1_anchor_union_bp",
        "mean_classified_acro1_union_bp",
        "mean_broad_acro1_annotation_bp",
        "mean_explicit_backbone_bp",
        "mean_unclassified_bp",
        "mean_non_anchor_remainder_bp",
    ]
    require_columns(clusters, required_cluster_columns, combined_path)

    manifest_ids = manifest.dataset_id.astype(str).tolist()
    if set(manifest_ids) != set(clusters.dataset_id.astype(str)):
        raise ValueError("Manifest and combined cluster-summary datasets differ")

    run_rows: list[dict[str, object]] = []
    long_rows: list[pd.DataFrame] = []
    for row in manifest.to_dict("records"):
        dataset = str(row["dataset_id"])
        prefix = resolved(row["ntriage_prefix"], manifest_path.parent)
        native_path = Path(f"{prefix}.clusters.tsv")
        native = pd.read_csv(
            native_path, sep="\t", usecols=["acre_id", "array_id", "length"]
        )
        native["length"] = pd.to_numeric(native.length, errors="raise").astype(int)
        if native.acre_id.duplicated().any():
            raise ValueError(f"{native_path}: duplicate acre_id")

        decision_path = root / "results" / dataset / "k_selection/selection_decision.json"
        decision = json.loads(decision_path.read_text())
        support_path = root / "results" / dataset / "k_selection/final_cluster_support.tsv"
        support = pd.read_csv(support_path, sep="\t")

        subset = clusters.loc[clusters.dataset_id.astype(str).eq(dataset)].copy()
        summary_units = int(subset.units.sum())
        if summary_units != len(native):
            raise ValueError(
                f"{dataset}: native units={len(native)} but cluster summaries={summary_units}"
            )
        if int(subset.selected_k.nunique()) != 1:
            raise ValueError(f"{dataset}: inconsistent selected_k values")
        authoritative_k = int(decision["authoritative_k"])
        if int(subset.selected_k.iloc[0]) != authoritative_k:
            raise ValueError(f"{dataset}: decision and cluster summary K differ")

        raw_total = int(subset.raw_acro1_anchor_total.sum())
        classified_total = int(subset.classified_acro1_total.sum())
        run_rows.append(
            {
                "dataset_id": dataset,
                "display_name": row["display_name"],
                "input_arrays": int(row["n_arrays"]),
                "arrays_with_complete_units": int(native.array_id.nunique()),
                "analyzed_input_bp": int(row["analyzed_bp"]),
                "complete_units": len(native),
                "mean_unit_length_bp": float(native.length.mean()),
                "min_unit_length_bp": int(native.length.min()),
                "max_unit_length_bp": int(native.length.max()),
                "units_ge_7000bp": int(native.length.ge(7000).sum()),
                "fraction_units_ge_7000bp": float(native.length.ge(7000).mean()),
                "units_ge_8000bp": int(native.length.ge(8000).sum()),
                "fraction_units_ge_8000bp": float(native.length.ge(8000).mean()),
                "min_period_bp": int(row["min_period_bp"]),
                "max_period_bp": int(row["max_period_bp"]),
                "period_target_bp": int(row["period_target_bp"]),
                "preferred_view": decision["preferred_view"],
                "preferred_k": int(decision["preferred_k"]),
                "authoritative_view": decision["authoritative_view"],
                "authoritative_k": authoritative_k,
                "fallback_used": bool(decision["fallback_used"]),
                "weakest_cluster_mean_agreement": float(
                    decision["weakest_authoritative_cluster_mean_agreement"]
                ),
                "minimum_cluster_units": int(support.units.min()),
                "minimum_cluster_arrays": int(support.arrays.min()),
                "raw_acro1_anchors_assigned": raw_total,
                "classified_s1_s2_s3_anchors": classified_total,
                "fraction_raw_anchors_classified": (
                    classified_total / raw_total if raw_total else 0.0
                ),
            }
        )

        unit_path = root / "results" / dataset / "cluster_analysis/unit_metrics.tsv.gz"
        units = pd.read_csv(
            unit_path,
            sep="\t",
            usecols=["cluster", "cluster_label", "length"],
        )
        units["units_ge_7000bp"] = units.length.ge(7000).astype(int)
        units["units_ge_8000bp"] = units.length.ge(8000).astype(int)
        long = (
            units.groupby(["cluster", "cluster_label"], as_index=False)
            .agg(
                units_from_unit_table=("length", "size"),
                units_ge_7000bp=("units_ge_7000bp", "sum"),
                units_ge_8000bp=("units_ge_8000bp", "sum"),
            )
        )
        long["dataset_id"] = dataset
        long_rows.append(long)

    run_summary = pd.DataFrame(run_rows)
    long_summary = pd.concat(long_rows, ignore_index=True)
    core = clusters.merge(
        long_summary,
        on=["dataset_id", "cluster", "cluster_label"],
        how="left",
        validate="one_to_one",
    )
    if not core.units.eq(core.units_from_unit_table).all():
        raise ValueError("Cluster unit counts disagree with unit_metrics")
    core["fraction_units_ge_7000bp"] = core.units_ge_7000bp / core.units
    core["fraction_units_ge_8000bp"] = core.units_ge_8000bp / core.units

    core_columns = [
        "dataset_id",
        "cluster",
        "cluster_label",
        "selected_k",
        "units",
        "fraction_of_units",
        "arrays",
        "mean_unit_length_bp",
        "median_unit_length_bp",
        "min_unit_length_bp",
        "max_unit_length_bp",
        "array_weighted_mean_unit_length_bp",
        "units_ge_7000bp",
        "fraction_units_ge_7000bp",
        "units_ge_8000bp",
        "fraction_units_ge_8000bp",
        "mean_acro1_copy_number",
        "median_acro1_copy_number",
        "array_weighted_mean_acro1_copy_number",
        "mean_classified_acro1_copy_number",
        "raw_acro1_anchor_total",
        "classified_acro1_total",
        "fraction_raw_acro1_anchors_classified",
        "subfamily_1_count",
        "subfamily_2_count",
        "subfamily_3_count",
        "fraction_subfamily_1_pooled",
        "fraction_subfamily_2_pooled",
        "fraction_subfamily_3_pooled",
        "mean_raw_acro1_anchor_union_bp",
        "mean_classified_acro1_union_bp",
        "mean_broad_acro1_annotation_bp",
        "mean_explicit_backbone_bp",
        "mean_unclassified_bp",
        "mean_non_anchor_remainder_bp",
    ]
    core = core[core_columns]

    output = combined_output
    output.mkdir(parents=True, exist_ok=True)
    run_path = output / "experiment_run_summary.tsv"
    core_path = output / "cluster_core_characteristics.tsv"
    run_summary.to_csv(run_path, sep="\t", index=False)
    core.to_csv(core_path, sep="\t", index=False)
    print(f"Wrote {len(run_summary)} run rows to {run_path}")
    print(f"Wrote {len(core)} cluster rows to {core_path}")


if __name__ == "__main__":
    main()
