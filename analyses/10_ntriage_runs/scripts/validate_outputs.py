#!/usr/bin/env python3
"""Validate completeness and internal consistency of the Nov-3 experiment.

This is an end-of-pipeline validator.  It checks input/run provenance, native
NTRIAGE tables, selected cluster assignments, 6--9 kb unit boundaries, ACRO1
and annotation coordinates, and the derived cluster-analysis tables.  A
machine-readable report is always written before a non-zero exit.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd

from analyze_clusters import (
    DEFAULT_MANIFEST,
    DEFAULT_ROOT,
    REQUIRED_DATASETS,
    dataset_column,
    load_manifest,
    normalize_assignments,
    output_prefix,
    selected_assignments_path,
)


NATIVE_SUFFIXES = [
    ".arrays.tsv",
    ".clusters.tsv",
    ".features.tsv",
    ".acres.bed",
    ".acres.fasta",
    ".annotations.bed",
    ".acro1_subfamilies.bed",
    ".acro1_composition.tsv",
]
ANALYSIS_FILES = [
    "unit_metrics.tsv.gz",
    "array_cluster_metrics.tsv.gz",
    "cluster_summary.tsv",
    "cluster_acro1_allele_composition.tsv",
    "cluster_unit_length_bins.tsv",
    "cluster_component_bp.tsv",
    "cluster_array_bootstrap_intervals.tsv",
    "cluster_descriptive_tests.tsv",
    "species_cluster_summary.tsv",
    "analysis_qc.json",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--datasets", nargs="*")
    parser.add_argument("--joint-id", default="joint_all")
    parser.add_argument(
        "--combined-output-dir",
        type=Path,
        default=Path("results/combined"),
        help=(
            "Directory containing aggregate analysis tables, relative to --root "
            "unless absolute (default: results/combined)"
        ),
    )
    parser.add_argument(
        "--validation-output-dir",
        type=Path,
        default=Path("validation"),
        help=(
            "Directory receiving validation_report.tsv and validation_summary.json, "
            "relative to --root unless absolute (default: validation)"
        ),
    )
    parser.add_argument(
        "--allow-partial",
        action="store_true",
        help="Do not require the canonical full-dataset manifest",
    )
    parser.add_argument(
        "--no-fail",
        action="store_true",
        help="Write errors but return exit status zero (useful while developing)",
    )
    return parser.parse_args()


class Report:
    def __init__(self) -> None:
        self.rows: list[dict[str, object]] = []

    def add(
        self,
        dataset: str,
        check: str,
        passed: bool,
        detail: str,
        *,
        severity: str = "error",
        observed: object = "",
        expected: object = "",
    ) -> None:
        status = "PASS" if passed else ("WARN" if severity == "warning" else "FAIL")
        self.rows.append(
            {
                "dataset_id": dataset,
                "check": check,
                "status": status,
                "severity_if_failed": severity,
                "observed": observed,
                "expected": expected,
                "detail": detail,
            }
        )

    def exception(self, dataset: str, check: str, error: Exception) -> None:
        self.add(dataset, check, False, f"{type(error).__name__}: {error}")

    def frame(self) -> pd.DataFrame:
        return pd.DataFrame(self.rows)


def check_exists(report: Report, dataset: str, check: str, path: Path) -> bool:
    exists = path.is_file() and path.stat().st_size > 0
    report.add(
        dataset,
        check,
        exists,
        str(path),
        observed=path.stat().st_size if exists else 0,
        expected=">0 bytes",
    )
    return exists


def root_relative(path: Path, root: Path) -> Path:
    path = path.expanduser()
    return path.resolve() if path.is_absolute() else (root / path).resolve()


def compare_id_sets(
    report: Report,
    dataset: str,
    check: str,
    reference: Sequence[object],
    observed: Sequence[object],
) -> None:
    left, right = set(map(str, reference)), set(map(str, observed))
    report.add(
        dataset,
        check,
        left == right,
        f"only_reference={len(left-right)}, only_observed={len(right-left)}",
        observed=len(right),
        expected=len(left),
    )


def validate_combined_outputs(
    report: Report,
    combined: Path,
    manifest: pd.DataFrame,
    id_col: str,
    manifest_path: Path,
    joint_id: str,
    root: Path,
) -> None:
    """Validate aggregate membership and the parameters that produced it."""
    expected_ids = manifest[id_col].astype(str).tolist()
    paths = {
        filename: combined / filename
        for filename in (
            "all_cluster_summaries.tsv",
            "all_species_cluster_summaries.tsv",
            "individual_vs_joint_cluster_crosswalk.tsv",
            "individual_vs_joint_unit_overlap.tsv",
            "analysis_parameters.json",
        )
    }
    present = {
        filename: check_exists(
            report, "GLOBAL", f"combined_file:{filename}", path
        )
        for filename, path in paths.items()
    }

    for filename in (
        "all_cluster_summaries.tsv",
        "all_species_cluster_summaries.tsv",
    ):
        if not present[filename]:
            continue
        try:
            table = pd.read_csv(paths[filename], sep="\t", dtype=str)
            if "dataset_id" not in table:
                report.add(
                    "GLOBAL",
                    f"combined_dataset_ids:{filename}",
                    False,
                    "aggregate table lacks dataset_id",
                    observed="missing",
                    expected="dataset_id column",
                )
            else:
                compare_id_sets(
                    report,
                    "GLOBAL",
                    f"combined_dataset_ids:{filename}",
                    expected_ids,
                    table.dataset_id.dropna().tolist(),
                )
        except Exception as error:
            report.exception("GLOBAL", f"combined_dataset_ids:{filename}", error)

    parameters_name = "analysis_parameters.json"
    if not present[parameters_name]:
        return
    try:
        parameters = json.loads(paths[parameters_name].read_text())
    except Exception as error:
        report.exception("GLOBAL", "combined_analysis_parameters_json", error)
        return

    recorded_manifest = parameters.get("manifest")
    if recorded_manifest is None:
        manifest_matches = False
        recorded_manifest_resolved = "missing"
    else:
        recorded_manifest_resolved = str(
            root_relative(Path(str(recorded_manifest)), root)
        )
        manifest_matches = Path(recorded_manifest_resolved) == manifest_path
    report.add(
        "GLOBAL",
        "combined_parameters_manifest",
        manifest_matches,
        "analysis_parameters manifest must match the requested validation manifest",
        observed=recorded_manifest_resolved,
        expected=str(manifest_path),
    )

    recorded_joint_id = parameters.get("joint_id")
    report.add(
        "GLOBAL",
        "combined_parameters_joint_id",
        str(recorded_joint_id) == joint_id,
        "analysis_parameters joint_id must match --joint-id",
        observed=recorded_joint_id if recorded_joint_id is not None else "missing",
        expected=joint_id,
    )

    recorded_datasets = parameters.get("datasets")
    if isinstance(recorded_datasets, list):
        compare_id_sets(
            report,
            "GLOBAL",
            "combined_parameters_datasets",
            expected_ids,
            recorded_datasets,
        )
    else:
        report.add(
            "GLOBAL",
            "combined_parameters_datasets",
            False,
            "analysis_parameters datasets must be a list",
            observed=type(recorded_datasets).__name__,
            expected="list",
        )

    recorded_output = parameters.get("combined_output_dir")
    if recorded_output is None:
        output_matches = False
        recorded_output_resolved = "missing"
    else:
        recorded_output_resolved = str(
            root_relative(Path(str(recorded_output)), root)
        )
        output_matches = Path(recorded_output_resolved) == combined
    report.add(
        "GLOBAL",
        "combined_parameters_output_dir",
        output_matches,
        "analysis_parameters combined_output_dir must match configured aggregate directory",
        observed=recorded_output_resolved,
        expected=str(combined),
    )


def read_bed_ids(path: Path) -> pd.DataFrame:
    try:
        result = pd.read_csv(
            path,
            sep="\t",
            header=None,
            usecols=[0, 1, 2, 3, 6],
            low_memory=False,
        )
        result.columns = ["chrom", "start", "end", "acre_id", "array_id"]
        result["start"] = pd.to_numeric(result.start, errors="coerce")
        result["end"] = pd.to_numeric(result.end, errors="coerce")
        result = result.dropna(subset=["start", "end"]).copy()
        result[["start", "end"]] = result[["start", "end"]].astype(int)
        return result
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=["chrom", "start", "end", "acre_id", "array_id"])


def interval_containment_failures(path: Path, units: pd.DataFrame, acre_column: int) -> tuple[int, int]:
    try:
        intervals = pd.read_csv(
            path,
            sep="\t",
            header=None,
            usecols=[0, 1, 2, acre_column],
            low_memory=False,
        )
    except pd.errors.EmptyDataError:
        return 0, 0
    intervals.columns = ["chrom", "start", "end", "acre_id"]
    intervals["start"] = pd.to_numeric(intervals.start, errors="coerce")
    intervals["end"] = pd.to_numeric(intervals.end, errors="coerce")
    intervals = intervals.dropna(subset=["start", "end"]).copy()
    intervals[["start", "end"]] = intervals[["start", "end"]].astype(int)
    bounds = units[["acre_id", "chrom", "start", "end"]].rename(
        columns={"chrom": "unit_chrom", "start": "unit_start", "end": "unit_end"}
    )
    joined = intervals.merge(bounds, on="acre_id", how="left", validate="many_to_one")
    failures = (
        joined.unit_start.isna()
        | joined.chrom.astype(str).ne(joined.unit_chrom.astype(str))
        | joined.start.lt(joined.unit_start)
        | joined.end.gt(joined.unit_end)
        | joined.end.le(joined.start)
    )
    return len(joined), int(failures.sum())


def complete_marker(root: Path, dataset: str) -> Path | None:
    candidates = [
        root / "runs" / dataset / ".complete.json",
        root / "runs" / dataset / "complete.json",
        root / "runs" / dataset / f"{dataset}.complete.json",
    ]
    return next((path for path in candidates if path.exists()), None)


def validate_dataset(
    report: Report, root: Path, dataset: str, row: Mapping[str, object]
) -> None:
    prefix = output_prefix(root, dataset, row)
    assignments_path = selected_assignments_path(root, dataset, row)

    missing_native = []
    for suffix in NATIVE_SUFFIXES:
        path = Path(f"{prefix}{suffix}")
        if not check_exists(report, dataset, f"native_file:{suffix}", path):
            missing_native.append(path)
    if not check_exists(report, dataset, "selected_assignments", assignments_path):
        return
    if missing_native:
        return

    marker = complete_marker(root, dataset)
    report.add(
        dataset,
        "run_complete_marker",
        marker is not None,
        str(marker) if marker else "No completion marker found",
        severity="warning",
        expected="runs/<dataset>/.complete.json",
    )
    command = root / "runs" / dataset / "command.json"
    if command.exists():
        text = command.read_text()
        has_range = (
            "--min-period" in text
            and "6000" in text
            and "--max-period" in text
            and "9000" in text
        )
        report.add(
            dataset,
            "command_uses_6_to_9kb_period_range",
            has_range,
            str(command),
            observed="6000/9000 present" if has_range else "range flags/values absent",
            expected="--min-period 6000 --max-period 9000",
        )
    else:
        report.add(
            dataset,
            "command_provenance",
            False,
            str(command),
            severity="warning",
            expected="command.json",
        )

    try:
        units = pd.read_csv(Path(f"{prefix}.clusters.tsv"), sep="\t")
        chosen = normalize_assignments(pd.read_csv(assignments_path, sep="\t"), dataset)
    except Exception as error:
        report.exception(dataset, "read_native_and_assignments", error)
        return

    report.add(
        dataset,
        "unit_ids_unique",
        not units.acre_id.duplicated().any(),
        "NTRIAGE clusters acre_id uniqueness",
        observed=int(units.acre_id.duplicated().sum()),
        expected=0,
    )
    report.add(
        dataset,
        "assignment_ids_unique",
        not chosen.acre_id.duplicated().any(),
        "Selected assignment acre_id uniqueness",
        observed=int(chosen.acre_id.duplicated().sum()),
        expected=0,
    )
    compare_id_sets(
        report, dataset, "assignments_cover_all_units", units.acre_id, chosen.acre_id
    )
    if not len(units):
        report.add(dataset, "nonempty_units", False, "No complete units")
        return

    length = pd.to_numeric(
        units.get("length", units.end.astype(int) - units.start.astype(int)), errors="coerce"
    )
    computed_length = units.end.astype(int) - units.start.astype(int)
    report.add(
        dataset,
        "coordinate_length_consistency",
        bool(length.eq(computed_length).all()),
        "end - start equals length",
        observed=int((length != computed_length).sum()),
        expected=0,
    )
    outside = int((length.lt(6000) | length.gt(9000) | length.isna()).sum())
    report.add(
        dataset,
        "complete_units_within_6_to_9kb",
        outside == 0,
        "Inclusive requested complete-unit range",
        observed=f"outside={outside}; min={length.min()}; max={length.max()}",
        expected="6000 <= length <= 9000",
    )
    selected_k_values = chosen.selected_k.unique()
    k = int(selected_k_values[0]) if len(selected_k_values) == 1 else -1
    clusters = sorted(chosen.selected_cluster.unique())
    report.add(
        dataset,
        "one_selected_k",
        len(selected_k_values) == 1,
        f"values={selected_k_values.tolist()}",
        observed=len(selected_k_values),
        expected=1,
    )
    report.add(
        dataset,
        "selected_k_matches_cluster_count",
        k == len(clusters),
        f"labels={clusters}",
        observed=len(clusters),
        expected=k,
    )
    counts = chosen.selected_cluster.value_counts()
    report.add(
        dataset,
        "selected_clusters_nonempty",
        len(counts) == k and counts.min() > 0,
        f"min_cluster_units={int(counts.min())}; max_cluster_units={int(counts.max())}",
        observed=len(counts),
        expected=k,
    )
    label_counts = chosen.groupby("selected_cluster").selected_cluster_label.nunique()
    report.add(
        dataset,
        "one_stable_label_per_cluster",
        label_counts.eq(1).all(),
        f"clusters_with_multiple_labels={int(label_counts.gt(1).sum())}",
        expected=0,
    )

    arrays = pd.read_csv(Path(f"{prefix}.arrays.tsv"), sep="\t")
    report.add(
        dataset,
        "unit_array_ids_present",
        set(units.array_id).issubset(set(arrays.array_id)),
        f"missing={len(set(units.array_id)-set(arrays.array_id))}",
        expected=0,
    )
    if "n_arrays" in row and str(row["n_arrays"]).strip():
        expected_arrays = int(float(str(row["n_arrays"])))
        report.add(
            dataset,
            "manifest_array_count",
            len(arrays) == expected_arrays,
            "NTRIAGE arrays.tsv rows versus manifest",
            observed=len(arrays),
            expected=expected_arrays,
        )

    features = pd.read_csv(
        Path(f"{prefix}.features.tsv"), sep="\t", usecols=["acre_id"]
    )
    compare_id_sets(report, dataset, "feature_ids_match_units", units.acre_id, features.acre_id)
    acres = read_bed_ids(Path(f"{prefix}.acres.bed"))
    compare_id_sets(report, dataset, "acre_bed_ids_match_units", units.acre_id, acres.acre_id)
    if len(acres) == len(units):
        coordinates = units[["acre_id", "chrom", "start", "end"]].merge(
            acres[["acre_id", "chrom", "start", "end"]],
            on="acre_id",
            suffixes=("_table", "_bed"),
            validate="one_to_one",
        )
        mismatch = (
            coordinates.chrom_table.astype(str).ne(coordinates.chrom_bed.astype(str))
            | coordinates.start_table.ne(coordinates.start_bed)
            | coordinates.end_table.ne(coordinates.end_bed)
        )
        report.add(
            dataset,
            "acre_bed_coordinates_match_units",
            not mismatch.any(),
            "clusters.tsv versus acres.bed",
            observed=int(mismatch.sum()),
            expected=0,
        )

    composition = pd.read_csv(Path(f"{prefix}.acro1_composition.tsv"), sep="\t")
    compare_id_sets(
        report, dataset, "composition_ids_match_units", units.acre_id, composition.acre_id
    )
    if all(column in composition for column in ("subfamily_1", "subfamily_2", "subfamily_3", "classified_total")):
        computed = composition[["subfamily_1", "subfamily_2", "subfamily_3"]].sum(axis=1)
        mismatch = computed.ne(composition.classified_total)
        report.add(
            dataset,
            "composition_total_consistency",
            not mismatch.any(),
            "S1 + S2 + S3 equals classified_total",
            observed=int(mismatch.sum()),
            expected=0,
        )

    for name, suffix, acre_column in (
        ("annotation", ".annotations.bed", 9),
        ("acro1_subfamily", ".acro1_subfamilies.bed", 9),
    ):
        try:
            intervals, failures = interval_containment_failures(
                Path(f"{prefix}{suffix}"), units, acre_column
            )
            report.add(
                dataset,
                f"{name}_intervals_within_units",
                failures == 0,
                f"checked {intervals} intervals",
                observed=failures,
                expected=0,
            )
        except Exception as error:
            report.exception(dataset, f"{name}_interval_containment", error)

    analysis_dir = root / "results" / dataset / "cluster_analysis"
    missing_analysis = []
    for filename in ANALYSIS_FILES:
        path = analysis_dir / filename
        if not check_exists(report, dataset, f"analysis_file:{filename}", path):
            missing_analysis.append(path)
    if missing_analysis:
        return

    try:
        metrics = pd.read_csv(analysis_dir / "unit_metrics.tsv.gz", sep="\t")
        summary = pd.read_csv(analysis_dir / "cluster_summary.tsv", sep="\t")
        bins = pd.read_csv(analysis_dir / "cluster_unit_length_bins.tsv", sep="\t")
        compare_id_sets(
            report, dataset, "derived_unit_ids_match_native", units.acre_id, metrics.acre_id
        )
        report.add(
            dataset,
            "cluster_summary_has_one_row_per_cluster",
            set(summary.cluster) == set(clusters) and not summary.cluster.duplicated().any(),
            f"summary_clusters={sorted(summary.cluster.tolist())}",
            observed=len(summary),
            expected=k,
        )
        report.add(
            dataset,
            "cluster_summary_unit_total",
            int(summary.units.sum()) == len(units),
            "Sum of per-cluster unit counts",
            observed=int(summary.units.sum()),
            expected=len(units),
        )
        report.add(
            dataset,
            "length_bin_unit_total",
            int(bins.units.sum()) == len(units),
            "Sum of length-bin counts",
            observed=int(bins.units.sum()),
            expected=len(units),
        )
        metric_bounds = {
            "raw_acro1_anchor_union_bp": metrics.length,
            "exact_acro1_union_bp": metrics.length,
            "classified_acro1_union_bp": metrics.length,
            "broad_acro1_annotation_bp": metrics.length,
            "explicit_backbone_bp": metrics.length,
            "unclassified_bp": metrics.length,
            "annotation_covered_bp": metrics.length,
            "non_anchor_remainder_bp": metrics.length,
            "non_acro1_remainder_bp": metrics.length,
        }
        for column, upper in metric_bounds.items():
            bad = metrics[column].lt(0) | metrics[column].gt(upper)
            report.add(
                dataset,
                f"derived_metric_bounds:{column}",
                not bad.any(),
                "Metric must be between zero and unit length",
                observed=int(bad.sum()),
                expected=0,
            )
        report.add(
            dataset,
            "legacy_exact_union_alias_matches_raw_anchor_union",
            metrics.exact_acro1_union_bp.equals(metrics.raw_acro1_anchor_union_bp),
            "Legacy exact_acro1_union_bp alias versus explicit raw-anchor union",
            observed=int(
                metrics.exact_acro1_union_bp.ne(
                    metrics.raw_acro1_anchor_union_bp
                ).sum()
            ),
            expected=0,
        )
        report.add(
            dataset,
            "legacy_non_acro1_alias_matches_non_anchor_remainder",
            metrics.non_acro1_remainder_bp.equals(metrics.non_anchor_remainder_bp),
            "Legacy non_acro1 remainder alias versus explicit non-anchor remainder",
            observed=int(
                metrics.non_acro1_remainder_bp.ne(
                    metrics.non_anchor_remainder_bp
                ).sum()
            ),
            expected=0,
        )
        accounting = metrics.annotation_accounting_difference_bp.abs()
        copy_anchor_difference = metrics.classified_copy_anchor_count_difference
        report.add(
            dataset,
            "classified_copy_count_matches_subfamily_anchor_count",
            copy_anchor_difference.eq(0).all(),
            "S1+S2+S3 composition counts versus subfamily BED anchor count",
            observed=int(copy_anchor_difference.ne(0).sum()),
            expected=0,
        )
        raw_required = {
            "raw_acro1_anchor_count",
            "raw_acro1_union_bp",
            "raw_acro1_anchors_per_kb",
            "classified_acro1_copy_number",
            "classified_acro1_anchor_count",
        }
        missing_raw = sorted(raw_required - set(metrics))
        report.add(
            dataset,
            "raw_and_classified_acro1_metrics_present",
            not missing_raw,
            f"missing={missing_raw}",
            observed=len(raw_required) - len(missing_raw),
            expected=len(raw_required),
        )
        if not missing_raw:
            report.add(
                dataset,
                "primary_copy_metric_is_raw_anchor_count",
                metrics.acro1_copy_number.equals(metrics.raw_acro1_anchor_count),
                "acro1_copy_number is the raw mapped-anchor proxy; classified count is separate",
                observed=int(metrics.acro1_copy_number.ne(metrics.raw_acro1_anchor_count).sum()),
                expected=0,
            )
        report.add(
            dataset,
            "annotation_category_accounting",
            accounting.eq(0).all(),
            "ACRO1 + explicit backbone + unclassified versus union annotation bp; "
            "non-zero can reflect category overlaps",
            severity="warning",
            observed=f"nonzero={int(accounting.ne(0).sum())}; max={int(accounting.max())}",
            expected=0,
        )
    except Exception as error:
        report.exception(dataset, "derived_analysis_consistency", error)


def main() -> None:
    args = parse_args()
    root = args.root.resolve()
    manifest_path = args.manifest.resolve()
    report = Report()
    try:
        manifest, id_col = load_manifest(manifest_path, args.datasets)
    except Exception as error:
        report.exception("GLOBAL", "load_manifest", error)
        manifest = pd.DataFrame()
        id_col = "dataset_id"

    if len(manifest):
        present = set(manifest[id_col])
        if args.allow_partial or args.datasets:
            report.add(
                "GLOBAL",
                "canonical_dataset_set",
                True,
                "Partial validation explicitly requested",
                severity="warning",
                observed=",".join(sorted(present)),
            )
        else:
            missing = sorted(REQUIRED_DATASETS - present)
            report.add(
                "GLOBAL",
                "canonical_dataset_set",
                not missing,
                f"missing={missing}",
                observed=len(present),
                expected=len(REQUIRED_DATASETS),
            )
        for row in manifest.to_dict("records"):
            dataset = str(row[id_col])
            try:
                validate_dataset(report, root, dataset, row)
            except Exception as error:
                report.exception(dataset, "unhandled_validation_exception", error)

    combined = root_relative(args.combined_output_dir, root)
    if len(manifest) and (not args.datasets or args.joint_id in set(manifest[id_col])):
        validate_combined_outputs(
            report,
            combined,
            manifest,
            id_col,
            manifest_path,
            args.joint_id,
            root,
        )

    output = root_relative(args.validation_output_dir, root)
    output.mkdir(parents=True, exist_ok=True)
    frame = report.frame()
    frame.to_csv(output / "validation_report.tsv", sep="\t", index=False)
    counts = frame.status.value_counts().to_dict() if len(frame) else {}
    failures = int((frame.status == "FAIL").sum()) if len(frame) else 1
    warnings = int((frame.status == "WARN").sum()) if len(frame) else 0
    summary: dict[str, Any] = {
        "manifest": str(manifest_path),
        "datasets_validated": manifest[id_col].tolist() if len(manifest) else [],
        "checks": len(frame),
        "status_counts": counts,
        "failures": failures,
        "warnings": warnings,
        "valid": failures == 0,
    }
    with (output / "validation_summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps(summary, indent=2), file=sys.stderr)
    if failures and not args.no_fail:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
