#!/usr/bin/env python3
"""Localize one strongest two-block ANI transition in each positive HPRC ACRO array.

The legacy find_recombination.py script stores only an array-level call based on
the fractions of high- and low-ANI ModDotPlot cells.  This script regenerates
4-kb ModDotPlot matrices and extends that logic by selecting the split whose
within-block versus cross-block ANI contrast has the largest Welch statistic.
The 30 retained legacy breakpoint fractions are used unchanged in the final
coordinate table and serve as an independent localization validation set.
"""

from __future__ import annotations

import argparse
import importlib.metadata
import json
import math
from pathlib import Path
import subprocess
import sys
import time

import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
OUT_ROOT = SCRIPT_DIR.parent
PROJECT_ROOT = OUT_ROOT.parent.parent
MANIFEST = PROJECT_ROOT / "NTRIAGE_OCT26" / "manifest.validated.tsv"
FULL_FASTA = (
    PROJECT_ROOT / "NTRIAGE_OCT26" / "manuscript_output"
    / "telomeric_ntriage_2794" / "inputs" / "telomeric_2794.fasta"
)
LEGACY_POINTS = PROJECT_ROOT / "sample_100" / "recombination_subset.tsv"
INPUT_DIR = OUT_ROOT / "inputs"
RAWDIR = OUT_ROOT / "moddotplot_w4000"
TABLE_DIR = OUT_ROOT / "tables"
LOG_DIR = OUT_ROOT / "logs"
POSITIVE_IDS = INPUT_DIR / "positive_array_ids.txt"
POSITIVE_FASTA = INPUT_DIR / "positive_arrays.fasta"
POINTS_OUT = TABLE_DIR / "positive_breakpoints_1275.tsv"
VALIDATION_OUT = TABLE_DIR / "localization_validation_30.tsv"
SUMMARY_OUT = TABLE_DIR / "localization_summary.tsv"
RUN_LOG = LOG_DIR / "localization_run.json"

WINDOW_BP = 4_000
HIGH_ANI = 98.0
LOW_ANI_MIN = 90.0
LOW_ANI_MAX = 96.0
MIN_HIGH_FRAC = 0.2
MIN_LOW_FRAC = 0.1


def run_checked(command: list[str]) -> float:
    started = time.perf_counter()
    subprocess.run(command, check=True)
    return time.perf_counter() - started


def prepare_positive_fasta(manifest: pd.DataFrame) -> float:
    INPUT_DIR.mkdir(parents=True, exist_ok=True)
    positives = manifest.loc[manifest.has_recombination == 1, "array_id"].tolist()
    POSITIVE_IDS.write_text("\n".join(positives) + "\n")
    extract_seconds = run_checked([
        "samtools", "faidx", str(FULL_FASTA), "-r", str(POSITIVE_IDS),
        "-o", str(POSITIVE_FASTA),
    ])
    run_checked(["samtools", "faidx", str(POSITIVE_FASTA)])
    with POSITIVE_FASTA.open() as handle:
        observed = sum(1 for line in handle if line.startswith(">"))
    if observed != len(positives):
        raise RuntimeError(f"Expected {len(positives)} FASTA records, observed {observed}")
    return extract_seconds


def run_moddotplot() -> float:
    RAWDIR.mkdir(parents=True, exist_ok=True)
    moddotplot = Path(sys.executable).with_name("moddotplot")
    if not moddotplot.exists():
        raise FileNotFoundError(f"ModDotPlot executable not found beside {sys.executable}")
    return run_checked([
        str(moddotplot), "static", "-f", str(POSITIVE_FASTA),
        "-w", str(WINDOW_BP), "--no-plot", "--no-hist",
        "--processes", "4", "-o", str(RAWDIR), "--quiet",
    ])


def read_matrix(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    data = pd.read_csv(path, sep="\t")
    starts = sorted(set(data.iloc[:, 1]).union(set(data.iloc[:, 4])))
    index = {value: i for i, value in enumerate(starts)}
    matrix = np.full((len(starts), len(starts)), np.nan)
    for left, right, value in zip(data.iloc[:, 1], data.iloc[:, 4], data.iloc[:, 6]):
        i, j = index[int(left)], index[int(right)]
        matrix[i, j] = matrix[j, i] = float(value)
    off_diagonal = data.loc[
        ~((data.iloc[:, 1] == data.iloc[:, 4]) & (data.iloc[:, 2] == data.iloc[:, 5])),
        data.columns[6],
    ].to_numpy(float)
    return np.asarray(starts, dtype=int), matrix, off_diagonal


def score_split(matrix: np.ndarray, split: int) -> dict[str, float] | None:
    left = matrix[:split, :split]
    right = matrix[split:, split:]
    cross = matrix[:split, split:]
    left_values = left[np.triu_indices(split, 1)]
    right_n = len(matrix) - split
    right_values = right[np.triu_indices(right_n, 1)]
    cross_values = cross.ravel()
    left_values = left_values[np.isfinite(left_values)]
    right_values = right_values[np.isfinite(right_values)]
    cross_values = cross_values[np.isfinite(cross_values)]
    if min(len(left_values), len(right_values), len(cross_values)) < 2:
        return None
    within_values = np.concatenate([left_values, right_values])
    within_mean = float(within_values.mean())
    cross_mean = float(cross_values.mean())
    within_variance = float(within_values.var(ddof=1))
    cross_variance = float(cross_values.var(ddof=1))
    standard_error = math.sqrt(
        within_variance / len(within_values) + cross_variance / len(cross_values)
    )
    welch_score = (within_mean - cross_mean) / standard_error if standard_error > 0 else float("nan")
    return {
        "welch_score": welch_score,
        "within_mean_ani": within_mean,
        "left_within_mean_ani": float(left_values.mean()),
        "right_within_mean_ani": float(right_values.mean()),
        "cross_mean_ani": cross_mean,
        "within_minus_cross_ani": within_mean - cross_mean,
        "equal_block_minus_cross_ani": (
            (float(left_values.mean()) + float(right_values.mean())) / 2 - cross_mean
        ),
        "within_cells": len(within_values),
        "cross_cells": len(cross_values),
    }


def localize_one(row: pd.Series, legacy: pd.DataFrame) -> dict[str, object]:
    array_id = str(row.array_id)
    bedpe = RAWDIR / array_id / f"{array_id}.bedpe"
    if not bedpe.exists():
        raise FileNotFoundError(bedpe)
    starts, matrix, off_diagonal = read_matrix(bedpe)
    core_start = int(row.original_start) - int(row.flanked_start) + 1
    core_length = int(row.original_end) - int(row.original_start)
    core_end = core_start + core_length
    candidates: list[dict[str, object]] = []
    for split in range(2, len(starts) - 2):
        position = int(starts[split])
        if not (core_start <= position <= core_end):
            continue
        values = score_split(matrix, split)
        if values is None or not np.isfinite(values["welch_score"]):
            continue
        candidates.append({"split": split, "position": position, **values})
    if not candidates:
        raise RuntimeError(f"No localizable split for {array_id}")
    inferred = max(candidates, key=lambda value: (float(value["welch_score"]), -int(value["split"])))

    legacy_row = legacy.loc[legacy.name == array_id]
    if len(legacy_row) == 1:
        legacy_fraction = float(legacy_row.iloc[0].breakpoint_frac)
        final_position = core_start + legacy_fraction * core_length
        source = "retained_legacy_coordinate"
        legacy_contrast = float(legacy_row.iloc[0].contrast)
    else:
        legacy_fraction = float("nan")
        final_position = float(inferred["position"])
        source = "inferred_4kb_moddotplot_two_block"
        legacy_contrast = float("nan")

    n_pairs = len(off_diagonal)
    high_fraction = float(np.mean(off_diagonal >= HIGH_ANI)) if n_pairs else float("nan")
    low_fraction = float(np.mean((off_diagonal >= LOW_ANI_MIN) & (off_diagonal <= LOW_ANI_MAX))) if n_pairs else float("nan")
    reproduced_flag = bool(
        n_pairs > 0 and high_fraction >= MIN_HIGH_FRAC and low_fraction >= MIN_LOW_FRAC
    )
    return {
        "array_id": array_id,
        "individual_id": row.individual_id,
        "haplotype": row.haplotype,
        "coordinate_source": source,
        "core_start_local": core_start,
        "core_end_local": core_end,
        "core_length_bp": core_length,
        "breakpoint_local_bp": final_position,
        "breakpoint_fraction_core": (final_position - core_start) / core_length,
        "inferred_breakpoint_local_bp": int(inferred["position"]),
        "inferred_breakpoint_fraction_core": (int(inferred["position"]) - core_start) / core_length,
        "legacy_breakpoint_fraction_core": legacy_fraction,
        "legacy_contrast": legacy_contrast,
        "inferred_split_after_window": int(inferred["split"]),
        "moddotplot_windows": len(starts),
        "moddotplot_pairs": n_pairs,
        "high_ani_fraction_ge_98": high_fraction,
        "low_ani_fraction_90_to_96": low_fraction,
        "legacy_binary_rule_reproduced": reproduced_flag,
        **{key: inferred[key] for key in (
            "welch_score", "within_mean_ani", "left_within_mean_ani",
            "right_within_mean_ani", "cross_mean_ani",
            "within_minus_cross_ani", "equal_block_minus_cross_ani",
            "within_cells", "cross_cells",
        )},
    }


def localize_all(manifest: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    legacy = pd.read_csv(LEGACY_POINTS, sep="\t")
    positives = manifest.loc[manifest.has_recombination == 1].copy()
    rows = []
    for index, (_, row) in enumerate(positives.iterrows(), start=1):
        rows.append(localize_one(row, legacy))
        if index % 100 == 0:
            print(f"Localized {index}/{len(positives)} arrays", flush=True)
    points = pd.DataFrame(rows)
    validation = points.loc[points.coordinate_source == "retained_legacy_coordinate"].copy()
    validation["absolute_error_bp"] = (
        validation.inferred_breakpoint_local_bp - validation.breakpoint_local_bp
    ).abs()
    validation["signed_error_bp"] = (
        validation.inferred_breakpoint_local_bp - validation.breakpoint_local_bp
    )
    validation["within_4kb"] = validation.absolute_error_bp <= WINDOW_BP
    validation["within_8kb"] = validation.absolute_error_bp <= 2 * WINDOW_BP
    summary = pd.DataFrame([
        ("full_manifest_arrays", len(manifest)),
        ("recombination_positive_arrays", len(positives)),
        ("retained_legacy_coordinates", int((points.coordinate_source == "retained_legacy_coordinate").sum())),
        ("newly_inferred_coordinates", int((points.coordinate_source == "inferred_4kb_moddotplot_two_block").sum())),
        ("legacy_binary_rule_reproduced_among_positives", int(points.legacy_binary_rule_reproduced.sum())),
        ("validation_median_absolute_error_bp", float(validation.absolute_error_bp.median())),
        ("validation_mean_absolute_error_bp", float(validation.absolute_error_bp.mean())),
        ("validation_within_4kb", int(validation.within_4kb.sum())),
        ("validation_within_8kb", int(validation.within_8kb.sum())),
        ("validation_maximum_absolute_error_bp", float(validation.absolute_error_bp.max())),
    ], columns=["metric", "value"])
    return points, validation, summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--skip-moddotplot", action="store_true")
    args = parser.parse_args()
    for directory in (INPUT_DIR, RAWDIR, TABLE_DIR, LOG_DIR):
        directory.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(MANIFEST, sep="\t", dtype={"haplotype": str})
    if len(manifest) != 2_794 or int(manifest.has_recombination.sum()) != 1_275:
        raise RuntimeError("Unexpected manifest counts")

    started = time.perf_counter()
    extract_seconds = prepare_positive_fasta(manifest)
    if args.skip_moddotplot:
        moddotplot_seconds = 0.0
    else:
        moddotplot_seconds = run_moddotplot()
    bedpe_count = len(list(RAWDIR.glob("*/*.bedpe")))
    if bedpe_count != 1_275:
        raise RuntimeError(f"Expected 1,275 ModDotPlot BEDPE files, observed {bedpe_count}")
    localization_started = time.perf_counter()
    points, validation, summary = localize_all(manifest)
    localization_seconds = time.perf_counter() - localization_started
    points.to_csv(POINTS_OUT, sep="\t", index=False, float_format="%.10g")
    validation.to_csv(VALIDATION_OUT, sep="\t", index=False, float_format="%.10g")
    summary.to_csv(SUMMARY_OUT, sep="\t", index=False, float_format="%.10g")
    record = {
        "window_bp": WINDOW_BP,
        "high_ani_threshold": HIGH_ANI,
        "low_ani_min": LOW_ANI_MIN,
        "low_ani_max": LOW_ANI_MAX,
        "minimum_high_fraction": MIN_HIGH_FRAC,
        "minimum_low_fraction": MIN_LOW_FRAC,
        "positive_arrays": 1_275,
        "moddotplot_version": importlib.metadata.version("ModDotPlot"),
        "extract_fasta_seconds": extract_seconds,
        "moddotplot_seconds": moddotplot_seconds,
        "localization_seconds": localization_seconds,
        "total_seconds": time.perf_counter() - started,
        "command": " ".join(sys.argv),
    }
    RUN_LOG.write_text(json.dumps(record, indent=2) + "\n")
    print(summary.to_string(index=False))
    print(f"Wrote {POINTS_OUT}")


if __name__ == "__main__":
    main()
