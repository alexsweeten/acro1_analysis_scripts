#!/usr/bin/env python3
"""Select a defensible cluster count from an NTRIAGE feature matrix.

The primary selection view caps the contribution of every parent array at four
complete ACRO units.  This prevents a handful of very long arrays from choosing
``K``.  Natural-abundance and (for the joint run) dataset-balanced views are
reported as sensitivity analyses.  The chosen K is then refit to *all* units
with 50 independent k-means++ starts.  If that full-data fit fails the cluster
support guard, progressively simpler sample-eligible K values are tested and
every rejected candidate is retained as an auditable diagnostic.

NTRIAGE's feature preprocessing is reproduced exactly: each numeric feature is
centered and divided by its sample standard deviation (ddof=1), constant
columns contribute zero, and each row is L2-normalized.  Scaling parameters are
estimated once from the complete run, before any balancing/subsampling.

Typical use
-----------
Run every row in the project manifest::

    python scripts/select_k.py --manifest inputs/dataset_manifest.tsv

Run one matrix directly::

    python scripts/select_k.py --dataset-id gorilla \
      --features runs/gorilla/gorilla.features.tsv \
      --arrays runs/gorilla/gorilla.arrays.tsv

The direct-mode output defaults to ``results/<dataset>/k_selection``.  In
manifest mode, the script accepts ``features_path``/``arrays_path`` directly or
derives them from ``ntriage_prefix`` (also ``output_prefix`` or ``run_prefix``).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import sys
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Sequence

# Avoid a noisy joblib physical-core probe on some macOS/Python builds.  This
# does not constrain BLAS/OpenMP threads; it only supplies the missing count.
os.environ.setdefault("LOKY_MAX_CPU_COUNT", str(os.cpu_count() or 1))

import numpy as np
import pandas as pd
import scipy
import sklearn
from scipy.optimize import linear_sum_assignment
from sklearn.cluster import KMeans
from sklearn.metrics import (
    adjusted_rand_score,
    calinski_harabasz_score,
    davies_bouldin_score,
    silhouette_score,
)


PROJECT = Path(__file__).resolve().parents[1]
METADATA_COLUMNS = ("acre_id", "array_id", "cluster")
DEFAULT_SEEDS = tuple(range(1, 11))
DEFAULT_FINAL_SEEDS = tuple(range(1001, 1051))
SELECTION_VERSION = 3
VISUALIZATION_VERSION = 2


@dataclass(frozen=True)
class DatasetSpec:
    dataset_id: str
    features: Path
    arrays: Path | None
    array_dataset_map: Path | None
    output_dir: Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Array-balanced multi-metric K selection and a 50-start final "
            "refit for NTRIAGE features."
        )
    )
    parser.add_argument("--manifest", type=Path, help="Dataset manifest TSV.")
    parser.add_argument(
        "--dataset",
        action="append",
        default=[],
        help="Dataset ID to run from the manifest (repeatable; default: all).",
    )
    parser.add_argument("--dataset-id", help="Dataset ID for direct mode.")
    parser.add_argument("--features", type=Path, help="NTRIAGE .features.tsv.")
    parser.add_argument("--arrays", type=Path, help="NTRIAGE .arrays.tsv.")
    parser.add_argument(
        "--array-dataset-map",
        type=Path,
        help=(
            "Optional TSV assigning array_id or chrom/source_name to dataset_id; "
            "used for the joint dataset-balanced sensitivity view."
        ),
    )
    parser.add_argument("--output-dir", type=Path, help="Direct-mode output directory.")
    parser.add_argument("--max-k", type=int, default=16)
    parser.add_argument("--sweep-starts", type=int, default=10)
    parser.add_argument(
        "--sensitivity-starts",
        type=int,
        default=5,
        help="Starts per K for non-primary sensitivity views (default: 5).",
    )
    parser.add_argument("--final-starts", type=int, default=50)
    parser.add_argument("--max-units-per-array", type=int, default=4)
    parser.add_argument(
        "--natural-sensitivity-cap",
        type=int,
        default=20000,
        help=(
            "Uniform natural-abundance sample cap for K sensitivity (0 means all). "
            "The final model is always fit to all units."
        ),
    )
    parser.add_argument("--silhouette-sample", type=int, default=3000)
    parser.add_argument("--min-cluster-units", type=int, default=5)
    parser.add_argument("--min-cluster-fraction", type=float, default=0.005)
    parser.add_argument("--min-cluster-arrays", type=int, default=3)
    parser.add_argument(
        "--min-cluster-mean-agreement",
        type=float,
        default=0.75,
        help=(
            "Minimum label-aligned mean start agreement required in every "
            "cluster, both during K selection and after the final refit."
        ),
    )
    parser.add_argument("--random-seed", type=int, default=1)
    parser.add_argument("--force", action="store_true", help="Replace complete outputs.")
    parser.add_argument("--no-plot", action="store_true")
    parser.add_argument(
        "--redraw-plots-only",
        action="store_true",
        help="Regenerate plots from saved k_selection.tsv files without fitting KMeans.",
    )
    args = parser.parse_args()

    if args.manifest is None and args.features is None:
        parser.error("provide --manifest or --features")
    if args.manifest is not None and args.features is not None:
        parser.error("--manifest and --features are mutually exclusive")
    if args.features is not None and not args.dataset_id:
        parser.error("direct mode requires --dataset-id")
    if args.max_k < 2:
        parser.error("--max-k must be at least 2")
    if args.sweep_starts < 2:
        parser.error("--sweep-starts must be at least 2 to estimate stability")
    if args.sensitivity_starts < 2:
        parser.error("--sensitivity-starts must be at least 2 to estimate stability")
    if args.final_starts < 2:
        parser.error("--final-starts must be at least 2 to estimate agreement")
    if args.max_units_per_array < 1:
        parser.error("--max-units-per-array must be positive")
    if not 0 <= args.min_cluster_fraction < 1:
        parser.error("--min-cluster-fraction must be in [0, 1)")
    if not 0 < args.min_cluster_mean_agreement <= 1:
        parser.error("--min-cluster-mean-agreement must be in (0, 1]")
    return args


def nonempty(value: object) -> bool:
    return value is not None and not pd.isna(value) and str(value).strip() != ""


def resolve_path(value: object, manifest: Path | None = None) -> Path | None:
    if not nonempty(value):
        return None
    path = Path(str(value)).expanduser()
    if path.is_absolute():
        return path
    candidates: list[Path] = []
    if manifest is not None:
        candidates.append(manifest.resolve().parent / path)
    candidates.append(PROJECT / path)
    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()
    # Output prefixes often do not exist when a manifest is first written.  A
    # project-relative path is the least surprising interpretation.
    return (PROJECT / path).resolve()


def first_value(row: pd.Series, names: Sequence[str]) -> object | None:
    for name in names:
        if name in row.index and nonempty(row[name]):
            return row[name]
    return None


def specs_from_manifest(path: Path, selected: set[str]) -> list[DatasetSpec]:
    path = path.resolve()
    frame = pd.read_csv(path, sep="\t", dtype=str)
    id_column = next(
        (name for name in ("dataset_id", "dataset", "name", "run_id") if name in frame),
        None,
    )
    if id_column is None:
        raise ValueError(f"{path}: no dataset_id/dataset/name/run_id column")

    specs: list[DatasetSpec] = []
    for _, row in frame.iterrows():
        dataset_id = str(row[id_column]).strip()
        if selected and dataset_id not in selected:
            continue
        feature_value = first_value(row, ("features_path", "features_file", "features"))
        prefix_value = first_value(
            row, ("ntriage_prefix", "output_prefix", "run_prefix", "prefix")
        )
        if feature_value is not None:
            features = resolve_path(feature_value, path)
        elif prefix_value is not None:
            prefix = resolve_path(prefix_value, path)
            assert prefix is not None
            features = Path(str(prefix) + ".features.tsv")
            if not features.exists() and Path(str(features) + ".gz").exists():
                features = Path(str(features) + ".gz")
        else:
            raise ValueError(
                f"manifest row {dataset_id!r} lacks features_path or ntriage_prefix"
            )
        assert features is not None

        arrays_value = first_value(row, ("arrays_path", "arrays_file", "arrays"))
        if arrays_value is not None:
            arrays = resolve_path(arrays_value, path)
        elif prefix_value is not None:
            prefix = resolve_path(prefix_value, path)
            arrays = Path(str(prefix) + ".arrays.tsv") if prefix is not None else None
            if arrays is not None and not arrays.exists():
                cluster_table = Path(str(prefix) + ".clusters.tsv")
                if cluster_table.exists():
                    arrays = cluster_table
        else:
            arrays = None

        map_value = first_value(
            row, ("array_dataset_map", "array_dataset_path", "cohort_map")
        )
        map_path = resolve_path(map_value, path)
        out_value = first_value(row, ("k_selection_dir", "results_dir"))
        output = (
            resolve_path(out_value, path)
            if out_value is not None
            else PROJECT / "results" / dataset_id / "k_selection"
        )
        assert output is not None
        specs.append(
            DatasetSpec(dataset_id, features, arrays, map_path, output.resolve())
        )

    missing = selected - {spec.dataset_id for spec in specs}
    if missing:
        raise ValueError(f"datasets absent from manifest: {', '.join(sorted(missing))}")
    return specs


def get_specs(args: argparse.Namespace) -> list[DatasetSpec]:
    if args.manifest is not None:
        return specs_from_manifest(args.manifest, set(args.dataset))
    output = args.output_dir or PROJECT / "results" / args.dataset_id / "k_selection"
    return [
        DatasetSpec(
            str(args.dataset_id),
            args.features.resolve(),
            args.arrays.resolve() if args.arrays else None,
            args.array_dataset_map.resolve() if args.array_dataset_map else None,
            output.resolve(),
        )
    ]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_json(value: object) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def build_selection_config(
    spec: DatasetSpec,
    args: argparse.Namespace,
    input_hashes: dict[str, str | None],
) -> dict[str, object]:
    """Return the scientific configuration used to validate resumability."""
    return {
        "dataset_id": spec.dataset_id,
        "input_features": str(spec.features.resolve()),
        "input_arrays": str(spec.arrays.resolve()) if spec.arrays else None,
        "array_dataset_map": (
            str(spec.array_dataset_map.resolve()) if spec.array_dataset_map else None
        ),
        "input_sha256": input_hashes,
        "max_k": args.max_k,
        "sweep_starts": args.sweep_starts,
        "sensitivity_starts": args.sensitivity_starts,
        "final_starts": args.final_starts,
        "max_units_per_array": args.max_units_per_array,
        "natural_sensitivity_cap": args.natural_sensitivity_cap,
        "silhouette_sample": args.silhouette_sample,
        "min_cluster_units": args.min_cluster_units,
        "min_cluster_fraction": args.min_cluster_fraction,
        "min_cluster_arrays": args.min_cluster_arrays,
        "min_cluster_mean_agreement": args.min_cluster_mean_agreement,
        "random_seed": args.random_seed,
        "preprocessing": "sample-SD standardization then row L2 normalization",
        "algorithm": "KMeans k-means++ Lloyd max_iter=300 tol=1e-4",
        "selection_version": SELECTION_VERSION,
    }


def scale_features(raw: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return NTRIAGE-equivalent standardized, row-L2-normalized features."""
    means = raw.mean(axis=0, dtype=np.float64)
    scales = raw.std(axis=0, ddof=1, dtype=np.float64)
    usable = np.isfinite(scales) & (scales > 1e-12)
    if not usable.any():
        raise ValueError("all feature columns are constant")
    x = (raw[:, usable] - means[usable]) / scales[usable]
    norms = np.linalg.norm(x, axis=1)
    nonzero = norms > 0
    x[nonzero] /= norms[nonzero, None]
    x[~nonzero] = 0.0
    if not np.isfinite(x).all():
        raise ValueError("non-finite values after feature preprocessing")
    return x, means, scales, usable


def stable_labels(labels: np.ndarray) -> tuple[np.ndarray, dict[int, int]]:
    """Relabel 0-based clusters 1..K by first input occurrence."""
    first = [(int(np.flatnonzero(labels == label)[0]), int(label)) for label in np.unique(labels)]
    remap = {old: new + 1 for new, (_, old) in enumerate(sorted(first))}
    stable = np.fromiter((remap[int(x)] for x in labels), dtype=np.int32, count=len(labels))
    return stable, remap


def deterministic_cap_by_group(
    groups: Sequence[object], cap: int, seed: int
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    positions = pd.Series(np.arange(len(groups), dtype=np.int64)).groupby(
        pd.Series(groups, dtype="object"), sort=False
    )
    chosen: list[np.ndarray] = []
    for _, series in positions:
        values = series.to_numpy(dtype=np.int64)
        if len(values) > cap:
            values = np.sort(rng.choice(values, size=cap, replace=False))
        chosen.append(values)
    return np.sort(np.concatenate(chosen)) if chosen else np.empty(0, dtype=np.int64)


def deterministic_natural_sample(n: int, cap: int, seed: int) -> np.ndarray:
    if cap <= 0 or n <= cap:
        return np.arange(n, dtype=np.int64)
    rng = np.random.default_rng(seed)
    return np.sort(rng.choice(n, size=cap, replace=False)).astype(np.int64)


def normalize_cohort(value: object) -> str:
    token = str(value).strip().split("|")[0].lower()
    token = token.split(":")[0]
    mappings = (
        (("hg002", "human"), "hg002"),
        (("mgorgor1", "gorilla"), "gorilla"),
        (("mpantro3", "chimp"), "chimpanzee"),
        (("mpanpan1", "bonobo"), "bonobo"),
        (("mponabe1", "sumatran"), "sumatran_orangutan"),
        (("mponpyg2", "bornean"), "bornean_orangutan"),
    )
    for needles, label in mappings:
        if any(needle in token for needle in needles):
            return label
    return "hprc"


def infer_array_cohorts(
    features: pd.DataFrame, arrays_path: Path | None, map_path: Path | None
) -> np.ndarray | None:
    """Infer one cohort per feature row for joint-run sensitivity analysis."""
    if arrays_path is None or not arrays_path.exists():
        return None
    arrays = pd.read_csv(arrays_path, sep="\t", dtype=str)
    if "array_id" not in arrays:
        return None
    # A .clusters.tsv is also accepted as the array lookup.  It has one row per
    # unit rather than per array, so collapse it before constructing mappings.
    arrays = arrays.drop_duplicates("array_id", keep="first").copy()

    cohort_column = next(
        (
            name
            for name in ("source_dataset", "dataset_id", "cohort", "species_id")
            if name in arrays
        ),
        None,
    )
    if cohort_column is not None:
        array_to_cohort = arrays.set_index("array_id")[cohort_column].to_dict()
    else:
        source_column = next(
            (name for name in ("source_name", "name", "chrom") if name in arrays), None
        )
        if source_column is None:
            return None
        arrays["_cohort"] = arrays[source_column].map(normalize_cohort)
        array_to_cohort = arrays.set_index("array_id")["_cohort"].to_dict()

    if map_path is not None and map_path.exists():
        mapping = pd.read_csv(map_path, sep="\t", dtype=str)
        cohort_map_column = next(
            (
                name
                for name in ("source_dataset", "dataset_id", "cohort", "species_id")
                if name in mapping
            ),
            None,
        )
        if cohort_map_column is not None:
            # Prepared metadata calls the input FASTA/BED sequence key
            # ``analysis_chrom`` and also has an ``array_id`` column containing
            # that external key.  NTRIAGE's internal IDs are ARRAY_000001, ...;
            # therefore chrom-based joins must take precedence over a direct
            # array_id join.
            if "analysis_chrom" in mapping and "chrom" in arrays:
                chrom_to_cohort = mapping.set_index("analysis_chrom")[
                    cohort_map_column
                ].to_dict()
                override = arrays["chrom"].map(chrom_to_cohort)
                array_to_cohort.update(
                    dict(zip(arrays.loc[override.notna(), "array_id"], override.dropna()))
                )
            elif "chrom" in mapping and "chrom" in arrays:
                chrom_to_cohort = mapping.set_index("chrom")[cohort_map_column].to_dict()
                override = arrays["chrom"].map(chrom_to_cohort)
                array_to_cohort.update(
                    dict(zip(arrays.loc[override.notna(), "array_id"], override.dropna()))
                )
            elif "source_name" in mapping and "source_name" in arrays:
                source_to_cohort = mapping.set_index("source_name")[cohort_map_column].to_dict()
                override = arrays["source_name"].map(source_to_cohort)
                array_to_cohort.update(
                    dict(zip(arrays.loc[override.notna(), "array_id"], override.dropna()))
                )
            elif "array_id" in mapping:
                internal_ids = set(arrays["array_id"])
                direct = mapping[mapping["array_id"].isin(internal_ids)]
                array_to_cohort.update(
                    direct.set_index("array_id")[cohort_map_column].to_dict()
                )

    cohort = features["array_id"].map(array_to_cohort)
    if cohort.isna().any():
        return None
    return cohort.astype(str).to_numpy()


def dataset_balanced_indices(
    base_indices: np.ndarray, cohorts: np.ndarray, seed: int
) -> np.ndarray | None:
    labels = pd.Series(cohorts[base_indices], dtype="object")
    groups = labels.groupby(labels, sort=True).groups
    if len(groups) < 2:
        return None
    minimum = min(len(indexes) for indexes in groups.values())
    if minimum < 2:
        return None
    rng = np.random.default_rng(seed)
    chosen: list[np.ndarray] = []
    for label in sorted(groups):
        relative = np.asarray(groups[label], dtype=np.int64)
        absolute = base_indices[relative]
        if len(absolute) > minimum:
            absolute = rng.choice(absolute, size=minimum, replace=False)
        chosen.append(np.sort(absolute))
    return np.sort(np.concatenate(chosen))


def minmax(series: pd.Series, higher_is_better: bool = True) -> pd.Series:
    numeric = pd.to_numeric(series, errors="coerce")
    valid = numeric[np.isfinite(numeric)]
    result = pd.Series(np.zeros(len(series), dtype=float), index=series.index)
    if valid.empty:
        return result
    low, high = float(valid.min()), float(valid.max())
    if math.isclose(low, high):
        result.loc[valid.index] = 1.0
    else:
        result.loc[valid.index] = (valid - low) / (high - low)
    if higher_is_better:
        return result
    inverted = pd.Series(np.zeros(len(series), dtype=float), index=series.index)
    inverted.loc[valid.index] = 1.0 - result.loc[valid.index]
    return inverted


def fit_once(x: np.ndarray, k: int, seed: int) -> KMeans:
    return KMeans(
        n_clusters=k,
        init="k-means++",
        n_init=1,
        max_iter=300,
        tol=1e-4,
        random_state=seed,
        algorithm="lloyd",
    ).fit(x)


def evaluate_view(
    x: np.ndarray,
    frame: pd.DataFrame,
    indices: np.ndarray,
    view: str,
    is_primary: bool,
    args: argparse.Namespace,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    sample = x[indices]
    sample_arrays = frame.iloc[indices]["array_id"].astype(str).to_numpy()
    max_k = min(args.max_k, len(sample) - 1)
    if max_k < 2:
        raise ValueError(f"{view}: fewer than three sampled units")
    starts = args.sweep_starts if is_primary else args.sensitivity_starts
    seeds = list(range(1, starts + 1))
    total_ss = float(np.square(sample - sample.mean(axis=0)).sum())
    summaries: list[dict[str, object]] = []
    runs: list[dict[str, object]] = []

    min_units_required = max(
        args.min_cluster_units,
        int(math.ceil(args.min_cluster_fraction * len(sample))),
    )
    min_arrays_required = min(args.min_cluster_arrays, len(np.unique(sample_arrays)))

    for k in range(2, max_k + 1):
        print(f"  {view}: K={k}, {len(seeds)} starts", flush=True)
        fitted: list[tuple[int, KMeans, np.ndarray]] = []
        for seed in seeds:
            model = fit_once(sample, k, seed)
            labels, _ = stable_labels(model.labels_)
            counts = np.bincount(labels, minlength=k + 1)[1:]
            actual_k = len(np.unique(labels))
            runs.append(
                {
                    "view": view,
                    "k": k,
                    "seed": seed,
                    "inertia": float(model.inertia_),
                    "iterations": int(model.n_iter_),
                    "actual_clusters": actual_k,
                    "min_cluster_size": int(counts.min()) if len(counts) else 0,
                    "max_cluster_size": int(counts.max()) if len(counts) else 0,
                }
            )
            fitted.append((seed, model, labels))

        best_seed, best_model, best_labels = min(fitted, key=lambda value: value[1].inertia_)
        counts = np.bincount(best_labels, minlength=k + 1)[1:]
        actual_k = len(np.unique(best_labels))
        array_counts = np.array(
            [len(np.unique(sample_arrays[best_labels == cluster])) for cluster in range(1, k + 1)],
            dtype=int,
        )
        pairwise_ari = [
            adjusted_rand_score(left[2], right[2])
            for left, right in combinations(fitted, 2)
        ]
        aligned_starts = np.column_stack(
            [
                align_to_reference(best_labels, model.labels_, k)
                for _, model, _ in fitted
            ]
        )
        unit_agreement = (aligned_starts == best_labels[:, None]).mean(axis=1)
        cluster_mean_agreement = np.array(
            [
                float(unit_agreement[best_labels == cluster].mean())
                if np.any(best_labels == cluster)
                else 0.0
                for cluster in range(1, k + 1)
            ],
            dtype=float,
        )
        for run_row, aligned_labels in zip(runs[-len(fitted):], aligned_starts.T):
            run_row["ari_to_best_partition"] = adjusted_rand_score(
                best_labels, aligned_labels
            )
            run_row["unit_agreement_to_best_partition"] = float(
                (aligned_labels == best_labels).mean()
            )
        valid = actual_k == k and len(counts) == k and counts.min() > 0
        if valid:
            silhouette = silhouette_score(
                sample,
                best_labels,
                metric="euclidean",
                sample_size=min(args.silhouette_sample, len(sample)),
                random_state=args.random_seed,
            )
            ch = calinski_harabasz_score(sample, best_labels)
            db = davies_bouldin_score(sample, best_labels)
        else:
            silhouette = ch = db = np.nan
        summaries.append(
            {
                "view": view,
                "sample_units": len(sample),
                "sample_arrays": len(np.unique(sample_arrays)),
                "k": k,
                "best_seed": best_seed,
                "inertia": float(best_model.inertia_),
                "explained_variance_fraction": (
                    1.0 - float(best_model.inertia_) / total_ss if total_ss > 0 else np.nan
                ),
                "silhouette": silhouette,
                "calinski_harabasz": ch,
                "davies_bouldin": db,
                "stability_ari_mean": float(np.mean(pairwise_ari)),
                "stability_ari_sd": (
                    float(np.std(pairwise_ari, ddof=1)) if len(pairwise_ari) > 1 else 0.0
                ),
                "stability_ari_min": float(np.min(pairwise_ari)),
                "mean_unit_start_agreement": float(unit_agreement.mean()),
                "weakest_cluster_mean_agreement": float(
                    cluster_mean_agreement.min()
                ),
                "actual_clusters": actual_k,
                "min_cluster_size": int(counts.min()) if len(counts) else 0,
                "max_cluster_size": int(counts.max()) if len(counts) else 0,
                "min_cluster_arrays": int(array_counts.min()) if len(array_counts) else 0,
                "cluster_size_cv": (
                    float(counts.std(ddof=1) / counts.mean()) if len(counts) > 1 else 0.0
                ),
                "required_min_cluster_units": min_units_required,
                "required_min_cluster_arrays": min_arrays_required,
                "required_min_cluster_mean_agreement": (
                    args.min_cluster_mean_agreement
                ),
            }
        )

    result = pd.DataFrame(summaries).sort_values("k").reset_index(drop=True)
    inertia = result["inertia"].to_numpy(dtype=float)
    relative_drop = np.full(len(result), np.nan)
    relative_drop[1:] = 1.0 - inertia[1:] / inertia[:-1]
    result["relative_inertia_drop"] = relative_drop
    elbow = np.zeros(len(result), dtype=float)
    for i in range(1, len(result) - 1):
        previous_gain = 1.0 - inertia[i] / inertia[i - 1]
        next_gain = 1.0 - inertia[i + 1] / inertia[i]
        elbow[i] = previous_gain - next_gain
    result["elbow_strength"] = elbow

    result["score_silhouette"] = minmax(result["silhouette"])
    result["score_calinski_harabasz"] = minmax(result["calinski_harabasz"])
    result["score_davies_bouldin"] = minmax(result["davies_bouldin"], False)
    result["score_pairwise_ari"] = minmax(result["stability_ari_mean"])
    result["score_weakest_cluster_agreement"] = minmax(
        result["weakest_cluster_mean_agreement"]
    )
    result["score_stability"] = 0.5 * (
        result["score_pairwise_ari"]
        + result["score_weakest_cluster_agreement"]
    )
    result["score_elbow"] = minmax(result["elbow_strength"])
    result["score_explained_variance"] = minmax(result["explained_variance_fraction"])
    result["robust_score"] = (
        0.25 * result["score_silhouette"]
        + 0.15 * result["score_calinski_harabasz"]
        + 0.15 * result["score_davies_bouldin"]
        + 0.25 * result["score_stability"]
        + 0.10 * result["score_elbow"]
        + 0.10 * result["score_explained_variance"]
    )
    result["support_eligible"] = (
        (result["actual_clusters"] == result["k"])
        & (result["min_cluster_size"] >= result["required_min_cluster_units"])
        & (result["min_cluster_arrays"] >= result["required_min_cluster_arrays"])
        & (
            result["weakest_cluster_mean_agreement"]
            >= result["required_min_cluster_mean_agreement"]
        )
    )
    eligible = result[result["support_eligible"]]
    selection_relaxed = False
    if eligible.empty:
        # Unit/array support may be relaxed for an unusually small input, but
        # cluster stability is a scientific validity guard and is never relaxed.
        eligible = result[
            (result["actual_clusters"] == result["k"])
            & (result["min_cluster_size"] >= 2)
            & (
                result["weakest_cluster_mean_agreement"]
                >= result["required_min_cluster_mean_agreement"]
            )
        ]
        selection_relaxed = True
    if eligible.empty:
        best_stability = float(result["weakest_cluster_mean_agreement"].max())
        raise RuntimeError(
            f"{view}: no K has weakest-cluster mean agreement >= "
            f"{args.min_cluster_mean_agreement:.3f}; best observed={best_stability:.3f}"
        )
    result["selection_eligible"] = result.index.isin(eligible.index)
    maximum = float(eligible["robust_score"].max())
    candidates = eligible[eligible["robust_score"] >= 0.95 * maximum]
    chosen_k = int(candidates["k"].min())
    result["view_selected"] = result["k"] == chosen_k
    result["within_95pct_best_score"] = (
        result["selection_eligible"] & (result["robust_score"] >= 0.95 * maximum)
    )
    result["selection_support_relaxed"] = selection_relaxed
    return result, pd.DataFrame(runs)


def align_to_reference(reference: np.ndarray, alternate: np.ndarray, k: int) -> np.ndarray:
    contingency = np.zeros((k, k), dtype=np.int64)
    np.add.at(contingency, (reference - 1, alternate), 1)
    rows, cols = linear_sum_assignment(-contingency)
    mapping = {int(column): int(row + 1) for row, column in zip(rows, cols)}
    # Degenerate starts can leave a label unmatched; map it to its largest-overlap
    # reference cluster so agreement remains defined and the problem is visible in
    # the per-start actual_clusters field.
    for label in np.unique(alternate):
        if int(label) not in mapping:
            mapping[int(label)] = int(np.argmax(contingency[:, int(label)]) + 1)
    return np.fromiter(
        (mapping[int(label)] for label in alternate), dtype=np.int32, count=len(alternate)
    )


def final_refit(
    x: np.ndarray,
    frame: pd.DataFrame,
    dataset_id: str,
    k: int,
    starts: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    fitted: list[tuple[int, KMeans]] = []
    for offset in range(starts):
        seed = 1001 + offset
        print(f"  final K={k}: start {offset + 1}/{starts}", flush=True)
        fitted.append((seed, fit_once(x, k, seed)))
    best_seed, best_model = min(fitted, key=lambda value: value[1].inertia_)
    reference, remap = stable_labels(best_model.labels_)

    aligned = np.empty((len(x), starts), dtype=np.int16)
    refit_rows: list[dict[str, object]] = []
    for column, (seed, model) in enumerate(fitted):
        aligned[:, column] = align_to_reference(reference, model.labels_, k)
        refit_rows.append(
            {
                "seed": seed,
                "k": k,
                "inertia": float(model.inertia_),
                "iterations": int(model.n_iter_),
                "actual_clusters": len(np.unique(model.labels_)),
                "selected_best_inertia": seed == best_seed,
                "ari_to_selected": adjusted_rand_score(reference, aligned[:, column]),
            }
        )

    agreement = (aligned == reference[:, None]).mean(axis=1)
    votes = np.zeros((len(x), k), dtype=np.int16)
    row_indexes = np.arange(len(x))
    for column in range(starts):
        votes[row_indexes, aligned[:, column] - 1] += 1
    ordered_votes = np.sort(votes, axis=1)
    vote_margin = (ordered_votes[:, -1] - ordered_votes[:, -2]) / starts

    # KMeans.transform avoids constructing an n x K x feature broadcast array.
    distances_raw = best_model.transform(x)
    distances = np.empty_like(distances_raw)
    for raw_label, stable_label in remap.items():
        distances[:, stable_label - 1] = distances_raw[:, raw_label]
    assigned_distance = distances[row_indexes, reference - 1]
    second_distance = np.partition(distances, kth=1, axis=1)[:, 1]
    distance_margin = second_distance - assigned_distance

    safe_id = re.sub(r"[^A-Za-z0-9]+", "_", dataset_id).strip("_").upper()
    assignments = pd.DataFrame(
        {
            "acre_id": frame["acre_id"].astype(str),
            "array_id": frame["array_id"].astype(str),
            "selected_k": k,
            "selected_cluster": reference,
            "selected_cluster_label": [
                f"{safe_id}_K{k:02d}_C{cluster:02d}" for cluster in reference
            ],
            "agreement": agreement,
            "vote_margin": vote_margin,
            "distance_to_centroid": assigned_distance,
            "second_distance": second_distance,
            "distance_margin": distance_margin,
            "n_final_starts": starts,
        }
    )

    cluster_rows: list[dict[str, object]] = []
    for cluster in range(1, k + 1):
        mask = reference == cluster
        cluster_rows.append(
            {
                "selected_cluster": cluster,
                "selected_cluster_label": f"{safe_id}_K{k:02d}_C{cluster:02d}",
                "units": int(mask.sum()),
                "arrays": int(frame.loc[mask, "array_id"].nunique()),
                "fraction_units": float(mask.mean()),
                "mean_agreement": float(agreement[mask].mean()),
                "agreement_q05": float(np.quantile(agreement[mask], 0.05)),
                "fraction_agreement_below_0.75": float((agreement[mask] < 0.75).mean()),
                "median_vote_margin": float(np.median(vote_margin[mask])),
                "median_distance_margin": float(np.median(distance_margin[mask])),
            }
        )
    return assignments, pd.DataFrame(refit_rows), pd.DataFrame(cluster_rows)


def assess_final_support(
    support: pd.DataFrame,
    total_units: int,
    args: argparse.Namespace,
) -> tuple[pd.DataFrame, bool, list[str]]:
    """Apply the same support guard to the authoritative full-data refit."""
    checked = support.copy()
    minimum_units = max(
        args.min_cluster_units,
        int(math.ceil(args.min_cluster_fraction * total_units)),
    )
    checked["required_min_units"] = minimum_units
    checked["required_min_arrays"] = args.min_cluster_arrays
    checked["required_min_mean_agreement"] = args.min_cluster_mean_agreement
    checked["passes_unit_support"] = checked["units"] >= minimum_units
    checked["passes_array_support"] = checked["arrays"] >= args.min_cluster_arrays
    checked["passes_mean_agreement"] = (
        checked["mean_agreement"] >= args.min_cluster_mean_agreement
    )
    checked["support_eligible"] = (
        checked["passes_unit_support"]
        & checked["passes_array_support"]
        & checked["passes_mean_agreement"]
    )
    reasons: list[str] = []
    for row in checked.loc[~checked["support_eligible"]].itertuples(index=False):
        failures: list[str] = []
        if not row.passes_unit_support:
            failures.append(f"units={row.units}<{minimum_units}")
        if not row.passes_array_support:
            failures.append(f"arrays={row.arrays}<{args.min_cluster_arrays}")
        if not row.passes_mean_agreement:
            failures.append(
                f"mean_agreement={row.mean_agreement:.3f}"
                f"<{args.min_cluster_mean_agreement:.3f}"
            )
        reasons.append(f"{row.selected_cluster_label}: {', '.join(failures)}")
    return checked, bool(checked["support_eligible"].all()), reasons


def write_plot(selection: pd.DataFrame, output: Path, dataset_id: str) -> None:
    required = {"preferred_selected", "authoritative_selected", "view_selected"}
    missing = sorted(required - set(selection.columns))
    if missing:
        raise ValueError(
            f"cannot draw authoritative K plot before decision columns exist: {missing}"
        )
    authoritative = selection[selection["authoritative_selected"].astype(bool)]
    preferred = selection[selection["preferred_selected"].astype(bool)]
    if len(authoritative) != 1 or len(preferred) != 1:
        raise ValueError(
            "k_selection.tsv must contain exactly one preferred and one "
            "authoritative selection"
        )
    authoritative_k = int(authoritative.iloc[0]["k"])
    preferred_k = int(preferred.iloc[0]["k"])
    authoritative_view = str(authoritative.iloc[0]["view"])
    preferred_view = str(preferred.iloc[0]["view"])

    decision_path = output.parent / "selection_decision.json"
    decision: dict[str, object] = {}
    if decision_path.exists():
        decision = json.loads(decision_path.read_text())
    rejected_pairs = {
        (str(candidate["view"]), int(candidate["k"]))
        for candidate in decision.get("rejected_candidates", [])
        if "view" in candidate and "k" in candidate
    }
    # Older decisions may omit the rejected list; a preferred result that is not
    # authoritative necessarily failed the full-data guard.
    if (preferred_view, preferred_k) != (authoritative_view, authoritative_k):
        rejected_pairs.add((preferred_view, preferred_k))

    os.environ.setdefault("MPLCONFIGDIR", str(output.parent / ".mplconfig"))
    os.environ.setdefault("XDG_CACHE_HOME", str(output.parent / ".cache"))
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    fig, axes = plt.subplots(2, 3, figsize=(14, 8), constrained_layout=True)
    colors = {
        "array_balanced": "#1261A0",
        "natural_abundance": "#D1495B",
        "dataset_balanced": "#2A9D8F",
    }
    panels = (
        ("robust_score", "Robust composite (higher is better)"),
        ("silhouette", "Silhouette (higher is better)"),
        ("stability_ari_mean", "Multi-start stability ARI"),
        ("weakest_cluster_mean_agreement", "Weakest-cluster mean agreement"),
        ("inertia", "Within-cluster sum of squares"),
        ("elbow_strength", "Elbow strength"),
    )
    for axis, (column, title) in zip(axes.flat, panels):
        for view, subset in selection.groupby("view", sort=False):
            axis.plot(
                subset["k"],
                subset[column],
                marker=".",
                markersize=3.5,
                linewidth=1.5,
                color=colors.get(view),
            )
            sample_choice = subset[subset["view_selected"].astype(bool)]
            if not sample_choice.empty:
                axis.scatter(
                    sample_choice["k"], sample_choice[column], s=48,
                    marker="x", color=colors.get(view), linewidths=1.7, zorder=4,
                )

        if np.isfinite(float(preferred.iloc[0][column])):
            axis.scatter(
                [preferred_k], [preferred.iloc[0][column]], s=90, marker="D",
                facecolors="none", edgecolors="#F4A261", linewidths=2.1, zorder=6,
            )
        for rejected_view, rejected_k in sorted(rejected_pairs):
            rejected = selection[
                (selection["view"] == rejected_view) & (selection["k"] == rejected_k)
            ]
            if not rejected.empty and np.isfinite(float(rejected.iloc[0][column])):
                axis.scatter(
                    [rejected_k], [rejected.iloc[0][column]], s=105, marker="X",
                    color="#C1121F", edgecolors="white", linewidths=0.7, zorder=7,
                )
        axis.axvline(
            authoritative_k, color="#111111", linestyle="--", linewidth=1.3,
            alpha=0.8, zorder=2,
        )
        if np.isfinite(float(authoritative.iloc[0][column])):
            axis.scatter(
                [authoritative_k], [authoritative.iloc[0][column]], s=145, marker="*",
                color="#FFD166", edgecolors="#111111", linewidths=1.0, zorder=8,
            )
        axis.set_title(title)
        axis.set_xlabel("K")
        axis.grid(alpha=0.25)
        if column == "weakest_cluster_mean_agreement":
            threshold = float(selection["required_min_cluster_mean_agreement"].iloc[0])
            axis.axhline(threshold, color="black", linestyle="--", linewidth=1)
    view_handles = [
        Line2D([0], [0], color=color, linewidth=2, label=view.replace("_", " "))
        for view, color in colors.items()
        if view in set(selection["view"])
    ]
    semantic_handles = [
        Line2D(
            [0], [0], marker="x", linestyle="none", color="#555555",
            markersize=7, label="sample-selected K (per view)",
        ),
        Line2D(
            [0], [0], marker="D", linestyle="none", markerfacecolor="none",
            markeredgecolor="#F4A261", markeredgewidth=2, markersize=7,
            label="preferred sample K",
        ),
    ]
    if rejected_pairs:
        semantic_handles.append(
            Line2D(
                [0], [0], marker="X", linestyle="none", color="#C1121F",
                markeredgecolor="white", markersize=8,
                label="rejected by full-data stability",
            )
        )
    semantic_handles.append(
        Line2D(
            [0], [0], marker="*", linestyle="--", color="#111111",
            markerfacecolor="#FFD166", markersize=11,
            label=f"authoritative K={authoritative_k}",
        )
    )
    fig.legend(
        handles=view_handles + semantic_handles,
        loc="outside lower center",
        ncol=min(4, len(view_handles + semantic_handles)),
        frameon=False,
        fontsize=8.5,
    )
    if (preferred_view, preferred_k) == (authoritative_view, authoritative_k):
        subtitle = f"authoritative K={authoritative_k} (preferred full-data-supported)"
    else:
        subtitle = (
            f"authoritative K={authoritative_k}; preferred {preferred_view} "
            f"K={preferred_k} rejected"
        )
    fig.suptitle(f"{dataset_id}: cluster-number sensitivity\n{subtitle}", fontsize=14)
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def redraw_saved_plots(specs: Sequence[DatasetSpec]) -> None:
    """Redraw K figures and migrate plot-only provenance without refitting."""
    current_script_sha = sha256_file(Path(__file__).resolve())
    for spec in specs:
        selection_path = spec.output_dir / "k_selection.tsv"
        parameters_path = spec.output_dir / "k_selection_parameters.json"
        if not selection_path.exists():
            raise FileNotFoundError(selection_path)
        selection = pd.read_csv(selection_path, sep="\t")
        output = spec.output_dir / "k_selection.png"
        write_plot(selection, output, spec.dataset_id)

        # Older version-3 configurations included the entire script SHA in the
        # scientific fingerprint.  A plot-only edit should not invalidate fitted
        # assignments, so migrate that SHA to explicit provenance and retain the
        # algorithm-version/input/parameter fingerprint as the resume key.
        if parameters_path.exists():
            parameters = json.loads(parameters_path.read_text())
            config = dict(parameters.get("selection_config", {}))
            fit_script_sha = config.pop("selector_script_sha256", None)
            if fit_script_sha is not None:
                parameters.setdefault(
                    "selection_script_sha256_at_fit", fit_script_sha
                )
            parameters["selection_config"] = config
            parameters["selection_config_sha256"] = sha256_json(config)
            parameters["visualization_version"] = VISUALIZATION_VERSION
            parameters["artifact_script_sha256_last_updated"] = current_script_sha
            parameters["plot_semantics"] = (
                "x=per-view sample selection; orange diamond=preferred sample K; "
                "red X=full-data rejection; gold star/dashed line=authoritative K"
            )
            parameters_path.write_text(json.dumps(parameters, indent=2) + "\n")
        print(f"{spec.dataset_id}: redrew {output}", flush=True)


def run_dataset(spec: DatasetSpec, args: argparse.Namespace) -> None:
    required_outputs = (
        spec.output_dir / "k_selection.tsv",
        spec.output_dir / "assignments.tsv.gz",
        spec.output_dir / "final_cluster_support.tsv",
        spec.output_dir / "selection_decision.json",
        spec.output_dir / "k_selection_parameters.json",
    )
    if not spec.features.exists():
        raise FileNotFoundError(spec.features)
    if spec.arrays is not None and not spec.arrays.exists():
        raise FileNotFoundError(spec.arrays)
    if spec.array_dataset_map is not None and not spec.array_dataset_map.exists():
        raise FileNotFoundError(spec.array_dataset_map)

    input_hashes: dict[str, str | None] = {
        "features": sha256_file(spec.features),
        "arrays": sha256_file(spec.arrays) if spec.arrays else None,
        "array_dataset_map": (
            sha256_file(spec.array_dataset_map) if spec.array_dataset_map else None
        ),
    }
    scientific_config = build_selection_config(spec, args, input_hashes)
    scientific_config_sha = sha256_json(scientific_config)
    existing = [path for path in required_outputs if path.exists()]
    if not args.force and existing:
        if len(existing) != len(required_outputs):
            raise RuntimeError(
                f"{spec.dataset_id}: partial K-selection outputs exist; use --force "
                "after inspecting them"
            )
        try:
            prior = json.loads(required_outputs[-1].read_text())
        except (OSError, json.JSONDecodeError) as error:
            raise RuntimeError(
                f"{spec.dataset_id}: unreadable prior parameters; use --force"
            ) from error
        if prior.get("selection_config_sha256") != scientific_config_sha:
            raise RuntimeError(
                f"{spec.dataset_id}: existing assignments were produced from a "
                "different feature/configuration fingerprint; use --force"
            )
        if not prior.get("final_refit_support", {}).get("pass", False):
            raise RuntimeError(
                f"{spec.dataset_id}: prior final refit did not pass support; use --force"
            )
        print(
            f"{spec.dataset_id}: verified feature/configuration fingerprint; skipping"
        )
        return
    spec.output_dir.mkdir(parents=True, exist_ok=True)
    if args.force:
        # These two files are completion markers.  Removing them up front keeps
        # a failed forced rerun from presenting stale assignments as current.
        for marker in (
            spec.output_dir / "k_selection_parameters.json",
            spec.output_dir / "selection_decision.json",
        ):
            if marker.exists():
                marker.unlink()
    print(f"{spec.dataset_id}: reading {spec.features}", flush=True)
    frame = pd.read_csv(spec.features, sep="\t")
    missing = [name for name in METADATA_COLUMNS if name not in frame]
    if missing:
        raise ValueError(f"{spec.features}: missing columns {missing}")
    if frame["acre_id"].duplicated().any():
        raise ValueError(f"{spec.features}: duplicate acre_id values")
    feature_names = list(frame.columns[3:])
    if not feature_names:
        raise ValueError(f"{spec.features}: no numeric feature columns after metadata")
    metadata_frame = frame.loc[:, list(METADATA_COLUMNS)].copy()
    raw = frame.iloc[:, 3:].to_numpy(dtype=np.float64, copy=True)
    # Release the duplicate pandas-owned feature block before allocating the
    # standardized matrix; this materially lowers peak memory on the joint run.
    del frame
    if not np.isfinite(raw).all():
        raise ValueError(f"{spec.features}: feature matrix contains non-finite values")
    x, means, scales, usable = scale_features(raw)
    del raw
    frame = metadata_frame
    pd.DataFrame(
        {
            "feature": feature_names,
            "mean": means,
            "sample_sd": scales,
            "used_nonconstant": usable,
        }
    ).to_csv(spec.output_dir / "feature_scaling.tsv", sep="\t", index=False)

    array_balanced = deterministic_cap_by_group(
        frame["array_id"].astype(str).to_numpy(),
        args.max_units_per_array,
        args.random_seed,
    )
    natural = deterministic_natural_sample(
        len(frame), args.natural_sensitivity_cap, args.random_seed
    )
    is_joint = spec.dataset_id in {
        "joint",
        "joint_all",
        "joint_primates",
        "all",
        "united",
    }
    views: list[tuple[str, np.ndarray]] = [
        ("array_balanced", array_balanced),
        ("natural_abundance", natural),
    ]
    cohorts = infer_array_cohorts(frame, spec.arrays, spec.array_dataset_map)
    dataset_balanced = None
    if cohorts is not None and is_joint:
        dataset_balanced = dataset_balanced_indices(
            array_balanced, cohorts, args.random_seed
        )
        if dataset_balanced is not None:
            views.append(("dataset_balanced", dataset_balanced))
    if is_joint and dataset_balanced is None:
        raise RuntimeError(
            "joint run requires a valid dataset-balanced view; check arrays_path "
            "and array_dataset_map cohort joins"
        )

    # The comparative joint question makes the equally weighted dataset view
    # the preferred selection.  Array-balanced remains a high-reliability
    # fallback because full-data stability can expose structures that exist only
    # in the balanced subsample.  Both receive the full sweep-start budget.
    preferred_view = "dataset_balanced" if is_joint else "array_balanced"
    fallback_view = "array_balanced" if is_joint else None
    full_sweep_views = {preferred_view}
    if fallback_view is not None:
        full_sweep_views.add(fallback_view)

    sample_rows: list[pd.DataFrame] = []
    selections: list[pd.DataFrame] = []
    all_runs: list[pd.DataFrame] = []
    for view, indices in views:
        if len(indices) < 3:
            continue
        sample_frame = frame.iloc[indices][["acre_id", "array_id"]].copy()
        sample_frame.insert(0, "view", view)
        if cohorts is not None:
            sample_frame["dataset_id"] = cohorts[indices]
        sample_rows.append(sample_frame)
        summary, runs = evaluate_view(
            x, frame, indices, view, view in full_sweep_views, args
        )
        selections.append(summary)
        all_runs.append(runs)

    selection = pd.concat(selections, ignore_index=True)
    view_choices: dict[str, int] = {}
    for view, subset in selection.groupby("view"):
        chosen = subset.loc[subset["view_selected"], "k"]
        if len(chosen) != 1:
            raise RuntimeError(f"{view} did not select exactly one K")
        view_choices[view] = int(chosen.iloc[0])
    preferred_k = view_choices[preferred_view]
    selection["preferred_selected"] = (
        (selection["view"] == preferred_view) & (selection["k"] == preferred_k)
    )
    pd.concat(all_runs, ignore_index=True).to_csv(
        spec.output_dir / "k_selection_runs.tsv.gz", sep="\t", index=False,
        compression="gzip",
    )
    pd.concat(sample_rows, ignore_index=True).to_csv(
        spec.output_dir / "view_samples.tsv.gz", sep="\t", index=False,
        compression="gzip",
    )

    print(
        f"{spec.dataset_id}: preferred {preferred_view} selection K={preferred_k}; "
        f"refitting {len(frame):,} units",
        flush=True,
    )
    def simpler_candidates(view: str, upper_k: int) -> list[tuple[str, int, str]]:
        """Order smaller eligible K by the documented defensible-band rule."""
        subset = selection[
            (selection["view"] == view)
            & selection["selection_eligible"]
            & (selection["k"] < upper_k)
        ].copy()
        if subset.empty:
            return []
        best_score = float(subset["robust_score"].max())
        near_best = subset[subset["robust_score"] >= 0.95 * best_score].sort_values(
            "k", ascending=True
        )
        remainder = subset[subset["robust_score"] < 0.95 * best_score].sort_values(
            ["robust_score", "k"], ascending=[False, True]
        )
        ordered = pd.concat([near_best, remainder], ignore_index=True)
        return [
            (view, int(row.k), "smaller_selection_eligible")
            for row in ordered.itertuples(index=False)
        ]

    candidate_queue: list[tuple[str, int, str]] = [
        (preferred_view, preferred_k, "preferred_view_selection")
    ]
    if fallback_view is not None:
        fallback_k = view_choices[fallback_view]
        if fallback_k < preferred_k:
            candidate_queue.append(
                (fallback_view, fallback_k, "fallback_view_selection")
            )
            candidate_queue.extend(simpler_candidates(fallback_view, fallback_k))
        else:
            candidate_queue.extend(simpler_candidates(preferred_view, preferred_k))
    else:
        candidate_queue.extend(simpler_candidates(preferred_view, preferred_k))

    rejected_candidates: list[dict[str, object]] = []
    candidate_attempts: list[dict[str, object]] = []
    attempted_k: set[int] = set()
    authoritative_view: str | None = None
    selected_k: int | None = None
    assignments: pd.DataFrame | None = None
    refit_runs: pd.DataFrame | None = None
    support: pd.DataFrame | None = None

    for candidate_view, candidate_k, candidate_basis in candidate_queue:
        # A full-data refit depends on K, not on the view that proposed it.
        # Never spend another 50 starts retrying a known failed numeric K.
        if candidate_k in attempted_k:
            continue
        attempted_k.add(candidate_k)
        if len(candidate_attempts) > 0:
            print(
                f"{spec.dataset_id}: trying simpler {candidate_view} K={candidate_k}",
                flush=True,
            )
        candidate_assignments, candidate_runs, support_raw = final_refit(
            x, frame, spec.dataset_id, candidate_k, args.final_starts
        )
        candidate_support, support_pass, support_reasons = assess_final_support(
            support_raw, len(frame), args
        )
        attempt = {
            "view": candidate_view,
            "k": candidate_k,
            "basis": candidate_basis,
            "pass": support_pass,
            "weakest_cluster_mean_agreement": float(
                candidate_support["mean_agreement"].min()
            ),
            "failure_details": support_reasons,
        }
        candidate_attempts.append(attempt)
        if support_pass:
            authoritative_view = candidate_view
            selected_k = candidate_k
            assignments = candidate_assignments
            refit_runs = candidate_runs
            support = candidate_support
            break

        rejected_dir = (
            spec.output_dir / "candidate_refits"
            / f"{candidate_view}_K{candidate_k:02d}_rejected"
        )
        rejected_dir.mkdir(parents=True, exist_ok=True)
        candidate_assignments.to_csv(
            rejected_dir / "assignments.tsv.gz", sep="\t", index=False,
            compression="gzip",
        )
        candidate_runs.to_csv(
            rejected_dir / "final_refit_runs.tsv", sep="\t", index=False
        )
        candidate_support.to_csv(
            rejected_dir / "final_cluster_support.tsv", sep="\t", index=False
        )
        rejection = {
            "view": candidate_view,
            "k": candidate_k,
            "basis": candidate_basis,
            "reason": "full-data final refit failed support",
            "details": support_reasons,
            "weakest_cluster_mean_agreement": float(
                candidate_support["mean_agreement"].min()
            ),
            "diagnostics_dir": str(rejected_dir),
        }
        rejected_candidates.append(rejection)
        (rejected_dir / "rejection.json").write_text(
            json.dumps(
                {**rejection, "selection_config_sha256": scientific_config_sha},
                indent=2,
            )
            + "\n"
        )
        print(
            f"{spec.dataset_id}: rejected {candidate_view} K={candidate_k}: "
            + "; ".join(support_reasons),
            flush=True,
        )

    if authoritative_view is None or selected_k is None:
        failure = {
            "dataset_id": spec.dataset_id,
            "preferred_view": preferred_view,
            "preferred_k": preferred_k,
            "candidate_attempts": candidate_attempts,
            "rejected_candidates": rejected_candidates,
            "selection_config_sha256": scientific_config_sha,
        }
        (spec.output_dir / "final_refit_support_failure.json").write_text(
            json.dumps(failure, indent=2) + "\n"
        )
        raise RuntimeError(
            f"{spec.dataset_id}: no simpler selection-eligible K passed the "
            "full-data support guard"
        )

    assert assignments is not None and refit_runs is not None and support is not None
    selection["authoritative_selected"] = (
        (selection["view"] == authoritative_view) & (selection["k"] == selected_k)
    )
    # Backward-compatible name used by downstream analysis; unlike
    # preferred_selected, this always denotes the full-data-supported result.
    selection["primary_selected"] = selection["authoritative_selected"]
    selection["preferred_view"] = preferred_view
    selection["authoritative_view"] = authoritative_view
    selection["authoritative_k"] = selected_k
    selection.to_csv(spec.output_dir / "k_selection.tsv", sep="\t", index=False)
    assignments.to_csv(
        spec.output_dir / "assignments.tsv.gz", sep="\t", index=False,
        compression="gzip",
    )
    refit_runs.to_csv(spec.output_dir / "final_refit_runs.tsv", sep="\t", index=False)
    support.to_csv(spec.output_dir / "final_cluster_support.tsv", sep="\t", index=False)

    decision = {
        "dataset_id": spec.dataset_id,
        "preferred_view": preferred_view,
        "preferred_k": preferred_k,
        "authoritative_view": authoritative_view,
        "authoritative_k": selected_k,
        "fallback_used": (
            authoritative_view != preferred_view or selected_k != preferred_k
        ),
        "candidate_attempts": candidate_attempts,
        "rejected_candidates": rejected_candidates,
        "final_refit_support_pass": True,
        "minimum_cluster_mean_agreement": args.min_cluster_mean_agreement,
        "weakest_authoritative_cluster_mean_agreement": float(
            support["mean_agreement"].min()
        ),
    }
    (spec.output_dir / "selection_decision.json").write_text(
        json.dumps(decision, indent=2) + "\n"
    )
    failure_marker = spec.output_dir / "final_refit_support_failure.json"
    if failure_marker.exists():
        failure_marker.unlink()

    if not args.no_plot:
        write_plot(selection, spec.output_dir / "k_selection.png", spec.dataset_id)

    sensitivity_choices = {
        view: int(subset.loc[subset["view_selected"], "k"].iloc[0])
        for view, subset in selection.groupby("view")
    }
    parameters = {
        "dataset_id": spec.dataset_id,
        "input_features": str(spec.features),
        "input_features_sha256": input_hashes["features"],
        "input_arrays": str(spec.arrays) if spec.arrays else None,
        "input_sha256": input_hashes,
        "selection_config": scientific_config,
        "selection_config_sha256": scientific_config_sha,
        "selection_script_sha256_at_fit": sha256_file(Path(__file__).resolve()),
        "visualization_version": VISUALIZATION_VERSION,
        "observations": len(frame),
        "parent_arrays": int(frame["array_id"].nunique()),
        "raw_features": len(feature_names),
        "nonconstant_features": int(usable.sum()),
        "preprocessing": (
            "full-matrix sample-SD column standardization (ddof=1), constant "
            "features set to zero, followed by row L2 normalization"
        ),
        "preferred_view": preferred_view,
        "preferred_k": preferred_k,
        "primary_view": authoritative_view,
        "authoritative_view": authoritative_view,
        "fallback_used": (
            authoritative_view != preferred_view or selected_k != preferred_k
        ),
        "candidate_attempts": candidate_attempts,
        "rejected_candidates": rejected_candidates,
        "max_units_per_array": args.max_units_per_array,
        "natural_sensitivity_cap": args.natural_sensitivity_cap,
        "view_sample_sizes": {view: len(indexes) for view, indexes in views},
        "sensitivity_selected_k": sensitivity_choices,
        "selected_k": selected_k,
        "k_values": list(range(2, min(args.max_k, len(array_balanced) - 1) + 1)),
        "sweep_starts": args.sweep_starts,
        "sweep_seeds": list(range(1, args.sweep_starts + 1)),
        "sensitivity_starts": args.sensitivity_starts,
        "sensitivity_seeds": list(range(1, args.sensitivity_starts + 1)),
        "final_starts": args.final_starts,
        "final_seeds": list(range(1001, 1001 + args.final_starts)),
        "algorithm": "scikit-learn KMeans; k-means++; Lloyd; n_init=1 per explicit start",
        "max_iter": 300,
        "tolerance": 1e-4,
        "silhouette_sample_size": args.silhouette_sample,
        "support_rule": {
            "minimum_units": args.min_cluster_units,
            "minimum_fraction_of_view": args.min_cluster_fraction,
            "minimum_parent_arrays": args.min_cluster_arrays,
            "minimum_cluster_mean_agreement": args.min_cluster_mean_agreement,
        },
        "score_weights": {
            "silhouette": 0.25,
            "calinski_harabasz": 0.15,
            "inverse_davies_bouldin": 0.15,
            "stability_composite_pairwise_ari_and_weakest_cluster_agreement": 0.25,
            "elbow_strength": 0.10,
            "explained_variance": 0.10,
        },
        "selection_rule": (
            "Among support-eligible K, choose the smallest K whose robust score "
            "is at least 95% of the maximum. Unit/array support may be relaxed "
            "for small inputs, but weakest-cluster mean agreement >= threshold "
            "is mandatory. If a full-data refit fails, try only smaller "
            "selection-eligible K values, ordering values within 95% of the "
            "remaining best score from smallest to largest and then by score. "
            "Joint data try the array-balanced selected K first."
        ),
        "final_refit_support": {
            "pass": True,
            "weakest_cluster_mean_agreement": float(
                support["mean_agreement"].min()
            ),
            "threshold": args.min_cluster_mean_agreement,
        },
        "cluster_label_rule": "1-based order of first unit occurrence within this run",
        "numpy_version": np.__version__,
        "pandas_version": pd.__version__,
        "scipy_version": scipy.__version__,
        "scikit_learn_version": sklearn.__version__,
        "command": " ".join(sys.argv),
    }
    (spec.output_dir / "k_selection_parameters.json").write_text(
        json.dumps(parameters, indent=2) + "\n"
    )
    print(
        f"{spec.dataset_id}: complete; K={selected_k}; "
        f"assignments={spec.output_dir / 'assignments.tsv.gz'}",
        flush=True,
    )


def main() -> None:
    args = parse_args()
    specs = get_specs(args)
    if not specs:
        raise ValueError("no datasets selected")
    if args.redraw_plots_only:
        redraw_saved_plots(specs)
        return
    failures: list[tuple[str, Exception]] = []
    for spec in specs:
        try:
            run_dataset(spec, args)
        except Exception as error:  # continue an all-dataset batch, then fail loudly
            failures.append((spec.dataset_id, error))
            print(f"ERROR {spec.dataset_id}: {error}", file=sys.stderr, flush=True)
            if len(specs) == 1:
                raise
    if failures:
        detail = "; ".join(f"{dataset}: {error}" for dataset, error in failures)
        raise RuntimeError(f"{len(failures)} dataset(s) failed: {detail}")


if __name__ == "__main__":
    main()
