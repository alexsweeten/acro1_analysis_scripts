#!/usr/bin/env python3
"""Test whether HG002 H1-H4 map onto the authoritative HPRC K=3 clusters.

Two complementary analyses are performed:

1. Fixed-centroid projection: HG002 is standardized with the saved HPRC
   scaling parameters and assigned to centroids reconstructed from the
   authoritative HPRC K=3 assignments.
2. Joint refit: the HPRC and HG002 feature matrices are standardized together,
   K=3 is refit, and the resulting labels are aligned to the original HPRC
   labels by maximum overlap among HPRC ACREs.

The joint fit is run both from the reconstructed HPRC centroids and with an
independent 10-start k-means++ fit. Agreement between those fits is retained as
an additional robustness check.
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent
OUT = HERE / "results"

HPRC_FEATURES = PROJECT / "runs/hprc/native/hprc.features.tsv"
HPRC_ASSIGNMENTS = PROJECT / "results/hprc/k_selection/assignments.tsv.gz"
HPRC_SCALING = PROJECT / "results/hprc/k_selection/feature_scaling.tsv"
HG002_FEATURES = (
    PROJECT / "hg002_acrocentric_k3/native/hg002_acrocentric_k3.features.tsv"
)
HG002_ASSIGNMENTS = (
    PROJECT
    / "hg002_acrocentric_k4_split_weakest/results/"
    "constrained_and_global_k4_unit_assignments.tsv.gz"
)

METADATA = ["acre_id", "array_id", "cluster"]
HPRC_LABELS = np.array([1, 2, 3], dtype=np.int32)


def read_features(path: Path) -> tuple[pd.DataFrame, np.ndarray, list[str]]:
    columns = pd.read_csv(path, sep="\t", nrows=0).columns.tolist()
    if columns[:3] != METADATA:
        raise ValueError(f"{path}: unexpected metadata columns {columns[:3]}")
    feature_names = columns[3:]
    dtypes: dict[str, object] = {name: np.float32 for name in feature_names}
    dtypes.update({name: str for name in METADATA})
    frame = pd.read_csv(path, sep="\t", dtype=dtypes)
    if frame[["acre_id", "array_id"]].duplicated().any():
        raise ValueError(f"{path}: duplicate ACRE/array keys")
    raw = frame.iloc[:, 3:].to_numpy(dtype=np.float32, copy=True)
    metadata = frame.loc[:, METADATA].copy()
    return metadata, raw, feature_names


def normalize_rows(x: np.ndarray) -> np.ndarray:
    norms = np.linalg.norm(x, axis=1)
    nonzero = norms > 0
    x[nonzero] /= norms[nonzero, None]
    x[~nonzero] = 0
    return x


def scale_saved(
    raw: np.ndarray, scaling: pd.DataFrame, feature_names: list[str]
) -> np.ndarray:
    if scaling.feature.tolist() != feature_names:
        raise ValueError("saved HPRC scaling features do not match feature matrix")
    usable = scaling.used_nonconstant.astype(str).str.lower().eq("true").to_numpy()
    means = scaling.loc[usable, "mean"].to_numpy(dtype=np.float32)
    scales = scaling.loc[usable, "sample_sd"].to_numpy(dtype=np.float32)
    x = (raw[:, usable] - means) / scales
    return normalize_rows(x)


def scale_joint(raw: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    means = raw.mean(axis=0, dtype=np.float64).astype(np.float32)
    scales = raw.std(axis=0, ddof=1, dtype=np.float64).astype(np.float32)
    usable = np.isfinite(scales) & (scales > 1e-12)
    x = (raw[:, usable] - means[usable]) / scales[usable]
    return normalize_rows(x), usable


def ordered_labels(
    metadata: pd.DataFrame, assignments: pd.DataFrame, label_column: str
) -> np.ndarray:
    merged = metadata[["acre_id", "array_id"]].merge(
        assignments[["acre_id", "array_id", label_column]],
        on=["acre_id", "array_id"],
        how="left",
        validate="one_to_one",
    )
    if merged[label_column].isna().any():
        raise ValueError(f"missing labels after merge: {label_column}")
    return merged[label_column].to_numpy()


def centroids(x: np.ndarray, labels: np.ndarray) -> np.ndarray:
    return np.vstack([x[labels == label].mean(axis=0) for label in HPRC_LABELS])


def distances_to_centroids(x: np.ndarray, centers: np.ndarray) -> np.ndarray:
    # Squared Euclidean distance; avoids a large three-dimensional broadcast.
    x2 = np.einsum("ij,ij->i", x, x)[:, None]
    c2 = np.einsum("ij,ij->i", centers, centers)[None, :]
    return np.maximum(x2 + c2 - 2 * (x @ centers.T), 0)


def align_to_hprc(reference: np.ndarray, candidate: np.ndarray) -> np.ndarray:
    candidate_values = np.unique(candidate)
    table = np.zeros((3, len(candidate_values)), dtype=np.int64)
    for i, ref in enumerate(HPRC_LABELS):
        for j, cand in enumerate(candidate_values):
            table[i, j] = np.sum((reference == ref) & (candidate == cand))
    rows, cols = linear_sum_assignment(-table)
    mapping = {candidate_values[col]: HPRC_LABELS[row] for row, col in zip(rows, cols)}
    return np.fromiter((mapping[value] for value in candidate), dtype=np.int32)


def crosswalk(source: np.ndarray, target: np.ndarray, source_name: str) -> pd.DataFrame:
    frame = pd.DataFrame({source_name: source, "hprc_cluster": target})
    result = (
        frame.value_counts(sort=False)
        .rename("units")
        .reset_index()
        .sort_values([source_name, "hprc_cluster"])
    )
    result["source_total"] = result.groupby(source_name).units.transform("sum")
    result["fraction_within_source"] = result.units / result.source_total
    return result


def compact_crosswalk(frame: pd.DataFrame, source_name: str) -> str:
    pivot = frame.pivot(index=source_name, columns="hprc_cluster", values="units").fillna(0)
    pivot = pivot.reindex(columns=HPRC_LABELS, fill_value=0).astype(int)
    lines = ["| Source | C1 | C2 | C3 | Dominant mapping |", "|---|---:|---:|---:|---|"]
    for source, row in pivot.iterrows():
        total = int(row.sum())
        best = int(row.idxmax())
        fraction = int(row.max()) / total
        lines.append(
            f"| {source} | {row.get(1, 0)} | {row.get(2, 0)} | {row.get(3, 0)} "
            f"| C{best} ({fraction:.1%}) |"
        )
    return "\n".join(lines)


def fit_joint(
    x: np.ndarray,
    hprc_n: int,
    hprc_reference: np.ndarray,
    init: np.ndarray | str,
    n_init: int,
    random_state: int,
) -> tuple[KMeans, np.ndarray]:
    model = KMeans(
        n_clusters=3,
        init=init,
        n_init=n_init,
        random_state=random_state,
        algorithm="lloyd",
        max_iter=300,
        tol=1e-4,
    ).fit(x)
    aligned = align_to_hprc(hprc_reference, model.labels_[:hprc_n])
    # Apply the same candidate->HPRC mapping to all rows by inferring it from
    # the already aligned HPRC subset.
    mapping: dict[int, int] = {}
    for raw_label in np.unique(model.labels_):
        values = aligned[model.labels_[:hprc_n] == raw_label]
        mapping[int(raw_label)] = int(np.bincount(values, minlength=4).argmax())
    all_aligned = np.fromiter(
        (mapping[int(value)] for value in model.labels_), dtype=np.int32
    )
    return model, all_aligned


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    hprc_meta, hprc_raw, hprc_features = read_features(HPRC_FEATURES)
    hg_meta, hg_raw, hg_features = read_features(HG002_FEATURES)
    if hprc_features != hg_features:
        raise ValueError("HPRC and HG002 feature columns differ")

    hprc_assign = pd.read_csv(HPRC_ASSIGNMENTS, sep="\t")
    hg_assign = pd.read_csv(HG002_ASSIGNMENTS, sep="\t")
    hprc_labels = ordered_labels(hprc_meta, hprc_assign, "selected_cluster").astype(np.int32)
    hg_labels = ordered_labels(hg_meta, hg_assign, "constrained_k4_cluster").astype(str)

    # Fixed HPRC-centroid projection.
    scaling = pd.read_csv(HPRC_SCALING, sep="\t")
    hprc_x = scale_saved(hprc_raw, scaling, hprc_features)
    hg_x = scale_saved(hg_raw, scaling, hprc_features)
    fixed_centers = centroids(hprc_x, hprc_labels)
    fixed_distances = distances_to_centroids(hg_x, fixed_centers)
    fixed_order = np.argsort(fixed_distances, axis=1)
    fixed_labels = fixed_order[:, 0].astype(np.int32) + 1
    fixed_margin = (
        fixed_distances[np.arange(len(hg_x)), fixed_order[:, 1]]
        - fixed_distances[np.arange(len(hg_x)), fixed_order[:, 0]]
    )
    fixed_crosswalk = crosswalk(hg_labels, fixed_labels, "hg002_cluster")
    fixed_crosswalk.to_csv(OUT / "fixed_projection_crosswalk.tsv", sep="\t", index=False)
    fixed_units = hg_meta[["acre_id", "array_id"]].copy()
    fixed_units["hg002_cluster"] = hg_labels
    fixed_units["projected_hprc_cluster"] = fixed_labels
    fixed_units["distance_to_nearest"] = np.sqrt(
        fixed_distances[np.arange(len(hg_x)), fixed_order[:, 0]]
    )
    fixed_units["distance_to_second"] = np.sqrt(
        fixed_distances[np.arange(len(hg_x)), fixed_order[:, 1]]
    )
    fixed_units["squared_distance_margin"] = fixed_margin
    fixed_units.to_csv(
        OUT / "fixed_projection_unit_assignments.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    margin_summary = (
        fixed_units.groupby("hg002_cluster", sort=True)
        .agg(
            units=("acre_id", "size"),
            median_distance_to_nearest=("distance_to_nearest", "median"),
            median_distance_to_second=("distance_to_second", "median"),
            mean_squared_distance_margin=("squared_distance_margin", "mean"),
            median_squared_distance_margin=("squared_distance_margin", "median"),
            q05_squared_distance_margin=(
                "squared_distance_margin", lambda values: values.quantile(0.05)
            ),
        )
        .reset_index()
    )
    margin_summary.to_csv(
        OUT / "fixed_projection_margin_summary.tsv", sep="\t", index=False
    )
    del hprc_x, hg_x, fixed_distances

    # Joint scaling and refits.
    hprc_n = len(hprc_raw)
    combined_raw = np.vstack([hprc_raw, hg_raw])
    del hprc_raw, hg_raw
    joint_x, usable = scale_joint(combined_raw)
    del combined_raw
    anchor_centers = centroids(joint_x[:hprc_n], hprc_labels)

    anchored_model, anchored_labels = fit_joint(
        joint_x, hprc_n, hprc_labels, anchor_centers, 1, 1
    )
    multistart_model, multistart_labels = fit_joint(
        joint_x, hprc_n, hprc_labels, "k-means++", 10, 1001
    )

    anchored_hprc = anchored_labels[:hprc_n]
    anchored_hg = anchored_labels[hprc_n:]
    multistart_hprc = multistart_labels[:hprc_n]
    multistart_hg = multistart_labels[hprc_n:]

    joint_crosswalk = crosswalk(hg_labels, anchored_hg, "hg002_cluster")
    joint_crosswalk.to_csv(OUT / "joint_refit_crosswalk.tsv", sep="\t", index=False)
    multistart_crosswalk = crosswalk(hg_labels, multistart_hg, "hg002_cluster")
    multistart_crosswalk.to_csv(
        OUT / "joint_multistart_crosswalk.tsv", sep="\t", index=False
    )
    hprc_retention = crosswalk(
        hprc_labels.astype(str), anchored_hprc, "original_hprc_cluster"
    )
    hprc_retention.to_csv(OUT / "hprc_retention_crosswalk.tsv", sep="\t", index=False)

    joint_units = hg_meta[["acre_id", "array_id"]].copy()
    joint_units["hg002_cluster"] = hg_labels
    joint_units["fixed_projection_hprc_cluster"] = fixed_labels
    joint_units["anchored_joint_hprc_cluster"] = anchored_hg
    joint_units["multistart_joint_hprc_cluster"] = multistart_hg
    joint_units.to_csv(
        OUT / "joint_unit_assignments.tsv.gz", sep="\t", index=False, compression="gzip"
    )

    summary = {
        "hprc_units": hprc_n,
        "hg002_units": len(hg_labels),
        "joint_nonconstant_features": int(usable.sum()),
        "fixed_projection_vs_anchored_joint_hg002_agreement": float(
            np.mean(fixed_labels == anchored_hg)
        ),
        "anchored_vs_multistart_joint_hg002_agreement": float(
            np.mean(anchored_hg == multistart_hg)
        ),
        "hprc_original_vs_anchored_joint_exact_retention": float(
            np.mean(hprc_labels == anchored_hprc)
        ),
        "hprc_original_vs_anchored_joint_units_retained": int(
            np.sum(hprc_labels == anchored_hprc)
        ),
        "hprc_original_vs_anchored_joint_ari": float(
            adjusted_rand_score(hprc_labels, anchored_hprc)
        ),
        "hprc_original_vs_multistart_joint_exact_retention": float(
            np.mean(hprc_labels == multistart_hprc)
        ),
        "hprc_original_vs_multistart_joint_units_retained": int(
            np.sum(hprc_labels == multistart_hprc)
        ),
        "hprc_original_vs_multistart_joint_ari": float(
            adjusted_rand_score(hprc_labels, multistart_hprc)
        ),
        "anchored_joint_inertia": float(anchored_model.inertia_),
        "multistart_joint_inertia": float(multistart_model.inertia_),
        "anchored_joint_iterations": int(anchored_model.n_iter_),
        "multistart_joint_iterations": int(multistart_model.n_iter_),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    report = f"""# HG002-to-HPRC NTRiage cluster crosswalk

This analysis tests whether the current HG002 H1-H4 assignments map onto the
authoritative merged-HPRC K=3 solution. The HPRC analysis contains
{hprc_n:,} ACREs; the acrocentric-only HG002 panel contributes {len(hg_labels):,}.

## Fixed HPRC-centroid projection

{compact_crosswalk(fixed_crosswalk, 'hg002_cluster')}

## Joint HPRC + HG002 K=3 refit

{compact_crosswalk(joint_crosswalk, 'hg002_cluster')}

The fixed projection and anchored joint refit agree for
{summary['fixed_projection_vs_anchored_joint_hg002_agreement']:.1%} of HG002
ACREs. Adding HG002 retains
{summary['hprc_original_vs_anchored_joint_units_retained']:,}/{hprc_n:,}
({summary['hprc_original_vs_anchored_joint_exact_retention']:.4%}) of the original
HPRC assignments (ARI={summary['hprc_original_vs_anchored_joint_ari']:.5f}).
The independent 10-start joint fit agrees with the anchored joint result for
{summary['anchored_vs_multistart_joint_hg002_agreement']:.1%} of HG002 ACREs and
retains {summary['hprc_original_vs_multistart_joint_units_retained']:,}/{hprc_n:,}
({summary['hprc_original_vs_multistart_joint_exact_retention']:.4%}) of the
original HPRC assignments
(ARI={summary['hprc_original_vs_multistart_joint_ari']:.5f}).

Cluster numbers refer to the existing HPRC labels: C1 is the mixed S1/S3
cluster, C2 is the S3-rich cluster, and C3 is the S2-rich cluster.
"""
    (HERE / "RESULTS.md").write_text(report)
    print(report)


if __name__ == "__main__":
    main()
