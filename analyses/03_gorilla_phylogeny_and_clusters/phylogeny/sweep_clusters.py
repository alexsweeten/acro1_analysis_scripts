#!/usr/bin/env python3
"""De novo multi-start k-means cluster-number sweep for gorilla ACREs."""

from __future__ import annotations

import json
import os
from itertools import combinations
from pathlib import Path

HERE = Path(__file__).resolve().parent
os.environ.setdefault("MPLCONFIGDIR", str(HERE / "matplotlib_cache"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sklearn
from sklearn.cluster import KMeans
from sklearn.metrics import (
    adjusted_rand_score,
    calinski_harabasz_score,
    davies_bouldin_score,
    silhouette_score,
)


FEATURES = HERE / "gorilla.features.tsv"
K_VALUES = list(range(2, 17))
SEEDS = list(range(1, 11))
CLUSTER_COLORS = [
    "#4E79A7", "#F28E2B", "#59A14F", "#E15759", "#B07AA1",
    "#76B7B2", "#EDC948", "#9C755F", "#FF9DA7", "#BAB0AC",
    "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02",
]
CHROMOSOME_ORDER = [
    "hsa2a", "hsa2b", "hsa8", "hsa9", "hsa13", "hsa15",
    "hsa18", "hsa19", "hsa21", "hsa22", "hsaY",
]
CHROMOSOME_COLORS = [
    "#1B9E77", "#66A61E", "#A6D854", "#E6AB02", "#D95F02", "#E7298A",
    "#7570B3", "#A6761D", "#1F78B4", "#6BAED6", "#202020",
]


def stable_labels(labels: np.ndarray) -> np.ndarray:
    first = [(np.flatnonzero(labels == label)[0], label) for label in np.unique(labels)]
    remap = {old: new + 1 for new, (_, old) in enumerate(sorted(first))}
    return np.fromiter((remap[x] for x in labels), dtype=np.int32, count=len(labels))


def scale_features(raw: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    means = raw.mean(axis=0)
    scales = raw.std(axis=0, ddof=1)
    usable = scales > 1e-12
    transformed = np.zeros_like(raw, dtype=np.float64)
    transformed[:, usable] = (raw[:, usable] - means[usable]) / scales[usable]
    norms = np.linalg.norm(transformed, axis=1)
    transformed[norms > 0] /= norms[norms > 0, None]
    return transformed, means, scales


def minmax(values: pd.Series, higher_is_better: bool = True) -> pd.Series:
    low, high = values.min(), values.max()
    result = pd.Series(np.ones(len(values)), index=values.index)
    if high != low:
        result = (values - low) / (high - low)
    return result if higher_is_better else 1.0 - result


def save_sweep_figure(summary: pd.DataFrame, best_k: int) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), sharex=True)
    panels = [
        ("silhouette", "Silhouette (higher is better)"),
        ("davies_bouldin", "Davies–Bouldin (lower is better)"),
        ("stability_ari_mean", "Mean multi-start stability ARI"),
        ("consensus_score", "Equal-weight consensus score"),
    ]
    for ax, (column, title) in zip(axes.flat, panels):
        ax.plot(summary["k"], summary[column], marker="o", linewidth=2, color="#4E79A7")
        ax.axvline(best_k, color="#D73027", linestyle="--", linewidth=1.5)
        ax.scatter(
            [best_k], summary.loc[summary["k"] == best_k, column],
            color="#D73027", s=70, zorder=4,
        )
        ax.set_title(title)
        ax.grid(alpha=0.25)
    for ax in axes[-1, :]:
        ax.set_xlabel("Number of clusters (k)")
    fig.suptitle(
        f"Gorilla-only de novo cluster sweep; selected k = {best_k}",
        fontsize=15,
        fontweight="bold",
    )
    fig.tight_layout()
    fig.savefig(HERE / "cluster_sweep_metrics.png", dpi=240, bbox_inches="tight", facecolor="white")
    fig.savefig(HERE / "cluster_sweep_metrics.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> None:
    frame = pd.read_csv(FEATURES, sep="\t")
    id_columns = ["acre_id", "array_id"]
    feature_names = [column for column in frame.columns if column not in id_columns]
    raw = frame[feature_names].to_numpy(dtype=np.float64, copy=True)
    x, means, scales = scale_features(raw)

    pd.DataFrame(
        {
            "feature": feature_names,
            "mean": means,
            "sample_sd": scales,
            "used": scales > 1e-12,
        }
    ).to_csv(HERE / "cluster_sweep_feature_scaling.tsv", sep="\t", index=False)

    total_ss = float(np.square(x - x.mean(axis=0)).sum())
    all_runs: list[dict[str, object]] = []
    summaries: list[dict[str, object]] = []
    assignments = frame[id_columns].copy()
    centroid_rows: list[dict[str, object]] = []

    for k in K_VALUES:
        print(f"k={k}: running {len(SEEDS)} independent starts", flush=True)
        results: list[tuple[int, KMeans, np.ndarray]] = []
        for seed in SEEDS:
            model = KMeans(
                n_clusters=k,
                init="k-means++",
                n_init=1,
                max_iter=300,
                tol=1e-4,
                random_state=seed,
                algorithm="lloyd",
            ).fit(x)
            labels = stable_labels(model.labels_)
            counts = np.bincount(labels, minlength=k + 1)[1:]
            all_runs.append(
                {
                    "k": k,
                    "seed": seed,
                    "inertia": model.inertia_,
                    "iterations": model.n_iter_,
                    "min_cluster_size": counts.min(),
                    "max_cluster_size": counts.max(),
                    "cluster_size_cv": counts.std(ddof=1) / counts.mean(),
                }
            )
            results.append((seed, model, labels))

        best_seed, best_model, best_labels = min(results, key=lambda item: item[1].inertia_)
        counts = np.bincount(best_labels, minlength=k + 1)[1:]
        pairwise_ari = [
            adjusted_rand_score(left[2], right[2])
            for left, right in combinations(results, 2)
        ]
        summaries.append(
            {
                "k": k,
                "best_seed": best_seed,
                "inertia": best_model.inertia_,
                "explained_variance_fraction": 1.0 - best_model.inertia_ / total_ss,
                "silhouette": silhouette_score(x, best_labels, metric="euclidean"),
                "calinski_harabasz": calinski_harabasz_score(x, best_labels),
                "davies_bouldin": davies_bouldin_score(x, best_labels),
                "stability_ari_mean": np.mean(pairwise_ari),
                "stability_ari_sd": np.std(pairwise_ari, ddof=1),
                "stability_ari_min": np.min(pairwise_ari),
                "min_cluster_size": counts.min(),
                "max_cluster_size": counts.max(),
                "cluster_size_cv": counts.std(ddof=1) / counts.mean(),
            }
        )
        assignments[f"k{k:02d}_cluster"] = best_labels
        for cluster in range(1, k + 1):
            members = x[best_labels == cluster]
            center = members.mean(axis=0)
            row: dict[str, object] = {
                "k": k,
                "cluster": cluster,
                "size": len(members),
            }
            row.update(dict(zip(feature_names, center)))
            centroid_rows.append(row)

    summary = pd.DataFrame(summaries).sort_values("k").reset_index(drop=True)
    inertia = summary["inertia"].to_numpy()
    relative_drop = np.full(len(summary), np.nan)
    relative_drop[1:] = 1.0 - inertia[1:] / inertia[:-1]
    summary["relative_inertia_drop"] = relative_drop
    elbow = np.full(len(summary), np.nan)
    for index in range(1, len(summary) - 1):
        prior_gain = 1.0 - summary.loc[index, "inertia"] / summary.loc[index - 1, "inertia"]
        next_gain = 1.0 - summary.loc[index + 1, "inertia"] / summary.loc[index, "inertia"]
        elbow[index] = prior_gain - next_gain
    summary["elbow_strength"] = elbow
    summary["score_silhouette"] = minmax(summary["silhouette"])
    summary["score_calinski_harabasz"] = minmax(summary["calinski_harabasz"])
    summary["score_davies_bouldin"] = minmax(summary["davies_bouldin"], False)
    summary["score_stability"] = minmax(summary["stability_ari_mean"])
    score_columns = [
        "score_silhouette",
        "score_calinski_harabasz",
        "score_davies_bouldin",
        "score_stability",
    ]
    summary["consensus_score"] = summary[score_columns].mean(axis=1)
    summary["consensus_rank"] = summary["consensus_score"].rank(
        method="min", ascending=False
    ).astype(int)
    best_k = int(
        summary.sort_values(["consensus_score", "k"], ascending=[False, True]).iloc[0]["k"]
    )

    runs = pd.DataFrame(all_runs)
    runs["selected_best_start"] = False
    for row in summaries:
        runs.loc[
            (runs["k"] == row["k"]) & (runs["seed"] == row["best_seed"]),
            "selected_best_start",
        ] = True
    runs.to_csv(HERE / "cluster_sweep_runs.tsv", sep="\t", index=False)
    summary.to_csv(HERE / "cluster_sweep_summary.tsv", sep="\t", index=False)
    assignments.to_csv(HERE / "cluster_sweep_assignments.tsv", sep="\t", index=False)
    pd.DataFrame(centroid_rows).to_csv(
        HERE / "cluster_sweep_centroids.tsv", sep="\t", index=False
    )

    labels = assignments[f"k{best_k:02d}_cluster"].to_numpy(dtype=int)
    chosen_centers = np.vstack([x[labels == cluster].mean(axis=0) for cluster in range(1, best_k + 1)])
    distances = np.linalg.norm(x - chosen_centers[labels - 1], axis=1)
    chosen = frame[id_columns].copy()
    chosen["de_novo_cluster"] = labels
    chosen["centroid_distance"] = distances
    chosen.to_csv(HERE / "gorilla_de_novo_clusters.tsv", sep="\t", index=False)

    metadata = pd.read_csv(HERE / "gorilla_tip_metadata.base.tsv", sep="\t")
    metadata = metadata.merge(chosen, on=["acre_id", "array_id"], validate="one_to_one")
    metadata.to_csv(HERE / "gorilla_tip_metadata.tsv", sep="\t", index=False)

    cluster_summary = (
        metadata.groupby("de_novo_cluster", observed=True)
        .agg(
            acres=("acre_id", "size"),
            parent_arrays=("array_id", "nunique"),
            chromosome_groups=("synteny_group", "nunique"),
            mean_length=("length", "mean"),
            mean_periodicity_score=("periodicity_score", "mean"),
            mean_component_coverage=("component_coverage", "mean"),
            mean_exclusive_fraction=("exclusive_fraction", "mean"),
        )
        .reset_index()
    )
    top_components = (
        metadata.groupby("de_novo_cluster")["top_component"]
        .agg(lambda values: values.value_counts().index[0])
        .rename("most_common_top_component")
    )
    cluster_summary = cluster_summary.merge(top_components, on="de_novo_cluster")
    cluster_summary.to_csv(HERE / "cluster_summary.tsv", sep="\t", index=False)
    pd.crosstab(metadata["de_novo_cluster"], metadata["synteny_group"]).to_csv(
        HERE / "cluster_chromosome_counts.tsv", sep="\t"
    )

    pd.DataFrame(
        {
            "de_novo_cluster": range(1, best_k + 1),
            "color": CLUSTER_COLORS[:best_k],
        }
    ).to_csv(HERE / "cluster_palette.tsv", sep="\t", index=False)
    present_chromosomes = [
        chromosome for chromosome in CHROMOSOME_ORDER
        if chromosome in set(metadata["synteny_group"])
    ]
    pd.DataFrame(
        {
            "synteny_group": present_chromosomes,
            "color": CHROMOSOME_COLORS[: len(present_chromosomes)],
        }
    ).to_csv(HERE / "chromosome_palette.tsv", sep="\t", index=False)

    parameters = {
        "input": FEATURES.name,
        "observations": len(frame),
        "raw_features": len(feature_names),
        "features_with_nonzero_variance": int((scales > 1e-12).sum()),
        "preprocessing": "gorilla-only sample-SD column standardization followed by row L2 normalization",
        "algorithm": "scikit-learn KMeans, k-means++, Lloyd",
        "k_values": K_VALUES,
        "seeds": SEEDS,
        "n_init_per_seed": 1,
        "max_iter": 300,
        "tolerance": 1e-4,
        "silhouette_sample_size": len(frame),
        "selection_rule": (
            "highest equal-weight mean of min-max-scaled silhouette, Calinski-Harabasz, "
            "inverse Davies-Bouldin, and mean pairwise ARI stability"
        ),
        "best_k": best_k,
        "numpy_version": np.__version__,
        "pandas_version": pd.__version__,
        "scikit_learn_version": sklearn.__version__,
    }
    (HERE / "cluster_sweep_parameters.json").write_text(json.dumps(parameters, indent=2) + "\n")
    save_sweep_figure(summary, best_k)
    best = summary.loc[summary["k"] == best_k].iloc[0]
    print(
        f"Selected gorilla-only k={best_k}: silhouette={best.silhouette:.4f}, "
        f"stability ARI={best.stability_ari_mean:.4f}, "
        f"explained variance={best.explained_variance_fraction:.2%}"
    )


if __name__ == "__main__":
    main()
