#!/usr/bin/env python3
"""De novo gorilla-only N-TRIAGE clustering and NOR composition figures."""

from __future__ import annotations

import json
from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import chi2_contingency
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.metrics import (
    adjusted_rand_score,
    calinski_harabasz_score,
    davies_bouldin_score,
    silhouette_score,
)


ROOT = Path(__file__).resolve().parents[1]
OUT = Path(__file__).resolve().parent
FEATURES = ROOT / "acro_gorilla_only/gorilla.features.tsv"
METADATA = ROOT / "acro_gorilla_only/gorilla_tip_metadata.base.tsv"
K_VALUES = list(range(2, 17))
SEEDS = list(range(1, 11))
SEED = 20260923
COLORS = {
    "G1": "#4E79A7",
    "G2": "#F28E2B",
    "G3": "#59A14F",
    "G4": "#E15759",
    "G5": "#B07AA1",
    "G6": "#76B7B2",
    "G7": "#EDC948",
    "G8": "#9C755F",
}


def stable_labels(labels):
    first = [(np.flatnonzero(labels == label)[0], label) for label in np.unique(labels)]
    remap = {old: new + 1 for new, (_, old) in enumerate(sorted(first))}
    return np.fromiter((remap[value] for value in labels), dtype=np.int32, count=len(labels))


def scale_features(raw):
    means = raw.mean(axis=0)
    scales = raw.std(axis=0, ddof=1)
    usable = scales > 1e-12
    transformed = np.zeros_like(raw, dtype=float)
    transformed[:, usable] = (raw[:, usable] - means[usable]) / scales[usable]
    norms = np.linalg.norm(transformed, axis=1)
    transformed[norms > 0] /= norms[norms > 0, None]
    return transformed, means, scales


def minmax(values, higher=True):
    low, high = values.min(), values.max()
    result = pd.Series(np.ones(len(values)), index=values.index)
    if high != low:
        result = (values - low) / (high - low)
    return result if higher else 1 - result


def cramer_v(table):
    chi2, p_value, dof, _ = chi2_contingency(table, correction=False)
    n = table.to_numpy().sum()
    denominator = n * min(table.shape[0] - 1, table.shape[1] - 1)
    return chi2, p_value, dof, np.sqrt(chi2 / denominator)


def plot_sweep(summary, best_k):
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.5), sharex=True, constrained_layout=True)
    panels = [
        ("silhouette", "Silhouette (higher is better)"),
        ("davies_bouldin", "Davies–Bouldin (lower is better)"),
        ("stability_ari_mean", "Multi-start stability ARI"),
        ("consensus_score", "Consensus selection score"),
    ]
    for axis, (column, title) in zip(axes.flat, panels):
        axis.plot(summary["k"], summary[column], marker="o", color="#4E79A7", linewidth=2)
        axis.axvline(best_k, color="#D1495B", linestyle="--", linewidth=1.5)
        chosen = summary[summary["k"] == best_k]
        axis.scatter(chosen["k"], chosen[column], color="#D1495B", s=75, zorder=4)
        axis.set_title(title, fontweight="bold")
        axis.grid(alpha=0.25)
    for axis in axes[-1]:
        axis.set_xlabel("Number of clusters (k)")
    fig.suptitle(f"Gorilla-only cluster-number sweep: selected k={best_k}", fontsize=16, fontweight="bold")
    fig.savefig(OUT / "gorilla_cluster_selection.png", dpi=280, bbox_inches="tight")
    fig.savefig(OUT / "gorilla_cluster_selection.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_embedding(frame, x, y, x_label, y_label, stem, title):
    fig, axis = plt.subplots(figsize=(10.8, 8.2), constrained_layout=True)
    order = sorted(frame["gorilla_cluster"].unique(), key=lambda value: int(value[1:]))
    for cluster in order:
        subset = frame[frame["gorilla_cluster"] == cluster]
        axis.scatter(
            subset[x],
            subset[y],
            s=25,
            alpha=0.72,
            color=COLORS[cluster],
            label=f"{cluster} (n={len(subset):,})",
            edgecolors="none",
        )
        axis.scatter(
            subset[x].median(),
            subset[y].median(),
            s=175,
            marker="*",
            color=COLORS[cluster],
            edgecolor="black",
            linewidth=0.9,
            zorder=5,
        )
        axis.annotate(
            cluster,
            (subset[x].median(), subset[y].median()),
            xytext=(0, 11),
            textcoords="offset points",
            fontsize=10,
            fontweight="bold",
            ha="center",
            va="bottom",
            zorder=6,
        )
    axis.set_xlabel(x_label)
    axis.set_ylabel(y_label)
    axis.set_title(title, fontsize=16, fontweight="bold", pad=12)
    axis.legend(title="Gorilla-only cluster", loc="upper right", frameon=True)
    axis.grid(alpha=0.2)
    sns.despine(ax=axis)
    fig.savefig(OUT / f"{stem}.png", dpi=280, bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_nor_composition(composition, clusters):
    fig, axis = plt.subplots(figsize=(13.8, 7.4), constrained_layout=True)
    statuses = ["NOR−", "NOR+"]
    left = np.zeros(2)
    for cluster in clusters:
        values = np.array(
            [
                composition.loc[(composition["nor_status"] == status) & (composition["gorilla_cluster"] == cluster), "fraction_of_status_acres"].iloc[0]
                for status in statuses
            ]
        )
        bars = axis.barh(
            statuses,
            values,
            left=left,
            color=COLORS[cluster],
            edgecolor="white",
            linewidth=1.5,
            label=cluster,
            height=0.56,
        )
        for bar, value, start in zip(bars, values, left):
            if value >= 0.05:
                axis.text(
                    start + value / 2,
                    bar.get_y() + bar.get_height() / 2,
                    f"{cluster}\n{100 * value:.1f}%",
                    ha="center",
                    va="center",
                    color="white",
                    fontsize=11,
                    fontweight="bold",
                )
        left += values
    axis.set_xlim(0, 1)
    axis.set_xlabel("Fraction of complete ACRE units")
    axis.set_ylabel("")
    axis.set_title(
        "Gorilla-only cluster composition of NOR+ and NOR− ACREs",
        fontsize=17,
        fontweight="bold",
        pad=14,
    )
    axis.legend(
        title="De novo gorilla cluster",
        loc="upper center",
        bbox_to_anchor=(0.5, -0.15),
        ncol=len(clusters),
        frameon=False,
    )
    axis.grid(axis="x", alpha=0.25)
    axis.set_axisbelow(True)
    sns.despine(ax=axis)
    fig.savefig(OUT / "gorilla_nor_cluster_composition.png", dpi=280, bbox_inches="tight")
    fig.savefig(OUT / "gorilla_nor_cluster_composition.pdf", bbox_inches="tight")
    plt.close(fig)


def main():
    sns.set_theme(style="whitegrid", context="talk")
    features = pd.read_csv(FEATURES, sep="\t")
    metadata = pd.read_csv(METADATA, sep="\t")
    if set(features["acre_id"]) != set(metadata["acre_id"]):
        raise ValueError("feature and metadata ACRE identifiers differ")
    metadata = features[["acre_id", "array_id"]].merge(
        metadata, on=["acre_id", "array_id"], validate="one_to_one"
    )
    feature_names = [column for column in features if column not in {"acre_id", "array_id"}]
    raw = features[feature_names].to_numpy(dtype=float)
    x, means, scales = scale_features(raw)
    pd.DataFrame(
        {"feature": feature_names, "mean": means, "sample_sd": scales, "used": scales > 1e-12}
    ).to_csv(OUT / "feature_scaling.tsv", sep="\t", index=False)

    total_ss = np.square(x - x.mean(axis=0)).sum()
    runs = []
    summaries = []
    assignments = features[["acre_id", "array_id"]].copy()
    best_models = {}
    for k in K_VALUES:
        print(f"k={k}: {len(SEEDS)} starts", flush=True)
        results = []
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
            results.append((seed, model, labels))
            runs.append(
                {
                    "k": k,
                    "seed": seed,
                    "inertia": model.inertia_,
                    "iterations": model.n_iter_,
                    "min_cluster_size": counts.min(),
                    "max_cluster_size": counts.max(),
                }
            )
        best_seed, best_model, best_labels = min(results, key=lambda value: value[1].inertia_)
        best_models[k] = (best_seed, best_model, best_labels)
        counts = np.bincount(best_labels, minlength=k + 1)[1:]
        ari = [adjusted_rand_score(a[2], b[2]) for a, b in combinations(results, 2)]
        summaries.append(
            {
                "k": k,
                "best_seed": best_seed,
                "inertia": best_model.inertia_,
                "explained_variance_fraction": 1 - best_model.inertia_ / total_ss,
                "silhouette": silhouette_score(x, best_labels),
                "calinski_harabasz": calinski_harabasz_score(x, best_labels),
                "davies_bouldin": davies_bouldin_score(x, best_labels),
                "stability_ari_mean": np.mean(ari),
                "stability_ari_sd": np.std(ari, ddof=1),
                "stability_ari_min": np.min(ari),
                "min_cluster_size": counts.min(),
                "max_cluster_size": counts.max(),
                "cluster_size_cv": counts.std(ddof=1) / counts.mean(),
            }
        )
        assignments[f"k{k:02d}_cluster"] = best_labels

    summary = pd.DataFrame(summaries)
    summary["score_silhouette"] = minmax(summary["silhouette"])
    summary["score_calinski_harabasz"] = minmax(summary["calinski_harabasz"])
    summary["score_davies_bouldin"] = minmax(summary["davies_bouldin"], False)
    summary["score_stability"] = minmax(summary["stability_ari_mean"])
    summary["consensus_score"] = summary[
        ["score_silhouette", "score_calinski_harabasz", "score_davies_bouldin", "score_stability"]
    ].mean(axis=1)
    summary["consensus_rank"] = summary["consensus_score"].rank(method="min", ascending=False).astype(int)
    best_k = int(summary.sort_values(["consensus_score", "k"], ascending=[False, True]).iloc[0]["k"])
    summary.to_csv(OUT / "cluster_sweep_summary.tsv", sep="\t", index=False)
    pd.DataFrame(runs).to_csv(OUT / "cluster_sweep_runs.tsv", sep="\t", index=False)
    assignments.to_csv(OUT / "cluster_sweep_assignments.tsv", sep="\t", index=False)
    plot_sweep(summary, best_k)

    labels = best_models[best_k][2]
    cluster_names = np.array([f"G{value}" for value in labels])
    pca = PCA(n_components=min(50, x.shape[0], x.shape[1]), random_state=SEED)
    pcs = pca.fit_transform(x)
    tsne = TSNE(
        n_components=2,
        perplexity=30,
        learning_rate="auto",
        init="pca",
        max_iter=2000,
        random_state=SEED,
    ).fit_transform(pcs)

    result = metadata.copy()
    result["gorilla_cluster"] = cluster_names
    result["selected_k"] = best_k
    result["PC1"] = pcs[:, 0]
    result["PC2"] = pcs[:, 1]
    result["TSNE1"] = tsne[:, 0]
    result["TSNE2"] = tsne[:, 1]
    result.to_csv(OUT / "gorilla_cluster_assignments.tsv.gz", sep="\t", index=False)

    joint = pd.read_csv(ROOT / "primate_full/six_primate_acre_metadata.tsv", sep="\t")
    joint = joint[joint["species_name"] == "gorilla"][["acre_id", "cluster8"]]
    crosswalk_data = result[["acre_id", "gorilla_cluster"]].merge(
        joint, on="acre_id", validate="one_to_one"
    )
    crosswalk = (
        crosswalk_data.groupby(["gorilla_cluster", "cluster8"])
        .size()
        .rename("acres")
        .reset_index()
    )
    crosswalk["fraction_of_gorilla_cluster"] = crosswalk["acres"] / crosswalk.groupby(
        "gorilla_cluster"
    )["acres"].transform("sum")
    crosswalk["fraction_of_joint_cluster"] = crosswalk["acres"] / crosswalk.groupby(
        "cluster8"
    )["acres"].transform("sum")
    crosswalk.to_csv(OUT / "gorilla_cluster_joint_C1_C8_crosswalk.tsv", sep="\t", index=False)
    pd.DataFrame(
        {
            "principal_component": np.arange(1, len(pca.explained_variance_ratio_) + 1),
            "explained_variance_ratio": pca.explained_variance_ratio_,
            "cumulative_explained_variance": np.cumsum(pca.explained_variance_ratio_),
        }
    ).to_csv(OUT / "pca_explained_variance.tsv", sep="\t", index=False)
    pd.DataFrame(
        pca.components_.T,
        index=feature_names,
        columns=[f"PC{i}" for i in range(1, len(pca.components_) + 1)],
    ).rename_axis("feature").reset_index().to_csv(OUT / "pca_loadings.tsv.gz", sep="\t", index=False)

    plot_embedding(
        result,
        "PC1",
        "PC2",
        f"PC1 ({100 * pca.explained_variance_ratio_[0]:.1f}% variance)",
        f"PC2 ({100 * pca.explained_variance_ratio_[1]:.1f}% variance)",
        "gorilla_optimal_clusters_pca",
        f"Gorilla-only N-TRIAGE clustering (optimal k={best_k}): PCA",
    )
    plot_embedding(
        result,
        "TSNE1",
        "TSNE2",
        "t-SNE 1",
        "t-SNE 2",
        "gorilla_optimal_clusters_tsne",
        f"Gorilla-only N-TRIAGE clustering (optimal k={best_k}): t-SNE",
    )

    cluster_summary = (
        result.groupby("gorilla_cluster")
        .agg(
            acres=("acre_id", "size"),
            parent_arrays=("array_id", "nunique"),
            synteny_groups=("synteny_group", "nunique"),
            median_length_bp=("length", "median"),
            mean_periodicity_score=("periodicity_score", "mean"),
            nor_positive_acres=("nor_status", lambda values: int((values == "NOR+").sum())),
            nor_negative_acres=("nor_status", lambda values: int((values == "NOR-").sum())),
        )
        .reset_index()
    )
    cluster_summary["nor_positive_fraction_current_acro"] = (
        cluster_summary["nor_positive_acres"]
        / (cluster_summary["nor_positive_acres"] + cluster_summary["nor_negative_acres"])
    )
    cluster_summary.to_csv(OUT / "gorilla_cluster_summary.tsv", sep="\t", index=False)

    current = result[result["nor_status"].isin(["NOR+", "NOR-"])].copy()
    current["nor_status"] = current["nor_status"].replace({"NOR-": "NOR−"})
    clusters = sorted(current["gorilla_cluster"].unique(), key=lambda value: int(value[1:]))
    rows = []
    for status in ("NOR−", "NOR+"):
        status_frame = current[current["nor_status"] == status]
        for cluster in clusters:
            subset = status_frame[status_frame["gorilla_cluster"] == cluster]
            rows.append(
                {
                    "nor_status": status,
                    "gorilla_cluster": cluster,
                    "cluster_acres": len(subset),
                    "status_acres": len(status_frame),
                    "fraction_of_status_acres": len(subset) / len(status_frame),
                    "arrays_with_cluster": subset["array_id"].nunique(),
                }
            )
    composition = pd.DataFrame(rows)
    composition.to_csv(OUT / "gorilla_nor_cluster_composition.tsv", sep="\t", index=False)
    plot_nor_composition(composition, clusters)

    acre_table = pd.crosstab(current["gorilla_cluster"], current["nor_status"])
    acre_chi2, acre_p, acre_dof, acre_v = cramer_v(acre_table)
    array_status = current.groupby("array_id")["nor_status"].first()
    array_counts = pd.crosstab(current["array_id"], current["gorilla_cluster"])
    dominant = array_counts.idxmax(axis=1)
    array_table = pd.crosstab(dominant, array_status.loc[dominant.index])
    array_chi2, array_p, array_dof, array_v = cramer_v(array_table)
    pd.DataFrame(
        [
            {
                "analysis_level": "ACRE",
                "observations": len(current),
                "chi_square": acre_chi2,
                "degrees_of_freedom": acre_dof,
                "p_value": acre_p,
                "cramers_v": acre_v,
            },
            {
                "analysis_level": "dominant_cluster_per_array",
                "observations": len(dominant),
                "chi_square": array_chi2,
                "degrees_of_freedom": array_dof,
                "p_value": array_p,
                "cramers_v": array_v,
            },
        ]
    ).to_csv(OUT / "gorilla_nor_cluster_association.tsv", sep="\t", index=False)

    parameters = {
        "input_features": str(FEATURES.relative_to(ROOT)),
        "input_metadata": str(METADATA.relative_to(ROOT)),
        "acres": len(features),
        "raw_features": len(feature_names),
        "features_with_nonzero_variance": int((scales > 1e-12).sum()),
        "preprocessing": "sample-SD feature standardization followed by row L2 normalization",
        "clustering": "KMeans on full normalized feature matrix",
        "k_values": K_VALUES,
        "seeds": SEEDS,
        "selection_rule": "equal-weight min-max consensus of silhouette, Calinski-Harabasz, inverse Davies-Bouldin, and multi-start ARI stability",
        "selected_k": best_k,
        "visualizations": {
            "PCA": "linear projection only; not used for clustering",
            "tSNE": "nonlinear projection only; not used for clustering",
        },
    }
    (OUT / "analysis_parameters.json").write_text(json.dumps(parameters, indent=2) + "\n")
    print(f"selected k={best_k}")


if __name__ == "__main__":
    main()
