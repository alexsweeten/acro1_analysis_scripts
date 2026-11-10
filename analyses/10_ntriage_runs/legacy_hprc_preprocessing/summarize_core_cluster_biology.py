#!/usr/bin/env python3
"""Summarize biological signatures of the corrected core-only N-TRIAGE clusters.

The diagnostic-repeat columns were inputs to N-TRIAGE clustering.  Their
cluster contrasts are therefore descriptive labels for the fitted clusters,
not independent tests of enrichment.  No p-values are reported here.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(os.environ.get("ACRO1_LEGACY_HPRC_ROOT", Path(__file__).resolve().parents[1]))
ANALYSIS = ROOT / "manuscript_output" / "telomeric_ntriage_2794"
PREFIX = ANALYSIS / "native" / "telomeric_2794_core"
TABLES = ANALYSIS / "tables"
MANIFEST = ROOT / "manifest.validated.tsv"

CLUSTERS = Path(f"{PREFIX}.clusters.tsv")
FEATURES = Path(f"{PREFIX}.features.tsv")
ACRO1 = Path(f"{PREFIX}.acro1_composition.tsv")

SCOPE_NOTE = (
    "descriptive_not_independent: diagnostic repeat features contributed to "
    "the clustering model"
)


def short_name(column: str) -> str:
    """Return a manuscript-readable repeat label."""
    name = column.removeprefix("diagnostic_frequency:")
    return name.removeprefix("COMP-subunit_ACRO_rnd-5_")


def standardized_difference(a: pd.Series, b: pd.Series) -> float:
    """Difference divided by the unweighted RMS within-group SD."""
    denominator = np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2.0)
    if denominator == 0 or not np.isfinite(denominator):
        return np.nan
    return float((a.mean() - b.mean()) / denominator)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    clusters = pd.read_csv(CLUSTERS, sep="\t")
    features = pd.read_csv(FEATURES, sep="\t")
    acro1 = pd.read_csv(ACRO1, sep="\t")
    manifest = pd.read_csv(MANIFEST, sep="\t", usecols=["array_id", "individual_id"])

    if len(clusters) != 77_276:
        raise ValueError(f"Expected 77,276 core units; found {len(clusters):,}")
    if not clusters.acre_id.equals(features.acre_id):
        raise ValueError("Cluster and feature rows are not in identical ACRE order")
    if clusters.acre_id.duplicated().any() or acro1.acre_id.duplicated().any():
        raise ValueError("ACRE identifiers must be unique")

    individual_for_array = manifest.set_index("array_id").individual_id
    clusters["individual_id"] = clusters.source_name.map(individual_for_array)
    if clusters.individual_id.isna().any():
        raise ValueError("Some cluster source names did not join to the manifest")

    # Give every complete unit one mutually exclusive position label.  A
    # one-unit array is classified as first; there were no such arrays here.
    last_index = clusters.groupby("array_id").unit_in_array.transform("max")
    clusters["position_class"] = np.select(
        [clusters.unit_in_array.eq(1), clusters.unit_in_array.eq(last_index)],
        ["first", "last"],
        default="interior",
    )

    # Long-form top-component composition.  This is the identity of the
    # single component with the greatest component coverage in each unit.
    component = (
        clusters.groupby(["cluster", "top_component"], observed=True)
        .size()
        .rename("units")
        .reset_index()
    )
    component["within_cluster_fraction"] = component.units / component.groupby(
        "cluster"
    ).units.transform("sum")
    component["rank_within_cluster"] = component.groupby("cluster")[
        "within_cluster_fraction"
    ].rank(method="first", ascending=False).astype(int)
    component = component.sort_values(
        ["cluster", "rank_within_cluster", "top_component"]
    )
    component.insert(0, "cluster_label", "C" + component.cluster.astype(str))
    component["top_component_definition"] = (
        "component with greatest annotated coverage for each complete ACRO unit"
    )
    component.to_csv(TABLES / "core_cluster_top_components.tsv", sep="\t", index=False)

    # Long-form, annotated diagnostic feature contrasts.  Retain presence
    # frequencies so sparse, highly specific annotations can be recognized.
    diagnostic_columns = [
        column for column in features if column.startswith("diagnostic_frequency:")
    ]
    diagnostic_rows: list[dict[str, object]] = []
    for cluster in sorted(features.cluster.unique()):
        inside = features.cluster.eq(cluster)
        for column in diagnostic_columns:
            a = features.loc[inside, column]
            b = features.loc[~inside, column]
            a_mean = float(a.mean())
            b_mean = float(b.mean())
            effect = standardized_difference(a, b)
            diagnostic_rows.append(
                {
                    "cluster_label": f"C{cluster}",
                    "cluster": int(cluster),
                    "feature": short_name(column),
                    "native_feature": column,
                    "cluster_mean_diagnostic_kmer_frequency": a_mean,
                    "other_clusters_mean_diagnostic_kmer_frequency": b_mean,
                    "absolute_difference": a_mean - b_mean,
                    "fold_enrichment": a_mean / b_mean if b_mean > 0 else np.inf,
                    "standardized_difference": effect,
                    "cluster_nonzero_unit_fraction": float(a.gt(0).mean()),
                    "other_clusters_nonzero_unit_fraction": float(b.gt(0).mean()),
                    "interpretation_scope": SCOPE_NOTE,
                }
            )
    diagnostic = pd.DataFrame(diagnostic_rows)
    # Annotations absent from every unit have undefined standardized
    # contrasts and contain no information about the fitted clusters.
    diagnostic = diagnostic.loc[
        np.isfinite(diagnostic.standardized_difference)
    ].copy()
    diagnostic["direction"] = np.where(
        diagnostic.standardized_difference.ge(0), "enriched", "depleted"
    )
    diagnostic["rank_within_direction"] = (
        diagnostic.assign(
            magnitude=np.where(
                diagnostic.direction.eq("enriched"),
                diagnostic.standardized_difference,
                -diagnostic.standardized_difference,
            )
        )
        .groupby(["cluster", "direction"]).magnitude.rank(
            method="first", ascending=False
        )
        .astype(int)
    )
    diagnostic = diagnostic.sort_values(
        ["cluster", "direction", "rank_within_direction", "feature"]
    )
    diagnostic.to_csv(
        TABLES / "core_cluster_all_diagnostic_feature_contrasts.tsv",
        sep="\t",
        index=False,
    )
    diagnostic.loc[diagnostic.rank_within_direction.le(5)].to_csv(
        TABLES / "core_cluster_top_diagnostic_feature_signatures.tsv",
        sep="\t",
        index=False,
    )

    # ACRO1 allele counts are summed before calculating fractions so units
    # containing more classified monomers contribute the appropriate weight.
    acro = clusters[["acre_id", "cluster"]].merge(
        acro1[["acre_id", "subfamily_1", "subfamily_2", "subfamily_3"]],
        on="acre_id",
        how="left",
        validate="one_to_one",
    )
    if acro[["subfamily_1", "subfamily_2", "subfamily_3"]].isna().any().any():
        raise ValueError("Some cluster units lack an ACRO1 composition row")
    acro_counts = acro.groupby("cluster")[[
        "subfamily_1",
        "subfamily_2",
        "subfamily_3",
    ]].sum()
    acro_fractions = acro_counts.div(acro_counts.sum(axis=1), axis=0)

    # A compact one-row-per-cluster table for direct manuscript use.
    position_counts = pd.crosstab(clusters.position_class, clusters.cluster).reindex(
        index=["first", "interior", "last"], columns=range(1, 9), fill_value=0
    )
    overall_cluster_fraction = clusters.cluster.value_counts(normalize=True)
    compact_rows: list[dict[str, object]] = []
    for cluster in range(1, 9):
        group = clusters.loc[clusters.cluster.eq(cluster)]
        top = component.loc[component.cluster.eq(cluster)].head(3)
        enriched = diagnostic.loc[
            diagnostic.cluster.eq(cluster) & diagnostic.direction.eq("enriched")
        ].head(3)
        depleted = diagnostic.loc[
            diagnostic.cluster.eq(cluster) & diagnostic.direction.eq("depleted")
        ].head(3)
        first_share = position_counts.loc["first", cluster] / position_counts.loc[
            "first"
        ].sum()
        last_share = position_counts.loc["last", cluster] / position_counts.loc[
            "last"
        ].sum()
        overall_share = float(overall_cluster_fraction[cluster])
        compact_rows.append(
            {
                "cluster_label": f"C{cluster}",
                "cluster": cluster,
                "complete_core_units": len(group),
                "fraction_of_all_core_units": len(group) / len(clusters),
                "arrays_with_cluster": group.source_name.nunique(),
                "individuals_with_cluster": group.individual_id.nunique(),
                "mean_unit_length_bp": group.length.mean(),
                "median_unit_length_bp": group.length.median(),
                "first_position_cluster_share": first_share,
                "first_position_enrichment_vs_overall": first_share / overall_share,
                "last_position_cluster_share": last_share,
                "last_position_enrichment_vs_overall": last_share / overall_share,
                "top_component_1": top.iloc[0].top_component,
                "top_component_1_fraction": top.iloc[0].within_cluster_fraction,
                "top_component_2": top.iloc[1].top_component if len(top) > 1 else "",
                "top_component_2_fraction": (
                    top.iloc[1].within_cluster_fraction if len(top) > 1 else 0.0
                ),
                "top_component_3": top.iloc[2].top_component if len(top) > 2 else "",
                "top_component_3_fraction": (
                    top.iloc[2].within_cluster_fraction if len(top) > 2 else 0.0
                ),
                "ACRO1_S1_canonical_fraction": acro_fractions.loc[
                    cluster, "subfamily_1"
                ],
                "ACRO1_S2_internal_CCC_deletion_fraction": acro_fractions.loc[
                    cluster, "subfamily_2"
                ],
                "ACRO1_S3_terminal_deletion_fraction": acro_fractions.loc[
                    cluster, "subfamily_3"
                ],
                "top_3_annotated_enrichments": "; ".join(
                    f"{row.feature} (d={row.standardized_difference:.3f})"
                    for row in enriched.itertuples()
                ),
                "top_3_annotated_depletions": "; ".join(
                    f"{row.feature} (d={row.standardized_difference:.3f})"
                    for row in depleted.itertuples()
                ),
                "interpretation_scope": SCOPE_NOTE,
            }
        )
    compact = pd.DataFrame(compact_rows)
    compact.to_csv(
        TABLES / "core_cluster_biological_summary.tsv", sep="\t", index=False
    )

    # Separate position table makes the end-polarity caveat transparent.
    position = position_counts.T.reset_index(names="cluster")
    position.insert(0, "cluster_label", "C" + position.cluster.astype(str))
    for label in ["first", "interior", "last"]:
        position[f"cluster_share_among_{label}_units"] = (
            position[label] / position[label].sum()
        )
    position["fraction_of_all_core_units"] = position.cluster.map(
        overall_cluster_fraction
    )
    position["orientation_note"] = (
        "first/last refer to input FASTA order; telomere-versus-centromere "
        "orientation has not been validated"
    )
    position.to_csv(TABLES / "core_cluster_position_summary.tsv", sep="\t", index=False)

    print(f"Wrote {len(compact)} cluster summaries and {len(diagnostic)} feature contrasts")


if __name__ == "__main__":
    main()
