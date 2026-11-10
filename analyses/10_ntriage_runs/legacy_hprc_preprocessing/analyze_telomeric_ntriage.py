#!/usr/bin/env python3
"""Array- and individual-aware analysis of the 2,794-array N-TRIAGE run."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import zlib
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import entropy
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score, silhouette_score
from sklearn.preprocessing import StandardScaler


SEED = 20261026
N_SIGNFLIP = 99_999
N_BOOTSTRAP = 20_000
N_POP_PERM = 19_999


def rng_for(label: str) -> np.random.Generator:
    return np.random.default_rng(np.random.SeedSequence([SEED, zlib.crc32(label.encode())]))


def bh_adjust(p_values: np.ndarray) -> np.ndarray:
    p = np.asarray(p_values, dtype=float)
    result = np.full_like(p, np.nan)
    valid = np.isfinite(p)
    pv = p[valid]
    if not len(pv):
        return result
    order = np.argsort(pv)
    ranked = pv[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.minimum(adjusted, 1.0)
    restored = np.empty_like(adjusted)
    restored[order] = adjusted
    result[valid] = restored
    return result


def signflip_scalar(values: np.ndarray, label: str, draws: int = N_SIGNFLIP) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    observed = float(np.mean(values))
    generator = rng_for(label)
    extreme = 0
    done = 0
    while done < draws:
        batch = min(2000, draws - done)
        signs = generator.choice(np.array([-1.0, 1.0]), size=(batch, len(values)))
        null = np.mean(signs * values[None, :], axis=1)
        extreme += int(np.sum(np.abs(null) >= abs(observed) - 1e-15))
        done += batch
    return observed, (extreme + 1) / (draws + 1)


def signflip_vector(values: np.ndarray, label: str, draws: int = N_SIGNFLIP) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    observed_vector = np.mean(values, axis=0)
    observed = float(np.sum(observed_vector ** 2))
    generator = rng_for(label)
    extreme = 0
    done = 0
    while done < draws:
        batch = min(2000, draws - done)
        signs = generator.choice(np.array([-1.0, 1.0]), size=(batch, len(values)))
        null_mean = signs @ values / len(values)
        null = np.sum(null_mean ** 2, axis=1)
        extreme += int(np.sum(null >= observed - 1e-15))
        done += batch
    return observed, (extreme + 1) / (draws + 1)


def bootstrap_mean_ci(values: np.ndarray, label: str, draws: int = N_BOOTSTRAP) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    generator = rng_for(label)
    means = np.empty(draws, dtype=float)
    done = 0
    while done < draws:
        batch = min(1000, draws - done)
        idx = generator.integers(0, len(values), size=(batch, len(values)))
        means[done : done + batch] = np.mean(values[idx], axis=1)
        done += batch
    return tuple(np.quantile(means, [0.025, 0.975]))


def matched_person_differences(frame: pd.DataFrame, value_cols: list[str]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for (individual, haplotype), group in frame.groupby(["individual_id", "haplotype"], sort=False):
        statuses = set(group["has_recombination"].astype(int))
        if statuses != {0, 1}:
            continue
        rec = group.loc[group["has_recombination"] == 1, value_cols].mean(axis=0)
        non = group.loc[group["has_recombination"] == 0, value_cols].mean(axis=0)
        row: dict[str, object] = {"individual_id": individual, "haplotype": haplotype}
        row.update({column: float(rec[column] - non[column]) for column in value_cols})
        rows.append(row)
    stratum = pd.DataFrame(rows)
    person = stratum.groupby("individual_id", as_index=False)[value_cols].mean()
    person["informative_haplotypes"] = stratum.groupby("individual_id").size().reindex(person["individual_id"]).to_numpy()
    return person


def haplotype_person_differences(frame: pd.DataFrame, value_cols: list[str]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for individual, group in frame.groupby("individual_id", sort=False):
        haps = set(group["haplotype"].astype(int))
        if haps != {1, 2}:
            continue
        h1 = group.loc[group["haplotype"] == 1, value_cols].mean(axis=0)
        h2 = group.loc[group["haplotype"] == 2, value_cols].mean(axis=0)
        row: dict[str, object] = {"individual_id": individual}
        row.update({column: float(h1[column] - h2[column]) for column in value_cols})
        rows.append(row)
    return pd.DataFrame(rows)


def qr_basis(design: np.ndarray) -> np.ndarray:
    q, r = np.linalg.qr(design)
    rank = int(np.linalg.matrix_rank(r))
    return q[:, :rank]


def residual_sse_batch(y: np.ndarray, q: np.ndarray) -> np.ndarray:
    # y: batch x observations x responses
    total = np.sum(y * y, axis=(1, 2))
    projection = np.einsum("np,bnd->bpd", q, y, optimize=True)
    return np.maximum(0.0, total - np.sum(projection * projection, axis=(1, 2)))


def freedman_lane_test(
    y: np.ndarray,
    x_reduced: np.ndarray,
    x_full: np.ndarray,
    label: str,
    draws: int = N_POP_PERM,
) -> dict[str, float]:
    y = np.asarray(y, dtype=float)
    if y.ndim == 1:
        y = y[:, None]
    q0 = qr_basis(x_reduced)
    q1 = qr_basis(x_full)
    fitted0 = q0 @ (q0.T @ y)
    residual0 = y - fitted0
    sse0 = float(np.sum(residual0 ** 2))
    residual1 = y - q1 @ (q1.T @ y)
    sse1 = float(np.sum(residual1 ** 2))
    df1 = q1.shape[1] - q0.shape[1]
    df2 = len(y) - q1.shape[1]
    observed = ((sse0 - sse1) / df1) / (sse1 / df2) if sse1 > 0 else math.inf
    generator = rng_for(label)
    extreme = 0
    done = 0
    while done < draws:
        batch = min(250, draws - done)
        permutations = np.vstack([generator.permutation(len(y)) for _ in range(batch)])
        y_perm = fitted0[None, :, :] + residual0[permutations, :]
        null_sse0 = residual_sse_batch(y_perm, q0)
        null_sse1 = residual_sse_batch(y_perm, q1)
        null_f = ((null_sse0 - null_sse1) / df1) / (null_sse1 / df2)
        extreme += int(np.sum(null_f >= observed - 1e-15))
        done += batch
    return {
        "pseudo_F": observed,
        "p_value": (extreme + 1) / (draws + 1),
        "partial_R2": (sse0 - sse1) / sse0 if sse0 > 0 else np.nan,
        "df_term": df1,
        "df_residual": df2,
        "permutations": draws,
    }


def dataframe_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    args = parser.parse_args()
    root = args.root
    native = root / "native"
    inputs = root / "inputs"
    tables = root / "tables"
    tables.mkdir(parents=True, exist_ok=True)
    # Analyze only the manifest-defined ACRO cores. The source FASTAs include
    # approximately 10-kb extraction flanks that must not enter N-TRIAGE.
    prefix = native / "telomeric_2794_core"

    manifest = pd.read_csv(args.manifest, sep="\t", na_values=["NA"], keep_default_na=True)
    arrays = pd.read_csv(f"{prefix}.arrays.tsv", sep="\t")
    clusters = pd.read_csv(f"{prefix}.clusters.tsv", sep="\t")
    acro = pd.read_csv(f"{prefix}.acro1_composition.tsv", sep="\t")
    switches = pd.read_csv(f"{prefix}.switches.tsv", sep="\t")
    audit = pd.read_csv(inputs / "input_audit.tsv", sep="\t")

    if len(manifest) != 2794 or len(arrays) != 2794:
        raise ValueError("Unexpected manifest or array row count")
    expected_units = int(arrays["complete_acres"].sum())
    if len(clusters) != expected_units:
        raise ValueError(
            f"Cluster table has {len(clusters):,} rows but array summary reports "
            f"{expected_units:,} complete ACRO units"
        )
    if set(arrays["source_name"]) != set(manifest["array_id"]):
        raise ValueError("N-TRIAGE array names do not match the manifest")

    arrays = arrays.rename(columns={"array_id": "ntriage_array_id"})
    array_frame = arrays.merge(manifest, left_on="source_name", right_on="array_id", validate="one_to_one")
    array_frame = array_frame.merge(audit[["array_id", "sequence_sha256"]], on="array_id", validate="one_to_one")

    cluster_counts = pd.crosstab(clusters["array_id"], clusters["cluster"]).reindex(columns=range(1, 9), fill_value=0)
    cluster_counts.columns = [f"cluster_{column}_count" for column in cluster_counts.columns]
    cluster_counts.index.name = "ntriage_array_id"
    cluster_counts = cluster_counts.reset_index()
    array_frame = array_frame.merge(cluster_counts, on="ntriage_array_id", validate="one_to_one")

    count_cols = [f"cluster_{i}_count" for i in range(1, 9)]
    frac_cols = [f"cluster_{i}_fraction" for i in range(1, 9)]
    clr_cols = [f"cluster_{i}_clr" for i in range(1, 9)]
    counts = array_frame[count_cols].to_numpy(float)
    fractions = counts / counts.sum(axis=1, keepdims=True)
    pseudocount = counts + 0.5
    logs = np.log(pseudocount / pseudocount.sum(axis=1, keepdims=True))
    clr = logs - logs.mean(axis=1, keepdims=True)
    array_frame[frac_cols] = fractions
    array_frame[clr_cols] = clr
    array_frame["cluster_entropy"] = np.apply_along_axis(entropy, 1, fractions)
    array_frame["dominant_cluster"] = np.argmax(fractions, axis=1) + 1

    acro_array = acro.groupby("array_id", as_index=False)[
        ["subfamily_1", "subfamily_2", "subfamily_3", "classified_total"]
    ].sum()
    acro_array = acro_array.rename(columns={
        "subfamily_1": "acro1_s1_count", "subfamily_2": "acro1_s2_count",
        "subfamily_3": "acro1_s3_count", "classified_total": "acro1_classified_total",
    })
    array_frame = array_frame.merge(acro_array, left_on="ntriage_array_id", right_on="array_id",
                                    suffixes=("", "_acro"), validate="one_to_one")
    s_count_cols = ["acro1_s1_count", "acro1_s2_count", "acro1_s3_count"]
    s_frac_cols = ["acro1_s1_fraction", "acro1_s2_fraction", "acro1_s3_fraction"]
    s_clr_cols = ["acro1_s1_clr", "acro1_s2_clr", "acro1_s3_clr"]
    s_counts = array_frame[s_count_cols].to_numpy(float)
    s_total = s_counts.sum(axis=1)
    array_frame[s_frac_cols] = s_counts / s_total[:, None]
    s_logs = np.log((s_counts + 0.5) / (s_total[:, None] + 1.5))
    array_frame[s_clr_cols] = s_logs - s_logs.mean(axis=1, keepdims=True)

    switch_counts = switches.groupby("array_id").size().rename("switch_count")
    array_frame = array_frame.merge(switch_counts, left_on="ntriage_array_id", right_index=True, how="left")
    array_frame["switch_count"] = array_frame["switch_count"].fillna(0).astype(int)
    array_frame["switches_per_100_units"] = 100 * array_frame["switch_count"] / array_frame["complete_acres"]
    array_frame["log_array_length"] = np.log(array_frame["input_length"])
    array_frame.to_csv(tables / "array_level_summary.tsv.gz", sep="\t", index=False, compression="gzip")

    # Basic run and cluster summaries.
    pd.DataFrame({
        "cluster": range(1, 9),
        "units": clusters.groupby("cluster").size().reindex(range(1, 9), fill_value=0).to_numpy(),
    }).assign(fraction=lambda d: d["units"] / len(clusters)).to_csv(
        tables / "native_cluster_counts.tsv", sep="\t", index=False
    )
    group_sizes = manifest.groupby(["has_recombination", "superpopulation"], dropna=False).size().reset_index(name="arrays")
    group_sizes.to_csv(tables / "manifest_group_sizes.tsv", sep="\t", index=False)

    # Recombination: individual/haplotype-matched, then individual-aggregated.
    recomb_clr_diff = matched_person_differences(array_frame, clr_cols)
    recomb_frac_diff = matched_person_differences(array_frame, frac_cols)
    recomb_clr_diff.to_csv(tables / "recombination_person_cluster_clr_differences.tsv", sep="\t", index=False)
    global_stat, global_p = signflip_vector(recomb_clr_diff[clr_cols].to_numpy(), "recomb_global_cluster_clr")
    global_recomb = pd.DataFrame([{
        "analysis": "cluster_composition_CLR_Aitchison",
        "contrast": "recombination_minus_no_recombination",
        "individuals": len(recomb_clr_diff),
        "informative_individual_haplotype_strata": int(recomb_clr_diff["informative_haplotypes"].sum()),
        "squared_mean_CLR_distance": global_stat,
        "permutations": N_SIGNFLIP,
        "two_sided_vector_signflip_p": global_p,
    }])
    global_recomb.to_csv(tables / "recombination_global_composition_test.tsv", sep="\t", index=False)

    descriptive = array_frame.groupby("has_recombination")[frac_cols].mean().T
    descriptive.columns = ["mean_fraction_no_recombination", "mean_fraction_recombination"]
    cluster_tests: list[dict[str, object]] = []
    for i, column in enumerate(frac_cols, start=1):
        values = recomb_frac_diff[column].to_numpy()
        estimate, p_value = signflip_scalar(values, f"recomb_cluster_{i}")
        ci_lo, ci_hi = bootstrap_mean_ci(values, f"recomb_cluster_{i}_bootstrap")
        cluster_tests.append({
            "cluster": i,
            "individuals": len(values),
            "mean_fraction_no_recombination": descriptive.loc[column, "mean_fraction_no_recombination"],
            "mean_fraction_recombination": descriptive.loc[column, "mean_fraction_recombination"],
            "mean_within_person_fraction_difference": estimate,
            "bootstrap_ci95_lo": ci_lo,
            "bootstrap_ci95_hi": ci_hi,
            "two_sided_signflip_p": p_value,
        })
    cluster_tests_df = pd.DataFrame(cluster_tests)
    cluster_tests_df["q_bh_across_8_clusters"] = bh_adjust(cluster_tests_df["two_sided_signflip_p"].to_numpy())
    cluster_tests_df.to_csv(tables / "recombination_cluster_fraction_tests.tsv", sep="\t", index=False)

    # Duplicate-sequence sensitivity.
    dedup = array_frame.sort_values("task_index").drop_duplicates("sequence_sha256", keep="first")
    dedup_clr = matched_person_differences(dedup, clr_cols)
    dedup_stat, dedup_p = signflip_vector(dedup_clr[clr_cols].to_numpy(), "recomb_global_deduplicated")
    pd.DataFrame([{
        "arrays_retained": len(dedup), "unique_sequences": dedup["sequence_sha256"].nunique(),
        "individuals": len(dedup_clr), "squared_mean_CLR_distance": dedup_stat,
        "two_sided_vector_signflip_p": dedup_p,
    }]).to_csv(tables / "recombination_exact_sequence_dedup_sensitivity.tsv", sep="\t", index=False)

    # ACRO1 structural composition and other array metrics.
    acro_clr_diff = matched_person_differences(array_frame, s_clr_cols)
    acro_frac_diff = matched_person_differences(array_frame, s_frac_cols)
    acro_stat, acro_p = signflip_vector(acro_clr_diff[s_clr_cols].to_numpy(), "recomb_acro1_global")
    acro_descriptive = array_frame.groupby("has_recombination")[s_frac_cols].mean().T
    acro_tests: list[dict[str, object]] = []
    for label, column in zip(["S1", "S2", "S3"], s_frac_cols):
        values = acro_frac_diff[column].to_numpy()
        estimate, p_value = signflip_scalar(values, f"recomb_acro1_{label}")
        ci_lo, ci_hi = bootstrap_mean_ci(values, f"recomb_acro1_{label}_bootstrap")
        acro_tests.append({
            "subfamily": label, "individuals": len(values),
            "mean_fraction_no_recombination": acro_descriptive.loc[column, 0],
            "mean_fraction_recombination": acro_descriptive.loc[column, 1],
            "mean_within_person_fraction_difference": estimate,
            "bootstrap_ci95_lo": ci_lo, "bootstrap_ci95_hi": ci_hi,
            "two_sided_signflip_p": p_value,
        })
    acro_tests_df = pd.DataFrame(acro_tests)
    acro_tests_df["q_bh_across_3_subfamilies"] = bh_adjust(acro_tests_df["two_sided_signflip_p"].to_numpy())
    acro_tests_df.to_csv(tables / "recombination_acro1_subfamily_tests.tsv", sep="\t", index=False)
    pd.DataFrame([{
        "analysis": "ACRO1_S1_S2_S3_CLR_Aitchison", "individuals": len(acro_clr_diff),
        "squared_mean_CLR_distance": acro_stat, "permutations": N_SIGNFLIP,
        "two_sided_vector_signflip_p": acro_p,
    }]).to_csv(tables / "recombination_acro1_global_test.tsv", sep="\t", index=False)

    metric_cols = ["switches_per_100_units", "switch_count", "complete_acres", "input_length",
                   "estimated_period", "periodicity_score", "cluster_entropy"]
    metric_diff = matched_person_differences(array_frame, metric_cols)
    metric_tests: list[dict[str, object]] = []
    for column in metric_cols:
        values = metric_diff[column].to_numpy()
        estimate, p_value = signflip_scalar(values, f"recomb_metric_{column}")
        ci_lo, ci_hi = bootstrap_mean_ci(values, f"recomb_metric_{column}_bootstrap")
        metric_tests.append({
            "metric": column, "individuals": len(values),
            "mean_no_recombination": array_frame.loc[array_frame.has_recombination == 0, column].mean(),
            "mean_recombination": array_frame.loc[array_frame.has_recombination == 1, column].mean(),
            "mean_within_person_difference": estimate,
            "bootstrap_ci95_lo": ci_lo, "bootstrap_ci95_hi": ci_hi,
            "two_sided_signflip_p": p_value,
        })
    metric_tests_df = pd.DataFrame(metric_tests)
    metric_tests_df["q_bh_across_metrics"] = bh_adjust(metric_tests_df["two_sided_signflip_p"].to_numpy())
    metric_tests_df.to_csv(tables / "recombination_array_metric_tests.tsv", sep="\t", index=False)

    # Haplotype 1 versus 2 is a phase-label QC test, not parent-of-origin inference.
    hap_clr_diff = haplotype_person_differences(array_frame, clr_cols)
    hap_frac_diff = haplotype_person_differences(array_frame, frac_cols)
    hap_stat, hap_p = signflip_vector(hap_clr_diff[clr_cols].to_numpy(), "haplotype_global")
    hap_rows: list[dict[str, object]] = []
    for i, column in enumerate(frac_cols, start=1):
        estimate, p_value = signflip_scalar(hap_frac_diff[column].to_numpy(), f"hap_cluster_{i}")
        hap_rows.append({"cluster": i, "individuals": len(hap_frac_diff),
                         "mean_haplotype1_minus_haplotype2": estimate,
                         "two_sided_signflip_p": p_value})
    hap_df = pd.DataFrame(hap_rows)
    hap_df["q_bh_across_8_clusters"] = bh_adjust(hap_df["two_sided_signflip_p"].to_numpy())
    hap_df.to_csv(tables / "haplotype_phase_cluster_tests.tsv", sep="\t", index=False)
    pd.DataFrame([{
        "analysis": "haplotype1_minus_haplotype2_CLR_Aitchison_QC", "individuals": len(hap_clr_diff),
        "squared_mean_CLR_distance": hap_stat, "permutations": N_SIGNFLIP,
        "two_sided_vector_signflip_p": hap_p,
    }]).to_csv(tables / "haplotype_phase_global_test.tsv", sep="\t", index=False)

    # Superpopulation and named-population analyses use one observation per individual.
    individual_rows: list[dict[str, object]] = []
    for individual, group in array_frame.groupby("individual_id", sort=False):
        superpops = group["superpopulation"].dropna().unique()
        populations = group["population"].dropna().unique()
        row: dict[str, object] = {
            "individual_id": individual,
            "superpopulation": superpops[0] if len(superpops) == 1 else np.nan,
            "population": populations[0] if len(populations) == 1 else np.nan,
            "arrays": len(group),
            "recombination_fraction": group["has_recombination"].mean(),
            "mean_log_array_length": group["log_array_length"].mean(),
            "mean_units": group["complete_acres"].mean(),
        }
        for column in frac_cols + clr_cols + s_frac_cols:
            row[column] = group[column].mean()
        individual_rows.append(row)
    individual = pd.DataFrame(individual_rows)
    individual.to_csv(tables / "individual_level_summary.tsv", sep="\t", index=False)

    def population_analysis(group_column: str, minimum_people: int, label_prefix: str) -> None:
        d = individual.dropna(subset=[group_column]).copy()
        counts_people = d.groupby(group_column).size()
        retained = counts_people[counts_people >= minimum_people].index
        d = d[d[group_column].isin(retained)].copy()
        levels = sorted(d[group_column].unique())
        baseline = levels[0]
        dummies = pd.get_dummies(d[group_column], dtype=float).reindex(columns=levels[1:])
        x0 = np.column_stack([
            np.ones(len(d)), d["recombination_fraction"].to_numpy(), d["mean_log_array_length"].to_numpy()
        ])
        x1 = np.column_stack([x0, dummies.to_numpy()])
        global_result = freedman_lane_test(d[clr_cols].to_numpy(), x0, x1, f"{label_prefix}_global")
        global_result.update({
            "group_variable": group_column, "individuals": len(d), "groups": len(levels),
            "minimum_people_per_group": minimum_people, "reference_level": baseline,
            "levels": ",".join(levels),
        })
        pd.DataFrame([global_result]).to_csv(tables / f"{label_prefix}_global_composition_test.tsv",
                                             sep="\t", index=False)
        rows: list[dict[str, object]] = []
        for i, column in enumerate(frac_cols, start=1):
            result = freedman_lane_test(d[column].to_numpy(), x0, x1, f"{label_prefix}_cluster_{i}")
            result.update({"cluster": i, "individuals": len(d), "groups": len(levels)})
            rows.append(result)
        tests = pd.DataFrame(rows)
        tests["q_bh_across_8_clusters"] = bh_adjust(tests["p_value"].to_numpy())
        tests.to_csv(tables / f"{label_prefix}_cluster_tests.tsv", sep="\t", index=False)
        means = d.groupby(group_column)[frac_cols + s_frac_cols].mean().reset_index()
        n_people = d.groupby(group_column).size().rename("individuals").reset_index()
        means.merge(n_people, on=group_column).to_csv(tables / f"{label_prefix}_individual_means.tsv",
                                                       sep="\t", index=False)

    population_analysis("superpopulation", 5, "superpopulation")
    population_analysis("population", 5, "population_exploratory")

    # Array-balanced PCA and feature signatures.
    features_path = Path(f"{prefix}.features.tsv")
    features = pd.read_csv(features_path, sep="\t")
    if len(features) != len(clusters) or not np.array_equal(features["acre_id"], clusters["acre_id"]):
        raise ValueError("Feature and cluster rows are not aligned")
    feature_names = [c for c in features.columns if c not in {"acre_id", "array_id", "cluster"}]
    x = features[feature_names].to_numpy(dtype=np.float32)
    by_array = features.groupby("array_id", sort=False).indices
    generator = rng_for("array_balanced_pca_sample")
    pca_sample_idx = np.sort(np.concatenate([
        generator.choice(indices, size=min(10, len(indices)), replace=False)
        for indices in by_array.values()
    ]))
    variance = np.var(x[pca_sample_idx], axis=0)
    variable = variance > 0
    x = x[:, variable]
    feature_names = list(np.asarray(feature_names)[variable])
    scaler = StandardScaler().fit(x[pca_sample_idx])
    x_sample_scaled = scaler.transform(x[pca_sample_idx]).astype(np.float32, copy=False)
    pca = PCA(n_components=20, svd_solver="randomized", random_state=SEED)
    pca.fit(x_sample_scaled)
    scores = np.empty((len(x), 20), dtype=np.float32)
    for start in range(0, len(x), 10_000):
        stop = min(len(x), start + 10_000)
        scores[start:stop] = pca.transform(scaler.transform(x[start:stop]))
    variance_table = pd.DataFrame({
        "principal_component": np.arange(1, 21),
        "explained_variance_ratio": pca.explained_variance_ratio_,
        "cumulative_explained_variance": np.cumsum(pca.explained_variance_ratio_),
    })
    variance_table.to_csv(tables / "pca_explained_variance.tsv", sep="\t", index=False)

    unit_meta = clusters[["acre_id", "array_id", "unit_in_array", "cluster", "source_name"]].copy()
    unit_meta = unit_meta.merge(manifest[["array_id", "individual_id", "haplotype", "population",
                                          "superpopulation", "has_recombination"]],
                                left_on="source_name", right_on="array_id", suffixes=("", "_manifest"),
                                validate="many_to_one")
    for pc in range(5):
        unit_meta[f"PC{pc + 1}"] = scores[:, pc]
    unit_meta.to_csv(tables / "unit_pca_metadata.tsv.gz", sep="\t", index=False, compression="gzip")

    sample_labels = features.loc[pca_sample_idx, "cluster"].to_numpy(int)
    signature_rows: list[dict[str, object]] = []
    centroid_rows: list[dict[str, object]] = []
    for cluster_id in range(1, 9):
        own = x_sample_scaled[sample_labels == cluster_id].mean(axis=0)
        other = x_sample_scaled[sample_labels != cluster_id].mean(axis=0)
        delta = own - other
        top = np.argsort(np.abs(delta))[::-1][:20]
        for rank, idx in enumerate(top, start=1):
            signature_rows.append({
                "cluster": cluster_id, "rank": rank, "feature": feature_names[idx],
                "standardized_centroid_difference": float(delta[idx]),
            })
        centroid_rows.append({"cluster": cluster_id, **{
            feature_names[j]: float(own[j]) for j in range(len(feature_names))
        }})
    pd.DataFrame(signature_rows).to_csv(tables / "cluster_top_feature_signatures.tsv", sep="\t", index=False)
    pd.DataFrame(centroid_rows).to_csv(tables / "cluster_standardized_centroids.tsv.gz", sep="\t",
                                          index=False, compression="gzip")

    # k sensitivity on at most four units per array. Reproduce N-TRIAGE's
    # clustering geometry: sample-SD standardization across every unit followed
    # by L2 normalization of each unit vector (src/core.cpp::kmeans).
    generator = rng_for("k_sensitivity_sample")
    k_idx = np.sort(np.concatenate([
        generator.choice(indices, size=min(4, len(indices)), replace=False)
        for indices in by_array.values()
    ]))
    native_mean = np.mean(x, axis=0, dtype=np.float64)
    native_sd = np.std(x, axis=0, ddof=1, dtype=np.float64)
    native_sd[native_sd <= 1e-12] = 1.0
    x_k = ((x[k_idx].astype(np.float64) - native_mean) / native_sd).astype(np.float32)
    row_norm = np.linalg.norm(x_k, axis=1, keepdims=True)
    row_norm[row_norm <= 0] = 1.0
    x_k /= row_norm
    native_labels = features.loc[k_idx, "cluster"].to_numpy(int)
    k_rows: list[dict[str, object]] = []
    for k in range(2, 15):
        model = KMeans(n_clusters=k, n_init=5, max_iter=300, random_state=SEED + k)
        labels = model.fit_predict(x_k)
        silhouette = silhouette_score(x_k, labels, metric="euclidean", sample_size=min(5000, len(x_k)),
                                      random_state=SEED)
        k_rows.append({
            "k": k, "sampled_units": len(x_k), "arrays_represented": len(by_array),
            "silhouette": silhouette, "inertia": model.inertia_,
            "adjusted_rand_vs_native_k8": adjusted_rand_score(native_labels, labels) if k == 8 else np.nan,
        })
    native_silhouette = silhouette_score(x_k, native_labels, metric="euclidean",
                                         sample_size=min(5000, len(x_k)), random_state=SEED)
    k_table = pd.DataFrame(k_rows)
    k_table["native_k8_silhouette_on_same_sample"] = native_silhouette
    k_table.to_csv(tables / "k_sensitivity_array_balanced_sample.tsv", sep="\t", index=False)

    provenance = {
        "seed": SEED,
        "arrays": len(array_frame),
        "individuals": int(array_frame["individual_id"].nunique()),
        "complete_ACRO_units": len(clusters),
        "total_input_bp": int(array_frame["input_length"].sum()),
        "native_clusters": 8,
        "ACRO1_classified_monomers": int(array_frame["acro1_classified_total"].sum()),
        "deep_switch_events": len(switches),
        "features_used_for_PCA": len(feature_names),
        "array_balanced_PCA_sample_units": len(pca_sample_idx),
        "signflip_draws": N_SIGNFLIP,
        "bootstrap_draws": N_BOOTSTRAP,
        "population_Freedman_Lane_permutations": N_POP_PERM,
        "ntriage_features_sha256": dataframe_sha256(features_path),
    }
    (tables / "analysis_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(json.dumps(provenance, indent=2))


if __name__ == "__main__":
    main()
