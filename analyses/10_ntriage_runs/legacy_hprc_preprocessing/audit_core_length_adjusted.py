#!/usr/bin/env python3
"""Independent audit of core-only N-TRIAGE composition associations.

All generated tables are explicitly named as core-only.  The script does not
read the earlier full-flank tables and never overwrites the main analysis
tables.
"""

from __future__ import annotations

import itertools
import math
import os
import zlib
from pathlib import Path

import numpy as np
import pandas as pd


SEED = 20261126
PRIMARY_PERMUTATIONS = 99_999
ANCESTRY_PERMUTATIONS = 19_999
LENGTH_ADJUSTED_PERMUTATIONS = 99_999
BATCH_SIZE = 250

ROOT = Path(os.environ.get("ACRO1_LEGACY_HPRC_ROOT", Path(__file__).resolve().parents[1]))
ANALYSIS = ROOT / "manuscript_output" / "telomeric_ntriage_2794"
PREFIX = ANALYSIS / "native" / "telomeric_2794_core"
TABLES = ANALYSIS / "tables"
MANIFEST = ROOT / "manifest.validated.tsv"
CLUSTERS = Path(f"{PREFIX}.clusters.tsv")
ARRAYS = Path(f"{PREFIX}.arrays.tsv")
CLUSTER_NAMES = [f"C{i}" for i in range(1, 9)]


def rng_for(label: str) -> np.random.Generator:
    return np.random.default_rng(
        np.random.SeedSequence([SEED, zlib.crc32(label.encode("utf-8"))])
    )


def bh_adjust(values: np.ndarray) -> np.ndarray:
    p = np.asarray(values, dtype=float)
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    result = np.empty_like(ranked)
    result[order] = np.minimum(ranked, 1.0)
    return result


def signflip(values: np.ndarray, label: str, vector: bool) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    observed_mean = values.mean(axis=0)
    observed = (
        float(np.sum(observed_mean**2)) if vector else float(observed_mean.squeeze())
    )
    generator = rng_for(label)
    extreme = 0
    done = 0
    while done < PRIMARY_PERMUTATIONS:
        batch = min(2000, PRIMARY_PERMUTATIONS - done)
        signs = generator.choice(np.array([-1.0, 1.0]), size=(batch, len(values)))
        null_mean = signs @ values / len(values)
        if vector:
            null = np.sum(null_mean**2, axis=1)
            extreme += int(np.sum(null >= observed - 1e-15))
        else:
            extreme += int(np.sum(np.abs(null_mean) >= abs(observed) - 1e-15))
        done += batch
    return observed, (extreme + 1) / (PRIMARY_PERMUTATIONS + 1)


def qr_basis(design: np.ndarray) -> np.ndarray:
    q, r = np.linalg.qr(design)
    return q[:, : np.linalg.matrix_rank(r)]


def residual_sse_by_response(y: np.ndarray, q: np.ndarray) -> np.ndarray:
    """Return SSE by response for y shaped batch x observations x responses."""
    total = np.sum(y * y, axis=1)
    projection = np.einsum("np,bnd->bpd", q, y, optimize=True)
    return np.maximum(0.0, total - np.sum(projection * projection, axis=1))


def freedman_lane(
    y: np.ndarray,
    x0: np.ndarray,
    x1: np.ndarray,
    label: str,
    draws: int,
) -> tuple[float, float, float, int, int]:
    y = np.asarray(y, dtype=float)
    if y.ndim == 1:
        y = y[:, None]
    q0, q1 = qr_basis(x0), qr_basis(x1)
    fitted0 = q0 @ (q0.T @ y)
    residual0 = y - fitted0
    sse0 = float(np.sum(residual0**2))
    residual1 = y - q1 @ (q1.T @ y)
    sse1 = float(np.sum(residual1**2))
    df1 = q1.shape[1] - q0.shape[1]
    df2 = len(y) - q1.shape[1]
    observed = ((sse0 - sse1) / df1) / (sse1 / df2)
    generator = rng_for(label)
    extreme = 0
    done = 0
    while done < draws:
        batch = min(BATCH_SIZE, draws - done)
        orders = np.vstack([generator.permutation(len(y)) for _ in range(batch)])
        yp = fitted0[None, :, :] + residual0[orders, :]
        null0 = residual_sse_by_response(yp, q0).sum(axis=1)
        null1 = residual_sse_by_response(yp, q1).sum(axis=1)
        null_f = ((null0 - null1) / df1) / (null1 / df2)
        extreme += int(np.sum(null_f >= observed - 1e-15))
        done += batch
    return observed, (sse0 - sse1) / sse0, (extreme + 1) / (draws + 1), df1, df2


def matched_person_differences(
    frame: pd.DataFrame, value_columns: list[str]
) -> tuple[pd.DataFrame, int]:
    rows: list[dict[str, object]] = []
    for (individual, haplotype), group in frame.groupby(
        ["individual_id", "haplotype"], sort=False
    ):
        if set(group["has_recombination"].astype(int)) != {0, 1}:
            continue
        rec = group.loc[group.has_recombination == 1, value_columns].mean()
        non = group.loc[group.has_recombination == 0, value_columns].mean()
        row: dict[str, object] = {"individual_id": individual, "haplotype": haplotype}
        row.update({column: float(rec[column] - non[column]) for column in value_columns})
        rows.append(row)
    strata = pd.DataFrame(rows)
    people = strata.groupby("individual_id", as_index=False)[value_columns].mean()
    return people, len(strata)


def block_center(values: np.ndarray, block_index: np.ndarray, sizes: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    result = values.copy()
    if values.ndim == 1:
        means = np.bincount(block_index, weights=values) / sizes
        return values - means[block_index]
    for column in range(values.shape[1]):
        means = np.bincount(block_index, weights=values[:, column]) / sizes
        result[:, column] -= means[block_index]
    return result


def main() -> None:
    clusters = pd.read_csv(CLUSTERS, sep="\t", usecols=["cluster", "source_name"])
    arrays = pd.read_csv(ARRAYS, sep="\t", usecols=["source_name", "input_length"])
    manifest = pd.read_csv(
        MANIFEST,
        sep="\t",
        dtype={"haplotype": str, "has_recombination": int},
        keep_default_na=False,
    )

    cluster_counts = (
        pd.crosstab(clusters.source_name, clusters.cluster)
        .reindex(columns=range(1, 9), fill_value=0)
        .rename(columns={i: f"C{i}" for i in range(1, 9)})
    )
    frame = (
        manifest.set_index("array_id")
        .join(arrays.set_index("source_name"), how="left")
        .join(cluster_counts, how="left")
        .reset_index()
    )
    if len(frame) != 2794 or frame[CLUSTER_NAMES + ["input_length"]].isna().any().any():
        raise ValueError("Core arrays, clusters, and manifest did not join one-to-one")

    counts = frame[CLUSTER_NAMES].to_numpy(dtype=float)
    unit_totals = counts.sum(axis=1)
    fractions = counts / unit_totals[:, None]
    pseudocount = (counts + 0.5) / (unit_totals[:, None] + 4.0)
    clr = np.log(pseudocount)
    clr -= clr.mean(axis=1, keepdims=True)
    hellinger = np.sqrt(fractions)
    transforms = {"CLR_count_plus_0.5": clr, "Hellinger": hellinger, "raw_fraction": fractions}
    for index, name in enumerate(CLUSTER_NAMES):
        frame[f"fraction_{name}"] = fractions[:, index]
        frame[f"clr_{name}"] = clr[:, index]

    # Core unit and cluster counts.
    count_table = pd.DataFrame(
        {
            "cluster": range(1, 9),
            "complete_core_units": [int((clusters.cluster == i).sum()) for i in range(1, 9)],
        }
    )
    count_table["fraction_of_core_units"] = count_table.complete_core_units / len(clusters)
    count_table["total_complete_core_units"] = len(clusters)
    count_table["arrays"] = len(frame)
    count_table["source_prefix"] = str(PREFIX)
    count_table.to_csv(TABLES / "audit_native_cluster_counts_core.tsv", sep="\t", index=False)

    # Primary individual-level matched recombination audit.
    clr_columns = [f"clr_{name}" for name in CLUSTER_NAMES]
    fraction_columns = [f"fraction_{name}" for name in CLUSTER_NAMES]
    person_clr, informative_strata = matched_person_differences(frame, clr_columns)
    person_fraction, _ = matched_person_differences(frame, fraction_columns)
    global_stat, global_p = signflip(
        person_clr[clr_columns].to_numpy(), "core_primary_recomb_global", vector=True
    )
    recombination_rows = [
        {
            "analysis": "global_CLR_Aitchison",
            "cluster": "all_8",
            "individuals": len(person_clr),
            "informative_individual_haplotype_strata": informative_strata,
            "effect": global_stat,
            "effect_definition": "squared_norm_of_mean_matched_CLR_difference",
            "p_value": global_p,
            "q_bh_across_cluster_rows": np.nan,
            "permutations": PRIMARY_PERMUTATIONS,
        }
    ]
    cluster_p_values: list[float] = []
    cluster_rows: list[dict[str, object]] = []
    for cluster, column in enumerate(fraction_columns, start=1):
        effect, p_value = signflip(
            person_fraction[[column]].to_numpy(), f"core_primary_recomb_C{cluster}", vector=False
        )
        cluster_p_values.append(p_value)
        cluster_rows.append(
            {
                "analysis": "cluster_fraction",
                "cluster": cluster,
                "individuals": len(person_fraction),
                "informative_individual_haplotype_strata": informative_strata,
                "effect": effect,
                "effect_definition": "mean_within_person_recombination_minus_nonrecombination_fraction",
                "p_value": p_value,
                "permutations": PRIMARY_PERMUTATIONS,
            }
        )
    q_values = bh_adjust(np.asarray(cluster_p_values))
    for row, q_value in zip(cluster_rows, q_values):
        row["q_bh_across_cluster_rows"] = q_value
    recombination_rows.extend(cluster_rows)
    pd.DataFrame(recombination_rows).to_csv(
        TABLES / "audit_recombination_primary_core.tsv", sep="\t", index=False
    )

    # Primary superpopulation audit: one independent row per person.
    individual_rows: list[dict[str, object]] = []
    for individual, group in frame.groupby("individual_id", sort=False):
        row: dict[str, object] = {
            "individual_id": individual,
            "superpopulation": group.superpopulation.iloc[0],
            "recombination_fraction": float(group.has_recombination.mean()),
            "mean_log_core_array_length": float(np.log(group.input_length).mean()),
        }
        for column in fraction_columns + clr_columns:
            row[column] = float(group[column].mean())
        individual_rows.append(row)
    individual = pd.DataFrame(individual_rows)
    individual = individual.loc[individual.superpopulation != "NA"].copy()
    retained = individual.groupby("superpopulation").size()
    retained = retained[retained >= 5].index
    individual = individual.loc[individual.superpopulation.isin(retained)].copy()
    levels = sorted(individual.superpopulation.unique())
    dummies = pd.get_dummies(individual.superpopulation, dtype=float).reindex(columns=levels[1:])
    x0 = np.column_stack(
        [
            np.ones(len(individual)),
            individual.recombination_fraction,
            individual.mean_log_core_array_length,
        ]
    )
    x1 = np.column_stack([x0, dummies.to_numpy()])
    ancestry_rows: list[dict[str, object]] = []
    result = freedman_lane(
        individual[clr_columns].to_numpy(),
        x0,
        x1,
        "core_primary_superpopulation_global",
        ANCESTRY_PERMUTATIONS,
    )
    ancestry_rows.append(
        {
            "analysis": "global_CLR_Aitchison",
            "cluster": "all_8",
            "pseudo_F": result[0],
            "partial_R2": result[1],
            "p_value": result[2],
            "q_bh_across_cluster_rows": np.nan,
            "df_term": result[3],
            "df_residual": result[4],
            "individuals": len(individual),
            "groups": len(levels),
            "levels": ",".join(levels),
            "covariates": "individual_recombination_fraction;mean_log_core_array_length",
            "permutations": ANCESTRY_PERMUTATIONS,
        }
    )
    ancestry_cluster_rows: list[dict[str, object]] = []
    ancestry_p: list[float] = []
    for cluster, column in enumerate(fraction_columns, start=1):
        result = freedman_lane(
            individual[column].to_numpy(),
            x0,
            x1,
            f"core_primary_superpopulation_C{cluster}",
            ANCESTRY_PERMUTATIONS,
        )
        ancestry_p.append(result[2])
        ancestry_cluster_rows.append(
            {
                "analysis": "cluster_fraction",
                "cluster": cluster,
                "pseudo_F": result[0],
                "partial_R2": result[1],
                "p_value": result[2],
                "df_term": result[3],
                "df_residual": result[4],
                "individuals": len(individual),
                "groups": len(levels),
                "levels": ",".join(levels),
                "covariates": "individual_recombination_fraction;mean_log_core_array_length",
                "permutations": ANCESTRY_PERMUTATIONS,
            }
        )
    for row, q_value in zip(ancestry_cluster_rows, bh_adjust(np.asarray(ancestry_p))):
        row["q_bh_across_cluster_rows"] = q_value
    ancestry_rows.extend(ancestry_cluster_rows)
    pd.DataFrame(ancestry_rows).to_csv(
        TABLES / "audit_superpopulation_primary_core.tsv", sep="\t", index=False
    )

    # Core-only, length-adjusted, individual×haplotype-blocked sensitivity.
    block_labels = (frame.individual_id.astype(str) + "#" + frame.haplotype.astype(str)).to_numpy()
    _, block_index = np.unique(block_labels, return_inverse=True)
    block_count = int(block_index.max() + 1)
    block_sizes = np.bincount(block_index)
    block_members = [np.flatnonzero(block_index == i) for i in range(block_count)]
    log_length = block_center(np.log(frame.input_length.to_numpy(float)), block_index, block_sizes)
    recombination = block_center(
        frame.has_recombination.to_numpy(float), block_index, block_sizes
    )
    design0 = log_length[:, None]
    design1 = np.column_stack([log_length, recombination])
    q0, q1 = qr_basis(design0), qr_basis(design1)
    df_term = q1.shape[1] - q0.shape[1]
    df_residual = len(frame) - block_count - q1.shape[1]

    observed: dict[str, dict[str, object]] = {}
    person_effects: dict[str, np.ndarray] = {}
    for name, values in transforms.items():
        centered = block_center(values, block_index, block_sizes)
        fitted0 = q0 @ (q0.T @ centered)
        residual0 = centered - fitted0
        residual1 = centered - q1 @ (q1.T @ centered)
        sse0_by_cluster = np.sum(residual0**2, axis=0)
        sse1_by_cluster = np.sum(residual1**2, axis=0)
        sse0, sse1 = float(sse0_by_cluster.sum()), float(sse1_by_cluster.sum())
        coefficient = np.linalg.lstsq(design1, centered, rcond=None)[0][1]
        global_f = ((sse0 - sse1) / df_term) / (sse1 / df_residual)
        cluster_f = ((sse0_by_cluster - sse1_by_cluster) / df_term) / (
            sse1_by_cluster / df_residual
        )
        observed[name] = {
            "fitted0": fitted0,
            "residual0": residual0,
            "global_f": global_f,
            "global_r2": (sse0 - sse1) / sse0,
            "cluster_f": cluster_f,
            "cluster_r2": (sse0_by_cluster - sse1_by_cluster) / sse0_by_cluster,
            "coefficient": coefficient,
            "global_extreme": 0,
            "cluster_extreme": np.zeros(8, dtype=int),
        }
        temp = frame[["individual_id", "haplotype", "has_recombination"]].copy()
        value_columns = [f"value_{i}" for i in range(8)]
        temp[value_columns] = values
        matched, _ = matched_person_differences(temp, value_columns)
        person_effects[name] = matched[value_columns].mean().to_numpy()

    generator = rng_for("core_length_adjusted_blocked")
    done = 0
    while done < LENGTH_ADJUSTED_PERMUTATIONS:
        batch = min(BATCH_SIZE, LENGTH_ADJUSTED_PERMUTATIONS - done)
        orders = np.empty((batch, len(frame)), dtype=np.int32)
        for members in block_members:
            local_order = np.argsort(generator.random((batch, len(members))), axis=1)
            orders[:, members] = members[local_order]
        for name in transforms:
            item = observed[name]
            permuted = item["fitted0"][None, :, :] + item["residual0"][orders, :]
            null0_by_cluster = residual_sse_by_response(permuted, q0)
            null1_by_cluster = residual_sse_by_response(permuted, q1)
            null0, null1 = null0_by_cluster.sum(axis=1), null1_by_cluster.sum(axis=1)
            null_global_f = ((null0 - null1) / df_term) / (null1 / df_residual)
            null_cluster_f = ((null0_by_cluster - null1_by_cluster) / df_term) / (
                null1_by_cluster / df_residual
            )
            item["global_extreme"] += int(
                np.sum(null_global_f >= item["global_f"] - 1e-15)
            )
            item["cluster_extreme"] += np.sum(
                null_cluster_f >= item["cluster_f"][None, :] - 1e-15, axis=0
            )
        done += batch

    global_rows: list[dict[str, object]] = []
    sensitivity_rows: list[dict[str, object]] = []
    for name, item in observed.items():
        global_rows.append(
            {
                "transformation": name,
                "arrays": len(frame),
                "complete_core_units": len(clusters),
                "individual_haplotype_blocks": block_count,
                "informative_individual_haplotype_blocks": informative_strata,
                "length_covariate": "log(core_input_length)",
                "core_length_source": f"{ARRAYS}:input_length",
                "pseudo_F": item["global_f"],
                "partial_R2": item["global_r2"],
                "df_term": df_term,
                "df_residual": df_residual,
                "permutations_within_individual_haplotype": LENGTH_ADJUSTED_PERMUTATIONS,
                "p_value": (item["global_extreme"] + 1)
                / (LENGTH_ADJUSTED_PERMUTATIONS + 1),
            }
        )
        cluster_p = (item["cluster_extreme"] + 1) / (
            LENGTH_ADJUSTED_PERMUTATIONS + 1
        )
        cluster_q = bh_adjust(cluster_p)
        for index, cluster in enumerate(range(1, 9)):
            unadjusted = float(person_effects[name][index])
            adjusted = float(item["coefficient"][index])
            sensitivity_rows.append(
                {
                    "transformation": name,
                    "cluster": cluster,
                    "unadjusted_person_matched_effect": unadjusted,
                    "length_adjusted_block_fixed_effect_beta": adjusted,
                    "effect_direction_stable": bool(
                        math.copysign(1.0, unadjusted) == math.copysign(1.0, adjusted)
                    ),
                    "pseudo_F": item["cluster_f"][index],
                    "partial_R2": item["cluster_r2"][index],
                    "p_value": cluster_p[index],
                    "q_bh_across_8_clusters": cluster_q[index],
                    "df_term": df_term,
                    "df_residual": df_residual,
                    "permutations_within_individual_haplotype": LENGTH_ADJUSTED_PERMUTATIONS,
                }
            )
    pd.DataFrame(global_rows).to_csv(
        TABLES / "recombination_length_adjusted_global_core.tsv", sep="\t", index=False
    )
    pd.DataFrame(sensitivity_rows).to_csv(
        TABLES / "recombination_length_adjusted_cluster_effects_core.tsv",
        sep="\t",
        index=False,
    )

    print(f"Audited {len(clusters):,} complete core units in {len(frame):,} arrays")
    print(f"Wrote five core-only audit tables to {TABLES}")


if __name__ == "__main__":
    main()
