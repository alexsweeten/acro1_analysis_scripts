#!/usr/bin/env python3
"""Reproducible switch/ACRO1 audit for the corrected core-only N-TRIAGE run.

This script deliberately writes separate ``core_switch_*`` tables and does not
alter the primary analysis tables.  Switch rates use the number of boundaries
that N-TRIAGE could evaluate as the exposure, rather than the number of units.
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
import statsmodels
from statsmodels.genmod.cov_struct import Exchangeable
from statsmodels.genmod.families import Binomial, Poisson
from statsmodels.genmod.generalized_estimating_equations import GEE
from statsmodels.formula.api import ols


SEED = 20261026
N_SIGNFLIP = 200_000
N_BOOTSTRAP = 20_000

ROOT = Path(os.environ.get("ACRO1_LEGACY_HPRC_ROOT", Path(__file__).resolve().parents[1]))
ANALYSIS = ROOT / "manuscript_output" / "telomeric_ntriage_2794"
PREFIX = ANALYSIS / "native" / "telomeric_2794_core"
TABLES = ANALYSIS / "tables"
MANIFEST = ROOT / "manifest.validated.tsv"

ARRAYS = Path(f"{PREFIX}.arrays.tsv")
CLUSTERS = Path(f"{PREFIX}.clusters.tsv")
ACRO1 = Path(f"{PREFIX}.acro1_composition.tsv")
SWITCHES = Path(f"{PREFIX}.switches.tsv")
BOUNDARIES = Path(f"{PREFIX}.switch_boundaries.tsv")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def as_bool(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)
    return series.astype(str).str.lower().eq("true")


def signflip_bootstrap(
    values: np.ndarray, rng: np.random.Generator
) -> tuple[float, float, float, float]:
    """Mean, two-sided Monte-Carlo sign-flip p, percentile bootstrap CI."""
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    observed = float(values.mean())
    extreme = 0
    done = 0
    while done < N_SIGNFLIP:
        batch = min(5_000, N_SIGNFLIP - done)
        signs = rng.choice(np.array([-1.0, 1.0]), size=(batch, len(values)))
        null_mean = np.mean(signs * values[None, :], axis=1)
        extreme += int(np.sum(np.abs(null_mean) >= abs(observed) - 1e-15))
        done += batch
    p_value = (extreme + 1) / (N_SIGNFLIP + 1)

    bootstrap_means = np.empty(N_BOOTSTRAP, dtype=float)
    done = 0
    while done < N_BOOTSTRAP:
        batch = min(2_000, N_BOOTSTRAP - done)
        indices = rng.integers(0, len(values), size=(batch, len(values)))
        bootstrap_means[done : done + batch] = values[indices].mean(axis=1)
        done += batch
    ci_lo, ci_hi = np.quantile(bootstrap_means, [0.025, 0.975])
    return observed, float(p_value), float(ci_lo), float(ci_hi)


def build_array_frame() -> tuple[
    pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame
]:
    manifest = pd.read_csv(MANIFEST, sep="\t", na_values=["NA"])
    arrays = pd.read_csv(ARRAYS, sep="\t")
    clusters = pd.read_csv(CLUSTERS, sep="\t")
    acro = pd.read_csv(ACRO1, sep="\t")
    switches = pd.read_csv(SWITCHES, sep="\t")
    boundaries = pd.read_csv(BOUNDARIES, sep="\t")
    boundaries["passes"] = as_bool(boundaries.passes_thresholds)

    acro_array = (
        acro.groupby("array_id")
        .agg(
            units=("acre_id", "size"),
            called_units=("classified_total", lambda values: int((values > 0).sum())),
            total=("classified_total", "sum"),
            s1=("subfamily_1", "sum"),
            s2=("subfamily_2", "sum"),
            s3=("subfamily_3", "sum"),
        )
        .reset_index()
    )
    for allele in (1, 2, 3):
        acro_array[f"fS{allele}"] = acro_array[f"s{allele}"] / acro_array.total
    acro_array["call_rate"] = acro_array.called_units / acro_array.units
    acro_array["calls_per_unit"] = acro_array.total / acro_array.units

    evaluable = boundaries.groupby("array_id").size().rename("n_eval").reset_index()
    raw_candidates = (
        boundaries.groupby("array_id").passes.sum().rename("raw_candidates").reset_index()
    )
    event_count = (
        switches.groupby("array_id").switch_id.size().rename("switch_count").reset_index()
    )

    data = (
        arrays.merge(acro_array, on="array_id", validate="one_to_one")
        .merge(evaluable, on="array_id", validate="one_to_one")
        .merge(raw_candidates, on="array_id", validate="one_to_one")
        .merge(event_count, on="array_id", how="left", validate="one_to_one")
        .merge(
            manifest,
            left_on="source_name",
            right_on="array_id",
            suffixes=("", "_manifest"),
            validate="one_to_one",
        )
    )
    data.switch_count = data.switch_count.fillna(0).astype(int)
    data["switch_rate100"] = 100 * data.switch_count / data.n_eval
    data["any_switch"] = data.switch_count.gt(0).astype(int)
    data["self_similarity"] = data.segmentation_method.eq("self_similarity").astype(int)
    data["period_z"] = (
        data.estimated_period - data.estimated_period.mean()
    ) / data.estimated_period.std(ddof=1)
    data["score_z"] = (
        data.periodicity_score - data.periodicity_score.mean()
    ) / data.periodicity_score.std(ddof=1)
    data["log_eval"] = np.log(data.n_eval)
    data["rate"] = data.switch_count / data.n_eval

    if len(data) != 2_794 or len(clusters) != 77_276:
        raise ValueError("Unexpected corrected-core array or unit count")
    if int(data.complete_acres.sum()) != len(clusters):
        raise ValueError("Array and unit counts disagree")
    expected_boundaries = int((data.complete_acres - 5).sum())
    if expected_boundaries != len(boundaries) or int(data.n_eval.sum()) != len(boundaries):
        raise ValueError("Evaluable-boundary count is not sum(complete units - 5)")
    if not np.array_equal(data.units.to_numpy(), data.complete_acres.to_numpy()):
        raise ValueError("ACRO1 and N-TRIAGE complete-unit counts disagree")

    return data, manifest, acro, switches, boundaries


def write_pooled_counts(
    data: pd.DataFrame,
    switches: pd.DataFrame,
    boundaries: pd.DataFrame,
) -> None:
    rows: list[dict[str, object]] = []
    groups: list[tuple[str, pd.DataFrame]] = [("all", data)]
    groups.extend(
        [
            ("no_recombination", data.loc[data.has_recombination.eq(0)]),
            ("recombination", data.loc[data.has_recombination.eq(1)]),
        ]
    )
    for label, group in groups:
        calls = float(group.total.sum())
        s1, s2, s3 = (int(group[f"s{i}"].sum()) for i in (1, 2, 3))
        rows.append(
            {
                "scope": label,
                "arrays": len(group),
                "individuals": group.individual_id.nunique(),
                "core_bp": int(group.input_length.sum()),
                "complete_core_units": int(group.complete_acres.sum()),
                "evaluable_boundaries": int(group.n_eval.sum()),
                "raw_passing_candidate_boundaries": int(group.raw_candidates.sum()),
                "collapsed_switch_events": int(group.switch_count.sum()),
                "arrays_with_at_least_one_switch": int(group.any_switch.sum()),
                "switch_events_per_100_evaluable_boundaries": (
                    100 * group.switch_count.sum() / group.n_eval.sum()
                ),
                "ACRO1_called_units": int(group.called_units.sum()),
                "ACRO1_unit_call_rate": group.called_units.sum() / group.units.sum(),
                "ACRO1_classified_monomers": int(calls),
                "ACRO1_S1_count": s1,
                "ACRO1_S2_count": s2,
                "ACRO1_S3_count": s3,
                "ACRO1_S1_fraction": s1 / calls,
                "ACRO1_S2_fraction": s2 / calls,
                "ACRO1_S3_fraction": s3 / calls,
                "self_similarity_arrays": int(group.self_similarity.sum()),
                "target_length_fallback_arrays": int((1 - group.self_similarity).sum()),
            }
        )
    pooled = pd.DataFrame(rows)
    pooled["candidate_collapse_definition"] = (
        "N-TRIAGE greedily retains passing boundaries at least 3 units apart"
    )
    pooled.to_csv(TABLES / "core_switch_acro1_pooled_counts.tsv", sep="\t", index=False)

    if pooled.loc[pooled.scope.eq("all"), "collapsed_switch_events"].iat[0] != len(
        switches
    ):
        raise ValueError("Pooled event count disagrees with switches.tsv")
    if pooled.loc[pooled.scope.eq("all"), "evaluable_boundaries"].iat[0] != len(
        boundaries
    ):
        raise ValueError("Pooled boundary count disagrees with switch_boundaries.tsv")


def write_matched_array_metrics(data: pd.DataFrame) -> None:
    # The fixed order and single RNG stream reproduce the original supplemental
    # audit exactly; changing the order changes Monte-Carlo draws.
    metric_order = [
        "input_length",
        "complete_acres",
        "estimated_period",
        "periodicity_score",
        "self_similarity",
        "call_rate",
        "calls_per_unit",
        "fS1",
        "fS2",
        "fS3",
        "switch_count",
        "switch_rate100",
        "any_switch",
    ]
    rng = np.random.default_rng(SEED)
    rows: list[dict[str, object]] = []
    for metric in metric_order:
        paired = (
            data.groupby(["individual_id", "has_recombination"])[metric]
            .mean()
            .unstack()
            .dropna()
        )
        values = (paired[1] - paired[0]).to_numpy()
        estimate, p_value, ci_lo, ci_hi = signflip_bootstrap(values, rng)
        rows.append(
            {
                "metric": metric,
                "matching_design": (
                    "direct donor match: mean across arrays within individual and status"
                ),
                "individuals": len(values),
                "mean_no_recombination_arrays": data.loc[
                    data.has_recombination.eq(0), metric
                ].mean(),
                "mean_recombination_arrays": data.loc[
                    data.has_recombination.eq(1), metric
                ].mean(),
                "mean_within_donor_recombination_minus_no_recombination": estimate,
                "bootstrap_ci95_lo": ci_lo,
                "bootstrap_ci95_hi": ci_hi,
                "two_sided_signflip_p": p_value,
                "signflip_draws": N_SIGNFLIP,
                "bootstrap_draws": N_BOOTSTRAP,
            }
        )
    pd.DataFrame(rows).to_csv(
        TABLES / "core_switch_acro1_donor_matched_metrics.tsv", sep="\t", index=False
    )


def write_models(data: pd.DataFrame) -> None:
    rhs = [
        ("unadjusted", "has_recombination"),
        (
            "QC_adjusted",
            "has_recombination+score_z+self_similarity+period_z",
        ),
        (
            "QC_and_ACRO1_adjusted",
            "has_recombination+score_z+self_similarity+period_z+fS1+fS2",
        ),
    ]
    rows: list[dict[str, object]] = []
    for label, predictors in rhs:
        formula = f"switch_count~{predictors}"
        fit = GEE.from_formula(
            formula,
            groups="individual_id",
            data=data,
            family=Poisson(),
            cov_struct=Exchangeable(),
            offset=data.log_eval,
        ).fit(maxiter=200)
        beta = float(fit.params["has_recombination"])
        se = float(fit.bse["has_recombination"])
        rows.append(
            {
                "model": f"Poisson_GEE_{label}",
                "outcome": "collapsed_switch_count",
                "effect_type": "rate_ratio_recombination_vs_no_recombination",
                "formula": formula,
                "exposure_or_offset": "offset=log(evaluable_boundaries)",
                "working_correlation": "Exchangeable within individual_id",
                "estimate": np.exp(beta),
                "ci95_lo": np.exp(beta - 1.96 * se),
                "ci95_hi": np.exp(beta + 1.96 * se),
                "coefficient_log_scale": beta,
                "robust_se_log_scale": se,
                "p_value": float(fit.pvalues["has_recombination"]),
                "arrays": len(data),
                "individuals": data.individual_id.nunique(),
            }
        )

        formula = f"any_switch~{predictors}"
        fit = GEE.from_formula(
            formula,
            groups="individual_id",
            data=data,
            family=Binomial(),
            cov_struct=Exchangeable(),
        ).fit(maxiter=200)
        beta = float(fit.params["has_recombination"])
        se = float(fit.bse["has_recombination"])
        rows.append(
            {
                "model": f"Binomial_GEE_{label}",
                "outcome": "any_collapsed_switch",
                "effect_type": "odds_ratio_recombination_vs_no_recombination",
                "formula": formula,
                "exposure_or_offset": "none",
                "working_correlation": "Exchangeable within individual_id",
                "estimate": np.exp(beta),
                "ci95_lo": np.exp(beta - 1.96 * se),
                "ci95_hi": np.exp(beta + 1.96 * se),
                "coefficient_log_scale": beta,
                "robust_se_log_scale": se,
                "p_value": float(fit.pvalues["has_recombination"]),
                "arrays": len(data),
                "individuals": data.individual_id.nunique(),
            }
        )

    formula = "rate~has_recombination+score_z+self_similarity+period_z+fS1+fS2"
    fit = ols(formula, data).fit(
        cov_type="cluster",
        cov_kwds={"groups": data.individual_id, "use_correction": True},
    )
    beta = float(fit.params["has_recombination"]) * 100
    se = float(fit.bse["has_recombination"]) * 100
    rows.append(
        {
            "model": "OLS_cluster_robust_QC_and_ACRO1_adjusted",
            "outcome": "switches_per_100_evaluable_boundaries",
            "effect_type": "additive_difference_recombination_minus_no_recombination",
            "formula": formula,
            "exposure_or_offset": "outcome already normalized by evaluable boundaries",
            "working_correlation": "cluster-robust SE by individual_id",
            "estimate": beta,
            "ci95_lo": beta - 1.96 * se,
            "ci95_hi": beta + 1.96 * se,
            "coefficient_log_scale": np.nan,
            "robust_se_log_scale": np.nan,
            "p_value": float(fit.pvalues["has_recombination"]),
            "arrays": len(data),
            "individuals": data.individual_id.nunique(),
        }
    )
    pd.DataFrame(rows).to_csv(
        TABLES / "core_switch_rate_models.tsv", sep="\t", index=False
    )


def build_switch_flank_table(
    data: pd.DataFrame,
    acro: pd.DataFrame,
    switches: pd.DataFrame,
    boundaries: pd.DataFrame,
) -> pd.DataFrame:
    composition = acro[
        [
            "array_id",
            "unit_in_array",
            "subfamily_1",
            "subfamily_2",
            "subfamily_3",
            "classified_total",
        ]
    ]
    rows: list[pd.DataFrame] = []
    key_columns = [
        "array_id",
        "boundary_after_unit",
        "boundary_before_unit",
        "boundary_position",
        "passes",
    ]
    for side, offsets, reference in [
        ("L", [-2, -1, 0], "boundary_after_unit"),
        ("R", [0, 1, 2], "boundary_before_unit"),
    ]:
        for offset in offsets:
            expanded = boundaries[key_columns].copy()
            expanded["side"] = side
            expanded["unit_in_array"] = expanded[reference] + offset
            rows.append(expanded)
    expanded = pd.concat(rows, ignore_index=True).merge(
        composition,
        on=["array_id", "unit_in_array"],
        how="left",
        validate="many_to_one",
    )
    if expanded.classified_total.isna().any():
        raise ValueError("A switch-flank unit did not join to ACRO1 composition")
    flank = (
        expanded.groupby(key_columns + ["side"])[
            ["subfamily_1", "subfamily_2", "subfamily_3", "classified_total"]
        ]
        .sum()
        .unstack("side")
    )
    flank.columns = [f"{name}_{side}" for name, side in flank.columns]
    flank = flank.reset_index()
    for allele in (1, 2, 3):
        flank[f"fS{allele}_L"] = (
            flank[f"subfamily_{allele}_L"] / flank.classified_total_L
        )
        flank[f"fS{allele}_R"] = (
            flank[f"subfamily_{allele}_R"] / flank.classified_total_R
        )
        flank[f"dS{allele}"] = flank[f"fS{allele}_R"] - flank[f"fS{allele}_L"]
    flank["TVD"] = 0.5 * sum(flank[f"dS{i}"].abs() for i in (1, 2, 3))
    flank["dominant_L"] = (
        flank[[f"fS{i}_L" for i in (1, 2, 3)]].to_numpy().argmax(axis=1) + 1
    )
    flank["dominant_R"] = (
        flank[[f"fS{i}_R" for i in (1, 2, 3)]].to_numpy().argmax(axis=1) + 1
    )
    flank["dominant_change"] = flank.dominant_L.ne(flank.dominant_R)

    primary_key = set(
        zip(
            switches.array_id,
            switches.boundary_after_unit,
            switches.boundary_before_unit,
        )
    )
    flank["primary"] = [
        (array, after, before) in primary_key
        for array, after, before in zip(
            flank.array_id,
            flank.boundary_after_unit,
            flank.boundary_before_unit,
        )
    ]
    metadata = data[
        ["array_id", "individual_id", "haplotype", "has_recombination"]
    ]
    flank = flank.merge(metadata, on="array_id", validate="many_to_one")
    primary = flank.loc[flank.primary & flank.TVD.notna()].copy()
    if len(primary) != len(switches):
        raise ValueError("Not every collapsed switch has an ACRO1 flank summary")
    return primary


def write_switch_flank_tests(primary: pd.DataFrame) -> None:
    rng = np.random.default_rng(SEED)
    metric_order = ["TVD", "abs_dS1", "abs_dS2", "abs_dS3"]
    for allele in (1, 2, 3):
        primary[f"abs_dS{allele}"] = primary[f"dS{allele}"].abs()
    rows: list[dict[str, object]] = []
    for metric in metric_order:
        paired = (
            primary.groupby(["individual_id", "has_recombination"])[metric]
            .mean()
            .unstack()
            .dropna()
        )
        values = (paired[1] - paired[0]).to_numpy()
        estimate, p_value, ci_lo, ci_hi = signflip_bootstrap(values, rng)
        rows.append(
            {
                "test": f"donor_matched_{metric}",
                "metric": metric,
                "comparison": "recombination_switches_minus_no_recombination_switches",
                "matching_design": "mean across collapsed switches within donor and status",
                "individuals": len(values),
                "events_no_recombination": int(
                    primary.has_recombination.eq(0).sum()
                ),
                "events_recombination": int(primary.has_recombination.eq(1).sum()),
                "event_mean_no_recombination": primary.loc[
                    primary.has_recombination.eq(0), metric
                ].mean(),
                "event_mean_recombination": primary.loc[
                    primary.has_recombination.eq(1), metric
                ].mean(),
                "effect": estimate,
                "bootstrap_ci95_lo": ci_lo,
                "bootstrap_ci95_hi": ci_hi,
                "two_sided_signflip_p": p_value,
                "signflip_draws": N_SIGNFLIP,
                "bootstrap_draws": N_BOOTSTRAP,
            }
        )

    # Direction tests ask whether left-to-right composition change among only
    # recombination-labelled switches differs from zero.  They are deliberately
    # separate from the absolute-change/TVD comparison above.
    for allele in (1, 2, 3):
        metric = f"dS{allele}"
        values = (
            primary.loc[primary.has_recombination.eq(1)]
            .groupby("individual_id")[metric]
            .mean()
            .to_numpy()
        )
        estimate, p_value, ci_lo, ci_hi = signflip_bootstrap(values, rng)
        rows.append(
            {
                "test": f"recombination_switch_signed_{metric}",
                "metric": metric,
                "comparison": "left_to_right_change_among_recombination_switches_vs_zero",
                "matching_design": "mean signed change across switches within donor",
                "individuals": len(values),
                "events_no_recombination": int(
                    primary.has_recombination.eq(0).sum()
                ),
                "events_recombination": int(primary.has_recombination.eq(1).sum()),
                "event_mean_no_recombination": primary.loc[
                    primary.has_recombination.eq(0), metric
                ].mean(),
                "event_mean_recombination": primary.loc[
                    primary.has_recombination.eq(1), metric
                ].mean(),
                "effect": estimate,
                "bootstrap_ci95_lo": ci_lo,
                "bootstrap_ci95_hi": ci_hi,
                "two_sided_signflip_p": p_value,
                "signflip_draws": N_SIGNFLIP,
                "bootstrap_draws": N_BOOTSTRAP,
            }
        )
    pd.DataFrame(rows).to_csv(
        TABLES / "core_switch_local_acro1_tests.tsv", sep="\t", index=False
    )

    dominant = (
        primary.groupby("has_recombination")
        .dominant_change.agg(["size", "sum", "mean"])
        .reset_index()
        .rename(
            columns={
                "size": "collapsed_switch_events",
                "sum": "dominant_ACRO1_allele_change_events",
                "mean": "dominant_ACRO1_allele_change_fraction",
            }
        )
    )
    dominant["scope"] = dominant.has_recombination.map(
        {0: "no_recombination", 1: "recombination"}
    )
    dominant["dominance_definition"] = (
        "argmax of pooled S1/S2/S3 monomer counts across 3 units per side; "
        "ties resolve to first maximum S1<S2<S3"
    )
    dominant.to_csv(
        TABLES / "core_switch_local_dominant_allele_changes.tsv",
        sep="\t",
        index=False,
    )


def write_provenance() -> None:
    provenance = {
        "analysis": "corrected core-only N-TRIAGE switch/ACRO1 audit",
        "prefix": str(PREFIX),
        "seed": SEED,
        "signflip_draws": N_SIGNFLIP,
        "bootstrap_draws": N_BOOTSTRAP,
        "software": {
            "python_numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "statsmodels": statsmodels.__version__,
        },
        "input_sha256": {
            str(path): sha256(path)
            for path in [MANIFEST, ARRAYS, CLUSTERS, ACRO1, SWITCHES, BOUNDARIES]
        },
        "definitions": {
            "evaluable_boundary": (
                "boundary represented in native switch_boundaries.tsv; for this run, "
                "sum(complete_core_units_per_array - 5)"
            ),
            "collapsed_switch": (
                "native switches.tsv event after N-TRIAGE greedily retains passing "
                "candidates separated by at least 3 units"
            ),
            "switch_flank": (
                "pooled raw ACRO1 monomer counts over 3 complete units left and 3 "
                "complete units right of each collapsed switch"
            ),
            "TVD": "0.5 * sum(abs(right_Sj_fraction - left_Sj_fraction), j=1..3)",
            "matched_switch_rate": (
                "direct donor match of mean array rates within individual and status; "
                "not donor-plus-phase matched"
            ),
        },
        "ntriage_switch_defaults": {
            "units_per_side": 3,
            "phase_sketch_size": 512,
            "minimum_within_flank_coherence": 0.90,
            "maximum_cross_flank_ANI": 0.92,
            "minimum_coherence_minus_cross_drop": 0.04,
            "maximum_adjacent_ANI": 0.92,
            "collapse_minimum_spacing_units": 3,
        },
        "interpretation": (
            "Associations are descriptive. Recombination labels derive from the same "
            "array sequence context and do not establish causal effects."
        ),
    }
    with (TABLES / "core_switch_acro1_audit_provenance.json").open("w") as handle:
        json.dump(provenance, handle, indent=2)
        handle.write("\n")


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    data, _, acro, switches, boundaries = build_array_frame()
    write_pooled_counts(data, switches, boundaries)
    write_matched_array_metrics(data)
    write_models(data)
    primary = build_switch_flank_table(data, acro, switches, boundaries)
    write_switch_flank_tests(primary)
    write_provenance()
    print(
        f"Audited {len(data):,} arrays, {int(data.complete_acres.sum()):,} core units, "
        f"{int(data.n_eval.sum()):,} evaluable boundaries, and {len(switches):,} switches"
    )


if __name__ == "__main__":
    main()
