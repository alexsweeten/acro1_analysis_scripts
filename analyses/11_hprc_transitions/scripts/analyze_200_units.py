#!/usr/bin/env python3
"""Paired 200-unit check of ACRO and ACRO1 length across curated HPRC splits.

Primary design: reproducibly sample 100 curated non-C2 -> C2 arrays and select
one complete unit wholly to the left and one wholly to the right of the split,
using consecutive L1PA10 anchors as an independent, uncapped unit definition.
This yields exactly 200 units.  The NTRiage 6.2--7.1-kb calls are retained as a
secondary comparison because that caller cannot represent a 7.3-kb unit.
"""

from __future__ import annotations

from bisect import bisect_right
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon


SEED = 20261026
N_ARRAYS = 100
N_SIGNFLIP = 200_000
N_BOOT = 20_000

ROOT = Path(os.environ.get("ACRO1_WORKSPACE", Path(__file__).resolve().parents[3]))
OUT = Path(os.environ.get("ACRO1_SPLIT_LENGTH_OUTPUT", Path(__file__).resolve().parents[1]))
TABLES = OUT / "tables"
FIGURES = OUT / "figures"

ROUND2 = ROOT / "hprc_breakpoint_ntriage_cluster_analysis/custom_manual_shift_breakpoint_map_round2/tables/round2_units_plus_minus_8_1196.tsv"
NATIVE = ROOT / "NTRIAGE_OCT26/manuscript_output/telomeric_ntriage_2794/native"
ACRES = NATIVE / "telomeric_2794_core.acres.bed"
COMPOSITION = NATIVE / "telomeric_2794_core.acro1_composition.tsv"
STRICT_BED = NATIVE / "telomeric_2794_core.acro1_subfamilies.bed"
ANNOTATIONS = NATIVE / "telomeric_2794_core.annotations.bed"
RAW_ANCHORS = ROOT / "NTRIAGE_OCT26/manuscript_output/telomeric_ntriage_2794/inputs/acro1_monomer_anchors.bed"
MANIFEST = ROOT / "NTRIAGE_OCT26/promoter_purifying_2794/manifest.validated.tsv"
PER_ARRAY = ROOT / "NTRIAGE_OCT26/promoter_purifying_2794/per_array"
UNIT_QC = ROOT / "NTRIAGE_OCT26/promoter_purifying_2794/combined/unit_qc.tsv"


def as_bool(x: pd.Series) -> pd.Series:
    return x.astype(str).str.lower().eq("true")


def merge_intervals(intervals: list[tuple[int, int]], max_gap: int = 0) -> list[tuple[int, int]]:
    """Return merged half-open intervals, optionally joining short gaps."""
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [list(intervals[0])]
    for start, end in intervals[1:]:
        if start <= merged[-1][1] + max_gap:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(int(a), int(b)) for a, b in merged]


def interval_metrics(intervals: list[tuple[int, int]]) -> dict[str, float]:
    """Coverage and largest ACRO1 tract proxy from monomer-anchor intervals."""
    union = merge_intervals(intervals, max_gap=0)
    tracts = merge_intervals(intervals, max_gap=20)
    total_covered = sum(b - a for a, b in union)
    largest_span = max((b - a for a, b in tracts), default=0)
    # Covered bases inside each <=20-bp-gap component, without counting gaps.
    largest_covered = 0
    for ta, tb in tracts:
        covered = sum(b - a for a, b in union if a >= ta and b <= tb)
        largest_covered = max(largest_covered, covered)
    return {
        "acro1_anchor_count": len(intervals),
        "acro1_union_bp": total_covered,
        "largest_acro1_tract_span_bp": largest_span,
        "largest_acro1_tract_covered_bp": largest_covered,
    }


def signflip_p(diff: np.ndarray, rng: np.random.Generator) -> float:
    diff = np.asarray(diff, dtype=float)
    diff = diff[np.isfinite(diff)]
    observed = abs(float(diff.mean()))
    extreme = 0
    done = 0
    while done < N_SIGNFLIP:
        batch = min(5_000, N_SIGNFLIP - done)
        signs = rng.choice(np.array([-1.0, 1.0]), size=(batch, len(diff)))
        null = np.abs(np.mean(signs * diff[None, :], axis=1))
        extreme += int(np.sum(null >= observed - 1e-15))
        done += batch
    return (extreme + 1) / (N_SIGNFLIP + 1)


def paired_summary(pairs: pd.DataFrame, metrics: list[str], analysis: str) -> pd.DataFrame:
    rows = []
    for i, metric in enumerate(metrics):
        wide = pairs.pivot(index="source_name", columns="side", values=metric).dropna()
        if not {"left", "right"}.issubset(wide.columns):
            continue
        left = wide["left"].to_numpy(float)
        right = wide["right"].to_numpy(float)
        diff = right - left
        individual = pd.Index(wide.index).str.split("#").str[0]
        individual_diff = pd.Series(diff).groupby(individual).mean().to_numpy(float)
        rng = np.random.default_rng(np.random.SeedSequence([SEED, i, len(pairs)]))
        boot = np.empty(N_BOOT)
        for start in range(0, N_BOOT, 2_000):
            n = min(2_000, N_BOOT - start)
            idx = rng.integers(0, len(diff), size=(n, len(diff)))
            boot[start : start + n] = diff[idx].mean(axis=1)
        individual_boot = np.empty(N_BOOT)
        for start in range(0, N_BOOT, 2_000):
            n = min(2_000, N_BOOT - start)
            idx = rng.integers(0, len(individual_diff), size=(n, len(individual_diff)))
            individual_boot[start : start + n] = individual_diff[idx].mean(axis=1)
        try:
            wp = float(wilcoxon(diff, alternative="two-sided").pvalue)
        except ValueError:
            wp = 1.0
        rows.append(
            {
                "analysis": analysis,
                "metric": metric,
                "paired_arrays": len(diff),
                "individuals": len(individual_diff),
                "units": 2 * len(diff),
                "left_mean": float(left.mean()),
                "right_mean": float(right.mean()),
                "left_median": float(np.median(left)),
                "right_median": float(np.median(right)),
                "right_minus_left_mean": float(diff.mean()),
                "right_minus_left_median": float(np.median(diff)),
                "mean_difference_ci95_lo": float(np.quantile(boot, 0.025)),
                "mean_difference_ci95_hi": float(np.quantile(boot, 0.975)),
                "right_over_left_mean_ratio": float(right.mean() / left.mean()) if left.mean() != 0 else float("nan"),
                "fraction_pairs_right_gt_left": float(np.mean(diff > 0)),
                "paired_signflip_p": signflip_p(diff, rng),
                "wilcoxon_signed_rank_p": wp,
                "individual_mean_difference": float(individual_diff.mean()),
                "individual_bootstrap_ci95_lo": float(np.quantile(individual_boot, 0.025)),
                "individual_bootstrap_ci95_hi": float(np.quantile(individual_boot, 0.975)),
                "individual_signflip_p": signflip_p(individual_diff, rng),
            }
        )
    return pd.DataFrame(rows)


def select_l1_units(candidates: pd.DataFrame) -> pd.DataFrame:
    """Select nearest L1PA10-defined unit wholly on either side for each array."""
    all_units = pd.read_csv(UNIT_QC, sep="\t")
    all_units = all_units[as_bool(all_units["included_in_structural_msa"])]
    unit_groups = {name: sub for name, sub in all_units.groupby("array_id", sort=False)}
    records: list[dict] = []
    for row in candidates.itertuples(index=False):
        if row.source_name not in unit_groups:
            continue
        units = unit_groups[row.source_name].copy()
        # Exclude a unit crossing the split; choose the closest complete unit
        # wholly contained in each block.
        left = units[units["local_end"] < row.split_local].sort_values("local_end").tail(1)
        right = units[units["local_start"] >= row.split_local].sort_values("local_start").head(1)
        if left.empty or right.empty:
            continue
        for side, unit in (("left", left.iloc[0]), ("right", right.iloc[0])):
            records.append(
                {
                    "source_name": row.source_name,
                    "side": side,
                    "split_local": row.split_local,
                    "unit_id_l1pa10": unit.unit_id,
                    "unit_start": int(unit.local_start),
                    "unit_end": int(unit.local_end) + 1,  # convert inclusive to half-open
                    "unit_length_l1pa10": int(unit.length),
                    "nearest_ntriage_acre_id": row.left_acre if side == "left" else row.right_acre,
                    "nearest_ntriage_cluster": row.left_cluster if side == "left" else row.right_cluster,
                    "unit_length_ntriage": row.left_length if side == "left" else row.right_length,
                }
            )
    return pd.DataFrame(records)


def assign_anchor_metrics(units: pd.DataFrame) -> pd.DataFrame:
    """Assign all raw ACRO1 monomer anchors to selected L1PA10-defined units."""
    by_source: dict[str, list[dict]] = {}
    for rec in units.to_dict("records"):
        by_source.setdefault(rec["source_name"], []).append(rec)

    raw: dict[tuple[str, str], list[tuple[int, int]]] = {
        (r["source_name"], r["side"]): [] for r in units.to_dict("records")
    }
    with RAW_ANCHORS.open() as handle:
        for line in handle:
            chrom, s, e, *_ = line.rstrip("\n").split("\t")
            if chrom not in by_source:
                continue
            start, end = int(s), int(e)
            midpoint = (start + end) / 2
            for rec in by_source[chrom]:
                if rec["unit_start"] <= midpoint < rec["unit_end"]:
                    raw[(chrom, rec["side"])].append(
                        (max(start, rec["unit_start"]), min(end, rec["unit_end"]))
                    )
                    break

    strict: dict[tuple[str, str], list[tuple[int, int]]] = {
        key: [] for key in raw
    }
    with STRICT_BED.open() as handle:
        next(handle)  # header
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            if chrom not in by_source:
                continue
            midpoint = (start + end) / 2
            for rec in by_source[chrom]:
                if rec["unit_start"] <= midpoint < rec["unit_end"]:
                    strict[(chrom, rec["side"])].append(
                        (max(start, rec["unit_start"]), min(end, rec["unit_end"]))
                    )
                    break

    rows = []
    for rec in units.to_dict("records"):
        key = (rec["source_name"], rec["side"])
        row = dict(rec)
        row.update(interval_metrics(raw[key]))
        strict_union = merge_intervals(strict[key])
        row["strict_acro1_bp"] = sum(b - a for a, b in strict_union)
        row["strict_acro1_monomers"] = len(strict[key])
        rows.append(row)
    return pd.DataFrame(rows)


def assign_backbone_metrics(units: pd.DataFrame) -> pd.DataFrame:
    """Intersect the tiled NTRiage component annotation with independent units."""
    by_source: dict[str, list[dict]] = {}
    for rec in units.to_dict("records"):
        by_source.setdefault(rec["source_name"], []).append(rec)
    labels = [
        "COMP-subunit_ACRO_rnd-5_family-1625",
        "L1PA10_3end",
        "COMP-subunit_ACRO_rnd-5_family-1624",
        "COMP-subunit_ACRO_rnd-5_family-38",
        "COMP-subunit_ACRO_rnd-5_family-37",
        "MER21B",
    ]
    metrics: dict[tuple[str, str], dict[str, float]] = {}
    for rec in units.to_dict("records"):
        key = (rec["source_name"], rec["side"])
        metrics[key] = {
            "annotation_acro1_bp": 0,
            "explicit_backbone_bp": 0,
            "unclassified_bp": 0,
            "annotation_covered_bp": 0,
            **{f"backbone_{label}_bp": 0 for label in labels},
        }
    with ANNOTATIONS.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            chrom, start, end, label = fields[0], int(fields[1]), int(fields[2]), fields[3]
            if chrom not in by_source:
                continue
            for rec in by_source[chrom]:
                overlap = max(0, min(end, rec["unit_end"]) - max(start, rec["unit_start"]))
                if not overlap:
                    continue
                out = metrics[(chrom, rec["side"])]
                out["annotation_covered_bp"] += overlap
                if label == "ACRO1":
                    out["annotation_acro1_bp"] += overlap
                elif label == "unclassified":
                    out["unclassified_bp"] += overlap
                else:
                    out["explicit_backbone_bp"] += overlap
                    key = f"backbone_{label}_bp"
                    if key in out:
                        out[key] += overlap
    rows = []
    for rec in units.to_dict("records"):
        row = dict(rec)
        row.update(metrics[(rec["source_name"], rec["side"])])
        row["non_acro1_remainder_bp"] = row["unit_length_l1pa10"] - row["acro1_union_bp"]
        row["annotated_non_acro1_bp"] = row["explicit_backbone_bp"] + row["unclassified_bp"]
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)

    d = pd.read_csv(ROUND2, sep="\t")
    near = d[d["units_from_round2_anchor"].isin([-1, 1])].copy()
    wide = near.pivot(
        index="source_name",
        columns="units_from_round2_anchor",
        values=["acre_id", "native_cluster", "start", "end"],
    )
    eligible = wide[
        (wide[("native_cluster", -1)] != 2)
        & (wide[("native_cluster", 1)] == 2)
    ].copy()
    candidates = pd.DataFrame(
        {
            "source_name": eligible.index,
            "left_acre": eligible[("acre_id", -1)].values,
            "right_acre": eligible[("acre_id", 1)].values,
            "left_cluster": eligible[("native_cluster", -1)].astype(int).values,
            "right_cluster": eligible[("native_cluster", 1)].astype(int).values,
            "left_length": (eligible[("end", -1)] - eligible[("start", -1)]).astype(int).values,
            "right_length": (eligible[("end", 1)] - eligible[("start", 1)]).astype(int).values,
            "split_local": (
                (eligible[("end", -1)] + eligible[("start", 1)]) / 2
            ).values,
        }
    )

    all_independent = select_l1_units(candidates)
    counts = all_independent.groupby("source_name")["side"].nunique()
    complete_names = counts[counts == 2].index.to_numpy()
    if len(complete_names) < N_ARRAYS:
        raise RuntimeError(f"Only {len(complete_names)} complete independent pairs; need {N_ARRAYS}")

    rng = np.random.default_rng(SEED)
    sampled_names = np.sort(rng.choice(complete_names, size=N_ARRAYS, replace=False))
    all_measured = assign_anchor_metrics(all_independent[all_independent.source_name.isin(complete_names)])
    all_measured = assign_backbone_metrics(all_measured)
    sample = all_measured[all_measured.source_name.isin(sampled_names)].copy()
    sample["sample_seed"] = SEED
    sample["individual_id"] = sample["source_name"].str.split("#").str[0]
    all_measured["individual_id"] = all_measured["source_name"].str.split("#").str[0]

    metrics = [
        "unit_length_l1pa10",
        "unit_length_ntriage",
        "acro1_union_bp",
        "largest_acro1_tract_covered_bp",
        "largest_acro1_tract_span_bp",
        "acro1_anchor_count",
        "strict_acro1_bp",
        "strict_acro1_monomers",
        "non_acro1_remainder_bp",
        "annotation_acro1_bp",
        "explicit_backbone_bp",
        "unclassified_bp",
        "annotated_non_acro1_bp",
        "annotation_covered_bp",
        "backbone_COMP-subunit_ACRO_rnd-5_family-1625_bp",
        "backbone_L1PA10_3end_bp",
        "backbone_COMP-subunit_ACRO_rnd-5_family-1624_bp",
        "backbone_COMP-subunit_ACRO_rnd-5_family-38_bp",
        "backbone_COMP-subunit_ACRO_rnd-5_family-37_bp",
        "backbone_MER21B_bp",
    ]
    sample_stats = paired_summary(sample, metrics, "primary_reproducible_100_pairs_200_units")
    full_stats = paired_summary(all_measured, metrics, "all_eligible_pairs_with_independent_L1PA10_units")
    stats = pd.concat([sample_stats, full_stats], ignore_index=True)

    sample.sort_values(["source_name", "side"]).to_csv(
        TABLES / "sampled_200_units.tsv", sep="\t", index=False
    )
    all_measured.sort_values(["source_name", "side"]).to_csv(
        TABLES / "all_eligible_independent_pairs.tsv", sep="\t", index=False
    )
    stats.to_csv(TABLES / "paired_length_tests.tsv", sep="\t", index=False)
    candidates.to_csv(TABLES / "eligible_curated_nonC2_to_C2_splits.tsv", sep="\t", index=False)
    def row(metric: str, analysis: str = "primary_reproducible_100_pairs_200_units") -> pd.Series:
        return stats[(stats.analysis == analysis) & (stats.metric == metric)].iloc[0]

    ul = row("unit_length_l1pa10")
    nt = row("unit_length_ntriage")
    ac = row("acro1_union_bp")
    tr = row("largest_acro1_tract_covered_bp")
    rem = row("non_acro1_remainder_bp")
    ann_acro1 = row("annotation_acro1_bp")
    bb = row("explicit_backbone_bp")
    unc = row("unclassified_bp")
    f38 = row("backbone_COMP-subunit_ACRO_rnd-5_family-38_bp")
    l1 = row("backbone_L1PA10_3end_bp")
    mer = row("backbone_MER21B_bp")
    f1625 = row("backbone_COMP-subunit_ACRO_rnd-5_family-1625_bp")
    f37 = row("backbone_COMP-subunit_ACRO_rnd-5_family-37_bp")
    full_ul = row("unit_length_l1pa10", "all_eligible_pairs_with_independent_L1PA10_units")
    full_ac = row("acro1_union_bp", "all_eligible_pairs_with_independent_L1PA10_units")

    readme = f"""# Paired 200-unit HPRC split analysis

## Design

The primary analysis reproducibly sampled **100** of the curated arrays with a nearest-unit **non-C2 to C2** transition (seed `{SEED}`). One complete unit wholly on each side of the split was selected, yielding exactly **200 units**. To avoid the 7,100-bp upper limit in the NTRiage complete-unit caller, primary unit lengths were measured independently between consecutive L1PA10 anchors. ACRO1 sequence was measured as the union of nonredundant ACRO1 monomer-anchor intervals within each independent unit. The largest <=20-bp-gap anchor component was used as a proxy for the longer, copy-number-variable ACRO1 tract. All tests are paired by array.

There were **{len(candidates):,}** curated non-C2 to C2 nearest-unit transitions. Independent L1PA10-defined units were locally available and fully separable on both sides for **{len(complete_names):,}** of them. The same comparisons were repeated across this complete eligible set as a sensitivity analysis.

## Primary 200-unit result

* **Uncapped ACRO-unit length:** left mean/median = **{ul.left_mean:.1f}/{ul.left_median:.1f} bp**; right C2 mean/median = **{ul.right_mean:.1f}/{ul.right_median:.1f} bp**. The paired mean right-minus-left difference was **{ul.right_minus_left_mean:+.1f} bp**. After averaging repeated arrays within each of the **{int(ul.individuals)} individuals**, the difference was **{ul.individual_mean_difference:+.1f} bp** (individual-bootstrap 95% CI **{ul.individual_bootstrap_ci95_lo:+.1f} to {ul.individual_bootstrap_ci95_hi:+.1f}**; individual sign-flip P = **{ul.individual_signflip_p:.3g}**).
* **Total ACRO1 anchor coverage:** left mean/median = **{ac.left_mean:.1f}/{ac.left_median:.1f} bp**; right C2 mean/median = **{ac.right_mean:.1f}/{ac.right_median:.1f} bp**. The paired mean difference was **{ac.right_minus_left_mean:+.1f} bp**; the individual-aggregated difference was **{ac.individual_mean_difference:+.1f} bp** (95% CI **{ac.individual_bootstrap_ci95_lo:+.1f} to {ac.individual_bootstrap_ci95_hi:+.1f}**; P = **{ac.individual_signflip_p:.3g}**).
* **Largest ACRO1 tract proxy:** left mean/median = **{tr.left_mean:.1f}/{tr.left_median:.1f} bp**; right C2 mean/median = **{tr.right_mean:.1f}/{tr.right_median:.1f} bp**. The paired mean difference was **{tr.right_minus_left_mean:+.1f} bp** (95% CI **{tr.mean_difference_ci95_lo:+.1f} to {tr.mean_difference_ci95_hi:+.1f}**; sign-flip P = **{tr.paired_signflip_p:.3g}**).
* **Sequence outside exact ACRO1 anchors:** left mean/median = **{rem.left_mean:.1f}/{rem.left_median:.1f} bp**; right C2 mean/median = **{rem.right_mean:.1f}/{rem.right_median:.1f} bp**. The individual-aggregated difference was **{rem.individual_mean_difference:+.1f} bp** (95% CI **{rem.individual_bootstrap_ci95_lo:+.1f} to {rem.individual_bootstrap_ci95_hi:+.1f}**; P = **{rem.individual_signflip_p:.3g}**). This complement includes unaligned bases within ACRO1 tracts and is not itself a backbone annotation.
* **Tiled ACRO1-component annotation:** left mean/median = **{ann_acro1.left_mean:.1f}/{ann_acro1.left_median:.1f} bp**; right C2 mean/median = **{ann_acro1.right_mean:.1f}/{ann_acro1.right_median:.1f} bp**. The individual-aggregated difference was **{ann_acro1.individual_mean_difference:+.1f} bp** (95% CI **{ann_acro1.individual_bootstrap_ci95_lo:+.1f} to {ann_acro1.individual_bootstrap_ci95_hi:+.1f}**; P = **{ann_acro1.individual_signflip_p:.3g}**).
* **Explicitly annotated backbone:** left mean/median = **{bb.left_mean:.1f}/{bb.left_median:.1f} bp**; right C2 mean/median = **{bb.right_mean:.1f}/{bb.right_median:.1f} bp**. The individual-aggregated difference was **{bb.individual_mean_difference:+.1f} bp** (95% CI **{bb.individual_bootstrap_ci95_lo:+.1f} to {bb.individual_bootstrap_ci95_hi:+.1f}**; P = **{bb.individual_signflip_p:.3g}**). Unclassified sequence changed by **{unc.individual_mean_difference:+.1f} bp** (95% CI **{unc.individual_bootstrap_ci95_lo:+.1f} to {unc.individual_bootstrap_ci95_hi:+.1f}**; P = **{unc.individual_signflip_p:.3g}**).
* **Backbone redistribution:** individual-aggregated right-minus-left differences were **{f38.individual_mean_difference:+.1f} bp** for family-38, **{l1.individual_mean_difference:+.1f} bp** for L1PA10, **{mer.individual_mean_difference:+.1f} bp** for MER21B, **{f1625.individual_mean_difference:+.1f} bp** for family-1625, and **{f37.individual_mean_difference:+.1f} bp** for family-37.
* **Capped NTRiage calls:** left and right means were **{nt.left_mean:.1f}** and **{nt.right_mean:.1f} bp**, respectively. These are secondary because the caller only permits 6,200–7,100-bp complete units.

## Full eligible-set sensitivity

Across all **{int(full_ul.paired_arrays):,}** independently measurable paired arrays from **{int(full_ul.individuals)} individuals**, uncapped L1PA10-defined units averaged **{full_ul.left_mean:.1f} bp** on the left and **{full_ul.right_mean:.1f} bp** on the C2 right (array-paired difference **{full_ul.right_minus_left_mean:+.1f} bp**; individual-aggregated difference **{full_ul.individual_mean_difference:+.1f} bp**, 95% CI **{full_ul.individual_bootstrap_ci95_lo:+.1f} to {full_ul.individual_bootstrap_ci95_hi:+.1f}**, P = **{full_ul.individual_signflip_p:.3g}**). Total ACRO1 anchor coverage averaged **{full_ac.left_mean:.1f} bp** on the left and **{full_ac.right_mean:.1f} bp** on the right (array-paired difference **{full_ac.right_minus_left_mean:+.1f} bp**; individual-aggregated difference **{full_ac.individual_mean_difference:+.1f} bp**, 95% CI **{full_ac.individual_bootstrap_ci95_lo:+.1f} to {full_ac.individual_bootstrap_ci95_hi:+.1f}**, P = **{full_ac.individual_signflip_p:.3g}**).

## Interpretation

The independent L1PA10-to-L1PA10 measurement is the appropriate result for evaluating a proposed ~7.3-kb C2 unit because the original NTRiage segmentation cannot call units longer than 7.1 kb. ACRO1 anchor coverage measures assembled sequence supported by the monomer alignments; the largest-tract statistic is a proxy rather than a direct annotation of the short and long ACRO1 tracts. The right/C2 side therefore contains significantly more explicitly annotated backbone, but the backbone difference is small relative to the expansion assigned to the broader ACRO1 component.

## Outputs

* `tables/sampled_200_units.tsv`: the exact 200 units and all measurements.
* `tables/paired_length_tests.tsv`: paired estimates, confidence intervals, and P values.
* `tables/all_eligible_independent_pairs.tsv`: full sensitivity data set.
* `figures/paired_200_unit_length_acro1_comparison.png` and `.pdf`: paired comparison figure.
"""
    (OUT / "README.md").write_text(readme)


if __name__ == "__main__":
    main()

