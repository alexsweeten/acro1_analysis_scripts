#!/usr/bin/env python3
"""Stratify fixed ACRO L1PA10-promoter divergence by native NTRIAGE cluster.

The unit of analysis is an array x NTRIAGE-cluster group.  HMM-defined complete
ACRO units from the promoter-purifying-selection pipeline are matched to native
NTRIAGE units by array-local coordinate overlap.  Divergence is then calculated
only among complete units assigned to the same cluster in the same array.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import math
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr


HERE = Path(__file__).resolve().parent
DATA_ROOT = Path(os.environ.get("ACRO1_WORKSPACE", HERE / "inputs"))
PIPE = Path(os.environ.get("ACRO1_PROMOTER_PIPELINE_DATA", DATA_ROOT / "promoter_purifying_2794"))
PIPE_SCRIPTS = HERE.parent / "pipeline"
sys.path.insert(0, str(PIPE_SCRIPTS))
import pp_common as C  # noqa: E402
import pp_core as K  # noqa: E402
from run_array import parse_tblout  # noqa: E402

NTRIAGE_UNITS = Path(os.environ.get("ACRO1_NATIVE_HPRC_UNITS", DATA_ROOT / "telomeric_2794_core.clusters.tsv"))
SIDE_METRICS = Path(os.environ.get("ACRO1_SPLIT_SIDE_METRICS", DATA_ROOT / "array_side_metrics.tsv.gz"))

OUT = Path(__file__).resolve().parent
TABLES = OUT / "tables"
FIGURES = OUT / "figures"
LOGS = OUT / "logs"
CACHE = OUT / "cache"

SEED = 20261027
MIN_CLUSTER_UNITS = 3
MIN_OVERLAP_FRAC = 0.50
N_SIGNFLIP = 200_000
N_BOOT = 20_000
MIN_PAIR_ARRAYS = 10
MIN_PAIR_INDIVIDUALS = 8


def deterministic_rng(label: str) -> np.random.Generator:
    digest = hashlib.sha256(f"{SEED}|{label}".encode()).digest()
    return np.random.default_rng(int.from_bytes(digest[:8], "little"))


def bh(values: pd.Series) -> pd.Series:
    x = values.to_numpy(float)
    ok = np.isfinite(x)
    out = np.full(len(x), np.nan)
    if not ok.any():
        return pd.Series(out, index=values.index)
    y = x[ok]
    order = np.argsort(y)
    ranked = y[order] * len(y) / np.arange(1, len(y) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty_like(ranked)
    q[order] = np.minimum(ranked, 1.0)
    out[np.flatnonzero(ok)] = q
    return pd.Series(out, index=values.index)


def signflip(values: np.ndarray, label: str, alternative: str = "two-sided") -> float:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not len(values):
        return math.nan
    observed = values.mean()
    rng = deterministic_rng(f"signflip|{label}")
    extreme = 0
    batch = 5_000
    for start in range(0, N_SIGNFLIP, batch):
        n = min(batch, N_SIGNFLIP - start)
        means = (rng.choice(np.array([-1.0, 1.0]), size=(n, len(values))) * values).mean(axis=1)
        if alternative == "greater":
            extreme += int(np.sum(means >= observed - 1e-15))
        elif alternative == "less":
            extreme += int(np.sum(means <= observed + 1e-15))
        else:
            extreme += int(np.sum(np.abs(means) >= abs(observed) - 1e-15))
    return (extreme + 1) / (N_SIGNFLIP + 1)


def bootstrap_mean(values: np.ndarray, label: str) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not len(values):
        return math.nan, math.nan
    rng = deterministic_rng(f"bootstrap|{label}")
    estimates = np.empty(N_BOOT)
    for index in range(N_BOOT):
        estimates[index] = rng.choice(values, size=len(values), replace=True).mean()
    lo, hi = np.quantile(estimates, [0.025, 0.975])
    return float(lo), float(hi)


def pooled_divergence(cs: dict[str, object], mask: np.ndarray) -> float:
    callable_columns = np.asarray(cs["ncall"]) >= C.MIN_CALLABLE_FRAC * int(cs["n"]) - 1e-9
    use = mask & callable_columns
    pairs = float(np.asarray(cs["pairs"])[use].sum())
    return float(np.asarray(cs["diffs"])[use].sum() / pairs) if pairs else math.nan


def local_l1pa10_test(encoded: np.ndarray) -> dict[str, object]:
    l1 = encoded[:, : C.L1PA10_END_FULL]
    cs = K.col_stats(l1)
    wm = K.window_metrics(cs)
    cpg = K.cpg_gapless(l1)
    starts = np.arange(1, len(wm["divergence"]) + 1)
    ends = starts + C.WIN - 1
    p0, p1 = C.PROMOTER_FULL
    promoter_index = p0 - 1
    is_promoter = starts == p0
    overlaps = ~((ends < p0) | (starts > p1)) & ~is_promoter
    reference = (wm["gc"][promoter_index], cpg[promoter_index], wm["occupancy"][promoter_index])
    candidates = (~is_promoter) & (~overlaps) & wm["eligible"] & np.isfinite(wm["divergence"])
    matched = candidates & K.match(wm["gc"], cpg, wm["occupancy"], reference, C.TOL_CPG)
    promoter = float(wm["divergence"][promoter_index])
    summary = K.summarise_vs_controls(promoter, wm["divergence"][matched])
    return {
        "promoter_divergence": promoter,
        "control_median": float(summary["control_median"]),
        "effect": float(summary["effect"]),
        "ratio": float(summary["ratio"]),
        "p_matched": float(summary["p"]),
        "n_controls": int(summary["n_controls"]),
        "status": "estimable" if summary["n_controls"] else "no_composition_matched_controls",
    }


def distal_backbone_test(encoded: np.ndarray, backbone_mask: pd.DataFrame) -> dict[str, object]:
    _, windows, promoter_index, _ = K.backbone_windows(encoded, backbone_mask)
    promoter = float(windows["divergence"][promoter_index])
    controls = windows["divergence"][windows["matched_primary_distal"]]
    summary = K.summarise_vs_controls(promoter, controls)
    return {
        "promoter_divergence": promoter,
        "control_median": float(summary["control_median"]),
        "effect": float(summary["effect"]),
        "ratio": float(summary["ratio"]),
        "p_matched": float(summary["p"]),
        "n_controls": int(summary["n_controls"]),
        "status": "estimable" if summary["n_controls"] else "no_composition_matched_controls",
    }


def cluster_metrics(sequences: list[str], backbone_mask: pd.DataFrame) -> dict[str, object]:
    encoded = K.encode(sequences)
    cs = K.col_stats(encoded)
    distal_mask = backbone_mask["primary_distal_backbone_ok"].astype(bool).to_numpy()
    promoter_mask = np.zeros(C.FULL_COLS, dtype=bool)
    promoter_mask[C.PROMOTER_FULL[0] - 1 : C.PROMOTER_FULL[1]] = True
    l1_mask = np.zeros(C.FULL_COLS, dtype=bool)
    l1_mask[: C.L1PA10_END_FULL] = True
    l1_mask[C.PROMOTER_FULL[0] - 1 : C.PROMOTER_FULL[1]] = False
    output: dict[str, object] = {
        "n_units": len(sequences),
        "raw_promoter_divergence": pooled_divergence(cs, promoter_mask),
        "distal_backbone_divergence": pooled_divergence(cs, distal_mask),
        "l1pa10_nonpromoter_divergence": pooled_divergence(cs, l1_mask),
    }
    for prefix, result in (
        ("local", local_l1pa10_test(encoded)),
        ("backbone", distal_backbone_test(encoded, backbone_mask)),
    ):
        output.update({f"{prefix}_{key}": value for key, value in result.items()})
    for background in ("distal_backbone_divergence", "l1pa10_nonpromoter_divergence"):
        value = float(output[background])
        output[f"promoter_depletion_scaled_by_{background}"] = (
            1.0 - float(output["raw_promoter_divergence"]) / value if value > 0 else math.nan
        )
    return output


def overlap_assignment(unit: pd.Series, native: pd.DataFrame) -> dict[str, object]:
    start, end = int(unit.local_start), int(unit.local_end)
    overlap = np.maximum(
        0,
        np.minimum(end, native["end"].to_numpy(int)) -
        np.maximum(start, native["start"].to_numpy(int)) + 1,
    )
    index = int(np.argmax(overlap))
    row = native.iloc[index]
    denom = min(end - start + 1, int(row.end) - int(row.start) + 1)
    fraction = float(overlap[index] / denom) if denom > 0 else 0.0
    return {
        "cluster": int(row.cluster),
        "native_unit_in_array": int(row.unit_in_array),
        "native_start": int(row.start),
        "native_end": int(row.end),
        "overlap_bp": int(overlap[index]),
        "overlap_fraction_shorter": fraction,
        "mapping_pass": bool(overlap[index] > 0 and fraction >= MIN_OVERLAP_FRAC),
    }


def read_one_fasta(path: Path) -> str:
    records = C.read_fasta(path)
    if len(records) != 1:
        raise ValueError(f"Expected one FASTA record in {path}, found {len(records)}")
    return records[0][1].upper()


def load_or_reconstruct_units(base: Path, safe_id: str, work: Path) -> tuple[pd.DataFrame, list[tuple[str, str]], str]:
    qc_path = base / "units_qc.tsv"
    units_path = base / "units_transcription.fa"
    if qc_path.exists() and units_path.exists():
        qc = pd.read_csv(qc_path, sep="\t")
        qc = qc[qc["included_in_structural_msa"].astype(str).str.lower().isin(["true", "1"])].copy()
        records = [
            (name.split()[0], sequence.upper()) for name, sequence in C.read_fasta(units_path)
        ]
        return qc, records, "saved_units"

    fasta = PIPE / "input_fasta" / f"{safe_id}.fa"
    tbl = work / "nhmmer.tbl"
    with open(work / "nhmmer.out", "w") as out, open(work / "nhmmer.err", "w") as err:
        subprocess.run(
            [
                "nhmmer", "--dna", "--tblout", str(tbl), "--noali", "-E", C.NHMMER_E,
                "--cpu", "1", str(C.HMM), str(fasta),
            ],
            stdout=out,
            stderr=err,
            check=True,
        )
    hmm_len = int(next(line.split()[1] for line in open(C.HMM) if line.startswith("LENG")))
    chains = K.build_chains(K.dedup_hits(parse_tblout(tbl, hmm_len)))
    sequence = read_one_fasta(fasta)
    rows, records = [], []
    index = 0
    for strand in "+-":
        if strand not in chains:
            continue
        for unit in K.build_units(sequence, chains[strand], strand):
            index += 1
            length = len(unit["seq"])
            if not C.UNIT_MIN <= length <= C.UNIT_MAX:
                continue
            unit_id = f"u{index:03d}{'p' if strand == '+' else 'm'}"
            rows.append({
                "unit_id": unit_id,
                "local_start": int(unit["local_start"]),
                "local_end": int(unit["local_end"]),
                "included_in_structural_msa": True,
            })
            records.append((unit_id, str(unit["seq"]).upper()))
    if len(records) < C.MIN_UNITS_BACKBONE:
        raise ValueError(f"Only {len(records)} complete units reconstructed")
    return pd.DataFrame(rows), records, "reconstructed_units"


def mapped_unit_sequences(
    base: Path, records: list[tuple[str, str]], work: Path
) -> tuple[dict[str, str], str]:
    existing = base / "backbone" / "units_reference_mapped.fa"
    if existing.exists():
        return {
            name.split()[0]: sequence.upper() for name, sequence in C.read_fasta(existing)
        }, "existing_mapped"
    if not records:
        raise ValueError("No complete unit sequences available")
    source = work / "units.fa"
    C.write_fasta(source, records)
    output = work / "units_plus_reference.fa"
    with open(output, "w") as out, open(work / "mafft.err", "w") as err:
        subprocess.run(
            [
                "mafft", "--quiet", "--thread", "1", "--retree", "1", "--maxiterate", "0",
                "--add", str(source), "--keeplength", str(C.REF_ALN9090),
            ],
            stdout=out,
            stderr=err,
            check=True,
        )
    mapped = {name.split()[0]: sequence.upper() for name, sequence in C.read_fasta(output)}
    selected = {name: mapped[name] for name, _ in records}
    if any(len(sequence) != C.FULL_COLS for sequence in selected.values()):
        raise ValueError("MAFFT did not retain the fixed 9,090-column reference")
    return selected, "mapped_from_saved_units"


def analyze_array(task: dict[str, object]) -> dict[str, object]:
    source = str(task["source_name"])
    safe_id = str(task["safe_id"])
    base = PIPE / "per_array" / safe_id
    result: dict[str, object] = {
        "source_name": source,
        "safe_id": safe_id,
        "individual_id": task["individual_id"],
        "haplotype": task["haplotype"],
        "has_recombination": task["has_recombination"],
        "status": "pending",
        "groups": [],
        "mappings": [],
    }
    cache = CACHE / f"{safe_id}.json"
    if cache.exists():
        return json.loads(cache.read_text())
    try:
        native = pd.DataFrame(task["native_units"])
        with tempfile.TemporaryDirectory(prefix="acro_cluster_promoter_") as temp:
            work = Path(temp)
            qc, records, unit_source = load_or_reconstruct_units(base, safe_id, work)
            mapped_fasta, alignment_source = mapped_unit_sequences(base, records, work)
            result["unit_source"] = unit_source
            result["alignment_source"] = alignment_source
            for _, unit in qc.iterrows():
                assignment = overlap_assignment(unit, native)
                record = {
                    "source_name": source,
                    "safe_id": safe_id,
                    "individual_id": task["individual_id"],
                    "unit_id": unit.unit_id,
                    "hmm_start": int(unit.local_start),
                    "hmm_end": int(unit.local_end),
                    **assignment,
                }
                result["mappings"].append(record)
            mapping = pd.DataFrame(result["mappings"])
            passed = mapping[mapping.mapping_pass].copy()
            result["n_hmm_units"] = len(qc)
            result["n_cluster_mapped_units"] = len(passed)
            result["mapping_rate"] = len(passed) / len(qc) if len(qc) else math.nan
            mask = pd.read_csv(C.BACKBONE_MASK, sep="\t")
            for cluster, group in passed.groupby("cluster"):
                if len(group) < MIN_CLUSTER_UNITS:
                    continue
                unit_ids = group.unit_id.astype(str).tolist()
                missing = [unit_id for unit_id in unit_ids if unit_id not in mapped_fasta]
                if missing:
                    raise ValueError(f"{len(missing)} mapped unit sequences absent, e.g. {missing[0]}")
                metrics = cluster_metrics([mapped_fasta[unit_id] for unit_id in unit_ids], mask)
                result["groups"].append({
                    "source_name": source,
                    "safe_id": safe_id,
                    "individual_id": task["individual_id"],
                    "haplotype": task["haplotype"],
                    "has_recombination": int(task["has_recombination"]),
                    "cluster": int(cluster),
                    **metrics,
                })
        result["status"] = "valid"
        temporary = cache.with_suffix(f".{os.getpid()}.tmp")
        temporary.write_text(json.dumps(result, allow_nan=True))
        os.replace(temporary, cache)
        return result
    except Exception as exc:
        result["status"] = "error"
        result["error"] = f"{type(exc).__name__}: {exc}"
        return result


def cluster_summary(frame: pd.DataFrame, cohort: str, metrics: list[str]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for cluster, sub in frame.groupby("cluster"):
        for metric in metrics:
            data = sub[["source_name", "individual_id", metric]].dropna()
            individual = data.groupby("individual_id", as_index=False)[metric].mean()
            ci_lo, ci_hi = bootstrap_mean(individual[metric].to_numpy(), f"{cohort}|C{cluster}|{metric}")
            rows.append({
                "cohort": cohort,
                "cluster": f"C{int(cluster)}",
                "metric": metric,
                "n_array_cluster_groups": len(data),
                "n_arrays": data.source_name.nunique(),
                "n_individuals": data.individual_id.nunique(),
                "array_weighted_mean": data[metric].mean(),
                "array_weighted_median": data[metric].median(),
                "array_weighted_q25": data[metric].quantile(0.25),
                "array_weighted_q75": data[metric].quantile(0.75),
                "individual_weighted_mean": individual[metric].mean(),
                "individual_boot_ci95_lo": ci_lo,
                "individual_boot_ci95_hi": ci_hi,
                "fraction_positive": float((data[metric] > 0).mean()),
            })
    return pd.DataFrame(rows)


def correlation_summary(frame: pd.DataFrame, cohort: str) -> pd.DataFrame:
    rows = []
    for cluster, sub in frame.groupby("cluster"):
        data = sub[["raw_promoter_divergence", "distal_backbone_divergence"]].dropna()
        rho, pvalue = (math.nan, math.nan)
        if len(data) >= 3:
            rho, pvalue = spearmanr(
                data.raw_promoter_divergence, data.distal_backbone_divergence
            )
        rows.append({
            "cohort": cohort,
            "cluster": f"C{int(cluster)}",
            "n_array_cluster_groups": len(data),
            "spearman_rho_promoter_vs_distal_backbone": rho,
            "spearman_p": pvalue,
        })
    out = pd.DataFrame(rows)
    out["q_bh_within_cohort"] = bh(out["spearman_p"])
    return out


def estimability_summary(frame: pd.DataFrame, cohort: str) -> pd.DataFrame:
    rows = []
    for cluster, sub in frame.groupby("cluster"):
        rows.append({
            "cohort": cohort,
            "cluster": f"C{int(cluster)}",
            "array_cluster_groups": len(sub),
            "local_l1pa10_estimable": int(sub.local_status.eq("estimable").sum()),
            "local_l1pa10_estimable_fraction": float(sub.local_status.eq("estimable").mean()),
            "distal_backbone_estimable": int(sub.backbone_status.eq("estimable").sum()),
            "distal_backbone_estimable_fraction": float(sub.backbone_status.eq("estimable").mean()),
        })
    return pd.DataFrame(rows)


def one_sample_effect_tests(frame: pd.DataFrame, cohort: str, metrics: list[str]) -> pd.DataFrame:
    rows = []
    for cluster, sub in frame.groupby("cluster"):
        for metric in metrics:
            data = sub[["individual_id", metric]].dropna()
            individual = data.groupby("individual_id", as_index=False)[metric].mean()
            ci_lo, ci_hi = bootstrap_mean(individual[metric].to_numpy(), f"one|{cohort}|C{cluster}|{metric}")
            rows.append({
                "cohort": cohort,
                "cluster": f"C{int(cluster)}",
                "metric": metric,
                "alternative": "greater_than_zero",
                "n_array_cluster_groups": len(data),
                "n_individuals": len(individual),
                "array_weighted_mean": data[metric].mean(),
                "individual_weighted_mean": individual[metric].mean(),
                "individual_boot_ci95_lo": ci_lo,
                "individual_boot_ci95_hi": ci_hi,
                "individual_mean_signflip_p": signflip(
                    individual[metric].to_numpy(), f"one|{cohort}|C{cluster}|{metric}", "greater"
                ),
            })
    out = pd.DataFrame(rows)
    out["q_bh_within_cohort_metric"] = out.groupby(["cohort", "metric"])[
        "individual_mean_signflip_p"
    ].transform(bh)
    return out


def pairwise_tests(frame: pd.DataFrame, cohort: str, metrics: list[str]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    clusters = sorted(frame.cluster.unique())
    for metric in metrics:
        wide = frame.pivot_table(
            index=["source_name", "individual_id"], columns="cluster", values=metric, aggfunc="first"
        )
        for i, a in enumerate(clusters):
            for b in clusters[i + 1 :]:
                if a not in wide or b not in wide:
                    continue
                paired = wide[[a, b]].dropna().reset_index()
                if len(paired) < MIN_PAIR_ARRAYS or paired.individual_id.nunique() < MIN_PAIR_INDIVIDUALS:
                    continue
                paired["difference"] = paired[a] - paired[b]
                individual = paired.groupby("individual_id", as_index=False).difference.mean()
                ci_lo, ci_hi = bootstrap_mean(
                    individual.difference.to_numpy(), f"pair|{cohort}|{metric}|C{a}-C{b}"
                )
                rows.append({
                    "cohort": cohort,
                    "metric": metric,
                    "cluster_a": f"C{int(a)}",
                    "cluster_b": f"C{int(b)}",
                    "n_paired_arrays": len(paired),
                    "n_individuals": paired.individual_id.nunique(),
                    "cluster_a_mean": paired[a].mean(),
                    "cluster_b_mean": paired[b].mean(),
                    "array_mean_a_minus_b": paired.difference.mean(),
                    "array_median_a_minus_b": paired.difference.median(),
                    "individual_mean_a_minus_b": individual.difference.mean(),
                    "individual_boot_ci95_lo": ci_lo,
                    "individual_boot_ci95_hi": ci_hi,
                    "individual_mean_signflip_p_two_sided": signflip(
                        individual.difference.to_numpy(), f"pair|{cohort}|{metric}|C{a}-C{b}", "two-sided"
                    ),
                })
    out = pd.DataFrame(rows)
    if len(out):
        out["q_bh_within_cohort_metric"] = out.groupby(["cohort", "metric"])[
            "individual_mean_signflip_p_two_sided"
        ].transform(bh)
    return out


def make_figure(frame: pd.DataFrame, cohort: str, output_stem: str) -> None:
    metrics = [
        ("raw_promoter_divergence", "Raw promoter divergence"),
        ("distal_backbone_divergence", "Distal-backbone divergence"),
        ("promoter_depletion_scaled_by_distal_backbone_divergence", "Promoter depletion vs distal backbone"),
        ("local_effect", "Local L1PA10 matched effect"),
    ]
    clusters = sorted(frame.cluster.unique())
    colors = plt.get_cmap("tab10")
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    for axis, (metric, title) in zip(axes.flat, metrics):
        positions, labels = [], []
        for position, cluster in enumerate(clusters, 1):
            selected = frame.loc[frame.cluster.eq(cluster), ["individual_id", metric]].dropna()
            if not len(selected):
                continue
            individual = selected.groupby("individual_id")[metric].mean().to_numpy(float)
            estimate = float(individual.mean())
            lo, hi = bootstrap_mean(individual, f"figure|{cohort}|C{cluster}|{metric}")
            positions.append(position)
            labels.append(f"C{int(cluster)}\n{len(selected):,} groups\n{len(individual):,} people")
            axis.errorbar(
                position,
                estimate,
                yerr=np.array([[estimate - lo], [hi - estimate]]),
                fmt="o",
                color=colors((cluster - 1) % 10),
                markeredgecolor="black",
                markeredgewidth=0.5,
                markersize=8,
                capsize=4,
                linewidth=2,
            )
        axis.axhline(0, color="black", linewidth=0.7, alpha=0.5)
        axis.set_xticks(positions, labels)
        axis.set_title(title, loc="left", fontweight="bold")
        axis.grid(axis="y", alpha=0.22)
    fig.suptitle(
        f"Fixed L1PA10-promoter divergence by within-array NTRIAGE cluster\n{cohort.replace('_', ' ')}",
        x=0.07, ha="left", fontsize=15, fontweight="bold",
    )
    fig.tight_layout(rect=(0, 0.01, 1, 0.93))
    fig.savefig(FIGURES / f"{output_stem}.png", dpi=300)
    fig.savefig(FIGURES / f"{output_stem}.pdf")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=min(8, os.cpu_count() or 1))
    args = parser.parse_args()
    for directory in (TABLES, FIGURES, LOGS, CACHE):
        directory.mkdir(parents=True, exist_ok=True)

    manifest = pd.read_csv(PIPE / "manifest.validated.tsv", sep="\t")
    native = pd.read_csv(NTRIAGE_UNITS, sep="\t")
    native_groups = {key: value for key, value in native.groupby("source_name", sort=False)}
    tasks = []
    for _, row in manifest.iterrows():
        source = row.array_id
        if source not in native_groups:
            continue
        tasks.append({
            "source_name": source,
            "safe_id": row.safe_id,
            "individual_id": row.individual_id,
            "haplotype": row.haplotype,
            "has_recombination": row.has_recombination,
            "native_units": native_groups[source][
                ["unit_in_array", "start", "end", "cluster"]
            ].to_dict("records"),
        })

    results: list[dict[str, object]] = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(analyze_array, task) for task in tasks]
        for index, future in enumerate(concurrent.futures.as_completed(futures), 1):
            results.append(future.result())
            if index % 200 == 0 or index == len(futures):
                errors = sum(row["status"] == "error" for row in results)
                print(f"completed {index}/{len(futures)}; errors={errors}", flush=True)

    errors = pd.DataFrame([
        {key: value for key, value in row.items() if key not in ("groups", "mappings")}
        for row in results if row["status"] == "error"
    ], columns=["source_name", "safe_id", "individual_id", "haplotype", "has_recombination", "status", "error"])
    errors.to_csv(LOGS / "errors.tsv", sep="\t", index=False)
    if len(errors):
        raise SystemExit(f"{len(errors)} arrays failed; see logs/errors.tsv")

    mappings = pd.DataFrame([record for row in results for record in row["mappings"]])
    groups = pd.DataFrame([record for row in results for record in row["groups"]])
    mappings.to_csv(TABLES / "hmm_to_ntriage_unit_mapping.tsv.gz", sep="\t", index=False, compression="gzip")
    groups.to_csv(TABLES / "array_cluster_promoter_metrics.tsv.gz", sep="\t", index=False, compression="gzip")

    status = pd.DataFrame([
        {key: value for key, value in row.items() if key not in ("groups", "mappings")}
        for row in results
    ])
    status.to_csv(TABLES / "array_processing_audit.tsv", sep="\t", index=False)

    side = pd.read_csv(SIDE_METRICS, sep="\t")
    prior_valid = set(side.loc[side.status.eq("valid"), "source_name"])
    canonical_valid = set(side.loc[
        side.status.eq("valid") & side.round2_transition.isin(["C1->C2", "C6->C2"]), "source_name"
    ])
    cohorts = {
        "all_2794_arrays": groups,
        "curated_switches_prior_valid_867": groups[groups.source_name.isin(prior_valid)].copy(),
        "canonical_C1_C6_to_C2_prior_valid_824": groups[groups.source_name.isin(canonical_valid)].copy(),
    }
    metrics = [
        "raw_promoter_divergence",
        "distal_backbone_divergence",
        "l1pa10_nonpromoter_divergence",
        "promoter_depletion_scaled_by_distal_backbone_divergence",
        "promoter_depletion_scaled_by_l1pa10_nonpromoter_divergence",
        "local_effect",
        "backbone_effect",
    ]
    summaries, pairs, one_samples, correlations, estimability = [], [], [], [], []
    for cohort, frame in cohorts.items():
        summaries.append(cluster_summary(frame, cohort, metrics))
        pairs.append(pairwise_tests(frame, cohort, metrics))
        one_samples.append(one_sample_effect_tests(frame, cohort, ["local_effect", "backbone_effect"]))
        correlations.append(correlation_summary(frame, cohort))
        estimability.append(estimability_summary(frame, cohort))
    summary = pd.concat(summaries, ignore_index=True)
    pairwise = pd.concat(pairs, ignore_index=True)
    one_sample = pd.concat(one_samples, ignore_index=True)
    summary.to_csv(TABLES / "cluster_summary_statistics.tsv", sep="\t", index=False)
    pairwise.to_csv(TABLES / "within_array_pairwise_cluster_tests.tsv", sep="\t", index=False)
    one_sample.to_csv(TABLES / "cluster_effect_vs_zero_tests.tsv", sep="\t", index=False)
    pd.concat(correlations, ignore_index=True).to_csv(
        TABLES / "promoter_background_correlations.tsv", sep="\t", index=False
    )
    pd.concat(estimability, ignore_index=True).to_csv(
        TABLES / "cluster_control_estimability.tsv", sep="\t", index=False
    )

    make_figure(
        cohorts["all_2794_arrays"], "all_2794_arrays", "promoter_divergence_by_cluster_all_2794"
    )
    make_figure(
        cohorts["curated_switches_prior_valid_867"],
        "curated_switches_prior_valid_867",
        "promoter_divergence_by_cluster_switch_arrays",
    )

    metadata = {
        "arrays_processed": len(results),
        "array_cluster_groups_min_3_units": len(groups),
        "hmm_units": len(mappings),
        "mapped_hmm_units": int(mappings.mapping_pass.sum()),
        "unit_mapping_rate": float(mappings.mapping_pass.mean()),
        "minimum_units_per_array_cluster_group": MIN_CLUSTER_UNITS,
        "minimum_overlap_fraction_of_shorter_unit": MIN_OVERLAP_FRAC,
        "signflip_draws": N_SIGNFLIP,
        "bootstrap_draws": N_BOOT,
        "seed": SEED,
        "fixed_promoter_full_alignment_columns_1based_inclusive": list(C.PROMOTER_FULL),
        "interpretation": "Divergence among assembled repeat copies is not a direct de novo mutation rate.",
    }
    (OUT / "analysis_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
