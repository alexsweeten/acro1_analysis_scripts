#!/usr/bin/env python3
"""Quantify family-38 abundance, unit length, and candidate PRDM9 motifs.

The annotation-based analysis measures family-38 sequence assigned by N-TRIAGE.
The alignment-based analysis maps complete ACRE sequences to the 1,300-bp
family-38 reference and measures near-full-length homologs directly.

Motif matches are sequence compatibility tests only.  They are not evidence of
PRDM9 occupancy, meiotic double-strand breaks, or recombination activity.
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import re
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import kruskal, mannwhitneyu


FAMILY38 = "COMP-subunit_ACRO_rnd-5_family-38"
SPECIES_ORDER = ["hg002", "mGorGor1", "mPanPan1", "mPanTro3", "mPonAbe1", "mPonPyg2"]
SPECIES_LABELS = {
    "hg002": "Human",
    "mGorGor1": "Gorilla",
    "mPanPan1": "Bonobo",
    "mPanTro3": "Chimpanzee",
    "mPonAbe1": "Sumatran\norangutan",
    "mPonPyg2": "Bornean\norangutan",
}
LINEAGE = {
    "hg002": "Human",
    "mGorGor1": "Gorilla",
    "mPanPan1": "Pan",
    "mPanTro3": "Pan",
    "mPonAbe1": "Orangutan",
    "mPonPyg2": "Orangutan",
}
COLORS = {
    "hg002": "#4C78A8",
    "mGorGor1": "#59A14F",
    "mPanPan1": "#E15759",
    "mPanTro3": "#F28E2B",
    "mPonAbe1": "#B07AA1",
    "mPonPyg2": "#9C755F",
}

# Published motif representations. N is represented by a dot in the regex.
MOTIFS = {
    "human_A_exact_13mer": {
        "display": "CCTCCCTNNCCAC",
        "regex": "CCTCCCT..CCAC",
        "interpretation": "specific human PRDM9-A-associated hotspot motif",
    },
    "human_A_broad_consensus": {
        "display": "NCCNCCNTNNCCNCN",
        "regex": ".CC.CC.T..CC.C.",
        "interpretation": "broad structural consensus compatible with human PRDM9-A",
    },
    "Pan_common_predicted": {
        "display": "AATTNNANTCNTCC",
        "regex": "AATT..A.TC.TCC",
        "interpretation": "computationally predicted common Pan submotif",
    },
    "gorilla_predicted": {
        "display": "CCNANNCCTC",
        "regex": "CC.A..CCTC",
        "interpretation": "computationally predicted gorilla submotif",
    },
}


def parse_args() -> argparse.Namespace:
    analysis_dir = Path(__file__).resolve().parent
    data_root = Path(
        os.environ.get("ACRO1_DATA_ROOT", analysis_dir / "data")
    )
    output_dir = Path(
        os.environ.get("ACRO1_OUTPUT_DIR", analysis_dir / "results")
    )
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root",
        "--data-root",
        dest="root",
        type=Path,
        default=data_root,
        help=(
            "Input-data root containing primate_full/ and references/ "
            "(default: $ACRO1_DATA_ROOT or ANALYSIS_DIR/data)."
        ),
    )
    parser.add_argument(
        "--metadata",
        type=Path,
        help="ACRE metadata TSV (default: ROOT/primate_full/six_primate_acre_metadata.tsv).",
    )
    parser.add_argument(
        "--annotations",
        type=Path,
        help="Component annotation BED (default: ROOT/primate_full/six_primate_full.annotations.bed).",
    )
    parser.add_argument(
        "--acre-fasta",
        type=Path,
        help="Complete-ACRE FASTA (default: ROOT/primate_full/six_primate_full.acres.fasta).",
    )
    parser.add_argument(
        "--component-fasta",
        type=Path,
        help="Component-reference FASTA (default: ROOT/references/dfam_2.fasta).",
    )
    parser.add_argument(
        "--paf",
        type=Path,
        help="Optional precomputed minimap2 PAF or PAF.gz. If omitted, mapping is run or reused from OUTDIR.",
    )
    parser.add_argument(
        "--minimap2",
        default=os.environ.get("MINIMAP2", shutil.which("minimap2") or "minimap2"),
        help="minimap2 executable (default: $MINIMAP2 or executable found on PATH).",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=output_dir,
        help="Output directory (default: $ACRO1_OUTPUT_DIR or ANALYSIS_DIR/results).",
    )
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--bootstrap", type=int, default=10_000)
    parser.add_argument("--seed", type=int, default=38)
    parser.add_argument("--skip-mapping", action="store_true")
    return parser.parse_args()


def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    name = None
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name is not None:
                seqs[name].append(line.upper())
    return {name: "".join(parts) for name, parts in seqs.items()}


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def overlapping_matches(seq: str, pattern: str) -> list[tuple[int, int, str, str]]:
    """Return zero-based half-open coordinates, strand, and matched sequence."""
    hits: list[tuple[int, int, str, str]] = []
    regex = re.compile(f"(?=({pattern}))", re.IGNORECASE)
    rc_pattern = revcomp(pattern.replace(".", "N")).replace("N", ".")
    rc_regex = re.compile(f"(?=({rc_pattern}))", re.IGNORECASE)
    for match in regex.finditer(seq):
        matched = match.group(1).upper()
        hits.append((match.start(), match.start() + len(matched), "+", matched))
    for match in rc_regex.finditer(seq):
        matched = match.group(1).upper()
        hits.append((match.start(), match.start() + len(matched), "-", matched))
    return sorted(set(hits))


def bootstrap_difference(x: np.ndarray, y: np.ndarray, n: int, rng: np.random.Generator) -> tuple[float, float, float]:
    observed = float(np.mean(x) - np.mean(y))
    draws = np.empty(n, dtype=float)
    for i in range(n):
        draws[i] = rng.choice(x, len(x), replace=True).mean() - rng.choice(y, len(y), replace=True).mean()
    low, high = np.quantile(draws, [0.025, 0.975])
    return observed, float(low), float(high)


def compare_groups(frame: pd.DataFrame, value: str, group: str, comparisons: list[tuple[str, list[str], str, list[str]]], bootstrap: int, seed: int) -> pd.DataFrame:
    rows = []
    rng = np.random.default_rng(seed)
    for left_name, left_values, right_name, right_values in comparisons:
        x = frame.loc[frame[group].isin(left_values), value].dropna().to_numpy(float)
        y = frame.loc[frame[group].isin(right_values), value].dropna().to_numpy(float)
        if not len(x) or not len(y):
            continue
        test = mannwhitneyu(x, y, alternative="two-sided")
        diff, lo, hi = bootstrap_difference(x, y, bootstrap, rng)
        rows.append(
            {
                "metric": value,
                "left_group": left_name,
                "right_group": right_name,
                "left_n": len(x),
                "right_n": len(y),
                "left_mean": np.mean(x),
                "right_mean": np.mean(y),
                "mean_difference": diff,
                "bootstrap_95ci_low": lo,
                "bootstrap_95ci_high": hi,
                "mann_whitney_u": test.statistic,
                "mann_whitney_p": test.pvalue,
            }
        )
    result = pd.DataFrame(rows)
    if not result.empty:
        # Benjamini-Hochberg correction within each four-contrast analysis.
        p = result["mann_whitney_p"].to_numpy(float)
        order = np.argsort(p)
        adjusted = np.empty_like(p)
        ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
        ranked = np.minimum.accumulate(ranked[::-1])[::-1]
        adjusted[order] = np.minimum(ranked, 1.0)
        result["mann_whitney_fdr_bh"] = adjusted
    return result


def run_mapping(
    component_fasta: Path,
    acre_fasta: Path,
    outdir: Path,
    minimap2: str,
    threads: int,
    skip: bool,
    precomputed_paf: Path | None,
) -> Path:
    if precomputed_paf is not None:
        if not precomputed_paf.exists():
            raise FileNotFoundError(f"Precomputed PAF not found: {precomputed_paf}")
        return precomputed_paf
    paf = outdir / "family38_all_mappings.paf.gz"
    if paf.exists():
        return paf
    if skip:
        raise FileNotFoundError(f"Mapping file not found: {paf}")
    resolved_minimap2 = shutil.which(minimap2)
    if resolved_minimap2 is None:
        candidate = Path(minimap2).expanduser()
        if not candidate.is_file():
            raise FileNotFoundError(
                f"minimap2 executable not found: {minimap2}. "
                "Use --minimap2 or set MINIMAP2."
            )
        resolved_minimap2 = str(candidate)
    command = [
        resolved_minimap2, "-x", "sr", "-k", "11", "-w", "5", "-N", "50",
        "--secondary=yes", "-c", "-t", str(threads),
        str(component_fasta),
        str(acre_fasta),
    ]
    completed = subprocess.run(command, check=True, capture_output=True, text=True)
    with gzip.open(paf, "wt") as handle:
        handle.write(completed.stdout)
    (outdir / "minimap2.log").write_text(completed.stderr)
    return paf


def parse_family38_paf(path: Path, metadata: pd.DataFrame) -> pd.DataFrame:
    meta = metadata.set_index("acre_id")
    records = []
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            if len(fields) < 12 or fields[5] != FAMILY38:
                continue
            acre_id = fields[0].split("|")[0]
            if acre_id not in meta.index:
                continue
            tags = {part[:2]: part[5:] for part in fields[12:] if len(part) > 5 and part[2:5] == ":A:"}
            record = {
                "acre_id": acre_id,
                "query_length": int(fields[1]),
                "query_start": int(fields[2]),
                "query_end": int(fields[3]),
                "strand": fields[4],
                "target_length": int(fields[6]),
                "target_start": int(fields[7]),
                "target_end": int(fields[8]),
                "matching_bases": int(fields[9]),
                "alignment_block_length": int(fields[10]),
                "mapq": int(fields[11]),
                "tp": tags.get("tp", ""),
            }
            record["query_span"] = record["query_end"] - record["query_start"]
            record["target_span"] = record["target_end"] - record["target_start"]
            record["target_coverage"] = record["target_span"] / record["target_length"]
            record["identity"] = record["matching_bases"] / record["alignment_block_length"]
            for column in ["array_id", "species_id", "species_name", "original_chrom"]:
                record[column] = meta.at[acre_id, column]
            records.append(record)
    return pd.DataFrame(records)


def deduplicate_mappings(hits: pd.DataFrame) -> pd.DataFrame:
    kept = []
    for _, group in hits.groupby("acre_id", sort=False):
        accepted: list[pd.Series] = []
        group = group.assign(_score=group["matching_bases"] + group["identity"] * 100).sort_values(
            ["_score", "query_span"], ascending=False
        )
        for _, row in group.iterrows():
            duplicate = False
            for prior in accepted:
                overlap = max(0, min(row.query_end, prior.query_end) - max(row.query_start, prior.query_start))
                if overlap / min(row.query_span, prior.query_span) >= 0.8:
                    duplicate = True
                    break
            if not duplicate:
                accepted.append(row)
        kept.extend(accepted)
    if not kept:
        return hits.iloc[0:0].copy()
    result = pd.DataFrame(kept).drop(columns=["_score"], errors="ignore")
    return result.sort_values(["acre_id", "query_start"]).reset_index(drop=True)


def annotation_analysis(
    annotations_path: Path,
    outdir: Path,
    metadata: pd.DataFrame,
    bootstrap: int,
    seed: int,
):
    annotations = pd.read_csv(
        annotations_path,
        sep="\t",
        header=None,
        names=[
            "synthetic_chrom", "start", "end", "component", "score", "strand",
            "thick_start", "thick_end", "rgb", "acre_id", "array_id",
            "raw_score", "normalized_score", "coverage",
        ],
    )
    family = annotations.loc[annotations.component.eq(FAMILY38)].copy()
    family["tract_bp"] = family.end - family.start
    family = family.merge(
        metadata[["acre_id", "species_id", "species_name", "original_chrom"]],
        on="acre_id", how="left", validate="many_to_one",
    )
    family.to_csv(outdir / "family38_annotation_tracts.tsv.gz", sep="\t", index=False, compression="gzip")

    per_acre_counts = family.groupby("acre_id").agg(
        family38_bp=("tract_bp", "sum"), family38_tracts=("tract_bp", "size")
    )
    per_acre = metadata[["acre_id", "array_id", "species_id", "species_name", "original_chrom", "length"]].merge(
        per_acre_counts, left_on="acre_id", right_index=True, how="left", validate="one_to_one"
    )
    per_acre[["family38_bp", "family38_tracts"]] = per_acre[["family38_bp", "family38_tracts"]].fillna(0).astype(int)
    per_acre["family38_fraction"] = per_acre.family38_bp / per_acre.length
    per_acre["family38_reference_copy_equivalents"] = per_acre.family38_bp / 1300.0
    per_acre["lineage"] = per_acre.species_id.map(LINEAGE)
    per_acre.to_csv(outdir / "family38_per_acre.tsv.gz", sep="\t", index=False, compression="gzip")

    per_species = per_acre.groupby(["species_id", "species_name"], sort=False).agg(
        acres=("acre_id", "size"),
        acres_with_family38=("family38_bp", lambda x: int((x > 0).sum())),
        total_acre_bp=("length", "sum"),
        total_family38_bp=("family38_bp", "sum"),
        mean_family38_bp_per_acre=("family38_bp", "mean"),
        median_family38_bp_per_acre=("family38_bp", "median"),
        mean_family38_fraction_per_acre=("family38_fraction", "mean"),
        median_family38_fraction_per_acre=("family38_fraction", "median"),
        total_family38_tracts=("family38_tracts", "sum"),
        mean_family38_tracts_per_acre=("family38_tracts", "mean"),
    ).reset_index()
    per_species["pooled_family38_fraction"] = per_species.total_family38_bp / per_species.total_acre_bp
    per_species["total_reference_copy_equivalents"] = per_species.total_family38_bp / 1300.0
    tract_stats = family.groupby("species_id").tract_bp.agg(["mean", "median", "min", "max", "std"]).reset_index()
    tract_stats = tract_stats.rename(columns={c: f"annotation_tract_{c}_bp" for c in ["mean", "median", "min", "max", "std"]})
    per_species = per_species.merge(tract_stats, on="species_id", how="left")
    per_species["_order"] = per_species.species_id.map({s: i for i, s in enumerate(SPECIES_ORDER)})
    per_species = per_species.sort_values("_order").drop(columns="_order")
    per_species.to_csv(outdir / "family38_annotation_by_species.tsv", sep="\t", index=False)

    per_array = per_acre.groupby(["species_id", "species_name", "array_id"], sort=False).agg(
        acres=("acre_id", "size"), total_acre_bp=("length", "sum"), family38_bp=("family38_bp", "sum")
    ).reset_index()
    per_array["family38_fraction"] = per_array.family38_bp / per_array.total_acre_bp
    per_array["mean_family38_bp_per_acre"] = per_array.family38_bp / per_array.acres
    per_array["lineage"] = per_array.species_id.map(LINEAGE)
    per_array.to_csv(outdir / "family38_per_array.tsv", sep="\t", index=False)

    comparisons = [
        ("Pan", ["mPanPan1", "mPanTro3"], "all non-Pan", ["hg002", "mGorGor1", "mPonAbe1", "mPonPyg2"]),
        ("Pan", ["mPanPan1", "mPanTro3"], "human", ["hg002"]),
        ("Pan", ["mPanPan1", "mPanTro3"], "gorilla", ["mGorGor1"]),
        ("Pan", ["mPanPan1", "mPanTro3"], "orangutans", ["mPonAbe1", "mPonPyg2"]),
    ]
    stats = compare_groups(per_array, "family38_fraction", "species_id", comparisons, bootstrap, seed)
    arrays = [per_array.loc[per_array.species_id.eq(s), "family38_fraction"].to_numpy() for s in SPECIES_ORDER]
    kw = kruskal(*arrays)
    stats["six_species_kruskal_h"] = kw.statistic
    stats["six_species_kruskal_p"] = kw.pvalue
    return family, per_acre, per_array, per_species, stats


def mapping_analysis(
    component_fasta: Path,
    acre_fasta: Path,
    outdir: Path,
    metadata: pd.DataFrame,
    minimap2: str,
    threads: int,
    skip: bool,
    precomputed_paf: Path | None,
    bootstrap: int,
    seed: int,
):
    paf = run_mapping(
        component_fasta,
        acre_fasta,
        outdir,
        minimap2,
        threads,
        skip,
        precomputed_paf,
    )
    all_hits = parse_family38_paf(paf, metadata)
    filtered = all_hits.loc[(all_hits.target_coverage >= 0.9) & (all_hits.identity >= 0.8)].copy()
    full_hits = deduplicate_mappings(filtered)
    full_hits.to_csv(outdir / "family38_full_hit_lengths.tsv.gz", sep="\t", index=False, compression="gzip")

    summary = full_hits.groupby(["species_id", "species_name"], sort=False).agg(
        full_hits=("acre_id", "size"),
        acres_with_full_hit=("acre_id", "nunique"),
        mean_query_span_bp=("query_span", "mean"),
        median_query_span_bp=("query_span", "median"),
        sd_query_span_bp=("query_span", "std"),
        min_query_span_bp=("query_span", "min"),
        max_query_span_bp=("query_span", "max"),
        mean_identity=("identity", "mean"),
    ).reset_index()
    summary["_order"] = summary.species_id.map({s: i for i, s in enumerate(SPECIES_ORDER)})
    summary = summary.sort_values("_order").drop(columns="_order")
    summary.to_csv(outdir / "family38_full_hit_summary.tsv", sep="\t", index=False)

    by_array = full_hits.groupby(["species_id", "species_name", "array_id"], sort=False).agg(
        full_hits=("acre_id", "size"), median_query_span_bp=("query_span", "median"), mean_query_span_bp=("query_span", "mean")
    ).reset_index()
    by_array["lineage"] = by_array.species_id.map(LINEAGE)
    by_array.to_csv(outdir / "family38_full_hits_per_array.tsv", sep="\t", index=False)
    comparisons = [
        ("Pan", ["mPanPan1", "mPanTro3"], "all non-Pan", ["hg002", "mGorGor1", "mPonAbe1", "mPonPyg2"]),
        ("Pan", ["mPanPan1", "mPanTro3"], "human", ["hg002"]),
        ("Pan", ["mPanPan1", "mPanTro3"], "gorilla", ["mGorGor1"]),
        ("Pan", ["mPanPan1", "mPanTro3"], "orangutans", ["mPonAbe1", "mPonPyg2"]),
    ]
    stats = compare_groups(by_array, "median_query_span_bp", "species_id", comparisons, bootstrap, seed + 1)
    arrays = [by_array.loc[by_array.species_id.eq(s), "median_query_span_bp"].to_numpy() for s in SPECIES_ORDER]
    kw = kruskal(*arrays)
    stats["six_species_kruskal_h"] = kw.statistic
    stats["six_species_kruskal_p"] = kw.pvalue
    return full_hits, summary, by_array, stats


def motif_analysis(
    component_fasta: Path,
    acre_fasta: Path,
    outdir: Path,
    metadata: pd.DataFrame,
    full_hits: pd.DataFrame,
):
    refs = read_fasta(component_fasta)
    acres_raw = read_fasta(acre_fasta)
    acres = {name.split("|")[0]: seq for name, seq in acres_raw.items()}

    reference_rows = []
    reference_site_rows = []
    for ref_name, seq in refs.items():
        for motif_name, motif in MOTIFS.items():
            hits = overlapping_matches(seq, motif["regex"])
            reference_rows.append(
                {
                    "reference": ref_name,
                    "reference_length": len(seq),
                    "is_family38": ref_name == FAMILY38,
                    "motif": motif_name,
                    "motif_display": motif["display"],
                    "match_count_both_strands": len(hits),
                }
            )
            for start, end, strand, matched in hits:
                reference_site_rows.append(
                    {
                        "reference": ref_name,
                        "is_family38": ref_name == FAMILY38,
                        "motif": motif_name,
                        "start_0based": start,
                        "end_0based_exclusive": end,
                        "strand": strand,
                        "matched_reference_sequence": matched,
                        "motif_oriented_sequence": matched if strand == "+" else revcomp(matched),
                    }
                )
    reference_scan = pd.DataFrame(reference_rows)
    reference_scan.to_csv(outdir / "family38_motif_reference_scan.tsv", sep="\t", index=False)
    pd.DataFrame(reference_site_rows).to_csv(outdir / "reference_motif_sites.tsv", sep="\t", index=False)

    site_rows = []
    hit_positive = []
    for hit_index, row in full_hits.reset_index(drop=True).iterrows():
        seq = acres[row.acre_id][int(row.query_start):int(row.query_end)]
        motif_counts = {}
        for motif_name, motif in MOTIFS.items():
            sites = overlapping_matches(seq, motif["regex"])
            motif_counts[motif_name] = len(sites)
            for start, end, strand, matched in sites:
                site_rows.append(
                    {
                        "hit_index": hit_index,
                        "acre_id": row.acre_id,
                        "array_id": row.array_id,
                        "species_id": row.species_id,
                        "motif": motif_name,
                        "hit_relative_start": start,
                        "hit_relative_end": end,
                        "acre_start": int(row.query_start) + start,
                        "acre_end": int(row.query_start) + end,
                        "strand_in_acre": strand,
                        "matched_sequence": matched,
                    }
                )
        hit_positive.append({"hit_index": hit_index, "acre_id": row.acre_id, "array_id": row.array_id, "species_id": row.species_id, **motif_counts})

    sites = pd.DataFrame(site_rows)
    sites.to_csv(outdir / "family38_motif_sites.tsv.gz", sep="\t", index=False, compression="gzip")
    positive = pd.DataFrame(hit_positive)
    positive.to_csv(outdir / "family38_motif_counts_per_full_hit.tsv.gz", sep="\t", index=False, compression="gzip")
    summary_rows = []
    for species in SPECIES_ORDER:
        group = positive.loc[positive.species_id.eq(species)]
        for motif_name, motif in MOTIFS.items():
            n_sites = int(group[motif_name].sum())
            n_positive = int((group[motif_name] > 0).sum())
            summary_rows.append(
                {
                    "species_id": species,
                    "species_name": metadata.loc[metadata.species_id.eq(species), "species_name"].iloc[0],
                    "motif": motif_name,
                    "motif_display": motif["display"],
                    "full_hits_scanned": len(group),
                    "motif_positive_hits": n_positive,
                    "motif_positive_fraction": n_positive / len(group) if len(group) else np.nan,
                    "total_motif_sites": n_sites,
                }
            )
    observed_summary = pd.DataFrame(summary_rows)
    observed_summary.to_csv(outdir / "family38_motif_observed_summary.tsv", sep="\t", index=False)
    return reference_scan, sites, positive, observed_summary


def create_figure(outdir: Path, per_array: pd.DataFrame, per_species: pd.DataFrame, by_array: pd.DataFrame, motif_summary: pd.DataFrame):
    sns.set_theme(style="whitegrid", context="talk")
    palette = [COLORS[s] for s in SPECIES_ORDER]
    labels = [SPECIES_LABELS[s] for s in SPECIES_ORDER]
    fig, axes = plt.subplots(2, 2, figsize=(17, 12), constrained_layout=True)

    sns.boxplot(data=per_array, x="species_id", y="family38_fraction", hue="species_id", order=SPECIES_ORDER, palette=palette, showfliers=False, legend=False, ax=axes[0, 0])
    sns.stripplot(data=per_array, x="species_id", y="family38_fraction", order=SPECIES_ORDER, color="black", alpha=0.25, size=2.5, ax=axes[0, 0])
    axes[0, 0].set(title="A. Family-38 abundance by array", xlabel="", ylabel="Fraction of annotated ACRE bp")
    axes[0, 0].set_xticks(np.arange(6), labels, rotation=20, ha="right")

    sns.boxplot(data=by_array, x="species_id", y="median_query_span_bp", hue="species_id", order=SPECIES_ORDER, palette=palette, showfliers=False, legend=False, ax=axes[0, 1])
    sns.stripplot(data=by_array, x="species_id", y="median_query_span_bp", order=SPECIES_ORDER, color="black", alpha=0.25, size=2.5, ax=axes[0, 1])
    axes[0, 1].axhline(1300, color="#555555", linestyle="--", linewidth=1.5, label="1,300-bp reference")
    axes[0, 1].legend(frameon=False, fontsize=11)
    axes[0, 1].set(title="B. Near-full-length family-38 homolog size", xlabel="", ylabel="Median aligned span per array (bp)")
    axes[0, 1].set_xticks(np.arange(6), labels, rotation=20, ha="right")

    plot_species = per_species.set_index("species_id").loc[SPECIES_ORDER]
    axes[1, 0].bar(np.arange(6), plot_species.pooled_family38_fraction, color=palette)
    axes[1, 0].set_xticks(np.arange(6), labels, rotation=20, ha="right")
    axes[1, 0].set(title="C. Pooled family-38 coverage", xlabel="", ylabel="Family-38 bp / complete ACRE bp")
    for i, (_, row) in enumerate(plot_species.iterrows()):
        axes[1, 0].text(i, row.pooled_family38_fraction + 0.006, f"{row.total_family38_bp/1e6:.2f} Mb", ha="center", va="bottom", fontsize=10)

    selected = ["human_A_exact_13mer", "human_A_broad_consensus", "Pan_common_predicted", "gorilla_predicted"]
    motif_plot = motif_summary.pivot(index="species_id", columns="motif", values="motif_positive_fraction").loc[SPECIES_ORDER, selected]
    x = np.arange(6)
    width = 0.19
    motif_colors = ["#D62728", "#1F77B4", "#2CA02C", "#9467BD"]
    motif_labels = ["Human A exact", "Human A broad", "Pan predicted", "Gorilla predicted"]
    for i, motif in enumerate(selected):
        axes[1, 1].bar(x + (i - 1.5) * width, motif_plot[motif], width, label=motif_labels[i], color=motif_colors[i])
    axes[1, 1].set_xticks(x, labels, rotation=20, ha="right")
    axes[1, 1].set_ylim(0, 1.08)
    axes[1, 1].set(title="D. Motifs in near-full-length homologs", xlabel="", ylabel="Fraction of homologs with ≥1 match")
    axes[1, 1].legend(frameon=False, fontsize=10, ncol=2)

    fig.suptitle("Family-38 abundance, homolog length, and candidate PRDM9 motifs", fontsize=22, weight="bold")
    fig.savefig(outdir / "family38_species_comparison.png", dpi=250, bbox_inches="tight")
    fig.savefig(outdir / "family38_species_comparison.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    root = args.root.expanduser().resolve()
    outdir = args.outdir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    metadata_path = (args.metadata or root / "primate_full/six_primate_acre_metadata.tsv").expanduser().resolve()
    annotations_path = (args.annotations or root / "primate_full/six_primate_full.annotations.bed").expanduser().resolve()
    acre_fasta = (args.acre_fasta or root / "primate_full/six_primate_full.acres.fasta").expanduser().resolve()
    component_fasta = (args.component_fasta or root / "references/dfam_2.fasta").expanduser().resolve()
    required = [metadata_path, annotations_path, acre_fasta, component_fasta]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError(
            "Missing required input(s): " + ", ".join(missing) +
            ". Supply --root/--data-root or the individual path options."
        )
    precomputed_paf = args.paf.expanduser().resolve() if args.paf else None
    metadata = pd.read_csv(metadata_path, sep="\t")

    family, per_acre, per_array, per_species, abundance_stats = annotation_analysis(
        annotations_path, outdir, metadata, args.bootstrap, args.seed
    )
    full_hits, hit_summary, by_array, length_stats = mapping_analysis(
        component_fasta,
        acre_fasta,
        outdir,
        metadata,
        args.minimap2,
        args.threads,
        args.skip_mapping,
        precomputed_paf,
        args.bootstrap,
        args.seed,
    )
    reference_scan, sites, positive, motif_summary = motif_analysis(
        component_fasta, acre_fasta, outdir, metadata, full_hits
    )
    create_figure(outdir, per_array, per_species, by_array, motif_summary)

    abundance_stats.to_csv(outdir / "family38_abundance_group_tests.tsv", sep="\t", index=False)
    length_stats.to_csv(outdir / "family38_length_group_tests.tsv", sep="\t", index=False)
    family38_ref = read_fasta(component_fasta)[FAMILY38]
    summary = {
        "family38_reference_length_bp": len(family38_ref),
        "complete_acres": int(len(metadata)),
        "complete_acre_bp": int(metadata.length.sum()),
        "family38_annotation_tracts": int(len(family)),
        "family38_annotation_bp": int(family.tract_bp.sum()),
        "near_full_length_hit_definition": "target coverage >= 0.90 and nucleotide identity >= 0.80, deduplicated at >=80% reciprocal query overlap",
        "near_full_length_family38_hits": int(len(full_hits)),
        "motif_matches_are_binding_evidence": False,
        "inputs": {
            "metadata": str(metadata_path),
            "annotations": str(annotations_path),
            "acre_fasta": str(acre_fasta),
            "component_references": str(component_fasta),
            "precomputed_paf": str(precomputed_paf) if precomputed_paf else None,
        },
    }
    (outdir / "analysis_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
