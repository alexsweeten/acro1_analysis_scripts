#!/usr/bin/env python3
"""Three divergence analyses for Supplementary Figure 2.

1. Array-balanced intra- versus inter-array divergence pooled across six
   great-ape assemblies, with inter-array comparisons restricted to the same
   species.
2. The same analysis displayed separately for every species.
3. Human HPRC split-array analysis comparing curated LHS blocks, curated RHS
   blocks, and arrays without a transition call.

Great-ape distances use exact, unsketched canonical 15-mer sets for complete
ACRO units. The HPRC analysis uses the same bounded 15-mer dissimilarity with
MinHash sketches: within-block divergence is averaged over unit pairs, whereas
between-block divergence is averaged over block-level sketches.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import math
import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial.distance import pdist, squareform
from scipy.stats import wilcoxon


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
FIGURES = ROOT / "figures"
TABLES = ROOT / "tables"
INTERMEDIATE = ROOT / "intermediate"

PRIMATE_ROOT = Path(os.environ.get("GREAT_APE_ACRE_ROOT", ROOT / "inputs" / "great_apes"))
PRIMATE_METADATA = PRIMATE_ROOT / "six_primate_acre_metadata.tsv"
PRIMATE_FASTA = PRIMATE_ROOT / "six_primate_full.acres.fasta"

WORKSPACE = Path(os.environ.get("ACRO1_WORKSPACE", ROOT / "inputs"))
HPRC_ROOT = WORKSPACE / "NTRIAGE_OCT26"
HPRC_MANIFEST = HPRC_ROOT / "manifest.validated.tsv"
HPRC_PIPE = HPRC_ROOT / "promoter_purifying_2794"
HPRC_ALL_RESULTS = HPRC_PIPE / "combined" / "all_array_results.tsv"
HPRC_SPLIT_METRICS = (
    HPRC_ROOT
    / "manuscript_output"
    / "split_side_divergence_selection"
    / "tables"
    / "array_side_metrics.tsv.gz"
)
HPRC_NTRIAGE = HPRC_ROOT / "manuscript_output" / "telomeric_ntriage_2794"
HPRC_CORE_FASTA = HPRC_NTRIAGE / "inputs" / "telomeric_2794.fasta"
HPRC_CLUSTERS = HPRC_NTRIAGE / "native" / "telomeric_2794_core.clusters.tsv"
HPRC_CURATED_AUDIT = (
    WORKSPACE
    / "hprc_breakpoint_ntriage_cluster_analysis"
    / "custom_manual_shift_breakpoint_map_round2"
    / "tables"
    / "round2_row_audit_1196.tsv"
)

KMER = 15
HPRC_UNIT_SKETCH = 512
HPRC_BLOCK_SKETCH = 2048
SEED = 20261029
N_SIGNFLIP = 100_000
BASES = {"A": 0, "C": 1, "G": 2, "T": 3}
BASE_ORDER = "ACGT"
INTRA_COLOR = "#1D395F"
INTER_COLOR = "#86CCF1"
HUMAN_COLORS = {
    "Split LHS": "#1D395F",
    "Split RHS": "#E8872D",
    "No split": "#86CCF1",
}
SPECIES_LABELS = {
    "human": "Human",
    "gorilla": "Gorilla",
    "chimpanzee": "Chimpanzee",
    "bonobo": "Bonobo",
    "sumatran_orangutan": "Sumatran\norangutan",
    "bornean_orangutan": "Bornean\norangutan",
}
SPECIES_ORDER = [
    "human",
    "gorilla",
    "chimpanzee",
    "bonobo",
    "sumatran_orangutan",
    "bornean_orangutan",
]


def configure_style() -> None:
    helvetica = Path("/System/Library/Fonts/Helvetica.ttc")
    if helvetica.exists():
        font_manager.fontManager.addfont(str(helvetica))
    mpl.rcParams.update(
        {
            "font.family": "Helvetica",
            "font.size": 15,
            "axes.titlesize": 22,
            "axes.labelsize": 18,
            "xtick.labelsize": 15,
            "ytick.labelsize": 14,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.linewidth": 1.3,
        }
    )


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    name: str | None = None
    pieces: list[str] = []
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(pieces).upper()
                name = line[1:].split()[0].split("|")[0]
                pieces = []
            else:
                pieces.append(line)
    if name is not None:
        records[name] = "".join(pieces).upper()
    return records


def read_selected_fasta(path: Path, selected: set[str]) -> dict[str, str]:
    records: dict[str, str] = {}
    name: str | None = None
    pieces: list[str] = []
    keep = False
    with open(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None and keep:
                    records[name] = "".join(pieces).upper()
                name = line[1:].split()[0].split("|")[0]
                keep = name in selected
                pieces = []
            elif keep:
                pieces.append(line)
    if name is not None and keep:
        records[name] = "".join(pieces).upper()
    return records


def distinct_canonical_kmers(sequence: str, k: int = KMER) -> np.ndarray:
    mask = (1 << (2 * k)) - 1
    shift = 2 * (k - 1)
    forward = reverse = valid = 0
    values: set[int] = set()
    for base in sequence:
        code = BASES.get(base)
        if code is None:
            forward = reverse = valid = 0
            continue
        forward = ((forward << 2) | code) & mask
        reverse = (reverse >> 2) | ((3 - code) << shift)
        valid += 1
        if valid >= k:
            values.add(min(forward, reverse))
    return np.fromiter(values, dtype=np.uint64, count=len(values))


def splitmix64(values: np.ndarray) -> np.ndarray:
    """Apply a deterministic bijective 64-bit mixer to encoded k-mers."""
    values = values.astype(np.uint64, copy=True)
    values += np.uint64(0x9E3779B97F4A7C15)
    values = (values ^ (values >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
    values = (values ^ (values >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
    return values ^ (values >> np.uint64(31))


def bottom_hash_sketch(sequence: str, size: int) -> np.ndarray:
    codes = distinct_canonical_kmers(sequence)
    if len(codes) == 0:
        return np.empty(0, dtype=np.uint64)
    hashed = splitmix64(codes)
    if len(hashed) > size:
        hashed = np.partition(hashed, size - 1)[:size]
    return np.sort(hashed)


def sketch_dissimilarity(left: np.ndarray, right: np.ndarray) -> float:
    denominator = min(len(left), len(right))
    if denominator == 0:
        return math.nan
    shared = np.intersect1d(left, right, assume_unique=True).size
    return 1.0 - shared / denominator


def iter_fasta(path: Path):
    name: str | None = None
    pieces: list[str] = []
    with open(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(pieces).upper()
                name = line[1:].split()[0]
                pieces = []
            else:
                pieces.append(line)
    if name is not None:
        yield name, "".join(pieces).upper()


def exact_kmer_distance_matrix(kmer_sets: list[np.ndarray], k: int = KMER) -> np.ndarray:
    sizes = np.asarray([len(values) for values in kmer_sets], dtype=np.int32)
    if np.any(sizes == 0):
        raise ValueError("A selected sequence contained no valid canonical k-mers")
    all_codes = np.concatenate(kmer_sets)
    rows = np.repeat(np.arange(len(kmer_sets), dtype=np.int32), sizes)
    _, columns = np.unique(all_codes, return_inverse=True)
    incidence = sparse.csr_matrix(
        (np.ones(len(rows), dtype=np.int16), (rows, columns)),
        shape=(len(kmer_sets), int(columns.max()) + 1),
    )
    shared = (incidence @ incidence.T).toarray().astype(np.int32, copy=False)
    union = sizes[:, None] + sizes[None, :] - shared
    jaccard = shared / union
    # A bounded dissimilarity is used because reference-poor ACRO classes can
    # have zero exact k-mer overlap, where Mash distance is saturated and no
    # longer interpretable as a nucleotide percentage.
    distance = 1.0 - jaccard
    np.fill_diagonal(distance, 0.0)
    return distance


def array_level_from_distance(frame: pd.DataFrame, distance: np.ndarray) -> pd.DataFrame:
    groups = {
        array_id: group.index.to_numpy(dtype=np.int32)
        for array_id, group in frame.groupby("array_id", sort=True)
    }
    rows: list[dict[str, object]] = []
    array_ids = sorted(groups)
    for array_id in array_ids:
        indices = groups[array_id]
        left, right = np.triu_indices(len(indices), 1)
        intra_values = 100.0 * distance[indices[left], indices[right]]
        between_array_means = []
        between_pair_count = 0
        for other_id in array_ids:
            if other_id == array_id:
                continue
            other = groups[other_id]
            block = 100.0 * distance[np.ix_(indices, other)]
            between_array_means.append(float(np.mean(block)))
            between_pair_count += int(block.size)
        source = frame.loc[indices]
        rows.append(
            {
                "species": source.species_name.iloc[0],
                "assembly": source.assembly.iloc[0],
                "array_id": array_id,
                "chromosome": source.original_chrom.iloc[0],
                "complete_units": len(indices),
                "intra_unit_pairs": len(intra_values),
                "inter_unit_pairs": between_pair_count,
                "intra_array_mean_divergence_percent": float(np.mean(intra_values)),
                "intra_array_median_divergence_percent": float(np.median(intra_values)),
                "inter_array_mean_divergence_percent": float(np.mean(between_array_means)),
                "inter_array_median_of_array_pair_means_percent": float(np.median(between_array_means)),
            }
        )
    return pd.DataFrame(rows)


def p_label(value: float) -> str:
    if value < 1e-4:
        return "$P$ < 10$^{-4}$"
    if value < 0.001:
        return f"$P$ = {value:.1e}"
    return f"$P$ = {value:.3f}"


def summarize_great_apes(array_level: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    summaries = []
    tests = []
    cohorts = [("All great apes", array_level)] + [
        (SPECIES_LABELS[species].replace("\n", " "), array_level[array_level.species == species])
        for species in SPECIES_ORDER
    ]
    for label, group in cohorts:
        intra = group.intra_array_mean_divergence_percent.to_numpy(float)
        inter = group.inter_array_mean_divergence_percent.to_numpy(float)
        for metric, values in (("Intra-array", intra), ("Inter-array", inter)):
            summaries.append(
                {
                    "cohort": label,
                    "metric": metric,
                    "arrays": len(values),
                    "mean_divergence_percent": np.mean(values),
                    "median_divergence_percent": np.median(values),
                    "q25_divergence_percent": np.quantile(values, 0.25),
                    "q75_divergence_percent": np.quantile(values, 0.75),
                }
            )
        test = wilcoxon(intra, inter, alternative="less", zero_method="wilcox")
        tests.append(
            {
                "cohort": label,
                "arrays": len(group),
                "mean_intra_minus_inter_percentage_points": float(np.mean(intra - inter)),
                "median_intra_minus_inter_percentage_points": float(np.median(intra - inter)),
                "paired_wilcoxon_statistic": float(test.statistic),
                "paired_wilcoxon_one_sided_p": float(test.pvalue),
                "alternative": "intra-array divergence < species-matched inter-array divergence",
            }
        )
    return pd.DataFrame(summaries), pd.DataFrame(tests)


def styled_boxplot(ax, data: list[np.ndarray], positions: list[float], colors: list[str], widths: float = 0.58) -> None:
    artists = ax.boxplot(
        data,
        positions=positions,
        widths=widths,
        patch_artist=True,
        showfliers=False,
        whis=1.5,
        medianprops={"color": "black", "linewidth": 2.2},
        boxprops={"edgecolor": "#333333", "linewidth": 1.4},
        whiskerprops={"color": "#555555", "linewidth": 1.4},
        capprops={"color": "#555555", "linewidth": 1.4},
    )
    for patch, color in zip(artists["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.92)


def save_figure(fig: plt.Figure, stem: str) -> None:
    for suffix, kwargs in (
        ("png", {"dpi": 300}),
        ("pdf", {}),
        ("svg", {}),
    ):
        fig.savefig(FIGURES / f"{stem}.{suffix}", bbox_inches="tight", facecolor="white", **kwargs)


def make_great_ape_figures(array_level: pd.DataFrame, tests: pd.DataFrame) -> None:
    configure_style()
    intra = array_level.intra_array_mean_divergence_percent.to_numpy()
    inter = array_level.inter_array_mean_divergence_percent.to_numpy()
    pooled_p = float(tests.loc[tests.cohort == "All great apes", "paired_wilcoxon_one_sided_p"].iloc[0])

    fig, ax = plt.subplots(figsize=(8.4, 7.2))
    styled_boxplot(ax, [intra, inter], [1, 2], [INTRA_COLOR, INTER_COLOR], 0.56)
    rng = np.random.default_rng(SEED)
    for pos, values, color in ((1, intra, INTRA_COLOR), (2, inter, INTER_COLOR)):
        ax.scatter(pos + rng.uniform(-0.12, 0.12, len(values)), values, s=11, color=color, alpha=0.20, linewidths=0)
    ax.set_xticks([1, 2], ["Intra-array", "Inter-array"])
    ax.set_ylabel("Exact canonical 15-mer dissimilarity (%)")
    ax.set_title("All great apes", fontweight="bold", loc="left")
    ax.text(0.02, 0.98, f"n = {len(array_level)} arrays; species-matched inter-array comparisons\n{p_label(pooled_p)}",
            transform=ax.transAxes, ha="left", va="top", color="#4F5D6B", fontsize=14)
    ax.grid(axis="y", color="#D9DDE2", linewidth=0.9)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("ACRO units are more similar within arrays than between arrays", fontsize=23, fontweight="bold", y=1.01)
    fig.tight_layout()
    save_figure(fig, "Supplementary_Figure_2A_all_great_apes_pooled")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(15.5, 7.5))
    positions = np.arange(len(SPECIES_ORDER), dtype=float) * 2.0 + 1.0
    plot_data: list[np.ndarray] = []
    plot_pos: list[float] = []
    colors: list[str] = []
    for center, species in zip(positions, SPECIES_ORDER):
        group = array_level[array_level.species == species]
        plot_data.extend(
            [
                group.intra_array_mean_divergence_percent.to_numpy(),
                group.inter_array_mean_divergence_percent.to_numpy(),
            ]
        )
        plot_pos.extend([center - 0.34, center + 0.34])
        colors.extend([INTRA_COLOR, INTER_COLOR])
    styled_boxplot(ax, plot_data, plot_pos, colors, 0.57)
    tick_labels = [
        f"{SPECIES_LABELS[species]}\nn={int((array_level.species == species).sum())}"
        for species in SPECIES_ORDER
    ]
    ax.set_xticks(positions, tick_labels)
    ax.set_ylabel("Exact canonical 15-mer dissimilarity (%)")
    ax.set_title("Great-ape assemblies analyzed separately", fontweight="bold", loc="left")
    handles = [
        mpl.patches.Patch(facecolor=INTRA_COLOR, edgecolor="#333333", label="Intra-array"),
        mpl.patches.Patch(facecolor=INTER_COLOR, edgecolor="#333333", label="Inter-array"),
    ]
    ax.legend(handles=handles, frameon=False, ncol=2, loc="upper left", fontsize=15)
    ax.grid(axis="y", color="#D9DDE2", linewidth=0.9)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    save_figure(fig, "Supplementary_Figure_2B_great_apes_by_species")
    plt.close(fig)


def encode_alignment(sequences: list[str], columns: np.ndarray | None = None) -> np.ndarray:
    length = len(sequences[0])
    if any(len(sequence) != length for sequence in sequences):
        raise ValueError("Alignment records have unequal lengths")
    if columns is None:
        columns = np.arange(length)
    width = int(np.sum(columns)) if np.asarray(columns).dtype == bool else len(columns)
    output = np.full((len(sequences), width), -1, dtype=np.int8)
    for row, sequence in enumerate(sequences):
        selected = np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)[columns]
        for base, code in BASES.items():
            output[row, selected == ord(base)] = code
    return output


def consensus_from_alignment(sequences: list[str]) -> np.ndarray:
    encoded = encode_alignment(sequences)
    output = np.full(encoded.shape[1], -1, dtype=np.int8)
    for base in range(4):
        counts = np.sum(encoded == base, axis=0)
        replace = counts > np.max(np.stack([np.sum(encoded == prior, axis=0) for prior in range(4)]), axis=0)
        output[replace] = base
    # Deterministic tie handling: np.argmax selects A, C, G, T in that order.
    count_matrix = np.stack([np.sum(encoded == base, axis=0) for base in range(4)])
    callable_columns = np.sum(count_matrix, axis=0) > 0
    output[callable_columns] = np.argmax(count_matrix[:, callable_columns], axis=0).astype(np.int8)
    return output


def pooled_nucleotide_divergence(sequences: list[str], columns: np.ndarray) -> tuple[float, int, int]:
    encoded = encode_alignment(sequences, columns)
    valid = np.sum(encoded >= 0, axis=0).astype(np.int64)
    pair_count = valid * (valid - 1) // 2
    same = np.zeros(encoded.shape[1], dtype=np.int64)
    for base in range(4):
        counts = np.sum(encoded == base, axis=0).astype(np.int64)
        same += counts * (counts - 1) // 2
    differences = pair_count - same
    usable = pair_count > 0
    denominator = int(np.sum(pair_count[usable]))
    numerator = int(np.sum(differences[usable]))
    return (numerator / denominator if denominator else math.nan, numerator, denominator)


def load_hprc_block(
    safe_id: str,
    source_name: str,
    individual_id: str,
    category: str,
    split_local_bp: float | None,
) -> dict[str, object] | None:
    directory = HPRC_PIPE / "per_array" / safe_id
    alignment_path = directory / "backbone" / "units_reference_mapped.fa"
    qc_path = directory / "units_qc.tsv"
    if not alignment_path.exists() or not qc_path.exists():
        return None
    aligned = read_fasta(alignment_path)
    qc = pd.read_csv(qc_path, sep="\t")
    qc = qc[qc.included_in_structural_msa.astype(str).str.lower().isin(["true", "1"])].copy()
    qc = qc[qc.unit_id.isin(aligned)].copy()
    if split_local_bp is not None:
        midpoint = (qc.local_start.astype(float) + qc.local_end.astype(float)) / 2.0
        if category == "Split LHS":
            qc = qc[midpoint <= split_local_bp]
        elif category == "Split RHS":
            qc = qc[midpoint > split_local_bp]
    qc = qc.sort_values(["local_start", "local_end", "unit_id"])
    if len(qc) < 3:
        return None
    sequences = [aligned[unit_id] for unit_id in qc.unit_id]
    if len({len(sequence) for sequence in sequences}) != 1:
        return None
    consensus = consensus_from_alignment(sequences)
    return {
        "source_name": source_name,
        "safe_id": safe_id,
        "individual_id": individual_id,
        "category": category,
        "unit_ids": qc.unit_id.tolist(),
        "unit_count": len(qc),
        "alignment_path": str(alignment_path),
        "consensus": consensus,
    }


def build_hprc_blocks() -> list[dict[str, object]]:
    split = pd.read_csv(HPRC_SPLIT_METRICS, sep="\t")
    split = split[
        (split.status == "valid")
        & split.round2_transition.isin(["C1->C2", "C6->C2"])
    ].copy()
    blocks: list[dict[str, object]] = []
    for row in split.itertuples(index=False):
        left = load_hprc_block(
            row.safe_id, row.source_name, row.individual_id, "Split LHS", float(row.split_local_bp)
        )
        right = load_hprc_block(
            row.safe_id, row.source_name, row.individual_id, "Split RHS", float(row.split_local_bp)
        )
        if left is not None and right is not None:
            blocks.extend([left, right])

    manifest = pd.read_csv(HPRC_MANIFEST, sep="\t")
    results = pd.read_csv(HPRC_ALL_RESULTS, sep="\t", low_memory=False)
    eligible = manifest[manifest.has_recombination == 0].merge(
        results[["array_id", "bb_status", "mixed_strand"]], on="array_id", how="left"
    )
    eligible = eligible[(eligible.bb_status == "estimable") & (~eligible.mixed_strand.fillna(True))]
    for row in eligible.itertuples(index=False):
        block = load_hprc_block(
            row.safe_id, row.array_id, row.individual_id, "No split", None
        )
        if block is not None:
            blocks.append(block)
    return blocks


def stable_seed(label: str) -> int:
    return int.from_bytes(hashlib.sha256(f"{SEED}|{label}".encode()).digest()[:8], "little")


def signflip_p(values: np.ndarray, alternative: str, label: str) -> float:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    observed = float(np.mean(values))
    rng = np.random.default_rng(stable_seed(label))
    extreme = 0
    batch = 5000
    for start in range(0, N_SIGNFLIP, batch):
        count = min(batch, N_SIGNFLIP - start)
        means = np.mean(rng.choice(np.array([-1.0, 1.0]), size=(count, len(values))) * values, axis=1)
        if alternative == "less":
            extreme += int(np.sum(means <= observed + 1e-15))
        elif alternative == "greater":
            extreme += int(np.sum(means >= observed - 1e-15))
        else:
            extreme += int(np.sum(np.abs(means) >= abs(observed) - 1e-15))
    return (extreme + 1) / (N_SIGNFLIP + 1)


def bh_adjust(values: pd.Series) -> pd.Series:
    raw = values.to_numpy(float)
    order = np.argsort(raw)
    ranked = raw[order] * len(raw) / np.arange(1, len(raw) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    output = np.empty_like(ranked)
    output[order] = np.minimum(ranked, 1.0)
    return pd.Series(output, index=values.index)


def paired_individual_test(frame: pd.DataFrame, metric: str, left: str, right: str) -> dict[str, object]:
    individual = frame.groupby(["individual_id", "category"])[metric].mean().unstack()
    paired = individual[[left, right]].dropna()
    differences = paired[left].to_numpy() - paired[right].to_numpy()
    return {
        "metric": metric,
        "left_group": left,
        "right_group": right,
        "n_individuals": len(paired),
        "left_individual_mean_percent": float(np.mean(paired[left])),
        "right_individual_mean_percent": float(np.mean(paired[right])),
        "mean_left_minus_right_percentage_points": float(np.mean(differences)),
        "median_left_minus_right_percentage_points": float(np.median(differences)),
        "two_sided_individual_signflip_p": signflip_p(
            differences, "two-sided", f"{metric}|{left}|{right}"
        ),
    }


def unit_count_matched_sensitivity(frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Match each split block to the closest-size no-split array from the same individual."""
    matches = []
    no_split = frame[frame.category == "No split"].copy()
    for category in ("Split LHS", "Split RHS"):
        for row in frame[frame.category == category].itertuples(index=False):
            candidates = no_split[no_split.individual_id == row.individual_id].copy()
            if candidates.empty:
                continue
            candidates["unit_count_difference"] = (candidates.unit_count - row.unit_count).abs()
            match = candidates.sort_values(["unit_count_difference", "source_name"]).iloc[0]
            matches.append(
                {
                    "split_category": category,
                    "individual_id": row.individual_id,
                    "split_source_name": row.source_name,
                    "no_split_source_name": match.source_name,
                    "split_units": row.unit_count,
                    "no_split_units": int(match.unit_count),
                    "unit_count_difference": int(match.unit_count - row.unit_count),
                    "intra_split_percent": row.intra_array_divergence_percent,
                    "intra_no_split_percent": match.intra_array_divergence_percent,
                    "intra_split_minus_no_split_percentage_points": (
                        row.intra_array_divergence_percent - match.intra_array_divergence_percent
                    ),
                    "inter_split_percent": row.inter_array_divergence_percent,
                    "inter_no_split_percent": match.inter_array_divergence_percent,
                    "inter_split_minus_no_split_percentage_points": (
                        row.inter_array_divergence_percent - match.inter_array_divergence_percent
                    ),
                }
            )
    match_frame = pd.DataFrame(matches)
    summaries = []
    for category, group in match_frame.groupby("split_category"):
        for metric, column in (
            ("intra_array_divergence_percent", "intra_split_minus_no_split_percentage_points"),
            ("inter_array_divergence_percent", "inter_split_minus_no_split_percentage_points"),
        ):
            individual = group.groupby("individual_id")[column].mean()
            summaries.append(
                {
                    "split_category": category,
                    "metric": metric,
                    "matched_split_blocks": len(group),
                    "individuals": len(individual),
                    "median_absolute_unit_count_difference": float(np.median(np.abs(group.unit_count_difference))),
                    "mean_split_minus_matched_no_split_percentage_points": float(group[column].mean()),
                    "individual_mean_split_minus_matched_no_split_percentage_points": float(individual.mean()),
                    "two_sided_individual_signflip_p": signflip_p(
                        individual.to_numpy(), "two-sided", f"matched|{category}|{metric}"
                    ),
                }
            )
    summary = pd.DataFrame(summaries)
    summary["bh_q_across_four_tests"] = bh_adjust(summary.two_sided_individual_signflip_p)
    return match_frame, summary


def block_sketch_from_unit_sketches(sketches: list[np.ndarray], size: int) -> np.ndarray:
    if not sketches:
        return np.empty(0, dtype=np.uint64)
    candidates = np.unique(np.concatenate(sketches))
    if len(candidates) > size:
        candidates = np.partition(candidates, size - 1)[:size]
    return np.sort(candidates)


def block_intra_dissimilarity(unit_sketches: list[np.ndarray]) -> tuple[float, int]:
    values = []
    for left in range(len(unit_sketches)):
        for right in range(left + 1, len(unit_sketches)):
            value = sketch_dissimilarity(unit_sketches[left], unit_sketches[right])
            if np.isfinite(value):
                values.append(value)
    return (float(np.mean(values)) if values else math.nan, len(values))


def block_inter_dissimilarity_matrix(sketches: list[np.ndarray]) -> np.ndarray:
    sizes = np.asarray([len(sketch) for sketch in sketches], dtype=np.int32)
    all_hashes = np.concatenate(sketches)
    rows = np.repeat(np.arange(len(sketches), dtype=np.int32), sizes)
    _, columns = np.unique(all_hashes, return_inverse=True)
    incidence = sparse.csr_matrix(
        (np.ones(len(rows), dtype=np.int16), (rows, columns)),
        shape=(len(sketches), int(columns.max()) + 1),
    )
    shared = (incidence @ incidence.T).toarray().astype(np.int32, copy=False)
    denominator = np.minimum(sizes[:, None], sizes[None, :])
    jaccard_estimate = shared / denominator
    distance = 1.0 - jaccard_estimate
    np.fill_diagonal(distance, 0.0)
    return distance


def analyze_hprc_blocks() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, int]:
    clusters = pd.read_csv(HPRC_CLUSTERS, sep="\t")
    clusters = clusters.sort_values(["source_name", "unit_in_array"])
    cluster_groups = {name: group.copy() for name, group in clusters.groupby("source_name", sort=False)}

    audit = pd.read_csv(HPRC_CURATED_AUDIT, sep="\t")
    audit = audit[audit.round2_transition.isin(["C1->C2", "C6->C2"])].copy()
    split_targets: dict[str, dict[str, object]] = {}
    for row in audit.itertuples(index=False):
        if row.source_name not in cluster_groups:
            continue
        total = len(cluster_groups[row.source_name])
        left = int(row.round2_units_left)
        if left >= 3 and total - left >= 3:
            split_targets[row.source_name] = {
                "individual_id": row.source_name.split("#", 1)[0],
                "left_units": left,
                "transition": row.round2_transition,
            }

    manifest = pd.read_csv(HPRC_MANIFEST, sep="\t")
    no_split = manifest[manifest.has_recombination == 0].copy()
    no_split_targets = {
        row.array_id: {"individual_id": row.individual_id}
        for row in no_split.itertuples(index=False)
        if row.array_id in cluster_groups and len(cluster_groups[row.array_id]) >= 3
    }
    target_names = set(split_targets) | set(no_split_targets)

    block_rows: list[dict[str, object]] = []
    block_sketches: list[np.ndarray] = []
    processed_targets = 0
    for source_name, sequence in iter_fasta(HPRC_CORE_FASTA):
        if source_name not in target_names:
            continue
        processed_targets += 1
        if processed_targets % 250 == 0:
            print(f"HPRC targets processed: {processed_targets}/{len(target_names)}", flush=True)
        units = cluster_groups[source_name]
        unit_large_sketches: list[np.ndarray] = []
        unit_small_sketches: list[np.ndarray] = []
        for row in units.itertuples(index=False):
            unit_sequence = sequence[int(row.start) : int(row.end)]
            large = bottom_hash_sketch(unit_sequence, HPRC_BLOCK_SKETCH)
            unit_large_sketches.append(large)
            unit_small_sketches.append(large[: min(HPRC_UNIT_SKETCH, len(large))])

        specifications: list[tuple[str, list[int], str, str]] = []
        if source_name in split_targets:
            target = split_targets[source_name]
            left_count = int(target["left_units"])
            specifications.extend(
                [
                    ("Split LHS", list(range(left_count)), str(target["individual_id"]), str(target["transition"])),
                    ("Split RHS", list(range(left_count, len(units))), str(target["individual_id"]), str(target["transition"])),
                ]
            )
        elif source_name in no_split_targets:
            target = no_split_targets[source_name]
            specifications.append(("No split", list(range(len(units))), str(target["individual_id"]), "none"))

        for category, selected, individual_id, transition in specifications:
            selected_small = [unit_small_sketches[index] for index in selected]
            selected_large = [unit_large_sketches[index] for index in selected]
            intra, pair_count = block_intra_dissimilarity(selected_small)
            block_sketch = block_sketch_from_unit_sketches(selected_large, HPRC_BLOCK_SKETCH)
            block_rows.append(
                {
                    "source_name": source_name,
                    "individual_id": individual_id,
                    "category": category,
                    "transition": transition,
                    "unit_count": len(selected),
                    "intra_unit_pairs": pair_count,
                    "intra_array_divergence_percent": 100.0 * intra,
                    "unit_sketch_size": HPRC_UNIT_SKETCH,
                    "block_sketch_size": len(block_sketch),
                }
            )
            block_sketches.append(block_sketch)

    frame = pd.DataFrame(block_rows)
    frame["inter_array_divergence_percent"] = np.nan
    for category, group in frame.groupby("category", sort=False):
        indices = group.index.to_numpy()
        distances = block_inter_dissimilarity_matrix([block_sketches[index] for index in indices])
        means = distances.sum(axis=1) / (len(indices) - 1)
        frame.loc[indices, "inter_array_divergence_percent"] = 100.0 * means

    summary_rows = []
    for category, group in frame.groupby("category", sort=False):
        for metric in ("intra_array_divergence_percent", "inter_array_divergence_percent"):
            values = group[metric].dropna().to_numpy()
            summary_rows.append(
                {
                    "category": category,
                    "metric": metric,
                    "blocks": len(values),
                    "individuals": group.individual_id.nunique(),
                    "mean_percent": float(np.mean(values)),
                    "median_percent": float(np.median(values)),
                    "q25_percent": float(np.quantile(values, 0.25)),
                    "q75_percent": float(np.quantile(values, 0.75)),
                }
            )
    summary = pd.DataFrame(summary_rows)

    test_rows = []
    for metric in ("intra_array_divergence_percent", "inter_array_divergence_percent"):
        for left, right in (
            ("Split LHS", "Split RHS"),
            ("Split RHS", "No split"),
            ("Split LHS", "No split"),
        ):
            test_rows.append(paired_individual_test(frame, metric, left, right))
    tests = pd.DataFrame(test_rows)
    tests["bh_q_across_six_tests"] = bh_adjust(tests.two_sided_individual_signflip_p)
    return frame, summary, tests, HPRC_BLOCK_SKETCH


def make_hprc_figure(frame: pd.DataFrame, tests: pd.DataFrame, block_sketch_size: int) -> None:
    configure_style()
    categories = ["Split LHS", "Split RHS", "No split"]
    metrics = [
        ("intra_array_divergence_percent", "Intra-array divergence", "Mean 15-mer dissimilarity among units within each block"),
        ("inter_array_divergence_percent", "Inter-array divergence", "Mean 15-mer dissimilarity to other blocks in the same group"),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(15.8, 7.5))
    for panel, (ax, (metric, title, subtitle)) in enumerate(zip(axes, metrics), start=1):
        values = [frame.loc[frame.category == category, metric].dropna().to_numpy() for category in categories]
        styled_boxplot(ax, values, [1, 2, 3], [HUMAN_COLORS[c] for c in categories], 0.58)
        rng = np.random.default_rng(SEED + panel)
        for position, category_values, category in zip([1, 2, 3], values, categories):
            take = category_values
            if len(take) > 700:
                take = rng.choice(take, 700, replace=False)
            ax.scatter(position + rng.uniform(-0.11, 0.11, len(take)), take, s=10,
                       color=HUMAN_COLORS[category], alpha=0.18, linewidths=0)
        category_labels = [
            f"{category}\nn={int((frame.category == category).sum())}"
            for category in categories
        ]
        ax.set_xticks([1, 2, 3], category_labels)
        ax.set_ylabel("Divergence (%)")
        ax.set_title(f"{'AB'[panel-1]}  {title}", fontweight="bold", loc="left")
        ax.text(0.02, 0.96, subtitle, transform=ax.transAxes, ha="left", va="top",
                fontsize=12.5, color="#4F5D6B")
        ax.grid(axis="y", color="#D9DDE2", linewidth=0.9)
        ax.set_axisbelow(True)
        ax.spines[["top", "right"]].set_visible(False)

    lhs_rhs = tests[(tests.left_group == "Split LHS") & (tests.right_group == "Split RHS")]
    for ax, metric in zip(axes, [item[0] for item in metrics]):
        row = lhs_rhs[lhs_rhs.metric == metric].iloc[0]
        ax.text(0.98, 0.96, f"LHS vs RHS: {p_label(float(row.two_sided_individual_signflip_p))}",
                transform=ax.transAxes, ha="right", va="top", fontsize=12.5,
                bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "edgecolor": "#B9C2CC"})
    fig.suptitle("Human ACRO divergence differs across split-array sequence blocks", fontsize=23,
                 fontweight="bold", y=1.01)
    fig.text(0.5, 0.005,
             f"Canonical 15-mer MinHash: {HPRC_UNIT_SKETCH} hashes per unit and {block_sketch_size:,} per block; statistics are individual-balanced.",
             ha="center", va="bottom", fontsize=12.5, color="#4F5D6B")
    fig.tight_layout(rect=(0, 0.035, 1, 0.98))
    save_figure(fig, "Supplementary_Figure_2C_human_split_block_comparison")
    plt.close(fig)


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    TABLES.mkdir(parents=True, exist_ok=True)
    INTERMEDIATE.mkdir(parents=True, exist_ok=True)

    metadata = pd.read_csv(PRIMATE_METADATA, sep="\t")
    sizes = metadata.groupby(["species_name", "array_id"]).size().rename("array_units")
    metadata = metadata.join(sizes, on=["species_name", "array_id"])
    metadata = metadata[metadata.array_units >= 2].copy()
    sequences = read_selected_fasta(PRIMATE_FASTA, set(metadata.acre_id))
    great_ape_parts = []
    excluded_units: list[dict[str, object]] = []
    for species in SPECIES_ORDER:
        print(f"Great-ape analysis: {species}", flush=True)
        frame = metadata[metadata.species_name == species].reset_index(drop=True)
        kmer_sets = [distinct_canonical_kmers(sequences[acre_id]) for acre_id in frame.acre_id]
        valid = np.asarray([len(values) > 0 for values in kmer_sets])
        if not np.all(valid):
            for row in frame.loc[~valid].itertuples(index=False):
                excluded_units.append(
                    {
                        "acre_id": row.acre_id,
                        "species": species,
                        "array_id": row.array_id,
                        "reason": "no_valid_canonical_15mers",
                    }
                )
            frame = frame.loc[valid].reset_index(drop=True)
            kmer_sets = [values for values, keep in zip(kmer_sets, valid) if keep]
            retained_sizes = frame.groupby("array_id").size()
            keep_arrays = retained_sizes[retained_sizes >= 2].index
            keep = frame.array_id.isin(keep_arrays).to_numpy()
            frame = frame.loc[keep].reset_index(drop=True)
            kmer_sets = [values for values, retain in zip(kmer_sets, keep) if retain]
        distance = exact_kmer_distance_matrix(kmer_sets)
        great_ape_parts.append(array_level_from_distance(frame, distance))
    great_ape = pd.concat(great_ape_parts, ignore_index=True)
    great_summary, great_tests = summarize_great_apes(great_ape)
    great_ape.to_csv(TABLES / "great_ape_array_level_divergence.tsv", sep="\t", index=False)
    great_summary.to_csv(TABLES / "great_ape_divergence_summary.tsv", sep="\t", index=False)
    great_tests.to_csv(TABLES / "great_ape_paired_tests.tsv", sep="\t", index=False)
    pd.DataFrame(excluded_units).to_csv(TABLES / "great_ape_excluded_units.tsv", sep="\t", index=False)
    make_great_ape_figures(great_ape, great_tests)

    human, human_summary, human_tests, block_sketch_size = analyze_hprc_blocks()
    human_matches, human_matched_tests = unit_count_matched_sensitivity(human)
    human.to_csv(TABLES / "human_split_block_array_level_divergence.tsv.gz", sep="\t", index=False, compression="gzip")
    human_summary.to_csv(TABLES / "human_split_block_divergence_summary.tsv", sep="\t", index=False)
    human_tests.to_csv(TABLES / "human_split_block_individual_tests.tsv", sep="\t", index=False)
    human_matches.to_csv(TABLES / "human_split_block_unit_count_matches.tsv.gz", sep="\t", index=False, compression="gzip")
    human_matched_tests.to_csv(TABLES / "human_split_block_unit_count_matched_tests.tsv", sep="\t", index=False)
    make_hprc_figure(human, human_tests, block_sketch_size)

    parameters = {
        "great_ape_source_metadata": str(PRIMATE_METADATA),
        "great_ape_source_fasta": str(PRIMATE_FASTA),
        "great_ape_kmer_size": KMER,
        "great_ape_kmer_method": "exact unsketched canonical k-mer Jaccard dissimilarity",
        "great_ape_inter_array_scope": "same species only; equal weight to every other array",
        "human_manifest": str(HPRC_MANIFEST),
        "human_curated_split_source": str(HPRC_SPLIT_METRICS),
        "human_split_transitions": ["C1->C2", "C6->C2"],
        "human_minimum_units_per_block": 3,
        "human_unit_sketch_size": HPRC_UNIT_SKETCH,
        "human_block_sketch_size": block_sketch_size,
        "human_intra_metric": "mean canonical 15-mer MinHash dissimilarity among complete units within each block",
        "human_inter_metric": "mean canonical 15-mer MinHash dissimilarity between a block and all other blocks in the same category",
        "signflip_permutations": N_SIGNFLIP,
        "seed": SEED,
    }
    (TABLES / "expanded_analysis_parameters.json").write_text(json.dumps(parameters, indent=2) + "\n")

    result = {
        "great_ape_arrays": int(len(great_ape)),
        "great_ape_complete_units": int(great_ape.complete_units.sum()),
        "human_blocks": int(len(human)),
        "human_block_sketch_size": block_sketch_size,
        "great_ape_summary": great_summary.to_dict(orient="records"),
        "great_ape_tests": great_tests.to_dict(orient="records"),
        "human_summary": human_summary.to_dict(orient="records"),
        "human_tests": human_tests.to_dict(orient="records"),
        "human_unit_count_matched_tests": human_matched_tests.to_dict(orient="records"),
    }
    (TABLES / "expanded_analysis_results.json").write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
