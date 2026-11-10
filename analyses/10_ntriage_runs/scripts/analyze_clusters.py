#!/usr/bin/env python3
"""Characterize the selected NTRIAGE clusters for every Nov-3 cohort.

The clustering features are used only to *define* clusters.  This script adds
array-aware descriptive summaries of unit length, ACRO1 allele/copy content,
and annotation-derived backbone content.  Tests are deliberately labelled
descriptive because tandem units from the same parent array are not
independent and some of the annotations also contributed to clustering.

Expected layout (created by the companion scripts)::

    inputs/dataset_manifest.tsv
    runs/<dataset>/native/<dataset>.<NTRIAGE outputs>
    results/<dataset>/k_selection/assignments.tsv.gz

Outputs are written below ``results/<dataset>/cluster_analysis``.  The script
does not hard-code a cluster count and is safe to rerun.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import os
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.stats import kruskal, mannwhitneyu
from sklearn.metrics import adjusted_rand_score


SCRIPT = Path(__file__).resolve()
DEFAULT_ROOT = SCRIPT.parents[1]
DEFAULT_MANIFEST = DEFAULT_ROOT / "inputs" / "dataset_manifest.tsv"
DEFAULT_PRIMATE_MAP = Path(
    os.environ.get(
        "ACRO1_PRIMATE_ARRAY_MAP",
        DEFAULT_ROOT / "inputs" / "six_primate_array_mapping.tsv",
    )
)
REQUIRED_DATASETS = {
    "hg002",
    "gorilla",
    "chimpanzee",
    "bonobo",
    "sumatran_orangutan",
    "bornean_orangutan",
    "joint_primates",
    "hprc",
    "joint_all",
}
PRIMARY_METRICS = [
    "length",
    "acro1_copy_number",
    "classified_acro1_copy_number",
    "raw_acro1_anchor_union_bp",
    "classified_acro1_union_bp",
    "broad_acro1_annotation_bp",
    "explicit_backbone_bp",
    "unclassified_bp",
    "non_anchor_remainder_bp",
]
ALLELE_COLUMNS = ["subfamily_1", "subfamily_2", "subfamily_3"]
SCOPE_NOTE = (
    "descriptive; tandem units are grouped by parent array, and ACRO/ACRO1 "
    "annotations overlap features used to define the clusters"
)


def sha256_file(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument(
        "--datasets",
        nargs="*",
        help="Dataset IDs to analyze (default: all rows in the manifest)",
    )
    parser.add_argument("--joint-id", default="joint_all")
    parser.add_argument(
        "--combined-output-dir",
        type=Path,
        default=Path("results/combined"),
        help=(
            "Directory for aggregate summaries/crosswalks, relative to --root "
            "unless absolute (default: results/combined)"
        ),
    )
    parser.add_argument("--primate-map", type=Path, default=DEFAULT_PRIMATE_MAP)
    parser.add_argument("--bootstrap", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--skip-plots", action="store_true")
    return parser.parse_args()


def resolve_path(value: object, base: Path) -> Path | None:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return None
    text = str(value).strip()
    if not text:
        return None
    path = Path(text).expanduser()
    return path if path.is_absolute() else (base / path).resolve()


def dataset_column(manifest: pd.DataFrame) -> str:
    for candidate in ("dataset_id", "dataset", "cohort", "name"):
        if candidate in manifest:
            return candidate
    raise ValueError("Manifest needs a dataset_id (or dataset/cohort/name) column")


def load_manifest(path: Path, selected: Sequence[str] | None) -> tuple[pd.DataFrame, str]:
    manifest = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    id_col = dataset_column(manifest)
    manifest[id_col] = manifest[id_col].str.strip()
    if manifest[id_col].duplicated().any():
        dup = manifest.loc[manifest[id_col].duplicated(), id_col].tolist()
        raise ValueError(f"Duplicate dataset IDs in manifest: {dup}")
    if selected:
        unknown = sorted(set(selected) - set(manifest[id_col]))
        if unknown:
            raise ValueError(f"Unknown --datasets: {', '.join(unknown)}")
        manifest = manifest.loc[manifest[id_col].isin(selected)].copy()
        order = {name: i for i, name in enumerate(selected)}
        manifest = manifest.sort_values(id_col, key=lambda s: s.map(order))
    return manifest, id_col


def output_prefix(root: Path, dataset: str, row: Mapping[str, object]) -> Path:
    for column in ("output_prefix", "ntriage_prefix", "run_prefix", "prefix"):
        if column in row:
            candidate = resolve_path(row[column], root)
            if candidate:
                text = str(candidate)
                for suffix in (".clusters.tsv", ".features.tsv", ".arrays.tsv"):
                    if text.endswith(suffix):
                        return Path(text[: -len(suffix)])
                return candidate
    return root / "runs" / dataset / "native" / dataset


def selected_assignments_path(root: Path, dataset: str, row: Mapping[str, object]) -> Path:
    for column in ("selected_assignments", "assignments", "cluster_assignments"):
        if column in row:
            candidate = resolve_path(row[column], root)
            if candidate:
                return candidate
    candidates = [
        root / "results" / dataset / "k_selection" / "assignments.tsv.gz",
        root / "results" / dataset / "k_selection" / "assignments.tsv",
        root / "results" / dataset / "selected_clusters.tsv.gz",
        root / "results" / dataset / "selected_clusters.tsv",
    ]
    for path in candidates:
        if path.exists():
            return path
    return candidates[0]


def normalize_assignments(frame: pd.DataFrame, dataset: str) -> pd.DataFrame:
    if "acre_id" not in frame:
        raise ValueError("Selected assignments need acre_id")
    cluster_source = next(
        (c for c in ("selected_cluster", "cluster", "best_cluster") if c in frame),
        None,
    )
    if cluster_source is None:
        raise ValueError("Selected assignments need selected_cluster (or cluster)")
    keep = ["acre_id", cluster_source]
    for optional in (
        "array_id",
        "selected_k",
        "selected_cluster_label",
        "agreement",
        "vote_margin",
        "distance_to_centroid",
        "second_distance",
        "distance_margin",
        "n_final_starts",
    ):
        if optional in frame and optional not in keep:
            keep.append(optional)
    result = frame[keep].copy()
    result = result.rename(columns={cluster_source: "selected_cluster"})
    result["selected_cluster"] = pd.to_numeric(
        result["selected_cluster"], errors="raise"
    ).astype(int)
    if result.acre_id.duplicated().any():
        raise ValueError(f"{dataset}: duplicate acre_id in selected assignments")
    selected_k = (
        int(pd.to_numeric(result.selected_k).iloc[0])
        if "selected_k" in result and len(result)
        else int(result.selected_cluster.nunique())
    )
    if "selected_cluster_label" not in result:
        tag = re.sub(r"[^A-Za-z0-9]+", "_", dataset).strip("_").upper()
        result["selected_cluster_label"] = result.selected_cluster.map(
            lambda c: f"{tag}_K{selected_k:02d}_C{c:02d}"
        )
    result["selected_k"] = selected_k
    return result


def load_units(prefix: Path, assignments: Path, dataset: str) -> pd.DataFrame:
    clusters_path = Path(f"{prefix}.clusters.tsv")
    units = pd.read_csv(clusters_path, sep="\t")
    required = {"acre_id", "array_id", "chrom", "start", "end"}
    missing = sorted(required - set(units))
    if missing:
        raise ValueError(f"{clusters_path}: missing columns {missing}")
    if units.acre_id.duplicated().any():
        raise ValueError(f"{dataset}: duplicate acre_id in NTRIAGE clusters")
    chosen = normalize_assignments(pd.read_csv(assignments, sep="\t"), dataset)
    if set(units.acre_id) != set(chosen.acre_id):
        only_units = len(set(units.acre_id) - set(chosen.acre_id))
        only_chosen = len(set(chosen.acre_id) - set(units.acre_id))
        raise ValueError(
            f"{dataset}: assignment/unit ID mismatch "
            f"({only_units} only in units; {only_chosen} only in assignments)"
        )
    if "cluster" in units:
        units = units.rename(columns={"cluster": "ntriage_initial_cluster"})
    # Do not duplicate the NTRIAGE array_id column from assignments.
    chosen = chosen.drop(columns=["array_id"], errors="ignore")
    units = units.merge(chosen, on="acre_id", how="left", validate="one_to_one")
    units["cluster"] = units.selected_cluster
    units["cluster_label"] = units.selected_cluster_label
    units["length"] = pd.to_numeric(
        units.get("length", units.end - units.start), errors="raise"
    ).astype(int)
    units["unit_key"] = (
        units.chrom.astype(str)
        + ":"
        + units.start.astype(str)
        + "-"
        + units.end.astype(str)
    )
    return units


def merged_lengths(frame: pd.DataFrame, key_columns: Sequence[str]) -> pd.DataFrame:
    """Return interval-union bp for every key without assuming non-overlap."""
    output_columns = list(key_columns) + ["union_bp"]
    if frame.empty:
        return pd.DataFrame(columns=output_columns)
    ordered = frame.sort_values(list(key_columns) + ["start", "end"])
    totals: list[tuple[object, ...]] = []
    old_key: tuple[object, ...] | None = None
    current_start = current_end = total = 0
    for row in ordered[list(key_columns) + ["start", "end"]].itertuples(
        index=False, name=None
    ):
        key = tuple(row[: len(key_columns)])
        start, end = int(row[-2]), int(row[-1])
        if old_key != key:
            if old_key is not None:
                total += current_end - current_start
                totals.append((*old_key, total))
            old_key = key
            current_start, current_end, total = start, end, 0
        elif start <= current_end:
            current_end = max(current_end, end)
        else:
            total += current_end - current_start
            current_start, current_end = start, end
    if old_key is not None:
        total += current_end - current_start
        totals.append((*old_key, total))
    return pd.DataFrame(totals, columns=output_columns)


def clip_intervals_to_units(
    intervals: pd.DataFrame, units: pd.DataFrame
) -> tuple[pd.DataFrame, dict[str, int]]:
    bounds = units[["acre_id", "start", "end"]].rename(
        columns={"start": "unit_start", "end": "unit_end"}
    )
    before = len(intervals)
    intervals = intervals.merge(bounds, on="acre_id", how="inner", validate="many_to_one")
    missing_unit = before - len(intervals)
    intervals["start"] = np.maximum(intervals.start, intervals.unit_start)
    intervals["end"] = np.minimum(intervals.end, intervals.unit_end)
    nonpositive = int((intervals.end <= intervals.start).sum())
    intervals = intervals.loc[intervals.end > intervals.start].drop(
        columns=["unit_start", "unit_end"]
    )
    return intervals, {
        "input_intervals": before,
        "intervals_without_matching_unit": missing_unit,
        "nonpositive_after_clipping": nonpositive,
        "retained_intervals": len(intervals),
    }


def load_annotation_metrics(prefix: Path, units: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, int]]:
    path = Path(f"{prefix}.annotations.bed")
    try:
        ann = pd.read_csv(
            path,
            sep="\t",
            header=None,
            usecols=[1, 2, 3, 9],
            low_memory=False,
        )
        ann.columns = ["start", "end", "family", "acre_id"]
        # Do not use pandas' ``comment='#'`` here: HPRC sequence names contain
        # embedded '#'.  Numeric coercion safely removes an optional BED header.
        ann["start"] = pd.to_numeric(ann.start, errors="coerce")
        ann["end"] = pd.to_numeric(ann.end, errors="coerce")
        ann = ann.dropna(subset=["start", "end"]).copy()
        ann[["start", "end"]] = ann[["start", "end"]].astype(int)
    except pd.errors.EmptyDataError:
        ann = pd.DataFrame(columns=["start", "end", "family", "acre_id"])
    ann, qc = clip_intervals_to_units(ann, units)
    if ann.empty:
        base = units[["acre_id"]].copy()
        for column in (
            "annotation_covered_bp",
            "broad_acro1_annotation_bp",
            "explicit_backbone_bp",
            "unclassified_bp",
        ):
            base[column] = 0
        return base, qc

    all_union = merged_lengths(ann, ["acre_id"]).rename(
        columns={"union_bp": "annotation_covered_bp"}
    )
    ann["category"] = np.select(
        [ann.family.eq("ACRO1"), ann.family.eq("unclassified")],
        ["broad_acro1_annotation_bp", "unclassified_bp"],
        default="explicit_backbone_bp",
    )
    category = merged_lengths(ann, ["acre_id", "category"])
    category = category.pivot(index="acre_id", columns="category", values="union_bp")
    category.columns.name = None
    category = category.reset_index()

    family = merged_lengths(ann, ["acre_id", "family"])
    family = family.pivot(index="acre_id", columns="family", values="union_bp")
    family.columns = [f"family_bp:{name}" for name in family.columns]
    family = family.reset_index()

    metrics = units[["acre_id"]].merge(all_union, on="acre_id", how="left")
    metrics = metrics.merge(category, on="acre_id", how="left")
    metrics = metrics.merge(family, on="acre_id", how="left")
    numeric = [c for c in metrics if c != "acre_id"]
    metrics[numeric] = metrics[numeric].fillna(0).astype(int)
    for column in (
        "annotation_covered_bp",
        "broad_acro1_annotation_bp",
        "explicit_backbone_bp",
        "unclassified_bp",
    ):
        if column not in metrics:
            metrics[column] = 0
    qc["annotation_families"] = int(ann.family.nunique())
    return metrics, qc


def load_anchor_metrics(prefix: Path, units: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, int]]:
    path = Path(f"{prefix}.acro1_subfamilies.bed")
    try:
        anchor = pd.read_csv(
            path,
            sep="\t",
            header=None,
            usecols=[1, 2, 9],
            low_memory=False,
        )
        anchor.columns = ["start", "end", "acre_id"]
        anchor["start"] = pd.to_numeric(anchor.start, errors="coerce")
        anchor["end"] = pd.to_numeric(anchor.end, errors="coerce")
        anchor = anchor.dropna(subset=["start", "end"]).copy()
        anchor[["start", "end"]] = anchor[["start", "end"]].astype(int)
    except pd.errors.EmptyDataError:
        anchor = pd.DataFrame(columns=["start", "end", "acre_id"])
    anchor, qc = clip_intervals_to_units(anchor, units)
    union = merged_lengths(anchor, ["acre_id"]).rename(
        columns={"union_bp": "classified_acro1_union_bp"}
    )
    counts = (
        anchor.groupby("acre_id").size().rename("classified_acro1_anchor_count").reset_index()
    )
    metrics = units[["acre_id"]].merge(union, on="acre_id", how="left")
    metrics = metrics.merge(counts, on="acre_id", how="left").fillna(0)
    metrics[["classified_acro1_union_bp", "classified_acro1_anchor_count"]] = metrics[
        ["classified_acro1_union_bp", "classified_acro1_anchor_count"]
    ].astype(int)
    return metrics, qc


def load_raw_anchor_metrics(
    path: Path | None, units: pd.DataFrame
) -> tuple[pd.DataFrame, dict[str, int]]:
    """Assign every input ACRO1 monomer anchor by midpoint to a complete unit.

    NTRIAGE's ``acro1_subfamilies.bed`` retains only anchors participating in
    exact 144/147-bp chains.  The raw input BED is therefore the less biased
    primary copy-number proxy; classified anchors remain available separately
    for the S1/S2/S3 composition analysis.
    """
    empty = units[["acre_id"]].copy()
    empty["raw_acro1_anchor_count"] = 0
    empty["raw_acro1_union_bp"] = 0
    if path is None or not path.exists():
        return empty, {
            "input_intervals": 0,
            "assigned_by_midpoint": 0,
            "unassigned_to_complete_unit": 0,
            "boundary_clipped": 0,
            "raw_anchor_input_missing": 1,
        }
    try:
        raw = pd.read_csv(
            path,
            sep="\t",
            header=None,
            usecols=[0, 1, 2],
            low_memory=False,
        )
    except pd.errors.EmptyDataError:
        return empty, {
            "input_intervals": 0,
            "assigned_by_midpoint": 0,
            "unassigned_to_complete_unit": 0,
            "boundary_clipped": 0,
            "raw_anchor_input_missing": 0,
        }
    raw.columns = ["chrom", "start", "end"]
    raw["start"] = pd.to_numeric(raw.start, errors="coerce")
    raw["end"] = pd.to_numeric(raw.end, errors="coerce")
    raw = raw.dropna(subset=["start", "end"]).copy()
    raw[["start", "end"]] = raw[["start", "end"]].astype(int)
    raw = raw.loc[raw.end > raw.start]
    input_count = len(raw)
    assigned_parts: list[pd.DataFrame] = []
    unit_groups = {str(chrom): group.sort_values("start") for chrom, group in units.groupby("chrom")}
    for chrom, anchors in raw.groupby(raw.chrom.astype(str), sort=False):
        target = unit_groups.get(str(chrom))
        if target is None or target.empty:
            continue
        starts = target.start.to_numpy(dtype=np.int64)
        ends = target.end.to_numpy(dtype=np.int64)
        mids = ((anchors.start.to_numpy(dtype=np.int64) + anchors.end.to_numpy(dtype=np.int64)) // 2)
        index = np.searchsorted(starts, mids, side="right") - 1
        valid = (index >= 0)
        valid_positions = np.flatnonzero(valid)
        valid[valid_positions] &= mids[valid_positions] < ends[index[valid_positions]]
        if not valid.any():
            continue
        anchor_part = anchors.iloc[np.flatnonzero(valid)].copy()
        target_positions = index[valid]
        anchor_part["acre_id"] = target.iloc[target_positions].acre_id.to_numpy()
        anchor_part["unit_start"] = starts[target_positions]
        anchor_part["unit_end"] = ends[target_positions]
        anchor_part["was_boundary_clipped"] = (
            anchor_part.start.lt(anchor_part.unit_start)
            | anchor_part.end.gt(anchor_part.unit_end)
        )
        anchor_part["start"] = np.maximum(anchor_part.start, anchor_part.unit_start)
        anchor_part["end"] = np.minimum(anchor_part.end, anchor_part.unit_end)
        anchor_part = anchor_part.loc[anchor_part.end > anchor_part.start]
        assigned_parts.append(anchor_part)
    if assigned_parts:
        assigned = pd.concat(assigned_parts, ignore_index=True)
        union = merged_lengths(assigned, ["acre_id"]).rename(
            columns={"union_bp": "raw_acro1_union_bp"}
        )
        counts = (
            assigned.groupby("acre_id").size().rename("raw_acro1_anchor_count").reset_index()
        )
        result = empty[["acre_id"]].merge(union, on="acre_id", how="left")
        result = result.merge(counts, on="acre_id", how="left").fillna(0)
        result[["raw_acro1_union_bp", "raw_acro1_anchor_count"]] = result[
            ["raw_acro1_union_bp", "raw_acro1_anchor_count"]
        ].astype(int)
        clipped = int(assigned.was_boundary_clipped.sum())
    else:
        assigned = pd.DataFrame()
        result = empty
        clipped = 0
    assigned_count = len(assigned)
    return result, {
        "input_intervals": input_count,
        "assigned_by_midpoint": assigned_count,
        "unassigned_to_complete_unit": input_count - assigned_count,
        "boundary_clipped": clipped,
        "raw_anchor_input_missing": 0,
    }


def load_composition(prefix: Path, units: pd.DataFrame) -> pd.DataFrame:
    path = Path(f"{prefix}.acro1_composition.tsv")
    composition = pd.read_csv(path, sep="\t")
    if "acre_id" not in composition:
        raise ValueError(f"{path}: missing acre_id")
    aliases = {
        "S1": "subfamily_1",
        "S2": "subfamily_2",
        "S3": "subfamily_3",
        "subfamily1": "subfamily_1",
        "subfamily2": "subfamily_2",
        "subfamily3": "subfamily_3",
    }
    composition = composition.rename(
        columns={k: v for k, v in aliases.items() if k in composition}
    )
    for column in ALLELE_COLUMNS:
        if column not in composition:
            composition[column] = 0
    keep = ["acre_id", *ALLELE_COLUMNS]
    for column in ("dominant_subfamily", "composition_class", "architecture_rle"):
        if column in composition:
            keep.append(column)
    composition = composition[keep]
    if composition.acre_id.duplicated().any():
        raise ValueError(f"{path}: duplicate acre_id")
    composition[ALLELE_COLUMNS] = composition[ALLELE_COLUMNS].fillna(0).astype(int)
    result = units[["acre_id"]].merge(
        composition, on="acre_id", how="left", validate="one_to_one"
    )
    result[ALLELE_COLUMNS] = result[ALLELE_COLUMNS].fillna(0).astype(int)
    for column in set(keep) - {"acre_id", *ALLELE_COLUMNS}:
        result[column] = result[column].fillna("none")
    return result


def load_species_map(
    manifest: pd.DataFrame, id_col: str, manifest_dir: Path, primate_map: Path
) -> dict[str, str]:
    mapping: dict[str, str] = {}
    candidate_paths: list[Path] = []
    if primate_map.exists():
        candidate_paths.append(primate_map)
    for column in ("array_metadata", "mapping", "array_mapping"):
        if column in manifest:
            for value in manifest[column].unique():
                path = resolve_path(value, manifest_dir)
                if path and path.exists() and path not in candidate_paths:
                    candidate_paths.append(path)
    for path in candidate_paths:
        try:
            table = pd.read_csv(path, sep="\t", dtype=str)
        except Exception:
            continue
        chrom_col = next(
            (c for c in ("synthetic_chrom", "chrom", "sequence_id") if c in table),
            None,
        )
        species_col = next(
            # ``source_dataset`` is the canonical cohort ID in the prepared
            # metadata (e.g. ``gorilla``); raw ``species_id`` may instead be
            # an assembly tag such as mGorGor1.
            (c for c in ("source_dataset", "dataset_id", "cohort", "species_id") if c in table),
            None,
        )
        if chrom_col and species_col:
            mapping.update(dict(zip(table[chrom_col], table[species_col])))
    return mapping


def add_provenance_columns(
    units: pd.DataFrame,
    dataset: str,
    species_map: Mapping[str, str],
    joint_id: str = "joint_all",
) -> pd.DataFrame:
    units = units.copy()
    # Both shipped aggregate datasets need their source species recovered from
    # synthetic sequence IDs. ``joint_id`` keeps the behavior extensible for a
    # caller-defined aggregate name.
    if dataset == joint_id or dataset in {"joint_all", "joint_primates"}:
        units["species_id"] = units.chrom.astype(str).map(species_map)
        units.loc[units.species_id.isna() & ~units.chrom.astype(str).str.startswith("PF"), "species_id"] = "hprc"
        units["species_id"] = units.species_id.fillna("unknown_primate")
    elif dataset == "hprc":
        units["species_id"] = "hprc"
    else:
        units["species_id"] = dataset
    source = units.get("source_name", units.chrom).astype(str)
    units["individual_id"] = np.where(
        units.species_id.eq("hprc"), source.str.split("#").str[0], units.species_id
    )
    return units


def build_unit_metrics(
    prefix: Path,
    units: pd.DataFrame,
    dataset: str,
    species_map: Mapping[str, str],
    raw_anchor_path: Path | None = None,
    joint_id: str = "joint_all",
) -> tuple[pd.DataFrame, dict[str, object]]:
    ann, ann_qc = load_annotation_metrics(prefix, units)
    classified_anchor, classified_anchor_qc = load_anchor_metrics(prefix, units)
    raw_anchor, raw_anchor_qc = load_raw_anchor_metrics(raw_anchor_path, units)
    composition = load_composition(prefix, units)
    metrics = units.merge(ann, on="acre_id", how="left", validate="one_to_one")
    metrics = metrics.merge(
        classified_anchor, on="acre_id", how="left", validate="one_to_one"
    )
    metrics = metrics.merge(raw_anchor, on="acre_id", how="left", validate="one_to_one")
    metrics = metrics.merge(composition, on="acre_id", how="left", validate="one_to_one")
    metrics = add_provenance_columns(metrics, dataset, species_map, joint_id=joint_id)
    metrics["classified_acro1_copy_number"] = metrics[ALLELE_COLUMNS].sum(axis=1)
    metrics["classified_total"] = metrics.classified_acro1_copy_number
    # Primary copy-number proxy: all raw mapped monomer anchors, including
    # anchors that do not participate in exact 144/147-bp allele chains.
    metrics["acro1_copy_number"] = metrics.raw_acro1_anchor_count
    metrics["acro1_anchor_count"] = metrics.raw_acro1_anchor_count
    # Preferred explicit name plus a backward-compatible legacy alias.  These
    # bases are the union of mapped raw-anchor intervals; they are not an exact
    # delineation of every ACRO1-derived base.
    metrics["raw_acro1_anchor_union_bp"] = metrics.raw_acro1_union_bp
    metrics["exact_acro1_union_bp"] = metrics.raw_acro1_anchor_union_bp
    metrics["classified_copy_anchor_count_difference"] = (
        metrics.classified_acro1_copy_number - metrics.classified_acro1_anchor_count
    )
    denominator = metrics.classified_total.replace(0, np.nan)
    for i, column in enumerate(ALLELE_COLUMNS, start=1):
        metrics[f"fraction_subfamily_{i}"] = (metrics[column] / denominator).fillna(0.0)
    metrics["acro1_anchors_per_kb"] = (
        1000.0 * metrics.raw_acro1_anchor_count / metrics.length
    )
    metrics["raw_acro1_anchors_per_kb"] = metrics.acro1_anchors_per_kb
    metrics["classified_acro1_anchors_per_kb"] = (
        1000.0 * metrics.classified_acro1_anchor_count / metrics.length
    )
    metrics["non_anchor_remainder_bp"] = (
        metrics.length - metrics.raw_acro1_anchor_union_bp
    )
    # Backward-compatible legacy alias.  The remainder is not proven to be
    # non-ACRO1, so new outputs and prose use ``non_anchor_remainder_bp``.
    metrics["non_acro1_remainder_bp"] = metrics.non_anchor_remainder_bp
    metrics["annotation_uncovered_bp"] = metrics.length - metrics.annotation_covered_bp
    metrics["annotation_category_sum_bp"] = (
        metrics.broad_acro1_annotation_bp
        + metrics.explicit_backbone_bp
        + metrics.unclassified_bp
    )
    metrics["annotation_accounting_difference_bp"] = (
        metrics.annotation_category_sum_bp - metrics.annotation_covered_bp
    )
    metrics["length_bin"] = pd.cut(
        metrics.length,
        bins=[5999, 6999, 7999, 9000],
        labels=["6_to_lt7kb", "7_to_lt8kb", "8_to_9kb"],
        include_lowest=True,
        right=True,
    ).astype("string")
    qc: dict[str, object] = {
        "dataset_id": dataset,
        "units": len(metrics),
        "arrays": int(metrics.array_id.nunique()),
        "clusters": int(metrics.cluster.nunique()),
        "unit_length_min": int(metrics.length.min()) if len(metrics) else None,
        "unit_length_max": int(metrics.length.max()) if len(metrics) else None,
        "annotation": ann_qc,
        "classified_acro1_anchors": classified_anchor_qc,
        "raw_acro1_anchors": raw_anchor_qc,
        "units_with_zero_classified_acro1": int(metrics.classified_total.eq(0).sum()),
        "units_with_zero_raw_acro1_anchors": int(metrics.raw_acro1_anchor_count.eq(0).sum()),
        "units_with_classified_copy_anchor_count_difference": int(
            metrics.classified_copy_anchor_count_difference.ne(0).sum()
        ),
        "units_with_unknown_primate_provenance": int(
            metrics.species_id.eq("unknown_primate").sum()
        ),
        "max_abs_annotation_accounting_difference_bp": int(
            metrics.annotation_accounting_difference_bp.abs().max()
        ) if len(metrics) else None,
        "interpretation_scope": SCOPE_NOTE,
    }
    return metrics, qc


def quantile(series: pd.Series, q: float) -> float:
    clean = pd.to_numeric(series, errors="coerce").dropna()
    return float(clean.quantile(q)) if len(clean) else math.nan


def array_cluster_table(units: pd.DataFrame) -> pd.DataFrame:
    keys = ["cluster", "cluster_label", "array_id", "species_id", "individual_id"]
    aggregation: dict[str, tuple[str, str]] = {"units": ("acre_id", "size")}
    for metric in PRIMARY_METRICS + ["annotation_covered_bp", "annotation_uncovered_bp"]:
        aggregation[f"mean_{metric}"] = (metric, "mean")
        aggregation[f"median_{metric}"] = (metric, "median")
    for allele in ALLELE_COLUMNS:
        aggregation[f"sum_{allele}"] = (allele, "sum")
    result = units.groupby(keys, observed=True).agg(**aggregation).reset_index()
    for statistic in ("mean", "median"):
        result[f"{statistic}_exact_acro1_union_bp"] = result[
            f"{statistic}_raw_acro1_anchor_union_bp"
        ]
        result[f"{statistic}_non_acro1_remainder_bp"] = result[
            f"{statistic}_non_anchor_remainder_bp"
        ]
    total = result[[f"sum_{c}" for c in ALLELE_COLUMNS]].sum(axis=1)
    for i, allele in enumerate(ALLELE_COLUMNS, start=1):
        result[f"fraction_subfamily_{i}"] = (
            result[f"sum_{allele}"] / total.replace(0, np.nan)
        )
    return result


def cluster_summary(units: pd.DataFrame, arrays: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for cluster, group in units.groupby("cluster", sort=True):
        label = str(group.cluster_label.iloc[0])
        ag = arrays.loc[arrays.cluster.eq(cluster)]
        counts = group[ALLELE_COLUMNS].sum()
        classified = int(counts.sum())
        raw_anchor_total = int(group.raw_acro1_anchor_count.sum())
        array_classification_fraction = (
            ag.mean_classified_acro1_copy_number
            / ag.mean_acro1_copy_number.replace(0, np.nan)
        )
        row: dict[str, object] = {
            "cluster": int(cluster),
            "cluster_label": label,
            "selected_k": int(group.selected_k.iloc[0]),
            "units": len(group),
            "fraction_of_units": len(group) / len(units),
            "arrays": int(group.array_id.nunique()),
            "individuals": int(group.individual_id.nunique()),
            "species": int(group.species_id.nunique()),
            "mean_unit_length_bp": float(group.length.mean()),
            "median_unit_length_bp": float(group.length.median()),
            "q25_unit_length_bp": quantile(group.length, 0.25),
            "q75_unit_length_bp": quantile(group.length, 0.75),
            "min_unit_length_bp": int(group.length.min()),
            "max_unit_length_bp": int(group.length.max()),
            "array_weighted_mean_unit_length_bp": float(ag.mean_length.mean()),
            "mean_acro1_copy_number": float(group.acro1_copy_number.mean()),
            "median_acro1_copy_number": float(group.acro1_copy_number.median()),
            "array_weighted_mean_acro1_copy_number": float(
                ag.mean_acro1_copy_number.mean()
            ),
            "mean_raw_acro1_anchor_count": float(group.raw_acro1_anchor_count.mean()),
            "median_raw_acro1_anchor_count": float(group.raw_acro1_anchor_count.median()),
            "mean_raw_acro1_anchors_per_kb": float(group.raw_acro1_anchors_per_kb.mean()),
            "units_without_raw_acro1_anchors": int(group.raw_acro1_anchor_count.eq(0).sum()),
            "fraction_units_without_raw_acro1_anchors": float(
                group.raw_acro1_anchor_count.eq(0).mean()
            ),
            "mean_classified_acro1_copy_number": float(
                group.classified_acro1_copy_number.mean()
            ),
            "median_classified_acro1_copy_number": float(
                group.classified_acro1_copy_number.median()
            ),
            "array_weighted_mean_classified_acro1_copy_number": float(
                ag.mean_classified_acro1_copy_number.mean()
            ),
            "units_without_classified_acro1": int(
                group.classified_acro1_copy_number.eq(0).sum()
            ),
            "fraction_units_without_classified_acro1": float(
                group.classified_acro1_copy_number.eq(0).mean()
            ),
            "subfamily_1_count": int(counts.subfamily_1),
            "subfamily_2_count": int(counts.subfamily_2),
            "subfamily_3_count": int(counts.subfamily_3),
            "classified_acro1_total": classified,
            "raw_acro1_anchor_total": raw_anchor_total,
            "unclassified_raw_acro1_anchor_total": max(
                raw_anchor_total - classified, 0
            ),
            "fraction_raw_acro1_anchors_classified": (
                classified / raw_anchor_total if raw_anchor_total else 0
            ),
            "fraction_raw_acro1_anchors_classified_array_weighted": float(
                array_classification_fraction.mean()
            ),
            "arrays_with_classified_acro1": int(ag.fraction_subfamily_1.notna().sum()),
            "fraction_subfamily_1_pooled": counts.subfamily_1 / classified if classified else 0,
            "fraction_subfamily_2_pooled": counts.subfamily_2 / classified if classified else 0,
            "fraction_subfamily_3_pooled": counts.subfamily_3 / classified if classified else 0,
            "fraction_subfamily_1_array_weighted": float(ag.fraction_subfamily_1.mean()),
            "fraction_subfamily_2_array_weighted": float(ag.fraction_subfamily_2.mean()),
            "fraction_subfamily_3_array_weighted": float(ag.fraction_subfamily_3.mean()),
            "interpretation_scope": SCOPE_NOTE,
        }
        for metric in (
            "raw_acro1_anchor_union_bp",
            "classified_acro1_union_bp",
            "broad_acro1_annotation_bp",
            "explicit_backbone_bp",
            "unclassified_bp",
            "non_anchor_remainder_bp",
            "annotation_covered_bp",
            "annotation_uncovered_bp",
        ):
            row[f"mean_{metric}"] = float(group[metric].mean())
            row[f"median_{metric}"] = float(group[metric].median())
            array_col = f"mean_{metric}"
            row[f"array_weighted_mean_{metric}"] = float(ag[array_col].mean())
        # Backward-compatible summary aliases.  New user-facing tables use the
        # explicit raw-anchor/non-anchor names above.
        for statistic in ("mean", "median", "array_weighted_mean"):
            row[f"{statistic}_exact_acro1_union_bp"] = row[
                f"{statistic}_raw_acro1_anchor_union_bp"
            ]
            row[f"{statistic}_non_acro1_remainder_bp"] = row[
                f"{statistic}_non_anchor_remainder_bp"
            ]
        rows.append(row)
    return pd.DataFrame(rows).sort_values("cluster")


def length_bin_table(units: pd.DataFrame) -> pd.DataFrame:
    bins = ["6_to_lt7kb", "7_to_lt8kb", "8_to_9kb", pd.NA]
    result = (
        units.groupby(["cluster", "cluster_label", "length_bin"], dropna=False, observed=False)
        .size()
        .rename("units")
        .reset_index()
    )
    totals = result.groupby("cluster").units.transform("sum")
    result["fraction_within_cluster"] = result.units / totals
    result["bin_definition"] = result.length_bin.map(
        {
            "6_to_lt7kb": "6000-6999 bp",
            "7_to_lt8kb": "7000-7999 bp",
            "8_to_9kb": "8000-9000 bp",
        }
    ).fillna("outside requested 6000-9000 bp range")
    return result


def component_tables(units: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    family_columns = [c for c in units if c.startswith("family_bp:")]
    family_rows: list[dict[str, object]] = []
    for cluster, group in units.groupby("cluster", sort=True):
        label = group.cluster_label.iloc[0]
        total_length = float(group.length.sum())
        for column in family_columns:
            bp = int(group[column].sum())
            family_rows.append(
                {
                    "cluster": int(cluster),
                    "cluster_label": label,
                    "component": column.removeprefix("family_bp:"),
                    "bp": bp,
                    "mean_bp_per_unit": float(group[column].mean()),
                    "fraction_of_cluster_unit_bp": bp / total_length if total_length else 0,
                }
            )
    family = pd.DataFrame(family_rows)
    if len(family):
        family["rank_by_bp_within_cluster"] = (
            family.groupby("cluster").bp.rank(method="first", ascending=False).astype(int)
        )
        family = family.sort_values(["cluster", "rank_by_bp_within_cluster"])

    if "top_component" in units:
        top = (
            units.groupby(["cluster", "cluster_label", "top_component"], observed=True)
            .size()
            .rename("units")
            .reset_index()
        )
        top["fraction_within_cluster"] = top.units / top.groupby("cluster").units.transform("sum")
        top["rank_within_cluster"] = (
            top.groupby("cluster").units.rank(method="first", ascending=False).astype(int)
        )
        top = top.sort_values(["cluster", "rank_within_cluster"])
    else:
        top = pd.DataFrame(
            columns=["cluster", "cluster_label", "top_component", "units", "fraction_within_cluster"]
        )
    return family, top


def bootstrap_array_intervals(
    arrays: pd.DataFrame, repetitions: int, seed: int
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    metrics = [f"mean_{metric}" for metric in PRIMARY_METRICS] + [
        "fraction_subfamily_1",
        "fraction_subfamily_2",
        "fraction_subfamily_3",
    ]
    rows: list[dict[str, object]] = []
    for cluster, group in arrays.groupby("cluster", sort=True):
        label = group.cluster_label.iloc[0]
        for metric in metrics:
            values = pd.to_numeric(group[metric], errors="coerce").dropna().to_numpy(float)
            if not len(values):
                continue
            estimate = float(values.mean())
            if len(values) == 1 or repetitions == 0:
                low = high = estimate
            else:
                # Chunking avoids a large B x n matrix for HPRC-wide clusters.
                boot = np.empty(repetitions, dtype=float)
                for begin in range(0, repetitions, 100):
                    count = min(100, repetitions - begin)
                    indices = rng.integers(0, len(values), size=(count, len(values)))
                    boot[begin : begin + count] = values[indices].mean(axis=1)
                low, high = np.quantile(boot, [0.025, 0.975])
            rows.append(
                {
                    "cluster": int(cluster),
                    "cluster_label": label,
                    "metric": metric,
                    "array_weighted_mean": estimate,
                    "bootstrap_95ci_low": float(low),
                    "bootstrap_95ci_high": float(high),
                    "arrays": len(values),
                    "bootstrap_repetitions": repetitions,
                    "resampling_unit": "parent array within cluster",
                }
            )
    return pd.DataFrame(rows)


def bh_adjust(pvalues: Iterable[float]) -> np.ndarray:
    p = np.asarray(list(pvalues), dtype=float)
    result = np.full(len(p), np.nan)
    finite = np.isfinite(p)
    if not finite.any():
        return result
    indices = np.where(finite)[0]
    order = indices[np.argsort(p[finite])]
    ranked = p[order] * len(order) / np.arange(1, len(order) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    result[order] = np.minimum(ranked, 1.0)
    return result


def descriptive_tests(arrays: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    metrics = [f"mean_{metric}" for metric in PRIMARY_METRICS]
    clusters = sorted(arrays.cluster.unique())
    for metric in metrics:
        groups = {
            int(c): pd.to_numeric(
                arrays.loc[arrays.cluster.eq(c), metric], errors="coerce"
            ).dropna().to_numpy(float)
            for c in clusters
        }
        usable = [v for v in groups.values() if len(v)]
        try:
            statistic, pvalue = kruskal(*usable) if len(usable) >= 2 else (np.nan, np.nan)
        except ValueError:
            statistic, pvalue = np.nan, np.nan
        rows.append(
            {
                "metric": metric,
                "test": "Kruskal-Wallis",
                "cluster_a": "all",
                "cluster_b": "all",
                "n_a": sum(map(len, usable)),
                "n_b": len(usable),
                "statistic": statistic,
                "median_difference_a_minus_b": np.nan,
                "p_value": pvalue,
                "multiple_testing_family": f"{metric}:omnibus",
                "interpretation_scope": SCOPE_NOTE,
            }
        )
        pair_indices: list[int] = []
        for a, b in itertools.combinations(clusters, 2):
            av, bv = groups[int(a)], groups[int(b)]
            try:
                statistic, pvalue = (
                    mannwhitneyu(av, bv, alternative="two-sided")
                    if len(av) and len(bv)
                    else (np.nan, np.nan)
                )
            except ValueError:
                statistic, pvalue = np.nan, np.nan
            pair_indices.append(len(rows))
            rows.append(
                {
                    "metric": metric,
                    "test": "Mann-Whitney U",
                    "cluster_a": int(a),
                    "cluster_b": int(b),
                    "n_a": len(av),
                    "n_b": len(bv),
                    "statistic": statistic,
                    "median_difference_a_minus_b": (
                        float(np.median(av) - np.median(bv)) if len(av) and len(bv) else np.nan
                    ),
                    "p_value": pvalue,
                    "multiple_testing_family": f"{metric}:all_pairwise_clusters",
                    "interpretation_scope": SCOPE_NOTE,
                }
            )
        adjusted = bh_adjust([rows[i]["p_value"] for i in pair_indices])
        for index, qvalue in zip(pair_indices, adjusted):
            rows[index]["bh_q_value"] = qvalue
    result = pd.DataFrame(rows)
    if "bh_q_value" not in result:
        result["bh_q_value"] = np.nan
    result.loc[result.test.eq("Kruskal-Wallis"), "bh_q_value"] = result.loc[
        result.test.eq("Kruskal-Wallis"), "p_value"
    ]
    return result


def species_cluster_summary(units: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for (species, cluster), group in units.groupby(["species_id", "cluster"], sort=True):
        counts = group[ALLELE_COLUMNS].sum()
        total = int(counts.sum())
        rows.append(
            {
                "species_id": species,
                "cluster": int(cluster),
                "cluster_label": group.cluster_label.iloc[0],
                "units": len(group),
                "arrays": int(group.array_id.nunique()),
                "mean_unit_length_bp": float(group.length.mean()),
                "median_unit_length_bp": float(group.length.median()),
                "mean_acro1_copy_number": float(group.acro1_copy_number.mean()),
                "mean_raw_acro1_anchor_count": float(group.raw_acro1_anchor_count.mean()),
                "mean_classified_acro1_copy_number": float(
                    group.classified_acro1_copy_number.mean()
                ),
                "subfamily_1_count": int(counts.subfamily_1),
                "subfamily_2_count": int(counts.subfamily_2),
                "subfamily_3_count": int(counts.subfamily_3),
                "fraction_subfamily_1": counts.subfamily_1 / total if total else 0,
                "fraction_subfamily_2": counts.subfamily_2 / total if total else 0,
                "fraction_subfamily_3": counts.subfamily_3 / total if total else 0,
                "mean_raw_acro1_anchor_union_bp": float(
                    group.raw_acro1_anchor_union_bp.mean()
                ),
                "mean_classified_acro1_union_bp": float(
                    group.classified_acro1_union_bp.mean()
                ),
                "mean_exact_acro1_union_bp": float(
                    group.raw_acro1_anchor_union_bp.mean()
                ),
                "mean_explicit_backbone_bp": float(group.explicit_backbone_bp.mean()),
            }
        )
    return pd.DataFrame(rows)


def make_plots(
    units: pd.DataFrame, arrays: pd.DataFrame, summary: pd.DataFrame, output: Path
) -> None:
    os.environ.setdefault("MPLCONFIGDIR", str(output.parent.parent / ".mpl_cache"))
    os.environ.setdefault("XDG_CACHE_HOME", str(output.parent.parent / ".cache"))
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    output.mkdir(parents=True, exist_ok=True)
    clusters = summary.cluster.astype(int).tolist()
    labels = summary.cluster_label.astype(str).tolist()
    x = np.arange(len(clusters))
    colors = ["#4E79A7", "#F28E2B", "#59A14F"]

    fig, axes = plt.subplots(2, 2, figsize=(max(10, len(clusters) * 1.15), 9))
    bottom = np.zeros(len(clusters))
    for i, (column, name) in enumerate(
        [
            ("fraction_subfamily_1_pooled", "S1 canonical"),
            ("fraction_subfamily_2_pooled", "S2 internal deletion"),
            ("fraction_subfamily_3_pooled", "S3 terminal deletion"),
        ]
    ):
        values = summary[column].to_numpy(float)
        axes[0, 0].bar(x, values, bottom=bottom, label=name, color=colors[i])
        bottom += values
    axes[0, 0].set_ylabel("Pooled classified ACRO1 fraction")
    axes[0, 0].set_ylim(0, 1)
    axes[0, 0].legend(frameon=False, fontsize=8)

    length_groups = [
        arrays.loc[arrays.cluster.eq(c), "mean_length"].dropna().to_numpy() for c in clusters
    ]
    axes[0, 1].boxplot(length_groups, showfliers=False)
    axes[0, 1].axhline(7000, color="#777777", linestyle="--", linewidth=0.8)
    axes[0, 1].axhline(8000, color="#777777", linestyle="--", linewidth=0.8)
    axes[0, 1].set_ylabel("Mean unit length per parent array (bp)")

    copy_groups = [
        arrays.loc[arrays.cluster.eq(c), "mean_acro1_copy_number"].dropna().to_numpy()
        for c in clusters
    ]
    axes[1, 0].boxplot(copy_groups, showfliers=False)
    axes[1, 0].set_ylabel(
        "Mean raw ACRO1 anchor count per unit\n(copy-number proxy; parent-array means)"
    )

    width = 0.36
    axes[1, 1].bar(
        x - width / 2,
        summary.array_weighted_mean_broad_acro1_annotation_bp,
        width,
        label="Broad ACRO1 annotation",
        color="#B07AA1",
    )
    axes[1, 1].bar(
        x + width / 2,
        summary.array_weighted_mean_explicit_backbone_bp,
        width,
        label="Explicit annotated backbone",
        color="#76B7B2",
    )
    axes[1, 1].set_ylabel("Array-weighted mean bp per unit")
    axes[1, 1].legend(frameon=False, fontsize=8)

    for axis in axes.flat:
        axis.set_xticks(x)
        axis.set_xticklabels(labels, rotation=45, ha="right", fontsize=7)
        axis.spines[["top", "right"]].set_visible(False)
    fig.suptitle("NTRIAGE cluster biology (array-aware summaries)")
    fig.tight_layout()
    fig.savefig(output / "cluster_biology_overview.png", dpi=220, bbox_inches="tight")
    fig.savefig(output / "cluster_biology_overview.pdf", bbox_inches="tight")
    plt.close(fig)

    counts = pd.crosstab(units.cluster, units.length_bin).reindex(index=clusters, fill_value=0)
    counts = counts.reindex(columns=["6_to_lt7kb", "7_to_lt8kb", "8_to_9kb"], fill_value=0)
    fractions = counts.div(counts.sum(axis=1), axis=0).fillna(0)
    fig, ax = plt.subplots(figsize=(max(7, len(clusters) * 0.9), 4.5))
    bottom = np.zeros(len(clusters))
    for column, color in zip(fractions.columns, ["#59A14F", "#F28E2B", "#E15759"]):
        values = fractions[column].to_numpy()
        ax.bar(x, values, bottom=bottom, label=column, color=color)
        bottom += values
    ax.set_xticks(x, labels, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Fraction of units")
    ax.set_ylim(0, 1)
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(output / "cluster_unit_length_bins.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def write_dataset_outputs(
    root: Path,
    dataset: str,
    units: pd.DataFrame,
    qc: dict[str, object],
    bootstrap: int,
    seed: int,
    plots: bool,
) -> dict[str, pd.DataFrame]:
    output = root / "results" / dataset / "cluster_analysis"
    figures = output / "figures"
    output.mkdir(parents=True, exist_ok=True)
    arrays = array_cluster_table(units)
    summary = cluster_summary(units, arrays)
    length_bins = length_bin_table(units)
    families, top = component_tables(units)
    intervals = bootstrap_array_intervals(arrays, bootstrap, seed)
    tests = descriptive_tests(arrays)
    species = species_cluster_summary(units)

    units.to_csv(output / "unit_metrics.tsv.gz", sep="\t", index=False)
    arrays.to_csv(output / "array_cluster_metrics.tsv.gz", sep="\t", index=False)
    summary.to_csv(output / "cluster_summary.tsv", sep="\t", index=False)
    allele_columns = [
        "cluster",
        "cluster_label",
        "units",
        "arrays",
        "subfamily_1_count",
        "subfamily_2_count",
        "subfamily_3_count",
        "classified_acro1_total",
        "raw_acro1_anchor_total",
        "unclassified_raw_acro1_anchor_total",
        "fraction_raw_acro1_anchors_classified",
        "fraction_raw_acro1_anchors_classified_array_weighted",
        "arrays_with_classified_acro1",
        "fraction_subfamily_1_pooled",
        "fraction_subfamily_2_pooled",
        "fraction_subfamily_3_pooled",
        "fraction_subfamily_1_array_weighted",
        "fraction_subfamily_2_array_weighted",
        "fraction_subfamily_3_array_weighted",
        "mean_acro1_copy_number",
        "median_acro1_copy_number",
        "array_weighted_mean_acro1_copy_number",
        "mean_raw_acro1_anchor_count",
        "median_raw_acro1_anchor_count",
        "mean_raw_acro1_anchors_per_kb",
        "units_without_raw_acro1_anchors",
        "fraction_units_without_raw_acro1_anchors",
        "mean_classified_acro1_copy_number",
        "median_classified_acro1_copy_number",
        "array_weighted_mean_classified_acro1_copy_number",
        "units_without_classified_acro1",
        "fraction_units_without_classified_acro1",
    ]
    summary[allele_columns].to_csv(
        output / "cluster_acro1_allele_composition.tsv", sep="\t", index=False
    )
    length_bins.to_csv(output / "cluster_unit_length_bins.tsv", sep="\t", index=False)
    families.to_csv(output / "cluster_component_bp.tsv", sep="\t", index=False)
    top.to_csv(output / "cluster_top_components.tsv", sep="\t", index=False)
    intervals.to_csv(output / "cluster_array_bootstrap_intervals.tsv", sep="\t", index=False)
    tests.to_csv(output / "cluster_descriptive_tests.tsv", sep="\t", index=False)
    species.to_csv(output / "species_cluster_summary.tsv", sep="\t", index=False)
    with (output / "analysis_qc.json").open("w") as handle:
        json.dump(qc, handle, indent=2, sort_keys=True)
        handle.write("\n")
    if plots:
        make_plots(units, arrays, summary, figures)
    return {"units": units, "summary": summary, "species": species}


def build_crosswalks(
    analyses: Mapping[str, dict[str, pd.DataFrame]], joint_id: str
) -> tuple[pd.DataFrame, pd.DataFrame]:
    columns = [
        "dataset_id",
        "individual_cluster",
        "individual_cluster_label",
        "joint_cluster",
        "joint_cluster_label",
        "intersection_units",
        "individual_cluster_units_shared",
        "joint_cluster_units_shared",
        "purity",
        "recall",
        "jaccard",
        "hungarian_match",
        "shared_units_total",
        "adjusted_rand_index_all_shared",
    ]
    overlap_columns = [
        "dataset_id",
        "individual_units",
        "joint_units",
        "shared_exact_coordinate_units",
        "fraction_individual_units_shared",
        "fraction_joint_units_matching_dataset",
        "adjusted_rand_index",
    ]
    if joint_id not in analyses:
        return pd.DataFrame(columns=columns), pd.DataFrame(columns=overlap_columns)
    joint = analyses[joint_id]["units"]
    rows: list[dict[str, object]] = []
    overlaps: list[dict[str, object]] = []
    for dataset, tables in analyses.items():
        if dataset == joint_id:
            continue
        individual = tables["units"]
        left = individual[["unit_key", "cluster", "cluster_label"]].rename(
            columns={"cluster": "individual_cluster", "cluster_label": "individual_cluster_label"}
        )
        right = joint[["unit_key", "cluster", "cluster_label", "species_id"]].rename(
            columns={"cluster": "joint_cluster", "cluster_label": "joint_cluster_label"}
        )
        shared = left.merge(right, on="unit_key", how="inner", validate="one_to_one")
        if shared.empty:
            overlaps.append(
                {
                    "dataset_id": dataset,
                    "individual_units": len(individual),
                    "joint_units": len(joint),
                    "shared_exact_coordinate_units": 0,
                    "fraction_individual_units_shared": 0.0,
                    "fraction_joint_units_matching_dataset": 0.0,
                    "adjusted_rand_index": np.nan,
                }
            )
            continue
        ari = adjusted_rand_score(shared.individual_cluster, shared.joint_cluster)
        contingency = pd.crosstab(shared.individual_cluster, shared.joint_cluster)
        row_index, col_index = linear_sum_assignment(-contingency.to_numpy())
        matched = {
            (int(contingency.index[i]), int(contingency.columns[j]))
            for i, j in zip(row_index, col_index)
        }
        individual_sizes = shared.individual_cluster.value_counts()
        joint_sizes = shared.joint_cluster.value_counts()
        individual_labels = shared.groupby("individual_cluster").individual_cluster_label.first()
        joint_labels = shared.groupby("joint_cluster").joint_cluster_label.first()
        for a in contingency.index:
            for b in contingency.columns:
                intersection = int(contingency.loc[a, b])
                if not intersection:
                    continue
                a_size, b_size = int(individual_sizes[a]), int(joint_sizes[b])
                rows.append(
                    {
                        "dataset_id": dataset,
                        "individual_cluster": int(a),
                        "individual_cluster_label": individual_labels[a],
                        "joint_cluster": int(b),
                        "joint_cluster_label": joint_labels[b],
                        "intersection_units": intersection,
                        "individual_cluster_units_shared": a_size,
                        "joint_cluster_units_shared": b_size,
                        "purity": intersection / a_size,
                        "recall": intersection / b_size,
                        "jaccard": intersection / (a_size + b_size - intersection),
                        "hungarian_match": (int(a), int(b)) in matched,
                        "shared_units_total": len(shared),
                        "adjusted_rand_index_all_shared": ari,
                    }
                )
        joint_dataset_units = int(joint.species_id.eq(dataset).sum())
        if dataset == "hprc":
            joint_dataset_units = int(joint.species_id.eq("hprc").sum())
        overlaps.append(
            {
                "dataset_id": dataset,
                "individual_units": len(individual),
                "joint_units": len(joint),
                "shared_exact_coordinate_units": len(shared),
                "fraction_individual_units_shared": len(shared) / len(individual),
                "fraction_joint_units_matching_dataset": (
                    len(shared) / joint_dataset_units if joint_dataset_units else np.nan
                ),
                "adjusted_rand_index": ari,
            }
        )
    return pd.DataFrame(rows, columns=columns), pd.DataFrame(overlaps, columns=overlap_columns)


def main() -> None:
    args = parse_args()
    root = args.root.resolve()
    manifest_path = args.manifest.resolve()
    manifest, id_col = load_manifest(manifest_path, args.datasets)
    species_map = load_species_map(manifest, id_col, manifest_path.parent, args.primate_map)
    analyses: dict[str, dict[str, pd.DataFrame]] = {}
    combined_summaries: list[pd.DataFrame] = []
    combined_species: list[pd.DataFrame] = []
    dataset_input_provenance: list[dict[str, object]] = []

    for row in manifest.to_dict("records"):
        dataset = str(row[id_col])
        prefix = output_prefix(root, dataset, row)
        assignments = selected_assignments_path(root, dataset, row)
        raw_anchor_path = resolve_path(row.get("acro1_monomers", ""), manifest_path.parent)
        native_clusters = Path(f"{prefix}.clusters.tsv")
        native_annotations = Path(f"{prefix}.annotations.bed")
        classified_anchors = Path(f"{prefix}.acro1_subfamilies.bed")
        provenance_files = {
            "selected_assignments": assignments,
            "native_clusters": native_clusters,
            "native_annotations": native_annotations,
            "classified_acro1_anchors": classified_anchors,
            "raw_acro1_anchor_input": raw_anchor_path,
        }
        dataset_input_provenance.append(
            {
                "dataset_id": dataset,
                "ntriage_prefix": str(prefix),
                "files": {
                    name: (
                        {
                            "path": str(path),
                            "sha256": sha256_file(path),
                            "bytes": path.stat().st_size,
                        }
                        if path is not None and path.exists()
                        else None
                    )
                    for name, path in provenance_files.items()
                },
            }
        )
        print(f"[{dataset}] loading {prefix}", file=sys.stderr, flush=True)
        units = load_units(prefix, assignments, dataset)
        units, qc = build_unit_metrics(
            prefix,
            units,
            dataset,
            species_map,
            raw_anchor_path=raw_anchor_path,
            joint_id=args.joint_id,
        )
        tables = write_dataset_outputs(
            root,
            dataset,
            units,
            qc,
            args.bootstrap,
            args.seed,
            not args.skip_plots,
        )
        analyses[dataset] = tables
        summary = tables["summary"].copy()
        summary.insert(0, "dataset_id", dataset)
        combined_summaries.append(summary)
        species = tables["species"].copy()
        species.insert(0, "dataset_id", dataset)
        combined_species.append(species)

    global_output = resolve_path(args.combined_output_dir, root)
    if global_output is None:  # argparse always supplies a value.
        raise ValueError("--combined-output-dir cannot be empty")
    global_output.mkdir(parents=True, exist_ok=True)
    pd.concat(combined_summaries, ignore_index=True).to_csv(
        global_output / "all_cluster_summaries.tsv", sep="\t", index=False
    )
    pd.concat(combined_species, ignore_index=True).to_csv(
        global_output / "all_species_cluster_summaries.tsv", sep="\t", index=False
    )
    crosswalk, overlap = build_crosswalks(analyses, args.joint_id)
    crosswalk.to_csv(
        global_output / "individual_vs_joint_cluster_crosswalk.tsv", sep="\t", index=False
    )
    overlap.to_csv(
        global_output / "individual_vs_joint_unit_overlap.tsv", sep="\t", index=False
    )
    parameters = {
        "manifest": str(manifest_path),
        "manifest_sha256": sha256_file(manifest_path),
        "analysis_script": str(SCRIPT),
        "analysis_script_sha256": sha256_file(SCRIPT),
        "command": [sys.executable, *sys.argv],
        "datasets": list(analyses),
        "dataset_input_provenance": dataset_input_provenance,
        "required_dataset_set_present": REQUIRED_DATASETS.issubset(analyses),
        "joint_id": args.joint_id,
        "combined_output_dir": str(global_output),
        "bootstrap_repetitions": args.bootstrap,
        "seed": args.seed,
        "primate_species_map": str(args.primate_map),
        "unit_length_bins_bp": [6000, 7000, 8000, 9000],
        "interpretation_scope": SCOPE_NOTE,
    }
    with (global_output / "analysis_parameters.json").open("w") as handle:
        json.dump(parameters, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(f"Wrote cluster analyses for {len(analyses)} datasets", file=sys.stderr)


if __name__ == "__main__":
    main()
