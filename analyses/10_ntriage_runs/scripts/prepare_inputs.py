#!/usr/bin/env python3
"""Build and validate inputs for the November comparative NTRIAGE experiment.

The six non-HPRC cohorts share a synthetic multi-FASTA. Their BED records are
split by species and their genomic ACRO1 monomer anchors are translated onto
the synthetic PF contigs. HPRC uses core (flank-excluded) intervals. Two joint
inputs combine either the six primates alone or all primate and HPRC records.
All BED coordinates are 0-based, half-open.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Iterator, Mapping, Sequence


DEFAULT_ROOT = Path(os.environ.get("ACRO1_NTRIAGE_ANALYSIS_ROOT", Path(__file__).resolve().parents[1]))
ACREAGE = Path(os.environ.get("ACREAGE_ROOT", DEFAULT_ROOT / "inputs"))
SIX_FASTA = ACREAGE / "primate_full/six_primate_arrays.fasta"
SIX_BED = ACREAGE / "primate_full/six_primate_arrays.bed"
SIX_MAPPING = ACREAGE / "primate_full/six_primate_array_mapping.tsv"
HPRC_INPUTS = Path(os.environ.get("HPRC_NTRIAGE_INPUTS", DEFAULT_ROOT / "hprc_inputs"))
HPRC_FASTA = HPRC_INPUTS / "telomeric_2794.fasta"
HPRC_BED = HPRC_INPUTS / "telomeric_2794.core.bed"
HPRC_ANCHORS = HPRC_INPUTS / "acro1_monomer_anchors.bed"
HPRC_AUDIT = HPRC_INPUTS / "input_audit.tsv"

DATASET_SPECS: tuple[tuple[str, str, str, str], ...] = (
    ("hg002", "HG002", "hg002", "HG002"),
    ("gorilla", "Gorilla", "mGorGor1", "mGorGor1"),
    ("chimpanzee", "Chimpanzee", "mPanTro3", "mPanTro3"),
    ("bonobo", "Bonobo", "mPanPan1", "mPanPan1"),
    ("sumatran_orangutan", "Sumatran orangutan", "mPonAbe1", "mPonAbe1"),
    ("bornean_orangutan", "Bornean orangutan", "mPonPyg2", "mPonPyg2"),
)

METADATA_FIELDS = (
    "array_id", "chrom", "source_name", "source_dataset", "cohort", "species_id", "species_name",
    "assembly", "individual_id", "haplotype", "analysis_chrom", "analysis_start",
    "analysis_end", "analysis_length_bp", "original_chrom", "original_start",
    "original_end", "original_bed_name",
)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def require_file(path: Path) -> None:
    if not path.is_file():
        raise FileNotFoundError(path)


def sha256_file(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Iterable[Mapping[str, object]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def iter_bed(path: Path) -> Iterator[list[str]]:
    with path.open() as handle:
        for line_no, raw in enumerate(handle, start=1):
            if not raw.strip() or raw.startswith(("#", "track", "browser")):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError(f"{path}:{line_no}: expected at least 3 BED columns")
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as error:
                raise ValueError(f"{path}:{line_no}: non-integer coordinates") from error
            if start < 0 or end <= start:
                raise ValueError(f"{path}:{line_no}: invalid interval {start}-{end}")
            yield fields


def write_bed(path: Path, records: Iterable[Sequence[object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for fields in records:
            handle.write("\t".join(map(str, fields)) + "\n")


def load_fai(fasta: Path) -> dict[str, int]:
    fai = Path(str(fasta) + ".fai")
    require_file(fai)
    lengths: dict[str, int] = {}
    with fai.open() as handle:
        for line_no, raw in enumerate(handle, start=1):
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise ValueError(f"{fai}:{line_no}: malformed index")
            if fields[0] in lengths:
                raise ValueError(f"{fai}: duplicate sequence {fields[0]}")
            lengths[fields[0]] = int(fields[1])
    return lengths


def ensure_fai(fasta: Path, samtools: str, force: bool = False) -> None:
    fai = Path(str(fasta) + ".fai")
    if force and fai.exists():
        fai.unlink()
    if not fai.exists():
        subprocess.run([samtools, "faidx", str(fasta)], check=True)


def validate_regions(path: Path, fai: Mapping[str, int]) -> tuple[int, int, int, int]:
    seen: set[str] = set()
    count = total = maximum = 0
    minimum: int | None = None
    for fields in iter_bed(path):
        chrom, start, end = fields[0], int(fields[1]), int(fields[2])
        if chrom not in fai:
            raise ValueError(f"{path}: {chrom!r} absent from FASTA index")
        if end > fai[chrom]:
            raise ValueError(f"{path}: {chrom}:{start}-{end} exceeds {fai[chrom]}")
        key = fields[3] if len(fields) >= 4 else f"{chrom}:{start}-{end}"
        if key in seen:
            raise ValueError(f"{path}: duplicate array ID {key!r}")
        seen.add(key)
        length = end - start
        count += 1
        total += length
        minimum = length if minimum is None else min(minimum, length)
        maximum = max(maximum, length)
    if not count:
        raise ValueError(f"{path}: no intervals")
    return count, total, minimum or 0, maximum


def anchor_qc(path: Path, fai: Mapping[str, int], regions: Path) -> dict[str, int]:
    interval_by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for fields in iter_bed(regions):
        interval_by_chrom[fields[0]].append((int(fields[1]), int(fields[2])))
    total = in_regions = 0
    anchor_chroms: set[str] = set()
    for fields in iter_bed(path):
        chrom, start, end = fields[0], int(fields[1]), int(fields[2])
        if chrom not in fai or end > fai[chrom]:
            raise ValueError(f"{path}: invalid anchor {chrom}:{start}-{end}")
        if len(fields) < 4 or "_monomer_" not in fields[3]:
            raise ValueError(f"{path}: anchor lacks *_monomer_* name: {fields}")
        total += 1
        anchor_chroms.add(chrom)
        if any(start >= left and end <= right for left, right in interval_by_chrom.get(chrom, [])):
            in_regions += 1
    return {
        "n_anchor_records": total,
        "n_anchor_chroms": len(anchor_chroms),
        "n_anchors_in_analysis_intervals": in_regions,
    }


def replace_or_keep(path: Path, force: bool, writer) -> None:
    if path.exists() and not force:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    try:
        writer(temporary)
        temporary.replace(path)
    finally:
        if temporary.exists():
            temporary.unlink()


def copy_fasta_pair(sources: Sequence[Path], destination: Path) -> None:
    with destination.open("wb") as out_handle:
        for source in sources:
            last_byte = b"\n"
            with source.open("rb") as in_handle:
                while True:
                    chunk = in_handle.read(8 * 1024 * 1024)
                    if not chunk:
                        break
                    out_handle.write(chunk)
                    last_byte = chunk[-1:]
            if last_byte not in (b"\n", b"\r"):
                out_handle.write(b"\n")


def primate_anchor_source(source_tag: str) -> Path:
    return (
        ACREAGE / "primate_acro_bedfiles" / source_tag
        / f"acro_monomers_full_{source_tag}.P.bed"
    )


def remap_primate_anchors(
    source: Path,
    destination: Path,
    mapping_by_name: Mapping[str, Mapping[str, str]],
    species_id: str,
) -> dict[str, int]:
    records: list[tuple[str, int, int, str, str, str]] = []
    excluded_outside = excluded_wrong_species = 0
    for fields in iter_bed(source):
        if len(fields) < 4:
            raise ValueError(f"{source}: expected named anchors")
        anchor_name = fields[3]
        array_name, separator, suffix = anchor_name.rpartition("_monomer_")
        if not separator or array_name not in mapping_by_name:
            raise ValueError(f"{source}: cannot map anchor {anchor_name!r}")
        row = mapping_by_name[array_name]
        if row["species_id"] != species_id:
            excluded_wrong_species += 1
            continue
        original_start = int(row["original_start"])
        input_length = int(row["input_length"])
        local_start = int(fields[1]) - original_start
        local_end = int(fields[2]) - original_start
        if local_start < 0 or local_end > input_length:
            excluded_outside += 1
            continue
        synthetic = row["synthetic_chrom"]
        records.append(
            (
                synthetic, local_start, local_end, f"{synthetic}_monomer_{suffix}",
                fields[4] if len(fields) >= 5 else ".",
                fields[5] if len(fields) >= 6 else ".",
            )
        )
    if excluded_wrong_species:
        raise ValueError(f"{source}: {excluded_wrong_species} anchors map to another species")
    records.sort(key=lambda item: (item[0], item[1], item[2], item[3]))
    write_bed(destination, records)
    return {
        "source_anchor_records": len(records) + excluded_outside,
        "retained_anchor_records": len(records),
        "excluded_outside_array": excluded_outside,
    }


def concatenate_text_files(sources: Sequence[Path], destination: Path) -> None:
    with destination.open("w") as out_handle:
        for source in sources:
            with source.open() as in_handle:
                for raw in in_handle:
                    if raw.strip() and not raw.startswith(("#", "track", "browser")):
                        out_handle.write(raw if raw.endswith("\n") else raw + "\n")


def filter_anchors_to_fasta(source: Path, destination: Path, fai: Mapping[str, int]) -> dict[str, int]:
    """Exclude full-length anchor hits that cross a FASTA record boundary."""
    retained: list[list[str]] = []
    excluded = 0
    total = 0
    for fields in iter_bed(source):
        total += 1
        chrom, start, end = fields[0], int(fields[1]), int(fields[2])
        if chrom not in fai or start < 0 or end > fai[chrom]:
            excluded += 1
            continue
        retained.append(fields)
    write_bed(destination, retained)
    return {
        "source_anchor_records": total,
        "retained_anchor_records": len(retained),
        "excluded_outside_array": excluded,
    }


def build_metadata(
    mapping_rows: Sequence[Mapping[str, str]],
    hprc_bed: Path | None = None,
    hprc_audit: Path | None = None,
) -> tuple[dict[str, list[dict[str, object]]], list[dict[str, object]]]:
    by_dataset: dict[str, list[dict[str, object]]] = defaultdict(list)
    species_to_dataset = {species: dataset for dataset, _, species, _ in DATASET_SPECS}
    for row in mapping_rows:
        dataset = species_to_dataset[row["species_id"]]
        by_dataset[dataset].append(
            {
                "array_id": row["synthetic_chrom"],
                "chrom": row["synthetic_chrom"],
                "source_name": row["synthetic_chrom"],
                "source_dataset": dataset,
                "cohort": "primate",
                "species_id": row["species_id"],
                "species_name": row["species_name"],
                "assembly": row["assembly"],
                "individual_id": "HG002" if dataset == "hg002" else row["species_id"],
                "haplotype": "",
                "analysis_chrom": row["synthetic_chrom"],
                "analysis_start": 0,
                "analysis_end": row["input_length"],
                "analysis_length_bp": row["input_length"],
                "original_chrom": row["original_chrom"],
                "original_start": row["original_start"],
                "original_end": row["original_end"],
                "original_bed_name": row["original_bed_name"],
            }
        )

    if (hprc_bed is None) != (hprc_audit is None):
        raise ValueError("hprc_bed and hprc_audit must be supplied together")
    if hprc_bed is not None and hprc_audit is not None:
        audit_by_array = {row["array_id"]: row for row in read_tsv(hprc_audit)}
        for fields in iter_bed(hprc_bed):
            array_id, start, end = fields[0], int(fields[1]), int(fields[2])
            audit = audit_by_array.get(array_id)
            if audit is None:
                raise ValueError(f"HPRC interval {array_id!r} absent from {hprc_audit}")
            by_dataset["hprc"].append(
                {
                    "array_id": array_id,
                    "chrom": array_id,
                    "source_name": array_id,
                    "source_dataset": "hprc",
                    "cohort": "hprc",
                    "species_id": "hprc",
                    "species_name": "human",
                    "assembly": "HPRC_release2",
                    "individual_id": audit["individual_id"],
                    "haplotype": audit["haplotype"],
                    "analysis_chrom": array_id,
                    "analysis_start": start,
                    "analysis_end": end,
                    "analysis_length_bp": end - start,
                    "original_chrom": "",
                    "original_start": "",
                    "original_end": "",
                    "original_bed_name": array_id,
                }
            )
    all_rows: list[dict[str, object]] = []
    for dataset, _, _, _ in DATASET_SPECS:
        all_rows.extend(by_dataset[dataset])
    if hprc_bed is not None:
        all_rows.extend(by_dataset["hprc"])
    return by_dataset, all_rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare the six primate, HPRC, and joint NTRIAGE inputs."
    )
    parser.add_argument("--analysis-root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--six-fasta", type=Path, default=SIX_FASTA)
    parser.add_argument("--six-bed", type=Path, default=SIX_BED)
    parser.add_argument("--six-mapping", type=Path, default=SIX_MAPPING)
    parser.add_argument("--hprc-fasta", type=Path, default=HPRC_FASTA)
    parser.add_argument("--hprc-bed", type=Path, default=HPRC_BED)
    parser.add_argument("--hprc-anchors", type=Path, default=HPRC_ANCHORS)
    parser.add_argument("--hprc-audit", type=Path, default=HPRC_AUDIT)
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument(
        "--six-primates-only",
        action="store_true",
        help=(
            "prepare only the six individual primate cohorts and their joint run; "
            "do not require or read HPRC inputs, materialize joint_all, or replace "
            "the full manifest/shared provenance"
        ),
    )
    parser.add_argument(
        "--force", action="store_true", help="replace deterministic derived BED/FASTA files"
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    root = args.analysis_root.resolve()
    inputs = root / "inputs"
    regions_dir = inputs / "regions"
    anchors_dir = inputs / "anchors"
    metadata_dir = inputs / "metadata"
    shared_provenance_dir = root / "provenance"
    provenance_dir = (
        shared_provenance_dir / "six_primates_no_hprc"
        if args.six_primates_only
        else shared_provenance_dir
    )
    for directory in (regions_dir, anchors_dir, metadata_dir, provenance_dir):
        directory.mkdir(parents=True, exist_ok=True)

    required: tuple[Path, ...] = (
        args.six_fasta, Path(str(args.six_fasta) + ".fai"), args.six_bed,
        args.six_mapping,
    )
    if not args.six_primates_only:
        required += (
            args.hprc_fasta,
            Path(str(args.hprc_fasta) + ".fai"),
            args.hprc_bed,
            args.hprc_anchors,
            args.hprc_audit,
        )
    for path in required:
        require_file(path)
    mapping_rows = read_tsv(args.six_mapping)
    expected_mapping_fields = {
        "synthetic_chrom", "species_id", "species_name", "assembly", "original_chrom",
        "original_start", "original_end", "input_length", "original_bed_name",
    }
    if not mapping_rows or not expected_mapping_fields.issubset(mapping_rows[0]):
        raise ValueError(f"{args.six_mapping}: missing required mapping columns")
    mapping_by_name = {row["original_bed_name"]: row for row in mapping_rows}
    mapping_by_chrom = {row["synthetic_chrom"]: row for row in mapping_rows}
    if len(mapping_by_name) != len(mapping_rows) or len(mapping_by_chrom) != len(mapping_rows):
        raise ValueError("Primate mapping identifiers are not unique")

    full_bed_by_chrom = {fields[0]: fields for fields in iter_bed(args.six_bed)}
    if set(full_bed_by_chrom) != set(mapping_by_chrom):
        missing = set(mapping_by_chrom) - set(full_bed_by_chrom)
        extra = set(full_bed_by_chrom) - set(mapping_by_chrom)
        raise ValueError(f"Six-primate BED/mapping mismatch; missing={missing}, extra={extra}")

    by_dataset_metadata, all_metadata = build_metadata(
        mapping_rows,
        None if args.six_primates_only else args.hprc_bed,
        None if args.six_primates_only else args.hprc_audit,
    )
    primate_metadata = [
        row
        for dataset, _, _, _ in DATASET_SPECS
        for row in by_dataset_metadata[dataset]
    ]
    display_names = {dataset: display for dataset, display, _, _ in DATASET_SPECS}
    display_names.update(
        {
            "joint_primates": "Six primates (joint; no HPRC)",
            "hprc": "HPRC",
            "joint_all": "All cohorts (joint)",
        }
    )

    anchor_audit_rows: list[dict[str, object]] = []
    primate_anchor_paths: list[Path] = []
    dataset_paths: dict[str, dict[str, object]] = {}
    for dataset, _, species_id, source_tag in DATASET_SPECS:
        bed_path = regions_dir / f"{dataset}.bed"
        anchor_path = anchors_dir / f"{dataset}.acro1_monomers.bed"
        metadata_path = metadata_dir / f"{dataset}.arrays.tsv"
        rows_for_species = [row for row in mapping_rows if row["species_id"] == species_id]
        expected_chroms = {row["synthetic_chrom"] for row in rows_for_species}

        def bed_writer(path: Path, chroms=expected_chroms) -> None:
            write_bed(path, (full_bed_by_chrom[chrom] for chrom in sorted(chroms)))

        replace_or_keep(bed_path, args.force, bed_writer)
        source_anchors = primate_anchor_source(source_tag)
        require_file(source_anchors)
        remap_stats: dict[str, int] = {}

        def anchor_writer(path: Path) -> None:
            remap_stats.update(
                remap_primate_anchors(source_anchors, path, mapping_by_name, species_id)
            )

        replace_or_keep(anchor_path, args.force, anchor_writer)
        if not remap_stats:
            source_count = sum(1 for _ in iter_bed(source_anchors))
            retained_count = sum(1 for _ in iter_bed(anchor_path))
            remap_stats = {
                "source_anchor_records": source_count,
                "retained_anchor_records": retained_count,
                "excluded_outside_array": source_count - retained_count,
            }
        anchor_audit_rows.append({"dataset_id": dataset, **remap_stats})
        write_tsv(metadata_path, by_dataset_metadata[dataset], METADATA_FIELDS)
        primate_anchor_paths.append(anchor_path)
        dataset_paths[dataset] = {
            "fasta": args.six_fasta.resolve(),
            "bed": bed_path.resolve(),
            "anchors": anchor_path.resolve(),
            "metadata": metadata_path.resolve(),
            "cohort": "primate",
            "source_bed": args.six_bed.resolve(),
            "source_anchors": source_anchors.resolve(),
        }

    joint_primates_bed = regions_dir / "joint_primates.bed"
    replace_or_keep(
        joint_primates_bed, args.force,
        lambda path: shutil.copyfile(args.six_bed, path),
    )
    joint_primates_anchors = anchors_dir / "joint_primates.acro1_monomers.bed"
    replace_or_keep(
        joint_primates_anchors, args.force,
        lambda path: concatenate_text_files(tuple(primate_anchor_paths), path),
    )
    joint_primates_metadata = metadata_dir / "joint_primates.arrays.tsv"
    write_tsv(joint_primates_metadata, primate_metadata, METADATA_FIELDS)
    dataset_paths["joint_primates"] = {
        "fasta": args.six_fasta.resolve(),
        "bed": joint_primates_bed.resolve(),
        "anchors": joint_primates_anchors.resolve(),
        "metadata": joint_primates_metadata.resolve(),
        "cohort": "joint",
        "source_bed": args.six_bed.resolve(),
        "source_anchors": "six remapped primate anchor sets",
    }

    if not args.six_primates_only:
        hprc_bed_local = regions_dir / "hprc.bed"
        replace_or_keep(
            hprc_bed_local,
            args.force,
            lambda path: shutil.copyfile(args.hprc_bed, path),
        )
        hprc_anchor_local = anchors_dir / "hprc.acro1_monomers.bed"
        hprc_anchor_stats: dict[str, int] = {}

        def hprc_anchor_writer(path: Path) -> None:
            hprc_anchor_stats.update(
                filter_anchors_to_fasta(
                    args.hprc_anchors, path, load_fai(args.hprc_fasta)
                )
            )

        replace_or_keep(hprc_anchor_local, args.force, hprc_anchor_writer)
        if not hprc_anchor_stats:
            source_count = sum(1 for _ in iter_bed(args.hprc_anchors))
            retained_count = sum(1 for _ in iter_bed(hprc_anchor_local))
            hprc_anchor_stats = {
                "source_anchor_records": source_count,
                "retained_anchor_records": retained_count,
                "excluded_outside_array": source_count - retained_count,
            }
        anchor_audit_rows.append({"dataset_id": "hprc", **hprc_anchor_stats})
        hprc_metadata_path = metadata_dir / "hprc.arrays.tsv"
        write_tsv(hprc_metadata_path, by_dataset_metadata["hprc"], METADATA_FIELDS)
        dataset_paths["hprc"] = {
            "fasta": args.hprc_fasta.resolve(),
            "bed": hprc_bed_local.resolve(),
            "anchors": hprc_anchor_local.resolve(),
            "metadata": hprc_metadata_path.resolve(),
            "cohort": "hprc",
            "source_bed": args.hprc_bed.resolve(),
            "source_anchors": args.hprc_anchors.resolve(),
        }

        joint_fasta = inputs / "joint_all.fasta"
        replace_or_keep(
            joint_fasta,
            args.force,
            lambda path: copy_fasta_pair((args.six_fasta, args.hprc_fasta), path),
        )
        ensure_fai(joint_fasta, args.samtools, force=args.force)
        joint_bed = regions_dir / "joint_all.bed"
        replace_or_keep(
            joint_bed,
            args.force,
            lambda path: concatenate_text_files((args.six_bed, hprc_bed_local), path),
        )
        joint_anchors = anchors_dir / "joint_all.acro1_monomers.bed"
        replace_or_keep(
            joint_anchors,
            args.force,
            lambda path: concatenate_text_files(
                (*primate_anchor_paths, hprc_anchor_local), path
            ),
        )
        joint_metadata = metadata_dir / "joint_all.arrays.tsv"
        write_tsv(joint_metadata, all_metadata, METADATA_FIELDS)
        dataset_paths["joint_all"] = {
            "fasta": joint_fasta.resolve(),
            "bed": joint_bed.resolve(),
            "anchors": joint_anchors.resolve(),
            "metadata": joint_metadata.resolve(),
            "cohort": "joint",
            "source_bed": "six_primate_arrays.bed + HPRC core BED",
            "source_anchors": "remapped primate anchors + HPRC anchors",
        }

        write_tsv(metadata_dir / "all_source_arrays.tsv", all_metadata, METADATA_FIELDS)
    write_tsv(
        provenance_dir / "primate_anchor_remapping_qc.tsv", anchor_audit_rows,
        ("dataset_id", "source_anchor_records", "retained_anchor_records", "excluded_outside_array"),
    )

    checksum_cache: dict[Path, str] = {}

    def cached_sha(path: Path) -> str:
        resolved = path.resolve()
        if resolved not in checksum_cache:
            checksum_cache[resolved] = sha256_file(resolved)
        return checksum_cache[resolved]

    manifest_rows: list[dict[str, object]] = []
    input_qc_rows: list[dict[str, object]] = []
    ordered_datasets = [dataset for dataset, _, _, _ in DATASET_SPECS] + [
        "joint_primates"
    ]
    if not args.six_primates_only:
        ordered_datasets.extend(["hprc", "joint_all"])
    for dataset in ordered_datasets:
        paths = dataset_paths[dataset]
        fasta = Path(paths["fasta"])
        bed = Path(paths["bed"])
        anchors = Path(paths["anchors"])
        metadata = Path(paths["metadata"])
        fai = load_fai(fasta)
        n_arrays, analyzed_bp, minimum, maximum = validate_regions(bed, fai)
        anchor_stats = anchor_qc(anchors, fai, bed)
        prefix = (root / "runs" / dataset / "native" / dataset).resolve()
        manifest_rows.append(
            {
                "dataset_id": dataset,
                "display_name": display_names[dataset],
                "cohort": paths["cohort"],
                "fasta": fasta,
                "bed": bed,
                "acro1_monomers": anchors,
                "array_metadata": metadata,
                "array_dataset_map": metadata,
                "ntriage_prefix": prefix,
                "features_path": str(prefix) + ".features.tsv",
                "arrays_path": str(prefix) + ".clusters.tsv",
                "n_arrays": n_arrays,
                "analyzed_bp": analyzed_bp,
                "n_anchor_records": anchor_stats["n_anchor_records"],
                "min_period_bp": 6000,
                "max_period_bp": 9000,
                "period_target_bp": 6650,
                "fasta_sha256": cached_sha(fasta),
                "bed_sha256": cached_sha(bed),
                "acro1_monomers_sha256": cached_sha(anchors),
                "source_bed": paths["source_bed"],
                "source_anchors": paths["source_anchors"],
            }
        )
        input_qc_rows.append(
            {
                "dataset_id": dataset,
                "n_fasta_records": len(fai),
                "n_arrays": n_arrays,
                "analyzed_bp": analyzed_bp,
                "min_array_bp": minimum,
                "max_array_bp": maximum,
                **anchor_stats,
            }
        )

    manifest_fields = (
        "dataset_id", "display_name", "cohort", "fasta", "bed", "acro1_monomers",
        "array_metadata", "array_dataset_map", "ntriage_prefix", "features_path",
        "arrays_path", "n_arrays", "analyzed_bp", "n_anchor_records", "min_period_bp",
        "max_period_bp", "period_target_bp", "fasta_sha256", "bed_sha256",
        "acro1_monomers_sha256", "source_bed", "source_anchors",
    )
    no_hprc_manifest = inputs / "six_primates_no_hprc_manifest.tsv"
    if args.six_primates_only:
        write_tsv(no_hprc_manifest, manifest_rows, manifest_fields)
    else:
        write_tsv(inputs / "dataset_manifest.tsv", manifest_rows, manifest_fields)
        no_hprc_ids = {
            "hg002",
            "gorilla",
            "chimpanzee",
            "bonobo",
            "sumatran_orangutan",
            "bornean_orangutan",
            "joint_primates",
        }
        write_tsv(
            no_hprc_manifest,
            (row for row in manifest_rows if row["dataset_id"] in no_hprc_ids),
            manifest_fields,
        )
    qc_fields = (
        "dataset_id", "n_fasta_records", "n_arrays", "analyzed_bp", "min_array_bp",
        "max_array_bp", "n_anchor_records", "n_anchor_chroms",
        "n_anchors_in_analysis_intervals",
    )
    write_tsv(provenance_dir / "input_qc.tsv", input_qc_rows, qc_fields)

    if args.six_primates_only:
        provenance = {
            "created_utc": utc_now(),
            "command": [sys.executable, *sys.argv],
            "coordinate_convention": "0-based half-open BED",
            "joint_primates_fasta_source": str(args.six_fasta.resolve()),
            "source_checksums": {
                "six_mapping_sha256": cached_sha(args.six_mapping),
            },
            "dataset_order": ordered_datasets,
            "notes": [
                "Primate monomer anchors crossing an extracted-array boundary are excluded, not clipped.",
                "joint_primates contains exactly the six primate cohorts and excludes HPRC.",
            ],
        }
    else:
        provenance = {
            "created_utc": utc_now(),
            "command": [sys.executable, *sys.argv],
            "coordinate_convention": "0-based half-open BED",
            "joint_primates_fasta_source": str(args.six_fasta.resolve()),
            "joint_fasta_sources": [
                str(args.six_fasta.resolve()),
                str(args.hprc_fasta.resolve()),
            ],
            "source_checksums": {
                "six_mapping_sha256": cached_sha(args.six_mapping),
                "hprc_input_audit_sha256": cached_sha(args.hprc_audit),
            },
            "dataset_order": ordered_datasets,
            "notes": [
                "Primate monomer anchors crossing an extracted-array boundary are excluded, not clipped.",
                "HPRC uses core intervals, excluding extraction flanks.",
                "joint_primates contains exactly the six primate cohorts and excludes HPRC.",
                "Joint FASTA record names are unique: PF* for primates and HPRC array IDs.",
            ],
        }
    (provenance_dir / "prepare_inputs.json").write_text(json.dumps(provenance, indent=2) + "\n")

    print(f"Prepared and validated {len(manifest_rows)} datasets under {root}")
    for row in input_qc_rows:
        print(
            f"  {row['dataset_id']}: {row['n_arrays']:,} arrays; "
            f"{row['analyzed_bp']:,} bp; {row['n_anchor_records']:,} ACRO1 anchors"
        )
    if not args.six_primates_only:
        print(f"Manifest: {inputs / 'dataset_manifest.tsv'}")
    print(f"Six-primate manifest: {no_hprc_manifest}")


if __name__ == "__main__":
    main()
