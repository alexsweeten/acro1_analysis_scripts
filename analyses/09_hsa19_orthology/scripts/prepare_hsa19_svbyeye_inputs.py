#!/usr/bin/env python3
"""Prepare real sequence alignments for an SVbyEye HSA19 synteny plot."""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
WORK = ROOT / "svbyeye"
SEQUENCES = WORK / "sequences"

ALL_ORDER = [
    "Human",
    "Chimpanzee",
    "Bonobo",
    "Gorilla",
    "Sumatran_orangutan",
    "Bornean_orangutan",
    "Siamang_hap1",
    "Siamang_hap2",
    "Hoolock_gibbon",
    "Nomascus_gibbon",
    "Rhesus_macaque",
    "Marmoset",
    "Golden_lion_tamarin",
    "Squirrel_monkey",
]

REPRESENTATIVE_ORDER = [
    "Human",
    "Chimpanzee",
    "Gorilla",
    "Sumatran_orangutan",
    "Siamang_hap1",
    "Nomascus_gibbon",
    "Rhesus_macaque",
    "Marmoset",
    "Squirrel_monkey",
]

HALF_WINDOW = 700_000


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def fetch_sequence(fasta: str, chrom: str, start0: int, end0: int) -> str:
    # samtools regions are one-based, inclusive; analysis coordinates are
    # zero-based, half-open.
    region = f"{chrom}:{start0 + 1}-{end0}"
    completed = subprocess.run(
        ["samtools", "faidx", fasta, region],
        check=True,
        text=True,
        capture_output=True,
    )
    lines = completed.stdout.splitlines()
    if not lines or not lines[0].startswith(">"):
        raise RuntimeError(f"No FASTA sequence returned for {fasta} {region}")
    return "".join(lines[1:]).upper()


def write_fasta(path: Path, name: str, sequence: str) -> None:
    with path.open("w") as handle:
        handle.write(f">{name}\n")
        for index in range(0, len(sequence), 80):
            handle.write(sequence[index : index + 80] + "\n")


def extract_sequences() -> pd.DataFrame:
    WORK.mkdir(parents=True, exist_ok=True)
    SEQUENCES.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(RESULTS / "assembly_manifest.tsv", sep="\t").set_index("assembly_key")
    loci = pd.read_csv(RESULTS / "orthologous_locus_matrix.tsv", sep="\t").set_index("assembly_key")
    rows = []
    combined_records = []
    for key in ALL_ORDER:
        locus = loci.loc[key]
        source = manifest.loc[key]
        center = (int(locus.locus_start) + int(locus.locus_end)) // 2
        start0 = center - HALF_WINDOW
        end0 = center + HALF_WINDOW
        if start0 < 0:
            raise ValueError(f"Extraction would cross contig start for {key}")
        sequence = fetch_sequence(str(source.fasta), str(locus.target_chrom), start0, end0)
        if len(sequence) != 2 * HALF_WINDOW:
            raise ValueError(f"{key}: expected {2 * HALF_WINDOW} bp, observed {len(sequence)}")
        reverse_complemented = str(locus.orientation_to_human) == "-"
        if reverse_complemented:
            sequence = reverse_complement(sequence)
        fasta_path = SEQUENCES / f"{key}.fa"
        write_fasta(fasta_path, key, sequence)
        combined_records.append((key, sequence))
        rows.append(
            {
                "assembly_key": key,
                "source_fasta": source.fasta,
                "source_chrom": locus.target_chrom,
                "source_start_0based": start0,
                "source_end_0based": end0,
                "sequence_length": len(sequence),
                "orientation_to_human": locus.orientation_to_human,
                "reverse_complemented": reverse_complemented,
                "extracted_fasta": str(fasta_path),
            }
        )
    combined_path = WORK / "HSA19_all14_oriented_1.4Mb.fa"
    with combined_path.open("w") as handle:
        for name, sequence in combined_records:
            handle.write(f">{name}\n")
            for index in range(0, len(sequence), 80):
                handle.write(sequence[index : index + 80] + "\n")
    extracted = pd.DataFrame(rows)
    extracted.to_csv(WORK / "extracted_sequence_manifest.tsv", sep="\t", index=False)
    return extracted


def run_progressive_alignments(order: list[str], output_name: str, threads: int) -> Path:
    output = WORK / output_name
    command_rows = []
    with output.open("w") as paf_handle:
        for query, target in zip(order[:-1], order[1:]):
            command = [
                "minimap2",
                "-t",
                str(threads),
                "-x",
                "asm20",
                "-c",
                "--eqx",
                "--secondary=no",
                str(SEQUENCES / f"{target}.fa"),
                str(SEQUENCES / f"{query}.fa"),
            ]
            completed = subprocess.run(command, check=True, text=True, stdout=paf_handle, stderr=subprocess.PIPE)
            command_rows.append(
                {
                    "query": query,
                    "target": target,
                    "command": " ".join(command),
                    "stderr": completed.stderr.strip(),
                }
            )
    pd.DataFrame(command_rows).to_csv(
        WORK / output_name.replace(".paf", ".commands.tsv"), sep="\t", index=False
    )
    return output


def build_annotations() -> None:
    loci = pd.read_csv(RESULTS / "orthologous_locus_matrix.tsv", sep="\t").set_index("assembly_key")
    segments = pd.read_csv(RESULTS / "HSA19_synteny_plot_segments.tsv", sep="\t")
    rows = []
    for key in ALL_ORDER:
        span = int(loci.loc[key].locus_span_bp)
        rows.append(
            {
                "seqnames": key,
                "start": int(HALF_WINDOW - span / 2) + 1,
                "end": int(HALF_WINDOW + span / 2),
                "feature": "Inferred orthologous interval",
                "ID": key,
            }
        )
    feature_names = {
        "ACRO1_HMM_tract": "ACRO1 profile-HMM tract",
        "ACRE_component": "Non-ACRO ACRE component",
        "mobile_element": "Associated mobile-element component",
    }
    for _, row in segments.iterrows():
        if row.feature_type not in feature_names:
            continue
        start = max(1, int(float(row.oriented_start_bp) + HALF_WINDOW) + 1)
        end = min(2 * HALF_WINDOW, int(float(row.oriented_end_bp) + HALF_WINDOW))
        if end < start:
            continue
        rows.append(
            {
                "seqnames": row.assembly_key,
                "start": start,
                "end": end,
                "feature": feature_names[row.feature_type],
                "ID": row.assembly_key,
            }
        )
    pd.DataFrame(rows).to_csv(WORK / "HSA19_svbyeye_annotations.tsv", sep="\t", index=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()
    extract_sequences()
    build_annotations()
    all_paf = run_progressive_alignments(ALL_ORDER, "HSA19_all14_progressive.paf", args.threads)
    representative_paf = run_progressive_alignments(
        REPRESENTATIVE_ORDER, "HSA19_representative9_progressive.paf", args.threads
    )
    print(f"Prepared SVbyEye inputs in {WORK}")
    print(f"All-assembly PAF: {all_paf}")
    print(f"Representative PAF: {representative_paf}")


if __name__ == "__main__":
    main()
