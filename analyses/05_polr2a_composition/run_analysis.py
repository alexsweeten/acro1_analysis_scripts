#!/usr/bin/env python3
"""Run the ACRO1--POLR2A CTD amino-acid-composition JSD analysis."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parent
INPUTS = ROOT / "inputs"
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
SOURCE = ROOT / "jsd_null_simulator.cpp"
BINARY = ROOT / "jsd_null_simulator"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run(command: list[str]) -> None:
    print("+", " ".join(command), flush=True)
    subprocess.run(command, check=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--replicates", type=int, default=1_000_000)
    parser.add_argument("--seed", type=int, default=20_261_103)
    parser.add_argument("--skip-simulation", action="store_true")
    parser.add_argument("--skip-plot", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.replicates < 1:
        raise SystemExit("--replicates must be positive")
    RESULTS.mkdir(exist_ok=True)
    FIGURES.mkdir(exist_ok=True)

    acro = INPUTS / "acro1_species_balanced_allele_consensus.tsv"
    polr2a = INPUTS / "polr2a_references.tsv"
    bases = INPUTS / "hg002_genome_base_composition.tsv"
    inputs = [acro, polr2a, bases]

    if not args.skip_simulation:
        compiler = shutil.which("c++") or shutil.which("clang++") or shutil.which("g++")
        if compiler is None:
            raise SystemExit("No C++17 compiler was found")
        run([
            compiler,
            "-O3",
            "-std=c++17",
            "-Wall",
            "-Wextra",
            "-pedantic",
            str(SOURCE),
            "-o",
            str(BINARY),
        ])
        run([
            str(BINARY),
            str(acro),
            str(polr2a),
            str(bases),
            str(args.replicates),
            str(args.seed),
            str(RESULTS / "jsd_null_summary.tsv"),
            str(RESULTS / "jsd_null_histogram.tsv"),
            str(RESULTS / "amino_acid_composition.tsv"),
            str(RESULTS / "observed_composition_summary.tsv"),
        ])

    if not args.skip_plot:
        rscript = shutil.which("Rscript")
        if rscript is None:
            raise SystemExit("Rscript is required to generate the figure")
        run([rscript, str(ROOT / "plot_results.R"), str(ROOT)])

    manifest = {
        "analysis": "ACRO1--POLR2A CTD codon-normalized amino-acid-composition JSD test",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "replicates": args.replicates,
        "seed": args.seed,
        "inputs": {str(path.relative_to(ROOT)): sha256(path) for path in inputs},
        "simulator_sha256": sha256(SOURCE),
        "notes": [
            "The three species-balanced ACRO1 allele consensuses are pooled with equal allele weight.",
            "Stop codons are excluded and the 20 sense-amino-acid frequencies are renormalized.",
            "The exact dinucleotide-count shuffle uses a randomized Euler trail and is not a uniform draw over all admissible sequences.",
            "Empirical P values are lower-tail plus-one estimates: (count[JSD_null <= JSD_observed] + 1)/(B + 1).",
        ],
    }
    with (RESULTS / "run_manifest.json").open("w") as handle:
        json.dump(manifest, handle, indent=2)
        handle.write("\n")


if __name__ == "__main__":
    main()
