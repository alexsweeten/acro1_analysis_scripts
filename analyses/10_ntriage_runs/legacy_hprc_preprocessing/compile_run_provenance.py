#!/usr/bin/env python3
"""Compile compact resource, cohort, and software provenance tables."""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
import sklearn


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def command_output(command: list[str]) -> str:
    completed = subprocess.run(command, check=True, text=True, capture_output=True)
    return completed.stdout.strip()


def local_hardware_summary() -> tuple[str, str]:
    """Return a sanitized machine summary without serial or device identifiers."""
    report = command_output(["system_profiler", "SPHardwareDataType"])

    def field(label: str, fallback: str = "not reported") -> str:
        match = re.search(rf"^\s*{re.escape(label)}:\s*(.+)$", report, flags=re.MULTILINE)
        return match.group(1).strip() if match else fallback

    chip = field("Chip")
    cores = field("Total Number of Cores")
    memory = field("Memory")
    model = field("Model Name")
    return f"{chip}; {cores}; {memory}", f"{model}; {platform.platform()}"


def resource_row(path: Path, stage: str, purpose: str, notes: str = "") -> dict[str, object]:
    record = json.loads(path.read_text())
    return {
        "stage": stage,
        "purpose": purpose,
        "wall_seconds": record["wall_seconds"],
        "user_cpu_seconds": record["user_cpu_seconds"],
        "system_cpu_seconds": record["system_cpu_seconds"],
        "peak_rss_bytes": record["peak_process_tree_rss_bytes"],
        "peak_rss_gib": record["peak_process_tree_rss_gib"],
        "memory_source": "50-ms psutil process RSS polling",
        "exit_code": record["exit_code"],
        "notes": notes,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--ntriage-bin", type=Path, required=True)
    parser.add_argument("--refs", type=Path, required=True)
    parser.add_argument("--ntriage-repo", type=Path, required=True)
    args = parser.parse_args()

    root = args.root.resolve()
    tables = root / "tables"
    logs = root / "logs"
    inputs = root / "inputs"
    native_prefix = root / "native" / "telomeric_2794_core"
    tables.mkdir(parents=True, exist_ok=True)

    resources = [
        resource_row(logs / "prepare_inputs_resources.json", "01_prepare_inputs",
                     "Validate 2,794 source FASTAs; write combined FASTA and core BED"),
        resource_row(logs / "samtools_faidx_core_resources.json", "02_fasta_index",
                     "Index the combined FASTA with samtools faidx"),
    ]

    minimap_log = (logs / "minimap2_acro1_cap3m.time_and_stderr.txt").read_text()
    match = re.search(r"Real time: ([0-9.]+) sec; CPU: ([0-9.]+) sec; Peak RSS: ([0-9.]+) GB", minimap_log)
    if not match:
        raise ValueError("Could not parse minimap2 resource line")
    minimap_wall, minimap_cpu, minimap_gb = map(float, match.groups())
    minimap_bytes = int(minimap_gb * 1_000_000_000)
    resources.append({
        "stage": "03_map_ACRO1_anchors", "purpose": "Sensitive minimap2 mapping of the 147-bp ACRO1 reference",
        "wall_seconds": minimap_wall, "user_cpu_seconds": minimap_cpu, "system_cpu_seconds": np.nan,
        "peak_rss_bytes": minimap_bytes, "peak_rss_gib": minimap_bytes / 1024 ** 3,
        "memory_source": "minimap2 built-in peak-RSS report", "exit_code": 0,
        "notes": "-x sr -k 9 -w 1 -f 0 -n 3 -m 20 -s 25 -p 0.50 -N 3000000; 8 threads",
    })
    resources.extend([
        resource_row(logs / "parse_acro1_core_resources.json", "04_filter_ACRO1_anchors",
                     "Filter and resolve nonredundant ACRO1 monomer anchors"),
        resource_row(logs / "ntriage_core_resources.json", "05_NTRIAGE_core",
                     "Core-only N-TRIAGE clustering and deep-subfamily analysis",
                     "fixed k=8; seed=1; 8 threads; 77,276 complete ACRO units"),
        resource_row(logs / "analysis_core_resources.json", "06_downstream_statistics",
                     "Individual-aware association tests, PCA, and sensitivity analyses",
                     "seed=20261026"),
    ])
    resource_table = pd.DataFrame(resources)
    resource_table.to_csv(tables / "runtime_memory_summary.tsv", sep="\t", index=False, na_rep="NA")
    total_wall = float(resource_table["wall_seconds"].sum())
    peak_gib = float(resource_table["peak_rss_gib"].max())
    pd.DataFrame([{
        "serial_stage_wall_seconds": total_wall,
        "serial_stage_wall_minutes": total_wall / 60,
        "maximum_observed_peak_rss_gib": peak_gib,
        "note": "Sum of measured stages; excludes figure rendering and the independent audit rerun.",
    }]).to_csv(tables / "runtime_memory_overall.tsv", sep="\t", index=False)

    manifest = pd.read_csv(args.manifest, sep="\t", na_values=["NA"])
    arrays = pd.read_csv(f"{native_prefix}.arrays.tsv", sep="\t")
    clusters = pd.read_csv(f"{native_prefix}.clusters.tsv", sep="\t")
    acro = pd.read_csv(f"{native_prefix}.acro1_composition.tsv", sep="\t")
    switches = pd.read_csv(f"{native_prefix}.switches.tsv", sep="\t")
    audit = pd.read_csv(inputs / "input_audit.tsv", sep="\t")
    cohort_rows = [
        ("arrays", len(manifest), "distinct manifest rows and core intervals"),
        ("individuals", manifest["individual_id"].nunique(), "independent people"),
        ("individual_haplotype_strata", manifest[["individual_id", "haplotype"]].drop_duplicates().shape[0],
         "haplotype 1/2 are phase labels, not parent-of-origin"),
        ("arrays_recombination", int(manifest["has_recombination"].sum()), "manifest label"),
        ("arrays_no_recombination", int((manifest["has_recombination"] == 0).sum()), "manifest label"),
        ("individuals_with_superpopulation", manifest.dropna(subset=["superpopulation"])["individual_id"].nunique(), ""),
        ("individuals_missing_superpopulation", manifest.loc[manifest["superpopulation"].isna(), "individual_id"].nunique(), ""),
        ("source_FASTA_bp_including_flanks", int(audit["sequence_length_bp"].sum()), "input extraction intervals"),
        ("manifest_defined_core_bp", int(audit["core_length_bp"].sum()), "sequence analyzed by N-TRIAGE"),
        ("complete_ACRO_units", len(clusters), "core-only N-TRIAGE segmentation"),
        ("mean_complete_units_per_array", len(clusters) / len(manifest), ""),
        ("median_complete_units_per_array", float(arrays["complete_acres"].median()), ""),
        ("ACRO1_monomers_classified", int(acro["classified_total"].sum()), "N-TRIAGE S1/S2/S3 assignments"),
        ("deep_subfamily_switch_events", len(switches), "nonredundant calls"),
    ]
    pd.DataFrame(cohort_rows, columns=["metric", "value", "definition"]).to_csv(
        tables / "cohort_run_summary.tsv", sep="\t", index=False
    )

    git_revision = command_output(["git", "-C", str(args.ntriage_repo), "rev-parse", "--short", "HEAD"])
    git_dirty = bool(command_output(["git", "-C", str(args.ntriage_repo), "status", "--porcelain", "--untracked-files=no"]))
    import psutil
    hardware_version, hardware_notes = local_hardware_summary()
    provenance = [
        ("benchmark hardware", hardware_version, "local machine", "NA", hardware_notes),
        ("N-TRIAGE", "0.2.0", str(args.ntriage_bin), sha256(args.ntriage_bin),
         f"git {git_revision}; working tree {'modified' if git_dirty else 'clean'}"),
        ("Dfam/component reference", "44 FASTA records", str(args.refs), sha256(args.refs), "dfam_2.fasta"),
        ("manifest", "2,794 rows", str(args.manifest.resolve()), sha256(args.manifest), "validated input manifest"),
        ("core BED", "2,794 intervals", str((inputs / "telomeric_2794.core.bed").resolve()),
         sha256(inputs / "telomeric_2794.core.bed"), "0-based half-open local core coordinates"),
        ("minimap2", command_output(["minimap2", "--version"]), command_output(["which", "minimap2"]), "NA", ""),
        ("samtools", command_output(["samtools", "--version"]).splitlines()[0].split()[-1],
         command_output(["which", "samtools"]), "NA", ""),
        ("Python", platform.python_version(), sys.executable, "NA", platform.platform()),
        ("pandas", pd.__version__, "Python package", "NA", ""),
        ("numpy", np.__version__, "Python package", "NA", ""),
        ("scipy", scipy.__version__, "Python package", "NA", ""),
        ("scikit-learn", sklearn.__version__, "Python package", "NA", ""),
        ("psutil", psutil.__version__, "Python package", "NA", "resource polling"),
    ]
    pd.DataFrame(provenance, columns=["software_or_resource", "version_or_size", "path_or_source",
                                      "sha256", "notes"]).to_csv(
        tables / "software_and_input_provenance.tsv", sep="\t", index=False
    )

    print(f"Wrote resource, cohort, and provenance tables under {tables}")


if __name__ == "__main__":
    main()
