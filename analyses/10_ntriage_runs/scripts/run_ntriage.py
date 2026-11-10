#!/usr/bin/env python3
"""Run NTRIAGE reproducibly for each prepared cohort and the joint input.

The native two-cluster assignment is deliberately provisional: downstream
``select_k.py`` chooses K with replicated clustering on the native feature
matrix. Segmentation, annotation, ACRO1 allele calls, and deep-subfamily output
are nevertheless produced once here and reused downstream.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import shlex
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Mapping, Sequence


DEFAULT_ROOT = Path(os.environ.get("ACRO1_NTRIAGE_ANALYSIS_ROOT", Path(__file__).resolve().parents[1]))
DEFAULT_BINARY = Path(os.environ.get("NTRIAGE_BINARY", "ntriage"))
DEFAULT_REFS = Path(os.environ.get("ACRO1_COMPONENT_REFERENCES", "component_references.fasta"))
REQUIRED_SUFFIXES = (
    ".arrays.tsv",
    ".acres.bed",
    ".acres.fasta",
    ".annotations.bed",
    ".clusters.tsv",
    ".features.tsv",
    ".acro1_subfamilies.bed",
    ".acro1_composition.tsv",
)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"dataset_id", "fasta", "bed", "acro1_monomers", "ntriage_prefix"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"{path}: missing columns {sorted(required - set(rows[0] if rows else []))}")
    if len({row["dataset_id"] for row in rows}) != len(rows):
        raise ValueError(f"{path}: dataset_id values are not unique")
    return rows


def config_hash(config: Mapping[str, object]) -> str:
    encoded = json.dumps(config, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def nonempty_outputs(prefix: Path) -> tuple[list[str], list[str]]:
    missing: list[str] = []
    empty: list[str] = []
    for suffix in REQUIRED_SUFFIXES:
        path = Path(str(prefix) + suffix)
        if not path.exists():
            missing.append(str(path))
        elif path.stat().st_size == 0:
            empty.append(str(path))
    return missing, empty


def build_command(
    row: Mapping[str, str],
    binary: Path,
    refs: Path,
    prefix: Path,
    args: argparse.Namespace,
) -> list[str]:
    return [
        str(binary),
        "--bed", row["bed"],
        "--fasta", row["fasta"],
        "--refs", str(refs),
        "--out", str(prefix),
        "--acro1-monomers", row["acro1_monomers"],
        "--kmer", str(args.kmer),
        "--bins", str(args.bins),
        "--threads", str(args.threads),
        "--clusters", str(args.provisional_clusters),
        "--min-period", str(args.min_period),
        "--max-period", str(args.max_period),
        "--period-target", str(args.period_target),
        "--seed", str(args.seed),
        "--deep-subfamilies",
        "--write-fasta",
    ]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the prepared comparative NTRIAGE inputs sequentially and resumably."
    )
    parser.add_argument("--analysis-root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument(
        "--summary-output",
        type=Path,
        default=Path("provenance/run_ntriage_summary.json"),
        help=(
            "Run-summary JSON, relative to --analysis-root unless absolute "
            "(default: provenance/run_ntriage_summary.json)"
        ),
    )
    parser.add_argument("--binary", type=Path, default=DEFAULT_BINARY)
    parser.add_argument("--refs", type=Path, default=DEFAULT_REFS)
    parser.add_argument(
        "--dataset", nargs="+", help="dataset IDs to run; default: every manifest row"
    )
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--kmer", type=int, default=15)
    parser.add_argument("--bins", type=int, default=256)
    parser.add_argument("--provisional-clusters", type=int, default=2)
    parser.add_argument("--min-period", type=int, default=6000)
    parser.add_argument("--max-period", type=int, default=9000)
    parser.add_argument("--period-target", type=int, default=6650)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--force", action="store_true", help="rerun completed datasets")
    parser.add_argument("--dry-run", action="store_true", help="print commands without running")
    parser.add_argument(
        "--continue-on-error", action="store_true", help="attempt later datasets after a failure"
    )
    parser.add_argument(
        "--skip-self-test", action="store_true", help="do not run the NTRIAGE built-in self-test"
    )
    parser.add_argument("--list", action="store_true", help="list manifest datasets and exit")
    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    if args.threads < 1:
        raise ValueError("--threads must be positive")
    if args.kmer < 3 or args.bins < 1 or args.provisional_clusters < 2:
        raise ValueError("invalid k-mer, feature-bin, or provisional-cluster setting")
    if not (1 <= args.min_period <= args.period_target <= args.max_period):
        raise ValueError("require min-period <= period-target <= max-period")


def run_one(
    row: Mapping[str, str], args: argparse.Namespace, binary_sha: str, refs_sha: str
) -> str:
    dataset = row["dataset_id"]
    prefix = Path(row["ntriage_prefix"])
    run_dir = prefix.parent.parent
    native_dir = prefix.parent
    logs_dir = run_dir / "logs"
    run_dir.mkdir(parents=True, exist_ok=True)
    native_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)
    command = build_command(row, args.binary.resolve(), args.refs.resolve(), prefix, args)
    run_config = {
        "dataset_id": dataset,
        "command": command,
        "binary_sha256": binary_sha,
        "refs_sha256": refs_sha,
        "input_sha256": {
            "fasta": row.get("fasta_sha256", ""),
            "bed": row.get("bed_sha256", ""),
            "acro1_monomers": row.get("acro1_monomers_sha256", ""),
        },
        "segmentation": {
            "min_period_bp": args.min_period,
            "max_period_bp": args.max_period,
            "period_target_bp": args.period_target,
        },
        "clustering_note": (
            "Native K=2 labels are provisional; authoritative K and labels are selected "
            "downstream from the NTRIAGE feature matrix."
        ),
    }
    fingerprint = config_hash(run_config)
    complete_path = run_dir / "complete.json"
    if complete_path.exists() and not args.force:
        completed = json.loads(complete_path.read_text())
        missing, empty = nonempty_outputs(prefix)
        if completed.get("config_sha256") == fingerprint and not missing and not empty:
            print(f"SKIP {dataset}: completed with matching configuration")
            return "skipped"
        raise RuntimeError(
            f"{dataset}: completion marker/configuration mismatch or missing outputs; "
            "inspect the run and use --force to rerun"
        )

    print(f"RUN  {dataset}: {shlex.join(command)}")
    if args.dry_run:
        return "dry-run"
    # A forced rerun must not retain an earlier success marker: if the new
    # process fails, a later invocation should never mistake stale outputs for
    # the requested configuration.
    if args.force and complete_path.exists():
        complete_path.unlink()
    started = utc_now()
    command_record = {
        **run_config,
        "config_sha256": fingerprint,
        "started_utc": started,
        "host": platform.node(),
        "platform": platform.platform(),
        "python": sys.version,
    }
    (run_dir / "command.json").write_text(json.dumps(command_record, indent=2) + "\n")
    stdout_path = logs_dir / "ntriage.stdout.log"
    stderr_path = logs_dir / "ntriage.stderr.log"
    environment = os.environ.copy()
    environment.setdefault("OMP_NUM_THREADS", str(args.threads))
    with stdout_path.open("w") as stdout_handle, stderr_path.open("w") as stderr_handle:
        completed = subprocess.run(
            command, stdout=stdout_handle, stderr=stderr_handle, env=environment, check=False
        )
    if completed.returncode != 0:
        failure = {
            "dataset_id": dataset,
            "returncode": completed.returncode,
            "started_utc": started,
            "failed_utc": utc_now(),
            "config_sha256": fingerprint,
            "stderr_log": str(stderr_path),
        }
        (run_dir / "failed.json").write_text(json.dumps(failure, indent=2) + "\n")
        raise RuntimeError(f"{dataset}: NTRIAGE failed ({completed.returncode}); see {stderr_path}")
    missing, empty = nonempty_outputs(prefix)
    if missing or empty:
        raise RuntimeError(f"{dataset}: output validation failed; missing={missing}; empty={empty}")
    completed_record = {
        "dataset_id": dataset,
        "status": "complete",
        "config_sha256": fingerprint,
        "started_utc": started,
        "completed_utc": utc_now(),
        "returncode": 0,
        "output_prefix": str(prefix),
        "validated_outputs": [str(prefix) + suffix for suffix in REQUIRED_SUFFIXES],
    }
    complete_path.write_text(json.dumps(completed_record, indent=2) + "\n")
    failed_path = run_dir / "failed.json"
    if failed_path.exists():
        failed_path.unlink()
    print(f"DONE {dataset}: {prefix}")
    return "complete"


def main() -> None:
    args = parse_args()
    validate_args(args)
    root = args.analysis_root.resolve()
    manifest = args.manifest or (root / "inputs/dataset_manifest.tsv")
    rows = read_manifest(manifest)
    available = [row["dataset_id"] for row in rows]
    if args.list:
        print("\n".join(available))
        return
    selected = args.dataset or available
    unknown = set(selected) - set(available)
    if unknown:
        raise ValueError(f"unknown dataset IDs: {sorted(unknown)}; available={available}")
    rows_by_id = {row["dataset_id"]: row for row in rows}
    for path in (args.binary, args.refs):
        if not path.is_file():
            raise FileNotFoundError(path)
    for dataset in selected:
        row = rows_by_id[dataset]
        for column in ("fasta", "bed", "acro1_monomers"):
            if not Path(row[column]).is_file():
                raise FileNotFoundError(f"{dataset} {column}: {row[column]}")

    binary_sha = sha256_file(args.binary)
    refs_sha = sha256_file(args.refs)
    if not args.skip_self_test and not args.dry_run:
        print("Running NTRIAGE built-in self-test")
        subprocess.run([str(args.binary), "--self-test"], check=True)
    statuses: dict[str, str] = {}
    failures: dict[str, str] = {}
    for dataset in selected:
        try:
            statuses[dataset] = run_one(rows_by_id[dataset], args, binary_sha, refs_sha)
        except Exception as error:
            failures[dataset] = str(error)
            print(f"ERROR {dataset}: {error}", file=sys.stderr)
            if not args.continue_on_error:
                raise
    summary = {
        "created_utc": utc_now(),
        "manifest": str(manifest.resolve()),
        "statuses": statuses,
        "failures": failures,
    }
    if not args.dry_run:
        summary_output = args.summary_output.expanduser()
        if not summary_output.is_absolute():
            summary_output = (root / summary_output).resolve()
        summary_output.parent.mkdir(parents=True, exist_ok=True)
        summary_output.write_text(json.dumps(summary, indent=2) + "\n")
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
