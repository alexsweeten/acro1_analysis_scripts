#!/usr/bin/env python3
"""Run one command and record wall time, CPU time, and peak process-tree RSS."""

from __future__ import annotations

import argparse
import datetime as dt
import json
import resource
import shlex
import subprocess
import time
from pathlib import Path

import psutil


def process_tree_rss(pid: int) -> int:
    try:
        root = psutil.Process(pid)
        try:
            processes = [root] + root.children(recursive=True)
        except (psutil.AccessDenied, PermissionError):
            # macOS sandboxing may prohibit enumerating the global process table.
            # N-TRIAGE is single-process/multithreaded, so root RSS is sufficient.
            processes = [root]
    except (psutil.NoSuchProcess, psutil.AccessDenied, PermissionError):
        return 0
    total = 0
    for process in processes:
        try:
            total += process.memory_info().rss
        except (psutil.NoSuchProcess, psutil.AccessDenied, PermissionError):
            pass
    return total


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--label", required=True)
    parser.add_argument("--metrics", type=Path, required=True)
    parser.add_argument("--stdout", type=Path, required=True)
    parser.add_argument("--stderr", type=Path, required=True)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command
    if command and command[0] == "--":
        command = command[1:]
    if not command:
        parser.error("a command is required after --")

    for path in (args.metrics, args.stdout, args.stderr):
        path.parent.mkdir(parents=True, exist_ok=True)

    start_iso = dt.datetime.now(dt.timezone.utc).isoformat()
    start = time.perf_counter()
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    peak_rss = 0
    with args.stdout.open("w") as stdout_handle, args.stderr.open("w") as stderr_handle:
        process = subprocess.Popen(command, stdout=stdout_handle, stderr=stderr_handle)
        while process.poll() is None:
            peak_rss = max(peak_rss, process_tree_rss(process.pid))
            time.sleep(0.05)
        peak_rss = max(peak_rss, process_tree_rss(process.pid))
        return_code = process.wait()
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    wall_seconds = time.perf_counter() - start
    end_iso = dt.datetime.now(dt.timezone.utc).isoformat()

    # ru_maxrss is bytes on macOS and KiB on Linux. psutil provides the primary,
    # platform-independent process-tree peak used in the report.
    record = {
        "label": args.label,
        "command": shlex.join(command),
        "start_utc": start_iso,
        "end_utc": end_iso,
        "exit_code": return_code,
        "wall_seconds": wall_seconds,
        "user_cpu_seconds": after.ru_utime - before.ru_utime,
        "system_cpu_seconds": after.ru_stime - before.ru_stime,
        "peak_process_tree_rss_bytes": peak_rss,
        "peak_process_tree_rss_gib": peak_rss / (1024 ** 3),
        "rusage_maxrss_platform_units": after.ru_maxrss,
        "poll_interval_seconds": 0.05,
    }
    args.metrics.write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record, indent=2))
    raise SystemExit(return_code)


if __name__ == "__main__":
    main()
