#!/usr/bin/env python3
"""List manifest task indices that have no valid COMPLETE sentinel (stdout), with a summary on stderr.

  find_incomplete.py            -> incomplete or never-run tasks
  find_incomplete.py --subset FILE   restrict to task indices in FILE (e.g. the pilot)
Completed arrays are never listed, so resubmitting this list cannot overwrite finished work.
"""
import argparse
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pp_common as C

ap = argparse.ArgumentParser()
ap.add_argument("--subset")
ap.add_argument("--min-len", type=int, default=0, help="keep arrays with FASTA length > this (resource class split)")
ap.add_argument("--max-len", type=int, default=10**12, help="keep arrays with FASTA length <= this")
ap.add_argument("--clean-stale-tmp", action="store_true", help="remove per_array/.tmp_* left by killed tasks (refuses while pp_ jobs are queued/running)")
a = ap.parse_args()
if a.clean_stale_tmp:
    import shutil
    import subprocess
    live = subprocess.run(["squeue", "-u", __import__("os").environ["USER"], "-h", "-o", "%j"], capture_output=True, text=True).stdout.split()
    if any(j.startswith("pp_") for j in live):
        sys.exit("pp_ jobs are still queued/running; not removing .tmp_* directories")
    stale = list(C.PER_ARRAY.glob(".tmp_*"))
    for d in stale:
        shutil.rmtree(d, ignore_errors=True)
    print(f"removed {len(stale)} stale tmp dirs", file=sys.stderr)
m = pd.read_csv(C.ROOT / "manifest.validated.tsv", sep="\t", usecols=["task_index", "safe_id", "array_id", "fasta_length"])
m = m[(m.fasta_length > a.min_len) & (m.fasta_length <= a.max_len)]
if a.subset:
    keep = {int(x) for x in open(a.subset).read().split()}
    m = m[m.task_index.isin(keep)]
todo, failed, never = [], [], 0
for t, sid in zip(m.task_index, m.safe_id):
    if (C.PER_ARRAY / sid / "COMPLETE").exists():
        continue
    todo.append(int(t))
    f = C.ROOT / "failed" / sid / "FAILED.json"
    if f.exists():
        failed.append((int(t), json.loads(f.read_text())["reason_code"]))
    else:
        never += 1
print("\n".join(map(str, todo)))
print(f"incomplete={len(todo)} (previous computational failures={len(failed)}, never_run_or_killed={never})", file=sys.stderr)
for t, r in failed[:20]:
    print(f"  failed task {t}: {r}", file=sys.stderr)
