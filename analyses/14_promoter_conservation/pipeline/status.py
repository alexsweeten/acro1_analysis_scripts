#!/usr/bin/env python3
"""Progress / failure monitor for the full run.

  status.py [--jobs JOBID ...]   (defaults to every job id in logs/submitted_jobs.txt)
"""
import argparse
import collections
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pp_common as C

ap = argparse.ArgumentParser()
ap.add_argument("--jobs", nargs="*")
a = ap.parse_args()
man = pd.read_csv(C.ROOT / "manifest.validated.tsv", sep="\t", usecols=["task_index", "safe_id"])
done = sum((C.PER_ARRAY / s / "COMPLETE").exists() for s in man.safe_id)
fails = list((C.ROOT / "failed").glob("*/FAILED.json")) if (C.ROOT / "failed").exists() else []
tmp = list(C.PER_ARRAY.glob(".tmp_*"))
print(f"COMPLETE {done}/{len(man)}   computational-failure records {len(fails)}   in-progress tmp dirs {len(tmp)}")
if fails:
    for k, v in collections.Counter(json.loads(f.read_text())["reason_code"] for f in fails).most_common():
        print(f"   failure reason {k}: {v}")
jobs = a.jobs or ([l.strip() for l in open(C.ROOT / "logs/submitted_jobs.txt")] if (C.ROOT / "logs/submitted_jobs.txt").exists() else [])
if jobs:
    out = subprocess.run(["sacct", "-j", ",".join(jobs), "-X", "-n", "-P", "--format=JobID,State,Elapsed,ReqMem,MaxRSS,ExitCode"],
                         capture_output=True, text=True).stdout.strip().splitlines()
    st = collections.Counter(l.split("|")[1].split()[0] for l in out if l)
    print("Slurm task states:", dict(st))
    bad = [l for l in out if l and l.split("|")[1].split()[0] in ("FAILED", "TIMEOUT", "OUT_OF_MEMORY", "NODE_FAIL", "CANCELLED")]
    for l in bad[:25]:
        print("   ", l)
    if len(bad) > 25:
        print(f"    ... {len(bad) - 25} more")
