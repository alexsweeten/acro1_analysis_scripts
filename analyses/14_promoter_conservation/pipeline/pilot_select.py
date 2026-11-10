#!/usr/bin/env python3
"""Deterministically choose a pilot spanning short / median / long input lengths and both recombination classes.

  pilot_select.py [--n 27] > pilot_tasks.txt
Arrays are ranked by FASTA length and split into three strata; within each stratum ranks are evenly spaced.
The longest array is always included (worst case for E-INS-i). Nothing about divergence is used.
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pp_common as C

ap = argparse.ArgumentParser()
ap.add_argument("--n", type=int, default=27)
a = ap.parse_args()
m = pd.read_csv(C.ROOT / "manifest.validated.tsv", sep="\t").sort_values(["fasta_length", "task_index"]).reset_index(drop=True)
per = a.n // 3
picks = []
edges = np.linspace(0, len(m), 4).astype(int)
for lo, hi in zip(edges[:-1], edges[1:]):
    picks += list(np.linspace(lo, hi - 1, per + 1)[:-1].round().astype(int) + 0)
picks[-1] = len(m) - 1                       # longest array
picks = list(dict.fromkeys(picks))
sel = m.loc[picks, ["task_index", "array_id", "fasta_length", "recombination_class"]]
sel["stratum"] = ["short"] * per + ["median"] * per + ["long"] * (len(sel) - 2 * per)
sel.to_csv(C.ROOT / "pilot_selection.tsv", sep="\t", index=False)
print("\n".join(str(i) for i in sel.task_index.tolist()))
print(sel.groupby(["stratum", "recombination_class"]).size().to_string(), file=sys.stderr)
