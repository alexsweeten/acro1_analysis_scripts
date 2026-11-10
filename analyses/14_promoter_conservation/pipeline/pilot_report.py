#!/usr/bin/env python3
"""Summarise pilot runtime / memory (own timers + sacct) and propose SLURM limits.

  pilot_report.py JOBID   -> prints tables, writes pilot_report.tsv and pilot_resources.tsv
Limits are chosen from the measured maxima with safety factors; E-INS-i time is extrapolated with a log-log fit on
the number of complete units, which is the only super-linear stage.
"""
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pp_common as C

job = sys.argv[1]
sel = pd.read_csv(C.ROOT / "pilot_selection.tsv", sep="\t")
man = pd.read_csv(C.ROOT / "manifest.validated.tsv", sep="\t")
rows = []
for _, s in sel.iterrows():
    sid = man.loc[man.task_index == s.task_index, "safe_id"].iloc[0]
    f = C.PER_ARRAY / sid / "result.tsv"
    if f.exists():
        r = pd.read_csv(f, sep="\t").iloc[0]
        rows.append(dict(task_index=s.task_index, stratum=s.stratum, fasta_length=r.fasta_length, n_anchors=r.n_anchors,
                         n_units=r.n_units_complete, status=r.primary_status, reason=r.primary_reason if isinstance(r.primary_reason, str) else "",
                         bb_status=r.bb_status, t_total_s=r.runtime_total_s, t_nhmmer_s=r.t_nhmmer_s, t_mafft_add_s=r.t_mafft_add_s,
                         t_einsi_s=r.t_mafft_einsi_s, t_backbone_s=r.t_mafft_backbone_s, peak_rss_mb=r.peak_rss_mb))
d = pd.DataFrame(rows).sort_values("fasta_length")
sa = subprocess.run(["sacct", "-j", job, "--format=JobID%20,State%12,Elapsed,MaxRSS,TotalCPU", "-n", "-P"], capture_output=True, text=True).stdout
acct = []
for l in sa.strip().splitlines():
    j, st, el, rss, cpu = l.split("|")
    if j.endswith(".batch"):
        acct.append(dict(array_task=int(j.split("_")[1].split(".")[0]), sacct_state=st, sacct_elapsed=el, sacct_maxrss_mb=float(rss.rstrip("K") or 0) / 1024))
acct = pd.DataFrame(acct)
pd.set_option("display.width", 220, "display.max_columns", 30)
print(d.round(1).to_string(index=False))
d.to_csv(C.ROOT / "pilot_report.tsv", sep="\t", index=False)
print("\nsacct MaxRSS (MB): max=%.0f median=%.0f" % (acct.sacct_maxrss_mb.max(), acct.sacct_maxrss_mb.median()))

ok = d[d.t_einsi_s > 1]
b, a = np.polyfit(np.log(ok.n_units), np.log(ok.t_einsi_s), 1)
print(f"\nE-INS-i time ~ {np.exp(a):.3g} * units^{b:.2f}   (n={len(ok)} arrays with >1 unit)")
mx = man.fasta_length.max()
per_kb = (d.n_units / d.fasta_length).median()
print("manifest length quantiles (kb):", {q: int(man.fasta_length.quantile(q) / 1000) for q in (.5, .9, .99, 1.0)})
pred = {q: np.exp(a) * (per_kb * man.fasta_length.quantile(q)) ** b for q in (.5, .9, .99, 1.0)}
print("predicted E-INS-i seconds at length quantiles:", {k: round(v) for k, v in pred.items()})
t_all = d.t_total_s.max()
rec = dict(max_pilot_total_s=t_all, max_pilot_peak_rss_mb=max(d.peak_rss_mb.max(), acct.sacct_maxrss_mb.max()))
pd.Series(rec).to_csv(C.ROOT / "pilot_resources.tsv", sep="\t", header=False)
print(rec)
