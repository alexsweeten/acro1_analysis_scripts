#!/usr/bin/env python3
"""Independent audit: recompute the primary and backbone results of finished arrays with naive loops
(no shared code with pp_core.py) and compare them to result.tsv.

  audit_array.py SAFE_ID [SAFE_ID ...]
"""
import itertools
import statistics
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pp_common as C   # constants + read_fasta only

W = 50


def col_pairs(seqs, c):
    col = [s[c] for s in seqs if s[c] in "ACGT"]
    n = len(col)
    d = sum(a != b for a, b in itertools.combinations(col, 2))
    return d, n * (n - 1) // 2, col


def window(seqs, s):          # s: 0-based start
    D = P = call = gc = 0
    for c in range(s, s + W):
        d, p, col = col_pairs(seqs, c)
        D, P, call, gc = D + d, P + p, call + len(col), gc + sum(x in "GC" for x in col)
    return D / P if P else float("nan"), gc / call if call else float("nan"), call / (len(seqs) * W)


def eligible(seqs, s):
    return all(sum(q[c] in "ACGT" for q in seqs) >= 0.8 * len(seqs) - 1e-9 for c in range(s, s + W))


def cpg_gapless(seqs, s):
    k = n = 0
    for q in seqs:
        f = [x for x in q[s:s + W] if x != "-"]
        for a, b in zip(f, f[1:]):
            if a in "ACGT" and b in "ACGT":
                n += 1
                k += a == "C" and b == "G"
    return k / n if n else float("nan")


def cpg_adj(seqs, s):
    k = n = 0
    for q in seqs:
        for c in range(s, s + W - 1):
            if q[c] in "ACGT" and q[c + 1] in "ACGT":
                n += 1
                k += q[c] == "C" and q[c + 1] == "G"
    return k / n if n else float("nan")


def stat(seqs, s, cpgfn):
    d, gc, occ = window(seqs, s)
    return d, gc, cpgfn(seqs, s), occ


def test(seqs, prom0, ctrl_ok, cpgfn, tol_cpg):
    pd_, pgc, pcpg, pocc = stat(seqs, prom0, cpgfn)
    divs = []
    for s in range(0, len(seqs[0]) - W + 1):
        if s == prom0 or not ctrl_ok(s):
            continue
        d, gc, cp, occ = stat(seqs, s, cpgfn)
        if abs(gc - pgc) <= C.TOL_GC + 1e-12 and abs(cp - pcpg) <= tol_cpg + 1e-12 and abs(occ - pocc) <= C.TOL_OCC + 1e-12:
            divs.append(d)
    n = len(divs)
    if not n:
        return dict(pd=pd_, n=0)
    return dict(pd=pd_, n=n, med=statistics.median(divs), effect=statistics.median(divs) - pd_,
                p=(1 + sum(x <= pd_ + 1e-12 for x in divs)) / (1 + n))


def derived_mask():
    rm = pd.read_csv(C.RM_TABLE, sep="\t")
    full = [s.upper() for _, s in C.read_fasta(C.REF_ALN9090)]
    N = C.FULL_COLS
    ok = [True] * (N + 1)
    for c in range(1, N + 1):
        if sum(q[c - 1] != "-" for q in full) / len(full) < 0.8:
            ok[c] = False
    banned = [False] * (N + 1)
    for _, r in rm.iterrows():
        nm, cls = r["repeat"].split("#")
        a, b = r.alignment_start_1based, r.alignment_end_1based_inclusive
        if cls.startswith("Simple_repeat") or cls.startswith("Satellite"):
            for c in range(a, b + 1):
                banned[c] = True
        if nm == "ACRO1":
            for c in range(max(1, a - 25), min(N, b + 25) + 1):
                banned[c] = True
    return ok, banned


def main():
    ok_ref, banned = derived_mask()
    N = C.FULL_COLS
    allbad = 0
    for sid in sys.argv[1:]:
        d = C.PER_ARRAY / sid
        r = pd.read_csv(d / "result.tsv", sep="\t").iloc[0]
        seqs = [s.upper() for _, s in C.read_fasta(d / "l1pa10_anchors.aligned.fa")]
        assert all(len(s) == 300 for s in seqs)
        a = test(seqs, 77, lambda s: eligible(seqs, s) and (s + W - 1 < 77 or s > 127 - 1) and not (s + W - 1 >= 77 and s <= 126), cpg_gapless, C.TOL_CPG)
        lines = [f"{r.array_id}: n_anchors={len(seqs)}"]
        chk = []
        chk.append(("primary promoter_div", a["pd"], r.promoter_divergence))
        chk.append(("primary n_controls", a["n"], r.n_matched_controls))
        if a["n"]:
            chk += [("primary control_median", a["med"], r.control_median_divergence), ("primary effect", a["effect"], r.effect), ("primary p", a["p"], r.p_matched)]
        if r.bb_status == "estimable":
            U = [s.upper() for _, s in C.read_fasta(d / "backbone" / "units_reference_mapped.fa")]
            assert all(len(s) == N for s in U)

            def prim_ok(s):
                cols = range(s + 1, s + W + 1)
                return eligible(U, s) and all(ok_ref[c] and not banned[c] and c > 296 for c in cols) and (s + W < 74 or s + 1 > 123)

            b = test(U, 73, prim_ok, cpg_adj, C.BACKBONE_TOL_CPG)
            chk += [("backbone promoter_div", b["pd"], r.bb_promoter_divergence), ("backbone n_controls", b["n"], r.bb_n_matched_controls)]
            if b["n"]:
                chk += [("backbone control_median", b["med"], r.bb_control_median_divergence), ("backbone effect", b["effect"], r.bb_effect), ("backbone p", b["p"], r.bb_p_matched)]
        bad = 0
        for name, mine, theirs in chk:
            same = abs(mine - theirs) < 1e-9 if isinstance(theirs, float) or isinstance(mine, float) else mine == theirs
            if not same:
                bad += 1
                lines.append(f"   MISMATCH {name}: audit={mine} pipeline={theirs}")
        lines.append(f"   {len(chk) - bad}/{len(chk)} quantities agree")
        allbad += bad
        print("\n".join(lines))
    sys.exit(1 if allbad else 0)


if __name__ == "__main__":
    main()
