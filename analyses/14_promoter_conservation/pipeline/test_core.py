#!/usr/bin/env python3
"""Unit tests: vectorised statistics vs naive brute force, chain DP vs exhaustive search, strand geometry."""
import itertools
import random
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pp_common as C
import pp_core as K

rng = random.Random(1)


def rand_aln(n, L, p_gap=0.1, p_amb=0.03):
    out = []
    for _ in range(n):
        s = []
        for _ in range(L):
            r = rng.random()
            s.append("-" if r < p_gap else "N" if r < p_gap + p_amb else rng.choice("ACGT"))
        out.append("".join(s))
    return out


def brute_window(seqs, s, w):
    diffs = pairs = 0
    gc = call = 0
    for c in range(s, s + w):
        col = [q[c] for q in seqs if q[c] in "ACGT"]
        call += len(col)
        gc += sum(b in "GC" for b in col)
        for a, b in itertools.combinations(col, 2):
            pairs += 1
            diffs += a != b
    cpg = dn = 0
    for q in seqs:
        f = [b for b in q[s:s + w] if b != "-"]
        for a, b in zip(f, f[1:]):
            if a in "ACGT" and b in "ACGT":
                dn += 1
                cpg += (a == "C" and b == "G")
    adj_cg = adj_n = 0
    for q in seqs:
        for c in range(s, s + w - 1):
            if q[c] in "ACGT" and q[c + 1] in "ACGT":
                adj_n += 1
                adj_cg += (q[c] == "C" and q[c + 1] == "G")
    return dict(div=diffs / pairs if pairs else np.nan, gc=gc / call if call else np.nan,
                occ=call / (len(seqs) * w), cpg=cpg / dn if dn else np.nan,
                cpg_adj=adj_cg / adj_n if adj_n else np.nan)


def test_window_stats():
    for n in (2, 5, 9):
        seqs = rand_aln(n, 120)
        X = K.encode(seqs)
        cs = K.col_stats(X)
        wm = K.window_metrics(cs)
        cg = K.cpg_gapless(X)
        ca = K.cpg_adjacent(X)
        for s in (0, 7, 33, 70):
            b = brute_window(seqs, s, C.WIN)
            assert np.isclose(wm["divergence"][s], b["div"], equal_nan=True), (n, s)
            assert np.isclose(wm["gc"][s], b["gc"], equal_nan=True)
            assert np.isclose(wm["occupancy"][s], b["occ"])
            assert np.isclose(cg[s], b["cpg"], equal_nan=True), (n, s, cg[s], b["cpg"])
            assert np.isclose(ca[s], b["cpg_adj"], equal_nan=True), (n, s)
    print("window statistics == brute force")


def test_not_mean_of_columns():
    # window with one gap-heavy column: pooled != mean of per-column divergence
    seqs = ["A" * 49 + "A", "A" * 49 + "C", "A" * 49 + "-", "A" * 49 + "-"]
    X = K.encode(seqs)
    cs = K.col_stats(X)
    pooled = K.window_metrics(cs)["divergence"][0]
    per_col = np.where(cs["pairs"] > 0, cs["diffs"] / np.maximum(cs["pairs"], 1), 0)
    assert abs(pooled - 1 / (49 * 6 + 1)) < 1e-12 and abs(pooled - per_col.mean()) > 1e-4
    print("pooled-not-averaged divergence OK")


def test_eligibility():
    seqs = ["A" * 100] * 5
    seqs = [list(s) for s in seqs]
    for i in range(2):          # 2/5 = 40% gaps at column 10 -> 60% callable < 80%
        seqs[i][10] = "-"
    seqs = ["".join(s) for s in seqs]
    wm = K.window_metrics(K.col_stats(K.encode(seqs)))
    assert not wm["eligible"][0] and not wm["eligible"][10] and wm["eligible"][11]
    seqs2 = ["A" * 100] * 5
    seqs2 = [list(s) for s in seqs2]
    seqs2[0][10] = "-"          # 4/5 = 80% exactly -> eligible
    wm = K.window_metrics(K.col_stats(K.encode(["".join(s) for s in seqs2])))
    assert wm["eligible"][0]
    print("80% eligibility boundary OK")


def test_p_and_effect():
    s = K.summarise_vs_controls(0.05, [0.04, 0.05, 0.06, 0.07])
    assert s["p"] == (1 + 2) / (1 + 4) and abs(s["effect"] - (0.055 - 0.05)) < 1e-12 and s["p_min"] == 1 / 5
    assert np.isnan(K.summarise_vs_controls(0.05, [])["effect"])
    print("empirical p / effect OK")


def test_chain_dp():
    for trial in range(300):
        n = rng.randint(1, 9)
        starts = sorted(rng.sample(range(0, 60000, 100), n))
        A = [dict(start=s, end=s + 300, score=rng.uniform(100, 400), strand="+", hit_id=i) for i, s in enumerate(starts)]
        got = K.best_chain(A)
        best = (0, 0.0)
        for r in range(1, n + 1):
            for comb in itertools.combinations(range(n), r):
                ok = all(C.SPACING_MIN <= A[comb[k + 1]]["start"] - A[comb[k]]["start"] <= C.SPACING_MAX for k in range(len(comb) - 1))
                if ok:
                    cand = (r, sum(A[i]["score"] for i in comb))
                    if cand > best:
                        best = cand
        assert len(got) == best[0] and abs(sum(a["score"] for a in got) - best[1]) < 1e-6, (trial, len(got), best)
    print("chain DP == exhaustive search (300 random cases)")


def test_dedup():
    hs = [dict(hit_id=0, start=1, end=300, score=200.0, evalue=1e-30, filter_status="pass"),
          dict(hit_id=1, start=151, end=450, score=150.0, evalue=1e-20, filter_status="pass"),   # 50% of shorter -> dup
          dict(hit_id=2, start=152, end=450, score=140.0, evalue=1e-20, filter_status="pass"),   # <50% -> kept
          dict(hit_id=3, start=1000, end=1300, score=90.0, evalue=1e-20, filter_status="fail_score")]
    kept = K.dedup_hits(hs)
    assert [h["hit_id"] for h in kept] == [0, 2] and hs[1]["dedup_status"] == "rejected_duplicate"
    print("dedup (>=50% of shorter, higher score kept) OK")


def test_units_geometry():
    rs = random.Random(3)
    cons = "".join(rs.choice("ACGT") for _ in range(400))
    spacer = lambda n: "".join(rs.choice("ACGT") for _ in range(n))
    # + strand: anchors at 1001, 6001, 12001 ; each anchor = cons
    seq = spacer(1000) + cons + spacer(5600) + cons + spacer(5600) + cons + spacer(500)
    st = [seq.find(cons)]
    for _ in range(2):
        st.append(seq.find(cons, st[-1] + 1))
    chain = [dict(start=s + 1, end=s + 400, score=300.0, strand="+") for s in st]
    units = K.build_units(seq, chain, "+")
    assert len(units) == 2 and all(u["seq"].startswith(cons) for u in units)
    assert units[0]["local_start"] == chain[0]["start"] and units[0]["local_end"] == chain[1]["start"] - 1
    # - strand: genome is the reverse complement of the above
    rc = C.revcomp(seq)
    st2 = [rc.find(C.revcomp(cons))]
    for _ in range(2):
        st2.append(rc.find(C.revcomp(cons), st2[-1] + 1))
    chain2 = [dict(start=s + 1, end=s + 400, score=300.0, strand="-") for s in st2]
    units2 = K.build_units(rc, chain2, "-")
    assert len(units2) == 2 and all(u["seq"].startswith(cons) for u in units2)
    assert [u["seq"] for u in units2] == [u["seq"] for u in units], "strand-symmetric units must be identical"
    for u in units2:
        assert K.oriented_anchor(rc, u["start_anchor"]) == cons
    print("+/- strand unit extraction is symmetric and begins at the anchor")


if __name__ == "__main__":
    for t in (test_window_stats, test_not_mean_of_columns, test_eligibility, test_p_and_effect,
              test_chain_dp, test_dedup, test_units_geometry):
        t()
    print("ALL TESTS PASSED")
