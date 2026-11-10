#!/usr/bin/env python3
import itertools
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pp_common as C
import pp_stats as S


def test_bh():
    q = S.bh([0.01, 0.04, 0.03, 0.005, np.nan])
    # sorted p: .005,.01,.03,.04 -> *4/rank: .02,.02,.04,.04 ; monotone
    assert np.allclose(q[:4], [0.02, 0.04, 0.04, 0.02]) and np.isnan(q[4]), q
    print("BH OK")


def test_signflip_exact_and_mc():
    x = np.array([0.3, 0.1, 0.2, -0.05, 0.15, 0.25, 0.4, -0.1])
    ex = S.signflip_mean_gt0(x, 10, "t")
    assert ex["method"] == "exact" and ex["draws"] == 256
    brute = np.mean([(np.array(s) @ x) / 8 >= x.mean() - 1e-12 for s in itertools.product((-1, 1), repeat=8)])
    assert abs(ex["p"] - brute) < 1e-12
    big = np.concatenate([x, x, x])        # n=24 -> Monte Carlo
    mc = S.signflip_mean_gt0(big, 200_000, "t2")
    assert mc["method"] == "monte_carlo" and mc["draws"] == 200_000
    # all-positive effects: observed is the max -> p = (1+1)/(D+1) ... only the all-plus pattern ties
    pos = np.full(30, 0.1) + np.arange(30) * 1e-3
    r = S.signflip_mean_gt0(pos, 200_000, "t3")
    assert r["p"] <= 3 / 200_001, r
    # symmetric data: p about 0.5
    sym = np.concatenate([np.linspace(0.01, 0.5, 15), -np.linspace(0.01, 0.5, 15)])
    assert 0.4 < S.signflip_mean_gt0(sym, 200_000, "t4")["p"] < 1.0
    again = S.signflip_mean_gt0(big, 200_000, "t2")
    assert again["p"] == mc["p"], "must be deterministic for a fixed seed/label"
    print("sign-flip exact==brute, MC plus-one + determinism OK")


def test_perm_matches_enumeration():
    rng = np.random.default_rng(0)
    x = rng.normal(size=10)
    g = np.array([1] * 4 + [0] * 6)
    obs = x[g == 1].mean() - x[g == 0].mean()
    ds = []
    for c in itertools.combinations(range(10), 4):
        m = np.zeros(10, bool)
        m[list(c)] = True
        ds.append(x[m].mean() - x[~m].mean())
    ds = np.array(ds)
    exact_two = np.mean(np.abs(ds) >= abs(obs) - 1e-12)
    r = S.perm_diff_means(x, g, 200_000, "perm")
    assert abs(r["p_two_sided"] - exact_two) < 0.01, (r, exact_two)
    assert abs(r["p_group1_greater"] - np.mean(ds >= obs - 1e-12)) < 0.01
    print("free label permutation ~ exact enumeration OK")


def test_stratified_perm():
    # labels only swap within clusters: with perfect within-cluster separation, the within-cluster test is significant
    # and a pure between-cluster offset (confound) is NOT.
    rng = np.random.default_rng(1)
    n_cl = 40
    cl = np.repeat(np.arange(n_cl), 4)
    off = np.repeat(rng.normal(0, 3, n_cl), 4)
    g = np.tile([1, 1, 0, 0], n_cl)
    null = off + rng.normal(0, 0.1, len(cl))
    eff = null + 0.5 * g
    p_null = S.perm_diff_means(null, g, 20_000, "s0", strata=cl)["p_two_sided"]
    p_eff = S.perm_diff_means(eff, g, 20_000, "s1", strata=cl)["p_two_sided"]
    assert p_eff < 1e-3 and p_null > 0.01, (p_null, p_eff)
    print("within-cluster stratified permutation OK")


def test_bootstrap():
    x = np.random.default_rng(2).normal(0.3, 1, 400)
    lo, hi = S.bootstrap_mean_ci(x, "b")
    se = x.std(ddof=1) / 20
    assert lo < x.mean() < hi and abs((hi - lo) - 2 * 1.96 * se) < 0.4 * se
    print("bootstrap CI OK")


if __name__ == "__main__":
    for t in (test_bh, test_signflip_exact_and_mc, test_perm_matches_enumeration, test_stratified_perm, test_bootstrap):
        t()
    print("ALL STAT TESTS PASSED")
