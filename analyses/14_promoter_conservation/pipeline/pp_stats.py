"""Randomisation inference. All routines are deterministic given (seed, label)."""
import itertools
import zlib

import numpy as np

import pp_common as C

CHUNK_ELEMS = 12_000_000     # bounds memory of the (draws x n) matrices


def rng_for(label, seed=C.SEED):
    """Independent, reproducible stream per test: SeedSequence([seed, crc32(label)])."""
    return np.random.default_rng(np.random.SeedSequence([seed, zlib.crc32(label.encode())]))


def bh(p):
    """Benjamini-Hochberg adjusted p-values (NaN-aware)."""
    p = np.asarray(p, dtype=float)
    q = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return q
    pv = p[ok]
    order = np.argsort(pv, kind="mergesort")
    ranked = pv[order] * len(pv) / (np.arange(len(pv)) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty_like(pv)
    out[order] = np.minimum(ranked, 1.0)
    q[ok] = out
    return q


def signflip_mean_gt0(x, draws, label):
    """One-sided H1: mean(x) > 0. Exact if n <= 20; else Monte Carlo with plus-one correction."""
    x = np.asarray(x, dtype=float)
    n = len(x)
    obs = x.mean()
    tol = 1e-12 * max(1.0, abs(obs))
    if n == 0:
        return dict(p=np.nan, method="none", draws=0)
    if n <= C.EXACT_MAX_N:
        signs = np.array(list(itertools.product((-1.0, 1.0), repeat=n)))
        extreme = int(((signs @ x) / n >= obs - tol).sum())
        return dict(p=extreme / 2 ** n, method="exact", draws=2 ** n)
    rng = rng_for(label)
    step = max(1, CHUNK_ELEMS // n)
    extreme, done = 0, 0
    while done < draws:
        b = min(step, draws - done)
        s = rng.integers(0, 2, size=(b, n), dtype=np.int8).astype(np.float64) * 2 - 1
        extreme += int(((s @ x) / n >= obs - tol).sum())
        done += b
    return dict(p=(extreme + 1) / (draws + 1), method="monte_carlo", draws=draws)


def bootstrap_mean_ci(x, label, nboot=C.N_BOOT, alpha=0.05):
    x = np.asarray(x, dtype=float)
    n = len(x)
    if n < 2:
        return (np.nan, np.nan)
    rng = rng_for(label + "|boot")
    step = max(1, CHUNK_ELEMS // n)
    means = []
    done = 0
    while done < nboot:
        b = min(step, nboot - done)
        means.append(x[rng.integers(0, n, size=(b, n))].mean(1))
        done += b
    m = np.concatenate(means)
    return (float(np.quantile(m, alpha / 2)), float(np.quantile(m, 1 - alpha / 2)))


def perm_diff_means(x, g, draws, label, strata=None):
    """Permute group labels; statistic = mean(x[g==1]) - mean(x[g==0]).

    strata=None: free permutation of array labels.  strata=int array: labels are shuffled only within strata
    (clusters), which conditions on cluster membership. Plus-one corrected p-values.
    """
    x = np.asarray(x, dtype=float)
    g = np.asarray(g).astype(bool)
    n, n1 = len(x), int(g.sum())
    n0 = n - n1
    obs = x[g].mean() - x[~g].mean()
    tol = 1e-12 * max(1.0, abs(obs))
    rng = rng_for(label)
    step = max(1, CHUNK_ELEMS // n)
    tot = x.sum()
    ge = le = ab = 0
    done = 0
    if strata is not None:
        order0 = np.argsort(strata, kind="mergesort")
        x, g = x[order0], g[order0]
        strata = np.asarray(strata)[order0]
        sid = np.unique(strata, return_inverse=True)[1].astype(np.float64)
        gf = g.astype(np.float64)
    while done < draws:
        b = min(step, draws - done)
        if strata is None:
            idx = np.broadcast_to(np.arange(n), (b, n)).copy()
            rng.permuted(idx, axis=1, out=idx)
            s1 = x[idx[:, :n1]].sum(1)
        else:
            perm = np.argsort(rng.random((b, n)) + sid[None, :], axis=1)
            s1 = (gf[perm] * x[None, :]).sum(1)
        d = s1 / n1 - (tot - s1) / n0
        ge += int((d >= obs - tol).sum())
        le += int((d <= obs + tol).sum())
        ab += int((np.abs(d) >= abs(obs) - tol).sum())
        done += b
    f = lambda k: (k + 1) / (draws + 1)
    return dict(diff=float(obs), p_two_sided=f(ab), p_group1_greater=f(ge), p_group1_less=f(le), draws=draws)
