"""Pure computation for the per-array worker (no I/O) so it can be unit-tested against brute force.

Conventions
-----------
* Alignment matrices are int8: 0..3 = A,C,G,T; 4 = ambiguous/other letter; 5 = gap.
* "Callable" means A/C/G/T. Gaps and ambiguous bases are ignored by every statistic.
* Divergence for a window = sum(differences) / sum(comparable pairs) over all its columns. Per-column
  divergence values are never averaged.
"""
import numpy as np

import pp_common as C

GAP = 5
_ENC = np.full(256, 4, dtype=np.int8)
for _i, _b in enumerate("ACGT"):
    _ENC[ord(_b)] = _i
    _ENC[ord(_b.lower())] = _i
_ENC[ord("-")] = GAP
_ENC[ord(".")] = GAP


def encode(seqs):
    L = {len(s) for s in seqs}
    if len(L) != 1:
        raise ValueError("sequences are not the same length")
    arr = np.frombuffer("".join(seqs).encode("ascii"), dtype=np.uint8).reshape(len(seqs), -1)
    return _ENC[arr]


# ------------------------------------------------------------------ phase 1: hits, chains, units
def filter_hits(hits, hmm_len):
    for h in hits:
        h["hmm_coverage"] = (h["hmm_to"] - h["hmm_from"] + 1) / hmm_len
        bad = []
        if h["hmm_coverage"] < C.MIN_HMM_COVERAGE:
            bad.append("coverage")
        if h["score"] < C.MIN_NHMMER_SCORE:
            bad.append("score")
        h["filter_status"] = "pass" if not bad else "fail_" + "_and_".join(bad)
    return hits


def dedup_hits(hits):
    """Highest score first; drop a hit overlapping a retained hit by >= 50% of the shorter one."""
    cand = sorted([h for h in hits if h["filter_status"] == "pass"],
                  key=lambda h: (-h["score"], h["evalue"], h["start"], h["end"]))
    kept = []
    for h in cand:
        h["dedup_status"], h["duplicate_of"], h["overlap_frac_shorter"] = "retained", "", 0.0
        for k in kept:
            ov = min(h["end"], k["end"]) - max(h["start"], k["start"]) + 1
            if ov <= 0:
                continue
            frac = ov / min(h["end"] - h["start"] + 1, k["end"] - k["start"] + 1)
            if frac >= C.DEDUP_OVERLAP_FRAC:
                h["dedup_status"], h["duplicate_of"], h["overlap_frac_shorter"] = "rejected_duplicate", k["hit_id"], frac
                break
        if h["dedup_status"] == "retained":
            kept.append(h)
    for h in hits:
        if h["filter_status"] != "pass":
            h["dedup_status"], h["duplicate_of"], h["overlap_frac_shorter"] = "not_evaluated", "", 0.0
    return kept


def best_chain(anchors):
    """Longest tandem chain (adjacent genomic-start spacing 3,500-20,000 bp); ties broken by summed score.

    `anchors` must be sorted by start. Returns the selected anchors (possibly a single one).
    """
    n = len(anchors)
    if n == 0:
        return []
    cnt = [1] * n
    sc = [round(a["score"], 6) for a in anchors]
    prev = [-1] * n
    for j in range(n):
        for i in range(j):
            d = anchors[j]["start"] - anchors[i]["start"]
            if C.SPACING_MIN <= d <= C.SPACING_MAX:
                cand = (cnt[i] + 1, sc[i] + round(anchors[j]["score"], 6))
                if cand > (cnt[j], sc[j]):
                    cnt[j], sc[j], prev[j] = cand[0], cand[1], i
    end = max(range(n), key=lambda j: (cnt[j], sc[j], -j))
    chain = []
    while end != -1:
        chain.append(anchors[end])
        end = prev[end]
    return chain[::-1]


def build_chains(retained):
    """Separate chain per strand; a chain needs >= 2 anchors."""
    out = {}
    for strand in "+-":
        sub = sorted([h for h in retained if h["strand"] == strand], key=lambda h: (h["start"], h["end"]))
        ch = best_chain(sub)
        if len(ch) >= 2:
            out[strand] = ch
    return out


def oriented_anchor(seq, h):
    s = seq[h["start"] - 1:h["end"]]
    return C.revcomp(s) if h["strand"] == "-" else s


def build_units(seq, chain, strand):
    """Adjacent-anchor intervals only, returned in transcription order and transcription orientation.

    '+' : [start_k, start_{k+1}) as-is.
    '-' : (end_k, end_{k+1}] reverse-complemented, so each unit begins with the reverse-complemented
          downstream-in-genome anchor (the transcriptionally current one).
    """
    units = []
    ks = range(len(chain) - 1) if strand == "+" else range(len(chain) - 2, -1, -1)
    for k in ks:
        a, b = chain[k], chain[k + 1]
        if strand == "+":
            s, e = a["start"], b["start"] - 1
            useq, first, nxt = seq[s - 1:e], a, b
        else:
            s, e = a["end"] + 1, b["end"]
            useq, first, nxt = C.revcomp(seq[s - 1:e]), b, a
        units.append(dict(strand=strand, local_start=s, local_end=e, seq=useq, start_anchor=first, next_anchor=nxt))
    return units


# ------------------------------------------------------------------ column / window statistics
def col_stats(X):
    n = X.shape[0]
    counts = np.stack([(X == b).sum(0) for b in range(4)]).astype(np.int64)
    ncall = counts.sum(0)
    pairs = ncall * (ncall - 1) // 2
    diffs = pairs - (counts * (counts - 1) // 2).sum(0)
    with np.errstate(divide="ignore", invalid="ignore"):
        p = counts / ncall
        H = -np.where(counts > 0, p * np.log2(p), 0.0).sum(0)
    H = np.nan_to_num(H)
    return dict(n=n, ncall=ncall, pairs=pairs, diffs=diffs, gc=counts[1] + counts[2],
                seg=(((counts > 0).sum(0)) >= 2).astype(np.int64), H=H)


def wsum(a, w):
    cs = np.concatenate([[0.0], np.cumsum(np.asarray(a, dtype=np.float64))])
    return cs[w:] - cs[:-w]


def _div(d, p):
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(p > 0, d / p, np.nan)


def window_metrics(cs, w=C.WIN):
    n = cs["n"]
    call = wsum(cs["ncall"], w)
    out = dict(
        divergence=_div(wsum(cs["diffs"], w), wsum(cs["pairs"], w)),
        gc=_div(wsum(cs["gc"], w), call),
        occupancy=call / (n * w),
        segregating_sites=wsum(cs["seg"], w),
        mean_entropy_bits=wsum(cs["H"], w) / w,
    )
    col_ok = cs["ncall"] >= C.MIN_CALLABLE_FRAC * n - 1e-9
    out["eligible"] = wsum((~col_ok).astype(float), w) == 0
    return out


def cpg_gapless(X, w=C.WIN):
    """Primary CpG: per sequence, drop gaps from the w-column fragment, then pool CpG / callable dinucleotides."""
    n, L = X.shape
    nw = L - w + 1
    cpg, dn = np.zeros(nw), np.zeros(nw)
    for s in range(nw):
        blk = X[:, s:s + w]
        for i in range(n):
            f = blk[i]
            f = f[f != GAP]
            if f.size < 2:
                continue
            a, b = f[:-1], f[1:]
            dn[s] += np.count_nonzero((a < 4) & (b < 4))
            cpg[s] += np.count_nonzero((a == 1) & (b == 2))
    return _div(cpg, dn)


def cpg_adjacent(X, w=C.WIN):
    """Backbone CpG: directly adjacent aligned callable dinucleotides within the window."""
    call = X < 4
    adj = (call[:, :-1] & call[:, 1:]).sum(0)
    cg = ((X[:, :-1] == 1) & (X[:, 1:] == 2)).sum(0)
    return _div(wsum(cg, w - 1), wsum(adj, w - 1))


def empirical_p(promoter_div, control_divs):
    c = np.asarray(control_divs, dtype=float)
    return (1 + int(np.sum(c <= promoter_div + 1e-12))) / (1 + len(c))


def summarise_vs_controls(promoter_div, control_divs):
    c = np.asarray(control_divs, dtype=float)
    n = len(c)
    if n == 0 or not np.isfinite(promoter_div):
        return dict(n_controls=n, control_median=np.nan, effect=np.nan, ratio=np.nan, p=np.nan,
                    p_min=(1 / (1 + n)) if n else np.nan)
    med = float(np.median(c))
    return dict(n_controls=n, control_median=med, effect=med - promoter_div,
                ratio=(promoter_div / med) if med > 0 else np.nan,
                p=empirical_p(promoter_div, c), p_min=1 / (1 + n))


def match(gc, cpg, occ, ref, tol_cpg):
    """Composition match to the promoter reference tuple (gc, cpg, occ)."""
    return (np.abs(gc - ref[0]) <= C.TOL_GC + 1e-12) & (np.abs(cpg - ref[1]) <= tol_cpg + 1e-12) & \
           (np.abs(occ - ref[2]) <= C.TOL_OCC + 1e-12)


def primary_windows(X):
    """Every 50-column window (starts 1..251) of the 300-column alignment with the matched-window test."""
    cs = col_stats(X)
    wm = window_metrics(cs)
    cpg = cpg_gapless(X)
    nw = len(wm["divergence"])
    start1 = np.arange(1, nw + 1)
    end1 = start1 + C.WIN - 1
    p0, p1 = C.PROMOTER_300
    is_prom = start1 == p0
    overlaps = ~((end1 < p0) | (start1 > p1)) & ~is_prom
    ps = p0 - 1
    ref = (wm["gc"][ps], cpg[ps], wm["occupancy"][ps])
    ctrl = (~is_prom) & (~overlaps) & wm["eligible"] & np.isfinite(wm["divergence"])
    with np.errstate(invalid="ignore"):
        m_primary = ctrl & match(wm["gc"], cpg, wm["occupancy"], ref, C.TOL_CPG)
        m_grid = ctrl & match(wm["gc"], cpg, wm["occupancy"], ref, C.TOL_CPG_GRID)
    cls = np.where(is_prom, "promoter", np.where(overlaps, "overlaps_promoter", "control"))
    table = dict(start_1based=start1, end_1based=end1, window_class=cls, eligible=wm["eligible"],
                 divergence=wm["divergence"], gc=wm["gc"], cpg=cpg, occupancy=wm["occupancy"],
                 segregating_sites=wm["segregating_sites"], mean_entropy_bits=wm["mean_entropy_bits"],
                 control_candidate=ctrl, matched_primary=m_primary, matched_cpg_grid_sensitivity=m_grid)
    return cs, table, ps, ref


def rolling_pooled(diffs, pairs, half):
    n = len(diffs)
    cd = np.concatenate([[0.0], np.cumsum(diffs, dtype=np.float64)])
    cp = np.concatenate([[0.0], np.cumsum(pairs, dtype=np.float64)])
    lo, hi = np.clip(np.arange(n) - half, 0, n), np.clip(np.arange(n) + half + 1, 0, n)
    return _div(cd[hi] - cd[lo], cp[hi] - cp[lo])


def rolling_mean(v, half):
    n = len(v)
    cv = np.concatenate([[0.0], np.cumsum(v, dtype=np.float64)])
    lo, hi = np.clip(np.arange(n) - half, 0, n), np.clip(np.arange(n) + half + 1, 0, n)
    return (cv[hi] - cv[lo]) / (hi - lo)


# ------------------------------------------------------------------ phase 4: conserved backbone
def backbone_windows(U, mask):
    """Matched distal-backbone test on the 9,090-column unit alignment (promoter = columns 74-123)."""
    cs = col_stats(U)
    wm = window_metrics(cs)
    cpg = cpg_adjacent(U)
    nw = len(wm["divergence"])
    start1 = np.arange(1, nw + 1)
    end1 = start1 + C.WIN - 1
    p0, p1 = C.PROMOTER_FULL
    is_prom = start1 == p0
    overlaps = ~((end1 < p0) | (start1 > p1)) & ~is_prom
    ps = p0 - 1
    ref = (wm["gc"][ps], cpg[ps], wm["occupancy"][ps])

    def all_true(flag):
        return wsum((~flag).astype(float), C.WIN) == 0

    prim_ok = all_true(mask["primary_distal_backbone_ok"].to_numpy())
    sens_ok = all_true(mask["all_non_acro1_backbone_ok"].to_numpy())
    base = wm["eligible"] & np.isfinite(wm["divergence"]) & ~is_prom & ~overlaps
    cand_p, cand_s = base & prim_ok, base & sens_ok
    with np.errstate(invalid="ignore"):
        mt = match(wm["gc"], cpg, wm["occupancy"], ref, C.BACKBONE_TOL_CPG)
    table = dict(start_1based=start1, end_1based=end1, window_class=np.where(is_prom, "promoter", np.where(overlaps, "overlaps_promoter", "other")),
                 eligible=wm["eligible"], divergence=wm["divergence"], gc=wm["gc"], cpg=cpg,
                 occupancy=wm["occupancy"], primary_distal_candidate=cand_p, all_non_acro1_candidate=cand_s,
                 matched_primary_distal=cand_p & mt, matched_all_non_acro1=cand_s & mt)
    return cs, table, ps, ref

