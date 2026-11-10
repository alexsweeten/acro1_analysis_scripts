#!/usr/bin/env python3
"""EXPLORATORY, POST-HOC (not prespecified): does the group contrast survive conditioning on overall divergence level?

Motivation: the absolute effect (control median - promoter divergence) is bounded by the array's overall divergence,
and the two groups differ ~10x in median control divergence. This does not replace any prespecified result.

Arrays are binned into quartiles of control-median divergence (estimable arrays with control median > 0). Reports per-stratum
group means / ratios and a label permutation of the group contrast conditional on stratum (labels shuffled only within strata).
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pp_common as C
import pp_stats as S

OUT = C.COMBINED / "exploratory_divergence_strata.tsv"
DRAWS = 200_000


def main():
    r = pd.read_csv(C.COMBINED / "all_array_results.tsv", sep="\t", dtype={"haplotype": str})
    est = r[(r.primary_status == "estimable") & (r.control_median_divergence > 0)].copy()
    est["stratum"] = pd.qcut(est.control_median_divergence, 4, labels=False)
    bounds = est.groupby("stratum").control_median_divergence.agg(["min", "max"])
    rows = []
    for (s, g), sub in est.groupby(["stratum", "recombination_class"]):
        rows.append(dict(row_type="stratum_group", stratum=int(s) + 1, control_median_min=bounds.at[s, "min"], control_median_max=bounds.at[s, "max"],
                         group=g, n=len(sub), mean_effect=sub.effect.mean(), median_ratio=sub.ratio.median(), frac_promoter_less_divergent=(sub.ratio < 1).mean(),
                         frac_promoter_more_divergent=(sub.ratio > 1).mean()))
    for g, sub in est.groupby("recombination_class"):
        rows.append(dict(row_type="all_informative", group=g, n=len(sub), mean_effect=sub.effect.mean(), median_ratio=sub.ratio.median(),
                         frac_promoter_less_divergent=(sub.ratio < 1).mean(), frac_promoter_more_divergent=(sub.ratio > 1).mean()))
    a, b = "telomeric_recombination", "telomeric_no_recombination"
    sel = est[est.recombination_class.isin([a, b])]
    raw = S.perm_diff_means(sel.effect.to_numpy(), (sel.recombination_class == a).to_numpy(), DRAWS, "exploratory|effect|unconditioned")
    # Stratum-adjusted difference T = sum_s w_s (mean_A,s - mean_B,s) / sum_s w_s, w_s = nA*nB/(nA+nB); null: labels shuffled within strata.
    strata = []
    for s, sub in sel.groupby("stratum"):
        x = sub.effect.to_numpy()
        g = (sub.recombination_class == a).to_numpy()
        if g.sum() and (~g).sum():
            strata.append((x, int(g.sum()), len(x) - int(g.sum()), g))
    w = np.array([na * nb / (na + nb) for _, na, nb, _ in strata])
    obs = float(sum(wi * (x[g].mean() - x[~g].mean()) for wi, (x, na, nb, g) in zip(w, strata)) / w.sum())
    rng = S.rng_for("exploratory|effect|stratum-adjusted")
    T = np.zeros(DRAWS)
    for wi, (x, na, nb, g) in zip(w, strata):
        done = 0
        step = max(1, 8_000_000 // len(x))
        tot = x.sum()
        while done < DRAWS:
            m = min(step, DRAWS - done)
            idx = np.broadcast_to(np.arange(len(x)), (m, len(x))).copy()
            rng.permuted(idx, axis=1, out=idx)
            sa = x[idx[:, :na]].sum(1)
            T[done:done + m] += wi * (sa / na - (tot - sa) / nb)
            done += m
    T /= w.sum()
    tol = 1e-12
    p_two = (1 + int((np.abs(T) >= abs(obs) - tol).sum())) / (DRAWS + 1)
    p_a_less = (1 + int((T <= obs + tol).sum())) / (DRAWS + 1)
    per = [(x[g].mean() - x[~g].mean()) for x, na, nb, g in strata]
    rows.append(dict(row_type="stratum_adjusted_contrast", group=f"{a} - {b}", n=len(sel), mean_effect=obs, p_two_sided=p_two, p_a_less_than_b=p_a_less,
                     permutations=DRAWS, note="CMH-style weights nA*nB/(nA+nB); labels permuted within divergence strata; EXPLORATORY, POST-HOC; per-stratum differences: "
                     + ", ".join(f"{v:+.4f}" for v in per)))
    rows.append(dict(row_type="unadjusted_contrast_same_arrays", group=f"{a} - {b}", n=len(sel), mean_effect=raw["diff"], p_two_sided=raw["p_two_sided"],
                     permutations=DRAWS, note="same arrays (control median > 0), free label permutation, for comparison"))
    pd.DataFrame(rows).to_csv(OUT, sep="\t", index=False, float_format="%.6g")
    print(pd.DataFrame(rows).drop(columns=["note"], errors="ignore").round(5).to_string(index=False))
    print("per-stratum differences (A - B):", [round(v, 5) for v in per])


if __name__ == "__main__":
    main()
