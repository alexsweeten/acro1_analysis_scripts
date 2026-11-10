#!/usr/bin/env python3
"""Aggregate per-array outputs into combined/. The ARRAY is the unit of inference throughout.

  aggregate.py [--outdir DIR] [--subset FILE] [--allow-partial] [--skip-r]

Without --allow-partial every manifest array must have a COMPLETE sentinel (a valid result or a documented
biological non-estimability), and all arrays must come from the same pipeline_version.
"""
import argparse
import subprocess
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pp_common as C
import pp_stats as S

HERE = Path(__file__).resolve().parent
GROUP_ORDER = ["telomeric_recombination", "telomeric_no_recombination", "non_telomeric_no_recombination"]

# (name, tier, effect column, p column, status column, extra mask fn, draws)
ANALYSES = [
    ("primary_300col_matched", "PRIMARY", "effect", "p_matched", "primary_status", None, C.N_SIGNFLIP_PRIMARY),
    ("backbone_distal_matched", "BACKBONE", "bb_effect", "bb_p_matched", "bb_status", None, C.N_SIGNFLIP_BACKBONE),
    ("primary_allcontrols_unmatched", "SENSITIVITY", "allctrl_unmatched_sens_effect", None, None, None, C.N_SIGNFLIP_PRIMARY),
    ("primary_cpg_tol_0.021", "SENSITIVITY", "cpg0.021_sens_effect", None, None, None, C.N_SIGNFLIP_PRIMARY),
    ("primary_excluding_mixed_strand", "SENSITIVITY", "effect", None, "primary_status", lambda d: ~d.mixed_strand.astype(bool), C.N_SIGNFLIP_PRIMARY),
    ("backbone_all_non_ACRO1", "SENSITIVITY", "bb_allnonacro1_sens_effect", None, None, None, C.N_SIGNFLIP_BACKBONE),
    ("backbone_excluding_mixed_strand", "SENSITIVITY", "bb_effect", None, "bb_status", lambda d: ~d.mixed_strand.astype(bool), C.N_SIGNFLIP_BACKBONE),
    ("backbone_excluding_mapping_flag", "SENSITIVITY", "bb_effect", None, "bb_status",
     lambda d: d.bb_review_flags.fillna("").eq(""), C.N_SIGNFLIP_BACKBONE),
]
CORE = ("primary_300col_matched", "backbone_distal_matched")


def load(args):
    man = pd.read_csv(C.ROOT / "manifest.validated.tsv", sep="\t", keep_default_na=False, dtype={"haplotype": str})
    if args.subset:
        keep = {int(x) for x in open(args.subset).read().split()}
        man = man[man.task_index.isin(keep)]
    rows, missing = [], []
    for _, m in man.iterrows():
        d = C.PER_ARRAY / m.safe_id
        if not (d / "COMPLETE").exists():
            missing.append(m.array_id)
            continue
        rows.append(pd.read_csv(d / "result.tsv", sep="\t", keep_default_na=True, dtype={"haplotype": str}))
    if missing and not args.allow_partial:
        sys.exit(f"{len(missing)} manifest arrays lack COMPLETE (e.g. {missing[:3]}); rerun them or use --allow-partial")
    res = pd.concat(rows, ignore_index=True)
    vers = res.pipeline_version.unique()
    if len(vers) != 1:
        sys.exit(f"mixed pipeline versions in results: {vers}")
    res = man[["task_index", "array_id"]].merge(res, on="array_id", how="inner")
    return man, res, missing


def apply_mapping_recheck(res, outdir):
    """Replace the stored promoter-mapping QC by the raw-coordinate recheck (stored value used aligned-ungapped indices).

    Only QC columns change; no divergence, effect or p-value depends on them. Original values are kept as *_aligned_index.
    """
    f = outdir / "promoter_mapping_recheck.tsv"
    if not f.exists():
        sys.exit(f"{f} missing: run scripts/recheck_mapping.py first")
    rc = pd.read_csv(f, sep="\t")
    if (rc.recheck_status != "ok").any():
        sys.exit("recheck contains arrays without a backbone mapping")
    res = res.merge(rc[["array_id", "recheck_concordance", "units_raw_prefix_violation", "units_with_different_promoter_base_sets"]],
                    on="array_id", how="left")
    res["bb_mapping_concordance_aligned_index"] = res["bb_mapping_concordance"]
    res["bb_anchor_not_unit_prefix_aligned_index"] = res["bb_anchor_not_unit_prefix"]
    res["bb_mapping_concordance"] = res["recheck_concordance"]
    res["bb_anchor_not_unit_prefix"] = res["units_raw_prefix_violation"]
    res["bb_units_with_different_promoter_base_sets"] = res["units_with_different_promoter_base_sets"]
    low = res["bb_mapping_concordance"] < C.MAPPING_CONCORDANCE_MIN
    res["bb_review_flags"] = np.where(low, "promoter_mapping_concordance_below_0.99", "")
    drop = {"promoter_mapping_concordance_below_0.99", "anchor_not_unit_prefix"}
    keep = res["qc_flags"].fillna("").map(lambda s: [x for x in s.split(";") if x and x not in drop])
    res["qc_flags"] = [";".join(k + (["promoter_mapping_concordance_below_0.99"] if l else [])) for k, l in zip(keep, low)]
    return res.drop(columns=["recheck_concordance", "units_raw_prefix_violation", "units_with_different_promoter_base_sets"])


def with_bh(res):
    res = res.copy()
    for pre, status, pcol in (("", "primary_status", "p_matched"), ("bb_", "bb_status", "bb_p_matched")):
        ok = res[status].eq("estimable")
        res[f"{pre}q_bh_within_group"] = np.nan
        res[f"{pre}q_bh_all_estimable"] = np.nan
        res.loc[ok, f"{pre}q_bh_all_estimable"] = S.bh(res.loc[ok, pcol].to_numpy())
        for g, idx in res[ok].groupby("recombination_class").groups.items():
            res.loc[idx, f"{pre}q_bh_within_group"] = S.bh(res.loc[idx, pcol].to_numpy())
    return res


def group_masks(df):
    masks = {}
    for g in GROUP_ORDER:
        mk = df.recombination_class.eq(g)
        if mk.any():
            masks[g] = mk
    return masks


def vec(df, col, status, extra):
    m = df[col].notna()
    if status:
        m &= df[status].eq("estimable")
    if extra is not None:
        m &= extra(df)
    return m


def summarise(x, draws, label):
    x = np.asarray(x, dtype=float)
    if len(x) == 0:
        return dict(n_estimable=0)
    sf = S.signflip_mean_gt0(x, draws, label)
    lo, hi = S.bootstrap_mean_ci(x, label)
    return dict(n_estimable=len(x), mean_effect=x.mean(), median_effect=float(np.median(x)), n_positive=int((x > 0).sum()),
                frac_positive=float((x > 0).mean()), boot_ci95_lo=lo, boot_ci95_hi=hi, signflip_p_one_sided_mean_gt0=sf["p"],
                signflip_method=sf["method"], signflip_draws=sf["draws"], seed_base=C.SEED, seed_stream=f"crc32('{label}')")


def array_weighted_group_summary(df):
    masks = group_masks(df)
    out = []
    for name, tier, ecol, pcol, scol, extra, draws in ANALYSES:
        v = vec(df, ecol, scol, extra)
        for g, gm in list(masks.items()) + [("ALL", pd.Series(True, index=df.index))]:
            sel = df[gm & v]
            row = dict(analysis=name, tier=tier, group=g, weighting="array", n_manifest_arrays=int(gm.sum()))
            row.update(summarise(sel[ecol].to_numpy(), draws, f"{name}|{g}|array"))
            row["n_nonestimable"] = row["n_manifest_arrays"] - row["n_estimable"]
            out.append(row)
    gs = pd.DataFrame(out)
    gs["q_bh_across_group_tests"] = np.nan
    for name in gs.analysis.unique():
        k = (gs.analysis == name) & (gs.group != "ALL")
        gs.loc[k, "q_bh_across_group_tests"] = S.bh(gs.loc[k, "signflip_p_one_sided_mean_gt0"].to_numpy())
    return gs


def pairwise_tests(df):
    masks = group_masks(df)
    names = [g for g in GROUP_ORDER if g in masks]
    pairs = [(a, b) for i, a in enumerate(names) for b in names[i + 1:]]
    out = []
    for name, tier, ecol, pcol, scol, extra, draws in ANALYSES[:2]:
        ok = vec(df, ecol, scol, extra)
        for a, b in pairs:
            sel = df[ok & (masks[a] | masks[b])]
            x, g = sel[ecol].to_numpy(), masks[a][sel.index].to_numpy()
            r = S.perm_diff_means(x, g, C.N_GROUP_PERM, f"{name}|{a}-vs-{b}|array")
            out.append(dict(analysis=name, group_a=a, group_b=b, n_a=int(g.sum()), n_b=int((~g).sum()),
                            mean_a=x[g].mean(), mean_b=x[~g].mean(), diff_a_minus_b=r["diff"],
                            p_two_sided=r["p_two_sided"], p_a_greater_than_b=r["p_group1_greater"], p_a_less_than_b=r["p_group1_less"],
                            permutations=r["draws"], seed_base=C.SEED,
                            directional_note="no prespecified directional hypothesis supplied; one-sided tails shown descriptively"))
    pt = pd.DataFrame(out)
    pt["q_bh_across_pairwise"] = np.nan
    for name in pt.analysis.unique() if len(pt) else []:
        k = pt.analysis == name
        pt.loc[k, "q_bh_across_pairwise"] = S.bh(pt.loc[k, "p_two_sided"].to_numpy())
    return pt


def cluster_aware(df, skip_r, outdir):
    masks = group_masks(df)
    out, model_in = [], []
    for name, tier, ecol, pcol, scol, extra, draws in ANALYSES[:2]:
        ok = vec(df, ecol, scol, extra)
        for unit, key in (("individual", df.individual_id), ("individual_haplotype", df.individual_id + "#" + df.haplotype.astype(str))):
            for g, gm in list(masks.items()) + [("ALL", pd.Series(True, index=df.index))]:
                sel = df[gm & ok]
                cm = sel.groupby(key[sel.index])[ecol].mean()
                row = dict(analysis=name, test="cluster_mean_signflip", cluster_unit=unit, group=g, n_arrays=len(sel), n_clusters=len(cm))
                row.update(summarise(cm.to_numpy(), draws, f"{name}|{g}|{unit}"))
                row["estimate"], row["detail"] = row.get("mean_effect"), "mean of cluster-mean effects"
                out.append(row)
        names = [g for g in GROUP_ORDER if g in masks]
        for i, a in enumerate(names):
            for b in names[i + 1:]:
                sel = df[ok & (masks[a] | masks[b])]
                x, g = sel[ecol].to_numpy(), masks[a][sel.index].to_numpy()
                ind = pd.factorize(sel.individual_id)[0]
                tab = pd.crosstab(ind, g)
                both = int(((tab > 0).all(axis=1)).sum())
                r = S.perm_diff_means(x, g, C.N_GROUP_PERM, f"{name}|{a}-vs-{b}|within-individual", strata=ind)
                out.append(dict(analysis=name, test="group_contrast_within_individual_stratified_permutation", cluster_unit="individual",
                                group=f"{a} - {b}", n_arrays=len(sel), n_clusters=int(tab.shape[0]), n_clusters_with_both_groups=both,
                                estimate=r["diff"], p_two_sided=r["p_two_sided"], p_a_greater_than_b=r["p_group1_greater"],
                                p_a_less_than_b=r["p_group1_less"], signflip_draws=r["draws"], seed_base=C.SEED,
                                detail="labels permuted only among arrays of the same individual; array-weighted difference of means"))
                model_in.append(pd.DataFrame(dict(analysis=name, array_id=sel.array_id, individual_id=sel.individual_id.to_numpy(),
                                                  group=np.where(g, a, b), reference=b, effect=x)))
    cl = pd.DataFrame(out)
    if model_in and not skip_r:
        mi = pd.concat(model_in)
        fin = outdir / "_cluster_model_input.tsv"
        mi.to_csv(fin, sep="\t", index=False)
        fout = outdir / "_cluster_model_output.tsv"
        rc = subprocess.run(["Rscript", str(HERE / "cluster_models.R"), str(fin), str(fout), str(C.SEED), str(C.N_CLUSTER_BOOT)],
                            capture_output=True, text=True)
        (outdir / "_cluster_models.log").write_text(rc.stdout + "\n" + rc.stderr)
        if rc.returncode:
            sys.exit(f"cluster_models.R failed; see {outdir}/_cluster_models.log")
        cl = pd.concat([cl, pd.read_csv(fout, sep="\t")], ignore_index=True)
        fin.unlink()
        fout.unlink()
    return cl


def position_tables(df, outdir):
    np.seterr(all="ignore")
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    n = len(df)
    div = np.full((n, C.ALN_COLS), np.nan)
    roll = np.full((n, C.ALN_COLS), np.nan)
    cnt = np.zeros((n, C.ALN_COLS))
    occ = np.full((n, C.ALN_COLS), np.nan)
    diffs = np.zeros(C.ALN_COLS)
    pairs = np.zeros(C.ALN_COLS)
    pro = None
    for i, sid in enumerate(df.safe_id):
        p = pd.read_csv(C.PER_ARRAY / sid / "per_position_profile.tsv", sep="\t")
        if len(p) != C.ALN_COLS:
            continue                      # non-estimable arrays never reached the alignment
        div[i], roll[i] = p.pairwise_divergence, p.rolling_11col_divergence
        cnt[i], occ[i] = p.anchor_count, p.occupancy
        diffs += p.pairwise_differences.fillna(0).to_numpy()
        pairs += p.comparable_pairs.fillna(0).to_numpy()
        pro = p
    callable_ = occ >= C.MIN_CALLABLE_FRAC
    dm = np.where(callable_, div, np.nan)
    rm = np.where(callable_, roll, np.nan)
    ref = pd.read_csv(C.PROSEQ, sep="\t")
    pos = np.arange(1, C.ALN_COLS + 1)
    with np.errstate(all="ignore"):
        t = pd.DataFrame(dict(
            alignment_column_1based=pos, n_arrays_callable_ge80pct=callable_.sum(0), n_arrays_with_alignment=int((~np.isnan(occ[:, 0])).sum()),
            mean_anchor_count=np.nanmean(cnt, 0), mean_occupancy=np.nanmean(occ, 0),
            array_mean_divergence=np.nanmean(dm, 0), array_median_divergence=np.nanmedian(dm, 0),
            array_q25_divergence=np.nanquantile(dm, .25, 0), array_q75_divergence=np.nanquantile(dm, .75, 0),
            array_mean_rolling11_divergence=np.nanmean(rm, 0), array_median_rolling11_divergence=np.nanmedian(rm, 0),
            projected_chr21_NEG_PROseq=ref.PROseq_mean, rolling_7col_PROseq=pro.rolling_7col_PROseq if pro is not None else np.nan,
            promoter_window=np.where((pos >= C.PROMOTER_300[0]) & (pos <= C.PROMOTER_300[1]), "yes", "no"),
            pooled_copy_divergence_NONINDEPENDENT_descriptive=np.where(pairs > 0, diffs / np.maximum(pairs, 1), np.nan)))
    t.to_csv(outdir / "per_position_array_weighted_profile.tsv", sep="\t", index=False, float_format="%.8g")
    hm = pd.DataFrame(roll, columns=[f"c{c}" for c in pos])
    hm.insert(0, "array_id", df.array_id.to_numpy())
    hm.to_csv(outdir / "array_by_position_rolling11_divergence.tsv.gz", sep="\t", index=False, float_format="%.5g")


def qc_tables(df, outdir):
    units, arows = [], []
    for _, r in df.iterrows():
        d = C.PER_ARRAY / r.safe_id
        u = pd.read_csv(d / "units_qc.tsv", sep="\t")
        if len(u):
            units.append(u)
        a = pd.read_csv(d / "l1pa10_anchors.tsv", sep="\t")
        sp = a.spacing_to_next.replace("", np.nan).astype(float).dropna() if len(a) else pd.Series(dtype=float)
        arows.append(dict(array_id=r.array_id, safe_id=r.safe_id, fasta_length=r.fasta_length, n_hits_total=r.n_hits_total,
                          n_hits_pass_filter=r.n_hits_pass_filter, n_hits_rejected_duplicate=r.n_hits_rejected_duplicate,
                          n_anchors=r.n_anchors, n_anchors_plus=r.n_anchors_plus, n_anchors_minus=r.n_anchors_minus,
                          mixed_strand=r.mixed_strand, n_units_adjacent=r.n_units_adjacent, n_units_complete=r.n_units_complete,
                          n_units_irregular=r.n_units_irregular, spacing_min=sp.min() if len(sp) else np.nan,
                          spacing_median=sp.median() if len(sp) else np.nan, spacing_max=sp.max() if len(sp) else np.nan,
                          mean_anchor_score=a.score.mean() if len(a) else np.nan, min_anchor_hmm_coverage=a.hmm_coverage.min() if len(a) else np.nan,
                          primary_status=r.primary_status, primary_reason=r.primary_reason, bb_status=r.bb_status, bb_reason=r.bb_reason,
                          qc_flags=r.qc_flags))
    pd.DataFrame(arows).to_csv(outdir / "anchor_qc.tsv", sep="\t", index=False)
    pd.concat(units, ignore_index=True).to_csv(outdir / "unit_qc.tsv", sep="\t", index=False) if units else None
    rows = []
    for a, rcol, scol, reasoncol in (("primary_300col", "primary_status", "primary_status", "primary_reason"),
                                      ("backbone", "bb_status", "bb_status", "bb_reason")):
        ne = df[df[scol].ne("estimable")]
        for _, r in ne.iterrows():
            rows.append(dict(array_id=r.array_id, analysis=a, reason_code=r[reasoncol], recombination_class=r.recombination_class,
                             individual_id=r.individual_id, fasta_length=r.fasta_length, n_anchors=r.n_anchors,
                             n_units_complete=r.n_units_complete, promoter_occupancy=r.promoter_occupancy,
                             promoter_divergence=r.promoter_divergence, n_matched_controls=r.n_matched_controls,
                             cpg0_021_sensitivity_n_controls=r["cpg0.021_sens_n_controls"],
                             cpg0_021_sensitivity_effect=r["cpg0.021_sens_effect"]))
    pd.DataFrame(rows, columns=["array_id", "analysis", "reason_code", "recombination_class", "individual_id", "fasta_length",
                                "n_anchors", "n_units_complete", "promoter_occupancy", "promoter_divergence", "n_matched_controls",
                                "cpg0_021_sensitivity_n_controls", "cpg0_021_sensitivity_effect"]).to_csv(
        outdir / "non_estimable_arrays.tsv", sep="\t", index=False)
    flagged = df[df.qc_flags.fillna("").ne("")][["array_id", "recombination_class", "qc_flags", "n_anchors", "n_units_complete",
                                                  "n_units_irregular", "bb_mapping_concordance", "fasta_length"]]
    flagged.to_csv(outdir / "qc_review_list.tsv", sep="\t", index=False)


def parameters(outdir):
    p = [("promoter_300col_start", C.PROMOTER_300[0]), ("promoter_300col_end", C.PROMOTER_300[1]), ("promoter_fullunit_start", C.PROMOTER_FULL[0]),
         ("promoter_fullunit_end", C.PROMOTER_FULL[1]), ("window_columns", C.WIN), ("nhmmer_E", C.NHMMER_E),
         ("min_hmm_coverage", C.MIN_HMM_COVERAGE), ("min_nhmmer_score", C.MIN_NHMMER_SCORE), ("dedup_overlap_frac_of_shorter", C.DEDUP_OVERLAP_FRAC),
         ("anchor_spacing_min", C.SPACING_MIN), ("anchor_spacing_max", C.SPACING_MAX), ("unit_length_min", C.UNIT_MIN), ("unit_length_max", C.UNIT_MAX),
         ("min_callable_fraction_per_column", C.MIN_CALLABLE_FRAC), ("tol_gc", C.TOL_GC), ("tol_cpg_primary", C.TOL_CPG),
         ("tol_cpg_grid_sensitivity", C.TOL_CPG_GRID), ("tol_occupancy", C.TOL_OCC), ("rolling_divergence_columns", C.ROLL_DIV),
         ("rolling_proseq_columns", C.ROLL_PRO), ("backbone_L1PA10_endpoint_column", C.L1PA10_END_FULL), ("backbone_acro1_buffer", C.ACRO1_BUFFER),
         ("backbone_tol_cpg", C.BACKBONE_TOL_CPG), ("backbone_min_units", C.MIN_UNITS_BACKBONE), ("mapping_concordance_min", C.MAPPING_CONCORDANCE_MIN),
         ("seed", C.SEED), ("signflip_draws_primary", C.N_SIGNFLIP_PRIMARY), ("signflip_draws_backbone", C.N_SIGNFLIP_BACKBONE),
         ("group_permutation_draws", C.N_GROUP_PERM), ("bootstrap_draws_mean_ci", C.N_BOOT), ("cluster_bootstrap_draws", C.N_CLUSTER_BOOT),
         ("signflip_exact_enumeration_max_n", C.EXACT_MAX_N)]
    pd.DataFrame(p, columns=["parameter", "value"]).to_csv(outdir / "analysis_parameters.tsv", sep="\t", index=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=str(C.COMBINED))
    ap.add_argument("--subset")
    ap.add_argument("--allow-partial", action="store_true")
    ap.add_argument("--skip-r", action="store_true")
    args = ap.parse_args()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    (out / "figures").mkdir(exist_ok=True)

    man, res, missing = load(args)
    res = apply_mapping_recheck(res, out)
    res = with_bh(res)
    res.to_csv(out / "all_array_results.tsv", sep="\t", index=False, float_format="%.10g")
    counts = dict(arrays_submitted=len(man), arrays_completed=len(res), arrays_primary_estimable=int(res.primary_status.eq("estimable").sum()),
                  arrays_primary_non_estimable=int(res.primary_status.ne("estimable").sum()),
                  arrays_backbone_estimable=int(res.bb_status.eq("estimable").sum()), arrays_backbone_non_estimable=int(res.bb_status.ne("estimable").sum()),
                  total_L1PA10_anchors=int(res.n_anchors.sum()), total_complete_units=int(res.n_units_complete.sum()),
                  total_adjacent_units=int(res.n_units_adjacent.sum()), arrays_mixed_strand=int(res.mixed_strand.astype(bool).sum()),
                  pipeline_version=res.pipeline_version.iloc[0])
    pd.Series(counts).rename("value").to_frame().to_csv(out / "run_counts.tsv", sep="\t", index_label="metric")
    print(pd.Series(counts).to_string())

    parameters(out)
    qc_tables(res, out)
    position_tables(res, out)
    gs = array_weighted_group_summary(res)
    gs.to_csv(out / "group_summary.tsv", sep="\t", index=False, float_format="%.8g")
    pt = pairwise_tests(res)
    pt.to_csv(out / "group_pairwise_tests.tsv", sep="\t", index=False, float_format="%.8g")
    cl = cluster_aware(res, args.skip_r, out)
    cl.to_csv(out / "cluster_aware_sensitivity.tsv", sep="\t", index=False, float_format="%.8g")

    mt = []
    for _, r in gs[gs.tier.isin(["PRIMARY", "BACKBONE"]) & (gs.group != "ALL")].iterrows():
        mt.append(dict(family=f"{r.analysis}: prespecified group sign-flip tests", test=r.group, p=r.signflip_p_one_sided_mean_gt0,
                       q_bh=r.q_bh_across_group_tests, n_tests_in_family=int(((gs.analysis == r.analysis) & (gs.group != 'ALL')).sum()), method="BH"))
    for _, r in pt.iterrows():
        mt.append(dict(family=f"{r.analysis}: pairwise group contrasts", test=f"{r.group_a} vs {r.group_b}", p=r.p_two_sided,
                       q_bh=r.q_bh_across_pairwise, n_tests_in_family=int((pt.analysis == r.analysis).sum()), method="BH"))
    for pre, status, label in (("", "primary_status", "primary_300col_matched"), ("bb_", "bb_status", "backbone_distal_matched")):
        ok = res[status].eq("estimable")
        for scope, k in [("ALL", ok)] + [(g, ok & res.recombination_class.eq(g)) for g in GROUP_ORDER if res.recombination_class.eq(g).any()]:
            q = res.loc[k, f"{pre}q_bh_all_estimable"]
            mt.append(dict(family=f"{label}: per-array matched-window p (BH across all estimable arrays)", test=f"{scope}: arrays with q<0.05",
                           p=np.nan, q_bh=np.nan, n_tests_in_family=int(ok.sum()), method=f"BH; n_q_lt_0.05={int((q < 0.05).sum())} of {int(k.sum())}"))
    pd.DataFrame(mt).to_csv(out / "multiple_testing.tsv", sep="\t", index=False, float_format="%.8g")
    if missing:
        pd.Series(missing, name="array_id_without_COMPLETE").to_csv(out / "arrays_without_complete.tsv", index=False)
    print("aggregation written to", out)


if __name__ == "__main__":
    main()
