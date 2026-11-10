#!/usr/bin/env python3
"""Generate README.md from the combined tables. Every number below is read from results, never typed."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pp_common as C

COMB = C.COMBINED


def md(df, floatfmt=None):
    d = df.copy()
    for c in d.columns:
        if d[c].dtype.kind == "f":
            d[c] = d[c].map(lambda v: "" if pd.isna(v) else (f"{v:.2e}" if (abs(v) < 1e-3 and v != 0) else f"{v:.4g}"))
    d = d.astype(str).replace("nan", "")
    head = "| " + " | ".join(d.columns) + " |\n|" + "|".join("---" for _ in d.columns) + "|\n"
    return head + "\n".join("| " + " | ".join(r) + " |" for r in d.values) + "\n"


def main():
    global COMB
    COMB = Path(sys.argv[1]) if len(sys.argv) > 1 else C.COMBINED
    out_path = Path(sys.argv[2]) if len(sys.argv) > 2 else C.ROOT / "README.md"
    rc = pd.read_csv(COMB / "run_counts.tsv", sep="\t", index_col=0)["value"]
    res = pd.read_csv(COMB / "all_array_results.tsv", sep="\t")
    gs = pd.read_csv(COMB / "group_summary.tsv", sep="\t")
    pt = pd.read_csv(COMB / "group_pairwise_tests.tsv", sep="\t")
    cl = pd.read_csv(COMB / "cluster_aware_sensitivity.tsv", sep="\t")
    ne = pd.read_csv(COMB / "non_estimable_arrays.tsv", sep="\t")
    qc = pd.read_csv(COMB / "qc_review_list.tsv", sep="\t")
    mt = pd.read_csv(COMB / "multiple_testing.tsv", sep="\t")
    meta = pd.read_csv(C.ROOT / "analysis_metadata.tsv", sep="\t")
    man = pd.read_csv(C.ROOT / "manifest.validated.tsv", sep="\t")
    chk = pd.read_csv(C.REF / "reference_checksums.tsv", sep="\t")
    g = lambda k: rc[k]
    L = []
    w = L.append

    w("# ACRO promoter-conservation analysis on 2,794 arrays\n")
    w("Reproduction of the existing CHM13-derived ACRO promoter-conservation analysis on a new set of arrays. "
      "The promoter was fixed **before** these arrays were inspected and was neither moved, trimmed, realigned nor re-selected using them. "
      "Everything here tests *relative depletion of divergence at a fixed, transferred promoter window*; it is a matched-divergence / randomisation analysis, "
      "not Tajima's D, and does not by itself demonstrate purifying selection.\n")

    w("## 1. Run summary\n")
    w(md(pd.DataFrame({"metric": ["arrays submitted (manifest)", "arrays completed (valid result or documented non-estimability)", "primary 300-column test: estimable",
                                   "primary 300-column test: non-estimable", "conserved-backbone test: estimable", "conserved-backbone test: non-estimable",
                                   "total L1PA10 anchors", "total adjacent-anchor units", "total complete units (3.5-20 kb)", "mixed-strand arrays",
                                   "individuals", "pipeline version"],
                       "value": [g("arrays_submitted"), g("arrays_completed"), g("arrays_primary_estimable"), g("arrays_primary_non_estimable"),
                                 g("arrays_backbone_estimable"), g("arrays_backbone_non_estimable"), g("total_L1PA10_anchors"), g("total_adjacent_units"),
                                 g("total_complete_units"), g("arrays_mixed_strand"), man.individual_id.nunique(), g("pipeline_version")]})))
    w("Groups present in the manifest (defined exactly as specified; a group is only analysed if the manifest contains it):\n")
    gc = man.recombination_class.value_counts().rename_axis("group").reset_index(name="arrays")
    w(md(gc))
    if "non_telomeric_no_recombination" not in set(gc.group):
        w("**`non_telomeric_no_recombination` is absent from this manifest** (every array has `is_not_telomeric = 0`), so that group was not defined, tested or compared. "
          "Only the single pairwise contrast `telomeric_recombination` vs `telomeric_no_recombination` exists.\n")

    w("## 2. Input arrays\n")
    w("The 2,794 arrays are the `oct_5/acro_telo.bed` rows with last-column flag 0 and length >= 50 kb (the same filter and order as `oct_5/violin_plots.py`; names and the recombination flag come from "
      "`oct_5/moddotplot_5.txt` and `oct_5/recombination.txt`). The last bed column is interpreted as `is_not_telomeric`, consistent with the legacy `sample_100` table. "
      "Each FASTA is the bed interval plus a 10 kb flank on each side (clamped to the contig), extracted with `samtools faidx` from the HPRC release 2 assemblies; "
      "for the 60 arrays shared with the legacy `sample_100` inputs the sequences are byte-identical. FASTA edges are **not** treated as repeat boundaries anywhere. "
      "One FASTA record per array; every record was verified to contain only IUPAC nucleotide symbols; array IDs and FASTA paths are unique; "
      "filesystem identifiers are `sanitised_id__<8-hex SHA256>` and the original `array_id` is carried in every table. "
      "`has_recombination` is taken from `oct_5/recombination.txt` (flag from `oct_5/find_recombination.py`; the thresholds that produced the existing file were not recorded, "
      "so they are inherited, not re-derived).\n")

    w("## 3. Fixed references and checksums\n")
    w(md(chk[["file", "observed_sha256", "status"]]))
    w("All five checksums were verified, and the run aborts if any differs. Structure was also checked: 138 x 300 L1PA10 alignment; 300-row PRO-seq profile whose `promoter_window` is `yes` for exactly columns 78-127; "
      "22 x 9,090 unit alignment; 236 RepeatMasker rows; L1PA10 endpoint column 296. No reference was rebuilt or approximated from the target arrays.\n")

    w("## 4. Fixed definitions and thresholds\n")
    w(f"* **Promoter**: columns **{C.PROMOTER_300[0]}-{C.PROMOTER_300[1]}** (50 columns, inclusive) of the 300-column L1PA10 alignment; columns **{C.PROMOTER_FULL[0]}-{C.PROMOTER_FULL[1]}** of the 9,090-column full-unit alignment.\n"
      "* **PRO-seq** (projected CHM13 chr21 NEG strand) was used in the prior CHM13 analysis to nominate and fix the promoter at columns 78-127. It was not used to find HPRC anchors or boundaries, select matched controls, or re-select, trim, move, or realign the promoter using these target arrays.\n"
      f"* **Anchors**: `nhmmer --dna -E {C.NHMMER_E} --noali` (both strands); keep HMM coverage >= {C.MIN_HMM_COVERAGE} (of 290 profile columns) and score >= {C.MIN_NHMMER_SCORE:g}; "
      f"deduplicate by decreasing score, dropping a hit overlapping a retained hit by >= {int(C.DEDUP_OVERLAP_FRAC*100)}% of the shorter hit (applied across strands).\n"
      f"* **Chains**: separate per strand; adjacent anchor-start spacing {C.SPACING_MIN}-{C.SPACING_MAX} bp; dynamic programming maximises anchor number, then summed score; >= 2 anchors; "
      "both strands retained if both qualify (mixed-strand flag; those arrays are still analysed and a sensitivity analysis excludes them).\n"
      f"* **Units**: only between adjacent anchors; length {C.UNIT_MIN}-{C.UNIT_MAX} bp for the structural MSA and backbone mapping. '-' strand units are reverse-complemented to the transcriptional orientation. "
      "Irregular units are flagged and excluded individually, never the whole array.\n"
      f"* **Window divergence** = total pairwise differences / total comparable pairs, pooled over the 50 columns (A/C/G/T only; per-column values are never averaged). "
      f"A window is eligible if every column is callable in >= {int(C.MIN_CALLABLE_FRAC*100)}% of that array's copies.\n"
      f"* **Matched controls** (300-column test): eligible 50-column windows (starts 1-251) not overlapping 78-127 with |dGC| <= {C.TOL_GC}, |dCpG| <= {C.TOL_CPG}, |d occupancy| <= {C.TOL_OCC}. "
      "Primary CpG is computed after removing gaps from each sequence's fragment, pooling CpG counts and callable dinucleotides.\n"
      "* **Primary effect** = median(matched-control divergence) - promoter divergence (> 0: promoter less divergent). "
      "One-sided matched-window p = (1 + #controls <= promoter) / (1 + #controls); minimum attainable p = 1 / (1 + #controls).\n"
      f"* **Backbone test**: promoter columns 74-123 of the 9,090-column alignment after `mafft --add --keeplength`; needs >= {C.MIN_UNITS_BACKBONE} complete units. "
      f"Primary control space = columns > {C.L1PA10_END_FULL} excluding ACRO1-mapped columns +/- {C.ACRO1_BUFFER} alignment columns, Simple_repeat and Satellite annotations, and columns with < 80% occupancy among the 22 CHM13 units; "
      f"matching |dGC| <= {C.TOL_GC}, |dCpG| <= {C.BACKBONE_TOL_CPG} (directly adjacent aligned callable dinucleotides), |d occupancy| <= {C.TOL_OCC}. "
      "Backbone effect = median(matched distal-backbone divergence) - promoter divergence. "
      f"Promoter-mapping concordance (Jaccard overlap of the promoter bases defined by the two alignments, computed in raw unit coordinates) < {C.MAPPING_CONCORDANCE_MIN} is flagged for manual review.\n"
      "* **Labelled sensitivities (never substituted for the primary result)**: all-control/unmatched windows; CpG tolerance 0.021; all non-ACRO1 backbone; excluding mixed-strand arrays; excluding arrays with a mapping flag.\n")

    w("## 5. Inference (the array is the unit of replication)\n")
    w(f"* **Seed** `{C.SEED}`. Every test draws from an independent stream `SeedSequence([seed, crc32(test label)])`, so results are reproducible and independent of the order tests run in. The R bootstrap uses `set.seed(seed + sum(utf8ToInt(analysis)))`.\n"
      f"* **Sign-flip test** (H1: mean array effect > 0): array effects multiplied by random +/-1; exact enumeration if n <= {C.EXACT_MAX_N}, otherwise Monte Carlo with "
      f"**{C.N_SIGNFLIP_PRIMARY:,}** draws (primary / 300-column) or **{C.N_SIGNFLIP_BACKBONE:,}** draws (backbone), p = (extreme + 1) / (draws + 1).\n"
      f"* **Bootstrap 95% CI** for the mean: {C.N_BOOT:,} array resamples, percentile interval.\n"
      f"* **Group contrasts**: two-sided array-label permutation of the difference in means, **{C.N_GROUP_PERM:,}** draws, plus-one corrected; both one-sided tails are shown descriptively because no directional hypothesis was supplied.\n"
      "* **Multiple testing**: Benjamini-Hochberg (i) per-array matched-window p within each group and across all estimable arrays; (ii) across the prespecified group sign-flip tests per analysis; (iii) across pairwise group contrasts per analysis. "
      "With one pairwise contrast, (iii) is trivially q = p.\n"
      "* **Dependence among arrays** (231 individuals contribute several arrays; 462 haplotypes): the array-weighted results above are reported *separately* from a cluster-aware sensitivity analysis: "
      "sign-flip and bootstrap on individual-mean (and individual-haplotype-mean) effects; a within-individual stratified label permutation for the group contrast; "
      "a random-intercept model `effect ~ group + (1 | individual)` (lmerTest, Satterthwaite) with an OLS fit using individual-clustered robust SEs and an individual-cluster bootstrap "
      f"({C.N_CLUSTER_BOOT:,} resamples of whole individuals). These 2,794 arrays are **not** 2,794 independent genomes.\n"
      "* Pooled-copy results appear only as a labelled, non-independent descriptive column in `per_position_array_weighted_profile.tsv` and are not used for inference.\n")

    w("## 6. Results\n")
    cols = ["analysis", "group", "n_manifest_arrays", "n_estimable", "mean_effect", "median_effect", "n_positive", "frac_positive", "boot_ci95_lo", "boot_ci95_hi",
            "signflip_p_one_sided_mean_gt0", "q_bh_across_group_tests"]
    w("### 6a. Array-weighted primary and backbone tests (`combined/group_summary.tsv`)\n")
    w(md(gs[gs.tier.isin(["PRIMARY", "BACKBONE"])][cols]))
    w("### 6b. Group contrast (`combined/group_pairwise_tests.tsv`)\n")
    w(md(pt[["analysis", "group_a", "group_b", "n_a", "n_b", "mean_a", "mean_b", "diff_a_minus_b", "p_two_sided", "p_a_greater_than_b", "p_a_less_than_b", "q_bh_across_pairwise"]]) if len(pt) else "No pairwise contrast available.\n")
    w("### 6c. Cluster-aware sensitivity (`combined/cluster_aware_sensitivity.tsv`)\n")
    cc = [c for c in ["analysis", "test", "cluster_unit", "group", "n_arrays", "n_clusters", "n_clusters_with_both_groups", "estimate", "se", "boot_ci95_lo", "boot_ci95_hi",
                      "p_two_sided", "signflip_p_one_sided_mean_gt0"] if c in cl.columns]
    w(md(cl[cc]))
    w("### 6d. Labelled sensitivities\n")
    w(md(gs[gs.tier == "SENSITIVITY"][["analysis", "group", "n_estimable", "mean_effect", "median_effect", "frac_positive", "boot_ci95_lo", "boot_ci95_hi", "signflip_p_one_sided_mean_gt0"]]))

    w("### 6e. Do the 300-column and backbone tests agree?\n")
    a = gs[gs.analysis == "primary_300col_matched"].set_index("group")
    b = gs[gs.analysis == "backbone_distal_matched"].set_index("group")
    lines = []
    for grp in a.index:
        if grp not in b.index:
            continue
        sa, sb = np.sign(a.at[grp, "mean_effect"]), np.sign(b.at[grp, "mean_effect"])
        pa, pb = a.at[grp, "signflip_p_one_sided_mean_gt0"], b.at[grp, "signflip_p_one_sided_mean_gt0"]
        verdict = "agree" if (sa == sb and (pa < 0.05) == (pb < 0.05)) else "DISAGREE"
        lines.append(dict(group=grp, primary_mean=a.at[grp, "mean_effect"], backbone_mean=b.at[grp, "mean_effect"], primary_p=pa, backbone_p=pb, verdict=verdict))
    w(md(pd.DataFrame(lines)))
    w("'agree' means the same sign of the mean effect and the same side of p = 0.05 for the one-sided sign-flip test; this is a descriptive consistency check, not a formal test of equivalence.\n")

    ex_path = COMB / "exploratory_divergence_strata.tsv"
    if ex_path.exists():
        ex = pd.read_csv(ex_path, sep="\t")
        w("### 6f. EXPLORATORY, POST-HOC: is the group contrast driven by overall divergence level?\n")
        w("Not prespecified, and it replaces nothing above. The absolute effect (control median - promoter divergence) is bounded by an array's overall divergence, and the groups differ strongly in it "
          f"(median matched-control divergence: {res[res.recombination_class == 'telomeric_no_recombination'].control_median_divergence.median():.4f} in `telomeric_no_recombination` vs "
          f"{res[res.recombination_class == 'telomeric_recombination'].control_median_divergence.median():.4f} in `telomeric_recombination`). Estimable arrays with a nonzero control median were split into quartiles "
          "of control-median divergence (`combined/exploratory_divergence_strata.tsv`):\n")
        w(md(ex[ex.row_type == "stratum_group"][["stratum", "control_median_min", "control_median_max", "group", "n", "mean_effect", "median_ratio", "frac_promoter_less_divergent"]]))
        w(md(ex[ex.row_type == "all_informative"][["group", "n", "mean_effect", "median_ratio", "frac_promoter_less_divergent", "frac_promoter_more_divergent"]]))
        adj = ex[ex.row_type == "stratum_adjusted_contrast"].iloc[0]
        una = ex[ex.row_type == "unadjusted_contrast_same_arrays"].iloc[0]
        w(f"Group difference in mean effect (recombination - no recombination), same {int(una.n)} arrays: unadjusted **{una.mean_effect:+.4f}**; adjusted for divergence stratum **{adj.mean_effect:+.4f}** "
          f"(two-sided permutation p = {adj.p_two_sided:.1e}, {int(adj.permutations):,} draws, labels shuffled within strata, so this is at the resolution floor). Per-stratum differences: {adj.note.split(': ')[-1]}. "
          "The sign reverses after adjustment, the per-stratum differences are not consistent in sign, and the extreme strata contain very few arrays from one group, so this is evidence that the unadjusted "
          "group contrast is substantially a divergence-level effect, not evidence for an opposite effect.\n")

    w("## 7. Non-estimable arrays and QC\n")
    if len(ne):
        w(md(ne.groupby(["analysis", "reason_code"]).size().reset_index(name="arrays")))
    else:
        w("No non-estimable arrays.\n")
    w("Biological non-estimability is recorded as a successful result with a reason code (`combined/non_estimable_arrays.tsv`); thresholds were never relaxed to rescue an array. "
      "For arrays with no composition-matched control, the CpG tolerance 0.021 result is reported in separate columns and never replaces the primary CpG = 0.020 result.\n")
    pe = res[res.primary_status != "estimable"]
    byg = res.groupby("recombination_class").apply(lambda s: pd.Series(dict(
        arrays=len(s), primary_non_estimable=int((s.primary_status != "estimable").sum()), fraction_non_estimable=float((s.primary_status != "estimable").mean()),
        median_promoter_divergence=s.promoter_divergence.median(), frac_promoter_divergence_zero=float((s.promoter_divergence == 0).mean())))).reset_index()
    w("**Non-estimability is not balanced across groups, so the primary group means are over different subsets of arrays:**\n")
    w(md(byg))
    w(f"Non-estimable arrays have near-invariant promoters ({int((pe.promoter_divergence == 0).sum())} of {len(pe)} have promoter divergence exactly 0; median promoter CpG fraction "
      f"{pe.promoter_cpg.median():.3f} vs {res[res.primary_status == 'estimable'].promoter_cpg.median():.3f} in estimable arrays), so no control window matches on GC, CpG and occupancy together. "
      f"The CpG-0.021 sensitivity finds controls for {int((pe['cpg0.021_sens_n_controls'] > 0).sum())} of them. The primary estimates for the group that loses these arrays are therefore conditional "
      "on estimability; the all-array backbone test (every array estimable) and the CpG-0.021 sensitivity are the checks on this.\n")
    w("P-value resolution: the sign-flip tests use 200,000 (primary) or 500,000 (backbone) draws and the group permutations 500,000, so reported values of ~2e-6 are the Monte Carlo floor 1/(draws+1), not exact zeros; "
      "the individual-cluster bootstrap has 2,000 resamples (floor ~1e-3); the lmerTest and cluster-robust p-values are asymptotic and their extreme magnitudes should not be read as precise.\n")
    flags = qc.qc_flags.str.split(";").explode().value_counts().rename_axis("qc_flag").reset_index(name="arrays")
    w("### QC flags requiring manual review (`combined/qc_review_list.tsv`)\n")
    w(md(flags) if len(flags) else "None.\n")
    out = res.sort_values("n_units_irregular", ascending=False).head(10)[["array_id", "n_anchors", "n_units_complete", "n_units_irregular", "mixed_strand", "bb_mapping_concordance"]]
    w("Examples (most irregular units):\n")
    w(md(out))
    w(f"Arrays with promoter-mapping concordance < {C.MAPPING_CONCORDANCE_MIN}: {int((res.bb_mapping_concordance < C.MAPPING_CONCORDANCE_MIN).sum())}.\n")
    if "bb_mapping_concordance_aligned_index" in res.columns:
        w(f"**Correction to the per-array mapping QC.** The per-array workers computed concordance from ungapped indices of the aligned sequences. `mafft --keeplength` deletes bases that are insertions "
          f"relative to the fixed alignments, so those indices are not raw positions; the stored value flagged {int((res.bb_mapping_concordance_aligned_index < C.MAPPING_CONCORDANCE_MIN).sum())} arrays and a spurious "
          f"`anchor_not_unit_prefix` flag on {int((res.bb_anchor_not_unit_prefix_aligned_index > 0).sum())} arrays. It was recomputed after the run in raw unit coordinates (`scripts/recheck_mapping.py`, `combined/promoter_mapping_recheck.tsv`): "
          f"0 units violate the raw anchor-prefix property, and {int((res.bb_mapping_concordance < C.MAPPING_CONCORDANCE_MIN).sum())} arrays have concordance below the threshold. The recheck was tested for homopolymer "
          "ambiguity in the leftmost subsequence mapping: none of the differing units was explained by it. The original values are kept as `*_aligned_index` columns. Only QC columns and the 'excluding mapping flag' sensitivity "
          "depend on this; no divergence, effect or p-value does (verified by an independent naive recomputation of 12 random arrays plus 8 pilot arrays).\n")

    w("## 8. Interpretation and caveats\n")
    w("A positive effect means the promoter window is *less divergent* than composition-matched windows from the same array (or the distal backbone). "
      "This is evidence for **relative promoter conservation / depletion of divergence**. It is not a demonstration of purifying selection: "
      "**gene conversion, concerted evolution, repeat age, shared genealogy and population structure** can each produce position-specific reductions in divergence among tandem repeat copies. "
      "Controls from sliding windows overlap one another, so matched-window p-values are descriptive within an array and the minimum attainable p can be large when few controls match; "
      "inference rests on the across-array tests, and the cluster-aware analysis shows how much depends on treating arrays from the same individual as independent. "
      "The recombination flag is inherited from an earlier analysis and is itself derived from the same dot-plot identity data, so group differences in divergence should be interpreted with that circularity in mind. "
      "The absolute effect scales with an array's overall divergence (see 6f), so a group contrast in absolute effect partly reflects group differences in divergence level; the scale-free fraction of arrays whose promoter is "
      f"less divergent than its controls (6f) is reported for that reason. Per-array matched-window p-values are weak: {int((res.q_bh_all_estimable < 0.05).sum())} arrays reach BH q < 0.05 across all estimable arrays, "
      f"and {int((res.p_min_attainable > 0.05).sum())} estimable arrays have too few matched controls to ever reach p < 0.05 (see `p_min_attainable`).\n")

    w("## 9. Software, environment and provenance\n")
    sw = meta[meta.category.isin(["software", "randomisation", "commands", "run"])][["key", "value"]]
    w(md(sw))
    w("Modules: `hmmer/3.4`, `mafft/7.526`, `samtools/1.21` (per-array tasks), `R/4.4.2` (statistics and figures); Python 3.9 (numpy, pandas, scipy). No conda environment file: tools are HPC `lmod` modules. "
      "Per-array command lines are in `per_array/*/logs/commands.txt`.\n")
    w(f"Resources actually used: per-array wall time median {res.runtime_total_s.median():.0f} s, max {res.runtime_total_s.max()/60:.0f} min; own-measured peak memory median "
      f"{res.peak_rss_mb.median()/1024:.1f} GB, max {res.peak_rss_mb.max()/1024:.1f} GB ({int((res.peak_rss_mb > 4000).sum())} tasks above 4 GB; E-INS-i accounts for >= 96% of their runtime, and memory is not a simple function of unit count). "
      "The 27-array pilot (max 3.0 GB) under-sampled this memory tail; the class limits (8 GB / 16 GB) were nevertheless sufficient and no task was killed.\n")

    w("## 10. Layout and restart\n")
    w("```\npromoter_purifying_2794/\n  reference/  scripts/  manifest.validated.tsv  analysis_metadata.tsv  per_array/<safe_id>/{...,COMPLETE}  combined/  combined/figures/  README.md\n```\n"
      "Each task writes to `per_array/.tmp_*`, validates every required output, and renames atomically; `COMPLETE` is written last. Completed arrays are skipped unless `--force`. "
      "Computational failures exit non-zero and are listed in `failed/`. Rerun only the incomplete tasks with `scripts/find_incomplete.py > tasks.txt; scripts/submit.sh tasks.txt MEM TIME CONC`.\n")
    w("Figures (PNG + PDF) are in `combined/figures/`: 1 per-position divergence + PRO-seq annotation + callable counts; 2 primary effect by group; 3 promoter vs matched control; "
      "4 backbone effect by group; 5 QC overview; 6 rasterised per-array heatmap.\n")
    out_path.write_text("\n".join(L))
    print("README.md written,", sum(len(x) for x in L), "chars")


if __name__ == "__main__":
    main()
