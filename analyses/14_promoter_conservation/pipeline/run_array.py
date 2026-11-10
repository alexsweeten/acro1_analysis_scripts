#!/usr/bin/env python3
"""Per-array worker: phases 1-4 for one ACRO array. Restartable and atomic.

  run_array.py --task-index N [--force] [--cpus K]

* Writes everything into per_array/.tmp_<safe_id>_<pid>/ and renames it to per_array/<safe_id>/ only after
  all required outputs validate; the COMPLETE sentinel is written last.
* Biological non-estimability is a normal result (exit 0, reason code in result.tsv).
* Computational failures (nhmmer/mafft errors, wrong alignment width, ...) exit non-zero and leave a
  FAILED.json under failed/<safe_id>/.
* Nothing here is tuned on the target arrays: thresholds and the promoter come from pp_common.py.
"""
import argparse
import hashlib
import json
import os
import resource
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pp_common as C
import pp_core as K

HERE = Path(__file__).resolve().parent
NAN = float("nan")


def pipeline_version():
    h = hashlib.sha256()
    for f in ("pp_common.py", "pp_core.py", "run_array.py"):
        h.update((HERE / f).read_bytes())
    return h.hexdigest()[:12]


class CompFailure(Exception):
    def __init__(self, code, msg):
        super().__init__(f"{code}: {msg}")
        self.code, self.msg = code, msg


REQUIRED = ["anchor_hits.tsv", "l1pa10_anchors.tsv", "l1pa10_anchors.fa", "l1pa10_anchors.aligned.fa",
            "units_all_adjacent.fa", "units_transcription.fa", "units_qc.tsv", "units_transcription.einsi.fa",
            "per_position_profile.tsv", "matched_window_test.tsv", "result.tsv",
            "backbone/units_reference_mapped.fa", "backbone/matched_backbone_windows.tsv", "backbone/result.tsv",
            "logs/anchor_discovery.log"]


class Ctx:
    def __init__(self, row, tmp, cpus):
        self.row, self.tmp, self.cpus = row, Path(tmp), cpus
        self.t = {}
        self.cmds = []
        self.logf = open(self.tmp / "logs" / "anchor_discovery.log", "w")
        self.scratch = self.tmp / "scratch"
        self.scratch.mkdir(exist_ok=True)

    def log(self, msg):
        self.logf.write(msg + "\n")
        self.logf.flush()

    def run(self, key, cmd, stdout=None, stderr=None, code="mafft_failure", env=None):
        self.cmds.append(" ".join(cmd) + (f" > {Path(stdout).name}" if stdout else ""))
        t0 = time.time()
        so = open(stdout, "w") if stdout else subprocess.DEVNULL
        se = open(stderr, "w") if stderr else subprocess.DEVNULL
        e = dict(os.environ, TMPDIR=str(self.scratch), **(env or {}))
        try:
            p = subprocess.run(cmd, stdout=so, stderr=se, env=e)
        finally:
            for f in (so, se):
                if hasattr(f, "close"):
                    f.close()
        self.t[key] = self.t.get(key, 0) + time.time() - t0
        if p.returncode != 0:
            raise CompFailure(code, f"{cmd[0]} exited {p.returncode}")


# ---------------------------------------------------------------- phase 1
def parse_tblout(path, hmm_len):
    hits = []
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        f = line.split(None, 15)
        a, b = int(f[6]), int(f[7])
        hits.append(dict(hit_id=len(hits) + 1, strand=f[11], start=min(a, b), end=max(a, b),
                         hmm_from=int(f[4]), hmm_to=int(f[5]), evalue=float(f[12]), score=float(f[13]),
                         bias=float(f[14])))
    return K.filter_hits(hits, hmm_len)


def phase1(ctx, seq):
    tmp = ctx.tmp
    hmm_len = int(next(l.split()[1] for l in open(C.HMM) if l.startswith("LENG")))
    ctx.run("t_nhmmer_s", ["nhmmer", "--dna", "--tblout", str(tmp / "logs/nhmmer.tbl"), "--noali", "-E", C.NHMMER_E,
                           "--cpu", str(ctx.cpus), str(C.HMM), ctx.row["fasta_path"]],
            stdout=tmp / "logs/nhmmer.out", stderr=tmp / "logs/nhmmer.err", code="nhmmer_failure")
    hits = parse_tblout(tmp / "logs/nhmmer.tbl", hmm_len)
    retained = K.dedup_hits(hits)
    chains = K.build_chains(retained)
    ctx.log(f"hits total={len(hits)} pass_filter={sum(h['filter_status']=='pass' for h in hits)} "
            f"retained_after_dedup={len(retained)} chains={ {s: len(c) for s, c in chains.items()} }")
    for h in hits:
        ctx.log(f"hit {h['hit_id']} {h['strand']} {h['start']}-{h['end']} cov={h['hmm_coverage']:.3f} "
                f"score={h['score']} {h['filter_status']} {h['dedup_status']}"
                + (f" dup_of={h['duplicate_of']} ov={h['overlap_frac_shorter']:.2f}" if h["duplicate_of"] else ""))
    # transcription order: '+' ascending, '-' descending genomic start
    anchors, k = [], 0
    for s in "+-":
        if s not in chains:
            continue
        ch = chains[s] if s == "+" else chains[s][::-1]
        for pos, h in enumerate(ch):
            k += 1
            h["anchor_id"] = f"a{k:03d}{'p' if s == '+' else 'm'}"
            h["chain_order"] = pos + 1
            h["chain_strand"] = s
            anchors.append(h)
    sel = {h["hit_id"] for h in anchors}
    for h in hits:
        h["selected"] = h["hit_id"] in sel
    fs = ctx.row["flanked_start"]
    rows = []
    for h in hits:
        rows.append(dict(array_id=ctx.row["array_id"], hit_id=h["hit_id"], strand=h["strand"], local_start=h["start"],
                         local_end=h["end"], contig_start=fs - 1 + h["start"], contig_end=fs - 1 + h["end"],
                         length=h["end"] - h["start"] + 1, hmm_from=h["hmm_from"], hmm_to=h["hmm_to"],
                         hmm_coverage=round(h["hmm_coverage"], 4), score=h["score"], evalue=h["evalue"], bias=h["bias"],
                         filter_status=h["filter_status"], dedup_status=h["dedup_status"], duplicate_of=h["duplicate_of"],
                         overlap_frac_shorter=round(h["overlap_frac_shorter"], 4),
                         selected_in_chain=h["selected"], anchor_id=h.get("anchor_id", "")))
    cols = ["array_id", "hit_id", "strand", "local_start", "local_end", "contig_start", "contig_end", "length", "hmm_from",
            "hmm_to", "hmm_coverage", "score", "evalue", "bias", "filter_status", "dedup_status", "duplicate_of",
            "overlap_frac_shorter", "selected_in_chain", "anchor_id"]
    pd.DataFrame(rows, columns=cols).to_csv(tmp / "anchor_hits.tsv", sep="\t", index=False)

    # anchors table, with inter-anchor start spacing along the chain (genomic, positive)
    arows = []
    for s in "+-":
        sub = [a for a in anchors if a["chain_strand"] == s]
        for i, a in enumerate(sub):
            nxt = sub[i + 1] if i + 1 < len(sub) else None
            prv = sub[i - 1] if i > 0 else None
            sp = lambda x, y: abs(x["start"] - y["start"]) if x is not None and y is not None else ""
            arows.append(dict(array_id=ctx.row["array_id"], anchor_id=a["anchor_id"], strand=s, chain_order=a["chain_order"],
                              hit_id=a["hit_id"], local_start=a["start"], local_end=a["end"],
                              contig_start=fs - 1 + a["start"], contig_end=fs - 1 + a["end"], length=a["end"] - a["start"] + 1,
                              hmm_coverage=round(a["hmm_coverage"], 4), score=a["score"], evalue=a["evalue"],
                              spacing_from_previous=sp(a, prv), spacing_to_next=sp(a, nxt)))
    acols = ["array_id", "anchor_id", "strand", "chain_order", "hit_id", "local_start", "local_end", "contig_start",
             "contig_end", "length", "hmm_coverage", "score", "evalue", "spacing_from_previous", "spacing_to_next"]
    pd.DataFrame(arows, columns=acols).to_csv(tmp / "l1pa10_anchors.tsv", sep="\t", index=False)
    C.write_fasta(tmp / "l1pa10_anchors.fa", [(a["anchor_id"], K.oriented_anchor(seq, a)) for a in anchors])

    # units: adjacent intervals only, in transcription orientation
    units, u = [], 0
    for s in "+-":
        if s in chains:
            for un in K.build_units(seq, chains[s], s):
                u += 1
                un["unit_id"] = f"u{u:03d}{'p' if s == '+' else 'm'}"
                un["length"] = len(un["seq"])
                un["in_range"] = C.UNIT_MIN <= un["length"] <= C.UNIT_MAX
                units.append(un)
    urows = []
    for un in units:
        nn = un["seq"].upper().count("N") / max(1, un["length"])
        reason = "" if un["in_range"] else ("unit_length_lt_3500" if un["length"] < C.UNIT_MIN else "unit_length_gt_20000")
        urows.append(dict(array_id=ctx.row["array_id"], unit_id=un["unit_id"], strand=un["strand"], local_start=un["local_start"],
                          local_end=un["local_end"], contig_start=fs - 1 + un["local_start"], contig_end=fs - 1 + un["local_end"],
                          length=un["length"], start_anchor_id=un["start_anchor"]["anchor_id"],
                          next_anchor_id=un["next_anchor"]["anchor_id"], n_fraction=round(nn, 5),
                          irregular=not un["in_range"], included_in_structural_msa=un["in_range"], exclusion_reason=reason))
    ucols = ["array_id", "unit_id", "strand", "local_start", "local_end", "contig_start", "contig_end", "length",
             "start_anchor_id", "next_anchor_id", "n_fraction", "irregular", "included_in_structural_msa", "exclusion_reason"]
    pd.DataFrame(urows, columns=ucols).to_csv(tmp / "units_qc.tsv", sep="\t", index=False)
    C.write_fasta(tmp / "units_all_adjacent.fa", [(un["unit_id"], un["seq"]) for un in units])
    C.write_fasta(tmp / "units_transcription.fa", [(un["unit_id"], un["seq"]) for un in units if un["in_range"]])
    return hits, anchors, chains, units


# ---------------------------------------------------------------- phase 2
def mafft_extract(path, wanted_ids, width, fail_code):
    recs = {n.split()[0]: s for n, s in C.read_fasta(path)}
    missing = [i for i in wanted_ids if i not in recs]
    if missing:
        raise CompFailure("mafft_failure", f"{len(missing)} target records missing from MAFFT output")
    out = [(i, recs[i].upper()) for i in wanted_ids]
    bad = [i for i, s in out if len(s) != width]
    if bad:
        raise CompFailure(fail_code, f"{len(bad)} records not {width} columns (e.g. {bad[0]}={len(dict(out)[bad[0]])})")
    return out, recs


def check_reference_unchanged(recs, ref_path):
    for n, s in C.read_fasta(ref_path):
        if n.split()[0] not in recs or recs[n.split()[0]].upper() != s.upper():
            raise CompFailure("mafft_failure", f"fixed reference alignment altered at {n}")


def primary_test(ctx, anchors):
    tmp = ctx.tmp
    ids = [a["anchor_id"] for a in anchors]
    ctx.run("t_mafft_add_s", ["mafft", "--quiet", "--thread", str(ctx.cpus), "--addfragments", str(tmp / "l1pa10_anchors.fa"),
                              "--keeplength", str(C.REF_ALN300)], stdout=tmp / "scratch/anchors_plus_ref.fa",
            stderr=tmp / "logs/mafft_addfragments.err")
    aligned, recs = mafft_extract(tmp / "scratch/anchors_plus_ref.fa", ids, C.ALN_COLS, "anchor_alignment_not_300_columns")
    check_reference_unchanged(recs, C.REF_ALN300)
    C.write_fasta(tmp / "l1pa10_anchors.aligned.fa", aligned)
    X = K.encode([s for _, s in aligned])
    cs, tab, ps, ref = K.primary_windows(X)
    pd.DataFrame(tab).to_csv(tmp / "matched_window_test.tsv", sep="\t", index=False, float_format="%.8g")

    # Per-position profile; prior CHM13 PRO-seq nominated the fixed promoter
    # but is not used in per-array anchor or divergence calculations.
    pro = pd.read_csv(C.PROSEQ, sep="\t")
    n = cs["n"]
    with np.errstate(divide="ignore", invalid="ignore"):
        coldiv = np.where(cs["pairs"] > 0, cs["diffs"] / cs["pairs"], np.nan)
    prof = pd.DataFrame(dict(
        alignment_column_1based=np.arange(1, C.ALN_COLS + 1), anchor_count=cs["ncall"], occupancy=cs["ncall"] / n,
        pairwise_divergence=coldiv, rolling_11col_divergence=K.rolling_pooled(cs["diffs"], cs["pairs"], C.ROLL_DIV // 2),
        projected_chr21_NEG_PROseq=pro["PROseq_mean"].to_numpy(),
        rolling_7col_PROseq=K.rolling_mean(pro["PROseq_mean"].to_numpy(), C.ROLL_PRO // 2),
        promoter_window=np.where((np.arange(1, C.ALN_COLS + 1) >= C.PROMOTER_300[0]) & (np.arange(1, C.ALN_COLS + 1) <= C.PROMOTER_300[1]), "yes", "no"),
        pairwise_differences=cs["diffs"], comparable_pairs=cs["pairs"]))
    prof.to_csv(tmp / "per_position_profile.tsv", sep="\t", index=False, float_format="%.8g")

    pdiv = tab["divergence"][ps]
    res = dict(n_aligned_anchors=n, promoter_divergence=pdiv, promoter_gc=tab["gc"][ps], promoter_cpg=tab["cpg"][ps],
               promoter_occupancy=tab["occupancy"][ps], promoter_segregating_sites=int(tab["segregating_sites"][ps]),
               promoter_mean_entropy_bits=tab["mean_entropy_bits"][ps], promoter_eligible=bool(tab["eligible"][ps]),
               n_eligible_windows=int(tab["eligible"].sum()), n_control_candidates=int(tab["control_candidate"].sum()))
    if n < 2:
        res.update(primary_status="non_estimable", primary_reason="fewer_than_two_aligned_anchors")
        return res
    if not res["promoter_eligible"]:
        res.update(primary_status="non_estimable", primary_reason="promoter_below_80pct_occupancy")
    ctrl_all = tab["divergence"][tab["control_candidate"]]
    ctrl_m = tab["divergence"][tab["matched_primary"]]
    ctrl_g = tab["divergence"][tab["matched_cpg_grid_sensitivity"]]
    if res["promoter_eligible"]:
        pm = K.summarise_vs_controls(pdiv, ctrl_m)
        res.update(n_matched_controls=pm["n_controls"], control_median_divergence=pm["control_median"], effect=pm["effect"],
                   ratio=pm["ratio"], p_matched=pm["p"], p_min_attainable=pm["p_min"])
        if pm["n_controls"] == 0:
            res.update(primary_status="non_estimable", primary_reason="no_composition_matched_controls")
        else:
            res.update(primary_status="estimable", primary_reason="")
        for tag, ctl in (("allctrl_unmatched_sens", ctrl_all), ("cpg0.021_sens", ctrl_g)):
            sm = K.summarise_vs_controls(pdiv, ctl)
            res.update({f"{tag}_n_controls": sm["n_controls"], f"{tag}_control_median": sm["control_median"],
                        f"{tag}_effect": sm["effect"], f"{tag}_p": sm["p"]})
    else:
        res.update(n_matched_controls=0)
    return res, aligned


# ---------------------------------------------------------------- phases 3 & 4
def structural_alignment(ctx, n_units, units):
    tmp = ctx.tmp
    out = tmp / "units_transcription.einsi.fa"
    if n_units == 0:
        out.write_text("")
        return "no_complete_units"
    if n_units == 1:
        u = next(x for x in units if x["in_range"])
        C.write_fasta(out, [(f"{u['unit_id']} single_unit_no_alignment_needed", u["seq"])])
        return "single_unit_no_alignment_needed"
    ctx.run("t_mafft_einsi_s", ["mafft", "--genafpair", "--maxiterate", "1000", "--ep", "0", "--thread", str(ctx.cpus),
                                str(tmp / "units_transcription.fa")], stdout=out, stderr=tmp / "logs/mafft_einsi.err")
    if out.stat().st_size == 0:
        raise CompFailure("mafft_failure", "E-INS-i produced empty output")
    return "einsi_aligned"


def promoter_concordance(units, anchors_aligned, mapped):
    """Do the bases in full-unit columns 74-123 coincide with those in 300-column anchor columns 78-127?"""
    a0, a1 = C.PROMOTER_300[0] - 1, C.PROMOTER_300[1]
    b0, b1 = C.PROMOTER_FULL[0] - 1, C.PROMOTER_FULL[1]
    amap, rows = dict(anchors_aligned), []
    inter = union = 0
    prefix_fail = 0
    for u in units:
        if not u["in_range"]:
            continue
        A, B = amap[u["start_anchor"]["anchor_id"]], dict(mapped)[u["unit_id"]]
        a_un = A.replace("-", "")
        ok = u["seq"].upper().startswith(a_un)
        prefix_fail += (not ok)

        def positions(aln, lo, hi):
            idx = np.cumsum([c != "-" for c in aln]) - 1
            return {int(idx[i]) for i in range(lo, hi) if aln[i] != "-"}

        sa, sb = positions(A, a0, a1), positions(B, b0, b1)
        i, un = len(sa & sb), len(sa | sb)
        inter, union = inter + i, union + un
        rows.append(dict(unit_id=u["unit_id"], anchor_id=u["start_anchor"]["anchor_id"], anchor_is_unit_prefix=ok,
                         anchor_window_bases=len(sa), unit_window_bases=len(sb), shared_bases=i,
                         anchor_window_gaps=(a1 - a0) - len(sa), unit_window_gaps=(b1 - b0) - len(sb)))
    conc = inter / union if union else NAN
    return conc, prefix_fail, pd.DataFrame(rows)


def backbone_test(ctx, units, anchors_aligned, res):
    tmp = ctx.tmp
    bdir = tmp / "backbone"
    ok_units = [u for u in units if u["in_range"]]
    out = dict(bb_n_units_mapped=0)
    if len(ok_units) < C.MIN_UNITS_BACKBONE:
        out.update(bb_status="non_estimable",
                   bb_reason="no_complete_units" if not ok_units else "fewer_than_three_complete_units_for_backbone")
        return out
    ids = [u["unit_id"] for u in ok_units]
    ctx.run("t_mafft_backbone_s", ["mafft", "--quiet", "--thread", "1", "--retree", "1", "--maxiterate", "0", "--add",
                                   str(tmp / "units_transcription.fa"), "--keeplength", str(C.REF_ALN9090)],
            stdout=tmp / "scratch/units_plus_ref.fa", stderr=tmp / "logs/mafft_backbone.err")
    mapped, recs = mafft_extract(tmp / "scratch/units_plus_ref.fa", ids, C.FULL_COLS, "unit_mapping_not_9090_columns")
    check_reference_unchanged(recs, C.REF_ALN9090)
    C.write_fasta(bdir / "units_reference_mapped.fa", mapped)
    conc, prefix_fail, qc = promoter_concordance(units, anchors_aligned, mapped)
    qc.to_csv(bdir / "promoter_mapping_qc.tsv", sep="\t", index=False)
    mask = pd.read_csv(C.BACKBONE_MASK, sep="\t")
    U = K.encode([s for _, s in mapped])
    cs, tab, ps, ref = K.backbone_windows(U, mask)
    pd.DataFrame(tab).to_csv(bdir / "matched_backbone_windows.tsv", sep="\t", index=False, float_format="%.8g")
    pdiv = tab["divergence"][ps]
    out.update(bb_n_units_mapped=len(mapped), bb_promoter_divergence=pdiv, bb_promoter_gc=tab["gc"][ps],
               bb_promoter_cpg_adjacent=tab["cpg"][ps], bb_promoter_occupancy=tab["occupancy"][ps],
               bb_promoter_eligible=bool(tab["eligible"][ps]), bb_mapping_concordance=conc,
               bb_anchor_not_unit_prefix=prefix_fail,
               bb_promoter_anchor_window_gaps=int(qc.anchor_window_gaps.sum()), bb_promoter_unit_window_gaps=int(qc.unit_window_gaps.sum()),
               bb_n_primary_distal_candidates=int(tab["primary_distal_candidate"].sum()))
    flags = []
    if not np.isfinite(conc) or conc < C.MAPPING_CONCORDANCE_MIN:
        flags.append("promoter_mapping_concordance_below_0.99")
    if prefix_fail:
        flags.append("anchor_not_unit_prefix")
    out["bb_review_flags"] = ";".join(flags)
    if not out["bb_promoter_eligible"]:
        out.update(bb_status="non_estimable", bb_reason="promoter_below_80pct_occupancy")
        return out
    pm = K.summarise_vs_controls(pdiv, tab["divergence"][tab["matched_primary_distal"]])
    out.update(bb_n_matched_controls=pm["n_controls"], bb_control_median_divergence=pm["control_median"], bb_effect=pm["effect"],
               bb_ratio=pm["ratio"], bb_p_matched=pm["p"], bb_p_min_attainable=pm["p_min"])
    sm = K.summarise_vs_controls(pdiv, tab["divergence"][tab["matched_all_non_acro1"]])
    out.update(bb_allnonacro1_sens_n_controls=sm["n_controls"], bb_allnonacro1_sens_control_median=sm["control_median"],
               bb_allnonacro1_sens_effect=sm["effect"], bb_allnonacro1_sens_p=sm["p"])
    if pm["n_controls"] == 0:
        out.update(bb_status="non_estimable", bb_reason="no_composition_matched_controls")
    else:
        out.update(bb_status="estimable", bb_reason="")
    return out


# ---------------------------------------------------------------- driver
RESULT_COLS = [
    "array_id", "safe_id", "individual_id", "haplotype", "population", "superpopulation", "species", "is_not_telomeric",
    "has_recombination", "recombination_class", "chromosome_or_contig", "original_start", "original_end", "fasta_length",
    "pipeline_version", "primary_status", "primary_reason", "n_hits_total", "n_hits_pass_filter", "n_hits_rejected_duplicate",
    "n_hits_retained", "n_anchors", "n_anchors_plus", "n_anchors_minus", "mixed_strand", "n_units_adjacent", "n_units_complete",
    "n_units_irregular", "structural_alignment", "promoter_divergence", "control_median_divergence", "effect", "ratio",
    "n_matched_controls", "p_matched", "p_min_attainable", "promoter_occupancy", "promoter_gc", "promoter_cpg",
    "promoter_segregating_sites", "promoter_mean_entropy_bits", "promoter_eligible", "n_aligned_anchors", "n_eligible_windows",
    "n_control_candidates", "allctrl_unmatched_sens_n_controls", "allctrl_unmatched_sens_control_median",
    "allctrl_unmatched_sens_effect", "allctrl_unmatched_sens_p", "cpg0.021_sens_n_controls", "cpg0.021_sens_control_median",
    "cpg0.021_sens_effect", "cpg0.021_sens_p", "bb_status", "bb_reason", "bb_n_units_mapped", "bb_promoter_divergence",
    "bb_control_median_divergence", "bb_effect", "bb_ratio", "bb_n_matched_controls", "bb_p_matched", "bb_p_min_attainable",
    "bb_promoter_occupancy", "bb_promoter_gc", "bb_promoter_cpg_adjacent", "bb_promoter_eligible", "bb_mapping_concordance",
    "bb_anchor_not_unit_prefix", "bb_promoter_anchor_window_gaps", "bb_promoter_unit_window_gaps",
    "bb_n_primary_distal_candidates", "bb_review_flags", "bb_allnonacro1_sens_n_controls",
    "bb_allnonacro1_sens_control_median", "bb_allnonacro1_sens_effect", "bb_allnonacro1_sens_p", "qc_flags",
    "runtime_total_s", "t_nhmmer_s", "t_mafft_add_s", "t_mafft_einsi_s", "t_mafft_backbone_s", "peak_rss_mb", "cpus", "host"]


def touch_placeholders(tmp):
    for f, hdr in (("matched_window_test.tsv", "start_1based\tend_1based\twindow_class\n"),
                   ("per_position_profile.tsv", "alignment_column_1based\tanchor_count\toccupancy\n"),
                   ("l1pa10_anchors.aligned.fa", ""), ("units_transcription.einsi.fa", ""),
                   ("backbone/units_reference_mapped.fa", ""),
                   ("backbone/matched_backbone_windows.tsv", "start_1based\tend_1based\twindow_class\n")):
        if not (tmp / f).exists():
            (tmp / f).write_text(hdr)


def process(row, tmp, cpus):
    ctx = Ctx(row, tmp, cpus)
    (tmp / "backbone").mkdir(exist_ok=True)
    t0 = time.time()
    recs = C.read_fasta(row["fasta_path"])
    if len(recs) != 1 or not C.IUPAC_RE.match(recs[0][1]):
        raise CompFailure("malformed_input_fasta", "expected exactly one IUPAC-only record")
    seq = recs[0][1].upper()
    hits, anchors, chains, units = phase1(ctx, seq)
    n_pass = sum(h["filter_status"] == "pass" for h in hits)
    res = dict(n_hits_total=len(hits), n_hits_pass_filter=n_pass,
               n_hits_rejected_duplicate=sum(h["dedup_status"] == "rejected_duplicate" for h in hits),
               n_hits_retained=sum(h["dedup_status"] == "retained" for h in hits), n_anchors=len(anchors),
               n_anchors_plus=len(chains.get("+", [])), n_anchors_minus=len(chains.get("-", [])),
               mixed_strand=len(chains) == 2, n_units_adjacent=len(units), n_units_complete=sum(u["in_range"] for u in units),
               n_units_irregular=sum(not u["in_range"] for u in units))
    flags = []
    if res["mixed_strand"]:
        flags.append("mixed_strand_chains")
    if res["n_units_irregular"]:
        flags.append("irregular_unit_spacing")
    bb = {}
    if not anchors:
        reason = "no_qualifying_L1PA10_hits" if n_pass == 0 else "no_tandem_chain_with_two_anchors"
        res.update(primary_status="non_estimable", primary_reason=reason, structural_alignment="not_run")
        bb = dict(bb_status="non_estimable", bb_reason=reason, bb_n_units_mapped=0)
        ctx.log(f"non-estimable: {reason}")
        touch_placeholders(tmp)
        structural_alignment(ctx, 0, units)
    else:
        out = primary_test(ctx, anchors)
        res.update(out[0] if isinstance(out, tuple) else out)
        aligned = out[1] if isinstance(out, tuple) else []
        if res["primary_status"] == "non_estimable":
            ctx.log(f"primary non-estimable: {res['primary_reason']}")
        touch_placeholders(tmp)
        bb = backbone_test(ctx, units, aligned, res) if aligned else dict(
            bb_status="non_estimable", bb_reason="fewer_than_two_aligned_anchors", bb_n_units_mapped=0)
        res["structural_alignment"] = structural_alignment(ctx, res["n_units_complete"], units)
    res.update(bb)
    if res.get("bb_review_flags"):
        flags.extend(res["bb_review_flags"].split(";"))
    for k in ("primary_reason", "bb_reason"):
        if res.get(k) == "no_complete_units":
            flags.append("no_complete_units")
    res["qc_flags"] = ";".join(dict.fromkeys(flags))
    ru = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    res.update(runtime_total_s=round(time.time() - t0, 2), peak_rss_mb=round(max(ru, resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) / 1024, 1),
               cpus=cpus, host=os.uname().nodename, pipeline_version=pipeline_version(), **{k: round(v, 2) for k, v in ctx.t.items()})
    for k in ("array_id", "safe_id", "individual_id", "haplotype", "population", "superpopulation", "species", "is_not_telomeric",
              "has_recombination", "recombination_class", "chromosome_or_contig", "original_start", "original_end", "fasta_length"):
        res[k] = row[k]
    df = pd.DataFrame([{c: res.get(c, NAN) for c in RESULT_COLS}])
    df.to_csv(tmp / "result.tsv", sep="\t", index=False, float_format="%.10g")
    bcols = [c for c in RESULT_COLS if c.startswith("bb_")]
    df[["array_id"] + bcols].to_csv(tmp / "backbone/result.tsv", sep="\t", index=False, float_format="%.10g")
    (tmp / "logs/commands.txt").write_text("\n".join(ctx.cmds) + "\n")
    ctx.logf.close()
    shutil.rmtree(ctx.scratch, ignore_errors=True)
    return res


def validate(tmp, res):
    miss = [f for f in REQUIRED if not (tmp / f).exists()]
    if miss:
        raise CompFailure("incomplete_outputs", f"missing {miss}")
    r = pd.read_csv(tmp / "result.tsv", sep="\t")
    if len(r) != 1 or r.at[0, "primary_status"] not in ("estimable", "non_estimable"):
        raise CompFailure("incomplete_outputs", "result.tsv malformed")
    if r.at[0, "primary_status"] == "estimable":
        prof = pd.read_csv(tmp / "per_position_profile.tsv", sep="\t")
        if len(prof) != C.ALN_COLS or {len(s) for _, s in C.read_fasta(tmp / "l1pa10_anchors.aligned.fa")} != {C.ALN_COLS}:
            raise CompFailure("incomplete_outputs", "primary outputs not 300 columns")
        if prof.promoter_window.eq("yes").sum() != C.WIN:
            raise CompFailure("incomplete_outputs", "promoter window is not 50 columns")
    if r.at[0, "bb_status"] == "estimable":
        if {len(s) for _, s in C.read_fasta(tmp / "backbone/units_reference_mapped.fa")} != {C.FULL_COLS}:
            raise CompFailure("incomplete_outputs", "backbone outputs not 9,090 columns")
    if not (tmp / "logs/anchor_discovery.log").stat().st_size:
        raise CompFailure("incomplete_outputs", "empty anchor_discovery.log")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--task-index", type=int, required=True)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--cpus", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 1)))
    a = ap.parse_args()
    man = pd.read_csv(C.ROOT / "manifest.validated.tsv", sep="\t", dtype=str, keep_default_na=False)
    row = man[man.task_index == str(a.task_index)]
    if len(row) != 1:
        sys.exit(f"task index {a.task_index} not in manifest")
    row = row.iloc[0].to_dict()
    for k in ("flanked_start", "flanked_end", "fasta_length", "original_start", "original_end"):
        row[k] = int(row[k])
    final = C.PER_ARRAY / row["safe_id"]
    if (final / "COMPLETE").exists() and not a.force:
        print(f"skip (COMPLETE): {row['array_id']}")
        return
    C.PER_ARRAY.mkdir(exist_ok=True)
    tmp = C.PER_ARRAY / f".tmp_{row['safe_id']}_{os.getpid()}"
    shutil.rmtree(tmp, ignore_errors=True)
    (tmp / "logs").mkdir(parents=True)
    try:
        res = process(row, tmp, a.cpus)
        validate(tmp, res)
        (tmp / "COMPLETE").write_text(f"pipeline_version={res['pipeline_version']}\nfinished={time.strftime('%Y-%m-%dT%H:%M:%S')}\n"
                                      f"status={res['primary_status']}\nreason={res['primary_reason']}\n")
        if final.exists():
            shutil.rmtree(final)         # only reached when incomplete/stale, or --force
        os.rename(tmp, final)
        print(f"OK {row['array_id']} primary={res['primary_status']} {res['primary_reason']} "
              f"anchors={res['n_anchors']} units={res['n_units_complete']} t={res['runtime_total_s']}s rss={res['peak_rss_mb']}MB")
    except CompFailure as e:
        fdir = C.ROOT / "failed" / row["safe_id"]
        shutil.rmtree(fdir, ignore_errors=True)
        fdir.parent.mkdir(exist_ok=True)
        (tmp / "FAILED.json").write_text(json.dumps(dict(array_id=row["array_id"], reason_code=e.code, message=e.msg), indent=1))
        shutil.rmtree(tmp / "scratch", ignore_errors=True)
        os.rename(tmp, fdir)
        print(f"FAILED {row['array_id']} {e}", file=sys.stderr)
        sys.exit(2)
    except Exception:
        shutil.rmtree(tmp, ignore_errors=True)
        raise


if __name__ == "__main__":
    main()
