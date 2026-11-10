#!/usr/bin/env python3
"""Recompute promoter-mapping concordance in RAW unit coordinates (post-hoc QC; no result is changed).

Why: `mafft --keeplength` deletes bases that are insertions relative to the reference alignment, so ungapped
indices of the aligned anchor / aligned unit are no longer raw positions. Here each aligned sequence is mapped
back to its raw sequence by leftmost order-preserving subsequence matching (aligned bases are a subsequence of the
raw bases), then the promoter base sets from the two alignments are compared as raw unit offsets.

  recheck_mapping.py [--out combined/promoter_mapping_recheck.tsv]
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pp_common as C


def raw_positions(aligned, raw):
    """For each non-gap column of `aligned`, the raw index of that base; None if not a subsequence."""
    pos, j, n = [], 0, len(raw)
    for ch in aligned:
        if ch == "-":
            pos.append(None)
            continue
        while j < n and raw[j] != ch:
            j += 1
        if j == n:
            return None
        pos.append(j)
        j += 1
    return pos


def window_set(pos, lo, hi):
    return {p for p in pos[lo:hi] if p is not None}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(C.COMBINED / "promoter_mapping_recheck.tsv"))
    a = ap.parse_args()
    man = pd.read_csv(C.ROOT / "manifest.validated.tsv", sep="\t", usecols=["array_id", "safe_id"])
    a0, a1 = C.PROMOTER_300[0] - 1, C.PROMOTER_300[1]
    b0, b1 = C.PROMOTER_FULL[0] - 1, C.PROMOTER_FULL[1]
    rows = []
    for _, m in man.iterrows():
        d = C.PER_ARRAY / m.safe_id
        res = pd.read_csv(d / "result.tsv", sep="\t").iloc[0]
        row = dict(array_id=m.array_id, orig_concordance=res.bb_mapping_concordance, orig_not_prefix=res.bb_anchor_not_unit_prefix)
        mapped = dict(C.read_fasta(d / "backbone/units_reference_mapped.fa"))
        if not mapped:
            rows.append(row | dict(recheck_status="no_backbone_mapping"))
            continue
        units = {n: s.upper() for n, s in C.read_fasta(d / "units_transcription.fa")}
        anchors = {n: s.upper() for n, s in C.read_fasta(d / "l1pa10_anchors.aligned.fa")}
        raw_anchor = {n: s.upper() for n, s in C.read_fasta(d / "l1pa10_anchors.fa")}
        uq = pd.read_csv(d / "units_qc.tsv", sep="\t")
        inter = union = n_units = n_fail = n_prefix_bad = n_changed = 0
        shifts = []
        for _, u in uq[uq.included_in_structural_msa].iterrows():
            uid, aid = u.unit_id, u.start_anchor_id
            n_units += 1
            raw_u, al_a, al_u = units[uid], anchors[aid].upper(), mapped[uid].upper()
            if not raw_u.startswith(raw_anchor[aid]):
                n_prefix_bad += 1
            pa, pu = raw_positions(al_a, raw_u[:len(raw_anchor[aid])]), raw_positions(al_u, raw_u)
            if pa is None or pu is None:
                n_fail += 1
                continue
            sa, sb = window_set(pa, a0, a1), window_set(pu, b0, b1)
            inter, union = inter + len(sa & sb), union + len(sa | sb)
            if sa != sb:
                n_changed += 1
                shifts.append(min(sa | sb) - min(sa & sb) if sa & sb else 999)
        row |= dict(recheck_status="ok", units_checked=n_units, units_raw_prefix_violation=n_prefix_bad, units_not_subsequence=n_fail,
                    units_with_different_promoter_base_sets=n_changed, recheck_concordance=inter / union if union else float("nan"))
        rows.append(row)
    out = pd.DataFrame(rows)
    out["recheck_flag_below_0.99"] = out.recheck_concordance < C.MAPPING_CONCORDANCE_MIN
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(a.out, sep="\t", index=False, float_format="%.6g")
    ok = out[out.recheck_status == "ok"]
    print(f"arrays rechecked {len(ok)}; raw-prefix violations {int(ok.units_raw_prefix_violation.sum())}; "
          f"units not a subsequence {int(ok.units_not_subsequence.sum())}")
    print(f"original concordance <0.99: {int((ok.orig_concordance < C.MAPPING_CONCORDANCE_MIN).sum())};  recheck <0.99: {int(ok['recheck_flag_below_0.99'].sum())}")
    print(f"original not-prefix flags: {int((ok.orig_not_prefix > 0).sum())};  recheck concordance min/median: {ok.recheck_concordance.min():.4f}/{ok.recheck_concordance.median():.4f}")


if __name__ == "__main__":
    main()
