#!/usr/bin/env python3
"""Stage and verify the five fixed reference files, then build the (target-independent) backbone mask.

Exits non-zero immediately if any checksum or structural expectation fails. No reference is rebuilt
or approximated from the target arrays.
"""
import shutil
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pp_common as C

OCT21 = Path("/data/sweetenap/hprc/oct_21")
ZIP = OCT21 / "reference.zip"


def fail(msg):
    sys.exit(f"REFERENCE CHECK FAILED: {msg}")


def stage():
    C.REF.mkdir(parents=True, exist_ok=True)
    for name in C.EXPECTED_SHA256:
        dst = C.REF / name
        if dst.exists():
            continue
        loose = OCT21 / name
        if loose.exists():
            shutil.copy2(loose, dst)
        elif ZIP.exists():
            with zipfile.ZipFile(ZIP) as z:
                member = f"reference/{name}"
                if member not in z.namelist():
                    fail(f"{name} not found in {ZIP}")
                with z.open(member) as src, open(dst, "wb") as out:
                    shutil.copyfileobj(src, out)
        else:
            fail(f"{name} not supplied (looked in {OCT21} and {ZIP}); please provide it")


def verify_checksums():
    rows = []
    for name, want in C.EXPECTED_SHA256.items():
        got = C.sha256_file(C.REF / name)
        rows.append((name, want, got, "OK" if got == want else "MISMATCH"))
        if got != want:
            fail(f"{name}: expected {want}, got {got}")
    pd.DataFrame(rows, columns=["file", "expected_sha256", "observed_sha256", "status"]).to_csv(
        C.REF / "reference_checksums.tsv", sep="\t", index=False)
    return rows


def verify_structure():
    aln = C.read_fasta(C.REF_ALN300)
    if len(aln) != C.N_REF_300 or {len(s) for _, s in aln} != {C.ALN_COLS}:
        fail("l1pa10_reference_alignment.fa must be 138 sequences x 300 columns")
    prof = pd.read_csv(C.PROSEQ, sep="\t")
    if len(prof) != C.ALN_COLS:
        fail("PRO-seq profile must have 300 rows")
    yes = prof.loc[prof.promoter_window == "yes", "reference_alignment_column_1based"]
    if (yes.min(), yes.max(), len(yes)) != (C.PROMOTER_300[0], C.PROMOTER_300[1], C.WIN):
        fail(f"promoter in PRO-seq profile is not 78-127: {yes.min()}-{yes.max()}")
    full = C.read_fasta(C.REF_ALN9090)
    if len(full) != C.N_REF_UNITS or {len(s) for _, s in full} != {C.FULL_COLS}:
        fail("acro_units_transcription.einsi.fa must be 22 sequences x 9,090 columns")
    rm = pd.read_csv(C.RM_TABLE, sep="\t")
    if len(rm) != 236:
        fail(f"RepeatMasker table must have 236 rows, found {len(rm)}")
    l1 = rm[rm["repeat"].str.startswith("L1PA10#")]
    if l1["alignment_end_1based_inclusive"].max() != C.L1PA10_END_FULL:
        fail("L1PA10 endpoint in RepeatMasker table is not column 296")
    return full, rm


def build_mask(full, rm):
    """Per-column flags for the 9,090-column unit alignment. Depends only on the 22 CHM13 units."""
    n = C.FULL_COLS
    occ = np.mean([[c not in C.GAP_CHARS for c in s] for _, s in full], axis=0)
    name = rm["repeat"].str.split("#").str[0]
    cls = rm["repeat"].str.split("#").str[1]
    a0, a1 = rm["alignment_start_1based"].to_numpy(), rm["alignment_end_1based_inclusive"].to_numpy()

    def paint(sel, buf=0):
        m = np.zeros(n + 2, dtype=bool)
        for s, e in zip(a0[sel], a1[sel]):
            m[max(1, s - buf):min(n, e + buf) + 1] = True
        return m[1:n + 1]

    acro1 = paint((name == "ACRO1").to_numpy())
    acro1_buf = paint((name == "ACRO1").to_numpy(), C.ACRO1_BUFFER)
    simple_sat = paint((cls.str.startswith("Simple_repeat") | cls.str.startswith("Satellite")).to_numpy())
    cols = np.arange(1, n + 1)
    low_occ = occ < C.BACKBONE_MASK_OCC
    l1pa10 = cols <= C.L1PA10_END_FULL
    base_ok = ~acro1_buf & ~simple_sat & ~low_occ
    mask = pd.DataFrame({
        "full_column_1based": cols,
        "reference_occupancy": occ.round(4),
        "in_L1PA10": l1pa10,
        "acro1_mapped": acro1,
        "acro1_buffered_25": acro1_buf,
        "simple_repeat_or_satellite": simple_sat,
        "below_80pct_reference_occupancy": low_occ,
        "primary_distal_backbone_ok": base_ok & ~l1pa10,
        "all_non_acro1_backbone_ok": base_ok,
    })
    mask.to_csv(C.BACKBONE_MASK, sep="\t", index=False)
    return mask


def main():
    stage()
    rows = verify_checksums()
    full, rm = verify_structure()
    mask = build_mask(full, rm)
    for r in rows:
        print(f"{r[3]:8s} {r[2]}  {r[0]}")
    print(f"mask: primary-ok columns={int(mask.primary_distal_backbone_ok.sum())}, "
          f"all-non-ACRO1-ok columns={int(mask.all_non_acro1_backbone_ok.sum())}, "
          f"ACRO1 columns={int(mask.acro1_mapped.sum())}")
    print("promoter: 300-col 78-127 (50), full-unit 74-123 (50); prior CHM13 PRO-seq nominated the fixed window")


if __name__ == "__main__":
    main()
