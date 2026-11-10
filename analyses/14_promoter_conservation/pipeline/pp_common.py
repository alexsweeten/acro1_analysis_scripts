"""Shared constants and helpers for the ACRO promoter-conservation pipeline.

Every fixed biological parameter lives here so that the worker, the aggregator and the README
generator cannot drift apart. None of these values may be tuned using the 2,794 target arrays.
"""
import hashlib
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REF = ROOT / "reference"
PER_ARRAY = ROOT / "per_array"
COMBINED = ROOT / "combined"
SRC = Path("/data/sweetenap/hprc")

# Fixed reference files: checked before anything runs.
EXPECTED_SHA256 = {
    "chr21_l1pa10.hmm": "fb15134b9446f087c8dcf27aea404758434875ce810a2e1798cb11ebfa4d0443",
    "l1pa10_reference_alignment.fa": "c47ce670e4df5d38a7bacaecf5f5fcb9440d1d0ef585cd494384cc3badeec092",
    "chr21_neg_proseq_projected_profile.tsv": "74286a783cf6d183cf556c1183d360e49d8d3da482926254c1f0d8b38ed31cf8",
    "acro_units_transcription.einsi.fa": "c0f09caf162d2bd0ff581b41f6449e07f87030bd020e9e418e62cecaaf10262b",
    "repeatmasker_breakpoints_einsi.tsv": "07ab4e62dd472a5203e839772511bd60612c1f65bb8a4038f3e31181f2a78d6d",
}
HMM = REF / "chr21_l1pa10.hmm"
REF_ALN300 = REF / "l1pa10_reference_alignment.fa"
PROSEQ = REF / "chr21_neg_proseq_projected_profile.tsv"
REF_ALN9090 = REF / "acro_units_transcription.einsi.fa"
RM_TABLE = REF / "repeatmasker_breakpoints_einsi.tsv"
BACKBONE_MASK = REF / "backbone_mask.tsv"

# ---- fixed promoter (defined before the target arrays were inspected) ----
ALN_COLS = 300
FULL_COLS = 9090
PROMOTER_300 = (78, 127)      # 1-based inclusive, 50 columns
PROMOTER_FULL = (74, 123)     # same promoter in the 9,090-column CHM13 unit alignment
WIN = 50
N_REF_300 = 138
N_REF_UNITS = 22

# ---- phase 1: anchors ----
NHMMER_E = "1e-10"
MIN_HMM_COVERAGE = 0.75
MIN_NHMMER_SCORE = 100.0
DEDUP_OVERLAP_FRAC = 0.5       # of the shorter hit
SPACING_MIN, SPACING_MAX = 3500, 20000
UNIT_MIN, UNIT_MAX = 3500, 20000

# ---- phase 2: primary matched-window test ----
MIN_CALLABLE_FRAC = 0.80
TOL_GC = 0.05
TOL_CPG = 0.02
TOL_CPG_GRID = 0.021           # labelled sensitivity only
TOL_OCC = 0.02
ROLL_DIV = 11
ROLL_PRO = 7

# ---- phase 4: conserved backbone ----
L1PA10_END_FULL = 296
ACRO1_BUFFER = 25
BACKBONE_TOL_CPG = 0.021
MIN_UNITS_BACKBONE = 3
MAPPING_CONCORDANCE_MIN = 0.99
BACKBONE_MASK_OCC = 0.80

# ---- statistics ----
SEED = 20261021
N_SIGNFLIP_PRIMARY = 200_000
N_SIGNFLIP_BACKBONE = 500_000
N_GROUP_PERM = 500_000
N_BOOT = 20_000
N_CLUSTER_BOOT = 2_000
EXACT_MAX_N = 20

IUPAC_RE = re.compile(r"^[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+$")
GAP_CHARS = "-."

_COMP = str.maketrans("ACGTRYSWKMBDHVNacgtryswkmbdhvn", "TGCAYRSWMKVHDBNtgcayrswmkvhdbn")


def revcomp(seq):
    return seq.translate(_COMP)[::-1]


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def safe_id(array_id):
    """Filesystem-safe id; the SHA256 suffix prevents collisions after sanitising."""
    base = re.sub(r"[^A-Za-z0-9._-]+", "_", array_id)
    return f"{base}__{hashlib.sha256(array_id.encode()).hexdigest()[:8]}"


def read_fasta(path):
    recs, name, buf = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n").rstrip("\r")
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(buf)))
                name, buf = line[1:].strip(), []
            elif line.strip():
                buf.append(line.strip())
    if name is not None:
        recs.append((name, "".join(buf)))
    return recs


def write_fasta(path, recs, width=80):
    with open(path, "w") as out:
        for name, seq in recs:
            out.write(f">{name}\n")
            for i in range(0, len(seq), width):
                out.write(seq[i:i + width] + "\n")


def fasta_text(recs, width=80):
    parts = []
    for name, seq in recs:
        parts.append(f">{name}\n")
        for i in range(0, len(seq), width):
            parts.append(seq[i:i + width] + "\n")
    return "".join(parts)
