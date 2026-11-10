#!/usr/bin/env python3
"""Edge-case integration tests running the real worker (nhmmer + MAFFT) on synthetic arrays built from pilot data.

Cases: mixed-strand array (units/anchors from both strands), no L1PA10 hits, a single anchor (no chain).
Run with the modules loaded. Does not touch per_array/.
"""
import random
import shutil
import sys
import tempfile
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pp_common as C
import run_array as RA

res = pd.read_csv("/tmp/pilot_combined/all_array_results.tsv", sep="\t")
man = pd.read_csv(C.ROOT / "manifest.validated.tsv", sep="\t")
rng = random.Random(7)


def seq_of(array_id):
    return C.read_fasta(man.loc[man.array_id == array_id, "fasta_path"].iloc[0])[0][1].upper()


minus = res[(res.n_anchors_minus > 0) & (res.n_anchors_plus == 0)].sort_values("fasta_length").iloc[0]
plus = res[(res.n_anchors_plus > 0) & (res.n_anchors_minus == 0)].sort_values("fasta_length").iloc[0]
rand = lambda n: "".join(rng.choice("ACGT") for _ in range(n))


def run_case(name, seq):
    work = Path(tempfile.mkdtemp(prefix=f"pp_edge_{name}_"))
    (work / "logs").mkdir()
    fa = work / "in.fa"
    C.write_fasta(fa, [(name, seq)])
    row = dict(man.iloc[0].to_dict(), array_id=name, safe_id=name, fasta_path=str(fa), flanked_start=1, flanked_end=len(seq), fasta_length=len(seq))
    for k in ("flanked_start", "flanked_end", "fasta_length", "original_start", "original_end"):
        row[k] = int(row[k])
    out = work / "out"
    (out / "logs").mkdir(parents=True)
    r = RA.process(row, out, 2)
    RA.validate(out, r)
    shutil.rmtree(work)
    return r


# 1. mixed strand: minus-strand array, spacer, plus-strand array
mixed = seq_of(minus.array_id) + rand(6000) + seq_of(plus.array_id)
r = run_case("mixed", mixed)
assert r["mixed_strand"] and r["n_anchors_plus"] == plus.n_anchors and r["n_anchors_minus"] == minus.n_anchors, r
assert r["n_units_adjacent"] == (plus.n_anchors - 1) + (minus.n_anchors - 1), r
assert r["primary_status"] in ("estimable", "non_estimable") and "mixed_strand_chains" in r["qc_flags"]
print(f"mixed-strand OK: anchors +{r['n_anchors_plus']} -{r['n_anchors_minus']}, units {r['n_units_complete']}, "
      f"primary={r['primary_status']}, backbone={r['bb_status']}, concordance={r['bb_mapping_concordance']}")

# 2. no hits
r = run_case("nohits", rand(60000))
assert r["primary_status"] == "non_estimable" and r["primary_reason"] == "no_qualifying_L1PA10_hits", r
assert r["bb_status"] == "non_estimable" and r["n_anchors"] == 0
print("no-hit array OK:", r["primary_reason"])

# 3. single anchor (one copy flanked by random sequence)
anc = pd.read_csv(C.PER_ARRAY / plus.safe_id / "l1pa10_anchors.tsv", sep="\t").iloc[0]
s = seq_of(plus.array_id)
single = rand(5000) + s[int(anc.local_start) - 1:int(anc.local_end)] + rand(5000)
r = run_case("single", single)
assert r["primary_status"] == "non_estimable" and r["primary_reason"] == "no_tandem_chain_with_two_anchors", r
print("single-anchor array OK:", r["primary_reason"], f"(hits pass={r['n_hits_pass_filter']})")

print("ALL EDGE-CASE TESTS PASSED")
