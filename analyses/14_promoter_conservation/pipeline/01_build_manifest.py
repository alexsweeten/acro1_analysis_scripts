#!/usr/bin/env python3
"""Build and validate manifest.validated.tsv for the 2,794 ACRO arrays.

Array definition (identical to oct_5/violin_plots.py):
  acro_telo.bed rows with last column == 0 and length >= 50 kb, named as in oct_5/moddotplot_5.txt.
Array IDs and FASTA regions carry the +/-10 kb flanks used by the legacy moddotplot / recombination runs,
so L1PA10 anchors adjacent to the ACRO1 interval are searchable. FASTA edges are not treated as repeat
boundaries anywhere downstream.
"""
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pp_common as C

OCT5 = C.SRC / "oct_5"
FLANK = 10_000
EXPECTED_N = 2794
INPUT_DIR = C.ROOT / "input_fasta"


_FAI = {}


def contig_length(contig):
    asm = assembly_for(contig)
    if asm not in _FAI:
        _FAI[asm] = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(str(asm) + ".fai")}
    return _FAI[asm][contig]


def load_set():
    bed = [l.rstrip("\n").split("\t") for l in open(OCT5 / "acro_telo.bed")]
    bed = [r for r in bed if r[-1] == "0" and int(r[2]) - int(r[1]) >= 50000]
    names = [l.split("-o ")[1].strip().split("/", 1)[1] for l in open(OCT5 / "moddotplot_5.txt")]
    if len(bed) != len(names):
        sys.exit(f"bed rows ({len(bed)}) != moddotplot_5 names ({len(names)})")
    recomb = dict(l.rstrip("\n").split("\t") for l in open(OCT5 / "recombination.txt"))
    rows = []
    for r, name in zip(bed, names):
        contig, s, e = r[0], int(r[1]), int(r[2])
        # Flanks are clamped to the contig, as in the legacy moddotplot region names.
        fs, fe = max(1, s - FLANK), min(e + FLANK, contig_length(contig))
        if name != f"{contig}_{fs}-{fe}":
            sys.exit(f"array name {name} != expected {contig}_{fs}-{fe}")
        rows.append(dict(array_id=name, chromosome_or_contig=contig, original_start=s, original_end=e,
                         flanked_start=fs, flanked_end=fe, is_not_telomeric=int(r[-1]),
                         has_recombination=int(recomb[name])))
    return pd.DataFrame(rows)


_CONTIG_ASM = {}


def assembly_for(contig):
    """Resolve by .fai membership; files are named hap1/hap2 for some samples and mat/pat for others."""
    sample = contig.split("#")[0]
    if sample not in _CONTIG_ASM:
        m = {}
        for fai in sorted((C.SRC / "assemblies" / sample).glob("*.fa.gz.fai")):
            for l in open(fai):
                m[l.split("\t")[0]] = Path(str(fai)[:-4])
        _CONTIG_ASM[sample] = m
    if contig not in _CONTIG_ASM[sample]:
        sys.exit(f"contig {contig} not found in any assembly index for {sample}")
    return _CONTIG_ASM[sample][contig]


def extract(df):
    INPUT_DIR.mkdir(exist_ok=True)
    by_asm = defaultdict(list)
    for i, r in df.iterrows():
        by_asm[assembly_for(r.chromosome_or_contig)].append(i)
    seqs = {}
    for asm, idx in by_asm.items():
        regions = "\n".join(f"{df.at[i, 'chromosome_or_contig']}:{df.at[i, 'flanked_start']}-{df.at[i, 'flanked_end']}"
                            for i in idx) + "\n"
        reg = INPUT_DIR / ".regions.txt"
        reg.write_text(regions)
        out = subprocess.run(["samtools", "faidx", str(asm), "-r", str(reg)], capture_output=True, text=True)
        if out.returncode:
            sys.exit(f"samtools faidx failed for {asm}: {out.stderr[:500]}")
        recs = []
        name, buf = None, []
        for line in out.stdout.splitlines():
            if line.startswith(">"):
                if name is not None:
                    recs.append("".join(buf))
                name, buf = line[1:], []
            else:
                buf.append(line.strip())
        if name is not None:
            recs.append("".join(buf))
        if len(recs) != len(idx):
            sys.exit(f"{asm}: asked {len(idx)} regions, got {len(recs)} records")
        for i, s in zip(idx, recs):
            seqs[i] = s
    (INPUT_DIR / ".regions.txt").unlink(missing_ok=True)
    return seqs


def main():
    df = load_set()
    if len(df) != EXPECTED_N or df.array_id.nunique() != EXPECTED_N:
        sys.exit(f"expected {EXPECTED_N} unique arrays, found {len(df)} / {df.array_id.nunique()}")
    df["safe_id"] = df.array_id.map(C.safe_id)
    if df.safe_id.nunique() != EXPECTED_N:
        sys.exit("safe_id collision")

    seqs = extract(df)
    paths, lens, shas, problems = [], [], [], []
    for i, r in df.iterrows():
        s = seqs[i]
        p = INPUT_DIR / f"{r.safe_id}.fa"
        C.write_fasta(p, [(r.array_id, s)])
        recs = C.read_fasta(p)
        if len(recs) != 1:
            problems.append((r.array_id, "not_one_record"))
        if not s:
            problems.append((r.array_id, "empty_sequence"))
        elif not C.IUPAC_RE.match(s):
            problems.append((r.array_id, "non_IUPAC_symbols"))
        paths.append(str(p))
        lens.append(len(s))
        shas.append(C.sha256_file(p))
    df["fasta_path"] = paths
    df["fasta_length"] = lens
    df["fasta_sha256"] = shas
    if df.fasta_path.duplicated().any():
        problems.append(("*", "duplicate_fasta_path"))
    if problems:
        pd.DataFrame(problems, columns=["array_id", "problem"]).to_csv(C.ROOT / "manifest.problems.tsv", sep="\t", index=False)
        sys.exit(f"{len(problems)} manifest problems; see manifest.problems.tsv")

    df["individual_id"] = df.chromosome_or_contig.str.split("#").str[0]
    df["haplotype"] = df.chromosome_or_contig.str.split("#").str[1]
    pop = pd.read_csv(C.SRC / "sample_to_superpopulation.tsv", sep="\t")
    pop = pop.rename(columns={"sample_id": "individual_id"})
    df = df.merge(pop, on="individual_id", how="left")
    df["population"] = df["population"].fillna("NA")
    df["superpopulation"] = df["superpopulation"].fillna("NA")
    df["species"] = "Homo sapiens"
    df["recombination_class"] = (df.is_not_telomeric.astype(str) + "_" + df.has_recombination.astype(str)).map({
        "0_1": "telomeric_recombination", "0_0": "telomeric_no_recombination",
        "1_0": "non_telomeric_no_recombination"}).fillna("undefined")
    df.insert(0, "task_index", range(1, len(df) + 1))
    cols = ["task_index", "array_id", "safe_id", "fasta_path", "individual_id", "haplotype", "population",
            "superpopulation", "species", "is_not_telomeric", "has_recombination", "recombination_class",
            "chromosome_or_contig", "original_start", "original_end", "flanked_start", "flanked_end",
            "fasta_length", "fasta_sha256"]
    df[cols].to_csv(C.ROOT / "manifest.validated.tsv", sep="\t", index=False)

    lines = [f"arrays\t{len(df)}", f"unique_array_id\t{df.array_id.nunique()}",
             f"unique_fasta_path\t{df.fasta_path.nunique()}", f"unique_safe_id\t{df.safe_id.nunique()}",
             f"individuals\t{df.individual_id.nunique()}",
             f"individual_haplotypes\t{(df.individual_id + '#' + df.haplotype).nunique()}",
             f"fasta_length_min/median/max\t{df.fasta_length.min()}/{int(df.fasta_length.median())}/{df.fasta_length.max()}",
             f"total_bp\t{int(df.fasta_length.sum())}",
             f"population_missing\t{(df.population == 'NA').sum()}"]
    for k, v in df.recombination_class.value_counts().items():
        lines.append(f"class:{k}\t{v}")
    (C.ROOT / "manifest.summary.tsv").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
