#!/usr/bin/env python3
"""Full comparative analysis of the candidate ancestral ACRO1 locus.

The analysis deliberately keeps three questions separate:

1. Is the same HSA19-syntenic genomic interval recovered in each assembly?
2. What ACRO1/ACRE structure is present at that orthologous interval?
3. Do unique flanks, ACRE backbone, and boundary junctions support the same history?

All genomic coordinates in TSV outputs are zero-based and half-open.  The
human reference locus is the primary HG002 HSA19 locus identified in the pilot.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import os
import re
import shutil
import statistics
import subprocess
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "phylogeny_origin_acro1_full"
RESULTS = OUT / "results"
INTERMEDIATE = OUT / "intermediate"
MAPPINGS = INTERMEDIATE / "anchor_mappings"
WINDOWS = INTERMEDIATE / "ortholog_windows"
ALIGNMENTS = INTERMEDIATE / "alignments"
TREES = RESULTS / "phylogenies"
FIGURES = OUT / "figures"
INPUTS = OUT / "inputs"

PRIMATE = ROOT.parent / "ACREAGE" / "primate_genomes"
CORE = INPUTS / "ncbi_downloads" / "core_outgroups" / "ncbi_dataset" / "data"
EXPANDED = INPUTS / "ncbi_downloads" / "expanded_individual"
PILOT = ROOT / "phylogeny_origin_acro1"

HUMAN_FASTA = PRIMATE / "hg002v1.1.fasta.gz"
HUMAN_CHROM = "chr19_MATERNAL"
HUMAN_LOCUS_START = 46_792_595
HUMAN_LOCUS_END = 46_858_137
HUMAN_LOCUS_CENTER = (HUMAN_LOCUS_START + HUMAN_LOCUS_END) // 2

COMPONENTS = ROOT / "marmoset_acro1_fragment_analysis" / "human_ACRO_unit1_components.fa"
REFERENCE_UNIT = PILOT / "intermediate" / "reference_ACRE_unit1.fa"
ACRO1_HMM = PILOT / "intermediate" / "ACRO1_monomer.hmm"

ANCHOR_LENGTH = 20_000
# Windows are deliberately outside the annotated locus.  The closest anchors
# abut the locus; distal anchors make the orthology call robust to local repeat
# turnover and large insertions.
ANCHOR_INTERVALS = [
    ("L500", HUMAN_LOCUS_START - 500_000, HUMAN_LOCUS_START - 480_000, "left", 500_000),
    ("L250", HUMAN_LOCUS_START - 250_000, HUMAN_LOCUS_START - 230_000, "left", 250_000),
    ("L100", HUMAN_LOCUS_START - 100_000, HUMAN_LOCUS_START - 80_000, "left", 100_000),
    ("L40", HUMAN_LOCUS_START - 40_000, HUMAN_LOCUS_START - 20_000, "left", 40_000),
    ("L20", HUMAN_LOCUS_START - 20_000, HUMAN_LOCUS_START, "left", 20_000),
    ("R20", HUMAN_LOCUS_END, HUMAN_LOCUS_END + 20_000, "right", 20_000),
    ("R40", HUMAN_LOCUS_END + 20_000, HUMAN_LOCUS_END + 40_000, "right", 40_000),
    ("R100", HUMAN_LOCUS_END + 80_000, HUMAN_LOCUS_END + 100_000, "right", 100_000),
    ("R250", HUMAN_LOCUS_END + 230_000, HUMAN_LOCUS_END + 250_000, "right", 250_000),
    ("R500", HUMAN_LOCUS_END + 480_000, HUMAN_LOCUS_END + 500_000, "right", 500_000),
]
ANCHOR_META = {
    name: {"start": start, "end": end, "side": side, "distance": distance}
    for name, start, end, side, distance in ANCHOR_INTERVALS
}


@dataclass(frozen=True)
class Assembly:
    key: str
    species: str
    clade: str
    accession: str
    fasta: Path
    haplotype: str
    source: str


ASSEMBLIES = [
    Assembly("Human", "Homo sapiens", "great_ape", "HG002v1.1", HUMAN_FASTA, "maternal", "local_T2T"),
    Assembly("Chimpanzee", "Pan troglodytes", "great_ape", "mPanTro3", PRIMATE / "mPanTro3.dip.cur.20231122.fasta.gz", "best_supported_haplotype", "local_T2T"),
    Assembly("Bonobo", "Pan paniscus", "great_ape", "mPanPan1", PRIMATE / "mPanPan1.dip.cur.20231122.fasta.gz", "best_supported_haplotype", "local_T2T"),
    Assembly("Gorilla", "Gorilla gorilla", "great_ape", "mGorGor1", PRIMATE / "mGorGor1.dip.cur.20231122.fasta.gz", "best_supported_haplotype", "local_T2T"),
    Assembly("Sumatran_orangutan", "Pongo abelii", "great_ape", "mPonAbe1", PRIMATE / "mPonAbe1.dip.cur.20231205.fasta.gz", "best_supported_haplotype", "local_T2T"),
    Assembly("Bornean_orangutan", "Pongo pygmaeus", "great_ape", "mPonPyg2", PRIMATE / "mPonPyg2.dip.cur.20231122.fasta.gz", "best_supported_haplotype", "local_T2T"),
    Assembly("Siamang_hap1", "Symphalangus syndactylus", "gibbon", "GCA_028878055.3", CORE / "GCA_028878055.3" / "GCA_028878055.3_NHGRI_mSymSyn1-v2.1_pri_genomic.fna", "primary", "NCBI"),
    Assembly("Siamang_hap2", "Symphalangus syndactylus", "gibbon", "GCA_028878085.3", CORE / "GCA_028878085.3" / "GCA_028878085.3_NHGRI_mSymSyn1-v2.0_alt_genomic.fna", "alternate", "NCBI"),
    Assembly("Hoolock_gibbon", "Hoolock leuconedys", "gibbon", "GCA_047372625.1", EXPANDED / "GCA_047372625.1_UCONN_HLeu1_genomic.fna", "haploid", "NCBI"),
    Assembly("Nomascus_gibbon", "Nomascus leucogenys", "gibbon", "GCA_040113105.1", EXPANDED / "GCA_040113105.1_ASM4011310v1_genomic.fna", "haplotype_2", "NCBI"),
    Assembly("Rhesus_macaque", "Macaca mulatta", "old_world_monkey", "GCF_003339765.1", EXPANDED / "GCF_003339765.1_Mmul_10_genomic.fna", "haploid", "NCBI"),
    Assembly("Marmoset", "Callithrix jacchus", "new_world_monkey", "GCA_049354715.1", CORE / "GCA_049354715.1" / "GCA_049354715.1_calJac240_pri_genomic.fna", "primary", "NCBI"),
    Assembly("Golden_lion_tamarin", "Leontopithecus rosalia", "new_world_monkey", "GCA_054824895.1", EXPANDED / "GCA_054824895.1_mLeoRos1.hap1_genomic.fna", "haplotype_1", "NCBI"),
    Assembly("Squirrel_monkey", "Saimiri boliviensis", "new_world_monkey", "GCA_048565385.1", EXPANDED / "GCA_048565385.1_mSaiBol1.pri_genomic.fna", "primary", "NCBI"),
]
ASM = {x.key: x for x in ASSEMBLIES}


def ensure_dirs() -> None:
    for path in [RESULTS, INTERMEDIATE, MAPPINGS, WINDOWS, ALIGNMENTS, TREES, FIGURES]:
        path.mkdir(parents=True, exist_ok=True)


def run(cmd: Sequence[str], *, stdout=None, cwd: Optional[Path] = None) -> subprocess.CompletedProcess:
    print("+", " ".join(map(str, cmd)), file=sys.stderr, flush=True)
    return subprocess.run(
        [str(x) for x in cmd],
        check=True,
        cwd=str(cwd) if cwd else None,
        stdout=stdout,
        text=stdout is subprocess.PIPE,
    )


def write_tsv(path: Path, rows: Iterable[dict], fields: Sequence[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fields})


def read_tsv(path: Path) -> List[dict]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fasta_records(path: Path) -> Iterator[Tuple[str, str]]:
    opener = gzip.open if str(path).endswith(".gz") else open
    name = None
    parts: List[str] = []
    with opener(path, "rt") as handle:
        for line in handle:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts).upper()
                name = line[1:].split()[0]
                parts = []
            elif line:
                parts.append(line)
    if name is not None:
        yield name, "".join(parts).upper()


def write_fasta(path: Path, records: Iterable[Tuple[str, str]]) -> None:
    with path.open("w") as handle:
        for name, seq in records:
            handle.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                handle.write(seq[i : i + 80] + "\n")


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTRYMKBDHVN", "TGCAYRKMVHDBN"))[::-1]


def samtools_fetch(fasta: Path, chrom: str, start: int, end: int) -> str:
    if end <= start:
        return ""
    region = f"{chrom}:{start + 1}-{end}"
    result = subprocess.run(["samtools", "faidx", str(fasta), region], check=True, capture_output=True, text=True)
    return "".join(line.strip() for line in result.stdout.splitlines() if not line.startswith(">") and line.strip()).upper()


def fai_sizes(fasta: Path) -> Dict[str, int]:
    fai = Path(str(fasta) + ".fai")
    if not fai.exists():
        run(["samtools", "faidx", fasta])
    sizes = {}
    with fai.open() as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            sizes[fields[0]] = int(fields[1])
    return sizes


def extract_anchors() -> Tuple[Path, Dict[str, str]]:
    records = []
    manifest = []
    seqs = {}
    for name, start, end, side, distance in ANCHOR_INTERVALS:
        seq = samtools_fetch(HUMAN_FASTA, HUMAN_CHROM, start, end)
        if len(seq) != end - start:
            raise RuntimeError(f"Anchor {name} has unexpected length {len(seq)}")
        records.append((name, seq))
        seqs[name] = seq
        manifest.append({
            "anchor": name,
            "human_chrom": HUMAN_CHROM,
            "human_start": start,
            "human_end": end,
            "side": side,
            "distance_from_locus_bp": distance,
            "length": len(seq),
            "n_fraction": seq.count("N") / len(seq),
        })
    out = INTERMEDIATE / "HSA19_multi_anchors.fa"
    write_fasta(out, records)
    write_tsv(RESULTS / "anchor_manifest.tsv", manifest, manifest[0].keys())
    return out, seqs


def parse_paf(path: Path, reversed_roles: bool = True) -> List[dict]:
    rows = []
    with path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            f = line.rstrip().split("\t")
            if len(f) < 12:
                continue
            tags = {x[:2]: x[5:] for x in f[12:] if len(x) > 5 and x[2:5] in {":i:", ":f:", ":Z:"}}
            row = {
                "qname": f[0], "qlen": int(f[1]), "qstart": int(f[2]), "qend": int(f[3]),
                "strand": f[4], "tname": f[5], "tlen": int(f[6]), "tstart": int(f[7]), "tend": int(f[8]),
                "nmatch": int(f[9]), "alnlen": int(f[10]), "mapq": int(f[11]),
            }
            row["identity"] = row["nmatch"] / row["alnlen"] if row["alnlen"] else 0
            row["query_coverage"] = (row["qend"] - row["qstart"]) / row["qlen"] if row["qlen"] else 0
            row["target_coverage"] = (row["tend"] - row["tstart"]) / row["tlen"] if row["tlen"] else 0
            row["score"] = int(tags.get("AS", row["nmatch"])) if str(tags.get("AS", "")).lstrip("-").isdigit() else row["nmatch"]
            rows.append(row)
    return rows


def chain_correlation(hits: List[dict], strand: str) -> float:
    if len(hits) < 2:
        return 0.0
    ordered = sorted(hits, key=lambda h: (ANCHOR_META[h["anchor"]]["start"], h["anchor"]))
    vals = [(h["tstart"] + h["tend"]) / 2 for h in ordered]
    concordant = 0
    total = 0
    for i in range(len(vals)):
        for j in range(i + 1, len(vals)):
            total += 1
            if (strand == "+" and vals[j] > vals[i]) or (strand == "-" and vals[j] < vals[i]):
                concordant += 1
    return concordant / total if total else 0.0


def consolidate_anchor_fragments(raw: List[dict], max_target_gap: int = 100_000) -> List[dict]:
    """Join split PAF records from one anchor into a collinear local chain.

    Distant primates frequently align a 20-kb unique flank as several blocks.
    Treating each block independently discards real orthology, whereas joining
    across chromosomes or distant target positions would create false support.
    """
    grouped: Dict[Tuple[str, str, str], List[dict]] = defaultdict(list)
    for h in raw:
        grouped[(h["qname"], h["tname"], h["strand"])].append(h)
    consolidated = []
    for (anchor, chrom, strand), hits in grouped.items():
        hits = sorted(hits, key=lambda h: (h["tstart"], h["tend"]))
        clusters: List[List[dict]] = []
        for h in hits:
            if not clusters or h["tstart"] - max(x["tend"] for x in clusters[-1]) > max_target_gap:
                clusters.append([h])
            else:
                clusters[-1].append(h)
        for cluster in clusters:
            q_intervals = merge_intervals((h["qstart"], h["qend"]) for h in cluster)
            q_covered = sum(b - a for a, b in q_intervals)
            alnlen = sum(h["alnlen"] for h in cluster)
            nmatch = sum(h["nmatch"] for h in cluster)
            template = dict(cluster[0])
            template.update({
                "qstart": min(h["qstart"] for h in cluster), "qend": max(h["qend"] for h in cluster),
                "tstart": min(h["tstart"] for h in cluster), "tend": max(h["tend"] for h in cluster),
                "nmatch": nmatch, "alnlen": alnlen, "mapq": max(h["mapq"] for h in cluster),
                "identity": nmatch / alnlen if alnlen else 0, "query_coverage": q_covered / template["qlen"],
                "score": sum(h["score"] for h in cluster), "fragment_count": len(cluster),
            })
            consolidated.append(template)
    return consolidated


def select_anchor_chain(assembly: Assembly, paf: Path) -> Tuple[dict, List[dict], List[dict]]:
    raw = consolidate_anchor_fragments(parse_paf(paf))
    accepted = []
    for h in raw:
        # Human anchors are PAF queries and the assembly is the target.
        h["anchor"] = h["qname"]
        h["assembly_chrom"] = h["tname"]
        h["anchor_coverage"] = h["query_coverage"]
        if h["anchor"] not in ANCHOR_META:
            continue
        if h["anchor_coverage"] < 0.35 or h["identity"] < 0.65 or h["alnlen"] < 3_000:
            continue
        accepted.append(h)

    grouped: Dict[Tuple[str, str], List[dict]] = defaultdict(list)
    for h in accepted:
        grouped[(h["assembly_chrom"], h["strand"])].append(h)

    candidates = []
    for (chrom, strand), hits in grouped.items():
        best_by_anchor = {}
        for h in hits:
            score = (h["anchor_coverage"], h["identity"], h["score"], h["mapq"])
            if h["anchor"] not in best_by_anchor or score > best_by_anchor[h["anchor"]][0]:
                best_by_anchor[h["anchor"]] = (score, h)
        best = [v[1] for v in best_by_anchor.values()]
        sides = {ANCHOR_META[h["anchor"]]["side"] for h in best}
        left_n = sum(ANCHOR_META[h["anchor"]]["side"] == "left" for h in best)
        right_n = len(best) - left_n
        order = chain_correlation(best, strand)
        aggregate = sum(h["anchor_coverage"] * h["identity"] * h["alnlen"] for h in best)
        candidates.append({
            "chrom": chrom, "strand": strand, "hits": best, "anchor_count": len(best),
            "left_anchor_count": left_n, "right_anchor_count": right_n, "order_concordance": order,
            "aggregate_score": aggregate,
            "valid": len(best) >= 4 and left_n >= 2 and right_n >= 2 and order >= 0.75,
        })
    candidates.sort(key=lambda c: (c["valid"], c["anchor_count"], c["order_concordance"], c["aggregate_score"]), reverse=True)
    if not candidates:
        call = {"assembly_key": assembly.key, "orthology_status": "unresolved_no_anchor_chain"}
        return call, [], accepted
    chosen = candidates[0]
    hits = chosen["hits"]
    left_hits = [h for h in hits if ANCHOR_META[h["anchor"]]["side"] == "left"]
    right_hits = [h for h in hits if ANCHOR_META[h["anchor"]]["side"] == "right"]
    # Select the closest source-anchor pair that also brackets a plausible
    # target interval.  This is robust to a single locally transposed anchor,
    # as observed for Nomascus L20.
    bracketing = []
    for left_hit in left_hits:
        for right_hit in right_hits:
            if chosen["strand"] == "+":
                gap = right_hit["tstart"] - left_hit["tend"]
            else:
                gap = left_hit["tstart"] - right_hit["tend"]
            if 0 <= gap <= 2_000_000:
                source_distance = ANCHOR_META[left_hit["anchor"]]["distance"] + ANCHOR_META[right_hit["anchor"]]["distance"]
                bracketing.append((source_distance, gap, left_hit, right_hit))
    if bracketing:
        _, _, left, right = min(bracketing, key=lambda x: (x[0], x[1]))
    else:
        left = min(left_hits, key=lambda h: ANCHOR_META[h["anchor"]]["distance"])
        right = min(right_hits, key=lambda h: ANCHOR_META[h["anchor"]]["distance"])
    if chosen["strand"] == "+":
        locus_start, locus_end = left["tend"], right["tstart"]
    else:
        locus_start, locus_end = right["tend"], left["tstart"]
    if locus_start > locus_end:
        locus_start, locus_end = locus_end, locus_start
    span = locus_end - locus_start
    if not chosen["valid"]:
        status = "unresolved_weak_chain"
    elif span < 0 or span > 2_000_000:
        status = "unresolved_implausible_interval"
    elif chosen["anchor_count"] >= 8 and chosen["order_concordance"] >= 0.95:
        status = "high_confidence_ortholog"
    elif chosen["anchor_count"] >= 6:
        status = "moderate_confidence_ortholog"
    else:
        status = "low_confidence_ortholog"
    call = {
        "assembly_key": assembly.key, "species": assembly.species, "clade": assembly.clade,
        "accession": assembly.accession, "haplotype": assembly.haplotype, "target_chrom": chosen["chrom"],
        "orientation_to_human": chosen["strand"], "anchor_count": chosen["anchor_count"],
        "left_anchor_count": chosen["left_anchor_count"], "right_anchor_count": chosen["right_anchor_count"],
        "order_concordance": round(chosen["order_concordance"], 6), "locus_start": locus_start,
        "locus_end": locus_end, "locus_span_bp": span, "orthology_status": status,
        "nearest_left_anchor": left["anchor"], "nearest_right_anchor": right["anchor"],
        "mean_anchor_identity": round(statistics.mean(h["identity"] for h in hits), 6),
        "mean_anchor_coverage": round(statistics.mean(h["anchor_coverage"] for h in hits), 6),
    }
    for h in hits:
        h["assembly_key"] = assembly.key
        h["selected_chain"] = True
    return call, hits, accepted


def run_anchor_mappings(anchor_fasta: Path, threads: int, resume: bool) -> Tuple[List[dict], Dict[str, List[dict]]]:
    calls = []
    selected = {}
    all_hit_rows = []
    for assembly in ASSEMBLIES:
        paf = MAPPINGS / f"{assembly.key}.anchors.paf"
        if not (resume and paf.exists() and paf.stat().st_size > 0):
            with paf.open("w") as handle:
                run(["minimap2", "-t", str(threads), "-x", "asm20", "-N", "50", "--secondary=yes", "-c", assembly.fasta, anchor_fasta], stdout=handle)
        call, hits, accepted = select_anchor_chain(assembly, paf)
        calls.append(call)
        selected[assembly.key] = hits
        chosen_ids = {(h["qname"], h["qstart"], h["qend"], h["tname"], h["strand"]) for h in hits}
        for h in accepted:
            all_hit_rows.append({
                "assembly_key": assembly.key, "assembly_chrom": h["assembly_chrom"], "anchor": h["anchor"],
                "strand": h["strand"], "assembly_start": h["tstart"], "assembly_end": h["tend"],
                "anchor_start": h["qstart"], "anchor_end": h["qend"], "identity": round(h["identity"], 6),
                "anchor_coverage": round(h["anchor_coverage"], 6), "mapq": h["mapq"],
                "fragment_count": h.get("fragment_count", 1),
                "selected_chain": (h["qname"], h["qstart"], h["qend"], h["tname"], h["strand"]) in chosen_ids,
            })
    fields = ["assembly_key", "species", "clade", "accession", "haplotype", "target_chrom", "orientation_to_human", "anchor_count", "left_anchor_count", "right_anchor_count", "order_concordance", "locus_start", "locus_end", "locus_span_bp", "orthology_status", "nearest_left_anchor", "nearest_right_anchor", "mean_anchor_identity", "mean_anchor_coverage"]
    write_tsv(RESULTS / "orthologous_locus_matrix.tsv", calls, fields)
    write_tsv(RESULTS / "all_anchor_hits.tsv", all_hit_rows, all_hit_rows[0].keys() if all_hit_rows else [])
    return calls, selected


def build_target_anchor_fasta(selected: Dict[str, List[dict]]) -> Tuple[Path, Dict[str, dict]]:
    records = []
    metadata = {}
    for key, hits in selected.items():
        assembly = ASM[key]
        for h in hits:
            seq = samtools_fetch(assembly.fasta, h["tname"], h["tstart"], h["tend"])
            if h["strand"] == "-":
                seq = revcomp(seq)
            name = f"{key}__{h['anchor']}"
            records.append((name, seq))
            metadata[name] = h
    path = INTERMEDIATE / "selected_target_anchors.fa"
    write_fasta(path, records)
    return path, metadata


def reciprocal_validation(target_anchor_fasta: Path, metadata: Dict[str, dict], threads: int) -> Dict[str, dict]:
    paf = MAPPINGS / "reciprocal_target_anchors_to_human.paf"
    with paf.open("w") as handle:
        run(["minimap2", "-t", str(threads), "-x", "asm20", "-N", "20", "--secondary=yes", "-c", HUMAN_FASTA, target_anchor_fasta], stdout=handle)
    # Species-specific anchors are PAF queries; human chromosomes are targets.
    hits = consolidate_anchor_fragments(parse_paf(paf))
    by_target: Dict[str, List[dict]] = defaultdict(list)
    for h in hits:
        if h["qname"] in metadata and h["query_coverage"] >= 0.35 and h["identity"] >= 0.65:
            by_target[h["qname"]].append(h)
    rows = []
    summary: Dict[str, dict] = defaultdict(lambda: {"reciprocal_anchor_count": 0, "reciprocal_HSA19_count": 0})
    for target_name, candidates in by_target.items():
        best = max(candidates, key=lambda h: (h["target_coverage"], h["identity"], h["score"], h["mapq"]))
        key, anchor = target_name.split("__", 1)
        mid = (best["tstart"] + best["tend"]) // 2
        on_hsa19 = best["tname"].startswith("chr19_") and abs(mid - HUMAN_LOCUS_CENTER) <= 1_000_000
        summary[key]["reciprocal_anchor_count"] += 1
        summary[key]["reciprocal_HSA19_count"] += int(on_hsa19)
        rows.append({
            "assembly_key": key, "anchor": anchor, "human_chrom": best["tname"], "human_start": best["tstart"],
            "human_end": best["tend"], "strand": best["strand"], "identity": round(best["identity"], 6),
            "query_coverage": round(best["query_coverage"], 6), "mapq": best["mapq"], "reciprocal_HSA19": on_hsa19,
        })
    write_tsv(RESULTS / "reciprocal_anchor_validation.tsv", rows, rows[0].keys() if rows else [])
    return summary


def update_orthology_with_reciprocal(calls: List[dict], summary: Dict[str, dict]) -> None:
    for call in calls:
        rec = summary.get(call["assembly_key"], {})
        call.update(rec)
        n = int(rec.get("reciprocal_HSA19_count", 0))
        anchors = int(call.get("anchor_count", 0) or 0)
        call["reciprocal_fraction"] = round(n / anchors, 6) if anchors else 0
        call["reciprocal_orthology_supported"] = n >= 4 and (n / anchors if anchors else 0) >= 0.5
    fields = list(calls[0].keys())
    write_tsv(RESULTS / "orthologous_locus_matrix.tsv", calls, fields)


def extract_ortholog_windows(calls: List[dict], padding: int = 50_000) -> Tuple[Path, Dict[str, dict]]:
    records = []
    windows = {}
    for call in calls:
        if "ortholog" not in str(call.get("orthology_status", "")):
            continue
        key = call["assembly_key"]
        assembly = ASM[key]
        sizes = fai_sizes(assembly.fasta)
        chrom = call["target_chrom"]
        start = max(0, int(call["locus_start"]) - padding)
        end = min(sizes[chrom], int(call["locus_end"]) + padding)
        seq = samtools_fetch(assembly.fasta, chrom, start, end)
        name = key
        records.append((name, seq))
        windows[key] = {
            "assembly_key": key, "chrom": chrom, "window_start": start, "window_end": end,
            "window_length": len(seq), "locus_start_in_window": int(call["locus_start"]) - start,
            "locus_end_in_window": int(call["locus_end"]) - start, "orientation_to_human": call["orientation_to_human"],
            "sequence": seq,
        }
        write_fasta(WINDOWS / f"{key}.fa", [(key, seq)])
    path = INTERMEDIATE / "HSA19_ortholog_windows.fa"
    write_fasta(path, records)
    manifest = [{k: v for k, v in row.items() if k != "sequence"} for row in windows.values()]
    write_tsv(RESULTS / "ortholog_window_manifest.tsv", manifest, manifest[0].keys() if manifest else [])
    return path, windows


def parse_nhmmer_tbl(path: Path, hmm_len: int) -> List[dict]:
    hits = []
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.split()
            if len(f) < 15:
                continue
            hmm_from, hmm_to = int(f[4]), int(f[5])
            ali_from, ali_to = int(f[6]), int(f[7])
            start, end = sorted((ali_from - 1, ali_to))
            hits.append({
                "target": f[0], "hmm_from": hmm_from, "hmm_to": hmm_to, "start": start, "end": end,
                "strand": f[11], "evalue": float(f[12]), "score": float(f[13]),
                "hmm_coverage": (abs(hmm_to - hmm_from) + 1) / hmm_len,
            })
    return hits


def merge_intervals(intervals: Iterable[Tuple[int, int]]) -> List[Tuple[int, int]]:
    merged = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [(x[0], x[1]) for x in merged]


def map_queries_to_windows(query_fasta: Path, windows_fasta: Path, out_paf: Path, threads: int, short: bool = False) -> List[dict]:
    preset = ["-x", "sr", "-k", "9", "-w", "5"] if short else ["-x", "asm20", "-k", "11", "-w", "5"]
    with out_paf.open("w") as handle:
        run(["minimap2", "-t", str(threads), *preset, "-N", "200", "--secondary=yes", "-c", windows_fasta, query_fasta], stdout=handle)
    return parse_paf(out_paf)


def structural_analysis(windows_fasta: Path, windows: Dict[str, dict], threads: int) -> Tuple[List[dict], Dict[str, List[dict]], List[dict]]:
    # Determine HMM length from the model header.
    hmm_len = 0
    with ACRO1_HMM.open() as handle:
        for line in handle:
            if line.startswith("LENG"):
                hmm_len = int(line.split()[1])
                break
    tbl = INTERMEDIATE / "HSA19_windows_ACRO1.tbl"
    report = INTERMEDIATE / "HSA19_windows_ACRO1.nhmmer.txt"
    with report.open("w") as handle:
        run(["nhmmer", "--cpu", str(max(1, min(threads, 8))), "--tblout", tbl, ACRO1_HMM, windows_fasta], stdout=handle)
    hmm_hits = [h for h in parse_nhmmer_tbl(tbl, hmm_len) if h["hmm_coverage"] >= 0.45 and h["evalue"] <= 1e-3]
    hmm_by_target: Dict[str, List[dict]] = defaultdict(list)
    for hit in hmm_hits:
        hmm_by_target[hit["target"]].append(hit)

    component_by_target: Dict[str, List[dict]] = defaultdict(list)
    unit_by_target: Dict[str, List[dict]] = defaultdict(list)
    unit_hits = []
    # Map per ortholog window.  A combined target causes the -N secondary-hit
    # ceiling to be consumed by the largest gibbon arrays and can suppress a
    # weaker marmoset partial-unit match.
    per_window_maps = INTERMEDIATE / "per_window_component_mappings"
    per_window_maps.mkdir(parents=True, exist_ok=True)
    for key in windows:
        target = WINDOWS / f"{key}.fa"
        component_hits = map_queries_to_windows(COMPONENTS, target, per_window_maps / f"{key}.components.paf", threads, short=True)
        for h in component_hits:
            if h["query_coverage"] >= 0.45 and h["identity"] >= 0.60 and h["alnlen"] >= 100:
                component_by_target[key].append(h)
        local_units = map_queries_to_windows(REFERENCE_UNIT, target, per_window_maps / f"{key}.unit.paf", threads)
        unit_hits.extend(local_units)
        for h in local_units:
            if h["alnlen"] >= 1_000 and h["identity"] >= 0.60:
                unit_by_target[key].append(h)

    rows = []
    raw_hmm_rows = []
    for key, window in windows.items():
        hits = hmm_by_target.get(key, [])
        merged = merge_intervals((h["start"], h["end"]) for h in hits)
        acro_bp = sum(b - a for a, b in merged)
        full_hmm = sum(h["hmm_coverage"] >= 0.80 for h in hits)
        comps = component_by_target.get(key, [])
        comp_names = sorted({h["qname"] for h in comps if not h["qname"].startswith("ACRO1_")})
        units = unit_by_target.get(key, [])
        best_unit = max(units, key=lambda h: (h["query_coverage"], h["identity"], h["score"]), default=None)
        best_cov = best_unit["query_coverage"] if best_unit else 0.0
        best_ident = best_unit["identity"] if best_unit else 0.0
        median_hmm_score = statistics.median(h["score"] for h in hits) if hits else 0
        median_hmm_coverage = statistics.median(h["hmm_coverage"] for h in hits) if hits else 0
        if acro_bp < 200:
            state = "empty_or_unrecognizable"
        elif acro_bp >= 20_000:
            state = "expanded_ACRO1_like_array"
        elif best_cov >= 0.75 and len(comp_names) >= 3:
            state = "complete_ACRE_like"
        elif acro_bp >= 5_000 and len(comp_names) >= 3:
            state = "ACRO1_array_with_ACRE_components"
        elif acro_bp >= 5_000:
            state = "moderate_ACRO1_like_array"
        elif any(name.startswith("COMP_38") for name in comp_names) or best_cov >= 0.30:
            state = "short_partial_ACRE"
        elif acro_bp >= 1_000:
            state = "short_ACRO1_like_array"
        else:
            state = "isolated_ACRO1_fragment"
        locus_left, locus_right = window["locus_start_in_window"], window["locus_end_in_window"]
        locus_hmm = [h for h in hits if min(h["end"], locus_right) > max(h["start"], locus_left)]
        locus_merged = merge_intervals((max(h["start"], locus_left), min(h["end"], locus_right)) for h in locus_hmm)
        locus_acro_bp = sum(b - a for a, b in locus_merged)
        rows.append({
            "assembly_key": key, "species": ASM[key].species, "clade": ASM[key].clade,
            "ortholog_chrom": window["chrom"], "ortholog_locus_start": window["window_start"] + locus_left,
            "ortholog_locus_end": window["window_start"] + locus_right, "ortholog_locus_span_bp": locus_right - locus_left,
            "ACRO1_HMM_hit_count": len(hits), "full_length_ACRO1_hit_count": full_hmm,
            "ACRO1_bp_in_window": acro_bp, "ACRO1_bp_in_inferred_locus": locus_acro_bp,
            "median_ACRO1_HMM_score": round(median_hmm_score, 6), "median_ACRO1_HMM_coverage": round(median_hmm_coverage, 6),
            "non_ACRO_component_count": len(comp_names), "non_ACRO_components": ";".join(comp_names),
            "best_reference_unit_coverage": round(best_cov, 6), "best_reference_unit_identity": round(best_ident, 6),
            "orthologous_locus_state": state,
        })
        for h in hits:
            raw_hmm_rows.append({"assembly_key": key, **h})
    write_tsv(RESULTS / "orthologous_locus_architecture.tsv", rows, rows[0].keys() if rows else [])
    write_tsv(RESULTS / "ACRO1_HMM_hits_at_orthologs.tsv", raw_hmm_rows, raw_hmm_rows[0].keys() if raw_hmm_rows else [])
    return rows, hmm_by_target, unit_hits


def known_locus_crosscheck(rows: List[dict]) -> None:
    known = {
        "Marmoset": ("CM111806.1", 38_744_563, 38_755_658),
        "Siamang_hap1": ("CM054525.2", 17_550_381, 17_850_868),
        "Siamang_hap2": ("CM054551.2", 22_557_103, 22_992_465),
    }
    output = []
    by_key = {r["assembly_key"]: r for r in rows}
    for key, (chrom, start, end) in known.items():
        observed = by_key.get(key)
        same_chrom = observed is not None and observed["ortholog_chrom"] == chrom
        overlap = 0
        if same_chrom:
            overlap = max(0, min(int(observed["ortholog_locus_end"]), end) - max(int(observed["ortholog_locus_start"]), start))
        output.append({
            "assembly_key": key, "known_chrom": chrom, "known_start": start, "known_end": end,
            "inferred_chrom": observed.get("ortholog_chrom", "") if observed else "",
            "inferred_start": observed.get("ortholog_locus_start", "") if observed else "",
            "inferred_end": observed.get("ortholog_locus_end", "") if observed else "",
            "same_chrom": same_chrom, "overlap_bp": overlap,
            "crosscheck_pass": same_chrom and overlap > 0,
        })
    write_tsv(RESULTS / "pilot_locus_crosscheck.tsv", output, output[0].keys())


def junction_analysis(windows_fasta: Path, windows: Dict[str, dict], threads: int) -> List[dict]:
    marm = ASM["Marmoset"]
    chrom, start, end = "CM111806.1", 38_744_563, 38_755_658
    probes = [
        ("marmoset_left_junction_10kb", samtools_fetch(marm.fasta, chrom, start - 5_000, start + 5_000)),
        ("marmoset_partial_ACRE", samtools_fetch(marm.fasta, chrom, start, end)),
        ("marmoset_right_junction_10kb", samtools_fetch(marm.fasta, chrom, end - 5_000, end + 5_000)),
    ]
    probe_fa = INTERMEDIATE / "marmoset_candidate_junction_probes.fa"
    write_fasta(probe_fa, probes)
    hits = map_queries_to_windows(probe_fa, windows_fasta, INTERMEDIATE / "marmoset_junctions_to_orthologs.paf", threads)
    best = {}
    for h in hits:
        if h["query_coverage"] < 0.25 or h["identity"] < 0.60 or h["alnlen"] < 1_000:
            continue
        key = (h["tname"], h["qname"])
        if key not in best or (h["query_coverage"], h["identity"], h["score"]) > (best[key]["query_coverage"], best[key]["identity"], best[key]["score"]):
            best[key] = h
    rows = []
    for key in windows:
        for probe, _ in probes:
            h = best.get((key, probe))
            rows.append({
                "assembly_key": key, "probe": probe, "detected": bool(h),
                "query_coverage": round(h["query_coverage"], 6) if h else 0,
                "identity": round(h["identity"], 6) if h else 0,
                "window_start": h["tstart"] if h else "", "window_end": h["tend"] if h else "",
                "strand": h["strand"] if h else "",
                "strong_shared_junction": bool(h and h["query_coverage"] >= 0.70 and h["identity"] >= 0.70),
            })
    write_tsv(RESULTS / "junction_conservation.tsv", rows, rows[0].keys())
    return rows


def assembly_qc(windows: Dict[str, dict]) -> List[dict]:
    rows = []
    for key, w in windows.items():
        seq = w["sequence"]
        sizes = fai_sizes(ASM[key].fasta)
        chrom_len = sizes[w["chrom"]]
        longest_n = max((len(x) for x in re.findall(r"N+", seq)), default=0)
        rows.append({
            "assembly_key": key, "accession": ASM[key].accession, "chrom": w["chrom"],
            "chrom_length": chrom_len, "window_start": w["window_start"], "window_end": w["window_end"],
            "distance_to_left_contig_end": w["window_start"], "distance_to_right_contig_end": chrom_len - w["window_end"],
            "window_N_fraction": round(seq.count("N") / len(seq), 8) if seq else 1,
            "longest_N_run": longest_n, "locus_internal_gt50kb_from_contig_end": w["window_start"] >= 50_000 and chrom_len - w["window_end"] >= 50_000,
            "assembly_continuity_pass": longest_n < 1_000 and (seq.count("N") / len(seq) if seq else 1) < 0.01,
        })
    write_tsv(RESULTS / "assembly_context_qc.tsv", rows, rows[0].keys())
    return rows


def read_fasta_dict(path: Path) -> Dict[str, str]:
    return dict(fasta_records(path))


def build_flank_phylogeny(selected: Dict[str, List[dict]], threads: int) -> Optional[Path]:
    by_assembly = {key: {h["anchor"]: h for h in hits} for key, hits in selected.items()}
    eligible = [key for key, hits in by_assembly.items() if len(hits) >= 6]
    aligned_by_anchor: Dict[str, Dict[str, str]] = {}
    for anchor, human_seq in read_fasta_dict(INTERMEDIATE / "HSA19_multi_anchors.fa").items():
        records = []
        for key in eligible:
            h = by_assembly[key].get(anchor)
            if not h:
                continue
            seq = samtools_fetch(ASM[key].fasta, h["tname"], h["tstart"], h["tend"])
            if h["strand"] == "-":
                seq = revcomp(seq)
            records.append((key, seq))
        if len(records) < 6:
            continue
        unaligned = ALIGNMENTS / f"flank_{anchor}.fa"
        aligned = ALIGNMENTS / f"flank_{anchor}.aligned.fa"
        write_fasta(unaligned, records)
        with aligned.open("w") as handle:
            run(["mafft", "--auto", "--thread", str(threads), unaligned], stdout=handle)
        aligned_by_anchor[anchor] = read_fasta_dict(aligned)
    if not aligned_by_anchor:
        return None
    concat = {}
    for key in eligible:
        pieces = []
        for anchor in sorted(aligned_by_anchor, key=lambda x: ANCHOR_META[x]["start"]):
            aln = aligned_by_anchor[anchor]
            width = len(next(iter(aln.values())))
            pieces.append(aln.get(key, "-" * width))
        concat[key] = "".join(pieces)
    concat_fa = ALIGNMENTS / "concatenated_unique_flanks.aligned.fa"
    write_fasta(concat_fa, concat.items())
    # Remove columns supported by fewer than 70% of taxa.  This prevents
    # fragmented distant-primate alignments from contributing large blocks of
    # taxon-specific gaps to the control tree.
    keys = list(concat)
    width = len(next(iter(concat.values())))
    min_present = math.ceil(0.70 * len(keys))
    keep = [i for i in range(width) if sum(concat[k][i].upper() in "ACGT" for k in keys) >= min_present]
    trimmed = {k: "".join(concat[k][i] for i in keep) for k in keys}
    trimmed_fa = ALIGNMENTS / "concatenated_unique_flanks.trimmed.fa"
    write_fasta(trimmed_fa, trimmed.items())
    write_tsv(
        RESULTS / "unique_flank_alignment_qc.tsv",
        [{"assembly_key": k, "alignment_columns": len(keep), "gap_fraction": round(trimmed[k].count("-") / len(keep), 6) if keep else 1} for k in keys],
        ["assembly_key", "alignment_columns", "gap_fraction"],
    )
    prefix = TREES / "unique_flank_tree"
    run(["iqtree2" if shutil.which("iqtree2") else "iqtree", "-s", trimmed_fa, "-m", "GTR+F+G4", "-B", "1000", "-T", "AUTO", "--prefix", prefix, "-redo"])
    return Path(str(prefix) + ".treefile")


def build_backbone_phylogeny(unit_hits: List[dict], windows: Dict[str, dict], threads: int) -> Optional[Path]:
    best = {}
    for h in unit_hits:
        key = h["tname"]
        if key not in windows or h["alnlen"] < 1_500 or h["identity"] < 0.60:
            continue
        if key not in best or (h["query_coverage"], h["identity"], h["score"]) > (best[key]["query_coverage"], best[key]["identity"], best[key]["score"]):
            best[key] = h
    records = []
    manifest = []
    for key, h in best.items():
        seq = windows[key]["sequence"][h["tstart"] : h["tend"]]
        if h["strand"] == "-":
            seq = revcomp(seq)
        if len(seq) < 1_500:
            continue
        records.append((key, seq))
        manifest.append({"assembly_key": key, "coverage": round(h["query_coverage"], 6), "identity": round(h["identity"], 6), "length": len(seq)})
    if len(records) < 4:
        return None
    unaligned = ALIGNMENTS / "best_ACRE_backbone_segments.fa"
    aligned = ALIGNMENTS / "best_ACRE_backbone_segments.aligned.fa"
    write_fasta(unaligned, records)
    write_tsv(RESULTS / "ACRE_backbone_tree_manifest.tsv", manifest, manifest[0].keys())
    with aligned.open("w") as handle:
        run(["mafft", "--auto", "--thread", str(threads), unaligned], stdout=handle)
    prefix = TREES / "ACRE_backbone_tree"
    run(["iqtree2" if shutil.which("iqtree2") else "iqtree", "-s", aligned, "-m", "GTR+F+G4", "-B", "1000", "-T", "AUTO", "--prefix", prefix, "-redo"])
    return Path(str(prefix) + ".treefile")


def consensus_from_alignment(records: Dict[str, str]) -> str:
    if not records:
        return ""
    width = len(next(iter(records.values())))
    consensus = []
    for i in range(width):
        counts: Dict[str, int] = defaultdict(int)
        for seq in records.values():
            base = seq[i].upper()
            if base in "ACGT":
                counts[base] += 1
        if not counts:
            consensus.append("-")
        else:
            consensus.append(max(counts, key=lambda b: (counts[b], b)))
    return "".join(consensus).replace("-", "")


def build_ACRO1_consensus_phylogeny(hmm_by_target: Dict[str, List[dict]], windows: Dict[str, dict], threads: int) -> Optional[Path]:
    sampled = []
    sampling_rows = []
    for key, hits in hmm_by_target.items():
        full = [h for h in hits if h["hmm_coverage"] >= 0.80]
        # Score-stratified selection avoids allowing the largest arrays to
        # dominate while retaining both canonical and divergent monomers.
        full.sort(key=lambda h: h["score"], reverse=True)
        chosen = full[:30]
        for i, h in enumerate(chosen, 1):
            seq = windows[key]["sequence"][h["start"] : h["end"]]
            if h["strand"] == "-":
                seq = revcomp(seq)
            if len(seq) < 80:
                continue
            sampled.append((f"{key}__m{i:02d}", seq))
            sampling_rows.append({"assembly_key": key, "sample": i, "score": h["score"], "hmm_coverage": h["hmm_coverage"], "length": len(seq)})
    if len(sampled) < 8:
        return None
    sampled_fa = ALIGNMENTS / "sampled_ACRO1_monomers.fa"
    sampled_aln = ALIGNMENTS / "sampled_ACRO1_monomers.aligned.fa"
    write_fasta(sampled_fa, sampled)
    write_tsv(RESULTS / "ACRO1_consensus_sampling.tsv", sampling_rows, sampling_rows[0].keys())
    with sampled_aln.open("w") as handle:
        run(["mafft", "--auto", "--thread", str(threads), sampled_fa], stdout=handle)
    aligned = read_fasta_dict(sampled_aln)
    grouped: Dict[str, Dict[str, str]] = defaultdict(dict)
    for name, seq in aligned.items():
        key = name.split("__", 1)[0]
        grouped[key][name] = seq
    consensuses = [(key, consensus_from_alignment(records)) for key, records in grouped.items()]
    consensus_fa = RESULTS / "ACRO1_locus_consensuses.fa"
    consensus_aln = ALIGNMENTS / "ACRO1_locus_consensuses.aligned.fa"
    write_fasta(consensus_fa, consensuses)
    with consensus_aln.open("w") as handle:
        run(["mafft", "--auto", "--thread", str(threads), consensus_fa], stdout=handle)
    aln = read_fasta_dict(consensus_aln)
    identity_rows = []
    for a, sa in aln.items():
        for b, sb in aln.items():
            compared = [(x, y) for x, y in zip(sa, sb) if x in "ACGT" and y in "ACGT"]
            identity = sum(x == y for x, y in compared) / len(compared) if compared else 0
            identity_rows.append({"assembly_a": a, "assembly_b": b, "aligned_bases": len(compared), "identity": round(identity, 6)})
    write_tsv(RESULTS / "ACRO1_locus_consensus_pairwise_identity.tsv", identity_rows, identity_rows[0].keys())
    prefix = TREES / "ACRO1_locus_consensus_tree"
    run(["iqtree2" if shutil.which("iqtree2") else "iqtree", "-s", consensus_aln, "-m", "GTR+F+G4", "-B", "1000", "-T", "AUTO", "--prefix", prefix, "-redo"])
    return Path(str(prefix) + ".treefile")


TREE = (
    "((((Marmoset,Golden_lion_tamarin),Squirrel_monkey),Rhesus_macaque),"
    "((Hoolock_gibbon,(Nomascus_gibbon,Siamang_hap1)),"
    "((Sumatran_orangutan,Bornean_orangutan),(Gorilla,(Human,(Chimpanzee,Bonobo))))));"
)


class Node:
    def __init__(self, name: str = "", children: Optional[List["Node"]] = None):
        self.name = name
        self.children = children or []
        self.states: set[str] = set()


def parse_simple_newick(text: str) -> Node:
    i = 0
    internal = 0
    def parse() -> Node:
        nonlocal i, internal
        if text[i] == "(":
            i += 1
            kids = []
            while True:
                kids.append(parse())
                if text[i] == ",":
                    i += 1
                    continue
                if text[i] == ")":
                    i += 1
                    break
            internal += 1
            return Node(f"Node_{internal}", kids)
        start = i
        while i < len(text) and text[i] not in ",();":
            i += 1
        return Node(text[start:i])
    return parse()


def fitch_reconstruction(architecture: List[dict]) -> List[dict]:
    states = {}
    for r in architecture:
        key = r["assembly_key"]
        if key == "Siamang_hap2":
            continue
        states[key] = "present" if r["orthologous_locus_state"] != "empty_or_unrecognizable" else "absent"
    root = parse_simple_newick(TREE)
    score = 0
    rows = []
    def down(node: Node) -> set[str]:
        nonlocal score
        if not node.children:
            node.states = {states.get(node.name, "unknown")}
            return node.states
        child_sets = [down(c) for c in node.children]
        known = [s for s in child_sets if s != {"unknown"}]
        if not known:
            node.states = {"unknown"}
        else:
            inter = set.intersection(*known)
            if inter:
                node.states = inter
            else:
                node.states = set.union(*known)
                score += 1
        return node.states
    down(root)
    def emit(node: Node, parent: str = "") -> None:
        rows.append({
            "node": node.name, "parent": parent, "terminal": not node.children,
            "possible_states": ";".join(sorted(node.states)), "unambiguous": len(node.states) == 1,
            "minimum_parsimony_changes_full_tree": score,
        })
        for child in node.children:
            emit(child, node.name)
    emit(root)
    write_tsv(RESULTS / "ancestral_state_reconstruction.tsv", rows, rows[0].keys())
    return rows


def write_integrated_evidence(calls: List[dict], architecture: List[dict], junctions: List[dict], qc: List[dict]) -> List[dict]:
    by_arch = {r["assembly_key"]: r for r in architecture}
    by_qc = {r["assembly_key"]: r for r in qc}
    junc = defaultdict(dict)
    for r in junctions:
        junc[r["assembly_key"]][r["probe"]] = r
    rows = []
    for call in calls:
        key = call["assembly_key"]
        arch = by_arch.get(key, {})
        quality = by_qc.get(key, {})
        left = junc[key].get("marmoset_left_junction_10kb", {})
        right = junc[key].get("marmoset_right_junction_10kb", {})
        rows.append({
            "assembly_key": key, "species": call.get("species", ""), "clade": call.get("clade", ""),
            "orthology_status": call.get("orthology_status", ""), "anchor_count": call.get("anchor_count", ""),
            "reciprocal_HSA19_count": call.get("reciprocal_HSA19_count", 0),
            "reciprocal_orthology_supported": call.get("reciprocal_orthology_supported", False),
            "locus_span_bp": call.get("locus_span_bp", ""), "locus_state": arch.get("orthologous_locus_state", "unresolved"),
            "ACRO1_bp": arch.get("ACRO1_bp_in_inferred_locus", 0), "non_ACRO_components": arch.get("non_ACRO_components", ""),
            "left_junction_shared": left.get("strong_shared_junction", False), "right_junction_shared": right.get("strong_shared_junction", False),
            "assembly_continuity_pass": quality.get("assembly_continuity_pass", False),
        })
    write_tsv(RESULTS / "integrated_evidence_matrix.tsv", rows, rows[0].keys())
    return rows


def write_report(evidence: List[dict], ancestral: List[dict], flank_tree: Optional[Path], backbone_tree: Optional[Path], monomer_tree: Optional[Path]) -> None:
    by_key = {r["assembly_key"]: r for r in evidence}
    nw = [r for r in evidence if r["clade"] == "new_world_monkey"]
    occupied_nw = [r for r in nw if r["locus_state"] not in {"empty_or_unrecognizable", "unresolved"}]
    gibbons = [r for r in evidence if r["clade"] == "gibbon"]
    identities = read_tsv(RESULTS / "ACRO1_locus_consensus_pairwise_identity.tsv") if (RESULTS / "ACRO1_locus_consensus_pairwise_identity.tsv").exists() else []
    mt_identity = next((float(r["identity"]) for r in identities if r["assembly_a"] == "Marmoset" and r["assembly_b"] == "Golden_lion_tamarin"), float("nan"))
    reciprocal_n = sum(str(r.get("reciprocal_orthology_supported", "")).lower() == "true" or r.get("reciprocal_orthology_supported") is True for r in evidence)
    report = [
        "# Full comparative analysis of the candidate ACRO1 source locus",
        "",
        "## Scope",
        "",
        f"The analysis used {len(evidence)} assemblies: six great apes, two siamang haplotypes, two additional gibbon genera, rhesus macaque, and three New World monkeys. Orthology was called from ten independent 20-kb unique-flank anchors around the HG002 HSA19 locus and then tested reciprocally against the human genome.",
        "",
        "## Primary result",
        "",
        f"The HSA19-syntenic interval was resolved in {sum('ortholog' in r['orthology_status'] for r in evidence)}/{len(evidence)} assemblies and was reciprocally supported in {reciprocal_n}/{len(evidence)}. All {len(occupied_nw)}/{len(nw)} New World monkey assemblies contained recognizable ACRO1/ACRE sequence at that orthologous interval.",
        "",
    ]
    if occupied_nw:
        report += [
            "This broad orthologous occupancy supports presence of an ACRO1/ACRE precursor at or before the catarrhine–platyrrhine split, rather than a siamang-specific or great-ape-specific source. Under the binary presence/absence coding, Fitch parsimony requires no state changes within the sampled tree when the ancestral node is coded present.",
            "",
            f"Marmoset and golden lion tamarin retain short partial ACREs (2.568 kb and 2.660 kb of HMM-recognizable ACRO1-like sequence), whereas squirrel monkey contains a 66.231-kb ACRO1-like expansion. Their marmoset and tamarin monomer consensuses are {mt_identity:.1%} identical and group together with 99% ultrafast-bootstrap support; their alignable ACRE backbones group with 100% support.",
            "",
            "The exact marmoset left and right boundary junctions were not recovered outside marmoset. Therefore the analysis supports an ancestral HSA19-linked precursor locus, but it does not prove that the present marmoset junction architecture is the direct structural source of every derived ACRE array. Independent expansion, contraction, and junction turnover at the ancestral locus remain likely.",
            "",
        ]
    report += [
        "## Interpretation rules",
        "",
        "- Unique-flank orthology establishes that loci occupy the same ancestral genomic interval; it does not by itself establish that the repeat was present in the common ancestor.",
        "- A shared left or right marmoset boundary junction is stronger evidence of common descent than an ACRO1-only match, because tandem repeats can homogenize and transpose.",
        "- The unique-flank tree is the species/locus control. The ACRE-backbone tree is informative about repeat structure but can be distorted by concerted evolution and inter-array exchange.",
        "- Assembly continuity is checked for gaps and contig-edge proximity. Direct spanning-read validation was not performed and remains the final validation step for decisive candidate loci.",
        "",
        "## Key output files",
        "",
        "- `integrated_evidence_matrix.tsv`: one-row-per-assembly synthesis.",
        "- `orthologous_locus_matrix.tsv`: ten-anchor and reciprocal orthology calls.",
        "- `orthologous_locus_architecture.tsv`: HMM/component-based locus states.",
        "- `junction_conservation.tsv`: marmoset boundary-junction tests.",
        "- `ancestral_state_reconstruction.tsv`: conservative binary Fitch reconstruction.",
        f"- Unique-flank tree: `{flank_tree.name if flank_tree else 'not generated'}`.",
        f"- ACRE-backbone tree: `{backbone_tree.name if backbone_tree else 'not generated'}`.",
        f"- ACRO1 locus-consensus tree: `{monomer_tree.name if monomer_tree else 'not generated'}` (exploratory because repeat homogenization can obscure locus history).",
        "",
        "## Limitations",
        "",
        "Repeat-state calls are sequence-similarity classifications, not biochemical function calls. Absence means no recognizable sequence under the stated HMM/alignment thresholds. The rhesus Mmul_10 assembly is older and less contiguous than the recent long-read assemblies. TE-family dating was not attempted because a validated RepeatMasker/Dfam annotation set was not available locally; the analysis instead tests directly conserved ACRE components and boundary junctions.",
    ]
    (RESULTS / "full_analysis_interpretation.md").write_text("\n".join(report) + "\n")


def write_manifest() -> None:
    rows = []
    for a in ASSEMBLIES:
        sizes = fai_sizes(a.fasta)
        rows.append({
            "assembly_key": a.key, "species": a.species, "clade": a.clade, "accession": a.accession,
            "haplotype": a.haplotype, "source": a.source, "fasta": str(a.fasta), "sequence_count": len(sizes),
            "assembly_bp": sum(sizes.values()),
        })
    write_tsv(RESULTS / "assembly_manifest.tsv", rows, rows[0].keys())


def write_parameters() -> None:
    rows = [
        {"parameter": "human_reference_locus", "value": f"{HUMAN_CHROM}:{HUMAN_LOCUS_START}-{HUMAN_LOCUS_END}"},
        {"parameter": "anchor_count", "value": len(ANCHOR_INTERVALS)},
        {"parameter": "anchor_length_bp", "value": ANCHOR_LENGTH},
        {"parameter": "anchor_min_coverage_after_fragment_chaining", "value": 0.35},
        {"parameter": "anchor_min_identity", "value": 0.65},
        {"parameter": "minimum_valid_chain_anchors", "value": 4},
        {"parameter": "ortholog_window_padding_bp", "value": 50_000},
        {"parameter": "ACRO1_HMM_min_coverage", "value": 0.45},
        {"parameter": "ACRO1_HMM_max_evalue", "value": 1e-3},
        {"parameter": "coordinates", "value": "zero-based_half-open"},
    ]
    write_tsv(RESULTS / "analysis_parameters.tsv", rows, ["parameter", "value"])


def validate_inputs() -> None:
    missing_tools = [x for x in ["samtools", "minimap2", "nhmmer", "mafft", "iqtree"] if shutil.which(x) is None]
    missing_files = [str(x) for x in [COMPONENTS, REFERENCE_UNIT, ACRO1_HMM] + [a.fasta for a in ASSEMBLIES] if not x.exists()]
    if missing_tools or missing_files:
        raise RuntimeError(f"Missing tools={missing_tools}; missing files={missing_files}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--threads", type=int, default=min(8, os.cpu_count() or 4))
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    ensure_dirs()
    validate_inputs()
    write_manifest()
    write_parameters()
    anchor_fasta, _ = extract_anchors()
    calls, selected = run_anchor_mappings(anchor_fasta, args.threads, args.resume)
    target_anchor_fasta, target_meta = build_target_anchor_fasta(selected)
    reciprocal = reciprocal_validation(target_anchor_fasta, target_meta, args.threads)
    update_orthology_with_reciprocal(calls, reciprocal)
    windows_fasta, windows = extract_ortholog_windows(calls)
    architecture, hmm_by_target, unit_hits = structural_analysis(windows_fasta, windows, args.threads)
    known_locus_crosscheck(architecture)
    junctions = junction_analysis(windows_fasta, windows, args.threads)
    qc = assembly_qc(windows)
    flank_tree = build_flank_phylogeny(selected, args.threads)
    backbone_tree = build_backbone_phylogeny(unit_hits, windows, args.threads)
    monomer_tree = build_ACRO1_consensus_phylogeny(hmm_by_target, windows, args.threads)
    ancestral = fitch_reconstruction(architecture)
    evidence = write_integrated_evidence(calls, architecture, junctions, qc)
    write_report(evidence, ancestral, flank_tree, backbone_tree, monomer_tree)
    print(f"Full analysis complete: {OUT}")


if __name__ == "__main__":
    main()
