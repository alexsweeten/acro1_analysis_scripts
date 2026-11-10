#!/usr/bin/env python3
"""De novo test for unusually long ORFs in primate ACRO1 monomer consensuses.

The analysis uses only CenSat ACRO1 interval labels and the underlying assembly
sequence. Tandem Repeats Finder is used to recover ~147-bp monomers de novo.
For every CenSat interval, the highest-scoring ACRO1-like TRF call supplies an
array-level monomer consensus and an observed genomic tract. The primary null
shuffles each monomer consensus while preserving its exact dinucleotide counts,
then reconstructs a ten-monomer tandem repeat and repeats the six-frame ORF scan.
"""

from __future__ import annotations

import argparse
import csv
import math
import random
import shutil
import subprocess
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PRIMATE = ROOT.parent / "ACREAGE" / "primate_genomes"


@dataclass(frozen=True)
class Assembly:
    key: str
    species: str
    assembly: str
    bed: Path
    fasta: Path | None
    fai: Path


ASSEMBLIES = [
    Assembly("CHM13", "Human", "T2T-CHM13v2.0", ROOT / "chm13" / "chm13v2.0.cenSatv2.1.bed", PRIMATE / "chm13v2.0.fa.gz", PRIMATE / "chm13v2.0.fa.gz.fai"),
    Assembly("HG002", "Human", "HG002v1.1", PRIMATE / "hg002v1.1.cenSatv2.0.bed", PRIMATE / "hg002v1.1.fasta.gz", PRIMATE / "hg002v1.1.fasta.gz.fai"),
    Assembly("Gorilla", "Gorilla", "mGorGor1", PRIMATE / "GorGor.cenSatv2.0.bed", PRIMATE / "mGorGor1.dip.cur.20231122.fasta.gz", PRIMATE / "mGorGor1.dip.cur.20231122.fasta.gz.fai"),
    Assembly("Chimpanzee", "Chimpanzee", "mPanTro3", PRIMATE / "PanTro.cenSatv2.0.bed", PRIMATE / "mPanTro3.dip.cur.20231122.fasta.gz", PRIMATE / "mPanTro3.dip.cur.20231122.fasta.gz.fai"),
    Assembly("Bonobo", "Bonobo", "mPanPan1", PRIMATE / "PanPan.cenSatv2.0.bed", PRIMATE / "mPanPan1.dip.cur.20231122.fasta.gz", PRIMATE / "mPanPan1.dip.cur.20231122.fasta.gz.fai"),
    Assembly("Sumatran_orangutan", "Sumatran orangutan", "mPonAbe1", PRIMATE / "PonAbe.cenSatv2.0.bed", PRIMATE / "mPonAbe1.dip.cur.20231205.fasta.gz", PRIMATE / "mPonAbe1.dip.cur.20231205.fasta.gz.fai"),
    Assembly("Bornean_orangutan", "Bornean orangutan", "mPonPyg2", PRIMATE / "PonPyg.cenSatv2.0.bed", PRIMATE / "mPonPyg2.dip.cur.20231122.fasta.gz", PRIMATE / "mPonPyg2.dip.cur.20231122.fasta.gz.fai"),
    Assembly("Gibbon", "Siamang gibbon", "mSymSyn1", PRIMATE / "SymSyn.cenSatv2.0.bed", None, PRIMATE / "mSymSyn1.analysis-dip.20240514.fasta.gz.fai"),
    Assembly("Marmoset", "Common marmoset", "calJac240_pri.v1.1", ROOT / "primate_acro_hsat1a_enrichment" / "eight_species_gibbon_marmoset" / "inputs" / "marmoset.cenSatv2.0.bed", None, PRIMATE / "calJac240_pri.v1.1.fa.fai.txt"),
]


CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
RC_TABLE = str.maketrans("ACGTN", "TGCAN")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(RC_TABLE)[::-1]


def translations(sequence: str) -> list[tuple[str, int, str]]:
    result = []
    for strand, oriented in [("+", sequence), ("-", reverse_complement(sequence))]:
        for frame in range(3):
            aa = "".join(
                CODON_TABLE.get(oriented[index : index + 3], "X")
                for index in range(frame, len(oriented) - 2, 3)
            )
            result.append((strand, frame, aa))
    return result


def longest_stop_free(aa: str) -> int:
    return max((len(part) for part in aa.split("*")), default=0)


def longest_atg_open(aa: str) -> int:
    best = 0
    active: list[int] = []
    for position, residue in enumerate(aa):
        if residue == "M":
            active.append(position)
        if residue == "*":
            if active:
                best = max(best, position - min(active))
            active = []
    if active:
        best = max(best, len(aa) - min(active))
    return best


def scan_orfs(sequence: str) -> dict[str, object]:
    scans = translations(sequence)
    stop_counts = [aa.count("*") for _, _, aa in scans]
    stopfree = [longest_stop_free(aa) for _, _, aa in scans]
    atg_open = [longest_atg_open(aa) for _, _, aa in scans]
    best_index = max(range(6), key=lambda index: (atg_open[index], stopfree[index], -stop_counts[index]))
    strand, frame, aa = scans[best_index]
    best_segment = max(aa.split("*"), key=len, default="")
    pst_fraction = (
        sum(residue in "PST" for residue in best_segment) / len(best_segment)
        if best_segment else math.nan
    )
    return {
        "max_stopfree_aa": max(stopfree),
        "max_atg_open_aa": max(atg_open),
        "zero_stop_frames": sum(count == 0 for count in stop_counts),
        "minimum_stop_count": min(stop_counts),
        "best_strand": strand,
        "best_frame": frame,
        "best_pst_fraction": pst_fraction,
    }


def scan_tandem_orf_lengths(monomer: str, copies: int = 10) -> dict[str, int]:
    """Fast six-frame ORF scan for an exact tandem reconstruction.

    When monomer length is divisible by three, the codon stream is periodic and
    can be evaluated from one monomer rather than materializing all ten copies.
    Other lengths fall back to the general scanner.
    """
    if len(monomer) % 3:
        result = scan_orfs(monomer * copies)
        return {
            "max_stopfree_aa": int(result["max_stopfree_aa"]),
            "max_atg_open_aa": int(result["max_atg_open_aa"]),
            "zero_stop_frames": int(result["zero_stop_frames"]),
        }
    period_codons = len(monomer) // 3
    best_stopfree = 0
    best_atg = 0
    zero_frames = 0
    for oriented in (monomer, reverse_complement(monomer)):
        circular = oriented + oriented[:2]
        for frame in range(3):
            codons = [circular[frame + 3 * index : frame + 3 * index + 3] for index in range(period_codons)]
            stops = [index for index, codon in enumerate(codons) if codon in {"TAA", "TAG", "TGA"}]
            starts = [index for index, codon in enumerate(codons) if codon == "ATG"]
            total_codons = (len(monomer) * copies - frame) // 3
            if not stops:
                zero_frames += 1
                best_stopfree = max(best_stopfree, total_codons)
                if starts:
                    best_atg = max(best_atg, total_codons - min(starts))
                continue
            cyclic_stopfree = max(
                (right - left - 1) % period_codons
                for left, right in zip(stops, stops[1:] + [stops[0]])
            )
            best_stopfree = max(best_stopfree, cyclic_stopfree)
            for start in starts:
                next_stop_distance = min((stop - start) % period_codons for stop in stops if stop != start)
                best_atg = max(best_atg, next_stop_distance)
    return {
        "max_stopfree_aa": best_stopfree,
        "max_atg_open_aa": best_atg,
        "zero_stop_frames": zero_frames,
    }


def dinucleotide_shuffle(sequence: str, rng: random.Random) -> str:
    """Return an Eulerian shuffle with the exact original dinucleotide counts."""
    if len(sequence) < 2:
        return sequence
    adjacency: dict[str, list[str]] = defaultdict(list)
    for left, right in zip(sequence, sequence[1:]):
        adjacency[left].append(right)
    for edges in adjacency.values():
        rng.shuffle(edges)
    stack = [sequence[0]]
    trail: list[str] = []
    while stack:
        node = stack[-1]
        if adjacency[node]:
            stack.append(adjacency[node].pop())
        else:
            trail.append(stack.pop())
    shuffled = "".join(reversed(trail))
    if len(shuffled) != len(sequence):
        raise RuntimeError("Dinucleotide shuffle did not preserve length")
    return shuffled


def read_fai(path: Path) -> dict[str, int]:
    sizes = {}
    if not path.exists():
        return sizes
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2 and fields[1].isdigit():
                sizes[fields[0]] = int(fields[1])
    return sizes


def read_acro_bed(path: Path) -> list[tuple[str, int, int, str]]:
    records = []
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith(("track", "browser", "#")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 4 and "ACRO1" in fields[3]:
                records.append((fields[0], int(fields[1]), int(fields[2]), fields[3]))
    return records


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    name: str | None = None
    parts: list[str] = []
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(parts).upper()
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
    if name is not None:
        records[name] = "".join(parts).upper()
    return records


def parse_trf(path: Path) -> dict[str, list[dict[str, object]]]:
    calls: dict[str, list[dict[str, object]]] = defaultdict(list)
    current: str | None = None
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if line.startswith("Sequence: "):
                current = line.split(maxsplit=1)[1]
                continue
            if not current or not line or not line[0].isdigit():
                continue
            fields = line.split()
            if len(fields) < 15:
                continue
            call = {
                "start": int(fields[0]),
                "end": int(fields[1]),
                "period": int(fields[2]),
                "copies": float(fields[3]),
                "consensus_size": int(fields[4]),
                "pct_match": int(fields[5]),
                "pct_indel": int(fields[6]),
                "score": int(fields[7]),
                "entropy": float(fields[12]),
                "consensus": fields[13].upper(),
            }
            if (
                135 <= call["period"] <= 159
                and call["copies"] >= 1.5
                and call["end"] - call["start"] + 1 >= 200
                and call["pct_match"] >= 65
                and call["entropy"] >= 1.5
                and set(call["consensus"]) <= set("ACGT")
            ):
                calls[current].append(call)
    return calls


def choose_nonoverlapping(calls: list[dict[str, object]]) -> list[dict[str, object]]:
    retained: list[dict[str, object]] = []
    for call in sorted(calls, key=lambda item: (-int(item["score"]), -int(item["end"]) + int(item["start"]))):
        start, end = int(call["start"]), int(call["end"])
        span = end - start + 1
        duplicate = False
        for old in retained:
            old_start, old_end = int(old["start"]), int(old["end"])
            overlap = max(0, min(end, old_end) - max(start, old_start) + 1)
            if overlap >= 0.5 * min(span, old_end - old_start + 1):
                duplicate = True
                break
        if not duplicate:
            retained.append(call)
    return sorted(retained, key=lambda item: int(item["start"]))


def empirical_p(observed: float, null: list[float]) -> float:
    return (1 + sum(value >= observed for value in null)) / (1 + len(null))


def fmt(value: object) -> object:
    if isinstance(value, float):
        if math.isnan(value):
            return "NA"
        return f"{value:.10g}"
    return value


parser = argparse.ArgumentParser()
parser.add_argument("--outdir", type=Path, default=ROOT / "primate_acro1_orf_test")
parser.add_argument("--permutations", type=int, default=999)
parser.add_argument("--seed", type=int, default=20261022)
parser.add_argument("--no-reuse-intermediate", action="store_true")
args = parser.parse_args()

if shutil.which("bedtools") is None or shutil.which("trf") is None:
    raise SystemExit("bedtools and trf are required")

args.outdir.mkdir(parents=True, exist_ok=True)
(args.outdir / "intermediate").mkdir(exist_ok=True)

audit_rows = []
array_rows: list[dict[str, object]] = []
tract_rows: list[dict[str, object]] = []
species_null_atg: dict[str, list[float]] = {}
species_null_stopfree: dict[str, list[float]] = {}
species_counts: Counter[str] = Counter()

for assembly_index, assembly in enumerate(ASSEMBLIES):
    records = read_acro_bed(assembly.bed) if assembly.bed.exists() else []
    sizes = read_fai(assembly.fai)
    fasta_available = assembly.fasta is not None and assembly.fasta.exists()
    coordinate_ok = bool(records) and all(chrom in sizes and end <= sizes[chrom] for chrom, start, end, label in records)
    status = "ready" if fasta_available and coordinate_ok else (
        "no_matching_fasta" if not fasta_available else "coordinate_mismatch"
    )
    audit_rows.append({
        "key": assembly.key,
        "species": assembly.species,
        "assembly": assembly.assembly,
        "bed": str(assembly.bed),
        "fasta": str(assembly.fasta) if assembly.fasta else "NA",
        "fai": str(assembly.fai),
        "acro_intervals": len(records),
        "acro_bp": sum(end - start for chrom, start, end, label in records),
        "coordinate_ok": coordinate_ok,
        "status": status,
    })
    if status != "ready":
        continue

    prefix = args.outdir / "intermediate" / assembly.key
    bed_path = prefix.with_suffix(".acro_arrays.bed")
    fasta_path = prefix.with_suffix(".acro_arrays.fa")
    if args.no_reuse_intermediate or not (bed_path.exists() and fasta_path.exists()):
        with bed_path.open("w") as handle:
            for index, (chrom, start, end, label) in enumerate(records, start=1):
                handle.write(f"{chrom}\t{start}\t{end}\t{assembly.key}_array{index:04d}\n")
        subprocess.run(
            ["bedtools", "getfasta", "-fi", str(assembly.fasta), "-bed", str(bed_path), "-nameOnly", "-fo", str(fasta_path)],
            check=True,
        )
    fasta_records = read_fasta(fasta_path)
    saved_dat = prefix.with_suffix(".trf.dat")
    if args.no_reuse_intermediate or not saved_dat.exists():
        with tempfile.TemporaryDirectory(prefix=f"acro1-trf-{assembly.key}-") as temp_name:
            temp = Path(temp_name)
            local_fasta = temp / fasta_path.name
            shutil.copyfile(fasta_path, local_fasta)
            subprocess.run(
                ["trf", local_fasta.name, "2", "7", "7", "80", "10", "30", "200", "-h", "-d"],
                check=True,
                cwd=temp,
                stdout=subprocess.DEVNULL,
            )
            dat_path = temp / f"{local_fasta.name}.2.7.7.80.10.30.200.dat"
            shutil.copyfile(dat_path, saved_dat)
    trf_calls = parse_trf(saved_dat)

    species_null_atg[assembly.key] = [0.0] * args.permutations
    species_null_stopfree[assembly.key] = [0.0] * args.permutations
    rng = random.Random(args.seed + assembly_index * 100003)

    for index, (chrom, start, end, label) in enumerate(records, start=1):
        array_id = f"{assembly.key}_array{index:04d}"
        sequence = fasta_records[array_id]
        retained = choose_nonoverlapping(trf_calls.get(array_id, []))
        for tract_index, call in enumerate(retained, start=1):
            tract = sequence[int(call["start"]) - 1 : int(call["end"])]
            tract_scan = scan_orfs(tract)
            tract_rows.append({
                "key": assembly.key,
                "species": assembly.species,
                "assembly": assembly.assembly,
                "array_id": array_id,
                "chrom": chrom,
                "array_start_bed": start,
                "array_end_bed": end,
                "tract_index": tract_index,
                "tract_start_1based_in_array": call["start"],
                "tract_end_1based_in_array": call["end"],
                "tract_length": len(tract),
                "trf_period": call["period"],
                "trf_copies": call["copies"],
                "trf_score": call["score"],
                "trf_pct_match": call["pct_match"],
                "actual_max_stopfree_aa": tract_scan["max_stopfree_aa"],
                "actual_max_atg_open_aa": tract_scan["max_atg_open_aa"],
                "actual_zero_stop_frames": tract_scan["zero_stop_frames"],
                "actual_best_strand": tract_scan["best_strand"],
                "actual_best_frame": tract_scan["best_frame"],
                "actual_best_pst_fraction": tract_scan["best_pst_fraction"],
            })
        if not retained:
            array_rows.append({
                "key": assembly.key,
                "species": assembly.species,
                "assembly": assembly.assembly,
                "array_id": array_id,
                "chrom": chrom,
                "array_start_bed": start,
                "array_end_bed": end,
                "array_length": end - start,
                "trf_acro1_tracts": 0,
                "status": "no_de_novo_147bp_call",
            })
            continue

        representative = max(retained, key=lambda item: (int(item["score"]), int(item["end"]) - int(item["start"])))
        consensus = str(representative["consensus"])
        synthetic = consensus * 10
        observed = scan_orfs(synthetic)
        actual_tract = sequence[int(representative["start"]) - 1 : int(representative["end"])]
        actual = scan_orfs(actual_tract)
        null_atg: list[float] = []
        null_stopfree: list[float] = []
        null_zero_frames: list[float] = []
        for permutation in range(args.permutations):
            shuffled = dinucleotide_shuffle(consensus, rng)
            null_scan = scan_tandem_orf_lengths(shuffled, 10)
            atg_value = float(null_scan["max_atg_open_aa"])
            stopfree_value = float(null_scan["max_stopfree_aa"])
            null_atg.append(atg_value)
            null_stopfree.append(stopfree_value)
            null_zero_frames.append(float(null_scan["zero_stop_frames"]))
            species_null_atg[assembly.key][permutation] += atg_value
            species_null_stopfree[assembly.key][permutation] += stopfree_value
        species_counts[assembly.key] += 1
        array_rows.append({
            "key": assembly.key,
            "species": assembly.species,
            "assembly": assembly.assembly,
            "array_id": array_id,
            "chrom": chrom,
            "array_start_bed": start,
            "array_end_bed": end,
            "array_length": end - start,
            "trf_acro1_tracts": len(retained),
            "status": "ok",
            "representative_period": representative["period"],
            "representative_copies": representative["copies"],
            "representative_score": representative["score"],
            "monomer_consensus": consensus,
            "consensus_length": len(consensus),
            "observed_max_stopfree_aa_10x": observed["max_stopfree_aa"],
            "observed_max_atg_open_aa_10x": observed["max_atg_open_aa"],
            "observed_zero_stop_frames_10x": observed["zero_stop_frames"],
            "observed_best_strand_10x": observed["best_strand"],
            "observed_best_frame_10x": observed["best_frame"],
            "observed_best_pst_fraction_10x": observed["best_pst_fraction"],
            "actual_tract_length": len(actual_tract),
            "actual_max_stopfree_aa": actual["max_stopfree_aa"],
            "actual_max_atg_open_aa": actual["max_atg_open_aa"],
            "actual_zero_stop_frames": actual["zero_stop_frames"],
            "actual_best_pst_fraction": actual["best_pst_fraction"],
            "null_mean_max_stopfree_aa_10x": sum(null_stopfree) / len(null_stopfree),
            "null_mean_max_atg_open_aa_10x": sum(null_atg) / len(null_atg),
            "p_max_stopfree_aa": empirical_p(float(observed["max_stopfree_aa"]), null_stopfree),
            "p_max_atg_open_aa": empirical_p(float(observed["max_atg_open_aa"]), null_atg),
            "p_zero_stop_frames": empirical_p(float(observed["zero_stop_frames"]), null_zero_frames),
        })


def write_dicts(path: Path, rows: list[dict[str, object]]) -> None:
    fields: list[str] = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: fmt(row.get(field, "NA")) for field in fields})


write_dicts(args.outdir / "input_audit.tsv", audit_rows)
write_dicts(args.outdir / "array_orf_results.tsv", array_rows)
write_dicts(args.outdir / "tract_orf_results.tsv", tract_rows)

summary_rows = []
for assembly in ASSEMBLIES:
    relevant = [row for row in array_rows if row["key"] == assembly.key]
    ok = [row for row in relevant if row.get("status") == "ok"]
    if not ok:
        summary_rows.append({
            "key": assembly.key,
            "species": assembly.species,
            "assembly": assembly.assembly,
            "arrays_tested": 0,
            "status": next((row["status"] for row in audit_rows if row["key"] == assembly.key), "not_tested"),
        })
        continue
    count = len(ok)
    observed_atg_mean = sum(float(row["observed_max_atg_open_aa_10x"]) for row in ok) / count
    observed_stopfree_mean = sum(float(row["observed_max_stopfree_aa_10x"]) for row in ok) / count
    null_atg_means = [value / count for value in species_null_atg[assembly.key]]
    null_stopfree_means = [value / count for value in species_null_stopfree[assembly.key]]
    actual_atg = sorted(float(row["actual_max_atg_open_aa"]) for row in ok)
    actual_stopfree = sorted(float(row["actual_max_stopfree_aa"]) for row in ok)
    summary_rows.append({
        "key": assembly.key,
        "species": assembly.species,
        "assembly": assembly.assembly,
        "arrays_total": len(relevant),
        "arrays_tested": count,
        "arrays_without_147bp_call": len(relevant) - count,
        "total_de_novo_tracts": sum(int(row["trf_acro1_tracts"]) for row in ok),
        "arrays_with_at_least_one_zero_stop_frame_10x": sum(int(row["observed_zero_stop_frames_10x"]) > 0 for row in ok),
        "fraction_arrays_with_zero_stop_frame_10x": sum(int(row["observed_zero_stop_frames_10x"]) > 0 for row in ok) / count,
        "mean_observed_atg_open_aa_10x": observed_atg_mean,
        "mean_null_atg_open_aa_10x": sum(null_atg_means) / len(null_atg_means),
        "aggregate_p_atg_open": empirical_p(observed_atg_mean, null_atg_means),
        "mean_observed_stopfree_aa_10x": observed_stopfree_mean,
        "mean_null_stopfree_aa_10x": sum(null_stopfree_means) / len(null_stopfree_means),
        "aggregate_p_stopfree": empirical_p(observed_stopfree_mean, null_stopfree_means),
        "median_actual_tract_atg_open_aa": actual_atg[len(actual_atg) // 2],
        "median_actual_tract_stopfree_aa": actual_stopfree[len(actual_stopfree) // 2],
        "status": "ok",
    })

write_dicts(args.outdir / "species_summary.tsv", summary_rows)

with (args.outdir / "analysis_parameters.tsv").open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["parameter", "value"])
    writer.writerows([
        ["analysis", "de novo ACRO1 monomer ORF test"],
        ["ACRO_interval_definition", "CenSat annotation label contains ACRO1"],
        ["monomer_discovery", "TRF 2 7 7 80 10 30 200"],
        ["accepted_period_bp", "135-159"],
        ["accepted_minimum_copy_number", "1.5"],
        ["representative_monomer", "highest-scoring accepted TRF call per CenSat interval"],
        ["ORF_scan", "all six frames; ten tandem copies of de novo monomer consensus"],
        ["primary_ORF_metric", "longest ATG-initiated segment before stop or sequence boundary"],
        ["secondary_ORF_metric", "longest stop-free amino-acid segment"],
        ["null", "exact dinucleotide-preserving shuffle of each monomer consensus, then ten-copy tandem reconstruction"],
        ["permutations", args.permutations],
        ["seed", args.seed],
        ["empirical_p", "(1 + null >= observed) / (1 + permutations)"],
    ])

print(args.outdir)
