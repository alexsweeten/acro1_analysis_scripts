#!/usr/bin/env python3
"""Build a compact synthetic FASTA/BED containing all six primate ACRE arrays."""

from __future__ import annotations

import csv
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
OUT = Path(__file__).resolve().parent

SAMPLES = [
    (
        "hg002",
        "human",
        "HG002v1.1",
        ROOT / "primate_genomes/hg002v1.1.fasta.gz",
        ROOT / "primate_acro_bedfiles/HG002/HG002_acro.named.bed",
    ),
    (
        "mPanPan1",
        "bonobo",
        "mPanPan1.dip.cur.20231122",
        ROOT / "primate_genomes/mPanPan1.dip.cur.20231122.fasta.gz",
        ROOT / "primate_acro_bedfiles/mPanPan1/mPanPan1_acro.named.bed",
    ),
    (
        "mPanTro3",
        "chimpanzee",
        "mPanTro3.dip.cur.20231122",
        ROOT / "primate_genomes/mPanTro3.dip.cur.20231122.fasta.gz",
        ROOT / "primate_acro_bedfiles/mPanTro3/mPanTro3_acro.named.bed",
    ),
    (
        "mGorGor1",
        "gorilla",
        "mGorGor1.dip.cur.20231122",
        ROOT / "primate_genomes/mGorGor1.dip.cur.20231122.fasta.gz",
        ROOT / "primate_acro_bedfiles/mGorGor1/mGorGor1_acro.named.bed",
    ),
    (
        "mPonPyg2",
        "bornean_orangutan",
        "mPonPyg2.dip.cur.20231122",
        ROOT / "primate_genomes/mPonPyg2.dip.cur.20231122.fasta.gz",
        ROOT / "primate_acro_bedfiles/mPonPyg2/mPonPyg2_acro.named.bed",
    ),
    (
        "mPonAbe1",
        "sumatran_orangutan",
        "mPonAbe1.dip.cur.20231205",
        ROOT / "primate_genomes/mPonAbe1.dip.cur.20231205.fasta.gz",
        ROOT / "primate_acro_bedfiles/mPonAbe1/mPonAbe1_acro.named.bed",
    ),
]


def read_bed(path: Path) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    with path.open() as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError(f"{path}:{line_number}: fewer than three BED fields")
            start, end = int(fields[1]), int(fields[2])
            if start < 0 or end <= start:
                raise ValueError(f"{path}:{line_number}: invalid interval")
            records.append(
                {
                    "chrom": fields[0],
                    "start": start,
                    "end": end,
                    "name": fields[3] if len(fields) >= 4 else f"bed_row_{line_number}",
                }
            )
    return records


def fasta_records(path: Path):
    name: str | None = None
    sequence: list[str] = []
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(sequence)
                name = line[1:].strip().split()[0]
                sequence = []
            else:
                sequence.append(line.strip())
    if name is not None:
        yield name, "".join(sequence)


def main() -> None:
    tmp = OUT / "extracted_arrays"
    tmp.mkdir(exist_ok=True)
    fasta_path = OUT / "six_primate_arrays.fasta"
    bed_path = OUT / "six_primate_arrays.bed"
    mapping_path = OUT / "six_primate_array_mapping.tsv"
    summary_path = OUT / "six_primate_input_summary.tsv"

    mapping_fields = [
        "input_index",
        "synthetic_chrom",
        "species_id",
        "species_name",
        "assembly",
        "original_chrom",
        "original_start",
        "original_end",
        "input_length",
        "original_bed_name",
        "source_bed",
        "source_fasta",
    ]
    input_index = 0
    summaries: list[tuple[str, str, int, int]] = []

    with (
        fasta_path.open("w") as fasta_out,
        bed_path.open("w") as bed_out,
        mapping_path.open("w", newline="") as mapping_out,
    ):
        writer = csv.DictWriter(mapping_out, fieldnames=mapping_fields, delimiter="\t")
        writer.writeheader()

        for species_id, species_name, assembly, genome, bed in SAMPLES:
            for required in (genome, Path(str(genome) + ".fai"), Path(str(genome) + ".gzi"), bed):
                if not required.exists():
                    raise FileNotFoundError(required)
            records = read_bed(bed)
            regions_path = tmp / f"{species_id}.regions.txt"
            extracted_path = tmp / f"{species_id}.fasta"
            with regions_path.open("w") as regions:
                for record in records:
                    # samtools faidx regions are one-based and end-inclusive.
                    regions.write(
                        f"{record['chrom']}:{int(record['start']) + 1}-{record['end']}\n"
                    )
            with extracted_path.open("w") as extracted:
                subprocess.run(
                    ["samtools", "faidx", "-r", str(regions_path), str(genome)],
                    stdout=extracted,
                    check=True,
                )

            extracted_records = list(fasta_records(extracted_path))
            if len(extracted_records) != len(records):
                raise RuntimeError(
                    f"{species_id}: requested {len(records)} arrays but extracted "
                    f"{len(extracted_records)}"
                )
            species_bp = 0
            for record, (_, sequence) in zip(records, extracted_records):
                expected = int(record["end"]) - int(record["start"])
                if len(sequence) != expected:
                    raise RuntimeError(
                        f"{species_id} {record['chrom']}:{record['start']}-{record['end']}: "
                        f"expected {expected} bp, received {len(sequence)}"
                    )
                input_index += 1
                species_bp += len(sequence)
                synthetic = f"PF{input_index:06d}"
                source_name = (
                    f"{species_id}|{record['chrom']}:{record['start']}-{record['end']}|"
                    f"{record['name']}"
                )
                fasta_out.write(f">{synthetic}\n")
                for offset in range(0, len(sequence), 80):
                    fasta_out.write(sequence[offset : offset + 80] + "\n")
                bed_out.write(f"{synthetic}\t0\t{len(sequence)}\t{source_name}\n")
                writer.writerow(
                    {
                        "input_index": input_index,
                        "synthetic_chrom": synthetic,
                        "species_id": species_id,
                        "species_name": species_name,
                        "assembly": assembly,
                        "original_chrom": record["chrom"],
                        "original_start": record["start"],
                        "original_end": record["end"],
                        "input_length": len(sequence),
                        "original_bed_name": record["name"],
                        "source_bed": bed.relative_to(ROOT),
                        "source_fasta": genome.relative_to(ROOT),
                    }
                )
            summaries.append((species_id, species_name, len(records), species_bp))

    with summary_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["species_id", "species_name", "input_arrays", "input_bp"])
        writer.writerows(summaries)
        writer.writerow(
            ["TOTAL", "all_species", sum(x[2] for x in summaries), sum(x[3] for x in summaries)]
        )

    subprocess.run(["samtools", "faidx", str(fasta_path)], check=True)
    print(f"Wrote {input_index} arrays and {sum(x[3] for x in summaries):,} bp")


if __name__ == "__main__":
    main()
