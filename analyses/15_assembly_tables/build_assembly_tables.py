#!/usr/bin/env python3
"""Build manuscript assembly tables from the verified local manifests."""

from __future__ import annotations

import csv
import os
from pathlib import Path


ROOT = Path(os.environ.get("ACRO1_WORKSPACE", Path(__file__).resolve().parents[2]))
OUT = Path(os.environ.get("ACRO1_ASSEMBLY_TABLE_OUTPUT", Path(__file__).resolve().parent))
HPRC_MANIFEST = ROOT / "NTRIAGE_OCT26/manifest.validated.tsv"
HPRC_INDEX = ROOT / (
    "NTRIAGE_OCT26/manuscript_output/hprc_telomeric_chromosome_assignment/"
    "inputs/assemblies_release2_v1.0.index.csv"
)


def ncbi(accession: str) -> str:
    return f"https://www.ncbi.nlm.nih.gov/datasets/genome/{accession}/"


comparative_rows = [
    {
        "dataset_group": "Comparative great apes; HSA19 orthology",
        "species": "Homo sapiens",
        "assembly": "T2T-HG002v1.1",
        "haplotype_or_phase": "Diploid; maternal haplotype used to define the HSA19 locus",
        "assembly_accession": "GCA_018852605.3 (paternal); GCA_018852615.3 (maternal)",
        "public_data_location": "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/hg002v1.1.fasta.gz",
        "local_fasta": "data/primate_genomes/hg002v1.1.fasta.gz",
        "local_annotation": "data/primate_genomes/hg002v1.1.cenSatv2.0.bed",
        "analysis_role": "Great-ape comparison, HSA19 locus definition, and human reference-array analyses",
    },
    {
        "dataset_group": "Comparative great apes; HSA19 orthology",
        "species": "Pan troglodytes",
        "assembly": "mPanTro3 (curated diploid release 2023-11-22)",
        "haplotype_or_phase": "Primary and alternate pseudohaplotypes",
        "assembly_accession": "GCA_028858775.2 (primary); GCA_028858805.2 (alternate)",
        "public_data_location": "https://genomeark.s3.amazonaws.com/species/Pan_troglodytes/mPanTro3/assembly_curated/mPanTro3.dip.cur.20231122.fasta.gz",
        "local_fasta": "data/primate_genomes/mPanTro3.dip.cur.20231122.fasta.gz",
        "local_annotation": "data/primate_genomes/PanTro.cenSatv2.0.bed",
        "analysis_role": "Great-ape comparison and HSA19 orthology",
    },
    {
        "dataset_group": "Comparative great apes; HSA19 orthology",
        "species": "Pan paniscus",
        "assembly": "mPanPan1 (curated diploid release 2023-11-22)",
        "haplotype_or_phase": "Maternal and paternal haplotypes",
        "assembly_accession": "GCA_028858845.2 (maternal); GCA_028858825.2 (paternal)",
        "public_data_location": "https://genomeark.s3.amazonaws.com/species/Pan_paniscus/mPanPan1/assembly_curated/mPanPan1.dip.cur.20231122.fasta.gz",
        "local_fasta": "data/primate_genomes/mPanPan1.dip.cur.20231122.fasta.gz",
        "local_annotation": "data/primate_genomes/PanPan.cenSatv2.0.bed",
        "analysis_role": "Great-ape comparison and HSA19 orthology",
    },
    {
        "dataset_group": "Comparative great apes; HSA19 orthology",
        "species": "Gorilla gorilla gorilla",
        "assembly": "mGorGor1 (curated diploid release 2023-11-22)",
        "haplotype_or_phase": "Maternal and paternal haplotypes",
        "assembly_accession": "GCA_028885495.2 (maternal); GCA_028885475.2 (paternal)",
        "public_data_location": "https://genomeark.s3.amazonaws.com/species/Gorilla_gorilla/mGorGor1/assembly_curated/mGorGor1.dip.cur.20231122.fasta.gz",
        "local_fasta": "data/primate_genomes/mGorGor1.dip.cur.20231122.fasta.gz",
        "local_annotation": "data/primate_genomes/GorGor.cenSatv2.0.bed",
        "analysis_role": "Great-ape comparison, gorilla phylogeny, and HSA19 orthology",
    },
    {
        "dataset_group": "Comparative great apes; HSA19 orthology",
        "species": "Pongo abelii",
        "assembly": "mPonAbe1 (curated diploid release 2023-12-05)",
        "haplotype_or_phase": "Primary and alternate pseudohaplotypes",
        "assembly_accession": "GCA_028885655.2 (primary); GCA_028885685.2 (alternate)",
        "public_data_location": "https://genomeark.s3.amazonaws.com/species/Pongo_abelii/mPonAbe1/assembly_curated/mPonAbe1.dip.cur.20231205.fasta.gz",
        "local_fasta": "data/primate_genomes/mPonAbe1.dip.cur.20231205.fasta.gz",
        "local_annotation": "data/primate_genomes/PonAbe.cenSatv2.0.bed",
        "analysis_role": "Great-ape comparison and HSA19 orthology",
    },
    {
        "dataset_group": "Comparative great apes; HSA19 orthology",
        "species": "Pongo pygmaeus",
        "assembly": "mPonPyg2 (curated diploid release 2023-11-22)",
        "haplotype_or_phase": "Primary and alternate pseudohaplotypes",
        "assembly_accession": "GCA_028885625.2 (primary); GCA_028885525.2 (alternate)",
        "public_data_location": "https://genomeark.s3.amazonaws.com/species/Pongo_pygmaeus/mPonPyg2/assembly_curated/mPonPyg2.dip.cur.20231122.fasta.gz",
        "local_fasta": "data/primate_genomes/mPonPyg2.dip.cur.20231122.fasta.gz",
        "local_annotation": "data/primate_genomes/PonPyg.cenSatv2.0.bed",
        "analysis_role": "Great-ape comparison and HSA19 orthology",
    },
    {
        "dataset_group": "HSA19 orthology",
        "species": "Symphalangus syndactylus",
        "assembly": "NHGRI_mSymSyn1-v2.1_pri",
        "haplotype_or_phase": "Primary",
        "assembly_accession": "GCA_028878055.3",
        "public_data_location": ncbi("GCA_028878055.3"),
        "local_fasta": str(ROOT / "phylogeny_origin_acro1_full/inputs/ncbi_downloads/core_outgroups/ncbi_dataset/data/GCA_028878055.3/GCA_028878055.3_NHGRI_mSymSyn1-v2.1_pri_genomic.fna"),
        "local_annotation": "Not required for the anchor-based HSA19 analysis",
        "analysis_role": "HSA19 orthology and non-acrocentric ACRO-locus analysis",
    },
    {
        "dataset_group": "HSA19 orthology",
        "species": "Symphalangus syndactylus",
        "assembly": "NHGRI_mSymSyn1-v2.0_alt",
        "haplotype_or_phase": "Alternate",
        "assembly_accession": "GCA_028878085.3",
        "public_data_location": ncbi("GCA_028878085.3"),
        "local_fasta": str(ROOT / "phylogeny_origin_acro1_full/inputs/ncbi_downloads/core_outgroups/ncbi_dataset/data/GCA_028878085.3/GCA_028878085.3_NHGRI_mSymSyn1-v2.0_alt_genomic.fna"),
        "local_annotation": "Not required for the anchor-based HSA19 analysis",
        "analysis_role": "HSA19 orthology and non-acrocentric ACRO-locus analysis",
    },
    {
        "dataset_group": "HSA19 orthology",
        "species": "Hoolock leuconedys",
        "assembly": "UCONN_HLeu1",
        "haplotype_or_phase": "Haploid representation",
        "assembly_accession": "GCA_047372625.1",
        "public_data_location": ncbi("GCA_047372625.1"),
        "local_fasta": str(ROOT / "phylogeny_origin_acro1_full/inputs/ncbi_downloads/expanded_individual/GCA_047372625.1_UCONN_HLeu1_genomic.fna"),
        "local_annotation": "Not required for the anchor-based HSA19 analysis",
        "analysis_role": "HSA19 orthology",
    },
    {
        "dataset_group": "HSA19 orthology",
        "species": "Nomascus leucogenys",
        "assembly": "ASM4011310v1",
        "haplotype_or_phase": "Haplotype 2",
        "assembly_accession": "GCA_040113105.1",
        "public_data_location": ncbi("GCA_040113105.1"),
        "local_fasta": str(ROOT / "phylogeny_origin_acro1_full/inputs/ncbi_downloads/expanded_individual/GCA_040113105.1_ASM4011310v1_genomic.fna"),
        "local_annotation": "Not required for the anchor-based HSA19 analysis",
        "analysis_role": "HSA19 orthology",
    },
    {
        "dataset_group": "HSA19 orthology",
        "species": "Macaca mulatta",
        "assembly": "Mmul_10",
        "haplotype_or_phase": "Haploid representation",
        "assembly_accession": "GCF_003339765.1",
        "public_data_location": ncbi("GCF_003339765.1"),
        "local_fasta": str(ROOT / "phylogeny_origin_acro1_full/inputs/ncbi_downloads/expanded_individual/GCF_003339765.1_Mmul_10_genomic.fna"),
        "local_annotation": "Not required for the anchor-based HSA19 analysis",
        "analysis_role": "HSA19 orthology",
    },
    {
        "dataset_group": "HSA19 orthology",
        "species": "Callithrix jacchus",
        "assembly": "calJac240_pri",
        "haplotype_or_phase": "Primary",
        "assembly_accession": "GCA_049354715.1",
        "public_data_location": ncbi("GCA_049354715.1"),
        "local_fasta": str(ROOT / "phylogeny_origin_acro1_full/inputs/ncbi_downloads/core_outgroups/ncbi_dataset/data/GCA_049354715.1/GCA_049354715.1_calJac240_pri_genomic.fna"),
        "local_annotation": "Not required for the anchor-based HSA19 analysis",
        "analysis_role": "HSA19 orthology and fragmented ACRO-locus analysis",
    },
    {
        "dataset_group": "HSA19 orthology",
        "species": "Leontopithecus rosalia",
        "assembly": "mLeoRos1.hap1",
        "haplotype_or_phase": "Haplotype 1",
        "assembly_accession": "GCA_054824895.1",
        "public_data_location": ncbi("GCA_054824895.1"),
        "local_fasta": str(ROOT / "phylogeny_origin_acro1_full/inputs/ncbi_downloads/expanded_individual/GCA_054824895.1_mLeoRos1.hap1_genomic.fna"),
        "local_annotation": "Not required for the anchor-based HSA19 analysis",
        "analysis_role": "HSA19 orthology",
    },
    {
        "dataset_group": "HSA19 orthology",
        "species": "Saimiri boliviensis",
        "assembly": "mSaiBol1.pri",
        "haplotype_or_phase": "Primary",
        "assembly_accession": "GCA_048565385.1",
        "public_data_location": ncbi("GCA_048565385.1"),
        "local_fasta": str(ROOT / "phylogeny_origin_acro1_full/inputs/ncbi_downloads/expanded_individual/GCA_048565385.1_mSaiBol1.pri_genomic.fna"),
        "local_annotation": "Not required for the anchor-based HSA19 analysis",
        "analysis_role": "HSA19 orthology",
    },
    {
        "dataset_group": "Human reference",
        "species": "Homo sapiens",
        "assembly": "T2T-CHM13v2.0",
        "haplotype_or_phase": "Effectively haploid complete reference",
        "assembly_accession": "GCA_009914755.4 (GenBank); GCF_009914755.1 (RefSeq)",
        "public_data_location": ncbi("GCA_009914755.4"),
        "local_fasta": "data/primate_genomes/chm13v2.0.fa.gz",
        "local_annotation": str(ROOT / "chm13/chm13v2.0.cenSatv2.1.bed"),
        "analysis_role": "Locus definition, PRO-seq projection, and reference-array analyses",
    },
]


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def markdown_escape(value: str) -> str:
    return value.replace("|", "\\|").replace("\n", "<br>")


def build_hprc_rows() -> list[dict[str, str]]:
    with HPRC_MANIFEST.open() as handle:
        arrays = list(csv.DictReader(handle, delimiter="\t"))
    used_pairs = sorted({(row["individual_id"], row["haplotype"]) for row in arrays})

    with HPRC_INDEX.open() as handle:
        release_rows = list(csv.DictReader(handle))
    by_pair = {(row["sample_id"], row["haplotype"]): row for row in release_rows}
    missing = [pair for pair in used_pairs if pair not in by_pair]
    if missing:
        raise RuntimeError(f"HPRC release index is missing {len(missing)} analyzed haplotypes: {missing[:10]}")

    rows = []
    for sample_id, haplotype in used_pairs:
        row = by_pair[(sample_id, haplotype)]
        rows.append(
            {
                "sample_id": sample_id,
                "haplotype": haplotype,
                "phasing": row["phasing"],
                "assembly_name": row["assembly_name"],
                "genbank_accession": row["genbank_accession"],
                "assembly_method": row["assembly_method"],
                "assembly_method_version": row["assembly_method_version"],
                "assembly_date": row["assembly_date"],
                "source": row["source"],
                "public_assembly_location": row["assembly"],
                "local_release_index": str(HPRC_INDEX),
                "local_array_manifest": str(HPRC_MANIFEST),
            }
        )
    if len(rows) != 462 or len({row["genbank_accession"] for row in rows}) != 462:
        raise RuntimeError("Expected 462 analyzed HPRC haplotypes with 462 unique GenBank accessions")
    return rows


def write_markdown(path: Path, hprc_rows: list[dict[str, str]]) -> None:
    columns = [
        ("Species", "species"),
        ("Assembly / release", "assembly"),
        ("Haplotype representation", "haplotype_or_phase"),
        ("Assembly accession(s)", "assembly_accession"),
        ("Public data location", "public_data_location"),
        ("Local FASTA used", "local_fasta"),
        ("Local annotation", "local_annotation"),
        ("Analysis role", "analysis_role"),
    ]
    lines = [
        "# Supplementary Table X. Genome assemblies and associated data sets",
        "",
        "Assembly accessions and data locations used in the ACRO analyses. For the five nonhuman great-ape diploid FASTAs, the listed GenBank accessions are the component primary/alternate or maternal/paternal assemblies corresponding to the combined GenomeArk FASTA. T2T-CHM13v2.0 was used as an additional human reference and is not counted among the 14 assemblies in the HSA19 comparative set.",
        "",
        "| " + " | ".join(title for title, _ in columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for row in comparative_rows:
        lines.append("| " + " | ".join(markdown_escape(row[key]) for _, key in columns) + " |")

    lines.extend(
        [
            "",
            "## HPRC Release 2 cohort",
            "",
            "| Cohort | Assemblies analyzed | Public accession and data location | Local analysis inputs | Use |",
            "| --- | ---: | --- | --- | --- |",
            "| HPRC Release 2 v1.0 | 462 phased haplotype assemblies from 231 individuals | 462 unique GenBank accessions and exact S3 FASTA locations are listed in `SUPPLEMENTARY_TABLE_HPRC_ASSEMBLIES.tsv`; umbrella BioProject PRJNA730822 | `NTRIAGE_OCT26/manuscript_output/hprc_telomeric_chromosome_assignment/inputs/assemblies_release2_v1.0.index.csv`; `NTRIAGE_OCT26/manifest.validated.tsv`; derived array FASTAs in `NTRIAGE_OCT26/telomeric_arrays/` | Population-level ACRO-array, NTRiage, transition, and promoter-conservation analyses |",
            "",
            f"All {len(hprc_rows)} analyzed sample–haplotype pairs matched the official Release 2 index, and all had a nonmissing, unique GenBank accession. Four index rows not analyzed here (CHM13, GRCh38, and the two HG002 reference haplotypes) were excluded from the cohort table.",
            "",
            "### Notes",
            "",
            "1. Local paths document the exact files used during analysis; public locations are the archival or consortium sources from which the corresponding assemblies can be obtained.",
            "2. The six great-ape entries are diploid analysis files and therefore list two component assembly accessions where applicable. Each entry counts as one comparative genome in the 14-assembly HSA19 analysis; the two siamang haplotypes were analyzed as separate assembly records.",
            "3. The HPRC table is kept as a separate machine-readable file because 462 assembly rows are not practical in the main Methods table.",
        ]
    )
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    comp_fields = list(comparative_rows[0])
    write_tsv(OUT / "SUPPLEMENTARY_TABLE_ASSEMBLIES.tsv", comparative_rows, comp_fields)
    hprc_rows = build_hprc_rows()
    hprc_fields = list(hprc_rows[0])
    write_tsv(OUT / "SUPPLEMENTARY_TABLE_HPRC_ASSEMBLIES.tsv", hprc_rows, hprc_fields)
    write_markdown(OUT / "SUPPLEMENTARY_TABLE_ASSEMBLIES.md", hprc_rows)


if __name__ == "__main__":
    main()
