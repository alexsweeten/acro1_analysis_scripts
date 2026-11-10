# ACRO1 manuscript analysis scripts

This repository contains the custom analysis code used for the ACRO1/ACRE
chapter, **Evolutionary analysis of the great ape ACRO1 satellite**. It is a
manuscript-specific companion to [NTRiage](https://github.com/alexsweeten/NTRiage),
which is maintained and versioned separately.

The repository intentionally contains code, compact provenance files, curated
transition audits, and summary results. It does not redistribute genome
assemblies, HPRC array FASTAs, alignment indexes, bigWig files, or other large
primary data. Assembly accessions are listed in
`analyses/15_assembly_tables/SUPPLEMENTARY_TABLE_ASSEMBLIES.tsv` and
`analyses/15_assembly_tables/SUPPLEMENTARY_TABLE_HPRC_ASSEMBLIES.tsv`.

## Analysis crosswalk

| Methods analysis | Directory | Principal entry point |
| --- | --- | --- |
| ACRE extraction, ACRO1 allele classification, consensuses, and marmoset component analysis | `analyses/01_acre_extraction_and_alleles` | `scripts/prepare_joint_input.py`; `scripts/label_six_primate_subfamilies.py` |
| Species-balanced ACRE phylogeny | `analyses/02_species_balanced_phylogeny` | `run_all.sh` |
| Gorilla phylogeny, G1--G6 clustering, and spatial association | `analyses/03_gorilla_phylogeny_and_clusters` | `phylogeny/prepare_gorilla.py`; `clustering/analyze_gorilla_clustering.py` |
| Six-frame ACRO1 ORF analysis | `analyses/04_orf_analysis` | `scripts/run_de_novo_orf_test.py` |
| ACRO1--POLR2A amino-acid composition test | `analyses/05_polr2a_composition` | `run_analysis.py` |
| Family-38 PRDM9-A motif scan | `analyses/06_prdm9_motif` | `analyze_family38.py` |
| ACRO--HSat1A abundance and boundary enrichment | `analyses/07_acro_hsat1a` | `run_six_species_analysis.R` |
| Telomere/rDNA positions and Pan reciprocal short-arm synteny | `analyses/08_telomere_rdna_pan_synteny` | `scripts/analyze_positions.R`; `scripts/build_all.sh` |
| HSA19 orthology and reciprocal-anchor analysis | `analyses/09_hsa19_orthology` | `run.sh` |
| HG002, gorilla, and HPRC NTRiage runs and crosswalks | `analyses/10_ntriage_runs` | `run_all.sh` |
| HPRC transition screening, localization, and manual audit | `analyses/11_hprc_transitions` | `scripts/find_recombination.py`; `scripts/localize_positive_breakpoints.py` |
| ACRE-scale HOR scoring | `analyses/12_hor_score` | `scripts/compute_hor_score.py` |
| Great-ape and HPRC intra-array divergence | `analyses/13_divergence` | analysis-specific Python entry points |
| Candidate L1PA10 interval conservation | `analyses/14_promoter_conservation` | `pipeline/run_array.py`; array- and ACRE-level scripts |
| Assembly and accession tables | `analyses/15_assembly_tables` | `build_assembly_tables.py` |

`METHODS_CROSSWALK.tsv` provides the same mapping in machine-readable form.

## Installation

Create the mixed Python/R/command-line environment with Conda or Mamba:

```bash
mamba env create -f environment.yml
mamba activate acro1-analysis
```

NTRiage v0.3.0 must be installed separately. The analyses also use MAFFT,
IQ-TREE, minimap2, HMMER, Mash, RepeatMasker, Tandem Repeat Finder, and
ModDotPlot as indicated in the Methods and individual scripts.

## Inputs

Large inputs are supplied by their original repositories rather than copied
here. Start from `config/input_paths.example.tsv`, replace each placeholder
with a local path, and consult the relevant analysis directory for the expected
columns and file formats. Coordinates use BED-style zero-based, half-open
intervals unless an analysis README states otherwise.

Some archived scripts preserve the exact multi-stage workflow used during the
analysis and therefore require their input variables or path arguments to be
set for a new installation. Files containing manual HPRC transition decisions
are retained under `analyses/11_hprc_transitions/audit` because those decisions
cannot be regenerated from software alone.

## Reproducibility status

See [REPRODUCIBILITY_STATUS.md](REPRODUCIBILITY_STATUS.md) for the distinction
between directly rerunnable workflows, archived provenance code, manual steps,
and unresolved source-data provenance. In particular, the PRO-seq track source
must be documented before a fully self-contained promoter-analysis release can
be claimed.

## Citation and license

Citation metadata are provided in `CITATION.cff`. Code in this repository is
released under the MIT License. Primary data retain the licenses and terms of
their original providers.

