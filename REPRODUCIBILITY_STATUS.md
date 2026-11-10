# Reproducibility status

This repository is a curated archive of the custom code underlying the Methods,
not a redistribution of the multi-gigabyte primary data.

## Directly parameterized workflows

- ACRO--HSat1A abundance and adjacency analysis.
- Family-38/PRDM9-A motif analysis.
- ACRO1--POLR2A JSD simulations.
- Generic ACRE-scale HOR scoring.

## Workflows requiring public input data and local configuration

- Great-ape and gorilla phylogenies.
- Pan reciprocal short-arm synteny.
- HSA19 reciprocal-anchor orthology.
- NTRiage cohort analyses.
- Great-ape and HPRC divergence analyses.
- Candidate-interval conservation analyses.

The original scripts for these analyses are included, but genome FASTAs,
repeat annotations, HPRC array FASTAs, large alignments, and software indexes
are not. Input accessions and compact provenance manifests are included where
available.

## Manual analyses

Transition localization combined an automated ModDotPlot screen and coordinate
localization with manual review. The row-level inclusion, removal, and boundary
shift records are included in `analyses/11_hprc_transitions/audit`. These audit
tables are part of the reproducible record and should be cited with the code.

## Reconstructed workflow

The original standalone HPRC HOR-score generator was not retained. The public
`compute_hor_score.py` is a generalized reconstruction of the preserved gorilla
implementation and the archived HPRC output schema. It implements the same
per-array MAFFT alignment, offset-block detection, distinct-partner counting,
array-size normalization, and 90%, 92%, and 96% identity thresholds. The
archived HPRC score table is retained for result verification.

## Unresolved provenance

The local promoter workflow used strand-specific CHM13 PRO-seq bigWig tracks
named `PROseq_k100_neg.bw` and `PROseq_k100_pos.bw`. Their accession, biological
source/cell type, and original publication were not recorded in the surviving
analysis directory. Those files are not redistributed here. This provenance
must be resolved before the manuscript describes the promoter workflow as
fully reproducible.

## External software

NTRiage source code is maintained at <https://github.com/alexsweeten/NTRiage>.
AniAnn's, RepeatMasker, MAFFT, IQ-TREE, minimap2, HMMER, Mash, TRF,
ModDotPlot, and SVbyEye are external dependencies and are not vendored here.

