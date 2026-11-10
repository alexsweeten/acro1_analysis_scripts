# Candidate L1PA10 interval conservation

`pipeline/` contains the HPRC per-array alignment, composition-matched control,
divergence, aggregation, and statistical-test code. `array_level/` contains
the strict-majority H1--H3 dominant-array analysis used for the main figure.
`acre_level/` contains the same-cluster, within-parent-array sensitivity
analysis.

Large input FASTAs, per-array caches, and bigWig tracks are intentionally
omitted. Paths can be supplied with the environment variables used in the
scripts, including `ACRO1_WORKSPACE`, `ACRO1_PROMOTER_PIPELINE_DATA`,
`ACRO1_HPRC_UNITS`, and `ACRO1_HPRC_ASSIGNMENTS`.

The biological accession, cell type, and publication source of the surviving
local `PROseq_k100_neg.bw` and `PROseq_k100_pos.bw` files have not yet been
recovered. The code is archived, but this provenance must be supplied before
claiming a fully self-contained reproduction of interval nomination.

