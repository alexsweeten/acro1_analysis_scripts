# ACRO--HSat1A association

This is the final six-species HG002 analysis used for the abundance correlation
and direct-boundary permutation tests. It supersedes the earlier three-species
exploratory script.

```bash
Rscript run_six_species_analysis.R \
  --input-dir /path/to/censat_inputs \
  --output-dir results \
  --permutations 100000 \
  --seed 20261021
```

The input directory must contain the assembly-matched CenSat BED files and
FASTA index files named in the script. The analysis collapses touching records
into runs, preserves their coordinates and gap structure during within-
chromosome label permutations, and reports exact species-level permutation
tests separately from the spatial null model.

