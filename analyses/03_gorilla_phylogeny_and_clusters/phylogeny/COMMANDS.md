# Commands used for the all-gorilla ACRE phylogeny

All commands were run from the ACREAGE project root:

```text
${ACREAGE_ROOT}
```

The commands below are the reproducible successful workflow. Input preparation
removes prior six-primate cluster labels; the de novo sweep refits scaling and
k-means using gorilla ACREs only.

## 1. Extract all gorilla ACRE sequences, features, and metadata

```bash
python3 \
  acro_gorilla_only/prepare_gorilla.py
```

## 2. Run the gorilla-only cluster-number sweep

```bash
OPENBLAS_NUM_THREADS=8 \
MPLCONFIGDIR=acro_gorilla_only/matplotlib_cache \
XDG_CACHE_HOME=acro_gorilla_only/font_cache \
python3 \
  acro_gorilla_only/sweep_clusters.py
```

This evaluates k=2–16 with ten k-means++ starts per k (seeds 1–10). The
selection criterion is the equal-weight consensus stored in
`cluster_sweep_parameters.json`; it selected k=6.

## 3. Align all 1,232 gorilla ACREs and normalize orientation

```bash
/usr/bin/time -l \
mafft \
  --thread 8 \
  --retree 1 \
  --maxiterate 0 \
  --adjustdirectionaccurately \
  acro_gorilla_only/gorilla_acres.fasta \
  > acro_gorilla_only/gorilla_acres.aligned.raw.fasta \
  2> acro_gorilla_only/mafft.log
```

```bash
python3 \
  acro_gorilla_only/finalize_alignment.py
```

`finalize_alignment.py` removes MAFFT `_R_` name prefixes, verifies the exact
tip set, records sequence orientation, and writes alignments retaining columns
with valid A/C/G/T bases in at least 2%, 5%, and 10% of tips.

## 4. Infer the maximum-likelihood tree

```bash
/usr/bin/time -l \
iqtree3 \
  -s acro_gorilla_only/gorilla_acres.aligned.occupancy10.fasta \
  -pre acro_gorilla_only/gorilla_acres \
  -m GTR+F+R4 \
  --fast \
  --alrt 1000 \
  -T 8 \
  --seed 1 \
  > acro_gorilla_only/iqtree.stdout.log \
  2> acro_gorilla_only/iqtree.stderr.log
```

The scalable `--fast` search and 1,000 SH-aLRT replicates retain all tips.
IQ-TREE ultrafast bootstrap is not combined with `--fast`; node labels in the
tree are therefore SH-aLRT support values.

## 5. Calculate tree statistics and draw four views

```bash
Rscript acro_gorilla_only/analyze_tree.R \
  > acro_gorilla_only/analyze_tree.log 2>&1
```

```bash
Rscript acro_gorilla_only/draw_trees.R \
  > acro_gorilla_only/draw_trees.log 2>&1
```

Both layouts show the same unrooted maximum-likelihood topology. For display,
the plotting script chooses the internal node minimizing the maximum topological
distance to all tips as a neutral center root; this changes no unrooted split or
support. The figures are cladograms for legibility, while inferred branch
lengths remain in `gorilla_acres.treefile`.

## 6. Validation and software-version capture

The exact validation and version-report commands, together with their results,
are recorded in `validation.log` and `software_versions.tsv` when the workflow
finishes.

```bash
python3 \
  -m py_compile \
  acro_gorilla_only/prepare_gorilla.py \
  acro_gorilla_only/sweep_clusters.py \
  acro_gorilla_only/finalize_alignment.py \
  acro_gorilla_only/validate_results.py
```

```bash
Rscript -e \
  'parse(file="acro_gorilla_only/analyze_tree.R"); parse(file="acro_gorilla_only/draw_trees.R")'
```

```bash
python3 \
  acro_gorilla_only/validate_results.py \
  > acro_gorilla_only/validation.log 2>&1
```

## Managed-environment note

IQ-TREE completed normally and recorded the final report, tree, and completion
timestamp. In this managed macOS sandbox, the surrounding `/usr/bin/time -l`
process subsequently returned status 1 because it could not query
`sysctl kern.clockrate`; this is visible in `iqtree.stderr.log` and did not
affect IQ-TREE or its outputs. The MAFFT timing wrapper ran outside that sandbox
restriction and completed normally.
