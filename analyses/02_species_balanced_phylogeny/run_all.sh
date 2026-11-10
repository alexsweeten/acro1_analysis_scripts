#!/usr/bin/env bash
set -euo pipefail

project="$(cd "$(dirname "$0")" && pwd)"
cd "$project"

"${PYTHON:-python3}" scripts/prepare_panel.py
mafft --auto --thread 4 inputs/panel_sequences.fasta \
  > results/panel_sequences.aligned.fasta \
  2> logs/mafft.log

"${IQTREE:-iqtree3}" \
  -s results/panel_sequences.aligned.fasta \
  -m GTR+F+G4 -bb 1000 -nt AUTO -redo \
  -pre results/panel_sequences

Rscript scripts/plot_trees.R
