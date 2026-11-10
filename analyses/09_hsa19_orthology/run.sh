#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."
python3 phylogeny_origin_acro1_full/scripts/run_full_analysis.py --threads 8 --resume
Rscript phylogeny_origin_acro1_full/scripts/plot_full_analysis.R phylogeny_origin_acro1_full

