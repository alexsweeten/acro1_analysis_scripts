#!/usr/bin/env bash
set -euo pipefail

ANALYSIS_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export MPLCONFIGDIR="${ANALYSIS_ROOT}/.mplconfig"
export XDG_CACHE_HOME="${ANALYSIS_ROOT}/.cache"
export LOKY_MAX_CPU_COUNT="${LOKY_MAX_CPU_COUNT:-8}"

python "${ANALYSIS_ROOT}/scripts/prepare_inputs.py" \
  --analysis-root "${ANALYSIS_ROOT}"
python "${ANALYSIS_ROOT}/scripts/run_ntriage.py" \
  --analysis-root "${ANALYSIS_ROOT}" \
  --continue-on-error
python "${ANALYSIS_ROOT}/scripts/select_k.py" \
  --manifest "${ANALYSIS_ROOT}/inputs/dataset_manifest.tsv"
python "${ANALYSIS_ROOT}/scripts/analyze_clusters.py" \
  --root "${ANALYSIS_ROOT}" \
  --manifest "${ANALYSIS_ROOT}/inputs/dataset_manifest.tsv"
python "${ANALYSIS_ROOT}/scripts/summarize_results.py" \
  --root "${ANALYSIS_ROOT}" \
  --manifest "${ANALYSIS_ROOT}/inputs/dataset_manifest.tsv"
python "${ANALYSIS_ROOT}/scripts/validate_outputs.py" \
  --root "${ANALYSIS_ROOT}" \
  --manifest "${ANALYSIS_ROOT}/inputs/dataset_manifest.tsv"
