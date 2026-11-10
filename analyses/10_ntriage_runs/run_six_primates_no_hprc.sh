#!/usr/bin/env bash
set -euo pipefail

ANALYSIS_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FOCUSED_MANIFEST="${ANALYSIS_ROOT}/inputs/six_primates_no_hprc_manifest.tsv"
FOCUSED_RESULTS="${ANALYSIS_ROOT}/results/six_primates_no_hprc"
FOCUSED_VALIDATION="${ANALYSIS_ROOT}/validation/joint_primates_no_hprc"

export MPLCONFIGDIR="${ANALYSIS_ROOT}/.mplconfig"
export XDG_CACHE_HOME="${ANALYSIS_ROOT}/.cache"
export LOKY_MAX_CPU_COUNT="${LOKY_MAX_CPU_COUNT:-8}"

python "${ANALYSIS_ROOT}/scripts/prepare_inputs.py" \
  --analysis-root "${ANALYSIS_ROOT}" \
  --six-primates-only
python "${ANALYSIS_ROOT}/scripts/run_ntriage.py" \
  --analysis-root "${ANALYSIS_ROOT}" \
  --manifest "${FOCUSED_MANIFEST}" \
  --summary-output "${ANALYSIS_ROOT}/provenance/six_primates_no_hprc/run_ntriage_summary.json" \
  --continue-on-error
python "${ANALYSIS_ROOT}/scripts/select_k.py" \
  --manifest "${FOCUSED_MANIFEST}"
python "${ANALYSIS_ROOT}/scripts/analyze_clusters.py" \
  --root "${ANALYSIS_ROOT}" \
  --manifest "${FOCUSED_MANIFEST}" \
  --joint-id joint_primates \
  --combined-output-dir "${FOCUSED_RESULTS}"
python "${ANALYSIS_ROOT}/scripts/summarize_results.py" \
  --root "${ANALYSIS_ROOT}" \
  --manifest "${FOCUSED_MANIFEST}" \
  --combined-output-dir "${FOCUSED_RESULTS}"
python "${ANALYSIS_ROOT}/scripts/plot_joint_primates_summary.py" \
  --root "${ANALYSIS_ROOT}" \
  --output-dir "${FOCUSED_RESULTS}/figures"
python "${ANALYSIS_ROOT}/scripts/validate_outputs.py" \
  --root "${ANALYSIS_ROOT}" \
  --manifest "${FOCUSED_MANIFEST}" \
  --joint-id joint_primates \
  --combined-output-dir "${FOCUSED_RESULTS}" \
  --validation-output-dir "${FOCUSED_VALIDATION}" \
  --allow-partial
