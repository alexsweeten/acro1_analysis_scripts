#!/bin/bash
# One manifest array per Slurm task. Reads the task's manifest index from $TASKLIST (line = SLURM_ARRAY_TASK_ID).
# Submit via submit.sh, which sets resources from the pilot and TASKLIST. Never use mafft --thread -1.
#SBATCH --job-name=pp2794
#SBATCH --partition=norm
#SBATCH --cpus-per-task=2
#SBATCH --output=/data/sweetenap/hprc/oct_21/promoter_purifying_2794/logs/slurm/%x_%A_%a.out
set -euo pipefail
ROOT=/data/sweetenap/hprc/oct_21/promoter_purifying_2794
module load hmmer/3.4 mafft/7.526 samtools/1.21
cd "$ROOT"
IDX=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${TASKLIST:?TASKLIST not set}")
echo "task ${SLURM_ARRAY_TASK_ID} -> manifest index ${IDX} on $(hostname), cpus=${SLURM_CPUS_PER_TASK}"
python3 scripts/run_array.py --task-index "$IDX" --cpus "${SLURM_CPUS_PER_TASK}" ${FORCE:+--force}
