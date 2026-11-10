#!/bin/bash
# Finalisation: metadata -> aggregation (array-weighted + cluster-aware inference) -> figures -> README.
# Run only after every manifest array has COMPLETE (aggregate.py refuses otherwise).
#SBATCH --job-name=pp_finalize
#SBATCH --partition=norm
#SBATCH --cpus-per-task=4
#SBATCH --mem=24g
#SBATCH --time=04:00:00
#SBATCH --output=/data/sweetenap/hprc/oct_21/promoter_purifying_2794/logs/slurm/finalize_%j.out
set -euo pipefail
ROOT=/data/sweetenap/hprc/oct_21/promoter_purifying_2794
module load hmmer/3.4 mafft/7.526 samtools/1.21 R/4.4.2
cd "$ROOT"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK} OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
python3 scripts/collect_metadata.py
python3 scripts/recheck_mapping.py
python3 scripts/aggregate.py ${ALLOW_PARTIAL:+--allow-partial}
Rscript scripts/make_figures.R combined
python3 scripts/make_readme.py
echo "finalize done"
