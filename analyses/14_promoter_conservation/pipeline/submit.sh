#!/bin/bash
# Usage: submit.sh TASKS_FILE MEM TIME CONCURRENCY [--force] [--tag NAME]
#   TASKS_FILE: one manifest task_index per line.  Splits into <=1000-task Slurm arrays (cluster convention).
# Prints the job ids. Rerun only failures with:  python3 scripts/find_incomplete.py > tasks.txt && submit.sh tasks.txt ...
set -euo pipefail
ROOT=/data/sweetenap/hprc/oct_21/promoter_purifying_2794
TASKS=$(readlink -f "$1"); MEM=$2; TIME=$3; CONC=$4; shift 4
FORCE=""; TAG="run"
while [ $# -gt 0 ]; do case $1 in --force) FORCE=1;; --tag) TAG=$2; shift;; esac; shift; done
mkdir -p "$ROOT/logs/slurm" "$ROOT/logs/tasklists"
STAMP=$(date +%Y%m%d_%H%M%S)
split -l 1000 -d -a 3 "$TASKS" "$ROOT/logs/tasklists/${TAG}_${STAMP}_"
for f in "$ROOT"/logs/tasklists/${TAG}_${STAMP}_*; do
  n=$(wc -l < "$f")
  sbatch --parsable --array=1-"$n"%"$CONC" --mem="$MEM" --time="$TIME" --job-name="pp_${TAG}" \
         --export=ALL,TASKLIST="$f"${FORCE:+,FORCE=1} "$ROOT/scripts/array_task.sbatch" | tee -a "$ROOT/logs/submitted_jobs.txt"
done
