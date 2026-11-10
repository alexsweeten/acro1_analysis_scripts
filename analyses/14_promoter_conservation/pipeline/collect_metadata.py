#!/usr/bin/env python3
"""Write analysis_metadata.tsv: software versions, modules, command lines, seeds, checksums, resources."""
import os
import platform
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pp_common as C
import run_array as RA


def sh(cmd):
    try:
        r = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=120)
        return (r.stdout + r.stderr).strip()
    except Exception as e:
        return f"unavailable: {e}"


def main():
    rows = []
    add = lambda k, v, cat="": rows.append((cat, k, str(v).replace("\n", " | ")))
    add("date", time.strftime("%Y-%m-%dT%H:%M:%S"), "run")
    add("host", platform.node(), "run")
    add("pipeline_version", RA.pipeline_version(), "run")
    for f in ("pp_common.py", "pp_core.py", "run_array.py", "pp_stats.py", "aggregate.py", "cluster_models.R", "make_figures.R",
              "00_stage_references.py", "01_build_manifest.py", "array_task.sbatch", "submit.sh"):
        add(f"sha256:{f}", C.sha256_file(C.ROOT / "scripts" / f)[:16], "scripts")
    add("python", sys.version.split()[0], "software")
    import numpy, pandas, scipy
    for m in (numpy, pandas, scipy):
        add(m.__name__, m.__version__, "software")
    add("nhmmer", sh("nhmmer -h | sed -n 2p"), "software")
    add("mafft", sh("mafft --version 2>&1 | head -1"), "software")
    add("samtools", sh("samtools --version | head -1"), "software")
    add("R", sh("R --version | head -1"), "software")
    add("R_packages", sh("""Rscript -e 'for(p in c("ggplot2","patchwork","scales","ragg","lme4","lmerTest","sandwich")) cat(p,as.character(packageVersion(p)),"; ")' 2>/dev/null"""), "software")
    add("modules_loaded", os.environ.get("LOADEDMODULES", sh("module list 2>&1 | tail -3")), "software")
    add("module_load_command", "module load hmmer/3.4 mafft/7.526 samtools/1.21 (per-array tasks); module load R/4.4.2 (stats/figures)", "software")
    add("environment_file", "none (HPC lmod modules; conda base python 3.9.15 on PATH)", "software")
    for _, r in pd.read_csv(C.REF / "reference_checksums.tsv", sep="\t").iterrows():
        add(f"sha256:{r.file}", f"{r.observed_sha256} ({r.status})", "reference")
    add("promoter_300col", "columns 78-127 inclusive (50)", "fixed_parameters")
    add("promoter_fullunit_9090col", "columns 74-123 inclusive (50)", "fixed_parameters")
    add("PROseq_role",
        "prior CHM13 PRO-seq nominated/fixed promoter columns 78-127; not used for HPRC anchors or boundaries, control-window selection, or target-array re-selection/repositioning",
        "fixed_parameters")
    for k, v in (("seed", C.SEED), ("signflip_draws_primary", C.N_SIGNFLIP_PRIMARY), ("signflip_draws_backbone", C.N_SIGNFLIP_BACKBONE),
                 ("group_permutation_draws", C.N_GROUP_PERM), ("bootstrap_draws_mean_ci", C.N_BOOT), ("cluster_bootstrap_draws", C.N_CLUSTER_BOOT),
                 ("rng", "numpy PCG64 via SeedSequence([seed, crc32(test label)]); R: set.seed(seed + sum(utf8ToInt(analysis)))")):
        add(k, v, "randomisation")
    cmds = next(iter(sorted(C.PER_ARRAY.glob("*/logs/commands.txt"))), None)
    if cmds:
        add("per_array_command_lines (example)", cmds.read_text().strip(), "commands")
    add("command:nhmmer", f"nhmmer --dna --tblout HITS.tbl --noali -E {C.NHMMER_E} --cpu $SLURM_CPUS_PER_TASK chr21_l1pa10.hmm ARRAY.fa", "commands")
    add("command:anchors_msa", "mafft --quiet --thread $SLURM_CPUS_PER_TASK --addfragments l1pa10_anchors.fa --keeplength l1pa10_reference_alignment.fa", "commands")
    add("command:structural_msa", "mafft --genafpair --maxiterate 1000 --ep 0 --thread $SLURM_CPUS_PER_TASK units_transcription.fa", "commands")
    add("command:backbone_map", "mafft --quiet --thread 1 --retree 1 --maxiterate 0 --add units_transcription.fa --keeplength acro_units_transcription.einsi.fa", "commands")
    jobs = C.ROOT / "logs" / "submitted_jobs.txt"
    if jobs.exists():
        add("slurm_jobs_submitted", jobs.read_text().split(), "slurm")
    pd.DataFrame(rows, columns=["category", "key", "value"]).to_csv(C.ROOT / "analysis_metadata.tsv", sep="\t", index=False)
    print(f"wrote {len(rows)} rows")


if __name__ == "__main__":
    main()
