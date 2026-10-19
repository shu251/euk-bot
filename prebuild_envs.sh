#!/bin/bash
# Pre-build euk-bot's conda envs on a TAMU HPRC login node. Grace compute nodes have no internet, so the
# head job can't create envs itself; it reuses the ones built here from the shared cache.
# Stub runs on the local executor: no real tools run, but Nextflow creates every env the run touches.
# Run twice to cover both assemblers and eggNOG-mapper. Takes a while; run it inside tmux.
#
# Usage (from the repo directory, on a login node):
#   bash prebuild_envs.sh
# Rerun after changing any module's conda directive; existing envs are reused.

set -euo pipefail

module purge
module load Nextflow/25.10.2
module load Anaconda3/2026.07-1

# Same cache locations as run_hprc.sbatch, so the head job finds these envs
USER_SCRATCH=${USER_SCRATCH:-/scratch/user/${USER}}
GROUP_SPACE=${GROUP_SPACE:-/scratch/group/hu-lab}
export NXF_HOME=${USER_SCRATCH}/.cache/nextflow
export XDG_CACHE_HOME=${USER_SCRATCH}/.cache
export NXF_CONDA_CACHEDIR=${GROUP_SPACE}/.cache/euk-bot/conda_envs
export CONDA_PKGS_DIRS=${GROUP_SPACE}/.cache/conda_pkgs
export NXF_OPTS='-Xms500m -Xmx2g'
umask 002
mkdir -p "${NXF_HOME}" "${NXF_CONDA_CACHEDIR}" "${CONDA_PKGS_DIRS}"

WORK_DIR=${USER_SCRATCH}/euk-bot/work_prebuild

echo "== Run 1: default branches (rnaSPAdes)"
nextflow run main.nf -profile test,conda -stub-run -work-dir "${WORK_DIR}" \
    --outdir "${USER_SCRATCH}/euk-bot/results_prebuild"

echo "== Run 2: MEGAHIT + eggNOG-mapper (data dir is a placeholder; stub only)"
nextflow run main.nf -profile test,conda -stub-run -work-dir "${WORK_DIR}" \
    --outdir "${USER_SCRATCH}/euk-bot/results_prebuild" \
    --assembler megahit --eggnog_data_dir tests/data

echo "== Envs in ${NXF_CONDA_CACHEDIR}:"
ls -1 "${NXF_CONDA_CACHEDIR}"
