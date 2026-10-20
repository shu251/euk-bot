#!/bin/bash
# Pre-build euk-bot's conda envs, or pull its Singularity images, into the shared caches so pipeline runs
# start without fetching anything. Stub runs on the local executor: no real tools run, but Nextflow creates
# every env / pulls every image the run touches. Run twice to cover the optional steps.
# HPRC policy: Singularity pulls must run on a compute node, so submit it with sbatch; compute nodes reach the
# internet through the WebProxy module.
#
# Usage (from the repo directory):
#   sbatch prebuild_envs.sh singularity   # Singularity images (compute node, required)
#   bash prebuild_envs.sh                 # conda envs (login node in tmux, or sbatch)
# Rerun after changing any module's conda/container directive; existing envs and images are reused.
#
#SBATCH --job-name=euk-bot-prebuild
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=30G
#SBATCH --output=euk-bot-prebuild_%j.log
#SBATCH --account=132704181195
#SBATCH --partition=medium

set -euo pipefail

PROFILE=${1:-conda}   # conda | singularity

module purge
module load Nextflow/25.10.2
module load Anaconda3/2026.07-1
module load WebProxy   # internet access from compute nodes

# Same cache locations as run_hprc.sbatch, so the head job finds these envs
USER_SCRATCH=${USER_SCRATCH:-/scratch/user/${USER}}
GROUP_SPACE=${GROUP_SPACE:-/scratch/group/hu-lab}
export NXF_HOME=${USER_SCRATCH}/.cache/nextflow
export XDG_CACHE_HOME=${USER_SCRATCH}/.cache
export NXF_CONDA_CACHEDIR=${GROUP_SPACE}/.cache/euk-bot/conda_envs
export CONDA_PKGS_DIRS=${USER_SCRATCH}/.cache/conda_pkgs   # many small files; keep off group quota
export NXF_SINGULARITY_CACHEDIR=${GROUP_SPACE}/.cache/euk-bot/singularity_images
export SINGULARITY_CACHEDIR=${TMPDIR:-${USER_SCRATCH}/.cache}/.singularity   # layer cache; HPRC: local /tmp, never $HOME
export SINGULARITY_TMPDIR=${TMPDIR:-${USER_SCRATCH}/.cache}
export NXF_OPTS='-Xms500m -Xmx2g'
umask 002
mkdir -p "${NXF_HOME}" "${NXF_CONDA_CACHEDIR}" "${CONDA_PKGS_DIRS}" "${NXF_SINGULARITY_CACHEDIR}" "${SINGULARITY_CACHEDIR}"

WORK_DIR=${USER_SCRATCH}/euk-bot/work_prebuild

echo "== Run 1: default branches"
nextflow run main.nf -profile test,${PROFILE} -stub-run -work-dir "${WORK_DIR}" \
    --outdir "${USER_SCRATCH}/euk-bot/results_prebuild"

echo "== Run 2: optional steps: eggNOG-mapper, assembly QC, read profiling (data dirs are placeholders; stub only)"
nextflow run main.nf -profile test,${PROFILE} -stub-run -work-dir "${WORK_DIR}" \
    --outdir "${USER_SCRATCH}/euk-bot/results_prebuild" \
    --eggnog_data_dir tests/data --run_assembly_qc --busco_db_dir tests/data \
    --run_read_profiling --kaiju_db tests/data

if [[ ${PROFILE} == singularity ]]; then
    echo "== Images in ${NXF_SINGULARITY_CACHEDIR}:"; ls -1 "${NXF_SINGULARITY_CACHEDIR}"
else
    echo "== Envs in ${NXF_CONDA_CACHEDIR}:"; ls -1 "${NXF_CONDA_CACHEDIR}"
fi
