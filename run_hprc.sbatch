#!/bin/bash
# Head job for euk-bot on TAMU HPRC. Nextflow itself submits each pipeline step as its own SLURM job;
# this job just keeps Nextflow alive, so it needs little CPU/memory but a long walltime.
# Compute nodes have no internet: build the conda envs first on a login node with prebuild_envs.sh.
#
# Usage (from the repo directory):
#   sbatch run_hprc.sbatch samplesheet_test3.csv /scratch/group/hu-lab/<project>/euk-bot_test3 [extra nextflow args...]
#
#SBATCH --job-name=euk-bot
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=euk-bot_%j.log
#SBATCH --account=132704181195
#SBATCH --partition=long

set -euo pipefail

SAMPLESHEET=${1:?"usage: sbatch run_hprc.sbatch <samplesheet.csv> <outdir> [nextflow args]"}
OUTDIR=${2:?"usage: sbatch run_hprc.sbatch <samplesheet.csv> <outdir> [nextflow args]"}
shift 2

module purge
module load Nextflow/25.10.2
module load Anaconda3/2026.07-1   # provides conda for -profile conda

# Storage layout (nothing under $HOME, which has a small quota on Grace):
#   personal scratch : Nextflow home, temp caches, work/ (large, disposable)
#   group space      : conda envs shared across the lab; put OUTDIR here too
USER_SCRATCH=${USER_SCRATCH:-/scratch/user/${USER}}
GROUP_SPACE=${GROUP_SPACE:-/scratch/group/hu-lab}
WORK_DIR=${WORK_DIR:-${USER_SCRATCH}/euk-bot/work}

export NXF_HOME=${USER_SCRATCH}/.cache/nextflow
export XDG_CACHE_HOME=${USER_SCRATCH}/.cache
export NXF_CONDA_CACHEDIR=${GROUP_SPACE}/.cache/euk-bot/conda_envs
export CONDA_PKGS_DIRS=${USER_SCRATCH}/.cache/conda_pkgs   # many small files; keep off group quota
umask 002   # keep shared envs group-writable
mkdir -p "${NXF_HOME}" "${NXF_CONDA_CACHEDIR}" "${CONDA_PKGS_DIRS}" "${WORK_DIR}"

export NXF_OPTS='-Xms1g -Xmx4g'

nextflow run main.nf \
    -profile slurm,conda \
    --input "${SAMPLESHEET}" \
    --outdir "${OUTDIR}" \
    --slurm_account "${SLURM_ACCOUNT:-132704181195}" \
    -work-dir "${WORK_DIR}" \
    -resume \
    "$@"
