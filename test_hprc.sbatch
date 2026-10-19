#!/bin/bash
# Smoke test for euk-bot on TAMU HPRC: checks the modules, compute-node internet access (needed to build
# conda envs), and that Nextflow can submit SLURM jobs. Runs the tiny test profile as a stub run, so no
# real tools or conda envs are used; it should finish in a few minutes.
#
# Usage (from the repo directory):
#   sbatch test_hprc.sbatch
#
#SBATCH --job-name=euk-bot-test
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=euk-bot-test_%j.log
#SBATCH --account=132704181195
#SBATCH --partition=short

set -uo pipefail

module purge
module load Nextflow/25.10.2
module load Miniforge3/26.3.2-2

USER_SCRATCH=${USER_SCRATCH:-/scratch/user/${USER}}
export NXF_HOME=${USER_SCRATCH}/.cache/nextflow
export XDG_CACHE_HOME=${USER_SCRATCH}/.cache
export NXF_OPTS='-Xms500m -Xmx2g'
mkdir -p "${NXF_HOME}"

echo "== Node: $(hostname)   Dir: $(pwd)"
echo "== Versions"
nextflow -version | grep -i version
mamba --version | head -1

echo "== Internet from compute node (conda envs are built here by the head job)"
for url in https://conda.anaconda.org https://repo.anaconda.com; do
    if curl -sS -m 15 -o /dev/null "${url}"; then echo "OK    ${url}"; else echo "FAIL  ${url}"; fi
done

echo "== Stub run of the test profile, tasks submitted to SLURM"
nextflow run main.nf \
    -profile test,slurm \
    -stub-run \
    --slurm_account 132704181195 \
    -work-dir "${USER_SCRATCH}/euk-bot/work_test"
status=$?

echo "== Nextflow exit status: ${status}"
exit ${status}
