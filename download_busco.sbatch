#!/bin/bash
# Download BUSCO lineage datasets for the optional assembly QC (--busco_db_dir). Uses the pipeline's BUSCO
# Singularity image, so pull images first (sbatch prebuild_envs.sh singularity). WebProxy gives the
# compute node internet access.
#
# Usage (from the repo directory):
#   sbatch download_busco.sbatch                                   # eukaryota_odb12
#   sbatch download_busco.sbatch eukaryota_odb12 alveolata_odb12   # any lineages; must match --busco_lineages
# Then run the pipeline with:
#   --run_assembly_qc --busco_db_dir /scratch/group/hu-lab/databases/busco
#
#SBATCH --job-name=busco-dl
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=busco-dl_%j.log
#SBATCH --account=132704181195
#SBATCH --partition=medium

set -euo pipefail

LINEAGES=("${@:-eukaryota_odb12}")
DATA_DIR=${DATA_DIR:-/scratch/group/hu-lab/databases/busco}
IMAGES=/scratch/group/hu-lab/.cache/euk-bot/singularity_images

module load WebProxy   # internet access from compute nodes

BUSCO_SIF=$(ls "${IMAGES}"/*busco-6.1.0*.img 2>/dev/null | head -1 || true)
[[ -n ${BUSCO_SIF} ]] || { echo "No BUSCO 6.1.0 image in ${IMAGES}; run prebuild_envs.sh singularity first" >&2; exit 1; }

umask 002   # keep shared databases group-readable
mkdir -p "${DATA_DIR}"
export SINGULARITY_CACHEDIR=${TMPDIR:-/scratch/user/${USER}/.cache}/.singularity

singularity exec -B "${DATA_DIR}" "${BUSCO_SIF}" busco --download_path "${DATA_DIR}" --download "${LINEAGES[@]}"

echo "== lineages in ${DATA_DIR}/lineages:"
ls "${DATA_DIR}/lineages"
