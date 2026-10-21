#!/bin/bash
# MMseqs2 protein DB of a built euk-bot reference version, for MetaEuk in the assembly branch (--metaeuk_db).
# Same proteins as the DIAMOND DB; built once with the pipeline's MetaEuk (Singularity image or conda env, from
# prebuild_envs.sh) into a read-only <version>_mmseqs folder beside the version folder.
#
# Usage (from the repo directory):
#   sbatch build_metaeuk_db.sbatch                  # v1
#   REF_VERSION=v2 sbatch build_metaeuk_db.sbatch
# Then run the pipeline with:
#   --metaeuk_db /scratch/group/hu-lab/databases/euk-bot-ref/v1_mmseqs
#
#SBATCH --job-name=euk-bot-mmseqsdb
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=euk-bot-mmseqsdb_%j.log
#SBATCH --account=132704181195
#SBATCH --partition=medium

set -euo pipefail

REF_VERSION=${REF_VERSION:-v1}
DB=/scratch/group/hu-lab/databases
PREFIX=${DB}/euk-bot-ref/${REF_VERSION}/euk-bot-ref.${REF_VERSION}
OUT_DIR=${DB}/euk-bot-ref/${REF_VERSION}_mmseqs
ENVS=/scratch/group/hu-lab/.cache/euk-bot/conda_envs
IMAGES=/scratch/group/hu-lab/.cache/euk-bot/singularity_images

# MetaEuk 7.bba0d80: the pipeline's Singularity image if pulled, else its conda env
METAEUK_SIF=$(ls "${IMAGES}"/*metaeuk-7.bba0d80*.img 2>/dev/null | head -1 || true)
METAEUK_ENV=$(dirname "$(dirname "$(grep -l "metaeuk=7.bba0d80" "${ENVS}"/env-*/conda-meta/history 2>/dev/null | head -1)")")
if [[ -n ${METAEUK_SIF} ]]; then
    METAEUK=(singularity exec --bind /scratch "${METAEUK_SIF}" metaeuk)
elif [[ -x ${METAEUK_ENV}/bin/metaeuk ]]; then
    METAEUK=("${METAEUK_ENV}/bin/metaeuk")
else
    echo "No MetaEuk 7.bba0d80 image in ${IMAGES} or env in ${ENVS}; run prebuild_envs.sh first" >&2; exit 1
fi

[[ -e ${OUT_DIR} ]] && { echo "${OUT_DIR} already exists" >&2; exit 1; }
umask 002   # keep shared databases group-readable
mkdir -p "${OUT_DIR}"

echo "== $(date +%T) metaeuk createdb ${PREFIX}.faa.gz"
"${METAEUK[@]}" createdb "${PREFIX}.faa.gz" "${OUT_DIR}/euk-bot-ref.${REF_VERSION}"

{
    echo "euk-bot reference DB ${REF_VERSION}: MMseqs2 protein DB for MetaEuk"
    echo "built:   $(date)  by ${USER}  (SLURM job ${SLURM_JOB_ID:-none})"
    echo "from:    ${PREFIX}.faa.gz (see ${DB}/euk-bot-ref/${REF_VERSION}/MANIFEST.txt)"
    echo "euk-bot: $(git rev-parse --short HEAD)$(git diff --quiet || echo ' + uncommitted changes')"
    echo "metaeuk: 7.bba0d80 createdb"
    echo
    echo "== outputs"
    ls -lh "${OUT_DIR}"
} > "${OUT_DIR}/MANIFEST.txt"
chmod -R a-w "${OUT_DIR}"   # a finished DB is never edited in place
echo "== $(date +%T) done"
cat "${OUT_DIR}/MANIFEST.txt"
