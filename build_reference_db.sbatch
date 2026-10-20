#!/bin/bash
# Build euk-bot's combined reference DB: MarFERReT v1.1.1 + EukProt v3 + PhyloDB prokaryotes, with exact
# duplicates across sources removed (see bin/build_reference_db.py), then a DIAMOND database from it using
# the pipeline's own DIAMOND (Singularity image or conda env) so the DB format matches. Run prebuild_envs.sh first.
#
# Each build goes in its own version folder with a MANIFEST.txt (sources, checksums, counts, git commit),
# and the job refuses to overwrite an existing version. Bump REF_VERSION when any source changes.
#
# Usage (from the repo directory):
#   sbatch build_reference_db.sbatch                       # builds REF_VERSION below
#   REF_VERSION=v2 sbatch build_reference_db.sbatch
# Then run the pipeline with:
#   --reference_dmnd     /scratch/group/hu-lab/databases/euk-bot-ref/v1/euk-bot-ref.v1.dmnd
#   --reference_taxonomy /scratch/group/hu-lab/databases/euk-bot-ref/v1/euk-bot-ref.v1.taxonomy.tsv.gz
#
#SBATCH --job-name=euk-bot-refdb
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=180G
#SBATCH --output=euk-bot-refdb_%j.log
#SBATCH --account=132704181195
#SBATCH --partition=medium

set -euo pipefail

REF_VERSION=${REF_VERSION:-v1}   # v1: MarFERReT v1.1.1 + EukProt v3 + PhyloDB 1.076 prokaryotes
DB=/scratch/group/hu-lab/databases
OUT_DIR=${DB}/euk-bot-ref/${REF_VERSION}
PREFIX=${OUT_DIR}/euk-bot-ref.${REF_VERSION}
MANIFEST=${OUT_DIR}/MANIFEST.txt
PHYLODB_FAA=${PHYLODB_FAA:-${DB}/phylodb_1.076/phylodb_1.076.pep.fa.gz}   # fresh copy; phylodb/reference.pep.fa.gz is corrupt
SOURCES=(
    "${DB}/marferret/data/MarFERReT.v1.1.1.proteins.faa.gz"
    "${DB}/marferret/data/MarFERReT.v1.1.1.metadata.csv"
    "${DB}/marferret/data/MarFERReT.v1.1.1.proteins_info.tab.gz"
    "${DB}/marferret/data/MarFERReT.v1.1.1.best_pfam_annotations.csv.gz"
    "${DB}/eukprot_v3/EukProt_proteins.v03.2021_11_22.tgz"
    "${DB}/eukprot_v3/EukProt_included_data_sets.v03.2021_11_22.txt"
    "${PHYLODB_FAA}"
    "${DB}/phylodb/prot-map.json"
    "${DB}/phylodb/taxonomy-table.txt"
    "${DB}/PR2/pr2-tax.tsv"
)
ENVS=/scratch/group/hu-lab/.cache/euk-bot/conda_envs
IMAGES=/scratch/group/hu-lab/.cache/euk-bot/singularity_images

module purge
module load Anaconda3/2026.07-1   # python3 for the merge script

# DIAMOND 2.1.10: the pipeline's Singularity image if pulled, else its conda env
DIAMOND_SIF=$(ls "${IMAGES}"/*diamond-2.1.10*.img 2>/dev/null | head -1 || true)
DIAMOND_ENV=$(dirname "$(dirname "$(grep -l "diamond=2.1.10" "${ENVS}"/env-*/conda-meta/history 2>/dev/null | head -1)")")
if [[ -n ${DIAMOND_SIF} ]]; then
    DIAMOND=(singularity exec "${DIAMOND_SIF}" diamond)
elif [[ -x ${DIAMOND_ENV}/bin/diamond ]]; then
    DIAMOND=("${DIAMOND_ENV}/bin/diamond")
else
    echo "No DIAMOND 2.1.10 image in ${IMAGES} or env in ${ENVS}; run prebuild_envs.sh first" >&2; exit 1
fi
gzip -t "${PHYLODB_FAA}" || { echo "${PHYLODB_FAA} is not a valid gzip file" >&2; exit 1; }

[[ -e ${OUT_DIR} ]] && { echo "${OUT_DIR} already exists; set a new REF_VERSION" >&2; exit 1; }
umask 002   # keep shared databases group-readable
mkdir -p "${OUT_DIR}"

{
    echo "euk-bot reference DB ${REF_VERSION}"
    echo "built:   $(date)  by ${USER}  (SLURM job ${SLURM_JOB_ID:-none})"
    echo "euk-bot: $(git rev-parse --short HEAD)$(git diff --quiet || echo ' + uncommitted changes')"
    echo "diamond: $("${DIAMOND[@]}" version)"
    echo
    echo "== sources (md5, size, path)"
    for f in "${SOURCES[@]}"; do
        echo "$(md5sum < "$f" | cut -d' ' -f1)  $(du -h "$f" | cut -f1)  $f"
    done
    echo
    echo "== merge counts"
} > "${MANIFEST}"

echo "== $(date +%T) merging sources"
python3 bin/build_reference_db.py \
    --marferret "${DB}/marferret/data" \
    --eukprot "${DB}/eukprot_v3" \
    --phylodb "${DB}/phylodb" \
    --phylodb-faa "${PHYLODB_FAA}" \
    --pr2-taxonomy "${DB}/PR2/pr2-tax.tsv" \
    --out "${PREFIX}" 2> >(tee -a "${MANIFEST}" >&2)

echo "== $(date +%T) diamond makedb ($("${DIAMOND[@]}" version))"
"${DIAMOND[@]}" makedb \
    --in "${PREFIX}.faa.gz" \
    --db "${PREFIX}" \
    --threads "${SLURM_CPUS_PER_TASK}"

{ echo; echo "== outputs"; ls -lh "${OUT_DIR}"; } >> "${MANIFEST}"
chmod -R a-w "${OUT_DIR}"   # a finished version is never edited in place
echo "== $(date +%T) done"
cat "${MANIFEST}"
