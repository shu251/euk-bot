#!/bin/bash
# EUKulele reference folder for a built euk-bot reference version, for consensus taxonomy of the assembly
# proteins (--eukulele_db). Writes EUKulele's tax-table.txt and prot-map.json (bin/make_eukulele_db.py: how PR2,
# PhyloDB and UniEuk lineages map onto EUKulele's ranks) into a read-only <version>_eukulele folder beside the
# version folder. The proteins and the DIAMOND DB are symlinks to the version folder, under the names EUKulele
# looks for (it then skips its own makedb); the pipeline gives EUKulele its own DIAMOND hits, so the DB is not
# searched a second time.
#
# Usage (from the repo directory):
#   sbatch build_eukulele_db.sbatch                  # v1
#   REF_VERSION=v2 sbatch build_eukulele_db.sbatch
# Then run the pipeline with:
#   --eukulele_db /scratch/group/hu-lab/databases/euk-bot-ref/v1_eukulele
#
#SBATCH --job-name=euk-bot-eukuleledb
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=euk-bot-eukuleledb_%j.log
#SBATCH --account=132704181195
#SBATCH --partition=medium

set -euo pipefail

REF_VERSION=${REF_VERSION:-v1}
DB=/scratch/group/hu-lab/databases
REF_DIR=${DB}/euk-bot-ref/${REF_VERSION}
NAME=euk-bot-ref.${REF_VERSION}
OUT_DIR=${DB}/euk-bot-ref/${REF_VERSION}_eukulele

module purge
module load Anaconda3/2026.07-1   # python3

[[ -e ${OUT_DIR} ]] && { echo "${OUT_DIR} already exists" >&2; exit 1; }
umask 002   # keep shared databases group-readable
mkdir -p "${OUT_DIR}/diamond"

echo "== $(date +%T) tax-table.txt + prot-map.json from ${REF_DIR}/${NAME}.taxonomy.tsv.gz"
python3 bin/make_eukulele_db.py "${REF_DIR}/${NAME}.taxonomy.tsv.gz" --outdir "${OUT_DIR}" | tee "${OUT_DIR}/rank_rules.txt"

# EUKulele finds the DIAMOND DB at diamond/<ref_fasta with leading/trailing '.', 'f', 'a' stripped>.dmnd;
# for euk-bot-ref.<version>.faa.gz nothing is stripped
ln -s "${REF_DIR}/${NAME}.faa.gz" "${OUT_DIR}/${NAME}.faa.gz"
ln -s "${REF_DIR}/${NAME}.dmnd" "${OUT_DIR}/diamond/${NAME}.faa.gz.dmnd"

{
    echo "euk-bot reference DB ${REF_VERSION}: EUKulele reference folder"
    echo "built:   $(date)  by ${USER}  (SLURM job ${SLURM_JOB_ID:-none})"
    echo "from:    ${REF_DIR}/${NAME}.taxonomy.tsv.gz (see ${REF_DIR}/MANIFEST.txt)"
    echo "euk-bot: $(git rev-parse --short HEAD)$(git diff --quiet || echo ' + uncommitted changes')"
    echo "ranks:   Domain Supergroup Division Class Order Family Genus Species (PR2 Subdivision dropped); see bin/make_eukulele_db.py"
    echo "EUKulele: --reference_dir ${OUT_DIR} --ref_fasta ${NAME}.faa.gz"
    echo
    echo "== rank rules"
    cat "${OUT_DIR}/rank_rules.txt"
    echo
    echo "== outputs"
    ls -lh "${OUT_DIR}" "${OUT_DIR}/diamond"
} > "${OUT_DIR}/MANIFEST.txt"
rm "${OUT_DIR}/rank_rules.txt"
chmod -R a-w "${OUT_DIR}"   # a finished DB is never edited in place
echo "== $(date +%T) done"
cat "${OUT_DIR}/MANIFEST.txt"
