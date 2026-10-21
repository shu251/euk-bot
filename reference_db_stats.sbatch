#!/bin/bash
# Summary tables for a built euk-bot reference DB (sizes, length distribution, per-organism taxonomy, Pfam
# coverage, and KEGG/eggNOG coverage once annotate_reference_db.sbatch has merged its output); see
# bin/reference_db_stats.py for the outputs. The version folder is read-only, so the tables go
# in the repo next to the manifest copy, small enough to keep in git.
#
# Usage (from the repo directory):
#   sbatch reference_db_stats.sbatch                    # v1
#   REF_VERSION=v2 sbatch reference_db_stats.sbatch
#
#SBATCH --job-name=euk-bot-refstats
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G   # eggNOG annotations are held in memory by protein ID
#SBATCH --output=euk-bot-refstats_%j.log
#SBATCH --account=132704181195
#SBATCH --partition=medium

set -euo pipefail

REF_VERSION=${REF_VERSION:-v1}
REF_DIR=/scratch/group/hu-lab/databases/euk-bot-ref/${REF_VERSION}
PREFIX=${REF_DIR}/euk-bot-ref.${REF_VERSION}
OUT_DIR=${OUT_DIR:-assets/reference_db/euk-bot-ref.${REF_VERSION}.stats}
EGGNOG=/scratch/group/hu-lab/databases/euk-bot-ref/${REF_VERSION}_eggnog/euk-bot-ref.${REF_VERSION}.emapper.annotations.tsv.gz

module purge
module load Anaconda3/2026.07-1   # python3

mkdir -p "${OUT_DIR}"
echo "== $(date +%T) reading ${PREFIX}.faa.gz + taxonomy"
EGGNOG_ARGS=()
if [[ -s ${EGGNOG} ]]; then
    EGGNOG_ARGS=(--eggnog "${EGGNOG}")
else
    echo "no eggNOG annotations at ${EGGNOG}; KEGG tables skipped"
fi
python3 bin/reference_db_stats.py \
    --faa "${PREFIX}.faa.gz" \
    --taxonomy "${PREFIX}.taxonomy.tsv.gz" \
    "${EGGNOG_ARGS[@]}" \
    --outdir "${OUT_DIR}"

# On-disk size of each file in the version folder
{ printf "file\tbytes\n"; for f in "${REF_DIR}"/*; do printf "%s\t%s\n" "$(basename "$f")" "$(stat -c %s "$f")"; done; } \
    > "${OUT_DIR}/files.tsv"

echo "== $(date +%T) done"
ls -lh "${OUT_DIR}"
cat "${OUT_DIR}/size_by_source.tsv"
