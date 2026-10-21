#!/bin/bash
# eggNOG-mapper annotation (KEGG KO/pathways, Pfam, COG, GO, ...) of a built euk-bot reference DB, run as a
# SLURM array over chunks of the protein FASTA. Same eggNOG-mapper version and settings as the pipeline's
# EGGNOG_MAPPER step. The version folder is read-only, so the result goes in a sibling <version>_eggnog folder;
# the .dmnd is unchanged. Needs the eggNOG data from download_eggnog.sbatch and the eggNOG-mapper image or env
# from prebuild_envs.sh.
#
# Three steps, from the repo directory:
#   sbatch -c 2 --mem=8G annotate_reference_db.sbatch split   # chunk the FASTA; prints the array commands
#   sbatch --array=0 annotate_reference_db.sbatch run         # one chunk first, to check time and memory
#   sbatch --array=1-<last>%8 annotate_reference_db.sbatch run
#   sbatch -c 2 --mem=8G annotate_reference_db.sbatch merge   # once every chunk is done
# (the #SBATCH resources below are for the run step)
# A chunk that finished is skipped if its array task is rerun. REF_VERSION=v2 for another version.
#
#SBATCH --job-name=euk-bot-refnog
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=120G
#SBATCH --output=euk-bot-refnog_%A_%a.log
#SBATCH --account=132704181195
#SBATCH --partition=medium

set -euo pipefail

STEP=${1:?"usage: sbatch annotate_reference_db.sbatch split|run|merge"}
REF_VERSION=${REF_VERSION:-v1}
DB=/scratch/group/hu-lab/databases
PREFIX=${DB}/euk-bot-ref/${REF_VERSION}/euk-bot-ref.${REF_VERSION}
EGGNOG_DATA=${EGGNOG_DATA:-${DB}/eggnog}
OUT_DIR=${DB}/euk-bot-ref/${REF_VERSION}_eggnog
WORK=${USER_SCRATCH:-/scratch/user/${USER}}/euk-bot/refdb_eggnog/${REF_VERSION}
CHUNK_PROTEINS=${CHUNK_PROTEINS:-2000000}
ENVS=/scratch/group/hu-lab/.cache/euk-bot/conda_envs
IMAGES=/scratch/group/hu-lab/.cache/euk-bot/singularity_images

umask 002   # keep shared databases group-readable

case ${STEP} in
split)
    [[ -e ${WORK}/chunks ]] && { echo "${WORK}/chunks already exists; remove it to re-split" >&2; exit 1; }
    mkdir -p "${WORK}/chunks" "${WORK}/out"
    # build_reference_db.py writes each record on two lines (header, sequence)
    zcat "${PREFIX}.faa.gz" | split -l $((CHUNK_PROTEINS * 2)) -d -a 3 --additional-suffix=.faa - "${WORK}/chunks/chunk_"
    for f in "${WORK}"/chunks/chunk_*.faa; do
        [[ $(head -c1 "$f") == ">" ]] || { echo "$f does not start with a FASTA header" >&2; exit 1; }
    done
    n=$(find "${WORK}/chunks" -name 'chunk_*.faa' | wc -l)
    echo "${n} chunks of up to ${CHUNK_PROTEINS} proteins in ${WORK}/chunks"
    echo "next: sbatch --array=0 annotate_reference_db.sbatch run, then --array=1-$((n - 1))%8"
    ;;

run)
    CHUNK=$(printf "chunk_%03d" "${SLURM_ARRAY_TASK_ID:?run as an array job: sbatch --array=...}")
    [[ -e ${WORK}/out/${CHUNK}.done ]] && { echo "${CHUNK} already done"; exit 0; }
    for f in eggnog.db eggnog_proteins.dmnd eggnog.taxa.db; do
        [[ -s ${EGGNOG_DATA}/${f} ]] || { echo "${EGGNOG_DATA}/${f} missing; run download_eggnog.sbatch" >&2; exit 1; }
    done

    # eggNOG-mapper 2.1.12: the pipeline's Singularity image if pulled, else its conda env
    EMAPPER_SIF=$(ls "${IMAGES}"/*eggnog-mapper-2.1.12*.img 2>/dev/null | head -1 || true)
    EMAPPER_ENV=$(dirname "$(dirname "$(grep -l "eggnog-mapper=2.1.12" "${ENVS}"/env-*/conda-meta/history 2>/dev/null | head -1)")")
    if [[ -n ${EMAPPER_SIF} ]]; then
        EMAPPER=(singularity exec --bind /scratch "${EMAPPER_SIF}" emapper.py)
    elif [[ -x ${EMAPPER_ENV}/bin/emapper.py ]]; then
        export PATH=${EMAPPER_ENV}/bin:${PATH}   # emapper.py calls diamond from the same env
        EMAPPER=(emapper.py)
    else
        echo "No eggNOG-mapper 2.1.12 image in ${IMAGES} or env in ${ENVS}; run prebuild_envs.sh first" >&2; exit 1
    fi

    TMP=${TMPDIR:-${WORK}/tmp}/${CHUNK}
    mkdir -p "${TMP}"
    echo "== $(date +%T) ${CHUNK} on $(hostname)"
    # --dbmem loads eggnog.db into RAM (~45 GB), much faster for millions of proteins
    "${EMAPPER[@]}" \
        -i "${WORK}/chunks/${CHUNK}.faa" --itype proteins \
        --data_dir "${EGGNOG_DATA}" \
        -m diamond \
        --dbmem \
        --cpu "${SLURM_CPUS_PER_TASK}" \
        --temp_dir "${TMP}" \
        --output_dir "${WORK}/out" \
        -o "${CHUNK}" \
        --override
    touch "${WORK}/out/${CHUNK}.done"
    rm -rf "${TMP}"
    echo "== $(date +%T) done"
    ;;

merge)
    n_chunks=$(find "${WORK}/chunks" -name 'chunk_*.faa' | wc -l)
    n_done=$(find "${WORK}/out" -name 'chunk_*.done' | wc -l)
    (( n_done == n_chunks )) || { echo "only ${n_done} of ${n_chunks} chunks done" >&2; exit 1; }
    [[ -e ${OUT_DIR} ]] && { echo "${OUT_DIR} already exists" >&2; exit 1; }
    mkdir -p "${OUT_DIR}"
    ANNOT=${OUT_DIR}/euk-bot-ref.${REF_VERSION}.emapper.annotations.tsv.gz

    # One header ("#query ..."), then the data rows of every chunk; "##" comment lines dropped
    first=$(ls "${WORK}"/out/chunk_*.emapper.annotations | head -1)
    { grep -m1 '^#query' "${first}"; grep -hv '^#' "${WORK}"/out/chunk_*.emapper.annotations; } | gzip -c > "${ANNOT}"

    n_proteins=$(grep -c '^>' "${WORK}"/chunks/chunk_*.faa | awk -F: '{n += $NF} END {print n}')
    n_annotated=$(zcat "${ANNOT}" | tail -n +2 | wc -l)
    {
        echo "euk-bot reference DB ${REF_VERSION}: eggNOG-mapper annotations"
        echo "merged:  $(date)  by ${USER}  (SLURM job ${SLURM_JOB_ID:-none})"
        echo "euk-bot: $(git rev-parse --short HEAD)$(git diff --quiet || echo ' + uncommitted changes')"
        echo "eggNOG-mapper version and settings: see the header of any ${WORK}/out/chunk_*.emapper.annotations"
        sed -n 's/^## //p' "${first}" | head -8
        echo
        echo "eggNOG data:"
        ls -l "${EGGNOG_DATA}"
        echo
        echo "proteins: ${n_proteins} in ${n_chunks} chunks; with an eggNOG hit: ${n_annotated}"
        echo
        ls -lh "${OUT_DIR}"
    } > "${OUT_DIR}/MANIFEST.txt"
    chmod -R a-w "${OUT_DIR}"
    cat "${OUT_DIR}/MANIFEST.txt"
    echo "chunks and per-chunk outputs can be removed once checked: rm -r ${WORK}"
    ;;

*)
    echo "unknown step '${STEP}': split, run or merge" >&2; exit 1
    ;;
esac
