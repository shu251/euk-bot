#!/bin/bash
# Download a prebuilt Kaiju index for --kaiju_db (read profiling). Default is nr_euk 2023-05-10: NCBI nr
# archaea, bacteria, viruses plus fungi and microbial eukaryotes, the newest prebuilt index that includes
# protists. About 82 GB to download; Kaiju needs about 204 GB RAM to load it.
# Each index goes in its own folder named after the tarball, so other versions can sit beside it.
# Compute nodes reach the internet through the WebProxy module. The download resumes, so rerunning after
# a failure picks up where it stopped.
#
# Usage (from the repo directory):
#   sbatch download_kaiju.sbatch                    # nr_euk 2023-05-10
#   sbatch download_kaiju.sbatch <tarball URL>      # any index from bioinformatics-centre.github.io/kaiju/downloads.html
# Then run the pipeline with:
#   --kaiju_db /scratch/group/hu-lab/databases/kaiju/kaiju_db_nr_euk_2023-05-10
#
#SBATCH --job-name=kaiju-dl
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=kaiju-dl_%j.log
#SBATCH --account=132704181195
#SBATCH --partition=medium

set -euo pipefail

module load WebProxy   # internet access from compute nodes

URL=${1:-https://kaiju-idx.s3.eu-central-1.amazonaws.com/2023/kaiju_db_nr_euk_2023-05-10.tgz}
TARBALL=$(basename "${URL}")
DATA_DIR=${DATA_DIR:-/scratch/group/hu-lab/databases/kaiju}/${TARBALL%.tgz}

if ! curl -sfI -m 30 "${URL}" > /dev/null; then
    echo "Can't reach ${URL} from $(hostname)." >&2
    exit 1
fi

umask 002   # keep shared databases group-readable
mkdir -p "${DATA_DIR}"
cd "${DATA_DIR}"

if ls ./*.fmi > /dev/null 2>&1 && [[ -s nodes.dmp && -s names.dmp ]]; then
    echo "== index already unpacked in ${DATA_DIR}"
else
    echo "== $(date +%T) downloading ${TARBALL}"
    wget --continue --no-verbose --tries=5 "${URL}"
    echo "== $(date +%T) unpacking"
    tar -xzf "${TARBALL}" && rm "${TARBALL}"
    ls ./*.fmi > /dev/null 2>&1 && [[ -s nodes.dmp && -s names.dmp ]] || {
        echo "Unpacked, but no *.fmi + nodes.dmp + names.dmp at the top of ${DATA_DIR}:" >&2; ls -lR >&2; exit 1; }
    {
        echo "Kaiju index ${TARBALL%.tgz}"
        echo "source:     ${URL}"
        echo "downloaded: $(date)  by ${USER}  (SLURM job ${SLURM_JOB_ID:-none})"
    } > MANIFEST.txt
    chmod -R a-w "${DATA_DIR}"   # shared read-only copy
fi

echo "== $(date +%T) done"
ls -lh "${DATA_DIR}"
