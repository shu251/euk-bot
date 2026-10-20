#!/bin/bash
# Download the eggNOG-mapper 2.1.x data (emapperdb 5.0.2) for --eggnog_data_dir: annotation DB, taxa DB and
# DIAMOND DB. About 15 GB to download and about 55 GB once unpacked, but only a handful of files.
# Grace compute nodes have had no internet, so the job checks first and exits right away if it can't reach
# the server. In that case run the same script on a login node inside tmux: bash download_eggnog.sbatch
# Downloads resume, so rerunning after a failure picks up where it stopped.
#
# Usage (from the repo directory):
#   sbatch download_eggnog.sbatch [data_dir]     # default: /scratch/group/hu-lab/databases/eggnog
#
#SBATCH --job-name=eggnog-dl
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=eggnog-dl_%j.log
#SBATCH --account=132704181195
#SBATCH --partition=medium

set -euo pipefail

DATA_DIR=${1:-/scratch/group/hu-lab/databases/eggnog}
BASE_URL=http://eggnog5.embl.de/download/emapperdb-5.0.2

if ! curl -sfI -m 30 "${BASE_URL}/eggnog.taxa.tar.gz" > /dev/null; then
    echo "Can't reach ${BASE_URL} from $(hostname)." >&2
    echo "Run it on a login node instead: tmux new -s eggnog; bash download_eggnog.sbatch ${DATA_DIR}" >&2
    exit 1
fi

umask 002   # keep shared databases group-readable
mkdir -p "${DATA_DIR}"
cd "${DATA_DIR}"

# Each file: skip if already unpacked by an earlier run, otherwise download (resuming) and unpack
for f in eggnog.db eggnog_proteins.dmnd; do
    [[ -s ${f} ]] && { echo "== ${f} already present"; continue; }
    echo "== $(date +%T) downloading ${f}.gz"
    wget --continue --no-verbose --tries=5 "${BASE_URL}/${f}.gz"
    echo "== $(date +%T) unpacking ${f}.gz"
    gunzip -f "${f}.gz"
done

if [[ -s eggnog.taxa.db ]]; then
    echo "== eggnog.taxa.db already present"
else
    wget --continue --no-verbose --tries=5 "${BASE_URL}/eggnog.taxa.tar.gz"
    tar -xzf eggnog.taxa.tar.gz && rm eggnog.taxa.tar.gz
fi

echo "== $(date +%T) done"
ls -lh "${DATA_DIR}"
