#!/bin/bash
# Stand-alone test, not part of the pipeline: which BBTools version runs BBNorm on paired gzipped reads without
# the "AssertionError: List size mismatch" crash seen in test3 (bbmap 40.02; worker threads die and BBNorm then
# hangs until the walltime). The input pairs are in sync (checked read names: 0 mismatched), so the reader is suspect.
# Each version runs the pipeline's BBNorm command on one sample's RiboDetector output for TEST_MIN minutes. The
# crash came within seconds, so a version still running at the timeout with no AssertionError passes.
#
# Images are pulled into the pipeline's Singularity cache under Nextflow's file names, so the version that passes
# is already cached when the BBNORM container is switched to it.
#
# Usage (from the repo directory; the work dir is the RIBODETECTOR task's, e.g. from the trace file hash):
#   sbatch benchmarks/bbnorm_version_test.sbatch /scratch/user/skhu/euk-bot/work/51/ead619*   # sample 140
#
#SBATCH --job-name=bbnorm-test
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
#SBATCH --output=bbnorm-test_%j.log
#SBATCH --account=132704181195
#SBATCH --partition=short

set -euo pipefail

RD_WORK=${1:?"usage: sbatch benchmarks/bbnorm_version_test.sbatch <work dir of a RIBODETECTOR task>"}
TEST_MIN=${TEST_MIN:-10}
VERSIONS=(40.02--h09cc210_0 39.13--he5f24ec_1 40.05--h09cc210_0)

USER_SCRATCH=${USER_SCRATCH:-/scratch/user/${USER}}
GROUP_SPACE=${GROUP_SPACE:-/scratch/group/hu-lab}
IMAGES=${GROUP_SPACE}/.cache/euk-bot/singularity_images
export SINGULARITY_CACHEDIR=${TMPDIR:-${USER_SCRATCH}/.cache}/.singularity   # HPRC: local /tmp, never $HOME
export SINGULARITY_TMPDIR=${TMPDIR:-${USER_SCRATCH}/.cache}
umask 002
mkdir -p "${SINGULARITY_CACHEDIR}"

R1=$(ls "${RD_WORK}"/*_non_rRNA_1.fq.gz)
R2=$(ls "${RD_WORK}"/*_non_rRNA_2.fq.gz)
SAMPLE=$(basename "${R1}" _non_rRNA_1.fq.gz)
OUT=${TMPDIR:-${USER_SCRATCH}/euk-bot/bbnorm_test}/bbnorm_test
mkdir -p "${OUT}"

module purge
module load WebProxy

echo "== $(date +%T) ${SAMPLE} on $(hostname); ${TEST_MIN} min per version"
xmx=$(( SLURM_MEM_PER_NODE * 85 / 100 / 1024 ))   # same 85% of task memory as the BBNORM module
for v in "${VERSIONS[@]}"; do
    img=${IMAGES}/quay.io-biocontainers-bbmap-${v}.img
    [[ -s ${img} ]] || singularity pull "${img}" "docker://quay.io/biocontainers/bbmap:${v}"

    echo "== $(date +%T) bbmap ${v}"
    rc=0
    timeout "${TEST_MIN}m" singularity exec --bind /scratch "${img}" \
        bbnorm.sh -Xmx${xmx}g \
            in="${R1}" in2="${R2}" \
            out="${OUT}/norm_1.fq.gz" out2="${OUT}/norm_2.fq.gz" \
            target=100 mindepth=2 prefilter=t tmpdir="${OUT}" overwrite=t \
            threads="${SLURM_CPUS_PER_TASK}" \
        2> "${OUT}/bbnorm_${v}.log" || rc=$?
    n_err=$(grep -c AssertionError "${OUT}/bbnorm_${v}.log" || true)
    if (( rc == 124 && n_err == 0 )); then verdict=PASS
    elif (( rc == 0 && n_err == 0 )); then verdict="PASS (finished)"
    else verdict=FAIL; fi
    echo "-- ${v}: exit ${rc} (124 = still running at ${TEST_MIN} min), AssertionErrors: ${n_err} -> ${verdict}"
    grep -m 3 -E 'AssertionError|Exception|Made prefilter|Made hash table' "${OUT}/bbnorm_${v}.log" || true
    rm -f "${OUT}"/norm_*.fq.gz
done
echo "== $(date +%T) done"
