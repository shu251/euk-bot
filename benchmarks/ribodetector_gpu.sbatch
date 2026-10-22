#!/bin/bash
# Stand-alone benchmark, not part of the pipeline: RiboDetector GPU mode on one sample, compared with the
# pipeline's CPU run (RIBODETECTOR, 12 cores, ~6.3 h and ~62 CPU-h per sample in test3).
# Input is the same BBDuk output the CPU task used, read from that task's work dir, with the same settings
# (ribodetector 0.3.0, -l 150, -e rrna, gzipped output) so only the device differs.
#
# The first run builds a CUDA conda env on personal scratch (~20 min, once). It is built on the GPU node so
# the CUDA build matches the node's driver.
#
# Usage (from the repo directory; the work dir is the CPU RIBODETECTOR task's, e.g. from the trace file hash):
#   sbatch benchmarks/ribodetector_gpu.sbatch /scratch/user/skhu/euk-bot/work/76/aa8a92*   # sample 137
#
#SBATCH --job-name=ribodetector-gpu
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=120G
#SBATCH --gres=gpu:1                  # any GPU type; e.g. gpu:a100:1 to pick one
#SBATCH --output=ribodetector-gpu_%j.log
#SBATCH --account=132704181195
#SBATCH --partition=gpu

set -euo pipefail

CPU_WORK=${1:?"usage: sbatch benchmarks/ribodetector_gpu.sbatch <work dir of the CPU RIBODETECTOR task>"}
READ_LEN=${READ_LEN:-150}       # params.ribodetector_read_len
CHUNK_SIZE=${CHUNK_SIZE:-256}

USER_SCRATCH=${USER_SCRATCH:-/scratch/user/${USER}}
ENV=${USER_SCRATCH}/.cache/euk-bot/ribodetector_gpu_env
export CONDA_PKGS_DIRS=${USER_SCRATCH}/.cache/conda_pkgs

R1=$(ls "${CPU_WORK}"/*_bbduk_1.fq.gz)
R2=$(ls "${CPU_WORK}"/*_bbduk_2.fq.gz)
SAMPLE=$(basename "${R1}" _bbduk_1.fq.gz)
OUT=${USER_SCRATCH}/euk-bot/ribodetector_gpu_test/${SAMPLE}
mkdir -p "${OUT}"

module purge
module load Anaconda3/2026.07-1
module load WebProxy

echo "== $(date +%T) ${SAMPLE} on $(hostname)"
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv
CUDA_DRIVER=$(nvidia-smi | grep -oP 'CUDA Version: \K[0-9.]+')
GPU_GB=$(nvidia-smi --query-gpu=memory.total --format=csv,noheader,nounits | head -1 | awk '{print int($1 / 1024)}')
GPU_MEM=$(( GPU_GB * 3 / 4 ))   # -m: GPU RAM RiboDetector may use; leave headroom

if [[ ! -x ${ENV}/bin/ribodetector ]]; then
    echo "== $(date +%T) building ${ENV} (CUDA ${CUDA_DRIVER})"
    CONDA_OVERRIDE_CUDA=${CUDA_DRIVER} conda create -y -p "${ENV}" --override-channels -c conda-forge -c bioconda \
        ribodetector=0.3.0 'pytorch=*=cuda*'
fi
"${ENV}/bin/python" -c 'import torch; assert torch.cuda.is_available(), "no CUDA in torch"; print("torch", torch.__version__, torch.cuda.get_device_name(0))'

echo "== $(date +%T) ribodetector GPU: -t ${SLURM_CPUS_PER_TASK} -m ${GPU_MEM} --chunk_size ${CHUNK_SIZE}"
start=$(date +%s)
"${ENV}/bin/ribodetector" \
    -t "${SLURM_CPUS_PER_TASK}" \
    -l "${READ_LEN}" \
    -i "${R1}" "${R2}" \
    -m "${GPU_MEM}" \
    -e rrna \
    --chunk_size "${CHUNK_SIZE}" \
    --log "${OUT}/${SAMPLE}.ribodetector_gpu.log" \
    -o "${OUT}/${SAMPLE}_non_rRNA_1.fq.gz" "${OUT}/${SAMPLE}_non_rRNA_2.fq.gz" \
    -r "${OUT}/${SAMPLE}.rRNA_1.fq.gz" "${OUT}/${SAMPLE}.rRNA_2.fq.gz"
elapsed=$(( $(date +%s) - start ))

# Same pair counts as the pipeline's <sample>.ribodetector.tsv, so the two runs line up
n_in=$(( $(zcat "${R1}" | wc -l) / 4 ))
n_kept=$(( $(zcat "${OUT}/${SAMPLE}_non_rRNA_1.fq.gz" | wc -l) / 4 ))
echo "== $(date +%T) done"
echo "-- GPU run"
awk -v s="${SAMPLE}" -v i="${n_in}" -v k="${n_kept}" -v t="${elapsed}" 'BEGIN {
    print "sample\tinput_pairs\tnon_rrna_pairs\tpct_rrna\tribodetector_min"
    printf "%s\t%d\t%d\t%.2f\t%.1f\n", s, i, k, (i > 0 ? 100 * (i - k) / i : 0), t / 60
}' | tee "${OUT}/${SAMPLE}.ribodetector_gpu.tsv"
echo "-- CPU run (pipeline)"
cat "${CPU_WORK}/${SAMPLE}.ribodetector.tsv"
echo "outputs in ${OUT}; remove when done: rm -r ${USER_SCRATCH}/euk-bot/ribodetector_gpu_test"
