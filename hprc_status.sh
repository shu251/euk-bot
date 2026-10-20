#!/bin/bash
# Collect everything needed to diagnose euk-bot on Grace into one file, to copy back and share.
# Read-only; safe to run on a login node. The PhyloDB gzip test reads 5.6 GB, so allow a few minutes.
#
# Usage (from the repo directory):
#   bash hprc_status.sh            # writes /scratch/user/$USER/hprc_status.txt
# Then on the Mac:
#   scp skhu@grace.hprc.tamu.edu:/scratch/user/skhu/hprc_status.txt ~/Desktop/Projects/euk-bot/fromhprc/

OUT=${1:-/scratch/user/${USER}/hprc_status.txt}
DB=/scratch/group/hu-lab/databases
ENVS=/scratch/group/hu-lab/.cache/euk-bot/conda_envs

{
    echo "== $(date)  $(hostname)"
    echo "== git"; git log --oneline -3; git status --short

    echo "== quota"; showquota 2>&1
    echo "== partitions (name, time limit, nodes, cpus, memory MB)"; sinfo -o "%P %l %D %c %m" 2>&1

    echo "== conda envs (built / with python)"
    for e in "${ENVS}"/env-*; do
        [[ -d $e ]] || continue
        spec=$(grep -h "# cmd:" "$e"/conda-meta/history 2>/dev/null | tail -1 | grep -o "[a-z-]*::[^ ]*" | tr '\n' ' ')
        echo "$(basename "$e")  ${spec:-<incomplete: no conda-meta history>}"
    done
    echo "== prebuild.log tail"; tail -25 prebuild.log 2>&1

    echo "== eggnog"; ls -lh "${DB}/eggnog" 2>&1

    echo "== phylodb"
    cat "${DB}/phylodb/DATABASE_INFO.txt" 2>&1
    ls -lh "${DB}/phylodb/diamond" 2>&1
    head -c 600 "${DB}/phylodb/prot-map.json"; echo
    gzip -t "${DB}/phylodb/reference.pep.fa.gz" 2>&1 && echo "gzip test: OK"

    echo "== eukprot"; ls -lh "${DB}/eukprot_v3" 2>&1
    first=$(find "${DB}/eukprot_v3" -name "*.fasta" 2>/dev/null | head -1)
    [[ -n $first ]] && { echo "$first"; grep -m 2 '>' "$first"; }
    head -3 "${DB}"/eukprot_v3/EukProt_included_data_sets*.txt 2>&1 | cut -c1-400

    echo "== raw reads: read length from first sample"
    r1=$(ls /scratch/group/hu-lab/data/frenemies-metatranscriptome/raw-reads/*_R1_001.fastq.gz 2>/dev/null | head -1)
    [[ -n $r1 ]] && zcat "$r1" | head -4000 | awk 'NR % 4 == 2 { n[length($0)]++ } END { for (l in n) print l, n[l] }' | sort -k2,2nr | head -3
} > "${OUT}" 2>&1

echo "Wrote ${OUT}"
