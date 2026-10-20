// Library strandedness from a subsample of R1 reads: DIAMOND blastx reports the frame each read hit in,
// so sense hits (+ frames) = forward-stranded (fr), antisense (- frames) = reverse-stranded (rf, e.g.
// dUTP), a mix = unstranded.
process STRANDEDNESS {
    tag "${meta.id}"
    label 'process_medium'
    conda 'bioconda::diamond=2.1.10'
    container 'quay.io/biocontainers/diamond:2.1.10--h5ca1c30_3'

    input:
    tuple val(meta), path(reads)
    path(db)

    output:
    path("${meta.id}.strandedness.tsv"), emit: summary

    script:
    def n_reads = params.strand_check_reads
    """
    zcat ${reads[0]} | head -n ${n_reads * 4} > r1_subsample.fq
    diamond blastx --db ${db} --query r1_subsample.fq --out hits.tsv \\
        --outfmt 6 qseqid qframe --max-target-seqs 1 --evalue 1e-5 --fast --threads ${task.cpus}
    awk -v s=${meta.id} -v n=${n_reads} '
        !seen[\$1]++ { if (\$2 > 0) fwd++; else rev++ }
        END {
            hits = fwd + rev; f = hits ? fwd / hits : 0
            call = hits < 1000 ? "too_few_hits" : f >= 0.8 ? "fr" : f <= 0.2 ? "rf" : "unstranded"
            print "sample\\treads_checked\\treads_with_hits\\tfrac_r1_sense\\tstrandedness"
            printf "%s\\t%d\\t%d\\t%.3f\\t%s\\n", s, n, hits, f, call
        }' hits.tsv > ${meta.id}.strandedness.tsv
    rm r1_subsample.fq hits.tsv
    """

    stub:
    """
    printf 'sample\\tstrandedness\\n${meta.id}\\tNA\\n' > ${meta.id}.strandedness.tsv
    """
}
