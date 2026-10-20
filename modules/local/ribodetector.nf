// Removes rRNA reads with a deep-learning classifier (no reference database needed).
// -e rrna: a pair is called rRNA if either mate is, so as many rRNA reads as possible are removed.
process RIBODETECTOR {
    tag "${meta.id}"
    label 'process_medium'
    conda 'bioconda::ribodetector=0.3.0'
    container 'quay.io/biocontainers/ribodetector:0.3.0--pyhdfd78af_0'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_non_rRNA_{1,2}.fq.gz"), emit: reads
    tuple val(meta), path("${meta.id}.ribodetector.tsv"),     emit: stats
    tuple val(meta), path("${meta.id}.ribodetector.log"),     emit: log

    script:
    """
    ribodetector_cpu \\
        -t ${task.cpus} \\
        -l ${params.ribodetector_read_len} \\
        -i ${reads[0]} ${reads[1]} \\
        -e rrna \\
        --chunk_size ${params.ribodetector_chunk_size} \\
        --log ${meta.id}.ribodetector.log \\
        -o ${meta.id}_non_rRNA_1.fq.gz ${meta.id}_non_rRNA_2.fq.gz

    # Per-sample rRNA fraction (read pairs in vs. kept)
    n_in=\$(( \$(zcat ${reads[0]} | wc -l) / 4 ))
    n_kept=\$(( \$(zcat ${meta.id}_non_rRNA_1.fq.gz | wc -l) / 4 ))
    awk -v s=${meta.id} -v i=\$n_in -v k=\$n_kept 'BEGIN {
        print "sample\\tinput_pairs\\tnon_rrna_pairs\\tpct_rrna"
        printf "%s\\t%d\\t%d\\t%.2f\\n", s, i, k, (i > 0 ? 100 * (i - k) / i : 0)
    }' > ${meta.id}.ribodetector.tsv
    """

    stub:
    """
    touch ${meta.id}_non_rRNA_1.fq.gz ${meta.id}_non_rRNA_2.fq.gz ${meta.id}.ribodetector.log
    printf 'sample\\tinput_pairs\\tnon_rrna_pairs\\tpct_rrna\\n' > ${meta.id}.ribodetector.tsv
    """
}
