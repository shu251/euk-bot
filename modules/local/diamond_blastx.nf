// Both mates go into one query; they share a read ID, so SUMMARIZE_HITS counts each pair once.
process DIAMOND_BLASTX {
    tag "${meta.id}"
    label 'process_high'
    conda 'bioconda::diamond=2.1.10'

    input:
    tuple val(meta), path(reads)
    path(db)

    output:
    tuple val(meta), path("${meta.id}.blastx.tsv"), emit: hits

    script:
    def args = task.ext.args ?: '--fast'
    """
    cat ${reads[0]} ${reads[1]} > query.fq.gz
    diamond blastx \\
        --db ${db} \\
        --query query.fq.gz \\
        --out ${meta.id}.blastx.tsv \\
        --outfmt 6 qseqid sseqid pident length evalue bitscore \\
        --evalue ${params.diamond_evalue} \\
        --max-target-seqs 1 \\
        --threads ${task.cpus} \\
        ${args}
    rm query.fq.gz
    """

    stub:
    """
    touch ${meta.id}.blastx.tsv
    """
}
