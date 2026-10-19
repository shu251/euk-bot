// Taxonomic/functional assignment of predicted proteins against the reference DB.
process DIAMOND_BLASTP {
    label 'process_high'
    conda 'bioconda::diamond=2.1.10'

    input:
    path(pep)
    path(db)

    output:
    path('proteins_vs_reference.tsv'), emit: hits

    script:
    def args = task.ext.args ?: '--sensitive'
    """
    diamond blastp \\
        --db ${db} \\
        --query ${pep} \\
        --out proteins_vs_reference.tsv \\
        --outfmt 6 qseqid sseqid pident length evalue bitscore \\
        --evalue ${params.diamond_evalue} \\
        --max-target-seqs 1 \\
        --threads ${task.cpus} \\
        ${args}
    """

    stub:
    """
    touch proteins_vs_reference.tsv
    """
}
