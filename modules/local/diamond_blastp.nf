// Taxonomic/functional assignment of predicted proteins against the reference DB.
process DIAMOND_BLASTP {
    tag "${name}"
    label 'process_high'
    conda 'bioconda::diamond=2.1.10'
    container 'quay.io/biocontainers/diamond:2.1.10--h5ca1c30_3'

    input:
    tuple val(name), path(pep)    // name: ORF caller (transdecoder / metaeuk)
    path(db)

    output:
    path("${name}.proteins_vs_reference.tsv"), emit: hits

    script:
    def args = task.ext.args ?: '--sensitive'
    """
    diamond blastp \\
        --db ${db} \\
        --query ${pep} \\
        --out ${name}.proteins_vs_reference.tsv \\
        --outfmt 6 qseqid sseqid pident length evalue bitscore \\
        --evalue ${params.diamond_evalue} \\
        --max-target-seqs 1 \\
        --threads ${task.cpus} \\
        ${args}
    """

    stub:
    """
    touch ${name}.proteins_vs_reference.tsv
    """
}
