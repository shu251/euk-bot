// Taxonomic/functional assignment of predicted proteins against the reference DB.
// Standard 12-column tabular output (outfmt 6 default), the format EUKulele reads. --max-target-seqs is
// task.ext.max_targets: 1 (best hit) by default, 100 when EUKulele uses these hits (conf/modules.config).
process DIAMOND_BLASTP {
    tag "${name}"
    label 'process_high'
    conda 'bioconda::diamond=2.1.10'
    container 'quay.io/biocontainers/diamond:2.1.10--h5ca1c30_3'

    input:
    tuple val(name), path(pep)    // name: ORF caller (transdecoder / metaeuk)
    path(db)

    output:
    tuple val(name), path("${name}.proteins_vs_reference.tsv"), emit: hits

    script:
    def args = task.ext.args ?: '--sensitive'
    def max_targets = task.ext.max_targets ?: 1
    """
    diamond blastp \\
        --db ${db} \\
        --query ${pep} \\
        --out ${name}.proteins_vs_reference.tsv \\
        --outfmt 6 \\
        --evalue ${params.diamond_evalue} \\
        --max-target-seqs ${max_targets} \\
        --threads ${task.cpus} \\
        ${args}
    """

    stub:
    """
    touch ${name}.proteins_vs_reference.tsv
    """
}
