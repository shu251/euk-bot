process SALMON_INDEX {
    label 'process_medium'
    conda 'bioconda::salmon=1.10.3'
    container 'quay.io/biocontainers/salmon:1.10.3--h45fbf2d_5'

    input:
    path(fasta)

    output:
    path('salmon_index'), emit: index

    script:
    """
    salmon index -t ${fasta} -i salmon_index -p ${task.cpus}
    """

    stub:
    """
    mkdir salmon_index
    """
}
