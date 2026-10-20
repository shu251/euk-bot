process SALMON_QUANT {
    tag "${meta.id}"
    label 'process_medium'
    conda 'bioconda::salmon=1.10.3'
    container 'quay.io/biocontainers/salmon:1.10.3--h45fbf2d_5'

    input:
    tuple val(meta), path(reads)
    path(index)

    output:
    tuple val(meta), path("${meta.id}"), emit: quant

    script:
    """
    salmon quant \\
        -i ${index} -l A \\
        -1 ${reads[0]} -2 ${reads[1]} \\
        --validateMappings \\
        -p ${task.cpus} \\
        -o ${meta.id}
    """

    stub:
    """
    mkdir ${meta.id}
    printf 'Name\\tLength\\tEffectiveLength\\tTPM\\tNumReads\\n' > ${meta.id}/quant.sf
    """
}
