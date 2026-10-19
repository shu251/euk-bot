process FASTQC {
    tag "${meta.id}"
    label 'process_low'
    conda 'bioconda::fastqc=0.12.1'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*_fastqc.zip'),  emit: zip
    tuple val(meta), path('*_fastqc.html'), emit: html

    script:
    """
    fastqc --threads ${task.cpus} ${reads}
    """

    stub:
    """
    touch ${meta.id}_1_fastqc.zip ${meta.id}_1_fastqc.html
    """
}
