process FASTQC {
    tag "${meta.id}"
    label 'process_low'
    conda 'bioconda::fastqc=0.12.1'
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'

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
    touch ${reads[0].simpleName}_fastqc.zip ${reads[0].simpleName}_fastqc.html
    """
}
