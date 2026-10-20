process MULTIQC {
    label 'process_low'
    conda 'bioconda::multiqc=1.25.1'
    container 'quay.io/biocontainers/multiqc:1.25.1--pyhdfd78af_0'

    input:
    path(files)

    output:
    path('multiqc_report.html'), emit: report
    path('multiqc_data'),        emit: data

    script:
    """
    multiqc --force .
    """

    stub:
    """
    touch multiqc_report.html
    mkdir multiqc_data
    """
}
