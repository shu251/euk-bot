process MULTIQC {
    label 'process_low'
    conda 'bioconda::multiqc=1.25.1'

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
