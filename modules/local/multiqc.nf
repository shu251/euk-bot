process MULTIQC {
    label 'process_low'
    conda 'bioconda::multiqc=1.25.1'
    container 'quay.io/biocontainers/multiqc:1.25.1--pyhdfd78af_0'

    tag "${name}"

    input:
    val(name)       // report name, e.g. raw / clean
    path(files)

    output:
    path("${name}_multiqc_report.html"),      emit: report
    path("${name}_multiqc_report_data"),      emit: data

    script:
    """
    multiqc --force --title "euk-bot: ${name} reads" --filename ${name}_multiqc_report.html .
    """

    stub:
    """
    touch ${name}_multiqc_report.html
    mkdir ${name}_multiqc_report_data
    """
}
