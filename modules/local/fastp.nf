// Adapter trimming (auto-detected for pairs), Trimmomatic-style sliding-window quality trimming (window 4,
// mean quality --trim_mean_quality, like SLIDINGWINDOW:4:15) and a minimum length (--min_read_length).
process FASTP {
    tag "${meta.id}"
    label 'process_medium'
    conda 'bioconda::fastp=0.23.4'
    container 'quay.io/biocontainers/fastp:0.23.4--h125f33a_5'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_trimmed_{1,2}.fq.gz"), emit: reads
    tuple val(meta), path("${meta.id}.fastp.json"),          emit: json
    tuple val(meta), path("${meta.id}.fastp.html"),          emit: html

    script:
    def args = task.ext.args ?: ''
    """
    fastp \\
        -i ${reads[0]} -I ${reads[1]} \\
        -o ${meta.id}_trimmed_1.fq.gz -O ${meta.id}_trimmed_2.fq.gz \\
        --detect_adapter_for_pe \\
        --cut_right --cut_window_size 4 --cut_mean_quality ${params.trim_mean_quality} \\
        --length_required ${params.min_read_length} \\
        --thread ${Math.min(task.cpus as int, 16)} \\
        --json ${meta.id}.fastp.json \\
        --html ${meta.id}.fastp.html \\
        ${args}
    """

    stub:
    """
    touch ${meta.id}_trimmed_1.fq.gz ${meta.id}_trimmed_2.fq.gz ${meta.id}.fastp.json ${meta.id}.fastp.html
    """
}
