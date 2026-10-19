process CAT_FASTQ {
    tag "${meta.id}"

    input:
    tuple val(meta), path(reads_1, stageAs: 'r1/*'), path(reads_2, stageAs: 'r2/*')

    output:
    tuple val(meta), path("${meta.id}_R{1,2}.merged.fastq.gz"), emit: reads

    script:
    // gzip members concatenate into a valid gzip stream, so no recompression needed
    """
    cat ${reads_1.join(' ')} > ${meta.id}_R1.merged.fastq.gz
    cat ${reads_2.join(' ')} > ${meta.id}_R2.merged.fastq.gz
    """

    stub:
    """
    touch ${meta.id}_R1.merged.fastq.gz ${meta.id}_R2.merged.fastq.gz
    """
}
