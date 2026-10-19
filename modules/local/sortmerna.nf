// Removes rRNA reads. --paired_in sends both mates to "rRNA" if either one hits.
process SORTMERNA {
    tag "${meta.id}"
    label 'process_high'
    conda 'bioconda::sortmerna=4.3.6'

    input:
    tuple val(meta), path(reads)
    path(dbs)

    output:
    tuple val(meta), path("${meta.id}_non_rRNA_{1,2}.fq.gz"), emit: reads
    tuple val(meta), path("${meta.id}.sortmerna.log"),        emit: log

    script:
    def refs = [dbs].flatten().collect { "--ref ${it}" }.join(' ')
    """
    sortmerna ${refs} \\
        --reads ${reads[0]} --reads ${reads[1]} \\
        --workdir smr_work \\
        --threads ${task.cpus} \\
        --fastx --paired_in --out2 \\
        --aligned rRNA --other non_rRNA

    mv non_rRNA_fwd.f*q.gz ${meta.id}_non_rRNA_1.fq.gz
    mv non_rRNA_rev.f*q.gz ${meta.id}_non_rRNA_2.fq.gz
    mv rRNA.log ${meta.id}.sortmerna.log
    rm -rf smr_work rRNA_fwd.* rRNA_rev.*
    """

    stub:
    """
    touch ${meta.id}_non_rRNA_1.fq.gz ${meta.id}_non_rRNA_2.fq.gz ${meta.id}.sortmerna.log
    """
}
