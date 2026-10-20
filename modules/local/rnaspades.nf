process RNASPADES {
    tag "${meta.id}"
    label 'process_high'
    label 'process_high_memory'
    conda 'bioconda::spades=4.0.0'
    container 'quay.io/biocontainers/spades:4.0.0--haf24da9_4'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.contigs.fa"), emit: contigs

    script:
    def args = task.ext.args ?: ''
    """
    rnaspades.py \\
        -1 ${reads[0]} -2 ${reads[1]} \\
        -o spades_out \\
        -t ${task.cpus} \\
        -m ${task.memory.toGiga()} \\
        ${args}

    # Prefix contig IDs with the sample so they stay unique after pooling
    awk -v s=${meta.id} '/^>/ { print ">" s "_" substr(\$1, 2); next } { print }' \\
        spades_out/transcripts.fasta > ${meta.id}.contigs.fa
    rm -rf spades_out
    """

    stub:
    """
    touch ${meta.id}.contigs.fa
    """
}
