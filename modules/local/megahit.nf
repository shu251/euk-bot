process MEGAHIT {
    tag "${meta.id}"
    label 'process_high'
    conda 'bioconda::megahit=1.2.9'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.contigs.fa"), emit: contigs

    script:
    def args = task.ext.args ?: ''
    """
    megahit \\
        -1 ${reads[0]} -2 ${reads[1]} \\
        -o megahit_out \\
        -t ${task.cpus} \\
        --memory ${task.memory.toBytes()} \\
        ${args}

    # Prefix contig IDs with the sample and drop MEGAHIT's header annotations
    awk -v s=${meta.id} '/^>/ { print ">" s "_" substr(\$1, 2); next } { print }' \\
        megahit_out/final.contigs.fa > ${meta.id}.contigs.fa
    rm -rf megahit_out
    """

    stub:
    """
    touch ${meta.id}.contigs.fa
    """
}
