// Pools all per-sample assemblies and clusters them into one non-redundant contig set.
process CDHIT_EST {
    label 'process_high'
    label 'process_high_memory'
    conda 'bioconda::cd-hit=4.8.1'
    container 'quay.io/biocontainers/cd-hit:4.8.1--h5ca1c30_13'

    input:
    path(contigs)

    output:
    path('clustered_contigs.fa'),         emit: fasta
    path('clustered_contigs.fa.clstr'),   emit: clstr

    script:
    """
    cat ${contigs} > pooled_contigs.fa
    cd-hit-est \\
        -i pooled_contigs.fa \\
        -o clustered_contigs.fa \\
        -c ${params.cdhit_identity} \\
        -n 10 -d 0 \\
        -M ${task.memory.toMega()} \\
        -T ${task.cpus}
    rm pooled_contigs.fa
    """

    stub:
    """
    touch clustered_contigs.fa clustered_contigs.fa.clstr
    """
}
