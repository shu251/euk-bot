// Contig statistics (count, total length, N50, GC...) for every assembly in one table.
process SEQKIT_STATS {
    label 'process_low'
    conda 'bioconda::seqkit=2.14.0'
    container 'quay.io/biocontainers/seqkit:2.14.0--hb192632_0'

    input:
    path(contigs)

    output:
    path('assembly_stats.tsv'), emit: stats

    script:
    """
    seqkit stats --all --tabular --threads ${task.cpus} ${contigs} > assembly_stats.tsv
    """

    stub:
    """
    printf 'file\\tnum_seqs\\tsum_len\\tN50\\n' > assembly_stats.tsv
    """
}
