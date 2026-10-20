// Pools all assemblies (every assembler, per-sample and co-assemblies) and clusters them into one
// non-redundant contig set. --cov-mode 1 -c: a contig joins a cluster only if the representative covers
// that fraction of it, so fragments contained in a longer contig collapse into it (as CD-HIT-EST did).
process MMSEQS_LINCLUST {
    label 'process_high'
    conda 'bioconda::mmseqs2=18.8cc5c'
    container 'quay.io/biocontainers/mmseqs2:18.8cc5c--hd6d6fdc_0'

    input:
    path(contigs)

    output:
    path('clustered_contigs.fa'),           emit: fasta
    path('clustered_contigs_clusters.tsv'), emit: clusters

    script:
    def args = task.ext.args ?: ''
    """
    cat ${contigs} > pooled_contigs.fa
    mmseqs easy-linclust \\
        pooled_contigs.fa clust tmp \\
        --min-seq-id ${params.cluster_identity} \\
        -c ${params.cluster_coverage} --cov-mode 1 \\
        --threads ${task.cpus} \\
        ${args}

    mv clust_rep_seq.fasta clustered_contigs.fa
    mv clust_cluster.tsv clustered_contigs_clusters.tsv
    rm -rf pooled_contigs.fa clust_all_seqs.fasta tmp
    """

    stub:
    """
    touch clustered_contigs.fa clustered_contigs_clusters.tsv
    """
}
