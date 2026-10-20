// Gene-content completeness of an assembly (transcriptome mode), for each lineage in --busco_lineages.
// For a mixed community this is a community-level score: useful for comparing assemblers and
// co-assembly on the same reads, not as completeness of any one organism.
process BUSCO {
    tag "${contigs.name - '.contigs.fa'}"
    label 'process_high'
    conda 'bioconda::busco=6.1.0'
    container 'quay.io/biocontainers/busco:6.1.0--pyhdfd78af_2'

    input:
    tuple val(meta), path(contigs)
    path(busco_db)

    output:
    path("*.busco.txt"), emit: summaries

    script:
    def name = contigs.name - '.contigs.fa'
    def lineages = params.busco_lineages.tokenize(',')*.trim()
    """
    for lineage in ${lineages.join(' ')}; do
        busco -i ${contigs} -m transcriptome -l \$lineage \\
            --download_path ${busco_db} --offline \\
            -c ${task.cpus} -o busco_\$lineage
        cp busco_\$lineage/short_summary.specific.\$lineage.busco_\$lineage.txt ${name}.\$lineage.busco.txt
        rm -rf busco_\$lineage
    done
    """

    stub:
    def name = contigs.name - '.contigs.fa'
    def lineages = params.busco_lineages.tokenize(',')*.trim()
    """
    for lineage in ${lineages.join(' ')}; do
        echo "C:0.0%[S:0.0%,D:0.0%],F:0.0%,M:100.0%,n:0" > ${name}.\$lineage.busco.txt
    done
    """
}
