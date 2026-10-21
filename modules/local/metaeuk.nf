// Homology-based protein prediction on the clustered contigs, run alongside TransDecoder: MetaEuk finds only
// proteins with a homolog in the reference (--metaeuk_db, an MMseqs2 DB from build_metaeuk_db.sbatch), but
// with more accurate gene structure. -e / --metaeuk-eval / --min-length as in the MetaEuk README example.
process METAEUK {
    label 'process_high'
    conda 'bioconda::metaeuk=7.bba0d80'
    container 'quay.io/biocontainers/metaeuk:7.bba0d80--pl5321hd6d6fdc_2'

    input:
    path(contigs)
    path(db_dir)    // folder holding one MMseqs2 protein DB

    output:
    path("${contigs.baseName}.metaeuk.fas"),            emit: pep
    path("${contigs.baseName}.metaeuk.codon.fas"),      emit: cds
    path("${contigs.baseName}.metaeuk.gff"),            emit: gff
    path("${contigs.baseName}.metaeuk.headersMap.tsv"), emit: headers

    script:
    def split_mem = Math.max(1, (task.memory.toGiga() * 0.7) as int)
    """
    db=\$(ls ${db_dir}/*.dbtype | grep -v -e '_h\\.dbtype\$' -e '\\.idx' | head -1)
    metaeuk easy-predict ${contigs} \${db%.dbtype} ${contigs.baseName}.metaeuk tmp \\
        --threads ${task.cpus} \\
        -e 100 --metaeuk-eval 0.0001 --min-length 40 \\
        --split-memory-limit ${split_mem}G
    rm -rf tmp
    """

    stub:
    """
    touch ${contigs.baseName}.metaeuk.fas ${contigs.baseName}.metaeuk.codon.fas \\
        ${contigs.baseName}.metaeuk.gff ${contigs.baseName}.metaeuk.headersMap.tsv
    """
}
