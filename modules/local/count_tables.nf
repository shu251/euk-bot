// Count tables (bin/build_count_tables.py): per-contig reads (+ CPM) with the contig's EUKulele consensus lineage
// and KEGG KO, plus long tables summed by KO, by lineage (all ranks and cut at --count_rank) and by lineage x KO.
// Library sizes come from the salmon quant dirs (NumReads summed = mapped; num_processed = input).
process COUNT_TABLES {
    label 'process_medium'
    conda 'conda-forge::python=3.9'
    container 'quay.io/biocontainers/python:3.9--1'   // has ps, which Nextflow needs; slim images don't

    input:
    path(contig_table)
    path(quant_dirs, stageAs: 'salmon/*')

    output:
    path('library_sizes.tsv'),   emit: library_sizes
    path('contig_counts.tsv'),   emit: contigs
    path('ko_counts.tsv'),       emit: ko
    path('taxon_counts*.tsv'),   emit: taxon
    path('taxon_ko_counts*.tsv'), emit: taxon_ko

    script:
    def args = task.ext.args ?: ''   // --rank / --libsize from params (conf/modules.config)
    """
    build_count_tables.py --contig-table ${contig_table} --quant-dirs ${quant_dirs} ${args}
    """

    stub:
    """
    touch library_sizes.tsv contig_counts.tsv ko_counts.tsv taxon_counts.tsv taxon_ko_counts.tsv \\
        taxon_counts.${params.count_rank}.tsv taxon_ko_counts.${params.count_rank}.tsv
    """
}
