// Contig-level table: each clustered contig's representative protein (best reference hit), its reference
// taxonomy and eggNOG annotation, EUKulele consensus taxonomy, the protein's own eggNOG annotation, and salmon
// NumReads per sample (bin/build_contig_table.py). Also a per-protein table with the same annotation columns.
process CONTIG_TABLE {
    label 'process_medium'
    conda 'conda-forge::python=3.9'
    container 'quay.io/biocontainers/python:3.9--1'   // has ps, which Nextflow needs; slim images don't

    input:
    path(numreads)
    path(annotation_files, stageAs: 'annotation/*')   // proteins, DIAMOND hits, eggNOG, EUKulele (by file name)
    path(taxonomy)                                     // optional: pass [] to skip
    path(reference_annotations)                        // optional: pass [] to skip

    output:
    path('contig_table.tsv'),  emit: contigs
    path('protein_table.tsv'), emit: proteins

    script:
    def tax = taxonomy ? "--taxonomy ${taxonomy}" : ''
    def ref = reference_annotations ? "--reference-annotations ${reference_annotations}" : ''
    """
    build_contig_table.py --numreads ${numreads} ${tax} ${ref} annotation/*
    """

    stub:
    """
    touch contig_table.tsv protein_table.tsv
    """
}
