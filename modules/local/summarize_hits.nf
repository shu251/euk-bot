process SUMMARIZE_HITS {
    label 'process_medium'
    conda 'conda-forge::python=3.9'
    container 'quay.io/biocontainers/python:3.9--1'   // has ps, which Nextflow needs; slim images don't

    input:
    path(hit_tables)
    path(taxonomy)    // optional: pass [] to skip

    output:
    path('reference_hit_counts.tsv'), emit: counts

    script:
    def tax = taxonomy ? "--taxonomy ${taxonomy}" : ''
    """
    summarize_diamond_hits.py ${tax} --out reference_hit_counts.tsv ${hit_tables}
    """

    stub:
    """
    touch reference_hit_counts.tsv
    """
}
