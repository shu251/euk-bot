// One row per assembly: contig stats, read mapping rate and BUSCO scores.
process ASSEMBLY_QC_SUMMARY {
    label 'process_low'
    conda 'conda-forge::python=3.9'
    container 'quay.io/biocontainers/python:3.9--1'   // has ps, which Nextflow needs; slim images don't

    input:
    path(stats)
    path(maprates)
    path(busco_summaries)    // optional: pass [] to skip

    output:
    path('assembly_qc.tsv'), emit: table

    script:
    def busco = busco_summaries ? "--busco ${busco_summaries}" : ''
    """
    summarize_assembly_qc.py --stats ${stats} --maprate ${maprates} ${busco} --out assembly_qc.tsv
    """

    stub:
    """
    touch assembly_qc.tsv
    """
}
