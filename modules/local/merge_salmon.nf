process MERGE_SALMON {
    label 'process_low'
    conda 'conda-forge::python=3.9'
    container 'quay.io/biocontainers/python:3.9--1'   // has ps, which Nextflow needs; slim images don't

    input:
    path(quant_dirs)

    output:
    path('salmon_numreads.tsv'), emit: numreads
    path('salmon_tpm.tsv'),      emit: tpm

    script:
    """
    merge_salmon.py --prefix salmon ${quant_dirs}
    """

    stub:
    """
    touch salmon_numreads.tsv salmon_tpm.tsv
    """
}
