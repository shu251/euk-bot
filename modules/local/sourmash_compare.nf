// Pairwise sample similarity (angular similarity of abundance-weighted k-mer sketches).
process SOURMASH_COMPARE {
    label 'process_low'
    conda 'bioconda::sourmash=4.9.4'
    container 'quay.io/biocontainers/sourmash:4.9.4--hdfd78af_0'

    input:
    path(sigs)

    output:
    path('sourmash_similarity.csv'), emit: csv

    script:
    """
    sourmash compare -k 31 --csv sourmash_similarity.csv ${sigs}
    """

    stub:
    """
    touch sourmash_similarity.csv
    """
}
