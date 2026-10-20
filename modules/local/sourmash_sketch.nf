// k-mer sketch of each sample's reads (abundance-weighted), for read-based sample similarity.
process SOURMASH_SKETCH {
    tag "${meta.id}"
    label 'process_low'
    conda 'bioconda::sourmash=4.9.4'
    container 'quay.io/biocontainers/sourmash:4.9.4--hdfd78af_0'

    input:
    tuple val(meta), path(reads)

    output:
    path("${meta.id}.sig"), emit: sig

    script:
    """
    sourmash sketch dna -p k=31,scaled=1000,abund --merge ${meta.id} -o ${meta.id}.sig ${reads}
    """

    stub:
    """
    touch ${meta.id}.sig
    """
}
