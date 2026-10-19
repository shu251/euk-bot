process TRANSDECODER {
    label 'process_medium'
    conda 'bioconda::transdecoder=5.7.1'

    input:
    path(fasta)

    output:
    path("${fasta}.transdecoder.pep"), emit: pep
    path("${fasta}.transdecoder.cds"), emit: cds
    path("${fasta}.transdecoder.gff3"), emit: gff

    script:
    """
    TransDecoder.LongOrfs -t ${fasta} -m ${params.transdecoder_min_len}
    TransDecoder.Predict  -t ${fasta} --single_best_only
    """

    stub:
    """
    touch ${fasta}.transdecoder.pep ${fasta}.transdecoder.cds ${fasta}.transdecoder.gff3
    """
}
