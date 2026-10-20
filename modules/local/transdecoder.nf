process TRANSDECODER {
    label 'process_medium'
    conda 'bioconda::transdecoder=5.7.1'
    container 'quay.io/biocontainers/transdecoder:5.7.1--pl5321hdfd78af_2'

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
