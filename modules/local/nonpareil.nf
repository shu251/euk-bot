// Sequencing coverage and diversity from read redundancy (k-mer mode, R1 only). Built for metagenomes,
// so for metatranscriptomes read it as relative: higher coverage = sample sequenced more completely.
process NONPAREIL {
    tag "${meta.id}"
    label 'process_medium'
    conda 'bioconda::nonpareil=3.5.5'
    container 'quay.io/biocontainers/nonpareil:3.5.5--r44h077b44d_2'

    input:
    tuple val(meta), path(reads)

    output:
    path("${meta.id}.npo"),           emit: npo
    path("${meta.id}.nonpareil.tsv"), emit: summary

    script:
    """
    zcat ${reads[0]} > r1.fq
    nonpareil -s r1.fq -T kmer -f fastq -b ${meta.id} -t ${task.cpus} -R ${task.memory.toMega()}
    rm r1.fq

    Rscript -e '
        library(Nonpareil)
        s <- tryCatch(summary(Nonpareil.set("${meta.id}.npo", plot = FALSE)), error = function(e) NULL)
        cols <- c("kappa", "C", "LR", "modelR", "LRstar", "diversity")
        vals <- if (is.null(s)) rep(NA, length(cols)) else as.numeric(s[1, cols])
        out <- data.frame(sample = "${meta.id}", t(setNames(vals, cols)), check.names = FALSE)
        names(out) <- c("sample", "np_redundancy", "np_coverage", "np_seq_effort_bp", "np_model_r", "np_effort_for_95pct_bp", "np_diversity")
        write.table(out, "${meta.id}.nonpareil.tsv", sep = "\\t", quote = FALSE, row.names = FALSE)
    '
    """

    stub:
    """
    touch ${meta.id}.npo
    printf 'sample\\tnp_coverage\\n${meta.id}\\tNA\\n' > ${meta.id}.nonpareil.tsv
    """
}
