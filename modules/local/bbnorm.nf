// Normalizes read depth by k-mer counts, for the assemblers only: highly expressed transcripts are downsampled
// to --bbnorm_target, which cuts assembly time and memory. Salmon and reference mapping use all clean reads.
// mindepth (--bbnorm_mindepth) is lowered from BBNorm's default of 5 so reads of rare transcripts are kept.
process BBNORM {
    tag "${meta.id}"
    label 'process_medium'
    conda 'bioconda::bbmap=40.02'
    container 'quay.io/biocontainers/bbmap:40.02--h09cc210_0'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_norm_{1,2}.fq.gz"), emit: reads
    tuple val(meta), path("${meta.id}.bbnorm.tsv"),        emit: stats
    tuple val(meta), path("${meta.id}.bbnorm.log"),        emit: log

    script:
    def xmx = Math.max(1, (task.memory.toGiga() * 0.85) as int)
    """
    bbnorm.sh -Xmx${xmx}g \\
        in=${reads[0]} in2=${reads[1]} \\
        out=${meta.id}_norm_1.fq.gz out2=${meta.id}_norm_2.fq.gz \\
        target=${params.bbnorm_target} mindepth=${params.bbnorm_mindepth} \\
        prefilter=t tmpdir=. \\
        threads=${task.cpus} \\
        2> ${meta.id}.bbnorm.log

    # Read pairs in vs. kept for assembly
    n_in=\$(( \$(zcat ${reads[0]} | wc -l) / 4 ))
    n_kept=\$(( \$(zcat ${meta.id}_norm_1.fq.gz | wc -l) / 4 ))
    printf 'sample\\tinput_pairs\\tnormalized_pairs\\n%s\\t%d\\t%d\\n' ${meta.id} \$n_in \$n_kept > ${meta.id}.bbnorm.tsv
    """

    stub:
    """
    touch ${meta.id}_norm_1.fq.gz ${meta.id}_norm_2.fq.gz ${meta.id}.bbnorm.log
    printf 'sample\\tinput_pairs\\tnormalized_pairs\\n' > ${meta.id}.bbnorm.tsv
    """
}
