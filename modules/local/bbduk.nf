// Removes PhiX spike-in and Illumina sequencing artifacts by k-mer match against BBTools' bundled references
// (--bbduk_ref keywords); a pair is removed if either mate matches.
process BBDUK {
    tag "${meta.id}"
    label 'process_medium'
    conda 'bioconda::bbmap=40.02'
    container 'quay.io/biocontainers/bbmap:40.02--h09cc210_0'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_bbduk_{1,2}.fq.gz"),     emit: reads
    tuple val(meta), path("${meta.id}.bbduk.tsv"),              emit: stats
    tuple val(meta), path("${meta.id}.bbduk.log"),              emit: log
    tuple val(meta), path("${meta.id}.bbduk_contaminants.txt"), emit: contaminants

    script:
    def xmx = Math.max(1, (task.memory.toGiga() * 0.85) as int)
    """
    bbduk.sh -Xmx${xmx}g \\
        in=${reads[0]} in2=${reads[1]} \\
        out=${meta.id}_bbduk_1.fq.gz out2=${meta.id}_bbduk_2.fq.gz \\
        ref=${params.bbduk_ref} k=31 hdist=1 \\
        stats=${meta.id}.bbduk_contaminants.txt \\
        threads=${task.cpus} \\
        2> ${meta.id}.bbduk.log

    # Read pairs in vs. kept
    n_in=\$(( \$(zcat ${reads[0]} | wc -l) / 4 ))
    n_kept=\$(( \$(zcat ${meta.id}_bbduk_1.fq.gz | wc -l) / 4 ))
    printf 'sample\\tinput_pairs\\tkept_pairs\\n%s\\t%d\\t%d\\n' ${meta.id} \$n_in \$n_kept > ${meta.id}.bbduk.tsv
    """

    stub:
    """
    touch ${meta.id}_bbduk_1.fq.gz ${meta.id}_bbduk_2.fq.gz ${meta.id}.bbduk.log ${meta.id}.bbduk_contaminants.txt
    printf 'sample\\tinput_pairs\\tkept_pairs\\n' > ${meta.id}.bbduk.tsv
    """
}
