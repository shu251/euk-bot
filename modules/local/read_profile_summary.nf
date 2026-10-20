// One row per sample with read QC, depth/diversity, composition, strandedness and read-based
// similarity, plus suggested co-assembly groups and empty columns to fill in by hand.
process READ_PROFILE_SUMMARY {
    label 'process_low'
    conda 'conda-forge::python=3.9'
    container 'quay.io/biocontainers/python:3.9--1'   // has ps, which Nextflow needs; slim images don't

    input:
    path(samplesheet)
    path(fastp_json)
    path(rrna_stats)
    path(similarity)
    path(nonpareil)
    path(strandedness)    // optional: pass [] to skip
    path(kaiju)           // optional: pass [] to skip

    output:
    path('read_profile.tsv'),           emit: table
    path('samplesheet_suggested.csv'),  emit: samplesheet

    script:
    def opt = [
        rrna_stats   ? "--rrna ${rrna_stats}" : '',
        strandedness ? "--strandedness ${strandedness}" : '',
        kaiju        ? "--kaiju ${kaiju}" : '',
    ].join(' ')
    """
    summarize_read_profile.py \\
        --samplesheet ${samplesheet} \\
        --fastp ${fastp_json} \\
        --similarity ${similarity} \\
        --nonpareil ${nonpareil} \\
        --group-similarity ${params.profile_group_similarity} \\
        ${opt}
    """

    stub:
    """
    touch read_profile.tsv samplesheet_suggested.csv
    """
}
