// One row per sample: raw read pairs split into removed by fastp / removed by BBDuk / rRNA / non-rRNA
// (counts and % of raw), pairs kept by BBNorm for assembly, fastp quality metrics and a read-loss flag.
process QC_SUMMARY {
    label 'process_low'
    conda 'conda-forge::python=3.9'
    container 'quay.io/biocontainers/python:3.9--1'   // has ps, which Nextflow needs; slim images don't

    input:
    path(fastp_json)
    path(bbduk_stats)     // optional: pass [] to skip
    path(rrna_stats)      // optional: pass [] to skip
    path(bbnorm_stats)    // optional: pass [] to skip

    output:
    path('preliminary_qc.tsv'), emit: table

    script:
    def opt = [
        bbduk_stats  ? "--bbduk ${bbduk_stats}" : '',
        rrna_stats   ? "--rrna ${rrna_stats}" : '',
        bbnorm_stats ? "--bbnorm ${bbnorm_stats}" : '',
    ].join(' ')
    """
    summarize_qc.py --fastp ${fastp_json} ${opt} \\
        --max-pct-removed-fastp ${params.max_pct_removed_fastp} \\
        --out preliminary_qc.tsv
    """

    stub:
    """
    touch preliminary_qc.tsv
    """
}
