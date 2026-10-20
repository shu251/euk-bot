include { FASTQC as FASTQC_RAW } from '../../modules/local/fastqc'
include { FASTP                } from '../../modules/local/fastp'
include { RIBODETECTOR         } from '../../modules/local/ribodetector'
include { MULTIQC              } from '../../modules/local/multiqc'

workflow QC_PREPROCESS {
    take:
    reads    // [meta, [R1, R2]]

    main:
    FASTQC_RAW(reads)
    FASTP(reads)

    ch_clean_reads   = FASTP.out.reads
    if (!params.skip_rrna_removal) {
        RIBODETECTOR(FASTP.out.reads)
        ch_clean_reads = RIBODETECTOR.out.reads
        // One table of per-sample rRNA fractions
        RIBODETECTOR.out.stats
            .map { _meta, tsv -> tsv }
            .collectFile(name: 'rrna_summary.tsv', keepHeader: true, skip: 1, sort: true,
                         storeDir: "${params.outdir}/qc/ribodetector")
    }

    ch_multiqc_files = FASTQC_RAW.out.zip
        .mix(FASTP.out.json)
        .map { _meta, f -> f }
        .collect()
    MULTIQC(ch_multiqc_files)

    emit:
    reads = ch_clean_reads    // [meta, [R1, R2]] trimmed, rRNA-depleted
}
