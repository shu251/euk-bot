include { FASTQC as FASTQC_RAW } from '../../modules/local/fastqc'
include { FASTP                } from '../../modules/local/fastp'
include { SORTMERNA            } from '../../modules/local/sortmerna'
include { MULTIQC              } from '../../modules/local/multiqc'

workflow QC_PREPROCESS {
    take:
    reads    // [meta, [R1, R2]]

    main:
    FASTQC_RAW(reads)
    FASTP(reads)

    ch_clean_reads   = FASTP.out.reads
    ch_sortmerna_log = channel.empty()
    if (!params.skip_rrna_removal) {
        ch_rrna_db = channel.fromPath(params.sortmerna_db, checkIfExists: true).collect()
        SORTMERNA(FASTP.out.reads, ch_rrna_db)
        ch_clean_reads   = SORTMERNA.out.reads
        ch_sortmerna_log = SORTMERNA.out.log
    }

    ch_multiqc_files = FASTQC_RAW.out.zip
        .mix(FASTP.out.json, ch_sortmerna_log)
        .map { _meta, f -> f }
        .collect()
    MULTIQC(ch_multiqc_files)

    emit:
    reads = ch_clean_reads    // [meta, [R1, R2]] trimmed, rRNA-depleted
}
