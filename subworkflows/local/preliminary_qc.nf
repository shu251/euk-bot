include { FASTQC as FASTQC_RAW     } from '../../modules/local/fastqc'
include { FASTQC as FASTQC_CLEAN   } from '../../modules/local/fastqc'
include { FASTQC as FASTQC_NORM    } from '../../modules/local/fastqc'
include { FASTP                    } from '../../modules/local/fastp'
include { BBDUK                    } from '../../modules/local/bbduk'
include { RIBODETECTOR             } from '../../modules/local/ribodetector'
include { BBNORM                   } from '../../modules/local/bbnorm'
include { MULTIQC as MULTIQC_RAW   } from '../../modules/local/multiqc'
include { MULTIQC as MULTIQC_CLEAN } from '../../modules/local/multiqc'
include { QC_SUMMARY               } from '../../modules/local/qc_summary'

/*
 * Step 1, preliminary QC:
 *   raw reads -> FastQC -> fastp (adapter + quality trimming, min length) -> BBDuk (PhiX, artifacts)
 *     -> RiboDetector (rRNA removal) = clean reads, used by every later step -> FastQC
 *   clean reads -> BBNorm (depth normalization, assembler input only) -> FastQC
 * MultiQC reports for the raw reads and for the clean + normalized reads, and a per-sample read accounting
 * table (qc/preliminary_qc.tsv). Later steps start only once both reports are done and the read-loss check
 * (checkReadLoss) has passed.
 */
workflow PRELIMINARY_QC {
    take:
    reads    // [meta, [R1, R2]]

    main:
    FASTQC_RAW(reads)
    FASTP(reads)
    ch_clean_reads = FASTP.out.reads
    ch_clean_logs  = channel.empty()    // tool logs for the clean MultiQC report

    ch_bbduk_stats = channel.empty()
    if (!params.skip_bbduk) {
        BBDUK(ch_clean_reads)
        ch_clean_reads = BBDUK.out.reads
        ch_bbduk_stats = BBDUK.out.stats
        ch_clean_logs  = ch_clean_logs.mix(BBDUK.out.log)
    }

    ch_rrna_reads = channel.empty()
    ch_rrna_stats = channel.empty()
    if (!params.skip_rrna_removal) {
        RIBODETECTOR(ch_clean_reads)
        ch_clean_reads = RIBODETECTOR.out.reads
        ch_rrna_reads  = RIBODETECTOR.out.rrna
        ch_rrna_stats  = RIBODETECTOR.out.stats
        // One table of per-sample rRNA fractions
        RIBODETECTOR.out.stats
            .map { _meta, tsv -> tsv }
            .collectFile(name: 'rrna_summary.tsv', keepHeader: true, skip: 1, sort: true,
                         storeDir: "${params.outdir}/qc/ribodetector")
    }

    FASTQC_CLEAN(ch_clean_reads)
    ch_clean_fastqc = FASTQC_CLEAN.out.zip

    // Normalized reads only feed the assemblers, so skip BBNorm when nothing gets assembled
    ch_assembly_reads = ch_clean_reads
    ch_norm_stats     = channel.empty()
    if (!params.skip_normalization && params.run_assembly && !params.stop_after_profiling) {
        BBNORM(ch_clean_reads)
        ch_assembly_reads = BBNORM.out.reads
        ch_norm_stats     = BBNORM.out.stats
        FASTQC_NORM(BBNORM.out.reads)
        ch_clean_fastqc   = ch_clean_fastqc.mix(FASTQC_NORM.out.zip)
    }

    MULTIQC_RAW('raw', FASTQC_RAW.out.zip.mix(FASTP.out.json).map { _meta, f -> f }.collect())
    MULTIQC_CLEAN('clean', ch_clean_fastqc.mix(ch_clean_logs).map { _meta, f -> f }.collect())

    QC_SUMMARY(
        FASTP.out.json.map { _meta, f -> f }.collect(),
        ch_bbduk_stats.map { _meta, f -> f }.collect().ifEmpty([]),
        ch_rrna_stats.map { _meta, f -> f }.collect().ifEmpty([]),
        ch_norm_stats.map { _meta, f -> f }.collect().ifEmpty([]),
    )

    // Emits once, after both MultiQC reports and the QC table exist and the read-loss check passed
    ch_qc_done = QC_SUMMARY.out.table
        .combine(MULTIQC_RAW.out.report)
        .combine(MULTIQC_CLEAN.out.report)
        .map { table, _raw, _clean -> checkReadLoss(table) }

    emit:
    reads          = ch_clean_reads.combine(ch_qc_done).map { meta, r, _ok -> [meta, r] }     // trimmed, filtered, rRNA-depleted
    assembly_reads = ch_assembly_reads.combine(ch_qc_done).map { meta, r, _ok -> [meta, r] }  // clean reads, normalized if BBNorm ran
    rrna           = ch_rrna_reads     // [meta, [R1, R2]] rRNA reads, for an optional 18S branch
    fastp_json     = FASTP.out.json    // [meta, json]
    rrna_stats     = ch_rrna_stats     // [meta, tsv] per-sample rRNA fraction
}

// Samples flagged in preliminary_qc.tsv lost more than --max_pct_removed_fastp of their read pairs in fastp,
// usually to the length filter. At a minimum length above 36 bp the run stops so it can be lowered first;
// at 36 bp or below it only warns.
def checkReadLoss(table) {
    if (table.size() == 0) {    // stub runs
        return true
    }
    def flagged = table.splitCsv(header: true, sep: '\t').findAll { row -> row.qc_flag }
    if (flagged) {
        def samples = flagged.collect { row -> "${row.sample} (${row.pct_removed_fastp}%)" }.join(', ')
        def msg = "Preliminary QC: ${flagged.size()} sample(s) lost more than ${params.max_pct_removed_fastp}% of read pairs " +
            "in fastp at --min_read_length ${params.min_read_length}: ${samples}. See qc/preliminary_qc.tsv and the MultiQC reports."
        if (params.min_read_length > 36) {
            error "${msg}\nRerun with --min_read_length 36 -resume, or with --max_pct_removed_fastp 100 to continue as is."
        }
        log.warn msg
    }
    return true
}
