include { SOURMASH_SKETCH      } from '../../modules/local/sourmash_sketch'
include { SOURMASH_COMPARE     } from '../../modules/local/sourmash_compare'
include { NONPAREIL            } from '../../modules/local/nonpareil'
include { KAIJU                } from '../../modules/local/kaiju'
include { READ_PROFILE_SUMMARY } from '../../modules/local/read_profile_summary'

/*
 * Step 2, read screening (optional, --run_read_profiling): profiles the QC'd reads before assembly, to judge
 * read quality and choose co-assembly groups. The strandedness calls come from main.nf (STRANDEDNESS runs there,
 * since the assembly also uses them); Kaiju runs when --kaiju_db is set.
 */
workflow READ_SCREENING {
    take:
    reads         // [meta, [R1, R2]] trimmed, rRNA-depleted
    fastp_json    // [meta, json]
    rrna_stats    // [meta, tsv]; empty with --skip_rrna_removal
    strandedness  // [meta, tsv]; empty without a protein DB for the check

    main:
    SOURMASH_SKETCH(reads)
    SOURMASH_COMPARE(SOURMASH_SKETCH.out.sig.collect())
    NONPAREIL(reads)

    ch_kaiju = channel.value([])
    if (params.kaiju_db) {
        KAIJU(reads, file(params.kaiju_db, checkIfExists: true))
        ch_kaiju = KAIJU.out.summary.collect()
    }

    READ_PROFILE_SUMMARY(
        file(params.input, checkIfExists: true),
        fastp_json.map { _meta, f -> f }.collect(),
        rrna_stats.map { _meta, f -> f }.collect().ifEmpty([]),
        SOURMASH_COMPARE.out.csv,
        NONPAREIL.out.summary.collect(),
        strandedness.map { _meta, f -> f }.collect().ifEmpty([]),
        ch_kaiju,
    )

    emit:
    table       = READ_PROFILE_SUMMARY.out.table
    samplesheet = READ_PROFILE_SUMMARY.out.samplesheet
}
