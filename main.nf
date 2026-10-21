#!/usr/bin/env nextflow
/*
 * euk-bot: metatranscriptome processing for paired-end FASTQ
 *
 *   QC_PREPROCESS : FastQC -> fastp -> RiboDetector (rRNA removal) -> MultiQC
 *   READ_PROFILE  : optional pre-assembly profile (sourmash, Nonpareil, strandedness, Kaiju) for grouping
 *   ASSEMBLY      : per-sample (rnaSPAdes + MEGAHIT) and per-group co-assembly -> MMseqs2 clustering -> TransDecoder
 *                   -> DIAMOND/eggNOG annotation -> salmon quantification
 *   REFERENCE_MAP : DIAMOND blastx of reads against a reference protein DB
 */

include { CAT_FASTQ     } from './modules/local/cat_fastq'
include { QC_PREPROCESS } from './subworkflows/local/qc_preprocess'
include { READ_PROFILE  } from './subworkflows/local/read_profile'
include { ASSEMBLY      } from './subworkflows/local/assembly'
include { REFERENCE_MAP } from './subworkflows/local/reference_map'

workflow {
    if (!params.input) {
        error "Provide a samplesheet with --input (columns: sample,fastq_1,fastq_2[,group])"
    }
    if (params.run_reference && !params.reference_dmnd && !params.stop_after_profiling) {
        error "The reference branch needs a DIAMOND database: --reference_dmnd (or --run_reference false)"
    }

    // Rows sharing a sample name (e.g. one per sequencing lane) are concatenated before QC.
    // The optional group column (e.g. vent site) defines co-assemblies.
    ch_samples = channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            [row.sample, row.group ?: '', [file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true)]]
        }
        .groupTuple()
        .map { sample, groups, pairs ->
            if (groups.toUnique().size() > 1) {
                error "Sample '${sample}' has rows with different groups: ${groups.toUnique().join(', ')}"
            }
            def meta = groups[0] ? [id: sample, group: groups[0]] : [id: sample]
            def sorted_pairs = pairs.sort { pair -> pair[0].name }
            [meta, sorted_pairs.collect { pair -> pair[0] }, sorted_pairs.collect { pair -> pair[1] }]
        }
        .branch { _meta, r1s, _r2s ->
            single: r1s.size() == 1
            multi:  true
        }

    CAT_FASTQ(ch_samples.multi)

    ch_reads = ch_samples.single
        .map { meta, r1s, r2s -> [meta, [r1s[0], r2s[0]]] }
        .mix(CAT_FASTQ.out.reads)

    QC_PREPROCESS(ch_reads)

    if (params.run_read_profiling || params.stop_after_profiling) {
        READ_PROFILE(QC_PREPROCESS.out.reads, QC_PREPROCESS.out.fastp_json, QC_PREPROCESS.out.rrna_stats)
    }
    if (params.stop_after_profiling) {
        return
    }

    if (params.run_assembly) {
        ASSEMBLY(QC_PREPROCESS.out.reads)
    }
    if (params.run_reference) {
        REFERENCE_MAP(QC_PREPROCESS.out.reads)
    }
}
