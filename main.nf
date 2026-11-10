#!/usr/bin/env nextflow
/*
 * euk-bot: metatranscriptome processing for paired-end FASTQ
 *
 *   PRELIMINARY_QC : FastQC -> fastp (trimming) -> BBDuk (PhiX) -> RiboDetector (rRNA) -> FastQC; BBNorm for the
 *                    assembler input; MultiQC before/after; per-sample QC table and read-loss check
 *   STRANDEDNESS   : library strandedness per sample (DIAMOND blastx of R1 reads), for the read profile and rnaSPAdes
 *   READ_SCREENING : optional pre-assembly screen (sourmash, Nonpareil, strandedness, Kaiju) and suggested
 *                    co-assembly groups
 *   ASSEMBLY       : per-sample (rnaSPAdes + MEGAHIT) and per-group co-assembly -> MMseqs2 clustering -> TransDecoder
 *                    (+ MetaEuk) -> DIAMOND/eggNOG annotation -> salmon quantification
 *   REFERENCE_MAP  : DIAMOND blastx of reads against a reference protein DB
 * Each run writes <outdir>/methods.txt: software, versions and key arguments of every step it uses.
 */

include { CAT_FASTQ      } from './modules/local/cat_fastq'
include { PRELIMINARY_QC } from './subworkflows/local/preliminary_qc'
include { STRANDEDNESS   } from './modules/local/strandedness'
include { READ_SCREENING } from './subworkflows/local/read_screening'
include { ASSEMBLY       } from './subworkflows/local/assembly'
include { REFERENCE_MAP  } from './subworkflows/local/reference_map'
include { writeMethods   } from './subworkflows/local/methods'

workflow {
    if (!params.input) {
        error "Provide a samplesheet with --input (columns: sample,fastq_1,fastq_2[,group])"
    }
    if (params.run_reference && !params.reference_dmnd && !params.stop_after_profiling) {
        error "The reference branch needs a DIAMOND database: --reference_dmnd (or --run_reference false)"
    }
    if (!(params.strandedness in ['auto', 'rf', 'fr', 'unstranded'])) {
        error "--strandedness must be auto, rf, fr or unstranded (got '${params.strandedness}')"
    }
    def strand_db = params.strand_check_dmnd ?: params.reference_dmnd
    def profiling = params.run_read_profiling || params.stop_after_profiling
    def assembling = params.run_assembly && !params.stop_after_profiling
    if (assembling && params.strandedness == 'auto' && !strand_db) {
        log.warn "--strandedness auto: no protein DB for the strandedness check (--strand_check_dmnd or --reference_dmnd), " +
            "so rnaSPAdes assembles every sample as unstranded"
    }
    if (assembling && !params.reference_dmnd && (params.eukulele_db || params.reference_annotations)) {
        log.warn "--eukulele_db and --reference_annotations use the proteins' DIAMOND hits, so they need --reference_dmnd; skipped"
    }
    writeMethods()

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

    PRELIMINARY_QC(ch_reads)

    // Strandedness per sample: reported in the read profile, and used by rnaSPAdes with --strandedness auto
    ch_strand = channel.empty()
    if (strand_db && (profiling || (assembling && params.strandedness == 'auto'))) {
        STRANDEDNESS(PRELIMINARY_QC.out.reads, file(strand_db, checkIfExists: true))
        ch_strand = STRANDEDNESS.out.summary
    }

    if (profiling) {
        READ_SCREENING(PRELIMINARY_QC.out.reads, PRELIMINARY_QC.out.fastp_json, PRELIMINARY_QC.out.rrna_stats, ch_strand)
    }
    if (params.stop_after_profiling) {
        return
    }

    if (params.run_assembly) {
        // [sample id, rf | fr | unstranded | too_few_hits] from each sample's strandedness table
        ch_strand_calls = ch_strand.map { meta, tsv -> [meta.id, tsv.splitCsv(header: true, sep: '\t')[0].strandedness] }
        ASSEMBLY(PRELIMINARY_QC.out.reads, PRELIMINARY_QC.out.assembly_reads, ch_strand_calls)
    }
    if (params.run_reference) {
        REFERENCE_MAP(PRELIMINARY_QC.out.reads)
    }
}
