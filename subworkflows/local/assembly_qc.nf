include { SEQKIT_STATS        } from '../../modules/local/seqkit_stats'
include { ASSEMBLY_MAPRATE    } from '../../modules/local/assembly_maprate'
include { BUSCO               } from '../../modules/local/busco'
include { ASSEMBLY_QC_SUMMARY } from '../../modules/local/assembly_qc_summary'

/*
 * Optional (--run_assembly_qc): compares every assembly before clustering, so assemblers and
 * co-assembly can be judged on the same reads. BUSCO runs only when --busco_db_dir is set.
 */
workflow ASSEMBLY_QC {
    take:
    contigs        // [meta, contigs.fa], one per sample/group x assembler
    reads          // [meta, [R1, R2]] for every sample and co-assembly group, matched on meta.id

    main:
    SEQKIT_STATS(contigs.map { _meta, fa -> fa }.collect())

    ch_maprate_in = contigs
        .map { meta, fa -> [meta.id, meta, fa] }
        .combine(reads.map { meta, pair -> [meta.id, pair] }, by: 0)
        .map { _id, meta, fa, pair -> [meta, fa, pair] }
    ASSEMBLY_MAPRATE(ch_maprate_in)

    ch_busco = channel.value([])
    if (params.busco_db_dir) {
        BUSCO(contigs, file(params.busco_db_dir, checkIfExists: true))
        ch_busco = BUSCO.out.summaries.flatten().collect()
    }

    ASSEMBLY_QC_SUMMARY(SEQKIT_STATS.out.stats, ASSEMBLY_MAPRATE.out.maprate.collect(), ch_busco)

    emit:
    table = ASSEMBLY_QC_SUMMARY.out.table
}
