include { DIAMOND_BLASTX } from '../../modules/local/diamond_blastx'
include { SUMMARIZE_HITS } from '../../modules/local/summarize_hits'

/*
 * Assembly-free branch: translate-search reads against a reference protein DB
 * and count best hits per reference protein.
 */
workflow REFERENCE_MAP {
    take:
    reads    // [meta, [R1, R2]]

    main:
    DIAMOND_BLASTX(reads, file(params.reference_dmnd, checkIfExists: true))

    taxonomy = params.reference_taxonomy ? file(params.reference_taxonomy, checkIfExists: true) : []
    SUMMARIZE_HITS(DIAMOND_BLASTX.out.hits.map { _meta, tsv -> tsv }.collect(), taxonomy)

    emit:
    hits   = DIAMOND_BLASTX.out.hits
    counts = SUMMARIZE_HITS.out.counts
}
