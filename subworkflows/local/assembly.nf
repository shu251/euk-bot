include { RNASPADES      } from '../../modules/local/rnaspades'
include { MEGAHIT        } from '../../modules/local/megahit'
include { CDHIT_EST      } from '../../modules/local/cdhit_est'
include { TRANSDECODER   } from '../../modules/local/transdecoder'
include { DIAMOND_BLASTP } from '../../modules/local/diamond_blastp'
include { EGGNOG_MAPPER  } from '../../modules/local/eggnog_mapper'
include { SALMON_INDEX   } from '../../modules/local/salmon_index'
include { SALMON_QUANT   } from '../../modules/local/salmon_quant'
include { MERGE_SALMON   } from '../../modules/local/merge_salmon'

/*
 * Per-sample assemblies are pooled and clustered into one non-redundant
 * reference, so every sample is quantified against the same contig set.
 */
workflow ASSEMBLY {
    take:
    reads    // [meta, [R1, R2]]

    main:
    if (params.assembler == 'rnaspades') {
        RNASPADES(reads)
        ch_contigs = RNASPADES.out.contigs
    } else if (params.assembler == 'megahit') {
        MEGAHIT(reads)
        ch_contigs = MEGAHIT.out.contigs
    } else {
        error "Unknown --assembler '${params.assembler}' (use 'rnaspades' or 'megahit')"
    }

    CDHIT_EST(ch_contigs.map { _meta, fa -> fa }.collect())
    TRANSDECODER(CDHIT_EST.out.fasta)

    if (params.reference_dmnd) {
        DIAMOND_BLASTP(TRANSDECODER.out.pep, file(params.reference_dmnd, checkIfExists: true))
    }
    if (params.eggnog_data_dir) {
        EGGNOG_MAPPER(TRANSDECODER.out.pep, file(params.eggnog_data_dir, checkIfExists: true))
    }

    SALMON_INDEX(CDHIT_EST.out.fasta)
    SALMON_QUANT(reads, SALMON_INDEX.out.index)
    MERGE_SALMON(SALMON_QUANT.out.quant.map { _meta, dir -> dir }.collect())

    emit:
    contigs = CDHIT_EST.out.fasta
    pep     = TRANSDECODER.out.pep
    counts  = MERGE_SALMON.out.numreads
    tpm     = MERGE_SALMON.out.tpm
}
