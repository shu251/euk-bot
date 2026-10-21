include { CAT_FASTQ as CAT_GROUP_READS } from '../../modules/local/cat_fastq'
include { RNASPADES       } from '../../modules/local/rnaspades'
include { MEGAHIT         } from '../../modules/local/megahit'
include { MMSEQS_LINCLUST } from '../../modules/local/mmseqs_linclust'
include { ASSEMBLY_QC     } from './assembly_qc'
include { TRANSDECODER   } from '../../modules/local/transdecoder'
include { METAEUK        } from '../../modules/local/metaeuk'
include { DIAMOND_BLASTP } from '../../modules/local/diamond_blastp'
include { EGGNOG_MAPPER  } from '../../modules/local/eggnog_mapper'
include { SALMON_INDEX   } from '../../modules/local/salmon_index'
include { SALMON_QUANT   } from '../../modules/local/salmon_quant'
include { MERGE_SALMON   } from '../../modules/local/merge_salmon'

/*
 * Each sample is assembled with every assembler in --assemblers, and samples sharing a samplesheet
 * group (e.g. a vent site) are also co-assembled with --coassembly_assemblers. All contigs are pooled
 * and clustered into one non-redundant reference, so every sample is quantified against the same set.
 */
workflow ASSEMBLY {
    take:
    reads            // [meta, [R1, R2]] clean reads, for quantification; meta.group is set when the samplesheet has a group column
    assembly_reads   // [meta, [R1, R2]] assembler input: the clean reads, normalized by BBNorm unless skipped

    main:
    def assemblers    = params.assemblers.tokenize(',')*.trim()
    def co_assemblers = params.coassembly_assemblers ? params.coassembly_assemblers.tokenize(',')*.trim() : []
    def unknown       = (assemblers + co_assemblers) - ['rnaspades', 'megahit']
    if (unknown) {
        error "Unknown assembler(s): ${unknown.join(', ')} (use rnaspades, megahit)"
    }

    ch_rnaspades_in = 'rnaspades' in assemblers ? assembly_reads : channel.empty()
    ch_megahit_in   = 'megahit' in assemblers ? assembly_reads : channel.empty()
    ch_assembled_reads = assembly_reads    // reads behind each assembly, for assembly QC

    // Co-assembly: pool the reads of each group with more than one sample
    if (co_assemblers) {
        ch_group_reads = assembly_reads
            .filter { meta, _reads -> meta.group }
            .map { meta, pair -> [meta.group, pair[0], pair[1]] }
            .groupTuple()
            .filter { _group, r1s, _r2s -> r1s.size() > 1 }
            .map { group, r1s, r2s -> [[id: "coasm_${group}", group: group], r1s, r2s] }
        CAT_GROUP_READS(ch_group_reads)
        ch_assembled_reads = ch_assembled_reads.mix(CAT_GROUP_READS.out.reads)
        if ('rnaspades' in co_assemblers) {
            ch_rnaspades_in = ch_rnaspades_in.mix(CAT_GROUP_READS.out.reads)
        }
        if ('megahit' in co_assemblers) {
            ch_megahit_in = ch_megahit_in.mix(CAT_GROUP_READS.out.reads)
        }
    }

    RNASPADES(ch_rnaspades_in)
    MEGAHIT(ch_megahit_in)
    ch_contigs = RNASPADES.out.contigs.mix(MEGAHIT.out.contigs)

    if (params.run_assembly_qc) {
        ASSEMBLY_QC(ch_contigs, ch_assembled_reads)
    }

    MMSEQS_LINCLUST(ch_contigs.map { _meta, fa -> fa }.collect())
    // ORFs: TransDecoder (ab initio), plus MetaEuk (reference homology) when --metaeuk_db is set; both
    // protein sets are kept and annotated, named by ORF caller
    TRANSDECODER(MMSEQS_LINCLUST.out.fasta)
    ch_pep = TRANSDECODER.out.pep.map { pep -> ['transdecoder', pep] }
    if (params.metaeuk_db) {
        METAEUK(MMSEQS_LINCLUST.out.fasta, file(params.metaeuk_db, checkIfExists: true))
        ch_pep = ch_pep.mix(METAEUK.out.pep.map { pep -> ['metaeuk', pep] })
    }

    if (params.reference_dmnd) {
        DIAMOND_BLASTP(ch_pep, file(params.reference_dmnd, checkIfExists: true))
    }
    if (params.eggnog_data_dir) {
        EGGNOG_MAPPER(ch_pep, file(params.eggnog_data_dir, checkIfExists: true))
    }

    SALMON_INDEX(MMSEQS_LINCLUST.out.fasta)
    SALMON_QUANT(reads, SALMON_INDEX.out.index)
    MERGE_SALMON(SALMON_QUANT.out.quant.map { _meta, dir -> dir }.collect())

    emit:
    contigs = MMSEQS_LINCLUST.out.fasta
    pep     = ch_pep    // [orf_caller, proteins.fa]
    counts  = MERGE_SALMON.out.numreads
    tpm     = MERGE_SALMON.out.tpm
}
