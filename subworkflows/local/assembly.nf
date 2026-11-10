include { CAT_FASTQ as CAT_GROUP_READS } from '../../modules/local/cat_fastq'
include { RNASPADES       } from '../../modules/local/rnaspades'
include { MEGAHIT         } from '../../modules/local/megahit'
include { MMSEQS_LINCLUST } from '../../modules/local/mmseqs_linclust'
include { ASSEMBLY_QC     } from './assembly_qc'
include { TRANSDECODER   } from '../../modules/local/transdecoder'
include { METAEUK        } from '../../modules/local/metaeuk'
include { DIAMOND_BLASTP } from '../../modules/local/diamond_blastp'
include { EGGNOG_MAPPER  } from '../../modules/local/eggnog_mapper'
include { EUKULELE       } from '../../modules/local/eukulele'
include { SALMON_INDEX   } from '../../modules/local/salmon_index'
include { SALMON_QUANT   } from '../../modules/local/salmon_quant'
include { MERGE_SALMON   } from '../../modules/local/merge_salmon'
include { CONTIG_TABLE   } from '../../modules/local/contig_table'

/*
 * Each sample is assembled with every assembler in --assemblers, and samples sharing a samplesheet
 * group (e.g. a vent site) are also co-assembled with --coassembly_assemblers. All contigs are pooled
 * and clustered into one non-redundant reference, so every sample is quantified against the same set.
 * rnaSPAdes is told the library type (meta.strandedness -> --ss); MEGAHIT has no stranded mode.
 */
workflow ASSEMBLY {
    take:
    reads            // [meta, [R1, R2]] clean reads, for quantification; meta.group is set when the samplesheet has a group column
    assembly_reads   // [meta, [R1, R2]] assembler input: the clean reads, normalized by BBNorm unless skipped
    strand_calls     // [sample id, call] from STRANDEDNESS; empty when it did not run

    main:
    def assemblers    = params.assemblers.tokenize(',')*.trim()
    def co_assemblers = params.coassembly_assemblers ? params.coassembly_assemblers.tokenize(',')*.trim() : []
    def unknown       = (assemblers + co_assemblers) - ['rnaspades', 'megahit']
    if (unknown) {
        error "Unknown assembler(s): ${unknown.join(', ')} (use rnaspades, megahit)"
    }

    // Library type per sample: --strandedness rf|fr|unstranded for all, or with auto the sample's own call
    // (unstranded when there is no rf/fr call, e.g. too_few_hits or no strandedness check)
    ch_assembly_reads = assembly_reads.map { meta, pair -> [meta + [strandedness: params.strandedness], pair] }
    if (params.strandedness == 'auto') {
        ch_assembly_reads = assembly_reads
            .map { meta, pair -> [meta.id, meta, pair] }
            .join(strand_calls, remainder: true)
            .filter { row -> row[1] instanceof Map }    // drop calls without assembler input (none expected)
            .map { _id, meta, pair, call -> [meta + [strandedness: call in ['rf', 'fr'] ? call : 'unstranded'], pair] }
    }

    ch_rnaspades_in = 'rnaspades' in assemblers ? ch_assembly_reads : channel.empty()
    ch_megahit_in   = 'megahit' in assemblers ? ch_assembly_reads : channel.empty()
    ch_assembled_reads = ch_assembly_reads    // reads behind each assembly, for assembly QC

    // Co-assembly: pool the reads of each group with more than one sample; stranded only if all its samples agree
    if (co_assemblers) {
        ch_group_reads = ch_assembly_reads
            .filter { meta, _reads -> meta.group }
            .map { meta, pair -> [meta.group, meta.strandedness, pair[0], pair[1]] }
            .groupTuple()
            .filter { _group, _strands, r1s, _r2s -> r1s.size() > 1 }
            .map { group, strands, r1s, r2s ->
                def strandedness = strands.toUnique().size() == 1 ? strands[0] : 'unstranded'
                [[id: "coasm_${group}", group: group, strandedness: strandedness], r1s, r2s]
            }
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

    // Annotation per ORF caller, all keyed by caller name: [orf_caller, file]
    ch_hits     = channel.empty()
    ch_eggnog   = channel.empty()
    ch_eukulele = channel.empty()
    if (params.reference_dmnd) {
        DIAMOND_BLASTP(ch_pep, file(params.reference_dmnd, checkIfExists: true))
        ch_hits = DIAMOND_BLASTP.out.hits
        // Consensus taxonomy from the top hits (DIAMOND_BLASTP keeps 100 per protein when --eukulele_db is set)
        if (params.eukulele_db) {
            EUKULELE(ch_pep.join(ch_hits), file(params.eukulele_db, checkIfExists: true))
            ch_eukulele = EUKULELE.out.taxonomy
        }
    }
    if (params.eggnog_data_dir) {
        EGGNOG_MAPPER(ch_pep, file(params.eggnog_data_dir, checkIfExists: true))
        ch_eggnog = EGGNOG_MAPPER.out.annotations
    }

    SALMON_INDEX(MMSEQS_LINCLUST.out.fasta)
    SALMON_QUANT(reads, SALMON_INDEX.out.index)
    MERGE_SALMON(SALMON_QUANT.out.quant.map { _meta, dir -> dir }.collect())

    // One table per contig: representative protein, taxonomy, function, counts
    ch_annotation_files = ch_pep.mix(ch_hits, ch_eggnog, ch_eukulele).map { _name, f -> f }.collect()
    CONTIG_TABLE(
        MERGE_SALMON.out.numreads,
        ch_annotation_files,
        params.reference_dmnd && params.reference_taxonomy ? file(params.reference_taxonomy, checkIfExists: true) : [],
        params.reference_dmnd && params.reference_annotations ? file(params.reference_annotations, checkIfExists: true) : []
    )

    emit:
    contigs = MMSEQS_LINCLUST.out.fasta
    pep     = ch_pep    // [orf_caller, proteins.fa]
    counts  = MERGE_SALMON.out.numreads
    tpm     = MERGE_SALMON.out.tpm
    table   = CONTIG_TABLE.out.contigs
}
