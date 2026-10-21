/*
 * <outdir>/methods.txt: each step this run uses, with software, version (read from each module's conda and
 * container lines, so they always match what runs) and key arguments. Written at launch, so it is there
 * while the run is going; a -resume run rewrites it.
 * The argument text mirrors the module scripts: update it here when a module's arguments change.
 */

// process name -> [conda: 'bioconda::fastp=0.23.4', container: 'quay.io/biocontainers/fastp:...']
def moduleSoftware() {
    def software = [:]
    files("${projectDir}/modules/local/*.nf").each { f ->
        def text = f.text
        def name = text.find(/process\s+(\w+)/) { match -> match[1] }
        if (name) {
            software[name] = [
                conda:     text.find(/conda\s+'([^']+)'/) { match -> match[1] } ?: '',
                container: text.find(/container\s+'([^']+)'/) { match -> match[1] } ?: '',
            ]
        }
    }
    return software
}

def gitCommit() {
    try {
        def dir = projectDir.toString()
        def commit = ['git', '-C', dir, 'rev-parse', '--short', 'HEAD'].execute().text.trim()
        def dirty = ['git', '-C', dir, 'status', '--porcelain', '--untracked-files=no'].execute().text.trim()
        return commit ? commit + (dirty ? ' + uncommitted changes' : '') : 'unknown (not a git checkout)'
    } catch (Exception _e) {
        return 'unknown (git not available)'
    }
}

// First lines of a database's MANIFEST.txt (name, build date, sources), if it has one
def manifestLines(dir, stopAt) {
    def manifest = file(dir).resolve('MANIFEST.txt')
    if (!manifest.exists()) {
        return []
    }
    return manifest.readLines().takeWhile { line -> !line.startsWith(stopAt) }.findAll { line -> line.trim() }
}

// methods.txt is built up in ctx: [out: lines, step: section number, software: moduleSoftware(), engine: ...]
def addSection(ctx, title) {
    ctx.step += 1
    def heading = "${ctx.step}. ${title}"
    ctx.out << '' << heading << ('-' * heading.size())
}

// One tool: "<label> <version>   [image or conda package]", then its key arguments / notes
def addTool(ctx, process, label, notes) {
    def s = ctx.software[process] ?: [conda: '', container: '']
    def version = s.conda ? s.conda.tokenize(' ')[0].tokenize('=')[-1] : ''
    def source = ctx.engine ? s.container : s.conda
    ctx.out << '' << "${label}${version ? ' ' + version : ''}${source ? '   [' + source + ']' : ''}"
    notes.each { line -> ctx.out << "    ${line}" }
}

// A euk-bot helper script (bin/), run in the pipeline's python image
def addScript(ctx, script, notes) {
    def image = ctx.engine ? '   [' + ctx.software['QC_SUMMARY'].container + ']' : ''
    ctx.out << '' << "euk-bot ${script}${image}"
    notes.each { line -> ctx.out << "    ${line}" }
}

def writeMethods() {
    def engine = workflow.containerEngine
    def ctx = [out: [], step: 0, software: moduleSoftware(), engine: engine]
    def out = ctx.out

    out << 'euk-bot methods'
    out << '==============='
    out << "Written at launch: the steps, software and settings this run uses."
    out << ''
    out << "date used:       ${workflow.start.format('yyyy-MM-dd')} (launched ${workflow.start.format('HH:mm')})"
    out << "project folder:  ${workflow.launchDir}"
    out << "results:         ${file(params.outdir)}"
    out << "samplesheet:     ${params.input}"
    out << "run:             ${workflow.runName}"
    out << "command:         ${workflow.commandLine}"
    out << "euk-bot:         ${gitCommit()}"
    out << "Nextflow:        ${nextflow.version}"
    out << "profile:         ${workflow.profile}"
    out << "software:        ${engine ? engine + ' images' : 'conda environments'} (shown in [ ] for each tool)"

    // 1. Preliminary QC
    def normalize = !params.skip_normalization && params.run_assembly && !params.stop_after_profiling
    addSection(ctx, 'Preliminary QC')
    out << ''
    out << 'Order: raw reads -> FastQC -> fastp -> ' + (params.skip_bbduk ? '' : 'BBDuk -> ') +
        (params.skip_rrna_removal ? '' : 'RiboDetector -> ') + 'clean reads (all later steps) -> FastQC' +
        (normalize ? '; clean reads -> BBNorm (assembler input only) -> FastQC.' : '.')
    out << ''
    out << 'cat'
    out << '    Samples sequenced on several lanes: lane FASTQs concatenated per sample before QC.'
    addTool(ctx, 'FASTQC', 'FastQC', ['Default settings; raw reads, clean reads' + (normalize ? ' and normalized reads.' : '.')])
    addTool(ctx, 'FASTP', 'fastp', [
        'Adapters auto-detected for paired-end reads (--detect_adapter_for_pe). Sliding-window quality trimming',
        "from the 3' end, like Trimmomatic SLIDINGWINDOW:4:${params.trim_mean_quality}: --cut_right --cut_window_size 4 --cut_mean_quality ${params.trim_mean_quality}.",
        "Pairs with a read shorter than ${params.min_read_length} bp after trimming removed (--length_required ${params.min_read_length}).",
        'Otherwise fastp defaults: a read with > 40% bases below Q15 or > 5 N bases is removed with its mate;',
        'poly-G tails trimmed automatically for NovaSeq/NextSeq data.',
        "QC check: a sample losing > ${params.max_pct_removed_fastp}% of its pairs in fastp stops the run (at a minimum length > 36 bp).",
    ])
    if (!params.skip_bbduk) {
        addTool(ctx, 'BBDUK', 'BBDuk (BBTools)', [
            "bbduk.sh ref=${params.bbduk_ref} k=31 hdist=1: pairs with a k-mer match to PhiX or Illumina sequencing",
            'artifacts (BBTools bundled references) removed.',
        ])
    }
    if (!params.skip_rrna_removal) {
        addTool(ctx, 'RIBODETECTOR', 'RiboDetector', [
            "ribodetector_cpu -l ${params.ribodetector_read_len} -e rrna --chunk_size ${params.ribodetector_chunk_size}",
            'Deep-learning rRNA classifier, no reference database. -e rrna: a pair is rRNA if either mate is.',
            "rRNA pairs ${params.save_rrna ? 'kept in reads/rRNA/' : 'discarded'}; non-rRNA pairs are the clean reads.",
        ])
    } else {
        out << '' << 'rRNA removal skipped (--skip_rrna_removal).'
    }
    if (normalize) {
        addTool(ctx, 'BBNORM', 'BBNorm (BBTools)', [
            "bbnorm.sh target=${params.bbnorm_target} mindepth=${params.bbnorm_mindepth} prefilter=t (2 passes, default):",
            "k-mer depth normalization of the clean reads, used only as assembler input. mindepth lowered from 5 so",
            'reads of rare transcripts are kept. Quantification and reference mapping use all clean reads.',
        ])
    }
    addTool(ctx, 'MULTIQC', 'MultiQC', [
        'qc/multiqc/raw_multiqc_report.html: FastQC of raw reads and fastp reports.',
        'qc/multiqc/clean_multiqc_report.html: FastQC of clean' + (normalize ? ' and normalized' : '') + ' reads' +
            (params.skip_bbduk ? '.' : ', BBDuk logs.'),
    ])
    addScript(ctx, 'bin/summarize_qc.py', [
        'Per-sample read accounting, qc/preliminary_qc.tsv: raw pairs removed by fastp / BBDuk / rRNA / non-rRNA,',
        'as counts and % of raw pairs (summing to 100), pairs kept by BBNorm, fastp quality metrics.',
    ])

    // 2. Read screening
    if (params.run_read_profiling || params.stop_after_profiling) {
        addSection(ctx, 'Read screening')
        addTool(ctx, 'SOURMASH_SKETCH', 'sourmash', [
            'sourmash sketch dna -p k=31,scaled=1000,abund: both mates of the QC\'d reads in one signature per sample.',
            'sourmash compare -k 31: pairwise abundance-weighted (angular) similarity, sourmash_similarity.csv.',
        ])
        addTool(ctx, 'NONPAREIL', 'Nonpareil', ['nonpareil -T kmer on R1 of the QC\'d reads: coverage, redundancy, diversity.'])
        def strand_db = params.strand_check_dmnd ?: params.reference_dmnd
        if (strand_db) {
            addTool(ctx, 'STRANDEDNESS', 'DIAMOND (strandedness)', [
                "blastx --fast --evalue 1e-5 --max-target-seqs 1 of the first ${params.strand_check_reads} R1 reads against",
                "${strand_db}.",
                "Library strandedness from the reading frame of each read's best hit:",
                'fr if >= 80% of reads with a hit are sense, rf if <= 20%, otherwise unstranded (< 1000 hits: not called).',
            ])
        }
        if (params.kaiju_db) {
            addTool(ctx, 'KAIJU', 'Kaiju', ["kaiju -a greedy -e 5, read pairs against ${params.kaiju_db}; kaiju2table -r superkingdom."])
            manifestLines(params.kaiju_db, 'downloaded').each { line -> out << "    ${line}" }
        }
        addScript(ctx, 'bin/summarize_read_profile.py', [
            'One row per sample (qc/read_profile/read_profile.tsv). Suggested co-assembly groups: average-linkage',
            "clustering of samples on sourmash similarity, cut at ${params.profile_group_similarity}.",
        ])
        if (params.stop_after_profiling) {
            out << '' << 'Run stops after read screening (--stop_after_profiling).'
        }
    }

    def full_run = !params.stop_after_profiling

    // Assembly branch
    if (full_run && params.run_assembly) {
        def assemblers    = params.assemblers.tokenize(',')*.trim()
        def co_assemblers = params.coassembly_assemblers ? params.coassembly_assemblers.tokenize(',')*.trim() : []
        addSection(ctx, 'Assembly')
        ['rnaspades', 'megahit'].each { a ->
            def scope = [assemblers.contains(a) ? 'each sample' : null,
                         co_assemblers.contains(a) ? 'each samplesheet group (co-assembly)' : null].findAll { x -> x }
            if (scope) {
                def process = a == 'rnaspades' ? 'RNASPADES' : 'MEGAHIT'
                def label   = a == 'rnaspades' ? 'rnaSPAdes (SPAdes)' : 'MEGAHIT'
                addTool(ctx, process, label, ["Default settings; assembles ${scope.join(' and ')}" +
                    (normalize ? ' from BBNorm-normalized reads.' : ' from the clean reads.')])
            }
        }
        if (co_assemblers) {
            out << '' << 'Co-assembly: the assembler input of all samples in a samplesheet group (groups with > 1 sample) pooled.'
        }

        if (params.run_assembly_qc) {
            addSection(ctx, 'Assembly QC')
            addTool(ctx, 'SEQKIT_STATS', 'SeqKit', ['seqkit stats --all: contig number, length, N50 per assembly.'])
            addTool(ctx, 'ASSEMBLY_MAPRATE', 'salmon (read mapping rate)', [
                'Assembler input reads mapped back to each assembly: salmon index; salmon quant -l A --validateMappings.',
            ])
            if (params.busco_db_dir) {
                addTool(ctx, 'BUSCO', 'BUSCO', [
                    "busco -m transcriptome -l ${params.busco_lineages} --offline; datasets in ${params.busco_db_dir}.",
                    'Community-level score, for comparing assemblers and co-assembly on the same reads.',
                ])
            }
            addScript(ctx, 'bin/summarize_assembly_qc.py', ['One row per assembly: assembly/qc/assembly_qc.tsv.'])
        }

        addSection(ctx, 'Clustering and ORF prediction')
        addTool(ctx, 'MMSEQS_LINCLUST', 'MMseqs2', [
            "Contigs of all assemblies pooled; mmseqs easy-linclust --min-seq-id ${params.cluster_identity}",
            "-c ${params.cluster_coverage} --cov-mode 1. Cluster representatives are the reference contig set.",
        ])
        addTool(ctx, 'TRANSDECODER', 'TransDecoder', [
            "TransDecoder.LongOrfs -m ${params.transdecoder_min_len}; TransDecoder.Predict --single_best_only.",
            'Ab initio ORFs on the clustered contigs (assembly/transdecoder/).',
        ])
        if (params.metaeuk_db) {
            addTool(ctx, 'METAEUK', 'MetaEuk', [
                'metaeuk easy-predict -e 100 --metaeuk-eval 0.0001 --min-length 40 (MetaEuk README example settings):',
                "homology-based proteins on the clustered contigs against ${params.metaeuk_db} (assembly/metaeuk/).",
                'Runs alongside TransDecoder; both protein sets are kept and annotated separately.',
            ])
            manifestLines(params.metaeuk_db, '== outputs').each { line -> out << "    ${line}" }
        }

        if (params.reference_dmnd || params.eggnog_data_dir) {
            addSection(ctx, 'Annotation')
            if (params.reference_dmnd) {
                addTool(ctx, 'DIAMOND_BLASTP', 'DIAMOND', [
                    "blastp --sensitive --evalue ${params.diamond_evalue} --max-target-seqs 1: predicted proteins",
                    "(${params.metaeuk_db ? 'TransDecoder and MetaEuk, separately' : 'TransDecoder'}) against ${params.reference_dmnd}.",
                ])
            }
            if (params.eggnog_data_dir) {
                addTool(ctx, 'EGGNOG_MAPPER', 'eggNOG-mapper', [
                    "emapper.py -m diamond (default sensitivity), data in ${params.eggnog_data_dir};",
                    "${params.metaeuk_db ? 'TransDecoder and MetaEuk proteins annotated separately.' : 'TransDecoder proteins.'}",
                ])
            }
        }

        addSection(ctx, 'Quantification')
        addTool(ctx, 'SALMON_QUANT', 'salmon', [
            'salmon index on the clustered contigs; salmon quant -l A --validateMappings per sample.',
        ])
        addScript(ctx, 'bin/merge_salmon.py', ['Merged NumReads and TPM tables: assembly/tables/.'])
    }

    // Reference branch
    if (full_run && params.run_reference) {
        addSection(ctx, 'Reference mapping')
        addTool(ctx, 'DIAMOND_BLASTX', 'DIAMOND', [
            "blastx --fast --evalue ${params.diamond_evalue} --max-target-seqs 1: R1 and R2 of the QC'd reads",
            "searched as single reads against ${params.reference_dmnd}.",
        ])
        manifestLines(file(params.reference_dmnd).parent.toString(), '== merge counts').each { line -> out << "    ${line}" }
        addScript(ctx, 'bin/summarize_diamond_hits.py', [
            "Best-hit read counts per reference protein${params.reference_taxonomy ? ', with taxonomy from ' + params.reference_taxonomy : ''}.",
        ])
    }

    def methods = file("${params.outdir}/methods.txt")
    methods.parent.mkdirs()
    methods.text = out.join('\n') + '\n'
}
