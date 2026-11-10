// Consensus taxonomy of the predicted proteins (EUKulele, mets mode) against the euk-bot reference, from the
// reference folder build_eukulele_db.sbatch makes (--eukulele_db). EUKulele's own DIAMOND step is skipped: the
// DIAMOND_BLASTP hits (top 100 per protein, the same -k as EUKulele) are placed where EUKulele writes its
// alignment, which it then reuses. BUSCO steps are off (--no_busco).
// The whole prot-map.json (one entry per reference protein) is loaded into memory, hence process_high_memory.
process EUKULELE {
    tag "${name}"
    label 'process_high'
    label 'process_high_memory'
    conda 'bioconda::eukulele=2.1.2'
    container 'quay.io/biocontainers/eukulele:2.1.2--pyhdfd78af_0'

    input:
    tuple val(name), path(pep), path(hits)    // name: ORF caller (transdecoder / metaeuk)
    path(db_dir)

    output:
    tuple val(name), path("${name}.eukulele_taxonomy.tsv"), emit: taxonomy
    path("${name}.eukulele_logs"),                         emit: logs

    script:
    def args = task.ext.args ?: ''
    """
    mkdir -p samples eukulele/mets_full/diamond
    ln -s ../${pep} samples/${name}.faa
    ref_fasta=\$(basename ${db_dir}/diamond/*.dmnd .dmnd)

    if [[ ! -s ${hits} ]]; then
        # No hits: EUKulele would otherwise run its own DIAMOND search
        printf '\\ttranscript_name\\tclassification_level\\tfull_classification\\tclassification\\tmax_score\\tambiguous\\tdatabase\\n' \\
            > ${name}.eukulele_taxonomy.tsv
        mkdir -p ${name}.eukulele_logs
        exit 0
    fi
    ln -s \$PWD/${hits} eukulele/mets_full/diamond/${name}.diamond.out

    EUKulele all \\
        --mets_or_mags mets \\
        --sample_dir samples \\
        --protein_extension .faa \\
        --reference_dir ${db_dir} \\
        --ref_fasta \${ref_fasta} \\
        --database euk-bot-ref \\
        --out_dir eukulele \\
        --scratch scratch \\
        --no_busco \\
        --CPUs ${task.cpus} \\
        ${args}

    mv eukulele/log ${name}.eukulele_logs
    # EUKulele reports a failed taxonomy step in its logs but still exits 0
    est=eukulele/taxonomy_estimation/${name}-estimated-taxonomy.out
    if [[ ! -s \${est} ]]; then
        echo "EUKulele wrote no taxonomy estimate; see ${name}.eukulele_logs/" >&2
        cat ${name}.eukulele_logs/tax_est_*.err >&2 || true
        exit 1
    fi
    mv \${est} ${name}.eukulele_taxonomy.tsv
    """

    stub:
    """
    touch ${name}.eukulele_taxonomy.tsv
    mkdir -p ${name}.eukulele_logs
    """
}
