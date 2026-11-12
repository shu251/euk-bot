// --dbmem loads eggnog.db (~45 GB) into memory for the annotation phase: test3 ran at 8.5 proteins/s without it
// (1.9M TransDecoder proteins would take ~60 h), vs ~500/s with it in annotate_reference_db.sbatch.
process EGGNOG_MAPPER {
    tag "${name}"
    label 'process_high'
    conda 'bioconda::eggnog-mapper=2.1.12'
    container 'quay.io/biocontainers/eggnog-mapper:2.1.12--pyhdfd78af_2'

    input:
    tuple val(name), path(pep)    // name: ORF caller (transdecoder / metaeuk)
    path(data_dir)

    output:
    tuple val(name), path("${name}.emapper.annotations"), emit: annotations

    script:
    """
    emapper.py \\
        -i ${pep} --itype proteins \\
        --data_dir ${data_dir} \\
        -m diamond \\
        --dbmem \\
        --cpu ${task.cpus} \\
        -o ${name}
    """

    stub:
    """
    touch ${name}.emapper.annotations
    """
}
