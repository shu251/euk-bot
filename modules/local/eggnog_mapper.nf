process EGGNOG_MAPPER {
    tag "${name}"
    label 'process_high'
    conda 'bioconda::eggnog-mapper=2.1.12'
    container 'quay.io/biocontainers/eggnog-mapper:2.1.12--pyhdfd78af_2'

    input:
    tuple val(name), path(pep)    // name: ORF caller (transdecoder / metaeuk)
    path(data_dir)

    output:
    path("${name}.emapper.annotations"), emit: annotations

    script:
    """
    emapper.py \\
        -i ${pep} --itype proteins \\
        --data_dir ${data_dir} \\
        -m diamond \\
        --cpu ${task.cpus} \\
        -o ${name}
    """

    stub:
    """
    touch ${name}.emapper.annotations
    """
}
