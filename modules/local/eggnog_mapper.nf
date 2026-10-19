process EGGNOG_MAPPER {
    label 'process_high'
    conda 'bioconda::eggnog-mapper=2.1.12'

    input:
    path(pep)
    path(data_dir)

    output:
    path('clustered.emapper.annotations'), emit: annotations

    script:
    """
    emapper.py \\
        -i ${pep} --itype proteins \\
        --data_dir ${data_dir} \\
        -m diamond \\
        --cpu ${task.cpus} \\
        -o clustered
    """

    stub:
    """
    touch clustered.emapper.annotations
    """
}
