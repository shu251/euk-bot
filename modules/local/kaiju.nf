// Read-level composition (protein-level, greedy mode) against a Kaiju database; summarized to
// superkingdom (Eukaryota / Bacteria / Archaea / Viruses / unclassified).
process KAIJU {
    tag "${meta.id}"
    label 'process_high'
    label 'process_high_memory'   // nr_euk index needs ~204 GB RAM
    conda 'bioconda::kaiju=1.10.3'
    container 'quay.io/biocontainers/kaiju:1.10.3--h5814d7d_0'

    input:
    tuple val(meta), path(reads)
    path(db)    // directory with one *.fmi index plus nodes.dmp and names.dmp

    output:
    path("${meta.id}.kaiju_superkingdom.tsv"), emit: summary

    script:
    """
    fmi=\$(ls ${db}/*.fmi | head -1)
    kaiju -t ${db}/nodes.dmp -f \$fmi -i ${reads[0]} -j ${reads[1]} -z ${task.cpus} -a greedy -e 5 -o ${meta.id}.kaiju.out
    # kaiju2table 1.10 no longer accepts -r superkingdom; roll up to superkingdom via nodes.dmp instead
    kaiju_superkingdom.sh ${db}/nodes.dmp ${meta.id}.kaiju.out > ${meta.id}.kaiju_superkingdom.tsv
    rm ${meta.id}.kaiju.out
    """

    stub:
    """
    printf 'file\\tpercent\\treads\\ttaxon_id\\ttaxon_name\\n' > ${meta.id}.kaiju_superkingdom.tsv
    """
}
