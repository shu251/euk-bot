// How well an assembly represents the reads it was built from: salmon mapping rate of those reads
// (the sample's reads for a per-sample assembly, the group's pooled reads for a co-assembly).
process ASSEMBLY_MAPRATE {
    tag "${contigs.name - '.contigs.fa'}"
    label 'process_medium'
    conda 'bioconda::salmon=1.10.3'
    container 'quay.io/biocontainers/salmon:1.10.3--h45fbf2d_5'

    input:
    tuple val(meta), path(contigs), path(reads)

    output:
    path("*.maprate.tsv"), emit: maprate

    script:
    def name = contigs.name - '.contigs.fa'
    """
    salmon index -t ${contigs} -i index -p ${task.cpus}
    salmon quant -i index -l A -1 ${reads[0]} -2 ${reads[1]} --validateMappings -p ${task.cpus} -o quant

    pct=\$(grep -o '"percent_mapped": [0-9.]*' quant/aux_info/meta_info.json | grep -o '[0-9.]*\$')
    n=\$(grep -o '"num_processed": [0-9]*' quant/aux_info/meta_info.json | grep -o '[0-9]*\$')
    printf 'assembly\\treads_processed\\tpct_reads_mapped\\n%s\\t%s\\t%s\\n' ${name} "\$n" "\$pct" > ${name}.maprate.tsv
    rm -rf index quant
    """

    stub:
    def name = contigs.name - '.contigs.fa'
    """
    printf 'assembly\\treads_processed\\tpct_reads_mapped\\n${name}\\t0\\t0\\n' > ${name}.maprate.tsv
    """
}
