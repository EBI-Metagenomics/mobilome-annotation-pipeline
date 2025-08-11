process AMRFINDER_PLUS {
    tag "$meta.id"
    label 'process_high'

    container 'quay.io/biocontainers/ncbi-amrfinderplus:3.11.4--h6e70893_0'

    input:
    tuple val(meta), path(fna), path(faa), path(gff)

    output:
    tuple val(meta), path("*_amrfinderplus.tsv"), emit: amrfinder_tsv

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """    
    amrfinder --plus \\
        -n ${fna} \\
        -p ${faa} \\
        -g ${gff} \\
        -d ${params.amrfinder_plus_db} \\
        -a prokka \\
        --output ${prefix}_amrfinderplus.tsv \\
        --threads ${task.cpus}
    """
}
