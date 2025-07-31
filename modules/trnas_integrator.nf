process TRNAS_INTEGRATOR {
    tag "$meta.id"
    label 'process_single'

    container 'quay.io/biocontainers/biopython:1.81'

    input:
    tuple val(meta), path(aragorn_file), path(prodigal_gff), path(prodigal_faa)

    output:
    tuple val(meta), path("*_merged.gff") , emit: merged_gff
    tuple val(meta), path("*_renamed.faa"), emit: merged_faa


    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    trnas_integrator.py \\
        --prodigal_gff ${prodigal_gff} \\
        --aragorn ${aragorn_file} \\
        --prodigal_faa ${prodigal_faa} \\
        --prefix ${prefix}
    """
}
