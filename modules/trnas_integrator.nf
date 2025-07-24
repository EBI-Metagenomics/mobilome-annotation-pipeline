process TRNAS_INTEGRATOR {
    tag "$meta.id"
    label 'process_single'

    container 'quay.io/biocontainers/biopython:1.81'

    input:
    tuple val(meta), path(aragorn_file), path(prodigal_gff)

    output:
    tuple val(meta), path("*_merged.gff"), emit: merged_gff
    tuple val(meta), path("names.map") , emit: merged_map


    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    trnas_integrator.py \\
        ${prodigal_gff} \\
        ${aragorn_file} \\
        --output-gff ${prefix}_merged.gff \\
        --locus-tag-prefix ${prefix}
    """
}
