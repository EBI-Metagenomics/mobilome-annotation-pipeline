process VALIDATE_ICE_ELEMENTS {
    tag "${meta.id}"
    label 'process_single'

    container 'quay.io/biocontainers/python:3.9--1'    
   
    input:
    tuple val(meta), path(macsyfinder_tsv), path(refined_tsv)
    
    output:
    tuple val(meta), path("*_refined.tsv") , emit: validated_ices
    tuple val(meta), path("*_rejected.tsv"), emit: rejected_ices, optional: true
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    validate_ices.py \\
        ${refined_tsv} \\
        --macsyfinder all_systems.tsv \\
        --output  ${prefix}_refined.tsv \\
        --rejected ${prefix}_rejected.tsv
    """
}
