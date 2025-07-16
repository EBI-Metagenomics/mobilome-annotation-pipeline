process REFINE_BOUNDARIES {
    tag "${meta.id}"
    label 'process_single'

    container 'quay.io/biocontainers/biopython:1.81'
    
    input:
    tuple val(meta), path(boundaries_tsv), path(trna_gff), path(vmatch_tsv)
    
    output:
    tuple val(meta), path("*_refined_data.tsv"), emit: refined_tsv
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ice_boundary_refinement.py \\
        --macsyfinder ${boundaries_tsv} \\
        --trna ${trna_gff} \\
        --drs ${vmatch_tsv} \\
        --output ${prefix}_refined_data.tsv -v
    """
}
