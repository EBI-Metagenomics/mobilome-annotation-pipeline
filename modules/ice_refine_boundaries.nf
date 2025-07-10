process REFINE_BOUNDARIES {
    tag "${meta.id}"
    

    conda "conda-forge::pandas=1.5.3"
    
    input:
    tuple val(meta), path(boundaries_tsv), path(trna_gff), path(dr_tsv)
    
    output:
    path "refined_boundaries.tsv", emit: refined_boundaries
    path "boundary_refinement_log.txt", emit: refinement_log
    path "versions.yml", emit: versions
    
    script:
    def args = task.ext.args ?: ''
    def min_ice_size = params.min_ice_size ?: 10000
    def max_ice_size = params.max_ice_size ?: 500000
    """
    refine_ice_boundaries.py \\
        ${ice_coordinates} \\
        ${trna_annotations} \\
        ${direct_repeats} \\
        refined_boundaries.tsv \\
        boundary_refinement_log.txt \\
        --min-ice-size ${min_ice_size} \\
        --max-ice-size ${max_ice_size} \\
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
    
    stub:
    """
    touch refined_boundaries.tsv
    touch boundary_refinement_log.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}
