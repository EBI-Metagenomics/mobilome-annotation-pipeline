process VALIDATE_ICE_ELEMENTS {
    tag "ice_validation"
    
    publishDir "${params.outdir}/final_results", mode: 'copy'

    conda "conda-forge::pandas=1.5.3"
    
    input:
    path refined_boundaries
    path trna_annotations
    
    output:
    path "validated_ice_elements.tsv", emit: validated_ice
    path "ice_summary_report.txt", emit: summary
    path "versions.yml", emit: versions
    
    script:
    def args = task.ext.args ?: ''
    def min_ice_size = params.min_ice_size ?: 10000
    def max_ice_size = params.max_ice_size ?: 500000
    def min_trna_count = params.min_trna_count ?: 1
    """
    validate_ice_elements.py \\
        ${refined_boundaries} \\
        ${trna_annotations} \\
        validated_ice_elements.tsv \\
        ice_summary_report.txt \\
        --min-ice-size ${min_ice_size} \\
        --max-ice-size ${max_ice_size} \\
        --min-trna-count ${min_trna_count} \\
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
    
    stub:
    """
    touch validated_ice_elements.tsv
    touch ice_summary_report.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """

}
