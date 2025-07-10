process ICE_QUALITY_ASSESSMENT {
    tag "$meta.id"
    
    publishDir "${params.outdir}/ice_quality_assessment", mode: 'copy'
    
    input:
    tuple val(meta), path(classified_ice_elements)
    tuple val(meta), path(integrated_ice_info)
    
    output:
    tuple val(meta), path("${prefix}_quality_assessed_ice.tsv"), emit: quality_assessed_ice
    tuple val(meta), path("${prefix}_quality_report.txt"), emit: quality_report
    tuple val(meta), path("${prefix}_failed_ice_elements.tsv"), emit: failed_elements
    path "versions.yml", emit: versions
    
    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def min_ice_size = params.min_ice_size ?: 10000
    def max_ice_size = params.max_ice_size ?: 500000
    def min_trna_count = params.min_trna_count ?: 1
    def min_confidence_score = params.min_confidence_score ?: 0.5
    def require_integrase = params.require_integrase ?: true
    def require_mobilization = params.require_mobilization ?: false
    
    """
    assess_ice_quality.py \\
        ${classified_ice_elements} \\
        ${integrated_ice_info} \\
        ${prefix} \\
        --min-size ${min_ice_size} \\
        --max-size ${max_ice_size} \\
        --min-trna-count ${min_trna_count} \\
        --min-confidence-score ${min_confidence_score} \\
        ${require_integrase ? '--require-integrase' : ''} \\
        ${require_mobilization ? '--require-mobilization' : ''} \\
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
    
    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_quality_assessed_ice.tsv
    touch ${prefix}_quality_report.txt
    touch ${prefix}_failed_ice_elements.tsv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g' || echo "3.9.0")
        pandas: \$(echo "1.5.3")
    END_VERSIONS
    """
}


