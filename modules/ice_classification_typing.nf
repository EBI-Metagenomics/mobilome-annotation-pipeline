process ICE_CLASSIFICATION_TYPING {
    tag "$meta.id"
    
    publishDir "${params.outdir}/ice_classification", mode: 'copy'
    
    input:
    tuple val(meta), path(integrated_ice_info)
    tuple val(meta), path(macsyfinder_gff3)
    
    output:
    tuple val(meta), path("${prefix}_classified_ice_elements.tsv"), emit: classified_ice_elements
    tuple val(meta), path("${prefix}_ice_typing_report.txt"), emit: typing_report
    path "versions.yml", emit: versions
    
    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    classify_ice_elements.py \\
        ${integrated_ice_info} \\
        ${prefix} \\
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
    touch ${prefix}_classified_ice_elements.tsv
    touch ${prefix}_ice_typing_report.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g' || echo "3.9.0")
        pandas: \$(echo "1.5.3")
    END_VERSIONS
    """
}
