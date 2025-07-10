process ICE_INFORMATION_INTEGRATION {
    tag "$meta.id"
    
    publishDir "${params.outdir}/ice_integration", mode: 'copy'

    conda "bioconda::biopython=1.81 conda-forge::pandas=1.5.3"

    
    input:
    tuple val(meta), path(validated_ice_elements)
    tuple val(meta), path(orit_results)
    tuple val(meta), path(ice_sequences)
    path genome_fasta
    tuple val(meta), path(macsyfinder_gff3)
    tuple val(meta), path(trna_annotations)
    tuple val(meta), path(direct_repeats)
    
    output:
    tuple val(meta), path("${prefix}_integrated_ice_info.json"), emit: integrated_ice_info
    tuple val(meta), path("${prefix}_ice_characteristics.tsv"), emit: ice_characteristics
    path "versions.yml", emit: versions
    
    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def gc_window_size = params.gc_window_size ?: 1000
    def gc_step_size = params.gc_step_size ?: 500
    
    """
    integrate_ice_information.py \\
        ${validated_ice_elements} \\
        ${macsyfinder_gff3} \\
        ${trna_annotations} \\
        ${direct_repeats} \\
        ${orit_results} \\
        ${genome_fasta} \\
        ${prefix} \\
        --gc-window-size ${gc_window_size} \\
        --gc-step-size ${gc_step_size} \\
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
    
    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_integrated_ice_info.json
    touch ${prefix}_ice_characteristics.tsv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g' || echo "3.9.0")
        pandas: \$(echo "1.5.3")
        biopython: \$(echo "1.81")
    END_VERSIONS
    """
}
