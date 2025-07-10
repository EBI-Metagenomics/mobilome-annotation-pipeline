process DEFENSEFINDER {
    conda "bioconda::defense-finder=2.0.1"
    
    input:
    tuple val(meta), path(input_file)
    path models_dir
    
    output:
    tuple val(meta), path("${meta.id}_defense_finder_systems.tsv"), emit: systems
    tuple val(meta), path("${meta.id}_defense_finder_genes.tsv"), emit: genes
    tuple val(meta), path("${meta.id}_defense_finder_hmmer.tsv"), emit: hmmer
    path "versions.yml", emit: versions
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def models_arg = models_dir ? "--models-dir ${models_dir}" : ""
    """
    # Create output directory
    mkdir -p defensefinder_output
    
    # Run DefenseFinder
    defense-finder run \\
        --out-dir defensefinder_output \\
        --workers ${task.cpus} \\
        --coverage 0.4 \\
        --db-type ordered_replicon \\
        ${models_arg} \\
        ${args} \\
        ${input_file}
    
    # Rename output files with sample prefix
    if [ -f "defensefinder_output/defense_finder_systems.tsv" ]; then
        mv defensefinder_output/defense_finder_systems.tsv ${prefix}_defense_finder_systems.tsv
    else
        touch ${prefix}_defense_finder_systems.tsv
    fi
    
    if [ -f "defensefinder_output/defense_finder_genes.tsv" ]; then
        mv defensefinder_output/defense_finder_genes.tsv ${prefix}_defense_finder_genes.tsv
    else
        touch ${prefix}_defense_finder_genes.tsv
    fi
    
    if [ -f "defensefinder_output/defense_finder_hmmer.tsv" ]; then
        mv defensefinder_output/defense_finder_hmmer.tsv ${prefix}_defense_finder_hmmer.tsv
    else
        touch ${prefix}_defense_finder_hmmer.tsv
    fi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        defense-finder: \$(defense-finder --version | grep -oP 'defense-finder \\K[0-9.]+')
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_defense_finder_systems.tsv
    touch ${prefix}_defense_finder_genes.tsv
    touch ${prefix}_defense_finder_hmmer.tsv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        defense-finder: \$(defense-finder --version | grep -oP 'defense-finder \\K[0-9.]+')
    END_VERSIONS
    """
}
