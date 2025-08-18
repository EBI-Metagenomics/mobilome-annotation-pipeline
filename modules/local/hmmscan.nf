process HMMSCAN {
    tag "${meta.id}"
    label 'process_medium'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/hmmer:3.4--hdbdd923_1'
        : 'biocontainers/hmmer:3.4--hdbdd923_1'}"

    input:
    tuple val(meta), path(faa_file)
    tuple val(meta2), path(ice_hmm_models)

    output:
    tuple val(meta), path("*_prescan.tbl"), emit: hmmscan_tbl
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    hmmscan \\
        --tblout ${prefix}_prescan.tbl \\
        ${ice_hmm_models} \\
        ${faa_file}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmer: \$(hmmscan -h | grep -o '^# HMMER [0-9.]*' | sed 's/^# HMMER *//')
    END_VERSIONS
    """
}
