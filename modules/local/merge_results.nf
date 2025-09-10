process MERGE_RESULTS {
    tag "${meta.id}"
    label 'process_single'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.9--1'
        : 'biocontainers/python:3.9--1'}"

    input:
    tuple val(meta), path(bed_files)

    output:
    tuple val(meta), path("*.merged.bed"), emit: bed
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    merge_results.py \\
        ${prefix}.merged.bed \\
        ${bed_files} \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
