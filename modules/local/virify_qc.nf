process VIRIFY_QC {
    tag "${meta.id}"

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.9--1'
        : 'biocontainers/python:3.9--1'}"

    input:
    tuple val(meta), path(virify_gff)

    output:
    tuple val(meta), path("*_virify_hq.gff"), emit: virify_hq
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    virify_qc.py \\
        --virify_gff ${virify_gff} \\
        --prefix ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
