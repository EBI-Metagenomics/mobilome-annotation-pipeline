process GFF_REDUCE {
    tag "${meta.id}"
    label 'process_single'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.9--1'
        : 'biocontainers/python:3.9--1'}"

    input:
    tuple val(meta), path(mobilome_prokka_gff)

    output:
    tuple val(meta), path("*_mobilome_clean.gff"), emit: mobilome_clean
    tuple val(meta), path("*_mobilome_extra.gff"), emit: mobilome_extra
    tuple val(meta), path("*_mobilome_nogenes.gff"), emit: mobilome_nogenes
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gff_minimizer.py \\
        --mobilome_prokka_gff ${mobilome_prokka_gff} \\
        --prefix ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
