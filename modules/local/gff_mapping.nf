process GFF_MAPPING {
    tag "${meta.id}"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.9--1'
        : 'biocontainers/python:3.9--1'}"

    input:
    tuple val(meta), path(mobilome_gff), path(genes_gff), path(combined_report), val(user_proteins)

    output:
    tuple val(meta), path("*_mobilome_clean.gff"), optional: true, emit: mobilome_clean_gff
    tuple val(meta), path("*_mobilome_extra.gff"), optional: true, emit: mobilome_extra_gff
    tuple val(meta), path("*_mobilome_full.gff"),  optional: true, emit: mobilome_full_gff
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def genes_arg = genes_gff ? "--user_gff ${genes_gff}" : ""
    def combined_report_arg = combined_report ? "--combined_report ${combined_report}" : ""
    def user_proteins_arg = user_proteins ? "--user_proteins" : ""
    """
    gff_mapping.py \\
        --prefix ${prefix} \\
        --mobilome_gff ${mobilome_gff} \\
        ${genes_arg} \\
        ${combined_report_arg} \\
        ${user_proteins_arg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
