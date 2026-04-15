process COMBINEREPORTER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1':
        'biocontainers/python:3.9--1'}"

    input:
    tuple val(meta), path(mobilome_gff), path(pathofact_gff), path(arg_gff), path(bgcs_gff)

    output:
    tuple val(meta), path("*_combined_report.tsv"), emit: tsv
    tuple val("${task.process}"), val('python'), eval("python --version"), topic: versions, emit: versions_combinereporter

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    combinereporter.py \\
        $args \\
        --mobilome ${mobilome} \\
        --pathofact_gff ${pathofact_gff} \\

        --ouput ${prefix}_combined_report.tsv
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args
    
    touch ${prefix}_combined_report.tsv
    """
}
