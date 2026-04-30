process COMBINEREPORTER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12--1':
        'biocontainers/python:3.12--1'}"

    input:
    tuple val(meta), path(mobilome_gff), path(pathofact_gff), path(arg_gff), path(bgcs_gff), path(ips_tsv)

    output:
    tuple val(meta), path("*_combined_report.tsv"), optional: true, emit: tsv
    tuple val("${task.process}"), val('python'), eval("python --version"), topic: versions, emit: versions_combinereporter

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mobilome_param = mobilome_gff ? "--mobilome ${mobilome_gff}" : ""
    def pathofact_param = pathofact_gff ? "--pathofact2 ${pathofact_gff}" : ""
    def arg_param = arg_gff ? "--amr ${arg_gff}" : ""
    def bgcs_param = bgcs_gff ? "--bgc ${bgcs_gff}" : ""
    def ips_param = ips_tsv ? "--interproscan ${ips_tsv}" : ""
    """
    pathofact2_report.py \\
        ${mobilome_param} \\
        ${pathofact_param} \\
        ${arg_param} \\
        ${bgcs_param} \\
        ${ips_param} \\
        --output ${prefix}_combined_report.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_combined_report.tsv
    """
}
