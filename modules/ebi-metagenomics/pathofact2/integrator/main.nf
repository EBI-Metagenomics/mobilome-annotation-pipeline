process PATHOFACT2_INTEGRATOR {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'biocontainers/python:3.12.12'}"

    input:
    tuple val(meta), path(gff), path(prot_annot), path(pathofact_support), val(annot_type)

    output:
    tuple val(meta), path("*_pathofact2.gff"), optional: true, emit: gff
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //g'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('pathofact2'), eval("echo ${VERSION}"), topic: versions, emit: versions_pathofact2

    when:
    task.ext.when == null || task.ext.when

    script:
    VERSION = '1.0.4' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    pathofact2_integrator.py \\
        -g ${gff} \\
        -a ${prot_annot} \\
        -s ${pathofact_support} \\
        -t ${annot_type} \\
        -o ${prefix}_pathofact2.gff
    """

    stub:
    VERSION = '1.0.4' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_pathofact2.gff
    """
}
