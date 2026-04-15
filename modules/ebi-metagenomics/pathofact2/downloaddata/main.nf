process PATHOFACT2_DOWNLOADDATA {
    tag '${task.process}'
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'biocontainers/python:3.12'}"

    input:
    val db_zenodo_id

    output:
    path "*.tar.gz", emit: zenodo_file
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //g'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def zenodo_id = db_zenodo_id ?: 18223764
    """
    download_zenodo.py \\
        ${args} \\
        --output-file Models.tar.gz \\
        $zenodo_id
    """

    stub:
    def args = task.ext.args ?: ''
    """
    echo $args
    echo "stub" | gzip > Models.tar.gz
    """
}
