process ZENODO_DOWNLOAD {
    tag "${meta.id}"
    label 'process_single'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:1.5.4--pyhdfd78af_0'
        : 'quay.io/microbiome-informatics/mgnify-pipelines-toolkit:1.5.4'}"

    input:
    // filename: the record file's name on Zenodo. download_zenodo selects files by index
    // (--file-index via ext.args), so declaring the expected name as the output makes the
    // task fail loudly if the index ever points at a different file.
    tuple val(meta), val(zenodo_id), val(filename)

    output:
    tuple val(meta), path("${filename}"), emit: file
    tuple val("${task.process}"), val('mgnify-pipelines-toolkit'), eval("get_mpt_version"), topic: versions, emit: versions_mgnifypipelinestoolkit

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    download_zenodo \\
        ${args} \\
        ${zenodo_id}
    """

    stub:
    // A real (empty-content) tarball for .tar.gz/.tgz targets, so downstream untar steps also stub-run.
    def is_tarball = filename.endsWith('.tar.gz') || filename.endsWith('.tgz')
    """
    if [ "${is_tarball}" == "true" ]; then
        mkdir stub && touch stub/stub.txt
        tar -czf ${filename} stub
    else
        echo "stub" | gzip > ${filename}
    fi
    """
}
