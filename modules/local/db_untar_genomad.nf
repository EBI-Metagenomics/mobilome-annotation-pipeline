process DB_UNTAR_GENOMAD {
    label 'process_single'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--h36e9172_9'
        : 'biocontainers/gnu-wget:1.18--h36e9172_9'}"

    input:
    path tarball

    output:
    path "genomad_db_v1.9/", emit: db
    path "versions.yml",     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    tar -xzf ${tarball}
    # The tarball's top-level directory is genomad_db/; only the archive name is versioned.
    mv genomad_db genomad_db_v1.9

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tar: \$(tar --version | head -1 | sed 's/tar (GNU tar) //')
    END_VERSIONS
    """

    stub:
    """
    mkdir genomad_db_v1.9
    touch genomad_db_v1.9/genomad_db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tar: stub
    END_VERSIONS
    """
}
