process DB_UNTAR_AMRFINDERPLUS {
    label 'process_single'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--h36e9172_9'
        : 'biocontainers/gnu-wget:1.18--h36e9172_9'}"

    input:
    path tarball

    output:
    path "amrfinderdb/", emit: db
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir amrfinderdb
    tar -xzf ${tarball} -C amrfinderdb/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tar: \$(tar --version | head -1 | sed 's/tar (GNU tar) //')
    END_VERSIONS
    """

    stub:
    """
    mkdir amrfinderdb
    touch amrfinderdb/dummy.db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tar: stub
    END_VERSIONS
    """
}
