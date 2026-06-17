process DB_UNTAR_PATHOFACT {
    label 'process_single'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--h36e9172_9'
        : 'biocontainers/gnu-wget:1.18--h36e9172_9'}"

    input:
    path tarball

    output:
    path "Models/",      emit: db
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    tar -xzf ${tarball}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tar: \$(tar --version | head -1 | sed 's/tar (GNU tar) //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p Models/TOX Models/VF
    touch Models/TOX/final_model.joblib
    touch Models/VF/final_model.joblib

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tar: stub
    END_VERSIONS
    """
}
