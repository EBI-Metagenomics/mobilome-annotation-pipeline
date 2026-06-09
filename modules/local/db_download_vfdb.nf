process DB_DOWNLOAD_VFDB {
    label 'process_medium'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/diamond:2.1.16--h13889ed_0'
        : 'biocontainers/diamond:2.1.16--h13889ed_0'}"

    output:
    path "VFDB_setB_pro.dmnd", emit: db
    path "versions.yml",       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    curl -L -o VFDB_setB_pro.fas.gz https://www.mgc.ac.cn/VFs/Down/VFDB_setB_pro.fas.gz

    diamond makedb \\
        --threads ${task.cpus} \\
        --in VFDB_setB_pro.fas.gz \\
        -d VFDB_setB_pro

    rm VFDB_setB_pro.fas.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: \$(diamond --version | sed 's/diamond version //g')
    END_VERSIONS
    """

    stub:
    """
    touch VFDB_setB_pro.dmnd

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: stub
    END_VERSIONS
    """
}
