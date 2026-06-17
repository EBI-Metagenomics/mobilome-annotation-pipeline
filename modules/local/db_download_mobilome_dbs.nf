process DB_DOWNLOAD_MOBILOME_DBS {
    label 'process_single'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--h36e9172_9'
        : 'biocontainers/gnu-wget:1.18--h36e9172_9'}"

    output:
    path "genomad_db_v1.9/", emit: genomad_db
    path "icf2_dbs/",        emit: icefinder_dbs
    path "versions.yml",     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    wget https://zenodo.org/records/14886553/files/genomad_db_v1.9.tar.gz
    tar -xzf genomad_db_v1.9.tar.gz
    rm genomad_db_v1.9.tar.gz

    wget ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/pipelines/tool-dbs/icefinder2lite/icf2_dbs.tar.gz
    tar -xzf icf2_dbs.tar.gz
    rm icf2_dbs.tar.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget --version | head -1 | cut -d ' ' -f 3)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p genomad_db_v1.9
    mkdir -p icf2_dbs/macsydata icf2_dbs/icehmm icf2_dbs/icefinder_prokka_uniprot
    touch icf2_dbs/icehmm/icescan.hmm

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: 1.18
    END_VERSIONS
    """
}
