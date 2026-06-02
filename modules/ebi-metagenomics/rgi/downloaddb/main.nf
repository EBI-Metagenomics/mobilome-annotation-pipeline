process RGI_DOWNLOADDB {
    label 'process_single'

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--h36e9172_9'
        : 'biocontainers/gnu-wget:1.18--h36e9172_9'}"

    output:
    path "card_dir/"   , emit: card_json
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir card_dir
    wget https://card.mcmaster.ca/latest/data
    tar -xvf data ./card.json
    mv card.json card_dir

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget --version | head -1 | cut -d ' ' -f 3)
        untar: \$(echo \$(tar --version 2>&1) | sed 's/^.*(GNU tar) //; s/ Copyright.*\$//')
    END_VERSIONS

    """
    stub:
    """
    mkdir card_dir
    touch card_dir/card.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget --version | head -1 | cut -d ' ' -f 3)
        untar: \$(echo \$(tar --version 2>&1) | sed 's/^.*(GNU tar) //; s/ Copyright.*\$//')
    END_VERSIONS
    """
}
