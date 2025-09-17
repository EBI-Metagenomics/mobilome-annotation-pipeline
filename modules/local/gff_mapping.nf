process GFF_MAPPING {
    tag "${meta.id}"
    label 'process_single'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.9--1'
        : 'biocontainers/python:3.9--1'}"

    input:
    tuple val(meta), path(mobilome_gff), path(map_file), path(cds_gff), val(gff_type)

    output:
    tuple val(meta), path("*_mobilome_clean.gff"), optional: true, emit: mobilome_clean_gff
    tuple val(meta), path("*_mobilome_extra.gff"), optional: true, emit: mobilome_extra_gff
    tuple val(meta), path("*_mobilome_full.gff"),  optional: true, emit: mobilome_full_gff
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gff_mapping.py \\
        --mobilome_gff ${mobilome_gff} \\
        --cds_gff ${cds_gff} \\
        --map_file ${map_file} \\
        --prefix ${prefix}_${gff_type} \\
        --mode $gff_type

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
