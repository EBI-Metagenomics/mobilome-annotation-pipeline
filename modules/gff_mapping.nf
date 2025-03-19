process GFF_MAPPING {
    tag "$meta.id"
    label 'process_single'

    // TODO: add the singularity container
    container 'quay.io/biocontainers/python:3.9--1'

    input:
    tuple val(meta), path(mobilome_clean), path(user_gff)

    output:
    tuple val(meta), path("${meta.id}_user_mobilome_clean.gff"), optional: true, emit: mobilome_clean_gff
    tuple val(meta), path("${meta.id}_user_mobilome_extra.gff"), optional: true, emit: mobilome_extra_gff
    tuple val(meta), path("${meta.id}_user_mobilome_full.gff"),  optional: true, emit: mobilome_full_gff

    script:
    def user_proteins_arg = (user_gff) ? "--user_gff ${user_gff}" : ""
    """
    gff_mapping.py \\
    --prefix ${meta.id} \\
    --mobilome_clean ${mobilome_clean} ${user_proteins_arg}
    """
}
