process GFF_MAPPING {
    tag "$meta.id"
    label 'process_single'

    // TODO: add the singularity container
    container 'quay.io/biocontainers/python:3.9--1'

    input:
    tuple val(meta), path(mobilome_clean), path(user_gff)

    output:
    tuple val(meta), path("${meta.id}_user_mobilome_extra.gff"), emit: mobilome_extra_gff
    tuple val(meta), path("${meta.id}_user_mobilome_full.gff"),  emit: mobilome_full_gff
    tuple val(meta), path("${meta.id}_user_mobilome_clean.gff"), emit: mobilome_clean_gff

    script:
    """
    gff_mapping.py \
    --prefix ${meta.id} \
    --mobilome_clean ${mobilome_clean} \
    --user_gff ${user_gff}
    """
}
