process AMRFINDER_REPORT {
    tag "$meta.id"
    label 'process_single'

    container 'quay.io/biocontainers/python:3.9--1'

    input:
    tuple val(meta), path(amrfinder_tsv), path(mobilome_gff), path(map_file), path(user_gff)

    output:
    tuple val(meta), path("${meta.id}_amr_location.tsv")

    script:
    def user_gff_arg = user_gff ? "--user_gff ${user_gff}" : ""
    """    
    amr_report.py \\
        --amr_out ${amrfinder_tsv} \\
        --mobilome ${mobilome_gff} \\
        --contigs_map ${map_file} ${user_gff_arg} \\
        --output ${meta.id}_amr_location.tsv
    """
}
