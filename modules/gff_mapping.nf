#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process GFF_MAPPING {
    publishDir "$params.outdir/gff_output_files", mode: 'copy'

    container 'quay.io/biocontainers/python:3.9--1'

    input:
    path mobilome_clean
    path user_gff

    output:
    path "user_mobilome_extra.gff"
    path "user_mobilome_full.gff"
    path "user_mobilome_clean.gff"

    script:
    """
    gff_mapping.py \
    --mobilome_clean ${mobilome_clean} \
    --user_gff ${user_gff}
    """
}

