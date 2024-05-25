#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process AMRFINDER_REPORT {
    publishDir "$params.outdir/func_annot", mode: 'copy'

    container 'quay.io/biocontainers/python:3.9--1'

    input:
    path amrfinder_tsv
    path mobilome_gff
    path map_file

    output:
    path("amr_location.txt")

    script:
    """    
    amr_report.py \
    --amr_out ${amrfinder_tsv} \
    --mobilome ${mobilome_gff} \
    --contigs_map ${map_file}
    """
}

