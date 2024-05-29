#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process AMRFINDER_REPORT {
    publishDir "$params.outdir/func_annot", mode: 'copy'

    container 'quay.io/biocontainers/python:3.9--1'

    input:
    path amrfinder_tsv
    path mobilome_gff
    path map_file
    path user_gff

    output:
    path("amr_location.tsv")

    script:
    def user_gff_arg = user_gff.name != "no_user_gff" ? "--user_gff ${user_gff}" : ""
    """    
    amr_report.py \
    --amr_out ${amrfinder_tsv} \
    --mobilome ${mobilome_gff} \
    --contigs_map ${map_file} $user_gff_arg
    """
}

