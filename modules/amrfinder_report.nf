#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process AMRFINDER_REPORT {
    publishDir "$params.outdir/func_annot", mode: 'copy'

    container 'quay.io/microbiome-informatics/virify-python3:1.2'

    input:
        path amrfinder_tsv
	path mobilome_gff

    output:
	path("amr_location.txt")

    when:
	amrfinder_tsv.fileExists()

    script:
    """    
    amr_report.py \
    --amr_out ${amrfinder_tsv} \
    --mobilome ${mobilome_gff} \
    """
}

