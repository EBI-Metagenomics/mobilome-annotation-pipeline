#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process VIRIFY_QC {

    // FIXME: this module doesn't produce a valid GFF
    // gt gff3validator: error: line 1 in file "virify_hq.gff" does not begin with "##gff-version" or "##gff-version"
    publishDir "$params.outdir/prediction", mode: 'copy'

    container 'quay.io/biocontainers/python:3.9--1'

    input:
    path vir_gff
    path vir_checkv

    output:
    path("virify_hq.gff"), emit: virify_hq

    script:
    if(vir_gff.exists())
        """
        virify_qc.py \
        --virify_gff ${vir_gff} \
        --checkv_out ${vir_checkv.join(' ')}
        """
    else
        """
        echo 'No input files for VIRify parsing... generating dummy files'
        touch virify_hq.gff
        """
}
