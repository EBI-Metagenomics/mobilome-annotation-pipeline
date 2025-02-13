#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process GFF_VALIDATOR {
    tag "$meta.id"
    label 'process_single'

    container 'quay.io/biocontainers/genometools-genometools:1.6.5--py310h3db02ab_0'

    input:
    path momo_gff

    script:
    """
    gt gff3validator ${momo_gff}
    """
}

