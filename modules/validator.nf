#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process GFF_VALIDATOR {

    container 'quay.io/microbiome-informatics/virify-python3:1.2'

    input:
        path momo_gff

    script:
    if (momo_gff.size() > 0)
        """
        gt gff3validator ${momo_gff}
        """
}

