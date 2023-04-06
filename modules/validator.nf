#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process gff_validator {
    memory "4 GB"
    cpus 1

    container "quay.io/microbiome-informatics/virify-python3:1.2"

    input:
        path momo_gff, name: 'momofy_predictions.gff'

    script:
    if(momo_gff.size() > 0)
        """
	gt gff3validator momofy_predictions.gff
        """
}

