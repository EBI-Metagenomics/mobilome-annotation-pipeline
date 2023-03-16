#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process momo_validator {
    publishDir "$launchDir/MoMofy_results", mode: 'copy'
    stageInMode = 'copy'

    memory "4 GB"
    cpus 4

    container "quay.io/microbiome-informatics/virify-python3:1.2"

    input:
        path momo_gff, name: 'momofy_predictions.gff'

    output:
        path 'validated_momofy_predictions.gff'

    script:
    if(momo_gff.size() > 0)
        """
	gt gff3 \
        -retainids yes \
        -checkids yes \
        -o validated_momofy_predictions.gff \
        momofy_predictions.gff
        """
}

