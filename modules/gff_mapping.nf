#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process GFF_MAPPING {
    publishDir "$launchDir/$params.outdir/gff_output_files"
    stageInMode = 'copy'

    container 'quay.io/microbiome-informatics/virify-python3:1.2'

    input:
        path mobilome_extra
	path user_gff

    output:
	path "user_mobilome_extra.gff"
	path "user_mobilome_full.gff"

    script:
    """
    gff_mapping.py \
    --mobilome_prokka_extra ${mobilome_extra} \
    --user_gff ${user_gff}
    """
}

