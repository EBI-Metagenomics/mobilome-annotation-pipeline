#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process GFF_REDUCE {
    publishDir "$launchDir/$params.outdir/gff_output_files"
    stageInMode = 'copy'

    container 'quay.io/microbiome-informatics/virify-python3:1.2'

    input:
        path mobilome_prokka_gff

    output:
        path "mobilome_clean.gff", emit: mobilome_clean
	path "mobilome_extra.gff", emit: mobilome_extra
	path "mobilome_nogenes.gff", emit: mobilome_nogenes

    script:
    """
    gff_minimizer.py --mobilome_prokka_gff ${mobilome_prokka_gff}
    """
}

