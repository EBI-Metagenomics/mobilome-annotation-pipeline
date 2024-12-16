process GFF_REDUCE {

    publishDir "${params.outdir}/gff_output_files", mode: 'copy'

    container 'quay.io/biocontainers/python:3.9--1'

    input:
    tuple val(meta), path(mobilome_prokka_gff)

    output:
    tuple val(meta), path("mobilome_clean.gff"), emit: mobilome_clean
    tuple val(meta), path("mobilome_extra.gff"), emit: mobilome_extra
    tuple val(meta), path("mobilome_nogenes.gff"), emit: mobilome_nogenes

    script:
    """
    gff_minimizer.py --mobilome_prokka_gff ${mobilome_prokka_gff} --prefix ${meta.id}
    """
}
