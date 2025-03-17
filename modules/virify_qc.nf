
process VIRIFY_QC {

    tag "${meta.id}"

    container 'quay.io/biocontainers/python:3.9--1'

    input:
    tuple val(meta), path(virify_gff)

    output:
    tuple val(meta), path("${meta.id}_virify_hq.gff"), emit: virify_hq

    script:
    """
    virify_qc.py \\
    --virify_gff ${virify_gff} \\
    --prefix ${meta.id}
    """
}
