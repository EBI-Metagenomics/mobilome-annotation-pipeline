process VIRIFY_QC {
    tag "${meta.id}"

    container 'quay.io/biocontainers/python:3.9--1'

    input:
    tuple val(meta), path(virify_gff)

    output:
    tuple val(meta), path("*_virify_hq.gff"), emit: virify_hq

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    virify_qc.py \\
        --virify_gff ${virify_gff} \\
        --prefix ${prefix}
    """
}
