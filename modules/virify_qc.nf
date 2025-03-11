
process VIRIFY_QC {

    // FIXME: this module doesn't produce a valid GFF
    // gt gff3validator: error: line 1 in file "virify_hq.gff" does not begin with "##gff-version" or "##gff-version"
    publishDir "$params.outdir/prediction", mode: 'copy'

    container 'quay.io/biocontainers/python:3.9--1'

    input:
    tuple val(meta), path(virify_gff), path(virify_checkv)

    output:
    tuple val(meta), path("virify_hq.gff"), emit: virify_hq

    script:
    """
    virify_qc.py \
    --virify_gff ${virify_gff} \
    --checkv_out ${virify_checkv.join(' ')}
    """
}
