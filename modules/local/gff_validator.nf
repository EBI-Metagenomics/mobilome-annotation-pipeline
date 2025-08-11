process GFF_VALIDATOR {
    tag "$meta.id"
    label 'process_single'

    container 'quay.io/biocontainers/genometools-genometools:1.6.5--py310h3db02ab_0'

    input:
    tuple val(meta), path(map_gff)

    script:
    """
    gt gff3validator ${map_gff}
    """
}

