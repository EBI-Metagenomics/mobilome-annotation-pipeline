process FASTA_WRITER {
    tag "$meta.id"
    label 'process_single'

    container 'quay.io/biocontainers/bedtools:2.23.0--h5b5514e_6'

    input:
    tuple val(meta), path(assembly), path(mobilome_nogenes)

    output:
    path "${meta.id}_mobilome.fasta"

    script:
    """
    awk -F '\t' '{ if (NF == 9 ) print \$1 "\t" \$4 "\t" \$5 "\t" \$9}' ${mobilome_nogenes} | cut -d';' -f1 > mobilome_nogenes.bed

    bedtools getfasta -fi ${assembly} -bed mobilome_nogenes.bed -name -fo ${meta.id}_mobilome.fasta

    sed -i 's/ID=//;s/::.*//' ${meta.id}_mobilome.fasta
    """
}
