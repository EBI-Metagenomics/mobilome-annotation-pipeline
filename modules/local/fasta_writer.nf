process FASTA_WRITER {
    // TODO: replace with a modified version of https://nf-co.re/modules/bedtools_getfasta/
    tag "${meta.id}"
    label 'process_single'

    container 'quay.io/biocontainers/bedtools:2.23.0--h5b5514e_6'

    input:
    tuple val(meta), path(assembly), path(mobilome_nogenes)

    output:
    path "*_mobilome.fasta"
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    awk -F '\t' '{ if (NF == 9 ) print \$1 "\t" \$4 "\t" \$5 "\t" \$9}' ${mobilome_nogenes} | cut -d';' -f1 > mobilome_nogenes.bed

    bedtools getfasta -fi ${assembly} -bed mobilome_nogenes.bed -name -fo ${prefix}_mobilome.fasta

    sed -i 's/ID=//;s/::.*//' ${prefix}_mobilome.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | head -1 | sed 's/bedtools v//g')
    END_VERSIONS
    """
}
