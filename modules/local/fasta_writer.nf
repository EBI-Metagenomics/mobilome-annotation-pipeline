process FASTA_WRITER {
    // TODO: replace with a modified version of https://nf-co.re/modules/bedtools_getfasta/
    tag "${meta.id}"
    label 'process_single'

    container 'quay.io/biocontainers/bedtools:2.23.0--h5b5514e_6'

    input:
    tuple val(meta), path(assembly), path(mobilome_nogenes)

    output:
    path "*_mobilome.fasta.gz"
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Assembly could be compressed (.gz)
    def is_assembly_compressed = assembly.name.endsWith(".gz")
    def assembly_file = assembly ? assembly.name.replace(".gz", "") : is_assembly_compressed
    """
    # GFF file is always compressed
    gzip -dc ${mobilome_nogenes} > ${prefix}_nogenes.gff

    awk -F '\t' '{ if (NF == 9 ) print \$1 "\t" \$4 "\t" \$5 "\t" \$9}' ${prefix}_nogenes.gff | cut -d';' -f1 > mobilome_nogenes.bed

    # Handle assembly file (compressed or uncompressed)
    if [ "${is_assembly_compressed}" == "true" ]; then
         gzip -c -d ${assembly} > ${assembly_file}
    fi

    bedtools getfasta -fi ${assembly_file} -bed mobilome_nogenes.bed -name -fo ${prefix}_mobilome.fasta

    sed -i 's/ID=//;s/::.*//' ${prefix}_mobilome.fasta

    gzip ${prefix}_mobilome.fasta


    # Clean up temporary files
    rm ${prefix}_nogenes.gff
    if [ "${is_assembly_compressed}" == "true" ]; then
         rm -f ${assembly_file}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | head -1 | sed 's/bedtools v//g')
    END_VERSIONS
    """
}
