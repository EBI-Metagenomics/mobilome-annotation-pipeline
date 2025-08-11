process PRODIGAL {
    tag "$meta.id"
    label 'process_single'

    container 'quay.io/biocontainers/prodigal:2.6.3--hec16e2b_4'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.gff")  , emit: annot
    tuple val(meta), path("*.faa")  , emit: faa

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def prodigal_args = task.ext.prodigal_args ?: ''
    """
    prodigal \\
        -i ${fasta} \\
        -o ${prefix}.gff \\
        -a ${prefix}.faa \\
        ${prodigal_args}
    """
}
