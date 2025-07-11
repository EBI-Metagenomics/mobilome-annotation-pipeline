process ARAGORN {
    tag "$meta.id"
    label 'process_single'

    container 'quay.io/biocontainers/aragorn:1.2.41--h7b50bb2_5'

    input:
    tuple val(meta), path(fasta_file)

    output:
    tuple val(meta), path("*_aragorn.gff"), emit: trna_gff

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def aragorn_args = task.ext.aragorn_args ?: ''
    """
    aragorn \\
        ${aragorn_args} \\
        -o ${prefix}_aragorn.tbl \\
        ${fasta_file}
    
    # Convert into GFF format 
    aragorn_to_gff.sh ${prefix}_aragorn.tbl ${prefix}_aragorn.gff
    """
}
