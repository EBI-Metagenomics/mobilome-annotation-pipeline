process ARAGORN {
    // TODO: Push to nf-core
    tag "${meta.id}"
    label 'process_single'

    container 'quay.io/biocontainers/aragorn:1.2.41--h7b50bb2_5'

    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path("*_aragorn.tbl"), emit: rnas_tbl
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    aragorn \\
        ${args} \\
        -o ${prefix}_aragorn.tbl \\
        ${contigs}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        aragorn: \$(aragorn -h | sed '2q;d' | sed 's/ARAGORN v//' | cut -d ' ' -f1)
    END_VERSIONS
    """
}
