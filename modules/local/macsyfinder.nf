process MACSYFINDER {
    tag "${meta.id}"
    label 'process_medium'

    container "quay.io/biocontainers/macsyfinder:2.1.4--pyhdfd78af_1"

    input:
    tuple val(meta), path(faa_file)
    path ice_models

    output:
    tuple val(meta), path("*_macsyfinder_results/all_systems.tsv"), emit: macsyfinder_tsv
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    macsyfinder \\
        --sequence-db ${faa_file} \\
        --models-dir ${ice_models} \\
        --out-dir ${prefix}_macsyfinder_results \\
        --worker ${task.cpus} ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        macsyfinder: \$(macsyfinder --version | head -n 1 | sed 's/Macsyfinder //g')
    END_VERSIONS
    """
}
