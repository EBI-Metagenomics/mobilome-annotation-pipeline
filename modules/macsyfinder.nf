process MACSYFINDER {
    tag "$meta.id"
    label 'process_medium'

    container "quay.io/biocontainers/macsyfinder:2.1.4--pyhdfd78af_1"
    
    input:
    tuple val(meta), path(faa_file)
    path ice_models

    output:
    tuple val(meta), path("${meta.id}_macsyfinder_results/all_systems.tsv"), emit: macsyfinder_tsv

    script:
    def macsyfinder_args = task.ext.macsyfinder_args ?: ''
    """
    # Run MacSyFinder for ICE detection
    macsyfinder \\
        --sequence-db ${faa_file} \\
        --models-dir ${ice_models} \\
        --out-dir ${meta.id}_macsyfinder_results \\
        --worker $task.cpus \\
        ${macsyfinder_args}
    """
}
