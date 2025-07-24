process HMMSCAN {
    tag "$meta.id"
    label 'process_medium'

    container 'quay.io/biocontainers/hmmer:3.0--h503566f_5'

    input:
    tuple val(meta), path(faa_file)
    tuple val(meta2), path(ice_hmm_models)

    output:
    tuple val(meta), path("*_prescan.tbl"), emit: hmmscan_tbl

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    hmmscan \\
        --tblout ${prefix}_prescan.tbl \\
        ${meta2.id} \\
        ${faa_file}
    """
}
