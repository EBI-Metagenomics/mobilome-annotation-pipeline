process VMATCH_PROCESS {
    tag "$meta.id"
    label 'process_single'
   
    container 'quay.io/microbiome-informatics/pandas-biopython:latest'

    input:
    tuple val(meta), path(vmatch_result)

    output:
    tuple val(meta), path("*_direct_repeats.tsv"), emit: dr_tsv

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    vmatch_processor.py \\
        ${vmatch_result} \\
        ${prefix}_direct_repeats.tsv
    """
}
