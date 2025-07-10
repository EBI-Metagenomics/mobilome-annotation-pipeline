process VMATCH {
    tag "$meta.id"
    label 'process_single'
    
    container 'quay.io/biocontainers/vmatch:2.3.0--h516909a_0'
    
    input:
    tuple val(meta), path(extended_regions)
    
    output:
    tuple val(meta), path("*_vmatch.out"), emit: vmatch_raw
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def index_args = task.ext.index_args ?: ''
    def vmatch_args = task.ext.vmatch_args ?: ''
    """
    # Create vmatch index
    mkvtree \\
        -db ${extended_regions} \\
        ${index_args} \\
        -indexname ${prefix}_vindex
    
    # Find direct repeats using vmatch
    vmatch \\
        ${vmatch_args} \\
        ${prefix}_vindex > ${prefix}_vmatch.out 
    """
}
