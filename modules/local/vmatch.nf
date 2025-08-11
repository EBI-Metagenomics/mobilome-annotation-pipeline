process VMATCH {
    tag "$meta.id"
    label 'process_single'
    
    container 'quay.io/biocontainers/vmatch:2.3.0--h516909a_0'
    
    input:
    tuple val(meta), path(contigs)
    
    output:
    tuple val(meta), path("*_vmatch.tsv"), emit: vmatch_tsv
    path "versions.yml"                 , emit: versions
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def index_args = task.ext.index_args ?: ''
    def vmatch_args = task.ext.vmatch_args ?: ''
    """
    # Create vmatch index
    mkvtree \\
        -db ${contigs} \\
        ${index_args} \\
        -indexname ${prefix}_vindex
    
    # Find direct repeats using vmatch
    vmatch \\
        ${vmatch_args} \\
        ${prefix}_vindex > ${prefix}_vmatch.out 

    # Transform into coordinates table
    vmatch_to_tsv.sh ${prefix}_vmatch.out ${prefix}_vmatch.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vmatch: \$(vmatch -version 2>&1 | grep -o '[0-9]\+\.[0-9]\+\.[0-9]\+')
    END_VERSIONS
    """
}
