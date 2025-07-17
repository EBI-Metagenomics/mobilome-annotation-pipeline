process VMATCH {
    tag "$meta.id"
    label 'process_single'
    
    container 'quay.io/biocontainers/vmatch:2.3.0--h516909a_0'
    
    input:
    tuple val(meta), path(macsy_boundaries), path(assembly)
    
    output:
    tuple val(meta), path("*_vmatch.tsv"), emit: vmatch_tsv
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def index_args = task.ext.index_args ?: ''
    def vmatch_args = task.ext.vmatch_args ?: ''
    """
    # Extract contigs of interest only
    macsy_to_fasta.sh \\
        ${macsy_boundaries} \\
        ${assembly} > ${prefix}_contigs.fasta

    # Create vmatch index
    mkvtree \\
        -db ${prefix}_contigs.fasta \\
        ${index_args} \\
        -indexname ${prefix}_vindex
    
    # Find direct repeats using vmatch
    vmatch \\
        ${vmatch_args} \\
        ${prefix}_vindex > ${prefix}_vmatch.out 

    # Transform into coordinates table
    vmatch_to_tsv.sh ${prefix}_vmatch.out ${prefix}_vmatch.tsv
    """
}
