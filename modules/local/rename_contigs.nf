process RENAME {
    tag "${meta.id}"
    label 'process_single'

    container 'quay.io/biocontainers/biopython:1.81'

    input:
    tuple val(meta), path(assembly_file)

    output:
    tuple val(meta), path('*_1kb_contigs.fasta'),   emit: contigs_1kb
    tuple val(meta), path('*_5kb_contigs.fasta'),   emit: contigs_5kb
    tuple val(meta), path('*_100kb_contigs.fasta'), emit: contigs_100kb
    tuple val(meta), path('*_contigID.map'),        emit: map_file

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    assembly_filter_rename.py \\
        --assembly ${assembly_file} \\
        --output ${prefix}
    """
}
