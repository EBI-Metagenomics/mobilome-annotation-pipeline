process RENAME {

    tag "${meta.id}"
    label 'process_single'

    container 'quay.io/biocontainers/biopython:1.75'

    input:
    tuple val(meta), path(assembly_file)

    output:
    tuple val(meta), path('1kb_contigs.fasta'), emit: contigs_1kb
    tuple val(meta), path('5kb_contigs.fasta'), emit: contigs_5kb
    tuple val(meta), path('contigID.map')     , emit: map_file

    script:
    """
    assembly_filter_rename.py --assembly ${assembly_file}
    """
}
