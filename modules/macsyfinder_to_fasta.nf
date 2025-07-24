process MACSYFINDER_TO_FASTA {
    tag "$meta.id"
    label 'process_single'
   

    container 'quay.io/microbiome-informatics/pandas-biopython:latest'

    input:
    tuple val(meta), path(macsyfinder_tsv), path(assembly_file), path(gff_file)

    output:
    tuple val(meta), path("*_contigs.fasta"), emit: macsy_contigs

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    macsy_to_fasta.py \\
        --macsy_input ${macsy_boundaries} \\
        --assembly ${assembly_file} \\
        --gff ${gff_file} \\
        --outout ${prefix}_contigs.fasta
    """
}
