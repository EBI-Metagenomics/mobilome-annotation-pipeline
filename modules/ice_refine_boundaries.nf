process REFINE_BOUNDARIES {
    tag "${meta.id}"
    label 'process_single'

    container 'quay.io/biocontainers/biopython:1.81'
    
    input:
    tuple val(meta), path(assembly), path(merged_gff), path(macsyfinder_tsv), path(uniprot_product_names), path(vmatch_tsv)
    
    output:
    tuple val(meta), path("*_ices.tsv")     , emit: ices_tsv
    tuple val(meta), path("*_ice_genes.tsv"), emit: ice_genes_tsv
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ice_boundary_refinement.py \\
        --assembly ${assembly} \\
        --gff_file ${merged_gff} \\
        --macsyfinder_out ${macsyfinder_tsv} \\
        --uniprot_annot ${uniprot_product_names} \\
        --drs_tsv ${vmatch_tsv} \\
        --prefix ${prefix}
    """
}
