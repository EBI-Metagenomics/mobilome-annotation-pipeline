process PROCESS_BLASTP_PROKKA {
    tag "$meta.id"
    label 'process_single'
   
    container 'quay.io/microbiome-informatics/pandas-biopython:latest'


    input:
    tuple val(meta), path(uniprot_tsv)

    output:
    tuple val(meta), path("*_uniprot_names.tsv"), emit: uniprot_product_names

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    prokka_blastp_processor.py \\
        --input ${uniprot_tsv} \\
        --output ${prefix}_uniprot_names.tsv
    """
}
