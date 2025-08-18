process PROCESS_BLASTP_PROKKA {
    tag "${meta.id}"
    label 'process_single'

    container 'quay.io/microbiome-informatics/pandas-biopython:latest'

    input:
    tuple val(meta), path(uniprot_tsv)

    output:
    tuple val(meta), path("*_uniprot_names.tsv"), emit: uniprot_product_names
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    prokka_compatible_blastp_processor.py \\
        --input ${uniprot_tsv} \\
        --output ${prefix}_uniprot_names.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}
