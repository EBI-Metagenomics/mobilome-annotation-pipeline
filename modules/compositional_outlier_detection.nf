process OUTLIER_FINDER {
    tag "${fasta.baseName}"
    label 'process_medium'

    container 'quay.io/microbiome-informatics/python_bio_numpy:v3.9.23'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.bed"), emit: bed

    script:
    def outliers_args = task.ext.outliers_args
    def prefix = task.ext.prefix ?: "${fasta.baseName}"
    """
    fast_composition_analyzer.py \\
        --input_fasta ${fasta} \\
        --output_bed ${prefix}.bed \\
        --threads ${task.cpus} $outliers_args
    """
}
