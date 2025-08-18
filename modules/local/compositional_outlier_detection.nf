process OUTLIER_FINDER {
    tag "${fasta.baseName}"
    label 'process_medium'

    container 'quay.io/microbiome-informatics/python_bio_numpy:v3.9.23'

    input:
    tuple val(meta), path(fasta)
    val score_threshold

    output:
    tuple val(meta), path("*.bed"), emit: bed
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${fasta.baseName}"
    """
    fast_composition_analyzer.py \\
        --input_fasta ${fasta} \\
        --output_bed ${prefix}.bed \\
        --score-threshold ${score_threshold} \\
        --threads ${task.cpus} ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        numpy: \$(python -c "import numpy; print(numpy.__version__)")
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
}
