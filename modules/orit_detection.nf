process ORIT_DETECTION {
    tag "$meta.id"
    label 'process_medium'
    
    conda "bioconda::blast=2.14.1 bioconda::python=3.9 bioconda::biopython=1.79"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.14.1--pl5321h6f7f691_0':
        'biocontainers/blast:2.14.1--pl5321h6f7f691_0' }"
    
    input:
    tuple val(meta), path(ice_sequences)    // ICE sequence files from validated ICE elements
    path orit_database                      // oriT reference database (BLAST format)
    
    output:
    tuple val(meta), path("${prefix}_orit_results.tsv"), emit: orit_results
    tuple val(meta), path("${prefix}_orit_sequences.fasta"), emit: orit_sequences, optional: true
    tuple val(meta), path("${prefix}_blast_output.txt"), emit: blast_output
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def evalue = task.ext.evalue ?: '0.01'
    def word_size = task.ext.word_size ?: '11'
    def identity_threshold = task.ext.identity_threshold ?: '49.0'
    def coverage_threshold = task.ext.coverage_threshold ?: '0.49'
    
    """
    detect_orit.py \\
        ${ice_sequences} \\
        ${orit_database} \\
        ${prefix} \\
        --evalue ${evalue} \\
        --word-size ${word_size} \\
        --identity-threshold ${identity_threshold} \\
        --coverage-threshold ${coverage_threshold} \\
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>&1 | sed 's/^.*blastn: //; s/ .*\$//')
        python: \$(python --version | sed 's/Python //g')
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
    
    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_orit_results.tsv
    touch ${prefix}_orit_sequences.fasta
    touch ${prefix}_blast_output.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>&1 | sed 's/^.*blastn: //; s/ .*\$//' || echo "2.14.1")
        python: \$(python --version | sed 's/Python //g' || echo "3.9.0")
        biopython: \$(echo "1.79")
        pandas: \$(echo "1.3.0")
    END_VERSIONS
    """
}
