process SYMBIOSIS_BLAST {
    tag "$meta.id"
    label 'process_medium'
    
    conda "bioconda::blast=2.14.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.14.1--h7d5a4b4_1':
        'biocontainers/blast:2.14.1--h7d5a4b4_1' }"

    input:
    tuple val(meta), path(fasta)
    path(symbiosis_db)

    output:
    tuple val(meta), path("*.sym.m8"), emit: blast_results
    tuple val(meta), path("*.sym.filtered"), emit: filtered_results
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    blastn \\
        -query ${fasta} \\
        -db ${symbiosis_db}/symbiosis \\
        -evalue 0.0001 \\
        -num_threads ${task.cpus} \\
        -max_hsps 1 \\
        -num_descriptions 1 \\
        -num_alignments 1 \\
        -outfmt "6 std slen stitle" \\
        -out ${prefix}.sym.m8 \\
        ${args}

    # Filter results with H-value >= 0.64
    filter_blast_results.py \\
        --input ${prefix}.sym.m8 \\
        --output ${prefix}.sym.filtered \\
        --threshold 0.64 \\
        --type nucleotide

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastp -version 2>&1 | sed 's/^.*blastp: //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sym.m8
    touch ${prefix}.sym.filtered
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastp -version 2>&1 | sed 's/^.*blastp: //; s/ .*\$//')
    END_VERSIONS
    """
}
