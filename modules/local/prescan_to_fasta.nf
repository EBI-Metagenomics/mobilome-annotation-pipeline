process PRESCAN_TO_FASTA {
    tag "${meta.id}"
    label 'process_single'

    container 'quay.io/biocontainers/biopython:1.81'

    input:
    tuple val(meta), path(hmmscan_tbl), path(faa), path(assembly), path(gff)
    val evalue_threshold

    output:
    tuple val(meta), path("*_candidates.fasta"), emit: candidates_fna, optional: true
    tuple val(meta), path("*_candidates.faa"), emit: candidates_faa, optional: true
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    prescan_to_fasta.py \\
        --hmm_out ${hmmscan_tbl} \\
        --proteins ${faa} \\
        --assembly ${assembly} \\
        --gff_file ${gff} \\
        --evalue_threshold ${evalue_threshold} \\
        --output ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
}
