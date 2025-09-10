process MACSYFINDER_TO_FASTA {
    tag "$meta.id"
    label 'process_single'

    container 'quay.io/microbiome-informatics/pandas-biopython:latest'

    input:
    tuple val(meta), path(macsyfinder_tsv), path(assembly_file), path(gff_file)

    output:
    tuple val(meta), path("*_contigs.fasta"), emit: macsy_contigs
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    macsy_to_fasta.py \\
        --macsy_input ${macsyfinder_tsv} \\
        --assembly ${assembly_file} \\
        --gff ${gff_file} \\
        --outout ${prefix}_contigs.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
}
