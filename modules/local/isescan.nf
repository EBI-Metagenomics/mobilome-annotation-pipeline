process ISESCAN {
    tag "${meta.id}"
    label 'process_high'

    container 'quay.io/biocontainers/isescan:1.7.3--h7b50bb2_0'

    input:
    tuple val(meta), path(assembly_file)

    output:
    tuple val(meta), path('*.fasta.tsv'), emit: iss_tsv
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    isescan.py \\
        --seqfile ${assembly_file} \\
        --output . \\
        --nthread ${task.cpus}

    if ls -l ${prefix}_1kb_contigs.fasta.tsv 2>/dev/null | grep -q .
    then
        echo 'ISEScan results exists'
    else
    echo 'ISEScan found 0 insertion sequences in input file... generating dummy files'
        touch ${prefix}_1kb_contigs.fasta.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isescan: \$(isescan.py --version | sed 's/isescan //g')
    END_VERSIONS
    """
}
