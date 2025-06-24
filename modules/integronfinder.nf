process INTEGRONFINDER {
    tag "$meta.id"
    label 'process_high'

    container 'quay.io/microbiome-informatics/integronfinder:71ee6e0'

    input:
    tuple val(meta), path(assembly_file)

    output:
    tuple val(meta), path("*.summary"), emit: contigs_summary
    tuple val(meta), path("*.gbk"),     emit: contigs_gbks

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    integron_finder --union-integrases \\
        --local-max \\
        --cpu ${task.cpus} \\
        --gbk \\
        ${assembly_file}

    if ls -l Results_Integron_Finder_${prefix}_5kb_contigs/*.gbk 2>/dev/null | grep -q .
    then
        echo 'IntegronFinder outputs complete'
        mv Results_Integron_Finder_${prefix}_5kb_contigs/*.gbk .
        mv Results_Integron_Finder_${prefix}_5kb_contigs/*.summary .
    else
        echo 'IntegronFinder found 0 integrons in assembly... generating dummy files'
        touch contig_dummy.gbk
        touch contig_dummy.summary
    fi
    """
}
