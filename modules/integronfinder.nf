process INTEGRONFINDER {
    tag "$meta.id"
    label 'process_high'

    container 'quay.io/microbiome-informatics/integronfinder:71ee6e0'

    input:
    tuple val(meta), path(assembly_file)

    output:
    tuple val(meta), path("Results_Integron_Finder_5kb_contigs/5kb_contigs.summary"), emit: contigs_summary
    tuple val(meta), path("Results_Integron_Finder_5kb_contigs/contig_*.gbk")       , emit: contigs_gbks

    script:
    """
    integron_finder --union-integrases \
    --local-max \
    --cpu ${task.cpus} \
    --gbk \
    ${assembly_file}

    if ls -l Results_Integron_Finder_5kb_contigs/contig_*.gbk 2>/dev/null | grep -q .
    then
        echo 'IntegronFinder outputs complete'
    else
        echo 'IntegronFinder found 0 integrons in assembly... generating dummy files'
        touch Results_Integron_Finder_5kb_contigs/contig_dummy.gbk
    fi
    """
}
