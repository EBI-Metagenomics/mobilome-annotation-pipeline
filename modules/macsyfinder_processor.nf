process MACSYFINDER_PROCESS {
    tag "$meta.id"
    label 'process_single'
   
    container 'quay.io/microbiome-informatics/pandas-biopython:latest'


    input:
    tuple val(meta), path(macsy_result), path(assembly_file), path(gff_file)

    output:
    tuple val(meta), path("*_ices_result/boundary_summary.json")        , emit: summary_json
    tuple val(meta), path("*_ices_result/boundaries.tsv")               , emit: boundaries_tsv
    tuple val(meta), path("*_ices_result/*_boundaries.json")            , emit: per_sys_json, optional: true
    tuple val(meta), path("*_ices_result/all_systems.fasta")            , emit: all_sys_fasta
    tuple val(meta), path("*_ices_result/all_systems_with_flanks.fasta"), emit: all_sys_flanks_fasta

    script:
    def macsy_process_args = task.ext.macsy_process_args
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    macsyfinder_processor.py \\
        --all-systems ${macsy_result} \\
        --genome ${assembly_file} \\
        --annotation ${gff_file} \\
        --output-dir ${prefix}_ices_result \\
        ${macsy_process_args}
    """
}
