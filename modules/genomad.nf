process GENOMAD {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/genomad:1.6.1--pyhdfd78af_0':
        'biocontainers/genomad:1.6.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(assembly_file)

    output:
    tuple val(meta), path("*_5kb_contigs_summary/*_5kb_contigs_virus_summary.tsv"),   emit: genomad_vir
    tuple val(meta), path("*_5kb_contigs_summary/*_5kb_contigs_plasmid_summary.tsv"), emit: genomad_plas

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """    
    genomad end-to-end \\
        --threads ${task.cpus} \\
        ${assembly_file} \\
        . \\
        ${params.genomad_db}
    """
}
