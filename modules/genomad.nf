process GENOMAD {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/genomad:1.6.1--pyhdfd78af_0':
        'biocontainers/genomad:1.6.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(assembly_file)

    output:
    tuple val(meta), path("genomad_out/5kb_contigs_summary/5kb_contigs_virus_summary.tsv"), emit: genomad_vir
    tuple val(meta), path("genomad_out/5kb_contigs_summary/5kb_contigs_plasmid_summary.tsv"), emit: genomad_plas

    script:
    """    
    genomad end-to-end ${assembly_file} \\
    --threads ${task.cpus} \\
    genomad_out ${params.genomad_db}
    """
}
