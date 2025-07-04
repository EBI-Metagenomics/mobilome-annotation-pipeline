process GENOMAD {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/genomad:1.11.1--pyhdfd78af_0':
        'biocontainers/genomad:1.11.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(assembly_file)

    output:
    tuple val(meta), path("*_5kb_contigs_summary/*_5kb_contigs_virus_summary.tsv"),   emit: genomad_vir
    tuple val(meta), path("*_5kb_contigs_summary/*_5kb_contigs_plasmid_summary.tsv"), emit: genomad_plas

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """    
    if [ -s ${assembly_file} ]; then
        genomad end-to-end ${assembly_file} \\
            --threads ${task.cpus} \\
            . ${params.genomad_db}
    else
        mkdir -p 5kb_contigs_summary
        touch ${prefix}_5kb_contigs_summary/${prefix}_5kb_contigs_virus_summary.tsv
        touch ${prefix}_5kb_contigs_summary/${prefix}_5kb_contigs_plasmid_summary.tsv
    fi
    """
}
