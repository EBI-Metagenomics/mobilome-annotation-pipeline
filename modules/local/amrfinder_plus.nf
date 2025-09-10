process AMRFINDER_PLUS {
    // FIXME: replace with the nf-core module
    tag "${meta.id}"
    label 'process_high'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/ncbi-amrfinderplus:3.12.8--h283d18e_0'
        : 'biocontainers/ncbi-amrfinderplus:3.12.8--h283d18e_0'}"

    input:
    tuple val(meta), path(fna), path(faa), path(gff)

    output:
    tuple val(meta), path("*_amrfinderplus.tsv"), emit: amrfinder_tsv
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """    
    amrfinder --plus \\
        -n ${fna} \\
        -p ${faa} \\
        -g ${gff} \\
        -d ${params.amrfinder_plus_db} \\
        -a prokka \\
        --output ${prefix}_amrfinderplus.tsv \\
        --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amrfinder: \$(amrfinder --version)
        amrfinderplus-database: \$(echo \$(amrfinder --database amrfinderdb --database_version 2> stdout) | rev | cut -f 1 -d ' ' | rev)
    END_VERSIONS
    """
}
