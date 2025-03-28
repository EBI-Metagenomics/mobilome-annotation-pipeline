process PROKKA {
    tag "$meta.id"
    label 'process_high'

   

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/prokka:1.14.6--pl526_0'
        : 'biocontainers/prokka:1.14.6--pl526_0'}"

    input:
    tuple val(meta), path(assembly_file)

    output:
    tuple val(meta), path("prokka_out/${meta.id}.gbk"), emit: prokka_gbk
    tuple val(meta), path("prokka_out/${meta.id}.gff"), emit: prokka_gff
    tuple val(meta), path("prokka_out/${meta.id}.faa"), emit: prokka_faa
    tuple val(meta), path("prokka_out/${meta.id}.fna"), emit: prokka_fna

    script:
    """
    # TMP folder issues in Prokka - https://github.com/tseemann/prokka/issues/402
    export TMPDIR="\$PWD/tmp"
    mkdir -p "\$PWD/tmp"
    # Disable the Java VM performane gathering tool, for improved performance
    export JAVA_TOOL_OPTIONS="-XX:-UsePerfData"

    prokka \\
        --outdir prokka_out \\
        --prefix ${meta.id} \\
        --cpus ${task.cpus} \\
        --metagenome \\
        ${assembly_file}
    """
}
