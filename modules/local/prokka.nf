process PROKKA {
    tag "$meta.id"
    label 'process_high'

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/prokka:1.14.6--pl526_0'
        : 'biocontainers/prokka:1.14.6--pl526_0'}"

    input:
    tuple val(meta), path(assembly_file)

    output:
    tuple val(meta), path("prokka_results/*.gbk"), emit: prokka_gbk
    tuple val(meta), path("prokka_results/*.gff"), emit: prokka_gff
    tuple val(meta), path("prokka_results/*.faa"), emit: prokka_faa
    tuple val(meta), path("prokka_results/*.fna"), emit: prokka_fna

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # TMP folder issues in Prokka - https://github.com/tseemann/prokka/issues/402
    export TMPDIR="\$PWD/tmp"
    mkdir -p "\$PWD/tmp"
    # Disable the Java VM performance gathering tool, for improved performance
    export JAVA_TOOL_OPTIONS="-XX:-UsePerfData"

    prokka \\
        --outdir prokka_results \\
        --prefix ${prefix} \\
        --cpus ${task.cpus} \\
        --metagenome \\
        ${assembly_file}
    """
}
