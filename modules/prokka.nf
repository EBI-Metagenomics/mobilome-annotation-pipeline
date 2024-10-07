#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process PROKKA {

    publishDir "$params.outdir/preprocessing", mode: 'copy'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/prokka:1.14.6--pl526_0' :
        'biocontainers/prokka:1.14.6--pl526_0' }"

    input:
    path assembly_file

    output:
    path 'prokka_out/contigs.gbk', emit: prokka_gbk
	path 'prokka_out/contigs.gff', emit: prokka_gff
	path 'prokka_out/contigs.faa', emit: prokka_faa
	path 'prokka_out/contigs.fna', emit: prokka_fna

    script:
    if(assembly_file.size() > 0)
        """
        # TMP folder issues in Prokka - https://github.com/tseemann/prokka/issues/402
        export TMPDIR="\$PWD/tmp"
        mkdir -p "\$PWD/tmp"
        # Disable the Java VM performane gathering tool, for improved performance
        export JAVA_TOOL_OPTIONS="-XX:-UsePerfData"

        prokka --outdir prokka_out \
        --prefix contigs \
        --cpus ${task.cpus} \
        --metagenome \
        ${assembly_file}
        """
    else
        """
        echo 'PROKKA dir empty due to empty input... generating dummy files'
        mkdir prokka_out
        touch prokka_out/contigs.gbk
        """
}


