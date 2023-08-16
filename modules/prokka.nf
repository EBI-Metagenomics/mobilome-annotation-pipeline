#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process PROKKA {
    publishDir "$launchDir/preprocessing"

    container 'quay.io/biocontainers/prokka:1.14.6--pl526_0'

    input:
        path assembly_file, name: 'contigs.fasta'

    output:
        path 'prokka_out/contigs.gbk', emit: prokka_gbk
	path 'prokka_out/contigs.gff', emit: prokka_gff
	path 'prokka_out/contigs.faa', emit: prokka_faa
	path 'prokka_out/contigs.fna', emit: prokka_fna

    script:
    if(assembly_file.size() > 0)
        """
        prokka --outdir prokka_out \
        --prefix contigs \
        --cpus ${task.cpus} \
        --metagenome \
        contigs.fasta
        """
    else
        """
        echo 'PROKKA dir empty due to empty input... generating dummy files'
        mkdir prokka_out
        touch prokka_out/contigs.gbk
        """
}


