#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process AMRFINDER_PLUS {
    publishDir "$params.outdir/func_annot"

    cpus 4
    memory { 8.GB * task.attempt }
    errorStrategy 'retry' 
    maxRetries 3

    container 'quay.io/biocontainers/ncbi_amrfinderplus:3.11.4'

    input:
        path fna
        path faa
        path gff

    output:
        path("amrfinderplus.tsv"), emit: amrfinder_tsv

    script:
    if(fna.size() > 0)
        """    
        amrfinder --plus \
        -n ${fna} \
        -p ${faa} \
        -g ${gff} \
        -d ${params.amrfinder_plus_db} \
        -a prokka \
        --output amrfinderplus.tsv \
        --threads ${task.cpus}
        """
    else
        """
        echo 'AMRFinderPlus output file empty due to empty input... generating dummy files'
        touch amrfinderplus.tsv
        """
}

