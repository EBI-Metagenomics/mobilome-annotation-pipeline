#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process AMRFINDER_PLUS {

    publishDir "$params.outdir/func_annot", mode: 'copy'

    container 'quay.io/biocontainers/ncbi-amrfinderplus:3.11.4--h6e70893_0'

    input:
        path fna
        path faa
        path gff

    output:
        path("amrfinderplus.tsv"), emit: amrfinder_tsv


    script:
    if (fna.size() > 0)
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
}

