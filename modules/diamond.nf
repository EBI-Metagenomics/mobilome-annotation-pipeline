#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process DIAMOND {
    publishDir "$launchDir/$params.outdir/func_annot"

    cpus 4
    memory { 8.GB * task.attempt }
    errorStrategy 'retry' 
    maxRetries 3

    container 'quay.io/biocontainers/diamond:2.0.12--hdcc8f71_0'

    input:
        path proteins_file
	path diamond_db 

    output:
        path 'mobileOG.tsv', emit:blast_out

    script:
    if(proteins_file.size() > 0)
        """
        diamond blastp \
        -q ${proteins_file} \
        --db ${diamond_db} \
        --outfmt 6 stitle qtitle pident bitscore slen evalue qlen sstart send qstart qend \
        -k 15 \
        -o mobileOG.tsv \
        -e 1e-20 \
        --query-cover 90 \
        --id 90 \
        --threads ${task.cpus}
        """
    else
        """
	echo 'No input files for diamond... generating dummy files'
        touch mobileOG.tsv
        """
}

