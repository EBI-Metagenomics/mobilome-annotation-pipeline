#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process diamond_mob {
    publishDir "$launchDir/mobileog_results"

    memory "8 GB"
    cpus 4

    container "quay.io/biocontainers/diamond:2.0.12--hdcc8f71_0"

    input:
        path proteins_file, name: 'proteins.faa'
	path diamond_db 

    output:
        path 'blastp_out.tsv', emit:blast_out

    script:
    if(proteins_file.size() > 0)
        """
        diamond blastp -q proteins.faa \
        --db ${diamond_db} \
        --outfmt 6 stitle qtitle pident bitscore slen evalue qlen sstart send qstart qend \
        -k 15 \
        -o blastp_out.tsv \
        -e 1e-20 \
        --query-cover 90 \
        --id 90 \
        --threads ${task.cpus}
        """
    else
        """
	echo 'No input files for diamond... generating dummy files'
        touch blastp_out.tsv
        """
}

