#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process VIRI_PARSER {
    publishDir "$launchDir/other_mge_results"

    container 'quay.io/microbiome-informatics/virify-python3:1.2'

    input:
        value vir_gff_path
        path assembly_file

    output:
        path("all_virify.fasta"), emit: viral_fasta
        path("plasmids.list"), emit: plasmids_list

    script:
    """   
    virify_parser.py \
    --virify_gff ${vir_gff_path} \
    --original_assem ${assembly_file}
    """
}

