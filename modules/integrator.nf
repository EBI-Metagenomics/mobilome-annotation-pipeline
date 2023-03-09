#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process integra {
    publishDir "$launchDir/MoMofy_results/", mode: 'copy'

    container "quay.io/microbiome-informatics/virify-python3:1.2"

    memory "4 GB"
    cpus 1

    input:
	path assembly
	path cds_gff
	path mapping_file
	path iss_fasta
	path iss_table
	path pal_fasta
	path pal_table
	path inf_table
	path inf_list
	path icf_table
	path icf_fasta
	path mog_table

    output:
	path 'test.out'
	path 'momofy_predictions.fna'
	path 'momofy_predictions.gff'
	path 'nested_integrons.txt'
	path 'discarded_mge.txt'
	path 'validated_momofy_predictions.gff'

    """    
    mge_integrator.py --assem ${assembly} \
    --cds_gff ${cds_gff} \
    --map ${mapping_file} \
    --iss_fa ${iss_fasta} \
    --iss_tsv ${iss_table} \
    --pal_fa ${pal_fasta} \
    --pal_tsv ${pal_table} \
    --inf_tsv ${inf_table} \
    --inf_gbks ${inf_list.join(' ')} \
    --icf_tsv ${icf_table} \
    --icf_fa ${icf_fasta} \
    --mog_tsv ${mog_table}

    gt gff3 \
    -retainids yes \
    -checkids yes \
    -o validated_momofy_predictions.gff \
    momofy_predictions.gff
    """
}

