#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process integra {
    publishDir "$launchDir/MoMofy_results/", mode: 'copy'

    container "quay.io/biocontainers/biopython:1.75"

    memory "4 GB"
    cpus 1

    input:
	path cgc_fasta
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
	//path 'test.out'
	path 'momofy_predictions.fna'
	path 'momofy_predictions.gff'
	path 'nested_integrons.txt'
	path 'discarded_iss.txt'


    """    
    mge_integrator.py --cgc_fa ${cgc_fasta} \
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
    """
}

//workflow {
//    integra(CGC_COOR, MAP_FILE, ISS_SEQS, ISS_PRED, PAL_SEQS, PAL_PRED, INF_PRED, INF_GBK, ICF_PRED, ICF_SEQS, MOG_ALN)
//}

