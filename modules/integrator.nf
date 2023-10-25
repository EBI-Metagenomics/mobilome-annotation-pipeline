#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process INTEGRATOR {
    publishDir "$launchDir/${params.outdir}/", mode: 'copy'
    stageOutMode = 'copy'

    container 'quay.io/microbiome-informatics/virify-python3:1.2'

    input:
	path prokka_gff
	path map_file
	path iss_tsv
	path pal_info
	path inf_summ
	path inf_gbk
	path icf_summ
	path icf_dr
	path mog_out
	path genomad_vir
	path genomad_plas
	path vir_results
	path crispr_tsv

    output:
	path 'mobilome_prokka.gff', emit: mobilome_prokka_gff
	path 'overlapping_integrons.txt'
	path 'discarded_mge.txt'

    script:
    	"""    
    	mge_integrator.py \
    	--pkka_gff ${prokka_gff} \
    	--map ${map_file} \
    	--iss_tsv ${iss_tsv} \
    	--pal_tsv ${pal_info} \
    	--inf_tsv ${inf_summ} \
    	--inf_gbks ${inf_gbk.join(' ')} \
    	--icf_tsv ${icf_summ} \
    	--icf_lim ${icf_dr} \
    	--mog_tsv ${mog_out} \
	--geno_out ${genomad_vir} \
	--geno_plas ${genomad_plas} \
        --virify_out ${vir_results} \
	--crispr_out ${crispr_tsv}
	"""
}

