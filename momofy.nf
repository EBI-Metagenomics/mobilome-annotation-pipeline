#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// MODULES
include { diamond_mob } from './modules/run_diamond'
include { rename } from './modules/rename_contigs'
include { prokka_annot } from './modules/run_prokka'
include { gbk_split } from './modules/gbk_splitter'
include { integronfinder } from './modules/run_integronfinder'
include { isescan } from './modules/run_isescan'
include { ice_finder } from './modules/run_icefinder'
include { integra } from './modules/integrator'

// INPUTS
assembly = Channel.fromPath( params.assembly, checkIfExists: true )
params.user_genes=false
params.palidis=false

workflow {
	rename(assembly)
	prokka_annot(rename.out.contigs_1kb)

	if (params.user_genes) {
		cds_gff = Channel.fromPath( params.prot_gff, checkIfExists: true )
		cds_faa = Channel.fromPath( params.prot_fasta, checkIfExists: true )
	}else{
		cds_gff = prokka_annot.out.prokka_gff
		cds_faa = prokka_annot.out.prokka_faa
	}

	if (params.palidis) {
	        pal_seq = Channel.fromPath( params.palidis_fasta, checkIfExists: true )
        	pal_info = Channel.fromPath( params.palidis_info, checkIfExists: true )
	}else{
		pal_seq = file('no_seqs')
		pal_info = file('no_info')
	}

	diamond_mob(cds_faa, params.mobileog_db)

	gbk_split(prokka_annot.out.prokka_gbk)

        ice_finder(gbk_split.out.gbks, file( "icefinder_results/gbk"), file( "icefinder_results/tmp"), file( "icefinder_results/result"))

	integronfinder(rename.out.contigs_5kb)
	isescan(rename.out.contigs_1kb)

	integra(assembly, cds_gff, rename.out.map_file, isescan.out.iss_fasta, isescan.out.iss_tsv, pal_seq, pal_info, integronfinder.out.inf_summ, integronfinder.out.inf_gbk.collect(), ice_finder.out.icf_summ_files, ice_finder.out.icf_fasta_files, diamond_mob.out.blast_out)	


}



