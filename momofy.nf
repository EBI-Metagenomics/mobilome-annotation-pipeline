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
include { subfastq } from './modules/splitfastq'
include { palidis } from './software/my_palidis/palidis'
include { integra } from './modules/integrator'

// INPUTS
cds_annot = Channel.fromPath( params.cds_annot, checkIfExists: true )
assembly = Channel.fromPath( params.assembly, checkIfExists: true )

params.palidis=false
pal_seq=file('no_seqs')
pal_info=file('no_info')
if (params.palidis) {
	read_1 = Channel.fromPath( params.read_1, checkIfExists: true )
	read_2 = Channel.fromPath( params.read_2, checkIfExists: true )
}

workflow {
	diamond_mob(cds_annot, params.mobileog_db)
	rename(assembly)
	prokka_annot(rename.out.contigs_5kb)
	gbk_split(prokka_annot.out.prokka_gbk)

        ice_finder(gbk_split.out.gbks, file( "icefinder_results/gbk"), file( "icefinder_results/tmp"), file( "icefinder_results/result"))

	integronfinder(rename.out.contigs_5kb)
	isescan(rename.out.contigs_1kb)

	if (params.palidis) {
		subfastq(read_1, read_2)
	        p1_ch = subfastq.out.pair_1.map { tuple('sample_1', it) }
                p2_ch = subfastq.out.pair_2.map { tuple('sample_1', it) }
               	palidis_input_reads = p1_ch.join(p2_ch)
	        contig_ch = rename.out.contigs_1kb.map { tuple('sample_1', it) }

        	palidis(palidis_input_reads, contig_ch)

		pal_seq=palidis.out.is_fasta_ch.ifEmpty(file('no_seqs'))
		pal_info=palidis.out.is_info_ch.ifEmpty(file('no_info'))
	}
	
	integra(cds_annot, rename.out.map_file, isescan.out.iss_fasta, isescan.out.iss_tsv, pal_seq, pal_info, integronfinder.out.inf_summ, integronfinder.out.inf_gbk.collect(), ice_finder.out.icf_summ_files, ice_finder.out.icf_fasta_files, diamond_mob.out.blast_out)	


}



