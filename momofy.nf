#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// MODULES
include { diamond_mob } from './modules/diamond'
include { rename } from './modules/rename_contigs'
include { prokka_annot } from './modules/prokka'
include { gbk_split } from './modules/gbk_splitter'
include { integronfinder } from './modules/integronfinder'
include { isescan } from './modules/isescan'
include { ice_finder } from './modules/icefinder'
include { integra } from './modules/integrator'
include { gff_validator } from './modules/validator'


def helpMessage() {
  log.info """
	MoMofy is a wraper that integrates the ouptput of different tools designed for the prediction of autonomous integrative mobile genetic elements in prokaryotic genomes and metagenomes.

        Usage:
         The basic command for running the pipeline is as follows:

         nextflow run momofy.nf --assembly contigs.fasta

         Mandatory arguments:
          --assembly                     (Meta)genomic assembly in fasta format (uncompress)

         Optional arguments:
          --user_genes                    User annotation files. See --prot_fasta and --prot_gff [false]
          --prot_gff                      Annotation file in GFF3 format. Mandatory with --user_genes true
          --prot_fasta                    Fasta file of aminoacids. Mandatory with --user_genes true
          --palidis                       Incorporate PaliDIS predictions to final output [false]
          --palidis_fasta                 Fasta file of PaliDIS insertion sequences. Mandatory with --palidis true
          --palidis_info                  Information file of PaliDIS insertion sequences. Mandatory with --palidis true
          --gff_validation                Run a step of format validation on the GFF3 file output [true]
          --outdir                        Output directory to place final MoMofy results [MoMofy_results]
          --help                          This usage statement [false]
        """
}

if (params.help) {
    helpMessage()
    exit 0
}

// Default options
params.user_genes = false
params.palidis = false
params.gff_validation = true
params.outdir = 'MoMofy_results'

workflow {

	assembly = Channel.fromPath( params.assembly, checkIfExists: true )
	rename( assembly )

	prokka_annot( rename.out.contigs_1kb )

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

	integra(cds_gff, rename.out.map_file, isescan.out.iss_fasta, isescan.out.iss_tsv, pal_seq, pal_info, integronfinder.out.inf_summ, integronfinder.out.inf_gbk.collect(), ice_finder.out.icf_summ_files, ice_finder.out.icf_fasta_files, ice_finder.out.icf_dr, diamond_mob.out.blast_out)	

	if (params.gff_validation) {
		gff_validator(integra.out.momo_gff)		
	}

}



