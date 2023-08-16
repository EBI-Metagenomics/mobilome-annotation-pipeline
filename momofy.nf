#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// MODULES
include { AMRFINDER_PLUS } from './modules/amrfinder_plus'
include { DIAMOND } from './modules/diamond'
include { GBK_SPLITTER } from './modules/gbk_splitter'
include { GFF_VALIDATOR } from './modules/validator'
include { INTEGRONFINDER } from './modules/integronfinder'
include { ISESCAN } from './modules/isescan'
include { ICEFINDER } from './modules/icefinder'
include { INTEGRATOR } from './modules/integrator'
include { PROKKA } from './modules/prokka'
include { RENAME } from './modules/rename_contigs'
include { VIRI_PARSER } from './modules/virify_parser'

def helpMessage() {
  log.info """
	MoMofy is a wraper that integrates the ouptput of different tools designed for the prediction of autonomous integrative mobile genetic elements in prokaryotic genomes and metagenomes.

        Usage:
         The basic command for running the pipeline is as follows:

         nextflow run momofy.nf --assembly contigs.fasta

         Mandatory arguments:
          --assembly                      (Meta)genomic assembly in fasta format (uncompress)

         Optional arguments:
          --user_genes                    Use the user annotation files. See --prot_fasta and --prot_gff [false]
          --prot_gff                      Annotation file in GFF3 format. Mandatory with --user_genes true
          --prot_fasta                    Fasta file of aminoacids. Mandatory with --user_genes true
          --virify                        Use VIRify results. See --vir_gff [false]
          --vir_gff                       The full path of VIRify results on GFF format. Mandatory with --vir_gff true
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
params.virify = false
params.palidis = false
params.gff_validation = true
params.outdir = 'MoMofy_results'

workflow {

	assembly = Channel.fromPath( params.assembly, checkIfExists: true )

	RENAME( assembly )

	if (params.virify){
		VIRI_PARSER( params.assembly, params.vir_gff )
	}

	PROKKA( RENAME.out.contigs_1kb )

	AMRFINDER_PLUS( PROKKA.out.prokka_fna, PROKKA.out.prokka_faa, PROKKA.out.prokka_gff )

	if ( params.user_genes ) {
		cds_gff = Channel.fromPath( params.prot_gff, checkIfExists: true )
		cds_faa = Channel.fromPath( params.prot_fasta, checkIfExists: true )
	}else{
		cds_gff = PROKKA.out.prokka_gff
		cds_faa = PROKKA.out.prokka_faa
	}

	if ( params.palidis ) {
	        pal_seq = Channel.fromPath( params.palidis_fasta, checkIfExists: true )
        	pal_info = Channel.fromPath( params.palidis_info, checkIfExists: true )
	}else{
		pal_seq = file('no_seqs')
		pal_info = file('no_info')
	}

	DIAMOND( cds_faa, params.mobileog_db )

	GBK_SPLITTER( PROKKA.out.prokka_gbk )

        ICEFINDER( GBK_SPLITTER.out.gbks, file( "icefinder_results/gbk"), file( "icefinder_results/tmp"), file( "icefinder_results/result" ) )

	INTEGRONFINDER( RENAME.out.contigs_5kb )

	ISESCAN( RENAME.out.contigs_1kb )

	INTEGRATOR(cds_gff, RENAME.out.map_file, ISESCAN.out.iss_fasta, ISESCAN.out.iss_tsv, pal_seq, pal_info, INTEGRONFINDER.out.inf_summ, INTEGRONFINDER.out.inf_gbk.collect(), ICEFINDER.out.icf_summ_files, ICEFINDER.out.icf_fasta_files, ICEFINDER.out.icf_dr, DIAMOND.out.blast_out)	

	if (params.gff_validation) {
		GFF_VALIDATOR(INTEGRATOR.out.momo_gff)		
	}

}



