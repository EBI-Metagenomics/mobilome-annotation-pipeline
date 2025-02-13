include { validateParameters ; paramsHelp ; paramsSummaryLog ; samplesheetToList ; paramsSummaryMap } from 'plugin/nf-schema'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Inputs preparing modules
include { RENAME           } from '../modules/rename_contigs'
include { GBK_SPLITTER     } from '../modules/gbk_splitter'

// Protein annotation modules
include { PROKKA           } from '../modules/prokka'
include { DIAMOND          } from '../modules/diamond'
include { AMRFINDER_PLUS   } from '../modules/amrfinder_plus'

// Mobile genetic elements annotation modules
include { INTEGRONFINDER   } from '../modules/integronfinder'
include { ISESCAN          } from '../modules/isescan'
include { ICEFINDER        } from '../modules/icefinder'
include { GENOMAD          } from '../modules/genomad'
//include { VIRIFY_QC } from './modules/virify_qc'
//include { CRISPR_FINDER } from './modules/crisprcas'

// Results integration and writing modules
include { AMRFINDER_REPORT } from '../modules/amrfinder_report'
include { FASTA_WRITER     } from '../modules/fasta_writer'
include { GFF_MAPPING      } from '../modules/gff_mapping'
include { GFF_REDUCE       } from '../modules/gff_reduce'
include { GFF_VALIDATOR    } from '../modules/validator'
include { INTEGRATOR       } from '../modules/integrator'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MOBILOMEANNOTATION {
    def summary_params = paramsSummaryMap(workflow)
    validateParameters()
    log.info(paramsSummaryLog(workflow))

    if (params.help) {
        log.info(paramsHelp("nextflow run ebi-metagenomics-mobilomeannotation --help"))
        exit(0)
    }


    ch_versions = Channel.empty()


    // ---- Combine input files into the inputs channel ---- //	
    def groupInputs = { list ->
        def meta = [id: list[0]]
        def assembly = list[1]
	def user_prots = list[2]
	def virify_gff = list[3]
        return tuple(meta, assembly, user_prots, virify_gff)
    }

    ch_inputs = Channel.fromList(samplesheetToList(params.input, "./assets/schema_input.json")).map(groupInputs) // [ meta, assembly, user_prots, virify_gff ]

    ch_inputs.toList().each { input ->
        println "Meta: ${input[0]}, Assembly: ${input[1]}, User Proteins: ${input[2]}, Virify GFF: ${input[3]}"
    }


    /*

    // PREPROCESSING
    RENAME( ch_inputs )

    PROKKA( RENAME.out.contigs_1kb )


    // PREDICTION
    GENOMAD( RENAME.out.contigs_5kb )
    GBK_SPLITTER( PROKKA.out.prokka_gbk )
	ICEFINDER( 
		GBK_SPLITTER.out.gbks,
		file( "${params.outdir}/prediction/icefinder_results/gbk"),
		file( "${params.outdir}/prediction/icefinder_results/tmp"),
		file( "${params.outdir}/prediction/icefinder_results/result")
	)
    // remove ICEFINDER tmp directory
    CLEANUP(
        file( "${params.outdir}/prediction/icefinder_results/tmp"),
        ICEFINDER.out.icf_summ_files
    )

    INTEGRONFINDER( RENAME.out.contigs_5kb )
    ISESCAN( RENAME.out.contigs_1kb )

    // ANNOTATION
    DIAMOND( PROKKA.out.prokka_faa, params.mobileog_db )

    if (params.virify){
        gff_input = Channel.fromPath( params.vir_gff, checkIfExists: true )
        checkv_input = file(params.vir_checkv, checkIfExists: true)
        VIRIFY_QC( gff_input, checkv_input )
        virify_results = VIRIFY_QC.out.virify_hq
    } else {
        virify_results = file('no_virify')
    }

    if (params.skip_crispr){
        crispr_tsv = file('no_crispr')
    } else {
        CRISPR_FINDER( RENAME.out.contigs_1kb )
        crispr_tsv = CRISPR_FINDER.out.crispr_report
    }

    if ( params.palidis ) {
        pal_info = Channel.fromPath( params.palidis_info, checkIfExists: true )
    }else{
        pal_info = file('no_info')
    }

    // INTEGRATION
    INTEGRATOR(
        PROKKA.out.prokka_gff, 
        RENAME.out.map_file, 
        ISESCAN.out.iss_tsv, 
        pal_info, 
        INTEGRONFINDER.out.inf_summ, 
        INTEGRONFINDER.out.inf_gbk.collect(), 
        ICEFINDER.out.icf_summ_files, 
        ICEFINDER.out.icf_dr, 
        DIAMOND.out.blast_out,
        GENOMAD.out.genomad_vir,
        GENOMAD.out.genomad_plas,
        virify_results,
        crispr_tsv,
    )


    // POSTPROCESSING
    GFF_REDUCE( INTEGRATOR.out.mobilome_prokka_gff )
    FASTA_WRITER( assembly, GFF_REDUCE.out.mobilome_nogenes )

    if ( params.user_genes ) {
        user_gff = Channel.fromPath( params.prot_gff, checkIfExists: true )
        GFF_MAPPING( GFF_REDUCE.out.mobilome_clean, user_gff )
    }

    if ( params.gff_validation ) {
        GFF_VALIDATOR( GFF_REDUCE.out.mobilome_nogenes )		
    }

    if ( !params.skip_amr ) {
        AMRFINDER_PLUS( PROKKA.out.prokka_fna, PROKKA.out.prokka_faa, PROKKA.out.prokka_gff )
	if ( params.user_genes ) {
            user_gff = Channel.fromPath( params.prot_gff, checkIfExists: true )
            AMRFINDER_REPORT( AMRFINDER_PLUS.out.amrfinder_tsv, INTEGRATOR.out.mobilome_prokka_gff, RENAME.out.map_file, user_gff )
	} else {
            user_gff = file('no_user_gff')
            AMRFINDER_REPORT( AMRFINDER_PLUS.out.amrfinder_tsv, INTEGRATOR.out.mobilome_prokka_gff, RENAME.out.map_file, user_gff )
	}
    }

    */

}
