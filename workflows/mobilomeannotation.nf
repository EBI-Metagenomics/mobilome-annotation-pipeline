include { validateParameters ; paramsHelp ; samplesheetToList } from 'plugin/nf-schema'

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

// include { VIRIFY_QC    } from './modules/virify_qc'
// include { CRISPR_FINDER } from './modules/crisprcas'

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
    
    validateParameters()

    def ch_inputs = Channel.fromList(samplesheetToList(params.input, "./assets/schema_input.json"))

    // PREPROCESSING
    RENAME( ch_inputs )

    PROKKA( RENAME.out.contigs_1kb )

    // PREDICTION
    GENOMAD( RENAME.out.contigs_5kb )
    
    GBK_SPLITTER( PROKKA.out.prokka_gbk )

    ICEFINDER( 
        GBK_SPLITTER.out.intput_list, GBK_SPLITTER.out.gbks
    )

    // remove ICEFINDER tmp directory
    // CLEANUP(
    //     file( "${params.outdir}/prediction/icefinder_results/tmp"),
    //     ICEFINDER.out.icf_summ_files
    // )

    INTEGRONFINDER( RENAME.out.contigs_5kb )
    ISESCAN( RENAME.out.contigs_1kb )

    // ANNOTATION
    DIAMOND( PROKKA.out.prokka_faa, params.mobileog_db )

    // TODO: disabled
    // if (params.virify){
    //     VIRIFY_QC(
    //         Channel.fromPath( params.vir_gff, checkIfExists: true ),
    //         file(params.vir_checkv, checkIfExists: true)
    //     )
    //     virify_results = VIRIFY_QC.out.virify_hq
    // } else {
    //     virify_results = file('no_virify')
    // }

    // if (params.skip_crispr) {
    //     crispr_tsv = file('no_crispr')
    // } else {
    //     CRISPR_FINDER( RENAME.out.contigs_1kb )
    //     crispr_tsv = CRISPR_FINDER.out.crispr_report
    // }

    // TODO: I've removed from the Integrator - palidis has to be in the samplesheet too
    // if ( params.palidis ) {
    //     pal_info = Channel.fromPath( params.palidis_info, checkIfExists: true )
    // }else{
    //     pal_info = 
    // }

    // INTEGRATION
    INTEGRATOR(
        PROKKA.out.prokka_gff.join(
            RENAME.out.map_file
        ).join(
            ISESCAN.out.iss_tsv
        ).join(
            INTEGRONFINDER.out.contigs_summary,
        ).join(
            INTEGRONFINDER.out.contigs_gbks.collect(), // TODO: I don't think this will work
        ).join(
            ICEFINDER.out.icf_summ_files
        ).join(
            ICEFINDER.out.icf_dr
        ).join(
            DIAMOND.out.blast_out
        ).join(
            GENOMAD.out.genomad_vir
        ).join(
            GENOMAD.out.genomad_vir
        ).join(
            ch_inputs.map { meta, _assembly, _user_proteins_gff, virify_gff -> {
                    [meta, virify_gff]
                }
            }
        ).join(
            ch_inputs.map { meta, _assembly, user_proteins_gff, _virify_gff -> {
                    [meta, user_proteins_gff]
                }
            }
        )
    )

    // POSTPROCESSING
    GFF_REDUCE( INTEGRATOR.out.mobilome_prokka_gff )

    FASTA_WRITER(
        ch_inputs.map {}.join( GFF_REDUCE.out.mobilome_nogenes )
    )

    if ( params.user_genes ) {
        user_gff = Channel.fromPath( params.prot_gff, checkIfExists: true )
        GFF_MAPPING(
            GFF_REDUCE.out.mobilome_clean.join( user_gff )
        )
    }

    if ( params.gff_validation ) {
        GFF_VALIDATOR( GFF_REDUCE.out.mobilome_nogenes )		
    }

    if ( !params.skip_amr ) {
        AMRFINDER_PLUS( PROKKA.out.prokka_fna.join( PROKKA.out.prokka_faa ).join( PROKKA.out.prokka_gff) )
        if ( params.user_genes ) {
            // TODO, this p
            AMRFINDER_REPORT(
                AMRFINDER_PLUS.out.amrfinder_tsv.join(
                    INTEGRATOR.out.mobilome_prokka_gff
                ).join(
                    RENAME.out.map_file
                ).join(
                    ch_inputs.map { meta, _assembly, user_proteins_gff, _virify_gff -> {
                            [meta, user_proteins_gff]
                        }
                    }
                )
            )
        } else {
            user_gff = file('no_user_gff')
            AMRFINDER_REPORT(
                AMRFINDER_PLUS.out.amrfinder_tsv.join(
                    INTEGRATOR.out.mobilome_prokka_gff
                ).join(
                    RENAME.out.map_file
                ).join(
                    ch_inputs.map { meta, _assembly, user_proteins_gff, _virify_gff -> {
                            [meta, user_proteins_gff]
                        }
                    }
                )
            )
        }
    }
}
