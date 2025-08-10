include { validateParameters ; paramsHelp ; samplesheetToList } from 'plugin/nf-schema'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Inputs preparing modules
include { RENAME           } from '../modules/rename_contigs'

// Annotation modules
include { PROKKA           } from '../modules/prokka'
include { AMRFINDER_PLUS   } from '../modules/amrfinder_plus'

// Mobile genetic elements prediction modules
include { INTEGRONFINDER   } from '../modules/integronfinder'
include { ISESCAN          } from '../modules/isescan'
include { GENOMAD          } from '../modules/genomad'
include { VIRIFY_QC        } from '../modules/virify_qc'

// Results integration and writing modules
include { AMRFINDER_REPORT } from '../modules/amrfinder_report'
include { FASTA_WRITER     } from '../modules/fasta_writer'
include { GFF_MAPPING      } from '../modules/gff_mapping'
include { GFF_REDUCE       } from '../modules/gff_reduce'
include { GFF_VALIDATOR    } from '../modules/gff_validator'
include { INTEGRATOR       } from '../modules/integrator'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { COMPOSITIONAL_OUTLIER_DETECTION } from '../subworkflows/compositional_outlier_detection'
include { ICEFINDER2_LITE                 } from '../subworkflows/icefinder2-lite'
// add blast annotation workflow from mobilome proteins after integration


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MOBILOMEANNOTATION WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MOBILOMEANNOTATION {

    validateParameters()

    def ch_inputs = Channel.fromList(samplesheetToList(params.input, "./assets/schema_input.json"))

    /*
    ******************************************************************************************************
    * The code below is transforming the input channels to handle optional inputs, such as the
    * user-provided GFF files for user_proteins and virify.
    * Nextflow doesn't handle optional inputs well, so we use a common hack to provide an empty
    * array ([]) when the input is missing.
    * For the user_proteins_gff, if the file is present, we emit a tuple with the metadata and the
    * file path. If the file is missing, we emit a tuple with the metadata and an empty array ([]).
    * Similarly, for the virify_gff, if the file is present, we emit a tuple with the metadata and
    * the file path. If the file is missing, we emit a tuple with the metadata and an empty array ([]).
    ******************************************************************************************************
    */

    def user_proteins_ch = ch_inputs.map { meta, _fasta, user_proteins_gff, _virify_gff -> {
            if ( user_proteins_gff ) {
                [meta, user_proteins_gff ]
            } else {
                [meta, []]
            }
        }
    }


    // PREPROCESSING
    RENAME( ch_inputs.map { meta, fasta, _user_proteins_gff, _virify_gff -> [meta, fasta] } )

    PROKKA( RENAME.out.contigs_1kb )

    // Parsing VIRify gff file when an input is provided
    def user_virify_gff_ch = ch_inputs.map { meta, _fasta, _user_proteins_gff, virify_gff -> {
           [meta, virify_gff]
        }
    }.filter { _meta, virify_gff -> virify_gff != [] }

    VIRIFY_QC( user_virify_gff_ch )

    // PREDICTION
    // Collecting ICEfinder2 databases
    db_ice_hmm_models = Channel.fromPath("${params.ice_hmm_models}.*", checkIfExists: true)
                        .collect()
                        .map { ice_db_files ->
                            [[id: file(params.ice_hmm_models).name ], ice_db_files]
                        }
    db_prokka_uniprot = Channel.fromPath("${params.prokka_uniprot_db}.*", checkIfExists: true)
                        .collect()
                        .map { uniprot_db_files ->
                            [[id: file(params.prokka_uniprot_db).name], uniprot_db_files]
                        }


    ICEFINDER2_LITE( 
        RENAME.out.contigs_5kb,
        db_ice_hmm_models, 
        params.ice_macsy_models, 
        db_prokka_uniprot
    )

    GENOMAD( RENAME.out.contigs_5kb )

    INTEGRONFINDER( RENAME.out.contigs_5kb )

    ISESCAN( RENAME.out.contigs_1kb )

    COMPOSITIONAL_OUTLIER_DETECTION( RENAME.out.contigs_100kb )


    /**********************************************************************************************
    * The INTEGRATOR step takes a bunch of outputs from the previous steps.
    * The following code is re-shaping the input to accommodate
    * optional inputs such as the user-provided GFF.
    * This is done this way because Nextflow doesn't handle optional inputs. One hack that the
    * community uses for inputs of type path is to provide an empty array ([]). So, we first
    * join with user-provided GFF with the remainder, try to get an empty element, and then we use map
    * to transform the null to [].
    ***********************************************************************************************/
    def integrator_ch = PROKKA.out.prokka_gff.join(
        RENAME.out.map_file
    ).join(
        ISESCAN.out.iss_tsv
    ).join(
        INTEGRONFINDER.out.contigs_summary
    ).join(
        INTEGRONFINDER.out.contigs_gbks
    ).join(
        ICEFINDER2_LITE.out.ices_tsv, remainder: true
    ).join(
        GENOMAD.out.genomad_vir
    ).join(
        GENOMAD.out.genomad_plas
    ).join(
        COMPOSITIONAL_OUTLIER_DETECTION.out.bed, remainder: true
    ).join(
        VIRIFY_QC.out.virify_hq, remainder: true
    )

    INTEGRATOR(
        integrator_ch.map {
            meta, prokka_gff, map_file, iss_tsv, contigs_summary, gbks, ices_tsv, genomad_vir, genomad_plas, compos_bed, virify_hq -> {
                [meta, prokka_gff, map_file, iss_tsv, contigs_summary, gbks, ices_tsv ? ices_tsv : [], genomad_vir, genomad_plas, compos_bed ? compos_bed : [], virify_hq ? virify_hq : [] ]
            }
        }
    )


    // POSTPROCESSING
    GFF_REDUCE( INTEGRATOR.out.mobilome_prokka_gff )

    FASTA_WRITER(
        ch_inputs.map { meta, fasta, _user_proteins_gff, _virify_gff -> [meta, fasta] } .join( GFF_REDUCE.out.mobilome_nogenes )
    )

    GFF_MAPPING(
        GFF_REDUCE.out.mobilome_clean.join( user_proteins_ch )
    )

    if ( params.gff_validation ) {
        GFF_VALIDATOR( GFF_REDUCE.out.mobilome_nogenes )		
    }
    
    // AMRFinder is optional. default skip_amr = FALSE
    def amr_finder_ch = PROKKA.out.prokka_fna.join( PROKKA.out.prokka_faa ).join( PROKKA.out.prokka_gff).filter({ it -> !it[0].skip_amrfinder_plus })

    AMRFINDER_PLUS( amr_finder_ch )

    AMRFINDER_REPORT(
        AMRFINDER_PLUS.out.amrfinder_tsv.join(
            INTEGRATOR.out.mobilome_prokka_gff
        ).join(
            RENAME.out.map_file
        ).join(
            user_proteins_ch
        )
    )
}
