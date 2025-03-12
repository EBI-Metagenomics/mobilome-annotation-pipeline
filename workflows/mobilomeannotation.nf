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
include { CRISPR_FINDER    } from '../modules/crisprcas'

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
    RENAME( ch_inputs.map { meta, fasta, _user_proteins_gff, _virify_gff -> [meta, fasta] } )

    PROKKA( RENAME.out.contigs_1kb )

    // PREDICTION
    GENOMAD( RENAME.out.contigs_5kb )
    
    GBK_SPLITTER( PROKKA.out.prokka_gbk )

    ICEFINDER( 
        GBK_SPLITTER.out.intput_list.join( GBK_SPLITTER.out.gbks )
    )

    INTEGRONFINDER( RENAME.out.contigs_5kb )

    ISESCAN( RENAME.out.contigs_1kb )

    // ANNOTATION
    DIAMOND( PROKKA.out.prokka_faa, file(params.mobileog_db, checkIfExists: true) )

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

    
    CRISPR_FINDER( RENAME.out.contigs_1kb.filter { it -> !it[0].skip_crispr_finder } )

    // TODO: I've removed from the Integrator - palidis has to be in the samplesheet too
    // if ( params.palidis ) {
    //     pal_info = Channel.fromPath( params.palidis_info, checkIfExists: true )
    // } else {
    //     pal_info = 
    // }
    
    def integrator_ch = PROKKA.out.prokka_gff.join(
        RENAME.out.map_file
    ).join(
        ISESCAN.out.iss_tsv
    ).join(
        INTEGRONFINDER.out.contigs_summary,
    ).join(
        INTEGRONFINDER.out.contigs_gbks.collect(),
    ).join(
        ICEFINDER.out.icf_summ_files
    ).join(
        ICEFINDER.out.icf_dr
    ).join(
        DIAMOND.out.blast_out
    ).join(
        GENOMAD.out.genomad_vir
    ).join(
        GENOMAD.out.genomad_plas
    )

    INTEGRATOR(
        integrator_ch,
        [[], []],
        [[], []]
    )

    // POSTPROCESSING
    GFF_REDUCE( INTEGRATOR.out.mobilome_prokka_gff )

    // TODO should the input of fasta writer be the renamed fasta files?
    FASTA_WRITER(
        ch_inputs.map { meta, fasta, _user_proteins_gff, _virify_gff -> [meta, fasta] } .join( GFF_REDUCE.out.mobilome_nogenes )
    )

    def user_proteins_ch = ch_inputs.map { meta, _fasta, user_proteins_gff, _virify_gff -> [meta, user_proteins_gff] }

    GFF_MAPPING(
        GFF_REDUCE.out.mobilome_clean,
        [[],[]] // TODO: user_proteins_ch
    )

    if ( params.gff_validation ) {
        GFF_VALIDATOR( GFF_REDUCE.out.mobilome_nogenes )		
    }
    
    // AMRFinder is optional
    def amr_finder_ch = PROKKA.out.prokka_fna.join( PROKKA.out.prokka_faa ).join( PROKKA.out.prokka_gff).filter({ it -> !it[0].skip_amrfinder_plus })

    AMRFINDER_PLUS( amr_finder_ch )

    AMRFINDER_REPORT(
        AMRFINDER_PLUS.out.amrfinder_tsv.join(
                INTEGRATOR.out.mobilome_prokka_gff
        ).join(
            RENAME.out.map_file
        ),
        [[],[]] // TODO: add the user proteins if present
    )
}
