include { PRODIGAL              } from '../modules/prodigal'
include { HMMSCAN as PRESCAN    } from '../modules/hmmscan'
include { PRESCAN_TO_FASTA      } from '../modules/prescan_to_fasta'
include { BLASTP_PROKKA         } from '../modules/blastp_prokka'
include { PROCESS_BLASTP_PROKKA } from '../modules/process_blastp_prokka'
include { ARAGORN               } from '../modules/aragorn'
include { TRNAS_INTEGRATOR      } from '../modules/trnas_integrator'
include { MACSYFINDER           } from '../modules/macsyfinder'
include { VMATCH                } from '../modules/vmatch'
include { REFINE_BOUNDARIES     } from '../modules/ice_refine_boundaries'

workflow ICEFINDER2_LITE {
    take:
    ch_assembly             // channel: tuple( val(meta), path(assembly_5kb) )
    ch_ice_hmm_models       // channel: tuple( val(db_id), path(ice_hmm_models_files))
    ch_ice_macsy_models     // channel: path(ice_macsy_models)
    ch_prokka_uniprot_db    // channel: tuple( val(db_id), path(prokka_uniprot_db))

    main:
    // Preannotation to detect candidate contigs
    PRODIGAL( ch_assembly )

    PRESCAN( PRODIGAL.out.faa, ch_ice_hmm_models)

    PRESCAN_TO_FASTA( PRESCAN.out.hmmscan_tbl.join( PRODIGAL.out.faa ).join( ch_assembly ) )

    // Actual annotation of candidate contigs 
    MACSYFINDER( PRESCAN_TO_FASTA.out.candidates_faa, ch_ice_macsy_models )

    BLASTP_PROKKA( PRESCAN_TO_FASTA.out.candidates_faa, ch_prokka_uniprot_db )

    PROCESS_BLASTP_PROKKA( BLASTP_PROKKA.out.uniprot_tsv )

    ARAGORN( PRESCAN_TO_FASTA.out.candidates_fna )

    TRNAS_INTEGRATOR( ARAGORN.out.rnas_tbl.join( PRODIGAL.out.annot )  )

    VMATCH( PRESCAN_TO_FASTA.out.candidates_fna )


    // Boundaries refinement of predicted ICEs
    REFINE_BOUNDARIES(
        ch_assembly.join(
        TRNAS_INTEGRATOR.out.merged_gff).join(
        MACSYFINDER.out.macsyfinder_tsv).join(
        PROCESS_BLASTP_PROKKA.out.uniprot_product_names).join(
        VMATCH.out.vmatch_tsv)
    )
    
    emit:
    ices_tsv = REFINE_BOUNDARIES.out.ices_tsv
 
}
