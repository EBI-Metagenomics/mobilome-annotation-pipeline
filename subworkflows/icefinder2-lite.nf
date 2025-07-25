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
    ch_ice_hmm_models       // channel: tuple( val(db_id), path(ice_hmm_models_files) )
    ch_ice_macsy_models     // channel: path(ice_macsy_models)
    ch_prokka_uniprot_db    // channel: tuple( val(db_id), path(prokka_uniprot_db) )

    main:
    // Prescanning for candidate contigs
    PRODIGAL( ch_assembly )
    PRESCAN( PRODIGAL.out.faa, ch_ice_hmm_models)
    PRESCAN_TO_FASTA( PRESCAN.out.hmmscan_tbl.join( PRODIGAL.out.faa ).join( ch_assembly ) )


    // Evaluating if candidates were found
    ch_candidates = PRESCAN_TO_FASTA.out.candidates_faa
        .join(PRESCAN_TO_FASTA.out.candidates_fna)
        .filter { meta, faa, fna ->
            faa.exists() && fna.exists()
        }

    // Only run downstream processes if candidates exist
    if (ch_candidates) {
        MACSYFINDER( ch_candidates.map { meta, faa, _fna -> tuple(meta, faa) }, ch_ice_macsy_models )
        BLASTP_PROKKA( ch_candidates.map { meta, faa, _fna -> tuple(meta, faa) }, ch_prokka_uniprot_db )
        PROCESS_BLASTP_PROKKA( BLASTP_PROKKA.out.uniprot_tsv )
        ARAGORN( ch_candidates.map { meta, _faa, fna -> tuple(meta, fna) } )
        TRNAS_INTEGRATOR( ARAGORN.out.rnas_tbl.join( PRODIGAL.out.annot ) )
        VMATCH( ch_candidates.map { meta, _faa, fna -> tuple(meta, fna) } )

        // REFINE_BOUNDARIES - this will only emit results for samples with ICEs
        REFINE_BOUNDARIES(
            ch_assembly.join(
            TRNAS_INTEGRATOR.out.merged_gff).join(
            MACSYFINDER.out.macsyfinder_tsv).join(
            PROCESS_BLASTP_PROKKA.out.uniprot_product_names).join(
            VMATCH.out.vmatch_tsv)
        )
        
        final_ices_tsv = REFINE_BOUNDARIES.out.ices_tsv
    } else {
        final_ices_tsv = Channel.empty()
    }

    emit:
    ices_tsv = final_ices_tsv  // This will be empty for samples without ICEs

}


