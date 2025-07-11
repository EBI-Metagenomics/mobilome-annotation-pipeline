include { MACSYFINDER           } from '../modules/macsyfinder'
include { MACSYFINDER_PROCESS   } from '../modules/macsyfinder_processor'
include { ARAGORN               } from '../modules/aragorn'
include { VMATCH                } from '../modules/vmatch'
include { REFINE_BOUNDARIES     } from '../modules/ice_refine_boundaries'
include { VALIDATE_ICE_ELEMENTS } from '../modules/validate_ice_elements'


workflow ICEFINDER2_LITE {
    take:
    ch_input_files     //  channel: [val(meta), path(fna), path(faa), path(gff)]
    ch_ice_models   // channel: path(ice_models)

    main:
    // Step 1: MacSyFinder ICE detection on proteins
    MACSYFINDER( ch_input_files.map{ meta, _fna, faa, _gff -> [meta, faa] },
        ch_ice_models
    )
 
    MACSYFINDER_PROCESS(MACSYFINDER.out.macsyfinder_tsv.join( ch_input_files.map{ meta, fna, _faa, gff -> [meta, fna, gff] } ))

    // Step 3: ICE boundary refinement using direct repeats and tRNA analysis
    ARAGORN(MACSYFINDER_PROCESS.out.all_sys_flanks_fasta)

    VMATCH(MACSYFINDER_PROCESS.out.all_sys_flanks_fasta)
   
    REFINE_BOUNDARIES(
        MACSYFINDER_PROCESS.out.boundaries_tsv.join(
        ARAGORN.out.trna_gff).join(
        VMATCH.out.vmatch_tsv)
    )
    

    // Step 4: Validate final ICE elements
    VALIDATE_ICE_ELEMENTS(
        REFINE_BOUNDARIES.out.refined_tsv
    )
    
    emit:
    ices_tsv = VALIDATE_ICE_ELEMENTS.out.validated_ices

}
