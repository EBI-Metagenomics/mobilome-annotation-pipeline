include { MACSYFINDER         } from '../modules/macsyfinder'
include { MACSYFINDER_PROCESS } from '../modules/macsyfinder_processor'
include { ARAGORN             } from '../modules/aragorn'
include { VMATCH              } from '../modules/vmatch'
include { VMATCH_PROCESS      } from '../modules/vmatch_processor'

include { REFINE_BOUNDARIES } from '../modules/ice_refine_boundaries'
include { VALIDATE_ICE_ELEMENTS } from '../modules/validate_ice_elements'
include { ORIT_DETECTION      } from '../modules/orit_detection'
include { ICE_INFORMATION_INTEGRATION } from '../modules/ice_information_integration'
include { ICE_CLASSIFICATION_TYPING } from '../modules/ice_classification_typing'
include { ICE_QUALITY_ASSESSMENT } from '../modules/ice_quality_assessment'
include { GENERATE_GFF3_OUTPUT } from '../modules/ice_gff3_output'


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
   
    VMATCH_PROCESS(VMATCH.out.vmatch_raw)
 

    

    REFINE_BOUNDARIES(
        MACSYFINDER_PROCESS.out.boundaries_tsv.join(
        ARAGORN.out.trna_gff).join(
        VMATCH_PROCESS.out.dr_tsv)
    )
    
    /*

    // Step 6: Validate final ICE elements
    VALIDATE_ICE_ELEMENTS(
        REFINE_BOUNDARIES.out.refined_boundaries,
        DETECT_TRNA.out.trna_annotations
    )
    

    // Step 4: oriT detection
    ORIT_DETECTION(
        VALIDATE_ICE_ELEMENTS.out.validated_ice,
    )


    // Step 5: Integrating results
    ICE_INFORMATION_INTEGRATION(
        validated_ice_elements,
        orit_results,
        ice_sequences,
        genome_fasta,
        macsyfinder_gff3,
        trna_annotations,
        direct_repeats
    )
    
    ICE_CLASSIFICATION_TYPING(
        ICE_INFORMATION_INTEGRATION.out.integrated_ice_info,
        macsyfinder_gff3
    )
    
    ICE_QUALITY_ASSESSMENT(
        ICE_CLASSIFICATION_TYPING.out.classified_ice_elements,
        ICE_INFORMATION_INTEGRATION.out.integrated_ice_info
    )
    
    // Step 6: Generate final GFF3 output
    GENERATE_GFF3_OUTPUT(
        ICE_QUALITY_ASSESSMENT.out.quality_assessed_ice,
        ICE_INFORMATION_INTEGRATION.out.integrated_ice_info,
        genome_fasta
    )
    
    


    emit:
    final_ice_gff3 = GENERATE_GFF3_OUTPUT.out.ice_gff3              // Final ICE annotations in GFF3 format
    ice_summary = ICE_QUALITY_ASSESSMENT.out.ice_summary            // Summary statistics
    classified_elements = ICE_CLASSIFICATION_TYPING.out.classified_ice_elements  // Classified ICE elements
    quality_report = ICE_QUALITY_ASSESSMENT.out.quality_report      // Quality assessment report
    integrated_info = ICE_INFORMATION_INTEGRATION.out.integrated_ice_info  // Comprehensive ICE information
    */


    emit:
    gff3_output = REFINE_BOUNDARIES.out.refined_boundaries



}
