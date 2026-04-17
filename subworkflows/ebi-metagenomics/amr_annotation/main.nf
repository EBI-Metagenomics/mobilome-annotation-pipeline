// Subworkflow to generate antimicrobial resistance annotation from protein and gene files
// Outputs are standardised and integrated into a single GFF3 format output
include { AMRFINDERPLUS_UPDATE             } from '../../../modules/nf-core/amrfinderplus/update/main'
include { AMRFINDERPLUS_RUN                } from '../../../modules/nf-core/amrfinderplus/run/main'
include { DEEPARG_DOWNLOADDATA             } from '../../../modules/nf-core/deeparg/downloaddata/main'
include { DEEPARG_PREDICT                  } from '../../../modules/nf-core/deeparg/predict/main'
include { RGI_DOWNLOADDB                   } from '../../../modules/ebi-metagenomics/rgi/downloaddb/main'
include { RGI_CARDANNOTATION               } from '../../../modules/nf-core/rgi/cardannotation/main'
include { RGI_MAIN                         } from '../../../modules/nf-core/rgi/main/main'
include { HAMRONIZATION_RGI                } from '../../../modules/nf-core/hamronization/rgi/main'
include { HAMRONIZATION_DEEPARG            } from '../../../modules/nf-core/hamronization/deeparg/main'
include { AMRINTEGRATOR                    } from '../../../modules/ebi-metagenomics/amrintegrator/main'

workflow AMR_ANNOTATION {
    take:
    ch_inputs                 // channel: tuple( val(meta), path(aminoacids), path(cds_gff) )
    ch_amrfinderplus_db       // channel: path( amrfinderplus_db )
    ch_deeparg_db             // channel: path( deeparg_db )
    ch_deeparg_db_version     // channel: val( deeparg_db_version )
    ch_deeparg_model          // channel: val( deeparg_model )
    ch_deeparg_tool_version   // channel: val( deeparg_tool_version )
    ch_rgi_db                 // channel: path( rgi_db )
    skip_amrfinderplus        // boolean 
    skip_deeparg              // boolean
    skip_rgi                  // boolean

    main:
    ch_versions = channel.empty()

    // Extract individual components from input channel
    ch_faa = ch_inputs.map{ meta, aminoacids, _cds_gff ->
        // Add is_proteins flag to meta for amrfinderplus nf-core module:
        // https://github.com/nf-core/modules/blob/e3da0d1bd481776ac3ddbd652e1cd72c4d7426d5/modules/nf-core/amrfinderplus/run/main.nf#L33
        def meta_with_protein_flag = meta + [is_proteins: true]
        tuple(meta_with_protein_flag, aminoacids)
    }
    ch_gff = ch_inputs.map{ meta, _aminoacids, cds_gff -> tuple(meta, cds_gff) }

    // RUNNING ANNOTATION TOOLS
    // AMRfinderplus run
    // Prepare channel for database
    if (!skip_amrfinderplus && ch_amrfinderplus_db) {
        ch_amrfinderplus_db = ch_amrfinderplus_db
    }
    else if (!skip_amrfinderplus && !ch_amrfinderplus_db) {
        AMRFINDERPLUS_UPDATE()
        ch_versions = ch_versions.mix(AMRFINDERPLUS_UPDATE.out.versions.first())
        ch_amrfinderplus_db = AMRFINDERPLUS_UPDATE.out.db
    }

    if (!skip_amrfinderplus) {
        AMRFINDERPLUS_RUN(ch_faa, ch_amrfinderplus_db)
        ch_versions = ch_versions.mix(AMRFINDERPLUS_RUN.out.versions.first())
        ch_amrfinderplus_results = AMRFINDERPLUS_RUN.out.report
    }

    // RGI run
    if (!skip_rgi) {
        if (!ch_rgi_db) {
            // Download and process CARD database
            RGI_DOWNLOADDB()
            ch_versions = ch_versions.mix(RGI_DOWNLOADDB.out.versions.first())
            RGI_CARDANNOTATION(RGI_DOWNLOADDB.out.card_json)
            card = RGI_CARDANNOTATION.out.db
            ch_versions = ch_versions.mix(RGI_CARDANNOTATION.out.versions.first())
        }
        else {
            // Use user-supplied database
            rgi_db = ch_rgi_db
            if (!rgi_db.contains("card_database_processed")) {
                RGI_CARDANNOTATION(rgi_db)
                card = RGI_CARDANNOTATION.out.db
                ch_versions = ch_versions.mix(RGI_CARDANNOTATION.out.versions.first())
            }
            else {
                card = rgi_db
            }
        }

        RGI_MAIN(ch_faa, card, [])
        ch_versions = ch_versions.mix(RGI_MAIN.out.versions.first())

        // Reporting
        HAMRONIZATION_RGI(RGI_MAIN.out.tsv, 'tsv', RGI_MAIN.out.tool_version, RGI_MAIN.out.db_version)
        ch_versions = ch_versions.mix(HAMRONIZATION_RGI.out.versions.first())
        ch_rgi_results = HAMRONIZATION_RGI.out.tsv
    }

    // DeepARG prepare download
    if (!skip_deeparg && ch_deeparg_db) {
        deeparg_db = ch_deeparg_db
    }
    else if (!skip_deeparg && !ch_deeparg_db) {
        // Download deeparg database
        DEEPARG_DOWNLOADDATA()
        ch_versions = ch_versions.mix(DEEPARG_DOWNLOADDATA.out.versions.first())
        deeparg_db = DEEPARG_DOWNLOADDATA.out.db
    }

    // DeepARG run
    if (!skip_deeparg) {
        ch_faa
           .map { meta, annotations ->
                def model = params.ch_deeparg_model
                [meta, annotations, model]
            }
            .set { ch_input_for_deeparg }

        DEEPARG_PREDICT(ch_input_for_deeparg, deeparg_db)
        ch_versions = ch_versions.mix(DEEPARG_PREDICT.out.versions.first())

        HAMRONIZATION_DEEPARG(DEEPARG_PREDICT.out.arg, 'tsv', ch_deeparg_tool_version, ch_deeparg_db_version)
        ch_versions = ch_versions.mix(HAMRONIZATION_DEEPARG.out.versions.first())
        ch_deeparg_results = HAMRONIZATION_DEEPARG.out.tsv
    }

    // Integrate and transform into a single GFF3 output file
    // Key insight: Tool results have meta with is_proteins:true, but ch_gff has original meta
    // We need to join by ID to avoid meta mismatch

    // Create keyed channels using meta.id for consistent joining
    ch_gff_keyed = ch_gff.map { meta, gff -> [meta.id, meta, gff] }

    // Create flag channels for tools that ran (keyed by meta.id)
    ch_deeparg_flags = skip_deeparg ? 
        channel.empty() : 
        ch_deeparg_results.map { meta, _file -> [meta.id, 'deeparg'] }
    
    ch_rgi_flags = skip_rgi ? 
        channel.empty() : 
        ch_rgi_results.map { meta, _file -> [meta.id, 'rgi'] }
    
    ch_amrfinder_flags = skip_amrfinderplus ? 
        channel.empty() : 
        ch_amrfinderplus_results.map { meta, _file -> [meta.id, 'amrfinder'] }

    // Combine all flags to identify samples that have at least one tool result
    ch_all_flags = ch_deeparg_flags
        .mix(ch_rgi_flags)
        .mix(ch_amrfinder_flags)
        .groupTuple()
        .map { id, tools -> [id, true] }  // Create a simple boolean flag

    // Filter GFF channel to only samples with results
    ch_gff_filtered = ch_gff_keyed
        .join(ch_all_flags, remainder: true)
        .filter { _id, _meta, _gff, has_results ->
            has_results == true
        }
        .map { id, meta, gff, _has_results ->
            [id, meta, gff]
        }

    // Build the input channel for AMRINTEGRATOR
    // Join tool results (keyed by meta.id) with the filtered GFF
    
    // Add deeparg results (or empty list if skipped)
    if (!skip_deeparg) {
        ch_deeparg_keyed = ch_deeparg_results.map { meta, file -> [meta.id, file] }
        ch_for_amrintegrator = ch_gff_filtered
            .join(ch_deeparg_keyed, remainder: true)
            .map { id, meta, gff, deeparg -> 
                [id, meta, gff, deeparg ?: []]  // Replace null with empty list
            }
    } else {
        ch_for_amrintegrator = ch_gff_filtered
            .map { id, meta, gff -> [id, meta, gff, []] }  // Empty list for skipped tool
    }

    // Add rgi results (or empty list if skipped)
    if (!skip_rgi) {
        ch_rgi_keyed = ch_rgi_results.map { meta, file -> [meta.id, file] }
        ch_for_amrintegrator = ch_for_amrintegrator
            .join(ch_rgi_keyed, remainder: true)
            .map { id, meta, gff, deeparg, rgi -> 
                [id, meta, gff, deeparg, rgi ?: []]  // Replace null with empty list
            }
    } else {
        ch_for_amrintegrator = ch_for_amrintegrator
            .map { id, meta, gff, deeparg -> 
                [id, meta, gff, deeparg, []]  // Empty list for skipped tool
            }
    }

    // Add amrfinderplus results (or empty list if skipped)
    if (!skip_amrfinderplus) {
        ch_amrfinder_keyed = ch_amrfinderplus_results.map { meta, file -> [meta.id, file] }
        ch_for_amrintegrator = ch_for_amrintegrator
            .join(ch_amrfinder_keyed, remainder: true)
            .map { id, meta, gff, deeparg, rgi, amrfinder -> 
                [id, meta, gff, deeparg, rgi, amrfinder ?: []]  // Replace null with empty list
            }
    } else {
        ch_for_amrintegrator = ch_for_amrintegrator
            .map { id, meta, gff, deeparg, rgi -> 
                [id, meta, gff, deeparg, rgi, []]  // Empty list for skipped tool
            }
    }

    // Final reordering to match AMRINTEGRATOR input: [meta, deeparg, rgi, amrfinder, gff]
    // Remove the id key and reorder
    ch_for_amrintegrator = ch_for_amrintegrator
        .map { id, meta, gff, deeparg, rgi, amrfinder ->
            [meta, deeparg, rgi, amrfinder, gff]
        }

    AMRINTEGRATOR(ch_for_amrintegrator)
    ch_versions = ch_versions.mix(AMRINTEGRATOR.out.versions.first())

    emit:
    gff      = AMRINTEGRATOR.out.gff            // channel: [ val(meta), [ gff ] ]
    versions = ch_versions                      // channel: [ versions.yml ]
}
