// Subworkflow to generate BGCs annotation using Sanntis, Gecco and Antismash
// Outputs are integrated into a single GFF3 format output
include { GFF2GBK                               } from '../../../modules/ebi-metagenomics/mgnifypipelinestoolkit/gff2gbk/main'
include { ANTISMASH_ANTISMASH                   } from '../../../modules/nf-core/antismash/antismash/main'
include { ANTISMASH_ANTISMASHDOWNLOADDATABASES  } from '../../../modules/nf-core/antismash/antismashdownloaddatabases/main'
include { ANTISMASH_JSON2GFF                    } from '../../../modules/ebi-metagenomics/antismash/json2gff/main'
include { SANNTIS                               } from '../../../modules/ebi-metagenomics/sanntis/main'
include { GECCO_RUN                             } from '../../../modules/nf-core/gecco/run/main'
include { GECCO_CONVERT                         } from '../../../modules/nf-core/gecco/convert/main'
include { BGCSMAPPER                            } from '../../../modules/ebi-metagenomics/mgnifypipelinestoolkit/bgcsmapper/main'
include { INTERPROSCAN                          } from '../../../modules/ebi-metagenomics/interproscan/main'

workflow BGC_ANNOTATION {

    take:
    ch_inputs             // channel: tuple( val(meta), path(contigs), path(gff), path(proteins), path(ips_annot) )
    ch_antismash_db       // channel: path(antismash_db)
    ch_ips_db             // channel: path(interproscan_db)
    skip_sanntis          // boolean
    skip_gecco            // boolean
    skip_antismash        // boolean
    ch_user_sanntis_gff   // channel: tuple( val(meta), path(gff) ) — catalogue result; empty if no manifest
    ch_user_gecco_gff     // channel: tuple( val(meta), path(gff) ) — catalogue result; empty if no manifest
    ch_user_antismash_gff // channel: tuple( val(meta), path(gff) ) — catalogue result; empty if no manifest
    manifest_provided     // boolean: true when --annotation_manifest is set

    main:
    ch_versions = channel.empty()

    // Extract individual inputs from input channel
    def inputs = ch_inputs.multiMap { meta, contigs, gff, proteins, ips_annot ->
        gff:         tuple(meta, gff)
        togbk_input: tuple(meta, contigs, gff, proteins)
        ips:         tuple(meta, ips_annot)
        prots:       tuple(meta, proteins)
    }

    // IPS channel emitted for downstream use (PATHOFACT2 / COMBINEREPORTER).
    // Populated inside the sanntis block when running internally; stays empty when
    // manifest is active (manifest IPS is routed in the main workflow, not here).
    ch_ips_out = channel.empty()

    // Declare result channels; populated below based on mode.
    ch_sanntis_results   = channel.empty()
    ch_gecco_results     = channel.empty()
    ch_antismash_results = channel.empty()

    if (manifest_provided) {
        // Pre-computed GFFs from the annotation manifest are used directly.
        // GFF2GBK, INTERPROSCAN, and all three BGC tool runs are skipped.
        ch_sanntis_results   = ch_user_sanntis_gff
        ch_gecco_results     = ch_user_gecco_gff
        ch_antismash_results = ch_user_antismash_gff
    } else {
        // Run tools internally.  All GFF2GBK-dependent code is nested here so that
        // GFF2GBK.out.gbk is never referenced when manifest_provided is true.
        GFF2GBK(inputs.togbk_input)

        // Run SanntiS
        if (!skip_sanntis) {

            // Treat [] (no IPS provided) the same as null: both mean IPS is absent.
            ch_ips_branched = inputs.ips.branch { meta, ips_annot ->
                provided: ips_annot
                missing:  !ips_annot
            }

            // Samples missing IPS: filter proteins to only those that need InterProScan.
            ch_prots_missing = inputs.prots
                .join(
                    ch_ips_branched.missing.map { meta, _ips -> tuple(meta, true) },
                    remainder: true
                )
                .filter { meta, _prots, needs_ips -> needs_ips == true }
                .map { meta, prots, _flag -> tuple(meta, prots) }

            INTERPROSCAN(ch_prots_missing, ch_ips_db)
            ch_versions = ch_versions.mix(INTERPROSCAN.out.versions)

            ips_tsv = ch_ips_branched.provided
                .mix(
                    INTERPROSCAN.out.tsv
                )
            ch_ips_out = ips_tsv

            ch_sanntis_input = ips_tsv
                .join(GFF2GBK.out.gbk, by: 0, remainder: true)
                .filter { meta, ips, gbk -> ips != null && gbk != null }
                .map { meta, ips, gbk -> tuple(meta, ips, gbk, []) }

            SANNTIS(ch_sanntis_input)
            ch_versions = ch_versions.mix(SANNTIS.out.versions)
            ch_sanntis_results = SANNTIS.out.gff
        }

        // Run GECCO
        if (!skip_gecco) {
            ch_gecco_input = GFF2GBK.out.gbk.map { meta, gbk -> tuple(meta, gbk, []) }
            GECCO_RUN(ch_gecco_input, [])

            ch_gecco_convert_input = GECCO_RUN.out.clusters
                .join(GECCO_RUN.out.gbk, by: 0, remainder: true)
                .filter { meta, clusters, gbk -> clusters != null && gbk != null }
                .map { meta, clusters, gbk -> tuple(meta, clusters, gbk) }

            GECCO_CONVERT(ch_gecco_convert_input, "clusters", "gff")
            ch_versions = ch_versions.mix(GECCO_CONVERT.out.versions)
            ch_gecco_results = GECCO_CONVERT.out.gff
        }

        // Run antiSMASH
        if (!skip_antismash && ch_antismash_db) {
            antismash_db = ch_antismash_db
        } else if (!skip_antismash && !ch_antismash_db) {
            ANTISMASH_ANTISMASHDOWNLOADDATABASES()
            ch_versions = ch_versions.mix(ANTISMASH_ANTISMASHDOWNLOADDATABASES.out.versions)
            antismash_db = ANTISMASH_ANTISMASHDOWNLOADDATABASES.out.database
        }

        if (!skip_antismash) {
            ANTISMASH_ANTISMASH(GFF2GBK.out.gbk, antismash_db, [])
            ANTISMASH_JSON2GFF(ANTISMASH_ANTISMASH.out.json_results)
            ch_antismash_results = ANTISMASH_JSON2GFF.out.gff
        }
    }

    // multiMap each result channel to broadcast to both consumers:
    //   flags  → ch_has_bgc_results (determines which samples reach BGCSMAPPER)
    //   bgcsmapper → ch_for_bgcsmapper join chain
    // A plain queue channel used in two operator chains would split items between them.
    ch_sanntis_split   = ch_sanntis_results.multiMap   { meta, gff -> flags: tuple(meta, gff); bgcsmapper: tuple(meta, gff) }
    ch_gecco_split     = ch_gecco_results.multiMap     { meta, gff -> flags: tuple(meta, gff); bgcsmapper: tuple(meta, gff) }
    ch_antismash_split = ch_antismash_results.multiMap { meta, gff -> flags: tuple(meta, gff); bgcsmapper: tuple(meta, gff) }

    // Per-tool flags: tuple(meta, true).
    // When manifest_provided, results come from user channels regardless of skip flags.
    ch_sanntis_flags = (skip_sanntis && !manifest_provided)
        ? channel.empty()
        : ch_sanntis_split.flags.map { meta, _file -> tuple(meta, true) }

    ch_gecco_flags = (skip_gecco && !manifest_provided)
        ? channel.empty()
        : ch_gecco_split.flags.map { meta, _file -> tuple(meta, true) }

    ch_antismash_flags = (skip_antismash && !manifest_provided)
        ? channel.empty()
        : ch_antismash_split.flags.map { meta, _file -> tuple(meta, true) }

    // Combine flags -> per-sample has_results boolean
    ch_has_bgc_results = ch_sanntis_flags
        .mix(ch_gecco_flags)
        .mix(ch_antismash_flags)
        .groupTuple()
        .map { meta, _vals -> tuple(meta, true) }

    // Filter base GFF to samples with >=1 tool output
    ch_gff_filtered = inputs.gff
        .join(ch_has_bgc_results, remainder: true)
        .filter { meta, gff, has_results -> has_results == true }
        .map { meta, gff, _has_results -> tuple(meta, gff) }

    /*
     * Build BGCSMAPPER input
     * Strategy: start from base GFF, join each optional tool result by meta key,
     * then apply placeholders in a single final map.
     * Expected BGCSMAPPER input:
     *   tuple(meta, base_gff, sanntis_gff, gecco_gff, antismash_gff)
     */

    ch_for_bgcsmapper = ch_gff_filtered

    if (manifest_provided || !skip_sanntis) {
        ch_for_bgcsmapper = ch_for_bgcsmapper
            .join(ch_sanntis_split.bgcsmapper, by: 0, remainder: true)
    } else {
        ch_for_bgcsmapper = ch_for_bgcsmapper
            .map { meta, gff -> tuple(meta, gff, null) }
    }

    if (manifest_provided || !skip_gecco) {
        ch_for_bgcsmapper = ch_for_bgcsmapper
            .join(ch_gecco_split.bgcsmapper, by: 0, remainder: true)
    } else {
        ch_for_bgcsmapper = ch_for_bgcsmapper
            .map { meta, gff, sanntis -> tuple(meta, gff, sanntis, null) }
    }

    if (manifest_provided || !skip_antismash) {
        ch_for_bgcsmapper = ch_for_bgcsmapper
            .join(ch_antismash_split.bgcsmapper, by: 0, remainder: true)
    } else {
        ch_for_bgcsmapper = ch_for_bgcsmapper
            .map { meta, gff, sanntis, gecco -> tuple(meta, gff, sanntis, gecco, null) }
    }

    // Single final map: apply [] placeholders for any null (skipped/missing) tool
    ch_for_bgcsmapper = ch_for_bgcsmapper
        .map { meta, gff, sanntis_gff, gecco_gff, antismash_gff ->
            tuple(
                meta,
                gff,
                sanntis_gff  ?: [],
                gecco_gff    ?: [],
                antismash_gff ?: []
            )
        }

    BGCSMAPPER(ch_for_bgcsmapper)

    ch_bgc_output = BGCSMAPPER.out.gff
        .join(BGCSMAPPER.out.json, by: 0, remainder: true)
        .map { meta, gff, json -> tuple(meta, gff, json ?: []) }

    emit:
    bgc_output = ch_bgc_output   // channel: [ val(meta), gff, json ]
    ips_tsv    = ch_ips_out      // channel: [ val(meta), tsv ] — provided + InterProScan-generated; empty if skip_sanntis
    versions   = ch_versions
}
