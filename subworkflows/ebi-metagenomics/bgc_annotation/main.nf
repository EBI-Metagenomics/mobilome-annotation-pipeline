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
    ch_inputs           // channel: tuple( val(meta), path(contigs), path(gff), path(proteins), path(ips_annot) )
    ch_antismash_db     // channel: path(antismash_db)
    ch_ips_db           // channel: path(interproscan_db)
    skip_sanntis        // boolean
    skip_gecco          // boolean
    skip_antismash      // boolean

    main:
    ch_versions = channel.empty()

    // Extract individual inputs from input channel

    def inputs = ch_inputs.multiMap { meta, contigs, gff, proteins, ips_annot ->
        gff:         tuple(meta, gff)
        togbk_input: tuple(meta, contigs, gff, proteins)
        ips:         tuple(meta, ips_annot)
        prots:       tuple(meta, proteins)
    }

    // Transform GFF into GBK
    GFF2GBK(inputs.togbk_input)

    // Run SanntiS
    if (!skip_sanntis) {

        // Split IPS input into provided vs missing
        ch_ips_branched = inputs.ips.branch { meta, ips_annot ->
            provided: ips_annot != null
            missing:  ips_annot == null
        }

        // Samples missing IPS: filter proteins to only those that need InterProScan.
        // Use remainder:true + filter to avoid a strict-mode join mismatch when all
        // samples already have IPS (missing branch is empty, inputs.prots is not).
        ch_prots_missing = inputs.prots
            .join(
                ch_ips_branched.missing.map { meta, _ips -> tuple(meta, true) },
                remainder: true
            )
            .filter { meta, _prots, needs_ips -> needs_ips == true }
            .map { meta, prots, _flag -> tuple(meta, prots) }

        /*
         * No automatic DB download anymore.
         * If IPS is missing and INTERPROSCAN is needed, ch_ips_db must be provided.
         * If the database is needed and not provided, the process will naturally fail
         */
        INTERPROSCAN(ch_prots_missing, ch_ips_db)
        ch_versions = ch_versions.mix(INTERPROSCAN.out.versions)

        /*
         * Build ips_tsv for all samples:
         *   - provided IPS for samples that had it
         *   - newly generated IPS for missing samples
         */
        ips_tsv = ch_ips_branched.provided.mix(INTERPROSCAN.out.tsv)

        // Run SanntiS
        ch_sanntis_input = ips_tsv
            .join(GFF2GBK.out.gbk, by: 0)
            .map { meta, ips, gbk ->
                tuple(meta, ips, gbk, [])
            }

        SANNTIS(ch_sanntis_input)
        ch_versions = ch_versions.mix(SANNTIS.out.versions)
        ch_sanntis_results = SANNTIS.out.gff
    }

    // Run GECCO
    if (!skip_gecco) {
        ch_gecco_input =  GFF2GBK.out.gbk.map { meta, gbk -> tuple(meta, gbk, []) }
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

    // Ensure each optional result channel exists even if tool is skipped
    ch_sanntis_results   = skip_sanntis   ? channel.empty() : ch_sanntis_results
    ch_gecco_results     = skip_gecco     ? channel.empty() : ch_gecco_results
    ch_antismash_results = skip_antismash ? channel.empty() : ch_antismash_results

    // Per-tool flags: tuple(meta, true)
    ch_sanntis_flags = skip_sanntis
        ? channel.empty()
        : ch_sanntis_results.map { meta, _file -> tuple(meta, true) }

    ch_gecco_flags = skip_gecco
        ? channel.empty()
        : ch_gecco_results.map { meta, _file -> tuple(meta, true) }

    ch_antismash_flags = skip_antismash
        ? channel.empty()
        : ch_antismash_results.map { meta, _file -> tuple(meta, true) }

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

    if (!skip_sanntis) {
        ch_for_bgcsmapper = ch_for_bgcsmapper
            .join(ch_sanntis_results, by: 0, remainder: true)
    } else {
        ch_for_bgcsmapper = ch_for_bgcsmapper
            .map { meta, gff -> tuple(meta, gff, null) }
    }

    if (!skip_gecco) {
        ch_for_bgcsmapper = ch_for_bgcsmapper
            .join(ch_gecco_results, by: 0, remainder: true)
    } else {
        ch_for_bgcsmapper = ch_for_bgcsmapper
            .map { meta, gff, sanntis -> tuple(meta, gff, sanntis, null) }
    }

    if (!skip_antismash) {
        ch_for_bgcsmapper = ch_for_bgcsmapper
            .join(ch_antismash_results, by: 0, remainder: true)
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
    versions   = ch_versions
}
