/* NF-CORE */
include { TABIX_BGZIP as TABIX_BGZIP_CLEAN           } from '../modules/nf-core/tabix/bgzip/main'
include { TABIX_BGZIP as TABIX_BGZIP_EXTRA           } from '../modules/nf-core/tabix/bgzip/main'
include { TABIX_BGZIP as TABIX_BGZIP_FULL            } from '../modules/nf-core/tabix/bgzip/main'
include { TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_CLEAN } from '../modules/nf-core/tabix/bgziptabix/main'
include { TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_EXTRA } from '../modules/nf-core/tabix/bgziptabix/main'
include { TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_FULL  } from '../modules/nf-core/tabix/bgziptabix/main'

/* LOCAL */
include { GFF_MAPPING } from '../modules/local/gff_mapping.nf'


workflow GFF_MAPPING_COMPRESSION_AND_INDEXING {
    take:
    ch_gff_inputs // tuple(meta, mobilome_gff, user_gff)

    main:
    /*
     * This subworkflow generates GFF files using the gff_mapping module, then produces
     * compressed and indexed versions (.gz and .csi). These files enable rendering GFFs
     * in genome browsers such as JBrowse 2.
     * 
     * The workflow applies bgzip compression and tabix indexing:
     * - .gz.gzi index: enables efficient random access to bgzipped GFF files
     * - .gz.csi index: genomic region index optimized for large contigs and metagenomes
    */

    ch_versions = Channel.empty()

    // Run GFF_MAPPING to create the three uncompressed GFF files
    GFF_MAPPING(ch_gff_inputs)
    ch_versions = ch_versions.mix(GFF_MAPPING.out.versions)

    /**************/
    /* GZI index */
    /*************/
    TABIX_BGZIP_CLEAN(
        GFF_MAPPING.out.mobilome_clean_gff
    )
    ch_versions = ch_versions.mix(TABIX_BGZIP_CLEAN.out.versions)

    TABIX_BGZIP_EXTRA(
        GFF_MAPPING.out.mobilome_extra_gff
    )
    ch_versions = ch_versions.mix(TABIX_BGZIP_EXTRA.out.versions)

    TABIX_BGZIP_FULL(
        GFF_MAPPING.out.mobilome_full_gff
    )
    ch_versions = ch_versions.mix(TABIX_BGZIP_FULL.out.versions)

    /**************/
    /* CSI index */
    /*************/
    TABIX_BGZIPTABIX_CLEAN(
        GFF_MAPPING.out.mobilome_clean_gff
    )
    ch_versions = ch_versions.mix(TABIX_BGZIPTABIX_CLEAN.out.versions)

    TABIX_BGZIPTABIX_EXTRA(
        GFF_MAPPING.out.mobilome_extra_gff
    )
    ch_versions = ch_versions.mix(TABIX_BGZIPTABIX_EXTRA.out.versions)

    TABIX_BGZIPTABIX_FULL(
        GFF_MAPPING.out.mobilome_full_gff
    )
    ch_versions = ch_versions.mix(TABIX_BGZIPTABIX_FULL.out.versions)

    emit:
    versions = ch_versions
}
