include { MOBILEOG_BLAST         } from '../modules/mobileog_blast'
include { RESFINDER_BLAST        } from '../modules/resfinder_blast'
include { VFDB_BLAST             } from '../modules/vfdb_blast'
include { ISFINDER_BLAST         } from '../modules/isfinder_blast'
include { METAL_RESISTANCE_BLAST } from '../modules/metal_resistance_blast'
include { DEGRADATION_BLAST      } from '../modules/degradation_blast'
include { SYMBIOSIS_BLAST        } from '../modules/symbiosis_blast'
include { INTEGRATE_ANNOTATIONS  } from '../modules/integrate_annotations'

workflow BLAST_ANNOTATIONS {
    take:
    ch_fasta      // channel: [meta, fasta_file]
    ch_faa        // channel: [meta, faa_file]
    ch_genes      // channel: [meta, gene_file]
    mobileog_db   // path: mobileog database
    resfinder_db  // path: resfinder database
    vfdb_db       // path: vfdb database
    isfinder_db   // path: isfinder database
    metal_db      // path: metal resistance database
    degradation_db // path: degradation database
    symbiosis_db  // path: symbiosis database

    main:
    ch_versions = Channel.empty()

    // Run all BLAST searches
    MOBILEOG_BLAST(ch_faa, mobileog_db)
    RESFINDER_BLAST(ch_fasta, resfinder_db)
    VFDB_BLAST(ch_faa, vfdb_db)
    ISFINDER_BLAST(ch_faa, isfinder_db)
    METAL_RESISTANCE_BLAST(ch_faa, metal_db)
    DEGRADATION_BLAST(ch_faa, degradation_db)
    SYMBIOSIS_BLAST(ch_fasta, symbiosis_db)

    // Collect versions
    ch_versions = ch_versions.mix(MOBILEOG_BLAST.out.versions)
    ch_versions = ch_versions.mix(RESFINDER_BLAST.out.versions)
    ch_versions = ch_versions.mix(VFDB_BLAST.out.versions)
    ch_versions = ch_versions.mix(ISFINDER_BLAST.out.versions)
    ch_versions = ch_versions.mix(METAL_RESISTANCE_BLAST.out.versions)
    ch_versions = ch_versions.mix(DEGRADATION_BLAST.out.versions)
    ch_versions = ch_versions.mix(SYMBIOSIS_BLAST.out.versions)

    // Combine all annotation results by meta
    ch_combined = ch_genes
        .join(MOBILEOG_BLAST.out.filtered_results, by: 0, remainder: true)
        .join(RESFINDER_BLAST.out.filtered_results, by: 0, remainder: true)
        .join(VFDB_BLAST.out.filtered_results, by: 0, remainder: true)
        .join(ISFINDER_BLAST.out.filtered_results, by: 0, remainder: true)
        .join(METAL_RESISTANCE_BLAST.out.filtered_results, by: 0, remainder: true)
        .join(DEGRADATION_BLAST.out.filtered_results, by: 0, remainder: true)
        .join(SYMBIOSIS_BLAST.out.filtered_results, by: 0, remainder: true)

    // Integrate all annotations
    INTEGRATE_ANNOTATIONS(ch_combined)
    ch_versions = ch_versions.mix(INTEGRATE_ANNOTATIONS.out.versions)

    emit:
    annotated_genes = INTEGRATE_ANNOTATIONS.out.annotated_genes
    versions = ch_versions
}
