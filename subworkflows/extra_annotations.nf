include { BLAST_BLASTP as MOBILEOG_BLAST         } from '../modules/nf-core/blast/blastp'
include { BLAST_BLASTN as RESFINDER_BLAST        } from '../modules/nf-core/blast/blastn'
include { BLAST_BLASTP as VFDB_BLAST             } from '../modules/nf-core/blast/blastp'
include { BLAST_BLASTP as ISFINDER_BLAST         } from '../modules/nf-core/blast/blastp'
include { BLAST_BLASTP as METAL_RESISTANCE_BLAST } from '../modules/nf-core/blast/blastp'
include { BLAST_BLASTP as DEGRADATION_BLAST      } from '../modules/nf-core/blast/blastp'
include { BLAST_BLASTN as SYMBIOSIS_BLAST        } from '../modules/nf-core/blast/blastn'
include { FILTER_BLAST_RESULTS as FILTER_BLAST_RESULTS_RESFINDER   } from '../modules/local/functional_annotation/filter_blast_results'
include { FILTER_BLAST_RESULTS as FILTER_BLAST_RESULTS_VFDB        } from '../modules/local/functional_annotation/filter_blast_results'
include { FILTER_BLAST_RESULTS as FILTER_BLAST_RESULTS_ISFINDER    } from '../modules/local/functional_annotation/filter_blast_results'
include { FILTER_BLAST_RESULTS as FILTER_BLAST_RESULTS_METAL       } from '../modules/local/functional_annotation/filter_blast_results'
include { FILTER_BLAST_RESULTS as FILTER_BLAST_RESULTS_DEGRADATION } from '../modules/local/functional_annotation/filter_blast_results'
include { FILTER_BLAST_RESULTS as FILTER_BLAST_RESULTS_SYMBIOSIS   } from '../modules/local/functional_annotation/filter_blast_results'
include { INTEGRATE_ANNOTATIONS                  } from '../modules/local/integrate_annotations'

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

    // Create database channels with metadata
    ch_mobileog_db = Channel.value([id: 'mobileog']).combine(Channel.fromPath(mobileog_db))
    ch_resfinder_db = Channel.value([id: 'resfinder']).combine(Channel.fromPath(resfinder_db))
    ch_vfdb_db = Channel.value([id: 'virulence']).combine(Channel.fromPath(vfdb_db))
    ch_isfinder_db = Channel.value([id: 'isfinder']).combine(Channel.fromPath(isfinder_db))
    ch_metal_db = Channel.value([id: 'metal']).combine(Channel.fromPath(metal_db))
    ch_degradation_db = Channel.value([id: 'degradation']).combine(Channel.fromPath(degradation_db))
    ch_symbiosis_db = Channel.value([id: 'symbiosis']).combine(Channel.fromPath(symbiosis_db))

    // Run all BLAST searches using nf-core modules
    MOBILEOG_BLAST(ch_faa, ch_mobileog_db, 'tsv')
    RESFINDER_BLAST(ch_fasta, ch_resfinder_db)
    VFDB_BLAST(ch_faa, ch_vfdb_db, 'tsv')
    ISFINDER_BLAST(ch_faa, ch_isfinder_db, 'tsv')
    METAL_RESISTANCE_BLAST(ch_faa, ch_metal_db, 'tsv')
    DEGRADATION_BLAST(ch_faa, ch_degradation_db, 'tsv')
    SYMBIOSIS_BLAST(ch_fasta, ch_symbiosis_db)

    // Filter BLAST results where needed
    FILTER_BLAST_RESULTS_RESFINDER(RESFINDER_BLAST.out.txt)
    FILTER_BLAST_RESULTS_VFDB(VFDB_BLAST.out.tsv)
    FILTER_BLAST_RESULTS_ISFINDER(ISFINDER_BLAST.out.tsv)
    FILTER_BLAST_RESULTS_METAL(METAL_RESISTANCE_BLAST.out.tsv)
    FILTER_BLAST_RESULTS_DEGRADATION(DEGRADATION_BLAST.out.tsv)
    FILTER_BLAST_RESULTS_SYMBIOSIS(SYMBIOSIS_BLAST.out.txt)

    // Collect versions
    ch_versions = ch_versions.mix(MOBILEOG_BLAST.out.versions)
    ch_versions = ch_versions.mix(RESFINDER_BLAST.out.versions)
    ch_versions = ch_versions.mix(VFDB_BLAST.out.versions)
    ch_versions = ch_versions.mix(ISFINDER_BLAST.out.versions)
    ch_versions = ch_versions.mix(METAL_RESISTANCE_BLAST.out.versions)
    ch_versions = ch_versions.mix(DEGRADATION_BLAST.out.versions)
    ch_versions = ch_versions.mix(SYMBIOSIS_BLAST.out.versions)
    ch_versions = ch_versions.mix(FILTER_BLAST_RESULTS_RESFINDER.out.versions)
    ch_versions = ch_versions.mix(FILTER_BLAST_RESULTS_VFDB.out.versions)
    ch_versions = ch_versions.mix(FILTER_BLAST_RESULTS_ISFINDER.out.versions)
    ch_versions = ch_versions.mix(FILTER_BLAST_RESULTS_METAL.out.versions)
    ch_versions = ch_versions.mix(FILTER_BLAST_RESULTS_DEGRADATION.out.versions)
    ch_versions = ch_versions.mix(FILTER_BLAST_RESULTS_SYMBIOSIS.out.versions)

    // Combine all annotation results by meta (using raw output for MOBILEOG as it doesn't need filtering)
    ch_combined = ch_genes
        .join(MOBILEOG_BLAST.out.tsv, by: 0, remainder: true)
        .join(FILTER_BLAST_RESULTS_RESFINDER.out.filtered_results, by: 0, remainder: true)
        .join(FILTER_BLAST_RESULTS_VFDB.out.filtered_results, by: 0, remainder: true)
        .join(FILTER_BLAST_RESULTS_ISFINDER.out.filtered_results, by: 0, remainder: true)
        .join(FILTER_BLAST_RESULTS_METAL.out.filtered_results, by: 0, remainder: true)
        .join(FILTER_BLAST_RESULTS_DEGRADATION.out.filtered_results, by: 0, remainder: true)
        .join(FILTER_BLAST_RESULTS_SYMBIOSIS.out.filtered_results, by: 0, remainder: true)

    // Integrate all annotations
    INTEGRATE_ANNOTATIONS(ch_combined)
    ch_versions = ch_versions.mix(INTEGRATE_ANNOTATIONS.out.versions)

    emit:
    annotated_genes = INTEGRATE_ANNOTATIONS.out.annotated_genes
    versions = ch_versions
}
