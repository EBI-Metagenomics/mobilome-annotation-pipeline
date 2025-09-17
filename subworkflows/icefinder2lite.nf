include { HMMSCAN                       } from '../modules/local/hmmscan'
include { PRESCAN_TO_FASTA              } from '../modules/local/prescan_to_fasta'
include { BLAST_BLASTP as BLASTP_PROKKA } from '../modules/nf-core/blast/blastp'
include { PROCESS_BLASTP_PROKKA         } from '../modules/local/process_blastp_prokka'
include { MACSYFINDER                   } from '../modules/local/macsyfinder'
include { VMATCH                        } from '../modules/local/vmatch'
include { REFINE_BOUNDARIES             } from '../modules/local/ice_refine_boundaries'

workflow ICEFINDER2_LITE {
    take:
    ch_icf2_inputs                  // channel: tuple( val(meta), path(assembly_5kb), path(proteins_faa), path(proteins_gff) )
    ch_icefinder_hmm_models         // channel: tuple( val(db_id), path(ice_hmm_models_files) )
    ch_icefinder_macsyfinder_models // channel: file(ice_macsy_models)
    ch_icefinder_prokka_uniprot_db  // channel: tuple( val(db_id), path(prokka_uniprot_db) )

    main:
    ch_versions = Channel.empty()

    ch_assem = ch_icf2_inputs.map{ meta, contigs, _faa, _gff -> tuple(meta, contigs) }
    ch_faa = ch_icf2_inputs.map{ meta, _contigs, faa, _gff -> tuple(meta, faa) }
    ch_gff = ch_icf2_inputs.map{ meta, _contigs, _faa, gff -> tuple(meta, gff) }

    // Prescanning for candidate contigs
    HMMSCAN( ch_faa, ch_ice_hmm_models)
    ch_versions = ch_versions.mix(HMMSCAN.out.versions)

    prescan_input_ch = HMMSCAN.out.hmmscan_tbl
        .join(ch_faa)
        .join(ch_assem)
        .join(ch_gff)

    PRESCAN_TO_FASTA(prescan_input_ch, params.prescan_evalue_threshold)
    ch_versions = ch_versions.mix(PRESCAN_TO_FASTA.out.versions)

    // Filter for samples with candidates (non-empty files)
    ch_candidates = PRESCAN_TO_FASTA.out.candidates_faa
        .join(PRESCAN_TO_FASTA.out.candidates_fna)
        .filter { _meta, faa, fna ->
            faa != null && fna != null
        }

    // Extract just the meta information from candidates for filtering
    ch_candidate_metas = ch_candidates.map { meta, _faa, _fna -> meta }

    /*
     * Filter down only the samples with candidates proceed to downstream analysis
     */
    ch_assembly_filtered = ch_assem
        .join(ch_candidate_metas, failOnMismatch: false)
    ch_gff_filtered = ch_gff
        .join(ch_candidate_metas, failOnMismatch: false)

    // Run downstream processes only on samples with candidates
    MACSYFINDER(
        ch_candidates.map { meta, faa, _fna -> tuple(meta, faa) },
        ch_icefinder_macsyfinder_models,
    )
    ch_versions = ch_versions.mix(MACSYFINDER.out.versions)

    BLASTP_PROKKA(
        ch_candidates.map { meta, faa, _fna -> tuple(meta, faa) },
        ch_icefinder_prokka_uniprot_db,
        'tsv',
    )
    ch_versions = ch_versions.mix(BLASTP_PROKKA.out.versions)

    PROCESS_BLASTP_PROKKA(BLASTP_PROKKA.out.tsv)
    ch_versions = ch_versions.mix(PROCESS_BLASTP_PROKKA.out.versions)

    VMATCH(
        ch_candidates.map { meta, _faa, fna -> tuple(meta, fna) }
    )
    ch_versions = ch_versions.mix(VMATCH.out.versions)

    // REFINE_BOUNDARIES - now only runs on samples that have candidates
    // All channels will have matching samples since we filtered ch_assembly
    REFINE_BOUNDARIES(ch_assembly_filtered
        .join(
            ch_gff_filtered
        ).join(
            MACSYFINDER.out.macsyfinder_tsv
        ).join(
            PROCESS_BLASTP_PROKKA.out.uniprot_product_names
        ).join(
            VMATCH.out.vmatch_tsv
        )
    )
    ch_versions = ch_versions.mix(REFINE_BOUNDARIES.out.versions)

    emit:
    ices_tsv = REFINE_BOUNDARIES.out.ices_tsv
    versions = ch_versions
}
