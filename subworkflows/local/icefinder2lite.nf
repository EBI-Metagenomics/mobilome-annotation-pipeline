include { HMMSCAN                       } from '../../modules/local/hmmscan'
include { PRESCAN_TO_FASTA              } from '../../modules/local/prescan_to_fasta'
include { BLAST_BLASTP as BLASTP_PROKKA } from '../../modules/nf-core/blast/blastp'
include { PROCESS_BLASTP_PROKKA         } from '../../modules/local/process_blastp_prokka'
include { MACSYFINDER                   } from '../../modules/local/macsyfinder'
include { VMATCH                        } from '../../modules/local/vmatch'
include { REFINE_BOUNDARIES             } from '../../modules/local/ice_refine_boundaries'

workflow ICEFINDER2_LITE {
    take:
    ch_icf2_inputs                  // channel: tuple( val(meta), path(assembly_5kb), path(proteins_faa), path(proteins_gff))
    ch_icefinder_hmm_models         // channel: tuple( val(db_id), path(ice_hmm_models_files) )
    ch_icefinder_macsyfinder_models // channel: file(ice_macsy_models)
    ch_icefinder_prokka_uniprot_db  // channel: tuple( val(db_id), path(prokka_uniprot_db) )

    main:
    ch_versions = Channel.empty()

    // multiMap is required here to broadcast all samples to all outputs;
    // multiple .map{} calls on the same queue channel would split items between operators.
    def ch_icf2_split = ch_icf2_inputs.multiMap { meta, contigs, faa, gff ->
        faa:         tuple(meta, faa)                   // for HMMSCAN
        for_prescan: tuple(meta, contigs, faa, gff)     // for prescan join (all fields at once)
        assem:       tuple(meta, contigs)               // for candidate filtering
        gff:         tuple(meta, gff)                   // for candidate filtering
    }

    // Prescanning for candidate contigs
    HMMSCAN(ch_icf2_split.faa, ch_icefinder_hmm_models)
    ch_versions = ch_versions.mix(HMMSCAN.out.versions)

    // Join HMMSCAN output with all input fields in one step to avoid reusing ch_faa
    prescan_input_ch = HMMSCAN.out.hmmscan_tbl
        .join(ch_icf2_split.for_prescan)
        .map { meta, hmmscan_tbl, contigs, faa, gff ->
            tuple(meta, hmmscan_tbl, faa, contigs, gff)
        }

    PRESCAN_TO_FASTA(prescan_input_ch, params.prescan_evalue_threshold)
    ch_versions = ch_versions.mix(PRESCAN_TO_FASTA.out.versions)


    // Filter for samples with candidates (non-empty files)
    ch_candidates = PRESCAN_TO_FASTA.out.candidates_faa
       .join(PRESCAN_TO_FASTA.out.candidates_fna)
        .filter { _meta, faa, fna ->
            faa != null && fna != null
        }

    // Split candidates once for all downstream uses;
    // multiple .map{} calls on the same queue channel would split items between operators.
    def ch_cand_split = ch_candidates.multiMap { meta, faa, fna ->
        for_metas_assem:  tuple(meta, true)  // for assembly_filtered join
        for_metas_gff:    tuple(meta, true)  // for gff_filtered join
        for_macsyfinder:  tuple(meta, faa)
        for_blastp:       tuple(meta, faa)
        for_vmatch:       tuple(meta, fna)
    }

    /*
     * Filter down only the samples with candidates proceed to downstream analysis
     */
    ch_assembly_filtered = ch_icf2_split.assem
        .join(ch_cand_split.for_metas_assem, remainder: true)
        .filter { _meta, _assembly, candidate_flag -> {
                candidate_flag == true
            }
        }
        .map { meta, assembly, _candidate_flag -> {
                [meta, assembly]
            }
        }

    ch_gff_filtered = ch_icf2_split.gff
        .join(ch_cand_split.for_metas_gff, remainder: true)
        .filter { _meta, _gff, candidate_flag -> {
                candidate_flag == true
            }
        }
        .map { meta, gff, _candidate_flag -> {
                [meta, gff]
            }
        }

    // Run downstream processes only on samples with candidates
    MACSYFINDER(
        ch_cand_split.for_macsyfinder,
        ch_icefinder_macsyfinder_models,
    )
    ch_versions = ch_versions.mix(MACSYFINDER.out.versions)

    BLASTP_PROKKA(
        ch_cand_split.for_blastp,
        ch_icefinder_prokka_uniprot_db,
        'tsv',
    )
    ch_versions = ch_versions.mix(BLASTP_PROKKA.out.versions)

    PROCESS_BLASTP_PROKKA(BLASTP_PROKKA.out.tsv)
    ch_versions = ch_versions.mix(PROCESS_BLASTP_PROKKA.out.versions)

    VMATCH(ch_cand_split.for_vmatch)
    ch_versions = ch_versions.mix(VMATCH.out.versions)

    // REFINE_BOUNDARIES - now only runs on samples that have candidates
    // All channels will have matching samples since we filtered ch_assembly
    def input_for_refine = ch_assembly_filtered
        .join(
            ch_gff_filtered
        ).join(
            MACSYFINDER.out.macsyfinder_tsv
        ).join(
            VMATCH.out.vmatch_tsv
        ).join(
            PROCESS_BLASTP_PROKKA.out.uniprot_product_names, remainder: true
        )

    REFINE_BOUNDARIES(
        // This shenanigans is to account for optional uniprot_product_names, moving to types in 25.10 should clean this
        input_for_refine.map {
            meta, assembly, gff_filtered, macsyfinder_tsv, vmatch_tsv, uniprot_product_names ->
            return [
                meta,
                assembly,
                gff_filtered,
                macsyfinder_tsv,
                vmatch_tsv,
                uniprot_product_names ?: []
            ]
        }
    )
    ch_versions = ch_versions.mix(REFINE_BOUNDARIES.out.versions)

    emit:
    ices_tsv = REFINE_BOUNDARIES.out.ices_tsv
    versions = ch_versions
}
