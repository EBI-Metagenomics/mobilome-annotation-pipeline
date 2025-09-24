include { PRODIGAL                      } from '../modules/nf-core/prodigal'
include { BLAST_BLASTP as BLASTP_PROKKA } from '../modules/nf-core/blast/blastp'

include { HMMSCAN                       } from '../modules/local/hmmscan'
include { PRESCAN_TO_FASTA              } from '../modules/local/prescan_to_fasta'
include { PROCESS_BLASTP_PROKKA         } from '../modules/local/process_blastp_prokka'
include { ARAGORN                       } from '../modules/local/aragorn'
include { TRNAS_INTEGRATOR              } from '../modules/local/trnas_integrator'
include { MACSYFINDER                   } from '../modules/local/macsyfinder'
include { VMATCH                        } from '../modules/local/vmatch'
include { REFINE_BOUNDARIES             } from '../modules/local/ice_refine_boundaries'

workflow ICEFINDER2_LITE {
    take:
    ch_assembly                      // channel: tuple( val(meta), path(assembly_5kb) )
    ch_icefinder_hmm_models         // channel: tuple( val(db_id), path(ice_hmm_models_files) )
    ch_icefinder_macsyfinder_models // channel: file(ice_macsy_models)
    ch_icefinder_prokka_uniprot_db  // channel: tuple( val(db_id), path(prokka_uniprot_db) )

    main:
    ch_versions = Channel.empty()

    // Prescanning for candidate contigs
    PRODIGAL(ch_assembly, 'gff')
    ch_versions = ch_versions.mix(PRODIGAL.out.versions)

    HMMSCAN(PRODIGAL.out.amino_acid_fasta, ch_icefinder_hmm_models)
    ch_versions = ch_versions.mix(HMMSCAN.out.versions)

    prescan_input_ch = HMMSCAN.out.hmmscan_tbl
        .join(PRODIGAL.out.amino_acid_fasta)
        .join(ch_assembly)

    PRESCAN_TO_FASTA(prescan_input_ch, params.prescan_evalue_threshold)
    ch_versions = ch_versions.mix(PRESCAN_TO_FASTA.out.versions)

    // Filter for samples with candidates (non-empty files)
    ch_candidates = PRESCAN_TO_FASTA.out.candidates_faa
        .join(PRESCAN_TO_FASTA.out.candidates_fna)
        .filter { _meta, faa, fna ->
            faa != null && fna != null
        }

    // Extract just the meta information from candidates for filtering
    ch_candidate_metas = ch_candidates.map { meta, _faa, _fna -> [meta, true] }

    /*
     * Filter down only the samples with candidates proceed to downstream analysis
     */
    ch_assembly_filtered = ch_assembly
        .join(ch_candidate_metas, remainder: true)
        .filter { _meta, _assembly, candidate_flag -> {
                candidate_flag == true
            }
        }
        .map { meta, assembly, _candidate_flag -> {
                [meta, assembly]
            }
        }

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

    ARAGORN(
        ch_candidates.map { meta, _faa, fna -> tuple(meta, fna) }
    )
    ch_versions = ch_versions.mix(ARAGORN.out.versions)

    /*
    * Filter down the PROKKA samples with those thare are relevant
    */
    ch_prokka_filtered = PRODIGAL.out.gene_annotations
        .join(ch_candidate_metas, remainder: true)
        .filter { _meta, _prokka, candidate_flag -> {
                candidate_flag == true
            }
        }
        .map { meta, prokka, _candidate_flag -> {
                [meta, prokka]
            }
        }

    TRNAS_INTEGRATOR(
        ARAGORN.out.rnas_tbl.join( ch_prokka_filtered )
    )
    ch_versions = ch_versions.mix(TRNAS_INTEGRATOR.out.versions)

    VMATCH(
        ch_candidates.map { meta, _faa, fna -> tuple(meta, fna) }
    )
    ch_versions = ch_versions.mix(VMATCH.out.versions)

    // REFINE_BOUNDARIES - now only runs on samples that have candidates
    // All channels will have matching samples since we filtered ch_assembly
    REFINE_BOUNDARIES(
        ch_assembly_filtered.join(
            TRNAS_INTEGRATOR.out.merged_gff
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
