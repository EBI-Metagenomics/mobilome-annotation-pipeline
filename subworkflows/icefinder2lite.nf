include { PRODIGAL                      } from '../modules/nf-core/prodigal'
include { BLAST_BLASTP as BLASTP_PROKKA } from '../modules/nf-core/blast/blastp'

include { HMMSCAN as PRESCAN            } from '../modules/local/hmmscan'
include { PRESCAN_TO_FASTA              } from '../modules/local/prescan_to_fasta'
include { PROCESS_BLASTP_PROKKA         } from '../modules/local/process_blastp_prokka'
include { ARAGORN                       } from '../modules/local/aragorn'
include { TRNAS_INTEGRATOR              } from '../modules/local/trnas_integrator'
include { MACSYFINDER                   } from '../modules/local/macsyfinder'
include { VMATCH                        } from '../modules/local/vmatch'
include { REFINE_BOUNDARIES             } from '../modules/local/ice_refine_boundaries'

workflow ICEFINDER2_LITE {
    take:
    ch_assembly          // channel: tuple( val(meta), path(assembly_5kb) )
    ch_ice_hmm_models    // channel: tuple( val(db_id), path(ice_hmm_models_files) )
    ch_ice_macsy_models  // channel: path(ice_macsy_models)
    ch_prokka_uniprot_db // channel: tuple( val(db_id), path(prokka_uniprot_db) )

    main:
    ch_versions = Channel.empty()

    // Prescanning for candidate contigs
    PRODIGAL(ch_assembly, 'gff')

    PRESCAN(PRODIGAL.out.amino_acid_fasta, ch_ice_hmm_models)

    prescan_input_ch = PRESCAN.out.hmmscan_tbl
        .join(PRODIGAL.out.amino_acid_fasta)
        .join(ch_assembly)

    PRESCAN_TO_FASTA(prescan_input_ch)

    // Evaluating if candidates were found
    ch_candidates = PRESCAN_TO_FASTA.out.candidates_faa
        .join(PRESCAN_TO_FASTA.out.candidates_fna)
        .filter { _meta, faa, fna ->
            faa.exists() && fna.exists()
        }

    // Only run downstream processes if candidates exist
    if (ch_candidates) {
        MACSYFINDER(ch_candidates.map { meta, faa, _fna -> tuple(meta, faa) }, ch_ice_macsy_models)
        BLASTP_PROKKA(ch_candidates.map { meta, faa, _fna -> tuple(meta, faa) }, ch_prokka_uniprot_db, 'tsv')
        PROCESS_BLASTP_PROKKA(BLASTP_PROKKA.out.tsv)
        ARAGORN(ch_candidates.map { meta, _faa, fna -> tuple(meta, fna) })
        TRNAS_INTEGRATOR(ARAGORN.out.rnas_tbl.join(PRODIGAL.out.gene_annotations))
        VMATCH(ch_candidates.map { meta, _faa, fna -> tuple(meta, fna) })

        // Collect versions
        ch_versions = ch_versions.mix(MACSYFINDER.out.versions)
        ch_versions = ch_versions.mix(BLASTP_PROKKA.out.versions)
        ch_versions = ch_versions.mix(ARAGORN.out.versions)
        ch_versions = ch_versions.mix(VMATCH.out.versions)

        // REFINE_BOUNDARIES - this will only emit results for samples with ICEs
        REFINE_BOUNDARIES(
            ch_assembly.join(
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
        final_ices_tsv = REFINE_BOUNDARIES.out.ices_tsv
    }
    else {
        final_ices_tsv = Channel.empty()
    }

    emit:
    ices_tsv = final_ices_tsv // This will be empty for samples without ICEs
    versions = ch_versions
}
