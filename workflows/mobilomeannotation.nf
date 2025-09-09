include { validateParameters ; paramsHelp ; samplesheetToList } from 'plugin/nf-schema'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Inputs preparing modules
include { RENAME                          } from '../modules/local/rename_contigs'

// Annotation modules
include { PROKKA                          } from '../modules/local/prokka'
include { AMRFINDER_PLUS                  } from '../modules/local/amrfinder_plus'

// Mobile genetic elements prediction modules
include { INTEGRONFINDER                  } from '../modules/local/integronfinder'
include { ISESCAN                         } from '../modules/local/isescan'
include { GENOMAD                         } from '../modules/local/genomad'
include { VIRIFY_QC                       } from '../modules/local/virify_qc'

// Results integration and writing modules
include { AMRFINDER_REPORT                } from '../modules/local/amrfinder_report'
include { FASTA_WRITER                    } from '../modules/local/fasta_writer'
include { GFF_MAPPING                     } from '../modules/local/gff_mapping'
include { GFF_REDUCE                      } from '../modules/local/gff_reduce'
include { GT_GFF3VALIDATOR                } from '../modules/nf-core/gt/gff3validator/main'
include { INTEGRATOR                      } from '../modules/local/integrator'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { COMPOSITIONAL_OUTLIER_DETECTION } from '../subworkflows/compositional_outlier_detection'
include { ICEFINDER2_LITE                 } from '../subworkflows/icefinder2lite'
include { CUSTOM_DUMPSOFTWAREVERSIONS     } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { MULTIQC                         } from '../modules/nf-core/multiqc/main'

// TODO: add blast annotation workflow from mobilome proteins after integration

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MOBILOMEANNOTATION WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MOBILOMEANNOTATION {
    main:

    validateParameters()

    def ch_inputs = Channel.fromList(samplesheetToList(params.input, "./assets/schema_input.json"))
    ch_versions = Channel.empty()


    /*
    ******************************************************************************************************
    * The code below is transforming the input channels to handle optional inputs, such as the
    * user-provided GFF files for user_proteins and virify.
    * Nextflow doesn't handle optional inputs well, so we use a common hack to provide an empty
    * array ([]) when the input is missing.
    * For the user_proteins_gff, if the file is present, we emit a tuple with the metadata and the
    * file path. If the file is missing, we emit a tuple with the metadata and an empty array ([]).
    * Similarly, for the virify_gff, if the file is present, we emit a tuple with the metadata and
    * the file path. If the file is missing, we emit a tuple with the metadata and an empty array ([]).
    ******************************************************************************************************
    */

    def user_proteins_ch = ch_inputs.map { meta, _fasta, user_proteins_gff, _virify_gff ->
        {
            if (user_proteins_gff) {
                [meta, user_proteins_gff]
            }
            else {
                [meta, []]
            }
        }
    }

    // PREPROCESSING
    RENAME(ch_inputs.map { meta, fasta, _user_proteins_gff, _virify_gff -> [meta, fasta] })
    ch_versions = ch_versions.mix(RENAME.out.versions)

    PROKKA(RENAME.out.contigs_1kb)
    ch_versions = ch_versions.mix(PROKKA.out.versions)

    // Parsing VIRify gff file when an input is provided
    def user_virify_gff_ch = ch_inputs
        .map { meta, _fasta, _user_proteins_gff, virify_gff ->
            {
                [meta, virify_gff]
            }
        }
        .filter { _meta, virify_gff -> virify_gff != [] }

    VIRIFY_QC(user_virify_gff_ch)
    ch_versions = ch_versions.mix(VIRIFY_QC.out.versions)

    // PREDICTION
    // Collecting ICEfinder2 databases
    db_icefinder_hmm_models = Channel.fromPath("${params.icefinder_hmm_models}.*", checkIfExists: true)
        .collect()
        .map { ice_db_files ->
            [[id: file(params.icefinder_hmm_models).name], ice_db_files]
        }

    db_icefinder_macsyfinder_models = file(params.icefinder_macsyfinder_models, checkIfExists: true)

    db_icefinder_prokka_uniprot = Channel.fromPath("${params.icefinder_prokka_uniprot_db}/*", checkIfExists: true)
        .collect()
        .map { uniprot_db_files ->
            [[id: "prokka_uniprot"], uniprot_db_files]
        }

    ICEFINDER2_LITE(
        RENAME.out.contigs_5kb,
        db_icefinder_hmm_models,
        db_icefinder_macsyfinder_models,
        db_icefinder_prokka_uniprot,
    )
    ch_versions = ch_versions.mix(ICEFINDER2_LITE.out.versions)

    GENOMAD(RENAME.out.contigs_5kb)
    ch_versions = ch_versions.mix(GENOMAD.out.versions)

    INTEGRONFINDER(RENAME.out.contigs_5kb)
    ch_versions = ch_versions.mix(INTEGRONFINDER.out.versions)

    ISESCAN(RENAME.out.contigs_1kb)
    ch_versions = ch_versions.mix(ISESCAN.out.versions)

    COMPOSITIONAL_OUTLIER_DETECTION(RENAME.out.contigs_100kb, params.outlier_score_threshold)
    ch_versions = ch_versions.mix(COMPOSITIONAL_OUTLIER_DETECTION.out.versions)

    /**********************************************************************************************
    * The INTEGRATOR step takes a bunch of outputs from the previous steps.
    * The following code is re-shaping the input to accommodate
    * optional inputs such as the user-provided GFF.
    * This is done this way because Nextflow doesn't handle optional inputs. One hack that the
    * community uses for inputs of type path is to provide an empty array ([]). So, we first
    * join with user-provided GFF with the remainder, try to get an empty element, and then we use map
    * to transform the null to [].
    ***********************************************************************************************/
    def integrator_ch = PROKKA.out.prokka_gff
        .join(
            RENAME.out.map_file
        )
        .join(
            ISESCAN.out.iss_tsv
        )
        .join(
            INTEGRONFINDER.out.contigs_summary
        )
        .join(
            INTEGRONFINDER.out.contigs_gbks
        )
        .join(
            ICEFINDER2_LITE.out.ices_tsv,
            remainder: true
        )
        .join(
            GENOMAD.out.genomad_vir
        )
        .join(
            GENOMAD.out.genomad_plas
        )
        .join(
            COMPOSITIONAL_OUTLIER_DETECTION.out.bed,
            remainder: true
        )
        .join(
            VIRIFY_QC.out.virify_hq,
            remainder: true
        )

    INTEGRATOR(
        integrator_ch.map { meta, prokka_gff, map_file, iss_tsv, contigs_summary, gbks, ices_tsv, genomad_vir, genomad_plas, compos_bed, virify_hq ->
            [meta, prokka_gff, map_file, iss_tsv, contigs_summary, gbks, ices_tsv ? ices_tsv : [], genomad_vir, genomad_plas, compos_bed ? compos_bed : [], virify_hq ? virify_hq : []]
        }
    )
    ch_versions = ch_versions.mix(INTEGRATOR.out.versions)


    // POSTPROCESSING
    GFF_REDUCE(INTEGRATOR.out.mobilome_prokka_gff)
    ch_versions = ch_versions.mix(GFF_REDUCE.out.versions)

    FASTA_WRITER(
        ch_inputs.map { meta, fasta, _user_proteins_gff, _virify_gff -> [meta, fasta] }.join(GFF_REDUCE.out.mobilome_nogenes)
    )
    ch_versions = ch_versions.mix(FASTA_WRITER.out.versions)

    GFF_MAPPING(
        GFF_REDUCE.out.mobilome_clean.join(user_proteins_ch)
    )
    ch_versions = ch_versions.mix(GFF_MAPPING.out.versions)

    if (params.gff_validation) {
        GT_GFF3VALIDATOR(GFF_REDUCE.out.mobilome_nogenes)
        ch_versions = ch_versions.mix(GT_GFF3VALIDATOR.out.versions)
    }

    // AMRFinder is optional. default skip_amr = FALSE
    def amr_finder_ch = PROKKA.out.prokka_fna.join(PROKKA.out.prokka_faa).join(PROKKA.out.prokka_gff).filter { it -> !it[0].skip_amrfinder_plus }

    AMRFINDER_PLUS(amr_finder_ch)
    ch_versions = ch_versions.mix(AMRFINDER_PLUS.out.versions)

    AMRFINDER_REPORT(
        AMRFINDER_PLUS.out.amrfinder_tsv.join(
            INTEGRATOR.out.mobilome_prokka_gff
        ).join(
            RENAME.out.map_file
        ).join(
            user_proteins_ch
        )
    )
    ch_versions = ch_versions.mix(AMRFINDER_REPORT.out.versions)

    //
    // Collate and save software versions
    //

    // Version collating //
    CUSTOM_DUMPSOFTWAREVERSIONS(
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config = Channel.fromPath(
        "${projectDir}/assets/multiqc_config.yml",
        checkIfExists: true
    )
    ch_multiqc_custom_config = params.multiqc_config
        ? Channel.fromPath(params.multiqc_config, checkIfExists: true)
        : Channel.empty()

    // ch_multiqc_logo = params.multiqc_logo
        // ? Channel.fromPath(params.multiqc_logo, checkIfExists: true)
        // : Channel.fromPath("${projectDir}/assets/mgnify_wordmark_dark_on_light.png", checkIfExists: true)

    MULTIQC(
        CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.first(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        [],
        [],
        [],
    )

    emit:
    versions = ch_versions
}
