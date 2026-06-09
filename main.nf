/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ebi-metagenomics/mobilome-annotation-pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/ebi-metagenomics/mobilome-annotation-pipeline
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MOBILOMEANNOTATION } from './workflows/mobilomeannotation'
include { DOWNLOAD_DATABASES } from './subworkflows/local/download_databases'

//
// WORKFLOW: Run main ebi-metagenomics/mobilome-annotation-pipeline analysis pipeline
//
workflow EBIMETAGENOMICS {
    MOBILOMEANNOTATION()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    if (params.download_dbs) {
        DOWNLOAD_DATABASES()

        workflow.onComplete {
            def db_dir = params.download_dbs
            if (workflow.success) {
                log.info """
                ============================================================
                 Database download complete!
                 Add the following block to your config file (e.g. my_paths.config)
                 and pass it with: nextflow run ... -c my_paths.config
                ============================================================
                params {
                    // Mobilome databases
                    genomad_db                   = "${db_dir}/genomad_db_v1.9"
                    icefinder_macsyfinder_models = "${db_dir}/icf2_dbs/macsydata"
                    icefinder_hmm_models         = "${db_dir}/icf2_dbs/icehmm/icescan"
                    icefinder_prokka_uniprot_db  = "${db_dir}/icf2_dbs/icefinder_prokka_uniprot"

                    // PATHOFACT2
                    pathofact_models             = "${db_dir}/Models.tar.gz"
                    virulencefactors_db          = "${db_dir}/VFDB_setB_pro.dmnd"
                    ncbi_cdd                     = "${db_dir}/database"

                    // AMR
                    amrfinderplus_db             = "${db_dir}/amrfinderdb"
                    deeparg_db                   = "${db_dir}/db"
                    rgi_db                       = "${db_dir}/card_dir"

                    // BGC
                    antismash_db                 = "${db_dir}/antismash_db"
                }
                ============================================================
                """.stripIndent()
            }
        }
    } else {
        EBIMETAGENOMICS()
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
