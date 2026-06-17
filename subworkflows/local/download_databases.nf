include { DB_DOWNLOAD_MOBILOME_DBS              } from '../../modules/local/db_download_mobilome_dbs'
include { DB_DOWNLOAD_VFDB                      } from '../../modules/local/db_download_vfdb'
include { DB_UNTAR_AMRFINDERPLUS                } from '../../modules/local/db_untar_amrfinderplus'
include { DB_UNTAR_PATHOFACT                    } from '../../modules/local/db_untar_pathofact'
include { PATHOFACT2_DOWNLOADDATA               } from '../../modules/ebi-metagenomics/pathofact2/downloaddata/main'
include { AMRFINDERPLUS_UPDATE                  } from '../../modules/nf-core/amrfinderplus/update/main'
include { DEEPARG_DOWNLOADDATA                  } from '../../modules/nf-core/deeparg/downloaddata/main'
include { ANTISMASH_ANTISMASHDOWNLOADDATABASES  } from '../../modules/nf-core/antismash/antismashdownloaddatabases/main'
include { LOCALCDSEARCH_DOWNLOAD                } from '../../modules/nf-core/localcdsearch/download/main'
include { RGI_DOWNLOADDB                        } from '../../modules/ebi-metagenomics/rgi/downloaddb/main'

workflow DOWNLOAD_DATABASES {
    main:

    DB_DOWNLOAD_MOBILOME_DBS()
    DB_DOWNLOAD_VFDB()
    PATHOFACT2_DOWNLOADDATA(Channel.of(params.zenodo_id ?: 18223764))
    DB_UNTAR_PATHOFACT(PATHOFACT2_DOWNLOADDATA.out.zenodo_file)
    LOCALCDSEARCH_DOWNLOAD(Channel.of(['cdd']))
    AMRFINDERPLUS_UPDATE()
    DB_UNTAR_AMRFINDERPLUS(AMRFINDERPLUS_UPDATE.out.db)
    DEEPARG_DOWNLOADDATA()
    RGI_DOWNLOADDB()
    ANTISMASH_ANTISMASHDOWNLOADDATABASES()
}
