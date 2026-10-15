include { DB_DOWNLOAD_ICEFINDER2                } from '../../modules/local/db_download_icefinder2'
include { ZENODO_DOWNLOAD as ZENODO_DOWNLOAD_GENOMAD   } from '../../modules/local/zenodo_download'
include { ZENODO_DOWNLOAD as ZENODO_DOWNLOAD_PATHOFACT } from '../../modules/local/zenodo_download'
include { UNTAR as UNTAR_GENOMAD                 } from '../../modules/nf-core/untar/main'
include { CHECKV_DOWNLOADDATABASE               } from '../../modules/nf-core/checkv/downloaddatabase/main'
include { DIAMOND_MAKEDB                        } from '../../modules/nf-core/diamond/makedb/main'
include { DB_UNTAR_AMRFINDERPLUS                } from '../../modules/local/db_untar_amrfinderplus'
include { DB_UNTAR_PATHOFACT                    } from '../../modules/local/db_untar_pathofact'
include { AMRFINDERPLUS_UPDATE                  } from '../../modules/nf-core/amrfinderplus/update/main'
include { DEEPARG_DOWNLOADDATA                  } from '../../modules/nf-core/deeparg/downloaddata/main'
include { ANTISMASH_ANTISMASHDOWNLOADDATABASES  } from '../../modules/nf-core/antismash/antismashdownloaddatabases/main'
include { LOCALCDSEARCH_DOWNLOAD                } from '../../modules/nf-core/localcdsearch/download/main'
include { RGI_DOWNLOADDB                        } from '../../modules/ebi-metagenomics/rgi/downloaddb/main'

workflow DOWNLOAD_DATABASES {
    main:

    // geNomad DB is pinned to a Zenodo record (file index set in conf/download_dbs.config)
    ZENODO_DOWNLOAD_GENOMAD(Channel.of([[id: 'genomad'], params.genomad_zenodo_id, "genomad_db_v${params.genomad_db_version}.tar.gz"]))
    UNTAR_GENOMAD(ZENODO_DOWNLOAD_GENOMAD.out.file)
    DB_DOWNLOAD_ICEFINDER2(params.icefinder2_db_url)
    CHECKV_DOWNLOADDATABASE()
    // Nextflow stages the remote FASTA itself: the diamond container ships no curl.
    // meta.id doubles as the DB prefix, giving VFDB_setB_pro.dmnd.
    DIAMOND_MAKEDB(
        Channel.of([[id: 'VFDB_setB_pro'], file(params.vfdb_url)]),
        [],
        [],
        []
    )
    ZENODO_DOWNLOAD_PATHOFACT(Channel.of([[id: 'pathofact'], params.pathofact_zenodo_id, 'Models.tar.gz']))
    DB_UNTAR_PATHOFACT(ZENODO_DOWNLOAD_PATHOFACT.out.file.map { _meta, tarball -> tarball })
    LOCALCDSEARCH_DOWNLOAD(Channel.of(['cdd']))
    AMRFINDERPLUS_UPDATE()
    DB_UNTAR_AMRFINDERPLUS(AMRFINDERPLUS_UPDATE.out.db)
    DEEPARG_DOWNLOADDATA()
    RGI_DOWNLOADDB()
    ANTISMASH_ANTISMASHDOWNLOADDATABASES()
}
