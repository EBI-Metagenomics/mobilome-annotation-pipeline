// Subworkflow to generate toxins and virulence factors annotation from protein sequences
// Outputs are filtered by threshold and integrated into a single GFF3 format output
include { PATHOFACT2_DOWNLOADDATA  } from '../../../modules/ebi-metagenomics/pathofact2/downloaddata/main'
include { PATHOFACT2_TOXINS        } from '../../../modules/ebi-metagenomics/pathofact2/toxins/main'
include { PATHOFACT2_VIRULENCE     } from '../../../modules/ebi-metagenomics/pathofact2/virulence/main'
include { PATHOFACT2_INTEGRATOR    } from '../../../modules/ebi-metagenomics/pathofact2/integrator/main'
include { PATHOFACT2_EXTRACTFASTA  } from '../../../modules/ebi-metagenomics/pathofact2/extractfasta/main'
include { LOCALCDSEARCH_ANNOTATE   } from '../../../modules/nf-core/localcdsearch/annotate/main'
include { LOCALCDSEARCH_DOWNLOAD   } from '../../../modules/nf-core/localcdsearch/download/main'
include { DIAMOND_BLASTP           } from '../../../modules/nf-core/diamond/blastp/main'
include { DIAMOND_MAKEDB           } from '../../../modules/nf-core/diamond/makedb/main'
include { WGET                     } from '../../../modules/nf-core/wget/main'

workflow PATHOFACT2 {
    take:
    ch_inputs          // channel: tuple( val(meta), path(aminoacids), path(cds_gff), path(ips_tsv) )
    ch_models          // channel: path( pathofact2_db )
    ch_vfdb            // channel: path( vfdb )
    ch_cdd             // channel: path( cdd_db )
    ch_zenodo_id        // channel: value( pathofact2_db_zenodo_id )
    ch_vfdb_url        // channel: tuple( val(meta2), val(vfdb_url) )

    main:
    ch_versions = channel.empty()

    // Extract individual components from input channel
    ch_faa = ch_inputs.map{ meta, aminoacids, _cds_gff, _ips_tsv -> tuple(meta, aminoacids) }
    ch_gff = ch_inputs.map{ meta, _aminoacids, cds_gff, _ips_tsv -> tuple(meta, cds_gff) }
    ch_ips = ch_inputs.map{ meta, _aminoacids, _cds_gff, ips_tsv -> tuple(meta, ips_tsv) }

    // Split inputs based on whether IPS annotation is provided
    ch_ips
        .branch { meta, ips_tsv ->
            with_ips: ips_tsv
                return tuple(meta, ips_tsv)
            without_ips: !ips_tsv
                return meta
        }
        .set { ch_ips_branched }

    ch_with_ips = ch_ips_branched.with_ips
    ch_without_ips = ch_ips_branched.without_ips

    // Preparing databases
    if (ch_models) {
        pathofact_models = ch_models
    } else {
        PATHOFACT2_DOWNLOADDATA(ch_zenodo_id)
        pathofact_models = PATHOFACT2_DOWNLOADDATA.out.zenodo_file
    }

    if (ch_vfdb) {
        vfdb_diamond_db = ch_vfdb
    } else {
        WGET(ch_vfdb_url)
        ch_versions = ch_versions.mix(WGET.out.versions.first())
        DIAMOND_MAKEDB(WGET.out.outfile, [], [], [])
        ch_versions = ch_versions.mix(DIAMOND_MAKEDB.out.versions.first())
        vfdb_diamond_db = DIAMOND_MAKEDB.out.db
    }

    // Prepare CDD database (will only be used if ch_without_ips has data)
    if (ch_cdd) {
        cdd_database = ch_cdd
    } else {
        LOCALCDSEARCH_DOWNLOAD(['cdd_ncbi'])
        cdd_database = LOCALCDSEARCH_DOWNLOAD.out.db
    }

    // Running prediction
    PATHOFACT2_TOXINS( ch_faa, pathofact_models )

    PATHOFACT2_VIRULENCE( ch_faa, pathofact_models )

    // Searching for hits in VFDB
    DIAMOND_BLASTP( ch_faa, vfdb_diamond_db, 6, 'qseqid sseqid pident length qlen slen evalue bitscore')
    ch_versions = ch_versions.mix(DIAMOND_BLASTP.out.versions.first())

    // Extracting positive matches
    ch_extractfasta_input = ch_faa
        .join(DIAMOND_BLASTP.out.txt)
        .join(PATHOFACT2_TOXINS.out.tsv)
        .join(PATHOFACT2_VIRULENCE.out.tsv)
    PATHOFACT2_EXTRACTFASTA(ch_extractfasta_input)

    // Running annotation using local-cd-search when ips_tsv is not provided
    ch_fasta_for_cdd = PATHOFACT2_EXTRACTFASTA.out.fasta
        .join(ch_without_ips)
    LOCALCDSEARCH_ANNOTATE(ch_fasta_for_cdd, cdd_database, false)

    // Combine IPS annotations with CDD annotations
    prot_annot = ch_with_ips.mix(LOCALCDSEARCH_ANNOTATE.out.result)

    // Set annotation type based on source
    annot_type = ch_with_ips
        .map { meta, _ips_tsv -> tuple(meta, 'ips') }
        .mix(
            LOCALCDSEARCH_ANNOTATE.out.result.map { meta, _annot -> tuple(meta, 'cdd') }
        )

    // Integrating results in a single gff file
    ch_for_integrator = ch_gff
        .join(prot_annot)
        .join(PATHOFACT2_EXTRACTFASTA.out.tsv)
        .join(annot_type)
    PATHOFACT2_INTEGRATOR(ch_for_integrator)

    // Handle cases where no predictions are made (integrator produces no output)
    ch_gff_output = PATHOFACT2_INTEGRATOR.out.gff.ifEmpty([])

    emit:
    gff  =  ch_gff_output            // channel: tuple( val(meta), path(gff) )
    versions = ch_versions           // channel: [ versions.yml ]

}
