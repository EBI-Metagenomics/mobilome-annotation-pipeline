// Parse the --annotation_manifest CSV and emit one channel per annotation type.
// Each emitted channel carries only the samples that have a non-empty entry for
// that tool; samples with a blank column are silently filtered out so downstream
// joins receive [] as the absent-input sentinel rather than null.
include { samplesheetToList } from 'plugin/nf-schema'

workflow PARSE_MANIFEST {

    take:
    manifest_path    // val: path string to the annotation manifest CSV
    samplesheet_ids  // val: Set<String> of sample IDs from the input samplesheet

    main:
    def manifest_list = samplesheetToList(manifest_path, "${projectDir}/assets/schema_manifest.json")
    def manifest_ids  = manifest_list.collect { row -> row[0].id } as Set
    def unknown_ids   = manifest_ids - samplesheet_ids
    if (unknown_ids) {
        error "Annotation manifest contains sample IDs not found in the samplesheet: " +
              "${unknown_ids.sort().join(', ')}. Check for typos in the manifest 'sample' column."
    }
    def ch_manifest = Channel.fromList(manifest_list)

    // multiMap broadcasts every row to all five named output channels.
    // Multiple .map{} calls on a single queue channel would round-robin items
    // between consumers, causing each tool channel to receive only a fraction.
    def ch_split = ch_manifest.multiMap { meta, ips_tsv, amrfinder_tsv, antismash_gff, gecco_gff, sanntis_gff ->
        ips_tsv:       tuple(meta, ips_tsv       ?: [])
        amrfinder_tsv: tuple(meta, amrfinder_tsv ?: [])
        antismash_gff: tuple(meta, antismash_gff ?: [])
        gecco_gff:     tuple(meta, gecco_gff     ?: [])
        sanntis_gff:   tuple(meta, sanntis_gff   ?: [])
    }

    emit:
    ips_tsv       = ch_split.ips_tsv.filter       { _meta, f -> f != [] }
    amrfinder_tsv = ch_split.amrfinder_tsv.filter { _meta, f -> f != [] }
    antismash_gff = ch_split.antismash_gff.filter { _meta, f -> f != [] }
    gecco_gff     = ch_split.gecco_gff.filter     { _meta, f -> f != [] }
    sanntis_gff   = ch_split.sanntis_gff.filter   { _meta, f -> f != [] }
}
