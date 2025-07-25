include { OUTLIER_FINDER } from '../modules/compositional_outlier_detection'
include { MERGE_RESULTS  } from '../modules/merge_results'

workflow COMPOSITIONAL_OUTLIER_DETECTION {
    take:
        ch_assembly  // tuple(meta, 100kb_contigs)

    main:
        // Split FASTA files into chunks with meta naming
        ch_chunks = ch_assembly
            .splitFasta( by:20, file:true )

        // Run compositional outlier detection on each chunk
        OUTLIER_FINDER (
            ch_chunks
        )

        // Group results while preserving the full meta map
        ch_grouped_results = OUTLIER_FINDER.out.bed.groupTuple(by: 0)

        // Merge results from all chunks for each sample
        MERGE_RESULTS (
            ch_grouped_results
        )

    emit:
        bed = MERGE_RESULTS.out.bed
}

