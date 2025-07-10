process ARAGORN {
    tag "$meta.id"
    label 'process_single'

    container 'quay.io/biocontainers/aragorn:1.2.41--h7b50bb2_5'

    input:
    tuple val(meta), path(fasta_file)

    output:
    path "*_aragorn.gff", emit: trna_gff
    path "*_aragorn.tbl", emit: trna_tbl

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def aragorn_args = task.ext.aragorn_args ?: ''
    """
    aragorn \\
        ${aragorn_args} \\
        -o ${prefix}_aragorn.tbl \\
        ${fasta_file}
    
    # Convert into GFF format 
    awk 'BEGIN { OFS="\t" } /^>/ {
        match(\$0, /contig=([^ ]+)/, m)
        contig = m[1]
        feature_count = 0
        next
    }
    /genes found/ { next }

    /^[0-9]+/ {
        feature_count++
        type = $2
        match(\$0, /\[([0-9]+),([0-9]+)\]/, coords)
        start = coords[1]
        end = coords[2]
        strand = (start <= end ? "+" : "-")
        print contig, "aragorn", "tRNA", start, end, ".", strand, ".", "ID=tRNA" feature_count ";product=" type
    }' ${prefix}_aragorn.tbl > ${prefix}_aragorn.gff
    """
}

