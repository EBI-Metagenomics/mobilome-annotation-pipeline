process DIAMOND {
    tag "$meta.id"
    label 'process_low'

    container 'quay.io/biocontainers/diamond:2.0.12--hdcc8f71_0'

    input:
    tuple val(meta), path(proteins_file)
    path diamond_db 

    output:
    tuple val(meta), path("${meta.id}_mobileog_hits.tsv"), emit: blast_out

    script:
    """
    diamond blastp \
    -q ${proteins_file} \
    --db ${diamond_db} \
    -k 15 \
    -e 1e-20 \
    --query-cover 90 \
    --id 90 \
    --threads ${task.cpus} \
    --outfmt 6 stitle qtitle pident bitscore slen evalue qlen sstart send qstart qend \
    -o ${meta.id}_mobileog_hits.tsv
    """
}

