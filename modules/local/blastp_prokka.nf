process BLASTP_PROKKA {
    tag "$meta.id"
    label 'process_medium'


    container 'https://depot.galaxyproject.org/singularity/blast:2.15.0--pl5321h6f7f691_1'


    input:
    tuple val(meta), path(faa)
    tuple val(meta2), path(database)

    output:
    tuple val(meta), path("*_blastp_uniprot.tsv"), emit: uniprot_tsv

    script:
    def blastp_args = task.ext.blastp_args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    blastp \\
        -query ${faa} \\
        -db ${meta2.id} \\
        -num_threads ${task.cpus} \\
        -out ${prefix}_blastp_uniprot.tsv \\
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle" \\
        ${blastp_args}
    """
}
