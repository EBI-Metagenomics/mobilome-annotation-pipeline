process GFF2GBK {
    tag "$meta.id"
    label 'process_single'

    container 'microbiome-informatics/mgnify-pipelines-toolkit:1.4.19'

    input:
    tuple val(meta), path(contigs), path(gff), path(proteins)

    output:
    tuple val(meta), path("*_from_gff.gbk"), emit: gbk
    tuple val("${task.process}"), val('mgnify-pipelines-toolkit'), eval("get_mpt_version"), topic: versions, emit: versions_mgnifypipelinestoolkit

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def proteins_arg = proteins.size() > 0 ? "--faa ${proteins}" : ""
    """
    gbk_generator \\
        --contigs ${contigs} \\
        --gff ${gff} \\
        --output_gbk ${prefix}_from_gff.gbk \\
        ${proteins_arg} \\
        ${args}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def proteins_arg = proteins.size() > 0 ? "--proteins ${proteins}" : ""
    """
    echo $args

    touch ${prefix}_from_gff.gbk
    """
}
