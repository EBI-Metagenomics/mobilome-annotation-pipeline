process ANTISMASH_JSON2GFF {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:1.4.16--pyhdfd78af_0' :
        'biocontainers/mgnify-pipelines-toolkit:1.4.16--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(antismash_json)

    output:
    tuple val(meta), path("*_antismash.gff"), emit: gff
    tuple val("${task.process}"), val('mgnify-pipelines-toolkit'), eval("get_mpt_version"), topic: versions, emit: versions_mgnifypipelinestoolkit

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    antismash_gff_builder \\
        ${args} \\
        --input ${antismash_json} \\
        --output ${prefix}_antismash.gff
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    touch ${prefix}_antismash.gff
    """
}
