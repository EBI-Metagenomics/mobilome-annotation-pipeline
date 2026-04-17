process BGCSMAPPER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:1.4.20--pyhdfd78af_0':
        'biocontainers/mgnify-pipelines-toolkit:1.4.20--pyhdfd78af_0' }"


    input:
    tuple val(meta), path(gff), path(sanntis), path(gecco), path(antismash)

    output:
    tuple val(meta), path("*_bgcs.gff"), optional: true, emit: gff
    tuple val(meta), path("*_bgcs.json"), optional: true, emit: json
    tuple val("${task.process}"), val('mgnify-pipelines-toolkit'), eval("get_mpt_version"), topic: versions, emit: versions_mgnifypipelinestoolkit

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sanntis_param = sanntis ? "--sanntis_gff ${sanntis}" : ""
    def gecco_param = gecco ? "--gecco_gff ${gecco}" : ""
    def antismash_param = antismash ? "--antismash_gff ${antismash}" : ""
    """
    bgc_mapper \\
        ${sanntis_param} \\
        ${gecco_param} \\
        ${antismash_param} \\
        --base_gff ${gff} \\
        --output_gff ${prefix}_bgcs.gff \\
        $args
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    touch ${prefix}_bgcs.gff
    touch ${prefix}_bgcs.json
    """
}
