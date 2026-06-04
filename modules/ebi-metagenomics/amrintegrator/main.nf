process AMRINTEGRATOR {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:1.4.6--pyhdfd78af_0':
        'biocontainers/mgnify-pipelines-toolkit:1.4.6--pyhdfd78af_0' }"

    input:
    tuple val(meta), 
     path(deeparg, stageAs: 'deeparg_hamr.tsv'), 
     path(rgi, stageAs: 'rgi_hamr.tsv'), 
     path(amrfp, stageAs: 'amrfinder.tsv'), 
     path(gff)

    output:
    tuple val(meta), path("*.gff"), optional: true, emit: gff
    path "versions.yml",                            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def deeparg_param = deeparg ? "--deeparg_hamr ${deeparg}" : ""
    def rgi_param = rgi ? "--rgi_hamr ${rgi}" : ""
    def amrfp_param = amrfp ? "--amrfp_out ${amrfp}" : ""
    """
    amr_integrator \\
        ${deeparg_param} \\
        ${rgi_param} \\
        ${amrfp_param} \\
        --cds_gff ${gff} \\
        --output integrated_${prefix}.gff \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: \$(get_mpt_version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args
    
    touch ${prefix}.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: \$(get_mpt_version)
    END_VERSIONS
    """
}
