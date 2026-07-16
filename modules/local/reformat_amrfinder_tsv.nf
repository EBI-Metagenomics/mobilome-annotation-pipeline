// Reformat an AMRFinderPlus TSV produced by the genomes-catalogue-pipeline so that
// its columns align with what amr_integrator expects.
//
// Catalogue column layout:  protein_id=1, drug_class=11, identity=17
// amr_integrator expected:  protein_id=1, drug_class=7,  identity=13
//
// The awk one-liner below extracts those three fields and places them at the
// positions the integrator reads, padding the unused columns with empty strings.
process REFORMAT_AMRFINDER_TSV {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12.12':
        'biocontainers/python:3.12.12' }"

    input:
    tuple val(meta), path(amrfinder_tsv)

    output:
    tuple val(meta), path("*_amrfinder_normalised.tsv"), emit: tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    awk -F"\\t" 'BEGIN{OFS="\\t"} {print \$1,"","","","","",\$11,"","","","","",\$17}' \\
        ${amrfinder_tsv} > ${prefix}_amrfinder_normalised.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_amrfinder_normalised.tsv
    """
}
