process INTEGRATOR {
    tag "${meta.id}"
    label 'process_single'

    container 'quay.io/biocontainers/biopython:1.81'

    input:
    tuple val(meta), path(gff_file), path(map_file), path(iss_tsv), path(inf_summ), path(inf_gbks), path(icf_tsv), path(genomad_vir), path(genomad_plas), path(compos_bed), path(vir_results)

    output:
    tuple val(meta), path("${meta.id}_mobilome.gff")       , emit: mobilome_gff
    tuple val(meta), path("${meta.id}_overlap_report.txt") , emit: overlapping_integrons_txt
    tuple val(meta), path("${meta.id}_discarded_mge.txt")  , emit: discarded_mge_txt
    path "versions.yml", emit: versions

    script:
    def virify_arg = vir_results ? "--virify_out ${vir_results}" : ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mge_integrator.py \\
        --gff_file ${gff_file} \\
        --map ${map_file} \\
        --iss_tsv ${iss_tsv} \\
        --inf_tsv ${inf_summ} \\
        --inf_gbks ${inf_gbks.join(' ')} \\
        --icf_tsv ${icf_tsv} \\
        --geno_out ${genomad_vir} \\
        --geno_plas ${genomad_plas} \\
        --comp_bed ${compos_bed} \\
        ${virify_arg} \\
        --prefix ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
}
