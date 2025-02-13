process INTEGRATOR {
    tag "$meta.id"
    label 'process_single'

    container 'quay.io/biocontainers/biopython:1.78'

    input:
    tuple val(meta), path(prokka_gff), path(map_file), path(iss_tsv), path(pal_info), path(inf_summ), path(inf_gbk), path(icf_summ), path(icf_dr), path(mog_out), path(genomad_vir), path(genomad_plas), path(vir_results), path(crispr_tsv)

    output:
    tuple val(meta), path("${meta.id}_mobilome_prokka.gff")       , emit: mobilome_prokka_gff
    tuple val(meta), path("${meta.id}_overlapping_integrons.txt") , emit: overlapping_integrons_txt
    tuple val(meta), path("${meta.id}_discarded_mge.txt")         , emit: discarded_mge_txt

    script:
    """
    mge_integrator.py \
    --pkka_gff ${prokka_gff} \
    --map ${map_file} \
    --iss_tsv ${iss_tsv} \
    --pal_tsv ${pal_info} \
    --inf_tsv ${inf_summ} \
    --inf_gbks ${inf_gbk.join(' ')} \
    --icf_tsv ${icf_summ} \
    --icf_lim ${icf_dr} \
    --mog_tsv ${mog_out} \
    --geno_out ${genomad_vir} \
    --geno_plas ${genomad_plas} \
    --virify_out ${vir_results} \
    --crispr_out ${crispr_tsv} \
    --prefix ${meta.id}
    """
}
