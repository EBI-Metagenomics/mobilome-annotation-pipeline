process INTEGRATOR {
    tag "$meta.id"
    label 'process_single'

    container 'quay.io/biocontainers/biopython:1.78'

    input:
    tuple val(meta), path(prokka_gff), path(map_file), path(iss_tsv), path(inf_summ), path(inf_gbks), path(icf_summ), path(icf_dr), path(mog_out), path(genomad_vir), path(genomad_plas), file(vir_results), file(crispr_tsv)
    // path(pal_info)

    output:
    tuple val(meta), path("${meta.id}_mobilome_prokka.gff")       , emit: mobilome_prokka_gff
    tuple val(meta), path("${meta.id}_overlapping_integrons.txt") , emit: overlapping_integrons_txt
    tuple val(meta), path("${meta.id}_discarded_mge.txt")         , emit: discarded_mge_txt

    script:
    // def palidis_arg = (pal_info) ? "--pal_tsv ${pal_info}" : ""
    // mge_integrator.py ${palidis_arg} \
    """
    mge_integrator.py \
    --pkka_gff ${prokka_gff} \
    --map ${map_file} \
    --iss_tsv ${iss_tsv} \
    --inf_tsv ${inf_summ} \
    --inf_gbks ${inf_gbks.join(' ')} \
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
