process CRISPR_FINDER {

    // TODO: publish dir statements go in the modules.config
    // publishDir "${params.outdir}/func_annot", mode: 'copy'

    container 'quay.io/microbiome-informatics/genomes-pipeline.crisprcasfinder:4.3.2patchedv5'

    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path("${meta.id}_crisprs_report.tsv "), emit: crispr_report

    script:
    """
    # CRISPRCasFinder doesn't like it if the folder is there already, which could happen
    # when retrying this process
    rm -rf crispr_results || true

    CRISPRCasFinder.pl \
    -i ${contigs} \
    -so ${params.crispr_so} \
    -def G \
    -drpt ${params.crispr_drpt} \
    -outdir crispr_results \
    -metagenome

    cp crispr_results/TSV/Crisprs_REPORT.tsv ${meta.id}_crisprs_report.tsv 
    """
}
