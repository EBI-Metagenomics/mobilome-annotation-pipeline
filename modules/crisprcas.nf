#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process CRISPR_FINDER {

    publishDir "$params.outdir/func_annot", mode: 'copy'

    container 'quay.io/microbiome-informatics/genomes-pipeline.crisprcasfinder:4.3.2'

    input:
        path contigs

    output:
        path("crispr_results/TSV/Crisprs_REPORT.tsv"), emit: crispr_report

    script:
    if(contigs.size() > 0)
       """
       CRISPRCasFinder.pl \
       -i ${contigs} \
       -so ${params.crispr_so} \
       -def G \
       -drpt ${params.crispr_drpt} \
       -outdir crispr_results \
       -metagenome
       """

    else
        """
        echo 'CRISPRCasFinder output file empty due to empty input... generating dummy files'
        mkdir -r crispr_results/TSV && touch crispr_results/TSV/Crisprs_REPORT.tsv 
        """

}
