#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process INTEGRONFINDER {

    publishDir "$params.outdir/prediction/integron_results", mode: 'copy'

    cpus 8
    memory { 8.GB * task.attempt }
    errorStrategy 'retry' 
    maxRetries 3

    container 'quay.io/microbiome-informatics/integronfinder:latest'

    input:
        path assembly_file

    output:
	path "Results_Integron_Finder_5kb_contigs/5kb_contigs.summary", emit: inf_summ
	path "Results_Integron_Finder_5kb_contigs/contig_*.gbk", emit: inf_gbk

    script:
    if(assembly_file.size() > 0)
        """
        integron_finder --union-integrases \
        --mute \
        --local-max \
        --cpu ${task.cpus} \
        --func-annot \
        --gbk \
        ${assembly_file}

	    if ls -l Results_Integron_Finder_5kb_contigs/contig_*.gbk 2>/dev/null | grep -q .
        then
	        echo 'IntegronFinder outputs complete'
        else
            echo 'IntegronFinder found 0 integrons in assembly... generating dummy files'
            touch Results_Integron_Finder_5kb_contigs/contig_dummy.gbk
        fi
        """
    else
        """
        echo 'IntegronFinder dir empty due to empty input... generating dummy files'
        mkdir Results_Integron_Finder_5kb_contigs
        touch Results_Integron_Finder_5kb_contigs/5kb_contigs.summary
        touch Results_Integron_Finder_5kb_contigs/contig_dummy.gbk
        """
}

