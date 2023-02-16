#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process integronfinder{
    publishDir "$launchDir/integron_results", mode: 'copy'
    stageInMode = 'copy'

    memory "8 GB"
    cpus 8

    container "quay.io/microbiome-informatics/integronfinder:latest"

    input:
        path assembly_file, name: 'contigs.fasta'

    output:
	path "Results_Integron_Finder_contigs/contigs.summary", emit: inf_summ
	path "Results_Integron_Finder_contigs/contig_*.gbk", emit: inf_gbk

    script:
    if(assembly_file.size() > 0)
        """
        integron_finder --union-integrases \
        --mute \
        --local-max \
        --cpu ${task.cpus} \
        --func-annot \
        --gbk \
        contigs.fasta

	if ls -l Results_Integron_Finder_contigs/contig_*.gbk 2>/dev/null | grep -q .
        then
	    echo 'IntegronFinder outputs complete'
        else
            echo 'IntegronFinder found 0 integrons in assembly... generating dummy files'
            touch Results_Integron_Finder_contigs/contig_dummy.gbk
        fi
        """
    else
        """
        echo 'IntegronFinder dir empty due to empty input... generating dummy files'
        mkdir Results_Integron_Finder_contigs
        touch Results_Integron_Finder_contigs/contigs.summary
        touch Results_Integron_Finder_contigs/contig_dummy.gbk
        """
}

