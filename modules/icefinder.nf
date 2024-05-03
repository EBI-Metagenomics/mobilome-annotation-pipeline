#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process ICEFINDER {

    publishDir "$params.outdir/prediction/icefinder_results", mode: 'copy'

    cpus 1
    memory { 16.GB * task.attempt }
    errorStrategy 'retry'
    maxRetries 3

    // This is required as it's the directory mounted in the singularity container
    beforeScript 'mkdir -p if_results'

    container {
        if ( params.icefinder_sif ) {
            return params.icefinder_sif
        }
        return "${params.singularity_cachedir}/icefinder-v1.0-local.sif"
    }

    containerOptions "--bind if_results:/install/ICEfinder_linux/result/"

    input:
    path input_list, stageAs: "/install/ICEfinder_linux/input.list"
    path gbks,       stageAs: "/install/ICEfinder_linux/gbk/*"

    output:
    path "if_results/icf_concat.summary", emit: icf_summ_files
    path "if_results/icf_dr.txt",         emit: icf_dr

    script:
    if(input_list.size() > 0)
        """
        mkdir -p /install/ICEfinder_linux/tmp

        cd /install/ICEfinder_linux/ && perl ./ICEfinder_local.pl input.list

        if ls -ld result/contig* 2>/dev/null | grep -q .
        then
            cat result/*/*summary.txt > result/icf_concat.summary
            grep 'DR:' result/*/ICE* > result/icf_dr.txt
        else
            echo 'ICEfinder found 0 ICE/IME in assembly... generating dummy files'
            touch result/icf_concat.summary
            touch result/icf_dr.txt
        fi
        """
    else
        """
        echo 'No input files for ICEfinder... generating dummy files'
        touch result/icf_concat.summary
        touch result/icf_dr.txt
        """
}
