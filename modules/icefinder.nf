#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process ICEFINDER {
    publishDir "$params.outdir/prediction/icefinder_results", mode: 'copy'

    cpus 1
    memory { 16.GB * task.attempt }
    errorStrategy 'retry'
    maxRetries 3

    container {
        if ( params.icefinder_sif ) {
            return params.icefinder_sif
        }
        return "${params.singularity_cachedir}/icefinder-v1.0-local.sif"
    }

    containerOptions {
        args = [
            "--bind ${input_list}:/install/ICEfinder_linux/input.list",
            "--bind ${gbk_folder}:/install/ICEfinder_linux/gbk/",
            "--bind ${tmp_folder}:/install/ICEfinder_linux/tmp/",
            "--bind ${res_folder}:/install/ICEfinder_linux/result/",
            "--pwd /install/ICEfinder_linux/"
        ]
        args.join(" ")
    }

    input:
    path input_list
    path gbk_folder, stageAs: "gbk"
    path tmp_folder, stageAs: "tmp"
    path res_folder, stageAs: "result"

    output:
    path "result/icf_concat.summary", emit: icf_summ_files
    path "result/icf_dr.txt", emit: icf_dr

    script:
    if(input_list.size() > 0)
        """
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
