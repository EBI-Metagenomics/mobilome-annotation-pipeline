process ICEFINDER {
    tag "${meta.id}"

    container {
        if (params.icefinder_sif) {
            return params.icefinder_sif
        }
        return "${params.singularity_cachedir}/icefinder-v1.0-local.sif"
    }

    containerOptions {
        def args = [
            "--bind ${input_list}:/install/ICEfinder_linux/input.list",
            "--bind ${task.workDir}/gbks:/install/ICEfinder_linux/gbk/",
            "--bind ${task.workDir}/tmp:/install/ICEfinder_linux/tmp/",
            "--bind ${task.workDir}/result:/install/ICEfinder_linux/result/",
            "--pwd /install/ICEfinder_linux/"
        ]
        args.join(" ")
    }

    input:
    tuple val(meta), path(input_list), path(gbk_files, stageAs: "gbks/*")
    // path tmp_folder, stageAs: "tmp"
    // path res_folder, stageAs: "result"

    output:
    tuple val(meta), path("result/icf_concat.summary"), emit: icf_summ_files
    tuple val(meta), path("result/icf_dr.txt"), emit: icf_dr

    script:
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
}
