process GBK_SPLITTER {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(gbk_file)

    output:
    tuple val(meta), path("*.gbk"),      emit: gbks
    tuple val(meta), path("input.list"), emit: intput_list

    script:
    """    
    gbk_splitter.pl ${gbk_file}

    if -s input.list
    then
        echo 'The file is not empty'
    else
        echo 'No contigs of size > 5kb ... generating dummy files'
        touch dummy.gbk
    fi
    """
}
