process GBK_SPLITTER {
    tag "$meta.id"
    label 'process_single'

    input:
    path gbk_file

    output:
    path "*.gbk"
    path "input.list", emit: gbks

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
