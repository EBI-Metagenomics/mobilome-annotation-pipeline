process GBK_SPLITTER {

    publishDir "${params.outdir}/prediction/icefinder_results/gbk", pattern: '*.gbk', mode: 'copy'
    publishDir "${params.outdir}/prediction/icefinder_results/", pattern: 'input.list', mode: 'copy'

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
