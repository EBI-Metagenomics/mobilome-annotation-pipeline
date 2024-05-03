#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process GBK_SPLITTER {

    // TODO: Ale, do we need this? Or with the icefinder publishDir bits is enough?
    // publishDir "$params.outdir/prediction/icefinder_results/gbk", pattern: '*.gbk', mode: 'copy'
    // publishDir "$params.outdir/prediction/icefinder_results/", pattern: 'input.list', mode: 'copy'

    input:
    path gbk_file
      
    output:
    path "input.list" , emit: input_list
    path "*.gbk"      , emit: gbks

    script:
    if (gbk_file.size() > 0)
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
    else
        """
        echo 'ICEfinder dir empty due to empty input... generating dummy files'
        touch input.list
        touch dummy.gbk
        """
}

