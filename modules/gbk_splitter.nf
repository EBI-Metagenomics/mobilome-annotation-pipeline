#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process GBK_SPLITTER {

    publishDir "$params.outdir/prediction/icefinder_results/gbk", pattern: '*.gbk', mode: 'copy'
    publishDir "$params.outdir/prediction/icefinder_results/", pattern: 'input.list', mode: 'copy'

    input:
    path gbk_file
      
    output:
    path "*.gbk"
    path "input.list" , emit: gbks

    // FIXME: this is no very canonical (in terms of nextflow pipelines)
    //        also, I don't think this is going to work in the cloud
    tmpDir = file("$params.outdir/prediction/icefinder_results/tmp")
    tmpDir.mkdirs()

    resDir = file("$params.outdir/prediction/icefinder_results/result")
    resDir.mkdirs()

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

