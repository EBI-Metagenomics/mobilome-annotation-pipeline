#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process ISESCAN {
    publishDir "$params.outdir/prediction"

    cpus 8
    memory { 8.GB * task.attempt }
    errorStrategy 'retry' 
    maxRetries 3

    container 'quay.io/microbiome-informatics/isescan-v1.7.2.3'

    input:
        path assembly_file

    output:
        path 'isescan_results/1kb_contigs.fasta.tsv', emit: iss_tsv

    script:
    if (assembly_file.size() > 0)
        """
        isescan.py --seqfile ${assembly_file} \
        --output isescan_results \
        --nthread ${task.cpus}

        if ls -l isescan_results/1kb_contigs.fasta.tsv 2>/dev/null | grep -q .
        then
            echo 'ISEScan results exists'
        else
	    echo 'ISEScan found 0 insertion sequences in input file... generating dummy files'
            touch isescan_results/1kb_contigs.fasta.tsv
        fi
        """
    else
        """
	echo 'ISEScan is not running due to empty input... generating dummy files'
        mkdir isescan_results
        touch isescan_results/1kb_contigs.fasta.tsv
        """
}

