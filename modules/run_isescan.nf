#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process isescan {
    publishDir "$launchDir/", mode: 'copy'
    stageInMode = 'copy'

    memory "8 GB"
    cpus 4

    container "quay.io/microbiome-informatics/isescan-v1.7.2.3"

    input:
        path assembly_file, name: 'contigs.fasta'

    output:
        path 'isescan_results/contigs.fasta.is.fna', emit: iss_fasta 
        path 'isescan_results/contigs.fasta.tsv', emit: iss_tsv

    script:
    if (assembly_file.size() > 0)
        """
        isescan.py --seqfile contigs.fasta \
        --output isescan_results \
        --nthread ${task.cpus}

        if ls -l isescan_results/contigs.fasta.tsv 2>/dev/null | grep -q .
        then
            echo 'ISEScan results exists'
        else
	    echo 'ISEScan found 0 insertion sequences in input file... generating dummy files'
            touch isescan_results/contigs.fasta.is.fna
            touch isescan_results/contigs.fasta.tsv
        fi
        """
    else
        """
	echo 'ISEScan is not running due to empty input... generating dummy files'
        mkdir isescan_results
	touch isescan_results/contigs.fasta.is.fna
        touch isescan_results/contigs.fasta.tsv
        """
}

