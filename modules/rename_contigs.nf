#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process RENAME {

    publishDir "$params.outdir/preprocessing", mode: 'copy'

    container 'quay.io/biocontainers/biopython:1.75'

    input:
    path assembly_file
      
    output:
    path '1kb_contigs.fasta', emit: contigs_1kb
    path '5kb_contigs.fasta', emit: contigs_5kb
    path 'contigID.map', emit: map_file

    script:
    """    
    assembly_filter_rename.py --assembly ${assembly_file}
    """
}

