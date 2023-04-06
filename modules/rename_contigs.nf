#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process rename {
    publishDir "$launchDir/preprocessing"

    container "quay.io/biocontainers/biopython:1.75"

    memory "1 GB"
    cpus 1

    input:
      path assembly_file, name: 'contigs.fasta'
      
    output:
      path '1kb_contigs.fasta', emit: contigs_1kb
      path '5kb_contigs.fasta', emit: contigs_5kb
      path 'contigID.map', emit: map_file
    
    """    
    assembly_filter_rename.py --assembly contigs.fasta
    """
}

