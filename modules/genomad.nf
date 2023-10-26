#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process GENOMAD {
    publishDir "$params.outdir/prediction/genomad_out"

    cpus 8
    memory { 32.GB * task.attempt }
    errorStrategy 'retry' 
    maxRetries 3

    container 'quay.io/biocontainers/genomad:1.6.1--pyhdfd78af_0'

    input:
        path assembly_file

    output:
        path("genomad_out/5kb_contigs_summary/5kb_contigs_virus_summary.tsv"), emit: genomad_vir
	path("genomad_out/5kb_contigs_summary/5kb_contigs_plasmid_summary.tsv"), emit: genomad_plas

    script:
    if(assembly_file.size() > 0)
        """    
	genomad end-to-end \
	--threads ${task.cpus} \
	${assembly_file} \
	genomad_out \
	${params.genomad_db}
        """
    else
        """
        echo 'geNomad output file empty due to empty input... generating dummy files'
	mkdir -p genomad_out/5kb_contigs_summary
        touch genomad_out/5kb_contigs_summary/5kb_contigs_virus_summary.tsv
	touch genomad_out/5kb_contigs_summary/5kb_contigs_plasmid_summary.tsv
        """
}

