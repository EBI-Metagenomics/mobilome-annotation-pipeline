#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process VIRIFY_QC {
    publishDir "$launchDir/$params.outdir/prediction"

    container 'quay.io/microbiome-informatics/virify-python3:1.2'

    input:
        path vir_gff
	path vir_checkv

    output:
        path("virify_hq.gff"), emit: virify_hq

    script:
    if(vir_gff.exists())
        """    
        virify_qc.py \
        --virify_gff ${vir_gff} \
	--checkv_out ${vir_checkv.join(' ')}
        """	
    else
        """
        echo 'No input files for VIRify parsing... generating dummy files'
        touch virify_hq.gff
        """

}




