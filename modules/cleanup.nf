#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process CLEANUP {
    input:
    path tmp_dir
    path icefinder_summary 

    script:
    """
    echo "Cleaning up ICEFINDER tmp directories..."

    real_tmp_dir=\$(readlink -f "${tmp_dir}")
    ls -la \${real_tmp_dir}
    rm -rf \${real_tmp_dir}
    """
}
